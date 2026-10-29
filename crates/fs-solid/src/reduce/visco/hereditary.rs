//! Complete equilibrium-mode material valves with actual constitutive memory.
//! The same tensor source, strain projection and exact Maxwell times survive;
//! no frequency-frozen C, scalar loss factor or independent-mode fit is used.

use super::{ModalStaticCompliance, field_identity, invalid, poll};
use crate::linear3::{
    TetAssemblyBudget,
    visco::{
        SolidMeshRef, SolidProjectionBudget, TetViscoField, ViscoElementError, ViscoElements,
        ViscoStateBudget,
    },
};
use crate::reduce::{PreparedValve, ReduceError, ValveCardRequest, ValveElementOrder};
use fs_blake3::{ContentHash, DomainHasher};
use fs_exec::{CancelGate, Cx};
use fs_material::visco_tensor::ProjectedTensorProny;
use fs_modal::ModePair;

mod memory;
pub use memory::HereditaryBranch;

/// Explicit cold assembly, source projection and exact history-compression caps.
#[derive(Clone, Copy, Debug)]
pub struct HereditaryValveControls {
    pub assembly: TetAssemblyBudget,
    pub internal_variables: ViscoStateBudget,
    pub projection: SolidProjectionBudget,
    /// Upper bound on 8*sum(max(source_rows,n)*n*n) for equal-time QR.
    pub maximum_factor_scalar_terms: usize,
    /// Admitted relative error in Phi^T M Phi=I and Phi^T K_eq Phi=diag(lambda).
    /// Finite in (0,1); this is a basis check, not a modal truncation certificate.
    pub maximum_basis_error: f64,
}

/// Full original equilibrium basis and projected tensor material. Metre-valued
/// coordinates x use q=sqrt(M_total)*x, M_x=M_total I, K_eq,x=M_total diag(lambda).
/// Equal source-time factors are orthogonally combined, retaining every loss
/// direction. This exactly supports wholly relaxed or wholly unrelaxed initial
/// histories; arbitrary independently authored quadrature memories are absent.
pub struct HereditaryModalValve {
    source: ViscoElements,
    modes: Vec<ModePair>,
    points: Vec<ProjectedTensorProny>,
    point_elements: Vec<usize>,
    viscosity: Vec<f64>,
    branches: Vec<HereditaryBranch>,
    pressure_projection: Vec<f64>,
    opening_projection: Vec<f64>,
    nodes: usize,
    rest_opening_m: f64,
    width_m: f64,
    total_mass_kg: f64,
    source_id: String,
    identity: ContentHash,
    compliance: ModalStaticCompliance,
    basis_error: f64,
    band: (f64, f64),
    instantaneous: f64,
    // Conservative per-coordinate source-strain ratios for the O(n) hot path.
    // Near a limit, evaluate the original quadrature maps without approximation.
    strain_bound: Vec<f64>,
}
impl HereditaryModalValve {
    pub fn source(&self) -> &ViscoElements {
        &self.source
    }
    /// Original shapes in canonical free-DOF order; phi^T M phi = I.
    pub fn modes(&self) -> &[ModePair] {
        &self.modes
    }
    /// phi^T load_per_pa [m2/sqrt(kg)], positive toward opening.
    /// Its transpose is the swept-volume observation on modal velocity.
    pub fn pressure_projection_m2_per_sqrt_kg(&self) -> &[f64] {
        &self.pressure_projection
    }
    /// Mean opening per mass-normalized coordinate [1/sqrt(kg)].
    /// Exact zero entries remain present and distinct from pressure coupling.
    pub fn opening_projection_per_sqrt_kg(&self) -> &[f64] {
        &self.opening_projection
    }
    pub const fn rest_opening_m(&self) -> f64 {
        self.rest_opening_m
    }
    pub const fn width_m(&self) -> f64 {
        self.width_m
    }
    pub const fn total_mass_kg(&self) -> f64 {
        self.total_mass_kg
    }
    pub fn source_id(&self) -> &str {
        &self.source_id
    }
    pub const fn identity(&self) -> ContentHash {
        self.identity
    }
    pub const fn compliance(&self) -> &ModalStaticCompliance {
        &self.compliance
    }
    pub const fn node_count(&self) -> usize {
        self.nodes
    }

    /// Exact original quadrature laws and strain maps in physical coordinates
    /// x [m], with q = sqrt(total_mass_kg)*x.
    pub fn material_points(&self) -> &[ProjectedTensorProny] {
        &self.points
    }
    pub fn point_elements(&self) -> &[usize] {
        &self.point_elements
    }
    /// Full Kelvin C_x [N s/m]; every off-diagonal entry is retained.
    pub fn viscosity_n_s_m(&self) -> &[f64] {
        &self.viscosity
    }
    pub fn branches(&self) -> &[HereditaryBranch] {
        &self.branches
    }
    pub fn memory_channel_count(&self) -> usize {
        self.branches.iter().map(|b| b.rows().len()).sum()
    }
    /// Largest admitted normalized mass/equilibrium-K discrepancy.
    pub const fn basis_error(&self) -> f64 {
        self.basis_error
    }
    /// Intersection of original source bands [rad/s], including DC.
    pub const fn frequency_band_rad_s(&self) -> (f64, f64) {
        self.band
    }
    /// Conservative bound on every instantaneous retained angular frequency.
    /// Uses the COMPLETE mass-normalized equilibrium plus branch stiffness.
    pub const fn instantaneous_frequency_bound_rad_s(&self) -> f64 {
        self.instantaneous
    }

    /// Original material strain validity at a supplied physical-coordinate
    /// endpoint. This endpoint check is not an enclosure of an entire sample.
    pub fn validate_displacement_m(&self, coordinates: &[f64]) -> Result<(), ReduceError> {
        if coordinates.len() != self.modes.len() || coordinates.iter().any(|v| !v.is_finite()) {
            return Err(invalid(
                "hereditary strain requires all finite physical coordinates",
            ));
        }
        let bound: f64 = coordinates
            .iter()
            .zip(&self.strain_bound)
            .map(|(x, ratio)| x.abs() * ratio)
            .sum();
        // Triangle inequality in the original Mandel/Frobenius norm. Keep
        // headroom for accumulation; uncertain cases use the original owner.
        let rounding = 64.0 * f64::EPSILON * (coordinates.len() + 6) as f64;
        if bound < 1.0 - rounding {
            return Ok(());
        }
        for point in &self.points {
            point
                .strain(coordinates)
                .map_err(ViscoElementError::Projection)?;
        }
        Ok(())
    }

    /// Project an arbitrary nodal force field with the SAME original shapes.
    /// Input N produces generalized force N/sqrt(kg). A dimensionless nodal
    /// force pattern per unit applied N produces the corresponding port map.
    /// Fixed components do no work. No mean-opening division is involved.
    pub fn project_nodal_force(
        &self,
        force_n: &[[f64; 3]],
        gate: &CancelGate,
    ) -> Result<Vec<f64>, ReduceError> {
        poll_gate(gate)?;
        if force_n.len() != self.nodes || force_n.iter().flatten().any(|v| !v.is_finite()) {
            return Err(invalid(
                "complete modal force needs one finite vector per source node",
            ));
        }
        let mut result = Vec::with_capacity(self.modes.len());
        for mode in &self.modes {
            poll_gate(gate)?;
            let value = mode
                .phi
                .iter()
                .zip(&self.source.equilibrium().free_dofs)
                .map(|(phi, &dof)| phi * force_n[dof / 3][dof % 3])
                .sum::<f64>();
            if !value.is_finite() {
                return Err(invalid("complete modal force projection overflow"));
            }
            result.push(value);
        }
        poll_gate(gate)?;
        Ok(result)
    }

    /// Expand original mass-normalized displacements or velocities into all
    /// source nodes. Fixed components remain zero. This is the exact transpose
    /// partner of project_nodal_force, not a separately normalized observer.
    pub fn expand_nodal(
        &self,
        coordinates: &[f64],
        gate: &CancelGate,
    ) -> Result<Vec<[f64; 3]>, ReduceError> {
        poll_gate(gate)?;
        if coordinates.len() != self.modes.len() || coordinates.iter().any(|v| !v.is_finite()) {
            return Err(invalid(
                "complete modal observation needs all finite original coordinates",
            ));
        }
        let mut result = vec![[0.0; 3]; self.nodes];
        for (row, &dof) in self.source.equilibrium().free_dofs.iter().enumerate() {
            poll_gate(gate)?;
            let value = self
                .modes
                .iter()
                .zip(coordinates)
                .map(|(mode, q)| mode.phi[row] * q)
                .sum::<f64>();
            if !value.is_finite() {
                return Err(invalid("complete modal observation overflow"));
            }
            result[dof / 3][dof % 3] = value;
        }
        poll_gate(gate)?;
        Ok(result)
    }
}

fn poll_gate(gate: &CancelGate) -> Result<(), ReduceError> {
    if gate.is_requested() {
        Err(crate::linear3::TetElasticError::Cancelled.into())
    } else {
        Ok(())
    }
}

/// Reduce the original equilibrium pencil, then project every source
/// tensor-Prony branch and Kelvin term on ALL selected modes, including modes
/// orthogonal to the mean opening. Original static compliance, topology,
/// pressure work, material bands, supports and cancellation remain mandatory.
/// The legacy frozen image and independent-mode conversion are separate APIs.
pub fn reduce_hereditary_pressure_loaded_valve(
    request: &ValveCardRequest<'_, TetViscoField<'_>>,
    nodal_load_per_pa_m2: &[[f64; 3]],
    controls: &HereditaryValveControls,
    cx: &Cx<'_>,
) -> Result<HereditaryModalValve, ReduceError> {
    poll(cx)?;
    if request.loss_factor != 0.0
        || !controls.maximum_basis_error.is_finite()
        || controls.maximum_basis_error <= 0.0
        || controls.maximum_basis_error >= 1.0
        || nodal_load_per_pa_m2.len() != request.nodes_m.len()
        || nodal_load_per_pa_m2
            .iter()
            .flatten()
            .any(|x| !x.is_finite())
    {
        return Err(invalid(
            "hereditary valve requires zero authored loss, finite source load and explicit basis accuracy",
        ));
    }
    crate::reduce::validate_element_topology(request)?;
    let mesh = match request.element_order {
        ValveElementOrder::Tet4 => SolidMeshRef::Tet4 {
            nodes_m: request.nodes_m,
            elements: request.tetrahedra,
        },
        ValveElementOrder::Tet10 {
            elements,
            quadrature,
        } => SolidMeshRef::Tet10 {
            nodes_m: request.nodes_m,
            elements,
            quadrature,
        },
    };
    let source = ViscoElements::new(
        mesh,
        request.material,
        request.fixed_dofs,
        controls.assembly,
        controls.internal_variables,
        cx,
    )?;
    // The equilibrium compliance and mode basis must lie in every original
    // source's domain. This query enforces DC validity without inventing K.
    source.complex_stiffness(0.0, cx)?;
    let mut reduced = crate::reduce::reduce_valve_impl(
        request,
        None,
        field_identity(request.material, request.tetrahedra.len(), 0.0, cx)?,
        cx,
        Some(nodal_load_per_pa_m2),
        Some(PreparedValve {
            assembly: source.equilibrium(),
            observation_nodes_m: request.nodes_m,
            static_stiffness: None,
            material_loss: None,
            complete_basis: true,
        }),
    )?;
    let modes = reduced
        .complete_material_modes
        .take()
        .ok_or_else(|| invalid("hereditary valve lost its complete equilibrium basis"))?;
    let n = modes.len();
    if !(1..=256).contains(&controls.projection.max_coordinates)
        || n == 0
        || n > controls.projection.max_coordinates
        || n.checked_mul(source.equilibrium().free_dofs.len())
            .is_none_or(|needed| needed > controls.projection.max_basis_scalars)
    {
        return Err(invalid(
            "hereditary complete basis exceeds its explicit capacity",
        ));
    }
    let mass = reduced.card.total_mass_kg;
    let root = mass.sqrt();
    if !root.is_finite() || root <= 0.0 {
        return Err(invalid("hereditary coordinate mass is unrepresentable"));
    }
    let columns = modes
        .iter()
        .map(|m| m.phi.iter().map(|v| v * root).collect::<Vec<_>>())
        .collect::<Vec<_>>();
    if columns
        .iter()
        .flatten()
        .zip(modes.iter().flat_map(|m| &m.phi))
        .any(|(actual, original)| !actual.is_finite() || (*original != 0.0 && *actual == 0.0))
    {
        return Err(invalid("hereditary physical basis is unrepresentable"));
    }
    let projected = source.project_basis(&columns, controls.projection, cx)?;
    let mut viscosity = vec![0.0; n * n];
    let mut equilibrium = vec![0.0; n * n];
    let mut band = (0.0_f64, f64::INFINITY);
    for point in projected.material_points() {
        poll(cx)?;
        let source_band = point.material().band_rad_s();
        band.0 = band.0.max(source_band.0);
        band.1 = band.1.min(source_band.1);
        for k in 0..n * n {
            viscosity[k] += point.viscosity_matrix()[k];
            equilibrium[k] += point.equilibrium_stiffness()[k];
        }
    }
    let mut basis_error = 0.0_f64;
    for i in 0..n {
        for j in 0..n {
            let target_m = if i == j { mass } else { 0.0 };
            let target_k = if i == j { mass * modes[i].lambda } else { 0.0 };
            let stiffness_scale = root * modes[i].lambda.sqrt() * root * modes[j].lambda.sqrt();
            let mass_error = ((projected.mass_matrix_kg()[i * n + j] - target_m) / mass).abs();
            let stiffness_error = ((equilibrium[i * n + j] - target_k) / stiffness_scale).abs();
            if !mass_error.is_finite()
                || !stiffness_error.is_finite()
                || !stiffness_scale.is_finite()
                || stiffness_scale <= 0.0
            {
                return Err(invalid("hereditary basis metric is unrepresentable"));
            }
            basis_error = basis_error.max(mass_error).max(stiffness_error);
        }
    }
    if basis_error > controls.maximum_basis_error || viscosity.iter().any(|x| !x.is_finite()) {
        return Err(invalid(
            "hereditary equilibrium basis or Kelvin projection failed admission",
        ));
    }
    let sources = projected
        .material_points()
        .iter()
        .flat_map(|p| p.branches())
        .map(|b| (b.relaxation_time_s(), b.rows()))
        .collect::<Vec<_>>();
    let branches = memory::compact(&sources, n, controls.maximum_factor_scalar_terms, || {
        cx.checkpoint().is_err()
    })
    .map_err(|error| match error {
        memory::CompressionError::Invalid(what) => invalid(what),
        memory::CompressionError::Cancelled => crate::linear3::TetElasticError::Cancelled.into(),
    })?;
    let mut tangent = vec![0.0; n * n];
    for i in 0..n {
        tangent[i * n + i] = modes[i].lambda;
    }
    for branch in &branches {
        poll(cx)?;
        for row in branch.rows() {
            for i in 0..n {
                for j in 0..n {
                    tangent[i * n + j] += (row[i] / root) * (row[j] / root);
                }
            }
        }
    }
    let instantaneous = tangent
        .chunks_exact(n)
        .map(|row| row.iter().map(|x| x.abs()).sum::<f64>())
        .fold(0.0_f64, f64::max)
        .sqrt();
    if tangent.iter().any(|x| !x.is_finite())
        || !instantaneous.is_finite()
        || instantaneous > band.1
    {
        return Err(invalid(
            "hereditary instantaneous modes exceed the original material band",
        ));
    }
    let mut strain_bound = vec![0.0_f64; n];
    for point in projected.material_points() {
        poll(cx)?;
        // Even an unbounded source keeps its original finite squared-norm
        // arithmetic. This numerical ceiling is below that overflow boundary.
        let ceiling = point
            .material()
            .maximum_strain_norm()
            .min(f64::MAX.sqrt() * 0.25);
        for (j, column) in point.strain_columns().iter().enumerate() {
            let mandel: [f64; 6] = core::array::from_fn(|i| {
                column[i]
                    * if i < 3 {
                        1.0
                    } else {
                        core::f64::consts::SQRT_2
                    }
            });
            let scale = mandel.iter().map(|x| x.abs()).fold(0.0_f64, f64::max);
            if scale != 0.0 {
                let norm = scale
                    * mandel
                        .iter()
                        .map(|x| (x / scale).powi(2))
                        .sum::<f64>()
                        .sqrt();
                let ratio = if norm.is_finite() {
                    norm / ceiling
                } else {
                    f64::INFINITY
                };
                strain_bound[j] = strain_bound[j].max(ratio);
            }
        }
    }
    let points = projected.material_points().to_vec();
    let point_elements = projected.point_elements().to_vec();
    drop(projected);
    let static_source = reduced
        .pressure_compliance
        .ok_or_else(|| invalid("hereditary valve lost its equilibrium compliance"))?;
    let axis = crate::reduce::unit(request.orifice.opening_axis)?;
    let free = &source.equilibrium().free_dofs;
    let mut pressure_projection = Vec::with_capacity(n);
    let mut opening_projection = Vec::with_capacity(n);
    let mut compliance = ModalStaticCompliance {
        exact_opening_force_m_per_n: static_source.exact_mean_force_m_per_n(),
        retained_opening_force_m_per_n: 0.0,
        exact_volume_m3_per_pa: static_source.exact_volume_m3_per_pa(),
        retained_volume_m3_per_pa: 0.0,
        exact_opening_m_per_pa: static_source.exact_opening_m_per_pa(),
        retained_opening_m_per_pa: 0.0,
    };
    for mode in &modes {
        poll(cx)?;
        let opening =
            crate::reduce::mean_face_translation(&mode.phi, free, request.orifice.face_nodes, axis);
        let pressure = crate::reduce::pressure_area(&mode.phi, free, nodal_load_per_pa_m2, 1.0);
        if !opening.is_finite() || !pressure.is_finite() {
            return Err(invalid(
                "hereditary valve has a nonfinite pressure or opening projection",
            ));
        }
        compliance.retained_opening_force_m_per_n += opening * (opening / mode.lambda);
        compliance.retained_volume_m3_per_pa += pressure * (pressure / mode.lambda);
        compliance.retained_opening_m_per_pa += opening * (pressure / mode.lambda);
        pressure_projection.push(pressure);
        opening_projection.push(opening);
    }
    if [
        compliance.retained_opening_force_m_per_n,
        compliance.retained_volume_m3_per_pa,
        compliance.retained_opening_m_per_pa,
    ]
    .iter()
    .any(|x| !x.is_finite())
    {
        return Err(invalid("hereditary valve compliance overflow"));
    }
    let mut hash = DomainHasher::new("org.frankensim.fs-solid.hereditary-modal-valve.v1");
    hash.update(reduced.card.identity.as_bytes());
    hash.update(&(request.nodes_m.len() as u64).to_le_bytes());
    hash.update(&(free.len() as u64).to_le_bytes());
    for &dof in free {
        hash.update(&(dof as u64).to_le_bytes());
    }
    hash.update(&(n as u64).to_le_bytes());
    for mode in &modes {
        hash.update(&mode.lambda.to_bits().to_le_bytes());
        for x in &mode.phi {
            hash.update(&x.to_bits().to_le_bytes());
        }
    }
    for x in viscosity
        .iter()
        .chain(&pressure_projection)
        .chain(&opening_projection)
        .chain(nodal_load_per_pa_m2.iter().flatten())
    {
        hash.update(&x.to_bits().to_le_bytes());
    }
    hash.update(&(points.len() as u64).to_le_bytes());
    for point in &points {
        hash.update(
            &point
                .material()
                .maximum_strain_norm()
                .to_bits()
                .to_le_bytes(),
        );
        hash.update(&point.volume_m3().to_bits().to_le_bytes());
        for x in point.strain_columns().iter().flatten() {
            hash.update(&x.to_bits().to_le_bytes());
        }
    }
    for branch in &branches {
        hash.update(&branch.relaxation_time_s().to_bits().to_le_bytes());
        hash.update(&(branch.rows().len() as u64).to_le_bytes());
        for x in branch.rows().iter().flatten() {
            hash.update(&x.to_bits().to_le_bytes());
        }
    }
    poll(cx)?;
    Ok(HereditaryModalValve {
        source,
        modes,
        points,
        point_elements,
        viscosity,
        branches,
        pressure_projection,
        opening_projection,
        nodes: request.nodes_m.len(),
        rest_opening_m: reduced.card.rest_gap_m,
        width_m: reduced.card.width_m,
        total_mass_kg: mass,
        source_id: request.source_id.to_string(),
        identity: hash.finalize(),
        compliance,
        basis_error,
        band,
        instantaneous,
        strain_bound,
    })
}
