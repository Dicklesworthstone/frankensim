//! Closed distributed 3D windings on the original source-bound Whitney field.
//!
//! A profile is J/terminal_current [1/m2], including turns and orientation in
//! the original mesh frame. The native vector owner admits each actual 1-A
//! field, its complete curl-kernel compatibility and current continuity. This
//! owner additionally requires closed current support, including on prescribed
//! exterior traces. The original exact Whitney load pairs A with J_unit.
//!
//! For fixed induced-linear material, L_ij=int(B_i.H_j dV) must agree with
//! int(A_j.J_i dV) for EVERY ordered pair. Only the small reciprocal terminal
//! matrix enters the existing mutual-inductance circuit assembler. Neither the
//! field solver, magnetic inverse, Kirchhoff equations nor time step is copied.
//! Full 3D end geometry is represented where meshed; any omitted series self
//! inductance is explicit. No saturation, remanence, thermal refresh, moving
//! winding, skin/proximity current or unbounded-exterior accuracy is inferred.

use super::{
    SpatialMagneticError, SpatialMagneticMaterial, SpatialVectorField, SpatialVectorMagnetostatics,
    VectorMagneticBoundary, VectorMagneticControls, apply, dot, poll,
};
use fs_phs::{
    StepRecord,
    circuit::{
        Branch, CircuitDae, CircuitError, CircuitGraph, MutualInductance,
        assemble_circuit_with_mutual_inductance,
    },
};

/// Bounded field-basis preparation; original field controls govern every solve.
#[derive(Clone, Copy, Debug)]
pub struct SpatialWindingConfig {
    /// Number of physical winding coordinates, from 1 through 16.
    pub max_windings: usize,
    /// Conservative numeric-vector count: windings*(3*edges+faces+13*cells).
    /// Includes authored profiles, exact source loads and native unit fields;
    /// original material/solver allocations and small terminal matrices are extra.
    pub max_basis_values: usize,
    /// Absolute source-pairing versus magnetic-energy agreement [H].
    pub pairing_tolerance_h: f64,
}

#[derive(Debug)]
pub enum SpatialWindingError {
    Magnetic(SpatialMagneticError),
    Circuit(CircuitError),
    Invalid(&'static str),
    Pairing {
        source: usize,
        receiver: usize,
        flux_h: f64,
        energy_h: f64,
    },
}
impl core::fmt::Display for SpatialWindingError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "spatial winding: {self:?}")
    }
}
impl std::error::Error for SpatialWindingError {}
impl From<SpatialMagneticError> for SpatialWindingError {
    fn from(e: SpatialMagneticError) -> Self {
        Self::Magnetic(e)
    }
}
impl From<CircuitError> for SpatialWindingError {
    fn from(e: CircuitError) -> Self {
        Self::Circuit(e)
    }
}
type Result<T> = core::result::Result<T, SpatialWindingError>;
fn invalid(message: &'static str) -> SpatialWindingError {
    SpatialWindingError::Invalid(message)
}
fn finite(value: f64) -> Result<f64> {
    Ok(super::finite(value)?)
}

fn pairing(a: &[f64], b: &[f64], cancelled: &mut impl FnMut() -> bool) -> Result<f64> {
    let mut value = 0.0;
    for (&a, &b) in a.iter().zip(b) {
        poll(cancelled)?;
        value = finite(a.mul_add(b, value))?;
    }
    Ok(value)
}

fn current_defect(
    problem: &SpatialVectorMagnetostatics,
    current: &[[f64; 3]],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<f64> {
    let mut maximum = 0.0_f64;
    for face in &problem.geometry.faces {
        poll(cancelled)?;
        let jump = if face.cells.len() == 2 {
            std::array::from_fn(|k| current[face.cells[0]][k] - current[face.cells[1]][k])
        } else {
            current[face.cells[0]]
        };
        maximum = maximum.max(finite(face.area_m2 * dot(jump, face.outward_normal))?.abs());
    }
    Ok(maximum)
}

/// Native unit fields and their source/energy-audited full 3D terminal matrix.
/// Positive terminal current follows each original authored profile; reversing
/// one profile reverses its mutual entries and reconstructed field contribution.
pub struct SpatialWindingBasis<'m> {
    problem: &'m SpatialVectorMagnetostatics,
    boundary: VectorMagneticBoundary,
    fields: Vec<SpatialVectorField<'m>>,
    loads_per_a: Vec<Vec<f64>>,
    active_inductance_h: Vec<f64>,
    inductance_h: Vec<f64>,
    external_inductance_h: Vec<f64>,
    maximum_pairing_defect_h: f64,
    controls: VectorMagneticControls,
}
impl<'m> SpatialWindingBasis<'m> {
    /// Prepare fixed induced-linear material and homogeneous tangential traces.
    /// Each profile must describe an entire closed distributed current; current
    /// cannot disappear through a prescribed boundary. The existing field owner
    /// still checks the complete relative harmonic/gauge compatibility.
    /// `external_inductance_h` names only explicitly omitted series self energy;
    /// zero is allowed and no leakage is invented for dependent windings.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        problem: &'m SpatialVectorMagnetostatics,
        profiles_per_m2: Vec<Vec<[f64; 3]>>,
        boundary: VectorMagneticBoundary,
        external_inductance_h: Vec<f64>,
        config: SpatialWindingConfig,
        controls: VectorMagneticControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self> {
        poll(cancelled)?;
        let count = profiles_per_m2.len();
        let mesh = problem.mesh();
        let values = mesh
            .edges
            .len()
            .checked_mul(3)
            .and_then(|v| v.checked_add(mesh.faces.len()))
            .and_then(|v| {
                mesh.tets
                    .len()
                    .checked_mul(13)
                    .and_then(|c| v.checked_add(c))
            })
            .and_then(|v| v.checked_mul(count));
        if count == 0
            || count > config.max_windings
            || config.max_windings > 16
            || values.is_none_or(|v| v > config.max_basis_values)
            || external_inductance_h.len() != count
            || external_inductance_h
                .iter()
                .any(|v| !v.is_finite() || *v < 0.0)
            || !config.pairing_tolerance_h.is_finite()
            || config.pairing_tolerance_h < 0.0
            || !controls.current_tolerance_a.is_finite()
            || controls.current_tolerance_a <= 0.0
            || profiles_per_m2.iter().any(|p| {
                p.len() != mesh.tets.len()
                    || p.iter().flatten().any(|v| !v.is_finite())
                    || p.iter().flatten().all(|&v| v == 0.0)
            })
            || boundary
                .potential
                .iter()
                .any(|face| face.edge_integrals_wb != [0.0; 3])
            || problem
                .materials()
                .iter()
                .any(|m| !matches!(m, SpatialMagneticMaterial::Linear(_)))
        {
            return Err(invalid(
                "bounded complete closed winding profiles, induced-linear material and homogeneous traces required",
            ));
        }
        let mut fields = Vec::with_capacity(count);
        let mut loads_per_a = Vec::with_capacity(count);
        for profile in &profiles_per_m2 {
            poll(cancelled)?;
            if current_defect(problem, profile, cancelled)? > controls.current_tolerance_a {
                return Err(invalid(
                    "distributed winding has an unresolved current discontinuity or exterior return",
                ));
            }
            let field = problem.solve(profile, &boundary, controls, cancelled)?;
            let mut load = vec![0.0; mesh.edges.len()];
            for (cell, vectors) in problem.cells.iter().enumerate() {
                poll(cancelled)?;
                for local in 0..6 {
                    let edge = vectors.edges[local];
                    load[edge] = finite(load[edge] + dot(profile[cell], vectors.integrals[local]))?;
                }
            }
            fields.push(field);
            loads_per_a.push(load);
        }
        let mut active_inductance_h = vec![0.0; count * count];
        let mut maximum_pairing_defect_h = 0.0_f64;
        for receiver in 0..count {
            for source in 0..=receiver {
                let mut energy = 0.0;
                for (cell, (&a, &b)) in fields[receiver]
                    .flux_density_t()
                    .iter()
                    .zip(fields[source].flux_density_t())
                    .enumerate()
                {
                    poll(cancelled)?;
                    energy = finite(
                        (problem.geometry.geometry.vol_signed[cell].abs()
                            / problem.permeability[cell])
                            .mul_add(dot(a, b), energy),
                    )?;
                }
                for (i, j) in [(receiver, source), (source, receiver)] {
                    let flux = pairing(&loads_per_a[i], fields[j].potential_wb(), cancelled)?;
                    let defect = finite(flux - energy)?.abs();
                    maximum_pairing_defect_h = maximum_pairing_defect_h.max(defect);
                    if defect > config.pairing_tolerance_h {
                        return Err(SpatialWindingError::Pairing {
                            source: j,
                            receiver: i,
                            flux_h: flux,
                            energy_h: energy,
                        });
                    }
                }
                active_inductance_h[receiver * count + source] = energy;
                active_inductance_h[source * count + receiver] = energy;
            }
        }
        let mut inductance_h = active_inductance_h.clone();
        for (i, &external) in external_inductance_h.iter().enumerate() {
            inductance_h[i * count + i] = finite(inductance_h[i * count + i] + external)?;
            if inductance_h[i * count + i] <= 0.0 {
                return Err(invalid(
                    "winding has no represented positive magnetic energy",
                ));
            }
        }
        poll(cancelled)?;
        Ok(Self {
            problem,
            boundary,
            fields,
            loads_per_a,
            active_inductance_h,
            inductance_h,
            external_inductance_h,
            maximum_pairing_defect_h,
            controls,
        })
    }
    pub const fn problem(&self) -> &'m SpatialVectorMagnetostatics {
        self.problem
    }
    pub const fn boundary(&self) -> &VectorMagneticBoundary {
        &self.boundary
    }
    pub fn winding_count(&self) -> usize {
        self.fields.len()
    }
    /// Original admitted 1-A field, including its complete source and solve reports.
    pub fn unit_field(&self, winding: usize) -> Option<&SpatialVectorField<'m>> {
        self.fields.get(winding)
    }
    pub fn active_inductance_h(&self) -> &[f64] {
        &self.active_inductance_h
    }
    pub fn inductance_h(&self) -> &[f64] {
        &self.inductance_h
    }
    pub fn external_inductance_h(&self) -> &[f64] {
        &self.external_inductance_h
    }
    pub const fn maximum_pairing_defect_h(&self) -> f64 {
        self.maximum_pairing_defect_h
    }

    fn admit_currents(&self, currents: &[f64]) -> Result<()> {
        if currents.len() != self.winding_count() || currents.iter().any(|v| !v.is_finite()) {
            return Err(invalid("complete finite 3D winding currents required"));
        }
        Ok(())
    }
    pub fn flux_linkages_wb(&self, currents: &[f64]) -> Result<Vec<f64>> {
        self.admit_currents(currents)?;
        self.inductance_h
            .chunks(self.winding_count())
            .map(|row| {
                row.iter()
                    .zip(currents)
                    .try_fold(0.0, |sum, (&l, &i)| finite(l.mul_add(i, sum)))
            })
            .collect()
    }
    /// Exact original Whitney source functional on a complete edge potential.
    /// This diagnoses gauge invariance; it does not admit an arbitrary potential
    /// as a physical field and excludes declared external series inductance.
    pub fn pair_potential_wb(
        &self,
        potential_wb: &[f64],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Vec<f64>> {
        poll(cancelled)?;
        if potential_wb.len() != self.problem.mesh().edges.len()
            || potential_wb.iter().any(|v| !v.is_finite())
        {
            return Err(invalid(
                "complete finite original Whitney potential required",
            ));
        }
        self.loads_per_a
            .iter()
            .map(|load| pairing(load, potential_wb, cancelled))
            .collect()
    }

    /// Optional actual-amplitude reconstruction; no field solve occurs here.
    /// Original curl/constitutive and source operators recover B/H, current and
    /// storage. Original absolute field/current/work allowances are rechecked;
    /// multiplying unit fields never silently multiplies their error allowance.
    pub fn field_for_currents(
        &self,
        currents: &[f64],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<SpatialWindingField<'m>> {
        poll(cancelled)?;
        self.admit_currents(currents)?;
        let problem = self.problem;
        let mut potential_wb = vec![0.0; problem.mesh().edges.len()];
        let mut loads = potential_wb.clone();
        let mut current_density_a_m2 = vec![[0.0; 3]; problem.mesh().tets.len()];
        let mut compatibility = 0.0;
        for (winding, (&current, field)) in currents.iter().zip(&self.fields).enumerate() {
            compatibility = finite(
                current
                    .abs()
                    .mul_add(field.source_compatibility_defect_a, compatibility),
            )?;
            for (edge, value) in potential_wb.iter_mut().enumerate() {
                poll(cancelled)?;
                *value = finite(current.mul_add(field.potential_wb()[edge], *value))?;
                loads[edge] =
                    finite(current.mul_add(self.loads_per_a[winding][edge], loads[edge]))?;
            }
            for (out, unit) in current_density_a_m2
                .iter_mut()
                .zip(field.impressed_current_density_a_m2())
            {
                poll(cancelled)?;
                for k in 0..3 {
                    out[k] = finite(current.mul_add(unit[k], out[k]))?;
                }
            }
        }
        let equation = apply(&problem.stiffness, &potential_wb, cancelled)?;
        let mut maximum_residual_a = 0.0_f64;
        for (edge, (&actual, &source)) in equation.iter().zip(&loads).enumerate() {
            poll(cancelled)?;
            if !self.fields[0].prescribed_edges()[edge] {
                maximum_residual_a = maximum_residual_a.max(finite(actual - source)?.abs());
            }
        }
        let current_defect = current_defect(problem, &current_density_a_m2, cancelled)?;
        if maximum_residual_a > self.controls.current_tolerance_a
            || current_defect > self.controls.current_tolerance_a
            || compatibility > self.controls.current_tolerance_a
        {
            return Err(invalid(
                "reconstructed winding current or field exceeds its original native allowance",
            ));
        }
        let face_flux_wb = apply(&problem.curl, &potential_wb, cancelled)?;
        let mut maximum_divergence_flux_wb = 0.0_f64;
        for row in problem.mesh().d2().rows {
            poll(cancelled)?;
            let flux = finite(
                row.iter()
                    .map(|&(face, sign)| f64::from(sign) * face_flux_wb[face])
                    .sum(),
            )?;
            maximum_divergence_flux_wb = maximum_divergence_flux_wb.max(flux.abs());
        }
        if maximum_divergence_flux_wb > self.controls.flux_tolerance_wb {
            return Err(invalid(
                "reconstructed magnetic flux exceeds its original div-curl allowance",
            ));
        }
        let mut flux_density_t = Vec::with_capacity(problem.cells.len());
        let mut field_strength_a_m = Vec::with_capacity(problem.cells.len());
        let mut cell_energy_j = Vec::with_capacity(problem.cells.len());
        let mut field_energy_j = 0.0;
        for (cell, vectors) in problem.cells.iter().enumerate() {
            poll(cancelled)?;
            let mut b = [0.0; 3];
            for local in 0..6 {
                for k in 0..3 {
                    b[k] = finite(
                        potential_wb[vectors.edges[local]].mul_add(vectors.curls[local][k], b[k]),
                    )?;
                }
            }
            let h = b.map(|v| v / problem.permeability[cell]);
            for value in h {
                finite(value)?;
            }
            let energy =
                finite(0.5 * problem.geometry.geometry.vol_signed[cell].abs() * dot(b, h))?;
            field_energy_j = finite(field_energy_j + energy)?;
            flux_density_t.push(b);
            field_strength_a_m.push(h);
            cell_energy_j.push(energy);
        }
        let mut maximum_recovered_flux_defect_wb = 0.0_f64;
        for (face_id, face) in problem.geometry.faces.iter().enumerate() {
            poll(cancelled)?;
            let [a, b, c] =
                problem.mesh().faces[face_id].map(|v| problem.positions_m()[v as usize]);
            let raw = super::cross(super::super::sub(b, a), super::super::sub(c, a));
            let sign = if dot(raw, face.outward_normal) > 0.0 {
                1.0
            } else {
                -1.0
            };
            for &cell in &face.cells {
                let flux = finite(face.area_m2 * dot(flux_density_t[cell], face.outward_normal))?;
                maximum_recovered_flux_defect_wb = maximum_recovered_flux_defect_wb
                    .max(finite(flux - sign * face_flux_wb[face_id])?.abs());
            }
        }
        if maximum_recovered_flux_defect_wb > self.controls.flux_tolerance_wb {
            return Err(invalid(
                "reconstructed cell B exceeds its original conforming face-flux allowance",
            ));
        }
        let flux_linkages_wb = self.flux_linkages_wb(currents)?;
        let mut external_energy_j = 0.0;
        for (&l, &i) in self.external_inductance_h.iter().zip(currents) {
            external_energy_j = finite((0.5 * l * i).mul_add(i, external_energy_j))?;
        }
        let source_pairing_j = pairing(&potential_wb, &loads, cancelled)?;
        let stored_energy_j = finite(field_energy_j + external_energy_j)?;
        let terminal_energy_j = finite(0.5 * pairing(currents, &flux_linkages_wb, cancelled)?)?;
        let energy_residual_j = finite(stored_energy_j - terminal_energy_j)?;
        if energy_residual_j.abs() > self.controls.energy_tolerance_j
            || finite(2.0 * field_energy_j - source_pairing_j)?.abs()
                > self.controls.energy_tolerance_j
        {
            return Err(invalid(
                "reconstructed 3D field, source pairing and native terminal energy disagree",
            ));
        }
        poll(cancelled)?;
        Ok(SpatialWindingField {
            problem,
            potential_wb,
            face_flux_wb,
            flux_density_t,
            field_strength_a_m,
            current_density_a_m2,
            cell_energy_j,
            flux_linkages_wb,
            field_energy_j,
            external_energy_j,
            stored_energy_j,
            maximum_residual_a,
            maximum_current_flux_defect_a: current_defect,
            maximum_divergence_flux_wb,
            maximum_recovered_flux_defect_wb,
            energy_residual_j,
        })
    }

    /// Replace selected original inductor placeholders with this complete L.
    /// Other inductors, capacitors, resistors, transformers and source port
    /// indices remain native. The original assembler admits the FULL matrix;
    /// dependent modes require explicit physical leakage or an ideal transformer.
    pub fn into_circuit(
        self,
        mut graph: CircuitGraph,
        branch_indices: Vec<usize>,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<SpatialWindingCircuit<'m>> {
        poll(cancelled)?;
        let n = self.winding_count();
        if graph.node_count > 32
            || graph.branches.len() > 64
            || graph.transformers.len() > 16
            || branch_indices.len() != n
        {
            return Err(invalid(
                "bounded complete original 3D winding circuit required",
            ));
        }
        let inductors: Vec<_> = graph
            .branches
            .iter()
            .enumerate()
            .filter_map(|(i, (_, _, branch))| {
                matches!(branch, Branch::Inductor { .. }).then_some(i)
            })
            .collect();
        let mut occupied = vec![false; graph.branches.len()];
        let mut winding_slots = Vec::with_capacity(n);
        for (winding, &branch) in branch_indices.iter().enumerate() {
            let slot = inductors
                .iter()
                .position(|&i| i == branch)
                .ok_or_else(|| invalid("each 3D winding must bind an original inductor branch"))?;
            if occupied[branch] {
                return Err(invalid("duplicate 3D winding circuit branch"));
            }
            occupied[branch] = true;
            winding_slots.push(slot);
            graph.branches[branch].2 = Branch::Inductor {
                henries: self.inductance_h[winding * n + winding],
            };
        }
        for &(_, _, branch) in &graph.branches {
            match branch {
                Branch::VoltageSource { port } | Branch::CurrentSource { port } if port >= 32 => {
                    return Err(invalid("3D winding circuit source port exceeds capacity"));
                }
                _ => {}
            }
        }
        let mut mutual = Vec::with_capacity(n * (n - 1) / 2);
        for i in 0..n {
            for j in 0..i {
                mutual.push(MutualInductance {
                    branches: (branch_indices[i], branch_indices[j]),
                    henries: self.inductance_h[i * n + j],
                });
            }
        }
        let dae = assemble_circuit_with_mutual_inductance(&graph, &mutual)?;
        if dae.system.state_dim() > 128 {
            return Err(invalid("3D winding native circuit state exceeds capacity"));
        }
        let winding_coordinates = winding_slots.iter().map(|&i| dae.flux_index[i]).collect();
        let remaining_coordinates = inductors
            .iter()
            .enumerate()
            .filter_map(|(i, &branch)| (!occupied[branch]).then_some(dae.flux_index[i]))
            .collect();
        poll(cancelled)?;
        Ok(SpatialWindingCircuit {
            basis: self,
            dae,
            graph,
            branch_indices,
            winding_slots,
            winding_coordinates,
            remaining_coordinates,
        })
    }
}

/// Reconstructed physical field, independent of optional circuit storage.
/// `stored_energy_j` includes the explicitly declared external inductances;
/// `field_energy_j` and cell energy include only the original meshed domain.
pub struct SpatialWindingField<'m> {
    pub problem: &'m SpatialVectorMagnetostatics,
    pub potential_wb: Vec<f64>,
    pub face_flux_wb: Vec<f64>,
    pub flux_density_t: Vec<[f64; 3]>,
    pub field_strength_a_m: Vec<[f64; 3]>,
    pub current_density_a_m2: Vec<[f64; 3]>,
    pub cell_energy_j: Vec<f64>,
    pub flux_linkages_wb: Vec<f64>,
    pub field_energy_j: f64,
    pub external_energy_j: f64,
    pub stored_energy_j: f64,
    pub maximum_residual_a: f64,
    pub maximum_current_flux_defect_a: f64,
    pub maximum_divergence_flux_wb: f64,
    pub maximum_recovered_flux_defect_wb: f64,
    pub energy_residual_j: f64,
}

/// Original native circuit plus the actual 3D source/field basis. Native flux
/// coordinates already contain ALL represented winding energy; reconstructed
/// field energy must never be added to the same circuit Hamiltonian again.
pub struct SpatialWindingCircuit<'m> {
    basis: SpatialWindingBasis<'m>,
    dae: CircuitDae,
    graph: CircuitGraph,
    branch_indices: Vec<usize>,
    winding_slots: Vec<usize>,
    winding_coordinates: Vec<usize>,
    remaining_coordinates: Vec<usize>,
}
impl<'m> SpatialWindingCircuit<'m> {
    pub const fn basis(&self) -> &SpatialWindingBasis<'m> {
        &self.basis
    }
    pub const fn dae(&self) -> &CircuitDae {
        &self.dae
    }
    pub const fn graph(&self) -> &CircuitGraph {
        &self.graph
    }
    pub fn branch_indices(&self) -> &[usize] {
        &self.branch_indices
    }
    pub fn flux_coordinates(&self) -> &[usize] {
        &self.winding_coordinates
    }

    /// Original consistent initialization preserves selected total linkages,
    /// other inductor fluxes and capacitor charges. Other storage uses original
    /// branch order; winding currents use the caller's explicit winding order.
    pub fn consistent_initial_state(
        &self,
        currents_a: &[f64],
        remaining_fluxes_wb: &[f64],
        charges_c: &[f64],
        inputs: &[f64],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Vec<f64>> {
        poll(cancelled)?;
        let linkages = self.basis.flux_linkages_wb(currents_a)?;
        if remaining_fluxes_wb.len() != self.remaining_coordinates.len()
            || charges_c.len() != self.dae.charge_index.len()
            || inputs.len() != self.dae.ports
            || remaining_fluxes_wb
                .iter()
                .chain(charges_c)
                .chain(inputs)
                .any(|v| !v.is_finite())
        {
            return Err(invalid(
                "complete finite original circuit storage and source inputs required",
            ));
        }
        let mut state = vec![0.0; self.dae.system.state_dim()];
        for (&index, &value) in self
            .winding_coordinates
            .iter()
            .zip(&linkages)
            .chain(self.remaining_coordinates.iter().zip(remaining_fluxes_wb))
            .chain(self.dae.charge_index.iter().zip(charges_c))
        {
            state[index] = value;
        }
        let result = self.dae.consistent_initial_state(&state, inputs)?;
        poll(cancelled)?;
        Ok(result)
    }
    pub fn currents_a(&self, state: &[f64]) -> Result<Vec<f64>> {
        let currents = self.dae.inductor_currents(state)?;
        Ok(self.winding_slots.iter().map(|&i| currents[i]).collect())
    }
    /// Native audited dynamics; fixed linear field extraction is a preparation
    /// cost, not a new field solve/factorization inside every circuit step.
    pub fn step_audited(
        &self,
        state: &[f64],
        inputs: &[f64],
        dt: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(StepRecord, f64)> {
        poll(cancelled)?;
        self.currents_a(state)?;
        let next = self.dae.step_audited(state, inputs, dt)?;
        poll(cancelled)?;
        Ok(next)
    }
    pub fn field_for_state(
        &self,
        state: &[f64],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<SpatialWindingField<'m>> {
        self.basis
            .field_for_currents(&self.currents_a(state)?, cancelled)
    }
}
