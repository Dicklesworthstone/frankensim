//! Planar magnetic vector potential with source-resolved linear material laws.
//!
//! A_z is in Wb/m, J_z in A/m2, coordinates/thickness in metres. Prescribed
//! boundary vertices fix A_z; unprescribed boundary edges carry the natural
//! zero tangential-H condition. Every connected mesh component needs a
//! complete prescribed exterior edge. The solution is a discrete P1 Galerkin estimate,
//! with a recomputed algebraic residual, not an error enclosure.

/// Retained planar induction and material eddy-current heat on this same mesh.
pub mod eddy;
/// Reversible nonlinear H(B) on the same mesh, Hodge and checked solver owners.
pub mod nonlinear;
/// Source-backed, bounded reversible permanent magnets on the same field owner.
pub mod permanent;

use fs_feec::{
    assembly::{incidence_to_csr, stiffness},
    weighted::{WeightedAssemblyLimits, WeightedError},
    weighted2::{CellWeight2, weighted_mass_matrix_2d},
};
use fs_material::{
    field_storage::{
        FieldStorageError, linear_magnetic_permeability, permanent::RecoilOperatingPoint,
    },
    state_point::ResolvedMaterialStatePoint,
};
use fs_qty::Dims;
use fs_rep_mesh::{Metric2, TriComplex2};
use fs_solver::{CheckedCgConfig, CheckedCgError, CsrOp, SolveReport, checked_cg};
use fs_sparse::{Csr, precond::IdentityPrecond};
use permanent::{MagneticEnergyReference, PlanarMagneticMaterial, PlanarRecoilMaterial};

/// Source, geometry, work, cancellation or actual residual refusal.
#[derive(Debug)]
pub enum MagneticError {
    /// Shared linear lossless permeability admission.
    Material(FieldStorageError),
    /// Shared bounded material Hodge assembly.
    Assembly(WeightedError),
    /// Invalid or unsupported physical/numerical input.
    Invalid(&'static str),
    /// Interrupted work publishes no partially accepted field.
    Cancelled,
    /// Checked CG exhausted its work budget or failed its true residual target.
    Nonconvergence(SolveReport),
}
impl core::fmt::Display for MagneticError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "planar magnetostatics: {self:?}")
    }
}
impl std::error::Error for MagneticError {}

pub(crate) fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), MagneticError> {
    if cancelled() {
        Err(MagneticError::Cancelled)
    } else {
        Ok(())
    }
}
fn finite(value: f64) -> Result<f64, MagneticError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(MagneticError::Invalid(
            "field arithmetic is not representable",
        ))
    }
}

/// Prepared immutable mesh, complete source snapshots and shared FEEC stiffness.
/// A different material/state creates a different problem; no coefficient cache
/// can survive a source replacement. Temperature is fixed within each solve.
pub struct PlanarMagnetostatics {
    mesh: TriComplex2,
    materials: Vec<ResolvedMaterialStatePoint>,
    permeability_h_m: Vec<f64>,
    recoil_materials: Vec<Option<PlanarRecoilMaterial>>,
    stiffness: Csr,
    thickness_m: f64,
}
impl PlanarMagnetostatics {
    /// Admit one material state per triangle, then assemble d0^T M(1/mu) d0.
    /// Absolute SI permeability is required. The shared circuit admission
    /// refuses field/frequency-dependent or bounded operating-field claims.
    /// Axisymmetric A_phi needs a different operator and is refused here.
    pub fn new(
        mesh: TriComplex2,
        materials: Vec<ResolvedMaterialStatePoint>,
        limits: WeightedAssemblyLimits,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, MagneticError> {
        Self::from_material_models(
            mesh,
            materials.into_iter().map(PlanarMagneticMaterial::Linear),
            limits,
            cancelled,
        )
    }

    /// Assemble the same operator from explicit per-cell constitutive choices.
    /// Recoil cells use their source-declared absolute recoil slope and retain
    /// Br, specimen direction and operating limits for every solve/recovery.
    /// Current density and remanence are independent sources. No coercivity or
    /// generic permeability is promoted to a reversible recoil coefficient.
    pub fn new_with_material_models(
        mesh: TriComplex2,
        models: Vec<PlanarMagneticMaterial>,
        limits: WeightedAssemblyLimits,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, MagneticError> {
        Self::from_material_models(mesh, models.into_iter(), limits, cancelled)
    }

    fn from_material_models(
        mesh: TriComplex2,
        models: impl ExactSizeIterator<Item = PlanarMagneticMaterial>,
        limits: WeightedAssemblyLimits,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, MagneticError> {
        poll(cancelled)?;
        let thickness_m = match mesh.metric() {
            Metric2::Planar { thickness } => thickness,
            Metric2::Axisymmetric { .. } => {
                return Err(MagneticError::Invalid("planar A_z requires planar metric"));
            }
        };
        if mesh.vertices().len() > limits.max_dofs
            || mesh.faces().len() > limits.max_cells
            || mesh.edges().len() > limits.max_dofs
            || mesh
                .faces()
                .len()
                .checked_mul(9)
                .is_none_or(|n| n > limits.max_triplets)
        {
            return Err(MagneticError::Invalid(
                "mesh exceeds declared assembly limits",
            ));
        }
        if models.len() != mesh.faces().len() {
            return Err(MagneticError::Invalid(
                "one source material is required per triangle",
            ));
        }
        let mut materials = Vec::with_capacity(models.len());
        let mut recoil_materials = Vec::with_capacity(models.len());
        let mut permeability_h_m = Vec::with_capacity(models.len());
        let mut weights = Vec::with_capacity(models.len());
        for model in models {
            poll(cancelled)?;
            let (material, recoil, mu) = match model {
                PlanarMagneticMaterial::Linear(material) => {
                    let mu =
                        linear_magnetic_permeability(&material).map_err(MagneticError::Material)?;
                    (material, None, mu)
                }
                PlanarMagneticMaterial::Recoil(recoil) => {
                    let mu = recoil.material().permeability_h_m();
                    (recoil.material().material().clone(), Some(recoil), mu)
                }
            };
            let nu = mu.recip();
            if !nu.is_finite() || nu <= 0.0 {
                return Err(MagneticError::Invalid(
                    "reluctivity must be representably positive",
                ));
            }
            permeability_h_m.push(mu);
            materials.push(material);
            recoil_materials.push(recoil);
            weights.push(CellWeight2::Scalar(nu));
        }
        let mass = weighted_mass_matrix_2d(
            &mesh,
            1,
            &weights,
            Dims([-1, -1, 2, 0, 2, 0]),
            limits,
            cancelled,
        )
        .map_err(|error| match error {
            WeightedError::Cancelled => MagneticError::Cancelled,
            error => MagneticError::Assembly(error),
        })?;
        poll(cancelled)?;
        let d0 = incidence_to_csr(&mesh.d0());
        let stiffness = stiffness(&d0, mass.matrix());
        poll(cancelled)?;
        for r in 0..stiffness.nrows() {
            poll(cancelled)?;
            if !stiffness.row(r).1.iter().all(|v| v.is_finite()) {
                return Err(MagneticError::Invalid(
                    "unrepresentable assembled stiffness",
                ));
            }
        }
        Ok(Self {
            mesh,
            materials,
            permeability_h_m,
            recoil_materials,
            stiffness,
            thickness_m,
        })
    }
    /// Exact geometry and durable feature identities used for the field.
    pub fn mesh(&self) -> &TriComplex2 {
        &self.mesh
    }
    /// Complete original per-cell query receipts, including fixed temperature.
    pub fn materials(&self) -> &[ResolvedMaterialStatePoint] {
        &self.materials
    }
    /// Absolute per-cell permeability [H/m].
    pub fn permeability_h_m(&self) -> &[f64] {
        &self.permeability_h_m
    }
    /// Original source law and orientation for every magnet cell; None marks
    /// induced-linear cells. Storage order is the original triangle order.
    pub fn recoil_materials(&self) -> &[Option<PlanarRecoilMaterial>] {
        &self.recoil_materials
    }
    /// Whether this owner carries an affine remanence source. Consumers using
    /// homogeneous inductance/induction energy must explicitly account for it.
    pub fn has_permanent_magnets(&self) -> bool {
        self.recoil_materials.iter().any(Option::is_some)
    }

    /// Solve for cellwise impressed winding current density, retained material
    /// remanence and prescribed A_z. Every magnet's resulting H must remain in
    /// its complete source-backed recoil envelope before a field is returned.
    /// Boundary entries are unique (vertex, Wb/m) pairs and must form complete
    /// exterior edges: isolated point constraints are refused. Unspecified
    /// exterior edges have zero natural tangential H; no open-boundary/infinite-domain
    /// approximation is inferred. Loads are integrated with exact P1 moments.
    /// Returned fields borrow this owner, retaining source/mesh correspondence.
    pub fn solve(
        &self,
        current_density_a_m2: &[f64],
        boundary_potential_wb_m: &[(usize, f64)],
        config: CheckedCgConfig,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<PlanarField<'_>, MagneticError> {
        self.solve_impl(
            current_density_a_m2,
            boundary_potential_wb_m,
            config,
            true,
            cancelled,
        )
    }

    // A derivative with respect to impressed current, not a realizable magnet
    // operating state. It shares the exact operator but has no remanence RHS.
    pub(crate) fn solve_response(
        &self,
        current_density_a_m2: &[f64],
        boundary_potential_wb_m: &[(usize, f64)],
        config: CheckedCgConfig,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<PlanarField<'_>, MagneticError> {
        self.solve_impl(
            current_density_a_m2,
            boundary_potential_wb_m,
            config,
            false,
            cancelled,
        )
    }

    fn solve_impl(
        &self,
        current_density_a_m2: &[f64],
        boundary_potential_wb_m: &[(usize, f64)],
        config: CheckedCgConfig,
        include_remanence: bool,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<PlanarField<'_>, MagneticError> {
        poll(cancelled)?;
        if current_density_a_m2.len() != self.mesh.faces().len()
            || !current_density_a_m2.iter().all(|v| v.is_finite())
            || !(config.tolerance.is_finite() && config.tolerance > 0.0 && config.tolerance < 1.0)
        {
            return Err(MagneticError::Invalid(
                "current field or residual tolerance is invalid",
            ));
        }
        let n = self.mesh.vertices().len();
        let prescribed = self.admit_boundary(boundary_potential_wb_m, cancelled)?;
        let mut rhs = vec![0.0; n];
        for (cell, face) in self.mesh.faces().iter().enumerate() {
            poll(cancelled)?;
            let volume = self
                .mesh
                .face_measure(cell)
                .ok_or(MagneticError::Invalid("missing triangle measure"))?;
            let load = finite(current_density_a_m2[cell] * (volume / 3.0))?;
            for &vertex in face {
                rhs[vertex as usize] = finite(rhs[vertex as usize] + load)?;
            }
            if let Some(recoil) = self.recoil_materials[cell]
                .as_ref()
                .filter(|_| include_remanence)
            {
                let br = recoil.remanence_t();
                let source_h = [
                    finite(br[0] / self.permeability_h_m[cell])?,
                    finite(br[1] / self.permeability_h_m[cell])?,
                ];
                let gradients = gradients(face.map(|v| self.mesh.vertices()[v as usize]))?;
                // Exact weak load ∫(Br/mu)·(d_y phi,-d_x phi) dV. Keeping the
                // whole P1 moment preserves both interface and exterior terms:
                // a constant Br has zero volume curl but a nonzero boundary load.
                for (local, &vertex) in face.iter().enumerate() {
                    let moment = finite(
                        source_h[0]
                            .mul_add(gradients[local][1], -source_h[1] * gradients[local][0]),
                    )?;
                    rhs[vertex as usize] = finite(volume.mul_add(moment, rhs[vertex as usize]))?;
                }
            }
        }
        let free: Vec<_> = (0..n).filter(|&v| prescribed[v].is_none()).collect();
        let mut index = vec![usize::MAX; n];
        for (i, &v) in free.iter().enumerate() {
            index[v] = i;
        }
        let mut row_ptr = vec![0];
        let mut columns = Vec::new();
        let mut values = Vec::new();
        let mut reduced_rhs = Vec::with_capacity(free.len());
        for &v in &free {
            poll(cancelled)?;
            let mut load = rhs[v];
            let (cols, vals) = self.stiffness.row(v);
            for (&w, &value) in cols.iter().zip(vals) {
                if let Some(potential) = prescribed[w] {
                    load = finite((-value).mul_add(potential, load))?;
                } else {
                    columns.push(index[w]);
                    values.push(value);
                }
            }
            row_ptr.push(columns.len());
            reduced_rhs.push(load);
        }
        let matrix = Csr::from_parts(free.len(), free.len(), row_ptr, columns, values);
        let mut potential: Vec<_> = prescribed
            .iter()
            .map(|value| value.unwrap_or(0.0))
            .collect();
        let report = if free.is_empty() {
            None
        } else {
            let op = CsrOp::symmetric(matrix);
            let result = checked_cg(&op, &IdentityPrecond, &reduced_rhs, config, |_| {
                poll(cancelled)
            })
            .map_err(|error| match error {
                CheckedCgError::InvalidInput(why) => MagneticError::Invalid(why),
                CheckedCgError::Interrupted(error) => error,
            })?;
            if !result.report.converged_euclidean() {
                return Err(MagneticError::Nonconvergence(result.report));
            }
            for (&v, &value) in free.iter().zip(&result.x) {
                potential[v] = value;
            }
            Some(result.report)
        };
        let mut field = self.recover_impl(potential, report, include_remanence, cancelled)?;
        field.impressed_current_density_a_m2 = Some(current_density_a_m2.to_vec());
        poll(cancelled)?;
        Ok(field)
    }

    fn admit_boundary(
        &self,
        boundary_potential_wb_m: &[(usize, f64)],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Vec<Option<f64>>, MagneticError> {
        Self::admit_mesh_boundary(&self.mesh, boundary_potential_wb_m, cancelled)
    }

    pub(crate) fn admit_mesh_boundary(
        mesh: &TriComplex2,
        boundary_potential_wb_m: &[(usize, f64)],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Vec<Option<f64>>, MagneticError> {
        let prescribed = Self::admit_mesh_traces(mesh, boundary_potential_wb_m, cancelled)?;
        Self::admit_mesh_gauges(mesh, &prescribed, cancelled)?;
        Ok(prescribed)
    }

    // Physical trace admission is shared with conduction interfaces, whose
    // positive surface conductance can connect otherwise separate mesh pieces.
    pub(crate) fn admit_mesh_traces(
        mesh: &TriComplex2,
        boundary_potential_wb_m: &[(usize, f64)],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Vec<Option<f64>>, MagneticError> {
        let n = mesh.vertices().len();
        let trace = mesh
            .boundary_trace()
            .map_err(|_| MagneticError::Invalid("mesh boundary trace is invalid"))?;
        let mut on_boundary = vec![false; n];
        for &vertex in &trace.vertices {
            on_boundary[vertex as usize] = true;
        }
        let mut prescribed = vec![None; n];
        for &(vertex, value) in boundary_potential_wb_m {
            poll(cancelled)?;
            if vertex >= n
                || !on_boundary[vertex]
                || !value.is_finite()
                || prescribed[vertex].is_some()
            {
                return Err(MagneticError::Invalid(
                    "boundary vertex/value is invalid or repeated",
                ));
            }
            prescribed[vertex] = Some(value);
        }
        // A point pin alone fixes a matrix gauge but cannot supply the
        // positive-length Dirichlet boundary of this physical model.
        let mut prescribed_edge = vec![false; n];
        for edge in &trace.edges {
            poll(cancelled)?;
            let [a, b] = edge.oriented_vertices.map(|v| v as usize);
            if prescribed[a].is_some() && prescribed[b].is_some() {
                prescribed_edge[a] = true;
                prescribed_edge[b] = true;
            }
        }
        if prescribed
            .iter()
            .zip(&prescribed_edge)
            .any(|(value, &on_edge)| value.is_some() && !on_edge)
        {
            return Err(MagneticError::Invalid(
                "prescribed potentials must cover exterior edges, not isolated points",
            ));
        }
        Ok(prescribed)
    }

    pub(crate) fn recover(
        &self,
        potential: Vec<f64>,
        report: Option<SolveReport>,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<PlanarField<'_>, MagneticError> {
        self.recover_impl(potential, report, true, cancelled)
    }

    /// Reconstruct an actual static field when the winding owner also retains
    /// the corresponding source current. Generic induction recovery deliberately
    /// does not establish this provenance for static force admission.
    pub(crate) fn recover_static(
        &self,
        potential: Vec<f64>,
        current_density_a_m2: Vec<f64>,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<PlanarField<'_>, MagneticError> {
        if current_density_a_m2.len() != self.mesh.faces().len()
            || current_density_a_m2.iter().any(|v| !v.is_finite())
        {
            return Err(MagneticError::Invalid(
                "static field recovery requires complete finite impressed currents",
            ));
        }
        let mut field = self.recover(potential, None, cancelled)?;
        field.impressed_current_density_a_m2 = Some(current_density_a_m2);
        poll(cancelled)?;
        Ok(field)
    }

    fn recover_impl(
        &self,
        potential: Vec<f64>,
        report: Option<SolveReport>,
        include_remanence: bool,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<PlanarField<'_>, MagneticError> {
        if potential.len() != self.mesh.vertices().len() || potential.iter().any(|v| !v.is_finite())
        {
            return Err(MagneticError::Invalid(
                "field recovery requires complete finite nodal potentials",
            ));
        }
        let mut magnetic_flux_density_t = Vec::with_capacity(self.mesh.faces().len());
        let mut magnetic_field_a_m = Vec::with_capacity(self.mesh.faces().len());
        let mut cell_energy_j = Vec::with_capacity(self.mesh.faces().len());
        let mut recoil_operating_points = Vec::with_capacity(self.mesh.faces().len());
        let mut energy_j = 0.0;
        for (cell, face) in self.mesh.faces().iter().enumerate() {
            poll(cancelled)?;
            let points = face.map(|v| self.mesh.vertices()[v as usize]);
            let gradients = gradients(points)?;
            let mut grad = [0.0; 2];
            // Subtract one local potential before differentiation: a constant
            // gauge must not create a field through cancellation of large values.
            for local in 1..3 {
                let delta = finite(potential[face[local] as usize] - potential[face[0] as usize])?;
                for a in 0..2 {
                    grad[a] = finite(delta.mul_add(gradients[local][a], grad[a]))?;
                }
            }
            let b = [grad[1], -grad[0]];
            let br = self.recoil_materials[cell]
                .as_ref()
                .filter(|_| include_remanence)
                .map_or([0.0; 2], PlanarRecoilMaterial::remanence_t);
            let induced_b = [finite(b[0] - br[0])?, finite(b[1] - br[1])?];
            let h = [
                finite(induced_b[0] / self.permeability_h_m[cell])?,
                finite(induced_b[1] / self.permeability_h_m[cell])?,
            ];
            let operating_point = if let Some(recoil) = self.recoil_materials[cell]
                .as_ref()
                .filter(|_| include_remanence)
            {
                let n = recoil.direction();
                let parallel = finite(n[0].mul_add(h[0], n[1] * h[1]))?;
                let transverse = finite(n[0].mul_add(h[1], -n[1] * h[0]))?.abs();
                Some(
                    recoil
                        .material()
                        .admit_operating_point(parallel, transverse)
                        .map_err(MagneticError::Material)?,
                )
            } else {
                None
            };
            let volume = self
                .mesh
                .face_measure(cell)
                .ok_or(MagneticError::Invalid("missing triangle measure"))?;
            let energy = finite(0.5 * volume * induced_b[0].mul_add(h[0], induced_b[1] * h[1]))?;
            if energy < 0.0 || (induced_b != [0.0; 2] && energy == 0.0) {
                return Err(MagneticError::Invalid(
                    "field energy is not representably positive",
                ));
            }
            energy_j = finite(energy_j + energy)?;
            magnetic_flux_density_t.push(b);
            magnetic_field_a_m.push(h);
            cell_energy_j.push(energy);
            recoil_operating_points.push(operating_point);
        }
        poll(cancelled)?;
        Ok(PlanarField {
            problem: self,
            impressed_current_density_a_m2: None,
            potential_wb_m: potential,
            magnetic_flux_density_t,
            magnetic_field_a_m,
            cell_energy_j,
            recoil_operating_points,
            energy_j,
            report,
            incremental_response: !include_remanence,
        })
    }

    fn admit_mesh_gauges(
        mesh: &TriComplex2,
        prescribed: &[Option<f64>],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), MagneticError> {
        let n = prescribed.len();
        let mut neighbors = vec![Vec::new(); n];
        for &[a, b] in mesh.edges() {
            poll(cancelled)?;
            neighbors[a as usize].push(b as usize);
            neighbors[b as usize].push(a as usize);
        }
        let mut visited = vec![false; n];
        let mut stack = Vec::new();
        for start in 0..n {
            if visited[start] {
                continue;
            }
            visited[start] = true;
            stack.push(start);
            let mut anchored = false;
            while let Some(v) = stack.pop() {
                poll(cancelled)?;
                anchored |= prescribed[v].is_some();
                for &w in &neighbors[v] {
                    if !visited[w] {
                        visited[w] = true;
                        stack.push(w);
                    }
                }
            }
            if !anchored {
                return Err(MagneticError::Invalid(
                    "every connected component needs a prescribed exterior edge",
                ));
            }
        }
        Ok(())
    }
}

pub(crate) fn gradients(points: [[f64; 2]; 3]) -> Result<[[f64; 2]; 3], MagneticError> {
    let [a, b, c] = points;
    let determinant = finite((b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0]))?;
    if determinant == 0.0 {
        return Err(MagneticError::Invalid("singular triangle"));
    }
    let gradients = [
        [(b[1] - c[1]) / determinant, (c[0] - b[0]) / determinant],
        [(c[1] - a[1]) / determinant, (a[0] - c[0]) / determinant],
        [(a[1] - b[1]) / determinant, (b[0] - a[0]) / determinant],
    ];
    if gradients.iter().flatten().all(|v| v.is_finite()) {
        Ok(gradients)
    } else {
        Err(MagneticError::Invalid("unrepresentable shape gradients"))
    }
}

/// A computed field tied to the exact immutable geometry and source states.
/// Cell fields are constant within each P1 triangle, not averaged across a
/// material interface. Energy is the discrete integral (B-Br).H/2 at fixed
/// state, where Br is zero for induced cells. The reference is explicit.
pub struct PlanarField<'a> {
    pub(crate) problem: &'a PlanarMagnetostatics,
    // Present after a static solve or a winding's static superposition.
    // Induction recovery alone cannot establish a current-free force contour.
    impressed_current_density_a_m2: Option<Vec<f64>>,
    potential_wb_m: Vec<f64>,
    magnetic_flux_density_t: Vec<[f64; 2]>,
    magnetic_field_a_m: Vec<[f64; 2]>,
    cell_energy_j: Vec<f64>,
    recoil_operating_points: Vec<Option<RecoilOperatingPoint>>,
    energy_j: f64,
    report: Option<SolveReport>,
    incremental_response: bool,
}
impl PlanarField<'_> {
    /// Original mesh, sources and coefficient state.
    pub fn problem(&self) -> &PlanarMagnetostatics {
        self.problem
    }
    /// Nodal A_z [Wb/m] in original mesh storage order.
    pub fn potential_wb_m(&self) -> &[f64] {
        &self.potential_wb_m
    }
    /// Cell B=(d_y A_z,-d_x A_z) [T] in the stored planar frame.
    pub fn magnetic_flux_density_t(&self) -> &[[f64; 2]] {
        &self.magnetic_flux_density_t
    }
    /// Physical H=(B-Br)/mu [A/m], or the derivative B/mu when
    /// `is_incremental_response()` is true. Source data stay with the problem.
    pub fn magnetic_field_a_m(&self) -> &[[f64; 2]] {
        &self.magnetic_field_a_m
    }
    /// Actual impressed J_z [A/m2] from a static solve or winding superposition,
    /// in triangle order.
    /// Actual fields also retain operating-admitted recoil sources through the
    /// original problem; zero J does not imply absence of permanent magnetization.
    /// Incremental responses omit remanence and do not establish a magnet's
    /// operating field even when this current snapshot is present.
    /// Absent for fields recovered from induction history; that history does
    /// not establish the static, current-free exterior used by Maxwell force.
    pub fn impressed_current_density_a_m2(&self) -> Option<&[f64]> {
        self.impressed_current_density_a_m2.as_deref()
    }
    /// Per-cell stored field energy [J] over actual planar thickness.
    pub fn cell_energy_j(&self) -> &[f64] {
        &self.cell_energy_j
    }
    /// Sum of stored field energy [J]. No mechanical force is inferred from it.
    pub fn energy_j(&self) -> f64 {
        self.energy_j
    }
    /// Reference of this field's constitutive storage. Remanence is held fixed;
    /// changing temperature or magnetization requires separate parameter work.
    pub fn energy_reference(&self) -> MagneticEnergyReference {
        if !self.incremental_response && self.problem.has_permanent_magnets() {
            MagneticEnergyReference::FixedRemanence
        } else {
            MagneticEnergyReference::ZeroFlux
        }
    }
    /// Admitted physical H components and remaining source-envelope margins
    /// per original triangle. Incremental responses have no operating-point
    /// entries. These are numerical estimates, not certified field bounds.
    pub fn recoil_operating_points(&self) -> &[Option<RecoilOperatingPoint>] {
        &self.recoil_operating_points
    }
    /// Whether this is a current derivative with remanence omitted. An
    /// incremental response has H=B/mu and no magnet operating-point claim;
    /// its source limits apply only after combining it with the bias field.
    pub fn is_incremental_response(&self) -> bool {
        self.incremental_response
    }
    /// Recomputed free-DOF residual report. Absent when all DOFs are prescribed
    /// or when a winding owner explicitly reconstructs accepted basis fields.
    pub fn report(&self) -> Option<&SolveReport> {
        self.report.as_ref()
    }
    pub(crate) fn thickness_m(&self) -> f64 {
        self.problem.thickness_m
    }
}
