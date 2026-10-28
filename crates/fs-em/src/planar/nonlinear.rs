//! Reversible H(B) and oriented recoil magnets on the planar P1/FEEC substrate.
//!
//! The residual is integral H(curl A).curl(test) - J*test. Its Newton tangent
//! has eigenvalues H/B and dH/dB, supplied by the SAME material law as energy.
//! Weighted Whitney assembly and checked CG remain the existing owners.
//! Bounded energy-decreasing backtracking never extrapolates a material law.
//! This is fixed-state magnetostatics: no hysteresis, heat, induction history,
//! motion or dynamic parameter-work interpretation is attached to this solve.

use super::{
    MagneticError, PlanarMagnetostatics, finite, gradients,
    permanent::{MagneticEnergyReference, PlanarRecoilMaterial},
    poll,
};
use fs_feec::{
    assembly::{incidence_to_csr, stiffness},
    weighted::{WeightedAssemblyLimits, WeightedError},
    weighted2::{CellWeight2, weighted_mass_matrix_2d},
};
use fs_material::{
    field_storage::{
        linear_magnetic_permeability,
        nonlinear::{ReversibleMagneticMaterial, ReversibleMagneticResponse},
        permanent::RecoilOperatingPoint,
    },
    state_point::ResolvedMaterialStatePoint,
};
use fs_qty::Dims;
use fs_rep_mesh::{Metric2, TriComplex2};
use fs_solver::{CheckedCgConfig, CheckedCgError, CsrOp, LinearOp, checked_cg};
use fs_sparse::{Csr, precond::IdentityPrecond};
use std::cell::{Cell, RefCell};

/// An explicit per-triangle physical law, retaining its original source owner.
#[derive(Clone, Debug, PartialEq)]
pub enum IsotropicMagneticMaterial {
    /// Existing source-resolved, fixed-state linear permeability.
    Linear(ResolvedMaterialStatePoint),
    /// Monotone reversible H(|B|), with zero remanence.
    Reversible(ReversibleMagneticMaterial),
    /// Source-backed B=mu_recoil H+Br with explicit specimen orientation and
    /// reversible operating envelope. Remanence is fixed during a trajectory.
    Recoil(PlanarRecoilMaterial),
}
impl IsotropicMagneticMaterial {
    /// Immutable source snapshot of this physical law.
    pub fn material(&self) -> &ResolvedMaterialStatePoint {
        match self {
            Self::Linear(material) => material,
            Self::Reversible(law) => law.material(),
            Self::Recoil(law) => law.material().material(),
        }
    }
    fn response(&self, b: f64) -> Result<ReversibleMagneticResponse, MagneticError> {
        match self {
            Self::Reversible(law) => law.evaluate(b).map_err(MagneticError::Material),
            Self::Linear(_) | Self::Recoil(_) => {
                let mu = match self {
                    Self::Linear(material) => {
                        linear_magnetic_permeability(material).map_err(MagneticError::Material)?
                    }
                    Self::Recoil(law) => law.material().permeability_h_m(),
                    Self::Reversible(_) => unreachable!(),
                };
                let h = finite(b / mu)?;
                let energy = finite((0.5 * b) * h)?;
                let reluctivity = finite(1.0 / mu)?;
                if reluctivity <= 0.0 || b > 0.0 && (h <= 0.0 || energy <= 0.0) {
                    return Err(MagneticError::Invalid(
                        "linear magnetic response is unrepresentable",
                    ));
                }
                Ok(ReversibleMagneticResponse {
                    field_strength_a_m: h,
                    secant_reluctivity_m_h: reluctivity,
                    differential_reluctivity_m_h: reluctivity,
                    energy_density_j_m3: energy,
                    coenergy_density_j_m3: energy,
                })
            }
        }
    }
}

/// Accuracy and finite work budgets; none is increased on nonconvergence.
#[derive(Clone, Copy, Debug)]
pub struct NonlinearMagneticControls {
    /// Euclidean norm of the free-coordinate physical residual [A m].
    pub absolute_residual_tolerance_a_m: f64,
    /// Relative target against max(initial residual, impressed-current load).
    pub relative_residual_tolerance: f64,
    /// Maximum accepted Newton updates. Zero permits checking an initial field.
    pub max_newton_iterations: usize,
    /// Maximum halvings of each Newton step after its full-step trial.
    pub max_backtracks: usize,
    /// Existing checked CG controls for EACH Newton direction.
    pub linear: CheckedCgConfig,
}

/// Recomputed nonlinear residual and actual bounded numerical work.
#[derive(Clone, Debug)]
pub struct NonlinearMagneticReport {
    /// Accepted nonlinear updates.
    pub newton_iterations: usize,
    /// Total Krylov iterations across all Newton directions.
    pub linear_iterations: usize,
    /// Total step halvings, including those in a failed line search.
    pub backtracks: usize,
    /// Recomputed Euclidean free-coordinate residual [A m].
    pub residual_norm_a_m: f64,
    /// Actual absolute acceptance target [A m].
    pub residual_target_a_m: f64,
}

/// Source/assembly refusal or an exhausted nonlinear solve; no field is returned.
#[derive(Debug)]
pub enum NonlinearMagneticError {
    /// Original source, geometry, checked linear solver or cancellation error.
    Magnetic(MagneticError),
    /// The initial field or a trial left one source's declared field envelope.
    FieldRange {
        /// Original mesh triangle.
        triangle: usize,
        /// Computed flux-density magnitude [T].
        observed_t: f64,
        /// Original source range's upper endpoint [T].
        maximum_t: f64,
    },
    /// The physical H left a source-backed reversible recoil envelope.
    RecoilRange {
        /// Original mesh triangle.
        triangle: usize,
        /// H projected onto the declared magnet direction [A/m].
        parallel_field_a_m: f64,
        /// Magnitude of the perpendicular component of H [A/m].
        transverse_field_a_m: f64,
        /// Original source interval for the parallel component [A/m].
        parallel_range_a_m: [f64; 2],
        /// Original source maximum transverse field [A/m].
        maximum_transverse_field_a_m: f64,
    },
    /// Newton or line-search work exhausted before the physical residual passed.
    Nonconvergence(NonlinearMagneticReport),
}
impl From<MagneticError> for NonlinearMagneticError {
    fn from(error: MagneticError) -> Self {
        Self::Magnetic(error)
    }
}
impl core::fmt::Display for NonlinearMagneticError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "nonlinear planar magnetostatics: {self:?}")
    }
}
impl std::error::Error for NonlinearMagneticError {}

/// Immutable geometry and material law collection for one static state.
pub struct PlanarNonlinearMagnetostatics {
    mesh: TriComplex2,
    materials: Vec<IsotropicMagneticMaterial>,
    incidence: Csr,
    limits: WeightedAssemblyLimits,
    thickness_m: f64,
}

struct Evaluation {
    b: Vec<[f64; 2]>,
    h: Vec<[f64; 2]>,
    energies: Vec<f64>,
    coenergies: Vec<f64>,
    recoil_operating_points: Vec<Option<RecoilOperatingPoint>>,
    tangent_weights: Vec<CellWeight2>,
    residual: Vec<f64>,
    energy: f64,
    coenergy: f64,
    objective: f64,
}

// A positive quadratic contribution supplied by another physical owner, e.g.
// the SAME conductor projection C_sigma/dt in a backward-Euler induction step.
// It contributes residual, tangent and objective together; the magnetic source
// law, FEEC assembly, Newton iteration and checked CG are not duplicated.
pub(crate) trait PositiveQuadratic {
    fn reference(&self) -> &[f64];
    fn apply(
        &self,
        x: &[f64],
        y: &mut [f64],
        cancelled: &dyn Fn() -> bool,
    ) -> Result<(), MagneticError>;
}
struct Tangent<'a> {
    magnetic: CsrOp,
    quadratic: Option<&'a dyn PositiveQuadratic>,
    free: &'a [usize],
    full_len: usize,
    cancelled: &'a dyn Fn() -> bool,
    failure: &'a RefCell<Option<MagneticError>>,
}
impl LinearOp for Tangent<'_> {
    fn n(&self) -> usize {
        self.free.len()
    }
    fn apply(&self, x: &[f64], y: &mut [f64]) {
        self.magnetic.apply(x, y);
        if let Some(quadratic) = self.quadratic {
            let mut full = vec![0.0; self.full_len];
            let mut applied = vec![0.0; self.full_len];
            for (&v, &value) in self.free.iter().zip(x) {
                full[v] = value;
            }
            if let Err(error) = quadratic.apply(&full, &mut applied, self.cancelled) {
                *self.failure.borrow_mut() = Some(error);
                y.fill(f64::NAN);
                return;
            }
            for (value, &v) in y.iter_mut().zip(self.free) {
                *value += applied[v];
            }
        }
    }
}

fn weighted_error(error: WeightedError) -> MagneticError {
    match error {
        WeightedError::Cancelled => MagneticError::Cancelled,
        other => MagneticError::Assembly(other),
    }
}
fn norm(values: impl Iterator<Item = f64>) -> f64 {
    values.fold(0.0, f64::hypot)
}

impl PlanarNonlinearMagnetostatics {
    /// Bind each original triangle to its actual induced or recoil source.
    /// Material names never select coefficients or a solver branch.
    pub fn new(
        mesh: TriComplex2,
        materials: Vec<IsotropicMagneticMaterial>,
        limits: WeightedAssemblyLimits,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, NonlinearMagneticError> {
        poll(cancelled)?;
        let thickness_m = match mesh.metric() {
            Metric2::Planar { thickness } => thickness,
            Metric2::Axisymmetric { .. } => {
                return Err(MagneticError::Invalid("planar A_z requires planar metric").into());
            }
        };
        if materials.len() != mesh.faces().len()
            || mesh.vertices().len() > limits.max_dofs
            || mesh.edges().len() > limits.max_dofs
            || mesh.faces().len() > limits.max_cells
            || mesh
                .faces()
                .len()
                .checked_mul(9)
                .is_none_or(|n| n > limits.max_triplets)
        {
            return Err(MagneticError::Invalid(
                "one material per triangle and sufficient assembly limits required",
            )
            .into());
        }
        let mut weights = Vec::with_capacity(materials.len());
        for material in &materials {
            poll(cancelled)?;
            weights.push(CellWeight2::Scalar(
                material.response(0.0)?.differential_reluctivity_m_h,
            ));
        }
        // The existing Hodge owner validates cell geometry and positivity even
        // when a later solve has every degree of freedom prescribed.
        weighted_mass_matrix_2d(
            &mesh,
            1,
            &weights,
            Dims([-1, -1, 2, 0, 2, 0]),
            limits,
            cancelled,
        )
        .map_err(weighted_error)?;
        let incidence = incidence_to_csr(&mesh.d0());
        poll(cancelled)?;
        Ok(Self {
            mesh,
            materials,
            incidence,
            limits,
            thickness_m,
        })
    }
    /// Original mesh, planar metric and durable geometry identities.
    pub fn mesh(&self) -> &TriComplex2 {
        &self.mesh
    }
    /// Complete source laws in original triangle order.
    pub fn materials(&self) -> &[IsotropicMagneticMaterial] {
        &self.materials
    }
    pub(crate) fn thickness_m(&self) -> f64 {
        self.thickness_m
    }

    fn evaluate(
        &self,
        potential: &[f64],
        rhs: &[f64],
        gauge: f64,
        quadratic: Option<&dyn PositiveQuadratic>,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Evaluation, NonlinearMagneticError> {
        let count = self.materials.len();
        let mut result = Evaluation {
            b: Vec::with_capacity(count),
            h: Vec::with_capacity(count),
            energies: Vec::with_capacity(count),
            coenergies: Vec::with_capacity(count),
            recoil_operating_points: Vec::with_capacity(count),
            tangent_weights: Vec::with_capacity(count),
            residual: rhs.iter().map(|x| -x).collect(),
            energy: 0.0,
            coenergy: 0.0,
            objective: 0.0,
        };
        for (triangle, (face, material)) in
            self.mesh.faces().iter().zip(&self.materials).enumerate()
        {
            poll(cancelled)?;
            let shape = gradients(face.map(|v| self.mesh.vertices()[v as usize]))?;
            let a0 = potential[face[0] as usize];
            let gradient = [0, 1].map(|axis| {
                (potential[face[1] as usize] - a0).mul_add(
                    shape[1][axis],
                    (potential[face[2] as usize] - a0) * shape[2][axis],
                )
            });
            let b = [finite(gradient[1])?, finite(-gradient[0])?];
            let remanence = match material {
                IsotropicMagneticMaterial::Recoil(law) => law.remanence_t(),
                _ => [0.0; 2],
            };
            let induced_b = [finite(b[0] - remanence[0])?, finite(b[1] - remanence[1])?];
            let magnitude = finite(induced_b[0].hypot(induced_b[1]))?;
            if let IsotropicMagneticMaterial::Reversible(law) = material
                && magnitude > law.maximum_flux_density_t()
            {
                return Err(NonlinearMagneticError::FieldRange {
                    triangle,
                    observed_t: magnitude,
                    maximum_t: law.maximum_flux_density_t(),
                });
            }
            let response = material.response(magnitude)?;
            let nu = response.secant_reluctivity_m_h;
            let differential = response.differential_reluctivity_m_h;
            let h = [finite(nu * induced_b[0])?, finite(nu * induced_b[1])?];
            let operating_point = if let IsotropicMagneticMaterial::Recoil(recoil) = material {
                let direction = recoil.direction();
                let parallel = finite(direction[0].mul_add(h[0], direction[1] * h[1]))?;
                let transverse = finite(direction[0].mul_add(h[1], -direction[1] * h[0]))?.abs();
                let source = recoil.material();
                let range = source.parallel_range_a_m();
                if parallel < range[0]
                    || parallel > range[1]
                    || transverse > source.maximum_transverse_a_m()
                {
                    return Err(NonlinearMagneticError::RecoilRange {
                        triangle,
                        parallel_field_a_m: parallel,
                        transverse_field_a_m: transverse,
                        parallel_range_a_m: range,
                        maximum_transverse_field_a_m: source.maximum_transverse_a_m(),
                    });
                }
                Some(
                    source
                        .admit_operating_point(parallel, transverse)
                        .map_err(MagneticError::Material)?,
                )
            } else {
                None
            };
            let volume = self
                .mesh
                .face_measure(triangle)
                .ok_or(MagneticError::Invalid("missing triangle measure"))?;
            let energy = finite(volume * response.energy_density_j_m3)?;
            // The Legendre dual of shifted recoil storage is
            // mu*|H|^2/2 + Br.H. It can be negative; B.H/2 is NOT its energy.
            let coenergy_density = finite(remanence[0].mul_add(
                h[0],
                remanence[1].mul_add(h[1], response.coenergy_density_j_m3),
            ))?;
            let coenergy = finite(volume * coenergy_density)?;
            result.energy = finite(result.energy + energy)?;
            result.coenergy = finite(result.coenergy + coenergy)?;
            for (local, &vertex) in face.iter().enumerate() {
                let pairing = h[0].mul_add(shape[local][1], -h[1] * shape[local][0]);
                result.residual[vertex as usize] =
                    finite(volume.mul_add(pairing, result.residual[vertex as usize]))?;
            }
            // Hodge acts on grad A, whereas the constitutive tangent acts on
            // B=rot(grad A). Rotate its positive tensor into the gradient frame.
            let weight = if magnitude == 0.0 || differential == nu {
                CellWeight2::Scalar(nu)
            } else {
                let direction = [gradient[0] / magnitude, gradient[1] / magnitude];
                let delta = differential - nu;
                let off_diagonal = delta * direction[0] * direction[1];
                CellWeight2::Tensor([
                    [delta.mul_add(direction[0] * direction[0], nu), off_diagonal],
                    [off_diagonal, delta.mul_add(direction[1] * direction[1], nu)],
                ])
            };
            result.b.push(b);
            result.h.push(h);
            result.energies.push(energy);
            result.coenergies.push(coenergy);
            result.recoil_operating_points.push(operating_point);
            result.tangent_weights.push(weight);
        }
        result.objective = result.energy;
        for (&load, &value) in rhs.iter().zip(potential) {
            result.objective = finite((-load).mul_add(value - gauge, result.objective))?;
        }
        if let Some(quadratic) = quadratic {
            let increment = potential
                .iter()
                .zip(quadratic.reference())
                .map(|(value, reference)| finite(value - reference))
                .collect::<Result<Vec<_>, _>>()?;
            let mut applied = vec![0.0; potential.len()];
            let cancel = RefCell::new(&mut *cancelled);
            quadratic.apply(&increment, &mut applied, &|| (*cancel.borrow_mut())())?;
            for (v, (&delta, &effort)) in increment.iter().zip(&applied).enumerate() {
                result.residual[v] = finite(result.residual[v] + effort)?;
                result.objective = finite((0.5 * delta).mul_add(effort, result.objective))?;
            }
        }
        Ok(result)
    }

    pub(crate) fn magnetic_load(
        &self,
        potential: &[f64],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Vec<f64>, NonlinearMagneticError> {
        Ok(self
            .evaluate(potential, &vec![0.0; potential.len()], 0.0, None, cancelled)?
            .residual)
    }

    fn assemble_tangent(
        &self,
        weights: &[CellWeight2],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Csr, NonlinearMagneticError> {
        let mass = weighted_mass_matrix_2d(
            &self.mesh,
            1,
            weights,
            Dims([-1, -1, 2, 0, 2, 0]),
            self.limits,
            cancelled,
        )
        .map_err(weighted_error)?;
        Ok(stiffness(&self.incidence, mass.matrix()))
    }

    // The moving-conductor Newton solve adds a genuinely nonsymmetric
    // transport derivative. Its magnetic residual and constitutive Hessian
    // remain exactly this owner's source-admitted FEEC formulation.
    pub(crate) fn magnetic_linearization(
        &self,
        potential: &[f64],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(Vec<f64>, Csr), NonlinearMagneticError> {
        let evaluation =
            self.evaluate(potential, &vec![0.0; potential.len()], 0.0, None, cancelled)?;
        let tangent = self.assemble_tangent(&evaluation.tangent_weights, cancelled)?;
        Ok((evaluation.residual, tangent))
    }

    pub(crate) fn recover(
        &self,
        potential: Vec<f64>,
        report: NonlinearMagneticReport,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<PlanarNonlinearField<'_>, NonlinearMagneticError> {
        let evaluation = self.evaluate(
            &potential,
            &vec![0.0; potential.len()],
            0.0,
            None,
            cancelled,
        )?;
        poll(cancelled)?;
        Ok(PlanarNonlinearField {
            problem: self,
            potential_wb_m: potential,
            magnetic_flux_density_t: evaluation.b,
            magnetic_field_a_m: evaluation.h,
            cell_energy_j: evaluation.energies,
            cell_coenergy_j: evaluation.coenergies,
            recoil_operating_points: evaluation.recoil_operating_points,
            energy_j: evaluation.energy,
            coenergy_j: evaluation.coenergy,
            report,
        })
    }

    /// Solve the physical nonlinear residual with source-derived Newton tangents.
    ///
    /// Boundary semantics are exactly those of `PlanarMagnetostatics`. The
    /// complete finite initial field must match prescribed potentials and stay
    /// within all source ranges; an invalid initial field is never clamped.
    /// All trial states are local. Failed convergence, material admission or
    /// cancellation returns no replacement. The final report recomputes the
    /// nonlinear residual; a small Krylov residual alone cannot accept a field.
    pub fn solve(
        &self,
        current_density_a_m2: &[f64],
        boundary_potential_wb_m: &[(usize, f64)],
        initial_potential_wb_m: &[f64],
        controls: NonlinearMagneticControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<PlanarNonlinearField<'_>, NonlinearMagneticError> {
        poll(cancelled)?;
        if current_density_a_m2.len() != self.mesh.faces().len()
            || current_density_a_m2.iter().any(|v| !v.is_finite())
        {
            return Err(
                MagneticError::Invalid("complete finite impressed current required").into(),
            );
        }
        let mut rhs = vec![0.0; self.mesh.vertices().len()];
        for (triangle, face) in self.mesh.faces().iter().enumerate() {
            poll(cancelled)?;
            let volume = self
                .mesh
                .face_measure(triangle)
                .ok_or(MagneticError::Invalid("missing triangle measure"))?;
            let load = finite(current_density_a_m2[triangle] * (volume / 3.0))?;
            for &v in face {
                rhs[v as usize] = finite(rhs[v as usize] + load)?;
            }
        }
        self.solve_loaded(
            &rhs,
            boundary_potential_wb_m,
            initial_potential_wb_m,
            None,
            controls,
            cancelled,
        )
    }

    pub(crate) fn solve_loaded(
        &self,
        rhs: &[f64],
        boundary_potential_wb_m: &[(usize, f64)],
        initial_potential_wb_m: &[f64],
        quadratic: Option<&dyn PositiveQuadratic>,
        controls: NonlinearMagneticControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<PlanarNonlinearField<'_>, NonlinearMagneticError> {
        poll(cancelled)?;
        let n = self.mesh.vertices().len();
        if rhs.len() != n
            || rhs.iter().any(|v| !v.is_finite())
            || initial_potential_wb_m.len() != n
            || initial_potential_wb_m.iter().any(|v| !v.is_finite())
            || quadratic.is_some_and(|q| {
                q.reference().len() != n || q.reference().iter().any(|v| !v.is_finite())
            })
            || !controls.absolute_residual_tolerance_a_m.is_finite()
            || controls.absolute_residual_tolerance_a_m <= 0.0
            || !controls.relative_residual_tolerance.is_finite()
            || !(0.0..1.0).contains(&controls.relative_residual_tolerance)
            || !controls.linear.tolerance.is_finite()
            || !(0.0..1.0).contains(&controls.linear.tolerance)
            || controls.linear.tolerance == 0.0
        {
            return Err(MagneticError::Invalid(
                "nonlinear field, current or residual controls are invalid",
            )
            .into());
        }
        let prescribed = PlanarMagnetostatics::admit_mesh_boundary(
            &self.mesh,
            boundary_potential_wb_m,
            cancelled,
        )?;
        if prescribed
            .iter()
            .zip(initial_potential_wb_m)
            .any(|(p, value)| p.is_some_and(|p| p != *value))
        {
            return Err(MagneticError::Invalid(
                "initial field must match every prescribed potential",
            )
            .into());
        }
        let gauge = prescribed.iter().flatten().next().copied().unwrap_or(0.0);
        let free: Vec<_> = (0..n).filter(|&v| prescribed[v].is_none()).collect();
        let mut index = vec![usize::MAX; n];
        for (i, &v) in free.iter().enumerate() {
            index[v] = i;
        }
        let mut potential = initial_potential_wb_m.to_vec();
        let mut evaluation = self.evaluate(&potential, rhs, gauge, quadratic, cancelled)?;
        let initial_norm = norm(free.iter().map(|&v| evaluation.residual[v]));
        let reference_norm = initial_norm.max(norm(free.iter().map(|&v| rhs[v])));
        let target = controls.absolute_residual_tolerance_a_m.max(finite(
            controls.relative_residual_tolerance * reference_norm,
        )?);
        let mut report = NonlinearMagneticReport {
            newton_iterations: 0,
            linear_iterations: 0,
            backtracks: 0,
            residual_norm_a_m: initial_norm,
            residual_target_a_m: target,
        };
        loop {
            poll(cancelled)?;
            report.residual_norm_a_m = finite(norm(free.iter().map(|&v| evaluation.residual[v])))?;
            if report.residual_norm_a_m <= target {
                return Ok(PlanarNonlinearField {
                    problem: self,
                    potential_wb_m: potential,
                    magnetic_flux_density_t: evaluation.b,
                    magnetic_field_a_m: evaluation.h,
                    cell_energy_j: evaluation.energies,
                    cell_coenergy_j: evaluation.coenergies,
                    recoil_operating_points: evaluation.recoil_operating_points,
                    energy_j: evaluation.energy,
                    coenergy_j: evaluation.coenergy,
                    report,
                });
            }
            if report.newton_iterations == controls.max_newton_iterations {
                return Err(NonlinearMagneticError::Nonconvergence(report));
            }
            let tangent = self.assemble_tangent(&evaluation.tangent_weights, cancelled)?;
            let mut row_ptr = vec![0];
            let mut columns = Vec::new();
            let mut values = Vec::new();
            for &v in &free {
                poll(cancelled)?;
                let (cols, vals) = tangent.row(v);
                for (&w, &value) in cols.iter().zip(vals) {
                    if index[w] != usize::MAX {
                        columns.push(index[w]);
                        values.push(finite(value)?);
                    }
                }
                row_ptr.push(columns.len());
            }
            let matrix = Csr::from_parts(free.len(), free.len(), row_ptr, columns, values);
            let direction_rhs: Vec<_> = free.iter().map(|&v| -evaluation.residual[v]).collect();
            let direction = {
                let cancel = RefCell::new(&mut *cancelled);
                let interrupted = Cell::new(false);
                let check = || {
                    if interrupted.get() || (*cancel.borrow_mut())() {
                        interrupted.set(true);
                        true
                    } else {
                        false
                    }
                };
                let failure = RefCell::new(None);
                let tangent = Tangent {
                    magnetic: CsrOp::symmetric(matrix),
                    quadratic,
                    free: &free,
                    full_len: n,
                    cancelled: &check,
                    failure: &failure,
                };
                let result = checked_cg(
                    &tangent,
                    &IdentityPrecond,
                    &direction_rhs,
                    controls.linear,
                    |_| {
                        if check() {
                            Err(MagneticError::Cancelled)
                        } else {
                            Ok(())
                        }
                    },
                );
                if let Some(error) = failure.into_inner() {
                    return Err(error.into());
                }
                result.map_err(|error| match error {
                    CheckedCgError::InvalidInput(why) => MagneticError::Invalid(why),
                    CheckedCgError::Interrupted(error) => error,
                })?
            };
            if !direction.report.converged_euclidean() {
                return Err(MagneticError::Nonconvergence(direction.report).into());
            }
            report.linear_iterations += direction.report.iters;
            let mut descent = 0.0;
            for (&v, &delta) in free.iter().zip(&direction.x) {
                descent = finite(evaluation.residual[v].mul_add(delta, descent))?;
            }
            if descent >= 0.0 {
                return Err(NonlinearMagneticError::Nonconvergence(report));
            }
            let mut alpha = 1.0_f64;
            let mut accepted = None;
            for backtrack in 0..=controls.max_backtracks {
                poll(cancelled)?;
                if backtrack > 0 {
                    report.backtracks += 1;
                }
                let mut trial = potential.clone();
                for (&v, &delta) in free.iter().zip(&direction.x) {
                    trial[v] = finite(alpha.mul_add(delta, potential[v]))?;
                }
                match self.evaluate(&trial, rhs, gauge, quadratic, cancelled) {
                    Ok(next) => {
                        let residual = norm(free.iter().map(|&v| next.residual[v]));
                        if residual <= target
                            || next.objective
                                <= (1e-4 * alpha).mul_add(descent, evaluation.objective)
                        {
                            accepted = Some((trial, next));
                            break;
                        }
                    }
                    Err(
                        NonlinearMagneticError::FieldRange { .. }
                        | NonlinearMagneticError::RecoilRange { .. },
                    ) => {}
                    Err(error) => return Err(error),
                }
                alpha *= 0.5;
            }
            let Some((next_potential, next_evaluation)) = accepted else {
                return Err(NonlinearMagneticError::Nonconvergence(report));
            };
            potential = next_potential;
            evaluation = next_evaluation;
            report.newton_iterations += 1;
        }
    }
}

/// Actual nonlinear field, bound to the source cards and geometry that produced it.
pub struct PlanarNonlinearField<'a> {
    problem: &'a PlanarNonlinearMagnetostatics,
    potential_wb_m: Vec<f64>,
    magnetic_flux_density_t: Vec<[f64; 2]>,
    magnetic_field_a_m: Vec<[f64; 2]>,
    cell_energy_j: Vec<f64>,
    cell_coenergy_j: Vec<f64>,
    recoil_operating_points: Vec<Option<RecoilOperatingPoint>>,
    energy_j: f64,
    coenergy_j: f64,
    report: NonlinearMagneticReport,
}
impl<'a> PlanarNonlinearField<'a> {
    /// Exact immutable geometry and source owner that produced this field.
    pub fn problem(&self) -> &'a PlanarNonlinearMagnetostatics {
        self.problem
    }
    /// Nodal A_z [Wb/m], preserving original vertex order.
    pub fn potential_wb_m(&self) -> &[f64] {
        &self.potential_wb_m
    }
    /// B=(d_y A_z,-d_x A_z) [T], one value per original triangle.
    pub fn magnetic_flux_density_t(&self) -> &[[f64; 2]] {
        &self.magnetic_flux_density_t
    }
    /// Actual constitutive H [A/m], without averaging material interfaces.
    pub fn magnetic_field_a_m(&self) -> &[[f64; 2]] {
        &self.magnetic_field_a_m
    }
    /// Source-integrated magnetic storage [J], including extrusion thickness.
    pub fn cell_energy_j(&self) -> &[f64] {
        &self.cell_energy_j
    }
    /// Per-cell Legendre-dual coenergy [J]; no motion derivative is implied.
    pub fn cell_coenergy_j(&self) -> &[f64] {
        &self.cell_coenergy_j
    }
    /// Sum of source-integrated energies [J] at the declared energy reference.
    pub fn energy_j(&self) -> f64 {
        self.energy_j
    }
    /// Sum of original cell magnetic coenergies [J].
    pub fn coenergy_j(&self) -> f64 {
        self.coenergy_j
    }
    /// Storage reference; recoil cells have zero energy at B=Br and induced
    /// cells retain their source-integrated zero-flux reference.
    pub fn energy_reference(&self) -> MagneticEnergyReference {
        if self.recoil_operating_points.iter().any(Option::is_some) {
            MagneticEnergyReference::FixedRemanence
        } else {
            MagneticEnergyReference::ZeroFlux
        }
    }
    /// Actual H components and remaining source-envelope margins, per triangle.
    /// These are numerical operating estimates, not certified field bounds.
    pub fn recoil_operating_points(&self) -> &[Option<RecoilOperatingPoint>] {
        &self.recoil_operating_points
    }
    /// Recomputed nonlinear residual and actual solve work.
    pub fn report(&self) -> &NonlinearMagneticReport {
        &self.report
    }
    pub(crate) fn thickness_m(&self) -> f64 {
        self.problem.thickness_m
    }
}
