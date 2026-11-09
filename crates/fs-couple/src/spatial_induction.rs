//! Reciprocal three-dimensional induction and finite rotor inertia.
//!
//! The original Whitney A–phi solve and original articulated impulse map use
//! one converged angular velocity. Conductor-frame Joule heat stays separate
//! from Lorentz mechanical work; the same local currents supply both accounts.
//! Every rejected speed trial restarts the accepted field and mechanical state.
//!
//! This is a homogeneous, isotropic circular cylinder on a fixed Eulerian mesh.
//! Its analytic occupancy is invariant under rotation. The admitted faceted
//! approximation is a geometry error, not a moving polyhedron or a sliding
//! interface. Uniform scalar permeability throughout the field domain excludes
//! unaccounted magnetic interface tractions. Hysteretic memory, nonaxisymmetric
//! material occupancy, deformation and transported material frames need their
//! own motion models and are not inferred from a prescribed velocity field.

mod geometry;
/// Source-pinned endpoint coefficients shared across rotor-speed trials.
pub mod thermal;

use crate::induction_rotor::{self, RotorStageError};
pub use crate::induction_rotor::{RotorBalance, RotorControls};
use fs_em::spatial::{
    SpatialMagneticError, SpatialMagneticMaterial,
    vector_potential::VectorMagneticBoundary,
    vector_potential::induction::{
        SpatialInduction, SpatialInductionCandidate, SpatialInductionControls,
        SpatialInductionStep,
        mechanical::{SpatialLorentzLoadError, SpatialLorentzLoadLimits},
    },
};
use fs_ga::GaError;
use fs_material::field_storage::spatial::SpatialMagneticLaw;
use fs_mbd::articulated::{ArticulatedError, ArticulatedModel, PrescribedBaseImpulseResponse};
use std::sync::Arc;

/// Analytic cylinder centered on the original revolute-z link's origin.
#[derive(Clone, Copy, Debug)]
pub struct CylindricalConductor {
    /// Zero denotes a solid cylinder; positive values admit a coaxial bore [m].
    pub inner_radius_m: f64,
    pub outer_radius_m: f64,
    pub half_length_m: f64,
    /// Maximum departure of the original boundary facets from the analytic
    /// cylinder/end planes [m]. This does not certify field convergence.
    pub maximum_boundary_error_m: f64,
}

#[derive(Debug)]
pub enum SpatialRotorError {
    Invalid(&'static str),
    Cancelled,
    Field(SpatialMagneticError),
    Mechanical(ArticulatedError),
    Load(SpatialLorentzLoadError),
    Geometry(GaError),
    Nonconvergence {
        iterations: usize,
        velocity_residual_rad_s: f64,
    },
    Balance(&'static str),
}
impl core::fmt::Display for SpatialRotorError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "spatial induction rotor refused: {self:?}")
    }
}
impl std::error::Error for SpatialRotorError {}
impl From<SpatialMagneticError> for SpatialRotorError {
    fn from(error: SpatialMagneticError) -> Self {
        match error {
            SpatialMagneticError::Cancelled => Self::Cancelled,
            other => Self::Field(other),
        }
    }
}
impl From<ArticulatedError> for SpatialRotorError {
    fn from(error: ArticulatedError) -> Self {
        Self::Mechanical(error)
    }
}
impl From<SpatialLorentzLoadError> for SpatialRotorError {
    fn from(error: SpatialLorentzLoadError) -> Self {
        match error {
            SpatialLorentzLoadError::Cancelled => Self::Cancelled,
            other => Self::Load(other),
        }
    }
}
impl From<RotorStageError> for SpatialRotorError {
    fn from(error: RotorStageError) -> Self {
        match error {
            RotorStageError::Invalid(reason) => Self::Invalid(reason),
            RotorStageError::Cancelled => Self::Cancelled,
            RotorStageError::Mechanical(error) => Self::Mechanical(error),
            RotorStageError::Load(error) => error.into(),
            RotorStageError::Geometry(error) => Self::Geometry(error),
            RotorStageError::Balance(reason) => Self::Balance(reason),
        }
    }
}
type Result<T> = core::result::Result<T, SpatialRotorError>;

/// Accepted field, charge-response history and finite-inertia rotor state.
/// Immutable geometry and the original native model are shared across steps.
#[derive(Clone)]
pub struct SpatialInductionRotor<'m> {
    field: SpatialInduction<'m>,
    binding: Arc<RotorBinding>,
    angle_rad: f64,
    angular_velocity_rad_s: f64,
}
struct RotorBinding {
    model: ArticulatedModel,
    cells: Vec<usize>,
    unit_velocity: Vec<[f64; 3]>,
    load_limits: SpatialLorentzLoadLimits,
    boundary_error_m: f64,
}

/// One publishable field/mechanical endpoint with original full-body impulse
/// and ideal-bearing reaction. No predictor load is substituted in the receipt.
pub struct SpatialRotorCandidate<'m> {
    state: SpatialInductionRotor<'m>,
    pub induction: SpatialInductionStep,
    pub mechanics: PrescribedBaseImpulseResponse,
    pub balance: RotorBalance,
}
impl<'m> SpatialRotorCandidate<'m> {
    pub fn state(&self) -> &SpatialInductionRotor<'m> {
        &self.state
    }
    pub fn into_state(self) -> SpatialInductionRotor<'m> {
        self.state
    }
    pub fn into_parts(
        self,
    ) -> (
        SpatialInductionRotor<'m>,
        SpatialInductionStep,
        PrescribedBaseImpulseResponse,
        RotorBalance,
    ) {
        (self.state, self.induction, self.mechanics, self.balance)
    }
}

fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<()> {
    if cancelled() {
        Err(SpatialRotorError::Cancelled)
    } else {
        Ok(())
    }
}
fn finite(value: f64) -> Result<f64> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(SpatialRotorError::Invalid("nonfinite rotor arithmetic"))
    }
}

impl<'m> SpatialInductionRotor<'m> {
    /// Bind all conducting tetrahedra to one homogeneous circular rotor.
    ///
    /// The native model has one centered axisymmetric revolute-z root link.
    /// All conductor cells retain the identical scalar-conductivity source;
    /// insulating cells may carry external impressed currents. The complete
    /// field domain must have the same scalar induced-linear permeability, so
    /// this conduction-current load does not omit a permeability-interface or
    /// permanent-magnet traction. No field mesh or source history is changed.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        field: SpatialInduction<'m>,
        model: ArticulatedModel,
        shape: CylindricalConductor,
        angle_rad: f64,
        angular_velocity_rad_s: f64,
        limits: SpatialLorentzLoadLimits,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self> {
        poll(cancelled)?;
        let problem = field.problem();
        if problem.mesh().vertex_count > limits.max_vertices
            || problem.mesh().tets.len() > limits.max_cells
        {
            return Err(SpatialRotorError::Invalid("rotor mesh exceeds load budget"));
        }
        let origin =
            induction_rotor::admit_model(&model, angle_rad, angular_velocity_rad_s, cancelled)?;
        let cells: Vec<_> = field
            .conductors()
            .iter()
            .enumerate()
            .filter_map(|(cell, source)| source.as_ref().map(|_| cell))
            .collect();
        let Some(&first_cell) = cells.first() else {
            return Err(SpatialRotorError::Invalid(
                "nonempty rotor conductor required",
            ));
        };
        let first = field.conductors()[first_cell].as_ref().unwrap();
        let first_mu = problem.permeability()[first_cell];
        if first_mu.isotropic_h_m().is_none() {
            return Err(SpatialRotorError::Invalid(
                "scalar rotor permeability required",
            ));
        }
        for (cell, material) in problem.materials().iter().enumerate() {
            poll(cancelled)?;
            if !matches!(material, SpatialMagneticMaterial::Linear(_))
                || problem.magnetic_laws()[cell] != SpatialMagneticLaw::Isotropic
                || problem.permeability()[cell] != first_mu
                || field.conductors()[cell].as_ref().is_some_and(|source| {
                    source != first || field.impressed_current_density_a_m2()[cell] != [0.0; 3]
                })
            {
                return Err(SpatialRotorError::Invalid(
                    "uniform scalar permeability and one homogeneous undriven conductor source required",
                ));
            }
        }
        let boundary_error_m = geometry::admit_cylinder(&field, &cells, origin, shape, cancelled)?;
        let mut unit_velocity = Vec::with_capacity(problem.mesh().vertex_count);
        for point in problem.positions_m() {
            poll(cancelled)?;
            unit_velocity.push([
                finite(origin.y - point[1])?,
                finite(point[0] - origin.x)?,
                0.0,
            ]);
        }
        poll(cancelled)?;
        Ok(Self {
            field,
            binding: Arc::new(RotorBinding {
                model,
                cells,
                unit_velocity,
                load_limits: limits,
                boundary_error_m,
            }),
            angle_rad,
            angular_velocity_rad_s,
        })
    }

    pub fn field(&self) -> &SpatialInduction<'m> {
        &self.field
    }
    pub fn model(&self) -> &ArticulatedModel {
        &self.binding.model
    }
    pub const fn angle_rad(&self) -> f64 {
        self.angle_rad
    }
    pub const fn angular_velocity_rad_s(&self) -> f64 {
        self.angular_velocity_rad_s
    }
    pub fn boundary_error_m(&self) -> f64 {
        self.binding.boundary_error_m
    }

    /// Solve motional EMF, charge continuity and finite rotor momentum at one
    /// shared midpoint speed. The native field's P1 current supplies the full
    /// mechanical wrench, its work and conductor-frame heat. Held external
    /// shaft torque is independent of the electromagnetic body load.
    ///
    /// Iterations are bounded and reuse this same accepted field and velocity.
    /// Failed field solves, budgets, balance or cancellation publish neither
    /// history. The existing field Peclet and Krylov limits remain in force.
    #[allow(clippy::too_many_arguments)]
    pub fn step_candidate(
        &self,
        target_s: f64,
        endpoint_current_a_m2: &[[f64; 3]],
        shaft_torque_nm: f64,
        field_controls: SpatialInductionControls,
        controls: RotorControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<SpatialRotorCandidate<'m>> {
        self.step_with(
            target_s,
            endpoint_current_a_m2,
            shaft_torque_nm,
            controls,
            |velocity, checkpoint| {
                self.field
                    .step_candidate(
                        endpoint_current_a_m2,
                        velocity,
                        target_s - self.field.time_s(),
                        field_controls,
                        &mut || checkpoint(),
                    )
                    .map_err(Into::into)
            },
            cancelled,
        )
    }

    /// Drive the same finite-inertia rotor with an applied magnetic history.
    ///
    /// The complete original tangential A trace is prescribed at the endpoint
    /// and varies affinely from the retained trace. Every speed trial includes
    /// the same boundary Faraday field and restarts the accepted history.
    /// The trace uses the original fixed WORLD mesh coordinates; body rotation
    /// enters through the native motional v cross B term at the solved speed.
    /// Boundary reaction work is external electromagnetic supply; it is not
    /// shaft work, Lorentz work, or an additional heat deposit. Native support,
    /// gauge and charge-continuity admission remain in force.
    #[allow(clippy::too_many_arguments)]
    pub fn step_with_boundary_candidate(
        &self,
        target_s: f64,
        endpoint_current_a_m2: &[[f64; 3]],
        endpoint_boundary: &VectorMagneticBoundary,
        shaft_torque_nm: f64,
        field_controls: SpatialInductionControls,
        controls: RotorControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<SpatialRotorCandidate<'m>> {
        self.step_with(
            target_s,
            endpoint_current_a_m2,
            shaft_torque_nm,
            controls,
            |velocity, checkpoint| {
                self.field
                    .step_with_boundary_candidate(
                        endpoint_current_a_m2,
                        velocity,
                        endpoint_boundary,
                        target_s - self.field.time_s(),
                        field_controls,
                        &mut || checkpoint(),
                    )
                    .map_err(Into::into)
            },
            cancelled,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn step_with(
        &self,
        target_s: f64,
        endpoint_current_a_m2: &[[f64; 3]],
        shaft_torque_nm: f64,
        controls: RotorControls,
        mut advance: impl FnMut(
            &[[f64; 3]],
            &mut dyn FnMut() -> bool,
        ) -> Result<SpatialInductionCandidate<'m>>,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<SpatialRotorCandidate<'m>> {
        poll(cancelled)?;
        let duration = finite(target_s - self.field.time_s())?;
        induction_rotor::admit_controls(duration, shaft_torque_nm, controls)?;
        if endpoint_current_a_m2.len() != self.field.problem().mesh().tets.len()
            || self
                .binding
                .cells
                .iter()
                .any(|&cell| endpoint_current_a_m2[cell] != [0.0; 3])
        {
            return Err(SpatialRotorError::Invalid(
                "complete external-only impressed-current source required",
            ));
        }
        self.binding
            .model
            .validate_applied_effort(&[shaft_torque_nm])?;
        let mut speed = self.angular_velocity_rad_s;
        let mut velocity = vec![[0.0; 3]; self.binding.unit_velocity.len()];
        let mut last_residual = f64::INFINITY;
        for iteration in 1..=controls.maximum_iterations {
            poll(cancelled)?;
            for (output, unit) in velocity.iter_mut().zip(&self.binding.unit_velocity) {
                poll(cancelled)?;
                *output = [finite(speed * unit[0])?, finite(speed * unit[1])?, 0.0];
            }
            let proposal = advance(&velocity, cancelled)?;
            let loads =
                proposal.lorentz_loads(&self.binding.cells, self.binding.load_limits, cancelled)?;
            let stage = induction_rotor::stage(
                &self.binding.model,
                self.angle_rad,
                self.angular_velocity_rad_s,
                speed,
                duration,
                shaft_torque_nm,
                &loads,
                proposal.step().mechanical_work_j,
                proposal.step().energy_residual_j,
                iteration,
                controls,
                cancelled,
            )?;
            last_residual = stage.balance.velocity_residual_rad_s;
            if let Some(endpoint) = stage.endpoint {
                let (field, induction) = proposal.into_parts();
                let state = Self {
                    field,
                    binding: Arc::clone(&self.binding),
                    angle_rad: endpoint.angle_rad,
                    angular_velocity_rad_s: endpoint.angular_velocity_rad_s,
                };
                poll(cancelled)?;
                return Ok(SpatialRotorCandidate {
                    state,
                    induction,
                    mechanics: stage.mechanics,
                    balance: stage.balance,
                });
            }
            let next = finite(
                speed
                    + controls.relaxation
                        * (stage.balance.mechanical_midpoint_velocity_rad_s - speed),
            )?;
            if next == speed {
                return Err(SpatialRotorError::Balance(
                    "coupled midpoint stalled before work/velocity acceptance",
                ));
            }
            speed = next;
        }
        Err(SpatialRotorError::Nonconvergence {
            iterations: controls.maximum_iterations,
            velocity_residual_rad_s: last_residual,
        })
    }
}
