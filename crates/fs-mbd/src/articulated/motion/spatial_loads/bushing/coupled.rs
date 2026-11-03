//! Full passive stiffness/damping on the original finite bushing coordinates.
//!
//! Cross-axis and translation/rotation coupling act on `[q, Log(R)]` and their
//! exact rates. The original geometry, reciprocal wrench lift, and simultaneous
//! implicit mechanics remain authoritative. These are explicitly supplied native
//! coefficients; no missing material tensor, thermal law or source is inferred.

mod matrix;
pub use matrix::{BushingCoordinateScales, BushingMatrixError, PassiveBushingMatrix};

use super::super::super::prescribed::{
    PrescribedBaseMotionError, PrescribedBaseMotionState, PrescribedBasePath,
};
use super::super::prescribed::PrescribedSpatialMotionStep;
use super::super::prescribed::implicit::{
    ImplicitPrescribedSpatialError, ImplicitPrescribedSpatialStep,
};
use super::*;

/// Fixed attachment geometry with independently prepared full passive K and C.
///
/// Constitutive effort is `K [q,r] + C [q_dot,r_dot]`. K stores half its
/// quadratic form and C dissipates its nonnegative quadratic form. Coupling
/// need not share principal axes. The native inverse-log Jacobian and A-side
/// translation couple preserve the original finite-rotation work conjugacy.
#[derive(Debug, Clone, PartialEq)]
pub struct CoupledSpatialBushing {
    a: MaterialFrame,
    b: MaterialFrame,
    rest_translation_a: Vec3,
    stiffness: PassiveBushingMatrix,
    damping: PassiveBushingMatrix,
}

impl CoupledSpatialBushing {
    /// Bind prepared native stiffness and damping to two fixed material frames.
    /// Matrix units and coordinate order are those of [`PassiveBushingMatrix`].
    /// Their numerical coordinate scales may differ; physical SI coordinates
    /// and the original supplied matrices remain unchanged.
    pub fn try_new(
        a: MaterialFrame,
        b: MaterialFrame,
        rest_translation_a: Vec3,
        stiffness: PassiveBushingMatrix,
        damping: PassiveBushingMatrix,
    ) -> Result<Self, FreeFloatingMotionError> {
        if a.link == b.link {
            return Err(FreeFloatingMotionError::InvalidLoad(
                "bushing attachment links",
            ));
        }
        validate_vec3(rest_translation_a, "bushing.rest_translation")?;
        Ok(Self {
            a,
            b,
            rest_translation_a,
            stiffness,
            damping,
        })
    }

    /// Represent an existing diagonal bushing in the full-law inventory.
    /// Exact zero directions remain zero; direct diagonal action is retained.
    pub fn try_from_diagonal(
        b: SpatialBushing,
        scales: BushingCoordinateScales,
    ) -> Result<Self, FreeFloatingMotionError> {
        let coefficients = |translation: Vec3, rotation: Vec3| {
            let values = coordinates(translation, rotation);
            core::array::from_fn(|i| core::array::from_fn(|j| if i == j { values[i] } else { 0.0 }))
        };
        let stiffness = PassiveBushingMatrix::try_new(
            coefficients(b.translation_stiffness, b.rotation_stiffness),
            scales,
        )
        .map_err(matrix_error)?;
        let damping = PassiveBushingMatrix::try_new(
            coefficients(b.translation_damping, b.rotation_damping),
            scales,
        )
        .map_err(matrix_error)?;
        Self::try_new(b.a, b.b, b.rest_translation_a, stiffness, damping)
    }

    pub const fn frames(&self) -> [MaterialFrame; 2] {
        [self.a, self.b]
    }
    pub const fn rest_translation_a(&self) -> Vec3 {
        self.rest_translation_a
    }
    pub const fn stiffness(&self) -> &PassiveBushingMatrix {
        &self.stiffness
    }
    pub const fn damping(&self) -> &PassiveBushingMatrix {
        &self.damping
    }
}

/// Evaluate the full passive laws using the original shape-relative geometry.
/// Work, storage and loss have the same native accounts as diagonal bushings.
pub fn evaluate_bushing_forces(
    model: &ArticulatedModel,
    state: &FreeFloatingMotionState,
    bushings: &[CoupledSpatialBushing],
    cancelled: impl FnMut() -> bool,
) -> Result<BushingForces, FreeFloatingMotionError> {
    super::evaluate_forces(model, state, bushings, cancelled)
}

/// Advance all full-matrix bushings in the original explicit Lie midpoint step.
/// The original timestep, chart, cancellation and work-defect semantics apply.
pub fn step_free_floating_bushings(
    model: &ArticulatedModel,
    state: &FreeFloatingMotionState,
    loads: SpatialMotionLoads<'_>,
    bushings: &[CoupledSpatialBushing],
    duration: f64,
    cancelled: impl FnMut() -> bool,
) -> Result<SpatialMotionStep, FreeFloatingMotionError> {
    super::step_with_bushings_and_effort(model, state, loads, bushings, duration, None, cancelled)
}

/// Advance full passive bushings and internal actuator feedback together.
///
/// The callback receives the complete actual free-carrier stage at elapsed
/// times zero and duration/2, with effort reset to the held values. Full
/// constitutive coupling, finite-body recoil and feedback use that same stage.
/// Actuator work remains separate from bushing storage and physical damping.
///
/// The original rate/effort limits, cancellation and chart admission apply.
/// Resolve the fastest attachment/controller mode with this explicit step.
/// Callback side effects remain caller-owned on any refusal.
#[allow(clippy::too_many_arguments)]
pub fn step_free_floating_bushings_with_effort(
    model: &ArticulatedModel,
    state: &FreeFloatingMotionState,
    loads: SpatialMotionLoads<'_>,
    bushings: &[CoupledSpatialBushing],
    duration: f64,
    mut effort: impl FnMut(
        f64,
        &FreeFloatingMotionState,
        &mut [f64],
    ) -> Result<(), FreeFloatingMotionError>,
    cancelled: impl FnMut() -> bool,
) -> Result<SpatialMotionStep, FreeFloatingMotionError> {
    super::step_with_bushings_and_effort(
        model,
        state,
        loads,
        bushings,
        duration,
        Some(&mut effort),
        cancelled,
    )
}

/// Advance full K/C coupling through the SAME native simultaneous implicit solve.
/// This supports stiff cross-coupled modes, base recoil, external loads and
/// axial connections, with their original work and physical-loss accounting.
/// Implicit midpoint remains non-L-stable; refine time for trajectory accuracy.
pub fn step_free_floating_bushings_implicit(
    model: &ArticulatedModel,
    state: &FreeFloatingMotionState,
    loads: SpatialMotionLoads<'_>,
    bushings: &[CoupledSpatialBushing],
    duration: f64,
    config: ImplicitSpatialConfig,
    cancelled: impl FnMut() -> bool,
) -> Result<ImplicitSpatialStep, ImplicitSpatialError> {
    super::super::implicit::step_with_bushings(
        model, state, loads, bushings, duration, config, cancelled,
    )
}

/// Solve full passive bushing coupling and smooth feedback at one midpoint.
///
/// Uses the existing native nonlinear owner, including physical base recoil,
/// all spatial loads and the independent final residual check. The callback
/// receives elapsed time duration/2 and the actual trial state; its effort
/// buffer is reset to held values before each call. Numerical derivative and
/// nonlinear trials may repeat in nonmonotone order, so the law must be smooth
/// and deterministic. Callback-owned effects cannot be rolled back.
///
/// Physical actuator work, bushing storage and damping retain separate accounts.
/// The original bounded solve, cancellation and geometric refusals apply.
/// Implicit midpoint is not L-stable; refine for phase and nonlinear geometry.
#[allow(clippy::too_many_arguments)]
pub fn step_free_floating_bushings_implicit_with_effort(
    model: &ArticulatedModel,
    state: &FreeFloatingMotionState,
    loads: SpatialMotionLoads<'_>,
    bushings: &[CoupledSpatialBushing],
    duration: f64,
    config: ImplicitSpatialConfig,
    mut effort: impl FnMut(
        f64,
        &FreeFloatingMotionState,
        &mut [f64],
    ) -> Result<(), FreeFloatingMotionError>,
    cancelled: impl FnMut() -> bool,
) -> Result<ImplicitSpatialStep, ImplicitSpatialError> {
    super::super::implicit::step_with_bushings_and_effort(
        model,
        state,
        loads,
        bushings,
        duration,
        config,
        Some(&mut effort),
        cancelled,
    )
}

/// Advance full passive K/C coupling on a fixed or prescribed moving support.
///
/// Uses the same finite geometry, reciprocal wrenches and native supported
/// midpoint lifecycle as diagonal bushings. All cross-axis forces enter both
/// joint dynamics and the complete support reaction; storage, damping loss
/// and signed support work retain their separate physical accounts.
///
/// This explicit step must resolve the fastest coupled attachment mode.
/// The caller supplies a consistent support path; chart admission and
/// cancellation have the same semantics as the original supported entry.
#[allow(clippy::too_many_arguments)]
pub fn step_prescribed_base_bushings(
    model: &ArticulatedModel,
    state: &PrescribedBaseMotionState,
    path: PrescribedBasePath,
    loads: SpatialMotionLoads<'_>,
    bushings: &[CoupledSpatialBushing],
    duration: f64,
    cancelled: impl FnMut() -> bool,
) -> Result<PrescribedSpatialMotionStep, PrescribedBaseMotionError> {
    super::super::prescribed::step_with_bushings_and_effort(
        model, state, path, loads, bushings, duration, None, cancelled,
    )
}

/// Advance supported full-matrix bushings with stage-dependent joint effort.
///
/// The callback receives the actual prescribed-base stage and elapsed time.
/// Its output starts from held effort at each stage. Feedback and all full
/// bushing couplings act on that same state, including support recovery.
/// Native inputs remain borrowed; callback-owned effects cannot be rolled back.
#[allow(clippy::too_many_arguments)]
pub fn step_prescribed_base_bushings_with_effort(
    model: &ArticulatedModel,
    state: &PrescribedBaseMotionState,
    path: PrescribedBasePath,
    loads: SpatialMotionLoads<'_>,
    bushings: &[CoupledSpatialBushing],
    duration: f64,
    mut effort: impl FnMut(
        f64,
        &PrescribedBaseMotionState,
        &mut [f64],
    ) -> Result<(), PrescribedBaseMotionError>,
    cancelled: impl FnMut() -> bool,
) -> Result<PrescribedSpatialMotionStep, PrescribedBaseMotionError> {
    super::super::prescribed::step_with_bushings_and_effort(
        model,
        state,
        path,
        loads,
        bushings,
        duration,
        Some(&mut effort),
        cancelled,
    )
}

/// Advance supported full K/C coupling with the native implicit midpoint solve.
///
/// All spatial loads, bushings and support dynamics use one simultaneous stage.
/// The nonlinear unknowns are original joint rates; native ABA and RNEA retain
/// the full mass distribution and prescribed support impulse/work. Implicit
/// midpoint is A-stable for a linear oscillator but not L-stable; resolve phase
/// and nonlinear geometry by timestep refinement.
#[allow(clippy::too_many_arguments)]
pub fn step_prescribed_base_bushings_implicit(
    model: &ArticulatedModel,
    state: &PrescribedBaseMotionState,
    path: PrescribedBasePath,
    loads: SpatialMotionLoads<'_>,
    bushings: &[CoupledSpatialBushing],
    duration: f64,
    config: ImplicitSpatialConfig,
    cancelled: impl FnMut() -> bool,
) -> Result<ImplicitPrescribedSpatialStep, ImplicitPrescribedSpatialError> {
    super::super::prescribed::implicit::step_with_bushings_and_effort(
        model, state, path, loads, bushings, duration, config, None, cancelled,
    )
}

/// Solve smooth actuator feedback and supported full-matrix bushings together.
///
/// Feedback receives each actual trial midpoint, including numerical derivative
/// stencils and independent final verification. The effort buffer is reset to
/// held effort before each call. The law must be deterministic and smooth over
/// admitted trials; callback-owned side effects cannot be rolled back.
#[allow(clippy::too_many_arguments)]
pub fn step_prescribed_base_bushings_implicit_with_effort(
    model: &ArticulatedModel,
    state: &PrescribedBaseMotionState,
    path: PrescribedBasePath,
    loads: SpatialMotionLoads<'_>,
    bushings: &[CoupledSpatialBushing],
    duration: f64,
    config: ImplicitSpatialConfig,
    mut effort: impl FnMut(
        f64,
        &PrescribedBaseMotionState,
        &mut [f64],
    ) -> Result<(), PrescribedBaseMotionError>,
    cancelled: impl FnMut() -> bool,
) -> Result<ImplicitPrescribedSpatialStep, ImplicitPrescribedSpatialError> {
    super::super::prescribed::implicit::step_with_bushings_and_effort(
        model,
        state,
        path,
        loads,
        bushings,
        duration,
        config,
        Some(&mut effort),
        cancelled,
    )
}

impl BushingLaw for CoupledSpatialBushing {
    fn geometry(&self) -> BushingGeometry {
        BushingGeometry {
            a: self.a,
            b: self.b,
            rest_translation_a: self.rest_translation_a,
            rotation_active: self.stiffness.rotation_active() || self.damping.rotation_active(),
        }
    }
    fn admit_coefficients(&self) -> Result<(), FreeFloatingMotionError> {
        Ok(())
    }
    fn response(
        &self,
        s: &Strain,
        translation_rate: Vec3,
        rotation_rate: Vec3,
    ) -> Result<ConstitutiveResponse, FreeFloatingMotionError> {
        let (elastic, energy) = self
            .stiffness
            .evaluate(coordinates(s.translation, s.rotation))
            .map_err(matrix_error)?;
        let (viscous, dissipation_power) = self
            .damping
            .evaluate(coordinates(translation_rate, rotation_rate))
            .map_err(matrix_error)?;
        Ok(ConstitutiveResponse {
            force: Vec3::new(
                elastic[0] + viscous[0],
                elastic[1] + viscous[1],
                elastic[2] + viscous[2],
            ),
            log_force: Vec3::new(
                elastic[3] + viscous[3],
                elastic[4] + viscous[4],
                elastic[5] + viscous[5],
            ),
            stored_energy: 0.5 * energy,
            dissipation_power,
        })
    }
    fn energy(&self, s: &Strain) -> Result<f64, FreeFloatingMotionError> {
        self.stiffness
            .evaluate(coordinates(s.translation, s.rotation))
            .map(|(_, energy)| 0.5 * energy)
            .map_err(matrix_error)
    }
}

fn coordinates(q: Vec3, r: Vec3) -> [f64; 6] {
    [q.x, q.y, q.z, r.x, r.y, r.z]
}
fn matrix_error(error: BushingMatrixError) -> FreeFloatingMotionError {
    FreeFloatingMotionError::InvalidLoad(error.0)
}
