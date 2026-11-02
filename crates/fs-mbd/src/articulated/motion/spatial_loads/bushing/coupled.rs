//! Full passive stiffness/damping on the original finite bushing coordinates.
//!
//! Cross-axis and translation/rotation coupling act on `[q, Log(R)]` and their
//! exact rates. The original geometry, reciprocal wrench lift, and simultaneous
//! implicit mechanics remain authoritative. These are explicitly supplied native
//! coefficients; no missing material tensor, thermal law or source is inferred.

mod matrix;
pub use matrix::{BushingCoordinateScales, BushingMatrixError, PassiveBushingMatrix};

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
    mut cancelled: impl FnMut() -> bool,
) -> Result<SpatialMotionStep, FreeFloatingMotionError> {
    poll(&mut cancelled)?;
    super::admit(model, loads, bushings, &mut cancelled)?;
    let evaluated = step_evaluated(
        model,
        state,
        loads.held,
        duration,
        &mut cancelled,
        |map, base, cancelled| super::evaluate(map, base, loads, bushings, cancelled),
        None,
    )?;
    finish_with_bushings(model, evaluated, loads, bushings, &mut cancelled)
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
