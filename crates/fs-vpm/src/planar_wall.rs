//! A stationary impermeable slip wall for native two-dimensional vortex motion.
//!
//! Each source has an opposite-circulation image reflected across one infinite
//! line. The same algebraic Biot–Savart core is used on both sides, so their
//! normal velocities cancel at the wall. Images are evaluated from the actual
//! source positions at every RK4 stage, including a vortex's own image. They
//! never become physical particles or enter the returned circulation inventory.
//!
//! This is an inviscid slip boundary: tangential velocity remains free. It does
//! not supply no-slip vorticity generation, viscous wall diffusion, wall motion,
//! corners or a mechanical wall-load model.
//! The algebraic core remains numerical regularization; core overlap with the
//! wall does not resolve a viscous boundary layer or certify fluid circulation
//! from the particle-strength inventory alone.

use crate::{
    VelocityKernel, VortexParticle, VpmBudget, VpmError, VpmRun, checked_induced_contribution,
    invalid_scalar, validate_particle,
};
use fs_exec::Cx;

/// One fixed planar wall in physical coordinates, with the fluid on its
/// positive-normal side. Distances are metres when particle coordinates are.
/// The normal is normalized once during construction; geometry is immutable.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PlanarWall2 {
    origin_m: [f64; 2],
    normal: [f64; 2],
}

impl PlanarWall2 {
    /// Construct a wall through `origin_m`, with a nonzero normal pointing into
    /// the fluid. The input normal need not have unit length.
    ///
    /// # Errors
    /// Refuses nonfinite coordinates or a zero/nonfinite normal.
    pub fn new(origin_m: [f64; 2], normal_into_fluid: [f64; 2]) -> Result<Self, VpmError> {
        for (name, values) in [
            ("wall origin", origin_m),
            ("wall normal", normal_into_fluid),
        ] {
            for value in values {
                if !value.is_finite() {
                    return Err(invalid_scalar(name, value, None));
                }
            }
        }
        let scale = normal_into_fluid[0].abs().max(normal_into_fluid[1].abs());
        if scale == 0.0 {
            return Err(invalid_scalar("wall normal length", scale, None));
        }
        let scaled = normal_into_fluid.map(|value| value / scale);
        let length = scaled[0].hypot(scaled[1]);
        Ok(Self {
            origin_m,
            normal: scaled.map(|value| value / length),
        })
    }

    /// A point on the wall in the particle coordinate system.
    #[must_use]
    pub const fn origin_m(self) -> [f64; 2] {
        self.origin_m
    }

    /// Unit normal pointing into the admitted fluid half-plane.
    #[must_use]
    pub const fn normal_into_fluid(self) -> [f64; 2] {
        self.normal
    }

    /// Evaluate the native real-plus-image velocity, without allocating images.
    /// Sources must lie strictly in the fluid; the target can be on the wall.
    /// The formula also defines its reflected extension outside the fluid.
    /// As in the original kernel, an exactly coincident real source contributes
    /// zero velocity; its distinct wall image still contributes.
    ///
    /// # Errors
    /// Refuses nonfinite inputs, negative core, a source at/beyond the wall,
    /// unrepresentable reflected geometry or a nonfinite velocity sum.
    pub fn induced_velocity(
        self,
        particles: &[VortexParticle],
        point: [f64; 2],
        core: f64,
    ) -> Result<[f64; 2], VpmError> {
        if !core.is_finite() || core < 0.0 {
            return Err(invalid_scalar("core", core, None));
        }
        for coordinate in point {
            if !coordinate.is_finite() {
                return Err(invalid_scalar("wall sample position", coordinate, None));
            }
        }
        let kernel = VelocityKernel::new(core);
        let mut velocity = [0.0; 2];
        const PHASE: &str = "fs-vpm.planar-wall-sample";
        // Keep physical sources first, matching the native time-step sweep.
        for (index, particle) in particles.iter().enumerate() {
            validate_particle(*particle, index)?;
            self.admit_position(particle.pos, index)?;
            checked_induced_contribution(&mut velocity, particle, point, kernel, PHASE, index)?;
        }
        for (index, particle) in particles.iter().enumerate() {
            checked_induced_contribution(
                &mut velocity,
                &self.image(particle),
                point,
                kernel,
                PHASE,
                index,
            )?;
        }
        Ok(velocity)
    }

    /// Sample the same real-plus-image velocity, full spatial Jacobian and
    /// analytic curl used by this wall's flow. Sources must lie strictly in the
    /// fluid; the probe may lie on the wall or in the reflected extension.
    ///
    /// Images are opposite-circulation boundary devices, never physical
    /// particles. Their strain and axial vorticity are included. Curl is summed
    /// from the original analytic kernel, retaining thin-core tails that a
    /// subtraction of rounded Jacobian entries could lose.
    ///
    /// This query uses O(N) work and O(1) stack storage with no image cloud.
    /// Like the original per-probe velocity query it has no execution receipt;
    /// callers must bound bulk sampling. It supplies a slip-wall field, not
    /// no-slip vorticity generation or viscous wall diffusion.
    ///
    /// # Errors
    /// Original scalar/clearance/reflection refusals, or a singular zero-core
    /// probe, unrepresentable source derivative or nonfinite accumulated field.
    /// A positive-core source center remains regular.
    pub fn induced_kinematics(
        self,
        particles: &[VortexParticle],
        point: [f64; 2],
        core: f64,
    ) -> Result<crate::ParticleKinematics2, VpmError> {
        if !core.is_finite() || core < 0.0 {
            return Err(invalid_scalar("core", core, None));
        }
        for coordinate in point {
            if !coordinate.is_finite() {
                return Err(invalid_scalar("wall sample position", coordinate, None));
            }
        }
        const PHASE: &str = "fs-vpm.planar-wall-kinematics";
        let mut result = crate::ParticleKinematics2::default();
        let mut add = |source: &VortexParticle, index| -> Result<(), VpmError> {
            let sample = crate::induced_kinematics(core::slice::from_ref(source), point, core)
                .ok_or(VpmError::NonFiniteResult {
                    phase: PHASE,
                    particle: index,
                })?;
            for i in 0..2 {
                result.velocity[i] += sample.velocity[i];
            }
            result.gradient[0][0] += sample.gradient[0][0];
            result.gradient[0][1] += sample.gradient[0][1];
            result.gradient[1][0] += sample.gradient[1][0];
            result.gradient[1][1] = -result.gradient[0][0];
            result.vorticity += sample.vorticity;
            if result
                .velocity
                .iter()
                .chain(result.gradient.iter().flatten())
                .chain(core::iter::once(&result.vorticity))
                .any(|value| !value.is_finite())
            {
                return Err(VpmError::NonFiniteResult {
                    phase: PHASE,
                    particle: index,
                });
            }
            Ok(())
        };
        // Match the original wall velocity's physical-first source ordering.
        for (index, source) in particles.iter().enumerate() {
            validate_particle(*source, index)?;
            self.admit_position(source.pos, index)?;
            add(source, index)?;
        }
        for (index, source) in particles.iter().enumerate() {
            add(&self.image(source), index)?;
        }
        Ok(result)
    }

    pub(super) fn admit_position(self, point: [f64; 2], particle: usize) -> Result<(), VpmError> {
        let (distance, image) = self.reflection(point);
        if !distance.is_finite() || image.iter().any(|value| !value.is_finite()) {
            return Err(VpmError::NonFiniteResult {
                phase: "fs-vpm.planar-wall-geometry",
                particle,
            });
        }
        if distance <= 0.0 {
            return Err(invalid_scalar("wall clearance", distance, Some(particle)));
        }
        Ok(())
    }

    // The original input/update units admit every physical position before
    // a velocity sweep. Images are constant-size stack values, not another
    // cloud, and repeated pair evaluation needs no second geometry admission.
    pub(super) fn image(self, particle: &VortexParticle) -> VortexParticle {
        VortexParticle::new(self.reflection(particle.pos).1, -particle.circulation)
    }

    fn reflection(self, point: [f64; 2]) -> (f64, [f64; 2]) {
        let mut distance = 0.0;
        for (axis, &coordinate) in point.iter().enumerate() {
            // A tangential wall-anchor coordinate is irrelevant; do not form
            // a potentially overflowing difference only to multiply by zero.
            if self.normal[axis] != 0.0 {
                distance += (coordinate - self.origin_m[axis]) * self.normal[axis];
            }
        }
        let image = core::array::from_fn(|axis| {
            let projection = distance * self.normal[axis];
            // Do not overflow 2*distance when the reflected coordinate fits.
            (point[axis] - projection) - projection
        });
        (distance, image)
    }
}

/// Evolve the physical vortices beside one stationary planar slip wall using
/// the original checked RK4 transaction. Every source, including the target
/// vortex itself, supplies a reflected opposite-circulation image at each stage.
///
/// All physical input, intermediate and endpoint positions must remain strictly
/// on the fluid side. A crossing refuses the complete run; it is not clipped or
/// repaired. Circulations retain their original bits. Wall-induced transport
/// can change vortex impulse, so free-space centroid/angular-impulse invariants
/// are not asserted. Negative time remains permitted for this inviscid path.
///
/// Images add `4*N*N*steps` pair/work units to the original RK4 plan, without
/// any additional particle or velocity buffers. The same local/ambient budgets,
/// fallible allocations and bounded checkpoints govern all image contributions.
///
/// # Errors
/// Original numerical/resource/execution refusals and invalid wall clearance
/// or reflected geometry at any physical stage. No partial cloud is returned.
pub fn simulate_with_cx(
    cx: &Cx<'_>,
    particles: &[VortexParticle],
    dt: f64,
    steps: usize,
    core: f64,
    wall: PlanarWall2,
    budget: VpmBudget,
) -> Result<VpmRun, VpmError> {
    crate::simulate_impl(cx, particles, dt, steps, core, budget, None, Some(wall))
}
