//! Objective transfer of the original cohesive law to rotating material frames.
//!
//! The frame maps local `[normal, shear_1, shear_2]` to world coordinates. With
//! attachment separation `d = x_b - x_a` and stress-free local offset `r`, the
//! physical jump is `g = R^T d - r`. Common rigid transport changes neither this
//! jump nor its original irreversible/frozen-BK history.
//!
//! For a physical trial, the original law supplies its exact straight-local-path
//! mean traction. The discrete pullback uses `Rbar = (R0 + R1)/2`, not a normalized
//! midpoint rotation. The spatial Cayley increment satisfies
//! `R1 - R0 = [c]x Rbar`. Thus `F = A Rbar tmean`, `M = dbar x F` obey
//! `F . (d1-d0) - M . c = A tmean . (g1-g0)`. Forces on attachments a/b are
//! `+F/-F`; **the frame-carrying side receives the additional couple `+M`**.
//! This cancels the force pair's midpoint moment, including a finite bond offset.
//!
//! These are discrete work efforts, not temporal averages of world traction.
//! A mechanics owner must use the same Cayley increment and midpoint attachment
//! moment arms, or independently audit its actual work. In particular `c` is
//! not the exponential rotation vector. This adapter does not make an arbitrary
//! rigid-body integrator energy conserving, choose a shared frame/torque split,
//! certify intrastep trajectory accuracy, or resolve a half-turn in one step.

use fs_material::cohesive::{CohesiveError, CohesivePoint, CohesiveRateStep};
use fs_mbd::{UnitQuaternion, Vec3};

/// Original interface-source retention around the numerical transport owner.
pub mod source;

/// Explicit rotation-chart and independent kinematic/work closure budgets.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ObjectiveCohesiveConfig {
    /// Maximum norm of the spatial Cayley increment `2 tan(angle/2)`, positive.
    /// Subdivide a physical step that exceeds this bound or reaches a half-turn.
    pub maximum_cayley_increment: f64,
    /// Absolute componentwise discrete kinematic closure tolerance [m].
    pub kinematic_tolerance_m: f64,
    /// Original material owner's absolute local path-work tolerance [J/m2].
    pub material_energy_tolerance_j_m2: f64,
    /// Independent absolute spatial-work closure tolerance [J].
    pub spatial_work_tolerance_j: f64,
}

/// A refusal leaves all accepted material and geometric state unchanged.
#[derive(Clone, Debug, PartialEq)]
pub enum ObjectiveCohesiveError {
    Invalid(&'static str),
    Material(CohesiveError),
    RotationLimit { cayley_norm: f64, maximum: f64 },
    KinematicDefect { residual_m: Vec3, tolerance_m: f64 },
    WorkDefect { defect_j: f64, tolerance_j: f64 },
}
impl From<CohesiveError> for ObjectiveCohesiveError {
    fn from(value: CohesiveError) -> Self {
        Self::Material(value)
    }
}
impl core::fmt::Display for ObjectiveCohesiveError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "objective cohesive trial refused: {self:?}")
    }
}
impl std::error::Error for ObjectiveCohesiveError {}

/// Attachment force pair and the balancing couple on the frame-carrying side.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ObjectiveCohesiveWrench {
    pub force_on_a_n: Vec3,
    pub force_on_b_n: Vec3,
    pub frame_couple_n_m: Vec3,
}

/// Original law/history with explicit finite geometry, all in SI units.
///
/// The reference offset is fixed in the material frame. Area is the original
/// bonded area, not the current projected/contact area. Construction does not
/// reset an existing damaged point or infer a new bond from geometric contact.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ObjectiveCohesivePoint<'a> {
    material: CohesivePoint<'a>,
    frame: UnitQuaternion,
    reference_offset_m: Vec3,
    separation_world_m: Vec3,
    area_m2: f64,
}

/// Immutable material successor and the compatible spatial work ports.
#[derive(Clone, Debug, PartialEq)]
pub struct ObjectiveCohesiveStep<'a> {
    pub point: ObjectiveCohesivePoint<'a>,
    /// Original intrinsic storage/fracture account and separate viscous loss.
    pub material: CohesiveRateStep<'a>,
    pub mean_traction_pa: [f64; 3],
    /// Work-conjugate discrete efforts, with midpoint attachment moment arms.
    pub discrete_wrench: ObjectiveCohesiveWrench,
    /// Actual endpoint traction/geometry; never substituted for discrete work.
    pub endpoint_wrench: ObjectiveCohesiveWrench,
    pub cayley_increment_world: Vec3,
    pub kinematic_residual_m: Vec3,
    /// Work supplied to the bond, the negative of the internal body/frame work.
    pub spatial_work_j: f64,
    pub work_defect_j: f64,
}

impl<'a> ObjectiveCohesivePoint<'a> {
    /// Attach an original accepted material point to its current material frame.
    /// The actual separation is `R * (reference_offset + original_local_jump)`.
    pub fn new(
        material: CohesivePoint<'a>,
        frame: UnitQuaternion,
        reference_offset_m: Vec3,
        area_m2: f64,
    ) -> Result<Self, ObjectiveCohesiveError> {
        if !area_m2.is_finite() || area_m2 <= 0.0 || !reference_offset_m.is_finite() {
            return Err(ObjectiveCohesiveError::Invalid(
                "finite offset and positive original area required",
            ));
        }
        let separation_world_m =
            frame.rotate_body_to_world(reference_offset_m.add(vector(material.jump_m())));
        finite_vector(separation_world_m)?;
        finite_scalar(area_m2 * material.stored_energy_j_m2())?;
        Ok(Self {
            material,
            frame,
            reference_offset_m,
            separation_world_m,
            area_m2,
        })
    }

    pub const fn material(&self) -> CohesivePoint<'a> {
        self.material
    }
    pub const fn frame(&self) -> UnitQuaternion {
        self.frame
    }
    pub const fn reference_offset_m(&self) -> Vec3 {
        self.reference_offset_m
    }
    pub const fn separation_world_m(&self) -> Vec3 {
        self.separation_world_m
    }
    pub const fn area_m2(&self) -> f64 {
        self.area_m2
    }

    /// Explicit common rigid transport, with no physical local deformation.
    ///
    /// Keeps the complete original material point bit-for-bit and rotates its
    /// geometry. This is not a mechanical timestep or impulse; it does not
    /// silently clamp a supplied strain, reproject history, evaluate viscosity,
    /// or turn a geometric coincidence into rebonding. Arbitrary finite frame
    /// changes are allowed here because no Cayley work step is being taken.
    pub fn transport_frame(self, frame: UnitQuaternion) -> Result<Self, ObjectiveCohesiveError> {
        let separation_world_m =
            frame.rotate_body_to_world(self.reference_offset_m.add(vector(self.material.jump_m())));
        finite_vector(separation_world_m)?;
        Ok(Self {
            frame,
            separation_world_m,
            ..self
        })
    }

    /// Trial from an actual new attachment separation in world coordinates.
    /// The new local physical jump is resolved from that geometry without any
    /// tolerance-based clamping. Use `transport_frame` for declared zero-strain
    /// transport, or `trial_local_jump` when the physical local jump is known.
    pub fn trial(
        &self,
        frame: UnitQuaternion,
        separation_world_m: Vec3,
        dt_s: f64,
        config: ObjectiveCohesiveConfig,
    ) -> Result<ObjectiveCohesiveStep<'a>, ObjectiveCohesiveError> {
        finite_vector(separation_world_m)?;
        let jump = frame
            .rotate_world_to_body(separation_world_m)
            .sub(self.reference_offset_m);
        self.trial_resolved(frame, separation_world_m, array(jump), dt_s, config)
    }

    /// Trial a declared straight physical local jump path while its material
    /// frame rotates. Keeping the same jump gives exactly zero local rate and
    /// damage increment, even when a stressed/viscous bond undergoes rotation.
    pub fn trial_local_jump(
        &self,
        frame: UnitQuaternion,
        jump_m: [f64; 3],
        dt_s: f64,
        config: ObjectiveCohesiveConfig,
    ) -> Result<ObjectiveCohesiveStep<'a>, ObjectiveCohesiveError> {
        let separation_world_m =
            frame.rotate_body_to_world(self.reference_offset_m.add(vector(jump_m)));
        finite_vector(separation_world_m)?;
        self.trial_resolved(frame, separation_world_m, jump_m, dt_s, config)
    }

    fn trial_resolved(
        &self,
        frame: UnitQuaternion,
        separation_world_m: Vec3,
        jump_m: [f64; 3],
        dt_s: f64,
        config: ObjectiveCohesiveConfig,
    ) -> Result<ObjectiveCohesiveStep<'a>, ObjectiveCohesiveError> {
        validate_config(config)?;
        let cayley = cayley_increment(self.frame, frame, config.maximum_cayley_increment)?;
        let delta_d = separation_world_m.sub(self.separation_world_m);
        let middle_d = self
            .separation_world_m
            .scale(0.5)
            .add(separation_world_m.scale(0.5));
        let co_rotating_delta = delta_d.sub(cayley.cross(middle_d));
        let geometric_delta = self
            .frame
            .rotate_world_to_body(co_rotating_delta)
            .scale(0.5)
            .add(frame.rotate_world_to_body(co_rotating_delta).scale(0.5));
        let residual = vector(jump_m)
            .sub(vector(self.material.jump_m()))
            .sub(geometric_delta);
        finite_vector(residual)?;
        if residual.x.abs().max(residual.y.abs()).max(residual.z.abs())
            > config.kinematic_tolerance_m
        {
            return Err(ObjectiveCohesiveError::KinematicDefect {
                residual_m: residual,
                tolerance_m: config.kinematic_tolerance_m,
            });
        }
        let (material, mean_traction_pa) = self.material.trial_over_time_with_mean_traction(
            jump_m,
            dt_s,
            config.material_energy_tolerance_j_m2,
        )?;
        let local_force = vector(mean_traction_pa).scale(self.area_m2);
        let force = self
            .frame
            .rotate_body_to_world(local_force)
            .scale(0.5)
            .add(frame.rotate_body_to_world(local_force).scale(0.5));
        let discrete_wrench = wrench(force, middle_d)?;
        let endpoint_force =
            frame.rotate_body_to_world(vector(material.traction_pa).scale(self.area_m2));
        let endpoint_wrench = wrench(endpoint_force, separation_world_m)?;
        let spatial_work_j =
            finite_scalar(force.dot(delta_d) - discrete_wrench.frame_couple_n_m.dot(cayley))?;
        let work_defect_j =
            finite_scalar(spatial_work_j - self.area_m2 * material.total_work_j_m2())?;
        if work_defect_j.abs() > config.spatial_work_tolerance_j {
            return Err(ObjectiveCohesiveError::WorkDefect {
                defect_j: work_defect_j,
                tolerance_j: config.spatial_work_tolerance_j,
            });
        }
        Ok(ObjectiveCohesiveStep {
            point: Self {
                material: material.equilibrium.point,
                frame,
                separation_world_m,
                ..*self
            },
            material,
            mean_traction_pa,
            discrete_wrench,
            endpoint_wrench,
            cayley_increment_world: cayley,
            kinematic_residual_m: residual,
            spatial_work_j,
            work_defect_j,
        })
    }
}

fn validate_config(config: ObjectiveCohesiveConfig) -> Result<(), ObjectiveCohesiveError> {
    if [
        config.maximum_cayley_increment,
        config.kinematic_tolerance_m,
        config.material_energy_tolerance_j_m2,
        config.spatial_work_tolerance_j,
    ]
    .into_iter()
    .any(|value| !value.is_finite() || value <= 0.0)
    {
        return Err(ObjectiveCohesiveError::Invalid(
            "positive finite rotation, kinematic and work budgets required",
        ));
    }
    Ok(())
}

pub(super) fn cayley_increment(
    old: UnitQuaternion,
    new: UnitQuaternion,
    maximum: f64,
) -> Result<Vec3, ObjectiveCohesiveError> {
    let [ow, ox, oy, oz] = old.components();
    let [nw, nx, ny, nz] = new.components();
    // Vector/scalar of q_new * conjugate(q_old); no second quaternion owner.
    // Differences make an identical frame produce bit-exact zero increment.
    let ov = Vec3::new(ox, oy, oz);
    let dv = Vec3::new(nx - ox, ny - oy, nz - oz);
    let v = dv.scale(ow).sub(ov.scale(nw - ow)).add(ov.cross(dv));
    let w = ow.mul_add(nw, ox.mul_add(nx, oy.mul_add(ny, oz * nz)));
    let c = v.scale(2.0 / w);
    let norm = c.x.hypot(c.y).hypot(c.z);
    if !c.is_finite() || !norm.is_finite() || norm > maximum {
        return Err(ObjectiveCohesiveError::RotationLimit {
            cayley_norm: norm,
            maximum,
        });
    }
    Ok(c)
}

fn wrench(
    force: Vec3,
    separation: Vec3,
) -> Result<ObjectiveCohesiveWrench, ObjectiveCohesiveError> {
    finite_vector(force)?;
    let frame_couple_n_m = separation.cross(force);
    finite_vector(frame_couple_n_m)?;
    Ok(ObjectiveCohesiveWrench {
        force_on_a_n: force,
        force_on_b_n: force.scale(-1.0),
        frame_couple_n_m,
    })
}
fn finite_scalar(value: f64) -> Result<f64, ObjectiveCohesiveError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(ObjectiveCohesiveError::Invalid(
            "non-finite cohesive geometry/work",
        ))
    }
}
fn finite_vector(value: Vec3) -> Result<(), ObjectiveCohesiveError> {
    if value.is_finite() {
        Ok(())
    } else {
        Err(ObjectiveCohesiveError::Invalid(
            "non-finite cohesive geometry/effort",
        ))
    }
}
fn vector(value: [f64; 3]) -> Vec3 {
    Vec3::new(value[0], value[1], value[2])
}
fn array(value: Vec3) -> [f64; 3] {
    [value.x, value.y, value.z]
}

#[cfg(test)]
mod tests;
