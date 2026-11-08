//! Local restoring stiffness of the original material/reservoir loads.

use fs_la::eigen::{SymmetricEigenvalues3Error, symmetric_eigenvalues3};
use fs_mbd::{RigidBodyState, Vec3};

use super::pose::HydrostaticPoseAxes;
use super::{
    DisplacementMesh, HydrostaticError, HydrostaticLoad, HydrostaticReservoir, finite, poll,
};
use crate::principal::PrincipalFrame;

/// Controls for a local, quasistatic three-degree-of-freedom analysis.
#[derive(Clone, Copy, Debug)]
pub struct HydrostaticStabilityControls {
    /// Projected perpendicular to gravity to define roll; pitch = up x roll.
    pub roll_direction_world: Vec3,
    /// Positive L in q = (heave/L, roll, pitch). All stiffness entries have units J.
    pub length_scale_m: f64,
    /// Positive physical heave perturbation; both this step and half are used.
    pub heave_step_m: f64,
    /// Angular perturbation in (0,0.25] radians; both this step and half are used.
    pub tilt_step_radians: f64,
    /// Maximum norm of actual buoyancy plus body weight before differentiation.
    pub maximum_force_defect_n: f64,
    /// Maximum norm of actual torque about the material COM before differentiation.
    pub maximum_torque_defect_nm: f64,
    /// Refinement, one-sided derivative jump and symmetry defects must each be no greater than
    /// neutral_stiffness_tolerance_j + this fraction times the largest matrix entry.
    /// In (0,1). Large defects refuse the analysis instead of silently symmetrizing.
    pub maximum_relative_stiffness_defect: f64,
    /// Positive absolute stiffness resolution. This is also the minimum neutral
    /// classification margin; measured numerical defects increase that margin.
    pub neutral_stiffness_tolerance_j: f64,
}

/// Numerical local stiffness classification, excluding horizontal drift and yaw.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum HydrostaticStability {
    /// Every heave/tilt direction has positive stiffness beyond estimated resolution.
    Restoring,
    /// At least one direction has negative stiffness beyond estimated resolution.
    Unstable,
    /// At least one direction is unresolved or neutral at the requested resolution.
    Marginal,
}

/// Complete result. No finite-difference sample or failed spectrum is published alone.
#[derive(Clone, Copy, Debug)]
pub struct HydrostaticStabilityReport {
    pub state: RigidBodyState,
    pub load: HydrostaticLoad,
    pub net_force_world_n: Vec3,
    pub torque_world_nm: Vec3,
    /// Orthonormal (up, roll, pitch) directions in the world frame.
    pub axes_world: [Vec3; 3],
    /// Symmetric fine-step estimate of -d(L*F_up, torque_roll, torque_pitch)/dq,
    /// where q = (heave/L, roll, pitch). Entries have units joules.
    pub restoring_stiffness_j: [[f64; 3]; 3],
    /// Ascending native symmetric eigenvalues; these are NOT squared frequencies.
    pub eigenvalues_j: [f64; 3],
    pub classification: HydrostaticStability,
    /// Absolute numerical margin used for classification, including observed
    /// refinement, symmetrization and floating-point eigenvalue resolution.
    /// An estimate, not an interval bound on derivatives or source uncertainty.
    pub estimated_resolution_j: f64,
    /// Maximum entrywise full-step/half-step difference divided by matrix scale.
    pub relative_refinement_defect: f64,
    /// Maximum entrywise antisymmetric difference divided by matrix scale.
    pub relative_symmetry_defect: f64,
    /// Extrapolated full/half-step one-sided derivative jump divided by matrix
    /// scale. Detects waterline-contact kinks that central differences can hide.
    pub relative_one_sided_defect: f64,
    /// Exactly 13 successful reservoir queries: center and +/- two steps per axis.
    /// Finite tanks also retain their own bounded inner waterline solves.
    pub evaluations: usize,
}

#[derive(Clone, Debug, PartialEq)]
pub enum HydrostaticStabilityError {
    Invalid(&'static str),
    Load(HydrostaticError),
    /// The requested state has not balanced its actual physical loads.
    NotEquilibrium {
        force_defect_n: f64,
        torque_defect_nm: f64,
    },
    /// The chosen perturbations did not resolve a consistent conservative tangent.
    UnresolvedStiffness {
        refinement_defect_j: f64,
        symmetry_defect_j: f64,
        one_sided_defect_j: f64,
        maximum_defect_j: f64,
    },
    Spectrum(SymmetricEigenvalues3Error),
}
impl core::fmt::Display for HydrostaticStabilityError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "hydrostatic stability refused: {self:?}")
    }
}
impl std::error::Error for HydrostaticStabilityError {}
impl From<HydrostaticError> for HydrostaticStabilityError {
    fn from(error: HydrostaticError) -> Self {
        Self::Load(error)
    }
}

/// Differentiate the SAME hydrostatic loads used by native rigid-body motion.
///
/// Requires a static force/torque-balanced state. Rotations are world-frame
/// exponential perturbations about its material COM. Each finite tank re-solves
/// its original conserved inventory at every perturbed pose; material sources
/// stay frozen. No displaced-volume sampling, water default or fixed-tank-level
/// approximation is introduced. Central differences at h and h/2 expose local
/// resolution; only an admitted symmetric tangent reaches the shared, bounded,
/// allocation-free 3x3 eigenvalue kernel.
///
/// Classification concerns local conservative heave and tilt only. Horizontal
/// translation and yaw are neutral in an unconstrained planar reservoir and are
/// excluded. This is not a finite-angle capsize, dynamic damping, added-mass,
/// wave or sloshing analysis. Eigenvalue magnitudes depend on the chosen L;
/// inertia and fluid dynamics are required to infer oscillation frequencies.
pub fn analyze_hydrostatic_stability<R: HydrostaticReservoir + ?Sized>(
    reservoir: &R,
    mesh: &DisplacementMesh,
    body: &PrincipalFrame,
    state: RigidBodyState,
    controls: HydrostaticStabilityControls,
    mut cancelled: impl FnMut() -> bool,
) -> Result<HydrostaticStabilityReport, HydrostaticStabilityError> {
    poll(&mut cancelled)?;
    validate(controls)?;
    if state.linear_momentum_world() != Vec3::ZERO || state.angular_momentum_body() != Vec3::ZERO {
        return Err(HydrostaticStabilityError::Invalid(
            "static zero-momentum state required",
        ));
    }
    let gravity = reservoir.gravity();
    let axes = HydrostaticPoseAxes::new(gravity, controls.roll_direction_world)?;
    let weight = gravity.acceleration_world().scale(body.properties().mass());
    if !weight.is_finite() {
        return Err(HydrostaticStabilityError::Invalid(
            "unrepresentable body weight",
        ));
    }
    let load = reservoir.load_frame(mesh, body, state, &mut cancelled)?;
    let net_force = load.wrench().force_world.add(weight);
    let torque = state
        .pose()
        .orientation()
        .rotate_body_to_world(load.wrench().torque_body);
    let force_defect = finite(magnitude(net_force))?;
    let torque_defect = finite(magnitude(torque))?;
    if force_defect > controls.maximum_force_defect_n
        || torque_defect > controls.maximum_torque_defect_nm
    {
        return Err(HydrostaticStabilityError::NotEquilibrium {
            force_defect_n: force_defect,
            torque_defect_nm: torque_defect,
        });
    }
    let center_effort = [
        finite(net_force.dot(axes.up) * controls.length_scale_m)?,
        finite(torque.dot(axes.roll))?,
        finite(torque.dot(axes.pitch))?,
    ];
    let mut matrices = [[[0.0; 3]; 3]; 2];
    let mut one_sided = [[[0.0; 3]; 3]; 2];
    let mut evaluations = 1;
    for (level, factor) in [1.0, 0.5].into_iter().enumerate() {
        for column in 0..3 {
            let step = if column == 0 {
                controls.heave_step_m
            } else {
                controls.tilt_step_radians
            } * factor;
            let coordinate_step = if column == 0 {
                step / controls.length_scale_m
            } else {
                step
            };
            let mut efforts = [[0.0; 3]; 2];
            for (side, sign) in [-1.0, 1.0].into_iter().enumerate() {
                poll(&mut cancelled)?;
                let mut tilt = [0.0; 2];
                let heave = if column == 0 {
                    sign * step
                } else {
                    tilt[column - 1] = sign * step;
                    0.0
                };
                let sample = axes.displaced_state(state, heave, tilt)?;
                if sample.pose() == state.pose() {
                    return Err(HydrostaticStabilityError::Invalid(
                        "pose perturbation is below floating-point resolution",
                    ));
                }
                if column == 0 {
                    let requested = axes.up.scale(heave);
                    let realized = sample
                        .pose()
                        .position_world()
                        .sub(state.pose().position_world());
                    let relative_pose_defect = magnitude(realized.sub(requested)) / step;
                    if !relative_pose_defect.is_finite()
                        || relative_pose_defect > 0.1 * controls.maximum_relative_stiffness_defect
                    {
                        return Err(HydrostaticStabilityError::Invalid(
                            "heave perturbation is not resolved in every world component",
                        ));
                    }
                }
                let sample_load = reservoir.load_frame(mesh, body, sample, &mut cancelled)?;
                evaluations += 1;
                let f = sample_load.wrench().force_world.add(weight);
                let t = sample
                    .pose()
                    .orientation()
                    .rotate_body_to_world(sample_load.wrench().torque_body);
                efforts[side] = [
                    finite(f.dot(axes.up) * controls.length_scale_m)?,
                    finite(t.dot(axes.roll))?,
                    finite(t.dot(axes.pitch))?,
                ];
            }
            for (row, (minus, plus)) in efforts[0].into_iter().zip(efforts[1]).enumerate() {
                // Divide by step before halving, avoiding overflow in 2*step.
                matrices[level][row][column] = finite((minus - plus) / coordinate_step * 0.5)?;
                one_sided[level][row][column] = finite(
                    ((center_effort[row] - plus) + (center_effort[row] - minus)) / coordinate_step,
                )?;
            }
        }
    }
    let mut scale = 0.0_f64;
    let mut refinement = 0.0_f64;
    let mut symmetry = 0.0_f64;
    let mut derivative_jump = 0.0_f64;
    for row in 0..3 {
        for column in 0..3 {
            scale = scale
                .max(matrices[0][row][column].abs())
                .max(matrices[1][row][column].abs());
            refinement =
                refinement.max(finite(matrices[0][row][column] - matrices[1][row][column])?.abs());
            symmetry =
                symmetry.max(finite(matrices[1][row][column] - matrices[1][column][row])?.abs());
            // Smooth one-sided curvature is O(h) and cancels. A persistent
            // derivative jump survives this extrapolation, even when both
            // central stencils agree (e.g. a just-submerged neutral body).
            derivative_jump = derivative_jump
                .max(finite(2.0 * one_sided[1][row][column] - one_sided[0][row][column])?.abs());
        }
    }
    let maximum_defect = finite(
        controls.neutral_stiffness_tolerance_j + controls.maximum_relative_stiffness_defect * scale,
    )?;
    if refinement > maximum_defect || symmetry > maximum_defect || derivative_jump > maximum_defect
    {
        return Err(HydrostaticStabilityError::UnresolvedStiffness {
            refinement_defect_j: refinement,
            symmetry_defect_j: symmetry,
            one_sided_defect_j: derivative_jump,
            maximum_defect_j: maximum_defect,
        });
    }
    let mut stiffness = matrices[1];
    for row in 0..3 {
        for column in 0..row {
            let value = stiffness[row][column] * 0.5 + stiffness[column][row] * 0.5;
            stiffness[row][column] = value;
            stiffness[column][row] = value;
        }
    }
    poll(&mut cancelled)?;
    let eigenvalues =
        symmetric_eigenvalues3(&stiffness).map_err(HydrostaticStabilityError::Spectrum)?;
    // For a 3x3 matrix, 3*max-entry bounds its spectral norm. The measured
    // differences are estimates of discretization error, not certified bounds.
    let resolution = finite(
        controls.neutral_stiffness_tolerance_j
            + 3.0 * (refinement + derivative_jump + 0.5 * symmetry + 64.0 * f64::EPSILON * scale),
    )?;
    let classification = if eigenvalues[0] > resolution {
        HydrostaticStability::Restoring
    } else if eigenvalues[0] < -resolution {
        HydrostaticStability::Unstable
    } else {
        HydrostaticStability::Marginal
    };
    poll(&mut cancelled)?;
    Ok(HydrostaticStabilityReport {
        state,
        load,
        net_force_world_n: net_force,
        torque_world_nm: torque,
        axes_world: [axes.up, axes.roll, axes.pitch],
        restoring_stiffness_j: stiffness,
        eigenvalues_j: eigenvalues,
        classification,
        estimated_resolution_j: resolution,
        relative_refinement_defect: if scale == 0.0 {
            0.0
        } else {
            refinement / scale
        },
        relative_symmetry_defect: if scale == 0.0 { 0.0 } else { symmetry / scale },
        relative_one_sided_defect: if scale == 0.0 {
            0.0
        } else {
            derivative_jump / scale
        },
        evaluations,
    })
}

fn magnitude(v: Vec3) -> f64 {
    v.x.hypot(v.y).hypot(v.z)
}

fn validate(c: HydrostaticStabilityControls) -> Result<(), HydrostaticStabilityError> {
    if [
        c.length_scale_m,
        c.heave_step_m,
        c.tilt_step_radians,
        c.maximum_force_defect_n,
        c.maximum_torque_defect_nm,
        c.maximum_relative_stiffness_defect,
        c.neutral_stiffness_tolerance_j,
    ]
    .iter()
    .any(|x| !x.is_finite() || *x <= 0.0)
        || c.tilt_step_radians > 0.25
        || c.maximum_relative_stiffness_defect >= 1.0
        || c.heave_step_m * 0.5 == 0.0
        || c.tilt_step_radians * 0.5 == 0.0
        || !(c.heave_step_m / c.length_scale_m).is_finite()
        || !(c.heave_step_m * 0.5 / c.length_scale_m).is_finite()
        || c.heave_step_m * 0.5 / c.length_scale_m == 0.0
    {
        return Err(HydrostaticStabilityError::Invalid(
            "finite positive stiffness scales and tolerances required",
        ));
    }
    Ok(())
}
