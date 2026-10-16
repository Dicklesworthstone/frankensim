//! Source-bounded multibody composition for Otis's US 31,128 hoisting apparatus.
//!
//! The 1861 grant prints the connected topology and operating sequence of
//! platform D; rope G; bar F; levers E; pawls f and hook racks C; drum H;
//! shaft I and pulleys J/K/L; power drum N and belts O/P; shipper S; hand and
//! stop ropes T/U/V; brake shoe Z; and the opposite-wound counterpoise Q/R.
//! It prints no mass, speed, force, spring rate, stopping distance, engagement
//! time, or power datum, so this module owns normalized joint coordinates and
//! claim predicates only.

use crate::articulated::{ArticulatedError, JointKind, JointModel};
use fs_ga::Vec3;

/// Declared display boundary at which the lower-travel stop is treated as reached.
pub const LOWER_LIMIT_NORMALIZED: f64 = 0.03;

/// Inputs to one deterministic US 31,128 topology query.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct OtisTopologyParams {
    /// Normalized platform-D travel in `[0, 1]`, bottom to top.
    pub platform_position_normalized: f64,
    /// Finite display phase for drums, pulleys, and belts.
    pub drive_phase_rad: f64,
    /// Requested direction: `-1` lower, `0` stop, `1` raise.
    pub drive_command: i8,
    /// Whether lifting rope G still connects bar F to drum H.
    pub rope_g_intact: bool,
    /// Whether stop rope U has been pulled.
    pub stop_rope_pulled: bool,
    /// Whether Claim 1's hook-rack and pawl lock remains present.
    pub claim_1_hook_lock_enabled: bool,
    /// Whether Claim 3's shipper-and-brake interlock remains present.
    pub claim_3_brake_interlock_enabled: bool,
    /// Whether Claim 4's opposite-drum counterpoise remains present.
    pub claim_4_counterpoise_enabled: bool,
}

/// Typed refusal from the source-bounded Otis composition.
#[derive(Debug, Clone, PartialEq)]
pub enum OtisTopologyError {
    /// A coordinate or discrete command was outside its declared domain.
    InvalidInput,
    /// The generic multibody joint owner refused the composition.
    Multibody(ArticulatedError),
}

impl From<ArticulatedError> for OtisTopologyError {
    fn from(value: ArticulatedError) -> Self {
        Self::Multibody(value)
    }
}

/// Generic-joint and source-sequence state consumed by browser renderers.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct OtisTopologyStep {
    /// Scalar joint coordinates in the composed mechanism.
    pub scalar_joint_coordinates: usize,
    /// Independent prescribed drives after source linkages are applied.
    pub independent_drive_dofs: usize,
    /// Platform-D prismatic axis.
    pub platform_axis: [f64; 3],
    /// Bar-F prismatic release axis.
    pub safety_bar_axis: [f64; 3],
    /// Lever-E and pawl-f revolute axis.
    pub safety_lever_axis: [f64; 3],
    /// Drum-H revolute axis.
    pub winding_drum_axis: [f64; 3],
    /// Shipper-S prismatic axis.
    pub shipper_axis: [f64; 3],
    /// Brake linkage X/Y revolute axis.
    pub brake_axis: [f64; 3],
    /// Counterpoise-R prismatic axis.
    pub counterpoise_axis: [f64; 3],
    /// Normalized platform-D coordinate.
    pub platform_position_normalized: f64,
    /// Opposed counterpoise coordinate when Claim 4 is present.
    pub counterpoise_position_normalized: f64,
    /// Wrapped display phase.
    pub drive_phase_rad: f64,
    /// Requested drive command.
    pub requested_drive_direction: i8,
    /// Motion direction admitted for platform D.
    pub platform_motion_direction: i8,
    /// Normalized shipper-S coordinate: `-1` O working, `0` idle, `1` P working.
    pub shipper_position_normalized: f64,
    /// Whether straight belt O is on working pulley L.
    pub straight_belt_o_working: bool,
    /// Whether cross-belt P is on working pulley L.
    pub cross_belt_p_working: bool,
    /// Whether both belts occupy idle pulleys J/K.
    pub both_belts_idle: bool,
    /// Whether shoe Z bears on working pulley L.
    pub brake_z_engaged: bool,
    /// Whether rope U and fork V have requested the source stop geometry.
    pub stop_rope_geometry_active: bool,
    /// Whether the lower-limit arm and projection have requested a stop.
    pub lower_limit_stop_active: bool,
    /// Whether rope G remains taut and connected.
    pub rope_g_taut: bool,
    /// Normalized release of bar F after rope-G failure.
    pub safety_bar_release_normalized: f64,
    /// Normalized opposed rotation of the two bent levers E.
    pub safety_lever_rotation_normalized: f64,
    /// Whether pawls f are in the hook racks C.
    pub pawls_f_engaged: bool,
    /// Whether Claim 1's self-locking hook geometry is satisfied.
    pub claim_1_hook_lock_satisfied: bool,
    /// Whether a removed Claim 1 admits a guided free-fall counterfactual.
    pub free_fall_counterfactual: bool,
    /// Whether Claim 3's simultaneous idle-belt/brake state is satisfied.
    pub claim_3_stop_interlock_satisfied: bool,
    /// Whether Q/R remain on the opposite side of H without touching the safety frame.
    pub claim_4_counterpoise_topology_satisfied: bool,
    /// Stable source-order mode label.
    pub mechanism_mode: &'static str,
}

fn revolute_axis(joint: &JointModel) -> [f64; 3] {
    let axis = joint.motion_subspace().angular;
    [axis.x, axis.y, axis.z]
}

fn prismatic_axis(joint: &JointModel) -> [f64; 3] {
    let axis = joint.motion_subspace().linear;
    [axis.x, axis.y, axis.z]
}

/// Compose the source-printed joints and evaluate one operating state.
///
/// # Errors
/// Refuses non-finite coordinates, platform travel outside `[0, 1]`, a drive
/// command other than `-1/0/1`, or a generic multibody-joint refusal.
pub fn step_otis_topology(
    params: OtisTopologyParams,
) -> Result<OtisTopologyStep, OtisTopologyError> {
    if !params.platform_position_normalized.is_finite()
        || !(0.0..=1.0).contains(&params.platform_position_normalized)
        || !params.drive_phase_rad.is_finite()
        || ![-1, 0, 1].contains(&params.drive_command)
    {
        return Err(OtisTopologyError::InvalidInput);
    }

    let platform = JointModel::prismatic(Vec3::new(0.0, 1.0, 0.0), None)?;
    let safety_bar = JointModel::prismatic(Vec3::new(0.0, 1.0, 0.0), None)?;
    let left_lever = JointModel::revolute(Vec3::new(0.0, 0.0, 1.0), None)?;
    let right_lever = JointModel::revolute(Vec3::new(0.0, 0.0, 1.0), None)?;
    let left_pawl = JointModel::revolute(Vec3::new(0.0, 0.0, 1.0), None)?;
    let right_pawl = JointModel::revolute(Vec3::new(0.0, 0.0, 1.0), None)?;
    let winding_drum = JointModel::revolute(Vec3::new(0.0, 0.0, 1.0), None)?;
    let drive_shaft = JointModel::revolute(Vec3::new(0.0, 0.0, 1.0), None)?;
    let power_drum = JointModel::revolute(Vec3::new(0.0, 0.0, 1.0), None)?;
    let shipper = JointModel::prismatic(Vec3::new(1.0, 0.0, 0.0), None)?;
    let brake = JointModel::revolute(Vec3::new(0.0, 0.0, 1.0), None)?;
    let counterpoise = JointModel::prismatic(Vec3::new(0.0, 1.0, 0.0), None)?;
    let joints = [
        platform,
        safety_bar,
        left_lever,
        right_lever,
        left_pawl,
        right_pawl,
        winding_drum,
        drive_shaft,
        power_drum,
        shipper,
        brake,
        counterpoise,
    ];
    debug_assert_eq!(platform.kind(), JointKind::Prismatic);
    debug_assert_eq!(safety_bar.kind(), JointKind::Prismatic);
    debug_assert_eq!(left_lever.kind(), JointKind::Revolute);
    debug_assert_eq!(winding_drum.kind(), JointKind::Revolute);
    debug_assert_eq!(shipper.kind(), JointKind::Prismatic);
    debug_assert_eq!(counterpoise.kind(), JointKind::Prismatic);

    let lower_limit_stop_active =
        params.drive_command < 0 && params.platform_position_normalized <= LOWER_LIMIT_NORMALIZED;
    let stop_requested =
        params.stop_rope_pulled || params.drive_command == 0 || lower_limit_stop_active;
    let straight_belt_o_working = params.drive_command > 0 && !stop_requested;
    let cross_belt_p_working = params.drive_command < 0 && !stop_requested;
    let both_belts_idle = stop_requested;
    let brake_z_engaged = stop_requested && params.claim_3_brake_interlock_enabled;
    let shipper_position_normalized = if straight_belt_o_working {
        -1.0
    } else if cross_belt_p_working {
        1.0
    } else {
        0.0
    };

    let rope_failed = !params.rope_g_intact;
    let pawls_f_engaged = rope_failed && params.claim_1_hook_lock_enabled;
    let free_fall_counterfactual = rope_failed && !params.claim_1_hook_lock_enabled;
    let platform_motion_direction = if rope_failed {
        if free_fall_counterfactual { -1 } else { 0 }
    } else if stop_requested {
        0
    } else {
        params.drive_command
    };
    let mechanism_mode = if pawls_f_engaged {
        "rope-failure-hook-lock"
    } else if free_fall_counterfactual {
        "claim-1-free-fall-counterfactual"
    } else if lower_limit_stop_active {
        "lower-limit-stop"
    } else if stop_requested {
        "service-stop"
    } else if platform_motion_direction > 0 {
        "raise"
    } else {
        "lower"
    };

    Ok(OtisTopologyStep {
        scalar_joint_coordinates: joints.iter().map(|joint| joint.dof_count()).sum(),
        independent_drive_dofs: 1,
        platform_axis: prismatic_axis(&platform),
        safety_bar_axis: prismatic_axis(&safety_bar),
        safety_lever_axis: revolute_axis(&left_lever),
        winding_drum_axis: revolute_axis(&winding_drum),
        shipper_axis: prismatic_axis(&shipper),
        brake_axis: revolute_axis(&brake),
        counterpoise_axis: prismatic_axis(&counterpoise),
        platform_position_normalized: params.platform_position_normalized,
        counterpoise_position_normalized: if params.claim_4_counterpoise_enabled {
            1.0 - params.platform_position_normalized
        } else {
            params.platform_position_normalized
        },
        drive_phase_rad: params.drive_phase_rad.rem_euclid(core::f64::consts::TAU),
        requested_drive_direction: params.drive_command,
        platform_motion_direction,
        shipper_position_normalized,
        straight_belt_o_working,
        cross_belt_p_working,
        both_belts_idle,
        brake_z_engaged,
        stop_rope_geometry_active: stop_requested,
        lower_limit_stop_active,
        rope_g_taut: params.rope_g_intact,
        safety_bar_release_normalized: if rope_failed { 1.0 } else { 0.0 },
        safety_lever_rotation_normalized: if pawls_f_engaged { 1.0 } else { 0.0 },
        pawls_f_engaged,
        claim_1_hook_lock_satisfied: pawls_f_engaged,
        free_fall_counterfactual,
        claim_3_stop_interlock_satisfied: !stop_requested || (both_belts_idle && brake_z_engaged),
        claim_4_counterpoise_topology_satisfied: params.claim_4_counterpoise_enabled,
        mechanism_mode,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn baseline() -> OtisTopologyParams {
        OtisTopologyParams {
            platform_position_normalized: 0.55,
            drive_phase_rad: 0.0,
            drive_command: 1,
            rope_g_intact: true,
            stop_rope_pulled: false,
            claim_1_hook_lock_enabled: true,
            claim_3_brake_interlock_enabled: true,
            claim_4_counterpoise_enabled: true,
        }
    }

    #[test]
    fn composes_the_complete_source_topology_from_generic_joints() {
        let step = step_otis_topology(baseline()).unwrap();
        assert_eq!(step.scalar_joint_coordinates, 12);
        assert_eq!(step.independent_drive_dofs, 1);
        assert_eq!(step.platform_axis, [0.0, 1.0, 0.0]);
        assert_eq!(step.shipper_axis, [1.0, 0.0, 0.0]);
        assert!(step.straight_belt_o_working);
        assert_eq!(step.mechanism_mode, "raise");
    }

    #[test]
    fn stop_idles_both_belts_and_applies_the_claim_3_brake() {
        let step = step_otis_topology(OtisTopologyParams {
            stop_rope_pulled: true,
            ..baseline()
        })
        .unwrap();
        assert!(step.both_belts_idle);
        assert!(step.brake_z_engaged);
        assert!(step.claim_3_stop_interlock_satisfied);
        assert_eq!(step.platform_motion_direction, 0);
    }

    #[test]
    fn rope_failure_locks_claim_1_or_exposes_the_counterfactual() {
        let caught = step_otis_topology(OtisTopologyParams {
            rope_g_intact: false,
            ..baseline()
        })
        .unwrap();
        assert!(caught.pawls_f_engaged);
        assert_eq!(caught.platform_motion_direction, 0);

        let removed = step_otis_topology(OtisTopologyParams {
            rope_g_intact: false,
            claim_1_hook_lock_enabled: false,
            ..baseline()
        })
        .unwrap();
        assert!(removed.free_fall_counterfactual);
        assert_eq!(removed.platform_motion_direction, -1);
    }

    #[test]
    fn lower_stop_and_counterpoise_follow_the_printed_arrangement() {
        let step = step_otis_topology(OtisTopologyParams {
            platform_position_normalized: 0.02,
            drive_command: -1,
            ..baseline()
        })
        .unwrap();
        assert!(step.lower_limit_stop_active);
        assert!(step.brake_z_engaged);
        assert_eq!(step.counterpoise_position_normalized, 0.98);
    }

    #[test]
    fn invalid_coordinates_and_commands_refuse() {
        assert_eq!(
            step_otis_topology(OtisTopologyParams {
                platform_position_normalized: 1.01,
                ..baseline()
            }),
            Err(OtisTopologyError::InvalidInput)
        );
        assert_eq!(
            step_otis_topology(OtisTopologyParams {
                drive_command: 2,
                ..baseline()
            }),
            Err(OtisTopologyError::InvalidInput)
        );
    }
}
