//! Source-bounded multibody composition for Nikolaus Otto's US 194,047 engine.
//!
//! The grant establishes one continuously driven crank, a piston constrained
//! to the cylinder axis, a connecting rod seated on the wrist and crank pins,
//! a half-speed side shaft, a sliding admission valve, and an exhaust-valve
//! linkage. It does not print build dimensions, masses, inertia, rpm, or valve
//! lift. This module therefore owns the generic joint topology, exact
//! slider-crank closure, the 2:1 shaft relation, and normalized valve timing.
//! Display geometry is supplied by the caller and is never presented as a
//! historical dimension.

use crate::articulated::{ArticulatedError, JointKind, JointModel};
use fs_ga::Vec3;

/// Inputs to one deterministic, display-scale mechanism query.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct OttoTopologyParams {
    /// Continuous four-stroke cycle coordinate in radians. Zero is intake TDC.
    pub crank_angle_rad: f64,
    /// Display crank radius in caller-declared model units.
    pub crank_radius: f64,
    /// Display wrist-to-crank pin distance in the same model units.
    pub connecting_rod_length: f64,
    /// Engine speed used only for the normalized governor presentation pose.
    pub engine_rpm: f64,
}

/// Typed refusal from the source-bounded Otto composition.
#[derive(Debug, Clone, PartialEq)]
pub enum OttoTopologyError {
    /// An input was non-finite, non-positive, or outside the declared domain.
    InvalidInput,
    /// The rod is too short to close over the selected crank radius.
    ImpossibleLinkage,
    /// The generic multibody joint owner refused the composition.
    Multibody(ArticulatedError),
}

impl From<ArticulatedError> for OttoTopologyError {
    fn from(value: ArticulatedError) -> Self {
        Self::Multibody(value)
    }
}

/// Generic-joint topology and one closed mechanism pose.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct OttoTopologyStep {
    /// Scalar joint coordinates in the composed topology.
    pub scalar_joint_coordinates: usize,
    /// Independent prescribed drives after the crank linkages are applied.
    pub independent_drive_dofs: usize,
    /// Main-crank revolute axis.
    pub crank_axis: [f64; 3],
    /// Piston prismatic axis.
    pub piston_axis: [f64; 3],
    /// Half-speed side-shaft revolute axis.
    pub side_shaft_axis: [f64; 3],
    /// Admission-slide prismatic axis.
    pub slide_valve_axis: [f64; 3],
    /// Exhaust-valve prismatic axis.
    pub exhaust_valve_axis: [f64; 3],
    /// Governor-spindle revolute axis.
    pub governor_axis: [f64; 3],
    /// Wrapped four-stroke cycle angle in `[0, 4π)`.
    pub cycle_angle_rad: f64,
    /// Crank pin x coordinate relative to the shaft center.
    pub crank_pin_x: f64,
    /// Crank pin y coordinate relative to the shaft center.
    pub crank_pin_y: f64,
    /// Piston wrist-pin x coordinate relative to the shaft center.
    pub piston_pin_x: f64,
    /// Piston wrist-pin y coordinate relative to the shaft center.
    pub piston_pin_y: f64,
    /// Connecting-rod angle from wrist pin toward crank pin.
    pub connecting_rod_angle_rad: f64,
    /// Recomputed wrist-to-crank distance for closure verification.
    pub connecting_rod_span: f64,
    /// Side-shaft angle at the source-required 2:1 ratio.
    pub side_shaft_angle_rad: f64,
    /// Normalized admission-slide displacement in `[-1, 1]`.
    pub slide_valve_normalized: f64,
    /// Normalized exhaust-valve lift in `[0, 1]`.
    pub exhaust_lift_normalized: f64,
    /// Normalized governor spread in `[0, 1]`; a declared display pose.
    pub governor_spread_normalized: f64,
    /// Stable four-stroke phase label.
    pub cycle_phase: &'static str,
}

fn revolute_axis(joint: &JointModel) -> [f64; 3] {
    let axis = joint.motion_subspace().angular;
    [axis.x, axis.y, axis.z]
}

fn prismatic_axis(joint: &JointModel) -> [f64; 3] {
    let axis = joint.motion_subspace().linear;
    [axis.x, axis.y, axis.z]
}

/// Compose the source-printed joints and close one slider-crank pose.
///
/// # Errors
/// Refuses non-finite inputs, non-positive geometry, rpm outside `[0, 600]`,
/// a rod no longer than the crank radius, or a generic joint refusal.
pub fn step_otto_topology(
    params: OttoTopologyParams,
) -> Result<OttoTopologyStep, OttoTopologyError> {
    if !params.crank_angle_rad.is_finite()
        || !params.crank_radius.is_finite()
        || !params.connecting_rod_length.is_finite()
        || !params.engine_rpm.is_finite()
        || params.crank_radius <= 0.0
        || params.connecting_rod_length <= 0.0
        || !(0.0..=600.0).contains(&params.engine_rpm)
    {
        return Err(OttoTopologyError::InvalidInput);
    }
    if params.connecting_rod_length <= params.crank_radius {
        return Err(OttoTopologyError::ImpossibleLinkage);
    }

    let crank = JointModel::revolute(Vec3::new(0.0, 0.0, 1.0), None)?;
    let piston = JointModel::prismatic(Vec3::new(1.0, 0.0, 0.0), None)?;
    let side_shaft = JointModel::revolute(Vec3::new(1.0, 0.0, 0.0), None)?;
    let slide_valve = JointModel::prismatic(Vec3::new(1.0, 0.0, 0.0), None)?;
    let exhaust_valve = JointModel::prismatic(Vec3::new(0.0, 1.0, 0.0), None)?;
    let exhaust_rocker = JointModel::revolute(Vec3::new(0.0, 0.0, 1.0), None)?;
    let governor = JointModel::revolute(Vec3::new(0.0, 1.0, 0.0), None)?;
    let governor_sleeve = JointModel::prismatic(Vec3::new(0.0, 1.0, 0.0), None)?;
    let joints = [
        crank,
        piston,
        side_shaft,
        slide_valve,
        exhaust_valve,
        exhaust_rocker,
        governor,
        governor_sleeve,
    ];
    debug_assert_eq!(crank.kind(), JointKind::Revolute);
    debug_assert_eq!(piston.kind(), JointKind::Prismatic);
    debug_assert_eq!(slide_valve.kind(), JointKind::Prismatic);
    debug_assert_eq!(exhaust_valve.kind(), JointKind::Prismatic);

    let cycle_angle = params
        .crank_angle_rad
        .rem_euclid(2.0 * core::f64::consts::TAU);
    // The public four-stroke coordinate starts the intake stroke at TDC. In
    // this left-facing cylinder projection, TDC puts the crank pin on -x.
    let crank_angle = (cycle_angle + core::f64::consts::PI).rem_euclid(core::f64::consts::TAU);
    let crank_pin_x = params.crank_radius * crank_angle.cos();
    let crank_pin_y = params.crank_radius * crank_angle.sin();
    let radial_y = params.crank_radius * crank_angle.sin();
    let closure = params
        .connecting_rod_length
        .mul_add(params.connecting_rod_length, -(radial_y * radial_y));
    if !closure.is_finite() || closure < 0.0 {
        return Err(OttoTopologyError::ImpossibleLinkage);
    }
    // The cylinder lies to the left of the crank in the museum projection.
    let piston_pin_x = crank_pin_x - closure.sqrt();
    let piston_pin_y = 0.0;
    let rod_dx = crank_pin_x - piston_pin_x;
    let rod_dy = crank_pin_y - piston_pin_y;
    let rod_span = rod_dx.hypot(rod_dy);
    let rod_angle = rod_dy.atan2(rod_dx);

    let side_shaft_angle = 0.5 * cycle_angle;
    let slide_valve_normalized = side_shaft_angle.sin();
    let exhaust_start = 3.0 * core::f64::consts::PI;
    let exhaust_lift_normalized = if cycle_angle >= exhaust_start {
        (cycle_angle - exhaust_start).sin().max(0.0)
    } else {
        0.0
    };
    let cycle_phase = if cycle_angle < core::f64::consts::PI {
        "intake"
    } else if cycle_angle < 2.0 * core::f64::consts::PI {
        "compression"
    } else if cycle_angle < 3.0 * core::f64::consts::PI {
        "power"
    } else {
        "exhaust"
    };

    Ok(OttoTopologyStep {
        scalar_joint_coordinates: joints.iter().map(|joint| joint.dof_count()).sum(),
        independent_drive_dofs: 1,
        crank_axis: revolute_axis(&crank),
        piston_axis: prismatic_axis(&piston),
        side_shaft_axis: revolute_axis(&side_shaft),
        slide_valve_axis: prismatic_axis(&slide_valve),
        exhaust_valve_axis: prismatic_axis(&exhaust_valve),
        governor_axis: revolute_axis(&governor),
        cycle_angle_rad: cycle_angle,
        crank_pin_x,
        crank_pin_y,
        piston_pin_x,
        piston_pin_y,
        connecting_rod_angle_rad: rod_angle,
        connecting_rod_span: rod_span,
        side_shaft_angle_rad: side_shaft_angle,
        slide_valve_normalized,
        exhaust_lift_normalized,
        governor_spread_normalized: (params.engine_rpm / 300.0).clamp(0.0, 1.0),
        cycle_phase,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn params(angle: f64) -> OttoTopologyParams {
        OttoTopologyParams {
            crank_angle_rad: angle,
            crank_radius: 0.65,
            connecting_rod_length: 2.4,
            engine_rpm: 180.0,
        }
    }

    #[test]
    fn composes_eight_coordinates_from_one_crank_drive() {
        let step = step_otto_topology(params(0.0)).expect("valid topology");
        assert_eq!(step.scalar_joint_coordinates, 8);
        assert_eq!(step.independent_drive_dofs, 1);
        assert_eq!(step.crank_axis, [0.0, 0.0, 1.0]);
        assert_eq!(step.piston_axis, [1.0, 0.0, 0.0]);
        assert_eq!(step.side_shaft_axis, [1.0, 0.0, 0.0]);
        assert_eq!(step.exhaust_valve_axis, [0.0, 1.0, 0.0]);
    }

    #[test]
    fn connecting_rod_closes_at_quadrature_and_dead_centers() {
        for angle in [
            0.0,
            0.5 * core::f64::consts::PI,
            core::f64::consts::PI,
            1.5 * core::f64::consts::PI,
        ] {
            let step = step_otto_topology(params(angle)).expect("valid mechanism pose");
            assert!((step.connecting_rod_span - 2.4).abs() < 1.0e-12);
            assert_eq!(step.piston_pin_y, 0.0);
        }
    }

    #[test]
    fn side_shaft_and_exhaust_follow_the_four_stroke_cycle() {
        let power = step_otto_topology(params(2.5 * core::f64::consts::PI)).unwrap();
        let exhaust = step_otto_topology(params(3.5 * core::f64::consts::PI)).unwrap();
        assert_eq!(power.cycle_phase, "power");
        assert_eq!(power.exhaust_lift_normalized, 0.0);
        assert_eq!(exhaust.cycle_phase, "exhaust");
        assert!((exhaust.exhaust_lift_normalized - 1.0).abs() < 1.0e-12);
        assert!((exhaust.side_shaft_angle_rad - 1.75 * core::f64::consts::PI).abs() < 1.0e-12);
    }

    #[test]
    fn four_stroke_labels_match_piston_direction() {
        let intake_tdc = step_otto_topology(params(0.0)).unwrap();
        let intake_bdc = step_otto_topology(params(core::f64::consts::PI)).unwrap();
        let compression_tdc = step_otto_topology(params(2.0 * core::f64::consts::PI)).unwrap();
        let power_bdc = step_otto_topology(params(3.0 * core::f64::consts::PI)).unwrap();

        assert_eq!(intake_tdc.cycle_phase, "intake");
        assert_eq!(intake_bdc.cycle_phase, "compression");
        assert_eq!(compression_tdc.cycle_phase, "power");
        assert_eq!(power_bdc.cycle_phase, "exhaust");
        assert!(intake_tdc.piston_pin_x < intake_bdc.piston_pin_x);
        assert!((compression_tdc.piston_pin_x - intake_tdc.piston_pin_x).abs() < 1.0e-12);
        assert!((power_bdc.piston_pin_x - intake_bdc.piston_pin_x).abs() < 1.0e-12);
    }

    #[test]
    fn invalid_and_impossible_geometry_refuses() {
        assert_eq!(
            step_otto_topology(OttoTopologyParams {
                crank_angle_rad: f64::NAN,
                ..params(0.0)
            }),
            Err(OttoTopologyError::InvalidInput)
        );
        assert_eq!(
            step_otto_topology(OttoTopologyParams {
                connecting_rod_length: 0.65,
                ..params(0.0)
            }),
            Err(OttoTopologyError::ImpossibleLinkage)
        );
        assert_eq!(
            step_otto_topology(OttoTopologyParams {
                engine_rpm: 600.1,
                ..params(0.0)
            }),
            Err(OttoTopologyError::InvalidInput)
        );
    }
}
