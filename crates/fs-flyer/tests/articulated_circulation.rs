//! Native nonlinear lift from actual link motion, including coupled joint work.
use std::cell::RefCell;

use fs_flyer::Refusal;
use fs_flyer::addedmass::Strip;
use fs_flyer::addedmass::articulated::{LinkedStrip, assemble_articulated_strips};
use fs_flyer::partitioned::articulated::aerodynamics::LinkedLiftingSurfaces;
use fs_flyer::partitioned::articulated::{
    ArticulatedLoads, ArticulatedState, ArticulatedStepConfig, JointImpedance, step,
};
use fs_ga::{Mat3, Se3, So3, So3Tangent, Twist, Vec3, Wrench};
use fs_mbd::articulated::generalized::{GeneralizedKinematics, free_floating_terms};
use fs_mbd::articulated::motion::FreeFloatingMotionState;
use fs_mbd::articulated::{
    ArticulatedModel, FreeFloatingBaseState, JointModel, Link, SpatialInertia,
};
use fs_wing::nonlinear::{
    InfluenceOperator, SectionFrame, StripRegime, StripSpec, solve_nonlinear_with_section_frames,
};
use fs_wing::wrench::CirculatoryWrenchAccumulator;
use fs_wing::{Panel, SurfaceId, flat_surface};

const RHO: f64 = 1.2;
const REFERENCE: [f64; 3] = [10.0, 0.0, 0.5];
fn v(x: [f64; 3]) -> Vec3 {
    Vec3::new(x[0], x[1], x[2])
}
fn a(x: Vec3) -> [f64; 3] {
    [x.x, x.y, x.z]
}
fn close(x: Vec3, y: Vec3, tolerance: f64) {
    assert!((x - y).norm() < tolerance, "{x:?} != {y:?}");
}
fn closure(_: usize, alpha: f64) -> (f64, StripRegime) {
    (core::f64::consts::TAU * alpha, StripRegime::Attached)
}
fn inertia(mass: f64, center: Vec3) -> SpatialInertia {
    SpatialInertia::new(
        mass,
        center,
        Mat3 {
            m: [1.0, 0.0, 0.0, 0.0, 1.1, 0.0, 0.0, 0.0, 1.2],
        },
    )
    .unwrap()
}
fn model(hinged: bool) -> ArticulatedModel {
    let mut links = vec![Link::new(
        "carrier",
        None,
        Se3::identity(),
        JointModel::FIXED,
        inertia(10.0, Vec3::default()),
    )];
    if hinged {
        links.push(Link::new(
            "pitching surface",
            Some(0),
            Se3::identity(),
            JointModel::revolute(Vec3::new(0.0, 1.0, 0.0), None).unwrap(),
            inertia(2.0, Vec3::new(0.5, 0.0, 0.0)),
        ));
    }
    ArticulatedModel::new(links).unwrap()
}
fn geometry() -> (Vec<Panel>, Vec<StripSpec>) {
    let panels = flat_surface(SurfaceId::WingLower, 2.0, 1.0, 1.0, 0.0, 2, 1).unwrap();
    let strips = (0..2)
        .map(|s| StripSpec {
            panel_indices: vec![s],
            chord_m: 1.0,
            twist_rad: 0.0,
        })
        .collect();
    (panels, strips)
}
fn state(hinged: bool, angle: f64, rate: f64) -> ArticulatedState {
    ArticulatedState {
        motion: FreeFloatingMotionState {
            base: FreeFloatingBaseState::new(
                Se3::identity(),
                Twist::new(Vec3::default(), v(REFERENCE)),
            ),
            joint_position: if hinged { vec![angle] } else { vec![] },
            joint_velocity: if hinged { vec![rate] } else { vec![] },
        },
        memory: vec![],
    }
}
fn map(model: &ArticulatedModel, state: &ArticulatedState) -> GeneralizedKinematics {
    free_floating_terms(
        model,
        state.motion.base,
        &state.motion.joint_position,
        &state.motion.joint_velocity,
        || false,
    )
    .unwrap()
    .kinematics
}
fn wind(point: [f64; 3]) -> Result<[f64; 3], Refusal> {
    // Smooth prescribed potential strain with a uniform translational part.
    Ok([0.1 + 0.03 * point[1], 0.03 * point[0], -0.05])
}

#[test]
fn frozen_root_recovers_original_nonlinear_circulation_and_wrench() {
    let model = model(false);
    let state = state(false, 0.0, 0.0);
    let kinematics = map(&model, &state);
    let (panels, strips) = geometry();
    let surfaces = LinkedLiftingSurfaces::new(&panels, &strips, &[0, 0]).unwrap();
    let flow = [0.3, -0.1, 0.05];
    let relative = a(v(REFERENCE) - v(flow));
    let actual = surfaces
        .evaluate(
            &kinematics,
            relative,
            RHO,
            &closure,
            |_| Ok(flow),
            None,
            || false,
        )
        .unwrap();
    let frames = panels
        .iter()
        .map(SectionFrame::from_panel)
        .collect::<Result<Vec<_>, _>>()
        .unwrap();
    let op = InfluenceOperator::build(&panels, relative, RHO).unwrap();
    let expected = solve_nonlinear_with_section_frames(
        &op,
        &panels,
        &strips,
        relative,
        RHO,
        &closure,
        &frames,
        &[relative; 2],
        None,
    )
    .unwrap();
    assert_eq!(actual.nonlinear, expected);
    assert_eq!(actual.strip_relative_velocity_world_mps, [relative; 2]);
    let mut wrench = CirculatoryWrenchAccumulator::new(RHO, [0.0; 3]).unwrap();
    for (panel, gamma) in panels.iter().zip(expected.gamma) {
        wrench.add_panel(panel, gamma, relative).unwrap();
    }
    assert_eq!(actual.world_wrenches, [wrench.total()]);
}

#[test]
fn joint_point_velocity_changes_real_loads_and_projects_the_same_virtual_work() {
    let model = model(true);
    let (panels, strips) = geometry();
    let surfaces = LinkedLiftingSurfaces::new(&panels, &strips, &[1, 1]).unwrap();
    let mut state = state(true, 0.12, 0.3);
    state.motion.base.twist_body.angular = Vec3::new(0.02, -0.03, 0.04);
    let kinematics = map(&model, &state);
    let actual = surfaces
        .evaluate(&kinematics, REFERENCE, RHO, &closure, wind, None, || false)
        .unwrap();
    let (sin, cos) = (0.12_f64.sin(), 0.12_f64.cos());
    let rotate = |p: Vec3| Vec3::new(cos * p.x + sin * p.z, p.y, -sin * p.x + cos * p.z);
    let angular = state.motion.base.twist_body.angular + Vec3::new(0.0, 0.3, 0.0);
    let mut point_power = 0.0;
    for (s, panel) in panels.iter().enumerate() {
        let point = rotate((v(panel.a) + v(panel.b)).scale(0.5));
        let material = v(REFERENCE) + angular.cross(point);
        let relative = material - v(wind(a(point)).unwrap());
        close(v(actual.strip_points_world_m[s]), point, 2e-14);
        close(
            v(actual.strip_relative_velocity_world_mps[s]),
            relative,
            2e-14,
        );
        let force = relative
            .cross(rotate(v(panel.b) - v(panel.a)))
            .scale(RHO * actual.nonlinear.gamma[s]);
        point_power += force.dot(material);
    }
    let generalized = kinematics
        .project_world_wrenches(&actual.world_wrenches)
        .unwrap();
    let mut velocity = state.motion.base.twist_body.to_array().to_vec();
    velocity.extend(&state.motion.joint_velocity);
    let projected_power: f64 = generalized
        .iter()
        .zip(velocity)
        .map(|(force, rate)| force * rate)
        .sum();
    assert!((projected_power - point_power).abs() < 2e-11);
    assert!(generalized[6].abs() > 0.1);
    let mut stopped = state.clone();
    stopped.motion.joint_velocity[0] = 0.0;
    let frozen = surfaces
        .evaluate(
            &map(&model, &stopped),
            REFERENCE,
            RHO,
            &closure,
            wind,
            None,
            || false,
        )
        .unwrap();
    assert!((actual.nonlinear.gamma[0] - frozen.nonlinear.gamma[0]).abs() > 1e-3);
    assert!((actual.world_wrenches[1].torque.y - frozen.world_wrenches[1].torque.y).abs() > 0.1);
}

#[test]
fn common_rotation_translation_and_wind_boost_preserve_physical_loading() {
    let model = model(true);
    let (panels, strips) = geometry();
    let surfaces = LinkedLiftingSurfaces::new(&panels, &strips, &[1, 1]).unwrap();
    let original = state(true, 0.08, 0.25);
    let expected = surfaces
        .evaluate(
            &map(&model, &original),
            REFERENCE,
            RHO,
            &closure,
            wind,
            None,
            || false,
        )
        .unwrap();
    let rotation = So3::exp(So3Tangent::new(Vec3::new(0.2, -0.3, 0.4))).unwrap();
    let position = Vec3::new(2.0, -3.0, 1.0);
    let pose = Se3::from_parts(rotation, position).unwrap();
    let boost = Vec3::new(100.0, -50.0, 30.0);
    let mut shifted = original;
    shifted.motion.base.world_from_base = pose;
    shifted.motion.base.twist_body.linear =
        v(REFERENCE) + rotation.inverse().rotate(boost).unwrap();
    let transformed_wind = |point| {
        let local = pose.inverse().unwrap().transform_point(v(point)).unwrap();
        Ok(a(rotation.rotate(v(wind(a(local))?)).unwrap() + boost))
    };
    let actual = surfaces
        .evaluate(
            &map(&model, &shifted),
            a(rotation.rotate(v(REFERENCE)).unwrap()),
            RHO,
            &closure,
            transformed_wind,
            None,
            || false,
        )
        .unwrap();
    for (a, b) in actual.nonlinear.gamma.iter().zip(&expected.nonlinear.gamma) {
        assert!((a - b).abs() < 2e-9);
    }
    for (actual, expected) in actual.world_wrenches.iter().zip(&expected.world_wrenches) {
        close(actual.force, rotation.rotate(expected.force).unwrap(), 2e-8);
        // Each torque remains about its own translated link origin.
        close(
            actual.torque,
            rotation.rotate(expected.torque).unwrap(),
            2e-8,
        );
    }
}

#[test]
fn actual_articulated_midpoint_couples_circulation_back_into_joint_velocity() {
    let model = model(true);
    let (panels, strips) = geometry();
    let surfaces = LinkedLiftingSurfaces::new(&panels, &strips, &[1, 1]).unwrap();
    let initial = state(true, 0.02, 0.0);
    let before = initial.clone();
    let fluid = [LinkedStrip::new(
        1,
        Strip {
            name: "wing",
            chord_m: 1.0,
            span_m: 2.0,
            position_m: [0.75, 0.0, 0.0],
            normal: [0.0, 0.0, 1.0],
            span_axis: [0.0, 1.0, 0.0],
            control_coord: None,
            control_gain: 0.0,
            control_angular_gain: [0.0; 3],
        },
    )
    .unwrap()];
    let impedance = [JointImpedance {
        damping: 0.0,
        stiffness: 0.0,
        rest_position: 0.0,
    }];
    let config = ArticulatedStepConfig {
        max_iterations: 64,
        linear_tolerance_mps: 1e-10,
        angular_tolerance_rad_s: 1e-10,
        joint_velocity_tolerances: &[1e-10],
    };
    let seen = RefCell::new(Vec::new());
    let advance = |circulation: bool| {
        step(
            &model,
            &impedance,
            &initial,
            &[],
            0.0,
            0.002,
            config,
            |_, stage, kinematics| {
                let loads = if circulation {
                    let result = surfaces.evaluate(
                        kinematics,
                        REFERENCE,
                        RHO,
                        &closure,
                        |_| Ok([0.0; 3]),
                        None,
                        || false,
                    )?;
                    seen.borrow_mut()
                        .push((stage.motion.joint_velocity[0], result.nonlinear.gamma[0]));
                    result.world_wrenches
                } else {
                    vec![Wrench::default(); kinematics.link_count()]
                };
                Ok(ArticulatedLoads {
                    fluid: assemble_articulated_strips(RHO, &fluid, kinematics)?,
                    world_wrenches: loads,
                    gravity_world: Vec3::default(),
                    joint_effort: vec![0.0],
                })
            },
            || false,
        )
        .unwrap()
    };
    let coupled = advance(true);
    let inertial_only = advance(false);
    assert!(
        (coupled.state.motion.joint_velocity[0] - inertial_only.state.motion.joint_velocity[0])
            .abs()
            > 1e-5
    );
    assert!(coupled.state.motion.joint_velocity[0].abs() > 1e-5);
    let samples = seen.borrow();
    assert!(samples.len() > 1);
    assert!(
        samples
            .windows(2)
            .any(|p| (p[1].0 - p[0].0).abs() > 1e-8 && (p[1].1 - p[0].1).abs() > 1e-8)
    );
    assert_eq!(initial, before);
}

#[test]
fn missing_links_invalid_flow_and_late_cancellation_return_no_partial_loads() {
    let model = model(true);
    let initial = state(true, 0.02, 0.1);
    let before = initial.clone();
    let kinematics = map(&model, &initial);
    let (panels, strips) = geometry();
    assert!(LinkedLiftingSurfaces::new(&panels, &[], &[]).is_err());
    let too_many = vec![strips[0].clone(); panels.len() + 1];
    assert!(LinkedLiftingSurfaces::new(&panels, &too_many, &[1, 1, 1]).is_err());
    let missing = LinkedLiftingSurfaces::new(&panels, &strips, &[2, 2]).unwrap();
    assert_eq!(
        missing
            .evaluate(&kinematics, REFERENCE, RHO, &closure, wind, None, || false)
            .unwrap_err()
            .code,
        "articulated-circulation-invalid"
    );
    let surfaces = LinkedLiftingSurfaces::new(&panels, &strips, &[1, 1]).unwrap();
    assert_eq!(
        surfaces
            .evaluate(
                &kinematics,
                REFERENCE,
                RHO,
                &closure,
                |_| Ok([30.0, 0.0, 0.0]),
                None,
                || false
            )
            .unwrap_err()
            .code,
        "incident-flow-invalid"
    );
    let mut polls = 0;
    surfaces
        .evaluate(&kinematics, REFERENCE, RHO, &closure, wind, None, || {
            polls += 1;
            false
        })
        .unwrap();
    let mut seen = 0;
    assert_eq!(
        surfaces
            .evaluate(&kinematics, REFERENCE, RHO, &closure, wind, None, || {
                seen += 1;
                seen == polls
            })
            .unwrap_err()
            .code,
        "articulated-circulation-cancelled"
    );
    assert_eq!(seen, polls);
    assert_eq!(initial, before);
}

#[test]
fn chord_rows_resolve_control_and_bound_motion_and_preserve_joint_virtual_work() {
    let model = model(true);
    let panels = flat_surface(SurfaceId::WingLower, 2.0, 1.0, 1.0, 0.0, 2, 3).unwrap();
    let strips: Vec<_> = (0..2)
        .map(|s| StripSpec {
            panel_indices: vec![s, 2 + s, 4 + s],
            chord_m: 1.0,
            twist_rad: 0.0,
        })
        .collect();
    let surfaces = LinkedLiftingSurfaces::new(&panels, &strips, &[1, 1]).unwrap();
    let mut original = state(true, 0.12, 0.3);
    original.motion.base.twist_body.angular = Vec3::new(0.02, -0.03, 0.04);
    let kinematics = map(&model, &original);
    let samples = RefCell::new(Vec::new());
    let actual = surfaces
        .evaluate(
            &kinematics,
            REFERENCE,
            RHO,
            &closure,
            |point| {
                samples.borrow_mut().push(point);
                wind(point)
            },
            None,
            || false,
        )
        .unwrap();
    assert_eq!(samples.borrow().len(), 2 * panels.len());
    let (sin, cos) = (0.12_f64.sin(), 0.12_f64.cos());
    let rotate = |p: Vec3| Vec3::new(cos * p.x + sin * p.z, p.y, -sin * p.x + cos * p.z);
    let angular = original.motion.base.twist_body.angular + Vec3::new(0.0, 0.3, 0.0);
    let mut force = Vec3::default();
    let mut torque = Vec3::default();
    let mut power = 0.0;
    for (j, panel) in panels.iter().enumerate() {
        let control = rotate(v(panel.ctrl));
        let point = rotate((v(panel.a) + v(panel.b)).scale(0.5));
        let material = v(REFERENCE) + angular.cross(point);
        let relative = material - v(wind(a(point)).unwrap());
        let control_relative = v(REFERENCE) + angular.cross(control) - v(wind(a(control)).unwrap());
        close(v(samples.borrow()[2 * j]), control, 2e-14);
        close(v(samples.borrow()[2 * j + 1]), point, 2e-14);
        close(
            v(actual.panel_control_relative_velocity_world_mps[j]),
            control_relative,
            2e-14,
        );
        close(
            v(actual.panel_bound_relative_velocity_world_mps[j]),
            relative,
            2e-14,
        );
        assert!((control_relative - relative).norm() > 0.01);
        let row_force = relative
            .cross(rotate(v(panel.b) - v(panel.a)))
            .scale(RHO * actual.nonlinear.gamma[j]);
        force = force + row_force;
        torque = torque + point.cross(row_force);
        power += row_force.dot(material);
    }
    close(actual.world_wrenches[1].force, force, 2e-11);
    close(actual.world_wrenches[1].torque, torque, 2e-11);
    close(v(actual.nonlinear.total_force_n), force, 2e-11);
    let generalized = kinematics
        .project_world_wrenches(&actual.world_wrenches)
        .unwrap();
    let mut rates = original.motion.base.twist_body.to_array().to_vec();
    rates.extend(&original.motion.joint_velocity);
    let projected: f64 = generalized.iter().zip(rates).map(|(f, u)| f * u).sum();
    assert!((projected - power).abs() < 3e-11);
    assert!(
        (v(actual.panel_bound_relative_velocity_world_mps[0])
            - v(actual.panel_bound_relative_velocity_world_mps[4]))
        .norm()
            > 0.1
    );
    // Actual native row loading must respond to the differing control-point
    // velocities. The old constant-reference row split cannot reproduce it.
    let world_panels: Vec<_> = panels
        .iter()
        .map(|p| Panel {
            a: a(rotate(v(p.a))),
            b: a(rotate(v(p.b))),
            ctrl: a(rotate(v(p.ctrl))),
            normal: a(rotate(v(p.normal))),
            ..*p
        })
        .collect();
    let frames: Vec<_> = strips
        .iter()
        .map(|s| SectionFrame::from_panel(&world_panels[s.panel_indices[0]]).unwrap())
        .collect();
    let frozen_operator = InfluenceOperator::build(&world_panels, REFERENCE, RHO).unwrap();
    let frozen = solve_nonlinear_with_section_frames(
        &frozen_operator,
        &world_panels,
        &strips,
        REFERENCE,
        RHO,
        &closure,
        &frames,
        &actual.strip_relative_velocity_world_mps,
        None,
    )
    .unwrap();
    assert!(
        actual
            .nonlinear
            .gamma
            .iter()
            .zip(&frozen.gamma)
            .any(|(a, b)| (a - b).abs() > 1e-4)
    );

    // Frame covariance includes both sampled locations and the resolved row
    // loads, with a common wind/body boost removed before rotational slip.
    let rotation = So3::exp(So3Tangent::new(Vec3::new(0.2, -0.3, 0.4))).unwrap();
    let pose = Se3::from_parts(rotation, Vec3::new(2.0, -3.0, 1.0)).unwrap();
    let boost = Vec3::new(100.0, -50.0, 30.0);
    let mut moved = original;
    moved.motion.base.world_from_base = pose;
    moved.motion.base.twist_body.linear = v(REFERENCE) + rotation.inverse().rotate(boost).unwrap();
    let transformed = surfaces
        .evaluate(
            &map(&model, &moved),
            a(rotation.rotate(v(REFERENCE)).unwrap()),
            RHO,
            &closure,
            |point| {
                let local = pose.inverse().unwrap().transform_point(v(point)).unwrap();
                Ok(a(rotation.rotate(v(wind(a(local))?)).unwrap() + boost))
            },
            None,
            || false,
        )
        .unwrap();
    for (a, b) in transformed
        .nonlinear
        .gamma
        .iter()
        .zip(&actual.nonlinear.gamma)
    {
        assert!((a - b).abs() < 3e-9);
    }
    close(
        transformed.world_wrenches[1].force,
        rotation.rotate(force).unwrap(),
        3e-8,
    );
    close(
        transformed.world_wrenches[1].torque,
        rotation.rotate(torque).unwrap(),
        3e-8,
    );
}
