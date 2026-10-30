//! G1/G3/G4: original rigid bodies, finite bonded geometry, fracture and retry.
use fs_couple::cohesive_contact::objective::{
    ObjectiveCohesiveConfig, ObjectiveCohesivePoint, source::SourceBoundObjectiveCohesivePoint,
};
use fs_couple::cohesive_contact::rigid::*;
use fs_evidence::ValidityDomain;
use fs_exec::CancelGate;
use fs_matdb::{
    ClaimSet, InterfaceSystemCard, InterpolationPolicy, MaterialStateId, PropertyClaim,
    PropertyKey, PropertyValue, Provenance, QueryPoint, SurfaceSpec, SystemContext,
    UncertaintyModel,
};
use fs_material::cohesive::{CohesivePoint, source::*};
use fs_material::state_point::{MaterialPropertySelection, interface_stiffness_per_area_dims};
use fs_mbd::{
    Gravity, MassProperties, Pose, RigidBodyIntegrator, RigidBodyState, UnitQuaternion, Vec3,
    Wrench,
};
use fs_qty::Pressure;

fn source(viscous: bool) -> ResolvedBilinearCohesive {
    let mut claims = ClaimSet::new();
    let mut coefficients = vec![
        (
            COHESIVE_STIFFNESS_PROPERTY,
            interface_stiffness_per_area_dims(),
            1e9,
        ),
        (COHESIVE_STRENGTH_PROPERTY, Pressure::DIMS, 2e5),
        (
            COHESIVE_FRACTURE_ENERGY_PROPERTY,
            FRACTURE_ENERGY_PER_AREA_DIMS,
            100.0,
        ),
        (
            COHESIVE_COMPRESSION_STIFFNESS_PROPERTY,
            interface_stiffness_per_area_dims(),
            2e9,
        ),
    ];
    if viscous {
        coefficients.extend([
            (
                COHESIVE_NORMAL_VISCOSITY_PROPERTY,
                COHESIVE_VISCOSITY_DIMS,
                2e5,
            ),
            (
                COHESIVE_TANGENTIAL_VISCOSITY_PROPERTY,
                COHESIVE_VISCOSITY_DIMS,
                3e5,
            ),
        ]);
    }
    for (name, dims, value) in coefficients {
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value: PropertyValue::Scalar { value, dims },
                validity: ValidityDomain::unconstrained().with("T", 280.0, 320.0),
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::ConstantWithinValidity,
                observations: vec![],
                provenance: Provenance {
                    source: "manufactured finite-rotation cohesive fixture".into(),
                    license: "CC0-1.0".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    let surface = |name: &str| SurfaceSpec {
        material: MaterialStateId {
            chemistry: name.into(),
            phase: "solid".into(),
            process: "original native rigid fixture".into(),
            revision: 1,
        },
        texture_frame: name.into(),
    };
    let card = InterfaceSystemCard::assemble(
        surface("A"),
        surface("B"),
        SystemContext {
            medium: "original bonded interface".into(),
            third_body: None,
            environment: "fixed source query".into(),
            history: "explicit original bond".into(),
        },
        claims,
        vec![],
    )
    .unwrap();
    let point = QueryPoint::new().with("T", 300.0).unwrap();
    if viscous {
        resolve_viscous_bilinear_cohesive(
            &card,
            &point,
            MaterialPropertySelection::SingleClaimOnly,
            CohesiveEnvelope::CommonMode,
            0.05,
            10.0,
        )
        .unwrap()
    } else {
        resolve_bilinear_cohesive(
            &card,
            &point,
            MaterialPropertySelection::SingleClaimOnly,
            0.05,
        )
        .unwrap()
    }
}
fn config() -> RigidCohesiveConfig {
    RigidCohesiveConfig {
        maximum_sites: 64,
        maximum_iterations: 40,
        maximum_evaluations: 4096,
        maximum_material_evaluations: 65536,
        maximum_backtracks: 16,
        maximum_step_s: 0.02,
        maximum_site_force_n: 1000.0,
        maximum_total_force_n: 1e4,
        maximum_total_torque_n_m: 1000.0,
        maximum_jump_increment_m: 0.01,
        force_absolute_tolerance_n: 1e-11,
        torque_absolute_tolerance_n_m: 1e-12,
        relative_tolerance: 1e-11,
        energy_absolute_tolerance_j: 2e-11,
        energy_relative_tolerance: 1e-10,
        linear_momentum_absolute_tolerance_n_s: 1e-12,
        angular_momentum_absolute_tolerance_n_m_s: 1e-12,
        momentum_relative_tolerance: 1e-10,
        frame_tolerance: 1e-10,
        objective: ObjectiveCohesiveConfig {
            maximum_cayley_increment: 1.0,
            kinematic_tolerance_m: 2e-12,
            material_energy_tolerance_j_m2: 2e-7,
            spatial_work_tolerance_j: 2e-12,
        },
    }
}
fn limits() -> RigidCohesiveBodyLimits {
    RigidCohesiveBodyLimits {
        maximum_linear_speed_m_s: 100.0,
        maximum_angular_speed_rad_s: 1e4,
        maximum_kinetic_energy_j: 1e6,
        maximum_step_displacement_m: 0.1,
    }
}
fn norm(v: Vec3) -> f64 {
    v.x.hypot(v.y).hypot(v.z)
}
fn close(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance,
        "{a:e} != {b:e}, tolerance {tolerance:e}"
    );
}
fn close_vec(a: Vec3, b: Vec3, tolerance: f64) {
    close(norm(a.sub(b)), 0.0, tolerance);
}
fn vector(x: [f64; 3]) -> Vec3 {
    Vec3::new(x[0], x[1], x[2])
}

fn bodies(
    source: &ResolvedBilinearCohesive,
    rotation: UnitQuaternion,
    translation: Vec3,
    velocity: [Vec3; 2],
    omega: [Vec3; 2],
) -> [RigidCohesiveBody; 2] {
    core::array::from_fn(|i| {
        let properties = MassProperties::from_principal_frame(
            [0.4, 0.7][i],
            Vec3::ZERO,
            if i == 0 {
                Vec3::new(0.002, 0.003, 0.004)
            } else {
                Vec3::new(0.003, 0.004, 0.005)
            },
            UnitQuaternion::from_axis_angle(Vec3::new(1.0, 2.0, -1.0), [0.4, -0.3][i]).unwrap(),
        )
        .unwrap();
        let x = rotation
            .rotate_body_to_world(Vec3::new([-0.001, 0.001][i], 0.0, 0.0))
            .add(translation);
        let state = RigidBodyState::new(
            Pose::new(x, rotation).unwrap(),
            rotation
                .rotate_body_to_world(velocity[i])
                .scale(properties.mass()),
            properties.angular_momentum_body_checked(omega[i]).unwrap(),
        )
        .unwrap();
        RigidCohesiveBody::new(
            source.surfaces()[i].material.clone(),
            properties,
            state,
            limits(),
        )
        .unwrap()
    })
}

fn site<'a>(
    source: &'a ResolvedBilinearCohesive,
    bodies: &[RigidCohesiveBody; 2],
    arms: [Vec3; 2],
    area: f64,
    history: CohesivePoint<'a>,
) -> RigidCohesiveSite<'a> {
    let frame = bodies[0].state().pose().orientation();
    let d = bodies[1]
        .state()
        .pose()
        .point_world_from_body(arms[1])
        .unwrap()
        .sub(
            bodies[0]
                .state()
                .pose()
                .point_world_from_body(arms[0])
                .unwrap(),
        );
    let reference = frame.rotate_world_to_body(d).sub(vector(history.jump_m()));
    let original = ObjectiveCohesivePoint::new(history, frame, reference, area).unwrap();
    let bound = SourceBoundObjectiveCohesivePoint::new(
        original,
        source,
        [bodies[0].material(), bodies[1].material()],
    )
    .unwrap();
    RigidCohesiveSite::new(bound, arms, UnitQuaternion::IDENTITY).unwrap()
}
fn assert_frame(frame: &RigidCohesiveFrame) {
    assert!(frame.energy_residual_j.abs() <= frame.energy_tolerance_j);
    assert!(norm(frame.linear_momentum_residual_n_s) <= frame.linear_momentum_tolerance_n_s);
    assert!(norm(frame.angular_momentum_residual_n_m_s) <= frame.angular_momentum_tolerance_n_m_s);
    assert!(frame.fracture_dissipation_j >= 0.0 && frame.viscous_dissipation_j >= 0.0);
}

#[test]
fn g1_actual_two_body_axial_response_matches_independent_midpoint_oscillator() {
    let source = source(false);
    let initial_jump = 3e-5;
    let initial = source
        .law()
        .bonded_point()
        .trial([initial_jump, 0.0, 0.0], 1e-7)
        .unwrap()
        .point;
    let bodies = bodies(
        &source,
        UnitQuaternion::IDENTITY,
        Vec3::ZERO,
        [Vec3::new(-0.004, 0.0, 0.0), Vec3::new(0.006, 0.0, 0.0)],
        [Vec3::ZERO; 2],
    );
    let old_p = bodies[0]
        .state()
        .linear_momentum_world()
        .add(bodies[1].state().linear_momentum_world());
    let original_mass = bodies.each_ref().map(|body| body.properties());
    let sites = vec![
        site(&source, &bodies, [Vec3::ZERO; 2], 0.7e-6, initial),
        site(&source, &bodies, [Vec3::ZERO; 2], 1.3e-6, initial),
    ];
    let h = 2e-4;
    let mut system = RigidCohesiveSystem::new(
        bodies,
        sites,
        Gravity::ZERO,
        h,
        config(),
        &CancelGate::new(),
    )
    .unwrap();
    let initial_energy = system.total_energy_j().unwrap();
    let omega2 = 2000.0 * (1.0 / 0.4 + 1.0 / 0.7);
    let mut expected_g = initial_jump;
    let mut expected_v = 0.01;
    for _ in 0..24 {
        let next_g = ((1.0 - 0.25 * h * h * omega2) * expected_g + h * expected_v)
            / (1.0 + 0.25 * h * h * omega2);
        expected_v -= h * omega2 * 0.5 * (expected_g + next_g);
        expected_g = next_g;
        let frame = *system.step([RigidCohesiveLoad::ZERO; 2]).unwrap();
        assert_frame(&frame);
        close(
            system.sites()[0].point().numerical().material().jump_m()[0],
            expected_g,
            2e-15,
        );
        close(
            system.sites()[1].point().numerical().material().jump_m()[0],
            expected_g,
            2e-15,
        );
        let actual_v = system.bodies()[1]
            .state()
            .center_of_mass_velocity_world(original_mass[1])
            .unwrap()
            .sub(
                system.bodies()[0]
                    .state()
                    .center_of_mass_velocity_world(original_mass[0])
                    .unwrap(),
            );
        close(actual_v.x, expected_v, 2e-12);
        close_vec(
            system.bodies()[0]
                .state()
                .linear_momentum_world()
                .add(system.bodies()[1].state().linear_momentum_world()),
            old_p,
            2e-14,
        );
        close(system.total_energy_j().unwrap(), initial_energy, 1e-13);
        for (body, original) in system.bodies().iter().zip(original_mass) {
            assert_eq!(body.properties(), original);
        }
    }
}

#[test]
fn g3_rotating_heterogeneous_patch_is_objective_and_keeps_full_force_moment_work() {
    let source = source(false);
    let rotation = UnitQuaternion::from_axis_angle(Vec3::new(1.0, -2.0, 0.4), 1.4).unwrap();
    let offset = Vec3::new(0.7, -0.4, 0.2);
    let velocity = [
        Vec3::new(-0.002, 0.001, -0.0005),
        Vec3::new(0.001, -0.0003, 0.0007),
    ];
    let omega = [Vec3::new(2.0, -0.5, 4.0), Vec3::new(-0.4, 0.8, -3.0)];
    let build = |rotation, translation| {
        let body = bodies(&source, rotation, translation, velocity, omega);
        let points = [(-0.012, 0.7e-6), (0.0, 1.1e-6), (0.009, 1.6e-6)]
            .into_iter()
            .map(|(y, area)| {
                let history = source
                    .law()
                    .bonded_point()
                    .trial([2e-5, 1e-5, -0.5e-5], 1e-7)
                    .unwrap()
                    .point;
                site(
                    &source,
                    &body,
                    [Vec3::new(0.0001, y, 0.002), Vec3::new(-0.0001, y, 0.002)],
                    area,
                    history,
                )
            })
            .collect();
        RigidCohesiveSystem::new(
            body,
            points,
            Gravity::ZERO,
            0.003,
            config(),
            &CancelGate::new(),
        )
        .unwrap()
    };
    let mut original = build(UnitQuaternion::IDENTITY, Vec3::ZERO);
    let mut transformed = build(rotation, offset);
    let mut had_couple = false;
    let mut differentiated_damage = false;
    for _ in 0..20 {
        let frame_a = *original.step([RigidCohesiveLoad::ZERO; 2]).unwrap();
        let frame_b = *transformed.step([RigidCohesiveLoad::ZERO; 2]).unwrap();
        assert_frame(&frame_a);
        assert_frame(&frame_b);
        close(
            frame_a.mechanical_energy_j,
            frame_b.mechanical_energy_j,
            2e-10,
        );
        for i in 0..2 {
            let a = original.bodies()[i].state();
            let b = transformed.bodies()[i].state();
            close_vec(
                b.pose().position_world(),
                rotation
                    .rotate_body_to_world(a.pose().position_world())
                    .add(offset),
                2e-10,
            );
            close_vec(
                b.linear_momentum_world(),
                rotation.rotate_body_to_world(a.linear_momentum_world()),
                2e-9,
            );
            close_vec(b.angular_momentum_body(), a.angular_momentum_body(), 2e-10);
            for basis in [
                Vec3::new(1.0, 0.0, 0.0),
                Vec3::new(0.0, 1.0, 0.0),
                Vec3::new(0.0, 0.0, 1.0),
            ] {
                close_vec(
                    b.pose().orientation().rotate_body_to_world(basis),
                    rotation
                        .rotate_body_to_world(a.pose().orientation().rotate_body_to_world(basis)),
                    2e-9,
                );
            }
        }
        for (a, b) in original.sites().iter().zip(transformed.sites()) {
            close_vec(
                vector(a.point().numerical().material().jump_m()),
                vector(b.point().numerical().material().jump_m()),
                2e-10,
            );
            close(
                a.point().numerical().material().damage(),
                b.point().numerical().material().damage(),
                2e-8,
            );
            assert!(core::ptr::eq(a.point().source(), &source));
        }
        had_couple |= original
            .site_frames()
            .iter()
            .flatten()
            .any(|f| norm(f.frame_couple_on_a_n_m) > 1e-7);
        let damage: Vec<_> = original
            .sites()
            .iter()
            .map(|site| site.point().numerical().material().damage())
            .collect();
        differentiated_damage |= (damage[0] - damage[2]).abs() > 0.05;
    }
    assert!(had_couple && differentiated_damage);
    assert!(
        original
            .sites()
            .iter()
            .any(|site| site.point().numerical().material().damage() > 0.8)
    );
}

#[test]
fn g1_original_viscosity_and_progressive_peeling_keep_separate_physical_losses() {
    let source = source(true);
    let body = bodies(
        &source,
        UnitQuaternion::IDENTITY,
        Vec3::ZERO,
        [Vec3::new(-0.08, 0.0, 0.0), Vec3::new(0.06, 0.0, 0.0)],
        [Vec3::new(0.2, 0.1, 0.6), Vec3::new(-0.1, 0.2, 2.0)],
    );
    let points = [-0.01, 0.0, 0.01]
        .into_iter()
        .map(|y| {
            site(
                &source,
                &body,
                [Vec3::new(0.0, y, 0.001); 2],
                1e-6,
                source.law().bonded_point(),
            )
        })
        .collect();
    let mut system = RigidCohesiveSystem::new(
        body,
        points,
        Gravity::ZERO,
        5e-4,
        config(),
        &CancelGate::new(),
    )
    .unwrap();
    let initial_energy = system.total_energy_j().unwrap();
    let mut fracture = 0.0;
    let mut viscous = 0.0;
    let mut previous_damage = [0.0; 3];
    let mut uneven = false;
    for _ in 0..65 {
        let frame = *system.step([RigidCohesiveLoad::ZERO; 2]).unwrap();
        assert_frame(&frame);
        fracture += frame.fracture_dissipation_j;
        viscous += frame.viscous_dissipation_j;
        let mut site_fracture = 0.0;
        let mut site_viscous = 0.0;
        for (i, value) in system.site_frames().iter().enumerate() {
            let value = value.unwrap();
            assert!(value.damage >= previous_damage[i]);
            previous_damage[i] = value.damage;
            site_fracture += value.fracture_dissipation_j;
            site_viscous += value.viscous_dissipation_j;
        }
        close(site_fracture, frame.fracture_dissipation_j, 1e-18);
        close(site_viscous, frame.viscous_dissipation_j, 1e-18);
        uneven |= (previous_damage[0] - previous_damage[2]).abs() > 0.03;
        close(
            system.total_energy_j().unwrap() + fracture + viscous,
            initial_energy,
            2e-10,
        );
    }
    assert!(uneven);
    assert!(previous_damage.into_iter().all(|d| d == 1.0));
    close(fracture, 3e-6 * source.law().fracture_energy_j_m2(), 1e-13);
    assert!(viscous > 1e-7);
    assert_eq!(system.steps(), 65);
}

#[test]
fn g1_failed_original_bonds_preserve_full_native_free_motion_and_reclosure_history() {
    let source = source(false);
    let failed = source
        .law()
        .bonded_point()
        .trial([0.0012, 0.0, 0.0], 1e-7)
        .unwrap()
        .point;
    assert_eq!(failed.damage(), 1.0);
    let opened = failed.trial([0.0005, 0.0, 0.0], 1e-7).unwrap().point;
    let body = bodies(
        &source,
        UnitQuaternion::IDENTITY,
        Vec3::ZERO,
        [
            Vec3::new(0.002, -0.001, 0.003),
            Vec3::new(0.004, 0.002, -0.001),
        ],
        [Vec3::new(1.0, 0.5, -0.4), Vec3::new(-0.6, 0.7, 0.8)],
    );
    let points = vec![site(&source, &body, [Vec3::ZERO; 2], 2e-6, opened)];
    let mut expected = body.each_ref().map(|b| b.state());
    let properties = body.each_ref().map(|b| b.properties());
    let mut free = RigidCohesiveSystem::new(
        body,
        points,
        Gravity::ZERO,
        5e-4,
        config(),
        &CancelGate::new(),
    )
    .unwrap();
    let native = RigidBodyIntegrator::new(Gravity::ZERO);
    for _ in 0..20 {
        for i in 0..2 {
            expected[i] = native
                .step_cayley(expected[i], properties[i], Wrench::ZERO, 5e-4)
                .unwrap()
                .state_after;
        }
        let frame = *free.step([RigidCohesiveLoad::ZERO; 2]).unwrap();
        assert_eq!(free.bodies().each_ref().map(|b| b.state()), expected);
        assert_eq!(frame.contact_loads, [RigidCohesiveLoad::ZERO; 2]);
        assert_eq!(frame.fracture_dissipation_j, 0.0);
        assert_eq!(free.sites()[0].point().numerical().material().damage(), 1.0);
    }
    let closed = failed.trial([-2e-5, 0.0, 0.0], 1e-7).unwrap().point;
    let body = bodies(
        &source,
        UnitQuaternion::IDENTITY,
        Vec3::ZERO,
        [Vec3::new(0.02, 0.0, 0.0), Vec3::new(-0.02, 0.0, 0.0)],
        [Vec3::ZERO; 2],
    );
    let points = vec![site(&source, &body, [Vec3::ZERO; 2], 2e-6, closed)];
    let mut reclosure = RigidCohesiveSystem::new(
        body,
        points,
        Gravity::ZERO,
        1e-4,
        config(),
        &CancelGate::new(),
    )
    .unwrap();
    let initial_energy = reclosure.total_energy_j().unwrap();
    let mut compressed = false;
    let mut released = false;
    for _ in 0..350 {
        let frame = *reclosure.step([RigidCohesiveLoad::ZERO; 2]).unwrap();
        assert_frame(&frame);
        let point = reclosure.sites()[0].point().numerical().material();
        assert_eq!(point.damage(), 1.0);
        assert_eq!(frame.fracture_dissipation_j, 0.0);
        assert_eq!(frame.viscous_dissipation_j, 0.0);
        compressed |= point.jump_m()[0] < 0.0 && frame.cohesive_stored_energy_j > 0.0;
        released |= point.jump_m()[0] > 0.0 && frame.cohesive_stored_energy_j == 0.0;
        close(reclosure.total_energy_j().unwrap(), initial_energy, 2e-12);
    }
    assert!(compressed && released);
}

#[test]
fn g1_external_world_force_couple_and_original_gravity_share_actual_work() {
    let source = source(true);
    let body = bodies(
        &source,
        UnitQuaternion::from_axis_angle(Vec3::new(1.0, 2.0, -1.0), 0.4).unwrap(),
        Vec3::new(0.2, -0.1, 0.3),
        [Vec3::new(0.002, -0.003, 0.001); 2],
        [Vec3::new(0.4, -0.2, 0.3), Vec3::new(-0.1, 0.6, 0.2)],
    );
    let points = vec![site(
        &source,
        &body,
        [Vec3::new(0.0, 0.01, 0.002); 2],
        1e-6,
        source.law().bonded_point(),
    )];
    let gravity = Gravity::new(Vec3::new(0.0, 0.0, -9.81)).unwrap();
    let mut system =
        RigidCohesiveSystem::new(body, points, gravity, 1e-3, config(), &CancelGate::new())
            .unwrap();
    let initial_energy = system.total_energy_j().unwrap();
    let loads = [
        RigidCohesiveLoad {
            force_world_n: Vec3::new(0.2, -0.1, 0.3),
            moment_world_n_m: Vec3::new(0.03, -0.02, 0.01),
        },
        RigidCohesiveLoad {
            force_world_n: Vec3::new(-0.05, 0.15, 0.2),
            moment_world_n_m: Vec3::new(-0.01, 0.03, 0.02),
        },
    ];
    let mut supplied = 0.0;
    let mut dissipated = 0.0;
    for _ in 0..20 {
        let old = system.bodies().each_ref().map(|b| b.state());
        let frame = *system.step(loads).unwrap();
        assert_frame(&frame);
        let mut independent_work = 0.0;
        for i in 0..2 {
            let after = system.bodies()[i].state();
            let lmid = old[i]
                .angular_momentum_body()
                .scale(0.5)
                .add(after.angular_momentum_body().scale(0.5));
            let c = old[i].pose().orientation().rotate_body_to_world(
                system.bodies()[i]
                    .properties()
                    .angular_velocity_body(lmid)
                    .scale(1e-3),
            );
            independent_work += loads[i].force_world_n.dot(
                after
                    .pose()
                    .position_world()
                    .sub(old[i].pose().position_world()),
            ) + loads[i].moment_world_n_m.dot(c);
        }
        close(frame.external_work_j, independent_work, 2e-16);
        supplied += independent_work;
        dissipated += frame.fracture_dissipation_j + frame.viscous_dissipation_j;
        close(
            system.total_energy_j().unwrap() + dissipated - supplied,
            initial_energy,
            2e-10,
        );
    }
}

#[test]
fn g4_late_body_guard_and_cancellation_preserve_bodies_bonds_observers_clock_and_retry() {
    let source = source(false);
    let build = || {
        let mut body = bodies(
            &source,
            UnitQuaternion::IDENTITY,
            Vec3::ZERO,
            [Vec3::ZERO; 2],
            [Vec3::ZERO; 2],
        );
        let b = &body[0];
        let mut restricted = b.limits();
        restricted.maximum_kinetic_energy_j = 1e-12;
        body[0] =
            RigidCohesiveBody::new(b.material().clone(), b.properties(), b.state(), restricted)
                .unwrap();
        let points = vec![site(
            &source,
            &body,
            [Vec3::ZERO; 2],
            1e-6,
            source.law().bonded_point(),
        )];
        RigidCohesiveSystem::new(
            body,
            points,
            Gravity::ZERO,
            1e-4,
            config(),
            &CancelGate::new(),
        )
        .unwrap()
    };
    let mut system = build();
    let original_bodies = system.bodies().clone();
    let original_sites = system.sites().to_vec();
    let cancelled = CancelGate::new();
    cancelled.request();
    assert!(matches!(
        system.step_under_gate([RigidCohesiveLoad::ZERO; 2], &cancelled),
        Err(RigidCohesiveError::Cancelled)
    ));
    let invalid = [
        RigidCohesiveLoad {
            force_world_n: Vec3::new(100.0, 0.0, 0.0),
            ..RigidCohesiveLoad::ZERO
        },
        RigidCohesiveLoad::ZERO,
    ];
    let refusal = system.step(invalid);
    assert!(
        matches!(
            refusal,
            Err(RigidCohesiveError::Range {
                what: "body kinetic energy",
                ..
            })
        ),
        "{refusal:?}"
    );
    assert_eq!(system.bodies(), &original_bodies);
    assert_eq!(system.sites(), &original_sites);
    assert!(system.site_frames().iter().all(Option::is_none));
    assert_eq!(system.steps(), 0);
    assert_eq!(system.time_s(), 0.0);
    assert!(system.last_frame().is_none());
    let valid = [
        RigidCohesiveLoad {
            force_world_n: Vec3::new(1e-6, 0.0, 0.0),
            ..RigidCohesiveLoad::ZERO
        },
        RigidCohesiveLoad::ZERO,
    ];
    let actual = *system.step(valid).unwrap();
    let mut fresh = build();
    assert_eq!(actual, *fresh.step(valid).unwrap());
    assert_eq!(system.bodies(), fresh.bodies());
    assert_eq!(system.sites(), fresh.sites());
    assert_eq!(system.site_frames(), fresh.site_frames());
}

#[test]
fn g0_g4_original_source_geometry_work_budget_and_strict_energy_are_enforced() {
    let source = source(false);
    let build = |cfg, preload| {
        let body = bodies(
            &source,
            UnitQuaternion::IDENTITY,
            Vec3::ZERO,
            [Vec3::ZERO; 2],
            [Vec3::ZERO; 2],
        );
        let point = source
            .law()
            .bonded_point()
            .trial([preload, 0.0, 0.0], 1e-7)
            .unwrap()
            .point;
        let points = vec![site(&source, &body, [Vec3::ZERO; 2], 1e-6, point)];
        RigidCohesiveSystem::new(body, points, Gravity::ZERO, 1e-4, cfg, &CancelGate::new())
            .unwrap()
    };
    let mut budget = config();
    budget.maximum_evaluations = 1;
    let mut system = build(budget, 0.0);
    let before = system.bodies().clone();
    let load = [
        RigidCohesiveLoad {
            force_world_n: Vec3::new(0.1, 0.0, 0.0),
            ..RigidCohesiveLoad::ZERO
        },
        RigidCohesiveLoad::ZERO,
    ];
    assert!(matches!(
        system.step(load),
        Err(RigidCohesiveError::Budget(
            "complete rigid patch evaluations"
        ))
    ));
    assert_eq!(system.bodies(), &before);
    assert_eq!(system.steps(), 0);
    let mut fresh = build(budget, 0.0);
    assert_eq!(
        *system.step([RigidCohesiveLoad::ZERO; 2]).unwrap(),
        *fresh.step([RigidCohesiveLoad::ZERO; 2]).unwrap()
    );

    let mut strict = config();
    strict.energy_absolute_tolerance_j = 1e-32;
    strict.energy_relative_tolerance = 1e-32;
    let mut system = build(strict, 3.7e-5);
    let before = system.bodies().clone();
    let points = system.sites().to_vec();
    assert!(matches!(
        system.step(load),
        Err(RigidCohesiveError::Balance {
            what: "actual rigid/cohesive energy",
            ..
        })
    ));
    assert_eq!(system.bodies(), &before);
    assert_eq!(system.sites(), &points);
    assert_eq!(system.steps(), 0);

    let body = bodies(
        &source,
        UnitQuaternion::IDENTITY,
        Vec3::ZERO,
        [Vec3::ZERO; 2],
        [Vec3::ZERO; 2],
    );
    let point = site(
        &source,
        &body,
        [Vec3::ZERO; 2],
        1e-6,
        source.law().bonded_point(),
    );
    let wrong = [
        RigidCohesiveBody::new(
            body[1].material().clone(),
            body[0].properties(),
            body[0].state(),
            limits(),
        )
        .unwrap(),
        RigidCohesiveBody::new(
            body[0].material().clone(),
            body[1].properties(),
            body[1].state(),
            limits(),
        )
        .unwrap(),
    ];
    assert!(matches!(
        RigidCohesiveSystem::new(
            wrong,
            vec![point],
            Gravity::ZERO,
            1e-4,
            config(),
            &CancelGate::new()
        ),
        Err(RigidCohesiveError::Source(_))
    ));
    let wrong_geometry = RigidCohesiveSite::new(
        point.point(),
        [Vec3::new(0.01, 0.0, 0.0), Vec3::ZERO],
        UnitQuaternion::IDENTITY,
    )
    .unwrap();
    assert!(matches!(
        RigidCohesiveSystem::new(
            body,
            vec![wrong_geometry],
            Gravity::ZERO,
            1e-4,
            config(),
            &CancelGate::new()
        ),
        Err(RigidCohesiveError::Balance {
            what: "initial attachment separation",
            ..
        })
    ));
}
