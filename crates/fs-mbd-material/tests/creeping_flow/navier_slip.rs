//! Physical interface slip reaches the original native force, work and motion.
use super::navier_slip_source as evidence;
use super::*;
use fs_mbd::{RigidBodyState, Wrench};
use fs_mbd_material::creeping_flow::navier_slip::NavierSlipBoundary;
use fs_mbd_material::creeping_flow::unsteady::BassetResolution;

fn initial(body: &fs_mbd_material::principal::PrincipalBody<'_>) -> RigidBodyState {
    body.state_from_kinematics(
        Se3::identity(),
        Twist::new(Vec3::new(0.0, 0.0, 1.0), Vec3::new(0.01, -0.02, 0.03)),
    )
    .unwrap()
}

#[test]
fn sourced_zero_slip_preserves_original_load_native_step_and_source_bits() {
    let source = sphere(20000.0);
    let body = source.principal_body(1e-12, || false).unwrap();
    let material = fluid(1.0, 1.0);
    let interface = evidence::constant(source.regions()[0].material().material(), 0.0, 400.0);
    let boundary = NavierSlipBoundary::new(&interface).unwrap();
    let make = || {
        CreepingSphere::new(
            &material,
            RADIUS,
            gravity(),
            fs_mbd::Vec3::new(0.004, 0.0, 0.0),
            limits(),
        )
        .unwrap()
    };
    let ordinary = make();
    let explicit = make().with_navier_slip(boundary).unwrap();
    assert!(core::ptr::eq(explicit.source(), &material));
    assert!(core::ptr::eq(
        explicit.navier_slip().unwrap().source(),
        &interface
    ));
    let state = initial(&body);
    assert_eq!(
        ordinary.load(&body, state, || false).unwrap(),
        explicit.load(&body, state, || false).unwrap()
    );
    assert_eq!(
        ordinary
            .step(&body, state, controls(0.01), || false)
            .unwrap(),
        explicit
            .step(&body, state, controls(0.01), || false)
            .unwrap()
    );
}

#[test]
fn sourced_zero_slip_preserves_selected_lift_in_both_builder_orders() {
    let source = sphere(20000.0);
    let body = source.principal_body(1e-12, || false).unwrap();
    let material = fluid(1.0, 1.0);
    let interface = evidence::constant(source.regions()[0].material().material(), 0.0, 400.0);
    let boundary = NavierSlipBoundary::new(&interface).unwrap();
    let make = || {
        CreepingSphere::new(
            &material,
            RADIUS,
            gravity(),
            fs_mbd::Vec3::new(0.004, 0.0, 0.0),
            limits(),
        )
        .unwrap()
    };
    let ordinary = make().with_rotational_lift();
    let lift_first = make()
        .with_rotational_lift()
        .with_navier_slip(boundary)
        .unwrap();
    let slip_first = make()
        .with_navier_slip(boundary)
        .unwrap()
        .with_rotational_lift();
    let state = initial(&body);
    let expected_load = ordinary.load(&body, state, || false).unwrap();
    assert_ne!(expected_load.lift_force_world_n, fs_mbd::Vec3::ZERO);
    let expected_step = ordinary
        .step(&body, state, controls(0.01), || false)
        .unwrap();
    for selected in [&lift_first, &slip_first] {
        assert!(selected.rotational_lift_enabled());
        assert!(core::ptr::eq(
            selected.navier_slip().unwrap().source(),
            &interface
        ));
        assert_eq!(
            selected.load(&body, state, || false).unwrap(),
            expected_load
        );
        assert_eq!(
            selected
                .step(&body, state, controls(0.01), || false)
                .unwrap(),
            expected_step
        );
    }
}

#[test]
fn positive_slip_and_selected_lift_refuse_in_both_builder_orders() {
    let source = sphere(20000.0);
    let body = source.principal_body(1e-12, || false).unwrap();
    let material = fluid(1.0, 1.0);
    let interface = evidence::constant(source.regions()[0].material().material(), RADIUS, 400.0);
    let boundary = NavierSlipBoundary::new(&interface).unwrap();
    let make = || {
        CreepingSphere::new(
            &material,
            RADIUS,
            Gravity::ZERO,
            fs_mbd::Vec3::ZERO,
            limits(),
        )
        .unwrap()
    };
    assert!(matches!(
        make().with_rotational_lift().with_navier_slip(boundary),
        Err(CreepingFlowError::Invalid(_))
    ));
    let selected = make()
        .with_navier_slip(boundary)
        .unwrap()
        .with_rotational_lift();
    assert!(selected.rotational_lift_enabled());
    let state = initial(&body);
    assert!(matches!(
        selected.load(&body, state, || false),
        Err(CreepingFlowError::Invalid(_))
    ));
    assert!(matches!(
        selected.step(&body, state, controls(0.01), || false),
        Err(CreepingFlowError::Invalid(_))
    ));
    assert!(matches!(
        selected.step_with_wrench(&body, state, Wrench::ZERO, controls(0.01), || false),
        Err(CreepingFlowError::Invalid(_))
    ));
}

#[test]
fn finite_interface_slip_changes_both_resistances_and_total_loss_without_changing_buoyancy() {
    let source = sphere(20000.0);
    let body = source.principal_body(1e-12, || false).unwrap();
    let material = fluid(1.0, 1.0);
    let interface = evidence::constant(source.regions()[0].material().material(), RADIUS, 400.0);
    let make =
        || CreepingSphere::new(&material, RADIUS, gravity(), fs_mbd::Vec3::ZERO, limits()).unwrap();
    let no_slip = make().load(&body, initial(&body), || false).unwrap();
    let slip = make()
        .with_navier_slip(NavierSlipBoundary::new(&interface).unwrap())
        .unwrap()
        .load(&body, initial(&body), || false)
        .unwrap();
    close(
        slip.drag_force_world_n.x,
        0.75 * no_slip.drag_force_world_n.x,
    );
    close(
        slip.rotary_torque_body_nm.z,
        0.25 * no_slip.rotary_torque_body_nm.z,
    );
    assert_eq!(slip.buoyancy_force_world_n, no_slip.buoyancy_force_world_n);
    assert_eq!(slip.translation_reynolds, no_slip.translation_reynolds);
    assert_eq!(slip.rotation_reynolds, no_slip.rotation_reynolds);
    let velocity = initial(&body)
        .center_of_mass_velocity_world(body.properties())
        .unwrap();
    close(
        slip.dissipated_power_w,
        -slip.drag_force_world_n.dot(velocity) - slip.rotary_torque_body_nm.z,
    );
    assert!(slip.dissipated_power_w > 0.0 && slip.dissipated_power_w < no_slip.dissipated_power_w);
}

#[test]
fn driven_partial_slip_uses_one_native_midpoint_and_closes_all_physical_work() {
    let source = sphere(20000.0);
    let body = source.principal_body(1e-12, || false).unwrap();
    let material = fluid(1.0, 1.0);
    let interface = evidence::constant(
        source.regions()[0].material().material(),
        2.0 * RADIUS,
        400.0,
    );
    let ambient = fs_mbd::Vec3::new(0.004, -0.005, 0.002);
    let medium = CreepingSphere::new(&material, RADIUS, Gravity::ZERO, ambient, limits())
        .unwrap()
        .with_navier_slip(NavierSlipBoundary::new(&interface).unwrap())
        .unwrap();
    let state = initial(&body);
    let applied = Wrench {
        force_world: fs_mbd::Vec3::new(1e-4, -2e-4, 3e-4),
        torque_body: fs_mbd::Vec3::new(0.0, 0.0, 1e-6),
    };
    let dt = 0.01;
    let result = medium
        .step_with_wrench(&body, state, applied, controls(dt), || false)
        .unwrap();
    let native = RigidBodyIntegrator::new(Gravity::ZERO)
        .step(
            state,
            body.properties(),
            Wrench {
                force_world: result
                    .load_midpoint
                    .wrench
                    .force_world
                    .add(applied.force_world),
                torque_body: result
                    .load_midpoint
                    .wrench
                    .torque_body
                    .add(applied.torque_body),
            },
            dt,
        )
        .unwrap();
    assert_eq!(result.dynamics, native);
    let mass = body.properties().mass();
    let inertia = body.properties().principal_inertia_body().z;
    let c = 6.0 * core::f64::consts::PI * RADIUS * 5.0 / 7.0;
    let rotary = 8.0 * core::f64::consts::PI * RADIUS.powi(3) / 7.0;
    let v0 = state
        .center_of_mass_velocity_world(body.properties())
        .unwrap();
    let v1 = native
        .state_after
        .center_of_mass_velocity_world(body.properties())
        .unwrap();
    let midpoint = (v0.x - ambient.x + 0.5 * dt * applied.force_world.x / mass)
        / (1.0 + 0.5 * dt * c / mass)
        + ambient.x;
    close(v1.x, 2.0 * midpoint - v0.x);
    let angular_mid =
        (1.0 + 0.5 * dt * applied.torque_body.z / inertia) / (1.0 + 0.5 * dt * rotary / inertia);
    close(
        native.state_after.angular_momentum_body().z,
        inertia + dt * (applied.torque_body.z - rotary * angular_mid),
    );
    close(
        result.viscous_heat_j,
        dt * result.load_midpoint.dissipated_power_w,
    );
    close(
        result.background_work_j,
        dt * result.load_midpoint.drag_force_world_n.dot(ambient),
    );
    close(
        result.stored_energy_change_j + result.viscous_heat_j
            - result.background_work_j
            - result.applied_work_j,
        0.0,
    );
    assert!(result.viscous_heat_j > 0.0);
}

#[test]
fn diffusion_and_relaxation_use_slip_resistances_but_reynolds_keeps_physical_radius() {
    let material = fluid(1.0, 1.0);
    for density in [100.0, 20000.0] {
        let source = sphere(density);
        let body = source.principal_body(1e-12, || false).unwrap();
        let interface =
            evidence::constant(source.regions()[0].material().material(), RADIUS, 400.0);
        let make = || {
            CreepingSphere::new(
                &material,
                RADIUS,
                Gravity::ZERO,
                fs_mbd::Vec3::ZERO,
                limits(),
            )
            .unwrap()
        };
        let no_slip = make();
        let slip = make()
            .with_navier_slip(NavierSlipBoundary::new(&interface).unwrap())
            .unwrap();
        let rest = body
            .state_from_kinematics(Se3::identity(), Twist::zero())
            .unwrap();
        let dt = if density == 100.0 { 0.001 } else { 0.1 };
        assert!(no_slip.step(&body, rest, controls(dt), || false).is_err());
        let accepted = slip.step(&body, rest, controls(dt), || false).unwrap();
        close(accepted.diffusion_ratio, 15.0 / density / 4.0);
        let fast = body
            .state_from_kinematics(
                Se3::identity(),
                Twist::new(Vec3::new(0.0, 0.0, 0.0), Vec3::new(10.0, 0.0, 0.0)),
            )
            .unwrap();
        assert!(matches!(
            slip.load(&body, fast, || false),
            Err(CreepingFlowError::Reynolds { .. })
        ));
    }
}

#[test]
fn no_slip_history_and_lift_cannot_be_silently_reused_for_partial_slip() {
    let source = sphere(20000.0);
    let body = source.principal_body(1e-12, || false).unwrap();
    let material = fluid(1.0, 1.0);
    let state = body
        .state_from_kinematics(Se3::identity(), Twist::zero())
        .unwrap();
    let resolution = BassetResolution {
        minimum_wavenumber_radius: 0.01,
        maximum_wavenumber_radius: 10.0,
        intervals: 8,
    };
    for length in [0.0, RADIUS] {
        let interface =
            evidence::constant(source.regions()[0].material().material(), length, 400.0);
        let medium = CreepingSphere::new(
            &material,
            RADIUS,
            Gravity::ZERO,
            fs_mbd::Vec3::ZERO,
            limits(),
        )
        .unwrap()
        .with_navier_slip(NavierSlipBoundary::new(&interface).unwrap())
        .unwrap();
        assert_eq!(
            medium
                .unsteady_motion(&body, state, resolution, || false)
                .is_ok(),
            length == 0.0
        );
        assert_eq!(
            medium
                .unsteady_rotating_motion(&body, state, resolution, || false)
                .is_ok(),
            length == 0.0
        );
        assert_eq!(
            medium
                .load_with_rotational_lift(&body, state, || false)
                .is_ok(),
            length == 0.0
        );
        assert_eq!(
            medium
                .step_with_rotational_lift(&body, state, Wrench::ZERO, controls(0.01), || false)
                .is_ok(),
            length == 0.0
        );
    }
}

#[test]
fn typed_temperature_roles_are_checked_for_both_interface_and_solid() {
    let reference = sphere(20000.0);
    let solid = reference.regions()[0].material().material();
    let material = fluid(1.0, 1.0);
    let supported = evidence::constant(solid, RADIUS, 400.0);
    let cases = [
        (None, true),
        (
            Some(QuantitySpec::dimensional(fs_qty::Temperature::DIMS)),
            true,
        ),
        (
            Some(QuantitySpec::semantic(SemanticType::new(
                QuantityKind::AbsoluteTemperature,
                ValueForm::Static,
            ))),
            true,
        ),
        (
            Some(QuantitySpec::semantic(SemanticType::new(
                QuantityKind::Pressure,
                ValueForm::Static,
            ))),
            false,
        ),
        (
            Some(QuantitySpec::semantic(SemanticType::new(
                QuantityKind::TemperatureDifference,
                ValueForm::Static,
            ))),
            false,
        ),
    ];
    for (quantity, admitted) in cases {
        let (validity, query) = match quantity {
            None => (
                ValidityDomain::unconstrained().with("T", 200.0, 400.0),
                point(300.0),
            ),
            Some(quantity) => (
                ValidityDomain::unconstrained().with_quantity("T", quantity, 200.0, 400.0),
                fs_matdb::QueryPoint::new()
                    .with_quantity("T", quantity, 300.0)
                    .unwrap(),
            ),
        };
        let interface = evidence::resolve(
            solid,
            RADIUS,
            fs_qty::Length::DIMS,
            validity.clone(),
            &query,
        );
        assert_eq!(
            NavierSlipBoundary::new(&interface).is_ok(),
            admitted,
            "interface {quantity:?}"
        );

        // The ordinary material resolver permits other physical coordinates
        // under arbitrary names. Only the Navier consumer interprets T as kelvin.
        let mut claims = ClaimSet::new();
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new("density", Density::DIMS),
                value: PropertyValue::Scalar {
                    value: 20000.0,
                    dims: Density::DIMS,
                },
                validity,
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::ConstantWithinValidity,
                observations: vec![],
                provenance: Provenance {
                    source: "synthetic typed solid coordinate".into(),
                    license: "CC0-1.0".into(),
                    artifact: None,
                },
            })
            .unwrap();
        let card = MaterialCard::assemble(solid.clone(), claims, vec![]).unwrap();
        let source = resolve_body(
            &[RegionBinding {
                key: "sphere",
                geometry: VolumeMoments::sphere(RADIUS).unwrap(),
                material: &card,
                point: &query,
                selection: MaterialPropertySelection::SingleClaimOnly,
            }],
            1,
            || false,
        )
        .unwrap();
        let body = source.principal_body(1e-12, || false).unwrap();
        let medium = CreepingSphere::new(
            &material,
            RADIUS,
            Gravity::ZERO,
            fs_mbd::Vec3::ZERO,
            limits(),
        )
        .unwrap()
        .with_navier_slip(NavierSlipBoundary::new(&supported).unwrap())
        .unwrap();
        assert_eq!(
            medium.load(&body, initial(&body), || false).is_ok(),
            admitted,
            "solid {quantity:?}"
        );
    }
}

#[test]
fn source_dimensions_passivity_static_support_and_solid_identity_are_enforced() {
    let source = sphere(20000.0);
    let solid = source.regions()[0].material().material();
    for (length, dims) in [
        (-RADIUS, fs_qty::Length::DIMS),
        (RADIUS, fs_qty::Dims::NONE),
    ] {
        let interface = evidence::resolve(
            solid,
            length,
            dims,
            ValidityDomain::unconstrained().with("T", 200.0, 400.0),
            &point(300.0),
        );
        assert!(NavierSlipBoundary::new(&interface).is_err());
    }
    let interface = evidence::resolve(
        solid,
        RADIUS,
        fs_qty::Length::DIMS,
        ValidityDomain::unconstrained()
            .with("T", 200.0, 400.0)
            .with("shear_rate", 0.0, 1.0),
        &point(300.0).with("shear_rate", 0.5).unwrap(),
    );
    assert!(NavierSlipBoundary::new(&interface).is_err());
    let mut foreign = solid.clone();
    foreign.revision += 1;
    let interface = evidence::constant(&foreign, RADIUS, 400.0);
    let material = fluid(1.0, 1.0);
    let medium = CreepingSphere::new(
        &material,
        RADIUS,
        Gravity::ZERO,
        fs_mbd::Vec3::ZERO,
        limits(),
    )
    .unwrap()
    .with_navier_slip(NavierSlipBoundary::new(&interface).unwrap())
    .unwrap();
    let body = source.principal_body(1e-12, || false).unwrap();
    assert!(medium.load(&body, initial(&body), || false).is_err());
}
