//! Same-source compliant pads in the actual articulated force and motion owners.
mod common;
use common::{close, density_card, point};
use fs_evidence::ValidityDomain;
use fs_ga::{Se3, So3, So3Tangent, Twist, Vec3, Wrench};
use fs_matdb::{
    ClaimSet, InterfaceSystemCard, InterpolationPolicy, MaterialCard, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, SurfaceSpec, SystemContext, UncertaintyModel,
};
use fs_material::state_point::{
    MaterialPropertySelection, ResolvedInterfaceStatePoint, ScalarAdmissibility,
    ScalarPropertyRequirement, interface_stiffness_per_area_dims,
    interface_viscosity_per_area_dims, resolve_interface_state_point,
};
use fs_mbd::articulated::motion::spatial_loads::SpatialMotionLoads;
use fs_mbd::articulated::motion::spatial_loads::bushing::{
    MaterialFrame, step_free_floating_bushings_implicit,
};
use fs_mbd::articulated::motion::{FreeFloatingMotionLoads, FreeFloatingMotionState};
use fs_mbd::articulated::{BaseState, FreeFloatingBaseState, JointModel, forward_kinematics};
use fs_mbd_material::articulated::{MaterialLink, ResolvedArticulated, resolve_articulated};
use fs_mbd_material::attachments::{
    AttachmentError, MaterialBushing, RectangularPad, step_bushings_implicit,
};
use fs_mbd_material::{RegionBinding, VolumeMoments};
use fs_qty::Temperature;

const ZERO: Vec3 = Vec3 {
    x: 0.0,
    y: 0.0,
    z: 0.0,
};
const PROPERTIES: [(&str, &str); 3] = [("kx", "cx"), ("ky", "cy"), ("kz", "cz")];

fn mechanism(a: &MaterialCard, b: &MaterialCard, joint: JointModel) -> ResolvedArticulated {
    let p = point(300.0);
    let regions = [a, b].map(|material| {
        [RegionBinding {
            key: "surface",
            material,
            point: &p,
            geometry: VolumeMoments::cuboid(Vec3::new(1.0, 1.0, 1.0)).unwrap(),
            selection: MaterialPropertySelection::SingleClaimOnly,
        }]
    });
    // Ideal rigid carriers with coincident attachment origins; collision
    // geometry is not part of this independent spring/mass reference.
    resolve_articulated(
        &[
            MaterialLink {
                name: "a",
                parent: None,
                parent_from_child_zero: Se3::identity(),
                joint: JointModel::FIXED,
                regions: &regions[0],
            },
            MaterialLink {
                name: "b",
                parent: Some(0),
                parent_from_child_zero: Se3::identity(),
                joint,
                regions: &regions[1],
            },
        ],
        2,
        || false,
    )
    .unwrap()
}

fn interface(a: &MaterialCard, b: &MaterialCard, k: f64, c: f64) -> InterfaceSystemCard {
    interface_with_gap(a, b, k, c, false)
}

fn interface_with_gap(
    a: &MaterialCard,
    b: &MaterialCard,
    k: f64,
    c: f64,
    limited: bool,
) -> InterfaceSystemCard {
    let mut claims = ClaimSet::new();
    for (axis, (stiffness, viscosity)) in PROPERTIES.into_iter().enumerate() {
        for (name, dims, value, thermal_scale) in [
            (
                stiffness,
                interface_stiffness_per_area_dims(),
                [100., 200., 300.][axis] * k,
                2.0,
            ),
            (
                viscosity,
                interface_viscosity_per_area_dims(),
                [5., 7., 11.][axis] * c,
                3.0,
            ),
        ] {
            let mut validity = ValidityDomain::unconstrained()
                .with("T", 300., 400.)
                .with("omega", 0., 1e4);
            if limited && name == "cz" {
                validity = validity.with("gap", -0.01, 0.01);
            }
            claims
                .insert_claim(PropertyClaim {
                    key: PropertyKey::new(name, dims),
                    value: PropertyValue::Curve {
                        abscissa: "T".into(),
                        abscissa_dims: Temperature::DIMS,
                        knots: vec![(300.0, value), (400.0, value * thermal_scale)],
                        dims,
                    },
                    validity,
                    uncertainty: UncertaintyModel::Unstated,
                    interpolation: InterpolationPolicy::LinearInside,
                    observations: vec![],
                    provenance: Provenance {
                        source: "synthetic bonded layer, not measured material data".into(),
                        license: "CC0-1.0".into(),
                        artifact: None,
                    },
                })
                .unwrap();
        }
    }
    InterfaceSystemCard::assemble(
        SurfaceSpec {
            material: a.id().clone(),
            texture_frame: "pad axes".into(),
        },
        SurfaceSpec {
            material: b.id().clone(),
            texture_frame: "pad axes".into(),
        },
        SystemContext {
            medium: "bonded fixture".into(),
            third_body: None,
            environment: "declared external temperature".into(),
            history: "undegraded".into(),
        },
        claims,
        vec![],
    )
    .unwrap()
}

fn source(card: &InterfaceSystemCard) -> ResolvedInterfaceStatePoint {
    source_at(card, &point(300.).with("omega", 0.).unwrap())
}

fn source_at(card: &InterfaceSystemCard, query: &QueryPoint) -> ResolvedInterfaceStatePoint {
    let requirements: Vec<_> = PROPERTIES
        .into_iter()
        .flat_map(|(k, c)| {
            [
                ScalarPropertyRequirement::try_new(
                    k,
                    interface_stiffness_per_area_dims(),
                    ScalarAdmissibility::NonNegative,
                )
                .unwrap(),
                ScalarPropertyRequirement::try_new(
                    c,
                    interface_viscosity_per_area_dims(),
                    ScalarAdmissibility::NonNegative,
                )
                .unwrap(),
            ]
        })
        .collect();
    resolve_interface_state_point(
        card,
        query,
        &requirements,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap()
}

fn geometry(rotation: Vec3) -> RectangularPad {
    RectangularPad {
        frames: [0, 1].map(|link| MaterialFrame {
            link,
            link_from_attachment: Se3::from_parts(
                So3::exp(So3Tangent::new(rotation)).unwrap(),
                ZERO,
            )
            .unwrap(),
        }),
        rest_translation_a: ZERO,
        size_m: [0.2, 0.3],
        maximum_displacement_m: 0.1,
        maximum_rotation_rad: 0.05,
    }
}

fn binding<'a>(
    mechanism: &'a ResolvedArticulated,
    card: &InterfaceSystemCard,
) -> MaterialBushing<'a> {
    MaterialBushing::rectangular_pad(
        mechanism,
        ["surface"; 2],
        geometry(ZERO),
        source(card),
        PROPERTIES,
        || false,
    )
    .unwrap()
}

fn state(q: f64, v: f64) -> FreeFloatingMotionState {
    FreeFloatingMotionState {
        base: FreeFloatingBaseState::stationary(Se3::identity()),
        joint_position: vec![q],
        joint_velocity: vec![v],
    }
}

#[test]
fn rectangular_source_integration_matches_surface_quadrature_and_actual_frame() {
    let a = density_card("a", 2., 2.);
    let b = density_card("b", 3., 3.);
    let card = interface(&a, &b, 1., 1.);
    for rotation in [false, true] {
        for axis in [
            Vec3::new(1., 0., 0.),
            Vec3::new(0., 1., 0.),
            Vec3::new(0., 0., 1.),
        ] {
            let joint = if rotation {
                JointModel::revolute(axis, None)
            } else {
                JointModel::prismatic(axis, None)
            }
            .unwrap();
            let model = mechanism(&a, &b, joint);
            let pad = binding(&model, &card);
            let s = state(0.01, 0.04);
            let response = pad.evaluate(&s, || false).unwrap();
            // Independent 2x2 Gauss integration of the original area traction
            // potential and rates, not the derived rotational coefficients.
            let mut energy = 0.0;
            let mut loss = 0.0;
            let mut power = 0.0;
            for sx in [-1.0, 1.0] {
                for sy in [-1.0, 1.0] {
                    let r = Vec3::new(sx * 0.1 / 3.0_f64.sqrt(), sy * 0.15 / 3.0_f64.sqrt(), 0.);
                    let direction = if rotation { axis.cross(r) } else { axis };
                    let u = direction.scale(s.joint_position[0]);
                    let v = direction.scale(s.joint_velocity[0]);
                    let elastic = Vec3::new(100. * u.x, 200. * u.y, 300. * u.z);
                    let viscous = Vec3::new(5. * v.x, 7. * v.y, 11. * v.z);
                    energy += 0.015 * 0.5 * u.dot(elastic);
                    loss += 0.015 * v.dot(viscous);
                    power -= 0.015 * v.dot(elastic + viscous);
                }
            }
            close(response.stored_energy, energy);
            close(response.dissipation_power, loss);
            close(response.internal_power, power);
            let fk = forward_kinematics(
                model.model(),
                BaseState::stationary(Se3::identity()),
                &s.joint_position,
                &s.joint_velocity,
            )
            .unwrap();
            let mut force = ZERO;
            let mut torque = ZERO;
            for (pose, wrench) in fk.world_from_link.iter().zip(&response.body_wrench) {
                let f = pose.rotation().rotate(wrench.force).unwrap();
                force = force + f;
                torque = torque
                    + pose.rotation().rotate(wrench.torque).unwrap()
                    + pose.translation().cross(f);
            }
            assert!(force.norm() < 1e-12 && torque.norm() < 1e-12);
        }
    }
    let model = mechanism(
        &a,
        &b,
        JointModel::prismatic(Vec3::new(1., 0., 0.), None).unwrap(),
    );
    let turned = MaterialBushing::rectangular_pad(
        &model,
        ["surface"; 2],
        geometry(Vec3::new(0., 0., core::f64::consts::FRAC_PI_2)),
        source(&card),
        PROPERTIES,
        || false,
    )
    .unwrap();
    let mut s = state(0.01, 0.04);
    let response = turned.evaluate(&s, || false).unwrap();
    close(response.stored_energy, 0.5 * 12.0 * 0.01_f64.powi(2)); // local y, not x
    s.base.world_from_base = Se3::from_parts(
        So3::exp(So3Tangent::new(Vec3::new(0.4, 0.3, -0.2))).unwrap(),
        Vec3::new(30., -20., 40.),
    )
    .unwrap();
    s.base.twist_body = Twist::new(Vec3::new(0.3, 0.1, -0.4), Vec3::new(10., 20., 30.));
    assert_eq!(turned.evaluate(&s, || false).unwrap(), response);
}

#[test]
fn source_only_stiffness_and_loss_changes_drive_the_native_implicit_trajectory() {
    let a = density_card("a", 2., 2.);
    let b = density_card("b", 3., 3.);
    let model = mechanism(
        &a,
        &b,
        JointModel::prismatic(Vec3::new(1., 0., 0.), None).unwrap(),
    );
    let initial = state(0.01, 0.04);
    let dt = 0.02;
    let zero_wrench = Wrench::new(ZERO, ZERO);
    let loads = SpatialMotionLoads {
        held: FreeFloatingMotionLoads {
            generalized_force: &[0.],
            external_body_wrench: &[zero_wrench; 2],
            gravity_world: ZERO,
        },
        world_points: &[],
        connections: &[],
    };
    let mut endpoints = Vec::new();
    for (k_scale, c_scale) in [(1., 1.), (4., 1.), (1., 9.)] {
        let card = interface(&a, &b, k_scale, c_scale);
        let pad = binding(&model, &card);
        let accepted =
            step_bushings_implicit(&[&pad], &initial, loads, dt, Default::default(), || false)
                .unwrap();
        let step = &accepted.spatial;
        let actual = &step.motion.state_after;
        pad.for_state(actual, || false).unwrap();
        // Independent reduced-mass oscillator and its 2x2 midpoint equations.
        let mass = 2.0 * 3.0 / 5.0;
        let k = 6.0 * k_scale;
        let c = 0.3 * c_scale;
        let h = dt / (2.0 * mass);
        let g = dt * dt * k / (4.0 * mass);
        let v = ((1.0 - h * c - g) * 0.04 - dt * k * 0.01 / mass) / (1.0 + h * c + g);
        let q = 0.01 + 0.5 * dt * (0.04 + v);
        assert!((accepted.midpoint_state.joint_position[0] - 0.5 * (0.01 + q)).abs() < 2e-11);
        assert!((actual.joint_velocity[0] - v).abs() < 2e-9);
        assert!((actual.joint_position[0] - q).abs() < 2e-11);
        assert!((step.dissipated_energy - c * dt * ((v + 0.04) / 2.0).powi(2)).abs() < 2e-11);
        assert!(step.mechanical_work_residual.abs() < 2e-10);
        assert!(
            (step.motion.momentum_after.linear_momentum_world
                - step.motion.momentum_before.linear_momentum_world)
                .norm()
                < 2e-11
        );
        endpoints.push(actual.joint_velocity[0]);
    }
    assert!((endpoints[1] - endpoints[0]).abs() > 1e-3);
    assert!((endpoints[2] - endpoints[0]).abs() > 1e-3);
    let card = interface(&a, &b, 1., 1.);
    let mut bounded = geometry(ZERO);
    bounded.maximum_displacement_m = 0.0106; // Initial/midpoint pass; endpoint exits.
    let pad = MaterialBushing::rectangular_pad(
        &model,
        ["surface"; 2],
        bounded,
        source(&card),
        PROPERTIES,
        || false,
    )
    .unwrap();
    let raw = step_free_floating_bushings_implicit(
        model.model(),
        &initial,
        loads,
        &[pad.for_state(&initial, || false).unwrap()],
        dt,
        Default::default(),
        || false,
    )
    .unwrap();
    assert!(pad.for_state(&raw.midpoint_state, || false).is_ok());
    assert!(
        pad.for_state(&raw.spatial.motion.state_after, || false)
            .is_err()
    );
    assert!(
        step_bushings_implicit(&[&pad], &initial, loads, dt, Default::default(), || false).is_err()
    );
    bounded.maximum_displacement_m = 0.0102; // Accepted native midpoint now exits too.
    let narrow = MaterialBushing::rectangular_pad(
        &model,
        ["surface"; 2],
        bounded,
        source(&card),
        PROPERTIES,
        || false,
    )
    .unwrap();
    assert!(narrow.for_state(&raw.midpoint_state, || false).is_err());
    assert!(
        step_bushings_implicit(&[&narrow], &initial, loads, dt, Default::default(), || {
            false
        })
        .is_err()
    );
    assert_eq!(initial, state(0.01, 0.04));
}

#[test]
fn source_rebinding_accounts_for_actual_storage_and_refuses_partial_updates() {
    let a = density_card("a", 2., 2.);
    let b = density_card("b", 3., 3.);
    let model = mechanism(
        &a,
        &b,
        JointModel::prismatic(Vec3::new(1., 0., 0.), None).unwrap(),
    );
    let card = interface(&a, &b, 1., 1.);
    let pad = binding(&model, &card);
    let s = state(0.01, 0.04);
    let original = pad.evaluate(&s, || false).unwrap();
    let original_sources = pad.axes().clone();
    let hot = point(400.).with("omega", 0.).unwrap();
    let mut polls = 0;
    let rebound = pad
        .rebind_at_state(&card, &hot, &s, 1.0, || {
            polls += 1;
            false
        })
        .unwrap();
    close(
        rebound.parameter_work_j,
        0.5 * (12.0 - 6.0) * 0.01_f64.powi(2),
    );
    let updated = rebound.binding.evaluate(&s, || false).unwrap();
    close(
        updated.stored_energy,
        original.stored_energy + rebound.parameter_work_j,
    );
    close(updated.dissipation_power, 3.0 * original.dissipation_power);
    for axis in rebound.binding.axes() {
        assert_eq!(axis.resolved().state_coordinate("T"), Some(400.));
        assert_eq!(
            axis.resolved().card_identity(),
            pad.axes()[0].resolved().card_identity()
        );
    }
    let cold = point(300.).with("omega", 0.).unwrap();
    close(
        rebound
            .binding
            .rebind_at_state(&card, &cold, &s, 1., || false)
            .unwrap()
            .parameter_work_j,
        -rebound.parameter_work_j,
    );
    assert!(pad.rebind_at_state(&card, &hot, &s, 0., || false).is_err());
    assert!(
        pad.rebind_at_state(
            &card,
            &point(401.).with("omega", 0.).unwrap(),
            &s,
            1.,
            || false
        )
        .is_err()
    );
    let changed = interface(&a, &b, 2., 1.);
    assert!(
        pad.rebind_at_state(&changed, &hot, &s, 1., || false)
            .is_err()
    );
    let mut n = 0;
    assert!(matches!(
        pad.rebind_at_state(&card, &hot, &s, 1., || {
            n += 1;
            n == polls
        }),
        Err(AttachmentError::Cancelled)
    ));
    assert_eq!(pad.axes(), &original_sources);
    assert_eq!(pad.evaluate(&s, || false).unwrap(), original);
    assert!(pad.require_angular_frequency_rad_s(1e4).is_ok());
    assert!(pad.require_angular_frequency_rad_s(1e4 + 1.0).is_err());
    assert!(pad.for_state(&state(0.2, 0.), || false).is_err());
    assert!(
        MaterialBushing::rectangular_pad(
            &model,
            ["missing", "surface"],
            geometry(ZERO),
            source(&card),
            PROPERTIES,
            || false
        )
        .is_err()
    );
    let reversed = interface(&b, &a, 1., 1.);
    assert!(
        MaterialBushing::rectangular_pad(
            &model,
            ["surface"; 2],
            geometry(ZERO),
            source(&reversed),
            PROPERTIES,
            || false
        )
        .is_err()
    );
    let mut invalid = geometry(ZERO);
    invalid.size_m[0] = 0.;
    assert!(
        MaterialBushing::rectangular_pad(
            &model,
            ["surface"; 2],
            invalid,
            source(&card),
            PROPERTIES,
            || false
        )
        .is_err()
    );
    let rotor = mechanism(
        &a,
        &b,
        JointModel::revolute(Vec3::new(0., 0., 1.), None).unwrap(),
    );
    assert!(
        binding(&rotor, &card)
            .for_state(&state(0.06, 0.), || false)
            .is_err()
    );
    let restricted = interface_with_gap(&a, &b, 1., 1., true);
    let query = point(300.)
        .with("omega", 0.)
        .unwrap()
        .with("gap", 0.)
        .unwrap();
    assert!(
        MaterialBushing::rectangular_pad(
            &model,
            ["surface"; 2],
            geometry(ZERO),
            source_at(&restricted, &query),
            PROPERTIES,
            || false
        )
        .is_err()
    );
}
