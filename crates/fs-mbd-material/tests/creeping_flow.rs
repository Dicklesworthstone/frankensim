//! G1/G3/G4: source swaps -> original inertia -> actual fluid-loaded motion.
//! Synthetic complete cards test the numerical model, not measured materials.
mod common;
use common::{close, density_card, point};
use fs_evidence::ValidityDomain;
use fs_ga::{Se3, Twist, Vec3};
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, UncertaintyModel,
};
use fs_material::fluid_transport::{FluidTransport, FluidTransportState};
use fs_material::liquid::{DYNAMIC_VISCOSITY_DIMS, ResolvedLiquidState, resolve_liquid_state};
use fs_material::state_point::{
    MaterialPropertySelection, SPECIFIC_HEAT_CAPACITY_DIMS, THERMAL_CONDUCTIVITY_DIMS,
};
use fs_mbd::{Gravity, RigidBodyIntegrator};
use fs_mbd_material::creeping_flow::{
    CreepingFlowError, CreepingFlowLimits, CreepingSphere, CreepingStepControls,
};
use fs_mbd_material::{RegionBinding, ResolvedBody, VolumeMoments, resolve_body};
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_qty::{Density, QuantitySpec};

const RADIUS: f64 = 0.01;
fn limits() -> CreepingFlowLimits {
    CreepingFlowLimits {
        maximum_reynolds: 0.05,
        maximum_diffusion_ratio: 0.05,
    }
}
fn controls(dt: f64) -> CreepingStepControls {
    CreepingStepControls {
        duration_s: dt,
        maximum_relaxation_fraction: 0.5,
        maximum_energy_defect_j: 1e-12,
    }
}
fn body(density: f64, geometry: VolumeMoments) -> ResolvedBody {
    let card = density_card("synthetic sphere", density, density);
    resolve_body(
        &[RegionBinding {
            key: "sphere",
            geometry,
            material: &card,
            point: &point(300.0),
            selection: MaterialPropertySelection::SingleClaimOnly,
        }],
        1,
        || false,
    )
    .unwrap()
}
fn sphere(density: f64) -> ResolvedBody {
    body(density, VolumeMoments::sphere(RADIUS).unwrap())
}
fn fluid_card(density: f64, viscosity: f64) -> MaterialCard {
    let mut claims = ClaimSet::new();
    for (name, dims, value) in [
        ("density", Density::DIMS, density),
        ("dynamic_viscosity", DYNAMIC_VISCOSITY_DIMS, viscosity),
        (
            "specific_heat_capacity",
            SPECIFIC_HEAT_CAPACITY_DIMS,
            1000.0,
        ),
        ("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS, 0.1),
    ] {
        claims
            .insert_claim(PropertyClaim {
                key: if name == "thermal_conductivity" {
                    PropertyKey::with_quantity(
                        name,
                        QuantitySpec::semantic(SemanticType::new(
                            QuantityKind::ThermalConductivity,
                            ValueForm::Static,
                        )),
                    )
                } else {
                    PropertyKey::new(name, dims)
                },
                value: PropertyValue::Scalar { value, dims },
                validity: ValidityDomain::unconstrained().with("T", 200.0, 400.0),
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::ConstantWithinValidity,
                observations: vec![],
                provenance: Provenance {
                    source: "synthetic creeping fluid".into(),
                    license: "CC0-1.0".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic fluid".into(),
            phase: "liquid".into(),
            process: "implementation fixture".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap()
}
fn fluid(density: f64, viscosity: f64) -> ResolvedLiquidState {
    resolve_liquid_state(
        &fluid_card(density, viscosity),
        &point(300.0),
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap()
}
fn gravity() -> Gravity {
    Gravity::new(fs_mbd::Vec3::new(0.0, 0.0, -0.01)).unwrap()
}

#[test]
fn sphere_geometry_and_source_density_supply_the_actual_native_inertia() {
    let density = 20000.0;
    let source = sphere(density);
    let volume = 4.0 * core::f64::consts::PI * RADIUS.powi(3) / 3.0;
    let inertia = source.inertia();
    close(inertia.mass(), density * volume);
    for (i, value) in inertia.inertia_com().m.into_iter().enumerate() {
        close(
            value,
            if i % 4 == 0 {
                2.0 * density * volume * RADIUS.powi(2) / 5.0
            } else {
                0.0
            },
        );
    }
    for radius in [0.0, -1.0, f64::NAN, f64::INFINITY, 1e200, 1e-200] {
        assert!(VolumeMoments::sphere(radius).is_err());
    }
}

#[test]
fn data_only_body_and_fluid_swaps_change_the_full_settling_trajectory() {
    for (density, fluid_density, viscosity) in [
        (10000.0, 1.0, 1.0),
        (20000.0, 1.0, 1.0),
        (20000.0, 1.0, 2.0),
        (20000.0, 2.0, 1.0),
    ] {
        let source = sphere(density);
        let body = source.principal_body(1e-12, || false).unwrap();
        let material = fluid(fluid_density, viscosity);
        let medium =
            CreepingSphere::new(&material, RADIUS, gravity(), fs_mbd::Vec3::ZERO, limits())
                .unwrap();
        assert!(core::ptr::eq(medium.source(), &material));
        let mut state = body
            .state_from_kinematics(Se3::identity(), Twist::zero())
            .unwrap();
        let mut heat = 0.0;
        let mut change = 0.0;
        let dt = 0.005;
        for _ in 0..200 {
            let result = medium.step(&body, state, controls(dt), || false).unwrap();
            let native = RigidBodyIntegrator::new(gravity())
                .step(state, body.properties(), result.load_midpoint.wrench, dt)
                .unwrap();
            assert_eq!(result.dynamics, native);
            assert_eq!(
                result.load_after,
                medium.load(&body, native.state_after, || false).unwrap()
            );
            heat += result.viscous_heat_j;
            change += result.stored_energy_change_j;
            state = result.dynamics.state_after;
        }
        // Independent continuum settling equation m*dv/dt=(rho_s-rho_f)*V*g-6*pi*mu*a*v.
        let tau = 2.0 * density * RADIUS.powi(2) / (9.0 * viscosity);
        let terminal = -2.0 * (density - fluid_density) * 0.01 * RADIUS.powi(2) / (9.0 * viscosity);
        let velocity = terminal * (1.0 - (-1.0 / tau).exp());
        let position = terminal * (1.0 - tau * (1.0 - (-1.0 / tau).exp()));
        let actual_velocity = state
            .center_of_mass_velocity_world(body.properties())
            .unwrap()
            .z;
        assert!((actual_velocity / velocity - 1.0).abs() < 2e-5);
        assert!((state.pose().position_world().z / position - 1.0).abs() < 2e-5);
        assert!(heat > 0.0);
        assert!((change + heat).abs() < 1e-15);
    }
}

fn rotating_trajectory(dt: f64, steps: usize) -> (f64, f64) {
    let source = sphere(20000.0);
    let body = source.principal_body(1e-12, || false).unwrap();
    let material = fluid(1.0, 1.0);
    let medium = CreepingSphere::new(
        &material,
        RADIUS,
        Gravity::ZERO,
        fs_mbd::Vec3::ZERO,
        limits(),
    )
    .unwrap();
    let mut state = body
        .state_from_kinematics(
            Se3::identity(),
            Twist::new(Vec3::new(0.0, 0.0, 1.0), Vec3::new(0.01, 0.0, 0.0)),
        )
        .unwrap();
    let energy_before = state.kinetic_energy(body.properties()).unwrap();
    let mut heat = 0.0;
    for _ in 0..steps {
        let result = medium.step(&body, state, controls(dt), || false).unwrap();
        heat += result.viscous_heat_j;
        state = result.dynamics.state_after;
    }
    close(
        state.kinetic_energy(body.properties()).unwrap() + heat,
        energy_before,
    );
    let velocity = state
        .center_of_mass_velocity_world(body.properties())
        .unwrap()
        .x;
    let omega = body
        .properties()
        .angular_velocity_body_checked(state.angular_momentum_body())
        .unwrap();
    let omega_world = state.pose().orientation().rotate_body_to_world(omega);
    (velocity, omega_world.z)
}

#[test]
fn native_translation_and_rotation_converge_to_independent_exponential_limits() {
    let duration: f64 = 0.4;
    let translation_exact = 0.01 * (-duration / (2.0 * 20000.0 * RADIUS.powi(2) / 9.0)).exp();
    let rotation_exact = (-duration / (20000.0 * RADIUS.powi(2) / 15.0)).exp();
    let mut errors = Vec::new();
    for (dt, steps) in [(0.02, 20), (0.01, 40), (0.005, 80)] {
        let (v, omega) = rotating_trajectory(dt, steps);
        errors.push([
            (v - translation_exact).abs(),
            (omega - rotation_exact).abs(),
        ]);
    }
    for pair in errors.windows(2) {
        for i in 0..2 {
            assert!(
                (3.9..4.1).contains(&(pair[0][i] / pair[1][i])),
                "{errors:?}"
            );
        }
    }
}

#[test]
fn maintained_flow_supplies_work_and_fluid_density_changes_only_the_correct_terms() {
    let source = sphere(20000.0);
    let body = source.principal_body(1e-12, || false).unwrap();
    let state = body
        .state_from_kinematics(Se3::identity(), Twist::zero())
        .unwrap();
    let a = fluid(1.0, 1.0);
    let b = fluid(2.0, 1.0);
    let flow = fs_mbd::Vec3::new(0.01, 0.0, 0.0);
    let medium = CreepingSphere::new(&a, RADIUS, Gravity::ZERO, flow, limits()).unwrap();
    let step = medium.step(&body, state, controls(0.01), || false).unwrap();
    assert!(step.dynamics.state_after.linear_momentum_world().x > 0.0);
    assert!(step.viscous_heat_j > 0.0);
    assert!(step.background_work_j > step.viscous_heat_j);
    close(
        step.stored_energy_change_j + step.viscous_heat_j,
        step.background_work_j,
    );
    let load = |material: &ResolvedLiquidState| {
        CreepingSphere::new(material, RADIUS, gravity(), flow, limits())
            .unwrap()
            .load(&body, state, || false)
            .unwrap()
    };
    let first = load(&a);
    let second = load(&b);
    assert_eq!(first.drag_force_world_n, second.drag_force_world_n);
    close(
        second.buoyancy_force_world_n.z,
        2.0 * first.buoyancy_force_world_n.z,
    );
    close(
        second.translation_reynolds,
        2.0 * first.translation_reynolds,
    );
}

#[test]
fn wrong_geometry_and_excluded_flow_regimes_refuse_before_publication() {
    let source = sphere(20000.0);
    let body = source.principal_body(1e-12, || false).unwrap();
    let state = body
        .state_from_kinematics(Se3::identity(), Twist::zero())
        .unwrap();
    let material = fluid(1.0, 1.0);
    let medium =
        CreepingSphere::new(&material, RADIUS, gravity(), fs_mbd::Vec3::ZERO, limits()).unwrap();
    let cuboid = equal_volume_cube();
    let wrong = cuboid.principal_body(1e-12, || false).unwrap();
    assert!(matches!(
        medium.load(&wrong, state, || false),
        Err(CreepingFlowError::Invalid(_))
    ));
    let dense_fluid = fluid(1000.0, 1.0);
    let unsteady = CreepingSphere::new(
        &dense_fluid,
        RADIUS,
        gravity(),
        fs_mbd::Vec3::ZERO,
        limits(),
    )
    .unwrap();
    assert!(matches!(
        unsteady.step(&body, state, controls(0.01), || false),
        Err(CreepingFlowError::Unsteady { .. })
    ));
    let fast = CreepingSphere::new(
        &material,
        RADIUS,
        Gravity::new(fs_mbd::Vec3::new(0.0, 0.0, -1000.0)).unwrap(),
        fs_mbd::Vec3::ZERO,
        limits(),
    )
    .unwrap();
    assert!(matches!(
        fast.step(&body, state, controls(0.01), || false),
        Err(CreepingFlowError::Reynolds { .. })
    ));
    assert!(medium.step(&body, state, controls(1.0), || false).is_err());
    for bad in [controls(0.0), controls(f64::NAN), controls(f64::INFINITY)] {
        assert!(medium.step(&body, state, bad, || false).is_err());
    }
    let mut bad = controls(0.01);
    bad.maximum_relaxation_fraction = 1.01;
    assert!(medium.step(&body, state, bad, || false).is_err());
    for maximum in [0.0, 0.11, f64::NAN, f64::INFINITY] {
        assert!(
            CreepingSphere::new(
                &material,
                RADIUS,
                gravity(),
                fs_mbd::Vec3::ZERO,
                CreepingFlowLimits {
                    maximum_reynolds: maximum,
                    ..limits()
                }
            )
            .is_err()
        );
    }
    assert_eq!(state.linear_momentum_world(), fs_mbd::Vec3::ZERO);
}
fn equal_volume_cube() -> ResolvedBody {
    // Equal volume does not make a cube an admissible spherical boundary.
    let side = (4.0 * core::f64::consts::PI / 3.0).cbrt() * RADIUS;
    body(
        20000.0,
        VolumeMoments::cuboid(Vec3::new(side, side, side)).unwrap(),
    )
}

#[test]
fn source_is_resolved_once_and_cancellation_preserves_the_entire_native_state() {
    struct Counted(std::cell::Cell<usize>);
    impl FluidTransport for Counted {
        fn transport_state(&self) -> Result<FluidTransportState, fs_material::MaterialError> {
            self.0.set(self.0.get() + 1);
            FluidTransportState::try_new(1.0, 1000.0, 0.1, 1.0)
        }
    }
    let fluid = Counted(std::cell::Cell::new(0));
    let medium =
        CreepingSphere::new(&fluid, RADIUS, gravity(), fs_mbd::Vec3::ZERO, limits()).unwrap();
    let source = sphere(20000.0);
    let body = source.principal_body(1e-12, || false).unwrap();
    let state = body
        .state_from_kinematics(Se3::identity(), Twist::zero())
        .unwrap();
    let saved = state;
    let mut calls = 0;
    let expected = medium
        .step(&body, state, controls(0.01), || {
            calls += 1;
            false
        })
        .unwrap();
    for stop in 1..=calls {
        let mut current = 0;
        assert!(matches!(
            medium.step(&body, state, controls(0.01), || {
                current += 1;
                current == stop
            }),
            Err(CreepingFlowError::Cancelled)
        ));
        assert_eq!(state, saved);
    }
    assert_eq!(
        medium.step(&body, state, controls(0.01), || false).unwrap(),
        expected
    );
    assert_eq!(fluid.0.get(), 1);
}
