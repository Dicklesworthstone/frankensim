//! Thermodynamic phase count and sourced composition reach actual buoyancy.

use super::*;
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::liquid::solution::{BinaryLiquidSolutionLaw, ResolvedBinaryLiquidEquilibrium};
use fs_material::liquid::{
    DYNAMIC_VISCOSITY_DIMS, DYNAMIC_VISCOSITY_PROPERTY, resolve_liquid_state,
};
use fs_material::state_point::{
    DENSITY_PROPERTY, SPECIFIC_HEAT_CAPACITY_DIMS, SPECIFIC_HEAT_CAPACITY_PROPERTY,
    THERMAL_CONDUCTIVITY_DIMS, THERMAL_CONDUCTIVITY_PROPERTY,
};
use fs_mbd_material::hydrostatic::{
    FiniteLayeredHydrostaticTank, FiniteTankControls, HydrostaticStepControls,
    HydrostaticTankGeometry,
};
use fs_qty::semantic::{CompositionBasis, QuantityKind, SemanticType, ValueForm};
use fs_qty::{Density, Dims, MolarMass, QuantitySpec, Temperature};
use fs_thermochem::SpeciesId;
use fs_thermochem::solution::{BinaryCoexistenceOptionsV1, BinaryRegularSolutionV1};

fn quantity(kind: QuantityKind) -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(kind, ValueForm::Static))
}

fn liquid_card() -> MaterialCard {
    // Explicit synthetic composition curves, not a pure-fluid mixing rule or
    // measured properties of named liquids. Density rises with species A, so
    // thermochemical composition order is opposite to hydrostatic layer order.
    let mut claims = ClaimSet::new();
    for (name, dims, values) in [
        (DENSITY_PROPERTY, Density::DIMS, [800.0, 1200.0]),
        (
            SPECIFIC_HEAT_CAPACITY_PROPERTY,
            SPECIFIC_HEAT_CAPACITY_DIMS,
            [2000.0, 3000.0],
        ),
        (
            THERMAL_CONDUCTIVITY_PROPERTY,
            THERMAL_CONDUCTIVITY_DIMS,
            [0.2, 0.6],
        ),
        (
            DYNAMIC_VISCOSITY_PROPERTY,
            DYNAMIC_VISCOSITY_DIMS,
            [0.002, 0.001],
        ),
    ] {
        let value_quantity = if name == THERMAL_CONDUCTIVITY_PROPERTY {
            quantity(QuantityKind::ThermalConductivity)
        } else {
            QuantitySpec::dimensional(dims)
        };
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::with_quantity(name, value_quantity),
                value: PropertyValue::Curve {
                    abscissa: "x_a".into(),
                    abscissa_dims: Dims([0; 6]),
                    knots: vec![(0.0, values[0]), (1.0, values[1])],
                    dims,
                },
                validity: ValidityDomain::unconstrained()
                    .with_quantity(
                        "T",
                        quantity(QuantityKind::AbsoluteTemperature),
                        300.0,
                        600.0,
                    )
                    .with_quantity(
                        "x_a",
                        quantity(QuantityKind::Composition(CompositionBasis::MoleFraction)),
                        0.0,
                        1.0,
                    )
                    .with("pressure", 100_000.0, 100_000.0),
                interpolation: InterpolationPolicy::LinearInside,
                uncertainty: UncertaintyModel::Unstated,
                observations: vec![],
                provenance: Provenance {
                    source:
                        "synthetic binary-liquid hydrostatics fixture; not measured material data"
                            .into(),
                    license: "CC0-1.0".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic A/B liquid".into(),
            phase: "liquid".into(),
            process: "numerical fixture".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap()
}

fn equilibrium(card: &MaterialCard, temperature: f64) -> ResolvedBinaryLiquidEquilibrium {
    let point = QueryPoint::new()
        .with_quantity("T", quantity(QuantityKind::AbsoluteTemperature), 300.0)
        .unwrap()
        .with_quantity(
            "x_a",
            quantity(QuantityKind::Composition(CompositionBasis::MoleFraction)),
            0.5,
        )
        .unwrap()
        .with("pressure", 100_000.0)
        .unwrap();
    let reference =
        resolve_liquid_state(card, &point, MaterialPropertySelection::SingleClaimOnly).unwrap();
    let model = BinaryRegularSolutionV1::new(
        [SpeciesId::new("A").unwrap(), SpeciesId::new("B").unwrap()],
        0.0,
        900.0,
        Temperature::new(300.0),
        Temperature::new(600.0),
    )
    .unwrap();
    BinaryLiquidSolutionLaw::try_new(
        model,
        [MolarMass::new(0.018), MolarMass::new(0.046)],
        card,
        reference,
        "T",
        "x_a",
    )
    .unwrap()
    .equilibrate(
        Temperature::new(temperature),
        [500_000.0, 500_000.0],
        BinaryCoexistenceOptionsV1::default(),
        || false,
    )
    .unwrap()
}

fn geometry() -> HydrostaticTankGeometry {
    HydrostaticTankGeometry {
        world_from_tank: Se3::identity(),
        inside_dimensions_m: Vec3::new(4.0, 4.0, 6.0),
    }
}

fn controls() -> FiniteTankControls {
    FiniteTankControls {
        maximum_volume_defect_m3: 1e-12,
        maximum_volume_evaluations: 64,
    }
}

#[test]
fn miscibility_changes_inventory_interfaces_and_native_buoyancy() {
    let card = liquid_card();
    let cold = equilibrium(&card, 300.0);
    let hot = equilibrium(&card, 600.0);
    assert_eq!(cold.phases().len(), 2);
    assert_eq!(hot.phases().len(), 1);
    assert_eq!(cold.total_component_moles(), hot.total_component_moles());
    assert_eq!(cold.total_mass_kg(), hot.total_mass_kg());
    close(cold.total_mass_kg(), 32_000.0);
    for result in [&cold, &hot] {
        for component in 0..2 {
            close(
                result
                    .phases()
                    .iter()
                    .map(|phase| phase.component_moles()[component])
                    .sum(),
                500_000.0,
            );
        }
    }
    #[cfg(feature = "electrothermal-coupling")]
    for result in [&cold, &hot] {
        for phase in result.phases() {
            let scaling =
                fs_lbm::material::plan_scaling_from_material(phase, 0.001, 0.001, 0.01).unwrap();
            assert_eq!(*scaling.transport(), phase.transport_state().unwrap());
            assert!(std::ptr::eq(scaling.material(), phase));
        }
    }
    let cold_tank = FiniteLayeredHydrostaticTank::from_binary_equilibrium(
        &cold,
        gravity(),
        geometry(),
        controls(),
        || false,
    )
    .unwrap();
    let hot_tank = FiniteLayeredHydrostaticTank::from_binary_equilibrium(
        &hot,
        gravity(),
        geometry(),
        controls(),
        || false,
    )
    .unwrap();
    let bottom = &cold.phases()[1];
    let top = &cold.phases()[0];
    assert!(bottom.liquid().density_kg_m3() > top.liquid().density_kg_m3());
    let retained: Vec<_> = cold_tank.layers().collect();
    assert!(std::ptr::addr_eq(
        retained[0].0.source,
        bottom as &dyn FluidTransport
    ));
    assert!(std::ptr::addr_eq(
        retained[1].0.source,
        top as &dyn FluidTransport
    ));
    assert_eq!(
        cold_tank.layer_masses_kg().collect::<Vec<_>>(),
        vec![bottom.mass_kg(), top.mass_kg()]
    );
    assert_eq!(cold_tank.liquid_mass_kg(), cold.total_mass_kg());
    close(cold_tank.liquid_volume_m3(), cold.total_volume_m3());
    // Equal phase amounts at z_A=1/2 do not imply equal phase masses/volumes.
    assert!((bottom.volume_m3() / cold.total_volume_m3() - 0.5).abs() > 0.1);

    let source = body(1000.0, Vec3::new(2.0, 2.0, 2.0));
    let body = source.principal_body(1e-12, || false).unwrap();
    let mesh = cube();
    let held = pose(1.25);
    let state = body.state_from_kinematics(held, Twist::zero()).unwrap();
    let interfaces = cold_tank.interfaces_at_pose(&mesh, held, || false).unwrap();
    // Cube extends from 0.25 to 2.25 m. Lower interface cuts its 4 m2 footprint;
    // the free surface lies above the complete 8 m3 displacement envelope.
    let lower = (bottom.volume_m3() - 1.0) / 12.0;
    assert!((0.25..2.25).contains(&lower));
    close(interfaces[0].surface_height_m, lower);
    close(
        interfaces[1].surface_height_m,
        (cold.total_volume_m3() + 8.0) / 16.0,
    );
    let lower_displacement = 4.0 * (lower - 0.25);
    let expected_mass = bottom.liquid().density_kg_m3() * lower_displacement
        + top.liquid().density_kg_m3() * (8.0 - lower_displacement);
    let cold_load = cold_tank.load(&mesh, &body, state, || false).unwrap();
    let hot_load = hot_tank.load(&mesh, &body, state, || false).unwrap();
    close(cold_load.wrench().force_world.z, 10.0 * expected_mass);
    close(hot_load.wrench().force_world.z, 80_000.0);
    assert!((cold_load.wrench().force_world.z - hot_load.wrench().force_world.z).abs() > 1000.0);
    assert_eq!(
        hot_tank
            .interfaces_at_pose(&mesh, held, || false)
            .unwrap()
            .len(),
        1
    );

    let step_controls = HydrostaticStepControls {
        duration_s: 1e-3,
        maximum_energy_defect_j: 1e-5,
    };
    let cold_step = cold_tank
        .step(&mesh, &body, state, step_controls, || false)
        .unwrap();
    let hot_step = hot_tank
        .step(&mesh, &body, state, step_controls, || false)
        .unwrap();
    assert!(
        cold_step
            .dynamics
            .state_after
            .linear_momentum_world()
            .z
            .abs()
            > 1.0
    );
    close(hot_step.dynamics.state_after.linear_momentum_world().z, 0.0);
    assert_eq!(cold_tank.liquid_mass_kg(), cold.total_mass_kg());
}

#[test]
fn phase_tank_construction_cancels_before_publishing_any_inventory() {
    let card = liquid_card();
    for temperature in [300.0, 600.0] {
        let result = equilibrium(&card, temperature);
        let before: Vec<_> = result
            .phases()
            .iter()
            .map(|phase| phase.component_moles())
            .collect();
        let mut polls = 0;
        let tank = FiniteLayeredHydrostaticTank::from_binary_equilibrium(
            &result,
            gravity(),
            geometry(),
            controls(),
            || {
                polls += 1;
                false
            },
        )
        .unwrap();
        for stop in 1..=polls {
            let mut count = 0;
            assert!(matches!(
                FiniteLayeredHydrostaticTank::from_binary_equilibrium(
                    &result,
                    gravity(),
                    geometry(),
                    controls(),
                    || {
                        count += 1;
                        count == stop
                    },
                ),
                Err(HydrostaticError::Cancelled)
            ));
        }
        assert_eq!(
            result
                .phases()
                .iter()
                .map(|phase| phase.component_moles())
                .collect::<Vec<_>>(),
            before
        );
        assert_eq!(tank.liquid_mass_kg(), result.total_mass_kg());
    }
}
