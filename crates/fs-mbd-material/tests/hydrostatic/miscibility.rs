//! Thermodynamic phase count and sourced composition reach actual buoyancy.

use super::*;
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSelection, ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim,
    PropertyKey, PropertyValue, Provenance, QueryPoint, SelectionPolicy, UncertaintyModel,
};
use fs_material::liquid::solution::caloric::BinaryLiquidCaloricLaw;
use fs_material::liquid::solution::{BinaryLiquidSolutionLaw, ResolvedBinaryLiquidEquilibrium};
use fs_material::liquid::{
    DYNAMIC_VISCOSITY_DIMS, DYNAMIC_VISCOSITY_PROPERTY, resolve_liquid_state,
};
use fs_material::state_point::{
    DENSITY_PROPERTY, SPECIFIC_HEAT_CAPACITY_DIMS, SPECIFIC_HEAT_CAPACITY_PROPERTY,
    THERMAL_CONDUCTIVITY_DIMS, THERMAL_CONDUCTIVITY_PROPERTY,
};
use fs_material::thermal::transport::isobaric::{
    IsobaricHeatCapacityLaw, IsobaricHeatCapacityQuery,
};
use fs_mbd_material::hydrostatic::{
    FiniteLayeredHydrostaticTank, FiniteTankControls, HydrostaticStepControls,
    HydrostaticTankGeometry,
};
use fs_qty::semantic::{CompositionBasis, QuantityKind, SemanticType, ValueForm};
use fs_qty::{Density, Dims, MolarMass, QuantitySpec, Temperature};
use fs_thermochem::solution::caloric::BinaryCaloricSolveV1;
use fs_thermochem::solution::{BinaryCoexistenceOptionsV1, BinaryRegularSolutionV1};
use fs_thermochem::{SpeciesId, UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R};

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

fn mixture_point(temperature: f64, composition: f64) -> QueryPoint {
    QueryPoint::new()
        .with_quantity(
            "T",
            quantity(QuantityKind::AbsoluteTemperature),
            temperature,
        )
        .unwrap()
        .with_quantity(
            "x_a",
            quantity(QuantityKind::Composition(CompositionBasis::MoleFraction)),
            composition,
        )
        .unwrap()
        .with_quantity("pressure", quantity(QuantityKind::Pressure), 100_000.0)
        .unwrap()
}

fn solution_law(card: &MaterialCard) -> BinaryLiquidSolutionLaw<'_> {
    let point = mixture_point(300.0, 0.5);
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
}

fn equilibrium(card: &MaterialCard, temperature: f64) -> ResolvedBinaryLiquidEquilibrium {
    solution_law(card)
        .equilibrate(
            Temperature::new(temperature),
            [500_000.0, 500_000.0],
            BinaryCoexistenceOptionsV1::default(),
            || false,
        )
        .unwrap()
}

fn pure_capacity(card: &MaterialCard, composition: f64) -> IsobaricHeatCapacityLaw {
    // The explicitly selected x_A=1 and x_A=0 endpoints are the synthetic
    // pure-liquid Cp sources; total-mixture Cp at x_A=1/2 is not substituted.
    let lower = mixture_point(300.0, composition);
    let upper = mixture_point(600.0, composition);
    let key = PropertyKey::new(SPECIFIC_HEAT_CAPACITY_PROPERTY, SPECIFIC_HEAT_CAPACITY_DIMS);
    IsobaricHeatCapacityLaw::try_new(
        card,
        IsobaricHeatCapacityQuery {
            specific_heat_cp: &key,
            temperature_axis: "T",
            pressure_axis: "pressure",
            lower: &lower,
            upper: &upper,
            selection: ClaimSelection::Policy(SelectionPolicy::SingleClaimOnly),
        },
        300.0,
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
fn sourced_heat_crosses_miscibility_gap_and_changes_native_buoyancy() {
    let card = liquid_card();
    let liquid = solution_law(&card);
    let caloric = BinaryLiquidCaloricLaw::try_new(
        &liquid,
        [
            (SpeciesId::new("A").unwrap(), pure_capacity(&card, 1.0)),
            (SpeciesId::new("B").unwrap(), pure_capacity(&card, 0.0)),
        ],
        "pressure",
    )
    .unwrap();
    let initial = caloric
        .equilibrate(
            Temperature::new(300.0),
            [500_000.0, 500_000.0],
            BinaryCoexistenceOptionsV1::default(),
            || false,
        )
        .unwrap();
    let cold = initial.equilibrium();
    assert_eq!(cold.phases().len(), 2);

    // Independent constant-Cp caloric oracle, not a prescribed final state.
    // At 500 K chi=1.8: one liquid. At 300 K the two phases have x and 1-x.
    // Both pure sensible increments and the change of mixing enthalpy are
    // necessary to determine the supplied heat and recover that temperature.
    let [x, y] = cold.phases()[0].mole_fractions();
    let sensible_heat = 500_000.0 * (0.018 * 3000.0 + 0.046 * 2000.0) * 200.0;
    let mixing_heat = 1_000_000.0 * R * 900.0 * (0.25 - x * y);
    let heat = sensible_heat + mixing_heat;
    assert!(mixing_heat > 0.05 * heat);
    let heated = caloric
        .add_heat(
            &initial,
            heat,
            BinaryCaloricSolveV1 {
                temperature_bracket: [Temperature::new(300.0), Temperature::new(600.0)],
                temperature_tolerance_k: 1e-8,
                energy_tolerance_j: 1e-3,
                maximum_iterations: 128,
                coexistence: BinaryCoexistenceOptionsV1::default(),
            },
            || false,
        )
        .unwrap();
    assert!((heated.equilibrium().temperature().value() - 500.0).abs() < 1e-8);
    assert!(heated.balance().energy_residual_j().abs() <= 1e-3);
    close(heated.balance().enthalpy_change_j(), heat);
    let hot = heated.equilibrium();
    assert_eq!(hot.phases().len(), 1);
    assert_eq!(hot.total_component_moles(), cold.total_component_moles());
    assert_eq!(hot.total_mass_kg(), cold.total_mass_kg());
    assert_eq!(hot.phases()[0].liquid().density_kg_m3(), 1000.0);

    // Each thermal endpoint supplies a frozen mechanical snapshot. This does
    // not claim a closed thermal/gravitational balance during interface motion.
    let source = body(1000.0, Vec3::new(2.0, 2.0, 2.0));
    let body = source.principal_body(1e-12, || false).unwrap();
    let mesh = cube();
    let held = pose(1.25);
    let state = body.state_from_kinematics(held, Twist::zero()).unwrap();
    let step_controls = HydrostaticStepControls {
        duration_s: 1e-3,
        maximum_energy_defect_j: 1e-5,
    };
    let mut momenta = [0.0; 2];
    for (index, equilibrium) in [cold, hot].into_iter().enumerate() {
        let tank = FiniteLayeredHydrostaticTank::from_binary_equilibrium(
            equilibrium,
            gravity(),
            geometry(),
            controls(),
            || false,
        )
        .unwrap();
        assert_eq!(tank.liquid_mass_kg(), 32_000.0);
        close(
            tank.layer_masses_kg().sum::<f64>(),
            equilibrium.total_mass_kg(),
        );
        let interfaces = tank.interfaces_at_pose(&mesh, held, || false).unwrap();
        assert_eq!(interfaces.len(), if index == 0 { 2 } else { 1 });
        close(
            interfaces.last().unwrap().surface_height_m,
            (equilibrium.total_volume_m3() + 8.0) / 16.0,
        );
        let expected_buoyancy = if index == 0 {
            let bottom = &cold.phases()[1];
            let top = &cold.phases()[0];
            let lower = (bottom.volume_m3() - 1.0) / 12.0;
            close(interfaces[0].surface_height_m, lower);
            let displaced_bottom = 4.0 * (lower - 0.25);
            10.0 * (bottom.liquid().density_kg_m3() * displaced_bottom
                + top.liquid().density_kg_m3() * (8.0 - displaced_bottom))
        } else {
            80_000.0
        };
        close(
            tank.load(&mesh, &body, state, || false)
                .unwrap()
                .wrench()
                .force_world
                .z,
            expected_buoyancy,
        );
        momenta[index] = tank
            .step(&mesh, &body, state, step_controls, || false)
            .unwrap()
            .dynamics
            .state_after
            .linear_momentum_world()
            .z;
    }
    assert!(momenta[0].abs() > 1.0);
    close(momenta[1], 0.0);
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
