//! Physical-time equilibrium propulsion: analytical references, not calibrated hydrogen.
use fs_material::thermochemical_gas::energy::parcel::transfer::nozzle::rocket::{
    RocketNozzle, equilibrium::{equilibrium_rocket_nozzle_flow, EquilibriumRocketNozzleSolve,
        EquilibriumRocketNozzleError, blowdown::{EquilibriumBlowdownSolve, EquilibriumGasReservoir,
        time::{EquilibriumRocketBlowdownError, EquilibriumRocketBlowdownSolve}}},
};
use fs_thermochem::{ChargeVector, ElementalMatrix, ElementId, ElementalReferenceIdV1,
    Nasa9StandardStateModelV1, Nasa9RegionV1, ReferenceEquationOfStateV1,
    StandardStateConventionV1, StandardStatePhaseV1, SpeciesId, required_nasa9_parameter_dims,
    equilibrium::composition::{EquilibriumEntropySolveV1, EquilibriumSolveV1, IdealGasEquilibriumV1}};
use fs_evidence::ValidityDomain;
use fs_matdb::{ConstitutiveModelCard, InitialStatePolicy, LawId, LawParameter, Provenance};
use fs_qty::{Amount, MolarMass, Temperature, Pressure};
use std::{collections::BTreeMap, sync::Arc};
const R: f64 = fs_thermochem::UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K;
fn id(s: &str) -> SpeciesId { SpeciesId::new(s).unwrap() }
fn model(name: &str, mass: f64, cp: f64, h: f64, s: f64) -> Nasa9StandardStateModelV1 {
    let mut parameters = BTreeMap::new();
    for ((name, dims), value) in required_nasa9_parameter_dims().into_iter().take(9)
        .zip([0.0, 0.0, cp, 0.0, 0.0, 0.0, 0.0, h, s]) {
        parameters.insert(name.into(), LawParameter { value, dims });
    }
    parameters.insert("reference_pressure".into(), LawParameter { value: 100000.0, dims: Pressure::DIMS });
    let card = ConstitutiveModelCard { law: LawId(fs_thermochem::NASA9_LAW_ID_V1.into()),
        law_version: fs_thermochem::NASA9_LAW_VERSION_V1, parameters,
        state_schema_version: fs_thermochem::NASA9_STATE_SCHEMA_VERSION_V1,
        initial_state: InitialStatePolicy::ZeroInternalState,
        validity: ValidityDomain::unconstrained().with("T", 200.0, 6000.0), sources: Vec::new(),
        provenance: Provenance { source: "manufactured reacting propulsion".into(),
            license: "test fixture".into(), artifact: None } };
    Nasa9StandardStateModelV1::new(id(name), MolarMass::new(mass),
        StandardStateConventionV1::new(StandardStatePhaseV1::Gas, ReferenceEquationOfStateV1::IdealGas,
            Pressure::new(100000.0), ElementalReferenceIdV1::new("propulsion-test").unwrap()).unwrap(),
        vec![Nasa9RegionV1::from_card(card).unwrap()]).unwrap()
}
fn chemistry(reacting: bool) -> Arc<IdealGasEquilibriumV1> {
    let models = if reacting { vec![model("H", 0.001, 2.5, 12000.0, 1.0), model("H2", 0.002, 3.5, 0.0, 0.0)] }
        else { vec![model("H2", 0.002, 3.5, 0.0, 0.0)] };
    let ids: Vec<_> = models.iter().map(|m| m.species().clone()).collect();
    let counts = if reacting { vec![1, 2] } else { vec![2] };
    Arc::new(IdealGasEquilibriumV1::new(
        ElementalMatrix::new(vec![ElementId::new("H").unwrap()], ids.clone(), vec![counts]).unwrap(),
        ChargeVector::new(ids.clone(), vec![0; ids.len()]).unwrap(), models).unwrap())
}
fn controls() -> EquilibriumRocketBlowdownSolve {
    EquilibriumRocketBlowdownSolve { reservoir: EquilibriumBlowdownSolve {
        minimum_pressure_pa: 1000.0, relative_density_tolerance: 2e-11,
        mass_tolerance_kg: 1e-13, energy_tolerance_j: 1e-5, maximum_density_iterations: 80 },
        nozzle: EquilibriumRocketNozzleSolve { thermochemistry: EquilibriumEntropySolveV1 {
            equilibrium: EquilibriumSolveV1 { relative_balance_tolerance: 1e-12,
                dimensionless_potential_tolerance: 1e-11, ..EquilibriumSolveV1::default() },
            minimum_temperature_k: 200.0, maximum_temperature_k: 6000.0,
            entropy_tolerance_j_k: 1e-9, maximum_temperature_iterations: 100 },
            minimum_pressure_pa: 1000.0, sonic_tolerance: 1e-9, relative_mass_flow_tolerance: 1e-9,
            energy_tolerance_j_kg: 1e-6, maximum_pressure_iterations: 100 },
        maximum_outgoing_fraction: 0.5, maximum_mass_iterations: 32 }
}
fn reservoir(reacting: bool) -> EquilibriumGasReservoir {
    EquilibriumGasReservoir::new(chemistry(reacting), vec![(id("H2"), Amount::new(1.0))],
        0.02, Temperature::new(2000.0), Pressure::new(100000.0), controls().nozzle.thermochemistry, || false).unwrap()
}
fn nozzle() -> RocketNozzle { RocketNozzle { throat_area_m2: 1e-5, exit_area_m2: 5e-5 } }
fn near(a: f64, b: f64, tolerance: f64) {
    assert!((a-b).abs() <= tolerance*b.abs().max(1e-30), "{a:.17e} != {b:.17e}");
}

#[test]
fn physical_time_depletion_and_impulse_are_second_order_against_perfect_gas_solution() {
    let initial = reservoir(false);
    let r = R / 0.002;
    let rate0 = 100000.0 * 1e-5 / (r * 2000.0).sqrt()
        * 1.4_f64.sqrt() * (5.0_f64 / 6.0).powi(3);
    let (mut lo, mut hi) = (1.0_f64, 10.0_f64);
    for _ in 0..100 {
        let m = (lo + hi) * 0.5;
        let area_ratio = ((1.0 + 0.2*m*m) / 1.2).powi(3) / m;
        if area_ratio < 5.0 { lo = m; } else { hi = m; }
    }
    let mach = (lo + hi) * 0.5;
    let exit_t = 2000.0 / (1.0 + 0.2*mach*mach);
    let exit_p = 100000.0 * (exit_t / 2000.0).powf(3.5);
    let force0 = rate0 * mach * (1.4*r*exit_t).sqrt() + exit_p*5e-5;
    let k = rate0 / initial.mass_kg();
    let duration = 0.25;
    let fraction = (1.0 + 0.2*k*duration).powf(-5.0);
    let exact_mass = initial.mass_kg()*fraction;
    let exact_impulse = force0/k * (1.0 - fraction.powf(1.2)) / 1.2;
    let mut errors = Vec::new();
    for n in [4, 8, 16] {
        let mut tank = initial.clone(); let mut impulse = 0.0; let mut energy = 0.0;
        for _ in 0..n {
            let old_mass = tank.mass_kg();
            let (endpoint, step) = tank.rocket_blowdown(0.0, nozzle(), duration/f64::from(n),
                controls(), || false).unwrap().into_parts();
            near(old_mass-endpoint.mass_kg(), step.gas.expelled_mass_kg, 1e-14);
            assert!(step.mass_quadrature_residual_kg.abs() <= controls().reservoir.mass_tolerance_kg);
            assert!(step.gas.specific_entropy_residual_j_kg_k.abs() <= 1e-9/0.002);
            assert!(step.mass_iterations <= 8);
            impulse += step.axial_impulse_n_s; energy += step.gas.exhausted_enthalpy_j; tank = endpoint;
        }
        near(energy, initial.internal_energy_j()-tank.internal_energy_j(), 2e-8);
        errors.push(((tank.mass_kg()-exact_mass).abs(), (impulse-exact_impulse).abs()));
    }
    for pair in errors.windows(2) {
        assert!((3.4..4.6).contains(&(pair[0].0/pair[1].0)), "mass {errors:?}");
        assert!((3.4..4.6).contains(&(pair[0].1/pair[1].1)), "impulse {errors:?}");
    }
    assert!(errors[2].0/initial.mass_kg() < 1e-5);
    assert!(errors[2].1/exact_impulse < 3e-5);
}

#[test]
fn reacting_time_step_reuses_the_original_nozzle_and_preserves_restart() {
    let initial = reservoir(true); let duration = 0.0625;
    let result = initial.rocket_blowdown(0.0, nozzle(), duration, controls(), || false).unwrap();
    let step = result.step(); let ch = step.midpoint_nozzle.chamber().equilibrium();
    let source = chemistry(true);
    let independent = equilibrium_rocket_nozzle_flow(&source, ch.temperature(), ch.pressure(),
        &[(id("H2"), Amount::new(1.0))], 0.0, nozzle(), controls().nozzle, || false).unwrap();
    near(step.midpoint_nozzle.mass_flow_kg_s, independent.mass_flow_kg_s, 1e-7);
    near(step.midpoint_nozzle.thrust_n, independent.thrust_n, 1e-7);
    near(step.axial_impulse_n_s, duration*independent.thrust_n, 1e-7);
    near(ch.mass().value(), 0.002, 1e-10);
    near(step.gas.exhausted_enthalpy_j, initial.internal_energy_j()-result.endpoint().internal_energy_j(), 1e-8);
    assert!(result.endpoint().thermodynamics().equilibrium().temperature().value() < 2000.0);
    assert!(result.endpoint().thermodynamics().equilibrium().pressure().value() < 100000.0);
    let atoms = |tank: &EquilibriumGasReservoir| tank.thermodynamics().equilibrium().species().iter()
        .find(|s| s.species == id("H")).unwrap().mole_fraction;
    assert!(atoms(result.endpoint()) < atoms(&initial));
    assert_eq!(initial.rocket_blowdown(0.0, nozzle(), duration, controls(), || false).unwrap(), result);
    let checkpoint = result.endpoint().clone();
    assert_eq!(result.endpoint().rocket_blowdown(0.0, nozzle(), duration, controls(), || false).unwrap(),
        checkpoint.rocket_blowdown(0.0, nozzle(), duration, controls(), || false).unwrap());
}

#[test]
fn time_inventory_regime_and_work_limits_refuse_without_spending_propellant() {
    let tank = reservoir(false); let before = tank.clone(); let c = controls();
    let source = chemistry(false);
    let initial = equilibrium_rocket_nozzle_flow(&source, Temperature::new(2000.0), Pressure::new(100000.0),
        &[(id("H2"), Amount::new(1.0))], 0.0, nozzle(), c.nozzle, || false).unwrap();
    let mut limited = c; limited.maximum_outgoing_fraction = 0.01;
    assert!(matches!(tank.rocket_blowdown(0.0, nozzle(), 1.0, limited, || false),
        Err(EquilibriumRocketBlowdownError::InventoryLimit { .. })));
    limited = c; limited.maximum_mass_iterations = 1;
    assert!(matches!(tank.rocket_blowdown(0.0, nozzle(), 0.0625, limited, || false),
        Err(EquilibriumRocketBlowdownError::NoConvergence { .. })));
    // Initial and midpoint are supported; final exit pressure is below ambient.
    let ambient = 0.91*initial.exit().sample().pressure.value();
    assert!(matches!(tank.rocket_blowdown(ambient, nozzle(), 0.1, c, || false),
        Err(EquilibriumRocketBlowdownError::Nozzle(EquilibriumRocketNozzleError::UnsupportedBackPressure { .. }))));
    limited = c; limited.nozzle.thermochemistry.entropy_tolerance_j_k *= 2.0;
    assert!(tank.rocket_blowdown(0.0, nozzle(), 0.01, limited, || false).is_err());
    for dt in [0.0, -1.0, f64::NAN, f64::INFINITY, 1e-300] {
        assert!(tank.rocket_blowdown(0.0, nozzle(), dt, c, || false).is_err());
    }
    let mut complete_polls = 0;
    tank.rocket_blowdown(0.0, nozzle(), 0.0625, c, || { complete_polls += 1; false }).unwrap();
    for limit in [0, complete_polls/2, complete_polls-1] {
        let mut polls = 0;
        assert!(tank.rocket_blowdown(0.0, nozzle(), 0.0625, c, || { polls += 1; polls > limit }).is_err());
    }
    assert_eq!(tank, before);
}
