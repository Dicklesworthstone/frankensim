use super::*;
use fs_evidence::ValidityDomain;
use fs_matdb::{ConstitutiveModelCard, InitialStatePolicy, LawId, LawParameter, Provenance};
use fs_qty::MolarMass;
use fs_thermochem::{ChargeVector, ElementalMatrix, ElementId, ElementalReferenceIdV1,
    Nasa9StandardStateModelV1, Nasa9RegionV1, ReferenceEquationOfStateV1,
    StandardStateConventionV1, StandardStatePhaseV1, required_nasa9_parameter_dims,
    equilibrium::composition::EquilibriumSolveV1};
use std::collections::BTreeMap;
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
        law_version: 1, parameters, state_schema_version: 0, initial_state: InitialStatePolicy::ZeroInternalState,
        validity: ValidityDomain::unconstrained().with("T", 200.0, 6000.0), sources: Vec::new(),
        provenance: Provenance { source: "manufactured dissociating nozzle, not physical hydrogen calibration".into(),
            license: "test fixture".into(), artifact: None } };
    Nasa9StandardStateModelV1::new(id(name), MolarMass::new(mass),
        StandardStateConventionV1::new(StandardStatePhaseV1::Gas, ReferenceEquationOfStateV1::IdealGas,
            Pressure::new(100000.0), ElementalReferenceIdV1::new("nozzle-test").unwrap()).unwrap(),
        vec![Nasa9RegionV1::from_card(card).unwrap()]).unwrap()
}
fn chemistry(reacting: bool) -> IdealGasEquilibriumV1 {
    let models = if reacting { vec![model("H", 0.001, 2.5, 12000.0, 1.0), model("H2", 0.002, 3.5, 0.0, 0.0)] }
        else { vec![model("H2", 0.002, 3.5, 0.0, 0.0)] };
    let ids: Vec<_> = models.iter().map(|m| m.species().clone()).collect();
    let counts = if reacting { vec![1, 2] } else { vec![2] };
    IdealGasEquilibriumV1::new(ElementalMatrix::new(vec![ElementId::new("H").unwrap()], ids.clone(), vec![counts]).unwrap(),
        ChargeVector::new(ids.clone(), vec![0; ids.len()]).unwrap(), models).unwrap()
}
fn controls(scale: f64) -> EquilibriumRocketNozzleSolve {
    EquilibriumRocketNozzleSolve { thermochemistry: EquilibriumEntropySolveV1 {
        equilibrium: EquilibriumSolveV1::default(), minimum_temperature_k: 200.0,
        maximum_temperature_k: 6000.0, entropy_tolerance_j_k: 1e-8 * scale, maximum_temperature_iterations: 100 },
        minimum_pressure_pa: 1000.0, sonic_tolerance: 1e-8, relative_mass_flow_tolerance: 1e-8,
        energy_tolerance_j_kg: 1e-6, maximum_pressure_iterations: 100 }
}
fn near(a: f64, b: f64, tolerance: f64) {
    assert!((a-b).abs() <= tolerance * b.abs().max(1e-30), "{a:.17e} != {b:.17e}");
}
fn run(o: &IdealGasEquilibriumV1, scale: f64, area_scale: f64, ambient: f64) -> EquilibriumRocketNozzleFlow {
    equilibrium_rocket_nozzle_flow(o, Temperature::new(2000.0), Pressure::new(100000.0),
        &[(id("H2"), Amount::new(scale))], ambient,
        RocketNozzle { throat_area_m2: 1e-5 * area_scale, exit_area_m2: 5e-5 * area_scale },
        controls(scale), || false).unwrap()
}

#[test]
fn g1_equilibrium_nozzle_reduces_to_analytic_mach_two_perfect_gas() {
    let o = chemistry(false);
    let mut c = controls(1.0); c.minimum_pressure_pa = 2000.0;
    let flow = equilibrium_rocket_nozzle_flow(&o, Temperature::new(1800.0), Pressure::new(2e6),
        &[(id("H2"), Amount::new(1.0))], 0.0,
        RocketNozzle { throat_area_m2: 0.001, exit_area_m2: 0.001 * 27.0 / 16.0 }, c, || false).unwrap();
    let r = R / 0.002;
    let rate = 2000.0 / (r * 1800.0).sqrt() * 1.4_f64.sqrt() * (5.0_f64 / 6.0).powi(3);
    let v = 2.0 * (1.4 * r * 1000.0).sqrt();
    let pe = 2e6 * (5.0_f64 / 9.0).powf(3.5);
    near(flow.throat().mach_number, 1.0, 1e-8);
    near(flow.throat().sample().temperature.value(), 1500.0, 2e-8);
    near(flow.exit().mach_number, 2.0, 3e-8);
    near(flow.exit().sample().temperature.value(), 1000.0, 3e-8);
    near(flow.exit().sample().pressure.value(), pe, 5e-8);
    near(flow.mass_flow_kg_s, rate, 3e-8);
    near(flow.thrust_n, rate * v + pe * 0.001 * 27.0 / 16.0, 5e-8);
    near(flow.specific_impulse_s, flow.thrust_n / rate / STANDARD_GRAVITY_M_S2, 5e-8);
}

#[test]
fn g1_reacting_nozzle_matches_independent_dissociation_and_entropy_roots() {
    let o = chemistry(true);
    let flow = run(&o, 1.0, 1.0, 0.0);
    // Independent reference: alpha=sqrt(K/(K+4p/p0)),
    // K=exp(1/2-24000/T+3/2 ln(T)); solve S(T,p)=S_chamber,
    // then v^2=c_equilibrium^2 and rho*v*A_exit=mdot.
    near(flow.throat().sample().temperature.value(), 1887.9155931346877, 1e-7);
    near(flow.throat().sample().pressure.value(), 56942.976229007836, 1e-7);
    near(flow.exit().sample().temperature.value(), 1468.8960900541344, 1e-7);
    near(flow.exit().sample().pressure.value(), 3518.813432226304, 1e-7);
    near(flow.exit().velocity_m_s, 7943.681969002135, 1e-7);
    near(flow.mass_flow_kg_s, 0.00018696926420027924, 1e-7);
    near(flow.thrust_n, 1.6611650443966697, 1e-7);
    near(flow.characteristic_velocity_m_s, 5348.472671576714, 1e-7);
    let amounts = |s: &EquilibriumThermodynamicStateV1| s.equilibrium().species().iter()
        .find(|s| s.species == id("H")).unwrap().amount.value();
    near(amounts(flow.throat().state()), 0.7791604251357318, 1e-7);
    near(amounts(flow.exit().state()), 0.4482322050939418, 1e-7);
    assert!(amounts(flow.exit().state()) < amounts(flow.chamber()));
    for station in [flow.throat(), flow.exit()] {
        let solution = station.state().equilibrium();
        let independent_drop = (flow.chamber().equilibrium().enthalpy().value() - solution.enthalpy().value()) / 0.002;
        near(station.velocity_m_s.powi(2) / 2.0, independent_drop, 1e-8);
        assert!(station.sample().entropy_change_j_k.abs() <= controls(1.0).thermochemistry.entropy_tolerance_j_k);
        assert!(solution.species().iter().all(|s| s.source.is_some()));
    }
    let a = amounts(flow.throat().state()) / 2.0;
    let frozen_cp = R * (3.5 + 1.5 * a) / 0.002;
    let gas_r = R * (1.0 + a) / 0.002;
    let frozen_c = (frozen_cp / (frozen_cp - gas_r) * gas_r * flow.throat().sample().temperature.value()).sqrt();
    assert!((flow.throat().velocity_m_s / frozen_c - 1.0).abs() > 0.01,
        "reacting choking cannot use frozen sound speed");
}

#[test]
fn g3_throat_area_and_feed_amount_have_distinct_physical_scalings() {
    let o = chemistry(true);
    let base = run(&o, 1.0, 1.0, 0.0);
    let large = run(&o, 1.0, 2.0, 0.0);
    near(large.mass_flow_kg_s, 2.0 * base.mass_flow_kg_s, 1e-7);
    near(large.thrust_n, 2.0 * base.thrust_n, 1e-7);
    near(large.exit().velocity_m_s, base.exit().velocity_m_s, 1e-7);
    for amount in [1e-6, 1e6] {
        let flow = run(&o, amount, 1.0, 0.0);
        near(flow.mass_flow_kg_s, base.mass_flow_kg_s, 1e-7);
        near(flow.thrust_n, base.thrust_n, 1e-7);
    }
    assert_eq!(run(&o, 1.0, 1.0, 0.0), base);
}

#[test]
fn g1_pressure_thrust_and_sonic_exit_are_explicit() {
    let o = chemistry(true);
    let vacuum = run(&o, 1.0, 1.0, 0.0);
    let pe = vacuum.exit().sample().pressure.value();
    let matched = run(&o, 1.0, 1.0, pe);
    near(vacuum.thrust_n - matched.thrust_n, pe * 5e-5, 1e-9);
    assert_eq!(matched.pressure_thrust_n, 0.0);
    assert_eq!(matched.mass_flow_kg_s, vacuum.mass_flow_kg_s);
    let sonic = equilibrium_rocket_nozzle_flow(&o, Temperature::new(2000.0), Pressure::new(100000.0),
        &[(id("H2"), Amount::new(1.0))], 0.0,
        RocketNozzle { throat_area_m2: 1e-5, exit_area_m2: 1e-5 }, controls(1.0), || false).unwrap();
    assert_eq!(sonic.throat(), sonic.exit());
    assert_eq!(sonic.relative_mass_flow_residual, 0.0);
}

#[test]
fn g3_g4_nozzle_limits_and_cancellation_return_no_partial_flow() {
    let o = chemistry(true); let before = o.clone();
    let feed = [(id("H2"), Amount::new(1.0))];
    let nozzle = RocketNozzle { throat_area_m2: 1e-5, exit_area_m2: 5e-5 };
    let call = |c, pa| equilibrium_rocket_nozzle_flow(&o, Temperature::new(2000.0), Pressure::new(100000.0),
        &feed, pa, nozzle, c, || false);
    assert!(matches!(call(controls(1.0), 90000.0), Err(Error::UnsupportedBackPressure { .. })));
    let mut c = controls(1.0); c.minimum_pressure_pa = 90000.0;
    assert!(matches!(call(c, 0.0), Err(Error::NotBracketed(_))));
    c = controls(1.0); c.maximum_pressure_iterations = 1;
    assert!(matches!(call(c, 0.0), Err(Error::NoConvergence { .. })));
    c = controls(1.0); c.thermochemistry.minimum_temperature_k = 100.0;
    assert!(call(c, 0.0).is_err());
    for limit in [0, 20, 200] {
        let mut polls = 0;
        assert!(equilibrium_rocket_nozzle_flow(&o, Temperature::new(2000.0), Pressure::new(100000.0),
            &feed, 0.0, nozzle, controls(1.0), || { polls += 1; polls > limit }).is_err());
    }
    assert_eq!(o, before);
}
