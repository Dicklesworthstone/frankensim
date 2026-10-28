//! Actual constant-pressure chemistry and shared endpoint heat feedback.
//! Manufactured source cards; these are not empirical combustion claims.
#[path = "material_exchange/gas.rs"]
#[allow(dead_code)]
mod source;
#[path = "reactive_isobaric_heat/body.rs"]
mod body;

use fs_couple::cavity_phs::piston::reactive::{GasReaction, GasReactionSolve};
use fs_couple::cavity_phs::piston::reactive::implicit::coupled::ImplicitMechanismSolve;
use fs_couple::cavity_phs::piston::reactive::implicit::coupled::isobaric::thermal::{
    IsobaricHeatSolve, IsobaricThermalError,
};
use fs_couple::cavity_phs::piston::thermal::{ConductanceSample, HeatExchangeSolve, ThermalPistonError};
use fs_couple::cavity_phs::piston::thermal::boundary::ports::ThermalPort;
use fs_material::thermochemical_gas::{ThermochemicalGasModel, energy::{GasEnergySolve, parcel::ClosedGasParcel}};
use fs_matdb::QueryPoint;
use fs_qty::{QuantitySpec, semantic::{QuantityKind, SemanticType, ValueForm}};
use fs_solver::NewtonKrylovConfig;
use fs_thermochem::*;
use fs_thermochem::equilibrium::*;
use fs_thermochem::mixture::FrozenIdealGasMixtureModelV1;
use std::convert::Infallible;
const R: f64 = UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K;

fn law(release: f64) -> GasReaction {
    let components = [source::component("A", 0.028, 0.0), source::component("B", 0.028, -release)];
    let ids = components.iter().map(|(s, _)| s.species().clone()).collect::<Vec<_>>();
    let id = ReactionId::new("ab").unwrap();
    let matrix = StoichiometricMatrix::new(ids.clone(), vec![id.clone()], vec![vec![-1], vec![1]]).unwrap();
    let elements = ElementalMatrix::new(vec![ElementId::new("X").unwrap()], ids.clone(), vec![vec![1, 1]]).unwrap();
    let certificate = verify_conservation(&elements, &matrix, &ChargeVector::new(ids, vec![0, 0]).unwrap()).unwrap();
    let equilibrium = IdealGasReactionEquilibriumV1::new(matrix, certificate, id,
        components.iter().map(|(s, _)| s.clone()).collect()).unwrap();
    let catalog = ThermochemicalGasModel::new(FrozenIdealGasMixtureModelV1::new(
        components.iter().map(|(s, _)| s.clone()).collect(),
        Composition::new(CompositionBasis::MoleFraction, vec![0.5, 0.5]).unwrap()).unwrap(),
        components.into_iter().map(|(_, t)| t).collect(), None).unwrap();
    let law = GasReaction::new(catalog, DeclaredStoichiometricMassActionV1::new(
        ThermodynamicReverseRateClosureV1::new(equilibrium,
            ReactionProgressRateScaleV1::new(ReactionProgressRateQuantityV1::new(20.0)).unwrap()),
        MassActionKineticOrderConventionV1::CallerDeclaredStoichiometricCoefficients)).unwrap();
    let policy = ImplicitMechanismSolve {
        extent_brackets_mol: law.channels().map(|(k, _)|
            (k.reverse_rate_closure().equilibrium().reaction().clone(), [-1.0, 1.0])).collect(),
        extent_tolerance_mol: 1e-11, interior_fraction: 0.05, jacobian_relative_step: 1e-5,
        maximum_iterations: 32, maximum_evaluations: 1024,
        newton: NewtonKrylovConfig { absolute_tolerance: 1e-12, relative_tolerance: 1e-12,
            linear_restart: 8, max_linear_cycles: 8, forcing_maximum: 1e-6,
            ..NewtonKrylovConfig::default() },
    };
    law.with_implicit_mechanism(policy).unwrap()
}
fn controls() -> GasReactionSolve {
    GasReactionSolve { gas: GasEnergySolve { minimum_temperature_k: 200.0,
        maximum_temperature_k: 1800.0, energy_tolerance_j_kg: 1e-7, maximum_iterations: 128 },
        maximum_consumed_fraction: 0.99999, species_tolerance_mol: 1e-11,
        entropy_tolerance_j_k: 1e-9 }
}
fn exchange() -> IsobaricHeatSolve {
    IsobaricHeatSolve { exchange: HeatExchangeSolve { heat_bracket_j: [-40.0, 40.0],
        source_tolerance_j: 1e-9, energy_tolerance_j: 1e-7, maximum_iterations: 128 },
        maximum_chemical_evaluations: 65536 }
}
fn gas(law: &GasReaction, fractions: [f64; 2], t: f64) -> ClosedGasParcel {
    let quantity = |kind| QuantitySpec::semantic(SemanticType::new(kind, ValueForm::Static));
    let point = QueryPoint::new().with_quantity("temperature", quantity(QuantityKind::AbsoluteTemperature), t).unwrap()
        .with_quantity("pressure", quantity(QuantityKind::Pressure), 100000.0).unwrap()
        .with("specimen-lot", 42.0).unwrap();
    ClosedGasParcel::new(law.catalog().with_composition(
        Composition::new(CompositionBasis::MoleFraction, fractions.to_vec()).unwrap()).unwrap(),
        point, 0.001).unwrap()
}
fn near(a: f64, b: f64, tol: f64) { assert!((a-b).abs() <= tol, "{a} != {b}"); }

#[test]
fn endpoint_heat_uses_cp_and_actual_pressure_work_without_a_frozen_temperature() {
    let law = law(0.0);
    let old = gas(&law, [0.5, 0.5], 400.0);
    let ports = [ThermalPort { id: 1, temperature_k: 700.0 }];
    let (dt, g, external) = (0.1, 0.4, 2.0);
    let result = law.isobaric_thermal_ports_candidate(&old, dt, external, controls(), exchange(),
        &ports, 1, |port, state| Ok::<_, Infallible>(ConductanceSample {
            conductance_w_k: g, source: (port.id, state.state().temperature) }), || false).unwrap();
    let capacity = old.mass_kg() / 0.028 * 3.5 * R;
    let q = dt * g * (700.0 - 400.0 - external / capacity) / (1.0 + dt * g / capacity);
    let t = 400.0 + (external + q) / capacity;
    near(result.step().exchanged_heat_j, q, 2e-8);
    near(result.gas().state().state().temperature, t, 2e-7);
    assert_eq!(result.gas().state().state().pressure, 100000.0);
    near(result.gas().volume_m3(), old.volume_m3()*t/400.0, 1e-12);
    assert_eq!(result.source()[0].source.1, result.gas().state().state().temperature);
    near(result.step().chemistry.boundary.boundary.work_by_gas_j,
        old.mass_kg()/0.028*R*(t-400.0), 2e-7);
    assert!(result.step().energy_residual_j.abs() <= exchange().exchange.energy_tolerance_j);
    assert!(result.step().enthalpy_residual_j.abs() <= exchange().exchange.energy_tolerance_j);
    assert!(result.step().contact_entropy_j_k > 0.0);
    assert_eq!(result.gas().point().axes().get("specimen-lot"), Some(&42.0));
}

#[test]
fn opposing_ports_keep_nonzero_transfers_and_entropy_at_zero_net_heat() {
    let law = law(0.0);
    let old = gas(&law, [0.5, 0.5], 400.0);
    let ports = [ThermalPort { id: 20, temperature_k: 600.0 },
        ThermalPort { id: 10, temperature_k: 200.0 }];
    let run = |ports: &[ThermalPort]| law.isobaric_thermal_ports_candidate(&old, 0.1, 0.0,
        controls(), exchange(), ports, 2, |port, _| Ok::<_, Infallible>(ConductanceSample {
            conductance_w_k: 0.02, source: port.id }), || false).unwrap();
    let result = run(&ports);
    assert_eq!(result, run(&[ports[1], ports[0]]));
    assert_eq!(result.step().exchanged_heat_j, 0.0);
    near(result.source()[0].heat_j, -0.4, 1e-14);
    near(result.source()[1].heat_j, 0.4, 1e-14);
    near(result.step().contact_entropy_j_k, 0.4*(1.0/200.0-1.0/600.0), 1e-15);
    assert_eq!(result.step().exchange_trials, 1);
    assert_eq!(result.gas(), &old);
}

#[test]
fn exothermic_endpoint_closes_burning_heat_and_expansion_at_the_same_state() {
    let law = law(6000.0);
    let old = gas(&law, [1.0, 0.0], 400.0);
    let ports = [ThermalPort { id: 3, temperature_k: 350.0 }];
    let result = law.isobaric_thermal_ports_candidate(&old, 0.1, 0.0, controls(), exchange(),
        &ports, 1, |_, state| Ok::<_, Infallible>(ConductanceSample {
            conductance_w_k: 0.03*(state.state().temperature/400.0),
            source: state.state().temperature }), || false).unwrap();
    assert!(result.step().exchanged_heat_j < 0.0);
    let end = result.gas();
    let q = result.step().exchanged_heat_j;
    let t = end.state().state().temperature;
    near(q, 0.1*0.03*(t/400.0)*(350.0-t), exchange().exchange.source_tolerance_j);
    let chemical = &result.step().chemistry.boundary.reaction;
    assert!(chemical.extent_mol > 0.0);
    near(chemical.extent_mol, 0.1*end.volume_m3()*chemical.rate.net_progress_rate().unwrap().value(), 1e-11);
    near(result.step().chemistry.enthalpy_change_j, q, 1e-7);
    let original = law.isobaric_implicit_candidate(&old, 0.1, q, controls(), || false).unwrap();
    assert_eq!(end, original.gas());
    assert_eq!(&result.step().chemistry, original.step());
    assert!(result.step().chemical_evaluations >= chemical.implicit_mechanism.as_ref().unwrap().evaluations);
}

#[test]
fn shared_chemical_work_limit_is_enforced_inside_the_last_heat_trial() {
    let law = law(6000.0);
    let old = gas(&law, [1.0, 0.0], 400.0);
    let ports = [ThermalPort { id: 1, temperature_k: 700.0 }];
    let run = |limits| law.isobaric_thermal_ports_candidate(&old, 0.1, 0.0, controls(), limits,
        &ports, 1, |_, _| Ok::<_, Infallible>(ConductanceSample { conductance_w_k: 0.1, source: () }), || false);
    let complete = run(exchange()).unwrap();
    let count = complete.step().chemical_evaluations;
    assert!(count > 1);
    let exact = IsobaricHeatSolve { maximum_chemical_evaluations: count, ..exchange() };
    assert_eq!(run(exact).unwrap(), complete);
    let limited = IsobaricHeatSolve { maximum_chemical_evaluations: count-1, ..exchange() };
    match run(limited) {
        Err(IsobaricThermalError::ChemicalEvaluationLimit { evaluations, maximum }) => {
            assert_eq!(maximum, count-1);
            assert!(evaluations <= maximum);
        },
        other => panic!("shared work must refuse before exceeding its cap: {other:?}"),
    }
}

#[test]
fn bad_ports_missing_sources_and_cancellation_publish_no_candidate() {
    let law = law(0.0);
    let old = gas(&law, [0.5, 0.5], 400.0);
    let saved = old.clone();
    let port = ThermalPort { id: 1, temperature_k: 700.0 };
    let result = law.isobaric_thermal_ports_candidate(&old, 0.1, 0.0, controls(), exchange(),
        &[port], 1, |_, _| Err::<ConductanceSample<()>, _>("missing source"), || false);
    assert!(matches!(result, Err(IsobaricThermalError::Thermal(ThermalPistonError::Interface("missing source")))));
    for ports in [vec![port, port], vec![ThermalPort { id: 1, temperature_k: f64::NAN }]] {
        assert!(law.isobaric_thermal_ports_candidate(&old, 0.1, 0.0, controls(), exchange(),
            &ports, 2, |_, _| Ok::<_, Infallible>(ConductanceSample { conductance_w_k: 0.0, source: () }), || false).is_err());
    }
    let mut count = 0;
    law.isobaric_thermal_ports_candidate(&old, 0.1, 0.0, controls(), exchange(),
        &[port], 1, |_, _| Ok::<_, Infallible>(ConductanceSample { conductance_w_k: 0.4, source: () }),
        || { count += 1; false }).unwrap();
    for stop in [1, count / 2, count] {
        let mut polls = 0;
        assert!(matches!(law.isobaric_thermal_ports_candidate(&old, 0.1, 0.0, controls(), exchange(),
            &[port], 1, |_, _| Ok::<_, Infallible>(ConductanceSample { conductance_w_k: 0.4, source: () }),
            || { polls += 1; polls == stop }), Err(IsobaricThermalError::Cancelled)));
    }
    assert_eq!(old, saved);
}
