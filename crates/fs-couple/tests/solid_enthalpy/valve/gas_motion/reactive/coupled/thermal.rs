//! G1/G3/G4 definitions: material heat feeds the same reaction/elastic endpoint.
//! Synthetic source data and independent caloric/kinetic/structural equations;
//! no empirical engine, flame or thermal-expansion qualification is implied.
use super::*;
use fs_couple::cavity_phs::modal::reactive::thermal::{
    CoupledReactiveModalThermalSolve, ReactiveModalThermalError,
};
use fs_couple::cavity_phs::piston::thermal::{ConductanceSample, HeatExchangeSolve};
use fs_couple::cavity_phs::piston::thermal::layers::{LayerConductionError, MaterialConductionLayers};
use fs_material::phase::UniformEnthalpyStepInput;
use fs_material::state_point::THERMAL_CONDUCTIVITY_DIMS;

const BODY_MASS: f64 = 0.001;
const BODY_T0: f64 = 500.0;
const GAS_HEAT: f64 = 0.1;
const BODY_HEAT: f64 = 0.2;

fn thermal_controls() -> CoupledReactiveModalThermalSolve {
    CoupledReactiveModalThermalSolve {
        coupled: solve(),
        exchange: HeatExchangeSolve { heat_bracket_j: [-2.0, 2.0],
            source_tolerance_j: 1e-10, energy_tolerance_j: 1e-8, maximum_iterations: 96 },
        entropy_tolerance_j_k: 1e-10,
    }
}

fn thermal_point(t: f64) -> QueryPoint {
    QueryPoint::new().with_quantity("T", QuantitySpec::dimensional(Temperature::DIMS), t)
        .unwrap().with("lot", 7.0).unwrap()
}

fn thermal_material(k: f64, cp: f64) -> (MaterialCard, PropertyKey, ResolvedHeatCapacityCurve) {
    let conductivity = PropertyKey::new("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS);
    let capacity = PropertyKey::new("specific_heat_capacity", SPECIFIC_HEAT_CAPACITY_DIMS);
    let density = PropertyKey::new("density", Density::DIMS);
    let mut claims = ClaimSet::new();
    for (key, lo, hi) in [(conductivity.clone(), k, 2.0 * k),
        (capacity.clone(), cp, cp), (density.clone(), 4000.0, 4000.0)] {
        let dims = key.quantity().dims();
        claims.insert_claim(PropertyClaim { key,
            value: PropertyValue::Curve { abscissa: "T".into(),
                abscissa_dims: Temperature::DIMS, knots: vec![(200.0, lo), (800.0, hi)], dims },
            validity: ValidityDomain::unconstrained().with_quantity("T",
                QuantitySpec::dimensional(Temperature::DIMS), 200.0, 800.0).with("lot", 7.0, 7.0),
            interpolation: InterpolationPolicy::LinearInside,
            uncertainty: UncertaintyModel::Unstated, observations: vec![],
            provenance: Provenance { source: "manufactured finite chemical-cavity heat store".into(),
                license: "CC0-1.0".into(), artifact: None },
        }).unwrap();
    }
    let card = MaterialCard::assemble(MaterialStateId { chemistry: "synthetic thermal store".into(),
        phase: "source chart".into(), process: "authored".into(), revision: 0 }, claims, vec![]).unwrap();
    let caloric = ResolvedHeatCapacityCurve::try_new(&card, HeatCapacitySourceQuery {
        heat_capacity: &capacity, density: &density, temperature_axis: "T",
        lower: &thermal_point(200.0), upper: &thermal_point(800.0),
        selection: SelectionPolicy::SingleClaimOnly,
    }, SolidLiquidPhase::Solid, 0.0, 1e-6).unwrap();
    (card, conductivity, caloric)
}

fn layers<'a>(card: &'a MaterialCard, key: &'a PropertyKey, point: &'a QueryPoint)
    -> MaterialConductionLayers<'a> {
    MaterialConductionLayers { material: card, context: point, temperature_axis: "T",
        liquid_fraction_axis: None, conductivity: key,
        selection: MaterialPropertySelection::SingleClaimOnly,
        area_m2: 0.2, gas_distance_m: 0.001, body_distance_m: 0.02 }
}

fn body(curve: &EquilibriumEnthalpyPhaseCurve, h: f64) -> UniformEnthalpyStepInput<'_> {
    UniformEnthalpyStepInput { curve, initial: curve.state_at_specific_enthalpy(h).unwrap(),
        mass_kg: BODY_MASS, volume_m3: BODY_MASS / 4000.0, surface_area_m2: 0.2,
        internal_heat_j: BODY_HEAT, duration_s: 1.0 / RATE as f64 }
}

// Eliminate gas reaction and held-force response with the independent parent
// oracle. Eliminate finite body storage with its independently authored Cp.
// Only the one scalar physical transfer Q remains. Both equal-species native
// transport sources have k_gas=0.03 W/(m K); no evaluated result is borrowed.
fn thermal_reference(k: f64, cp: f64) -> (f64, f64, JointReference) {
    let sample = |q: f64| {
        let gas = joint_reference(3000.0, 1000.0, GAS_HEAT + q);
        let tb = BODY_T0 + (BODY_HEAT - q) / (BODY_MASS * cp);
        let kb = k * (1.0 + (tb - 200.0) / 600.0);
        let conductance = 0.2 / (0.001 / 0.03 + 0.02 / kb);
        (q - conductance * (tb - gas.temperature) / RATE as f64, tb, gas)
    };
    let q = bracketed_root(-2.0, 2.0, |q| sample(q).0);
    let (_, tb, gas) = sample(q);
    (q, tb, gas)
}

#[test]
fn source_conductivity_and_capacity_change_the_same_chemical_motion_endpoint() {
    with_cx(|cx| {
        let mut pressures = Vec::new();
        let mut heats = Vec::new();
        for (k, cp) in [(0.2, 50.0), (2.0, 50.0), (2.0, 100.0)] {
            let (reaction, gas) = reacting_gas(3000.0, [1.0, 0.0, 0.0]);
            let reaction = implicit(reaction);
            let run = GasValveMotion::new(gas,
                motion(1000.0, 0.0, 0.0, COUPLED_AREA, RATE, cx), cx).unwrap();
            let (card, key, caloric) = thermal_material(k, cp);
            let point = thermal_point(BODY_T0);
            let input = body(caloric.curve(), cp * (BODY_T0 - 200.0));
            let next = run.coupled_reactive_thermal_layers_candidate(
                &reaction, GAS_HEAT, OUTSIDE, input, chemistry(), thermal_controls(),
                layers(&card, &key, &point), cx,
            ).unwrap();
            let (q, tb, expected) = thermal_reference(k, cp);
            let report = next.step();
            close(report.exchanged_heat_j, q, 2e-9);
            close(next.thermal().state.temperature_k(), tb, 1e-7);
            close(next.motion().gas().state().state().temperature, expected.temperature, 3e-7);
            close(next.motion().gas().state().state().pressure, expected.pressure, 1e-5);
            close(report.coupled.reaction.reaction.extent_mol, expected.extent, 3e-13);
            close(next.motion().motion().nodal_motion(cx).unwrap().velocity_m_s[7][0],
                expected.velocity, 1e-9);
            close(report.body_entropy_change_j_k, BODY_MASS * cp * (tb / BODY_T0).ln(), 1e-10);
            close(report.interface_entropy_j_k, q * (1.0 / expected.temperature - 1.0 / tb), 1e-10);
            close(report.total_entropy_residual_j_k, expected.entropy
                + BODY_MASS * cp * (tb / BODY_T0).ln()
                - GAS_HEAT / expected.temperature - BODY_HEAT / tb, 1e-8);
            close(report.total_energy_residual_j, 0.0, 1e-8);
            assert_eq!(next.thermal().external_heat_j, -report.exchanged_heat_j);
            assert_eq!(next.thermal().report,
                layers(&card, &key, &point).resolve(next.motion().gas().state(), next.thermal().state)
                    .unwrap().source);
            assert_eq!(next.motion().motion().accepted_steps(), 1);
            assert_eq!(run.motion().accepted_steps(), 0);
            assert!(report.exchange_trials >= 3);
            assert!(report.coupled_trials > report.coupled.trials);
            pressures.push(expected.pressure);
            heats.push(report.exchanged_heat_j);
        }
        assert!(heats[1] > 2.0 * heats[0]);
        assert!((pressures[1] - pressures[0]).abs() > 10.0);
        assert!((pressures[2] - pressures[1]).abs() > 0.1);
    });
}

#[test]
fn latent_heat_store_and_native_reaction_motion_continue_on_one_clock() {
    with_cx(|cx| {
        let (reaction, gas) = reacting_gas(3000.0, [1.0, 0.0, 0.0]);
        let reaction = implicit(reaction);
        let mut run = GasValveMotion::new(gas,
            motion(1000.0, 0.0, 0.05, COUPLED_AREA, RATE, cx), cx).unwrap();
        let initial_u = run.gas().state().specific_internal_energy_j_kg();
        let initial_identity = run.motion().valve().identity();
        let (card, key, _) = thermal_material(2.0, 50.0);
        let knot = |h, t, f| EnthalpyPhaseKnot { specific_enthalpy_j_kg: h,
            temperature_k: t, liquid_mass_fraction: f, bulk_density_kg_m3: 4000.0 };
        let curve = EquilibriumEnthalpyPhaseCurve::try_new(card.content_hash(), vec![
            knot(0.0, 200.0, 0.0), knot(15_000.0, 500.0, 0.0),
            knot(25_000.0, 500.0, 1.0), knot(40_000.0, 800.0, 1.0),
        ]).unwrap();
        let mut input = body(&curve, 20_000.0);
        input.internal_heat_j = 0.0;
        let point = thermal_point(BODY_T0);
        let mut outgoing = 0.0;
        let mut heat = 0.0;
        for tick in 1..=6 {
            let previous_fraction = input.initial.liquid_mass_fraction();
            let candidate = run.coupled_reactive_thermal_layers_candidate(
                &reaction, GAS_HEAT, OUTSIDE, input, chemistry(), thermal_controls(),
                layers(&card, &key, &point), cx,
            ).unwrap();
            let (motion, thermal, step) = candidate.into_parts();
            close(thermal.state.temperature_k(), BODY_T0, 1e-12);
            assert!(thermal.state.liquid_mass_fraction() < previous_fraction);
            close(step.body_entropy_change_j_k, -step.exchanged_heat_j / BODY_T0, 1e-11);
            outgoing += step.coupled.motion.modal.viscous_dissipation_j
                + step.coupled.exterior_pressure_work_j;
            heat += step.exchanged_heat_j;
            input.initial = thermal.state;
            run = motion;
            assert_eq!(run.motion().accepted_steps(), tick);
            assert_eq!(run.motion().valve().identity(), initial_identity);
            close(run.gas().mass_kg() * (run.gas().state().specific_internal_energy_j_kg() - initial_u)
                + BODY_MASS * (input.initial.specific_enthalpy_j_kg() - 20_000.0)
                + energy(run.motion()) + outgoing - tick as f64 * GAS_HEAT, 0.0, 1e-7);
        }
        close(BODY_MASS * (20_000.0 - input.initial.specific_enthalpy_j_kg()), heat, 1e-9);
    });
}

#[test]
fn late_interface_failure_clock_source_and_budget_refusals_leave_all_histories_retryable() {
    with_cx(|cx| {
        let (reaction, gas) = reacting_gas(3000.0, [1.0, 0.0, 0.0]);
        let reaction = implicit(reaction);
        let run = GasValveMotion::new(gas,
            motion(1000.0, 0.0, 0.0, COUPLED_AREA, RATE, cx), cx).unwrap();
        let original = run.clone();
        let (card, key, caloric) = thermal_material(2.0, 50.0);
        let point = thermal_point(BODY_T0);
        let input = body(caloric.curve(), 15_000.0);
        let mut calls = 0;
        let failed = run.coupled_reactive_thermal_candidate(
            &reaction, GAS_HEAT, OUTSIDE, input, chemistry(), thermal_controls(), |g, b| {
                calls += 1;
                if calls == 2 { return Err("late interface refusal"); }
                let source = layers(&card, &key, &point).resolve(g, b).unwrap();
                Ok(ConductanceSample { conductance_w_k: source.conductance_w_k, source: () })
            }, cx,
        );
        assert!(matches!(failed, Err(ReactiveModalThermalError::Interface("late interface refusal"))));
        assert_eq!(calls, 2);
        let mut wrong_clock = input;
        wrong_clock.duration_s *= 2.0;
        assert!(run.coupled_reactive_thermal_layers_candidate(&reaction, GAS_HEAT, OUTSIDE,
            wrong_clock, chemistry(), thermal_controls(), layers(&card, &key, &point), cx).is_err());
        let other = thermal_material(0.2, 50.0).0;
        assert!(run.coupled_reactive_thermal_layers_candidate(&reaction, GAS_HEAT, OUTSIDE,
            input, chemistry(), thermal_controls(), layers(&other, &key, &point), cx).is_err());
        let mut short = thermal_controls();
        short.coupled.mechanical.maximum_trials = 1;
        assert!(run.coupled_reactive_thermal_layers_candidate(&reaction, GAS_HEAT, OUTSIDE,
            input, chemistry(), short, layers(&card, &key, &point), cx).is_err());
        assert_eq!(run.gas(), original.gas());
        assert_eq!(run.motion().model().states(), original.motion().model().states());
        assert_eq!(run.motion().accepted_steps(), 0);
        assert_eq!(input.initial, caloric.curve().state_at_specific_enthalpy(15_000.0).unwrap());
        let retry = || run.coupled_reactive_thermal_layers_candidate(&reaction, GAS_HEAT, OUTSIDE,
            input, chemistry(), thermal_controls(), layers(&card, &key, &point), cx).unwrap();
        let first = retry();
        let second = retry();
        assert_eq!(first.motion().gas(), second.motion().gas());
        assert_eq!(first.motion().motion().model().states(), second.motion().motion().model().states());
        assert_eq!(first.thermal(), second.thermal());
    });
}

#[test]
fn cancellation_after_real_reaction_motion_and_source_evaluation_publishes_nothing() {
    let gate = CancelGate::new();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        let cx = Cx::new(&gate, arena, StreamKey { seed: 3, kernel_id: 17, tile: 0, iteration: 0 },
            Budget::INFINITE, ExecMode::Deterministic);
        let (reaction, gas) = reacting_gas(3000.0, [1.0, 0.0, 0.0]);
        let reaction = implicit(reaction);
        let run = GasValveMotion::new(gas,
            motion(1000.0, 0.0, 0.0, COUPLED_AREA, RATE, &cx), &cx).unwrap();
        let original = run.clone();
        let (card, key, caloric) = thermal_material(2.0, 50.0);
        let point = thermal_point(BODY_T0);
        let input = body(caloric.curve(), 15_000.0);
        let mut calls = 0;
        let failed = run.coupled_reactive_thermal_candidate(
            &reaction, GAS_HEAT, OUTSIDE, input, chemistry(), thermal_controls(), |g, b| {
                calls += 1;
                let sample = layers(&card, &key, &point).resolve(g, b)?;
                gate.request();
                Ok::<_, LayerConductionError>(sample)
            }, &cx,
        );
        assert!(matches!(failed, Err(ReactiveModalThermalError::Cancelled)));
        assert_eq!(calls, 1);
        assert_eq!(run.gas(), original.gas());
        assert_eq!(run.motion().model().states(), original.motion().model().states());
        assert_eq!(run.motion().accepted_steps(), 0);
        assert_eq!(input.initial.specific_enthalpy_j_kg(), 15_000.0);
    });
}
