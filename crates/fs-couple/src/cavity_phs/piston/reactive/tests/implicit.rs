//! Analytic isomer relaxation and actual material/piston consumers. All
//! thermodynamic and kinetic data here are synthetic implementation fixtures.

use super::*;
use crate::cavity_phs::piston::reactive::implicit::ImplicitReactionSolve;
use crate::cavity_phs::piston::reactive::material_rate::{
    MaterialRateScale, ACTIVITY_RATE_SCALE_DIMS, ACTIVITY_RATE_SCALE_PROPERTY,
};
use crate::cavity_phs::piston::reactive::subcycling::ReactionSubcycling;

fn policy() -> ImplicitReactionSolve {
    ImplicitReactionSolve {
        extent_bracket_mol: [-1.0, 1.0],
        extent_tolerance_mol: 1e-11,
        maximum_iterations: 128,
    }
}
fn controls() -> GasReactionSolve {
    GasReactionSolve { maximum_consumed_fraction: 0.99, ..chemical_controls() }
}

#[test]
fn stiff_reversible_isomers_match_backward_euler_in_both_directions() {
    // A <-> B: equal molar masses, cp and formation references give constant
    // T and K=1. The independently reduced ODE is xi' = k*(n_A0-n_B0-2xi),
    // k = scale*R*T/p0; V cancels because activities are partial pressures.
    let explicit = reaction(0.0, 3.5, 1e6);
    let law = explicit.clone().with_implicit(policy()).unwrap();
    let duration = 0.01;
    for fractions in [[0.75, 0.0, 0.25], [0.0, 0.75, 0.25]] {
        let before = gas(&law, fractions);
        let original = before.clone();
        assert!(matches!(explicit.parcel_candidate(&before, duration, controls(), || false),
            Err(GasReactionError::ConsumptionLimit { .. })));
        let result = law.parcel_candidate(&before, duration, controls(), || false).unwrap();
        let report = result.step();
        let equation = report.implicit.unwrap();
        let z = duration * 1e6 * R * 600.0 / 100000.0;
        assert!(z > 256.0);
        let difference = report.species[0].initial_mol - report.species[1].initial_mol;
        let expected = z * difference / (1.0 + 2.0*z);
        near(report.extent_mol, expected, 1e-12);
        let independently_measured = report.extent_mol - duration * before.volume_m3()
            * report.rate.net_progress_rate().unwrap().value();
        near(equation.residual_mol, independently_measured, 1e-16);
        assert!(equation.residual_mol.abs() <= policy().extent_tolerance_mol);
        assert!(equation.trials <= policy().maximum_iterations + 2);
        assert!(report.species.iter().all(|s| s.final_mol >= 0.0));
        assert!(report.entropy_change_j_k >= 0.0);
        assert!(report.maximum_consumed_fraction < 0.5);
        near(result.gas().state().state().temperature, 600.0, 1e-7);
        assert_eq!(result.gas().mass_kg(), before.mass_kg());
        assert_eq!(result.gas().volume_m3(), before.volume_m3());
        assert_eq!(before, original);
        assert_eq!(result, law.parcel_candidate(&before, duration, controls(), || false).unwrap());
    }
}

fn temperature_rate(scale: f64) -> MaterialRateScale {
    let mut claims = ClaimSet::new();
    claims.insert_claim(PropertyClaim {
        key: PropertyKey::new(ACTIVITY_RATE_SCALE_PROPERTY, ACTIVITY_RATE_SCALE_DIMS),
        value: PropertyValue::Curve {
            abscissa: "temperature".into(), abscissa_dims: Temperature::DIMS,
            knots: vec![(200.0, 0.25*scale), (600.0, scale), (1800.0, 4.0*scale)],
            dims: ACTIVITY_RATE_SCALE_DIMS,
        },
        validity: ValidityDomain::unconstrained().with_quantity("temperature",
            quantity(QuantityKind::AbsoluteTemperature), 200.0, 1800.0),
        uncertainty: UncertaintyModel::Unstated,
        interpolation: InterpolationPolicy::LinearInside,
        observations: vec![], provenance: provenance(),
    }).unwrap();
    let material = MaterialCard::assemble(MaterialStateId {
        chemistry: "synthetic-isomer-kinetic-law".into(), phase: "gas".into(),
        process: "manufactured-temperature-law".into(), revision: 0,
    }, claims, vec![]).unwrap();
    MaterialRateScale::new(material, MaterialPropertySelection::SingleClaimOnly,
        Pressure::new(100000.0)).unwrap()
}

#[test]
fn endpoint_formation_energy_and_source_rate_change_actual_pressure_and_sound() {
    let fixed = reaction(4000.0, 4.0, 1e6).with_implicit(policy()).unwrap();
    let sourced = reaction(4000.0, 4.0, 1e6)
        .with_material_rate(temperature_rate(1e6)).unwrap()
        .with_implicit(policy()).unwrap();
    let before = gas(&sourced, [0.75, 0.0, 0.25]);
    let original_u = before.mass_kg() * before.state().specific_internal_energy_j_kg();
    let result = sourced.parcel_candidate(&before, 5e-5, controls(), || false).unwrap();
    let constant = fixed.parcel_candidate(&before, 5e-5, controls(), || false).unwrap();
    let report = result.step();
    let n: Vec<_> = report.species.iter().map(|s| s.final_mol).collect();
    let expected_t = (original_u + 4000.0*n[1]) / (R * (2.5*n[0] + 3.0*n[1] + 1.5*n[2]));
    near(result.gas().state().state().temperature, expected_t, 1e-7);
    near(result.gas().state().state().pressure,
        n.iter().sum::<f64>() * R * expected_t / before.volume_m3(), 1e-6);
    let scale = report.rate_material.as_ref().unwrap()
        .property(ACTIVITY_RATE_SCALE_PROPERTY).unwrap().value_si();
    near(scale, 1e6 * (1.0 + (expected_t - 600.0)/400.0), 1e-5);
    assert!(scale > 1e6);
    assert!(report.extent_mol > constant.step().extent_mol);
    assert!(result.gas().state().state().sound_speed > before.state().state().sound_speed);
    assert_ne!(result.gas().state().state().thermal_conductivity,
        before.state().state().thermal_conductivity);
    assert!(report.stored_energy_change_j.abs() <= before.mass_kg()*controls().gas.energy_tolerance_j_kg);
    assert!(report.entropy_change_j_k > 0.0);
}

#[test]
fn existing_temporal_refinement_reuses_implicit_steps_and_matches_relaxation() {
    let law = reaction(0.0, 3.5, 1000.0).with_implicit(policy()).unwrap()
        .with_subcycling(ReactionSubcycling {
            absolute_amount_mol: 1e-6, absolute_temperature_k: 1e-5,
            absolute_pressure_pa: 0.01, relative_tolerance: 0.0,
            maximum_substeps: 256, maximum_evaluations: 511,
        }).unwrap();
    let before = gas(&law, [0.75, 0.0, 0.25]);
    let duration = 0.002;
    let result = law.parcel_candidate(&before, duration, controls(), || false).unwrap();
    let report = result.step();
    let refined = report.subcycling.as_ref().unwrap();
    assert!(report.implicit.is_none());
    assert!(refined.steps.iter().all(|step| step.reaction.implicit.is_some()));
    assert!(refined.steps.len() > 2);
    assert!(refined.estimated_error_ratio <= 1.0);
    assert_eq!(refined.steps.last().unwrap().end_s, duration);
    let k = 1000.0 * R * 600.0 / 100000.0;
    let exact = 0.5 * report.species[0].initial_mol * (1.0 - (-2.0*k*duration).exp());
    near(report.species[1].final_mol, exact, 2e-6);
}

#[test]
fn absent_reagents_equilibrium_and_failed_brackets_preserve_inventory() {
    let law = reaction(0.0, 3.5, 1e6).with_implicit(policy()).unwrap();
    for fractions in [[0.0, 0.0, 1.0], [0.5, 0.5, 0.0]] {
        let before = gas(&law, fractions);
        let result = law.parcel_candidate(&before, 1.0, controls(), || false).unwrap();
        assert_eq!(result.gas(), &before);
        assert_eq!(result.step().extent_mol, 0.0);
        assert_eq!(result.step().implicit.unwrap().trials, 1);
    }
    let constrained = reaction(0.0, 3.5, 1e6).with_implicit(ImplicitReactionSolve {
        extent_bracket_mol: [-1e-8, 1e-8], ..policy()
    }).unwrap();
    let before = gas(&constrained, [1.0, 0.0, 0.0]);
    let original = before.clone();
    assert!(matches!(constrained.parcel_candidate(&before, 0.01, controls(), || false),
        Err(GasReactionError::ImplicitOutsideBracket { .. })));
    assert_eq!(before, original);
    assert!(law.clone().with_implicit(policy()).is_err());
    assert!(GasReaction::combine(vec![law.clone(), law]).is_err());
}

#[test]
fn source_envelope_and_nonlinear_work_limits_are_not_bypassed() {
    let law = reaction(4000.0, 4.0, 1e6).with_implicit(policy()).unwrap();
    let before = gas(&law, [0.75, 0.0, 0.25]);
    let original = before.clone();
    let restricted = GasReactionSolve {
        gas: GasEnergySolve { maximum_temperature_k: 600.01, ..controls().gas },
        ..controls()
    };
    assert!(matches!(law.parcel_candidate(&before, 5e-5, restricted, || false),
        Err(GasReactionError::Gas(_))));
    let limited = reaction(4000.0, 4.0, 1e6).with_implicit(ImplicitReactionSolve {
        maximum_iterations: 1, extent_tolerance_mol: 1e-16, ..policy()
    }).unwrap();
    assert!(matches!(limited.parcel_candidate(&before, 5e-5, controls(), || false),
        Err(GasReactionError::ImplicitNotConverged { .. })));
    assert_eq!(before, original);
}

#[test]
fn original_material_piston_consumes_implicit_gas_and_cancels_transactionally() {
    let law = reaction(4000.0, 4.0, 1e6).with_implicit(policy()).unwrap();
    let piston = MaterialGasPiston::new(gas(&law, [0.75, 0.0, 0.25]),
        body(8000.0), 0.01, 0.0).unwrap();
    let original = piston.clone();
    let mut calls = 0;
    let result = piston.reactive_candidate(&law, drive(), controls(), mechanical_controls(), || {
        calls += 1; false
    }).unwrap();
    assert!(result.step().reaction.implicit.is_some());
    assert!(result.piston().piston().velocity_m_s() > 0.0);
    assert!(result.step().total_energy_residual_j.abs() <= mechanical_controls().energy_tolerance_j);
    assert_eq!(piston, original);
    for cancel_at in [1, calls/2, calls] {
        let mut observed = 0;
        assert!(matches!(piston.reactive_candidate(&law, drive(), controls(), mechanical_controls(), || {
            observed += 1; observed == cancel_at
        }), Err(GasReactionError::Cancelled)));
        assert_eq!(piston, original);
        assert_eq!(piston.piston().ports.next_tick(), original.piston().ports.next_tick());
    }
}
