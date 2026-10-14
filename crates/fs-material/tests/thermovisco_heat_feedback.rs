//! G0/G1/G3: original Maxwell history, physical heat and caloric feedback.
//! Synthetic coefficients; algebraic references are not measured material data.
use core::convert::Infallible;
use std::cell::Cell;

use fs_blake3::ContentHash;
use fs_material::phase::{EnthalpyPhaseKnot, EquilibriumEnthalpyPhaseCurve,
    UniformEnthalpyStep, UniformEnthalpyStepInput};
use fs_material::thermovisco::{Limits, ThermalKnot, ThermoMaxwell};
use fs_material::thermovisco::thermal::{CaloricPoint, FeedbackConfig, FeedbackError};

const VOLUME: f64 = 1e-6;
const MASS: f64 = 1e-3;
fn config() -> FeedbackConfig {
    FeedbackConfig { max_iterations: 100, temperature_tolerance_k: 1e-10,
        energy_tolerance_j: 1e-9, relaxation: 1.0 }
}
fn law(hot_relaxation: f64) -> ThermoMaxwell {
    ThermoMaxwell::new(1e8, vec![(3e8, 1.0)], vec![
        ThermalKnot { temperature_k: 280.0, modulus_scale: 1.0, relaxation_scale: 1.0, free_strain: 0.0 },
        ThermalKnot { temperature_k: 300.0, modulus_scale: 1.0, relaxation_scale: 1.0, free_strain: 0.0 },
        ThermalKnot { temperature_k: 320.0, modulus_scale: 1.0, relaxation_scale: hot_relaxation, free_strain: 0.0 },
    ], Limits { max_strain: 0.04, max_strain_rate_s: 10.0 }).unwrap()
}
fn chart(cp: f64, tag: u8) -> EquilibriumEnthalpyPhaseCurve {
    EquilibriumEnthalpyPhaseCurve::try_new(ContentHash([tag; 32]), [280.0, 300.0, 600.0]
        .into_iter().map(|t| EnthalpyPhaseKnot {
            specific_enthalpy_j_kg: (t - 300.0) * cp, temperature_k: t,
            liquid_mass_fraction: 0.0, bulk_density_kg_m3: 1000.0,
        }).collect()).unwrap()
}
fn adiabatic(input: UniformEnthalpyStepInput<'_>) -> Result<UniformEnthalpyStep<()>, Infallible> {
    let h = input.initial.specific_enthalpy_j_kg() + input.internal_heat_j / input.mass_kg;
    Ok(UniformEnthalpyStep { state: input.curve.state_at_specific_enthalpy(h).unwrap(),
        external_heat_j: 0.0, energy_residual_tolerance_j: 1e-10, report: () })
}
fn close(a: f64, b: f64, tolerance: f64) {
    assert!((a-b).abs() <= tolerance, "actual={a}, expected={b}");
}

#[test]
fn only_physical_dissipation_heats_and_the_original_endpoint_balance_closes() {
    let material = law(1.0);
    let curve = chart(100.0, 1);
    let old = CaloricPoint::try_new(material.birth(300.0, 0.0).unwrap(),
        &curve, 0.0, MASS, VOLUME, 0.001).unwrap();
    for dt in [0.01, 0.2, 2.0] {
        let step = old.advance(0.02, dt, config(), adiabatic, || false).unwrap();
        let a = 1.0 / (1.0 + dt);
        let q = a * 0.02;
        let physical = VOLUME * 3e8 * dt * q * q;
        let numerical = VOLUME * 0.5 * (1e8*0.02_f64.powi(2) + 3e8*q*q);
        close(step.next.mechanical().branch_strains()[0], q, 1e-15);
        close(MASS*step.next.caloric_state().specific_enthalpy_j_kg(), physical, 1e-12);
        close(step.applied_heat_j, physical, 1e-12);
        close(VOLUME*step.mechanical_energy.numerical_dissipation, numerical, 1e-12);
        assert!(numerical > 0.0);
        assert!(step.combined_residual_j.abs() < 1e-10);
        assert_eq!(step.next.mechanical().temperature_k(), step.next.caloric_state().temperature_k());
        assert_eq!(step.next.time_s(), dt);
    }
    assert_eq!(old.time_s(), 0.0);
    assert_eq!(old.mechanical().branch_strains(), &[0.0]);
}

#[test]
fn caloric_data_alone_changes_temperature_relaxation_and_stress() {
    let material = law(0.05);
    let a = chart(5.0, 2);
    let b = chart(10.0, 3);
    let low = CaloricPoint::try_new(material.birth(300.0, 0.0).unwrap(),
        &a, 0.0, MASS, VOLUME, 0.001).unwrap()
        .advance(0.02, 0.2, config(), adiabatic, || false).unwrap();
    let high = CaloricPoint::try_new(material.birth(300.0, 0.0).unwrap(),
        &b, 0.0, MASS, VOLUME, 0.001).unwrap()
        .advance(0.02, 0.2, config(), adiabatic, || false).unwrap();
    assert!(low.next.caloric_state().temperature_k() > high.next.caloric_state().temperature_k());
    assert!(low.next.mechanical().stress_pa() < high.next.mechanical().stress_pa());
    let frozen = material.birth(300.0, 0.0).unwrap().advance(300.0, 0.02, 0.2, 0.001).unwrap();
    assert!(high.next.mechanical().branch_strains()[0] < frozen.point.branch_strains()[0]);
    for step in [&low, &high] {
        assert!(step.iterations > 1);
        assert!(step.heat_defect_j < config().energy_tolerance_j);
        assert!(step.thermal_residual_j.abs() < config().energy_tolerance_j);
    }
}

#[test]
fn nonlinear_trials_never_age_the_old_state_and_keep_nonclone_transport_reports() {
    struct Report(usize);
    let material = law(0.05);
    let curve = chart(5.0, 4);
    let loaded = material.birth(300.0, 0.0).unwrap().advance(300.0, 0.01, 0.01, 0.001).unwrap().point;
    let old = CaloricPoint::try_new(loaded.clone(), &curve, 0.0, MASS, VOLUME, 0.001).unwrap();
    let mut calls = 0;
    let step = old.advance(0.02, 0.2, config(), |input| {
        calls += 1;
        assert_eq!(input.initial, *old.caloric_state());
        assert_eq!(input.duration_s, 0.2);
        assert_eq!(input.mass_kg, MASS);
        let t = adiabatic(input).unwrap();
        Ok::<_, Infallible>(UniformEnthalpyStep { state: t.state, external_heat_j: t.external_heat_j,
            energy_residual_tolerance_j: t.energy_residual_tolerance_j, report: Report(calls) })
    }, || false).unwrap();
    assert_eq!(calls, step.iterations);
    assert_eq!(step.transport.report.0, calls);
    let expected = loaded.advance(step.next.caloric_state().temperature_k(), 0.02, 0.2,
        config().energy_tolerance_j/VOLUME).unwrap();
    assert_eq!(step.next.mechanical(), &expected.point);
    assert_eq!(old.mechanical(), &loaded);
}

#[test]
fn lossless_spring_does_not_turn_numerical_dissipation_into_heat() {
    let material = ThermoMaxwell::new(1e8, vec![], vec![
        ThermalKnot { temperature_k: 280.0, modulus_scale: 1.0, relaxation_scale: 1.0, free_strain: 0.0 },
        ThermalKnot { temperature_k: 600.0, modulus_scale: 1.0, relaxation_scale: 1.0, free_strain: 0.0 },
    ], Limits { max_strain: 0.04, max_strain_rate_s: 10.0 }).unwrap();
    let curve = chart(100.0, 5);
    let old = CaloricPoint::try_new(material.birth(300.0, 0.0).unwrap(), &curve, 0.0, MASS, VOLUME, 0.001).unwrap();
    let step = old.advance(0.02, 0.2, config(), adiabatic, || false).unwrap();
    assert!(step.mechanical_energy.numerical_dissipation > 0.0);
    assert_eq!(step.applied_heat_j, 0.0);
    assert_eq!(step.next.caloric_state().temperature_k(), 300.0);
    assert_eq!(step.iterations, 1);
}

#[test]
fn foreign_states_heat_creation_and_transport_failures_cannot_publish() {
    let material = law(1.0);
    let curve = chart(100.0, 6);
    let other = chart(100.0, 7);
    let old = CaloricPoint::try_new(material.birth(300.0, 0.0).unwrap(), &curve, 0.0, MASS, VOLUME, 0.001).unwrap();
    let saved = old.clone();
    let bad = old.advance(0.01, 0.2, config(), |input| {
        let mut proposal = adiabatic(input).unwrap();
        proposal.state = other.state_at_specific_enthalpy(proposal.state.specific_enthalpy_j_kg()).unwrap();
        Ok::<_, Infallible>(proposal)
    }, || false);
    assert!(matches!(bad, Err(FeedbackError::Invalid(_))));
    let invented = old.advance(0.01, 0.2, config(), |input| {
        let mut proposal = adiabatic(input).unwrap();
        proposal.external_heat_j = 1.0;
        Ok::<_, Infallible>(proposal)
    }, || false);
    assert!(matches!(invented, Err(FeedbackError::Invalid(_))));
    let failed = old.advance::<(), _>(0.01, 0.2, config(), |_| Err("original transport failure"), || false);
    assert!(matches!(failed, Err(FeedbackError::Transport("original transport failure"))));
    assert_eq!(old, saved);
}

#[test]
fn cancellation_work_exhaustion_and_late_material_failure_leave_both_histories_intact() {
    let material = law(0.05);
    let curve = chart(5.0, 8);
    let old = CaloricPoint::try_new(material.birth(300.0, 0.0).unwrap(), &curve, 0.0, MASS, VOLUME, 0.001).unwrap();
    let saved = old.clone();
    let flag = Cell::new(false);
    let cancelled = old.advance(0.02, 0.2, config(), |input| {
        flag.set(true); adiabatic(input)
    }, || flag.get());
    assert!(matches!(cancelled, Err(FeedbackError::Cancelled)));
    let stopped = old.advance(0.02, 0.2, FeedbackConfig { max_iterations: 1, ..config() }, adiabatic, || false);
    assert!(matches!(stopped, Err(FeedbackError::NotConverged { iterations: 1, .. })));
    let late = old.advance(0.02, 0.2, config(), |input| {
        let h = 100.0 * 5.0;
        Ok::<_, Infallible>(UniformEnthalpyStep { state: input.curve.state_at_specific_enthalpy(h).unwrap(),
            external_heat_j: MASS*h - input.internal_heat_j,
            energy_residual_tolerance_j: 1e-10, report: () })
    }, || false);
    assert!(matches!(late, Err(FeedbackError::Mechanical(_))));
    assert_eq!(old, saved);
}

#[test]
fn unsupported_caloric_physics_is_refused_instead_of_inventing_missing_work() {
    let material = law(1.0);
    let curve = chart(100.0, 9);
    for (modulus, expansion) in [(0.5, 0.0), (1.0, 0.001)] {
        let unsupported = ThermoMaxwell::new(1e8, vec![], vec![
            ThermalKnot { temperature_k: 300.0, modulus_scale: 1.0, relaxation_scale: 1.0, free_strain: 0.0 },
            ThermalKnot { temperature_k: 600.0, modulus_scale: modulus, relaxation_scale: 1.0, free_strain: expansion },
        ], Limits { max_strain: 0.04, max_strain_rate_s: 10.0 }).unwrap();
        assert!(CaloricPoint::try_new(unsupported.birth(300.0, 0.0).unwrap(), &curve, 0.0, MASS, VOLUME, 0.001).is_err());
    }
    for (fraction, density) in [(1.0, 1000.0), (0.0, 900.0)] {
        let mut knots = curve.knots().to_vec();
        let last = knots.last_mut().unwrap();
        last.liquid_mass_fraction = fraction;
        last.bulk_density_kg_m3 = density;
        let unsupported = EquilibriumEnthalpyPhaseCurve::try_new(ContentHash([10;32]), knots).unwrap();
        assert!(CaloricPoint::try_new(material.birth(300.0, 0.0).unwrap(), &unsupported, 0.0, MASS, VOLUME, 0.001).is_err());
    }
    assert!(CaloricPoint::try_new(material.birth(300.0, 0.0).unwrap(), &curve, 0.0, 2.0*MASS, VOLUME, 0.001).is_err());
}
