//! G0/G1/G3: original Maxwell history, physical heat and caloric feedback.
//! Synthetic coefficients; algebraic references are not measured material data.
use core::convert::Infallible;
use std::cell::Cell;

use fs_blake3::ContentHash;
use fs_material::phase::{EnthalpyPhaseKnot, EquilibriumEnthalpyPhaseCurve,
    UniformEnthalpyStep, UniformEnthalpyStepInput};
use fs_material::thermovisco::{Limits, ThermalKnot, ThermoMaxwell};
use fs_material::thermovisco::thermal::{CaloricPoint, FeedbackConfig, FeedbackError, MechanicalControl};

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
    let material = affine_law(-0.005, vec![(3e8, 1.0)]);
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
    for (middle_modulus, expansion) in [(0.75, 0.0), (1.0, 0.001)] {
        let unsupported = ThermoMaxwell::new(1e8, vec![], vec![
            ThermalKnot { temperature_k: 300.0, modulus_scale: 1.0, relaxation_scale: 1.0, free_strain: 0.0 },
            ThermalKnot { temperature_k: 400.0, modulus_scale: middle_modulus, relaxation_scale: 1.0, free_strain: 0.0 },
            ThermalKnot { temperature_k: 600.0, modulus_scale: 1.0, relaxation_scale: 1.0, free_strain: expansion },
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

fn stress(value: f64) -> MechanicalControl {
    MechanicalControl::Stress { stress_pa: value, tolerance_pa: 1e-6 }
}

#[test]
fn sustained_signed_load_matches_discrete_creep_and_closes_prescribed_work() {
    let material = law(1.0);
    let curve = chart(100.0, 11);
    for load in [-1e6, 1e6] {
        let mut point = CaloricPoint::try_new(material.birth(300.0, 0.0).unwrap(),
            &curve, 0.0, MASS, VOLUME, 0.001).unwrap();
        let mut work = 0.0;
        let mut numerical = 0.0;
        for n in 1..=12 {
            let previous = point.mechanical().total_strain();
            let step = point.advance_controlled(stress(load), 0.2, config(), adiabatic, || false).unwrap();
            // Standard-linear-solid BE creep: tau_creep=tau*(E_inf+E_1)/E_inf.
            let expected = load/1e8 + (load/4e8 - load/1e8) / 1.05_f64.powi(n);
            close(step.next.mechanical().total_strain(), expected, 1e-13);
            close(step.next.mechanical().stress_pa(), load, 1e-6);
            close(step.external_mechanical_work_j,
                VOLUME*load*(step.next.mechanical().total_strain()-previous), 1e-16);
            assert_eq!(step.control, stress(load));
            assert!(step.stress_residual_pa.unwrap().abs() <= 1e-6);
            work += step.external_mechanical_work_j;
            numerical += VOLUME*step.mechanical_energy.numerical_dissipation;
            point = step.next;
        }
        close(VOLUME*point.mechanical().stored_energy_j_m3()
            + MASS*point.caloric_state().specific_enthalpy_j_kg() + numerical, work, 1e-10);
        assert!(point.caloric_state().temperature_k() > 300.0);
    }
}

#[test]
fn caloric_data_changes_actual_strain_under_the_same_sustained_load() {
    let material = law(0.05);
    let curves = [chart(5.0, 12), chart(10.0, 13)];
    let mut results = Vec::new();
    for curve in &curves {
        let mut point = CaloricPoint::try_new(material.birth(300.0, 0.0).unwrap(),
            curve, 0.0, MASS, VOLUME, 0.001).unwrap();
        for _ in 0..12 {
            let step = point.advance_controlled(stress(2e6), 0.2, config(), adiabatic, || false).unwrap();
            assert!(step.combined_residual_j.abs() < 1e-9);
            point = step.next;
        }
        results.push(point);
    }
    assert!(results[0].caloric_state().temperature_k() > results[1].caloric_state().temperature_k()+0.5);
    assert!(results[0].mechanical().total_strain() > results[1].mechanical().total_strain()+1e-4);
    for result in &results { close(result.mechanical().stress_pa(), 2e6, 1e-6); }
}

#[test]
fn stress_control_and_strain_control_share_the_exact_constitutive_endpoint() {
    let material = law(0.05);
    let curve = chart(100.0, 14);
    let loaded = material.birth(300.0, 0.001).unwrap()
        .advance(300.0, 0.011, 0.1, 1e-3).unwrap().point;
    let old = CaloricPoint::try_new(loaded.clone(), &curve, 0.0, MASS, VOLUME, 0.001).unwrap();
    let thermostat = |input: UniformEnthalpyStepInput<'_>| {
        Ok::<_, Infallible>(UniformEnthalpyStep { state: input.initial,
            external_heat_j: -input.internal_heat_j, energy_residual_tolerance_j: 1e-10, report: () })
    };
    let a = old.advance_controlled(stress(-1e6), 0.2, config(), thermostat, || false).unwrap();
    let b = old.advance(a.next.mechanical().total_strain(), 0.2, config(), thermostat, || false).unwrap();
    assert_eq!(a.next.mechanical(), b.next.mechanical());
    assert_eq!(a.next.caloric_state(), b.next.caloric_state());
    assert_eq!(a.applied_heat_j, b.applied_heat_j);
    assert_eq!(b.stress_residual_pa, None);
    close(a.external_mechanical_work_j, b.external_mechanical_work_j, 1e-12);
    assert_eq!(old.mechanical(), &loaded);
}

#[test]
fn load_removal_preserves_memory_and_generates_heat_without_external_work() {
    let material = law(0.05);
    let curve = chart(5.0, 15);
    let mut point = CaloricPoint::try_new(material.birth(300.0, 0.0).unwrap(),
        &curve, 0.0, MASS, VOLUME, 0.001).unwrap();
    for _ in 0..5 {
        point = point.advance_controlled(stress(2e6), 0.2, config(), adiabatic, || false).unwrap().next;
    }
    let before = point.clone();
    let unloaded = point.advance_controlled(stress(0.0), 0.2, config(), adiabatic, || false).unwrap();
    assert_eq!(unloaded.external_mechanical_work_j, 0.0);
    assert!(unloaded.next.mechanical().total_strain() < point.mechanical().total_strain());
    assert!(unloaded.next.mechanical().branch_strains().iter().any(|q| q.abs() > 1e-5));
    assert!(unloaded.applied_heat_j > 0.0);
    assert!(unloaded.next.caloric_state().temperature_k() > point.caloric_state().temperature_k());
    assert!(unloaded.combined_residual_j.abs() < 1e-9);
    assert_eq!(point, before);
}

#[test]
fn invalid_or_unattainable_loads_refuse_without_changing_history() {
    let material = law(0.05);
    let curve = chart(5.0, 16);
    let old = CaloricPoint::try_new(material.birth(300.0, 0.0).unwrap(),
        &curve, 0.0, MASS, VOLUME, 0.001).unwrap();
    let saved = old.clone();
    for control in [stress(f64::NAN), stress(f64::INFINITY),
        MechanicalControl::Stress { stress_pa: 1.0, tolerance_pa: 0.0 },
        MechanicalControl::Stress { stress_pa: 1.0, tolerance_pa: f64::NAN },
        MechanicalControl::TotalStrain(f64::NAN)] {
        let called = Cell::new(false);
        assert!(matches!(old.advance_controlled(control, 0.2, config(), |input| {
            called.set(true); adiabatic(input)
        }, || false), Err(FeedbackError::Invalid(_))));
        assert!(!called.get());
    }
    assert!(matches!(old.advance_controlled(stress(1e12), 0.2, config(), adiabatic, || false),
        Err(FeedbackError::Mechanical(_))));
    let cancelled = Cell::new(false);
    assert!(matches!(old.advance_controlled(stress(2e6), 0.2, config(), |input| {
        cancelled.set(true); adiabatic(input)
    }, || cancelled.get()), Err(FeedbackError::Cancelled)));
    assert_eq!(old, saved);
}

fn affine_law(slope: f64, branches: Vec<(f64, f64)>) -> ThermoMaxwell {
    ThermoMaxwell::new(1e8, branches, [280.0, 300.0, 320.0].into_iter()
        .map(|t| ThermalKnot { temperature_k: t, modulus_scale: 1.0+slope*(t-300.0),
            relaxation_scale: 1.0, free_strain: 0.0 }).collect(),
        Limits { max_strain: 0.04, max_strain_rate_s: 10.0 }).unwrap()
}

#[test]
fn affine_elastic_heat_has_the_independent_closed_form_and_preserves_its_sign() {
    let curve = chart(100.0, 17);
    for slope in [-0.005, 0.005] {
        let material = affine_law(slope, vec![]);
        let loaded = material.birth(300.0, 0.0).unwrap()
            .advance(300.0, 0.01, 0.1, 1e-5).unwrap().point;
        let old = CaloricPoint::try_new(loaded, &curve, 0.0, MASS, VOLUME, 0.001).unwrap();
        let step = old.advance(0.02, 0.2, config(), adiabatic, || false).unwrap();
        let xi0 = 0.5*VOLUME*1e8*0.01_f64.powi(2);
        let xi1 = 0.5*VOLUME*1e8*0.02_f64.powi(2);
        let capacity = MASS*100.0;
        let t = 300.0/(1.0-slope*(xi1-xi0)/capacity);
        close(step.next.caloric_state().temperature_k(), t, 1e-9);
        close(step.reversible_heat_j, slope*t*(xi1-xi0), 1e-12);
        close(step.applied_heat_j, capacity*(t-300.0), 1e-10);
        assert_eq!(step.mechanical_energy.dissipated, 0.0);
        assert_eq!(step.reversible_heat_j.is_sign_positive(), slope.is_sign_positive());
        assert!(step.mechanical_energy.numerical_dissipation > 0.0);
        close(step.next.mechanical_entropy_j_k(), -slope*xi1, 1e-16);
        close(step.next.mechanical_internal_energy_j(), (1.0-300.0*slope)*xi1, 1e-15);
        close(step.mechanical_internal_energy_change_j, (1.0-300.0*slope)*(xi1-xi0), 1e-15);
        close(step.reversible_heat_j + VOLUME*step.mechanical_energy.modulus_work
            + step.mechanical_entropy_energy_change_j, 0.0, 1e-15);
        assert!(step.combined_residual_j.abs() < 1e-10);
    }
}

#[test]
fn heating_held_strain_changes_free_energy_without_inventing_internal_work() {
    let material = affine_law(-0.005, vec![]);
    let curve = chart(100.0, 18);
    let loaded = material.birth(300.0, 0.0).unwrap()
        .advance(300.0, 0.01, 0.1, 1e-5).unwrap().point;
    let old = CaloricPoint::try_new(loaded, &curve, 0.0, MASS, VOLUME, 0.001).unwrap();
    let step = old.advance(0.01, 1.0, config(), |input| {
        let state = input.curve.state_at_specific_enthalpy(100.0*20.0).unwrap();
        Ok::<_, Infallible>(UniformEnthalpyStep { state,
            external_heat_j: MASS*2000.0-input.internal_heat_j,
            energy_residual_tolerance_j: 1e-10, report: () })
    }, || false).unwrap();
    close(step.next.mechanical().stress_pa(), 0.9*old.mechanical().stress_pa(), 1e-8);
    close(step.mechanical_internal_energy_change_j, 0.0, 1e-16);
    close(step.reversible_heat_j, 0.0, 1e-16);
    assert!(step.mechanical_energy.modulus_work < 0.0);
    assert!(step.mechanical_entropy_energy_change_j > 0.0);
    assert_eq!(step.external_mechanical_work_j, 0.0);
    assert_eq!(step.mechanical_energy.dissipated, 0.0);
    close(step.next.mechanical_internal_energy_j(), old.mechanical_internal_energy_j(), 1e-16);
    close(step.next.mechanical_entropy_j_k(), old.mechanical_entropy_j_k(), 1e-16);
}

#[test]
fn affine_hereditary_history_matches_independent_cooling_relaxation_and_unloading_accounts() {
    let slope = -0.005;
    let intercept = 1.0-300.0*slope;
    let material = affine_law(slope, vec![(3e8, 1.0)]);
    let curve = chart(100.0, 19);
    let mut old = CaloricPoint::try_new(material.birth(300.0, 0.0).unwrap(),
        &curve, 0.0, MASS, VOLUME, 0.001).unwrap();
    let (mut q, mut strain, mut temperature, mut work, mut numerical) = (0.0, 0.0, 300.0, 0.0, 0.0);
    let capacity = MASS*100.0;
    let dt = 0.2;
    for next_strain in [0.02, 0.02, 0.0, -0.01, 0.0] {
        let next_q = (q+next_strain-strain)/(1.0+dt);
        let xi0 = 0.5*VOLUME*(1e8*strain*strain + 3e8*q*q);
        let xi1 = 0.5*VOLUME*(1e8*next_strain*next_strain + 3e8*next_q*next_q);
        let base_loss = VOLUME*3e8*dt*next_q*next_q;
        let next_t = (capacity*temperature + intercept*base_loss)
            /(capacity-slope*(base_loss+xi1-xi0));
        let scale = intercept+slope*next_t;
        let step = old.advance(next_strain, dt, config(), adiabatic, || false).unwrap();
        close(step.next.mechanical().branch_strains()[0], next_q, 1e-15);
        close(step.next.caloric_state().temperature_k(), next_t, 1e-9);
        close(VOLUME*step.mechanical_energy.dissipated, scale*base_loss, 1e-12);
        close(step.reversible_heat_j, slope*next_t*(xi1-xi0), 1e-12);
        close(step.next.mechanical_internal_energy_j(), intercept*xi1, 1e-14);
        assert!(step.mechanical_energy.dissipated >= 0.0);
        assert!(step.combined_residual_j.abs() < 1e-9);
        work += step.external_mechanical_work_j;
        numerical += VOLUME*step.mechanical_energy.numerical_dissipation;
        close(step.next.mechanical_internal_energy_j()
            + MASS*step.next.caloric_state().specific_enthalpy_j_kg() + numerical, work, 1e-9);
        q = next_q; strain = next_strain; temperature = next_t; old = step.next;
    }
}
