//! Analytic caloric/reaction fixtures test physical coupling, not empirical rates.
use super::*;
use crate::material::compression::{
    MAX_COMPRESSION_SUBSTEPS, MaterialCompressionError, MaterialCompressionSolve,
    advance_closed_parcel_compression,
};
use fs_material::thermochemical_gas::energy::entropy::GasEntropySolve;

fn compression_controls(substeps: usize) -> MaterialCompressionSolve {
    let mut kinetics = controls();
    kinetics.chemistry.method = crate::chemistry::IsochoricKineticsMethod::ExactFirstOrder {
        limits: crate::chemistry::linear::FirstOrderLimits::default(),
    };
    kinetics.chemistry.maximum_step_s = 1.0;
    MaterialCompressionSolve {
        kinetics,
        volume: GasEntropySolve {
            minimum_temperature_k: 200.0,
            maximum_temperature_k: 1900.0,
            entropy_tolerance_j_kg_k: 1e-8,
            maximum_iterations: 96,
        },
        substeps,
    }
}

#[test]
fn equal_calorics_recover_adiabatic_pressure_temperature_and_work() {
    for final_volume in [0.005, 0.02] {
        let (bank, mechanism, mut parcel) = fixture(0.0, 2000.0);
        let before = parcel.clone();
        let report = advance_closed_parcel_compression(
            &mut parcel,
            &bank,
            &mechanism,
            Volume::new(final_volume),
            Time::new(0.1),
            compression_controls(8),
            || false,
        )
        .unwrap();
        let ratio = before.volume_m3() / final_volume;
        let expected_t = 400.0 * ratio.powf(1.0 / 3.0);
        let expected_p = before.state().state().pressure * ratio.powf(4.0 / 3.0);
        let expected_work = 0.01 * 3.0 * R * (expected_t - 400.0);
        assert!((parcel.state().state().temperature - expected_t).abs() < 2e-7);
        assert!((parcel.state().state().pressure / expected_p - 1.0).abs() < 1e-9);
        assert!((report.boundary_work_on_gas_j - expected_work).abs() < 1e-7);
        assert!(report.energy_residual_j.abs() < 1e-8);
        assert!(report.entropy_change_j_k > 0.0); // Chemical mixing; volume maps are isentropic.
        assert_eq!(parcel.mass_kg().to_bits(), before.mass_kg().to_bits());
        assert_eq!(parcel.volume_m3(), final_volume);
        let product = parcel
            .model()
            .thermodynamics()
            .mole_composition()
            .fractions()[1];
        assert!((product - (1.0 - (-0.4_f64).exp())).abs() < 1e-11);
        parcel.state().state().validate_coherence().unwrap();
    }
}

#[test]
fn exothermic_compression_converges_to_independent_coupled_solution() {
    // A -> B, k=4/s, equal Cv=3R, formation(B)=-Q. Mole count is constant.
    // With V(t)=V0 exp(-3*a*t), the first law is
    // dT/dt = a*T + (Q/3R)*k*exp(-k*t), not reaction followed by one stroke.
    let duration = 0.1;
    let final_volume = 0.005;
    let release = 24000.0;
    let a = (0.01_f64 / final_volume).ln() / (3.0 * duration);
    let exact_t = 400.0 * (a * duration).exp()
        + release / (3.0 * R) * 4.0 * ((a * duration).exp() - (-4.0 * duration).exp()) / (a + 4.0);
    let mut errors = Vec::new();
    for substeps in [2, 4, 8, 16] {
        let (bank, mechanism, mut parcel) = fixture(-release, 2000.0);
        let before = parcel.clone();
        let report = advance_closed_parcel_compression(
            &mut parcel,
            &bank,
            &mechanism,
            Volume::new(final_volume),
            Time::new(duration),
            compression_controls(substeps),
            || false,
        )
        .unwrap();
        errors.push((parcel.state().state().temperature - exact_t).abs());
        let actual_delta = parcel.mass_kg()
            * (parcel.state().specific_internal_energy_j_kg()
                - before.state().specific_internal_energy_j_kg());
        assert!((actual_delta - report.internal_energy_change_j).abs() < 2e-7);
        assert!((actual_delta - report.boundary_work_on_gas_j).abs() < 2e-7);
        assert!(report.energy_residual_j.abs() < 2e-8);
        assert!(report.boundary_work_on_gas_j > 0.0);
        assert_eq!(report.substeps, substeps);
        assert!(report.chemical_accepted_steps >= substeps);
        assert_eq!(
            parcel.model().transport_materials(),
            bank.transport_materials()
        );
        assert!(
            parcel.state().state().dynamic_viscosity > before.state().state().dynamic_viscosity
        );
        assert!(parcel.state().state().sound_speed > before.state().state().sound_speed);
        parcel.state().state().validate_coherence().unwrap();
    }
    for pair in errors.windows(2) {
        assert!(
            pair[0] / pair[1] > 3.8,
            "second-order splitting: {errors:?}"
        );
    }
    assert!(errors[3] < 0.04, "analytic endpoint: {errors:?}");
}

#[test]
fn fixed_volume_is_the_original_rigid_reactor_and_zero_time_is_exact() {
    let (bank, mechanism, mut parcel) = fixture(-24000.0, 2000.0);
    let mut reference = parcel.clone();
    let solve = compression_controls(12);
    let native = advance_closed_parcel_adiabatic(
        &mut reference,
        &bank,
        &mechanism,
        Time::new(0.05),
        solve.kinetics,
        || false,
    )
    .unwrap();
    let report = advance_closed_parcel_compression(
        &mut parcel,
        &bank,
        &mechanism,
        Volume::new(0.01),
        Time::new(0.05),
        solve,
        || false,
    )
    .unwrap();
    assert_eq!(parcel, reference);
    assert_eq!(report.boundary_work_on_gas_j, 0.0);
    assert_eq!(report.energy_residual_j, native.balance.energy_residual_j);
    assert_eq!(report.substeps, 1);
    let before = parcel.clone();
    let zero = advance_closed_parcel_compression(
        &mut parcel,
        &bank,
        &mechanism,
        Volume::new(0.01),
        Time::new(0.0),
        solve,
        || false,
    )
    .unwrap();
    assert_eq!(parcel, before);
    assert_eq!(zero.boundary_work_on_gas_j, 0.0);
    assert_eq!(zero.chemical_accepted_steps, 0);
}

#[test]
fn late_cancellation_and_intermediate_source_failure_are_atomic() {
    let (bank, mechanism, original) = fixture(-24000.0, 2000.0);
    let solve = compression_controls(4);
    let mut probe = original.clone();
    let mut total_polls = 0;
    advance_closed_parcel_compression(
        &mut probe,
        &bank,
        &mechanism,
        Volume::new(0.005),
        Time::new(0.1),
        solve,
        || {
            total_polls += 1;
            false
        },
    )
    .unwrap();
    let mut parcel = original.clone();
    let mut polls = 0;
    let error = advance_closed_parcel_compression(
        &mut parcel,
        &bank,
        &mechanism,
        Volume::new(0.005),
        Time::new(0.1),
        solve,
        || {
            polls += 1;
            polls == total_polls
        },
    )
    .unwrap_err();
    assert_eq!(error, MaterialCompressionError::Cancelled);
    assert_eq!(parcel, original);
    let (restricted, mechanism, mut parcel) = fixture(-24000.0, 450.0);
    let before = parcel.clone();
    assert!(
        advance_closed_parcel_compression(
            &mut parcel,
            &restricted,
            &mechanism,
            Volume::new(0.005),
            Time::new(0.1),
            solve,
            || false,
        )
        .is_err()
    );
    assert_eq!(parcel, before);
}

#[test]
fn invalid_strokes_and_work_budgets_do_not_publish() {
    let (bank, mechanism, original) = fixture(0.0, 2000.0);
    for (volume, duration, steps) in [
        (0.005, 0.0, 4),
        (0.005, 0.1, 0),
        (0.005, 0.1, MAX_COMPRESSION_SUBSTEPS + 1),
        (0.0, 0.1, 4),
        (0.005, -0.1, 4),
    ] {
        let mut parcel = original.clone();
        assert!(matches!(
            advance_closed_parcel_compression(
                &mut parcel,
                &bank,
                &mechanism,
                Volume::new(volume),
                Time::new(duration),
                compression_controls(steps),
                || false,
            ),
            Err(MaterialCompressionError::Invalid(_))
        ));
        assert_eq!(parcel, original);
    }
}
