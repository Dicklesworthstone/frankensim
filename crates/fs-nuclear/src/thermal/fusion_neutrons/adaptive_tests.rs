//! G1: analytic absorption followed by source-off heating and material cooling.
use super::*;

#[test]
fn adaptive_fusion_capture_resolves_material_cooling_and_shutoff() {
    let mut receiver = cell(20, 1, 2.0, 1400.0, 0.0);
    receiver.boundary.htc_w_m2_k = 4.0;
    let initial = thermal(vec![receiver], &[(MaterialId(1), &HOT)], vec![]);
    let state = capture_state();
    let pulse = pulse(&state, 0.0, 1.0, 10.0);
    let heating = bound(&initial, Arc::clone(&state), Arc::clone(&pulse));
    // C=8 J/K, G=4 W/K. During the pulse, N'=s-.7N and
    // Pheat=3*.7N. Afterwards the resident neutrons continue to absorb.
    let count = pulse.source().emission().total_particle_count();
    let lambda = 0.7_f64;
    let cooling = 0.5_f64;
    let convolution = |t: f64| ((-lambda * t).exp() - (-cooling * t).exp()) / (cooling - lambda);
    let at_shutoff = 50.0 * (-cooling).exp()
        + 3.0 * count / 8.0 * (-(-cooling).exp_m1() / cooling - convolution(1.0));
    let exact_temperature = 300.0
        + at_shutoff * (-cooling).exp()
        + 3.0 * count * -(-lambda).exp_m1() / 8.0 * convolution(1.0);
    let exact_neutrons = count / lambda * -(-lambda).exp_m1() * (-lambda).exp();
    let mut errors = Vec::new();
    for temperature in [1e-3, 1e-5] {
        let controls = NuclearAdaptiveControls {
            step: StepConfig {
                max_dt_s: 1.0,
                energy_tolerance_j: 1e-6,
            },
            initial_dt_s: 1.0,
            minimum_dt_s: 1e-8,
            maximum_trials: 20_000,
            maximum_substeps: 20_000,
            tolerances: NuclearTimeTolerances {
                temperature_k: temperature,
                specific_energy_j_kg: 4.0 * temperature,
                nuclei: 1e-6,
                nuclear_energy_j: 1e-6,
                neutron_population: 1e-6,
                neutron_count: 1e-6,
                relative: 0.0,
            },
        };
        let (material, evolved, report) = heating
            .advance_adaptive_to(&initial, 2.0, controls, || Ok(()))
            .unwrap();
        errors.push((material.temperature_k(CellId(20)).unwrap() - exact_temperature).abs());
        assert!(report.accuracy_rejections > 0);
        assert!(report.maximum_accepted_normalized_error <= 1.0);
        assert_eq!(evolved.audit().substeps, 2 * report.accepted_comparisons);
        let source = evolved.source(CellId(20)).unwrap();
        close(
            source.state().neutrons().neutron_counts()[0],
            exact_neutrons,
            1e-10,
        );
        close(
            evolved.cumulative_released().unwrap().deposited_j,
            3.0 * (count - exact_neutrons),
            1e-9,
        );
        assert_eq!(
            source
                .forcing()
                .pulse()
                .unwrap()
                .rates_at(Time::new(2.0))
                .unwrap(),
            &[0.0]
        );
        assert!(
            evolved
                .last_step(CellId(20))
                .unwrap()
                .spatial()
                .released
                .deposited_j
                > 0.0
        );
        assert!(report.thermal.residual_j.abs() < controls.step.energy_tolerance_j);
    }
    assert!(errors[1] < errors[0] / 5.0, "temperature errors {errors:?}");
    assert!(errors[1] < 0.02, "temperature errors {errors:?}");
    assert_eq!(initial.time_s(), 0.0);
    assert_eq!(state.neutrons().neutron_counts(), &[0.0]);
}
