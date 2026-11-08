//! G1/G3 physical harmonic measurements: analytic amplitudes, phases and ratios.

use fs_fft::C64;
use fs_psycho::harmonic::{
    FundamentalStatus, HarmonicAnalysisPlan, HarmonicFitBudget, HarmonicFitError, HarmonicLevelRatio,
};

fn budget() -> HarmonicFitBudget {
    HarmonicFitBudget {
        max_workspace_bytes: 4 * 1024 * 1024,
        min_reciprocal_condition: 1e-10,
    }
}

fn signal(
    count: usize,
    timestep: f64,
    fundamental: f64,
    orders: &[usize],
    phasors: &[C64],
    dc: f64,
    start: f64,
) -> Vec<f64> {
    (0..count)
        .map(|j| {
            let t = start + j as f64 * timestep;
            dc + orders
                .iter()
                .zip(phasors)
                .map(|(&order, phasor)| {
                    let (sine, cosine) =
                        (std::f64::consts::TAU * fundamental * order as f64 * t).sin_cos();
                    phasor.re * cosine + phasor.im * sine
                })
                .sum::<f64>()
        })
        .collect()
}

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() < tolerance,
        "actual {actual:e}, expected {expected:e}, tolerance {tolerance:e}"
    );
}

fn db(level: HarmonicLevelRatio) -> f64 {
    let HarmonicLevelRatio::Decibels(value) = level else {
        panic!("expected finite level, got {level:?}")
    };
    value
}

#[test]
fn partial_period_record_recovers_physical_harmonics_and_even_deficits() {
    let count = 997;
    let timestep = 1.0 / 48_000.0;
    let fundamental = 423.37;
    let orders = [1, 2, 3, 4, 5, 6, 7];
    let phasors = [
        C64::new(3.0, 4.0),
        C64::new(0.3, 0.4),
        C64::new(-0.6, 0.8),
        C64::new(-0.24, 0.32),
        C64::new(0.0, 0.2),
        C64::new(0.15, -0.2),
        C64::new(0.1, 0.0),
    ];
    let amplitudes = [5.0_f64, 0.5, 1.0, 0.4, 0.2, 0.25, 0.1];
    let plan =
        HarmonicAnalysisPlan::new(count, timestep, fundamental, &orders, 1e-10, budget()).unwrap();
    let samples = signal(count, timestep, fundamental, &orders, &phasors, 7.25, 0.0);
    let mut scratch = vec![0.0; count];
    let mut fitted = [C64::default(); 7];
    let analysis = plan
        .analyze_into(&samples, &mut scratch, &mut fitted)
        .unwrap();

    assert_eq!(analysis.harmonics(), orders);
    assert_eq!(analysis.fundamental_status(), FundamentalStatus::Resolved);
    close(analysis.dc(), 7.25, 2e-12);
    assert!(analysis.residual_rms() < 2e-12);
    let arithmetic_mean = samples.iter().sum::<f64>() / count as f64;
    assert!((arithmetic_mean - analysis.dc()).abs() > 1e-3);
    for (index, &order) in orders.iter().enumerate() {
        let component = analysis.component(order).unwrap();
        assert_eq!(component.order, order);
        assert_eq!(component.frequency_hz, fundamental * order as f64);
        close(component.peak_phasor.re, phasors[index].re, 2e-12);
        close(component.peak_phasor.im, phasors[index].im, 2e-12);
        close(
            db(component.relative_level.unwrap()),
            20.0 * (amplitudes[index] / 5.0).log10(),
            2e-10,
        );
    }
    assert!(analysis.component(8).is_none());
    let distortion = analysis.distortion().unwrap().unwrap();
    let even_square = 0.5_f64.powi(2) + 0.4_f64.powi(2) + 0.25_f64.powi(2);
    let odd_square = 1.0 + 0.2_f64.powi(2) + 0.1_f64.powi(2);
    close(
        distortion.total_ratio.unwrap(),
        (even_square + odd_square).sqrt() / 5.0,
        2e-13,
    );
    close(
        distortion.even_ratio.unwrap(),
        even_square.sqrt() / 5.0,
        2e-13,
    );
    close(
        distortion.odd_ratio.unwrap(),
        odd_square.sqrt() / 5.0,
        2e-13,
    );
    assert!(distortion.complete_band);
    assert_eq!(distortion.highest_order, 7);
    for order in [2, 4, 6] {
        let expected = 20.0 * amplitudes[order - 1].log10()
            - 10.0 * (amplitudes[order - 2] * amplitudes[order]).log10();
        close(
            db(analysis.even_harmonic_deficit(order).unwrap()),
            expected,
            2e-10,
        );
    }
    assert!(analysis.even_harmonic_deficit(0).is_none());
    assert!(analysis.even_harmonic_deficit(3).is_none());
    assert!(analysis.even_harmonic_deficit(8).is_none());
}

#[test]
fn changing_gain_and_phase_origin_preserves_harmonic_ratios() {
    let count = 513;
    let timestep = 1.0 / 8192.0;
    let fundamental = 147.31;
    let orders = [1, 2, 3];
    let original = [
        C64::new(0.6, -0.8),
        C64::new(-0.3, 0.4),
        C64::new(0.24, 0.18),
    ];
    let plan =
        HarmonicAnalysisPlan::new(count, timestep, fundamental, &orders, 0.0, budget()).unwrap();
    let mut scratch = vec![0.0; count];
    let mut fitted = [C64::default(); 3];
    let scratch_ptr = scratch.as_ptr();
    let fitted_ptr = fitted.as_ptr();
    for scale in [1e-200, 1.0, 1e200] {
        for start in [0.0, 0.00137] {
            let authored: Vec<_> = original
                .iter()
                .map(|p| C64::new(scale * p.re, scale * p.im))
                .collect();
            let samples = signal(
                count,
                timestep,
                fundamental,
                &orders,
                &authored,
                0.7 * scale,
                start,
            );
            let analysis = plan
                .analyze_into(&samples, &mut scratch, &mut fitted)
                .unwrap();
            let distortion = analysis.distortion().unwrap().unwrap();
            close(
                distortion.total_ratio.unwrap(),
                (0.25_f64 + 0.09).sqrt(),
                2e-13,
            );
            close(distortion.even_ratio.unwrap(), 0.5, 2e-13);
            close(distortion.odd_ratio.unwrap(), 0.3, 2e-13);
            close(
                db(analysis.even_harmonic_deficit(2).unwrap()),
                20.0 * 0.5_f64.log10() - 10.0 * 0.3_f64.log10(),
                2e-10,
            );
            close(analysis.dc() / scale, 0.7, 2e-13);
            assert!(analysis.residual_rms() / scale < 2e-12);
            for (index, &order) in orders.iter().enumerate() {
                let (sine, cosine) =
                    (std::f64::consts::TAU * fundamental * order as f64 * start).sin_cos();
                let fitted = analysis.peak_phasors()[index];
                close(
                    fitted.re / scale,
                    original[index].re * cosine + original[index].im * sine,
                    2e-12,
                );
                close(
                    fitted.im / scale,
                    original[index].im * cosine - original[index].re * sine,
                    2e-12,
                );
            }
        }
    }
    assert_eq!(scratch.as_ptr(), scratch_ptr);
    assert_eq!(fitted.as_ptr(), fitted_ptr);
}

#[test]
fn unmodeled_inharmonic_energy_is_residual_not_harmonic_distortion() {
    // Integer DFT bins make the unmodeled tone orthogonal to the harmonic
    // columns. It is not an integer multiple of the physical fundamental.
    let count = 1024;
    let timestep = 1.0 / 8192.0;
    let plan =
        HarmonicAnalysisPlan::new(count, timestep, 128.0, &[1, 2], 1e-9, budget()).unwrap();
    let mut samples = signal(
        count,
        timestep,
        128.0,
        &[1, 2],
        &[C64::new(1.0, 0.0), C64::new(0.0, 0.1)],
        3.0,
        0.0,
    );
    for (j, value) in samples.iter_mut().enumerate() {
        *value += 0.7 * (std::f64::consts::TAU * 312.0 * j as f64 * timestep).cos();
    }
    let mut scratch = vec![0.0; count];
    let mut phasors = [C64::default(); 2];
    let analysis = plan
        .analyze_into(&samples, &mut scratch, &mut phasors)
        .unwrap();
    let distortion = analysis.distortion().unwrap().unwrap();
    close(distortion.total_ratio.unwrap(), 0.1, 3e-13);
    close(distortion.even_ratio.unwrap(), 0.1, 3e-13);
    assert_eq!(distortion.odd_ratio, None);
    close(analysis.residual_rms(), 0.7 / 2.0_f64.sqrt(), 3e-13);
    close(analysis.dc(), 3.0, 3e-13);
}

#[test]
fn absent_orders_zero_reference_and_declared_resolution_are_distinct() {
    let count = 257;
    let timestep = 1.0 / 8192.0;
    let mut scratch = vec![0.0; count];
    let mut phasors = [C64::default(); 3];
    let no_fundamental =
        HarmonicAnalysisPlan::new(count, timestep, 193.7, &[2, 3, 4], 1e-8, budget()).unwrap();
    let samples = signal(
        count,
        timestep,
        193.7,
        &[2, 3, 4],
        &[C64::new(0.4, 0.0); 3],
        0.0,
        0.0,
    );
    let analysis = no_fundamental
        .analyze_into(&samples, &mut scratch, &mut phasors)
        .unwrap();
    assert_eq!(analysis.fundamental_status(), FundamentalStatus::NotRetained);
    assert_eq!(analysis.distortion().unwrap(), None);
    assert!(analysis.component(1).is_none());
    assert_eq!(analysis.component(2).unwrap().relative_level, None);
    assert_eq!(analysis.even_harmonic_deficit(2), None);

    let sparse =
        HarmonicAnalysisPlan::new(count, timestep, 193.7, &[1, 3, 5], 1e-8, budget()).unwrap();
    let samples = signal(
        count,
        timestep,
        193.7,
        &[1, 3, 5],
        &[C64::new(0.4, 0.0); 3],
        0.0,
        0.0,
    );
    let analysis = sparse
        .analyze_into(&samples, &mut scratch, &mut phasors)
        .unwrap();
    let distortion = analysis.distortion().unwrap().unwrap();
    assert!(!distortion.complete_band);
    assert_eq!(distortion.highest_order, 5);
    assert_eq!(distortion.even_ratio, None);
    close(distortion.odd_ratio.unwrap(), 2.0_f64.sqrt(), 1e-12);
    assert_eq!(analysis.even_harmonic_deficit(2), None);

    let dense =
        HarmonicAnalysisPlan::new(count, timestep, 193.7, &[1, 2, 3], 1e-4, budget()).unwrap();
    let samples = signal(
        count,
        timestep,
        193.7,
        &[1, 2, 3],
        &[
            C64::new(1e-7, 0.0),
            C64::new(0.4, 0.0),
            C64::new(0.2, 0.0),
        ],
        1.0,
        0.0,
    );
    let analysis = dense
        .analyze_into(&samples, &mut scratch, &mut phasors)
        .unwrap();
    assert_eq!(analysis.fundamental_status(), FundamentalStatus::BelowMinimum);
    assert_eq!(analysis.distortion().unwrap(), None);
    assert_eq!(analysis.component(2).unwrap().relative_level, None);
    assert!(analysis.component(1).unwrap().peak_phasor.re > 0.0);
    assert!(matches!(
        analysis.even_harmonic_deficit(2),
        Some(HarmonicLevelRatio::Decibels(_))
    ));

    let analysis = dense
        .analyze_into(&vec![7.0; count], &mut scratch, &mut phasors)
        .unwrap();
    assert_eq!(analysis.fundamental_status(), FundamentalStatus::Zero);
    assert_eq!(analysis.distortion().unwrap(), None);
    assert_eq!(
        analysis.even_harmonic_deficit(2),
        Some(HarmonicLevelRatio::Undefined)
    );
    assert_eq!(analysis.dc(), 7.0);
    assert_eq!(analysis.residual_rms(), 0.0);

    let single =
        HarmonicAnalysisPlan::new(count, timestep, 193.7, &[1], 1e-8, budget()).unwrap();
    let samples = signal(
        count,
        timestep,
        193.7,
        &[1],
        &[C64::new(0.4, 0.0)],
        0.0,
        0.0,
    );
    let analysis = single
        .analyze_into(&samples, &mut scratch, &mut phasors[..1])
        .unwrap();
    let distortion = analysis.distortion().unwrap().unwrap();
    assert_eq!(distortion.total_ratio, None);
    assert_eq!(distortion.even_ratio, None);
    assert_eq!(distortion.odd_ratio, None);
}

#[test]
fn analysis_keeps_the_shared_clock_budget_and_transactional_fit_refusals() {
    for floor in [-1.0, f64::NAN, f64::INFINITY] {
        assert!(matches!(
            HarmonicAnalysisPlan::new(257, 1.0 / 8192.0, 193.7, &[1, 2, 3], floor, budget()),
            Err(HarmonicFitError::InvalidParameter { .. })
        ));
    }
    assert!(
        HarmonicAnalysisPlan::new(257, 1.0 / 8192.0, 1024.0, &[1, 4], 0.0, budget()).is_err()
    );
    let mut too_small = budget();
    too_small.max_workspace_bytes = 1;
    assert!(matches!(
        HarmonicAnalysisPlan::new(257, 1.0 / 8192.0, 193.7, &[1, 2, 3], 0.0, too_small),
        Err(HarmonicFitError::BudgetExceeded { .. })
    ));
    let plan =
        HarmonicAnalysisPlan::new(257, 1.0 / 8192.0, 193.7, &[1, 2, 3], 1e-8, budget()).unwrap();
    assert_eq!(plan.fit_plan().sample_count(), 257);
    assert_eq!(plan.minimum_fundamental_peak(), 1e-8);
    let mut scratch = vec![0.0; 257];
    let mut phasors = [C64::new(7.0, 8.0); 3];
    let before = phasors;
    assert!(
        plan.analyze_into(&vec![f64::NAN; 257], &mut scratch, &mut phasors)
            .is_err()
    );
    assert_eq!(phasors, before);
    assert!(
        plan.analyze_into(&vec![0.0; 256], &mut scratch, &mut phasors)
            .is_err()
    );
    assert_eq!(phasors, before);
}
