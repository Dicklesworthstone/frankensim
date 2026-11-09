//! G1/G3: coherent tensor-Prony work and heat against independent real ODE fields.

use fs_fft::harmonic_window::{HarmonicWindow, HarmonicWindowProfile};
use fs_material::visco::{
    MAX_HARMONIC_MATERIAL_LINES, Mandel6, PronyBranch, TensorProny, TensorPronyHarmonicWindowPoint,
    ViscoError,
};
use fs_math::c64::C64;

fn outer(diagonal: f64, a: [f64; 6]) -> Mandel6 {
    core::array::from_fn(|i| {
        core::array::from_fn(|j| a[i] * a[j] + if i == j { diagonal } else { 0.0 })
    })
}
fn law() -> TensorProny {
    TensorProny::new_mandel(
        outer(8.0, [1.0, -0.5, 0.2, 0.4, -0.3, 0.6]),
        vec![
            PronyBranch::new(outer(0.0, [1.0, 2.0, -0.5, 0.8, -0.4, 0.3]), 0.2).unwrap(),
            PronyBranch::new(outer(0.4, [-0.2, 0.3, 1.0, -0.6, 0.7, 0.8]), 1.4).unwrap(),
        ],
        0.05,
    )
    .unwrap()
    .with_kelvin_viscosity_mandel(outer(0.01, [0.1, -0.2, 0.05, 0.08, -0.03, 0.04]))
    .unwrap()
    .with_band((0.1, 100.0))
    .unwrap()
}
fn strain() -> [C64; 6] {
    core::array::from_fn(|i| C64::new((i as f64 - 2.0) * 0.001, (3.0 - i as f64) * 0.0007))
}
fn apply(matrix: &Mandel6, vector: &[f64; 6]) -> [f64; 6] {
    core::array::from_fn(|i| (0..6).map(|j| matrix[i][j] * vector[j]).sum())
}
fn dot(a: &[f64; 6], b: &[f64; 6]) -> f64 {
    a.iter().zip(b).map(|(x, y)| x * y).sum()
}
fn near(actual: f64, expected: f64, scale: f64, relative: f64) {
    assert!(
        (actual - expected).abs() <= relative * scale,
        "actual={actual:e}, expected={expected:e}, scale={scale:e}"
    );
}

// Solve the actual Maxwell ODE in real cosine/sine coordinates. Sum the
// physical fields before squaring. No production spectral/Gram utility enters.
fn instantaneous(law: &TensorProny, frequencies: &[f64], strains: &[[C64; 6]], t: f64) -> [f64; 4] {
    instantaneous_with_bias(law, frequencies, strains, &[0.0; 6], t)
}
fn instantaneous_with_bias(
    law: &TensorProny,
    frequencies: &[f64],
    strains: &[[C64; 6]],
    bias: &[f64; 6],
    t: f64,
) -> [f64; 4] {
    let (mut eps, mut rate) = (*bias, [0.0; 6]);
    for (&w, e) in frequencies.iter().zip(strains) {
        let (s, c) = (w * t).sin_cos();
        for i in 0..6 {
            eps[i] += e[i].re * c + e[i].im * s;
            rate[i] += w * (e[i].im * c - e[i].re * s);
        }
    }
    let mut stress = apply(law.c_eq(), &eps);
    let mut storage = 0.5 * dot(&eps, &stress);
    let mut kelvin = 0.0;
    if let Some(viscosity) = law.kv_pa_s() {
        let viscous = apply(viscosity, &rate);
        kelvin = dot(&rate, &viscous);
        for i in 0..6 {
            stress[i] += viscous[i];
        }
    }
    let mut maxwell = 0.0;
    for branch in law.branches() {
        let (mut q, mut q_rate) = ([0.0; 6], [0.0; 6]);
        for (&w, e) in frequencies.iter().zip(strains) {
            let (s, c) = (w * t).sin_cos();
            let x = w * branch.tau_s;
            let denominator = 1.0 + x * x;
            for i in 0..6 {
                let a = (x * e[i].im + x * x * e[i].re) / denominator;
                let b = (x * x * e[i].im - x * e[i].re) / denominator;
                q[i] += a * c + b * s;
                q_rate[i] += w * (b * c - a * s);
            }
        }
        for i in 0..6 {
            near(
                q_rate[i] + q[i] / branch.tau_s,
                rate[i],
                rate[i].abs().max(0.1),
                3e-14,
            );
        }
        let branch_stress = apply(&branch.c_pa, &q);
        let square = dot(&q, &branch_stress);
        storage += 0.5 * square;
        maxwell += square / branch.tau_s;
        for i in 0..6 {
            stress[i] += branch_stress[i];
        }
    }
    [storage, kelvin, maxwell, dot(&stress, &rate)]
}
fn integrate(
    law: &TensorProny,
    frequencies: &[f64],
    strains: &[[C64; 6]],
    window: HarmonicWindow,
) -> [f64; 3] {
    integrate_with_bias(law, frequencies, strains, &[0.0; 6], window)
}
fn integrate_with_bias(
    law: &TensorProny,
    frequencies: &[f64],
    strains: &[[C64; 6]],
    bias: &[f64; 6],
    window: HarmonicWindow,
) -> [f64; 3] {
    let panels = 8192;
    let mut result = [0.0; 3];
    for i in 0..=panels {
        let t = window.start_s() + window.duration_s() * i as f64 / panels as f64;
        let value = instantaneous_with_bias(law, frequencies, strains, bias, t);
        let weight = if i == 0 || i == panels {
            1.0
        } else if i % 2 == 1 {
            4.0
        } else {
            2.0
        };
        for k in 0..3 {
            result[k] += weight * value[k + 1];
        }
    }
    result.map(|v| v * window.duration_s() / (3.0 * panels as f64))
}
fn quantities(p: &TensorPronyHarmonicWindowPoint) -> [f64; 8] {
    [
        p.start_stored_energy_j_m3,
        p.end_stored_energy_j_m3,
        p.stored_energy_change_j_m3,
        p.kelvin_heat_j_m3,
        p.maxwell_heat_j_m3,
        p.material_heat_j_m3,
        p.stress_work_j_m3,
        p.energy_balance_residual_j_m3,
    ]
}

#[test]
fn multitone_tensor_heat_and_work_match_independent_real_field_integration() {
    let law = law();
    let frequencies = [2.7, 3.1, 11.0];
    let strains = [
        strain(),
        strain().map(|z| z * C64::new(0.6, -0.4)),
        strain().map(|z| z * C64::new(-0.2, 0.1)),
    ];
    let window = HarmonicWindow::new(0.37, 1.9).unwrap();
    let spectrum = law.prepare_harmonic_spectrum(&frequencies).unwrap();
    let result = spectrum
        .prepare_window(window)
        .unwrap()
        .evaluate(&strains)
        .unwrap();
    let expected = integrate(&law, &frequencies, &strains, window);
    let start = instantaneous(&law, &frequencies, &strains, window.start_s())[0];
    let end = instantaneous(
        &law,
        &frequencies,
        &strains,
        window.start_s() + window.duration_s(),
    )[0];
    let scale = start.max(end).max(expected[2].abs());
    near(result.start_stored_energy_j_m3, start, scale, 3e-13);
    near(result.end_stored_energy_j_m3, end, scale, 3e-13);
    near(result.stored_energy_change_j_m3, end - start, scale, 3e-13);
    near(result.kelvin_heat_j_m3, expected[0], scale, 2e-11);
    near(result.maxwell_heat_j_m3, expected[1], scale, 2e-11);
    near(result.stress_work_j_m3, expected[2], scale, 2e-11);
    near(result.energy_balance_residual_j_m3, 0.0, scale, 4e-13);
}

#[test]
fn finite_beat_heating_is_not_the_sum_of_separate_tone_heat() {
    let law = law();
    let frequencies = [10.0, 10.2];
    let strains = [strain(); 2];
    let window = HarmonicWindow::new(0.0, 0.7).unwrap();
    let spectrum = law.prepare_harmonic_spectrum(&frequencies).unwrap();
    let coherent = spectrum
        .prepare_window(window)
        .unwrap()
        .evaluate(&strains)
        .unwrap();
    let diagonal: f64 = frequencies
        .iter()
        .map(|&w| {
            law.prepare_harmonic_material(w)
                .unwrap()
                .observe_window(&strain(), window)
                .unwrap()
                .material_heat_j_m3
        })
        .sum();
    assert!(coherent.material_heat_j_m3 > 1.9 * diagonal);
    let expected = integrate(&law, &frequencies, &strains, window);
    near(
        coherent.material_heat_j_m3,
        expected[0] + expected[1],
        coherent.material_heat_j_m3,
        2e-12,
    );
}

#[test]
fn equal_frequency_cancellation_precedes_source_admission_and_physical_products() {
    let law = law();
    let mut tiny = [C64::ZERO; 6];
    tiny[0] = C64::new(0.001, -0.0007);
    let mut large = [C64::ZERO; 6];
    large[0] = C64::new(1e16, -2e16);
    let strains = [large, tiny, large.map(|z| -z)];
    let spectrum = law.prepare_harmonic_spectrum(&[2.7; 3]).unwrap();
    spectrum.validate_strains(&strains).unwrap();
    let window = HarmonicWindow::new(0.31, 2.1).unwrap();
    let result = spectrum
        .prepare_window(window)
        .unwrap()
        .evaluate(&strains)
        .unwrap();
    let reference = law
        .prepare_harmonic_material(2.7)
        .unwrap()
        .observe_window(&tiny, window)
        .unwrap();
    for (actual, expected) in quantities(&result).into_iter().zip(quantities(&reference)) {
        near(
            actual,
            expected,
            reference
                .material_heat_j_m3
                .max(reference.start_stored_energy_j_m3),
            2e-13,
        );
    }
    let cancelled = spectrum
        .prepare_window(window)
        .unwrap()
        .evaluate(&[large, [C64::ZERO; 6], large.map(|z| -z)])
        .unwrap();
    assert_eq!(quantities(&cancelled), [0.0; 8]);
}

#[test]
fn one_tone_matches_existing_original_material_window_at_every_quantity() {
    let law = law();
    let spectrum = law.prepare_harmonic_spectrum(&[2.7]).unwrap();
    let single = law.prepare_harmonic_material(2.7).unwrap();
    for (start, duration) in [(0.0, 2.0), (0.4, 1e-14), (1e18, 0.3)] {
        let window = HarmonicWindow::new(start, duration).unwrap();
        let actual = spectrum
            .prepare_window(window)
            .unwrap()
            .evaluate(&[strain()])
            .unwrap();
        let expected = single.observe_window(&strain(), window).unwrap();
        for (a, e) in quantities(&actual).into_iter().zip(quantities(&expected)) {
            near(
                a,
                e,
                expected
                    .start_stored_energy_j_m3
                    .max(expected.material_heat_j_m3),
                2e-13,
            );
        }
    }
}

#[test]
fn a_source_clock_shift_and_coherent_phases_preserve_every_energy() {
    let law = law();
    let frequencies = [2.7, 3.1, 11.0];
    let spectrum = law.prepare_harmonic_spectrum(&frequencies).unwrap();
    let strains = [
        strain(),
        strain().map(|z| z * C64::new(-0.4, 0.7)),
        strain().map(|z| z.scale(0.3)),
    ];
    let shift = 1.4;
    let shifted = core::array::from_fn::<_, 3, _>(|i| {
        let (s, c) = (-frequencies[i] * shift).sin_cos();
        strains[i].map(|z| z * C64::new(c, s))
    });
    let original = spectrum
        .prepare_window(HarmonicWindow::new(2.3, 0.9).unwrap())
        .unwrap()
        .evaluate(&strains)
        .unwrap();
    let translated = spectrum
        .prepare_window(HarmonicWindow::new(2.3 - shift, 0.9).unwrap())
        .unwrap()
        .evaluate(&shifted)
        .unwrap();
    for (a, e) in quantities(&translated)
        .into_iter()
        .zip(quantities(&original))
    {
        near(
            a,
            e,
            original
                .material_heat_j_m3
                .max(original.start_stored_energy_j_m3),
            5e-13,
        );
    }
}

#[test]
fn complete_waveform_domain_and_original_frequency_domain_are_preserved() {
    let law = law();
    let mut e = [C64::ZERO; 6];
    e[0] = C64::from_re(0.04);
    let spectrum = law.prepare_harmonic_spectrum(&[2.7, 3.1]).unwrap();
    assert!(spectrum.validate_strains(&[e, e]).is_err());
    assert!(spectrum.validate_strains(&[e]).is_err());
    let mut invalid = e;
    invalid[2].im = f64::NAN;
    assert!(spectrum.validate_strains(&[e, invalid]).is_err());
    for frequency in [0.0, -1.0, 0.05, 101.0, f64::NAN, f64::INFINITY] {
        assert!(law.prepare_harmonic_spectrum(&[frequency]).is_err());
    }
    assert!(
        law.prepare_harmonic_spectrum(&vec![2.7; MAX_HARMONIC_MATERIAL_LINES + 1])
            .is_err()
    );
    let shaped = HarmonicWindow::with_profile(0.0, 1.0, HarmonicWindowProfile::Triangular).unwrap();
    assert!(spectrum.prepare_window(shaped).is_err());
    let empty = law.prepare_harmonic_spectrum(&[]).unwrap();
    let zero = empty
        .prepare_window(HarmonicWindow::new(0.0, 1.0).unwrap())
        .unwrap()
        .evaluate(&[])
        .unwrap();
    assert_eq!(quantities(&zero), [0.0; 8]);
}

#[test]
fn material_channel_checkpoint_preserves_the_callers_cancellation_and_reuse() {
    #[derive(Debug)]
    enum Error {
        Material,
        Cancelled,
    }
    impl From<ViscoError> for Error {
        fn from(_: ViscoError) -> Self {
            Self::Material
        }
    }
    let law = law();
    let spectrum = law.prepare_harmonic_spectrum(&[2.7, 3.1]).unwrap();
    let prepared = spectrum
        .prepare_window(HarmonicWindow::new(0.3, 1.7).unwrap())
        .unwrap();
    let mut checkpoints = 0;
    let result = prepared.evaluate_with_checkpoint(&[strain(); 2], || {
        checkpoints += 1;
        if checkpoints == 4 {
            Err(Error::Cancelled)
        } else {
            Ok(())
        }
    });
    assert!(matches!(result, Err(Error::Cancelled)));
    let a = prepared.evaluate(&[strain(); 2]).unwrap();
    let b = prepared.evaluate(&[strain(); 2]).unwrap();
    assert_eq!(a, b);
}

#[test]
fn relaxed_bias_keeps_real_storage_and_work_cross_terms_without_extra_heat() {
    let law = law().with_band((0.0, 100.0)).unwrap();
    let frequencies = [2.7, 3.1, 11.0];
    let strains = [
        strain(),
        strain().map(|z| z * C64::new(0.6, -0.4)),
        strain().map(|z| z.scale(0.2)),
    ];
    let bias = [0.006, -0.004, 0.003, -0.002, 0.001, 0.003];
    let window = HarmonicWindow::new(0.37, 1.9).unwrap();
    let spectrum = law.prepare_harmonic_spectrum(&frequencies).unwrap();
    spectrum
        .validate_strains_with_bias(&strains, &bias)
        .unwrap();
    let prepared = spectrum.prepare_window(window).unwrap();
    let point = prepared.evaluate_with_bias(&strains, &bias).unwrap();
    let ac = prepared.evaluate(&strains).unwrap();
    let expected = integrate_with_bias(&law, &frequencies, &strains, &bias, window);
    let start = instantaneous_with_bias(&law, &frequencies, &strains, &bias, window.start_s())[0];
    let end = instantaneous_with_bias(
        &law,
        &frequencies,
        &strains,
        &bias,
        window.start_s() + window.duration_s(),
    )[0];
    let scale = start.max(end).max(point.material_heat_j_m3);
    near(point.start_stored_energy_j_m3, start, scale, 4e-13);
    near(point.end_stored_energy_j_m3, end, scale, 4e-13);
    near(point.stored_energy_change_j_m3, end - start, scale, 4e-13);
    near(point.kelvin_heat_j_m3, expected[0], scale, 2e-11);
    near(point.maxwell_heat_j_m3, expected[1], scale, 2e-11);
    near(point.stress_work_j_m3, expected[2], scale, 2e-11);
    near(point.energy_balance_residual_j_m3, 0.0, scale, 4e-13);
    assert!((point.stress_work_j_m3 - ac.stress_work_j_m3).abs() > 1e-5 * scale);
    assert_eq!(point.kelvin_heat_j_m3, ac.kelvin_heat_j_m3);
    assert_eq!(point.maxwell_heat_j_m3, ac.maxwell_heat_j_m3);
    assert_eq!(point.material_heat_j_m3, ac.material_heat_j_m3);
}

#[test]
fn bias_work_survives_a_window_smaller_than_endpoint_energy_resolution() {
    let law = TensorProny::new_mandel(outer(8.0, [0.0; 6]), vec![], 0.05).unwrap();
    let spectrum = law.prepare_harmonic_spectrum(&[3.3]).unwrap();
    let mut strain = [C64::ZERO; 6];
    strain[0] = C64::new(0.0, 0.002);
    let bias = [0.01, 0.0, 0.0, 0.0, 0.0, 0.0];
    let duration = 1e-20;
    let point = spectrum
        .prepare_window(HarmonicWindow::new(0.0, duration).unwrap())
        .unwrap()
        .evaluate_with_bias(&[strain], &bias)
        .unwrap();
    let work = 8.0 * bias[0] * 3.3 * strain[0].im * duration;
    assert_eq!(point.start_stored_energy_j_m3, point.end_stored_energy_j_m3);
    assert!(point.stored_energy_change_j_m3 > 0.0);
    near(point.stored_energy_change_j_m3, work, work, 4e-14);
    near(point.stress_work_j_m3, work, work, 4e-14);
    assert_eq!(point.material_heat_j_m3, 0.0);
}

#[test]
fn relaxed_dc_only_stores_equilibrium_energy_and_preserves_the_full_domain() {
    let law = law().with_band((0.0, 100.0)).unwrap();
    let bias = [0.006, -0.004, 0.003, -0.002, 0.001, 0.003];
    let spectrum = law.prepare_harmonic_spectrum(&[]).unwrap();
    let point = spectrum
        .prepare_window(HarmonicWindow::new(1e15, 1e-8).unwrap())
        .unwrap()
        .evaluate_with_bias(&[], &bias)
        .unwrap();
    let expected = 0.5 * dot(&bias, &apply(law.c_eq(), &bias));
    near(point.start_stored_energy_j_m3, expected, expected, 4e-14);
    assert_eq!(point.start_stored_energy_j_m3, point.end_stored_energy_j_m3);
    assert_eq!(point.stored_energy_change_j_m3, 0.0);
    assert_eq!(point.stress_work_j_m3, 0.0);
    assert_eq!(point.material_heat_j_m3, 0.0);
    let mut e = [C64::ZERO; 6];
    e[0] = C64::from_re(0.03);
    let bias = [0.03, 0.0, 0.0, 0.0, 0.0, 0.0];
    let ac = law.prepare_harmonic_spectrum(&[3.1]).unwrap();
    ac.validate_strains(&[e]).unwrap();
    assert!(ac.validate_strains_with_bias(&[e], &bias).is_err());
}

#[test]
fn zero_bias_keeps_ac_band_and_cancellation_while_nonzero_bias_requires_dc() {
    let law = law();
    let spectrum = law.prepare_harmonic_spectrum(&[3.1; 3]).unwrap();
    let large = [C64::new(1e200, -1e200); 6];
    let fields = [large, strain(), large.map(|z| -z)];
    let prepared = spectrum
        .prepare_window(HarmonicWindow::new(0.3, 1.2).unwrap())
        .unwrap();
    let ac = prepared.evaluate(&fields).unwrap();
    let zero = prepared.evaluate_with_bias(&fields, &[0.0; 6]).unwrap();
    assert_eq!(quantities(&zero), quantities(&ac));
    let bias = [0.001; 6];
    assert!(matches!(
        prepared.evaluate_with_bias(&fields, &bias),
        Err(ViscoError::OutOfBand { omega: 0.0, .. })
    ));
    assert!(spectrum.validate_strains_with_bias(&fields, &bias).is_err());
    assert!(
        prepared
            .evaluate_with_bias(&fields, &[f64::NAN; 6])
            .is_err()
    );
    let cancelled = prepared.evaluate_with_bias_and_checkpoint(&fields, &[0.0; 6], || {
        Err(ViscoError::Parameters {
            what: "cancelled biased observation",
        })
    });
    assert!(matches!(
        cancelled,
        Err(ViscoError::Parameters {
            what: "cancelled biased observation"
        })
    ));
    assert_eq!(
        quantities(&prepared.evaluate(&fields).unwrap()),
        quantities(&ac)
    );
}
