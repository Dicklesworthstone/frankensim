//! Coherent thin-body energy against original 3D real fields and transverse constraints.

use fs_fft::harmonic_window::{HarmonicWindow, HarmonicWindowProfile};
use fs_material::{
    plane_stress_viscosity::MaterialPlane,
    visco::{
        MAX_HARMONIC_MATERIAL_LINES, Mandel6, PlaneStressHarmonic, PlaneStressHarmonicSpectrum,
        PlaneStressSpectrumWindow, PronyBranch, TensorProny, TensorPronyHarmonicWindowPoint,
        ViscoError,
    },
};
use fs_math::c64::C64;

fn isotropic(k: f64, g: f64) -> Mandel6 {
    core::array::from_fn(|i| {
        core::array::from_fn(|j| {
            if i < 3 && j < 3 {
                k - 2.0 * g / 3.0 + if i == j { 2.0 * g } else { 0.0 }
            } else if i == j {
                2.0 * g
            } else {
                0.0
            }
        })
    })
}
fn law() -> TensorProny {
    TensorProny::new_mandel(
        isotropic(4.0, 3.0),
        vec![
            PronyBranch::new(isotropic(20.0, 0.1), 0.08).unwrap(),
            PronyBranch::new(isotropic(0.1, 4.0), 2.0).unwrap(),
        ],
        0.05,
    )
    .unwrap()
    .with_kelvin_viscosity_mandel(isotropic(0.02, 0.2))
    .unwrap()
    .with_band((0.1, 100.0))
    .unwrap()
}
fn in_plane() -> [C64; 3] {
    [
        C64::new(0.003, 0.001),
        C64::new(-0.001, 0.002),
        C64::new(0.002, -0.001),
    ]
}
fn near(actual: f64, expected: f64, scale: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance * scale,
        "actual={actual:e}, expected={expected:e}, scale={scale:e}"
    );
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
fn compare(actual: &TensorPronyHarmonicWindowPoint, expected: &TensorPronyHarmonicWindowPoint) {
    let scale = expected
        .start_stored_energy_j_m3
        .max(expected.material_heat_j_m3);
    for (a, b) in quantities(actual).into_iter().zip(quantities(expected)) {
        near(a, b, scale, 8e-13);
    }
}

// Independent isotropic closed form: lambda(omega)/(lambda(omega)+2G(omega)).
// Bulk and shear memories have different times, so a fixed real Poisson ratio
// cannot reconstruct these fields. No production constraint solver is used.
fn closed_form_strain(omega: f64, e: [C64; 3]) -> [C64; 6] {
    let (mut k, mut g) = (C64::new(4.0, -0.02 * omega), C64::new(3.0, -0.2 * omega));
    for (bulk, shear, tau) in [(20.0, 0.1, 0.08), (0.1, 4.0, 2.0)] {
        let x = omega * tau;
        let transfer = C64::new(x * x, -x).scale(1.0 / (1.0 + x * x));
        k = k + transfer.scale(bulk);
        g = g + transfer.scale(shear);
    }
    let lambda = k - g.scale(2.0 / 3.0);
    let normal = -(lambda / (lambda + g.scale(2.0))) * (e[0] + e[1]);
    [
        e[0],
        e[1],
        normal,
        e[2].scale(core::f64::consts::FRAC_1_SQRT_2),
        C64::ZERO,
        C64::ZERO,
    ]
}
fn apply(c: &Mandel6, e: &[f64; 6]) -> [f64; 6] {
    core::array::from_fn(|i| (0..6).map(|j| c[i][j] * e[j]).sum())
}
fn dot(a: &[f64; 6], b: &[f64; 6]) -> f64 {
    a.iter().zip(b).map(|(x, y)| x * y).sum()
}
fn real_fields(law: &TensorProny, frequencies: &[f64], strains: &[[C64; 6]], t: f64) -> [f64; 4] {
    let (mut e, mut rate) = ([0.0; 6], [0.0; 6]);
    for (&w, z) in frequencies.iter().zip(strains) {
        let (s, c) = (w * t).sin_cos();
        for i in 0..6 {
            e[i] += z[i].re * c + z[i].im * s;
            rate[i] += w * (z[i].im * c - z[i].re * s);
        }
    }
    let mut stress = apply(law.c_eq(), &e);
    let mut stored = 0.5 * dot(&e, &stress);
    let viscous = apply(law.kv_pa_s().unwrap(), &rate);
    let kelvin = dot(&rate, &viscous);
    for i in 0..6 {
        stress[i] += viscous[i];
    }
    let mut heat = 0.0;
    for branch in law.branches() {
        let mut q = [0.0; 6];
        for (&w, z) in frequencies.iter().zip(strains) {
            let (s, c) = (w * t).sin_cos();
            let x = w * branch.tau_s;
            for i in 0..6 {
                let a = (x * x * z[i].re + x * z[i].im) / (1.0 + x * x);
                let b = (x * x * z[i].im - x * z[i].re) / (1.0 + x * x);
                q[i] += a * c + b * s;
            }
        }
        let branch_stress = apply(&branch.c_pa, &q);
        let square = dot(&q, &branch_stress);
        stored += 0.5 * square;
        heat += square / branch.tau_s;
        for i in 0..6 {
            stress[i] += branch_stress[i];
        }
    }
    for i in [2, 4, 5] {
        near(stress[i], 0.0, 0.1, 3e-14);
    }
    // Only work-conjugate in-plane tractions are supplied to a thin body.
    let work = stress[0] * rate[0] + stress[1] * rate[1] + stress[3] * rate[3];
    [stored, kelvin, heat, work]
}

#[test]
fn dispersive_transverse_relaxation_matches_independent_real_fields() {
    let law = law();
    let frequencies = [1.3, 4.1, 9.7];
    let e = [
        in_plane(),
        in_plane().map(|z| z * C64::new(0.6, -0.3)),
        in_plane().map(|z| z.scale(0.2)),
    ];
    let full = core::array::from_fn::<_, 3, _>(|i| closed_form_strain(frequencies[i], e[i]));
    let ratios = [1.3, 9.7].map(|w| closed_form_strain(w, [C64::ONE, C64::ZERO, C64::ZERO])[2]);
    assert!(
        (ratios[0] - ratios[1]).abs() > 0.05,
        "test must exercise dispersive normal strain"
    );
    let window = HarmonicWindow::new(0.37, 1.7).unwrap();
    let spectrum = law
        .prepare_plane_stress_harmonic_spectrum(MaterialPlane::Axes12, &frequencies)
        .unwrap();
    let point = spectrum
        .prepare_window(window)
        .unwrap()
        .evaluate(&e)
        .unwrap();
    let panels = 4096;
    let mut integral = [0.0; 3];
    for i in 0..=panels {
        let t = window.start_s() + window.duration_s() * i as f64 / panels as f64;
        let value = real_fields(&law, &frequencies, &full, t);
        let weight = if i == 0 || i == panels {
            1.0
        } else if i % 2 == 1 {
            4.0
        } else {
            2.0
        };
        for k in 0..3 {
            integral[k] += weight * value[k + 1];
        }
    }
    let integral = integral.map(|x| x * window.duration_s() / (3.0 * panels as f64));
    let start = real_fields(&law, &frequencies, &full, window.start_s())[0];
    let end = real_fields(
        &law,
        &frequencies,
        &full,
        window.start_s() + window.duration_s(),
    )[0];
    let scale = start.max(end).max(point.material_heat_j_m3);
    near(point.start_stored_energy_j_m3, start, scale, 8e-13);
    near(point.end_stored_energy_j_m3, end, scale, 8e-13);
    near(point.stored_energy_change_j_m3, end - start, scale, 8e-13);
    near(point.kelvin_heat_j_m3, integral[0], scale, 3e-11);
    near(point.maxwell_heat_j_m3, integral[1], scale, 3e-11);
    near(point.stress_work_j_m3, integral[2], scale, 3e-11);
    near(point.energy_balance_residual_j_m3, 0.0, scale, 8e-13);
}

#[test]
fn anisotropic_ordered_planes_retain_original_tensor_energy_and_beats() {
    let outer = |d: f64, a: [f64; 6]| {
        core::array::from_fn(|i| {
            core::array::from_fn(|j| a[i] * a[j] + if i == j { d } else { 0.0 })
        })
    };
    let law = TensorProny::new_mandel(
        outer(8.0, [1.0, -0.5, 0.2, 0.4, -0.3, 0.6]),
        vec![PronyBranch::new(outer(0.2, [1.0, 2.0, -0.5, 0.8, -0.4, 0.3]), 0.2).unwrap()],
        0.05,
    )
    .unwrap();
    let frequencies = [10.0, 10.2];
    let window = HarmonicWindow::new(0.0, 0.7).unwrap();
    for plane in [
        MaterialPlane::Axes12,
        MaterialPlane::Axes21,
        MaterialPlane::Axes13,
        MaterialPlane::Axes23,
    ] {
        let spectrum = law
            .prepare_plane_stress_harmonic_spectrum(plane, &frequencies)
            .unwrap();
        let actual = spectrum
            .prepare_window(window)
            .unwrap()
            .evaluate(&[in_plane(); 2])
            .unwrap();
        let full = frequencies.map(|w| {
            law.prepare_plane_stress_harmonic(plane, w)
                .unwrap()
                .recover_strain(&in_plane())
                .unwrap()
        });
        let expected = law
            .prepare_harmonic_spectrum(&frequencies)
            .unwrap()
            .prepare_window(window)
            .unwrap()
            .evaluate(&full)
            .unwrap();
        compare(&actual, &expected);
        let diagonal: f64 = frequencies
            .iter()
            .zip(full)
            .map(|(&w, e)| {
                law.prepare_harmonic_material(w)
                    .unwrap()
                    .observe_window(&e, window)
                    .unwrap()
                    .material_heat_j_m3
            })
            .sum();
        assert!(actual.material_heat_j_m3 > 1.9 * diagonal);
    }
}

#[test]
fn duplicate_coordinates_cancel_before_transverse_recovery_and_domain_admission() {
    let law = law();
    let e = in_plane();
    let large = [C64::new(1e250, -1e250); 3];
    let spectrum = law
        .prepare_plane_stress_harmonic_spectrum(MaterialPlane::Axes12, &[4.1; 3])
        .unwrap();
    let coordinates = [large, e, large.map(|z| -z)];
    spectrum.validate_strains(&coordinates).unwrap();
    let window = HarmonicWindow::new(0.3, 1.1).unwrap();
    let actual = spectrum
        .prepare_window(window)
        .unwrap()
        .evaluate(&coordinates)
        .unwrap();
    let full = law
        .prepare_plane_stress_harmonic(MaterialPlane::Axes12, 4.1)
        .unwrap()
        .recover_strain(&e)
        .unwrap();
    let expected = law
        .prepare_harmonic_material(4.1)
        .unwrap()
        .observe_window(&full, window)
        .unwrap();
    compare(&actual, &expected);
    let cancelled = law
        .prepare_plane_stress_harmonic_spectrum(MaterialPlane::Axes12, &[4.1; 2])
        .unwrap()
        .prepare_window(window)
        .unwrap()
        .evaluate(&[large, large.map(|z| -z)])
        .unwrap();
    assert_eq!(quantities(&cancelled), [0.0; 8]);
}

#[test]
fn source_clock_and_phase_translation_preserve_thin_body_energy() {
    let law = law();
    let frequencies = [2.7, 3.1, 11.0];
    let spectrum = law
        .prepare_plane_stress_harmonic_spectrum(MaterialPlane::Axes12, &frequencies)
        .unwrap();
    let e = [
        in_plane(),
        in_plane().map(|z| z * C64::new(-0.4, 0.7)),
        in_plane().map(|z| z.scale(0.3)),
    ];
    let shift = 1.4;
    let shifted = core::array::from_fn::<_, 3, _>(|i| {
        let (s, c) = (-frequencies[i] * shift).sin_cos();
        e[i].map(|z| z * C64::new(c, s))
    });
    let original = spectrum
        .prepare_window(HarmonicWindow::new(2.3, 0.9).unwrap())
        .unwrap()
        .evaluate(&e)
        .unwrap();
    let translated = spectrum
        .prepare_window(HarmonicWindow::new(2.3 - shift, 0.9).unwrap())
        .unwrap()
        .evaluate(&shifted)
        .unwrap();
    compare(&translated, &original);
}

#[test]
fn original_full_strain_frequency_window_and_work_bounds_are_retained() {
    let stiff = TensorProny::new_mandel(isotropic(100.0, 1.0), vec![], 0.05).unwrap();
    let spectrum = stiff
        .prepare_plane_stress_harmonic_spectrum(MaterialPlane::Axes12, &[1.0])
        .unwrap();
    let in_plane = [C64::from_re(0.03), C64::from_re(0.03), C64::ZERO];
    assert!(0.03_f64.hypot(0.03) < 0.05);
    assert!(
        spectrum.validate_strains(&[in_plane]).is_err(),
        "recovered normal strain must enter the source domain"
    );
    let law = law();
    for w in [0.0, -1.0, 0.05, 101.0, f64::NAN, f64::INFINITY] {
        assert!(
            law.prepare_plane_stress_harmonic_spectrum(MaterialPlane::Axes12, &[w])
                .is_err()
        );
    }
    assert!(
        law.prepare_plane_stress_harmonic_spectrum(
            MaterialPlane::Axes12,
            &vec![2.7; MAX_HARMONIC_MATERIAL_LINES + 1]
        )
        .is_err()
    );
    let frequencies = [2.7, 3.1, 11.0];
    let spectrum = law
        .prepare_plane_stress_harmonic_spectrum(MaterialPlane::Axes12, &frequencies)
        .unwrap();
    assert!(spectrum.validate_strains(&[]).is_err());
    assert!(
        spectrum
            .validate_strains(&[[C64::new(f64::NAN, 0.0); 3]; 3])
            .is_err()
    );
    let shaped = HarmonicWindow::with_profile(0.0, 1.0, HarmonicWindowProfile::Triangular).unwrap();
    assert!(spectrum.prepare_window(shaped).is_err());
    assert_eq!(spectrum.plane(), MaterialPlane::Axes12);
    assert_eq!(spectrum.angular_frequencies_rad_s(), frequencies);
    assert!(std::ptr::eq(spectrum.law(), &law));
    let f = frequencies.len();
    assert!(
        spectrum.numerical_payload_bytes()
            <= core::mem::size_of::<PlaneStressHarmonicSpectrum<'_>>()
                + 128 * f * (law.branches().len() + 1)
                + f * core::mem::size_of::<PlaneStressHarmonic>()
    );
    let window = spectrum
        .prepare_window(HarmonicWindow::new(0.0, 1.0).unwrap())
        .unwrap();
    assert!(
        window.numerical_payload_bytes()
            <= core::mem::size_of::<PlaneStressSpectrumWindow<'_, '_>>() + 256 * f + 64 * f * f
    );
    let empty = law
        .prepare_plane_stress_harmonic_spectrum(MaterialPlane::Axes12, &[])
        .unwrap()
        .prepare_window(HarmonicWindow::new(0.0, 1.0).unwrap())
        .unwrap()
        .evaluate(&[])
        .unwrap();
    assert_eq!(quantities(&empty), [0.0; 8]);
}

#[test]
fn cancelled_observation_keeps_caller_error_and_prepared_source_reusable() {
    #[derive(Debug, PartialEq)]
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
    let spectrum = law
        .prepare_plane_stress_harmonic_spectrum(MaterialPlane::Axes12, &[2.7, 3.1])
        .unwrap();
    let window = spectrum
        .prepare_window(HarmonicWindow::new(0.3, 1.7).unwrap())
        .unwrap();
    let mut checkpoints = 0;
    let failed = window.evaluate_with_checkpoint(&[in_plane(); 2], || {
        checkpoints += 1;
        if checkpoints == 4 {
            Err(Error::Cancelled)
        } else {
            Ok(())
        }
    });
    assert_eq!(failed.unwrap_err(), Error::Cancelled);
    assert!(
        window
            .evaluate(&[in_plane(); 2])
            .unwrap()
            .material_heat_j_m3
            > 0.0
    );
}

#[test]
fn relaxed_plane_bias_uses_its_dc_constraint_and_keeps_ac_heat_unchanged() {
    let law = law().with_band((0.0, 100.0)).unwrap();
    let frequencies = [1.3, 4.1, 9.7];
    let coordinates = [
        in_plane(),
        in_plane().map(|z| z * C64::new(0.6, -0.3)),
        in_plane().map(|z| z.scale(0.2)),
    ];
    let bias = [0.006, 0.002, -0.001];
    let full =
        core::array::from_fn::<_, 3, _>(|i| closed_form_strain(frequencies[i], coordinates[i]));
    let full_bias = closed_form_strain(0.0, bias.map(C64::from_re)).map(|z| z.re);
    let ac_bias = closed_form_strain(frequencies[1], bias.map(C64::from_re));
    assert!((ac_bias[2] - C64::from_re(full_bias[2])).abs() > 0.001);
    let spectrum = law
        .prepare_plane_stress_harmonic_spectrum(MaterialPlane::Axes12, &frequencies)
        .unwrap();
    spectrum
        .validate_strains_with_bias(&coordinates, &bias)
        .unwrap();
    let window = HarmonicWindow::new(0.37, 1.7).unwrap();
    let prepared = spectrum.prepare_window(window).unwrap();
    let actual = prepared.evaluate_with_bias(&coordinates, &bias).unwrap();
    let expected = law
        .prepare_harmonic_spectrum(&frequencies)
        .unwrap()
        .prepare_window(window)
        .unwrap()
        .evaluate_with_bias(&full, &full_bias)
        .unwrap();
    compare(&actual, &expected);
    let ac = prepared.evaluate(&coordinates).unwrap();
    assert_eq!(actual.material_heat_j_m3, ac.material_heat_j_m3);
    assert_ne!(actual.start_stored_energy_j_m3, ac.start_stored_energy_j_m3);
    assert_ne!(actual.stress_work_j_m3, ac.stress_work_j_m3);
}

#[test]
fn plane_bias_retains_full_static_strain_domain_and_optional_dc_applicability() {
    let law = law();
    let spectrum = law
        .prepare_plane_stress_harmonic_spectrum(MaterialPlane::Axes12, &[3.1])
        .unwrap();
    let window = HarmonicWindow::new(0.3, 1.7).unwrap();
    let prepared = spectrum.prepare_window(window).unwrap();
    let ordinary = prepared.evaluate(&[in_plane()]).unwrap();
    let zero = prepared
        .evaluate_with_bias(&[in_plane()], &[0.0; 3])
        .unwrap();
    assert_eq!(quantities(&zero), quantities(&ordinary));
    assert!(matches!(
        prepared.evaluate_with_bias(&[in_plane()], &[0.001; 3]),
        Err(ViscoError::OutOfBand { omega: 0.0, .. })
    ));
    assert!(
        spectrum
            .validate_strains_with_bias(&[in_plane()], &[0.001; 3])
            .is_err()
    );
    assert!(
        prepared
            .evaluate_with_bias(&[in_plane()], &[f64::NAN; 3])
            .is_err()
    );
    let stiff = TensorProny::new_mandel(isotropic(100.0, 1.0), vec![], 0.05).unwrap();
    let static_source = stiff
        .prepare_plane_stress_harmonic_spectrum(MaterialPlane::Axes12, &[])
        .unwrap();
    assert!(
        static_source
            .validate_strains_with_bias(&[], &[0.03, 0.03, 0.0])
            .is_err()
    );
    let bias = [0.005, -0.003, 0.001];
    let point = static_source
        .prepare_window(window)
        .unwrap()
        .evaluate_with_bias(&[], &bias)
        .unwrap();
    assert!(point.start_stored_energy_j_m3 > 0.0);
    assert_eq!(point.start_stored_energy_j_m3, point.end_stored_energy_j_m3);
    assert_eq!(point.stored_energy_change_j_m3, 0.0);
    assert_eq!(point.stress_work_j_m3, 0.0);
    assert_eq!(point.material_heat_j_m3, 0.0);
}
