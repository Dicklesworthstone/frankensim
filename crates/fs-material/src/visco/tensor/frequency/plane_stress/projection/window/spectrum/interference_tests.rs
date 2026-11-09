//! G1/G3 checks of physical interval products at coherent material nulls.

use super::*;

fn close(actual: f64, expected: f64, relative: f64) {
    assert!(
        (actual - expected).abs() <= relative * expected.abs(),
        "{actual:e} != {expected:e}"
    );
}

fn one_channel(scale: f64) -> [[f64; 6]; 6] {
    let mut gram = [[0.0; 6]; 6];
    gram[0][0] = scale;
    gram
}

fn imaginary_strains(scale: f64) -> [StrainParts; 3] {
    [1.0, -1.0, 0.25].map(|imaginary| {
        let mut strain = [[Dd::ZERO; 2]; 6];
        strain[0][1] = Dd::from_f64(imaginary * scale);
        strain
    })
}

#[test]
fn g1_kelvin_heat_resolves_a_real_three_tone_rate_null() {
    let duration = 1e-10;
    let temporal = HarmonicEnergySpectrum::new(
        &[1.0, 2.0, 4.0],
        HarmonicWindow::new(-0.5 * duration, duration).unwrap(),
    )
    .unwrap();
    let strains = imaginary_strains(1.0);
    // epsilon_dot = cos(t)-2*cos(2t)+cos(4t) = -9*t^2/2+O(t^4).
    // Physical Kelvin heat is 81*T^5/320+O(T^7), including sum tones.
    let expected = (81.0 / 320.0) * duration.powi(5);
    close(
        temporal.kelvin_heat(&one_channel(1.0), &strains).unwrap(),
        expected,
        2e-10,
    );
    let rates = [1.0, 2.0, 4.0].map(|w| [(0.0, 0), (w, 0)]);
    let mut work = Sum::default();
    temporal.add_stress_work(&one_channel(1.0), &strains, &rates, &mut work);
    close(work.value(), expected, 2e-10);
}

#[test]
fn g1_coherent_storage_increment_survives_fourth_order_energy() {
    let duration = 1e-10;
    let temporal = HarmonicEnergySpectrum::new(
        &[1.0, 2.0, 4.0],
        HarmonicWindow::new(0.0, duration).unwrap(),
    )
    .unwrap();
    let strains: [StrainParts; 3] = [1.0, -2.0, 1.0].map(|real| {
        let mut strain = [[Dd::ZERO; 2]; 6];
        strain[0][0] = Dd::from_f64(real);
        strain
    });
    let response = temporal
        .spring(
            &one_channel(1.0),
            &strains,
            &[[(1.0, 0), (0.0, 0)]; 3],
            None,
        )
        .unwrap();
    // q(0)=0 and q(T)=-9*T^2/2+O(T^4); no rounded energy subtraction.
    assert_eq!(response[0], 0.0);
    close(response[2], (81.0 / 8.0) * duration.powi(4), 2e-10);
}

#[test]
fn g3_original_viscosity_and_rate_scales_precede_physical_work() {
    let frequencies = [1.0, 2.0, 4.0].map(|w| w * 2_f64.powi(300));
    let duration = 2_f64.powi(-400);
    let temporal = HarmonicEnergySpectrum::new(
        &frequencies,
        HarmonicWindow::new(-0.5 * duration, duration).unwrap(),
    )
    .unwrap();
    let strains = imaginary_strains(2_f64.powi(-100));
    let expected = (81.0 / 320.0) * 2_f64.powi(600);
    close(
        temporal
            .kelvin_heat(&one_channel(2_f64.powi(500)), &strains)
            .unwrap(),
        expected,
        2e-12,
    );
    // Each original viscosity-times-rate phasor exceeds f64, while the
    // coherent stress-work integral is representable in physical units.
    assert!((2_f64.powi(1000) * 2_f64.powi(200)).is_infinite());
    let rates = frequencies.map(|w| [(0.0, 0), (w, 0)]);
    let mut work = Sum::default();
    temporal.add_stress_work(&one_channel(2_f64.powi(1000)), &strains, &rates, &mut work);
    close(work.value(), expected, 2e-12);
}

#[test]
fn g3_original_maxwell_tau_restores_heat_from_subnormal_memory() {
    let frequencies = [1.0, 2.0, 4.0];
    let duration = 2_f64.powi(-20);
    let tau = 2_f64.powi(-1000);
    let temporal = HarmonicEnergySpectrum::new(
        &frequencies,
        HarmonicWindow::new(-0.5 * duration, duration).unwrap(),
    )
    .unwrap();
    let strains = imaginary_strains(2_f64.powi(-200));
    // For x=omega*tau, x^2/(1+x^2)-i*x/(1+x^2) has these binary
    // factors; the neglected relative x^2 is far below DD precision.
    let spectrum = frequencies.map(|w| [(w * w, -2000), (w, -1000)]);
    let response = temporal
        .spring(
            &one_channel(2_f64.powi(500)),
            &strains,
            &spectrum,
            Some(tau),
        )
        .unwrap();
    let expected = (81.0 / 320.0) * 2_f64.powi(-500);
    assert_eq!(tau * 2_f64.powi(-200), 0.0);
    close(response[3], expected, 2e-11);
    let mut work = Sum::default();
    temporal.add_stress_work(
        &one_channel(2_f64.powi(1000)),
        &strains,
        &spectrum,
        &mut work,
    );
    close(work.value(), expected, 2e-11);
}

#[test]
fn g1_anisotropic_maxwell_and_biased_equilibrium_close_independent_real_odes() {
    let frequencies = [1.0, 2.0, 4.0];
    let tau = 0.3;
    let start = 0.17;
    let duration = 0.2;
    let strains: [StrainParts; 3] = core::array::from_fn(|f| {
        core::array::from_fn(|i| {
            [
                Dd::from_f64((f + 1) as f64 * (i as f64 - 2.0) * 0.001),
                Dd::from_f64((2.0 - f as f64) * (3.0 - i as f64) * 0.0007),
            ]
        })
    });
    let bias = [0.01, -0.003, 0.004, 0.002, -0.001, 0.005].map(Dd::from_f64);
    let mut gram = [[0.0; 6]; 6];
    gram[0] = [1.0, -2.0, 0.3, 0.0, 0.5, -0.1];
    gram[1] = [0.2, 0.0, -0.6, 0.7, 0.1, -0.2];
    let matrix = core::array::from_fn(|i| {
        core::array::from_fn(|j| (0..6).map(|k| gram[k][i] * gram[k][j]).sum())
    });
    let spectrum = frequencies.map(|w| {
        let x = w * tau;
        [(x * x / (1.0 + x * x), 0), (x / (1.0 + x * x), 0)]
    });
    let temporal =
        HarmonicEnergySpectrum::new(&frequencies, HarmonicWindow::new(start, duration).unwrap())
            .unwrap();
    let unit = [[(1.0, 0), (0.0, 0)]; 3];
    let equilibrium = temporal
        .spring_with_bias(&gram, &strains, &unit, None, &bias)
        .unwrap();
    let maxwell = temporal
        .spring(&gram, &strains, &spectrum, Some(tau))
        .unwrap();
    let mut work = Sum::default();
    temporal.add_stress_work_with_bias(&matrix, &strains, &unit, &bias, &mut work);
    temporal.add_stress_work(&matrix, &strains, &spectrum, &mut work);
    // Independently solve tau*q_dot+q=tau*epsilon_dot in real components,
    // then quadrature the original stress work and q^T*C*q/tau heat.
    let physical = |t: f64| {
        let mut strain = bias.map(Dd::to_f64);
        let mut rate = [0.0; 6];
        let mut memory = [0.0; 6];
        for f in 0..3 {
            let w = frequencies[f];
            let x = w * tau;
            let (s, c) = (w * t).sin_cos();
            for i in 0..6 {
                let re = strains[f][i][0].to_f64();
                let im = strains[f][i][1].to_f64();
                strain[i] += re * c + im * s;
                rate[i] += w * (im * c - re * s);
                memory[i] +=
                    ((x * x * re + x * im) * c + (x * x * im - x * re) * s) / (1.0 + x * x);
            }
        }
        let mut value = [0.0; 3];
        for row in gram {
            let e = (0..6).map(|i| row[i] * strain[i]).sum::<f64>();
            let q = (0..6).map(|i| row[i] * memory[i]).sum::<f64>();
            value[0] += 0.5 * (e * e + q * q);
            value[1] += q * q / tau;
        }
        for i in 0..6 {
            let stress = (0..6)
                .map(|j| matrix[i][j] * (strain[j] + memory[j]))
                .sum::<f64>();
            value[2] += stress * rate[i];
        }
        value
    };
    let mut integral = [0.0; 3];
    let steps = 2048;
    for k in 0..=steps {
        let weight = if k == 0 || k == steps {
            1.0
        } else if k % 2 == 0 {
            2.0
        } else {
            4.0
        };
        let value = physical(start + duration * k as f64 / steps as f64);
        for i in 0..3 {
            integral[i] += weight * value[i] * duration / (3 * steps) as f64;
        }
    }
    close(equilibrium[0] + maxwell[0], physical(start)[0], 3e-12);
    close(
        equilibrium[1] + maxwell[1],
        physical(start + duration)[0],
        3e-12,
    );
    close(maxwell[3], integral[1], 3e-11);
    close(work.value(), integral[2], 3e-11);
    close(
        equilibrium[2] + maxwell[2],
        physical(start + duration)[0] - physical(start)[0],
        3e-11,
    );
    close(
        work.value(),
        equilibrium[2] + maxwell[2] + maxwell[3],
        3e-11,
    );
}

#[test]
fn g3_retained_source_clock_agrees_with_translated_original_fields() {
    use fs_fft::harmonic_window::phase_rad_s_from_clock_dd;
    let frequencies = [1.0, 2.0, 4.0];
    let clock = Dd::from_pair(2_f64.powi(40), 2_f64.powi(-13));
    let duration = 2_f64.powi(-20);
    let strains = imaginary_strains(0.01);
    let translated: [StrainParts; 3] = core::array::from_fn(|f| {
        let rotation = phase_rad_s_from_clock_dd(frequencies[f], clock).unwrap();
        strains[f].map(|[re, im]| {
            [
                re * rotation[0] - im * rotation[1],
                re * rotation[1] + im * rotation[0],
            ]
        })
    });
    let actual = HarmonicEnergySpectrum::new_at_clock(&frequencies, clock, duration).unwrap();
    let shifted =
        HarmonicEnergySpectrum::new(&frequencies, HarmonicWindow::new(0.0, duration).unwrap())
            .unwrap();
    let rounded = HarmonicEnergySpectrum::new(
        &frequencies,
        HarmonicWindow::new(clock.hi, duration).unwrap(),
    )
    .unwrap();
    let heat = actual.kelvin_heat(&one_channel(1.0), &strains).unwrap();
    close(
        heat,
        shifted.kelvin_heat(&one_channel(1.0), &translated).unwrap(),
        3e-15,
    );
    assert!((heat - rounded.kelvin_heat(&one_channel(1.0), &strains).unwrap()).abs() > 1e-8 * heat);
}

#[test]
fn g3_shared_centered_coefficients_preserve_profile_parity_and_shape() {
    use fs_fft::harmonic_window::{HarmonicBandWorkspace, HarmonicWindowProfile};
    let frequencies = [0.5, 1.0, 2.0];
    let fields = [[0.3, -0.7], [-0.8, 0.5], [0.9, 0.2]].map(|z| z.map(Dd::from_f64));
    for profile in [
        HarmonicWindowProfile::Uniform,
        HarmonicWindowProfile::Triangular,
    ] {
        let window = HarmonicWindow::with_profile(0.31, 0.2, profile).unwrap();
        let band = window
            .prepare_real_band_rad_s(&frequencies)
            .unwrap()
            .unwrap();
        let mut workspace = HarmonicBandWorkspace::new();
        let borrowed = window
            .prepare_real_band_rad_s_in(&frequencies, &mut workspace)
            .unwrap()
            .unwrap();
        let mut projected = vec![[Dd::ZERO; 2]; band.mode_count()];
        band.project_parts_dd(&fields, &mut projected).unwrap();
        for (mode, value) in projected.iter().enumerate() {
            let mut direct = Dd::ZERO;
            for source in 0..3 {
                let coefficient = band.centered_coefficient_dd(mode, source).unwrap();
                assert_eq!(
                    coefficient,
                    borrowed.centered_coefficient_dd(mode, source).unwrap()
                );
                let phase = HarmonicEnergyWindow::new(frequencies[source], window)
                    .unwrap()
                    .center;
                let [re, im] = fields[source];
                let centered = [re * phase[0] - im * phase[1], re * phase[1] + im * phase[0]];
                direct = direct + coefficient * centered[mode % 2];
            }
            assert_eq!(value[1 - mode % 2], Dd::ZERO);
            assert!((direct - value[mode % 2]).to_f64().abs() < 1e-29);
        }
        assert_eq!(
            band.centered_coefficient_dd(band.mode_count(), 0),
            Err(WindowError::InvalidBandShape)
        );
        assert_eq!(
            borrowed.centered_coefficient_dd(0, 3),
            Err(WindowError::InvalidBandShape)
        );
    }
}
