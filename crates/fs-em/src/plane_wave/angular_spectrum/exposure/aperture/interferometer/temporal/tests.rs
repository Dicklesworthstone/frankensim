use super::*;
use fs_fft::harmonic_window::HarmonicWindowProfile;

fn wave(a: C64) -> [C64; 6] {
    [a, C64::ZERO, C64::ZERO, C64::ZERO, a.scale(-1.0), C64::ZERO]
}

#[allow(clippy::too_many_arguments)]
fn owner(
    frequencies: &[f64],
    kx: &[f64],
    ky: &[f64],
    values: &[Vec<[C64; 6]>],
    arms: usize,
    aperture: AngularSpectrumAperture,
    window: HarmonicWindow,
    proper: bool,
    positive: bool,
) -> PreparedAngularSpectrumApertureInterferometer {
    let grid = AngularSpectrumGrid {
        width: kx.len(),
        height: ky.len(),
        spacing_m: [1.0; 2],
    };
    let n = kx.len() * ky.len();
    let f = frequencies.len();
    let mut fields = Vec::new();
    for (frequency, field) in values.iter().enumerate() {
        for arm in 0..arms {
            let mut input: Vec<FieldParts> = field[arm * n..(arm + 1) * n]
                .iter()
                .map(|v| v.map(|z| [Dd::from_f64(z.re), Dd::from_f64(z.im)]))
                .collect();
            if !proper {
                center_fields(
                    &mut input,
                    frequencies[frequency],
                    |i| kx[i % kx.len()],
                    |i| ky[i / kx.len()],
                    aperture,
                    window,
                    &mut || false,
                )
                .unwrap();
            }
            fields.extend(input);
        }
    }
    let factor =
        positive.then(|| Positive::new(window, frequencies, proper, &mut || false).unwrap());
    let mode_count = factor.as_ref().map_or(0, |v| v.factor.mode_count());
    let temporal = if proper {
        Temporal::Proper(pair_kernels(f, window, |i| -frequencies[i], &mut || false).unwrap())
    } else {
        Temporal::Real(pair_kernels(f, window, |i| -frequencies[i], &mut || false).unwrap())
    };
    PreparedAngularSpectrumApertureInterferometer {
        arm_count: arms,
        modes: n,
        side: AngularSpectrumIlluminationSide::TotalLeft,
        aperture,
        window,
        banks: vec![Bank {
            frequencies: frequencies.to_vec(),
            fields,
            temporal,
            positive: factor,
        }],
        spatial_real: (!proper).then(|| {
            Spatial::from_axes(grid, aperture, |i| kx[i], |i| ky[i], &mut || false).unwrap()
        }),
        spatial_proper: proper.then(|| {
            Spatial::from_axes(grid, aperture, |i| kx[i], |i| ky[i], &mut || false).unwrap()
        }),
        gram: GramWorkspace::new(n, f, if proper { 0 } else { f }, &mut || false).unwrap(),
        projected: vec![[[Dd::ZERO; 2]; 6]; mode_count * n],
    }
}

#[test]
fn real_repeated_two_tone_readout_keeps_dark_residual_and_signed_power() {
    let duration: f64 = 1e-12;
    let area = AngularSpectrumAperture::try_new([0.0; 2], [0.7, 0.9]).unwrap();
    let values = [vec![wave(C64::ONE)], vec![wave(C64::from_re(-1.0))]];
    for (profile, coefficient) in [
        (HarmonicWindowProfile::Uniform, 9.0 / 320.0),
        (HarmonicWindowProfile::Triangular, 9.0 / 960.0),
    ] {
        let window = HarmonicWindow::with_profile(-duration / 2.0, duration, profile).unwrap();
        let mut owner = owner(
            &[1.0, 2.0],
            &[0.0],
            &[0.0],
            &values,
            1,
            area,
            window,
            false,
            true,
        );
        let expected = coefficient * duration.powi(4);
        let pointer = owner.projected.as_ptr();
        let read = owner
            .measure_moments(&[C64::ONE], &[0.0], || false)
            .unwrap();
        assert!((read.electric_second_moment_v2_m2()[0][0] / expected - 1.0).abs() < 2e-12);
        assert!((read.reading().normal_power_w / (expected * area.area_m2()) + 1.0).abs() < 2e-12);
        assert_eq!(
            read,
            owner
                .measure_moments(&[C64::ONE], &[0.0], || false)
                .unwrap()
        );
        assert_eq!(pointer, owner.projected.as_ptr());
    }
}

#[test]
fn proper_optical_null_adds_independent_banks_without_inventing_coherence() {
    let duration: f64 = 1e-12;
    let area = AngularSpectrumAperture::try_new([0.0; 2], [1.0; 2]).unwrap();
    let values: Vec<_> = [1.0, -2.0, 1.0].map(|a| vec![wave(C64::from_re(a))]).into();
    let window =
        HarmonicWindow::with_profile(-duration / 2.0, duration, HarmonicWindowProfile::Triangular)
            .unwrap();
    let frequencies = [1e15, 1e15 + 1.0, 1e15 + 2.0];
    let mut first = owner(
        &frequencies,
        &[0.0],
        &[0.0],
        &values,
        1,
        area,
        window,
        true,
        true,
    );
    let mut second = owner(
        &frequencies,
        &[0.0],
        &[0.0],
        &values,
        1,
        area,
        window,
        true,
        true,
    );
    first.banks.append(&mut second.banks);
    let read = first
        .measure_moments(&[C64::ONE], &[0.0], || false)
        .unwrap();
    let expected = 2.0 * duration.powi(4) / 480.0;
    assert!((read.electric_second_moment_v2_m2()[0][0] / expected - 1.0).abs() < 2e-12);
}

#[test]
fn all_moments_preserve_complex_arm_gains_signed_delays_and_spatial_parity() {
    let (frequencies, kx, ky) = ([1.3, 2.1], [0.2, -0.7], [0.3, -0.5]);
    let values: Vec<_> = (0..2)
        .map(|f| {
            (0..8)
                .map(|m| {
                    core::array::from_fn(|c| {
                        C64::new(
                            ((f + 1) * (c + 1)) as f64 / (7.0 * (m + 1) as f64),
                            (c as f64 - 2.0) * (m as f64 - 1.5) / (11.0 * (f + 1) as f64),
                        )
                    })
                })
                .collect()
        })
        .collect();
    let gains = [C64::new(0.8, 0.3), C64::new(-0.6, 0.1)];
    let delays = [0.17, -0.23];
    let area = AngularSpectrumAperture::try_new([-0.19, 0.23], [0.7, 0.9]).unwrap();
    for profile in [
        HarmonicWindowProfile::Uniform,
        HarmonicWindowProfile::Triangular,
    ] {
        let window = HarmonicWindow::with_profile(0.17, 0.43, profile).unwrap();
        for proper in [false, true] {
            let mut stable = owner(
                &frequencies,
                &kx,
                &ky,
                &values,
                2,
                area,
                window,
                proper,
                true,
            );
            let mut original = owner(
                &frequencies,
                &kx,
                &ky,
                &values,
                2,
                area,
                window,
                proper,
                false,
            );
            let actual = stable.measure_moments(&gains, &delays, || false).unwrap();
            let expected = original.measure_moments(&gains, &delays, || false).unwrap();
            for (a, b) in [
                actual.electric_second_moment_v2_m2(),
                actual.magnetic_second_moment_a2_m2(),
                actual.electric_magnetic_second_moment_w_m2(),
            ]
            .into_iter()
            .flatten()
            .flatten()
            .zip(
                [
                    expected.electric_second_moment_v2_m2(),
                    expected.magnetic_second_moment_a2_m2(),
                    expected.electric_magnetic_second_moment_w_m2(),
                ]
                .into_iter()
                .flatten()
                .flatten(),
            ) {
                assert!((a - b).abs() < 3e-14, "{a:e} != {b:e}");
            }
        }
    }
}

#[test]
fn coherent_arm_null_and_cancellation_leave_cached_observer_reusable() {
    let area = AngularSpectrumAperture::try_new([0.0; 2], [1.0; 2]).unwrap();
    let window = HarmonicWindow::new(0.0, 0.1).unwrap();
    let values = [vec![wave(C64::ONE); 2], vec![wave(C64::new(0.3, 0.2)); 2]];
    for proper in [false, true] {
        let mut owner = owner(
            &[1.0, 2.0],
            &[0.0],
            &[0.0],
            &values,
            2,
            area,
            window,
            proper,
            true,
        );
        let null = owner
            .measure_moments(&[C64::ONE, C64::from_re(-1.0)], &[0.0; 2], || false)
            .unwrap();
        assert_eq!(null.reading().normal_power_w, 0.0);
        let gains = [C64::ONE, C64::new(0.1, -0.3)];
        let expected = owner
            .measure_moments(&gains, &[0.2, -0.1], || false)
            .unwrap();
        let mut polls = 0;
        assert!(matches!(
            owner.measure_moments(&gains, &[0.2, -0.1], || {
                polls += 1;
                polls == 10
            }),
            Err(AngularSpectrumError::Cancelled)
        ));
        assert_eq!(
            expected,
            owner
                .measure_moments(&gains, &[0.2, -0.1], || false)
                .unwrap()
        );
    }
}

#[test]
fn cached_factor_bounds_count_retained_modes_rotations_and_shared_scratch_once() {
    let grid = AngularSpectrumGrid {
        width: 2,
        height: 3,
        spacing_m: [1.0; 2],
    };
    let window = HarmonicWindow::new(0.0, 1e-12).unwrap();
    let first = Plan::new(grid, window, &[1.0, 2.0], false, 0)
        .unwrap()
        .unwrap();
    let next = Plan::new(grid, window, &[1.0, 2.0], false, first.modes)
        .unwrap()
        .unwrap();
    assert_eq!(
        next.extra_bytes,
        2 * (first.modes + 2) * core::mem::size_of::<Dd>()
    );
    assert_eq!(
        first.extra_bytes - next.extra_bytes,
        FIXED_BYTES + 6 * first.modes * core::mem::size_of::<FieldParts>()
    );
    assert_eq!(
        first.extra_work,
        (6 * first.modes * (2 + 2 + 3 + 1)).saturating_sub(6 * 2 * (2 + 3 + 2))
    );
    assert!(Plan::new(grid, window, &[1.0], false, 0).unwrap().is_none());
    assert!(
        Plan::new(
            grid,
            HarmonicWindow::new(0.0, 2.0).unwrap(),
            &[1.0, 2.0],
            false,
            0
        )
        .unwrap()
        .is_none()
    );
}

#[test]
fn delayed_optical_dark_envelope_retains_shifted_fourth_moment() {
    let duration: f64 = 1e-12;
    let area = AngularSpectrumAperture::try_new([0.0; 2], [1.0; 2]).unwrap();
    let values: Vec<_> = [1.0, -2.0, 1.0]
        .map(|a| vec![wave(C64::from_re(a)); 2])
        .into();
    let frequencies = [1e15, 1e15 + 1.0, 1e15 + 2.0];
    for (profile, variance, fourth) in [
        (HarmonicWindowProfile::Uniform, 1.0 / 12.0, 1.0 / 80.0),
        (HarmonicWindowProfile::Triangular, 1.0 / 24.0, 1.0 / 240.0),
    ] {
        let window = HarmonicWindow::with_profile(-duration / 2.0, duration, profile).unwrap();
        for delay in [8e-16_f64, 1e-13, 1e-12] {
            let expected = 0.5
                * (delay.powi(4)
                    + 6.0 * delay.powi(2) * variance * duration.powi(2)
                    + fourth * duration.powi(4));
            // A(t-delay)=[1-exp(i*(t-delay))]^2 times a common carrier;
            // its leading proper real-field power is <(t-delay)^4>/2.
            // The first inactive arm must not determine its physical phase.
            let mut observed = owner(
                &frequencies,
                &[0.0],
                &[0.0],
                &values,
                2,
                area,
                window,
                true,
                true,
            );
            let actual = observed
                .measure_moments(&[C64::ZERO, C64::ONE], &[0.17, delay], || false)
                .unwrap()
                .electric_second_moment_v2_m2()[0][0];
            assert!(
                (actual / expected - 1.0).abs() < 3e-12,
                "profile={profile:?}, delay={delay:e}, actual={actual:e}, expected={expected:e}"
            );
            let coherent = observed
                .measure_moments(&[C64::ONE, C64::from_re(-0.5)], &[delay; 2], || false)
                .unwrap()
                .electric_second_moment_v2_m2()[0][0];
            assert!((coherent / (0.25 * expected) - 1.0).abs() < 3e-12);
        }
    }
}

#[test]
fn proper_delay_projection_admits_added_carrier_phase_work() {
    assert_eq!(delay_projection_work(3, 2, true).unwrap(), 8);
    assert_eq!(delay_projection_work(3, 2, false).unwrap(), 0);
    assert!(delay_projection_work(usize::MAX, 2, true).is_err());
}

#[test]
fn overflowing_relative_carrier_keeps_valid_independent_delay_phases() {
    let area = AngularSpectrumAperture::try_new([0.0; 2], [1.0; 2]).unwrap();
    let window = HarmonicWindow::new(0.0, 0.1).unwrap();
    let values = [vec![wave(C64::ONE); 2], vec![wave(C64::new(0.3, 0.2)); 2]];
    let mut stable = owner(
        &[2.0, 3.0],
        &[0.0],
        &[0.0],
        &values,
        2,
        area,
        window,
        true,
        true,
    );
    let mut original = owner(
        &[2.0, 3.0],
        &[0.0],
        &[0.0],
        &values,
        2,
        area,
        window,
        true,
        false,
    );
    let delay = 0.3 * f64::MAX;
    let gains = [C64::ONE, C64::from_re(0.2)];
    let a = stable
        .measure_moments(&gains, &[delay, -delay], || false)
        .unwrap();
    let b = original
        .measure_moments(&gains, &[delay, -delay], || false)
        .unwrap();
    let (a, b) = (
        a.electric_second_moment_v2_m2()[0][0],
        b.electric_second_moment_v2_m2()[0][0],
    );
    assert!((a - b).abs() < 3e-14, "{a:e} != {b:e}");
}
