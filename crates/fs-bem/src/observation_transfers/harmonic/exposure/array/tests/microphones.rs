use super::*;

fn patterns() -> [FirstOrder; 2] {
    [
        FirstOrder::new(0.25, [0.6, 0.8, 0.0]).unwrap(),
        FirstOrder::new(0.7, [-2.0 / 3.0, 1.0 / 3.0, 2.0 / 3.0]).unwrap(),
    ]
}

fn fixture() -> Fixture {
    let mut fixture = Fixture::new();
    for (f, owner) in fixture.owners.iter_mut().enumerate() {
        let z = owner.medium.density * owner.medium.sound_speed;
        for (r, field) in owner.transfers.iter_mut().enumerate() {
            for c in 0..3 {
                field[c + 1] = C64::new(
                    (0.3 + 0.2 * r as f64 - 0.1 * c as f64) / z,
                    (-0.6 + 0.4 * f as f64 + 0.3 * c as f64) / z,
                );
            }
        }
    }
    fixture
}

fn mixed_exposure(fixture: &Fixture) -> PreparedHarmonicExposure<'_> {
    let banks: [[HarmonicExposureLine<'_>; 2]; 3] = core::array::from_fn(|b| {
        [
            line(&fixture.owners[0], &fixture.amplitudes[b][..1]),
            line(&fixture.owners[1], &fixture.amplitudes[b][1..]),
        ]
    });
    PreparedHarmonicExposure::new_with_statistics(
        &fixture.owners[0],
        "directional-clock",
        &banks[0],
        &[&banks[1], &banks[2]],
        &[
            HarmonicExposureFactorStatistics::Real,
            HarmonicExposureFactorStatistics::ProperComplex,
        ],
        preparation(),
    )
    .unwrap()
}

// Independent real-time first-order observation and Simpson integration over
// the exact +/-1 real and four-phase proper ensemble. Projection is applied
// to real p/u at every delayed sample, not to prepared channel phasors.
fn quadrature(
    fixture: &Fixture,
    patterns: &[FirstOrder],
    channels: &[HarmonicPressureChannel],
    window: HarmonicWindow,
) -> (Vec<HarmonicPressureChannelMoments>, Vec<f64>) {
    let proper = [C64::ONE, C64::new(0.0, 1.0), -C64::ONE, C64::new(0.0, -1.0)];
    let mut means = vec![vec![0.0; channels.len()]; 8];
    let mut moments = vec![HarmonicPressureChannelMoments::default(); channels.len()];
    let steps = 2048;
    for (state, average) in means.iter_mut().enumerate() {
        let real = if state % 2 == 0 { 1.0 } else { -1.0 };
        for sample in 0..=steps {
            let u = sample as f64 / steps as f64;
            let time = window.start_s() + window.duration_s() * u;
            let simpson = if sample == 0 || sample == steps {
                1.0
            } else if sample % 2 == 0 {
                2.0
            } else {
                4.0
            };
            let sensitivity = match window.profile() {
                HarmonicWindowProfile::Uniform => 1.0,
                HarmonicWindowProfile::Triangular => 2.0 * (1.0 - (2.0 * u - 1.0).abs()),
            };
            let weight = simpson * sensitivity / (3.0 * steps as f64);
            for (c, channel) in channels.iter().enumerate() {
                let mut y = 0.0;
                for tap in &channel.taps {
                    let microphone = patterns[tap.receiver];
                    for f in 0..2 {
                        let owner = &fixture.owners[f];
                        let coefficient = fixture.amplitudes[0][f]
                            + fixture.amplitudes[1][f].scale(real)
                            + fixture.amplitudes[2][f] * proper[state / 2];
                        let angle = owner.omega * (time - tap.delay_s);
                        let (sin, cos) = angle.sin_cos();
                        let values = owner.transfers[tap.receiver].map(|transfer| {
                            let amplitude = transfer * coefficient;
                            amplitude.re * cos + amplitude.im * sin
                        });
                        let axial: f64 = microphone
                            .front_axis()
                            .iter()
                            .zip(&values[1..])
                            .map(|(n, velocity)| n * velocity)
                            .sum();
                        y += tap.weight
                            * (microphone.pressure_fraction() * values[0]
                                - (1.0 - microphone.pressure_fraction())
                                    * owner.medium.density
                                    * owner.medium.sound_speed
                                    * axial);
                    }
                }
                average[c] += weight * y;
                moments[c].mean_square_pressure_pa2 += weight * y * y / 8.0;
            }
        }
    }
    for (c, moment) in moments.iter_mut().enumerate() {
        moment.mean_pressure_pa = means.iter().map(|mean| mean[c] / 8.0).sum();
    }
    let n = channels.len();
    let covariance = (0..n * n)
        .map(|index| {
            let row = index / n;
            let col = index % n;
            means
                .iter()
                .map(|v| {
                    (v[row] - moments[row].mean_pressure_pa)
                        * (v[col] - moments[col].mean_pressure_pa)
                        / 8.0
                })
                .sum()
        })
        .collect();
    (moments, covariance)
}

#[test]
fn directional_arrays_match_delayed_mixed_ensemble_time_quadrature() {
    let fixture = fixture();
    let exposure = mixed_exposure(&fixture);
    let channels = channels();
    let patterns = patterns();
    let array = PreparedHarmonicPressureArray::new_with_microphones(
        &exposure,
        &patterns,
        &channels,
        limits(),
    )
    .unwrap();
    drop(exposure);
    assert_eq!(array.microphone_patterns(), Some(patterns.as_slice()));
    for profile in [
        HarmonicWindowProfile::Uniform,
        HarmonicWindowProfile::Triangular,
    ] {
        let window = HarmonicWindow::with_profile(0.17, 0.43, profile).unwrap();
        let observed = array.observe_window(window).unwrap();
        let (expected, covariance) = quadrature(&fixture, &patterns, &channels, window);
        for (actual, expected) in observed.moments().iter().zip(expected) {
            close(actual.mean_pressure_pa, expected.mean_pressure_pa, 3e-11);
            close(
                actual.mean_square_pressure_pa2,
                expected.mean_square_pressure_pa2,
                3e-11,
            );
        }
        for (&actual, expected) in observed.averaged_covariance_pa2().iter().zip(covariance) {
            close(actual, expected, 3e-11);
        }
        assert_eq!(observed.microphone_patterns(), Some(patterns.as_slice()));
        assert_eq!(observed.physical_receivers(), fixture.owners[0].receivers);
        assert_eq!(observed.source_clock(), "directional-clock");
    }
}

#[test]
fn a_standing_wave_pressure_node_still_drives_an_axial_velocity_microphone() {
    let mut reference = owner(1.0, &[C64::ZERO; 3]);
    reference.medium.density = 1.0;
    reference.medium.sound_speed = 1.0;
    for (receiver, field) in reference.receivers.iter_mut().zip(&mut reference.transfers) {
        // At k*x=pi/2, two unit counterpropagating waves have p=0, u_x=2i/Z.
        receiver.position_m = [core::f64::consts::FRAC_PI_2, 0.0, 0.0];
        field[1] = C64::new(0.0, 2.0);
    }
    let unit = [C64::ONE];
    let mean = [line(&reference, &unit)];
    let exposure =
        PreparedHarmonicExposure::new(&reference, "standing-wave", &mean, &[], preparation())
            .unwrap();
    let channels: Vec<_> = (0..3)
        .map(|r| channel(&format!("pickup-{r}"), &[tap(r, 1.0, 0.0)]))
        .collect();
    let patterns = [
        FirstOrder::default(),
        FirstOrder::new(0.0, [-1.0, 0.0, 0.0]).unwrap(),
        FirstOrder::new(0.0, [0.0, 1.0, 0.0]).unwrap(),
    ];
    let array = PreparedHarmonicPressureArray::new_with_microphones(
        &exposure,
        &patterns,
        &channels,
        limits(),
    )
    .unwrap();
    let result = array
        .observe_window(HarmonicWindow::new(0.0, core::f64::consts::FRAC_PI_2).unwrap())
        .unwrap();
    assert_eq!(
        result.moments()[0],
        HarmonicPressureChannelMoments::default()
    );
    assert_eq!(
        result.moments()[2],
        HarmonicPressureChannelMoments::default()
    );
    close(
        result.moments()[1].mean_pressure_pa,
        4.0 / core::f64::consts::PI,
        1e-14,
    );
    close(result.moments()[1].mean_square_pressure_pa2, 2.0, 1e-14);
}

#[test]
fn microphone_rotation_and_explicit_omni_preserve_the_original_field() {
    let fixture = fixture();
    let exposure = mixed_exposure(&fixture);
    let channels = channels();
    let patterns = patterns();
    let array = PreparedHarmonicPressureArray::new_with_microphones(
        &exposure,
        &patterns,
        &channels,
        limits(),
    )
    .unwrap();
    let mut rotated = self::fixture();
    for owner in &mut rotated.owners {
        for field in &mut owner.transfers {
            let [pressure, x, y, z] = *field;
            *field = [pressure, -y, x, z];
        }
        for receiver in &mut owner.receivers {
            let [x, y, z] = receiver.position_m;
            receiver.position_m = [-y, x, z];
        }
    }
    let rotated_patterns = patterns.map(|pattern| {
        let [x, y, z] = pattern.front_axis();
        FirstOrder::new(pattern.pressure_fraction(), [-y, x, z]).unwrap()
    });
    let rotated_exposure = mixed_exposure(&rotated);
    let rotated_array = PreparedHarmonicPressureArray::new_with_microphones(
        &rotated_exposure,
        &rotated_patterns,
        &channels,
        limits(),
    )
    .unwrap();
    let window =
        HarmonicWindow::with_profile(-0.13, 0.71, HarmonicWindowProfile::Triangular).unwrap();
    let baseline = array.observe_window(window).unwrap();
    let result = rotated_array.observe_window(window).unwrap();
    for (a, b) in result.moments().iter().zip(baseline.moments()) {
        close(a.mean_pressure_pa, b.mean_pressure_pa, 1e-14);
        close(
            a.mean_square_pressure_pa2,
            b.mean_square_pressure_pa2,
            1e-14,
        );
    }
    for (&a, &b) in result
        .averaged_covariance_pa2()
        .iter()
        .zip(baseline.averaged_covariance_pa2())
    {
        close(a, b, 1e-14);
    }
    let pressure = PreparedHarmonicPressureArray::new(&exposure, &channels, limits()).unwrap();
    let omni = PreparedHarmonicPressureArray::new_with_microphones(
        &exposure,
        &[FirstOrder::default(); 2],
        &channels,
        limits(),
    )
    .unwrap();
    let a = pressure.observe_window(window).unwrap();
    let b = omni.observe_window(window).unwrap();
    assert_eq!(a.moments(), b.moments());
    assert_eq!(a.averaged_covariance_pa2(), b.averaged_covariance_pa2());
    assert_eq!(a.microphone_patterns(), None);
}

#[test]
fn oblique_velocity_projection_keeps_a_weak_survivor_between_large_components() {
    let mut reference = owner(1.0, &[C64::ZERO]);
    reference.medium.density = 1.0;
    reference.medium.sound_speed = 1.0;
    reference.transfers[0][1..].copy_from_slice(&[
        C64::from_re(1e16),
        C64::ONE,
        C64::from_re(-1e16),
    ]);
    let unit = [C64::ONE];
    let mean = [line(&reference, &unit)];
    let exposure =
        PreparedHarmonicExposure::new(&reference, "vector-null", &mean, &[], preparation())
            .unwrap();
    let axis = 1.0 / 3.0_f64.sqrt();
    let patterns = [FirstOrder::new(0.0, [axis; 3]).unwrap()];
    let array = PreparedHarmonicPressureArray::new_with_microphones(
        &exposure,
        &patterns,
        &[channel("oblique", &[tap(0, 1.0, 0.0)])],
        limits(),
    )
    .unwrap();
    let result = array
        .observe_window(HarmonicWindow::new(-0.5, 1.0).unwrap())
        .unwrap();
    close(
        result.moments()[0].mean_pressure_pa,
        -2.0 * axis * 0.5_f64.sin(),
        3e-15,
    );
    close(
        result.moments()[0].mean_square_pressure_pa2,
        0.5 * axis * axis * (1.0 + 1.0_f64.sin()),
        3e-15,
    );
}

#[test]
fn microphone_admission_budgets_and_cancellation_leave_reusable_inputs() {
    let fixture = fixture();
    let exposure = mixed_exposure(&fixture);
    let channels = channels();
    let patterns = patterns();
    for count in [0, 1, 3] {
        assert!(matches!(
            PreparedHarmonicPressureArray::new_with_microphones(
                &exposure,
                &vec![FirstOrder::default(); count],
                &channels,
                limits()
            ),
            Err(TransferError::Invalid(_))
        ));
    }
    let array = PreparedHarmonicPressureArray::new_with_microphones(
        &exposure,
        &patterns,
        &channels,
        limits(),
    )
    .unwrap();
    let mut budget = limits();
    budget.maximum_workspace_bytes = array.workspace_bytes();
    assert!(
        PreparedHarmonicPressureArray::new_with_microphones(
            &exposure, &patterns, &channels, budget
        )
        .is_ok()
    );
    budget.maximum_workspace_bytes -= 1;
    assert!(matches!(
        PreparedHarmonicPressureArray::new_with_microphones(
            &exposure, &patterns, &channels, budget
        ),
        Err(TransferError::Cap { .. })
    ));
    budget = limits();
    // Pressure-only preparation is admitted at exactly its own product count; the full p/u work is not free.
    // array.rs `prepare_products` = 3 lines taps + nonempty_banks taps + lines bank_count + microphone_products.
    // Here lines 6 (three banks of two), taps 7 (2 + 2 + 3), nonempty banks 3, bank_count 3 (two factors and
    // the reference) and no microphone term without patterns: 3*6*7 + 3*7 + 6*3 = 165. The p/u preparation adds
    // 2 microphones * (4 lines + 8) = 64, so it needs 229. (54a084112a moved the per-tap factor from 2 to 3 and
    // added the bank term; the earlier literal 128 sat between the previous 102 and 166.)
    budget.maximum_prepare_products = 165;
    assert!(PreparedHarmonicPressureArray::new(&exposure, &channels, budget).is_ok());
    assert!(matches!(
        PreparedHarmonicPressureArray::new_with_microphones(
            &exposure, &patterns, &channels, budget
        ),
        Err(TransferError::Cap { .. })
    ));
    let mut checks = 0;
    assert!(matches!(
        PreparedHarmonicPressureArray::new_with_microphones_cancellable(
            &exposure,
            &patterns,
            &channels,
            limits(),
            || {
                checks += 1;
                checks == 40
            }
        ),
        Err(TransferError::Cancelled)
    ));
    let replay = PreparedHarmonicPressureArray::new_with_microphones(
        &exposure,
        &patterns,
        &channels,
        limits(),
    )
    .unwrap();
    let window = HarmonicWindow::new(0.0, 0.1).unwrap();
    assert_eq!(
        array.observe_window(window).unwrap().moments(),
        replay.observe_window(window).unwrap().moments()
    );
    let dark = PreparedHarmonicPressureArray::new_with_microphones(
        &exposure,
        &patterns,
        &[channel("dark", &[])],
        limits(),
    )
    .unwrap();
    assert_eq!(
        dark.observe_window(window).unwrap().moments(),
        &[HarmonicPressureChannelMoments::default()]
    );
}
