//! Original Lorentz P/J energy oracles for the shared finite-window contraction.

use super::*;
use fs_material::field_storage::dielectric::lorentz::LorentzPole;

type Electric = [[C64; 3]; 2];

fn electric() -> Electric {
    [
        [C64::new(0.8, 0.2), C64::new(-0.1, 0.3), C64::new(0.2, -0.4)],
        [C64::new(-0.3, 0.6), C64::new(0.4, 0.1), C64::new(-0.2, 0.1)],
    ]
}

fn near(actual: f64, expected: f64, scale: f64) {
    assert!(
        (actual - expected).abs() <= 5e-11 * scale.max(1e-30),
        "actual {actual:.17e}, expected {expected:.17e}, scale {scale:.17e}"
    );
}

fn phase(omega: f64, time: f64) -> C64 {
    let (sine, cosine) = (omega * time).sin_cos();
    C64::new(cosine, sine)
}

fn filtered(
    pole: LorentzPole,
    frequencies: [f64; 2],
    electric: Electric,
) -> (Vec<[C64; 6]>, Vec<[C64; 3]>) {
    let mut fields = Vec::new();
    let mut kinetic = Vec::new();
    for i in 0..2 {
        let transfer = PoleTransfer::lorentz(pole, frequencies[i]).unwrap();
        fields.push(core::array::from_fn(|c| {
            electric[i][c % 3]
                * if c < 3 {
                    transfer.heat
                } else {
                    transfer.storage
                }
        }));
        kinetic.push(electric[i].map(|e| e * transfer.kinetic_storage));
    }
    (fields, kinetic)
}

fn real_kernels(
    frequencies: [f64; 2],
    window: HarmonicWindow,
) -> Vec<MomentKernel<HarmonicProduct>> {
    frequencies
        .into_iter()
        .flat_map(|a| {
            frequencies.into_iter().map(move |b| MomentKernel {
                product: window.product_kernel(a, b).unwrap(),
                storage_frequency: Dd::from_f64(a),
            })
        })
        .collect()
}

// Independent moderate-scale original ODE, without the normalized transfer:
// P = delta*omega0^2 E/(omega0^2-omega^2+i*gamma*omega), J=i*omega*P.
fn original_phasors(
    pole: LorentzPole,
    frequencies: [f64; 2],
    electric: Electric,
) -> (Electric, Electric) {
    let natural = pole.angular_frequency_rad_s();
    let polarization = core::array::from_fn(|i| {
        let omega = frequencies[i];
        let response = C64::from_re(pole.permittivity_increment_f_m() * natural * natural)
            / C64::new(
                natural * natural - omega * omega,
                pole.damping_rate_per_s() * omega,
            );
        electric[i].map(|e| e * response)
    });
    let current =
        core::array::from_fn(|i| polarization[i].map(|p| rotate_i(p).scale(frequencies[i])));
    (polarization, current)
}

fn original_energy(
    pole: LorentzPole,
    frequencies: [f64; 2],
    polarization: Electric,
    current: Electric,
    time: f64,
) -> [f64; 2] {
    let mut energy = [0.0; 2];
    for c in 0..3 {
        let p: f64 = (0..2)
            .map(|i| (polarization[i][c] * phase(frequencies[i], time)).re)
            .sum();
        let j: f64 = (0..2)
            .map(|i| (current[i][c] * phase(frequencies[i], time)).re)
            .sum();
        energy[0] += p * p / (2.0 * pole.permittivity_increment_f_m());
        energy[1] += (j / pole.angular_frequency_rad_s()).powi(2)
            / (2.0 * pole.permittivity_increment_f_m());
    }
    energy
}

fn mean_product(window: HarmonicWindow, frequencies: [f64; 2], a: Electric, b: Electric) -> f64 {
    let mut sum = Dd::ZERO;
    for i in 0..2 {
        for j in 0..2 {
            let kernel = window
                .product_kernel(frequencies[i], frequencies[j])
                .unwrap();
            for c in 0..3 {
                // Actual shared owner uses exp(-i omega t); the original
                // physical phasors in this test use exp(+i omega t).
                sum = sum
                    + kernel
                        .mean_real_product_dd(a[i][c].conj(), b[j][c].conj())
                        .unwrap();
            }
        }
    }
    sum.to_f64()
}

#[test]
fn lorentz_short_real_windows_retain_both_endpoint_stores_and_original_work() {
    let pole = LorentzPole::try_new(0.7, 7.0, 0.4).unwrap();
    let frequencies = [6.8, 7.15];
    let electric = electric();
    let (fields, kinetic) = filtered(pole, frequencies, electric);
    let (polarization, current) = original_phasors(pole, frequencies, electric);
    for (start, duration) in [(0.17, 0.012), (0.17, 0.12), (0.41, 0.079)] {
        let window = HarmonicWindow::new(start, duration).unwrap();
        let kernels = real_kernels(frequencies, window);
        let (heat, storage) =
            pole_moments(&kernels, &fields, Some(&kinetic), product, &mut || false).unwrap();
        let (_, restoring_only) =
            pole_moments(&kernels, &fields, None, product, &mut || false).unwrap();
        let initial = original_energy(pole, frequencies, polarization, current, start);
        let final_energy =
            original_energy(pole, frequencies, polarization, current, start + duration);
        let rates: [f64; 2] = core::array::from_fn(|c| (final_energy[c] - initial[c]) / duration);
        let work = mean_product(window, frequencies, electric, current);
        let heat_reference = mean_product(window, frequencies, current, current)
            * pole.damping_rate_per_s()
            / (pole.permittivity_increment_f_m() * pole.angular_frequency_rad_s().powi(2));
        let scale = work
            .abs()
            .max(heat_reference)
            .max(rates[0].abs() + rates[1].abs());
        assert!(heat.to_f64() > 0.0 && rates[1].abs() > 1e-3);
        near(heat.to_f64(), heat_reference, scale);
        near(restoring_only.to_f64(), rates[0], scale);
        near(storage.to_f64(), rates[0] + rates[1], scale);
        near((heat + storage).to_f64(), work, scale);

        let first = 0.37 * duration;
        let mut split = [0.0; 2];
        for (part_start, part_duration) in [(start, first), (start + first, duration - first)] {
            let part = real_kernels(
                frequencies,
                HarmonicWindow::new(part_start, part_duration).unwrap(),
            );
            let (q, rate) =
                pole_moments(&part, &fields, Some(&kinetic), product, &mut || false).unwrap();
            split[0] += q.to_f64() * part_duration;
            split[1] += rate.to_f64() * part_duration;
        }
        near(split[0] / duration, heat.to_f64(), scale);
        near(split[1] / duration, storage.to_f64(), scale);
    }
}

#[test]
fn lorentz_proper_covariance_retains_kinetic_slow_beats_without_carrier_derivatives() {
    let window = HarmonicWindow::new(0.17, 0.12).unwrap();
    for carrier in [0.0, 1_125_899_906_842_624.0] {
        let frequencies = [carrier + 2.0, carrier + 2.25];
        let pole = LorentzPole::try_new(0.7, 2.0 * (carrier + 4.0), 0.4 * (carrier + 1.0)).unwrap();
        let (fields, kinetic) = filtered(pole, frequencies, electric());
        let mut kernels = Vec::new();
        for a in frequencies {
            for b in frequencies {
                kernels.push(MomentKernel {
                    product: window.proper_product_kernel(a, b).unwrap(),
                    storage_frequency: (Dd::from_f64(a) - Dd::from_f64(b)) * Dd::from_f64(0.5),
                });
            }
        }
        let (heat, storage) = pole_moments(
            &kernels,
            &fields,
            Some(&kinetic),
            proper_product,
            &mut || false,
        )
        .unwrap();
        let (_, restoring_only) =
            pole_moments(&kernels, &fields, None, proper_product, &mut || false).unwrap();
        // Exact proper-factor covariance has only the slow difference term.
        // E[U]=(|R(t)|^2+|K(t)|^2)/4, so the cross derivative is
        // Re(i*delta*(R0*conj(R1)+K0*conj(K1))*mean_phase(delta))/2.
        let delta = frequencies[0] - frequencies[1];
        let half_width = delta * 0.5 * window.duration_s();
        let beat = phase(delta, window.start_s() + 0.5 * window.duration_s())
            .scale(half_width.sin() / half_width);
        let mut restoring_cross = C64::ZERO;
        let mut kinetic_cross = C64::ZERO;
        let mut heat_cross = C64::ZERO;
        let mut diagonal_heat = 0.0;
        for c in 0..3 {
            restoring_cross = restoring_cross + fields[0][c + 3] * fields[1][c + 3].conj();
            kinetic_cross = kinetic_cross + kinetic[0][c] * kinetic[1][c].conj();
            heat_cross = heat_cross + fields[0][c] * fields[1][c].conj();
            diagonal_heat += 0.5 * (fields[0][c].norm_sq() + fields[1][c].norm_sq());
        }
        let restoring_rate = 0.5 * delta * (rotate_i(restoring_cross) * beat).re;
        let kinetic_rate = 0.5 * delta * (rotate_i(kinetic_cross) * beat).re;
        let expected_heat = diagonal_heat + (heat_cross * beat).re;
        assert!(kinetic_rate.abs() > 1e-5);
        near(
            storage.to_f64(),
            restoring_rate + kinetic_rate,
            restoring_rate.abs() + kinetic_rate.abs(),
        );
        near(
            restoring_only.to_f64(),
            restoring_rate,
            restoring_rate.abs(),
        );
        near(heat.to_f64(), expected_heat, expected_heat);
        assert!((expected_heat - diagonal_heat).abs() > 0.01 * expected_heat);
    }
}

#[test]
fn lorentz_lossless_windows_keep_reversible_work_and_honor_cancellation() {
    let pole = LorentzPole::try_new(0.7, 7.0, 0.0).unwrap();
    let frequencies = [4.0, 11.0];
    let window = HarmonicWindow::new(0.11, 0.037).unwrap();
    let electric = electric();
    let (fields, kinetic) = filtered(pole, frequencies, electric);
    let retained = (fields.clone(), kinetic.clone());
    let kernels = real_kernels(frequencies, window);
    let (heat, storage) =
        pole_moments(&kernels, &fields, Some(&kinetic), product, &mut || false).unwrap();
    let (_, current) = original_phasors(pole, frequencies, electric);
    let work = mean_product(window, frequencies, electric, current);
    assert_eq!(heat.to_f64(), 0.0);
    assert!(work.abs() > 1e-3);
    near(storage.to_f64(), work, work.abs());
    assert!(PoleTransfer::lorentz(pole, 7.0).is_err());

    let mut polls = 0;
    let cancelled = pole_moments(&kernels, &fields, Some(&kinetic), product, &mut || {
        polls += 1;
        polls == 2
    });
    assert!(matches!(cancelled, Err(MomentError::Cancelled)));
    assert_eq!(polls, 2);
    assert_eq!((fields.clone(), kinetic.clone()), retained);
    assert!(matches!(
        pole_moments(&kernels, &fields, Some(&kinetic[..1]), product, &mut || {
            false
        }),
        Err(MomentError::Unrepresentable)
    ));
}
