//! Full-matrix power admission for the common real modal wave reference.
//!
//! P = a + b, U = (a - b) / zc0 gives net mean input power
//! (||a||² - ||b||²) / (2 zc0), so every coherent input is passive exactly
//! when the largest singular value of R is at most one. Entrywise magnitudes
//! and eigenvalues cannot establish this for a coupled reflection operator.

use super::{C64, DuctError, Fft, FftC64, MmLineError};
use fs_la::eigen_complex::{ComplexSvdError, svd_complex_with_cancel};

const MAX_TAPS: usize = 65_536;
const POLL_TILE: usize = 256;

pub(super) fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), MmLineError> {
    if cancelled() {
        Err(MmLineError::Duct(DuctError::Cancelled {
            completed_steps: 0,
        }))
    } else {
        Ok(())
    }
}

/// Admit one complete reflection operator without clipping any mode.
pub(super) fn reflection_gain(
    reflection: &[C64],
    n: usize,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<f64, MmLineError> {
    poll(cancelled)?;
    if n == 0 || n > fs_duct::modal::MAX_MODES || reflection.len() != n * n {
        return Err(MmLineError::Invalid {
            what: "matrix reflection needs a supported square modal operator",
        });
    }
    if reflection
        .iter()
        .any(|z| !z.re.is_finite() || !z.im.is_finite())
    {
        return Err(MmLineError::Invalid {
            what: "matrix reflection contains a non-finite coefficient",
        });
    }
    let diagonal = reflection.iter().enumerate().all(|(entry, z)| {
        entry / n == entry % n || (z.re == 0.0 && z.im == 0.0)
    });
    let gain = if diagonal {
        // Uniform ducts often stay exactly diagonal: avoid a decomposition
        // and its small work allocations when the singular values are known.
        (0..n).map(|mode| reflection[mode * n + mode].abs()).fold(0.0, f64::max)
    } else {
        // The shared SVD works on columns directly rather than squaring the
        // condition number in R^H R. Nonconvergence is a refusal, never a
        // partially converged estimate admitted as passive.
        let decomposition = svd_complex_with_cancel(
            reflection,
            n,
            n,
            32.0 * f64::EPSILON,
            64,
            cancelled,
        )
        .map_err(|error| match error {
            ComplexSvdError::Cancelled => MmLineError::Duct(DuctError::Cancelled {
                completed_steps: 0,
            }),
            ComplexSvdError::BadInput { what } => MmLineError::Invalid { what },
            ComplexSvdError::NotConverged { .. } => MmLineError::Invalid {
                what: "matrix reflection singular-value admission did not converge",
            },
        })?;
        decomposition.singular_values[0]
    };
    if !gain.is_finite() || gain > 1.0 + crate::driving_point::PASSIVITY_TOL {
        return Err(MmLineError::Invalid {
            what: "matrix reflection operator amplifies a coherent modal input",
        });
    }
    Ok(gain)
}

/// Peak singular value on exactly the existing 4N-point realization grid.
///
/// Four modulated N-point FFTs produce F_4N[4m+q] without retaining 4N
/// samples for every entry. Only the nonnegative half of each interlaced
/// spectrum is kept; real taps make the other half conjugate. The transform,
/// modulation, spectra and assembled matrix buffers are reused. This is a
/// sampled-grid admission check, not a bound between sampled frequencies.
pub(super) fn realized_gain(
    taps: &[Vec<f64>],
    n: usize,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<f64, MmLineError> {
    poll(cancelled)?;
    if n == 0 || n > fs_duct::modal::MAX_MODES || taps.len() != n * n {
        return Err(MmLineError::Invalid {
            what: "matrix FIR needs a supported square modal operator",
        });
    }
    let length = taps[0].len();
    if length < 2
        || length > MAX_TAPS
        || !length.is_power_of_two()
        || taps.iter().any(|entry| entry.len() != length)
    {
        return Err(MmLineError::Invalid {
            what: "matrix FIR entries need one power-of-two length from 2 through 65536",
        });
    }
    let fft = Fft::new(length);
    let zero = FftC64::new(0.0, 0.0);
    let mut buffer = vec![zero; length];
    let mut scratch = vec![zero; length];
    let mut modulation = vec![zero; length];
    let mut spectra = vec![vec![zero; length / 2 + 1]; n * n];
    let mut reflection = vec![C64::ZERO; n * n];
    let mut peak = 0.0_f64;
    for offset in 0..4 {
        poll(cancelled)?;
        for (time, value) in modulation.iter_mut().enumerate() {
            if time.is_multiple_of(POLL_TILE) {
                poll(cancelled)?;
            }
            let phase = -core::f64::consts::TAU * offset as f64 * time as f64
                / (4 * length) as f64;
            let (sin, cos) = fs_math::det::sin_cos(phase);
            *value = FftC64::new(cos, sin);
        }
        // offset=0 includes DC and Nyquist; the other three grids stop
        // immediately before Nyquist. There are exactly 2N+1 bins in all.
        let bins = length / 2 + usize::from(offset == 0);
        for (entry, spectrum) in taps.iter().zip(&mut spectra) {
            poll(cancelled)?;
            for (time, ((value, phase), &tap)) in buffer
                .iter_mut()
                .zip(&modulation)
                .zip(entry)
                .enumerate()
            {
                if time.is_multiple_of(POLL_TILE) {
                    poll(cancelled)?;
                }
                if !tap.is_finite() {
                    return Err(MmLineError::Invalid {
                        what: "matrix FIR contains a non-finite tap",
                    });
                }
                *value = FftC64::new(tap * phase.re, tap * phase.im);
            }
            fft.forward(&mut buffer, &mut scratch);
            spectrum[..bins].copy_from_slice(&buffer[..bins]);
        }
        for bin in 0..bins {
            for (value, spectrum) in reflection.iter_mut().zip(&spectra) {
                *value = C64::new(spectrum[bin].re, spectrum[bin].im);
            }
            peak = peak.max(reflection_gain(&reflection, n, cancelled)?);
        }
    }
    Ok(peak)
}

/// Squared relative Frobenius change of the common causal prefix over ALL
/// entries. Exclude the preceding rung's fourth quarter, already checked by
/// the separate tail criterion. Scaling before squaring preserves very weak
/// finite matrices and prevents overflow; both norms use the same scale.
pub(super) fn prefix_change_energy(
    previous: Option<&[Vec<f64>]>,
    current: &[Vec<f64>],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<f64, MmLineError> {
    poll(cancelled)?;
    let Some(previous) = previous else {
        return Ok(f64::INFINITY);
    };
    if previous.is_empty() || previous.len() != current.len() {
        return Err(MmLineError::Invalid {
            what: "matrix prefix comparison needs matching nonempty operators",
        });
    }
    let old_length = previous[0].len();
    let prefix = 3 * (old_length / 4);
    if prefix == 0
        || previous.iter().any(|entry| entry.len() != old_length)
        || current.iter().any(|entry| entry.len() < prefix)
    {
        return Err(MmLineError::Invalid {
            what: "matrix prefix comparison needs a shared causal prefix",
        });
    }
    let mut scale = 0.0_f64;
    for (old, new) in previous.iter().zip(current) {
        for (old_tile, new_tile) in old[..prefix]
            .chunks(POLL_TILE)
            .zip(new[..prefix].chunks(POLL_TILE))
        {
            poll(cancelled)?;
            for (&a, &b) in old_tile.iter().zip(new_tile) {
                if !a.is_finite() || !b.is_finite() {
                    return Err(MmLineError::Invalid {
                        what: "matrix prefix contains a non-finite tap",
                    });
                }
                scale = scale.max(a.abs()).max(b.abs());
            }
        }
    }
    if scale == 0.0 {
        return Ok(0.0);
    }
    let (mut old_energy, mut new_energy, mut change) = (0.0_f64, 0.0_f64, 0.0_f64);
    for (old, new) in previous.iter().zip(current) {
        for (old_tile, new_tile) in old[..prefix]
            .chunks(POLL_TILE)
            .zip(new[..prefix].chunks(POLL_TILE))
        {
            poll(cancelled)?;
            for (&a, &b) in old_tile.iter().zip(new_tile) {
                let a = a / scale;
                let b = b / scale;
                old_energy += a * a;
                new_energy += b * b;
                change += (a - b) * (a - b);
            }
        }
    }
    Ok(change / old_energy.max(new_energy))
}

#[cfg(test)]
mod tests {
    use super::*;

    // H diag(a z^-1, b z^-3) H^T, H the orthogonal two-mode Hadamard
    // mixer. Its singular values are exactly a and b at EVERY frequency.
    fn mixed_delays(a: f64, b: f64) -> Vec<Vec<f64>> {
        let mut taps = vec![vec![0.0; 8]; 4];
        for (entry, ir) in taps.iter_mut().enumerate() {
            ir[1] = 0.5 * a;
            ir[3] = if entry == 0 || entry == 3 {
                0.5 * b
            } else {
                -0.5 * b
            };
        }
        taps
    }

    #[test]
    fn coherent_gain_is_not_an_entry_or_eigenvalue_bound() {
        let coupled = vec![C64::new(0.6, 0.0); 4];
        assert!(coupled.iter().all(|z| z.abs() < 1.0));
        // Each individually driven column even has norm sqrt(0.72) < 1.
        assert!(2.0 * 0.6_f64.powi(2) < 1.0);
        assert!(reflection_gain(&coupled, 2, &mut || false).is_err());
        // Eigenvalues alone also miss amplification in a nonnormal operator.
        let nonnormal = [C64::ZERO, C64::new(1.1, 0.0), C64::ZERO, C64::ZERO];
        assert!(reflection_gain(&nonnormal, 2, &mut || false).is_err());
        let lossless = [
            C64::ZERO,
            C64::new(0.0, 1.0),
            C64::new(0.0, -1.0),
            C64::ZERO,
        ];
        assert!((reflection_gain(&lossless, 2, &mut || false).unwrap() - 1.0).abs() < 1.0e-14);
    }

    #[test]
    fn coupled_fir_power_gain_checks_all_incoming_modes() {
        let passive = mixed_delays(0.95, 0.7);
        let peak = realized_gain(&passive, 2, &mut || false).unwrap();
        assert!((peak - 0.95).abs() < 2.0e-13);
        let active = mixed_delays(1.2, 0.3);
        // The old plane-only gate sees at most (1.2 + 0.3)/2 = 0.75.
        assert!(active[0].iter().map(|h| h.abs()).sum::<f64>() < 1.0);
        assert!(realized_gain(&active, 2, &mut || false).is_err());
        // An equal coherent two-mode unit-energy input excites the a branch.
        let unit = core::f64::consts::FRAC_1_SQRT_2;
        let returned_energy: f64 = (0..8)
            .flat_map(|time| {
                let active = &active;
                (0..2).map(move |row| {
                    let incoming = unit * (active[2 * row][time] + active[2 * row + 1][time]);
                    incoming * incoming
                })
            })
            .sum();
        assert!((returned_energy - 1.44).abs() < 1.0e-14);
    }

    #[test]
    fn interlaced_fft_matches_direct_oversampled_matrix_dtft() {
        let n = 3;
        let length = 16;
        let taps: Vec<Vec<f64>> = (0..n * n)
            .map(|entry| {
                (0..length)
                    .map(|time| (((11 * entry + 7 * time + 3) % 17) as f64 - 8.0) * 0.002)
                    .collect()
            })
            .collect();
        let mut direct = 0.0_f64;
        for bin in 0..=2 * length {
            let mut reflection = vec![C64::ZERO; n * n];
            for (value, entry) in reflection.iter_mut().zip(&taps) {
                for (time, &tap) in entry.iter().enumerate() {
                    let phase = -core::f64::consts::TAU * bin as f64 * time as f64
                        / (4 * length) as f64;
                    let (sin, cos) = fs_math::det::sin_cos(phase);
                    *value = *value + C64::new(cos * tap, sin * tap);
                }
            }
            direct = direct.max(reflection_gain(&reflection, n, &mut || false).unwrap());
        }
        let interlaced = realized_gain(&taps, n, &mut || false).unwrap();
        assert!((direct - interlaced).abs() < 2.0e-13);
    }

    #[test]
    fn prefix_convergence_detects_hidden_mode_change_at_every_scale() {
        for scale in [1.0e-200, 1.0, 1.0e200] {
            let mut previous = vec![vec![0.0; 8]; 4];
            previous[0][1] = scale * 0.5;
            previous[1][2] = scale * 0.25;
            let mut current: Vec<Vec<f64>> = previous
                .iter()
                .map(|entry| {
                    let mut longer = entry.clone();
                    longer.resize(16, 0.0);
                    longer
                })
                .collect();
            assert_eq!(
                prefix_change_energy(Some(&previous), &current, &mut || false).unwrap(),
                0.0
            );
            // R00 and all fourth-quarter tails are unchanged; an off-diagonal
            // causal reflection has not converged.
            current[1][2] = scale * 0.5;
            let change =
                prefix_change_energy(Some(&previous), &current, &mut || false).unwrap();
            assert!((change - 0.125).abs() < 1.0e-14);
        }
        let zero = vec![vec![0.0; 8]; 4];
        assert_eq!(prefix_change_energy(Some(&zero), &zero, &mut || false).unwrap(), 0.0);
        assert!(prefix_change_energy(None, &zero, &mut || false).unwrap().is_infinite());
    }

    #[test]
    fn fixed_rolloff_keeps_passive_scalar_prefix_convergence_meaningful() {
        let reference = crate::driving_point::reflectance_rolloff(65_536);
        let source = |length: usize| -> Vec<C64> {
            (0..=length / 2).map(|bin| {
                let phase = core::f64::consts::TAU * bin as f64 * 16.0 / length as f64;
                let (sin, cos) = fs_math::det::sin_cos(phase);
                C64::new(cos, sin).scale((0.7 - 1.0) / (0.7 + 1.0))
            }).collect()
        };
        let realize = |length: usize, fixed: bool| -> Vec<f64> {
            let spectrum = source(length);
            if fixed {
                let rolloff: Vec<_> = reference.iter().step_by(65_536 / length).copied().collect();
                crate::driving_point::realize_reflectance_ir_with_rolloff(
                    &spectrum, &rolloff, length, true,
                )
            } else {
                crate::driving_point::realize_reflectance_ir_from_positive_spectrum(
                    &spectrum, length, 16.0, true,
                )
            }
        };
        let previous = vec![realize(1024, true)];
        let current = vec![realize(2048, true)];
        let change = prefix_change_energy(Some(&previous), &current, &mut || false).unwrap();
        assert!(change < 1.0e-10);
        assert!(realized_gain(&current, 1, &mut || false).unwrap() < 0.18);
        assert!((current[0].iter().sum::<f64>() - (0.7 - 1.0) / (0.7 + 1.0)).abs() < 2.0e-13);

        // Recomputing the phase taper changes the numerical target between
        // rungs, even though this physical response is just a passive delay.
        let old_previous = vec![realize(1024, false)];
        let old_current = vec![realize(2048, false)];
        assert!(
            prefix_change_energy(Some(&old_previous), &old_current, &mut || false).unwrap()
                > 1.0e-8
        );
    }

    #[test]
    fn admitted_coupled_fir_push_preserves_coherent_modal_power_and_delay() {
        use super::super::{MatrixFirLine, MatrixLineReceipt};
        for sign in [1.0, -1.0] {
            let taps = mixed_delays(0.95, 0.7);
            let peak = realized_gain(&taps, 2, &mut || false).unwrap();
            // Synthetic exact FIR control for the actual runtime convolution;
            // it is not a duct-realization accuracy receipt.
            let mut line = MatrixFirLine {
                n_modes: 2,
                n_fft: 8,
                taps,
                history: vec![vec![0.0; 8]; 2],
                write: 0,
                incoming: vec![0.0; 2],
                energy: vec![0.0; 2],
                zc0: 2.0,
                receipt: MatrixLineReceipt {
                    n_fft: 8,
                    rungs_tried: 2,
                    tail_energy: 0.0,
                    entry_tail_energy: vec![0.0; 4],
                    prefix_change_energy: 0.0,
                    peak_shift_cents: 0.0,
                    peak_reflectance: 0.825,
                    peak_operator_norm: peak,
                },
            };
            let unit = core::f64::consts::FRAC_1_SQRT_2;
            let (delay, amplitude) = if sign > 0.0 { (1, 0.95) } else { (3, 0.7) };
            for time in 0..32 {
                let outgoing = if time == 0 { [unit, sign * unit] } else { [0.0; 2] };
                let incoming = line.push(&outgoing).unwrap();
                let expected = if time == delay { amplitude * unit } else { 0.0 };
                assert!((incoming[0] - expected).abs() < 1.0e-14);
                assert!((incoming[1] - sign * expected).abs() < 1.0e-14);
            }
            let returned: f64 = line.mode_energy().iter().sum();
            assert!((returned - amplitude * amplitude).abs() < 1.0e-14);
            assert!(returned <= 1.0);
        }
    }

    #[test]
    fn admission_rejects_nonfinite_data_and_honors_cancellation() {
        let mut taps = mixed_delays(0.9, 0.7);
        taps[3][4] = f64::NAN;
        assert!(realized_gain(&taps, 2, &mut || false).is_err());
        assert!(matches!(
            realized_gain(&mixed_delays(0.9, 0.7), 2, &mut || true),
            Err(MmLineError::Duct(DuctError::Cancelled { .. }))
        ));
        let mut polls = 0;
        assert!(matches!(
            reflection_gain(&[C64::new(0.25, 0.0); 4], 2, &mut || {
                polls += 1;
                polls >= 3
            }),
            Err(MmLineError::Duct(DuctError::Cancelled { .. }))
        ));
    }
}
