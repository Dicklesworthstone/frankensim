//! Positive original receiver/reservoir contraction; no new acoustic solver.
use super::*;
use fs_bem::observation_transfers::HarmonicField;
use fs_fft::harmonic_window::{
    AngularHarmonicSpectrumTerm, HarmonicProductTone, HarmonicSpectrumTermLimits,
};
use fs_math::{dd::Dd, det};

fn invalid(reason: &'static str) -> MolecularWindowError {
    MolecularWindowError::Invalid(reason)
}
fn budget() -> MolecularWindowError {
    invalid("coherent molecular preparation budget")
}
fn finite(z: C64) -> bool {
    z.re.is_finite() && z.im.is_finite()
}

pub(super) fn spectrum(
    binding: &MolecularGasSpectrumHeating<'_>,
    rates: &[&[C64]],
    clock: &str,
    limits: MolecularWindowLimits,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(PreparedHarmonicProductSpectrum, usize), MolecularWindowError> {
    let banks = binding.absorptions();
    let bins = banks.len();
    let first = &banks[0];
    let receivers = first.receivers().len();
    if limits.maximum_bins == 0
        || limits.maximum_bins > 64
        || bins > limits.maximum_bins
        || limits.maximum_receivers == 0
        || limits.maximum_receivers > 65_536
        || receivers > limits.maximum_receivers
        || limits.maximum_products == 0
        || limits.maximum_products > 16_777_216
        || limits.maximum_workspace_bytes == 0
        || limits.maximum_workspace_bytes > 64 * 1024 * 1024
    {
        return Err(budget());
    }
    if clock.trim().is_empty()
        || clock.len() > 4096
        || rates.len() != bins
        || !limits.maximum_constitutive_relative_error.is_finite()
        || limits.maximum_constitutive_relative_error <= 0.0
        || limits.maximum_constitutive_relative_error >= 1.0
    {
        return Err(invalid(
            "original coherent rates, clock and positive loss tolerance required",
        ));
    }
    let mut rate_entries = 0usize;
    for (bank, &row) in banks.iter().zip(rates) {
        poll(cancelled)?;
        if bank.receivers() != first.receivers() || bank.quadrature() != first.quadrature() {
            return Err(invalid(
                "coherent frequencies require identical original receiver/volume geometry",
            ));
        }
        if row.len() != bank.sources().len() || row.iter().any(|&z| !finite(z)) {
            return Err(invalid(
                "coherent peak rates must match each original finite source bank",
            ));
        }
        rate_entries = rate_entries.checked_add(row.len()).ok_or_else(budget)?;
    }
    let pairs = bins * (bins + 1) / 2;
    let terms = 2 * pairs;
    // Native material construction has its own 32-reservoir hard bound. Use
    // this upper bound BEFORE any new field/product workspace allocation.
    let max_arms = first.source().processes().len();
    let products = rate_entries
        .checked_mul(4)
        .and_then(|n| {
            (2 * bins + 4 * pairs)
                .checked_mul(max_arms)
                .and_then(|m| n.checked_add(m))
        })
        .and_then(|n| n.checked_mul(receivers))
        .ok_or_else(budget)?;
    if products > limits.maximum_products {
        return Err(budget());
    }
    let native_bytes = terms * core::mem::size_of::<HarmonicProductTone>();
    let mut bytes = clock.len();
    for (count, size) in [
        (rate_entries, core::mem::size_of::<C64>()),
        (bins, core::mem::size_of::<Vec<C64>>()),
        (
            bins,
            core::mem::size_of::<super::super::spectrum::MolecularGasSpectrumBinPower<'_>>()
                + core::mem::size_of::<usize>(),
        ),
        (bins * receivers, core::mem::size_of::<C64>()),
        (bins * max_arms, core::mem::size_of::<C64>()),
        (
            receivers,
            core::mem::size_of::<HarmonicField>() + core::mem::size_of::<Dd>(),
        ),
        (bins, core::mem::size_of::<[Dd; 2]>()),
        (pairs, core::mem::size_of::<[Dd; 4]>()),
        (terms, core::mem::size_of::<AngularHarmonicSpectrumTerm>()),
        (1, native_bytes),
    ] {
        bytes = count
            .checked_mul(size)
            .and_then(|n| bytes.checked_add(n))
            .ok_or_else(budget)?;
    }
    for bank in banks {
        bytes = bytes
            .checked_add(bank.workspace_bytes())
            .ok_or_else(budget)?;
    }
    if bytes > limits.maximum_workspace_bytes {
        return Err(budget());
    }

    let bulk = first.source().bulk_relaxation()?;
    poll(cancelled)?;
    let arms = bulk.arms();
    let mut coefficients = vec![C64::ZERO; bins * arms.len()];
    for (bin, bank) in banks.iter().enumerate() {
        poll(cancelled)?;
        let omega = bank.angular_frequency_rad_s();
        let modulus = bank.source().caloric_response(omega)?.bulk_modulus_pa;
        let compliance = modulus.recip();
        let root_loss = det::sqrt(omega) * det::sqrt(compliance.im);
        if !root_loss.is_finite() || (arms.is_empty() != (root_loss == 0.0)) {
            return Err(invalid(
                "original molecular heat coefficient is unrepresentable",
            ));
        }
        let mut relative_norm = Dd::ZERO;
        for (arm_index, arm) in arms.iter().enumerate() {
            // r_j/sqrt(delta_K_j*tau_j), per original pressure. Separate
            // roots avoid requiring delta_K*tau to fit before normalization.
            let coefficient = (arm.response(omega)? / modulus)
                .scale(1.0 / det::sqrt(arm.delta_bulk_modulus_pa))
                .scale(1.0 / det::sqrt(arm.relaxation_time_s));
            if !finite(coefficient) || coefficient == C64::ZERO {
                return Err(invalid(
                    "molecular reservoir heat response is unrepresentable",
                ));
            }
            coefficients[bin * arms.len() + arm_index] = coefficient;
            let re = Dd::from_f64(coefficient.re / root_loss);
            let im = Dd::from_f64(coefficient.im / root_loss);
            relative_norm = relative_norm + re * re + im * im;
        }
        let error = if arms.is_empty() {
            0.0
        } else {
            (relative_norm - Dd::ONE).to_f64().abs()
        };
        if !error.is_finite() || error > limits.maximum_constitutive_relative_error {
            return Err(MolecularWindowError::ConstitutiveLoss {
                bin,
                relative_error: error,
                maximum: limits.maximum_constitutive_relative_error,
            });
        }
    }

    let mut pressures = vec![C64::ZERO; bins * receivers];
    let mut fields = vec![HarmonicField::default(); receivers];
    for (bin, (bank, &row)) in banks.iter().zip(rates).enumerate() {
        poll(cancelled)?;
        bank.response()
            .harmonic()
            .evaluate_coherent(row, &mut fields, &mut *cancelled)?;
        for (receiver, field) in fields.iter().enumerate() {
            pressures[bin * receivers + receiver] = field.pressure_pa;
        }
    }
    let mut weights = vec![Dd::ZERO; receivers];
    for point in first.quadrature() {
        weights[point.receiver] = weights[point.receiver] + Dd::from_f64(point.volume_m3);
    }
    let mut moments = vec![[Dd::ZERO; 4]; pairs];
    let mut amplitudes = vec![[Dd::ZERO; 2]; bins];
    for (receiver, weight) in weights.iter().enumerate() {
        poll(cancelled)?;
        let volume = weight.to_f64();
        if volume == 0.0 {
            continue;
        }
        let root_volume = det::sqrt(volume);
        for arm in 0..arms.len() {
            poll(cancelled)?;
            for (bin, amplitude) in amplitudes.iter_mut().enumerate() {
                let coefficient = coefficients[bin * arms.len() + arm].scale(root_volume);
                if !finite(coefficient) || coefficient == C64::ZERO {
                    return Err(invalid(
                        "volume-weighted molecular heat response is unrepresentable",
                    ));
                }
                let p = pressures[bin * receivers + receiver];
                let [pr, pi, cr, ci] =
                    [p.re, p.im, coefficient.re, coefficient.im].map(Dd::from_f64);
                *amplitude = [pr * cr - pi * ci, pr * ci + pi * cr];
                let squared = amplitude[0] * amplitude[0] + amplitude[1] * amplitude[1];
                if p != C64::ZERO && squared.to_f64() == 0.0 {
                    return Err(invalid("nonzero molecular reservoir heat underflowed"));
                }
            }
            let mut pair = 0;
            for (i, a) in amplitudes.iter().enumerate() {
                for (j, b) in amplitudes[..=i].iter().enumerate() {
                    let factor = Dd::from_f64(if i == j { 0.5 } else { 1.0 });
                    let [rr, ii, ri, ir] = [a[0] * b[0], a[1] * b[1], a[0] * b[1], a[1] * b[0]];
                    for (sum, value) in
                        moments[pair]
                            .iter_mut()
                            .zip([rr + ii, ir - ri, rr - ii, ri + ir])
                    {
                        *sum = *sum + factor * value;
                    }
                    pair += 1;
                }
            }
        }
    }
    let mut output = Vec::with_capacity(terms);
    let mut pair = 0;
    for (i, a) in banks.iter().enumerate() {
        poll(cancelled)?;
        let wa = Dd::from_f64(a.angular_frequency_rad_s());
        for b in &banks[..=i] {
            let wb = Dd::from_f64(b.angular_frequency_rad_s());
            let [dr, di, sr, si] = moments[pair];
            for (omega, amplitude) in [(wa - wb, [dr, di]), (wa + wb, [sr, si])] {
                if amplitude.iter().any(|v| v.hi != 0.0 || v.lo != 0.0) {
                    output.push(AngularHarmonicSpectrumTerm {
                        angular_frequency_rad_s: omega,
                        peak_amplitude: amplitude,
                    });
                }
            }
            pair += 1;
        }
    }
    let spectrum = PreparedHarmonicProductSpectrum::from_angular_terms_cancellable(
        output,
        HarmonicSpectrumTermLimits {
            maximum_terms: terms,
            maximum_workspace_bytes: native_bytes,
        },
        &mut *cancelled,
    )?;
    poll(cancelled)?;
    Ok((spectrum, bytes))
}
