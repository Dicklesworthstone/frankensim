//! Positive spatial contraction into the original exact temporal spectrum.
use super::*;
use fs_fft::harmonic_window::{
    AngularHarmonicSpectrumTerm, HarmonicProductTone, HarmonicSpectrumTermLimits,
};
use fs_math::det;

const MAX_LINES: usize = 64;
const MAX_POINTS: usize = 1_048_576;
const MAX_PAIRS: usize = 16_777_216;
const MAX_WORKSPACE: usize = 64 * 1024 * 1024;
// Positive four-point Gauss-Legendre rule on [-1,1]. Panel subdivision is an
// explicit caller resolution axis, independent of all thermal step sizes.
const GAUSS: [(f64, f64); 4] = [
    (-0.861_136_311_594_052_6, 0.347_854_845_137_453_8),
    (-0.339_981_043_584_856_3, 0.652_145_154_862_546_1),
    (0.339_981_043_584_856_3, 0.652_145_154_862_546_1),
    (0.861_136_311_594_052_6, 0.347_854_845_137_453_8),
];

pub(super) fn admit(
    responses: &[BulkTractionAbsorptionResponse],
    clock: &str,
    cells: usize,
    limits: BulkViscousWindowLimits,
) -> Result<(), BulkViscousWindowError> {
    if responses.is_empty()
        || responses.len() > limits.maximum_lines
        || limits.maximum_lines == 0
        || limits.maximum_lines > MAX_LINES
        || cells == 0
        || cells > limits.maximum_cells
        || limits.maximum_cells == 0
        || limits.maximum_cells > super::super::MAX_BULK_HEATING_CELLS
        || limits.maximum_spatial_points == 0
        || limits.maximum_spatial_points > MAX_POINTS
        || limits.maximum_pair_samples == 0
        || limits.maximum_pair_samples > MAX_PAIRS
        || limits.maximum_workspace_bytes == 0
        || limits.maximum_workspace_bytes > MAX_WORKSPACE
    {
        return Err(budget());
    }
    if clock.trim().is_empty()
        || !limits.maximum_panel_exponent.is_finite()
        || limits.maximum_panel_exponent <= 0.0
        || limits.maximum_panel_exponent > 4.0
        || !limits.maximum_diagonal_relative_error.is_finite()
        || limits.maximum_diagonal_relative_error <= 0.0
        || limits.maximum_diagonal_relative_error >= 1.0
    {
        return Err(invalid(
            "source clock and finite positive spatial accuracy controls required",
        ));
    }
    let first = &responses[0];
    for (index, response) in responses.iter().enumerate() {
        if response.equilibrium() != first.equilibrium()
            || response.area_m2() != first.area_m2()
            || response.length_m() != first.length_m()
        {
            return Err(invalid(
                "all coherent tones require the original cylinder and complete equilibrium",
            ));
        }
        for previous in &responses[..index] {
            if previous.angular_frequency_rad_s() == response.angular_frequency_rad_s() {
                let a = previous.response();
                let b = response.response();
                if [
                    a.inlet_traction_impedance_pa_s_m3,
                    a.inlet_pressure_over_flow_pa_s_m3,
                    a.outlet_traction_over_inlet_flow_pa_s_m3,
                    a.outlet_pressure_over_inlet_flow_pa_s_m3,
                    a.outlet_over_inlet_flow,
                ] != [
                    b.inlet_traction_impedance_pa_s_m3,
                    b.inlet_pressure_over_flow_pa_s_m3,
                    b.outlet_traction_over_inlet_flow_pa_s_m3,
                    b.outlet_pressure_over_inlet_flow_pa_s_m3,
                    b.outlet_over_inlet_flow,
                ] {
                    return Err(invalid(
                        "equal-frequency tones require the same original loaded unit response",
                    ));
                }
            }
        }
    }
    Ok(())
}

pub(super) fn cells(
    responses: &[BulkTractionAbsorptionResponse],
    binding: &BulkViscousHeating<'_>,
    limits: BulkViscousWindowLimits,
    checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
) -> Result<Vec<BulkViscousCellSpectrum>, BulkViscousWindowError> {
    let mut indices: Vec<_> = (0..responses.len()).collect();
    indices.sort_unstable_by(|&a, &b| {
        responses[a]
            .angular_frequency_rad_s()
            .total_cmp(&responses[b].angular_frequency_rad_s())
            .then_with(|| a.cmp(&b))
    });
    let mut frequencies = Vec::new();
    let mut group_of = vec![0; responses.len()];
    let mut maximum_rate: f64 = 0.0;
    for &index in &indices {
        checkpoint()?;
        let omega = responses[index].angular_frequency_rad_s();
        if frequencies.last() != Some(&omega) {
            frequencies.push(omega);
        }
        group_of[index] = frequencies.len() - 1;
        let wave = responses[index]
            .equilibrium()
            .absorption()
            .thermoviscous_wave(omega)?;
        maximum_rate = maximum_rate.max(wave.wavenumber_rad_m().abs());
    }
    let pairs = frequencies.len() * (frequencies.len() + 1) / 2;
    let terms = 2 * pairs;
    let mut points = 0usize;
    let mut panels = Vec::with_capacity(binding.cells.len());
    for cell in &binding.cells {
        let count = (maximum_rate * (cell.end_distance_m - cell.start_distance_m)
            / limits.maximum_panel_exponent)
            .ceil()
            .max(1.0);
        if !count.is_finite() || count > (limits.maximum_spatial_points / GAUSS.len()) as f64 {
            return Err(budget());
        }
        let count = count as usize;
        points = points.checked_add(count * GAUSS.len()).ok_or_else(budget)?;
        panels.push(count);
    }
    if points > limits.maximum_spatial_points
        || points.checked_mul(pairs).ok_or_else(budget)? > limits.maximum_pair_samples
    {
        return Err(budget());
    }
    let native_bytes = terms * core::mem::size_of::<HarmonicProductTone>();
    let retained = binding
        .cells
        .len()
        .checked_mul(
            native_bytes
                + core::mem::size_of::<BulkViscousCellSpectrum>()
                + core::mem::size_of::<super::super::BulkViscousCellHeating>()
                + core::mem::size_of::<(CellId, f64)>(),
        )
        .ok_or_else(budget)?;
    let scratch = native_bytes
        + terms * core::mem::size_of::<AngularHarmonicSpectrumTerm>()
        + pairs * core::mem::size_of::<[Dd; 4]>()
        + responses.len() * (core::mem::size_of::<Dd>() + 2 * core::mem::size_of::<usize>())
        + frequencies.len() * (core::mem::size_of::<[Dd; 2]>() + core::mem::size_of::<f64>())
        + panels.len() * core::mem::size_of::<usize>();
    if retained.checked_add(scratch).ok_or_else(budget)? > limits.maximum_workspace_bytes {
        return Err(budget());
    }
    let material = responses[0].equilibrium().absorption();
    let eta = (4.0 / 3.0) * material.state().transport().dynamic_viscosity_pa_s()
        + material.bulk_viscosity_pa_s();
    let root_eta = det::sqrt(eta);
    let mut result = Vec::with_capacity(binding.cells.len());
    for (cell, panels) in binding.cells.iter().zip(panels) {
        checkpoint()?;
        let mut moments = vec![[Dd::ZERO; 4]; pairs];
        let mut diagonal = vec![Dd::ZERO; responses.len()];
        let mut amplitudes = vec![[Dd::ZERO; 2]; frequencies.len()];
        let width = (cell.end_distance_m - cell.start_distance_m) / panels as f64;
        for panel in 0..panels {
            let left = cell.start_distance_m + panel as f64 * width;
            let right = if panel + 1 == panels {
                cell.end_distance_m
            } else {
                cell.start_distance_m + (panel + 1) as f64 * width
            };
            if right <= left {
                return Err(DepositionError::Arithmetic.into());
            }
            for (node, weight) in GAUSS {
                checkpoint()?;
                let x = left + 0.5 * (right - left) * (1.0 + node);
                let volume = 0.5 * (right - left) * weight * responses[0].area_m2();
                if !volume.is_finite() || volume <= 0.0 {
                    return Err(DepositionError::Arithmetic.into());
                }
                // Apply physical roots before products; no unweighted square
                // is required to fit in binary64. Coherent fields combine first.
                let scale = root_eta * det::sqrt(volume);
                if !scale.is_finite() || (eta > 0.0 && scale == 0.0) {
                    return Err(DepositionError::Arithmetic.into());
                }
                amplitudes.fill([Dd::ZERO; 2]);
                for &index in &indices {
                    let gradient = responses[index].at_distance(x)?.velocity_gradient_per_s;
                    let raw = [Dd::from_f64(gradient.re), Dd::from_f64(gradient.im)];
                    for (sum, value) in amplitudes[group_of[index]].iter_mut().zip(raw) {
                        *sum = *sum + value;
                    }
                    let weighted = raw.map(|v| v * Dd::from_f64(scale));
                    diagonal[index] = diagonal[index]
                        + Dd::from_f64(0.5)
                            * (weighted[0] * weighted[0] + weighted[1] * weighted[1]);
                }
                for value in &mut amplitudes {
                    for part in value {
                        *part = *part * Dd::from_f64(scale);
                    }
                }
                let mut pair = 0;
                for (i, a) in amplitudes.iter().enumerate() {
                    for (j, b) in amplitudes[..=i].iter().enumerate() {
                        let factor = Dd::from_f64(if i == j { 0.5 } else { 1.0 });
                        let rr = a[0] * b[0];
                        let ii = a[1] * b[1];
                        let ri = a[0] * b[1];
                        let ir = a[1] * b[0];
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
        let mut maximum_error: f64 = 0.0;
        for (line, (response, integral)) in responses.iter().zip(diagonal).enumerate() {
            checkpoint()?;
            let exact = response
                .absorption_between(cell.start_distance_m, cell.end_distance_m)?
                .viscous_heat_w;
            let difference = (integral - Dd::from_f64(exact)).to_f64().abs();
            let relative_error = if exact == 0.0 && difference == 0.0 {
                0.0
            } else {
                difference / exact
            };
            if !relative_error.is_finite()
                || relative_error > limits.maximum_diagonal_relative_error
            {
                return Err(BulkViscousWindowError::SpatialResolution {
                    cell: cell.cell,
                    line,
                    relative_error,
                    maximum: limits.maximum_diagonal_relative_error,
                });
            }
            maximum_error = maximum_error.max(relative_error);
        }
        let mut output = Vec::with_capacity(terms);
        let mut pair = 0;
        for (i, &omega_a) in frequencies.iter().enumerate() {
            for &omega_b in &frequencies[..=i] {
                let [dr, di, sr, si] = moments[pair];
                for (omega, amplitude) in [
                    (Dd::from_f64(omega_a) - Dd::from_f64(omega_b), [dr, di]),
                    (Dd::from_f64(omega_a) + Dd::from_f64(omega_b), [sr, si]),
                ] {
                    // Exact zeros need no phase calculation. This is not a
                    // small-amplitude cutoff; every retained low part survives.
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
        let spectrum = spectral_call(checkpoint, |cancelled| {
            PreparedHarmonicProductSpectrum::from_angular_terms_cancellable(
                output,
                HarmonicSpectrumTermLimits {
                    maximum_terms: terms,
                    maximum_workspace_bytes: native_bytes,
                },
                cancelled,
            )
        })?;
        result.push(BulkViscousCellSpectrum {
            cell: cell.cell,
            start_distance_m: cell.start_distance_m,
            end_distance_m: cell.end_distance_m,
            spatial_points: panels * GAUSS.len(),
            maximum_diagonal_relative_error: maximum_error,
            spectrum,
        });
    }
    Ok(result)
}
