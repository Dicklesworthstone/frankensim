//! Broadband body-frame radiation filters from solver-neutral samples.
//!
//! This module is deliberately an offline bridge: a radiation solver supplies
//! complex spherical-harmonic transfer samples plus direct directional
//! held-out samples, complex harmonics are changed exactly into a real
//! tesseral channel basis, and each channel is fit and Tustin-discretized into
//! a causal real state-space filter. The producer may be BEM, an analytic
//! model, or a future solver; this L3 composition layer does not depend on it.
//! It does not move the source, add propagation delay, or synthesize listener
//! pressure.

use fs_math::c64::C64;
use fs_math::det;
use fs_vfit::discretize::{DiscreteStateSpace, DiscretizeError, bilinear_state_space};
use fs_vfit::{FitOptions, VfError, WeightPreset, vector_fit};

/// Maximum frequencies in either offline solve grid.
pub const MAX_BROADBAND_FREQUENCIES: usize = 512;
/// Maximum input-channel-frequency cells retained for fitting.
pub const MAX_BROADBAND_TRANSFER_SAMPLES: usize = 4_194_304;
/// Maximum scalar inputs in one neutral sample bundle.
pub const MAX_BROADBAND_INPUTS: usize = 256;
/// Maximum independent directions in the held-out far-field gate.
pub const MAX_VALIDATION_DIRECTIONS: usize = 4096;
/// Maximum SH degree supported by the normalized recurrence in this module.
pub const MAX_RADIATION_SH_DEGREE: usize = 64;
/// Evidence authority of the produced artifact.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BroadbandRadiationAuthority {
    /// Numerically gated estimate, not an enclosure or physical certificate.
    EstimateOnly,
}

/// Harmonic sign carried by every neutral sample bundle.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum HarmonicTimeConvention {
    /// Physical fields use `Re(H(omega) exp(-i omega t))`.
    ExpNegativeIOmegaT,
}

/// Worst-case solver diagnostics for one frequency's complete input batch.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct RadiationSampleDiagnostics {
    /// Fraction of direct far-field power captured by the SH truncation.
    pub captured_fraction: f64,
    /// Producer-reported surface panels per acoustic wavelength.
    pub panels_per_wavelength: f64,
    /// Producer-reported condition-number lower bound.
    pub condition_lower_bound: f64,
}

/// One training frequency in canonical complex-SH coefficient ordering.
#[derive(Debug, Clone, PartialEq)]
pub struct ComplexShTrainingSample {
    /// Physical angular frequency [rad/s].
    pub omega_rad_s: f64,
    /// Input-major `a_lm`, indexed `l(l+1)+m`, length `(l_max+1)^2`.
    pub coefficients_by_input: Vec<Vec<C64>>,
    /// Worst-case producer diagnostics across this frequency's input batch.
    pub diagnostics: RadiationSampleDiagnostics,
}

/// One disjoint held-out frequency with direct far-field references.
#[derive(Debug, Clone, PartialEq)]
pub struct DirectFarFieldHeldOutSample {
    /// Physical angular frequency [rad/s].
    pub omega_rad_s: f64,
    /// Nonzero finite body-frame directions, independent of SH projection.
    pub directions: Vec<[f64; 3]>,
    /// Input-major direct complex far-field amplitudes, aligned to directions.
    pub far_field_by_input: Vec<Vec<C64>>,
    /// Worst-case producer diagnostics across this frequency's input batch.
    pub diagnostics: RadiationSampleDiagnostics,
}

/// Radiation-solver-neutral offline input to the fitting bridge.
#[derive(Debug, Clone, PartialEq)]
pub struct SampledRadiationData {
    /// Stable producer/artifact identity retained in the report.
    pub source_id: String,
    /// Input harmonic sign; currently the one fs-bem and this bridge document.
    pub harmonic_time_convention: HarmonicTimeConvention,
    /// Highest degree present in every complex-SH training row.
    pub l_max: usize,
    /// Stable scalar-input transfer identities in the common row order. The
    /// source identity binds their units and producer semantics.
    pub input_ids: Vec<String>,
    /// Strictly increasing physical training grid.
    pub training: Vec<ComplexShTrainingSample>,
    /// Strictly increasing, disjoint physical withheld grid.
    pub held_out: Vec<DirectFarFieldHeldOutSample>,
}

/// Offline fit, discretization, and withheld-validation controls.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct BroadbandRadiationControls {
    /// Runtime sample rate [Hz].
    pub sample_rate_hz: f64,
    /// Minimum accepted producer-reported SH captured fraction.
    pub minimum_captured_fraction: f64,
    /// Vector-fit pole count per real tesseral channel.
    pub fit_order: usize,
    /// Vector-fit pole-relocation iteration cap.
    pub fit_iterations: usize,
    /// Vector-fit row weighting.
    pub fit_weights: WeightPreset,
    /// Absolute far-field amplitude floor in pointwise relative errors.
    pub far_field_signal_floor: f64,
    /// Maximum accepted pointwise normalized complex error.
    pub maximum_normalized_error: f64,
    /// Maximum accepted RMS normalized complex error.
    pub rms_normalized_error: f64,
}

/// One canonical real tesseral channel.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct RealTesseralChannel {
    /// Spherical-harmonic degree.
    pub l: usize,
    /// Signed order: zero is zonal, `+m` is `sqrt(2) Re(Y_lm)`, and
    /// `-m` is `sqrt(2) Im(Y_lm)`.
    pub signed_m: i64,
}

/// Filters for one caller-identified scalar input, in artifact channel order.
#[derive(Debug, Clone, PartialEq)]
pub struct BroadbandRadiationInputArtifact {
    /// Stable scalar-input transfer identity.
    pub id: String,
    /// Scalar-input to body-frame far-field real-SH coefficients.
    pub filters: Vec<DiscreteStateSpace>,
}

/// Withheld evidence and approximation boundaries for an artifact.
#[derive(Debug, Clone, PartialEq)]
pub struct BroadbandRadiationReport {
    /// This slice can publish estimates only.
    pub authority: BroadbandRadiationAuthority,
    /// Stable identity of the neutral radiation sample source.
    pub source_id: String,
    /// Harmonic sign of the accepted source samples.
    pub harmonic_time_convention: HarmonicTimeConvention,
    /// Exact physical training angular-frequency grid [rad/s].
    pub training_omega_rad_s: Vec<f64>,
    /// Exact disjoint held-out physical angular-frequency grid [rad/s].
    pub held_out_omega_rad_s: Vec<f64>,
    /// Maximum pointwise normalized complex far-field error.
    pub maximum_normalized_complex_error: f64,
    /// RMS pointwise normalized complex far-field error.
    pub rms_normalized_complex_error: f64,
    /// Complete admitted fit, budget, and error controls.
    pub controls: BroadbandRadiationControls,
    /// Worst-case producer diagnostics across every accepted sample.
    pub sampling_diagnostics: RadiationSampleDiagnostics,
    /// Worst normalized training residual reported by any channel fit.
    pub maximum_channel_weighted_rms: f64,
    /// Total direct input-frequency-direction comparisons in the held-out gate.
    pub held_out_comparison_count: usize,
    /// Explicit approximation boundary carried with every artifact.
    pub no_claims: &'static str,
}

/// Offline body-frame far-field scalar-input filter bank.
#[derive(Debug, Clone, PartialEq)]
pub struct BroadbandRadiationArtifact {
    /// Sample interval shared by every discrete realization [s].
    pub sample_interval_s: f64,
    /// Highest retained SH degree.
    pub l_max: usize,
    /// Canonical channel metadata; each input's filters match this order.
    pub channels: Vec<RealTesseralChannel>,
    /// One filter bank per scalar input, in caller order.
    pub inputs: Vec<BroadbandRadiationInputArtifact>,
    /// Fit and withheld evidence.
    pub report: BroadbandRadiationReport,
}

/// Stable refusal surface for artifact construction and basis conversion.
#[derive(Debug)]
pub enum BroadbandRadiationError {
    /// Caller admission failure carrying a stable rule description.
    InvalidInput(&'static str),
    /// Vector-fit failure: `(input index, real-SH channel, source)`.
    Fit(usize, usize, VfError),
    /// Tustin failure: `(input index, real-SH channel, source)`.
    Discretize(usize, usize, DiscretizeError),
    /// SH capture failure: `(measured minimum, required minimum)`.
    CapturedFraction(f64, f64),
    /// Direct held-out failure: `(maximum, rms, maximum limit, rms limit)`.
    HeldOutError(f64, f64, f64, f64),
    /// Failed properness, stability, or evaluation invariant.
    NumericalFailure(&'static str),
}

/// Honest boundary attached verbatim to every accepted report.
pub const BROADBAND_RADIATION_NO_CLAIMS: &str = "Estimate-only sampled linear exterior-acoustic transfer about an undeformed stationary reference body. The artifact maps each source-defined scalar input to body-frame far-field amplitude before 1/r propagation; source_id, not this bridge, owns the input units and producer semantics. It does not enclose source-transfer, BEM, spherical-harmonic, vector-fit, discretization, or above-band errors beyond reported sampled metrics; and it does not model moving-surface/FW-H terms, deformation-updated geometry or normals, fluid loading back-coupling, propagation delay, Doppler/convective effects, near-field pressure, impacts, contact or air-film sound, supports, rooms, heads, or nonlinear acoustics.";

/// Canonical channel order: degree-major, then zonal, cosine/sine pairs.
fn real_tesseral_channels(l_max: usize) -> Vec<RealTesseralChannel> {
    let mut out = Vec::with_capacity((l_max + 1) * (l_max + 1));
    for l in 0..=l_max {
        out.push(RealTesseralChannel { l, signed_m: 0 });
        for m in 1..=l {
            let m = i64::try_from(m).expect("degree is bounded");
            out.push(RealTesseralChannel { l, signed_m: m });
            out.push(RealTesseralChannel { l, signed_m: -m });
        }
    }
    out
}

/// Exact unitary change from canonical complex SH coefficients to real
/// tesseral coefficients. Coefficient values remain complex frequency
/// responses.
pub fn complex_sh_to_real_tesseral(
    l_max: usize,
    coefficients: &[C64],
) -> Result<Vec<C64>, BroadbandRadiationError> {
    validate_coefficients(l_max, coefficients, "complex SH")?;
    let inv_sqrt_two = 1.0 / det::sqrt(2.0);
    let mut out = Vec::with_capacity(coefficients.len());
    for l in 0..=l_max {
        out.push(coefficients[complex_sh_index(l, 0)]);
        for m in 1..=l {
            let mi = i64::try_from(m).expect("m is bounded by MAX_RADIATION_SH_DEGREE");
            let positive = coefficients[complex_sh_index(l, mi)];
            let signed_negative =
                coefficients[complex_sh_index(l, -mi)].scale(if m % 2 == 0 { 1.0 } else { -1.0 });
            out.push((positive + signed_negative).scale(inv_sqrt_two));
            let difference = positive - signed_negative;
            out.push(C64::new(-difference.im, difference.re).scale(inv_sqrt_two));
        }
    }
    Ok(out)
}

/// Inverse unitary change into canonical `l(l+1)+m` complex ordering.
pub fn real_tesseral_to_complex_sh(
    l_max: usize,
    coefficients: &[C64],
) -> Result<Vec<C64>, BroadbandRadiationError> {
    validate_coefficients(l_max, coefficients, "real tesseral")?;
    let count = coefficients.len();
    let inv_sqrt_two = 1.0 / det::sqrt(2.0);
    let mut out = vec![C64::ZERO; count];
    let mut at = 0usize;
    for l in 0..=l_max {
        out[l * (l + 1)] = coefficients[at];
        at += 1;
        for m in 1..=l {
            let cosine = coefficients[at];
            let sine = coefficients[at + 1];
            at += 2;
            let i_sine = C64::new(-sine.im, sine.re);
            out[l * (l + 1) + m] = (cosine - i_sine).scale(inv_sqrt_two);
            out[l * (l + 1) - m] = (cosine + i_sine).scale(if m % 2 == 0 {
                inv_sqrt_two
            } else {
                -inv_sqrt_two
            });
        }
    }
    Ok(out)
}

/// Evaluate complex response coefficients on the documented real tesseral
/// basis at a nonzero body-frame direction.
pub fn evaluate_real_tesseral(
    l_max: usize,
    coefficients: &[C64],
    direction: [f64; 3],
) -> Result<C64, BroadbandRadiationError> {
    validate_coefficients(l_max, coefficients, "real tesseral")?;
    let norm_squared = direction.iter().map(|value| value * value).sum::<f64>();
    if !(norm_squared > 0.0 && norm_squared.is_finite()) {
        return Err(BroadbandRadiationError::InvalidInput(
            "evaluation direction must be nonzero and finite",
        ));
    }
    let x = direction[2] / det::sqrt(norm_squared);
    let phi = det::atan2(direction[1], direction[0]);
    let pbar = normalized_associated_legendre(l_max, x);
    let triangular = |l: usize, m: usize| l * (l + 1) / 2 + m;
    let mut result = C64::ZERO;
    let mut at = 0usize;
    for l in 0..=l_max {
        result = result + coefficients[at].scale(pbar[triangular(l, 0)]);
        at += 1;
        for m in 1..=l {
            let scale = det::sqrt(2.0) * pbar[triangular(l, m)];
            result = result
                + coefficients[at].scale(scale * det::cos(m as f64 * phi))
                + coefficients[at + 1].scale(scale * det::sin(m as f64 * phi));
            at += 2;
        }
    }
    Ok(result)
}

/// Build and withhold-gate a reusable body-frame broadband radiation bank from
/// solver-neutral samples.
///
/// Input samples use `exp(-i omega t)`, whereas fs-vfit uses `s=+i omega`;
/// every real-SH response is therefore conjugated before fitting. Physical
/// frequency `omega` is presented to the continuous fit at
/// `(2/T) tan(omega T/2)`, so an unprewarped Tustin realization samples the
/// intended physical response at the original digital frequency.
pub fn build_broadband_radiation_artifact(
    samples: &SampledRadiationData,
    controls: BroadbandRadiationControls,
) -> Result<BroadbandRadiationArtifact, BroadbandRadiationError> {
    validate_inputs(samples, controls)?;
    let sampling_diagnostics = summarize_diagnostics(samples);
    if sampling_diagnostics.captured_fraction < controls.minimum_captured_fraction {
        return Err(BroadbandRadiationError::CapturedFraction(
            sampling_diagnostics.captured_fraction,
            controls.minimum_captured_fraction,
        ));
    }
    let t_s = 1.0 / controls.sample_rate_hz;
    let channels = real_tesseral_channels(samples.l_max);
    let channel_count = channels.len();
    let mut training_samples =
        vec![vec![vec![C64::ZERO; samples.training.len()]; channel_count]; samples.input_ids.len()];

    for (frequency_index, sample) in samples.training.iter().enumerate() {
        for (input_index, coefficients) in sample.coefficients_by_input.iter().enumerate() {
            let real = complex_sh_to_real_tesseral(samples.l_max, coefficients)?;
            for channel in 0..channel_count {
                training_samples[input_index][channel][frequency_index] = real[channel];
            }
        }
    }

    let fit_omega: Vec<f64> = samples
        .training
        .iter()
        .map(|sample| 2.0 / t_s * det::tan(sample.omega_rad_s * t_s / 2.0))
        .collect();
    let fit_options = FitOptions {
        order: controls.fit_order,
        iterations: controls.fit_iterations,
        weights: controls.fit_weights,
        fit_e: false,
        fit_d: false,
    };
    let mut maximum_fit_rms = 0.0_f64;
    let mut artifact_inputs = Vec::with_capacity(samples.input_ids.len());
    for (input_index, id) in samples.input_ids.iter().enumerate() {
        let mut filters = Vec::with_capacity(channel_count);
        for channel in 0..channel_count {
            let response: Vec<C64> = training_samples[input_index][channel]
                .iter()
                .map(|sample| sample.conj())
                .collect();
            if response
                .iter()
                .all(|sample| sample.abs() <= controls.far_field_signal_floor)
            {
                filters.push(zero_filter(t_s));
                continue;
            }
            let fit = vector_fit(&fit_omega, &response, &fit_options)
                .map_err(|source| BroadbandRadiationError::Fit(input_index, channel, source))?;
            if !fit.model.is_stable() || fit.model.e != 0.0 || fit.model.d != 0.0 {
                return Err(BroadbandRadiationError::NumericalFailure(
                    "fit must be stable and strictly proper",
                ));
            }
            if !fit.report.weighted_rms.is_finite() {
                return Err(BroadbandRadiationError::NumericalFailure(
                    "fit diagnostic must be finite",
                ));
            }
            maximum_fit_rms = maximum_fit_rms.max(fit.report.weighted_rms);
            let discrete = bilinear_state_space(&fit.model, t_s, 0.0).map_err(|source| {
                BroadbandRadiationError::Discretize(input_index, channel, source)
            })?;
            if discrete.e_leftover != 0.0 {
                return Err(BroadbandRadiationError::NumericalFailure(
                    "proper fit produced a nonzero discrete improper remainder",
                ));
            }
            filters.push(discrete);
        }
        artifact_inputs.push(BroadbandRadiationInputArtifact {
            id: id.clone(),
            filters,
        });
    }

    let mut normalized_errors = Vec::new();
    for sample in &samples.held_out {
        for (input_index, expected_by_direction) in sample.far_field_by_input.iter().enumerate() {
            let mut predicted_real = Vec::with_capacity(channel_count);
            for filter in &artifact_inputs[input_index].filters {
                let fitted = if filter.n == 0 {
                    C64::ZERO
                } else {
                    filter.eval(sample.omega_rad_s).map_err(|_| {
                        BroadbandRadiationError::NumericalFailure(
                            "discrete realization evaluation failed",
                        )
                    })?
                };
                if !fitted.re.is_finite() || !fitted.im.is_finite() {
                    return Err(BroadbandRadiationError::NumericalFailure(
                        "discrete realization evaluation must be finite",
                    ));
                }
                predicted_real.push(fitted.conj());
            }
            for (&expected, &direction) in expected_by_direction.iter().zip(&sample.directions) {
                let actual = evaluate_real_tesseral(samples.l_max, &predicted_real, direction)?;
                normalized_errors.push(
                    (actual - expected).abs() / expected.abs().max(controls.far_field_signal_floor),
                );
            }
        }
    }

    let (maximum_error, rms_error) = error_summary(&normalized_errors);
    if maximum_error > controls.maximum_normalized_error
        || rms_error > controls.rms_normalized_error
    {
        return Err(BroadbandRadiationError::HeldOutError(
            maximum_error,
            rms_error,
            controls.maximum_normalized_error,
            controls.rms_normalized_error,
        ));
    }

    Ok(BroadbandRadiationArtifact {
        sample_interval_s: t_s,
        l_max: samples.l_max,
        channels,
        inputs: artifact_inputs,
        report: BroadbandRadiationReport {
            authority: BroadbandRadiationAuthority::EstimateOnly,
            source_id: samples.source_id.clone(),
            harmonic_time_convention: samples.harmonic_time_convention,
            training_omega_rad_s: samples
                .training
                .iter()
                .map(|sample| sample.omega_rad_s)
                .collect(),
            held_out_omega_rad_s: samples
                .held_out
                .iter()
                .map(|sample| sample.omega_rad_s)
                .collect(),
            maximum_normalized_complex_error: maximum_error,
            rms_normalized_complex_error: rms_error,
            controls,
            sampling_diagnostics,
            maximum_channel_weighted_rms: maximum_fit_rms,
            held_out_comparison_count: normalized_errors.len(),
            no_claims: BROADBAND_RADIATION_NO_CLAIMS,
        },
    })
}

fn validate_inputs(
    samples: &SampledRadiationData,
    controls: BroadbandRadiationControls,
) -> Result<(), BroadbandRadiationError> {
    let invalid = BroadbandRadiationError::InvalidInput;
    if samples.source_id.is_empty()
        || samples.input_ids.is_empty()
        || samples.input_ids.len() > MAX_BROADBAND_INPUTS
        || samples.input_ids.iter().any(String::is_empty)
        || samples
            .input_ids
            .iter()
            .enumerate()
            .any(|(index, id)| samples.input_ids[..index].contains(id))
        || samples.l_max > MAX_RADIATION_SH_DEGREE
    {
        return Err(invalid(
            "source/input IDs must be nonempty and unique, the input batch bounded, and l_max supported",
        ));
    }
    if !(controls.sample_rate_hz > 0.0 && controls.sample_rate_hz.is_finite())
        || controls.fit_order == 0
        || controls.fit_iterations == 0
        || !(0.0..=1.0).contains(&controls.minimum_captured_fraction)
        || controls.minimum_captured_fraction == 0.0
        || !(controls.far_field_signal_floor > 0.0 && controls.far_field_signal_floor.is_finite())
        || !(controls.maximum_normalized_error > 0.0
            && controls.maximum_normalized_error.is_finite())
        || !(controls.rms_normalized_error > 0.0 && controls.rms_normalized_error.is_finite())
    {
        return Err(invalid(
            "fit, sample-rate, capture, signal-floor, and error controls must be finite and within documented bounds",
        ));
    }
    let needed = controls
        .fit_order
        .checked_add(2)
        .ok_or_else(|| invalid("fit order is too large"))?;
    if samples.training.len() < needed
        || samples.training.len() > MAX_BROADBAND_FREQUENCIES
        || samples.held_out.is_empty()
        || samples.held_out.len() > MAX_BROADBAND_FREQUENCIES
    {
        return Err(invalid(
            "training needs order+2 samples and both grids must be nonempty and bounded",
        ));
    }
    let nyquist = core::f64::consts::PI * controls.sample_rate_hz;
    for (index, sample) in samples.training.iter().enumerate() {
        let omega = sample.omega_rad_s;
        let warped =
            2.0 * controls.sample_rate_hz * det::tan(omega / (2.0 * controls.sample_rate_hz));
        if !(omega > 0.0 && omega.is_finite() && omega < nyquist && warped.is_finite())
            || index > 0 && samples.training[index - 1].omega_rad_s >= omega
        {
            return Err(invalid(
                "frequency grids must be strictly increasing, positive, finite, Tustin-warpable, and below Nyquist",
            ));
        }
    }
    for (index, sample) in samples.held_out.iter().enumerate() {
        let omega = sample.omega_rad_s;
        let warped =
            2.0 * controls.sample_rate_hz * det::tan(omega / (2.0 * controls.sample_rate_hz));
        if !(omega > 0.0 && omega.is_finite() && omega < nyquist && warped.is_finite())
            || index > 0 && samples.held_out[index - 1].omega_rad_s >= omega
            || samples
                .training
                .iter()
                .any(|training| training.omega_rad_s == omega)
        {
            return Err(invalid(
                "held-out frequencies must be increasing, Tustin-warpable, below Nyquist, and disjoint from training",
            ));
        }
    }
    let coefficient_count = (samples.l_max + 1) * (samples.l_max + 1);
    let transfer_cells = samples
        .input_ids
        .len()
        .checked_mul(coefficient_count)
        .and_then(|count| count.checked_mul(samples.training.len()))
        .ok_or_else(|| invalid("transfer sample count overflowed"))?;
    if transfer_cells > MAX_BROADBAND_TRANSFER_SAMPLES {
        return Err(invalid("transfer sample budget exceeded"));
    }
    for sample in &samples.training {
        if sample.coefficients_by_input.len() != samples.input_ids.len()
            || sample.coefficients_by_input.iter().any(|row| {
                row.len() != coefficient_count
                    || row
                        .iter()
                        .any(|value| !value.re.is_finite() || !value.im.is_finite())
            })
            || !valid_diagnostic(sample.diagnostics)
        {
            return Err(invalid(
                "every training frequency needs finite canonical SH rows and valid diagnostics for every input",
            ));
        }
    }
    for sample in &samples.held_out {
        if !(8..=MAX_VALIDATION_DIRECTIONS).contains(&sample.directions.len())
            || sample.far_field_by_input.len() != samples.input_ids.len()
            || sample.directions.iter().any(|direction| {
                let norm_squared = direction.iter().map(|value| value * value).sum::<f64>();
                !(norm_squared > 0.0 && norm_squared.is_finite())
            })
            || sample.far_field_by_input.iter().any(|row| {
                row.len() != sample.directions.len()
                    || row
                        .iter()
                        .any(|value| !value.re.is_finite() || !value.im.is_finite())
            })
            || !valid_diagnostic(sample.diagnostics)
        {
            return Err(invalid(
                "each held-out frequency needs bounded finite directions, direct far fields, and valid diagnostics for every input",
            ));
        }
    }
    Ok(())
}

fn validate_coefficients(
    l_max: usize,
    coefficients: &[C64],
    family: &'static str,
) -> Result<(), BroadbandRadiationError> {
    if l_max > MAX_RADIATION_SH_DEGREE
        || coefficients.len() != (l_max + 1) * (l_max + 1)
        || coefficients
            .iter()
            .any(|value| !value.re.is_finite() || !value.im.is_finite())
    {
        return Err(BroadbandRadiationError::InvalidInput(match family {
            "complex SH" => "complex SH coefficients must be finite with count (l_max + 1)^2",
            _ => "real tesseral coefficients must be finite with count (l_max + 1)^2",
        }));
    }
    Ok(())
}

fn complex_sh_index(l: usize, m: i64) -> usize {
    let offset = usize::try_from(i64::try_from(l).expect("degree is bounded") + m)
        .expect("l + m is nonnegative");
    l * l + offset
}

fn normalized_associated_legendre(l_max: usize, x: f64) -> Vec<f64> {
    let sine = det::sqrt((1.0 - x * x).max(0.0));
    let mut values = vec![0.0; (l_max + 1) * (l_max + 2) / 2];
    let index = |l: usize, m: usize| l * (l + 1) / 2 + m;
    values[0] = 0.5 / det::sqrt(core::f64::consts::PI);
    for m in 1..=l_max {
        let mf = m as f64;
        values[index(m, m)] =
            -det::sqrt((2.0 * mf + 1.0) / (2.0 * mf)) * sine * values[index(m - 1, m - 1)];
    }
    for m in 0..l_max {
        let mf = m as f64;
        values[index(m + 1, m)] = det::sqrt(2.0 * mf + 3.0) * x * values[index(m, m)];
    }
    for m in 0..=l_max {
        for l in (m + 2)..=l_max {
            let lf = l as f64;
            let mf = m as f64;
            let a = det::sqrt((4.0 * lf * lf - 1.0) / (lf * lf - mf * mf));
            let previous = lf - 1.0;
            let b = det::sqrt((previous * previous - mf * mf) / (4.0 * previous * previous - 1.0));
            values[index(l, m)] = a * (x * values[index(l - 1, m)] - b * values[index(l - 2, m)]);
        }
    }
    values
}

fn valid_diagnostic(diagnostic: RadiationSampleDiagnostics) -> bool {
    diagnostic.captured_fraction > 0.0
        && diagnostic.captured_fraction <= 1.0 + 64.0 * f64::EPSILON
        && diagnostic.captured_fraction.is_finite()
        && diagnostic.panels_per_wavelength > 0.0
        && diagnostic.panels_per_wavelength.is_finite()
        && diagnostic.condition_lower_bound >= 1.0
        && diagnostic.condition_lower_bound.is_finite()
}

fn summarize_diagnostics(samples: &SampledRadiationData) -> RadiationSampleDiagnostics {
    let mut summary = RadiationSampleDiagnostics {
        captured_fraction: f64::INFINITY,
        panels_per_wavelength: f64::INFINITY,
        condition_lower_bound: 0.0,
    };
    for diagnostic in samples
        .training
        .iter()
        .map(|sample| sample.diagnostics)
        .chain(samples.held_out.iter().map(|sample| sample.diagnostics))
    {
        summary.captured_fraction = summary.captured_fraction.min(diagnostic.captured_fraction);
        summary.panels_per_wavelength = summary
            .panels_per_wavelength
            .min(diagnostic.panels_per_wavelength);
        summary.condition_lower_bound = summary
            .condition_lower_bound
            .max(diagnostic.condition_lower_bound);
    }
    summary
}

fn zero_filter(t_s: f64) -> DiscreteStateSpace {
    DiscreteStateSpace {
        n: 0,
        a: Vec::new(),
        b: Vec::new(),
        c: Vec::new(),
        d: 0.0,
        e_leftover: 0.0,
        t_s,
    }
}

fn error_summary(errors: &[f64]) -> (f64, f64) {
    let maximum = errors.iter().copied().fold(0.0_f64, f64::max);
    let mean_square = errors.iter().map(|error| error * error).sum::<f64>() / errors.len() as f64;
    (maximum, det::sqrt(mean_square))
}

#[cfg(test)]
mod tests {
    use super::*;
    use fs_bem::helmholtz::DirectivityTable;

    #[test]
    fn g0_real_tesseral_change_reconstructs_complex_basis_and_signs() {
        let l_max: usize = 3;
        let coefficients: Vec<C64> = (0..(l_max + 1).pow(2))
            .map(|index| C64::new(index as f64 * 0.13 - 0.7, 0.4 - index as f64 * 0.09))
            .collect();
        let complex = DirectivityTable {
            k: 2.0,
            l_max,
            coefficients: coefficients.clone(),
            captured_fraction: 1.0,
        };
        let real = complex_sh_to_real_tesseral(l_max, &complex.coefficients).unwrap();
        let round_trip = real_tesseral_to_complex_sh(l_max, &real).unwrap();
        for (actual, expected) in round_trip.iter().zip(&coefficients) {
            assert!((*actual - *expected).abs() < 2.0e-15);
        }
        for direction in [[0.3, -0.4, 0.8], [-0.9, 0.2, 0.1], [0.1, 0.7, -0.6]] {
            let actual = evaluate_real_tesseral(l_max, &real, direction).unwrap();
            assert!((actual - complex.evaluate(direction)).abs() < 2.0e-14);
        }
    }

    fn controls() -> BroadbandRadiationControls {
        BroadbandRadiationControls {
            sample_rate_hz: 48_000.0,
            minimum_captured_fraction: 0.75,
            fit_order: 1,
            fit_iterations: 8,
            fit_weights: WeightPreset::Uniform,
            far_field_signal_floor: 1.0e-12,
            maximum_normalized_error: 1.0e-7,
            rms_normalized_error: 1.0e-7,
        }
    }

    fn diagnostic() -> RadiationSampleDiagnostics {
        RadiationSampleDiagnostics {
            captured_fraction: 1.0,
            panels_per_wavelength: 12.0,
            condition_lower_bound: 2.0,
        }
    }

    fn warped_one_pole_bem_coefficient(omega_rad_s: f64) -> C64 {
        let sample_rate_hz = controls().sample_rate_hz;
        let warped = 2.0 * sample_rate_hz * det::tan(omega_rad_s / (2.0 * sample_rate_hz));
        // H_-(omega) = conj(H_+(i omega_warp)), H_+(s)=1/(s+1200).
        C64::new(1_200.0, -warped).recip()
    }

    fn neutral_one_pole_samples() -> SampledRadiationData {
        let directions = vec![
            [1.0, 0.0, 0.0],
            [-1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, -1.0, 0.0],
            [0.0, 0.0, 1.0],
            [0.0, 0.0, -1.0],
            [1.0, 1.0, 1.0],
            [-1.0, 1.0, -1.0],
        ];
        let y00 = 0.5 / det::sqrt(core::f64::consts::PI);
        let training = [400.0, 700.0, 1_100.0, 1_700.0, 2_600.0, 3_900.0]
            .into_iter()
            .map(|omega_rad_s| ComplexShTrainingSample {
                omega_rad_s,
                coefficients_by_input: vec![vec![warped_one_pole_bem_coefficient(omega_rad_s)]],
                diagnostics: diagnostic(),
            })
            .collect();
        let held_out = [550.0, 900.0, 1_400.0, 2_100.0, 3_200.0]
            .into_iter()
            .map(|omega_rad_s| {
                let direct = warped_one_pole_bem_coefficient(omega_rad_s).scale(y00);
                DirectFarFieldHeldOutSample {
                    omega_rad_s,
                    directions: directions.clone(),
                    far_field_by_input: vec![vec![direct; directions.len()]],
                    diagnostics: diagnostic(),
                }
            })
            .collect();
        SampledRadiationData {
            source_id: "analytic-one-pole-oracle".to_owned(),
            harmonic_time_convention: HarmonicTimeConvention::ExpNegativeIOmegaT,
            l_max: 0,
            input_ids: vec!["uniform-pulsation".to_owned()],
            training,
            held_out,
        }
    }

    #[test]
    fn g0_neutral_train_withhold_is_deterministic_and_strictly_proper() {
        let samples = neutral_one_pole_samples();
        let first = build_broadband_radiation_artifact(&samples, controls()).unwrap();
        let second = build_broadband_radiation_artifact(&samples, controls()).unwrap();
        assert_eq!(first, second);
        assert!(first.channels.len() == 1 && first.inputs[0].filters.len() == 1);
        assert_eq!(
            first.report.authority,
            BroadbandRadiationAuthority::EstimateOnly
        );
        assert_eq!(first.report.held_out_comparison_count, 40);
        assert!(first.report.maximum_normalized_complex_error < 1.0e-7);
        assert!(
            first.inputs[0]
                .filters
                .iter()
                .all(|filter| filter.e_leftover == 0.0)
        );
    }

    #[test]
    fn g0_admission_and_complete_withheld_gate_catch_mutations() {
        let mut overlap = neutral_one_pole_samples();
        overlap.held_out[0].omega_rad_s = overlap.training[0].omega_rad_s;
        assert!(validate_inputs(&overlap, controls()).is_err());

        let mut nyquist = neutral_one_pole_samples();
        nyquist.training[5].omega_rad_s = core::f64::consts::PI * 48_000.0;
        assert!(validate_inputs(&nyquist, controls()).is_err());

        let mut phase_mutation = neutral_one_pole_samples();
        phase_mutation.held_out[2].far_field_by_input[0][4] =
            -phase_mutation.held_out[2].far_field_by_input[0][4];
        assert!(matches!(
            build_broadband_radiation_artifact(&phase_mutation, controls()),
            Err(BroadbandRadiationError::HeldOutError(..))
        ));
    }
}
