//! Cold construction of a mutual load from the original baffled-piston owner.
//!
//! This is a finite-band proper rational approximation of stationary, uniform
//! circular openings in one infinite rigid baffle and homogeneous inviscid gas.
//! The matrix retains every self and mutual term. Propagation is approximated
//! by the shared poles; no separate delay is extracted or silently discarded.
//! Three requested orders are inspected, with withheld physical samples and
//! both continuous and actual interval-mean ZOH passivity admission. These
//! sampled errors do not certify unsampled frequencies, closed-loop accuracy,
//! discrete stored energy, or the complete M2.5-W3b instrument qualification.
//! A fit that is not continuously passive gets the M2.5-W3b residue repair
//! (poles, D and E fixed). Its in-band change must stay within a tenth of the
//! error allowance, and the repaired fit passes the same error and passivity
//! gates. An order that passes unchanged is still preferred to any repaired
//! order. No clipping, added damping or global gain is applied.

/// Mixed rectangular and circular flow sources through the same fitting owner.
pub mod mixed;
/// Finite material-body diffraction through the same rational realization owner.
pub mod boundary;

use super::{MimoRadiationLimits, MimoRadiationModel};
use crate::pcm_wav::baffled::openings::CircularOpening;
use fs_bem::helmholtz::{HelmholtzError, Medium};
use fs_bem::helmholtz::aperture::{
    BaffledPiston, baffled_piston_array_impedance, reciprocity_and_passivity,
};
use fs_blake3::DomainHasher;
use fs_exec::CancelGate;
use fs_material::gas::GasState;
use fs_math::{c64::C64, det};
use fs_vfit::{FitOptions, PoleTerm, WeightPreset};
use fs_vfit::causality::{CausalityError, CausalityOptions, check_mimo_causality};
use fs_vfit::mimo::{
    MimoModel, MimoPassivityReport, check_mimo_passivity, mimo_vector_fit,
    repair_mimo_passivity_with_budget, shared_pole_relocation_scratch_bytes,
};

/// A complete bounded DG side-port set, its axial outlet and an independently
/// driven resolved mouth. Every opening retains all mutual terms; the shared
/// fitter eliminates local numerators before solving its common denominator
/// for larger arrays. Existing fit-memory and runtime-state limits still apply.
pub const MAX_BAFFLED_MIMO_PORTS: usize = fs_duct::dg::MAX_DG_SIDE_PORTS + 2;
/// Maximum rows in each fitting, causality or withheld frequency sweep.
pub const MAX_BAFFLED_MIMO_FREQUENCIES: usize = 513;
/// Maximum admitted allowance for the dense fitting and certification scratch.
pub const MAX_BAFFLED_MIMO_FIT_MEMORY_BYTES: usize = 512 * 1024 * 1024;

/// Explicit numerical choices; no instrument name chooses a coefficient.
#[derive(Clone, Copy, Debug)]
pub struct BaffledMimoBakeOptions {
    /// Log-spaced training rows, including both declared band endpoints.
    /// Withheld geometric midpoints use `frequency_samples - 1` distinct rows.
    /// An additional uniform half-step sweep checks source-data causality.
    pub frequency_samples: usize,
    /// Exactly three strictly increasing shared-pole orders, each in 1..=64.
    pub orders: [usize; 3],
    /// Shared VF relocation cap per order, in 1..=32.
    pub fit_iterations: usize,
    /// Maximum matrix-relative Frobenius error on both training and held-out rows.
    pub maximum_relative_error: f64,
    /// Conservative scratch admission, checked before sampling or fitting.
    pub max_fit_memory_bytes: usize,
    /// Explicit assumptions for the original dispersion-relation consistency test.
    pub causality: CausalityOptions,
}

/// One original fit attempt; refusals never become an accepted lower-quality load.
#[derive(Clone, Debug)]
pub struct BaffledMimoBakeRung {
    /// Requested shared-pole order.
    pub order: usize,
    /// Recomputed error against the original training matrix, when finite.
    pub fit_relative_error: Option<f64>,
    /// Error at independently sampled interleaved frequencies, when finite.
    pub held_out_relative_error: Option<f64>,
    /// True only if the continuous owner admitted this unchanged fit.
    pub continuous_passive: bool,
    /// Largest in-band matrix-relative change of an applied residue repair.
    /// The errors above and the sampled result are then the repaired fit's.
    pub repair_relative_change: Option<f64>,
    /// True only if the actual interval-mean sampled owner admitted this fit.
    pub sampled_passive: bool,
    /// Original fit, error, discretization or passivity refusal, if any.
    pub refusal: Option<String>,
}

/// Numerical diagnostics of the selected physical-source fit (repaired only
/// when its rung records a residue repair).
#[derive(Clone, Debug)]
pub struct BaffledMimoBakeReport {
    /// Smallest requested order passing every gate unchanged, after inspecting
    /// all three; only if none does, the smallest passing after residue repair.
    pub selected_order: usize,
    /// Number of original physical training samples.
    pub frequency_samples: usize,
    /// Number of original physical samples withheld from identification.
    pub held_out_samples: usize,
    /// Selected model's largest training matrix-relative error.
    pub maximum_fit_relative_error: f64,
    /// Selected model's largest withheld matrix-relative error.
    pub maximum_held_out_relative_error: f64,
    /// Worst dispersion-relation miss divided by its declared error bound.
    pub maximum_causality_ratio: f64,
    /// Selected model's continuous passivity result; its certificate class is retained.
    pub continuous: MimoPassivityReport,
    /// Actual mean-output ZOH grid result; this remains `GridOnly`.
    pub sampled: MimoPassivityReport,
    /// All three requested attempts, including refusals and nonmonotone errors.
    pub rungs: Vec<BaffledMimoBakeRung>,
}

/// Cold-build failure; no partially built model or runtime is published.
#[derive(Debug)]
pub enum BaffledMimoBakeError {
    /// Invalid geometry, gas, clock, numerical option or resource declaration.
    Invalid(&'static str),
    /// Original shared surface geometry or quadrature admission refused.
    Surface(String),
    /// A declared working or retained memory allowance is insufficient.
    Memory { /// Which separate allowance failed.
        kind: &'static str, /// Conservative required bytes.
        required: usize, /// Caller supplied bytes.
        allowed: usize },
    /// Original physical quadrature refusal at a specific frequency.
    Sample { /// Physical frequency [Hz].
        frequency_hz: f64, /// Original BEM/piston owner's error.
        source: HelmholtzError },
    /// Original source-data causality consistency test refused.
    Causality(CausalityError),
    /// None of the bounded three orders passed, with or without residue repair.
    NoConvergedFit { /// Original diagnostics in requested order.
        rungs: Vec<BaffledMimoBakeRung> },
    /// Cancellation before publishing an admitted result.
    Cancelled,
}
impl core::fmt::Display for BaffledMimoBakeError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Invalid(what) => write!(f, "baffled MIMO bake: {what}"),
            Self::Surface(what) => write!(f, "baffled MIMO surface: {what}"),
            Self::Memory { kind, required, allowed } => write!(f,
                "baffled MIMO {kind} needs {required} bytes, allowance is {allowed}"),
            Self::Sample { frequency_hz, source } => write!(f,
                "baffled MIMO sample at {frequency_hz} Hz: {source}"),
            Self::Causality(source) => write!(f, "baffled MIMO source: {source}"),
            Self::NoConvergedFit { rungs } => write!(f,
                "baffled MIMO order ladder refused: {rungs:?}"),
            Self::Cancelled => f.write_str("baffled MIMO bake cancelled"),
        }
    }
}
impl std::error::Error for BaffledMimoBakeError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Sample { source, .. } => Some(source),
            Self::Causality(source) => Some(source),
            _ => None,
        }
    }
}
fn invalid(what: &'static str) -> BaffledMimoBakeError { BaffledMimoBakeError::Invalid(what) }
fn poll(gate: &CancelGate) -> Result<(), BaffledMimoBakeError> {
    if gate.is_requested() { Err(BaffledMimoBakeError::Cancelled) } else { Ok(()) }
}

/// Derive one matching load from actual ordered aperture geometry and gas data.
///
/// The band must include the complete mechanics Nyquist limit, irrespective of
/// microphone output rate. Source samples use the existing adaptive disc rule;
/// force/velocity is divided by BOTH aperture areas and conjugated once from
/// `exp(-i omega t)` to the fitter's `s=+i omega` convention. Gas substitutions
/// rerun the same owners, rather than rescaling a previous realization.
///
/// Cancellation is checked between physical rows, causality, and fit orders.
/// Existing individual quadrature, VF and eigensolves have no mid-call poll.
///
/// # Errors
/// Invalid source/budgets, unresolved physical quadrature, source causality
/// refusal, exhausted accuracy/passivity ladder, or cancellation.
#[allow(clippy::too_many_arguments)]
pub fn bake_baffled_mimo(
    exterior: GasState,
    apertures: &[CircularOpening],
    mechanics_rate_hz: u32,
    band_hz: [f64; 2],
    options: BaffledMimoBakeOptions,
    limits: MimoRadiationLimits,
    gate: &CancelGate,
) -> Result<(MimoRadiationModel, BaffledMimoBakeReport), BaffledMimoBakeError> {
    poll(gate)?;
    admit(&exterior, apertures, mechanics_rate_hz, band_hz, options, limits)?;
    let n = apertures.len();
    let pistons: Vec<_> = apertures.iter().map(|a| BaffledPiston {
        center_m: a.center_m, radius_m: a.radius_m,
    }).collect();
    let (impedance, report) = fit_source(
        n, mechanics_rate_hz, band_hz, options, limits, gate,
        |hz, gate| sample_rows(&pistons, &exterior, hz, gate),
    )?;
    let declaration_identity = declaration(&exterior, apertures, mechanics_rate_hz, band_hz);
    let realization_identity = realization(declaration_identity, &impedance, options);
    Ok((MimoRadiationModel { declaration_identity, realization_identity, exterior,
        apertures: apertures.to_vec(), band_hz, impedance }, report))
}

// One unchanged fitting/certification owner for every physical source sampler.
// Geometry/material/quadrature admission belongs to each caller before this sweep.
#[allow(clippy::too_many_arguments)]
fn fit_source<E: From<BaffledMimoBakeError>>(
    n: usize, mechanics_rate_hz: u32, band_hz: [f64; 2],
    options: BaffledMimoBakeOptions, limits: MimoRadiationLimits, gate: &CancelGate,
    sample: impl Fn(&[f64], &CancelGate) -> Result<Vec<Vec<C64>>, E>,
) -> Result<(MimoModel, BaffledMimoBakeReport), E> {
    let count = options.frequency_samples;
    let span = det::ln(band_hz[1] / band_hz[0]);
    let mut train_hz: Vec<_> = (0..count).map(|i|
        band_hz[0] * det::exp(span * i as f64 / (count - 1) as f64)).collect();
    train_hz[0] = band_hz[0];
    train_hz[count - 1] = band_hz[1];
    let held_hz: Vec<_> = train_hz.windows(2).map(|f| f[0] * det::sqrt(f[1] / f[0])).collect();
    let training = sample(&train_hz, gate)?;
    let held_out = sample(&held_hz, gate)?;
    // The causality owner's half-step grid never asks the strictly-positive
    // piston owner for DC. The highest ORIGINAL sample is exactly band_hz[1].
    let causal_hz: Vec<_> = (0..count).map(|i|
        band_hz[1] * (i as f64 + 0.5) / (count as f64 - 0.5)).collect();
    let causal_data = sample(&causal_hz, gate)?;
    let omega = |frequencies: &[f64]| -> Vec<f64> {
        frequencies.iter().map(|f| core::f64::consts::TAU * f).collect()
    };
    let train_omega = omega(&train_hz);
    let held_omega = omega(&held_hz);
    poll(gate)?;
    let causality = check_mimo_causality(&omega(&causal_hz), &causal_data, n,
        &vec![0.0; n * n], &options.causality).map_err(BaffledMimoBakeError::Causality)?;
    let maximum_causality_ratio = causality.iter().map(|r| r.worst_ratio).fold(0.0, f64::max);
    let mut rungs = Vec::with_capacity(3);
    let (mut selected, mut repaired) = (None, None);
    let band = (core::f64::consts::TAU * band_hz[0], core::f64::consts::TAU * band_hz[1]);
    for order in options.orders {
        poll(gate)?;
        let mut rung = BaffledMimoBakeRung { order,
            fit_relative_error: None, held_out_relative_error: None, continuous_passive: false,
            repair_relative_change: None, sampled_passive: false, refusal: None };
        let result = (|| -> Result<_, String> {
            let fit = FitOptions { order, iterations: options.fit_iterations,
                weights: WeightPreset::InverseMagnitude, fit_d: true, fit_e: false };
            let (model, _) = mimo_vector_fit(&train_omega, &training, n, &fit, true)
                .map_err(|e| e.to_string())?;
            if !finite_model(&model) || model.order() != order {
                return Err("shared fit must retain its requested stable finite proper order".into());
            }
            let training_error = relative_error(&model, &train_omega, &training)?;
            let held_error = relative_error(&model, &held_omega, &held_out)?;
            rung.fit_relative_error = Some(training_error);
            rung.held_out_relative_error = Some(held_error);
            if training_error > options.maximum_relative_error || held_error > options.maximum_relative_error {
                return Err("unchanged fit exceeds training or withheld matrix error allowance".into());
            }
            let mut continuous = check_mimo_passivity(&model, band).map_err(|e| e.to_string())?;
            rung.continuous_passive = continuous.passive;
            let (model, training_error, held_error) = if continuous.passive {
                (model, training_error, held_error)
            } else {
                // M2.5-W3b residue repair, counted in the operator-error budget.
                let available = options.max_fit_memory_bytes
                    .saturating_sub(source_and_passivity_scratch_bytes(n, options));
                let (repaired, repair) = repair_mimo_passivity_with_budget(
                    &model, band, available).map_err(|e|
                    format!("unchanged continuous fit is not passive and its residue repair refused: {e}"))?;
                rung.repair_relative_change = Some(repair.max_relative_change);
                if !finite_model(&repaired)
                    || !(repair.max_relative_change <= 0.1 * options.maximum_relative_error)
                {
                    return Err("residue repair exceeds a tenth of the matrix error allowance".into());
                }
                let training_error = relative_error(&repaired, &train_omega, &training)?;
                let held_error = relative_error(&repaired, &held_omega, &held_out)?;
                rung.fit_relative_error = Some(training_error);
                rung.held_out_relative_error = Some(held_error);
                if training_error > options.maximum_relative_error || held_error > options.maximum_relative_error {
                    return Err("repaired fit exceeds training or withheld matrix error allowance".into());
                }
                continuous = check_mimo_passivity(&repaired, band).map_err(|e| e.to_string())?;
                if !continuous.passive { return Err("repaired continuous fit is not passive".into()); }
                (repaired, training_error, held_error)
            };
            let sampled = model.zoh(1.0 / f64::from(mechanics_rate_hz))
                .map_err(|e| e.to_string())?.mean_sampled_passivity(limits.sampled_passivity_points);
            rung.sampled_passive = sampled.passive;
            if !sampled.passive { return Err("actual interval-mean ZOH fit is not sampled-passive".into()); }
            Ok((model, continuous, sampled, training_error, held_error))
        })();
        poll(gate)?;
        match result {
            Ok(candidate) if rung.repair_relative_change.is_none() => { selected.get_or_insert(candidate); }
            Ok(candidate) => { repaired.get_or_insert(candidate); }
            Err(reason) => rung.refusal = Some(reason),
        }
        rungs.push(rung);
    }
    // An order that passes unchanged is preferred to any repaired order.
    let Some((impedance, continuous, sampled, maximum_fit_relative_error,
        maximum_held_out_relative_error)) = selected.or(repaired) else {
        return Err(BaffledMimoBakeError::NoConvergedFit { rungs }.into());
    };
    let report = BaffledMimoBakeReport { selected_order: impedance.order(),
        frequency_samples: count, held_out_samples: held_hz.len(),
        maximum_fit_relative_error, maximum_held_out_relative_error, maximum_causality_ratio,
        continuous, sampled, rungs };
    poll(gate)?;
    Ok((impedance, report))
}

fn admit(gas: &GasState, apertures: &[CircularOpening], rate: u32, band: [f64; 2],
    opts: BaffledMimoBakeOptions, limits: MimoRadiationLimits) -> Result<(), BaffledMimoBakeError>
{
    admit_common(gas, apertures.len(), rate, band, opts, limits)?;
    admit_apertures(apertures)
}

fn admit_apertures(apertures: &[CircularOpening]) -> Result<(), BaffledMimoBakeError> {
    for (i, a) in apertures.iter().enumerate() {
        let area = core::f64::consts::PI * a.radius_m * a.radius_m;
        if !a.radius_m.is_finite() || a.radius_m <= 0.0
            || !area.is_finite() || area <= 0.0 || a.center_m.iter().any(|x| !x.is_finite())
            || apertures[..i].iter().any(|b| {
                let distance = det::hypot(a.center_m[0] - b.center_m[0], a.center_m[1] - b.center_m[1]);
                !distance.is_finite() || distance <= a.radius_m + b.radius_m
            })
        {
            return Err(invalid("circular apertures must be finite, positive and separated"));
        }
    }
    Ok(())
}

fn admit_common(gas: &GasState, n: usize, rate: u32, band: [f64; 2],
    opts: BaffledMimoBakeOptions, limits: MimoRadiationLimits) -> Result<(), BaffledMimoBakeError>
{
    gas.validate_coherence().map_err(|_| invalid("source gas must be coherent"))?;
    let order = opts.orders[2];
    if !(1..=MAX_BAFFLED_MIMO_PORTS).contains(&n) || rate == 0
        || !band.iter().all(|f| f.is_finite() && *f > 0.0)
        || band[0] >= band[1] || !(band[1] / band[0]).is_finite()
        || !(core::f64::consts::TAU * band[1]).is_finite()
        || band[1] < 0.5 * f64::from(rate)
        || !(32..=MAX_BAFFLED_MIMO_FREQUENCIES).contains(&opts.frequency_samples)
        || opts.orders.iter().any(|o| !(1..=64).contains(o))
        || !opts.orders.windows(2).all(|o| o[0] < o[1])
        || opts.frequency_samples < order + 2
        || !(1..=32).contains(&opts.fit_iterations)
        || !opts.maximum_relative_error.is_finite()
        || !(0.0 < opts.maximum_relative_error && opts.maximum_relative_error <= 1.0)
        || opts.max_fit_memory_bytes == 0 || opts.max_fit_memory_bytes > MAX_BAFFLED_MIMO_FIT_MEMORY_BYTES
        || !(1..=8).contains(&opts.causality.subtractions)
        || opts.frequency_samples < 8 * (opts.causality.subtractions + 1)
        || !opts.causality.tail_bound_ratio.is_finite() || opts.causality.tail_bound_ratio < 1.0
        || !opts.causality.relative_data_tolerance.is_finite()
        || opts.causality.relative_data_tolerance < 0.0
        || opts.causality.relative_data_tolerance > opts.maximum_relative_error
        || limits.max_states > 256 || order * n > limits.max_states
        || !(17..=65_537).contains(&limits.sampled_passivity_points)
    {
        return Err(invalid("source, full-clock band, three-order ladder or explicit resource limits are inadmissible"));
    }
    // The shared VF keeps every upper-triangle response. Account for its
    // channel elimination, equilibrated QR, vectors and dense eigensolve scratch
    // BEFORE either the source sweep or the fitter can allocate it. Hard caps
    // above keep every product representable on wasm32.
    let states = order * n;
    let required = fit_scratch_bytes(n, opts);
    if required > opts.max_fit_memory_bytes {
        return Err(BaffledMimoBakeError::Memory { kind: "fit scratch", required,
            allowed: opts.max_fit_memory_bytes });
    }
    // Conservative term count = order; the actual host later checks its exact
    // real/complex decomposition with the same retained-payload accounting.
    let required = n * n * (64 * order + 64) + states * 64 + n * 256
        + (limits.sampled_passivity_points + 6400 + order * 128) * 8;
    if required > limits.max_memory_bytes {
        return Err(BaffledMimoBakeError::Memory { kind: "runtime payload", required,
            allowed: limits.max_memory_bytes });
    }
    Ok(())
}

fn fit_scratch_bytes(n: usize, opts: BaffledMimoBakeOptions) -> usize {
    let order = opts.orders[2];
    let responses = n * (n + 1) / 2;
    let fit = FitOptions { order, iterations: opts.fit_iterations,
        weights: WeightPreset::InverseMagnitude, fit_d: true, fit_e: false };
    shared_pole_relocation_scratch_bytes(opts.frequency_samples, responses, &fit)
        .and_then(|relocation| relocation.checked_add(source_and_passivity_scratch_bytes(n, opts)))
        .unwrap_or(usize::MAX)
}

fn source_and_passivity_scratch_bytes(n: usize, opts: BaffledMimoBakeOptions) -> usize {
    let states = opts.orders[2] * n;
    8 * (32 * states * states + 16 * opts.frequency_samples * n * n) + 1024 * 1024
}

fn sample_rows(pistons: &[BaffledPiston], gas: &GasState, hz: &[f64], gate: &CancelGate)
    -> Result<Vec<Vec<C64>>, BaffledMimoBakeError>
{
    let n = pistons.len();
    let areas: Vec<_> = pistons.iter().map(|p| core::f64::consts::PI * p.radius_m * p.radius_m).collect();
    let mut rows = Vec::with_capacity(hz.len());
    for &frequency_hz in hz {
        poll(gate)?;
        let force_velocity = baffled_piston_array_impedance(pistons,
            core::f64::consts::TAU * frequency_hz / gas.sound_speed, Medium::from(gas))
            .map_err(|source| BaffledMimoBakeError::Sample { frequency_hz, source })?;
        let row: Vec<_> = force_velocity.iter().enumerate().map(|(i, z)|
            z.conj().scale(1.0 / areas[i / n] / areas[i % n])).collect();
        if row.iter().any(|z| !z.re.is_finite() || !z.im.is_finite()) {
            return Err(invalid("physical acoustic matrix conversion is not finite"));
        }
        let (reciprocity, margin) = reciprocity_and_passivity(&row, n);
        if !reciprocity.is_finite() || !margin.is_finite() || reciprocity > 1e-12 || margin < -1e-12 {
            return Err(invalid("original physical matrix fails reciprocity or sampled passivity"));
        }
        rows.push(row);
    }
    Ok(rows)
}

fn finite_model(m: &MimoModel) -> bool {
    m.is_stable() && m.d.iter().all(|v| v.is_finite()) && m.e.iter().all(|v| *v == 0.0)
        && m.residues.iter().flatten().all(|z| z.re.is_finite() && z.im.is_finite())
        && m.poles.iter().all(|pole| match pole {
            PoleTerm::Real { pole, .. } => pole.is_finite(),
            PoleTerm::Pair { pole, .. } => pole.re.is_finite() && pole.im.is_finite() && pole.im > 0.0,
        })
}
fn relative_error(model: &MimoModel, omega: &[f64], data: &[Vec<C64>]) -> Result<f64, String> {
    let mut maximum = 0.0_f64;
    for (&w, row) in omega.iter().zip(data) {
        let norm = det::sqrt(row.iter().map(|z| z.norm_sq()).sum());
        let gap = det::sqrt(model.eval_iw(w).iter().zip(row).map(|(a, b)| (*a - *b).norm_sq()).sum());
        let error = gap / norm;
        if !error.is_finite() { return Err("fit error is not finite against the physical source matrix".into()); }
        maximum = maximum.max(error);
    }
    Ok(maximum)
}

fn scalar(hash: &mut DomainHasher, value: f64) { hash.update(&value.to_bits().to_le_bytes()); }
fn declaration(gas: &GasState, apertures: &[CircularOpening], rate: u32, band: [f64; 2]) -> [u8; 32] {
    let mut hash = DomainHasher::new("org.frankensim.fs-couple.baffled-mimo-source.v1");
    for value in [gas.temperature, gas.pressure, gas.density, gas.sound_speed, gas.dynamic_viscosity,
        gas.thermal_conductivity, gas.gamma, gas.specific_gas_constant, gas.specific_heat_cp,
        gas.prandtl, gas.characteristic_impedance, gas.water_mole_fraction] { scalar(&mut hash, value); }
    hash.update(&rate.to_le_bytes());
    for f in band { scalar(&mut hash, f); }
    hash.update(&(apertures.len() as u64).to_le_bytes());
    for a in apertures { for v in [a.center_m[0], a.center_m[1], a.radius_m] { scalar(&mut hash, v); } }
    hash.finalize().0
}
fn realization(source: [u8; 32], model: &MimoModel, opts: BaffledMimoBakeOptions) -> [u8; 32] {
    let mut hash = DomainHasher::new("org.frankensim.fs-couple.baffled-mimo-realization.v1");
    hash.update(&source);
    for value in [opts.frequency_samples, opts.orders[0], opts.orders[1], opts.orders[2],
        opts.fit_iterations, opts.max_fit_memory_bytes, opts.causality.subtractions, model.poles.len()] {
        hash.update(&(value as u64).to_le_bytes());
    }
    for value in [opts.maximum_relative_error, opts.causality.tail_bound_ratio,
        opts.causality.relative_data_tolerance] { scalar(&mut hash, value); }
    for pole in &model.poles {
        match pole {
            PoleTerm::Real { pole, .. } => { hash.update(&[0]); scalar(&mut hash, *pole); },
            PoleTerm::Pair { pole, .. } => {
                hash.update(&[1]); scalar(&mut hash, pole.re); scalar(&mut hash, pole.im);
            }
        }
    }
    for z in model.residues.iter().flatten() { scalar(&mut hash, z.re); scalar(&mut hash, z.im); }
    for value in model.d.iter().chain(&model.e) { scalar(&mut hash, *value); }
    hash.finalize().0
}

#[cfg(test)]
mod tests;
