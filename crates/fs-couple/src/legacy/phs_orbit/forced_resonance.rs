//! Local finite-amplitude resonance of the original complete forced pHS.
//!
//! The shared pHS survey owns peak discovery/refinement. Each physical query
//! corrects a complete-state predictor with the original harmonic-balance solve,
//! then observes the original nonlinear port law and its analytic frequency
//! derivative. A fixed anchor makes query order irrelevant to seed selection.

use super::{
    PeriodicPortInput, PhsOrbitError,
    sideband::{PhsPortSidebandControls, PhsPortSidebandError, PhsPortSidebandSource},
    solve_forced_harmonic_balance_scaled_seeded,
    spectrum::{PortSpectrumControls, PortSpectrumError, observe_port_spectrum},
};
use fs_math::{c64::C64, det};
use fs_orbit::{HbAnchor, HbBudget, HbOrbit, HbScales, sideband::SidebandScales};
use fs_phs::{
    PortHamiltonian,
    harmonic::resonance::{
        ResonanceSearchError, ResonanceSurvey, ResonanceSurveyControls, ResonanceSurveyError,
        survey_resonances,
    },
};

/// Existing physical/numerical owners plus an explicit local-branch admission.
#[derive(Clone, Debug)]
pub struct PhsForcedResonanceControls {
    /// Fixed retained state order and original bounded Newton/AFT solve.
    /// The order must equal `tangent.sideband.harmonics` and retain the input.
    pub harmonic_balance: HbBudget,
    pub harmonic_balance_scales: HbScales,
    pub tangent: PhsPortSidebandControls,
    pub tangent_scales: SidebandScales,
    /// Original physical output spectrum. Its output order must equal the
    /// tangent's output order, independently of the retained state order.
    pub spectrum: PortSpectrumControls,
    /// Positive maximum |ln(omega / anchor.omega)|, at most one. The complete
    /// grid must fit; no continuation steps or extrapolated brackets are added.
    pub maximum_log_frequency_distance: f64,
    /// Positive absolute RMS distance from the analytic anchor predictor,
    /// per ORIGINAL state coordinate, including every material-memory state.
    /// DC has weight one and each conjugate harmonic pair has weight two.
    /// This is a numerical branch-distance gate, not a branch-uniqueness proof.
    pub predictor_rms_tolerances: Vec<f64>,
    /// One total query allowance for the grid and all shared peak refinements.
    pub survey: ResonanceSurveyControls,
}

/// Original physical peaks and the largest admitted predictor corrections.
#[derive(Clone, Debug)]
pub struct PhsForcedResonanceSurvey {
    /// Frequencies and derivatives refer to the DRIVE FUNDAMENTAL [rad/s].
    /// The selected physical output line is `output_harmonic * omega`.
    /// Observations are peak phasors under exp(-i*k*omega*t), with the authored
    /// input phase retained. The reported backward residual is the original
    /// differentiated-equation solve's residual, not a forward-error estimate.
    pub resonances: ResonanceSurvey,
    /// Per-state physical RMS correction maxima over all original HB queries.
    pub maximum_predictor_correction_rms: Vec<f64>,
}

/// Native failures remain distinguishable from leaving the local predictor.
#[derive(Clone, Debug, PartialEq)]
pub enum PhsForcedResonanceError {
    Invalid(&'static str),
    Orbit(PhsOrbitError),
    Tangent(PhsPortSidebandError),
    Spectrum(PortSpectrumError),
    PredictorDisagreement {
        coordinate: usize,
        correction_rms: f64,
        tolerance: f64,
    },
    Cancelled,
}
impl core::fmt::Display for PhsForcedResonanceError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "forced pHS resonance query refused: {self:?}")
    }
}
impl std::error::Error for PhsForcedResonanceError {}

type Result<T> = core::result::Result<T, PhsForcedResonanceError>;

/// Survey finite-amplitude peaks of one original physical output harmonic.
///
/// The original Hamiltonian, J/R/G, complete state (including explicit memory),
/// and authored input DC/complex phasors stay fixed. The input and supplied
/// anchor must share a clock. Preparation independently rechecks the anchor's
/// complete nonlinear equations and computes `z = dX/dln(omega)` with the
/// existing exact tangent. Every grid/refinement query starts from the SAME
/// `X_anchor + ln(omega/omega_anchor) * z_anchor`; it never inherits the previous
/// query's answer. The original seeded, physically scaled HB owner corrects it.
/// A per-coordinate RMS correction gate refuses distant solutions.
///
/// Every corrected orbit is independently admitted by the original tangent
/// owner, including base balance, quadrature and retained/omitted derivative
/// equations. The physical observer evaluates `G^T grad H`, not a state proxy;
/// its derivative is the original `G^T Hxx z`. Both keep their shifted-grid
/// output-resolution gates. Their positive Fourier coefficients become peak
/// phasors by `2*conj(c_k)`; the log derivative is divided by the fundamental
/// frequency, not by k*omega. The shared pHS owner locates the amplitude peaks.
///
/// These are local, fixed-order numerical results. Correction limits do not
/// prove a unique branch or exclude two nearby solutions. Near a fold this
/// frequency graph may cease to exist; native singularities/correction limits
/// refuse without regularization or automatic branch switching. Use the
/// existing pseudo-arclength continuation for folds and the separate original
/// Floquet owner for stability. Repeat state order and grid resolution as
/// appropriate: no infinite-tail, complete-peak, or stability claim is added.
/// A phase-fixed survey does not infer frequency-dependent drive filters or
/// time delays; those require their original physical memory/source law.
///
/// Work is one anchor tangent plus at most `survey.maximum_evaluations` original
/// HB/spectrum/tangent queries, each under the supplied owner's caps. Only one
/// query is live; the additional retained anchor derivative has at most 1024
/// complex entries (the tangent owner's hard cap). No orbit cache accumulates.
/// Cancellation is polled between owners and within their existing cancellable
/// operations; one synchronous bounded HB correction cannot be interrupted.
///
/// # Errors
/// Invalid controls, clock/source mismatch, unresolved physical state/output,
/// native solve refusal, predictor disagreement, exhausted query capacity or
/// cancellation. No partial peak list is returned and the anchor is immutable.
#[allow(clippy::too_many_arguments)]
pub fn survey_forced_port_resonances(
    system: &PortHamiltonian,
    input: &PeriodicPortInput,
    anchor: &HbOrbit,
    output_port: usize,
    output_harmonic: usize,
    frequencies_rad_s: &[f64],
    controls: &PhsForcedResonanceControls,
    mut cancelled: impl FnMut() -> bool,
) -> core::result::Result<PhsForcedResonanceSurvey, ResonanceSurveyError<PhsForcedResonanceError>> {
    let source_error = |error| ResonanceSurveyError::Search(ResonanceSearchError::Source(error));
    poll(&mut cancelled).map_err(source_error)?;
    admit(
        system,
        input,
        anchor,
        output_port,
        output_harmonic,
        frequencies_rad_s,
        controls,
    )
    .map_err(source_error)?;
    // Lazy preparation lets the shared survey admit its complete grid/search
    // controls before any original storage or nonlinear-solver work is done.
    let mut anchor_derivative: Option<Vec<C64>> = None;
    let mut maximum_predictor_correction_rms = vec![0.0_f64; system.state_dim()];
    let resonances = survey_resonances(
        frequencies_rad_s,
        controls.survey,
        &mut cancelled,
        |omega, _need_derivative, cancelled| {
            poll(cancelled)?;
            if anchor_derivative.is_none() {
                let mut source =
                    PhsPortSidebandSource::forced(system, input, anchor, &controls.tangent)
                        .map_err(PhsForcedResonanceError::Tangent)?;
                let mut prepared = source
                    .prepare(&controls.tangent_scales, &mut *cancelled)
                    .map_err(PhsForcedResonanceError::Tangent)?;
                anchor_derivative = Some(
                    prepared
                        .forced_frequency_tangent(&mut *cancelled)
                        .map_err(PhsForcedResonanceError::Tangent)?
                        .state_log_frequency_derivative
                        .coefficients,
                );
            }
            let predictor = predict(
                anchor,
                anchor_derivative
                    .as_ref()
                    .ok_or(PhsForcedResonanceError::Invalid(
                        "anchor tangent was not prepared",
                    ))?,
                omega,
                controls.harmonic_balance.harmonics,
            )?;
            let drive = input
                .at_frequency(omega)
                .map_err(PhsForcedResonanceError::Orbit)?;
            poll(cancelled)?;
            let solution = solve_forced_harmonic_balance_scaled_seeded(
                system,
                &drive,
                predictor.packed_seed(&HbAnchor::Forced { omega }),
                &controls.harmonic_balance,
                &controls.harmonic_balance_scales,
            )
            .map_err(PhsForcedResonanceError::Orbit)?;
            poll(cancelled)?;
            assess_correction(
                &predictor,
                &solution.orbit,
                &controls.predictor_rms_tolerances,
                &mut maximum_predictor_correction_rms,
            )?;
            let value =
                observe_port_spectrum(system, &solution.orbit, &controls.spectrum, &mut *cancelled)
                    .map_err(PhsForcedResonanceError::Spectrum)?
                    .peak_phasors[output_port][output_harmonic - 1];
            let mut source =
                PhsPortSidebandSource::forced(system, &drive, &solution.orbit, &controls.tangent)
                    .map_err(PhsForcedResonanceError::Tangent)?;
            let mut prepared = source
                .prepare(&controls.tangent_scales, &mut *cancelled)
                .map_err(PhsForcedResonanceError::Tangent)?;
            let tangent = prepared
                .forced_frequency_tangent(&mut *cancelled)
                .map_err(PhsForcedResonanceError::Tangent)?;
            let log_derivative = tangent.port_log_frequency_derivative
                [(tangent.output_harmonics + output_harmonic) * system.port_dim() + output_port];
            let derivative = C64::new(
                2.0 * (log_derivative.re / omega),
                -2.0 * (log_derivative.im / omega),
            );
            poll(cancelled)?;
            Ok((
                value,
                derivative,
                tangent.state_log_frequency_derivative.relative_residual,
            ))
        },
    )?;
    poll(&mut cancelled).map_err(source_error)?;
    Ok(PhsForcedResonanceSurvey {
        resonances,
        maximum_predictor_correction_rms,
    })
}

fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<()> {
    if cancelled() {
        Err(PhsForcedResonanceError::Cancelled)
    } else {
        Ok(())
    }
}

#[allow(clippy::too_many_arguments)]
fn admit(
    system: &PortHamiltonian,
    input: &PeriodicPortInput,
    anchor: &HbOrbit,
    port: usize,
    harmonic: usize,
    frequencies: &[f64],
    c: &PhsForcedResonanceControls,
) -> Result<()> {
    if anchor.omega.to_bits() != input.omega_rad_s().to_bits()
        || port >= system.port_dim()
        || harmonic == 0
        || harmonic > c.spectrum.harmonics
        || c.spectrum.harmonics != c.tangent.output_harmonics
        || c.harmonic_balance.harmonics != c.tangent.sideband.harmonics
        || c.harmonic_balance.harmonics == 0
        || !c.maximum_log_frequency_distance.is_finite()
        || c.maximum_log_frequency_distance <= 0.0
        || c.maximum_log_frequency_distance > 1.0
        || c.predictor_rms_tolerances.len() != system.state_dim()
        || c.predictor_rms_tolerances
            .iter()
            .any(|v| !v.is_finite() || *v <= 0.0)
    {
        return Err(PhsForcedResonanceError::Invalid(
            "one original positive output harmonic, common source clock/orders and positive physical branch limits required",
        ));
    }
    if frequencies.iter().any(|&omega| {
        let distance = det::ln(omega / anchor.omega);
        !distance.is_finite() || distance.abs() > c.maximum_log_frequency_distance
    }) {
        return Err(PhsForcedResonanceError::Invalid(
            "the complete frequency grid must stay within the authored anchor distance",
        ));
    }
    Ok(())
}

fn predict(anchor: &HbOrbit, derivative: &[C64], omega: f64, h: usize) -> Result<HbOrbit> {
    let n = anchor.coeffs.len();
    let distance = det::ln(omega / anchor.omega);
    // Solver metadata is neither copied nor trusted as evidence for a seed.
    let mut result = HbOrbit {
        omega,
        coeffs: anchor
            .coeffs
            .iter()
            .map(|row| row.iter().take(h + 1).copied().collect())
            .collect(),
        residual_trace: Vec::new(),
        residual: 0.0,
        aft_samples: 0,
        sampling_defect: 0.0,
    };
    for (j, row) in result.coeffs.iter_mut().enumerate() {
        // Native anchor admission already refused any active omitted band.
        row.resize(h + 1, C64::ZERO);
        for (k, coefficient) in row.iter_mut().enumerate() {
            let z = derivative[(h + k) * n + j];
            *coefficient = *coefficient + z.scale(distance);
            if k == 0 {
                // The real source has conjugate symmetry; the tiny imaginary
                // zero-band roundoff of its complex solve is not a real state.
                coefficient.im = 0.0;
            }
            if !coefficient.re.is_finite() || !coefficient.im.is_finite() {
                return Err(PhsForcedResonanceError::Invalid(
                    "anchor predictor overflow",
                ));
            }
        }
    }
    Ok(result)
}

fn assess_correction(
    predictor: &HbOrbit,
    orbit: &HbOrbit,
    tolerances: &[f64],
    maxima: &mut [f64],
) -> Result<()> {
    for (j, (predicted, corrected)) in predictor.coeffs.iter().zip(&orbit.coeffs).enumerate() {
        let mut rms = 0.0;
        for (k, (a, b)) in predicted.iter().zip(corrected).enumerate() {
            let difference = (*b - *a).abs();
            rms = det::hypot(rms, difference);
            if k != 0 {
                rms = det::hypot(rms, difference);
            }
        }
        if !rms.is_finite() || rms > tolerances[j] {
            return Err(PhsForcedResonanceError::PredictorDisagreement {
                coordinate: j,
                correction_rms: rms,
                tolerance: tolerances[j],
            });
        }
        maxima[j] = maxima[j].max(rms);
    }
    Ok(())
}
