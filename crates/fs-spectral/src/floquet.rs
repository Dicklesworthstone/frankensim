//! Complete numerical Floquet spectra and original-flow phase neutrality.
//!
//! This owner interprets an explicitly supplied dimensionless monodromy. It
//! neither constructs an evolution model nor establishes that its base waveform
//! is a periodic orbit. Integration, temporal convergence and physical source
//! equivalence belong to the caller. All eigenpairs and the full basis are
//! checked by the shared nonsymmetric dense service, with no dropped modes.

use crate::dense::{
    DenseSpectrum, DenseSpectrumError, DenseSpectrumOptions, eigenpair_row_residuals,
    nonsymmetric_eigenbasis, nonsymmetric_spectrum,
};
use fs_math::{c64::C64, det};

/// Complete spectral acceptance plus an explicit unit-circle exclusion margin.
#[derive(Clone, Copy, Debug)]
pub struct FloquetControls {
    pub spectrum: DenseSpectrumOptions,
    /// Strictly between zero and one. This numerical guard is not an enclosure
    /// of integration error, nonnormal spectral sensitivity or nonlinear motion.
    pub unit_circle_guard: f64,
    /// Maximum distance to +1 for a candidate autonomous phase multiplier.
    /// Positive and strictly smaller than `unit_circle_guard`, which separates
    /// ALL other multipliers from that same +1 point.
    pub phase_multiplier_tolerance: f64,
    /// Original-row backward residual of M*g=g for the supplied complete scaled
    /// physical flow direction; positive and below one. No nearest root alone
    /// establishes phase neutrality.
    pub maximum_phase_residual: f64,
}

/// Classification of this supplied, numerically resolved finite map only.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum FloquetClassification {
    Stable,
    Unstable,
    Unresolved,
}

#[derive(Clone, Debug, PartialEq)]
pub enum FloquetError {
    Invalid(&'static str),
    Spectrum(DenseSpectrumError),
    Cancelled,
}
impl core::fmt::Display for FloquetError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "Floquet assessment refused: {self:?}")
    }
}
impl std::error::Error for FloquetError {}
impl From<DenseSpectrumError> for FloquetError {
    fn from(error: DenseSpectrumError) -> Self {
        match error {
            DenseSpectrumError::Cancelled => Self::Cancelled,
            other => Self::Spectrum(other),
        }
    }
}
type Result<T> = core::result::Result<T, FloquetError>;

/// Every multiplier and original-coordinate right vector, including phase.
#[derive(Clone, Debug)]
pub struct FloquetSpectrum {
    pub spectrum: DenseSpectrum,
    /// Maximum modulus over ALL original state directions.
    pub spectral_radius: f64,
    /// Index into `spectrum.values`, present only for an independently checked,
    /// isolated autonomous phase mode. The mode itself remains in the spectrum.
    pub phase_multiplier_index: Option<usize>,
    /// Measured original-row neutral-flow residual. None means that no phase
    /// quotient was requested, as required for an externally forced system.
    pub phase_residual: Option<f64>,
    /// Full radius when no quotient was requested, otherwise transverse radius
    /// only if the physical phase direction was resolved and isolated.
    pub effective_spectral_radius: Option<f64>,
    /// Natural logarithm of the effective radius / supplied physical period.
    /// None for an unresolved phase, zero radius, or unrepresentable growth.
    /// No complex logarithm branch, decay-energy Q or modal damping fit is inferred.
    pub dominant_growth_rate_per_s: Option<f64>,
    pub classification: FloquetClassification,
}

/// Compute and assess the complete spectrum of an ACTUAL supplied monodromy.
///
/// The supplied period is in seconds; the matrix is dimensionless in the
/// caller's disclosed state chart, and real growth rates are in inverse seconds.
/// These numerical results mint no admission or spectral-truth authority.
///
/// `phase_flow=None` keeps every multiplier in the stability radius, including
/// for a zero-amplitude externally forced response. `Some(g)` declares an
/// autonomous nonstationary physical flow in exactly the same scaled state
/// chart as the matrix. Every component, including constitutive memory, belongs
/// in g. The full M*g=g equation, one near-+1 multiplier and separation of all
/// other roots must pass before a transverse radius is reported. A stationary
/// zero flow is invalid; multiple neutral directions remain unresolved.
///
/// The shared dense owner's payload/work caps also cover this O(n^2) phase
/// assessment and O(n) normalized flow scratch. Input matrix/flow and caller
/// storage are excluded. Complete-basis conditioning and backward residuals
/// are numerical gates, not a forward-error or invariant-subspace enclosure.
///
/// # Errors
/// Invalid matrix/physical period/controls/phase declaration, shared spectral
/// or resource refusal, unrepresentable arithmetic, or cancellation. No partial
/// accepted spectrum is returned on failure.
pub fn analyze_monodromy(
    matrix: &[C64],
    dimension: usize,
    period: f64,
    phase_flow: Option<&[C64]>,
    controls: FloquetControls,
    mut cancelled: impl FnMut() -> bool,
) -> Result<FloquetSpectrum> {
    poll(&mut cancelled)?;
    admit(dimension, period, phase_flow, controls)?;
    let spectrum = nonsymmetric_spectrum(matrix, dimension, controls.spectrum, &mut cancelled)?;
    assess(
        matrix,
        dimension,
        period,
        phase_flow,
        controls,
        spectrum,
        &mut cancelled,
    )
}

/// Assess previously computed COMPLETE monodromy multipliers without repeating
/// QR discovery. Original eigenvector residuals and complete-basis conditioning
/// are still checked by `dense::nonsymmetric_eigenbasis` through the same native
/// SVD owner. Input order is immaterial; the returned spectrum is canonical and
/// keeps multiplicity. Supplying an incomplete/false list does not bypass the
/// original matrix checks. Other semantics and refusals match `analyze_monodromy`.
#[allow(clippy::too_many_arguments)]
pub fn assess_monodromy_with_multipliers(
    matrix: &[C64],
    dimension: usize,
    period: f64,
    multipliers: &[C64],
    phase_flow: Option<&[C64]>,
    controls: FloquetControls,
    mut cancelled: impl FnMut() -> bool,
) -> Result<FloquetSpectrum> {
    poll(&mut cancelled)?;
    admit(dimension, period, phase_flow, controls)?;
    let spectrum = nonsymmetric_eigenbasis(
        matrix,
        dimension,
        multipliers,
        controls.spectrum,
        &mut cancelled,
    )?;
    assess(
        matrix,
        dimension,
        period,
        phase_flow,
        controls,
        spectrum,
        &mut cancelled,
    )
}

fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<()> {
    if cancelled() {
        Err(FloquetError::Cancelled)
    } else {
        Ok(())
    }
}
fn finite(z: C64) -> bool {
    z.re.is_finite() && z.im.is_finite()
}

fn admit(n: usize, period: f64, phase: Option<&[C64]>, controls: FloquetControls) -> Result<()> {
    if n == 0
        || n > 1024
        || !period.is_finite()
        || period <= 0.0
        || !controls.unit_circle_guard.is_finite()
        || controls.unit_circle_guard <= 0.0
        || controls.unit_circle_guard >= 1.0
        || !controls.phase_multiplier_tolerance.is_finite()
        || controls.phase_multiplier_tolerance <= 0.0
        || controls.phase_multiplier_tolerance >= controls.unit_circle_guard
        || !controls.maximum_phase_residual.is_finite()
        || controls.maximum_phase_residual <= 0.0
        || controls.maximum_phase_residual >= 1.0
        || phase.is_some_and(|g| {
            g.len() != n || g.iter().any(|&z| !finite(z)) || g.iter().all(|&z| z == C64::ZERO)
        })
    {
        return Err(FloquetError::Invalid(
            "positive period, guarded Floquet controls and complete nonzero phase flow required",
        ));
    }
    Ok(())
}

#[allow(clippy::too_many_arguments)]
fn assess(
    matrix: &[C64],
    n: usize,
    period: f64,
    phase_flow: Option<&[C64]>,
    controls: FloquetControls,
    spectrum: DenseSpectrum,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<FloquetSpectrum> {
    let spectral_radius = spectrum.values.iter().map(|z| z.abs()).fold(0.0, f64::max);
    if !spectral_radius.is_finite() {
        return Err(FloquetError::Invalid("Floquet radius is not representable"));
    }
    let mut phase_residual = None;
    let mut phase_multiplier_index = None;
    let effective_spectral_radius = if let Some(flow) = phase_flow {
        // Scale the flow BEFORE the residual products so tiny physical units do
        // not disappear and large but finite flow components do not overflow.
        let scale = flow
            .iter()
            .map(|z| z.re.abs().max(z.im.abs()))
            .fold(0.0, f64::max);
        let mut normalized = Vec::with_capacity(n);
        for &value in flow {
            poll(cancelled)?;
            let value_scaled = C64::new(value.re / scale, value.im / scale);
            if !finite(value_scaled)
                || (value.re != 0.0 && value_scaled.re == 0.0)
                || (value.im != 0.0 && value_scaled.im == 0.0)
            {
                return Err(FloquetError::Invalid(
                    "normalized full phase flow is not representable",
                ));
            }
            normalized.push(value_scaled);
        }
        let residual =
            eigenpair_row_residuals(matrix, n, C64::from_re(1.0), &normalized, &mut *cancelled)?
                .into_iter()
                .fold(0.0, f64::max);
        phase_residual = Some(residual);
        let mut near = None;
        let mut count = 0;
        for (index, &value) in spectrum.values.iter().enumerate() {
            poll(cancelled)?;
            if (value - C64::from_re(1.0)).abs() <= controls.phase_multiplier_tolerance {
                count += 1;
                near = Some(index);
            }
        }
        if residual <= controls.maximum_phase_residual
            && count == 1
            && spectrum.values.iter().enumerate().all(|(i, &value)| {
                Some(i) == near || (value - C64::from_re(1.0)).abs() > controls.unit_circle_guard
            })
        {
            phase_multiplier_index = near;
            Some(
                spectrum
                    .values
                    .iter()
                    .enumerate()
                    .filter(|(i, _)| Some(*i) != near)
                    .map(|(_, z)| z.abs())
                    .fold(0.0, f64::max),
            )
        } else {
            None
        }
    } else {
        Some(spectral_radius)
    };
    let classification = match effective_spectral_radius {
        Some(radius) if radius < 1.0 - controls.unit_circle_guard => FloquetClassification::Stable,
        Some(radius) if radius > 1.0 + controls.unit_circle_guard => {
            FloquetClassification::Unstable
        }
        _ => FloquetClassification::Unresolved,
    };
    let dominant_growth_rate_per_s = effective_spectral_radius.and_then(|radius| {
        if radius <= 0.0 {
            return None;
        }
        let growth = det::ln(radius) / period;
        growth.is_finite().then_some(growth)
    });
    poll(cancelled)?;
    Ok(FloquetSpectrum {
        spectrum,
        spectral_radius,
        phase_multiplier_index,
        phase_residual,
        effective_spectral_radius,
        dominant_growth_rate_per_s,
        classification,
    })
}
