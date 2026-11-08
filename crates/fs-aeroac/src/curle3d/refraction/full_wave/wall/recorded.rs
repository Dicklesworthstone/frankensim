//! Recorded mixed-source sound at one frequency-dependent locally reacting wall.
//!
//! All source mechanisms share the original Welch clock, window and factors.
//! Full pressure and vector velocity interfere before quadratic observations.
//! Frequency-band observation adds resolved bin powers, not complex amplitudes
//! at unrelated frequencies, and retains only one frequency's transfer banks.

use super::super::array::{FullWaveCompactSourceLimits, remaining_budget};
use super::super::quadrupole::array::FullWaveQuadrupoleSourceLimits;
use super::super::{QuadratureLimits, admit, force, poll};
use super::{
    ImpedanceWallFieldKind, PlanarImpedanceWall, PreparedImpedanceWallCompactSourceObserver,
    PreparedImpedanceWallQuadrupoleObserver,
};
use crate::AeroacError;
use crate::curle3d::statistics::{StatisticsSum, invalid};
use crate::curle3d::{CoherenceBudget, CompactFieldStatistics};
use core::ops::Range;
use fs_fft::welch::WelchCrossSpectrum;
use fs_math::c64::C64;

const REFLECTED_SCOPE: &str = "Original synchronized stationary compact volume flows \
[m³/s], forces ON fluid [N] and symmetric stress moments [N m]; reflected-only full \
spherical pressure and vector velocity at one passive frequency-dependent locally \
reacting infinite wall in a homogeneous lossless fluid. Original source/component \
cross terms survive within every Welch factor. Bin-integrated second moments and \
stationary resolved-bin sums retain sampling, window, leakage and frequency-resolution \
limitations. Native quadrature estimates are not interval certificates or relative \
guarantees at interference nulls. No direct field, deterministic cross-frequency phases, \
transient waveform, lossless compliant pole, inferred independence, unresolved SGS \
sound, finite edges, mean flow, bulk loss, source feedback or measured-SPL qualification.";

const TOTAL_SCOPE: &str = "Original synchronized stationary compact volume flows \
[m³/s], forces ON fluid [N] and symmetric stress moments [N m]; coherent direct-plus- \
reflected spherical pressure and vector velocity at one passive frequency-dependent \
locally reacting infinite wall in a homogeneous lossless fluid. Every source, tensor \
and path cross term survives within its original Welch factor. Bin-integrated second \
moments and stationary resolved-bin sums retain sampling, window, leakage and frequency- \
resolution limitations. No relative accuracy guarantee at interference nulls, singular \
total self-field, interval certificate, deterministic cross-frequency phases, transient \
waveform, lossless compliant pole, inferred independence, unresolved SGS sound, finite \
edges, mean flow, bulk loss, source feedback or measured-SPL qualification.";

pub(super) const fn scope(field_kind: ImpedanceWallFieldKind) -> &'static str {
    match field_kind {
        ImpedanceWallFieldKind::Reflected => REFLECTED_SCOPE,
        ImpedanceWallFieldKind::Total => TOTAL_SCOPE,
    }
}

pub(super) fn same_geometry(a: PlanarImpedanceWall, b: PlanarImpedanceWall) -> bool {
    a.fluid() == b.fluid() && a.point_m() == b.point_m() && a.outward_normal() == b.outward_normal()
}

/// Native per-bin source limits plus aggregate work for the entire band.
#[derive(Clone, Copy, Debug)]
pub struct ImpedanceWallWelchBandLimits {
    pub compact: FullWaveCompactSourceLimits,
    pub quadrupole: FullWaveQuadrupoleSourceLimits,
    pub maximum_bins: usize,
    /// Positive cap shared by every source mechanism and frequency. Exact
    /// rigid/release images consume no spectral quadrature evaluations.
    pub maximum_total_evaluations: usize,
    /// `max_modes` bounds bins * recorded segments; `max_source_mode_pairs`
    /// bounds (compact sources + stress cells) * bins * recorded segments.
    /// Overlapping segments are algebraic factors, not independent samples.
    pub coherence: CoherenceBudget,
}

/// Stationary power and full p/u moments integrated over selected original bins.
#[derive(Clone, Debug)]
pub struct ImpedanceWallWelchBand {
    /// Included native bin indices; DC and an even-record Nyquist bin refuse.
    pub bins: Range<usize>,
    pub bin_width_hz: f64,
    /// Center frequencies of the first and last included bins, not band edges.
    pub center_frequency_range_hz: [f64; 2],
    pub statistics: CompactFieldStatistics,
    pub evaluations: usize,
}

/// Propagate a joint original spectrum through its wall law at every selected
/// frequency. The callback is invoked once per bin, in increasing frequency,
/// and returns the original owner's specific impedance [Pa s/m], fluid and
/// plane. Only impedance may vary: a changed fluid or plane refuses. Native
/// wall admission rejects active loads and unsupported lossless surface poles.
/// No interpolation, material fit, damping floor or causality claim is added.
///
/// Channels are 4N `[Q,Fx,Fy,Fz]` blocks followed by 6M physical symmetric
/// `[xx,yy,zz,xy,xz,yz]` stress blocks. Either group may be empty; their combined
/// channel count must match the one synchronized spectrum. The supplied source
/// arrays bind its original Cartesian positions and order; source mechanisms
/// must be disjoint, with no permeable-surface double counting.
///
/// Every bin uses its exact original frequency and its own complete wall Q/F/
/// stress kernels. All complex fields add within each recorded factor BEFORE
/// pressure power, velocity covariance and signed active/reactive intensity.
/// Distinct bins add as stationary resolved-spectrum moments; they are never
/// assigned deterministic relative phase. This is not a finite-time waveform
/// or an exact finite-window beat calculation. Resolve both source spectrum
/// and wall dispersion with the original segment length: rapidly varying
/// transfer inside one bin and spectral leakage require frequency refinement.
/// DC, Nyquist, missing channels and frequencies outside the record refuse.
///
/// Only one frequency's cached transfers are retained at a time. All native
/// per-source/per-bin caps and the combined quadrature and factor budgets
/// apply. Cancellation or any refusal returns no partial band. The immutable
/// recording is not transformed, resampled, rewindowed or copied.
#[allow(clippy::too_many_arguments)]
pub fn observe_welch_band(
    spectrum: &WelchCrossSpectrum,
    bins: Range<usize>,
    compact_source_positions_m: &[[f64; 3]],
    stress_source_positions_m: &[[f64; 3]],
    observer_m: [f64; 3],
    field_kind: ImpedanceWallFieldKind,
    mut wall_at_frequency: impl FnMut(f64) -> Result<PlanarImpedanceWall, AeroacError>,
    limits: ImpedanceWallWelchBandLimits,
    mut cancelled: impl FnMut() -> bool,
) -> Result<ImpedanceWallWelchBand, AeroacError> {
    poll(&mut cancelled)?;
    admit(limits.compact.monopole)?;
    force::admit(limits.compact.force)?;
    QuadratureLimits::from(limits.quadrupole.quadrupole).admit()?;
    let count = bins
        .end
        .checked_sub(bins.start)
        .filter(|&n| n > 0 && n <= limits.maximum_bins)
        .ok_or_else(|| invalid("recorded wall band needs a bounded nonempty bin range"))?;
    let compact_count = compact_source_positions_m.len();
    let stress_count = stress_source_positions_m.len();
    let sources = compact_count
        .checked_add(stress_count)
        .ok_or_else(|| invalid("recorded wall band source count overflow"))?;
    let channels = compact_count
        .checked_mul(4)
        .and_then(|n| stress_count.checked_mul(6).and_then(|m| n.checked_add(m)))
        .ok_or_else(|| invalid("recorded wall band channel count overflow"))?;
    let factors = count
        .checked_mul(spectrum.segment_count())
        .ok_or_else(|| invalid("recorded wall band factor count overflow"))?;
    if sources == 0
        || channels != spectrum.channel_count()
        || compact_count > limits.compact.max_sources
        || stress_count > limits.quadrupole.max_sources
        || limits.maximum_total_evaluations == 0
        || limits.compact.maximum_total_evaluations == 0
        || limits.quadrupole.maximum_total_evaluations == 0
        || factors > limits.coherence.max_modes
        || sources
            .checked_mul(factors)
            .is_none_or(|n| n > limits.coherence.max_source_mode_pairs)
    {
        return Err(invalid(
            "recorded wall band source, quadrature or combined factor budget",
        ));
    }
    // Contiguous, increasing bins need only endpoint admission. This excludes
    // DC/Nyquist before invoking the material law or preparing any transfer.
    let mut frequencies = [0.0; 2];
    for (value, bin) in frequencies.iter_mut().zip([bins.start, bins.end - 1]) {
        *value = spectrum.frequency_hz(bin).map_err(AeroacError::Welch)?;
        spectrum
            .peak_bin_hz(bin, *value)
            .map_err(AeroacError::Welch)?;
    }
    let mut statistics = StatisticsSum::default();
    let mut evaluations = 0usize;
    let mut geometry = None;
    for bin in bins.clone() {
        poll(&mut cancelled)?;
        let frequency_hz = spectrum.frequency_hz(bin).map_err(AeroacError::Welch)?;
        let wall = wall_at_frequency(frequency_hz)?;
        if geometry.is_some_and(|original| !same_geometry(original, wall)) {
            return Err(invalid(
                "recorded wall band changed the original fluid or plane",
            ));
        }
        geometry = Some(wall);
        let exact = wall.impedance_pa_s_m().is_none() || wall.impedance_pa_s_m() == Some(C64::ZERO);
        let mut compact_limits = limits.compact;
        compact_limits.maximum_total_evaluations = remaining_budget(
            compact_limits.maximum_total_evaluations,
            limits.maximum_total_evaluations - evaluations,
            exact || compact_count == 0,
        )?;
        let compact = match field_kind {
            ImpedanceWallFieldKind::Reflected => {
                PreparedImpedanceWallCompactSourceObserver::new_reflected(
                    wall,
                    frequency_hz,
                    observer_m,
                    compact_source_positions_m,
                    compact_limits,
                    &mut cancelled,
                )
            }
            ImpedanceWallFieldKind::Total => PreparedImpedanceWallCompactSourceObserver::new_total(
                wall,
                frequency_hz,
                observer_m,
                compact_source_positions_m,
                compact_limits,
                &mut cancelled,
            ),
        }?;
        evaluations = spent(
            evaluations,
            compact.evaluations(),
            limits.maximum_total_evaluations,
        )?;
        let mut stress_limits = limits.quadrupole;
        stress_limits.maximum_total_evaluations = remaining_budget(
            stress_limits.maximum_total_evaluations,
            limits.maximum_total_evaluations - evaluations,
            exact || stress_count == 0,
        )?;
        let stress = match field_kind {
            ImpedanceWallFieldKind::Reflected => {
                PreparedImpedanceWallQuadrupoleObserver::new_reflected(
                    wall,
                    frequency_hz,
                    observer_m,
                    stress_source_positions_m,
                    stress_limits,
                    &mut cancelled,
                )
            }
            ImpedanceWallFieldKind::Total => PreparedImpedanceWallQuadrupoleObserver::new_total(
                wall,
                frequency_hz,
                observer_m,
                stress_source_positions_m,
                stress_limits,
                &mut cancelled,
            ),
        }?;
        evaluations = spent(
            evaluations,
            stress.evaluations(),
            limits.maximum_total_evaluations,
        )?;
        compact.accumulate_welch_bin(
            Some(&stress),
            spectrum,
            bin,
            limits.coherence,
            &mut statistics,
            &mut cancelled,
        )?;
    }
    poll(&mut cancelled)?;
    Ok(ImpedanceWallWelchBand {
        bins,
        bin_width_hz: spectrum.bin_width_hz(),
        center_frequency_range_hz: frequencies,
        statistics: statistics.finish(scope(field_kind))?,
        evaluations,
    })
}

fn spent(previous: usize, current: usize, maximum: usize) -> Result<usize, AeroacError> {
    previous
        .checked_add(current)
        .filter(|&n| n <= maximum)
        .ok_or_else(|| invalid("recorded wall band exhausted its shared quadrature budget"))
}
