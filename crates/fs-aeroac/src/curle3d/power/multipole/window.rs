//! Real outgoing radiation on an explicit retarded-time interval.
//!
//! At x = center + R n, let tau = t - R/c on the original source clock.
//! In the radiation-zone limit the parent's sqrt(watt) source amplitude A
//! gives P(tau) = 2 <Re sum_l [-i A_l(n) exp(-i omega_l tau)] ^ 2>.
//! The leading radiative velocity is p/(rho c); this is the limiting outward
//! pressure/velocity flux, not a local near-field p^2 intensity model.
//!
//! Difference- and sum-frequency products have different spatial phases:
//! k_a(s_a-center) - k_b(s_b-center) and
//! k_a(s_a-center) + k_b(s_b-center), respectively. Dropping the second
//! product, or replacing both by one same-frequency separation, is incorrect
//! for a finite window. Both use the original analytic angular moments.
//! The shared harmonic-window owner supplies the real temporal product.

use super::{
    AeroacError, CompactMultipolePhasor, ComplexDd, Dd, FluidAcousticState, MultipoleSum,
    PreparedMultipoleRadiatedPower, RadiatingMultipole, finite_dd, invalid, nonfinite,
    validate_power,
};
use fs_fft::harmonic_window::{HarmonicWindow, HarmonicWindowProfile, WindowError, angular_frequency};
use std::ops::Range;

mod kernel;
use kernel::{AngularKernel, SpatialCorrection, colocated_product, prepare_correction, transform};

/// Authority of the full-sphere real-field energy and outgoing wave momentum.
pub const MULTIPOLE_WINDOW_SCOPE: &str = "Finite uniform retarded-time integral of prescribed coherent Q/F/Lighthill-stress radiation in the original homogeneous stationary lossless free space, with all spatial, beat and real-field sum-frequency terms; radiation-zone enclosing-sphere flux, not local near-field intensity, finite-radius storage balance, mechanical impulse, source feedback, stochastic PSD inference, moving sources, mean flow, boundaries or measured-SPL qualification.";

/// One original physical source bank at its owner's positive frequency.
///
/// Every bank shares one Cartesian frame and emission clock. Sources retain
/// peak Q [m^3/s], force ON the fluid [N] and symmetric stress [N m]. Positions,
/// material, normalization and frequency come only from the borrowed owner.
/// The caller must not count a physical mechanism twice in different banks.
#[derive(Clone, Copy, Debug)]
pub struct MultipoleRadiationLine<'p> {
    pub operator: &'p PreparedMultipoleRadiatedPower,
    pub sources: &'p [CompactMultipolePhasor],
}

/// Preparation and repeated-window work admitted before allocation.
#[derive(Clone, Copy, Debug)]
pub struct RadiatedEnergyBudget {
    /// Supplied banks, including repeated frequencies and silent banks.
    pub max_lines: usize,
    /// Sum of the original source counts in all supplied banks.
    pub max_source_terms: usize,
    /// L^2 + P^2 for L supplied banks and P original source terms. Bounds
    /// frequency/source grouping, angular pairs and every window application.
    pub max_pair_work: usize,
    /// Conservative peak bound for owned buffers and angular/product scratch.
    /// Uses supplied counts before coalescing. Borrowed
    /// source owners/banks, allocator overhead and call stacks are separate.
    pub max_workspace_bytes: usize,
}

/// The radiation-zone flux integrated over one physical sphere-time window.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct WindowedRadiatedEnergy {
    /// tau = t - R/c for spheres centered at the prepared center.
    pub window: HarmonicWindow,
    pub mean_power_w: f64,
    pub mean_power_first_moment_w: [f64; 3],
    pub radiated_energy_j: f64,
    /// Integral of direction times differential outgoing energy [J].
    /// Its norm is at most radiated_energy_j in the exact model.
    pub energy_first_moment_j: [f64; 3],
    /// energy_first_moment_j / c [N s = kg m/s], not mechanical recoil.
    pub wave_momentum_n_s: [f64; 3],
    pub scope: &'static str,
}

/// A source-bound, deterministic finite-spectrum radiation observer.
///
/// Equal frequencies AND positions coalesce physical source terms before
/// normalization. A global colocated baseline and separate angular remainders
/// retain weak cancellation. Preparation costs O(L^2 + P^2), contracts spatial
/// moments once, and retains O(L + F^2) values for F distinct frequencies.
/// Repeated windows cost O(F^2), allocate nothing, and reuse one temporal
/// kernel across energy and all three directional components.
///
/// The explicit center defines the enclosing spheres, not a new source clock.
/// Translating sources and center together preserves the result. Moving just
/// the sources changes arrival times across the sphere and can change a finite
/// window's energy. No start-up transient or stochastic source law is inferred.
pub struct PreparedMultipoleRadiatedEnergy<'p> {
    lines: Vec<MultipoleRadiationLine<'p>>,
    frequencies: Vec<Frequency>,
    pairs: Vec<SpatialCorrection>,
    acoustic: FluidAcousticState,
    sphere_center_m: [f64; 3],
    pair_work: usize,
    workspace_bytes: usize,
}

struct FrequencyBuild<'p> {
    operator: &'p PreparedMultipoleRadiatedPower,
    sum: MultipoleSum,
}

struct SourceBuild {
    frequency: usize,
    position_m: [f64; 3],
    sum: MultipoleSum,
}

struct Frequency {
    hz: f64,
    total: RadiatingMultipole,
    atoms: Range<usize>,
}

struct Atom {
    wavevector: [Dd; 3],
    source: RadiatingMultipole,
}

impl<'p> PreparedMultipoleRadiatedEnergy<'p> {
    /// Bind original source banks and the center of the observing spheres.
    ///
    /// At least one bank supplies the original medium and frequency clock;
    /// an empty or zero-amplitude source bank represents darkness. Every
    /// owner must have exactly the same material and free-space domain.
    pub fn new(
        lines: &[MultipoleRadiationLine<'p>],
        sphere_center_m: [f64; 3],
        budget: RadiatedEnergyBudget,
    ) -> Result<Self, AeroacError> {
        Self::new_cancellable(lines, sphere_center_m, budget, || false)
    }

    pub fn new_cancellable(
        lines: &[MultipoleRadiationLine<'p>],
        sphere_center_m: [f64; 3],
        budget: RadiatedEnergyBudget,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, AeroacError> {
        poll(&mut cancelled)?;
        if lines.is_empty()
            || lines.len() > budget.max_lines
            || sphere_center_m.iter().any(|v| !v.is_finite())
        {
            return Err(invalid("radiated energy needs an admitted original line bank and sphere center"));
        }
        let acoustic = lines[0].operator.compact.acoustic;
        let mut source_terms = 0_usize;
        for line in lines {
            poll(&mut cancelled)?;
            line.operator.check_count(line.sources)?;
            if line.operator.compact.image_group_size != 1
                || line.operator.compact.acoustic != acoustic
            {
                return Err(invalid("radiated energy requires the same original lossless free-space medium"));
            }
            source_terms = source_terms.checked_add(line.sources.len())
                .ok_or(invalid("radiated-energy source count overflows"))?;
        }
        let line_pairs = lines.len().checked_mul(lines.len())
            .ok_or(invalid("radiated-energy line-pair work overflows"))?;
        let pair_work = source_terms.checked_mul(source_terms)
            .and_then(|n| n.checked_add(line_pairs))
            .ok_or(invalid("radiated-energy source-pair work overflows"))?;
        if source_terms > budget.max_source_terms || pair_work > budget.max_pair_work {
            return Err(invalid("radiated energy exceeds its source/pair budget"));
        }
        let pair_capacity = line_pairs.checked_add(lines.len())
            .map(|n| n / 2).ok_or(invalid("radiated-energy pair storage overflows"))?;
        let workspace_bytes = workspace_bytes(lines.len(), source_terms, pair_capacity)?;
        if workspace_bytes > budget.max_workspace_bytes {
            return Err(invalid("radiated energy exceeds its workspace budget"));
        }

        let mut retained = reserve(lines.len())?;
        retained.extend_from_slice(lines);
        let mut frequency_build: Vec<FrequencyBuild<'p>> = reserve(lines.len())?;
        let mut source_build: Vec<SourceBuild> = reserve(source_terms)?;
        for line in lines {
            poll(&mut cancelled)?;
            let hz = line.operator.compact.frequency_hz;
            let mut found = None;
            for (i, frequency) in frequency_build.iter().enumerate() {
                poll(&mut cancelled)?;
                if frequency.operator.compact.frequency_hz == hz {
                    found = Some(i);
                    break;
                }
            }
            let frequency = found.unwrap_or(frequency_build.len());
            if found.is_none() {
                frequency_build.push(FrequencyBuild { operator: line.operator, sum: MultipoleSum::default() });
            }
            // Keep the original source order in the global baseline, even
            // when position grouping has a different traversal order.
            for &source in line.sources {
                poll(&mut cancelled)?;
                frequency_build[frequency].sum.add(source)?;
            }
            for group in &line.operator.compact.groups {
                let mut found = None;
                for (i, atom) in source_build.iter().enumerate() {
                    poll(&mut cancelled)?;
                    if atom.frequency == frequency && atom.position_m == group.position_m {
                        found = Some(i);
                        break;
                    }
                }
                let atom = found.unwrap_or(source_build.len());
                if found.is_none() {
                    source_build.push(SourceBuild {
                        frequency,
                        position_m: group.position_m,
                        sum: MultipoleSum::default(),
                    });
                }
                for &source_index in &line.operator.compact.members[group.members.clone()] {
                    poll(&mut cancelled)?;
                    source_build[atom].sum.add(line.sources[source_index])?;
                }
            }
        }

        let mut frequencies = reserve(frequency_build.len())?;
        let mut atoms = reserve(source_build.len())?;
        for (index, frequency) in frequency_build.iter().enumerate() {
            poll(&mut cancelled)?;
            let hz = frequency.operator.compact.frequency_hz;
            // Difference-frequency spatial phases retain their authored Hz
            // before the common 2*pi/c scaling; they are not differences of
            // already rounded binary64 wave numbers.
            let k = angular_frequency(hz).map_err(window_error)?
                / Dd::from_f64(acoustic.sound_speed_m_s());
            let begin = atoms.len();
            for source in &source_build {
                poll(&mut cancelled)?;
                if source.frequency != index {
                    continue;
                }
                let wavevector = core::array::from_fn(|axis| {
                    (Dd::from_f64(source.position_m[axis])
                        - Dd::from_f64(sphere_center_m[axis])) * k
                });
                if wavevector.iter().any(|v| !finite_dd(*v)) {
                    return Err(nonfinite("radiated-energy spatial phase"));
                }
                atoms.push(Atom {
                    wavevector,
                    source: transform(frequency.operator.normalized(source.sum)?, ComplexDd::negative_i),
                });
            }
            frequencies.push(Frequency {
                hz,
                total: transform(frequency.operator.normalized(frequency.sum)?, ComplexDd::negative_i),
                atoms: begin..atoms.len(),
            });
        }
        let actual_pairs = frequencies.len() * (frequencies.len() + 1) / 2;
        let mut pairs = reserve(actual_pairs)?;
        for (i, a) in frequencies.iter().enumerate() {
            for (j, b) in frequencies[..=i].iter().enumerate() {
                poll(&mut cancelled)?;
                pairs.push(prepare_correction(
                    &atoms[a.atoms.clone()], &atoms[b.atoms.clone()], i == j, &mut cancelled,
                )?);
            }
        }
        poll(&mut cancelled)?;
        Ok(Self {
            lines: retained,
            frequencies,
            pairs,
            acoustic,
            sphere_center_m,
            pair_work,
            workspace_bytes,
        })
    }

    pub fn source_lines(&self) -> &[MultipoleRadiationLine<'p>] { &self.lines }
    pub const fn acoustic_state(&self) -> FluidAcousticState { self.acoustic }
    pub const fn sphere_center_m(&self) -> [f64; 3] { self.sphere_center_m }
    pub fn frequency_count(&self) -> usize { self.frequencies.len() }
    pub const fn preparation_pair_work(&self) -> usize { self.pair_work }
    /// Conservative preallocation bound, including temporary preparation data.
    pub const fn workspace_bytes(&self) -> usize { self.workspace_bytes }

    /// Integrate real outgoing energy, preserving both sum and beat terms.
    /// Only uniform windows are unweighted physical energy intervals.
    pub fn observe(&self, window: HarmonicWindow) -> Result<WindowedRadiatedEnergy, AeroacError> {
        self.observe_cancellable(window, || false)
    }

    pub fn observe_cancellable(
        &self,
        window: HarmonicWindow,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<WindowedRadiatedEnergy, AeroacError> {
        poll(&mut cancelled)?;
        if window.profile() != HarmonicWindowProfile::Uniform {
            return Err(invalid("radiated energy requires a uniform physical time interval"));
        }
        let mut mean = [Dd::ZERO; 4];
        let mut pair = 0;
        for (i, a) in self.frequencies.iter().enumerate() {
            for (j, b) in self.frequencies[..=i].iter().enumerate() {
                poll(&mut cancelled)?;
                let temporal = window.product_kernel_hz(a.hz, b.hz).map_err(window_error)?;
                let baseline = colocated_product(temporal, a.total, b.total)?;
                let correction = self.pairs[pair].apply(temporal)?;
                // The real pressure normalization contributes 2. A distinct
                // frequency pair contributes both orders of the real product.
                let factor = Dd::from_f64(if i == j { 2.0 } else { 4.0 });
                for axis in 0..4 {
                    mean[axis] = mean[axis] + factor * (baseline[axis] + correction[axis]);
                }
                pair += 1;
            }
        }
        poll(&mut cancelled)?;
        finish(window, mean, self.acoustic.sound_speed_m_s())
    }
}

fn workspace_bytes(lines: usize, sources: usize, pairs: usize) -> Result<usize, AeroacError> {
    let mut bytes = core::mem::size_of::<PreparedMultipoleRadiatedEnergy<'_>>()
        + core::mem::size_of::<AngularKernel>()
        + 2 * core::mem::size_of::<SpatialCorrection>()
        + core::mem::size_of::<fs_fft::harmonic_window::HarmonicProduct>()
        + core::mem::size_of::<WindowedRadiatedEnergy>();
    for (count, width) in [
        (lines, core::mem::size_of::<MultipoleRadiationLine<'_>>()),
        (lines, core::mem::size_of::<FrequencyBuild<'_>>()),
        (lines, core::mem::size_of::<Frequency>()),
        (sources, core::mem::size_of::<SourceBuild>()),
        (sources, core::mem::size_of::<Atom>()),
        (pairs, core::mem::size_of::<SpatialCorrection>()),
    ] {
        bytes = count.checked_mul(width).and_then(|n| bytes.checked_add(n))
            .ok_or(invalid("radiated-energy workspace overflows"))?;
    }
    Ok(bytes)
}

fn reserve<T>(count: usize) -> Result<Vec<T>, AeroacError> {
    let mut result = Vec::new();
    result.try_reserve_exact(count)
        .map_err(|_| invalid("radiated-energy allocation exceeds available memory"))?;
    Ok(result)
}

fn finish(window: HarmonicWindow, mean: [Dd; 4], c: f64) -> Result<WindowedRadiatedEnergy, AeroacError> {
    validate_power(mean[0])?;
    let duration = Dd::from_f64(window.duration_s());
    let energy = mean.map(|value| value * duration);
    let momentum = [energy[1], energy[2], energy[3]].map(|value| value / Dd::from_f64(c));
    if mean.iter().chain(&energy).chain(&momentum).any(|v| !finite_dd(*v))
        || (Dd::ZERO.lt(mean[0]) && !Dd::ZERO.lt(energy[0]))
    {
        return Err(nonfinite("finite-window radiated energy or wave momentum"));
    }
    let result = WindowedRadiatedEnergy {
        window,
        mean_power_w: mean[0].to_f64(),
        mean_power_first_moment_w: [mean[1], mean[2], mean[3]].map(Dd::to_f64),
        radiated_energy_j: energy[0].to_f64(),
        energy_first_moment_j: [energy[1], energy[2], energy[3]].map(Dd::to_f64),
        wave_momentum_n_s: momentum.map(Dd::to_f64),
        scope: MULTIPOLE_WINDOW_SCOPE,
    };
    if [result.mean_power_w, result.radiated_energy_j].iter()
        .chain(&result.mean_power_first_moment_w)
        .chain(&result.energy_first_moment_j)
        .chain(&result.wave_momentum_n_s).any(|v| !v.is_finite())
    {
        return Err(nonfinite("rounded finite-window radiation result"));
    }
    Ok(result)
}

fn window_error(error: WindowError) -> AeroacError {
    match error {
        WindowError::InvalidWindow => invalid("invalid retarded-time radiation window"),
        WindowError::InvalidFrequency => invalid("invalid radiated-energy frequency"),
        WindowError::NonFiniteAmplitude | WindowError::Unrepresentable =>
            nonfinite("radiated-energy temporal product"),
    }
}

fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), AeroacError> {
    if cancelled() { Err(AeroacError::Cancelled { what: "finite-window multipole radiation" }) }
    else { Ok(()) }
}
