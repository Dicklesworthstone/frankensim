//! Whole-sphere radiated power of coherent compact 3D source arrays.
//!
//! In the existing negative-time convention, the far pressure is
//! -i k exp(i k R)/(4 pi R) sum_j (rho c Q_j + n.F_j) exp(-i k n.s_j).
//! Integrating its outward energy flux gives
//! P = k^2/(8 pi rho c) <|sum_j (rho c Q_j + n.F_j) exp(-i k n.s_j)|^2>_sphere.
//! Scalar, vector and tensor spherical moments retain every separated-source
//! and monopole/dipole interference term without angular quadrature.
//!
//! This is cycle-average power in a homogeneous LOSSLESS exterior. It is not
//! local pressure squared, the sum of individual source powers, or a lossy
//! medium's emitted power. Prescribed compact-source and common-clock
//! assumptions remain those of the original field observer.

use super::statistics::{finite, invalid, nonfinite};
use super::{CompactSourcePhasor, PreparedCompactSourceObserver};
use crate::AeroacError;
use fs_material::fluid_acoustics::FluidAcousticState;
use fs_math::special::{SphericalBesselMoments, spherical_bessel_moments};
use fs_math::{c64::C64, dd::Dd, det};

/// Authority of the homogeneous, lossless compact-source power integral.
pub const COMPACT_POWER_SCOPE: &str = "Cycle-average full-sphere power of prescribed compact 3D volume flow [m^3/s] and force on fluid [N]; homogeneous stationary lossless free space, with no finite-source validation, moving sources, scattering, absorption, mutual loading or measured-SPL qualification.";

/// Admission bounds for geometry preparation and each coherent evaluation.
#[derive(Clone, Copy, Debug)]
pub struct RadiatedPowerBudget {
    /// Maximum original source count, including colocated or silent sources.
    pub max_sources: usize,
    /// Maximum square of the original source count, checked before allocation.
    /// Bounds grouping and the cached pair operator, including repeated group sums.
    pub max_source_pairs: usize,
}

/// Additional work admission for an explicit source-coherence factor bank.
#[derive(Clone, Copy, Debug)]
pub struct RadiatedPowerModeBudget {
    /// Maximum authored factors, including zero factors.
    pub max_modes: usize,
    /// Maximum original source-count squared times factor count.
    pub max_source_pair_mode_products: usize,
}

/// Outward power integrated over the complete sphere, with peak normalization.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CompactRadiatedPower {
    /// Time-average watts, including all declared interference.
    pub power_w: f64,
    /// The original prescribed-source and exterior-model limitations.
    pub scope: &'static str,
}

#[derive(Clone, Debug)]
struct PositionGroup {
    position_m: [f64; 3],
    members: core::ops::Range<usize>,
}

#[derive(Clone, Copy, Debug)]
struct PairMoment {
    // Direction from the earlier group to the later group.
    direction: [f64; 3],
    moments: SphericalBesselMoments,
}

/// Cached analytic angular integral for actual compact-source positions.
///
/// Exactly equal positions are grouped, preserving the original input order
/// and component phases. Evaluation first forms a field-summed colocated
/// baseline, then adds analytic separation corrections. The scalar and tensor
/// moment decrements remain separate from one and one-third, so a weak dipole
/// or monopole pair is not lost by subtracting two rounded O(1) powers.
///
/// For N supplied sources and M distinct positions, preparation is O(N^2)
/// and retains O(N + M^2) data. Immutable application allocates nothing, uses
/// no trigonometry, and takes O(N M + M^2) work: group sums are recomputed
/// when needed instead of requiring mutable or caller-owned scratch.
/// Double-double field and quadratic reductions retain small cancellation
/// residuals; this is not a certified enclosure or a relative guarantee near
/// every null. Numerical negativity refuses without a power clamp.
#[derive(Clone, Debug)]
pub struct PreparedCompactRadiatedPower {
    acoustic: FluidAcousticState,
    frequency_hz: f64,
    source_count: usize,
    source_pairs: usize,
    groups: Vec<PositionGroup>,
    members: Vec<usize>,
    pairs: Vec<PairMoment>,
    flow_gain: Dd,
    force_gain: Dd,
}

impl PreparedCompactRadiatedPower {
    /// Prepare one frequency and the source positions in one original frame.
    ///
    /// Q is physical volume injection [m^3/s], F is total force ON the fluid [N].
    /// Each source must admit the compact model; the entire array need not be
    /// compact. Source phases refer to one emission clock. No observer delay
    /// is applied to the input phasors before this source-domain integral.
    pub fn new(
        acoustic: FluidAcousticState,
        frequency_hz: f64,
        source_positions_m: &[[f64; 3]],
        budget: RadiatedPowerBudget,
    ) -> Result<Self, AeroacError> {
        Self::new_cancellable(acoustic, frequency_hz, source_positions_m, budget, || false)
    }

    /// Prepare with cancellation during source validation, grouping and pair work.
    /// No partially prepared operator is published after a refusal.
    pub fn new_cancellable(
        acoustic: FluidAcousticState,
        frequency_hz: f64,
        source_positions_m: &[[f64; 3]],
        budget: RadiatedPowerBudget,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, AeroacError> {
        Self::prepare(
            acoustic,
            frequency_hz,
            source_positions_m.len(),
            budget,
            |index| Ok(source_positions_m[index]),
            &mut cancelled,
        )
    }

    /// Preserve the exact source geometry, frequency and material of an observer.
    ///
    /// The observer position does not enter total radiated power. A molecular
    /// observer explicitly refuses: attenuation changes flux with radius, so
    /// neither lossless power nor power at infinity is its emitted source power.
    /// No temporary position array is allocated.
    pub fn from_observer(
        observer: &PreparedCompactSourceObserver,
        budget: RadiatedPowerBudget,
    ) -> Result<Self, AeroacError> {
        Self::from_observer_cancellable(observer, budget, || false)
    }

    /// Preserve an existing observer's source authority with cancellable preparation.
    pub fn from_observer_cancellable(
        observer: &PreparedCompactSourceObserver,
        budget: RadiatedPowerBudget,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, AeroacError> {
        poll(&mut cancelled)?;
        if observer.molecular_state().is_some() {
            return Err(invalid("full-sphere compact power requires a lossless exterior"));
        }
        Self::prepare(
            observer.acoustic(),
            observer.frequency_hz(),
            observer.source_count(),
            budget,
            |index| {
                observer.source_position_m(index)
                    .ok_or_else(|| invalid("compact power source index is unavailable"))
            },
            &mut cancelled,
        )
    }

    #[allow(clippy::float_cmp)] // Only exactly equal coordinates permit field coalescing.
    fn prepare(
        acoustic: FluidAcousticState,
        frequency_hz: f64,
        source_count: usize,
        budget: RadiatedPowerBudget,
        position: impl Fn(usize) -> Result<[f64; 3], AeroacError>,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, AeroacError> {
        poll(cancelled)?;
        let source_pairs = source_count.checked_mul(source_count)
            .ok_or_else(|| invalid("compact power source-pair count overflows"))?;
        if source_count > budget.max_sources || source_pairs > budget.max_source_pairs {
            return Err(invalid("compact power exceeds source/pair budget"));
        }
        let omega = core::f64::consts::TAU * frequency_hz;
        let k = omega / acoustic.sound_speed_m_s();
        if [frequency_hz, omega, k].iter().any(|value| !value.is_finite() || *value <= 0.0) {
            return Err(invalid("compact power frequency and wave number must be positive finite"));
        }
        let root_impedance = Dd::from_f64(acoustic.characteristic_impedance_pa_s_m()).sqrt();
        let base_gain = Dd::from_f64(k) / Dd::from_f64(8.0 * core::f64::consts::PI).sqrt();
        // These gains turn source coefficients into sqrt(watt) angular amplitudes.
        let flow_gain = base_gain * root_impedance;
        let force_gain = base_gain / root_impedance;
        if [flow_gain, force_gain].iter().any(|value| !finite_dd(*value) || value.hi <= 0.0) {
            return Err(invalid("compact power normalization is not representable"));
        }
        for index in 0..source_count {
            poll(cancelled)?;
            if position(index)?.iter().any(|value| !value.is_finite()) {
                return Err(nonfinite("compact power source position"));
            }
        }
        let mut groups = Vec::<PositionGroup>::new();
        let mut members = Vec::new();
        groups.try_reserve_exact(source_count)
            .map_err(|_| invalid("compact power group allocation exceeds available memory"))?;
        members.try_reserve_exact(source_count)
            .map_err(|_| invalid("compact power member allocation exceeds available memory"))?;
        for first in 0..source_count {
            poll(cancelled)?;
            let point = position(first)?;
            if groups.iter().any(|group| group.position_m == point) {
                continue;
            }
            let start = members.len();
            for index in first..source_count {
                poll(cancelled)?;
                if position(index)? == point {
                    members.push(index);
                }
            }
            groups.push(PositionGroup { position_m: point, members: start..members.len() });
        }
        // The admitted N^2 bounds this triangular M^2 count without overflow.
        let pair_count = groups.len() * groups.len().saturating_sub(1) / 2;
        let mut pairs = Vec::new();
        pairs.try_reserve_exact(pair_count)
            .map_err(|_| invalid("compact power pair allocation exceeds available memory"))?;
        for (i, a) in groups.iter().enumerate() {
            poll(cancelled)?;
            for b in &groups[..i] {
                poll(cancelled)?;
                let displacement: [f64; 3] = core::array::from_fn(|axis| a.position_m[axis] - b.position_m[axis]);
                let distance = det::hypot(det::hypot(displacement[0], displacement[1]), displacement[2]);
                let phase_distance = k * distance;
                if [distance, phase_distance].iter().any(|value| !value.is_finite() || *value <= 0.0) {
                    return Err(invalid("compact power source separation is not representable"));
                }
                let moments = spherical_bessel_moments(phase_distance)
                    .ok_or_else(|| invalid("compact power spherical moments are not representable"))?;
                pairs.push(PairMoment { direction: displacement.map(|value| value / distance), moments });
            }
        }
        poll(cancelled)?;
        Ok(Self { acoustic, frequency_hz, source_count, source_pairs, groups, members, pairs, flow_gain, force_gain })
    }

    /// Original fluid-acoustic snapshot, without a substituted medium.
    #[must_use]
    pub const fn acoustic_state(&self) -> FluidAcousticState {
        self.acoustic
    }

    /// One physical source frequency in Hz.
    #[must_use]
    pub const fn frequency_hz(&self) -> f64 {
        self.frequency_hz
    }

    /// Number of source vectors required in the original input order.
    #[must_use]
    pub const fn source_count(&self) -> usize {
        self.source_count
    }

    /// Number of exactly distinct source positions; close positions stay distinct.
    #[must_use]
    pub fn position_count(&self) -> usize {
        self.groups.len()
    }

    /// Total cycle-average power for coherent peak source phasors.
    pub fn evaluate(
        &self,
        sources: &[CompactSourcePhasor],
    ) -> Result<CompactRadiatedPower, AeroacError> {
        self.evaluate_cancellable(sources, || false)
    }

    /// Coherent power with source/pair cancellation polls and no partial result.
    pub fn evaluate_cancellable(
        &self,
        sources: &[CompactSourcePhasor],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<CompactRadiatedPower, AeroacError> {
        let power = self.coherent_power(sources, &mut cancelled)?;
        poll(&mut cancelled)?;
        finish(power)
    }

    /// Power of explicit source factors with E[xi_r conj(xi_s)] = delta_rs.
    ///
    /// Each factor supplies one Q/F vector per original source. Sources and
    /// mechanisms interfere within a factor; uncorrelated factors add powers.
    /// These are peak-unit spectral-line factors, not unscaled PSD/Hz. For a
    /// one-sided PSD bin of width df, factor the covariance 2 df S.
    /// This cycle-average quantity does not require a proper or Gaussian
    /// distribution, and does not infer correlation from unknown phases.
    pub fn evaluate_coherence_modes(
        &self,
        modes: &[&[CompactSourcePhasor]],
        budget: RadiatedPowerModeBudget,
    ) -> Result<CompactRadiatedPower, AeroacError> {
        self.evaluate_coherence_modes_cancellable(modes, budget, || false)
    }

    /// Bounded factor-bank evaluation with cancellation before publication.
    pub fn evaluate_coherence_modes_cancellable(
        &self,
        modes: &[&[CompactSourcePhasor]],
        budget: RadiatedPowerModeBudget,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<CompactRadiatedPower, AeroacError> {
        poll(&mut cancelled)?;
        let work = self.source_pairs.checked_mul(modes.len());
        if modes.len() > budget.max_modes
            || work.is_none_or(|count| count > budget.max_source_pair_mode_products)
        {
            return Err(invalid("compact power factors exceed the mode/pair budget"));
        }
        if modes.iter().any(|mode| mode.len() != self.source_count) {
            return Err(invalid("each compact power factor requires one vector per original source"));
        }
        let mut power = Dd::ZERO;
        for mode in modes {
            let contribution = self.coherent_power(mode, &mut cancelled)?;
            validate_power(contribution)?;
            power = power + contribution;
        }
        poll(&mut cancelled)?;
        finish(power)
    }

    /// O(N) specialization for caller-declared mutually uncorrelated sources.
    ///
    /// Colocated sources are still uncorrelated here: their fields must not be
    /// coalesced. Coherent Q/F components within one source are retained; their
    /// monopole/dipole cross term integrates to zero over a complete sphere.
    pub fn evaluate_independent_sources(
        &self,
        sources: &[CompactSourcePhasor],
    ) -> Result<CompactRadiatedPower, AeroacError> {
        self.evaluate_independent_sources_cancellable(sources, || false)
    }

    /// Independent-source power with cancellation and no partial publication.
    pub fn evaluate_independent_sources_cancellable(
        &self,
        sources: &[CompactSourcePhasor],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<CompactRadiatedPower, AeroacError> {
        poll(&mut cancelled)?;
        self.check_count(sources)?;
        let mut power = Dd::ZERO;
        for source in sources {
            poll(&mut cancelled)?;
            let mut sum = SourceSum::default();
            sum.add(*source)?;
            power = power + self.normalized(sum)?.colocated_power();
        }
        poll(&mut cancelled)?;
        finish(power)
    }

    fn coherent_power(
        &self,
        sources: &[CompactSourcePhasor],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Dd, AeroacError> {
        poll(cancelled)?;
        self.check_count(sources)?;
        let mut total = SourceSum::default();
        for source in sources {
            poll(cancelled)?;
            total.add(*source)?;
        }
        let mut power = self.normalized(total)?.colocated_power();
        let mut pair_index = 0;
        for (i, group) in self.groups.iter().enumerate().skip(1) {
            let a = self.group_sum(group, sources, cancelled)?;
            for other in &self.groups[..i] {
                poll(cancelled)?;
                let b = self.group_sum(other, sources, cancelled)?;
                let correction = self.pairs[pair_index].correction(a, b);
                power = power + correction + correction;
                pair_index += 1;
            }
        }
        validate_power(power)?;
        Ok(power)
    }

    fn group_sum(
        &self,
        group: &PositionGroup,
        sources: &[CompactSourcePhasor],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<RadiatingSource, AeroacError> {
        let mut sum = SourceSum::default();
        for &index in &self.members[group.members.clone()] {
            poll(cancelled)?;
            sum.add(sources[index])?;
        }
        self.normalized(sum)
    }

    fn normalized(&self, sum: SourceSum) -> Result<RadiatingSource, AeroacError> {
        let source = RadiatingSource {
            volume: sum.volume.scale(self.flow_gain),
            force: sum.force.map(|value| value.scale(self.force_gain)),
        };
        if !source.volume.finite() || source.force.iter().any(|value| !value.finite()) {
            return Err(nonfinite("coalesced compact radiating amplitudes"));
        }
        Ok(source)
    }

    fn check_count(&self, sources: &[CompactSourcePhasor]) -> Result<(), AeroacError> {
        if sources.len() == self.source_count {
            Ok(())
        } else {
            Err(invalid("compact power needs one vector per original source"))
        }
    }
}

impl PairMoment {
    fn correction(self, a: RadiatingSource, b: RadiatingSource) -> Dd {
        let an = a.project(self.direction);
        let bn = b.project(self.direction);
        let force_product = a.force.iter().zip(b.force)
            .fold(Dd::ZERO, |sum, (x, y)| sum + x.real_product(y));
        let mixed = a.volume.imaginary_product(bn) + an.imaginary_product(b.volume);
        self.moments.j0_minus_one * a.volume.real_product(b.volume)
            + self.moments.j1_over_x_minus_one_third * force_product
            - self.moments.j2 * an.real_product(bn)
            + self.moments.j1 * mixed
    }
}

#[derive(Clone, Copy, Default)]
struct ComplexDd {
    re: Dd,
    im: Dd,
}

impl ComplexDd {
    fn add(&mut self, value: C64) {
        self.re = self.re + Dd::from_f64(value.re);
        self.im = self.im + Dd::from_f64(value.im);
    }

    fn scale(self, value: Dd) -> Self {
        Self { re: self.re * value, im: self.im * value }
    }

    fn real_product(self, other: Self) -> Dd {
        self.re * other.re + self.im * other.im
    }

    fn imaginary_product(self, other: Self) -> Dd {
        self.im * other.re - self.re * other.im
    }

    fn finite(self) -> bool {
        finite_dd(self.re) && finite_dd(self.im)
    }
}

#[derive(Clone, Copy, Default)]
struct SourceSum {
    volume: ComplexDd,
    force: [ComplexDd; 3],
}

impl SourceSum {
    fn add(&mut self, source: CompactSourcePhasor) -> Result<(), AeroacError> {
        if !finite(source.volume_flow_m3_s) || source.force_on_fluid_n.iter().any(|value| !finite(*value)) {
            return Err(nonfinite("compact power source amplitude"));
        }
        self.volume.add(source.volume_flow_m3_s);
        for (sum, value) in self.force.iter_mut().zip(source.force_on_fluid_n) {
            sum.add(value);
        }
        Ok(())
    }
}

#[derive(Clone, Copy)]
struct RadiatingSource {
    volume: ComplexDd,
    force: [ComplexDd; 3],
}

impl RadiatingSource {
    fn colocated_power(self) -> Dd {
        self.volume.real_product(self.volume)
            + self.force.iter().fold(Dd::ZERO, |sum, value| sum + value.real_product(*value))
                / Dd::from_f64(3.0)
    }

    fn project(self, direction: [f64; 3]) -> ComplexDd {
        let mut result = ComplexDd::default();
        for (value, weight) in self.force.iter().zip(direction) {
            let scaled = value.scale(Dd::from_f64(weight));
            result.re = result.re + scaled.re;
            result.im = result.im + scaled.im;
        }
        result
    }
}

fn finite_dd(value: Dd) -> bool {
    value.hi.is_finite() && value.lo.is_finite()
}

fn validate_power(power: Dd) -> Result<(), AeroacError> {
    if !finite_dd(power) {
        Err(nonfinite("compact radiated power"))
    } else if power.lt(Dd::ZERO) {
        Err(invalid("compact radiated power lost numerical positivity"))
    } else {
        Ok(())
    }
}

fn finish(power: Dd) -> Result<CompactRadiatedPower, AeroacError> {
    validate_power(power)?;
    Ok(CompactRadiatedPower { power_w: power.to_f64(), scope: COMPACT_POWER_SCOPE })
}

fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), AeroacError> {
    if cancelled() {
        Err(AeroacError::Cancelled { what: "compact radiated power" })
    } else {
        Ok(())
    }
}
