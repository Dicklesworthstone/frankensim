//! First angular moment of the ORIGINAL coherent outgoing power pattern.
//!
//! With the parent's sqrt(watt) amplitudes, A(n) = sum_j (a_j+n.f_j)
//! exp(-i k n.s_j), power is <|A|^2> and its first moment is <n |A|^2>.
//! Dividing that moment by the original c gives outgoing wave momentum flux.
//! It is not the mechanical force on a prescribed emitter: that requires the
//! source loading, mean flow and mechanical momentum balance.
//!
//! For e along a source separation and x=k*d, the extra angular identity is
//! <n_i n_j n_l exp(-i x n.e)> = -i [t*(delta_ij e_l + delta_il e_j
//! + delta_jl e_i) - j3(x)*e_i e_j e_l], t=(j1+j3)/5 = j2/x.
//! This follows from the Legendre integral and spherical-Bessel series in
//! <https://dlmf.nist.gov/10.54.E2> and <https://dlmf.nist.gov/10.53.E1>.
//! The colocated Q/F cross term is retained before separation corrections;
//! unlike scalar sphere power, it does not vanish in this vector integral.

use super::{
    CompactSourcePhasor, Dd, PairMoment, PreparedCompactRadiatedPower, RadiatedPowerBudget,
    RadiatedPowerModeBudget, RadiatingSource, SourceSum, finite_dd, invalid, nonfinite,
    validate_power,
};
use crate::AeroacError;
use fs_material::fluid_acoustics::FluidAcousticState;
use fs_math::det;

/// Authority of the original free-space power pattern's first angular moment.
pub const COMPACT_MOMENTUM_SCOPE: &str = "Cycle-average full-sphere outgoing wave momentum flux of prescribed compact 3D volume flow [m^3/s] and force on fluid [N], with all coherent cross terms in a homogeneous stationary lossless exterior; not a mechanical recoil force, mutual loading, mean-flow momentum balance, wall-domain result, finite-source validation or measured-SPL qualification.";

/// Cycle-average power and directional emission in the original Cartesian frame.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CompactRadiatedMomentum {
    /// Total outward power [W], including every declared interference term.
    pub power_w: f64,
    /// Integral of outward direction times differential radiated power [W].
    /// Dividing by nonzero `power_w` gives the power-weighted mean direction,
    /// whose norm is at most one. This is zero for an isotropic radiator.
    pub power_first_moment_w: [f64; 3],
    /// Outgoing wave momentum per time [N], `power_first_moment_w / c`.
    /// A signed radiation-pattern observable, not an applied mechanical force.
    pub wave_momentum_flux_n: [f64; 3],
    pub scope: &'static str,
}

#[derive(Clone, Copy, Debug)]
pub(super) struct ThirdMoment {
    pub(super) j3: Dd,
    pub(super) isotropic: Dd,
}

impl ThirdMoment {
    pub(super) fn new(
        x: f64,
        lower: fs_math::special::SphericalBesselMoments,
    ) -> Result<Self, AeroacError> {
        let (j3, isotropic) = if x <= 4.0 {
            // Upward recurrence would cancel tiny j3. The fixed
            // series has truncation below 1e-39 through x = 4.
            let xd = Dd::from_f64(x);
            let squared = xd * xd;
            let mut term = squared * xd / Dd::from_f64(105.0);
            let mut sum = term;
            for step in 1..=24_u32 {
                let n = f64::from(step);
                term = -term * squared / Dd::from_f64(2.0 * n * (2.0 * n + 7.0));
                sum = sum + term;
            }
            (sum, (lower.j1 + sum) / Dd::from_f64(5.0))
        } else {
            let isotropic = lower.j2 / Dd::from_f64(x);
            // At large separation j3 approximately equals -j1.
            // Keep j2/x directly, rather than cancelling j1+j3.
            (Dd::from_f64(5.0) * isotropic - lower.j1, isotropic)
        };
        if !finite_dd(j3) || !finite_dd(isotropic) {
            return Err(nonfinite("compact momentum angular moments"));
        }
        Ok(Self { j3, isotropic })
    }
}

/// Original compact-array power plus cached rank-three angular moments.
///
/// Reuses the source grouping, material, phases and pair budget of the power
/// owner. Preparation adds two DD coefficients per distinct-position pair.
/// Coherent evaluation computes power and its vector moment in one traversal,
/// with no allocation or transcendental calls. Its work is O(N M + M^2) for
/// N sources at M positions. Independent-source evaluation is O(N).
///
/// Only the complete free-space sphere is admitted. Full-image scalar power
/// can be divided by image-group order, but the corresponding VECTOR integral
/// cannot: a reflecting wall can receive normal momentum. No such symmetry
/// shortcut or a wall-force claim is inferred here.
#[derive(Clone, Debug)]
pub struct PreparedCompactRadiatedMomentum {
    compact: PreparedCompactRadiatedPower,
    moments: Vec<ThirdMoment>,
}

impl PreparedCompactRadiatedPower {
    /// Extend this original free-space source owner with directional emission.
    pub fn with_momentum(self) -> Result<PreparedCompactRadiatedMomentum, AeroacError> {
        self.with_momentum_cancellable(|| false)
    }

    /// Cache the additional angular moments with bounded pair cancellation.
    pub fn with_momentum_cancellable(
        self,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<PreparedCompactRadiatedMomentum, AeroacError> {
        poll(&mut cancelled)?;
        if self.image_group_size != 1 {
            return Err(invalid(
                "compact momentum requires free space; image power symmetry is not vector symmetry",
            ));
        }
        let mut moments = Vec::new();
        moments
            .try_reserve_exact(self.pairs.len())
            .map_err(|_| invalid("compact momentum pair allocation exceeds available memory"))?;
        let k = core::f64::consts::TAU * self.frequency_hz / self.acoustic.sound_speed_m_s();
        for (i, a) in self.groups.iter().enumerate().skip(1) {
            for b in &self.groups[..i] {
                poll(&mut cancelled)?;
                let d: [f64; 3] =
                    core::array::from_fn(|axis| a.position_m[axis] - b.position_m[axis]);
                let x = k * det::hypot(det::hypot(d[0], d[1]), d[2]);
                let pair = self.pairs[moments.len()];
                moments.push(ThirdMoment::new(x, pair.moments)?);
            }
        }
        poll(&mut cancelled)?;
        Ok(PreparedCompactRadiatedMomentum {
            compact: self,
            moments,
        })
    }
}

impl PreparedCompactRadiatedMomentum {
    /// Prepare physical Q/F sources sharing one frequency, frame and source clock.
    pub fn new(
        acoustic: FluidAcousticState,
        frequency_hz: f64,
        positions_m: &[[f64; 3]],
        budget: RadiatedPowerBudget,
    ) -> Result<Self, AeroacError> {
        Self::new_cancellable(acoustic, frequency_hz, positions_m, budget, || false)
    }

    pub fn new_cancellable(
        acoustic: FluidAcousticState,
        frequency_hz: f64,
        positions_m: &[[f64; 3]],
        budget: RadiatedPowerBudget,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, AeroacError> {
        PreparedCompactRadiatedPower::new_cancellable(
            acoustic,
            frequency_hz,
            positions_m,
            budget,
            &mut cancelled,
        )?
        .with_momentum_cancellable(cancelled)
    }

    /// Original power owner, including the retained material and source identity.
    pub const fn power_operator(&self) -> &PreparedCompactRadiatedPower {
        &self.compact
    }

    pub fn evaluate(
        &self,
        sources: &[CompactSourcePhasor],
    ) -> Result<CompactRadiatedMomentum, AeroacError> {
        self.evaluate_cancellable(sources, || false)
    }

    /// Sum coherent fields before quadratic angular integration; no partial result.
    pub fn evaluate_cancellable(
        &self,
        sources: &[CompactSourcePhasor],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<CompactRadiatedMomentum, AeroacError> {
        let result = self.coherent(sources, &mut cancelled)?;
        poll(&mut cancelled)?;
        self.finish(result)
    }

    /// Explicit peak-unit source factors with E[xi_r conj(xi_s)] = delta_rs.
    /// Each factor interferes coherently; different factors add both power and
    /// vector momentum. No independence, Gaussian law or PSD/Hz units are inferred.
    pub fn evaluate_coherence_modes(
        &self,
        modes: &[&[CompactSourcePhasor]],
        budget: RadiatedPowerModeBudget,
    ) -> Result<CompactRadiatedMomentum, AeroacError> {
        self.evaluate_coherence_modes_cancellable(modes, budget, || false)
    }

    pub fn evaluate_coherence_modes_cancellable(
        &self,
        modes: &[&[CompactSourcePhasor]],
        budget: RadiatedPowerModeBudget,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<CompactRadiatedMomentum, AeroacError> {
        poll(&mut cancelled)?;
        if modes.len() > budget.max_modes
            || self
                .compact
                .source_pairs
                .checked_mul(modes.len())
                .is_none_or(|n| n > budget.max_source_pair_mode_products)
            || modes.iter().any(|mode| mode.len() != self.compact.source_count)
        {
            return Err(invalid("compact momentum factors exceed shape or work budget"));
        }
        let mut total = (Dd::ZERO, [Dd::ZERO; 3]);
        for mode in modes {
            let (power, moment) = self.coherent(mode, &mut cancelled)?;
            total.0 = total.0 + power;
            for (sum, value) in total.1.iter_mut().zip(moment) {
                *sum = *sum + value;
            }
        }
        poll(&mut cancelled)?;
        self.finish(total)
    }

    /// Caller-declared independent physical sources; each source's Q/F cross
    /// term still contributes directional momentum even at a shared position.
    pub fn evaluate_independent_sources(
        &self,
        sources: &[CompactSourcePhasor],
    ) -> Result<CompactRadiatedMomentum, AeroacError> {
        self.evaluate_independent_sources_cancellable(sources, || false)
    }

    pub fn evaluate_independent_sources_cancellable(
        &self,
        sources: &[CompactSourcePhasor],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<CompactRadiatedMomentum, AeroacError> {
        poll(&mut cancelled)?;
        self.compact.check_count(sources)?;
        let mut total = (Dd::ZERO, [Dd::ZERO; 3]);
        for source in sources {
            poll(&mut cancelled)?;
            let mut sum = SourceSum::default();
            sum.add(*source)?;
            let source = self.compact.normalized(sum)?;
            total.0 = total.0 + source.colocated_power();
            for (sum, value) in total.1.iter_mut().zip(colocated_moment(source)) {
                *sum = *sum + value;
            }
        }
        poll(&mut cancelled)?;
        self.finish(total)
    }

    fn coherent(
        &self,
        sources: &[CompactSourcePhasor],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(Dd, [Dd; 3]), AeroacError> {
        poll(cancelled)?;
        self.compact.check_count(sources)?;
        let mut total = SourceSum::default();
        for source in sources {
            poll(cancelled)?;
            total.add(*source)?;
        }
        let source = self.compact.normalized(total)?;
        let mut power = source.colocated_power();
        let mut moment = colocated_moment(source);
        let mut index = 0;
        for (i, group) in self.compact.groups.iter().enumerate().skip(1) {
            let a = self.compact.group_sum(group, sources, cancelled)?;
            for other in &self.compact.groups[..i] {
                poll(cancelled)?;
                let b = self.compact.group_sum(other, sources, cancelled)?;
                let pair = self.compact.pairs[index];
                let correction = pair.correction(a, b);
                power = power + correction + correction;
                let correction = first_moment_correction(pair, self.moments[index], a, b);
                for (sum, value) in moment.iter_mut().zip(correction) {
                    *sum = *sum + value + value;
                }
                index += 1;
            }
        }
        validate_power(power)?;
        Ok((power, moment))
    }

    fn finish(&self, (power, moment): (Dd, [Dd; 3])) -> Result<CompactRadiatedMomentum, AeroacError> {
        finish_momentum(
            power,
            moment,
            self.compact.acoustic.sound_speed_m_s(),
            COMPACT_MOMENTUM_SCOPE,
        )
    }
}

pub(super) fn colocated_moment(source: RadiatingSource) -> [Dd; 3] {
    let two_thirds = Dd::from_f64(2.0) / Dd::from_f64(3.0);
    source.force.map(|force| two_thirds * source.volume.real_product(force))
}

pub(super) fn first_moment_correction(
    pair: PairMoment,
    third: ThirdMoment,
    a: RadiatingSource,
    b: RadiatingSource,
) -> [Dd; 3] {
    let an = a.project(pair.direction);
    let bn = b.project(pair.direction);
    let force_dot = a.force.iter().zip(b.force).fold(Dd::ZERO, |sum, (a, b)| {
        sum + a.imaginary_product(b)
    });
    let axial = pair.moments.j1 * a.volume.imaginary_product(b.volume)
        - pair.moments.j2 * (a.volume.real_product(bn) + an.real_product(b.volume))
        + third.isotropic * force_dot
        - third.j3 * an.imaginary_product(bn);
    core::array::from_fn(|axis| {
        Dd::from_f64(pair.direction[axis]) * axial
            + pair.moments.j1_over_x_minus_one_third
                * (a.volume.real_product(b.force[axis]) + a.force[axis].real_product(b.volume))
            + third.isotropic
                * (a.force[axis].imaginary_product(bn) + an.imaginary_product(b.force[axis]))
    })
}

pub(super) fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), AeroacError> {
    if cancelled() {
        Err(AeroacError::Cancelled {
            what: "compact radiated momentum",
        })
    } else {
        Ok(())
    }
}

pub(super) fn finish_momentum(
    power: Dd,
    moment: [Dd; 3],
    sound_speed_m_s: f64,
    scope: &'static str,
) -> Result<CompactRadiatedMomentum, AeroacError> {
    validate_power(power)?;
    let c = Dd::from_f64(sound_speed_m_s);
    let flux = moment.map(|value| value / c);
    if moment.iter().chain(&flux).any(|value| !finite_dd(*value)) {
        return Err(nonfinite("compact radiated momentum"));
    }
    let result = CompactRadiatedMomentum {
        power_w: power.to_f64(),
        power_first_moment_w: moment.map(Dd::to_f64),
        wave_momentum_flux_n: flux.map(Dd::to_f64),
        scope,
    };
    if !result.power_w.is_finite()
        || result
            .power_first_moment_w
            .iter()
            .chain(&result.wave_momentum_flux_n)
            .any(|value| !value.is_finite())
    {
        return Err(nonfinite("compact radiated momentum"));
    }
    Ok(result)
}
