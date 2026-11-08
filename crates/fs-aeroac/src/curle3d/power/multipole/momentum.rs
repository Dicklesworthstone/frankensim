//! Directional radiation of the original Q/F/Lighthill-stress source owner.
//!
//! Write the normalized amplitude as a + n.f + n.T.n with traceless T.
//! The original power operator has already absorbed trace(T) into a.
//! Its first angular moment needs the same moments through rank four and
//! <n_i n_j n_k n_l n_m exp(-i x n.e)> =
//! -i [a5 sym(delta delta e) + b5 sym(delta e e e) + j5 e_i e_j e_k e_l e_m],
//! a5 = j3/x^2, b5 = -j4/x. The sym sums have 15 and 10 distinct terms.
//! This is the derivative of the rank-four Legendre integral (DLMF 10.54.2).
//! All contractions below are for T and U traceless; the source trace itself
//! is retained in the original scalar amplitude, not discarded.

use super::super::momentum::{
    CompactRadiatedMomentum, ThirdMoment, colocated_moment, finish_momentum,
    first_moment_correction, poll,
};
use super::{
    AeroacError, CompactMultipolePhasor, ComplexDd, Dd, FluidAcousticState, MultipoleSum,
    PairMoment, PreparedMultipoleRadiatedPower, RadiatedPowerBudget, RadiatedPowerModeBudget,
    RadiatingMultipole, TensorPairMoment, det, finite_dd, imaginary_dot, invalid, nonfinite,
    project, real_dot, validate_power,
};

/// Authority of full-sphere directional radiation including resolved volume stress.
pub const MULTIPOLE_MOMENTUM_SCOPE: &str = "Cycle-average full-sphere outgoing wave momentum flux of prescribed compact 3D volume flow [m^3/s], force on fluid [N] and symmetric Lighthill stress [N m], including every coherent cross term in the original homogeneous stationary lossless exterior; not mechanical recoil, mutual loading, mean-flow momentum balance, wall-domain momentum, finite-source validation, unresolved turbulence or measured-SPL qualification.";

/// The original multipole power owner with directional emission moments.
///
/// Retains the original material, positions, trace decomposition and source
/// identity. Preparation adds three DD coefficients per distinct-position pair.
/// Coherent queries return power and its vector moment in one O(N M + M^2)
/// traversal for N sources at M positions, without allocation or trigonometry.
/// Independent-source queries take O(N); explicit coherence factors retain
/// all within-factor Q/F/T cross terms and add their signed vector moments.
///
/// Only free space is admitted. The scalar image-group symmetry shortcut does
/// not give the physical-domain vector integral or a wall momentum balance.
#[derive(Clone, Debug)]
pub struct PreparedMultipoleRadiatedMomentum {
    multipole: PreparedMultipoleRadiatedPower,
    moments: Vec<FifthMoment>,
}

impl PreparedMultipoleRadiatedPower {
    /// Extend this exact original free-space multipole owner with momentum flux.
    pub fn with_momentum(self) -> Result<PreparedMultipoleRadiatedMomentum, AeroacError> {
        self.with_momentum_cancellable(|| false)
    }

    /// Prepare rank-five coefficients under the original source/pair budget.
    pub fn with_momentum_cancellable(
        self,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<PreparedMultipoleRadiatedMomentum, AeroacError> {
        poll(&mut cancelled)?;
        if self.compact.image_group_size != 1 {
            return Err(invalid(
                "multipole momentum requires free space; image power symmetry is not vector symmetry",
            ));
        }
        let mut moments = Vec::new();
        moments
            .try_reserve_exact(self.moments.len())
            .map_err(|_| invalid("multipole momentum pair allocation exceeds available memory"))?;
        let k = core::f64::consts::TAU * self.compact.frequency_hz
            / self.compact.acoustic.sound_speed_m_s();
        for (i, a) in self.compact.groups.iter().enumerate().skip(1) {
            for b in &self.compact.groups[..i] {
                poll(&mut cancelled)?;
                let d: [f64; 3] =
                    core::array::from_fn(|axis| a.position_m[axis] - b.position_m[axis]);
                let x = k * det::hypot(det::hypot(d[0], d[1]), d[2]);
                moments.push(FifthMoment::new(x, self.moments[moments.len()])?);
            }
        }
        poll(&mut cancelled)?;
        Ok(PreparedMultipoleRadiatedMomentum {
            multipole: self,
            moments,
        })
    }
}

impl PreparedMultipoleRadiatedMomentum {
    /// Prepare physical Q/F/T sources with one frequency, frame and emission clock.
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
        PreparedMultipoleRadiatedPower::new_cancellable(
            acoustic,
            frequency_hz,
            positions_m,
            budget,
            &mut cancelled,
        )?
        .with_momentum_cancellable(cancelled)
    }

    /// Original power, material, source geometry and admitted model by reference.
    pub const fn power_operator(&self) -> &PreparedMultipoleRadiatedPower {
        &self.multipole
    }

    pub fn evaluate(
        &self,
        sources: &[CompactMultipolePhasor],
    ) -> Result<CompactRadiatedMomentum, AeroacError> {
        self.evaluate_cancellable(sources, || false)
    }

    /// Coherent application with no partial result on cancellation or refusal.
    pub fn evaluate_cancellable(
        &self,
        sources: &[CompactMultipolePhasor],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<CompactRadiatedMomentum, AeroacError> {
        let result = self.coherent(sources, &mut cancelled)?;
        poll(&mut cancelled)?;
        self.finish(result)
    }

    /// Peak-unit factors with E[xi_r conj(xi_s)] = delta_rs.
    ///
    /// Q, F and T remain coherent within each factor. No independence,
    /// pseudocovariance, Gaussian law or PSD-per-hertz units are inferred.
    pub fn evaluate_coherence_modes(
        &self,
        modes: &[&[CompactMultipolePhasor]],
        budget: RadiatedPowerModeBudget,
    ) -> Result<CompactRadiatedMomentum, AeroacError> {
        self.evaluate_coherence_modes_cancellable(modes, budget, || false)
    }

    pub fn evaluate_coherence_modes_cancellable(
        &self,
        modes: &[&[CompactMultipolePhasor]],
        budget: RadiatedPowerModeBudget,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<CompactRadiatedMomentum, AeroacError> {
        poll(&mut cancelled)?;
        if modes.len() > budget.max_modes
            || self
                .multipole
                .compact
                .source_pairs
                .checked_mul(modes.len())
                .is_none_or(|n| n > budget.max_source_pair_mode_products)
            || modes.iter().any(|mode| mode.len() != self.multipole.source_count())
        {
            return Err(invalid("multipole momentum factors exceed shape or work budget"));
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

    /// Caller-declared independent physical sources; each source's Q/F/T
    /// cross terms still direct radiation even at a shared position.
    pub fn evaluate_independent_sources(
        &self,
        sources: &[CompactMultipolePhasor],
    ) -> Result<CompactRadiatedMomentum, AeroacError> {
        self.evaluate_independent_sources_cancellable(sources, || false)
    }

    pub fn evaluate_independent_sources_cancellable(
        &self,
        sources: &[CompactMultipolePhasor],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<CompactRadiatedMomentum, AeroacError> {
        poll(&mut cancelled)?;
        self.multipole.check_count(sources)?;
        let mut total = (Dd::ZERO, [Dd::ZERO; 3]);
        for &source in sources {
            poll(&mut cancelled)?;
            let mut sum = MultipoleSum::default();
            sum.add(source)?;
            let source = self.multipole.normalized(sum)?;
            total.0 = total.0 + source.colocated_power();
            for (sum, value) in total.1.iter_mut().zip(colocated(source)) {
                *sum = *sum + value;
            }
        }
        poll(&mut cancelled)?;
        self.finish(total)
    }

    fn coherent(
        &self,
        sources: &[CompactMultipolePhasor],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(Dd, [Dd; 3]), AeroacError> {
        poll(cancelled)?;
        self.multipole.check_count(sources)?;
        let mut total = MultipoleSum::default();
        for &source in sources {
            poll(cancelled)?;
            total.add(source)?;
        }
        let source = self.multipole.normalized(total)?;
        let mut power = source.colocated_power();
        let mut moment = colocated(source);
        let mut pair_index = 0;
        for (i, group) in self.multipole.compact.groups.iter().enumerate().skip(1) {
            let a = self.multipole.group_sum(group, sources, cancelled)?;
            for other in &self.multipole.compact.groups[..i] {
                poll(cancelled)?;
                let b = self.multipole.group_sum(other, sources, cancelled)?;
                let pair = self.multipole.compact.pairs[pair_index];
                let tensor = self.multipole.moments[pair_index];
                let delta_power = tensor.correction(pair, a, b);
                power = power + delta_power + delta_power;
                let delta_moment = correction(pair, tensor, self.moments[pair_index], a, b);
                for (sum, value) in moment.iter_mut().zip(delta_moment) {
                    *sum = *sum + value + value;
                }
                pair_index += 1;
            }
        }
        validate_power(power)?;
        Ok((power, moment))
    }

    fn finish(&self, (power, moment): (Dd, [Dd; 3])) -> Result<CompactRadiatedMomentum, AeroacError> {
        finish_momentum(
            power,
            moment,
            self.multipole.compact.acoustic.sound_speed_m_s(),
            MULTIPOLE_MOMENTUM_SCOPE,
        )
    }
}

#[derive(Clone, Copy, Debug)]
struct FifthMoment {
    isotropic: Dd,
    directional: Dd,
    j5: Dd,
}

impl FifthMoment {
    fn new(x: f64, tensor: TensorPairMoment) -> Result<Self, AeroacError> {
        let (j5, isotropic, directional) = if x <= 4.0 {
            let xd = Dd::from_f64(x);
            let squared = xd * xd;
            let mut term = squared * squared * xd / Dd::from_f64(10395.0);
            let mut j5 = term;
            for step in 1..=24_u32 {
                let n = f64::from(step);
                term = -term * squared / Dd::from_f64(2.0 * n * (2.0 * n + 11.0));
                j5 = j5 + term;
            }
            // These recurrences for the QUOTIENTS avoid underflow of
            // j3 or j4 before division, even when j3/x^2 remains finite.
            let directional = -(tensor.j3 + j5) / Dd::from_f64(9.0);
            let isotropic =
                (tensor.mixed_third / Dd::from_f64(2.0) - directional) / Dd::from_f64(7.0);
            (j5, isotropic, directional)
        } else {
            let xd = Dd::from_f64(x);
            let directional = -tensor.j4 / xd;
            (
                -Dd::from_f64(9.0) * directional - tensor.j3,
                (tensor.j3 / xd) / xd,
                directional,
            )
        };
        if [j5, isotropic, directional].iter().any(|value| !finite_dd(*value)) {
            return Err(nonfinite("multipole momentum angular moments"));
        }
        Ok(Self { isotropic, directional, j5 })
    }
}

fn colocated(source: RadiatingMultipole) -> [Dd; 3] {
    let compact = colocated_moment(source.compact);
    let tensor_force = tensor_vector(source.tensor, source.compact.force, false);
    core::array::from_fn(|axis| {
        compact[axis] + Dd::from_f64(4.0) * tensor_force[axis] / Dd::from_f64(15.0)
    })
}

fn correction(
    pair: PairMoment,
    tensor: TensorPairMoment,
    fifth: FifthMoment,
    a: RadiatingMultipole,
    b: RadiatingMultipole,
) -> [Dd; 3] {
    let two = Dd::from_f64(2.0);
    let four = Dd::from_f64(4.0);
    let third = ThirdMoment { j3: tensor.j3, isotropic: tensor.mixed_third / two };
    let compact = first_moment_correction(pair, third, a.compact, b.compact);
    let at = a.tensor_direction(pair.direction);
    let bt = b.tensor_direction(pair.direction);
    let aq = project(at, pair.direction);
    let bq = project(bt, pair.direction);
    let af = a.compact.project(pair.direction);
    let bf = b.compact.project(pair.direction);
    let a_tensor_force = tensor_vector(a.tensor, b.compact.force, false);
    let b_tensor_force = tensor_vector(b.tensor, a.compact.force, false);
    let a_tensor_cross = tensor_vector(a.tensor, bt, true);
    let b_tensor_cross = tensor_vector(b.tensor, at, true);
    let fourth_directional = tensor.fourth_directional / four;

    // Rank three: scalar/stress. Rank four: force/stress, subtracting
    // its nonzero colocated limit. Rank five: stress/stress.
    let axial = -tensor.j3
        * (a.compact.volume.imaginary_product(bq) + aq.imaginary_product(b.compact.volume))
        + two * fourth_directional
            * (real_dot(a.compact.force, bt) + real_dot(at, b.compact.force))
        + tensor.j4 * (af.real_product(bq) + aq.real_product(bf))
        + two * fifth.isotropic * imaginary_tensor_product(a.tensor, b.tensor)
        + four * fifth.directional * imaginary_dot(at, bt)
        + fifth.j5 * aq.imaginary_product(bq);

    core::array::from_fn(|axis| {
        compact[axis] + Dd::from_f64(pair.direction[axis]) * axial
            + tensor.mixed_third
                * (a.compact.volume.imaginary_product(bt[axis])
                    + at[axis].imaginary_product(b.compact.volume))
            + tensor.fourth_isotropic_delta * (a_tensor_force[axis] + b_tensor_force[axis])
            + fourth_directional
                * (a.compact.force[axis].real_product(bq)
                    + aq.real_product(b.compact.force[axis])
                    + two * (af.real_product(bt[axis]) + at[axis].real_product(bf)))
            + four * fifth.isotropic * (a_tensor_cross[axis] - b_tensor_cross[axis])
            + two * fifth.directional
                * (at[axis].imaginary_product(bq) + aq.imaginary_product(bt[axis]))
    })
}

// Real or imaginary part of T conj(v), with the physical shear convention.
fn tensor_vector(t: [ComplexDd; 6], v: [ComplexDd; 3], imaginary: bool) -> [Dd; 3] {
    [[t[0], t[3], t[4]], [t[3], t[1], t[5]], [t[4], t[5], t[2]]].map(|row| {
        if imaginary { imaginary_dot(row, v) } else { real_dot(row, v) }
    })
}

fn imaginary_tensor_product(a: [ComplexDd; 6], b: [ComplexDd; 6]) -> Dd {
    a.into_iter().zip(b).enumerate().fold(Dd::ZERO, |sum, (index, (a, b))| {
        let term = a.imaginary_product(b);
        sum + if index < 3 { term } else { term + term }
    })
}
