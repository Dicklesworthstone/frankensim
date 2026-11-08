//! Radiated power including symmetric Lighthill stress moments.
//!
//! With peak phasors exp(-i omega t), the far pressure is
//! -i k exp(i k R)/(4 pi R) sum_j [rho c Q_j + n.F_j - i k n.T_j.n]
//! exp(-i k n.s_j). Thus the original compact power normalization also
//! applies to T = integral T_ij dV [N m], after multiplying it by -i k.
//! No density factor is applied to an already physical stress moment.
//!
//! The complete scalar through rank-four angular integral follows from
//! the Legendre identity DLMF 10.54.2 and the series DLMF 10.53.1:
//! <https://dlmf.nist.gov/10.54.E2>, <https://dlmf.nist.gov/10.53.E1>.
//! Tensor traces are absorbed into the scalar amplitude before squaring,
//! leaving a traceless tensor and preserving monopole/trace cancellation.
//! As in the original owner, pair corrections retain small remainders
//! separately from their colocated limits.

use super::{
    CompactRadiatedPower, CompactSourcePhasor, ComplexDd, Dd, PairMoment, PositionGroup,
    PreparedCompactRadiatedPower, RadiatedPowerBudget, RadiatedPowerModeBudget, RadiatingSource,
    SourceImage, SourceSum, finish, finite, finite_dd, invalid, nonfinite, poll, validate_power,
};
use super::momentum::ThirdMoment;
use crate::AeroacError;
use crate::curle3d::quadrupole::QuadrupolePhasor;
use fs_material::fluid_acoustics::FluidAcousticState;
use fs_math::{c64::C64, det, special::SphericalBesselMoments};

mod momentum;
pub use momentum::{MULTIPOLE_MOMENTUM_SCOPE, PreparedMultipoleRadiatedMomentum};

mod window;
pub use window::{
    MULTIPOLE_WINDOW_SCOPE, MultipoleRadiationLine, PreparedMultipoleRadiatedEnergy,
    RadiatedEnergyBudget, WindowedRadiatedEnergy,
};

/// Authority of the lossless complete compact-multipole power integral.
pub const MULTIPOLE_POWER_SCOPE: &str = "Cycle-average full-sphere power of prescribed compact 3D volume flow [m^3/s], force on fluid [N] and symmetric Lighthill stress moments [N m]; homogeneous stationary lossless free space with all coherent cross terms, but no source-size validation, unresolved turbulence, moving sources, scattering, absorption, mutual loading or measured-SPL qualification.";

/// Authority of a complete multipole image group in its physical fluid domain.
pub const MULTIPOLE_IDEAL_BOUNDARY_POWER_SCOPE: &str = "Cycle-average power of prescribed compact 3D volume flow [m^3/s], force on fluid [N] and symmetric Lighthill stress moments [N m] into the original homogeneous stationary lossless ideal-plane, perpendicular-corner or pi/m-wedge fluid domain; complete coherent image/tensor interference, but no finite wall impedance, general diffraction, absorption, unresolved turbulence, source feedback or measured-SPL qualification.";

/// The physical source terms at one original compact-source position.
///
/// All entries use the SAME emission clock and complex peak convention.
/// The stress is the original integrated symmetric tensor, including its
/// trace. Its xy/xz/yz entries are physical shear moments, not doubled.
/// The user owns source partitioning: do not supply the same physical term
/// both as an equivalent compact force/volume source and as volume stress.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct CompactMultipolePhasor {
    /// Physical volume injection and total force on the fluid.
    pub compact: CompactSourcePhasor,
    /// Resolved or otherwise explicitly supplied Lighthill stress moment.
    pub quadrupole: QuadrupolePhasor,
}

impl From<CompactSourcePhasor> for CompactMultipolePhasor {
    fn from(compact: CompactSourcePhasor) -> Self {
        Self {
            compact,
            quadrupole: QuadrupolePhasor::default(),
        }
    }
}

impl From<QuadrupolePhasor> for CompactMultipolePhasor {
    fn from(quadrupole: QuadrupolePhasor) -> Self {
        Self {
            compact: CompactSourcePhasor::default(),
            quadrupole,
        }
    }
}

/// Prepared analytic power of coherent monopoles, dipoles and quadrupoles.
///
/// The original geometry, material, exact-position groups, budgets and
/// optional ideal-wall images are retained by value. Only this extension
/// stores the additional rank-three/four moments; existing compact-only
/// preparation and evaluation remain unchanged. Preparation is O(P^2), with
/// O(P + M^2) retained data for P paths and M distinct positions. Evaluation
/// is allocation-free, contains no trigonometry and takes O(P M + M^2) work.
///
/// Every returned result is integrated far-field radiated power in watts,
/// not local p^2/(rho c), and not a sum of independently radiating images.
/// No partial owner/result is published on cancellation or numerical refusal.
#[derive(Clone, Debug)]
pub struct PreparedMultipoleRadiatedPower {
    compact: PreparedCompactRadiatedPower,
    stress_gain: Dd,
    moments: Vec<TensorPairMoment>,
}

impl PreparedCompactRadiatedPower {
    /// Extend this admitted geometry with full symmetric stress radiation.
    ///
    /// Image forces retain their original vector transform; a stress image
    /// uses sign * R T R^T. Source vectors still name only physical sources.
    pub fn with_quadrupoles(self) -> Result<PreparedMultipoleRadiatedPower, AeroacError> {
        self.with_quadrupoles_cancellable(|| false)
    }

    /// Extend with bounded rank-three/four preparation and cancellation polls.
    pub fn with_quadrupoles_cancellable(
        self,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<PreparedMultipoleRadiatedPower, AeroacError> {
        poll(&mut cancelled)?;
        let k = (core::f64::consts::TAU * self.frequency_hz) / self.acoustic.sound_speed_m_s();
        let stress_gain = self.force_gain * Dd::from_f64(k);
        if !finite_dd(stress_gain) || stress_gain.hi <= 0.0 {
            return Err(invalid(
                "multipole power normalization is not representable",
            ));
        }
        let mut moments = Vec::new();
        moments
            .try_reserve_exact(self.pairs.len())
            .map_err(|_| invalid("multipole power moment allocation exceeds available memory"))?;
        let mut pair = 0;
        for (index, a) in self.groups.iter().enumerate().skip(1) {
            for b in &self.groups[..index] {
                poll(&mut cancelled)?;
                let d: [f64; 3] =
                    core::array::from_fn(|axis| a.position_m[axis] - b.position_m[axis]);
                let x = k * det::hypot(det::hypot(d[0], d[1]), d[2]);
                moments.push(TensorPairMoment::new(x, self.pairs[pair].moments)?);
                pair += 1;
            }
        }
        poll(&mut cancelled)?;
        Ok(PreparedMultipoleRadiatedPower {
            compact: self,
            stress_gain,
            moments,
        })
    }
}

impl PreparedMultipoleRadiatedPower {
    /// Prepare one lossless frequency and the original physical source positions.
    pub fn new(
        acoustic: FluidAcousticState,
        frequency_hz: f64,
        source_positions_m: &[[f64; 3]],
        budget: RadiatedPowerBudget,
    ) -> Result<Self, AeroacError> {
        Self::new_cancellable(acoustic, frequency_hz, source_positions_m, budget, || false)
    }

    /// Prepare with the same source/pair budget and cancellation contract.
    pub fn new_cancellable(
        acoustic: FluidAcousticState,
        frequency_hz: f64,
        source_positions_m: &[[f64; 3]],
        budget: RadiatedPowerBudget,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, AeroacError> {
        PreparedCompactRadiatedPower::new_cancellable(
            acoustic,
            frequency_hz,
            source_positions_m,
            budget,
            &mut cancelled,
        )?
        .with_quadrupoles_cancellable(&mut cancelled)
    }

    /// Exact original geometry and low-order normalization, available without copying.
    #[must_use]
    pub const fn compact_power(&self) -> &PreparedCompactRadiatedPower {
        &self.compact
    }

    /// Original physical sources, excluding their images.
    #[must_use]
    pub const fn source_count(&self) -> usize {
        self.compact.source_count
    }

    /// Physical integration domain and multipole-model authority.
    #[must_use]
    pub const fn scope(&self) -> &'static str {
        if self.compact.image_group_size == 1 {
            MULTIPOLE_POWER_SCOPE
        } else {
            MULTIPOLE_IDEAL_BOUNDARY_POWER_SCOPE
        }
    }

    /// Full coherent power, including volume/force/stress cross terms.
    pub fn evaluate(
        &self,
        sources: &[CompactMultipolePhasor],
    ) -> Result<CompactRadiatedPower, AeroacError> {
        self.evaluate_cancellable(sources, || false)
    }

    /// Coherent power with cancellation before any result is published.
    pub fn evaluate_cancellable(
        &self,
        sources: &[CompactMultipolePhasor],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<CompactRadiatedPower, AeroacError> {
        let power = self.coherent_power(sources, &mut cancelled)?;
        poll(&mut cancelled)?;
        finish(power, self.scope())
    }

    /// Explicit source factors with E[xi_r conj(xi_s)] = delta_rs.
    ///
    /// Volume, force and stress components of each factor remain coherent,
    /// including across sources and all images of the same physical source.
    /// The result does not infer independent sources from unknown phases.
    pub fn evaluate_coherence_modes(
        &self,
        modes: &[&[CompactMultipolePhasor]],
        budget: RadiatedPowerModeBudget,
    ) -> Result<CompactRadiatedPower, AeroacError> {
        self.evaluate_coherence_modes_cancellable(modes, budget, || false)
    }

    /// Bounded factor evaluation with no partial result on cancellation.
    pub fn evaluate_coherence_modes_cancellable(
        &self,
        modes: &[&[CompactMultipolePhasor]],
        budget: RadiatedPowerModeBudget,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<CompactRadiatedPower, AeroacError> {
        poll(&mut cancelled)?;
        let work = self.compact.source_pairs.checked_mul(modes.len());
        if modes.len() > budget.max_modes
            || work.is_none_or(|count| count > budget.max_source_pair_mode_products)
        {
            return Err(invalid(
                "multipole power factors exceed the mode/pair budget",
            ));
        }
        if modes.iter().any(|mode| mode.len() != self.source_count()) {
            return Err(invalid(
                "each multipole factor requires one vector per original source",
            ));
        }
        let mut power = Dd::ZERO;
        for mode in modes {
            let contribution = self.coherent_power(mode, &mut cancelled)?;
            validate_power(contribution)?;
            power = power + contribution;
        }
        poll(&mut cancelled)?;
        finish(power, self.scope())
    }

    /// Power for caller-declared mutually uncorrelated physical sources.
    ///
    /// Q, F and T within one source retain all of their physical cross terms.
    /// Each source remains coherent with its own images. Free-space work is
    /// O(N); an image group of order g takes O(N g^2), without allocating
    /// sparse factors or mixing independent colocated sources.
    pub fn evaluate_independent_sources(
        &self,
        sources: &[CompactMultipolePhasor],
    ) -> Result<CompactRadiatedPower, AeroacError> {
        self.evaluate_independent_sources_cancellable(sources, || false)
    }

    /// Independent-source power with cancellation before publication.
    pub fn evaluate_independent_sources_cancellable(
        &self,
        sources: &[CompactMultipolePhasor],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<CompactRadiatedPower, AeroacError> {
        poll(&mut cancelled)?;
        self.check_count(sources)?;
        let mut power = Dd::ZERO;
        if self.compact.image_group_size == 1 {
            for &source in sources {
                poll(&mut cancelled)?;
                let mut sum = MultipoleSum::default();
                sum.add(source)?;
                power = power + self.normalized(sum)?.colocated_power();
            }
        } else {
            for (&source, images) in sources.iter().zip(
                self.compact
                    .images
                    .chunks_exact(self.compact.image_group_size),
            ) {
                poll(&mut cancelled)?;
                let mut total = MultipoleSum::default();
                for &image in images {
                    poll(&mut cancelled)?;
                    total.add_image(source, image)?;
                }
                let mut contribution = self.normalized(total)?.colocated_power();
                for (index, &a) in images.iter().enumerate().skip(1) {
                    let mut sum = MultipoleSum::default();
                    sum.add_image(source, a)?;
                    let left = self.normalized(sum)?;
                    for &b in &images[..index] {
                        poll(&mut cancelled)?;
                        if a.group == b.group {
                            continue;
                        }
                        let mut sum = MultipoleSum::default();
                        sum.add_image(source, b)?;
                        let right = self.normalized(sum)?;
                        let hi = a.group.max(b.group);
                        let lo = a.group.min(b.group);
                        let pair = hi * (hi - 1) / 2 + lo;
                        let correction = if a.group > b.group {
                            self.correction(pair, left, right)
                        } else {
                            self.correction(pair, right, left)
                        };
                        contribution = contribution + correction + correction;
                    }
                }
                validate_power(contribution)?;
                power = power + contribution;
            }
            power = power / Dd::from_f64(self.compact.image_group_size as f64);
        }
        poll(&mut cancelled)?;
        finish(power, self.scope())
    }

    fn coherent_power(
        &self,
        sources: &[CompactMultipolePhasor],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Dd, AeroacError> {
        poll(cancelled)?;
        self.check_count(sources)?;
        let mut total = MultipoleSum::default();
        if self.compact.image_group_size == 1 {
            for &source in sources {
                poll(cancelled)?;
                total.add(source)?;
            }
        } else {
            for &image in &self.compact.images {
                poll(cancelled)?;
                total.add_image(sources[image.source], image)?;
            }
        }
        let mut power = self.normalized(total)?.colocated_power();
        let mut pair = 0;
        for (index, group) in self.compact.groups.iter().enumerate().skip(1) {
            let a = self.group_sum(group, sources, cancelled)?;
            for other in &self.compact.groups[..index] {
                poll(cancelled)?;
                let b = self.group_sum(other, sources, cancelled)?;
                let correction = self.correction(pair, a, b);
                power = power + correction + correction;
                pair += 1;
            }
        }
        if self.compact.image_group_size != 1 {
            power = power / Dd::from_f64(self.compact.image_group_size as f64);
        }
        validate_power(power)?;
        Ok(power)
    }

    fn group_sum(
        &self,
        group: &PositionGroup,
        sources: &[CompactMultipolePhasor],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<RadiatingMultipole, AeroacError> {
        let mut sum = MultipoleSum::default();
        for &index in &self.compact.members[group.members.clone()] {
            poll(cancelled)?;
            if self.compact.image_group_size == 1 {
                sum.add(sources[index])?;
            } else {
                let image = self.compact.images[index];
                sum.add_image(sources[image.source], image)?;
            }
        }
        self.normalized(sum)
    }

    fn normalized(&self, sum: MultipoleSum) -> Result<RadiatingMultipole, AeroacError> {
        let mut compact = self.compact.normalized(sum.compact)?;
        // The -i factor is essential: stress and force have different orders
        // of time differentiation in their radiative pressures.
        let mut tensor = sum
            .stress
            .map(|value| value.scale(self.stress_gain).negative_i());
        let third = Dd::ONE / Dd::from_f64(3.0);
        let mean = tensor[0].plus(tensor[1]).plus(tensor[2]).scale(third);
        compact.volume = compact.volume.plus(mean);
        // Differences preserve exact isotropy; enforce zero trace in the same
        // summation order used by the analytic contraction.
        let xx = tensor[0]
            .minus(tensor[1])
            .plus(tensor[0].minus(tensor[2]))
            .scale(third);
        let yy = tensor[1]
            .minus(tensor[0])
            .plus(tensor[1].minus(tensor[2]))
            .scale(third);
        tensor[0] = xx;
        tensor[1] = yy;
        tensor[2] = xx.plus(yy).scale(-Dd::ONE);
        if !compact.volume.finite() || tensor.iter().any(|value| !value.finite()) {
            return Err(nonfinite("coalesced multipole radiating amplitudes"));
        }
        Ok(RadiatingMultipole { compact, tensor })
    }

    fn correction(&self, pair: usize, a: RadiatingMultipole, b: RadiatingMultipole) -> Dd {
        self.moments[pair].correction(self.compact.pairs[pair], a, b)
    }

    fn check_count(&self, sources: &[CompactMultipolePhasor]) -> Result<(), AeroacError> {
        if sources.len() == self.source_count() {
            Ok(())
        } else {
            Err(invalid(
                "multipole power needs one vector per original source",
            ))
        }
    }
}

#[derive(Clone, Copy, Debug)]
struct TensorPairMoment {
    j3: Dd,
    j4: Dd,
    mixed_third: Dd,
    fourth_isotropic_delta: Dd,
    fourth_directional: Dd,
}

impl TensorPairMoment {
    fn new(x: f64, lower: SphericalBesselMoments) -> Result<Self, AeroacError> {
        let third = ThirdMoment::new(x, lower)?;
        let j3 = third.j3;
        // Retain the rank-four series instead of a cancelling recurrence at
        // tiny separation. Rank three is shared with the original Q/F owner.
        let j4 = if x <= 4.0 {
            let squared = Dd::from_f64(x) * Dd::from_f64(x);
            let mut term = squared * squared / Dd::from_f64(945.0);
            let mut sum = term;
            for step in 1..=24_u32 {
                let n = f64::from(step);
                term = -term * squared / Dd::from_f64(2.0 * n * (2.0 * n + 9.0));
                sum = sum + term;
            }
            sum
        } else {
            Dd::from_f64(7.0) * j3 / Dd::from_f64(x) - lower.j2
        };
        let result = Self {
            j3,
            j4,
            mixed_third: Dd::from_f64(2.0) * third.isotropic,
            fourth_isotropic_delta: Dd::from_f64(2.0)
                * (lower.j0_minus_one / Dd::from_f64(15.0)
                    + Dd::from_f64(2.0) * lower.j2 / Dd::from_f64(21.0)
                    + j4 / Dd::from_f64(35.0)),
            fourth_directional: if x <= 4.0 {
                -Dd::from_f64(4.0) * (lower.j2 + j4) / Dd::from_f64(7.0)
            } else {
                // j4 approaches -j2: retain j3/x directly at large x.
                -Dd::from_f64(4.0) * j3 / Dd::from_f64(x)
            },
        };
        if [
            result.j3,
            result.j4,
            result.mixed_third,
            result.fourth_isotropic_delta,
            result.fourth_directional,
        ]
        .iter()
        .any(|value| !finite_dd(*value))
        {
            return Err(invalid("multipole spherical moments are not representable"));
        }
        Ok(result)
    }

    fn correction(self, pair: PairMoment, a: RadiatingMultipole, b: RadiatingMultipole) -> Dd {
        let at = a.tensor_direction(pair.direction);
        let bt = b.tensor_direction(pair.direction);
        let aq = project(at, pair.direction);
        let bq = project(bt, pair.direction);
        let af = a.compact.project(pair.direction);
        let bf = b.compact.project(pair.direction);
        let scalar_tensor = a.compact.volume.real_product(bq) + aq.real_product(b.compact.volume);
        let force_tensor = imaginary_dot(a.compact.force, bt) + imaginary_dot(at, b.compact.force);
        let axial_third = af.imaginary_product(bq) + aq.imaginary_product(bf);
        // The trace-free tensors eliminate every trace contraction. The
        // coefficients are those of delta-delta, e-e-delta and e-e-e-e.
        pair.correction(a.compact, b.compact) - pair.moments.j2 * scalar_tensor
            + self.mixed_third * force_tensor
            - self.j3 * axial_third
            + self.fourth_isotropic_delta * tensor_product(a.tensor, b.tensor)
            + self.fourth_directional * real_dot(at, bt)
            + self.j4 * aq.real_product(bq)
    }
}

#[derive(Clone, Copy, Default)]
struct MultipoleSum {
    compact: SourceSum,
    stress: [ComplexDd; 6],
}

const TENSOR_COMPONENTS: [(usize, usize); 6] = [(0, 0), (1, 1), (2, 2), (0, 1), (0, 2), (1, 2)];

impl MultipoleSum {
    fn add(&mut self, source: CompactMultipolePhasor) -> Result<(), AeroacError> {
        self.compact.add(source.compact)?;
        self.check_stress(source.quadrupole)?;
        for (sum, value) in self
            .stress
            .iter_mut()
            .zip(source.quadrupole.stress_moment_n_m)
        {
            sum.add(value);
        }
        Ok(())
    }

    fn add_image(
        &mut self,
        source: CompactMultipolePhasor,
        image: SourceImage,
    ) -> Result<(), AeroacError> {
        self.compact.add_image(source.compact, image)?;
        self.check_stress(source.quadrupole)?;
        let rotation = image.force_transform.map(|row| row.map(Dd::from_f64));
        for (target, &(row, column)) in TENSOR_COMPONENTS.iter().enumerate() {
            for (component, &(i, j)) in TENSOR_COMPONENTS.iter().enumerate() {
                let value = source.quadrupole.stress_moment_n_m[component];
                if value == C64::ZERO {
                    continue;
                }
                let mut coefficient = rotation[row][i] * rotation[column][j];
                if i != j {
                    coefficient = coefficient + rotation[row][j] * rotation[column][i];
                }
                coefficient = coefficient * Dd::from_f64(image.sign);
                self.stress[target].re =
                    self.stress[target].re + Dd::from_f64(value.re) * coefficient;
                self.stress[target].im =
                    self.stress[target].im + Dd::from_f64(value.im) * coefficient;
            }
        }
        Ok(())
    }

    fn check_stress(&self, source: QuadrupolePhasor) -> Result<(), AeroacError> {
        if source.stress_moment_n_m.iter().any(|value| !finite(*value)) {
            Err(nonfinite("multipole power stress amplitude"))
        } else {
            Ok(())
        }
    }
}

#[derive(Clone, Copy)]
struct RadiatingMultipole {
    compact: RadiatingSource,
    // Traceless xx, yy, zz, xy, xz, yz, in normalized sqrt(watt) units.
    tensor: [ComplexDd; 6],
}

impl RadiatingMultipole {
    fn colocated_power(self) -> Dd {
        self.compact.colocated_power()
            + Dd::from_f64(2.0) * tensor_product(self.tensor, self.tensor) / Dd::from_f64(15.0)
    }

    fn tensor_direction(self, direction: [f64; 3]) -> [ComplexDd; 3] {
        let t = self.tensor;
        [
            project([t[0], t[3], t[4]], direction),
            project([t[3], t[1], t[5]], direction),
            project([t[4], t[5], t[2]], direction),
        ]
    }
}

impl ComplexDd {
    fn plus(self, other: Self) -> Self {
        Self {
            re: self.re + other.re,
            im: self.im + other.im,
        }
    }
    fn minus(self, other: Self) -> Self {
        Self {
            re: self.re - other.re,
            im: self.im - other.im,
        }
    }
    fn negative_i(self) -> Self {
        Self {
            re: self.im,
            im: -self.re,
        }
    }
}

fn project(values: [ComplexDd; 3], direction: [f64; 3]) -> ComplexDd {
    values
        .into_iter()
        .zip(direction)
        .fold(ComplexDd::default(), |sum, (value, weight)| {
            sum.plus(value.scale(Dd::from_f64(weight)))
        })
}

fn real_dot(a: [ComplexDd; 3], b: [ComplexDd; 3]) -> Dd {
    a.into_iter()
        .zip(b)
        .fold(Dd::ZERO, |sum, (a, b)| sum + a.real_product(b))
}

fn imaginary_dot(a: [ComplexDd; 3], b: [ComplexDd; 3]) -> Dd {
    a.into_iter()
        .zip(b)
        .fold(Dd::ZERO, |sum, (a, b)| sum + a.imaginary_product(b))
}

fn tensor_product(a: [ComplexDd; 6], b: [ComplexDd; 6]) -> Dd {
    a.into_iter()
        .zip(b)
        .enumerate()
        .fold(Dd::ZERO, |sum, (index, (a, b))| {
            let term = a.real_product(b);
            sum + if index < 3 { term } else { term + term }
        })
}
