//! Total perturbation-energy transport of fixed monopoles in uniform flow.
//!
//! This integrates <p u> + U <E>, not pressure work alone. With M=U/c,
//! beta²=1-|M|² and e parallel to M, put A=beta I+(beta²-beta)e e^T.
//! On an enclosing ellipsoid x=A r n the outgoing Green function has phase
//! k(r-M.r n). A source at s contributes Q exp[i k M.A^-1 s]
//! exp[-i k n.A^-1 s]. The flux through that ellipsoid is therefore
//!
//! P=rho omega²/(8 pi c beta^4) <(1-M.n)|sum_j Q_j
//!       exp[-i k(n-M).A^-1 s_j]|²>_sphere.
//!
//! For d=A^-1(s_j-s_l), x=k|d|, a=M.d/|d|, the exact pair kernel is
//! exp(i a x)[j0(x)+i a j1(x)]. This follows by differentiating the scalar
//! sphere integral <exp(-i n.x)>=j0(|x|), DLMF 10.54.2:
//! <https://dlmf.nist.gov/10.54.E2>. The convected Green function and source
//! convention are those of [`crate::curle3d::array::convected`], retaining
//! Wang, *Flow effects in computational aeroacoustics* (2020), equations
//! 3.28–3.30: <https://pure.tue.nl/ws/portalfiles/portal/149260007/20200511_Z.Wang.pdf>.
//! No extra source force, entropy, or mean-flow energy extraction is inferred.
//!
//! For the original observer's perpendicular ideal planes or pi/m wedge,
//! U is tangent to every wall. Reflection R with parity s gives p(Rx)=s p(x),
//! u(Rx)=s R u(x), E(Rx)=E(x), and J(Rx)=R J(x). The energy flux through each
//! reflected sector of an enclosing sphere is equal, and wall flux vanishes.
//! Whole-image power divided by the group order (2^planes or 2m for a wedge)
//! is thus power into the physical domain.

use super::{
    CompactRadiatedPower, ComplexDd, PositionGroup, RadiatedPowerBudget, RadiatedPowerModeBudget,
    finish, finite_dd, poll,
};
use crate::AeroacError;
use crate::curle3d::array::convected::PreparedConvectedMonopoleObserver;
use crate::curle3d::array::wedge::WedgeImages;
use crate::curle3d::statistics::{finite, invalid, nonfinite};
use fs_material::fluid_acoustics::FluidAcousticState;
use fs_math::{c64::C64, dd::Dd, det, special::spherical_bessel_moments};

/// Physical authority of the integrated quantity and the prescribed source law.
pub const CONVECTED_MONOPOLE_POWER_SCOPE: &str = "Cycle-average total outward quadratic acoustic perturbation-energy transport <p u> + U <E> of fixed prescribed compact volume sources injecting at the uniform mean velocity; homogeneous subsonic lossless free space, with coherent interference; no extra source force/entropy, nonuniform flow, absorption, scattering, source feedback, finite-source validation or measured-SPL qualification.";

/// Authority of power into an existing observer's ideal reflecting domain.
pub const CONVECTED_MONOPOLE_IDEAL_BOUNDARY_POWER_SCOPE: &str = "Cycle-average total outward quadratic acoustic perturbation-energy transport <p u> + U <E> of fixed prescribed compact volume sources injecting at the uniform mean velocity into the original homogeneous subsonic lossless domain of one to three perpendicular infinite ideal rigid or pressure-release planes, or an admitted ideal pi/m wedge; mean flow is tangent to every plane and each physical source retains the complete coherent image group. No extra source force/entropy, nonuniform flow, wall layers, finite impedance, general diffraction, absorption, source feedback, finite-source validation or measured-SPL qualification.";

#[derive(Clone, Copy, Debug)]
struct ConvectedPair {
    real_minus_one: Dd,
    imaginary: Dd,
}

#[derive(Clone, Copy, Debug)]
struct ConvectedImage {
    source: usize,
    position_m: [f64; 3],
    sign: f64,
}

/// Analytic radiated power of coherent fixed volume sources in uniform flow.
///
/// For P direct/image paths, preparation caches the exact scalar/vector
/// spherical moments in O(P²) work. Evaluation is allocation-free O(P²),
/// without angular quadrature or trigonometry. Exactly coincident positions
/// share a field sum; all source
/// phases remain in the original common clock. A cancellation-preserving
/// coincident baseline plus pair decrements retains weak oppositely driven
/// sources even when their separation is far below a wavelength.
///
/// The power is independent of the chosen enclosing surface: away from the
/// sources the cycle-average energy equation has div(<p u>+U<E>)=0.
#[derive(Clone, Debug)]
pub struct PreparedConvectedMonopolePower {
    acoustic: FluidAcousticState,
    flow_velocity_m_s: [f64; 3],
    frequency_hz: f64,
    source_count: usize,
    source_pairs: usize,
    groups: Vec<PositionGroup>,
    members: Vec<usize>,
    pairs: Vec<ConvectedPair>,
    // Physical volume flow times this gain has units sqrt(watt).
    gain: Dd,
    // Empty in free space; otherwise every path maps to its original source.
    images: Vec<ConvectedImage>,
    image_group_size: usize,
}

impl PreparedConvectedMonopolePower {
    /// Prepare physical source positions and an explicit inclusive Mach margin.
    /// `maximum_mach` must be in [0,1); zero admits only a stationary medium.
    /// Source vectors contain peak volume-flow phasors [m³/s], in exp(-i omega t).
    pub fn new(
        acoustic: FluidAcousticState,
        flow_velocity_m_s: [f64; 3],
        frequency_hz: f64,
        source_positions_m: &[[f64; 3]],
        maximum_mach: f64,
        budget: RadiatedPowerBudget,
    ) -> Result<Self, AeroacError> {
        Self::new_cancellable(
            acoustic,
            flow_velocity_m_s,
            frequency_hz,
            source_positions_m,
            maximum_mach,
            budget,
            || false,
        )
    }

    /// Cancellable preparation; refusal never exposes a partially prepared bank.
    #[allow(clippy::too_many_arguments)]
    pub fn new_cancellable(
        acoustic: FluidAcousticState,
        flow_velocity_m_s: [f64; 3],
        frequency_hz: f64,
        source_positions_m: &[[f64; 3]],
        maximum_mach: f64,
        budget: RadiatedPowerBudget,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, AeroacError> {
        Self::prepare(
            acoustic,
            flow_velocity_m_s,
            frequency_hz,
            source_positions_m.len(),
            |index| Ok(source_positions_m[index]),
            maximum_mach,
            budget,
            &mut cancelled,
        )
    }

    /// Preserve an observer's material, actual flow, frequency and source
    /// positions. Its receiver position never enters total power. Original
    /// perpendicular ideal planes or pi/m wedges retain their complete coherent
    /// image group; power is integrated only over the physical fluid sector. Geometry,
    /// strict source half-spaces and tangent flow were admitted by that owner.
    /// Budgets count EVERY direct/image path and its square, while later
    /// amplitude vectors still contain only the original physical sources.
    pub fn from_observer(
        observer: &PreparedConvectedMonopoleObserver,
        maximum_mach: f64,
        budget: RadiatedPowerBudget,
    ) -> Result<Self, AeroacError> {
        Self::from_observer_cancellable(observer, maximum_mach, budget, || false)
    }

    /// Poll during image expansion and the original bounded pair preparation.
    /// Cancellation or refusal never publishes an incomplete image group.
    pub fn from_observer_cancellable(
        observer: &PreparedConvectedMonopoleObserver,
        maximum_mach: f64,
        budget: RadiatedPowerBudget,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, AeroacError> {
        poll(&mut cancelled)?;
        let planes = observer.reflection_planes();
        if planes.is_empty() {
            return Self::prepare(
                observer.acoustic(),
                observer.flow_velocity_m_s(),
                observer.frequency_hz(),
                observer.source_count(),
                |index| {
                    observer
                        .source_position_m(index)
                        .ok_or_else(|| invalid("convected power source index is unavailable"))
                },
                maximum_mach,
                budget,
                &mut cancelled,
            );
        }
        // Retain the original observer's geometry owner. A two-wall pi/m
        // wedge has 2m paths, not the four paths of perpendicular planes.
        // The shared visitor also retains translated affine reflections and
        // the even-order restriction for mixed rigid/release wedge walls.
        let wedge = match observer.reflection_wedge_order() {
            Some(order) => Some(WedgeImages::new(
                planes.try_into().map_err(|_| {
                    invalid("convected wedge power requires the original two planes")
                })?,
                order,
            )?),
            None => None,
        };
        let image_group_size = wedge.map_or(1usize << planes.len(), WedgeImages::path_count);
        let paths = observer
            .source_count()
            .checked_mul(image_group_size)
            .ok_or_else(|| invalid("convected power image path count overflows"))?;
        let pairs = paths
            .checked_mul(paths)
            .ok_or_else(|| invalid("convected power image pair count overflows"))?;
        if paths > budget.max_sources || pairs > budget.max_source_pairs {
            return Err(invalid("convected images exceed source/pair budget"));
        }
        let mut images = Vec::new();
        images
            .try_reserve_exact(paths)
            .map_err(|_| invalid("convected power image allocation exceeds available memory"))?;
        for source in 0..observer.source_count() {
            poll(&mut cancelled)?;
            let position = observer
                .source_position_m(source)
                .ok_or_else(|| invalid("convected power source index is unavailable"))?;
            if let Some(wedge) = wedge {
                // As in the convected field owner, H^-1 enumerates the entire
                // source group and chi(H^-1)=chi(H). Preserve the direct point
                // exactly and use the SAME affine map for every other image.
                wedge.visit_inverse(|index, image, sign| {
                    poll(&mut cancelled)?;
                    images.push(ConvectedImage {
                        source,
                        position_m: if index == 0 {
                            position
                        } else {
                            image.position(position)?
                        },
                        sign,
                    });
                    Ok(())
                })?;
            } else {
                for mask in 0..image_group_size {
                    poll(&mut cancelled)?;
                    let mut image = ConvectedImage {
                        source,
                        position_m: position,
                        sign: 1.0,
                    };
                    for (index, plane) in planes.iter().enumerate() {
                        if mask & (1 << index) != 0 {
                            // Same DD reflection as the admitted observer and the
                            // stationary power owner; no new geometry is inferred.
                            image.position_m = plane.reflected_position(image.position_m)?;
                            image.sign *= plane.boundary().sign();
                        }
                    }
                    images.push(image);
                }
            }
        }
        let mut prepared = Self::prepare(
            observer.acoustic(),
            observer.flow_velocity_m_s(),
            observer.frequency_hz(),
            images.len(),
            |index| Ok(images[index].position_m),
            maximum_mach,
            budget,
            &mut cancelled,
        )?;
        prepared.source_count = observer.source_count();
        prepared.images = images;
        prepared.image_group_size = image_group_size;
        poll(&mut cancelled)?;
        Ok(prepared)
    }

    #[allow(clippy::too_many_arguments)]
    fn prepare(
        acoustic: FluidAcousticState,
        flow_velocity_m_s: [f64; 3],
        frequency_hz: f64,
        source_count: usize,
        position: impl Fn(usize) -> Result<[f64; 3], AeroacError>,
        maximum_mach: f64,
        budget: RadiatedPowerBudget,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, AeroacError> {
        poll(cancelled)?;
        let source_pairs = source_count
            .checked_mul(source_count)
            .ok_or_else(|| invalid("convected power source-pair count overflows"))?;
        if source_count > budget.max_sources || source_pairs > budget.max_source_pairs {
            return Err(invalid("convected power exceeds source/pair budget"));
        }
        let c = acoustic.sound_speed_m_s();
        let omega = core::f64::consts::TAU * frequency_hz;
        let k = omega / c;
        let mach = flow_velocity_m_s.map(|v| v / c);
        let speed = norm(mach);
        if !maximum_mach.is_finite()
            || !(0.0..1.0).contains(&maximum_mach)
            || !speed.is_finite()
            || speed > maximum_mach
            || [frequency_hz, omega, k]
                .iter()
                .any(|v| !v.is_finite() || *v <= 0.0)
        {
            return Err(invalid(
                "convected power frequency and subsonic flow margin",
            ));
        }
        let beta_squared = (1.0 - speed) * (1.0 + speed);
        let beta = det::sqrt(beta_squared);
        let direction = if speed == 0.0 {
            [0.0; 3]
        } else {
            mach.map(|v| v / speed)
        };
        let gain = Dd::from_f64(k)
            * (Dd::from_f64(acoustic.characteristic_impedance_pa_s_m())
                / Dd::from_f64(8.0 * core::f64::consts::PI))
            .sqrt()
            / Dd::from_f64(beta_squared);
        if !finite_dd(gain) || gain.hi <= 0.0 {
            return Err(invalid(
                "convected power normalization is not representable",
            ));
        }
        for index in 0..source_count {
            poll(cancelled)?;
            if position(index)?.iter().any(|v| !v.is_finite()) {
                return Err(nonfinite("convected power source position"));
            }
        }
        let mut groups = Vec::<PositionGroup>::new();
        let mut members = Vec::new();
        groups
            .try_reserve_exact(source_count)
            .map_err(|_| invalid("convected power group allocation exceeds available memory"))?;
        members
            .try_reserve_exact(source_count)
            .map_err(|_| invalid("convected power member allocation exceeds available memory"))?;
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
            groups.push(PositionGroup {
                position_m: point,
                members: start..members.len(),
            });
        }
        let pair_count = groups.len() * groups.len().saturating_sub(1) / 2;
        let mut pairs = Vec::new();
        pairs
            .try_reserve_exact(pair_count)
            .map_err(|_| invalid("convected power pair allocation exceeds available memory"))?;
        for (index, a) in groups.iter().enumerate() {
            for b in &groups[..index] {
                poll(cancelled)?;
                let separation: [f64; 3] =
                    core::array::from_fn(|i| a.position_m[i] - b.position_m[i]);
                // A^-1 d = d/beta + e(e.d)(1/beta²-1/beta).
                let parallel = dot(direction, separation);
                let stretch = parallel * (1.0 - beta) / beta_squared;
                let transformed =
                    core::array::from_fn(|i| separation[i] / beta + direction[i] * stretch);
                let range = norm(transformed);
                let x = k * range;
                if [range, x].iter().any(|v| !v.is_finite() || *v <= 0.0) {
                    return Err(invalid(
                        "convected power transformed separation is not representable",
                    ));
                }
                let projection = dot(mach, transformed.map(|v| v / range));
                let phase = projection * x;
                let moments = spherical_bessel_moments(x).ok_or_else(|| {
                    invalid("convected power sphere moments are not representable")
                })?;
                let (sin, cos) = det::sin_cos(phase);
                let sin = Dd::from_f64(sin);
                // Retain cos(phase)-1 when phase is too small for f64 cos to
                // differ from one. All three leading real decrements have
                // the same sign, including the convective vector moment.
                let cos_minus_one = if phase.abs() <= 1.0 {
                    let half_sin = Dd::from_f64(det::sin(0.5 * phase));
                    Dd::from_f64(-2.0) * half_sin * half_sin
                } else {
                    Dd::from_f64(cos) - Dd::ONE
                };
                let vector = Dd::from_f64(projection) * moments.j1;
                let pair = ConvectedPair {
                    real_minus_one: moments.j0_minus_one
                        + cos_minus_one
                        + moments.j0_minus_one * cos_minus_one
                        - vector * sin,
                    imaginary: sin + moments.j0_minus_one * sin + vector + vector * cos_minus_one,
                };
                if !finite_dd(pair.real_minus_one) || !finite_dd(pair.imaginary) {
                    return Err(nonfinite("convected power interference moment"));
                }
                pairs.push(pair);
            }
        }
        poll(cancelled)?;
        Ok(Self {
            acoustic,
            flow_velocity_m_s,
            frequency_hz,
            source_count,
            source_pairs,
            groups,
            members,
            pairs,
            gain,
            images: Vec::new(),
            image_group_size: 1,
        })
    }

    /// Original coherent density and adiabatic sound speed, without a substituted medium.
    #[must_use]
    pub const fn acoustic_state(&self) -> FluidAcousticState {
        self.acoustic
    }
    /// Uniform fluid velocity relative to the fixed laboratory sources [m/s].
    #[must_use]
    pub const fn flow_velocity_m_s(&self) -> [f64; 3] {
        self.flow_velocity_m_s
    }
    /// Physical source frequency in Hz; the fixed laboratory receiver sees the same frequency.
    #[must_use]
    pub const fn frequency_hz(&self) -> f64 {
        self.frequency_hz
    }
    /// Required number of volume-flow phasors in the original source order.
    #[must_use]
    pub const fn source_count(&self) -> usize {
        self.source_count
    }
    /// Exactly distinct direct/image positions; nearby sources remain distinct.
    #[must_use]
    pub fn position_count(&self) -> usize {
        self.groups.len()
    }

    /// Coherent transported power for the original peak volume-flow phasors.
    pub fn evaluate(&self, sources: &[C64]) -> Result<CompactRadiatedPower, AeroacError> {
        self.evaluate_cancellable(sources, || false)
    }

    /// Allocation-free coherent power with source and pair cancellation polls.
    pub fn evaluate_cancellable(
        &self,
        sources: &[C64],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<CompactRadiatedPower, AeroacError> {
        let power = self.coherent_power(sources, &mut cancelled)?;
        poll(&mut cancelled)?;
        finish(power, self.scope())
    }

    /// Explicit uncorrelated source factors, with E[xi_r conj(xi_s)]=delta_rs.
    /// Coherence is preserved within each factor; factors add powers. These
    /// are peak-unit spectral-line factors, not unscaled PSD/Hz values.
    pub fn evaluate_coherence_modes(
        &self,
        modes: &[&[C64]],
        budget: RadiatedPowerModeBudget,
    ) -> Result<CompactRadiatedPower, AeroacError> {
        self.evaluate_coherence_modes_cancellable(modes, budget, || false)
    }

    /// Bounded source factor bank; no partial total is returned on refusal.
    pub fn evaluate_coherence_modes_cancellable(
        &self,
        modes: &[&[C64]],
        budget: RadiatedPowerModeBudget,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<CompactRadiatedPower, AeroacError> {
        poll(&mut cancelled)?;
        let work = self.source_pairs.checked_mul(modes.len());
        if modes.len() > budget.max_modes
            || work.is_none_or(|work| work > budget.max_source_pair_mode_products)
        {
            return Err(invalid("convected power coherence work exceeds budget"));
        }
        let mut total = Dd::ZERO;
        for mode in modes {
            let power = self.coherent_power(mode, &mut cancelled)?;
            super::validate_power(power)?;
            total = total + power;
        }
        poll(&mut cancelled)?;
        finish(total, self.scope())
    }

    fn scope(&self) -> &'static str {
        if self.image_group_size > 1 {
            CONVECTED_MONOPOLE_IDEAL_BOUNDARY_POWER_SCOPE
        } else {
            CONVECTED_MONOPOLE_POWER_SCOPE
        }
    }

    fn coherent_power(
        &self,
        sources: &[C64],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Dd, AeroacError> {
        poll(cancelled)?;
        if sources.len() != self.source_count {
            return Err(invalid(
                "convected power source count differs from prepared geometry",
            ));
        }
        let mut baseline = ComplexDd::default();
        for &source in sources {
            poll(cancelled)?;
            if !finite(source) {
                return Err(nonfinite("convected power volume-flow amplitude"));
            }
            if self.image_group_size == 1 {
                baseline.add(source);
            }
        }
        if self.image_group_size > 1 {
            for image in &self.images {
                poll(cancelled)?;
                baseline.add(sources[image.source].scale(image.sign));
            }
        }
        baseline = baseline.scale(self.gain);
        let mut power = baseline.real_product(baseline);
        let mut pair_index = 0;
        for (index, a) in self.groups.iter().enumerate() {
            let a = self.group_sum(a, sources, cancelled)?;
            for b in &self.groups[..index] {
                poll(cancelled)?;
                let b = self.group_sum(b, sources, cancelled)?;
                let moment = self.pairs[pair_index];
                power = power
                    + Dd::from_f64(2.0)
                        * (a.real_product(b) * moment.real_minus_one
                            - a.imaginary_product(b) * moment.imaginary);
                pair_index += 1;
            }
        }
        Ok(if self.image_group_size == 1 {
            power
        } else {
            power / Dd::from_f64(self.image_group_size as f64)
        })
    }

    fn group_sum(
        &self,
        group: &PositionGroup,
        sources: &[C64],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<ComplexDd, AeroacError> {
        let mut result = ComplexDd::default();
        for &index in &self.members[group.members.clone()] {
            poll(cancelled)?;
            if self.image_group_size == 1 {
                result.add(sources[index]);
            } else {
                let image = self.images[index];
                result.add(sources[image.source].scale(image.sign));
            }
        }
        result = result.scale(self.gain);
        if !result.finite() {
            return Err(nonfinite("convected power grouped amplitude"));
        }
        Ok(result)
    }
}

fn norm(vector: [f64; 3]) -> f64 {
    det::hypot(det::hypot(vector[0], vector[1]), vector[2])
}
fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    (0..3).map(|i| a[i] * b[i]).sum()
}
