//! Causal pressure and particle velocity from stationary volume stresses.
//!
//! A compact symmetric stress moment Q [N m] gives
//! `p(x,t) = d_i d_j [Q_ij(t-r/c)/(4 pi r)]`. The derivatives act on
//! direction and spreading as well as retarded time. Consequently pressure
//! retains Q/r^3, Qdot/(c r^2) and Qddot/(c^2 r); the first two are not
//! radiation-zone terms. This is the time-domain form of the Lighthill
//! Green-function solution used by [`super::PreparedQuadrupoleObserver`].
//! Reference: Lighthill (1952), *On sound generated aerodynamically I*,
//! <https://doi.org/10.1098/rspa.1952.0060>.
//!
//! Momentum, `rho du/dt = -grad(p)`, also requires the causal integral of Q.
//! The resulting full vector velocity retains transverse near-field motion.
//! Source derivatives, initial history and interpolation are caller-owned:
//! the observer neither differentiates a noisy record nor fabricates a
//! stationary tone from a transient. All sources must be sampled at their
//! OWN retarded times for the SAME receiver time before coherent evaluation.
//! Exact ideal-wall images retain their own retarded times and the full
//! tensor/vector transformation. The original stationary, homogeneous,
//! lossless point-source limits apply; no finite-cell quadrature error,
//! motion, general scattering or feedback is inferred.

use super::{QUADRUPOLE_SCOPE, invalid, nonfinite, poll};
use crate::AeroacError;
use fs_material::fluid_acoustics::FluidAcousticState;
use fs_math::{dd::Dd, det};

mod reflection;
/// Doppler-correct pressure of uniform and accelerating physical stress moments.
pub mod moving;

/// One source's physical symmetric stress history at one retarded time.
///
/// Every array uses xx, yy, zz, xy, xz, yz; off-diagonal values are physical
/// tensor entries, not doubled engineering-shear entries. The four arrays
/// describe the SAME source clock and frame. They may come from an analytic
/// history or a separately qualified reconstruction of accepted flow states.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct TransientQuadrupoleSource {
    /// Integral of Q [N m s]. For an initially resting causal source use
    /// `integral(-infinity, retarded_time) Q(s) ds`. A different integration
    /// constant declares a different initial near field. Never reset at an
    /// audio callback boundary. A settled harmonic uses its matching
    /// zero-mean harmonic antiderivative instead.
    pub stress_integral_n_m_s: [f64; 6],
    /// Integrated physical stress Q [N m], including its trace.
    pub stress_moment_n_m: [f64; 6],
    /// First time derivative of Q [N m/s].
    pub stress_rate_n_m_s: [f64; 6],
    /// Second time derivative of Q [N m/s^2].
    pub stress_acceleration_n_m_s2: [f64; 6],
}

/// Coherent real field at one receiver time, including signed energy flux.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RetardedQuadrupoleField {
    /// Pressure perturbation [Pa], not a peak phasor or RMS value.
    pub pressure_pa: f64,
    /// Full Cartesian particle velocity [m/s], including induction terms.
    pub particle_velocity_m_s: [f64; 3],
    /// Instantaneous p u [W/m^2], formed AFTER summing the source fields.
    /// A pressure null can carry transverse velocity. This is not p^2/Z,
    /// a harmonic half-product, or the sum of individual source intensities.
    pub instantaneous_intensity_w_m2: [f64; 3],
    /// Original prescribed-source and propagation-model boundaries.
    pub scope: &'static str,
}

#[derive(Clone, Debug)]
struct Transfer {
    source_m: [f64; 3],
    delay_s: f64,
    direction: [f64; 3],
    // The three radial pressure orders: 1/(4 pi r^3), 1/(4 pi c r^2),
    // 1/(4 pi c^2 r). Velocity has one additional inverse length and rho.
    pressure_scales: [f64; 3],
    velocity_scales: [f64; 4],
    // A virtual receiver H^-1 x uses the original source tensor. Its velocity
    // returns to the physical frame through (H^-1)^T. None is the direct path.
    inverse_image: Option<[[f64; 3]; 3]>,
}

impl Transfer {
    fn new(
        acoustic: FluidAcousticState,
        observer_m: [f64; 3],
        source_m: [f64; 3],
    ) -> Result<Self, AeroacError> {
        if source_m.iter().any(|x| !x.is_finite()) {
            return Err(nonfinite("transient quadrupole source position"));
        }
        let displacement = core::array::from_fn::<_, 3, _>(|i| observer_m[i] - source_m[i]);
        let range = det::hypot(
            det::hypot(displacement[0], displacement[1]),
            displacement[2],
        );
        let c = acoustic.sound_speed_m_s();
        let rho = acoustic.density_kg_m3();
        let green_scale = (0.25 / core::f64::consts::PI) / range;
        let pressure_scales = [
            (green_scale / range) / range,
            (green_scale / range) / c,
            (green_scale / c) / c,
        ];
        let velocity_scales = [
            (pressure_scales[0] / range) / rho,
            (pressure_scales[1] / range) / rho,
            (pressure_scales[2] / range) / rho,
            (pressure_scales[2] / c) / rho,
        ];
        let delay_s = range / c;
        if [range, delay_s]
            .iter()
            .chain(&pressure_scales)
            .chain(&velocity_scales)
            .any(|x| !x.is_finite() || *x <= 0.0)
        {
            return Err(invalid(
                "transient quadrupole separation or field coefficients are not representable",
            ));
        }
        Ok(Self {
            source_m,
            delay_s,
            direction: displacement.map(|x| x / range),
            pressure_scales,
            velocity_scales,
            inverse_image: None,
        })
    }
}

/// Prepared coherent observation of stationary transient volume-stress cells.
///
/// Only each cell is compact; the whole source array may span wavelengths.
/// Preparation allocates once under an explicit source-count budget. Each
/// evaluation is allocation-free O(P), where P is the admitted physical-source
/// path count, with cached geometry/material scales, no transcendentals and
/// fixed-order double-double field sums. Exact ideal-wall images require
/// distinct retarded histories, so they are never collapsed to one effective
/// delay. Empty source arrays are silence. Neither the preparation budget nor
/// a successful evaluation validates source compactness, bandwidth or physics.
#[derive(Clone, Debug)]
pub struct PreparedTransientQuadrupoleObserver {
    acoustic: FluidAcousticState,
    observer_m: [f64; 3],
    // Direct paths come first; each image then repeats physical source order.
    transfers: Vec<Transfer>,
    physical_source_count: usize,
    reflection_planes: Vec<super::IdealAcousticPlane>,
    reflection_wedge_order: Option<usize>,
}

impl PreparedTransientQuadrupoleObserver {
    /// Prepare physical source and receiver positions in one Cartesian frame.
    pub fn new(
        acoustic: FluidAcousticState,
        observer_m: [f64; 3],
        source_positions_m: &[[f64; 3]],
        max_sources: usize,
    ) -> Result<Self, AeroacError> {
        Self::new_cancellable(
            acoustic,
            observer_m,
            source_positions_m,
            max_sources,
            || false,
        )
    }

    /// Poll before each source and publication. Singular geometry, invalid
    /// coefficients, allocation failure or cancellation publishes no observer.
    pub fn new_cancellable(
        acoustic: FluidAcousticState,
        observer_m: [f64; 3],
        source_positions_m: &[[f64; 3]],
        max_sources: usize,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, AeroacError> {
        poll(&mut cancelled)?;
        if source_positions_m.len() > max_sources {
            return Err(invalid(
                "transient quadrupole source count exceeds preparation budget",
            ));
        }
        if observer_m.iter().any(|x| !x.is_finite()) {
            return Err(nonfinite("transient quadrupole observer position"));
        }
        let mut transfers = Vec::new();
        transfers
            .try_reserve_exact(source_positions_m.len())
            .map_err(|_| invalid("transient quadrupole allocation exceeds available memory"))?;
        for source in source_positions_m {
            poll(&mut cancelled)?;
            transfers.push(Transfer::new(acoustic, observer_m, *source)?);
        }
        poll(&mut cancelled)?;
        Ok(Self {
            acoustic,
            observer_m,
            transfers,
            physical_source_count: source_positions_m.len(),
            reflection_planes: Vec::new(),
            reflection_wedge_order: None,
        })
    }

    /// Original acoustic snapshot shared by every prepared path.
    pub const fn acoustic(&self) -> FluidAcousticState {
        self.acoustic
    }

    /// Fixed receiver position [m] in the original Cartesian frame.
    pub const fn observer_position_m(&self) -> [f64; 3] {
        self.observer_m
    }

    /// Number of original physical sources, unchanged by ideal-wall images.
    pub fn source_count(&self) -> usize {
        self.physical_source_count
    }

    /// Number of independently retarded paths, including every direct path.
    pub fn path_count(&self) -> usize {
        self.transfers.len()
    }

    /// Original source position [m], retaining the caller's source order.
    pub fn source_position_m(&self, index: usize) -> Option<[f64; 3]> {
        (index < self.source_count()).then(|| self.transfers[index].source_m)
    }

    /// Direct-path travel time r/c [s] for the original physical source.
    /// Reflections retain additional distinct delays through `path_delay_s`.
    pub fn source_delay_s(&self, index: usize) -> Option<f64> {
        (index < self.source_count()).then(|| self.transfers[index].delay_s)
    }

    /// Original source index for a path in image-major/source-major order.
    /// The first `source_count()` paths are the direct group.
    pub fn path_source_index(&self, index: usize) -> Option<usize> {
        self.transfers
            .get(index)
            .map(|_| index % self.source_count())
    }

    /// Travel time [s] of one exact direct or reflected path.
    pub fn path_delay_s(&self, index: usize) -> Option<f64> {
        self.transfers.get(index).map(|t| t.delay_s)
    }

    /// Minimum and maximum delay [s] over ALL paths of one physical source.
    /// Use these bounds to admit complete source-history windows before
    /// observation; no echo is discarded or replaced by an effective delay.
    pub fn source_delay_bounds_s(&self, index: usize) -> Option<[f64; 2]> {
        let direct = self.source_delay_s(index)?;
        Some(
            self.transfers
                .iter()
                .skip(index)
                .step_by(self.source_count())
                .fold([direct; 2], |[lo, hi], path| {
                    [lo.min(path.delay_s), hi.max(path.delay_s)]
                }),
        )
    }

    /// Physical ideal boundaries, in their original preparation order.
    pub fn reflection_planes(&self) -> &[super::IdealAcousticPlane] {
        &self.reflection_planes
    }

    /// Exact wedge opening pi/m, or None for free space/perpendicular walls.
    pub const fn reflection_wedge_order(&self) -> Option<usize> {
        self.reflection_wedge_order
    }

    /// Original source and propagation-model boundaries of every field.
    pub const fn scope(&self) -> &'static str {
        if self.reflection_planes.is_empty() {
            QUADRUPOLE_SCOPE
        } else if self.reflection_wedge_order.is_some() {
            reflection::WEDGE_SCOPE
        } else {
            reflection::ORTHOGONAL_SCOPE
        }
    }

    /// Sum source fields at one receiver time before computing energy flux.
    ///
    /// Every input is already at its own path's retarded time. Supply
    /// `path_count()` samples in the order exposed by `path_source_index` and
    /// `path_delay_s`; reflected paths must not reuse a direct-path sample.
    /// The caller owns interpolation, antialiasing and source bandwidth.
    pub fn evaluate_at_retarded(
        &self,
        sources: &[TransientQuadrupoleSource],
    ) -> Result<RetardedQuadrupoleField, AeroacError> {
        self.evaluate_at_retarded_cancellable(sources, || false)
    }

    /// Allocation-free coherent evaluation; poll at most every 256 paths.
    /// Refusal returns no partial field and does not mutate any source state.
    pub fn evaluate_at_retarded_cancellable(
        &self,
        sources: &[TransientQuadrupoleSource],
        cancelled: impl FnMut() -> bool,
    ) -> Result<RetardedQuadrupoleField, AeroacError> {
        if sources.len() != self.path_count() {
            return Err(invalid(
                "one retarded history is required per transient quadrupole path",
            ));
        }
        self.evaluate_sources(|path, _, _| Ok(sources[path]), cancelled)
    }

    /// Observe histories at one receiver time, applying each path delay once.
    ///
    /// The callback receives (physical source index, emission time in seconds)
    /// in image-major/source-major order. Reflections observe the SAME source
    /// at each path's distinct retarded time; they do not advance or duplicate
    /// its state. The callback supplies a physical value and its matching
    /// derivatives/integral from that source's own history. The observer
    /// performs no interpolation, extrapolation, differentiation or source
    /// advancement. Callbacks should only observe history; their own side
    /// effects cannot be undone if a later callback or cancellation refuses.
    ///
    /// Emission times are rounded to f64 before invoking the callback. A
    /// positive rounded delay does not certify phase accuracy: at a large
    /// epoch, partially resolved timestamps can perturb source interference.
    /// Use [`Self::evaluate_at_clock_dd`] to retain short path delays and
    /// sample intervals on the original double-double source clock.
    pub fn evaluate_at_time(
        &self,
        receiver_time_s: f64,
        source_at: impl FnMut(usize, f64) -> Result<TransientQuadrupoleSource, AeroacError>,
    ) -> Result<RetardedQuadrupoleField, AeroacError> {
        self.evaluate_at_time_cancellable(receiver_time_s, source_at, || false)
    }

    /// Allocation-free history observation with bounded cancellation polls.
    /// Retarded timestamps that overflow or round back to the receiver time
    /// refuse without a field. The f64 precision limits of
    /// [`Self::evaluate_at_time`] still apply to admitted timestamps.
    pub fn evaluate_at_time_cancellable(
        &self,
        receiver_time_s: f64,
        mut source_at: impl FnMut(usize, f64) -> Result<TransientQuadrupoleSource, AeroacError>,
        cancelled: impl FnMut() -> bool,
    ) -> Result<RetardedQuadrupoleField, AeroacError> {
        if !receiver_time_s.is_finite() {
            return Err(nonfinite("transient quadrupole receiver time"));
        }
        self.evaluate_sources(
            |_, index, delay_s| {
                let emission_time_s =
                    (Dd::from_f64(receiver_time_s) - Dd::from_f64(delay_s)).to_f64();
                if !emission_time_s.is_finite() || emission_time_s >= receiver_time_s {
                    return Err(invalid(
                        "transient quadrupole retarded time is not representable",
                    ));
                }
                source_at(index, emission_time_s)
            },
            cancelled,
        )
    }

    /// Observe histories on the original double-double source clock.
    ///
    /// Both the receiver time and callback emission times retain their low
    /// components in seconds. This preserves relative propagation phase when
    /// a short delay or sample interval lies below the f64 spacing of a large
    /// epoch. It composes directly with `HarmonicSampleClock::time_at`; no
    /// phase origin, source alignment or history interpolation is inferred.
    ///
    /// The callback receives (original source index, receiver time - r/c).
    /// Keep that timestamp in DD until subtracting its history epoch or
    /// evaluating its phase; rounding the absolute time to f64 loses the
    /// information again. Integral, derivatives and all history ownership
    /// requirements of [`Self::evaluate_at_time`] remain unchanged.
    pub fn evaluate_at_clock_dd(
        &self,
        receiver_time_s: Dd,
        source_at: impl FnMut(usize, Dd) -> Result<TransientQuadrupoleSource, AeroacError>,
    ) -> Result<RetardedQuadrupoleField, AeroacError> {
        self.evaluate_at_clock_dd_cancellable(receiver_time_s, source_at, || false)
    }

    /// Allocation-free DD-clock observation with the original cancellation
    /// bounds. Nonfinite or DD-unresolved retarded times refuse without a
    /// field. As with the f64 callback, user callback side effects cannot be
    /// undone if a later source fails or cancellation is requested.
    pub fn evaluate_at_clock_dd_cancellable(
        &self,
        receiver_time_s: Dd,
        mut source_at: impl FnMut(usize, Dd) -> Result<TransientQuadrupoleSource, AeroacError>,
        cancelled: impl FnMut() -> bool,
    ) -> Result<RetardedQuadrupoleField, AeroacError> {
        if !receiver_time_s.hi.is_finite() || !receiver_time_s.lo.is_finite() {
            return Err(nonfinite("transient quadrupole receiver clock"));
        }
        let receiver_time_s = Dd::from_pair(receiver_time_s.hi, receiver_time_s.lo);
        if !receiver_time_s.hi.is_finite() || !receiver_time_s.lo.is_finite() {
            return Err(nonfinite("transient quadrupole receiver clock"));
        }
        self.evaluate_sources(
            |_, index, delay_s| {
                let emission_time_s = receiver_time_s - Dd::from_f64(delay_s);
                if !emission_time_s.hi.is_finite()
                    || !emission_time_s.lo.is_finite()
                    || !emission_time_s.lt(receiver_time_s)
                {
                    return Err(invalid(
                        "transient quadrupole retarded time is unresolved on the DD clock",
                    ));
                }
                source_at(index, emission_time_s)
            },
            cancelled,
        )
    }

    fn evaluate_sources(
        &self,
        mut source_at: impl FnMut(usize, usize, f64) -> Result<TransientQuadrupoleSource, AeroacError>,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<RetardedQuadrupoleField, AeroacError> {
        poll(&mut cancelled)?;
        let mut sums = [Dd::ZERO; 4];
        for (index, transfer) in self.transfers.iter().enumerate() {
            if index.is_multiple_of(256) {
                poll(&mut cancelled)?;
            }
            let source = source_at(index, index % self.source_count(), transfer.delay_s)?;
            if source
                .stress_integral_n_m_s
                .iter()
                .chain(&source.stress_moment_n_m)
                .chain(&source.stress_rate_n_m_s)
                .chain(&source.stress_acceleration_n_m_s2)
                .any(|x| !x.is_finite())
            {
                return Err(nonfinite("transient quadrupole source history"));
            }
            let direction = transfer.direction.map(Dd::from_f64);
            let integral = Projection::new(source.stress_integral_n_m_s, direction);
            let moment = Projection::new(source.stress_moment_n_m, direction);
            let rate = Projection::new(source.stress_rate_n_m_s, direction);
            let acceleration = Projection::new(source.stress_acceleration_n_m_s2, direction);
            let p = transfer.pressure_scales.map(Dd::from_f64);
            sums[0] = sums[0]
                + moment.pressure_near() * p[0]
                + rate.pressure_near() * p[1]
                + acceleration.radial * p[2];
            let v = transfer.velocity_scales.map(Dd::from_f64);
            let mut velocity = if transfer.inverse_image.is_none() {
                [sums[1], sums[2], sums[3]]
            } else {
                [Dd::ZERO; 3]
            };
            for axis in 0..3 {
                let n = direction[axis];
                // -grad of the retarded time integral of pressure. Direction
                // derivatives are essential: Qn is not generally parallel to n.
                velocity[axis] = velocity[axis]
                    + integral.velocity_near(n, axis) * v[0]
                    + moment.velocity_near(n, axis) * v[1]
                    + ((rate.radial * Dd::from_f64(6.0) - rate.trace) * n
                        - rate.vector[axis] * Dd::from_f64(2.0))
                        * v[2]
                    + acceleration.radial * n * v[3];
            }
            if let Some(inverse) = transfer.inverse_image {
                for axis in 0..3 {
                    for virtual_axis in 0..3 {
                        sums[axis + 1] = sums[axis + 1]
                            + Dd::from_f64(inverse[virtual_axis][axis]) * velocity[virtual_axis];
                    }
                }
            } else {
                for axis in 0..3 {
                    sums[axis + 1] = velocity[axis];
                }
            }
        }
        poll(&mut cancelled)?;
        let [pressure_pa, ux, uy, uz] = sums.map(Dd::to_f64);
        let particle_velocity_m_s = [ux, uy, uz];
        let instantaneous_intensity_w_m2 = particle_velocity_m_s.map(|u| pressure_pa * u);
        if !pressure_pa.is_finite()
            || particle_velocity_m_s
                .iter()
                .chain(&instantaneous_intensity_w_m2)
                .any(|x| !x.is_finite())
        {
            return Err(nonfinite(
                "accumulated transient quadrupole field or intensity",
            ));
        }
        Ok(RetardedQuadrupoleField {
            pressure_pa,
            particle_velocity_m_s,
            instantaneous_intensity_w_m2,
            scope: self.scope(),
        })
    }
}

struct Projection {
    trace: Dd,
    radial: Dd,
    vector: [Dd; 3],
}

impl Projection {
    fn new(tensor: [f64; 6], n: [Dd; 3]) -> Self {
        let [xx, yy, zz, xy, xz, yz] = tensor.map(Dd::from_f64);
        let vector = [
            xx * n[0] + xy * n[1] + xz * n[2],
            xy * n[0] + yy * n[1] + yz * n[2],
            xz * n[0] + yz * n[1] + zz * n[2],
        ];
        Self {
            trace: xx + yy + zz,
            radial: n[0] * vector[0] + n[1] * vector[1] + n[2] * vector[2],
            vector,
        }
    }

    fn pressure_near(&self) -> Dd {
        self.radial * Dd::from_f64(3.0) - self.trace
    }

    fn velocity_near(&self, n: Dd, axis: usize) -> Dd {
        (self.radial * Dd::from_f64(15.0) - self.trace * Dd::from_f64(3.0)) * n
            - self.vector[axis] * Dd::from_f64(6.0)
    }
}
