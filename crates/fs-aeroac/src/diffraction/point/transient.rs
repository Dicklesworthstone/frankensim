//! Causal spherical diffraction for the original real volume-flow clock.
//!
//! The Biot--Tolstoy half-plane kernel is the same one used by the harmonic
//! point owner (Ouis 2002, equations 3--6, doi:10.1006/jsvi.2000.3553).
//! Put f(s,t)=Q(t-sqrt(L²+s²)/c)/sqrt(L²+s²). Each direct/image potential is
//! f(0)/(8*pi)+H(q)*(Q(t-D/c)/D-f(0))/(4*pi)
//! - integral_0^infinity partial_s(f)*atan(q/s) ds/(4*pi²).
//! Pressure is rho*partial_t of this expression and velocity its negative
//! spatial gradient. The REAL causal integral retains the complete edge tail.
//! Analytic derivatives, including the q=0 limit, share the harmonic geometry.
//!
//! Original Q and dQ/dt at uniform knots define a C1 cubic Hermite history.
//! The prepared taps integrate these basis functions, NOT sampled harmonic
//! transfers. Eight cached real convolution kernels yield p and all u components.
//! Zero padding in the original FFT performs linear convolution without a
//! periodic source, acausal pre-echo or an imposed frequency-grid propagator.

use super::{
    AeroacError, HalfPlaneFrame, HalfPlaneMonopoleLimits, HalfPlaneReceiver, IdealPlaneBoundary,
    PointGeometry, QuadratureLimits, Spectrum, integrate_spectrum, invalid, poll,
};
use crate::curle3d::moving::TimedCompactSource;
use core::f64::consts::PI;
use fs_material::fluid_acoustics::FluidAcousticState;
use fs_math::{c64::C64, dd::Dd, det};

mod convolution;
use convolution::Convolution;

const MAX_SAMPLES: usize = 65_536;
const MAX_TOTAL_EVALUATIONS: usize = 16_777_216;
const SCOPE: &str = "Causal stationary volume-flow history diffracted by one infinite \
ideal rigid or pressure-release half-plane in an original homogeneous lossless fluid. \
C1 cubic interpolation of original Q and dQ/dt; complete edge memory within the \
declared observation horizon. Full pressure and velocity, instantaneous pressure \
work and perturbation energy. Quadrature estimates exclude interpolation, geometry \
and floating-point error. No source motion, mean flow, force injection, finite \
barrier, wall impedance, absorption, feedback or measured-SPL qualification.";

/// Bounded causal preparation on one explicit uniform physical clock.
#[derive(Clone, Copy, Debug)]
pub struct HalfPlaneTransientLimits {
    /// Each tap acts on [Q, dt*dQ/dt], both in m³/s. The two pressure
    /// coefficients and six velocity coefficients use the original norm gates.
    pub quadrature: HalfPlaneMonopoleLimits,
    /// Total integrand evaluations across every tap, in 1..=16_777_216.
    pub maximum_total_evaluations: usize,
    /// Conservative retained + preparation/application peak numerical payload.
    /// The caller's borrowed source record and allocator metadata are excluded.
    pub maximum_payload_bytes: usize,
}

/// Numerical preparation work; no source interpolation or geometry certificate.
#[derive(Clone, Debug)]
pub struct HalfPlaneTransientDiagnostics {
    pub evaluations: usize,
    pub maximum_depth_reached: usize,
    pub maximum_retained_panels: usize,
    pub minimum_edge_delay_s: f64,
    pub minimum_arrival_delay_s: f64,
    pub fft_length: usize,
    pub peak_payload_bytes: usize,
    /// Sums of per-tap embedded error estimates. Layout matches coefficients().
    /// Multiply Q channels by max|Q| and derivative channels by max|dt*dQ/dt|.
    pub quadrature_error_l1: [f64; 8],
}

/// A real instantaneous observation on the source's original clock.
#[derive(Clone, Copy, Debug)]
pub struct HalfPlaneTransientField {
    pub receiver_time_s: f64,
    pub pressure_pa: f64,
    pub particle_velocity_m_s: [f64; 3],
    /// p*u, without a harmonic factor of one half. Sum coherent fields first.
    pub instantaneous_intensity_w_m2: [f64; 3],
    /// p²/(2K) + rho*|u|²/2 in the original stationary medium.
    pub perturbation_energy_j_m3: f64,
    pub pressure_quadrature_error_pa: f64,
    pub velocity_quadrature_error_m_s: [f64; 3],
    pub scope: &'static str,
}

/// Finite observation horizon, with no dropping of old source samples.
///
/// Every supported output uses all accepted source knots since the quiescent
/// start. A longer record is refused; this is not an indefinitely running
/// finite-memory approximation to the diffraction tail.
#[derive(Clone, Debug)]
pub struct PreparedHalfPlaneTransient {
    acoustic: FluidAcousticState,
    frame: HalfPlaneFrame,
    boundary: IdealPlaneBoundary,
    source_position_m: [f64; 3],
    receiver: HalfPlaneReceiver,
    sample_interval_s: f64,
    coefficients: Vec<[f64; 8]>,
    errors: Vec<[f64; 8]>,
    convolution: Convolution,
    diagnostics: HalfPlaneTransientDiagnostics,
}

impl PreparedHalfPlaneTransient {
    /// maximum_samples includes the original zero-strength start knot.
    /// The source must be in open fluid; receivers may select an exact sheet face.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        acoustic: FluidAcousticState,
        frame: HalfPlaneFrame,
        boundary: IdealPlaneBoundary,
        source_position_m: [f64; 3],
        receiver: HalfPlaneReceiver,
        sample_interval_s: f64,
        maximum_samples: usize,
        limits: HalfPlaneTransientLimits,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, AeroacError> {
        poll(&mut cancelled)?;
        let controls = QuadratureLimits::from(limits.quadrature);
        controls.admit()?;
        if !sample_interval_s.is_finite()
            || sample_interval_s <= 0.0
            || !(1..=MAX_SAMPLES).contains(&maximum_samples)
            || !(1..=MAX_TOTAL_EVALUATIONS).contains(&limits.maximum_total_evaluations)
        {
            return Err(invalid(
                "causal diffraction needs a positive clock and bounded horizon/work",
            ));
        }
        let step_length = acoustic.sound_speed_m_s() * sample_interval_s;
        let horizon = sample_interval_s * maximum_samples as f64;
        if !step_length.is_finite() || step_length <= 0.0 || !horizon.is_finite() {
            return Err(invalid(
                "causal diffraction time or propagation span is unrepresentable",
            ));
        }
        let (fft_length, fft_workspace) = Convolution::layout(maximum_samples)?;
        let peak = fft_length
            .checked_mul(11 * 16)
            .and_then(|v| {
                maximum_samples
                    .checked_mul(512)
                    .and_then(|n| v.checked_add(n))
            })
            .and_then(|v| v.checked_add(fft_workspace))
            .and_then(|v| {
                (controls.maximum_evaluations / 30 + 6)
                    .checked_mul(512)
                    .and_then(|q| v.checked_add(q))
            })
            .and_then(|v| v.checked_add(4096))
            .ok_or_else(|| invalid("causal diffraction payload overflow"))?;
        if peak > limits.maximum_payload_bytes {
            return Err(invalid(
                "causal diffraction exceeds the declared numerical payload",
            ));
        }
        let source = frame.local_receiver(HalfPlaneReceiver::WorldPosition(source_position_m))?;
        let point = frame.local_receiver(receiver)?;
        let geometry = PointGeometry::new(source, point)?;
        let edge_delay = geometry.distance / acoustic.sound_speed_m_s();
        let mut arrival = edge_delay;
        for j in 0..2 {
            if geometry.q[j] > 0.0 {
                arrival = arrival.min(geometry.direct_distance[j] / acoustic.sound_speed_m_s());
            }
        }
        if !edge_delay.is_finite() || edge_delay <= 0.0 || !arrival.is_finite() || arrival <= 0.0 {
            return Err(invalid("causal diffraction delays are unrepresentable"));
        }
        let mut coefficients = Vec::new();
        let mut errors = Vec::new();
        coefficients
            .try_reserve_exact(maximum_samples)
            .map_err(|_| invalid("causal diffraction coefficient allocation failed"))?;
        errors
            .try_reserve_exact(maximum_samples)
            .map_err(|_| invalid("causal diffraction error allocation failed"))?;
        let mut evaluations = 0usize;
        let mut maximum_depth = 0usize;
        let mut maximum_panels = 0usize;
        let mut error_l1 = [Dd::ZERO; 8];
        for lag in 0..maximum_samples {
            poll(&mut cancelled)?;
            let kernel = TapKernel::new(
                &geometry,
                frame.axes(),
                boundary.sign(),
                acoustic,
                sample_interval_s,
                lag,
            )?;
            let (values, error) = if kernel.regions() == 0 {
                (kernel.base, [0.0; 8])
            } else {
                let remaining = limits
                    .maximum_total_evaluations
                    .checked_sub(evaluations)
                    .ok_or_else(|| invalid("causal diffraction total work exhausted"))?;
                if remaining == 0 {
                    return Err(invalid("causal diffraction total work exhausted"));
                }
                let active = QuadratureLimits {
                    maximum_evaluations: controls.maximum_evaluations.min(remaining),
                    ..controls
                };
                let result = integrate_spectrum::<8, 2>(&kernel, [0.0; 8], active, &mut cancelled)?;
                evaluations += result.evaluations;
                maximum_depth = maximum_depth.max(result.maximum_depth_reached);
                maximum_panels = maximum_panels.max(result.retained_panels);
                (result.value.map(|v| v.re), result.error)
            };
            if values.iter().any(|v| !v.is_finite()) {
                return Err(invalid("causal diffraction coefficient is unrepresentable"));
            }
            coefficients.push(values);
            errors.push(error);
            for j in 0..8 {
                error_l1[j] = error_l1[j] + Dd::from_f64(error[j]);
            }
        }
        let error_l1 = error_l1.map(Dd::to_f64);
        if error_l1.iter().any(|v| !v.is_finite()) {
            return Err(invalid("causal diffraction error sum is unrepresentable"));
        }
        let convolution =
            Convolution::prepare(&coefficients, fft_length, fft_workspace, &mut cancelled)?;
        poll(&mut cancelled)?;
        Ok(Self {
            acoustic,
            frame,
            boundary,
            source_position_m,
            receiver,
            sample_interval_s,
            coefficients,
            errors,
            convolution,
            diagnostics: HalfPlaneTransientDiagnostics {
                evaluations,
                maximum_depth_reached: maximum_depth,
                maximum_retained_panels: maximum_panels,
                minimum_edge_delay_s: edge_delay,
                minimum_arrival_delay_s: arrival,
                fft_length,
                peak_payload_bytes: peak,
                quadrature_error_l1: error_l1,
            },
        })
    }
    pub const fn acoustic(&self) -> FluidAcousticState {
        self.acoustic
    }
    pub const fn frame(&self) -> HalfPlaneFrame {
        self.frame
    }
    pub const fn boundary(&self) -> IdealPlaneBoundary {
        self.boundary
    }
    pub const fn source_position_m(&self) -> [f64; 3] {
        self.source_position_m
    }
    pub const fn receiver(&self) -> HalfPlaneReceiver {
        self.receiver
    }
    pub const fn sample_interval_s(&self) -> f64 {
        self.sample_interval_s
    }
    pub fn maximum_samples(&self) -> usize {
        self.coefficients.len()
    }
    pub fn diagnostics(&self) -> &HalfPlaneTransientDiagnostics {
        &self.diagnostics
    }
    /// Original real tap layout: p/Q, p/(dt*Qdot), ux/Q, uy/Q, uz/Q,
    /// ux/(dt*Qdot), uy/(dt*Qdot), uz/(dt*Qdot). Taps carry no source amplitude.
    pub fn coefficients(&self) -> &[[f64; 8]] {
        &self.coefficients
    }

    /// Allocation-free direct evaluation, useful for individual receivers/times.
    /// Uses the same taps as the O(N log N) complete-record operation.
    pub fn observe_at(
        &self,
        source: &[TimedCompactSource],
        index: usize,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<HalfPlaneTransientField, AeroacError> {
        self.admit_record(source, &mut cancelled)?;
        if index >= source.len() {
            return Err(invalid(
                "causal observation index is outside its source record",
            ));
        }
        let mut value = [Dd::ZERO; 4];
        let mut error = [Dd::ZERO; 4];
        for lag in 0..=index {
            if lag % 64 == 0 {
                poll(&mut cancelled)?;
            }
            let q = source[index - lag].strength.volume_flow_m3_s;
            let d = self.sample_interval_s * source[index - lag].strength.volume_flow_rate_m3_s2;
            for component in 0..4 {
                let (a, b) = channels(component);
                value[component] = value[component]
                    + Dd::from_f64(self.coefficients[lag][a]) * Dd::from_f64(q)
                    + Dd::from_f64(self.coefficients[lag][b]) * Dd::from_f64(d);
                error[component] = error[component]
                    + Dd::from_f64(self.errors[lag][a]) * Dd::from_f64(q.abs())
                    + Dd::from_f64(self.errors[lag][b]) * Dd::from_f64(d.abs());
            }
        }
        poll(&mut cancelled)?;
        self.field(
            source[index].source_time_s,
            value.map(Dd::to_f64),
            error.map(Dd::to_f64),
        )
    }

    /// Observe every accepted time by cached, zero-padded LINEAR convolution.
    ///
    /// The first source knot must have Q=Qdot=0, declaring a quiescent past.
    /// Input force components must be exactly zero, so they cannot be silently
    /// dropped from the original CompactSourceStrength. No future sample enters
    /// an output. A record extending beyond the prepared horizon is refused.
    /// Appending explicitly supplied zero source knots observes the remaining
    /// tail within that same horizon; no automatic termination is invented.
    pub fn observe_record(
        &self,
        source: &[TimedCompactSource],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Vec<HalfPlaneTransientField>, AeroacError> {
        let amplitudes = self.admit_record(source, &mut cancelled)?;
        let values =
            self.convolution
                .apply(source, self.sample_interval_s, amplitudes, &mut cancelled)?;
        // Hermite support starts in the interval BEFORE the first nonzero
        // knot. A known quiescent prefix must retain exact causal zeros even
        // when the record-wide FFT has floating-point spectral roundoff.
        let support_start = source
            .iter()
            .position(|sample| {
                sample.strength.volume_flow_m3_s != 0.0
                    || sample.strength.volume_flow_rate_m3_s2 != 0.0
            })
            .map(|index| index.saturating_sub(1));
        // A conservative record-wide bound on the embedded quadrature estimate.
        // FFT and other floating-point roundoff remain outside this estimate.
        let errors = core::array::from_fn(|component| {
            let (a, b) = channels(component);
            (Dd::from_f64(self.diagnostics.quadrature_error_l1[a]) * Dd::from_f64(amplitudes[0])
                + Dd::from_f64(self.diagnostics.quadrature_error_l1[b])
                    * Dd::from_f64(amplitudes[1]))
            .to_f64()
        });
        let mut out = Vec::new();
        out.try_reserve_exact(source.len())
            .map_err(|_| invalid("causal diffraction output allocation failed"))?;
        for (index, (sample, value)) in source.iter().zip(values).enumerate() {
            if index % 64 == 0 {
                poll(&mut cancelled)?;
            }
            // Suppress only the analytically empty prefix, never a physical
            // edge tail or a later zero in the emitted source history.
            let before_arrival = support_start.is_none_or(|start| {
                // The taps use the declared relative clock. Subtracting
                // rounded absolute epochs could move the first arrival.
                ((Dd::from_f64(index as f64) - Dd::from_f64(start as f64))
                    * Dd::from_f64(self.sample_interval_s))
                .to_f64()
                    < self.diagnostics.minimum_arrival_delay_s
            });
            let value = if before_arrival { [0.0; 4] } else { value };
            out.push(self.field(sample.source_time_s, value, errors)?);
        }
        poll(&mut cancelled)?;
        Ok(out)
    }

    fn admit_record(
        &self,
        source: &[TimedCompactSource],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<[f64; 2], AeroacError> {
        poll(cancelled)?;
        if source.is_empty() || source.len() > self.maximum_samples() {
            return Err(invalid(
                "causal source record is empty or exceeds its complete-memory horizon",
            ));
        }
        let start = source[0].source_time_s;
        if !start.is_finite()
            || source[0].strength.volume_flow_m3_s != 0.0
            || source[0].strength.volume_flow_rate_m3_s2 != 0.0
        {
            return Err(invalid(
                "causal source must declare its quiescent zero-strength start",
            ));
        }
        let mut amplitudes = [0.0_f64; 2];
        for (i, sample) in source.iter().enumerate() {
            if i % 64 == 0 {
                poll(cancelled)?;
            }
            let strength = sample.strength;
            if !sample.source_time_s.is_finite()
                || !strength.volume_flow_m3_s.is_finite()
                || !strength.volume_flow_rate_m3_s2.is_finite()
                || strength
                    .force_on_fluid_n
                    .iter()
                    .chain(&strength.force_rate_on_fluid_n_s)
                    .any(|v| *v != 0.0)
            {
                return Err(invalid(
                    "causal half-plane input requires finite pure volume-flow history",
                ));
            }
            let expected = (Dd::from_f64(start)
                + Dd::from_f64(self.sample_interval_s) * Dd::from_f64(i as f64))
            .to_f64();
            // No absolute-epoch tolerance: physical dt must remain resolvable.
            if sample.source_time_s != expected
                || (i > 0 && sample.source_time_s <= source[i - 1].source_time_s)
            {
                return Err(invalid(
                    "causal source knots must retain the prepared uniform physical clock",
                ));
            }
            let derivative = self.sample_interval_s * strength.volume_flow_rate_m3_s2;
            if !derivative.is_finite() {
                return Err(invalid("scaled source derivative is unrepresentable"));
            }
            amplitudes[0] = amplitudes[0].max(strength.volume_flow_m3_s.abs());
            amplitudes[1] = amplitudes[1].max(derivative.abs());
        }
        Ok(amplitudes)
    }

    fn field(
        &self,
        time: f64,
        value: [f64; 4],
        error: [f64; 4],
    ) -> Result<HalfPlaneTransientField, AeroacError> {
        let pressure = value[0];
        let velocity = [value[1], value[2], value[3]];
        let intensity = velocity.map(|v| (Dd::from_f64(pressure) * Dd::from_f64(v)).to_f64());
        let mut energy = Dd::from_f64(pressure)
            * Dd::from_f64(pressure)
            * Dd::from_f64(0.5 * self.acoustic.adiabatic_compressibility_per_pa());
        for v in velocity {
            energy = energy
                + Dd::from_f64(v)
                    * Dd::from_f64(v)
                    * Dd::from_f64(0.5 * self.acoustic.density_kg_m3());
        }
        let energy = energy.to_f64();
        if value
            .iter()
            .chain(&error)
            .chain(&intensity)
            .any(|v| !v.is_finite())
            || !energy.is_finite()
            || energy < 0.0
            || error.iter().any(|v| *v < 0.0)
        {
            return Err(invalid("causal diffracted field/work is unrepresentable"));
        }
        Ok(HalfPlaneTransientField {
            receiver_time_s: time,
            pressure_pa: pressure,
            particle_velocity_m_s: velocity,
            instantaneous_intensity_w_m2: intensity,
            perturbation_energy_j_m3: energy,
            pressure_quadrature_error_pa: error[0],
            velocity_quadrature_error_m_s: [error[1], error[2], error[3]],
            scope: SCOPE,
        })
    }
}

fn channels(component: usize) -> (usize, usize) {
    if component == 0 {
        (0, 1)
    } else {
        (1 + component, 4 + component)
    }
}

/// Value, first and second derivatives with respect to the dimensionless knot
/// coordinate. The second derivative may jump at knots; all integration
/// regions split there. Q and Qdot themselves are continuous.
fn hermite(x: f64) -> [[f64; 3]; 2] {
    let a = x.abs();
    if a >= 1.0 {
        return [[0.0; 3]; 2];
    }
    let b = 1.0 - a;
    [
        [
            b * b * (1.0 + 2.0 * a),
            x.signum() * (-6.0 * a * b),
            -6.0 + 12.0 * a,
        ],
        [
            x * b * b,
            b * (1.0 - 3.0 * a),
            x.signum() * (-4.0 + 6.0 * a),
        ],
    ]
}

struct TapKernel<'a> {
    geometry: &'a PointGeometry,
    axes: [[f64; 3]; 3],
    sign: f64,
    gamma: f64,
    lag: usize,
    pressure_scale: f64,
    velocity_scale: f64,
    base: [f64; 8],
    boundaries: [f64; 6],
    count: usize,
}
impl<'a> TapKernel<'a> {
    fn new(
        geometry: &'a PointGeometry,
        axes: [[f64; 3]; 3],
        sign: f64,
        acoustic: FluidAcousticState,
        dt: f64,
        lag: usize,
    ) -> Result<Self, AeroacError> {
        let l = geometry.distance;
        let step = acoustic.sound_speed_m_s() * dt;
        let gamma = l / step;
        let pressure_scale = acoustic.density_kg_m3() / l / dt / (4.0 * PI * PI);
        let velocity_scale = (1.0 / l) / l / (4.0 * PI * PI);
        if [gamma, pressure_scale, velocity_scale]
            .iter()
            .any(|v| !v.is_finite() || *v <= 0.0)
        {
            return Err(invalid("causal diffraction basis scale is unrepresentable"));
        }
        let at = |radius: f64, direction: [f64; 3]| {
            let x = (Dd::from_f64(lag as f64) - Dd::from_f64(radius) / Dd::from_f64(step)).to_f64();
            let shape = hermite(x);
            let mut result = [0.0; 8];
            for j in 0..2 {
                result[j] = acoustic.density_kg_m3() / radius / dt * shape[j][1] / (4.0 * PI);
                let radial =
                    (shape[j][0] + radius / step * shape[j][1]) / radius / radius / (4.0 * PI);
                for i in 0..3 {
                    result[2 + 3 * j + i] = radial * direction[i];
                }
            }
            result
        };
        let edge = at(l, geometry.edge_gradient);
        let mut base = [Dd::ZERO; 8];
        for j in 0..2 {
            let factor = if j == 0 { 1.0 } else { sign };
            let mut term = edge.map(|v| Dd::from_f64(v) * Dd::from_f64(0.5));
            if geometry.q[j] > 0.0 {
                let d = geometry.direct_distance[j];
                let direct = at(d, geometry.direct_delta[j].map(|v| v / d));
                for i in 0..8 {
                    term[i] = term[i] + Dd::from_f64(direct[i]) - Dd::from_f64(edge[i]);
                }
            }
            for i in 0..8 {
                base[i] = base[i] + term[i] * Dd::from_f64(factor);
            }
        }
        let mut base = base.map(Dd::to_f64);
        for j in 0..2 {
            let local = [base[2 + 3 * j], base[3 + 3 * j], base[4 + 3 * j]];
            let world = rotate_real(axes, local);
            base[2 + 3 * j..5 + 3 * j].copy_from_slice(&world);
        }
        // Only delays inside this compact temporal basis can contribute.
        // sqrt((R-L)*(R+L))/L retains the narrow first-arrival interval.
        let to_w = |multiple: f64| {
            let ratio = multiple / gamma;
            if ratio <= 1.0 {
                0.0
            } else {
                det::sqrt(ratio - 1.0) * det::sqrt(ratio + 1.0)
            }
        };
        let lower = to_w((lag as f64 - 1.0).max(0.0));
        let upper = to_w(lag as f64 + 1.0);
        let mut boundaries = [0.0; 6];
        let mut count = 0;
        if upper > lower {
            let mut candidates = [
                lower,
                upper,
                to_w(lag as f64),
                (geometry.q[0] / l).abs(),
                (geometry.q[1] / l).abs(),
            ];
            candidates.sort_by(f64::total_cmp);
            for w in candidates {
                if w >= lower && w <= upper && (count == 0 || w > boundaries[count - 1]) {
                    boundaries[count] = w;
                    count += 1;
                }
            }
        }
        if boundaries[..count]
            .iter()
            .chain(&base)
            .any(|v| !v.is_finite())
        {
            return Err(invalid(
                "causal diffraction interval/field is unrepresentable",
            ));
        }
        Ok(Self {
            geometry,
            axes,
            sign,
            gamma,
            lag,
            pressure_scale,
            velocity_scale,
            base,
            boundaries,
            count,
        })
    }
}
impl Spectrum<8> for TapKernel<'_> {
    fn regions(&self) -> usize {
        self.count.saturating_sub(1)
    }
    fn extent(&self) -> (f64, f64) {
        (0.0, 0.0)
    }
    fn resolved(&self, _region: usize, _lower: f64, _upper: f64) -> bool {
        true
    }
    fn sample(&self, region: usize, x: f64) -> Result<[C64; 8], AeroacError> {
        let low = self.boundaries[region];
        let width = self.boundaries[region + 1] - low;
        let w = (Dd::from_f64(low) + Dd::from_f64(width) * Dd::from_f64(x)).to_f64();
        let radius = det::hypot(1.0, w);
        let coordinate = (Dd::from_f64(self.lag as f64)
            - Dd::from_f64(self.gamma) * Dd::from_f64(radius))
        .to_f64();
        let shape = hermite(coordinate);
        let inv = 1.0 / radius;
        let inv3 = inv * inv * inv;
        let inv5 = inv3 * inv * inv;
        let mut values = [Dd::ZERO; 8];
        for j in 0..2 {
            let factor = if j == 0 { 1.0 } else { self.sign };
            let eta = self.geometry.q[j] / self.geometry.distance;
            let angle = det::atan2(eta, w);
            let scale = eta.abs().max(w);
            let quotient = if scale == 0.0 {
                1.0
            } else {
                let a = eta / scale;
                let b = w / scale;
                b * b / (a * a + b * b)
            };
            for basis in 0..2 {
                let [value, first, second] = shape[basis];
                let a = value + self.gamma * radius * first;
                let p =
                    self.pressure_scale * w * inv3 * (first + self.gamma * radius * second) * angle;
                values[basis] = values[basis] + Dd::from_f64(p * factor);
                let radial = w * (self.gamma * self.gamma * second * inv3 + 3.0 * a * inv5) * angle;
                let angular = -a * inv3 * quotient;
                for axis in 0..3 {
                    let v = self.velocity_scale
                        * (radial * self.geometry.edge_gradient[axis]
                            + angular * self.geometry.q_gradient[j][axis]);
                    let at = 2 + 3 * basis + axis;
                    values[at] = values[at] + Dd::from_f64(v * factor);
                }
            }
        }
        let mut values = values.map(Dd::to_f64);
        for basis in 0..2 {
            let v = rotate_real(
                self.axes,
                [
                    values[2 + 3 * basis],
                    values[3 + 3 * basis],
                    values[4 + 3 * basis],
                ],
            );
            values[2 + 3 * basis..5 + 3 * basis].copy_from_slice(&v);
        }
        let result = core::array::from_fn(|i| {
            C64::from_re(
                (Dd::from_f64(values[i]) * Dd::from_f64(width)
                    + Dd::from_f64(self.base[i]) / Dd::from_f64(self.regions() as f64))
                .to_f64(),
            )
        });
        if result.iter().any(|v| !v.re.is_finite()) {
            return Err(invalid(
                "causal diffraction quadrature sample is unrepresentable",
            ));
        }
        Ok(result)
    }
}
fn rotate_real(axes: [[f64; 3]; 3], vector: [f64; 3]) -> [f64; 3] {
    core::array::from_fn(|i| {
        (0..3)
            .fold(Dd::ZERO, |sum, j| {
                sum + Dd::from_f64(axes[j][i]) * Dd::from_f64(vector[j])
            })
            .to_f64()
    })
}

#[cfg(test)]
mod tests;
