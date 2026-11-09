//! Causal physical force diffraction, including original force impulse.
//!
//! For a force ON the fluid, p = grad_source G[F] and
//! u = -grad_receiver grad_source G[integral(F dt)] / rho, away from the
//! source. G is the SAME real Biot--Tolstoy Green function as the volume
//! owner. Both derivatives precede the s=L*w substitution. No displaced
//! source pair, harmonic synthesis or pressure/impedance velocity is used.
//!
//! Original F and Fdot define a C1 Hermite force. Its primitive is integrated
//! analytically, including the nonzero impulse after emission. The mixed
//! derivative at q=0 retains the same half-weighted geometrical term as the
//! harmonic physical-force owner; logarithmic panels resolve narrow shadows.

use super::super::{HalfPlaneForceLimits, force::mixed_geometry};
use super::{
    AeroacError, C64, Convolution, Dd, FluidAcousticState, HalfPlaneFrame, HalfPlaneReceiver,
    HalfPlaneTransientField, IdealPlaneBoundary, MAX_SAMPLES, MAX_TOTAL_EVALUATIONS, PI,
    PointGeometry, QuadratureLimits, Spectrum, TimedCompactSource, det, hermite,
    integrate_spectrum, invalid, poll, rotate_real,
};

const SCOPE: &str = "Original stationary physical force ON fluid diffracted by one \
infinite ideal rigid or pressure-release half-plane in a homogeneous lossless fluid. \
C1 original F/Fdot interpolation, analytic force impulse, complete causal edge memory \
within the declared horizon, and full pressure/velocity. Quadrature estimates exclude \
interpolation, geometry and floating-point error. No source motion, mean flow, volume \
injection, finite barrier, absorption, feedback or measured-SPL qualification.";

/// Accuracy is applied to all three force directions and both Hermite bases.
#[derive(Clone, Copy, Debug)]
pub struct HalfPlaneForceTransientLimits {
    /// Per-tap physical pressure-vector and velocity-tensor tolerances.
    pub quadrature: HalfPlaneForceLimits,
    /// Shared finite preparation work over every lag and force direction.
    pub maximum_total_evaluations: usize,
    /// Peak numerical payload, including cached spectra and convolution work.
    pub maximum_payload_bytes: usize,
}

#[derive(Clone, Debug)]
pub struct HalfPlaneForceTransientDiagnostics {
    pub evaluations: usize,
    pub maximum_depth_reached: usize,
    pub maximum_retained_panels: usize,
    pub minimum_edge_delay_s: f64,
    pub minimum_arrival_delay_s: f64,
    pub fft_length: usize,
    pub peak_payload_bytes: usize,
    /// Sum of embedded tap estimates, in the same layout as coefficients().
    pub quadrature_error_l1: [f64; 24],
}

/// Original source position, receiver, material and finite uniform clock.
/// Prepared taps retain the ENTIRE force history, including its impulse.
#[derive(Clone, Debug)]
pub struct PreparedHalfPlaneForceTransient {
    acoustic: FluidAcousticState,
    frame: HalfPlaneFrame,
    boundary: IdealPlaneBoundary,
    source_position_m: [f64; 3],
    receiver: HalfPlaneReceiver,
    sample_interval_s: f64,
    coefficients: Vec<[f64; 24]>,
    errors: Vec<[f64; 24]>,
    convolution: Convolution<24>,
    diagnostics: HalfPlaneForceTransientDiagnostics,
}

impl PreparedHalfPlaneForceTransient {
    /// Bind one stationary source, ideal sheet and receiver to a finite uniform
    /// clock. Later records supply original force [N] and force rate [N/s].
    /// They must start quiescent, contain no volume injection, and retain zero
    /// force knots after emission for the desired tail observation interval.
    /// All three source directions and velocity components share quadrature.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        acoustic: FluidAcousticState,
        frame: HalfPlaneFrame,
        boundary: IdealPlaneBoundary,
        source_position_m: [f64; 3],
        receiver: HalfPlaneReceiver,
        sample_interval_s: f64,
        maximum_samples: usize,
        limits: HalfPlaneForceTransientLimits,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, AeroacError> {
        poll(&mut cancelled)?;
        let controls = QuadratureLimits::from(limits.quadrature);
        controls.admit()?;
        let dt = sample_interval_s;
        let step = acoustic.sound_speed_m_s() * dt;
        if !dt.is_finite()
            || dt <= 0.0
            || !step.is_finite()
            || step <= 0.0
            || !(1..=MAX_SAMPLES).contains(&maximum_samples)
            || !(dt * maximum_samples as f64).is_finite()
            || !(1..=MAX_TOTAL_EVALUATIONS).contains(&limits.maximum_total_evaluations)
        {
            return Err(invalid(
                "causal force requires a representable clock and bounded horizon/work",
            ));
        }
        let (fft_length, workspace) = Convolution::<24>::layout(maximum_samples)?;
        let peak = fft_length
            .checked_mul(31 * 16)
            .and_then(|v| {
                maximum_samples
                    .checked_mul(1300)
                    .and_then(|n| v.checked_add(n))
            })
            .and_then(|v| v.checked_add(workspace))
            .and_then(|v| {
                (controls.maximum_evaluations / 30 + 8)
                    .checked_mul(1536)
                    .and_then(|n| v.checked_add(n))
            })
            .and_then(|v| v.checked_add(8192))
            .ok_or_else(|| invalid("causal force payload overflow"))?;
        if peak > limits.maximum_payload_bytes {
            return Err(invalid(
                "causal force exceeds its declared numerical payload",
            ));
        }
        let source = frame.local_receiver(HalfPlaneReceiver::WorldPosition(source_position_m))?;
        let point = frame.local_receiver(receiver)?;
        let geometry = PointGeometry::new(source, point)?;
        let reciprocal = PointGeometry::new(point, source)?;
        let mixed = mixed_geometry(
            source,
            point,
            geometry.distance,
            geometry.edge_gradient,
            reciprocal.edge_gradient,
        );
        if mixed
            .0
            .iter()
            .flatten()
            .chain(mixed.1.iter().flatten().flatten())
            .any(|x| !x.is_finite())
        {
            return Err(invalid(
                "causal force geometry derivatives are unrepresentable",
            ));
        }
        let edge_delay = geometry.distance / acoustic.sound_speed_m_s();
        let mut arrival = edge_delay;
        for term in 0..2 {
            if geometry.q[term] > 0.0 {
                arrival = arrival.min(geometry.direct_distance[term] / acoustic.sound_speed_m_s());
            }
        }
        if !arrival.is_finite() || arrival <= 0.0 || !edge_delay.is_finite() {
            return Err(invalid("causal force arrival is unrepresentable"));
        }
        let mut coefficients = Vec::new();
        let mut errors = Vec::new();
        coefficients
            .try_reserve_exact(maximum_samples)
            .map_err(|_| invalid("causal force tap allocation failed"))?;
        errors
            .try_reserve_exact(maximum_samples)
            .map_err(|_| invalid("causal force estimate allocation failed"))?;
        let mut evaluations = 0;
        let mut depth = 0;
        let mut panels = 0;
        let mut error_l1 = [Dd::ZERO; 24];
        for lag in 0..maximum_samples {
            poll(&mut cancelled)?;
            let kernel = ForceTap::new(
                &geometry,
                &reciprocal,
                &mixed,
                frame.axes(),
                boundary.sign(),
                acoustic,
                dt,
                lag,
            )?;
            let (value, error) = if kernel.regions() == 0 {
                (kernel.base, [0.0; 24])
            } else {
                let remaining = limits
                    .maximum_total_evaluations
                    .checked_sub(evaluations)
                    .filter(|&v| v > 0)
                    .ok_or_else(|| invalid("causal force total quadrature work exhausted"))?;
                let result = integrate_spectrum::<24, 6>(
                    &kernel,
                    [0.0; 24],
                    QuadratureLimits {
                        maximum_evaluations: controls.maximum_evaluations.min(remaining),
                        ..controls
                    },
                    &mut cancelled,
                )?;
                evaluations += result.evaluations;
                depth = depth.max(result.maximum_depth_reached);
                panels = panels.max(result.retained_panels);
                (result.value.map(|v| v.re), result.error)
            };
            if value.iter().chain(&error).any(|v| !v.is_finite()) {
                return Err(invalid("causal force coefficients are unrepresentable"));
            }
            for j in 0..24 {
                error_l1[j] = error_l1[j] + Dd::from_f64(error[j]);
            }
            coefficients.push(value);
            errors.push(error);
        }
        let error_l1 = error_l1.map(Dd::to_f64);
        if error_l1.iter().any(|v| !v.is_finite()) {
            return Err(invalid("causal force estimate sum is unrepresentable"));
        }
        let convolution =
            Convolution::prepare(&coefficients, fft_length, workspace, &mut cancelled)?;
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
            diagnostics: HalfPlaneForceTransientDiagnostics {
                evaluations,
                maximum_depth_reached: depth,
                maximum_retained_panels: panels,
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
    pub fn diagnostics(&self) -> &HalfPlaneForceTransientDiagnostics {
        &self.diagnostics
    }
    /// Four output rows [p,ux,uy,uz], each with six columns
    /// [Fx,Fy,Fz,dt*Fxdot,dt*Fydot,dt*Fzdot], in original world axes.
    pub fn coefficients(&self) -> &[[f64; 24]] {
        &self.coefficients
    }

    /// Allocation-free evaluation of one accepted output knot.
    pub fn observe_at(
        &self,
        source: &[TimedCompactSource],
        index: usize,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<HalfPlaneTransientField, AeroacError> {
        self.admit(source, &mut cancelled)?;
        if index >= source.len() {
            return Err(invalid("causal force index is outside its record"));
        }
        let mut values = [Dd::ZERO; 4];
        let mut errors = [Dd::ZERO; 4];
        for lag in 0..=index {
            if lag % 64 == 0 {
                poll(&mut cancelled)?;
            }
            let sample = &source[index - lag];
            for output in 0..4 {
                for input in 0..6 {
                    let a = self.input(sample, input);
                    let k = 6 * output + input;
                    values[output] =
                        values[output] + Dd::from_f64(self.coefficients[lag][k]) * Dd::from_f64(a);
                    errors[output] =
                        errors[output] + Dd::from_f64(self.errors[lag][k]) * Dd::from_f64(a.abs());
                }
            }
        }
        poll(&mut cancelled)?;
        self.field(
            source[index].source_time_s,
            values.map(Dd::to_f64),
            errors.map(Dd::to_f64),
        )
    }
    /// Cached zero-padded linear FFT convolution of every original force channel.
    /// Supply zero force knots after emission to observe the complete edge tail
    /// and retained impulse within this horizon. No future knot enters an output.
    pub fn observe_record(
        &self,
        source: &[TimedCompactSource],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Vec<HalfPlaneTransientField>, AeroacError> {
        let amplitudes = self.admit(source, &mut cancelled)?;
        let values = self.convolution.apply_channels(
            source.len(),
            |i, j| self.input(&source[i], j),
            amplitudes,
            core::array::from_fn(|i| core::array::from_fn(|j| 6 * i + j)),
            &mut cancelled,
        )?;
        let error = core::array::from_fn(|i| {
            (0..6)
                .fold(Dd::ZERO, |s, j| {
                    s + Dd::from_f64(self.diagnostics.quadrature_error_l1[6 * i + j])
                        * Dd::from_f64(amplitudes[j])
                })
                .to_f64()
        });
        let support = source
            .iter()
            .position(|s| (0..6).any(|j| self.input(s, j) != 0.0))
            .map(|i| i.saturating_sub(1));
        let mut result = Vec::new();
        result
            .try_reserve_exact(source.len())
            .map_err(|_| invalid("causal force result allocation failed"))?;
        for (i, value) in values.into_iter().enumerate() {
            if i % 64 == 0 {
                poll(&mut cancelled)?;
            }
            let empty = support.is_none_or(|start| {
                ((Dd::from_f64(i as f64) - Dd::from_f64(start as f64))
                    * Dd::from_f64(self.sample_interval_s))
                .to_f64()
                    < self.diagnostics.minimum_arrival_delay_s
            });
            result.push(self.field(
                source[i].source_time_s,
                if empty { [0.0; 4] } else { value },
                error,
            )?);
        }
        poll(&mut cancelled)?;
        Ok(result)
    }
    fn input(&self, sample: &TimedCompactSource, channel: usize) -> f64 {
        if channel < 3 {
            sample.strength.force_on_fluid_n[channel]
        } else {
            self.sample_interval_s * sample.strength.force_rate_on_fluid_n_s[channel - 3]
        }
    }
    fn admit(
        &self,
        source: &[TimedCompactSource],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<[f64; 6], AeroacError> {
        poll(cancelled)?;
        if source.is_empty() || source.len() > self.maximum_samples() {
            return Err(invalid(
                "causal force record is empty or exceeds its complete-memory horizon",
            ));
        }
        let start = source[0].source_time_s;
        if !start.is_finite()
            || source[0]
                .strength
                .force_on_fluid_n
                .iter()
                .chain(&source[0].strength.force_rate_on_fluid_n_s)
                .any(|&v| v != 0.0)
        {
            return Err(invalid(
                "causal force requires its quiescent zero-strength start",
            ));
        }
        let mut amplitude = [0.0_f64; 6];
        for (i, sample) in source.iter().enumerate() {
            if i % 64 == 0 {
                poll(cancelled)?;
            }
            let expected = (Dd::from_f64(start)
                + Dd::from_f64(self.sample_interval_s) * Dd::from_f64(i as f64))
            .to_f64();
            if !sample.source_time_s.is_finite()
                || sample.source_time_s != expected
                || (i > 0 && sample.source_time_s <= source[i - 1].source_time_s)
                || sample.strength.volume_flow_m3_s != 0.0
                || sample.strength.volume_flow_rate_m3_s2 != 0.0
            {
                return Err(invalid(
                    "causal force requires pure original force on its prepared uniform clock",
                ));
            }
            for j in 0..6 {
                let value = self.input(sample, j);
                if !value.is_finite()
                    || (j >= 3
                        && (!sample.strength.force_rate_on_fluid_n_s[j - 3].is_finite()
                            || (sample.strength.force_rate_on_fluid_n_s[j - 3] != 0.0
                                && value == 0.0)))
                {
                    return Err(invalid("causal force or scaled rate is unrepresentable"));
                }
                amplitude[j] = amplitude[j].max(value.abs());
            }
        }
        Ok(amplitude)
    }
    fn field(
        &self,
        time: f64,
        value: [f64; 4],
        error: [f64; 4],
    ) -> Result<HalfPlaneTransientField, AeroacError> {
        let velocity = [value[1], value[2], value[3]];
        let energy =
            crate::curle3d::energy::instantaneous(self.acoustic, [0.0; 3], value[0], velocity)?;
        if error.iter().any(|v| !v.is_finite() || *v < 0.0) {
            return Err(invalid("causal force field estimate is unrepresentable"));
        }
        Ok(HalfPlaneTransientField {
            receiver_time_s: time,
            pressure_pa: value[0],
            particle_velocity_m_s: velocity,
            instantaneous_intensity_w_m2: energy.pressure_work_flux_w_m2,
            perturbation_energy_j_m3: energy.energy_density_j_m3,
            instantaneous_energy_flux_w_m2: energy.energy_flux_w_m2,
            pressure_quadrature_error_pa: error[0],
            velocity_quadrature_error_m_s: [error[1], error[2], error[3]],
            scope: SCOPE,
        })
    }
}

/// Integral from the quiescent past of the two ORIGINAL Hermite basis functions.
/// Beyond support the value basis retains unit impulse; the rate basis has zero.
fn primitive(x: f64) -> [f64; 2] {
    if x <= -1.0 {
        return [0.0; 2];
    }
    if x >= 1.0 {
        return [1.0, 0.0];
    }
    let a = 1.0 - x.abs();
    let a3 = a * a * a;
    let half = a3 * (1.0 - 0.5 * a);
    [
        if x < 0.0 { half } else { 1.0 - half },
        a3 * (0.25 * a - 1.0 / 3.0),
    ]
}

struct ForceTap<'a> {
    geometry: &'a PointGeometry,
    reciprocal: &'a PointGeometry,
    mixed: &'a ([[f64; 3]; 3], [[[f64; 3]; 3]; 2]),
    axes: [[f64; 3]; 3],
    sign: f64,
    gamma: f64,
    lag: usize,
    pressure_scale: f64,
    velocity_scale: f64,
    base: [f64; 24],
    boundaries: [f64; 6],
    count: usize,
}
impl<'a> ForceTap<'a> {
    #[allow(clippy::too_many_arguments)]
    fn new(
        geometry: &'a PointGeometry,
        reciprocal: &'a PointGeometry,
        mixed: &'a ([[f64; 3]; 3], [[[f64; 3]; 3]; 2]),
        axes: [[f64; 3]; 3],
        sign: f64,
        medium: FluidAcousticState,
        dt: f64,
        lag: usize,
    ) -> Result<Self, AeroacError> {
        let l = geometry.distance;
        let gamma = l / (medium.sound_speed_m_s() * dt);
        let pressure_scale = (1.0 / l) / l / (4.0 * PI * PI);
        let velocity_scale = pressure_scale / l * (dt / medium.density_kg_m3());
        if [gamma, pressure_scale, velocity_scale]
            .iter()
            .any(|v| !v.is_finite() || *v <= 0.0)
        {
            return Err(invalid("causal force basis scale is unrepresentable"));
        }
        let at = |radius: f64,
                  receiver_gradient: [f64; 3],
                  source_gradient: [f64; 3],
                  hessian: [[f64; 3]; 3]| {
            let z = radius / (medium.sound_speed_m_s() * dt);
            let x = (Dd::from_f64(lag as f64) - Dd::from_f64(z)).to_f64();
            let shape = hermite(x);
            let integral = primitive(x);
            let mut result = [0.0; 24];
            for basis in 0..2 {
                let [f, fd, _] = shape[basis];
                let a = integral[basis] + z * f;
                let b = 2.0 * a + z * z * fd;
                let p = -(f + z * fd) / radius / radius / (4.0 * PI);
                let v = dt / medium.density_kg_m3() / radius / radius / radius / (4.0 * PI);
                for j in 0..3 {
                    result[3 * basis + j] = p * source_gradient[j];
                    for i in 0..3 {
                        result[6 * (i + 1) + 3 * basis + j] = -v
                            * (b * receiver_gradient[i] * source_gradient[j] - a * hessian[i][j]);
                    }
                }
            }
            result
        };
        let edge = at(l, geometry.edge_gradient, reciprocal.edge_gradient, mixed.0);
        let mut base = [Dd::ZERO; 24];
        for term in 0..2 {
            let factor = if term == 0 { 1.0 } else { sign };
            let visible = if geometry.q[term] > 0.0 {
                1.0
            } else if geometry.q[term] < 0.0 {
                0.0
            } else {
                0.5
            };
            let mut addition = edge.map(|v| Dd::from_f64(0.5) * Dd::from_f64(v));
            if visible > 0.0 {
                let d = geometry.direct_distance[term];
                if d <= 0.0 {
                    return Err(invalid("causal force coincides with a visible source"));
                }
                let normal = geometry.direct_delta[term].map(|v| v / d);
                let reflection = [1.0, if term == 0 { 1.0 } else { -1.0 }, 1.0];
                let source_gradient = core::array::from_fn(|j| -reflection[j] * normal[j]);
                let hessian = core::array::from_fn(|i| {
                    core::array::from_fn(|j| {
                        -reflection[j] * ((if i == j { 1.0 } else { 0.0 }) - normal[i] * normal[j])
                    })
                });
                let direct = at(d, normal, source_gradient, hessian);
                for k in 0..24 {
                    addition[k] = addition[k]
                        + Dd::from_f64(visible) * (Dd::from_f64(direct[k]) - Dd::from_f64(edge[k]));
                }
            }
            for k in 0..24 {
                base[k] = base[k] + Dd::from_f64(factor) * addition[k];
            }
        }
        let base = rotate(axes, base.map(Dd::to_f64));
        let to_w = |multiple: f64| {
            let ratio = multiple / gamma;
            if ratio <= 1.0 {
                0.0
            } else {
                det::sqrt(ratio - 1.0) * det::sqrt(ratio + 1.0)
            }
        };
        let upper = to_w(lag as f64 + 1.0);
        let mut candidates = [
            0.0,
            upper,
            to_w((lag as f64 - 1.0).max(0.0)),
            to_w(lag as f64),
            (geometry.q[0] / l).abs(),
            (geometry.q[1] / l).abs(),
        ];
        candidates.sort_by(f64::total_cmp);
        let mut boundaries = [0.0; 6];
        let mut count = 0;
        if upper > 0.0 {
            for w in candidates {
                if w <= upper && (count == 0 || w > boundaries[count - 1]) {
                    boundaries[count] = w;
                    count += 1;
                }
            }
        }
        if base.iter().chain(&boundaries).any(|v| !v.is_finite()) {
            return Err(invalid(
                "causal force base or integration extent is unrepresentable",
            ));
        }
        Ok(Self {
            geometry,
            reciprocal,
            mixed,
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
impl Spectrum<24> for ForceTap<'_> {
    fn regions(&self) -> usize {
        self.count.saturating_sub(1)
    }
    fn extent(&self) -> (f64, f64) {
        (0.0, 0.0)
    }
    fn resolved(&self, _region: usize, _lower: f64, _upper: f64) -> bool {
        true
    }
    fn sample(&self, region: usize, x: f64) -> Result<[C64; 24], AeroacError> {
        let lo = self.boundaries[region];
        let hi = self.boundaries[region + 1];
        let (w, jac) = if lo == 0.0 || hi <= 2.0 * lo {
            let width = hi - lo;
            (
                (Dd::from_f64(lo) + Dd::from_f64(width) * Dd::from_f64(x)).to_f64(),
                width,
            )
        } else {
            let width = det::ln(hi) - det::ln(lo);
            let w = det::exp(det::ln(lo) + width * x);
            (w, width * w)
        };
        let radius = det::hypot(1.0, w);
        let z = self.gamma * radius;
        let coordinate = (Dd::from_f64(self.lag as f64)
            - Dd::from_f64(self.gamma) * Dd::from_f64(radius))
        .to_f64();
        let shape = hermite(coordinate);
        let integral = primitive(coordinate);
        let inv = 1.0 / radius;
        let inv3 = inv * inv * inv;
        let inv5 = inv3 * inv * inv;
        let inv7 = inv5 * inv * inv;
        let mut result = [Dd::ZERO; 24];
        for term in 0..2 {
            let sign = if term == 0 { 1.0 } else { self.sign };
            let eta = self.geometry.q[term] / self.geometry.distance;
            let angle = det::atan2(eta, w);
            let scale = eta.abs().max(w);
            let a = eta / scale;
            let b = w / scale;
            let den = a * a + b * b;
            let quotient = b * b / den;
            for basis in 0..2 {
                let [f, fd, fdd] = shape[basis];
                let p0 = f + z * fd;
                let pr = w * (3.0 * p0 + z * z * fdd) * inv5 * angle * jac;
                let pq = -p0 * inv3 * quotient * jac;
                let h = integral[basis] + z * f;
                let first = 3.0 * h + z * z * fd;
                let radial = w * first * inv5 * angle * jac;
                let second = w
                    * (first * inv5 - (15.0 * h + 6.0 * z * z * fd + z * z * z * fdd) * inv7)
                    * angle
                    * jac;
                let cross = first * inv5 * quotient * jac;
                // jac/scale keeps the finite integrated narrow-shadow weight.
                let angular2 = 2.0 * h * inv3 * quotient * (a / den) * (jac / scale);
                let angular = -h * inv3 * quotient * jac;
                for j in 0..3 {
                    let p = -self.pressure_scale
                        * (pr * self.reciprocal.edge_gradient[j]
                            + pq * self.reciprocal.q_gradient[term][j]);
                    result[3 * basis + j] = result[3 * basis + j] + Dd::from_f64(sign * p);
                    for i in 0..3 {
                        let v = self.velocity_scale
                            * (second
                                * self.geometry.edge_gradient[i]
                                * self.reciprocal.edge_gradient[j]
                                + radial * self.mixed.0[i][j]
                                + cross
                                    * (self.geometry.edge_gradient[i]
                                        * self.reciprocal.q_gradient[term][j]
                                        + self.geometry.q_gradient[term][i]
                                            * self.reciprocal.edge_gradient[j])
                                + angular2
                                    * self.geometry.q_gradient[term][i]
                                    * self.reciprocal.q_gradient[term][j]
                                + angular * self.mixed.1[term][i][j]);
                        let at = 6 * (i + 1) + 3 * basis + j;
                        result[at] = result[at] + Dd::from_f64(sign * v);
                    }
                }
            }
        }
        let world = rotate(self.axes, result.map(Dd::to_f64));
        let value = core::array::from_fn(|k| {
            C64::from_re(
                (Dd::from_f64(world[k])
                    + Dd::from_f64(self.base[k]) / Dd::from_f64(self.regions() as f64))
                .to_f64(),
            )
        });
        if value.iter().any(|v| !v.re.is_finite()) {
            return Err(invalid("causal force quadrature sample is unrepresentable"));
        }
        Ok(value)
    }
}
fn rotate(axes: [[f64; 3]; 3], local: [f64; 24]) -> [f64; 24] {
    let mut world = [0.0; 24];
    for basis in 0..2 {
        world[3 * basis..3 * basis + 3].copy_from_slice(&rotate_real(
            axes,
            [local[3 * basis], local[3 * basis + 1], local[3 * basis + 2]],
        ));
        let columns: [[f64; 3]; 3] = core::array::from_fn(|j| {
            rotate_real(
                axes,
                [
                    local[6 + 3 * basis + j],
                    local[12 + 3 * basis + j],
                    local[18 + 3 * basis + j],
                ],
            )
        });
        for i in 0..3 {
            let row = rotate_real(axes, [columns[0][i], columns[1][i], columns[2][i]]);
            world[6 * (i + 1) + 3 * basis..6 * (i + 1) + 3 * basis + 3].copy_from_slice(&row);
        }
    }
    world
}

#[cfg(test)]
mod tests;
