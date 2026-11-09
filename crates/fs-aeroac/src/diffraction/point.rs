//! Spherical source diffraction by the original ideal straight edge.
//!
//! The exact half-plane Green function follows the Biot--Tolstoy impulse
//! solution, Ouis (2002), equations 3--6, doi:10.1006/jsvi.2000.3553.
//! With L²=(r+rs)²+dz², q=2 sqrt(r rs) cos((theta +/- theta_s)/2),
//! D²=L²-q², R²=L²+s² and f(s)=exp(i k R)/R, each direct/image term is
//! H(q)*exp(i k D)/(4 pi D) - q/(4 pi²) integral f(s)/(s²+q²) ds.
//! Integration by parts isolates the geometrical-boundary cancellation:
//! f(0)/(8 pi) + H(q)*(exp(i k D)/D-f(0))/(4 pi)
//! - integral f'(s)*atan(q/s) ds/(4 pi²).
//! The last integral runs along s=exp(i pi/4)*t, t>=0. No pole or square-root
//! cut is crossed; outgoing waves decay exponentially on this contour.
//! Differentiating this SAME integral supplies all three velocity components.
//! At q=0 the residual integral vanishes while its q derivative remains finite.
//!
//! This includes spherical spreading, near fields and both sheet faces without
//! a plane-wave, high-frequency, ray-amplitude or finite-edge approximation.
//! Embedded quadrature errors are estimates, not interval certificates.

use super::{HalfPlaneFrame, HalfPlaneReceiver, cis, finite, half_angle, invalid, poll};
use crate::AeroacError;
use crate::curle3d::refraction::full_wave::{
    QuadratureLimits, Spectrum, integrate_spectrum, log_weighted,
};
use crate::curle3d::{CoherentCompactField, IdealPlaneBoundary};
use crate::observation::HarmonicFieldLine;
use core::f64::consts::{FRAC_1_SQRT_2, FRAC_PI_2, PI, SQRT_2, TAU};
use fs_material::fluid_acoustics::FluidAcousticState;
use fs_math::{c64::C64, dd::Dd, det};

mod flow;

/// The same physical p/Q and u/Q tolerances and bounded adaptive work controls
/// used by the existing spherical-wave interface propagator.
pub use crate::curle3d::refraction::full_wave::FullWaveMonopoleLimits as HalfPlaneMonopoleLimits;

const SCOPE: &str = "Prescribed stationary point volume source diffracted by one \
infinite straight-edged rigid or pressure-release sheet in a stationary homogeneous \
lossless fluid. Full spherical pressure and analytic Euler velocity; estimated \
quadrature error, without spectral truncation. No flow, wall impedance, finite \
barrier, source feedback or measured-SPL qualification.";

/// Numerical estimates for the original unit volume-flow transfer.
/// Multiply absolute errors by |Q| for a source phasor. Geometry, material and
/// floating-point errors are outside the embedded quadrature estimates.
#[derive(Clone, Copy, Debug)]
pub struct HalfPlaneMonopoleDiagnostics {
    pub evaluations: usize,
    pub maximum_depth_reached: usize,
    pub retained_panels: usize,
    /// Shortest source-edge-receiver path [m].
    pub minimum_edge_path_m: f64,
    pub pressure_quadrature_error_pa_per_m3_s: f64,
    /// Cartesian world components, each in [1/m²].
    pub velocity_quadrature_error_per_m2: [f64; 3],
}

/// One original source, sheet, receiver, medium and frequency.
/// Preparation resolves the exact outgoing Green function; later source
/// amplitudes require only four complex products and no allocation.
#[derive(Clone, Copy, Debug)]
pub struct PreparedHalfPlaneMonopole {
    acoustic: FluidAcousticState,
    frequency_hz: f64,
    edge_flow_velocity_m_s: f64,
    frame: HalfPlaneFrame,
    boundary: IdealPlaneBoundary,
    source_position_m: [f64; 3],
    receiver: HalfPlaneReceiver,
    transfer: [C64; 4],
    diagnostics: HalfPlaneMonopoleDiagnostics,
}

impl PreparedHalfPlaneMonopole {
    /// Prepare a physical source volume-flow transfer with exp(-i omega t)
    /// peak phasors. Q is referenced to the source itself; all propagation
    /// phases are included once. The source must lie in the open fluid.
    /// Exact receiver sheet traces use HalfPlaneReceiver::SheetFace.
    ///
    /// # Errors
    /// Invalid frequency/geometry, an edge or source singularity, unresolved
    /// numerical work, or cancellation returns no prepared observation.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        acoustic: FluidAcousticState,
        frequency_hz: f64,
        frame: HalfPlaneFrame,
        boundary: IdealPlaneBoundary,
        source_position_m: [f64; 3],
        receiver: HalfPlaneReceiver,
        limits: HalfPlaneMonopoleLimits,
        cancelled: impl FnMut() -> bool,
    ) -> Result<Self, AeroacError> {
        Self::with_edge_parallel_flow(
            acoustic,
            frequency_hz,
            frame,
            boundary,
            source_position_m,
            receiver,
            0.0,
            limits,
            cancelled,
        )
    }

    /// A stationary physical volume source in uniform subsonic mean flow
    /// parallel to the sheet edge. The source and observation frequency remain
    /// on the laboratory clock; the supplied velocity is the FLUID velocity.
    ///
    /// The source enters continuity as Q*delta, with no additional force or
    /// vorticity injection. The original potential solves the convected wave
    /// equation; p=-rho*(partial_t+U*partial_z)*phi and u=grad(phi).
    /// A finite source contains all edge wave numbers, so its pressure and
    /// velocity require the complete Lorentz transform, not one intrinsic
    /// plane-wave frequency. The original material state is retained.
    ///
    /// # Errors
    /// Supersonic/nonfinite flow, invalid original source geometry, unresolved
    /// transformed coefficients/tolerances, exhausted work or cancellation.
    #[allow(clippy::too_many_arguments)]
    pub fn with_edge_parallel_flow(
        acoustic: FluidAcousticState,
        frequency_hz: f64,
        frame: HalfPlaneFrame,
        boundary: IdealPlaneBoundary,
        source_position_m: [f64; 3],
        receiver: HalfPlaneReceiver,
        edge_flow_velocity_m_s: f64,
        limits: HalfPlaneMonopoleLimits,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, AeroacError> {
        poll(&mut cancelled)?;
        let controls = QuadratureLimits::from(limits);
        controls.admit()?;
        let omega = TAU * frequency_hz;
        let k = omega / acoustic.sound_speed_m_s();
        let pressure_scale = omega * acoustic.density_kg_m3();
        if !frequency_hz.is_finite()
            || frequency_hz <= 0.0
            || !k.is_finite()
            || k <= 0.0
            || !pressure_scale.is_finite()
            || pressure_scale <= 0.0
        {
            return Err(invalid(
                "half-plane point source needs a resolved positive frequency",
            ));
        }
        let source = frame.local_receiver(HalfPlaneReceiver::WorldPosition(source_position_m))?;
        let point = frame.local_receiver(receiver)?;
        let (result, minimum_edge_path_m) = if edge_flow_velocity_m_s == 0.0 {
            // Keep the original zero-flow arithmetic and original world axes.
            let kernel = PointKernel::new(
                source,
                point,
                k,
                pressure_scale,
                boundary,
                frame.axes(),
                controls,
            )?;
            let result = integrate_spectrum::<4, 1>(&kernel, [0.0; 4], controls, &mut cancelled)?;
            (result, kernel.distance)
        } else {
            flow::prepare(
                acoustic,
                k,
                pressure_scale,
                source,
                point,
                boundary,
                frame.axes(),
                edge_flow_velocity_m_s,
                controls,
                &mut cancelled,
            )?
        };
        poll(&mut cancelled)?;
        Ok(Self {
            acoustic,
            frequency_hz,
            edge_flow_velocity_m_s,
            frame,
            boundary,
            source_position_m,
            receiver,
            transfer: result.value,
            diagnostics: HalfPlaneMonopoleDiagnostics {
                evaluations: result.evaluations,
                maximum_depth_reached: result.maximum_depth_reached,
                retained_panels: result.retained_panels,
                minimum_edge_path_m,
                pressure_quadrature_error_pa_per_m3_s: result.error[0],
                velocity_quadrature_error_per_m2: [
                    result.error[1],
                    result.error[2],
                    result.error[3],
                ],
            },
        })
    }

    pub const fn acoustic(&self) -> FluidAcousticState {
        self.acoustic
    }
    pub const fn frequency_hz(&self) -> f64 {
        self.frequency_hz
    }
    /// Original uniform fluid velocity along the placed edge [m/s].
    pub const fn edge_flow_velocity_m_s(&self) -> f64 {
        self.edge_flow_velocity_m_s
    }
    pub fn mean_flow_velocity_m_s(&self) -> [f64; 3] {
        self.frame.axes()[2].map(|axis| axis * self.edge_flow_velocity_m_s)
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
    pub const fn diagnostics(&self) -> &HalfPlaneMonopoleDiagnostics {
        &self.diagnostics
    }

    /// Coherent pressure and full Cartesian near/far particle velocity.
    pub fn evaluate_harmonic_line(
        &self,
        volume_flow_m3_s: C64,
    ) -> Result<HarmonicFieldLine<3>, AeroacError> {
        if !finite(volume_flow_m3_s) {
            return Err(invalid("half-plane source volume flow must be finite"));
        }
        let values = self.transfer.map(|v| v * volume_flow_m3_s);
        if values.iter().any(|&v| !finite(v)) {
            return Err(invalid("half-plane point-source field is unrepresentable"));
        }
        Ok(HarmonicFieldLine {
            frequency_hz: self.frequency_hz,
            pressure_pa: values[0],
            particle_velocity_m_s: [values[1], values[2], values[3]],
        })
    }

    /// Pressure, physical velocity, and signed active/reactive pressure work.
    pub fn evaluate(&self, volume_flow_m3_s: C64) -> Result<CoherentCompactField, AeroacError> {
        let line = self.evaluate_harmonic_line(volume_flow_m3_s)?;
        let flux = line
            .particle_velocity_m_s
            .map(|v| line.pressure_pa.scale(0.5) * v.conj());
        if flux.iter().any(|&v| !finite(v)) {
            return Err(invalid(
                "half-plane point-source intensity is unrepresentable",
            ));
        }
        Ok(CoherentCompactField {
            pressure_pa: line.pressure_pa,
            particle_velocity_m_s: line.particle_velocity_m_s,
            active_intensity_w_m2: flux.map(|v| v.re),
            reactive_intensity_w_m2: flux.map(|v| v.im),
            scope: if self.edge_flow_velocity_m_s == 0.0 {
                SCOPE
            } else {
                flow::SCOPE
            },
        })
    }
    /// Mean perturbation energy and laboratory transport from the same full
    /// coherent field. The flux includes mean-flow advection U*<E> in addition
    /// to the pressure work reported by evaluate().
    pub fn evaluate_energy_transport(
        &self,
        volume_flow_m3_s: C64,
    ) -> Result<crate::curle3d::energy::AcousticEnergyTransport, AeroacError> {
        super::energy_transport(
            self.acoustic,
            self.mean_flow_velocity_m_s(),
            self.evaluate_harmonic_line(volume_flow_m3_s)?,
        )
    }
}

struct PointKernel {
    distance: f64,
    kappa: Dd,
    phase: C64,
    eta: [f64; 2],
    q_gradient: [[f64; 3]; 2],
    edge_gradient: [f64; 3],
    axes: [[f64; 3]; 3],
    sign: f64,
    pressure_scale: f64,
    velocity_scale: f64,
    base: [C64; 4],
    breaks: [f64; 3],
    count: usize,
    scale: f64,
    absolute: [f64; 2],
}

impl PointKernel {
    #[allow(clippy::too_many_arguments)]
    fn new(
        source: [f64; 3],
        point: [f64; 3],
        k: f64,
        omega_rho: f64,
        boundary: IdealPlaneBoundary,
        axes: [[f64; 3]; 3],
        controls: QuadratureLimits,
    ) -> Result<Self, AeroacError> {
        let [xs, ys, zs] = source;
        let [x, y, z] = point;
        let rs = det::hypot(xs, ys);
        let r = det::hypot(x, y);
        let dz = (Dd::from_f64(z) - Dd::from_f64(zs)).to_f64();
        let radial_sum = r + rs;
        let distance = det::hypot(radial_sum, dz);
        if !r.is_finite()
            || r <= 0.0
            || !rs.is_finite()
            || rs <= 0.0
            || !distance.is_finite()
            || distance <= 0.0
        {
            return Err(invalid(
                "half-plane point source and receiver must resolve the edge",
            ));
        }
        let kappa = Dd::from_f64(k) * Dd::from_f64(distance);
        let velocity_scale = (1.0 / distance) / distance / (4.0 * PI * PI);
        let pressure_scale = omega_rho / distance / (4.0 * PI * PI);
        if !kappa.to_f64().is_finite()
            || kappa.to_f64() <= 0.0
            || !velocity_scale.is_finite()
            || velocity_scale <= 0.0
            || !pressure_scale.is_finite()
            || pressure_scale <= 0.0
        {
            return Err(invalid(
                "half-plane point-source transfer scale is unrepresentable",
            ));
        }
        let [ct, st] = half_angle(x / r, y / r);
        let [cs, ss] = half_angle(xs / rs, ys / rs);
        let cosine = [ct * cs + st * ss, ct * cs - st * ss];
        let sine = [st * cs - ct * ss, st * cs + ct * ss];
        let root = det::sqrt(r) * det::sqrt(rs);
        let q = cosine.map(|c| 2.0 * root * c);
        let eta = q.map(|q| q / distance);
        let radial = [x / r, y / r, 0.0];
        let tangential = [-y / r, x / r, 0.0];
        let ratio = det::sqrt(rs) / det::sqrt(r);
        let q_gradient = core::array::from_fn(|j| {
            core::array::from_fn(|i| ratio * (cosine[j] * radial[i] - sine[j] * tangential[i]))
        });
        let edge_gradient = [
            (radial_sum / distance) * radial[0],
            (radial_sum / distance) * radial[1],
            dz / distance,
        ];
        let phase = cis(kappa);
        let edge = phase.scale(1.0 / distance);
        let edge_derivative = phase * C64::new(-1.0 / distance, k).scale(1.0 / distance);
        let sign = boundary.sign();
        let mut green = C64::ZERO;
        let mut gradient = [C64::ZERO; 3];
        for j in 0..2 {
            let factor = if j == 0 { 1.0 } else { sign };
            let delta = [
                (Dd::from_f64(x) - Dd::from_f64(xs)).to_f64(),
                (Dd::from_f64(y) - Dd::from_f64(if j == 0 { ys } else { -ys })).to_f64(),
                dz,
            ];
            let d = det::hypot(det::hypot(delta[0], delta[1]), delta[2]);
            if !d.is_finite() || (q[j] > 0.0 && d <= 0.0) {
                return Err(invalid(
                    "half-plane receiver coincides with a visible point-source singularity",
                ));
            }
            // An occluded image may coincide with the receiver. Its Green
            // term is absent; the remaining diffraction field stays finite.
            let mut value = edge.scale(1.0 / (8.0 * PI));
            let mut derivative = edge_gradient.map(|g| edge_derivative.scale(g / (8.0 * PI)));
            if q[j] > 0.0 {
                // D-L = -q²/(D+L) retains the small transition difference.
                let difference = -(q[j] / (d + distance)) * q[j];
                let phase_difference = Dd::from_f64(k) * Dd::from_f64(difference);
                let half_phase = phase_difference * Dd::from_f64(0.5);
                let (s, _) = det::sin_cos_dd(half_phase);
                let em1 = C64::new(
                    (-Dd::from_f64(2.0) * s * s).to_f64(),
                    det::sin_cos_dd(phase_difference).0.to_f64(),
                );
                let diff =
                    phase * (em1.scale(1.0 / d) + C64::from_re((-difference / d) / distance));
                value = super::add(value, diff.scale(1.0 / (4.0 * PI)));
                let direct = cis(Dd::from_f64(k) * Dd::from_f64(d)).scale(1.0 / d);
                let direct_derivative = direct * C64::new(-1.0 / d, k);
                for i in 0..3 {
                    derivative[i] = super::add(
                        derivative[i],
                        (direct_derivative.scale(delta[i] / d)
                            - edge_derivative.scale(edge_gradient[i]))
                        .scale(1.0 / (4.0 * PI)),
                    );
                }
            }
            green = super::add(green, value.scale(factor));
            for i in 0..3 {
                gradient[i] = super::add(gradient[i], derivative[i].scale(factor));
            }
        }
        let world_gradient = rotate(axes, gradient);
        let base = [
            C64::new(green.im, -green.re).scale(omega_rho),
            -world_gradient[0],
            -world_gradient[1],
            -world_gradient[2],
        ];
        if base.iter().any(|&v| !finite(v))
            || q_gradient.iter().flatten().any(|v| !v.is_finite())
            || eta.iter().any(|v| !v.is_finite())
        {
            return Err(invalid(
                "half-plane source geometry or transfer is unrepresentable",
            ));
        }
        // Resolve both arctangent transition widths explicitly, including
        // arbitrarily narrow shadow boundaries. Zero-width regions are absent.
        let lo = eta[0].abs().min(eta[1].abs());
        let hi = eta[0].abs().max(eta[1].abs());
        let mut breaks = [0.0; 3];
        let mut count = 1;
        for value in [lo, hi] {
            if value > breaks[count - 1] {
                breaks[count] = value;
                count += 1;
            }
        }
        Ok(Self {
            distance,
            kappa,
            phase,
            eta,
            q_gradient,
            edge_gradient,
            axes,
            sign,
            pressure_scale,
            velocity_scale,
            base,
            breaks,
            count,
            scale: 1.0 / det::sqrt(kappa.to_f64().max(1.0)),
            absolute: [controls.absolute_pressure, controls.absolute_velocity],
        })
    }

    fn coordinate(&self, region: usize, x: f64) -> (f64, f64) {
        let lower = self.breaks[region];
        if region + 1 < self.count {
            let width = self.breaks[region + 1] - lower;
            (lower + width * x, width)
        } else {
            let rest = 1.0 - x;
            (lower + self.scale * x / rest, self.scale / rest / rest)
        }
    }

    // Conservative algebraic contour envelopes used only to decide whether
    // phase-resolution splitting is necessary. The tail is still integrated,
    // never dropped or reported as an interval certificate.
    fn tail_is_small(&self, w: f64) -> bool {
        // |sqrt(1+i*w²)| >= max(1,w). Integrating the unit interval
        // separately gives the same conservative constants below for w<1.
        let inverse = 1.0 / w.max(1.0);
        let k = self.kappa.to_f64();
        let attenuation = self.phase_extent(w).1;
        let common = k * inverse + inverse * inverse;
        let radial = k * k * inverse * inverse + k * inverse * inverse * inverse + inverse.powi(4);
        let mut pressure = 0.0;
        let mut velocity = 0.0;
        for j in 0..2 {
            let a = self.eta[j].abs();
            let q_norm = self.q_gradient[j]
                .iter()
                .fold(0.0, |n, &v| det::hypot(n, v));
            pressure += 4.0 * a * common;
            velocity += 8.0 * a * radial + 2.0 * common * q_norm;
        }
        let pressure_bound = log_weighted(
            C64::ONE,
            det::ln(self.pressure_scale) - attenuation,
            pressure,
        )
        .re;
        let velocity_bound = log_weighted(
            C64::ONE,
            det::ln(self.velocity_scale) - attenuation,
            velocity,
        )
        .re;
        pressure_bound <= self.absolute[0] * 0.01 && velocity_bound <= self.absolute[1] * 0.01
    }

    fn phase_extent(&self, w: f64) -> (f64, f64) {
        if !w.is_finite() {
            return (f64::INFINITY, f64::INFINITY);
        }
        let r = C64::new(1.0, w * w).sqrt();
        let difference = C64::new(0.0, w * w) / (r + C64::ONE);
        (
            self.kappa.to_f64() * difference.re,
            self.kappa.to_f64() * r.im,
        )
    }
}

impl Spectrum<4> for PointKernel {
    fn regions(&self) -> usize {
        self.count
    }
    fn extent(&self) -> (f64, f64) {
        (0.0, 0.0)
    }

    fn resolved(&self, region: usize, lower: f64, upper: f64) -> bool {
        let (a, _) = self.coordinate(region, lower);
        let (b, _) = self.coordinate(region, upper);
        if self.tail_is_small(a) {
            return true;
        }
        let (pa, _) = self.phase_extent(a);
        let (pb, _) = self.phase_extent(b);
        // Avoid accepting an aliased outgoing oscillation. The endpoint at
        // infinity is refined until its entire remaining tail is attenuated.
        pb - pa <= FRAC_PI_2
    }

    fn sample(&self, region: usize, x: f64) -> Result<[C64; 4], AeroacError> {
        let (w, jacobian) = self.coordinate(region, x);
        let r = C64::new(1.0, w * w).sqrt();
        let delta = C64::new(0.0, w * w) / (r + C64::ONE);
        let attenuation = (self.kappa * Dd::from_f64(r.im)).to_f64();
        if !w.is_finite() || !jacobian.is_finite() || !attenuation.is_finite() {
            return Err(invalid(
                "half-plane contour coordinate or attenuation is unrepresentable",
            ));
        }
        let logarithmic = attenuation > 700.0;
        let decay = if logarithmic {
            1.0
        } else {
            det::exp(-attenuation)
        };
        let base = self.base.map(|v| v.scale(1.0 / self.count as f64));
        let e = cis(self.kappa * Dd::from_f64(delta.re)).scale(decay);
        let kr = r.scale(self.kappa.to_f64());
        let ikr = C64::new(-kr.im, kr.re);
        let inverse = r.recip();
        let inverse3 = inverse * inverse * inverse;
        let g = e * (ikr - C64::ONE) * inverse3;
        let dg = e * (C64::from_re(3.0) - ikr.scale(3.0) - kr * kr) * inverse3 * inverse * inverse;
        let mut value = C64::ZERO;
        let mut gradient = [C64::ZERO; 3];
        for j in 0..2 {
            let factor = if j == 0 { 1.0 } else { self.sign };
            let a = self.eta[j];
            let arctangent = contour_arctangent(a, w);
            let angular_scale = a.abs().max(w);
            let ratio = if angular_scale == 0.0 {
                C64::ONE
            } else {
                let aw = a / angular_scale;
                let ww = w / angular_scale;
                C64::new(0.0, ww * ww) / C64::new(aw * aw, ww * ww)
            };
            let c = C64::new(0.0, w) * g * arctangent;
            let cr = C64::new(0.0, w) * dg * arctangent;
            let cq = C64::new(FRAC_1_SQRT_2, FRAC_1_SQRT_2) * g * ratio;
            value = super::add(value, c.scale(factor));
            for i in 0..3 {
                gradient[i] = super::add(
                    gradient[i],
                    (cr.scale(self.edge_gradient[i]) + cq.scale(self.q_gradient[j][i]))
                        .scale(factor),
                );
            }
        }
        if !finite(value) || gradient.iter().any(|&v| !finite(v)) {
            return Err(invalid(
                "half-plane undecayed contour field is unrepresentable",
            ));
        }
        // Reuse the original propagator's scale-aware exponential owner.
        // A tiny contour exponential can still carry resolved physical work.
        let scale = |unit: C64, physical: f64| {
            if logarithmic {
                log_weighted(
                    unit,
                    det::ln(physical) + det::ln(jacobian) - attenuation,
                    1.0,
                )
            } else {
                unit.scale(physical * jacobian)
            }
        };
        let pressure = scale(
            self.phase * C64::new(-value.im, value.re),
            self.pressure_scale,
        );
        let velocity =
            rotate(self.axes, gradient).map(|v| scale(self.phase * v, self.velocity_scale));
        let correction = [pressure, velocity[0], velocity[1], velocity[2]];
        let result = core::array::from_fn(|i| super::add(base[i], correction[i]));
        if result.iter().any(|&v| !finite(v)) {
            return Err(invalid("half-plane outgoing contour is unrepresentable"));
        }
        Ok(result)
    }
}

// atan(q/(exp(i pi/4)*w)), without a large q/w or squaring physical lengths.
// The logarithm ratio is bounded away from its cut. At small ratio its odd
// series avoids subtracting nearly equal logarithms.
fn contour_arctangent(q: f64, w: f64) -> C64 {
    let scale = q.abs().max(w);
    if scale == 0.0 {
        return C64::ZERO;
    }
    let a = q / scale;
    let b = w / scale;
    let cross = SQRT_2 * a * b;
    let sum = a * a + b * b;
    let z = cross / sum;
    let imaginary = if z.abs() < 0.05 {
        let square = z * z;
        -0.5 * z
            * (1.0
                + square
                    * (1.0 / 3.0
                        + square
                            * (1.0 / 5.0
                                + square * (1.0 / 7.0 + square * (1.0 / 9.0 + square / 11.0)))))
    } else {
        0.25 * det::ln((1.0 - z) / (1.0 + z))
    };
    C64::new(0.5 * det::atan2(cross, b * b - a * a), imaginary)
}

fn rotate(axes: [[f64; 3]; 3], vector: [C64; 3]) -> [C64; 3] {
    core::array::from_fn(|world| {
        let mut re = Dd::ZERO;
        let mut im = Dd::ZERO;
        for local in 0..3 {
            re = re + Dd::from_f64(axes[local][world]) * Dd::from_f64(vector[local].re);
            im = im + Dd::from_f64(axes[local][world]) * Dd::from_f64(vector[local].im);
        }
        C64::new(re.to_f64(), im.to_f64())
    })
}

#[cfg(test)]
mod tests;
