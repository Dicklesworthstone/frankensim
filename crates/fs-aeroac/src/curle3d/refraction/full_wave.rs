//! Spherical-wave transmission through one original stationary fluid plane.
//!
//! This integrates the outgoing Sommerfeld spectrum, including its two critical
//! circles and evanescent tail. For clearances a,b, transverse range R, original
//! densities rho1,rho2, and `zeta_i = sqrt(k_i²-q²)` with nonnegative real/imaginary
//! parts, put `D=rho2*zeta1+rho1*zeta2` and `E=exp(i*zeta1*a+i*zeta2*b)`. Then
//! `p/Q = omega*rho1*rho2/(2*pi) integral q*J0(qR)*E/D dq`.
//! Differentiating this SAME field gives radial and normal particle velocity;
//! neither is replaced by pressure divided by a scalar acoustic impedance.
//!
//! The scalar outgoing identity and branch convention are given by W. C. Chew,
//! Purdue ECE 604, lecture 35, equations (35.1.11)--(35.1.14):
//! <https://engineering.purdue.edu/wcchew/ece604s20/Lecture%20Notes/Lect35.pdf>.
//! The original pressure/normal-velocity interface conditions supply
//! `Tp=2*rho2*zeta1/D`, cancelling the apparent incident grazing singularity.
//! All phasors use exp(-i omega t), with Q the prescribed source volume flow.
//!
//! Three branch-regularizing substitutions feed bounded adaptive Gauss--Kronrod
//! quadrature. Its error estimates are NOT certified enclosures. The omitted
//! evanescent tail has a separate analytic bound, evaluated in floating point.
//! Geometry, Bessel and other floating-point errors are outside both quantities.
//! The existing Bessel owner documents <5e-13 error relative to its envelope
//! through argument 3000, with increasing phase error at much larger arguments;
//! diagnostics retain the largest argument so this limit remains visible.

use super::{PlanarFluidInterface, norm, poll, positive};
use crate::AeroacError;
use crate::bessel::{j0, j1};
use crate::curle3d::CoherentCompactField;
use crate::curle3d::statistics::{ComplexSum, finite, invalid, nonfinite};
use crate::observation::HarmonicFieldLine;
use core::cmp::Ordering;
use core::f64::consts::{FRAC_PI_2, PI, TAU};
use fs_math::{c64::C64, dd::Dd, det};
use std::collections::BinaryHeap;

/// Coherent and statistically specified arrays of physical Q/F sources.
pub mod array;
/// Spherical-wave transmission through one finite original fluid layer.
pub mod layer;
/// Reflected and coherent total fields in the original incident fluid.
pub mod reflection;
/// Spherical reflection from a passive, locally reacting impedance plane.
pub mod wall;
/// Original symmetric volume-stress sources through the refracting interface.
pub mod quadrupole;

mod force;
pub use force::{
    FULL_WAVE_FORCE_SCOPE, FullWaveForceDiagnostics, FullWaveForceLimits, PreparedFullWaveForce,
};

/// Hard bound on integrand evaluations, including discarded parent panels.
pub const MAX_FULL_WAVE_EVALUATIONS: usize = 1_048_576;

/// Applicability retained by every coherent field returned by this owner.
pub const FULL_WAVE_MONOPOLE_SCOPE: &str = "Prescribed stationary compact 3D volume source; \
outgoing spherical-wave transmission across one infinite planar boundary between two \
stationary homogeneous lossless fluids, including reactive and evanescent fields. \
Estimated finite quadrature error and a separate analytic tail bound; no interval \
certificate, flow, absorption, finite boundaries, source feedback or measured-SPL qualification.";

/// Explicit numerical accuracy and work controls; no high-frequency admission.
#[derive(Clone, Copy, Debug)]
pub struct FullWaveMonopoleLimits {
    /// Positive absolute tolerance on p/Q [Pa / (m³/s)].
    pub absolute_pressure_tolerance_pa_per_m3_s: f64,
    /// Positive absolute tolerance on the complex velocity-vector norm / Q [1/m²].
    pub absolute_velocity_tolerance_per_m2: f64,
    /// Relative tolerance in [0,1). Each acceptance target is the larger of its
    /// absolute tolerance and this fraction of the computed transfer magnitude.
    pub relative_tolerance: f64,
    /// Integrand evaluations, in 1..=MAX_FULL_WAVE_EVALUATIONS. Each embedded
    /// panel costs 15; splitting it costs 30 more. An exact equal-speed field
    /// needs no quadrature. Exhaustion returns no prepared observer.
    pub maximum_evaluations: usize,
    /// Bisection depth per transformed region, in 1..=48. Root panels have depth 0.
    pub maximum_depth: usize,
}

/// Numerical evidence for a unit volume-flow transfer, excluding roundoff.
/// Multiply absolute errors/bounds by |Q| for a particular source amplitude.
#[derive(Clone, Copy, Debug)]
pub struct FullWaveMonopoleDiagnostics {
    pub evaluations: usize,
    pub maximum_depth_reached: usize,
    pub retained_panels: usize,
    /// Upper transverse wavenumber included by quadrature [rad/m]. Zero for
    /// the exact equal-speed expression, which has no spectral truncation.
    pub transverse_wavenumber_cutoff_rad_m: f64,
    /// Largest q*R in the retained spectrum; see the Bessel accuracy scope.
    pub maximum_bessel_argument: f64,
    /// Sum of embedded-panel pressure error estimates; not a certified bound.
    pub pressure_quadrature_error_pa_per_m3_s: f64,
    /// Separate radial/normal estimates [1/m²], in the original plane frame.
    pub velocity_quadrature_error_per_m2: [f64; 2],
    /// Analytic bound on omitted pressure transfer, separate from quadrature.
    pub pressure_tail_bound_pa_per_m3_s: f64,
    /// Analytic radial/normal tail bounds [1/m²]. Radial is exactly zero on axis.
    pub velocity_tail_bound_per_m2: [f64; 2],
    /// True only when the ORIGINAL sound speeds are exactly equal.
    pub exact_equal_speed: bool,
}

/// A complete transmitted source-to-microphone field on the original source clock.
/// Preparation has bounded temporary storage/work; application allocates nothing
/// and reuses four complex coefficients without quadrature or transcendental work.
#[derive(Clone, Debug)]
pub struct PreparedFullWaveMonopole {
    interface: PlanarFluidInterface,
    frequency_hz: f64,
    source_m: [f64; 3],
    observer_m: [f64; 3],
    transfer: [C64; 4],
    diagnostics: FullWaveMonopoleDiagnostics,
}

impl PreparedFullWaveMonopole {
    /// Prepare one physical volume source and one observer strictly on opposite
    /// sides. Near-interface and critical contributions are retained; very small
    /// clearances or highly oscillatory fields may exhaust the explicit budget.
    /// No geometrical-acoustics fallback or artificial spectral damping is used.
    pub fn new(
        interface: PlanarFluidInterface,
        frequency_hz: f64,
        source_m: [f64; 3],
        observer_m: [f64; 3],
        limits: FullWaveMonopoleLimits,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, AeroacError> {
        poll(cancelled)?;
        admit(limits)?;
        let omega = TAU * frequency_hz;
        if !positive(omega) {
            return Err(invalid("full-wave monopole frequency"));
        }
        if source_m.iter().chain(&observer_m).any(|v| !v.is_finite()) {
            return Err(nonfinite("full-wave monopole positions"));
        }
        let a = -interface.distance(source_m);
        let b = interface.distance(observer_m);
        let h = (a + b).to_f64();
        if !positive(a.to_f64()) || !positive(b.to_f64()) || !positive(h) {
            return Err(invalid(
                "full-wave source and observer must lie strictly on their fluid sides",
            ));
        }
        let lateral = core::array::from_fn(|i| {
            (Dd::from_f64(observer_m[i])
                - Dd::from_f64(source_m[i])
                - (a + b) * Dd::from_f64(interface.normal[i]))
            .to_f64()
        });
        let range = norm(lateral);
        if !range.is_finite() {
            return Err(nonfinite("full-wave lateral separation"));
        }
        let tangent = if range == 0.0 {
            [0.0; 3]
        } else {
            lateral.map(|v| v / range)
        };
        let c1 = interface.incident.sound_speed_m_s();
        let c2 = interface.receiving.sound_speed_m_s();
        let rho1 = interface.incident.density_kg_m3();
        let rho2 = interface.receiving.density_kg_m3();
        let density_scale = rho1.max(rho2);
        let density = [rho1 / density_scale, rho2 / density_scale];
        if density.iter().any(|v| !positive(*v)) {
            return Err(nonfinite("full-wave density contrast"));
        }
        let (field, diagnostics) = if c1 == c2 {
            let distance = det::hypot(range, h);
            let phase = Dd::from_f64(omega) * Dd::from_f64(distance) / Dd::from_f64(c1);
            if !positive(distance)
                || !positive(omega / c1)
                || !phase.hi.is_finite()
                || !phase.lo.is_finite()
            {
                return Err(nonfinite("equal-speed full-wave distance or phase"));
            }
            let (sin, cos) = det::sin_cos_dd(phase);
            let green = C64::new(cos.to_f64(), sin.to_f64()).scale((0.25 / PI) / distance);
            let rho_effective = rho1 * (density[1] / (density[0] + density[1])) * 2.0;
            let pressure = minus_i(green).scale(omega * rho_effective);
            // Differentiate the spherical Green function before dividing by
            // receiving density. This retains induction even when k*r is tiny.
            let radial = green
                * C64::new(1.0 / distance, -omega / c1)
                    .scale(2.0 * (density[0] / (density[0] + density[1])));
            (
                [
                    pressure,
                    radial.scale(range / distance),
                    radial.scale(h / distance),
                ],
                FullWaveMonopoleDiagnostics {
                    evaluations: 0,
                    maximum_depth_reached: 0,
                    retained_panels: 0,
                    transverse_wavenumber_cutoff_rad_m: 0.0,
                    maximum_bessel_argument: 0.0,
                    pressure_quadrature_error_pa_per_m3_s: 0.0,
                    velocity_quadrature_error_per_m2: [0.0; 2],
                    pressure_tail_bound_pa_per_m3_s: 0.0,
                    velocity_tail_bound_per_m2: [0.0; 2],
                    exact_equal_speed: true,
                },
            )
        } else {
            let k = [omega / c1, omega / c2];
            if k.iter().any(|v| !positive(*v)) || k[0] == k[1] {
                return Err(nonfinite("full-wave distinct fluid wavenumbers"));
            }
            let kmin = k[0].min(k[1]);
            let kmax = k[0].max(k[1]);
            // Avoid subtracting rounded squared wavenumbers near equal speeds.
            let gap = det::sqrt(kmax - kmin) * det::sqrt(kmax) * det::sqrt(1.0 + kmin / kmax);
            let pressure_scale = (omega / (2.0 * PI)) * (rho1 * density[1]);
            let velocity_scale = density[0] / (2.0 * PI);
            if !positive(gap) || !positive(pressure_scale) || !positive(velocity_scale) {
                return Err(nonfinite("full-wave spectral scales"));
            }
            let mut kernel = Kernel {
                a,
                b,
                range,
                density,
                kmin,
                kmax,
                gap,
                incident_has_larger_k: k[0] > k[1],
                pressure_scale,
                velocity_scale,
                tail_extent: 0.0,
            };
            let tail = kernel.prepare_tail(limits, cancelled)?;
            integrate(&kernel, tail, limits, cancelled)?
        };
        let transfer = [
            field[0],
            field[1].scale(tangent[0]) + field[2].scale(interface.normal[0]),
            field[1].scale(tangent[1]) + field[2].scale(interface.normal[1]),
            field[1].scale(tangent[2]) + field[2].scale(interface.normal[2]),
        ];
        if transfer.iter().any(|v| !finite(*v)) {
            return Err(nonfinite("full-wave field transfer"));
        }
        poll(cancelled)?;
        Ok(Self {
            interface,
            frequency_hz,
            source_m,
            observer_m,
            transfer,
            diagnostics,
        })
    }

    pub const fn diagnostics(&self) -> &FullWaveMonopoleDiagnostics {
        &self.diagnostics
    }
    pub const fn interface(&self) -> PlanarFluidInterface {
        self.interface
    }
    pub const fn frequency_hz(&self) -> f64 {
        self.frequency_hz
    }
    pub const fn source_position_m(&self) -> [f64; 3] {
        self.source_m
    }
    pub const fn observer_position_m(&self) -> [f64; 3] {
        self.observer_m
    }

    fn amplitudes(&self, volume_flow_m3_s: C64) -> Result<[C64; 4], AeroacError> {
        if !finite(volume_flow_m3_s) {
            return Err(nonfinite("full-wave volume-flow phasor"));
        }
        let field = self.transfer.map(|v| v * volume_flow_m3_s);
        if field.iter().any(|v| !finite(*v)) {
            return Err(nonfinite("full-wave field amplitude"));
        }
        Ok(field)
    }

    /// Coherent peak Q [m³/s] in exp(-i omega t) convention. The transfer already
    /// contains all propagation phase; active/reactive intensity uses p*conj(u)/2.
    pub fn evaluate(&self, volume_flow_m3_s: C64) -> Result<CoherentCompactField, AeroacError> {
        let [pressure, ux, uy, uz] = self.amplitudes(volume_flow_m3_s)?;
        let velocity = [ux, uy, uz];
        let intensity = velocity.map(|v| pressure.scale(0.5) * v.conj());
        if intensity.iter().any(|v| !finite(*v)) {
            return Err(nonfinite("full-wave intensity"));
        }
        Ok(CoherentCompactField {
            pressure_pa: pressure,
            particle_velocity_m_s: velocity,
            active_intensity_w_m2: intensity.map(|v| v.re),
            reactive_intensity_w_m2: intensity.map(|v| v.im),
            scope: FULL_WAVE_MONOPOLE_SCOPE,
        })
    }

    /// Feed the original finite-window observer with no second propagation delay.
    pub fn evaluate_harmonic_line(
        &self,
        volume_flow_m3_s: C64,
    ) -> Result<HarmonicFieldLine<3>, AeroacError> {
        let [pressure, ux, uy, uz] = self.amplitudes(volume_flow_m3_s)?;
        Ok(HarmonicFieldLine {
            frequency_hz: self.frequency_hz,
            pressure_pa: pressure,
            particle_velocity_m_s: [ux, uy, uz],
        })
    }
}

fn admit(limits: FullWaveMonopoleLimits) -> Result<(), AeroacError> {
    QuadratureLimits::from(limits).admit()
}

fn minus_i(z: C64) -> C64 {
    C64::new(z.im, -z.re)
}

struct Kernel {
    a: Dd,
    b: Dd,
    range: f64,
    density: [f64; 2],
    kmin: f64,
    kmax: f64,
    gap: f64,
    incident_has_larger_k: bool,
    pressure_scale: f64,
    velocity_scale: f64,
    tail_extent: f64,
}

#[derive(Clone, Copy)]
struct SpectralPoint {
    q: f64,
    zeta: [C64; 2],
    measure: [f64; 2],
}

// The original two critical circles and tail share this branch map across
// transmission and reflection. The measure factors multiply to q*dq/dx.
fn spectral_point(
    kmin: f64,
    kmax: f64,
    gap: f64,
    incident_has_larger_k: bool,
    tail_extent: f64,
    region: usize,
    x: f64,
) -> SpectralPoint {
    let (q, low, high, measure) = if region < 2 {
        let (sin, cos) = det::sin_cos(FRAC_PI_2 * x);
        if region == 0 {
            let q = kmin * sin;
            let zlow = kmin * cos;
            (
                q,
                C64::from_re(zlow),
                C64::from_re(det::hypot(gap, zlow)),
                [q, zlow * FRAC_PI_2],
            )
        } else {
            let ilo = gap * sin;
            let rhi = gap * cos;
            (
                det::hypot(kmin, ilo),
                C64::new(0.0, ilo),
                C64::from_re(rhi),
                [ilo, rhi * FRAC_PI_2],
            )
        }
    } else {
        let t = tail_extent * x;
        (
            det::hypot(kmax, t),
            C64::new(0.0, det::hypot(gap, t)),
            C64::new(0.0, t),
            [t, tail_extent],
        )
    };
    let zeta = if incident_has_larger_k {
        [high, low]
    } else {
        [low, high]
    };
    SpectralPoint { q, zeta, measure }
}

impl Kernel {
    fn point(&self, region: usize, x: f64) -> SpectralPoint {
        spectral_point(
            self.kmin,
            self.kmax,
            self.gap,
            self.incident_has_larger_k,
            self.tail_extent,
            region,
            x,
        )
    }

    fn sample(&self, region: usize, x: f64) -> Result<[C64; 3], AeroacError> {
        let SpectralPoint { q, zeta, measure } = self.point(region, x);
        let phase = self.a * Dd::from_f64(zeta[0].re) + self.b * Dd::from_f64(zeta[1].re);
        let attenuation =
            (self.a * Dd::from_f64(zeta[0].im) + self.b * Dd::from_f64(zeta[1].im)).to_f64();
        let argument = q * self.range;
        if !phase.hi.is_finite()
            || !phase.lo.is_finite()
            || !attenuation.is_finite()
            || !argument.is_finite()
        {
            return Err(nonfinite("full-wave spectral phase"));
        }
        let (sin, cos) = det::sin_cos_dd(phase);
        let phasor = C64::new(cos.to_f64(), sin.to_f64());
        let denominator = zeta[0].scale(self.density[1]) + zeta[1].scale(self.density[0]);
        let jzero = j0(argument);
        let jone = j1(argument);
        let values = if attenuation <= 700.0 {
            let exponential = phasor.scale(det::exp(-attenuation));
            let weight = exponential * (C64::from_re(measure[0]) / denominator).scale(measure[1]);
            let radial = weight.scale((self.velocity_scale * q) * jone);
            [
                weight.scale(self.pressure_scale * jzero),
                C64::new(-radial.im, radial.re),
                (weight * zeta[1]).scale(self.velocity_scale * jzero),
            ]
        } else {
            // Raw exp(-attenuation) can underflow before the original material,
            // frequency and spectral measure restore a representable field.
            // Normalize D before its magnitude: even |D| need not fit in f64.
            let scale = denominator.re.abs().max(denominator.im.abs());
            if !positive(scale) {
                return Err(nonfinite("strongly evanescent interface denominator"));
            }
            let normalized = C64::new(denominator.re / scale, denominator.im / scale);
            let magnitude = normalized.abs();
            let unit = phasor * normalized.conj().scale(1.0 / magnitude);
            let log_weight = det::ln(measure[0]) + det::ln(measure[1])
                - det::ln(scale)
                - det::ln(magnitude)
                - attenuation;
            let normal_magnitude = zeta[1].abs();
            let normal_unit = if normal_magnitude == 0.0 {
                C64::ZERO
            } else {
                C64::new(zeta[1].re / normal_magnitude, zeta[1].im / normal_magnitude)
            };
            [
                log_weighted(unit, log_weight + det::ln(self.pressure_scale), jzero),
                log_weighted(
                    C64::new(-unit.im, unit.re),
                    log_weight + det::ln(self.velocity_scale) + det::ln(q),
                    jone,
                ),
                log_weighted(
                    unit * normal_unit,
                    log_weight + det::ln(self.velocity_scale) + det::ln(normal_magnitude),
                    jzero,
                ),
            ]
        };
        if values.iter().any(|v| !finite(*v)) {
            return Err(nonfinite("full-wave spectral field"));
        }
        Ok(values)
    }

    // Each q, Re(zeta), Im(zeta) is monotone on its region. These endpoint
    // variations therefore bound total phase/envelope variation independently
    // of sampled GK agreement, preventing unresolved oscillations from passing.
    fn resolved(&self, region: usize, lower: f64, upper: f64) -> bool {
        let lo = self.point(region, lower);
        let hi = self.point(region, upper);
        let a = self.a.to_f64();
        let b = self.b.to_f64();
        let phase_span = (hi.q - lo.q).abs() * self.range
            + a * (hi.zeta[0].re - lo.zeta[0].re).abs()
            + b * (hi.zeta[1].re - lo.zeta[1].re).abs();
        let decay_span =
            a * (hi.zeta[0].im - lo.zeta[0].im).abs() + b * (hi.zeta[1].im - lo.zeta[1].im).abs();
        phase_span <= FRAC_PI_2 && decay_span <= 2.0
    }

    fn prepare_tail(
        &mut self,
        limits: FullWaveMonopoleLimits,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<[f64; 3], AeroacError> {
        let h = (self.a + self.b).to_f64();
        let hk = h * self.kmax;
        if !hk.is_finite() {
            return Err(nonfinite("full-wave clearance wavenumber"));
        }
        // For q²=kmax²+t²: q*dq=t*dt, |D|>=(rho1+rho2)*t,
        // |E|<=exp(-h*t), |J0|,|J1|<=1, q<=t+kmax.
        // Density has already been scaled out of D and both prefactors.
        let log_density = det::ln(self.density[0] + self.density[1]);
        let log_h = det::ln(h);
        let log_pressure = det::ln(self.pressure_scale) - log_density - log_h;
        let log_velocity = det::ln(self.velocity_scale) - log_density - 2.0 * log_h;
        let goals = [
            det::ln(limits.absolute_pressure_tolerance_pa_per_m3_s) - det::ln(4.0),
            det::ln(limits.absolute_velocity_tolerance_per_m2) - det::ln(4.0 * det::sqrt(2.0)),
        ];
        let logs = |s: f64| [log_pressure - s, log_velocity - s + det::ln(s + hk + 1.0)];
        let mut low = 0.0;
        let mut high = 1.0;
        let mut admitted = false;
        for _ in 0..32 {
            poll(cancelled)?;
            let bound = logs(high);
            if bound[0] <= goals[0] && bound[1] <= goals[1] {
                admitted = true;
                break;
            }
            low = high;
            high *= 2.0;
        }
        if !admitted {
            return Err(invalid(
                "full-wave evanescent tail cannot meet the requested tolerance",
            ));
        }
        for _ in 0..24 {
            poll(cancelled)?;
            let middle = low + (high - low) * 0.5;
            let bound = logs(middle);
            if bound[0] <= goals[0] && bound[1] <= goals[1] {
                high = middle;
            } else {
                low = middle;
            }
        }
        self.tail_extent = high / h;
        if !positive(self.tail_extent) {
            return Err(nonfinite("full-wave evanescent cutoff"));
        }
        let bound = logs(high).map(det::exp);
        Ok([
            bound[0],
            if self.range == 0.0 { 0.0 } else { bound[1] },
            bound[1],
        ])
    }
}

// Construct real and imaginary components separately, preserving exact zeros
// and signs without first materializing a potentially overflowing magnitude.
fn log_weighted(unit: C64, log_scale: f64, real_weight: f64) -> C64 {
    if real_weight == 0.0 {
        return C64::ZERO;
    }
    let logarithm = log_scale + det::ln(real_weight.abs());
    let component = |value: f64| {
        if value == 0.0 {
            0.0
        } else {
            det::exp(logarithm + det::ln(value.abs())) * value.signum() * real_weight.signum()
        }
    };
    C64::new(component(unit.re), component(unit.im))
}

/// Units belong to the public source-specific limits. The shared integration
/// engine only sees the pressure/velocity transfer norms and work budgets.
#[derive(Clone, Copy)]
struct QuadratureLimits {
    absolute_pressure: f64,
    absolute_velocity: f64,
    relative: f64,
    maximum_evaluations: usize,
    maximum_depth: usize,
}

impl From<FullWaveMonopoleLimits> for QuadratureLimits {
    fn from(limits: FullWaveMonopoleLimits) -> Self {
        Self {
            absolute_pressure: limits.absolute_pressure_tolerance_pa_per_m3_s,
            absolute_velocity: limits.absolute_velocity_tolerance_per_m2,
            relative: limits.relative_tolerance,
            maximum_evaluations: limits.maximum_evaluations,
            maximum_depth: limits.maximum_depth,
        }
    }
}

impl QuadratureLimits {
    fn admit(self) -> Result<(), AeroacError> {
        if !positive(self.absolute_pressure)
            || !positive(self.absolute_velocity)
            || !self.relative.is_finite()
            || !(0.0..1.0).contains(&self.relative)
            || !(1..=MAX_FULL_WAVE_EVALUATIONS).contains(&self.maximum_evaluations)
            || !(1..=48).contains(&self.maximum_depth)
        {
            return Err(invalid(
                "full-wave tolerances or explicit evaluation/depth budgets",
            ));
        }
        Ok(())
    }
}

/// Statically sized samples keep the original monopole's three-channel work
/// unchanged; force derivatives share the same nodes with seven channels.
trait Spectrum<const N: usize> {
    fn regions(&self) -> usize {
        3
    }
    fn extent(&self) -> (f64, f64);
    fn resolved(&self, region: usize, lower: f64, upper: f64) -> bool;
    fn sample(&self, region: usize, x: f64) -> Result<[C64; N], AeroacError>;
}

impl Spectrum<3> for Kernel {
    fn extent(&self) -> (f64, f64) {
        let cutoff = det::hypot(self.kmax, self.tail_extent);
        (cutoff, cutoff * self.range)
    }
    fn resolved(&self, region: usize, lower: f64, upper: f64) -> bool {
        Kernel::resolved(self, region, lower, upper)
    }
    fn sample(&self, region: usize, x: f64) -> Result<[C64; 3], AeroacError> {
        Kernel::sample(self, region, x)
    }
}

struct IntegratedSpectrum<const N: usize> {
    value: [C64; N],
    error: [f64; N],
    evaluations: usize,
    maximum_depth_reached: usize,
    retained_panels: usize,
    cutoff: f64,
    maximum_argument: f64,
}

fn real_norm(values: &[f64]) -> f64 {
    let Some((&first, remaining)) = values.split_first() else {
        return 0.0;
    };
    remaining
        .iter()
        .fold(first.abs(), |sum, v| det::hypot(sum, *v))
}

fn complex_norm(values: &[C64]) -> f64 {
    let Some((&first, remaining)) = values.split_first() else {
        return 0.0;
    };
    remaining
        .iter()
        .fold(first.abs(), |sum, v| det::hypot(sum, v.abs()))
}

// Embedded 7-point Gauss / 15-point Kronrod rule on [-1,1]. Shared nodes
// integrate every complex transfer together, with DD signed reductions.
const NODES: [f64; 8] = [
    0.991_455_371_120_812_6,
    0.949_107_912_342_758_5,
    0.864_864_423_359_769_1,
    0.741_531_185_599_394_5,
    0.586_087_235_467_691_1,
    0.405_845_151_377_397_2,
    0.207_784_955_007_898_5,
    0.0,
];
const KRONROD: [f64; 8] = [
    0.022_935_322_010_529_225,
    0.063_092_092_629_978_55,
    0.104_790_010_322_250_19,
    0.140_653_259_715_525_92,
    0.169_004_726_639_267_9,
    0.190_350_578_064_785_42,
    0.204_432_940_075_298_89,
    0.209_482_141_084_727_83,
];
const GAUSS: [f64; 4] = [
    0.129_484_966_168_869_7,
    0.279_705_391_489_276_64,
    0.381_830_050_505_118_9,
    0.417_959_183_673_469_4,
];

struct Panel<const N: usize> {
    region: usize,
    lower: f64,
    upper: f64,
    depth: usize,
    value: [C64; N],
    error: [f64; N],
    unresolved: bool,
    priority: f64,
    serial: usize,
}
impl<const N: usize> PartialEq for Panel<N> {
    fn eq(&self, rhs: &Self) -> bool {
        self.serial == rhs.serial
    }
}
impl<const N: usize> Eq for Panel<N> {}
impl<const N: usize> PartialOrd for Panel<N> {
    fn partial_cmp(&self, rhs: &Self) -> Option<Ordering> {
        Some(self.cmp(rhs))
    }
}
impl<const N: usize> Ord for Panel<N> {
    fn cmp(&self, rhs: &Self) -> Ordering {
        self.unresolved
            .cmp(&rhs.unresolved)
            .then_with(|| self.priority.total_cmp(&rhs.priority))
            .then_with(|| rhs.serial.cmp(&self.serial))
    }
}

struct Work {
    evaluations: usize,
    maximum_depth_reached: usize,
    serial: usize,
}

impl Work {
    #[allow(clippy::too_many_arguments)]
    fn panel<const N: usize, const PRESSURE: usize>(
        &mut self,
        kernel: &impl Spectrum<N>,
        region: usize,
        lower: f64,
        upper: f64,
        depth: usize,
        limits: QuadratureLimits,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Panel<N>, AeroacError> {
        if self.evaluations + 15 > limits.maximum_evaluations {
            return Err(invalid("full-wave quadrature evaluation budget exhausted"));
        }
        self.maximum_depth_reached = self.maximum_depth_reached.max(depth);
        let midpoint = lower + (upper - lower) * 0.5;
        let half_width = (upper - lower) * 0.5;
        let mut samples = [[C64::ZERO; N]; 15];
        let mut kronrod = [ComplexSum::default(); N];
        let mut gauss = [ComplexSum::default(); N];
        let mut absolute = [Dd::ZERO; N];
        for node in 0..15 {
            poll(cancelled)?;
            let index = node / 2;
            let x = if node == 14 {
                midpoint
            } else {
                midpoint
                    + if node % 2 == 0 {
                        -half_width * NODES[index]
                    } else {
                        half_width * NODES[index]
                    }
            };
            let values = kernel.sample(region, x)?;
            self.evaluations += 1;
            samples[node] = values;
            for component in 0..N {
                kronrod[component].add(values[component].scale(KRONROD[index]))?;
                absolute[component] =
                    absolute[component] + Dd::from_f64(values[component].abs() * KRONROD[index]);
                if index % 2 == 1 {
                    gauss[component].add(values[component].scale(GAUSS[index / 2]))?;
                }
            }
        }
        let result = kronrod.map(ComplexSum::value);
        let gauss = gauss.map(ComplexSum::value);
        let mut error = [0.0; N];
        for component in 0..N {
            let mean = result[component].scale(0.5);
            let mut deviation = Dd::ZERO;
            for (node, values) in samples.iter().enumerate() {
                deviation =
                    deviation + Dd::from_f64((values[component] - mean).abs() * KRONROD[node / 2]);
            }
            let resasc = deviation.to_f64() * half_width;
            let difference = (result[component] - gauss[component]).abs() * half_width;
            let estimated = if resasc > 0.0 && difference > 0.0 {
                let ratio = (200.0 * (difference / resasc)).min(1.0);
                resasc * ratio * det::sqrt(ratio)
            } else {
                difference
            };
            error[component] =
                estimated.max(50.0 * f64::EPSILON * absolute[component].to_f64() * half_width);
        }
        let value = result.map(|v| v.scale(half_width));
        if value.iter().any(|v| !finite(*v)) || error.iter().any(|v| !v.is_finite() || *v < 0.0) {
            return Err(nonfinite("full-wave quadrature panel"));
        }
        let priority = (real_norm(&error[..PRESSURE]) / limits.absolute_pressure)
            .max(real_norm(&error[PRESSURE..]) / limits.absolute_velocity);
        let panel = Panel {
            region,
            lower,
            upper,
            depth,
            value,
            error,
            unresolved: !kernel.resolved(region, lower, upper),
            priority,
            serial: self.serial,
        };
        self.serial += 1;
        Ok(panel)
    }
}

fn integrate(
    kernel: &Kernel,
    tail: [f64; 3],
    limits: FullWaveMonopoleLimits,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<([C64; 3], FullWaveMonopoleDiagnostics), AeroacError> {
    let result = integrate_spectrum::<3, 1>(kernel, tail, limits.into(), cancelled)?;
    Ok((
        result.value,
        FullWaveMonopoleDiagnostics {
            evaluations: result.evaluations,
            maximum_depth_reached: result.maximum_depth_reached,
            retained_panels: result.retained_panels,
            transverse_wavenumber_cutoff_rad_m: result.cutoff,
            maximum_bessel_argument: result.maximum_argument,
            pressure_quadrature_error_pa_per_m3_s: result.error[0],
            velocity_quadrature_error_per_m2: [result.error[1], result.error[2]],
            pressure_tail_bound_pa_per_m3_s: tail[0],
            velocity_tail_bound_per_m2: [tail[1], tail[2]],
            exact_equal_speed: false,
        },
    ))
}

fn integrate_spectrum<const N: usize, const PRESSURE: usize>(
    kernel: &impl Spectrum<N>,
    tail: [f64; N],
    limits: QuadratureLimits,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<IntegratedSpectrum<N>, AeroacError> {
    let (cutoff, maximum_argument) = kernel.extent();
    if !cutoff.is_finite() || !maximum_argument.is_finite() {
        return Err(nonfinite("full-wave retained spectral extent"));
    }
    let mut panels = BinaryHeap::new();
    // Each transformed-region root costs 15 evaluations; each split adds one
    // retained panel and costs another 30. No allocation occurs while sampling.
    panels
        .try_reserve_exact(limits.maximum_evaluations / 30 + kernel.regions())
        .map_err(|_| invalid("full-wave quadrature allocation exceeds available memory"))?;
    let mut sums = [ComplexSum::default(); N];
    let mut errors = [Dd::ZERO; N];
    let mut unresolved = 0_usize;
    let mut work = Work {
        evaluations: 0,
        maximum_depth_reached: 0,
        serial: 0,
    };
    for region in 0..kernel.regions() {
        let panel = work.panel::<N, PRESSURE>(kernel, region, 0.0, 1.0, 0, limits, cancelled)?;
        for component in 0..N {
            sums[component].add(panel.value[component])?;
            errors[component] = errors[component] + Dd::from_f64(panel.error[component]);
        }
        unresolved += usize::from(panel.unresolved);
        panels.push(panel);
    }
    loop {
        poll(cancelled)?;
        let value = sums.map(ComplexSum::value);
        let error = errors.map(Dd::to_f64);
        if value.iter().any(|v| !finite(*v)) || error.iter().any(|v| !v.is_finite() || *v < 0.0) {
            return Err(nonfinite("full-wave accumulated field or error estimates"));
        }
        let pressure_target = limits
            .absolute_pressure
            .max(limits.relative * complex_norm(&value[..PRESSURE]));
        let velocity_target = limits
            .absolute_velocity
            .max(limits.relative * complex_norm(&value[PRESSURE..]));
        let total_error = core::array::from_fn::<_, N, _>(|i| error[i] + tail[i]);
        let pressure_error = real_norm(&total_error[..PRESSURE]);
        let velocity_error = real_norm(&total_error[PRESSURE..]);
        if unresolved == 0 && pressure_error <= pressure_target && velocity_error <= velocity_target
        {
            return Ok(IntegratedSpectrum {
                value,
                error,
                evaluations: work.evaluations,
                maximum_depth_reached: work.maximum_depth_reached,
                retained_panels: panels.len(),
                cutoff,
                maximum_argument,
            });
        }
        let panel = panels
            .pop()
            .ok_or_else(|| invalid("full-wave quadrature lost all panels"))?;
        if panel.depth >= limits.maximum_depth {
            return Err(AeroacError::NotConverged {
                what: "full-wave quadrature depth or phase resolution",
                residual: (pressure_error / pressure_target).max(velocity_error / velocity_target),
            });
        }
        let midpoint = panel.lower + (panel.upper - panel.lower) * 0.5;
        if midpoint == panel.lower || midpoint == panel.upper {
            return Err(invalid(
                "full-wave quadrature interval is no longer representable",
            ));
        }
        // Construct both children before replacing a published parent sum.
        let left = work.panel::<N, PRESSURE>(
            kernel,
            panel.region,
            panel.lower,
            midpoint,
            panel.depth + 1,
            limits,
            cancelled,
        )?;
        let right = work.panel::<N, PRESSURE>(
            kernel,
            panel.region,
            midpoint,
            panel.upper,
            panel.depth + 1,
            limits,
            cancelled,
        )?;
        for component in 0..N {
            sums[component].add(-panel.value[component])?;
            sums[component].add(left.value[component])?;
            sums[component].add(right.value[component])?;
            errors[component] = errors[component] - Dd::from_f64(panel.error[component])
                + Dd::from_f64(left.error[component])
                + Dd::from_f64(right.error[component]);
        }
        unresolved = unresolved - usize::from(panel.unresolved)
            + usize::from(left.unresolved)
            + usize::from(right.unresolved);
        panels.push(left);
        panels.push(right);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn physical_scales_recover_fields_after_raw_evanescent_decay_underflows() {
        // The first case revives pressure through the original density; the
        // second revives all three components through frequency/spectral scale.
        // Evaluate the independent on-axis/in-plane tail algebra with two
        // balanced half-decays, so the oracle never needs log reconstruction.
        for (scale, density, pressure_scale, range) in [
            (1.0, [1.0, 1.0], 1e300, 0.0),
            (1e160, [1e-200, 1.0], 1e60, 1e-162),
        ] {
            let kernel = Kernel {
                a: Dd::from_f64(0.5 / scale),
                b: Dd::from_f64(0.5 / scale),
                range,
                density,
                kmin: scale,
                kmax: 2.0 * scale,
                gap: det::sqrt(3.0) * scale,
                incident_has_larger_k: true,
                pressure_scale,
                velocity_scale: density[0] / (2.0 * PI),
                tail_extent: 1000.0 * scale,
            };
            let point = kernel.point(2, 0.8);
            let attenuation = (kernel.a * Dd::from_f64(point.zeta[0].im)
                + kernel.b * Dd::from_f64(point.zeta[1].im))
            .to_f64();
            assert_eq!(det::exp(-attenuation), 0.0);
            let half_decay = det::exp(-0.5 * attenuation);
            let denominator = density[1] * point.zeta[0].im + density[0] * point.zeta[1].im;
            let measure = (point.measure[0] / denominator) * (point.measure[1] * half_decay);
            let argument = point.q * range;
            let expected = [
                C64::new(0.0, -measure * (pressure_scale * half_decay) * j0(argument)),
                C64::from_re(
                    measure * ((kernel.velocity_scale * point.q) * half_decay) * j1(argument),
                ),
                C64::from_re(
                    measure
                        * ((kernel.velocity_scale * point.zeta[1].im) * half_decay)
                        * j0(argument),
                ),
            ];
            let actual = kernel.sample(2, 0.8).unwrap();
            assert!(actual[0].abs() > 0.0);
            if range != 0.0 {
                assert!(actual[1].abs() > 0.0 && actual[2].abs() > 0.0);
            }
            for (actual, expected) in actual.into_iter().zip(expected) {
                if expected == C64::ZERO {
                    assert_eq!(actual, C64::ZERO);
                } else {
                    assert!(
                        (actual - expected).abs() <= 2e-12 * expected.abs(),
                        "actual={actual:?}, balanced-decay reference={expected:?}"
                    );
                }
            }
        }
    }
}
