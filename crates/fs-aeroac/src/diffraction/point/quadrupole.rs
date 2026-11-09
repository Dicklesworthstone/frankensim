//! Physical symmetric volume-stress diffraction by the original half-plane.
//!
//! p=M_ij d_si d_sj G and u_l=d_xl p/(i omega rho), with the original
//! Lighthill stress moment [N m], including its trace and physical shear.
//! All 24 coefficients differentiate the parent's Biot--Tolstoy Green function
//! analytically BEFORE contour substitution; no displaced-source stencil.
//!
//! Third derivatives require a stronger shadow-boundary regularization. With
//! f(s)=exp(ik sqrt(L²+s²))/sqrt(L²+s²), A=f'(R=L)/L, subtract
//! r(s)=A*s/(1+s²/L²) from the integration-by-parts integrand. Its primitive is
//! integral r(s)*atan(q/s) ds = A*L²*pi/2*sign(q)*ln(1+|q|/L).
//! The explicit q|q| term cancels the geometrical wave's second-derivative
//! jump. The remaining integrand is O(s³)*atan(q/s), so three source/receiver
//! derivatives are regular even at q=0. The complete contour is integrated;
//! the subtraction is an exact algebraic identity, not a low-frequency model.
//! Green-function reference: Ouis (2002), doi:10.1006/jsvi.2000.3553.
//! Stress convention: Lighthill (1952), doi:10.1098/rspa.1952.0060.

use super::*;
use crate::curle3d::quadrupole::QuadrupolePhasor;
pub use crate::curle3d::refraction::full_wave::quadrupole::FullWaveQuadrupoleLimits as HalfPlaneQuadrupoleLimits;
mod derivatives;
use derivatives::Derivatives as Jet;

const SCOPE: &str = "Prescribed stationary symmetric Lighthill stress moment [N m] \
diffracted by one infinite ideal rigid or pressure-release straight-edged sheet \
in a stationary homogeneous lossless fluid; full source Hessian and analytic \
Euler velocity with original trace, shear and source clock. Quadrature errors \
are estimates. No finite source-size qualification, unresolved stress completion, \
finite barrier, absorption, flow, source feedback or measured-SPL qualification.";

/// Absolute errors for the physical Cartesian six-column map. Column order is
/// xx,yy,zz,xy,xz,yz; shear multiplicity is already included, as in the map.
/// For one source, sum each error times the corresponding |M_j|. Geometry and
/// floating-point errors are outside these embedded quadrature estimates.
#[derive(Clone, Copy, Debug)]
pub struct HalfPlaneQuadrupoleDiagnostics {
    pub evaluations: usize,
    pub maximum_depth_reached: usize,
    pub retained_panels: usize,
    pub minimum_edge_path_m: f64,
    pub pressure_quadrature_error_pa_per_n_m: [f64; 6],
    pub velocity_quadrature_error_m_s_per_n_m: [[f64; 6]; 3],
}

/// One original physical stress position, material, edge, receiver and clock.
/// A reusable 4-by-6 map retains every tensor cross term. Warm evaluation
/// requires 24 complex products and fixed-order DD sums, without allocation,
/// quadrature, transcendental functions or another propagation phase.
#[derive(Clone, Copy, Debug)]
pub struct PreparedHalfPlaneQuadrupole {
    acoustic: FluidAcousticState,
    frequency_hz: f64,
    frame: HalfPlaneFrame,
    boundary: IdealPlaneBoundary,
    source_position_m: [f64; 3],
    receiver: HalfPlaneReceiver,
    transfer: [[C64; 6]; 4],
    diagnostics: HalfPlaneQuadrupoleDiagnostics,
}
impl PreparedHalfPlaneQuadrupole {
    /// Resolve the complete pressure/velocity transfer with the SAME bounded
    /// adaptive quadrature as volume and force sources. Limits apply to full
    /// Hessian/third-derivative Frobenius norms: the 24 integrated channels
    /// use sqrt(2) for symmetric shear. They do not promise relative accuracy
    /// at a particular source-pattern or coherent interference null.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        acoustic: FluidAcousticState,
        frequency_hz: f64,
        frame: HalfPlaneFrame,
        boundary: IdealPlaneBoundary,
        source_position_m: [f64; 3],
        receiver: HalfPlaneReceiver,
        limits: HalfPlaneQuadrupoleLimits,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, AeroacError> {
        poll(&mut cancelled)?;
        let controls = QuadratureLimits::from(limits);
        controls.admit()?;
        let omega = TAU * frequency_hz;
        let k = omega / acoustic.sound_speed_m_s();
        let rho_omega = omega * acoustic.density_kg_m3();
        if !frequency_hz.is_finite()
            || frequency_hz <= 0.0
            || !k.is_finite()
            || k <= 0.0
            || !rho_omega.is_finite()
            || rho_omega <= 0.0
        {
            return Err(invalid(
                "half-plane stress needs a resolved positive frequency",
            ));
        }
        let source = frame.local_receiver(HalfPlaneReceiver::WorldPosition(source_position_m))?;
        let observer = frame.local_receiver(receiver)?;
        let kernel = StressKernel::new(
            source,
            observer,
            k,
            rho_omega,
            boundary,
            frame.axes(),
            controls,
        )?;
        let result = integrate_spectrum::<24, 6>(&kernel, [0.0; 24], controls, &mut cancelled)?;
        let scale = |j| if j < 3 { 1.0 } else { SQRT_2 };
        let transfer = core::array::from_fn(|i| {
            core::array::from_fn(|j| result.value[6 * i + j].scale(scale(j)))
        });
        let diagnostics = HalfPlaneQuadrupoleDiagnostics {
            evaluations: result.evaluations,
            maximum_depth_reached: result.maximum_depth_reached,
            retained_panels: result.retained_panels,
            minimum_edge_path_m: kernel.point.distance,
            pressure_quadrature_error_pa_per_n_m: core::array::from_fn(|j| {
                result.error[j] * scale(j)
            }),
            velocity_quadrature_error_m_s_per_n_m: core::array::from_fn(|i| {
                core::array::from_fn(|j| result.error[6 * (i + 1) + j] * scale(j))
            }),
        };
        if transfer.iter().flatten().any(|&v| !finite(v)) {
            return Err(invalid("half-plane physical stress map is unrepresentable"));
        }
        poll(&mut cancelled)?;
        Ok(Self {
            acoustic,
            frequency_hz,
            frame,
            boundary,
            source_position_m,
            receiver,
            transfer,
            diagnostics,
        })
    }
    pub const fn acoustic(&self) -> FluidAcousticState {
        self.acoustic
    }
    pub const fn frequency_hz(&self) -> f64 {
        self.frequency_hz
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
    /// Rows p,ux,uy,uz; columns physical Mxx,Myy,Mzz,Mxy,Mxz,Myz.
    pub const fn transfer(&self) -> &[[C64; 6]; 4] {
        &self.transfer
    }
    pub const fn diagnostics(&self) -> &HalfPlaneQuadrupoleDiagnostics {
        &self.diagnostics
    }
    pub fn evaluate_harmonic_line(
        &self,
        source: QuadrupolePhasor,
    ) -> Result<HarmonicFieldLine<3>, AeroacError> {
        if source.stress_moment_n_m.into_iter().any(|v| !finite(v)) {
            return Err(invalid(
                "half-plane stress moment must have finite physical entries",
            ));
        }
        let mut values = [C64::ZERO; 4];
        for (i, row) in self.transfer.iter().enumerate() {
            let (mut re, mut im) = (Dd::ZERO, Dd::ZERO);
            for (&h, m) in row.iter().zip(source.stress_moment_n_m) {
                let (a, b, c, d) = (
                    Dd::from_f64(h.re),
                    Dd::from_f64(h.im),
                    Dd::from_f64(m.re),
                    Dd::from_f64(m.im),
                );
                re = re + a * c - b * d;
                im = im + a * d + b * c;
            }
            values[i] = C64::new(re.to_f64(), im.to_f64());
        }
        if values.into_iter().any(|v| !finite(v)) {
            return Err(invalid("half-plane stress field is unrepresentable"));
        }
        Ok(HarmonicFieldLine {
            frequency_hz: self.frequency_hz,
            pressure_pa: values[0],
            particle_velocity_m_s: [values[1], values[2], values[3]],
        })
    }
    pub fn evaluate(&self, source: QuadrupolePhasor) -> Result<CoherentCompactField, AeroacError> {
        let line = self.evaluate_harmonic_line(source)?;
        let flux = line
            .particle_velocity_m_s
            .map(|v| line.pressure_pa.scale(0.5) * v.conj());
        if flux.into_iter().any(|v| !finite(v)) {
            return Err(invalid("half-plane stress intensity is unrepresentable"));
        }
        Ok(CoherentCompactField {
            pressure_pa: line.pressure_pa,
            particle_velocity_m_s: line.particle_velocity_m_s,
            active_intensity_w_m2: flux.map(|v| v.re),
            reactive_intensity_w_m2: flux.map(|v| v.im),
            scope: SCOPE,
        })
    }
    pub fn evaluate_energy_transport(
        &self,
        source: QuadrupolePhasor,
    ) -> Result<crate::curle3d::energy::AcousticEnergyTransport, AeroacError> {
        super::super::energy_transport(
            self.acoustic,
            [0.0; 3],
            self.evaluate_harmonic_line(source)?,
        )
    }
}

// exp(z), with a bounded outgoing real part. The original outer Green phase
// is evaluated separately in DD; this only evaluates the relative contour.
fn exponential(z: C64) -> C64 {
    let (s, c) = det::sin_cos(z.im);
    C64::new(c, s).scale(det::exp(z.re))
}
// phi_n(z)=(exp(z)-sum_{j<n} z^j/j!)/z^n and its first three derivatives.
// The convergent series avoids all small-argument divided differences.
fn phi(z: Jet, order: usize) -> Jet {
    let x = z.value;
    let mut f = [C64::ZERO; 4];
    if x.abs() < 2.0 {
        for (d, value) in f.iter_mut().enumerate() {
            let mut coefficient = 1.0;
            for j in 1..=order + d {
                coefficient /= j as f64;
            }
            for j in 1..=d {
                coefficient *= j as f64;
            }
            let mut term = C64::from_re(coefficient);
            let mut sum = term;
            for n in 1..48 {
                term = (term * x).scale((n + d) as f64 / (n * (n + order + d)) as f64);
                sum = sum + term;
                if term.abs() < 2e-17 * sum.abs() {
                    break;
                }
            }
            *value = sum;
        }
    } else {
        let e = exponential(x);
        let inverse = x.recip();
        f[0] = (e - C64::ONE) * inverse;
        for d in 1..4 {
            f[d] = (e - f[d - 1].scale(d as f64)) * inverse;
        }
        if order == 2 {
            f[0] = (f[0] - C64::ONE) * inverse;
            for d in 1..4 {
                f[d] = (f[d] - f[d - 1].scale(d as f64)) * inverse;
            }
        }
    }
    z.chain(f[0], f[1], f[2], f[3])
}
// phi(x)=sign(x)*ln(1+|x|)+x|x|/2 is C3 at zero. These exact derivatives
// retain the original cancellation without a 1/q or a finite difference.
fn shadow_primitive(x: Jet) -> Jet {
    let value = x.value.re;
    let a = value.abs();
    let f = if a < 0.125 {
        let mut sum = 0.0;
        let mut power = 1.0;
        for n in 0..24 {
            sum += power / (n + 3) as f64;
            power *= -a;
        }
        value * (1.0 + a * a * sum)
    } else {
        value.signum() * (det::ln(1.0 + a) + 0.5 * a * a)
    };
    let inv = 1.0 / (1.0 + a);
    x.chain(
        C64::from_re(f),
        C64::from_re(1.0 + a * a * inv),
        C64::from_re(value * (2.0 + a) * inv * inv),
        C64::from_re(2.0 * inv * inv * inv),
    )
}
fn phase_jet(l: Jet, k: f64) -> Jet {
    let ik = C64::new(0.0, k);
    let phase = cis(Dd::from_f64(k) * Dd::from_f64(l.value.re - 1.0));
    l.chain(phase, phase * ik, phase * ik * ik, phase * ik * ik * ik)
}

struct StressKernel {
    point: PointKernel,
    length: Jet,
    angle: [Jet; 2],
    base: [C64; 24],
    pressure_scale: f64,
    velocity_scale: f64,
    phase_floor: f64,
    breaks: [f64; 4],
    count: usize,
}
impl StressKernel {
    #[allow(clippy::too_many_arguments)]
    fn new(
        source: [f64; 3],
        receiver: [f64; 3],
        k: f64,
        rho_omega: f64,
        boundary: IdealPlaneBoundary,
        axes: [[f64; 3]; 3],
        controls: QuadratureLimits,
    ) -> Result<Self, AeroacError> {
        let point = PointKernel::new(source, receiver, k, 1.0, boundary, axes, controls)?;
        let l = point.distance;
        let pressure_scale = 1.0 / l / l / l;
        let velocity_scale = pressure_scale / l / rho_omega;
        if !pressure_scale.is_finite()
            || pressure_scale <= 0.0
            || !velocity_scale.is_finite()
            || velocity_scale <= 0.0
        {
            return Err(invalid(
                "half-plane stress transfer scales are unrepresentable",
            ));
        }
        let sx = Jet::coordinate(source[0] / l, axes[0], [0.0; 3]);
        let sy = Jet::coordinate(source[1] / l, axes[1], [0.0; 3]);
        let rx = Jet::coordinate(receiver[0] / l, [0.0; 3], axes[0]);
        let ry = Jet::coordinate(receiver[1] / l, [0.0; 3], axes[1]);
        let dz = (Dd::from_f64(receiver[2]) - Dd::from_f64(source[2])).to_f64() / l;
        let rz = Jet::coordinate(dz, axes[2].map(|a| -a), axes[2]);
        let rs = (sx * sx + sy * sy).sqrt();
        let r = (rx * rx + ry * ry).sqrt();
        let mut length = ((r + rs) * (r + rs) + rz * rz).sqrt();
        length.value = C64::ONE;
        let ts = (sx + sy.times(C64::new(0.0, 1.0))).argument_derivatives();
        let tr = (rx + ry.times(C64::new(0.0, 1.0))).argument_derivatives();
        let [ct, st] = half_angle(
            receiver[0] / det::hypot(receiver[0], receiver[1]),
            receiver[1] / det::hypot(receiver[0], receiver[1]),
        );
        let [cs, ss] = half_angle(
            source[0] / det::hypot(source[0], source[1]),
            source[1] / det::hypot(source[0], source[1]),
        );
        let cosine = [ct * cs + st * ss, ct * cs - st * ss];
        let sine = [st * cs - ct * ss, st * cs + ct * ss];
        let root = (r * rs).sqrt().scale(2.0);
        let angle = core::array::from_fn(|t| {
            let phi = (tr + ts.scale(if t == 0 { -1.0 } else { 1.0 })).scale(0.5);
            let mut q = root
                * phi.chain(
                    C64::from_re(cosine[t]),
                    C64::from_re(-sine[t]),
                    C64::from_re(-cosine[t]),
                    C64::from_re(sine[t]),
                );
            q.value = C64::from_re(point.eta[t]);
            q
        });
        let kap = point.kappa.to_f64();
        let one = Jet::real(1.0);
        let inv = length.recip();
        let e = phase_jet(length, kap);
        let edge = e * inv;
        let coefficient = edge * (length.times(C64::new(0.0, kap)) - one) * inv * inv;
        let mut base_jet = Jet::real(0.0);
        for t in 0..2 {
            let q = angle[t];
            let h = q * inv;
            let mut term = edge.scale(1.0 / (8.0 * PI))
                - (coefficient * q * q).scale(1.0 / (16.0 * PI))
                - (coefficient * length * length * shadow_primitive(h)).scale(1.0 / (8.0 * PI));
            if q.value.re > 0.0 {
                let residual = if h.value.re < 0.5 {
                    let delta = (h * h).scale(-1.0) * ((one - h * h).sqrt() + one).recip();
                    let ikl = length.times(C64::new(0.0, kap));
                    edge * delta
                        * delta
                        * (one + delta).recip()
                        * (ikl * ikl * phi(ikl * delta, 2)
                            - (ikl - one) * (Jet::real(1.5) + delta.scale(0.5)))
                } else {
                    let reflection = if t == 0 { 1.0 } else { -1.0 };
                    let dx = rx - sx;
                    let dy = ry - sy.scale(reflection);
                    let d = (dx * dx + dy * dy + rz * rz).sqrt();
                    phase_jet(d, kap) * d.recip() - edge + (coefficient * q * q).scale(0.5)
                };
                term = term + residual.scale(1.0 / (4.0 * PI));
            }
            base_jet = base_jet + term.scale(if t == 0 { 1.0 } else { point.sign });
        }
        if !length.all_finite()
            || angle.into_iter().any(|v| !v.all_finite())
            || !base_jet.all_finite()
        {
            return Err(invalid(
                "half-plane stress geometry derivatives are unrepresentable",
            ));
        }
        let geometry = 1.0 + 1.0 / r.value.re + 1.0 / rs.value.re;
        // A conservative derivative/geometry envelope only turns OFF phase
        // pre-splitting in an exponentially small tail; every rational and
        // outgoing term still passes the embedded quadrature error test.
        let phase_floor = det::ln(100_000.0)
            + 4.0 * det::ln(1.0 + kap)
            + 3.0 * det::ln(geometry)
            + (det::ln(pressure_scale) - det::ln(controls.absolute_pressure))
                .max(det::ln(velocity_scale) - det::ln(controls.absolute_velocity));
        let mut candidates = [0.0, point.eta[0].abs(), point.eta[1].abs(), point.scale];
        candidates.sort_by(f64::total_cmp);
        let mut breaks = [0.0; 4];
        let mut count = 1;
        for value in candidates {
            if value > breaks[count - 1] {
                breaks[count] = value;
                count += 1;
            }
        }
        let mut result = Self {
            point,
            length,
            angle,
            base: [C64::ZERO; 24],
            pressure_scale,
            velocity_scale,
            phase_floor,
            breaks,
            count,
        };
        result.base = result.channels(base_jet);
        if result.base.into_iter().any(|v| !finite(v)) || !phase_floor.is_finite() {
            return Err(invalid(
                "half-plane stress base or phase bound is unrepresentable",
            ));
        }
        Ok(result)
    }
    fn channels(&self, value: Jet) -> [C64; 24] {
        let mut result = [C64::ZERO; 24];
        for j in 0..6 {
            let shear = if j < 3 { 1.0 } else { SQRT_2 };
            result[j] = (self.point.phase * value.source2[j]).scale(self.pressure_scale * shear);
            for i in 0..3 {
                let v = self.point.phase * value.mixed2[i][j];
                result[6 * (i + 1) + j] = C64::new(v.im, -v.re).scale(self.velocity_scale * shear);
            }
        }
        result
    }
    fn coordinate(&self, region: usize, x: f64) -> (f64, f64) {
        let lower = self.breaks[region];
        if region + 1 < self.count {
            let upper = self.breaks[region + 1];
            if lower == 0.0 || upper <= 2.0 * lower {
                let width = upper - lower;
                (lower + width * x, width)
            } else {
                let start = det::ln(lower);
                let width = det::ln(upper) - start;
                let value = det::exp(start + width * x);
                (value, width * value)
            }
        } else {
            let rest = 1.0 - x;
            (
                lower + self.point.scale * x / rest,
                self.point.scale / rest / rest,
            )
        }
    }
}
impl Spectrum<24> for StressKernel {
    fn regions(&self) -> usize {
        self.count
    }
    fn extent(&self) -> (f64, f64) {
        (0.0, 0.0)
    }
    fn resolved(&self, region: usize, lower: f64, upper: f64) -> bool {
        let (a, _) = self.coordinate(region, lower);
        let (b, _) = self.coordinate(region, upper);
        (a >= 1.0 && self.point.phase_extent(a).1 >= self.phase_floor.max(0.0))
            || self.point.phase_extent(b).0 - self.point.phase_extent(a).0 <= FRAC_PI_2
    }
    fn sample(&self, region: usize, x: f64) -> Result<[C64; 24], AeroacError> {
        let (w, jac) = self.coordinate(region, x);
        let s = C64::new(w * FRAC_1_SQRT_2, w * FRAC_1_SQRT_2);
        let square = C64::new(0.0, w * w);
        let l = self.length;
        let r = (l * l + Jet::constant(square)).sqrt();
        let delta = Jet::constant(square) * (r + l).recip();
        let ik = C64::new(0.0, self.point.kappa.to_f64());
        // (f'(s)-r(s))/s³. All factors remain regular at the contour origin.
        let remainder = phase_jet(l, self.point.kappa.to_f64())
            * r.recip()
            * r.recip()
            * r.recip()
            * (r + l).recip()
            * (phi(delta.times(ik), 1) * (r.times(ik) - Jet::real(1.0)).times(ik) + l.recip());
        let mut integral = Jet::real(0.0);
        for t in 0..2 {
            let q = self.angle[t];
            let a = q.value.re;
            let scale = a.abs().max(w);
            let sn = s.scale(1.0 / scale);
            let an = a / scale;
            let ss = sn * sn;
            let den = (ss + C64::from_re(an * an)).recip();
            let common = ss * ss * den;
            // s³*atan(q/s) and its first three q derivatives. Scale factors
            // cancel BEFORE forming the third derivative at a narrow shadow.
            let v = s * s * s * contour_arctangent(a, w);
            let first = common.scale(scale * scale);
            let second = (common * den).scale(-2.0 * an * scale);
            let third = common * den * den * (C64::from_re(3.0 * an * an) - ss).scale(2.0);
            let angular = q.chain(v, first, second, third);
            integral =
                integral + (remainder * angular).scale(if t == 0 { 1.0 } else { self.point.sign });
        }
        let contour = integral
            .times(C64::new(FRAC_1_SQRT_2, FRAC_1_SQRT_2))
            .scale(-jac / (4.0 * PI * PI));
        let correction = self.channels(contour);
        let result = core::array::from_fn(|i| {
            super::super::add(self.base[i].scale(1.0 / self.count as f64), correction[i])
        });
        if !w.is_finite() || !jac.is_finite() || result.into_iter().any(|v| !finite(v)) {
            return Err(invalid("half-plane stress contour is unrepresentable"));
        }
        Ok(result)
    }
}
#[cfg(test)]
mod tests;
