//! Exact plane and spherical source diffraction by one ideal half-plane.
//!
//! The sheet occupies x >= 0, y = 0, with its straight edge on the z axis.
//! Coordinates and output velocity use that same orthonormal physical frame.
//! Both faces are retained: y = +0 selects the upper face and y = -0 the
//! lower face. No geometric visibility cutoff discards the shadow field.
//! [`PreparedHalfPlaneSpectrum`] places this same exact edge in world space
//! and retains coherent interference between supplied incident plane waves.
//! [`PreparedHalfPlaneMonopole`] resolves the exact finite-distance spherical
//! Green function and analytic velocity with bounded contour quadrature.
//! [`PreparedHalfPlaneForce`] differentiates the same Green function in source
//! and receiver coordinates for physical point forces and full Euler velocity.
//!
//! For exp(-i omega t), arrival azimuth alpha and receiver azimuth theta,
//! the transverse solution is
//! `exp(-i kr cos(theta-alpha)) F(sqrt(2kr) cos((theta-alpha)/2))`
//! plus/minus its alpha -> -alpha reflection, for rigid/release boundaries.
//! `F(t) = exp(-i pi/4)/sqrt(pi) integral(-infinity..t) exp(i s²) ds`.
//! This is Sommerfeld's exact Helmholtz solution, not the paraxial knife-edge
//! approximation. See A3--A4 in
//! <https://pmc.ncbi.nlm.nih.gov/articles/PMC4844525/>. The edge-parallel
//! wave-number component separates as exp(i kz z), with transverse k = sqrt(k²-kz²).
//! [`HalfPlaneDiffraction::with_edge_parallel_flow`] retains that separation
//! for uniform subsonic flow along the edge. Each incident direction has
//! intrinsic frequency `Omega = omega / (1 + U n_z/c)` and `k = Omega/c`;
//! the observed frequency remains the original laboratory `omega`.
//!
//! Analytic derivatives retain the singular 1/sqrt(r) edge velocity and
//! signed active/reactive flux. The edge itself therefore refuses. This
//! stationary, zero-thickness, infinite, ideal sheet does not approximate
//! finite barriers, impedance coatings, arbitrary wedges,
//! absorption, flow across the edge, or reciprocal source loading.

use crate::AeroacError;
use crate::curle3d::{CoherentCompactField, IdealPlaneBoundary};
use crate::observation::HarmonicFieldLine;
use fs_material::fluid_acoustics::FluidAcousticState;
use fs_math::{c64::C64, dd::Dd, det};

mod point;
mod spectrum;
pub use point::{
    HalfPlaneForceDiagnostics, HalfPlaneForceLimits, HalfPlaneMonopoleDiagnostics,
    HalfPlaneMonopoleLimits, PreparedHalfPlaneForce, PreparedHalfPlaneMonopole,
};
pub use spectrum::{
    HalfPlaneFace, HalfPlaneFrame, HalfPlaneIncidentWave, HalfPlaneReceiver,
    HalfPlaneSpectrumBudget, MAX_HALF_PLANE_COMPONENTS, PreparedHalfPlaneSpectrum,
};

const INV_SQRT_TWO_PI: f64 = 0.398_942_280_401_432_7;
const SCOPE: &str = "Prescribed plane-wave diffraction by one infinite straight-edged, \
zero-thickness rigid or pressure-release half-plane in a stationary homogeneous lossless \
fluid; exact Sommerfeld pressure and Euler velocity, including shadow and transition fields. \
No finite barrier, finite-impedance wall, finite-source, arbitrary-wedge, flow, absorption, \
source-feedback or measured-SPL qualification.";

const FLOW_SCOPE: &str = "Prescribed plane-wave diffraction by one infinite straight-edged, \
zero-thickness rigid or pressure-release half-plane with uniform subsonic flow parallel \
to its edge; exact separated linearized-Euler pressure and perturbation velocity. \
Active intensity is pressure work, excluding mean-flow energy advection. No cross-edge \
flow, shear, finite barrier, wall impedance, absorption, source feedback or measured-SPL qualification.";

/// One material and incident direction bound to the exact half-plane solution.
///
/// The complex input amplitude is pressure [Pa] at the edge-frame origin on
/// the source clock, before diffraction. The unit direction points along
/// propagation, toward the sheet. Oblique incidence along its edge is retained;
/// exactly grazing incidence parallel to the sheet refuses rather than choosing
/// an unspecified limiting face. The infinite geometry requires no mesh.
#[derive(Clone, Copy, Debug)]
pub struct HalfPlaneDiffraction {
    acoustic: FluidAcousticState,
    frequency_hz: f64,
    intrinsic_omega: f64,
    edge_flow_velocity_m_s: f64,
    k: f64,
    transverse_k: f64,
    direction: [f64; 3],
    arrival_half: [f64; 2],
    boundary: IdealPlaneBoundary,
}

impl HalfPlaneDiffraction {
    /// Bind an admitted homogeneous material, frequency [Hz], unit propagation
    /// direction and ideal boundary. No density or sound speed is defaulted.
    pub fn new(
        acoustic: FluidAcousticState,
        frequency_hz: f64,
        propagation_direction: [f64; 3],
        boundary: IdealPlaneBoundary,
    ) -> Result<Self, AeroacError> {
        Self::with_edge_parallel_flow(acoustic, frequency_hz, propagation_direction, boundary, 0.0)
    }

    /// Bind a stationary sheet in uniform subsonic mean flow along its z edge.
    /// The direction is the incident WAVE NORMAL, not its laboratory ray/group
    /// velocity (which is `c*n + U*e_z`). The supplied frequency and incident
    /// pressure refer to the laboratory clock and the edge-frame origin.
    ///
    /// Translational symmetry retains kz in every diffracted term. Thus the
    /// material derivative is exactly `-i*(omega-U*kz)`, and the original
    /// Sommerfeld pressure AND velocity solve the linearized Euler system at
    /// that intrinsic frequency. No effective sound speed or stationary-fluid
    /// pressure-to-velocity conversion is used. Rigid/release boundary data
    /// remain valid because the mean flow is tangent to both faces and edge.
    /// Cross-edge flow, shear and a Kutta/wake model are outside this owner.
    pub fn with_edge_parallel_flow(
        acoustic: FluidAcousticState,
        frequency_hz: f64,
        propagation_direction: [f64; 3],
        boundary: IdealPlaneBoundary,
        edge_flow_velocity_m_s: f64,
    ) -> Result<Self, AeroacError> {
        let [dx, dy, dz] = propagation_direction;
        let transverse = det::hypot(dx, dy);
        let length = det::hypot(transverse, dz);
        let omega = core::f64::consts::TAU * frequency_hz;
        let mach = edge_flow_velocity_m_s / acoustic.sound_speed_m_s();
        let intrinsic_omega = omega / (1.0 + mach * dz);
        let k = intrinsic_omega / acoustic.sound_speed_m_s();
        let transverse_k = k * transverse;
        if !frequency_hz.is_finite()
            || frequency_hz <= 0.0
            || !edge_flow_velocity_m_s.is_finite()
            || !mach.is_finite()
            || mach.abs() >= 1.0
            || propagation_direction.iter().any(|v| !v.is_finite())
            || (length - 1.0).abs() > 64.0 * f64::EPSILON
            || dy == 0.0
            || !k.is_finite()
            || k <= 0.0
            || !transverse_k.is_finite()
            || transverse_k <= 0.0
            || !(intrinsic_omega * acoustic.density_kg_m3()).is_finite()
            || intrinsic_omega * acoustic.density_kg_m3() <= 0.0
        {
            return Err(invalid(
                "half-plane needs a finite non-grazing unit wave, subsonic edge flow and positive resolved frequency",
            ));
        }
        Ok(Self {
            acoustic,
            frequency_hz,
            intrinsic_omega,
            edge_flow_velocity_m_s,
            k,
            transverse_k,
            direction: propagation_direction,
            arrival_half: half_angle(-dx / transverse, -dy / transverse),
            boundary,
        })
    }

    /// The immutable original medium used for wavelength and particle motion.
    pub const fn acoustic(&self) -> FluidAcousticState {
        self.acoustic
    }

    /// Original laboratory frequency [Hz], unchanged by a stationary sheet.
    pub const fn laboratory_frequency_hz(&self) -> f64 {
        self.frequency_hz
    }

    /// Material-frame angular frequency `omega - U*kz` [rad/s].
    pub const fn intrinsic_angular_frequency_rad_s(&self) -> f64 {
        self.intrinsic_omega
    }

    /// Incident spatial wave-number magnitude [rad/m], relative to the fluid.
    pub const fn wavenumber_rad_m(&self) -> f64 {
        self.k
    }

    /// Original uniform mean-fluid velocity [m/s] in the edge frame.
    pub const fn mean_flow_velocity_m_s(&self) -> [f64; 3] {
        [0.0, 0.0, self.edge_flow_velocity_m_s]
    }

    /// Coherent total pressure, Cartesian perturbation velocity and pressure work.
    /// In flow, active intensity is `<p u>`; transported perturbation energy
    /// additionally contains `U*<E>` (see [`Self::evaluate_energy_transport`]).
    ///
    /// Coordinates are metres in the documented edge frame. Signed y = +/-0
    /// selects a sheet face for x > 0. At x < 0 both signed zeros describe the
    /// same open fluid. r = 0 refuses because the edge velocity is singular.
    /// Evaluation is allocation-free with at most two 256-term transition
    /// evaluations; cancellation is polled at most eight iterations apart.
    pub fn evaluate(
        &self,
        position_m: [f64; 3],
        incident_pressure_pa: C64,
        cancelled: impl FnMut() -> bool,
    ) -> Result<CoherentCompactField, AeroacError> {
        let [pressure, ux, uy, uz] =
            self.amplitudes(position_m, incident_pressure_pa, cancelled)?;
        let velocity = [ux, uy, uz];
        let flux = velocity.map(|value| pressure.scale(0.5) * value.conj());
        if flux.iter().any(|&value| !finite(value)) {
            return Err(invalid("half-plane intensity is not representable"));
        }
        Ok(CoherentCompactField {
            pressure_pa: pressure,
            particle_velocity_m_s: velocity,
            active_intensity_w_m2: flux.map(|value| value.re),
            reactive_intensity_w_m2: flux.map(|value| value.im),
            scope: if self.edge_flow_velocity_m_s == 0.0 {
                SCOPE
            } else {
                FLOW_SCOPE
            },
        })
    }

    /// Time-averaged quadratic energy and its laboratory flux, including mean
    /// advection. Reduces the complete coherent field, preserving near-edge
    /// motion; this is not a plane-wave pressure-squared approximation.
    pub fn evaluate_energy_transport(
        &self,
        position_m: [f64; 3],
        incident_pressure_pa: C64,
        cancelled: impl FnMut() -> bool,
    ) -> Result<crate::curle3d::energy::AcousticEnergyTransport, AeroacError> {
        let line = self.evaluate_harmonic_line(position_m, incident_pressure_pa, cancelled)?;
        energy_transport(self.acoustic, self.mean_flow_velocity_m_s(), line)
    }

    /// Feed the original finite-window/sampled acoustic observation pipeline.
    /// All direct, reflected and diffracted phases are already applied; no
    /// second delay or plane-wave pressure-to-velocity approximation is needed.
    pub fn evaluate_harmonic_line(
        &self,
        position_m: [f64; 3],
        incident_pressure_pa: C64,
        cancelled: impl FnMut() -> bool,
    ) -> Result<HarmonicFieldLine<3>, AeroacError> {
        let [pressure, ux, uy, uz] =
            self.amplitudes(position_m, incident_pressure_pa, cancelled)?;
        Ok(HarmonicFieldLine {
            frequency_hz: self.frequency_hz,
            pressure_pa: pressure,
            particle_velocity_m_s: [ux, uy, uz],
        })
    }

    fn amplitudes(
        &self,
        [x, y, z]: [f64; 3],
        amplitude: C64,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<[C64; 4], AeroacError> {
        poll(&mut cancelled)?;
        let radius = det::hypot(x, y);
        let kr = Dd::from_f64(self.transverse_k) * Dd::from_f64(radius);
        if [x, y, z].iter().any(|v| !v.is_finite())
            || !finite(amplitude)
            || radius <= 0.0
            || !kr.to_f64().is_finite()
            || kr.to_f64() <= 0.0
            || !(2.0 * kr.to_f64()).is_finite()
        {
            return Err(invalid(
                "half-plane receiver must resolve a finite distance from the edge",
            ));
        }
        let [ct, st] = half_angle(x / radius, y / radius);
        let [ca, sa] = self.arrival_half;
        let minus = [ct * ca + st * sa, st * ca - ct * sa];
        let plus = [ct * ca - st * sa, st * ca + ct * sa];
        let k = Dd::from_f64(self.k);
        let axial_phase = k * Dd::from_f64(self.direction[2]) * Dd::from_f64(z);
        let phase_x = k * Dd::from_f64(self.direction[0]) * Dd::from_f64(x);
        let phase_y = k * Dd::from_f64(self.direction[1]) * Dd::from_f64(y);
        let edge_phase = cis(kr + axial_phase);
        let root = det::sqrt(2.0 * kr.to_f64());
        let direct = term(
            root * minus[0],
            cis(phase_x + phase_y + axial_phase),
            edge_phase,
            &mut cancelled,
        )?;
        let reflected = term(
            root * plus[0],
            cis(phase_x - phase_y + axial_phase),
            edge_phase,
            &mut cancelled,
        )?;
        let sign = self.boundary.sign();
        let pressure = add(direct, reflected.scale(sign));
        let difference = add(direct, reflected.scale(-sign));
        // exp(i phase) F'(t) has phase kr + kz*z - pi/4 for BOTH terms.
        // Keep that identity before evaluating trig, especially in the shadow
        // where separately rounded large phases would corrupt a weak field.
        let derivative = edge_phase
            * C64::new(INV_SQRT_TWO_PI, -INV_SQRT_TWO_PI)
                .scale(core::f64::consts::FRAC_1_SQRT_2 * det::sqrt(self.transverse_k / radius));
        let inverse_rho_omega = 1.0 / (self.intrinsic_omega * self.acoustic.density_kg_m3());
        let edge_velocity = C64::new(derivative.im, -derivative.re).scale(inverse_rho_omega);
        let inverse_impedance = 1.0 / self.acoustic.characteristic_impedance_pa_s_m();
        let velocity = [
            add(
                pressure.scale(self.direction[0] * inverse_impedance),
                edge_velocity.scale(plus[0] + sign * minus[0]),
            ),
            add(
                difference.scale(self.direction[1] * inverse_impedance),
                edge_velocity.scale(plus[1] + sign * minus[1]),
            ),
            pressure.scale(self.direction[2] * inverse_impedance),
        ];
        let result = [pressure, velocity[0], velocity[1], velocity[2]].map(|v| amplitude * v);
        if result.iter().any(|&value| !finite(value)) {
            return Err(invalid(
                "half-plane pressure or edge velocity is not representable",
            ));
        }
        poll(&mut cancelled)?;
        Ok(result)
    }
}

// A peak harmonic's mean quadratic energy is the half-sum of its two real
// quadratures. Reuse the original energy owner so both material modulus and
// mean-flow transport have the same semantics as the other acoustic solvers.
fn energy_transport(
    acoustic: FluidAcousticState,
    flow: [f64; 3],
    line: HarmonicFieldLine<3>,
) -> Result<crate::curle3d::energy::AcousticEnergyTransport, AeroacError> {
    use crate::curle3d::energy;
    let real = energy::instantaneous(
        acoustic,
        flow,
        line.pressure_pa.re,
        line.particle_velocity_m_s.map(|v| v.re),
    )?;
    let imag = energy::instantaneous(
        acoustic,
        flow,
        line.pressure_pa.im,
        line.particle_velocity_m_s.map(|v| v.im),
    )?;
    let half_sum =
        |a, b| (Dd::from_f64(a) * Dd::from_f64(0.5) + Dd::from_f64(b) * Dd::from_f64(0.5)).to_f64();
    Ok(energy::AcousticEnergyTransport {
        energy_density_j_m3: half_sum(real.energy_density_j_m3, imag.energy_density_j_m3),
        pressure_work_flux_w_m2: core::array::from_fn(|i| {
            half_sum(
                real.pressure_work_flux_w_m2[i],
                imag.pressure_work_flux_w_m2[i],
            )
        }),
        energy_flux_w_m2: core::array::from_fn(|i| {
            half_sum(real.energy_flux_w_m2[i], imag.energy_flux_w_m2[i])
        }),
    })
}

// theta is in [0, 2pi], so sin(theta/2) >= 0 and the sign of cos(theta/2)
// distinguishes both sheet faces. Stable half-angle identities avoid atan2
// subtraction and preserve the lower face's signed zero exactly.
fn half_angle(cos: f64, sin: f64) -> [f64; 2] {
    if cos >= 0.0 {
        let half_cos = det::sqrt(0.5 * (1.0 + cos)).copysign(sin);
        [half_cos, sin.abs() / (2.0 * half_cos.abs())]
    } else {
        let half_sin = det::sqrt(0.5 * (1.0 - cos));
        [sin / (2.0 * half_sin), half_sin]
    }
}

fn cis(phase: Dd) -> C64 {
    let (sin, cos) = det::sin_cos_dd(phase);
    C64::new(cos.to_f64(), sin.to_f64())
}

fn term(
    t: f64,
    plane: C64,
    edge: C64,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<C64, AeroacError> {
    if t.abs() < 1.5 {
        return Ok(plane * transition_series(t, cancelled)?);
    }
    let tail = edge * scaled_tail(t.abs(), cancelled)?;
    Ok(if t < 0.0 { tail } else { plane - tail })
}

// Defining convergent Fresnel series, with DD reduction. Only |t| < 1.5
// enters; there is no large-argument subtractive cancellation in this sum.
fn transition_series(t: f64, cancelled: &mut impl FnMut() -> bool) -> Result<C64, AeroacError> {
    let square = Dd::from_f64(t) * Dd::from_f64(t);
    let mut value = Dd::from_f64(t);
    let mut re = Dd::ZERO;
    let mut im = Dd::ZERO;
    for n in 0..64 {
        if n % 8 == 0 {
            poll(cancelled)?;
        }
        match n % 4 {
            0 => re = re + value,
            1 => im = im + value,
            2 => re = re - value,
            _ => im = im - value,
        }
        if value.to_f64().abs() < 1e-33 {
            return Ok(C64::new(
                (Dd::from_f64(0.5) + (re + im) * Dd::from_f64(INV_SQRT_TWO_PI)).to_f64(),
                ((im - re) * Dd::from_f64(INV_SQRT_TWO_PI)).to_f64(),
            ));
        }
        value = value * square * Dd::from_f64(f64::from(2 * n + 1))
            / Dd::from_f64(f64::from((n + 1) * (2 * n + 3)));
    }
    Err(invalid("half-plane transition series did not converge"))
}

// exp(-i x²) F(-x), x >= 1.5, from the incomplete-gamma continued fraction
// Gamma(1/2,-i x²). Contracting DLMF 8.9.2 gives b_n=z+1/2+2n,
// a_n=-n(n-1/2). The scaled tail avoids subtracting Fresnel C/S from 1/2
// in deep shadow; its phase is combined analytically by term above.
fn scaled_tail(x: f64, cancelled: &mut impl FnMut() -> bool) -> Result<C64, AeroacError> {
    let mut b = C64::new(0.5, -x * x);
    let mut c = C64::from_re(1e100);
    let mut d = b.recip();
    let mut h = d;
    for n in 1..=256 {
        if n % 8 == 1 {
            poll(cancelled)?;
        }
        let n = f64::from(n);
        let a = -n * (n - 0.5);
        b.re += 2.0;
        d = (d.scale(a) + b).recip();
        c = b + C64::from_re(a) / c;
        let delta = c * d;
        h = h * delta;
        if (delta - C64::ONE).abs() <= 4.0 * f64::EPSILON {
            return Ok(h * C64::new(1.0, -1.0).scale(0.5 * INV_SQRT_TWO_PI * x));
        }
    }
    Err(invalid(
        "half-plane transition continued fraction did not converge",
    ))
}

fn add(a: C64, b: C64) -> C64 {
    C64::new(
        (Dd::from_f64(a.re) + Dd::from_f64(b.re)).to_f64(),
        (Dd::from_f64(a.im) + Dd::from_f64(b.im)).to_f64(),
    )
}

fn finite(value: C64) -> bool {
    value.re.is_finite() && value.im.is_finite()
}

fn invalid(what: &'static str) -> AeroacError {
    AeroacError::InvalidParameter { what }
}

fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), AeroacError> {
    if cancelled() {
        Err(AeroacError::Cancelled {
            what: "half-plane diffraction",
        })
    } else {
        Ok(())
    }
}

#[cfg(test)]
mod tests;

#[cfg(test)]
mod flow_tests;
