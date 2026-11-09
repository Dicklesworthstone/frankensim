//! Exact uniform-flow image of the ORIGINAL finite-source Green function.
//!
//! Let M=U/c, beta²=1-M², z'=z/beta and alpha=k*M/beta². The
//! convected fundamental solution is
//! G_U(x,s)=exp(-i*alpha*(z-zs))*G_0(x',s';k/beta)/beta.
//! This is the Lorentz transform in Griesmaier/Raumer, §2.1, Proposition 2.1
//! and (2.2): <https://www.waves.kit.edu/downloads/CRC1173_Preprint_2022-16.pdf>.
//! The half-plane and both boundary conditions survive because U is tangent
//! to the sheet AND edge. A transverse or nonuniform flow does not.
//!
//! For the original continuity source Q, phi=-Q*G_U,
//! p/Q=rho*(-i*omega+U*partial_z)*G_U, and u/Q=-grad(G_U).
//! If p0=-i*omega*rho*G_0 and u0=-grad'(G_0), the physical fields are
//! p=phase*(p0/beta³-rho*U*u0z/beta²),
//! ux,y=phase*u0x,y/beta,
//! uz=phase*(u0z/beta²-U*p0/(K*beta³)).
//! The original omega and material coefficients are never replaced by a
//! scalar effective sound speed. Every quadrature sample is transformed
//! before its physical pressure/vector-velocity acceptance check.

use super::{PointKernel, cis, finite, invalid, rotate};
use crate::AeroacError;
use crate::curle3d::IdealPlaneBoundary;
use crate::curle3d::refraction::full_wave::{
    IntegratedSpectrum, QuadratureLimits, Spectrum, integrate_spectrum,
};
use fs_material::fluid_acoustics::FluidAcousticState;
use fs_math::{c64::C64, dd::Dd, det};

pub(super) const SCOPE: &str = "Prescribed stationary point volume source diffracted by \
one infinite straight-edged rigid or pressure-release sheet in uniform subsonic fluid \
flow parallel to its edge. Exact convected potential Green function with full physical \
pressure and perturbation velocity on the laboratory source clock; estimated quadrature \
error. Active intensity is pressure work; mean energy transport additionally includes \
fluid advection. No cross-edge/sheared flow, source motion, impedance coating, absorption, \
finite barrier, independent momentum injection, source feedback or measured-SPL qualification.";

struct FlowImage {
    phase: C64,
    axes: [[f64; 3]; 3],
    inverse_beta: f64,
    inverse_beta_squared: f64,
    pressure_gain: f64,
    pressure_from_velocity: f64,
    velocity_from_pressure: f64,
}

impl FlowImage {
    fn apply(&self, values: [C64; 4]) -> Result<[C64; 4], AeroacError> {
        // DD products retain the physical pressure/edge-velocity cancellation.
        let pair = |a: C64, x: f64, b: C64, y: f64| {
            C64::new(
                (Dd::from_f64(a.re) * Dd::from_f64(x) + Dd::from_f64(b.re) * Dd::from_f64(y))
                    .to_f64(),
                (Dd::from_f64(a.im) * Dd::from_f64(x) + Dd::from_f64(b.im) * Dd::from_f64(y))
                    .to_f64(),
            )
        };
        let pressure = pair(
            values[0],
            self.pressure_gain,
            values[3],
            -self.pressure_from_velocity,
        );
        let local = [
            values[1].scale(self.inverse_beta),
            values[2].scale(self.inverse_beta),
            pair(
                values[3],
                self.inverse_beta_squared,
                values[0],
                -self.velocity_from_pressure,
            ),
        ];
        let velocity = rotate(self.axes, local);
        let result =
            [pressure, velocity[0], velocity[1], velocity[2]].map(|value| self.phase * value);
        if result.iter().any(|&value| !finite(value)) {
            return Err(invalid(
                "convected half-plane pressure or velocity is unrepresentable",
            ));
        }
        Ok(result)
    }

    fn contour_limits(&self, controls: QuadratureLimits) -> Result<QuadratureLimits, AeroacError> {
        // These tightened bounds ONLY decide whether the original tail still
        // needs phase refinement. The shared integrator accepts the transformed
        // physical world-space samples using the caller's original controls.
        let pressure = (0.5 * (controls.absolute_pressure / self.pressure_gain))
            .min(0.5 * (controls.absolute_velocity / self.velocity_from_pressure.abs()));
        let velocity = (0.5 * (controls.absolute_pressure / self.pressure_from_velocity.abs()))
            .min(0.5 * (controls.absolute_velocity / self.inverse_beta_squared));
        let limits = QuadratureLimits {
            absolute_pressure: pressure,
            absolute_velocity: velocity,
            ..controls
        };
        limits.admit()?;
        Ok(limits)
    }
}

struct ConvectedKernel {
    original: PointKernel,
    image: FlowImage,
}
impl Spectrum<4> for ConvectedKernel {
    fn regions(&self) -> usize {
        self.original.regions()
    }
    fn extent(&self) -> (f64, f64) {
        self.original.extent()
    }
    fn resolved(&self, region: usize, lower: f64, upper: f64) -> bool {
        self.original.resolved(region, lower, upper)
    }
    fn sample(&self, region: usize, x: f64) -> Result<[C64; 4], AeroacError> {
        self.image.apply(self.original.sample(region, x)?)
    }
}

#[allow(clippy::too_many_arguments)]
pub(super) fn prepare(
    acoustic: FluidAcousticState,
    k: f64,
    omega_rho: f64,
    source: [f64; 3],
    point: [f64; 3],
    boundary: IdealPlaneBoundary,
    axes: [[f64; 3]; 3],
    velocity: f64,
    controls: QuadratureLimits,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(IntegratedSpectrum<4>, f64), AeroacError> {
    let mach = velocity / acoustic.sound_speed_m_s();
    if !velocity.is_finite() || !mach.is_finite() || mach.abs() >= 1.0 {
        return Err(invalid(
            "finite-source edge flow must be uniform and strictly subsonic",
        ));
    }
    let beta_squared = (1.0 - mach) * (1.0 + mach);
    let beta = det::sqrt(beta_squared);
    let inverse_beta = 1.0 / beta;
    let inverse_beta_squared = 1.0 / beta_squared;
    let pressure_gain = inverse_beta * inverse_beta_squared;
    let pressure_from_velocity = (Dd::from_f64(acoustic.density_kg_m3())
        * Dd::from_f64(velocity)
        * Dd::from_f64(inverse_beta_squared))
    .to_f64();
    let velocity_from_pressure = (Dd::from_f64(velocity) * Dd::from_f64(pressure_gain)
        / Dd::from_f64(acoustic.adiabatic_bulk_modulus_pa()))
    .to_f64();
    let dz = (Dd::from_f64(point[2]) - Dd::from_f64(source[2])).to_f64();
    let phase_argument =
        -(Dd::from_f64(k) * Dd::from_f64(mach) * Dd::from_f64(dz) / Dd::from_f64(beta_squared));
    let transformed_k = k * inverse_beta;
    let axial = dz * inverse_beta;
    if [
        pressure_gain,
        pressure_from_velocity,
        velocity_from_pressure,
        phase_argument.to_f64(),
        transformed_k,
        axial,
    ]
    .iter()
    .any(|x| !x.is_finite())
        || transformed_k <= 0.0
        || (velocity != 0.0 && (pressure_from_velocity == 0.0 || velocity_from_pressure == 0.0))
    {
        return Err(invalid(
            "finite-source Lorentz transform is unrepresentable",
        ));
    }
    let image = FlowImage {
        phase: cis(phase_argument),
        axes,
        inverse_beta,
        inverse_beta_squared,
        pressure_gain,
        pressure_from_velocity,
        velocity_from_pressure,
    };
    let contour_limits = image.contour_limits(controls)?;
    // Only relative axial position enters the straight-edge owner. Translate
    // first, then stretch, retaining small separations at large world coordinates.
    let original = PointKernel::new(
        [source[0], source[1], 0.0],
        [point[0], point[1], axial],
        transformed_k,
        omega_rho,
        boundary,
        [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
        contour_limits,
    )?;
    let physical_edge_path = det::hypot(
        det::hypot(source[0], source[1]) + det::hypot(point[0], point[1]),
        dz,
    );
    let kernel = ConvectedKernel { original, image };
    let result = integrate_spectrum::<4, 1>(&kernel, [0.0; 4], controls, cancelled)?;
    Ok((result, physical_edge_path))
}

#[cfg(test)]
mod tests;
