//! Physical force sources inside the original finite layer, observed outside.
//!
//! Let WA/WB be the parent's reciprocal pressure/normal-velocity weights.
//! Differentiating the original INTERIOR source gives A_s=i*rho2*B, hence
//! pF=[i*q*WA*J1/rho2, -WB*J0]/(2*pi). Exterior momentum uses rho1 and
//! d_observer_normal exp(i*zeta1*a)=-i*zeta1 exp(i*zeta1*a).
//! In the physical receiver radial basis, ordered [tt,ss,tn,nt,nn], velocity is
//! [q² WA (J0-J1/x)/rho2, q² WA J1/x/rho2, -i q WB J1,
//!  -i q zeta1 WA J1/rho2, zeta1 WB J0]/(2*pi*omega*rho1).
//! This radial basis is the NEGATIVE of the reciprocal interior kernel's
//! tangent. Normal always means the original first-face normal into the layer.

use super::super::{InteriorKernel, interior_coordinates, reduces_to_interface};
use crate::AeroacError;
use crate::curle3d::CoherentCompactField;
use crate::curle3d::refraction::PlanarFluidInterface;
use crate::curle3d::refraction::full_wave::force::{admit, j1_over_argument};
use crate::curle3d::refraction::full_wave::layer::{
    FullWaveLayerForceDiagnostics, PlanarFluidLayer,
};
use crate::curle3d::refraction::full_wave::reflection::{
    apply_force, field, force_cartesian, line, prepare_reflected_tail,
};
use crate::curle3d::refraction::full_wave::{
    FullWaveForceLimits, PreparedFullWaveForce, Spectrum, finite, integrate_spectrum, j0, j1,
    nonfinite, poll,
};
use crate::observation::HarmonicFieldLine;
use core::f64::consts::{PI, TAU};
use fs_math::{c64::C64, det};

/// Applicability of the original layer-force source and exterior observation.
pub const FULL_WAVE_INTERIOR_SOURCE_FORCE_SCOPE: &str = "Prescribed stationary compact \
3D force ON the original layer fluid [N], strictly between both faces; complete pressure \
and particle velocity in the incident exterior of one infinite planar homogeneous lossless \
fluid layer. Both internal directions, coherent returns, critical, reactive, transverse \
and evanescent fields are retained. Slower-than-both guided layers refuse. Numerical \
quadrature estimates and analytic tail bounds are not interval certificates; no interior \
receiver, source loading, absorption, moving fluids, finite edges or measured-SPL qualification.";

/// Cached Cartesian response to an original force inside the actual layer.
/// A 4-by-3 transfer retains every source direction and phase; evaluation uses
/// fixed DD accumulation with no allocation, quadrature or extra propagation delay.
#[derive(Clone, Debug)]
pub struct PreparedInteriorSourceFullWaveForce {
    layer: PlanarFluidLayer,
    frequency_hz: f64,
    source_m: [f64; 3],
    observer_m: [f64; 3],
    transfer: [[C64; 3]; 4],
    diagnostics: FullWaveLayerForceDiagnostics,
}

impl PreparedInteriorSourceFullWaveForce {
    /// Source strictly inside; receiver strictly before the first face.
    /// To observe the other exterior, explicitly reverse the same original
    /// layer's exteriors, first-face point and normal. Force remains expressed
    /// in the original Cartesian frame and acts ON its actual source fluid.
    pub fn new(
        layer: PlanarFluidLayer,
        frequency_hz: f64,
        source_m: [f64; 3],
        observer_m: [f64; 3],
        limits: FullWaveForceLimits,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, AeroacError> {
        poll(cancelled)?;
        admit(limits)?;
        let (coordinates, depth, _) =
            interior_coordinates(layer, frequency_hz, observer_m, source_m)?;
        let (transfer, diagnostics) = if reduces_to_interface(layer) {
            let interface = PlanarFluidInterface::new(
                layer.layer_state(),
                layer.incident_state(),
                layer.first_face_point_m(),
                layer.receiving_normal().map(|v| -v),
            )?;
            let original = PreparedFullWaveForce::new(
                interface,
                frequency_hz,
                source_m,
                observer_m,
                limits,
                cancelled,
            )?;
            (original.transfer, (*original.diagnostics()).into())
        } else {
            let mut interior = InteriorKernel::new(layer, frequency_hz, observer_m, source_m)?;
            let r1 = det::ln(layer.incident_state().density_kg_m3());
            let r2 = det::ln(layer.layer_state().density_kg_m3());
            let w = det::ln(coordinates.omega);
            // |A/L|<=2rho1 exp(-alpha2*s)/alpha1; |B/L|<=2 exp(-alpha2*s).
            // The interior force's tangential pressure therefore carries rho1/rho2;
            // its normal pressure has unit density factor. Receiver derivatives
            // give the same max(1/rho1,1/rho2) degree-2 velocity envelope.
            let (extent, bounds) = prepare_reflected_tail(
                (coordinates.a + depth).to_f64(),
                interior.common.critical[interior.common.critical_count - 1],
                limits.into(),
                [1, 2],
                [
                    (r1 - r2).max(0.0) - det::ln(PI),
                    (-r1).max(-r2) - det::ln(PI) - w,
                ],
                [2, 5],
                cancelled,
            )?;
            interior.common.tail_extent = extent;
            let radial = coordinates.range != 0.0;
            let tail = [
                if radial { bounds[0] } else { 0.0 },
                bounds[0],
                bounds[1],
                0.5 * bounds[1],
                if radial { bounds[1] } else { 0.0 },
                if radial { bounds[1] } else { 0.0 },
                bounds[1],
            ];
            let kernel = SourceForceKernel {
                interior,
                log_pressure: -det::ln(TAU),
                log_rho2: r2,
                log_velocity: -det::ln(TAU) - w - r1,
            };
            let result = integrate_spectrum::<7, 2>(&kernel, tail, limits.into(), cancelled)?;
            (
                force_cartesian(
                    result.value,
                    coordinates.tangent.map(|v| -v),
                    layer.receiving_normal(),
                )?,
                FullWaveLayerForceDiagnostics {
                    evaluations: result.evaluations,
                    maximum_depth_reached: result.maximum_depth_reached,
                    retained_panels: result.retained_panels,
                    transverse_wavenumber_cutoff_rad_m: result.cutoff,
                    maximum_bessel_argument: result.maximum_argument,
                    pressure_quadrature_error_pa_per_n: [result.error[0], result.error[1]],
                    velocity_quadrature_error_m_s_per_n: core::array::from_fn(|i| {
                        result.error[i + 2]
                    }),
                    pressure_tail_bound_pa_per_n: [tail[0], tail[1]],
                    velocity_tail_bound_m_s_per_n: core::array::from_fn(|i| tail[i + 2]),
                    reduced_to_single_interface: false,
                },
            )
        };
        poll(cancelled)?;
        Ok(Self {
            layer,
            frequency_hz,
            source_m,
            observer_m,
            transfer,
            diagnostics,
        })
    }

    pub const fn layer(&self) -> PlanarFluidLayer {
        self.layer
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
    pub const fn diagnostics(&self) -> &FullWaveLayerForceDiagnostics {
        &self.diagnostics
    }

    /// Complex peak force ON the layer fluid [N], with exp(-i omega t).
    pub fn evaluate(&self, force: [C64; 3]) -> Result<CoherentCompactField, AeroacError> {
        field(
            apply_force(&self.transfer, force)?,
            FULL_WAVE_INTERIOR_SOURCE_FORCE_SCOPE,
        )
    }

    pub fn evaluate_harmonic_line(
        &self,
        force: [C64; 3],
    ) -> Result<HarmonicFieldLine<3>, AeroacError> {
        Ok(line(self.frequency_hz, apply_force(&self.transfer, force)?))
    }
}

struct SourceForceKernel {
    interior: InteriorKernel,
    log_pressure: f64,
    log_rho2: f64,
    log_velocity: f64,
}

impl Spectrum<7> for SourceForceKernel {
    fn regions(&self) -> usize {
        self.interior.common.critical_count + 1
    }
    fn extent(&self) -> (f64, f64) {
        self.interior.common.extent()
    }
    fn resolved(&self, region: usize, lower: f64, upper: f64) -> bool {
        self.interior.resolved(region, lower, upper)
    }
    fn sample(&self, region: usize, x: f64) -> Result<[C64; 7], AeroacError> {
        let weight = self.interior.sample_weight(region, x)?;
        let spectrum = &weight.spectrum;
        let argument = spectrum.q * self.interior.common.coordinates.range;
        let zero = j0(argument);
        let one = j1(argument);
        let ratio = j1_over_argument(argument, one);
        let q = det::ln(spectrum.q);
        let vertical = spectrum.zeta[0].abs();
        let z = det::ln(vertical);
        let direction = if vertical == 0.0 {
            C64::ZERO
        } else {
            C64::new(
                spectrum.zeta[0].re / vertical,
                spectrum.zeta[0].im / vertical,
            )
        };
        let i = C64::new(0.0, 1.0);
        let a = |unit: C64, scale: f64, value: f64| {
            spectrum.scaled(
                weight.pressure.unit * unit,
                weight.pressure.log_magnitude + scale,
                value,
            )
        };
        let b = |unit: C64, scale: f64, value: f64| {
            spectrum.scaled(
                weight.normal_velocity.unit * unit,
                weight.normal_velocity.log_magnitude + scale,
                value,
            )
        };
        let tangent_velocity = self.log_velocity - self.log_rho2;
        let values = [
            a(i, self.log_pressure - self.log_rho2 + q, one),
            b(C64::from_re(-1.0), self.log_pressure, zero),
            a(C64::ONE, tangent_velocity + 2.0 * q, zero - ratio),
            a(C64::ONE, tangent_velocity + 2.0 * q, ratio),
            b(-i, self.log_velocity + q, one),
            a(-i * direction, tangent_velocity + q + z, one),
            b(direction, self.log_velocity + z, zero),
        ];
        if values.iter().any(|v| !finite(*v)) {
            return Err(nonfinite("interior-source force spectrum"));
        }
        Ok(values)
    }
}
