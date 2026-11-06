//! Original volume stresses inside a finite layer, radiating into an exterior.
//!
//! The reciprocal interior kernel supplies pressure/normal weights WA/WB.
//! H=i/(2*pi*rho2) integral WA*J0 and p=M:d_s d_s H. Original interior
//! source differentiation gives A_s=i*rho2*B and A_ss=-zeta2^2*A. Thus a
//! mixed normal/tangential source derivative uses B, while two normal source
//! derivatives use the ORIGINAL layer wavenumber, never the exterior one.
//! Exterior receiver momentum uses rho1 and d_rn=-i*zeta1; the physical
//! receiving tangent is the negative of the reciprocal kernel's tangent.
//! All fourteen channels retain both internal directions and coherent returns.

use super::super::{InteriorKernel, interior_coordinates, reduces_to_interface};
use crate::AeroacError;
use crate::curle3d::CoherentCompactField;
use crate::curle3d::quadrupole::QuadrupolePhasor;
use crate::curle3d::refraction::PlanarFluidInterface;
use crate::curle3d::refraction::full_wave::force::j1_over_argument;
use crate::curle3d::refraction::full_wave::layer::{
    FullWaveLayerQuadrupoleDiagnostics, PlanarFluidLayer,
};
use crate::curle3d::refraction::full_wave::quadrupole::{
    FullWaveQuadrupoleLimits, PreparedFullWaveQuadrupole, cartesian, j2_over_argument,
};
use crate::curle3d::refraction::full_wave::reflection::{
    apply_quadrupole, field, line, prepare_reflected_tail,
};
use crate::curle3d::refraction::full_wave::{
    QuadratureLimits, Spectrum, finite, integrate_spectrum, j0, j1, nonfinite, poll,
};
use crate::observation::HarmonicFieldLine;
use core::f64::consts::{PI, SQRT_2, TAU};
use fs_math::{c64::C64, det};

/// Applicability of the original physical layer stress and exterior field.
pub const FULL_WAVE_INTERIOR_SOURCE_QUADRUPOLE_SCOPE: &str = "Prescribed stationary \
symmetric 3D Lighthill stress moment [N m] strictly inside the original finite \
fluid layer; full spherical-wave pressure and vector velocity in its incident \
exterior. Both internal directions, all coherent returns, physical trace/shear, \
near, transverse, critical and evanescent fields are retained. Infinite planar \
homogeneous lossless nontrapping fluids only. Numerical quadrature estimates and \
analytic tails are not interval certificates. No interior receiver, source-size \
or SGS validation, absorption, moving media, finite edges, source loading or \
measured-SPL qualification.";

/// Complete cached Cartesian response to an original interior stress moment.
/// The 4-by-6 map retains original tensor entries and source phase. Warm
/// application uses the shared fixed-order DD contraction with no allocation,
/// quadrature, finite differences or additional propagation delay.
#[derive(Clone, Debug)]
pub struct PreparedInteriorSourceFullWaveQuadrupole {
    layer: PlanarFluidLayer,
    frequency_hz: f64,
    source_m: [f64; 3],
    observer_m: [f64; 3],
    transfer: [[C64; 6]; 4],
    diagnostics: FullWaveLayerQuadrupoleDiagnostics,
}

impl PreparedInteriorSourceFullWaveQuadrupole {
    /// Source strictly between both original faces; observer strictly before
    /// the first face. To observe the other exterior, explicitly reverse the
    /// same physical layer's exteriors, first-face point and normal. The stress
    /// remains in the original Cartesian frame and original source fluid.
    ///
    /// All channels share native branch maps, the coherent layer denominator,
    /// resonance resolution and bounded quadrature. Tolerances control the
    /// pressure-Hessian and velocity-tensor Frobenius norms per N m, before
    /// applying a particular source tensor. Failure publishes no partial field.
    pub fn new(
        layer: PlanarFluidLayer,
        frequency_hz: f64,
        source_m: [f64; 3],
        observer_m: [f64; 3],
        limits: FullWaveQuadrupoleLimits,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, AeroacError> {
        poll(cancelled)?;
        QuadratureLimits::from(limits).admit()?;
        let (coordinates, depth, _) =
            interior_coordinates(layer, frequency_hz, observer_m, source_m)?;
        let (transfer, diagnostics) = if reduces_to_interface(layer) {
            // An absent original second face leaves the physical layer source
            // radiating across the first face into the original exterior.
            let interface = PlanarFluidInterface::new(
                layer.layer_state(),
                layer.incident_state(),
                layer.first_face_point_m(),
                layer.receiving_normal().map(|v| -v),
            )?;
            let original = PreparedFullWaveQuadrupole::new(
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
            // Above all critical circles, |A/L|<=2rho1 e^(-alpha2*s)/alpha1
            // and |B/L|<=2e^(-alpha2*s). Put q*dq=t*dt and alpha_j>=t.
            // Pressure has degree 2 with max(rho1/rho2,1); receiving velocity
            // has degree 3 with max(1/rho1,1/rho2)/omega. Both use the actual
            // exterior clearance plus source depth, not the full layer width.
            let (extent, [p, v]) = prepare_reflected_tail(
                (coordinates.a + depth).to_f64(),
                interior.common.critical[interior.common.critical_count - 1],
                limits.into(),
                [2, 3],
                [
                    (r1 - r2).max(0.0) - det::ln(PI),
                    (-r1).max(-r2) - det::ln(PI) - w,
                ],
                [5, 14],
                cancelled,
            )?;
            interior.common.tail_extent = extent;
            let odd_p = if coordinates.range == 0.0 { 0.0 } else { p };
            let odd_v = if coordinates.range == 0.0 { 0.0 } else { v };
            let tail = [
                p,
                p,
                SQRT_2 * odd_p,
                p,
                odd_v,
                odd_v,
                SQRT_2 * v,
                odd_v,
                SQRT_2 * odd_v,
                SQRT_2 * v,
                v,
                v,
                SQRT_2 * odd_v,
                v,
            ];
            let kernel = SourceQuadrupoleKernel {
                interior,
                log_pressure: -det::ln(TAU),
                log_rho2: r2,
                log_velocity: -det::ln(TAU) - w - r1,
            };
            let result = integrate_spectrum::<14, 4>(&kernel, tail, limits.into(), cancelled)?;
            (
                cartesian(
                    result.value,
                    coordinates.tangent.map(|v| -v),
                    layer.receiving_normal(),
                )?,
                FullWaveLayerQuadrupoleDiagnostics {
                    evaluations: result.evaluations,
                    maximum_depth_reached: result.maximum_depth_reached,
                    retained_panels: result.retained_panels,
                    transverse_wavenumber_cutoff_rad_m: result.cutoff,
                    maximum_bessel_argument: result.maximum_argument,
                    pressure_quadrature_error_pa_per_n_m: core::array::from_fn(|i| result.error[i]),
                    velocity_quadrature_error_m_s_per_n_m: core::array::from_fn(|i| {
                        result.error[i + 4]
                    }),
                    pressure_tail_bound_pa_per_n_m: core::array::from_fn(|i| tail[i]),
                    velocity_tail_bound_m_s_per_n_m: core::array::from_fn(|i| tail[i + 4]),
                    reduced_to_single_interface: false,
                },
            )
        };
        if transfer.iter().flatten().any(|value| !finite(*value)) {
            return Err(nonfinite("interior-source quadrupole Cartesian transfer"));
        }
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
    /// Four pressure and ten velocity Frobenius channel estimates. The
    /// reduction flag denotes an absent second face (layer equals rear fluid).
    pub const fn diagnostics(&self) -> &FullWaveLayerQuadrupoleDiagnostics {
        &self.diagnostics
    }

    /// Original symmetric peak stress [N m], with exp(-i omega t).
    pub fn evaluate(&self, source: QuadrupolePhasor) -> Result<CoherentCompactField, AeroacError> {
        field(
            apply_quadrupole(&self.transfer, source)?,
            FULL_WAVE_INTERIOR_SOURCE_QUADRUPOLE_SCOPE,
        )
    }

    pub fn evaluate_harmonic_line(
        &self,
        source: QuadrupolePhasor,
    ) -> Result<HarmonicFieldLine<3>, AeroacError> {
        Ok(line(
            self.frequency_hz,
            apply_quadrupole(&self.transfer, source)?,
        ))
    }
}

struct SourceQuadrupoleKernel {
    interior: InteriorKernel,
    log_pressure: f64,
    log_rho2: f64,
    log_velocity: f64,
}

impl Spectrum<14> for SourceQuadrupoleKernel {
    fn regions(&self) -> usize {
        <InteriorKernel as Spectrum<3>>::regions(&self.interior)
    }
    fn extent(&self) -> (f64, f64) {
        <InteriorKernel as Spectrum<3>>::extent(&self.interior)
    }
    fn resolved(&self, region: usize, lower: f64, upper: f64) -> bool {
        <InteriorKernel as Spectrum<3>>::resolved(&self.interior, region, lower, upper)
    }
    fn sample(&self, region: usize, x: f64) -> Result<[C64; 14], AeroacError> {
        let weight = self.interior.sample_weight(region, x)?;
        let spectrum = &weight.spectrum;
        let argument = spectrum.q * self.interior.common.coordinates.range;
        if !argument.is_finite() {
            return Err(nonfinite("interior-source quadrupole Bessel argument"));
        }
        let zero = j0(argument);
        let one = j1(argument);
        let ratio = j1_over_argument(argument, one);
        let derivative = zero - ratio;
        let second_ratio = j2_over_argument(argument, zero, ratio);
        let third = one - second_ratio;
        let q = det::ln(spectrum.q);
        let magnitude = [spectrum.zeta[0].abs(), spectrum.zeta[1].abs()];
        let z = magnitude.map(det::ln);
        let direction = core::array::from_fn::<_, 2, _>(|index| {
            if magnitude[index] == 0.0 {
                C64::ZERO
            } else {
                C64::new(
                    spectrum.zeta[index].re / magnitude[index],
                    spectrum.zeta[index].im / magnitude[index],
                )
            }
        });
        let square = direction[1] * direction[1];
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
        // Keep independent A/B zeros and all derivative magnitudes balanced
        // with native attenuation; no raw q^3 or zeta1*zeta2^2 is formed.
        let p = self.log_pressure - self.log_rho2;
        let v = self.log_velocity - self.log_rho2;
        let values = [
            a(-i, p + 2.0 * q, derivative),
            a(-i, p + 2.0 * q, ratio),
            b(-C64::ONE, self.log_pressure + q, SQRT_2 * one),
            a(-i * square, p + 2.0 * z[1], zero),
            a(C64::ONE, v + 3.0 * q, third),
            a(C64::ONE, v + 3.0 * q, second_ratio),
            b(i, self.log_velocity + 2.0 * q, SQRT_2 * derivative),
            a(square, v + q + 2.0 * z[1], one),
            a(C64::ONE, v + 3.0 * q, SQRT_2 * second_ratio),
            b(i, self.log_velocity + 2.0 * q, SQRT_2 * ratio),
            a(i * direction[0], v + 2.0 * q + z[0], derivative),
            a(i * direction[0], v + 2.0 * q + z[0], ratio),
            b(direction[0], self.log_velocity + q + z[0], SQRT_2 * one),
            a(i * direction[0] * square, v + z[0] + 2.0 * z[1], zero),
        ];
        if values.iter().any(|value| !finite(*value)) {
            return Err(nonfinite("interior-source quadrupole spectrum"));
        }
        Ok(values)
    }
}
