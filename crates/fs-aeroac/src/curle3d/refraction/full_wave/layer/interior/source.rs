//! Radiation from an original volume source inside a finite fluid layer.
//!
//! For the exp(-i omega t) clock, pressure obeys
//! `div(grad(p)/rho) + omega² p/(rho c²) = i omega Q delta`.
//! Its reciprocal Green function makes p/Q invariant under source/receiver
//! exchange. The exterior velocity is nevertheless a DIFFERENT derivative:
//! it differentiates the original exterior coordinate, giving `-zeta_1 p /
//! (omega rho_1)` in the face-normal direction. Copying the interior receiver's
//! velocity would violate momentum and interface continuity.
//!
//! The original interior kernel supplies the common spectrum, all coherent
//! returns, critical branch maps and scaled numerator. Only the receiving
//! derivatives and their physical tail estimates change here.

use super::super::{
    FullWaveLayerMonopoleDiagnostics, FullWaveMonopoleLimits, PlanarFluidInterface,
    PlanarFluidLayer, PreparedFullWaveMonopole, Spectrum, admit, finite, integrate_spectrum, j0,
    j1, nonfinite, poll,
};
use super::{InteriorKernel, apply, field, interior_coordinates, line};
use crate::AeroacError;
use crate::curle3d::CoherentCompactField;
use crate::observation::HarmonicFieldLine;
use core::f64::consts::TAU;
use fs_math::{c64::C64, det};

/// Coherent original layer volume-flow and force arrays observed in an exterior.
pub mod array;
mod force;
pub use force::{FULL_WAVE_INTERIOR_SOURCE_FORCE_SCOPE, PreparedInteriorSourceFullWaveForce};

/// Applicability of the original interior source and exterior observation.
pub const FULL_WAVE_INTERIOR_SOURCE_SCOPE: &str = "Prescribed stationary compact 3D \
volume source strictly inside one infinite planar homogeneous lossless fluid layer; \
full-wave radiation into the original incident exterior, with coherent internal returns, \
critical, reactive and evanescent fields. Slower-than-both trapping layers refuse. \
Numerical quadrature estimates and analytic spectral-tail bounds are not interval \
certificates; no source loading, absorption, moving fluids, finite edges or measured-SPL qualification.";

/// A physical layer source observed in the half-space before its first face.
/// Preparation is bounded; application reuses four cached complex coefficients.
#[derive(Clone, Debug)]
pub struct PreparedInteriorSourceFullWaveMonopole {
    layer: PlanarFluidLayer,
    frequency_hz: f64,
    source_m: [f64; 3],
    observer_m: [f64; 3],
    transfer: [C64; 4],
    diagnostics: FullWaveLayerMonopoleDiagnostics,
}

impl PreparedInteriorSourceFullWaveMonopole {
    /// The source lies strictly between the original faces and the observer
    /// strictly in the incident exterior. To observe the other exterior,
    /// explicitly supply the same physical layer with its two exteriors and
    /// face orientation reversed. No material identity or side is inferred.
    pub fn new(
        layer: PlanarFluidLayer,
        frequency_hz: f64,
        source_m: [f64; 3],
        observer_m: [f64; 3],
        limits: FullWaveMonopoleLimits,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, AeroacError> {
        poll(cancelled)?;
        admit(limits)?;
        interior_coordinates(layer, frequency_hz, observer_m, source_m)?;
        let (transfer, diagnostics) = if layer.layer_state() == layer.receiving_state() {
            // The second face is physically absent. Keep the original exact
            // single-interface limits, including reactive spherical velocity.
            let interface = PlanarFluidInterface::new(
                layer.layer_state(),
                layer.incident_state(),
                layer.first_face_point_m(),
                layer.receiving_normal().map(|x| -x),
            )?;
            let original = PreparedFullWaveMonopole::new(
                interface,
                frequency_hz,
                source_m,
                observer_m,
                limits,
                cancelled,
            )?;
            (original.transfer, (*original.diagnostics()).into())
        } else {
            let mut kernel = ExteriorKernel {
                interior: InteriorKernel::new(layer, frequency_hz, observer_m, source_m)?,
            };
            let tail = kernel.interior.prepare_tail(
                limits,
                layer.incident_state().density_kg_m3(),
                cancelled,
            )?;
            let result = integrate_spectrum::<3, 1>(&kernel, tail, limits.into(), cancelled)?;
            let [pressure, radial, normal] = result.value;
            let tangent = kernel.interior.common.coordinates.tangent;
            let transfer = core::array::from_fn(|i| {
                if i == 0 {
                    pressure
                } else {
                    // The shared coordinate is exterior -> interior; the
                    // physical receiving radial direction points the other way.
                    radial.scale(-tangent[i - 1]) + normal.scale(layer.receiving_normal()[i - 1])
                }
            });
            (
                transfer,
                FullWaveLayerMonopoleDiagnostics {
                    evaluations: result.evaluations,
                    maximum_depth_reached: result.maximum_depth_reached,
                    retained_panels: result.retained_panels,
                    transverse_wavenumber_cutoff_rad_m: result.cutoff,
                    maximum_bessel_argument: result.maximum_argument,
                    pressure_quadrature_error_pa_per_m3_s: result.error[0],
                    velocity_quadrature_error_per_m2: [result.error[1], result.error[2]],
                    pressure_tail_bound_pa_per_m3_s: tail[0],
                    velocity_tail_bound_per_m2: [tail[1], tail[2]],
                    reduced_to_single_interface: false,
                },
            )
        };
        if transfer.iter().any(|v| !finite(*v)) {
            return Err(nonfinite("interior-source exterior transfer"));
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
    pub const fn diagnostics(&self) -> &FullWaveLayerMonopoleDiagnostics {
        &self.diagnostics
    }

    /// Complex peak Q [m³/s] on the original source clock; no extra ray delay.
    pub fn evaluate(&self, volume_flow_m3_s: C64) -> Result<CoherentCompactField, AeroacError> {
        field(
            apply(self.transfer, volume_flow_m3_s)?,
            FULL_WAVE_INTERIOR_SOURCE_SCOPE,
        )
    }

    pub fn evaluate_harmonic_line(
        &self,
        volume_flow_m3_s: C64,
    ) -> Result<HarmonicFieldLine<3>, AeroacError> {
        Ok(line(
            self.frequency_hz,
            apply(self.transfer, volume_flow_m3_s)?,
        ))
    }
}

struct ExteriorKernel {
    interior: InteriorKernel,
}

impl Spectrum<3> for ExteriorKernel {
    fn regions(&self) -> usize {
        self.interior.common.critical_count + 1
    }
    fn extent(&self) -> (f64, f64) {
        self.interior.common.extent()
    }
    fn resolved(&self, region: usize, lower: f64, upper: f64) -> bool {
        self.interior.resolved(region, lower, upper)
    }
    fn sample(&self, region: usize, x: f64) -> Result<[C64; 3], AeroacError> {
        let weight = self.interior.sample_weight(region, x)?;
        let spectrum = &weight.spectrum;
        let numerator = &weight.pressure;
        let common = &self.interior.common;
        let argument = spectrum.q * common.coordinates.range;
        let zero = j0(argument);
        let log_velocity = numerator.log_magnitude
            - det::ln(TAU)
            - det::ln(common.layer.incident_state().density_kg_m3());
        let vertical = spectrum.zeta[0].abs();
        let direction = if vertical == 0.0 {
            C64::ZERO
        } else {
            C64::new(
                -spectrum.zeta[0].re / vertical,
                -spectrum.zeta[0].im / vertical,
            )
        };
        let values = [
            spectrum.scaled(
                numerator.unit,
                numerator.log_magnitude + det::ln(common.coordinates.omega) - det::ln(TAU),
                zero,
            ),
            spectrum.scaled(
                C64::new(-numerator.unit.im, numerator.unit.re),
                log_velocity + det::ln(spectrum.q),
                j1(argument),
            ),
            spectrum.scaled(
                direction * numerator.unit,
                log_velocity + det::ln(vertical),
                zero,
            ),
        ];
        if values.iter().any(|v| !finite(*v)) {
            return Err(nonfinite("interior-source spectral field"));
        }
        Ok(values)
    }
}
