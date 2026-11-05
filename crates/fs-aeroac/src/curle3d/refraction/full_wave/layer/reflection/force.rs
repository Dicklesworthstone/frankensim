//! Original force-source derivatives of the complete finite-layer reflection.
//!
//! The parent retains the full reflected Green weight. Its two pressure and
//! five velocity derivatives use exactly the shared single-interface force
//! expressions, including opposite mixed signs and normal image-force reversal.

use super::super::super::force::admit;
use super::super::super::reflection::{
    PreparedReflectedFullWaveForce, ReflectionCoordinates, apply_force, field, force_cartesian,
    line, sample_force_spectrum, singular_total,
};
use super::super::super::{FullWaveForceLimits, Spectrum, integrate_spectrum, poll};
use super::super::{FullWaveLayerForceDiagnostics, PlanarFluidLayer};
use super::{LayerReflectionKernel, add_direct, reduces_to_first_interface};
use crate::AeroacError;
use crate::curle3d::{CoherentCompactField, CompactSourcePhasor, PreparedCompactSourceObserver};
use crate::observation::HarmonicFieldLine;
use core::f64::consts::PI;
use fs_math::{c64::C64, det};

pub const FULL_WAVE_LAYER_REFLECTED_FORCE_SCOPE: &str = "Prescribed stationary compact \
3D force ON fluid [N]; reflected spherical-wave field from one infinite planar \
homogeneous lossless fluid layer, including coherent internal returns, transverse \
and reactive motion and evanescence. Slower-than-both guided layers refuse. Estimated \
quadrature error and analytic tail bounds are not interval certificates; no source \
feedback, absorption, moving fluids, finite edges or measured-SPL qualification.";

pub const FULL_WAVE_LAYER_TOTAL_FORCE_SCOPE: &str = "Prescribed stationary compact \
3D force ON fluid [N]; coherent direct plus finite-fluid-layer reflected spherical \
field in the incident fluid. Reflected-transfer error estimates give no relative \
guarantee at interference nulls. Slower-than-both guided layers and singular direct \
self-fields refuse; no source feedback, absorption, moving fluids, finite edges, \
interval certificate or measured-SPL qualification.";

type Transfer = [[C64; 3]; 4];

/// Two cached 4-by-3 Cartesian transfers retain all original force directions
/// and phases. Applying them performs no allocation, quadrature, finite
/// differences or transcendental work. Total self-fields are undefined.
#[derive(Clone, Debug)]
pub struct PreparedLayeredReflectedFullWaveForce {
    layer: PlanarFluidLayer,
    frequency_hz: f64,
    source_m: [f64; 3],
    observer_m: [f64; 3],
    reflected: Transfer,
    total: Option<Transfer>,
    diagnostics: FullWaveLayerForceDiagnostics,
}

impl PreparedLayeredReflectedFullWaveForce {
    /// Prepare an original force and receiver strictly before the first face.
    /// Pressure-vector and velocity-tensor norm tolerances control the reflected
    /// transfer per newton; all seven channels share bounded spectral nodes.
    /// The total inherits absolute reflected estimates, not relative accuracy
    /// near a direct/returned cancellation. Failure publishes no partial field.
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
        let coordinates =
            ReflectionCoordinates::new(layer.exterior, frequency_hz, source_m, observer_m)?;
        let (reflected, diagnostics) = if reduces_to_first_interface(layer) {
            let original = PreparedReflectedFullWaveForce::new(
                layer.exterior,
                frequency_hz,
                source_m,
                observer_m,
                limits,
                cancelled,
            )?;
            let mut transfer = [[C64::ZERO; 3]; 4];
            for column in 0..3 {
                let force =
                    core::array::from_fn(|i| if i == column { C64::ONE } else { C64::ZERO });
                let original = original.evaluate_reflected_harmonic_line(force)?;
                transfer[0][column] = original.pressure_pa;
                for row in 0..3 {
                    transfer[row + 1][column] = original.particle_velocity_m_s[row];
                }
            }
            (transfer, (*original.diagnostics()).into())
        } else {
            let log_pressure = -det::ln(4.0 * PI);
            let log_velocity = log_pressure
                - det::ln(coordinates.omega)
                - det::ln(layer.incident_state().density_kg_m3());
            let mut common = LayerReflectionKernel::new(layer, coordinates)?;
            let bound = common.prepare_tail(
                limits.into(),
                [1, 2],
                [log_pressure, log_velocity],
                [2, 5],
                cancelled,
            )?;
            let radial = coordinates.range != 0.0;
            let tail = [
                if radial { bound[0] } else { 0.0 },
                bound[0],
                bound[1],
                0.5 * bound[1],
                if radial { bound[1] } else { 0.0 },
                if radial { bound[1] } else { 0.0 },
                bound[1],
            ];
            let kernel = LayerReflectedForceKernel {
                common,
                log_pressure,
                log_velocity,
            };
            let result = integrate_spectrum::<7, 2>(&kernel, tail, limits.into(), cancelled)?;
            let transfer =
                force_cartesian(result.value, coordinates.tangent, layer.receiving_normal())?;
            (
                transfer,
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
        let total = if source_m == observer_m {
            None
        } else {
            let direct = PreparedCompactSourceObserver::new(
                layer.incident_state(),
                frequency_hz,
                observer_m,
                &[source_m],
                1,
            )?;
            let mut transfer = [[C64::ZERO; 3]; 4];
            for column in 0..3 {
                poll(cancelled)?;
                let force =
                    core::array::from_fn(|i| if i == column { C64::ONE } else { C64::ZERO });
                let direct =
                    direct.evaluate_harmonic_line(&[CompactSourcePhasor::force_dipole(force)])?;
                let direct = [
                    direct.pressure_pa,
                    direct.particle_velocity_m_s[0],
                    direct.particle_velocity_m_s[1],
                    direct.particle_velocity_m_s[2],
                ];
                let summed =
                    add_direct(direct, core::array::from_fn(|row| reflected[row][column]))?;
                for row in 0..4 {
                    transfer[row][column] = summed[row];
                }
            }
            Some(transfer)
        };
        poll(cancelled)?;
        Ok(Self {
            layer,
            frequency_hz,
            source_m,
            observer_m,
            reflected,
            total,
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
    pub const fn has_total_field(&self) -> bool {
        self.total.is_some()
    }
    /// Reflected pressure-vector and velocity-tensor error estimates per N.
    pub const fn diagnostics(&self) -> &FullWaveLayerForceDiagnostics {
        &self.diagnostics
    }

    /// Complex peak physical force ON the original fluid [N], exp(-i omega t).
    pub fn evaluate_reflected(
        &self,
        force_on_fluid_n: [C64; 3],
    ) -> Result<CoherentCompactField, AeroacError> {
        field(
            apply_force(&self.reflected, force_on_fluid_n)?,
            FULL_WAVE_LAYER_REFLECTED_FORCE_SCOPE,
        )
    }
    pub fn evaluate_total(
        &self,
        force_on_fluid_n: [C64; 3],
    ) -> Result<CoherentCompactField, AeroacError> {
        field(
            apply_force(
                self.total.as_ref().ok_or_else(singular_total)?,
                force_on_fluid_n,
            )?,
            FULL_WAVE_LAYER_TOTAL_FORCE_SCOPE,
        )
    }
    pub fn evaluate_reflected_harmonic_line(
        &self,
        force_on_fluid_n: [C64; 3],
    ) -> Result<HarmonicFieldLine<3>, AeroacError> {
        Ok(line(
            self.frequency_hz,
            apply_force(&self.reflected, force_on_fluid_n)?,
        ))
    }
    pub fn evaluate_total_harmonic_line(
        &self,
        force_on_fluid_n: [C64; 3],
    ) -> Result<HarmonicFieldLine<3>, AeroacError> {
        Ok(line(
            self.frequency_hz,
            apply_force(
                self.total.as_ref().ok_or_else(singular_total)?,
                force_on_fluid_n,
            )?,
        ))
    }
}

struct LayerReflectedForceKernel {
    common: LayerReflectionKernel,
    log_pressure: f64,
    log_velocity: f64,
}

impl Spectrum<7> for LayerReflectedForceKernel {
    fn regions(&self) -> usize {
        <LayerReflectionKernel as Spectrum<3>>::regions(&self.common)
    }
    fn extent(&self) -> (f64, f64) {
        <LayerReflectionKernel as Spectrum<3>>::extent(&self.common)
    }
    fn resolved(&self, region: usize, lower: f64, upper: f64) -> bool {
        <LayerReflectionKernel as Spectrum<3>>::resolved(&self.common, region, lower, upper)
    }
    fn sample(&self, region: usize, x: f64) -> Result<[C64; 7], AeroacError> {
        sample_force_spectrum(
            self.common.sample_weight(region, x)?,
            self.common.common.coordinates.range,
            self.log_pressure,
            self.log_velocity,
        )
    }
}
