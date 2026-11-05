//! Original source-force derivatives of the internal Green function.
//!
//! For exp(-i omega t), pF=+i/(omega*rho_layer) grad_source(pQ).
//! Each native returning family has source-normal derivative i*sigma*zeta2,
//! where sigma=+1 at the first face and -1 at the second. Reuse its original
//! pressure numerator C and receiving-velocity numerator V. With the parent's
//! common weight K, the pressure channels are K*[i*q*C*J1,-sigma*zeta2*C*J0]
//! /(4*pi), and receiving momentum gives the five velocity channels
//! K*[q²*C*(J0-J1/x),q²*C*J1/x,-i*sigma*q*zeta2*C*J1,
//! i*sigma*q*rho2*V*J1,-zeta2*rho2*V*J0]/(4*pi*omega*rho2).
//! The mixed channels have distinct source/receiver signs; the normal-normal
//! channel has no remaining face sign. No finite-difference source is used.

use super::{FaceKernel, FaceWeight, InternalCoordinates, LayerKernel, unit};
use crate::AeroacError;
use crate::curle3d::refraction::full_wave::force::{admit, j1_over_argument};
use crate::curle3d::refraction::full_wave::layer::{
    FullWaveLayerForceDiagnostics, PlanarFluidLayer,
};
use crate::curle3d::refraction::full_wave::reflection::{
    PreparedReflectedFullWaveForce, apply_force, field, force_cartesian, line,
    prepare_reflected_tail, singular_total,
};
use crate::curle3d::refraction::full_wave::{
    FullWaveForceLimits, Spectrum, integrate_spectrum, j0, j1, poll,
};
use crate::curle3d::refraction::{PlanarFluidInterface, positive};
use crate::curle3d::statistics::{ComplexSum, finite, invalid, nonfinite};
use crate::curle3d::{CoherentCompactField, CompactSourcePhasor, PreparedCompactSourceObserver};
use crate::observation::HarmonicFieldLine;
use core::f64::consts::PI;
use fs_math::{c64::C64, det};

pub const FULL_WAVE_LAYER_INTERNAL_REFLECTED_FORCE_SCOPE: &str = "Prescribed stationary compact \
3D force ON the original layer fluid [N] and receiver strictly inside the same infinite \
planar homogeneous lossless layer; complete coherent boundary-return pressure and actual \
velocity, including transverse, reactive, critical and evanescent motion. The analytic \
direct field is excluded. Slower-than-both guided layers refuse. Numerical quadrature \
estimates and analytic tail bounds are not interval certificates; no source loading, \
absorption, moving fluids, finite edges or measured-SPL qualification.";

pub const FULL_WAVE_LAYER_INTERNAL_TOTAL_FORCE_SCOPE: &str = "Prescribed stationary compact \
3D force ON the original layer fluid [N] and receiver strictly inside the same infinite \
planar homogeneous lossless layer; coherent analytic direct plus complete boundary-return \
pressure and vector velocity. Reflected-transfer errors do not guarantee relative accuracy \
at interference nulls. Coincident direct fields and slower-than-both guided layers refuse; \
no source loading, absorption, moving fluids, finite edges, interval certificate or measured-SPL qualification.";

type Transfer = [[C64; 3]; 4];

/// Cached boundary and coherent total Cartesian responses. Arbitrary complex
/// physical forces apply with fixed DD accumulation and no further quadrature.
#[derive(Clone, Debug)]
pub struct PreparedLayerInternalFullWaveForce {
    layer: PlanarFluidLayer,
    frequency_hz: f64,
    source_m: [f64; 3],
    observer_m: [f64; 3],
    reflected: Transfer,
    total: Option<Transfer>,
    diagnostics: FullWaveLayerForceDiagnostics,
}

impl PreparedLayerInternalFullWaveForce {
    /// Source and receiver lie strictly between the original faces. The seven
    /// boundary-return channels share bounded adaptive samples per returning
    /// family. Active families divide the absolute tolerances, use zero relative
    /// tolerance and share one evaluation budget. The singular direct self-field
    /// is absent from the optional total at coincident points. An absent second
    /// face delegates to the original reflected-force owner with the original
    /// limits, including its unequal-speed and exact equal-speed cases.
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
        let geometry = InternalCoordinates::new(layer, frequency_hz, source_m, observer_m)?;
        if !positive(geometry.spectral.omega / layer.layer_state().sound_speed_m_s()) {
            return Err(nonfinite("internal force wavenumber"));
        }
        let (reflected, diagnostics) = if layer.layer_state() == layer.receiving_state() {
            let interface = PlanarFluidInterface::new(
                layer.layer_state(),
                layer.incident_state(),
                layer.first_face_point_m(),
                layer.receiving_normal().map(|v| -v),
            )?;
            let original = PreparedReflectedFullWaveForce::new(
                interface,
                frequency_hz,
                source_m,
                observer_m,
                limits,
                cancelled,
            )?;
            let mut reflected = [[C64::ZERO; 3]; 4];
            for column in 0..3 {
                poll(cancelled)?;
                let force =
                    core::array::from_fn(|i| if i == column { C64::ONE } else { C64::ZERO });
                let value = original.evaluate_reflected_harmonic_line(force)?;
                reflected[0][column] = value.pressure_pa;
                for axis in 0..3 {
                    reflected[axis + 1][column] = value.particle_velocity_m_s[axis];
                }
            }
            (reflected, (*original.diagnostics()).into())
        } else {
            reflected_transfer(layer, geometry, limits, cancelled)?
        };
        let total = if source_m == observer_m {
            None
        } else {
            let direct = PreparedCompactSourceObserver::new(
                layer.layer_state(),
                frequency_hz,
                observer_m,
                &[source_m],
                1,
            )?;
            let mut transfer = [[C64::ZERO; 3]; 4];
            for column in 0..3 {
                poll(cancelled)?;
                let unit = core::array::from_fn(|i| if i == column { C64::ONE } else { C64::ZERO });
                let response =
                    direct.evaluate_harmonic_line(&[CompactSourcePhasor::force_dipole(unit)])?;
                let direct = [
                    response.pressure_pa,
                    response.particle_velocity_m_s[0],
                    response.particle_velocity_m_s[1],
                    response.particle_velocity_m_s[2],
                ];
                for row in 0..4 {
                    let mut sum = ComplexSum::default();
                    sum.add(direct[row])?;
                    sum.add(reflected[row][column])?;
                    transfer[row][column] = sum.value();
                }
            }
            Some(transfer)
        };
        if reflected
            .iter()
            .flatten()
            .chain(total.iter().flatten().flatten())
            .any(|v| !finite(*v))
        {
            return Err(nonfinite("internal layer force transfer"));
        }
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
    /// Summed absolute estimates and omitted-tail bounds per unit force. Work
    /// and panels add across families; depth, cutoff and Bessel argument are
    /// their maxima. Direct-field arithmetic is outside these diagnostics.
    pub const fn diagnostics(&self) -> &FullWaveLayerForceDiagnostics {
        &self.diagnostics
    }
    pub const fn total_field_available(&self) -> bool {
        self.total.is_some()
    }
    pub fn evaluate_reflected(&self, force: [C64; 3]) -> Result<CoherentCompactField, AeroacError> {
        field(
            apply_force(&self.reflected, force)?,
            FULL_WAVE_LAYER_INTERNAL_REFLECTED_FORCE_SCOPE,
        )
    }
    pub fn evaluate_total(&self, force: [C64; 3]) -> Result<CoherentCompactField, AeroacError> {
        field(
            apply_force(self.total.as_ref().ok_or_else(singular_total)?, force)?,
            FULL_WAVE_LAYER_INTERNAL_TOTAL_FORCE_SCOPE,
        )
    }
    pub fn evaluate_reflected_harmonic_line(
        &self,
        force: [C64; 3],
    ) -> Result<HarmonicFieldLine<3>, AeroacError> {
        Ok(line(
            self.frequency_hz,
            apply_force(&self.reflected, force)?,
        ))
    }
    pub fn evaluate_total_harmonic_line(
        &self,
        force: [C64; 3],
    ) -> Result<HarmonicFieldLine<3>, AeroacError> {
        Ok(line(
            self.frequency_hz,
            apply_force(self.total.as_ref().ok_or_else(singular_total)?, force)?,
        ))
    }
}

fn reflected_transfer(
    layer: PlanarFluidLayer,
    geometry: InternalCoordinates,
    limits: FullWaveForceLimits,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(Transfer, FullWaveLayerForceDiagnostics), AeroacError> {
    let media = [
        layer.incident_state(),
        layer.layer_state(),
        layer.receiving_state(),
    ];
    let active = [media[0] != media[1], media[2] != media[1]];
    let count = active.iter().filter(|v| **v).count();
    let mut part_limits = FullWaveForceLimits {
        absolute_pressure_tolerance_pa_per_n: limits.absolute_pressure_tolerance_pa_per_n
            / count as f64,
        absolute_velocity_tolerance_m_s_per_n: limits.absolute_velocity_tolerance_m_s_per_n
            / count as f64,
        relative_tolerance: 0.0,
        ..limits
    };
    admit(part_limits)?;
    let mut common = LayerKernel::new(layer, geometry.spectral)?;
    let rho_log = media.map(|fluid| det::ln(fluid.density_kg_m3()));
    let log_pressure = -det::ln(4.0 * PI);
    let log_velocity = log_pressure - det::ln(geometry.spectral.omega) - rho_log[1];
    let radial = geometry.spectral.range != 0.0;
    let mut sums = [ComplexSum::default(); 7];
    let mut diagnostics = FullWaveLayerForceDiagnostics {
        evaluations: 0,
        maximum_depth_reached: 0,
        retained_panels: 0,
        transverse_wavenumber_cutoff_rad_m: 0.0,
        maximum_bessel_argument: 0.0,
        pressure_quadrature_error_pa_per_n: [0.0; 2],
        velocity_quadrature_error_m_s_per_n: [0.0; 5],
        pressure_tail_bound_pa_per_n: [0.0; 2],
        velocity_tail_bound_m_s_per_n: [0.0; 5],
        reduced_to_single_interface: false,
    };
    for face in 0..2 {
        if !active[face] {
            continue;
        }
        poll(cancelled)?;
        if diagnostics.evaluations == limits.maximum_evaluations {
            return Err(invalid(
                "internal layer force exhausted its shared evaluation budget",
            ));
        }
        part_limits.maximum_evaluations = limits.maximum_evaluations - diagnostics.evaluations;
        let return_depth = geometry.source[face] + geometry.observer[face];
        let depth = return_depth.to_f64();
        if !positive(depth) {
            return Err(nonfinite("internal layer force return clearance"));
        }
        // Both the pressure bracket and its normalized receiver derivative
        // are <=4*exp(-t*depth) once t*h>=1/2, t²=q²-kmax². Source/receiver
        // derivatives add one/two powers, each bounded by t+kmax. Keep each
        // face's actual return length instead of imposing the shorter tail
        // on a distant return. The enlarged cutoff preserves returned bounds.
        let (extent, bound) = prepare_reflected_tail(
            depth,
            common.critical[common.critical_count - 1],
            part_limits.into(),
            [1, 2],
            [log_pressure + det::ln(4.0), log_velocity + det::ln(4.0)],
            [2, 5],
            cancelled,
        )?;
        common.tail_extent = extent.max(0.5 / layer.thickness_m());
        if !positive(common.tail_extent) {
            return Err(nonfinite("internal layer force tail cutoff"));
        }
        let tail = [
            if radial { bound[0] } else { 0.0 },
            bound[0],
            bound[1],
            0.5 * bound[1],
            if radial { bound[1] } else { 0.0 },
            if radial { bound[1] } else { 0.0 },
            bound[1],
        ];
        let kernel = FaceForceKernel {
            common: FaceKernel {
                common: &common,
                face,
                source_depth: geometry.source[face],
                remaining_depth: geometry.observer[1 - face],
                return_depth,
                rho_log,
                log_pressure: log_pressure + det::ln(geometry.spectral.omega) + rho_log[1],
                log_velocity: log_pressure,
            },
            log_pressure,
            log_velocity,
        };
        let result = integrate_spectrum::<7, 2>(&kernel, tail, part_limits.into(), cancelled)?;
        for (sum, value) in sums.iter_mut().zip(result.value) {
            sum.add(value)?;
        }
        diagnostics.evaluations += result.evaluations;
        diagnostics.maximum_depth_reached = diagnostics
            .maximum_depth_reached
            .max(result.maximum_depth_reached);
        diagnostics.retained_panels += result.retained_panels;
        diagnostics.transverse_wavenumber_cutoff_rad_m = diagnostics
            .transverse_wavenumber_cutoff_rad_m
            .max(result.cutoff);
        diagnostics.maximum_bessel_argument = diagnostics
            .maximum_bessel_argument
            .max(result.maximum_argument);
        for axis in 0..2 {
            diagnostics.pressure_quadrature_error_pa_per_n[axis] += result.error[axis];
            diagnostics.pressure_tail_bound_pa_per_n[axis] += tail[axis];
        }
        for axis in 0..5 {
            diagnostics.velocity_quadrature_error_m_s_per_n[axis] += result.error[axis + 2];
            diagnostics.velocity_tail_bound_m_s_per_n[axis] += tail[axis + 2];
        }
    }
    Ok((
        force_cartesian(
            sums.map(ComplexSum::value),
            geometry.spectral.tangent,
            layer.receiving_normal(),
        )?,
        diagnostics,
    ))
}

struct FaceForceKernel<'a> {
    common: FaceKernel<'a>,
    log_pressure: f64,
    log_velocity: f64,
}

impl Spectrum<7> for FaceForceKernel<'_> {
    fn regions(&self) -> usize {
        <FaceKernel<'_> as Spectrum<3>>::regions(&self.common)
    }
    fn extent(&self) -> (f64, f64) {
        <FaceKernel<'_> as Spectrum<3>>::extent(&self.common)
    }
    fn resolved(&self, region: usize, lower: f64, upper: f64) -> bool {
        <FaceKernel<'_> as Spectrum<3>>::resolved(&self.common, region, lower, upper)
    }
    fn sample(&self, region: usize, x: f64) -> Result<[C64; 7], AeroacError> {
        let FaceWeight {
            spectrum,
            pressure,
            velocity,
            direction,
            log_coefficient,
        } = self.common.sample_weight(region, x)?;
        let argument = spectrum.q * self.common.common.coordinates.range;
        let zero = j0(argument);
        let one = j1(argument);
        let ratio = j1_over_argument(argument, one);
        let q = det::ln(spectrum.q);
        let z = unit(spectrum.zeta[1]);
        let log_z = det::ln(spectrum.zeta[1].abs());
        let sigma = if self.common.face == 0 { 1.0 } else { -1.0 };
        let i = C64::new(0.0, 1.0);
        let pressure_direction = direction * pressure.unit;
        let pressure_log = log_coefficient + pressure.log_magnitude;
        let velocity_direction = direction * velocity.unit;
        let velocity_log = log_coefficient + velocity.log_magnitude + self.common.rho_log[1];
        let values = [
            spectrum.scaled(
                i * pressure_direction,
                pressure_log + self.log_pressure + q,
                one,
            ),
            spectrum.scaled(
                z * pressure_direction,
                pressure_log + self.log_pressure + log_z,
                -sigma * zero,
            ),
            spectrum.scaled(
                pressure_direction,
                pressure_log + self.log_velocity + 2.0 * q,
                zero - ratio,
            ),
            spectrum.scaled(
                pressure_direction,
                pressure_log + self.log_velocity + 2.0 * q,
                ratio,
            ),
            spectrum.scaled(
                -i * z * pressure_direction,
                pressure_log + self.log_velocity + q + log_z,
                sigma * one,
            ),
            spectrum.scaled(
                i * velocity_direction,
                velocity_log + self.log_velocity + q,
                sigma * one,
            ),
            spectrum.scaled(
                z * velocity_direction,
                velocity_log + self.log_velocity + log_z,
                -zero,
            ),
        ];
        if values.iter().any(|v| !finite(*v)) {
            return Err(nonfinite("internal force spectrum"));
        }
        Ok(values)
    }
}
