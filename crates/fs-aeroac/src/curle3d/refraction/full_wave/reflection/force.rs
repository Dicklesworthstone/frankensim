//! Original physical force derivatives of the reflected spherical Green function.
//!
//! For W=q*dq*R*exp(i*zeta1*(a+b))/zeta1, the two pressure coefficients
//! are integral W*[i*q*J1,zeta1*J0]/(4*pi). The five velocity coefficients
//! (tt,ss,tn,nt,nn) are integral W*[q²*(J0-J1/x),q²*J1/x,
//! i*q*zeta1*J1,-i*q*zeta1*J1,-zeta1²*J0]/(4*pi*omega*rho1).
//! The opposite mixed signs and negative normal-normal term reflect the
//! original SOURCE normal derivative. An image force reverses its normal
//! component; treating it as an unchanged force gives the wrong reflected field.

use super::super::force::{admit, equal_speed, j1_over_argument};
use super::super::quadrupole::j2_over_argument;
use super::super::{
    FullWaveForceDiagnostics, FullWaveForceLimits, PlanarFluidInterface, Spectrum, finite,
    integrate_spectrum, j0, j1, nonfinite, poll, positive,
};
use super::{
    ReflectionCoordinates, ReflectionKernel, ReflectionSpectralWeight, exact_reflection, field,
    line, singular_total,
};
use crate::AeroacError;
use crate::curle3d::statistics::ComplexSum;
use crate::curle3d::{CoherentCompactField, CompactSourcePhasor, PreparedCompactSourceObserver};
use crate::observation::HarmonicFieldLine;
use core::f64::consts::PI;
use fs_math::{c64::C64, det};

pub const FULL_WAVE_REFLECTED_FORCE_SCOPE: &str = "Prescribed stationary compact 3D \
force ON fluid [N]; reflected-only spherical-wave field on the incident side of \
one infinite plane between original stationary homogeneous lossless fluids, \
including reactive, transverse and evanescent motion. Estimated quadrature error \
and analytic spectral-tail bounds are not interval certificates; no direct self-field, \
source feedback, absorption, finite edges, mean flow or measured-SPL qualification.";

pub const FULL_WAVE_TOTAL_FORCE_SCOPE: &str = "Prescribed stationary compact 3D force \
ON fluid [N]; coherent direct-plus-reflected spherical-wave field on the incident \
side of one infinite original lossless fluid interface. Relative quadrature control \
applies to the reflected transfer; only its absolute error estimate carries to total \
fields near interference nulls. No interval certificate, singular direct self-field, \
source feedback, absorption, finite edges, mean flow or measured-SPL qualification.";

type Transfer = [[C64; 3]; 4];

/// Original incident-side force response with separate reflected and total fields.
/// Preparation caches two Cartesian 4-by-3 matrices. Every later evaluation
/// preserves complex force phases and requires no allocation or integration.
///
/// Coincident source and observer positions retain a finite reflected field;
/// total-field requests then refuse the singular direct point-source field.
#[derive(Clone, Debug)]
pub struct PreparedReflectedFullWaveForce {
    interface: PlanarFluidInterface,
    frequency_hz: f64,
    source_m: [f64; 3],
    observer_m: [f64; 3],
    reflected: Transfer,
    total: Option<Transfer>,
    diagnostics: FullWaveForceDiagnostics,
}

impl PreparedReflectedFullWaveForce {
    /// Prepare a physical force and an observer strictly in the incident fluid.
    /// Limits control reflected pressure-vector and velocity-tensor norms per N.
    /// The total field inherits those absolute estimates, not a relative error
    /// guarantee at a coherent null. Original homogeneous direct-field arithmetic
    /// and coordinate/Bessel roundoff remain outside quadrature diagnostics.
    pub fn new(
        interface: PlanarFluidInterface,
        frequency_hz: f64,
        source_m: [f64; 3],
        observer_m: [f64; 3],
        limits: FullWaveForceLimits,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, AeroacError> {
        poll(cancelled)?;
        admit(limits)?;
        let coordinates =
            ReflectionCoordinates::new(interface, frequency_hz, source_m, observer_m)?;
        let incident = interface.incident_state();
        let receiving = interface.receiving_state();
        if !positive(coordinates.omega / incident.sound_speed_m_s()) {
            return Err(nonfinite("reflected force wavenumber"));
        }
        let (cylindrical, diagnostics) =
            if incident.sound_speed_m_s() == receiving.sound_speed_m_s() {
                let reflection = exact_reflection(interface);
                let mut field = [C64::ZERO; 7];
                if reflection != 0.0 {
                    field = equal_speed(
                        coordinates.omega,
                        incident.sound_speed_m_s(),
                        coordinates.range,
                        coordinates.height.to_f64(),
                        incident.density_kg_m3(),
                        [1.0, 1.0],
                    )?;
                    for (index, value) in field.iter_mut().enumerate() {
                        *value = value.scale(if index >= 5 { -reflection } else { reflection });
                    }
                }
                (
                    field,
                    FullWaveForceDiagnostics {
                        evaluations: 0,
                        maximum_depth_reached: 0,
                        retained_panels: 0,
                        transverse_wavenumber_cutoff_rad_m: 0.0,
                        maximum_bessel_argument: 0.0,
                        pressure_quadrature_error_pa_per_n: [0.0; 2],
                        velocity_quadrature_error_m_s_per_n: [0.0; 5],
                        pressure_tail_bound_pa_per_n: [0.0; 2],
                        velocity_tail_bound_m_s_per_n: [0.0; 5],
                        exact_equal_speed: true,
                    },
                )
            } else {
                let log_pressure = -det::ln(4.0 * PI);
                let log_velocity =
                    log_pressure - det::ln(coordinates.omega) - det::ln(incident.density_kg_m3());
                let mut common = ReflectionKernel::new(interface, coordinates)?;
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
                let kernel = ReflectionForceKernel {
                    common,
                    log_pressure,
                    log_velocity,
                };
                let result = integrate_spectrum::<7, 2>(&kernel, tail, limits.into(), cancelled)?;
                (
                    result.value,
                    FullWaveForceDiagnostics {
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
                        exact_equal_speed: false,
                    },
                )
            };
        let reflected = cartesian(
            cylindrical,
            coordinates.tangent,
            interface.receiving_normal(),
        )?;
        let total = if source_m == observer_m {
            None
        } else {
            let direct = PreparedCompactSourceObserver::new(
                incident,
                frequency_hz,
                observer_m,
                &[source_m],
                1,
            )?;
            let mut field = reflected;
            for column in 0..3 {
                poll(cancelled)?;
                let force =
                    core::array::from_fn(|axis| if axis == column { C64::ONE } else { C64::ZERO });
                let line =
                    direct.evaluate_harmonic_line(&[CompactSourcePhasor::force_dipole(force)])?;
                let values = [
                    line.pressure_pa,
                    line.particle_velocity_m_s[0],
                    line.particle_velocity_m_s[1],
                    line.particle_velocity_m_s[2],
                ];
                for row in 0..4 {
                    let mut sum = ComplexSum::default();
                    sum.add(values[row])?;
                    sum.add(reflected[row][column])?;
                    field[row][column] = sum.value();
                }
            }
            if field.iter().flatten().any(|v| !finite(*v)) {
                return Err(nonfinite("direct-plus-reflected force transfer"));
            }
            Some(field)
        };
        poll(cancelled)?;
        Ok(Self {
            interface,
            frequency_hz,
            source_m,
            observer_m,
            reflected,
            total,
            diagnostics,
        })
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
    /// Numerical errors apply to reflected transfer coefficients per newton.
    pub const fn diagnostics(&self) -> &FullWaveForceDiagnostics {
        &self.diagnostics
    }
    pub const fn has_total_field(&self) -> bool {
        self.total.is_some()
    }

    pub fn evaluate_reflected(
        &self,
        force_on_fluid_n: [C64; 3],
    ) -> Result<CoherentCompactField, AeroacError> {
        field(
            apply(&self.reflected, force_on_fluid_n)?,
            FULL_WAVE_REFLECTED_FORCE_SCOPE,
        )
    }

    /// Combine direct and reflected amplitudes before forming intensity.
    pub fn evaluate_total(
        &self,
        force_on_fluid_n: [C64; 3],
    ) -> Result<CoherentCompactField, AeroacError> {
        field(
            apply(self.total_transfer()?, force_on_fluid_n)?,
            FULL_WAVE_TOTAL_FORCE_SCOPE,
        )
    }

    pub fn evaluate_reflected_harmonic_line(
        &self,
        force_on_fluid_n: [C64; 3],
    ) -> Result<HarmonicFieldLine<3>, AeroacError> {
        Ok(line(
            self.frequency_hz,
            apply(&self.reflected, force_on_fluid_n)?,
        ))
    }

    pub fn evaluate_total_harmonic_line(
        &self,
        force_on_fluid_n: [C64; 3],
    ) -> Result<HarmonicFieldLine<3>, AeroacError> {
        Ok(line(
            self.frequency_hz,
            apply(self.total_transfer()?, force_on_fluid_n)?,
        ))
    }

    fn total_transfer(&self) -> Result<&Transfer, AeroacError> {
        self.total.as_ref().ok_or_else(singular_total)
    }
}

pub(in super::super) fn cartesian(
    value: [C64; 7],
    tangent: [f64; 3],
    normal: [f64; 3],
) -> Result<Transfer, AeroacError> {
    let [pt, pn, vtt, vss, vtn, vnt, vnn] = value;
    let mut field = [[C64::ZERO; 3]; 4];
    for column in 0..3 {
        field[0][column] = pt.scale(tangent[column]) + pn.scale(normal[column]);
        for row in 0..3 {
            let nn = normal[row] * normal[column];
            let identity = if row == column { 1.0 } else { 0.0 };
            let mut sum = ComplexSum::default();
            sum.add(vss.scale(identity - nn))?;
            sum.add((vtt - vss).scale(tangent[row] * tangent[column]))?;
            sum.add(vtn.scale(tangent[row] * normal[column]))?;
            sum.add(vnt.scale(normal[row] * tangent[column]))?;
            sum.add(vnn.scale(nn))?;
            field[row + 1][column] = sum.value();
        }
    }
    if field.iter().flatten().any(|v| !finite(*v)) {
        return Err(nonfinite("reflected force Cartesian transfer"));
    }
    Ok(field)
}

pub(in super::super) fn apply(
    transfer: &Transfer,
    force: [C64; 3],
) -> Result<[C64; 4], AeroacError> {
    if force.iter().any(|v| !finite(*v)) {
        return Err(nonfinite("reflected force phasor"));
    }
    let mut sums = [ComplexSum::default(); 4];
    for (column, amplitude) in force.into_iter().enumerate() {
        if amplitude != C64::ZERO {
            for (row, sum) in transfer.iter().zip(sums.iter_mut()) {
                sum.add(row[column] * amplitude)?;
            }
        }
    }
    let values = sums.map(ComplexSum::value);
    if values.iter().any(|v| !finite(*v)) {
        return Err(nonfinite("reflected force field"));
    }
    Ok(values)
}

struct ReflectionForceKernel {
    common: ReflectionKernel,
    log_pressure: f64,
    log_velocity: f64,
}

impl Spectrum<7> for ReflectionForceKernel {
    fn regions(&self) -> usize {
        <ReflectionKernel as Spectrum<3>>::regions(&self.common)
    }
    fn extent(&self) -> (f64, f64) {
        <ReflectionKernel as Spectrum<3>>::extent(&self.common)
    }
    fn resolved(&self, region: usize, lower: f64, upper: f64) -> bool {
        <ReflectionKernel as Spectrum<3>>::resolved(&self.common, region, lower, upper)
    }
    fn sample(&self, region: usize, x: f64) -> Result<[C64; 7], AeroacError> {
        let weight = self.common.sample_weight(region, x)?;
        sample_force_spectrum(
            weight,
            self.common.coordinates.range,
            self.log_pressure,
            self.log_velocity,
        )
    }
}

/// Analytic original-source and receiving derivatives of a reflected weight.
/// Interfaces and finite layers differ only in that Green-function weight.
pub(in super::super) fn sample_force_spectrum(
    weight: ReflectionSpectralWeight,
    range: f64,
    log_pressure: f64,
    log_velocity: f64,
) -> Result<[C64; 7], AeroacError> {
    force_channels(weight, range, log_pressure, log_velocity, false)
}

/// Derivative of each angular Bessel factor with respect to q*range, retaining
/// the original physical source/receiver powers and signs. The impedance-wall
/// pole subtraction uses this to evaluate its removable singularity exactly.
pub(in super::super) fn sample_force_angular_derivative(
    weight: ReflectionSpectralWeight,
    range: f64,
    log_pressure: f64,
    log_velocity: f64,
) -> Result<[C64; 7], AeroacError> {
    force_channels(weight, range, log_pressure, log_velocity, true)
}

fn force_channels(
    weight: ReflectionSpectralWeight,
    range: f64,
    log_pressure: f64,
    log_velocity: f64,
    angular_derivative: bool,
) -> Result<[C64; 7], AeroacError> {
    let argument = weight.q * range;
    if !argument.is_finite() {
        return Err(nonfinite("reflected force Bessel argument"));
    }
    let jzero = j0(argument);
    let jone = j1(argument);
    let ratio = j1_over_argument(argument, jone);
    let (jzero, jone, ratio) = if angular_derivative {
        (
            -jone,
            jzero - ratio,
            -j2_over_argument(argument, jzero, ratio),
        )
    } else {
        (jzero, jone, ratio)
    };
    let q_log = det::ln(weight.q);
    let magnitude = weight.zeta.abs();
    let direction = if magnitude == 0.0 {
        C64::ZERO
    } else {
        C64::new(weight.zeta.re / magnitude, weight.zeta.im / magnitude)
    };
    let z_log = det::ln(magnitude);
    let i = C64::new(0.0, 1.0);
    let values = [
        weight.scaled(i, log_pressure + q_log, jone),
        weight.scaled(direction, log_pressure + z_log, jzero),
        weight.scaled(C64::ONE, log_velocity + 2.0 * q_log, jzero - ratio),
        weight.scaled(C64::ONE, log_velocity + 2.0 * q_log, ratio),
        weight.scaled(i * direction, log_velocity + q_log + z_log, jone),
        weight.scaled(-i * direction, log_velocity + q_log + z_log, jone),
        weight.scaled(-(direction * direction), log_velocity + 2.0 * z_log, jzero),
    ];
    if values.iter().any(|v| !finite(*v)) {
        return Err(nonfinite("reflected force spectral field"));
    }
    Ok(values)
}
