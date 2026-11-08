//! Physical force derivatives of the SAME finite-impedance reflected Green field.
//!
//! A force ON fluid [N] differentiates the original source position, not the
//! image position: p_F = i/(omega*rho) F.grad_source(p_Q). The shared seven
//! reflected-force channels retain the opposite mixed velocity signs, normal
//! image-force reversal, and finite transverse velocity on the symmetry axis.

use super::super::force::{admit, equal_speed};
use super::super::reflection::{
    ReflectionCoordinates, apply_force, field, force_cartesian, line,
    sample_force_angular_derivative, sample_force_spectrum, singular_total,
};
use super::super::{FullWaveForceLimits, Spectrum, finite, nonfinite, poll};
use super::{PlanarImpedanceWall, WallKernel, integrate_wall_spectrum};
use crate::AeroacError;
use crate::curle3d::statistics::ComplexSum;
use crate::curle3d::{CoherentCompactField, CompactSourcePhasor, PreparedCompactSourceObserver};
use crate::observation::HarmonicFieldLine;
use core::f64::consts::PI;
use fs_math::{c64::C64, det};

const REFLECTED_SCOPE: &str = "Prescribed stationary compact 3D force ON fluid [N]; \
reflected spherical-wave pressure and velocity in one homogeneous lossless fluid \
at an infinite passive locally reacting impedance wall. Full propagating, evanescent \
and transverse near fields; numerical quadrature estimates and analytic tail bounds, \
not interval certificates. Includes outgoing lossless compliant surface waves. No finite edges, mean flow, bulk \
absorption, source feedback or measured-SPL qualification.";

const TOTAL_SCOPE: &str = "Prescribed stationary compact 3D force ON fluid [N]; \
coherent direct-plus-reflected pressure and velocity in one homogeneous lossless fluid \
at an infinite passive locally reacting impedance wall. Reflected-transfer absolute \
errors carry to the total; relative accuracy is not guaranteed at interference nulls. \
Includes outgoing lossless compliant surface waves. No singular direct self-field, finite edges, mean flow, bulk \
absorption, source feedback, interval certificate or measured-SPL qualification.";

type Transfer = [[C64; 3]; 4];

/// Errors for unit force transfers before applying the original complex force.
/// Pressure entries are (t,n), velocity entries (tt,ss,tn,nt,nn). Their Euclidean
/// norms bound the pressure-vector and velocity-tensor estimates respectively;
/// multiply by ||F||_2 for a particular force. Geometry, Bessel, coordinate and
/// floating-point errors are outside these estimates and analytic tail bounds.
#[derive(Clone, Copy, Debug)]
pub struct ImpedanceWallForceDiagnostics {
    pub evaluations: usize,
    pub maximum_depth_reached: usize,
    pub retained_panels: usize,
    pub transverse_wavenumber_cutoff_rad_m: f64,
    pub maximum_bessel_argument: f64,
    pub pressure_quadrature_error_pa_per_n: [f64; 2],
    pub velocity_quadrature_error_m_s_per_n: [f64; 5],
    pub pressure_tail_bound_pa_per_n: [f64; 2],
    pub velocity_tail_bound_m_s_per_n: [f64; 5],
    /// Exact rigid/release image forces use no quadrature or spectral cutoff.
    pub exact_ideal_boundary: bool,
}

/// Two cached Cartesian 4-by-3 field matrices for one source and receiver.
/// Application preserves all complex force phases with compensated reductions,
/// without allocation, quadrature, numerical differentiation or trigonometry.
#[derive(Clone, Debug)]
pub struct PreparedImpedanceWallForce {
    wall: PlanarImpedanceWall,
    frequency_hz: f64,
    source_m: [f64; 3],
    observer_m: [f64; 3],
    reflected: Transfer,
    total: Option<Transfer>,
    diagnostics: ImpedanceWallForceDiagnostics,
}

impl PreparedImpedanceWallForce {
    /// The source lies strictly in the fluid; receivers may lie on the wall.
    /// A coincident receiver can query the finite reflected field only. Limits
    /// govern the reflected transfer per newton, not relative error at a force
    /// or interference null. The wall impedance belongs to this frequency.
    pub fn new(
        wall: PlanarImpedanceWall,
        frequency_hz: f64,
        source_m: [f64; 3],
        observer_m: [f64; 3],
        limits: FullWaveForceLimits,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, AeroacError> {
        poll(cancelled)?;
        admit(limits)?;
        let coordinates = ReflectionCoordinates::with_boundary_receiver(
            wall.geometry,
            frequency_hz,
            source_m,
            observer_m,
            true,
        )?;
        let mut diagnostics = ImpedanceWallForceDiagnostics {
            evaluations: 0,
            maximum_depth_reached: 0,
            retained_panels: 0,
            transverse_wavenumber_cutoff_rad_m: 0.0,
            maximum_bessel_argument: 0.0,
            pressure_quadrature_error_pa_per_n: [0.0; 2],
            velocity_quadrature_error_m_s_per_n: [0.0; 5],
            pressure_tail_bound_pa_per_n: [0.0; 2],
            velocity_tail_bound_m_s_per_n: [0.0; 5],
            exact_ideal_boundary: false,
        };
        let cylindrical = if wall.normalized.is_none() || wall.normalized == Some(C64::ZERO) {
            let mut values = equal_speed(
                coordinates.omega,
                wall.fluid().sound_speed_m_s(),
                coordinates.range,
                coordinates.height.to_f64(),
                wall.fluid().density_kg_m3(),
                [1.0, 1.0],
            )?;
            let sign = if wall.normalized.is_none() { 1.0 } else { -1.0 };
            // Source-normal reversal changes nt and nn relative to the
            // positive-clearance free-space tensor. Both mixed signs matter.
            for (index, value) in values.iter_mut().enumerate() {
                *value = value.scale(if index >= 5 { -sign } else { sign });
            }
            diagnostics.exact_ideal_boundary = true;
            values
        } else {
            let log_pressure = -det::ln(4.0 * PI);
            let log_velocity =
                log_pressure - det::ln(coordinates.omega) - det::ln(wall.fluid().density_kg_m3());
            let mut common = WallKernel::new(wall, coordinates)?;
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
            let kernel = WallForceKernel {
                common,
                log_pressure,
                log_velocity,
            };
            let result = integrate_wall_spectrum::<7, 2>(
                &kernel,
                &kernel.common,
                tail,
                limits.into(),
                [[1, 0], [0, 1], [2, 0], [2, 0], [1, 1], [1, 1], [0, 2]],
                |region, x| {
                    sample_force_angular_derivative(
                        kernel.common.sample_weight(region, x)?,
                        coordinates.range,
                        log_pressure,
                        log_velocity,
                    )
                },
                cancelled,
            )?;
            diagnostics.evaluations = result.evaluations;
            diagnostics.maximum_depth_reached = result.maximum_depth_reached;
            diagnostics.retained_panels = result.retained_panels;
            diagnostics.transverse_wavenumber_cutoff_rad_m = result.cutoff;
            diagnostics.maximum_bessel_argument = result.maximum_argument;
            diagnostics.pressure_quadrature_error_pa_per_n = [result.error[0], result.error[1]];
            diagnostics.velocity_quadrature_error_m_s_per_n =
                core::array::from_fn(|i| result.error[i + 2]);
            diagnostics.pressure_tail_bound_pa_per_n = [tail[0], tail[1]];
            diagnostics.velocity_tail_bound_m_s_per_n = core::array::from_fn(|i| tail[i + 2]);
            result.value
        };
        let reflected = force_cartesian(cylindrical, coordinates.tangent, wall.outward_normal())?;
        let total = if source_m == observer_m {
            None
        } else {
            let direct = PreparedCompactSourceObserver::new(
                wall.fluid(),
                frequency_hz,
                observer_m,
                &[source_m],
                1,
            )?;
            let mut total = [[C64::ZERO; 3]; 4];
            for column in 0..3 {
                let unit = core::array::from_fn(|i| if i == column { C64::ONE } else { C64::ZERO });
                let value =
                    direct.evaluate_harmonic_line(&[CompactSourcePhasor::force_dipole(unit)])?;
                let values = [
                    value.pressure_pa,
                    value.particle_velocity_m_s[0],
                    value.particle_velocity_m_s[1],
                    value.particle_velocity_m_s[2],
                ];
                for row in 0..4 {
                    let mut sum = ComplexSum::default();
                    sum.add(values[row])?;
                    sum.add(reflected[row][column])?;
                    total[row][column] = sum.value();
                }
            }
            if total.iter().flatten().any(|z| !finite(*z)) {
                return Err(nonfinite("impedance wall total force transfer"));
            }
            Some(total)
        };
        poll(cancelled)?;
        Ok(Self {
            wall,
            frequency_hz,
            source_m,
            observer_m,
            reflected,
            total,
            diagnostics,
        })
    }

    pub const fn wall(&self) -> PlanarImpedanceWall {
        self.wall
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
    pub const fn diagnostics(&self) -> &ImpedanceWallForceDiagnostics {
        &self.diagnostics
    }
    pub const fn total_field_available(&self) -> bool {
        self.total.is_some()
    }

    pub fn evaluate_reflected(
        &self,
        force_on_fluid_n: [C64; 3],
    ) -> Result<CoherentCompactField, AeroacError> {
        field(
            apply_force(&self.reflected, force_on_fluid_n)?,
            REFLECTED_SCOPE,
        )
    }

    /// Sum the fields before forming active/reactive intensity, including
    /// direct/reflected and cross-axis force interference.
    pub fn evaluate_total(
        &self,
        force_on_fluid_n: [C64; 3],
    ) -> Result<CoherentCompactField, AeroacError> {
        field(
            apply_force(
                self.total.as_ref().ok_or_else(singular_total)?,
                force_on_fluid_n,
            )?,
            TOTAL_SCOPE,
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

struct WallForceKernel {
    common: WallKernel,
    log_pressure: f64,
    log_velocity: f64,
}

impl Spectrum<7> for WallForceKernel {
    fn regions(&self) -> usize {
        <WallKernel as Spectrum<3>>::regions(&self.common)
    }
    fn extent(&self) -> (f64, f64) {
        <WallKernel as Spectrum<3>>::extent(&self.common)
    }
    fn resolved(&self, region: usize, lower: f64, upper: f64) -> bool {
        <WallKernel as Spectrum<3>>::resolved(&self.common, region, lower, upper)
    }
    fn sample(&self, region: usize, x: f64) -> Result<[C64; 7], AeroacError> {
        sample_force_spectrum(
            self.common.sample_weight(region, x)?,
            self.common.coordinates.range,
            self.log_pressure,
            self.log_velocity,
        )
    }
}
