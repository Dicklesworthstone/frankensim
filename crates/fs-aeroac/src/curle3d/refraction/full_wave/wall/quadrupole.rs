//! Physical volume-stress Hessian of the finite-impedance reflected Green field.
//!
//! A stationary symmetric Lighthill moment [N m] differentiates the original
//! source coordinates twice. The existing reflected quadrupole channels retain
//! trace, both physical shear contractions and all receiver velocity derivatives;
//! the original wall kernel supplies the passive spectral reflection, including
//! its evanescent amplification and lossy or outgoing lossless surface mode.

use super::super::quadrupole::{FullWaveQuadrupoleLimits, cartesian};
use super::super::reflection::{
    ReflectionCoordinates, apply_quadrupole, field, line, quadrupole_angular_derivative_channels,
    quadrupole_channels, quadrupole_image_transfer, singular_total,
};
use super::super::{QuadratureLimits, Spectrum, finite, nonfinite, poll};
use super::{PlanarImpedanceWall, WallKernel, integrate_wall_spectrum};
use crate::AeroacError;
use crate::curle3d::CoherentCompactField;
use crate::curle3d::quadrupole::{PreparedQuadrupoleObserver, QuadrupolePhasor};
use crate::curle3d::statistics::ComplexSum;
use crate::observation::HarmonicFieldLine;
use core::f64::consts::{PI, SQRT_2};
use fs_math::{c64::C64, det};

mod array;
pub use array::{
    IMPEDANCE_WALL_REFLECTED_QUADRUPOLE_ARRAY_SCOPE, IMPEDANCE_WALL_TOTAL_QUADRUPOLE_ARRAY_SCOPE,
    PreparedImpedanceWallQuadrupoleObserver,
};

const REFLECTED_SCOPE: &str = "Prescribed stationary symmetric 3D Lighthill stress \
moment [N m]; reflected spherical-wave pressure and vector velocity in one homogeneous \
lossless fluid at an infinite passive locally reacting impedance wall. Full trace, \
physical shear, propagating, evanescent and transverse near fields; numerical quadrature \
estimates and analytic tail bounds, not interval certificates. Includes outgoing lossless \
compliant surface waves. No finite edges, mean flow, bulk absorption, unresolved SGS sound, source feedback \
or measured-SPL qualification.";

const TOTAL_SCOPE: &str = "Prescribed stationary symmetric 3D Lighthill stress \
moment [N m]; coherent direct-plus-reflected pressure and vector velocity in one \
homogeneous lossless fluid at an infinite passive locally reacting impedance wall. \
Reflected-transfer absolute errors carry to the total; relative accuracy is not \
guaranteed at interference nulls. Includes outgoing lossless compliant surface waves. \
No singular direct self-field, finite edges, mean flow, bulk absorption, unresolved SGS sound, source feedback, \
interval certificate or measured-SPL qualification.";

type Transfer = [[C64; 6]; 4];

/// Errors for unit symmetric stress transfers before applying the original
/// complex moment. Pressure channels are (tt,ss,sqrt(2) tn,nn); velocity channels
/// are (ttt,tss,sqrt(2) ttn,tnn,sqrt(2) sts,sqrt(2) ssn,ntt,nss,sqrt(2) ntn,nnn),
/// with receiver derivative first. Their Euclidean norms bound the respective
/// Frobenius transfer estimates; multiply by ||M||_F, counting physical shear
/// twice, for a particular source. Geometry, Bessel and floating-point errors
/// are outside these estimates and analytic tail bounds.
#[derive(Clone, Copy, Debug)]
pub struct ImpedanceWallQuadrupoleDiagnostics {
    pub evaluations: usize,
    pub maximum_depth_reached: usize,
    pub retained_panels: usize,
    pub transverse_wavenumber_cutoff_rad_m: f64,
    pub maximum_bessel_argument: f64,
    pub pressure_quadrature_error_pa_per_n_m: [f64; 4],
    pub velocity_quadrature_error_m_s_per_n_m: [f64; 10],
    pub pressure_tail_bound_pa_per_n_m: [f64; 4],
    pub velocity_tail_bound_m_s_per_n_m: [f64; 10],
    /// Exact rigid/release image stresses need no quadrature or cutoff.
    pub exact_ideal_boundary: bool,
}

/// Cached Cartesian 4-by-6 reflected and coherent total field matrices.
/// The physical symmetric stress order is (xx,yy,zz,xy,xz,yz), with no Voigt
/// engineering-shear multiplier. Application preserves all source phases with
/// compensated reductions, without allocation, quadrature or differentiation.
#[derive(Clone, Debug)]
pub struct PreparedImpedanceWallQuadrupole {
    wall: PlanarImpedanceWall,
    frequency_hz: f64,
    source_m: [f64; 3],
    observer_m: [f64; 3],
    reflected: Transfer,
    total: Option<Transfer>,
    diagnostics: ImpedanceWallQuadrupoleDiagnostics,
}

impl PreparedImpedanceWallQuadrupole {
    /// The source lies strictly in the fluid; the receiver may lie on the wall.
    /// A coincident receiver retains only the finite reflected response. Limits
    /// control the complete reflected stress transfer per newton metre, before
    /// applying a source. The wall impedance belongs to this frequency.
    pub fn new(
        wall: PlanarImpedanceWall,
        frequency_hz: f64,
        source_m: [f64; 3],
        observer_m: [f64; 3],
        limits: FullWaveQuadrupoleLimits,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, AeroacError> {
        poll(cancelled)?;
        QuadratureLimits::from(limits).admit()?;
        let coordinates = ReflectionCoordinates::with_boundary_receiver(
            wall.geometry,
            frequency_hz,
            source_m,
            observer_m,
            true,
        )?;
        let mut diagnostics = ImpedanceWallQuadrupoleDiagnostics {
            evaluations: 0,
            maximum_depth_reached: 0,
            retained_panels: 0,
            transverse_wavenumber_cutoff_rad_m: 0.0,
            maximum_bessel_argument: 0.0,
            pressure_quadrature_error_pa_per_n_m: [0.0; 4],
            velocity_quadrature_error_m_s_per_n_m: [0.0; 10],
            pressure_tail_bound_pa_per_n_m: [0.0; 4],
            velocity_tail_bound_m_s_per_n_m: [0.0; 10],
            exact_ideal_boundary: false,
        };
        let reflected = if wall.normalized.is_none() || wall.normalized == Some(C64::ZERO) {
            diagnostics.exact_ideal_boundary = true;
            quadrupole_image_transfer(
                wall.geometry,
                frequency_hz,
                source_m,
                observer_m,
                if wall.normalized.is_none() { 1.0 } else { -1.0 },
                cancelled,
            )?
        } else {
            let log_pressure = -det::ln(4.0 * PI);
            let log_velocity =
                log_pressure - det::ln(coordinates.omega) - det::ln(wall.fluid().density_kg_m3());
            let mut common = WallKernel::new(wall, coordinates)?;
            // Derivative orders and Frobenius weights are the same as for
            // two-fluid reflection. The wall owner additionally bounds the
            // full passive evanescent amplification instead of clipping |R|.
            let [p, v] = common.prepare_tail(
                limits.into(),
                [2, 3],
                [log_pressure, log_velocity],
                [5, 14],
                cancelled,
            )?;
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
            let kernel = WallQuadrupoleKernel {
                common,
                log_pressure,
                log_velocity,
            };
            let result = integrate_wall_spectrum::<14, 4>(
                &kernel,
                &kernel.common,
                tail,
                limits.into(),
                [
                    [2, 0],
                    [2, 0],
                    [1, 1],
                    [0, 2],
                    [3, 0],
                    [3, 0],
                    [2, 1],
                    [1, 2],
                    [3, 0],
                    [2, 1],
                    [2, 1],
                    [2, 1],
                    [1, 2],
                    [0, 3],
                ],
                |region, x| {
                    quadrupole_angular_derivative_channels(
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
            diagnostics.pressure_quadrature_error_pa_per_n_m =
                core::array::from_fn(|i| result.error[i]);
            diagnostics.velocity_quadrature_error_m_s_per_n_m =
                core::array::from_fn(|i| result.error[i + 4]);
            diagnostics.pressure_tail_bound_pa_per_n_m = core::array::from_fn(|i| tail[i]);
            diagnostics.velocity_tail_bound_m_s_per_n_m = core::array::from_fn(|i| tail[i + 4]);
            cartesian(result.value, coordinates.tangent, wall.outward_normal())?
        };
        if reflected.iter().flatten().any(|z| !finite(*z)) {
            return Err(nonfinite("impedance wall reflected quadrupole transfer"));
        }
        let total = if source_m == observer_m {
            None
        } else {
            let direct = PreparedQuadrupoleObserver::new_cancellable(
                wall.fluid(),
                frequency_hz,
                observer_m,
                &[source_m],
                1,
                &mut *cancelled,
            )?;
            let mut total = reflected;
            for column in 0..6 {
                poll(cancelled)?;
                let mut source = QuadrupolePhasor::default();
                source.stress_moment_n_m[column] = C64::ONE;
                let value = direct.evaluate_harmonic_line(&[source])?;
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
                return Err(nonfinite("impedance wall total quadrupole transfer"));
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
    pub const fn diagnostics(&self) -> &ImpedanceWallQuadrupoleDiagnostics {
        &self.diagnostics
    }
    pub const fn has_total_field(&self) -> bool {
        self.total.is_some()
    }

    /// Apply the full physical symmetric moment on its exp(-i omega t) clock.
    pub fn evaluate_reflected(
        &self,
        source: QuadrupolePhasor,
    ) -> Result<CoherentCompactField, AeroacError> {
        field(apply_quadrupole(&self.reflected, source)?, REFLECTED_SCOPE)
    }

    /// Sum complex pressure and velocity before forming the energy flux.
    pub fn evaluate_total(
        &self,
        source: QuadrupolePhasor,
    ) -> Result<CoherentCompactField, AeroacError> {
        field(
            apply_quadrupole(self.total.as_ref().ok_or_else(singular_total)?, source)?,
            TOTAL_SCOPE,
        )
    }

    pub fn evaluate_reflected_harmonic_line(
        &self,
        source: QuadrupolePhasor,
    ) -> Result<HarmonicFieldLine<3>, AeroacError> {
        Ok(line(
            self.frequency_hz,
            apply_quadrupole(&self.reflected, source)?,
        ))
    }

    pub fn evaluate_total_harmonic_line(
        &self,
        source: QuadrupolePhasor,
    ) -> Result<HarmonicFieldLine<3>, AeroacError> {
        Ok(line(
            self.frequency_hz,
            apply_quadrupole(self.total.as_ref().ok_or_else(singular_total)?, source)?,
        ))
    }
}

struct WallQuadrupoleKernel {
    common: WallKernel,
    log_pressure: f64,
    log_velocity: f64,
}

impl Spectrum<14> for WallQuadrupoleKernel {
    fn regions(&self) -> usize {
        <WallKernel as Spectrum<3>>::regions(&self.common)
    }
    fn extent(&self) -> (f64, f64) {
        <WallKernel as Spectrum<3>>::extent(&self.common)
    }
    fn resolved(&self, region: usize, lower: f64, upper: f64) -> bool {
        <WallKernel as Spectrum<3>>::resolved(&self.common, region, lower, upper)
    }
    fn sample(&self, region: usize, x: f64) -> Result<[C64; 14], AeroacError> {
        quadrupole_channels(
            self.common.sample_weight(region, x)?,
            self.common.coordinates.range,
            self.log_pressure,
            self.log_velocity,
        )
    }
}
