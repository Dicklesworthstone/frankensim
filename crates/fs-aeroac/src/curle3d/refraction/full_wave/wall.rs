//! Spherical-wave reflection from an infinite locally reacting impedance wall.
//!
//! For exp(-i omega t), the outward normal points from fluid into wall and
//! `p = Z u_n`. The outgoing Weyl spectrum therefore has
//! `R(q) = (Z zeta - rho omega)/(Z zeta + rho omega)`, with
//! `zeta = sqrt(k²-q²)` on the outgoing branch. Applying a normal-incidence
//! coefficient to an image would omit both angular dependence and the
//! evanescent near field. This owner integrates that same spectrum and its
//! spatial derivatives through the parent's bounded Gauss--Kronrod owner.
//!
//! A passive compliant wall can amplify evanescent components without
//! generating power. For normalized z=Z/(rho c), their maximum amplitude is
//! `(abs(z)+Im(z))/Re(z)` when Im(z)>0, and at most one otherwise. That factor
//! enters the analytic tail bound; it is never clipped to one. The spectrum
//! is split at its evanescent maximum to resolve a narrow lossy surface mode.
//! Purely reactive compliance has a real-contour pole and explicitly refuses;
//! a lossless surface-wave contour prescription is outside this implementation.

use super::reflection::{
    self, ReflectionCoordinates, ReflectionSpectralWeight, prepare_reflected_tail,
};
use super::*;
use crate::curle3d::{CompactSourcePhasor, PreparedCompactSourceObserver};
use fs_material::fluid_acoustics::FluidAcousticState;

mod array;
mod force;
mod quadrupole;
pub use array::{
    IMPEDANCE_WALL_REFLECTED_ARRAY_SCOPE, IMPEDANCE_WALL_TOTAL_ARRAY_SCOPE,
    ImpedanceWallCompactSourceDiagnostics, ImpedanceWallFieldKind,
    PreparedImpedanceWallCompactSourceObserver,
};
pub use force::{ImpedanceWallForceDiagnostics, PreparedImpedanceWallForce};
pub use quadrupole::{ImpedanceWallQuadrupoleDiagnostics, PreparedImpedanceWallQuadrupole};

const SCOPE: &str = "Prescribed stationary compact 3D volume source in one homogeneous \
lossless fluid, spherical reflection from an infinite passive locally reacting impedance \
plane. Full propagating/evanescent pressure and velocity, coherent direct/reflected \
interference, estimated quadrature error and separate analytic tail bound. No lossless \
compliant surface-wave pole, nonlocal wall response, finite edges, mean flow, bulk \
absorption, source feedback or measured-SPL qualification.";

/// Original fluid, wall geometry and specific impedance [Pa s/m], at the
/// frequency supplied to preparation. Frequency-dependent material laws must
/// be evaluated by their original owner for each requested frequency.
#[derive(Clone, Copy, Debug)]
pub struct PlanarImpedanceWall {
    geometry: PlanarFluidInterface,
    impedance_pa_s_m: Option<C64>,
    normalized: Option<C64>,
}

impl PlanarImpedanceWall {
    /// `outward_normal` points from the fluid INTO the wall. `None` denotes
    /// exact rigidity and zero impedance denotes exact pressure release.
    /// Finite impedances require nonnegative resistance. Pure compliance
    /// (zero resistance and positive reactance in exp(-i omega t)) refuses
    /// because its surface mode lies on the integration contour.
    pub fn new(
        fluid: FluidAcousticState,
        point_m: [f64; 3],
        outward_normal: [f64; 3],
        impedance_pa_s_m: Option<C64>,
    ) -> Result<Self, AeroacError> {
        let geometry = PlanarFluidInterface::new(fluid, fluid, point_m, outward_normal)?;
        let normalized = match impedance_pa_s_m {
            None => None,
            Some(z) => {
                if !finite(z) || z.re < 0.0 || (z.re == 0.0 && z.im > 0.0) {
                    return Err(invalid(
                        "passive wall impedance excludes a lossless compliant pole",
                    ));
                }
                let reference = fluid.density_kg_m3() * fluid.sound_speed_m_s();
                if !positive(reference) {
                    return Err(nonfinite("wall acoustic impedance reference"));
                }
                let normalized = C64::new(z.re / reference, z.im / reference);
                if !finite(normalized) || (normalized == C64::ZERO && z != C64::ZERO) {
                    return Err(nonfinite("normalized wall impedance"));
                }
                Some(normalized)
            }
        };
        Ok(Self {
            geometry,
            impedance_pa_s_m,
            normalized,
        })
    }

    pub const fn fluid(self) -> FluidAcousticState {
        self.geometry.incident_state()
    }
    pub const fn point_m(self) -> [f64; 3] {
        self.geometry.point_m()
    }
    pub const fn outward_normal(self) -> [f64; 3] {
        self.geometry.receiving_normal()
    }
    pub const fn impedance_pa_s_m(self) -> Option<C64> {
        self.impedance_pa_s_m
    }
}

/// Numerical evidence for reflected unit-volume-flow pressure and velocity.
/// The same absolute errors apply to the total field; relative accuracy at a
/// direct/reflected cancellation is not guaranteed. These are estimates plus
/// a floating-point tail bound, not interval certificates or geometry errors.
#[derive(Clone, Copy, Debug)]
pub struct ImpedanceWallDiagnostics {
    pub evaluations: usize,
    pub maximum_depth_reached: usize,
    pub retained_panels: usize,
    pub transverse_wavenumber_cutoff_rad_m: f64,
    pub maximum_bessel_argument: f64,
    pub pressure_quadrature_error_pa_per_m3_s: f64,
    pub velocity_quadrature_error_per_m2: [f64; 2],
    pub pressure_tail_bound_pa_per_m3_s: f64,
    pub velocity_tail_bound_per_m2: [f64; 2],
    /// Rigid/release walls use their exact image, with no quadrature.
    pub exact_ideal_boundary: bool,
}

/// Cached spherical response; applying a source phasor allocates nothing and
/// performs no quadrature or transcendental evaluation.
#[derive(Clone, Debug)]
pub struct PreparedImpedanceWallMonopole {
    wall: PlanarImpedanceWall,
    frequency_hz: f64,
    reflected: [C64; 4],
    total: Option<[C64; 4]>,
    diagnostics: ImpedanceWallDiagnostics,
}

impl PreparedImpedanceWallMonopole {
    /// Sources lie strictly in the fluid; receivers may lie on the plane.
    /// Coincident source/receiver positions retain the finite reflected field
    /// but refuse the singular direct-plus-reflected field. Cancellation and
    /// work exhaustion return no partially prepared response.
    pub fn new(
        wall: PlanarImpedanceWall,
        frequency_hz: f64,
        source_m: [f64; 3],
        observer_m: [f64; 3],
        limits: FullWaveMonopoleLimits,
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
        let mut diagnostics = ImpedanceWallDiagnostics {
            evaluations: 0,
            maximum_depth_reached: 0,
            retained_panels: 0,
            transverse_wavenumber_cutoff_rad_m: 0.0,
            maximum_bessel_argument: 0.0,
            pressure_quadrature_error_pa_per_m3_s: 0.0,
            velocity_quadrature_error_per_m2: [0.0; 2],
            pressure_tail_bound_pa_per_m3_s: 0.0,
            velocity_tail_bound_per_m2: [0.0; 2],
            exact_ideal_boundary: false,
        };
        let reflected = if wall.normalized.is_none() || wall.normalized == Some(C64::ZERO) {
            let distance = wall.geometry.distance(source_m);
            let image = core::array::from_fn(|axis| {
                (Dd::from_f64(source_m[axis])
                    - distance * Dd::from_f64(2.0 * wall.outward_normal()[axis]))
                .to_f64()
            });
            let sign = if wall.normalized.is_none() { 1.0 } else { -1.0 };
            let observer = PreparedCompactSourceObserver::new(
                wall.fluid(),
                frequency_hz,
                observer_m,
                &[image],
                1,
            )?;
            let value = observer
                .evaluate_harmonic_line(&[CompactSourcePhasor::monopole(C64::from_re(sign))])?;
            diagnostics.exact_ideal_boundary = true;
            [
                value.pressure_pa,
                value.particle_velocity_m_s[0],
                value.particle_velocity_m_s[1],
                value.particle_velocity_m_s[2],
            ]
        } else {
            let mut kernel = WallKernel::new(wall, coordinates)?;
            let bound = kernel.prepare_tail(
                limits.into(),
                [0, 1],
                [kernel.log_pressure, -det::ln(4.0 * PI)],
                [1, 2],
                cancelled,
            )?;
            let tail = [
                bound[0],
                if coordinates.range == 0.0 {
                    0.0
                } else {
                    bound[1]
                },
                bound[1],
            ];
            let result = integrate_spectrum::<3, 1>(&kernel, tail, limits.into(), cancelled)?;
            diagnostics.evaluations = result.evaluations;
            diagnostics.maximum_depth_reached = result.maximum_depth_reached;
            diagnostics.retained_panels = result.retained_panels;
            diagnostics.transverse_wavenumber_cutoff_rad_m = result.cutoff;
            diagnostics.maximum_bessel_argument = result.maximum_argument;
            diagnostics.pressure_quadrature_error_pa_per_m3_s = result.error[0];
            diagnostics.velocity_quadrature_error_per_m2 = [result.error[1], result.error[2]];
            diagnostics.pressure_tail_bound_pa_per_m3_s = tail[0];
            diagnostics.velocity_tail_bound_per_m2 = [tail[1], tail[2]];
            let [pressure, radial, normal] = result.value;
            let n = wall.outward_normal();
            let t = coordinates.tangent;
            [
                pressure,
                radial.scale(t[0]) + normal.scale(n[0]),
                radial.scale(t[1]) + normal.scale(n[1]),
                radial.scale(t[2]) + normal.scale(n[2]),
            ]
        };
        if reflected.iter().any(|z| !finite(*z)) {
            return Err(nonfinite("wall reflected field"));
        }
        let total = if source_m == observer_m {
            None
        } else {
            let direct = PreparedCompactSourceObserver::new(
                wall.fluid(),
                frequency_hz,
                observer_m,
                &[source_m],
                1,
            )?
            .evaluate_harmonic_line(&[CompactSourcePhasor::monopole(C64::ONE)])?;
            let direct = [
                direct.pressure_pa,
                direct.particle_velocity_m_s[0],
                direct.particle_velocity_m_s[1],
                direct.particle_velocity_m_s[2],
            ];
            let mut total = [C64::ZERO; 4];
            for row in 0..4 {
                let mut sum = ComplexSum::default();
                sum.add(direct[row])?;
                sum.add(reflected[row])?;
                total[row] = sum.value();
            }
            if total.iter().any(|z| !finite(*z)) {
                return Err(nonfinite("wall total field"));
            }
            Some(total)
        };
        poll(cancelled)?;
        Ok(Self {
            wall,
            frequency_hz,
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
    pub const fn diagnostics(&self) -> &ImpedanceWallDiagnostics {
        &self.diagnostics
    }
    pub const fn total_field_available(&self) -> bool {
        self.total.is_some()
    }

    pub fn evaluate_reflected(
        &self,
        volume_flow_m3_s: C64,
    ) -> Result<CoherentCompactField, AeroacError> {
        reflection::field(reflection::apply(self.reflected, volume_flow_m3_s)?, SCOPE)
    }
    pub fn evaluate_total(
        &self,
        volume_flow_m3_s: C64,
    ) -> Result<CoherentCompactField, AeroacError> {
        reflection::field(
            reflection::apply(
                self.total.ok_or_else(reflection::singular_total)?,
                volume_flow_m3_s,
            )?,
            SCOPE,
        )
    }
    pub fn evaluate_total_harmonic_line(
        &self,
        volume_flow_m3_s: C64,
    ) -> Result<HarmonicFieldLine<3>, AeroacError> {
        Ok(reflection::line(
            self.frequency_hz,
            reflection::apply(
                self.total.ok_or_else(reflection::singular_total)?,
                volume_flow_m3_s,
            )?,
        ))
    }
    pub fn evaluate_reflected_harmonic_line(
        &self,
        volume_flow_m3_s: C64,
    ) -> Result<HarmonicFieldLine<3>, AeroacError> {
        Ok(reflection::line(
            self.frequency_hz,
            reflection::apply(self.reflected, volume_flow_m3_s)?,
        ))
    }
}

struct WallKernel {
    coordinates: ReflectionCoordinates,
    k: f64,
    z: C64,
    unit_z: C64,
    log_z_over_k: f64,
    log_pressure: f64,
    log_gain: f64,
    tail_extent: f64,
    split: Option<f64>,
}

impl WallKernel {
    fn new(
        wall: PlanarImpedanceWall,
        coordinates: ReflectionCoordinates,
    ) -> Result<Self, AeroacError> {
        let k = coordinates.omega / wall.fluid().sound_speed_m_s();
        let z = wall
            .normalized
            .ok_or_else(|| invalid("rigid wall uses its exact image"))?;
        if !positive(k) || !positive(z.abs()) {
            return Err(nonfinite("wall spectral scale"));
        }
        let log_gain = if z.im > 0.0 {
            det::ln(z.abs()) + det::ln(1.0 + z.im / z.abs()) - det::ln(z.re)
        } else {
            0.0
        };
        Ok(Self {
            coordinates,
            k,
            z,
            unit_z: C64::new(z.re / z.abs(), z.im / z.abs()),
            log_z_over_k: det::ln(z.abs()) - det::ln(k),
            log_pressure: det::ln(coordinates.omega) + det::ln(wall.fluid().density_kg_m3())
                - det::ln(4.0 * PI),
            log_gain,
            tail_extent: 0.0,
            split: None,
        })
    }

    // q, zeta and the absolute transformed measure q*dq/|zeta|.
    fn point(&self, region: usize, x: f64) -> (f64, C64, f64) {
        if region == 0 {
            let (sin, cos) = det::sin_cos(FRAC_PI_2 * x);
            let q = self.k * sin;
            let zeta = if x == 1.0 { 0.0 } else { self.k * cos };
            (q, C64::from_re(zeta), FRAC_PI_2 * q)
        } else {
            let (lower, upper) = match (region, self.split) {
                (1, Some(split)) => (0.0, split),
                (_, Some(split)) => (split, self.tail_extent),
                (_, None) => (0.0, self.tail_extent),
            };
            let t = lower + (upper - lower) * x;
            (det::hypot(self.k, t), C64::new(0.0, t), upper - lower)
        }
    }

    fn coefficient(&self, zeta: C64) -> C64 {
        let magnitude = zeta.abs();
        if magnitude == 0.0 {
            return -C64::ONE;
        }
        let product_log = self.log_z_over_k + det::ln(magnitude);
        let scale = product_log.max(0.0);
        let direction = self.unit_z * C64::new(zeta.re / magnitude, zeta.im / magnitude);
        let product = direction.scale(det::exp(product_log - scale));
        let reference = C64::from_re(det::exp(-scale));
        (product - reference) / (product + reference)
    }

    fn prepare_tail(
        &mut self,
        limits: QuadratureLimits,
        degrees: [u32; 2],
        log_scales: [f64; 2],
        channels: [usize; 2],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<[f64; 2], AeroacError> {
        let (extent, bounds) = prepare_reflected_tail(
            self.coordinates.height.to_f64(),
            self.k,
            limits,
            degrees,
            log_scales.map(|scale| scale + self.log_gain),
            channels,
            cancelled,
        )?;
        self.tail_extent = extent;
        let peak = self.k / self.z.abs();
        self.split = (peak > 0.0 && peak < extent).then_some(peak);
        Ok(bounds)
    }

    // Share the original outgoing Green weight with analytic source/receiver
    // derivatives. The physical prefactors belong to each source's channels.
    fn sample_weight(
        &self,
        region: usize,
        x: f64,
    ) -> Result<ReflectionSpectralWeight, AeroacError> {
        let (q, zeta, measure) = self.point(region, x);
        let r = self.coefficient(zeta);
        let magnitude = r.abs();
        if !magnitude.is_finite() {
            return Err(nonfinite("wall spectral reflection"));
        }
        if magnitude == 0.0 {
            return Ok(ReflectionSpectralWeight {
                q,
                zeta,
                unit: C64::ZERO,
                log_magnitude: 0.0,
            });
        }
        let phase = self.coordinates.height * Dd::from_f64(zeta.re);
        let decay = (self.coordinates.height * Dd::from_f64(zeta.im)).to_f64();
        if !phase.hi.is_finite() || !phase.lo.is_finite() || !decay.is_finite() {
            return Err(nonfinite("wall spectral propagation"));
        }
        let (sin, cos) = det::sin_cos_dd(phase);
        let direction = if region == 0 {
            C64::ONE
        } else {
            C64::new(0.0, -1.0)
        };
        Ok(ReflectionSpectralWeight {
            q,
            zeta,
            unit: C64::new(cos.to_f64(), sin.to_f64())
                * direction
                * C64::new(r.re / magnitude, r.im / magnitude),
            log_magnitude: det::ln(measure) + det::ln(magnitude) - decay,
        })
    }
}

impl Spectrum<3> for WallKernel {
    fn regions(&self) -> usize {
        if self.split.is_some() { 3 } else { 2 }
    }
    fn extent(&self) -> (f64, f64) {
        let cutoff = det::hypot(self.k, self.tail_extent);
        (cutoff, cutoff * self.coordinates.range)
    }
    fn resolved(&self, region: usize, lower: f64, upper: f64) -> bool {
        let (q0, z0, _) = self.point(region, lower);
        let (q1, z1, _) = self.point(region, upper);
        let (_, zm, _) = self.point(region, lower + 0.5 * (upper - lower));
        let r0 = self.coefficient(z0);
        let r1 = self.coefficient(z1);
        let rm = self.coefficient(zm);
        // The evanescent branches split at the unique amplitude maximum, so
        // their samples cannot straddle an unseen narrow surface-mode peak.
        // Resolve fractional changes there: an absolute amplitude threshold
        // would require work proportional to the (unbounded) passive gain.
        // The embedded quadrature still controls the absolute field errors.
        let scale = r0.abs().max(r1.abs()).max(rm.abs()).max(1.0);
        let h = self.coordinates.height.to_f64();
        scale.is_finite()
            && (q1 - q0).abs() * self.coordinates.range + h * (z1.re - z0.re).abs() <= FRAC_PI_2
            && h * (z1.im - z0.im).abs() <= 2.0
            && (r0 - rm).abs() <= 0.125 * scale
            && (r1 - rm).abs() <= 0.125 * scale
    }
    fn sample(&self, region: usize, x: f64) -> Result<[C64; 3], AeroacError> {
        let weight = self.sample_weight(region, x)?;
        let arg = weight.q * self.coordinates.range;
        let normal_abs = weight.zeta.abs();
        let normal = if normal_abs == 0.0 {
            C64::ZERO
        } else {
            C64::new(-weight.zeta.re / normal_abs, -weight.zeta.im / normal_abs)
        };
        let log_velocity = -det::ln(4.0 * PI);
        let values = [
            weight.scaled(C64::ONE, self.log_pressure, j0(arg)),
            weight.scaled(
                C64::new(0.0, 1.0),
                log_velocity + det::ln(weight.q),
                j1(arg),
            ),
            weight.scaled(normal, log_velocity + det::ln(normal_abs), j0(arg)),
        ];
        if values.iter().any(|v| !finite(*v)) {
            return Err(nonfinite("wall reflected spectrum"));
        }
        Ok(values)
    }
}
