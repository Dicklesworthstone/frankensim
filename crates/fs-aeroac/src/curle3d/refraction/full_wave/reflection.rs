//! Incident-side spherical reflection at the original two-fluid plane.
//!
//! For source/receiver clearances a,b in the incident fluid, H=a+b and
//! `zeta_j=sqrt(k_j²-q²)` on the outgoing branch, pressure/normal-velocity
//! continuity gives `R=(rho2*zeta1-rho1*zeta2)/(rho2*zeta1+rho1*zeta2)`.
//! With `W=q*dq*R*exp(i*zeta1*H)/zeta1`, the reflected monopole coefficients
//! are `omega*rho1/(4*pi) integral W*J0`, `i/(4*pi) integral W*q*J1`, and
//! `-1/(4*pi) integral W*zeta1*J0` for pressure, radial and normal velocity.
//! The minus normal sign is physical: the returning field travels away from
//! the plane. The parent cites the outgoing Weyl identity and interface law.
//!
//! All propagating, critical, total-internal-reflection and evanescent
//! components are retained. The incident grazing quotient is cancelled in
//! each branch map before evaluation. The original homogeneous direct field
//! is added coherently only for total observations; intensity is formed after
//! adding pressure and velocity, including their interference cross terms.
//!
//! Numerical tolerances and diagnostics control the REFLECTED transfer. The
//! total inherits its absolute error estimates, but relative accuracy is not
//! guaranteed at a direct/reflected cancellation. Roundoff, geometry and the
//! native Bessel accuracy scope remain outside these numerical estimates.

use super::{
    ComplexSum, Dd, FullWaveMonopoleDiagnostics, FullWaveMonopoleLimits, PlanarFluidInterface,
    QuadratureLimits, SpectralPoint, Spectrum, admit, finite, integrate_spectrum, invalid, j0, j1,
    log_weighted, minus_i, nonfinite, norm, poll, positive, spectral_point,
};
use crate::AeroacError;
use crate::curle3d::{CoherentCompactField, CompactSourcePhasor, PreparedCompactSourceObserver};
use crate::observation::HarmonicFieldLine;
use core::f64::consts::{FRAC_PI_2, PI, TAU};
use fs_math::{c64::C64, det};

mod force;
pub use force::{
    FULL_WAVE_REFLECTED_FORCE_SCOPE, FULL_WAVE_TOTAL_FORCE_SCOPE, PreparedReflectedFullWaveForce,
};
pub(super) use force::{apply as apply_force, cartesian as force_cartesian, sample_force_spectrum};

mod quadrupole;
pub use quadrupole::{
    FULL_WAVE_REFLECTED_QUADRUPOLE_SCOPE, FULL_WAVE_TOTAL_QUADRUPOLE_SCOPE,
    PreparedReflectedFullWaveQuadrupole,
};
pub(super) use quadrupole::{apply as apply_quadrupole, quadrupole_channels};

/// The returned field excludes the direct source contribution.
pub const FULL_WAVE_REFLECTED_MONOPOLE_SCOPE: &str = "Prescribed stationary compact 3D volume \
source; reflected spherical-wave field in the incident fluid at one infinite planar \
boundary between original stationary homogeneous lossless fluids, including critical, \
reactive and evanescent components. Numerical quadrature estimates and analytic spectral \
tail bounds are not interval certificates. No source feedback, finite edges, absorption, \
moving fluids or measured-SPL qualification.";

/// Pressure and velocity are added before computing the total energy flux.
pub const FULL_WAVE_TOTAL_MONOPOLE_SCOPE: &str = "Prescribed stationary compact 3D volume \
source; coherent direct plus reflected spherical-wave field in the incident fluid at one \
infinite planar boundary between original stationary homogeneous lossless fluids. \
Reflected-transfer error estimates do not guarantee relative accuracy at interference \
nulls. No singular direct self-field, source feedback, finite edges, absorption, moving \
fluids, interval certificate or measured-SPL qualification.";

/// Separate cached reflected and coherent total responses on the source clock.
/// Preparation bounds quadrature work; evaluation has no allocation, quadrature
/// or transcendental work. At a coincident source/receiver only the reflected
/// field exists. It is not a self-consistent source-loading or feedback model.
#[derive(Clone, Debug)]
pub struct PreparedReflectedFullWaveMonopole {
    interface: PlanarFluidInterface,
    frequency_hz: f64,
    source_m: [f64; 3],
    observer_m: [f64; 3],
    reflected: [C64; 4],
    total: Option<[C64; 4]>,
    diagnostics: FullWaveMonopoleDiagnostics,
}

impl PreparedReflectedFullWaveMonopole {
    /// Source and observer must lie strictly in the original incident fluid.
    /// The original monopole limits apply to the reflected transfer; the direct
    /// homogeneous contribution is analytic. Exact equal sound speeds use the
    /// physical density-jump image, including its full induction field.
    pub fn new(
        interface: PlanarFluidInterface,
        frequency_hz: f64,
        source_m: [f64; 3],
        observer_m: [f64; 3],
        limits: FullWaveMonopoleLimits,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, AeroacError> {
        poll(cancelled)?;
        admit(limits)?;
        let coordinates =
            ReflectionCoordinates::new(interface, frequency_hz, source_m, observer_m)?;
        let rho = interface.incident_state().density_kg_m3();
        let (cylindrical, diagnostics) = if interface.incident_state().sound_speed_m_s()
            == interface.receiving_state().sound_speed_m_s()
        {
            (
                equal_speed_monopole(interface, coordinates)?,
                FullWaveMonopoleDiagnostics {
                    evaluations: 0,
                    maximum_depth_reached: 0,
                    retained_panels: 0,
                    transverse_wavenumber_cutoff_rad_m: 0.0,
                    maximum_bessel_argument: 0.0,
                    pressure_quadrature_error_pa_per_m3_s: 0.0,
                    velocity_quadrature_error_per_m2: [0.0; 2],
                    pressure_tail_bound_pa_per_m3_s: 0.0,
                    velocity_tail_bound_per_m2: [0.0; 2],
                    exact_equal_speed: true,
                },
            )
        } else {
            let mut kernel = ReflectionKernel::new(interface, coordinates)?;
            let bound = kernel.prepare_tail(
                limits.into(),
                [0, 1],
                [
                    det::ln(coordinates.omega) + det::ln(rho) - det::ln(4.0 * PI),
                    -det::ln(4.0 * PI),
                ],
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
            (
                result.value,
                FullWaveMonopoleDiagnostics {
                    evaluations: result.evaluations,
                    maximum_depth_reached: result.maximum_depth_reached,
                    retained_panels: result.retained_panels,
                    transverse_wavenumber_cutoff_rad_m: result.cutoff,
                    maximum_bessel_argument: result.maximum_argument,
                    pressure_quadrature_error_pa_per_m3_s: result.error[0],
                    velocity_quadrature_error_per_m2: [result.error[1], result.error[2]],
                    pressure_tail_bound_pa_per_m3_s: tail[0],
                    velocity_tail_bound_per_m2: [tail[1], tail[2]],
                    exact_equal_speed: false,
                },
            )
        };
        let [pressure, radial, normal] = cylindrical;
        let n = interface.receiving_normal();
        let t = coordinates.tangent;
        let reflected = [
            pressure,
            radial.scale(t[0]) + normal.scale(n[0]),
            radial.scale(t[1]) + normal.scale(n[1]),
            radial.scale(t[2]) + normal.scale(n[2]),
        ];
        if reflected.iter().any(|v| !finite(*v)) {
            return Err(nonfinite("reflected monopole transfer"));
        }
        let total = if source_m == observer_m {
            None
        } else {
            let original = PreparedCompactSourceObserver::new(
                interface.incident_state(),
                frequency_hz,
                observer_m,
                &[source_m],
                1,
            )?;
            let direct =
                original.evaluate_harmonic_line(&[CompactSourcePhasor::monopole(C64::ONE)])?;
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
            if total.iter().any(|v| !finite(*v)) {
                return Err(nonfinite("total incident-side monopole transfer"));
            }
            Some(total)
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
    /// Numerical evidence for the reflected transfer only; see the module scope.
    pub const fn diagnostics(&self) -> &FullWaveMonopoleDiagnostics {
        &self.diagnostics
    }
    /// False only for coincident original source and observer positions.
    pub const fn total_field_available(&self) -> bool {
        self.total.is_some()
    }

    /// Original complex peak volume flow Q [m³/s], exp(-i omega t).
    pub fn evaluate_reflected(
        &self,
        volume_flow_m3_s: C64,
    ) -> Result<CoherentCompactField, AeroacError> {
        field(
            apply(self.reflected, volume_flow_m3_s)?,
            FULL_WAVE_REFLECTED_MONOPOLE_SCOPE,
        )
    }

    /// Sum the original direct and reflected fields before computing intensity.
    /// Refuses the singular direct field at a coincident source/receiver.
    pub fn evaluate_total(
        &self,
        volume_flow_m3_s: C64,
    ) -> Result<CoherentCompactField, AeroacError> {
        field(
            apply(self.total.ok_or_else(singular_total)?, volume_flow_m3_s)?,
            FULL_WAVE_TOTAL_MONOPOLE_SCOPE,
        )
    }

    pub fn evaluate_reflected_harmonic_line(
        &self,
        volume_flow_m3_s: C64,
    ) -> Result<HarmonicFieldLine<3>, AeroacError> {
        Ok(line(
            self.frequency_hz,
            apply(self.reflected, volume_flow_m3_s)?,
        ))
    }

    pub fn evaluate_total_harmonic_line(
        &self,
        volume_flow_m3_s: C64,
    ) -> Result<HarmonicFieldLine<3>, AeroacError> {
        Ok(line(
            self.frequency_hz,
            apply(self.total.ok_or_else(singular_total)?, volume_flow_m3_s)?,
        ))
    }
}

pub(super) fn singular_total() -> AeroacError {
    invalid("total incident-side field is singular at the prescribed point source")
}

pub(super) fn apply(transfer: [C64; 4], amplitude: C64) -> Result<[C64; 4], AeroacError> {
    if !finite(amplitude) {
        return Err(nonfinite("reflected volume-flow phasor"));
    }
    let values = transfer.map(|value| value * amplitude);
    if values.iter().any(|v| !finite(*v)) {
        return Err(nonfinite("reflected monopole field"));
    }
    Ok(values)
}

pub(super) fn field(
    amplitudes: [C64; 4],
    scope: &'static str,
) -> Result<CoherentCompactField, AeroacError> {
    let [pressure, ux, uy, uz] = amplitudes;
    let velocity = [ux, uy, uz];
    let intensity = velocity.map(|v| pressure.scale(0.5) * v.conj());
    if intensity.iter().any(|v| !finite(*v)) {
        return Err(nonfinite("incident-side acoustic intensity"));
    }
    Ok(CoherentCompactField {
        pressure_pa: pressure,
        particle_velocity_m_s: velocity,
        active_intensity_w_m2: intensity.map(|v| v.re),
        reactive_intensity_w_m2: intensity.map(|v| v.im),
        scope,
    })
}

pub(super) fn line(frequency_hz: f64, amplitudes: [C64; 4]) -> HarmonicFieldLine<3> {
    let [pressure, ux, uy, uz] = amplitudes;
    HarmonicFieldLine {
        frequency_hz,
        pressure_pa: pressure,
        particle_velocity_m_s: [ux, uy, uz],
    }
}

#[derive(Clone, Copy)]
pub(super) struct ReflectionCoordinates {
    pub(super) omega: f64,
    pub(super) height: Dd,
    pub(super) range: f64,
    pub(super) tangent: [f64; 3],
}

impl ReflectionCoordinates {
    pub(super) fn new(
        interface: PlanarFluidInterface,
        frequency_hz: f64,
        source: [f64; 3],
        observer: [f64; 3],
    ) -> Result<Self, AeroacError> {
        let omega = TAU * frequency_hz;
        if !positive(omega) {
            return Err(invalid("reflected monopole frequency"));
        }
        if source.iter().chain(&observer).any(|v| !v.is_finite()) {
            return Err(nonfinite("incident-side acoustic positions"));
        }
        let a = -interface.distance(source);
        let b = -interface.distance(observer);
        let height = a + b;
        if !positive(a.to_f64()) || !positive(b.to_f64()) || !positive(height.to_f64()) {
            return Err(invalid(
                "source and observer must lie strictly in the incident fluid",
            ));
        }
        let lateral = core::array::from_fn(|i| {
            (Dd::from_f64(observer[i])
                - Dd::from_f64(source[i])
                - (a - b) * Dd::from_f64(interface.receiving_normal()[i]))
            .to_f64()
        });
        let range = norm(lateral);
        if !range.is_finite() {
            return Err(nonfinite("incident-side lateral separation"));
        }
        let tangent = if range == 0.0 {
            [0.0; 3]
        } else {
            lateral.map(|v| v / range)
        };
        Ok(Self {
            omega,
            height,
            range,
            tangent,
        })
    }
}

fn exact_reflection(interface: PlanarFluidInterface) -> f64 {
    let rho1 = interface.incident_state().density_kg_m3();
    let rho2 = interface.receiving_state().density_kg_m3();
    let scale = rho1.max(rho2);
    let a = rho1 / scale;
    let b = rho2 / scale;
    (b - a) / (b + a)
}

fn equal_speed_monopole(
    interface: PlanarFluidInterface,
    coordinates: ReflectionCoordinates,
) -> Result<[C64; 3], AeroacError> {
    let k = coordinates.omega / interface.incident_state().sound_speed_m_s();
    if !positive(k) {
        return Err(nonfinite("reflected acoustic wavenumber"));
    }
    let reflection = exact_reflection(interface);
    if reflection == 0.0 {
        return Ok([C64::ZERO; 3]);
    }
    let h = coordinates.height.to_f64();
    let distance = det::hypot(coordinates.range, h);
    let phase = Dd::from_f64(coordinates.omega) * Dd::from_f64(distance)
        / Dd::from_f64(interface.incident_state().sound_speed_m_s());
    if !positive(distance) || !phase.hi.is_finite() || !phase.lo.is_finite() {
        return Err(nonfinite("reflected image distance or phase"));
    }
    let (sin, cos) = det::sin_cos_dd(phase);
    let phasor = C64::new(cos.to_f64(), sin.to_f64());
    let log_green = -det::ln(4.0 * PI) - det::ln(distance);
    let pressure = log_weighted(
        minus_i(phasor),
        log_green
            + det::ln(coordinates.omega)
            + det::ln(interface.incident_state().density_kg_m3()),
        reflection,
    );
    let green = phasor.scale((0.25 / PI) / distance);
    let radial = green * C64::new(1.0 / distance, -k).scale(reflection);
    Ok([
        pressure,
        radial.scale(coordinates.range / distance),
        radial.scale(-h / distance),
    ])
}

/// q*dq*R*exp(i*zeta1*H)/zeta1, with physical scale retained before derivatives.
pub(super) struct ReflectionSpectralWeight {
    pub(super) q: f64,
    pub(super) zeta: C64,
    pub(super) unit: C64,
    pub(super) log_magnitude: f64,
}

impl ReflectionSpectralWeight {
    pub(super) fn scaled(&self, direction: C64, log_scale: f64, real_weight: f64) -> C64 {
        log_weighted(
            self.unit * direction,
            self.log_magnitude + log_scale,
            real_weight,
        )
    }
}

struct ReflectionKernel {
    coordinates: ReflectionCoordinates,
    kmin: f64,
    kmax: f64,
    gap: f64,
    incident_has_larger_k: bool,
    log_density: [f64; 2],
    tail_extent: f64,
}

impl ReflectionKernel {
    fn new(
        interface: PlanarFluidInterface,
        coordinates: ReflectionCoordinates,
    ) -> Result<Self, AeroacError> {
        let k = [interface.incident_state(), interface.receiving_state()]
            .map(|fluid| coordinates.omega / fluid.sound_speed_m_s());
        if k.iter().any(|v| !positive(*v)) || k[0] == k[1] {
            return Err(nonfinite("distinct reflected acoustic wavenumbers"));
        }
        let kmin = k[0].min(k[1]);
        let kmax = k[0].max(k[1]);
        let gap = det::sqrt(kmax - kmin) * det::sqrt(kmax) * det::sqrt(1.0 + kmin / kmax);
        if !positive(gap) {
            return Err(nonfinite("reflected critical-wavenumber gap"));
        }
        Ok(Self {
            coordinates,
            kmin,
            kmax,
            gap,
            incident_has_larger_k: k[0] > k[1],
            log_density: [
                det::ln(interface.incident_state().density_kg_m3()),
                det::ln(interface.receiving_state().density_kg_m3()),
            ],
            tail_extent: 0.0,
        })
    }

    fn point(&self, region: usize, x: f64) -> SpectralPoint {
        spectral_point(
            self.kmin,
            self.kmax,
            self.gap,
            self.incident_has_larger_k,
            self.tail_extent,
            region,
            x,
        )
    }

    // Positive original densities and outgoing zeta make these admittance
    // products lie in the closed first quadrant. Scale their logarithms before
    // addition: neither a huge product nor a tiny density ratio must materialize.
    fn coefficient(&self, zeta: [C64; 2]) -> C64 {
        let logs = [
            self.log_density[1] + det::ln(zeta[0].abs()),
            self.log_density[0] + det::ln(zeta[1].abs()),
        ];
        let scale = logs[0].max(logs[1]);
        let terms = core::array::from_fn::<_, 2, _>(|i| {
            let direction = if zeta[i].re > 0.0 {
                C64::ONE
            } else if zeta[i].im > 0.0 {
                C64::new(0.0, 1.0)
            } else {
                C64::ZERO
            };
            direction.scale(det::exp(logs[i] - scale))
        });
        (terms[0] - terms[1]) / (terms[0] + terms[1])
    }

    fn sample_weight(
        &self,
        region: usize,
        x: f64,
    ) -> Result<ReflectionSpectralWeight, AeroacError> {
        let point = self.point(region, x);
        let reflection = self.coefficient(point.zeta);
        let magnitude = reflection.abs();
        if !magnitude.is_finite() {
            return Err(nonfinite("fluid reflection coefficient"));
        }
        if magnitude == 0.0 {
            return Ok(ReflectionSpectralWeight {
                q: point.q,
                zeta: point.zeta[0],
                unit: C64::ZERO,
                log_magnitude: 0.0,
            });
        }
        let phase = self.coordinates.height * Dd::from_f64(point.zeta[0].re);
        let attenuation = (self.coordinates.height * Dd::from_f64(point.zeta[0].im)).to_f64();
        if !phase.hi.is_finite() || !phase.lo.is_finite() || !attenuation.is_finite() {
            return Err(nonfinite("reflected spectral phase"));
        }
        let (sin, cos) = det::sin_cos_dd(phase);
        let phasor = C64::new(cos.to_f64(), sin.to_f64());
        // Cancel q*dq/zeta1 analytically at BOTH original critical circles.
        // All remaining scalar factors are nonnegative on the transformed roots.
        let (direction, log_measure) = match (region, self.incident_has_larger_k) {
            (0, false) => (C64::ONE, det::ln(FRAC_PI_2) + det::ln(point.q)),
            (0, true) => (
                C64::ONE,
                det::ln(FRAC_PI_2) + det::ln(point.q) + det::ln(point.zeta[1].abs())
                    - det::ln(point.zeta[0].abs()),
            ),
            (1, false) => (
                C64::new(0.0, -1.0),
                det::ln(FRAC_PI_2) + det::ln(point.zeta[1].abs()),
            ),
            (1, true) => (C64::ONE, det::ln(FRAC_PI_2) + det::ln(point.zeta[1].abs())),
            (_, true) => (C64::new(0.0, -1.0), det::ln(self.tail_extent)),
            (_, false) => (
                C64::new(0.0, -1.0),
                det::ln(self.tail_extent) + det::ln(point.measure[0])
                    - det::ln(point.zeta[0].abs()),
            ),
        };
        Ok(ReflectionSpectralWeight {
            q: point.q,
            zeta: point.zeta[0],
            unit: phasor
                * direction
                * C64::new(reflection.re / magnitude, reflection.im / magnitude),
            log_magnitude: log_measure + det::ln(magnitude) - attenuation,
        })
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
            self.kmax,
            limits,
            degrees,
            log_scales,
            channels,
            cancelled,
        )?;
        self.tail_extent = extent;
        Ok(bounds)
    }
}

/// Bound a positive spectral tail C*integral exp(-h*t)*(t+kmax)^degree dt.
/// Each caller supplies its original attenuation clearance and prefactor.
/// Above every critical circle, a reflected lossless-fluid amplitude has
/// magnitude at most one; transmitted kernels include their admittance bound
/// in the prefactor. Degrees zero through three cover Q, force and stress.
pub(super) fn prepare_reflected_tail(
    h: f64,
    kmax: f64,
    limits: QuadratureLimits,
    degrees: [u32; 2],
    log_scales: [f64; 2],
    channels: [usize; 2],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(f64, [f64; 2]), AeroacError> {
    if degrees.iter().any(|degree| *degree > 3) {
        return Err(invalid("evanescent tail polynomial degree is unsupported"));
    }
    let hk = h * kmax;
    if !hk.is_finite() {
        return Err(nonfinite("reflected clearance wavenumber"));
    }
    // For q²=kmax²+t², q*dq/|zeta1| <= dt, |R|<=1,
    // exp(-H*Im(zeta1)) <= exp(-H*t), and q,|zeta1| <= t+kmax.
    let log_h = det::ln(h);
    let goals = [limits.absolute_pressure, limits.absolute_velocity].map(det::ln);
    let goals = core::array::from_fn::<_, 2, _>(|i| {
        goals[i] - det::ln(4.0 * det::sqrt(channels[i] as f64))
    });
    let logs = |s: f64| {
        core::array::from_fn::<_, 2, _>(|i| {
            let x = s + hk + 1.0;
            let polynomial = match degrees[i] {
                0 => 0.0,
                1 => det::ln(x),
                // (s+hk)^2 + 2(s+hk) + 2 = x²+1; no x² overflow.
                2 => 2.0 * det::ln(x) + det::ln(1.0 + (1.0 / x) * (1.0 / x)),
                // (s+hk)^3 + 3(s+hk)^2 + 6(s+hk) + 6 = x³+3x+2.
                // Ratios retain the cubic envelope without forming x³.
                3 => {
                    let inverse = 1.0 / x;
                    3.0 * det::ln(x) + det::ln(1.0 + inverse * inverse * (3.0 + 2.0 * inverse))
                }
                _ => unreachable!("tail degrees admitted above"),
            };
            log_scales[i] - (f64::from(degrees[i]) + 1.0) * log_h - s + polynomial
        })
    };
    let (mut low, mut high) = (0.0, 1.0);
    let mut admitted = false;
    for _ in 0..32 {
        poll(cancelled)?;
        let bound = logs(high);
        if bound[0] <= goals[0] && bound[1] <= goals[1] {
            admitted = true;
            break;
        }
        low = high;
        high *= 2.0;
    }
    if !admitted {
        return Err(invalid(
            "reflected evanescent tail cannot meet requested tolerances",
        ));
    }
    for _ in 0..24 {
        poll(cancelled)?;
        let middle = low + (high - low) * 0.5;
        let bound = logs(middle);
        if bound[0] <= goals[0] && bound[1] <= goals[1] {
            high = middle;
        } else {
            low = middle;
        }
    }
    let tail_extent = high / h;
    if !positive(tail_extent) {
        return Err(nonfinite("reflected evanescent cutoff"));
    }
    Ok((tail_extent, logs(high).map(det::exp)))
}

impl Spectrum<3> for ReflectionKernel {
    fn extent(&self) -> (f64, f64) {
        let cutoff = det::hypot(self.kmax, self.tail_extent);
        (cutoff, cutoff * self.coordinates.range)
    }

    fn resolved(&self, region: usize, lower: f64, upper: f64) -> bool {
        let lo = self.point(region, lower);
        let hi = self.point(region, upper);
        let h = self.coordinates.height.to_f64();
        let phase = (hi.q - lo.q).abs() * self.coordinates.range
            + h * (hi.zeta[0].re - lo.zeta[0].re).abs();
        let decay = h * (hi.zeta[0].im - lo.zeta[0].im).abs();
        // R is monotone real on the two outer regions and traverses a monotone
        // semicircle between the critical circles. Endpoint separation bounds
        // all within-panel coefficient variation, including narrow transitions.
        phase <= FRAC_PI_2
            && decay <= 2.0
            && (self.coefficient(hi.zeta) - self.coefficient(lo.zeta)).abs() <= 0.25
    }

    fn sample(&self, region: usize, x: f64) -> Result<[C64; 3], AeroacError> {
        let weight = self.sample_weight(region, x)?;
        let argument = weight.q * self.coordinates.range;
        if !argument.is_finite() {
            return Err(nonfinite("reflected Bessel argument"));
        }
        let jzero = j0(argument);
        let log_velocity = -det::ln(4.0 * PI);
        let log_pressure = det::ln(self.coordinates.omega) + self.log_density[0] + log_velocity;
        let magnitude = weight.zeta.abs();
        let normal = if magnitude == 0.0 {
            C64::ZERO
        } else {
            C64::new(-weight.zeta.re / magnitude, -weight.zeta.im / magnitude)
        };
        let values = [
            weight.scaled(C64::ONE, log_pressure, jzero),
            weight.scaled(
                C64::new(0.0, 1.0),
                log_velocity + det::ln(weight.q),
                j1(argument),
            ),
            weight.scaled(normal, log_velocity + det::ln(magnitude), jzero),
        ];
        if values.iter().any(|v| !finite(*v)) {
            return Err(nonfinite("reflected monopole spectrum"));
        }
        Ok(values)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use fs_material::fluid_acoustics::FluidAcousticState;

    fn kernel(rho: [f64; 2], speed: [f64; 2], hz: f64, height: f64) -> ReflectionKernel {
        let interface = PlanarFluidInterface::new(
            FluidAcousticState::try_new(rho[0], speed[0]).unwrap(),
            FluidAcousticState::try_new(rho[1], speed[1]).unwrap(),
            [0.0; 3],
            [0.0, 0.0, 1.0],
        )
        .unwrap();
        let coordinates = ReflectionCoordinates::new(
            interface,
            hz,
            [0.0, 0.0, -height * 0.5],
            [0.0, 0.0, -height * 0.5],
        )
        .unwrap();
        let mut kernel = ReflectionKernel::new(interface, coordinates).unwrap();
        kernel.tail_extent = 400.0;
        kernel
    }

    #[test]
    fn critical_quotients_and_strong_evanescence_retain_representable_fields() {
        for speed in [[300.0, 600.0], [600.0, 300.0]] {
            let kernel = kernel([1.2, 2.4], speed, 700.0, 0.05);
            for region in 0..3 {
                for x in [0.0, 1.0] {
                    assert!(kernel.sample(region, x).unwrap().iter().all(|v| finite(*v)));
                }
            }
        }
        let attenuated = kernel([1e250, 2e250], [1.0, 2.0], 1.0, 2.0);
        assert_eq!(det::exp(-800.0), 0.0);
        let pressure = attenuated.sample(2, 1.0).unwrap()[0];
        assert!(pressure.abs().is_finite() && pressure.abs() > 0.0);
    }

    #[test]
    fn near_equal_speeds_still_resolve_the_original_critical_reflection_transition() {
        let kernel = kernel([1.2, 2.4], [343.0, 343.01], 1000.0, 0.08);
        let lo = kernel.point(1, 0.0);
        let hi = kernel.point(1, 1.0);
        assert!(
            kernel.coordinates.height.to_f64() * (hi.zeta[0].re - lo.zeta[0].re).abs() < FRAC_PI_2
        );
        assert!(kernel.coordinates.height.to_f64() * (hi.zeta[0].im - lo.zeta[0].im).abs() < 2.0);
        assert!(!kernel.resolved(1, 0.0, 1.0));
        assert!((kernel.coefficient(hi.zeta) - kernel.coefficient(lo.zeta)).abs() > 1.9);
    }

    #[test]
    fn cubic_tail_matches_independent_integration_by_parts() {
        let limits = QuadratureLimits {
            absolute_pressure: 1e-7,
            absolute_velocity: 1e-9,
            relative: 0.0,
            maximum_evaluations: 1024,
            maximum_depth: 24,
        };
        let scales: [f64; 2] = [0.7, 1.3];
        let channels: [usize; 2] = [5, 14];
        for h in [0.005, 0.2, 3.0] {
            let kmax = 80.0;
            let (cutoff, bounds) = prepare_reflected_tail(
                h,
                kmax,
                limits,
                [2, 3],
                scales.map(f64::ln),
                channels,
                &mut || false,
            )
            .unwrap();
            // Integration by parts, directly in physical t units; this does
            // not use the dimensionless shifted/log polynomial above.
            let boundary = (-h * cutoff).exp() / h;
            let mut integral = boundary;
            let mut power = 1.0;
            for degree in 1..=3 {
                power *= cutoff + kmax;
                integral = boundary * power + f64::from(degree) * integral / h;
                if degree >= 2 {
                    let index = (degree - 2) as usize;
                    let expected = scales[index] * integral;
                    assert!((bounds[index] - expected).abs() <= 3e-12 * expected);
                }
            }
            for (i, tolerance) in [limits.absolute_pressure, limits.absolute_velocity]
                .into_iter()
                .enumerate()
            {
                assert!(bounds[i] * 4.0 * (channels[i] as f64).sqrt() <= tolerance);
            }
        }
        assert!(
            prepare_reflected_tail(0.1, 80.0, limits, [2, 4], [0.0; 2], channels, &mut || false,)
                .is_err()
        );
    }
}
