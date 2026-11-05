//! Full source spectrum through one finite homogeneous fluid layer.
//!
//! Let zeta_j=sqrt(k_j²-q²) have nonnegative real/imaginary parts and
//! Y_j=zeta_j/rho_j. Pressure and normal-velocity continuity give
//! L=(Y_1+Y_3) cos(zeta_2 h)
//!   -i h (zeta_2²/rho_2 + rho_2 Y_1 Y_3) sinc(zeta_2 h).
//! With exterior E=exp(i zeta_1 a+i zeta_3 b),
//! p/Q=omega/(2 pi) integral q J0(qR) E/L dq. Its analytic receiver
//! derivatives supply BOTH velocity components. This sums every coherent
//! internal return and retains source evanescence and frustrated transmission.
//!
//! The original continuity equations are C. C. Mei, MIT Wave Propagation,
//! chapter 3, equations (1.18)--(1.25):
//! <https://web.mit.edu/1.138j/www/material/chap-3.pdf>.
//! The parent documents the outgoing Sommerfeld identity and source clock.
//! The sinc denominator is regular at internal grazing; multiplying it by
//! exp(-Im(zeta_2) h) avoids growing transfer matrices. An additional panel
//! resolution condition bounds denominator variation, so narrow coherent
//! resonances cannot pass solely through coincident quadrature samples.
//!
//! A positive-density layer slower than BOTH exteriors can support real guided
//! poles. That class is explicitly refused, not damped or silently sampled.
//! Finite lateral edges, absorption, moving fluids and source loading are absent.

use super::{
    Dd, FullWaveMonopoleDiagnostics, FullWaveMonopoleLimits, PlanarFluidInterface,
    PreparedFullWaveMonopole, QuadratureLimits, Spectrum, admit, finite, integrate_spectrum,
    invalid, j0, j1, log_weighted, minus_i, nonfinite, norm, poll, positive,
};
use crate::AeroacError;
use crate::curle3d::CoherentCompactField;
use crate::observation::HarmonicFieldLine;
use core::f64::consts::{FRAC_PI_2, TAU};
use fs_material::fluid_acoustics::FluidAcousticState;
use fs_math::{c64::C64, det};

mod force;
pub use force::{
    FULL_WAVE_LAYER_FORCE_SCOPE, FullWaveLayerForceDiagnostics, PreparedLayeredFullWaveForce,
};
pub mod reflection;

/// Applicability of the original finite-layer source-to-receiver transfer.
pub const FULL_WAVE_LAYER_MONOPOLE_SCOPE: &str = "Prescribed stationary compact 3D volume \
source; coherent spherical-wave transmission through one infinite planar homogeneous \
lossless fluid layer between original stationary lossless fluids, including internal \
returns, reactive fields and evanescence. Slower-than-both trapping layers refuse. \
Numerical quadrature estimates and analytic spectral-tail bounds are not interval \
certificates; no absorption, moving media, finite edges, source feedback or measured-SPL qualification.";

/// Two parallel faces and three original stationary fluid snapshots.
/// The first face passes through `first_face_point_m`; the receiving normal
/// points into the layer and toward the second face a distance h away.
#[derive(Clone, Copy, Debug)]
pub struct PlanarFluidLayer {
    exterior: PlanarFluidInterface,
    layer: FluidAcousticState,
    thickness_m: f64,
}

impl PlanarFluidLayer {
    /// Admit a nonnegative finite thickness. A nonzero layer must not be slower
    /// than both exteriors: c_layer >= min(c_incident,c_receiving). This excludes
    /// real guided-mode poles without inserting any fictitious attenuation.
    /// A zero-thickness layer is an exact original two-fluid interface.
    pub fn new(
        incident: FluidAcousticState,
        layer: FluidAcousticState,
        receiving: FluidAcousticState,
        first_face_point_m: [f64; 3],
        receiving_normal: [f64; 3],
        thickness_m: f64,
    ) -> Result<Self, AeroacError> {
        if !thickness_m.is_finite() || thickness_m < 0.0 {
            return Err(invalid(
                "fluid-layer thickness must be finite and nonnegative",
            ));
        }
        if thickness_m > 0.0
            && layer.sound_speed_m_s() < incident.sound_speed_m_s().min(receiving.sound_speed_m_s())
        {
            return Err(invalid(
                "slower-than-both fluid layer requires a guided-pole solver",
            ));
        }
        Ok(Self {
            exterior: PlanarFluidInterface::new(
                incident,
                receiving,
                first_face_point_m,
                receiving_normal,
            )?,
            layer,
            thickness_m,
        })
    }

    pub const fn incident_state(self) -> FluidAcousticState {
        self.exterior.incident_state()
    }
    pub const fn layer_state(self) -> FluidAcousticState {
        self.layer
    }
    pub const fn receiving_state(self) -> FluidAcousticState {
        self.exterior.receiving_state()
    }
    pub const fn first_face_point_m(self) -> [f64; 3] {
        self.exterior.point_m()
    }
    pub const fn receiving_normal(self) -> [f64; 3] {
        self.exterior.receiving_normal()
    }
    pub const fn thickness_m(self) -> f64 {
        self.thickness_m
    }

    fn reduces_to_interface(self) -> bool {
        self.thickness_m == 0.0
            || (self.incident_state() == self.layer && self.receiving_state() == self.layer)
    }
}

/// Numerical estimates and analytic omitted-tail bounds for unit volume flow.
/// Multiply absolute entries by |Q|. Numerical quadrature is estimated, not
/// certified; geometry, Bessel and floating-point errors are separate.
#[derive(Clone, Copy, Debug)]
pub struct FullWaveLayerMonopoleDiagnostics {
    pub evaluations: usize,
    pub maximum_depth_reached: usize,
    pub retained_panels: usize,
    pub transverse_wavenumber_cutoff_rad_m: f64,
    pub maximum_bessel_argument: f64,
    pub pressure_quadrature_error_pa_per_m3_s: f64,
    /// Radial and normal entries [1/m²] in the original face frame.
    pub velocity_quadrature_error_per_m2: [f64; 2],
    pub pressure_tail_bound_pa_per_m3_s: f64,
    pub velocity_tail_bound_per_m2: [f64; 2],
    /// An absent or completely identical layer used the original interface owner.
    /// Its original unequal-speed interface may still require quadrature.
    pub reduced_to_single_interface: bool,
}

impl From<FullWaveMonopoleDiagnostics> for FullWaveLayerMonopoleDiagnostics {
    fn from(value: FullWaveMonopoleDiagnostics) -> Self {
        Self {
            evaluations: value.evaluations,
            maximum_depth_reached: value.maximum_depth_reached,
            retained_panels: value.retained_panels,
            transverse_wavenumber_cutoff_rad_m: value.transverse_wavenumber_cutoff_rad_m,
            maximum_bessel_argument: value.maximum_bessel_argument,
            pressure_quadrature_error_pa_per_m3_s: value.pressure_quadrature_error_pa_per_m3_s,
            velocity_quadrature_error_per_m2: value.velocity_quadrature_error_per_m2,
            pressure_tail_bound_pa_per_m3_s: value.pressure_tail_bound_pa_per_m3_s,
            velocity_tail_bound_per_m2: value.velocity_tail_bound_per_m2,
            reduced_to_single_interface: true,
        }
    }
}

/// A complete transmitted spherical-wave field with four cached coefficients.
/// Preparation owns all bounded quadrature; applying Q or obtaining a native
/// harmonic observation line allocates nothing and performs fixed work.
#[derive(Clone, Debug)]
pub struct PreparedLayeredFullWaveMonopole {
    layer: PlanarFluidLayer,
    frequency_hz: f64,
    source_m: [f64; 3],
    observer_m: [f64; 3],
    transfer: [C64; 4],
    diagnostics: FullWaveLayerMonopoleDiagnostics,
}

impl PreparedLayeredFullWaveMonopole {
    /// Source and observer must be strictly outside opposite faces. The limits
    /// retain their original physical pressure/velocity units and shared node
    /// count. Internal histories are summed coherently in the Green function,
    /// never truncated into a finite list of rays or regularized by damping.
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
        let coordinates = LayerCoordinates::new(layer, frequency_hz, source_m, observer_m)?;
        let (transfer, diagnostics) = if layer.reduces_to_interface() {
            let original = PreparedFullWaveMonopole::new(
                layer.exterior,
                frequency_hz,
                source_m,
                observer_m,
                limits,
                cancelled,
            )?;
            (original.transfer, (*original.diagnostics()).into())
        } else {
            let mut kernel = LayerKernel::new(layer, coordinates)?;
            let rho3 = layer.receiving_state().density_kg_m3();
            let bound = kernel.prepare_tail(
                limits.into(),
                [0, 1],
                [
                    det::ln(coordinates.omega) - det::ln(TAU),
                    -det::ln(TAU) - det::ln(rho3),
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
            let [pressure, radial, normal] = result.value;
            let transfer = core::array::from_fn(|i| {
                if i == 0 {
                    pressure
                } else {
                    radial.scale(coordinates.tangent[i - 1])
                        + normal.scale(layer.receiving_normal()[i - 1])
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
            return Err(nonfinite("full-wave layered monopole transfer"));
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

    fn amplitudes(&self, volume_flow_m3_s: C64) -> Result<[C64; 4], AeroacError> {
        if !finite(volume_flow_m3_s) {
            return Err(nonfinite("layered volume-flow phasor"));
        }
        let values = self.transfer.map(|v| v * volume_flow_m3_s);
        if values.iter().any(|v| !finite(*v)) {
            return Err(nonfinite("layered monopole field"));
        }
        Ok(values)
    }

    /// Complex peak Q [m³/s] on the original exp(-i omega t) source clock.
    pub fn evaluate(&self, volume_flow_m3_s: C64) -> Result<CoherentCompactField, AeroacError> {
        let [pressure, ux, uy, uz] = self.amplitudes(volume_flow_m3_s)?;
        let velocity = [ux, uy, uz];
        let intensity = velocity.map(|v| pressure.scale(0.5) * v.conj());
        if intensity.iter().any(|v| !finite(*v)) {
            return Err(nonfinite("layered monopole intensity"));
        }
        Ok(CoherentCompactField {
            pressure_pa: pressure,
            particle_velocity_m_s: velocity,
            active_intensity_w_m2: intensity.map(|v| v.re),
            reactive_intensity_w_m2: intensity.map(|v| v.im),
            scope: FULL_WAVE_LAYER_MONOPOLE_SCOPE,
        })
    }

    /// All layer and exterior phases are already present; apply no extra delay.
    pub fn evaluate_harmonic_line(
        &self,
        volume_flow_m3_s: C64,
    ) -> Result<HarmonicFieldLine<3>, AeroacError> {
        let [pressure, ux, uy, uz] = self.amplitudes(volume_flow_m3_s)?;
        Ok(HarmonicFieldLine {
            frequency_hz: self.frequency_hz,
            pressure_pa: pressure,
            particle_velocity_m_s: [ux, uy, uz],
        })
    }
}

#[derive(Clone, Copy)]
struct LayerCoordinates {
    omega: f64,
    a: Dd,
    b: Dd,
    range: f64,
    tangent: [f64; 3],
}

impl LayerCoordinates {
    fn new(
        layer: PlanarFluidLayer,
        frequency_hz: f64,
        source: [f64; 3],
        observer: [f64; 3],
    ) -> Result<Self, AeroacError> {
        let omega = TAU * frequency_hz;
        if !positive(omega) {
            return Err(invalid("full-wave layer frequency"));
        }
        if source.iter().chain(&observer).any(|v| !v.is_finite()) {
            return Err(nonfinite("full-wave layer positions"));
        }
        let a = -layer.exterior.distance(source);
        let b = layer.exterior.distance(observer) - Dd::from_f64(layer.thickness_m);
        let normal_distance = a + b + Dd::from_f64(layer.thickness_m);
        if !positive(a.to_f64()) || !positive(b.to_f64()) || !positive(normal_distance.to_f64()) {
            return Err(invalid(
                "full-wave source and observer must lie strictly outside opposite layer faces",
            ));
        }
        let lateral = core::array::from_fn(|i| {
            (Dd::from_f64(observer[i])
                - Dd::from_f64(source[i])
                - normal_distance * Dd::from_f64(layer.receiving_normal()[i]))
            .to_f64()
        });
        let range = norm(lateral);
        if !range.is_finite() {
            return Err(nonfinite("full-wave layer lateral distance"));
        }
        let tangent = if range == 0.0 {
            [0.0; 3]
        } else {
            lateral.map(|v| v / range)
        };
        Ok(Self {
            omega,
            a,
            b,
            range,
            tangent,
        })
    }
}

#[derive(Clone, Copy)]
struct LayerPoint {
    q: f64,
    zeta: [C64; 3],
    // Product equals q*dq/dx; retain factors before logarithmic scaling.
    measure: [f64; 2],
}

/// Shared physical q*dq*E/L weight, retaining its scale before derivatives.
/// All logarithms belong to preparation; the final observation is cached.
struct LayerSpectralWeight {
    q: f64,
    zeta: [C64; 3],
    unit: C64,
    log_magnitude: f64,
}

impl LayerSpectralWeight {
    /// `direction` is a unit complex phase (or exact zero); `log_scale` keeps
    /// physical coefficients and derivative magnitudes attached until exp.
    fn scaled(&self, direction: C64, log_scale: f64, real_weight: f64) -> C64 {
        log_weighted(
            self.unit * direction,
            self.log_magnitude + log_scale,
            real_weight,
        )
    }
}

struct LayerKernel {
    layer: PlanarFluidLayer,
    coordinates: LayerCoordinates,
    k: [f64; 3],
    critical: [f64; 3],
    critical_count: usize,
    // Absolute square-root differences between medium k and each critical k.
    gaps: [[f64; 3]; 3],
    // Scaled coefficients 1/rho1,1/rho3,1/rho2,rho2/(rho1rho3).
    coefficients: [f64; 4],
    log_density_scale: f64,
    common_speed: bool,
    tail_extent: f64,
}

fn wave_gap(a: f64, b: f64) -> f64 {
    let high = a.max(b);
    let low = a.min(b);
    det::sqrt(high - low) * det::sqrt(high) * det::sqrt(1.0 + low / high)
}

impl LayerKernel {
    fn new(layer: PlanarFluidLayer, coordinates: LayerCoordinates) -> Result<Self, AeroacError> {
        let media = [layer.incident_state(), layer.layer, layer.receiving_state()];
        let k = media.map(|fluid| coordinates.omega / fluid.sound_speed_m_s());
        if k.iter().any(|v| !positive(*v)) {
            return Err(nonfinite("layered acoustic wavenumbers"));
        }
        for i in 0..3 {
            for j in i + 1..3 {
                if k[i] == k[j] && media[i].sound_speed_m_s() != media[j].sound_speed_m_s() {
                    return Err(nonfinite(
                        "distinct original layer wavenumbers round together",
                    ));
                }
            }
        }
        let mut ordered = k;
        ordered.sort_by(f64::total_cmp);
        let mut critical = [0.0; 3];
        let mut critical_count = 0;
        for value in ordered {
            if critical_count == 0 || value != critical[critical_count - 1] {
                critical[critical_count] = value;
                critical_count += 1;
            }
        }
        let gaps = core::array::from_fn(|i| {
            core::array::from_fn(|j| {
                if j < critical_count {
                    wave_gap(k[i], critical[j])
                } else {
                    0.0
                }
            })
        });
        if gaps.iter().flatten().any(|v| !v.is_finite()) {
            return Err(nonfinite("layer critical-wavenumber gaps"));
        }
        let rho_log = media.map(|fluid| det::ln(fluid.density_kg_m3()));
        let logs = [
            -rho_log[0],
            -rho_log[2],
            -rho_log[1],
            rho_log[1] - rho_log[0] - rho_log[2],
        ];
        let log_density_scale = logs.into_iter().fold(f64::NEG_INFINITY, f64::max);
        let coefficients = logs.map(|v| det::exp(v - log_density_scale));
        if coefficients.iter().any(|v| !positive(*v)) {
            return Err(nonfinite(
                "layer density contrast exceeds representable coefficient ratios",
            ));
        }
        Ok(Self {
            layer,
            coordinates,
            k,
            critical,
            critical_count,
            gaps,
            coefficients,
            log_density_scale,
            common_speed: critical_count == 1,
            tail_extent: 0.0,
        })
    }

    fn point(&self, region: usize, x: f64) -> LayerPoint {
        if region == self.critical_count {
            let t = self.tail_extent * x;
            return LayerPoint {
                q: det::hypot(self.critical[self.critical_count - 1], t),
                zeta: core::array::from_fn(|i| {
                    C64::new(0.0, det::hypot(self.gaps[i][self.critical_count - 1], t))
                }),
                measure: [t, self.tail_extent],
            };
        }
        let (sin, cos) = det::sin_cos(FRAC_PI_2 * x);
        if region == 0 {
            let q = self.critical[0] * sin;
            let vertical = self.critical[0] * cos;
            return LayerPoint {
                q,
                zeta: core::array::from_fn(|i| C64::from_re(det::hypot(self.gaps[i][0], vertical))),
                measure: [q, vertical * FRAC_PI_2],
            };
        }
        let low = self.critical[region - 1];
        let high = self.critical[region];
        let gap = wave_gap(low, high);
        let imaginary = gap * sin;
        let real = gap * cos;
        LayerPoint {
            q: det::hypot(low, imaginary),
            zeta: core::array::from_fn(|i| {
                if self.k[i] <= low {
                    C64::new(0.0, det::hypot(self.gaps[i][region - 1], imaginary))
                } else {
                    C64::from_re(det::hypot(self.gaps[i][region], real))
                }
            }),
            measure: [imaginary, real * FRAC_PI_2],
        }
    }

    // Return exp(-Im(delta))*cos(delta), exp(-Im(delta))*h*sinc(delta).
    // Both are real for this lossless layer. No growing hyperbolic quantity is formed.
    fn layer_functions(&self, vertical: C64) -> (f64, f64) {
        let h = self.layer.thickness_m;
        if vertical.im > 0.0 {
            let eta = vertical.im * h;
            let decay = det::exp(-2.0 * eta);
            let ratio = if eta == 0.0 {
                1.0
            } else if eta > 350.0 {
                0.5 / eta
            } else {
                -det::expm1(-2.0 * eta) / (2.0 * eta)
            };
            (0.5 * (1.0 + decay), h * ratio)
        } else {
            let delta = Dd::from_f64(vertical.re) * Dd::from_f64(h);
            let (sin, cos) = det::sin_cos_dd(delta);
            let phase = delta.to_f64();
            let sinc = if phase.abs() < 1e-3 {
                let x2 = phase * phase;
                1.0 - x2 * (1.0 / 6.0 - x2 * (1.0 / 120.0 - x2 / 5040.0))
            } else {
                sin.to_f64() / phase
            };
            (cos.to_f64(), h * sinc)
        }
    }

    // Common-cutoff cancellation removes L's original common zeta analytically.
    // Otherwise the supplied W is an arbitrary common positive wavenumber scale.
    fn denominator(&self, point: LayerPoint, wave_scale: f64) -> C64 {
        let [a1, a3, b2, b13] = self.coefficients;
        let (cosine, sine_over_zeta) = self.layer_functions(point.zeta[1]);
        if self.common_speed {
            C64::from_re((a1 + a3) * cosine)
                + minus_i(point.zeta[1]).scale((b2 + b13) * sine_over_zeta)
        } else {
            let z = point
                .zeta
                .map(|v| C64::new(v.re / wave_scale, v.im / wave_scale));
            (z[0].scale(a1) + z[2].scale(a3)).scale(cosine)
                + minus_i(z[1] * z[1]).scale(b2 * (wave_scale * sine_over_zeta))
                + minus_i(z[0] * z[2]).scale(b13 * (wave_scale * sine_over_zeta))
        }
    }

    fn sample_weight(&self, region: usize, x: f64) -> Result<LayerSpectralWeight, AeroacError> {
        let point = self.point(region, x);
        let wave_scale = point.q.max(self.critical[self.critical_count - 1]);
        let mut denominator = self.denominator(point, wave_scale);
        let divisor = if self.common_speed {
            let magnitude = point.zeta[0].abs();
            if !positive(magnitude) {
                return Err(nonfinite("layer common-cutoff factor"));
            }
            denominator =
                denominator * C64::new(point.zeta[0].re / magnitude, point.zeta[0].im / magnitude);
            magnitude
        } else {
            wave_scale
        };
        let scale = denominator.re.abs().max(denominator.im.abs());
        if !positive(scale) {
            return Err(nonfinite("layer Green denominator"));
        }
        let normalized = C64::new(denominator.re / scale, denominator.im / scale);
        let magnitude = normalized.abs();
        let phase = self.coordinates.a * Dd::from_f64(point.zeta[0].re)
            + self.coordinates.b * Dd::from_f64(point.zeta[2].re);
        let attenuation = (self.coordinates.a * Dd::from_f64(point.zeta[0].im)
            + self.coordinates.b * Dd::from_f64(point.zeta[2].im)
            + Dd::from_f64(self.layer.thickness_m) * Dd::from_f64(point.zeta[1].im))
        .to_f64();
        if !phase.hi.is_finite()
            || !phase.lo.is_finite()
            || !attenuation.is_finite()
            || !positive(divisor)
            || !(point.q * self.coordinates.range).is_finite()
        {
            return Err(nonfinite("layer spectral phase or derivative scale"));
        }
        let (sin, cos) = det::sin_cos_dd(phase);
        let unit = C64::new(cos.to_f64(), sin.to_f64()) * normalized.conj().scale(1.0 / magnitude);
        let log_magnitude = det::ln(point.measure[0]) + det::ln(point.measure[1])
            - det::ln(divisor)
            - det::ln(scale)
            - det::ln(magnitude)
            - self.log_density_scale
            - attenuation;
        if log_magnitude.is_nan() || !finite(unit) {
            return Err(nonfinite("layer spectral scale"));
        }
        Ok(LayerSpectralWeight {
            q: point.q,
            zeta: point.zeta,
            unit,
            log_magnitude,
        })
    }

    /// Shared degree-0/1/2 absolute evanescent envelopes. Scales are logarithms
    /// of physical field prefactors BEFORE multiplying E/L by q or zeta.
    fn prepare_tail(
        &mut self,
        limits: QuadratureLimits,
        degrees: [u32; 2],
        log_scales: [f64; 2],
        channels: [usize; 2],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<[f64; 2], AeroacError> {
        let h = (self.coordinates.a + self.coordinates.b + Dd::from_f64(self.layer.thickness_m))
            .to_f64();
        let hk = h * self.critical[self.critical_count - 1];
        if !hk.is_finite() || degrees.iter().any(|d| *d > 2) || channels.contains(&0) {
            return Err(nonfinite("layer tail clearance or derivative degree"));
        }
        let rho1 = self.layer.incident_state().density_kg_m3();
        let rho3 = self.layer.receiving_state().density_kg_m3();
        let la = -det::ln(rho1);
        let lb = -det::ln(rho3);
        let largest = la.max(lb);
        let log_density =
            det::ln(2.0) - largest - det::ln(det::exp(la - largest) + det::exp(lb - largest));
        let goals = [limits.absolute_pressure, limits.absolute_velocity].map(det::ln);
        let goals = core::array::from_fn::<_, 2, _>(|i| {
            goals[i] - det::ln(4.0) - 0.5 * det::ln(channels[i] as f64)
        });
        let logs = |s: f64| {
            core::array::from_fn::<_, 2, _>(|i| {
                let shifted = s + hk + 1.0;
                let polynomial = match degrees[i] {
                    0 => 0.0,
                    1 => det::ln(shifted),
                    _ => 2.0 * det::ln(shifted) + det::ln(1.0 + (1.0 / shifted) * (1.0 / shifted)),
                };
                log_scales[i] + log_density - f64::from(degrees[i] + 1) * det::ln(h) - s
                    + polynomial
            })
        };
        let mut low = 0.0;
        let mut high = 1.0;
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
            return Err(invalid("layer evanescent tail cannot meet tolerance"));
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
        self.tail_extent = high / h;
        if !positive(self.tail_extent) {
            return Err(nonfinite("layer evanescent cutoff"));
        }
        Ok(logs(high).map(det::exp))
    }

    fn resolved(&self, region: usize, lower: f64, upper: f64) -> bool {
        let lo = self.point(region, lower);
        let hi = self.point(region, upper);
        let mid = self.point(region, lower + (upper - lower) * 0.5);
        let lengths = [
            self.coordinates.a.to_f64(),
            self.layer.thickness_m,
            self.coordinates.b.to_f64(),
        ];
        let mut phase = (hi.q - lo.q).abs() * self.coordinates.range;
        let mut decay = 0.0;
        for i in 0..3 {
            phase += lengths[i] * (hi.zeta[i].re - lo.zeta[i].re).abs();
            decay += lengths[i] * (hi.zeta[i].im - lo.zeta[i].im).abs();
        }
        if phase > FRAC_PI_2 || decay > 2.0 {
            return false;
        }
        let wave_scale = hi.q.max(self.critical[self.critical_count - 1]);
        let midpoint = self.denominator(mid, wave_scale).abs();
        if !positive(midpoint) {
            return false;
        }
        let [a1, a3, b2, b13] = self.coefficients;
        let delta = self.layer.thickness_m * (hi.zeta[1] - lo.zeta[1]).abs();
        let (cl, sl) = self.layer_functions(lo.zeta[1]);
        let (ch, sh) = self.layer_functions(hi.zeta[1]);
        if self.common_speed {
            let variation = if lo.zeta[1].im == 0.0 && hi.zeta[1].im == 0.0 {
                (a1 + a3 + b2 + b13) * delta.min(2.0)
            } else {
                let sinh_l = (lo.zeta[1].im * sl).abs();
                let sinh_h = (hi.zeta[1].im * sh).abs();
                (a1 + a3) * (ch - cl).abs() + (b2 + b13) * (sinh_h - sinh_l).abs()
            };
            return variation.is_finite() && variation <= 0.25 * midpoint;
        }
        let zlo = lo
            .zeta
            .map(|v| C64::new(v.re / wave_scale, v.im / wave_scale));
        let zhi = hi
            .zeta
            .map(|v| C64::new(v.re / wave_scale, v.im / wave_scale));
        let maximum = core::array::from_fn::<_, 3, _>(|i| zlo[i].abs().max(zhi[i].abs()));
        let variation = core::array::from_fn::<_, 3, _>(|i| (zhi[i] - zlo[i]).abs());
        let a_max = a1 * maximum[0] + a3 * maximum[2];
        let a_delta = a1 * variation[0] + a3 * variation[2];
        let b_max = b2 * maximum[1] * maximum[1] + b13 * maximum[0] * maximum[2];
        let b_delta = b2 * ((zhi[1] - zlo[1]) * (zhi[1] + zlo[1])).abs()
            + b13 * (maximum[0] * variation[2] + maximum[2] * variation[0]);
        let (c_max, c_delta, s_max, s_delta) = if lo.zeta[1].im == 0.0 && hi.zeta[1].im == 0.0 {
            let minimum_phase = (lo.zeta[1].re.min(hi.zeta[1].re) * self.layer.thickness_m).abs();
            let envelope = if minimum_phase > 1.0 {
                1.0 / minimum_phase
            } else {
                1.0
            };
            let derivative = if minimum_phase > 1.0 {
                (1.0 / minimum_phase + 1.0 / (minimum_phase * minimum_phase)).min(0.5)
            } else {
                0.5
            };
            (
                1.0,
                delta.min(2.0),
                wave_scale * self.layer.thickness_m * envelope,
                wave_scale * self.layer.thickness_m * (derivative * delta).min(2.0),
            )
        } else {
            (
                cl.abs().max(ch.abs()),
                (ch - cl).abs(),
                wave_scale * sl.abs().max(sh.abs()),
                wave_scale * (sh - sl).abs(),
            )
        };
        // Triangle inequality on A*cos - i B*h*sinc, using monotone zeta
        // endpoints and |sinc'|<=1/2. A peak cannot hide between sampled nodes
        // when the denominator varies by at most a quarter of its midpoint.
        let bound = a_delta * c_max + a_max * c_delta + b_delta * s_max + b_max * s_delta;
        bound.is_finite() && bound <= 0.25 * midpoint
    }
}

impl Spectrum<3> for LayerKernel {
    fn regions(&self) -> usize {
        self.critical_count + 1
    }
    fn extent(&self) -> (f64, f64) {
        let cutoff = det::hypot(self.critical[self.critical_count - 1], self.tail_extent);
        (cutoff, cutoff * self.coordinates.range)
    }
    fn resolved(&self, region: usize, lower: f64, upper: f64) -> bool {
        LayerKernel::resolved(self, region, lower, upper)
    }
    fn sample(&self, region: usize, x: f64) -> Result<[C64; 3], AeroacError> {
        let weight = self.sample_weight(region, x)?;
        let argument = weight.q * self.coordinates.range;
        let zero = j0(argument);
        let one = j1(argument);
        let log_velocity = -det::ln(TAU) - det::ln(self.layer.receiving_state().density_kg_m3());
        let vertical = weight.zeta[2].abs();
        let direction = if vertical == 0.0 {
            C64::ZERO
        } else {
            C64::new(weight.zeta[2].re / vertical, weight.zeta[2].im / vertical)
        };
        let values = [
            weight.scaled(
                C64::ONE,
                det::ln(self.coordinates.omega) - det::ln(TAU),
                zero,
            ),
            weight.scaled(C64::new(0.0, 1.0), log_velocity + det::ln(weight.q), one),
            weight.scaled(direction, log_velocity + det::ln(vertical), zero),
        ];
        if values.iter().any(|v| !finite(*v)) {
            return Err(nonfinite("layered acoustic spectral field"));
        }
        Ok(values)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use core::f64::consts::PI;

    fn kernel(density: f64) -> LayerKernel {
        let incident = FluidAcousticState::try_new(1.2, 300.0).unwrap();
        let internal = FluidAcousticState::try_new(density, 600.0).unwrap();
        let layer = PlanarFluidLayer::new(
            incident,
            internal,
            incident,
            [0.0; 3],
            [0.0, 0.0, 1.0],
            0.31,
        )
        .unwrap();
        let coordinates =
            LayerCoordinates::new(layer, 1500.0, [0.0, 0.0, -0.12], [0.09, 0.0, 0.39]).unwrap();
        LayerKernel::new(layer, coordinates).unwrap()
    }

    #[test]
    fn internal_grazing_and_large_evanescent_depth_keep_regular_denominators() {
        let mut kernel = kernel(2.4);
        let external = wave_gap(kernel.k[0], kernel.k[1]);
        let point = LayerPoint {
            q: kernel.k[1],
            zeta: [C64::from_re(external), C64::ZERO, C64::from_re(external)],
            measure: [0.0; 2],
        };
        let scale = kernel.k[0];
        let actual = kernel
            .denominator(point, scale)
            .scale(scale * det::exp(kernel.log_density_scale));
        let admittance = external / kernel.layer.incident_state().density_kg_m3();
        let expected = C64::new(
            2.0 * admittance,
            -kernel.layer.thickness_m
                * kernel.layer.layer_state().density_kg_m3()
                * admittance
                * admittance,
        );
        assert!((actual - expected).abs() <= 2e-14 * expected.abs());
        let alpha = 1e6;
        let (scaled_cos, scaled_sine_over_zeta) = kernel.layer_functions(C64::new(0.0, alpha));
        assert_eq!(scaled_cos, 0.5);
        assert!((scaled_sine_over_zeta * alpha - 0.5).abs() < 1e-15);
        kernel.tail_extent = alpha;
        let point = kernel.point(kernel.critical_count, 1.0);
        let denominator = kernel.denominator(point, point.q);
        assert!(finite(denominator) && denominator.abs() > 0.0);
        assert_eq!(denominator.re, 0.0);
    }

    #[test]
    fn denominator_resolution_catches_resonances_that_phase_span_cannot() {
        let kernel = kernel(1000.0);
        let centre = (PI / (kernel.k[1] * kernel.layer.thickness_m)).acos() / FRAC_PI_2;
        let lower = centre - 0.01;
        let upper = centre + 0.01;
        let lo = kernel.point(0, lower);
        let hi = kernel.point(0, upper);
        let phase = (hi.q - lo.q).abs() * kernel.coordinates.range
            + (hi.zeta[0].re - lo.zeta[0].re).abs() * kernel.coordinates.a.to_f64()
            + (hi.zeta[1].re - lo.zeta[1].re).abs() * kernel.layer.thickness_m
            + (hi.zeta[2].re - lo.zeta[2].re).abs() * kernel.coordinates.b.to_f64();
        assert!(phase < FRAC_PI_2);
        assert!(!kernel.resolved(0, lower, upper));
        assert!(kernel.resolved(0, centre - 1e-9, centre + 1e-9));
    }
}
