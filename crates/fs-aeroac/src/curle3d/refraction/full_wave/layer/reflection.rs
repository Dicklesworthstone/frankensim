//! Complete incident-side return from one finite original fluid layer.
//!
//! With the parent's admittances, C=cos(zeta2*h), S=h*sinc(zeta2*h),
//! pressure and normal-velocity continuity give R=N/L, where
//! N=(Y1-Y3)C+i(zeta2²/rho2-rho2*Y1*Y3)S and L is the original
//! transmission denominator. This retains all internal returns coherently.
//! The returning source spectrum is q*dq*R*exp(i*zeta1*(a+b))/zeta1.
//! No finite image truncation, artificial attenuation or source feedback is
//! introduced. The parent's slower-than-both guided-pole refusal still applies.
//!
//! Branch maps, bounded quadrature, density normalization, scaled C/S and
//! resonance resolution are shared with transmission. The incident grazing
//! quotient is cancelled before evaluation. Numerical tolerances control the
//! reflected transfer; a total-field interference null has no relative-error
//! guarantee. Above all three critical circles, |R|<=1 supplies the original
//! reflection tail envelope with the exterior return clearance a+b.

use super::super::reflection::{
    PreparedReflectedFullWaveMonopole, ReflectionCoordinates, ReflectionSpectralWeight, apply,
    field, line, prepare_reflected_tail, singular_total,
};
use super::super::{
    ComplexSum, FullWaveMonopoleLimits, QuadratureLimits, Spectrum, admit, finite,
    integrate_spectrum, j0, j1, minus_i, nonfinite, poll, positive,
};
use super::{
    FullWaveLayerMonopoleDiagnostics, LayerCoordinates, LayerKernel, LayerPoint, PlanarFluidLayer,
};
use crate::AeroacError;
use crate::curle3d::{CoherentCompactField, CompactSourcePhasor, PreparedCompactSourceObserver};
use crate::observation::HarmonicFieldLine;
use core::f64::consts::{FRAC_PI_2, PI};
use fs_math::{c64::C64, dd::Dd, det};

/// Applicability of the returning field, excluding the original direct field.
pub const FULL_WAVE_LAYER_REFLECTED_MONOPOLE_SCOPE: &str = "Prescribed stationary compact \
3D volume source; reflected spherical-wave field from one infinite planar homogeneous \
lossless fluid layer between original stationary lossless fluids, including coherent \
internal returns, reactive fields and evanescence. Slower-than-both guided layers refuse. \
Numerical quadrature estimates and analytic tail bounds are not interval certificates; \
no absorption, moving fluids, finite edges, source feedback or measured-SPL qualification.";

/// The source-side direct and returned amplitudes interfere before intensity.
pub const FULL_WAVE_LAYER_TOTAL_MONOPOLE_SCOPE: &str = "Prescribed stationary compact \
3D volume source; coherent direct plus finite-fluid-layer reflected spherical-wave field \
in the original incident fluid. Reflected-transfer error estimates do not guarantee \
relative accuracy at interference nulls. Slower-than-both guided layers and singular \
direct self-fields refuse; no absorption, moving fluids, finite edges, source feedback, \
interval certificate or measured-SPL qualification.";

/// Separate cached reflected and coherent total responses on the source clock.
/// Applying the original complex flow performs fixed work with no allocation,
/// quadrature or transcendental evaluation. Only the reflected field is defined
/// at a coincident source/observer; it does not feed back onto the source.
#[derive(Clone, Debug)]
pub struct PreparedLayeredReflectedFullWaveMonopole {
    layer: PlanarFluidLayer,
    frequency_hz: f64,
    source_m: [f64; 3],
    observer_m: [f64; 3],
    reflected: [C64; 4],
    total: Option<[C64; 4]>,
    diagnostics: FullWaveLayerMonopoleDiagnostics,
}

impl PreparedLayeredReflectedFullWaveMonopole {
    /// Source and observer must be strictly before the first original face.
    /// Limits apply to the reflected transfer; the direct field is analytic.
    /// A zero-thickness layer, or one identical to its receiving half-space,
    /// uses the original single-interface reflection owner exactly.
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
        let coordinates =
            ReflectionCoordinates::new(layer.exterior, frequency_hz, source_m, observer_m)?;
        let (reflected, diagnostics) = if reduces_to_first_interface(layer) {
            let original = PreparedReflectedFullWaveMonopole::new(
                layer.exterior,
                frequency_hz,
                source_m,
                observer_m,
                limits,
                cancelled,
            )?;
            let response = original.evaluate_reflected_harmonic_line(C64::ONE)?;
            (
                [
                    response.pressure_pa,
                    response.particle_velocity_m_s[0],
                    response.particle_velocity_m_s[1],
                    response.particle_velocity_m_s[2],
                ],
                (*original.diagnostics()).into(),
            )
        } else {
            let mut kernel = LayerReflectionKernel::new(layer, coordinates)?;
            let log_velocity = -det::ln(4.0 * PI);
            let bound = kernel.prepare_tail(
                limits.into(),
                [0, 1],
                [
                    det::ln(coordinates.omega)
                        + det::ln(layer.incident_state().density_kg_m3())
                        + log_velocity,
                    log_velocity,
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
            let n = layer.receiving_normal();
            let t = coordinates.tangent;
            let transfer = [
                pressure,
                radial.scale(t[0]) + normal.scale(n[0]),
                radial.scale(t[1]) + normal.scale(n[1]),
                radial.scale(t[2]) + normal.scale(n[2]),
            ];
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
        if reflected.iter().any(|v| !finite(*v)) {
            return Err(nonfinite("layer reflected monopole transfer"));
        }
        let total = if source_m == observer_m {
            None
        } else {
            let direct = PreparedCompactSourceObserver::new(
                layer.incident_state(),
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
            Some(add_direct(direct, reflected)?)
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
    /// Estimates and analytic omitted-tail bounds for the reflected transfer.
    pub const fn diagnostics(&self) -> &FullWaveLayerMonopoleDiagnostics {
        &self.diagnostics
    }
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
            FULL_WAVE_LAYER_REFLECTED_MONOPOLE_SCOPE,
        )
    }
    pub fn evaluate_total(
        &self,
        volume_flow_m3_s: C64,
    ) -> Result<CoherentCompactField, AeroacError> {
        field(
            apply(self.total.ok_or_else(singular_total)?, volume_flow_m3_s)?,
            FULL_WAVE_LAYER_TOTAL_MONOPOLE_SCOPE,
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

fn reduces_to_first_interface(layer: PlanarFluidLayer) -> bool {
    layer.thickness_m == 0.0 || layer.layer_state() == layer.receiving_state()
}

fn add_direct(direct: [C64; 4], reflected: [C64; 4]) -> Result<[C64; 4], AeroacError> {
    let mut total = [C64::ZERO; 4];
    for row in 0..4 {
        let mut sum = ComplexSum::default();
        sum.add(direct[row])?;
        sum.add(reflected[row])?;
        total[row] = sum.value();
    }
    if total.iter().any(|v| !finite(*v)) {
        return Err(nonfinite("total layer-reflected field"));
    }
    Ok(total)
}

struct LayerReflectionKernel {
    common: LayerKernel,
}

impl LayerReflectionKernel {
    fn new(
        layer: PlanarFluidLayer,
        coordinates: ReflectionCoordinates,
    ) -> Result<Self, AeroacError> {
        // Reuse the transmission kernel's phase-resolution guard with the
        // actual exterior return H*zeta1; zeta3 carries no propagation phase.
        Ok(Self {
            common: LayerKernel::new(
                layer,
                LayerCoordinates {
                    omega: coordinates.omega,
                    a: coordinates.height,
                    b: Dd::from_f64(0.0),
                    range: coordinates.range,
                    tangent: coordinates.tangent,
                },
            )?,
        })
    }

    fn numerator(&self, point: LayerPoint, wave_scale: f64) -> C64 {
        let [a1, a3, b2, b13] = self.common.coefficients;
        let (cosine, sine_over_zeta) = self.common.layer_functions(point.zeta[1]);
        if self.common.common_speed {
            C64::from_re((a1 - a3) * cosine)
                - minus_i(point.zeta[1]).scale((b2 - b13) * sine_over_zeta)
        } else {
            let z = point
                .zeta
                .map(|v| C64::new(v.re / wave_scale, v.im / wave_scale));
            (z[0].scale(a1) - z[2].scale(a3)).scale(cosine)
                - minus_i(z[1] * z[1]).scale(b2 * (wave_scale * sine_over_zeta))
                + minus_i(z[0] * z[2]).scale(b13 * (wave_scale * sine_over_zeta))
        }
    }

    fn sample_weight(
        &self,
        region: usize,
        x: f64,
    ) -> Result<ReflectionSpectralWeight, AeroacError> {
        let common = &self.common;
        let point = common.point(region, x);
        let kmax = common.critical[common.critical_count - 1];
        let wave_scale = point.q.max(kmax);
        let numerator = self.numerator(point, wave_scale);
        let denominator = common.denominator(point, wave_scale);
        let ns = numerator.re.abs().max(numerator.im.abs());
        let ds = denominator.re.abs().max(denominator.im.abs());
        if !ns.is_finite() || !positive(ds) {
            return Err(nonfinite("layer reflection numerator or denominator"));
        }
        if ns == 0.0 {
            return Ok(ReflectionSpectralWeight {
                q: point.q,
                zeta: point.zeta[0],
                unit: C64::ZERO,
                log_magnitude: 0.0,
            });
        }
        let numerator = C64::new(numerator.re / ns, numerator.im / ns);
        let denominator = C64::new(denominator.re / ds, denominator.im / ds);
        let na = numerator.abs();
        let da = denominator.abs();
        let reflection_unit = numerator.scale(1.0 / na) * denominator.conj().scale(1.0 / da);
        let reflection_log = det::ln(ns) - det::ln(ds) + det::ln(na) - det::ln(da);
        let phase = common.coordinates.a * Dd::from_f64(point.zeta[0].re);
        let attenuation = (common.coordinates.a * Dd::from_f64(point.zeta[0].im)).to_f64();
        if !phase.hi.is_finite()
            || !phase.lo.is_finite()
            || !attenuation.is_finite()
            || !(point.q * common.coordinates.range).is_finite()
        {
            return Err(nonfinite("layer reflected spectral phase"));
        }
        let negative_i = C64::new(0.0, -1.0);
        // Cancel q*dq/zeta1 analytically wherever the incident root vanishes.
        let (direction, log_measure) = if region == common.critical_count && common.k[0] == kmax {
            (negative_i, det::ln(common.tail_extent))
        } else if region == 0 && common.k[0] == common.critical[0] {
            (C64::ONE, det::ln(point.q) + det::ln(FRAC_PI_2))
        } else if region > 0
            && region < common.critical_count
            && common.k[0] == common.critical[region - 1]
        {
            (negative_i, det::ln(point.measure[1]))
        } else if region > 0
            && region < common.critical_count
            && common.k[0] == common.critical[region]
        {
            (C64::ONE, det::ln(point.measure[0]) + det::ln(FRAC_PI_2))
        } else {
            (
                if point.zeta[0].im > 0.0 {
                    negative_i
                } else {
                    C64::ONE
                },
                det::ln(point.measure[0]) + det::ln(point.measure[1])
                    - det::ln(point.zeta[0].abs()),
            )
        };
        let (sin, cos) = det::sin_cos_dd(phase);
        let unit = C64::new(cos.to_f64(), sin.to_f64()) * direction * reflection_unit;
        let log_magnitude = log_measure + reflection_log - attenuation;
        if log_magnitude.is_nan() || !finite(unit) {
            return Err(nonfinite("layer reflected spectral scale"));
        }
        Ok(ReflectionSpectralWeight {
            q: point.q,
            zeta: point.zeta[0],
            unit,
            log_magnitude,
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
            self.common.coordinates.a.to_f64(),
            self.common.critical[self.common.critical_count - 1],
            limits,
            degrees,
            log_scales,
            channels,
            cancelled,
        )?;
        self.common.tail_extent = extent;
        Ok(bounds)
    }
}

impl Spectrum<3> for LayerReflectionKernel {
    fn regions(&self) -> usize {
        <LayerKernel as Spectrum<3>>::regions(&self.common)
    }
    fn extent(&self) -> (f64, f64) {
        <LayerKernel as Spectrum<3>>::extent(&self.common)
    }
    fn resolved(&self, region: usize, lower: f64, upper: f64) -> bool {
        // The existing variation bound also bounds N: its terms differ from
        // L only in signs. No relative-N guard may reject a physical R=0.
        self.common.resolved(region, lower, upper)
    }
    fn sample(&self, region: usize, x: f64) -> Result<[C64; 3], AeroacError> {
        let weight = self.sample_weight(region, x)?;
        let argument = weight.q * self.common.coordinates.range;
        let zero = j0(argument);
        let log_velocity = -det::ln(4.0 * PI);
        let log_pressure = log_velocity
            + det::ln(self.common.coordinates.omega)
            + det::ln(self.common.layer.incident_state().density_kg_m3());
        let magnitude = weight.zeta.abs();
        let normal = if magnitude == 0.0 {
            C64::ZERO
        } else {
            C64::new(-weight.zeta.re / magnitude, -weight.zeta.im / magnitude)
        };
        let values = [
            weight.scaled(C64::ONE, log_pressure, zero),
            weight.scaled(
                C64::new(0.0, 1.0),
                log_velocity + det::ln(weight.q),
                j1(argument),
            ),
            weight.scaled(normal, log_velocity + det::ln(magnitude), zero),
        ];
        if values.iter().any(|v| !finite(*v)) {
            return Err(nonfinite("layer reflected monopole spectrum"));
        }
        Ok(values)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use fs_material::fluid_acoustics::FluidAcousticState;

    fn kernel(rho: [f64; 3], c: [f64; 3], hz: f64, h: f64, height: f64) -> LayerReflectionKernel {
        let fluids =
            core::array::from_fn::<_, 3, _>(|i| FluidAcousticState::try_new(rho[i], c[i]).unwrap());
        let layer = PlanarFluidLayer::new(
            fluids[0],
            fluids[1],
            fluids[2],
            [0.0; 3],
            [0.0, 0.0, 1.0],
            h,
        )
        .unwrap();
        let coordinates = ReflectionCoordinates::new(
            layer.exterior,
            hz,
            [0.0, 0.0, -height * 0.5],
            [0.03, 0.0, -height * 0.5],
        )
        .unwrap();
        let mut kernel = LayerReflectionKernel::new(layer, coordinates).unwrap();
        kernel.common.tail_extent = 400.0;
        kernel
    }

    #[test]
    fn all_incident_grazing_positions_and_large_scaled_fields_stay_finite() {
        for c in [
            [300.0, 600.0, 400.0],
            [600.0, 450.0, 300.0],
            [350.0; 3],
            [300.0, 300.0, 450.0],
            [300.0, 600.0, 300.0],
        ] {
            let kernel = kernel([1.2, 2.4, 1.8], c, 700.0, 0.018, 0.05);
            for region in 0..kernel.regions() {
                for x in [0.0, 1.0] {
                    assert!(kernel.sample(region, x).unwrap().iter().all(|v| finite(*v)));
                }
            }
        }
        let large = kernel([1e250, 2e250, 1.5e250], [1.0, 2.0, 1.5], 1.0, 0.2, 2.0);
        assert_eq!(det::exp(-800.0), 0.0);
        let p = large.sample(large.common.critical_count, 1.0).unwrap()[0];
        assert!(finite(p) && p.abs() > 0.0);
    }

    #[test]
    fn fabry_perot_reflection_zero_does_not_require_relative_numerator_resolution() {
        let hz = 700.0;
        let c = 350.0;
        let kernel = kernel([1.2, 5.0, 1.2], [c; 3], hz, c / (2.0 * hz), 0.05);
        let point = kernel.common.point(0, 0.0);
        let ratio = kernel.numerator(point, kernel.common.k[0])
            / kernel.common.denominator(point, kernel.common.k[0]);
        assert!(ratio.abs() < 1e-13);
        assert!(kernel.resolved(0, 0.0, 1e-5));
        assert!(kernel.sample(0, 0.0).unwrap().iter().all(|v| finite(*v)));
    }
}
