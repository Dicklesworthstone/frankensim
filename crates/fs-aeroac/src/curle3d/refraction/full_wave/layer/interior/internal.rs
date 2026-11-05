//! Original volume sources and receivers strictly inside a finite fluid layer.
//!
//! The direct spherical field is analytic. Only the reflected field is
//! integrated, so equal source/receiver depths do not create a nondecaying
//! direct Sommerfeld tail. With source/receiver depths s,z, Yj=zeta_j/rho_j,
//! Cj(d)=cos(zeta2*d)-i*rho2*Yj*sin(zeta2*d)/zeta2, and the original layer
//! denominator L, its reflected spectral bracket is
//! H=[(Y2-Y1)exp(i*zeta2*s)C3(h-z)
//!   +(Y2-Y3)exp(i*zeta2*(h-s))C1(z)]/L.
//! This regular form retains all internal returns without subtracting a large
//! direct spectrum or forming the grazing 0/0 of reflection-coefficient ratios.
//! Differentiating C3(h-z) and C1(z) gives +i*rho2*V3 and -i*rho2*V1;
//! receiving velocity is the derivative of this same pressure field.

use super::super::{LayerCoordinates, LayerKernel, LayerSpectralWeight, wave_gap};
use super::ScaledNumerator;
use crate::AeroacError;
use crate::curle3d::refraction::full_wave::layer::{
    FullWaveLayerMonopoleDiagnostics, PlanarFluidLayer,
};
use crate::curle3d::refraction::full_wave::reflection::{
    PreparedReflectedFullWaveMonopole, apply, field, line, prepare_reflected_tail, singular_total,
};
use crate::curle3d::refraction::full_wave::{
    FullWaveMonopoleLimits, Spectrum, admit, integrate_spectrum, j0, j1, poll,
};
use crate::curle3d::refraction::{PlanarFluidInterface, norm, positive};
use crate::curle3d::statistics::{ComplexSum, finite, invalid, nonfinite};
use crate::curle3d::{CoherentCompactField, CompactSourcePhasor, PreparedCompactSourceObserver};
use crate::observation::HarmonicFieldLine;
use core::f64::consts::{FRAC_PI_2, PI, TAU};
use fs_math::{c64::C64, dd::Dd, det};

pub mod array;
mod force;
pub use force::{
    FULL_WAVE_LAYER_INTERNAL_REFLECTED_FORCE_SCOPE, FULL_WAVE_LAYER_INTERNAL_TOTAL_FORCE_SCOPE,
    PreparedLayerInternalFullWaveForce,
};

pub const FULL_WAVE_LAYER_INTERNAL_REFLECTED_SCOPE: &str = "Prescribed stationary \
compact volume source and receiver strictly inside one infinite planar homogeneous \
lossless fluid layer; complete coherent reflected spherical field, including both \
return directions, critical, reactive and evanescent components. The finite returned \
self-field does not imply source loading. Slower-than-both guided layers refuse. \
Numerical estimates and analytic tail bounds are not interval certificates; no \
absorption, moving fluids, finite edges or measured-SPL qualification.";

pub const FULL_WAVE_LAYER_INTERNAL_TOTAL_SCOPE: &str = "Prescribed stationary compact \
volume source and receiver strictly inside one infinite planar homogeneous lossless \
fluid layer; exact direct spherical wave plus every coherent reflected return and \
complete particle velocity. Singular total self-fields and slower-than-both guided \
layers refuse. Reflected-transfer estimates give no relative guarantee at total-field \
interference nulls; no source loading, absorption, moving fluids, finite edges, \
interval certificate or measured-SPL qualification.";

/// Cached reflected and total p/u transfers for an original source inside the layer.
/// Application has fixed work, no allocation, quadrature or extra propagation delay.
#[derive(Clone, Debug)]
pub struct PreparedLayerInternalFullWaveMonopole {
    layer: PlanarFluidLayer,
    frequency_hz: f64,
    source_m: [f64; 3],
    observer_m: [f64; 3],
    reflected: [C64; 4],
    total: Option<[C64; 4]>,
    diagnostics: FullWaveLayerMonopoleDiagnostics,
}

impl PreparedLayerInternalFullWaveMonopole {
    /// Both points lie strictly between the original faces. A coincident pair
    /// admits its finite reflected field; every total-field method then refuses.
    /// Limits control the reflected transfer. The two returning families share
    /// one evaluation budget and divide the absolute tolerances; they use zero
    /// relative tolerance so cancellation cannot relax the combined absolute
    /// target. An absent second face uses the original interface owner and its
    /// original limits, including exact homogeneous and equal-speed reductions.
    /// No partial observer is published on budget exhaustion or cancellation.
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
        let geometry = InternalCoordinates::new(layer, frequency_hz, source_m, observer_m)?;
        let (reflected, diagnostics) = if layer.layer_state() == layer.receiving_state() {
            let interface = PlanarFluidInterface::new(
                layer.layer_state(),
                layer.incident_state(),
                layer.first_face_point_m(),
                layer.receiving_normal().map(|v| -v),
            )?;
            let original = PreparedReflectedFullWaveMonopole::new(
                interface,
                frequency_hz,
                source_m,
                observer_m,
                limits,
                cancelled,
            )?;
            let value = original.evaluate_reflected_harmonic_line(C64::ONE)?;
            (
                [
                    value.pressure_pa,
                    value.particle_velocity_m_s[0],
                    value.particle_velocity_m_s[1],
                    value.particle_velocity_m_s[2],
                ],
                (*original.diagnostics()).into(),
            )
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
            Some(total)
        };
        if reflected
            .iter()
            .chain(total.iter().flatten())
            .any(|v| !finite(*v))
        {
            return Err(nonfinite("internal layer monopole transfer"));
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
    /// Summed absolute estimates/bounds for the reflected transfer per unit Q.
    /// For two families the cutoff and Bessel argument are their maxima; work
    /// and retained panels are summed, while maximum subdivision depth is not.
    pub const fn diagnostics(&self) -> &FullWaveLayerMonopoleDiagnostics {
        &self.diagnostics
    }
    pub const fn total_field_available(&self) -> bool {
        self.total.is_some()
    }

    /// Original complex peak Q [m³/s], on the exp(-i omega t) source clock.
    pub fn evaluate_reflected(&self, q: C64) -> Result<CoherentCompactField, AeroacError> {
        field(
            apply(self.reflected, q)?,
            FULL_WAVE_LAYER_INTERNAL_REFLECTED_SCOPE,
        )
    }
    pub fn evaluate_total(&self, q: C64) -> Result<CoherentCompactField, AeroacError> {
        field(
            apply(self.total.ok_or_else(singular_total)?, q)?,
            FULL_WAVE_LAYER_INTERNAL_TOTAL_SCOPE,
        )
    }
    pub fn evaluate_reflected_harmonic_line(
        &self,
        q: C64,
    ) -> Result<HarmonicFieldLine<3>, AeroacError> {
        Ok(line(self.frequency_hz, apply(self.reflected, q)?))
    }
    pub fn evaluate_total_harmonic_line(
        &self,
        q: C64,
    ) -> Result<HarmonicFieldLine<3>, AeroacError> {
        Ok(line(
            self.frequency_hz,
            apply(self.total.ok_or_else(singular_total)?, q)?,
        ))
    }
}

struct InternalCoordinates {
    spectral: LayerCoordinates,
    source: [Dd; 2],
    observer: [Dd; 2],
}

impl InternalCoordinates {
    fn new(
        layer: PlanarFluidLayer,
        hz: f64,
        source: [f64; 3],
        observer: [f64; 3],
    ) -> Result<Self, AeroacError> {
        let omega = TAU * hz;
        if !positive(omega) {
            return Err(invalid("internal full-wave frequency"));
        }
        if source.iter().chain(&observer).any(|v| !v.is_finite()) {
            return Err(nonfinite("internal full-wave positions"));
        }
        let h = Dd::from_f64(layer.thickness_m());
        let s = layer.exterior.distance(source);
        let z = layer.exterior.distance(observer);
        let source_depths = [s, h - s];
        let observer_depths = [z, h - z];
        if source_depths
            .iter()
            .chain(&observer_depths)
            .any(|v| !positive(v.to_f64()))
        {
            return Err(invalid(
                "source and observer must lie strictly inside the original layer",
            ));
        }
        let lateral = core::array::from_fn(|axis| {
            (Dd::from_f64(observer[axis])
                - Dd::from_f64(source[axis])
                - (z - s) * Dd::from_f64(layer.receiving_normal()[axis]))
            .to_f64()
        });
        let range = norm(lateral);
        if !range.is_finite() {
            return Err(nonfinite("internal layer lateral separation"));
        }
        Ok(Self {
            spectral: LayerCoordinates {
                omega,
                a: Dd::ZERO,
                b: Dd::ZERO,
                range,
                tangent: if range == 0.0 {
                    [0.0; 3]
                } else {
                    lateral.map(|v| v / range)
                },
            },
            source: source_depths,
            observer: observer_depths,
        })
    }
}

fn reflected_transfer(
    layer: PlanarFluidLayer,
    geometry: InternalCoordinates,
    limits: FullWaveMonopoleLimits,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<([C64; 4], FullWaveLayerMonopoleDiagnostics), AeroacError> {
    let media = [
        layer.incident_state(),
        layer.layer_state(),
        layer.receiving_state(),
    ];
    let active = [media[0] != media[1], media[2] != media[1]];
    let count = active.iter().filter(|v| **v).count();
    let mut part_limits = FullWaveMonopoleLimits {
        absolute_pressure_tolerance_pa_per_m3_s: limits.absolute_pressure_tolerance_pa_per_m3_s
            / count as f64,
        absolute_velocity_tolerance_per_m2: limits.absolute_velocity_tolerance_per_m2
            / count as f64,
        relative_tolerance: 0.0,
        ..limits
    };
    admit(part_limits)?;
    let mut common = LayerKernel::new(layer, geometry.spectral)?;
    let rho_log = media.map(|fluid| det::ln(fluid.density_kg_m3()));
    let log_velocity = -det::ln(4.0 * PI);
    let log_pressure = log_velocity + det::ln(geometry.spectral.omega) + rho_log[1];
    let mut sums = [ComplexSum::default(); 3];
    let mut diagnostics = FullWaveLayerMonopoleDiagnostics {
        evaluations: 0,
        maximum_depth_reached: 0,
        retained_panels: 0,
        transverse_wavenumber_cutoff_rad_m: 0.0,
        maximum_bessel_argument: 0.0,
        pressure_quadrature_error_pa_per_m3_s: 0.0,
        velocity_quadrature_error_per_m2: [0.0; 2],
        pressure_tail_bound_pa_per_m3_s: 0.0,
        velocity_tail_bound_per_m2: [0.0; 2],
        reduced_to_single_interface: false,
    };
    for face in 0..2 {
        if !active[face] {
            continue;
        }
        poll(cancelled)?;
        if diagnostics.evaluations == limits.maximum_evaluations {
            return Err(invalid(
                "internal layer exhausted its shared evaluation budget",
            ));
        }
        part_limits.maximum_evaluations = limits.maximum_evaluations - diagnostics.evaluations;
        let return_depth = geometry.source[face] + geometry.observer[face];
        let depth = return_depth.to_f64();
        if !positive(depth) {
            return Err(nonfinite("internal layer return clearance"));
        }
        // Above every critical circle, |r21|,|r23|<=1. Each family is
        // r_face exp(-alpha2*depth)(1 +/- r_other exp(-2alpha2*d)) /
        // (1-r21*r23 exp(-2alpha2*h)); the minus sign gives its normal derivative.
        // At t*h>=1/2, both magnitudes are <=4exp(-t*depth), t²=q²-kmax².
        // Reuse the original absolute tail integrals with that physical factor.
        let (extent, bound) = prepare_reflected_tail(
            depth,
            common.critical[common.critical_count - 1],
            part_limits.into(),
            [0, 1],
            [log_pressure + det::ln(4.0), log_velocity + det::ln(4.0)],
            [1, 2],
            cancelled,
        )?;
        common.tail_extent = extent.max(0.5 / layer.thickness_m());
        if !positive(common.tail_extent) {
            return Err(nonfinite("internal layer tail cutoff"));
        }
        // Increasing the retained cutoff preserves the conservative omitted-tail
        // bounds returned by the common tail owner.
        let tail = [
            bound[0],
            if geometry.spectral.range == 0.0 {
                0.0
            } else {
                bound[1]
            },
            bound[1],
        ];
        let kernel = FaceKernel {
            common: &common,
            face,
            source_depth: geometry.source[face],
            remaining_depth: geometry.observer[1 - face],
            return_depth,
            rho_log,
            log_pressure,
            log_velocity,
        };
        let result = integrate_spectrum::<3, 1>(&kernel, tail, part_limits.into(), cancelled)?;
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
        diagnostics.pressure_quadrature_error_pa_per_m3_s += result.error[0];
        diagnostics.pressure_tail_bound_pa_per_m3_s += tail[0];
        for axis in 0..2 {
            diagnostics.velocity_quadrature_error_per_m2[axis] += result.error[axis + 1];
            diagnostics.velocity_tail_bound_per_m2[axis] += tail[axis + 1];
        }
    }
    let [pressure, radial, normal] = sums.map(ComplexSum::value);
    let n = layer.receiving_normal();
    let t = geometry.spectral.tangent;
    Ok((
        [
            pressure,
            radial.scale(t[0]) + normal.scale(n[0]),
            radial.scale(t[1]) + normal.scale(n[1]),
            radial.scale(t[2]) + normal.scale(n[2]),
        ],
        diagnostics,
    ))
}

struct FaceKernel<'a> {
    common: &'a LayerKernel,
    face: usize,
    source_depth: Dd,
    remaining_depth: Dd,
    return_depth: Dd,
    rho_log: [f64; 3],
    log_pressure: f64,
    log_velocity: f64,
}

/// Original per-face pressure and receiving-velocity basis, shared with the
/// physical source-force derivative without another boundary solve.
struct FaceWeight {
    spectrum: LayerSpectralWeight,
    pressure: ScaledNumerator,
    velocity: ScaledNumerator,
    direction: C64,
    log_coefficient: f64,
}

fn unit(value: C64) -> C64 {
    let length = value.abs();
    if length == 0.0 {
        C64::ZERO
    } else {
        C64::new(value.re / length, value.im / length)
    }
}

impl FaceKernel<'_> {
    fn admittance_difference(&self, z: [C64; 3]) -> Result<ScaledNumerator, AeroacError> {
        let outer = 2 * self.face;
        let fluids = [
            self.common.layer.incident_state(),
            self.common.layer.layer_state(),
            self.common.layer.receiving_state(),
        ];
        let density_difference = fluids[outer].density_kg_m3() - fluids[1].density_kg_m3();
        // Y2-Yj = zeta2*(rhoj-rho2)/(rho2*rhoj) + (zeta2-zetaj)/rhoj.
        // Rationalize the second difference to preserve small original contrast
        // in a large evanescent wavenumber; never form squared wavenumbers.
        let (wave_unit, wave_log) = if self.common.k[1] == self.common.k[outer] {
            (C64::ZERO, f64::NEG_INFINITY)
        } else {
            let scale = z[1].abs().max(z[outer].abs());
            let sum = C64::new(
                z[1].re / scale + z[outer].re / scale,
                z[1].im / scale + z[outer].im / scale,
            );
            (
                unit(sum.conj()).scale((self.common.k[1] - self.common.k[outer]).signum()),
                2.0 * det::ln(wave_gap(self.common.k[1], self.common.k[outer]))
                    - det::ln(scale)
                    - det::ln(sum.abs())
                    - self.rho_log[outer],
            )
        };
        ScaledNumerator::sum(
            unit(z[1]).scale(density_difference.signum()),
            det::ln(z[1].abs()) + det::ln(density_difference.abs())
                - self.rho_log[1]
                - self.rho_log[outer],
            wave_unit,
            wave_log,
        )
    }

    fn sample_weight(&self, region: usize, x: f64) -> Result<FaceWeight, AeroacError> {
        let spectrum = self
            .common
            .sample_weight_at_depth(region, x, self.return_depth)?;
        let z = spectrum.zeta;
        let outer = 2 * (1 - self.face);
        let coefficient = self.admittance_difference(z)?;
        let (cosine, sine) = LayerKernel::layer_functions_at(z[1], self.remaining_depth);
        let lc = det::ln(cosine.abs());
        let ls = det::ln(sine.abs());
        let l2 = det::ln(z[1].abs());
        let lj = det::ln(z[outer].abs());
        let negative_i = C64::new(0.0, -1.0);
        let pressure = ScaledNumerator::sum(
            C64::from_re(cosine.signum()),
            lc,
            (negative_i * unit(z[outer])).scale(sine.signum()),
            lj + ls + self.rho_log[1] - self.rho_log[outer],
        )?;
        let velocity = ScaledNumerator::sum(
            unit(z[outer]).scale(cosine.signum()),
            lj + lc - self.rho_log[outer],
            (negative_i * unit(z[1]) * unit(z[1])).scale(sine.signum()),
            2.0 * l2 + ls - self.rho_log[1],
        )?;
        let phase = self.source_depth * Dd::from_f64(z[1].re);
        if !phase.hi.is_finite() || !phase.lo.is_finite() {
            return Err(nonfinite("internal layer source phase"));
        }
        let (sin, cos) = det::sin_cos_dd(phase);
        let direction = C64::new(cos.to_f64(), sin.to_f64()) * coefficient.unit * unit(z[1]).conj();
        let log_coefficient = coefficient.log_magnitude - l2;
        Ok(FaceWeight {
            spectrum,
            pressure,
            velocity,
            direction,
            log_coefficient,
        })
    }
}

impl Spectrum<3> for FaceKernel<'_> {
    fn regions(&self) -> usize {
        self.common.critical_count + 1
    }
    fn extent(&self) -> (f64, f64) {
        self.common.extent()
    }
    fn resolved(&self, region: usize, lower: f64, upper: f64) -> bool {
        let lo = self.common.point(region, lower);
        let hi = self.common.point(region, upper);
        let longest = (Dd::from_f64(self.common.layer.thickness_m()) + self.source_depth).to_f64();
        let phase = (hi.q - lo.q).abs() * self.common.coordinates.range
            + longest * (hi.zeta[1].re - lo.zeta[1].re).abs();
        let decay_depth = if lo.zeta[1].im * self.remaining_depth.to_f64() > 350.0 {
            self.return_depth.to_f64()
        } else {
            longest
        };
        phase <= FRAC_PI_2
            && self
                .common
                .resolved_at_depth(region, lower, upper, decay_depth)
    }
    fn sample(&self, region: usize, x: f64) -> Result<[C64; 3], AeroacError> {
        let FaceWeight {
            spectrum,
            pressure,
            velocity,
            direction,
            log_coefficient,
        } = self.sample_weight(region, x)?;
        let argument = spectrum.q * self.common.coordinates.range;
        let zero = j0(argument);
        let values = [
            spectrum.scaled(
                direction * pressure.unit,
                log_coefficient + pressure.log_magnitude + self.log_pressure,
                zero,
            ),
            spectrum.scaled(
                C64::new(0.0, 1.0) * direction * pressure.unit,
                log_coefficient + pressure.log_magnitude + self.log_velocity + det::ln(spectrum.q),
                j1(argument),
            ),
            spectrum.scaled(
                direction * velocity.unit,
                log_coefficient + velocity.log_magnitude + self.log_velocity + self.rho_log[1],
                if self.face == 0 { zero } else { -zero },
            ),
        ];
        if values.iter().any(|v| !finite(*v)) {
            return Err(nonfinite("internal layer reflected spectrum"));
        }
        Ok(values)
    }
}
