//! Spherical-wave fields inside the original finite fluid layer.
//!
//! At depth z, put d=h-z, Y3=zeta3/rho3 and retain the parent's full
//! denominator L. The interior pressure numerator is
//! C=cos(zeta2*d)-i*rho2*Y3*d*sinc(zeta2*d); normal velocity uses
//! V=Y3*cos(zeta2*d)-i*zeta2²/rho2*d*sinc(zeta2*d).
//! These are the same pressure/normal-velocity continuity equations as the
//! exterior owner, not a local plane-wave p/Z approximation. Both traveling
//! directions, reactive fields and every coherent internal return remain.

use super::super::reflection::{apply, field, line};
use super::{
    AeroacError, CoherentCompactField, Dd, FullWaveLayerMonopoleDiagnostics,
    FullWaveMonopoleLimits, HarmonicFieldLine, LayerCoordinates, LayerKernel, LayerSpectralWeight,
    PlanarFluidInterface, PlanarFluidLayer, PreparedFullWaveMonopole, Spectrum, admit, finite,
    integrate_spectrum, invalid, j0, j1, minus_i, nonfinite, norm, poll, positive,
};
use core::f64::consts::{PI, TAU};
use fs_math::{c64::C64, det};

/// Original sources and receivers both inside a finite refracting layer.
pub mod internal;
/// Original volume sources inside the layer radiating into its incident exterior.
pub mod source;

/// Coherent original exterior volume-flow and force arrays observed in the layer.
pub mod array;
mod force;
pub use force::{FULL_WAVE_LAYER_INTERIOR_FORCE_SCOPE, PreparedInteriorFullWaveForce};
mod quadrupole;
pub use quadrupole::{
    FULL_WAVE_LAYER_INTERIOR_QUADRUPOLE_SCOPE, PreparedInteriorFullWaveQuadrupole,
};

pub(super) fn reduces_to_interface(layer: PlanarFluidLayer) -> bool {
    layer.layer_state() == layer.receiving_state()
}

/// Applicability of an exterior physical volume source observed in the layer.
pub const FULL_WAVE_LAYER_INTERIOR_MONOPOLE_SCOPE: &str = "Prescribed stationary \
compact 3D volume source in the incident exterior; full spherical-wave pressure \
and velocity strictly inside one infinite planar homogeneous lossless fluid layer. \
Both internal directions, evanescence and coherent returns are retained. Trapping \
layers refuse; no absorption, moving media, finite edges, source loading, interval \
certificate or measured-SPL qualification.";

/// Cached interior pressure and all three velocity components per original Q.
/// Preparation shares the exterior owner's branches, resonance resolution and
/// bounded quadrature. Applying the original source-clock amplitude is allocation
/// free and needs four complex multiplies, with no new boundary solve.
#[derive(Clone, Debug)]
pub struct PreparedLayerInteriorFullWaveMonopole {
    layer: PlanarFluidLayer,
    frequency_hz: f64,
    source_m: [f64; 3],
    observer_m: [f64; 3],
    transfer: [C64; 4],
    diagnostics: FullWaveLayerMonopoleDiagnostics,
}

impl PreparedLayerInteriorFullWaveMonopole {
    /// Source strictly before the first face; receiver strictly inside the
    /// original layer. An identical layer/receiving fluid delegates the original
    /// single-interface owner, including its exact homogeneous spherical field.
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
        interior_coordinates(layer, frequency_hz, source_m, observer_m)?;
        let (transfer, diagnostics) = if layer.layer_state() == layer.receiving_state() {
            let interface = PlanarFluidInterface::new(
                layer.incident_state(),
                layer.layer_state(),
                layer.first_face_point_m(),
                layer.receiving_normal(),
            )?;
            let original = PreparedFullWaveMonopole::new(
                interface,
                frequency_hz,
                source_m,
                observer_m,
                limits,
                cancelled,
            )?;
            (original.transfer, (*original.diagnostics()).into())
        } else {
            let mut kernel = InteriorKernel::new(layer, frequency_hz, source_m, observer_m)?;
            let tail =
                kernel.prepare_tail(limits, layer.layer_state().density_kg_m3(), cancelled)?;
            let result = integrate_spectrum::<3, 1>(&kernel, tail, limits.into(), cancelled)?;
            let [pressure, radial, normal] = result.value;
            let tangent = kernel.common.coordinates.tangent;
            let transfer = core::array::from_fn(|i| {
                if i == 0 {
                    pressure
                } else {
                    radial.scale(tangent[i - 1]) + normal.scale(layer.receiving_normal()[i - 1])
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
            return Err(nonfinite("interior layered monopole transfer"));
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

    /// Complex peak Q [m³/s] on the original exp(-i omega t) source clock.
    pub fn evaluate(&self, volume_flow_m3_s: C64) -> Result<CoherentCompactField, AeroacError> {
        field(
            apply(self.transfer, volume_flow_m3_s)?,
            FULL_WAVE_LAYER_INTERIOR_MONOPOLE_SCOPE,
        )
    }

    /// Propagation and standing-wave phases are already present; add no delay.
    pub fn evaluate_harmonic_line(
        &self,
        volume_flow_m3_s: C64,
    ) -> Result<HarmonicFieldLine<3>, AeroacError> {
        Ok(line(
            self.frequency_hz,
            apply(self.transfer, volume_flow_m3_s)?,
        ))
    }
}

fn interior_coordinates(
    layer: PlanarFluidLayer,
    frequency_hz: f64,
    source: [f64; 3],
    observer: [f64; 3],
) -> Result<(LayerCoordinates, Dd, Dd), AeroacError> {
    let omega = TAU * frequency_hz;
    if !positive(omega) {
        return Err(invalid("interior full-wave frequency"));
    }
    if source.iter().chain(&observer).any(|v| !v.is_finite()) {
        return Err(nonfinite("interior full-wave positions"));
    }
    let a = -layer.exterior.distance(source);
    let depth = layer.exterior.distance(observer);
    let remaining = Dd::from_f64(layer.thickness_m()) - depth;
    if !positive(a.to_f64()) || !positive(depth.to_f64()) || !positive(remaining.to_f64()) {
        return Err(invalid(
            "source must be before the first face and observer strictly inside the layer",
        ));
    }
    let lateral = core::array::from_fn(|axis| {
        (Dd::from_f64(observer[axis])
            - Dd::from_f64(source[axis])
            - (a + depth) * Dd::from_f64(layer.receiving_normal()[axis]))
        .to_f64()
    });
    let range = norm(lateral);
    if !range.is_finite() {
        return Err(nonfinite("interior full-wave lateral separation"));
    }
    let tangent = if range == 0.0 {
        [0.0; 3]
    } else {
        lateral.map(|v| v / range)
    };
    Ok((
        LayerCoordinates {
            omega,
            a,
            b: Dd::from_f64(0.0),
            range,
            tangent,
        },
        depth,
        remaining,
    ))
}

/// Keep coefficient magnitudes logarithmic until combined with physical decay.
/// This avoids forming rho2/rho3 or zeta2² at large original contrasts.
struct ScaledNumerator {
    unit: C64,
    log_magnitude: f64,
}

impl ScaledNumerator {
    fn sum(a: C64, log_a: f64, b: C64, log_b: f64) -> Result<Self, AeroacError> {
        let scale = log_a.max(log_b);
        if scale == f64::NEG_INFINITY {
            return Ok(Self {
                unit: C64::ZERO,
                log_magnitude: scale,
            });
        }
        if !scale.is_finite() {
            return Err(nonfinite("interior acoustic numerator scale"));
        }
        let value = a.scale(det::exp(log_a - scale)) + b.scale(det::exp(log_b - scale));
        let magnitude = value.abs();
        Ok(Self {
            unit: if magnitude == 0.0 {
                C64::ZERO
            } else {
                value.scale(1.0 / magnitude)
            },
            log_magnitude: scale + det::ln(magnitude),
        })
    }
}

struct InteriorWeight {
    spectrum: LayerSpectralWeight,
    pressure: ScaledNumerator,
    normal_velocity: ScaledNumerator,
}

struct InteriorKernel {
    common: LayerKernel,
    // Preserve original face clearances through oscillatory phase and decay.
    depth_m: Dd,
    remaining_m: Dd,
}

impl InteriorKernel {
    fn new(
        layer: PlanarFluidLayer,
        frequency_hz: f64,
        source: [f64; 3],
        observer: [f64; 3],
    ) -> Result<Self, AeroacError> {
        let (coordinates, depth_m, remaining_m) =
            interior_coordinates(layer, frequency_hz, source, observer)?;
        Ok(Self {
            common: LayerKernel::new(layer, coordinates)?,
            depth_m,
            remaining_m,
        })
    }

    fn sample_weight(&self, region: usize, x: f64) -> Result<InteriorWeight, AeroacError> {
        let spectrum = self
            .common
            .sample_weight_at_depth(region, x, self.depth_m)?;
        let [_, z2, z3] = spectrum.zeta;
        let (cosine, sine) = LayerKernel::layer_functions_at(z2, self.remaining_m);
        let unit = |v: C64| {
            let length = v.abs();
            if length == 0.0 {
                C64::ZERO
            } else {
                C64::new(v.re / length, v.im / length)
            }
        };
        let rho2 = det::ln(self.common.layer.layer_state().density_kg_m3());
        let rho3 = det::ln(self.common.layer.receiving_state().density_kg_m3());
        let lc = det::ln(cosine.abs());
        let ls = det::ln(sine.abs());
        let l2 = det::ln(z2.abs());
        let l3 = det::ln(z3.abs());
        let pressure = ScaledNumerator::sum(
            C64::from_re(cosine.signum()),
            lc,
            minus_i(unit(z3)).scale(sine.signum()),
            l3 + ls + rho2 - rho3,
        )?;
        let normal_velocity = ScaledNumerator::sum(
            unit(z3).scale(cosine.signum()),
            l3 + lc - rho3,
            minus_i(unit(z2) * unit(z2)).scale(sine.signum()),
            2.0 * l2 + ls - rho2,
        )?;
        Ok(InteriorWeight {
            spectrum,
            pressure,
            normal_velocity,
        })
    }

    fn resolved(&self, region: usize, lower: f64, upper: f64) -> bool {
        let alpha = self.common.point(region, lower).zeta[1].im;
        // Once the remaining-depth numerator has reached its scaled
        // evanescent limit, so has the full-layer denominator. Their large
        // exponents already cancel analytically: only the actual interior
        // depth controls the remaining exponential variation. Keep the full
        // phase and denominator-variation checks, and retain the original
        // conservative depth before this common saturation point.
        let decay_depth = if alpha * self.remaining_m.to_f64() > 350.0 {
            self.depth_m.to_f64()
        } else {
            self.common.layer.thickness_m()
        };
        self.common
            .resolved_at_depth(region, lower, upper, decay_depth)
    }

    fn prepare_tail(
        &mut self,
        limits: FullWaveMonopoleLimits,
        receiving_density_kg_m3: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<[f64; 3], AeroacError> {
        // Above every critical circle, put yj=alpha_j/rho_j > 0 and
        // L=iD, D=y1*C(h)+N(h), C=cosh+(y3/y2)sinh, N=y3*cosh+y2*sinh.
        // Positive addition gives C(h)>=C(d)cosh(alpha2*z) and likewise N.
        // Thus |C(d)/L|<=2*rho1*exp(-alpha2*z)/alpha1 and
        // |V(d)/L|<=2*exp(-alpha2*z). The reciprocal exterior normal
        // derivative has the same bound: |Y1*C(d)/L|<=2*exp(-alpha2*z).
        // Only radial momentum changes with the actual receiving density.
        // With t=sqrt(q²-kmax²), alpha_j>=t,
        // q*dq=t*dt and q<=t+kmax, the integrated envelopes below follow.
        let h = (self.common.coordinates.a + self.depth_m).to_f64();
        let hk = h * self.common.critical[self.common.critical_count - 1];
        if !positive(h) || !hk.is_finite() {
            return Err(nonfinite("interior acoustic tail clearance"));
        }
        let lh = det::ln(h);
        let rho1 = det::ln(self.common.layer.incident_state().density_kg_m3());
        let receiving_density = det::ln(receiving_density_kg_m3);
        let scales = [
            det::ln(self.common.coordinates.omega) + rho1 - det::ln(PI) - lh,
            rho1 - receiving_density - det::ln(PI) - 2.0 * lh,
            -det::ln(PI) - 2.0 * lh,
        ];
        let logs = |s: f64| {
            [
                scales[0] - s,
                scales[1] - s + det::ln(s + hk + 1.0),
                scales[2] - s + det::ln(s + 1.0),
            ]
        };
        let goals = [
            det::ln(limits.absolute_pressure_tolerance_pa_per_m3_s) - det::ln(4.0),
            det::ln(limits.absolute_velocity_tolerance_per_m2) - det::ln(4.0) - 0.5 * det::ln(2.0),
        ];
        let accepted = |s| {
            let values = logs(s);
            values[0] <= goals[0] && values[1] <= goals[1] && values[2] <= goals[1]
        };
        let (mut low, mut high) = (0.0, 1.0);
        let mut admitted = false;
        for _ in 0..32 {
            poll(cancelled)?;
            if accepted(high) {
                admitted = true;
                break;
            }
            low = high;
            high *= 2.0;
        }
        if !admitted {
            return Err(invalid("interior evanescent tail cannot meet tolerance"));
        }
        for _ in 0..24 {
            poll(cancelled)?;
            let middle = low + (high - low) * 0.5;
            if accepted(middle) {
                high = middle;
            } else {
                low = middle;
            }
        }
        self.common.tail_extent = high / h;
        if !positive(self.common.tail_extent) {
            return Err(nonfinite("interior evanescent cutoff"));
        }
        let mut bound = logs(high).map(det::exp);
        if self.common.coordinates.range == 0.0 {
            bound[1] = 0.0;
        }
        Ok(bound)
    }
}

impl Spectrum<3> for InteriorKernel {
    fn regions(&self) -> usize {
        self.common.regions()
    }
    fn extent(&self) -> (f64, f64) {
        self.common.extent()
    }
    fn resolved(&self, region: usize, lower: f64, upper: f64) -> bool {
        InteriorKernel::resolved(self, region, lower, upper)
    }
    fn sample(&self, region: usize, x: f64) -> Result<[C64; 3], AeroacError> {
        let InteriorWeight {
            spectrum,
            pressure,
            normal_velocity,
        } = self.sample_weight(region, x)?;
        let argument = spectrum.q * self.common.coordinates.range;
        let rho2 = det::ln(self.common.layer.layer_state().density_kg_m3());
        let values = [
            spectrum.scaled(
                pressure.unit,
                pressure.log_magnitude + det::ln(self.common.coordinates.omega) - det::ln(TAU),
                j0(argument),
            ),
            spectrum.scaled(
                C64::new(-pressure.unit.im, pressure.unit.re),
                pressure.log_magnitude + det::ln(spectrum.q) - rho2 - det::ln(TAU),
                j1(argument),
            ),
            spectrum.scaled(
                normal_velocity.unit,
                normal_velocity.log_magnitude - det::ln(TAU),
                j0(argument),
            ),
        ];
        if values.iter().any(|v| !finite(*v)) {
            return Err(nonfinite("interior acoustic spectral field"));
        }
        Ok(values)
    }
}
