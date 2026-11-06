//! Original exterior volume stress observed in the layer's standing-wave field.
//!
//! With the parent's pressure numerator A and normal-velocity numerator B,
//! H=i/(2*pi*rho1) integral W*A*J0. The physical stress gives p=M:d_s d_s H.
//! Tangential receiving derivatives use A/(2*pi*omega*rho1*rho2), while
//! normal velocity uses B/(2*pi*omega*rho1). Differentiating the ORIGINAL
//! source twice leaves A/B intact and uses zeta1, never the layer wavenumber.
//! Thus each normal-velocity channel retains both internal directions and
//! cannot be reconstructed by dividing pressure by a traveling-wave impedance.

use super::super::super::force::j1_over_argument;
use super::super::super::quadrupole::{
    FullWaveQuadrupoleLimits, PreparedFullWaveQuadrupole, cartesian, j2_over_argument,
};
use super::super::super::reflection::{apply_quadrupole, field, line, prepare_reflected_tail};
use super::super::super::{
    QuadratureLimits, Spectrum, finite, integrate_spectrum, j0, j1, nonfinite, poll,
};
use super::super::{FullWaveLayerQuadrupoleDiagnostics, PlanarFluidLayer};
use super::{InteriorKernel, interior_coordinates, reduces_to_interface};
use crate::AeroacError;
use crate::curle3d::CoherentCompactField;
use crate::curle3d::quadrupole::QuadrupolePhasor;
use crate::observation::HarmonicFieldLine;
use core::f64::consts::{PI, SQRT_2, TAU};
use fs_math::{c64::C64, det};

/// Scope shared by the original source's coherent pressure/velocity observations.
pub const FULL_WAVE_LAYER_INTERIOR_QUADRUPOLE_SCOPE: &str = "Prescribed stationary \
symmetric 3D Lighthill stress moment [N m] in the original incident exterior; full \
spherical-wave pressure and vector velocity strictly inside one infinite planar \
homogeneous lossless fluid layer. Both internal directions, coherent returns, trace, \
physical shear, near, transverse, critical and evanescent fields are retained. \
Slower-than-both trapping layers refuse. Numerical quadrature estimates and analytic \
tail bounds are not interval certificates. No interior source, source-size or SGS \
validation, absorption, moving fluids, finite edges, source feedback or measured-SPL \
qualification.";

/// Cached 4-by-6 Cartesian transfer from original physical stress to a receiver
/// in the layer. Preparation shares native branches, resonance resolution and
/// cancellation. Warm observation uses the original fixed-order DD tensor
/// reduction without allocation, quadrature or another propagation delay.
#[derive(Clone, Debug)]
pub struct PreparedInteriorFullWaveQuadrupole {
    layer: PlanarFluidLayer,
    frequency_hz: f64,
    source_m: [f64; 3],
    observer_m: [f64; 3],
    transfer: [[C64; 6]; 4],
    diagnostics: FullWaveLayerQuadrupoleDiagnostics,
}

impl PreparedInteriorFullWaveQuadrupole {
    /// The source strictly precedes the first face; the receiver lies strictly
    /// inside the original finite layer. All fourteen physical tensor channels
    /// share quadrature nodes and the native standing-wave pressure/velocity
    /// numerators. Tolerances control complete per-N-m transfer norms before
    /// applying a particular source stress. Failure publishes no partial field.
    pub fn new(
        layer: PlanarFluidLayer,
        frequency_hz: f64,
        source_m: [f64; 3],
        observer_m: [f64; 3],
        limits: FullWaveQuadrupoleLimits,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, AeroacError> {
        poll(cancelled)?;
        QuadratureLimits::from(limits).admit()?;
        let (coordinates, depth, _) =
            interior_coordinates(layer, frequency_hz, source_m, observer_m)?;
        let (transfer, diagnostics) = if reduces_to_interface(layer) {
            // The original rear face disappears only when its two original
            // fluid states agree. The first face remains fully refracting.
            let original = PreparedFullWaveQuadrupole::new(
                layer.exterior,
                frequency_hz,
                source_m,
                observer_m,
                limits,
                cancelled,
            )?;
            (original.transfer, (*original.diagnostics()).into())
        } else {
            let r1 = det::ln(layer.incident_state().density_kg_m3());
            let r2 = det::ln(layer.layer_state().density_kg_m3());
            let w = det::ln(coordinates.omega);
            let log_pressure = -det::ln(TAU) - r1;
            let log_tangent = log_pressure - w - r2;
            let log_normal = log_pressure - w;
            let mut common = InteriorKernel::new(layer, frequency_hz, source_m, observer_m)?;
            // Above all critical circles, |A/L| <= 2rho1 exp(-alpha2*z)/alpha1
            // and |B/L| <= 2 exp(-alpha2*z). With q*dq=t*dt, alpha1>=t and
            // q,alpha1<=t+kmax, pressure has degree 2 and velocity degree 3.
            // Tangential velocity uses 1/rho2; normal velocity uses 1/rho1.
            // Their larger envelope safely controls all TEN velocity entries.
            let (extent, bound) = prepare_reflected_tail(
                (coordinates.a + depth).to_f64(),
                common.common.critical[common.common.critical_count - 1],
                limits.into(),
                [2, 3],
                [-det::ln(PI), -det::ln(PI) - w + (-r1).max(-r2)],
                [5, 14],
                cancelled,
            )?;
            common.common.tail_extent = extent;
            let [p, v] = bound;
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
            let kernel = InteriorQuadrupoleKernel {
                common,
                log_pressure,
                log_tangent,
                log_normal,
            };
            let result = integrate_spectrum::<14, 4>(&kernel, tail, limits.into(), cancelled)?;
            (
                cartesian(result.value, coordinates.tangent, layer.receiving_normal())?,
                FullWaveLayerQuadrupoleDiagnostics {
                    evaluations: result.evaluations,
                    maximum_depth_reached: result.maximum_depth_reached,
                    retained_panels: result.retained_panels,
                    transverse_wavenumber_cutoff_rad_m: result.cutoff,
                    maximum_bessel_argument: result.maximum_argument,
                    pressure_quadrature_error_pa_per_n_m: core::array::from_fn(|i| result.error[i]),
                    velocity_quadrature_error_m_s_per_n_m: core::array::from_fn(|i| {
                        result.error[i + 4]
                    }),
                    pressure_tail_bound_pa_per_n_m: core::array::from_fn(|i| tail[i]),
                    velocity_tail_bound_m_s_per_n_m: core::array::from_fn(|i| tail[i + 4]),
                    reduced_to_single_interface: false,
                },
            )
        };
        if transfer.iter().flatten().any(|value| !finite(*value)) {
            return Err(nonfinite("interior quadrupole Cartesian transfer"));
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
    /// Four pressure and ten velocity Frobenius entries in the original
    /// shared channel ordering. The reduction flag here denotes an absent
    /// original rear face (layer fluid equals receiving-exterior fluid).
    pub const fn diagnostics(&self) -> &FullWaveLayerQuadrupoleDiagnostics {
        &self.diagnostics
    }

    /// Original physical symmetric peak stress [N m], exp(-i omega t).
    /// All tensor terms interfere before signed active/reactive flux is formed.
    pub fn evaluate(&self, source: QuadrupolePhasor) -> Result<CoherentCompactField, AeroacError> {
        field(
            apply_quadrupole(&self.transfer, source)?,
            FULL_WAVE_LAYER_INTERIOR_QUADRUPOLE_SCOPE,
        )
    }

    pub fn evaluate_harmonic_line(
        &self,
        source: QuadrupolePhasor,
    ) -> Result<HarmonicFieldLine<3>, AeroacError> {
        Ok(line(
            self.frequency_hz,
            apply_quadrupole(&self.transfer, source)?,
        ))
    }
}

struct InteriorQuadrupoleKernel {
    common: InteriorKernel,
    log_pressure: f64,
    log_tangent: f64,
    log_normal: f64,
}

impl Spectrum<14> for InteriorQuadrupoleKernel {
    fn regions(&self) -> usize {
        <InteriorKernel as Spectrum<3>>::regions(&self.common)
    }
    fn extent(&self) -> (f64, f64) {
        <InteriorKernel as Spectrum<3>>::extent(&self.common)
    }
    fn resolved(&self, region: usize, lower: f64, upper: f64) -> bool {
        <InteriorKernel as Spectrum<3>>::resolved(&self.common, region, lower, upper)
    }
    fn sample(&self, region: usize, x: f64) -> Result<[C64; 14], AeroacError> {
        let weight = self.common.sample_weight(region, x)?;
        let argument = weight.spectrum.q * self.common.common.coordinates.range;
        if !argument.is_finite() {
            return Err(nonfinite("interior quadrupole Bessel argument"));
        }
        let zero = j0(argument);
        let one = j1(argument);
        let ratio = j1_over_argument(argument, one);
        let derivative = zero - ratio;
        let second_ratio = j2_over_argument(argument, zero, ratio);
        let third = one - second_ratio;
        let q = det::ln(weight.spectrum.q);
        let magnitude = weight.spectrum.zeta[0].abs();
        let z = det::ln(magnitude);
        let direction = if magnitude == 0.0 {
            C64::ZERO
        } else {
            C64::new(
                weight.spectrum.zeta[0].re / magnitude,
                weight.spectrum.zeta[0].im / magnitude,
            )
        };
        let square = direction * direction;
        let i = C64::new(0.0, 1.0);
        let pressure = |direction: C64, scale: f64, value: f64| {
            weight.spectrum.scaled(
                weight.pressure.unit * direction,
                weight.pressure.log_magnitude + scale,
                value,
            )
        };
        let normal = |direction: C64, scale: f64, value: f64| {
            weight.spectrum.scaled(
                weight.normal_velocity.unit * direction,
                weight.normal_velocity.log_magnitude + scale,
                value,
            )
        };
        // Form source/receiver derivative magnitudes in log space before the
        // common exponential. Both A and B retain their independent zeros and
        // cancellations, including normal motion at a pressure node.
        let p = self.log_pressure;
        let v = self.log_tangent;
        let n = self.log_normal;
        let values = [
            pressure(-i, p + 2.0 * q, derivative),
            pressure(-i, p + 2.0 * q, ratio),
            pressure(direction, p + q + z, SQRT_2 * one),
            pressure(-i * square, p + 2.0 * z, zero),
            pressure(C64::ONE, v + 3.0 * q, third),
            pressure(C64::ONE, v + 3.0 * q, second_ratio),
            pressure(-i * direction, v + 2.0 * q + z, SQRT_2 * derivative),
            pressure(square, v + q + 2.0 * z, one),
            pressure(C64::ONE, v + 3.0 * q, SQRT_2 * second_ratio),
            pressure(-i * direction, v + 2.0 * q + z, SQRT_2 * ratio),
            normal(-i, n + 2.0 * q, derivative),
            normal(-i, n + 2.0 * q, ratio),
            normal(direction, n + q + z, SQRT_2 * one),
            normal(-i * square, n + 2.0 * z, zero),
        ];
        if values.iter().any(|value| !finite(*value)) {
            return Err(nonfinite("interior quadrupole spectrum"));
        }
        Ok(values)
    }
}
