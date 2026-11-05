//! Analytic original-force derivatives of the interior standing-wave field.
//!
//! With the parent's WA/WB weights, pF=[i q WA J1,zeta1 WA J0]/(2pi rho1).
//! The tt,ss,tn velocity terms use WA/(2pi omega rho1 rho2); nt,nn use
//! WB/(2pi omega rho1). These are +i/(omega rho1) SOURCE derivatives of
//! the original monopole field, retaining both internal propagation directions.
use super::super::super::force::{admit, j1_over_argument};
use super::super::super::reflection::{
    apply_force, field, force_cartesian, line, prepare_reflected_tail,
};
use super::super::super::{
    FullWaveForceLimits, PreparedFullWaveForce, Spectrum, finite, integrate_spectrum, j0, j1,
    nonfinite, poll,
};
use super::super::{FullWaveLayerForceDiagnostics, PlanarFluidLayer};
use super::{InteriorKernel, interior_coordinates, reduces_to_interface};
use crate::AeroacError;
use crate::curle3d::CoherentCompactField;
use crate::observation::HarmonicFieldLine;
use core::f64::consts::{PI, TAU};
use fs_math::{c64::C64, det};

pub const FULL_WAVE_LAYER_INTERIOR_FORCE_SCOPE: &str = "Prescribed stationary compact 3D \
force ON the original incident fluid [N]; complete pressure and particle velocity \
strictly inside one infinite planar homogeneous lossless fluid layer, including coherent \
internal returns, standing, transverse, reactive and evanescent fields. Slower-than-both \
guided layers refuse. Numerical quadrature estimates and analytic tail bounds are not \
interval certificates; no interior source, absorption, moving fluids, finite edges, \
source feedback or measured-SPL qualification.";

/// A cached 4-by-3 Cartesian source transfer in the original layer fluid.
/// Every force orientation/phase uses fixed DD accumulation without allocation,
/// quadrature, finite differences or additional propagation delay.
#[derive(Clone, Debug)]
pub struct PreparedInteriorFullWaveForce {
    layer: PlanarFluidLayer,
    frequency_hz: f64,
    source_m: [f64; 3],
    observer_m: [f64; 3],
    transfer: [[C64; 3]; 4],
    diagnostics: FullWaveLayerForceDiagnostics,
}
impl PreparedInteriorFullWaveForce {
    /// Source strictly precedes the first face; observer is strictly inside.
    /// Pressure-vector and velocity-tensor norms retain the original per-N
    /// limits. All seven coefficients share spectral work and resonance checks.
    pub fn new(
        layer: PlanarFluidLayer,
        frequency_hz: f64,
        source_m: [f64; 3],
        observer_m: [f64; 3],
        limits: FullWaveForceLimits,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, AeroacError> {
        poll(cancelled)?;
        admit(limits)?;
        let (coordinates, depth, _) =
            interior_coordinates(layer, frequency_hz, source_m, observer_m)?;
        let (transfer, diagnostics) = if reduces_to_interface(layer) {
            let original = PreparedFullWaveForce::new(
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
            // Above all critical circles, |A/L| <= 2rho1 exp(-alpha2*s)/alpha1
            // and |B/L| <= 2 exp(-alpha2*s). Source and receiving derivatives
            // give degree-1 pressure and degree-2 velocity envelopes at a+s.
            let (extent, bound) = prepare_reflected_tail(
                (coordinates.a + depth).to_f64(),
                common.common.critical[common.common.critical_count - 1],
                limits.into(),
                [1, 2],
                [-det::ln(PI), -det::ln(PI) - w + (-r1).max(-r2)],
                [2, 5],
                cancelled,
            )?;
            common.common.tail_extent = extent;
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
            let kernel = InteriorForceKernel {
                common,
                log_pressure,
                log_tangent,
                log_normal,
            };
            let result = integrate_spectrum::<7, 2>(&kernel, tail, limits.into(), cancelled)?;
            (
                force_cartesian(result.value, coordinates.tangent, layer.receiving_normal())?,
                FullWaveLayerForceDiagnostics {
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
                    reduced_to_single_interface: false,
                },
            )
        };
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
    pub const fn diagnostics(&self) -> &FullWaveLayerForceDiagnostics {
        &self.diagnostics
    }
    /// Complex peak force ON the original incident fluid [N], exp(-i omega t).
    pub fn evaluate(&self, force: [C64; 3]) -> Result<CoherentCompactField, AeroacError> {
        field(
            apply_force(&self.transfer, force)?,
            FULL_WAVE_LAYER_INTERIOR_FORCE_SCOPE,
        )
    }
    pub fn evaluate_harmonic_line(
        &self,
        force: [C64; 3],
    ) -> Result<HarmonicFieldLine<3>, AeroacError> {
        Ok(line(self.frequency_hz, apply_force(&self.transfer, force)?))
    }
}
struct InteriorForceKernel {
    common: InteriorKernel,
    log_pressure: f64,
    log_tangent: f64,
    log_normal: f64,
}
impl Spectrum<7> for InteriorForceKernel {
    fn regions(&self) -> usize {
        <InteriorKernel as Spectrum<3>>::regions(&self.common)
    }
    fn extent(&self) -> (f64, f64) {
        <InteriorKernel as Spectrum<3>>::extent(&self.common)
    }
    fn resolved(&self, r: usize, a: f64, b: f64) -> bool {
        <InteriorKernel as Spectrum<3>>::resolved(&self.common, r, a, b)
    }
    fn sample(&self, region: usize, x: f64) -> Result<[C64; 7], AeroacError> {
        let w = self.common.sample_weight(region, x)?;
        let arg = w.spectrum.q * self.common.common.coordinates.range;
        let zero = j0(arg);
        let one = j1(arg);
        let ratio = j1_over_argument(arg, one);
        let q = det::ln(w.spectrum.q);
        let z = det::ln(w.spectrum.zeta[0].abs());
        let dz = direction(w.spectrum.zeta[0]);
        let i = C64::new(0.0, 1.0);
        let pressure = |direction: C64, scale: f64, value: f64| {
            w.spectrum.scaled(
                w.pressure.unit * direction,
                w.pressure.log_magnitude + scale,
                value,
            )
        };
        let normal = |direction: C64, scale: f64, value: f64| {
            w.spectrum.scaled(
                w.normal_velocity.unit * direction,
                w.normal_velocity.log_magnitude + scale,
                value,
            )
        };
        let values = [
            pressure(i, self.log_pressure + q, one),
            pressure(dz, self.log_pressure + z, zero),
            pressure(C64::ONE, self.log_tangent + 2.0 * q, zero - ratio),
            pressure(C64::ONE, self.log_tangent + 2.0 * q, ratio),
            pressure(i * dz, self.log_tangent + q + z, one),
            normal(i, self.log_normal + q, one),
            normal(dz, self.log_normal + z, zero),
        ];
        if values.iter().any(|v| !finite(*v)) {
            return Err(nonfinite("interior force spectrum"));
        }
        Ok(values)
    }
}

fn direction(z: C64) -> C64 {
    let magnitude = z.abs();
    if magnitude == 0.0 {
        C64::ZERO
    } else {
        C64::new(z.re / magnitude, z.im / magnitude)
    }
}
