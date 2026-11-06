//! Original symmetric volume stresses returned by a finite refracting layer.
//!
//! The parent supplies W=q*dq*R_layer*exp(i*zeta1*(a+b))/zeta1. The reflected
//! scalar Green function is i/(4*pi) integral W*J0, exactly as for one fluid
//! plane. Its source Hessian and returning receiver derivative therefore use
//! the original fourteen angular channels unchanged. Only the native scalar
//! reflection weight contains the layer, including every coherent return.
//! Source and receiver normal covectors have opposite signs; physical trace,
//! shear, near and reactive fields survive the Cartesian contraction.

use super::super::super::quadrupole::{FullWaveQuadrupoleLimits, cartesian};
use super::super::super::reflection::{
    PreparedReflectedFullWaveQuadrupole, ReflectionCoordinates, apply_quadrupole, field, line,
    quadrupole_channels, singular_total,
};
use super::super::super::{
    QuadratureLimits, Spectrum, finite, integrate_spectrum, nonfinite, poll,
};
use super::super::{FullWaveLayerQuadrupoleDiagnostics, PlanarFluidLayer};
use super::{LayerReflectionKernel, add_direct, reduces_to_first_interface};
use crate::AeroacError;
use crate::curle3d::CoherentCompactField;
use crate::curle3d::quadrupole::{PreparedQuadrupoleObserver, QuadrupolePhasor};
use crate::observation::HarmonicFieldLine;
use core::f64::consts::{PI, SQRT_2};
use fs_math::{c64::C64, det};

/// Returned field only, retaining the original complete symmetric stress.
pub const FULL_WAVE_LAYER_REFLECTED_QUADRUPOLE_SCOPE: &str = "Prescribed stationary \
symmetric 3D Lighthill stress moment [N m]; reflected spherical-wave field from \
one infinite homogeneous lossless fluid layer, including all coherent internal \
returns, physical trace/shear, transverse, near, critical and evanescent fields. \
Slower-than-both guided layers refuse. Quadrature estimates and analytic tails \
are not interval certificates. No source-size or SGS validation, absorption, \
mean flow, finite edges, source loading or measured-SPL qualification.";

/// Original direct and returned pressure/velocity interfere before energy.
pub const FULL_WAVE_LAYER_TOTAL_QUADRUPOLE_SCOPE: &str = "Prescribed stationary \
symmetric 3D Lighthill stress moment [N m]; coherent direct-plus-finite-layer \
reflected spherical-wave field in the original incident fluid. Numerical relative \
control applies to the reflected transfer; its absolute estimates carry to total \
fields without relative guarantees at interference nulls. Guided layers and \
singular direct self-fields refuse. No source-size or SGS validation, absorption, \
mean flow, finite edges, source loading, interval certificate or measured-SPL \
qualification.";

type Transfer = [[C64; 6]; 4];

/// Two cached Cartesian maps preserve all original tensor entries and phases.
/// Warm application uses the original fixed-order DD contraction without
/// allocation, quadrature, transcendental calls or another propagation delay.
/// A coincident source/receiver retains its finite returned field; the total
/// prescribed point-source field remains singular and explicitly refuses.
#[derive(Clone, Debug)]
pub struct PreparedLayeredReflectedFullWaveQuadrupole {
    layer: PlanarFluidLayer,
    frequency_hz: f64,
    source_m: [f64; 3],
    observer_m: [f64; 3],
    reflected: Transfer,
    total: Option<Transfer>,
    diagnostics: FullWaveLayerQuadrupoleDiagnostics,
}

impl PreparedLayeredReflectedFullWaveQuadrupole {
    /// Source and receiver lie strictly before the original first face.
    /// Limits control reflected pressure-Hessian and velocity-tensor Frobenius
    /// norms per N m. All fourteen channels share the native layer denominator,
    /// critical maps, resonance guard and bounded quadrature. A zero-thickness
    /// layer, or one matching the receiving half-space, delegates exactly to
    /// the original single-plane owner, including its direct/total convention.
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
        let coordinates =
            ReflectionCoordinates::new(layer.exterior, frequency_hz, source_m, observer_m)?;
        let (reflected, total, diagnostics) = if reduces_to_first_interface(layer) {
            let original = PreparedReflectedFullWaveQuadrupole::new(
                layer.exterior,
                frequency_hz,
                source_m,
                observer_m,
                limits,
                cancelled,
            )?;
            let mut reflected = [[C64::ZERO; 6]; 4];
            let mut total = original.has_total_field().then_some([[C64::ZERO; 6]; 4]);
            for column in 0..6 {
                poll(cancelled)?;
                let source = basis(column);
                let values = amplitudes(original.evaluate_reflected_harmonic_line(source)?);
                for row in 0..4 {
                    reflected[row][column] = values[row];
                }
                if let Some(total) = &mut total {
                    let values = amplitudes(original.evaluate_total_harmonic_line(source)?);
                    for row in 0..4 {
                        total[row][column] = values[row];
                    }
                }
            }
            (reflected, total, (*original.diagnostics()).into())
        } else {
            let log_pressure = -det::ln(4.0 * PI);
            let log_velocity = log_pressure
                - det::ln(coordinates.omega)
                - det::ln(layer.incident_state().density_kg_m3());
            let mut common = LayerReflectionKernel::new(layer, coordinates)?;
            // Above all three critical circles, |R_layer|<=1. Quadratic
            // source and cubic receiving derivatives use the actual exterior
            // return clearance a+b, with physical shear Frobenius weights.
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
            let kernel = LayerReflectedQuadrupoleKernel {
                common,
                log_pressure,
                log_velocity,
            };
            let result = integrate_spectrum::<14, 4>(&kernel, tail, limits.into(), cancelled)?;
            let reflected = cartesian(result.value, coordinates.tangent, layer.receiving_normal())?;
            let total = if source_m == observer_m {
                None
            } else {
                let direct = PreparedQuadrupoleObserver::new_cancellable(
                    layer.incident_state(),
                    frequency_hz,
                    observer_m,
                    &[source_m],
                    1,
                    &mut *cancelled,
                )?;
                let mut total = [[C64::ZERO; 6]; 4];
                for column in 0..6 {
                    poll(cancelled)?;
                    let values = amplitudes(direct.evaluate_harmonic_line(&[basis(column)])?);
                    let values =
                        add_direct(values, core::array::from_fn(|row| reflected[row][column]))?;
                    for row in 0..4 {
                        total[row][column] = values[row];
                    }
                }
                Some(total)
            };
            (
                reflected,
                total,
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
        if reflected.iter().flatten().any(|value| !finite(*value))
            || total
                .as_ref()
                .is_some_and(|t| t.iter().flatten().any(|value| !finite(*value)))
        {
            return Err(nonfinite(
                "layer reflected/total quadrupole Cartesian transfer",
            ));
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
    pub const fn has_total_field(&self) -> bool {
        self.total.is_some()
    }
    /// Original four-pressure/ten-velocity Frobenius diagnostics for the
    /// returned transfer. Total fields inherit absolute reflected errors;
    /// a direct/return interference null has no relative-error guarantee.
    pub const fn diagnostics(&self) -> &FullWaveLayerQuadrupoleDiagnostics {
        &self.diagnostics
    }

    /// Original physical symmetric stress [N m], exp(-i omega t).
    pub fn evaluate_reflected(
        &self,
        source: QuadrupolePhasor,
    ) -> Result<CoherentCompactField, AeroacError> {
        field(
            apply_quadrupole(&self.reflected, source)?,
            FULL_WAVE_LAYER_REFLECTED_QUADRUPOLE_SCOPE,
        )
    }

    /// Coherent direct/returned pressure and velocity, then signed intensity.
    pub fn evaluate_total(
        &self,
        source: QuadrupolePhasor,
    ) -> Result<CoherentCompactField, AeroacError> {
        field(
            apply_quadrupole(self.total.as_ref().ok_or_else(singular_total)?, source)?,
            FULL_WAVE_LAYER_TOTAL_QUADRUPOLE_SCOPE,
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

fn basis(column: usize) -> QuadrupolePhasor {
    let mut source = QuadrupolePhasor::default();
    source.stress_moment_n_m[column] = C64::ONE;
    source
}

fn amplitudes(line: HarmonicFieldLine<3>) -> [C64; 4] {
    let [x, y, z] = line.particle_velocity_m_s;
    [line.pressure_pa, x, y, z]
}

struct LayerReflectedQuadrupoleKernel {
    common: LayerReflectionKernel,
    log_pressure: f64,
    log_velocity: f64,
}

impl Spectrum<14> for LayerReflectedQuadrupoleKernel {
    fn regions(&self) -> usize {
        <LayerReflectionKernel as Spectrum<3>>::regions(&self.common)
    }
    fn extent(&self) -> (f64, f64) {
        <LayerReflectionKernel as Spectrum<3>>::extent(&self.common)
    }
    fn resolved(&self, region: usize, lower: f64, upper: f64) -> bool {
        <LayerReflectionKernel as Spectrum<3>>::resolved(&self.common, region, lower, upper)
    }
    fn sample(&self, region: usize, x: f64) -> Result<[C64; 14], AeroacError> {
        quadrupole_channels(
            self.common.sample_weight(region, x)?,
            self.common.common.coordinates.range,
            self.log_pressure,
            self.log_velocity,
        )
    }
}
