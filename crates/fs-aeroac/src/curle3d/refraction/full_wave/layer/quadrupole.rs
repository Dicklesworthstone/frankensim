//! Original volume-stress Hessian through a coherent finite fluid layer.
//!
//! With the parent's w=q*dq*E/L, the outgoing scalar Green function is
//! H=i/(2*pi*rho1) integral w*J0. The physical source is p=M_ij d_si d_sj H,
//! and u_l=d_xl p/(i*omega*rho3). Thus source-normal derivatives use zeta1
//! and receiver-normal derivatives use zeta3; the intermediate layer changes
//! the common weight, retaining every coherent internal return. The original
//! quadrupole owner supplies angular moments and Cartesian tensor contraction.
//! No finite differences, source remeshing, pressure/impedance approximation,
//! artificial damping, or finite image truncation enters the production field.

use super::super::force::j1_over_argument;
use super::super::quadrupole::{
    FullWaveQuadrupoleDiagnostics, FullWaveQuadrupoleLimits, PreparedFullWaveQuadrupole, cartesian,
    j2_over_argument,
};
use super::super::reflection::{field, line};
use super::super::{
    QuadratureLimits, Spectrum, finite, integrate_spectrum, j0, j1, nonfinite, poll,
};
use super::{LayerCoordinates, LayerKernel, PlanarFluidLayer};
use crate::AeroacError;
use crate::curle3d::CoherentCompactField;
use crate::curle3d::quadrupole::QuadrupolePhasor;
use crate::curle3d::statistics::ComplexSum;
use crate::observation::HarmonicFieldLine;
use core::f64::consts::{SQRT_2, TAU};
use fs_math::{c64::C64, det};

pub mod array;

/// Physical and numerical scope retained by every coherent observation.
pub const FULL_WAVE_LAYER_QUADRUPOLE_SCOPE: &str = "Prescribed stationary symmetric \
3D Lighthill stress moment [N m] in the original incident fluid; coherent spherical \
full-wave transmission through one infinite homogeneous lossless fluid layer, \
including every internal return, trace, physical shear, transverse, near, critical \
and evanescent fields. Slower-than-both trapping layers refuse. Estimated quadrature \
error and analytic tail bounds are not interval certificates. No source-size \
validation, unresolved SGS sound, absorption, mean flow, finite edges, source \
feedback or measured-SPL qualification.";

/// Original fourteen-channel Frobenius transfer errors per physical N m.
/// Pressure channels are (tt,ss,sqrt(2)tn,nn); velocity channels are
/// (ttt,tss,sqrt(2)ttn,tnn,sqrt(2)sts,sqrt(2)ssn,ntt,nss,sqrt(2)ntn,nnn).
/// Multiply their Euclidean norms by the source tensor's Frobenius norm,
/// counting off-diagonal stress entries twice, for a particular field estimate.
/// Numerical quadrature estimates and analytic omitted-spectrum bounds exclude
/// geometry, native Bessel and other floating-point errors.
#[derive(Clone, Copy, Debug)]
pub struct FullWaveLayerQuadrupoleDiagnostics {
    pub evaluations: usize,
    pub maximum_depth_reached: usize,
    pub retained_panels: usize,
    pub transverse_wavenumber_cutoff_rad_m: f64,
    pub maximum_bessel_argument: f64,
    pub pressure_quadrature_error_pa_per_n_m: [f64; 4],
    pub velocity_quadrature_error_m_s_per_n_m: [f64; 10],
    pub pressure_tail_bound_pa_per_n_m: [f64; 4],
    pub velocity_tail_bound_m_s_per_n_m: [f64; 10],
    /// An exact physical layer reduction used the original interface owner.
    /// Each observation owner documents its reduction domain. Unequal original
    /// speeds can still require quadrature after an interface reduction.
    pub reduced_to_single_interface: bool,
}

impl From<FullWaveQuadrupoleDiagnostics> for FullWaveLayerQuadrupoleDiagnostics {
    fn from(value: FullWaveQuadrupoleDiagnostics) -> Self {
        Self {
            evaluations: value.evaluations,
            maximum_depth_reached: value.maximum_depth_reached,
            retained_panels: value.retained_panels,
            transverse_wavenumber_cutoff_rad_m: value.transverse_wavenumber_cutoff_rad_m,
            maximum_bessel_argument: value.maximum_bessel_argument,
            pressure_quadrature_error_pa_per_n_m: value.pressure_quadrature_error_pa_per_n_m,
            velocity_quadrature_error_m_s_per_n_m: value.velocity_quadrature_error_m_s_per_n_m,
            pressure_tail_bound_pa_per_n_m: value.pressure_tail_bound_pa_per_n_m,
            velocity_tail_bound_m_s_per_n_m: value.velocity_tail_bound_m_s_per_n_m,
            reduced_to_single_interface: true,
        }
    }
}

/// Cached 4-by-6 Cartesian response to an original physical stress moment.
/// All six entries, including its trace and physical shears, remain coherent
/// on their original source clock. Warm application allocates nothing and
/// performs no quadrature, transcendental calls or extra propagation delay.
#[derive(Clone, Debug)]
pub struct PreparedLayeredFullWaveQuadrupole {
    layer: PlanarFluidLayer,
    frequency_hz: f64,
    source_m: [f64; 3],
    observer_m: [f64; 3],
    transfer: [[C64; 6]; 4],
    diagnostics: FullWaveLayerQuadrupoleDiagnostics,
}

impl PreparedLayeredFullWaveQuadrupole {
    /// Prepare an incident-side source and a receiver beyond the second face.
    /// All fourteen channels share the original branch-regularized layer
    /// quadrature and its resonance-resolution checks. Pressure-Hessian and
    /// velocity-tensor tolerances control the complete physical transfer before
    /// applying a particular source tensor. Failure publishes no partial field.
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
        let coordinates = LayerCoordinates::new(layer, frequency_hz, source_m, observer_m)?;
        let (transfer, diagnostics) = if layer.reduces_to_interface() {
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
            let log_pressure = -det::ln(TAU) - det::ln(layer.incident_state().density_kg_m3());
            let log_velocity = log_pressure
                - det::ln(coordinates.omega)
                - det::ln(layer.receiving_state().density_kg_m3());
            let mut common = LayerKernel::new(layer, coordinates)?;
            // The source Hessian and receiving third derivative require
            // quadratic/cubic high-q envelopes. Squared sqrt(2) channel
            // weights retain the complete symmetric-tensor Frobenius norms.
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
            let kernel = LayerQuadrupoleKernel {
                common,
                log_pressure,
                log_velocity,
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
            return Err(nonfinite("layered quadrupole Cartesian transfer"));
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
    pub const fn diagnostics(&self) -> &FullWaveLayerQuadrupoleDiagnostics {
        &self.diagnostics
    }

    fn amplitudes(&self, source: QuadrupolePhasor) -> Result<[C64; 4], AeroacError> {
        if source.stress_moment_n_m.iter().any(|value| !finite(*value)) {
            return Err(nonfinite("layered quadrupole stress moment"));
        }
        let mut sums = [ComplexSum::default(); 4];
        for (column, moment) in source.stress_moment_n_m.into_iter().enumerate() {
            if moment != C64::ZERO {
                for (row, sum) in self.transfer.iter().zip(&mut sums) {
                    sum.add(row[column] * moment)?;
                }
            }
        }
        let values = sums.map(ComplexSum::value);
        if values.iter().any(|value| !finite(*value)) {
            return Err(nonfinite("layered quadrupole field amplitude"));
        }
        Ok(values)
    }

    /// Physical peak moments [N m], with the original exp(-i omega t) clock.
    /// Signed active/reactive intensity is formed after coherent p/u reduction.
    pub fn evaluate(&self, source: QuadrupolePhasor) -> Result<CoherentCompactField, AeroacError> {
        field(self.amplitudes(source)?, FULL_WAVE_LAYER_QUADRUPOLE_SCOPE)
    }

    /// Native line for finite-window or sampled observation; propagation phase
    /// is already present and the source clock receives no additional delay.
    pub fn evaluate_harmonic_line(
        &self,
        source: QuadrupolePhasor,
    ) -> Result<HarmonicFieldLine<3>, AeroacError> {
        Ok(line(self.frequency_hz, self.amplitudes(source)?))
    }
}

struct LayerQuadrupoleKernel {
    common: LayerKernel,
    log_pressure: f64,
    log_velocity: f64,
}

impl Spectrum<14> for LayerQuadrupoleKernel {
    fn regions(&self) -> usize {
        <LayerKernel as Spectrum<3>>::regions(&self.common)
    }
    fn extent(&self) -> (f64, f64) {
        <LayerKernel as Spectrum<3>>::extent(&self.common)
    }
    fn resolved(&self, region: usize, lower: f64, upper: f64) -> bool {
        <LayerKernel as Spectrum<3>>::resolved(&self.common, region, lower, upper)
    }
    fn sample(&self, region: usize, x: f64) -> Result<[C64; 14], AeroacError> {
        let weight = self.common.sample_weight(region, x)?;
        let argument = weight.q * self.common.coordinates.range;
        if !argument.is_finite() {
            return Err(nonfinite("layered quadrupole Bessel argument"));
        }
        let zero = j0(argument);
        let one = j1(argument);
        let ratio = j1_over_argument(argument, one);
        let derivative = zero - ratio;
        let second_ratio = j2_over_argument(argument, zero, ratio);
        let third = one - second_ratio;
        let q_log = det::ln(weight.q);
        let magnitude = [weight.zeta[0].abs(), weight.zeta[2].abs()];
        let direction = core::array::from_fn::<_, 2, _>(|index| {
            let z = weight.zeta[2 * index];
            if magnitude[index] == 0.0 {
                C64::ZERO
            } else {
                C64::new(z.re / magnitude[index], z.im / magnitude[index])
            }
        });
        let [source_log, receiver_log] = magnitude.map(det::ln);
        let [source, receiver] = direction;
        let source_squared = source * source;
        let i = C64::new(0.0, 1.0);
        let p = self.log_pressure;
        let v = self.log_velocity;
        // Attach derivative magnitudes before exponentiation: a q^3 or
        // zeta1^2*zeta3 factor can restore a finite field after bare E underflows.
        let values = [
            weight.scaled(-i, p + 2.0 * q_log, derivative),
            weight.scaled(-i, p + 2.0 * q_log, ratio),
            weight.scaled(source, p + q_log + source_log, SQRT_2 * one),
            weight.scaled(-i * source_squared, p + 2.0 * source_log, zero),
            weight.scaled(C64::ONE, v + 3.0 * q_log, third),
            weight.scaled(C64::ONE, v + 3.0 * q_log, second_ratio),
            weight.scaled(
                -i * source,
                v + 2.0 * q_log + source_log,
                SQRT_2 * derivative,
            ),
            weight.scaled(source_squared, v + q_log + 2.0 * source_log, one),
            weight.scaled(C64::ONE, v + 3.0 * q_log, SQRT_2 * second_ratio),
            weight.scaled(-i * source, v + 2.0 * q_log + source_log, SQRT_2 * ratio),
            weight.scaled(-i * receiver, v + 2.0 * q_log + receiver_log, derivative),
            weight.scaled(-i * receiver, v + 2.0 * q_log + receiver_log, ratio),
            weight.scaled(
                source * receiver,
                v + q_log + source_log + receiver_log,
                SQRT_2 * one,
            ),
            weight.scaled(
                -i * source_squared * receiver,
                v + 2.0 * source_log + receiver_log,
                zero,
            ),
        ];
        if values.iter().any(|value| !finite(*value)) {
            return Err(nonfinite("layered quadrupole spectral field"));
        }
        Ok(values)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use fs_material::fluid_acoustics::FluidAcousticState;

    #[test]
    fn cubic_layer_tail_matches_independent_integration_by_parts() {
        let fluid = |rho, speed| FluidAcousticState::try_new(rho, speed).unwrap();
        let layer = PlanarFluidLayer::new(
            fluid(1.2, 300.0),
            fluid(2.4, 600.0),
            fluid(1.8, 400.0),
            [0.0; 3],
            [0.0, 0.0, 1.0],
            0.018,
        )
        .unwrap();
        let coordinates =
            LayerCoordinates::new(layer, 700.0, [0.0, 0.0, -0.025], [0.04, 0.0, 0.041]).unwrap();
        let omega = coordinates.omega;
        let mut kernel = LayerKernel::new(layer, coordinates).unwrap();
        let limits = QuadratureLimits::from(FullWaveQuadrupoleLimits {
            absolute_pressure_tolerance_pa_per_n_m: 2e-6,
            absolute_velocity_tolerance_m_s_per_n_m: 2e-9,
            relative_tolerance: 1e-10,
            maximum_evaluations: 262_144,
            maximum_depth: 40,
        });
        let scales = [1.0 / (TAU * 1.2), 1.0 / (TAU * omega * 1.2 * 1.8)];
        let actual = kernel
            .prepare_tail(limits, [2, 3], scales.map(det::ln), [5, 14], &mut || false)
            .unwrap();
        let h = (coordinates.a
            + coordinates.b
            + super::super::super::Dd::from_f64(layer.thickness_m()))
        .to_f64();
        let t = kernel.tail_extent;
        let k = kernel.critical[kernel.critical_count - 1];
        // I_d(T) = exp(-h*T)*(T+k)^d/h + d*I_(d-1)(T)/h.
        // This recurrence does not use the shifted polynomial in production.
        let decay = (-h * t).exp();
        let mut integral = decay / h;
        for degree in 1..=3 {
            integral = decay * (t + k).powi(degree) / h + f64::from(degree) * integral / h;
            if degree >= 2 {
                let channel = (degree - 2) as usize;
                let expected = scales[channel] * (2.0 / (1.0 / 1.2 + 1.0 / 1.8)) * integral;
                assert!((actual[channel] - expected).abs() <= 8e-14 * expected);
            }
        }
        assert!(
            kernel
                .prepare_tail(limits, [2, 4], scales.map(det::ln), [5, 14], &mut || false)
                .is_err()
        );
    }
}
