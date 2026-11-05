//! Physical force derivatives of the original coherent finite-layer Green function.
//!
//! With the parent's physical weight w=q*dq*E/L and exp(-i omega t),
//! p/F = 1/(2*pi*rho1) integral w*[i*q*J1, zeta1*J0].
//! The five cylindrical velocity entries (tt,ss,tn,nt,nn) use
//! 1/(2*pi*omega*rho1*rho3) integral w*
//! [q²*(J0-J1/x), q²*J1/x, i*q*zeta1*J1, i*q*zeta3*J1, zeta1*zeta3*J0].
//! These are analytic source and receiver derivatives of the same Green
//! function; the finite layer's resonances and every internal return remain
//! coherent. No finite-difference source, pressure/impedance approximation,
//! artificial damping or truncation of internal reflections is introduced.

use super::super::force::{admit, j1_over_argument};
use super::super::{
    FullWaveForceDiagnostics, FullWaveForceLimits, PreparedFullWaveForce, Spectrum, finite,
    integrate_spectrum, j0, j1, nonfinite, poll,
};
use super::{LayerCoordinates, LayerKernel, PlanarFluidLayer};
use crate::AeroacError;
use crate::curle3d::CoherentCompactField;
use crate::curle3d::statistics::ComplexSum;
use crate::observation::HarmonicFieldLine;
use core::f64::consts::TAU;
use fs_math::{c64::C64, det};

/// Physical and numerical scope attached to every returned coherent field.
pub const FULL_WAVE_LAYER_FORCE_SCOPE: &str = "Prescribed stationary compact 3D force \
ON fluid [N]; coherent spherical-wave transmission through one infinite planar \
homogeneous lossless fluid layer between original stationary lossless fluids, \
including internal returns, transverse motion, reactive fields and evanescence. \
Slower-than-both trapping layers refuse. Numerical quadrature estimates and analytic \
spectral-tail bounds are not interval certificates; no absorption, moving media, \
finite edges, source feedback or measured-SPL qualification.";

/// Unit-force transfer estimates and analytic omitted-spectrum bounds.
/// Pressure entries are radial/normal; velocity entries are tt,ss,tn,nt,nn.
/// Their Euclidean norms respectively bound the pressure vector and velocity
/// tensor's Frobenius norm. Multiply by ||F|| for a prescribed force amplitude.
/// Embedded quadrature errors are estimates; geometry, Bessel, coordinate and
/// floating-point errors are outside both the estimates and the tail bounds.
#[derive(Clone, Copy, Debug)]
pub struct FullWaveLayerForceDiagnostics {
    pub evaluations: usize,
    pub maximum_depth_reached: usize,
    pub retained_panels: usize,
    pub transverse_wavenumber_cutoff_rad_m: f64,
    pub maximum_bessel_argument: f64,
    pub pressure_quadrature_error_pa_per_n: [f64; 2],
    pub velocity_quadrature_error_m_s_per_n: [f64; 5],
    pub pressure_tail_bound_pa_per_n: [f64; 2],
    pub velocity_tail_bound_m_s_per_n: [f64; 5],
    /// An absent or completely identical layer used the original force owner.
    /// An unequal-speed collapsed interface can still require quadrature.
    pub reduced_to_single_interface: bool,
}

impl From<FullWaveForceDiagnostics> for FullWaveLayerForceDiagnostics {
    fn from(value: FullWaveForceDiagnostics) -> Self {
        Self {
            evaluations: value.evaluations,
            maximum_depth_reached: value.maximum_depth_reached,
            retained_panels: value.retained_panels,
            transverse_wavenumber_cutoff_rad_m: value.transverse_wavenumber_cutoff_rad_m,
            maximum_bessel_argument: value.maximum_bessel_argument,
            pressure_quadrature_error_pa_per_n: value.pressure_quadrature_error_pa_per_n,
            velocity_quadrature_error_m_s_per_n: value.velocity_quadrature_error_m_s_per_n,
            pressure_tail_bound_pa_per_n: value.pressure_tail_bound_pa_per_n,
            velocity_tail_bound_m_s_per_n: value.velocity_tail_bound_m_s_per_n,
            reduced_to_single_interface: true,
        }
    }
}

/// Cached complete pressure/velocity transfer for an original physical force.
/// Twelve complex Cartesian coefficients serve every force orientation and
/// phase. Applying them allocates nothing and evaluates no transcendental
/// functions, quadrature, finite differences or additional propagation delay.
#[derive(Clone, Debug)]
pub struct PreparedLayeredFullWaveForce {
    layer: PlanarFluidLayer,
    frequency_hz: f64,
    source_m: [f64; 3],
    observer_m: [f64; 3],
    transfer: [[C64; 3]; 4],
    diagnostics: FullWaveLayerForceDiagnostics,
}

impl PreparedLayeredFullWaveForce {
    /// Prepare a force in the incident fluid and a receiver strictly beyond
    /// the second layer face. The original force limits control pressure-vector
    /// and velocity-tensor norms, preserving their different physical units.
    /// All seven channels share bounded quadrature nodes and the parent's
    /// resonance-resolution and outgoing-branch rules. Failure publishes no
    /// partial transfer. Guided-pole regimes are refused by PlanarFluidLayer.
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
        let coordinates = LayerCoordinates::new(layer, frequency_hz, source_m, observer_m)?;
        let (transfer, diagnostics) = if layer.reduces_to_interface() {
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
            let rho1 = layer.incident_state().density_kg_m3();
            let rho3 = layer.receiving_state().density_kg_m3();
            let log_pressure = -det::ln(TAU) - det::ln(rho1);
            let log_velocity = log_pressure - det::ln(coordinates.omega) - det::ln(rho3);
            let mut common = LayerKernel::new(layer, coordinates)?;
            let bound = common.prepare_tail(
                limits.into(),
                [1, 2],
                [log_pressure, log_velocity],
                [2, 5],
                cancelled,
            )?;
            let radial = coordinates.range != 0.0;
            let tail = [
                if radial { bound[0] } else { 0.0 },
                bound[0],
                bound[1],
                bound[1] * 0.5,
                if radial { bound[1] } else { 0.0 },
                if radial { bound[1] } else { 0.0 },
                bound[1],
            ];
            let kernel = LayerForceKernel {
                common,
                log_pressure,
                log_velocity,
            };
            let result = integrate_spectrum::<7, 2>(&kernel, tail, limits.into(), cancelled)?;
            let [pt, pn, vtt, vss, vtn, vnt, vnn] = result.value;
            let tangent = coordinates.tangent;
            let normal = layer.receiving_normal();
            let mut field = [[C64::ZERO; 3]; 4];
            for column in 0..3 {
                field[0][column] = pt.scale(tangent[column]) + pn.scale(normal[column]);
                for row in 0..3 {
                    let nn = normal[row] * normal[column];
                    let identity = if row == column { 1.0 } else { 0.0 };
                    let mut sum = ComplexSum::default();
                    sum.add(vss.scale(identity - nn))?;
                    sum.add((vtt - vss).scale(tangent[row] * tangent[column]))?;
                    sum.add(vtn.scale(tangent[row] * normal[column]))?;
                    sum.add(vnt.scale(normal[row] * tangent[column]))?;
                    sum.add(vnn.scale(nn))?;
                    field[row + 1][column] = sum.value();
                }
            }
            (
                field,
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
        if transfer.iter().flatten().any(|v| !finite(*v)) {
            return Err(nonfinite("full-wave layered force transfer"));
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
    pub const fn diagnostics(&self) -> &FullWaveLayerForceDiagnostics {
        &self.diagnostics
    }

    fn amplitudes(&self, force_on_fluid_n: [C64; 3]) -> Result<[C64; 4], AeroacError> {
        if force_on_fluid_n.iter().any(|v| !finite(*v)) {
            return Err(nonfinite("full-wave layered force phasor"));
        }
        let mut sums = [ComplexSum::default(); 4];
        for (column, amplitude) in force_on_fluid_n.into_iter().enumerate() {
            if amplitude != C64::ZERO {
                for (row, sum) in self.transfer.iter().zip(sums.iter_mut()) {
                    sum.add(row[column] * amplitude)?;
                }
            }
        }
        let values = sums.map(ComplexSum::value);
        if values.iter().any(|v| !finite(*v)) {
            return Err(nonfinite("full-wave layered force field"));
        }
        Ok(values)
    }

    /// Complex peak force ON the original fluid in newtons, exp(-i omega t).
    /// The full field retains its signed active and reactive intensity.
    pub fn evaluate(
        &self,
        force_on_fluid_n: [C64; 3],
    ) -> Result<CoherentCompactField, AeroacError> {
        let [pressure, ux, uy, uz] = self.amplitudes(force_on_fluid_n)?;
        let velocity = [ux, uy, uz];
        let intensity = velocity.map(|v| pressure.scale(0.5) * v.conj());
        if intensity.iter().any(|v| !finite(*v)) {
            return Err(nonfinite("full-wave layered force intensity"));
        }
        Ok(CoherentCompactField {
            pressure_pa: pressure,
            particle_velocity_m_s: velocity,
            active_intensity_w_m2: intensity.map(|v| v.re),
            reactive_intensity_w_m2: intensity.map(|v| v.im),
            scope: FULL_WAVE_LAYER_FORCE_SCOPE,
        })
    }

    /// Native finite-window observer line with all propagation phase included.
    pub fn evaluate_harmonic_line(
        &self,
        force_on_fluid_n: [C64; 3],
    ) -> Result<HarmonicFieldLine<3>, AeroacError> {
        let [pressure, ux, uy, uz] = self.amplitudes(force_on_fluid_n)?;
        Ok(HarmonicFieldLine {
            frequency_hz: self.frequency_hz,
            pressure_pa: pressure,
            particle_velocity_m_s: [ux, uy, uz],
        })
    }
}

struct LayerForceKernel {
    common: LayerKernel,
    log_pressure: f64,
    log_velocity: f64,
}

impl Spectrum<7> for LayerForceKernel {
    fn regions(&self) -> usize {
        <LayerKernel as Spectrum<3>>::regions(&self.common)
    }
    fn extent(&self) -> (f64, f64) {
        <LayerKernel as Spectrum<3>>::extent(&self.common)
    }
    fn resolved(&self, region: usize, lower: f64, upper: f64) -> bool {
        <LayerKernel as Spectrum<3>>::resolved(&self.common, region, lower, upper)
    }
    fn sample(&self, region: usize, x: f64) -> Result<[C64; 7], AeroacError> {
        let weight = self.common.sample_weight(region, x)?;
        let argument = weight.q * self.common.coordinates.range;
        if !argument.is_finite() {
            return Err(nonfinite("layered force Bessel argument"));
        }
        let jzero = j0(argument);
        let jone = j1(argument);
        let ratio = j1_over_argument(argument, jone);
        let q_log = det::ln(weight.q);
        let magnitude = [weight.zeta[0].abs(), weight.zeta[2].abs()];
        let direction = core::array::from_fn::<_, 2, _>(|i| {
            let z = weight.zeta[2 * i];
            if magnitude[i] == 0.0 {
                C64::ZERO
            } else {
                C64::new(z.re / magnitude[i], z.im / magnitude[i])
            }
        });
        let vertical_log = magnitude.map(det::ln);
        let i = C64::new(0.0, 1.0);
        let values = [
            weight.scaled(i, self.log_pressure + q_log, jone),
            weight.scaled(direction[0], self.log_pressure + vertical_log[0], jzero),
            weight.scaled(C64::ONE, self.log_velocity + 2.0 * q_log, jzero - ratio),
            weight.scaled(C64::ONE, self.log_velocity + 2.0 * q_log, ratio),
            weight.scaled(
                i * direction[0],
                self.log_velocity + q_log + vertical_log[0],
                jone,
            ),
            weight.scaled(
                i * direction[1],
                self.log_velocity + q_log + vertical_log[1],
                jone,
            ),
            weight.scaled(
                direction[0] * direction[1],
                self.log_velocity + vertical_log[0] + vertical_log[1],
                jzero,
            ),
        ];
        if values.iter().any(|v| !finite(*v)) {
            return Err(nonfinite("layered force spectral field"));
        }
        Ok(values)
    }
}
