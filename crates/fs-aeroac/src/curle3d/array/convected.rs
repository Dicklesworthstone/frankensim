//! Coherent monopoles in a homogeneous fluid with uniform subsonic velocity U.
//!
//! For exp(-i omega t), M=U/c, beta²=1-|M|² and r=x-y, the outgoing
//! convected Green function is G=exp(i k d)/(4 pi R), where
//! R=sqrt(beta² |r|²+(M.r)²), d=(R-M.r)/beta². A physical volume injection
//! Q [m³/s] has potential phi=-Q G, velocity u=grad(phi), and pressure
//! p=-rho(-i omega+U.grad)phi. Pressure and travel time are supplied by the
//! original [`UniformFlowCompactObserver`], which already implements this
//! source law through the Galilean moving-source owner. This array adds the
//! potential's particle velocity, coherent summation and transported energy.
//! See Wang, *Flow effects in computational aeroacoustics* (2020),
//! equations 3.28–3.30: <https://pure.tue.nl/ws/portalfiles/portal/149260007/20200511_Z.Wang.pdf>.
//!
//! Q injects fluid at the mean velocity, with no added perturbation force or
//! entropy. A compact force in flow also produces a convected vortical wake;
//! replacing its velocity by grad(p)/(i omega rho) is not admitted here.
//! Source and receiver positions are fixed in the laboratory frame. Density
//! and sound speed come from one original checked material snapshot.
//!
//! Away from sources the linearized Euler equations give
//! d_t E + div(U E + p u)=0, E=p²/(2 rho c²)+rho |u|²/2. We report this
//! positive quadratic perturbation energy and its transport separately from
//! the pressure-work intensity p u. These are not mean-flow energy change,
//! a turbulence model or a source-validation claim. The ideal-plane adapter
//! adds exact image reflections when this SAME flow is tangent to every wall.

use super::super::statistics::{ComplexSum, StatisticsSum, finite, invalid, nonfinite};
use super::{CoherenceBudget, CoherentCompactField, CompactFieldStatistics};
use crate::AeroacError;
use crate::curle3d::{
    CompactSourcePhasor,
    moving::{MovingCompactLimits, flow::UniformFlowCompactObserver},
};
use fs_material::fluid_acoustics::FluidAcousticState;
use fs_math::{c64::C64, dd::Dd, det};

mod reflection;

/// Applicability shared by deterministic and statistical observations.
pub const CONVECTED_MONOPOLE_SCOPE: &str = "Prescribed compact volume injection at the uniform mean velocity; \
homogeneous subsonic lossless linear fluid, fixed sources and receivers; no extra force/entropy, \
source-size validation, nonuniform/vortical flow, scattering, absorption, mutual loading or measured-SPL qualification.";

/// Explicit distance and subsonic margins for every original source path.
#[derive(Clone, Copy, Debug)]
pub struct ConvectedMonopoleLimits {
    /// Inclusive speed/c limit in [0,1); zero admits only quiescent fluid.
    pub maximum_mach: f64,
    /// Minimum source/receiver separation [m], positive and finite.
    pub minimum_distance_m: f64,
}

/// Coherent field and the transported quadratic perturbation energy.
#[derive(Clone, Debug)]
pub struct ConvectedMonopoleField {
    /// The actual pressure, velocity and pressure-work intensity. Its active
    /// intensity is only <p u>, not the total transport below.
    pub field: CoherentCompactField,
    /// <p²>/(2 rho c²) + rho <|u|²>/2 [J/m³].
    pub mean_perturbation_energy_density_j_m3: f64,
    /// <p u> + U <E> [W/m²], in the original laboratory coordinates.
    pub mean_perturbation_energy_flux_w_m2: [f64; 3],
}

/// Explicit source second moments with the same uniform-flow energy balance.
#[derive(Clone, Debug)]
pub struct ConvectedMonopoleStatistics {
    /// Original spectral-line pressure/velocity covariance and <p u>.
    pub field: CompactFieldStatistics,
    pub mean_perturbation_energy_density_j_m3: f64,
    pub mean_perturbation_energy_flux_w_m2: [f64; 3],
}

#[derive(Clone, Debug)]
struct Path {
    source_m: [f64; 3],
    delay_s: f64,
    // Pressure, ux, uy, uz per unit volume flow. Phase is already included.
    field: [C64; 4],
}

/// Prepared linear field from fixed physical monopoles in uniform mean flow.
///
/// Preparation stores four complex coefficients per source. Application is
/// allocation-free O(N), retains each source's physical phase, and uses the
/// same double-double field sums and coherence statistics as the quiescent
/// source owner. No source covariance matrix or boundary solver is introduced.
#[derive(Clone, Debug)]
pub struct PreparedConvectedMonopoleObserver {
    acoustic: FluidAcousticState,
    flow_velocity_m_s: [f64; 3],
    observer_m: [f64; 3],
    frequency_hz: f64,
    paths: Vec<Path>,
    limits: ConvectedMonopoleLimits,
    // Images are folded into each original source, not statistical sources.
    reflection_planes: Vec<super::IdealAcousticPlane>,
}

impl PreparedConvectedMonopoleObserver {
    /// Prepare original material, physical mean flow and all source positions.
    pub fn new(
        acoustic: FluidAcousticState,
        flow_velocity_m_s: [f64; 3],
        frequency_hz: f64,
        observer_m: [f64; 3],
        source_positions_m: &[[f64; 3]],
        max_sources: usize,
        limits: ConvectedMonopoleLimits,
    ) -> Result<Self, AeroacError> {
        Self::new_cancellable(
            acoustic,
            flow_velocity_m_s,
            frequency_hz,
            observer_m,
            source_positions_m,
            max_sources,
            limits,
            || false,
        )
    }

    /// Bounded preparation; cancellation or refusal publishes no partial bank.
    #[allow(clippy::too_many_arguments)]
    pub fn new_cancellable(
        acoustic: FluidAcousticState,
        flow_velocity_m_s: [f64; 3],
        frequency_hz: f64,
        observer_m: [f64; 3],
        source_positions_m: &[[f64; 3]],
        max_sources: usize,
        limits: ConvectedMonopoleLimits,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, AeroacError> {
        poll(&mut cancelled)?;
        let c = acoustic.sound_speed_m_s();
        let omega = core::f64::consts::TAU * frequency_hz;
        let k = omega / c;
        let mach = flow_velocity_m_s.map(|v| v / c);
        let speed = norm(mach);
        if source_positions_m.len() > max_sources
            || !limits.maximum_mach.is_finite()
            || !(0.0..1.0).contains(&limits.maximum_mach)
            || !limits.minimum_distance_m.is_finite()
            || limits.minimum_distance_m <= 0.0
            || !speed.is_finite()
            || speed > limits.maximum_mach
            || [frequency_hz, omega, k, acoustic.density_kg_m3() * omega]
                .iter()
                .any(|v| !v.is_finite() || *v <= 0.0)
            || observer_m.iter().any(|v| !v.is_finite())
        {
            return Err(invalid(
                "convected monopole frequency, subsonic flow, geometry or source budget",
            ));
        }
        let beta_squared = (1.0 - speed) * (1.0 + speed);
        let mut paths = Vec::new();
        paths
            .try_reserve_exact(source_positions_m.len())
            .map_err(|_| invalid("convected monopole allocation exceeds available memory"))?;
        for &source_m in source_positions_m {
            poll(&mut cancelled)?;
            if source_m.iter().any(|v| !v.is_finite()) {
                return Err(nonfinite("convected monopole source position"));
            }
            let pressure_owner = UniformFlowCompactObserver::new(
                acoustic,
                flow_velocity_m_s,
                observer_m,
                source_m,
                MovingCompactLimits {
                    maximum_mach: limits.maximum_mach,
                    minimum_distance_m: limits.minimum_distance_m,
                },
            )?;
            let displacement = core::array::from_fn::<_, 3, _>(|i| observer_m[i] - source_m[i]);
            let range = norm(displacement);
            if !range.is_finite() || range < limits.minimum_distance_m {
                return Err(invalid("convected monopole source/receiver separation"));
            }
            let direction = displacement.map(|v| v / range);
            let projection = dot(mach, direction);
            let stretch = det::hypot(projection, det::sqrt(beta_squared));
            let delay_s = pressure_owner.propagation_delay_s();
            // The original pressure owner resolves the convected arrival; the
            // new velocity gradient uses exactly that same physical phase.
            let delay_factor = delay_s / (range / c);
            let inverse_r = (1.0 / range) / stretch;
            let phase = omega * delay_s;
            if [delay_s, inverse_r, phase]
                .iter()
                .any(|v| !v.is_finite() || *v <= 0.0)
            {
                return Err(invalid(
                    "convected monopole propagation is not representable",
                ));
            }
            let (sin, cos) = det::sin_cos(phase);
            let green = C64::new(cos, sin).scale((0.25 / core::f64::consts::PI) * inverse_r);
            let mut field = [C64::ZERO; 4];
            // grad(R) and grad(d) are dimensionless. The latter expression
            // avoids subtracting nearly identical gradients and dividing by beta².
            for axis in 0..3 {
                let grad_r = (beta_squared * direction[axis] + projection * mach[axis]) / stretch;
                let grad_d = (direction[axis] - mach[axis] * delay_factor) / stretch;
                field[axis + 1] = green * C64::new(grad_r * inverse_r, -k * grad_d);
            }
            field[0] = pressure_owner
                .pressure_phasor(frequency_hz, CompactSourcePhasor::monopole(C64::ONE))?;
            if field.iter().any(|&value| !finite(value))
                || green == C64::ZERO
                || field[0] == C64::ZERO
                || field[1..].iter().all(|&value| value == C64::ZERO)
            {
                return Err(nonfinite("convected monopole field transfer"));
            }
            paths.push(Path {
                source_m,
                delay_s,
                field,
            });
        }
        poll(&mut cancelled)?;
        Ok(Self {
            acoustic,
            flow_velocity_m_s,
            observer_m,
            frequency_hz,
            paths,
            limits,
            reflection_planes: Vec::new(),
        })
    }

    /// Original coherent density and adiabatic sound speed.
    pub const fn acoustic(&self) -> FluidAcousticState {
        self.acoustic
    }
    /// Actual uniform laboratory-frame fluid velocity [m/s].
    pub const fn flow_velocity_m_s(&self) -> [f64; 3] {
        self.flow_velocity_m_s
    }
    pub const fn observer_position_m(&self) -> [f64; 3] {
        self.observer_m
    }
    pub const fn frequency_hz(&self) -> f64 {
        self.frequency_hz
    }
    pub fn source_count(&self) -> usize {
        self.paths.len()
    }
    pub fn source_position_m(&self, index: usize) -> Option<[f64; 3]> {
        self.paths.get(index).map(|path| path.source_m)
    }
    /// Physical DIRECT-path travel time d/c, including directional convection.
    /// This phase is already included in the prepared coefficients.
    pub fn source_delay_s(&self, index: usize) -> Option<f64> {
        self.paths.get(index).map(|path| path.delay_s)
    }

    /// Original admitted boundaries; empty for free-space propagation.
    pub fn reflection_planes(&self) -> &[super::IdealAcousticPlane] {
        &self.reflection_planes
    }

    /// Physical propagation model retained by every field/statistical query.
    pub const fn scope(&self) -> &'static str {
        if self.reflection_planes.is_empty() {
            CONVECTED_MONOPOLE_SCOPE
        } else {
            reflection::CONVECTED_IDEAL_BOUNDARY_SCOPE
        }
    }

    /// Coherently propagate peak volume flows [m³/s] on one source clock.
    pub fn evaluate(&self, flows: &[C64]) -> Result<ConvectedMonopoleField, AeroacError> {
        self.evaluate_cancellable(flows, || false)
    }

    /// Allocation-free source summation with a cancellation poll per path.
    pub fn evaluate_cancellable(
        &self,
        flows: &[C64],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<ConvectedMonopoleField, AeroacError> {
        let values = self.amplitudes(flows, &mut cancelled)?;
        let mut sum = StatisticsSum::default();
        sum.add(values)?;
        let moments = self.finish(sum)?;
        poll(&mut cancelled)?;
        Ok(ConvectedMonopoleField {
            field: CoherentCompactField {
                pressure_pa: values[0],
                particle_velocity_m_s: [values[1], values[2], values[3]],
                active_intensity_w_m2: moments.field.active_intensity_w_m2,
                reactive_intensity_w_m2: moments.field.reactive_intensity_w_m2,
                scope: self.scope(),
            },
            mean_perturbation_energy_density_j_m3: moments.mean_perturbation_energy_density_j_m3,
            mean_perturbation_energy_flux_w_m2: moments.mean_perturbation_energy_flux_w_m2,
        })
    }

    /// Source factors Q=B xi, E[xi xi*]=I, in peak m³/s units. Phases and
    /// cross-source interference survive within each explicitly supplied factor.
    /// This is a cycle-average second moment, not an inferred random phase law.
    pub fn evaluate_coherence_modes(
        &self,
        modes: &[&[C64]],
        budget: CoherenceBudget,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<ConvectedMonopoleStatistics, AeroacError> {
        poll(&mut cancelled)?;
        if modes.len() > budget.max_modes
            || modes
                .len()
                .checked_mul(self.source_count())
                .is_none_or(|work| work > budget.max_source_mode_pairs)
        {
            return Err(invalid("convected monopole coherence exceeds work budget"));
        }
        let mut sum = StatisticsSum::default();
        for &mode in modes {
            sum.add(self.amplitudes(mode, &mut cancelled)?)?;
        }
        poll(&mut cancelled)?;
        self.finish(sum)
    }

    /// Propagate one measured interior Welch bin of physical volume flows.
    /// There is exactly one Q channel [m³/s] per original source, in prepared
    /// position order; synchronized recorded phases are preserved. No force
    /// dipole, mass-flow conversion or unmeasured source correlation is added.
    ///
    /// The original uniform-flow p/u paths propagate every factor coherently.
    /// One-sided /Hz factors become peak bin factors through sqrt(2 df), then
    /// the existing statistics and perturbation-energy owner form the results.
    /// Pressure work remains distinct from flow-transported energy. Outputs
    /// integrate one bin; they are neither PSD/Hz nor coherent-gain tone fits.
    /// Welch leakage/resolution remain, and no cross-frequency coherence,
    /// segment independence or statistical confidence is inferred.
    ///
    /// Allocation-free O(sources * segments), with bounded cancellation and
    /// no repeated propagation kernels or dense covariance construction.
    pub fn evaluate_welch_bin(
        &self,
        spectrum: &fs_fft::welch::WelchCrossSpectrum,
        bin: usize,
        budget: CoherenceBudget,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<ConvectedMonopoleStatistics, AeroacError> {
        poll(&mut cancelled)?;
        let measured = super::welch::admit(
            spectrum,
            bin,
            self.frequency_hz,
            self.source_count(),
            1,
            1,
            budget,
        )?;
        let mut sum = StatisticsSum::default();
        for factor in measured.factors() {
            poll(&mut cancelled)?;
            let mut sums = [ComplexSum::default(); 4];
            for (path, flow) in self.paths.iter().zip(factor) {
                poll(&mut cancelled)?;
                let peak = measured
                    .peak_factor(*flow)
                    .map_err(|_| invalid("recorded volume-flow peak factor is unrepresentable"))?;
                path.accumulate(C64::new(peak.re, peak.im), &mut sums)?;
            }
            sum.add(sums.map(ComplexSum::value))?;
        }
        poll(&mut cancelled)?;
        self.finish(sum)
    }

    /// Feed original propagated p/u into the existing finite-window observer.
    /// Its reported <p u> remains pressure work; flow-carried energy is separate.
    pub fn evaluate_harmonic_line(
        &self,
        flows: &[C64],
    ) -> Result<crate::observation::HarmonicFieldLine<3>, AeroacError> {
        let [pressure, ux, uy, uz] = self.amplitudes(flows, &mut || false)?;
        Ok(crate::observation::HarmonicFieldLine {
            frequency_hz: self.frequency_hz,
            pressure_pa: pressure,
            particle_velocity_m_s: [ux, uy, uz],
        })
    }

    fn amplitudes(
        &self,
        flows: &[C64],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<[C64; 4], AeroacError> {
        poll(cancelled)?;
        if flows.len() != self.source_count() {
            return Err(invalid("one volume flow is required per convected source"));
        }
        let mut sums = [ComplexSum::default(); 4];
        for (path, &flow) in self.paths.iter().zip(flows) {
            poll(cancelled)?;
            path.accumulate(flow, &mut sums)?;
        }
        let values = sums.map(ComplexSum::value);
        if values.iter().any(|&value| !finite(value)) {
            return Err(nonfinite("accumulated convected monopole field"));
        }
        poll(cancelled)?;
        Ok(values)
    }

    fn finish(&self, sum: StatisticsSum) -> Result<ConvectedMonopoleStatistics, AeroacError> {
        let field = sum.finish(self.scope())?;
        let mut energy = Dd::from_f64(
            (0.5 * field.mean_square_pressure_pa2) / self.acoustic.adiabatic_bulk_modulus_pa(),
        );
        for axis in 0..3 {
            energy = energy
                + Dd::from_f64(0.5 * self.acoustic.density_kg_m3())
                    * Dd::from_f64(field.velocity_covariance_m2_s2[axis][axis].re);
        }
        let mean_perturbation_energy_density_j_m3 = energy.to_f64();
        let mean_perturbation_energy_flux_w_m2 = core::array::from_fn(|axis| {
            (Dd::from_f64(field.active_intensity_w_m2[axis])
                + Dd::from_f64(self.flow_velocity_m_s[axis]) * energy)
                .to_f64()
        });
        if !mean_perturbation_energy_density_j_m3.is_finite()
            || mean_perturbation_energy_flux_w_m2
                .iter()
                .any(|v| !v.is_finite())
        {
            return Err(nonfinite("convected perturbation energy or transport"));
        }
        Ok(ConvectedMonopoleStatistics {
            field,
            mean_perturbation_energy_density_j_m3,
            mean_perturbation_energy_flux_w_m2,
        })
    }
}

impl Path {
    fn accumulate(&self, flow: C64, sums: &mut [ComplexSum; 4]) -> Result<(), AeroacError> {
        if !finite(flow) {
            return Err(nonfinite("convected volume flow"));
        }
        for (sum, &transfer) in sums.iter_mut().zip(&self.field) {
            sum.add(transfer * flow)?;
        }
        Ok(())
    }
}

fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), AeroacError> {
    if cancelled() {
        Err(AeroacError::Cancelled {
            what: "convected monopole observation",
        })
    } else {
        Ok(())
    }
}
fn norm(value: [f64; 3]) -> f64 {
    det::hypot(det::hypot(value[0], value[1]), value[2])
}
fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a.into_iter()
        .zip(b)
        .fold(Dd::ZERO, |sum, (x, y)| {
            sum + Dd::from_f64(x) * Dd::from_f64(y)
        })
        .to_f64()
}

#[cfg(test)]
#[path = "convected/tests.rs"]
mod tests;
