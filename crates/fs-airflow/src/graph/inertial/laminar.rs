//! Circular-pipe Stokes diffusion on the original passive-system owner.
//!
//! The incompressible, axially developed Newtonian equation gives
//! `Z(s) = L0 [s + 8/tau + sum(g_j s/(s+p_j))]`, with
//! `L0=rho length/area` and `tau=rho radius²/mu`. The original positive
//! [`ZkDiffusiveBasis`] supplies the viscous poles and their sampled band
//! admission. No acoustic compliance or thermal shunt is introduced.
//! See Womersley (1955), <https://doi.org/10.1113/jphysiol.1955.sp005276>,
//! and the continued fraction <https://dlmf.nist.gov/10.10.E1>.
//!
//! Mean flow Q and memory flows z_j obey `z'_j=p_j(Q-z_j)`. Their storage is
//! `H=L0 Q²/2 + sum(L0 g_j z_j²/(2 p_j))`; irreversible power is
//! `L0 [8 Q²/tau + sum(g_j (Q-z_j)²)]`. Memory pressure times mean flow is
//! therefore NOT all heat. At steady flow z_j=Q, this storage recovers the
//! circular parabolic profile's kinetic factor 4/3, at every pole count.
//!
//! Walls are stationary, rigid and no-slip. Radius, density and viscosity stay
//! fixed over this owner's history. Entrance flow, profile-dependent convection,
//! compressibility, cavitation and transition prediction remain outside it.
//! Gas use requires caller-qualified low-Mach/frozen-density support: a transport
//! snapshot alone does not supply sound-speed or absolute-pressure admission.
//! The frequency-domain error is sampled over an explicit band; arbitrary
//! startup pulses do not inherit a time-domain error certificate.

use std::sync::Arc;

use fs_exec::Cx;
use fs_material::MaterialError;
use fs_material::fluid_transport::{FluidTransport, FluidTransportState};
use fs_math::{c64::C64, det};
use fs_phs::linear_exact::{LinearSplitBudget, PreparedExactZoh};
use fs_phs::{
    PhsError, PortHamiltonian, PreparedStepRecord, QuadraticStorage, ResistanceWitness,
    ZkDiffusive, ZkDiffusiveBasis,
};
use fs_qty::{Length, Pressure, VolumetricFlowRate};

use super::{GraphError, checkpoint};
use crate::graph::thermal::medium::series::hydraulic::{MaterialFlowResistance, ResistanceError};

/// Geometry and physical approximation choices, all explicit at construction.
#[derive(Clone, Copy, Debug)]
pub struct LaminarPipeSpec {
    pub length: Length,
    pub radius: Length,
    /// Mean-speed Reynolds ceiling in (0,2300]; not a transition prediction.
    pub maximum_reynolds: f64,
    /// Positive, increasing angular-frequency band [rad/s] for wall memory.
    pub band_rad_s: [f64; 2],
    /// Sampled relative complex AND real wall-loss error, strictly in (0,1).
    pub maximum_relative_wall_error: f64,
}

#[derive(Debug)]
pub enum LaminarPipeError {
    Invalid(&'static str),
    Material(MaterialError),
    Resistance(ResistanceError),
    Native(PhsError),
    Execution(GraphError),
}
impl core::fmt::Display for LaminarPipeError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "laminar pipe: {self:?}")
    }
}
impl std::error::Error for LaminarPipeError {}
impl From<MaterialError> for LaminarPipeError {
    fn from(value: MaterialError) -> Self {
        Self::Material(value)
    }
}
impl From<ResistanceError> for LaminarPipeError {
    fn from(value: ResistanceError) -> Self {
        Self::Resistance(value)
    }
}
impl From<PhsError> for LaminarPipeError {
    fn from(value: PhsError) -> Self {
        Self::Native(value)
    }
}
impl From<GraphError> for LaminarPipeError {
    fn from(value: GraphError) -> Self {
        Self::Execution(value)
    }
}

/// Original material owner, frozen physical coefficients, and native dynamics.
/// A borrowed original source is supported through `FluidTransport for &T`.
pub struct LaminarPipe<M> {
    source: M,
    transport: FluidTransportState,
    law: MaterialFlowResistance,
    memory: ZkDiffusive,
    inertance: f64,
    steady_resistance: f64,
    flow_scales: Vec<f64>,
    system: PortHamiltonian,
    identity: Arc<()>,
}

/// Complete private viscous history. Cloning produces an independent proposal.
/// Coordinates are energy-normalized flows, not samples of radial velocity.
#[derive(Clone, Debug)]
pub struct LaminarPipeState {
    coordinates: Vec<f64>,
    mean_flow_scale: f64,
    accepted_time_s: f64,
    identity: Arc<()>,
}
impl LaminarPipeState {
    pub const fn accepted_time_s(&self) -> f64 {
        self.accepted_time_s
    }
    pub fn volume_flow_m3_s(&self) -> f64 {
        self.coordinates[0] / self.mean_flow_scale
    }
    pub fn mean_kinetic_energy_j(&self) -> f64 {
        0.5 * self.coordinates[0].powi(2)
    }
    pub fn profile_memory_energy_j(&self) -> f64 {
        energy(&self.coordinates[1..])
    }
    pub fn kinetic_energy_j(&self) -> f64 {
        energy(&self.coordinates)
    }
    /// Original native state order for explicit composition through fs-phs.
    pub fn coordinates(&self) -> &[f64] {
        &self.coordinates
    }
}

/// A held-pressure interval, retaining actual integrated flow and native work.
#[derive(Clone, Copy, Debug)]
pub struct LaminarPipeStep {
    pub duration_s: f64,
    pub initial_volume_flow_m3_s: f64,
    pub final_volume_flow_m3_s: f64,
    /// Time integral / duration, not endpoint flow or endpoint averaging.
    pub average_volume_flow_m3_s: f64,
    pub displaced_volume_m3: f64,
    pub reynolds_envelope: f64,
    pub initial_mean_kinetic_energy_j: f64,
    pub final_mean_kinetic_energy_j: f64,
    pub initial_profile_memory_energy_j: f64,
    pub final_profile_memory_energy_j: f64,
    /// `dissipated` alone is irreversible heat. `delta_h` includes memory.
    /// Exact ZOH independently integrates pressure work and resistance loss.
    pub native: PreparedStepRecord,
}

impl<M: FluidTransport> LaminarPipe<M> {
    /// Bind one actual source snapshot to the canonical circular wall spectrum.
    /// Reuse a basis across many pipes; no eigenproblem enters warm stepping.
    pub fn try_new(
        cx: &Cx<'_>,
        source: M,
        spec: LaminarPipeSpec,
        basis: &ZkDiffusiveBasis,
    ) -> Result<Self, LaminarPipeError> {
        checkpoint(cx)?;
        let law = MaterialFlowResistance::circular_laminar(
            spec.length,
            Length::new(2.0 * spec.radius.value()),
            spec.maximum_reynolds,
        )?;
        let transport = source.transport_state()?;
        let rho = transport.density_kg_m3();
        let mu = transport.dynamic_viscosity_pa_s();
        // Unit gamma/Pr select ONLY the existing viscous branch: these are
        // not inferred caloric properties of the retained fluid. That branch
        // is independent of compressibility, Cp, conductivity and Prandtl.
        let memory = basis.bind(
            spec.radius.value(),
            rho,
            mu,
            1.0,
            1.0,
            spec.band_rad_s,
            spec.maximum_relative_wall_error,
        )?;
        let inertance = rho * (spec.length.value() / law.flow_area().value());
        let steady_resistance = inertance * memory.series_drag_per_s();
        let scale = det::sqrt(inertance);
        if ![inertance, steady_resistance, scale, 1.0 / scale]
            .into_iter()
            .all(positive)
        {
            return Err(LaminarPipeError::Invalid(
                "unrepresentable pipe inertia or resistance",
            ));
        }
        let n = 1 + memory.series_terms().len();
        let mut flow_scales = vec![scale];
        let mut resistance = vec![0.0; n * n];
        let mut squares = vec![0.0; n * n];
        let mut q = vec![0.0; n * n];
        let mut g = vec![0.0; 2 * n];
        resistance[0] = memory.series_drag_per_s();
        squares[0] = det::sqrt(resistance[0]);
        g[0] = 1.0 / scale;
        g[1] = -g[0];
        q[0] = 1.0;
        for (j, term) in memory.series_terms().iter().enumerate() {
            let k = j + 1;
            let a = det::sqrt(term.gain_per_s);
            let b = det::sqrt(term.rate_per_s);
            let coupling = a * b;
            let arm_scale = scale * det::sqrt(term.gain_per_s / term.rate_per_s);
            if ![a, b, coupling, arm_scale, 1.0 / arm_scale]
                .into_iter()
                .all(positive)
            {
                return Err(LaminarPipeError::Invalid(
                    "unrepresentable viscous memory scale",
                ));
            }
            flow_scales.push(arm_scale);
            resistance[0] += term.gain_per_s;
            resistance[k] = -coupling;
            resistance[k * n] = -coupling;
            resistance[k * n + k] = term.rate_per_s;
            squares[k * n] = a;
            squares[k * n + k] = -b;
            q[k * n + k] = 1.0;
        }
        let witness = ResistanceWitness::try_from_source_rows(n, resistance.clone(), squares)?;
        let system = PortHamiltonian::new(
            n,
            2,
            vec![0.0; n * n],
            resistance,
            g,
            Box::new(QuadraticStorage::new(q, n)?),
        )?
        .with_resistance_witness(witness)?;
        checkpoint(cx)?;
        Ok(Self {
            source,
            transport,
            law,
            memory,
            inertance,
            steady_resistance,
            flow_scales,
            system,
            identity: Arc::new(()),
        })
    }

    pub const fn source(&self) -> &M {
        &self.source
    }
    pub const fn transport(&self) -> FluidTransportState {
        self.transport
    }
    pub const fn steady_law(&self) -> MaterialFlowResistance {
        self.law
    }
    pub fn wall_memory(&self) -> &ZkDiffusive {
        &self.memory
    }
    pub const fn mean_inertance_pa_s2_m3(&self) -> f64 {
        self.inertance
    }
    pub const fn steady_resistance_pa_s_m3(&self) -> f64 {
        self.steady_resistance
    }
    /// Native pressure inputs `[p_from,p_to]` have volume-flow outputs `[Q,-Q]`.
    /// The native harmonic and coupled-system owners consume these same matrices.
    pub fn system(&self) -> &PortHamiltonian {
        &self.system
    }

    /// Transfer the system and ORIGINAL material owner to explicit PHS assembly.
    /// Before consuming, retain needed geometry/band diagnostics. That assembly
    /// owns the resulting history, source and Reynolds/forcing admission; the
    /// standalone pipe's prepared-state checks cannot qualify a different solve.
    pub fn into_native(self) -> (M, PortHamiltonian) {
        (self.source, self.system)
    }

    /// Full finite-pole impedance in the original `exp(+i omega t)` convention.
    /// Outside the sampled physical band the original wall owner refuses.
    pub fn impedance_at(&self, omega_rad_s: f64) -> Result<C64, LaminarPipeError> {
        let value = self
            .memory
            .response_at(omega_rad_s)?
            .0
            .scale(self.inertance);
        if !positive(value.re) || !positive(value.im) {
            return Err(LaminarPipeError::Invalid("unrepresentable pipe impedance"));
        }
        Ok(value)
    }

    /// Initialize a fully developed steady profile, including all viscous memory.
    /// Zero flow is rest; nonzero mean flow with zero memory is NOT Poiseuille.
    pub fn steady_state(
        &self,
        volume_flow: VolumetricFlowRate,
        accepted_time_s: f64,
    ) -> Result<LaminarPipeState, LaminarPipeError> {
        if !accepted_time_s.is_finite() || accepted_time_s < 0.0 {
            return Err(LaminarPipeError::Invalid(
                "finite nonnegative accepted time required",
            ));
        }
        self.admit_flow(volume_flow.value().abs())?;
        let coordinates: Vec<_> = self
            .flow_scales
            .iter()
            .map(|s| s * volume_flow.value())
            .collect();
        if coordinates.iter().any(|x| !x.is_finite())
            || !energy(&coordinates).is_finite()
            || (volume_flow.value() != 0.0
                && (coordinates.iter().any(|x| *x == 0.0) || energy(&coordinates) == 0.0))
        {
            return Err(LaminarPipeError::Invalid(
                "unrepresentable initial pipe state",
            ));
        }
        Ok(LaminarPipeState {
            coordinates,
            mean_flow_scale: self.flow_scales[0],
            accepted_time_s,
            identity: Arc::clone(&self.identity),
        })
    }

    /// Prepare the original exact held-input propagator for a fixed interval.
    /// Dense exponentials and heat integrals are cold work at the bounded pole
    /// count. Warm `advance` and `average_flow_relation` allocate no memory.
    pub fn prepare(
        &self,
        cx: &Cx<'_>,
        duration_s: f64,
    ) -> Result<PreparedLaminarPipe<'_, M>, LaminarPipeError> {
        checkpoint(cx)?;
        if !positive(duration_s) {
            return Err(LaminarPipeError::Invalid(
                "positive finite interval required",
            ));
        }
        let native = PreparedExactZoh::new(
            &self.system,
            duration_s,
            LinearSplitBudget {
                max_modes: self.flow_scales.len(),
                max_omega_dt: 1.0,
            },
        )?;
        let pending = vec![0.0; self.flow_scales.len()];
        checkpoint(cx)?;
        Ok(PreparedLaminarPipe {
            pipe: self,
            native,
            pending,
        })
    }

    fn admit_flow(&self, maximum: f64) -> Result<f64, LaminarPipeError> {
        let speed = maximum / self.law.flow_area().value();
        if maximum > 0.0 && speed == 0.0 {
            return Err(LaminarPipeError::Invalid("unrepresentable mean speed"));
        }
        let reynolds = self
            .transport
            .reynolds(speed, self.law.reynolds_length().value())?;
        if reynolds > self.law.reynolds_range().1 {
            return Err(ResistanceError::Reynolds {
                value: reynolds,
                admitted: self.law.reynolds_range(),
            }
            .into());
        }
        Ok(reynolds)
    }
}

/// One reusable native trajectory/heat cache. No state is changed on refusal.
pub struct PreparedLaminarPipe<'a, M> {
    pipe: &'a LaminarPipe<M>,
    native: PreparedExactZoh<'a>,
    pending: Vec<f64>,
}
impl<M: FluidTransport> PreparedLaminarPipe<'_, M> {
    pub fn duration_s(&self) -> f64 {
        self.native.duration()
    }

    /// `Q_bar = history + mobility * (p_from-p_to)` for THIS state and clock.
    /// This is the actual native time-integrated port relation for a coupled
    /// pressure solve, not a quasi-steady resistance or endpoint approximation.
    /// Final pressure/flow applicability is admitted by `advance`.
    pub fn average_flow_relation(
        &self,
        state: &LaminarPipeState,
    ) -> Result<(f64, f64), LaminarPipeError> {
        self.admit_owner(state)?;
        let (history, mobility) =
            self.native
                .average_relation(&state.coordinates, &[0.0, 0.0], 0)?;
        if !positive(mobility) {
            return Err(LaminarPipeError::Invalid(
                "unrepresentable pressure-flow mobility",
            ));
        }
        Ok((history, mobility))
    }

    /// Evolve both pressure ports held over the prepared interval, atomically.
    /// A comparison bound admits every interior mean flow, including reversals:
    /// the cooperative Q/z equations preserve the extrema of initial flows and
    /// the eventual steady flow; `initial_max + dt*|dp|/L0` is a second bound.
    /// Their minimum is conservative, not a measured temporal error estimate.
    pub fn advance(
        &mut self,
        cx: &Cx<'_>,
        state: &mut LaminarPipeState,
        pressures: [Pressure; 2],
    ) -> Result<LaminarPipeStep, LaminarPipeError> {
        checkpoint(cx)?;
        self.admit_owner(state)?;
        let duration_s = self.native.duration();
        let next_time = state.accepted_time_s + duration_s;
        let physical_pressures = pressures.map(|p| p.value());
        let dp = physical_pressures[0] - physical_pressures[1];
        if physical_pressures.iter().any(|p| !p.is_finite())
            || !dp.is_finite()
            || !next_time.is_finite()
            || next_time <= state.accepted_time_s
        {
            return Err(LaminarPipeError::Invalid(
                "finite pressures and representable clock required",
            ));
        }
        // Common pressure is an exact null direction of an incompressible
        // two-port. Remove it before dense maps to avoid catastrophic common-
        // mode cancellation; the total work remains (p_from-p_to)*integral Q.
        let input = [dp, 0.0];
        let initial_max = state
            .coordinates
            .iter()
            .zip(&self.pipe.flow_scales)
            .map(|(x, scale)| (x / scale).abs())
            .fold(0.0_f64, f64::max);
        let steady_max = dp.abs() / self.pipe.steady_resistance;
        let forced_increment = (dp.abs() / self.pipe.inertance) * duration_s;
        let maximum = initial_max
            .max(steady_max)
            .min(initial_max + forced_increment);
        if !maximum.is_finite() || (initial_max == 0.0 && dp != 0.0 && maximum == 0.0) {
            return Err(LaminarPipeError::Invalid("unrepresentable flow envelope"));
        }
        let reynolds_envelope = self.pipe.admit_flow(maximum)?;
        let initial_volume_flow_m3_s = state.volume_flow_m3_s();
        let initial_mean_kinetic_energy_j = state.mean_kinetic_energy_j();
        let initial_profile_memory_energy_j = state.profile_memory_energy_j();
        let mut output = [0.0; 2];
        let native =
            self.native
                .step_into(&state.coordinates, &input, &mut self.pending, &mut output)?;
        let final_volume_flow_m3_s = self.pending[0] / self.pipe.flow_scales[0];
        let final_mean_kinetic_energy_j = 0.5 * self.pending[0].powi(2);
        let final_profile_memory_energy_j = energy(&self.pending[1..]);
        let displaced_volume_m3 = duration_s * output[0];
        if ![
            final_volume_flow_m3_s,
            final_mean_kinetic_energy_j,
            final_profile_memory_energy_j,
            displaced_volume_m3,
        ]
        .into_iter()
        .all(f64::is_finite)
            || (output[0] != 0.0 && displaced_volume_m3 == 0.0)
            || output[1] != -output[0]
        {
            return Err(LaminarPipeError::Invalid(
                "unrepresentable flow, energy or port continuity",
            ));
        }
        checkpoint(cx)?;
        state.coordinates.copy_from_slice(&self.pending);
        state.accepted_time_s = next_time;
        Ok(LaminarPipeStep {
            duration_s,
            initial_volume_flow_m3_s,
            final_volume_flow_m3_s,
            average_volume_flow_m3_s: output[0],
            displaced_volume_m3,
            reynolds_envelope,
            initial_mean_kinetic_energy_j,
            final_mean_kinetic_energy_j,
            initial_profile_memory_energy_j,
            final_profile_memory_energy_j,
            native,
        })
    }

    fn admit_owner(&self, state: &LaminarPipeState) -> Result<(), LaminarPipeError> {
        if !Arc::ptr_eq(&state.identity, &self.pipe.identity) {
            return Err(LaminarPipeError::Invalid(
                "viscous history belongs to a different pipe",
            ));
        }
        Ok(())
    }
}

fn positive(value: f64) -> bool {
    value.is_finite() && value > 0.0
}
fn energy(values: &[f64]) -> f64 {
    values.iter().map(|x| 0.5 * x * x).sum()
}

#[cfg(test)]
mod tests;
