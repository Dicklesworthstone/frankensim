//! Unsteady developed pipes between the same two pressure reservoirs.
//!
//! Every branch retains its original [`LaminarPipe`] equation, source and
//! viscous-profile memory. Their common pressure difference drives independent
//! exact trajectories, so `Q_total(t) = sum(Q_i(t))` at every time, rather than
//! splitting flow by steady resistances during startup or reversal. Pressure
//! work, kinetic storage and irreversible heat add over those same trajectories.
//!
//! This is a parallel connection of the original developed-flow models, not a
//! mixing-junction, entrance or manifold-flow model. All branches share one
//! fixed reference density. Each branch retains its own frequency band and
//! Reynolds admission. Viscosity changes use the original state-preserving
//! material continuation; temperature evolution, density changes and reservoir
//! storage remain with their physical owners.

use std::sync::Arc;

use fs_exec::Cx;
use fs_material::fluid_transport::FluidTransport;
use fs_math::c64::C64;
use fs_phs::PreparedStepRecord;
use fs_qty::{Pressure, VolumetricFlowRate};

use super::{
    LaminarPipe, LaminarPipeError, LaminarPipeState, LaminarPipeStep, PreparedLaminarPipe,
    checkpoint, positive,
};

/// Original source-bearing pipes, all oriented between the same two reservoirs.
pub struct LaminarParallel<M> {
    branches: Vec<LaminarPipe<M>>,
    identity: Arc<()>,
}

/// Independent viscous histories with exactly one accepted physical clock.
/// Cloning produces an independent proposed state without resampling a source.
#[derive(Clone, Debug)]
pub struct LaminarParallelState {
    branches: Vec<LaminarPipeState>,
    identity: Arc<()>,
}

impl LaminarParallelState {
    pub fn accepted_time_s(&self) -> f64 {
        self.branches[0].accepted_time_s()
    }
    /// Branch order is the original authored order, including reversed flows.
    pub fn branches(&self) -> &[LaminarPipeState] {
        &self.branches
    }
    pub fn volume_flow_m3_s(&self) -> f64 {
        self.branches
            .iter()
            .map(LaminarPipeState::volume_flow_m3_s)
            .sum()
    }
    pub fn mean_kinetic_energy_j(&self) -> f64 {
        self.branches
            .iter()
            .map(LaminarPipeState::mean_kinetic_energy_j)
            .sum()
    }
    pub fn profile_memory_energy_j(&self) -> f64 {
        self.branches
            .iter()
            .map(LaminarPipeState::profile_memory_energy_j)
            .sum()
    }
    pub fn kinetic_energy_j(&self) -> f64 {
        self.branches
            .iter()
            .map(LaminarPipeState::kinetic_energy_j)
            .sum()
    }
}

/// Whole-bundle sums of the original branch flow and energy records.
pub type LaminarParallelStep = LaminarPipeStep;

impl<M: FluidTransport> LaminarParallel<M> {
    /// Connect already admitted pipes without resampling their material sources.
    pub fn try_new(cx: &Cx<'_>, branches: Vec<LaminarPipe<M>>) -> Result<Self, LaminarPipeError> {
        checkpoint(cx)?;
        let Some(first) = branches.first() else {
            return Err(LaminarPipeError::Invalid(
                "a laminar parallel bundle needs a branch",
            ));
        };
        let density = first.transport().density_kg_m3();
        for pipe in &branches {
            checkpoint(cx)?;
            if pipe.transport().density_kg_m3() != density {
                return Err(LaminarPipeError::Invalid(
                    "parallel branches require exactly common reference density",
                ));
            }
        }
        Ok(Self {
            branches,
            identity: Arc::new(()),
        })
    }

    pub fn branches(&self) -> &[LaminarPipe<M>] {
        &self.branches
    }

    /// Continue every branch under its own new viscosity, retaining all accepted
    /// profiles and their clock. Geometry, density and band admission remain the
    /// responsibility of the original pipe owner; no branch is reset or dropped.
    pub fn with_sources<N: FluidTransport>(
        &self,
        cx: &Cx<'_>,
        sources: Vec<N>,
    ) -> Result<LaminarParallel<N>, LaminarPipeError> {
        checkpoint(cx)?;
        if sources.len() != self.branches.len() {
            return Err(LaminarPipeError::Invalid(
                "one source per original parallel branch required",
            ));
        }
        let branches = self
            .branches
            .iter()
            .zip(sources)
            .map(|(pipe, source)| pipe.with_source(cx, source))
            .collect::<Result<Vec<_>, _>>()?;
        let mut next = LaminarParallel::try_new(cx, branches)?;
        next.identity = Arc::clone(&self.identity);
        Ok(next)
    }

    /// Original harmonic impedances combine as `1/Z = sum(1/Z_i)`.
    /// Every branch must admit this frequency; no branch is bypassed or replaced
    /// by its steady resistance outside its sampled physical frequency band.
    pub fn impedance_at(&self, omega_rad_s: f64) -> Result<C64, LaminarPipeError> {
        let mut admittance = C64::ZERO;
        for pipe in &self.branches {
            admittance = admittance + pipe.impedance_at(omega_rad_s)?.recip();
        }
        let impedance = admittance.recip();
        if !positive(admittance.re)
            || !positive(-admittance.im)
            || !positive(impedance.re)
            || !positive(impedance.im)
        {
            return Err(LaminarPipeError::Invalid(
                "unrepresentable parallel impedance",
            ));
        }
        Ok(impedance)
    }

    /// Initialize each physical Poiseuille profile at the same pressure drop.
    /// Zero pressure initializes rest. A transient generally has a different
    /// split and must retain its original branch states instead of calling this.
    pub fn steady_state(
        &self,
        pressure_drop: Pressure,
        accepted_time_s: f64,
    ) -> Result<LaminarParallelState, LaminarPipeError> {
        let dp = pressure_drop.value();
        if !dp.is_finite() {
            return Err(LaminarPipeError::Invalid(
                "finite parallel pressure difference required",
            ));
        }
        let branches = self
            .branches
            .iter()
            .map(|pipe| {
                let flow = dp / pipe.steady_resistance_pa_s_m3();
                if !flow.is_finite() || (dp != 0.0 && flow == 0.0) {
                    return Err(LaminarPipeError::Invalid(
                        "unrepresentable steady branch flow",
                    ));
                }
                pipe.steady_state(VolumetricFlowRate::new(flow), accepted_time_s)
            })
            .collect::<Result<Vec<_>, _>>()?;
        self.state_from_branches(branches)
    }

    /// Assemble existing original histories without destroying any branch's
    /// nonsteady profile. All histories must belong to their respective pipes
    /// and have the same accepted clock. They need not have the same flow sign.
    pub fn state_from_branches(
        &self,
        branches: Vec<LaminarPipeState>,
    ) -> Result<LaminarParallelState, LaminarPipeError> {
        if branches.len() != self.branches.len() {
            return Err(LaminarPipeError::Invalid(
                "one history per original parallel branch required",
            ));
        }
        let time = branches[0].accepted_time_s();
        if branches.iter().zip(&self.branches).any(|(state, pipe)| {
            !Arc::ptr_eq(&state.identity, &pipe.identity) || state.accepted_time_s() != time
        }) {
            return Err(LaminarPipeError::Invalid(
                "original branch histories at one common clock required",
            ));
        }
        let state = LaminarParallelState {
            branches,
            identity: Arc::clone(&self.identity),
        };
        if !state.volume_flow_m3_s().is_finite() || !state.kinetic_energy_j().is_finite() {
            return Err(LaminarPipeError::Invalid(
                "unrepresentable parallel flow or storage",
            ));
        }
        Ok(state)
    }

    /// Prepare each original exact propagator independently. This avoids a dense
    /// exponential on the combined bundle and retains per-branch irreversible
    /// heat directly. Cold preparation allocates; warm advancement reuses buffers.
    pub fn prepare(
        &self,
        cx: &Cx<'_>,
        duration_s: f64,
    ) -> Result<PreparedLaminarParallel<'_, M>, LaminarPipeError> {
        checkpoint(cx)?;
        let prepared = self
            .branches
            .iter()
            .map(|pipe| pipe.prepare(cx, duration_s))
            .collect::<Result<Vec<_>, _>>()?;
        let pending = self
            .branches
            .iter()
            .map(|pipe| pipe.steady_state(VolumetricFlowRate::new(0.0), 0.0))
            .collect::<Result<Vec<_>, _>>()?;
        checkpoint(cx)?;
        Ok(PreparedLaminarParallel {
            bundle: self,
            prepared,
            pending,
            records: Vec::with_capacity(self.branches.len()),
            pending_records: Vec::with_capacity(self.branches.len()),
        })
    }

    fn admit_owner(&self, state: &LaminarParallelState) -> Result<(), LaminarPipeError> {
        if !Arc::ptr_eq(&state.identity, &self.identity) {
            return Err(LaminarPipeError::Invalid(
                "viscous histories belong to a different parallel bundle",
            ));
        }
        Ok(())
    }
}

/// Original branch propagators with reusable transactional state and report buffers.
pub struct PreparedLaminarParallel<'a, M> {
    bundle: &'a LaminarParallel<M>,
    prepared: Vec<PreparedLaminarPipe<'a, M>>,
    pending: Vec<LaminarPipeState>,
    records: Vec<LaminarPipeStep>,
    pending_records: Vec<LaminarPipeStep>,
}

impl<M: FluidTransport> PreparedLaminarParallel<'_, M> {
    pub fn duration_s(&self) -> f64 {
        self.prepared[0].duration_s()
    }

    /// Records from the last SUCCESSFUL bundle step in authored branch order.
    /// Empty before the first success; a refused step preserves these records.
    /// Each record retains its original native numerical admission and actual
    /// irreversible heat, which a thermal consumer can deposit in that branch.
    pub fn branch_steps(&self) -> &[LaminarPipeStep] {
        &self.records
    }

    /// Actual integrated relation `Q_total_bar = history + mobility * dp`.
    /// Histories are summed with their signs; opposing branch flows are retained.
    pub fn average_flow_relation(
        &self,
        state: &LaminarParallelState,
    ) -> Result<(f64, f64), LaminarPipeError> {
        self.bundle.admit_owner(state)?;
        let mut history = 0.0;
        let mut mobility = 0.0;
        for (prepared, branch) in self.prepared.iter().zip(&state.branches) {
            let relation = prepared.average_flow_relation(branch)?;
            history += relation.0;
            mobility += relation.1;
        }
        if !history.is_finite() || !positive(mobility) {
            return Err(LaminarPipeError::Invalid(
                "unrepresentable parallel flow relation",
            ));
        }
        Ok((history, mobility))
    }

    /// Find a HELD pressure difference for a prescribed interval-average total
    /// flow, using the original condensed trajectories. This does not prescribe
    /// instantaneous total flow or enforce a sealed, zero-storage junction.
    /// `advance` still admits every branch's complete Reynolds envelope.
    pub fn pressure_for_average_flow(
        &self,
        state: &LaminarParallelState,
        total_flow: VolumetricFlowRate,
    ) -> Result<Pressure, LaminarPipeError> {
        let (history, mobility) = self.average_flow_relation(state)?;
        let increment = total_flow.value() - history;
        let pressure = increment / mobility;
        if !total_flow.value().is_finite()
            || !pressure.is_finite()
            || (increment != 0.0 && pressure == 0.0)
        {
            return Err(LaminarPipeError::Invalid(
                "unrepresentable held pressure for total mean flow",
            ));
        }
        Ok(Pressure::new(pressure))
    }

    /// Advance all branches between the SAME held pressures. Their exact
    /// trajectories define total flow at both reservoirs; no steady flow split
    /// or endpoint-average approximation enters the connection.
    ///
    /// The returned native record sums actual work, storage change and heat.
    /// Its solver certificate is intentionally absent: original independent
    /// solver admissions remain on `branch_steps`, not a fabricated joint solve.
    /// The Reynolds envelope is the maximum across branches, each of which
    /// checks its OWN ceiling. Refusal preserves every accepted state and report.
    pub fn advance(
        &mut self,
        cx: &Cx<'_>,
        state: &mut LaminarParallelState,
        pressures: [Pressure; 2],
    ) -> Result<LaminarParallelStep, LaminarPipeError> {
        checkpoint(cx)?;
        self.bundle.admit_owner(state)?;
        self.pending_records.clear();
        for ((prepared, pending), accepted) in self
            .prepared
            .iter_mut()
            .zip(&mut self.pending)
            .zip(&state.branches)
        {
            pending.coordinates.copy_from_slice(&accepted.coordinates);
            pending.accepted_time_s = accepted.accepted_time_s;
            self.pending_records
                .push(prepared.advance(cx, pending, pressures)?);
        }
        let mut sum = LaminarPipeStep {
            duration_s: self.duration_s(),
            initial_volume_flow_m3_s: 0.0,
            final_volume_flow_m3_s: 0.0,
            average_volume_flow_m3_s: 0.0,
            displaced_volume_m3: 0.0,
            reynolds_envelope: 0.0,
            initial_mean_kinetic_energy_j: 0.0,
            final_mean_kinetic_energy_j: 0.0,
            initial_profile_memory_energy_j: 0.0,
            final_profile_memory_energy_j: 0.0,
            native: PreparedStepRecord {
                delta_h: 0.0,
                dissipated: 0.0,
                supplied: 0.0,
                newton_iters: 0,
                solver_residual: 0.0,
                solver_certificate: None,
            },
        };
        for (index, branch) in self.pending_records.iter().enumerate() {
            if index % 64 == 0 {
                checkpoint(cx)?;
            }
            sum.initial_volume_flow_m3_s += branch.initial_volume_flow_m3_s;
            sum.final_volume_flow_m3_s += branch.final_volume_flow_m3_s;
            sum.average_volume_flow_m3_s += branch.average_volume_flow_m3_s;
            sum.displaced_volume_m3 += branch.displaced_volume_m3;
            sum.reynolds_envelope = sum.reynolds_envelope.max(branch.reynolds_envelope);
            sum.initial_mean_kinetic_energy_j += branch.initial_mean_kinetic_energy_j;
            sum.final_mean_kinetic_energy_j += branch.final_mean_kinetic_energy_j;
            sum.initial_profile_memory_energy_j += branch.initial_profile_memory_energy_j;
            sum.final_profile_memory_energy_j += branch.final_profile_memory_energy_j;
            sum.native.delta_h += branch.native.delta_h;
            sum.native.dissipated += branch.native.dissipated;
            sum.native.supplied += branch.native.supplied;
            sum.native.newton_iters = sum
                .native
                .newton_iters
                .checked_add(branch.native.newton_iters)
                .ok_or(LaminarPipeError::Invalid(
                    "parallel native iteration count overflow",
                ))?;
            sum.native.solver_residual = sum
                .native
                .solver_residual
                .max(branch.native.solver_residual);
        }
        if ![
            sum.initial_volume_flow_m3_s,
            sum.final_volume_flow_m3_s,
            sum.average_volume_flow_m3_s,
            sum.displaced_volume_m3,
            sum.initial_mean_kinetic_energy_j,
            sum.final_mean_kinetic_energy_j,
            sum.initial_profile_memory_energy_j,
            sum.final_profile_memory_energy_j,
            sum.native.delta_h,
            sum.native.dissipated,
            sum.native.supplied,
            sum.final_mean_kinetic_energy_j + sum.final_profile_memory_energy_j,
        ]
        .into_iter()
        .all(f64::is_finite)
        {
            return Err(LaminarPipeError::Invalid(
                "unrepresentable parallel flow, work or storage",
            ));
        }
        checkpoint(cx)?;
        core::mem::swap(&mut state.branches, &mut self.pending);
        core::mem::swap(&mut self.records, &mut self.pending_records);
        Ok(sum)
    }
}

#[cfg(test)]
mod tests;
