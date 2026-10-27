//! Low-induction material MHD on the existing electrical, MAC and enthalpy owners.
//!
//! The native moving-conductor current operator supplies BOTH I²R and its
//! reciprocal Lorentz force. A bounded implicit iteration evaluates that force
//! at the same endpoint velocity as native backward-Euler mechanical work.
//! Neither transfer defects nor backward-Euler numerical loss become heat.

use super::*;
use crate::conjugate::enthalpy::{
    JouleEnthalpyError, JouleEnthalpyInput, StepConfig, admit_map as admit_electrical_map,
    checkpoint, electrode_power, finite,
};
use crate::conjugate::{ElectricMotionInput, MovingConductorSolution, UnsteadyHeatingAudit};
use fs_conduction::transient::deposition::EnergyAudit;

/// Electrical topology is an explicit choice, not inferred from a missing drive.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum MhdCurrentBoundary {
    /// Every connected conductor needs the original owner's potential reference.
    Electrodes,
    /// No electrodes; retain circulating currents and check every cell's KCL.
    Isolated { current_tolerance_a: f64 },
}

/// Fixed imposed field and bounded current/momentum iteration for one advance.
#[derive(Clone, Copy, Debug)]
pub struct MhdEnthalpyInput<'a> {
    /// One stationary imposed B [T] per voxel. The caller owns solenoidal-field
    /// validity and the low magnetic Reynolds number approximation.
    pub magnetic_flux_density_t: &'a [[f64; 3]],
    pub current_boundary: MhdCurrentBoundary,
    /// Maximum current/momentum trials per accepted thermal interval. Original
    /// whole-advance electrical AND material-flow solve ceilings also apply.
    pub max_coupling_iterations: usize,
    /// Picard velocity relaxation in (0, 1]. Strong magnetic damping may need a
    /// smaller value. Failure to converge refuses the complete advance.
    pub relaxation: f64,
    /// Maximum componentwise difference between electrical evaluation velocity
    /// and native endpoint velocity [m/s]. Mechanical power is checked separately
    /// against JouleEnthalpyInput::power_tolerance_w, even at small velocities.
    pub velocity_tolerance_m_s: f64,
}

/// The exact source epoch and velocity used by one native current solve.
#[derive(Clone, Debug, PartialEq)]
pub struct MhdCurrentState {
    pub time_s: f64,
    pub temperature_k: Vec<f64>,
    pub specific_enthalpy_j_kg: Vec<f64>,
    /// None denotes insulation or an ideal electrode with no bulk getter.
    pub resistivity_ohm_m: Vec<Option<f64>>,
    pub velocity_m_s: Vec<[f64; 3]>,
    /// Anchor-on-material force [N] per original voxel, equal and opposite to
    /// the Lorentz load in an exact frozen cell. The force on the anchor is
    /// its negative. Mobile cells have zero reaction. This direct magnetic
    /// reaction excludes fluid/solid traction and wall reactions.
    pub frozen_lorentz_reaction_n: Vec<[f64; 3]>,
    pub solution: MovingConductorSolution,
}

/// Accepted physical and numerical accounts, all over actual thermal durations.
#[derive(Clone, Debug, PartialEq)]
pub struct MhdEnthalpyReport {
    /// Contains BOTH deposited Joule and physical viscous/Darcy heat already.
    pub thermal: EnergyAudit,
    pub mechanical: UnsteadyHeatingAudit,
    pub electrode_energy_j: f64,
    pub joule_energy_j: f64,
    pub deposited_joule_j: f64,
    /// Signed force work on the native accepted endpoint velocity; negative
    /// for magnetic braking. This is never an additional heat source.
    pub lorentz_work_j: f64,
    /// Work computed by the reciprocal current owner on its sampled velocity.
    pub sampled_lorentz_work_j: f64,
    /// Accepted anchor-on-material Lorentz reaction impulse [N s], integrated
    /// per original voxel over actual thermal intervals, before endpoint phase
    /// changes. The impulse delivered to the external anchor is its negative.
    pub frozen_lorentz_reaction_impulse_n_s: Vec<[f64; 3]>,
    /// Exactly zero: every supported cell and native constrained face must
    /// have zero velocity. A moving support is outside this explicit model.
    pub frozen_support_work_j: f64,
    pub electrical_defect_j: f64,
    pub heat_transfer_defect_j: f64,
    /// Actual native-endpoint Lorentz work minus electrically sampled work.
    pub mechanical_transfer_defect_j: f64,
    /// Sum of absolute interval transfer defects; opposite signs cannot hide
    /// an unconverged time-centering mismatch.
    pub absolute_mechanical_transfer_defect_j: f64,
    /// Native face body work minus independently lifted cell force work.
    pub native_body_transfer_defect_j: f64,
    /// Original external body force work, excluding Lorentz and buoyancy.
    pub other_body_work_j: f64,
    /// Original thermal input with deposited electrical/mechanical heat removed.
    pub other_thermal_input_j: f64,
    /// Delta thermal + delta kinetic + BE numerical loss minus incoming
    /// enthalpy, original thermal inputs, electrode/wall/body/buoyancy work and
    /// the separately reported native pressure-constraint work. The numerical
    /// loss and pressure work retain their own mechanical accounts.
    pub combined_residual_j: f64,
    pub current_solves: usize,
    pub momentum_trials: usize,
    pub last_applied: Option<MhdCurrentState>,
    pub last_duration_s: Option<f64>,
    pub last_velocity_defect_m_s: Option<f64>,
    /// Fresh material/current support at final h/T and actual native velocity.
    /// This diagnostic solve does not deposit another heat increment.
    pub endpoint: MhdCurrentState,
}

#[derive(Debug)]
pub enum MaterialFlowMhdError<E> {
    Flow(MaterialFlowError<E>),
    Electrical(JouleEnthalpyError),
    CouplingDidNotConverge {
        iterations: usize,
        velocity_defect_m_s: f64,
        mechanical_power_defect_w: f64,
    },
}
impl<E> From<DepositionError> for MaterialFlowMhdError<E> {
    fn from(error: DepositionError) -> Self {
        Self::Electrical(error.into())
    }
}
impl<E> From<JouleEnthalpyError> for MaterialFlowMhdError<E> {
    fn from(error: JouleEnthalpyError) -> Self {
        Self::Electrical(error)
    }
}
impl<E> From<MaterialFlowError<E>> for MaterialFlowMhdError<E> {
    fn from(error: MaterialFlowError<E>) -> Self {
        Self::Flow(error)
    }
}
impl<E: core::fmt::Display> core::fmt::Display for MaterialFlowMhdError<E> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Flow(error) => write!(f, "material MHD flow: {error}"),
            Self::Electrical(error) => write!(f, "material MHD current: {error}"),
            Self::CouplingDidNotConverge {
                iterations,
                velocity_defect_m_s,
                mechanical_power_defect_w,
            } => write!(
                f,
                "material MHD failed after {iterations} trials: velocity defect {velocity_defect_m_s} m/s, mechanical work defect {mechanical_power_defect_w} W"
            ),
        }
    }
}
impl<E: std::error::Error + 'static> std::error::Error for MaterialFlowMhdError<E> {}

type ElectricalResponse<'a, F, T> =
    dyn Fn(&FlowMaterialSample<F, T>) -> Result<Option<f64>, ChtError> + 'a;

#[derive(Clone)]
struct Pending {
    state: MhdCurrentState,
    start_s: f64,
    dt_s: f64,
    end_s: f64,
    mapped_w: f64,
    electrode_w: f64,
    lorentz_w: f64,
    other_body_w: f64,
    body_defect_w: f64,
    velocity_defect_m_s: f64,
}

#[derive(Clone, Default)]
struct Work {
    electrode_j: f64,
    joule_j: f64,
    mapped_j: f64,
    lorentz_j: f64,
    sampled_lorentz_j: f64,
    frozen_lorentz_reaction_impulse_n_s: Vec<[f64; 3]>,
    other_body_j: f64,
    body_defect_j: f64,
    absolute_transfer_defect_j: f64,
    absolute_electrical_defect_j: f64,
    absolute_body_defect_j: f64,
    absolute_heat_defect_j: f64,
}

#[derive(Clone)]
struct MhdFlow<'a, 'f, F, R, O, T> {
    flow: MaterialFlowFeedback<'f, F, R, O, T>,
    input: JouleEnthalpyInput<'a>,
    motion: MhdEnthalpyInput<'a>,
    response: &'a ElectricalResponse<'a, F, T>,
    current_solves: usize,
    momentum_trials: usize,
    pending: Option<Pending>,
    last: Option<Pending>,
    work: Work,
}

fn invalid(reason: &'static str) -> JouleEnthalpyError {
    JouleEnthalpyError::Invalid(reason)
}

#[allow(clippy::too_many_arguments)]
fn evaluate<F, T: FlowTransportSnapshot>(
    input: JouleEnthalpyInput<'_>,
    motion: MhdEnthalpyInput<'_>,
    gate: &CancelGate,
    solves: &mut usize,
    thermal: &DepositionSimulation<'_>,
    samples: &[FlowMaterialSample<F, T>],
    velocity: &[[f64; 3]],
    response: &ElectricalResponse<'_, F, T>,
) -> Result<MhdCurrentState, JouleEnthalpyError> {
    admit_electrical_map(thermal, input, gate)?;
    if *solves >= input.max_current_solves {
        return Err(JouleEnthalpyError::CurrentBudget);
    }
    if samples.len() != input.cell_ids.len() || velocity.len() != samples.len() {
        return Err(invalid("MHD requires every original flow material sample"));
    }
    let mut temperature_k = Vec::with_capacity(samples.len());
    let mut specific_enthalpy_j_kg = Vec::with_capacity(samples.len());
    for (voxel, (sample, &id)) in samples.iter().zip(input.cell_ids).enumerate() {
        checkpoint(gate)?;
        if sample.voxel != voxel
            || sample.cell != id
            || sample.specific_enthalpy_j_kg.to_bits()
                != thermal.cells()[&id].specific_h_j_kg.to_bits()
            || sample.temperature_k.to_bits() != thermal.temperature_k(id)?.to_bits()
        {
            return Err(invalid(
                "MHD current must use the exact receiving material epoch",
            ));
        }
        temperature_k.push(sample.temperature_k);
        specific_enthalpy_j_kg.push(sample.specific_enthalpy_j_kg);
    }
    let mut resistivity_ohm_m = vec![None; samples.len()];
    let resolve = |c: usize, _, _| {
        let value = response(&samples[c])?;
        resistivity_ohm_m[c] = value;
        Ok(value)
    };
    let native = ElectricMotionInput {
        velocity_m_s: velocity,
        magnetic_flux_density_t: motion.magnetic_flux_density_t,
        max_cells: input.max_cells,
        power_tolerance_w: input.power_tolerance_w,
    };
    *solves += 1;
    let solution = match motion.current_boundary {
        MhdCurrentBoundary::Electrodes => input.setup.solve_with_cell_motion(
            input.domain,
            &temperature_k,
            input.contacts,
            native,
            resolve,
            gate,
        )?,
        MhdCurrentBoundary::Isolated {
            current_tolerance_a,
        } => input.setup.solve_isolated_cell_motion(
            input.domain,
            &temperature_k,
            input.contacts,
            native,
            current_tolerance_a,
            resolve,
            gate,
        )?,
    };
    // Electrical cells retain their original material identity through phase
    // changes. Frozen conductors still carry current and dissipate I²R; their
    // Lorentz load is carried by the explicitly stationary external anchor.
    let mut frozen_lorentz_reaction_n = vec![[0.0; 3]; samples.len()];
    for (c, sample) in samples.iter().enumerate() {
        checkpoint(gate)?;
        if matches!(
            sample.transport.phase_mobility(),
            Some(PhaseMobility::ImmobileSolid { .. })
        ) {
            if velocity[c].iter().any(|&value| value != 0.0) {
                return Err(invalid("anchored MHD requires exactly stationary frozen material"));
            }
            for axis in 0..3 {
                frozen_lorentz_reaction_n[c][axis] = finite(
                    -solution.lorentz_force_density_n_m3[c][axis] * input.domain.volume(c),
                )?;
            }
        }
    }
    Ok(MhdCurrentState {
        time_s: thermal.time_s(),
        temperature_k,
        specific_enthalpy_j_kg,
        resistivity_ohm_m,
        velocity_m_s: velocity.to_vec(),
        frozen_lorentz_reaction_n,
        solution,
    })
}

// This is the transpose of the existing MAC cell-to-face force lift: each cell
// contributes V/2 to either normal face and cell u is their arithmetic average.
// It is independent of the native assembled momentum RHS/work reduction.
fn cell_work(
    domain: &VoxelDomain,
    force: &[[f64; 3]],
    velocity: &[[f64; 3]],
    gate: &CancelGate,
) -> Result<f64, JouleEnthalpyError> {
    if force.len() != domain.cell_count() || velocity.len() != force.len() {
        return Err(invalid("MHD force/velocity map changed"));
    }
    let mut power = 0.0;
    for (c, (force, velocity)) in force.iter().zip(velocity).enumerate() {
        checkpoint(gate)?;
        for axis in 0..3 {
            power = finite(power + force[axis] * velocity[axis] * domain.volume(c))?;
        }
    }
    Ok(power)
}

impl<F: Clone, R, O: DepositionObserver, T: FlowTransportSnapshot> DepositionObserver
    for MhdFlow<'_, '_, F, R, O, T>
where
    R: Clone + FnMut(usize, &DepositedCell, f64) -> Result<F, MaterialError>,
{
    type Error = MaterialFlowMhdError<O::Error>;

    fn requires_step_interval(&self) -> bool {
        true
    }

    fn begin(&mut self) -> Result<(), Self::Error> {
        self.flow.begin()?;
        self.current_solves = 0;
        self.momentum_trials = 0;
        self.pending = None;
        self.last = None;
        self.work = Work {
            frozen_lorentz_reaction_impulse_n_s: vec![[0.0; 3]; self.input.cell_ids.len()],
            ..Work::default()
        };
        Ok(())
    }

    fn sources(
        &mut self,
        thermal: &DepositionSimulation<'_>,
        check: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<BTreeMap<CellId, f64>, Self::Error> {
        check()?;
        if self.pending.is_some() {
            return Err(invalid("unaccepted MHD interval").into());
        }
        if self.current_solves >= self.input.max_current_solves {
            return Err(JouleEnthalpyError::CurrentBudget.into());
        }
        // Retained momentum sources only prepare material/advection; its solve
        // happens below. Count each ACTUAL trial, rather than this preparation.
        let sources = self.flow.sources(thermal, check)?;
        self.flow.solves = self.momentum_trials;
        Ok(sources)
    }

    fn step_sources(
        &mut self,
        thermal: &DepositionSimulation<'_>,
        start_s: f64,
        dt_s: f64,
        end_s: f64,
        check: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<BTreeMap<CellId, f64>, Self::Error> {
        check()?;
        if self.pending.is_some() {
            return Err(invalid("MHD interval already prepared").into());
        }
        let epoch = self
            .flow
            .applied()
            .ok_or(invalid("missing MHD source epoch"))?;
        if epoch.time_s.to_bits() != start_s.to_bits()
            || thermal.time_s().to_bits() != start_s.to_bits()
        {
            return Err(invalid("stale MHD source time").into());
        }
        let mut velocity = self
            .flow
            .unsteady_flow()
            .ok_or(invalid("missing native MHD momentum"))?
            .velocity_m_s
            .clone();
        // An all-solid interval has an exact zero-dimensional momentum system.
        // It still solves current and deposits heat, but consumes no native
        // momentum trial, exactly as the original phase owner specifies.
        let native_solve = epoch.samples.iter()
            .any(|sample| sample.transport.mobile_transport().is_some());
        let mut velocity_defect = f64::INFINITY;
        let mut work_defect = f64::INFINITY;
        for _ in 0..self.motion.max_coupling_iterations {
            check()?;
            checkpoint(self.flow.gate)?;
            if native_solve && self.momentum_trials >= self.flow.max_solves {
                return Err(DepositionError::Budget("live material MHD momentum trials").into());
            }
            let current = evaluate(
                self.input,
                self.motion,
                self.flow.gate,
                &mut self.current_solves,
                thermal,
                &epoch.samples,
                &velocity,
                self.response,
            )?;
            // Every trial starts from the SAME accepted face momentum and inner
            // history. Rejected Picard iterates cannot age either owner.
            let mut candidate = self.flow.clone();
            let mut force = current.solution.lorentz_force_density_n_m3.clone();
            for (c, sample) in epoch.samples.iter().enumerate() {
                check()?;
                if matches!(
                    sample.transport.phase_mobility(),
                    Some(PhaseMobility::ImmobileSolid { .. })
                ) {
                    // Keep this load in the explicit support receipt, not in
                    // the fluid RHS. The native solid-face constraints remain
                    // untouched. Other solid loading still needs its own owner.
                    if epoch.external_force_n_m3.as_ref()
                        .is_some_and(|forces| forces[c].iter().any(|&f| f != 0.0))
                        || epoch.buoyancy_force_n_m3.as_ref()
                            .is_some_and(|forces| forces[c].iter().any(|&f| f != 0.0))
                    {
                        return Err(invalid(
                            "anchored MHD support accounts only for frozen Lorentz loads",
                        ).into());
                    }
                    force[c] = [0.0; 3];
                }
            }
            if let Some(external) = &epoch.external_force_n_m3 {
                for (total, old) in force.iter_mut().zip(external) {
                    check()?;
                    for axis in 0..3 {
                        total[axis] = finite(total[axis] + old[axis])?;
                    }
                }
            }
            candidate
                .applied
                .as_mut()
                .expect("checked source epoch")
                .external_force_n_m3 = Some(force);
            self.momentum_trials += usize::from(native_solve);
            candidate.solves = self.momentum_trials;
            let mut sources = candidate.step_sources(thermal, start_s, dt_s, end_s, check)?;
            let next_velocity = &candidate
                .unsteady_flow()
                .expect("admitted momentum")
                .velocity_m_s;
            for (sample, actual) in epoch.samples.iter().zip(next_velocity) {
                check()?;
                if matches!(
                    sample.transport.phase_mobility(),
                    Some(PhaseMobility::ImmobileSolid { .. })
                ) && actual.iter().any(|&v| v != 0.0)
                {
                    return Err(invalid(
                        "native momentum moved an externally anchored frozen cell",
                    ).into());
                }
            }
            velocity_defect = 0.0_f64;
            for (actual, sampled) in next_velocity.iter().zip(&velocity) {
                check()?;
                for axis in 0..3 {
                    velocity_defect =
                        velocity_defect.max(finite(actual[axis] - sampled[axis])?.abs());
                }
            }
            let lorentz_w = cell_work(
                self.input.domain,
                &current.solution.lorentz_force_density_n_m3,
                next_velocity,
                self.flow.gate,
            )?;
            let other_body_w = match &epoch.external_force_n_m3 {
                Some(force) => cell_work(self.input.domain, force, next_velocity, self.flow.gate)?,
                None => 0.0,
            };
            let buoyancy_w = match &epoch.buoyancy_force_n_m3 {
                Some(force) => cell_work(self.input.domain, force, next_velocity, self.flow.gate)?,
                None => 0.0,
            };
            let mechanical = candidate
                .unsteady_mechanical_step()
                .ok_or(invalid("MHD needs native work receipt"))?;
            let body_defect_w = finite(
                mechanical.power.supplied_body_power_w - lorentz_w - other_body_w - buoyancy_w,
            )?;
            if body_defect_w.abs() > self.input.power_tolerance_w
                || (mechanical.supplied_buoyancy_power_w - buoyancy_w).abs()
                    > self.input.power_tolerance_w
            {
                return Err(
                    invalid("native MHD force lift does not preserve mechanical work").into(),
                );
            }
            work_defect = finite(lorentz_w - current.solution.mechanical_power_w)?;
            let electric = &current.solution.electric;
            let electrode_w = electric.electrodes.iter().try_fold(0.0, |sum, electrode| {
                finite(sum + electrode_power(electrode)?)
            })?;
            let combined_power_defect = finite(electrode_w - electric.total_joule_w - lorentz_w)?;
            if velocity_defect <= self.motion.velocity_tolerance_m_s
                && work_defect.abs() <= self.input.power_tolerance_w
                && combined_power_defect.abs() <= self.input.power_tolerance_w
            {
                let mut mapped_w = 0.0;
                if electric.joule_w.len() != self.input.cell_ids.len() {
                    return Err(invalid("MHD Joule source map changed").into());
                }
                for (&cell, &watts) in self.input.cell_ids.iter().zip(&electric.joule_w) {
                    check()?;
                    if !watts.is_finite() || watts < 0.0 {
                        return Err(invalid("invalid MHD irreversible heat").into());
                    }
                    mapped_w = finite(mapped_w + watts)?;
                    let value = sources.entry(cell).or_insert(0.0);
                    *value = finite(*value + watts)?;
                }
                if (mapped_w - electric.total_joule_w).abs() > self.input.power_tolerance_w
                    || (electrode_w - electric.delivered_w).abs() > self.input.power_tolerance_w
                {
                    return Err(invalid("MHD electrode/Joule deposition handoff failed").into());
                }
                self.pending = Some(Pending {
                    state: current,
                    start_s,
                    dt_s,
                    end_s,
                    mapped_w,
                    electrode_w,
                    lorentz_w,
                    other_body_w,
                    body_defect_w,
                    velocity_defect_m_s: velocity_defect,
                });
                self.flow = candidate;
                return Ok(sources);
            }
            for (sampled, actual) in velocity.iter_mut().zip(next_velocity) {
                check()?;
                for axis in 0..3 {
                    sampled[axis] = finite(
                        sampled[axis] + self.motion.relaxation * (actual[axis] - sampled[axis]),
                    )?;
                }
            }
        }
        Err(MaterialFlowMhdError::CouplingDidNotConverge {
            iterations: self.motion.max_coupling_iterations,
            velocity_defect_m_s: velocity_defect,
            mechanical_power_defect_w: work_defect,
        })
    }

    fn conductances(
        &mut self,
        thermal: &DepositionSimulation<'_>,
        check: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<BTreeMap<u64, f64>, Self::Error> {
        self.flow.conductances(thermal, check).map_err(Into::into)
    }
    fn heat_transfers(
        &mut self,
        thermal: &DepositionSimulation<'_>,
        check: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<Vec<HeatTransfer>, Self::Error> {
        self.flow.heat_transfers(thermal, check).map_err(Into::into)
    }
    fn limit_timestep_s(
        &mut self,
        thermal: &DepositionSimulation<'_>,
        proposed_s: f64,
        check: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<f64, Self::Error> {
        self.flow
            .limit_timestep_s(thermal, proposed_s, check)
            .map_err(Into::into)
    }
    fn event(
        &mut self,
        event: &Event,
        thermal: &DepositionSimulation<'_>,
    ) -> Result<(), Self::Error> {
        self.flow.event(event, thermal)?;
        admit_electrical_map(thermal, self.input, self.flow.gate)?;
        Ok(())
    }
    fn step(
        &mut self,
        start_s: f64,
        dt_s: f64,
        old: &BTreeMap<CellId, f64>,
        thermal: &DepositionSimulation<'_>,
        check: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), Self::Error> {
        let pending = self.pending.take().ok_or(invalid("missing MHD interval"))?;
        if pending.start_s.to_bits() != start_s.to_bits()
            || pending.dt_s.to_bits() != dt_s.to_bits()
            || pending.end_s.to_bits() != thermal.time_s().to_bits()
        {
            return Err(invalid("MHD work interval changed before acceptance").into());
        }
        self.flow.step(start_s, dt_s, old, thermal, check)?;
        let state = &pending.state.solution;
        let w = &mut self.work;
        w.electrode_j = finite(w.electrode_j + pending.electrode_w * dt_s)?;
        w.joule_j = finite(w.joule_j + state.electric.total_joule_w * dt_s)?;
        w.mapped_j = finite(w.mapped_j + pending.mapped_w * dt_s)?;
        w.lorentz_j = finite(w.lorentz_j + pending.lorentz_w * dt_s)?;
        w.sampled_lorentz_j = finite(w.sampled_lorentz_j + state.mechanical_power_w * dt_s)?;
        for (impulse, force) in w.frozen_lorentz_reaction_impulse_n_s.iter_mut()
            .zip(&pending.state.frozen_lorentz_reaction_n)
        {
            check()?;
            for axis in 0..3 {
                impulse[axis] = finite(impulse[axis] + force[axis] * dt_s)?;
            }
        }
        w.other_body_j = finite(w.other_body_j + pending.other_body_w * dt_s)?;
        w.body_defect_j = finite(w.body_defect_j + pending.body_defect_w * dt_s)?;
        w.absolute_transfer_defect_j = finite(
            w.absolute_transfer_defect_j
                + (pending.lorentz_w - state.mechanical_power_w).abs() * dt_s,
        )?;
        w.absolute_electrical_defect_j = finite(
            w.absolute_electrical_defect_j
                + (pending.electrode_w - state.electric.total_joule_w - state.mechanical_power_w)
                    .abs()
                    * dt_s,
        )?;
        w.absolute_body_defect_j =
            finite(w.absolute_body_defect_j + pending.body_defect_w.abs() * dt_s)?;
        w.absolute_heat_defect_j = finite(
            w.absolute_heat_defect_j
                + (state.electric.total_joule_w - pending.mapped_w).abs() * dt_s,
        )?;
        if [
            w.absolute_transfer_defect_j,
            w.absolute_electrical_defect_j,
            w.absolute_body_defect_j,
            w.absolute_heat_defect_j,
        ]
        .iter()
        .any(|&defect| defect > self.input.energy_tolerance_j)
        {
            return Err(invalid(
                "accumulated MHD electrical/mechanical/heat transfer exceeds tolerance",
            )
            .into());
        }
        self.last = Some(pending);
        Ok(())
    }
}

/// Advance source-resolved moving current, native fluid momentum and total h.
///
/// Attach MaterialFlowFeedback::with_unsteady_viscous_heating first. Every actual
/// thermal interval retains one starting material epoch for resistivity, local
/// viscosity, conductivity, buoyancy and any mobile mushy drag. Old accepted
/// native velocity advects h and sets its CFL. The native electrical and Stokes
/// owners then iterate to a common endpoint velocity; I²R and physical
/// molecular/Darcy losses heat the original cells exactly once. The electrical
/// getter must read the supplied retained sample, not re-resolve material or
/// mutate external physical history. None explicitly means insulation.
///
/// The native force lift and cell velocity average are reciprocal. Their work
/// is nevertheless checked independently on the actual accepted velocity, with
/// a separate temporal transfer defect. Electrical dissipation cannot hide
/// incomplete iteration, pressure projection error or BE numerical damping.
/// The latter retains its nonnegative native account and is never heat.
/// Whole-advance current limits include every trial and the final current-only
/// support solve; the flow limit counts every actual momentum trial. An inner
/// duration-dependent observer is cloned from the same accepted history on
/// every trial and advances physically only once at acceptance.
///
/// This is fixed-volume, low-induction, closed unsteady creeping flow with a
/// caller-imposed stationary magnetic field and fixed original contacts. The
/// imposed field source is stationary and supplies no extra mechanical work;
/// electrode power, signed Lorentz transfer, heat and kinetic storage close
/// explicitly. There is no induced-field energy store, Hall term or moving
/// geometry. Solenoidal B and low-induction applicability remain caller claims.
/// Native Reynolds and mechanical-work admission remain in force. The original
/// external body force and Boussinesq source remain composable and keep distinct
/// work accounts. This entry refuses exact solid phase locking: reaction on a
/// conducting frozen skeleton needs an explicit support owner. Select
/// advance_mhd_enthalpy_with_anchored_phase for stationary external anchors.
/// Mobile mushy
/// resistance still uses the original h-dependent material closure.
///
/// The full electrical/thermal/flow map must be identical. Resume with returned
/// thermal state and flow together; a changed drive or B is an explicit new
/// advance. Changing B establishes a newly imposed static field; magnetic
/// switching/Faraday transients are outside this model. Thermal constitutive
/// feedback remains first-order in time and needs
/// refinement. Every electrical, material, convergence, budget, endpoint or
/// cancellation refusal publishes neither state. The report's combined energy
/// residual must satisfy the supplied electrical energy tolerance as well as
/// each original owner's independent gates.
#[allow(clippy::too_many_arguments, clippy::type_complexity)]
pub fn advance_mhd_enthalpy_with_material_flow<'t, 'f, F, R, O, T>(
    thermal: &DepositionSimulation<'t>,
    target_s: f64,
    config: StepConfig,
    input: JouleEnthalpyInput<'_>,
    motion: MhdEnthalpyInput<'_>,
    flow: &MaterialFlowFeedback<'f, F, R, O, T>,
    resistivity: impl Fn(&FlowMaterialSample<F, T>) -> Result<Option<f64>, ChtError>,
) -> Result<
    (
        DepositionSimulation<'t>,
        MaterialFlowFeedback<'f, F, R, O, T>,
        MhdEnthalpyReport,
    ),
    MaterialFlowMhdError<O::Error>,
>
where
    F: Clone,
    T: FlowTransportSnapshot,
    R: Clone + FnMut(usize, &DepositedCell, f64) -> Result<F, MaterialError>,
    O: DepositionObserver,
{
    advance_mhd_enthalpy_impl(
        thermal, target_s, config, input, motion, flow, false, resistivity,
    )
}

/// Continue conducting solid, mushy and liquid material with stationary anchors.
///
/// This explicit boundary model supports the Lorentz load of every exact solid
/// cell with an equal and opposite external force. The force on the material and
/// its accepted impulse remain visible per original voxel; the external anchor
/// receives their negatives. The anchor's velocity and work are exactly zero.
/// This accounts only for direct magnetic loads, not native wall/fluid traction.
///
/// Construct the flow with MaterialFlowFeedback::try_new_phase and attach native
/// unsteady viscous heating. Every original cell remains in the same electrical,
/// material and total-enthalpy map: frozen cells can conduct and heat without
/// liquid-only viscosity or heat capacity. Melting releases the original native
/// faces at zero speed; mobile cells receive the unchanged reciprocal Lorentz
/// force lift and physical molecular/Darcy heat. All-solid intervals solve
/// current without inventing a momentum solve. Phase and transport are sampled
/// at the actual interval start, so timestep refinement remains necessary.
///
/// The original native phase owner still refuses freezing a moving face. This
/// API supplies no capture impulse/heat, moving support, deforming solid,
/// unequal-density phase change or free-surface transition. Nonmagnetic loading
/// on constrained cells still refuses. All low-induction, fixed-field, work,
/// convergence and transactional requirements of
/// advance_mhd_enthalpy_with_material_flow also apply.
#[allow(clippy::too_many_arguments, clippy::type_complexity)]
pub fn advance_mhd_enthalpy_with_anchored_phase<'t, 'f, F, R, O>(
    thermal: &DepositionSimulation<'t>,
    target_s: f64,
    config: StepConfig,
    input: JouleEnthalpyInput<'_>,
    motion: MhdEnthalpyInput<'_>,
    flow: &MaterialFlowFeedback<'f, F, R, O, PhaseFlowTransport>,
    resistivity: impl Fn(&FlowMaterialSample<F, PhaseFlowTransport>) -> Result<Option<f64>, ChtError>,
) -> Result<
    (
        DepositionSimulation<'t>,
        MaterialFlowFeedback<'f, F, R, O, PhaseFlowTransport>,
        MhdEnthalpyReport,
    ),
    MaterialFlowMhdError<O::Error>,
>
where
    F: Clone,
    R: Clone + FnMut(usize, &DepositedCell, f64) -> Result<F, MaterialError>,
    O: DepositionObserver,
{
    advance_mhd_enthalpy_impl(
        thermal, target_s, config, input, motion, flow, true, resistivity,
    )
}

#[allow(clippy::too_many_arguments, clippy::type_complexity)]
fn advance_mhd_enthalpy_impl<'t, 'f, F, R, O, T>(
    thermal: &DepositionSimulation<'t>,
    target_s: f64,
    config: StepConfig,
    input: JouleEnthalpyInput<'_>,
    motion: MhdEnthalpyInput<'_>,
    flow: &MaterialFlowFeedback<'f, F, R, O, T>,
    anchored_skeleton: bool,
    resistivity: impl Fn(&FlowMaterialSample<F, T>) -> Result<Option<f64>, ChtError>,
) -> Result<
    (
        DepositionSimulation<'t>,
        MaterialFlowFeedback<'f, F, R, O, T>,
        MhdEnthalpyReport,
    ),
    MaterialFlowMhdError<O::Error>,
>
where
    F: Clone,
    T: FlowTransportSnapshot,
    R: Clone + FnMut(usize, &DepositedCell, f64) -> Result<F, MaterialError>,
    O: DepositionObserver,
{
    let gate = flow.gate;
    checkpoint(gate)?;
    admit_electrical_map(thermal, input, gate)?;
    flow.check_map(thermal, &mut || Ok(()))?;
    if !core::ptr::eq(input.domain, flow.domain)
        || flow.domain.fluid_count() != flow.domain.cell_count()
        || input.cell_ids.len() != flow.input.cell_ids.len()
        || input
            .cell_ids
            .iter()
            .zip(flow.input.cell_ids)
            .any(|(&e, &f)| f != Some(e))
        || flow.exact_phase != anchored_skeleton
        || !flow
            .unsteady
            .as_ref()
            .is_some_and(MaterialMomentum::has_heating)
    {
        return Err(invalid(
            "MHD needs the same complete material domain/map, explicit phase support and native unsteady physical heating",
        )
        .into());
    }
    if motion.magnetic_flux_density_t.len() != input.cell_ids.len()
        || motion.max_coupling_iterations == 0
        || !motion.relaxation.is_finite()
        || motion.relaxation <= 0.0
        || motion.relaxation > 1.0
        || !motion.velocity_tolerance_m_s.is_finite()
        || motion.velocity_tolerance_m_s <= 0.0
    {
        return Err(invalid("invalid MHD field or coupling limits").into());
    }
    for field in motion.magnetic_flux_density_t {
        checkpoint(gate)?;
        if field.iter().any(|value| !value.is_finite()) {
            return Err(invalid("nonfinite imposed magnetic field").into());
        }
    }
    if let MhdCurrentBoundary::Isolated {
        current_tolerance_a,
    } = motion.current_boundary
        && (!current_tolerance_a.is_finite()
            || current_tolerance_a < 0.0
            || !input.setup.electrodes.is_empty())
    {
        return Err(
            invalid("isolated MHD needs no electrodes and a finite nonnegative KCL limit").into(),
        );
    }
    let source = MhdFlow {
        flow: flow.clone(),
        input,
        motion,
        response: &resistivity,
        current_solves: 0,
        momentum_trials: 0,
        pending: None,
        last: None,
        work: Work {
            frozen_lorentz_reaction_impulse_n_s: vec![[0.0; 3]; input.cell_ids.len()],
            ..Work::default()
        },
    };
    let (next, mut accepted, thermal_audit) =
        thermal.advance_observed(&source, target_s, config, || checkpoint(gate))?;
    let endpoint = evaluate(
        input,
        motion,
        gate,
        &mut accepted.current_solves,
        &next,
        accepted.flow.endpoint(),
        &accepted
            .flow
            .unsteady_flow()
            .expect("admitted momentum")
            .velocity_m_s,
        accepted.response,
    )?;
    let mechanical = accepted
        .flow
        .unsteady_heating_audit()
        .ok_or(invalid("missing accepted MHD mechanical account"))?;
    let work = &accepted.work;
    let other_thermal_input_j =
        finite(thermal_audit.external_heat_j - work.mapped_j - mechanical.generated_heat_j)?;
    let combined_residual_j = finite(
        thermal_audit.stored_energy_change_j + mechanical.final_kinetic_energy_j
            - mechanical.initial_kinetic_energy_j
            + mechanical.backward_euler_loss_j
            - thermal_audit.incoming_energy_j
            - other_thermal_input_j
            - work.electrode_j
            - mechanical.supplied_wall_work_j
            - work.other_body_j
            - mechanical.supplied_buoyancy_work_j
            - mechanical.pressure_constraint_work_j,
    )?;
    if combined_residual_j.abs() > input.energy_tolerance_j {
        return Err(invalid("combined MHD thermal/kinetic/work balance exceeds tolerance").into());
    }
    let report = MhdEnthalpyReport {
        thermal: thermal_audit,
        mechanical,
        electrode_energy_j: work.electrode_j,
        joule_energy_j: work.joule_j,
        deposited_joule_j: work.mapped_j,
        lorentz_work_j: work.lorentz_j,
        sampled_lorentz_work_j: work.sampled_lorentz_j,
        frozen_lorentz_reaction_impulse_n_s: work.frozen_lorentz_reaction_impulse_n_s.clone(),
        frozen_support_work_j: 0.0,
        electrical_defect_j: finite(work.electrode_j - work.joule_j - work.sampled_lorentz_j)?,
        heat_transfer_defect_j: finite(work.joule_j - work.mapped_j)?,
        mechanical_transfer_defect_j: finite(work.lorentz_j - work.sampled_lorentz_j)?,
        absolute_mechanical_transfer_defect_j: work.absolute_transfer_defect_j,
        native_body_transfer_defect_j: work.body_defect_j,
        other_body_work_j: work.other_body_j,
        other_thermal_input_j,
        combined_residual_j,
        current_solves: accepted.current_solves,
        momentum_trials: accepted.momentum_trials,
        last_applied: accepted.last.as_ref().map(|last| last.state.clone()),
        last_duration_s: accepted.last.as_ref().map(|last| last.dt_s),
        last_velocity_defect_m_s: accepted.last.as_ref().map(|last| last.velocity_defect_m_s),
        endpoint,
    };
    checkpoint(gate)?;
    Ok((next, accepted.flow, report))
}
