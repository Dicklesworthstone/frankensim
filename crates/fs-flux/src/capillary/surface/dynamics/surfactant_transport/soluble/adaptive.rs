//! Coupled local time-error control for native soluble capillary motion.
//!
//! One complete original chemical/mechanical map is compared with two half
//! maps. Backward-Euler adsorption makes this composition first order. The
//! unextrapolated fine endpoint must resolve position, velocity, each surface
//! amount, the finite bulk amount and its original nonideal species activity.
//! Constant topology, reference masses, coefficients, external nodal forces,
//! solvent volume and temperature retain the native endpoint's physical scope.

use super::*;
use crate::capillary::surface::surfactant::{
    soluble::adaptive::{
        AdaptiveSolubleSurfaceControls, activity_error, amount_error, log_activity,
    },
    time_control::{self, Limits, Order, Schedule},
};

/// Original chemical error scales plus physical mechanical accuracy and work.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AdaptiveSolubleSurfaceDynamicsControls {
    /// The original amount/activity and time controls. Its relative tolerance
    /// also scales mechanical errors by changes from the accepted state.
    /// The physical transport maximum_step_s caps the full comparison span;
    /// each chemical half of each original map executes exactly one substep.
    pub transport: AdaptiveSolubleSurfaceControls,
    pub absolute_position_tolerance_m: f64,
    pub absolute_velocity_tolerance_m_s: f64,
    /// All original mechanical proposals, including coarse and rejected work.
    /// Every comparison reserves three proposals. Each retains its original
    /// bounded Newton/Krylov policy from SurfaceDynamicsControls.
    pub maximum_mechanical_steps: usize,
}

/// Whole-duration physical accounts, separate from all attempted solver work.
#[derive(Clone, Debug)]
pub struct AdaptiveSolubleSurfaceDynamicsReport {
    pub interval_s: [f64; 2],
    pub trials: usize,
    pub accepted_intervals: usize,
    pub rejected_intervals: usize,
    pub accepted_mechanical_steps: usize,
    pub evaluated_mechanical_steps: usize,
    pub accepted_chemical_substeps: usize,
    pub evaluated_chemical_substeps: usize,
    pub accepted_root_evaluations: usize,
    pub evaluated_root_evaluations: usize,
    pub accepted_surface_pair_exchanges: usize,
    pub evaluated_surface_pair_exchanges: usize,
    pub maximum_accepted_normalized_error: f64,
    pub suggested_next_step_s: f64,
    pub total_amount_mol: [f64; 2],
    pub amount_defect_mol: f64,
    pub gross_adsorption_mol: f64,
    pub gross_desorption_mol: f64,
    pub kinetic_energy_j: [f64; 2],
    /// Reconstructed from original and final momentum, without subtracting
    /// two large absolute kinetic-energy sums.
    pub kinetic_change_j: f64,
    /// Original-to-final bulk mixing plus actual area/amount surface change.
    pub free_energy_change_j: f64,
    pub chemical_free_energy_change_j: f64,
    pub geometric_surface_energy_change_j: f64,
    pub exchange_free_energy_dissipation_j: f64,
    pub surface_diffusion_free_energy_dissipation_j: f64,
    pub surface_viscous_dissipation_j: f64,
    /// Held forces contracted with the complete original-to-final displacement.
    pub external_work_j: f64,
    /// Whole free-energy change minus accepted chemical and geometric changes.
    pub energy_account_residual_j: f64,
    /// Whole kinetic/free-energy change plus physical losses minus held work.
    pub energy_residual_j: f64,
    pub maximum_impulse_residual_n_s: f64,
    pub maximum_kinematic_residual_m: f64,
}

/// No refusal publishes a private prefix, candidate bath or probe history.
#[derive(Debug)]
pub enum AdaptiveSolubleSurfaceDynamicsError {
    Dynamics(SolubleSurfaceDynamicsError),
    TrialBudget { completed_trials: usize },
    MechanicalStepBudget { completed: usize, required: usize },
    ChemicalSubstepBudget { completed: usize, required: usize },
    PairExchangeBudget { completed: usize, required: usize },
    TemporalAccuracy { step_s: f64, normalized_error: f64 },
    TimeResolution { time_s: f64 },
}
impl From<SolubleSurfaceDynamicsError> for AdaptiveSolubleSurfaceDynamicsError {
    fn from(error: SolubleSurfaceDynamicsError) -> Self {
        Self::Dynamics(error)
    }
}
impl From<SurfaceDynamicsError> for AdaptiveSolubleSurfaceDynamicsError {
    fn from(error: SurfaceDynamicsError) -> Self {
        Self::Dynamics(error.into())
    }
}
impl From<SolubleSurfaceError> for AdaptiveSolubleSurfaceDynamicsError {
    fn from(error: SolubleSurfaceError) -> Self {
        Self::Dynamics(SolubleSurfaceDynamicsError::Account(error))
    }
}
impl core::fmt::Display for AdaptiveSolubleSurfaceDynamicsError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "adaptive moving soluble surface: {self:?}")
    }
}
impl std::error::Error for AdaptiveSolubleSurfaceDynamicsError {}
type Result<T> = core::result::Result<T, AdaptiveSolubleSurfaceDynamicsError>;

#[derive(Default)]
struct Work {
    mechanics: usize,
    chemical: usize,
    roots: usize,
    pairs: usize,
}

#[derive(Default)]
struct Ledger {
    work: Work,
    chemical: Sum,
    geometric: Sum,
    exchange: Sum,
    diffusion: Sum,
    viscous: Sum,
    adsorbed: Sum,
    desorbed: Sum,
    maximum_impulse: f64,
    maximum_kinematic: f64,
}
impl Ledger {
    fn accept(&mut self, report: &SolubleSurfaceDynamicsReport) {
        self.work.mechanics += 1;
        self.work.chemical += report.total_substeps;
        self.work.roots += report.total_root_evaluations;
        self.work.pairs += report.total_surface_pair_exchanges;
        self.chemical.add(report.chemical_free_energy_change_j);
        self.geometric.add(report.mechanics.surface.energy_change_j);
        self.exchange.add(report.exchange_free_energy_dissipation_j);
        self.diffusion
            .add(report.surface_diffusion_free_energy_dissipation_j);
        self.viscous
            .add(report.mechanics.surface_viscous_dissipation_j);
        for chemical in [report.before_motion, report.after_motion] {
            self.adsorbed.add(chemical.gross_adsorption_mol);
            self.desorbed.add(chemical.gross_desorption_mol);
        }
        self.maximum_impulse = self
            .maximum_impulse
            .max(report.mechanics.impulse_residual_n_s);
        self.maximum_kinematic = self
            .maximum_kinematic
            .max(report.mechanics.kinematic_residual_m);
    }
}

struct Candidate {
    state: SecondOrderState,
    bulk: AqueousSolution,
    amounts: Vec<f64>,
}
impl Candidate {
    fn new(
        state: &SecondOrderState,
        bulk: &AqueousSolution,
        amounts: &[f64],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self> {
        Ok(Self {
            state: SecondOrderState {
                t: state.t,
                q: copied(&state.q, cancelled)?,
                v: copied(&state.v, cancelled)?,
                a: copied(&state.a, cancelled)?,
                steps: state.steps,
                history: Vec::new(),
            },
            bulk: bulk.clone(),
            amounts: copied(amounts, cancelled)?,
        })
    }

    fn reset(&mut self, original: &Self, cancelled: &mut impl FnMut() -> bool) -> Result<()> {
        for (from, into) in [
            (&original.state.q, &mut self.state.q),
            (&original.state.v, &mut self.state.v),
            (&original.state.a, &mut self.state.a),
            (&original.amounts, &mut self.amounts),
        ] {
            for (i, (&value, out)) in from.iter().zip(into).enumerate() {
                poll_at(cancelled, i)?;
                *out = value;
            }
        }
        self.state.t = original.state.t;
        self.state.steps = original.state.steps;
        self.state.history.clear();
        self.bulk
            .replace_inventories(
                original.bulk.solvent_volume_m3(),
                original.bulk.solute_moles(),
            )
            .map_err(SolubleSurfaceError::from)?;
        Ok(())
    }
}

impl CapillarySurfaceDynamics<'_> {
    /// Complete a moving soluble-surface duration to an exact absolute target.
    ///
    /// Native coupled endpoints supply every coarse and fine proposal. Their
    /// chemistry uses actual initial/solved areas; mechanics retains original
    /// masses, forces, viscosity, packing and independent physical gates. The
    /// first-order fine estimate compares position, velocity, every face amount
    /// and the original finite bulk amount/activity. Actual representable half
    /// times determine the correction, including on a nonzero accepted clock.
    ///
    /// Coarse and rejected proposals share the whole-call mechanical, chemical,
    /// root and pair caps. Accepted fine receipts alone supply physical losses
    /// and telemetry. Only endpoint arrays are copied for trials; original
    /// accumulated history is neither copied nor published as a private prefix.
    ///
    /// Independent original-to-final amount, geometric/chemical and combined
    /// mechanical/free-energy accounts must pass the unchanged fixed tolerances
    /// before all three states and accepted history publish once. Only excessive
    /// temporal error is retried; physical/source/solver refusals propagate.
    /// Target time must be strictly later. Held forces, coefficients, topology,
    /// reference mass, temperature and solvent volume must be cut at events.
    /// Local temporal estimates do not certify global or spatial accuracy.
    #[allow(clippy::too_many_arguments, clippy::too_many_lines)]
    pub fn advance_soluble_surfactant_adaptive_to(
        &self,
        state: &mut SecondOrderState,
        bulk: &mut AqueousSolution,
        amounts_mol: &mut [f64],
        transport: SolubleSurfaceTransport<'_>,
        viscosity: SurfaceViscosityField<'_>,
        external_forces_n: &[Vector],
        target_s: f64,
        mechanics: SurfaceDynamicsControls,
        accuracy: AdaptiveSolubleSurfaceDynamicsControls,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<AdaptiveSolubleSurfaceDynamicsReport> {
        poll(&mut cancelled)?;
        let temporal = accuracy.transport;
        let mut schedule = Schedule::interval(
            state.t,
            target_s,
            Limits {
                initial_step_s: temporal.initial_step_s,
                minimum_step_s: temporal.minimum_step_s,
                maximum_step_s: transport.controls.maximum_step_s,
                maximum_trials: temporal.maximum_trials,
            },
            Order::First,
        )
        .map_err(schedule_error)?;
        if target_s <= state.t
            || state.q.len() != 3 * self.masses_kg.len()
            || state.v.len() != state.q.len()
            || state.a.len() != state.q.len()
            || amounts_mol.len() != self.faces.len()
            || !vectors(external_forces_n, self.masses_kg.len())
            || [
                accuracy.absolute_position_tolerance_m,
                accuracy.absolute_velocity_tolerance_m_s,
                temporal.absolute_surface_amount_tolerance_mol,
                temporal.absolute_bulk_amount_tolerance_mol,
                temporal.absolute_bulk_activity_tolerance,
            ]
            .iter()
            .any(|value| !value.is_finite() || *value <= 0.0)
            || !temporal.relative_tolerance.is_finite()
            || !(0.0..1.0).contains(&temporal.relative_tolerance)
        {
            return Err(SurfaceDynamicsError::Invalid(
                "complete moving state, positive duration and finite positive physical error scales required",
            ).into());
        }
        if self.faces.len() > transport.controls.maximum_faces
            || transport.edges.len() > transport.controls.maximum_edges
        {
            return Err(SolubleSurfaceError::WorkBudget.into());
        }
        let start_s = state.t;
        let before_areas = face_areas(self.faces, &state.q, &mut cancelled)?;
        let mut current = Candidate::new(state, bulk, amounts_mol, &mut cancelled)?;
        // Preserve native source, layout, packing and graph admission even
        // before a work-budget refusal. This performs no physical time step.
        let initial = advance_soluble_surface(
            &mut current.bulk,
            transport.laws,
            &before_areas,
            &mut current.amounts,
            transport.kinetics,
            transport.edges,
            0.0,
            transport.controls,
            &mut cancelled,
        )
        .map_err(|error| stage_error(SurfaceDiffusionPhase::BeforeMotion, error))?;
        let mut coarse = Candidate::new(
            &current.state,
            &current.bulk,
            &current.amounts,
            &mut cancelled,
        )?;
        let mut fine = Candidate::new(
            &current.state,
            &current.bulk,
            &current.amounts,
            &mut cancelled,
        )?;
        let mut current_activity = log_activity(&current.bulk, transport.kinetics.bulk_species)?;
        let mut work = Work::default();
        let mut ledger = Ledger::default();
        let mut history = Vec::new();
        let pairs_per_trial = if transport.edges.is_empty() {
            0
        } else {
            transport
                .edges
                .len()
                .checked_mul(2)
                .and_then(|n| n.checked_sub(1))
                .and_then(|n| n.checked_mul(12))
                .ok_or(SolubleSurfaceError::WorkBudget)?
        };
        while !schedule.complete() {
            poll(&mut cancelled)?;
            schedule.check_trial_budget().map_err(schedule_error)?;
            if 3 > accuracy
                .maximum_mechanical_steps
                .saturating_sub(work.mechanics)
            {
                return Err(AdaptiveSolubleSurfaceDynamicsError::MechanicalStepBudget {
                    completed: work.mechanics,
                    required: 3,
                });
            }
            if 6 > transport
                .controls
                .maximum_substeps
                .saturating_sub(work.chemical)
            {
                return Err(AdaptiveSolubleSurfaceDynamicsError::ChemicalSubstepBudget {
                    completed: work.chemical,
                    required: 6,
                });
            }
            if pairs_per_trial
                > transport
                    .controls
                    .maximum_surface_pair_exchanges
                    .saturating_sub(work.pairs)
            {
                return Err(AdaptiveSolubleSurfaceDynamicsError::PairExchangeBudget {
                    completed: work.pairs,
                    required: pairs_per_trial,
                });
            }
            let trial = schedule.begin().map_err(schedule_error)?;
            coarse.reset(&current, &mut cancelled)?;
            fine.reset(&current, &mut cancelled)?;
            self.adaptive_soluble_step(
                &mut coarse,
                transport,
                viscosity,
                external_forces_n,
                trial.span_s,
                mechanics,
                &mut work,
                &mut cancelled,
            )?;
            let first = self.adaptive_soluble_step(
                &mut fine,
                transport,
                viscosity,
                external_forces_n,
                trial.halves_s[0],
                mechanics,
                &mut work,
                &mut cancelled,
            )?;
            let second = self.adaptive_soluble_step(
                &mut fine,
                transport,
                viscosity,
                external_forces_n,
                trial.halves_s[1],
                mechanics,
                &mut work,
                &mut cancelled,
            )?;
            if coarse.state.t.to_bits() != trial.end_s.to_bits()
                || fine.state.t.to_bits() != trial.end_s.to_bits()
            {
                return Err(AdaptiveSolubleSurfaceDynamicsError::TimeResolution {
                    time_s: current.state.t,
                });
            }
            let coarse_activity = log_activity(&coarse.bulk, transport.kinetics.bulk_species)?;
            let fine_activity = log_activity(&fine.bulk, transport.kinetics.bulk_species)?;
            let mut error = activity_error(
                current_activity,
                coarse_activity,
                fine_activity,
                temporal.absolute_bulk_activity_tolerance,
                temporal.relative_tolerance,
                trial.fine_error_factor,
            );
            for (before, full, halves, allowance) in [
                (
                    &current.state.q,
                    &coarse.state.q,
                    &fine.state.q,
                    accuracy.absolute_position_tolerance_m,
                ),
                (
                    &current.state.v,
                    &coarse.state.v,
                    &fine.state.v,
                    accuracy.absolute_velocity_tolerance_m_s,
                ),
                (
                    &current.amounts,
                    &coarse.amounts,
                    &fine.amounts,
                    temporal.absolute_surface_amount_tolerance_mol,
                ),
            ] {
                for (i, ((&start, &a), &b)) in before.iter().zip(full).zip(halves).enumerate() {
                    poll_at(&mut cancelled, i)?;
                    error = error.max(amount_error(
                        start,
                        a,
                        b,
                        allowance,
                        temporal.relative_tolerance,
                        trial.fine_error_factor,
                    ));
                }
            }
            let species = transport.kinetics.bulk_species;
            error = error.max(amount_error(
                current.bulk.solute_moles()[species],
                coarse.bulk.solute_moles()[species],
                fine.bulk.solute_moles()[species],
                temporal.absolute_bulk_amount_tolerance_mol,
                temporal.relative_tolerance,
                trial.fine_error_factor,
            ));
            if !schedule.consider(trial, error).map_err(schedule_error)? {
                continue;
            }
            ledger.accept(&first);
            ledger.accept(&second);
            history
                .try_reserve(fine.state.history.len())
                .map_err(|_| SurfaceDynamicsError::Allocation)?;
            history.append(&mut fine.state.history);
            core::mem::swap(&mut current, &mut fine);
            current_activity = fine_activity;
        }
        let species = transport.kinetics.bulk_species;
        let mut defect = Sum::default();
        let mut surface_after = Sum::default();
        defect.add(current.bulk.solute_moles()[species] - bulk.solute_moles()[species]);
        for (i, (&old, &new)) in amounts_mol.iter().zip(&current.amounts).enumerate() {
            poll_at(&mut cancelled, i)?;
            defect.add(new - old);
            surface_after.add(new);
        }
        let defect = finite(defect.value())?;
        if defect.abs() > transport.controls.amount_tolerance_mol {
            return Err(SolubleSurfaceError::AmountBalance { defect_mol: defect }.into());
        }
        let after_areas = face_areas(self.faces, &current.state.q, &mut cancelled)?;
        let free_change = energy_change(
            bulk,
            species,
            bulk.solute_moles(),
            current.bulk.solute_moles(),
            transport.laws,
            &before_areas,
            &after_areas,
            amounts_mol,
            &current.amounts,
            transport.kinetics.standard_potential_offset_j_mol,
            &mut cancelled,
        )?;
        let chemical = finite(ledger.chemical.value())?;
        let geometric = finite(ledger.geometric.value())?;
        if chemical > transport.controls.energy_tolerance_j {
            return Err(SolubleSurfaceError::FreeEnergy { change_j: chemical }.into());
        }
        let account = sum(&[free_change, -chemical, -geometric])?;
        if account.abs() > transport.controls.energy_tolerance_j {
            return Err(SolubleSurfaceError::EnergyAccount {
                residual_j: account,
            }
            .into());
        }
        let mut kinetic_change = Sum::default();
        let mut external_work = Sum::default();
        for i in 0..state.q.len() {
            poll_at(&mut cancelled, i)?;
            kinetic_change.add(
                self.masses_kg[i / 3]
                    * (current.state.v[i] - state.v[i])
                    * (0.5 * current.state.v[i] + 0.5 * state.v[i]),
            );
            external_work.add(external_forces_n[i / 3][i % 3] * (current.state.q[i] - state.q[i]));
        }
        let kinetic_change = finite(kinetic_change.value())?;
        let external_work = finite(external_work.value())?;
        let exchange = finite(ledger.exchange.value())?;
        let diffusion = finite(ledger.diffusion.value())?;
        let viscous = finite(ledger.viscous.value())?;
        let energy = sum(&[
            kinetic_change,
            free_change,
            exchange,
            diffusion,
            viscous,
            -external_work,
        ])?;
        if energy.abs() > mechanics.energy_tolerance_j {
            return Err(SurfaceDynamicsError::Balance {
                impulse_n_s: ledger.maximum_impulse,
                kinematic_m: ledger.maximum_kinematic,
                energy_j: energy,
            }
            .into());
        }
        let statistics = schedule.statistics;
        let report = AdaptiveSolubleSurfaceDynamicsReport {
            interval_s: [start_s, target_s],
            trials: statistics.trials,
            accepted_intervals: statistics.accepted_intervals,
            rejected_intervals: statistics.rejected_intervals,
            accepted_mechanical_steps: ledger.work.mechanics,
            evaluated_mechanical_steps: work.mechanics,
            accepted_chemical_substeps: ledger.work.chemical,
            evaluated_chemical_substeps: work.chemical,
            accepted_root_evaluations: ledger.work.roots,
            evaluated_root_evaluations: work.roots,
            accepted_surface_pair_exchanges: ledger.work.pairs,
            evaluated_surface_pair_exchanges: work.pairs,
            maximum_accepted_normalized_error: statistics.maximum_accepted_normalized_error,
            suggested_next_step_s: statistics.suggested_next_step_s,
            total_amount_mol: [
                sum(&[
                    initial.bulk_amount_before_mol,
                    initial.surface_amount_before_mol,
                ])?,
                sum(&[current.bulk.solute_moles()[species], surface_after.value()])?,
            ],
            amount_defect_mol: defect,
            gross_adsorption_mol: finite(ledger.adsorbed.value())?,
            gross_desorption_mol: finite(ledger.desorbed.value())?,
            kinetic_energy_j: [
                kinetic(self.masses_kg, &state.v)?,
                kinetic(self.masses_kg, &current.state.v)?,
            ],
            kinetic_change_j: kinetic_change,
            free_energy_change_j: free_change,
            chemical_free_energy_change_j: chemical,
            geometric_surface_energy_change_j: geometric,
            exchange_free_energy_dissipation_j: exchange,
            surface_diffusion_free_energy_dissipation_j: diffusion,
            surface_viscous_dissipation_j: viscous,
            external_work_j: external_work,
            energy_account_residual_j: account,
            energy_residual_j: energy,
            maximum_impulse_residual_n_s: ledger.maximum_impulse,
            maximum_kinematic_residual_m: ledger.maximum_kinematic,
        };
        current.state.history = history;
        // Reserve the only extension of original history before changing any
        // physical value. Trials never copy or mutate that accepted prefix.
        state
            .history
            .try_reserve(current.state.history.len())
            .map_err(|_| SurfaceDynamicsError::Allocation)?;
        poll(&mut cancelled)?;
        publish_state(state, current.state);
        *bulk = current.bulk;
        amounts_mol.copy_from_slice(&current.amounts);
        Ok(report)
    }

    #[allow(clippy::too_many_arguments)]
    fn adaptive_soluble_step(
        &self,
        candidate: &mut Candidate,
        transport: SolubleSurfaceTransport<'_>,
        viscosity: SurfaceViscosityField<'_>,
        external_forces_n: &[Vector],
        dt_s: f64,
        mechanics: SurfaceDynamicsControls,
        work: &mut Work,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<SolubleSurfaceDynamicsReport> {
        let roots = transport
            .controls
            .maximum_root_evaluations
            .saturating_sub(work.roots);
        if roots == 0 {
            return Err(stage_error(
                SurfaceDiffusionPhase::BeforeMotion,
                SolubleSurfaceError::RootBudget {
                    completed_evaluations: work.roots,
                },
            )
            .into());
        }
        let controls = SolubleSurfaceControls {
            maximum_step_s: 0.5 * dt_s,
            maximum_substeps: 2,
            maximum_root_evaluations: roots,
            maximum_surface_pair_exchanges: transport
                .controls
                .maximum_surface_pair_exchanges
                .saturating_sub(work.pairs),
            ..transport.controls
        };
        let report = self
            .advance_soluble_surfactant(
                &mut candidate.state,
                &mut candidate.bulk,
                &mut candidate.amounts,
                SolubleSurfaceTransport {
                    controls,
                    ..transport
                },
                viscosity,
                external_forces_n,
                dt_s,
                mechanics,
                cancelled,
            )
            .map_err(|error| match error {
                SolubleSurfaceDynamicsError::Transport {
                    phase,
                    error:
                        SolubleSurfaceError::RootBudget {
                            completed_evaluations,
                        },
                } => stage_error(
                    phase,
                    SolubleSurfaceError::RootBudget {
                        completed_evaluations: work.roots + completed_evaluations,
                    },
                ),
                other => other,
            })?;
        work.mechanics += 1;
        work.chemical += report.total_substeps;
        work.roots += report.total_root_evaluations;
        work.pairs += report.total_surface_pair_exchanges;
        Ok(report)
    }
}

fn copied(source: &[f64], cancelled: &mut impl FnMut() -> bool) -> Result<Vec<f64>> {
    let mut output = Vec::new();
    output
        .try_reserve_exact(source.len())
        .map_err(|_| SurfaceDynamicsError::Allocation)?;
    for (i, &value) in source.iter().enumerate() {
        poll_at(cancelled, i)?;
        output.push(value);
    }
    Ok(output)
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<()> {
    if cancelled() {
        Err(SurfaceDynamicsError::Cancelled.into())
    } else {
        Ok(())
    }
}
fn poll_at(cancelled: &mut impl FnMut() -> bool, index: usize) -> Result<()> {
    if index % 192 == 0 {
        poll(cancelled)?;
    }
    Ok(())
}
fn schedule_error(error: time_control::Error) -> AdaptiveSolubleSurfaceDynamicsError {
    match error {
        time_control::Error::Invalid => SurfaceDynamicsError::Invalid(
            "finite forward target and positive ordered time/work limits required",
        )
        .into(),
        time_control::Error::TrialBudget { completed_trials } => {
            AdaptiveSolubleSurfaceDynamicsError::TrialBudget { completed_trials }
        }
        time_control::Error::TemporalAccuracy {
            step_s,
            normalized_error,
        } => AdaptiveSolubleSurfaceDynamicsError::TemporalAccuracy {
            step_s,
            normalized_error,
        },
        time_control::Error::TimeResolution { elapsed_s } => {
            AdaptiveSolubleSurfaceDynamicsError::TimeResolution { time_s: elapsed_s }
        }
    }
}

#[cfg(test)]
mod tests;
