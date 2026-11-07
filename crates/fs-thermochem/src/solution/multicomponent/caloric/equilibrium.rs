//! Finite heat transactions on an accepted local liquid phase branch.
//!
//! The original caloric owner supplies every pure Cp and Cp/T integral. The
//! original equilibrium solver redistributes the SAME conserved inventory at
//! the SAME phase count. Its unshifted analytic caloric response predicts each
//! bounded temperature continuation step and supplies the shared caloric
//! inverse's Newton tangent. Actual enthalpy increments decide convergence.
//!
//! This is bounded seeded numerical continuation, not a mathematical branch
//! identity certificate or a global Gibbs envelope. A nonzero heat transaction
//! requires the complete offered thermal bracket to be reachable without phase
//! loss, coincident phases or an unresolved joint Hessian. Zero heat with the
//! original temperature inside that bracket re-admits the original response
//! directly, with no phase solves or branch claim elsewhere in the bracket.
//! No phase insertion, latent-heat smearing,
//! concentration floor or regularized physical heat capacity is introduced.
//! The explicit [`discovery`] child separately admits changing phase counts
//! through the original finite search without altering fixed-branch continuation.

/// Conserved thermal endpoints with the original bounded phase-count discovery.
pub mod discovery;
/// Conserved relaxation and mixing of actual homogeneous or separated inputs.
pub mod transactions;

use super::{
    MulticomponentCaloricErrorV1, MulticomponentCaloricSolveV1, MulticomponentLiquidCaloricModelV1,
    Owner, R, Sum, checked, log_ratio, nonnegative, positive, scaled,
};
use crate::solution::caloric::temperature::{self, TemperatureControl, TemperatureSample};
use crate::solution::caloric::{BinaryCaloricErrorV1, PureLiquidHeatCapacityV1};
use crate::solution::ln_complement;
use crate::solution::multicomponent::discovery::MulticomponentPhaseDiscoveryErrorV1;
use crate::solution::multicomponent::equilibrium::caloric::{
    MulticomponentCaloricResponseControlV1, MulticomponentCaloricResponseErrorV1,
};
use crate::solution::multicomponent::equilibrium::{
    MulticomponentEquilibriumCandidateV1, MulticomponentEquilibriumErrorV1,
    MulticomponentEquilibriumSeedV1,
};
use fs_math::det::ln;
use fs_qty::Temperature;
use std::sync::Arc;

/// Existing thermal tolerances plus bounded local continuation work.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MulticomponentEquilibriumCaloricSolveV1 {
    pub thermal: MulticomponentCaloricSolveV1,
    pub response: MulticomponentCaloricResponseControlV1,
    /// Positive maximum temperature advance per original phase solve [K].
    /// The analytic predictor also limits every positive phase-component
    /// amount's predicted relative change to one quarter per step.
    pub maximum_temperature_step_k: f64,
    /// Whole-call phase-solve bound, 1..=65_536 across ALL root samples.
    /// Each phase solve retains the initial candidate's original inner limits.
    pub maximum_phase_solves: usize,
}

/// A refused branch, source or numerical budget publishes no thermal endpoint.
#[derive(Clone, Debug, PartialEq)]
pub enum MulticomponentEquilibriumCaloricErrorV1 {
    Invalid(&'static str),
    Caloric(MulticomponentCaloricErrorV1),
    Equilibrium(MulticomponentEquilibriumErrorV1),
    Response(MulticomponentCaloricResponseErrorV1),
    Discovery(MulticomponentPhaseDiscoveryErrorV1),
    /// Remaining whole-call work cannot admit another complete, unchanged
    /// configured discovery policy. No shortened final search is substituted.
    DiscoverySolveBudget {
        completed_candidate_solves: usize,
        required_trial_allowance: usize,
    },
    PhaseSolveBudget {
        completed_phase_solves: usize,
    },
    /// Actual finite H/S increments violate the positive-Cp branch's endpoint
    /// Clausius bounds after the caller's entropy/energy allowances.
    ClausiusBounds {
        entropy_change_j_per_k: f64,
        bounds_j_per_k: [f64; 2],
        tolerance_j_per_k: f64,
    },
    Unrepresentable(&'static str),
    Cancelled,
}
impl core::fmt::Display for MulticomponentEquilibriumCaloricErrorV1 {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "multicomponent local equilibrium calorics: {self:?}")
    }
}
impl std::error::Error for MulticomponentEquilibriumCaloricErrorV1 {}
impl From<MulticomponentCaloricErrorV1> for MulticomponentEquilibriumCaloricErrorV1 {
    fn from(error: MulticomponentCaloricErrorV1) -> Self {
        match error {
            MulticomponentCaloricErrorV1::Cancelled => Self::Cancelled,
            other => Self::Caloric(other),
        }
    }
}
impl From<BinaryCaloricErrorV1> for MulticomponentEquilibriumCaloricErrorV1 {
    fn from(error: BinaryCaloricErrorV1) -> Self {
        match error {
            BinaryCaloricErrorV1::Cancelled => Self::Cancelled,
            other => Self::Caloric(other.into()),
        }
    }
}
impl From<MulticomponentEquilibriumErrorV1> for MulticomponentEquilibriumCaloricErrorV1 {
    fn from(error: MulticomponentEquilibriumErrorV1) -> Self {
        match error {
            MulticomponentEquilibriumErrorV1::Cancelled => Self::Cancelled,
            other => Self::Equilibrium(other),
        }
    }
}
impl From<MulticomponentCaloricResponseErrorV1> for MulticomponentEquilibriumCaloricErrorV1 {
    fn from(error: MulticomponentCaloricResponseErrorV1) -> Self {
        match error {
            MulticomponentCaloricResponseErrorV1::Cancelled => Self::Cancelled,
            other => Self::Response(other),
        }
    }
}
impl From<MulticomponentPhaseDiscoveryErrorV1> for MulticomponentEquilibriumCaloricErrorV1 {
    fn from(error: MulticomponentPhaseDiscoveryErrorV1) -> Self {
        match error {
            MulticomponentPhaseDiscoveryErrorV1::Cancelled => Self::Cancelled,
            other => Self::Discovery(other),
        }
    }
}
type Result<T> = core::result::Result<T, MulticomponentEquilibriumCaloricErrorV1>;

/// Exact accepted candidate bound to the original immutable pure-Cp owner.
/// The stored tangent retains that candidate's numerical/local qualifications.
#[derive(Clone, Debug, PartialEq)]
pub struct MulticomponentEquilibriumCaloricStateV1 {
    owner: Arc<Owner>,
    candidate: MulticomponentEquilibriumCandidateV1,
    enthalpy_j: f64,
    frozen_capacity_j_per_k: f64,
    redistribution_capacity_j_per_k: f64,
    capacity_j_per_k: f64,
    phase_component_derivatives: Box<[f64]>,
    response_controls: MulticomponentCaloricResponseControlV1,
    maximum_scaled_response_residual: f64,
}
impl MulticomponentEquilibriumCaloricStateV1 {
    pub const fn candidate(&self) -> &MulticomponentEquilibriumCandidateV1 {
        &self.candidate
    }
    pub const fn temperature(&self) -> Temperature {
        self.candidate.feed().temperature()
    }
    pub fn component_moles(&self) -> &[f64] {
        self.candidate.component_moles()
    }
    pub const fn total_moles(&self) -> f64 {
        self.candidate.total_moles()
    }
    pub fn heat_capacities(&self) -> &[PureLiquidHeatCapacityV1] {
        &self.owner.heat_capacity
    }
    /// Pure sensible enthalpy plus the actual accepted phases' mixing enthalpy.
    /// Heat transactions use direct increments, never subtraction of this datum.
    pub const fn enthalpy_j(&self) -> f64 {
        self.enthalpy_j
    }
    pub const fn frozen_heat_capacity_j_per_k(&self) -> f64 {
        self.frozen_capacity_j_per_k
    }
    pub const fn redistribution_heat_capacity_j_per_k(&self) -> f64 {
        self.redistribution_capacity_j_per_k
    }
    pub const fn isobaric_heat_capacity_j_per_k(&self) -> f64 {
        self.capacity_j_per_k
    }
    pub const fn response_controls(&self) -> MulticomponentCaloricResponseControlV1 {
        self.response_controls
    }
    pub const fn maximum_scaled_response_residual(&self) -> f64 {
        self.maximum_scaled_response_residual
    }
}

/// Conserved signed heat and system entropy between actual admitted states.
#[derive(Clone, Debug, PartialEq)]
pub struct MulticomponentEquilibriumCaloricBalanceV1 {
    final_state: MulticomponentEquilibriumCaloricStateV1,
    heat_input_j: f64,
    enthalpy_change_j: f64,
    entropy_change_j_per_k: f64,
    residual_j: f64,
    temperature_bracket: [Temperature; 2],
    iterations: usize,
    completed_phase_solves: usize,
}
impl MulticomponentEquilibriumCaloricBalanceV1 {
    pub const fn final_state(&self) -> &MulticomponentEquilibriumCaloricStateV1 {
        &self.final_state
    }
    pub const fn heat_input_j(&self) -> f64 {
        self.heat_input_j
    }
    pub const fn enthalpy_change_j(&self) -> f64 {
        self.enthalpy_change_j
    }
    /// System entropy change; reservoir entropy is not inferred from heat alone.
    pub const fn entropy_change_j_per_k(&self) -> f64 {
        self.entropy_change_j_per_k
    }
    pub const fn residual_j(&self) -> f64 {
        self.residual_j
    }
    pub const fn temperature_bracket(&self) -> [Temperature; 2] {
        self.temperature_bracket
    }
    pub const fn iterations(&self) -> usize {
        self.iterations
    }
    pub const fn completed_phase_solves(&self) -> usize {
        self.completed_phase_solves
    }
}

impl MulticomponentLiquidCaloricModelV1 {
    /// Bind an already accepted candidate without replacing or re-solving it.
    /// The complete thermodynamic owner must match, and the original joint
    /// Hessian must admit the physical local Cp without regularization.
    pub fn equilibrium_caloric_state(
        &self,
        candidate: MulticomponentEquilibriumCandidateV1,
        response_controls: MulticomponentCaloricResponseControlV1,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<MulticomponentEquilibriumCaloricStateV1> {
        poll(&mut cancelled)?;
        if candidate.feed().model() != self.thermodynamics() {
            return Err(MulticomponentEquilibriumCaloricErrorV1::Invalid(
                "equilibrium requires the same complete original caloric potential",
            ));
        }
        let temperature = candidate.feed().temperature();
        self.admit_temperature(temperature)?;
        let response = candidate.isobaric_caloric_response(
            self.heat_capacities(),
            response_controls,
            &mut cancelled,
        )?;
        let mut derivatives =
            Vec::with_capacity(candidate.phases().len() * candidate.component_moles().len());
        for phase in 0..candidate.phases().len() {
            poll(&mut cancelled)?;
            derivatives.extend_from_slice(
                response
                    .phase_component_moles_temperature_derivatives_mol_per_k(phase)
                    .expect("original response has every accepted phase"),
            );
        }
        let mut enthalpy = Sum::default();
        for (&amount, source) in candidate
            .component_moles()
            .iter()
            .zip(self.heat_capacities())
        {
            poll(&mut cancelled)?;
            if amount != 0.0 {
                enthalpy.add(scaled(
                    amount,
                    source
                        .enthalpy_change_j_per_mol(source.reference_temperature(), temperature)?,
                    "equilibrium pure-liquid enthalpy",
                )?);
            }
        }
        for phase in candidate.phases() {
            poll(&mut cancelled)?;
            enthalpy.add(scaled(
                phase.total_moles(),
                phase.thermodynamics().mixing_enthalpy_j_per_mol(),
                "equilibrium phase mixing enthalpy",
            )?);
        }
        let frozen_capacity_j_per_k = response.frozen_heat_capacity_j_per_k();
        let redistribution_capacity_j_per_k = response.redistribution_heat_capacity_j_per_k();
        let capacity_j_per_k = response.isobaric_heat_capacity_j_per_k();
        let maximum_scaled_response_residual = response.maximum_scaled_residual();
        poll(&mut cancelled)?;
        Ok(MulticomponentEquilibriumCaloricStateV1 {
            owner: Arc::clone(&self.owner),
            candidate,
            enthalpy_j: checked(enthalpy.value(), "total equilibrium liquid enthalpy")?,
            frozen_capacity_j_per_k,
            redistribution_capacity_j_per_k,
            capacity_j_per_k,
            phase_component_derivatives: derivatives.into_boxed_slice(),
            response_controls,
            maximum_scaled_response_residual,
        })
    }

    /// Continue the original accepted phase branch to a prescribed temperature.
    /// The same predictor, species conservation, unregularized response and
    /// whole-call phase-solve budget used by finite heating remain in force.
    /// Both temperatures must belong to the offered thermal bracket.
    ///
    /// This supplies a caloric query, not a prescribed heat transaction. The
    /// actual H/S increments must satisfy the positive-Cp branch's endpoint
    /// Clausius bounds. A disappearing phase or unresolved branch refuses;
    /// discovery cannot silently replace the original branch.
    pub fn equilibrium_at_temperature(
        &self,
        initial: &MulticomponentEquilibriumCaloricStateV1,
        target: Temperature,
        controls: MulticomponentEquilibriumCaloricSolveV1,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<MulticomponentEquilibriumCaloricStateV1> {
        poll(&mut cancelled)?;
        self.admit_equilibrium_caloric(initial)?;
        let [lo, hi] = self.admit_equilibrium_solve(0.0, controls)?;
        if !target.value().is_finite()
            || !(lo..=hi).contains(&target.value())
            || !(lo..=hi).contains(&initial.temperature().value())
        {
            return Err(MulticomponentEquilibriumCaloricErrorV1::Invalid(
                "branch query temperatures outside the offered thermal bracket",
            ));
        }
        let admitted = self.equilibrium_caloric_state(
            initial.candidate.clone(),
            controls.response,
            &mut cancelled,
        )?;
        let mut completed = 0;
        let final_state = self.continue_equilibrium_temperature(
            admitted,
            target,
            controls,
            &mut completed,
            &mut cancelled,
        )?;
        let heat = self.equilibrium_caloric_change(initial, &final_state, false, &mut cancelled)?;
        let entropy =
            self.equilibrium_caloric_change(initial, &final_state, true, &mut cancelled)?;
        let delta_t = target.value() - initial.temperature().value();
        if delta_t != 0.0
            && (heat == 0.0
                || heat.signum() != delta_t.signum()
                || entropy == 0.0
                || entropy.signum() != delta_t.signum())
        {
            return Err(MulticomponentEquilibriumCaloricErrorV1::Unrepresentable(
                "temperature query does not resolve positive heat and entropy response",
            ));
        }
        admit_clausius_change(initial, &final_state, heat, entropy, controls.thermal)?;
        poll(&mut cancelled)?;
        Ok(final_state)
    }

    /// Add signed heat while continuing the original accepted phase count and
    /// support. Each local predictor preserves every species column, then the
    /// original equilibrium solver corrects it. Every corrected point requires
    /// a resolved unshifted Cp and remains inside the predictor neighborhood.
    ///
    /// The original candidate's phase controls stay unchanged. An unsupported
    /// bracket or disappearing/singular phase branch refuses; this operation
    /// does not fall back to a new phase count or a frozen-composition Cp.
    /// Reported finite entropy must lie between DeltaH/T_initial and
    /// DeltaH/T_final, with the declared entropy allowance, energy allowance
    /// divided by the smaller endpoint T, and 64 eps of arithmetic scale.
    pub fn add_heat_equilibrium(
        &self,
        initial: &MulticomponentEquilibriumCaloricStateV1,
        heat_input_j: f64,
        controls: MulticomponentEquilibriumCaloricSolveV1,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<MulticomponentEquilibriumCaloricBalanceV1> {
        poll(&mut cancelled)?;
        self.admit_equilibrium_caloric(initial)?;
        let thermal = controls.thermal;
        let [lo, hi] = self.admit_equilibrium_solve(heat_input_j, controls)?;
        // Re-admit the physical response under this call's exact controls; no
        // new equilibrium solve or numerical promotion of the candidate occurs.
        let admitted = self.equilibrium_caloric_state(
            initial.candidate.clone(),
            controls.response,
            &mut cancelled,
        )?;
        if heat_input_j == 0.0 && (lo..=hi).contains(&initial.temperature().value()) {
            poll(&mut cancelled)?;
            return Ok(MulticomponentEquilibriumCaloricBalanceV1 {
                final_state: admitted,
                heat_input_j,
                enthalpy_change_j: 0.0,
                entropy_change_j_per_k: 0.0,
                residual_j: 0.0,
                temperature_bracket: [initial.temperature(); 2],
                iterations: 0,
                completed_phase_solves: 0,
            });
        }
        let (root, completed_phase_solves) =
            self.solve_equilibrium_heat(admitted, controls, &mut cancelled, |state, cancelled| {
                let change_j = self.equilibrium_caloric_change(initial, state, false, cancelled)?;
                let residual_j =
                    checked(change_j - heat_input_j, "local equilibrium heat residual")?;
                Ok((change_j, residual_j))
            })?;
        let entropy_change_j_per_k =
            self.equilibrium_caloric_change(initial, &root.sample.state, true, &mut cancelled)?;
        let delta_t = root.sample.state.temperature().value() - initial.temperature().value();
        if heat_input_j != 0.0
            && (delta_t == 0.0
                || delta_t.signum() != heat_input_j.signum()
                || root.sample.change_j == 0.0
                || root.sample.change_j.signum() != heat_input_j.signum()
                || entropy_change_j_per_k == 0.0
                || entropy_change_j_per_k.signum() != heat_input_j.signum())
        {
            return Err(MulticomponentEquilibriumCaloricErrorV1::Unrepresentable(
                "nonzero heat does not resolve the positive local temperature/entropy response",
            ));
        }
        admit_clausius_change(
            initial,
            &root.sample.state,
            root.sample.change_j,
            entropy_change_j_per_k,
            thermal,
        )?;
        poll(&mut cancelled)?;
        Ok(MulticomponentEquilibriumCaloricBalanceV1 {
            final_state: root.sample.state,
            heat_input_j,
            enthalpy_change_j: root.sample.change_j,
            entropy_change_j_per_k,
            residual_j: root.sample.residual_j,
            temperature_bracket: root.bracket.map(Temperature::new),
            iterations: root.iterations,
            completed_phase_solves,
        })
    }

    /// Direct H change of the same inventory and phase count; original pure
    /// datums cancel before arithmetic. Phase amounts AND compositions change.
    pub fn equilibrium_enthalpy_change_j(
        &self,
        initial: &MulticomponentEquilibriumCaloricStateV1,
        final_state: &MulticomponentEquilibriumCaloricStateV1,
    ) -> Result<f64> {
        self.equilibrium_caloric_change(initial, final_state, false, &mut || false)
    }

    /// Direct Cp/T and ideal/nonideal phase entropy change for the same feed.
    /// This does not infer a reservoir temperature or entropy production.
    pub fn equilibrium_entropy_change_j_per_k(
        &self,
        initial: &MulticomponentEquilibriumCaloricStateV1,
        final_state: &MulticomponentEquilibriumCaloricStateV1,
    ) -> Result<f64> {
        self.equilibrium_caloric_change(initial, final_state, true, &mut || false)
    }

    fn admit_equilibrium_caloric(
        &self,
        state: &MulticomponentEquilibriumCaloricStateV1,
    ) -> Result<()> {
        if !Arc::ptr_eq(&self.owner, &state.owner) && self.owner != state.owner {
            return Err(MulticomponentEquilibriumCaloricErrorV1::Invalid(
                "equilibrium parcel requires its original complete caloric owner",
            ));
        }
        Ok(())
    }

    fn admit_equilibrium_solve(
        &self,
        heat_input_j: f64,
        controls: MulticomponentEquilibriumCaloricSolveV1,
    ) -> Result<[f64; 2]> {
        let thermal = controls.thermal;
        let [lo, hi] = thermal.temperature_bracket.map(|t| t.value());
        if !heat_input_j.is_finite()
            || !lo.is_finite()
            || !hi.is_finite()
            || lo >= hi
            || lo < self.temperature_domain()[0].value()
            || hi > self.temperature_domain()[1].value()
            || !positive(thermal.temperature_tolerance_k)
            || !nonnegative(thermal.energy_tolerance_j)
            || !nonnegative(thermal.entropy_tolerance_j_per_k)
            || !(1..=256).contains(&thermal.maximum_iterations)
            || !positive(controls.maximum_temperature_step_k)
            || !(1..=65_536).contains(&controls.maximum_phase_solves)
        {
            return Err(MulticomponentEquilibriumCaloricErrorV1::Invalid(
                "equilibrium thermal controls or finite signed heat input",
            ));
        }
        Ok([lo, hi])
    }

    // ONE continuation and thermal inverse for heating, relaxation and mixing.
    // Callers retain their own physical input account. In particular a small Q
    // need not be rounded into a large artificial heat offset before inversion.
    fn solve_equilibrium_heat<C: FnMut() -> bool>(
        &self,
        admitted: MulticomponentEquilibriumCaloricStateV1,
        controls: MulticomponentEquilibriumCaloricSolveV1,
        cancelled: &mut C,
        mut account: impl FnMut(&MulticomponentEquilibriumCaloricStateV1, &mut C) -> Result<(f64, f64)>,
    ) -> Result<(
        temperature::TemperatureRoot<MulticomponentEquilibriumCaloricStateV1>,
        usize,
    )> {
        let thermal = controls.thermal;
        let mut continuation = Some(admitted);
        let mut completed_phase_solves = 0;
        let root = temperature::solve(
            TemperatureControl {
                bracket: thermal.temperature_bracket.map(|t| t.value()),
                temperature_tolerance_k: thermal.temperature_tolerance_k,
                energy_tolerance_j: thermal.energy_tolerance_j,
                maximum_iterations: thermal.maximum_iterations,
            },
            |value| -> Result<_> {
                poll(cancelled)?;
                let state = self.continue_equilibrium_temperature(
                    continuation
                        .take()
                        .expect("last complete continuation state"),
                    Temperature::new(value),
                    controls,
                    &mut completed_phase_solves,
                    cancelled,
                )?;
                let (change_j, residual_j) = account(&state, cancelled)?;
                continuation = Some(state.clone());
                Ok(TemperatureSample {
                    state,
                    change_j,
                    residual_j,
                })
            },
            |state| Some(state.isobaric_heat_capacity_j_per_k()),
        )?;
        Ok((root, completed_phase_solves))
    }

    fn equilibrium_caloric_change(
        &self,
        initial: &MulticomponentEquilibriumCaloricStateV1,
        final_state: &MulticomponentEquilibriumCaloricStateV1,
        entropy: bool,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<f64> {
        self.admit_equilibrium_caloric(initial)?;
        self.admit_equilibrium_caloric(final_state)?;
        if initial.component_moles() != final_state.component_moles()
            || initial.candidate.phases().len() != final_state.candidate.phases().len()
        {
            return Err(MulticomponentEquilibriumCaloricErrorV1::Invalid(
                "local caloric change requires the same component inventory and phase count",
            ));
        }
        let mut change = Sum::default();
        for (&amount, source) in initial.component_moles().iter().zip(self.heat_capacities()) {
            poll(cancelled)?;
            if amount != 0.0 {
                let increment = if entropy {
                    source.entropy_change_j_per_mol_k(
                        initial.temperature(),
                        final_state.temperature(),
                    )?
                } else {
                    source.enthalpy_change_j_per_mol(
                        initial.temperature(),
                        final_state.temperature(),
                    )?
                };
                change.add(scaled(
                    amount,
                    increment,
                    "equilibrium pure-source increment",
                )?);
            }
        }
        let n = initial.component_moles().len();
        let coefficients = if entropy {
            self.thermodynamics().chi_constants()
        } else {
            self.thermodynamics().chi_inverse_temperature_k()
        };
        let sign = if entropy { -1.0 } else { 1.0 };
        for (old, new) in initial
            .candidate
            .phases()
            .iter()
            .zip(final_state.candidate.phases())
        {
            let before = old.thermodynamics().mole_fractions();
            let after = new.thermodynamics().mole_fractions();
            let delta_amount = new.total_moles() - old.total_moles();
            let middle_amount = 0.5 * new.total_moles() + 0.5 * old.total_moles();
            for i in 0..n {
                poll(cancelled)?;
                if entropy && before[i] != 0.0 {
                    // Symmetric difference of n_i*ln(x_i), with a complement
                    // logarithm for a dominant fraction beside a real trace.
                    let delta_n = new.component_moles()[i] - old.component_moles()[i];
                    let middle_n = 0.5 * new.component_moles()[i] + 0.5 * old.component_moles()[i];
                    let middle_log = 0.5 * fraction_log(before, i) + 0.5 * fraction_log(after, i);
                    change.add(scaled(
                        delta_n,
                        -R * middle_log,
                        "phase amount ideal entropy",
                    )?);
                    change.add(scaled(
                        middle_n,
                        -R * fraction_log_change(before, after, i),
                        "phase composition ideal entropy",
                    )?);
                }
                for j in i + 1..n {
                    let coefficient = coefficients[i * n + j];
                    let middle_product = 0.5 * ((coefficient * before[i]) * before[j])
                        + 0.5 * ((coefficient * after[i]) * after[j]);
                    let product_change = (coefficient * (after[i] - before[i]))
                        * (0.5 * after[j] + 0.5 * before[j])
                        + (coefficient * (after[j] - before[j]))
                            * (0.5 * after[i] + 0.5 * before[i]);
                    change.add(scaled(
                        delta_amount,
                        sign * R * middle_product,
                        "phase amount excess increment",
                    )?);
                    change.add(scaled(
                        middle_amount,
                        sign * R * product_change,
                        "phase composition excess increment",
                    )?);
                }
            }
        }
        Ok(checked(change.value(), "equilibrium caloric increment")?)
    }

    fn continue_equilibrium_temperature(
        &self,
        mut state: MulticomponentEquilibriumCaloricStateV1,
        target: Temperature,
        controls: MulticomponentEquilibriumCaloricSolveV1,
        completed: &mut usize,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<MulticomponentEquilibriumCaloricStateV1> {
        while state.temperature() != target {
            poll(cancelled)?;
            if *completed >= controls.maximum_phase_solves {
                return Err(MulticomponentEquilibriumCaloricErrorV1::PhaseSolveBudget {
                    completed_phase_solves: *completed,
                });
            }
            let remaining = target.value() - state.temperature().value();
            let mut distance = remaining.abs().min(controls.maximum_temperature_step_k);
            for (phase, derivatives) in state.candidate.phases().iter().zip(
                state
                    .phase_component_derivatives
                    .chunks_exact(state.component_moles().len()),
            ) {
                poll(cancelled)?;
                for (&amount, &derivative) in phase.component_moles().iter().zip(derivatives) {
                    if derivative != 0.0 {
                        distance = distance.min(0.25 * (amount / derivative.abs()));
                    }
                }
            }
            let next_temperature = if distance == remaining.abs() {
                target
            } else {
                Temperature::new(state.temperature().value() + remaining.signum() * distance)
            };
            if !positive(distance) || next_temperature == state.temperature() {
                return Err(MulticomponentEquilibriumCaloricErrorV1::Unrepresentable(
                    "local phase predictor cannot resolve a positive temperature step",
                ));
            }
            let step = next_temperature.value() - state.temperature().value();
            let seed = predicted_seed(&state, step, cancelled)?;
            let candidate = self.thermodynamics().equilibrium_candidate(
                next_temperature,
                state.component_moles(),
                &seed,
                state.candidate.controls(),
                &mut *cancelled,
            )?;
            *completed += 1;
            let next =
                self.equilibrium_caloric_state(candidate, controls.response, &mut *cancelled)?;
            // Bounded prediction does not prove branch identity. Refuse a
            // large corrected jump rather than silently relabelling phases.
            for (old, new) in state.candidate.phases().iter().zip(next.candidate.phases()) {
                poll(cancelled)?;
                if old
                    .component_moles()
                    .iter()
                    .zip(new.component_moles())
                    .any(|(&a, &b)| a > 0.0 && ((b - a) / a).abs() > 0.5)
                {
                    return Err(MulticomponentEquilibriumCaloricErrorV1::Unrepresentable(
                        "corrected candidate left the local phase predictor neighborhood",
                    ));
                }
            }
            let increment = self.equilibrium_caloric_change(&state, &next, false, cancelled)?;
            if increment == 0.0 || increment.signum() != step.signum() {
                return Err(MulticomponentEquilibriumCaloricErrorV1::Unrepresentable(
                    "continued branch does not resolve positive isobaric heat capacity",
                ));
            }
            state = next;
        }
        poll(cancelled)?;
        Ok(state)
    }
}

fn admit_clausius_change(
    initial: &MulticomponentEquilibriumCaloricStateV1,
    final_state: &MulticomponentEquilibriumCaloricStateV1,
    enthalpy_change_j: f64,
    entropy_change_j_per_k: f64,
    controls: MulticomponentCaloricSolveV1,
) -> Result<()> {
    let endpoint_entropy = [
        checked(
            enthalpy_change_j / initial.temperature().value(),
            "initial Clausius bound",
        )?,
        checked(
            enthalpy_change_j / final_state.temperature().value(),
            "final Clausius bound",
        )?,
    ];
    let bounds_j_per_k = [
        endpoint_entropy[0].min(endpoint_entropy[1]),
        endpoint_entropy[0].max(endpoint_entropy[1]),
    ];
    let scale = entropy_change_j_per_k
        .abs()
        .max(bounds_j_per_k[0].abs())
        .max(bounds_j_per_k[1].abs());
    let tolerance_j_per_k = checked(
        controls.entropy_tolerance_j_per_k
            + controls.energy_tolerance_j
                / initial
                    .temperature()
                    .value()
                    .min(final_state.temperature().value())
            + (64.0 * f64::EPSILON) * scale,
        "finite Clausius allowance",
    )?;
    if bounds_j_per_k[0] - entropy_change_j_per_k > tolerance_j_per_k
        || entropy_change_j_per_k - bounds_j_per_k[1] > tolerance_j_per_k
    {
        return Err(MulticomponentEquilibriumCaloricErrorV1::ClausiusBounds {
            entropy_change_j_per_k,
            bounds_j_per_k,
            tolerance_j_per_k,
        });
    }
    Ok(())
}

fn predicted_seed(
    state: &MulticomponentEquilibriumCaloricStateV1,
    step_k: f64,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<MulticomponentEquilibriumSeedV1> {
    let n = state.component_moles().len();
    let p = state.candidate.phases().len();
    let mut amounts = vec![vec![0.0; n]; p];
    for (i, &feed) in state.component_moles().iter().enumerate() {
        poll(cancelled)?;
        if feed == 0.0 {
            continue;
        }
        let anchor = (0..p)
            .max_by(|&a, &b| {
                state.candidate.phases()[a].component_moles()[i]
                    .total_cmp(&state.candidate.phases()[b].component_moles()[i])
                    .then_with(|| b.cmp(&a))
            })
            .expect("accepted positive phase count");
        let mut others = Sum::default();
        for (phase, row) in amounts.iter_mut().enumerate() {
            if phase == anchor {
                continue;
            }
            let value = state.candidate.phases()[phase].component_moles()[i]
                + step_k * state.phase_component_derivatives[phase * n + i];
            if !positive(value) {
                return Err(MulticomponentEquilibriumCaloricErrorV1::Unrepresentable(
                    "local phase predictor loses a positive component",
                ));
            }
            row[i] = value;
            others.add(value);
        }
        amounts[anchor][i] = feed - others.value();
        if !positive(amounts[anchor][i]) {
            return Err(MulticomponentEquilibriumCaloricErrorV1::Unrepresentable(
                "local phase predictor cannot conserve its largest component allocation",
            ));
        }
    }
    let mut weights = Vec::with_capacity(p);
    for row in &mut amounts {
        poll(cancelled)?;
        let mut sum = Sum::default();
        for &amount in row.iter() {
            sum.add(amount);
        }
        let total = sum.value();
        if !positive(total) {
            return Err(MulticomponentEquilibriumCaloricErrorV1::Unrepresentable(
                "local phase predictor amount",
            ));
        }
        weights.push(total);
        for amount in row {
            let fraction = *amount / total;
            if *amount > 0.0 && !positive(fraction) {
                return Err(MulticomponentEquilibriumCaloricErrorV1::Unrepresentable(
                    "local phase predictor loses a real trace in normalization",
                ));
            }
            *amount = fraction;
        }
    }
    Ok(MulticomponentEquilibriumSeedV1 {
        phase_mole_fractions: amounts,
        phase_weights: weights,
    })
}

fn other_fractions(fractions: &[f64], index: usize) -> f64 {
    let mut others = Sum::default();
    for (i, &value) in fractions.iter().enumerate() {
        if i != index {
            others.add(value);
        }
    }
    others.value()
}

fn fraction_log(fractions: &[f64], index: usize) -> f64 {
    if fractions[index] > 0.5 {
        ln_complement(other_fractions(fractions, index))
    } else {
        ln(fractions[index])
    }
}

fn fraction_log_change(before: &[f64], after: &[f64], index: usize) -> f64 {
    if before[index] > 0.5 && after[index] > 0.5 {
        let old = other_fractions(before, index);
        let new = other_fractions(after, index);
        if old == new {
            return 0.0;
        }
        let (lo, hi, sign) = if new < old {
            (new, old, 1.0)
        } else {
            (old, new, -1.0)
        };
        -sign * ln_complement((hi - lo) / (1.0 - lo))
    } else {
        log_ratio(after[index], before[index])
    }
}

fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<()> {
    if cancelled() {
        Err(MulticomponentEquilibriumCaloricErrorV1::Cancelled)
    } else {
        Ok(())
    }
}

#[cfg(test)]
mod tests;
