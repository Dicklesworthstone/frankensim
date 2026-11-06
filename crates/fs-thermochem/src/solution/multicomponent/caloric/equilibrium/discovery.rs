//! Finite heat and mixing endpoints with bounded original phase discovery.
//!
//! Each temperature sample runs the SAME declared deterministic discovery
//! policy, then binds its selected candidate to the original pure-liquid Cp
//! owner and unregularized local response. The shared caloric inverse accepts
//! actual enthalpy increments, including every initial and final phase. Phase
//! count may change between samples without continuing a vanishing phase.
//!
//! These are bounded-search endpoints, not a globally proved Gibbs envelope,
//! a continuous phase trajectory or a nucleation law. An unresolved enthalpy
//! jump is not interpolated into an invented coexistence/latent-heat plateau.
//! Failed discovery, unresolved response and exhausted work return no endpoint.

use super::super::{MulticomponentCaloricErrorV1, MulticomponentCaloricSolveV1};
use super::transactions::MulticomponentCaloricInventoryRefV1;
use super::{
    MulticomponentCaloricResponseControlV1, MulticomponentEquilibriumCaloricBalanceV1,
    MulticomponentEquilibriumCaloricErrorV1, MulticomponentEquilibriumCaloricStateV1,
    MulticomponentLiquidCaloricModelV1, Result, Sum, Temperature, TemperatureControl,
    TemperatureSample, checked, poll, temperature,
};
use crate::solution::multicomponent::discovery::{
    MulticomponentPhaseDiscoveryControlV1, MulticomponentPhaseDiscoveryReportV1,
};

/// Thermal tolerances and one unchanged discovery policy for every temperature.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MulticomponentDiscoveredCaloricSolveV1 {
    pub thermal: MulticomponentCaloricSolveV1,
    pub discovery: MulticomponentPhaseDiscoveryControlV1,
    pub response: MulticomponentCaloricResponseControlV1,
    /// Whole-call attempted candidate solves, 1..=1_048_576. Before each
    /// temperature sample, the complete `discovery.maximum_seed_trials`
    /// allowance must fit; only actual attempts are charged afterward.
    pub maximum_candidate_solves: usize,
}

/// Actual thermal balance and the original final candidate's search evidence.
#[derive(Clone, Debug, PartialEq)]
pub struct MulticomponentDiscoveredCaloricBalanceV1 {
    balance: MulticomponentEquilibriumCaloricBalanceV1,
    discovery: MulticomponentPhaseDiscoveryReportV1,
    controls: MulticomponentDiscoveredCaloricSolveV1,
    completed_searches: usize,
}
impl MulticomponentDiscoveredCaloricBalanceV1 {
    pub const fn final_state(&self) -> &MulticomponentEquilibriumCaloricStateV1 {
        self.balance.final_state()
    }
    /// `completed_phase_solves` counts all candidate attempts, including
    /// rejected starts, across every sampled temperature.
    pub const fn balance(&self) -> &MulticomponentEquilibriumCaloricBalanceV1 {
        &self.balance
    }
    pub const fn discovery_report(&self) -> &MulticomponentPhaseDiscoveryReportV1 {
        &self.discovery
    }
    pub const fn controls(&self) -> MulticomponentDiscoveredCaloricSolveV1 {
        self.controls
    }
    pub const fn completed_searches(&self) -> usize {
        self.completed_searches
    }
}

impl MulticomponentLiquidCaloricModelV1 {
    /// Apply signed heat to an actual homogeneous or separated liquid parcel,
    /// allowing the original bounded discovery to select a different phase
    /// count at the endpoint. No phase anchor or homogeneous fallback is used.
    ///
    /// The selected candidate carries only the original finite search's local
    /// qualifications. Every sample must admit physical, unregularized Cp.
    /// At Q=0, actual input-to-final entropy may not decrease beyond tolerance.
    /// For nonzero Q, reservoir entropy remains the caller's separate account.
    pub fn add_heat_discovered<'a>(
        &self,
        initial: impl Into<MulticomponentCaloricInventoryRefV1<'a>>,
        heat_input_j: f64,
        controls: MulticomponentDiscoveredCaloricSolveV1,
        cancelled: impl FnMut() -> bool,
    ) -> Result<MulticomponentDiscoveredCaloricBalanceV1> {
        self.equilibrate_discovered(&[initial.into()], heat_input_j, controls, cancelled)
    }

    /// Mix actual homogeneous or separated parcels at possibly different
    /// temperatures, conserving every original component and its full H/S.
    /// Every temperature trial independently uses the declared bounded search;
    /// it does not reuse a different parcel as a source of heat or phase mass.
    pub fn mix_discovered<'a>(
        &self,
        left: impl Into<MulticomponentCaloricInventoryRefV1<'a>>,
        right: impl Into<MulticomponentCaloricInventoryRefV1<'a>>,
        heat_input_j: f64,
        controls: MulticomponentDiscoveredCaloricSolveV1,
        cancelled: impl FnMut() -> bool,
    ) -> Result<MulticomponentDiscoveredCaloricBalanceV1> {
        self.equilibrate_discovered(
            &[left.into(), right.into()],
            heat_input_j,
            controls,
            cancelled,
        )
    }

    fn equilibrate_discovered(
        &self,
        inputs: &[MulticomponentCaloricInventoryRefV1<'_>],
        heat_input_j: f64,
        controls: MulticomponentDiscoveredCaloricSolveV1,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<MulticomponentDiscoveredCaloricBalanceV1> {
        poll(&mut cancelled)?;
        let bracket = self.admit_solve(heat_input_j, controls.thermal)?;
        if !(1..=1_048_576).contains(&controls.maximum_candidate_solves)
            || controls.discovery.maximum_seed_trials > controls.maximum_candidate_solves
        {
            return Err(MulticomponentEquilibriumCaloricErrorV1::Invalid(
                "whole-call discovery candidate budget",
            ));
        }
        for &input in inputs {
            match input {
                MulticomponentCaloricInventoryRefV1::Homogeneous(state) => self.admit(state)?,
                MulticomponentCaloricInventoryRefV1::Equilibrium(state) => {
                    self.admit_equilibrium_caloric(state)?;
                }
            }
        }
        let amounts = self.combine_component_moles(inputs.iter().map(|x| x.component_moles()))?;
        let n = self.heat_capacities().len();
        let (_, fractions) = self.inventory(&amounts[..n])?;
        let feed = self
            .thermodynamics()
            .evaluate(inputs[0].temperature(), &fractions[..n], &mut cancelled)
            .map_err(MulticomponentCaloricErrorV1::from)?;
        // One common feed cancels analytically between all initial and final
        // phases. Its H/S mixing terms are T-independent for original A+B/T.
        let input_to_feed = self.mixing_changes(
            inputs.iter().flat_map(|input| input.phases()),
            feed.mole_fractions(),
            &mut cancelled,
        )?;
        let mut completed_candidate_solves = 0;
        let mut completed_searches = 0;
        let thermal = controls.thermal;
        let root = temperature::solve(
            TemperatureControl {
                bracket,
                temperature_tolerance_k: thermal.temperature_tolerance_k,
                energy_tolerance_j: thermal.energy_tolerance_j,
                maximum_iterations: thermal.maximum_iterations,
            },
            |value| -> Result<_> {
                poll(&mut cancelled)?;
                if controls.maximum_candidate_solves - completed_candidate_solves
                    < controls.discovery.maximum_seed_trials
                {
                    return Err(
                        MulticomponentEquilibriumCaloricErrorV1::DiscoverySolveBudget {
                            completed_candidate_solves,
                            required_trial_allowance: controls.discovery.maximum_seed_trials,
                        },
                    );
                }
                let discovered = self.thermodynamics().discover_equilibrium_candidate(
                    Temperature::new(value),
                    &amounts[..n],
                    controls.discovery,
                    &mut cancelled,
                )?;
                completed_candidate_solves += discovered.report().attempted_candidates();
                completed_searches += 1;
                let report = discovered.report().clone();
                let state = self.equilibrium_caloric_state(
                    discovered.into_candidate(),
                    controls.response,
                    &mut cancelled,
                )?;
                let output_to_feed = self.mixing_changes(
                    MulticomponentCaloricInventoryRefV1::from(&state).phases(),
                    feed.mole_fractions(),
                    &mut cancelled,
                )?;
                let mut change = Sum::default();
                change.add(
                    self.sensible_inventory_change(
                        inputs
                            .iter()
                            .map(|x| (x.temperature(), x.component_moles())),
                        state.temperature(),
                        false,
                        &mut cancelled,
                    )?,
                );
                change.add(input_to_feed[0]);
                change.add(-output_to_feed[0]);
                let change_j = checked(change.value(), "discovered phase enthalpy change")?;
                // Preserve a small heat request separately from potentially
                // large phase storage until the compensated residual forms.
                change.add(-heat_input_j);
                Ok(TemperatureSample {
                    state: (state, report),
                    change_j,
                    residual_j: checked(change.value(), "discovered phase heat residual")?,
                })
            },
            |state| Some(state.0.isobaric_heat_capacity_j_per_k()),
        )?;
        let (final_state, discovery) = root.sample.state;
        let output_to_feed = self.mixing_changes(
            MulticomponentCaloricInventoryRefV1::from(&final_state).phases(),
            feed.mole_fractions(),
            &mut cancelled,
        )?;
        let mut entropy = Sum::default();
        entropy.add(
            self.sensible_inventory_change(
                inputs
                    .iter()
                    .map(|x| (x.temperature(), x.component_moles())),
                final_state.temperature(),
                true,
                &mut cancelled,
            )?,
        );
        entropy.add(input_to_feed[1]);
        entropy.add(-output_to_feed[1]);
        let entropy_change_j_per_k = checked(entropy.value(), "discovered phase entropy change")?;
        if heat_input_j == 0.0 && entropy_change_j_per_k < -thermal.entropy_tolerance_j_per_k {
            return Err(MulticomponentCaloricErrorV1::EntropyDecrease {
                change_j_per_k: entropy_change_j_per_k,
            }
            .into());
        }
        poll(&mut cancelled)?;
        Ok(MulticomponentDiscoveredCaloricBalanceV1 {
            balance: MulticomponentEquilibriumCaloricBalanceV1 {
                final_state,
                heat_input_j,
                enthalpy_change_j: root.sample.change_j,
                entropy_change_j_per_k,
                residual_j: root.sample.residual_j,
                temperature_bracket: root.bracket.map(Temperature::new),
                iterations: root.iterations,
                completed_phase_solves: completed_candidate_solves,
            },
            discovery,
            controls,
            completed_searches,
        })
    }
}
