//! Continuing simultaneous chemical/pressure/heat feedback on the existing
//! material assembly. Move the original accepted candidates, preserving every
//! gas species/source, circuit state, thermal event cursor and port history.

use super::super::exchange::reactive::{
    CoupledMaterialExchangeCandidate, CoupledMaterialExchangeError,
};
use super::*;
use crate::cavity_phs::piston::thermal::boundary::ports::network::reactive::{
    CoupledChamberReaction, CoupledChamberStep,
};
use crate::circuit_enthalpy::CircuitEnthalpyError;

/// The exact native report after its complete physical endpoint moves into the
/// assembly. Ordinary midpoint and reacting endpoint accounts remain distinct.
pub struct CoupledExchangeChamberStep<I = ResolvedMaterialStatePoint> {
    pub id: u64,
    pub surfaces: Vec<ThermalPortSample<(Contact, I)>>,
    pub boundary: CoupledChamberStep,
}

/// Actual material/contact/electrical sources and all same-interval chamber
/// reports. Reaction progress lives only inside its original boundary account.
pub struct CoupledExchangeAdvanceSources<I = ResolvedMaterialStatePoint> {
    pub material: ContactThermalSample<MaterialExchangeSample>,
    pub chambers: Vec<CoupledExchangeChamberStep<I>>,
    pub transferred_heat_j: f64,
}

/// Complete original finite-network and all-owner balances for an accepted
/// interval. This is a report of the state advance, not a parallel history.
pub struct CoupledMaterialExchangeAdvance<I = ResolvedMaterialStatePoint> {
    pub report: CircuitEnthalpyReport<CoupledExchangeAdvanceSources<I>>,
    pub step: NetworkStep,
}

fn admission_error<E>(error: StateAdmissionError) -> CoupledMaterialExchangeError<E> {
    match error {
        StateAdmissionError::Invalid(message) => CircuitEnthalpyError::Invalid(message),
        StateAdmissionError::Cancelled => CircuitEnthalpyError::Cancelled,
    }
}

impl<'m> MaterialExchangeState<'m> {
    /// Preview one complete same-endpoint reactive assembly interval.
    ///
    /// The thermal timeline supplies the only duration. Every live chamber
    /// needs exactly one drive; optional reaction bindings select the implicit
    /// policy without replacing its accepted species or pressure history.
    /// The original prepared network resolves chemistry, gas-wall heat,
    /// electrical/contact losses, radiation and finite-cell enthalpy together.
    /// Each resulting gas/volume/velocity and circuit/thermal endpoint moves
    /// into the returned state, with its exact native reports retained.
    ///
    /// An empty reaction list preserves ordinary chamber integration. Empty
    /// chambers retain electrical/thermal evolution. All fixed topology,
    /// source-domain, independent-wall and first-order limitations of the
    /// prepared candidate remain; no flow or chemistry-only subcycling is
    /// inserted here. Any refusal leaves the live owner entirely unchanged.
    #[allow(clippy::too_many_arguments)]
    pub fn coupled_reactive_candidate(
        &self,
        prepared: &PreparedMaterialExchangeNetwork<'_, '_>,
        cx: &Cx<'_>,
        target_s: f64,
        ports: &[f64],
        drives: &[ChamberDrive],
        reactions: &[CoupledChamberReaction<'_>],
        gas_patches: &[MaterialContact],
        controls: NetworkSolve,
        guess_j: &[f64],
    ) -> Result<(Self, CoupledMaterialExchangeAdvance), CoupledMaterialExchangeError> {
        let chambers = self
            .driven_chambers(&self.chambers, cx, target_s, drives, controls)
            .map_err(admission_error)?;
        let proposal = prepared.coupled_reactive_step_candidate(
            cx,
            &self.thermal,
            target_s,
            &self.circuit,
            ports,
            &chambers,
            reactions,
            gas_patches,
            controls,
            guess_j,
        )?;
        self.collect_coupled_candidate(cx, proposal)
    }

    /// Advance all continuing owners atomically after the final checkpoint.
    /// Chemical progress and contact/electrical heats are already included;
    /// consuming the returned records must not apply any of them again.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_coupled_reactive(
        &mut self,
        prepared: &PreparedMaterialExchangeNetwork<'_, '_>,
        cx: &Cx<'_>,
        target_s: f64,
        ports: &[f64],
        drives: &[ChamberDrive],
        reactions: &[CoupledChamberReaction<'_>],
        gas_patches: &[MaterialContact],
        controls: NetworkSolve,
        guess_j: &[f64],
    ) -> Result<CoupledMaterialExchangeAdvance, CoupledMaterialExchangeError> {
        let (next, report) = self.coupled_reactive_candidate(
            prepared,
            cx,
            target_s,
            ports,
            drives,
            reactions,
            gas_patches,
            controls,
            guess_j,
        )?;
        cx.checkpoint()
            .map_err(|_| CircuitEnthalpyError::Cancelled)?;
        *self = next;
        Ok(report)
    }

    /// Preview the SAME continuing assembly with a caller-selected gas heat law.
    ///
    /// Reacting contacts see the actual gas endpoint, ordinary contacts their
    /// original gas-energy midpoint, and both see the same proposed thermal
    /// cell/material snapshot. The callback returns a conductance and its
    /// original source; it must be repeatable and leave accepted state untouched.
    /// Surface admission, receiving-cell debit and all-owner balance remain
    /// owned by the original prepared network, with no second thermal solve.
    #[allow(clippy::too_many_arguments, clippy::type_complexity)]
    pub fn coupled_reactive_with_gas_law_candidate<I, E>(
        &self,
        prepared: &PreparedMaterialExchangeNetwork<'_, '_>,
        cx: &Cx<'_>,
        target_s: f64,
        ports: &[f64],
        drives: &[ChamberDrive],
        reactions: &[CoupledChamberReaction<'_>],
        surfaces: &[GasSurface],
        controls: NetworkSolve,
        guess_j: &[f64],
        resolve: impl Fn(
            Contact,
            &DepositionSimulation<'m>,
            &MaterialExchangeSample,
            &ResolvedThermochemicalGasState,
        ) -> Result<ConductanceSample<I>, E>,
    ) -> Result<(Self, CoupledMaterialExchangeAdvance<I>), CoupledMaterialExchangeError<E>> {
        let chambers = self
            .driven_chambers(&self.chambers, cx, target_s, drives, controls)
            .map_err(admission_error)?;
        let proposal = prepared.coupled_reactive_with_gas_law_candidate(
            cx,
            &self.thermal,
            target_s,
            &self.circuit,
            ports,
            &chambers,
            reactions,
            surfaces,
            controls,
            guess_j,
            resolve,
        )?;
        self.collect_coupled_candidate(cx, proposal)
    }

    /// Publish every custom heat-law endpoint together. Failed source queries,
    /// chemical/pressure/heat iterations and late cancellation preserve the
    /// complete accepted thermal clock, circuit, species and port-tick prefix.
    #[allow(clippy::too_many_arguments, clippy::type_complexity)]
    pub fn advance_coupled_reactive_with_gas_law<I, E>(
        &mut self,
        prepared: &PreparedMaterialExchangeNetwork<'_, '_>,
        cx: &Cx<'_>,
        target_s: f64,
        ports: &[f64],
        drives: &[ChamberDrive],
        reactions: &[CoupledChamberReaction<'_>],
        surfaces: &[GasSurface],
        controls: NetworkSolve,
        guess_j: &[f64],
        resolve: impl Fn(
            Contact,
            &DepositionSimulation<'m>,
            &MaterialExchangeSample,
            &ResolvedThermochemicalGasState,
        ) -> Result<ConductanceSample<I>, E>,
    ) -> Result<CoupledMaterialExchangeAdvance<I>, CoupledMaterialExchangeError<E>> {
        let (next, report) = self.coupled_reactive_with_gas_law_candidate(
            prepared, cx, target_s, ports, drives, reactions, surfaces, controls, guess_j, resolve,
        )?;
        cx.checkpoint()
            .map_err(|_| CircuitEnthalpyError::Cancelled)?;
        *self = next;
        Ok(report)
    }

    // State assembly is shared by sourced conduction and custom heat laws.
    // It moves existing owner results and retains their distinct accounts; no
    // chemical, circuit, caloric or contact operation is rerun on consumption.
    fn collect_coupled_candidate<I, E>(
        &self,
        cx: &Cx<'_>,
        proposal: CoupledMaterialExchangeCandidate<'m, I>,
    ) -> Result<(Self, CoupledMaterialExchangeAdvance<I>), CoupledMaterialExchangeError<E>> {
        let (thermal, circuit, report, step) = proposal.into_parts();
        let mut next_chambers = BTreeMap::new();
        let mut accounts = Vec::with_capacity(report.source.chambers.len());
        for chamber in report.source.chambers {
            cx.checkpoint()
                .map_err(|_| CircuitEnthalpyError::Cancelled)?;
            let id = chamber.id;
            let (piston, surfaces, boundary) = chamber.boundary.into_parts();
            if !self.chambers.contains_key(&id) || next_chambers.insert(id, piston).is_some() {
                return Err(CircuitEnthalpyError::Invalid(
                    "coupled candidate chamber identities differ from live assembly",
                ));
            }
            accounts.push(CoupledExchangeChamberStep {
                id,
                surfaces,
                boundary,
            });
        }
        if next_chambers.len() != self.chambers.len() {
            return Err(CircuitEnthalpyError::Invalid(
                "coupled candidate omitted a live chamber",
            ));
        }
        let source = CoupledExchangeAdvanceSources {
            material: report.source.material,
            chambers: accounts,
            transferred_heat_j: report.source.transferred_heat_j,
        };
        let report = CircuitEnthalpyReport {
            source,
            electrical: report.electrical,
            thermal: report.thermal,
            nonlinear: report.nonlinear,
            maximum_cell_residual_j: report.maximum_cell_residual_j,
            circuit_trials: report.circuit_trials,
            deposited_joule_j: report.deposited_joule_j,
            non_electrical_heat_j: report.non_electrical_heat_j,
            total_energy_residual_j: report.total_energy_residual_j,
        };
        cx.checkpoint()
            .map_err(|_| CircuitEnthalpyError::Cancelled)?;
        Ok((
            Self {
                thermal,
                circuit,
                chambers: next_chambers,
            },
            CoupledMaterialExchangeAdvance { report, step },
        ))
    }
}
