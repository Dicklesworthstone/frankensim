//! Prescribed exterior entropy in the original atomic reacting composition.

use super::*;
use crate::multicomponent::{ExteriorMixtureEntropyReport, MixtureGas};
use crate::multicomponent::dissipation::ReservoirMolecularEntropyReport;
use crate::multicomponent::dissipation::entropy::ExteriorTransportInventory;
pub use crate::multicomponent::dissipation::entropy::reactive::ReactiveExteriorEntropyReport;

/// Original-source entropy receipts of both complete molecular bath intervals.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SplitExteriorMolecularEntropyReport {
    /// Molecular evolution before chemistry and inviscid transport.
    pub first: ReservoirMolecularEntropyReport,
    /// Molecular evolution after chemistry and inviscid transport.
    pub second: ReservoirMolecularEntropyReport,
}

/// Complete reacting evolution with actual prescribed-exterior entropy supply.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ReactiveExteriorEntropyStepReport {
    /// Original chemistry, transport, common duration and conservative balances.
    pub coupled: ReactiveStepReport,
    /// Direct original-to-final source entropy minus actual boundary supply,
    /// with actual chemical species changes retained in the source account.
    pub entropy: ReactiveExteriorEntropyReport,
    /// The middle inviscid update in its own initial-state entropy basis.
    /// This component receipt does not replace complete-composition admission.
    pub transport_entropy: ExteriorMixtureEntropyReport,
    /// Actual complete molecular halves when explicitly selected. These
    /// component receipts do not replace the original-to-final entropy gate.
    pub dissipation_entropy: Option<SplitExteriorMolecularEntropyReport>,
}

impl ReactiveExteriorEntropyStepReport {
    /// Closed entry points keep their original algorithms and exact reports.
    /// Only their report shape changes; there is no external entropy input.
    pub(super) fn from_closed(closed: ReactiveEntropyStepReport) -> Self {
        Self {
            coupled: closed.coupled,
            entropy: ReactiveExteriorEntropyReport {
                source_entropy_change_j_k: closed.entropy.source_entropy_change_j_k,
                boundary_entropy_input_j_k: [0.0; 2],
                entropy_generation_j_k: closed.entropy.entropy_change_j_k,
                source_entropy_generation_j_k: closed.entropy.source_entropy_change_j_k,
                reference_entropy_residual_j_k: closed.entropy.reference_entropy_residual_j_k,
                maximum_relative_species_residual: closed.entropy.maximum_relative_species_residual,
                rejected_entropy_steps: closed.entropy.rejected_entropy_steps,
            },
            transport_entropy: ExteriorMixtureEntropyReport {
                source_entropy_change_j_k: closed.transport_entropy.source_entropy_change_j_k,
                boundary_entropy_input_j_k: [0.0; 2],
                entropy_generation_j_k: closed.transport_entropy.entropy_change_j_k,
                source_entropy_generation_j_k: closed.transport_entropy.source_entropy_change_j_k,
                reference_entropy_residual_j_k: closed.transport_entropy.reference_entropy_residual_j_k,
                minimum_step_entropy_generation_j_k: closed.transport_entropy.minimum_step_entropy_change_j_k,
                maximum_relative_species_residual: closed.transport_entropy.maximum_relative_species_residual,
                rejected_entropy_steps: closed.transport_entropy.rejected_entropy_steps,
            },
            dissipation_entropy: closed.dissipation_entropy.map(|halves| {
                SplitExteriorMolecularEntropyReport {
                    first: closed_molecular_report(halves.first),
                    second: closed_molecular_report(halves.second),
                }
            }),
        }
    }
}

fn closed_molecular_report(report: MixtureEntropyReport) -> ReservoirMolecularEntropyReport {
    ReservoirMolecularEntropyReport {
        source_entropy_change_j_k: report.source_entropy_change_j_k,
        boundary_entropy_input_j_k: [0.0; 2],
        entropy_generation_j_k: report.entropy_change_j_k,
        source_entropy_generation_j_k: report.source_entropy_change_j_k,
        reference_entropy_residual_j_k: report.reference_entropy_residual_j_k,
        minimum_step_entropy_generation_j_k: report.minimum_step_entropy_change_j_k,
        maximum_relative_species_residual: report.maximum_relative_species_residual,
        rejected_entropy_steps: report.rejected_entropy_steps,
    }
}

impl ReactingEuler1d {
    /// Advance original fixed-volume chemistry around entropy-consistent
    /// inviscid transport through prescribed exterior gas states.
    ///
    /// The native C(dt/2), A(dt), C(dt/2) owner retains its reactor, accepted and
    /// trial workspaces, selected chemistry method and bounded restart policy.
    /// Every active exterior is fixed throughout both flow RK stages and every
    /// whole-composition retry, and may have nonzero axial velocity. A shorter
    /// middle flow step restarts BOTH chemistry halves at the shorter duration.
    ///
    /// Each Some replaces that end's reflecting or transmissive declaration;
    /// a None physical end must be a stationary reflecting wall. Both None
    /// execute the original closed entropy composition exactly. Exterior cells
    /// must use the original gas bank's canonical species order and referenced
    /// energy. They are prescribed Riemann traces, not the stationary stagnation
    /// conditions of the legacy compact-reservoir characteristic closure.
    ///
    /// The complete account starts BEFORE chemistry. It evaluates the actual
    /// final source-state entropy, observes the actual per-cell chemical species
    /// changes, and contracts the actual accepted SSP end transfers against
    /// the prescribed exterior entropy efforts in that SAME immutable basis.
    /// Species, momentum, kinetic energy and advective entropy potential enter
    /// the boundary supply. Summed component entropy reports are not used.
    ///
    /// Only the source-reference contribution of numerical species continuity
    /// drift is removed. Each drift is bounded independently against original
    /// inventory plus absolute chemical and boundary throughput. Real chemical
    /// production and physical entropy transported by each species are retained.
    /// The original reactor keeps its own source, energy, positivity, entropy
    /// and work checks; no additional reaction heat is introduced.
    ///
    /// First-order space and fixed constant-area geometry are required. Common
    /// positive species support, or common absence, is required on every
    /// transport face. No trace seeding, clipping or flux fallback is performed.
    /// This entry includes no molecular halves or finite exterior storage.
    /// The caller's entropy allowance applies separately to the middle flow
    /// and to the complete actual-state composition. Neither admission nor
    /// symmetric splitting proves unconditional positivity or temporal order.
    /// Refusal or cancellation leaves the accepted field and clock unchanged.
    pub fn step_entropy_consistent_with_exterior_states(
        &mut self,
        maximum_duration: Time,
        controls: ReactiveStepControl,
        entropy: MixtureEntropyControl,
        exterior_states: [Option<Cell<'_>>; 2],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<ReactiveExteriorEntropyStepReport, ReactiveError> {
        if exterior_states.iter().all(Option::is_none) {
            let closed = self.step_entropy_consistent(
                maximum_duration, controls, entropy, cancelled,
            )?;
            return Ok(ReactiveExteriorEntropyStepReport::from_closed(closed));
        }
        let mut account = CoupledEntropy::new_with_exterior_states(
            &self.accepted, controls.transport, entropy, exterior_states, &mut cancelled,
        )?;
        let (coupled, _) = self.step_inner::<ConstantTransport>(
            maximum_duration,
            controls,
            None,
            stabilized::MolecularMethod::Explicit(DissipationControl::default()),
            exterior_states,
            None,
            Some(&mut account),
            &mut cancelled,
        )?;
        Ok(account.complete_exterior(coupled))
    }
}

/// One optional reusable receipt slot; failed or shortened attempts are cleared
/// before the original chemistry owner restarts. The private transport inventory
/// is constructed only by the native SSP boundary-quadrature owner.
pub(super) struct ExteriorCoupling {
    transport: Option<(ExteriorMixtureEntropyReport, ExteriorTransportInventory)>,
    completed: Option<(
        ReactiveExteriorEntropyReport,
        ExteriorMixtureEntropyReport,
        Option<SplitExteriorMolecularEntropyReport>,
    )>,
}

impl ExteriorCoupling {
    pub(super) fn reset(&mut self) {
        self.transport = None;
        self.completed = None;
    }
}

impl CoupledEntropy {
    pub(super) fn new_with_exterior_states(
        flow: &MixtureEuler1d,
        transport: StepControl,
        controls: MixtureEntropyControl,
        exterior_states: [Option<Cell<'_>>; 2],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, ReactiveError> {
        poll(cancelled)?;
        if transport.reconstruction != Reconstruction::FirstOrder {
            return Err(FlowError::EntropyDomain(
                "entropy-consistent reacting flow requires first-order reconstruction",
            ).into());
        }
        if flow.mesh.uniform_cross_section().is_none() {
            return Err(FlowError::EntropyDomain(
                "entropy-consistent reacting flow requires fixed constant-area geometry",
            ).into());
        }
        let Boundaries::Ends { left, right } = flow.boundaries else {
            return Err(FlowError::EntropyDomain(
                "prescribed reacting exterior states require physical ends",
            ).into());
        };
        for (end, boundary) in [left, right].into_iter().enumerate() {
            if !matches!(boundary, Boundary::Reflecting | Boundary::Transmissive) {
                return Err(FlowError::EntropyDomain(
                    "reacting exterior entropy requires reflecting or transmissive ends",
                ).into());
            }
            if exterior_states[end].is_none() && boundary != Boundary::Reflecting {
                return Err(FlowError::EntropyDomain(
                    "an open reacting end requires its exterior thermodynamic state",
                ).into());
            }
        }
        let account = ReactiveEntropyAccount::new_with_exterior_states(
            flow, controls, exterior_states, cancelled,
        )?;
        let mut exterior = Vec::new();
        exterior.try_reserve_exact(1).map_err(|_| FlowError::Allocation)?;
        exterior.push(ExteriorCoupling { transport: None, completed: None });
        poll(cancelled)?;
        Ok(Self { controls, account, completed: None, exterior })
    }

    pub(in crate::multicomponent::reactive) fn has_exterior_states(&self) -> bool {
        !self.exterior.is_empty()
    }

    pub(in crate::multicomponent::reactive) fn advance_transport(
        &mut self,
        flow: &mut MixtureEuler1d,
        duration: Time,
        controls: StepControl,
        exterior_states: [Option<Cell<'_>>; 2],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(TransportStepReport, Option<MixtureEntropyReport>), ReactiveError> {
        if let Some(exterior) = self.exterior.first_mut() {
            let (transport, report, inventory) =
                flow.step_entropy_consistent_with_exterior_inventory(
                    duration, controls, self.controls, exterior_states, cancelled,
                )?;
            exterior.transport = Some((report, inventory));
            Ok((transport, None))
        } else {
            let (transport, report) = flow.step_entropy_consistent(
                duration, controls, self.controls, cancelled,
            )?;
            Ok((transport, Some(report)))
        }
    }

    /// Observe one complete accepted molecular half without retaining a second
    /// opaque stage inventory. The whole source account owns transactional
    /// accumulation in its immutable pre-split basis.
    pub(in crate::multicomponent::reactive) fn observe_molecular_inventory(
        &mut self,
        gas: &MixtureGas,
        inventory: &ExteriorTransportInventory,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), ReactiveError> {
        self.account.observe_molecular_inventory(gas, inventory, cancelled)?;
        Ok(())
    }

    pub(in crate::multicomponent::reactive) fn check_exterior_candidate(
        &mut self,
        initial: &MixtureEuler1d,
        candidate: &MixtureEuler1d,
        temperature_bounds: Option<[f64; 2]>,
        molecular: [Option<ReservoirMolecularEntropyReport>; 2],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), ReactiveError> {
        let dissipation_entropy = match molecular {
            [None, None] => None,
            [Some(first), Some(second)] => Some(SplitExteriorMolecularEntropyReport {
                first, second,
            }),
            _ => return Err(ReactiveError::Invalid(
                "incomplete reacting exterior molecular entropy receipts",
            )),
        };
        let exterior = self.exterior.first_mut().ok_or(ReactiveError::Invalid(
            "missing reacting exterior entropy account",
        ))?;
        let (transport, inventory) = exterior.transport.as_ref().ok_or(ReactiveError::Invalid(
            "missing reacting exterior transport inventory",
        ))?;
        let entropy = self.account.check_candidate_with_exterior(
            initial, candidate, temperature_bounds, inventory, cancelled,
        )?;
        exterior.completed = Some((entropy, *transport, dissipation_entropy));
        Ok(())
    }

    /// Every fallible contraction and report reduction precedes native swap.
    pub(super) fn complete_exterior(self, coupled: ReactiveStepReport) -> ReactiveExteriorEntropyStepReport {
        let (entropy, transport_entropy, dissipation_entropy) = self.exterior.into_iter().next()
            .and_then(|exterior| exterior.completed)
            .expect("accepted reacting exterior step retained complete entropy admission");
        ReactiveExteriorEntropyStepReport { coupled, entropy, transport_entropy, dissipation_entropy }
    }
}
