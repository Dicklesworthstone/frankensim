//! Closed original-source entropy admission for the native reacting composition.
//!
//! Every stage remains with its original owner. The whole-step account observes
//! actual chemical species changes and measures the original-to-final source
//! entropy before publication. Genuine chemical entropy is retained; only the
//! arbitrary-reference effect of numerical species-source drift is separated.

use super::*;
use crate::multicomponent::dissipation::entropy::reactive::ReactiveEntropyAccount;
pub use crate::multicomponent::dissipation::entropy::reactive::ReactiveEntropyReport;
use crate::multicomponent::{
    Boundaries, Boundary, MixtureEntropyControl, MixtureEntropyReport, Reconstruction,
};

/// Original-source entropy receipts of the two complete molecular half steps.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SplitMolecularEntropyReport {
    /// Molecular evolution before chemistry and advection.
    pub first: MixtureEntropyReport,
    /// Molecular evolution after chemistry and advection.
    pub second: MixtureEntropyReport,
}

/// One accepted closed reacting step with component and whole-state entropy.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ReactiveEntropyStepReport {
    /// Native chemistry, physical transport, duration and conservative balances.
    pub coupled: ReactiveStepReport,
    /// Direct original-to-final entropy and actual chemical-source continuity.
    /// This is the admission account for the COMPLETE accepted composition.
    pub entropy: ReactiveEntropyReport,
    /// Actual middle inviscid step; excludes both chemical and molecular halves.
    pub transport_entropy: MixtureEntropyReport,
    /// Actual molecular halves when explicitly selected.
    pub dissipation_entropy: Option<SplitMolecularEntropyReport>,
}

impl ReactingEuler1d {
    /// Advance original chemistry around source-conjugate inviscid transport,
    /// admitting the actual entropy of the complete closed composition.
    ///
    /// The original fixed-U,V reactor retains its source bank, kinetic laws,
    /// selected backward-Euler or SDIRK method, temperature interval and all
    /// physical/error/work checks. Chemical heat release enters through the
    /// existing species-energy gauge change exactly once. The middle native
    /// SSP-RK2 flow step uses the original variable-caloric entropy-stable face.
    ///
    /// The composition C(dt/2), A(dt), C(dt/2) uses one actual duration. A shorter
    /// acoustic or entropy-admitted flow step restarts BOTH chemical halves from
    /// accepted state. A refused complete entropy/source account also retries
    /// the whole composition under the existing bounded halving policy.
    ///
    /// First-order reconstruction, fixed constant area, and periodic or two
    /// stationary reflecting ends are required. The actual complete source
    /// entropy includes chemical production. Only the entropy-reference effect
    /// of numerical species-source drift is removed, with each such drift
    /// independently admitted against original mass plus absolute chemical
    /// throughput. Species created by chemistry are not mistaken for mass drift.
    ///
    /// The same caller-selected entropy allowance applies to the middle flow
    /// step and the direct original-to-final composition. The original reactor
    /// retains its own declared entropy tolerance. Neither split symmetry nor
    /// component admission establishes an unconditional entropy/positivity or
    /// global temporal-order theorem. Refusal/cancellation leaves the complete
    /// accepted flow and clock unchanged.
    pub fn step_entropy_consistent(
        &mut self,
        maximum_duration: Time,
        controls: ReactiveStepControl,
        entropy: MixtureEntropyControl,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<ReactiveEntropyStepReport, ReactiveError> {
        let mut account =
            CoupledEntropy::new(&self.accepted, controls.transport, entropy, &mut cancelled)?;
        let (coupled, _) = self.step_inner::<ConstantTransport>(
            maximum_duration,
            controls,
            None,
            stabilized::MolecularMethod::Explicit(DissipationControl::default()),
            [None; 2],
            None,
            Some(&mut account),
            &mut cancelled,
        )?;
        Ok(account.complete(coupled))
    }

    /// Add original source-conjugate SSP-RK2 molecular half steps to the same
    /// closed chemistry/inviscid composition and complete entropy admission.
    ///
    /// The stages are D(dt/2), C(dt/2), A(dt), C(dt/2), D(dt/2). Both molecular
    /// halves retain the original Maxwell-Stefan, Fourier, viscous and reciprocal
    /// Soret/Dufour owners and must complete their exact half duration. Their
    /// heat and species enthalpy enter the existing conservative energy flux,
    /// with no added reaction heat or separate entropy-generated heat source.
    ///
    /// Molecular endpoint/source refusal or acoustic shortening replays all
    /// five stages from accepted state. The final complete endpoint is checked
    /// against the reactor temperature interval and one whole-state entropy
    /// allowance before publication. Component entropy receipts are retained
    /// separately and do not replace the direct complete-state account.
    ///
    /// Every transport face requires each species positive on both traces or
    /// absent on both. One-sided absence refuses without trace seeding, clipping
    /// or fallback. Open reservoirs, moving/variable-area geometry, reconstruction
    /// above first order and external forcing are outside this closed entry.
    pub fn step_entropy_consistent_with_dissipation<L: TransportLaw>(
        &mut self,
        operator: &mut DissipationOperator<L>,
        maximum_duration: Time,
        controls: ReactiveStepControl,
        dissipation_controls: DissipationControl,
        entropy: MixtureEntropyControl,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<ReactiveEntropyStepReport, ReactiveError> {
        let mut account =
            CoupledEntropy::new(&self.accepted, controls.transport, entropy, &mut cancelled)?;
        let (coupled, _) = self.step_inner(
            maximum_duration,
            controls,
            Some(operator),
            stabilized::MolecularMethod::EntropyExplicit {
                controls: dissipation_controls,
                entropy,
            },
            [None; 2],
            None,
            Some(&mut account),
            &mut cancelled,
        )?;
        Ok(account.complete(coupled))
    }
}

/// One immutable initial anchor with retry-local chemical transfer accounting.
/// Component receipts and complete entropy arithmetic finish before the native
/// accepted/trial swap. Legacy compositions do not construct this account.
pub(super) struct CoupledEntropy {
    controls: MixtureEntropyControl,
    account: ReactiveEntropyAccount,
    completed: Option<(
        ReactiveEntropyReport,
        MixtureEntropyReport,
        Option<SplitMolecularEntropyReport>,
    )>,
}

impl CoupledEntropy {
    pub(super) fn new(
        flow: &MixtureEuler1d,
        transport: StepControl,
        controls: MixtureEntropyControl,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, ReactiveError> {
        poll(cancelled)?;
        if transport.reconstruction != Reconstruction::FirstOrder {
            return Err(FlowError::EntropyDomain(
                "entropy-consistent reacting flow requires first-order reconstruction",
            )
            .into());
        }
        if flow.mesh.uniform_cross_section().is_none() {
            return Err(FlowError::EntropyDomain(
                "entropy-consistent reacting flow requires fixed constant-area geometry",
            )
            .into());
        }
        if !matches!(
            flow.boundaries,
            Boundaries::Periodic
                | Boundaries::Ends {
                    left: Boundary::Reflecting,
                    right: Boundary::Reflecting,
                }
        ) {
            return Err(FlowError::EntropyDomain(
                "entropy-consistent reacting flow requires periodic or stationary reflecting ends",
            )
            .into());
        }
        Ok(Self {
            controls,
            account: ReactiveEntropyAccount::new(flow, controls, cancelled)?,
            completed: None,
        })
    }

    pub(super) const fn controls(&self) -> MixtureEntropyControl {
        self.controls
    }

    pub(super) fn begin_attempt(&mut self) {
        self.account.begin_attempt();
        self.completed = None;
    }

    pub(super) fn observe_chemistry_row(
        &mut self,
        before: &[f64],
        after: &[f64],
        volume: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), ReactiveError> {
        self.account
            .observe_chemistry_row(before, after, volume, cancelled)?;
        Ok(())
    }

    #[allow(clippy::too_many_arguments)]
    pub(super) fn check_candidate(
        &mut self,
        initial: &MixtureEuler1d,
        candidate: &MixtureEuler1d,
        temperature_bounds: Option<[f64; 2]>,
        transport: MixtureEntropyReport,
        molecular: [Option<MixtureEntropyReport>; 2],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), ReactiveError> {
        let dissipation = match molecular {
            [None, None] => None,
            [Some(first), Some(second)] => Some(SplitMolecularEntropyReport { first, second }),
            _ => {
                return Err(ReactiveError::Invalid(
                    "incomplete reacting molecular entropy receipts",
                ));
            }
        };
        let entropy =
            self.account
                .check_candidate(initial, candidate, temperature_bounds, cancelled)?;
        self.completed = Some((entropy, transport, dissipation));
        Ok(())
    }

    /// No fallible report arithmetic follows native publication.
    pub(super) fn complete(self, coupled: ReactiveStepReport) -> ReactiveEntropyStepReport {
        let (entropy, transport_entropy, dissipation_entropy) = self
            .completed
            .expect("accepted reacting step retained complete entropy admission");
        ReactiveEntropyStepReport {
            coupled,
            entropy,
            transport_entropy,
            dissipation_entropy,
        }
    }
}
