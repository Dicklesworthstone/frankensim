//! One original chemical/source entropy account for the native network transaction.

use super::*;
use crate::multicomponent::dissipation::entropy::reactive::ReactiveEntropyAccount;
pub use crate::multicomponent::reactive::{ReactiveEntropyReport, SplitMolecularEntropyReport};
use crate::multicomponent::{MixtureEntropyControl, MixtureEntropyReport, cell};

/// Complete accepted network evolution with component and whole-state entropy.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ReactiveNetworkEntropyReport {
    /// Original chemistry/transport operators, actual duration and energy balances.
    pub coupled: ReactiveNetworkReport,
    /// Direct original-to-final entropy, retaining actual chemical source entropy.
    /// Component entropy reports do not replace this complete-composition gate.
    pub entropy: ReactiveEntropyReport,
    /// Native measured inviscid update between the chemical half steps.
    pub transport_entropy: MixtureEntropyReport,
    /// Both complete source-conjugate molecular halves, when requested.
    pub dissipation_entropy: Option<SplitMolecularEntropyReport>,
}

impl ReactingNetwork1d {
    /// Advance closed chemistry/advection with one complete source-entropy gate.
    ///
    /// The existing fixed-U,V reactor advances every branch cell and finite
    /// plenum around the original entropy-admitted network advection step.
    /// Actual represented per-cell chemical species changes enter the entropy
    /// account, including the original species entropy constants. Species
    /// created or consumed by chemistry are physical sources; only the reference
    /// effect of numerical species-source continuity defects is removed from
    /// admission. Each defect is bounded against original inventory plus the
    /// absolute chemical throughput of both halves and every storage owner.
    ///
    /// Acoustic shortening or entropy refusal restarts both chemical halves
    /// from accepted state under the original bounded retry policy. The final
    /// complete original-to-final source state receives one entropy allowance,
    /// with no summed component tolerances or added reaction heat. All fields,
    /// clocks and accepted receipts publish together after final admission.
    ///
    /// External ends must be stationary reflecting walls; isolated periodic
    /// branches are admitted. This closed fixed-volume entry excludes baths,
    /// wall/valve sources and moving plenums. The native engineering advection
    /// flux is measured, not replaced: neither entropy admission nor symmetric
    /// splitting establishes a fully discrete stability or temporal-order theorem.
    pub fn step_with_entropy(
        &mut self,
        maximum_duration: Time,
        controls: ReactiveStepControl,
        entropy: MixtureEntropyControl,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<ReactiveNetworkEntropyReport, ReactiveNetworkError> {
        let mut account = CoupledEntropy::new(&self.accepted, entropy, &mut cancelled)?;
        let (coupled, _) = self.step_method::<ConstantTransport>(
            maximum_duration,
            controls,
            None,
            MolecularMethod::Explicit(DissipationControl::default()),
            NetworkSources::default(),
            None,
            None,
            Some(&mut account),
            &mut cancelled,
        )?;
        Ok(account.complete(coupled))
    }

    /// Compose D/2, C/2, A, C/2, D/2 in the same complete entropy transaction.
    ///
    /// Both molecular halves use the original source-conjugate network SSPRK2
    /// operator, including finite plenum storage, authored port mixing lengths,
    /// Maxwell--Stefan diffusion, divided source enthalpies, Fourier heat and
    /// viscous work. They complete the actual half duration and their own
    /// substep/whole-interval entropy gates before composition can succeed.
    /// Internal molecular ports are not an external entropy supply.
    ///
    /// Every molecular face requires positive species on both traces or absence
    /// on both. No concentration floor is used. The final molecular state must
    /// remain in the original reactor's temperature interval. A refusal replays
    /// all five stages without publishing intermediate molecular port reports.
    /// The complete gate retains the real chemical entropy; molecular heat or
    /// reaction energy is never added again as an independent source.
    ///
    /// Closed stationary constant-area molecular geometry is required. Baths,
    /// moving volumes, other sources and stabilized/adaptive molecular methods
    /// remain outside this entry. Component gates are measured admissions, not
    /// a proof of fully discrete coupled entropy stability.
    pub fn step_with_dissipation_and_entropy<L: TransportLaw>(
        &mut self,
        operator: &mut NetworkDissipationOperator<L>,
        maximum_duration: Time,
        controls: ReactiveStepControl,
        dissipation_controls: DissipationControl,
        entropy: MixtureEntropyControl,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<ReactiveNetworkEntropyReport, ReactiveNetworkError> {
        let mut account = CoupledEntropy::new(&self.accepted, entropy, &mut cancelled)?;
        let (coupled, _) = self.step_method(
            maximum_duration,
            controls,
            Some(operator),
            MolecularMethod::EntropyExplicit {
                controls: dissipation_controls,
                entropy,
            },
            NetworkSources::default(),
            None,
            None,
            Some(&mut account),
            &mut cancelled,
        )?;
        Ok(account.complete(coupled))
    }
}

/// Retry-local source transfers in one immutable original inventory basis.
pub(super) struct CoupledEntropy {
    controls: MixtureEntropyControl,
    account: ReactiveEntropyAccount,
    transport: Option<MixtureEntropyReport>,
    completed: Option<(ReactiveEntropyReport, Option<SplitMolecularEntropyReport>)>,
}

impl CoupledEntropy {
    fn new(
        network: &MixtureNetwork1d,
        controls: MixtureEntropyControl,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, ReactiveNetworkError> {
        Ok(Self {
            controls,
            account: ReactiveEntropyAccount::from_tracker(
                network.entropy_account(controls, None, cancelled)?,
            ),
            transport: None,
            completed: None,
        })
    }

    pub(super) fn begin_attempt(&mut self) {
        self.account.begin_attempt();
        self.transport = None;
        self.completed = None;
    }

    pub(super) fn observe_chemistry_row(
        &mut self,
        before: &[f64],
        after: &[f64],
        volume: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), ReactiveNetworkError> {
        self.account
            .observe_chemistry_row(before, after, volume, cancelled)?;
        Ok(())
    }

    pub(super) fn advance_transport(
        &mut self,
        network: &mut MixtureNetwork1d,
        duration: Time,
        controls: crate::multicomponent::StepControl,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<NetworkStepReport, ReactiveNetworkError> {
        let (transport, entropy) =
            network.step_with_entropy(duration, controls, self.controls, cancelled)?;
        self.transport = Some(entropy);
        Ok(transport)
    }

    pub(super) fn finish(
        &mut self,
        initial: &MixtureNetwork1d,
        candidate: &MixtureNetwork1d,
        chemistry: IsochoricKineticsSolve,
        molecular: Option<SplitMolecularEntropyReport>,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), ReactiveNetworkError> {
        let n = initial.gas().species_count();
        let m = n + 2;
        let original = initial
            .branches
            .iter()
            .flat_map(|branch| {
                branch
                    .states
                    .chunks_exact(m)
                    .zip(branch.mesh.volumes())
                    .map(move |(q, &v)| (cell(q, n), v))
            })
            .chain(
                initial
                    .plenums
                    .iter()
                    .map(|plenum| (plenum.state(), plenum.volume_m3())),
            );
        let final_cells = candidate
            .branches
            .iter()
            .flat_map(|branch| branch.states.chunks_exact(m).map(move |q| cell(q, n)))
            .chain(candidate.plenums.iter().map(Plenum::state));
        let cells = original
            .zip(final_cells)
            .map(|((old, volume), new)| (old, new, volume));
        let report = self.account.check_candidate_cells(
            candidate.gas(),
            cells,
            Some([
                chemistry.minimum_temperature_k,
                chemistry.maximum_temperature_k,
            ]),
            cancelled,
        )?;
        self.completed = Some((report, molecular));
        Ok(())
    }

    fn complete(self, coupled: ReactiveNetworkReport) -> ReactiveNetworkEntropyReport {
        let (entropy, dissipation_entropy) =
            self.completed.expect("admitted complete reacting entropy");
        ReactiveNetworkEntropyReport {
            coupled,
            entropy,
            transport_entropy: self.transport.expect("admitted middle network entropy"),
            dissipation_entropy,
        }
    }
}
