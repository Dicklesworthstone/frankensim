//! Original-source entropy with native prescribed-reservoir network transport.

use super::*;
use crate::multicomponent::dissipation::entropy::EntropyTracker;
pub use crate::multicomponent::dissipation::entropy::reactive::ReactiveExteriorEntropyReport;
use crate::multicomponent::{ExteriorMixtureEntropyReport, StepControl};

/// Complete reacting network state and its actual prescribed entropy supply.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ReactiveNetworkExteriorEntropyReport {
    /// Native chemistry, transport, actual duration and conservative balances.
    pub coupled: ReactiveNetworkReport,
    /// Direct original-to-final source entropy minus the actual exterior supply.
    /// Real chemical conversion remains in the admitted entropy generation.
    pub entropy: ReactiveExteriorEntropyReport,
    /// Middle native network update in its own initial-state entropy basis.
    /// This component report does not replace complete-step admission.
    pub transport_entropy: ExteriorMixtureEntropyReport,
}

impl ReactiveNetworkExteriorEntropyReport {
    fn from_closed(closed: ReactiveNetworkEntropyReport) -> Self {
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
                reference_entropy_residual_j_k: closed
                    .transport_entropy
                    .reference_entropy_residual_j_k,
                minimum_step_entropy_generation_j_k: closed
                    .transport_entropy
                    .minimum_step_entropy_change_j_k,
                maximum_relative_species_residual: closed
                    .transport_entropy
                    .maximum_relative_species_residual,
                rejected_entropy_steps: closed.transport_entropy.rejected_entropy_steps,
            },
        }
    }
}

impl ReactingNetwork1d {
    /// Advance C(dt/2), A(dt), C(dt/2) with prescribed-reservoir entropy supply.
    ///
    /// Every branch cell and finite plenum retains the original fixed-U,V
    /// reactor and caloric source. The middle transport uses the existing
    /// characteristic pressure-reservoir closure and measured native network
    /// entropy gate. One optional stationary reservoir is supplied per physical
    /// branch end; a finite-plenum port cannot be replaced by a reservoir.
    /// Unowned exterior ends must be stationary reflecting walls. All absent
    /// reservoirs execute the closed entropy composition with identical physics.
    ///
    /// The complete account is anchored before chemistry. Each actual accepted
    /// SSP species, momentum and referenced-energy transfer is contracted with
    /// that end's prescribed thermodynamic efforts in this SAME source basis.
    /// Advected enthalpy already includes pressure-flow work. Receiving-reservoir
    /// efforts also apply to outflow, so total generation includes the reservoir
    /// coupling's irreversibility. No heat, chemical energy or p*dV is added twice.
    /// Internal branch/plenum transfers cancel inside the state inventory.
    ///
    /// Real chemical species entropy remains physical. Only the source-reference
    /// effect of numerical species-source drift is removed, after bounding each
    /// drift against original inventory plus absolute chemical and both-stage
    /// boundary throughput. A zero reservoir activity refuses any nonzero species
    /// transfer at either SSP stage; opposite transfers cannot conceal it.
    ///
    /// Acoustic shortening or entropy refusal restarts the entire composition.
    /// Failure or final cancellation preserves all accepted fields, clocks and
    /// receipts. This fixed-network entry excludes molecular baths, wall/valve
    /// sources and moving volumes. Measured admission of the engineering flux
    /// and symmetric splitting establish no fully discrete stability/order theorem.
    pub fn step_with_reservoirs_and_entropy(
        &mut self,
        maximum_duration: Time,
        controls: ReactiveStepControl,
        reservoirs: &[[Option<Cell<'_>>; 2]],
        entropy: MixtureEntropyControl,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<ReactiveNetworkExteriorEntropyReport, ReactiveNetworkError> {
        let active = self.accepted.admit_reservoirs(reservoirs, &mut cancelled)?;
        if !active {
            return self
                .step_with_entropy(maximum_duration, controls, entropy, cancelled)
                .map(ReactiveNetworkExteriorEntropyReport::from_closed);
        }
        let mut account = CoupledEntropy::new_with_reservoirs(
            &self.accepted,
            entropy,
            reservoirs,
            &mut cancelled,
        )?;
        let (coupled, _) = self.step_method::<ConstantTransport>(
            maximum_duration,
            controls,
            None,
            MolecularMethod::Explicit(DissipationControl::default()),
            NetworkSources::default(),
            None,
            Some(reservoirs),
            Some(&mut account),
            &mut cancelled,
        )?;
        Ok(account.complete_exterior(coupled))
    }
}

/// Only an open call retains its middle account, including opaque physical SSP
/// inventories. Closed calls allocate nothing for this optional state.
pub(super) struct ExteriorCoupling {
    transport: Option<EntropyTracker>,
    completed: Option<ReactiveExteriorEntropyReport>,
}

impl ExteriorCoupling {
    pub(super) fn reset(&mut self) {
        self.transport = None;
        self.completed = None;
    }
}

impl CoupledEntropy {
    fn new_with_reservoirs(
        network: &MixtureNetwork1d,
        controls: MixtureEntropyControl,
        reservoirs: &[[Option<Cell<'_>>; 2]],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, ReactiveNetworkError> {
        let account = ReactiveEntropyAccount::from_tracker(network.entropy_account(
            controls,
            Some(reservoirs),
            cancelled,
        )?);
        let mut exterior = Vec::new();
        exterior
            .try_reserve_exact(1)
            .map_err(|_| FlowError::Allocation)?;
        exterior.push(ExteriorCoupling {
            transport: None,
            completed: None,
        });
        Ok(Self {
            controls,
            account,
            transport: None,
            completed: None,
            exterior,
        })
    }

    pub(super) fn advance_exterior_transport(
        &mut self,
        network: &mut MixtureNetwork1d,
        duration: Time,
        controls: StepControl,
        reservoirs: Option<&[[Option<Cell<'_>>; 2]]>,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<NetworkStepReport, ReactiveNetworkError> {
        let reservoirs = reservoirs.ok_or(ReactiveNetworkError::Invalid(
            "open reacting entropy requires its prescribed reservoir states",
        ))?;
        let (report, account) = network.step_with_reservoirs_and_entropy_account(
            duration,
            controls,
            reservoirs,
            self.controls,
            cancelled,
        )?;
        self.exterior[0].transport = Some(account);
        Ok(report)
    }

    pub(super) fn finish_exterior<'a, 'b>(
        &mut self,
        gas: &MixtureGas,
        cells: impl IntoIterator<Item = (Cell<'a>, Cell<'b>, f64)>,
        temperature_bounds: Option<[f64; 2]>,
        molecular: Option<SplitMolecularEntropyReport>,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), ReactiveNetworkError> {
        if molecular.is_some() {
            return Err(ReactiveNetworkError::Invalid(
                "open reacting entropy requires actual molecular bath transfer inventories",
            ));
        }
        let exterior = &mut self.exterior[0];
        let transport = exterior
            .transport
            .as_ref()
            .ok_or(ReactiveNetworkError::Invalid(
                "open reacting entropy requires an admitted middle network transport",
            ))?;
        let report = self.account.check_candidate_with_exterior_cells(
            gas,
            cells,
            temperature_bounds,
            transport,
            cancelled,
        )?;
        exterior.completed = Some(report);
        Ok(())
    }

    fn complete_exterior(
        self,
        coupled: ReactiveNetworkReport,
    ) -> ReactiveNetworkExteriorEntropyReport {
        let exterior = self
            .exterior
            .into_iter()
            .next()
            .expect("accepted reacting network retained its exterior account");
        ReactiveNetworkExteriorEntropyReport {
            coupled,
            entropy: exterior
                .completed
                .expect("complete reacting network entropy admitted"),
            transport_entropy: exterior
                .transport
                .expect("middle reservoir transport admitted")
                .exterior_report(),
        }
    }
}
