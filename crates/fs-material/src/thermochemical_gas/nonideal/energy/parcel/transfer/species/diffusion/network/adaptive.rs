//! Error-controlled complete windows on the original real-gas inventory owner.
//!
//! The shared gas controller compares one full step with two conservative half
//! steps. It accepts only the fine endpoints; no Richardson extrapolation,
//! trace clipping or surrogate EOS enters a physical state. Every probe samples
//! its own actual PR fugacities, density, qualified binary data and partial
//! enthalpies. Source receipts remain those of the last accepted half step.

use super::{
    Error as StepError, GasEnergyError, GasEnergySolve, IsochoricRealGasParcel,
    RealGasDiffusionNetworkCandidate, RealGasDiffusionNetworkLink, RealGasDiffusionNetworkSolve,
    RealGasDiffusionNetworkStep, audit, positive, step_counted_candidate,
};
use crate::thermochemical_gas::energy::parcel::transfer::{
    inventory,
    species::{
        diffusion::network::{
            DiffusionNetworkSolve,
            adaptive::{
                Candidate, DiffusionNetworkAdvanceCandidate, DiffusionNetworkAdvanceReport,
                Inventory, Transfer, TrialFailure, advance_with_kernel,
            },
        },
        network::transport::GasTransportWindowAudit,
    },
};
use fs_thermochem::SpeciesId;
use std::collections::BTreeMap;

pub use crate::thermochemical_gas::energy::parcel::transfer::species::diffusion::network::adaptive::{
    DiffusionNetworkAdvanceControls,
    DiffusionNetworkAdvanceError as RealGasDiffusionNetworkAdvanceError,
    DiffusionNetworkIntegratedTransfer,
};
use RealGasDiffusionNetworkAdvanceError as Error;

/// The same endpoint accuracy criteria as the ideal/thermal gas owners, plus a
/// cumulative budget for actual PR/Maxwell--Stefan face calls. One active face
/// call counts once even for a multicomponent catalog; its original <=32-species
/// and complete-pair bounds remain. Closed faces consume no face calls.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RealGasDiffusionNetworkAdvanceControls {
    pub temporal: DiffusionNetworkAdvanceControls,
    /// All coarse, fine, rejected and failed calls; 1..=1,000,000.
    pub maximum_face_evaluations: usize,
}

/// Shared accepted species/enthalpy ledger with an actual PR final-state audit.
/// `pair_evaluations` is zero here; `face_evaluations` counts original real-gas
/// constitutive calls. No intermediate rejected state is retained.
pub type RealGasDiffusionNetworkAdvanceReport =
    DiffusionNetworkAdvanceReport<RealGasDiffusionNetworkStep>;
pub type RealGasDiffusionNetworkAdvanceCandidate =
    DiffusionNetworkAdvanceCandidate<RealGasDiffusionNetworkStep, IsochoricRealGasParcel>;

impl Inventory for IsochoricRealGasParcel {
    fn species_inventory(&self) -> Result<BTreeMap<SpeciesId, f64>, GasEnergyError> {
        inventory(self.model().source(), self.mass_kg())
    }
    fn temperature_k(&self) -> f64 {
        self.state().temperature_k()
    }
    fn audit_window(
        before: &BTreeMap<u64, Self>,
        after: &BTreeMap<u64, Self>,
        physical: DiffusionNetworkSolve,
        _elapsed_s: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<GasTransportWindowAudit, Error> {
        audit(
            before,
            after,
            RealGasDiffusionNetworkSolve {
                exchange: physical.exchange,
                entropy_tolerance_j_k: physical.entropy_tolerance_j_k,
            },
            cancelled,
        )
        .map_err(|error| {
            if error == StepError::Cancelled {
                Error::Cancelled
            } else {
                Error::RealGasWindow(error)
            }
        })
    }
}

impl Candidate for RealGasDiffusionNetworkCandidate {
    type Parcel = IsochoricRealGasParcel;
    type Receipt = RealGasDiffusionNetworkStep;
    fn parcels(&self) -> &BTreeMap<u64, IsochoricRealGasParcel> {
        RealGasDiffusionNetworkCandidate::parcels(self)
    }
    fn entropy_change(&self) -> f64 {
        self.step().entropy_change_j_k
    }
    fn energy_residual(&self) -> f64 {
        self.step().exchange.total_energy_residual_j
    }
    fn transfers(&self) -> impl Iterator<Item = Transfer<'_>> {
        self.step().exchange.fluxes.iter().map(|row| Transfer {
            donor: row.transfer.donor,
            receiver: row.transfer.receiver,
            species: &row.transfer.species,
            mass_kg: row.transfer.mass_kg,
            enthalpy_j: row.carried_enthalpy_j,
        })
    }
    fn into_parts(self) -> (BTreeMap<u64, IsochoricRealGasParcel>, Self::Receipt) {
        RealGasDiffusionNetworkCandidate::into_parts(self)
    }
}

/// Advance exactly the requested positive duration, publishing every endpoint
/// and accepted gross transfer together or returning no candidate.
///
/// Local error is the full-minus-two-half difference in EACH actual species
/// inventory and original PR temperature (first-order divisor 2^1-1 = 1).
/// Only native gross-withdrawal, endpoint caloric-bracket and actual entropy
/// refusals may trigger a shorter step; source/domain, phase/mechanical,
/// conservation and nonconvergence failures propagate without reinterpretation.
/// Rejected probes spend cumulative work but no accepted mass or enthalpy.
///
/// The accepted enthalpy residual ledger and independent original-to-final
/// species/U/S audit each use the same fixed physical budgets, never budgets
/// multiplied by the number of substeps. PR departure energy and density-chart
/// entropy remain in their original audit; formation references cancel before
/// arithmetic. Arrivals can diffuse onward only during later physical substeps.
///
/// This estimates local temporal error, not global/spatial accuracy, implicit
/// stability or phase equilibrium. The original isothermal-face approximation,
/// actual source identities, positive common species support and T/p-only
/// binary-coefficient qualification remain unchanged throughout the window.
#[allow(clippy::too_many_arguments)]
pub fn advance_candidate(
    original: &BTreeMap<u64, IsochoricRealGasParcel>,
    links: &[RealGasDiffusionNetworkLink<'_>],
    duration_s: f64,
    endpoints: &BTreeMap<u64, GasEnergySolve>,
    physical: RealGasDiffusionNetworkSolve,
    controls: RealGasDiffusionNetworkAdvanceControls,
    mut cancelled: impl FnMut() -> bool,
) -> Result<RealGasDiffusionNetworkAdvanceCandidate, Error> {
    if cancelled() {
        return Err(Error::Cancelled);
    }
    if ![
        physical.exchange.maximum_species_fraction,
        physical.exchange.mass_tolerance_kg,
        physical.exchange.energy_tolerance_j,
        physical.entropy_tolerance_j_k,
    ]
    .into_iter()
    .all(positive)
        || physical.exchange.maximum_species_fraction >= 1.0
        || !(1..=1_000_000).contains(&controls.maximum_face_evaluations)
    {
        return Err(Error::Invalid(
            "positive original real-gas physical and bounded face-work controls required",
        ));
    }
    let mut workspaces = BTreeMap::new();
    advance_with_kernel(
        original,
        links.len(),
        duration_s,
        endpoints,
        DiffusionNetworkSolve {
            exchange: physical.exchange,
            entropy_tolerance_j_k: physical.entropy_tolerance_j_k,
        },
        controls.temporal,
        None,
        Some(controls.maximum_face_evaluations),
        |state, dt, work, cancelled| {
            step_counted_candidate(
                state,
                links,
                dt,
                endpoints,
                physical,
                &mut workspaces,
                &mut work.face_evaluations,
                controls.maximum_face_evaluations,
                cancelled,
            )
            .map_err(|error| match error {
                StepError::FaceWorkLimit { .. } => TrialFailure::FaceWork,
                other => TrialFailure::RealGas(other),
            })
        },
        cancelled,
    )
}
