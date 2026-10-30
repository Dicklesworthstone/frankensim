//! One native gas-transport snapshot drives complete continuing acoustic loads.
//!
//! This is the same transaction as the source-wall DG transport adapter. Only
//! the acoustic owner differs: each node retains its bore, finite branches,
//! terminal-channel momentum, exterior poles and aggregate physical accounts.

/// Wall ducts, passive complete loads and moving sources on one gas network.
pub mod assembly;

use super::{FiniteGasDgBackground, FiniteGasDgError, FiniteGasDgPortStep, admit_initial};
use crate::bernoulli_aperture::wall::material::{
    FiniteGasDgTransportCandidate, FiniteGasDgTransportStep, FiniteGasTransportNode,
    advance_finite_gas_transport_with_paths, finite_gas_transport_candidate_with_paths,
};
pub use crate::bernoulli_aperture::wall::material::{
    FiniteGasDgTransportError, FiniteGasDgTransportInput, FiniteGasDgTransportSolve,
};
use crate::dg_port::DgPortLoad;
use fs_material::thermochemical_gas::energy::parcel::ClosedGasParcel;
use fs_material::thermochemical_gas::energy::parcel::composition::GasCompositionStep;
use fs_material::thermochemical_gas::energy::parcel::transfer::{
    nozzle::network::snapshot::SnapshotNozzleLink,
    species::{diffusion::network::DiffusionNetworkLink, network::transport::GasTransportPaths},
};
use std::collections::BTreeMap;

/// Original transport and complete-load reports, with one enclosing mass/energy account.
pub type FiniteGasDgPortTransportStep =
    FiniteGasDgTransportStep<FiniteGasDgPortStep<GasCompositionStep>>;

/// Immutable preview of every finite parcel and complete acoustic endpoint.
pub type FiniteGasDgPortTransportCandidate =
    FiniteGasDgTransportCandidate<DgPortLoad, FiniteGasDgPortStep<GasCompositionStep>>;

impl FiniteGasTransportNode for DgPortLoad {
    type Step<S> = FiniteGasDgPortStep<S>;
    type Input = FiniteGasDgTransportInput;

    fn transport_input(
        input: &Self::Input,
        _bound: bool,
    ) -> Result<&FiniteGasDgTransportInput, &'static str> {
        Ok(input)
    }

    fn dt(&self) -> f64 {
        self.propagator().dt()
    }
    fn step_count(&self) -> u64 {
        self.state().step_count()
    }
    fn stored_energy_j(&self) -> f64 {
        DgPortLoad::stored_energy_j(self)
    }

    fn admit<E, C: FnMut() -> bool>(
        &self,
        gas: &ClosedGasParcel,
        input: &FiniteGasDgTransportInput,
        cancelled: &mut C,
    ) -> Result<(), FiniteGasDgError<E>> {
        if input.ports.len() != 1
            || input.ports[0].id != self.inlet()
            || !input.ports[0].effort.is_finite()
        {
            return Err(FiniteGasDgError::Invalid(
                "complete load transport requires exactly its original inlet and finite held flow",
            ));
        }
        admit_initial(
            self,
            gas,
            input.acoustic.expect("validated acoustic input"),
            cancelled,
        )
    }

    fn advance<S, E, C: FnMut() -> bool>(
        &mut self,
        gas: &mut ClosedGasParcel,
        input: &FiniteGasDgTransportInput,
        molecular_heat_j: f64,
        evolve: impl FnMut(
            &ClosedGasParcel,
            f64,
            f64,
            &mut C,
        ) -> Result<FiniteGasDgBackground<S>, E>,
        cancelled: &mut C,
    ) -> Result<Self::Step<S>, FiniteGasDgError<E>> {
        if molecular_heat_j != 0.0 {
            return Err(FiniteGasDgError::Invalid(
                "complete passive loads do not own the source-wall molecular heat integral",
            ));
        }
        self.step_with_finite_gas_controlled(
            gas,
            input.ports[0].effort,
            input.acoustic.expect("validated acoustic input"),
            evolve,
            cancelled,
        )
    }

    fn accounts<S>(step: &Self::Step<S>) -> [f64; 4] {
        [
            step.acoustic.ledger.port_work_step_j,
            step.acoustic_dissipation_j,
            step.refresh.parameter_work_j,
            0.0,
        ]
    }

    fn source<S>(step: &Self::Step<S>) -> &S { &step.background }
}

/// Advance native nozzle/diffusion transport and complete acoustic loads together.
///
/// Map keys are the original finite parcel identities. Each bound parcel owns
/// the uniform circular main bore PLUS its main-gas finite terminal passages;
/// independent branch gases remain their original sources. Admission checks
/// this actual volume, complete source, work controls and common acoustic clock
/// before the transport snapshot. A load input contains exactly one `PortInput`
/// for its original inlet, with held volume flow in `effort`; a reservoir has
/// no acoustic controls and no ports.
///
/// Original species/partial-enthalpy transfer runs once. Each load then pays its
/// actual bore/channel material work using that unchanged transported candidate.
/// The enclosing energy account includes every represented branch and channel,
/// original loss and net acoustic port work. Rational exterior work is already
/// in that net port account and is never debited twice. Source and endpoint
/// driving-force gates precede joint publication; failure changes no owner.
///
/// These remain separate well-mixed volumes linked by actual material transport,
/// not one continuous heterogeneous acoustic junction. No mean-pressure force,
/// implicit thermodynamic transport from acoustic flow, jet noise or deposition
/// of acoustic/numerical dissipation into gas heat is inferred. Original
/// fixed-geometry, first-order partition and source-domain limits still apply.
#[allow(clippy::too_many_arguments)]
pub fn advance_finite_gas_port_transport<C: FnMut() -> bool>(
    parcels: &mut BTreeMap<u64, ClosedGasParcel>,
    nodes: &mut BTreeMap<u64, DgPortLoad>,
    inputs: &BTreeMap<u64, FiniteGasDgTransportInput>,
    nozzles: &[SnapshotNozzleLink],
    diffusion: &[DiffusionNetworkLink<'_>],
    duration_s: f64,
    controls: FiniteGasDgTransportSolve,
    cancelled: &mut C,
) -> Result<FiniteGasDgPortTransportStep, FiniteGasDgTransportError> {
    advance_finite_gas_port_transport_with_paths(
        parcels,
        nodes,
        inputs,
        GasTransportPaths {
            nozzles,
            diffusion,
            ..Default::default()
        },
        duration_s,
        controls,
        cancelled,
    )
}

/// The same transaction including original porous and membrane transport paths.
/// Every mechanism shares one unchanged gross-species inventory and every final
/// driving-force gate is rechecked after all heating and acoustic work.
#[allow(clippy::too_many_arguments)]
pub fn advance_finite_gas_port_transport_with_paths<C: FnMut() -> bool>(
    parcels: &mut BTreeMap<u64, ClosedGasParcel>,
    nodes: &mut BTreeMap<u64, DgPortLoad>,
    inputs: &BTreeMap<u64, FiniteGasDgTransportInput>,
    paths: GasTransportPaths<'_>,
    duration_s: f64,
    controls: FiniteGasDgTransportSolve,
    cancelled: &mut C,
) -> Result<FiniteGasDgPortTransportStep, FiniteGasDgTransportError> {
    advance_finite_gas_transport_with_paths(
        parcels, nodes, inputs, paths, duration_s, controls, cancelled,
    )
}

/// Preview the exact same transaction without publishing any parcel or load.
/// The returned candidate exposes both endpoint maps and the original physical
/// reports; it contains all channel, branch, exterior and source histories.
#[allow(clippy::too_many_arguments)]
pub fn finite_gas_port_transport_candidate_with_paths<C: FnMut() -> bool>(
    parcels: &BTreeMap<u64, ClosedGasParcel>,
    nodes: &BTreeMap<u64, DgPortLoad>,
    inputs: &BTreeMap<u64, FiniteGasDgTransportInput>,
    paths: GasTransportPaths<'_>,
    duration_s: f64,
    controls: FiniteGasDgTransportSolve,
    cancelled: &mut C,
) -> Result<FiniteGasDgPortTransportCandidate, FiniteGasDgTransportError> {
    finite_gas_transport_candidate_with_paths(
        parcels, nodes, inputs, paths, duration_s, controls, cancelled,
    )
}
