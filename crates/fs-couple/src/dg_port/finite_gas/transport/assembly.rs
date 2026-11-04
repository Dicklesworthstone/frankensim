//! Heterogeneous original acoustic owners on one conservative gas transaction.
//!
//! This module dispatches only ownership and physical input types. Transport,
//! caloric feedback, final source gates, compensated accounts and joint
//! publication remain in the same kernel as every homogeneous entry point.

/// Native finite-rate chemistry coupled to transport and acoustic work.
pub mod reactive;

pub use super::{FiniteGasDgTransportError, FiniteGasDgTransportInput, FiniteGasDgTransportSolve};
use crate::bernoulli_aperture::dg::LoadedApertureDg;
use crate::bernoulli_aperture::dg::finite_gas::FiniteGasApertureDgStep;
use crate::bernoulli_aperture::dg::finite_gas::transport::FiniteGasApertureTransportInput;
use crate::bernoulli_aperture::wall::material::{
    FiniteGasDgNode, FiniteGasDgStep, FiniteGasDgTransportCandidate, FiniteGasDgTransportStep,
    FiniteGasDgHeatBudget, FiniteGasDgInternalHeat, FiniteGasTransportAccounts,
    FiniteGasTransportNode, advance_finite_gas_transport_with_paths,
    finite_gas_transport_candidate_with_paths,
};
use crate::dg_port::DgPortLoad;
use crate::dg_port::finite_gas::{
    FiniteGasDgBackground, FiniteGasDgError, FiniteGasDgPortStep,
    PreparedFiniteGasDgPortHeat,
};
use fs_material::thermochemical_gas::energy::parcel::ClosedGasParcel;
use fs_material::thermochemical_gas::energy::parcel::composition::GasCompositionStep;
use fs_material::thermochemical_gas::energy::parcel::transfer::species::network::transport::GasTransportPaths;
use fs_phs::linear_exact::DenseExponentialBudget;
use std::collections::BTreeMap;

/// An original accepted physical owner, keyed by its finite gas parcel identity.
/// Cloning preserves every evolving history and never reconstructs a source.
#[derive(Clone)]
pub enum FiniteGasAcousticNode {
    /// Original uniform or regional source-wall binding, DG state and propagator.
    Wall(FiniteGasDgNode),
    /// Original complete passive bore, branches, finite channels and terminals.
    Port(DgPortLoad),
    /// Original moving mechanics and its complete downstream acoustic load.
    Aperture(LoadedApertureDg),
}

impl From<FiniteGasDgNode> for FiniteGasAcousticNode {
    fn from(value: FiniteGasDgNode) -> Self {
        Self::Wall(value)
    }
}
impl From<DgPortLoad> for FiniteGasAcousticNode {
    fn from(value: DgPortLoad) -> Self {
        Self::Port(value)
    }
}
impl From<LoadedApertureDg> for FiniteGasAcousticNode {
    fn from(value: LoadedApertureDg) -> Self {
        Self::Aperture(value)
    }
}
impl FiniteGasAcousticNode {
    /// Original fixed physical sample duration [s].
    pub fn time_step_s(&self) -> f64 {
        FiniteGasTransportNode::dt(self)
    }
    /// Accepted physical samples of the original owner.
    pub fn accepted_steps(&self) -> u64 {
        FiniteGasTransportNode::step_count(self)
    }
    /// All represented mechanical, material and acoustic storage [J].
    pub fn stored_energy_j(&self) -> f64 {
        FiniteGasTransportNode::stored_energy_j(self)
    }
}

/// Native physical controls; input kind must match its actual bound owner.
#[derive(Clone, Debug)]
pub enum FiniteGasAcousticInput {
    /// A wall's held native ports, a complete load's one original inlet, or an
    /// unbound reservoir's caloric controls with no acoustic controls or ports.
    Passive(FiniteGasDgTransportInput),
    /// A source-wall DG volume returns its native molecular relaxation heat to
    /// its own finite gas. The original slip-square integral is prepared once
    /// from accepted history and held ports, then reused by every coupled trial.
    /// Other native losses stay exported. Requires the original molecular law,
    /// unchanged species fractions and reservoir support; complete port loads,
    /// moving apertures and unbound reservoirs cannot use this input variant.
    MolecularWall {
        thermal: FiniteGasDgTransportInput,
        heat_budget: DenseExponentialBudget,
    },
    /// Return both original molecular and fluid-viscosity heat from one accepted
    /// source-wall DG sample. Each channel has its own complete preparation
    /// allowance; both integrals are reused across every caloric, chemistry and
    /// thermal/work trial. Thermal-wall, solid and terminal losses stay exported.
    /// Source/clock/volume restrictions are the same as `MolecularWall`; a wall
    /// without molecular or viscous memory returns zero for that absent source.
    ViscousMolecularWall {
        thermal: FiniteGasDgTransportInput,
        heat_budgets: FiniteGasDgHeatBudget,
    },
    /// Return molecular and viscous heat belonging to a complete load's main
    /// finite gas inventory, including its resolved terminal passages. One
    /// native acoustic trial determines the actual implicit terminal and branch
    /// inputs before transport; that SAME trial and its heat are reused by all
    /// caloric, chemistry and thermal/work endpoint solves. Independent branch
    /// gas, exterior, solid, thermal-wall and numerical losses stay exported.
    ViscousMolecularPort {
        thermal: FiniteGasDgTransportInput,
        heat_budgets: FiniteGasDgHeatBudget,
    },
    /// Held drives of an actual moving source, with no prescribed acoustic flow.
    Aperture(FiniteGasApertureTransportInput),
}

/// Native per-owner reports without erasing their distinct physical accounts.
#[derive(Debug)]
pub enum FiniteGasAcousticStep<S = GasCompositionStep> {
    /// Source-wall DG and finite gas endpoint.
    Wall(FiniteGasDgStep<S>),
    /// Complete passive load and finite gas endpoint.
    Port(FiniteGasDgPortStep<S>),
    /// Original moving mechanics, complete load and finite gas endpoint.
    Aperture(FiniteGasApertureDgStep<S>),
}

/// One original transport snapshot and one enclosing mass/energy account.
pub type FiniteGasAcousticTransportStep = FiniteGasDgTransportStep<FiniteGasAcousticStep>;
/// Unpublished endpoints of every original gas and heterogeneous physical owner.
pub type FiniteGasAcousticTransportCandidate =
    FiniteGasDgTransportCandidate<FiniteGasAcousticNode, FiniteGasAcousticStep>;

fn mismatch<E>() -> FiniteGasDgError<E> {
    FiniteGasDgError::Invalid("gas transport input kind must match its original acoustic owner")
}

/// Heat and its original solved complete-load trajectory remain one preparation.
/// It is private to this transaction and never replaces accepted node history.
#[derive(Default)]
pub(crate) struct PreparedFiniteGasAcousticHeat {
    heat: FiniteGasDgInternalHeat,
    port: Option<PreparedFiniteGasDgPortHeat>,
}

impl FiniteGasTransportNode for FiniteGasAcousticNode {
    type Input = FiniteGasAcousticInput;
    type Step<S> = FiniteGasAcousticStep<S>;
    type PreparedHeat = PreparedFiniteGasAcousticHeat;

    fn internal_heat(prepared: &Self::PreparedHeat) -> FiniteGasDgInternalHeat {
        prepared.heat
    }

    fn transport_input(
        input: &Self::Input,
        bound: bool,
    ) -> Result<&FiniteGasDgTransportInput, &'static str> {
        match input {
            FiniteGasAcousticInput::Passive(input) => Ok(input),
            FiniteGasAcousticInput::MolecularWall { thermal, .. }
            | FiniteGasAcousticInput::ViscousMolecularWall { thermal, .. }
            | FiniteGasAcousticInput::ViscousMolecularPort { thermal, .. } if bound => Ok(thermal),
            FiniteGasAcousticInput::Aperture(input) if bound => {
                LoadedApertureDg::transport_input(input, bound)
            }
            FiniteGasAcousticInput::Aperture(_)
            | FiniteGasAcousticInput::MolecularWall { .. }
            | FiniteGasAcousticInput::ViscousMolecularWall { .. }
            | FiniteGasAcousticInput::ViscousMolecularPort { .. } => {
                Err("unbound gas reservoirs require plain caloric controls")
            }
        }
    }

    fn dt(&self) -> f64 {
        match self {
            Self::Wall(node) => FiniteGasTransportNode::dt(node),
            Self::Port(node) => FiniteGasTransportNode::dt(node),
            Self::Aperture(node) => FiniteGasTransportNode::dt(node),
        }
    }

    fn step_count(&self) -> u64 {
        match self {
            Self::Wall(node) => FiniteGasTransportNode::step_count(node),
            Self::Port(node) => FiniteGasTransportNode::step_count(node),
            Self::Aperture(node) => FiniteGasTransportNode::step_count(node),
        }
    }

    fn stored_energy_j(&self) -> f64 {
        match self {
            Self::Wall(node) => FiniteGasTransportNode::stored_energy_j(node),
            Self::Port(node) => FiniteGasTransportNode::stored_energy_j(node),
            Self::Aperture(node) => FiniteGasTransportNode::stored_energy_j(node),
        }
    }

    fn admit<E, C: FnMut() -> bool>(
        &self,
        gas: &ClosedGasParcel,
        input: &Self::Input,
        cancelled: &mut C,
    ) -> Result<(), FiniteGasDgError<E>> {
        // The shared kernel runs this exact-kind gate for EVERY node before
        // evaluating any material path. No drive can be silently ignored.
        match (self, input) {
            (Self::Wall(node), FiniteGasAcousticInput::Passive(input))
            | (Self::Wall(node), FiniteGasAcousticInput::MolecularWall { thermal: input, .. })
            | (Self::Wall(node), FiniteGasAcousticInput::ViscousMolecularWall { thermal: input, .. }) => {
                node.admit(gas, input, cancelled)
            }
            (Self::Port(node), FiniteGasAcousticInput::Passive(input))
            | (Self::Port(node), FiniteGasAcousticInput::ViscousMolecularPort { thermal: input, .. }) => {
                node.admit(gas, input, cancelled)
            }
            (Self::Aperture(node), FiniteGasAcousticInput::Aperture(input)) => {
                node.admit(gas, input, cancelled)
            }
            _ => Err(mismatch()),
        }
    }

    fn prepare_internal_heat<E, C: FnMut() -> bool>(
        &self,
        gas: &ClosedGasParcel,
        input: &Self::Input,
        cancelled: &mut C,
    ) -> Result<Self::PreparedHeat, FiniteGasDgError<E>> {
        match (self, input) {
            (Self::Wall(node), FiniteGasAcousticInput::MolecularWall { thermal, heat_budget }) => {
                Ok(PreparedFiniteGasAcousticHeat {
                    heat: FiniteGasDgInternalHeat {
                        molecular_heat_j: node.prepare_molecular_heat(thermal, *heat_budget, cancelled)?,
                        viscous_heat_j: 0.0,
                    },
                    port: None,
                })
            }
            (Self::Wall(node), FiniteGasAcousticInput::ViscousMolecularWall { thermal, heat_budgets }) => {
                Ok(PreparedFiniteGasAcousticHeat {
                    heat: node.prepare_viscous_molecular_heat(thermal, *heat_budgets, cancelled)?,
                    port: None,
                })
            }
            (Self::Port(node), FiniteGasAcousticInput::ViscousMolecularPort { thermal, heat_budgets }) => {
                let prepared = node.prepare_finite_gas_source_heat(
                    gas, thermal.ports[0].effort,
                    thermal.acoustic.expect("admitted acoustic controls"),
                    *heat_budgets, cancelled,
                )?;
                Ok(PreparedFiniteGasAcousticHeat {
                    heat: prepared.heat(),
                    port: Some(prepared),
                })
            }
            (_, FiniteGasAcousticInput::MolecularWall { .. }
                | FiniteGasAcousticInput::ViscousMolecularWall { .. }
                | FiniteGasAcousticInput::ViscousMolecularPort { .. }) => Err(mismatch()),
            _ => Ok(PreparedFiniteGasAcousticHeat::default()),
        }
    }

    fn advance<S, E, C: FnMut() -> bool>(
        &mut self,
        gas: &mut ClosedGasParcel,
        input: &Self::Input,
        prepared: &Self::PreparedHeat,
        evolve: impl FnMut(
            &ClosedGasParcel,
            f64,
            f64,
            &mut C,
        ) -> Result<FiniteGasDgBackground<S>, E>,
        cancelled: &mut C,
    ) -> Result<Self::Step<S>, FiniteGasDgError<E>> {
        match (self, input) {
            (Self::Wall(node), FiniteGasAcousticInput::Passive(input))
            | (Self::Wall(node), FiniteGasAcousticInput::MolecularWall { thermal: input, .. })
            | (Self::Wall(node), FiniteGasAcousticInput::ViscousMolecularWall { thermal: input, .. }) => node
                .advance(gas, input, &prepared.heat, evolve, cancelled)
                .map(FiniteGasAcousticStep::Wall),
            (Self::Port(node), FiniteGasAcousticInput::Passive(input)) => node
                .advance(gas, input, &prepared.heat, evolve, cancelled)
                .map(FiniteGasAcousticStep::Port),
            (Self::Port(node), FiniteGasAcousticInput::ViscousMolecularPort { thermal, .. }) => {
                let original = prepared.port.as_ref().ok_or_else(mismatch)?;
                let candidate = node.complete_finite_gas_source_heat(
                    original, gas, thermal.acoustic.expect("admitted acoustic controls"),
                    evolve, cancelled,
                )?;
                *node = candidate.load;
                *gas = candidate.gas;
                Ok(FiniteGasAcousticStep::Port(candidate.step))
            }
            (Self::Aperture(node), FiniteGasAcousticInput::Aperture(input)) => node
                .advance(gas, input, &prepared.heat, evolve, cancelled)
                .map(FiniteGasAcousticStep::Aperture),
            _ => Err(mismatch()),
        }
    }

    fn accounts<S>(step: &Self::Step<S>) -> FiniteGasTransportAccounts {
        match step {
            FiniteGasAcousticStep::Wall(step) => FiniteGasDgNode::accounts(step),
            FiniteGasAcousticStep::Port(step) => DgPortLoad::accounts(step),
            FiniteGasAcousticStep::Aperture(step) => LoadedApertureDg::accounts(step),
        }
    }

    fn source<S>(step: &Self::Step<S>) -> &S {
        match step {
            FiniteGasAcousticStep::Wall(step) => &step.background,
            FiniteGasAcousticStep::Port(step) => &step.background,
            FiniteGasAcousticStep::Aperture(step) => &step.gas.background,
        }
    }
}

/// Advance wall ducts, passive loads, moving sources and finite reservoirs on
/// ONE native material-transport snapshot and physical clock.
///
/// Each owner keeps its original source, volume, loss and physical input gates.
/// Its own solver advances its complete state; this adapter supplies no new
/// transport, caloric or mechanical law. Native external work accounts compose
/// exactly once: valve-to-duct work stays internal, and complete-load external
/// port work already includes signed exterior exchange. Material paths retain
/// their own explicit nozzle, diffusion, porous and membrane models.
///
/// The original fixed-volume first-order partition and all source-domain limits
/// apply. Distinct chambers do not become one continuous acoustic junction.
/// An input-kind mismatch, final driving-force or energy refusal, or cancellation
/// leaves every accepted parcel and every original physical owner unchanged.
#[allow(clippy::too_many_arguments)]
pub fn advance_finite_gas_acoustic_transport_with_paths<C: FnMut() -> bool>(
    parcels: &mut BTreeMap<u64, ClosedGasParcel>,
    nodes: &mut BTreeMap<u64, FiniteGasAcousticNode>,
    inputs: &BTreeMap<u64, FiniteGasAcousticInput>,
    paths: GasTransportPaths<'_>,
    duration_s: f64,
    controls: FiniteGasDgTransportSolve,
    cancelled: &mut C,
) -> Result<FiniteGasAcousticTransportStep, FiniteGasDgTransportError> {
    advance_finite_gas_transport_with_paths(
        parcels, nodes, inputs, paths, duration_s, controls, cancelled,
    )
}

/// Preview exactly the same mixed-owner transaction with no accepted mutation.
#[allow(clippy::too_many_arguments)]
pub fn finite_gas_acoustic_transport_candidate_with_paths<C: FnMut() -> bool>(
    parcels: &BTreeMap<u64, ClosedGasParcel>,
    nodes: &BTreeMap<u64, FiniteGasAcousticNode>,
    inputs: &BTreeMap<u64, FiniteGasAcousticInput>,
    paths: GasTransportPaths<'_>,
    duration_s: f64,
    controls: FiniteGasDgTransportSolve,
    cancelled: &mut C,
) -> Result<FiniteGasAcousticTransportCandidate, FiniteGasDgTransportError> {
    finite_gas_transport_candidate_with_paths(
        parcels, nodes, inputs, paths, duration_s, controls, cancelled,
    )
}
