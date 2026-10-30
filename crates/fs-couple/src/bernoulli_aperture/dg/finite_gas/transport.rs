//! Original moving sources on the shared finite-gas transport transaction.
//!
//! The material owner transports species and partial enthalpies once. Each
//! original aperture then advances its mechanics and complete acoustic load,
//! paying native gas-parameter work from that same transported candidate. No
//! valve flow is silently converted into thermodynamic mass transport.

use super::{FiniteGasApertureDgStep, acoustic, admit};
use crate::bernoulli_aperture::dg::LoadedApertureDg;
use crate::bernoulli_aperture::tube::TubeDrive;
use crate::bernoulli_aperture::two_sided::LoadedAperture;
use crate::bernoulli_aperture::wall::material::{
    FiniteGasDgTransportCandidate, FiniteGasDgTransportStep, FiniteGasTransportNode,
    advance_finite_gas_transport_with_paths, finite_gas_transport_candidate_with_paths,
};
pub use crate::dg_port::finite_gas::transport::{
    FiniteGasDgTransportError, FiniteGasDgTransportInput, FiniteGasDgTransportSolve,
};
use crate::dg_port::finite_gas::{FiniteGasDgBackground, FiniteGasDgError};
use crate::dg_port::sum_terminal_work;
use fs_material::thermochemical_gas::energy::parcel::ClosedGasParcel;
use fs_material::thermochemical_gas::energy::parcel::composition::GasCompositionStep;
use fs_material::thermochemical_gas::energy::parcel::transfer::{
    nozzle::network::snapshot::SnapshotNozzleLink,
    species::{diffusion::network::DiffusionNetworkLink, network::transport::GasTransportPaths},
};
use std::collections::BTreeMap;

/// Actual held physical drives of the original nonlinear mechanical source.
#[derive(Clone, Copy, Debug)]
pub struct FiniteGasApertureDrive {
    /// Prescribed upstream pressure and independently supplied body volume flow.
    pub drive: TubeDrive,
    /// Original blade forces or material-modal force-projection inputs [N].
    /// A scalar aperture requires the second force to be exactly zero.
    pub forces_n: [f64; 2],
    /// Source-bound modal pressure tractions [Pa]; zero for non-modal sources.
    /// Their work is already included in the original mechanical work account.
    pub pressures_pa: [f64; 2],
}

/// One finite parcel's caloric controls and optional moving-source drive.
#[derive(Clone, Debug)]
pub struct FiniteGasApertureTransportInput {
    /// Same native heat/work controls as the wall and passive-load adapters.
    /// `ports` must be empty: the original valve computes its own load flow.
    pub thermal: FiniteGasDgTransportInput,
    /// Present exactly for parcels with a moving source; absent for reservoirs.
    pub source: Option<FiniteGasApertureDrive>,
}

/// Original species transport, moving-source reports and one enclosing account.
/// The inherited `acoustic_*` accounts include the original mechanical owner.
pub type FiniteGasApertureTransportStep =
    FiniteGasDgTransportStep<FiniteGasApertureDgStep<GasCompositionStep>>;

/// Private gas, mechanical, material, acoustic and terminal endpoint histories.
pub type FiniteGasApertureTransportCandidate =
    FiniteGasDgTransportCandidate<LoadedApertureDg, FiniteGasApertureDgStep<GasCompositionStep>>;

impl FiniteGasTransportNode for LoadedApertureDg {
    type Step<S> = FiniteGasApertureDgStep<S>;
    type Input = FiniteGasApertureTransportInput;

    fn transport_input(
        input: &Self::Input,
        bound: bool,
    ) -> Result<&FiniteGasDgTransportInput, &'static str> {
        if input.source.is_some() != bound || !input.thermal.ports.is_empty() {
            return Err(
                "finite aperture transport requires exact source bindings and no prescribed acoustic ports",
            );
        }
        Ok(&input.thermal)
    }

    fn dt(&self) -> f64 {
        self.load().propagator().dt()
    }

    fn step_count(&self) -> u64 {
        self.duct_state().step_count()
    }

    fn stored_energy_j(&self) -> f64 {
        LoadedApertureDg::stored_energy_j(self)
    }

    fn admit<E, C: FnMut() -> bool>(
        &self,
        gas: &ClosedGasParcel,
        input: &Self::Input,
        cancelled: &mut C,
    ) -> Result<(), FiniteGasDgError<E>> {
        let source = input.source.ok_or(FiniteGasDgError::Invalid(
            "finite aperture transport is missing its original source drive",
        ))?;
        if self.aperture().accepted_steps() != self.duct_state().step_count()
            || self.aperture().time_step_s().to_bits() != self.dt().to_bits()
            || !source.drive.upstream_pressure_pa.is_finite()
            || !source.drive.body_flow_m3_s.is_finite()
            || !source.forces_n.iter().all(|v| v.is_finite())
            || !source.pressures_pa.iter().all(|v| v.is_finite())
            || matches!(self.aperture(), LoadedAperture::Single(_) if source.forces_n[1] != 0.0)
        {
            return Err(FiniteGasDgError::Invalid(
                "finite aperture transport requires one mechanical/acoustic clock and admitted finite source drives",
            ));
        }
        self.aperture()
            .validate_pressures(source.pressures_pa)
            .map_err(acoustic)?;
        admit(
            self,
            gas,
            input.thermal.acoustic.ok_or(FiniteGasDgError::Invalid(
                "finite aperture transport is missing acoustic work controls",
            ))?,
            cancelled,
        )
    }

    fn advance<S, E, C: FnMut() -> bool>(
        &mut self,
        gas: &mut ClosedGasParcel,
        input: &Self::Input,
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
                "moving apertures do not own the source-wall molecular heat integral",
            ));
        }
        let source = input.source.expect("admitted moving-source input");
        self.step_with_finite_gas_controlled(
            gas,
            source.drive,
            source.forces_n,
            source.pressures_pa,
            input.thermal.acoustic.expect("admitted acoustic controls"),
            evolve,
            cancelled,
        )
    }

    fn accounts<S>(step: &Self::Step<S>) -> [f64; 4] {
        let frame = &step.coupled;
        // The valve/duct pressure-flow work is INTERNAL to this node. Only
        // original prescribed drives and signed exterior exchange cross its
        // boundary. Modal pressure work is already in mechanical_work_j.
        [
            frame.upstream_work_j + frame.body_work_j + frame.aperture.mechanical_work_j
                - sum_terminal_work(&frame.terminals),
            frame.inferred_total_dissipation_j,
            frame.parameter_work_j,
            0.0,
        ]
    }

    fn source<S>(step: &Self::Step<S>) -> &S { &step.gas.background }
}

/// Advance native nozzle/diffusion transport and complete moving sources together.
#[allow(clippy::too_many_arguments)]
pub fn advance_finite_gas_aperture_transport<C: FnMut() -> bool>(
    parcels: &mut BTreeMap<u64, ClosedGasParcel>,
    nodes: &mut BTreeMap<u64, LoadedApertureDg>,
    inputs: &BTreeMap<u64, FiniteGasApertureTransportInput>,
    nozzles: &[SnapshotNozzleLink],
    diffusion: &[DiffusionNetworkLink<'_>],
    duration_s: f64,
    controls: FiniteGasDgTransportSolve,
    cancelled: &mut C,
) -> Result<FiniteGasApertureTransportStep, FiniteGasDgTransportError> {
    advance_finite_gas_aperture_transport_with_paths(
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

/// Publish one common endpoint for native material transport, every original
/// moving aperture, its complete DG load and all finite gas parcels.
///
/// Each node owns its original scalar, independent-blade or material-modal
/// source and one downstream complete load. A parcel must cover the actual
/// main bore and its main-gas terminal passages. Independent branch gases keep
/// their original sources. All source and acoustic clocks must agree exactly.
/// Unbound reservoirs have no source, acoustic controls or physical ports.
///
/// Nozzle, diffusion, porous and membrane paths share one unchanged gross
/// species inventory. Each native source advances once in its initial medium;
/// its repeatable caloric callback pays the actual complete-load parameter
/// work from the transported candidate. Final driving-force and enclosing
/// mass/energy gates precede publication of every participant together.
///
/// This retains the existing first-order, fixed-volume partition and native
/// validity limits. Material paths are explicit physical connections, separate
/// from the valve's oscillatory acoustic flow. Mean-pressure prestress, automatic
/// aperture mass transport, two-sided density-discontinuous valve physics and
/// thermal deposition of acoustic/mechanical losses are not supplied here.
/// An existing upstream acoustic tract is refused rather than assigned the
/// downstream inventory. Any refusal or cancellation preserves every owner.
#[allow(clippy::too_many_arguments)]
pub fn advance_finite_gas_aperture_transport_with_paths<C: FnMut() -> bool>(
    parcels: &mut BTreeMap<u64, ClosedGasParcel>,
    nodes: &mut BTreeMap<u64, LoadedApertureDg>,
    inputs: &BTreeMap<u64, FiniteGasApertureTransportInput>,
    paths: GasTransportPaths<'_>,
    duration_s: f64,
    controls: FiniteGasDgTransportSolve,
    cancelled: &mut C,
) -> Result<FiniteGasApertureTransportStep, FiniteGasDgTransportError> {
    advance_finite_gas_transport_with_paths(
        parcels, nodes, inputs, paths, duration_s, controls, cancelled,
    )
}

/// Preview the exact same transaction while retaining every accepted gas,
/// source, contact, material-memory, wave and terminal history unchanged.
#[allow(clippy::too_many_arguments)]
pub fn finite_gas_aperture_transport_candidate_with_paths<C: FnMut() -> bool>(
    parcels: &BTreeMap<u64, ClosedGasParcel>,
    nodes: &BTreeMap<u64, LoadedApertureDg>,
    inputs: &BTreeMap<u64, FiniteGasApertureTransportInput>,
    paths: GasTransportPaths<'_>,
    duration_s: f64,
    controls: FiniteGasDgTransportSolve,
    cancelled: &mut C,
) -> Result<FiniteGasApertureTransportCandidate, FiniteGasDgTransportError> {
    finite_gas_transport_candidate_with_paths(
        parcels, nodes, inputs, paths, duration_s, controls, cancelled,
    )
}
