//! One native molecular snapshot and final source admission for both heat laws.

use super::*;
use fs_gas::material::transport::parcel::ParcelMolecularError;
use fs_material::thermochemical_gas::energy::parcel::transfer::snapshot::MAX_SNAPSHOT_PARCELS;
use fs_material::thermochemical_gas::energy::parcel::transfer::species::network::transport::GasTransportError;
use fs_material::thermochemical_gas::energy::{GasEnergySolve, parcel::ClosedGasParcel};

pub(super) enum MolecularStageError {
    Invalid(&'static str),
    Molecular(ParcelMolecularError),
    Transport(GasTransportError),
    Cancelled,
}
impl From<ParcelMolecularError> for MolecularStageError {
    fn from(error: ParcelMolecularError) -> Self {
        match error {
            ParcelMolecularError::Cancelled => Self::Cancelled,
            other => Self::Molecular(other),
        }
    }
}
impl From<GasTransportError> for MolecularStageError {
    fn from(error: GasTransportError) -> Self {
        match error {
            GasTransportError::Cancelled => Self::Cancelled,
            other => Self::Transport(other),
        }
    }
}
impl From<MolecularStageError> for MaterialFlowError {
    fn from(error: MolecularStageError) -> Self {
        match error {
            MolecularStageError::Invalid(reason) => Self::Invalid(reason),
            MolecularStageError::Molecular(error) => Self::Molecular(error),
            MolecularStageError::Transport(error) => Self::Transport(error),
            MolecularStageError::Cancelled => Self::Cancelled,
        }
    }
}
impl<E> From<MolecularStageError> for MaterialGasLawFlowError<E> {
    fn from(error: MolecularStageError) -> Self {
        match error {
            MolecularStageError::Invalid(reason) => Self::Invalid(reason),
            MolecularStageError::Molecular(error) => Self::Molecular(error),
            MolecularStageError::Transport(error) => Self::Transport(error),
            MolecularStageError::Cancelled => Self::Cancelled,
        }
    }
}
fn check(cx: &Cx<'_>) -> Result<(), MolecularStageError> {
    cx.checkpoint().map_err(|_| MolecularStageError::Cancelled)
}

pub(super) struct MolecularInputs {
    pub(super) duration_s: f64,
    pub(super) parcels: BTreeMap<u64, ClosedGasParcel>,
    pub(super) endpoints: BTreeMap<u64, GasEnergySolve>,
}

pub(super) fn inputs(
    state: &MaterialExchangeState<'_>,
    cx: &Cx<'_>,
    target_s: f64,
    drives: &[ChamberDrive],
    controls: NetworkSolve,
) -> Result<MolecularInputs, MolecularStageError> {
    check(cx)?;
    let duration_s = target_s - state.time_s();
    if !duration_s.is_finite()
        || duration_s <= 0.0
        || state.chambers.is_empty()
        || state.chambers.len() > controls.maximum_chambers
        || state.chambers.len() > MAX_SNAPSHOT_PARCELS
        || drives.len() != state.chambers.len()
    {
        return Err(MolecularStageError::Invalid(
            "one forward clock and one gas solver per bounded live chamber required",
        ));
    }
    let mut endpoints = BTreeMap::new();
    for drive in drives {
        check(cx)?;
        if !state.chambers.contains_key(&drive.chamber)
            || endpoints
                .insert(drive.chamber, drive.piston_solve.gas)
                .is_some()
        {
            return Err(MolecularStageError::Invalid(
                "duplicate drive or unknown live chamber",
            ));
        }
    }
    let mut parcels = BTreeMap::new();
    for (&id, chamber) in &state.chambers {
        check(cx)?;
        parcels.insert(id, chamber.gas().clone());
    }
    check(cx)?;
    Ok(MolecularInputs {
        duration_s,
        parcels,
        endpoints,
    })
}

impl<'m> MaterialExchangeState<'m> {
    #[allow(clippy::too_many_arguments)]
    pub(super) fn molecular_stage<L: TransportLaw>(
        &self,
        cx: &Cx<'_>,
        target_s: f64,
        drives: &[ChamberDrive],
        molecular: &mut ParcelMolecularTransport<L>,
        links: &[ParcelMolecularLink],
        paths: GasTransportPaths<'_>,
        molecular_controls: ParcelMolecularSolve,
        controls: NetworkSolve,
    ) -> Result<FlowStage<ParcelMolecularStep>, MolecularStageError> {
        let MolecularInputs {
            duration_s,
            parcels,
            endpoints,
        } = inputs(self, cx, target_s, drives, controls)?;
        let (parcels, flow) = molecular
            .step_candidate(
                &parcels,
                links,
                paths,
                duration_s,
                &endpoints,
                molecular_controls,
                || cx.checkpoint().is_err(),
            )?
            .into_parts();
        let mut moved: BTreeSet<_> = flow
            .transport
            .exchange
            .parcels
            .iter()
            .filter(|p| p.incoming_mass_kg > 0.0 || p.outgoing_mass_kg > 0.0)
            .map(|p| p.id)
            .collect();
        // Net cancellation cannot restore the closed two-port contract when
        // one of the original faces carries an additional heat stream.
        for link in &flow.links {
            check(cx)?;
            if link.heat_correction_left_to_right_j != 0.0 {
                moved.extend(link.parcels);
            }
        }
        for heat in &flow.transport.exchange.heat {
            check(cx)?;
            if heat.heat_left_to_right_j != 0.0 {
                moved.extend(heat.parcels);
            }
        }
        check(cx)?;
        Ok(FlowStage {
            parcels,
            moved,
            report: flow,
            mass_tolerance_kg: molecular_controls.physical.exchange.mass_tolerance_kg,
        })
    }

    /// Validate completed feedback without sampling or spending another flux.
    pub(super) fn admit_molecular_endpoint<L: TransportLaw>(
        &self,
        cx: &Cx<'_>,
        molecular: &mut ParcelMolecularTransport<L>,
        report: &ParcelMolecularStep,
        controls: ParcelMolecularSolve,
    ) -> Result<(), MolecularStageError> {
        let mut endpoints = BTreeMap::new();
        for (&id, chamber) in &self.chambers {
            check(cx)?;
            endpoints.insert(id, chamber.gas().clone());
        }
        report
            .transport
            .admit_endpoint_driving_forces(&endpoints, || cx.checkpoint().is_err())?;
        let molecular_endpoints: BTreeSet<_> =
            report.links.iter().flat_map(|link| link.parcels).collect();
        for id in molecular_endpoints {
            check(cx)?;
            let parcel = endpoints.get(&id).ok_or(MolecularStageError::Invalid(
                "missing final molecular gas owner",
            ))?;
            molecular.qualify_parcel(parcel, controls.face, || cx.checkpoint().is_err())?;
        }
        check(cx)
    }
}
