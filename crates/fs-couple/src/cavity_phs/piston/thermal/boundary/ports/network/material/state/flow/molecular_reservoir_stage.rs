//! One original bath/internal transport stage and final source gate for both heat laws.

use super::super::stage::{MolecularStageError, inputs};
use super::*;

fn check(cx: &Cx<'_>) -> Result<(), MolecularStageError> {
    cx.checkpoint().map_err(|_| MolecularStageError::Cancelled)
}

impl<'m> MaterialExchangeState<'m> {
    #[allow(clippy::too_many_arguments)]
    pub(super) fn molecular_reservoir_stage<L: TransportLaw>(
        &self,
        cx: &Cx<'_>,
        target_s: f64,
        drives: &[ChamberDrive],
        molecular: &mut ParcelMolecularTransport<L>,
        internal_links: &[ParcelMolecularLink],
        reservoirs: &[ReservoirMolecularLink<'_>],
        molecular_controls: ParcelMolecularSolve,
        controls: NetworkSolve,
    ) -> Result<(BTreeMap<u64, GasPiston>, ReservoirMolecularStep), MolecularStageError> {
        let inputs = inputs(self, cx, target_s, drives, controls)?;
        let (parcels, flow) = molecular
            .network_candidate(
                &inputs.parcels,
                internal_links,
                reservoirs,
                inputs.duration_s,
                &inputs.endpoints,
                molecular_controls,
                || cx.checkpoint().is_err(),
            )?
            .into_parts();
        let mut staged = BTreeMap::new();
        for (id, parcel) in parcels {
            check(cx)?;
            let mut chamber = self
                .chambers
                .get(&id)
                .ok_or(MolecularStageError::Invalid(
                    "reservoir transport returned an unknown chamber",
                ))?
                .clone();
            chamber.gas = parcel;
            staged.insert(id, chamber);
        }
        Ok((staged, flow))
    }

    pub(super) fn admit_molecular_reservoir_endpoint<L: TransportLaw>(
        &self,
        cx: &Cx<'_>,
        molecular: &mut ParcelMolecularTransport<L>,
        flow: &ReservoirMolecularStep,
        controls: ParcelMolecularSolve,
    ) -> Result<(), MolecularStageError> {
        let connected: BTreeSet<_> = flow
            .internal_links
            .iter()
            .flat_map(|link| link.parcels)
            .chain(flow.links.iter().map(|link| link.parcel))
            .collect();
        for id in connected {
            check(cx)?;
            let gas = self.chambers.get(&id).ok_or(MolecularStageError::Invalid(
                "missing final reservoir-connected molecular owner",
            ))?;
            molecular.qualify_parcel(gas.gas(), controls.face, || cx.checkpoint().is_err())?;
        }
        check(cx)
    }
}
