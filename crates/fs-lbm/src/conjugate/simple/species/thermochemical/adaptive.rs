//! Refreshed gas mixing through the original material time-integration owner.

use super::*;
use fs_material::thermochemical_gas::energy::parcel::transfer::species::network::{
    SpeciesNetworkAdvanceCandidate, SpeciesNetworkAdvanceControls, SpeciesNetworkAdvanceError,
    transport::GasTransportSolve,
};

/// Complete source-bearing endpoints over an error-controlled physical window.
#[derive(Clone, Debug, PartialEq)]
pub struct TurbulentGasSpeciesAdvance {
    /// Original temporal owner: accepted-only species/heat ledgers, local error
    /// ratios, fixed whole-window conservation and source-relative endpoint audit.
    pub exchange: SpeciesNetworkAdvanceCandidate,
    /// Native species-conducting faces held fixed over this window. Additional
    /// heat paths remain separately identified in the original exchange report.
    pub conducting_faces: usize,
    /// Largest final departure from the SAME initial accepted carrier [kg].
    pub maximum_carrier_mass_defect_kg: f64,
    /// Sum of absolute final departures from that original carrier [kg].
    pub absolute_carrier_mass_defect_kg: f64,
}

/// Shared by species-only and combined heat/species public entry points.
#[allow(clippy::too_many_arguments)]
pub(super) fn advance_prepared(
    prepared: &PreparedTurbulentGasMixing<'_>,
    original: &BTreeMap<u64, ClosedGasParcel>,
    duration_s: f64,
    endpoints: &BTreeMap<u64, GasEnergySolve>,
    physical: GasTransportSolve,
    temporal: SpeciesNetworkAdvanceControls,
    mut heat: impl FnMut(
        &BTreeMap<u64, ClosedGasParcel>,
        f64,
        &mut dyn FnMut() -> bool,
    ) -> Result<Vec<SpeciesNetworkHeatTransfer>, SpeciesNetworkError>,
    gate: &CancelGate,
) -> Result<TurbulentGasSpeciesAdvance, ChtError> {
    poll(gate)?;
    let exchange = network::advance_candidate_adaptive(
        original,
        duration_s,
        endpoints,
        physical,
        temporal,
        |state, dt, cancelled| {
            let heat = heat(state, dt, cancelled)?;
            prepared
                .step_candidate(state, dt, endpoints, physical.exchange, &heat, cancelled)
                .map(|step| step.exchange)
        },
        || gate.is_requested(),
    )
    .map_err(|error| match error {
        SpeciesNetworkAdvanceError::Cancelled => ChtError::Cancelled,
        other => refuse(&format!("original adaptive gas owner refused: {other}")),
    })?;
    let (maximum_carrier_mass_defect_kg, absolute_carrier_mass_defect_kg) = prepared
        .audit_carrier(
            exchange.parcels(),
            physical.exchange.mass_tolerance_kg,
            &mut || gate.is_requested(),
        )
        .map_err(material_error)?;
    poll(gate)?;
    Ok(TurbulentGasSpeciesAdvance {
        exchange,
        conducting_faces: prepared.faces.len(),
        maximum_carrier_mass_defect_kg,
        absolute_carrier_mass_defect_kg,
    })
}

impl FvFlow {
    /// Advance turbulent source-owned gas mixing through a complete interval.
    ///
    /// Native geometry, accepted carrier masses and this flow's `nu_t / Sc_t`
    /// stay fixed for the window. Every full/half probe resamples actual species
    /// composition and donor partial enthalpy from its own original gas states.
    /// Each conservative stage performs one original caloric/EOS/transport
    /// reconstruction per affected cell. Newly arrived species can move onward
    /// during later accepted time, and no original source or flow is mutated.
    ///
    /// The EXISTING material controller compares one full step with two half
    /// steps, controlling each constituent's inventory and source temperature.
    /// Only the two fine candidates can be accepted. There is no extrapolated
    /// composition or energy state and no second-order/global-error claim. Typed
    /// gross-withdrawal limits, source thermal bracket overshoot, local error
    /// and actual endpoint entropy admission can shorten the trial. Arbitrary
    /// material/source failures remain fatal. Step/evaluation limits bound work.
    ///
    /// Accepted partial-enthalpy transfers use the existing compensated ledger.
    /// The original-to-final material audit separately checks each species,
    /// total mass, source-relative internal energy and actual endpoint entropy
    /// under the SAME physical budgets, never multiplied by the number of steps.
    /// Carrier drift is likewise measured against the original accepted masses.
    /// Formation references are canceled by the original caloric owner before
    /// arithmetic, not subtracted from large displayed enthalpy totals.
    ///
    /// The source snapshot's closed fixed-volume and 256-parcel/4096-stream
    /// limits remain. Temporal refinement does not qualify the turbulence law,
    /// replace spatial refinement, or refresh the momentum closure. The caller
    /// owns the next flow epoch. Species-only mixing adds no separate heat flux;
    /// thermal transport retains its explicit coupled entry point.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_turbulent_gas_species(
        &self,
        domain: &VoxelDomain,
        carrier: &MaterialAdvectionStep,
        original: &BTreeMap<u64, ClosedGasParcel>,
        turbulent_schmidt: f64,
        duration_s: f64,
        endpoints: &BTreeMap<u64, GasEnergySolve>,
        physical: GasTransportSolve,
        temporal: SpeciesNetworkAdvanceControls,
        gate: &CancelGate,
    ) -> Result<TurbulentGasSpeciesAdvance, ChtError> {
        let prepared = PreparedTurbulentGasMixing::new(
            self,
            domain,
            carrier,
            original,
            turbulent_schmidt,
            endpoints,
            physical.exchange,
            gate,
        )?;
        advance_prepared(
            &prepared,
            original,
            duration_s,
            endpoints,
            physical,
            temporal,
            |_, _, _| Ok(Vec::new()),
            gate,
        )
    }
}
