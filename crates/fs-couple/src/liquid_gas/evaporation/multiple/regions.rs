//! Distinct finite atmospheres on one original native liquid transaction.
//!
//! Each region supplies its actual temperature, pressure, vapor inventory and
//! original gas calorics. All of its native cells share that region's donor
//! account. Regions prepare against the SAME unchanged native state, followed
//! by one batched mass/enthalpy transaction and one native physical flow step.
//! No intermediate liquid endpoint or averaged atmosphere is introduced.

/// Regional post-exchange EOS pressure on the original native gas links.
pub mod pressure;

use super::{
    BTreeMap, CancelGate, ClosedGasParcel, DistributedLiquidGasSolve, Error,
    FreeSurfaceEnthalpyReport, FreeSurfaceMassExchangeReport, GasEvaporationError,
    GasFilmHeatSolve, GasLiquidChamberHeatSolve, GasLiquidChamberHeatStep, GasLiquidChamberStep,
    LiquidGasEvaporationSurface, LiquidGasHeatSurface, MaterialFreeSurfaceEnthalpy3,
    NativeLiquidGasSurfaceStep, ThermochemicalGasModel, Total, account_native_phase, poll,
    prepare_native_phase,
};

/// An original measured gas face and the finite parcel occupying its region.
#[derive(Clone, Copy)]
pub struct RegionalLiquidGasSurface<'a> {
    pub region_id: u64,
    pub surface: LiquidGasEvaporationSurface<'a>,
}

/// An original phase/sensible-heat face and its finite gas region.
#[derive(Clone, Copy)]
pub struct RegionalLiquidGasHeatSurface<'a> {
    pub region_id: u64,
    pub surface: LiquidGasHeatSurface<'a>,
}

/// One region's actual phase, heat and native-storage account BEFORE flow.
/// Native flow may subsequently carry liquid enthalpy between these regions;
/// its single complete account appears on `RegionalLiquidGasReport`.
#[derive(Clone, Debug)]
pub struct RegionalLiquidGasRegionReport {
    pub phase: GasLiquidChamberStep,
    pub exchange: GasLiquidChamberStep,
    pub heat: Option<GasLiquidChamberHeatStep>,
    /// One original inventory per liquid cell, in Cartesian coordinate order.
    pub surfaces: Vec<NativeLiquidGasSurfaceStep>,
    pub native_phase_caloric_defect_j: f64,
    pub mass_residual_kg: f64,
    pub volatile_mass_residual_kg: f64,
    pub absolute_mass_defect_kg: f64,
    pub energy_residual_j: f64,
    pub absolute_energy_defect_j: f64,
}

/// All active regional source accounts and the sole native physical step.
#[derive(Clone, Debug)]
pub struct RegionalLiquidGasReport {
    pub duration_s: f64,
    /// Only referenced regions appear here; unused input parcels are retained.
    pub regions: BTreeMap<u64, RegionalLiquidGasRegionReport>,
    pub mass_exchange: FreeSurfaceMassExchangeReport,
    pub flow: FreeSurfaceEnthalpyReport,
    /// Actual regional EOS loads sampled by the original native flow owner.
    /// Absent on the original unpressurized entry points.
    pub pressure: Option<pressure::RegionalGasPressureReport>,
    pub native_liquid_enthalpy_change_j: f64,
    pub native_phase_caloric_defect_j: f64,
    pub mass_residual_kg: f64,
    pub volatile_mass_residual_kg: f64,
    pub absolute_mass_defect_kg: f64,
    pub energy_residual_j: f64,
    pub absolute_energy_defect_j: f64,
}

pub struct RegionalLiquidGasCandidate<'s, 'm, 'catalog> {
    liquid: MaterialFreeSurfaceEnthalpy3<'s, 'm>,
    parcels: BTreeMap<u64, ClosedGasParcel>,
    catalog: &'catalog ThermochemicalGasModel,
    report: RegionalLiquidGasReport,
}

impl<'s, 'm, 'catalog> RegionalLiquidGasCandidate<'s, 'm, 'catalog> {
    pub const fn liquid(&self) -> &MaterialFreeSurfaceEnthalpy3<'s, 'm> {
        &self.liquid
    }
    pub const fn parcels(&self) -> &BTreeMap<u64, ClosedGasParcel> {
        &self.parcels
    }
    pub const fn catalog(&self) -> &'catalog ThermochemicalGasModel {
        self.catalog
    }
    pub const fn report(&self) -> &RegionalLiquidGasReport {
        &self.report
    }
    pub fn into_parts(
        self,
    ) -> (
        MaterialFreeSurfaceEnthalpy3<'s, 'm>,
        BTreeMap<u64, ClosedGasParcel>,
        &'catalog ThermochemicalGasModel,
        RegionalLiquidGasReport,
    ) {
        (self.liquid, self.parcels, self.catalog, self.report)
    }
}

/// Exchange each native interface with its actual finite atmosphere, then
/// advance original native flow once. All faces of one native cell must name
/// ONE region, and a gas coordinate cannot belong to two different regions.
/// One original volatile identity/datum applies across this pure-liquid owner.
///
/// The original multiface phase owner keeps distinct nonlinear film rates;
/// cells never acquire duplicate inventories. Every region uses the unchanged
/// input liquid, its own original gas and one joint regional donor/caloric
/// solve. Source and native defects are bounded per region AND over the whole
/// interval, so opposing defects cannot cancel between regional atmospheres.
/// The total face count is bounded by `controls.chamber.maximum_surfaces`.
///
/// This is an explicit first-order phase/flow split with prescribed regional
/// membership and fixed external gas volumes. The original native boundary
/// loading remains unchanged: this does not infer regional pressure traction,
/// inter-region gas transport, gas momentum, volume work or spatial gas CFD.
/// Interface exhaustion/overflow and unsupported attached native owners retain
/// their original refusals. No partial region, gas map or liquid is published.
#[allow(clippy::too_many_arguments)]
pub fn advance_liquid_regions_with_gas<'s, 'm, 'catalog>(
    liquid: &MaterialFreeSurfaceEnthalpy3<'s, 'm>,
    parcels: &BTreeMap<u64, ClosedGasParcel>,
    catalog: &'catalog ThermochemicalGasModel,
    surfaces: &[RegionalLiquidGasSurface<'_>],
    controls: DistributedLiquidGasSolve,
    gate: &CancelGate,
) -> Result<RegionalLiquidGasCandidate<'s, 'm, 'catalog>, Error> {
    advance_regions(
        liquid,
        parcels,
        catalog,
        surfaces
            .iter()
            .map(|row| (row.region_id, row.surface, None)),
        controls,
        None,
        None,
        gate,
    )
}

/// Reciprocal sensible heat and phase transfer use each region's actual gas.
/// Each original regional chamber first exchanges its sensible heat, then
/// evaluates phase at those heated states. Only the resulting combined native
/// increments enter the single batch; heat is not applied again during flow.
#[allow(clippy::too_many_arguments)]
pub fn advance_liquid_regions_with_gas_and_heat<'s, 'm, 'catalog>(
    liquid: &MaterialFreeSurfaceEnthalpy3<'s, 'm>,
    parcels: &BTreeMap<u64, ClosedGasParcel>,
    catalog: &'catalog ThermochemicalGasModel,
    surfaces: &[RegionalLiquidGasHeatSurface<'_>],
    controls: DistributedLiquidGasSolve,
    heat_controls: GasLiquidChamberHeatSolve,
    gate: &CancelGate,
) -> Result<RegionalLiquidGasCandidate<'s, 'm, 'catalog>, Error> {
    advance_regions(
        liquid,
        parcels,
        catalog,
        surfaces.iter().map(|row| {
            (
                row.region_id,
                row.surface.surface,
                Some((
                    row.surface.heat_transfer_coefficient_w_m2_k,
                    row.surface.heat_controls,
                )),
            )
        }),
        controls,
        Some(heat_controls),
        None,
        gate,
    )
}

#[allow(clippy::too_many_arguments, clippy::type_complexity)]
fn advance_regions<'s, 'm, 'catalog, 'surface>(
    liquid: &MaterialFreeSurfaceEnthalpy3<'s, 'm>,
    parcels: &BTreeMap<u64, ClosedGasParcel>,
    catalog: &'catalog ThermochemicalGasModel,
    surfaces: impl ExactSizeIterator<
        Item = (
            u64,
            LiquidGasEvaporationSurface<'surface>,
            Option<(f64, GasFilmHeatSolve)>,
        ),
    >,
    controls: DistributedLiquidGasSolve,
    heat_controls: Option<GasLiquidChamberHeatSolve>,
    pressure_load: Option<pressure::RegionalGasPressure<'_>>,
    gate: &CancelGate,
) -> Result<RegionalLiquidGasCandidate<'s, 'm, 'catalog>, Error> {
    poll(gate)?;
    if let Some(load) = pressure_load {
        pressure::admit(liquid, parcels, load, gate)?;
    }
    if surfaces.len() == 0
        || surfaces.len() > controls.chamber.maximum_surfaces
        || [controls.mass_tolerance_kg, controls.energy_tolerance_j]
            .into_iter()
            .any(|x| !x.is_finite() || x < 0.0)
    {
        return Err(Error::Invalid(
            "regional phase needs a bounded positive face roster and finite balance budgets",
        ));
    }
    let mut regions = BTreeMap::<u64, Vec<_>>::new();
    let mut cell_regions = BTreeMap::new();
    let mut gas_regions = BTreeMap::new();
    let mut binding = None;
    for (id, surface, heat) in surfaces {
        poll(gate)?;
        if !parcels.contains_key(&id) {
            return Err(Error::Invalid(
                "native face names an absent finite gas region",
            ));
        }
        if let Some(load) = pressure_load {
            pressure::admit_surface(load, id, &surface)?;
        }
        let common = *binding.get_or_insert(surface.binding);
        if surface.binding != common {
            return Err(Error::Invalid(
                "native regions require one volatile species and liquid/gas datum binding",
            ));
        }
        if cell_regions
            .insert(surface.liquid, id)
            .is_some_and(|old| old != id)
            || gas_regions
                .insert(surface.gas, id)
                .is_some_and(|old| old != id)
        {
            return Err(Error::Invalid(
                "each native liquid cell and exposed gas coordinate requires one region",
            ));
        }
        regions.entry(id).or_default().push((surface, heat));
    }
    let binding = binding.expect("positive regional face count admitted");
    let duration_s = liquid.flow().scaling().timestep_s();
    let mut proposals = BTreeMap::new();
    let mut transfers = Vec::with_capacity(cell_regions.len());
    for (id, surfaces) in regions {
        poll(gate)?;
        let proposal = prepare_native_phase(
            liquid,
            &parcels[&id],
            catalog,
            surfaces.into_iter(),
            controls,
            heat_controls,
            gate,
        )?;
        transfers.extend_from_slice(&proposal.0.transfers);
        proposals.insert(id, proposal);
    }

    // No region has changed the native owner. Apply all selected transfers on
    // one private clone, with one material refresh and one native donor per cell.
    poll(gate)?;
    let mut next_liquid = liquid.clone();
    let mass_exchange = next_liquid.exchange_interface_mass(&transfers, controls.native, gate)?;
    let mut samples = BTreeMap::new();
    for sample in &mass_exchange.samples {
        poll(gate)?;
        samples.insert(sample.coordinates, *sample);
    }
    if samples.len() != cell_regions.len() {
        return Err(Error::Invalid(
            "regional native transaction omitted a liquid cell",
        ));
    }
    let mut next_parcels = BTreeMap::new();
    for (&id, gas) in parcels {
        poll(gate)?;
        if !proposals.contains_key(&id) {
            next_parcels.insert(id, gas.clone());
        }
    }
    let mut reports = BTreeMap::new();
    let mut gas_mass_change = Total::default();
    let mut volatile_mass_change = Total::default();
    let mut absolute_mass = Total::default();
    let mut energy = Total::default();
    let mut absolute_energy = Total::default();
    let mut native_caloric = Total::default();
    absolute_mass.add(mass_exchange.absolute_mass_defect_kg)?;
    for (id, (mut prepared, gas)) in proposals {
        poll(gate)?;
        let account = account_native_phase(&mut prepared, &mut samples, |_, _| Ok(()), gate)?;
        let mass_change = gas.mass_kg() - parcels[&id].mass_kg();
        let vapor_change = gas
            .species_mass_kg(&binding.species)
            .map_err(GasEvaporationError::from)?
            - parcels[&id]
                .species_mass_kg(&binding.species)
                .map_err(GasEvaporationError::from)?;
        let local_mass = account.native_mass_change_kg + mass_change;
        let local_vapor = account.native_mass_change_kg + vapor_change;
        let source_mass_defect = prepared.gas_mass_defect.abs()
            + prepared.gas_species_defect.abs()
            + account.absolute_mass_consistency_kg;
        let local_absolute_mass = source_mass_defect + account.native_projection_defect_kg;
        let local_energy = prepared.exchange.energy_residual_j + account.native_caloric_defect_j;
        let local_absolute_energy =
            prepared.exchange.sum_absolute_energy_residual_j + account.absolute_caloric_defect_j;
        gate_balances(
            local_mass,
            local_vapor,
            local_absolute_mass,
            local_energy,
            local_absolute_energy,
            controls,
        )?;
        gas_mass_change.add(mass_change)?;
        volatile_mass_change.add(vapor_change)?;
        absolute_mass.add(source_mass_defect)?;
        energy.add(local_energy)?;
        absolute_energy.add(local_absolute_energy)?;
        native_caloric.add(account.native_caloric_defect_j)?;
        reports.insert(
            id,
            RegionalLiquidGasRegionReport {
                phase: prepared.phase,
                exchange: prepared.exchange,
                heat: prepared.film_heat,
                surfaces: account.surfaces,
                native_phase_caloric_defect_j: account.native_caloric_defect_j,
                mass_residual_kg: local_mass,
                volatile_mass_residual_kg: local_vapor,
                absolute_mass_defect_kg: local_absolute_mass,
                energy_residual_j: local_energy,
                absolute_energy_defect_j: local_absolute_energy,
            },
        );
        next_parcels.insert(id, gas);
    }
    if !samples.is_empty() {
        return Err(Error::Invalid(
            "regional source omitted a native transaction receipt",
        ));
    }

    poll(gate)?;
    let (flow, pressure) = if let Some(load) = pressure_load {
        let (flow, report) =
            pressure::advance_flow(&mut next_liquid, parcels, &next_parcels, load, gate)?;
        (flow, Some(report))
    } else {
        (next_liquid.step(gate)?, None)
    };
    let flow_mass = flow.mass_after_kg - flow.mass_before_kg;
    let native_mass = mass_exchange.mass_change_kg + flow_mass;
    let mass_residual_kg = native_mass + gas_mass_change.value();
    let volatile_mass_residual_kg = native_mass + volatile_mass_change.value();
    absolute_mass.add(flow_mass.abs())?;
    energy.add(flow.balance_defect_j)?;
    absolute_energy.add(flow.balance_defect_j.abs())?;
    let energy_residual_j = energy.value();
    let absolute_energy_defect_j = absolute_energy.value();
    let absolute_mass_defect_kg = absolute_mass.value();
    gate_balances(
        mass_residual_kg,
        volatile_mass_residual_kg,
        absolute_mass_defect_kg,
        energy_residual_j,
        absolute_energy_defect_j,
        controls,
    )?;
    let native_liquid_enthalpy_change_j =
        mass_exchange.enthalpy_change_j + (flow.enthalpy_after_j - flow.enthalpy_before_j);
    let native_phase_caloric_defect_j = native_caloric.value();
    if !native_liquid_enthalpy_change_j.is_finite() || !native_phase_caloric_defect_j.is_finite() {
        return Err(Error::Invalid(
            "regional native energy account is unrepresentable",
        ));
    }
    poll(gate)?;
    Ok(RegionalLiquidGasCandidate {
        liquid: next_liquid,
        parcels: next_parcels,
        catalog,
        report: RegionalLiquidGasReport {
            duration_s,
            regions: reports,
            mass_exchange,
            flow,
            pressure,
            native_liquid_enthalpy_change_j,
            native_phase_caloric_defect_j,
            mass_residual_kg,
            volatile_mass_residual_kg,
            absolute_mass_defect_kg,
            energy_residual_j,
            absolute_energy_defect_j,
        },
    })
}

fn gate_balances(
    mass: f64,
    volatile: f64,
    absolute_mass: f64,
    energy: f64,
    absolute_energy: f64,
    controls: DistributedLiquidGasSolve,
) -> Result<(), Error> {
    if ![mass, volatile, absolute_mass, energy, absolute_energy]
        .into_iter()
        .all(f64::is_finite)
    {
        return Err(Error::Invalid("regional phase balance is unrepresentable"));
    }
    let worst_mass = mass.abs().max(volatile.abs()).max(absolute_mass);
    if worst_mass > controls.mass_tolerance_kg {
        return Err(Error::MassBalance {
            residual_kg: worst_mass,
            tolerance_kg: controls.mass_tolerance_kg,
        });
    }
    let worst_energy = energy.abs().max(absolute_energy);
    if worst_energy > controls.energy_tolerance_j {
        return Err(Error::EnergyBalance {
            residual_j: worst_energy,
            tolerance_j: controls.energy_tolerance_j,
        });
    }
    Ok(())
}
