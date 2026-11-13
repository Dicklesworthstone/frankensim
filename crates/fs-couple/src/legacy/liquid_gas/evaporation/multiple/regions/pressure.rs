//! Original regional EOS pressure at every actual native gas link.

use std::cell::Cell;

use super::super::{FreeSurfaceFiniteGasPressureLoad, admit_pressure};
use super::{
    BTreeMap, CancelGate, ClosedGasParcel, DistributedLiquidGasSolve, Error,
    FreeSurfaceEnthalpyReport, GasLiquidChamberHeatSolve, LiquidGasEvaporationSurface,
    MaterialFreeSurfaceEnthalpy3, RegionalLiquidGasCandidate, RegionalLiquidGasHeatSurface,
    RegionalLiquidGasSurface, ThermochemicalGasModel, advance_regions, poll,
};

/// Explicit gas-region membership for the original native pressure boundary.
///
/// Keys are integer native gas-cell coordinates. The map must cover EVERY gas
/// neighbor sampled by the native D3Q19 pressure owner, including diagonal
/// links and exposed liquid cells without an authored phase-transfer patch.
/// Extra in-domain entries are allowed, but every named parcel must exist.
/// The original native owner decides which entries are actual exposed links;
/// Cartesian phase patches are not a substitute for that pressure boundary.
///
/// Supply a current map on each call. Native flow may change the exposed links;
/// membership is prescribed, never inferred from volume fractions or averaged
/// between distinct finite atmospheres. Different gas neighbors of one native
/// cell may apply different pressures. Every authored phase face must agree
/// with this map at its own gas coordinate.
#[derive(Clone, Copy, Debug)]
pub struct RegionalGasPressure<'a> {
    pub gas_regions: &'a BTreeMap<[usize; 3], u64>,
    /// Bound on the original combined gas/capillary lattice-density excursion.
    pub maximum_density_excursion: f64,
}

/// Actual regional pressure loads used by the one native physical flow step.
#[derive(Clone, Debug, PartialEq)]
pub struct RegionalGasPressureReport {
    /// Only regions sampled at actual native gas links appear here. The
    /// initial pressure precedes phase/sensible-heat exchange; the applied
    /// pressure is that same parcel's actual post-exchange EOS value.
    pub regions: BTreeMap<u64, FreeSurfaceFiniteGasPressureLoad>,
    /// Number of original native D3Q19 pressure-link callbacks for each region.
    /// These are streaming characteristics, not measured interfacial areas.
    pub sampled_links: BTreeMap<u64, usize>,
}

/// Exchange regional phase mass/enthalpy once, then apply each resulting gas
/// EOS pressure through the original native link-pressure boundary once.
///
/// The fixed original pressure reference, stationary solid geometry, fixed
/// external parcel volumes and combined density-excursion bound are required.
/// Actual post-exchange pressure is frozen during the one native flow step.
/// Parcels without phase patches can still load other interfaces through the
/// complete pressure map. Missing or contradictory coverage refuses the whole
/// private candidate without an atmospheric fallback or partial gas update.
///
/// This adds mechanical boundary loading, not gas-volume work or pressure heat.
/// The gas volumes are NOT the geometric complement of moving liquid. No gas
/// momentum, gas topology evolution, pressure-work debit or closure of total
/// mechanical plus thermal energy is inferred. Thermal/source balances retain
/// their original meaning and never turn flow truncation into caloric heat.
#[allow(clippy::too_many_arguments)]
pub fn advance_liquid_regions_with_gas_and_pressure<'s, 'm, 'catalog>(
    liquid: &MaterialFreeSurfaceEnthalpy3<'s, 'm>,
    parcels: &BTreeMap<u64, ClosedGasParcel>,
    catalog: &'catalog ThermochemicalGasModel,
    surfaces: &[RegionalLiquidGasSurface<'_>],
    controls: DistributedLiquidGasSolve,
    pressure: RegionalGasPressure<'_>,
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
        Some(pressure),
        gate,
    )
}

/// Reciprocal regional sensible heat and phase transfer followed by each
/// original post-heat/phase EOS pressure on its actual native gas links.
/// Sensible/phase heat is already in native enthalpy and is not supplied again
/// during pressure loading. All admission and no-work boundaries of
/// [`advance_liquid_regions_with_gas_and_pressure`] apply unchanged.
#[allow(clippy::too_many_arguments)]
pub fn advance_liquid_regions_with_gas_and_heat_and_pressure<'s, 'm, 'catalog>(
    liquid: &MaterialFreeSurfaceEnthalpy3<'s, 'm>,
    parcels: &BTreeMap<u64, ClosedGasParcel>,
    catalog: &'catalog ThermochemicalGasModel,
    surfaces: &[RegionalLiquidGasHeatSurface<'_>],
    controls: DistributedLiquidGasSolve,
    heat_controls: GasLiquidChamberHeatSolve,
    pressure: RegionalGasPressure<'_>,
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
        Some(pressure),
        gate,
    )
}

pub(super) fn admit(
    liquid: &MaterialFreeSurfaceEnthalpy3<'_, '_>,
    parcels: &BTreeMap<u64, ClosedGasParcel>,
    pressure: RegionalGasPressure<'_>,
    gate: &CancelGate,
) -> Result<(), Error> {
    admit_pressure(liquid, Some(pressure.maximum_density_excursion))?;
    if pressure.gas_regions.is_empty() {
        return Err(Error::Invalid(
            "regional gas pressure requires a complete native gas-coordinate map",
        ));
    }
    let domain = liquid.flow().lattice().domain_cells();
    for (coordinates, id) in pressure.gas_regions {
        poll(gate)?;
        if !parcels.contains_key(id) {
            return Err(Error::Invalid(
                "regional gas pressure map names an absent finite gas parcel",
            ));
        }
        if (0..3).any(|axis| {
            coordinates[axis] >= domain[axis] || i64::try_from(coordinates[axis]).is_err()
        }) {
            return Err(Error::Invalid(
                "regional gas pressure coordinates must be representable and in-domain",
            ));
        }
    }
    Ok(())
}

pub(super) fn admit_surface(
    pressure: RegionalGasPressure<'_>,
    region_id: u64,
    surface: &LiquidGasEvaporationSurface<'_>,
) -> Result<(), Error> {
    match pressure.gas_regions.get(&surface.gas) {
        Some(&id) if id == region_id => Ok(()),
        Some(_) => Err(Error::Invalid(
            "regional phase face and native gas pressure map name different parcels",
        )),
        None => Err(Error::Invalid(
            "regional gas pressure map omits an authored phase face",
        )),
    }
}

struct SampledPressure {
    load: FreeSurfaceFiniteGasPressureLoad,
    links: Cell<usize>,
}

pub(super) fn advance_flow(
    liquid: &mut MaterialFreeSurfaceEnthalpy3<'_, '_>,
    initial_parcels: &BTreeMap<u64, ClosedGasParcel>,
    parcels: &BTreeMap<u64, ClosedGasParcel>,
    pressure: RegionalGasPressure<'_>,
    gate: &CancelGate,
) -> Result<(FreeSurfaceEnthalpyReport, RegionalGasPressureReport), Error> {
    poll(gate)?;
    if initial_parcels.len() != parcels.len() {
        return Err(Error::Invalid(
            "regional pressure requires the original fixed gas parcel roster",
        ));
    }
    for (id, initial) in initial_parcels {
        poll(gate)?;
        if parcels
            .get(id)
            .is_none_or(|gas| gas.volume_m3() != initial.volume_m3())
        {
            return Err(Error::Invalid(
                "regional phase exchange changed an original fixed gas volume",
            ));
        }
    }
    let mut loads = BTreeMap::new();
    for &id in pressure.gas_regions.values() {
        poll(gate)?;
        if loads.contains_key(&id) {
            continue;
        }
        let gas = parcels.get(&id).ok_or(Error::Invalid(
            "regional pressure lost an original mapped gas parcel",
        ))?;
        let initial = &initial_parcels[&id];
        loads.insert(
            id,
            SampledPressure {
                load: FreeSurfaceFiniteGasPressureLoad {
                    initial_pressure_pa: initial.state().state().pressure,
                    applied_pressure_pa: gas.state().state().pressure,
                    volume_m3: gas.volume_m3(),
                    maximum_density_excursion: pressure.maximum_density_excursion,
                },
                links: Cell::new(0),
            },
        );
    }
    // The original material callback accepts f64, not Result. NaN makes that
    // original owner refuse atomically; retain our precise coverage refusal
    // for the caller instead of replacing a missing region with ambient gas.
    let coverage_error = Cell::new(None);
    let stepped = liquid.step_with_heat_and_gas_link_pressure(
        &[],
        gate,
        pressure.maximum_density_excursion,
        |link| {
            let [Ok(x), Ok(y), Ok(z)] = link.link.gas_cell.map(usize::try_from) else {
                coverage_error.set(Some(
                    "native regional pressure gas coordinate is unrepresentable",
                ));
                return f64::NAN;
            };
            let Some(id) = pressure.gas_regions.get(&[x, y, z]) else {
                coverage_error.set(Some(
                    "regional gas pressure map omits an exposed native gas link",
                ));
                return f64::NAN;
            };
            let sample = &loads[id];
            let Some(count) = sample.links.get().checked_add(1) else {
                coverage_error.set(Some(
                    "regional native pressure-link count is unrepresentable",
                ));
                return f64::NAN;
            };
            sample.links.set(count);
            sample.load.applied_pressure_pa
        },
    );
    poll(gate)?;
    if let Some(reason) = coverage_error.get() {
        return Err(Error::Invalid(reason));
    }
    let flow = stepped?;
    let mut regions = BTreeMap::new();
    let mut sampled_links = BTreeMap::new();
    for (id, sample) in loads {
        poll(gate)?;
        let count = sample.links.get();
        if count != 0 {
            regions.insert(id, sample.load);
            sampled_links.insert(id, count);
        }
    }
    Ok((
        flow,
        RegionalGasPressureReport {
            regions,
            sampled_links,
        },
    ))
}
