//! Gas-side heat laws share the original circuit/contact/radiation assembly.
//! This boundary changes only the constitutive response, not the integrator,
//! source selection, heat destinations, reaction schedule or energy accounting.
use super::*;
use crate::cavity_phs::piston::thermal::ConductanceSample;
use crate::cavity_phs::piston::thermal::boundary::ports::network::{Contact, NetworkStep};
use fs_material::thermochemical_gas::ResolvedThermochemicalGasState;
use fs_thermal::enthalpy::CellId;
use std::collections::{BTreeMap, BTreeSet};

/// One explicitly exposed patch, independent of its heat-transfer law.
/// Areas of ALL gas patches on a cell are admitted together, across chambers.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct GasSurface {
    pub contact: Contact,
    pub area_m2: f64,
}

impl PreparedMaterialExchangeNetwork<'_, '_> {
    /// Compose a declared gas heat law with existing electrical losses, ordered
    /// thermal contacts, radiation, finite enthalpy cells and wall motion.
    ///
    /// `resolve` receives the actual gas energy midpoint and the complete cell
    /// material sample at the same nonlinear trial. Return TOTAL conductance
    /// in W/K and the original source record. The owner computes signed heat,
    /// sums it once per gas and debits each receiving cell separately. Do not
    /// deposit heat, advance chemistry, or mutate physical state in this pure
    /// callback. A zero-area surface must return zero conductance; zero area
    /// still evaluates its material law and does not hide invalid source data.
    ///
    /// Surface area must agree with the geometry used by the law. Aggregate
    /// area admission is not a spatial non-overlap proof. Gas conduction and
    /// convection may coexist on DISJOINT patches. Radiation may share them.
    /// Nonzero authored ambient convection on a bound cell refuses.
    ///
    /// This is the original uniform-cell/gas, independent-wall, first-order
    /// coupling, not a resolved boundary layer or a guarantee of timestep
    /// accuracy. Every accepted endpoint is returned together, or none is.
    #[allow(clippy::too_many_arguments, clippy::type_complexity)]
    pub fn step_with_gas_law_candidate<'m, I, E>(
        &self, cx: &Cx<'_>, thermal: &DepositionSimulation<'m>, target_s: f64,
        circuit: &MaterialCircuit, ports: &[f64], chambers: &[Chamber<'_>],
        surfaces: &[GasSurface], controls: NetworkSolve, guess_j: &[f64],
        resolve: impl Fn(Contact, &DepositionSimulation<'m>, &MaterialExchangeSample,
            &ResolvedThermochemicalGasState) -> Result<ConductanceSample<I>, E>,
    ) -> Result<NetworkCandidate<'m, ContactThermalSample<MaterialExchangeSample>, I>,
        CircuitEnthalpyError<NetworkError<MaterialExchangeError, E>>> {
        let admitted = GasSurfaceSet::admit(cx, thermal, chambers, surfaces, controls)?;
        let binding = self.binding.electrical;
        let routes = binding.heat_destinations(thermal).map_err(|e| coupled_error(e.into()))?;
        let destinations = binding.bulk.resistors.iter().map(|r| r.destination).collect::<Vec<_>>();
        let result = super::super::super::step_with_contacts_candidate(thermal, target_s,
            CircuitEnthalpyInput { circuit, ports, destinations: &destinations }, &routes,
            chambers, &admitted.contacts, controls, guess_j,
            |state| {
                // Recheck the original owner's exposure at every actual trial.
                admitted.check_trial(cx, state)?;
                self.resolve(cx, state)
            },
            |contact, state, source, gas| resolve(contact, state, source, gas),
            || cx.checkpoint().is_err());
        let mut candidate = match result {
            Err(CircuitEnthalpyError::Source(NetworkError::Material(e))) => return Err(coupled_error(e)),
            other => other?,
        };
        admitted.check_samples(candidate.chambers().iter()
            .flat_map(|chamber| chamber.boundary.source())
            .map(|sample| (sample.port.id, sample.conductance_w_k)))?;
        let radiation_j = candidate.material().material.radiation.as_ref()
            .map_or(0.0, |r| (target_s - thermal.time_s()) * r.net_heat_w());
        close_radiation(cx, &mut candidate.step, radiation_j, controls)?;
        Ok(candidate)
    }

    /// Reuse the original series gas/body conduction law in a mixed-law
    /// assembly. The COMPLETE endpoint material bundle remains authoritative,
    /// including any electrical properties and pinned source selections.
    /// The surrounding assembly must include this patch in aggregate coverage.
    pub fn resolve_conductive_gas_contact(
        &self, patch: MaterialContact, thermal: &DepositionSimulation<'_>,
        source: &MaterialExchangeSample, gas: &ResolvedThermochemicalGasState,
    ) -> Result<ConductanceSample<ResolvedMaterialStatePoint>, MaterialNetworkError> {
        resolve_contact(self.binding.electrical.bulk, patch.contact.cell, thermal,
            &source.bulk, gas, patch.geometry)
    }
}

// One surface admission and radiation correction for midpoint and reacting
// endpoint drivers. This owns no gas, thermal or electrical state.
pub(super) struct GasSurfaceSet {
    mapped: BTreeMap<u64, GasSurface>,
    areas: BTreeMap<CellId, f64>,
    pub(super) contacts: Vec<Contact>,
}
impl GasSurfaceSet {
    pub(super) fn admit<B>(
        cx: &Cx<'_>, thermal: &DepositionSimulation<'_>, chambers: &[Chamber<'_>],
        surfaces: &[GasSurface], controls: NetworkSolve,
    ) -> Result<Self, CircuitEnthalpyError<DrivenNetworkError<MaterialExchangeError, B>>> {
        checkpoint(cx).map_err(coupled_error)?;
        if surfaces.len() > controls.maximum_contacts || chambers.len() > controls.maximum_chambers {
            return Err(CircuitEnthalpyError::Invalid("gas surface/chamber capacity exceeded"));
        }
        let ids = chambers.iter().map(|chamber| chamber.id).collect::<BTreeSet<_>>();
        let mut mapped = BTreeMap::new();
        for &surface in surfaces {
            checkpoint(cx).map_err(coupled_error)?;
            if !surface.area_m2.is_finite() || surface.area_m2 < 0.0
                || !ids.contains(&surface.contact.chamber)
                || mapped.insert(surface.contact.id, surface).is_some() {
                return Err(CircuitEnthalpyError::Invalid("unique gas surfaces, existing chambers and finite nonnegative areas required"));
            }
        }
        let mut areas = BTreeMap::<CellId, f64>::new();
        for surface in mapped.values() {
            let total = areas.entry(surface.contact.cell).or_default();
            *total += surface.area_m2;
            if !total.is_finite() {
                return Err(CircuitEnthalpyError::Invalid("aggregate gas patch area is unrepresentable"));
            }
        }
        let contacts = mapped.values().map(|surface| surface.contact).collect();
        let admitted = Self { mapped, areas, contacts };
        admitted.check_trial(cx, thermal).map_err(coupled_error)?;
        Ok(admitted)
    }

    pub(super) fn check_trial(&self, cx: &Cx<'_>, thermal: &DepositionSimulation<'_>)
        -> Result<(), MaterialExchangeError> {
        for (&cell, &area) in &self.areas {
            checkpoint(cx)?;
            admit_area(thermal, cell, area)?;
        }
        Ok(())
    }

    pub(super) fn check_samples<B>(&self, samples: impl Iterator<Item = (u64, f64)>)
        -> Result<(), CircuitEnthalpyError<DrivenNetworkError<MaterialExchangeError, B>>> {
        for (id, conductance) in samples {
            let surface = self.mapped.get(&id).ok_or(CircuitEnthalpyError::Invalid("gas result names an absent surface"))?;
            if surface.area_m2 == 0.0 && conductance != 0.0 {
                return Err(CircuitEnthalpyError::Invalid("zero-area gas surface returned nonzero conductance"));
            }
        }
        Ok(())
    }
}

pub(super) fn close_radiation<B>(
    cx: &Cx<'_>, step: &mut NetworkStep, radiation_j: f64, controls: NetworkSolve,
) -> Result<(), CircuitEnthalpyError<DrivenNetworkError<MaterialExchangeError, B>>> {
    let external = step.external_network_heat_j-radiation_j;
    let residual = step.total_energy_residual_j+radiation_j;
    if !radiation_j.is_finite() || !external.is_finite() || !residual.is_finite()
        || residual.abs() > controls.total_energy_tolerance_j {
        return Err(coupled_error(MaterialExchangeError::Invalid("closed radiation cannot supply combined physical energy")));
    }
    checkpoint(cx).map_err(coupled_error)?;
    step.external_network_heat_j = external;
    step.total_energy_residual_j = residual;
    Ok(())
}

fn admit_area(thermal: &DepositionSimulation<'_>, id: CellId, area: f64)
    -> Result<(), MaterialExchangeError> {
    let invalid = MaterialExchangeError::Invalid;
    let cell = thermal.cells().get(&id).ok_or(invalid("gas surface names an inactive cell"))?;
    let exposed = thermal.exposed_area_m2(id).map_err(|e|
        MaterialExchangeError::Electrical(ElectricalContactNetworkError::Thermal(e)))?;
    if cell.boundary.htc_w_m2_k != 0.0 || area > exposed {
        return Err(invalid("combined gas patches exceed exposed surface or duplicate ambient convection"));
    }
    Ok(())
}
