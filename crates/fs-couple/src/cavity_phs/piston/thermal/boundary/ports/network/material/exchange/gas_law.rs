//! Gas-side heat laws share the original circuit/contact/radiation assembly.
//! This boundary changes only the constitutive response, not the integrator,
//! source selection, heat destinations, reaction schedule or energy accounting.
use super::*;
use crate::cavity_phs::piston::thermal::ConductanceSample;
use crate::cavity_phs::piston::thermal::boundary::ports::network::Contact;
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
        checkpoint(cx).map_err(coupled_error)?;
        if surfaces.len() > controls.maximum_contacts || chambers.len() > controls.maximum_chambers {
            return Err(CircuitEnthalpyError::Invalid("gas surface/chamber capacity exceeded"));
        }
        let ids = chambers.iter().map(|c| c.id).collect::<BTreeSet<_>>();
        let mut mapped = BTreeMap::new();
        for surface in surfaces {
            checkpoint(cx).map_err(coupled_error)?;
            if !surface.area_m2.is_finite() || surface.area_m2 < 0.0
                || !ids.contains(&surface.contact.chamber)
                || mapped.insert(surface.contact.id, surface).is_some() {
                return Err(CircuitEnthalpyError::Invalid("unique gas surfaces, existing chambers and finite nonnegative areas required"));
            }
        }
        let mut areas = BTreeMap::<CellId, f64>::new();
        // Stable physical contact-ID ordering, even for area accumulation.
        for surface in mapped.values() {
            let total = areas.entry(surface.contact.cell).or_default();
            *total += surface.area_m2;
            if !total.is_finite() {
                return Err(CircuitEnthalpyError::Invalid("aggregate gas patch area is unrepresentable"));
            }
        }
        for (&cell, &area) in &areas {
            checkpoint(cx).map_err(coupled_error)?;
            admit_area(thermal, cell, area).map_err(coupled_error)?;
        }
        let binding = self.binding.electrical;
        let routes = binding.heat_destinations(thermal).map_err(|e| coupled_error(e.into()))?;
        let destinations = binding.bulk.resistors.iter().map(|r| r.destination).collect::<Vec<_>>();
        let contacts = mapped.values().map(|s| s.contact).collect::<Vec<_>>();
        let result = super::super::super::step_with_contacts_candidate(thermal, target_s,
            CircuitEnthalpyInput { circuit, ports, destinations: &destinations }, &routes,
            chambers, &contacts, controls, guess_j,
            |state| {
                // Recheck the original owner's exposure at every actual trial.
                for (&cell, &area) in &areas {
                    checkpoint(cx)?;
                    admit_area(state, cell, area)?;
                }
                self.resolve(cx, state)
            },
            |contact, state, source, gas| resolve(contact, state, source, gas),
            || cx.checkpoint().is_err());
        let mut candidate = match result {
            Err(CircuitEnthalpyError::Source(NetworkError::Material(e))) => return Err(coupled_error(e)),
            other => other?,
        };
        // An admitted zero-area law cannot inject a hidden heat path. This is
        // checked before publication without replacing the callback's payload.
        for chamber in candidate.chambers() {
            for sample in chamber.boundary.source() {
                if mapped[&sample.port.id].area_m2 == 0.0 && sample.conductance_w_k != 0.0 {
                    return Err(CircuitEnthalpyError::Invalid("zero-area gas surface returned nonzero conductance"));
                }
            }
        }
        let radiation_j = candidate.material().material.radiation.as_ref()
            .map_or(0.0, |r| (target_s - thermal.time_s()) * r.net_heat_w());
        let external = candidate.step().external_network_heat_j - radiation_j;
        let residual = candidate.step().total_energy_residual_j + radiation_j;
        if !radiation_j.is_finite() || !external.is_finite() || !residual.is_finite()
            || residual.abs() > controls.total_energy_tolerance_j {
            return Err(coupled_error(MaterialExchangeError::Invalid("closed radiation cannot supply combined physical energy")));
        }
        checkpoint(cx).map_err(coupled_error)?;
        candidate.step.external_network_heat_j = external;
        candidate.step.total_energy_residual_j = residual;
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
