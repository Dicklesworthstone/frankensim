//! Finite adsorbent storage on the original heterogeneous pore timeline.
//!
//! Every accepted interval samples the original dusty-gas throats and original
//! competitive site laws on one gas snapshot. Their existing joint owner spends
//! gross gas inventories once, updates every finite solid temperature, and
//! reconstructs each gas with its original species/energy law. The next interval
//! sees those accepted inventories and temperatures, including adsorption heat.
//!
//! Geometry stays rigid: bound molecules are site occupancy, not a prescribed
//! liquid saturation. No sorption-to-permeability law, capillary work, moving
//! pore volume or separate Darcy heating is inferred. Explicit skeleton heat
//! paths retain the original dry caloric, conductivity and interface sources.

/// Componentwise temporal accuracy for the complete original gas/site/heat system.
pub mod adaptive;

use super::*;
use crate::adsorption::{
    AdsorptionError,
    competitive::{
        CompetitiveSolid, CompetitiveSolve,
        network::{self, AdsorptionNetworkBinding, AdsorptionNetworkError, AdsorptionNetworkStep,
            conduction::AdsorbentConductionLink, gas_heat::AdsorbentGasHeatLink},
    },
};

/// One independently stored competitive site population per selected pore body.
/// Cell indices are zero-based and cannot bind the two external reservoirs.
/// Dry mass and source charts belong to each original solid; geometric porosity
/// does not silently supply a skeleton density or a capacity calibration.
#[derive(Clone, Debug, PartialEq)]
pub struct AdsorbingPoreDomain {
    domain: PoreDomain,
    solids: BTreeMap<usize, CompetitiveSolid>,
    catalog: ThermochemicalGasModel,
}

/// Gross surface streams survive adsorption/desorption reversals in the window.
#[derive(Clone, Debug, PartialEq)]
pub struct PoreAdsorptionTransfer {
    pub cell: usize,
    pub species: SpeciesId,
    pub gross_adsorption_mol: f64,
    pub gross_desorption_mol: f64,
}

/// Gross internal heat streams survive direction changes across the window.
/// Original source identities and the final actual layer/interface sample are
/// retained by `AdsorbingPoreStep::pore.last_substep.conduction`.
#[derive(Clone, Debug, PartialEq)]
pub struct PoreConductionTransfer {
    pub id: u64,
    pub cells: [usize; 2],
    pub gross_a_to_b_j: f64,
    pub gross_b_to_a_j: f64,
}

/// Gross sensible heat on a sourced gas/solid path over the complete window.
/// Both directions survive reversals; the last native sample retains the
/// original source queries. This heat is separate from adsorbed-species enthalpy.
#[derive(Clone, Debug, PartialEq)]
pub struct PoreGasHeatTransfer {
    pub id: u64,
    pub cell: usize,
    pub gross_gas_to_solid_j: f64,
    pub gross_solid_to_gas_j: f64,
}

/// Bounded window accounting, retaining the actual final native source report.
#[derive(Clone, Debug, PartialEq)]
pub struct AdsorbingPoreStep {
    /// Whole-window gas+bound-species mass and gas+solid energy/entropy admission.
    /// `pore_mass_change_kg` remains the change of free gas in pore bodies.
    pub pore: PoreDomainStep<AdsorptionNetworkStep>,
    pub adsorption: Vec<PoreAdsorptionTransfer>,
    pub conduction: Vec<PoreConductionTransfer>,
    pub gas_heat: Vec<PoreGasHeatTransfer>,
    pub adsorbed_mass_change_kg: f64,
    /// Sum of the original dry+bound sensible and binding-energy increments.
    pub solid_energy_change_j: f64,
}

#[derive(Clone, Debug, PartialEq)]
pub struct AdsorbingPoreCandidate {
    domain: AdsorbingPoreDomain,
    step: AdsorbingPoreStep,
}
impl AdsorbingPoreCandidate {
    pub fn domain(&self) -> &AdsorbingPoreDomain {
        &self.domain
    }
    pub fn step(&self) -> &AdsorbingPoreStep {
        &self.step
    }
    pub fn into_parts(self) -> (AdsorbingPoreDomain, AdsorbingPoreStep) {
        (self.domain, self.step)
    }
}

impl PoreSubstep for AdsorptionNetworkStep {
    fn exchange(&self) -> &SpeciesNetworkStep {
        &self.exchange
    }
}

impl AdsorbingPoreDomain {
    /// Attach retained finite site populations to physical pore bodies.
    /// Every actual update also checks the original complete thermochemical
    /// catalog and each source/domain through the native adsorption owner.
    pub fn try_new(
        domain: PoreDomain,
        solids: BTreeMap<usize, CompetitiveSolid>,
        catalog: ThermochemicalGasModel,
    ) -> Result<Self, PoreDomainError> {
        if solids.is_empty() || solids.keys().any(|&cell| cell >= domain.cells.len()) {
            return Err(Error::Invalid(
                "finite adsorbents must bind existing pore bodies exactly once",
            ));
        }
        Ok(Self {
            domain,
            solids,
            catalog,
        })
    }
    pub fn pore_domain(&self) -> &PoreDomain {
        &self.domain
    }
    pub fn solids(&self) -> &BTreeMap<usize, CompetitiveSolid> {
        &self.solids
    }
    pub fn catalog(&self) -> &ThermochemicalGasModel {
        &self.catalog
    }

    /// Evolve transport, competitive loading and finite thermal storage together.
    /// The SAME native bounded pore loop handles stepping, typed retries, work
    /// limits, exact face throughput and cancellation. No accepted input changes
    /// until the complete gas+solid window passes its single balance allowance.
    ///
    /// Energy defects and complete-system entropy increments are compensated
    /// sums of the native reference-cancelled substeps. Actual final gas+bound
    /// species inventories are independently checked across the entire window.
    /// Per-step mass/energy/entropy limits are not multiplied by the step count.
    #[allow(clippy::too_many_arguments)]
    pub fn step_candidate(
        &self,
        duration_s: f64,
        faces: &[PoreFace<'_>],
        endpoints: &[GasEnergySolve],
        controls: PoreDomainSolve,
        adsorption_controls: &BTreeMap<usize, CompetitiveSolve>,
        cancelled: impl FnMut() -> bool,
    ) -> Result<AdsorbingPoreCandidate, PoreDomainError> {
        self.step_candidate_with_conduction(duration_s, faces, endpoints, controls,
            adsorption_controls, &[], cancelled)
    }

    /// Evolve explicitly authored source-backed skeleton heat paths on the same
    /// pore timeline. Link solid IDs are the zero-based pore-cell indices.
    /// Flux is frozen at each accepted starting state; later intervals observe
    /// the resulting temperatures through the original adsorption/flow owners.
    /// Every source, energy, entropy and cancellation refusal preserves the
    /// original gas, site, thermal and clock histories.
    #[allow(clippy::too_many_arguments)]
    pub fn step_candidate_with_conduction(
        &self,
        duration_s: f64,
        faces: &[PoreFace<'_>],
        endpoints: &[GasEnergySolve],
        controls: PoreDomainSolve,
        adsorption_controls: &BTreeMap<usize, CompetitiveSolve>,
        links: &[AdsorbentConductionLink<'_, '_>],
        cancelled: impl FnMut() -> bool,
    ) -> Result<AdsorbingPoreCandidate, PoreDomainError> {
        self.step_candidate_with_heat(duration_s, faces, endpoints, controls,
            adsorption_controls, links, &[], cancelled)
    }

    /// Evolve explicitly sourced gas/solid sensible heat on the original pore
    /// timeline, together with skeleton heat, adsorption and dusty-gas flow.
    /// Gas-link solid IDs are zero-based pore-cell indices. Every interval uses
    /// its actual accepted gas temperature/pressure and finite solid temperature;
    /// resulting temperatures feed the next original kinetics and transport solve.
    /// The complete window retains one absolute energy and entropy allowance.
    #[allow(clippy::too_many_arguments)]
    pub fn step_candidate_with_heat(
        &self,
        duration_s: f64,
        faces: &[PoreFace<'_>],
        endpoints: &[GasEnergySolve],
        controls: PoreDomainSolve,
        adsorption_controls: &BTreeMap<usize, CompetitiveSolve>,
        links: &[AdsorbentConductionLink<'_, '_>],
        gas_links: &[AdsorbentGasHeatLink<'_>],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<AdsorbingPoreCandidate, PoreDomainError> {
        check(&mut cancelled)?;
        if !self.solids.keys().eq(adsorption_controls.keys()) {
            return Err(Error::Invalid(
                "every retained pore adsorbent needs its exact local controls",
            ));
        }
        let entropy_limit = adsorption_controls
            .values()
            .fold(controls.transport.entropy_tolerance_j_k, |limit, solve| {
                limit.min(solve.entropy_tolerance_j_k)
            });
        let mut solids = self.solids.clone();
        let mut energy = Sum::default();
        let mut entropy = Sum::default();
        let mut solid_energy = Sum::default();
        let mut surface_totals: BTreeMap<(usize, SpeciesId), [Sum; 2]> = BTreeMap::new();
        let mut heat_totals: BTreeMap<u64, ([usize; 2], [Sum; 2])> = BTreeMap::new();
        let mut gas_heat_totals: BTreeMap<u64, (usize, [Sum; 2])> = BTreeMap::new();
        let run = self.domain.run_window(
            duration_s,
            faces,
            endpoints,
            controls,
            &mut cancelled,
            |parcels, paths, dt, endpoints, transport, cancelled| {
                let bindings: Vec<_> = solids
                    .iter()
                    .map(|(&cell, solid)| AdsorptionNetworkBinding {
                        id: cell as u64,
                        parcel: cell as u64 + 1,
                        solid,
                        solve: adsorption_controls[&cell],
                    })
                    .collect();
                let candidate = network::step_candidate_with_heat(
                    parcels,
                    paths,
                    &bindings,
                    &self.catalog,
                    dt,
                    endpoints,
                    transport,
                    links,
                    gas_links,
                    &mut *cancelled,
                )
                .map_err(PoreAdvanceError::Adsorption)?;
                drop(bindings);
                let (next, updated, step) = candidate.into_parts();
                add(&mut energy, step.total_energy_residual_j)?;
                add(&mut entropy, step.entropy_change_j_k)?;
                for row in &step.solids {
                    if cancelled() {
                        return Err(PoreAdvanceError::Adsorption(
                            AdsorptionNetworkError::Cancelled,
                        ));
                    }
                    add(&mut solid_energy, row.solid_energy_change_j)?;
                    for (i, species) in row.species.iter().enumerate() {
                        let totals = surface_totals
                            .entry((row.id as usize, species.species.clone()))
                            .or_insert_with(|| [Sum::default(), Sum::default()]);
                        add(&mut totals[0], row.sites.gross_adsorption_mol[i])?;
                        add(&mut totals[1], row.sites.gross_desorption_mol[i])?;
                    }
                }
                for sample in &step.conduction {
                    if cancelled() {
                        return Err(PoreAdvanceError::Adsorption(
                            AdsorptionNetworkError::Cancelled,
                        ));
                    }
                    let (_, totals) = heat_totals.entry(sample.id).or_insert_with(|| (
                        sample.solids.map(|id| id as usize),
                        [Sum::default(), Sum::default()],
                    ));
                    if sample.energy_a_to_b_j >= 0.0 {
                        add(&mut totals[0], sample.energy_a_to_b_j)?;
                    } else {
                        add(&mut totals[1], -sample.energy_a_to_b_j)?;
                    }
                }
                for sample in &step.gas_heat {
                    if cancelled() {
                        return Err(PoreAdvanceError::Adsorption(AdsorptionNetworkError::Cancelled));
                    }
                    let (_, totals) = gas_heat_totals.entry(sample.id).or_insert_with(|| (
                        sample.solid as usize, [Sum::default(), Sum::default()],
                    ));
                    if sample.energy_gas_to_solid_j >= 0.0 {
                        add(&mut totals[0], sample.energy_gas_to_solid_j)?;
                    } else {
                        add(&mut totals[1], -sample.energy_gas_to_solid_j)?;
                    }
                }
                solids = updated
                    .into_iter()
                    .map(|(id, solid)| (id as usize, solid))
                    .collect();
                Ok((next, step))
            },
        )?;
        let (mass, species, pore_mass, bound_mass) =
            inventories(self, &run.domain, &solids, &mut cancelled)?;
        if mass.abs() > controls.transport.exchange.mass_tolerance_kg
            || species > controls.transport.exchange.mass_tolerance_kg
            || energy.value.abs() > controls.transport.exchange.energy_tolerance_j
            || entropy.value < -entropy_limit
        {
            return Err(Error::Invalid(
                "complete adsorbing pore window exceeded its species/mass/energy/entropy budget",
            ));
        }
        let adsorption = surface_totals
            .into_iter()
            .map(|((cell, species), totals)| PoreAdsorptionTransfer {
                cell,
                species,
                gross_adsorption_mol: totals[0].value,
                gross_desorption_mol: totals[1].value,
            })
            .collect();
        let conduction = heat_totals.into_iter().map(|(id, (cells, totals))|
            PoreConductionTransfer {
                id, cells, gross_a_to_b_j: totals[0].value, gross_b_to_a_j: totals[1].value,
            }).collect();
        let gas_heat = gas_heat_totals.into_iter().map(|(id, (cell, totals))|
            PoreGasHeatTransfer { id, cell,
                gross_gas_to_solid_j: totals[0].value, gross_solid_to_gas_j: totals[1].value,
            }).collect();
        check(&mut cancelled)?;
        let (domain, pore) = run.finish((mass, species, energy.value, entropy.value, pore_mass));
        Ok(AdsorbingPoreCandidate {
            domain: Self {
                domain,
                solids,
                catalog: self.catalog.clone(),
            },
            step: AdsorbingPoreStep {
                pore,
                adsorption,
                conduction,
                gas_heat,
                adsorbed_mass_change_kg: bound_mass,
                solid_energy_change_j: solid_energy.value,
            },
        })
    }

    /// Publish the gas, every finite adsorbent and physical clock atomically.
    #[allow(clippy::too_many_arguments)]
    pub fn advance(
        &mut self,
        duration_s: f64,
        faces: &[PoreFace<'_>],
        endpoints: &[GasEnergySolve],
        controls: PoreDomainSolve,
        adsorption_controls: &BTreeMap<usize, CompetitiveSolve>,
        cancelled: impl FnMut() -> bool,
    ) -> Result<AdsorbingPoreStep, PoreDomainError> {
        self.advance_with_conduction(duration_s, faces, endpoints, controls,
            adsorption_controls, &[], cancelled)
    }

    /// Atomically publish the same native gas/site/thermal candidate with
    /// explicitly sourced internal heat paths and full-window heat receipts.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_with_conduction(
        &mut self,
        duration_s: f64,
        faces: &[PoreFace<'_>],
        endpoints: &[GasEnergySolve],
        controls: PoreDomainSolve,
        adsorption_controls: &BTreeMap<usize, CompetitiveSolve>,
        links: &[AdsorbentConductionLink<'_, '_>],
        cancelled: impl FnMut() -> bool,
    ) -> Result<AdsorbingPoreStep, PoreDomainError> {
        self.advance_with_heat(duration_s, faces, endpoints, controls,
            adsorption_controls, links, &[], cancelled)
    }

    /// Atomically publish the same native candidate including finite gas/solid
    /// sensible heat. A late source, work, balance or cancellation refusal leaves
    /// every original gas, site population, temperature and clock unchanged.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_with_heat(
        &mut self,
        duration_s: f64,
        faces: &[PoreFace<'_>],
        endpoints: &[GasEnergySolve],
        controls: PoreDomainSolve,
        adsorption_controls: &BTreeMap<usize, CompetitiveSolve>,
        links: &[AdsorbentConductionLink<'_, '_>],
        gas_links: &[AdsorbentGasHeatLink<'_>],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<AdsorbingPoreStep, PoreDomainError> {
        let candidate = self.step_candidate_with_heat(
            duration_s,
            faces,
            endpoints,
            controls,
            adsorption_controls,
            links,
            gas_links,
            &mut cancelled,
        )?;
        check(&mut cancelled)?;
        let (next, step) = candidate.into_parts();
        *self = next;
        Ok(step)
    }
}

fn add(sum: &mut Sum, value: f64) -> Result<(), PoreAdvanceError> {
    sum.add(value).map_err(|error| {
        PoreAdvanceError::Adsorption(AdsorptionNetworkError::Thermodynamic(error.into()))
    })
}

// Do not retry source or model-domain errors based on their display text.
pub(super) fn retryable(error: &AdsorptionNetworkError) -> bool {
    match error {
        AdsorptionNetworkError::Transport(error) => super::retryable(error),
        AdsorptionNetworkError::Gas(error) => {
            super::retryable(&GasTransportError::Exchange(error.clone()))
        }
        AdsorptionNetworkError::Adsorption { error, .. }
        | AdsorptionNetworkError::Thermodynamic(error) => matches!(
            error,
            AdsorptionError::InventoryLimit { .. }
                | AdsorptionError::ThermalRange
                | AdsorptionError::NotConverged { .. }
                | AdsorptionError::EntropyDecreased { .. }
                | AdsorptionError::Gas(
                    GasEnergyError::OutsideBracket { .. } | GasEnergyError::NotConverged { .. }
                )
        ),
        _ => false,
    }
}

fn inventories(
    old: &AdsorbingPoreDomain,
    next: &PoreDomain,
    solids: &BTreeMap<usize, CompetitiveSolid>,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(f64, f64, f64, f64), PoreDomainError> {
    if !old.solids.keys().eq(solids.keys()) || !old.domain.parcels.keys().eq(next.parcels.keys()) {
        return Err(Error::Invalid(
            "adsorbing pore update changed the retained body/site roster",
        ));
    }
    let mut total = Sum::default();
    let mut bound_mass = Sum::default();
    let mut maximum = 0.0_f64;
    for source in old.catalog.thermodynamics().components() {
        let mut species = Sum::default();
        for (&id, before) in &old.domain.parcels {
            check(cancelled)?;
            species
                .add(
                    next.parcels[&id]
                        .species_mass_kg(source.species())
                        .map_err(gas_error)?
                        - before
                            .species_mass_kg(source.species())
                            .map_err(gas_error)?,
                )
                .map_err(gas_error)?;
        }
        for (&cell, before) in &old.solids {
            check(cancelled)?;
            let after = &solids[&cell];
            if before.model() != after.model() || before.dry_mass_kg() != after.dry_mass_kg() {
                return Err(Error::Invalid(
                    "adsorbing pore update replaced an original solid source or dry mass",
                ));
            }
            let delta = (after.adsorbed_moles(source.species()).unwrap_or(0.0)
                - before.adsorbed_moles(source.species()).unwrap_or(0.0))
                * source.molar_mass().value();
            species.add(delta).map_err(gas_error)?;
            bound_mass.add(delta).map_err(gas_error)?;
        }
        maximum = maximum.max(species.value.abs());
        total.add(species.value).map_err(gas_error)?;
    }
    let mut pore_mass = Sum::default();
    for cell in 0..old.domain.cells.len() {
        check(cancelled)?;
        let id = cell as u64 + 1;
        pore_mass
            .add(next.parcels[&id].mass_kg() - old.domain.parcels[&id].mass_kg())
            .map_err(gas_error)?;
    }
    Ok((total.value, maximum, pore_mass.value, bound_mass.value))
}
