//! Conservative absorption and emission between native opaque coating bodies.
//!
//! Each link is a distinct pair of parallel, coextensive patches with unit
//! view factor across a phase-averaged far-field vacuum gap. Each side is
//! internally isothermal; different sides retain independent native T/h laws.
//! The original optical response owns polarized multiple reflections. Native
//! `HeatTransfer` applies identical signed watts to both caloric stores, without
//! external heat, a conductive contact, or a second thermal integration.

use super::*;
use crate::thermal_surface::{CoatingExchangeRates, CoatingGapResponse};
use fs_thermal::enthalpy::HeatTransfer;

/// One opaque specimen side, resolved from its original native caloric state.
/// `stack` must be `IsothermalCoatingStack::Opaque`; all its films and substrate
/// share this receiver's temperature, specific enthalpy and reference mass.
#[derive(Clone, Debug)]
pub struct CoatingGapSurface {
    pub source_identity: ContentHash,
    pub receiver: SurfaceHeatingReceiver,
    pub stack: IsothermalCoatingStack,
}

/// One physical gap. Signed transfer is from `surfaces[0]` to `surfaces[1]`.
/// Its positive physical area occupies that much currently exposed area on
/// EACH side. Distinct links must describe disjoint patches, even when they
/// share a receiver. Geometry/non-overlap and full view factor are declarations;
/// no distance-based view factor, occlusion, or diffraction is inferred.
#[derive(Clone, Debug)]
pub struct CoatingGapLink {
    pub link_id: u64,
    pub surfaces: [CoatingGapSurface; 2],
    pub area_m2: f64,
    pub quadrature: CoatingExchangeQuadrature,
}

#[derive(Clone, Copy, Debug)]
pub struct NativeCoatingGapBudget {
    pub max_links: usize,
    /// Distinct original optical-law rows, shared across both sides/all links.
    pub max_total_optical_knots: usize,
    /// Link rate resolutions across one advance, including endpoint checks.
    pub max_source_evaluations: usize,
    /// Original angle/wavelength solves; preparation charges two per node.
    pub max_optical_evaluations: usize,
}

/// Accepted internal radiation, already applied once by the native heat owner.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct NativeCoatingGapEnergy {
    pub duration_s: f64,
    pub by_link_first_to_second_j: BTreeMap<u64, f64>,
    pub by_cell_net_j: BTreeMap<CellId, f64>,
    /// Sum of signed per-cell internal radiation, an arithmetic closure witness.
    pub balance_residual_j: f64,
}

#[derive(Debug)]
struct GapModel {
    links: Vec<CoatingGapLink>,
    cells: BTreeMap<CellId, CellBinding>,
    budget: NativeCoatingGapBudget,
}
#[derive(Clone, Debug)]
struct PreparedGap {
    states: [[u64; 2]; 2],
    emitters: [CoatingThermalEmission; 2],
    response: Arc<CoatingGapResponse>,
    rates: CoatingExchangeRates,
}
#[derive(Clone, Debug, Default)]
struct GapEnergySums {
    duration: CheckedSum,
    links: Vec<CheckedSum>,
    cells: BTreeMap<CellId, CheckedSum>,
}

/// A network of independent physical vacuum gaps on ONE native thermal clock.
///
/// Use `DepositionSimulation::advance_observed`, optionally composed with
/// illumination, conduction and other existing observers. Each endpoint uses
/// its own native caloric law and temperature/enthalpy-dependent optics. The
/// prepared response is shared across clones and reused while both optical
/// identities remain unchanged. A failed advance publishes neither candidate.
/// This owner must be the sole radiative owner of its patches/common bands.
#[derive(Clone, Debug)]
pub struct NativeCoatingGapExchange {
    model: Arc<GapModel>,
    identity: ContentHash,
    time_s: f64,
    prepared_at_s: f64,
    gaps: Vec<PreparedGap>,
    transfers: Vec<HeatTransfer>,
    conductances_w_k: BTreeMap<CellId, f64>,
    work: NativeCoatingExchangeWork,
    sums: GapEnergySums,
    energy: NativeCoatingGapEnergy,
}

impl NativeCoatingGapExchange {
    pub fn try_new(
        thermal: &DepositionSimulation<'_>,
        links: &[CoatingGapLink],
        budget: NativeCoatingGapBudget,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<Self, NativeCoatingExchangeError> {
        checkpoint()?;
        if links.is_empty()
            || links.len() > budget.max_links
            || links.len() > 4096
            || budget.max_source_evaluations == 0
        {
            return Err(NativeCoatingExchangeError::Budget);
        }
        for link in links {
            for side in &link.surfaces {
                if side.stack.films().len() > MAX_COHERENT_STACK_LAYERS {
                    return Err(NativeCoatingExchangeError::Budget);
                }
            }
        }
        let mut links = links.to_vec();
        links.sort_by_key(|link| link.link_id);
        let mut cells = BTreeMap::new();
        let mut laws = BTreeSet::new();
        let mut rows = 0;
        for (index, link) in links.iter().enumerate() {
            checkpoint()?;
            if index > 0 && links[index - 1].link_id == link.link_id {
                return Err(invalid(
                    "each physical coating gap requires a unique link ID",
                ));
            }
            if !positive(link.area_m2)
                || link.surfaces[0].receiver.thermal_cell == link.surfaces[1].receiver.thermal_cell
            {
                return Err(invalid(
                    "a coating gap needs positive area and distinct native receivers",
                ));
            }
            for side in &link.surfaces {
                if side.source_identity == ContentHash([0; 32]) {
                    return Err(invalid(
                        "each coating gap side requires a nonzero source identity",
                    ));
                }
                if !matches!(side.stack, IsothermalCoatingStack::Opaque { .. }) {
                    return Err(CoatingExchangeError::OpaqueRequired.into());
                }
                let id = side.receiver.thermal_cell;
                if let std::collections::btree_map::Entry::Vacant(entry) = cells.entry(id) {
                    entry.insert(CellBinding::capture(thermal, id)?);
                }
                admit_law_rows(
                    &side.stack,
                    &mut laws,
                    &mut rows,
                    budget.max_total_optical_knots,
                )?;
            }
        }
        let mut owner = Self {
            sums: GapEnergySums {
                links: vec![CheckedSum::default(); links.len()],
                ..GapEnergySums::default()
            },
            model: Arc::new(GapModel {
                links,
                cells,
                budget,
            }),
            identity: ContentHash([0; 32]),
            time_s: thermal.time_s(),
            prepared_at_s: thermal.time_s(),
            gaps: Vec::new(),
            transfers: Vec::new(),
            conductances_w_k: BTreeMap::new(),
            work: NativeCoatingExchangeWork::default(),
            energy: NativeCoatingGapEnergy::default(),
        };
        owner.refresh(thermal, &mut checkpoint)?;
        let mut hash = DomainHasher::new("org.frankensim.render.native-coating-gap-exchange.v1");
        hash.update(&thermal.time_s().to_bits().to_le_bytes());
        for (&id, binding) in &owner.model.cells {
            binding.hash(id, &mut hash);
        }
        for (link, gap) in owner.model.links.iter().zip(&owner.gaps) {
            hash.update(&link.link_id.to_le_bytes());
            hash.update(gap.response.identity().as_bytes());
            for (side, emitter) in link.surfaces.iter().zip(&gap.emitters) {
                hash.update(&side.receiver.thermal_cell.0.to_le_bytes());
                hash.update(&side.receiver.material_id.0.to_le_bytes());
                for t in [
                    side.receiver.temperature_bounds_k.0,
                    side.receiver.temperature_bounds_k.1,
                ] {
                    hash.update(&t.to_bits().to_le_bytes());
                }
                hash.update(emitter.identity().as_bytes());
                for film in side.stack.films() {
                    if let CoatingOpticsModel::Law(law) = film {
                        hash.update(law.identity().as_bytes());
                    }
                }
                if let IsothermalCoatingStack::Opaque {
                    substrate: SubstrateOpticsModel::Law(law),
                    ..
                } = &side.stack
                {
                    hash.update(law.identity().as_bytes());
                }
            }
        }
        owner.identity = hash.finalize();
        checkpoint()?;
        Ok(owner)
    }

    pub const fn identity(&self) -> ContentHash {
        self.identity
    }
    pub const fn time_s(&self) -> f64 {
        self.time_s
    }
    pub const fn work(&self) -> NativeCoatingExchangeWork {
        self.work
    }
    pub fn heat_transfers(&self) -> &[HeatTransfer] {
        &self.transfers
    }
    pub fn conductance_bounds_w_k(&self) -> &BTreeMap<CellId, f64> {
        &self.conductances_w_k
    }
    pub fn energy(&self) -> &NativeCoatingGapEnergy {
        &self.energy
    }
    /// Latest resolved side, in the link's authored first/second order.
    pub fn emitter(&self, link_id: u64, side: usize) -> Option<&CoatingThermalEmission> {
        let index = self
            .model
            .links
            .binary_search_by_key(&link_id, |link| link.link_id)
            .ok()?;
        self.gaps[index].emitters.get(side)
    }

    fn validate(
        &self,
        thermal: &DepositionSimulation<'_>,
        check: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<BTreeMap<CellId, f64>, NativeCoatingExchangeError> {
        let mut temperatures = BTreeMap::new();
        for (&id, binding) in &self.model.cells {
            check()?;
            temperatures.insert(id, binding.temperature(thermal, id)?);
        }
        let selected: Vec<_> = self.model.cells.keys().copied().collect();
        let available = thermal.exposed_areas_m2(&selected, &mut *check)?;
        let mut occupied = BTreeMap::<CellId, CheckedSum>::new();
        for link in &self.model.links {
            for side in &link.surfaces {
                check()?;
                let receiver = side.receiver;
                let (lo, hi) = receiver.temperature_bounds_k;
                if !positive(lo)
                    || !positive(hi)
                    || hi < lo
                    || self.model.cells[&receiver.thermal_cell].material != receiver.material_id
                {
                    return Err(invalid(
                        "gap receiver requires its original material and positive ordered temperature support",
                    ));
                }
                let t = temperatures[&receiver.thermal_cell];
                if !positive(t) || t < lo || t > hi {
                    return Err(NativeCoatingEmissionError::Binding(
                        SurfaceHeatingError::TemperatureDomain {
                            thermal_cell: receiver.thermal_cell,
                            temperature_k: t,
                            minimum_k: lo,
                            maximum_k: hi,
                        },
                    )
                    .into());
                }
                occupied
                    .entry(receiver.thermal_cell)
                    .or_default()
                    .add(link.area_m2)?;
            }
        }
        for (cell, area) in occupied {
            if area.value()? > available[&cell] {
                return Err(invalid("coating gaps exceed current native exposed area"));
            }
        }
        Ok(temperatures)
    }

    fn refresh(
        &mut self,
        thermal: &DepositionSimulation<'_>,
        check: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), NativeCoatingExchangeError> {
        let temperatures = self.validate(thermal, check)?;
        let mut gaps = Vec::with_capacity(self.model.links.len());
        let mut transfers = Vec::with_capacity(self.model.links.len());
        let mut conductances = BTreeMap::<CellId, CheckedSum>::new();
        for (index, link) in self.model.links.iter().enumerate() {
            check()?;
            self.work.source_evaluations = self
                .work
                .source_evaluations
                .checked_add(1)
                .ok_or(NativeCoatingExchangeError::Budget)?;
            if self.work.source_evaluations > self.model.budget.max_source_evaluations {
                return Err(NativeCoatingExchangeError::Budget);
            }
            let ids = link
                .surfaces
                .each_ref()
                .map(|side| side.receiver.thermal_cell);
            let states = ids.map(|id| {
                [
                    temperatures[&id].to_bits(),
                    thermal.cells()[&id].specific_h_j_kg.to_bits(),
                ]
            });
            let old = self.gaps.get(index);
            let mut resolve = |side: usize| {
                if let Some(old) = old.filter(|old| old.states[side] == states[side]) {
                    Ok(old.emitters[side].clone())
                } else {
                    let s = &link.surfaces[side];
                    resolve_stack_emitter(
                        &s.stack,
                        s.receiver,
                        s.source_identity,
                        thermal,
                        self.model.budget.max_total_optical_knots,
                        check,
                    )
                }
            };
            let emitters = [resolve(0)?, resolve(1)?];
            let optical_ids = emitters
                .each_ref()
                .map(CoatingThermalEmission::optical_identity);
            let response = if let Some(old) =
                old.filter(|old| old.response.optical_identities() == optical_ids)
            {
                Arc::clone(&old.response)
            } else {
                let evaluations = link
                    .quadrature
                    .wavelength_bins
                    .checked_mul(link.quadrature.cosine_bins)
                    .and_then(|n| n.checked_mul(2))
                    .ok_or(NativeCoatingExchangeError::Budget)?;
                self.work.optical_evaluations = self
                    .work
                    .optical_evaluations
                    .checked_add(evaluations)
                    .ok_or(NativeCoatingExchangeError::Budget)?;
                if self.work.optical_evaluations > self.model.budget.max_optical_evaluations {
                    return Err(NativeCoatingExchangeError::Budget);
                }
                Arc::new(optical(check, |checkpoint| {
                    CoatingGapResponse::from_emitters(
                        &emitters[0],
                        &emitters[1],
                        link.area_m2,
                        link.quadrature,
                        CoatingExchangeBudget {
                            max_spectral_nodes: link.quadrature.wavelength_bins,
                            max_optical_evaluations: evaluations,
                        },
                        checkpoint,
                    )
                })?)
            };
            let rates = if let Some(old) = old.filter(|old| old.states == states) {
                old.rates
            } else {
                optical(check, |checkpoint| {
                    response.evaluate(temperatures[&ids[0]], temperatures[&ids[1]], checkpoint)
                })?
            };
            transfers.push(HeatTransfer {
                a: ids[0],
                b: ids[1],
                power_w: rates.outward_power_w,
            });
            for id in ids {
                conductances.entry(id).or_default().add(rates.secant_w_k)?;
            }
            gaps.push(PreparedGap {
                states,
                emitters,
                response,
                rates,
            });
        }
        self.conductances_w_k = conductances
            .into_iter()
            .map(|(id, sum)| Ok((id, sum.value()?)))
            .collect::<Result<_, OpticalHeatingError>>()?;
        self.gaps = gaps;
        self.transfers = transfers;
        self.prepared_at_s = thermal.time_s();
        check()?;
        Ok(())
    }
}

impl DepositionObserver for NativeCoatingGapExchange {
    type Error = NativeCoatingExchangeError;
    fn begin(&mut self) -> Result<(), Self::Error> {
        self.work = NativeCoatingExchangeWork::default();
        Ok(())
    }
    fn heat_transfers(
        &mut self,
        thermal: &DepositionSimulation<'_>,
        check: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<Vec<HeatTransfer>, Self::Error> {
        if thermal.time_s() != self.time_s {
            return Err(invalid("coating gap clock mismatch"));
        }
        self.refresh(thermal, check)?;
        Ok(self.transfers.clone())
    }
    fn limit_timestep_s(
        &mut self,
        thermal: &DepositionSimulation<'_>,
        proposed_s: f64,
        check: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<f64, Self::Error> {
        if thermal.time_s() != self.time_s || thermal.time_s() != self.prepared_at_s {
            return Err(invalid(
                "gap timestep requires its prepared native start state",
            ));
        }
        Ok(thermal.limit_timestep_for_additional_conductances(
            proposed_s,
            &self.conductances_w_k,
            check,
        )?)
    }
    fn event(&mut self, _: &Event, thermal: &DepositionSimulation<'_>) -> Result<(), Self::Error> {
        if thermal.time_s() != self.time_s {
            return Err(invalid("coating gap event clock mismatch"));
        }
        self.refresh(thermal, &mut || Ok(()))
    }
    fn step(
        &mut self,
        start_s: f64,
        duration_s: f64,
        _: &BTreeMap<CellId, f64>,
        thermal: &DepositionSimulation<'_>,
        check: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), Self::Error> {
        if start_s != self.time_s
            || self.prepared_at_s != start_s
            || !positive(duration_s)
            || thermal.time_s() <= start_s
        {
            return Err(invalid(
                "gap energy requires its matching accepted native interval",
            ));
        }
        let mut sums = self.sums.clone();
        sums.duration.add(duration_s)?;
        for (sum, transfer) in sums.links.iter_mut().zip(&self.transfers) {
            check()?;
            let energy = finite(transfer.power_w * duration_s)?;
            sum.add(energy)?;
            sums.cells.entry(transfer.a).or_default().add(-energy)?;
            sums.cells.entry(transfer.b).or_default().add(energy)?;
        }
        // Account the actual accepted-start pair watts before endpoint optics.
        self.refresh(thermal, check)?;
        let by_cell_net_j = sums
            .cells
            .iter()
            .map(|(&id, sum)| Ok((id, sum.value()?)))
            .collect::<Result<BTreeMap<_, _>, OpticalHeatingError>>()?;
        let mut residual = CheckedSum::default();
        for &energy in by_cell_net_j.values() {
            residual.add(energy)?;
        }
        let energy = NativeCoatingGapEnergy {
            duration_s: sums.duration.value()?,
            by_cell_net_j,
            balance_residual_j: residual.value()?,
            by_link_first_to_second_j: self
                .model
                .links
                .iter()
                .zip(&sums.links)
                .map(|(link, sum)| Ok((link.link_id, sum.value()?)))
                .collect::<Result<_, OpticalHeatingError>>()?,
        };
        check()?;
        self.sums = sums;
        self.energy = energy;
        self.time_s = thermal.time_s();
        Ok(())
    }
}
