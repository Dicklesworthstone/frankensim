//! Several liquid surfaces sharing one original backward-Euler heat graph.
//!
//! One interface query per selected cell and two original heat-rate evaluations
//! resolve the complete endpoint, independent of the number of phase surfaces.
//! Each surface retains its own law, film, liquid domain and vapor receipt.
//! All proposed masses are visible to one source callback and installed only
//! after the shared heat balance and every local phase gate have passed.

use super::super::{LiquidBranch, ThermalEvaporationError};
use crate::enthalpy::backward_euler::{
    BackwardEulerStep, BackwardEulerTrial, BoundaryExchangeStepCoefficients,
    BoundaryStepCoefficients, StepCoefficients, TransportStepCoefficients,
};
use crate::enthalpy::observe::EnergySum;
use crate::enthalpy::{CellId, DepositionError, DepositionSimulation, finite, nonnegative};
use crate::evaporation::{
    EvaporationError, EvaporationLaw, GasFilm, LiquidSurfaceState, SurfaceStep, flux, properties,
};
use std::cell::RefCell;
use std::collections::BTreeMap;

/// Original liquid identity and caloric support for a source-resolved gas film.
/// The film is queried from the actual initial/endpoint state rather than stored
/// here. No pressure, humidity, gas temperature or transport coefficient is inferred.
pub struct EvaporatingCellSource<'l, L: EvaporationLaw + ?Sized = dyn EvaporationLaw + 'l> {
    pub id: CellId,
    pub law: &'l L,
    pub pure_liquid_enthalpy_range_j_kg: [f64; 2],
}
impl<L: EvaporationLaw + ?Sized> Copy for EvaporatingCellSource<'_, L> {}
impl<L: EvaporationLaw + ?Sized> Clone for EvaporatingCellSource<'_, L> {
    fn clone(&self) -> Self {
        *self
    }
}

/// One existing thermal cell's phase law and prescribed far-field boundary.
/// Mass, geometry, caloric identity and stored enthalpy stay in the timeline.
/// The default trait object admits heterogeneous laws in the same binding slice;
/// the generic form also retains unsized single-law callers without boxing.
pub struct EvaporatingCellBinding<'l, L: EvaporationLaw + ?Sized = dyn EvaporationLaw + 'l> {
    pub id: CellId,
    pub law: &'l L,
    pub film: GasFilm,
    pub pure_liquid_enthalpy_range_j_kg: [f64; 2],
}
impl<L: EvaporationLaw + ?Sized> Copy for EvaporatingCellBinding<'_, L> {}
impl<L: EvaporationLaw + ?Sized> Clone for EvaporatingCellBinding<'_, L> {
    fn clone(&self) -> Self {
        *self
    }
}

/// Physical open-system account on every selected cell's original caloric datum.
#[derive(Clone, Debug, PartialEq)]
pub struct EvaporatingNetworkBackwardEulerAudit {
    /// Original external sources and thermal baths, excluding internal links [J].
    pub external_heat_j: f64,
    /// Actual network storage change, including every changed liquid mass [J].
    pub stored_energy_change_j: f64,
    /// Network mass change plus all represented vapor transfers [kg].
    pub mass_balance_residual_kg: f64,
    /// Storage plus full carried vapor enthalpy minus external heat [J],
    /// evaluated after cancelling each carried-liquid reference datum.
    pub energy_residual_j: f64,
    /// Canonical CellId order. Each receipt uses its own film and latent heat;
    /// signed opposite transfers never replace the individual acceptance gates.
    pub surfaces: BTreeMap<CellId, SurfaceStep>,
}

/// Original thermal sources plus complete vapor-stream enthalpy supplied by
/// an external owner. Every selected cell must appear, including zero excess.
pub struct VaporEnthalpyStepCoefficients<T> {
    /// Physical heat sources only; exclude the vapor-carried energy below.
    pub heat: BoundaryExchangeStepCoefficients<T>,
    /// Complete signed-stream specific enthalpy minus `h_liquid + L` [J/kg],
    /// with all terms expressed in that cell's retained liquid datum. The
    /// value may have either sign. For condensation it includes any upstream
    /// vapor sensible-enthalpy difference; source compatibility remains with
    /// the gas/material owner rather than being inferred by this adapter.
    pub vapor_enthalpy_excess_j_kg: BTreeMap<CellId, f64>,
}

/// Complete vapor-enthalpy account, deliberately distinct from the existing
/// prescribed-film timeline account. Latent heat remains one subaccount of
/// each full vapor receipt; the correction is never added to latent heat.
#[derive(Clone, Debug, PartialEq)]
pub struct EvaporatingNetworkVaporEnthalpyAudit {
    /// Full corrected local/global energy residuals and original mass gates.
    /// Each surface's vapor_energy_j includes its correction exactly once.
    pub thermal: EvaporatingNetworkBackwardEulerAudit,
    /// Signed `transferred_mass * vapor_enthalpy_excess` [J] per selected cell.
    /// Retained directly before adding an arbitrary liquid caloric reference,
    /// rather than recovered by subtracting large vapor-enthalpy totals.
    pub vapor_enthalpy_correction_j: BTreeMap<CellId, f64>,
}

struct PreparedSurface<'l, L: EvaporationLaw + ?Sized> {
    binding: EvaporatingCellBinding<'l, L>,
    branch: LiquidBranch,
    area_m2: f64,
    initial_transfer_kg: f64,
}

/// One event-free interval over the original heat network and selected liquids.
/// Energy unknowns include ALL active heat cells in the native CellId order.
pub struct EvaporatingNetworkBackwardEulerStep<
    's,
    'm,
    'l,
    L: EvaporationLaw + ?Sized = dyn EvaporationLaw + 'l,
> {
    heat: BackwardEulerStep<'s, 'm>,
    initial: &'s DepositionSimulation<'m>,
    surfaces: BTreeMap<CellId, PreparedSurface<'l, L>>,
}

impl<'s, 'm, 'l, L: EvaporationLaw + ?Sized> EvaporatingNetworkBackwardEulerStep<'s, 'm, 'l, L> {
    /// Bind active, distinct liquid cells and freeze their current exposed area.
    /// Empty bindings explicitly select heat-only transport. Start-time events
    /// must already be processed; an interior event refuses and endpoint events
    /// remain pending. Unsupported or initially boiling surfaces refuse.
    pub fn new(
        initial: &'s DepositionSimulation<'m>,
        bindings: &[EvaporatingCellBinding<'l, L>],
        target_s: f64,
        checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<Self, ThermalEvaporationError> {
        Self::prepare(initial, bindings, target_s, None, checkpoint)
    }

    /// Prepare the same heat/phase owner with source-resolved INITIAL films.
    /// Every selected cell must be active and appear exactly once in the map.
    /// The returned step retains these films as its initial physical reference;
    /// use `evaluate_with_films` to resolve later speculative endpoints.
    pub fn new_with_films(
        initial: &'s DepositionSimulation<'m>,
        sources: &[EvaporatingCellSource<'l, L>],
        target_s: f64,
        resolve_films: impl FnOnce(
            &DepositionSimulation<'m>,
        )
            -> Result<BTreeMap<CellId, GasFilm>, ThermalEvaporationError>,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<Self, ThermalEvaporationError> {
        checkpoint()?;
        if sources.len() > initial.cells.len() {
            return Err(DepositionError::InvalidInput("too many evaporating cell sources").into());
        }
        let heat = BackwardEulerStep::new(initial, target_s, &mut checkpoint)?;
        let mut selected = BTreeMap::new();
        for (index, source) in sources.iter().enumerate() {
            if index % 256 == 0 {
                checkpoint()?;
            }
            LiquidBranch::new(source.id, source.pure_liquid_enthalpy_range_j_kg)?.check(initial)?;
            if !initial.cells.contains_key(&source.id) {
                return Err(DepositionError::InvalidInput(
                    "evaporation requires an active liquid cell",
                )
                .into());
            }
            if selected.insert(source.id, ()).is_some() {
                return Err(
                    DepositionError::InvalidInput("duplicate evaporating cell source").into(),
                );
            }
        }
        checkpoint()?;
        let films = resolve_films(initial)?;
        validate_films(selected.keys().copied(), &films, &mut checkpoint)?;
        let bindings: Vec<_> = sources
            .iter()
            .map(|source| EvaporatingCellBinding {
                id: source.id,
                law: source.law,
                film: films[&source.id],
                pure_liquid_enthalpy_range_j_kg: source.pure_liquid_enthalpy_range_j_kg,
            })
            .collect();
        Self::prepare(initial, &bindings, target_s, Some(heat), checkpoint)
    }

    fn prepare(
        initial: &'s DepositionSimulation<'m>,
        bindings: &[EvaporatingCellBinding<'l, L>],
        target_s: f64,
        prepared_heat: Option<BackwardEulerStep<'s, 'm>>,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<Self, ThermalEvaporationError> {
        checkpoint()?;
        if bindings.len() > initial.cells.len() {
            return Err(DepositionError::InvalidInput("too many evaporating cell bindings").into());
        }
        let mut surfaces = BTreeMap::new();
        for (index, binding) in bindings.iter().enumerate() {
            if index % 256 == 0 {
                checkpoint()?;
            }
            binding.film.validate()?;
            let branch = LiquidBranch::new(binding.id, binding.pure_liquid_enthalpy_range_j_kg)?;
            branch.check(initial)?;
            let cell = initial
                .cells
                .get(&binding.id)
                .ok_or(DepositionError::InvalidInput(
                    "evaporation requires an active liquid cell",
                ))?;
            let area_m2 = finite(
                cell.geometry.surface_area_m2 - cell.boundary.reservoir.map_or(0.0, |r| r.area_m2),
            )?;
            if !nonnegative(area_m2) {
                return Err(DepositionError::InvalidInput("contacts exceed surface area").into());
            }
            if surfaces
                .insert(
                    binding.id,
                    PreparedSurface {
                        binding: *binding,
                        branch,
                        area_m2,
                        initial_transfer_kg: 0.0,
                    },
                )
                .is_some()
            {
                return Err(
                    DepositionError::InvalidInput("duplicate evaporating cell binding").into(),
                );
            }
        }
        // Same ordered area subtraction and authoritative asymmetric patches as
        // exposed_area_m2, but visit each original contact only once for all cells.
        for (index, link) in initial.links.values().enumerate() {
            if index % 256 == 0 {
                checkpoint()?;
            }
            let covered = initial.contact_areas(link);
            for (id, area) in [(link.a, covered[0]), (link.b, covered[1])] {
                if let Some(surface) = surfaces.get_mut(&id) {
                    surface.area_m2 = finite(surface.area_m2 - area)?;
                    if !nonnegative(surface.area_m2) {
                        return Err(
                            DepositionError::InvalidInput("contacts exceed surface area").into(),
                        );
                    }
                }
            }
        }
        let heat = match prepared_heat {
            Some(heat) => heat,
            None => BackwardEulerStep::new(initial, target_s, &mut checkpoint)?,
        };
        for surface in surfaces.values_mut() {
            checkpoint()?;
            let binding = surface.binding;
            let p = properties(binding.law, initial.temperature_k(binding.id)?)?;
            let rate = flux(p, binding.film)?.mass_flux_kg_m2_s;
            surface.initial_transfer_kg = finite(rate * surface.area_m2 * heat.duration_s())?;
        }
        checkpoint()?;
        Ok(Self {
            heat,
            initial,
            surfaces,
        })
    }

    pub fn cell_ids(&self) -> &[CellId] {
        self.heat.cell_ids()
    }
    pub const fn duration_s(&self) -> f64 {
        self.heat.duration_s()
    }
    pub const fn target_s(&self) -> f64 {
        self.heat.target_s()
    }
    pub fn resolvable_increments_j(&self, energy_increments_j: &[f64]) -> Vec<f64> {
        self.heat.resolvable_increments_j(energy_increments_j)
    }

    /// Resolve all liquid inventories before calling the endpoint source once.
    /// `thermal` retains original masses for the BE storage coordinates; the
    /// canonical liquid map supplies every proposed new mass. Return independent
    /// sources/boundary exchanges/transfers only; phase heat is added here.
    /// Every evaluation uses one physical heat-rate pass and one native residual
    /// pass over the shared graph, with no per-surface heat solve or split.
    pub fn evaluate<T, E>(
        &self,
        energy_increments_j: &[f64],
        resolve: impl FnOnce(
            &DepositionSimulation<'m>,
            &BTreeMap<CellId, LiquidSurfaceState>,
        ) -> Result<BoundaryExchangeStepCoefficients<T>, E>,
        checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<EvaporatingNetworkBackwardEulerTrial<'m, T>, E>
    where
        E: From<ThermalEvaporationError> + From<DepositionError>,
    {
        self.evaluate_resolved(
            energy_increments_j,
            |_| Ok(None),
            |thermal, liquids, _| resolve(thermal, liquids).map(|heat| (heat, None)),
            checkpoint,
        )
    }

    /// Resolve all endpoint films BEFORE evaluating any interface rate. The
    /// heat callback sees the exact validated films and proposed liquid masses
    /// used in this residual. Returned maps must match the selected cells
    /// exactly, including an empty map for a heat-only step.
    /// The thermal view carries endpoint time/h/T and the interval's INITIAL
    /// masses for energy coordinates; proposed masses are in the liquid map.
    ///
    /// Film sources may depend on time and all endpoint caloric states; their
    /// accepted history must remain immutable during speculative calls. This
    /// remains a prescribed-film thermal balance: a finite-gas owner must also
    /// close its species, energy and upstream-vapor enthalpy accounts.
    pub fn evaluate_with_films<T, E>(
        &self,
        energy_increments_j: &[f64],
        resolve_films: impl FnOnce(&DepositionSimulation<'m>) -> Result<BTreeMap<CellId, GasFilm>, E>,
        resolve_heat: impl FnOnce(
            &DepositionSimulation<'m>,
            &BTreeMap<CellId, LiquidSurfaceState>,
            &BTreeMap<CellId, GasFilm>,
        ) -> Result<BoundaryExchangeStepCoefficients<T>, E>,
        checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<EvaporatingNetworkBackwardEulerTrial<'m, T>, E>
    where
        E: From<ThermalEvaporationError> + From<DepositionError>,
    {
        self.evaluate_resolved(
            energy_increments_j,
            |thermal| resolve_films(thermal).map(Some),
            |thermal, liquids, films| {
                resolve_heat(thermal, liquids, films.expect("resolved endpoint films"))
                    .map(|heat| (heat, None))
            },
            checkpoint,
        )
    }

    /// Evaluate original heat and phase mass transfer with source-owned
    /// complete vapor enthalpy. Films are resolved first, then the callback
    /// sees all proposed liquid masses and those exact films. It supplies
    /// physical heat and a finite excess for EVERY selected cell, including
    /// explicit zero values. The thermal view retains initial masses for the
    /// native energy coordinates; no gas inventory or datum is guessed here.
    ///
    /// The reduced cell balance is `m0*(h1-h0) + dm*L + dm*excess = Q`.
    /// Positive dm exports vapor and negative dm imports it. Latent heat stays
    /// `dm*L`; full carried enthalpy becomes `dm*(h1+L+excess)`. Independent
    /// heat-rate evaluation and every original local/global admission gate
    /// are retained. No separate thermal solve or source query is added.
    ///
    /// This event-free trial has a distinct audit and acceptance type. Its
    /// correction must enter an enclosing finite-gas balance explicitly;
    /// existing prescribed-film adaptive timelines do not accept this type.
    pub fn evaluate_with_vapor_enthalpy<T, E>(
        &self,
        energy_increments_j: &[f64],
        resolve_films: impl FnOnce(&DepositionSimulation<'m>) -> Result<BTreeMap<CellId, GasFilm>, E>,
        resolve_heat_and_vapor: impl FnOnce(
            &DepositionSimulation<'m>,
            &BTreeMap<CellId, LiquidSurfaceState>,
            &BTreeMap<CellId, GasFilm>,
        ) -> Result<VaporEnthalpyStepCoefficients<T>, E>,
        checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<EvaporatingNetworkVaporEnthalpyTrial<'m, T>, E>
    where
        E: From<ThermalEvaporationError> + From<DepositionError>,
    {
        let inner = self.evaluate_resolved(
            energy_increments_j,
            |thermal| resolve_films(thermal).map(Some),
            |thermal, liquids, films| {
                let resolved = resolve_heat_and_vapor(
                    thermal,
                    liquids,
                    films.expect("resolved endpoint films"),
                )?;
                Ok((resolved.heat, Some(resolved.vapor_enthalpy_excess_j_kg)))
            },
            checkpoint,
        )?;
        Ok(EvaporatingNetworkVaporEnthalpyTrial { inner })
    }

    fn evaluate_resolved<T, E>(
        &self,
        energy_increments_j: &[f64],
        resolve_films: impl FnOnce(
            &DepositionSimulation<'m>,
        ) -> Result<Option<BTreeMap<CellId, GasFilm>>, E>,
        resolve: impl FnOnce(
            &DepositionSimulation<'m>,
            &BTreeMap<CellId, LiquidSurfaceState>,
            Option<&BTreeMap<CellId, GasFilm>>,
        ) -> Result<
            (
                BoundaryExchangeStepCoefficients<T>,
                Option<BTreeMap<CellId, f64>>,
            ),
            E,
        >,
        checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<EvaporatingNetworkBackwardEulerTrial<'m, T>, E>
    where
        E: From<ThermalEvaporationError> + From<DepositionError>,
    {
        let checkpoint = RefCell::new(checkpoint);
        let dt = self.duration_s();
        let heat = self.heat.evaluate_with_boundary_exchanges::<_, E>(
            energy_increments_j,
            |thermal| {
                let mut poll = || checkpoint.borrow_mut()();
                poll()?;
                let films = resolve_films(thermal)?;
                if let Some(films) = &films {
                    validate_films(self.surfaces.keys().copied(), films, &mut poll)?;
                }
                let mut liquids = BTreeMap::new();
                let mut receipts = BTreeMap::new();
                for (&id, prepared) in &self.surfaces {
                    poll()?;
                    prepared.branch.check(thermal)?;
                    let old = &self.initial.cells[&id];
                    let h = thermal.cells[&id].specific_h_j_kg;
                    let temperature_k = thermal.temperature_k(id)?;
                    let p = properties(prepared.binding.law, temperature_k)
                        .map_err(ThermalEvaporationError::from)?;
                    let film = films
                        .as_ref()
                        .map_or(prepared.binding.film, |films| films[&id]);
                    let rate = flux(p, film).map_err(ThermalEvaporationError::from)?;
                    let raw_dm = finite(rate.mass_flux_kg_m2_s * prepared.area_m2 * dt)?;
                    let mass = finite(old.mass_kg - raw_dm)?;
                    if mass <= 0.0 {
                        return Err(ThermalEvaporationError::from(
                            EvaporationError::InventoryLimit {
                                maximum_duration_s: (dt
                                    * (old.mass_kg / prepared.initial_transfer_kg.max(raw_dm)))
                                .next_down()
                                .max(0.0),
                            },
                        )
                        .into());
                    }
                    let dm = finite(old.mass_kg - mass)?;
                    if raw_dm != 0.0 && dm == 0.0 {
                        return Err(ThermalEvaporationError::from(
                            EvaporationError::Unrepresentable(
                                "nonzero transfer rounds away in liquid mass",
                            ),
                        )
                        .into());
                    }
                    let liquid = LiquidSurfaceState {
                        mass_kg: mass,
                        specific_enthalpy_j_kg: h,
                        area_m2: prepared.area_m2,
                    };
                    let mut proposed_cell = thermal.cells[&id].clone();
                    proposed_cell.mass_kg = mass;
                    thermal.materials[&proposed_cell.material_id].validate_cell(&proposed_cell)?;
                    let latent_energy_j = finite(dm * p.latent_heat_j_kg)?;
                    if latent_energy_j != 0.0 && finite(latent_energy_j / dt)? == 0.0 {
                        return Err(ThermalEvaporationError::from(
                            EvaporationError::Unrepresentable(
                                "nonzero latent transfer underflows to zero thermal power",
                            ),
                        )
                        .into());
                    }
                    liquids.insert(id, liquid);
                    receipts.insert(
                        id,
                        SurfaceStep {
                            liquid,
                            transferred_mass_kg: dm,
                            latent_energy_j,
                            vapor_energy_j: finite(finite(dm * h)? + latent_energy_j)?,
                            supplied_heat_j: 0.0,
                            energy_residual_j: 0.0,
                            mass_rate_residual_kg: finite(dm - raw_dm)?,
                            recession_m: if prepared.area_m2 == 0.0 {
                                0.0
                            } else {
                                finite(dm / prepared.area_m2 / p.liquid_density_kg_m3)?
                            },
                            surface_temperature_k: temperature_k,
                        },
                    );
                }
                let (resolved_heat, excess) = resolve(thermal, &liquids, films.as_ref())?;
                let mut correction_energy = BTreeMap::new();
                if let Some(excess) = &excess {
                    validate_enthalpy_excess(self.surfaces.keys().copied(), excess, &mut poll)?;
                    for (&id, surface) in &mut receipts {
                        poll()?;
                        let correction_j = finite(surface.transferred_mass_kg * excess[&id])?;
                        surface.vapor_energy_j = finite(surface.vapor_energy_j + correction_j)?;
                        correction_energy.insert(id, correction_j);
                    }
                }
                let BoundaryExchangeStepCoefficients {
                    transport:
                        TransportStepCoefficients {
                            heat:
                                BoundaryStepCoefficients {
                                    bulk:
                                        StepCoefficients {
                                            mut sources_w,
                                            conductances_w_k,
                                            payload,
                                        },
                                    boundaries,
                                },
                            transfers,
                        },
                    exchanges,
                } = resolved_heat;
                let mut physical = thermal.clone();
                physical.set_step_conductances(&conductances_w_k, &mut poll)?;
                let mut rates = physical.heat_rates_with_boundary_exchanges(
                    f64::MAX,
                    &boundaries,
                    &exchanges,
                    &mut poll,
                )?;
                physical.add_rate_sources(&mut rates, &sources_w, &mut poll)?;
                physical.add_step_transfers(&mut rates, &transfers, &mut poll)?;
                let external_heat_j = finite(dt * rates.external_w.value()?)?;
                for (&id, surface) in &mut receipts {
                    poll()?;
                    let old = &self.initial.cells[&id];
                    surface.supplied_heat_j = finite(dt * rates.watts[&id].value()?)?;
                    let mut residual = EnergySum::default();
                    residual.add(finite(
                        old.mass_kg
                            * finite(surface.liquid.specific_enthalpy_j_kg - old.specific_h_j_kg)?,
                    )?)?;
                    residual.add(surface.latent_energy_j)?;
                    if let Some(&correction_j) = correction_energy.get(&id) {
                        residual.add(correction_j)?;
                    }
                    residual.add(-surface.supplied_heat_j)?;
                    surface.energy_residual_j = residual.value()?;
                    physical
                        .cells
                        .get_mut(&id)
                        .expect("active liquid cell")
                        .mass_kg = surface.liquid.mass_kg;
                    let phase_energy_j = match correction_energy.get(&id) {
                        Some(&correction_j) => finite(surface.latent_energy_j + correction_j)?,
                        None => surface.latent_energy_j,
                    };
                    let phase_power_w = finite(phase_energy_j / dt)?;
                    if phase_energy_j != 0.0 && phase_power_w == 0.0 {
                        return Err(ThermalEvaporationError::from(
                            EvaporationError::Unrepresentable(
                                "nonzero vapor enthalpy transfer underflows to zero thermal power",
                            ),
                        )
                        .into());
                    }
                    if phase_power_w != 0.0 {
                        let source = sources_w.entry(id).or_insert(0.0);
                        *source = finite(*source - phase_power_w)?;
                    }
                }
                let mut stored_change = EnergySum::default();
                let mut energy_residual = EnergySum::default();
                let mut mass_residual = EnergySum::default();
                for (index, (&id, cell)) in physical.cells.iter().enumerate() {
                    if index % 256 == 0 {
                        poll()?;
                    }
                    let old = &self.initial.cells[&id];
                    let dh = finite(cell.specific_h_j_kg - old.specific_h_j_kg)?;
                    stored_change.add(finite(cell.mass_kg * dh)?)?;
                    energy_residual.add(finite(old.mass_kg * dh)?)?;
                    mass_residual.add(finite(cell.mass_kg - old.mass_kg)?)?;
                }
                for (&id, surface) in &receipts {
                    poll()?;
                    stored_change.add(-finite(
                        surface.transferred_mass_kg * self.initial.cells[&id].specific_h_j_kg,
                    )?)?;
                    energy_residual.add(surface.latent_energy_j)?;
                    if let Some(&correction_j) = correction_energy.get(&id) {
                        energy_residual.add(correction_j)?;
                    }
                    mass_residual.add(surface.transferred_mass_kg)?;
                }
                energy_residual.add(-external_heat_j)?;
                let audit = EvaporatingNetworkBackwardEulerAudit {
                    external_heat_j,
                    stored_energy_change_j: stored_change.value()?,
                    mass_balance_residual_kg: mass_residual.value()?,
                    energy_residual_j: energy_residual.value()?,
                    surfaces: receipts,
                };
                poll()?;
                Ok(BoundaryExchangeStepCoefficients {
                    transport: TransportStepCoefficients {
                        heat: BoundaryStepCoefficients {
                            bulk: StepCoefficients {
                                sources_w,
                                conductances_w_k,
                                payload: Endpoint {
                                    payload,
                                    audit: EvaporatingNetworkVaporEnthalpyAudit {
                                        thermal: audit,
                                        vapor_enthalpy_correction_j: correction_energy,
                                    },
                                },
                            },
                            boundaries,
                        },
                        transfers,
                    },
                    exchanges,
                })
            },
            || checkpoint.borrow_mut()(),
        )?;
        checkpoint.borrow_mut()()?;
        Ok(EvaporatingNetworkBackwardEulerTrial { heat })
    }
}

fn validate_enthalpy_excess(
    ids: impl ExactSizeIterator<Item = CellId>,
    excess: &BTreeMap<CellId, f64>,
    checkpoint: &mut impl FnMut() -> Result<(), DepositionError>,
) -> Result<(), ThermalEvaporationError> {
    if excess.len() != ids.len() {
        return Err(DepositionError::InvalidInput(
            "vapor enthalpy excess must exactly match selected liquid cells",
        )
        .into());
    }
    for (index, id) in ids.enumerate() {
        if index % 256 == 0 {
            checkpoint()?;
        }
        if !excess.get(&id).is_some_and(|value| value.is_finite()) {
            return Err(DepositionError::InvalidInput(
                "finite vapor enthalpy excess required for every selected liquid cell",
            )
            .into());
        }
    }
    checkpoint()?;
    Ok(())
}

fn validate_films(
    ids: impl ExactSizeIterator<Item = CellId>,
    films: &BTreeMap<CellId, GasFilm>,
    checkpoint: &mut impl FnMut() -> Result<(), DepositionError>,
) -> Result<(), ThermalEvaporationError> {
    if films.len() != ids.len() {
        return Err(DepositionError::InvalidInput(
            "resolved films must exactly match selected liquid cells",
        )
        .into());
    }
    for (index, id) in ids.enumerate() {
        if index % 256 == 0 {
            checkpoint()?;
        }
        films
            .get(&id)
            .ok_or(DepositionError::InvalidInput(
                "resolved films must exactly match selected liquid cells",
            ))?
            .validate()?;
    }
    checkpoint()?;
    Ok(())
}

struct Endpoint<T> {
    payload: T,
    audit: EvaporatingNetworkVaporEnthalpyAudit,
}

/// One speculative heat and multi-surface endpoint, with no published history.
pub struct EvaporatingNetworkBackwardEulerTrial<'m, T> {
    heat: BackwardEulerTrial<'m, Endpoint<T>>,
}
impl<'m, T> EvaporatingNetworkBackwardEulerTrial<'m, T> {
    pub fn residual_j(&self) -> &[f64] {
        self.heat.residual_j()
    }
    pub fn maximum_cell_residual_j(&self) -> f64 {
        self.audit()
            .surfaces
            .values()
            .fold(self.heat.maximum_cell_residual_j(), |maximum, surface| {
                maximum.max(surface.energy_residual_j.abs())
            })
    }
    pub const fn audit(&self) -> &EvaporatingNetworkBackwardEulerAudit {
        &self.heat.payload().audit.thermal
    }
    pub fn payload(&self) -> &T {
        &self.heat.payload().payload
    }

    /// Preserve native local/total energy gates, then separately gate each
    /// surface's energy and represented mass rate and the network's open-system
    /// balances. All original material validations precede the mass installation;
    /// final cancellation still discards the complete candidate.
    pub fn accept(
        self,
        cell_tolerance_j: f64,
        total_tolerance_j: f64,
        mass_tolerance_kg: f64,
        checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<
        (
            DepositionSimulation<'m>,
            T,
            EvaporatingNetworkBackwardEulerAudit,
        ),
        ThermalEvaporationError,
    > {
        let (next, payload, audit) = self.accept_complete(
            cell_tolerance_j,
            total_tolerance_j,
            mass_tolerance_kg,
            checkpoint,
        )?;
        Ok((next, payload, audit.thermal))
    }

    fn accept_complete(
        self,
        cell_tolerance_j: f64,
        total_tolerance_j: f64,
        mass_tolerance_kg: f64,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<
        (
            DepositionSimulation<'m>,
            T,
            EvaporatingNetworkVaporEnthalpyAudit,
        ),
        ThermalEvaporationError,
    > {
        checkpoint()?;
        if !nonnegative(cell_tolerance_j)
            || !nonnegative(total_tolerance_j)
            || !nonnegative(mass_tolerance_kg)
        {
            return Err(DepositionError::InvalidInput(
                "implicit evaporation tolerances must be finite nonnegative joules and kilograms",
            )
            .into());
        }
        let audit = self.audit();
        for surface in audit.surfaces.values() {
            checkpoint()?;
            if surface.energy_residual_j.abs() > cell_tolerance_j {
                return Err(DepositionError::EnergyBalance(surface.energy_residual_j).into());
            }
            if surface.mass_rate_residual_kg.abs() > mass_tolerance_kg {
                return Err(ThermalEvaporationError::MassBalance {
                    residual_kg: surface.mass_rate_residual_kg,
                    tolerance_kg: mass_tolerance_kg,
                });
            }
        }
        if audit.energy_residual_j.abs() > total_tolerance_j {
            return Err(DepositionError::EnergyBalance(audit.energy_residual_j).into());
        }
        if audit.mass_balance_residual_kg.abs() > mass_tolerance_kg {
            return Err(ThermalEvaporationError::MassBalance {
                residual_kg: audit.mass_balance_residual_kg,
                tolerance_kg: mass_tolerance_kg,
            });
        }
        let (mut next, endpoint, _) =
            self.heat
                .accept(cell_tolerance_j, total_tolerance_j, &mut checkpoint)?;
        for (&id, surface) in &endpoint.audit.thermal.surfaces {
            checkpoint()?;
            let mut proposed = next.cells[&id].clone();
            proposed.mass_kg = surface.liquid.mass_kg;
            next.materials[&proposed.material_id].validate_cell(&proposed)?;
        }
        checkpoint()?;
        for (&id, surface) in &endpoint.audit.thermal.surfaces {
            next.cells
                .get_mut(&id)
                .expect("validated active liquid cell")
                .mass_kg = surface.liquid.mass_kg;
        }
        checkpoint()?;
        Ok((next, endpoint.payload, endpoint.audit))
    }
}

/// Complete source-owned vapor enthalpy on the same speculative heat endpoint.
/// This cannot be passed to the existing prescribed-film timeline machinery.
pub struct EvaporatingNetworkVaporEnthalpyTrial<'m, T> {
    inner: EvaporatingNetworkBackwardEulerTrial<'m, T>,
}
impl<'m, T> EvaporatingNetworkVaporEnthalpyTrial<'m, T> {
    pub fn residual_j(&self) -> &[f64] {
        self.inner.residual_j()
    }
    pub fn maximum_cell_residual_j(&self) -> f64 {
        self.inner.maximum_cell_residual_j()
    }
    pub const fn audit(&self) -> &EvaporatingNetworkVaporEnthalpyAudit {
        &self.inner.heat.payload().audit
    }
    pub fn payload(&self) -> &T {
        self.inner.payload()
    }

    /// Preserve every native and corrected local/global balance, mass-rate,
    /// material and cancellation gate before installing any proposed mass.
    pub fn accept(
        self,
        cell_tolerance_j: f64,
        total_tolerance_j: f64,
        mass_tolerance_kg: f64,
        checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<
        (
            DepositionSimulation<'m>,
            T,
            EvaporatingNetworkVaporEnthalpyAudit,
        ),
        ThermalEvaporationError,
    > {
        self.inner.accept_complete(
            cell_tolerance_j,
            total_tolerance_j,
            mass_tolerance_kg,
            checkpoint,
        )
    }
}
