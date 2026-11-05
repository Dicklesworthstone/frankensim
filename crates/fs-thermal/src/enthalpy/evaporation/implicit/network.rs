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
        let heat = BackwardEulerStep::new(initial, target_s, &mut checkpoint)?;
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
        let checkpoint = RefCell::new(checkpoint);
        let dt = self.duration_s();
        let heat = self.heat.evaluate_with_boundary_exchanges::<_, E>(
            energy_increments_j,
            |thermal| {
                let mut poll = || checkpoint.borrow_mut()();
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
                    let rate =
                        flux(p, prepared.binding.film).map_err(ThermalEvaporationError::from)?;
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
                } = resolve(thermal, &liquids)?;
                let mut physical = thermal.clone();
                physical.set_step_conductances(&conductances_w_k, &mut poll)?;
                physical.add_step_sources(&sources_w, &mut poll)?;
                let mut rates = physical.heat_rates_with_boundary_exchanges(
                    f64::MAX,
                    &boundaries,
                    &exchanges,
                    &mut poll,
                )?;
                physical.add_step_transfers(&mut rates, &transfers, &mut poll)?;
                let external_heat_j = finite(dt * rates.external_w)?;
                for (&id, surface) in &mut receipts {
                    poll()?;
                    let old = &self.initial.cells[&id];
                    surface.supplied_heat_j = finite(dt * rates.watts[&id])?;
                    let mut residual = EnergySum::default();
                    residual.add(finite(
                        old.mass_kg
                            * finite(surface.liquid.specific_enthalpy_j_kg - old.specific_h_j_kg)?,
                    )?)?;
                    residual.add(surface.latent_energy_j)?;
                    residual.add(-surface.supplied_heat_j)?;
                    surface.energy_residual_j = residual.value()?;
                    physical
                        .cells
                        .get_mut(&id)
                        .expect("active liquid cell")
                        .mass_kg = surface.liquid.mass_kg;
                    let phase_power_w = finite(surface.latent_energy_j / dt)?;
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
                                payload: Endpoint { payload, audit },
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

struct Endpoint<T> {
    payload: T,
    audit: EvaporatingNetworkBackwardEulerAudit,
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
        &self.heat.payload().audit
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
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<
        (
            DepositionSimulation<'m>,
            T,
            EvaporatingNetworkBackwardEulerAudit,
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
        for (&id, surface) in &endpoint.audit.surfaces {
            checkpoint()?;
            let mut proposed = next.cells[&id].clone();
            proposed.mass_kg = surface.liquid.mass_kg;
            next.materials[&proposed.material_id].validate_cell(&proposed)?;
        }
        checkpoint()?;
        for (&id, surface) in &endpoint.audit.surfaces {
            next.cells
                .get_mut(&id)
                .expect("validated active liquid cell")
                .mass_kg = surface.liquid.mass_kg;
        }
        checkpoint()?;
        Ok((next, endpoint.payload, endpoint.audit))
    }
}
