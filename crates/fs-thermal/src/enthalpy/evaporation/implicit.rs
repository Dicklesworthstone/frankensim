//! Simultaneous endpoint heat transport and one liquid/vapor mass exchange.
//!
//! The original backward-Euler owner evaluates every cell, contact and thermal
//! boundary at one endpoint. For the evaporating cell its storage equation is
//! `m0*(h1-h0) + dm*L(T1) = dt*P1`. The retained mass `m0` is a coordinate for
//! this reduced open-system balance, not the proposed endpoint inventory.
//! Acceptance installs `m1 = m0-dm` and accounts for the full vapor enthalpy
//! `dm*(h1+L1)` once. There is no sensible-only intermediate temperature.
//!
//! This adapter supplies residuals, not a nonlinear solver or an error estimate.
//! Callers own bounded nonlinear work, timestep refinement and the prescribed
//! film's mass/energy account. Geometry and exposure stay fixed over the step;
//! endpoint events remain pending until all coupled owners accept the interval.

use super::{LiquidBranch, ThermalEvaporationError};
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

/// Bounded automatic nonlinear solves of this original endpoint residual.
#[path = "implicit/solve.rs"]
pub mod solve;

#[path = "implicit/timeline.rs"]
pub mod timeline;
pub use timeline::{EvaporationImplicitTimelineControls, EvaporationImplicitTimelineReport};
pub use solve::{EvaporationImplicitControls, EvaporationImplicitError, EvaporationImplicitReport};

/// Physical open-system account for one interval, before endpoint events.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct EvaporatingBackwardEulerAudit {
    /// Original sources and prescribed thermal baths, positive inward [J].
    /// Internal contacts/transfers and phase enthalpy are excluded.
    pub external_heat_j: f64,
    /// Actual liquid/network storage change, including the changed mass [J].
    pub stored_energy_change_j: f64,
    /// Final network mass plus vapor transfer minus initial network mass [kg].
    pub mass_balance_residual_kg: f64,
    /// Storage change plus full vapor enthalpy minus external heat [J].
    /// Evaluated after analytically cancelling the common carried-liquid datum;
    /// summing rounded report fields can lose accuracy at large datum offsets.
    pub energy_residual_j: f64,
    /// Local phase receipt. Its heat includes the cell's signed contact heat;
    /// `external_heat_j` above excludes those conservative internal exchanges.
    pub surface: SurfaceStep,
}

/// An event-free step on the original thermal topology and caloric laws.
/// Unknowns are `m0*(h1-h0)` in `cell_ids()` order, including every active cell.
pub struct EvaporatingBackwardEulerStep<'s, 'm, 'l, L: EvaporationLaw + ?Sized> {
    heat: BackwardEulerStep<'s, 'm>,
    initial: &'s DepositionSimulation<'m>,
    law: &'l L,
    film: GasFilm,
    branch: LiquidBranch,
    area_m2: f64,
    initial_transfer_kg: f64,
}

impl<'s, 'm, 'l, L: EvaporationLaw + ?Sized> EvaporatingBackwardEulerStep<'s, 'm, 'l, L> {
    /// Prepare one positive-duration interval after processing start-time events.
    /// The caller explicitly admits the pure-liquid caloric interval. A missing
    /// liquid cell, initially boiling source or unsupported original law refuses.
    /// Events strictly inside the interval refuse; endpoint events stay pending.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        initial: &'s DepositionSimulation<'m>,
        id: CellId,
        law: &'l L,
        film: GasFilm,
        target_s: f64,
        pure_liquid_enthalpy_range_j_kg: [f64; 2],
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<Self, ThermalEvaporationError> {
        checkpoint()?;
        film.validate()?;
        let branch = LiquidBranch::new(id, pure_liquid_enthalpy_range_j_kg)?;
        branch.check(initial)?;
        let area_m2 = initial.exposed_area_m2(id)?;
        let heat = BackwardEulerStep::new(initial, target_s, &mut checkpoint)?;
        let p = properties(law, initial.temperature_k(id)?)?;
        let rate = flux(p, film)?.mass_flux_kg_m2_s;
        let initial_transfer_kg = finite(rate * area_m2 * heat.duration_s())?;
        checkpoint()?;
        Ok(Self {
            heat,
            initial,
            law,
            film,
            branch,
            area_m2,
            initial_transfer_kg,
        })
    }

    /// Stable ordering shared by energy coordinates and local joule residuals.
    pub fn cell_ids(&self) -> &[CellId] {
        self.heat.cell_ids()
    }

    pub const fn duration_s(&self) -> f64 {
        self.heat.duration_s()
    }

    pub const fn target_s(&self) -> f64 {
        self.heat.target_s()
    }

    /// Use the original owner's representability projection for solver iterates.
    /// It does not discard any heat from the residual or weaken acceptance.
    pub fn resolvable_increments_j(&self, energy_increments_j: &[f64]) -> Vec<f64> {
        self.heat.resolvable_increments_j(energy_increments_j)
    }

    /// Evaluate one tentative endpoint, with no accepted history changes.
    ///
    /// `resolve` receives original endpoint calorics/topology with the INITIAL
    /// mass retained for BE coordinates, plus the explicit proposed liquid mass
    /// and enthalpy. Mass-dependent sources must use that proposed liquid state.
    /// Return only independently owned heat: the adapter adds phase transfer.
    /// Source, link, bath and conservative-transfer coefficients retain the
    /// meanings of the original `BoundaryExchangeStepCoefficients` ingress.
    ///
    /// Each call performs one interface-property query and two original thermal
    /// rate evaluations: one independent physical heat account and one native
    /// residual including latent exchange. Work is bounded by active topology;
    /// the caller must bound the number of nonlinear evaluations.
    pub fn evaluate<T, E>(
        &self,
        energy_increments_j: &[f64],
        resolve: impl FnOnce(
            &DepositionSimulation<'m>,
            LiquidSurfaceState,
        ) -> Result<BoundaryExchangeStepCoefficients<T>, E>,
        checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<EvaporatingBackwardEulerTrial<'m, T>, E>
    where
        E: From<ThermalEvaporationError> + From<DepositionError>,
    {
        let checkpoint = RefCell::new(checkpoint);
        let id = self.branch.id;
        let old = &self.initial.cells[&id];
        let dt = self.duration_s();
        let heat = self.heat.evaluate_with_boundary_exchanges::<_, E>(
            energy_increments_j,
            |thermal| {
                let mut poll = || checkpoint.borrow_mut()();
                poll()?;
                self.branch.check(thermal)?;
                let temperature_k = thermal.temperature_k(id)?;
                let p =
                    properties(self.law, temperature_k).map_err(ThermalEvaporationError::from)?;
                let rate = flux(p, self.film).map_err(ThermalEvaporationError::from)?;
                let raw_dm = finite(rate.mass_flux_kg_m2_s * self.area_m2 * dt)?;
                let mass = finite(old.mass_kg - raw_dm)?;
                if mass <= 0.0 {
                    return Err(
                        ThermalEvaporationError::from(EvaporationError::InventoryLimit {
                            maximum_duration_s: (dt
                                * (old.mass_kg / self.initial_transfer_kg.max(raw_dm)))
                            .next_down()
                            .max(0.0),
                        })
                        .into(),
                    );
                }
                let dm = finite(old.mass_kg - mass)?;
                if raw_dm != 0.0 && dm == 0.0 {
                    return Err(
                        ThermalEvaporationError::from(EvaporationError::Unrepresentable(
                            "nonzero transfer rounds away in liquid mass",
                        ))
                        .into(),
                    );
                }
                let h = thermal.cells[&id].specific_h_j_kg;
                let liquid = LiquidSurfaceState {
                    mass_kg: mass,
                    specific_enthalpy_j_kg: h,
                    area_m2: self.area_m2,
                };
                let mut proposed_cell = thermal.cells[&id].clone();
                proposed_cell.mass_kg = mass;
                thermal.materials[&proposed_cell.material_id].validate_cell(&proposed_cell)?;
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
                } = resolve(thermal, liquid)?;
                // This pass measures original physical supply directly; it does
                // not recover it by subtracting a latent sink from mixed totals.
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
                let supplied_heat_j = finite(dt * rates.watts[&id])?;
                let external_heat_j = finite(dt * rates.external_w)?;
                let latent_energy_j = finite(dm * p.latent_heat_j_kg)?;
                let phase_power_w = finite(latent_energy_j / dt)?;
                if latent_energy_j != 0.0 && phase_power_w == 0.0 {
                    return Err(
                        ThermalEvaporationError::from(EvaporationError::Unrepresentable(
                            "nonzero latent transfer underflows to zero thermal power",
                        ))
                        .into(),
                    );
                }
                let mut surface_residual = EnergySum::default();
                surface_residual.add(finite(old.mass_kg * finite(h - old.specific_h_j_kg)?)?)?;
                surface_residual.add(latent_energy_j)?;
                surface_residual.add(-supplied_heat_j)?;
                let surface = SurfaceStep {
                    liquid,
                    transferred_mass_kg: dm,
                    latent_energy_j,
                    vapor_energy_j: finite(finite(dm * h)? + latent_energy_j)?,
                    supplied_heat_j,
                    energy_residual_j: surface_residual.value()?,
                    mass_rate_residual_kg: finite(dm - raw_dm)?,
                    recession_m: if self.area_m2 == 0.0 {
                        0.0
                    } else {
                        finite(dm / self.area_m2 / p.liquid_density_kg_m3)?
                    },
                    surface_temperature_k: temperature_k,
                };
                physical
                    .cells
                    .get_mut(&id)
                    .expect("active liquid cell")
                    .mass_kg = mass;
                let mut stored_change = EnergySum::default();
                let mut energy_residual = EnergySum::default();
                let mut mass_residual = EnergySum::default();
                for (index, (&cell_id, cell)) in physical.cells.iter().enumerate() {
                    if index % 256 == 0 {
                        poll()?;
                    }
                    let old_cell = &self.initial.cells[&cell_id];
                    let dh = finite(cell.specific_h_j_kg - old_cell.specific_h_j_kg)?;
                    stored_change.add(finite(cell.mass_kg * dh)?)?;
                    energy_residual.add(finite(old_cell.mass_kg * dh)?)?;
                    mass_residual.add(finite(cell.mass_kg - old_cell.mass_kg)?)?;
                }
                stored_change.add(-finite(dm * old.specific_h_j_kg)?)?;
                energy_residual.add(latent_energy_j)?;
                energy_residual.add(-external_heat_j)?;
                mass_residual.add(dm)?;
                let audit = EvaporatingBackwardEulerAudit {
                    external_heat_j,
                    stored_energy_change_j: stored_change.value()?,
                    mass_balance_residual_kg: mass_residual.value()?,
                    energy_residual_j: energy_residual.value()?,
                    surface,
                };
                if phase_power_w != 0.0 {
                    let source = sources_w.entry(id).or_insert(0.0);
                    *source = finite(*source - phase_power_w)?;
                }
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
        let audit = heat.payload().audit;
        checkpoint.borrow_mut()()?;
        Ok(EvaporatingBackwardEulerTrial { heat, id, audit })
    }
}

struct Endpoint<T> {
    payload: T,
    audit: EvaporatingBackwardEulerAudit,
}

/// A complete tentative endpoint. No mass or history is published by evaluation.
pub struct EvaporatingBackwardEulerTrial<'m, T> {
    heat: BackwardEulerTrial<'m, Endpoint<T>>,
    id: CellId,
    audit: EvaporatingBackwardEulerAudit,
}

impl<'m, T> EvaporatingBackwardEulerTrial<'m, T> {
    /// Original local energy residuals, including represented latent transfer.
    pub fn residual_j(&self) -> &[f64] {
        self.heat.residual_j()
    }

    /// Includes the independent physical surface balance, so cancellation in
    /// mixed source arithmetic cannot hide a local defect.
    pub fn maximum_cell_residual_j(&self) -> f64 {
        self.heat
            .maximum_cell_residual_j()
            .max(self.audit.surface.energy_residual_j.abs())
    }

    pub const fn audit(&self) -> &EvaporatingBackwardEulerAudit {
        &self.audit
    }

    pub fn payload(&self) -> &T {
        &self.heat.payload().payload
    }

    /// Accept only after native local/total energy gates, the independent open
    /// energy account, represented-rate and whole-network mass gates all pass.
    /// Tolerances are finite nonnegative absolute joules and kilograms. Events
    /// at the target remain pending; reference geometry and authored inputs are
    /// unchanged. Algebraic acceptance does not establish timestep accuracy.
    pub fn accept(
        self,
        cell_tolerance_j: f64,
        total_tolerance_j: f64,
        mass_tolerance_kg: f64,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<(DepositionSimulation<'m>, T, EvaporatingBackwardEulerAudit), ThermalEvaporationError>
    {
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
        if self.audit.surface.energy_residual_j.abs() > cell_tolerance_j {
            return Err(
                DepositionError::EnergyBalance(self.audit.surface.energy_residual_j).into(),
            );
        }
        if self.audit.energy_residual_j.abs() > total_tolerance_j {
            return Err(DepositionError::EnergyBalance(self.audit.energy_residual_j).into());
        }
        for residual_kg in [
            self.audit.mass_balance_residual_kg,
            self.audit.surface.mass_rate_residual_kg,
        ] {
            if residual_kg.abs() > mass_tolerance_kg {
                return Err(ThermalEvaporationError::MassBalance {
                    residual_kg,
                    tolerance_kg: mass_tolerance_kg,
                });
            }
        }
        let (mut next, endpoint, _) =
            self.heat
                .accept(cell_tolerance_j, total_tolerance_j, &mut checkpoint)?;
        let cell = next
            .cells
            .get_mut(&self.id)
            .expect("validated active liquid cell");
        cell.mass_kg = endpoint.audit.surface.liquid.mass_kg;
        next.materials[&cell.material_id].validate_cell(cell)?;
        checkpoint()?;
        Ok((next, endpoint.payload, self.audit))
    }
}

#[cfg(test)]
#[path = "implicit/tests.rs"]
mod tests;
