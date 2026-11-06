//! Native deposition events around the shared adaptive finite-gas march.
//!
//! This module owns surface activation and native event admission only. Heat,
//! phase exchange, gas reconstruction, time accuracy and whole-window accounts
//! remain with the existing coupled endpoint and adaptive owners.

use std::collections::{BTreeMap, BTreeSet};

use super::super::admit_liquid_source;
use super::adaptive::{
    self, AdaptiveImplicitGasCandidate, AdaptiveImplicitGasControls, AdaptiveImplicitGasError,
};
use super::*;
use crate::thermochemical_gas::energy::parcel::transfer::species::boundary::admit_catalog;
use fs_thermal::enthalpy::{
    DepositionObserver, EnergyAudit, Event, StepConfig, ThermalStorageBasis,
};

/// Advance the joint heat/phase/gas solution through native deposition events.
///
/// The shared adaptive owner clips comparisons at each event boundary. Native
/// births, contact/covering changes, separations and boundary changes apply in
/// their authored order, including events at the initial and target times.
/// Selected interfaces activate only after birth, and each interval retains
/// its starting exposure. Existing heat cells advance before the first liquid
/// birth while the gas remains idle; empty geometry advances through the native
/// event owner without consuming thermal substeps.
/// The heat callback's liquid map contains only active selected cells and may
/// be empty while unselected native heat cells are present.
///
/// Only accepted fine intervals contribute physical phase exchange. Their
/// numerical work and conservation allowances span the complete request, not
/// individual events. Source discontinuities must coincide with authored events
/// or explicit advance boundaries. At an event endpoint the heat callback sees
/// the pre-event native state and must retain its source's left limit; native
/// Boundary events provide this ordering directly. The returned payload belongs to the
/// last accepted fine interval before any events at its endpoint; a window with
/// no heat interval has no payload. Refusal leaves both original inputs intact.
#[allow(clippy::too_many_arguments, clippy::type_complexity)]
pub fn advance_deposition_with_gas_timeline<'m, 'catalog, T>(
    thermal: &DepositionSimulation<'m>,
    surfaces: &[ImplicitGasSurface<'_>],
    gas: &ClosedGasParcel,
    catalog: &'catalog ThermochemicalGasModel,
    target_s: f64,
    controls: AdaptiveImplicitGasControls,
    resolve_heat: impl FnMut(
        &DepositionSimulation<'m>,
        &BTreeMap<CellId, LiquidSurfaceState>,
        &ClosedGasParcel,
    ) -> Result<BoundaryExchangeStepCoefficients<T>, ImplicitGasError>,
    checkpoint: impl FnMut() -> Result<(), DepositionError>,
) -> Result<AdaptiveImplicitGasCandidate<'m, 'catalog, T>, AdaptiveImplicitGasError> {
    adaptive::advance_window(
        thermal,
        surfaces,
        gas,
        catalog,
        target_s,
        controls,
        true,
        resolve_heat,
        checkpoint,
    )
}

/// Immutable declarations include future births without querying their laws.
pub(super) struct Roster<'l> {
    surfaces: BTreeMap<CellId, ImplicitGasSurface<'l>>,
}

/// Native birth coordinates, retained before later heat or mass exchange.
#[derive(Clone, Copy, Debug)]
pub(super) struct Birth {
    pub mass_kg: f64,
    pub liquid_enthalpy_offset_j_kg: f64,
}

impl<'l> Roster<'l> {
    pub(super) fn new(
        thermal: &DepositionSimulation<'_>,
        surfaces: &[ImplicitGasSurface<'l>],
        gas: &ClosedGasParcel,
        catalog: &ThermochemicalGasModel,
        controls: ImplicitGasControls,
        checkpoint: &mut impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<Self, AdaptiveImplicitGasError> {
        checkpoint()?;
        admit_gas(gas, catalog, controls)?;
        let mut known: BTreeSet<_> = thermal.cells().keys().copied().collect();
        for (index, event) in thermal.pending_events().iter().enumerate() {
            if index % 256 == 0 {
                checkpoint()?;
            }
            if let Event::Birth(cell) = &event.event {
                known.insert(cell.id);
            }
        }
        let mut ordered = BTreeMap::new();
        let mut species = BTreeSet::new();
        for surface in surfaces {
            checkpoint()?;
            let [lo, hi] = surface.pure_liquid_enthalpy_range_j_kg;
            if !known.contains(&surface.cell_id)
                || !lo.is_finite()
                || !hi.is_finite()
                || lo > hi
                || !surface.binding.liquid_enthalpy_offset_j_kg.is_finite()
                || !nonnegative(surface.mass_transfer_coefficient_m_s)
                || !nonnegative(surface.heat_transfer_coefficient_w_m2_k)
                || ordered.insert(surface.cell_id, *surface).is_some()
            {
                return Err(AdaptiveImplicitGasError::Invalid(
                    "unique known liquid cells and finite admitted film/enthalpy data required",
                ));
            }
            species.insert(&surface.binding.species);
        }
        for species in species {
            if !catalog
                .thermodynamics()
                .components()
                .iter()
                .any(|source| source.species() == species)
            {
                return Err(GasEvaporationError::Invalid(
                    "evaporating species requires its complete original gas source",
                )
                .into());
            }
        }
        // Admit every declaration before querying any selected phase law.
        // Future births retain their declaration but have no physical state yet.
        for surface in ordered.values() {
            checkpoint()?;
            if thermal.cells().contains_key(&surface.cell_id) {
                admit_surface(*surface, thermal, catalog, controls)?;
            }
        }
        checkpoint()?;
        Ok(Self { surfaces: ordered })
    }

    pub(super) fn active(
        &self,
        thermal: &DepositionSimulation<'_>,
        checkpoint: &mut impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<Vec<ImplicitGasSurface<'l>>, AdaptiveImplicitGasError> {
        let mut active = Vec::new();
        for (index, (&id, &surface)) in self.surfaces.iter().enumerate() {
            if index % 256 == 0 {
                checkpoint()?;
            }
            if thermal.cells().contains_key(&id) {
                active.push(surface);
            }
        }
        checkpoint()?;
        Ok(active)
    }

    #[allow(clippy::too_many_arguments)]
    pub(super) fn advance_events<'m>(
        &self,
        thermal: &DepositionSimulation<'m>,
        gas: &ClosedGasParcel,
        catalog: &ThermochemicalGasModel,
        target_s: f64,
        controls: ImplicitGasControls,
        checkpoint: &mut impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<
        (
            DepositionSimulation<'m>,
            BTreeMap<CellId, Birth>,
            EnergyAudit,
        ),
        AdaptiveImplicitGasError,
    > {
        checkpoint()?;
        admit_gas(gas, catalog, controls)?;
        let observer = BirthObserver {
            roster: self,
            catalog,
            controls,
            births: BTreeMap::new(),
        };
        let (thermal, observer, audit) = thermal.advance_observed_limited(
            &observer,
            target_s,
            StepConfig {
                max_dt_s: 1.0,
                energy_tolerance_j: controls.thermal.total_energy_tolerance_j,
            },
            0,
            &mut *checkpoint,
        )?;
        checkpoint()?;
        Ok((thermal, observer.births, audit))
    }
}

#[derive(Clone)]
struct BirthObserver<'a, 'l> {
    roster: &'a Roster<'l>,
    catalog: &'a ThermochemicalGasModel,
    controls: ImplicitGasControls,
    births: BTreeMap<CellId, Birth>,
}

impl DepositionObserver for BirthObserver<'_, '_> {
    type Error = AdaptiveImplicitGasError;

    fn event(
        &mut self,
        event: &Event,
        thermal: &DepositionSimulation<'_>,
    ) -> Result<(), Self::Error> {
        if let Event::Birth(cell) = event {
            let offset = if let Some(&surface) = self.roster.surfaces.get(&cell.id) {
                admit_surface(surface, thermal, self.catalog, self.controls)?;
                surface.binding.liquid_enthalpy_offset_j_kg
            } else {
                0.0
            };
            self.births.insert(
                cell.id,
                Birth {
                    mass_kg: cell.mass_kg,
                    liquid_enthalpy_offset_j_kg: offset,
                },
            );
        }
        Ok(())
    }

    fn step(
        &mut self,
        _: f64,
        _: f64,
        _: &BTreeMap<CellId, f64>,
        _: &DepositionSimulation<'_>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), Self::Error> {
        checkpoint()?;
        Ok(())
    }
}

fn admit_surface(
    surface: ImplicitGasSurface<'_>,
    thermal: &DepositionSimulation<'_>,
    catalog: &ThermochemicalGasModel,
    controls: ImplicitGasControls,
) -> Result<(), AdaptiveImplicitGasError> {
    let cell = thermal
        .cells()
        .get(&surface.cell_id)
        .ok_or(DepositionError::InvalidInput(
            "selected liquid cell must be active",
        ))?;
    if thermal.material_storage_basis(cell.material_id)? != ThermalStorageBasis::Enthalpy {
        return Err(DepositionError::InvalidInput(
            "finite-gas evaporation requires original liquid enthalpy storage",
        )
        .into());
    }
    let [lo, hi] = surface.pure_liquid_enthalpy_range_j_kg;
    if !(lo..=hi).contains(&cell.specific_h_j_kg) {
        return Err(DepositionError::MaterialDomain.into());
    }
    admit_liquid_source(
        surface.law,
        surface.binding,
        catalog,
        thermal.temperature_k(surface.cell_id)?,
        cell.specific_h_j_kg,
        controls.phase,
    )?;
    Ok(())
}

fn admit_gas(
    gas: &ClosedGasParcel,
    catalog: &ThermochemicalGasModel,
    controls: ImplicitGasControls,
) -> Result<(), AdaptiveImplicitGasError> {
    admit_catalog(gas, catalog)?;
    let state = gas.state().state();
    if !positive(gas.mass_kg() * state.specific_heat_cp / state.gamma)
        || !(controls.phase.gas.gas.minimum_temperature_k
            ..=controls.phase.gas.gas.maximum_temperature_k)
            .contains(&state.temperature)
    {
        return Err(ImplicitGasError::Invalid(
            "initial gas must lie inside its admitted caloric bracket",
        )
        .into());
    }
    Ok(())
}
