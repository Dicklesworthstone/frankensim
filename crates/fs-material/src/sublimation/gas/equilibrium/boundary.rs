//! Conservative completion of melting on the original molecular interface.
//!
//! The original frozen streams give `m1=m0-dm` and
//! `m1*(h1-h0)=dm*(h0-h_s-L_sub-alignment+R*T_s/2)-m_in*delta_e-Q_g`.
//! Replacing h0 by the original liquid fusion knot h_l on the right gives the
//! event denominator. Its ratio to `m0*(h_l-h0)` selects the same streams' time
//! of complete melting, including deposition and reflected heat together.
//!
//! A bounded probe exposure locates that event before an oversized requested
//! interval or its gross withdrawals are evaluated. The original owners then
//! evaluate the accepted streams and reflection, check their initial donor
//! inventories, and reconstruct the original gas once. The exact source knot
//! is accepted only within the original represented mass and energy budgets.
//!
//! Completion returns the whole molten body and unconsumed time. The caller
//! selects its next liquid interface, exposed area and transport regime; a
//! continuum evaporation law is never inferred from molecular sublimation.
//! Retain the same source chart, gas species and caloric-datum offset when
//! handing off this body. Fusion is already in bulk enthalpy, not extra heat.
//! These event times belong to the explicit frozen-rate discretization and
//! require timestep refinement to resolve evolving gas and surface conditions.
//! If external heat is split over the interval, retry the complete coupled
//! interval at an earlier event; do not shorten the clock of an already
//! accepted full-duration heat substep.

use super::super::{
    SurfaceSnapshot, condensed_thermal_target, molecular_stream_energy, nonnegative,
};
use super::{
    ClosedGasParcel, EquilibriumSublimation, EquilibriumSurfaceCandidate, GasChamberSurface,
    GasSurfaceBinding, GasSurfaceCandidate, GasSurfaceError, GasSurfaceSolve, PhaseInventoryStep,
    ReflectedHeatStep, SolidSurfaceState, SublimationError, SublimationLaw, ThermalAccommodation,
    ThermochemicalGasModel, admit_surface_inventory, finish_equilibrium_surface, invalid, mixed,
    poll, positive, species_snapshot, validate_controls,
};

/// Accepted molecular time and the terminal solid-interface decision.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MeltingBoundaryStep {
    /// Physical duration consumed by both molecular streams and reflection [s].
    pub elapsed_duration_s: f64,
    /// Requested duration not consumed before completion [s].
    pub remaining_duration_s: f64,
    /// Exactly zero solid at the original fully liquid fusion knot.
    pub melting_completed: bool,
    /// Represented condensed storage minus the original transfer target [J].
    /// Included in exchange().step().energy_residual_j, never added as heat.
    pub condensed_energy_residual_j: f64,
}

/// Original scalar exchange and phase inventories at the same accepted time.
#[derive(Clone, Debug, PartialEq)]
pub struct EquilibriumMeltingCandidate<'curve, 'catalog> {
    endpoint: EquilibriumSurfaceCandidate<'curve, 'catalog>,
    boundary: MeltingBoundaryStep,
}
impl<'curve, 'catalog> EquilibriumMeltingCandidate<'curve, 'catalog> {
    /// Shared original-owner endpoint; complete melting is terminal for its
    /// solid interface. `exchange().solid()` carries the TOTAL condensed body.
    #[must_use]
    pub const fn endpoint(&self) -> &EquilibriumSurfaceCandidate<'curve, 'catalog> {
        &self.endpoint
    }
    #[must_use]
    pub const fn exchange(&self) -> &GasSurfaceCandidate<'catalog> {
        self.endpoint.exchange()
    }
    #[must_use]
    pub const fn phase(&self) -> &PhaseInventoryStep {
        self.endpoint.phase()
    }
    #[must_use]
    pub const fn reflection(&self) -> Option<&ReflectedHeatStep> {
        self.endpoint.reflection()
    }
    /// The same admitted source; its ordinary state method still requires solid.
    #[must_use]
    pub const fn inventory(&self) -> EquilibriumSublimation<'curve> {
        self.endpoint.inventory()
    }
    #[must_use]
    pub const fn boundary(&self) -> &MeltingBoundaryStep {
        &self.boundary
    }
    /// Consume the consistent condensed/gas endpoint and its handoff time.
    #[must_use]
    pub fn into_parts(
        self,
    ) -> (
        EquilibriumSurfaceCandidate<'curve, 'catalog>,
        MeltingBoundaryStep,
    ) {
        (self.endpoint, self.boundary)
    }
}

/// Advance at most the requested duration, stopping at complete melting.
///
/// An initial solid seed is required. Without completion, the accepted endpoint
/// agrees with the existing scalar interval. An event clips BOTH mass streams
/// and every reflected species' heat before the original gross donor gates.
/// Earlier donor exhaustion, unsupported sources or an unrepresentable event
/// still refuse atomically; the event cannot relax those physical budgets.
#[allow(clippy::too_many_arguments)]
pub fn advance_to_melting<'curve, 'catalog>(
    inventory: EquilibriumSublimation<'curve>,
    law: &dyn SublimationLaw,
    binding: &GasSurfaceBinding,
    bulk: SolidSurfaceState,
    gas: &ClosedGasParcel,
    catalog: &'catalog ThermochemicalGasModel,
    accommodation: Option<&ThermalAccommodation>,
    maximum_duration_s: f64,
    controls: GasSurfaceSolve,
    mut cancelled: impl FnMut() -> bool,
) -> Result<EquilibriumMeltingCandidate<'curve, 'catalog>, GasSurfaceError> {
    poll(&mut cancelled)?;
    validate_controls(controls, maximum_duration_s)?;
    inventory.state(bulk)?;
    let caloric = |h| {
        inventory
            .curve()
            .temperature_at_specific_enthalpy(h)
            .map_err(|_| {
                SublimationError::Invalid("original solid/liquid caloric chart refused state")
            })
    };
    let surface = GasChamberSurface {
        id: 0,
        law,
        binding,
        solid: bulk,
        temperature_from_enthalpy: &caloric,
    };
    let gas_snapshot = species_snapshot(gas, catalog, &binding.species)?;
    let mut evaluate =
        |dt| -> Result<(SurfaceSnapshot, Option<ReflectedHeatStep>), GasSurfaceError> {
            poll(&mut cancelled)?;
            let snapshot = mixed::snapshot(
                &surface,
                inventory.source,
                gas_snapshot,
                catalog,
                dt,
                controls,
            )?;
            let reflection = match accommodation {
                Some(law) => Some(law.evaluate(
                    binding,
                    gas,
                    snapshot.temperature_k,
                    bulk.area_m2,
                    dt,
                    &mut cancelled,
                )?),
                None => None,
            };
            Ok((snapshot, reflection))
        };
    // Bound geometric exposure, not time: a tiny area can require a long but
    // representable interval. This also admits an earlier finite event when
    // area * requested_duration would overflow before any donor check.
    let probe_duration_s = if bulk.area_m2 > 0.0 {
        maximum_duration_s.min(1.0 / bulk.area_m2)
    } else {
        maximum_duration_s
    };
    let (mut snapshot, mut reflection) = evaluate(probe_duration_s)?;
    let reflected_heat = reflection.as_ref().map_or(0.0, |row| row.heat_into_gas_j);
    let energy = molecular_stream_energy(snapshot, reflected_heat)?;
    let boundary_h = inventory.source.liquid_at_melting_h_j_kg();
    let drive =
        condensed_thermal_target(snapshot, boundary_h, energy.difference_j_kg, reflected_heat);
    if !drive.is_finite() {
        return Err(invalid("terminal melting balance is unrepresentable"));
    }
    let mut elapsed_duration_s = maximum_duration_s;
    let mut terminal_enthalpy_j_kg = None;
    if drive > 0.0 {
        let headroom_j = bulk.mass_kg * (boundary_h - bulk.specific_enthalpy_j_kg);
        if !positive(headroom_j) {
            return Err(invalid("terminal melting headroom is unrepresentable"));
        }
        // Retain the ordinary ratio ordering whenever representable. Alternate
        // orders rescue a finite event hidden by an intermediate overflow or
        // underflow; the accepted endpoint still measures the extensive defect.
        let times = [
            probe_duration_s * (headroom_j / drive),
            headroom_j * (probe_duration_s / drive),
            (probe_duration_s * headroom_j) / drive,
        ];
        let event_duration_s = times.into_iter().find(|t| positive(*t)).unwrap_or(times[0]);
        if event_duration_s == 0.0 || event_duration_s.is_nan() {
            return Err(invalid("terminal melting duration is unrepresentable"));
        }
        if event_duration_s <= maximum_duration_s {
            elapsed_duration_s = event_duration_s;
            terminal_enthalpy_j_kg = Some(boundary_h);
        }
    }
    if elapsed_duration_s != probe_duration_s {
        (snapshot, reflection) = evaluate(elapsed_duration_s)?;
    }
    admit_surface_inventory(snapshot, elapsed_duration_s, controls)?;
    let prepared = mixed::prepare_with_boundary(
        &surface,
        inventory.source,
        snapshot,
        catalog,
        controls,
        reflection.as_ref().map_or(0.0, |row| row.heat_into_gas_j),
        terminal_enthalpy_j_kg,
    )?;
    let remaining_duration_s = maximum_duration_s - elapsed_duration_s;
    if !nonnegative(remaining_duration_s) {
        return Err(invalid(
            "terminal melting cannot represent its remaining duration",
        ));
    }
    // A root can round just beyond the requested time while the ORIGINAL
    // represented update already lands exactly on h_l. That is completion too;
    // it introduces neither a temperature epsilon nor an extra energy change.
    let melting_completed =
        terminal_enthalpy_j_kg.is_some() || prepared.solid.specific_enthalpy_j_kg == boundary_h;
    let boundary = MeltingBoundaryStep {
        elapsed_duration_s,
        remaining_duration_s,
        melting_completed,
        condensed_energy_residual_j: prepared.energy_residual_j,
    };
    let endpoint = finish_equilibrium_surface(
        inventory,
        binding,
        bulk,
        gas,
        catalog,
        gas_snapshot,
        prepared,
        reflection,
        accommodation,
        controls,
        melting_completed,
        &mut cancelled,
    )?;
    if melting_completed {
        let absolute_residual = boundary.condensed_energy_residual_j.abs()
            + endpoint
                .exchange()
                .gas_exchange()
                .step()
                .energy_residual_j
                .abs();
        if !absolute_residual.is_finite() || absolute_residual > controls.gas.energy_tolerance_j {
            return Err(GasSurfaceError::EnergyBalance {
                residual_j: absolute_residual,
                tolerance_j: controls.gas.energy_tolerance_j,
            });
        }
    }
    poll(&mut cancelled)?;
    Ok(EquilibriumMeltingCandidate { endpoint, boundary })
}

#[cfg(test)]
mod tests;
