//! Liquid-vapor exchange from an explicitly resolved solid/liquid inventory.
//!
//! The original equilibrium chart owns bulk sensible and fusion storage. The
//! original film law owns LIQUID-to-vapor latent heat. Its caloric compatibility
//! therefore uses pure liquid h at the melting plateau, not bulk mixture h.
//! The shared extensive update retains their difference exactly once: cooling
//! can freeze liquid while T stays at melting, and condensation can melt solid.
//!
//! Only liquid present at the start of the vapor substep can evaporate. The
//! optional heat-first split supplies its accepted post-heat inventory; liquid
//! formed later within vapor transfer cannot finance that frozen rate.
//! The caller supplies exposed LIQUID area;
//! phase mass fraction does not determine surface coverage. A remaining liquid
//! seed is required by the default stepping path. The opt-in [`boundary`]
//! owner can return exact freezing completion and unconsumed physical time;
//! selecting the subsequent interface and transport regime remains explicit.
//!
//! This is a sharp, pure-substance, local-equilibrium reduction on the original
//! chart. Its liquid sensible branch remains available after complete melting.
//! Mushy alloys, separate phase temperatures, moving geometry and inferred
//! pressure-dependent phase boundaries are outside this reduction. The parent
//! film's frozen-rate, fixed-gas-volume and caloric-account assumptions apply.

/// Conservative completion of the liquid interface at its original fusion knot.
pub mod boundary;
/// Reciprocal gas-film sensible heat on the same original fusion inventory.
pub mod heat;

use super::{
    ClosedGasParcel, EvaporationLaw, GasEvaporationCandidate, GasEvaporationError,
    GasEvaporationSolve, GasLiquidBinding, LiquidCaloricState, LiquidSurfaceState,
    ThermochemicalGasModel, finish_surface, nonnegative, positive, prepare_surface_with_caloric,
};
use crate::phase::{EquilibriumEnthalpyPhaseCurve, EquilibriumPhaseState};
use crate::sublimation::gas::chamber::mixed::SharpMeltingCaloric;

fn invalid(message: &'static str) -> GasEvaporationError {
    GasEvaporationError::Invalid(message)
}

/// Borrow the original sharp fusion chart without copying or refitting its knots.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct EquilibriumEvaporation<'a> {
    source: SharpMeltingCaloric<'a>,
    pub(in crate::evaporation::gas) solid_at_melting_h_j_kg: f64,
    liquid_at_melting_h_j_kg: f64,
    liquid_at_melting_density_kg_m3: f64,
}

/// Both constituent masses resolved from the same original bulk enthalpy.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PhaseInventoryState {
    /// Original source identity, temperature, density and equilibrium fractions.
    pub equilibrium: EquilibriumPhaseState,
    /// Total condensed mass, including solid that cannot evaporate as liquid.
    pub mass_kg: f64,
    /// Actual solid inventory [kg].
    pub solid_mass_kg: f64,
    /// Actual liquid inventory [kg].
    pub liquid_mass_kg: f64,
}

/// Fusion and vapor transfer are distinct parts of one actual mass/energy update.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PhaseInventoryStep {
    /// Accepted initial condensed partition.
    pub initial: PhaseInventoryState,
    /// Proposed partition at the same endpoint as the original gas exchange.
    pub final_state: PhaseInventoryState,
    /// Positive liquid -> solid; negative melting. This is not vapor condensation.
    pub freezing_mass_kg: f64,
    /// Signed fusion release [J], already contained in bulk enthalpy. Never add
    /// this diagnostic as another heat source.
    pub fusion_energy_j: f64,
    /// Delta liquid mass + mass into vapor + freezing mass [kg].
    pub liquid_partition_residual_kg: f64,
}

impl<'a> EquilibriumEvaporation<'a> {
    /// Admit the same sharp fusion source used by the solid-interface owner,
    /// with a strictly increasing, all-liquid sensible branch above melting.
    pub fn new(curve: &'a EquilibriumEnthalpyPhaseCurve) -> Result<Self, GasEvaporationError> {
        let source = SharpMeltingCaloric::new(curve).map_err(|_| {
            invalid("evaporation requires an original sharp solid/liquid fusion chart")
        })?;
        let knots = curve.knots();
        let end = knots
            .iter()
            .position(|knot| knot.liquid_mass_fraction == 1.0)
            .ok_or(invalid(
                "sharp fusion requires its original liquid endpoint",
            ))?;
        let solid_at_melting_h_j_kg = knots[..end]
            .iter()
            .rfind(|knot| knot.liquid_mass_fraction == 0.0)
            .ok_or(invalid("sharp fusion requires its original solid endpoint"))?
            .specific_enthalpy_j_kg;
        for pair in knots[end..].windows(2) {
            if pair[0].liquid_mass_fraction != 1.0
                || pair[1].liquid_mass_fraction != 1.0
                || pair[0].temperature_k >= pair[1].temperature_k
            {
                return Err(invalid(
                    "liquid sensible branch must be single phase and strictly increasing",
                ));
            }
        }
        Ok(Self {
            source,
            solid_at_melting_h_j_kg,
            liquid_at_melting_h_j_kg: knots[end].specific_enthalpy_j_kg,
            liquid_at_melting_density_kg_m3: knots[end].bulk_density_kg_m3,
        })
    }

    /// The unmodified material-bound bulk enthalpy source.
    #[must_use]
    pub const fn curve(self) -> &'a EquilibriumEnthalpyPhaseCurve {
        self.source.curve()
    }

    /// Original plateau energy difference, separate from vaporization latent heat.
    #[must_use]
    pub fn fusion_latent_heat_j_kg(self) -> f64 {
        self.source.fusion_latent_heat_j_kg()
    }

    // Numeric original-source evaluation avoids constructing/hashing a phase
    // identity during the shared film preparation. Reports retain full states.
    pub(in crate::evaporation::gas) fn caloric(
        self,
        h: f64,
    ) -> Result<LiquidCaloricState, GasEvaporationError> {
        self.caloric_with_completion(h, false)
    }

    // Only an event-owned endpoint may use the exact original solid fusion
    // knot; the liquid constituent still supplies the limiting vaporization h.
    pub(in crate::evaporation::gas) fn caloric_with_completion(
        self,
        h: f64,
        allow_completed: bool,
    ) -> Result<LiquidCaloricState, GasEvaporationError> {
        let state = self.curve().values_at_specific_enthalpy(h).map_err(|_| {
            invalid("original solid/liquid enthalpy chart refused the evaporation state")
        })?;
        if state.liquid_mass_fraction == 0.0
            && !(allow_completed && h == self.solid_at_melting_h_j_kg)
        {
            return Err(invalid(
                "no liquid interface remains; freezing completion requires its own event owner",
            ));
        }
        let fully_liquid = state.liquid_mass_fraction == 1.0;
        Ok(LiquidCaloricState {
            temperature_k: state.temperature_k,
            liquid_enthalpy_j_kg: if fully_liquid {
                h
            } else {
                self.liquid_at_melting_h_j_kg
            },
            liquid_mass_fraction: state.liquid_mass_fraction,
            liquid_density_kg_m3: Some(if fully_liquid {
                state.bulk_density_kg_m3
            } else {
                self.liquid_at_melting_density_kg_m3
            }),
        })
    }

    /// Resolve the original phase partition with a remaining liquid seed.
    /// `bulk.mass_kg` and h describe the entire condensed body; area is exposed
    /// liquid area, not total body area or a mass-fraction-derived estimate.
    pub fn state(
        self,
        bulk: LiquidSurfaceState,
    ) -> Result<PhaseInventoryState, GasEvaporationError> {
        self.state_with_completion(bulk, false)
    }

    pub(in crate::evaporation::gas) fn state_with_completion(
        self,
        bulk: LiquidSurfaceState,
        allow_completed: bool,
    ) -> Result<PhaseInventoryState, GasEvaporationError> {
        if !positive(bulk.mass_kg) || !nonnegative(bulk.area_m2) {
            return Err(invalid(
                "positive total condensed mass and nonnegative exposed liquid area required",
            ));
        }
        let equilibrium = self
            .curve()
            .state_at_specific_enthalpy(bulk.specific_enthalpy_j_kg)
            .map_err(|_| invalid("original solid/liquid caloric state refused"))?;
        let solid_mass_kg = bulk.mass_kg * equilibrium.solid_mass_fraction();
        let liquid_mass_kg = bulk.mass_kg * equilibrium.liquid_mass_fraction();
        let completed = allow_completed
            && bulk.specific_enthalpy_j_kg == self.solid_at_melting_h_j_kg
            && equilibrium.is_fully_solid();
        if !nonnegative(liquid_mass_kg)
            || (liquid_mass_kg == 0.0 && !completed)
            || (equilibrium.liquid_mass_fraction() > 0.0 && liquid_mass_kg == 0.0)
            || !nonnegative(solid_mass_kg)
            || (equilibrium.solid_mass_fraction() > 0.0 && solid_mass_kg == 0.0)
        {
            return Err(invalid(
                "liquid phase is absent or a condensed constituent inventory is unrepresentable",
            ));
        }
        Ok(PhaseInventoryState {
            equilibrium,
            mass_kg: bulk.mass_kg,
            solid_mass_kg,
            liquid_mass_kg,
        })
    }

    /// Reconstruct fusion from actual endpoints, independently of the film rate.
    pub fn partition_step(
        self,
        initial: LiquidSurfaceState,
        final_bulk: LiquidSurfaceState,
        transferred_mass_kg: f64,
        mass_tolerance_kg: f64,
    ) -> Result<PhaseInventoryStep, GasEvaporationError> {
        self.partition_step_with_completion(
            initial,
            final_bulk,
            transferred_mass_kg,
            mass_tolerance_kg,
            false,
        )
    }

    pub(in crate::evaporation::gas) fn partition_step_with_completion(
        self,
        initial: LiquidSurfaceState,
        final_bulk: LiquidSurfaceState,
        transferred_mass_kg: f64,
        mass_tolerance_kg: f64,
        allow_completed: bool,
    ) -> Result<PhaseInventoryStep, GasEvaporationError> {
        if !positive(mass_tolerance_kg) || !transferred_mass_kg.is_finite() {
            return Err(invalid(
                "finite transfer and positive phase-partition mass tolerance required",
            ));
        }
        let initial = self.state(initial)?;
        let final_state = self.state_with_completion(final_bulk, allow_completed)?;
        let freezing_mass_kg = final_state.solid_mass_kg - initial.solid_mass_kg;
        let residual = (final_state.liquid_mass_kg - initial.liquid_mass_kg)
            + transferred_mass_kg
            + freezing_mass_kg;
        let total_residual = (final_state.mass_kg - initial.mass_kg) + transferred_mass_kg;
        let fusion_energy_j = freezing_mass_kg * self.fusion_latent_heat_j_kg();
        if ![residual, total_residual, fusion_energy_j]
            .into_iter()
            .all(f64::is_finite)
            || residual.abs().max(total_residual.abs()) > mass_tolerance_kg
        {
            return Err(invalid(
                "actual liquid/solid/vapor inventories failed their mass balance",
            ));
        }
        Ok(PhaseInventoryStep {
            initial,
            final_state,
            freezing_mass_kg,
            fusion_energy_j,
            liquid_partition_residual_kg: residual,
        })
    }
}

/// One original gas endpoint and the same condensed body's fusion inventory.
#[derive(Clone, Debug, PartialEq)]
pub struct EquilibriumSurfaceCandidate<'curve, 'catalog> {
    inventory: EquilibriumEvaporation<'curve>,
    exchange: GasEvaporationCandidate<'catalog>,
    phase: PhaseInventoryStep,
}
impl<'curve, 'catalog> EquilibriumSurfaceCandidate<'curve, 'catalog> {
    /// Underlying `liquid()` and `liquid_enthalpy_change_j` describe the WHOLE
    /// condensed body here. Recession is a vapor-transfer-equivalent thickness;
    /// it does not represent geometry changes caused by freezing or melting.
    #[must_use]
    pub const fn exchange(&self) -> &GasEvaporationCandidate<'catalog> {
        &self.exchange
    }
    /// Actual original-source phase inventory and its fusion subaccount.
    #[must_use]
    pub const fn phase(&self) -> &PhaseInventoryStep {
        &self.phase
    }
    /// Retained binding to the original condensed caloric source.
    #[must_use]
    pub const fn inventory(&self) -> EquilibriumEvaporation<'curve> {
        self.inventory
    }
    /// Consume the mutually consistent exchange and phase proposals together.
    #[must_use]
    pub fn into_parts(self) -> (GasEvaporationCandidate<'catalog>, PhaseInventoryStep) {
        (self.exchange, self.phase)
    }
}

/// Evaporate/condense through an authored liquid surface of a melting body.
///
/// The shared original film and finite gas own rates, species, carried energy,
/// inventory limits and endpoint checks. Bulk h owns fusion; no second fusion
/// heat is supplied. Condensation carries the original upstream vapor enthalpy.
/// Both caller-owned inputs stay unchanged on every failure and cancellation.
#[allow(clippy::too_many_arguments)]
pub fn advance_equilibrium_surface<'curve, 'catalog, L: EvaporationLaw + ?Sized>(
    inventory: EquilibriumEvaporation<'curve>,
    law: &L,
    binding: &GasLiquidBinding,
    bulk: LiquidSurfaceState,
    gas: &ClosedGasParcel,
    catalog: &'catalog ThermochemicalGasModel,
    mass_transfer_coefficient_m_s: f64,
    duration_s: f64,
    controls: GasEvaporationSolve,
    mut cancelled: impl FnMut() -> bool,
) -> Result<EquilibriumSurfaceCandidate<'curve, 'catalog>, GasEvaporationError> {
    let prepared = prepare_surface_with_caloric(
        law,
        binding,
        bulk,
        gas,
        catalog,
        |h| inventory.caloric(h),
        mass_transfer_coefficient_m_s,
        duration_s,
        controls,
        Ok,
        &mut cancelled,
    )?;
    let phase = inventory.partition_step(
        bulk,
        prepared.liquid,
        prepared.dm,
        controls.gas.mass_tolerance_kg,
    )?;
    let exchange = finish_surface(
        prepared,
        law,
        binding,
        gas,
        catalog,
        mass_transfer_coefficient_m_s,
        duration_s,
        controls,
        cancelled,
    )?;
    Ok(EquilibriumSurfaceCandidate {
        inventory,
        exchange,
        phase,
    })
}

#[cfg(test)]
pub(in crate::evaporation::gas) mod tests;
