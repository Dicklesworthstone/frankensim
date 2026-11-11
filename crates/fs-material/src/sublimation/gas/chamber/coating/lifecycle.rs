//! Molecular turnover through supported-film formation and complete extinction.
//!
//! The original thermal turnover owner integrates finite-vapor backpressure at
//! frozen support and gas temperatures. The SAME integrated opposing streams
//! then enter the original supported-storage and source-relative gas balances.
//! Gas composition, pressure, transport and caloric state are reconstructed once.
//! The nonvolatile support never supplies vapor, even when the film vanishes.
//!
//! This opt-in model admits immediate heterogeneous deposition. At zero film
//! inventory, an undersaturated interface re-emits sticking incoming molecules
//! at the support temperature without storing a negative film. Supersaturation
//! starts a positive film with concurrent emission. No numerical mass seed,
//! nucleation barrier or finite adsorption residence time is inferred.
//!
//! Only one exchanging supported surface is admitted: independent analytical
//! reservoirs would double-spend gas shared by several surfaces. Both body and
//! gas temperatures evolve after the frozen-rate interval, so its extinction
//! time is a local first-order event estimate, not an exact nonisothermal event.
//! Nonsticking molecules and carriers reflect specularly. Optional dilute film
//! resistance has the original prescribed closure, without film energy storage.

use super::super::{
    GasChamberCandidate, GasChamberSolve, InventoryAdmission, advance_chamber_with_inventory,
};
use super::{
    GasSupportedChamberCandidate, GasSupportedChamberSolve, GasSupportedSurface,
    GasSupportedSurfaceStep, SupportedSurfaceState, prepare, snapshot_with_film, validate_storage,
};
use crate::sublimation::gas::{
    ClosedGasParcel, GasSurfaceError, GasSurfaceSolve, ThermochemicalGasModel, invalid,
    nonnegative, phase_at, poll, positive, validate_controls,
};
use crate::thermochemical_gas::energy::parcel::transfer::species::boundary::{
    GasSpeciesBoundarySolve, admit_catalog,
};
use fs_thermal::sublimation::molecular::turnover::{
    MolecularTurnover, integrate_molecular_turnover,
};
use fs_thermal::sublimation::{SolverConfig, SurfaceKinetics, VaporReservoir};

/// Original caloric solvers and explicit physical admission budgets.
/// Gross turnover is integrated, so no old-inventory fraction is supplied.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct GasSupportedLifecycleSolve {
    pub gas: GasSpeciesBoundarySolve,
    pub storage: SolverConfig,
    pub latent_compatibility_tolerance_j_kg: f64,
    /// Deposit volume / fixed accessible gas volume, in (0, 0.01].
    /// The gas volume excludes fixed supports; the cap bounds omitted deposit
    /// displacement, not pressure-work error or the thin-film approximation.
    pub maximum_deposit_volume_fraction: f64,
}
impl GasSupportedLifecycleSolve {
    pub(crate) fn validate(self) -> Result<(), GasSurfaceError> {
        validate_controls(self.chamber_controls().chamber.surface, 0.0)?;
        validate_storage(self.storage)?;
        if !positive(self.latent_compatibility_tolerance_j_kg)
            || !positive(self.maximum_deposit_volume_fraction)
            || self.maximum_deposit_volume_fraction > 0.01
        {
            return Err(invalid(
                "invalid supported lifecycle source or volume budget",
            ));
        }
        Ok(())
    }

    // The existing caloric/storage owners receive their original controls.
    // The legacy fraction remains a valid internal value for shared admission;
    // the explicit IntegratedSingleSurface policy, never this fraction, admits
    // the analytically solved turnover. It is not a user-facing physical knob.
    pub(crate) fn chamber_controls(self) -> GasSupportedChamberSolve {
        GasSupportedChamberSolve {
            chamber: GasChamberSolve {
                surface: GasSurfaceSolve {
                    gas: self.gas,
                    maximum_gross_inventory_fraction: 0.5,
                    latent_compatibility_tolerance_j_kg: self.latent_compatibility_tolerance_j_kg,
                },
                maximum_surfaces: 1,
            },
            storage: self.storage,
        }
    }
}

/// One complete original gas/support transaction and its kinetic event receipt.
/// All inputs remain unchanged on failure or cancellation.
#[derive(Clone, Debug, PartialEq)]
pub struct GasSupportedLifecycleCandidate<'catalog> {
    endpoint: GasSupportedChamberCandidate<'catalog>,
    turnover: MolecularTurnover,
    deposit_volume_fraction: [f64; 2],
}
impl<'catalog> GasSupportedLifecycleCandidate<'catalog> {
    #[must_use]
    pub const fn endpoint(&self) -> &GasSupportedChamberCandidate<'catalog> {
        &self.endpoint
    }
    #[must_use]
    pub const fn gas(&self) -> &ClosedGasParcel {
        self.endpoint.gas()
    }
    #[must_use]
    pub const fn catalog(&self) -> &'catalog ThermochemicalGasModel {
        self.endpoint.catalog()
    }
    #[must_use]
    pub fn state(&self) -> SupportedSurfaceState {
        self.endpoint.surfaces()[0].state
    }
    /// The reservoir here is the frozen-temperature kinetic subproblem, not
    /// the final thermochemical gas. Use `gas()` for the accepted physical gas.
    #[must_use]
    pub const fn turnover(&self) -> &MolecularTurnover {
        &self.turnover
    }
    #[must_use]
    pub const fn deposit_volume_fraction(&self) -> [f64; 2] {
        self.deposit_volume_fraction
    }
    #[must_use]
    pub fn into_endpoint(self) -> GasSupportedChamberCandidate<'catalog> {
        self.endpoint
    }
    #[must_use]
    pub fn into_parts(
        self,
    ) -> (
        GasSupportedChamberCandidate<'catalog>,
        MolecularTurnover,
        [f64; 2],
    ) {
        (self.endpoint, self.turnover, self.deposit_volume_fraction)
    }
}

fn volume_fraction(
    mass_kg: f64,
    density_kg_m3: f64,
    gas_volume_m3: f64,
    maximum: f64,
) -> Result<f64, GasSurfaceError> {
    let fraction = (mass_kg / density_kg_m3) / gas_volume_m3;
    if !nonnegative(fraction) || fraction > maximum {
        return Err(invalid(
            "supported film exceeds the declared fixed-gas-volume approximation",
        ));
    }
    Ok(fraction)
}

/// Whether the represented net transfer keeps the inventory class of the
/// integrated tail. An inventory is negative, exactly zero or positive; a
/// non-finite one is never encoded. The tail is nonnegative, so a positive tail
/// whose old inventory plus net transfer rounds to zero OR BELOW is a positive
/// tail that the net-transfer representation cannot carry, and an exactly zero
/// tail needs an exactly zero sum. Both kinds refuse, whichever way the last bit
/// of the sum rounds.
fn net_transfer_encodes_tail(sum: f64, tail: f64) -> bool {
    fn class(value: f64) -> Option<u8> {
        if !value.is_finite() {
            None
        } else if value > 0.0 {
            Some(2)
        } else if value == 0.0 {
            Some(1)
        } else {
            Some(0)
        }
    }
    matches!((class(sum), class(tail)), (Some(a), Some(b)) if a == b)
}

/// Advance one internally isothermal supported deposit through birth/extinction.
///
/// `film` retains the original linear series resistance: `None` is the local
/// free-molecular limit, `Some(0.0)` seals molecular exchange, and a positive
/// value is the caller's admitted coefficient in m/s. It is frozen together
/// with area and both initial temperatures during the analytical mass interval.
/// The full original caloric laws then receive the actual opposing molecular
/// energies, including heat exchanged by equal gross streams on a bare support.
/// No latent source is added again and no gas heat capacity is substituted.
///
/// Mass may be emitted, mixed and recondensed within this interval. Exact zero
/// coating and exact zero individual gas species are valid; total gas vacuum
/// remains the original gas owner's refusal. Both endpoint source compatibility
/// and every molecular/net inventory residual use the supplied absolute budgets.
/// A zero duration performs full source/geometry admission and preserves state.
/// For evolving temperatures, refine the physical time interval; the analytic
/// mass substep alone does not establish coupled temporal accuracy or stability.
#[allow(clippy::too_many_arguments)]
pub fn advance_supported_lifecycle<'catalog>(
    surface: &GasSupportedSurface<'_>,
    film: Option<f64>,
    gas: &ClosedGasParcel,
    catalog: &'catalog ThermochemicalGasModel,
    duration_s: f64,
    controls: GasSupportedLifecycleSolve,
    mut cancelled: impl FnMut() -> bool,
) -> Result<GasSupportedLifecycleCandidate<'catalog>, GasSurfaceError> {
    poll(&mut cancelled)?;
    controls.validate()?;
    if !nonnegative(duration_s) || film.is_some_and(|value| !nonnegative(value)) {
        return Err(invalid(
            "invalid supported lifecycle duration or gas-film coefficient",
        ));
    }
    // Source identity precedes phase rates and heat solves, including at dt=0.
    admit_catalog(gas, catalog)?;
    let chamber_controls = controls.chamber_controls();
    let mut integrated = None;
    let mut storage = None;
    let mut initial_fraction = 0.0;
    let mut final_fraction = 0.0;
    let (candidate, _) = advance_chamber_with_inventory(
        std::slice::from_ref(surface),
        None,
        gas,
        catalog,
        duration_s,
        chamber_controls.chamber,
        InventoryAdmission::IntegratedSingleSurface,
        || Ok(()),
        |surface, context| {
            // Zero exposure retains the exact original caloric/kinetic source
            // admission without forming an explicit gross-withdrawal proposal.
            let mut snapshot =
                snapshot_with_film(surface, context, catalog, film, 0.0, chamber_controls)?;
            let initial_phase = phase_at(surface.law, snapshot.temperature_k, context.molar_mass)?;
            initial_fraction = volume_fraction(
                surface.state.deposit_mass_kg,
                initial_phase.solid_density_kg_m3,
                gas.volume_m3(),
                controls.maximum_deposit_volume_fraction,
            )?;
            let turnover = integrate_molecular_turnover(
                surface.law,
                snapshot.temperature_k,
                SurfaceKinetics {
                    accommodation_coefficient: surface.binding.accommodation_coefficient,
                    mass_transfer_coefficient_m_s: film,
                },
                snapshot.solid,
                VaporReservoir::Finite {
                    vapor_mass_kg: context.vapor_mass_kg,
                    volume_m3: gas.volume_m3(),
                    temperature_k: context.temperature_k,
                },
                duration_s,
            )?;
            for residual in [
                turnover.stream_mass_residual_kg,
                turnover.deposit_mass_residual_kg,
                turnover.vapor_mass_residual_kg,
                (snapshot.solid.mass_kg - turnover.transferred_mass_kg) - turnover.deposit_mass_kg,
            ] {
                if !residual.is_finite() || residual.abs() > controls.gas.mass_tolerance_kg {
                    return Err(invalid(
                        "integrated supported turnover exceeds its mass budget",
                    ));
                }
            }
            let VaporReservoir::Finite { vapor_mass_kg, .. } = turnover.reservoir else {
                return Err(invalid(
                    "integrated supported turnover lost its finite reservoir",
                ));
            };
            if !net_transfer_encodes_tail(
                snapshot.solid.mass_kg - turnover.transferred_mass_kg,
                turnover.deposit_mass_kg,
            ) || !net_transfer_encodes_tail(
                context.vapor_mass_kg + turnover.transferred_mass_kg,
                vapor_mass_kg,
            ) {
                return Err(invalid(
                    "integrated molecular net transfer cannot represent a positive inventory tail",
                ));
            }
            snapshot.outgoing_mass_kg = turnover.outgoing_mass_kg;
            snapshot.incoming_mass_kg = turnover.incoming_mass_kg;
            snapshot.transferred_mass_kg = turnover.transferred_mass_kg;
            snapshot.incoming_specific_energy_j_kg = if turnover.incoming_mass_kg == 0.0 {
                None
            } else {
                // An initially absent species can be emitted and return in the
                // SAME interval. Its energy is still the retained original law
                // at the frozen gas temperature, never a guessed latent value.
                let energy = catalog.species_specific_enthalpy_j_kg(
                    &surface.binding.species,
                    context.temperature_k,
                )? - 0.5 * context.specific_gas_constant * context.temperature_k;
                if !energy.is_finite() {
                    return Err(invalid(
                        "integrated incoming molecular energy is unrepresentable",
                    ));
                }
                Some(energy)
            };
            snapshot.partial_enthalpy_change_j_kg =
                if turnover.outgoing_mass_kg == 0.0 && turnover.incoming_mass_kg == 0.0 {
                    0.0
                } else {
                    catalog.species_specific_enthalpy_change_j_kg(
                        &surface.binding.species,
                        context.temperature_k,
                        snapshot.temperature_k,
                    )?
                };
            integrated = Some(turnover);
            Ok(snapshot)
        },
        |surface, snapshot, reflected_heat, cancelled| {
            let (prepared, retained) = prepare(
                surface,
                snapshot,
                catalog,
                chamber_controls,
                reflected_heat,
                cancelled,
            )?;
            let final_phase = phase_at(
                surface.law,
                prepared.surface_temperature_k[1],
                snapshot.gas.molar_mass,
            )?;
            final_fraction = volume_fraction(
                prepared.solid.mass_kg,
                final_phase.solid_density_kg_m3,
                gas.volume_m3(),
                controls.maximum_deposit_volume_fraction,
            )?;
            storage = Some(retained);
            Ok(prepared)
        },
        |_, _| Ok(()),
        &mut cancelled,
    )?;
    let turnover = integrated.ok_or_else(|| invalid("missing integrated supported turnover"))?;
    let stored = storage.ok_or_else(|| invalid("missing original supported storage"))?;
    let GasChamberCandidate {
        gas,
        surfaces,
        step,
    } = candidate;
    let VaporReservoir::Finite { vapor_mass_kg, .. } = turnover.reservoir else {
        return Err(invalid(
            "integrated supported turnover lost its finite reservoir",
        ));
    };
    let actual_vapor_mass = gas.gas().species_mass_kg(&surface.binding.species)?;
    if (actual_vapor_mass == 0.0) != (vapor_mass_kg == 0.0)
        || (stored.state.deposit_mass_kg == 0.0) != (turnover.deposit_mass_kg == 0.0)
        || (actual_vapor_mass - vapor_mass_kg).abs() > controls.gas.mass_tolerance_kg
    {
        return Err(invalid(
            "original gas reconstruction changed the integrated inventory topology",
        ));
    }
    let mut exchanges = surfaces.into_iter();
    let exchange = exchanges
        .next()
        .ok_or_else(|| invalid("missing supported lifecycle endpoint"))?;
    if exchanges.next().is_some() {
        return Err(invalid(
            "supported lifecycle produced more than one surface",
        ));
    }
    let endpoint = GasSupportedChamberCandidate {
        gas,
        surfaces: vec![GasSupportedSurfaceStep {
            id: exchange.id,
            state: stored.state,
            exchange,
            substrate_energy_change_j: stored.substrate_energy_change_j,
            deposit_energy_change_j: stored.deposit_energy_change_j,
            reflected_heat_into_gas_j: 0.0,
        }],
        step,
    };
    poll(&mut cancelled)?;
    Ok(GasSupportedLifecycleCandidate {
        endpoint,
        turnover,
        deposit_volume_fraction: [initial_fraction, final_fraction],
    })
}

#[cfg(test)]
mod tests;
