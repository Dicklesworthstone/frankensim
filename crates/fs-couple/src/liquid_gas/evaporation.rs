//! Original finite-rate phase transfer on native moving-liquid inventory.
//!
//! The existing Stefan-film and original gas caloric owners determine one
//! signed transfer. Native representability is selected before that solve;
//! latent and upstream-vapor sensible heat reach the actual liquid h once.
//! Native VOF transport then advances the same physical interval. The returned
//! candidate publishes the gas and native liquid together, without a second
//! film inventory or a reset of the atmosphere on subsequent calls.

use fs_exec::CancelGate;
use fs_lbm::material::free_surface::enthalpy::{
    FreeSurfaceEnthalpyError, FreeSurfaceEnthalpyReport, FreeSurfaceMassExchange,
    FreeSurfaceMassExchangeControls, FreeSurfaceMassExchangeReport, MaterialFreeSurfaceEnthalpy3,
};
use fs_material::evaporation::gas::{
    GasEvaporationError, GasEvaporationSolve, GasEvaporationStep, GasLiquidBinding,
    advance_surface_with_gas_mass_selection,
};
use fs_material::thermochemical_gas::{ThermochemicalGasModel, energy::parcel::ClosedGasParcel};
use fs_thermal::evaporation::{EvaporationError, EvaporationLaw, LiquidSurfaceState};
use std::cell::RefCell;

/// One measured Cartesian liquid/gas patch at the beginning of the interval.
#[derive(Clone, Copy)]
pub struct LiquidGasEvaporationSurface<'a> {
    pub liquid: [usize; 3],
    /// Face-adjacent in-domain gas; walls and occupied neighbors refuse.
    pub gas: [usize; 3],
    pub area_m2: f64,
    pub mass_transfer_coefficient_m_s: f64,
    pub law: &'a dyn EvaporationLaw,
    /// Explicit volatile identity and constant native-liquid/gas datum offset.
    pub binding: &'a GasLiquidBinding,
}

#[derive(Clone, Copy, Debug)]
pub struct EvaporatingLiquidGasSolve {
    pub phase: GasEvaporationSolve,
    pub native: FreeSurfaceMassExchangeControls,
    /// Whole physical interval, including actual native transport drift [kg].
    pub mass_tolerance_kg: f64,
    /// Whole source-relative thermal defect after native transport [J].
    pub energy_tolerance_j: f64,
}

#[derive(Debug)]
pub enum EvaporatingLiquidGasError {
    Phase(GasEvaporationError),
    Liquid(FreeSurfaceEnthalpyError),
    Invalid(&'static str),
    MassBalance { residual_kg: f64, tolerance_kg: f64 },
    EnergyBalance { residual_j: f64, tolerance_j: f64 },
    Cancelled,
}
impl core::fmt::Display for EvaporatingLiquidGasError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "finite-gas/native-liquid phase transfer: {self:?}")
    }
}
impl std::error::Error for EvaporatingLiquidGasError {}
impl From<FreeSurfaceEnthalpyError> for EvaporatingLiquidGasError {
    fn from(value: FreeSurfaceEnthalpyError) -> Self {
        Self::Liquid(value)
    }
}
impl From<GasEvaporationError> for EvaporatingLiquidGasError {
    fn from(value: GasEvaporationError) -> Self {
        if value == GasEvaporationError::Cancelled {
            Self::Cancelled
        } else {
            Self::Phase(value)
        }
    }
}

/// Original owner reports and independent whole-interval conservation checks.
#[derive(Clone, Debug)]
pub struct EvaporatingLiquidGasReport {
    pub duration_s: f64,
    pub phase: GasEvaporationStep,
    pub mass_exchange: FreeSurfaceMassExchangeReport,
    pub flow: FreeSurfaceEnthalpyReport,
    /// Actual native liquid change in its retained caloric datum [J].
    pub native_liquid_enthalpy_change_j: f64,
    /// Native phase endpoint minus the original phase candidate's sensible
    /// storage, evaluated before adding any formation reference [J].
    pub native_phase_caloric_defect_j: f64,
    /// Actual gas + native phase + subsequent native-flow mass change [kg].
    pub mass_residual_kg: f64,
    /// Same account using the actual volatile species instead of total gas [kg].
    pub volatile_mass_residual_kg: f64,
    /// Original phase defect + native caloric adjustment + native flow defect.
    /// Gas-species representation rounding remains separately mass-budgeted,
    /// as in the original source-relative boundary owner.
    pub energy_residual_j: f64,
}

pub struct EvaporatingLiquidGasCandidate<'s, 'm, 'catalog> {
    liquid: MaterialFreeSurfaceEnthalpy3<'s, 'm>,
    gas: ClosedGasParcel,
    catalog: &'catalog ThermochemicalGasModel,
    report: EvaporatingLiquidGasReport,
}
impl<'s, 'm, 'catalog> EvaporatingLiquidGasCandidate<'s, 'm, 'catalog> {
    pub const fn liquid(&self) -> &MaterialFreeSurfaceEnthalpy3<'s, 'm> {
        &self.liquid
    }
    pub const fn gas(&self) -> &ClosedGasParcel {
        &self.gas
    }
    pub const fn catalog(&self) -> &'catalog ThermochemicalGasModel {
        self.catalog
    }
    pub const fn report(&self) -> &EvaporatingLiquidGasReport {
        &self.report
    }
    pub fn into_parts(
        self,
    ) -> (
        MaterialFreeSurfaceEnthalpy3<'s, 'm>,
        ClosedGasParcel,
        &'catalog ThermochemicalGasModel,
        EvaporatingLiquidGasReport,
    ) {
        (self.liquid, self.gas, self.catalog, self.report)
    }
}

fn poll(gate: &CancelGate) -> Result<(), EvaporatingLiquidGasError> {
    if gate.is_requested() {
        Err(EvaporatingLiquidGasError::Cancelled)
    } else {
        Ok(())
    }
}

/// Advance phase transfer and original native flow over ONE retained timestep.
///
/// The phase owner supplies all Stefan-film, latent, vapor-sensible and gas
/// caloric physics. Its h-to-T callback is the actual native liquid chart, and
/// its mass selector is the actual native VOF representation. The selected dm
/// debits the finite gas and credits native liquid once; no intermediate film
/// survives. The native source owner refreshes viscosity, buoyancy and supported
/// capillary bindings before its normal momentum/enthalpy transport.
/// Every refusal or cancellation discards the complete private candidate.
///
/// This is a first-order frozen-film/flow split on one caller-measured patch;
/// refine the original physical timestep for accuracy. Gas volume is a fixed
/// connected control volume, not the complement of the changing liquid fill.
/// The native constant-density volume law remains authoritative; phase-law
/// recession is a diagnostic and does not move a second interface.
/// No sealed-vessel displacement, gas traction/recoil, added sensible film
/// conduction, nucleation, or interface overflow/exhaustion remap is inferred.
/// Material crosses at the original liquid's local velocity. Existing solved
/// surface-stress/contact histories remain subject to native exchange refusal.
/// Attached electrical cores require their original conduction timeline and
/// therefore refuse this direct native-flow entry point.
#[allow(clippy::too_many_arguments)]
pub fn advance_liquid_with_gas<'s, 'm, 'catalog>(
    liquid: &MaterialFreeSurfaceEnthalpy3<'s, 'm>,
    gas: &ClosedGasParcel,
    catalog: &'catalog ThermochemicalGasModel,
    surface: LiquidGasEvaporationSurface<'_>,
    controls: EvaporatingLiquidGasSolve,
    gate: &CancelGate,
) -> Result<EvaporatingLiquidGasCandidate<'s, 'm, 'catalog>, EvaporatingLiquidGasError> {
    type Error = EvaporatingLiquidGasError;
    poll(gate)?;
    if [controls.mass_tolerance_kg, controls.energy_tolerance_j]
        .iter()
        .any(|x| !x.is_finite() || *x < 0.0)
    {
        return Err(Error::Invalid(
            "finite nonnegative whole-interval mass and energy budgets required",
        ));
    }
    liquid.admit_interface_gas_patch(surface.liquid, surface.gas, surface.area_m2)?;
    let initial = liquid.project_interface_mass_exchange(surface.liquid, 0.0)?;
    let h0 = liquid
        .specific_enthalpy_at(surface.liquid)
        .ok_or(Error::Invalid(
            "native phase interface has no retained enthalpy",
        ))?;
    let duration_s = liquid.flow().scaling().timestep_s();
    let chart_error = RefCell::new(None);
    let mut mass_error = None;
    let mut selected = None;
    let proposal = advance_surface_with_gas_mass_selection(
        surface.law,
        surface.binding,
        LiquidSurfaceState {
            mass_kg: initial.initial_mass_kg,
            specific_enthalpy_j_kg: h0,
            area_m2: surface.area_m2,
        },
        gas,
        catalog,
        |h| {
            liquid
                .liquid_temperature_at_specific_enthalpy_k(h)
                .map_err(|error| {
                    *chart_error.borrow_mut() = Some(error);
                    EvaporationError::Invalid("original native liquid caloric branch refused")
                })
        },
        surface.mass_transfer_coefficient_m_s,
        duration_s,
        controls.phase,
        |dm| match liquid.project_interface_mass_exchange(surface.liquid, -dm) {
            Ok(projection) => {
                selected = Some(projection);
                Ok(-projection.mass_change_kg)
            }
            Err(error) => {
                mass_error = Some(error);
                Err(GasEvaporationError::Invalid(
                    "original native interface mass selection refused",
                ))
            }
        },
        || gate.is_requested(),
    );
    if let Some(error) = chart_error.into_inner().or(mass_error) {
        return Err(error.into());
    }
    let phase = proposal?;
    let projection = selected.ok_or(Error::Invalid("phase owner did not select native mass"))?;
    let (phase_liquid, next_gas, catalog, phase) = phase.into_parts();
    // The native post-transfer mass owns its actual caloric update. Any last-bit
    // difference from SI subtraction in the phase owner is audited below.
    let target_delta_h = phase_liquid.specific_enthalpy_j_kg - h0;
    let heat = projection.final_mass_kg * target_delta_h;
    if !heat.is_finite() || (target_delta_h != 0.0 && heat == 0.0) {
        return Err(Error::Invalid(
            "native phase caloric increment is unrepresentable",
        ));
    }
    let mut next_liquid = liquid.clone();
    let mass_exchange = next_liquid.exchange_interface_mass(
        &[FreeSurfaceMassExchange {
            coordinates: surface.liquid,
            mass_kg: projection.mass_change_kg,
            specific_enthalpy_j_kg: h0,
            heat_j: heat,
        }],
        controls.native,
        gate,
    )?;
    if mass_exchange.mass_change_kg != -phase.transferred_mass_kg {
        return Err(Error::Invalid(
            "native phase transfer differs from its selected gas debit",
        ));
    }
    let stored_delta_h = mass_exchange.samples[0].specific_enthalpy_j_kg[1] - h0;
    let native_phase_caloric_defect_j =
        projection.final_mass_kg * stored_delta_h - phase_liquid.mass_kg * target_delta_h;
    poll(gate)?;
    let flow = next_liquid.step(gate)?;
    let native_transport_mass = flow.mass_after_kg - flow.mass_before_kg;
    let native_mass = mass_exchange.mass_change_kg + native_transport_mass;
    let mass_residual_kg = native_mass + (next_gas.mass_kg() - gas.mass_kg());
    let volatile_mass_residual_kg = native_mass
        + (next_gas
            .species_mass_kg(&surface.binding.species)
            .map_err(GasEvaporationError::from)?
            - gas
                .species_mass_kg(&surface.binding.species)
                .map_err(GasEvaporationError::from)?);
    let energy_residual_j =
        phase.energy_residual_j + native_phase_caloric_defect_j + flow.balance_defect_j;
    let native_liquid_enthalpy_change_j =
        mass_exchange.enthalpy_change_j + (flow.enthalpy_after_j - flow.enthalpy_before_j);
    if ![
        native_phase_caloric_defect_j,
        mass_residual_kg,
        volatile_mass_residual_kg,
        energy_residual_j,
        native_liquid_enthalpy_change_j,
    ]
    .iter()
    .all(|x| x.is_finite())
    {
        return Err(Error::Invalid(
            "combined native-liquid/gas account is unrepresentable",
        ));
    }
    let worst_mass = if mass_residual_kg.abs() >= volatile_mass_residual_kg.abs() {
        mass_residual_kg
    } else {
        volatile_mass_residual_kg
    };
    if worst_mass.abs() > controls.mass_tolerance_kg {
        return Err(Error::MassBalance {
            residual_kg: worst_mass,
            tolerance_kg: controls.mass_tolerance_kg,
        });
    }
    if energy_residual_j.abs() > controls.energy_tolerance_j {
        return Err(Error::EnergyBalance {
            residual_j: energy_residual_j,
            tolerance_j: controls.energy_tolerance_j,
        });
    }
    poll(gate)?;
    Ok(EvaporatingLiquidGasCandidate {
        liquid: next_liquid,
        gas: next_gas,
        catalog,
        report: EvaporatingLiquidGasReport {
            duration_s,
            phase,
            mass_exchange,
            flow,
            native_liquid_enthalpy_change_j,
            native_phase_caloric_defect_j,
            mass_residual_kg,
            volatile_mass_residual_kg,
            energy_residual_j,
        },
    })
}
