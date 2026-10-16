//! Source-resolved humid air using the existing binary gas-mixture evaluator.
//!
//! Component query points describe the dry carrier and the dilute water-vapor
//! model at a common temperature and total mixture pressure. They do not assert
//! that pure water vapor is a stable phase at that total pressure. Relative
//! humidity is an explicit mixture input, separate from either component card.

use fs_matdb::{MaterialCard, QueryPoint};
use fs_qty::{Dims, QuantitySpec};

use crate::gas::{ConductivityModel, GasState, ResolvedGasState, resolve_sutherland_gas_state};
use crate::state_point::{MaterialPropertySelection, MaterialStatePointError};

/// A humid gas state together with both source-resolved component bundles.
///
/// The two parameter identities retain their own source receipts. Neither one
/// alone identifies the mixture: relative humidity and the explicit fixed
/// USSA/Eucken/Wilke/WMS model choices also determine the returned state.
#[derive(Clone, Debug, PartialEq)]
pub struct ResolvedMoistAirState {
    state: GasState,
    dry_air: ResolvedGasState,
    water_vapor: ResolvedGasState,
    relative_humidity: f64,
    water_mass_fraction: f64,
}

impl ResolvedMoistAirState {
    /// Derived medium for the existing acoustic and transport consumers.
    #[must_use]
    pub const fn state(&self) -> &GasState {
        &self.state
    }

    /// Dry-carrier parameters and their source usage receipts.
    #[must_use]
    pub const fn dry_air(&self) -> &ResolvedGasState {
        &self.dry_air
    }

    /// Dilute-vapor parameters and their source usage receipts.
    #[must_use]
    pub const fn water_vapor(&self) -> &ResolvedGasState {
        &self.water_vapor
    }

    /// Relative humidity as a fraction, using the Buck liquid-water basis.
    #[must_use]
    pub const fn relative_humidity(&self) -> f64 {
        self.relative_humidity
    }

    /// Water mass divided by total gas mass, distinct from RH and mole fraction.
    #[must_use]
    pub const fn water_mass_fraction(&self) -> f64 {
        self.water_mass_fraction
    }
}

/// Resolve the declared dry-air and water-vapor component models, then mix them.
///
/// Each component keeps its complete query and selection plan. The dry query
/// must explicitly select `source-composition-ussa1976 = 1`; the vapor query
/// must select `source-component-water-vapor = 1`, both dimensionless. Source
/// validity additionally enforces every card-specific condition. The dry-card
/// zero-RH condition applies to the component, not to the humid mixture.
/// Both queries must supply the same typed temperature and total pressure.
///
/// Conductivity uses USSA for the carrier and Eucken for the vapor; viscosity
/// uses Wilke mixing and conductivity uses Wassiljewa–Mason–Saxena mixing.
/// These are dilute-gas engineering approximations, with no condensation,
/// enhancement-factor, arbitrary-composition or measured-accuracy claim.
/// Evaluation is fixed-cost synchronous arithmetic after bounded resolution.
///
/// # Errors
/// Refuses missing/ambiguous/incompatible source claims, unsupported component
/// conditions, mismatched temperatures or pressures, invalid RH, and every
/// refusal of the existing moist gas evaluator. No partial result is returned.
pub fn resolve_moist_air_state(
    dry_card: &MaterialCard,
    dry_point: &QueryPoint,
    dry_selection: MaterialPropertySelection,
    vapor_card: &MaterialCard,
    vapor_point: &QueryPoint,
    vapor_selection: MaterialPropertySelection,
    relative_humidity: f64,
) -> Result<ResolvedMoistAirState, MaterialStatePointError> {
    for (point, flag) in [
        (dry_point, "source-composition-ussa1976"),
        (vapor_point, "source-component-water-vapor"),
    ] {
        if point.axes().get(flag) != Some(&1.0)
            || point.axis_quantities().get(flag) != Some(&QuantitySpec::dimensional(Dims::NONE))
        {
            return Err(MaterialStatePointError::InvalidRequirement {
                property: flag.into(),
                reason: "humid-air component identity needs its explicit dimensionless source flag",
            });
        }
    }
    let dry_air = resolve_sutherland_gas_state(
        dry_card,
        dry_point,
        ConductivityModel::Ussa1976AirFit,
        dry_selection,
    )?;
    let water_vapor = resolve_sutherland_gas_state(
        vapor_card,
        vapor_point,
        ConductivityModel::Eucken,
        vapor_selection,
    )?;
    let state =
        crate::gas::mix_moist_air(*dry_air.state(), *water_vapor.state(), relative_humidity)
            .map_err(|_| MaterialStatePointError::InvalidDerived {
                quantity: "humid-air component states, relative humidity or mixture model domain",
            })?;
    // w = x M_water / M_mix; M = R_universal / R_specific.
    let water_mass_fraction = state.water_mole_fraction * state.specific_gas_constant
        / water_vapor.state().specific_gas_constant;
    Ok(ResolvedMoistAirState {
        state,
        dry_air,
        water_vapor,
        relative_humidity,
        water_mass_fraction,
    })
}
