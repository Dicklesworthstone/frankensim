//! Finite Henry dissolution/exsolution with the original thermochemical gas.
//!
//! C* = H_s^cp p_i uses the SELECTED species partial pressure and the IUPAC
//! concentration/pressure solubility convention, mol/(m³ Pa). A declared
//! quasi-steady liquid film gives dn_liquid/dt = A K_l(C* - C). Without gas-side
//! resistance K_l = D/delta. An explicit pressure-based gas coefficient k_g
//! instead gives 1/K_l = delta/D + H/k_g: the two physical film resistances
//! are in series, with Henry equilibrium at their common interface.
//! With fixed liquid/gas volumes, T, H, D and k_g, the native ideal-gas relation
//! closes an affine finite-inventory ODE. Its exponential update includes the
//! changing gas partial pressure, conserves available solute, and approaches
//! equilibrium without an infinite gas reservoir or an explicit-Euler clamp.
//!
//! Other gas species retain their original inventories. The existing species
//! boundary owner reconstructs composition, pressure, transport and acoustics;
//! its source-relative energy input cancels flow work to hold T exactly.
//! Original NASA-9 gas calorics and an explicitly aligned dissolved PARTIAL
//! enthalpy chart determine gas and liquid thermostat heat separately. The
//! solvent stays at fixed T and volume. Negligible solute partial molar volume
//! identifies dissolved h with u; this assumption is not a guessed liquid EOS.
//! Streams use gas-side partial h in both directions; dissolution heat is
//! assigned to the solution thermostat. The combined thermostat heat is invariant.
//!
//! An optional source-qualified gas film supplies k_g = D_G/(R T delta_G):
//! K_L = 1/(delta_L/D_L + H_s^cp R T delta_G/D_G). Both film traces use
//! one common molar flux; the bulk gas pressure is not the interface pressure.
//! Supply either that original diffusivity source or the explicit coefficient,
//! never both for the same film.
//!
//! Applicability: dilute ideal solution, one exchanging noncondensable species,
//! existing gas volume and already quasi-steady films. Gas-side resistance is
//! negligible by default; an explicit gas film additionally assumes dilute
//! Fick diffusion in its declared carrier, without Stefan-flow correction.
//! A supplied coefficient must apply to the entire interval; no pressure- or
//! composition-dependent coefficient is silently frozen. A bare coefficient
//! has no gas diffusion time; only a gas diffusivity and thickness admit that
//! additional timescale check. Neither check certifies a gas-film transient.
//! The authored film-relaxation/concentration bounds
//! do not validate those assumptions. This model supplies no nucleation,
//! moving bubble radius, transient diffusion field, evaporation, gas reaction,
//! or adiabatic thermal evolution. No entropy-production theorem is claimed.
//!
//! Sources: IUPAC Gold Book, H_s^cp, https://goldbook.iupac.org/terms/view/09024 ;
//! MIT 1.85, gas transfer through a water-side film,
//! https://ocw.mit.edu/courses/1-85-water-and-wastewater-treatment-engineering-spring-2006/ef48e6d3cd622b8fc2c0adcc652a0817_midterm_review.pdf .
//! Two-film interface balances and series resistance: ETH Zurich, film model,
//! https://www.hyper-tvt.ethz.ch/fundamentals-masstransfer-filmodel.html .

/// Simultaneous absorption/exsolution of independent dilute species into one liquid.
pub mod mixture;

mod gas_film;
use gas_film::{FilmResponse, ResolvedGasFilm};
pub use gas_film::{GAS_FILM_DIFFUSIVITY_PROPERTY, GAS_FILM_MOLE_FRACTION_AXIS, GasFilmSource};

use crate::phase::{PhaseStateError, ResolvedHeatCapacityCurve};
use crate::state_point::{
    MaterialPropertySelection, MaterialStatePointError, ResolvedMaterialStatePoint,
    ScalarAdmissibility, ScalarPropertyRequirement, resolve_material_state_point,
};
use crate::thermochemical_gas::ThermochemicalGasModel;
use crate::thermochemical_gas::energy::{
    GasEnergyError,
    parcel::{
        ClosedGasParcel,
        transfer::species::boundary::{
            GasSpeciesBoundarySolve, GasSpeciesBoundaryStep, SpeciesBoundaryRelativeTransfer,
        },
    },
};
use fs_matdb::{MaterialCard, QueryPoint};
use fs_math::det;
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_qty::{Dims, Pressure, QuantitySpec, Temperature};
use fs_thermochem::{SpeciesId, UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R};

// Reuse the existing pressure-normalized molar-flux quantity. This coefficient
// carries no membrane model, source authority or implicit STP conversion.
pub use crate::thermochemical_gas::energy::parcel::transfer::species::network::transport::membrane::MolarPermeance;

pub const HENRY_SOLUBILITY_PROPERTY: &str = "henry_molar_solubility";
pub const HENRY_SOLUBILITY_DIMS: Dims = Dims([-2, -1, 2, 0, 0, 1]);
pub const DISSOLVED_DIFFUSIVITY_PROPERTY: &str = "dissolved_gas_diffusivity";
pub const DIFFUSIVITY_DIMS: Dims = Dims([2, 0, -1, 0, 0, 0]);
pub const CONCENTRATION_DIMS: Dims = Dims([-3, 0, 0, 0, 0, 1]);

/// The card must describe this solute, solvent and actual mixture conditions.
/// H and D may depend on canonical temperature; they must remain constant as
/// pressure and concentration change at the selected isothermal state.
#[derive(Clone, Debug, PartialEq)]
pub struct HenrySource {
    pub material: MaterialCard,
    pub point: QueryPoint,
    pub selection: MaterialPropertySelection,
}

/// Original authored partial specific enthalpy of the DISSOLVED species in
/// this solution, not the neat solute's heat capacity or a guessed latent heat.
/// The explicit offset aligns its datum with the retained NASA-9 gas species.
/// Density carried by the generic chart is not used to infer solution volume.
#[derive(Clone, Debug, PartialEq)]
pub struct DissolvedGasCalorics {
    pub partial_enthalpy: ResolvedHeatCapacityCurve,
    pub enthalpy_offset_j_kg: f64,
}

#[derive(Clone, Debug, PartialEq)]
pub enum DissolvedGasError {
    Invalid(&'static str),
    Material(MaterialStatePointError),
    Caloric(PhaseStateError),
    Gas(GasEnergyError),
    Concentration {
        actual_mol_m3: f64,
        maximum_mol_m3: f64,
    },
    FilmRelaxation {
        ratio: f64,
        maximum: f64,
    },
    /// Gas diffusion cannot be treated as quasi-steady at the authored limit.
    GasFilmRelaxation {
        /// Gas-film diffusion time divided by the common bulk relaxation time.
        ratio: f64,
        /// Authored maximum of that dimensionless ratio.
        maximum: f64,
    },
    /// A gas bulk/interface or simultaneous-mixture dilution bound is exceeded.
    GasFilmConcentration {
        /// Selected or conservatively summed exchanging-species mole fraction.
        actual_mole_fraction: f64,
        /// Authored dilute-Fick upper limit, strictly below one.
        maximum_mole_fraction: f64,
    },
    Balance {
        quantity: &'static str,
        residual: f64,
        tolerance: f64,
    },
    Cancelled,
}
impl From<MaterialStatePointError> for DissolvedGasError {
    fn from(e: MaterialStatePointError) -> Self {
        Self::Material(e)
    }
}
impl From<PhaseStateError> for DissolvedGasError {
    fn from(e: PhaseStateError) -> Self {
        Self::Caloric(e)
    }
}
impl From<GasEnergyError> for DissolvedGasError {
    fn from(e: GasEnergyError) -> Self {
        if e == GasEnergyError::Cancelled {
            Self::Cancelled
        } else {
            Self::Gas(e)
        }
    }
}
impl core::fmt::Display for DissolvedGasError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "finite Henry gas exchange: {self:?}")
    }
}
impl std::error::Error for DissolvedGasError {}
type Result<T> = core::result::Result<T, DissolvedGasError>;
fn positive(x: f64) -> bool {
    x.is_finite() && x > 0.0
}
fn nonnegative(x: f64) -> bool {
    x.is_finite() && x >= 0.0
}
fn invalid(s: &'static str) -> DissolvedGasError {
    DissolvedGasError::Invalid(s)
}
fn finite(x: f64) -> Result<f64> {
    if x.is_finite() {
        Ok(x)
    } else {
        Err(invalid("nonfinite dissolved-gas arithmetic"))
    }
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<()> {
    if cancelled() {
        Err(DissolvedGasError::Cancelled)
    } else {
        Ok(())
    }
}
fn budget(quantity: &'static str, residual: f64, tolerance: f64) -> Result<()> {
    if residual.is_finite() && residual.abs() <= tolerance {
        Ok(())
    } else {
        Err(DissolvedGasError::Balance {
            quantity,
            residual,
            tolerance,
        })
    }
}
fn typed(kind: QuantityKind) -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(kind, ValueForm::Static))
}

fn admit_frozen_sources(resolved: &ResolvedMaterialStatePoint) -> Result<()> {
    for property in resolved.properties() {
        if !property.is_constant_at_fixed_temperature()
            || !property.is_static_at_fixed_external_conditions()
            || property
                .source_axes()
                .iter()
                .any(|(axis, _)| axis != "temperature")
            || property
                .answer()
                .evidence
                .model
                .validity
                .bounds()
                .keys()
                .any(|axis| {
                    let q = property.answer().receipt.axis_quantities.get(axis);
                    matches!(axis.as_str(), "T" | "p" | "P")
                        || (axis != "temperature"
                            && q.is_some_and(|q| q.dims() == Temperature::DIMS))
                        || (axis != "pressure" && q.is_some_and(|q| q.dims() == Pressure::DIMS))
                })
        {
            return Err(invalid(
                "isothermal Henry exchange needs constant-at-fixed-T sources and canonical axes",
            ));
        }
    }
    Ok(())
}

/// Original constitutive sources frozen at one explicitly authored temperature.
#[derive(Clone, Debug, PartialEq)]
pub struct HenryDissolution {
    species: SpeciesId,
    catalog: ThermochemicalGasModel,
    source: HenrySource,
    resolved: ResolvedMaterialStatePoint,
    calorics: DissolvedGasCalorics,
    temperature_k: f64,
    molar_mass_kg_mol: f64,
    henry_mol_m3_pa: f64,
    diffusivity_m2_s: f64,
    dissolved_enthalpy_j_kg: f64,
    gas_film: Option<ResolvedGasFilm>,
}

impl HenryDissolution {
    pub fn new(
        catalog: &ThermochemicalGasModel,
        species: SpeciesId,
        source: HenrySource,
        calorics: DissolvedGasCalorics,
    ) -> Result<Self> {
        let gas = catalog
            .thermodynamics()
            .components()
            .iter()
            .find(|s| s.species() == &species)
            .ok_or(invalid(
                "selected solute requires its original gas species catalog",
            ))?;
        if catalog.water_species() == Some(&species) {
            return Err(invalid(
                "water vapor phase change is not noncondensable Henry exchange",
            ));
        }
        let temperature_k = *source
            .point
            .axes()
            .get("temperature")
            .ok_or(invalid("Henry source requires canonical temperature"))?;
        if !positive(temperature_k)
            || !calorics.enthalpy_offset_j_kg.is_finite()
            || source.point.axis_quantities().get("temperature")
                != Some(&typed(QuantityKind::AbsoluteTemperature))
            || source.point.axis_quantities().get("pressure")
                != Some(&typed(QuantityKind::Pressure))
        {
            return Err(invalid(
                "typed canonical T/p and an explicit caloric datum required",
            ));
        }
        let requirements = [
            ScalarPropertyRequirement::try_new(
                HENRY_SOLUBILITY_PROPERTY,
                HENRY_SOLUBILITY_DIMS,
                ScalarAdmissibility::StrictlyPositive,
            )?,
            ScalarPropertyRequirement::try_new(
                DISSOLVED_DIFFUSIVITY_PROPERTY,
                DIFFUSIVITY_DIMS,
                ScalarAdmissibility::NonNegative,
            )?,
        ];
        let resolved = resolve_material_state_point(
            &source.material,
            &source.point,
            &requirements,
            source.selection.clone(),
        )?;
        admit_frozen_sources(&resolved)?;
        let henry = resolved
            .property(HENRY_SOLUBILITY_PROPERTY)
            .unwrap()
            .value_si();
        let diffusivity = resolved
            .property(DISSOLVED_DIFFUSIVITY_PROPERTY)
            .unwrap()
            .value_si();
        let dissolved_h = finite(
            calorics
                .partial_enthalpy
                .state_at_temperature(temperature_k)?
                .specific_enthalpy_j_kg()
                + calorics.enthalpy_offset_j_kg,
        )?;
        catalog.species_specific_enthalpy_j_kg(&species, temperature_k)?;
        Ok(Self {
            species,
            catalog: catalog.clone(),
            source,
            resolved,
            calorics,
            temperature_k,
            molar_mass_kg_mol: gas.molar_mass().value(),
            henry_mol_m3_pa: henry,
            diffusivity_m2_s: diffusivity,
            dissolved_enthalpy_j_kg: dissolved_h,
            gas_film: None,
        })
    }
    pub fn species(&self) -> &SpeciesId {
        &self.species
    }
    pub const fn source(&self) -> &HenrySource {
        &self.source
    }
    pub const fn calorics(&self) -> &DissolvedGasCalorics {
        &self.calorics
    }
    pub const fn catalog(&self) -> &ThermochemicalGasModel {
        &self.catalog
    }
    pub const fn temperature_k(&self) -> f64 {
        self.temperature_k
    }
    pub const fn henry_mol_m3_pa(&self) -> f64 {
        self.henry_mol_m3_pa
    }
    pub const fn diffusivity_m2_s(&self) -> f64 {
        self.diffusivity_m2_s
    }

    fn sample(
        &self,
        gas: &ClosedGasParcel,
        liquid: DissolvedGasInventory,
        controls: DissolvedGasControls,
        film: FilmResponse,
    ) -> Result<DissolutionSample> {
        if gas.state().state().temperature != self.temperature_k {
            return Err(invalid(
                "both phases must retain the exact isothermal source temperature",
            ));
        }
        let pressure = gas.species_partial_pressure_pa(&self.species)?;
        let concentration = finite(liquid.dissolved_moles / liquid.liquid_volume_m3)?;
        let equilibrium = finite(self.henry_mol_m3_pa * pressure)?;
        let interface = finite(film.liquid_weight * equilibrium + film.gas_weight * concentration)?;
        let interface_pressure = if film.gas_weight == 0.0 {
            pressure
        } else {
            finite(interface / self.henry_mol_m3_pa)?
        };
        if interface > 0.0 && interface_pressure == 0.0 {
            return Err(invalid("positive Henry interface pressure underflows"));
        }
        // A conducting constant-total-pressure gas film cannot have a species
        // partial pressure larger than its total pressure. A sealed path has
        // only a Henry-equivalent pressure on its disconnected liquid side.
        if film.conductance_m3_s > 0.0
            && interface_pressure > gas.state().state().pressure * (1.0 + 32.0 * f64::EPSILON)
        {
            return Err(invalid(
                "gas-film interfacial partial pressure exceeds total pressure",
            ));
        }
        for c in [concentration, interface] {
            if !nonnegative(c) || c > controls.maximum_concentration_mol_m3 {
                return Err(DissolvedGasError::Concentration {
                    actual_mol_m3: c,
                    maximum_mol_m3: controls.maximum_concentration_mol_m3,
                });
            }
        }
        let point = |c| {
            self.source
                .point
                .clone()
                .with_quantity(
                    "pressure",
                    typed(QuantityKind::Pressure),
                    gas.state().state().pressure,
                )
                .and_then(|p| {
                    p.with_quantity(
                        "dissolved_concentration",
                        typed(QuantityKind::AmountConcentration),
                        c,
                    )
                })
                .map_err(|_| invalid("unrepresentable actual Henry source conditions"))
        };
        // Both liquid bulk and interfacial equilibrium must lie in the SAME
        // original claims' support; a supported dilute bulk cannot hide a
        // nondilute interfacial concentration.
        let bulk = self
            .resolved
            .at_query_point(&self.source.material, &point(concentration)?)?;
        let interface_material = self
            .resolved
            .at_query_point(&self.source.material, &point(interface)?)?;
        for at in [&bulk, &interface_material] {
            if at.property(HENRY_SOLUBILITY_PROPERTY).unwrap().value_si() != self.henry_mol_m3_pa
                || at
                    .property(DISSOLVED_DIFFUSIVITY_PROPERTY)
                    .unwrap()
                    .value_si()
                    != self.diffusivity_m2_s
            {
                return Err(invalid(
                    "frozen Henry/diffusivity source changed during an isothermal step",
                ));
            }
        }
        let gas_film_material = self
            .gas_film
            .as_ref()
            .map(|source| -> Result<_> {
                Ok([
                    source.sample(gas.state().state().pressure, pressure)?,
                    source.sample(gas.state().state().pressure, interface_pressure)?,
                ])
            })
            .transpose()?;
        let flux = finite(film.coefficient_m_s * (equilibrium - concentration))?;
        if film.coefficient_m_s > 0.0 && equilibrium != concentration && flux == 0.0 {
            return Err(invalid("nonzero film flux underflows"));
        }
        Ok(DissolutionSample {
            partial_pressure_pa: pressure,
            concentration_mol_m3: concentration,
            equilibrium_concentration_mol_m3: equilibrium,
            interface_partial_pressure_pa: interface_pressure,
            interface_concentration_mol_m3: interface,
            molar_flux_mol_m2_s: flux,
            material: [bulk, interface_material],
            gas_film_material,
        })
    }

    /// Exact fixed-coefficient finite-pair relaxation, followed by ONE original
    /// gas reconstruction. Both histories publish together after all source,
    /// inventory, thermal and final cancellation checks. The gas has no separate
    /// clock; the supplied finite-liquid clock owns this exchange interval.
    ///
    /// A and delta are physical authored film geometry, not a fitted rate.
    /// delta²/D divided by the bulk relaxation time must meet the declared
    /// quasi-steady limit. A new initially uniform diffusion layer is outside
    /// this model even if the timescale ratio is small.
    pub fn advance(
        &self,
        gas: &mut ClosedGasParcel,
        liquid: &mut DissolvedGasInventory,
        film: DissolutionFilm,
        end_time_s: f64,
        controls: DissolvedGasControls,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<DissolutionStep> {
        let prepared = self.prepare(gas, *liquid, film, end_time_s, controls, &mut cancelled)?;
        let candidate = gas.exchange_species_boundary_relative(
            &self.catalog,
            &[prepared.transfer(self)],
            controls.gas,
            &mut cancelled,
        )?;
        let exchange = self.finish(prepared, candidate.gas(), controls)?;
        let energy_residual = finite(
            (candidate.step().stored_energy_change_j + exchange.dissolved_internal_energy_change_j)
                - exchange.thermostat_heat_into_system_j,
        )?;
        budget(
            "combined internal energy",
            energy_residual,
            controls.gas.energy_tolerance_j,
        )?;
        let (next_gas, _, gas_step) = candidate.into_parts();
        let step = DissolutionStep {
            initial: exchange.initial,
            state: exchange.state,
            gas: gas_step,
            sources: exchange.sources,
            dissolved_moles_change: exchange.dissolved_moles_change,
            relaxation_rate_per_s: exchange.relaxation_rate_per_s,
            film_relaxation_ratio: exchange.film_relaxation_ratio,
            gas_film_relaxation_ratio: exchange.gas_film_relaxation_ratio,
            partial_enthalpy_into_liquid_j: exchange.partial_enthalpy_into_liquid_j,
            dissolved_internal_energy_change_j: exchange.dissolved_internal_energy_change_j,
            thermostat_heat_into_gas_j: exchange.thermostat_heat_into_gas_j,
            thermostat_heat_into_liquid_j: exchange.thermostat_heat_into_liquid_j,
            thermostat_heat_into_system_j: exchange.thermostat_heat_into_system_j,
            mass_residual_kg: exchange.mass_residual_kg,
            energy_residual_j: energy_residual,
        };
        poll(&mut cancelled)?;
        *gas = next_gas;
        *liquid = step.state;
        Ok(step)
    }

    // Both the scalar and simultaneous owners prepare from the unchanged gas.
    // No intermediate composition, clock, source resolve or energy can publish.
    fn prepare(
        &self,
        gas: &ClosedGasParcel,
        initial: DissolvedGasInventory,
        film: DissolutionFilm,
        end_time_s: f64,
        controls: DissolvedGasControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<PreparedDissolution> {
        poll(cancelled)?;
        let dt = end_time_s - initial.time_s;
        controls.gas.gas.validate()?;
        if !positive(dt)
            || !end_time_s.is_finite()
            || !positive(controls.maximum_step_s)
            || dt > controls.maximum_step_s
            || !positive(controls.maximum_concentration_mol_m3)
            || !positive(controls.maximum_film_relaxation_ratio)
            || controls.maximum_film_relaxation_ratio >= 1.0
            || !positive(controls.gas.mass_tolerance_kg)
            || !positive(controls.gas.energy_tolerance_j)
            || !nonnegative(film.area_m2)
            || !positive(film.thickness_m)
        {
            return Err(invalid(
                "positive bounded time, dilute/film limits, budgets and physical film geometry required",
            ));
        }
        let initial_mass = gas.species_mass_kg(&self.species)?;
        let capacity = finite(
            1.0 / initial.liquid_volume_m3
                + self.henry_mol_m3_pa * (R * self.temperature_k / gas.volume_m3()),
        )?;
        if !positive(capacity) {
            return Err(invalid("positive finite liquid/gas capacity required"));
        }
        let response = self.film_response(film, capacity)?;
        let before = self.sample(gas, initial, controls, response)?;
        let rate = finite(response.conductance_m3_s * capacity)?;
        let ratio = response.liquid_relaxation_ratio;
        if rate > 0.0 && ratio > controls.maximum_film_relaxation_ratio {
            return Err(DissolvedGasError::FilmRelaxation {
                ratio,
                maximum: controls.maximum_film_relaxation_ratio,
            });
        }
        if let Some(gas_ratio) = response.gas_relaxation_ratio
            && rate > 0.0
            && gas_ratio > controls.maximum_film_relaxation_ratio
        {
            return Err(DissolvedGasError::GasFilmRelaxation {
                ratio: gas_ratio,
                maximum: controls.maximum_film_relaxation_ratio,
            });
        }
        if response.conductance_m3_s > 0.0 && !positive(rate) {
            return Err(invalid("positive film relaxation underflows or overflows"));
        }
        let drive = before.equilibrium_concentration_mol_m3 - before.concentration_mol_m3;
        let exponent = finite(rate * dt)?;
        let fraction = -det::expm1(-exponent);
        let requested = finite((drive / capacity) * fraction)?;
        if (rate > 0.0 && (exponent == 0.0 || fraction == 0.0))
            || (drive != 0.0 && fraction != 0.0 && requested == 0.0)
        {
            return Err(invalid("positive finite diffusion relaxation underflows"));
        }
        let next_moles = finite(initial.dissolved_moles + requested)?;
        let dn = next_moles - initial.dissolved_moles;
        let dm = finite(dn * self.molar_mass_kg_mol)?;
        if !nonnegative(next_moles)
            || (initial.dissolved_moles > 0.0 && next_moles == 0.0)
            || (initial_mass > 0.0 && initial_mass - dm == 0.0)
            || (requested != 0.0 && (dn == 0.0 || dm == 0.0))
            || initial_mass - dm < 0.0
        {
            return Err(invalid(
                "finite solute transfer depletes or loses a represented inventory",
            ));
        }
        budget(
            "represented kinetic transfer",
            (dn - requested) * self.molar_mass_kg_mol,
            controls.gas.mass_tolerance_kg,
        )?;
        let next_liquid = DissolvedGasInventory {
            dissolved_moles: next_moles,
            time_s: end_time_s,
            ..initial
        };
        // Match the original boundary owner's exact operation order for
        // h_i-u_i=R*T/M_i. Their sum is zero sensible heat, so T is unchanged.
        let gas_bath_heat = finite(dm * (R * self.temperature_k / self.molar_mass_kg_mol))?;
        let hg = self
            .catalog
            .species_specific_enthalpy_j_kg(&self.species, self.temperature_k)?;
        let stream_h = finite(dm * hg)?;
        let dissolved_change = finite(dm * self.dissolved_enthalpy_j_kg)?;
        let liquid_bath_heat = finite(dm * (self.dissolved_enthalpy_j_kg - hg))?;
        let bath_heat = finite(gas_bath_heat + liquid_bath_heat)?;
        Ok(PreparedDissolution {
            initial_gas_mass_kg: initial_mass,
            mass_into_liquid_kg: dm,
            capacity,
            drive,
            film: response,
            exchange: DissolutionExchange {
                species: self.species.clone(),
                initial,
                state: next_liquid,
                sources: [before.clone(), before],
                dissolved_moles_change: dn,
                relaxation_rate_per_s: rate,
                film_relaxation_ratio: ratio,
                gas_film_relaxation_ratio: response.gas_relaxation_ratio,
                partial_enthalpy_into_liquid_j: stream_h,
                dissolved_internal_energy_change_j: dissolved_change,
                thermostat_heat_into_gas_j: gas_bath_heat,
                thermostat_heat_into_liquid_j: liquid_bath_heat,
                thermostat_heat_into_system_j: bath_heat,
                mass_residual_kg: 0.0,
            },
        })
    }

    fn finish(
        &self,
        mut prepared: PreparedDissolution,
        gas: &ClosedGasParcel,
        controls: DissolvedGasControls,
    ) -> Result<DissolutionExchange> {
        let after = self.sample(gas, prepared.exchange.state, controls, prepared.film)?;
        let next_mass = gas.species_mass_kg(&self.species)?;
        let mass_residual =
            finite((next_mass - prepared.initial_gas_mass_kg) + prepared.mass_into_liquid_kg)?;
        budget(
            "combined selected-species mass",
            mass_residual,
            controls.gas.mass_tolerance_kg,
        )?;
        let final_drive = after.equilibrium_concentration_mol_m3 - after.concentration_mol_m3;
        // A direction reversal larger than admitted inventory roundoff is not
        // accepted as an exact relaxation; no pressure or loading is clipped.
        let drive_roundoff =
            prepared.capacity * controls.gas.mass_tolerance_kg / self.molar_mass_kg_mol;
        if (prepared.drive > 0.0 && final_drive < -drive_roundoff)
            || (prepared.drive < 0.0 && final_drive > drive_roundoff)
        {
            return Err(invalid("finite Henry relaxation crossed its equilibrium"));
        }
        prepared.exchange.sources[1] = after;
        prepared.exchange.mass_residual_kg = mass_residual;
        Ok(prepared.exchange)
    }
}

struct PreparedDissolution {
    exchange: DissolutionExchange,
    initial_gas_mass_kg: f64,
    mass_into_liquid_kg: f64,
    capacity: f64,
    drive: f64,
    film: FilmResponse,
}
impl PreparedDissolution {
    fn transfer(&self, law: &HenryDissolution) -> SpeciesBoundaryRelativeTransfer {
        SpeciesBoundaryRelativeTransfer {
            species: law.species.clone(),
            mass_into_gas_kg: -self.mass_into_liquid_kg,
            relative_energy_into_gas_j: self.exchange.thermostat_heat_into_gas_j,
        }
    }
}

/// One species' original finite-film exchange within a shared gas reconstruction.
/// Gas storage and the combined energy residual belong to the enclosing step.
#[derive(Clone, Debug, PartialEq)]
pub struct DissolutionExchange {
    pub species: SpeciesId,
    pub initial: DissolvedGasInventory,
    pub state: DissolvedGasInventory,
    pub sources: [DissolutionSample; 2],
    /// Positive for absorption into the liquid, negative for exsolution.
    pub dissolved_moles_change: f64,
    pub relaxation_rate_per_s: f64,
    /// Liquid-film diffusion time divided by the common bulk relaxation time.
    pub film_relaxation_ratio: f64,
    /// The same ratio for a source-qualified gas film with known diffusivity
    /// and thickness; absent for a bare gas coefficient or liquid-only model.
    pub gas_film_relaxation_ratio: Option<f64>,
    pub partial_enthalpy_into_liquid_j: f64,
    pub dissolved_internal_energy_change_j: f64,
    pub thermostat_heat_into_gas_j: f64,
    pub thermostat_heat_into_liquid_j: f64,
    pub thermostat_heat_into_system_j: f64,
    pub mass_residual_kg: f64,
}

/// Only dissolved solute is stored here; the original finite gas remains its
/// own owner. Solvent mass/temperature are fixed and no solvent evaporates.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DissolvedGasInventory {
    liquid_volume_m3: f64,
    dissolved_moles: f64,
    time_s: f64,
}
impl DissolvedGasInventory {
    pub fn new(liquid_volume_m3: f64, dissolved_moles: f64, time_s: f64) -> Result<Self> {
        if !positive(liquid_volume_m3) || !nonnegative(dissolved_moles) || !nonnegative(time_s) {
            return Err(invalid(
                "positive liquid volume, nonnegative finite solute and clock required",
            ));
        }
        Ok(Self {
            liquid_volume_m3,
            dissolved_moles,
            time_s,
        })
    }
    pub const fn liquid_volume_m3(self) -> f64 {
        self.liquid_volume_m3
    }
    pub const fn dissolved_moles(self) -> f64 {
        self.dissolved_moles
    }
    pub const fn time_s(self) -> f64 {
        self.time_s
    }
}
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DissolutionFilm {
    /// Zero explicitly closes the existing path.
    pub area_m2: f64,
    pub thickness_m: f64,
    /// Fixed pressure-normalized gas-side coefficient [mol/(m² s Pa)].
    /// None preserves negligible gas-side resistance unless the law retains a
    /// source-qualified gas film; supplying both descriptions is ambiguous.
    /// Some(0) seals the path. This is not a velocity coefficient [m/s] or a
    /// volumetric kLa [1/s]. Its applicability throughout the isothermal
    /// pressure/composition change must be supplied independently; it carries
    /// no inferred material authority or gas-film diffusion time.
    /// Two simultaneously sealed films refuse an undetermined interface state.
    pub gas_side_transfer: Option<MolarPermeance>,
}
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DissolvedGasControls {
    pub gas: GasSpeciesBoundarySolve,
    pub maximum_step_s: f64,
    /// Applies to BOTH liquid bulk and actual liquid interface concentration.
    pub maximum_concentration_mol_m3: f64,
    /// Authored ceiling in (0,1) on lambda*delta²/D for the liquid film and any
    /// source-qualified gas film. A bare gas coefficient supplies no diffusion
    /// time. Not an accuracy theorem; a sealed path has no bulk relaxation.
    pub maximum_film_relaxation_ratio: f64,
}
#[derive(Clone, Debug, PartialEq)]
pub struct DissolutionSample {
    pub partial_pressure_pa: f64,
    pub concentration_mol_m3: f64,
    /// Henry concentration implied by the BULK gas, not the actual interface.
    pub equilibrium_concentration_mol_m3: f64,
    /// Henry-equivalent interfacial partial pressure [Pa]. For a sealed gas
    /// film it is a disconnected equilibrium potential, not a bulk gas state.
    pub interface_partial_pressure_pa: f64,
    /// Liquid-side concentration in Henry equilibrium with that interface [mol/m³].
    pub interface_concentration_mol_m3: f64,
    /// Signed common quasi-steady film flux, positive into the liquid.
    pub molar_flux_mol_m2_s: f64,
    /// Original claim bundle at [bulk liquid, liquid-side interface].
    pub material: [ResolvedMaterialStatePoint; 2],
    /// Original gas diffusivity source at [bulk gas, gas-side interface].
    pub gas_film_material: Option<[ResolvedMaterialStatePoint; 2]>,
}
#[derive(Clone, Debug, PartialEq)]
pub struct DissolutionStep {
    pub initial: DissolvedGasInventory,
    pub state: DissolvedGasInventory,
    pub gas: GasSpeciesBoundaryStep,
    pub sources: [DissolutionSample; 2],
    /// Positive for dissolution, negative for exsolution.
    pub dissolved_moles_change: f64,
    pub relaxation_rate_per_s: f64,
    pub film_relaxation_ratio: f64,
    /// Gas-film diffusion/bulk relaxation ratio when its original diffusivity
    /// source and physical thickness are supplied, otherwise absent.
    pub gas_film_relaxation_ratio: Option<f64>,
    /// Signed original gas partial enthalpy entering the liquid, including flow work.
    pub partial_enthalpy_into_liquid_j: f64,
    pub dissolved_internal_energy_change_j: f64,
    /// Heat required to keep the finite GAS at its original T.
    pub thermostat_heat_into_gas_j: f64,
    /// Heat required to keep the finite SOLUTION at that T, including dissolution.
    pub thermostat_heat_into_liquid_j: f64,
    /// Required external heat; this is not an adiabatic evolution.
    pub thermostat_heat_into_system_j: f64,
    pub mass_residual_kg: f64,
    pub energy_residual_j: f64,
}

#[cfg(test)]
#[path = "dissolved_gas/tests.rs"]
mod tests;
