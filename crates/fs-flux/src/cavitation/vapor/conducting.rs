//! Resolved liquid conduction and Stefan transport on the native vapor clock.
//!
//! The spherical liquid grid moves with the geometric bubble surface while
//! each shell retains its volume. Liquid flows relative to that grid during
//! phase exchange. The original field operator carries this same signed mass
//! flux and its source enthalpy through every face. Its interface Fourier heat
//! drives the existing saturated-vapor balance at EACH adaptive stage.
//!
//! Radius, liquid interface velocity, vapor temperature, all liquid-cell
//! temperatures and their physical quadratures use one fixed-stack DP5 solve.
//! There is no lagged heat callback, independently advanced liquid clock or
//! reset thermal layer. Actual radial viscous loss in the resolved shells is
//! returned to their heat inventory; loss outside the finite grid is reported
//! by the field operator without being deposited in these cells.
//!
//! The common fixed enthalpy reference gives the discrete thermal account
//! `Delta(U_v+U_g-h_ref m_v)+Delta H_liquid+W_contents+Q_outer+H_outer-Q_visc=0`.
//! Outer conductive heat and advected enthalpy are positive OUT of the grid.
//! They connect to the prescribed-temperature exterior reservoir; the finite
//! shell domain is not a sealed finite-liquid vessel. W_contents integrates
//! (p_v+p_g)dV and remains in the same mechanical/thermal ownership. The optional
//! conserved dry gas shares the native temperature and volume; its original
//! energy increment is independently checked against integral m_g Cv Tdot.
//!
//! The field uses its explicitly bounded saturation-branch liquid caloric
//! approximation, with pressure-path, property-drift and pressure-heating
//! screens. The original vapor/local-equilibrium and incompressible-liquid
//! assumptions still apply. Spatial resolution, outer-domain extent and
//! temporal refinement remain necessary; a small energy residual alone does
//! not establish an accurate diffusion layer or bubble trajectory.

use super::liquid::{LiquidThermalRates, MAX_LIQUID_THERMAL_CELLS, SphericalLiquidField};
use super::{
    BubbleContentsLoad, SaturatedVaporBubble, VaporBubbleError, VaporBubbleNumerics,
    VaporBubbleObservation, VaporBubbleState, VaporBubbleStep, stage_derivatives,
    validate_numerics,
};
use crate::cavitation::{CavitationError, positive, stepper};
#[cfg(feature = "caloric-cavitation")]
use fs_material::cavitation::vapor::gas::ConservedBubbleGas;
#[cfg(feature = "caloric-cavitation")]
use fs_material::thermochemical_gas::energy::GasEnergyError;

/// Conserved noncondensable inventory and explicit wet-mixture validity bounds.
#[cfg(feature = "caloric-cavitation")]
#[path = "conducting/gas.rs"]
pub mod gas;

const PHASE_COORDINATES: usize = 8;
const FIELD_END: usize = PHASE_COORDINATES + MAX_LIQUID_THERMAL_CELLS;
const OUTER_HEAT: usize = FIELD_END;
const OUTER_ENTHALPY: usize = FIELD_END + 1;
const VISCOUS_HEAT: usize = FIELD_END + 2;
const GAS_ENERGY: usize = FIELD_END + 3;
const GAS_WORK: usize = FIELD_END + 4;
const COORDINATES: usize = FIELD_END + 5;

/// Local error scales and independent whole-window thermal balances.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ConductingVaporBubbleNumerics {
    pub bubble: VaporBubbleNumerics,
    /// Relative control applies to each liquid cell's temperature CHANGE
    /// during the requested window, so an absolute Kelvin datum cannot hide
    /// errors in a thin, weakly superheated or subcooled boundary layer.
    pub liquid_temperature_absolute_tolerance_k: f64,
    pub liquid_energy_balance_tolerance_j: f64,
    pub total_energy_balance_tolerance_j: f64,
}

#[derive(Clone, Debug, PartialEq)]
pub enum ConductingVaporBubbleError {
    Physics(VaporBubbleError),
    #[cfg(feature = "caloric-cavitation")]
    Gas(GasEnergyError),
    InvalidInput(&'static str),
    OutsideDomain {
        quantity: &'static str,
        value: f64,
        minimum: f64,
        maximum: f64,
    },
    SourceMismatch,
    Balance {
        quantity: &'static str,
        residual_j: f64,
        tolerance_j: f64,
    },
    Cancelled,
}
impl core::fmt::Display for ConductingVaporBubbleError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "resolved-liquid vapor bubble: {self:?}")
    }
}
impl std::error::Error for ConductingVaporBubbleError {}
impl From<VaporBubbleError> for ConductingVaporBubbleError {
    fn from(error: VaporBubbleError) -> Self {
        match error {
            VaporBubbleError::Cancelled => Self::Cancelled,
            other => Self::Physics(other),
        }
    }
}
impl From<CavitationError> for ConductingVaporBubbleError {
    fn from(error: CavitationError) -> Self {
        VaporBubbleError::from(error).into()
    }
}

#[cfg(feature = "caloric-cavitation")]
impl From<GasEnergyError> for ConductingVaporBubbleError {
    fn from(error: GasEnergyError) -> Self {
        if error == GasEnergyError::Cancelled {
            Self::Cancelled
        } else {
            Self::Gas(error)
        }
    }
}

/// Accepted/stage gas contribution and explicitly estimated mixture limits.
/// The temperature and volume are the SAME ones in the native bubble state.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct GasVaporObservation {
    pub mass_kg: f64,
    pub density_kg_m3: f64,
    pub partial_pressure_pa: f64,
    pub cv_j_kg_k: f64,
    pub energy_rate_w: f64,
    pub pressure_work_power_w: f64,
    /// Heat INTO the conserved dry component, m_g Cv Tdot + p_g Vdot.
    pub heat_power_w: f64,
    pub mixture_density_kg_m3: f64,
    /// Frozen-composition speed of the declared additive mixture approximation.
    pub mixture_sound_speed_m_s: f64,
    pub vapor_mass_fraction: f64,
    /// Fourier power entering the gas-vapor interior: Q_liquid-L*mdot_v.
    pub internal_heat_power_w: f64,
    pub estimated_internal_temperature_difference_k: f64,
    pub relative_composition_gradient_scale: f64,
    pub thermal_response_ratio: f64,
    pub diffusion_response_ratio: f64,
    pub pressure_relaxation_ratio: f64,
}

/// Native physical endpoint and the same resolved liquid thermal response.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ConductingVaporBubbleObservation {
    pub bubble: VaporBubbleObservation,
    pub liquid: LiquidThermalRates,
    pub gas: Option<GasVaporObservation>,
}

/// One admitted window of radial, phase and liquid thermal evolution.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ConductingVaporBubbleStep {
    pub bubble: VaporBubbleStep,
    pub liquid: LiquidThermalRates,
    pub gas: Option<GasVaporObservation>,
    /// Independent integral m_g Cv(T) Tdot dt, zero without gas [J].
    pub integrated_noncondensable_energy_j: f64,
    /// Integral p_g Vdot dt [J]. Already INCLUDED in bubble.pressure_work_j.
    pub noncondensable_pressure_work_j: f64,
    /// Integral (m_g Cv Tdot+p_g Vdot) dt [J], heat into the dry component.
    pub noncondensable_heat_j: f64,
    /// Original m_g Delta u_g minus its independently integrated coordinate.
    pub noncondensable_energy_residual_j: f64,
    /// Actual source enthalpy change in the fixed-volume liquid shells [J].
    pub liquid_energy_change_j: f64,
    /// Fourier heat leaving the finite thermal grid into the exterior bath [J].
    pub outer_conductive_heat_j: f64,
    /// Enthalpy leaving the moving grid, relative to its fixed reference [J].
    pub outer_advected_enthalpy_j: f64,
    /// Original radial viscous dissipation inside the resolved shells [J].
    pub viscous_heat_j: f64,
    /// Delta H_liq + Q_interface + H_interface + Q_outer + H_outer - Q_visc.
    pub liquid_energy_residual_j: f64,
    /// Delta(U_v+U_g-h_ref m_v) + Delta H_liq + W_contents + Q_outer + H_outer - Q_visc.
    pub total_energy_residual_j: f64,
}

/// One authoritative radial clock and retained liquid temperature profile.
pub struct ConductingVaporBubble {
    model: SaturatedVaporBubble,
    liquid: SphericalLiquidField,
    temperatures_k: [f64; MAX_LIQUID_THERMAL_CELLS],
    #[cfg(feature = "caloric-cavitation")]
    gas: Option<gas::GasBinding>,
}

/// Complete uncommitted physical window. Observation may still refuse it.
#[must_use = "inspect and commit to publish the complete phase and liquid field"]
pub struct ConductingVaporBubbleProposal<'a> {
    target: &'a mut ConductingVaporBubble,
    step: ConductingVaporBubbleStep,
    temperatures_k: [f64; MAX_LIQUID_THERMAL_CELLS],
}
impl ConductingVaporBubbleProposal<'_> {
    pub const fn step(&self) -> &ConductingVaporBubbleStep {
        &self.step
    }
    pub fn liquid_temperatures_k(&self) -> &[f64] {
        &self.temperatures_k[..self.target.liquid.grid().cell_count()]
    }
    pub fn commit(self) -> ConductingVaporBubbleStep {
        self.target.model.state = self.step.bubble.observation.state;
        self.target.temperatures_k = self.temperatures_k;
        self.step
    }
}

fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), ConductingVaporBubbleError> {
    if cancelled() {
        Err(ConductingVaporBubbleError::Cancelled)
    } else {
        Ok(())
    }
}
fn budget(
    quantity: &'static str,
    residual_j: f64,
    tolerance_j: f64,
) -> Result<(), ConductingVaporBubbleError> {
    if !residual_j.is_finite() || residual_j.abs() > tolerance_j {
        Err(ConductingVaporBubbleError::Balance {
            quantity,
            residual_j,
            tolerance_j,
        })
    } else {
        Ok(())
    }
}

impl ConductingVaporBubble {
    /// Retain an explicitly sampled liquid temperature for every prepared
    /// shell. Both owners must use the same original water and interface data.
    /// The grid's reference radii define fixed shell volumes; its actual radii
    /// are evaluated at the native bubble's current radius, without remapping.
    pub fn new(
        model: SaturatedVaporBubble,
        liquid: SphericalLiquidField,
        temperatures_k: &[f64],
    ) -> Result<Self, ConductingVaporBubbleError> {
        if model.source() != liquid.source() {
            return Err(ConductingVaporBubbleError::SourceMismatch);
        }
        if temperatures_k.len() != liquid.grid().cell_count() {
            return Err(ConductingVaporBubbleError::InvalidInput(
                "one liquid temperature per prepared shell",
            ));
        }
        let mut temperatures = [0.0; MAX_LIQUID_THERMAL_CELLS];
        temperatures[..temperatures_k.len()].copy_from_slice(temperatures_k);
        let (thermal, _) = model.thermal_state(model.state())?;
        // Validate only geometry/source/profile here. Actual pressure, bath
        // forcing and phase-flow admission belong to observe/advance.
        liquid.prepare(
            &thermal,
            &temperatures,
            liquid.source().state().temperature_k(),
        )?;
        Ok(Self {
            model,
            liquid,
            temperatures_k: temperatures,
            #[cfg(feature = "caloric-cavitation")]
            gas: None,
        })
    }

    /// Add an explicitly initialized noncondensable inventory at the SAME
    /// accepted T and spherical volume. Existing inventory cannot be replaced.
    /// No heat, mass or pressure is created to force an equilibrium radius.
    #[cfg(feature = "caloric-cavitation")]
    pub fn with_noncondensable_gas(
        mut self,
        inventory: ConservedBubbleGas,
        limits: gas::GasVaporLimits,
    ) -> Result<Self, ConductingVaporBubbleError> {
        if self.gas.is_some() {
            return Err(ConductingVaporBubbleError::InvalidInput(
                "a conserved noncondensable inventory is already bound",
            ));
        }
        let (thermal, _) = self.model.thermal_state(self.state())?;
        self.gas = Some(gas::GasBinding::new(inventory, limits, &thermal)?);
        Ok(self)
    }

    /// Immutable original inventory/source; its reference parcel is NOT an
    /// independently evolved gas state. Current T/V belong to Self::state.
    #[cfg(feature = "caloric-cavitation")]
    pub fn gas(&self) -> Option<&ConservedBubbleGas> {
        self.gas.as_ref().map(|binding| &binding.inventory)
    }

    /// Original water-phase owner. Use Self::observe for the composed dynamics
    /// when a conserved gas inventory is present.
    pub const fn model(&self) -> &SaturatedVaporBubble {
        &self.model
    }
    pub const fn liquid(&self) -> &SphericalLiquidField {
        &self.liquid
    }
    pub const fn state(&self) -> VaporBubbleState {
        self.model.state()
    }
    pub fn liquid_temperatures_k(&self) -> &[f64] {
        &self.temperatures_k[..self.liquid.grid().cell_count()]
    }

    /// Read the native dynamics with heat from the CURRENT resolved profile.
    /// The far temperature defines the exterior conductive and advective bath.
    pub fn observe(
        &self,
        pressure_pa: f64,
        bath_temperature_k: f64,
    ) -> Result<ConductingVaporBubbleObservation, ConductingVaporBubbleError> {
        self.evaluate(
            self.state(),
            &self.temperatures_k,
            pressure_pa,
            bath_temperature_k,
        )
    }

    fn evaluate(
        &self,
        state: VaporBubbleState,
        temperatures_k: &[f64; MAX_LIQUID_THERMAL_CELLS],
        pressure_pa: f64,
        bath_temperature_k: f64,
    ) -> Result<ConductingVaporBubbleObservation, ConductingVaporBubbleError> {
        let (thermal, variation) = self.model.thermal_state(state)?;
        let stage = self
            .liquid
            .prepare(&thermal, temperatures_k, bath_temperature_k)?;
        #[cfg(feature = "caloric-cavitation")]
        let sample = self
            .gas
            .as_ref()
            .map(|binding| binding.sample(state.temperature_k(), thermal.volume_m3))
            .transpose()?;
        #[cfg(feature = "caloric-cavitation")]
        let load = sample
            .as_ref()
            .map_or_else(BubbleContentsLoad::default, |sample| BubbleContentsLoad {
                gas_pressure_pa: sample.pressure_pa,
                gas_heat_capacity_j_k: sample.heat_capacity_j_k,
            });
        #[cfg(not(feature = "caloric-cavitation"))]
        let load = BubbleContentsLoad::default();
        let bubble = self.model.observe_thermal_loaded(
            thermal,
            variation,
            pressure_pa,
            stage.interface_heat_power_w(),
            load,
        )?;
        let liquid = self.liquid.finish(stage, &bubble, pressure_pa)?;
        #[cfg(feature = "caloric-cavitation")]
        let gas = match (self.gas.as_ref(), sample) {
            (Some(binding), Some(sample)) => Some(binding.observe(
                sample,
                &bubble,
                self.model.source().state().liquid_density_kg_m3(),
            )?),
            _ => None,
        };
        #[cfg(not(feature = "caloric-cavitation"))]
        let gas = None;
        Ok(ConductingVaporBubbleObservation {
            bubble,
            liquid,
            gas,
        })
    }

    /// One fully coupled physical advance. Pressure and bath-temperature laws
    /// must be deterministic and resolved by the supplied maximum time step.
    pub fn advance_controlled<P, B, C>(
        &mut self,
        end_time_s: f64,
        numerics: &ConductingVaporBubbleNumerics,
        pressure: &P,
        bath_temperature: &B,
        cancelled: &mut C,
    ) -> Result<ConductingVaporBubbleStep, ConductingVaporBubbleError>
    where
        P: Fn(f64) -> f64,
        B: Fn(f64) -> f64,
        C: FnMut() -> bool,
    {
        Ok(self
            .prepare_advance_controlled(
                end_time_s,
                numerics,
                pressure,
                bath_temperature,
                cancelled,
            )?
            .commit())
    }

    pub fn prepare_advance_controlled<'a, P, B, C>(
        &'a mut self,
        end_time_s: f64,
        numerics: &ConductingVaporBubbleNumerics,
        pressure: &P,
        bath_temperature: &B,
        cancelled: &mut C,
    ) -> Result<ConductingVaporBubbleProposal<'a>, ConductingVaporBubbleError>
    where
        P: Fn(f64) -> f64,
        B: Fn(f64) -> f64,
        C: FnMut() -> bool,
    {
        poll(cancelled)?;
        validate_numerics(&numerics.bubble)?;
        if ![
            numerics.liquid_temperature_absolute_tolerance_k,
            numerics.liquid_energy_balance_tolerance_j,
            numerics.total_energy_balance_tolerance_j,
        ]
        .into_iter()
        .all(positive)
        {
            return Err(ConductingVaporBubbleError::InvalidInput(
                "positive liquid-temperature and independent thermal energy budgets",
            ));
        }
        let start = self.state();
        let initial = self.observe(pressure(start.time_s()), bath_temperature(start.time_s()))?;
        let (initial_phase, _) = self.model.thermal_state(start)?;
        let anchor = self.liquid.reference_enthalpy_j_kg();
        let mut values = [0.0; COORDINATES];
        values[0] = start.radius_m();
        values[1] = start.liquid_velocity_m_s();
        values[2] = start.temperature_k();
        // Integrate temperature increments around the retained accepted
        // profile. The shared relative-error controller therefore measures
        // resolved changes instead of the arbitrary absolute Kelvin datum.
        let mut absolute = [numerics.liquid_temperature_absolute_tolerance_k; COORDINATES];
        absolute[..PHASE_COORDINATES].copy_from_slice(&[
            numerics.bubble.radial.radius_absolute_tolerance_m,
            numerics.bubble.radial.velocity_absolute_tolerance_m_s,
            numerics.bubble.temperature_absolute_tolerance_k,
            numerics.bubble.energy_absolute_tolerance_j,
            numerics.bubble.energy_absolute_tolerance_j,
            numerics.bubble.energy_absolute_tolerance_j,
            numerics.bubble.mass_absolute_tolerance_kg,
            numerics.bubble.energy_absolute_tolerance_j,
        ]);
        absolute[OUTER_HEAT..].fill(numerics.bubble.energy_absolute_tolerance_j);
        let result = stepper::advance_array(
            start.time_s(),
            values,
            end_time_s,
            &numerics.bubble.radial,
            absolute,
            &mut |time, values: &[f64; COORDINATES]| {
                let state = VaporBubbleState::try_new(time, values[0], values[1], values[2])?;
                let mut temperatures = [0.0; MAX_LIQUID_THERMAL_CELLS];
                for (index, temperature) in temperatures.iter_mut().enumerate() {
                    *temperature = self.temperatures_k[index] + values[PHASE_COORDINATES + index];
                }
                let observation =
                    self.evaluate(state, &temperatures, pressure(time), bath_temperature(time))?;
                let mut rates = [0.0; COORDINATES];
                rates[..PHASE_COORDINATES]
                    .copy_from_slice(&stage_derivatives(&observation.bubble, anchor));
                rates[PHASE_COORDINATES..FIELD_END]
                    .copy_from_slice(&observation.liquid.temperature_rate_k_s);
                rates[OUTER_HEAT] = observation.liquid.outer_conductive_power_w;
                rates[OUTER_ENTHALPY] = observation.liquid.outer_advected_power_w;
                rates[VISCOUS_HEAT] = observation.liquid.viscous_heat_power_w;
                rates[GAS_ENERGY] = observation.bubble.noncondensable_heat_capacity_j_k
                    * observation.bubble.temperature_rate_k_s;
                rates[GAS_WORK] = observation.bubble.noncondensable_pressure_pa
                    * observation.bubble.volume_rate_m3_s;
                Ok(rates)
            },
            cancelled,
            &|error: &ConductingVaporBubbleError| {
                matches!(
                    error,
                    ConductingVaporBubbleError::OutsideDomain { .. }
                        | ConductingVaporBubbleError::Physics(
                            VaporBubbleError::OutsideModel { .. }
                                | VaporBubbleError::TemperatureOutsideDomain { .. }
                                | VaporBubbleError::Dynamics(
                                    CavitationError::RadiusOutsideDomain { .. }
                                        | CavitationError::NonFiniteDerived(_)
                                )
                        )
                )
            },
        )?;
        poll(cancelled)?;
        let next = VaporBubbleState::try_new(
            result.time_s,
            result.values[0],
            result.values[1],
            result.values[2],
        )?;
        let mut temperatures = [0.0; MAX_LIQUID_THERMAL_CELLS];
        for (index, temperature) in temperatures.iter_mut().enumerate() {
            *temperature = self.temperatures_k[index] + result.values[PHASE_COORDINATES + index];
        }
        let endpoint = self.evaluate(
            next,
            &temperatures,
            pressure(result.time_s),
            bath_temperature(result.time_s),
        )?;
        let mut phase_values = [0.0; PHASE_COORDINATES];
        phase_values.copy_from_slice(&result.values[..PHASE_COORDINATES]);
        let phase_result = stepper::ArrayStep {
            time_s: result.time_s,
            values: phase_values,
            accepted_steps: result.accepted_steps,
            rejected_steps: result.rejected_steps,
            rhs_evaluations: result.rhs_evaluations,
            maximum_normalized_local_error: result.maximum_normalized_local_error,
        };
        #[cfg(feature = "caloric-cavitation")]
        let gas_energy_change = self
            .gas
            .as_ref()
            .map(|binding| {
                binding
                    .inventory
                    .energy_change_j(start.temperature_k(), next.temperature_k())
            })
            .transpose()?
            .unwrap_or(0.0);
        #[cfg(not(feature = "caloric-cavitation"))]
        let gas_energy_change = 0.0;
        let gas_energy_residual = gas_energy_change - result.values[GAS_ENERGY];
        budget(
            "original noncondensable caloric source",
            gas_energy_residual,
            numerics.bubble.energy_balance_tolerance_j,
        )?;
        let bubble = self.model.finish_step(
            initial_phase,
            anchor,
            &phase_result,
            endpoint.bubble,
            &numerics.bubble,
            gas_energy_change,
        )?;
        let liquid_change = endpoint.liquid.energy_change_from(&initial.liquid)?;
        let outer_heat = result.values[OUTER_HEAT];
        let outer_enthalpy = result.values[OUTER_ENTHALPY];
        let viscous = result.values[VISCOUS_HEAT];
        let liquid_residual = liquid_change
            + bubble.supplied_heat_j
            + bubble.reference_invariant_transferred_enthalpy_j
            + outer_heat
            + outer_enthalpy
            - viscous;
        let total_residual = liquid_change
            + bubble.reference_invariant_energy_change_j
            + bubble.pressure_work_j
            + outer_heat
            + outer_enthalpy
            - viscous;
        budget(
            "resolved liquid first law",
            liquid_residual,
            numerics.liquid_energy_balance_tolerance_j,
        )?;
        budget(
            "combined phase/liquid first law",
            total_residual,
            numerics.total_energy_balance_tolerance_j,
        )?;
        if !viscous.is_finite() || viscous < 0.0 {
            return Err(ConductingVaporBubbleError::InvalidInput(
                "integrated resolved viscous heat",
            ));
        }
        // Only qualified complete candidates may publish. The retained gas
        // reference remains immutable; current T/V have one authoritative owner.
        #[cfg(feature = "caloric-cavitation")]
        if let Some(binding) = &self.gas {
            poll(cancelled)?;
            binding
                .inventory
                .qualify(next.temperature_k(), endpoint.bubble.volume_m3)?;
        }
        poll(cancelled)?;
        let gas_heat = result.values[GAS_ENERGY] + result.values[GAS_WORK];
        if !gas_heat.is_finite() {
            return Err(ConductingVaporBubbleError::InvalidInput(
                "represented gas heat integral",
            ));
        }
        let step = ConductingVaporBubbleStep {
            bubble,
            liquid: endpoint.liquid,
            gas: endpoint.gas,
            integrated_noncondensable_energy_j: result.values[GAS_ENERGY],
            noncondensable_pressure_work_j: result.values[GAS_WORK],
            noncondensable_heat_j: gas_heat,
            noncondensable_energy_residual_j: gas_energy_residual,
            liquid_energy_change_j: liquid_change,
            outer_conductive_heat_j: outer_heat,
            outer_advected_enthalpy_j: outer_enthalpy,
            viscous_heat_j: viscous,
            liquid_energy_residual_j: liquid_residual,
            total_energy_residual_j: total_residual,
        };
        Ok(ConductingVaporBubbleProposal {
            target: self,
            step,
            temperatures_k: temperatures,
        })
    }
}

#[cfg(test)]
#[path = "conducting/tests.rs"]
mod tests;
