//! Equilibrium miscibility energy on the shared conservative thermal timeline.
//!
//! The fixed component inventory defines composition and a mass normalization;
//! timeline cells may have different masses of that same composition. Closed
//! storage uses zero specific enthalpy at its cold bracket endpoint. Explicit
//! stream storage instead retains the original common pure-species datum, so
//! mass exchange between compositions carries their actual mixing enthalpy.
//! Every temperature inversion uses the original caloric equilibrium solver,
//! including the energy
//! of phase redistribution, without a sampled enthalpy chart or another Cp law.
//!
//! This is thermal storage at externally maintained pressure. It does not move
//! species, infer interphase conductances, or change reference cell geometry.
//! Numeric thermal trials do not refresh mixture transport: publish resolved
//! phases through the original source-bound law when those properties are used.
//! Entropy uses the same pure Cp/T integrals and accepted mixing phases.
//! Only the explicit stream constructor supplies an advection binding and
//! caloric identity. Chemical entropy transported between different compositions
//! remains outside the existing stationary-species flow entropy capability.

use super::{
    BinaryCaloricSolveV1, BinaryCaloricStateV1, BinaryLiquidCaloricError, BinaryLiquidCaloricLaw,
    IsobaricHeatCapacityLaw, SourceOwner, poll,
};
use fs_qty::Temperature;
use fs_thermal::enthalpy::{DepositionError, ThermalLaw};
use std::sync::Arc;

/// Equilibrium storage and explicitly declared bulk conductivity on the native
/// implicit thermal solver, including phase-redistribution heat capacity.
pub mod conduction;

/// Common-datum liquid storage for the existing stationary-species flow graph.
pub mod stream;

/// Immutable fixed-composition caloric law, retaining both original Cp sources.
#[derive(Clone, Debug, PartialEq)]
pub struct BinaryLiquidThermalStorage {
    owner: Arc<SourceOwner>,
    endpoints: [BinaryCaloricStateV1; 2],
    mass_kg: f64,
    lower_specific_h_j_kg: f64,
    upper_specific_h_j_kg: f64,
    maximum_slope: f64,
    controls: BinaryCaloricSolveV1,
    stream: Option<stream::StreamBinding>,
}

impl BinaryLiquidCaloricLaw<'_, '_> {
    /// Prepare one immutable thermal domain and prove a reusable global CFL
    /// bound. Both endpoint phases must admit the original material sources.
    /// The original caloric owner bounds dH/dT from its source knots, including
    /// a signed heat capacity of mixing. Equilibrium redistribution adds a
    /// nonnegative response, so dT/dh <= mass / minimum_capacity. The bound
    /// covers every source segment and both supplied thermal endpoints.
    pub fn thermal_storage(
        &self,
        component_moles: [f64; 2],
        controls: BinaryCaloricSolveV1,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<BinaryLiquidThermalStorage, BinaryLiquidCaloricError> {
        poll(&mut cancelled)?;
        let cold = self.equilibrate(
            controls.temperature_bracket[0],
            component_moles,
            controls.coexistence,
            &mut cancelled,
        )?;
        // Reuse the thermal owner's full control/domain admission, including
        // finite iteration and convergence budgets, rather than duplicating it.
        self.owner
            .model
            .add_heat(cold.caloric_state(), 0.0, controls, &mut cancelled)?;
        let hot = self.equilibrate(
            controls.temperature_bracket[1],
            component_moles,
            controls.coexistence,
            &mut cancelled,
        )?;
        let mass_kg = cold.equilibrium().total_mass_kg();
        let span = self
            .owner
            .model
            .enthalpy_change_j(cold.caloric_state(), hot.caloric_state())?;
        let upper_specific_h_j_kg = span / mass_kg;
        let minimum_capacity = self
            .owner
            .model
            .isobaric_heat_capacity_lower_bound_j_per_k(
                component_moles,
                controls.temperature_bracket,
            )?;
        // Products and sum round downward; the reciprocal rounds upward. This
        // also handles subnormal inventories without a relative-epsilon guess.
        let maximum_slope = (mass_kg / minimum_capacity).next_up();
        if !positive(upper_specific_h_j_kg)
            || !positive(minimum_capacity)
            || !positive(maximum_slope)
        {
            return Err(BinaryLiquidCaloricError::Invalid(
                "unrepresentable thermal-storage enthalpy span or heat-capacity bound",
            ));
        }
        poll(&mut cancelled)?;
        Ok(BinaryLiquidThermalStorage {
            owner: Arc::clone(&self.owner),
            endpoints: [cold.caloric_state().clone(), hot.caloric_state().clone()],
            mass_kg,
            lower_specific_h_j_kg: 0.0,
            upper_specific_h_j_kg,
            maximum_slope,
            controls,
            stream: None,
        })
    }
}

impl BinaryLiquidThermalStorage {
    pub const fn temperature_domain(&self) -> [Temperature; 2] {
        self.controls.temperature_bracket
    }
    /// Inclusive specific-enthalpy bounds [J/kg]. Ordinary thermal storage
    /// uses zero at the cold endpoint; stream storage retains the original
    /// common pure-species datum, including the equilibrium mixing enthalpy.
    pub const fn specific_enthalpy_domain_j_kg(&self) -> [f64; 2] {
        [self.lower_specific_h_j_kg, self.upper_specific_h_j_kg]
    }
    pub fn sources(&self) -> &[IsobaricHeatCapacityLaw; 2] {
        &self.owner.sources
    }

    /// Equilibrium specific heat at the original fixed pressure and overall
    /// composition [J/(kg K)]. The phase amounts and compositions may change
    /// with temperature; this includes their enthalpy response as well as the
    /// two original pure-component Cp sources. It is independent of the amount
    /// of material chosen to normalize this storage law.
    ///
    /// A phase-onset corner has no unique two-sided Cp and refuses explicitly.
    /// The implicit conduction adapter separately selects a documented
    /// one-sided Newton tangent at those corners.
    pub fn specific_heat_cp_j_kg_k(
        &self,
        temperature: Temperature,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<f64, BinaryLiquidCaloricError> {
        poll(&mut cancelled)?;
        let [lower, upper] = self.temperature_domain();
        if !temperature.value().is_finite()
            || temperature.value() < lower.value()
            || temperature.value() > upper.value()
        {
            return Err(BinaryLiquidCaloricError::Invalid(
                "temperature outside thermal-storage domain",
            ));
        }
        let state = self.owner.model.equilibrate(
            temperature,
            self.endpoints[0].component_moles(),
            self.controls.coexistence,
            &mut cancelled,
        )?;
        let cp = state.isobaric_heat_capacity_j_per_k()? / self.mass_kg;
        if !positive(cp) {
            return Err(BinaryLiquidCaloricError::Invalid(
                "unrepresentable equilibrium specific heat",
            ));
        }
        poll(&mut cancelled)?;
        Ok(cp)
    }

    /// Query the original equilibrium calorics; never extrapolate the prepared
    /// thermal domain or subtract two large datum-dependent total enthalpies.
    pub fn specific_enthalpy_j_kg(
        &self,
        temperature: Temperature,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<f64, BinaryLiquidCaloricError> {
        poll(&mut cancelled)?;
        let [lower, upper] = self.temperature_domain();
        if !temperature.value().is_finite()
            || temperature.value() < lower.value()
            || temperature.value() > upper.value()
        {
            return Err(BinaryLiquidCaloricError::Invalid(
                "temperature outside thermal-storage domain",
            ));
        }
        if temperature == lower {
            return Ok(self.lower_specific_h_j_kg);
        }
        if temperature == upper {
            return Ok(self.upper_specific_h_j_kg);
        }
        let state = self.owner.model.equilibrate(
            temperature,
            self.endpoints[0].component_moles(),
            self.controls.coexistence,
            &mut cancelled,
        )?;
        let increment = self
            .owner
            .model
            .enthalpy_change_j(&self.endpoints[0], &state)?
            / self.mass_kg;
        if !positive(increment) {
            return Err(BinaryLiquidCaloricError::Invalid(
                "positive thermal-storage enthalpy increment is unrepresentable",
            ));
        }
        let h = self.lower_specific_h_j_kg + increment;
        if self.stream.is_some()
            && (h <= self.lower_specific_h_j_kg || h >= self.upper_specific_h_j_kg)
        {
            return Err(BinaryLiquidCaloricError::Invalid(
                "interior stream enthalpy is not representable in the common datum",
            ));
        }
        self.admit(h)?;
        poll(&mut cancelled)?;
        Ok(h)
    }

    /// Recover the same equilibrium phase selection as the closed heat solver.
    pub fn caloric_state_at_specific_enthalpy(
        &self,
        specific_h_j_kg: f64,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<BinaryCaloricStateV1, BinaryLiquidCaloricError> {
        poll(&mut cancelled)?;
        self.admit(specific_h_j_kg)?;
        // Preserve exact endpoints through the J/kg-to-J round trip.
        if specific_h_j_kg == self.lower_specific_h_j_kg {
            return Ok(self.endpoints[0].clone());
        }
        if specific_h_j_kg == self.upper_specific_h_j_kg {
            return Ok(self.endpoints[1].clone());
        }
        let heat_j = (specific_h_j_kg - self.lower_specific_h_j_kg) * self.mass_kg;
        if !positive(heat_j) {
            return Err(BinaryLiquidCaloricError::Invalid(
                "positive thermal-storage heat input is unrepresentable",
            ));
        }
        let balance =
            self.owner
                .model
                .add_heat(&self.endpoints[0], heat_j, self.controls, &mut cancelled)?;
        Ok(balance.final_state().clone())
    }

    fn admit(&self, h: f64) -> Result<(), BinaryLiquidCaloricError> {
        if !h.is_finite() || h < self.lower_specific_h_j_kg || h > self.upper_specific_h_j_kg {
            Err(BinaryLiquidCaloricError::Invalid(
                "enthalpy outside thermal-storage domain",
            ))
        } else {
            Ok(())
        }
    }
}

impl ThermalLaw for BinaryLiquidThermalStorage {
    fn caloric_identity(&self) -> Option<[u8; 32]> {
        self.stream
            .as_ref()
            .map(stream::StreamBinding::caloric_identity)
    }
    fn advected_composition(
        &self,
    ) -> Option<fs_thermal::enthalpy::advection::species::CaloricComposition<'_>> {
        self.stream.as_ref().map(stream::StreamBinding::composition)
    }
    fn temperature_k(&self, h: f64) -> Result<f64, DepositionError> {
        self.caloric_state_at_specific_enthalpy(h, || false)
            .map(|state| state.temperature().value())
            .map_err(|_| DepositionError::MaterialDomain)
    }
    fn maximum_temperature_slope(&self) -> f64 {
        self.maximum_slope
    }
    fn specific_entropy_change_j_kg_k(
        &self,
        initial_h_j_kg: f64,
        final_h_j_kg: f64,
    ) -> Result<Option<f64>, DepositionError> {
        self.admit(initial_h_j_kg)
            .and_then(|()| self.admit(final_h_j_kg))
            .map_err(|_| DepositionError::MaterialDomain)?;
        if initial_h_j_kg == final_h_j_kg {
            return Ok(Some(0.0));
        }
        let initial = self
            .caloric_state_at_specific_enthalpy(initial_h_j_kg, || false)
            .map_err(|_| DepositionError::MaterialDomain)?;
        let final_state = self
            .caloric_state_at_specific_enthalpy(final_h_j_kg, || false)
            .map_err(|_| DepositionError::MaterialDomain)?;
        let change = self
            .owner
            .model
            .entropy_change_j_per_k(&initial, &final_state)
            .map_err(|_| DepositionError::MaterialDomain)?
            / self.mass_kg;
        // The accepted numerical inverse must resolve the requested interval.
        // Never turn an unresolved positive heat interval into zero entropy.
        if !change.is_finite()
            || (final_h_j_kg > initial_h_j_kg && change <= 0.0)
            || (final_h_j_kg < initial_h_j_kg && change >= 0.0)
        {
            return Err(DepositionError::MaterialDomain);
        }
        Ok(Some(change))
    }
}

fn positive(value: f64) -> bool {
    value.is_finite() && value > 0.0
}
