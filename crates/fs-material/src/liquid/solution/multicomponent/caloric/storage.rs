//! Fixed-composition homogeneous liquid storage on the native thermal owner.
//!
//! Every inverse and entropy change delegates to the original multicomponent
//! caloric model. Its pure Cp and Cp/T sources retain every original knot;
//! the original A+B/T mixing enthalpy is constant along this homogeneous path.
//! Closed storage measures h from the cold endpoint, avoiding large datum
//! subtraction. Explicit stream storage instead retains the common pure-species
//! datum, including nonideal mixing heat and entropy. Cell mass need not equal
//! the inventory used for normalization.
//!
//! Pressure, composition and reference geometry stay fixed. This deliberately
//! includes constrained metastable homogeneous liquids; no equilibrium split,
//! phase conductivity, density evolution or momentum closure is inferred.
//! Numeric thermal trials do not refresh other mixture transport properties.
//! The conduction adapter separately admits the complete original k(T) path.

/// Original mixture conductivity and analytic implicit thermal tangents.
pub mod conduction;
/// Common-datum homogeneous liquids on stationary species/enthalpy flow graphs.
pub mod stream;

use super::{
    IsobaricHeatCapacityLaw, MulticomponentCaloricSolveV1, MulticomponentCaloricStateV1,
    MulticomponentLiquidCaloricError, MulticomponentLiquidCaloricLaw, Result, SourceOwner, poll,
};
use fs_matdb::QueryPoint;
use fs_qty::Temperature;
use fs_thermal::enthalpy::{DepositionError, ThermalLaw};
use std::sync::Arc;

/// Immutable original homogeneous calorics at one conserved composition.
#[derive(Clone, Debug, PartialEq)]
pub struct MulticomponentLiquidThermalStorage {
    owner: Arc<SourceOwner>,
    endpoints: [MulticomponentCaloricStateV1; 2],
    source_endpoints: [QueryPoint; 2],
    mass_kg: f64,
    lower_specific_h_j_kg: f64,
    upper_specific_h_j_kg: f64,
    maximum_slope: f64,
    controls: MulticomponentCaloricSolveV1,
    stream: Option<stream::StreamBinding>,
}

impl MulticomponentLiquidCaloricLaw<'_, '_> {
    /// Admit the original material bundle at both thermal endpoints and retain
    /// one mass-normalized homogeneous caloric law. The original heat solver
    /// admits all root controls; unsupported material endpoints refuse atomically.
    ///
    /// The reusable CFL bound sums downward-rounded positive source minima.
    /// Each minimum covers the original source segments intersecting the
    /// prepared temperature bracket, including clipped endpoint segments.
    /// It is conservative; the physical local Cp always uses the exact source
    /// interpolation, independently of this bound.
    pub fn homogeneous_thermal_storage(
        &self,
        component_moles: &[f64],
        controls: MulticomponentCaloricSolveV1,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<MulticomponentLiquidThermalStorage> {
        poll(&mut cancelled)?;
        let cold = self.homogeneous(
            controls.temperature_bracket[0],
            component_moles,
            &mut cancelled,
        )?;
        self.model()
            .add_heat(cold.caloric_state(), 0.0, controls, &mut cancelled)?;
        let hot = self.homogeneous(
            controls.temperature_bracket[1],
            component_moles,
            &mut cancelled,
        )?;
        let mass_kg = cold.liquid().mass_kg();
        let span = self
            .model()
            .enthalpy_change_j(cold.caloric_state(), hot.caloric_state())?
            / mass_kg;
        let mut minimum_capacity = 0.0;
        let [lower, upper] = controls.temperature_bracket.map(|t| t.value());
        for ((&amount, source), pure) in component_moles
            .iter()
            .zip(self.sources())
            .zip(self.model().heat_capacities())
        {
            poll(&mut cancelled)?;
            if amount == 0.0 {
                continue;
            }
            let mut minimum = f64::INFINITY;
            for segment in source.source_line().knots.windows(2) {
                if segment[1].0 < lower || segment[0].0 > upper {
                    continue;
                }
                for (temperature, _) in segment {
                    // Evaluate the original molar owner at its exact knots,
                    // retaining its mass-to-molar rounding and source units.
                    let cp = pure
                        .heat_capacity_at(Temperature::new(*temperature))
                        .map_err(super::MulticomponentCaloricErrorV1::from)?
                        .value();
                    minimum = minimum.min(cp);
                }
            }
            let contribution = (amount * minimum).next_down().max(0.0);
            minimum_capacity = (minimum_capacity + contribution).next_down();
        }
        let maximum_slope = (mass_kg / minimum_capacity).next_up();
        if !positive(span) || !positive(minimum_capacity) || !positive(maximum_slope) {
            return Err(MulticomponentLiquidCaloricError::Invalid(
                "unrepresentable homogeneous storage span or heat-capacity bound",
            ));
        }
        poll(&mut cancelled)?;
        Ok(MulticomponentLiquidThermalStorage {
            owner: Arc::clone(&self.owner),
            endpoints: [cold.caloric_state().clone(), hot.caloric_state().clone()],
            source_endpoints: [
                cold.liquid().liquid().material().query(),
                hot.liquid().liquid().material().query(),
            ],
            mass_kg,
            lower_specific_h_j_kg: 0.0,
            upper_specific_h_j_kg: span,
            maximum_slope,
            controls,
            stream: None,
        })
    }
}

impl MulticomponentLiquidThermalStorage {
    pub const fn temperature_domain(&self) -> [Temperature; 2] {
        self.controls.temperature_bracket
    }
    /// Inclusive h bounds [J/kg]. Closed storage uses zero at the cold endpoint;
    /// explicit stream storage retains the original common pure-species datum.
    pub const fn specific_enthalpy_domain_j_kg(&self) -> [f64; 2] {
        [self.lower_specific_h_j_kg, self.upper_specific_h_j_kg]
    }
    pub fn sources(&self) -> &[IsobaricHeatCapacityLaw] {
        &self.owner.sources
    }
    pub fn component_moles(&self) -> &[f64] {
        self.endpoints[0].component_moles()
    }

    /// Exact homogeneous Cp from the original pure-liquid sources [J/(kg K)].
    pub fn specific_heat_cp_j_kg_k(
        &self,
        temperature: Temperature,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<f64> {
        poll(&mut cancelled)?;
        self.admit_temperature(temperature)?;
        let state =
            self.owner
                .model
                .homogeneous(temperature, self.component_moles(), &mut cancelled)?;
        let cp = state.isobaric_heat_capacity_j_per_k() / self.mass_kg;
        if !positive(cp) {
            return Err(MulticomponentLiquidCaloricError::Invalid(
                "unrepresentable homogeneous specific heat",
            ));
        }
        Ok(cp)
    }

    /// Original sensible enthalpy increment, without datum-dependent subtraction.
    pub fn specific_enthalpy_j_kg(
        &self,
        temperature: Temperature,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<f64> {
        poll(&mut cancelled)?;
        self.admit_temperature(temperature)?;
        if temperature == self.temperature_domain()[0] {
            return Ok(self.lower_specific_h_j_kg);
        }
        if temperature == self.temperature_domain()[1] {
            return Ok(self.upper_specific_h_j_kg);
        }
        let state =
            self.owner
                .model
                .homogeneous(temperature, self.component_moles(), &mut cancelled)?;
        let increment = self
            .owner
            .model
            .enthalpy_change_j(&self.endpoints[0], &state)?
            / self.mass_kg;
        let h = self.lower_specific_h_j_kg + increment;
        if !positive(increment)
            || h <= self.lower_specific_h_j_kg
            || h >= self.upper_specific_h_j_kg
        {
            return Err(MulticomponentLiquidCaloricError::Invalid(
                "interior homogeneous storage enthalpy is unrepresentable",
            ));
        }
        poll(&mut cancelled)?;
        Ok(h)
    }

    /// Reuse the original safeguarded caloric inverse and its bounded work.
    pub fn caloric_state_at_specific_enthalpy(
        &self,
        h: f64,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<MulticomponentCaloricStateV1> {
        poll(&mut cancelled)?;
        self.admit(h)?;
        if h == self.lower_specific_h_j_kg {
            return Ok(self.endpoints[0].clone());
        }
        if h == self.upper_specific_h_j_kg {
            return Ok(self.endpoints[1].clone());
        }
        let heat_j = (h - self.lower_specific_h_j_kg) * self.mass_kg;
        if !positive(heat_j) {
            return Err(MulticomponentLiquidCaloricError::Invalid(
                "positive homogeneous storage heat is unrepresentable",
            ));
        }
        let balance =
            self.owner
                .model
                .add_heat(&self.endpoints[0], heat_j, self.controls, &mut cancelled)?;
        Ok(balance.final_state().clone())
    }

    fn admit_temperature(&self, temperature: Temperature) -> Result<()> {
        let [lower, upper] = self.temperature_domain();
        if !temperature.value().is_finite()
            || temperature.value() < lower.value()
            || temperature.value() > upper.value()
        {
            Err(MulticomponentLiquidCaloricError::Invalid(
                "temperature outside homogeneous storage domain",
            ))
        } else {
            Ok(())
        }
    }
    fn admit(&self, h: f64) -> Result<()> {
        if !h.is_finite() || h < self.lower_specific_h_j_kg || h > self.upper_specific_h_j_kg {
            Err(MulticomponentLiquidCaloricError::Invalid(
                "enthalpy outside homogeneous storage domain",
            ))
        } else {
            Ok(())
        }
    }
}

impl ThermalLaw for MulticomponentLiquidThermalStorage {
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
    fn advected_entropy_anchor(
        &self,
    ) -> Option<fs_thermal::enthalpy::advection::species::CaloricEntropyAnchor<'_>> {
        self.stream
            .as_ref()
            .map(stream::StreamBinding::entropy_anchor)
    }
    fn temperature_k(&self, h: f64) -> core::result::Result<f64, DepositionError> {
        self.caloric_state_at_specific_enthalpy(h, || false)
            .map(|state| state.temperature().value())
            .map_err(|_| DepositionError::MaterialDomain)
    }
    fn maximum_temperature_slope(&self) -> f64 {
        self.maximum_slope
    }
    fn specific_entropy_change_j_kg_k(
        &self,
        initial_h: f64,
        final_h: f64,
    ) -> core::result::Result<Option<f64>, DepositionError> {
        self.admit(initial_h)
            .and_then(|()| self.admit(final_h))
            .map_err(|_| DepositionError::MaterialDomain)?;
        if initial_h == final_h {
            return Ok(Some(0.0));
        }
        let initial = self
            .caloric_state_at_specific_enthalpy(initial_h, || false)
            .map_err(|_| DepositionError::MaterialDomain)?;
        let final_state = self
            .caloric_state_at_specific_enthalpy(final_h, || false)
            .map_err(|_| DepositionError::MaterialDomain)?;
        let change = self
            .owner
            .model
            .entropy_change_j_per_k(&initial, &final_state)
            .map_err(|_| DepositionError::MaterialDomain)?
            / self.mass_kg;
        if !change.is_finite()
            || (final_h > initial_h && change <= 0.0)
            || (final_h < initial_h && change >= 0.0)
        {
            return Err(DepositionError::MaterialDomain);
        }
        Ok(Some(change))
    }
}

fn positive(value: f64) -> bool {
    value.is_finite() && value > 0.0
}
