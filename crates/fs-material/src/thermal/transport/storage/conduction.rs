//! Native implicit conduction with additive, amount-bound sensible storage.
//!
//! The original assembly stores e(T) = integral C_total(T) dT / m_ref.
//! Therefore dT/de = m_ref/C_total(T) and dk/de = (dk/dT)*dT/de. The existing
//! thermal graph owns heat fluxes, energy residuals, time stepping and acceptance.
//! No constituent capacity is frozen or converted to a thermodynamic Cp or Cv.

use super::super::{
    ThermalTransportError, admit_temperature_point, conductivity::PreparedScalarConductivity,
    positive,
};
use super::HeatStorageThermalLaw;
use fs_blake3::{ContentHash, DomainHasher};
use fs_matdb::{ClaimSelection, MaterialCard, PropertyKey, PropertyLineSlice, QueryPoint};
use fs_thermal::enthalpy::conduction::ConductivityLaw;
use fs_thermal::enthalpy::conduction::implicit::{
    ImplicitConductionLaw, ImplicitConductionObservables, ImplicitConductionState,
};
use fs_thermal::enthalpy::{DepositionError, ThermalLaw};

/// One exact storage assembly and an explicitly selected effective conductivity.
///
/// The caller owns the conductivity model's suitability for this assembly,
/// direction and fixed geometry. It is never averaged from constituent heat
/// capacities. The conductivity card may therefore differ from the constituent
/// cards; all original identities and fixed source conditions remain retained.
/// The co-temperature approximation, invariant total mass, optional sensible
/// entropy and `DeclaredSensibleStorage` basis remain with the borrowed law.
#[derive(Clone, Debug)]
pub struct HeatStorageConductionLaw<'law, 'storage> {
    caloric: &'law HeatStorageThermalLaw<'storage>,
    conductivity: PreparedScalarConductivity,
    conductivity_card: ContentHash,
    identity: ContentHash,
}

impl<'law, 'storage> HeatStorageConductionLaw<'law, 'storage> {
    /// Prepare source conductivity over the entire original storage interval.
    ///
    /// Both endpoints use positive absolute kelvin and must exactly match the
    /// storage span. Remaining typed source coordinates stay fixed. Canonical
    /// source-line admission retains every k(T) breakpoint and refuses holes,
    /// discrete-only support, changing claims and nonpositive conductivity.
    pub fn try_new(
        card: &MaterialCard,
        caloric: &'law HeatStorageThermalLaw<'storage>,
        conductivity: &PropertyKey,
        temperature_axis: &str,
        lower: &QueryPoint,
        upper: &QueryPoint,
        selection: ClaimSelection,
    ) -> Result<Self, ThermalTransportError> {
        for point in [lower, upper] {
            admit_temperature_point(point, temperature_axis)?;
        }
        let span = caloric.storage().temperature_span_k();
        if lower.axes()[temperature_axis] != span[0] || upper.axes()[temperature_axis] != span[1] {
            return Err(ThermalTransportError::Invalid(
                "assembly conductivity must cover the complete original storage interval",
            ));
        }
        let conductivity = PreparedScalarConductivity::try_new(
            card,
            conductivity,
            temperature_axis,
            lower,
            upper,
            selection,
        )?;
        let conductivity_card = card.content_hash();
        let mut hash = DomainHasher::new("org.frankensim.fs-material.heat-storage-conduction.v1");
        hash.update(caloric.identity().as_bytes());
        hash.update(conductivity_card.as_bytes());
        conductivity.append_identity(&mut hash);
        Ok(Self {
            caloric,
            conductivity,
            conductivity_card,
            identity: hash.finalize(),
        })
    }

    /// The exact fixed-mass caloric object registered on the thermal timeline.
    #[must_use]
    pub const fn caloric(&self) -> &'law HeatStorageThermalLaw<'storage> {
        self.caloric
    }

    /// Original conductivity card, independently bound to this assembly model.
    #[must_use]
    pub const fn conductivity_card_identity(&self) -> ContentHash {
        self.conductivity_card
    }

    /// Original conductivity support and source receipts, without resampling.
    #[must_use]
    pub const fn source_line(&self) -> &PropertyLineSlice {
        self.conductivity.source_line()
    }

    pub(in crate::thermal::transport) fn integrated_conductivity_curve(
        &self,
    ) -> Result<crate::thermal::conduction::ConductivityCurve, DepositionError> {
        self.conductivity.integrated_curve()
    }

    #[must_use]
    pub const fn identity(&self) -> ContentHash {
        self.identity
    }

    /// One original caloric inverse and same-state analytic tangents.
    ///
    /// Successful evaluation allocates nothing and performs no source query.
    /// Interior k knots use the right derivative; the upper endpoint uses the
    /// left derivative. The assembled total capacity remains continuous.
    pub fn state_at_specific_storage(
        &self,
        storage_j_kg: f64,
    ) -> Result<ImplicitConductionState, ThermalTransportError> {
        let invalid = ThermalTransportError::Invalid;
        let temperature_k = self
            .caloric
            .temperature_k(storage_j_kg)
            .map_err(|_| invalid("storage outside original heat-capacity assembly"))?;
        let capacity = self.caloric.storage().heat_capacity_j_k(temperature_k)?;
        let (conductivity, dk_dt) = self.conductivity.at_temperature(temperature_k)?;
        let temperature_slope = self.caloric.reference_mass_kg() / capacity;
        let conductivity_slope = dk_dt * temperature_slope;
        if !positive(temperature_slope)
            || !conductivity_slope.is_finite()
            || (dk_dt != 0.0 && conductivity_slope == 0.0)
        {
            return Err(invalid(
                "assembly conduction tangents are not representable",
            ));
        }
        Ok(ImplicitConductionState {
            temperature_k,
            temperature_slope_kg_k_j: temperature_slope,
            conductivity_w_m_k: conductivity,
            conductivity_slope_kg_w_m_k_j: conductivity_slope,
        })
    }
}

impl ConductivityLaw for HeatStorageConductionLaw<'_, '_> {
    fn conductivity_w_m_k(&self, storage_j_kg: f64) -> Result<f64, DepositionError> {
        let temperature = self.caloric.temperature_k(storage_j_kg)?;
        self.conductivity
            .at_temperature(temperature)
            .map(|(conductivity, _)| conductivity)
            .map_err(|_| DepositionError::MaterialDomain)
    }
}

impl ImplicitConductionLaw for HeatStorageConductionLaw<'_, '_> {
    fn caloric_law(&self) -> &dyn ThermalLaw {
        self.caloric
    }

    fn conductivity_law(&self) -> &dyn ConductivityLaw {
        self
    }

    fn enthalpy_bounds_j_kg(&self) -> [f64; 2] {
        self.caloric.specific_storage_span_j_kg()
    }

    fn state_at_specific_enthalpy(
        &self,
        storage_j_kg: f64,
    ) -> Result<ImplicitConductionState, DepositionError> {
        self.state_at_specific_storage(storage_j_kg)
            .map_err(|_| DepositionError::MaterialDomain)
    }

    fn time_observables(
        &self,
        storage_j_kg: f64,
    ) -> Result<ImplicitConductionObservables, DepositionError> {
        Ok(ImplicitConductionObservables {
            temperature_k: self.caloric.temperature_k(storage_j_kg)?,
            phase_mass_fraction: None,
        })
    }
}
