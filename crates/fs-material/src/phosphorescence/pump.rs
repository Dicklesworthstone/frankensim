//! Original molecular pump cross sections at an exact material condition.
//!
//! Each value is square metres per ground-state emitter, ordered `[G -> S,
//! G -> T]`. Both channels may absorb the same wavelength; transport must use
//! their sum for extinction and retain their separate excitation destinations.
//! These are neither bulk attenuation nor mass absorption coefficients. No
//! number density, oscillator strength or molecular data is inferred.
//!
//! Wavelength is vacuum SI metres. Native source queries own interpolation,
//! selected claims and validity. An explicitly owned absolute-temperature axis
//! may refresh with the thermal cell; all other external coordinates remain
//! fixed in the convenience refresh. This snapshot evolves no populations.

use super::{
    ContentHash, Dims, DomainHasher, MaterialCard, MaterialPropertySelection, PhosphorescenceError,
    QuantitySpec, QueryPoint, ResolvedMaterialStatePoint, ScalarAdmissibility,
    ScalarPropertyRequirement, resolve_material_state_point,
};
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};

/// Original molecular cross sections, ordered `[ground -> singlet, ground -> triplet]`.
pub const GROUND_STATE_PUMP_CROSS_SECTION_PROPERTIES: [&str; 2] = [
    "luminescence_ground_singlet_absorption_cross_section",
    "luminescence_ground_triplet_absorption_cross_section",
];
/// SI square metres per emitter, not square metres per kilogram or per mole.
pub const GROUND_STATE_PUMP_CROSS_SECTION_DIMS: Dims = Dims([2, 0, 0, 0, 0, 0]);
const WAVELENGTH_QUANTITY: QuantitySpec = QuantitySpec::dimensional(Dims([1, 0, 0, 0, 0, 0]));
const TEMPERATURE_QUANTITY: QuantitySpec = QuantitySpec::semantic(SemanticType::new(
    QuantityKind::AbsoluteTemperature,
    ValueForm::Static,
));

/// Immutable original source bundle and the two cross sections at its point.
/// Identities attribute the source and state; they do not certify measurements.
#[derive(Clone, Debug, PartialEq)]
pub struct GroundStatePumpStatePoint {
    resolved: ResolvedMaterialStatePoint,
    wavelength_axis: String,
    temperature_axis: Option<String>,
    wavelength_m: f64,
    temperature_k: Option<f64>,
    cross_sections_m2: [f64; 2],
    identity: ContentHash,
}

impl GroundStatePumpStatePoint {
    /// Resolve both explicit transition roles. Zero disables a channel;
    /// missing, ambiguous or negative coefficients never receive defaults.
    pub fn resolve(
        card: &MaterialCard,
        point: &QueryPoint,
        selection: MaterialPropertySelection,
        wavelength_axis: &str,
        temperature_axis: Option<&str>,
    ) -> Result<Self, PhosphorescenceError> {
        let requirements = GROUND_STATE_PUMP_CROSS_SECTION_PROPERTIES
            .map(|name| {
                ScalarPropertyRequirement::try_new(
                    name,
                    GROUND_STATE_PUMP_CROSS_SECTION_DIMS,
                    ScalarAdmissibility::NonNegative,
                )
            })
            .into_iter()
            .collect::<Result<Vec<_>, _>>()?;
        Self::try_from_material(
            resolve_material_state_point(card, point, &requirements, selection)?,
            wavelength_axis,
            temperature_axis,
        )
    }

    /// Retain a complete original bundle while rechecking SI quantities,
    /// nonnegative values and the explicitly owned wavelength/temperature axes.
    /// Frequency, field, stress, strain and history response or bounds refuse.
    pub fn try_from_material(
        resolved: ResolvedMaterialStatePoint,
        wavelength_axis: &str,
        temperature_axis: Option<&str>,
    ) -> Result<Self, PhosphorescenceError> {
        if temperature_axis == Some(wavelength_axis) {
            return Err(invalid(
                "wavelength and temperature require distinct coordinates",
            ));
        }
        let wavelength_m = coordinate(&resolved, wavelength_axis, WAVELENGTH_QUANTITY)?;
        let temperature_k = temperature_axis
            .map(|axis| coordinate(&resolved, axis, TEMPERATURE_QUANTITY))
            .transpose()?;
        let owned = [
            (wavelength_axis, WAVELENGTH_QUANTITY),
            (temperature_axis.unwrap_or(""), TEMPERATURE_QUANTITY),
        ];
        let owned = &owned[..1 + usize::from(temperature_axis.is_some())];
        let mut cross_sections_m2 = [0.0; 2];
        for (target, name) in GROUND_STATE_PUMP_CROSS_SECTION_PROPERTIES
            .into_iter()
            .enumerate()
        {
            let property = resolved
                .property(name)
                .ok_or(PhosphorescenceError::Property {
                    name,
                    reason: "missing or ambiguous original molecular cross section",
                })?;
            let requirement = property.requirement();
            let value = property.value_si();
            if requirement.quantity()
                != QuantitySpec::dimensional(GROUND_STATE_PUMP_CROSS_SECTION_DIMS)
                || requirement.hardness_test().is_some()
                || requirement.elastic_component().is_some()
                || requirement.strain_component().is_some()
                || requirement.stress_component().is_some()
                || !value.is_finite()
                || value < 0.0
                || !property.has_only_owned_coordinates_and_fixed_external_conditions(owned)
            {
                return Err(PhosphorescenceError::Property {
                    name,
                    reason: "requires a nonnegative molecular area at owned wavelength/temperature and fixed external conditions",
                });
            }
            cross_sections_m2[target] = value;
        }
        let mut hash = DomainHasher::new("org.frankensim.material.ground-state-pump.v1");
        hash.update(resolved.identity().as_bytes());
        hash.update(&(wavelength_axis.len() as u64).to_le_bytes());
        hash.update(wavelength_axis.as_bytes());
        hash.update(&[u8::from(temperature_axis.is_some())]);
        if let Some(axis) = temperature_axis {
            hash.update(&(axis.len() as u64).to_le_bytes());
            hash.update(axis.as_bytes());
        }
        Ok(Self {
            resolved,
            wavelength_axis: wavelength_axis.to_owned(),
            temperature_axis: temperature_axis.map(str::to_owned),
            wavelength_m,
            temperature_k,
            cross_sections_m2,
            identity: hash.finalize(),
        })
    }

    /// Exact original card, complete query, requirements and selected receipts.
    pub const fn resolved(&self) -> &ResolvedMaterialStatePoint {
        &self.resolved
    }
    /// Identity of the complete source snapshot and owned coordinate roles.
    pub const fn identity(&self) -> ContentHash {
        self.identity
    }
    /// Native `[G -> S, G -> T]` cross sections, square metres per ground emitter.
    pub const fn cross_sections_m2(&self) -> [f64; 2] {
        self.cross_sections_m2
    }
    /// Original vacuum wavelength in SI metres.
    pub const fn wavelength_m(&self) -> f64 {
        self.wavelength_m
    }
    /// Explicitly owned absolute temperature, when declared, in kelvin.
    pub const fn temperature_k(&self) -> Option<f64> {
        self.temperature_k
    }
    /// Exact owned wavelength coordinate name; no wavelength unit is inferred.
    pub fn wavelength_axis(&self) -> &str {
        &self.wavelength_axis
    }
    /// Exact owned absolute-temperature coordinate, when declared.
    pub fn temperature_axis(&self) -> Option<&str> {
        self.temperature_axis.as_deref()
    }

    /// Refresh the complete point using the same card, quantity conventions and
    /// exact selected claims. Unsupported wavelengths or conditions refuse.
    pub fn at_query_point(
        &self,
        card: &MaterialCard,
        point: &QueryPoint,
    ) -> Result<Self, PhosphorescenceError> {
        Self::try_from_material(
            self.resolved.at_query_point(card, point)?,
            &self.wavelength_axis,
            self.temperature_axis.as_deref(),
        )
    }

    /// Evaluate one actual pump node and temperature without changing any
    /// other external condition. `Some(T)` is required exactly when a thermal
    /// axis was bound; a supplied temperature is never silently ignored.
    pub fn at_wavelength_temperature(
        &self,
        card: &MaterialCard,
        wavelength_m: f64,
        temperature_k: Option<f64>,
    ) -> Result<Self, PhosphorescenceError> {
        if !wavelength_m.is_finite()
            || wavelength_m <= 0.0
            || temperature_k.is_some() != self.temperature_axis.is_some()
            || temperature_k.is_some_and(|value| !value.is_finite() || value <= 0.0)
        {
            return Err(invalid(
                "positive wavelength and matching owned temperature required",
            ));
        }
        let mut point = self
            .resolved
            .query()
            .with(self.wavelength_axis.as_str(), wavelength_m)
            .map_err(|_| invalid("invalid pump wavelength coordinate"))?;
        if let (Some(axis), Some(value)) = (&self.temperature_axis, temperature_k) {
            point = point
                .with(axis.as_str(), value)
                .map_err(|_| invalid("invalid pump temperature coordinate"))?;
        }
        self.at_query_point(card, &point)
    }
}

fn coordinate(
    resolved: &ResolvedMaterialStatePoint,
    axis: &str,
    quantity: QuantitySpec,
) -> Result<f64, PhosphorescenceError> {
    let value = resolved.state_coordinate(axis);
    if resolved.axis_quantities().get(axis) != Some(&quantity)
        || !value.is_some_and(|value| value.is_finite() && value > 0.0)
    {
        return Err(invalid(
            "positive explicitly typed SI wavelength or absolute temperature required",
        ));
    }
    Ok(value.expect("validated pump coordinate"))
}

fn invalid(reason: &'static str) -> PhosphorescenceError {
    PhosphorescenceError::Property {
        name: "ground-state molecular pump",
        reason,
    }
}
