//! Source-resolved bulk and uniform DC conductor primitives.
//!
//! Bulk resistivity belongs to the material state, not to a chosen geometry.
//! Spatial solvers consume [`ResolvedConductorState`]; lumped conductors lower
//! that same state with explicit length and cross-section. Neither path chooses
//! a metal, adds contact resistance, or infers a temperature-evolution law.

/// Shared planar scalar/anisotropic conductivity with fixed material frames.
pub mod planar;
/// Shared 3D scalar/orthotropic conductivity with explicit principal directions.
pub mod spatial;
/// Pinned scalar transport at explicitly owned electric/current magnitudes.
pub mod response;
/// Source-bound electrothermal law nodes for the shared constitutive graph.
pub mod graph;

use core::fmt;

use fs_matdb::{MaterialCard, PropertyKey, QueryPoint};
use fs_qty::{Dims, QuantitySpec};

use crate::state_point::{
    MaterialPropertySelection, MaterialStatePointError, ResolvedMaterialStatePoint,
    ResolvedScalarProperty, ScalarAdmissibility, ScalarPropertyRequirement,
    resolve_material_state_point,
};

/// Canonical source property for bulk electrical resistivity.
pub const ELECTRICAL_RESISTIVITY_PROPERTY: &str = "electrical_resistivity";

/// Electrical resistivity dimensions `[L, M, t, T, I, N]` for ohm metre.
pub const ELECTRICAL_RESISTIVITY_DIMS: Dims = Dims([3, 1, -3, 0, -2, 0]);

/// Canonical source property for isotropic DC electrical conductivity.
pub const ELECTRICAL_CONDUCTIVITY_PROPERTY: &str = "electrical_conductivity";

/// Electrical conductivity dimensions for siemens per metre.
pub const ELECTRICAL_CONDUCTIVITY_DIMS: Dims = Dims([-3, -1, 3, 0, 2, 0]);

/// Independent scalar source quantity for the same isotropic DC Ohmic law.
///
/// Conductivity and resistivity are reciprocals only under this stated scalar
/// law. This does not invert anisotropic tensors, AC admittances, or infer a
/// relation between independently measured curves at different conditions.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ElectricalPropertyBasis {
    /// A measured/declared resistivity is used directly.
    Resistivity,
    /// Resistivity is derived as `1 / electrical_conductivity` at the query point.
    Conductivity,
}

impl ElectricalPropertyBasis {
    /// Exact canonical property name supplied by the retained source bundle.
    #[must_use]
    pub const fn property_name(self) -> &'static str {
        match self {
            Self::Resistivity => ELECTRICAL_RESISTIVITY_PROPERTY,
            Self::Conductivity => ELECTRICAL_CONDUCTIVITY_PROPERTY,
        }
    }

    /// Dimensional source requirement; no quantity semantics are inferred.
    #[must_use]
    pub const fn dimensions(self) -> Dims {
        match self {
            Self::Resistivity => ELECTRICAL_RESISTIVITY_DIMS,
            Self::Conductivity => ELECTRICAL_CONDUCTIVITY_DIMS,
        }
    }

    /// Requirement for inclusion in an existing multiphysics material bundle.
    pub fn requirement(self) -> Result<ScalarPropertyRequirement, MaterialStatePointError> {
        ScalarPropertyRequirement::try_new(
            self.property_name(), self.dimensions(), ScalarAdmissibility::StrictlyPositive,
        )
    }

    /// Bind the independent electrical quantity to an exact source key. Names
    /// are never normalized or guessed; the scalar DC dimensions remain fixed.
    pub fn requirement_for(
        self,
        key: &PropertyKey,
    ) -> Result<ScalarPropertyRequirement, ConductorError> {
        if key.quantity() != QuantitySpec::dimensional(self.dimensions()) {
            return Err(ConductorError::InvalidInput {
                quantity: "exact electrical source key must match the selected scalar DC basis",
            });
        }
        Ok(ScalarPropertyRequirement::try_with_key(
            key,
            ScalarAdmissibility::StrictlyPositive,
        )?)
    }
}

/// Ideal quasi-static electrical boundary condition for a uniform resistor.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum OhmicDrive {
    /// Prescribed signed current [A].
    Current(f64),
    /// Prescribed signed potential difference [V], with `I = V/R`.
    Voltage(f64),
}

/// Geometry-independent isotropic DC conductor state with source-use evidence.
///
/// Resolve again at each changed physical state. This object is immutable and
/// does not extrapolate a temperature coefficient from one sampled resistivity.
#[derive(Clone, Debug, PartialEq)]
pub struct ResolvedConductorState {
    material: ResolvedMaterialStatePoint,
    source_basis: ElectricalPropertyBasis,
    source_key: PropertyKey,
    resistivity_ohm_m: f64,
    liquid_mass_fraction_axis: Option<String>,
}

impl ResolvedConductorState {
    /// Adopt an already resolved multiphysics bundle without querying a second
    /// material snapshot or discarding its mechanical/thermal source receipts.
    /// The complete original bundle and query conditions are retained intact.
    /// This does not refresh or expand the source's domain of validity.
    ///
    /// # Errors
    /// The requested bundle must contain exactly one electrical source basis:
    /// resistivity or conductivity. Thus changing that source requirement is a
    /// data change for consumers of this method, including layered conductors.
    /// A bundle containing both bases requires an explicit independent-variable
    /// choice through [`Self::try_from_material_with_basis`]; conflicting data is
    /// never averaged, repaired or silently used as a fallback.
    ///
    /// Missing/ambiguous sources, incompatible quantities, nonpositive values
    /// and an unrepresentable reciprocal refuse. Extra properties are retained.
    pub fn try_from_material(
        material: ResolvedMaterialStatePoint,
    ) -> Result<Self, ConductorError> {
        let has = |name| material.properties().iter().any(|p| p.requirement().name() == name);
        let basis = match (has(ELECTRICAL_RESISTIVITY_PROPERTY), has(ELECTRICAL_CONDUCTIVITY_PROPERTY)) {
            (true, false) => ElectricalPropertyBasis::Resistivity,
            (false, true) => ElectricalPropertyBasis::Conductivity,
            (true, true) => return Err(ConductorError::InvalidInput {
                quantity: "both electrical source bases are present; select the independent quantity explicitly",
            }),
            (false, false) => return Err(ConductorError::InvalidInput {
                quantity: "resolved electrical_resistivity or electrical_conductivity",
            }),
        };
        Self::try_from_material_with_basis(material, basis)
    }

    /// Adopt one explicitly selected independent quantity from a complete
    /// source bundle. Other observations remain available unchanged; they are
    /// not treated as independent inputs to, or silently reconciled with, this
    /// scalar Ohmic law. Source uncertainty remains on the original property;
    /// taking a reciprocal does not turn it into an independent measurement.
    ///
    /// # Errors
    /// Missing/ambiguous selected properties, wrong quantities, invalid values
    /// or an unrepresentable conductivity-to-resistivity conversion.
    pub fn try_from_material_with_basis(
        material: ResolvedMaterialStatePoint,
        source_basis: ElectricalPropertyBasis,
    ) -> Result<Self, ConductorError> {
        let key = PropertyKey::new(source_basis.property_name(), source_basis.dimensions());
        Self::try_from_material_with_property(material, source_basis, &key)
    }

    /// Adopt an exact electrical property from the original complete bundle.
    /// The caller declares whether that property is scalar DC resistivity or
    /// conductivity. Its key, conditions, support and receipts survive all
    /// lowering and refresh; an absent key never falls back to a similar name.
    ///
    /// # Errors
    /// A key with the wrong quantity, a missing exact property, a nonpositive
    /// source value, or an unrepresentable conductivity reciprocal.
    pub fn try_from_material_with_property(
        material: ResolvedMaterialStatePoint,
        source_basis: ElectricalPropertyBasis,
        source_key: &PropertyKey,
    ) -> Result<Self, ConductorError> {
        source_basis.requirement_for(source_key)?;
        let property = material.property_by_key(source_key)
            .ok_or(ConductorError::InvalidInput {
                quantity: source_basis.property_name(),
            })?;
        let value = property.value_si();
        if property.requirement().quantity()
            != QuantitySpec::dimensional(source_basis.dimensions())
            || !value.is_finite() || value <= 0.0
        {
            return Err(ConductorError::InvalidInput {
                quantity: source_basis.property_name(),
            });
        }
        let resistivity_ohm_m = match source_basis {
            ElectricalPropertyBasis::Resistivity => value,
            ElectricalPropertyBasis::Conductivity => {
                let reciprocal = 1.0 / value;
                if !reciprocal.is_finite() || reciprocal <= 0.0 {
                    return Err(ConductorError::InvalidDerived {
                        quantity: "resistivity from electrical conductivity",
                    });
                }
                reciprocal
            }
        };
        Ok(Self {
            material,
            source_basis,
            source_key: source_key.clone(),
            resistivity_ohm_m,
            liquid_mass_fraction_axis: None,
        })
    }

    /// Declare one phase coordinate supplied by a physical enthalpy owner.
    /// The already-resolved coordinate must be a dimensional scalar in [0,1].
    /// This admits a stationary source snapshot; it neither evolves phase nor
    /// invents a solid/liquid mixing rule. Every other source axis must retain
    /// its existing fixed-external meaning. Coupled owners must bind this exact
    /// axis to their actual liquid MASS fraction on every refresh.
    ///
    /// # Errors
    /// Missing, mistyped or out-of-range phase coordinates, a changed declared
    /// axis, or another source-dependent coordinate without a static owner.
    pub fn with_liquid_mass_fraction_axis(
        self,
        axis: impl Into<String>,
    ) -> Result<Self, ConductorError> {
        self.with_phase_axis_and_coordinates(axis.into(), &[])
    }

    // Only a physical response owner may admit additional changing coordinates.
    // The ordinary phase snapshot and public refresh remain strictly static.
    fn with_phase_axis_and_coordinates(
        mut self,
        axis: String,
        additional: &[(&str, QuantitySpec)],
    ) -> Result<Self, ConductorError> {
        validate_liquid_mass_fraction_axis(&self.material, &axis)?;
        let mut owned = vec![(axis.as_str(), QuantitySpec::dimensional(Dims([0; 6])))];
        owned.extend_from_slice(additional);
        if self
            .liquid_mass_fraction_axis
            .as_ref()
            .is_some_and(|old| old != &axis)
            || !self.source_property()
                .has_only_owned_coordinates_and_fixed_external_conditions(&owned)
        {
            return Err(ConductorError::InvalidInput {
                quantity: "phase conductor requires one retained phase axis and admitted owned source conditions",
            });
        }
        self.liquid_mass_fraction_axis = Some(axis);
        Ok(self)
    }

    /// Explicit stationary liquid-mass-fraction coordinate, when selected.
    pub fn liquid_mass_fraction_axis(&self) -> Option<&str> {
        self.liquid_mass_fraction_axis.as_deref()
    }

    /// Exact source-resolved material bundle, including its usage receipt.
    #[must_use]
    pub const fn material(&self) -> &ResolvedMaterialStatePoint {
        &self.material
    }

    /// Independent source quantity and its explicit algebraic interpretation.
    #[must_use]
    pub const fn source_basis(&self) -> ElectricalPropertyBasis { self.source_basis }

    /// Complete original source key, including its quantity declaration.
    #[must_use]
    pub const fn source_key(&self) -> &PropertyKey { &self.source_key }

    /// Rebuild this same electrical law around a refreshed complete bundle.
    /// Source selection and domain checks belong to the material refresh owner;
    /// this method retains the exact key, basis and any declared phase axis.
    pub fn with_material(
        &self,
        material: ResolvedMaterialStatePoint,
    ) -> Result<Self, ConductorError> {
        self.with_material_and_coordinates(material, &[])
    }

    fn with_material_and_coordinates(
        &self,
        material: ResolvedMaterialStatePoint,
        additional: &[(&str, QuantitySpec)],
    ) -> Result<Self, ConductorError> {
        let next = Self::try_from_material_with_property(
            material, self.source_basis, &self.source_key,
        )?;
        match self.liquid_mass_fraction_axis() {
            Some(axis) => next.with_phase_axis_and_coordinates(axis.to_owned(), additional),
            None => Ok(next),
        }
    }

    /// Original selected property, including support, uncertainty and receipt.
    /// A derived resistivity is never inserted into the material bundle.
    #[must_use]
    pub fn source_property(&self) -> &ResolvedScalarProperty {
        self.material.property_by_key(&self.source_key)
            .expect("the immutable conductor retains its admitted source property")
    }

    /// Isotropic bulk resistivity [ohm m], suitable for spatial current balance.
    #[must_use]
    pub const fn resistivity_ohm_m(&self) -> f64 {
        self.resistivity_ohm_m
    }

    /// Lower this same resolved state to a uniform lumped resistor.
    ///
    /// Geometry is supplied at the resolved state; no thermal expansion or
    /// contact resistance is inferred. The source bundle is moved intact.
    ///
    /// # Errors
    /// Non-finite/non-positive geometry or unrepresentable resistance.
    pub fn into_uniform_conductor(
        self,
        length_m: f64,
        area_m2: f64,
    ) -> Result<ResolvedConductor, ConductorError> {
        validate_geometry(length_m, area_m2)?;
        let resistance_ohm =
            uniform_resistance_ohm(self.resistivity_ohm_m, length_m, area_m2)?;
        Ok(ResolvedConductor {
            material: self.material,
            source_basis: self.source_basis,
            source_key: self.source_key,
            length_m,
            area_m2,
            resistance_ohm,
        })
    }
}

fn validate_liquid_mass_fraction_axis(
    material: &ResolvedMaterialStatePoint,
    axis: &str,
) -> Result<(), ConductorError> {
    let point = material.query();
    if axis.is_empty()
        || point.axis_quantities().get(axis) != Some(&QuantitySpec::dimensional(Dims([0; 6])))
        || !point
            .axes()
            .get(axis)
            .is_some_and(|fraction| fraction.is_finite() && (0.0..=1.0).contains(fraction))
    {
        return Err(ConductorError::InvalidInput {
            quantity: "liquid mass fraction needs one existing dimensional scalar coordinate in [0,1]",
        });
    }
    Ok(())
}

/// Resolve bulk resistivity without inventing a unit wire or solver geometry.
///
/// Only `electrical_resistivity` is requested. Its dimensions, numerical domain,
/// claim selection, complete query point and source support are checked by the
/// common material state-point owner. A new temperature, pressure or composition
/// requires a new query; missing support is a refusal, not a fallback metal.
///
/// # Errors
/// Missing, ambiguous, wrong-unit, inadmissible or out-of-domain properties, or
/// an incomplete/foreign property pin plan. No partial state is returned.
pub fn resolve_conductor_state(
    card: &MaterialCard,
    point: &QueryPoint,
    selection: MaterialPropertySelection,
) -> Result<ResolvedConductorState, ConductorError> {
    resolve_conductor_state_with_basis(card, point, selection, ElectricalPropertyBasis::Resistivity)
}

/// Resolve an explicitly chosen resistivity or conductivity source at ONE
/// common material state. For conductivity, interpolate/evaluate the original
/// source first, then derive its reciprocal. Interpolating reciprocal samples
/// would be a different material law and is never done here.
///
/// The exact source requirements, pins and uncertainty are retained. An absent
/// or inadmissible chosen source cannot fall back to the other representation.
///
/// # Errors
/// All source/selection refusals from [`resolve_conductor_state`] and an
/// unrepresentable reciprocal. No synthetic property claims are produced.
pub fn resolve_conductor_state_with_basis(
    card: &MaterialCard,
    point: &QueryPoint,
    selection: MaterialPropertySelection,
    source_basis: ElectricalPropertyBasis,
) -> Result<ResolvedConductorState, ConductorError> {
    let key = PropertyKey::new(source_basis.property_name(), source_basis.dimensions());
    resolve_conductor_state_with_property(card, point, selection, source_basis, &key)
}

/// Resolve one explicitly named electrical source without rewriting its card.
/// The full caller-supplied point and selection policy are checked by the same
/// state-point owner as canonical sources. Conductivity is evaluated before
/// taking its reciprocal; no temperature law or missing condition is inferred.
///
/// # Errors
/// Wrong basis/key quantities, source/selection/domain failures or an
/// unrepresentable reciprocal. No partial state or synthetic claim is returned.
pub fn resolve_conductor_state_with_property(
    card: &MaterialCard,
    point: &QueryPoint,
    selection: MaterialPropertySelection,
    source_basis: ElectricalPropertyBasis,
    source_key: &PropertyKey,
) -> Result<ResolvedConductorState, ConductorError> {
    let requirement = source_basis.requirement_for(source_key)?;
    let material = resolve_material_state_point(card, point, &[requirement], selection)?;
    ResolvedConductorState::try_from_material_with_property(material, source_basis, source_key)
}

/// A uniform conductor resolved from one source-supported material state.
#[derive(Clone, Debug, PartialEq)]
pub struct ResolvedConductor {
    material: ResolvedMaterialStatePoint,
    source_basis: ElectricalPropertyBasis,
    source_key: PropertyKey,
    length_m: f64,
    area_m2: f64,
    resistance_ohm: f64,
}

impl ResolvedConductor {
    /// Resolve the signed current for an ideal current or voltage source.
    /// The complementary variable is algebraic; no circuit dynamics are solved.
    pub fn current_a_for_drive(&self, drive: OhmicDrive) -> Result<f64, ConductorError> {
        current_a_for_drive(self.resistance_ohm, drive)
    }
    /// Exact source-resolved material bundle, including its usage receipt.
    #[must_use]
    pub const fn material(&self) -> &ResolvedMaterialStatePoint {
        &self.material
    }

    /// Independent source basis used to derive this conductor's resistance.
    #[must_use]
    pub const fn source_basis(&self) -> ElectricalPropertyBasis { self.source_basis }

    /// Complete original source key, retained through geometric lowering.
    #[must_use]
    pub const fn source_key(&self) -> &PropertyKey { &self.source_key }

    /// Refresh the electrical source while retaining this exact key and the
    /// caller-owned geometry. No thermal expansion or phase evolution is added.
    pub fn with_material(
        &self,
        material: ResolvedMaterialStatePoint,
    ) -> Result<Self, ConductorError> {
        ResolvedConductorState::try_from_material_with_property(
            material, self.source_basis, &self.source_key,
        )?.into_uniform_conductor(self.length_m, self.area_m2)
    }

    /// The original resistivity or conductivity source, never a synthetic claim.
    #[must_use]
    pub fn source_property(&self) -> &ResolvedScalarProperty {
        self.material.property_by_key(&self.source_key)
            .expect("the immutable conductor retains its admitted source property")
    }

    /// Uniform DC resistance [ohm].
    #[must_use]
    pub const fn resistance_ohm(&self) -> f64 {
        self.resistance_ohm
    }

    /// Caller-owned uniform conductor length [m].
    #[must_use]
    pub const fn length_m(&self) -> f64 {
        self.length_m
    }

    /// Caller-owned uniform conductor cross-sectional area [m²].
    #[must_use]
    pub const fn area_m2(&self) -> f64 {
        self.area_m2
    }

    /// Resistive Joule power `I²R` [W] for a finite signed current [A].
    ///
    /// Current direction does not affect the dissipated power. This is a
    /// scalar constitutive evaluation, not a circuit or energy-balance solve.
    pub fn joule_power_w(&self, current_a: f64) -> Result<f64, ConductorError> {
        joule_power_w(self.resistance_ohm, current_a)
    }
}

/// Refusal from source-resolved conductor construction or evaluation.
#[derive(Clone, Debug, PartialEq)]
pub enum ConductorError {
    /// The material card could not supply one admissible resistivity claim.
    MaterialState(MaterialStatePointError),
    /// A caller-owned geometric or electrical input was non-finite or invalid.
    InvalidInput {
        /// Stable input name.
        quantity: &'static str,
    },
    /// A finite input combination could not produce a finite physical result.
    InvalidDerived {
        /// Stable result name.
        quantity: &'static str,
    },
}

impl fmt::Display for ConductorError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::MaterialState(source) => {
                write!(formatter, "conductor material refusal: {source}")
            }
            Self::InvalidInput { quantity } => {
                write!(formatter, "invalid conductor input: {quantity}")
            }
            Self::InvalidDerived { quantity } => {
                write!(formatter, "unrepresentable conductor result: {quantity}")
            }
        }
    }
}

impl std::error::Error for ConductorError {}

impl From<MaterialStatePointError> for ConductorError {
    fn from(source: MaterialStatePointError) -> Self {
        Self::MaterialState(source)
    }
}

fn validate_geometry(length_m: f64, area_m2: f64) -> Result<(), ConductorError> {
    if !(length_m.is_finite() && length_m > 0.0) {
        return Err(ConductorError::InvalidInput {
            quantity: "length_m",
        });
    }
    if !(area_m2.is_finite() && area_m2 > 0.0) {
        return Err(ConductorError::InvalidInput {
            quantity: "area_m2",
        });
    }
    Ok(())
}

fn uniform_resistance_ohm(
    resistivity_ohm_m: f64,
    length_m: f64,
    area_m2: f64,
) -> Result<f64, ConductorError> {
    let direct = resistivity_ohm_m * length_m / area_m2;
    let resistance_ohm = if direct.is_finite() && direct > 0.0 {
        direct
    } else {
        resistivity_ohm_m * (length_m / area_m2)
    };
    if !(resistance_ohm.is_finite() && resistance_ohm > 0.0) {
        return Err(ConductorError::InvalidDerived {
            quantity: "resistance_ohm",
        });
    }
    Ok(resistance_ohm)
}

fn joule_power_w(resistance_ohm: f64, current_a: f64) -> Result<f64, ConductorError> {
    if !current_a.is_finite() {
        return Err(ConductorError::InvalidInput {
            quantity: "current_a",
        });
    }
    let direct = current_a * current_a * resistance_ohm;
    let power_w = if direct.is_finite() && (direct > 0.0 || current_a == 0.0) {
        direct
    } else {
        (current_a * resistance_ohm) * current_a
    };
    if !(power_w.is_finite() && (power_w > 0.0 || current_a == 0.0)) {
        return Err(ConductorError::InvalidDerived {
            quantity: "joule_power_w",
        });
    }
    Ok(power_w)
}

fn current_a_for_drive(resistance_ohm: f64, drive: OhmicDrive) -> Result<f64, ConductorError> {
    let (value, quantity) = match drive {
        OhmicDrive::Current(value) => (value, "current_a"),
        OhmicDrive::Voltage(value) => (value, "voltage_v"),
    };
    if !value.is_finite() {
        return Err(ConductorError::InvalidInput { quantity });
    }
    let current = match drive {
        OhmicDrive::Current(value) => value,
        OhmicDrive::Voltage(value) => value / resistance_ohm,
    };
    if !current.is_finite() || (value != 0.0 && current == 0.0) {
        return Err(ConductorError::InvalidDerived {
            quantity: "current_a",
        });
    }
    Ok(current)
}

/// Resolve a uniform isotropic DC conductor from one resistivity claim.
///
/// The card must provide `electrical_resistivity` in ohm metre dimensions at
/// the complete caller-supplied query point. Geometry is the actual length and
/// area at that state, not reference geometry to be thermally expanded here;
/// no contact resistance, current distribution, circuit topology, temperature
/// update, or transient energy balance is inferred.
///
/// # Errors
/// Refuses invalid geometry, missing/ambiguous/wrong-unit/out-of-domain source
/// claims, incomplete pin plans, and non-finite resistance.
pub fn resolve_uniform_conductor(
    card: &MaterialCard,
    point: &QueryPoint,
    selection: MaterialPropertySelection,
    length_m: f64,
    area_m2: f64,
) -> Result<ResolvedConductor, ConductorError> {
    // Preserve the established geometry-first admission order while delegating
    // all material resolution and lowering to the common bulk-state path.
    validate_geometry(length_m, area_m2)?;
    resolve_conductor_state(card, point, selection)?.into_uniform_conductor(length_m, area_m2)
}

/// Resolve and geometrically lower an explicitly selected electrical source.
/// Shares exactly the same geometry, ideal-drive and Joule-power owners as
/// [`resolve_uniform_conductor`]. Only the independent source quantity differs.
///
/// # Errors
/// Invalid geometry, source/selection failures or unrepresentable derivation.
pub fn resolve_uniform_conductor_with_basis(
    card: &MaterialCard,
    point: &QueryPoint,
    selection: MaterialPropertySelection,
    source_basis: ElectricalPropertyBasis,
    length_m: f64,
    area_m2: f64,
) -> Result<ResolvedConductor, ConductorError> {
    validate_geometry(length_m, area_m2)?;
    resolve_conductor_state_with_basis(card, point, selection, source_basis)?
        .into_uniform_conductor(length_m, area_m2)
}

/// Resolve an exact source key and use the existing uniform resistance and
/// Joule-power owners at the caller's physical geometry.
///
/// # Errors
/// Invalid geometry, wrong source quantity, material query refusals or an
/// unrepresentable resistance. No contact or temperature model is inferred.
pub fn resolve_uniform_conductor_with_property(
    card: &MaterialCard,
    point: &QueryPoint,
    selection: MaterialPropertySelection,
    source_basis: ElectricalPropertyBasis,
    source_key: &PropertyKey,
    length_m: f64,
    area_m2: f64,
) -> Result<ResolvedConductor, ConductorError> {
    validate_geometry(length_m, area_m2)?;
    resolve_conductor_state_with_property(card, point, selection, source_basis, source_key)?
        .into_uniform_conductor(length_m, area_m2)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn g1_ideal_voltage_and_current_drives_obey_ohms_law_and_passivity() {
        for resistance in [0.25, 10.0, 1000.0] {
            for voltage in [-12.0, 0.0, 12.0] {
                let current =
                    current_a_for_drive(resistance, OhmicDrive::Voltage(voltage)).unwrap();
                assert_eq!(current, voltage / resistance);
                assert_eq!(
                    current_a_for_drive(resistance, OhmicDrive::Current(current)).unwrap(),
                    current
                );
                let power = joule_power_w(resistance, current).unwrap();
                assert!(power >= 0.0);
                assert!((power - voltage * current).abs() <= 4.0 * f64::EPSILON * power);
            }
        }
    }

    #[test]
    fn g0_voltage_drive_refuses_nonfinite_and_unrepresentable_current() {
        for voltage in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            assert!(matches!(
                current_a_for_drive(1.0, OhmicDrive::Voltage(voltage)),
                Err(ConductorError::InvalidInput {
                    quantity: "voltage_v"
                })
            ));
        }
        for (resistance, voltage) in [(f64::MIN_POSITIVE, f64::MAX), (f64::MAX, f64::MIN_POSITIVE)]
        {
            assert!(matches!(
                current_a_for_drive(resistance, OhmicDrive::Voltage(voltage)),
                Err(ConductorError::InvalidDerived {
                    quantity: "current_a"
                })
            ));
        }
    }

    #[test]
    fn g0_resistance_recovers_after_intermediate_product_overflow() {
        let resistance = uniform_resistance_ohm(1.0e200, 1.0e200, 1.0e200).unwrap();
        assert_eq!(resistance, 1.0e200);
    }

    #[test]
    fn g0_joule_power_recovers_after_current_square_underflow() {
        let power = joule_power_w(1.0e200, 1.0e-200).unwrap();
        assert_eq!(power, 1.0e-200);
    }

    #[test]
    fn g0_joule_power_refuses_unrepresentable_nonzero_result() {
        assert_eq!(
            joule_power_w(f64::MIN_POSITIVE, f64::MIN_POSITIVE),
            Err(ConductorError::InvalidDerived {
                quantity: "joule_power_w"
            })
        );
    }
}
