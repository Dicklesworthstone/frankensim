//! Reversible thermal response of source-backed linear field storage.
//!
//! At fixed geometry and charge q or flux linkage lambda, K is respectively
//! C(T) or L(T). The electrical owner uses the Helmholtz contribution
//! psi = z^2/(2 K). This module supplies s = -dpsi/dT and U = psi + T s;
//! using psi alone as total stored energy would omit electro/magnetocaloric
//! exchange. For an affine positive K(T), s = psi K'/K and
//! (dU/dT)_z = -2 T psi (K'/K)^2. The CALORIC plus field heat capacity,
//! not either contribution alone, must remain positive in a coupled owner.
//!
//! This is the fixed-geometry, isotropic, nondispersive, linear reversible rung.
//! No Curie transition, saturation, hysteresis or loss is inferred from a
//! permeability/permittivity curve. Derivative jumps require entropy jumps and
//! a different caloric chart; globally non-affine source paths are refused.
//! Standard thermodynamic convention: U = F + TS, S = -(dF/dT) at fixed work
//! coordinates (MIT 2.43, lecture 6, Helmholtz free energy).
//! <https://web.mit.edu/beretta/www/AdvancedThermodynamics_2_43_Lecture6.htm>

use fs_matdb::MaterialCard;
use fs_qty::{
    QuantitySpec,
    semantic::{QuantityKind, SemanticType, ValueForm},
};

use super::{
    ELECTRIC_PERMITTIVITY_PROPERTY, FieldStorageError, MAGNETIC_PERMEABILITY_PROPERTY,
    ResolvedCapacitor, ResolvedInductor, energy, positive,
};
use crate::state_point::{
    ResolvedMaterialStatePoint, affine_temperature_slope, pinned_temperature_line,
};

/// Spatial material field thermodynamics, using the same affine source admission
/// and potential derivatives as lumped C(T)/L(T).
pub mod magnetic;
/// Spatial scalar or oriented dielectric thermodynamics at fixed displacement.
pub mod electric;

/// A native material storage element. Geometry and source receipts are retained;
/// this is not another circuit state, graph branch or material evaluator.
#[derive(Clone, Debug, PartialEq)]
pub enum ThermalFieldElement {
    /// q is signed charge [C], K is capacitance [F].
    Capacitor(ResolvedCapacitor),
    /// z is signed flux LINKAGE [Wb], K is inductance [H].
    Inductor(ResolvedInductor),
}
impl ThermalFieldElement {
    /// Complete material snapshot shared with the original field owner.
    #[must_use]
    pub fn material(&self) -> &ResolvedMaterialStatePoint {
        match self {
            Self::Capacitor(value) => value.material(),
            Self::Inductor(value) => value.material(),
        }
    }
    /// Capacitance [F] or inductance [H], according to the enum variant.
    #[must_use]
    pub fn coefficient(&self) -> f64 {
        match self {
            Self::Capacitor(value) => value.capacitance_f(),
            Self::Inductor(value) => value.inductance_h(),
        }
    }
    fn property_name(&self) -> &'static str {
        match self {
            Self::Capacitor(_) => ELECTRIC_PERMITTIVITY_PROPERTY,
            Self::Inductor(_) => MAGNETIC_PERMEABILITY_PROPERTY,
        }
    }
    fn with_material(
        &self,
        material: ResolvedMaterialStatePoint,
    ) -> Result<Self, FieldStorageError> {
        match self {
            Self::Capacitor(value) => {
                ResolvedCapacitor::try_from_material(material, value.area_m2(), value.gap_m())
                    .map(Self::Capacitor)
            }
            Self::Inductor(value) => ResolvedInductor::try_from_material(
                material,
                value.area_m2(),
                value.path_length_m(),
                value.turns(),
            )
            .map(Self::Inductor),
        }
    }
}

/// Immutable affine-in-temperature coefficient law on a complete supported path.
/// The original material card, selected claim, nonthermal coordinates, and
/// geometry are fixed. A different material or physical context is a new law,
/// not an unaccounted replacement of live electrical or caloric history.
#[derive(Clone, Debug, PartialEq)]
pub struct ThermalFieldStorage {
    card: MaterialCard,
    original: ThermalFieldElement,
    temperature_axis: String,
    temperature_range_k: [f64; 2],
    original_temperature_k: f64,
    coefficient_slope_per_k: f64,
    minimum_coefficient: f64,
}

/// One exact source refresh and its derivative in the admitted affine model.
#[derive(Clone, Debug, PartialEq)]
pub struct ResolvedThermalFieldStorage {
    element: ThermalFieldElement,
    temperature_k: f64,
    coefficient_slope_per_k: f64,
}

/// Thermal potentials of one scalar field element, or a documented signed
/// branch partition supplied by a coupled owner at fixed total coordinates.
/// The signed entropy and internal-energy contributions are relative to the
/// zero-field material caloric chart; they are not separately positive stores.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FieldThermalState {
    /// Electrical Helmholtz contribution [J]: z^2/(2 K) for a scalar element;
    /// coupled branch partitions can be signed and must be summed by group.
    pub free_energy_j: f64,
    /// Field entropy contribution -dpsi/dT at fixed z [J/K].
    pub entropy_j_k: f64,
    /// Field internal-energy contribution psi + T s [J].
    pub internal_energy_j: f64,
    /// Fixed-z field heat-capacity contribution dU/dT [J/K], nonpositive for
    /// this scalar affine law. Coupled consumers may report signed branch
    /// partitions; only their complete group's sum has that interpretation.
    pub heat_capacity_j_k: f64,
}

fn finite(value: f64, what: &'static str) -> Result<f64, FieldStorageError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(FieldStorageError::Invalid(what))
    }
}

// Shared source-path admission. The model derivative comes from the exact
// selected source path, never from finite differences of unrelated queries.
fn affine_path(
    card: &MaterialCard,
    material: &ResolvedMaterialStatePoint,
    property_name: &str,
    temperature_axis: &str,
    temperature_range_k: [f64; 2],
) -> Result<(f64, f64, f64), FieldStorageError> {
    let (reference, slope, endpoints) = affine_source_path(
        card, material, property_name, temperature_axis, temperature_range_k, true,
    )?;
    Ok((reference, slope, endpoints[0].min(endpoints[1])))
}

/// Exact affine source-path admission shared by scalar storage and coupled
/// active-material blocks. Signed cross coefficients need finite support;
/// diagonal storage coefficients additionally require positivity.
pub(crate) fn affine_source_path(
    card: &MaterialCard,
    material: &ResolvedMaterialStatePoint,
    property_name: &str,
    temperature_axis: &str,
    temperature_range_k: [f64; 2],
    strictly_positive: bool,
) -> Result<(f64, f64, [f64; 2]), FieldStorageError> {
    let invalid = FieldStorageError::Invalid;
    let [low, high] = temperature_range_k;
    if !low.is_finite() || !high.is_finite() || low <= 0.0 || low >= high {
        return Err(invalid("field thermal range must be finite positive increasing kelvin"));
    }
    let context = material.query();
    let absolute = QuantitySpec::semantic(SemanticType::new(
        QuantityKind::AbsoluteTemperature, ValueForm::Static,
    ));
    if context.axis_quantities().get(temperature_axis) != Some(&absolute) {
        return Err(invalid("field thermal coordinate needs explicit absolute-temperature semantics"));
    }
    let reference_k = *context.axes().get(temperature_axis)
        .ok_or(invalid("field thermal coordinate is missing"))?;
    if !(low..=high).contains(&reference_k) {
        return Err(invalid("original field temperature lies outside its thermal range"));
    }
    material.at_query_point(card, &context)?;
    let property = material.property(property_name)
        .ok_or(invalid("field thermal property is missing or ambiguous"))?;
    if !property.has_only_fixed_external_conditions(Some(temperature_axis), false) {
        return Err(invalid("field law requires one owned temperature and fixed external source coordinates"));
    }
    let point = |temperature| context.clone().with(temperature_axis, temperature)
        .map_err(|_| invalid("field thermal query is unrepresentable"));
    let claim = card.claims().claim(property.answer().receipt.selected)
        .ok_or(invalid("field source claim is absent from original card"))?;
    let knots = pinned_temperature_line(card, &claim.key, temperature_axis,
        &point(low)?, &point(high)?, property.answer().receipt.selected)?;
    if knots.first().map(|v| v.0) != Some(low)
        || knots.last().map(|v| v.0) != Some(high)
        || knots.iter().any(|(_, value)| !value.is_finite() || strictly_positive && *value <= 0.0)
    {
        return Err(invalid("admissible field source must continuously cover the complete thermal range"));
    }
    // Positive storage retains its original relative admission. Signed cross
    // coefficients can cancel exactly at the reference temperature; admit their
    // reference against the source-path scale rather than divide by that zero.
    let (slope_reference_k, slope_reference_value) = if strictly_positive {
        (reference_k, property.value_si())
    } else { knots[0] };
    let slope = affine_temperature_slope(&knots, slope_reference_k, slope_reference_value).map_err(invalid)?;
    if !strictly_positive {
        let expected = finite(slope.mul_add(reference_k-knots[0].0,knots[0].1),
            "signed affine source reference is unrepresentable")?;
        let scale = knots.iter().fold(property.value_si().abs().max(expected.abs()),
            |scale,(_,value)| scale.max(value.abs()));
        if (expected-property.value_si()).abs()>128.0*f64::EPSILON*scale {
            return Err(invalid("signed coefficient reference differs from its admitted affine path"));
        }
    }
    Ok((reference_k, slope, [knots.first().unwrap().1, knots.last().unwrap().1]))
}

fn thermal_state(
    free_energy_j: f64,
    temperature_k: f64,
    ratio: f64,
) -> Result<FieldThermalState, FieldStorageError> {
    if ratio == 0.0 || free_energy_j == 0.0 {
        return Ok(FieldThermalState {
            free_energy_j, entropy_j_k: 0.0, internal_energy_j: free_energy_j,
            heat_capacity_j_k: 0.0,
        });
    }
    let entropy_j_k = finite(free_energy_j * ratio, "field entropy is unrepresentable")?;
    let internal_energy_j = finite(temperature_k.mul_add(entropy_j_k, free_energy_j),
        "field internal energy is unrepresentable")?;
    let heat_capacity_j_k = finite((-2.0*temperature_k*entropy_j_k)*ratio,
        "field heat capacity is unrepresentable")?;
    if entropy_j_k == 0.0 || heat_capacity_j_k == 0.0 {
        return Err(FieldStorageError::Invalid("nonzero field thermal contribution underflows"));
    }
    Ok(FieldThermalState { free_energy_j, entropy_j_k, internal_energy_j, heat_capacity_j_k })
}

pub(crate) fn capacity_bound(
    worst_energy_j: f64,
    maximum_temperature_k: f64,
    maximum_ratio: f64,
) -> Result<f64, FieldStorageError> {
    if worst_energy_j == 0.0 || maximum_ratio == 0.0 { return Ok(0.0); }
    let magnitude = finite((2.0*maximum_temperature_k*worst_energy_j*maximum_ratio.abs())
        *maximum_ratio.abs(), "field capacity bound is unrepresentable")?;
    if magnitude == 0.0 {
        return Err(FieldStorageError::Invalid("nonzero field capacity bound underflows"));
    }
    finite(-magnitude*(1.0+2048.0*f64::EPSILON), "field capacity bound is unrepresentable")
}

impl ThermalFieldStorage {
    /// Bind a native capacitor to its original card over a positive temperature
    /// interval. The original absolute-temperature coordinate must be explicitly
    /// typed, and its value must lie in the supplied interval.
    pub fn capacitor(
        card: &MaterialCard,
        original: &ResolvedCapacitor,
        temperature_axis: &str,
        temperature_range_k: [f64; 2],
    ) -> Result<Self, FieldStorageError> {
        Self::new(
            card,
            ThermalFieldElement::Capacitor(original.clone()),
            temperature_axis,
            temperature_range_k,
        )
    }

    /// Bind the native uniform magnetic-path element. This does not extend its
    /// geometric reduction to moving gaps, leakage or spatial field solutions.
    pub fn inductor(
        card: &MaterialCard,
        original: &ResolvedInductor,
        temperature_axis: &str,
        temperature_range_k: [f64; 2],
    ) -> Result<Self, FieldStorageError> {
        Self::new(
            card,
            ThermalFieldElement::Inductor(original.clone()),
            temperature_axis,
            temperature_range_k,
        )
    }

    fn new(
        card: &MaterialCard,
        original: ThermalFieldElement,
        temperature_axis: &str,
        temperature_range_k: [f64; 2],
    ) -> Result<Self, FieldStorageError> {
        let invalid = FieldStorageError::Invalid;
        let material = original.material();
        let (reference_k, source_slope, _) = affine_path(card, material,
            original.property_name(), temperature_axis, temperature_range_k)?;
        let [low, high] = temperature_range_k;
        let context = material.query();
        let property = material
            .property(original.property_name())
            .ok_or(invalid("field thermal property is missing or ambiguous"))?;
        let point = |temperature| {
            context
                .clone()
                .with(temperature_axis, temperature)
                .map_err(|_| invalid("field thermal query is unrepresentable"))
        };
        let lower = point(low)?;
        let upper = point(high)?;
        // Geometry is fixed, so the logarithmic derivative of K equals that
        // of its source epsilon or mu; no finite-difference source queries.
        let scale = positive(
            original.coefficient() / property.value_si(),
            "field geometry scale is unrepresentable",
        )?;
        let coefficient_slope_per_k = finite(
            source_slope * scale,
            "field coefficient derivative is unrepresentable",
        )?;
        if source_slope != 0.0 && coefficient_slope_per_k == 0.0 {
            return Err(invalid("field coefficient derivative underflows"));
        }
        let at_lower = original.with_material(material.at_query_point(card, &lower)?)?;
        let at_upper = original.with_material(material.at_query_point(card, &upper)?)?;
        let minimum_coefficient = at_lower.coefficient().min(at_upper.coefficient());
        Ok(Self {
            card: card.clone(),
            original,
            temperature_axis: temperature_axis.to_owned(),
            temperature_range_k,
            original_temperature_k: reference_k,
            coefficient_slope_per_k,
            minimum_coefficient,
        })
    }

    /// Original native geometry and complete field material snapshot.
    #[must_use]
    pub const fn original(&self) -> &ThermalFieldElement {
        &self.original
    }
    /// Same immutable material card used for every temperature query.
    #[must_use]
    pub const fn card(&self) -> &MaterialCard {
        &self.card
    }
    /// Complete original physical context and selected source receipts.
    #[must_use]
    pub fn material(&self) -> &ResolvedMaterialStatePoint {
        self.original.material()
    }
    /// The one changing physical coordinate.
    #[must_use]
    pub fn temperature_axis(&self) -> &str {
        &self.temperature_axis
    }
    /// Entire admitted continuous source interval [K].
    #[must_use]
    pub const fn temperature_range_k(&self) -> [f64; 2] {
        self.temperature_range_k
    }
    /// Absolute temperature of the retained original snapshot [K].
    #[must_use]
    pub const fn original_temperature_k(&self) -> f64 {
        self.original_temperature_k
    }
    /// Derived dC/dT [F/K] or dL/dT [H/K] with fixed geometry.
    #[must_use]
    pub const fn coefficient_slope_per_k(&self) -> f64 {
        self.coefficient_slope_per_k
    }

    /// Re-evaluate the ENTIRE original source bundle at one supported temperature.
    /// All claims, quantity conventions and other physical coordinates remain
    /// pinned. The result supplies the native C/L and the thermal potentials.
    pub fn at_temperature(
        &self,
        temperature_k: f64,
    ) -> Result<ResolvedThermalFieldStorage, FieldStorageError> {
        let invalid = FieldStorageError::Invalid;
        let [low, high] = self.temperature_range_k;
        if !temperature_k.is_finite() || !(low..=high).contains(&temperature_k) {
            return Err(invalid(
                "field temperature is outside its admitted source interval",
            ));
        }
        let point = self
            .material()
            .query()
            .with(&self.temperature_axis, temperature_k)
            .map_err(|_| invalid("field thermal query is unrepresentable"))?;
        let element = self
            .original
            .with_material(self.material().at_query_point(&self.card, &point)?)?;
        let expected = finite(
            self.coefficient_slope_per_k.mul_add(
                temperature_k - self.original_temperature_k,
                self.original.coefficient(),
            ),
            "affine field coefficient is unrepresentable",
        )?;
        if (element.coefficient() - expected).abs()
            > 512.0 * f64::EPSILON * element.coefficient().abs().max(expected.abs())
        {
            return Err(invalid(
                "refreshed source left the admitted affine field law",
            ));
        }
        Ok(ResolvedThermalFieldStorage {
            element,
            temperature_k,
            coefficient_slope_per_k: self.coefficient_slope_per_k,
        })
    }

    /// Conservative lower bound on the field contribution to dU/dT [J/K]
    /// throughout the admitted interval for |q| or |lambda| <= the declared bound.
    /// Uses Tmax and Kmin separately, so it can be stricter than the actual
    /// minimum. The coupled owner adds ALL fields in a cell to its caloric
    /// capacity, requires a positive total, and enforces the coordinate bound.
    /// This is analytic model admission with roundoff padding, not an interval
    /// certificate or a bound on nonlinear/time-discretization error.
    pub fn minimum_heat_capacity_j_k(
        &self,
        maximum_absolute_coordinate: f64,
    ) -> Result<f64, FieldStorageError> {
        let invalid = FieldStorageError::Invalid;
        if !maximum_absolute_coordinate.is_finite() || maximum_absolute_coordinate < 0.0 {
            return Err(invalid(
                "field coordinate bound must be finite and nonnegative",
            ));
        }
        if maximum_absolute_coordinate == 0.0 || self.coefficient_slope_per_k == 0.0 {
            return Ok(0.0);
        }
        let worst_energy = energy(maximum_absolute_coordinate, self.minimum_coefficient)?;
        let ratio = self.coefficient_slope_per_k / self.minimum_coefficient;
        if ratio == 0.0 || !ratio.is_finite() {
            return Err(invalid("field capacity derivative ratio is unrepresentable"));
        }
        capacity_bound(worst_energy, self.temperature_range_k[1], ratio)
    }
}

impl ResolvedThermalFieldStorage {
    /// Native C/L snapshot for the original circuit owner.
    #[must_use]
    pub const fn element(&self) -> &ThermalFieldElement {
        &self.element
    }
    /// Actual absolute temperature [K].
    #[must_use]
    pub const fn temperature_k(&self) -> f64 {
        self.temperature_k
    }
    /// Actual field Helmholtz, entropy, internal-energy and capacity contributions.
    pub fn state(&self, coordinate: f64) -> Result<FieldThermalState, FieldStorageError> {
        let free_energy_j = energy(coordinate, self.element.coefficient())?;
        let ratio = finite(
            self.coefficient_slope_per_k / self.element.coefficient(),
            "field thermal derivative ratio is unrepresentable",
        )?;
        if self.coefficient_slope_per_k != 0.0 && ratio == 0.0 {
            return Err(FieldStorageError::Invalid(
                "nonzero field thermal derivative underflows",
            ));
        }
        thermal_state(free_energy_j, self.temperature_k, ratio)
    }
}

