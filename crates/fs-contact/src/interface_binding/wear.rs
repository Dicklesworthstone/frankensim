//! Original ordered-interface wear coefficients and per-surface hardness.
//! The caller supplies actual unsigned sliding distance and compressive load.
//! Temperature and claims remain with the original material owners; this layer
//! advances native Archard volume without inventing recession or an energy sink.

use super::BoundDryInterfaceSystem;
use super::thermal_friction::elliptic::SourceThermalEllipticCoulomb;
use fs_blake3::ContentHash;
use fs_matdb::{MaterialCard, PropertyKey, QueryPoint};
use fs_material::phase::EquilibriumPhaseState;
use fs_material::state_point::{ResolvedInterfaceStatePoint, ResolvedMaterialStatePoint};
use fs_material::thermal::source::{ThermalMaterialError, ThermalMaterialPoint};
use fs_qty::semantic::{QuantityKind, ValueForm};
use fs_qty::{Dims, Pressure, QuantitySpec};
use fs_tribo::{ArchardLaw, TriboError, WearAdvance, WearState};

#[derive(Clone, Debug, PartialEq)]
pub enum SourceWearError {
    Invalid(&'static str),
    Material(ThermalMaterialError),
    Tribology(TriboError),
}
impl From<ThermalMaterialError> for SourceWearError {
    fn from(error: ThermalMaterialError) -> Self {
        Self::Material(error)
    }
}
impl From<TriboError> for SourceWearError {
    fn from(error: TriboError) -> Self {
        Self::Tribology(error)
    }
}
impl core::fmt::Display for SourceWearError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "source Archard wear: {self:?}")
    }
}
impl std::error::Error for SourceWearError {}

/// Retained wear and original selected claims, in immutable surface A/B order.
/// Coefficients are dimensionless Archard k; hardness is actual indentation
/// hardness in Pa, not a tensile strength or an unconverted hardness number.
#[derive(Clone, Debug, PartialEq)]
pub struct SourceArchardWear {
    friction_owner: ContentHash,
    interface: BoundDryInterfaceSystem,
    coefficients: ResolvedInterfaceStatePoint,
    coefficient_keys: [PropertyKey; 2],
    hardness: [ResolvedMaterialStatePoint; 2],
    hardness_keys: [PropertyKey; 2],
    temperature_axes: [String; 2],
    volumes: [WearState; 2],
}

impl SourceArchardWear {
    /// Bind actual receipts to the SAME source friction owner and ordered bulk
    /// materials. Both roles may explicitly select the same interface claim.
    /// Hardness test context is retained in the complete property keys. Only
    /// declared friction coordinates and each bulk temperature evolve; other
    /// supported external conditions stay fixed. No source authority is raised.
    pub fn new(
        friction: &SourceThermalEllipticCoulomb,
        coefficients: ResolvedInterfaceStatePoint,
        coefficient_keys: [PropertyKey; 2],
        hardness: [ResolvedMaterialStatePoint; 2],
        hardness_keys: [PropertyKey; 2],
        temperature_axes: [String; 2],
        volumes: [WearState; 2],
    ) -> Result<Self, SourceWearError> {
        let original = &friction.principal_sources()[0];
        let interface = original.interface();
        let query = original.initial_source().query();
        if coefficients.card_identity() != interface.card_identity()
            || coefficients.surface_materials() != interface.surface_materials()
            || coefficients.axis_quantities() != query.axis_quantities()
        {
            return Err(SourceWearError::Invalid(
                "wear requires the same ordered interface and query schema",
            ));
        }
        let owned_names: Vec<_> = std::iter::once(original.axes().temperature.as_str())
            .chain(original.axes().nominal_pressure.as_deref())
            .chain(original.query_axes().iter().map(|axis| axis.name.as_str()))
            .collect();
        let owned: Vec<_> = owned_names
            .iter()
            .map(|name| {
                query
                    .axis_quantities()
                    .get(*name)
                    .copied()
                    .map(|q| (*name, q))
                    .ok_or(SourceWearError::Invalid(
                        "wear requires typed owned coordinates",
                    ))
            })
            .collect::<Result<_, _>>()?;
        if coefficients.query_point().iter().any(|(name, value)| {
            !owned_names.contains(&name.as_str()) && query.axes().get(name) != Some(value)
        }) {
            return Err(SourceWearError::Invalid(
                "wear cannot substitute fixed interface conditions",
            ));
        }
        for side in 0..2 {
            let k = coefficients
                .property_by_key(&coefficient_keys[side])
                .ok_or(SourceWearError::Invalid(
                    "missing exact wear coefficient key",
                ))?;
            let h = hardness[side].property_by_key(&hardness_keys[side]).ok_or(
                SourceWearError::Invalid("missing exact hardness key and test context"),
            )?;
            let hq = hardness_keys[side].quantity();
            let physical_hardness = hq.dims() == Pressure::DIMS
                && hq.semantic_type().is_none_or(|ty| {
                    matches!(ty.kind(), QuantityKind::Hardness(_)) && ty.form() == ValueForm::Static
                });
            if coefficient_keys[side].quantity() != QuantitySpec::dimensional(Dims::NONE)
                || !k.value_si().is_finite()
                || k.value_si() < 0.0
                || !k.has_only_owned_coordinates_and_fixed_external_conditions(&owned)
                || hardness[side].material() != &interface.surface_materials()[side]
                || !physical_hardness
                || !h.value_si().is_finite()
                || h.value_si() <= 0.0
                || !h.has_only_fixed_external_conditions(Some(&temperature_axes[side]), false)
                || !hardness[side]
                    .state_coordinate(&temperature_axes[side])
                    .is_some_and(|t| t.is_finite() && t > 0.0)
            {
                return Err(SourceWearError::Invalid(
                    "wear needs nonnegative dimensionless k, original positive Pa hardness and owned source conditions",
                ));
            }
        }
        Ok(Self {
            friction_owner: friction.identity(),
            interface: interface.clone(),
            coefficients,
            coefficient_keys,
            hardness,
            hardness_keys,
            temperature_axes,
            volumes,
        })
    }

    /// Cumulative removed volume in original A/B roles [m³].
    pub fn volumes_m3(&self) -> [f64; 2] {
        self.volumes.map(|v| v.volume_m3())
    }
    /// Current exact selected interface claims and physical query.
    pub const fn coefficient_source(&self) -> &ResolvedInterfaceStatePoint {
        &self.coefficients
    }
    /// Current original bulk hardness claims, including test contexts.
    pub const fn hardness_sources(&self) -> &[ResolvedMaterialStatePoint; 2] {
        &self.hardness
    }
    /// Check source ownership before doing mechanical or thermal work.
    pub fn admit_friction(
        &self,
        friction: &SourceThermalEllipticCoulomb,
    ) -> Result<(), SourceWearError> {
        if self.friction_owner != friction.identity() {
            return Err(SourceWearError::Invalid(
                "wear checkpoint belongs to another friction source",
            ));
        }
        Ok(())
    }

    /// Resolve the SAME claims at the actual endpoint and advance native wear.
    /// `distance_m` is unsigned sliding arc length, never net displacement or
    /// an impulse-derived distance. N and endpoint k/H are held over this phase;
    /// changing conditions require temporal refinement. Failure mutates nothing.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_candidate(
        &self,
        friction: &SourceThermalEllipticCoulomb,
        interface_query: &QueryPoint,
        materials: [&MaterialCard; 2],
        phases: [EquilibriumPhaseState; 2],
        normal_n: f64,
        distance_m: f64,
    ) -> Result<(Self, [WearAdvance; 2]), SourceWearError> {
        self.admit_friction(friction)?;
        let original = &friction.principal_sources()[0];
        let initial = original.initial_source().query();
        if interface_query.axes().iter().any(|(name, value)| {
            name != &original.axes().temperature
                && original.axes().nominal_pressure.as_ref() != Some(name)
                && !original.query_axes().iter().any(|axis| &axis.name == name)
                && initial.axes().get(name) != Some(value)
        }) {
            return Err(SourceWearError::Invalid(
                "wear cannot evolve unowned interface conditions",
            ));
        }
        let mut next = self.clone();
        next.coefficients = self
            .coefficients
            .at_query_point(self.interface.card(), interface_query)
            .map_err(ThermalMaterialError::Material)?;
        for side in 0..2 {
            if phases[side].liquid_mass_fraction() != 0.0 {
                return Err(SourceWearError::Invalid(
                    "dry Archard wear requires solid surface material",
                ));
            }
            let context = self.hardness[side].query();
            let point = ThermalMaterialPoint {
                material: materials[side],
                context: &context,
                temperature_axis: &self.temperature_axes[side],
                liquid_fraction_axis: None,
            }
            .at(phases[side])?;
            next.hardness[side] = self.hardness[side]
                .at_query_point(materials[side], &point)
                .map_err(ThermalMaterialError::Material)?;
        }
        let mut advance = |side: usize| {
            let k = next
                .coefficients
                .property_by_key(&self.coefficient_keys[side])
                .ok_or(SourceWearError::Invalid(
                    "wear coefficient disappeared on refresh",
                ))?
                .value_si();
            let h = next.hardness[side]
                .property_by_key(&self.hardness_keys[side])
                .ok_or(SourceWearError::Invalid("hardness disappeared on refresh"))?
                .value_si();
            Ok::<_, SourceWearError>(
                ArchardLaw {
                    wear_coefficient: k,
                }
                .advance(
                    self.interface.interface(),
                    &mut next.volumes[side],
                    normal_n,
                    distance_m,
                    h,
                )?,
            )
        };
        let reports = [advance(0)?, advance(1)?];
        Ok((next, reports))
    }
}
