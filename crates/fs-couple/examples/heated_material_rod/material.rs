//! Existing canonical pack ingress and shared caloric/material evaluators.
use std::io::Read;
use fs_evidence::ValidityDomain;
use fs_matdb::{ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId,
    NormalizedMaterialCardPack, PropertyClaim, PropertyKey, PropertyValue,
    Provenance, QueryPoint, SelectionPolicy, UncertaintyModel};
use fs_material::phase::{HeatCapacitySourceQuery, ResolvedHeatCapacityCurve, SolidLiquidPhase};
use fs_material::state_point::{DENSITY_PROPERTY, INVERSE_TEMPERATURE_DIMS,
    LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY, MaterialPropertySelection,
    ResolvedMaterialStatePoint, SPECIFIC_HEAT_CAPACITY_DIMS, SPECIFIC_HEAT_CAPACITY_PROPERTY,
    ScalarAdmissibility, ScalarPropertyRequirement, THERMAL_CONDUCTIVITY_DIMS,
    THERMAL_CONDUCTIVITY_PROPERTY, resolve_material_state_point};
use fs_qty::{Density, Pressure, QuantitySpec, Temperature};
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use super::Result;

pub(super) struct Material {
    pub card: MaterialCard,
    pub caloric: ResolvedHeatCapacityCurve,
    k_requirement: ScalarPropertyRequirement,
    k_selection: MaterialPropertySelection,
}
pub(super) fn point(t: f64) -> Result<QueryPoint> {
    Ok(QueryPoint::new().with("T", t)?)
}
impl Material {
    pub fn new(card: MaterialCard) -> Result<Self> {
        let density = PropertyKey::new(DENSITY_PROPERTY, Density::DIMS);
        let cp = PropertyKey::new(SPECIFIC_HEAT_CAPACITY_PROPERTY, SPECIFIC_HEAT_CAPACITY_DIMS);
        let caloric = ResolvedHeatCapacityCurve::try_new(&card, HeatCapacitySourceQuery {
            heat_capacity: &cp, density: &density, temperature_axis: "T",
            lower: &point(300.0)?, upper: &point(500.0)?, selection: SelectionPolicy::SingleClaimOnly,
        }, SolidLiquidPhase::Solid, 0.0, 1e-4)?;
        let dimensional = QuantitySpec::dimensional(THERMAL_CONDUCTIVITY_DIMS);
        let typed = QuantitySpec::semantic(SemanticType::new(QuantityKind::ThermalConductivity, ValueForm::Static));
        let quantity = card.claims().registered_quantity(THERMAL_CONDUCTIVITY_PROPERTY)
            .ok_or("material is missing thermal conductivity")?;
        if quantity != dimensional && quantity != typed {
            return Err("thermal conductivity needs the existing static conductivity quantity".into());
        }
        let k_requirement = ScalarPropertyRequirement::try_with_quantity(
            THERMAL_CONDUCTIVITY_PROPERTY, quantity, ScalarAdmissibility::StrictlyPositive)?;
        let state = resolve_material_state_point(&card, &point(300.0)?,
            std::slice::from_ref(&k_requirement), MaterialPropertySelection::SingleClaimOnly)?;
        let k_selection = MaterialPropertySelection::PinnedByProperty(vec![(
            THERMAL_CONDUCTIVITY_PROPERTY.into(), state.properties()[0].answer().receipt.selected,
        )]);
        Ok(Self { card, caloric, k_requirement, k_selection })
    }
    pub fn conductivity_at(&self, t: f64) -> Result<ResolvedMaterialStatePoint> {
        let state = resolve_material_state_point(&self.card, &point(t)?,
            std::slice::from_ref(&self.k_requirement), self.k_selection.clone())?;
        if !state.properties()[0].is_constant_at_fixed_temperature() {
            return Err("heat transport requires static coefficients or a T-only conductivity law".into());
        }
        if state.properties()[0].angular_frequency_band_rad_s()?.is_some_and(|(lo, _)| lo > 0.0) {
            return Err("conductivity source does not support static heat transport".into());
        }
        Ok(state)
    }
    pub fn load(path: &str) -> Result<Self> {
        const MAX_BYTES: u64 = 16*1024*1024;
        let mut bytes = Vec::new();
        std::fs::File::open(path)?.take(MAX_BYTES+1).read_to_end(&mut bytes)?;
        if bytes.len() as u64 > MAX_BYTES { return Err("material pack exceeds 16 MiB".into()); }
        let pack = NormalizedMaterialCardPack::from_bytes(&bytes)?;
        Self::new(pack.card().clone())
    }
}

/// An explicit synthetic dataset, NOT a metal or a measured specimen.
pub(super) fn demo(cp: f64, k: f64, elastic_upper_k: f64) -> Result<Material> {
    let mut claims = ClaimSet::new();
    for (name, dims, a, b, upper) in [
        (DENSITY_PROPERTY, Density::DIMS, 1000.0, 1000.0, 500.0),
        (SPECIFIC_HEAT_CAPACITY_PROPERTY, SPECIFIC_HEAT_CAPACITY_DIMS, cp, cp, 500.0),
        (THERMAL_CONDUCTIVITY_PROPERTY, THERMAL_CONDUCTIVITY_DIMS, k, k, 500.0),
        ("young_modulus", Pressure::DIMS, 1e6, 5e5, elastic_upper_k),
        ("shear_modulus", Pressure::DIMS, 4e5, 2e5, elastic_upper_k),
        (LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY, INVERSE_TEMPERATURE_DIMS, 1e-5, 3e-5, 500.0),
    ] {
        claims.insert_claim(PropertyClaim {
            key: PropertyKey::new(name, dims),
            value: PropertyValue::Curve { abscissa: "T".into(), abscissa_dims: Temperature::DIMS,
                knots: vec![(300.0, a), (upper, b)], dims },
            validity: ValidityDomain::unconstrained().with("T", 300.0, upper),
            uncertainty: UncertaintyModel::Unstated, interpolation: InterpolationPolicy::LinearInside,
            observations: vec![], provenance: Provenance { source: "synthetic thermoelastic rod data".into(),
                license: "CC0-1.0".into(), artifact: None },
        })?;
    }
    Material::new(MaterialCard::assemble(MaterialStateId { chemistry: "synthetic".into(),
        phase: "solid".into(), process: "illustrative data, not an identified specimen".into(), revision: 0 },
        claims, vec![])? )
}
