//! Original material state and physical geometry for a small-shear third body.
//!
//! Bulk simple shear gives F = (A G/h)e and a plastic shear-rate law with
//! resistance A eta_plastic/h above yield force A tau_y. The plastic viscosity
//! is an independently selected constitutive property. An oscillatory loss
//! viscosity, a friction coefficient or an instrument spectrum cannot supply it.
//! Time integration remains in fs-tribo's shared tangential-law interface.

use super::{MaterialStatePointError, ResolvedMaterialStatePoint, affine_temperature_slope,
    pinned_temperature_line};
use fs_matdb::MaterialCard;
use fs_qty::{Dims, DynViscosity, Pressure, QuantitySpec, Temperature};
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};

pub const SHEAR_LAYER_MODULUS_PROPERTY: &str = "shear_modulus";
pub const SHEAR_LAYER_YIELD_PROPERTY: &str = "shear_yield_stress";
pub const PLASTIC_SHEAR_VISCOSITY_PROPERTY: &str = "plastic_shear_viscosity";

fn invalid(quantity: &'static str) -> MaterialStatePointError {
    MaterialStatePointError::InvalidDerived { quantity }
}
fn positive(value: f64, quantity: &'static str) -> Result<f64, MaterialStatePointError> {
    if value.is_finite() && value > 0.0 { Ok(value) } else { Err(invalid(quantity)) }
}

/// Source-bound, pressure-independent isothermal elastic-viscoplastic shear.
/// Area and thickness describe a uniform, flat layer with prescribed retained
/// contact; this owner does not infer area from load or erase strain on opening.
/// The material and all query/claim receipts remain attached to the coefficients.
#[derive(Clone, Debug, PartialEq)]
pub struct ResolvedShearLayer {
    resolved: ResolvedMaterialStatePoint,
    area_m2: f64,
    thickness_m: f64,
    maximum_elastic_shear_strain: f64,
    maximum_elastic_slip_m: f64,
    stiffness_n_m: f64,
    yield_force_n: f64,
    plastic_resistance_n_s_m: f64,
    temperature_kelvin: f64,
    omega_band_rad_s: (f64, f64),
}
impl ResolvedShearLayer {
    /// Require positive G and eta_plastic, nonnegative shear yield stress, and
    /// explicit source support from DC. No pressure/strain/rate-dependent
    /// coefficient or applicability axis is frozen while contact loads change.
    /// Temperature-only curves may be evaluated at the explicit fixed T.
    ///
    /// The caller's maximum elastic engineering shear must lie in (0, 0.1].
    /// This is a declared small-strain model admission, not a measured error
    /// guarantee. Consumers must enforce h times that limit on stored slip and
    /// check their complete represented frequency band before propagation.
    pub fn new(resolved: ResolvedMaterialStatePoint, area_m2: f64, thickness_m: f64,
        maximum_elastic_shear_strain: f64) -> Result<Self, MaterialStatePointError> {
        positive(area_m2, "shear layer requires positive physical area")?;
        positive(thickness_m, "shear layer requires positive physical thickness")?;
        if !maximum_elastic_shear_strain.is_finite()
            || maximum_elastic_shear_strain <= 0.0 || maximum_elastic_shear_strain > 0.1 {
            return Err(invalid("shear layer requires an explicit small elastic shear limit in (0,0.1]"));
        }
        let mut temperature = None;
        for (axis, value) in resolved.query_point() {
            let quantity = resolved.axis_quantities().get(axis).copied();
            let absolute = quantity == Some(QuantitySpec::semantic(SemanticType::new(
                QuantityKind::AbsoluteTemperature, ValueForm::Static)));
            let legacy = axis == "T" && quantity.is_none_or(|q|
                q == QuantitySpec::dimensional(Temperature::DIMS));
            if absolute || legacy {
                positive(*value, "shear layer requires positive absolute temperature")?;
                if temperature.replace((axis.as_str(), *value)).is_some() {
                    return Err(invalid("shear layer needs exactly one absolute-temperature coordinate"));
                }
            }
        }
        let (temperature_axis, temperature_kelvin) = temperature
            .ok_or_else(|| invalid("shear layer requires an explicit absolute-temperature coordinate"))?;
        let mut band = (0.0_f64, f64::MAX);
        let mut coefficient = |name: &str, dims: Dims, allow_zero: bool|
            -> Result<f64, MaterialStatePointError> {
            let property = resolved.property(name).ok_or_else(||
                MaterialStatePointError::InvalidRequirement { property: name.into(),
                    reason: "missing unambiguous shear-layer constitutive property" })?;
            if property.requirement().quantity() != QuantitySpec::dimensional(dims) {
                return Err(invalid("shear-layer coefficient has incompatible quantity semantics"));
            }
            let value = property.value_si();
            if !value.is_finite() || value < 0.0 || (!allow_zero && value == 0.0) {
                return Err(invalid("shear-layer modulus/viscosity must be positive and yield nonnegative"));
            }
            if !property.is_constant_at_fixed_temperature() {
                return Err(invalid("shear-layer law needs constant or temperature-only coefficients"));
            }
            let validity = &property.answer().evidence.model.validity;
            if validity.bound(temperature_axis).is_none() {
                return Err(invalid("each shear-layer coefficient needs explicit absolute-temperature support"));
            }
            for axis in validity.bounds().keys() {
                if axis == temperature_axis { continue; }
                let quantity = resolved.axis_quantities().get(axis).copied();
                let frequency = quantity.and_then(QuantitySpec::semantic_type).is_some_and(|ty|
                    ty.form() == ValueForm::Static && matches!(ty.kind(), QuantityKind::Frequency(_)));
                let legacy = axis == "omega" && quantity.is_none_or(|q|
                    q == QuantitySpec::dimensional(fs_qty::Frequency::DIMS));
                if !frequency && !legacy {
                    return Err(invalid("isothermal shear layer cannot freeze pressure, strain, rate or other source coordinates"));
                }
            }
            let (lo, hi) = property.angular_frequency_band_rad_s()?
                .ok_or_else(|| invalid("each shear-layer coefficient needs explicit static/frequency applicability"))?;
            if lo != 0.0 || hi <= 0.0 {
                return Err(invalid("shear-layer coefficient must cover static storage and a positive frequency band"));
            }
            band = (band.0.max(lo), band.1.min(hi));
            Ok(value)
        };
        let modulus = coefficient(SHEAR_LAYER_MODULUS_PROPERTY, Pressure::DIMS, false)?;
        let yield_stress = coefficient(SHEAR_LAYER_YIELD_PROPERTY, Pressure::DIMS, true)?;
        let viscosity = coefficient(PLASTIC_SHEAR_VISCOSITY_PROPERTY, DynViscosity::DIMS, false)?;
        let geometry = positive(area_m2 / thickness_m, "shear-layer area/thickness is unrepresentable")?;
        let stiffness_n_m = positive(geometry * modulus, "shear-layer stiffness is unrepresentable")?;
        let plastic_resistance_n_s_m = positive(geometry * viscosity,
            "shear-layer plastic resistance is unrepresentable")?;
        let yield_force_n = area_m2 * yield_stress;
        if !yield_force_n.is_finite() || (yield_stress > 0.0 && yield_force_n == 0.0) {
            return Err(invalid("shear-layer yield force is unrepresentable"));
        }
        let maximum_elastic_slip_m = positive(thickness_m * maximum_elastic_shear_strain,
            "shear-layer small-strain slip limit is unrepresentable")?;
        Ok(Self { resolved, area_m2, thickness_m, maximum_elastic_shear_strain,
            maximum_elastic_slip_m, stiffness_n_m, yield_force_n,
            plastic_resistance_n_s_m, temperature_kelvin, omega_band_rad_s: band })
    }
    /// Complete original material identity, coefficients, state and usage receipts.
    #[must_use]
    pub const fn resolved(&self) -> &ResolvedMaterialStatePoint { &self.resolved }
    #[must_use]
    pub const fn area_m2(&self) -> f64 { self.area_m2 }
    #[must_use]
    pub const fn thickness_m(&self) -> f64 { self.thickness_m }
    #[must_use]
    pub const fn maximum_elastic_shear_strain(&self) -> f64 { self.maximum_elastic_shear_strain }
    #[must_use]
    pub const fn maximum_elastic_slip_m(&self) -> f64 { self.maximum_elastic_slip_m }
    #[must_use]
    pub const fn stiffness_n_m(&self) -> f64 { self.stiffness_n_m }
    #[must_use]
    pub const fn yield_force_n(&self) -> f64 { self.yield_force_n }
    #[must_use]
    pub const fn plastic_resistance_n_s_m(&self) -> f64 { self.plastic_resistance_n_s_m }
    #[must_use]
    pub const fn temperature_kelvin(&self) -> f64 { self.temperature_kelvin }

    /// Admit `G(T)` as one differentiable affine law over the complete thermal
    /// owner's temperature interval, retaining this exact modulus claim.
    ///
    /// The returned `dk/dT = (A/h) dG/dT` belongs to the free energy
    /// `psi = k(T)|e|^2/2`. Its entropy is `-k'|e|^2/2`; its internal energy is
    /// `psi + T s`, not `psi`. A coupled owner must carry both coefficient work
    /// and reversible heat when continuing stored strain at another T.
    /// Affine G has zero second T derivative, so it introduces no additional
    /// fixed-strain heat capacity. This does not admit expansion, phase change,
    /// nonlinear G(T), or derivative jumps hidden between two endpoints.
    pub fn thermal_stiffness_slope_n_m_k(
        &self,
        card: &MaterialCard,
        temperature_axis: &str,
        temperature_span_k: (f64, f64),
    ) -> Result<f64, MaterialStatePointError> {
        let (lo, hi) = temperature_span_k;
        let absolute = QuantitySpec::semantic(SemanticType::new(
            QuantityKind::AbsoluteTemperature, ValueForm::Static));
        let quantity = self.resolved.axis_quantities().get(temperature_axis).copied();
        let temperature = quantity == Some(absolute) || temperature_axis == "T"
            && quantity.is_none_or(|q| q == QuantitySpec::dimensional(Temperature::DIMS));
        if card.content_hash() != self.resolved.card_identity() || !temperature
            || !lo.is_finite() || !hi.is_finite() || lo <= 0.0 || hi <= lo
            || self.temperature_kelvin < lo || self.temperature_kelvin > hi
            || self.resolved.query_point().iter().find(|(name, _)| name == temperature_axis)
                .map(|(_, value)| *value)
                != Some(self.temperature_kelvin) {
            return Err(invalid("shear thermal path needs the original card, temperature axis and complete positive interval"));
        }
        let property = self.resolved.property(SHEAR_LAYER_MODULUS_PROPERTY)
            .expect("the immutable layer retains its modulus");
        let claim = card.claims().claim(property.answer().receipt.selected)
            .ok_or_else(|| invalid("shear thermal modulus claim is absent"))?;
        let lower = self.resolved.query().with(temperature_axis, lo)
            .map_err(|_| invalid("shear thermal path start violates its original schema"))?;
        let upper = self.resolved.query().with(temperature_axis, hi)
            .map_err(|_| invalid("shear thermal path end violates its original schema"))?;
        let knots = pinned_temperature_line(card, &claim.key, temperature_axis, &lower,
            &upper, property.answer().receipt.selected)?;
        if knots.iter().any(|(_, g)| !g.is_finite() || *g <= 0.0) {
            return Err(invalid("shear thermal modulus must be positive over the complete source path"));
        }
        let slope = affine_temperature_slope(&knots, self.temperature_kelvin,
            property.value_si()).map_err(invalid)?;
        let stiffness_slope = (self.area_m2 / self.thickness_m) * slope;
        if !stiffness_slope.is_finite() || slope != 0.0 && stiffness_slope == 0.0 {
            return Err(invalid("shear thermal stiffness derivative is unrepresentable"));
        }
        Ok(stiffness_slope)
    }
    #[must_use]
    pub const fn omega_band_rad_s(&self) -> (f64, f64) { self.omega_band_rad_s }
    pub fn require_angular_frequency_rad_s(&self, omega: f64) -> Result<(), MaterialStatePointError> {
        if !omega.is_finite() || omega < self.omega_band_rad_s.0 || omega > self.omega_band_rad_s.1 {
            return Err(invalid("shear-layer frequency exceeds original source applicability"));
        }
        Ok(())
    }
}


#[cfg(test)]
mod tests {
    use super::*;
    use crate::state_point::{MaterialPropertySelection, ScalarAdmissibility,
        ScalarPropertyRequirement, resolve_material_state_point};
    use fs_evidence::ValidityDomain;
    use fs_matdb::{ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId,
        PropertyClaim, PropertyKey, PropertyValue, Provenance, QueryPoint, UncertaintyModel};

    fn source(modulus: f64, frequency_lower: f64, pressure_bound: bool,
        legacy_viscosity: bool, temperature_curve: bool) -> ResolvedMaterialStatePoint {
        let mut claims = ClaimSet::new();
        let mut requirements = Vec::new();
        for (name, dims, value, admit) in [
            (SHEAR_LAYER_MODULUS_PROPERTY, Pressure::DIMS, modulus, ScalarAdmissibility::StrictlyPositive),
            (SHEAR_LAYER_YIELD_PROPERTY, Pressure::DIMS, 1000.0, ScalarAdmissibility::NonNegative),
            (if legacy_viscosity { "oscillatory_shear_viscosity" } else { PLASTIC_SHEAR_VISCOSITY_PROPERTY },
                DynViscosity::DIMS, 2000.0, ScalarAdmissibility::StrictlyPositive),
        ] {
            let mut validity = ValidityDomain::unconstrained()
                .with("T", 280.0, 320.0).with("omega", frequency_lower, 200_000.0);
            if pressure_bound { validity = validity.with("pressure", 0.0, 1e6); }
            let curve = temperature_curve && name == SHEAR_LAYER_MODULUS_PROPERTY;
            claims.insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value: if curve { PropertyValue::Curve {
                    abscissa: "T".into(), abscissa_dims: Temperature::DIMS,
                    knots: vec![(280.0, 0.5*value), (320.0, 1.5*value)], dims,
                } } else { PropertyValue::Scalar { value, dims } },
                validity, uncertainty: UncertaintyModel::Unstated,
                interpolation: if curve { InterpolationPolicy::LinearInside } else { InterpolationPolicy::ConstantWithinValidity },
                observations: vec![], provenance: Provenance {
                    source: "manufactured shear-layer constitutive fixture".into(),
                    license: "test-only".into(), artifact: None },
            }).unwrap();
            requirements.push(ScalarPropertyRequirement::try_new(name, dims, admit).unwrap());
        }
        let card = MaterialCard::assemble(MaterialStateId {
            chemistry: "manufactured".into(), phase: "solid".into(), process: "shear-test".into(), revision: 0,
        }, claims, vec![]).unwrap();
        let mut query = QueryPoint::new().with("T", 300.0).unwrap().with("omega", 1000.0).unwrap();
        if pressure_bound { query = query.with("pressure", 1e5).unwrap(); }
        resolve_material_state_point(&card, &query, &requirements,
            MaterialPropertySelection::SingleClaimOnly).unwrap()
    }

    #[test]
    fn g1_geometry_and_material_properties_drive_distinct_shear_storage_yield_and_flow() {
        let material = source(1e6, 0.0, false, false, false);
        let layer = ResolvedShearLayer::new(material.clone(), 1e-6, 1e-5, 0.05).unwrap();
        assert_eq!(layer.resolved(), &material);
        let area = ResolvedShearLayer::new(material.clone(), 2e-6, 1e-5, 0.05).unwrap();
        let thickness = ResolvedShearLayer::new(material, 1e-6, 2e-5, 0.05).unwrap();
        let stiffer = ResolvedShearLayer::new(source(2e6, 0.0, false, false, false),
            1e-6, 1e-5, 0.05).unwrap();
        assert_eq!(area.stiffness_n_m(), 2.0*layer.stiffness_n_m());
        assert_eq!(area.yield_force_n(), 2.0*layer.yield_force_n());
        assert_eq!(area.plastic_resistance_n_s_m(), 2.0*layer.plastic_resistance_n_s_m());
        assert_eq!(thickness.stiffness_n_m(), 0.5*layer.stiffness_n_m());
        assert_eq!(thickness.yield_force_n(), layer.yield_force_n());
        assert_eq!(thickness.maximum_elastic_slip_m(), 2.0*layer.maximum_elastic_slip_m());
        assert_eq!(stiffer.stiffness_n_m(), 2.0*layer.stiffness_n_m());
        assert_eq!(stiffer.yield_force_n(), layer.yield_force_n());
        assert_eq!(stiffer.plastic_resistance_n_s_m(), layer.plastic_resistance_n_s_m());
        assert_ne!(stiffer.resolved().identity(), layer.resolved().identity());
        // Volume integral of 1/2 G gamma_e^2 must equal 1/2 k e^2.
        let strain = 0.01_f64;
        let volume_energy = 0.5*1e6*strain*strain*layer.area_m2()*layer.thickness_m();
        let e = strain*layer.thickness_m();
        let spring_energy = 0.5*layer.stiffness_n_m()*e*e;
        assert!((spring_energy/volume_energy-1.0).abs() < 1e-14);
    }

    #[test]
    fn g0_static_source_semantics_refuse_pressure_or_oscillatory_viscosity_substitution() {
        for material in [
            source(1e6, 1.0, false, false, false),
            source(1e6, 0.0, true, false, false),
            source(1e6, 0.0, false, true, false),
        ] {
            assert!(ResolvedShearLayer::new(material, 1e-6, 1e-5, 0.05).is_err());
        }
        for (area, thickness, limit) in [(0.0, 1e-5, 0.05),
            (1e-6, 0.0, 0.05), (1e-6, 1e-5, 0.0), (1e-6, 1e-5, 0.11)] {
            assert!(ResolvedShearLayer::new(source(1e6, 0.0, false, false, false),
                area, thickness, limit).is_err());
        }
        let layer = ResolvedShearLayer::new(source(1e6, 0.0, false, false, false),
            1e-6, 1e-5, 0.05).unwrap();
        assert_eq!(layer.omega_band_rad_s(), (0.0, 200_000.0));
        assert!(layer.require_angular_frequency_rad_s(0.0).is_ok());
        assert!(layer.require_angular_frequency_rad_s(150_000.0).is_ok());
        assert!(layer.require_angular_frequency_rad_s(200_001.0).is_err());
        assert!(layer.require_angular_frequency_rad_s(f64::NAN).is_err());
    }

    #[test]
    fn g0_temperature_curve_retains_original_source_at_the_frozen_physical_state() {
        let layer = ResolvedShearLayer::new(source(1e6, 0.0, false, false, true),
            1e-6, 1e-5, 0.05).unwrap();
        assert_eq!(layer.temperature_kelvin(), 300.0);
        let property = layer.resolved().property(SHEAR_LAYER_MODULUS_PROPERTY).unwrap();
        assert_eq!(property.value_si(), 1e6);
        assert_eq!(property.source_axes(), &[("T".into(), Temperature::DIMS)]);
    }

    fn thermal_layer(knots: Vec<(f64, f64)>, axis: &str)
        -> (MaterialCard, ResolvedShearLayer) {
        let temperature = QuantitySpec::semantic(SemanticType::new(
            QuantityKind::AbsoluteTemperature, ValueForm::Static));
        let mut claims = ClaimSet::new();
        let mut requirements = Vec::new();
        for (name, dims, value) in [
            (SHEAR_LAYER_MODULUS_PROPERTY, Pressure::DIMS, 1e6),
            (SHEAR_LAYER_YIELD_PROPERTY, Pressure::DIMS, 1000.0),
            (PLASTIC_SHEAR_VISCOSITY_PROPERTY, DynViscosity::DIMS, 2000.0),
        ] {
            let modulus = name == SHEAR_LAYER_MODULUS_PROPERTY;
            claims.insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value: if modulus { PropertyValue::Curve {
                    abscissa: axis.into(), abscissa_dims: Temperature::DIMS,
                    knots: knots.clone(), dims,
                } } else { PropertyValue::Scalar { value, dims } },
                validity: ValidityDomain::unconstrained()
                    .with_quantity(axis, temperature, 280.0, 320.0)
                    .with("omega", 0.0, 200_000.0),
                uncertainty: UncertaintyModel::Unstated,
                interpolation: if modulus { InterpolationPolicy::LinearInside }
                    else { InterpolationPolicy::ConstantWithinValidity },
                observations: vec![], provenance: Provenance {
                    source: "manufactured affine shear thermal fixture".into(),
                    license: "test-only".into(), artifact: None },
            }).unwrap();
            requirements.push(ScalarPropertyRequirement::try_new(name, dims,
                ScalarAdmissibility::StrictlyPositive).unwrap());
        }
        let card = MaterialCard::assemble(MaterialStateId {
            chemistry: "manufactured".into(), phase: "solid".into(),
            process: "thermal-shear-test".into(), revision: 0,
        }, claims, vec![]).unwrap();
        let point = QueryPoint::new().with_quantity(axis, temperature, 300.0).unwrap()
            .with("omega", 1000.0).unwrap();
        let resolved = resolve_material_state_point(&card, &point, &requirements,
            MaterialPropertySelection::SingleClaimOnly).unwrap();
        let layer = ResolvedShearLayer::new(resolved, 2e-6, 1e-5, 0.05).unwrap();
        (card, layer)
    }

    #[test]
    fn g1_affine_shear_temperature_path_retains_typed_axis_and_exact_source_slope() {
        for axis in ["T", "film_temperature"] {
            let (card, layer) = thermal_layer(
                vec![(280.0, 1.2e6), (300.0, 1e6), (320.0, 0.8e6)], axis);
            let slope = layer.thermal_stiffness_slope_n_m_k(&card, axis,
                (280.0, 320.0)).unwrap();
            assert!((slope + 2000.0).abs() < 1e-12);
            assert!(layer.thermal_stiffness_slope_n_m_k(&card, axis,
                (279.0, 320.0)).is_err());
            assert!(layer.thermal_stiffness_slope_n_m_k(&card, "other_temperature",
                (280.0, 320.0)).is_err());
            let (constant_card, constant) = thermal_layer(
                vec![(280.0, 1e6), (320.0, 1e6)], axis);
            assert_eq!(constant.thermal_stiffness_slope_n_m_k(&constant_card, axis,
                (280.0, 320.0)).unwrap(), 0.0);
            assert!(layer.thermal_stiffness_slope_n_m_k(&constant_card, axis,
                (280.0, 320.0)).is_err());
        }
    }

    #[test]
    fn g0_thermal_shear_rejects_modulus_kinks_hidden_between_endpoint_values() {
        let (card, layer) = thermal_layer(
            vec![(280.0, 1.2e6), (300.0, 1.1e6), (320.0, 0.8e6)], "T");
        assert!(layer.thermal_stiffness_slope_n_m_k(&card, "T",
            (280.0, 320.0)).is_err());
    }
}
