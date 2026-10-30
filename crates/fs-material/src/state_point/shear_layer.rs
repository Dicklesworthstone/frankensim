//! Original material state and physical geometry for a small-shear third body.
//!
//! Bulk simple shear gives F = (A G/h)e and a plastic shear-rate law with
//! resistance A eta_plastic/h above yield force A tau_y. The plastic viscosity
//! is an independently selected constitutive property. An oscillatory loss
//! viscosity, a friction coefficient or an instrument spectrum cannot supply it.
//! Time integration remains in fs-tribo's shared tangential-law interface.

use super::{MaterialStatePointError, ResolvedMaterialStatePoint, affine_temperature_slope,
    pinned_temperature_line};
use fs_matdb::{MaterialCard, QueryPoint};
use fs_qty::{Dims, DynViscosity, Pressure, QuantitySpec, Temperature, Time};
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};

pub const SHEAR_LAYER_MODULUS_PROPERTY: &str = "shear_modulus";
pub const SHEAR_LAYER_YIELD_PROPERTY: &str = "shear_yield_stress";
pub const PLASTIC_SHEAR_VISCOSITY_PROPERTY: &str = "plastic_shear_viscosity";
/// Bounded source histories per retained material layer, independent of mesh.
pub const MAX_SHEAR_RELAXATION_ARMS: usize = 64;

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
    relaxation_properties: Vec<(String, String)>,
    relaxation: Vec<(f64, f64)>,
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
            plastic_resistance_n_s_m, temperature_kelvin, omega_band_rad_s: band,
            relaxation_properties: Vec::new(), relaxation: Vec::new() })
    }

    /// Bind Maxwell shear arms from the SAME original resolved material. Each
    /// pair names (DeltaG [Pa], tau [s]); `shear_modulus` remains G_infinity.
    /// Every arm shares the original physical area/thickness and temperature.
    /// Source applicability must cover DC and the consumer's actual band;
    /// no oscillatory modulus fit, missing coefficient or relaxation is inferred.
    pub fn with_relaxation(mut self, properties: Vec<(String, String)>)
        -> Result<Self, MaterialStatePointError> {
        if properties.len() > MAX_SHEAR_RELAXATION_ARMS || !self.relaxation_properties.is_empty() {
            return Err(invalid("shear relaxation needs a bounded one-time source binding"));
        }
        let temperature_axis = self.resolved.query_point().iter().find(|(axis, _)| {
            self.resolved.axis_quantities().get(axis).copied().is_some_and(|quantity|
                quantity.semantic_type().is_some_and(|ty|
                    ty.kind() == QuantityKind::AbsoluteTemperature && ty.form() == ValueForm::Static))
                || axis == "T"
        }).map(|(axis, _)| axis.as_str()).expect("admitted unique temperature");
        let mut band = self.omega_band_rad_s;
        let mut relaxation = Vec::with_capacity(properties.len());
        for (index, (modulus, time)) in properties.iter().enumerate() {
            if modulus == SHEAR_LAYER_MODULUS_PROPERTY
                || properties[..index].contains(&(modulus.clone(), time.clone())) {
                return Err(invalid("shear relaxation cannot repeat an arm or count equilibrium modulus twice"));
            }
            let mut values = [0.0; 2];
            for (slot, (name, dims)) in [(modulus, Pressure::DIMS), (time, Time::DIMS)].iter().enumerate() {
                let property = self.resolved.property(name).ok_or_else(||
                    MaterialStatePointError::InvalidRequirement { property: (*name).clone(),
                        reason: "missing selected shear relaxation source property" })?;
                if property.requirement().quantity() != QuantitySpec::dimensional(*dims)
                    || !property.is_constant_at_fixed_temperature() {
                    return Err(invalid("shear memory needs typed, fixed-state or temperature-only coefficients"));
                }
                values[slot] = positive(property.value_si(), "shear memory modulus and time must be positive")?;
                let validity = &property.answer().evidence.model.validity;
                if validity.bound(temperature_axis).is_none() {
                    return Err(invalid("shear memory requires explicit temperature applicability"));
                }
                for axis in validity.bounds().keys() {
                    if axis == temperature_axis { continue; }
                    let quantity = self.resolved.axis_quantities().get(axis).copied();
                    let frequency = quantity.and_then(QuantitySpec::semantic_type).is_some_and(|ty|
                        ty.form() == ValueForm::Static && matches!(ty.kind(), QuantityKind::Frequency(_)));
                    let legacy = axis == "omega" && quantity.is_none_or(|q|
                        q == QuantitySpec::dimensional(fs_qty::Frequency::DIMS));
                    if !frequency && !legacy {
                        return Err(invalid("shear memory cannot freeze another constitutive state coordinate"));
                    }
                }
                let (lo, hi) = property.angular_frequency_band_rad_s()?
                    .ok_or_else(|| invalid("shear memory requires a declared frequency band"))?;
                if lo != 0.0 || hi <= 0.0 {
                    return Err(invalid("shear memory must include static storage and positive bandwidth"));
                }
                band = (band.0.max(lo), band.1.min(hi));
            }
            let stiffness = positive(self.area_m2/self.thickness_m*values[0],
                "shear memory stiffness is unrepresentable")?;
            relaxation.push((stiffness, values[1]));
        }
        self.relaxation = relaxation;
        self.relaxation_properties = properties;
        self.omega_band_rad_s = band;
        Ok(self)
    }
    /// Original source roles in stable constitutive-history order.
    pub fn relaxation_properties(&self) -> &[(String, String)] { &self.relaxation_properties }
    /// Physical (A*DeltaG/h [N/m], tau [s]) for each original source arm.
    pub fn relaxation(&self) -> &[(f64, f64)] { &self.relaxation }
    /// Re-evaluate the original complete source without losing a memory arm,
    /// changing geometry or substituting a property selection policy.
    pub fn at_query_point(&self, card: &MaterialCard, query: &QueryPoint)
        -> Result<Self, MaterialStatePointError> {
        Self::new(self.resolved.at_query_point(card, query)?, self.area_m2,
            self.thickness_m, self.maximum_elastic_shear_strain)?
            .with_relaxation(self.relaxation_properties.clone())
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
        // The shared film thermal owner accounts affine equilibrium entropy.
        // Maxwell arms store their own strain; an unaccounted d(DeltaG)/dT
        // would omit their entropy. Permit temperature-dependent tau, whose
        // change alters dissipation without changing the elastic potential.
        for (name, _) in &self.relaxation_properties {
            let arm = self.resolved.property(name).expect("bound shear memory property");
            let claim = card.claims().claim(arm.answer().receipt.selected)
                .ok_or_else(|| invalid("shear memory thermal source claim is absent"))?;
            let values = pinned_temperature_line(card, &claim.key, temperature_axis,
                &lower, &upper, arm.answer().receipt.selected)?;
            if values.iter().any(|(_, value)| *value != arm.value_si()) {
                return Err(invalid("thermal shear memory requires temperature-independent arm moduli"));
            }
        }
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
            ("memory_modulus", Pressure::DIMS, 5e5),
            ("memory_time", Time::DIMS, 0.01),
            ("thermal_memory_modulus", Pressure::DIMS, 1e6),
        ] {
            let modulus = name == SHEAR_LAYER_MODULUS_PROPERTY || name == "thermal_memory_modulus";
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
    fn g1_shear_memory_retains_original_material_geometry_and_thermal_roles() {
        let (card, original) = thermal_layer(vec![(280.0, 1.2e6), (320.0, 0.8e6)], "T");
        let names = vec![("memory_modulus".into(), "memory_time".into())];
        let layer = original.clone().with_relaxation(names.clone()).unwrap();
        assert_eq!(layer.resolved(), original.resolved());
        assert_eq!(layer.relaxation_properties(), names);
        assert!((layer.relaxation()[0].0/1e5-1.0).abs() < 1e-14);
        assert_eq!(layer.relaxation()[0].1, 0.01);
        assert_eq!(layer.stiffness_n_m(), original.stiffness_n_m());
        let warmer = layer.at_query_point(&card,
            &layer.resolved().query().with("T", 310.0).unwrap()).unwrap();
        assert_eq!(warmer.temperature_kelvin(), 310.0);
        assert_eq!(warmer.relaxation_properties(), layer.relaxation_properties());
        assert_eq!(warmer.relaxation(), layer.relaxation());
        assert!(warmer.stiffness_n_m() < layer.stiffness_n_m());
        assert!(layer.thermal_stiffness_slope_n_m_k(&card, "T", (280.0, 320.0)).is_ok());
        let half = ResolvedShearLayer::new(original.resolved().clone(),
            original.area_m2()/2.0, original.thickness_m(), 0.05).unwrap()
            .with_relaxation(names).unwrap();
        assert_eq!(half.relaxation()[0], (0.5*layer.relaxation()[0].0, 0.01));
        for names in [
            vec![("missing".into(), "memory_time".into())],
            vec![("memory_modulus".into(), "memory_modulus".into())],
            vec![(SHEAR_LAYER_MODULUS_PROPERTY.into(), "memory_time".into())],
            vec![("memory_modulus".into(), "memory_time".into()); 2],
        ] { assert!(original.clone().with_relaxation(names).is_err()); }
        let variable = original.with_relaxation(
            vec![("thermal_memory_modulus".into(), "memory_time".into())]).unwrap();
        assert!(variable.thermal_stiffness_slope_n_m_k(&card, "T", (280.0, 320.0)).is_err(),
            "a changing arm modulus cannot omit its mechanical entropy");
    }

    #[test]
    fn g0_thermal_shear_rejects_modulus_kinks_hidden_between_endpoint_values() {
        let (card, layer) = thermal_layer(
            vec![(280.0, 1.2e6), (300.0, 1.1e6), (320.0, 0.8e6)], "T");
        assert!(layer.thermal_stiffness_slope_n_m_k(&card, "T",
            (280.0, 320.0)).is_err());
    }
}
