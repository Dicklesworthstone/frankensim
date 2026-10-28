//! Intrinsic pore permeability and quasi-steady, isothermal gas permeation.
//!
//! Darcy superficial velocity is u = -K grad(p)/mu. Ideal-gas conservation
//! makes rho*u constant through a slab, not volume flow. Integrating
//! k_app = k*(1+b/p) gives
//! j*mu*R*T*sum(L/k) = (p_left-p_right)*(p_mean+b).
//! This is exact for constant temperature/composition and common-b layers.
//! Factoring the pressure difference avoids cancellation of pressure squares.
//! With different layer factors, each layer instead obeys its own integrated
//! potential. A bounded scalar flux solve propagates from low to high pressure;
//! the factors are never averaged into an apparent common material.
//! Formulation references:
//! <https://doc.comsol.com/6.3/doc/com.comsol.help.echem/echem_ug_fluidflow.08.111.html>
//! <https://documentation.pflotran.org/user_guide/cards/subsurface/klinkenberg_effect_card.html>
//!
//! Intrinsic permeability [m2], hydraulic conductivity, membrane permeability
//! and permeance are distinct quantities. This is a rigid, dry, gas-filled
//! skeleton, not a selective membrane, Darcy-Forchheimer inertia, pore-storage
//! model or Knudsen flow. Explicit source support and speed limits bound the
//! linear/slip regime. No universal gas slip coefficient is supplied.

use crate::{MaterialError, gas::GasState, state_point::ResolvedMaterialStatePoint,
    thermochemical_gas::ThermochemicalGasModel};
use fs_qty::{Dims, QuantitySpec};

pub const INTRINSIC_PERMEABILITY_DIMS: Dims = Dims([2, 0, 0, 0, 0, 0]);
/// Maximum compiled laminate size. Common-factor queries take constant work.
pub const MAX_POROUS_LAYERS: usize = 64;
/// Maximum scalar iterations for a heterogeneous-factor laminate.
pub const MAX_LAYERED_DARCY_ITERATIONS: usize = 128;

fn invalid(what: &str) -> MaterialError {
    MaterialError::Parameters { what: what.into() }
}
fn positive(value: f64) -> bool { value.is_finite() && value > 0.0 }
fn nonnegative(value: f64) -> bool { value.is_finite() && value >= 0.0 }
fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0]*b[0] + a[1]*b[1] + a[2]*b[2]
}
fn unit(direction: [f64; 3]) -> bool {
    direction.iter().all(|x| x.is_finite()) && (dot(direction, direction)-1.0).abs() <= 2e-12
}

/// Nonnegative principal intrinsic permeabilities with explicit world axes.
/// Zero is sealed along that material axis; it is never replaced by a floor.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct IntrinsicPermeability {
    principal_m2: [f64; 3],
    axes: [[f64; 3]; 3],
}
impl IntrinsicPermeability {
    pub fn isotropic(value_m2: f64) -> Result<Self, MaterialError> {
        Self::orthotropic([value_m2; 3], [[1.,0.,0.], [0.,1.,0.], [0.,0.,1.]])
    }
    /// Rows are mutually orthogonal unit material axes in world coordinates,
    /// in exactly the order of the three principal permeabilities.
    pub fn orthotropic(principal_m2: [f64; 3], axes: [[f64; 3]; 3]) -> Result<Self, MaterialError> {
        if !principal_m2.iter().all(|&k| nonnegative(k)) || !axes.iter().all(|&a| unit(a))
            || dot(axes[0], axes[1]).abs() > 2e-12
            || dot(axes[0], axes[2]).abs() > 2e-12
            || dot(axes[1], axes[2]).abs() > 2e-12 {
            return Err(invalid("intrinsic permeability needs nonnegative m2 and explicit orthonormal axes"));
        }
        Ok(Self { principal_m2, axes })
    }
    pub const fn principal_m2(&self) -> &[f64; 3] { &self.principal_m2 }
    pub const fn material_axes_in_world(&self) -> &[[f64; 3]; 3] { &self.axes }
    /// n.K.n for a planar normal pressure gradient in a wide layer; tangential
    /// flow is unconstrained. A laterally confined anisotropic plug instead
    /// requires 1/(n.K^-1.n), a different physical boundary condition.
    pub fn normal_permeability(&self, normal: [f64; 3]) -> Result<f64, MaterialError> {
        if !unit(normal) { return Err(invalid("permeability normal must be a finite unit vector")); }
        let mut value = 0.0;
        let mut conducting = false;
        for (k, axis) in self.principal_m2.iter().zip(self.axes) {
            let projection = dot(normal, axis);
            conducting |= *k > 0.0 && projection != 0.0;
            value += k * projection * projection;
        }
        if !nonnegative(value) || conducting && value == 0.0 {
            return Err(invalid("normal permeability is unrepresentable"));
        }
        Ok(value)
    }
}

/// Explicit source/model domain. Speed is superficial (Darcy) velocity, not
/// pore velocity, and must follow the supplied linear-regime evidence.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PorousGasSupport {
    pub temperature_k: (f64, f64),
    pub pressure_pa: (f64, f64),
    pub maximum_superficial_speed_m_s: f64,
}
impl PorousGasSupport {
    fn validate(self) -> Result<(), MaterialError> {
        if ![self.temperature_k.0, self.temperature_k.1, self.pressure_pa.0,
            self.pressure_pa.1, self.maximum_superficial_speed_m_s].into_iter().all(positive)
            || self.temperature_k.0 > self.temperature_k.1 || self.pressure_pa.0 > self.pressure_pa.1 {
            return Err(invalid("porous gas support needs positive ordered T/p bounds and a finite positive speed limit"));
        }
        Ok(())
    }
    fn admit(self, gas: &GasState, left: f64, right: f64) -> Result<(), MaterialError> {
        gas.validate_coherence()?;
        if gas.temperature < self.temperature_k.0 || gas.temperature > self.temperature_k.1
            || [left, right, gas.pressure].iter().any(|&p| !p.is_finite()
                || p < self.pressure_pa.0 || p > self.pressure_pa.1) {
            return Err(invalid("gas permeability query lies outside source temperature/absolute-pressure support"));
        }
        Ok(())
    }
}

/// One unchanged skeleton layer. A material bundle retains exact source claims,
/// units, query coordinates and uncertainty; a citation/model-input coefficient
/// is never promoted to a measured material.
#[derive(Clone, Debug, PartialEq)]
pub struct PorousGasLayer {
    permeability: IntrinsicPermeability,
    thickness_m: f64,
    klinkenberg_slip_pa: f64,
    provenance: String,
    source: Option<ResolvedMaterialStatePoint>,
    wet_state: Option<super::wet::WetPoreGasState>,
}
impl PorousGasLayer {
    /// Explicit data/model-input constructor. The caller separately supplies
    /// the physical applicability domain; a citation is not qualification.
    pub fn try_new(permeability: IntrinsicPermeability, thickness_m: f64,
        klinkenberg_slip_pa: f64, provenance: impl Into<String>) -> Result<Self, MaterialError> {
        let provenance = provenance.into();
        if !positive(thickness_m) || !nonnegative(klinkenberg_slip_pa)
            || provenance.trim().is_empty() || provenance.len() > 4096 {
            return Err(invalid("porous layer needs positive thickness, nonnegative slip and bounded explicit provenance"));
        }
        Ok(Self { permeability, thickness_m, klinkenberg_slip_pa, provenance,
            source: None, wet_state: None })
    }
    /// Bind source-resolved intrinsic permeability [m2], and optionally an
    /// explicitly selected Klinkenberg factor [Pa]. Only validity-wide scalar
    /// constants are admitted: a sampled apparent-permeability curve is not k0.
    pub fn from_material_isotropic(source: ResolvedMaterialStatePoint,
        permeability_property: &str, thickness_m: f64, slip_property: Option<&str>)
        -> Result<Self, MaterialError> {
        let k = source_value(&source, permeability_property, INTRINSIC_PERMEABILITY_DIMS)?;
        Self::bind_source(source, IntrinsicPermeability::isotropic(k)?, thickness_m, slip_property)
    }
    pub fn from_material_orthotropic(source: ResolvedMaterialStatePoint,
        properties: [&str; 3], axes: [[f64; 3]; 3], thickness_m: f64,
        slip_property: Option<&str>) -> Result<Self, MaterialError> {
        if properties[0] == properties[1] || properties[0] == properties[2] || properties[1] == properties[2] {
            return Err(invalid("orthotropic permeability needs three distinct source properties"));
        }
        let mut k = [0.0; 3];
        for i in 0..3 { k[i] = source_value(&source, properties[i], INTRINSIC_PERMEABILITY_DIMS)?; }
        Self::bind_source(source, IntrinsicPermeability::orthotropic(k, axes)?, thickness_m, slip_property)
    }
    fn bind_source(source: ResolvedMaterialStatePoint, permeability: IntrinsicPermeability,
        thickness_m: f64, slip_property: Option<&str>) -> Result<Self, MaterialError> {
        let slip = slip_property.map(|name| source_value(&source, name, fs_qty::Pressure::DIMS))
            .transpose()?.unwrap_or(0.0);
        let mut layer = Self::try_new(permeability, thickness_m, slip,
            format!("source-bound permeability material {:?}", source.card_identity()))?;
        layer.source = Some(source);
        Ok(layer)
    }
    pub const fn permeability(&self) -> &IntrinsicPermeability { &self.permeability }
    pub const fn thickness_m(&self) -> f64 { self.thickness_m }
    pub const fn klinkenberg_slip_pa(&self) -> f64 { self.klinkenberg_slip_pa }
    pub fn provenance(&self) -> &str { &self.provenance }
    pub const fn source(&self) -> Option<&ResolvedMaterialStatePoint> { self.source.as_ref() }
    /// Current prescribed saturation/geometry response; intrinsic k and its
    /// original material source remain unchanged. No liquid inventory advances.
    pub const fn wet_state(&self) -> Option<&super::wet::WetPoreGasState> { self.wet_state.as_ref() }
    pub fn gas_relative_permeability(&self) -> f64 {
        self.wet_state.as_ref().map_or(1.0, |state| state.gas_relative_permeability())
    }
    /// Replace, rather than multiply, the previous wet response. A dry-gas
    /// slip calibration does not qualify wet pores and therefore refuses.
    pub fn with_wet_state(&self, state: super::wet::WetPoreGasState) -> Result<Self, MaterialError> {
        if self.klinkenberg_slip_pa != 0.0 {
            return Err(invalid("wet-pore mobility requires no-slip intrinsic permeability; dry slip calibration does not qualify wet pores"));
        }
        let mut layer = self.clone();
        layer.wet_state = Some(state);
        Ok(layer)
    }
}
fn source_value(source: &ResolvedMaterialStatePoint, name: &str, dims: Dims) -> Result<f64, MaterialError> {
    let property = source.property(name).ok_or_else(|| invalid("missing selected permeability/slip source property"))?;
    if property.requirement().quantity() != QuantitySpec::dimensional(dims)
        || !property.source_axes().is_empty()
        || property.answer().receipt.decision != fs_matdb::EvaluationDecision::ConstantWithinValidity
        || !nonnegative(property.value_si()) {
        return Err(invalid("permeability/slip source must be a nonnegative validity-wide constant in its exact quantity"));
    }
    let answer = property.answer();
    for axis in answer.evidence.model.validity.bounds().keys() {
        let quantity = answer.receipt.axis_quantities.get(axis).copied();
        let value = source.query_point().iter().find(|(name, _)| name == axis).map(|(_, value)| *value);
        if source_axis(axis, quantity)?.is_some() {
            if !value.is_some_and(positive) {
                return Err(invalid("permeability source needs positive absolute T/p coordinates"));
            }
        } else if !crate::state_point::is_fixed_external_coordinate(axis, quantity, value) {
            return Err(invalid("permeability source has an unsupported evolving condition"));
        }
    }
    Ok(property.value_si())
}

#[derive(Clone, Copy)]
enum SourceAxis { Temperature, Pressure }
// Only named total-state coordinates belong to this constitutive owner.
// A pressure-dimensional specimen stress, partial pressure or reference value
// cannot silently become the actual absolute pore-gas pressure.
fn source_axis(axis: &str, quantity: Option<QuantitySpec>) -> Result<Option<SourceAxis>, MaterialError> {
    use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
    let thermal = quantity.is_none_or(|q| q == QuantitySpec::dimensional(fs_qty::Temperature::DIMS)
        || q == QuantitySpec::semantic(SemanticType::new(QuantityKind::AbsoluteTemperature, ValueForm::Static)));
    let pressure = quantity.is_none_or(|q| q == QuantitySpec::dimensional(fs_qty::Pressure::DIMS)
        || q == QuantitySpec::semantic(SemanticType::new(QuantityKind::Pressure, ValueForm::Static)));
    if axis == "T" || axis == "temperature" {
        return if thermal { Ok(Some(SourceAxis::Temperature)) }
            else { Err(invalid("permeability temperature axis has incompatible quantity semantics")) };
    }
    if axis == "pressure" {
        return if pressure { Ok(Some(SourceAxis::Pressure)) }
            else { Err(invalid("permeability pressure axis has incompatible quantity semantics")) };
    }
    if quantity.is_some_and(|q| q.dims() == fs_qty::Temperature::DIMS || q.dims() == fs_qty::Pressure::DIMS) {
        return Err(invalid("permeability source has an ambiguous temperature/pressure alias"));
    }
    Ok(None)
}

/// Precompiled series geometry and calibrated per-layer slip factors.
/// Common-factor flow is analytic; heterogeneous factors use a bounded,
/// allocation-free interface-pressure solve in the original layer order.
#[derive(Clone, Debug, PartialEq)]
pub struct DarcyGasBarrier {
    layers: Vec<PorousGasLayer>,
    normal: [f64; 3],
    support: PorousGasSupport,
    /// None is exactly sealed, not an artificially enormous finite resistance.
    geometric_resistance_per_m: Option<f64>,
    resistance_fractions: [f64; MAX_POROUS_LAYERS],
    common_slip_pa: Option<f64>,
    minimum_slip_pa: f64,
    maximum_slip_pa: f64,
    slip_gas_model: Option<ThermochemicalGasModel>,
}

/// Signed left-to-right flow. Each face uses its own density: mass is conserved
/// even when compression makes the two volume flow rates different.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DarcyGasFlow {
    pub mass_rate_left_to_right_kg_s: f64,
    pub mass_flux_left_to_right_kg_m2_s: f64,
    pub left_volume_rate_m3_s: f64,
    pub right_volume_rate_m3_s: f64,
    /// Mean-density volume flow divided by pressure difference; at equilibrium
    /// this is the infinitesimal conductance, including every layer's slip.
    pub mean_volume_conductance_m3_pa_s: f64,
    /// A*integral(u*(-dp)), excluding thermal storage. Diagnostic, not extra heat.
    pub drag_dissipation_w: f64,
    pub temperature_k: f64,
    pub left_pressure_pa: f64,
    pub right_pressure_pa: f64,
}

/// Optional resolved laminate profile; normal flow queries do not allocate it.
#[derive(Clone, Debug, PartialEq)]
pub struct DarcyGasFlowProfile {
    pub flow: DarcyGasFlow,
    /// Absolute pressures at all n+1 layer faces, in original left-to-right
    /// order, including the prescribed endpoints. None for a sealed path or
    /// zero exposed area: those inputs do not determine the interior pressure.
    /// Very small layer drops can round to equal absolute face pressures.
    pub interface_pressures_pa: Option<Vec<f64>>,
    /// Scalar iterations; zero for the exact common-factor/equilibrium paths.
    pub iterations: usize,
    /// Reconstructed total pressure rise minus the prescribed rise, divided
    /// by that rise, evaluated in scaled increments before absolute-pressure
    /// rounding. Zero at equilibrium. This is a floating-point residual,
    /// not an interval enclosure or independent experimental validation.
    pub normalized_pressure_residual: f64,
}

#[derive(Clone, Copy)]
struct LayeredRoot { scaled_flux: f64, iterations: usize }

#[derive(Clone, Copy)]
struct LayeredQuery {
    low: f64,
    drop: f64,
    scale: f64,
    low_scaled: f64,
    drop_scaled: f64,
    reverse_layers: bool,
}

impl LayeredQuery {
    fn new(left: f64, right: f64, maximum_slip: f64) -> Result<Self, MaterialError> {
        let low = left.min(right);
        let drop = (left-right).abs();
        let scale = left.max(right).max(maximum_slip);
        let low_scaled = low/scale;
        let drop_scaled = drop/scale;
        if ![low_scaled, drop_scaled].into_iter().all(positive) {
            return Err(invalid("laminate pressure scaling is unrepresentable"));
        }
        Ok(Self { low, drop, scale, low_scaled, drop_scaled, reverse_layers: left > right })
    }
}

#[derive(Clone, Copy)]
struct PressureMarch { fraction: f64, derivative: f64 }

impl DarcyGasBarrier {
    /// Nonzero slip requires the gas-model constructor because GasState has
    /// no chemical identity and cannot qualify gas-specific slip on its own.
    pub fn try_new(layers: Vec<PorousGasLayer>, normal: [f64; 3], support: PorousGasSupport)
        -> Result<Self, MaterialError> {
        Self::compile(layers, normal, support, None)
    }
    /// Bind slip to the complete immutable gas model, including its canonical
    /// composition and source laws. No molecular-mass or name inference.
    pub fn try_new_for_gas_model(layers: Vec<PorousGasLayer>, normal: [f64; 3],
        support: PorousGasSupport, model: &ThermochemicalGasModel) -> Result<Self, MaterialError> {
        Self::compile(layers, normal, support, Some(model.clone()))
    }
    fn compile(layers: Vec<PorousGasLayer>, normal: [f64; 3], mut support: PorousGasSupport,
        model: Option<ThermochemicalGasModel>) -> Result<Self, MaterialError> {
        support.validate()?;
        if layers.is_empty() || layers.len() > MAX_POROUS_LAYERS || !unit(normal) {
            return Err(invalid("Darcy barrier needs 1..64 layers and an explicit unit normal"));
        }
        let first_slip = layers[0].klinkenberg_slip_pa;
        let common_slip = layers.iter().all(|layer| layer.klinkenberg_slip_pa == first_slip)
            .then_some(first_slip);
        let minimum_slip = layers.iter().map(|layer| layer.klinkenberg_slip_pa).fold(f64::INFINITY, f64::min);
        let maximum_slip = layers.iter().map(|layer| layer.klinkenberg_slip_pa).fold(0.0, f64::max);
        if maximum_slip > 0.0 && model.is_none() {
            return Err(invalid("nonzero Klinkenberg slip requires its explicit calibrated gas model"));
        }
        let mut permeabilities = [0.0; MAX_POROUS_LAYERS];
        let mut sealed = false;
        for (i, layer) in layers.iter().enumerate() {
            let intrinsic = layer.permeability.normal_permeability(normal)?;
            let relative = layer.gas_relative_permeability();
            permeabilities[i] = intrinsic*relative;
            if intrinsic > 0.0 && relative > 0.0 && permeabilities[i] == 0.0 {
                return Err(invalid("nonzero wet-pore Darcy permeability is unrepresentable"));
            }
            sealed |= permeabilities[i] == 0.0;
            if let Some(source) = &layer.source {
                for property in source.properties() {
                    let answer = property.answer();
                    for (axis, &(lo, hi)) in answer.evidence.model.validity.bounds() {
                        let quantity = answer.receipt.axis_quantities.get(axis).copied();
                        match source_axis(axis, quantity)? {
                            Some(SourceAxis::Temperature) => support.temperature_k =
                                (support.temperature_k.0.max(lo), support.temperature_k.1.min(hi)),
                            Some(SourceAxis::Pressure) => support.pressure_pa =
                                (support.pressure_pa.0.max(lo), support.pressure_pa.1.min(hi)),
                            None => {}
                        }
                    }
                }
            }
        }
        let mut resistance = 0.0;
        let mut resistance_fractions = [0.0; MAX_POROUS_LAYERS];
        if !sealed {
            for (i, (layer, k)) in layers.iter().zip(permeabilities).enumerate() {
                let part = layer.thickness_m/k;
                if !positive(part) || !positive(resistance+part) {
                    return Err(invalid("porous layer series resistance is unrepresentable"));
                }
                resistance_fractions[i] = part;
                resistance += part;
            }
            for fraction in &mut resistance_fractions[..layers.len()] {
                *fraction /= resistance;
                if common_slip.is_none() && !positive(*fraction) {
                    return Err(invalid("heterogeneous laminate resistance fractions are unrepresentable"));
                }
            }
        }
        support.validate()?;
        Ok(Self { layers, normal, support,
            geometric_resistance_per_m: if sealed { None } else { Some(resistance) },
            resistance_fractions, common_slip_pa: common_slip,
            minimum_slip_pa: minimum_slip, maximum_slip_pa: maximum_slip, slip_gas_model: model })
    }
    pub fn layers(&self) -> &[PorousGasLayer] { &self.layers }
    pub const fn normal(&self) -> [f64; 3] { self.normal }
    pub const fn support(&self) -> PorousGasSupport { self.support }
    pub const fn geometric_resistance_per_m(&self) -> Option<f64> { self.geometric_resistance_per_m }
    /// The factor only when every layer has exactly the same calibrated b.
    /// None means heterogeneous factors, never an averaged or inferred factor.
    pub const fn klinkenberg_slip_pa(&self) -> Option<f64> { self.common_slip_pa }
    pub const fn calibrated_gas_model(&self) -> Option<&ThermochemicalGasModel> { self.slip_gas_model.as_ref() }
    /// Compile one current saturation response per original layer. The original
    /// ordered geometry, intrinsic tensors, sources and gas support are retained.
    /// Updating a wet barrier never applies its previous relative factor again.
    pub fn with_wet_states(&self, states: &[super::wet::WetPoreGasState]) -> Result<Self, MaterialError> {
        if states.len() != self.layers.len() {
            return Err(invalid("wet Darcy response must cover every original layer exactly once"));
        }
        let layers = self.layers.iter().zip(states)
            .map(|(layer, state)| layer.with_wet_state(state.clone()))
            .collect::<Result<Vec<_>, _>>()?;
        Self::compile(layers, self.normal, self.support, self.slip_gas_model.clone())
    }

    /// Frozen, low-pressure-drop resistance [Pa s/m3]. None means sealed.
    /// The consumer owns its low-drop bound and must enforce the solved
    /// superficial-speed support. Slip refuses, rather than being frozen.
    pub fn linear_volume_resistance(&self, gas: &GasState, area_m2: f64)
        -> Result<Option<f64>, MaterialError> {
        self.support.admit(gas, gas.pressure, gas.pressure)?;
        if !nonnegative(area_m2) || self.maximum_slip_pa != 0.0 {
            return Err(invalid("linear Darcy resistance needs nonnegative area and continuum no-slip permeability"));
        }
        let Some(geometry) = self.geometric_resistance_per_m.filter(|_| area_m2 > 0.0) else { return Ok(None); };
        let resistance = gas.dynamic_viscosity * geometry / area_m2;
        if !positive(resistance) { return Err(invalid("Darcy volume resistance is unrepresentable")); }
        Ok(Some(resistance))
    }

    pub fn mass_flow(&self, gas: &GasState, left_pressure_pa: f64,
        right_pressure_pa: f64, area_m2: f64) -> Result<DarcyGasFlow, MaterialError> {
        if self.maximum_slip_pa > 0.0 {
            return Err(invalid("gas-specific Klinkenberg flow requires mass_flow_for_model"));
        }
        self.evaluate(gas, left_pressure_pa, right_pressure_pa, area_m2).map(|(flow, _)| flow)
    }
    /// State and model must be supplied from the same gas owner. That owner
    /// admits gas T/p/source/phase; this adds the permeability-domain checks.
    pub fn mass_flow_for_model(&self, gas: &GasState, model: &ThermochemicalGasModel,
        left_pressure_pa: f64, right_pressure_pa: f64, area_m2: f64)
        -> Result<DarcyGasFlow, MaterialError> {
        if self.maximum_slip_pa > 0.0 && self.slip_gas_model.as_ref() != Some(model) {
            return Err(invalid("Klinkenberg calibration does not cover this gas model/composition"));
        }
        self.evaluate(gas, left_pressure_pa, right_pressure_pa, area_m2).map(|(flow, _)| flow)
    }

    /// Flow plus interface pressures. This diagnostic path allocates at most
    /// n+1 pressures; ordinary mass-flow queries remain allocation-free.
    pub fn mass_flow_profile(&self, gas: &GasState, left: f64, right: f64, area: f64)
        -> Result<DarcyGasFlowProfile, MaterialError> {
        if self.maximum_slip_pa > 0.0 {
            return Err(invalid("gas-specific Klinkenberg profile requires mass_flow_profile_for_model"));
        }
        let evaluated = self.evaluate(gas, left, right, area)?;
        self.profile(left, right, area, evaluated)
    }

    /// Gas-calibrated flow and interface pressures, using the same accepted
    /// flux as the ordinary query. No averaged layer slip factor is introduced.
    pub fn mass_flow_profile_for_model(&self, gas: &GasState, model: &ThermochemicalGasModel,
        left: f64, right: f64, area: f64) -> Result<DarcyGasFlowProfile, MaterialError> {
        if self.maximum_slip_pa > 0.0 && self.slip_gas_model.as_ref() != Some(model) {
            return Err(invalid("Klinkenberg calibration does not cover this gas model/composition"));
        }
        let evaluated = self.evaluate(gas, left, right, area)?;
        self.profile(left, right, area, evaluated)
    }

    fn profile(&self, left: f64, right: f64, area: f64,
        (flow, root): (DarcyGasFlow, Option<LayeredRoot>)) -> Result<DarcyGasFlowProfile, MaterialError> {
        let mut result = DarcyGasFlowProfile { flow, interface_pressures_pa: None,
            iterations: root.map_or(0, |root| root.iterations), normalized_pressure_residual: 0.0 };
        if self.geometric_resistance_per_m.is_none() || area == 0.0 { return Ok(result); }
        let mut pressures = vec![left; self.layers.len()+1];
        if left != right {
            let query = LayeredQuery::new(left, right, self.maximum_slip_pa)?;
            let scaled_flux = if let Some(root) = root { root.scaled_flux } else {
                (0.5*left/query.scale + 0.5*right/query.scale)
                    + self.common_slip_pa.ok_or_else(|| invalid("missing laminate flux solve"))?/query.scale
            };
            let marched = self.march(query, scaled_flux, false, Some(&mut pressures))?
                .ok_or_else(|| invalid("laminate pressure reconstruction failed"))?;
            result.normalized_pressure_residual = marched.fraction-1.0;
            if result.normalized_pressure_residual.abs() > self.pressure_residual_tolerance() {
                return Err(invalid("laminate interface-pressure residual exceeds its bound"));
            }
            // Preserve the exact prescribed endpoints; the separately retained
            // scaled residual discloses reconstruction roundoff before this.
            pressures[0] = left;
            pressures[self.layers.len()] = right;
        }
        result.interface_pressures_pa = Some(pressures);
        Ok(result)
    }

    fn pressure_residual_tolerance(&self) -> f64 {
        64.0*f64::EPSILON*self.layers.len() as f64
    }

    // lambda = |j| mu R T sum(L/k) / (delta_p * pressure_scale).
    // Its exact bracket uses min/max b, since sum(b_i delta_p_i) lies
    // between min(b)*delta_p and max(b)*delta_p for monotone pressure.
    fn solve_layered(&self, query: LayeredQuery, budget: usize) -> Result<LayeredRoot, MaterialError> {
        let mean_scaled = query.low_scaled+0.5*query.drop_scaled;
        let mut lower = mean_scaled+self.minimum_slip_pa/query.scale;
        let mut upper = mean_scaled+self.maximum_slip_pa/query.scale;
        let mut candidate = if upper > 4.0*lower {
            fs_math::det::sqrt(lower)*fs_math::det::sqrt(upper)
        } else { 0.5*lower+0.5*upper };
        for iteration in 1..=budget.min(MAX_LAYERED_DARCY_ITERATIONS) {
            let marched = self.march(query, candidate, true, None)?;
            let mut newton = f64::NAN;
            if let Some(marched) = marched {
                if marched.fraction == 1.0 {
                    return Ok(LayeredRoot { scaled_flux: candidate, iterations: iteration });
                }
                lower = candidate;
                if positive(marched.derivative) {
                    newton = candidate+(1.0-marched.fraction)/marched.derivative;
                }
            } else { upper = candidate; }
            if (upper-lower)/lower <= 64.0*f64::EPSILON {
                let solved = 0.5*lower+0.5*upper;
                let residual = self.march(query, solved, false, None)?
                    .ok_or_else(|| invalid("laminate flux reconstruction failed"))?.fraction-1.0;
                if residual.abs() <= self.pressure_residual_tolerance() {
                    return Ok(LayeredRoot { scaled_flux: solved, iterations: iteration });
                }
                return Err(invalid("laminate scalar flux bracket does not close its pressure residual"));
            }
            candidate = if upper > 4.0*lower {
                // Geometric bisection handles extreme positive slip contrasts
                // without hundreds of additive halvings near a tiny root.
                fs_math::det::sqrt(lower)*fs_math::det::sqrt(upper)
            } else {
                let width = upper-lower;
                if newton > lower+0.25*width && newton < upper-0.25*width { newton }
                else { 0.5*lower+0.5*upper }
            };
            if !(candidate > lower && candidate < upper) {
                return Err(invalid("laminate scalar flux solve stagnated before acceptance"));
            }
        }
        Err(invalid("laminate scalar flux solve exhausted its bounded iteration budget"))
    }

    // Integrate from the lower-pressure face. The rationalized positive root
    // solves c = a*ds + drop_scaled*ds^2/2, without pressure-square subtraction
    // or a negative radicand. Trial overshoots stop before leaving the physical
    // pressure interval; only an accepted flux is reconstructed in full.
    fn march(&self, query: LayeredQuery, scaled_flux: f64, stop_at_high: bool,
        mut interfaces: Option<&mut [f64]>) -> Result<Option<PressureMarch>, MaterialError> {
        let mut fraction = 0.0;
        let mut correction = 0.0;
        let mut derivative = 0.0;
        for ordinal in 0..self.layers.len() {
            let i = if query.reverse_layers { self.layers.len()-1-ordinal } else { ordinal };
            let weight = self.resistance_fractions[i];
            let c = scaled_flux*weight;
            let a = query.low_scaled+self.layers[i].klinkenberg_slip_pa/query.scale
                + query.drop_scaled*fraction;
            if stop_at_high {
                let remaining = (1.0-fraction).max(0.0);
                if c > remaining*(a+0.5*query.drop_scaled*remaining) { return Ok(None); }
            }
            if weight > 0.0 && !positive(c) {
                return Err(invalid("laminate layer pressure increment underflowed"));
            }
            let cross = fs_math::det::sqrt(query.drop_scaled)*fs_math::det::sqrt(c);
            let scale = a.max(cross);
            let root = scale*fs_math::det::sqrt((a/scale)*(a/scale)+2.0*(cross/scale)*(cross/scale));
            let increment = c/(0.5*a+0.5*root);
            if !nonnegative(increment) || !positive(a) {
                return Err(invalid("laminate layer pressure propagation is unrepresentable"));
            }
            let end_a = a+query.drop_scaled*increment;
            derivative = (a/end_a)*derivative+weight/end_a;
            let adjusted = increment-correction;
            let next = fraction+adjusted;
            correction = (next-fraction)-adjusted;
            fraction = next;
            if let Some(pressures) = interfaces.as_deref_mut() {
                let face = if query.reverse_layers { i } else { i+1 };
                pressures[face] = query.low+query.drop*fraction.clamp(0.0,1.0);
            }
        }
        if !nonnegative(fraction) {
            return Err(invalid("laminate total pressure propagation is unrepresentable"));
        }
        if stop_at_high && fraction > 1.0 { return Ok(None); }
        Ok(Some(PressureMarch { fraction, derivative }))
    }

    fn evaluate(&self, gas: &GasState, left: f64, right: f64, area: f64)
        -> Result<(DarcyGasFlow, Option<LayeredRoot>), MaterialError> {
        self.support.admit(gas, left, right)?;
        if !nonnegative(area) { return Err(invalid("Darcy flow area must be finite and nonnegative")); }
        let mut flow = DarcyGasFlow { mass_rate_left_to_right_kg_s: 0.0,
            mass_flux_left_to_right_kg_m2_s: 0.0, left_volume_rate_m3_s: 0.0,
            right_volume_rate_m3_s: 0.0, mean_volume_conductance_m3_pa_s: 0.0,
            drag_dissipation_w: 0.0, temperature_k: gas.temperature,
            left_pressure_pa: left, right_pressure_pa: right };
        let Some(geometry) = self.geometric_resistance_per_m.filter(|_| area > 0.0) else { return Ok((flow, None)); };
        let mean = 0.5*left + 0.5*right;
        let rt = gas.specific_gas_constant*gas.temperature;
        let hydraulic = gas.dynamic_viscosity*geometry;
        let conductance = 1.0/hydraulic;
        if ![rt, hydraulic, conductance].into_iter().all(positive) {
            return Err(invalid("gas permeation coefficients are unrepresentable"));
        }
        let difference = left-right;
        let mut solved = None;
        let (enhancement, pressure_over_rt) = if let Some(slip) = self.common_slip_pa {
            (1.0+slip/mean, (mean+slip)/rt)
        } else if difference == 0.0 {
            // At equilibrium every interface has the same pressure, so the
            // exact infinitesimal series resistance is sum[r_i/(1+b_i/p)].
            let scale = mean.max(self.maximum_slip_pa);
            let p = mean/scale;
            let fraction: f64 = self.layers.iter().zip(&self.resistance_fractions).map(|(layer, &weight)|
                weight*(p/(p+layer.klinkenberg_slip_pa/scale))).sum();
            (1.0/fraction, 0.0)
        } else {
            let query = LayeredQuery::new(left, right, self.maximum_slip_pa)?;
            let root = self.solve_layered(query, MAX_LAYERED_DARCY_ITERATIONS)?;
            solved = Some(root);
            (root.scaled_flux*(query.scale/mean), root.scaled_flux*(query.scale/rt))
        };
        if !positive(enhancement) || difference != 0.0 && !positive(pressure_over_rt) {
            return Err(invalid("gas permeation slip coefficients are unrepresentable"));
        }
        flow.mean_volume_conductance_m3_pa_s = area*conductance*enhancement;
        if difference != 0.0 {
            let flux = (difference*conductance)*pressure_over_rt;
            let left_speed = flux*rt/left;
            let right_speed = flux*rt/right;
            if !flux.is_finite() || flux == 0.0
                || !left_speed.is_finite() || !right_speed.is_finite()
                || left_speed.abs().max(right_speed.abs()) > self.support.maximum_superficial_speed_m_s {
                return Err(invalid("gas permeation exceeds its linear-speed support or numerical range"));
            }
            flow.mass_flux_left_to_right_kg_m2_s = flux;
            flow.mass_rate_left_to_right_kg_s = area*flux;
            flow.left_volume_rate_m3_s = area*left_speed;
            flow.right_volume_rate_m3_s = area*right_speed;
            // u = j RT/p. This drag-work diagnostic is not an extra heat input.
            let (high, low) = (left.max(right), left.min(right));
            let ratio = (high-low)/low;
            let logarithm = if ratio < 1e-4 {
                // log(1+x), relative remainder < 2e-21 on this interval.
                ratio*(1.0+ratio*(-0.5+ratio*(1.0/3.0+ratio*(-0.25+ratio*0.2))))
            } else if ratio.is_finite() { fs_math::det::ln(1.0+ratio) }
            else { fs_math::det::ln(high)-fs_math::det::ln(low) };
            flow.drag_dissipation_w = flow.mass_rate_left_to_right_kg_s.abs()*rt*logarithm;
            if flow.mass_rate_left_to_right_kg_s == 0.0 {
                return Err(invalid("nonzero gas mass transport underflowed"));
            }
        }
        if ![flow.mass_rate_left_to_right_kg_s, flow.mass_flux_left_to_right_kg_m2_s,
            flow.left_volume_rate_m3_s, flow.right_volume_rate_m3_s,
            flow.mean_volume_conductance_m3_pa_s, flow.drag_dissipation_w].iter().all(|x| x.is_finite())
            || flow.mean_volume_conductance_m3_pa_s <= 0.0 || flow.drag_dissipation_w < 0.0 {
            return Err(invalid("gas permeation flow/power is unrepresentable"));
        }
        Ok((flow, solved))
    }
}

#[cfg(test)]
mod tests;
