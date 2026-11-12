//! Original two-voltage channel-current surfaces and explicit thermal scaling.
//!
//! Current is drain-to-source positive. The source is channel-only: body-diode
//! conduction remains with its separate original device owner. Original affine
//! source cells supply both current and analytical partial derivatives; no
//! nominal drain voltage, I/V quotient, fitted grid or extrapolation is used.

use super::*;
use fs_matdb::{
    EnvelopeAnswer, InterpolationPolicy, PreparedSurfaceResponse, PropertyClaim, PropertyValue,
    QueryPoint,
};
use fs_phs::PhsError;
use fs_phs::circuit::converter::commutation::gate::{GateChannelLaw, GateChannelPoint};

/// Explicit positive dimensionless multiplier of the complete current surface.
/// This is an authored separable law, not an inferred temperature correction.
/// Both this claim and the current claim must apply at the actual temperature.
#[derive(Clone, Debug)]
pub struct GateCurrentThermalScale<'a> {
    pub key: &'a PropertyKey,
    pub selection: MaterialPropertySelection,
}

/// Complete channel-only I(Vds,Vgs), retaining both signed voltage domains.
/// The original material context must already own both distinct typed volt axes.
/// An authored zero-drain mesh line and same-sign cells certify passive current
/// throughout the original surface. Cutoff, saturation and passive negative
/// differential conductance are admitted.
#[derive(Clone, Debug)]
pub struct GateCurrentSurfaceBinding<'a> {
    pub cell: CellId,
    pub key: &'a PropertyKey,
    pub drain_voltage_axis: &'a str,
    pub gate_voltage_axis: &'a str,
    pub drain_voltage_range_v: [f64; 2],
    pub gate_voltage_range_v: [f64; 2],
    pub thermal_scale: Option<GateCurrentThermalScale<'a>>,
    pub selection: MaterialPropertySelection,
}

#[derive(Clone, Debug)]
struct CurrentSurfaceData {
    claim: PropertyClaim,
    response: PreparedSurfaceResponse,
    drain_index: usize,
    drain_range_v: [f64; 2],
    gate_range_v: [f64; 2],
}
impl CurrentSurfaceData {
    fn same_response(&self, other: &Self) -> bool {
        self.drain_index == other.drain_index
            && self.drain_range_v.map(f64::to_bits) == other.drain_range_v.map(f64::to_bits)
            && self.gate_range_v.map(f64::to_bits) == other.gate_range_v.map(f64::to_bits)
            && self.response.axes() == other.response.axes()
            && self.response.samples() == other.response.samples()
            && self.response.triangles() == other.response.triangles()
    }
}
impl PartialEq for CurrentSurfaceData {
    fn eq(&self, other: &Self) -> bool {
        self.claim == other.claim && self.same_response(other)
    }
}

#[derive(Clone, Debug)]
struct ScaleSource {
    key: PropertyKey,
    source: ResolvedMaterialStatePoint,
    value: f64,
}
/// Equal = the same material state (card, material, axes, query point, requirement, source
/// axes, selected claim, value bits) behind the same key and multiplier. The receipt's
/// selection-policy tag and hash are left out on purpose, because a refresh at the same state
/// re-resolves under `PinnedByProperty`.
impl PartialEq for ScaleSource {
    fn eq(&self, other: &Self) -> bool {
        self.key == other.key
            && self.value == other.value
            && same_material(&self.source, &other.source)
    }
}

/// Fresh original source receipts plus one immutable prepared source surface.
/// The surface and its triangulation stay fixed along the history. Only an
/// explicitly bound positive thermal multiplier changes the current law.
#[derive(Clone, Debug)]
pub struct GateCurrentSurfaceMaterial {
    cell: CellId,
    key: PropertyKey,
    source: ResolvedMaterialStatePoint,
    support: EnvelopeAnswer,
    data: Arc<CurrentSurfaceData>,
    scale: Option<ScaleSource>,
    law: GateChannelLaw,
}
/// Equal = the same material state (card, material, axes, query point, requirement, source
/// axes, selected claim, value bits) behind the same cell, key, support, multiplier source and
/// surface data. The receipt's selection-policy tag and hash are left out of the source on
/// purpose, because a refresh at the same state re-resolves under `PinnedByProperty`.
impl PartialEq for GateCurrentSurfaceMaterial {
    fn eq(&self, other: &Self) -> bool {
        self.cell == other.cell
            && self.key == other.key
            && same_material(&self.source, &other.source)
            && self.support == other.support
            && self.scale == other.scale
            && (Arc::ptr_eq(&self.data, &other.data) || self.data == other.data)
    }
}
impl GateCurrentSurfaceMaterial {
    pub const fn cell(&self) -> CellId {
        self.cell
    }
    pub const fn key(&self) -> &PropertyKey {
        &self.key
    }
    pub const fn source(&self) -> &ResolvedMaterialStatePoint {
        &self.source
    }
    pub fn source_claim(&self) -> &PropertyClaim {
        &self.data.claim
    }
    pub fn response(&self) -> &PreparedSurfaceResponse {
        &self.data.response
    }
    pub const fn support(&self) -> &EnvelopeAnswer {
        &self.support
    }
    pub fn drain_voltage_axis(&self) -> &str {
        &self.data.response.axes()[self.data.drain_index]
    }
    pub fn gate_voltage_axis(&self) -> &str {
        &self.data.response.axes()[1 - self.data.drain_index]
    }
    pub fn drain_voltage_range_v(&self) -> [f64; 2] {
        self.data.drain_range_v
    }
    pub fn gate_voltage_range_v(&self) -> [f64; 2] {
        self.data.gate_range_v
    }
    pub fn thermal_scale_source(&self) -> Option<&ResolvedMaterialStatePoint> {
        self.scale.as_ref().map(|scale| &scale.source)
    }
    pub fn thermal_scale(&self) -> f64 {
        self.scale.as_ref().map_or(1.0, |scale| scale.value)
    }
}

/// Numerical evaluation owns the same prepared material-data table. Its
/// original source receipts are refreshed by the thermal owner before use.
#[derive(Clone, Debug)]
struct CurrentResponse {
    data: Arc<CurrentSurfaceData>,
    scale: f64,
}
impl GateChannelResponse for CurrentResponse {
    fn gate_voltage_range_v(&self) -> [f64; 2] {
        self.data.gate_range_v
    }
    fn evaluate(&self, vds: f64, vgs: f64) -> Result<GateChannelPoint, PhsError> {
        if !vds.is_finite()
            || !vgs.is_finite()
            || !(self.data.drain_range_v[0]..=self.data.drain_range_v[1]).contains(&vds)
            || !(self.data.gate_range_v[0]..=self.data.gate_range_v[1]).contains(&vgs)
        {
            return Err(PhsError::Unadmitted {
                what: "actual channel voltages left the original two-voltage source support",
            });
        }
        let drain = self.data.drain_index;
        let mut at = [0.0; 2];
        at[drain] = vds;
        at[1 - drain] = vgs;
        let value = self
            .data
            .response
            .evaluate(at)
            .map_err(|_| PhsError::Unadmitted {
                what: "actual channel state has no original surface-cell response and tangent",
            })?;
        let originals = [
            value.value,
            value.gradient[drain],
            value.gradient[1 - drain],
        ];
        let scaled = originals.map(|value| value * self.scale);
        if scaled
            .iter()
            .zip(originals)
            .any(|(&after, before)| !after.is_finite() || (before != 0.0 && after == 0.0))
        {
            return Err(PhsError::Unadmitted {
                what: "thermally scaled channel current or source-cell derivative is unrepresentable",
            });
        }
        // The canonical evaluator preserves an exact constant-valued edge,
        // including the certified I=0 row. No tolerance or clipping is needed.
        Ok(GateChannelPoint {
            current_a: scaled[0],
            drain_derivative_s: scaled[1],
            gate_derivative_s: scaled[2],
        })
    }
}

impl channel_sealed::Material for GateCurrentSurfaceMaterial {}
impl GateChannelMaterial for GateCurrentSurfaceMaterial {
    type Response = GateChannelLaw;
    fn native_response(&self) -> Self::Response {
        self.law.clone()
    }
    fn same_electrical_response(&self, other: &Self) -> bool {
        self.thermal_scale().to_bits() == other.thermal_scale().to_bits()
            && (Arc::ptr_eq(&self.data, &other.data) || self.data.same_response(&other.data))
    }
}
impl channel_sealed::Binding for GateCurrentSurfaceBinding<'_> {}
impl GateChannelBinding for GateCurrentSurfaceBinding<'_> {
    type Material = GateCurrentSurfaceMaterial;
    fn cell(&self) -> CellId {
        self.cell
    }
    fn resolve(
        cells: &[MaterialCell<'_>],
        bindings: &[Self; 4],
        thermal: &DepositionSimulation<'_>,
        sample: &MaterialNetworkSample,
        retained: Option<&[Self::Material; 4]>,
        cancel: &RefCell<impl FnMut() -> bool>,
    ) -> Result<[Self::Material; 4], Error> {
        let mut resolved = Vec::with_capacity(4);
        for (index, binding) in bindings.iter().enumerate() {
            resolved.push(resolve_one(
                cells,
                binding,
                thermal,
                sample,
                retained.map(|old| &old[index]),
                cancel,
            )?);
        }
        check(cancel)?;
        resolved
            .try_into()
            .map_err(|_| invalid("four complete current-surface sources required"))
    }
}

fn same_held_context(
    old: &ResolvedMaterialStatePoint,
    query: &QueryPoint,
    body: &MaterialCell<'_>,
    owned: &[&str],
) -> Result<(), Error> {
    for (axis, value) in old.query_point() {
        if axis != body.point.temperature_axis
            && Some(axis.as_str()) != body.point.liquid_fraction_axis
            && !owned.contains(&axis.as_str())
            && query.axes().get(axis).map(|v| v.to_bits()) != Some(value.to_bits())
        {
            return Err(invalid(
                "channel current changed its held original source context",
            ));
        }
    }
    Ok(())
}

fn actual_source(
    body: &MaterialCell<'_>,
    key: &PropertyKey,
    query: &QueryPoint,
    admissibility: ScalarAdmissibility,
    selection: &MaterialPropertySelection,
    retained: Option<&ResolvedMaterialStatePoint>,
    owned: &[&str],
) -> Result<ResolvedMaterialStatePoint, Error> {
    if let Some(old) = retained {
        same_held_context(old, query, body, owned)?;
        old.at_query_point(body.point.material, query)
            .map_err(source)
    } else {
        let requirement =
            ScalarPropertyRequirement::try_with_key(key, admissibility).map_err(source)?;
        resolve_material_state_point(
            body.point.material,
            query,
            &[requirement],
            selection.clone(),
        )
        .map_err(source)
    }
}

fn owned_thermal<'a>(
    body: &'a MaterialCell<'_>,
    query: &QueryPoint,
) -> Vec<(&'a str, QuantitySpec)> {
    let mut owned = Vec::with_capacity(2);
    if let Some(quantity) = query.axis_quantities().get(body.point.temperature_axis) {
        owned.push((body.point.temperature_axis, *quantity));
    }
    if let Some(axis) = body.point.liquid_fraction_axis {
        owned.push((axis, QuantitySpec::dimensional(Dims::NONE)));
    }
    owned
}

fn resolve_scale(
    body: &MaterialCell<'_>,
    binding: Option<&GateCurrentThermalScale<'_>>,
    query: &QueryPoint,
    retained: Option<&ScaleSource>,
) -> Result<Option<ScaleSource>, Error> {
    let Some(binding) = binding else {
        if retained.is_some() {
            return Err(invalid(
                "channel current cannot drop its original thermal multiplier",
            ));
        }
        return Ok(None);
    };
    if binding.key.quantity() != QuantitySpec::dimensional(Dims::NONE)
        || retained.is_some_and(|old| &old.key != binding.key)
    {
        return Err(invalid(
            "channel thermal multiplier needs its original dimensionless property",
        ));
    }
    let source = actual_source(
        body,
        binding.key,
        query,
        ScalarAdmissibility::StrictlyPositive,
        &binding.selection,
        retained.map(|old| &old.source),
        &[],
    )?;
    let property = source
        .property_by_key(binding.key)
        .ok_or_else(|| invalid("channel current lost its thermal multiplier source"))?;
    body.point.admit_static(property).map_err(point)?;
    if !property
        .has_only_owned_coordinates_and_fixed_external_conditions(&owned_thermal(body, query))
    {
        return Err(invalid(
            "separable channel multiplier can vary only with owned thermal coordinates or fixed external conditions",
        ));
    }
    let value = property.value_si();
    Ok(Some(ScaleSource {
        key: binding.key.clone(),
        source,
        value,
    }))
}

fn prepare_data(
    claim: &PropertyClaim,
    binding: &GateCurrentSurfaceBinding<'_>,
    cancel: &RefCell<impl FnMut() -> bool>,
) -> Result<CurrentSurfaceData, Error> {
    let PropertyValue::Surface2D {
        axes,
        axis_dims,
        samples,
        triangles,
        dims,
    } = &claim.value
    else {
        return Err(invalid(
            "channel current requires one complete original two-voltage surface",
        ));
    };
    let drain = axes
        .iter()
        .position(|axis| axis == binding.drain_voltage_axis)
        .ok_or_else(|| invalid("channel current source lacks its owned drain voltage"))?;
    if axes[1 - drain] != binding.gate_voltage_axis
        || *axis_dims != [Voltage::DIMS; 2]
        || *dims != Current::DIMS
        || claim.interpolation != InterpolationPolicy::LinearInside
    {
        return Err(invalid(
            "channel current needs original ampere/volt affine cells on both owned voltage axes",
        ));
    }
    let response = claim.value.prepare_surface_response().map_err(|error| {
        source(fs_material::state_point::MaterialStatePointError::Query {
            property: binding.key.name().to_owned(),
            source: error,
        })
    })?;
    let mut incident = vec![false; samples.len()];
    for triangle in triangles {
        check(cancel)?;
        let mut low = f64::INFINITY;
        let mut high = f64::NEG_INFINITY;
        for &index in triangle {
            incident[index as usize] = true;
            let vds = samples[index as usize].0[drain];
            low = low.min(vds);
            high = high.max(vds);
        }
        if low < 0.0 && high > 0.0 {
            return Err(invalid(
                "passive channel source cells must meet on an authored zero-drain mesh; crossing cells are not split or inferred",
            ));
        }
    }
    for (index, &(at, current)) in samples.iter().enumerate() {
        check(cancel)?;
        let vds = at[drain];
        if (vds == 0.0 && current != 0.0)
            || (vds > 0.0 && current < 0.0)
            || (vds < 0.0 && current > 0.0)
        {
            return Err(invalid(
                "original channel source violates signed passive current or exact zero-drain current",
            ));
        }
        if !incident[index]
            && (binding.drain_voltage_range_v[0]..=binding.drain_voltage_range_v[1]).contains(&vds)
            && (binding.gate_voltage_range_v[0]..=binding.gate_voltage_range_v[1])
                .contains(&at[1 - drain])
        {
            return Err(invalid(
                "channel operating rectangle contains an isolated source sample without a tangent",
            ));
        }
    }
    Ok(CurrentSurfaceData {
        claim: claim.clone(),
        response,
        drain_index: drain,
        drain_range_v: binding.drain_voltage_range_v,
        gate_range_v: binding.gate_voltage_range_v,
    })
}

fn resolve_one(
    cells: &[MaterialCell<'_>],
    binding: &GateCurrentSurfaceBinding<'_>,
    thermal: &DepositionSimulation<'_>,
    sample: &MaterialNetworkSample,
    retained: Option<&GateCurrentSurfaceMaterial>,
    cancel: &RefCell<impl FnMut() -> bool>,
) -> Result<GateCurrentSurfaceMaterial, Error> {
    check(cancel)?;
    let voltage = QuantitySpec::dimensional(Voltage::DIMS);
    let [drain_low, drain_high] = binding.drain_voltage_range_v;
    let [gate_low, gate_high] = binding.gate_voltage_range_v;
    if binding.key.quantity() != QuantitySpec::dimensional(Current::DIMS)
        || ![drain_low, drain_high, gate_low, gate_high]
            .iter()
            .all(|v| v.is_finite())
        || drain_low >= 0.0
        || drain_high <= 0.0
        || gate_low >= gate_high
        || binding.drain_voltage_axis == binding.gate_voltage_axis
    {
        return Err(invalid(
            "channel current needs an ampere property, two volt axes and finite signed-drain/gate support",
        ));
    }
    let body = cells
        .iter()
        .find(|cell| cell.cell == binding.cell)
        .ok_or_else(|| invalid("channel current needs its original semiconductor package"))?;
    let cell = thermal
        .cells()
        .get(&binding.cell)
        .ok_or_else(|| invalid("channel current package is inactive"))?;
    if cell.material_id != body.material_id {
        return Err(invalid(
            "channel current changed its original caloric material",
        ));
    }
    let axes = [binding.drain_voltage_axis, binding.gate_voltage_axis];
    if axes.iter().any(|&axis| {
        axis == body.point.temperature_axis
            || Some(axis) == body.point.liquid_fraction_axis
            || body.point.context.axis_quantities().get(axis) != Some(&voltage)
            || !body.point.context.axes().contains_key(axis)
    }) {
        return Err(invalid(
            "channel current requires two distinct existing typed voltage coordinates",
        ));
    }
    if let Some(old) = retained {
        if old.cell != binding.cell
            || &old.key != binding.key
            || old.drain_voltage_axis() != binding.drain_voltage_axis
            || old.gate_voltage_axis() != binding.gate_voltage_axis
            || old.drain_voltage_range_v().map(f64::to_bits)
                != binding.drain_voltage_range_v.map(f64::to_bits)
            || old.gate_voltage_range_v().map(f64::to_bits)
                != binding.gate_voltage_range_v.map(f64::to_bits)
            || old.scale.is_some() != binding.thermal_scale.is_some()
        {
            return Err(invalid(
                "channel current must retain its original package, key, two axes, support and multiplier",
            ));
        }
    }
    let phase = sample
        .cells
        .iter()
        .find(|row| row.cell == binding.cell)
        .ok_or_else(|| invalid("channel current lacks its admitted caloric phase"))?
        .phase;
    let query = body
        .point
        .at(phase)
        .map_err(point)?
        .with(binding.drain_voltage_axis, 0.0)
        .and_then(|q| q.with(binding.gate_voltage_axis, gate_low))
        .map_err(|_| invalid("channel source reference voltages are unrepresentable"))?;
    let material = actual_source(
        body,
        binding.key,
        &query,
        ScalarAdmissibility::Finite,
        &binding.selection,
        retained.map(|old| &old.source),
        &axes,
    )?;
    let property = material
        .property_by_key(binding.key)
        .ok_or_else(|| invalid("channel current lost its original source property"))?;
    let voltages = [
        (binding.drain_voltage_axis, voltage),
        (binding.gate_voltage_axis, voltage),
    ];
    body.point
        .admit_static_with_axes(property, &voltages)
        .map_err(point)?;
    let mut owned = owned_thermal(body, &query);
    owned.extend(voltages);
    if !property.has_only_owned_coordinates_and_fixed_external_conditions(&owned) {
        return Err(invalid(
            "channel surface cannot freeze another evolving bias, frequency or response coordinate",
        ));
    }
    let selected = property.answer().receipt.selected;
    let claim = body
        .point
        .material
        .claims()
        .claim(selected)
        .ok_or_else(|| invalid("original channel surface claim is absent"))?;
    let lower = query
        .clone()
        .with(binding.drain_voltage_axis, drain_low)
        .map_err(|_| invalid("channel lower drain voltage is unrepresentable"))?;
    let upper = query
        .clone()
        .with(binding.drain_voltage_axis, drain_high)
        .and_then(|q| q.with(binding.gate_voltage_axis, gate_high))
        .map_err(|_| invalid("channel upper source voltages are unrepresentable"))?;
    // The original envelope owner proves continuous support over the whole
    // rectangle, including its interior; endpoint samples alone are insufficient.
    let support = body
        .point
        .material
        .claims()
        .query_envelope_typed(
            binding.key,
            &lower,
            &upper,
            ClaimSelection::Pinned(selected),
        )
        .map_err(|error| {
            source(fs_material::state_point::MaterialStatePointError::Support {
                property: binding.key.name().to_owned(),
                source: Box::new(error),
            })
        })?;
    let data = if let Some(old) = retained {
        if &old.data.claim != claim {
            return Err(invalid(
                "channel current changed its original full source surface",
            ));
        }
        old.data.clone()
    } else {
        Arc::new(prepare_data(claim, binding, cancel)?)
    };
    let scale = resolve_scale(
        body,
        binding.thermal_scale.as_ref(),
        &query,
        retained.and_then(|old| old.scale.as_ref()),
    )?;
    check(cancel)?;
    let factor = scale.as_ref().map_or(1.0, |scale| scale.value);
    let law = if let Some(old) = retained
        && old.thermal_scale().to_bits() == factor.to_bits()
    {
        old.law.clone()
    } else {
        GateChannelLaw::new(CurrentResponse {
            data: data.clone(),
            scale: factor,
        })
        .map_err(model)?
    };
    Ok(GateCurrentSurfaceMaterial {
        cell: binding.cell,
        key: binding.key.clone(),
        source: material,
        support,
        data,
        scale,
        law,
    })
}
