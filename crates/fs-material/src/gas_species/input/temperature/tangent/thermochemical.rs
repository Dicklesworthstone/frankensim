//! Analytic tangents of the actual selected thermochemical source snapshot.

use super::{MaterialError, Temperature, invalid};
use crate::state_point::{ResolvedMaterialStatePoint, THERMAL_CONDUCTIVITY_PROPERTY};
use crate::thermochemical_gas::{
    ResolvedThermochemicalGasState, SpeciesTransportMaterial, ThermochemicalGasModel,
};
use fs_matdb::{
    ClaimId, InterpolationPolicy, PreparedSurfaceResponse, PropertyClaim, PropertyValue,
};
use fs_qty::Pressure;
use fs_thermochem::transport::MAX_WILKE_GAS_COMPONENTS_V1;

const TRANSPORT_PROPERTIES: [&str; 2] = ["dynamic_viscosity", THERMAL_CONDUCTIVITY_PROPERTY];

#[derive(Debug)]
pub(super) struct PreparedTransport {
    properties: Vec<[Option<PreparedSurface>; 2]>,
}

#[derive(Debug)]
struct PreparedSurface {
    selected: ClaimId,
    response: PreparedSurfaceResponse,
}

impl PreparedTransport {
    pub(super) fn try_new(
        model: &ThermochemicalGasModel,
        snapshot: &ResolvedThermochemicalGasState,
    ) -> Result<Option<Self>, MaterialError> {
        let count = model.transport_materials().len();
        if count > MAX_WILKE_GAS_COMPONENTS_V1 || snapshot.components().len() != count {
            return Err(invalid(
                "gas tangent preparation needs its original bounded sources",
            ));
        }
        // Do not allocate a cache for the existing scalar/curve-only path.
        let mut properties = Vec::new();
        for (i, (source, resolved)) in model
            .transport_materials()
            .iter()
            .zip(snapshot.components())
            .enumerate()
        {
            let mut entries = [None, None];
            for (entry, name) in entries.iter_mut().zip(TRANSPORT_PROPERTIES) {
                let (selected, claim) = selected_claim(source, resolved, name)?;
                if matches!(&claim.value, PropertyValue::Surface2D { .. }) {
                    surface_temperature_axis(claim)?;
                    *entry = Some(PreparedSurface {
                        selected,
                        response: claim
                            .value
                            .prepare_surface_response()
                            .map_err(surface_error)?,
                    });
                }
            }
            if entries.iter().any(Option::is_some) && properties.is_empty() {
                properties.reserve_exact(count);
                properties.resize_with(i, || [None, None]);
            }
            if properties.capacity() != 0 {
                properties.push(entries);
            }
        }
        Ok((!properties.is_empty()).then_some(Self { properties }))
    }
}

pub(super) fn derivatives(
    model: &ThermochemicalGasModel,
    snapshot: &ResolvedThermochemicalGasState,
    prepared: Option<&PreparedTransport>,
) -> Result<[f64; 3], MaterialError> {
    let count = model.transport_materials().len();
    if count > MAX_WILKE_GAS_COMPONENTS_V1
        || snapshot.components().len() != count
        || prepared.is_some_and(|cache| cache.properties.len() != count)
    {
        return Err(invalid(
            "thermochemical tangent needs its original bounded component sources",
        ));
    }
    let t = snapshot.state().temperature;
    let mut viscosity_derivatives = [0.0; MAX_WILKE_GAS_COMPONENTS_V1];
    let mut conductivity_derivatives = [0.0; MAX_WILKE_GAS_COMPONENTS_V1];
    for (i, (source, resolved)) in model
        .transport_materials()
        .iter()
        .zip(snapshot.components())
        .enumerate()
    {
        viscosity_derivatives[i] = source_derivative(
            source,
            resolved,
            TRANSPORT_PROPERTIES[0],
            t,
            prepared.and_then(|cache| cache.properties[i][0].as_ref()),
        )?;
        conductivity_derivatives[i] = source_derivative(
            source,
            resolved,
            TRANSPORT_PROPERTIES[1],
            t,
            prepared.and_then(|cache| cache.properties[i][1].as_ref()),
        )?;
    }
    let wilke = snapshot
        .viscosity()
        .temperature_tangent(&viscosity_derivatives[..count])
        .map_err(|error| MaterialError::State {
            what: format!("thermochemical Wilke tangent: {error}"),
        })?;
    let fractions = model.thermodynamics().mole_composition().fractions();
    let denominators = snapshot.viscosity().receipt().denominator_bits();
    let mut cp_derivative = 0.0;
    let mut conductivity_derivative = 0.0;
    for (i, (species, resolved)) in model
        .thermodynamics()
        .components()
        .iter()
        .zip(snapshot.components())
        .enumerate()
    {
        let cp = species
            .evaluate_heat_capacity_temperature_derivative(Temperature::new(t))
            .map_err(|error| MaterialError::State {
                what: format!("thermochemical Cp tangent: {error}"),
            })?;
        cp_derivative += fractions[i] * cp.value();
        let pure = resolved
            .property(THERMAL_CONDUCTIVITY_PROPERTY)
            .expect("original gas resolved conductivity")
            .value_si();
        let denominator = f64::from_bits(denominators[i]);
        // The numerator is the ORIGINAL selected k source; the denominator and
        // its derivative are outputs of the ORIGINAL Wilke interaction owner.
        // No independent WMS interaction, source query or pair matrix is built.
        let contribution = fractions[i] * (pure / denominator);
        conductivity_derivative += fractions[i] * (conductivity_derivatives[i] / denominator)
            - contribution * (wilke.denominator_derivatives_per_k()[i] / denominator);
    }
    cp_derivative /= snapshot.thermodynamics().properties().molar_mass().value();
    if !cp_derivative.is_finite() || !conductivity_derivative.is_finite() {
        return Err(invalid(
            "thermochemical heat-capacity or conductivity tangent is unrepresentable",
        ));
    }
    Ok([
        cp_derivative,
        wilke.dynamic_viscosity_derivative_pa_s_k,
        conductivity_derivative,
    ])
}

// The complete value query has ALREADY selected and admitted this claim at
// the actual pressure/composition/context. Differentiate its original source
// coefficients without another material query, hashing pass or source switch.
fn source_derivative(
    source: &SpeciesTransportMaterial,
    resolved: &ResolvedMaterialStatePoint,
    name: &str,
    t: f64,
    prepared: Option<&PreparedSurface>,
) -> Result<f64, MaterialError> {
    let (selected, claim) = selected_claim(source, resolved, name)?;
    match &claim.value {
        PropertyValue::Scalar { .. }
            if claim.interpolation == InterpolationPolicy::ConstantWithinValidity
                && claim
                    .validity
                    .bound("temperature")
                    .is_some_and(|(lo, hi)| hi > lo) =>
        {
            Ok(0.0)
        }
        PropertyValue::Curve {
            abscissa,
            abscissa_dims,
            knots,
            ..
        } if abscissa == "temperature"
            && *abscissa_dims == Temperature::DIMS
            && claim.interpolation == InterpolationPolicy::LinearInside
            && knots.len() >= 2 =>
        {
            let (mut lo, mut hi) = (knots[0].0, knots[knots.len() - 1].0);
            if let Some((lower, upper)) = claim.validity.bound("temperature") {
                lo = lo.max(lower);
                hi = hi.min(upper);
            }
            if !lo.is_finite() || !hi.is_finite() || hi <= lo || t < lo || t > hi {
                return Err(invalid(
                    "gas source tangent needs a nonzero supported temperature interval",
                ));
            }
            let upper = knots.partition_point(|&(coordinate, _)| coordinate <= t);
            let mut index = (upper - 1).min(knots.len() - 2);
            if t == hi && index > 0 && t == knots[index].0 {
                index -= 1;
            }
            let derivative = segment_slope(knots[index], knots[index + 1])?;
            if t > lo && t < hi && index > 0 && t == knots[index].0 {
                let left = segment_slope(knots[index - 1], knots[index])?;
                if (derivative - left).abs()
                    > 128.0 * f64::EPSILON * derivative.abs().max(left.abs())
                {
                    return Err(invalid(
                        "gas transport source has no unique temperature derivative at this knot",
                    ));
                }
            }
            Ok(derivative)
        }
        PropertyValue::Surface2D { .. } => {
            let axis = surface_temperature_axis(claim)?;
            let pressure = resolved
                .query_point()
                .iter()
                .find_map(|(axis, value)| (axis == "pressure").then_some(*value))
                .ok_or_else(|| invalid("original gas surface pressure is absent"))?;
            let point = if axis == 0 {
                [t, pressure]
            } else {
                [pressure, t]
            };
            let allowed = claim
                .validity
                .bound("temperature")
                .map_or([true; 2], |(lo, hi)| [t > lo, t < hi]);
            let uncached;
            let response = if let Some(prepared) = prepared {
                if prepared.selected != selected {
                    return Err(invalid(
                        "prepared gas surface lost its original selected claim",
                    ));
                }
                &prepared.response
            } else {
                // The convenience API supports the same surfaces. Iterative
                // consumers retain a prepared view and avoid this cold work.
                uncached = claim
                    .value
                    .prepare_surface_response()
                    .map_err(surface_error)?;
                &uncached
            };
            response
                .inward_axis_derivative(point, axis, allowed)
                .map_err(surface_error)
        }
        _ => Err(invalid(
            "gas transport tangent requires a thermally supported constant, linear temperature curve or temperature/pressure surface",
        )),
    }
}

fn selected_claim<'s>(
    source: &'s SpeciesTransportMaterial,
    resolved: &ResolvedMaterialStatePoint,
    name: &str,
) -> Result<(ClaimId, &'s PropertyClaim), MaterialError> {
    let selected = resolved
        .property(name)
        .ok_or_else(|| invalid("original gas transport property is absent"))?
        .answer()
        .receipt
        .selected;
    let claim = source
        .material
        .claims()
        .claim(selected)
        .ok_or_else(|| invalid("original selected gas transport claim is absent"))?;
    Ok((selected, claim))
}

fn surface_temperature_axis(claim: &PropertyClaim) -> Result<usize, MaterialError> {
    if let PropertyValue::Surface2D {
        axes, axis_dims, ..
    } = &claim.value
        && claim.interpolation == InterpolationPolicy::LinearInside
    {
        for thermal in 0..2 {
            let pressure = 1 - thermal;
            if axes[thermal] == "temperature"
                && axis_dims[thermal] == Temperature::DIMS
                && axes[pressure] == "pressure"
                && axis_dims[pressure] == Pressure::DIMS
            {
                return Ok(thermal);
            }
        }
    }
    Err(invalid(
        "gas surface tangent requires original linear temperature/pressure cells",
    ))
}

fn surface_error(error: fs_matdb::MatDbError) -> MaterialError {
    MaterialError::State {
        what: format!("original gas transport surface tangent: {error}"),
    }
}

fn segment_slope(a: (f64, f64), b: (f64, f64)) -> Result<f64, MaterialError> {
    let width = b.0 - a.0;
    let difference = b.1 - a.1;
    let derivative = difference / width;
    if !width.is_finite()
        || width <= 0.0
        || !difference.is_finite()
        || !derivative.is_finite()
        || difference != 0.0 && derivative == 0.0
    {
        return Err(invalid(
            "original gas transport source tangent is unrepresentable",
        ));
    }
    Ok(derivative)
}
