//! Analytic tangents of the actual selected thermochemical source snapshot.

use super::{MaterialError, Temperature, invalid};
use crate::state_point::{ResolvedMaterialStatePoint, THERMAL_CONDUCTIVITY_PROPERTY};
use crate::thermochemical_gas::{
    ResolvedThermochemicalGasState, SpeciesTransportMaterial, ThermochemicalGasModel,
};
use fs_matdb::{InterpolationPolicy, PropertyValue};
use fs_thermochem::transport::MAX_WILKE_GAS_COMPONENTS_V1;

pub(super) fn derivatives(
    model: &ThermochemicalGasModel,
    snapshot: &ResolvedThermochemicalGasState,
) -> Result<[f64; 3], MaterialError> {
    let count = model.transport_materials().len();
    if count > MAX_WILKE_GAS_COMPONENTS_V1 || snapshot.components().len() != count {
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
        viscosity_derivatives[i] = source_derivative(source, resolved, "dynamic_viscosity", t)?;
        conductivity_derivatives[i] =
            source_derivative(source, resolved, THERMAL_CONDUCTIVITY_PROPERTY, t)?;
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
) -> Result<f64, MaterialError> {
    let property = resolved
        .property(name)
        .ok_or(invalid("original gas transport property is absent"))?;
    let claim = source
        .material
        .claims()
        .claim(property.answer().receipt.selected)
        .ok_or(invalid("original selected gas transport claim is absent"))?;
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
        _ => Err(invalid(
            "gas transport tangent requires a constant with thermal validity or an original linear temperature curve",
        )),
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
