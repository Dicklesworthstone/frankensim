//! Observation derivatives are separate from transport derivatives. In
//! particular T(h) may be smooth while the independently authored f_liquid(h)
//! has a corner. Do not assign a classical derivative by selecting one side.

use super::{EnthalpyAdjointError, ScanSensitivityError, poll};
use crate::transient::enthalpy::heterogeneous::HeterogeneousEnthalpyBackwardEuler;
use fs_exec::Cx;
use fs_material::phase::EquilibriumEnthalpyPhaseCurve;

#[derive(Debug)]
pub(super) enum FractionDerivative {
    Smooth(f64),
    Corner,
    Endpoint,
    Unrepresentable,
}

// Retain refusal metadata instead of refusing every existing h/T derivative.
// Only a requested nonzero fraction seed consumes this extra derivative.
fn derivative(curve: &EquilibriumEnthalpyPhaseCurve, h: f64) -> FractionDerivative {
    let Ok(slope) = curve.liquid_mass_fraction_derivative_at_specific_enthalpy(h) else {
        return FractionDerivative::Unrepresentable;
    };
    if !slope.is_finite() { return FractionDerivative::Unrepresentable; }
    let knots = curve.knots();
    if let Ok(index) = knots.binary_search_by(|k| k.specific_enthalpy_j_kg.partial_cmp(&h)
        .expect("finite admitted chart state")) {
        if index == 0 || index + 1 == knots.len() { return FractionDerivative::Endpoint; }
        let Ok(left) = curve.liquid_mass_fraction_derivative_at_specific_enthalpy(
            knots[index-1].specific_enthalpy_j_kg) else {
            return FractionDerivative::Unrepresentable;
        };
        if !left.is_finite() { return FractionDerivative::Unrepresentable; }
        // Exact comparison matches the transport owner's classical-derivative
        // admission: a tolerance would conceal a real constitutive corner.
        if left != slope { return FractionDerivative::Corner; }
    }
    FractionDerivative::Smooth(slope)
}

pub(super) fn bind(cx: &Cx<'_>, storage: &HeterogeneousEnthalpyBackwardEuler<'_, '_>, h: &[f64])
    -> Result<Vec<FractionDerivative>, ScanSensitivityError>
{
    poll(cx, 0)?;
    let mut out = Vec::new();
    out.try_reserve_exact(h.len()).map_err(|_| EnthalpyAdjointError::InvalidInput(
        "scan liquid-fraction derivative allocation refused"))?;
    for (vertex, &value) in h.iter().enumerate() {
        if vertex % 256 == 0 { poll(cx, vertex)?; }
        let curve = storage.phase_curve_at(vertex).ok_or(EnthalpyAdjointError::InvalidInput(
            "scan liquid-fraction derivative needs the original per-node chart"))?;
        out.push(derivative(curve, value));
    }
    poll(cx, h.len())?;
    Ok(out)
}

pub(super) fn pullback(cx: &Cx<'_>, slopes: &[FractionDerivative], h: &[f64], seed: &[f64])
    -> Result<Vec<f64>, ScanSensitivityError>
{
    poll(cx, 0)?;
    if slopes.len() != seed.len() || h.len() != seed.len() {
        return Err(EnthalpyAdjointError::InvalidInput("liquid-fraction seed must use original vertex order").into());
    }
    let mut out = Vec::new();
    out.try_reserve_exact(seed.len()).map_err(|_| EnthalpyAdjointError::InvalidInput(
        "scan liquid-fraction seed allocation refused"))?;
    for (vertex, (&weight, slope)) in seed.iter().zip(slopes).enumerate() {
        if vertex % 256 == 0 { poll(cx, vertex)?; }
        if !weight.is_finite() { return Err(EnthalpyAdjointError::NonFiniteArithmetic.into()); }
        let value = if weight == 0.0 { 0.0 } else {
            match slope {
                FractionDerivative::Smooth(slope) => weight * slope,
                FractionDerivative::Corner => return Err(EnthalpyAdjointError::ChartKink {
                    vertex, specific_enthalpy_j_kg: h[vertex],
                }.into()),
                FractionDerivative::Endpoint => return Err(EnthalpyAdjointError::ChartEndpoint {
                    vertex, specific_enthalpy_j_kg: h[vertex],
                }.into()),
                FractionDerivative::Unrepresentable => return Err(EnthalpyAdjointError::NonFiniteArithmetic.into()),
            }
        };
        if !value.is_finite() { return Err(EnthalpyAdjointError::NonFiniteArithmetic.into()); }
        out.push(value);
    }
    poll(cx, seed.len())?;
    Ok(out)
}
