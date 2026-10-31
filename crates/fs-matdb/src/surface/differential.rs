//! Source-cell tangents for the original sampled scalar surface evaluator.

use super::{coordinate_admitted, expanded_determinant, validate, vertices};
use crate::{EvaluationDecision, InterpolationPolicy, MatDbError, PropertyValue};

/// The canonical value and one original source cell's affine gradient.
///
/// Gradient coordinates have the same order and units as the source axes.
/// At a shared edge the ordinary evaluator selects the first containing cell.
/// At an exact sample the first incident source cell supplies the tangent while
/// the value/decision remain the original exact sample. These are explicit
/// one-sided derivatives, not a claim of differentiability at a source corner.
#[derive(Clone, Debug, PartialEq)]
pub struct SurfaceDifferential {
    pub value: f64,
    pub gradient: [f64; 2],
    pub triangle: u32,
    pub evaluation: EvaluationDecision,
}

/// Immutable numeric preparation of the SAME declared source cells.
///
/// This owns source coordinates/ordinates and cached analytical cell gradients.
/// It performs no source selection, validity-box extension, triangulation,
/// finite-difference query, or derivative uncertainty inference. The consumer
/// must retain its selected claim/receipt, admit LinearInside interpolation,
/// and check all actual applicability coordinates at each physical trial.
/// Values use the original exact-predicate/barycentric evaluator, including its
/// robust expansion fallback. No source query or identity hash occurs here.
#[derive(Clone, Debug)]
pub struct PreparedSurfaceResponse {
    axes: [String; 2],
    samples: Vec<([f64; 2], f64)>,
    triangles: Vec<[u32; 3]>,
    gradients: Vec<[f64; 2]>,
}

impl PropertyValue {
    /// Prepare analytical tangents on this payload's original surface cells.
    /// Source interpolation policy and physical applicability stay with the
    /// original claim owner. Preparation does not prove continuous box support.
    pub fn prepare_surface_response(&self) -> Result<PreparedSurfaceResponse, MatDbError> {
        let Self::Surface2D {
            axes,
            samples,
            triangles,
            ..
        } = self
        else {
            return Err(unsupported(
                "a differential surface requires a Surface2D payload",
            ));
        };
        validate(axes, samples, triangles)?;
        if triangles.is_empty() {
            return Err(unsupported(
                "a differential surface requires declared affine cells",
            ));
        }
        let gradients = triangles
            .iter()
            .map(|&triangle| {
                gradient(
                    vertices(samples, triangle),
                    triangle.map(|index| samples[index as usize].1),
                )
            })
            .collect::<Result<Vec<_>, _>>()?;
        Ok(PreparedSurfaceResponse {
            axes: axes.clone(),
            samples: samples.clone(),
            triangles: triangles.clone(),
            gradients,
        })
    }
}

impl PreparedSurfaceResponse {
    pub const fn axes(&self) -> &[String; 2] {
        &self.axes
    }
    pub fn samples(&self) -> &[([f64; 2], f64)] {
        &self.samples
    }
    pub fn triangles(&self) -> &[[u32; 3]] {
        &self.triangles
    }

    /// Evaluate the canonical source value with its selected cell tangent.
    /// Holes and isolated samples do not gain derivative support. Robust
    /// predicate/evaluation fallbacks retain their original allocation behavior.
    pub fn evaluate(&self, point: [f64; 2]) -> Result<SurfaceDifferential, MatDbError> {
        let (value, evaluation) = super::evaluate(
            &self.axes,
            &self.samples,
            &self.triangles,
            point,
            InterpolationPolicy::LinearInside,
        )?;
        let triangle = match &evaluation {
            EvaluationDecision::LinearInsideTriangle { triangle } => *triangle,
            EvaluationDecision::ExactSurfaceSample { sample } => self
                .triangles
                .iter()
                .position(|triangle| triangle.contains(sample))
                .map(|index| index as u32)
                .ok_or_else(|| unsupported("an isolated sample has no source-cell tangent"))?,
            _ => return Err(unsupported("unexpected canonical surface decision")),
        };
        Ok(SurfaceDifferential {
            value,
            gradient: self.gradients[triangle as usize],
            triangle,
            evaluation,
        })
    }
}

fn unsupported(reason: &'static str) -> MatDbError {
    MatDbError::UnsupportedEvaluation { reason }
}

fn gradient(points: [[f64; 2]; 3], values: [f64; 3]) -> Result<[f64; 2], MatDbError> {
    if values.iter().all(|value| *value == values[0]) {
        return Ok([0.0; 2]);
    }
    // Exact power-of-two normalization keeps the ordinate inside the same
    // expansion arithmetic domain as source coordinates. No ordinates change
    // in the value evaluator. Unrepresentable tangent scales refuse explicitly.
    let maximum = values.iter().map(|v| v.abs()).fold(0.0_f64, f64::max);
    let exponent = ((maximum.to_bits() >> 52) & 0x7ff).max(1);
    let scale = f64::from_bits(exponent << 52);
    let normalized = values.map(|value| value / scale);
    if normalized
        .iter()
        .zip(values)
        .any(|(&value, original)| !coordinate_admitted(value) || value == 0.0 && original != 0.0)
    {
        return Err(unsupported(
            "surface tangent ordinate scale is numerically unresolved",
        ));
    }
    let denominator = expanded_determinant(points[0], points[1], points[2]);
    let x = [0, 1, 2].map(|i| [normalized[i], points[i][1]]);
    let y = [0, 1, 2].map(|i| [points[i][0], normalized[i]]);
    let numerators = [
        expanded_determinant(x[0], x[1], x[2]),
        expanded_determinant(y[0], y[1], y[2]),
    ];
    let result = numerators.map(|numerator| (numerator / denominator) * scale);
    if !denominator.is_finite()
        || denominator == 0.0
        || result
            .iter()
            .zip(numerators)
            .any(|(&value, numerator)| !value.is_finite() || value == 0.0 && numerator != 0.0)
    {
        return Err(unsupported(
            "surface tangent is not representable as a finite scalar",
        ));
    }
    Ok(result)
}
