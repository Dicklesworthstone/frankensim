//! Source-cell tangents for the original sampled scalar surface evaluator.

use std::collections::BTreeMap;

use fs_ivl::{Sign, orient2d};

use super::{contains, coordinate_admitted, expanded_determinant, validate, vertices};
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
    boundary_edges: Vec<[u32; 2]>,
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
        // Source admission already rules out crossing/overlapping cells and
        // mismatched coincident vertices. An edge used once is exposed; cache
        // those original vertex identities rather than rebuild topology at
        // every nonlinear contact face.
        let mut edge_uses = BTreeMap::<[u32; 2], usize>::new();
        for triangle in triangles {
            for edge in 0..3 {
                let a = triangle[edge];
                let b = triangle[(edge + 1) % 3];
                *edge_uses.entry([a.min(b), a.max(b)]).or_default() += 1;
            }
        }
        Ok(PreparedSurfaceResponse {
            axes: axes.clone(),
            samples: samples.clone(),
            triangles: triangles.clone(),
            gradients,
            boundary_edges: edge_uses.into_iter()
                .filter_map(|(edge, uses)| (uses == 1).then_some(edge)).collect(),
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

    /// Classical derivative on the explicitly varying source coordinates.
    ///
    /// All incident source cells that contribute to the requested coordinate
    /// subspace must have the same projected gradient. Two varying coordinates
    /// require an open source neighbourhood; one varying coordinate requires
    /// both directions of its fixed-coordinate slice. Inactive components are
    /// zero. Thus an external-pressure edge can support a temperature tangent
    /// without pretending to support pressure variation. No cells, samples or
    /// interpolation values are added, and claim applicability remains with
    /// the caller. Work is linear in the admitted source-cell count.
    pub fn gradient_for_axes(
        &self, point: [f64; 2], varying_axes: [bool; 2],
    ) -> Result<[f64; 2], MatDbError> {
        if varying_axes == [false; 2] {
            super::evaluate(&self.axes, &self.samples, &self.triangles, point,
                InterpolationPolicy::LinearInside)?;
            return Ok([0.0; 2]);
        }
        self.evaluate(point)?;
        let full = varying_axes.iter().all(|&varying| varying);
        if full && self.boundary_edges.iter().any(|&[a, b]| {
            on_edge(self.samples[a as usize].0, self.samples[b as usize].0, point)
        }) {
            return Err(unsupported("surface derivative lacks an open two-coordinate neighbourhood"));
        }
        let varying = usize::from(!varying_axes[0]);
        let mut supported = [false; 2];
        let mut answer: Option<[f64; 2]> = None;
        for (&triangle, &gradient) in self.triangles.iter().zip(&self.gradients) {
            let positions = vertices(&self.samples, triangle);
            if !contains(positions, point) {
                continue;
            }
            let directions = if full { [true; 2] }
                else { slice_directions(positions, point, varying) };
            // A wedge touching a fixed-coordinate line only at this vertex
            // cannot change the derivative of that line.
            if !directions.iter().any(|&direction| direction) {
                continue;
            }
            for side in 0..2 { supported[side] |= directions[side]; }
            let projected = [0, 1].map(|axis| if varying_axes[axis] { gradient[axis] } else { 0.0 });
            if answer.is_some_and(|previous| previous != projected) {
                return Err(unsupported("surface has unequal incident gradients in a varying coordinate"));
            }
            answer = Some(projected);
        }
        if !supported.iter().all(|&side| side) {
            return Err(unsupported("surface derivative lacks a two-sided fixed-coordinate neighbourhood"));
        }
        answer.ok_or_else(|| unsupported("surface point has no contributing derivative cell"))
    }
}

fn on_edge(a: [f64; 2], b: [f64; 2], point: [f64; 2]) -> bool {
    (0..2).all(|axis| a[axis].min(b[axis]) <= point[axis]
        && point[axis] <= a[axis].max(b[axis]))
        && orient2d(a, b, point) == Sign::Zero
}

// Exact local tangent cone of a containing triangle along one coordinate.
// Only edges through the actual point constrain an infinitesimal direction.
// Their directional orientation signs are simple transverse-coordinate
// comparisons, avoiding a rounded or artificially displaced query point.
fn slice_directions(triangle: [[f64; 2]; 3], point: [f64; 2], varying: usize) -> [bool; 2] {
    let winding = orient2d(triangle[0], triangle[1], triangle[2]);
    let mut supported = [true; 2]; // negative, positive coordinate directions
    for edge in 0..3 {
        let a = triangle[edge];
        let b = triangle[(edge + 1) % 3];
        if orient2d(a, b, point) != Sign::Zero { continue; }
        let [first, second] = if varying == 0 { [a[1], b[1]] } else { [b[0], a[0]] };
        let sign = if first < second { Sign::Negative }
            else if first > second { Sign::Positive } else { Sign::Zero };
        supported[0] &= sign == Sign::Zero || sign.flip() == winding;
        supported[1] &= sign == Sign::Zero || sign == winding;
    }
    supported
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
