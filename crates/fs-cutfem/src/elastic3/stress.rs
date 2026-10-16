//! Bulk stress from the retained Q1 field and the assembly's reference tensor.
//!
//! Voigt entries are `[xx, yy, zz, xy, yz, xz]`, with physical shear stresses
//! (no engineering-strain factor in the output). These are numerical samples,
//! not recovered nodal stresses or bounds on a continuum maximum. Ghost and
//! Nitsche stabilization affect the solved field but are not material stress.
use super::*;

/// One original bulk quadrature point; point order is cell order followed by
/// that cell's retained rule order. No geometry is regenerated or reclassified.
#[derive(Debug, Clone, PartialEq)]
pub struct BulkStressPoint3 {
    /// Index in the operator's `scales()` and `volumes()` arrays.
    pub cell: usize,
    pub position: [f64; 3],
    /// Positive reference-volume quadrature weight, not normalized.
    pub weight: f64,
    /// `C_reference : eps(u)` for the original admitted material.
    pub reference_stress: [f64; 6],
    /// Current physical stress: `scale[cell] * reference_stress`.
    pub stress: [f64; 6],
}

/// Partial derivatives of a scalar functional of the physical stress samples.
/// These do not include the effect of equilibrium changing with material.
#[derive(Debug, Clone, PartialEq)]
pub struct BulkStressPullback3 {
    /// `dJ/du` in the operator's independent displacement coordinates.
    pub displacement: Vec<f64>,
    /// `dJ/dscale` at FIXED u, one value per retained active cell.
    pub scales: Vec<f64>,
}

fn poll(checkpoint: &mut impl FnMut() -> ControlFlow<()>) -> Result<(), ElasticityError3> {
    if checkpoint().is_break() {
        Err(ElasticityError3::Cancelled)
    } else {
        Ok(())
    }
}

/// Stress caused by one Cartesian displacement component with gradient g.
fn basis(lame: [f64; 2], g: [f64; 3], component: usize) -> [f64; 6] {
    let [lambda, mu] = lame;
    let mut s = [lambda * g[component]; 6];
    s[component] += 2.0 * mu * g[component];
    s[3] = mu
        * match component {
            0 => g[1],
            1 => g[0],
            _ => 0.0,
        };
    s[4] = mu
        * match component {
            1 => g[2],
            2 => g[1],
            _ => 0.0,
        };
    s[5] = mu
        * match component {
            0 => g[2],
            2 => g[0],
            _ => 0.0,
        };
    s
}

fn stress(lame: [f64; 2], gradients: &[[f64; 3]; 8], local: &[f64; 24]) -> [f64; 6] {
    let mut value = [0.0; 6];
    for i in 0..24 {
        let b = basis(lame, gradients[i / 3], i % 3);
        for k in 0..6 {
            value[k] += b[k] * local[i];
        }
    }
    value
}

impl CutElasticity3 {
    /// Exact `z^T (dK/dscale_c) u`, including each cell's original Nitsche
    /// block and half of each incident ghost term. No current scale factor is
    /// applied. Direct contraction avoids cancellation in polarization when
    /// primal and adjoint magnitudes differ substantially.
    pub fn scale_bilinear_forms(
        &self,
        z: &[f64],
        u: &[f64],
        mut checkpoint: impl FnMut() -> ControlFlow<()>,
    ) -> Result<Vec<f64>, ElasticityError3> {
        self.admit_stress_state(Some(u), &mut checkpoint)?;
        self.admit_stress_state(Some(z), &mut checkpoint)?;
        let mut result = vec![0.0; self.cells.len()];
        for (id, cell) in self.cells.iter().enumerate() {
            poll(&mut checkpoint)?;
            let local = self.stress_local(cell, u);
            let dual = self.stress_local(cell, z);
            for i in 0..24 {
                let applied: f64 = cell.stiffness[i]
                    .iter()
                    .zip(&local)
                    .map(|(k, u)| k * u)
                    .sum();
                result[id] += dual[i] * applied;
            }
        }
        for face in &self.ghosts {
            poll(&mut checkpoint)?;
            let mut value = 0.0;
            for c in 0..3 {
                let jump = |x: &[f64]| -> f64 {
                    face.nodes
                        .iter()
                        .zip(&face.jump)
                        .map(|(&n, &j)| if self.fixed[n] { 0.0 } else { j * x[3 * n + c] })
                        .sum()
                };
                value += face.weight * jump(z) * jump(u);
            }
            for &cell in &face.cells {
                result[cell] += 0.5 * value;
            }
        }
        if !result.iter().all(|v| v.is_finite()) {
            return Err(ElasticityError3::Invalid(
                "bilinear scale contraction overflow",
            ));
        }
        poll(&mut checkpoint)?;
        Ok(result)
    }

    fn stress_point_count(
        &self,
        checkpoint: &mut impl FnMut() -> ControlFlow<()>,
    ) -> Result<usize, ElasticityError3> {
        let mut count = 0usize;
        for cell in &self.cells {
            poll(checkpoint)?;
            count = count
                .checked_add(cell.rules.bulk().len())
                .ok_or(ElasticityError3::Invalid(
                    "bulk stress point count overflow",
                ))?;
        }
        Ok(count)
    }

    fn admit_stress_state(
        &self,
        u: Option<&[f64]>,
        checkpoint: &mut impl FnMut() -> ControlFlow<()>,
    ) -> Result<(), ElasticityError3> {
        poll(checkpoint)?;
        if u.is_some_and(|v| v.len() != self.n()) {
            return Err(ElasticityError3::Invalid(
                "invalid bulk stress displacement shape",
            ));
        }
        if let Some(u) = u {
            for chunk in u.chunks(256) {
                poll(checkpoint)?;
                if chunk.iter().any(|v| !v.is_finite()) {
                    return Err(ElasticityError3::Invalid(
                        "nonfinite bulk stress displacement",
                    ));
                }
            }
        }
        for chunk in self.scales.chunks(256) {
            poll(checkpoint)?;
            if chunk
                .iter()
                .any(|v| !v.is_finite() || *v <= 0.0 || *v > 1.0)
            {
                return Err(ElasticityError3::Invalid(
                    "invalid bulk stress material scale",
                ));
            }
        }
        Ok(())
    }

    fn stress_local(&self, cell: &Cell3, u: &[f64]) -> [f64; 24] {
        std::array::from_fn(|i| {
            let node = cell.nodes[i / 3];
            if self.fixed[node] {
                0.0
            } else {
                u[3 * node + i % 3]
            }
        })
    }

    /// Evaluate the actual bulk quadrature stress with the assembly's Lamé
    /// tensor and current cell scales. Clamped coefficients are masked exactly
    /// as in the operator. `max_points` admits the output allocation before any
    /// point is evaluated. Every point and publication boundary polls the caller.
    /// No state changes on success, refusal or cancellation.
    pub fn bulk_stress(
        &self,
        u: &[f64],
        max_points: usize,
        mut checkpoint: impl FnMut() -> ControlFlow<()>,
    ) -> Result<Vec<BulkStressPoint3>, ElasticityError3> {
        self.admit_stress_state(Some(u), &mut checkpoint)?;
        let count = self.stress_point_count(&mut checkpoint)?;
        if count > max_points {
            return Err(ElasticityError3::Invalid(
                "bulk stress point allowance exhausted",
            ));
        }
        let mut points = Vec::with_capacity(count);
        for (id, cell) in self.cells.iter().enumerate() {
            let local = self.stress_local(cell, u);
            for &(position, weight) in cell.rules.bulk() {
                poll(&mut checkpoint)?;
                let (_, gradients) = q1(cell.bounds, position);
                let reference_stress = stress(self.lame, &gradients, &local);
                let stress = reference_stress.map(|s| self.scales[id] * s);
                if !reference_stress
                    .iter()
                    .chain(&stress)
                    .all(|s| s.is_finite())
                {
                    return Err(ElasticityError3::Invalid("bulk stress overflow"));
                }
                points.push(BulkStressPoint3 {
                    cell: id,
                    position,
                    weight,
                    reference_stress,
                    stress,
                });
            }
        }
        poll(&mut checkpoint)?;
        Ok(points)
    }

    /// Exact VJP of the PHYSICAL stresses returned by `bulk_stress`.
    /// `derivatives[i]` is dJ/dsigma_i in the exact retained point order;
    /// quadrature weights are NOT inserted here. Include them in the supplied
    /// cotangents when differentiating a volume integral. The direct scale
    /// derivative holds u fixed. An equilibrium adjoint must separately include
    /// the operator's complete scale derivative, including stabilization.
    pub fn bulk_stress_pullback(
        &self,
        u: &[f64],
        derivatives: &[[f64; 6]],
        checkpoint: impl FnMut() -> ControlFlow<()>,
    ) -> Result<BulkStressPullback3, ElasticityError3> {
        self.stress_pullback(Some(u), derivatives, checkpoint)
    }

    /// Exact VJP of `reference_stress`, without the current stiffness scale.
    /// This is the linear map `(C_reference B)^T`, useful for explicitly
    /// declared stress-relaxation models. It has no direct scale derivative.
    /// Point ordering and unweighted cotangent conventions match the physical
    /// pullback. The reference material is the admitted card, not unit Young's
    /// modulus. Neither entry point differentiates geometry or quadrature.
    pub fn reference_bulk_stress_pullback(
        &self,
        derivatives: &[[f64; 6]],
        checkpoint: impl FnMut() -> ControlFlow<()>,
    ) -> Result<Vec<f64>, ElasticityError3> {
        Ok(self
            .stress_pullback(None, derivatives, checkpoint)?
            .displacement)
    }

    fn stress_pullback(
        &self,
        u: Option<&[f64]>,
        derivatives: &[[f64; 6]],
        mut checkpoint: impl FnMut() -> ControlFlow<()>,
    ) -> Result<BulkStressPullback3, ElasticityError3> {
        self.admit_stress_state(u, &mut checkpoint)?;
        if derivatives.len() != self.stress_point_count(&mut checkpoint)? {
            return Err(ElasticityError3::Invalid(
                "one stress cotangent per bulk point required",
            ));
        }
        for d in derivatives {
            poll(&mut checkpoint)?;
            if !d.iter().all(|v| v.is_finite()) {
                return Err(ElasticityError3::Invalid("nonfinite bulk stress cotangent"));
            }
        }
        let mut displacement = vec![0.0; self.n()];
        let mut scales = if u.is_some() {
            vec![0.0; self.cells.len()]
        } else {
            Vec::new()
        };
        let mut index = 0;
        for (id, cell) in self.cells.iter().enumerate() {
            let local = u.map(|u| self.stress_local(cell, u));
            let factor = if u.is_some() { self.scales[id] } else { 1.0 };
            for &(position, _) in cell.rules.bulk() {
                poll(&mut checkpoint)?;
                let (_, gradients) = q1(cell.bounds, position);
                let d = &derivatives[index];
                index += 1;
                if let Some(local) = &local {
                    let value = stress(self.lame, &gradients, local);
                    scales[id] += d.iter().zip(value).map(|(a, b)| a * b).sum::<f64>();
                }
                for i in 0..24 {
                    let node = cell.nodes[i / 3];
                    if self.fixed[node] {
                        continue;
                    }
                    let b = basis(self.lame, gradients[i / 3], i % 3);
                    displacement[3 * node + i % 3] +=
                        factor * d.iter().zip(b).map(|(a, b)| a * b).sum::<f64>();
                }
            }
        }
        if !displacement.iter().chain(&scales).all(|v| v.is_finite()) {
            return Err(ElasticityError3::Invalid("bulk stress pullback overflow"));
        }
        poll(&mut checkpoint)?;
        Ok(BulkStressPullback3 {
            displacement,
            scales,
        })
    }
}
