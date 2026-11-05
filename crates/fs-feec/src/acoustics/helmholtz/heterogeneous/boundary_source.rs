//! Prepared quadratic normal motion on the original physical boundary faces.
//!
//! A dimensionless six-node shape `g` driven by an inward normal velocity
//! amplitude `u` [m/s] contributes `b_i = integral(phi_i g) u` [m^3/s].
//! The pressure equation receives `i omega b`; its conjugate work is
//! `integral(p g) u`. Shape preparation retains the original H1 pressure basis,
//! including edge and face tests, without replacing midside motion by corner
//! values, dividing by area, or normalizing the signed shape.
//!
//! The volume owner's admitted affine or curved map supplies the surface
//! measure. A quadratic velocity trace does not alter that geometry or infer
//! a structural displacement, material law, or radiation model.

use super::super::tri_quadrature;
use super::{HeterogeneousError, HeterogeneousHelmholtzFe, MAX_QUADRATURE_ORDER, invalid, poll};
use fs_math::{c64::C64, dd::Dd};
use std::collections::BTreeSet;

/// Bound on simultaneously driven prepared boundary-source columns. Exact
/// trajectory preparations additionally retain their own augmented-state cap.
pub const MAX_BOUNDARY_SOURCES: usize = 128;

pub(super) fn admit_count(count: usize) -> Result<(), HeterogeneousError> {
    if count > MAX_BOUNDARY_SOURCES {
        return Err(HeterogeneousError::Limit {
            what: "simultaneous prepared boundary acoustic sources",
            requested: count,
            maximum: MAX_BOUNDARY_SOURCES,
        });
    }
    Ok(())
}

/// Frozen quadratic inward-normal velocity shape on one original fluid model.
///
/// Values are dimensionless. A velocity amplitude [m/s] multiplies the stored
/// load [m^2], giving the original pressure-space source [m^3/s]. Preparation
/// integrates the signed shape as supplied, with no area normalization. The
/// model borrow prevents reusing a column with another material/geometry owner.
/// Repeated drives reuse this projection and its work-conjugate observation.
pub struct BoundaryVelocitySource<'m, 'f, 'c> {
    pub(super) model: &'m HeterogeneousHelmholtzFe<'f, 'c>,
    pub(super) load: Vec<f64>,
    pub(super) faces: Vec<[u32; 3]>,
    integral_m2: f64,
}

impl BoundaryVelocitySource<'_, '_, '_> {
    /// Canonically sorted, unique original boundary faces supporting the shape.
    /// Faces with identically zero supplied values remain present so physical
    /// boundary-law admission can inspect the complete declared support.
    #[must_use]
    pub fn support_faces(&self) -> &[[u32; 3]] {
        &self.faces
    }

    /// Whether the original pressure-space projection has any nonzero entry.
    /// A signed shape can have zero net area integral and still couple to the
    /// pressure field. This checks the cached load itself, without dividing by
    /// area or normalizing the shape.
    #[must_use]
    pub fn has_nonzero_projection(&self) -> bool {
        self.load.iter().any(|&weight| weight != 0.0)
    }

    /// Signed area integral of the dimensionless shape [m^2]. This can be zero
    /// for a nonzero source; neither projection nor work uses it as a divisor.
    #[must_use]
    pub const fn integrated_shape_m2(&self) -> f64 {
        self.integral_m2
    }

    pub(super) fn validate(
        &self,
        model: &HeterogeneousHelmholtzFe<'_, '_>,
    ) -> Result<(), HeterogeneousError> {
        if !core::ptr::eq(self.model, model) {
            return Err(
                invalid("boundary source belongs to another spatial/material operator").into(),
            );
        }
        Ok(())
    }

    pub(super) fn harmonic_pressure(&self, pressure: &[C64]) -> C64 {
        let mut real = Dd::from_f64(0.0);
        let mut imaginary = Dd::from_f64(0.0);
        for (&p, &weight) in pressure.iter().zip(&self.load) {
            real = real + Dd::from_f64(p.re) * Dd::from_f64(weight);
            imaginary = imaginary + Dd::from_f64(p.im) * Dd::from_f64(weight);
        }
        C64::new(real.to_f64(), imaginary.to_f64())
    }

    pub(super) fn midpoint_pressure(&self, initial: &[f64], next: &[f64]) -> f64 {
        let mut value = Dd::from_f64(0.0);
        for ((&a, &b), &weight) in initial.iter().zip(next).zip(&self.load) {
            value = value + Dd::from_f64(f64::midpoint(a, b)) * Dd::from_f64(weight);
        }
        value.to_f64()
    }
}

impl<'f, 'c> HeterogeneousHelmholtzFe<'f, 'c> {
    /// Project quadratic normal motion onto every original pressure trace test.
    ///
    /// Each entry supplies an original boundary triangle and six dimensionless
    /// values in `corner0, corner1, corner2, edge01, edge12, edge20` order. Edge
    /// values are the values AT the reference-edge midpoints mapped into the
    /// physical face, not hierarchical coefficients. Positive amplitude means
    /// inward velocity along the local normal; winding does
    /// not change this scalar convention. All values follow the supplied face
    /// order, including when that order is reversed or cyclically permuted.
    ///
    /// The original face rule uses at least `pressure_order + 2` Gauss points
    /// per Duffy axis, and at least the curved geometry's declared order. This
    /// integrates the quadratic shape against the full pressure trace and the
    /// pointwise physical measure. The curved metric needs its own quadrature
    /// convergence; use [`Self::quadratic_velocity_source_with_quadrature`] to
    /// vary that axis. No material, face-average or corner-only approximation
    /// is used.
    /// Faces must be nonempty, unique and on the real exterior boundary. Zero
    /// and signed shapes are allowed. Separate prepared sources may overlap
    /// and interfere; each subsequent solver admits its pressure-release and
    /// other physical boundary-law conflicts before using any amplitude.
    /// Cancellation, invalid support or nonfinite projection publishes no
    /// partially prepared source.
    pub fn quadratic_velocity_source<'m>(
        &'m self,
        traces: &[([u32; 3], [f64; 6])],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<BoundaryVelocitySource<'m, 'f, 'c>, HeterogeneousError> {
        self.quadratic_velocity_source_with_quadrature(
            traces,
            self.geometry.face_quadrature_order(self.geometry.order + 2),
            cancelled,
        )
    }

    /// Prepare the same physical source with an independently selected face
    /// quadrature order, between `pressure_order + 2` and 16 Gauss points per
    /// Duffy axis. This varies source/curved-metric integration without changing
    /// the original pressure basis, geometry, material or source amplitudes.
    pub fn quadratic_velocity_source_with_quadrature<'m>(
        &'m self,
        traces: &[([u32; 3], [f64; 6])],
        quadrature_order: usize,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<BoundaryVelocitySource<'m, 'f, 'c>, HeterogeneousError> {
        poll(cancelled)?;
        if traces.is_empty() {
            return Err(invalid("nonempty quadratic boundary source required").into());
        }
        if quadrature_order < self.geometry.order + 2 || quadrature_order > MAX_QUADRATURE_ORDER {
            return Err(invalid("quadratic source quadrature order outside admitted range").into());
        }
        let quadrature = tri_quadrature(quadrature_order);
        let mut load = vec![Dd::from_f64(0.0); self.ndof()];
        let mut integral = Dd::from_f64(0.0);
        let mut faces = BTreeSet::new();
        for &(face, values) in traces {
            poll(cancelled)?;
            let mut key = face;
            key.sort_unstable();
            if !faces.insert(key) || values.iter().any(|value| !value.is_finite()) {
                return Err(invalid(
                    "quadratic boundary source requires unique faces and finite values",
                )
                .into());
            }
            let rule = self.geometry.face_rule(face, &quadrature)?;
            if !rule.area.is_finite() || rule.area <= 0.0 {
                return Err(invalid("unrepresentable quadratic source face geometry").into());
            }
            let mut input_to_local = [0; 3];
            for (input, vertex) in face.iter().enumerate() {
                input_to_local[input] = rule
                    .vertices
                    .iter()
                    .position(|v| v == vertex)
                    .ok_or_else(|| invalid("quadratic source vertex mismatch"))?;
            }
            for (q, &(local_barycentric, _)) in quadrature.iter().enumerate() {
                poll(cancelled)?;
                let [a, b, c] = input_to_local.map(|local| local_barycentric[local]);
                let shapes = [
                    a * (2.0 * a - 1.0),
                    b * (2.0 * b - 1.0),
                    c * (2.0 * c - 1.0),
                    4.0 * a * b,
                    4.0 * b * c,
                    4.0 * c * a,
                ];
                let mut value = Dd::from_f64(0.0);
                for (&nodal, shape) in values.iter().zip(shapes) {
                    value = value + Dd::from_f64(nodal) * Dd::from_f64(shape);
                }
                let weighted = value * Dd::from_f64(rule.weights[q]);
                if !value.to_f64().is_finite() || !weighted.to_f64().is_finite() {
                    return Err(invalid("quadratic boundary source projection overflow").into());
                }
                integral = integral + weighted;
                for (local, &global) in rule.dofs.iter().enumerate() {
                    load[global] = load[global] + weighted * Dd::from_f64(rule.phi[local][q]);
                }
            }
        }
        let load: Vec<_> = load.into_iter().map(Dd::to_f64).collect();
        let integral_m2 = integral.to_f64();
        if !integral_m2.is_finite() || load.iter().any(|value| !value.is_finite()) {
            return Err(invalid("quadratic boundary source projection overflow").into());
        }
        poll(cancelled)?;
        Ok(BoundaryVelocitySource {
            model: self,
            load,
            faces: faces.into_iter().collect(),
            integral_m2,
        })
    }
}

#[cfg(test)]
mod tests;
