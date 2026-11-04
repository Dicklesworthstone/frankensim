//! Signed acoustic power through original, explicitly sided element faces.
//!
//! A shared material interface can be measured from either tetrahedron. Each
//! side keeps its own density and FE gradient; their discrepancy is observable,
//! not corrected into artificial continuity. Reactive power is not dissipation.

use super::super::super::tri_quadrature;
use super::super::time::HeterogeneousTimeState;
use super::super::{
    AcousticFeError, HeterogeneousError, HeterogeneousHelmholtzFe, MAX_QUADRATURE_ORDER,
    Medium, finite, invalid, poll,
};
use fs_math::{c64::C64, dd::Dd};
use std::collections::BTreeSet;

/// One original mesh face, with its explicitly selected one-sided volume trace.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct ElementFace {
    /// Original tetrahedron supplying the material, pressure gradient and normal.
    pub element: usize,
    /// Three distinct vertices of that tetrahedron; winding does not set the sign.
    pub vertices: [u32; 3],
}

/// Explicit work/storage admission for a reusable surface observation.
#[derive(Clone, Copy, Debug)]
pub struct SurfaceProbeLimits {
    /// Maximum distinct physical triangles.
    pub max_faces: usize,
    /// Maximum retained quadrature points. Each stores one local FE stencil;
    /// storage grows with the polynomial order as well as this point count.
    pub max_quadrature_points: usize,
}

impl Default for SurfaceProbeLimits {
    fn default() -> Self {
        Self {
            max_faces: 4096,
            max_quadrature_points: 65_536,
        }
    }
}

/// Integral of `p conj(v_out)/2` in the original `exp(+i omega t)` convention.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct HarmonicSurfacePower {
    /// Signed NET transported power [W], positive out of the selected elements.
    /// Coherent incident/reflected contributions are combined before integration.
    pub active_power_w: f64,
    /// Signed reactive power [W]; retained separately from transported power/heat.
    pub reactive_power_w: f64,
}

struct NormalShape {
    dof: usize,
    value: f64,
    normal_derivative_m_inv: f64,
}

struct SurfaceSample {
    shapes: Vec<NormalShape>,
    inverse_density: f64,
    weight_m2: f64,
}

/// Prepared high-order power measurement on a selected set of original faces.
///
/// Preparation reuses the original point-probe stencils and material samples,
/// then retains only pressure and NORMAL-derivative weights. Repeated harmonic
/// or transient fields need no allocations, geometry rebuilding or source queries.
/// Integrals use fixed-order double-double accumulation for cancelling fluxes.
/// Finite-element gradient error and smooth-material quadrature error still need
/// independent convergence; no exact coarse-mesh conservation is asserted.
pub struct HeterogeneousSurfaceProbe<'m, 'f, 'c> {
    model: &'m HeterogeneousHelmholtzFe<'f, 'c>,
    faces: Vec<ElementFace>,
    samples: Vec<SurfaceSample>,
    area_m2: f64,
    quadrature_order: usize,
}

impl<'f, 'c> HeterogeneousHelmholtzFe<'f, 'c> {
    /// Prepare a surface-power measurement from the ORIGINAL frozen source.
    ///
    /// Normals point outward from each explicitly selected tetrahedron, even
    /// for internal faces. A physical face may appear only once; measure opposite
    /// sides separately instead of silently cancelling duplicate faces.
    ///
    /// All faces and budgets are checked before any source query. The original
    /// triangle Duffy rule uses `quadrature_order` points per axis, between
    /// `geometry.order()+1` and `MAX_QUADRATURE_ORDER`. This integrates the
    /// pressure-times-gradient polynomial for constant media; smoothly varying
    /// density requires its own quadrature study. Source errors propagate intact.
    /// As in `new`, source-generic preparation reports cancellation and invalid
    /// limits through `AcousticFeError::InvalidInput`; no partial probe is returned.
    pub fn surface_probe<'m, E: From<AcousticFeError>>(
        &'m self,
        faces: &[ElementFace],
        quadrature_order: usize,
        limits: SurfaceProbeLimits,
        material_at: &mut impl FnMut(usize, [f64; 3]) -> Result<Medium, E>,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<HeterogeneousSurfaceProbe<'m, 'f, 'c>, E> {
        let checkpoint = |cancelled: &mut dyn FnMut() -> bool| {
            if cancelled() {
                Err(invalid("acoustic surface preparation cancelled"))
            } else {
                Ok(())
            }
        };
        checkpoint(cancelled)?;
        if faces.is_empty()
            || limits.max_faces == 0
            || faces.len() > limits.max_faces
            || limits.max_quadrature_points == 0
            || quadrature_order < self.geometry.order + 1
            || quadrature_order > MAX_QUADRATURE_ORDER
        {
            return Err(invalid("acoustic surface faces, quadrature or limits outside admitted range").into());
        }
        let point_count = faces.len().checked_mul(quadrature_order * quadrature_order)
            .ok_or_else(|| invalid("acoustic surface quadrature count overflow"))?;
        if point_count > limits.max_quadrature_points {
            return Err(invalid("acoustic surface exceeds its quadrature-point budget").into());
        }
        let mut seen = BTreeSet::new();
        let mut admitted = Vec::with_capacity(faces.len());
        let mut area = Dd::ZERO;
        for &face in faces {
            checkpoint(cancelled)?;
            let mut key = face.vertices;
            key.sort_unstable();
            if !seen.insert(key) {
                return Err(invalid("acoustic surface physical faces must be unique").into());
            }
            let (_, normal, face_area) = self.geometry.element_face_geometry(face.element, face.vertices)?;
            let tet = self.geometry.tetrahedra()[face.element];
            let mut local = [0; 3];
            for (i, vertex) in face.vertices.iter().enumerate() {
                local[i] = tet.iter().position(|v| v == vertex)
                    .ok_or_else(|| invalid("surface face is not owned by the selected tetrahedron"))?;
            }
            area = area + Dd::from_f64(face_area);
            admitted.push((face, normal, face_area, local));
        }
        let area_m2 = area.to_f64();
        if !area_m2.is_finite() || area_m2 <= 0.0 {
            return Err(invalid("unrepresentable acoustic surface area").into());
        }
        let quadrature = tri_quadrature(quadrature_order);
        let mut samples = Vec::with_capacity(point_count);
        for (face, normal, face_area, local) in admitted {
            for &(barycentric, weight) in &quadrature {
                checkpoint(cancelled)?;
                let mut point = [0.0; 4];
                for i in 0..3 {
                    point[local[i]] = barycentric[i];
                }
                let probe = self.probe(face.element, point, material_at)?;
                let weight_m2 = (2.0 * weight) * face_area;
                if !weight_m2.is_finite() || weight_m2 <= 0.0 {
                    return Err(invalid("unrepresentable acoustic surface quadrature weight").into());
                }
                let mut shapes = Vec::with_capacity(probe.shapes.len());
                // Keep zero pressure traces too: their normal derivatives can
                // contribute to velocity on the face (including interior bubbles).
                for shape in probe.shapes {
                    let normal_derivative_m_inv = shape.gradient_m_inv.iter().zip(normal)
                        .map(|(gradient, n)| gradient * n).sum::<f64>();
                    if !normal_derivative_m_inv.is_finite() {
                        return Err(invalid("unrepresentable acoustic surface derivative").into());
                    }
                    shapes.push(NormalShape {
                        dof: shape.dof,
                        value: shape.value,
                        normal_derivative_m_inv,
                    });
                }
                samples.push(SurfaceSample {
                    shapes,
                    inverse_density: probe.inverse_density,
                    weight_m2,
                });
            }
        }
        checkpoint(cancelled)?;
        Ok(HeterogeneousSurfaceProbe {
            model: self,
            faces: faces.to_vec(),
            samples,
            area_m2,
            quadrature_order,
        })
    }
}

impl HeterogeneousSurfaceProbe<'_, '_, '_> {
    /// Original face/element selections, in their retained integration order.
    #[must_use]
    pub fn faces(&self) -> &[ElementFace] {
        &self.faces
    }

    /// Total physical surface area [m^2], without signed cancellation.
    #[must_use]
    pub const fn area_m2(&self) -> f64 {
        self.area_m2
    }

    /// Original triangle Gauss-point count on each Duffy axis.
    #[must_use]
    pub const fn quadrature_order(&self) -> usize {
        self.quadrature_order
    }

    /// Number of cached local pressure/normal-gradient stencils.
    #[must_use]
    pub fn sample_count(&self) -> usize {
        self.samples.len()
    }

    /// Integrate coherent pressure and normal velocity at their physical
    /// frequency. Use the same original coefficients/frequency as the solve.
    /// A bare harmonic slice carries no independent model identity.
    pub fn harmonic(
        &self,
        pressure_pa: &[C64],
        frequency_hz: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<HarmonicSurfacePower, HeterogeneousError> {
        poll(cancelled)?;
        let omega = core::f64::consts::TAU * frequency_hz;
        if pressure_pa.len() != self.model.ndof() || !omega.is_finite() || omega <= 0.0 {
            return Err(invalid("surface power requires original coefficients and a finite positive frequency").into());
        }
        let mut active = Dd::ZERO;
        let mut reactive = Dd::ZERO;
        for sample in &self.samples {
            poll(cancelled)?;
            let factor = sample.inverse_density / omega;
            if !factor.is_finite() || factor <= 0.0 {
                return Err(invalid("unrepresentable acoustic surface velocity factor").into());
            }
            let mut pressure = C64::ZERO;
            let mut derivative = C64::ZERO;
            for shape in &sample.shapes {
                let coefficient = pressure_pa[shape.dof];
                pressure = pressure + coefficient.scale(shape.value);
                derivative = derivative + coefficient.scale(shape.normal_derivative_m_inv);
            }
            let velocity = C64::new(-derivative.im, derivative.re).scale(factor);
            let intensity = pressure.scale(0.5) * velocity.conj();
            if !finite(pressure) || !finite(derivative) || !finite(velocity) || !finite(intensity) {
                return Err(invalid("nonfinite harmonic acoustic surface field").into());
            }
            active = active + Dd::from_f64(intensity.re) * Dd::from_f64(sample.weight_m2);
            reactive = reactive + Dd::from_f64(intensity.im) * Dd::from_f64(sample.weight_m2);
        }
        poll(cancelled)?;
        let power = HarmonicSurfacePower {
            active_power_w: active.to_f64(),
            reactive_power_w: reactive.to_f64(),
        };
        if !power.active_power_w.is_finite() || !power.reactive_power_w.is_finite() {
            return Err(invalid("acoustic surface power overflow").into());
        }
        Ok(power)
    }

    /// Instantaneous signed outward power [W], `integral(p v_out)`, on the
    /// exact original transient operator. There is no harmonic one-half factor.
    pub fn transient(
        &self,
        state: &HeterogeneousTimeState<'_, '_, '_>,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<f64, HeterogeneousError> {
        poll(cancelled)?;
        if !core::ptr::eq(self.model, state.model) {
            return Err(invalid("surface state belongs to another spatial/material operator").into());
        }
        let mut power = Dd::ZERO;
        for sample in &self.samples {
            poll(cancelled)?;
            let mut pressure = 0.0;
            let mut derivative = 0.0;
            for shape in &sample.shapes {
                pressure += state.pressure_pa()[shape.dof] * shape.value;
                derivative += state.pressure_potential_pa_s()[shape.dof] * shape.normal_derivative_m_inv;
            }
            let velocity = -derivative * sample.inverse_density;
            if !pressure.is_finite() || !derivative.is_finite() || !velocity.is_finite() {
                return Err(invalid("nonfinite transient acoustic surface field").into());
            }
            power = power + Dd::from_f64(pressure) * Dd::from_f64(velocity)
                * Dd::from_f64(sample.weight_m2);
        }
        poll(cancelled)?;
        let power = power.to_f64();
        if !power.is_finite() {
            return Err(invalid("transient acoustic surface power overflow").into());
        }
        Ok(power)
    }
}

#[cfg(test)]
mod tests;
