//! Spatially varying passive walls on the original pressure traces.
//!
//! The outward velocity law is v = Y(x, omega) p, with exp(+i omega t).
//! Both real and imaginary parts are integrated against the FULL original
//! pressure trace: averaging Y over a face loses its correlation with p.
//! Positive quadrature weights and pointwise Re(Y) >= 0 give a passive
//! discrete conductance form. This does not certify the unsampled continuous
//! material law, its frequency dependence, or a causal time-domain realization.

use super::super::geometry::GeometryFacePoint;
use super::super::tri_quadrature;
use super::{
    AcousticFeError, HeterogeneousError, HeterogeneousHelmholtzFe, MAX_QUADRATURE_ORDER, finite,
    invalid,
};
use fs_math::{c64::C64, dd::Dd};
use fs_sparse::{Coo, Csr};
use std::collections::BTreeSet;

pub(super) mod relaxation;

/// A sampled specific admittance [m/(Pa s)] on original boundary faces.
///
/// The exact spatial/material owner AND frequency are retained. Subsequent
/// harmonic right-hand sides reuse these forms without resampling a material.
/// This is a locally reacting, prescribed wall law; it does not supply
/// exterior radiation, nonlocal wall coupling or a transient memory model.
pub struct SpatialBoundaryAdmittance<'m, 'f, 'c> {
    model: &'m HeterogeneousHelmholtzFe<'f, 'c>,
    frequency_hz: f64,
    quadrature_order: usize,
    pub(super) faces: Vec<[u32; 3]>,
    pub(super) conductance: Csr,
    pub(super) susceptance: Csr,
}

/// A prescribed, frequency-independent local law `v_out = Y(x) p`.
///
/// `Y` is a finite nonnegative specific conductance [m/(Pa s)]. The original
/// surface integral is retained, so nonuniform pressure and wall properties
/// remain correlated. This law is causal and memoryless: it is deliberately
/// distinct from a frequency sample of [`SpatialBoundaryAdmittance`]. A real
/// value at one frequency alone does not establish a broadband wall law.
pub struct SpatialResistiveBoundary<'m, 'f, 'c> {
    model: &'m HeterogeneousHelmholtzFe<'f, 'c>,
    quadrature_order: usize,
    pub(super) faces: Vec<[u32; 3]>,
    pub(super) conductance: Csr,
}

impl<'m, 'f, 'c> SpatialResistiveBoundary<'m, 'f, 'c> {
    /// Original boundary triangles in canonical order.
    #[must_use]
    pub fn faces(&self) -> &[[u32; 3]] {
        &self.faces
    }

    /// Gauss points per surface Duffy axis used for this material field.
    #[must_use]
    pub const fn quadrature_order(&self) -> usize {
        self.quadrature_order
    }

    /// Use the SAME memoryless physical law in a harmonic solve.
    ///
    /// The already sampled original trace form is cloned without calling the
    /// material again. This direction is physically defined because the law is
    /// frequency independent; converting an arbitrary harmonic sample into a
    /// transient boundary is deliberately not provided.
    pub fn at_frequency(
        &self,
        frequency_hz: f64,
    ) -> Result<SpatialBoundaryAdmittance<'m, 'f, 'c>, HeterogeneousError> {
        let omega = core::f64::consts::TAU * frequency_hz;
        if !omega.is_finite()
            || omega <= 0.0
            || !(omega * omega).is_finite()
            || omega * omega == 0.0
        {
            return Err(invalid("positive representable harmonic frequency required").into());
        }
        Ok(SpatialBoundaryAdmittance {
            model: self.model,
            frequency_hz,
            quadrature_order: self.quadrature_order,
            faces: self.faces.clone(),
            conductance: self.conductance.clone(),
            susceptance: Coo::new(self.model.ndof(), self.model.ndof()).assemble(),
        })
    }

    pub(super) fn validate(
        &self,
        model: &HeterogeneousHelmholtzFe<'_, '_>,
    ) -> Result<(), HeterogeneousError> {
        if !core::ptr::eq(self.model, model) {
            return Err(
                invalid("spatial resistance belongs to another spatial/material operator").into(),
            );
        }
        Ok(())
    }
}

impl SpatialBoundaryAdmittance<'_, '_, '_> {
    /// Exact frequency [Hz] at which the original wall callback was evaluated.
    #[must_use]
    pub const fn frequency_hz(&self) -> f64 {
        self.frequency_hz
    }

    /// Gauss points per surface Duffy axis; material/geometry convergence is
    /// independent of the pressure-basis order.
    #[must_use]
    pub const fn quadrature_order(&self) -> usize {
        self.quadrature_order
    }

    /// Canonically sorted original boundary triangles, without orientation.
    #[must_use]
    pub fn faces(&self) -> &[[u32; 3]] {
        &self.faces
    }

    pub(super) fn validate(
        &self,
        model: &HeterogeneousHelmholtzFe<'_, '_>,
        frequency_hz: f64,
    ) -> Result<(), HeterogeneousError> {
        if !core::ptr::eq(self.model, model) {
            return Err(
                invalid("spatial admittance belongs to another spatial/material operator").into(),
            );
        }
        if self.frequency_hz != frequency_hz {
            return Err(invalid("spatial admittance was sampled at a different frequency").into());
        }
        Ok(())
    }
}

impl<'f, 'c> HeterogeneousHelmholtzFe<'f, 'c> {
    /// Integrate a spatial wall law on the ORIGINAL physical boundary.
    ///
    /// The callback receives frequency [Hz], a canonical original face, and
    /// the original pointwise position, outward normal and surface metric.
    /// It returns specific admittance [m/(Pa s)] in exp(+i omega t), with
    /// finite components and Re(Y) >= 0. Its error propagates unchanged.
    /// Each physical quadrature point is sampled exactly once. Curved faces
    /// retain all original midside coordinates; no corner-plane or averaged
    /// material is substituted.
    ///
    /// Faces must be nonempty, unique exterior triangles. Quadrature is
    /// explicit: at least pressure_order + 1 (and the curved geometry's
    /// declared minimum), at most 16 Gauss points per axis. That lower bound
    /// integrates a CONSTANT law on affine faces; spatial variation and curved
    /// metrics require their own quadrature convergence. Material jumps must
    /// be resolved by the chosen face partition/rule. Cancellation is polled
    /// at face, quadrature and assembly boundaries; no partial owner escapes.
    pub fn spatial_boundary_admittance<'m, E: From<AcousticFeError>>(
        &'m self,
        frequency_hz: f64,
        faces: &[[u32; 3]],
        quadrature_order: usize,
        admittance_at: &mut impl FnMut(f64, [u32; 3], GeometryFacePoint) -> Result<C64, E>,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<SpatialBoundaryAdmittance<'m, 'f, 'c>, E> {
        let poll = |cancelled: &mut dyn FnMut() -> bool| {
            if cancelled() {
                Err(AcousticFeError::Cancelled)
            } else {
                Ok(())
            }
        };
        poll(cancelled)?;
        let omega = core::f64::consts::TAU * frequency_hz;
        if !omega.is_finite()
            || omega <= 0.0
            || !(omega * omega).is_finite()
            || omega * omega == 0.0
        {
            return Err(invalid("positive representable harmonic frequency required").into());
        }
        let (faces, conductance, susceptance) = self.sample_boundary_admittance(
            faces,
            quadrature_order,
            &mut |face, point| admittance_at(frequency_hz, face, point),
            cancelled,
        )?;
        Ok(SpatialBoundaryAdmittance {
            model: self,
            frequency_hz,
            quadrature_order,
            faces,
            conductance,
            susceptance,
        })
    }

    /// Prepare a spatially varying, passive MEMORYLESS boundary law once.
    ///
    /// The callback returns `Y(x)` [m/(Pa s)] at each original curved or affine
    /// face quadrature point. It must describe a frequency-independent local
    /// resistance, not the real part of a dispersive impedance sample. There
    /// is no face averaging, representative normal or replacement geometry.
    /// Sampling and quadrature admission match [`Self::spatial_boundary_admittance`].
    /// Subsequent midpoint/exact steps reuse the original sparse form without
    /// material callbacks or added time-step work. A callback refusal propagates
    /// unchanged and cancellation returns no partially prepared boundary.
    pub fn spatial_resistive_boundary<'m, E: From<AcousticFeError>>(
        &'m self,
        faces: &[[u32; 3]],
        quadrature_order: usize,
        conductance_at: &mut impl FnMut([u32; 3], GeometryFacePoint) -> Result<f64, E>,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<SpatialResistiveBoundary<'m, 'f, 'c>, E> {
        let (faces, conductance, _) = self.sample_boundary_admittance(
            faces,
            quadrature_order,
            &mut |face, point| conductance_at(face, point).map(C64::from_re),
            cancelled,
        )?;
        Ok(SpatialResistiveBoundary {
            model: self,
            quadrature_order,
            faces,
            conductance,
        })
    }

    fn sample_boundary_admittance<E: From<AcousticFeError>>(
        &self,
        faces: &[[u32; 3]],
        quadrature_order: usize,
        admittance_at: &mut impl FnMut([u32; 3], GeometryFacePoint) -> Result<C64, E>,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(Vec<[u32; 3]>, Csr, Csr), E> {
        let poll = |cancelled: &mut dyn FnMut() -> bool| {
            if cancelled() {
                Err(AcousticFeError::Cancelled)
            } else {
                Ok(())
            }
        };
        poll(cancelled)?;
        let minimum = self.geometry.face_quadrature_order(self.geometry.order + 1);
        if quadrature_order < minimum || quadrature_order > MAX_QUADRATURE_ORDER {
            return Err(
                invalid("spatial admittance quadrature order outside admitted range").into(),
            );
        }
        if faces.is_empty() {
            return Err(invalid("spatial admittance requires nonempty boundary faces").into());
        }
        let mut canonical = BTreeSet::new();
        for &face in faces {
            poll(cancelled)?;
            let mut key = face;
            key.sort_unstable();
            if !canonical.insert(key) {
                return Err(invalid("spatial admittance boundary faces must be unique").into());
            }
            self.geometry.boundary_owner(key)?;
        }
        let quadrature = tri_quadrature(quadrature_order);
        let mut real = Coo::new(self.ndof(), self.ndof());
        let mut imaginary = Coo::new(self.ndof(), self.ndof());
        for &face in &canonical {
            poll(cancelled)?;
            let (element, omitted) = self.geometry.boundary_owner(face)?;
            let local: Vec<_> = (0..4).filter(|&a| a != omitted).collect();
            let rule = self.geometry.face_rule(face, &quadrature)?;
            let count = rule.dofs.len();
            let mut re = vec![Dd::from_f64(0.0); count * count];
            let mut im = vec![Dd::from_f64(0.0); count * count];
            for (q, &(face_barycentric, _)) in quadrature.iter().enumerate() {
                poll(cancelled)?;
                let mut barycentric = [0.0; 4];
                for i in 0..3 {
                    barycentric[local[i]] = face_barycentric[i];
                }
                let point = self
                    .geometry
                    .geometry_face_point(element, omitted, barycentric)?;
                let y = admittance_at(face, point)?;
                let weight = rule.weights[q];
                if !finite(y) || y.re < 0.0 || !weight.is_finite() || weight <= 0.0 {
                    return Err(invalid(
                        "spatial admittance requires finite passive samples and positive measure",
                    )
                    .into());
                }
                for i in 0..count {
                    for j in i..count {
                        let product = Dd::from_f64(weight)
                            * Dd::from_f64(rule.phi[i][q])
                            * Dd::from_f64(rule.phi[j][q]);
                        re[i * count + j] = re[i * count + j] + product * Dd::from_f64(y.re);
                        im[i * count + j] = im[i * count + j] + product * Dd::from_f64(y.im);
                    }
                }
            }
            for (i, &row) in rule.dofs.iter().enumerate() {
                poll(cancelled)?;
                for (j, &col) in rule.dofs.iter().enumerate().skip(i) {
                    for (coo, value) in [
                        (&mut real, re[i * count + j].to_f64()),
                        (&mut imaginary, im[i * count + j].to_f64()),
                    ] {
                        if !value.is_finite() {
                            return Err(invalid("spatial admittance assembly overflow").into());
                        }
                        if value != 0.0 {
                            coo.push(row, col, value);
                            if i != j {
                                coo.push(col, row, value);
                            }
                        }
                    }
                }
            }
        }
        poll(cancelled)?;
        let conductance = real.assemble();
        let susceptance = imaginary.assemble();
        for row in 0..self.ndof() {
            poll(cancelled)?;
            if conductance
                .row(row)
                .1
                .iter()
                .chain(susceptance.row(row).1)
                .any(|v| !v.is_finite())
            {
                return Err(invalid("spatial admittance global assembly overflow").into());
            }
        }
        Ok((canonical.into_iter().collect(), conductance, susceptance))
    }
}

#[cfg(test)]
mod tests;
