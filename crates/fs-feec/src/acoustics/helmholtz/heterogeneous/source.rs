//! Distributed volume injection in the original heterogeneous pressure space.
//!
//! With a prescribed volume-rate density `s` [1/s], linear continuity is
//! `beta p_dot + div(v) = s`. Thus the pressure-potential equation receives
//! `b_i = integral(phi_i s)`, and the harmonic pressure equation receives
//! `i omega b`. Its conjugate work is `integral(p s)`, with no extra density,
//! speed, cell-volume or port-area factor. This is a prescribed acoustic
//! source, not an inferred mass, heat-release, phase-change or chemistry law.
//! The quiescent first-order pressure/velocity source convention is described
//! in the COMSOL Acoustics Module User's Guide, "Pressure Acoustics, Time
//! Explicit Model":
//! <https://doc.comsol.com/6.4/doc/com.comsol.help.aco/aco_ug_pressure.05.100.html>.

use super::{
    HeterogeneousError, HeterogeneousHelmholtzFe, MAX_QUADRATURE_ORDER, duffy_quadrature, invalid,
    poll,
};
use fs_math::{c64::C64, dd::Dd};

/// Bound on simultaneously driven volume-source columns. The exact trajectory
/// owner's augmented-state limit applies additionally to its preparations.
pub const MAX_VOLUME_SOURCES: usize = 128;

pub(super) fn admit_count(count: usize) -> Result<(), HeterogeneousError> {
    if count > MAX_VOLUME_SOURCES {
        return Err(HeterogeneousError::Limit {
            what: "simultaneous volume acoustic sources",
            requested: count,
            maximum: MAX_VOLUME_SOURCES,
        });
    }
    Ok(())
}

/// Frozen spatial distribution of volume injection on one original fluid mesh.
///
/// Shape density is in m^-3. Multiplying by a volume-flow amplitude [m^3/s]
/// gives the local volume-rate density [1/s]. The signed shape is integrated
/// as supplied: it is not normalized, clipped or replaced by a point source.
/// Preparation samples the callback once; every subsequent drive reuses the
/// same original-basis load. Its model borrow prevents cross-wiring a source
/// onto a different material/geometry owner.
pub struct VolumeVelocitySource<'m, 'f, 'c> {
    pub(super) model: &'m HeterogeneousHelmholtzFe<'f, 'c>,
    pub(super) load: Vec<f64>,
    integral: f64,
}

impl VolumeVelocitySource<'_, '_, '_> {
    /// Dimensionless signed volume integral of the supplied shape.
    #[must_use]
    pub const fn integrated_shape(&self) -> f64 {
        self.integral
    }

    pub(super) fn validate(
        &self,
        model: &HeterogeneousHelmholtzFe<'_, '_>,
    ) -> Result<(), HeterogeneousError> {
        if !core::ptr::eq(self.model, model) {
            return Err(
                invalid("volume source belongs to another spatial/material operator").into(),
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
    /// Project a distributed source onto EVERY original hierarchical test.
    ///
    /// The callback receives the original tetrahedron index and SI position,
    /// and returns a finite signed density [m^-3]. Quadrature has its own
    /// explicit order, between `geometry.order()+2` and 16. Localized sources
    /// require independent source-quadrature and h/p convergence. Curved maps
    /// use their original physical positions and pointwise Jacobian, so the
    /// source rule also resolves the geometry. No material is queried or
    /// reassembled. Source errors propagate unchanged, and
    /// cancellation or nonfinite data publishes no partially prepared source.
    pub fn volume_velocity_source<'m, E: From<HeterogeneousError>>(
        &'m self,
        quadrature_order: usize,
        density_per_m3: &mut impl FnMut(usize, [f64; 3]) -> Result<f64, E>,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<VolumeVelocitySource<'m, 'f, 'c>, E> {
        poll(cancelled)?;
        if quadrature_order < self.geometry.order + 2 || quadrature_order > MAX_QUADRATURE_ORDER {
            return Err(HeterogeneousError::from(invalid(
                "volume source quadrature order outside admitted range",
            ))
            .into());
        }
        // Preserve the existing affine projection; a curved element must
        // never use its corner tetrahedron's constant volume measure.
        let metrics = if self.geometry.is_curved() {
            None
        } else {
            Some(&self.geometry.affine_geometry)
        };
        let quadrature = duffy_quadrature(quadrature_order);
        let mut load = vec![Dd::from_f64(0.0); self.ndof()];
        let mut integral = Dd::from_f64(0.0);
        for (element, vertices) in self.geometry.tetrahedra().iter().enumerate() {
            poll(cancelled)?;
            let dofs = self.geometry.space.element_dofs(element);
            for &(barycentric, weight) in &quadrature {
                poll(cancelled)?;
                let (position, jacobian) = if let Some(metrics) = &metrics {
                    let mut position = [0.0; 3];
                    for (&vertex, &bary) in vertices.iter().zip(&barycentric) {
                        for (axis, coordinate) in position.iter_mut().enumerate() {
                            *coordinate += bary * self.geometry.positions[vertex as usize][axis];
                        }
                    }
                    (position, 6.0 * metrics.vol_signed[element].abs())
                } else {
                    let point = self
                        .geometry
                        .geometry_point(element, barycentric)
                        .map_err(HeterogeneousError::from)?;
                    (point.position_m, point.determinant_m3)
                };
                let density = density_per_m3(element, position)?;
                let weighted = density * weight * jacobian;
                if !density.is_finite() || !weighted.is_finite() {
                    return Err(HeterogeneousError::from(invalid(
                        "nonfinite volume source density",
                    ))
                    .into());
                }
                integral = integral + Dd::from_f64(weighted);
                for &(global, shape) in &dofs {
                    load[global] = load[global]
                        + Dd::from_f64(weighted)
                            * Dd::from_f64(shape.eval(barycentric, self.geometry.order));
                }
            }
        }
        let load: Vec<_> = load.into_iter().map(Dd::to_f64).collect();
        let integral = integral.to_f64();
        if !integral.is_finite() || load.iter().any(|value| !value.is_finite()) {
            return Err(
                HeterogeneousError::from(invalid("volume source projection overflow")).into(),
            );
        }
        poll(cancelled)?;
        Ok(VolumeVelocitySource {
            model: self,
            load,
            integral,
        })
    }
}

#[cfg(test)]
mod tests;
