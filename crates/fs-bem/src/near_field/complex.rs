//! The original adaptive receiver rows for a source-bound molecular field.
//!
//! Outer owners bind material, frequency and surface identity. This internal
//! carrier shares geometry and quadrature with real-k receivers and admits
//! physical momentum independently from the propagation wavenumber.
use super::*;
use crate::helmholtz::laplace::DampedWaveNumber;

fn finite(value: C64) -> bool {
    value.re.is_finite() && value.im.is_finite()
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), HelmholtzError> {
    if cancelled() {
        Err(HelmholtzError::Cancelled)
    } else {
        Ok(())
    }
}

/// Logical row payload, including row descriptors and the prepared owner.
/// Geometry admission scratch, allocator metadata and stacks remain excluded.
pub(crate) fn required_row_bytes(panels: usize, receivers: usize) -> Option<usize> {
    panels
        .checked_mul(core::mem::size_of::<Integral>())?
        .checked_add(core::mem::size_of::<Vec<Integral>>())?
        .checked_mul(receivers)?
        .checked_add(core::mem::size_of::<Prepared<'_>>())
}

/// Direct pressure/gradient contractions and propagated quadrature discrepancy.
/// The discrepancies exclude spatial discretization, solve and material errors.
pub(crate) struct Fields {
    pub(crate) pressure_pa: Vec<C64>,
    pub(crate) gradient_pa_per_m: Vec<[C64; 3]>,
    pub(crate) quadrature_error_estimate_pa: Vec<f64>,
    pub(crate) quadrature_error_estimate_pa_per_m: Vec<[f64; 3]>,
}

pub(crate) struct Prepared<'a> {
    surface: &'a SpherePanels,
    rho_omega: f64,
    rows: Vec<Vec<Integral>>,
    kernel_evaluations: usize,
}

impl<'a> Prepared<'a> {
    /// Geometry has already admitted exterior or interior winding/clearance.
    /// Complex |k| resolves both phase and attenuation on each adaptive tile.
    /// No centroid fallback or partial prepared bank is returned on refusal.
    pub(crate) fn new(
        geometry: &Geometry<'a>,
        wave: DampedWaveNumber,
        rho_omega: f64,
        options: Options,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, HelmholtzError> {
        poll(cancelled)?;
        options.validate()?;
        if !finite(wave.k)
            || wave.k.re <= 0.0
            || wave.k.im < 0.0
            || wave.phase_origin_m != 0.0
            || !rho_omega.is_finite()
            || rho_omega <= 0.0
        {
            return Err(bad(
                "complex receiver needs outgoing unshifted k and real positive rho*omega",
            ));
        }
        let surface = geometry.surface;
        let triangles = surface
            .triangles()
            .ok_or_else(|| bad("missing retained receiver triangles"))?;
        if triangles
            .len()
            .checked_mul(geometry.points.len())
            .and_then(|count| count.checked_mul(80))
            .is_none_or(|count| count > options.maximum_kernel_evaluations)
            || required_row_bytes(triangles.len(), geometry.points.len()).is_none()
        {
            return Err(bad(
                "near-field kernel-evaluation or row-payload budget exhausted",
            ));
        }
        let mut work = Work {
            options,
            evaluations: 0,
        };
        let mut rows = Vec::with_capacity(geometry.points.len());
        for &point in &geometry.points {
            poll(cancelled)?;
            let mut row = Vec::with_capacity(triangles.len());
            for (panel, &triangle) in triangles.iter().enumerate() {
                row.push(integrate_fields_cancellable(
                    wave,
                    point,
                    triangle,
                    surface.normals()[panel],
                    0,
                    &mut work,
                    true,
                    cancelled,
                )?);
            }
            rows.push(row);
        }
        poll(cancelled)?;
        Ok(Self {
            surface,
            rho_omega,
            rows,
            kernel_evaluations: work.evaluations,
        })
    }

    pub(crate) fn kernel_evaluations(&self) -> usize {
        self.kernel_evaluations
    }

    /// Contract original P0 physical traces. The outer owner checks source,
    /// geometry and frequency binding; this layer checks extents and values.
    pub(crate) fn evaluate_fields(
        &self,
        pressure: &[C64],
        actual_velocity: &[C64],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Fields, HelmholtzError> {
        poll(cancelled)?;
        let panels = self.surface.areas().len();
        if pressure.len() != panels || actual_velocity.len() != panels {
            return Err(bad(
                "complex receiver trace extent differs from original surface",
            ));
        }
        for (panel, (&p, &v)) in pressure.iter().zip(actual_velocity).enumerate() {
            if panel % 128 == 0 {
                poll(cancelled)?;
            }
            if !finite(p) || !finite(v) {
                return Err(bad("nonfinite complex receiver trace"));
            }
        }
        let receivers = self.rows.len();
        let mut result = Fields {
            pressure_pa: Vec::with_capacity(receivers),
            gradient_pa_per_m: Vec::with_capacity(receivers),
            quadrature_error_estimate_pa: Vec::with_capacity(receivers),
            quadrature_error_estimate_pa_per_m: Vec::with_capacity(receivers),
        };
        for row in &self.rows {
            poll(cancelled)?;
            let mut p_sum = ComplexSum::default();
            let mut gradient = [ComplexSum::default(); 3];
            let mut p_error = 0.0;
            let mut gradient_error = [0.0; 3];
            for (panel, ((cell, &p), &v)) in
                row.iter().zip(pressure).zip(actual_velocity).enumerate()
            {
                if panel % 128 == 0 {
                    poll(cancelled)?;
                }
                let q = v * C64::new(0.0, self.rho_omega);
                p_sum.add(cell.d * p - cell.s * q);
                p_error += p.abs() * cell.error_d + q.abs() * cell.error_s;
                for axis in 0..3 {
                    gradient[axis].add(cell.gradient.d[axis] * p - cell.gradient.s[axis] * q);
                    gradient_error[axis] += p.abs() * cell.gradient.error_d[axis]
                        + q.abs() * cell.gradient.error_s[axis];
                }
            }
            let p = p_sum.finish();
            let gradient = gradient.map(ComplexSum::finish);
            if !finite(p)
                || gradient.iter().any(|&value| !finite(value))
                || !p_error.is_finite()
                || gradient_error.iter().any(|value| !value.is_finite())
            {
                return Err(bad("nonfinite complex receiver field or discrepancy"));
            }
            result.pressure_pa.push(p);
            result.gradient_pa_per_m.push(gradient);
            result.quadrature_error_estimate_pa.push(p_error);
            result
                .quadrature_error_estimate_pa_per_m
                .push(gradient_error);
        }
        poll(cancelled)?;
        Ok(result)
    }

    /// Reuse the direct gradient buffers for Euler particle velocity. Neither
    /// rho*c_eq*k nor a plane-wave pressure-to-velocity shortcut is introduced.
    pub(crate) fn evaluate(
        &self,
        pressure: &[C64],
        actual_velocity: &[C64],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<VectorObservation, HelmholtzError> {
        let mut fields = self.evaluate_fields(pressure, actual_velocity, cancelled)?;
        let inverse = C64::new(0.0, -1.0 / self.rho_omega);
        for (velocity, error) in fields
            .gradient_pa_per_m
            .iter_mut()
            .zip(&mut fields.quadrature_error_estimate_pa_per_m)
        {
            poll(cancelled)?;
            for axis in 0..3 {
                velocity[axis] = velocity[axis] * inverse;
                error[axis] /= self.rho_omega;
                if !finite(velocity[axis]) || !error[axis].is_finite() {
                    return Err(bad("nonfinite complex receiver particle velocity"));
                }
            }
        }
        poll(cancelled)?;
        Ok(VectorObservation {
            scalar: Observation {
                pressure: fields.pressure_pa,
                quadrature_error_estimate_pa: fields.quadrature_error_estimate_pa,
            },
            particle_velocity_m_s: fields.gradient_pa_per_m,
            quadrature_error_estimate_m_s: fields.quadrature_error_estimate_pa_per_m,
        })
    }
}

#[cfg(test)]
mod tests;
