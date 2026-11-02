//! Reusable P1 pressure-jump receiver rows for a passive molecular exterior.
//!
//! The open-screen owner supplies the topology and binds these rows to its
//! immutable solved source. The existing P1 triangle integrator supplies the
//! outgoing pressure and analytic gradient, with separate componentwise
//! quadrature gates. Momentum keeps the original real `rho*omega`.
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

/// Logical retained row payload. Geometry scratch, returned observations,
/// allocator metadata and recursive integration stacks are separate.
pub(crate) fn required_row_bytes(panels: usize, receivers: usize) -> Option<usize> {
    panels
        .checked_mul(receivers)?
        .checked_mul(core::mem::size_of::<Moments>())?
        .checked_add(core::mem::size_of::<OpenRows>())
}

/// Flat receiver-major rows; each original facet retains all three P1 moments.
pub(crate) struct OpenRows {
    rows: Vec<Moments>,
    panels: usize,
    rho_omega: f64,
    kernel_evaluations: usize,
}

impl OpenRows {
    /// Admit the complete bank before allocation or quadrature. Topology is
    /// owned by the original screen; this layer checks exact face/edge distance,
    /// representable geometry, derivative work and retained row payload.
    #[allow(clippy::too_many_arguments)]
    pub(crate) fn new(
        triangles: &[Triangle],
        normals: &[Point],
        wave: DampedWaveNumber,
        rho_omega: f64,
        points: &[Point],
        separation: f64,
        options: Options,
        maximum_row_bytes: usize,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, HelmholtzError> {
        poll(cancelled)?;
        options.validate()?;
        let panels = triangles.len();
        if panels == 0
            || panels > MAX_PANELS
            || normals.len() != panels
            || points.len() > MAX_RECEIVERS
            || !finite(wave.k)
            || wave.k.re <= 0.0
            || wave.k.im < 0.0
            || wave.phase_origin_m != 0.0
            || !rho_omega.is_finite()
            || rho_omega <= 0.0
            || !rho_omega.recip().is_finite()
            || !separation.is_finite()
            || separation <= 0.0
        {
            return Err(bad("invalid molecular open-screen receiver request"));
        }
        let entries = panels
            .checked_mul(points.len())
            .ok_or_else(|| bad("open-screen receiver row extent overflow"))?;
        if entries
            .checked_mul(80)
            .is_none_or(|count| count > options.maximum_kernel_evaluations)
            || required_row_bytes(panels, points.len())
                .is_none_or(|bytes| bytes > maximum_row_bytes)
        {
            return Err(bad(
                "open-screen receiver rows exceed their payload or quadrature budget",
            ));
        }
        let mut coordinate_scale = 1.0_f64;
        for (i, (&triangle, &normal)) in triangles.iter().zip(normals).enumerate() {
            if i % 32 == 0 {
                poll(cancelled)?;
            }
            for value in triangle.iter().flatten() {
                if !value.is_finite() || value.abs() > 1e9 {
                    return Err(bad("unrepresentable open-screen receiver geometry"));
                }
                coordinate_scale = coordinate_scale.max(value.abs());
            }
            let area_vector = cross(sub(triangle[1], triangle[0]), sub(triangle[2], triangle[0]));
            let area2 = norm(area_vector);
            if !area2.is_finite()
                || area2 <= 0.0
                || normal.iter().any(|v| !v.is_finite())
                || (norm(normal) - 1.0).abs() > 1e-10
                || (dot(area_vector.map(|v| v / area2), normal) - 1.0).abs() > 1e-10
            {
                return Err(bad("invalid oriented open-screen receiver triangle"));
            }
        }
        for &point in points {
            poll(cancelled)?;
            for value in point {
                if !value.is_finite() || value.abs() > 1e9 {
                    return Err(bad("unrepresentable open-screen receiver position"));
                }
                coordinate_scale = coordinate_scale.max(value.abs());
            }
        }
        if separation <= 64.0 * f64::EPSILON * coordinate_scale {
            return Err(bad("open-screen clearance is below coordinate resolution"));
        }
        for &point in points {
            for (panel, &triangle) in triangles.iter().enumerate() {
                if panel % 32 == 0 {
                    poll(cancelled)?;
                }
                let clearance = distance(point, triangle);
                if !clearance.is_finite() || clearance <= separation {
                    return Err(bad("open-screen receiver touches surface"));
                }
            }
        }
        let mut rows = Vec::new();
        rows.try_reserve_exact(entries)
            .map_err(|_| HelmholtzError::AllocationFailed {
                what: "molecular open-screen receiver moments",
            })?;
        let mut work = Work {
            options,
            evaluations: 0,
        };
        for &point in points {
            poll(cancelled)?;
            for (&triangle, &normal) in triangles.iter().zip(normals) {
                rows.push(integrate_linear(
                    wave,
                    point,
                    triangle,
                    normal,
                    [[1., 0., 0.], [0., 1., 0.], [0., 0., 1.]],
                    0,
                    &mut work,
                    true,
                    cancelled,
                )?);
            }
        }
        poll(cancelled)?;
        Ok(Self {
            rows,
            panels,
            rho_omega,
            kernel_evaluations: work.evaluations,
        })
    }

    pub(crate) fn kernel_evaluations(&self) -> usize {
        self.kernel_evaluations
    }

    /// Contract the original facet-corner jump without making a second copy
    /// of it. The source owner validates nodal extents and identity. Source
    /// values are also checked when the receiver list is empty.
    pub(crate) fn evaluate(
        &self,
        source: impl Fn(usize) -> [C64; 3],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<VectorObservation, HelmholtzError> {
        poll(cancelled)?;
        for panel in 0..self.panels {
            if panel % 128 == 0 {
                poll(cancelled)?;
            }
            if source(panel).iter().any(|&v| !finite(v)) {
                return Err(bad("nonfinite open-screen pressure jump"));
            }
        }
        let receivers = self.rows.len() / self.panels;
        let mut pressure = Vec::with_capacity(receivers);
        let mut pressure_error = Vec::with_capacity(receivers);
        let mut velocity = Vec::with_capacity(receivers);
        let mut velocity_error = Vec::with_capacity(receivers);
        let inverse = C64::new(0.0, -1.0 / self.rho_omega);
        for row in self.rows.chunks_exact(self.panels) {
            poll(cancelled)?;
            let mut sums = [ComplexSum::default(); 4];
            let mut errors = [0.0; 4];
            let mut corrections = [0.0; 4];
            for (panel, moments) in row.iter().enumerate() {
                if panel % 128 == 0 {
                    poll(cancelled)?;
                }
                for (cell, jump) in moments.iter().zip(source(panel)) {
                    sums[0].add(cell.d * jump);
                    compensated(
                        &mut errors[0],
                        &mut corrections[0],
                        cell.error_d * jump.abs(),
                    );
                    for axis in 0..3 {
                        sums[axis + 1].add(cell.gradient.d[axis] * jump);
                        compensated(
                            &mut errors[axis + 1],
                            &mut corrections[axis + 1],
                            cell.gradient.error_d[axis] * jump.abs(),
                        );
                    }
                }
            }
            let fields = sums.map(ComplexSum::finish);
            let vector = core::array::from_fn(|axis| fields[axis + 1] * inverse);
            let estimates = core::array::from_fn(|axis| errors[axis + 1] / self.rho_omega);
            if fields.iter().chain(&vector).any(|&v| !finite(v))
                || errors.iter().chain(&estimates).any(|v| !v.is_finite())
            {
                return Err(bad("nonfinite open-screen field or quadrature discrepancy"));
            }
            pressure.push(fields[0]);
            pressure_error.push(errors[0]);
            velocity.push(vector);
            velocity_error.push(estimates);
        }
        poll(cancelled)?;
        Ok(VectorObservation {
            scalar: Observation {
                pressure,
                quadrature_error_estimate_pa: pressure_error,
            },
            particle_velocity_m_s: velocity,
            quadrature_error_estimate_m_s: velocity_error,
        })
    }
}

#[cfg(test)]
mod tests;
