//! Both-sided pressure and Euler velocity from the original adaptive P0 rows.
use super::*;
use crate::{helmholtz::bc::IncidentSample, near_field};

/// Source-bound single-layer rows and the matched incident pressure/velocity.
/// No Green integration or source sampling is repeated during evaluation.
pub struct PreparedSoftScreenReceivers<'b, 'p, 's, 'i> {
    source: &'b SoftScreenScattering<'p, 's, 'i>,
    points: Vec<Point>,
    incident: Vec<IncidentSample>,
    rows: Vec<near_field::RayleighSingleLayer>,
    kernel_evaluations: usize,
}

impl<'p, 's, 'i> SoftScreenScattering<'p, 's, 'i> {
    /// Prepare adaptive G and analytic grad(G) at points on EITHER screen side.
    /// Every gradient component has its own quadrature acceptance gate. Exact
    /// face/edge clearance and the original receiver/work ceilings apply.
    ///
    /// `maximum_row_bytes` covers retained raw rows, points, incident samples
    /// and owner; it excludes borrowed geometry/solution, temporary triangles
    /// (72 bytes/facet), output fields and the original integration stack.
    /// Discrepancies bound neither source/mesh error nor relative dark-field
    /// cancellation error. An ideal pressure-release sheet has zero thickness.
    pub fn prepare_receivers<'b>(
        &'b self,
        points: &[Point],
        separation: f64,
        options: near_field::Options,
        maximum_row_bytes: usize,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<PreparedSoftScreenReceivers<'b, 'p, 's, 'i>, HelmholtzError> {
        poll(&mut cancelled)?;
        options.validate()?;
        let n = self.prepared.unknowns();
        let bytes = n
            .checked_mul(points.len())
            .and_then(|v| v.checked_mul(core::mem::size_of::<near_field::RayleighSingleLayer>()))
            .and_then(|v| {
                v.checked_add(
                    points.len()
                        * (core::mem::size_of::<Point>() + core::mem::size_of::<IncidentSample>()),
                )
            })
            .and_then(|v| {
                v.checked_add(core::mem::size_of::<
                    PreparedSoftScreenReceivers<'_, '_, '_, '_>,
                >())
            });
        if points.len() > near_field::MAX_RECEIVERS
            || bytes.is_none_or(|v| v > maximum_row_bytes)
            || n.checked_mul(points.len())
                .and_then(|v| v.checked_mul(80))
                .is_none_or(|v| v > options.maximum_kernel_evaluations)
            || self.incident.source_count() as u128 * points.len() as u128
                > self.maximum_source_point_products as u128
        {
            return Err(bad(
                "soft-screen receiver rows exceed storage or source/quadrature work",
            ));
        }
        self.prepared
            .admit_points(points, separation, &mut cancelled)?;
        let rho_omega = self.k * self.medium.sound_speed * self.medium.density;
        let mut incident = Vec::with_capacity(points.len());
        for &point in points {
            // Refuse a source/receiver coincidence before adaptive quadrature.
            incident.push(self.incident.sample_wave(
                self.k,
                rho_omega,
                Some(self.medium.density * self.medium.sound_speed),
                point,
                &mut cancelled,
            )?);
        }
        let triangles: Vec<_> = self.prepared.facets.iter().map(|f| f.triangle).collect();
        let (rows, kernel_evaluations) =
            near_field::single_layer_rows(self.k, &triangles, points, options, &mut cancelled)?;
        poll(&mut cancelled)?;
        Ok(PreparedSoftScreenReceivers {
            source: self,
            points: points.to_vec(),
            incident,
            rows,
            kernel_evaluations,
        })
    }
}

impl PreparedSoftScreenReceivers<'_, '_, '_, '_> {
    pub fn points_m(&self) -> &[Point] {
        &self.points
    }
    pub fn kernel_evaluations(&self) -> usize {
        self.kernel_evaluations
    }

    pub(super) fn incident_samples(&self) -> &[IncidentSample] {
        &self.incident
    }

    /// Incident-plus-scattered p [Pa] and Euler velocity [m/s]. The incident
    /// field is added once, retaining interference and reactive near fields.
    /// `near_field::FirstOrder` consumes these same physical p/u components.
    pub fn evaluate(
        &self,
        cancelled: impl FnMut() -> bool,
    ) -> Result<near_field::VectorObservation, HelmholtzError> {
        self.evaluate_fields(true, cancelled)
    }

    /// Outgoing scattered field only, on the same original receiver rows.
    pub fn evaluate_scattered(
        &self,
        cancelled: impl FnMut() -> bool,
    ) -> Result<near_field::VectorObservation, HelmholtzError> {
        self.evaluate_fields(false, cancelled)
    }

    fn evaluate_fields(
        &self,
        with_incident: bool,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<near_field::VectorObservation, HelmholtzError> {
        poll(&mut cancelled)?;
        let n = self.source.density.len();
        let rho_omega = self.source.k * self.source.medium.sound_speed * self.source.medium.density;
        let inverse = C64::new(0.0, -1.0 / rho_omega);
        let mut pressure = Vec::with_capacity(self.points.len());
        let mut pressure_error = Vec::with_capacity(self.points.len());
        let mut velocity = Vec::with_capacity(self.points.len());
        let mut velocity_error = Vec::with_capacity(self.points.len());
        for (index, row) in self.rows.chunks_exact(n).enumerate() {
            poll(&mut cancelled)?;
            let mut sums = [ComplexSum::default(); 4];
            let mut errors = [fs_math::dd::Dd::ZERO; 4];
            for (panel, (cell, &density)) in row.iter().zip(&self.source.density).enumerate() {
                if panel % 128 == 0 {
                    poll(&mut cancelled)?;
                }
                sums[0].add_product(cell.value, density);
                errors[0] = errors[0] + fs_math::dd::Dd::from_f64(cell.error_value * density.abs());
                for axis in 0..3 {
                    sums[axis + 1].add_product(cell.gradient[axis], density);
                    errors[axis + 1] = errors[axis + 1]
                        + fs_math::dd::Dd::from_f64(cell.error_gradient[axis] * density.abs());
                }
            }
            // Add the incident pressure before rounding the coherent sum.
            if with_incident {
                sums[0].add(self.incident[index].pressure_pa);
            }
            let p = sums[0].finish();
            let mut u = [C64::ZERO; 3];
            for axis in 0..3 {
                let mut sum = ComplexSum::default();
                sum.add(sums[axis + 1].finish() * inverse);
                if with_incident {
                    sum.add(self.incident[index].velocity_m_s[axis]);
                }
                u[axis] = sum.finish();
            }
            let errors = errors.map(|e| e.to_f64());
            let ue = core::array::from_fn(|axis| errors[axis + 1] / rho_omega);
            if !finite(p)
                || u.iter().any(|&v| !finite(v))
                || errors.iter().chain(&ue).any(|v| !v.is_finite())
            {
                return Err(bad("soft-screen receiver field or discrepancy overflows"));
            }
            pressure.push(p);
            pressure_error.push(errors[0]);
            velocity.push(u);
            velocity_error.push(ue);
        }
        poll(&mut cancelled)?;
        Ok(near_field::VectorObservation {
            scalar: near_field::Observation {
                pressure,
                quadrature_error_estimate_pa: pressure_error,
            },
            particle_velocity_m_s: velocity,
            quadrature_error_estimate_m_s: velocity_error,
        })
    }
}
