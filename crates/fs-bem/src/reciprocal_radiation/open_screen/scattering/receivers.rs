//! Cached analytic illumination plus the original adaptive scattered field.
use super::*;
use crate::{helmholtz::bc::IncidentSample, near_field};

/// One P1 quadrature bank shared by every incident case, retaining the exact
/// source/geometry identity. Analytic incident p/u is cached at preparation.
pub struct PreparedScreenScatteringReceivers<'b, 'p, 's, 'i> {
    source: &'b ScreenScattering<'p, 's, 'i>,
    rows: PreparedScreenReceivers<'b, 'p, 's>,
    incident: Vec<Vec<IncidentSample>>,
}

impl<'p, 's, 'i> ScreenScattering<'p, 's, 'i> {
    /// Prepare total pressure and Cartesian particle velocity for all cases.
    /// `maximum_row_bytes` includes the original P1 rows, points, every cached
    /// analytic incident sample and this owner. Other exclusions are those of
    /// JumpSolution::prepare_receivers. The original source-point work ceiling
    /// also bounds this complete analytic receiver bank.
    ///
    /// Finite-emitter singularities refuse before receiver quadrature. Error
    /// estimates cover integration of the solved scattered jump, not incident
    /// source error, boundary discretization or relative cancellation error.
    pub fn prepare_receivers<'b>(
        &'b self,
        points_m: &[Point],
        separation_m: f64,
        options: near_field::Options,
        maximum_row_bytes: usize,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<PreparedScreenScatteringReceivers<'b, 'p, 's, 'i>, HelmholtzError> {
        poll(&mut cancelled)?;
        options.validate()?;
        let cases = self.cases();
        let terms: u128 = self
            .incidents
            .iter()
            .map(|source| source.source_count() as u128)
            .sum();
        let extra = cases
            .checked_mul(points_m.len())
            .and_then(|n| n.checked_mul(core::mem::size_of::<IncidentSample>()))
            .and_then(|bytes| {
                cases
                    .checked_mul(core::mem::size_of::<Vec<IncidentSample>>())
                    .and_then(|vectors| bytes.checked_add(vectors))
            })
            .and_then(|bytes| {
                bytes.checked_add(
                    core::mem::size_of::<PreparedScreenScatteringReceivers<'_, '_, '_, '_>>()
                        - core::mem::size_of::<PreparedScreenReceivers<'_, '_, '_>>(),
                )
            })
            .ok_or_else(|| bad("open screen total receiver payload overflows"))?;
        if points_m.len() > near_field::MAX_RECEIVERS
            || terms * points_m.len() as u128 > self.maximum_source_point_products as u128
            || extra > maximum_row_bytes
        {
            return Err(bad(
                "open screen total receivers exceed source work or row payload",
            ));
        }
        let source = &self.solutions[0];
        let rho_c = source.medium.density * source.medium.sound_speed;
        let rho_omega = source.k.re * source.medium.sound_speed * source.medium.density;
        let mut incident = Vec::with_capacity(cases);
        for field in &self.incidents {
            poll(&mut cancelled)?;
            let mut samples = Vec::with_capacity(points_m.len());
            for &point in points_m {
                samples.push(field.sample_wave(
                    source.k.re,
                    rho_omega,
                    Some(rho_c),
                    point,
                    &mut cancelled,
                )?);
            }
            incident.push(samples);
        }
        let rows = source.prepare_receivers(
            points_m,
            separation_m,
            options,
            maximum_row_bytes - extra,
            &mut cancelled,
        )?;
        poll(&mut cancelled)?;
        Ok(PreparedScreenScatteringReceivers {
            source: self,
            rows,
            incident,
        })
    }
}

impl PreparedScreenScatteringReceivers<'_, '_, '_, '_> {
    pub fn points_m(&self) -> &[Point] {
        self.rows.points_m()
    }
    pub fn kernel_evaluations(&self) -> usize {
        self.rows.kernel_evaluations()
    }

    /// Physical incident-plus-scattered p [Pa] and Euler velocity [m/s]. Both
    /// components add their matched analytic incident field exactly once.
    /// This repeats no Green kernels, incident sampling or boundary solve.
    pub fn evaluate(
        &self,
        case: usize,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<near_field::VectorObservation, HelmholtzError> {
        poll(&mut cancelled)?;
        let jump = self.source.scattered(case)?.pressure_jump_pa();
        let mut total = self.rows.evaluate_jump(jump, &mut cancelled)?;
        for (index, incident) in self.incident[case].iter().enumerate() {
            poll(&mut cancelled)?;
            total.scalar.pressure[index] = total.scalar.pressure[index] + incident.pressure_pa;
            for axis in 0..3 {
                total.particle_velocity_m_s[index][axis] =
                    total.particle_velocity_m_s[index][axis] + incident.velocity_m_s[axis];
            }
            if !finite(total.scalar.pressure[index])
                || total.particle_velocity_m_s[index]
                    .iter()
                    .any(|&value| !finite(value))
            {
                return Err(bad("open screen total receiver field overflows"));
            }
        }
        poll(&mut cancelled)?;
        Ok(total)
    }
}
