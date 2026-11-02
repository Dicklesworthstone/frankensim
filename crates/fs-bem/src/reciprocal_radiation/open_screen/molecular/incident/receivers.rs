//! One original P1 row bank and cached analytic incident field for every case.
use super::*;
use crate::helmholtz::bc::IncidentSample;
use crate::near_field;

/// Original geometry/material/source binding, with no receiver rebinding.
/// Incident terms are sampled once at preparation; case evaluation contracts
/// one coherent total jump and adds the matched analytic field exactly once.
pub struct PreparedMolecularScreenIlluminationReceivers<'b, 'p, 's, 'i> {
    source: &'b MolecularScreenIllumination<'p, 's, 'i>,
    rows: PreparedMolecularScreenReceivers<'b, 'p, 's>,
    incident: Vec<Vec<IncidentSample>>,
}

impl<'p, 's, 'i> MolecularScreenIllumination<'p, 's, 'i> {
    /// Prepare shared adaptive rows and each original analytic incident field.
    /// The declared row-byte limit includes the original row bank, receiver
    /// positions, all cached incident samples and this owner. The original
    /// source budget separately bounds total-jump/work scratch. Observation
    /// output is 96 bytes/receiver plus vector descriptors. These are logical
    /// payload bounds, excluding allocator metadata and recursion/thread stacks.
    ///
    /// Emitter singularities refuse before expensive receiver quadrature.
    /// Quadrature error estimates concern the combined outgoing jump only;
    /// they do not bound source model, mesh, solve or analytic evaluation error.
    pub fn prepare_receivers<'b>(
        &'b self,
        points_m: &[Point],
        separation_m: f64,
        options: near_field::Options,
        maximum_row_bytes: usize,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<PreparedMolecularScreenIlluminationReceivers<'b, 'p, 's, 'i>, HelmholtzError> {
        poll(&mut cancelled)?;
        options.validate()?;
        let cases = self.incidents.len();
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
                    core::mem::size_of::<
                        PreparedMolecularScreenIlluminationReceivers<'_, '_, '_, '_>,
                    >() - core::mem::size_of::<PreparedMolecularScreenReceivers<'_, '_, '_>>(),
                )
            })
            .ok_or_else(|| bad("molecular illuminated receiver payload overflows"))?;
        if points_m.len() > near_field::MAX_RECEIVERS
            || terms * points_m.len() as u128 > self.maximum_source_point_products as u128
            || extra > maximum_row_bytes
        {
            return Err(bad(
                "molecular illuminated receivers exceed source work or retained payload",
            ));
        }
        let mut incident = Vec::with_capacity(cases);
        for source in &self.incidents {
            poll(&mut cancelled)?;
            incident.push(source.pressure_velocity_at_many(
                points_m,
                self.maximum_source_point_products,
                &mut cancelled,
            )?);
        }
        let rows = self.motion.prepare_receivers(
            points_m,
            separation_m,
            options,
            maximum_row_bytes - extra,
            &mut cancelled,
        )?;
        poll(&mut cancelled)?;
        Ok(PreparedMolecularScreenIlluminationReceivers {
            source: self,
            rows,
            incident,
        })
    }
}

impl PreparedMolecularScreenIlluminationReceivers<'_, '_, '_, '_> {
    pub fn points_m(&self) -> &[Point] {
        self.rows.points_m()
    }
    pub fn kernel_evaluations(&self) -> usize {
        self.rows.kernel_evaluations()
    }

    /// Total physical p/u for the selected incident and actual source rates.
    /// This evaluates no Green kernels and repeats no boundary factorization.
    pub fn evaluate(
        &self,
        case: usize,
        rates: &[C64],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<near_field::VectorObservation, HelmholtzError> {
        poll(&mut cancelled)?;
        let jump = self
            .source
            .total_pressure_jump(case, rates, &mut cancelled)?;
        let mut total = self.rows.evaluate_jump(&jump, &mut cancelled)?;
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
                    .any(|&v| !finite(v))
            {
                return Err(bad("molecular illuminated total receiver overflows"));
            }
        }
        poll(&mut cancelled)?;
        Ok(total)
    }
}
