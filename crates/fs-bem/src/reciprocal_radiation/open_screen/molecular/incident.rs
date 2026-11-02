//! One pressure-jump solve bank for shell motion and coherent illumination.
//!
//! The incident field is regular across a zero-thickness screen and has no
//! pressure jump of its own. A blocked screen therefore solves
//! H*p_blocked = -i*rho*omega*integral(phi*u_incident.n). Its force on the
//! original structural coordinates is -G*p_blocked. Subsequent motion adds
//! the original unit-source jump columns, never a duplicated incident load.
use super::*;
use crate::helmholtz::{bc::IncidentField, molecular::PreparedMolecularIncident};

mod receivers;
pub use receivers::PreparedMolecularScreenIlluminationReceivers;

/// Full illumination limits in addition to the original mechanical source.
#[derive(Debug, Clone, Copy)]
pub struct MolecularScreenIlluminationLimits {
    /// Port/motion jump limits retain their original meaning. Its total work
    /// ceiling now admits BOTH motion and blocked columns with one H/LU.
    pub boundary: MolecularScreenPortLimits,
    /// Complete independent incident cases, in 1..=256; no case is truncated.
    pub maximum_incident_cases: usize,
    /// Every retained blocked nodal jump, cases*original_vertices.
    pub maximum_incident_jump_entries: usize,
    /// All analytic source terms times (cached quadrature points + facets),
    /// admitting both original sampled drives and finite-emitter clearance.
    /// Each later receiver preparation uses the same ceiling for terms*points.
    pub maximum_source_point_products: u64,
    /// Exact distance from every finite emitter to each triangle face/edge.
    /// Either side is admitted; this is not a closed-body exterior test.
    pub minimum_source_clearance_m: f64,
    /// Full logical live payload: max(original solve payload, declared operator
    /// workspace) plus 16*c*(vertices+ports)+64*n+32*(vertices+ports)+24*t
    /// +256*c+4096. n is active pressure count, c cases, t all source terms.
    /// This includes retained incident owners, blocked jumps/forces, source
    /// position scratch, complex weak drive and later total-jump/work scratch.
    /// Geometry ownership and separately admitted receiver rows/outputs,
    /// allocator metadata and thread stacks are excluded.
    pub maximum_live_bytes: usize,
}

/// Original motion bank, blocked sound fields and exact mechanical projection.
/// Incident declarations stay borrowed and the evaluated gas stays immutable.
pub struct MolecularScreenIllumination<'p, 's, 'i> {
    motion: MolecularScreenVelocityPortSolution<'p, 's>,
    incidents: Vec<PreparedMolecularIncident<'i>>,
    blocked: Vec<Vec<C64>>,
    forces: Vec<Vec<C64>>,
    pairing: Vec<Vec<f64>>,
    blocked_residuals: Vec<f64>,
    frequency_work: u64,
    maximum_source_point_products: u64,
}

impl<'s> PreparedScreen<'s> {
    /// Illuminate both sides while retaining all original structural unit
    /// sources. ONE original H and LU serve every motion and incident column.
    /// Incident normal velocity is integrated at the original cached P1
    /// triangle quadrature points; it is not sampled only at centroids or
    /// clamped at free-edge vertices. That declared quadrature/mesh must be
    /// refined for close or rapidly varying sources independently of residuals.
    ///
    /// Analytic incident pressure is continuous across this single-jump sheet.
    /// Blocked pressure jump supplies -G*p_blocked force on the structure.
    /// The companion structural owner must add that force before solving its
    /// loaded motion. Receiver evaluation then adds incident sound once to the
    /// coherent blocked-plus-motion double-layer field.
    #[allow(clippy::too_many_arguments)]
    pub fn solve_molecular_illumination<'p, 'i>(
        &'p self,
        source: &MolecularExterior,
        omega: f64,
        rows: &[&[[f64; 3]]],
        incidents: &[IncidentField<'i>],
        limits: MolecularScreenIlluminationLimits,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<MolecularScreenIllumination<'p, 's, 'i>, MolecularError> {
        poll(&mut cancelled)?;
        let cases = incidents.len();
        if cases == 0
            || cases > limits.maximum_incident_cases
            || limits.maximum_incident_cases == 0
            || limits.maximum_incident_cases > crate::helmholtz::MAX_RADIATION_FIELDS_PER_BATCH
            || !limits.minimum_source_clearance_m.is_finite()
            || limits.minimum_source_clearance_m <= 0.0
            || limits.maximum_source_point_products == 0
        {
            return Err(bad("molecular screen illumination needs bounded cases and positive source clearance/work").into());
        }
        let terms: u128 = incidents
            .iter()
            .map(|field| match field {
                IncidentField::PlaneWave(_) => 1,
                IncidentField::Monopoles(sources) => sources.len() as u128,
                IncidentField::ForceDipoles(sources) => sources.len() as u128,
                IncidentField::Sources(sources) => sources.len() as u128,
            })
            .sum();
        let (c, p, n, v, f) = (
            cases as u128,
            rows.len() as u128,
            self.n as u128,
            self.surface.vertices.len() as u128,
            self.facets.len() as u128,
        );
        let samples: u128 = self.samples.iter().map(|row| row.len() as u128).sum();
        let base_work = self.molecular_velocity_port_work(rows.len())? as u128;
        let total_work =
            base_work + 128 * c * n * n + 256 * c * samples + 64 * c * v + 64 * c * p * v;
        let base_bytes = self
            .molecular_port_payload(rows.len())?
            .max(limits.boundary.operator.maximum_workspace_bytes as u128);
        let bytes =
            base_bytes + 16 * c * (v + p) + 64 * n + 32 * (v + p) + 24 * terms + 256 * c + 4096;
        if c * v > limits.maximum_incident_jump_entries as u128
            || total_work > limits.boundary.ports.maximum_work as u128
            || total_work > self.options.maximum_work as u128
            || terms * (samples + f) > limits.maximum_source_point_products as u128
            || bytes > limits.maximum_live_bytes as u128
        {
            return Err(
                bad("complete molecular screen illumination exceeds declared resources").into(),
            );
        }
        let frequency_work = u64::try_from(total_work)
            .map_err(|_| bad("molecular screen illumination work overflows"))?;
        let mut prepared = Vec::new();
        prepared
            .try_reserve_exact(cases)
            .map_err(|_| HelmholtzError::AllocationFailed {
                what: "molecular screen incident owners",
            })?;
        let first = source.prepare_incident(omega, incidents[0], &mut cancelled)?;
        for &field in &incidents[1..] {
            prepared.push(first.with_field(field, &mut cancelled)?);
        }
        prepared.insert(0, first);
        // No source enclosing sphere or winding rule is valid for an open
        // screen: use its original triangles, including their edges.
        let mut scale = 1.0_f64;
        for (vertex, point) in self.surface.vertices.iter().enumerate() {
            if vertex % 256 == 0 {
                poll(&mut cancelled)?;
            }
            for value in point {
                scale = scale.max(value.abs());
            }
        }
        for field in &prepared {
            poll(&mut cancelled)?;
            for point in field.source_positions() {
                for value in point {
                    scale = scale.max(value.abs());
                }
                if limits.minimum_source_clearance_m <= 64.0 * f64::EPSILON * scale {
                    return Err(bad(
                        "molecular screen source clearance is below coordinate resolution",
                    )
                    .into());
                }
                for (index, facet) in self.facets.iter().enumerate() {
                    if index % 32 == 0 {
                        poll(&mut cancelled)?;
                    }
                    let distance = crate::near_field::distance(point, facet.triangle);
                    if !distance.is_finite() || distance <= limits.minimum_source_clearance_m {
                        return Err(
                            bad("molecular screen finite source touches its face or edge").into(),
                        );
                    }
                }
            }
        }
        let (motion, (blocked, forces, blocked_residuals, pairing)) = self
            .solve_molecular_ports_with(
                source,
                omega,
                rows,
                limits.boundary,
                Some(&prepared[0]),
                |pairing, h, lu, rho_omega, cancelled| {
                    let mut blocked = Vec::with_capacity(cases);
                    let mut forces = Vec::with_capacity(cases);
                    let mut residuals = Vec::with_capacity(cases);
                    for incident in &prepared {
                        poll(cancelled)?;
                        let mut weak: Vec<Sum> = (0..self.n).map(|_| Sum::default()).collect();
                        for (facet, samples) in self.facets.iter().zip(&self.samples) {
                            poll(cancelled)?;
                            for sample in samples {
                                let field =
                                    incident.pressure_velocity_at(sample.point, &mut *cancelled)?;
                                let mut normal = Sum::default();
                                for axis in 0..3 {
                                    normal.add(field.velocity_m_s[axis].scale(facet.normal[axis]));
                                }
                                let drive = -normal.value();
                                for a in 0..3 {
                                    if let Some(u) = self.dof[facet.ids[a]] {
                                        weak[u].add(drive.scale(sample.weight * sample.shape[a]));
                                    }
                                }
                            }
                        }
                        let mut rhs = Vec::with_capacity(self.n);
                        for (i, value) in weak.into_iter().enumerate() {
                            if i % 128 == 0 {
                                poll(cancelled)?;
                            }
                            let value = value.value();
                            let force = C64::new(-rho_omega * value.im, rho_omega * value.re);
                            if !finite(force)
                                || (value.re != 0.0 && force.im == 0.0)
                                || (value.im != 0.0 && force.re == 0.0)
                            {
                                return Err(bad(
                                    "molecular blocked physical drive is unrepresentable",
                                )
                                .into());
                            }
                            rhs.push(force);
                        }
                        let (jump, residual) = self.solve_molecular_rhs(&rhs, h, lu, cancelled)?;
                        let mut force = Vec::with_capacity(pairing.len());
                        for row in &pairing {
                            poll(cancelled)?;
                            force.push(-self.molecular_jump_effort(row, &jump, cancelled)?);
                        }
                        blocked.push(jump);
                        forces.push(force);
                        residuals.push(residual);
                    }
                    Ok((blocked, forces, residuals, pairing))
                },
                &mut cancelled,
            )?;
        poll(&mut cancelled)?;
        Ok(MolecularScreenIllumination {
            motion,
            incidents: prepared,
            blocked,
            forces,
            pairing,
            blocked_residuals,
            frequency_work,
            maximum_source_point_products: limits.maximum_source_point_products,
        })
    }
}

impl<'p, 's, 'i> MolecularScreenIllumination<'p, 's, 'i> {
    pub fn motion(&self) -> &MolecularScreenVelocityPortSolution<'p, 's> {
        &self.motion
    }
    pub fn blocked_forces(&self) -> &[Vec<C64>] {
        &self.forces
    }
    pub fn blocked_pressure_jumps_pa(&self) -> &[Vec<C64>] {
        &self.blocked
    }
    pub fn blocked_backward_residuals(&self) -> &[f64] {
        &self.blocked_residuals
    }
    pub fn frequency_work(&self) -> u64 {
        self.frequency_work
    }
    pub fn maximum_backward_residual(&self) -> f64 {
        self.blocked_residuals
            .iter()
            .copied()
            .fold(self.motion.maximum_backward_residual(), f64::max)
    }
    fn admit_case(&self, case: usize) -> Result<(), HelmholtzError> {
        if case >= self.blocked.len() {
            Err(bad("molecular screen incident case is out of range"))
        } else {
            Ok(())
        }
    }
    /// Total original jump; the analytic incident field has no sheet jump.
    /// This combines blocked and actual motion fields before observation/work.
    pub fn total_pressure_jump(
        &self,
        case: usize,
        rates: &[C64],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Vec<C64>, HelmholtzError> {
        poll(&mut cancelled)?;
        self.admit_case(case)?;
        let mut jump = self.motion.pressure_jump(rates, &mut cancelled)?;
        for (vertex, (value, &blocked)) in jump.iter_mut().zip(&self.blocked[case]).enumerate() {
            if vertex % 128 == 0 {
                poll(&mut cancelled)?;
            }
            *value = *value + blocked;
            if !finite(*value) {
                return Err(bad("molecular screen total pressure jump overflows"));
            }
        }
        poll(&mut cancelled)?;
        Ok(jump)
    }
    /// Direct -G*total_jump projection, in the original source coordinates.
    /// Its signed work can be negative: incident sound can drive the structure.
    /// Temporary total jump and modal output fit the admitted boundary payload.
    pub fn total_fluid_reaction(
        &self,
        case: usize,
        rates: &[C64],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Vec<C64>, HelmholtzError> {
        let jump = self.total_pressure_jump(case, rates, &mut cancelled)?;
        let mut reaction = Vec::with_capacity(self.pairing.len());
        for row in &self.pairing {
            poll(&mut cancelled)?;
            reaction.push(-self.motion.prepared.molecular_jump_effort(
                row,
                &jump,
                &mut cancelled,
            )?);
        }
        poll(&mut cancelled)?;
        Ok(reaction)
    }
}

#[cfg(test)]
mod tests;
