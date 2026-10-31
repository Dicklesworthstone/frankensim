//! Outgoing angular fields of the retained material-boundary response.
//!
//! Finite radiator solutions already contain the complete outgoing field.
//! Point-source solutions contain scattered traces, so their analytic outgoing
//! volume monopole or force dipole is added exactly once. Complex source drives are summed before
//! angular projection; summing individual source powers would lose interference.

use super::IterativeBoundaryResponse;
use crate::helmholtz::{
    self, DirectivityProjectionOptions, DirectivityTable, HelmholtzError, bc::IncidentSource,
};
use crate::observation_transfers::{IterativeObservationError, SourceId, TransferError, cap};
use fs_blake3::ContentHash;
use fs_exec::CancelGate;
use fs_math::{c64::C64, det};
use fs_vfit::observation::MAX_OBSERVATION_ENTRIES;

/// Source-major outgoing amplitudes with `p ~ F exp(ikr)/r`.
/// Each column retains its source-coordinate units: multiply by that source's
/// complex drive (including flow [m³/s] or force [N]) before combining fields.
/// These are not finite-distance pressures.
#[derive(Clone, Debug)]
pub struct IterativeFarFieldSample {
    pub omega_rad_s: f64,
    pub operator_hash: ContentHash,
    pub sources: Vec<SourceId>,
    pub directions: Vec<[f64; 3]>,
    pub amplitudes: Vec<C64>,
    /// Triangle/centroid phase evaluations plus direct compact-source phases.
    pub kernel_evaluations: usize,
}

/// Independent angular, kernel-work and temporary-storage budgets.
#[derive(Clone, Copy, Debug)]
pub struct IterativeDirectivityOptions {
    pub l_max: usize,
    pub angular: DirectivityProjectionOptions,
    /// Total phase evaluations across every angular refinement, not per batch.
    pub maximum_kernel_evaluations: usize,
    /// Complex entries for the combined pressure AND velocity boundary arrays.
    pub maximum_combined_entries: usize,
}

/// A coherent driven field projected onto the existing complex SH basis.
#[derive(Clone, Debug)]
pub struct IterativeCoherentDirectivity {
    pub omega_rad_s: f64,
    pub operator_hash: ContentHash,
    /// Complex drives in the response's original source order and SI units.
    pub drives: Vec<C64>,
    /// `captured_fraction` measures omitted SH power against independently
    /// refined angular power. BEM/FMM and panel quadrature errors remain separate.
    pub table: DirectivityTable,
    pub kernel_evaluations: usize,
}

impl IterativeBoundaryResponse<'_, '_> {
    /// Observe every retained source in explicit unit directions, without
    /// another boundary solve. Work is admitted for the entire result before
    /// its first evaluation; cancellation never returns a partial source table.
    pub fn far_field(
        &self,
        directions: &[[f64; 3]],
        maximum_kernel_evaluations: usize,
        gate: &CancelGate,
    ) -> Result<IterativeFarFieldSample, IterativeObservationError> {
        poll(gate)?;
        if directions.is_empty() {
            return Err(invalid("far-field directions must not be empty"));
        }
        for (index, &direction) in directions.iter().enumerate() {
            if index % 256 == 0 {
                poll(gate)?;
            }
            let norm = det::sqrt(direction.iter().map(|v| v * v).sum::<f64>());
            if !norm.is_finite() || (norm - 1.0).abs() > 1e-10 {
                return Err(invalid("far-field directions must be finite unit vectors"));
            }
        }
        let entries = directions
            .len()
            .checked_mul(self.sources().len())
            .ok_or_else(|| invalid("far-field output shape overflow"))?;
        cap(entries, MAX_OBSERVATION_ENTRIES)?;
        let panel_work = panel_work(self, directions.len())?;
        let mut work = KernelWork::new(maximum_kernel_evaluations);
        for index in 0..self.sources().len() {
            work.charge(panel_work)?;
            if self.source_incident(index)?.is_some() {
                work.charge(directions.len())?;
            }
        }
        let mut amplitudes = reserve(entries)?;
        for (index, solution) in self.solutions().iter().enumerate() {
            poll(gate)?;
            let mut field = helmholtz::far_field_integrated_cancellable(
                self.surface(),
                solution,
                self.medium(),
                directions,
                panel_work,
                || gate.is_requested(),
            )?;
            if let Some(source) = self.source_incident(index)? {
                add_compact_sources(
                    &mut field,
                    directions,
                    &[source],
                    solution.k,
                    self.omega_rad_s() * self.medium().density,
                    gate,
                )?;
            }
            amplitudes.extend(field);
        }
        poll(gate)?;
        Ok(IterativeFarFieldSample {
            omega_rad_s: self.omega_rad_s(),
            operator_hash: self.operator_hash(),
            sources: self.sources().to_vec(),
            directions: directions.to_vec(),
            amplitudes,
            kernel_evaluations: work.used,
        })
    }

    /// Project the actual complex superposition of source fields. Combining
    /// boundary traces once costs O(sources * panels); each angular evaluation
    /// then visits the surface once, independent of the number of radiators.
    /// The existing angular refinement resolves total field variation even
    /// when it lies outside the requested SH degree. Neither interference
    /// nulls nor negative mutual terms are clipped or replaced by source powers.
    pub fn coherent_directivity(
        &self,
        drives: &[C64],
        options: IterativeDirectivityOptions,
        gate: &CancelGate,
    ) -> Result<IterativeCoherentDirectivity, IterativeObservationError> {
        poll(gate)?;
        if drives.len() != self.sources().len() || drives.iter().any(|&v| !finite(v)) {
            return Err(invalid(
                "coherent directivity requires one finite drive per source",
            ));
        }
        let entries = self
            .surface()
            .areas()
            .len()
            .checked_mul(2)
            .ok_or_else(|| invalid("combined far-field trace shape overflow"))?;
        cap(
            entries,
            options
                .maximum_combined_entries
                .min(MAX_OBSERVATION_ENTRIES),
        )?;
        let first = self
            .solutions()
            .first()
            .ok_or_else(|| invalid("coherent directivity requires a retained source"))?;
        // Only this private pressure/velocity image is used by the far-field
        // owner. Individual solve/power diagnostics stay on the original response.
        let mut combined = first.clone();
        for (panel, (pressure, velocity)) in combined
            .pressure
            .iter_mut()
            .zip(&mut combined.velocity)
            .enumerate()
        {
            if panel % 256 == 0 {
                poll(gate)?;
            }
            *pressure = *pressure * drives[0];
            *velocity = *velocity * drives[0];
        }
        for (solution, &drive) in self.solutions().iter().zip(drives).skip(1) {
            for panel in 0..combined.pressure.len() {
                if panel % 256 == 0 {
                    poll(gate)?;
                }
                combined.pressure[panel] =
                    combined.pressure[panel] + solution.pressure[panel] * drive;
                combined.velocity[panel] =
                    combined.velocity[panel] + solution.velocity[panel] * drive;
            }
        }
        if combined
            .pressure
            .iter()
            .chain(&combined.velocity)
            .any(|&v| !finite(v))
        {
            return Err(invalid(
                "coherent boundary trace exceeds finite numerical range",
            ));
        }
        let mut direct_sources = Vec::new();
        for (index, &drive) in drives.iter().enumerate() {
            if let Some(mut source) = self.source_incident(index)? {
                match &mut source {
                    IncidentSource::Monopole(source) => {
                        source.volume_velocity_m3_s = source.volume_velocity_m3_s * drive;
                    }
                    IncidentSource::ForceDipole(source) => {
                        for component in &mut source.force_on_fluid_n {
                            *component = *component * drive;
                        }
                    }
                    IncidentSource::PlaneWave(_) => {
                        return Err(invalid("a plane wave has no outgoing total far field"));
                    }
                }
                direct_sources.push(source);
            }
        }
        let mut work = KernelWork::new(options.maximum_kernel_evaluations);
        let table = helmholtz::directivity_sh_table_from_evaluator_with_options(
            combined.k,
            options.l_max,
            options.angular,
            |directions| {
                let boundary_work = panel_work(self, directions.len())?;
                let direct_work = directions
                    .len()
                    .checked_mul(direct_sources.len())
                    .ok_or_else(|| budget_overflow(work.maximum))?;
                work.charge(boundary_work)?;
                work.charge(direct_work)?;
                let mut field = helmholtz::far_field_integrated_cancellable(
                    self.surface(),
                    &combined,
                    self.medium(),
                    directions,
                    boundary_work,
                    || gate.is_requested(),
                )?;
                add_compact_sources(
                    &mut field,
                    directions,
                    &direct_sources,
                    combined.k,
                    self.omega_rad_s() * self.medium().density,
                    gate,
                )?;
                Ok(field)
            },
            || gate.is_requested(),
        )?;
        poll(gate)?;
        Ok(IterativeCoherentDirectivity {
            omega_rad_s: self.omega_rad_s(),
            operator_hash: self.operator_hash(),
            drives: drives.to_vec(),
            table,
            kernel_evaluations: work.used,
        })
    }
}

fn panel_work(
    response: &IterativeBoundaryResponse<'_, '_>,
    directions: usize,
) -> Result<usize, HelmholtzError> {
    let nodes = if response.surface().triangles().is_some() {
        16
    } else {
        1
    };
    directions
        .checked_mul(response.surface().areas().len())
        .and_then(|value| value.checked_mul(nodes))
        .ok_or_else(|| budget_overflow(usize::MAX))
}

fn add_compact_sources(
    field: &mut [C64],
    directions: &[[f64; 3]],
    sources: &[IncidentSource],
    k: f64,
    omega_rho: f64,
    gate: &CancelGate,
) -> Result<(), HelmholtzError> {
    let flow_factor = C64::new(0.0, -omega_rho / (4.0 * core::f64::consts::PI));
    let force_factor = C64::new(0.0, -k / (4.0 * core::f64::consts::PI));
    for (value, direction) in field.iter_mut().zip(directions) {
        poll(gate)?;
        let direction = helmholtz::far_field_unit_direction(*direction);
        for source in sources {
            let (position, amplitude) = match source {
                IncidentSource::Monopole(source) => {
                    (source.position_m, flow_factor * source.volume_velocity_m3_s)
                }
                IncidentSource::ForceDipole(source) => {
                    let axial = (0..3).fold(C64::ZERO, |sum, axis| {
                        sum + source.force_on_fluid_n[axis].scale(direction[axis])
                    });
                    (source.position_m, force_factor * axial)
                }
                IncidentSource::PlaneWave(_) => {
                    return Err(HelmholtzError::BadParameter {
                        what: "a plane wave has no outgoing total far field",
                    });
                }
            };
            let phase = -k
                * direction
                    .iter()
                    .zip(position)
                    .map(|(a, b)| a * b)
                    .sum::<f64>();
            *value = *value + amplitude * C64::new(det::cos(phase), det::sin(phase));
        }
        if !finite(*value) {
            return Err(HelmholtzError::BadParameter {
                what: "total outgoing field exceeds finite numerical range",
            });
        }
    }
    Ok(())
}

struct KernelWork {
    used: usize,
    maximum: usize,
}
impl KernelWork {
    fn new(maximum: usize) -> Self {
        Self { used: 0, maximum }
    }
    fn charge(&mut self, count: usize) -> Result<(), HelmholtzError> {
        let required = self.used.checked_add(count);
        if required.is_none_or(|value| value > self.maximum) {
            return Err(HelmholtzError::FarFieldBudget {
                required_evaluations: required.unwrap_or(usize::MAX),
                maximum_evaluations: self.maximum,
            });
        }
        self.used = required.unwrap();
        Ok(())
    }
}
fn budget_overflow(maximum_evaluations: usize) -> HelmholtzError {
    HelmholtzError::FarFieldBudget {
        required_evaluations: usize::MAX,
        maximum_evaluations,
    }
}
fn reserve<T>(count: usize) -> Result<Vec<T>, HelmholtzError> {
    let mut output = Vec::new();
    output
        .try_reserve_exact(count)
        .map_err(|_| HelmholtzError::AllocationFailed {
            what: "iterative far-field samples",
        })?;
    Ok(output)
}
fn invalid(what: &'static str) -> IterativeObservationError {
    TransferError::Invalid(what).into()
}
fn finite(value: C64) -> bool {
    value.re.is_finite() && value.im.is_finite()
}
fn poll(gate: &CancelGate) -> Result<(), HelmholtzError> {
    if gate.is_requested() {
        Err(HelmholtzError::Cancelled)
    } else {
        Ok(())
    }
}

#[cfg(test)]
mod tests;
