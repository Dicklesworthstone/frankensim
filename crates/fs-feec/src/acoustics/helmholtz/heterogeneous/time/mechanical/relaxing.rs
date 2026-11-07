//! Original spatial admittance memories shared by the coupled wall integrators.

use super::super::super::SpatialRelaxingBoundary;
use super::{
    BoundaryAdmittance, C64, Coo, Csr, Dd, HeterogeneousError, HeterogeneousHelmholtzFe,
    energy_change, invalid, paired_work, poll, positive_form, relative,
};

/// Raw pressure memory and its original geometric trace; no residue inverse.
pub(super) struct TraceMemory {
    pub(super) patch: usize,
    pub(super) term: usize,
    pub(super) dofs: Vec<usize>,
    pub(super) conductance: Csr,
    pub(super) tau: f64,
}

pub(super) fn prepare(
    model: &HeterogeneousHelmholtzFe<'_, '_>,
    laws: &[&SpatialRelaxingBoundary<'_, '_, '_>],
    max_memory_dofs: usize,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Vec<TraceMemory>, HeterogeneousError> {
    let mut arms = Vec::new();
    let mut count = 0usize;
    for (patch, law) in laws.iter().enumerate() {
        poll(cancelled)?;
        law.instantaneous.validate(model)?;
        // Geometric support remains live even where every residue is zero.
        let (trace, _) = model.admittance_matrices(
            &[BoundaryAdmittance {
                faces: law.faces().to_vec(),
                admittance_m_pa_s: C64::ONE,
            }],
            cancelled,
        )?;
        let dofs = (0..model.ndof())
            .filter(|&row| {
                let (columns, values) = trace.row(row);
                columns.binary_search(&row).is_ok_and(|i| values[i] > 0.0)
            })
            .collect::<Vec<_>>();
        count = dofs
            .len()
            .checked_mul(law.terms.len())
            .and_then(|n| count.checked_add(n))
            .ok_or_else(|| invalid("coupled spatial wall memory extent overflow"))?;
        if count > max_memory_dofs {
            return Err(HeterogeneousError::Limit {
                what: "coupled spatial admittance memory degrees of freedom",
                requested: count,
                maximum: max_memory_dofs,
            });
        }
        for (term, original) in law.terms.iter().enumerate() {
            poll(cancelled)?;
            let mut local = Coo::new(dofs.len(), dofs.len());
            for (row, &global) in dofs.iter().enumerate() {
                poll(cancelled)?;
                let (columns, values) = original.conductance.row(global);
                for (&column, &value) in columns.iter().zip(values) {
                    if value != 0.0 {
                        let column = dofs
                            .binary_search(&column)
                            .map_err(|_| invalid("incomplete coupled admittance trace support"))?;
                        local.push(row, column, value);
                    }
                }
            }
            arms.push(TraceMemory {
                patch,
                term,
                dofs: dofs.clone(),
                conductance: local.assemble(),
                tau: original.relaxation_time_s,
            });
        }
    }
    Ok(arms)
}

pub(super) fn energy(
    arms: &[TraceMemory],
    memories: &[Vec<f64>],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(f64, f64), HeterogeneousError> {
    let mut stored = Dd::from_f64(0.0);
    let mut roundoff = 0.0;
    for (arm, memory) in arms.iter().zip(memories) {
        let (value, error) = positive_form(&arm.conductance, memory, cancelled)?;
        stored = stored + Dd::from_f64(0.5) * Dd::from_f64(arm.tau) * Dd::from_f64(value);
        roundoff += 0.5 * arm.tau * error;
    }
    if !stored.to_f64().is_finite() || !roundoff.is_finite() {
        return Err(invalid("coupled admittance stored energy overflow").into());
    }
    Ok((stored.to_f64(), roundoff))
}

pub(super) fn subtract_history(
    arms: &[TraceMemory],
    weights: &[(f64, f64)],
    memories: &[Vec<f64>],
    load: &mut [f64],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(), HeterogeneousError> {
    for ((arm, &(history, _)), memory) in arms.iter().zip(weights).zip(memories) {
        poll(cancelled)?;
        let mut projected = vec![0.0; arm.dofs.len()];
        arm.conductance.spmv(memory, &mut projected);
        for (&dof, value) in arm.dofs.iter().zip(projected) {
            load[dof] =
                (Dd::from_f64(load[dof]) - Dd::from_f64(history) * Dd::from_f64(value)).to_f64();
            if !load[dof].is_finite() {
                return Err(invalid("coupled admittance history load overflow").into());
            }
        }
    }
    Ok(())
}

#[derive(Default)]
pub(super) struct MidpointUpdate {
    pub(super) memories: Vec<Vec<f64>>,
    pub(super) outward: Vec<f64>,
    pub(super) energy_change: f64,
    pub(super) heat: f64,
    pub(super) work: f64,
    pub(super) heat_roundoff: f64,
    pub(super) work_roundoff: f64,
    pub(super) backward_error: f64,
}

pub(super) fn midpoint(
    arms: &[TraceMemory],
    weights: &[(f64, f64)],
    initial: &[Vec<f64>],
    pressure: &[f64],
    duration_s: f64,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<MidpointUpdate, HeterogeneousError> {
    let mut result = MidpointUpdate::default();
    if arms.is_empty() {
        return Ok(result);
    }
    result.outward = vec![0.0; pressure.len()];
    result.memories.reserve(arms.len());
    let mut change = Dd::from_f64(0.0);
    let mut heat = Dd::from_f64(0.0);
    let mut work = Dd::from_f64(0.0);
    let h = Dd::from_f64(duration_s);
    for ((arm, &(history, driven)), q0) in arms.iter().zip(weights).zip(initial) {
        poll(cancelled)?;
        let mut endpoint = Vec::with_capacity(arm.dofs.len());
        let mut middle = Vec::with_capacity(arm.dofs.len());
        let mut trace = Vec::with_capacity(arm.dofs.len());
        for (&dof, &before) in arm.dofs.iter().zip(q0) {
            poll(cancelled)?;
            let p = pressure[dof];
            let qm = (Dd::from_f64(history) * Dd::from_f64(before)
                + Dd::from_f64(driven) * Dd::from_f64(p))
            .to_f64();
            let after = (Dd::from_f64(2.0) * Dd::from_f64(qm) - Dd::from_f64(before)).to_f64();
            let defect = (Dd::from_f64(arm.tau) * (Dd::from_f64(after) - Dd::from_f64(before))
                + h * (Dd::from_f64(qm) - Dd::from_f64(p)))
            .to_f64();
            let scale = arm.tau * (before.abs() + after.abs()) + duration_s * (qm.abs() + p.abs());
            if ![qm, after, defect, scale].into_iter().all(f64::is_finite) {
                return Err(invalid("coupled admittance memory candidate overflow").into());
            }
            result.backward_error = result.backward_error.max(relative(defect.abs(), scale));
            endpoint.push(after);
            middle.push(qm);
            trace.push(p);
        }
        let mut projected = vec![0.0; arm.dofs.len()];
        arm.conductance.spmv(&middle, &mut projected);
        for (&dof, &value) in arm.dofs.iter().zip(&projected) {
            result.outward[dof] =
                (Dd::from_f64(result.outward[dof]) + Dd::from_f64(value)).to_f64();
        }
        let (value, error) = paired_work(duration_s, &trace, &projected);
        work = work + Dd::from_f64(value);
        result.work_roundoff += error;
        let (value, error) = positive_form(&arm.conductance, &middle, cancelled)?;
        heat = heat + h * Dd::from_f64(value);
        result.heat_roundoff += duration_s * error;
        change = change
            + Dd::from_f64(arm.tau)
                * Dd::from_f64(energy_change(&arm.conductance, q0, &endpoint, cancelled)?);
        result.memories.push(endpoint);
    }
    result.energy_change = change.to_f64();
    result.heat = heat.to_f64();
    result.work = work.to_f64();
    Ok(result)
}
