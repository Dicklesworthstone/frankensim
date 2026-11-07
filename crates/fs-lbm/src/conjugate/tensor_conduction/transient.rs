//! Backward-Euler storage on the same conservative full-tensor face operator.

use fs_exec::CancelGate;

use super::{ConductivityTensor, TensorConductionSolution, solve_tensor_conduction_inner};
use crate::conjugate::energy::{EnergyConfig, PseudoStep, ThermalSetup};
use crate::conjugate::transient::{TransientRecord, record_conduction_step};
use crate::conjugate::{ChtError, SolidMaterial, VoxelDomain, finite, finite_positive, poll};

/// One accepted interval, retaining its actual tensor face temperatures and heat.
#[derive(Debug, Clone, PartialEq)]
pub struct TensorConductionStep {
    /// Endpoint temperature, physical face traces/fluxes and sparse-solve evidence.
    pub solution: TensorConductionSolution,
    /// Independently computed storage/source/reservoir heat ledger in joules.
    pub record: TransientRecord,
}

/// Advance frozen full-tensor conduction with one physical backward-Euler step.
///
/// `heat_capacity_j_k` contains the complete, strictly positive sensible heat
/// capacity C for each cell, in J/K. It stays fixed over this interval; it is
/// not a volumetric density and is not multiplied by volume again. For an
/// ordinary material it is rho*c*V. No phase change, inferred latent capacity,
/// material selection or time-varying tensor law is introduced by this entry.
///
/// The existing condensed face operator receives only the cell terms
/// C/dt and (C/dt)*T_old. Faces and zero-thickness contacts store no energy.
/// The immutable previous field is the physical start of the interval, never
/// a solver trial. The supplied setup is the endpoint source/boundary sample.
/// Callers may compose intervals and change that sample at known forcing times.
///
/// Positive storage also admits an entirely insulated body. A prescribed cell
/// retains its declared endpoint temperature; its external reservoir supplies
/// both its conductive reaction and C*(T_new-T_old), exactly once. The record
/// uses the same physical heat ledger as the native conjugate transient owner.
///
/// Backward Euler is first order in time. Full SPD tensor diffusion is
/// dissipative; no discrete maximum principle for arbitrary rotated anisotropy
/// is claimed. Native true-residual and local heat-balance checks include
/// storage, while the retained face heat contains conduction alone.
///
/// # Errors
/// The steady tensor admission rules still apply, including its supported
/// boundary/contact scope. Mismatched fields, nonpositive capacities or Kelvin
/// temperatures, unrepresentable time/storage, allocation failure, an
/// unconverged endpoint, and cancellation publish no partial step.
#[allow(clippy::too_many_arguments)]
pub fn advance_tensor_conduction(
    domain: &VoxelDomain,
    solids: &[SolidMaterial],
    setup: &ThermalSetup,
    tensors: &[ConductivityTensor],
    heat_capacity_j_k: &[f64],
    previous_temperature_k: &[f64],
    start_time_s: f64,
    time_step_s: f64,
    config: &EnergyConfig,
    gate: &CancelGate,
) -> Result<TensorConductionStep, ChtError> {
    let prepared = prepare_tensor_step(
        domain,
        tensors.len(),
        heat_capacity_j_k,
        previous_temperature_k,
        start_time_s,
        time_step_s,
        gate,
    )?;
    let step = prepared.storage();
    let solution =
        solve_tensor_conduction_inner(domain, solids, setup, tensors, config, Some(&step), gate)?;
    poll(gate)?;
    let record = prepared.record(domain, setup, &solution)?;
    poll(gate)?;
    Ok(TensorConductionStep { solution, record })
}

/// One admitted physical sensible-storage interval. Frozen and original-law
/// tensor drivers borrow the same immutable old state and use the same ledger.
pub(super) struct PreparedTensorStep<'a> {
    coefficient_w_k: Vec<f64>,
    heat_capacity_j_k: &'a [f64],
    previous_temperature_k: &'a [f64],
    time_s: f64,
    time_step_s: f64,
}

impl PreparedTensorStep<'_> {
    /// Borrow cell-only storage; no trace or contact gains artificial capacity.
    pub(super) fn storage(&self) -> PseudoStep<'_> {
        PseudoStep {
            coefficient_w_k: &self.coefficient_w_k,
            previous: self.previous_temperature_k,
            junction_coefficient_w_k: &[],
        }
    }

    /// Account actual final heat and original-state storage, including the
    /// prescribed-cell thermostat contribution through the existing owner.
    pub(super) fn record(
        &self,
        domain: &VoxelDomain,
        setup: &ThermalSetup,
        solution: &TensorConductionSolution,
    ) -> Result<TransientRecord, ChtError> {
        let solids_present = (0..domain.cell_count()).any(|cell| !domain.is_fluid(cell));
        record_conduction_step(
            domain,
            setup,
            &solution.temperature,
            solution.report.balance,
            solution.report.iterations,
            |cell| {
                self.heat_capacity_j_k[cell]
                    * (solution.temperature[cell] - self.previous_temperature_k[cell])
            },
            (self.time_s, self.time_step_s, solids_present),
            0.0,
            0.0,
        )
    }
}

/// Admit the same physical time, field shape and C/dt for either tensor owner.
/// Conductivity entries are frozen tensors or immutable original laws.
pub(super) fn prepare_tensor_step<'a>(
    domain: &VoxelDomain,
    conductivity_entries: usize,
    heat_capacity_j_k: &'a [f64],
    previous_temperature_k: &'a [f64],
    start_time_s: f64,
    time_step_s: f64,
    gate: &CancelGate,
) -> Result<PreparedTensorStep<'a>, ChtError> {
    poll(gate)?;
    finite("tensor.transient.start_time_s", start_time_s)?;
    finite_positive("tensor.transient.time_step_s", time_step_s)?;
    let time_s = start_time_s + time_step_s;
    finite("tensor.transient.time_s", time_s)?;
    if start_time_s < 0.0 || time_s <= start_time_s {
        return Err(ChtError::InvalidInput {
            field: "tensor.transient.time",
            reason: "a nonnegative start and representable forward interval are required".into(),
        });
    }
    let cells = domain.cell_count();
    if heat_capacity_j_k.len() != cells
        || previous_temperature_k.len() != cells
        || conductivity_entries != cells
    {
        return Err(ChtError::InvalidInput {
            field: "tensor.transient.fields",
            reason: "one capacity, previous temperature and conductivity entry per domain cell are required"
                .into(),
        });
    }
    cells
        .checked_mul(std::mem::size_of::<f64>())
        .filter(|&bytes| isize::try_from(bytes).is_ok())
        .ok_or_else(|| ChtError::InvalidInput {
            field: "tensor.transient.storage",
            reason: "storage allocation exceeds the addressable entry count".into(),
        })?;
    let mut coefficient_w_k = Vec::new();
    coefficient_w_k
        .try_reserve_exact(cells)
        .map_err(|_| ChtError::InvalidInput {
            field: "tensor.transient.storage",
            reason: "storage allocation refused".into(),
        })?;
    for cell in 0..cells {
        if cell.is_multiple_of(4096) {
            poll(gate)?;
        }
        finite_positive(
            "tensor.transient.heat_capacity_j_k",
            heat_capacity_j_k[cell],
        )?;
        finite_positive(
            "tensor.transient.previous_temperature_k",
            previous_temperature_k[cell],
        )?;
        let coefficient = heat_capacity_j_k[cell] / time_step_s;
        finite_positive("tensor.transient.storage_w_k", coefficient)?;
        coefficient_w_k.push(coefficient);
    }
    Ok(PreparedTensorStep {
        coefficient_w_k,
        heat_capacity_j_k,
        previous_temperature_k,
        time_s,
        time_step_s,
    })
}
