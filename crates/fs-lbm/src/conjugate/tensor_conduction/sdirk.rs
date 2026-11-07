//! Two-stage, second-order physical time integration on the native tensor operator.

use fs_exec::CancelGate;

use super::{
    ConductivityTensor, TensorConductionSolution, advance_tensor_conduction,
    solve_tensor_conduction_stage,
};
use crate::conjugate::energy::{EnergyConfig, PseudoStep, ThermalSetup};
use crate::conjugate::transient::{TransientRecord, record_conduction_step};
use crate::conjugate::{ChtError, SolidMaterial, VoxelDomain, finite, finite_positive, poll};

/// Diagonal coefficient of the two-stage, L-stable SDIRK2 tableau.
pub const TENSOR_SDIRK2_GAMMA: f64 = 1.0 - std::f64::consts::FRAC_1_SQRT_2;

/// An accepted SDIRK2 interval with both actual spatial solves retained.
#[derive(Debug, Clone, PartialEq)]
pub struct TensorConductionSdirk2Step {
    /// Physical sample times: start + gamma*dt, then start + dt.
    pub stage_time_s: [f64; 2],
    /// First stage and stiffly accurate endpoint, including their physical
    /// face temperatures, conductive heat and native sparse-solve evidence.
    pub stages: [TensorConductionSolution; 2],
    /// Whole-interval storage/source/reservoir heat in joules, using the
    /// tableau weights. Iterations sum all Krylov work across both stages.
    pub record: TransientRecord,
}

/// Advance frozen full-tensor conduction with physical two-stage SDIRK2.
///
/// The tableau has rows `[gamma, 0]`, `[1-gamma, gamma]`, weights
/// `[1-gamma, gamma]`, and `gamma = 1 - 1/sqrt(2)`. Both stages use
/// the same condensed native tensor face operator. For smooth forcing and
/// fixed constraints, the method is second order in time and L-stable.
/// This does not imply monotonicity for arbitrarily large steps or a discrete
/// maximum principle for rotated anisotropy.
///
/// `setups[0]` and `setups[1]` are independently supplied physical source,
/// boundary, sink and contact samples at the returned stage times. Evaluate
/// time-dependent forcing at those times; repeating one endpoint sample
/// generally loses second-order accuracy. Resolve forcing discontinuities
/// with interval boundaries rather than assuming order across a jump.
///
/// `heat_capacity_j_k` is complete sensible capacity in J/K per cell and
/// `tensors` is one frozen SPD conductivity tensor per cell. Neither is
/// multiplied by cell volume again or changed between stages. This entry
/// does not infer phase change, latent heat, or temperature-dependent laws.
/// Faces and zero-thickness contacts have no artificial heat capacity.
///
/// The immutable `previous_temperature_k` is the physical start state in
/// BOTH stage equations. Stage 2 carries stage 1's actual source minus
/// conductive face/sink outflow as an internal heat rate, never as external
/// power or a fabricated previous temperature. Prescribed cells carry their
/// stage-1 reservoir storage rate. The final ledger integrates each stage's
/// actual external heat, including the stage-2 reservoir correction for
/// storage already carried from stage 1. This also admits changing prescribed
/// masks without counting that carried heat twice.
///
/// `config.max_iterations` is one shared Krylov budget for the entire
/// interval, including native refinements in both stages. Cancellation or
/// refusal of either stage returns no partial interval and changes no input.
///
/// # Errors
/// The frozen tensor admission rules apply at each stage. Invalid fields,
/// nonpositive capacities or temperatures, unrepresentable ordered stage
/// times/storage/heat, allocation refusal, an exhausted shared iteration
/// budget, an unconverged physical stage, and cancellation return an error.
#[allow(clippy::too_many_arguments, clippy::too_many_lines)]
pub fn advance_tensor_conduction_sdirk2(
    domain: &VoxelDomain,
    solids: &[SolidMaterial],
    setups: [&ThermalSetup; 2],
    tensors: &[ConductivityTensor],
    heat_capacity_j_k: &[f64],
    previous_temperature_k: &[f64],
    start_time_s: f64,
    time_step_s: f64,
    config: &EnergyConfig,
    gate: &CancelGate,
) -> Result<TensorConductionSdirk2Step, ChtError> {
    poll(gate)?;
    let (stage_interval_s, stage_time_s) =
        stage_clock(start_time_s, time_step_s).ok_or_else(|| invalid(
            "a finite nonnegative start and representable strictly ordered stage times are required",
        ))?;
    let first = advance_tensor_conduction(
        domain,
        solids,
        setups[0],
        tensors,
        heat_capacity_j_k,
        previous_temperature_k,
        start_time_s,
        stage_interval_s,
        config,
        gate,
    )?;
    let cells = domain.cell_count();
    // The first stage has admitted all field sizes, capacities and old
    // temperatures. Retain just two additional cell fields for stage 2.
    let mut coefficient_w_k = storage(cells)?;
    let mut carried_power_w = storage(cells)?;
    for cell in 0..cells {
        if cell.is_multiple_of(4096) {
            poll(gate)?;
        }
        let coefficient = heat_capacity_j_k[cell] / stage_interval_s;
        finite_positive("tensor.sdirk2.storage_w_k", coefficient)?;
        coefficient_w_k.push(coefficient);
        let power = setups[0].power_w.get(cell).copied().unwrap_or(0.0);
        let rate = first.solution.outward_heat_w[cell]
            .iter()
            .fold(power, |inward, outward| inward - outward);
        carried_power_w.push(rate);
    }
    for (index, sink) in setups[0].cell_sinks.iter().enumerate() {
        if index.is_multiple_of(4096) {
            poll(gate)?;
        }
        let temperature_difference = first.solution.temperature[sink.cell] - sink.temperature_k;
        carried_power_w[sink.cell] =
            (-sink.conductance_w_k).mul_add(temperature_difference, carried_power_w[sink.cell]);
    }
    // A fixed cell's total rate includes its thermostat. Native source
    // powers at fixed cells are not counted a second time.
    for (index, &(cell, _)) in setups[0].fixed_temperature.iter().enumerate() {
        if index.is_multiple_of(4096) {
            poll(gate)?;
        }
        carried_power_w[cell] = coefficient_w_k[cell]
            * (first.solution.temperature[cell] - previous_temperature_k[cell]);
    }
    let carry_weight = (1.0 - TENSOR_SDIRK2_GAMMA) / TENSOR_SDIRK2_GAMMA;
    for (cell, carried) in carried_power_w.iter_mut().enumerate() {
        if cell.is_multiple_of(4096) {
            poll(gate)?;
        }
        *carried *= carry_weight;
        finite("tensor.sdirk2.carried_power_w", *carried)?;
    }
    let remaining = config
        .max_iterations
        .checked_sub(first.solution.report.iterations)
        .ok_or_else(|| invalid("first stage exceeded the shared Krylov iteration budget"))?;
    let second_config = EnergyConfig {
        max_iterations: remaining,
        ..*config
    };
    let storage_step = PseudoStep {
        coefficient_w_k: &coefficient_w_k,
        previous: previous_temperature_k,
        junction_coefficient_w_k: &[],
    };
    let second = solve_tensor_conduction_stage(
        domain,
        solids,
        setups[1],
        tensors,
        &second_config,
        Some(&storage_step),
        Some(&carried_power_w),
        gate,
    )?;
    poll(gate)?;
    let iterations = first
        .solution
        .report
        .iterations
        .checked_add(second.report.iterations)
        .filter(|&count| count <= config.max_iterations)
        .ok_or_else(|| invalid("stages exceeded the shared Krylov iteration budget"))?;
    let solids_present = (0..cells).any(|cell| !domain.is_fluid(cell));
    // Reuse the physical ledger owner to account for actual stage-2 face
    // heat, sinks and constrained storage. This is an internal accounting
    // record, not a second backward-Euler interval.
    let mut record = record_conduction_step(
        domain,
        setups[1],
        &second.temperature,
        second.report.balance,
        iterations,
        |cell| heat_capacity_j_k[cell] * (second.temperature[cell] - previous_temperature_k[cell]),
        (stage_time_s[1], stage_interval_s, solids_present),
        0.0,
        0.0,
    )?;
    let mut fixed_carried_power_w = 0.0;
    // The second native solve has admitted unique, in-range fixed cells.
    for (index, &(cell, _)) in setups[1].fixed_temperature.iter().enumerate() {
        if index.is_multiple_of(4096) {
            poll(gate)?;
        }
        fixed_carried_power_w += carried_power_w[cell];
        finite("tensor.sdirk2.fixed_carried_power_w", fixed_carried_power_w)?;
    }
    (record.source_j, record.boundary_outflow_j) = integrated_heat(
        first.record.source_j,
        first.record.boundary_outflow_j,
        record.source_j,
        record.boundary_outflow_j,
        stage_interval_s,
        fixed_carried_power_w,
    );
    record.closure_j =
        record.stored_energy_change_j - (record.source_j - record.boundary_outflow_j);
    for value in [record.source_j, record.boundary_outflow_j, record.closure_j] {
        finite("tensor.sdirk2.energy_record", value)?;
    }
    poll(gate)?;
    Ok(TensorConductionSdirk2Step {
        stage_time_s,
        stages: [first.solution, second],
        record,
    })
}

fn invalid(reason: &str) -> ChtError {
    ChtError::InvalidInput {
        field: "tensor.sdirk2",
        reason: reason.into(),
    }
}

fn storage(cells: usize) -> Result<Vec<f64>, ChtError> {
    cells
        .checked_mul(std::mem::size_of::<f64>())
        .filter(|&bytes| isize::try_from(bytes).is_ok())
        .ok_or_else(|| invalid("cell storage exceeds the addressable entry count"))?;
    let mut field = Vec::new();
    field
        .try_reserve_exact(cells)
        .map_err(|_| invalid("cell storage allocation refused"))?;
    Ok(field)
}

// Pure temporal arithmetic shared by the production step and direct numerical
// tests. Spatial heat remains owned by the native face solve.
fn stage_clock(start: f64, dt: f64) -> Option<(f64, [f64; 2])> {
    let a = TENSOR_SDIRK2_GAMMA * dt;
    let times = [start + a, start + dt];
    (start.is_finite()
        && start >= 0.0
        && dt.is_finite()
        && dt > 0.0
        && a > 0.0
        && times[0].is_finite()
        && times[1].is_finite()
        && start < times[0]
        && times[0] < times[1])
        .then_some((a, times))
}

fn integrated_heat(
    first_source_j: f64,
    first_outflow_j: f64,
    second_source_j: f64,
    second_outflow_j: f64,
    stage_interval_s: f64,
    fixed_carried_power_w: f64,
) -> (f64, f64) {
    let weight = (1.0 - TENSOR_SDIRK2_GAMMA) / TENSOR_SDIRK2_GAMMA;
    let source_j = weight.mul_add(first_source_j, second_source_j);
    // The stage-2 thermostat supplies (C*dT - b*rate_1), not C*dT.
    // Its already-carried part therefore adds to net reservoir outflow.
    let boundary_outflow_j = stage_interval_s.mul_add(
        fixed_carried_power_w,
        weight.mul_add(first_outflow_j, second_outflow_j),
    );
    (source_j, boundary_outflow_j)
}

#[cfg(test)]
mod tests {
    use super::{TENSOR_SDIRK2_GAMMA, integrated_heat, stage_clock};

    fn close(actual: f64, expected: f64, tolerance: f64) {
        assert!(
            (actual - expected).abs() < tolerance,
            "{actual} != {expected}"
        );
    }

    #[test]
    fn clock_admits_only_resolved_physical_stage_times() {
        let (a, times) = stage_clock(2.0, 4.0).unwrap();
        close(a, 4.0 * TENSOR_SDIRK2_GAMMA, 1e-15);
        assert_eq!(times, [2.0 + a, 6.0]);
        for (start, dt) in [
            (-1.0, 1.0),
            (0.0, 0.0),
            (0.0, -1.0),
            (f64::NAN, 1.0),
            (0.0, f64::INFINITY),
            (1e16, 2.0),
            (f64::MAX, f64::MAX),
        ] {
            assert!(stage_clock(start, dt).is_none(), "{start} {dt}");
        }
    }

    #[test]
    fn stage_heat_quadrature_exactly_integrates_linear_forcing() {
        let start = 3.0;
        let dt = 2.0;
        let (a, time) = stage_clock(start, dt).unwrap();
        let power = |t| 5.0 + 7.0 * t;
        let (source, outflow) =
            integrated_heat(a * power(time[0]), 0.0, a * power(time[1]), 0.0, a, 0.0);
        close(
            source,
            5.0 * dt + 3.5 * ((start + dt).powi(2) - start.powi(2)),
            1e-12,
        );
        assert_eq!(outflow, 0.0);
    }

    #[test]
    fn original_state_tableau_has_second_order_diffusive_decay() {
        let gamma = TENSOR_SDIRK2_GAMMA;
        let march = |steps: usize| {
            let dt = 1.0 / steps as f64;
            let (a, _) = stage_clock(0.0, dt).unwrap();
            let mut old = 1.0;
            for _ in 0..steps {
                let first = old / (1.0 + a);
                let carried = (1.0 - gamma) / gamma * -first;
                // Actual stage-2 cell equation: (T2-T0)/a + T2 = carry.
                let second = (old + a * carried) / (1.0 + a);
                let z = dt;
                let amplification = (1.0 + (2.0 * gamma - 1.0) * z) / (1.0 + gamma * z).powi(2);
                close(second, amplification * old, 1e-14);
                let (_, outflow) = integrated_heat(0.0, a * first, 0.0, a * second, a, 0.0);
                close(second - old, -outflow, 1e-14);
                old = second;
            }
            old
        };
        let exact = (-1.0_f64).exp();
        let coarse = (march(10) - exact).abs();
        let fine = (march(20) - exact).abs();
        assert!((3.9..4.1).contains(&(coarse / fine)), "{coarse} {fine}");
    }

    #[test]
    fn changing_prescribed_mask_retains_original_storage_and_reservoir_heat() {
        let (a, _) = stage_clock(0.0, 2.0).unwrap();
        let r = (1.0 - TENSOR_SDIRK2_GAMMA) / TENSOR_SDIRK2_GAMMA;
        let capacity = [2.0, 3.0];
        let old = [300.0, 310.0];
        let g = 2.0 / 3.0;
        // First cell prescribed in stage 1, second cell free.
        let first = [
            350.0,
            (capacity[1] * old[1] + a * g * 350.0) / (capacity[1] + a * g),
        ];
        let q1 = g * (first[0] - first[1]);
        let first_rate = [capacity[0] * (first[0] - old[0]) / a, q1];
        // Reverse the mask in stage 2. T0 remains the reference in both rows.
        let second = [
            (capacity[0] * old[0] + a * (r * first_rate[0] + g * 330.0)) / (capacity[0] + a * g),
            330.0,
        ];
        let q2 = g * (second[0] - second[1]);
        let first_outflow = -a * q1 - capacity[0] * (first[0] - old[0]);
        let second_outflow = a * q2 - capacity[1] * (second[1] - old[1]);
        let (source, outflow) = integrated_heat(
            0.0,
            first_outflow,
            0.0,
            second_outflow,
            a,
            r * first_rate[1],
        );
        let stored = capacity[0] * (second[0] - old[0]) + capacity[1] * (second[1] - old[1]);
        assert_eq!(source, 0.0);
        close(stored, -outflow, 1e-12);
        let omitted_correction = r * first_outflow + second_outflow;
        assert!((stored + omitted_correction).abs() > 20.0);
    }

    #[test]
    fn repeated_fixed_cell_storage_is_counted_once() {
        let (a, _) = stage_clock(1.0, 3.0).unwrap();
        let r = (1.0 - TENSOR_SDIRK2_GAMMA) / TENSOR_SDIRK2_GAMMA;
        let c = 7.0;
        let first_storage = c * (340.0 - 300.0);
        let final_storage = c * (360.0 - 300.0);
        let (_, outflow) = integrated_heat(
            0.0,
            -first_storage,
            0.0,
            -final_storage,
            a,
            r * first_storage / a,
        );
        close(outflow, -final_storage, 1e-12);
    }
}
