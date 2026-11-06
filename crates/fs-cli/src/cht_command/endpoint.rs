//! Explicit work bounds and acceptance reports for native endpoint coupling.

use fs_exec::CancelGate;
use fs_lbm::conjugate::{ChtError, MaterialEndpointConfig, MaterialEndpointMarch};
use std::fmt::Write as _;

use super::{J, Result, bad, num, number, solver_failure};

pub(super) const HELP: &str = "\nImplicit convection coupling: transient.endpoint_coupling is an object with\nmaximum_passes_per_step, maximum_momentum_solves, maximum_energy_solves,\nmaximum_source_refreshes, temperature_tolerance_k,\ntransport_relative_tolerance, force_tolerance_n_m3 and temperature_relaxation.\nAll fields are required; solve and refresh budgets cover the whole march.\nIt requires flow \"unsteady\" and energy \"march\". Native constant-property\nbuoyancy supports laminar or explicitly selected transient SGS transport;\nevolving source materials retain their admitted SGS, source-enthalpy\nand solid-storage combinations. Each coupling trial reuses\nthe same accepted physical history. Results report actual work and endpoint\ndefects; implicit coupling does not raise the thermal time integrator's order.\n";

pub(super) fn parse(value: &J, steps: usize) -> Result<MaterialEndpointConfig> {
    const AT: &str = "transient.endpoint_coupling";
    if value.as_object().is_none() {
        return Err(bad(format!("{AT} must be an object")));
    }
    let count = |key: &str, maximum: usize| -> Result<usize> {
        let n = number(value, key, AT)?;
        if !(1.0..=maximum as f64).contains(&n) || n.fract() != 0.0 {
            return Err(bad(format!(
                "{AT}.{key} must be a whole number in 1..={maximum}"
            )));
        }
        Ok(n as usize)
    };
    let config = MaterialEndpointConfig {
        maximum_passes_per_step: count("maximum_passes_per_step", 1024)?,
        maximum_momentum_solves: count("maximum_momentum_solves", 1_000_000)?,
        maximum_energy_solves: count("maximum_energy_solves", 1_000_000)?,
        maximum_source_refreshes: count("maximum_source_refreshes", 2_000_001)?,
        temperature_tolerance_k: number(value, "temperature_tolerance_k", AT)?,
        transport_relative_tolerance: number(value, "transport_relative_tolerance", AT)?,
        force_tolerance_n_m3: number(value, "force_tolerance_n_m3", AT)?,
        temperature_relaxation: number(value, "temperature_relaxation", AT)?,
    };
    if [
        config.temperature_tolerance_k,
        config.transport_relative_tolerance,
        config.force_tolerance_n_m3,
    ]
    .iter()
    .any(|v| !v.is_finite() || *v < 0.0)
        || config.transport_relative_tolerance >= 1.0
        || !config.temperature_relaxation.is_finite()
        || config.temperature_relaxation <= 0.0
        || config.temperature_relaxation > 1.0
    {
        return Err(bad(format!(
            "{AT} requires finite nonnegative tolerances, relative tolerance below one and temperature_relaxation in (0,1]"
        )));
    }
    if config.maximum_momentum_solves < steps
        || config.maximum_energy_solves < steps
        || config.maximum_source_refreshes <= steps
    {
        return Err(bad(format!(
            "{AT} whole-march budgets must cover every step and its endpoint, plus the opening source refresh"
        )));
    }
    Ok(config)
}

/// Serialize observations of the accepted native march without source queries.
pub(super) fn take<M>(run: MaterialEndpointMarch<M>, gate: &CancelGate) -> Result<(M, String)> {
    let c = run.config;
    let mut maxima = [0.0_f64; 3];
    for (index, record) in run.records.iter().enumerate() {
        if index % 256 == 0 && gate.is_requested() {
            return Err(solver_failure(&ChtError::Cancelled));
        }
        for (i, value) in [
            record.temperature_defect_k,
            record.transport_relative_defect,
            record.force_defect_n_m3,
        ]
        .into_iter()
        .enumerate()
        {
            num(value)?;
            maxima[i] = maxima[i].max(value);
        }
    }
    let mut out = format!(
        concat!(
            "{{\"mode\":\"implicit-endpoint\",\"accepted_steps\":{},",
            "\"momentum_solves\":{},\"energy_solves\":{},\"source_refreshes\":{},",
            "\"maximum_temperature_defect_k\":{},\"maximum_transport_relative_defect\":{},",
            "\"maximum_force_defect_n_m3\":{},\"controls\":{{",
            "\"maximum_passes_per_step\":{},\"maximum_momentum_solves\":{},",
            "\"maximum_energy_solves\":{},\"maximum_source_refreshes\":{},",
            "\"temperature_tolerance_k\":{},\"transport_relative_tolerance\":{},",
            "\"force_tolerance_n_m3\":{},\"temperature_relaxation\":{}}},\"records\":["
        ),
        run.records.len(),
        run.momentum_solves,
        run.energy_solves,
        run.source_refreshes,
        num(maxima[0])?,
        num(maxima[1])?,
        num(maxima[2])?,
        c.maximum_passes_per_step,
        c.maximum_momentum_solves,
        c.maximum_energy_solves,
        c.maximum_source_refreshes,
        num(c.temperature_tolerance_k)?,
        num(c.transport_relative_tolerance)?,
        num(c.force_tolerance_n_m3)?,
        num(c.temperature_relaxation)?,
    );
    let every = run.records.len().div_ceil(200).max(1);
    let mut first = true;
    for (index, record) in run.records.iter().enumerate() {
        if (index + 1) % every != 0 && index + 1 != run.records.len() {
            continue;
        }
        if gate.is_requested() {
            return Err(solver_failure(&ChtError::Cancelled));
        }
        if !first {
            out.push(',');
        }
        first = false;
        let _ = write!(
            out,
            concat!(
                "{{\"time_s\":{},\"passes\":{},\"temperature_defect_k\":{},",
                "\"transport_relative_defect\":{},\"force_defect_n_m3\":{}}}"
            ),
            num(record.time_s)?,
            record.passes,
            num(record.temperature_defect_k)?,
            num(record.transport_relative_defect)?,
            num(record.force_defect_n_m3)?,
        );
    }
    out.push_str("]}");
    Ok((run.solution, out))
}
