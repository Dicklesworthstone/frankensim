//! The same data, geometry and source laws through the library implicit owner.
use super::*;
use fs_lbm::conjugate::enthalpy::transport::interfaces::implicit::{
    ImplicitJouleConfig,step_joule_enthalpy_implicit,
};
use fs_solver::NewtonKrylovConfig;

/// Advance one output interval without emitting partially converged substeps.
/// The explicit and implicit routes share all geometry, tables and physical
/// laws. This presentation adapter owns no current/heat/phase residual.
pub(super) fn advance<'m>(
    old:&DepositionSimulation<'m>,target:f64,max_dt:f64,input:JouleEnthalpyInput<'_>,
    table:&Table,contact:Option<&ContactTable>,gate:&CancelGate,difference_step_j:f64,
)->Result<IntervalResult<'m>,Box<dyn Error>> {
    let duration=target-old.time_s();
    if !duration.is_finite() || duration<=0.0 || !max_dt.is_finite() || max_dt<=0.0
        || !difference_step_j.is_finite() || difference_step_j<=0.0 {
        return Err("invalid implicit duration, maximum step or numerical stencil".into());
    }
    // Equal subintervals avoid a floating-point tail of nearly zero duration.
    // This partition is never allowed to coarsen the caller's maximum step.
    let count=(duration/max_dt).ceil().max(1.0);
    if !count.is_finite() || count>20_000.0 {return Err("implicit output interval exceeds 20000 steps".into());}
    let count=count as usize;
    let tolerance=1e-8/count as f64;
    let controls=ImplicitJouleConfig {
        max_unknowns:3,max_newton_iterations:24,difference_step_j,
        cell_energy_tolerance_j:tolerance,total_energy_tolerance_j:tolerance,
        newton:NewtonKrylovConfig {absolute_tolerance:tolerance/10.0,relative_tolerance:1e-10,
            forcing_minimum:1e-4,linear_restart:3,..NewtonKrylovConfig::default()},
    };
    let mut current=old.clone();let mut last=None;
    let (mut solves,mut iterations,mut electrode_j,mut external_j,mut worst)=(0usize,0usize,0.0,0.0,0.0_f64);
    for index in 1..=count {
        let stop=if index==count {target} else {old.time_s()+duration*(index as f64/count as f64)};
        // All probes and final current checks share the WHOLE output interval's
        // original budget; it is not silently renewed at every nonlinear step.
        let remaining=input.max_current_solves.checked_sub(solves).ok_or("implicit current budget")?;
        if remaining==0 {return Err("implicit current budget exhausted".into());}
        let guess=initial_guess(&current,input,table,difference_step_j);
        let (next,report)=step_joule_enthalpy_implicit(&current,stop,
            JouleEnthalpyInput {max_current_solves:remaining,..input},&[0,1],gate,controls,&guess,
            |_,cell,_|table.transport(cell.specific_h_j_kg),thermal_link,
            |_,face,sides|contact.ok_or(JouleEnthalpyError::Invalid("unexpected electrical interface"))?.response(face,sides),
            |_|Ok(std::collections::BTreeMap::new()))?;
        electrode_j+=report.electrode_energy_j;external_j+=report.thermal.external_heat_j;
        solves+=report.current_solves;iterations+=report.nonlinear.iterations;
        worst=worst.max(report.maximum_cell_residual_j);
        current=next;last=Some(report);
    }
    let report=last.ok_or("implicit interval has no accepted step")?;
    // Independent interval storage check before the existing observation stream
    // receives this endpoint. Individual substeps already passed local gates.
    let residual=current.energy_j()-old.energy_j()-external_j;
    if ![electrode_j,external_j,residual].iter().all(|x|x.is_finite()) || residual.abs()>1e-8 {
        return Err("implicit observation interval energy does not close".into());
    }
    Ok(IntervalResult {
        next:current,electrode_j,external_j,solves,
        current_a:report.endpoint.electric.electrodes[0].current_a,
        conductivity:report.material.cells[1].response.conductivity_w_m_k(),
        contact_resistance:report.interfaces.contacts.first().map(|s|s.response.resistance_ohm_m2()),
        newton_iterations:iterations,maximum_cell_residual_j:worst,
        temporal_indicator_j:0.0,rejected_intervals:0,
    })
}

/// Numerical interior guess shared by fixed and adaptive paths. Never preheats
/// accepted matter; the residual always starts from the original enthalpy.
fn initial_guess(state:&DepositionSimulation<'_>,input:JouleEnthalpyInput<'_>,table:&Table,scale:f64)->Vec<f64> {
    let powered=input.setup.electrodes.iter().any(|e|match e.drive {
        ElectrodeDrive::Voltage(v)=>v!=0.0,ElectrodeDrive::Current(i)=>i!=0.0,
    });
    state.cells().values().map(|c| {
        if powered && c.specific_h_j_kg==table.electrical[0][0] {8.0*scale} else {0.0}
    }).collect()
}

/// One requested output interval, using the library controller and unchanged
/// physical material/contact tables. A refused interval emits no partial row.
#[allow(clippy::too_many_arguments)]
pub(super) fn adaptive_advance<'m>(
    old:&DepositionSimulation<'m>,target:f64,max_dt:f64,input:JouleEnthalpyInput<'_>,
    table:&Table,contact:Option<&ContactTable>,gate:&CancelGate,difference_step_j:f64,
    cell_error_budget_j:f64,minimum_step_s:f64,
)->Result<IntervalResult<'m>,Box<dyn Error>> {
    use fs_lbm::conjugate::enthalpy::transport::interfaces::implicit::adaptive::{
        AdaptiveJouleConfig,advance_joule_enthalpy_adaptive,
    };
    let config=AdaptiveJouleConfig {
        step:ImplicitJouleConfig {
            max_unknowns:3,max_newton_iterations:24,difference_step_j,
            cell_energy_tolerance_j:1e-8,total_energy_tolerance_j:1e-8,
            newton:NewtonKrylovConfig {absolute_tolerance:1e-9,relative_tolerance:1e-10,
                forcing_minimum:1e-4,linear_restart:3,..NewtonKrylovConfig::default()},
        },
        maximum_step_s:max_dt,minimum_step_s,max_step_attempts:20_000,cell_error_budget_j,
    };
    let (next,report)=advance_joule_enthalpy_adaptive(old,target,input,&[0,1],gate,config,
        |state,_|Ok(initial_guess(state,input,table,difference_step_j)),
        |_,cell,_|table.transport(cell.specific_h_j_kg),thermal_link,
        |_,face,sides|contact.ok_or(JouleEnthalpyError::Invalid("unexpected electrical interface"))?.response(face,sides),
        |_|Ok(std::collections::BTreeMap::new()))?;
    Ok(IntervalResult {
        next,electrode_j:report.electrode_energy_j,external_j:report.external_heat_j,
        solves:report.current_solves,current_a:report.last_step.endpoint.electric.electrodes[0].current_a,
        conductivity:report.last_step.material.cells[1].response.conductivity_w_m_k(),
        contact_resistance:report.last_step.interfaces.contacts.first().map(|s|s.response.resistance_ohm_m2()),
        newton_iterations:report.accepted_newton_iterations,maximum_cell_residual_j:report.maximum_cell_residual_j,
        temporal_indicator_j:report.estimated_cell_errors_j.iter().map(|(_,e)|*e).fold(0.0,f64::max),
        rejected_intervals:report.rejected_intervals,
    })
}

#[cfg(test)]
#[path="implicit/tests.rs"]
mod tests;
