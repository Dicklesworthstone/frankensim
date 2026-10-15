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
        let mut guess=vec![0.0;current.cells().len()];
        // At cold source turn-on, mixed Krylov directions need an interior
        // initial guess for the bounded h chart. This is only a Newton guess:
        // initial storage and every residual still use the ORIGINAL cold state.
        // An exactly unpowered cold system retains its exact zero-residual guess.
        let powered=input.setup.electrodes.iter().any(|e|match e.drive {
            ElectrodeDrive::Voltage(v)=>v!=0.0,ElectrodeDrive::Current(i)=>i!=0.0,
        });
        if powered {
            for (i,cell) in current.cells().values().enumerate() {
                if cell.specific_h_j_kg==table.electrical[0][0] {
                    guess[i]=8.0*difference_step_j;
                }
            }
        }
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
    })
}

#[cfg(test)]
#[path="implicit/tests.rs"]
mod tests;
