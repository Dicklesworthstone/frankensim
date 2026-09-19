//! Static contact-set iteration over the shared bounded scalar reaction solver.
//! Dynamics retain their existing joint normal/friction sweep; this stationary
//! solve has no time-dependent loss or tangential sliding state.
use super::*;
use super::super::contact::{ModalContactConfig, force_tolerance, solve_reaction};
use super::super::contact::multiple::MultiContactConfig;

pub(super) fn admit_setup(
    network: &CoupledModalSystem, p: usize, config: MultiContactConfig,
) -> Result<(), ModalCouplingError> {
    if network.samples_rendered() != 0 || !(1..=32).contains(&config.max_contacts)
        || p == 0 || p > config.max_contacts || !(1..=128).contains(&config.max_sweeps) {
        return Err(invalid("multiple contacts require a sample-zero network and bounded nonempty contact/sweep counts"));
    }
    let n = network.mode_count();
    let k = network.columns.len();
    let terms = n.checked_mul(k+2).and_then(|v| v.checked_add((k+1)*(k+1)))
        .and_then(|v| v.checked_mul(p))
        .and_then(|v| n.checked_mul(p*p).and_then(|w| v.checked_add(w)))
        .ok_or_else(|| invalid("multi-contact setup work overflow"))?;
    if terms > config.max_setup_terms {
        return Err(invalid("multi-contact setup exceeds max_setup_terms"));
    }
    Ok(())
}

// The caller supplies the original static potential derivative. No contact
// force, gap or stiffness is changed to accelerate coordinate convergence.
pub(super) fn solve_joint(
    compliance: &[f64], free_x: &[f64], config_at: impl Fn(usize) -> ModalContactConfig,
    law_at: impl Fn(usize, f64) -> Result<f64, ModalCouplingError>,
    max_sweeps: usize, reactions: &mut [f64], gate: Option<&CancelGate>,
) -> Result<usize, ModalCouplingError> {
    let p = reactions.len();
    reactions.fill(0.0);
    let mut worst = (0.0_f64, 1.0_f64);
    for sweep in 1..=max_sweeps {
        for i in 0..p {
            poll(gate)?;
            let mut free = free_x[i];
            for j in 0..p {
                if i != j { free = finite(free - compliance[i*p+j]*reactions[j])?; }
            }
            let mut local = config_at(i);
            if local.force_absolute_tolerance_n*0.125 > 0.0 { local.force_absolute_tolerance_n *= 0.125; }
            if local.force_relative_tolerance*0.125 > 0.0 { local.force_relative_tolerance *= 0.125; }
            let evaluate = |reaction: f64| -> Result<(f64, f64), ModalCouplingError> {
                poll(gate)?;
                let x = finite(free - compliance[i*p+i]*reaction)?;
                let expected = law_at(i, x)?;
                Ok((finite(reaction - expected)?, force_tolerance(reaction, expected, local)?))
            };
            reactions[i] = match solve_reaction(evaluate, local) {
                Ok((r, _)) => r,
                // Another contact may relieve this bounded scratch iterate.
                // It is NOT accepted unless the final joint equations pass.
                Err(ModalCouplingError::Budget { what: "contact normal force", .. }) => local.maximum_force_n,
                Err(error) => return Err(error),
            };
        }
        worst = (0.0, 1.0);
        for i in 0..p {
            poll(gate)?;
            let mut x = free_x[i];
            for j in 0..p { x = finite(x - compliance[i*p+j]*reactions[j])?; }
            let expected = law_at(i, x)?;
            let residual = finite(reactions[i] - expected)?;
            let tolerance = force_tolerance(reactions[i], expected, config_at(i))?;
            if residual.abs()/tolerance > worst.0.abs()/worst.1 { worst = (residual, tolerance); }
        }
        if worst.0.abs() <= worst.1 { return Ok(sweep); }
    }
    Err(ModalCouplingError::ContactSolve {
        residual_n: worst.0, tolerance_n: worst.1, iterations: max_sweeps,
    })
}
