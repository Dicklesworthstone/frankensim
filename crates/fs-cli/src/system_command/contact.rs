//! `contact_interface`: a normal contact between two scene ends,
//! `fs-dcontact`'s [`ContactInterface`], on the kernel's load line.
//!
//! The scene authors one power law `R = K [x - gap]_+^exponent` on the
//! closure `x`, the displacement of end `a` towards end `b`, and the start
//! of `x`. The owner steps the reaction as the discrete gradient of the
//! law's potential over each step's move and settles it on the load line by
//! bracketed bisection, inside the budgets the scene gives. An optional
//! Hunt–Crossley loss uses the same owner's native law and exports the lost
//! mechanical energy through its heat port. Nothing here evaluates the law
//! or solves anything.

use fs_dcontact::Obstacle;
use fs_dcontact::component::{ContactInterface, ContactSolve, ContactState};
use fs_ports::component::{Component, PortSpec};
use fs_system::{NodeId, SystemBuilder};

use super::friction::reads_ends;
use super::{J, Refusal, fields, number, number_or};

/// The kernel holds a contact's closure and its last reaction.
pub(super) type State = ContactState;

// The defaults of the reaction search: the budgets of fs-dcontact's own
// component tests (tests/component.rs `solve`). The force and penetration
// ceilings have no default: they are the scene's scale.
const MAX_ITERATIONS: f64 = 100.0;
const FORCE_ABSOLUTE_TOLERANCE_N: f64 = 1.0e-9;
const FORCE_RELATIVE_TOLERANCE: f64 = 1.0e-12;

/// The reaction search's budgets: the iterations and tolerances with their
/// defaults, the two ceilings as given.
fn solve(definition: &J, what: &str) -> Result<ContactSolve, String> {
    let inner = format!("{what}.solve");
    let solve = definition
        .get("solve")
        .ok_or_else(|| format!("{what}: missing solve (maximum_force_n, maximum_penetration_m)"))?;
    fields(
        solve,
        &[
            "max_iterations",
            "maximum_force_n",
            "maximum_penetration_m",
            "force_absolute_tolerance_n",
            "force_relative_tolerance",
        ],
        &inner,
    )?;
    let iterations = number_or(solve, "max_iterations", MAX_ITERATIONS, &inner)?;
    if !(1.0..=128.0).contains(&iterations) || iterations.fract() != 0.0 {
        return Err(format!(
            "{inner}: max_iterations is a whole number from 1 to 128"
        ));
    }
    // A whole number from 1 to 128, checked above.
    #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
    let max_iterations = iterations as usize;
    Ok(ContactSolve {
        max_iterations,
        maximum_force_n: number(solve, "maximum_force_n", &inner)?,
        maximum_penetration_m: number(solve, "maximum_penetration_m", &inner)?,
        force_absolute_tolerance_n: number_or(
            solve,
            "force_absolute_tolerance_n",
            FORCE_ABSOLUTE_TOLERANCE_N,
            &inner,
        )?,
        force_relative_tolerance: number_or(
            solve,
            "force_relative_tolerance",
            FORCE_RELATIVE_TOLERANCE,
            &inner,
        )?,
    })
}

/// A `contact_interface` element: the owner's interface on the kernel's load
/// line, starting at the closure the scene declares.
pub(super) fn interface(
    builder: &mut SystemBuilder,
    name: &str,
    definition: &J,
    what: &str,
) -> Result<(NodeId, Vec<PortSpec>), Refusal> {
    fields(
        definition,
        &[
            "kind",
            "stiffness",
            "exponent",
            "internal_loss_s_per_m",
            "gap_m",
            "closure_m",
            "solve",
            "reads_ends",
        ],
        what,
    )?;
    let ends = reads_ends(definition, what)?;
    let gap_m = number_or(definition, "gap_m", 0.0, what)?;
    let closure_m = number_or(definition, "closure_m", 0.0, what)?;
    let internal_loss = number_or(definition, "internal_loss_s_per_m", 0.0, what)?;
    if !closure_m.is_finite() {
        return Err(format!("{what}: closure_m must be finite").into());
    }
    // One point on the unit opening map, weight one: the closure itself.
    let law = Obstacle::new(
        vec![-1.0],
        1,
        1,
        vec![gap_m],
        vec![1.0],
        number(definition, "stiffness", what)?,
        number(definition, "exponent", what)?,
        format!("frankensim system scene element {name}: an authored power law"),
    )
    .and_then(|law| law.with_internal_loss(internal_loss))
    .map_err(|error| format!("{what}: {error}"))?;
    let mut component = ContactInterface::new(law, solve(definition, what)?)
        .map_err(|error| format!("{what}: {error}"))?;
    // The kernel wires each end input read here, `a.end` or `b.end`, to the
    // end its neighbours tell.
    if ends.contains(&true) {
        component = component
            .reading_ends(ends)
            .map_err(|error| format!("{what}: {error}"))?;
    }
    let ports = component.ports().to_vec();
    let node = builder
        .add_on_line(name, component, ContactState::at(closure_m))
        .map_err(|error| Refusal::Admission(error.to_string()))?;
    Ok((node, ports))
}

/// The closure at the end of the last step and the reaction it took.
pub(super) fn values(state: Option<&State>) -> Vec<(&'static str, f64)> {
    vec![
        ("closure_m", state.map_or(f64::NAN, |state| state.closure_m)),
        (
            "reaction_n",
            state.map_or(f64::NAN, |state| state.reaction_n),
        ),
    ]
}
