//! Implicit node-temperature derivatives of the same closed circulation.
//! Local branch source derivatives use the original segment kernels; the
//! frozen node-mixing matrix and its transpose close the circulating feedback.

use super::*;

pub(super) fn temperatures(linear: &TransportLinearization<'_, '_>, cx: &Cx<'_>,
    direction: &TransportDirection) -> Result<Vec<f64>, TransportError> {
    let mut rhs = vec![0.0; linear.mixers.len()];
    for (node, mixer) in linear.mixers.iter().enumerate() {
        poll(cx)?;
        for &(branch, fraction) in &mixer.incoming {
            let mut outlet = 0.0;
            for row in linear.network.row_range(branch) {
                poll(cx)?;
                outlet = dot(linear.segments[row].outlet,
                    [outlet, direction.walls_k[row], direction.log_conductances[row]])?;
            }
            add(&mut rhs[node], fraction * outlet)?;
        }
    }
    linear.network.circulation.as_ref().expect("closed circulation").solve(cx, &rhs, false)
}

pub(super) fn pullback(linear: &TransportLinearization<'_, '_>, cx: &Cx<'_>,
    objective: &TransportObjective, mut nodes: Vec<f64>, mut gradient: TransportGradient,
    heat_weight: f64) -> Result<TransportGradient, TransportError> {
    // Differentiate the explicit output map, temporarily holding every node
    // temperature independent. Its node weights form the implicit adjoint RHS.
    for (branch, &weight) in objective.branch_outlets.iter().enumerate() {
        poll(cx)?;
        if linear.primal.branches[branch].outlet_temperature_k.is_none() { continue; }
        let inlet = reverse_branch(linear, cx, branch, weight, Some(&objective.references),
            heat_weight, &mut gradient)?;
        add(&mut nodes[linear.upstream[branch]], inlet)?;
    }
    let adjoint = linear.network.circulation.as_ref().expect("closed circulation")
        .solve(cx, &nodes, true)?;
    // Each node residual receives the wall/hA forcing of its incoming branches.
    // The implicit inlet dependence is already in the transposed matrix; do not
    // propagate it again as a DAG sweep or double-count a physical circulation.
    for (node, mixer) in linear.mixers.iter().enumerate() {
        poll(cx)?;
        for &(branch, fraction) in &mixer.incoming {
            reverse_branch(linear, cx, branch, fraction * adjoint[node], None, 0.0,
                &mut gradient)?;
        }
    }
    poll(cx)?;
    Ok(gradient)
}

#[allow(clippy::too_many_arguments)]
fn reverse_branch(linear: &TransportLinearization<'_, '_>, cx: &Cx<'_>, branch: usize,
    mut weight: f64, references: Option<&[f64]>, heat_weight: f64,
    gradient: &mut TransportGradient) -> Result<f64, TransportError> {
    for row in linear.network.row_range(branch).rev() {
        poll(cx)?;
        let kernel = linear.segments[row];
        let weights = [weight, references.map_or(0.0, |values| values[row]), heat_weight];
        add(&mut gradient.walls[row],
            dot([kernel.outlet[1], kernel.reference[1], kernel.heat[1]], weights)?)?;
        add(&mut gradient.log_conductances[row],
            dot([kernel.outlet[2], kernel.reference[2], kernel.heat[2]], weights)?)?;
        weight = dot([kernel.outlet[0], kernel.reference[0], kernel.heat[0]], weights)?;
    }
    Ok(weight)
}
