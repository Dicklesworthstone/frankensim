//! Implicit tangents/adjoints of the imposed return-temperature solve.
//! Do not differentiate elimination iterations and do not freeze the intake.
use super::*;
use super::super::recirculation::solve_feedback;

pub(super) fn apply(lin: &TransportLinearization<'_, '_>, cx: &Cx<'_>,
    direction: &TransportDirection) -> Result<TransportDifferential, TransportError> {
    lin.check_direction(cx, direction)?;
    let feedback = lin.network.feedback.as_ref().expect("feedback dispatch");
    let mut mixed = direction.clone();
    for &source in &feedback.sources { mixed.inlets_k[source] = 0.0; }
    let base = lin.apply_open(cx, &mixed)?;
    let mut rhs: Vec<f64> = feedback.sources.iter().zip(&feedback.fractions)
        .map(|(&node, &fraction)| (1.0 - fraction) * direction.inlets_k[node]).collect();
    for link in &feedback.links {
        poll(cx)?;
        let i = feedback.sources.binary_search(&link.supply_node).expect("admitted source");
        let value = base.node_temperatures_k[link.return_node]
            .ok_or(TransportError::InvalidInput("return tangent absent"))?;
        rhs[i] = finite(rhs[i] + link.fraction * value, "return tangent right-hand side")?;
    }
    let solution = solve_feedback(cx, &feedback.matrix, &rhs, false)?;
    for (&node, &value) in feedback.sources.iter().zip(&solution) { mixed.inlets_k[node] = value; }
    lin.apply_open(cx, &mixed)
}

pub(super) fn pullback(lin: &TransportLinearization<'_, '_>, cx: &Cx<'_>,
    objective: &TransportObjective) -> Result<TransportGradient, TransportError> {
    let feedback = lin.network.feedback.as_ref().expect("feedback dispatch");
    let mut result = lin.pullback_open(cx, objective)?;
    let rhs: Vec<_> = feedback.sources.iter().map(|&node| result.inlets[node]).collect();
    let lambda = solve_feedback(cx, &feedback.matrix, &rhs, true)?;
    let mut returned = lin.zero_objective();
    for link in &feedback.links {
        poll(cx)?;
        let i = feedback.sources.binary_search(&link.supply_node).expect("admitted source");
        add(&mut returned.node_temperatures[link.return_node], link.fraction * lambda[i])?;
    }
    let extra = lin.pullback_open(cx, &returned)?;
    for (out, extra) in result.walls.iter_mut().zip(&extra.walls) { poll(cx)?; add(out, *extra)?; }
    for (out, extra) in result.log_conductances.iter_mut().zip(&extra.log_conductances) { poll(cx)?; add(out, *extra)?; }
    for (out, extra) in result.inlets.iter_mut().zip(&extra.inlets) { poll(cx)?; add(out, *extra)?; }
    // These controls are FRESH temperatures, not the eliminated mixed values.
    for (i, &node) in feedback.sources.iter().enumerate() {
        result.inlets[node] = finite((1.0 - feedback.fractions[i]) * lambda[i], "fresh inlet adjoint")?;
    }
    poll(cx)?;
    Ok(result)
}
