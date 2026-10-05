//! Physical air-capacity and convection controls of a complete thermal adjoint.
//! The solid/radiation/wall feedback has ALREADY been solved by the caller.

use super::{AirPath, BTreeSet, Cx, Result, admitted_exchange_terms, bad, finite, poll, zeros};

/// Physical controls, not derivatives of a frozen effective air reference.
#[derive(Debug, Clone, PartialEq)]
pub struct AirTransportGradient {
    /// dJ/dln(mass flow * specific heat), in independent path order.
    /// Geometry and every convective coefficient are held fixed.
    pub log_capacity_rates: Vec<f64>,
    /// dJ/dln(h), in region order, including h's effect on downstream air.
    /// The direct consistent Robin contribution supplied by the caller is
    /// included exactly once. Wetted areas and capacity rates are fixed.
    pub log_htc: Vec<f64>,
}

/// Compose a complete solid/air dual with the physical exponential air law.
///
/// `reference_bars` differentiates independent additive changes to the ORIGINAL
/// convective Robin references. `direct_log_htc_bars` differentiates the same
/// convective coefficients at FIXED references. Both must come from the SAME
/// full-system nodal dual, including any smooth material or radiation feedback.
/// Neither frozen-solid bars nor already composed air-h derivatives may be
/// supplied. The caller owns this binding and checks its physical primal/dual.
///
/// Reverse each actual path march at fixed wall means. For z=h*A/C,
/// d(epsilon)/dln(z)=z*exp(-z) and d(g)/dln(z)=exp(-z)-g, g=epsilon/z.
/// Every segment's change reaches all downstream references on its own path.
/// Hence a fan-affinity control has dJ/dln(s)=sum(log_capacity_rates)
/// + sum(log_htc * dln(h)/dln(s)), NOT just the direct Robin h contraction.
/// The hydraulic caller must independently establish its flow scaling law.
///
/// This is an Estimated local derivative of the admitted exchange model, not
/// a hydraulic solve, geometry/flow-direction derivative, inverse certificate,
/// uncertainty interval or derivative of a moving temperature maximum. O(ports)
/// work/storage; no new solid solve or dense Jacobian. `max_ports` explicitly
/// bounds all regions before allocation. No partial result on refusal.
///
/// # Errors
/// Wrong lengths/order, duplicate trace ownership, invalid wall temperatures,
/// exceeded port allowance, nonfinite arithmetic, or cancellation.
#[allow(clippy::too_many_arguments)]
pub fn pullback_transport_controls(
    cx: &Cx<'_>, paths: &[AirPath], regions: &[&str], walls_k: &[f64],
    reference_bars: &[f64], direct_log_htc_bars: &[f64], max_ports: usize,
) -> Result<AirTransportGradient> {
    poll(cx)?;
    if paths.is_empty() || paths.len() > max_ports {
        return Err(bad("air transport controls require nonempty paths within the port budget"));
    }
    let mut n = 0_usize;
    for path in paths {
        poll(cx)?;
        n = n.checked_add(path.segments().len()).ok_or_else(|| bad("air transport port count overflow"))?;
        if n > max_ports { return Err(bad("air transport control port budget exhausted")); }
    }
    if [regions.len(), walls_k.len(), reference_bars.len(), direct_log_htc_bars.len()]
        .iter().any(|&len| len != n) {
        return Err(bad("one named wall and both direct dual controls per air segment are required"));
    }
    let mut seen = BTreeSet::new();
    for (i, segment) in paths.iter().flat_map(|p| p.segments()).enumerate() {
        poll(cx)?;
        if segment.region() != regions[i] || !seen.insert(regions[i]) {
            return Err(bad("air transport controls require exact ordering and unique trace ownership"));
        }
        finite(walls_k[i])?;
        finite(reference_bars[i])?;
        finite(direct_log_htc_bars[i])?;
    }
    let mut capacity = zeros(cx, paths.len())?;
    let mut htc = zeros(cx, n)?;
    htc.copy_from_slice(direct_log_htc_bars);
    let mut start = 0;
    for (branch, path) in paths.iter().enumerate() {
        poll(cx)?;
        let end = start + path.segments().len();
        let march = path.march(&walls_k[start..end])?;
        let mut outlet_bar = 0.0;
        for (j, segment) in path.segments().iter().enumerate().rev() {
            poll(cx)?;
            let i = start+j;
            let (z, _, g) = admitted_exchange_terms(segment, path.capacity_rate_w_per_k())?;
            let carry = fs_math::det::exp(-z);
            let dg = log_reference_weight_slope(z, g, carry);
            let delta = finite(walls_k[i] - march.segments[j].inlet_temperature_k)?;
            let outgoing = finite(outlet_bar * finite(z*carry)?)?;
            let local = finite(reference_bars[i]*dg)?;
            let ntu_bar = finite(delta * finite(outgoing-local)?)?;
            htc[i] = finite(htc[i]+ntu_bar)?;
            capacity[branch] = finite(capacity[branch]-ntu_bar)?;
            outlet_bar = finite(g.mul_add(reference_bars[i], finite(carry*outlet_bar)?))?;
        }
        start = end;
    }
    poll(cx)?;
    Ok(AirTransportGradient { log_capacity_rates: capacity, log_htc: htc })
}

// Avoid cancellation in exp(-z)-expm1(-z)/(-z) for weak exchange.
// Through z^6; the first omitted term has magnitude <= z^7/5760 here.
fn log_reference_weight_slope(z: f64, g: f64, carry: f64) -> f64 {
    if z < 1e-3 {
        z * (-0.5 + z*(1.0/3.0 + z*(-1.0/8.0
            + z*(1.0/30.0 + z*(-1.0/144.0 + z/840.0)))))
    } else { carry-g }
}

#[cfg(test)]
mod tests;
