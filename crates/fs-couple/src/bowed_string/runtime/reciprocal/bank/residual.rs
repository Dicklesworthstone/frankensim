//! Residual flexibility: the static compliance of the bending modes a string
//! does NOT retain, at its point supports (pexfp.299, MagentaBadger 7592).
//!
//! A point support seen through N full-span sines cannot form its kink, so the
//! retained basis alone is a Rayleigh-Ritz upper bound whose pitch error falls
//! as 1/N. The omitted modes still yield, almost statically, under a support
//! force. Their compliance is the exact static Green's function of the same
//! tension and EI pinned span minus the retained modes' static sum. All the
//! supports on one string share it as one matrix, cross terms included:
//! - A spring bridge condenses exactly to k/(1 + k G_bb). Its dashpot condenses
//!   to first order in c G omega, c/(1 + k G_bb)^2.
//! - A stop's closure gains the bridge cross term, minus k_eff G_sb times the
//!   bridge link's own column, and the stops' compliance among themselves
//!   becomes G_ss' - k_eff G_sb G_s'b. The obstacle owner books that storage
//!   and the declared numerical loss of its step-held quasi-static response.
//! - `MultiBowStringBank::with_bow_residual_flexibility` uses this same Schur
//!   complement for fixed bow maps and their joint omitted-mode compliance.

use super::*;
use fs_math::det;

/// Exact static deflection [m/N] at span fraction `x` under a unit transverse
/// force at `y`: EI w'''' - T w'' = delta(y), pinned at both ends of the span.
/// EI = 0 is the taut string, x(L - y)/(T L) for x <= y.
#[must_use]
pub fn static_green_m_n(card: &BowedStringCard, x: f64, y: f64) -> f64 {
    let (l, t) = (card.length_m, card.tension_n);
    let (a, b) = (x.min(y)*l, x.max(y)*l);
    let string = a*(l-b)/l;
    if !(card.bending_stiffness_n_m2 > 0.0) { return string/t; }
    let k = det::sqrt(t/card.bending_stiffness_n_m2);
    // sinh(k a) sinh(k (L - b)) / (k sinh(k L)), with every exponential at most one.
    let bending = det::exp(k*(a-b))*det::expm1(-2.0*k*a)*det::expm1(-2.0*k*(l-b))
        / (-2.0*det::expm1(-2.0*k*l)) / k;
    (string-bending)/t
}

/// The retained modes' share of `static_green_m_n`: over the card's
/// `mode_count` sines, the sum of 2 sin(n pi x) sin(n pi y) / (L (T k^2 + EI k^4)),
/// k = n pi / L. That is phi(x) phi(y) / omega^2 for each mass-normalized mode.
#[must_use]
pub fn retained_compliance_m_n(card: &BowedStringCard, x: f64, y: f64) -> f64 {
    let (lo, hi) = (x.min(y), x.max(y));
    let mut sum = 0.0;
    for n in 1..=card.mode_count {
        let order = n as f64*core::f64::consts::PI;
        let k = order/card.length_m;
        let stiffness = card.tension_n*k*k + card.bending_stiffness_n_m2*k*k*k*k;
        sum += 2.0*det::sin(order*lo)*det::sin(order*hi)/(card.length_m*stiffness);
    }
    sum
}

/// What the omitted modes add at (x, y) [m/N]: `static_green_m_n` minus
/// `retained_compliance_m_n`. Symmetric in x and y.
///
/// # Errors
/// Stations outside the open span, or a span without positive tension.
pub fn residual_compliance_m_n(card: &BowedStringCard, x: f64, y: f64) -> Result<f64, ReciprocalStringError> {
    if !interior(x) || !interior(y) || !(card.tension_n > 0.0) || !(card.length_m > 0.0)
        || !(card.bending_stiffness_n_m2 >= 0.0) {
        return Err(invalid("residual compliance needs interior stations on a tensioned span"));
    }
    let g = static_green_m_n(card, x, y) - retained_compliance_m_n(card, x, y);
    if !g.is_finite() { return Err(invalid("residual compliance is not representable")); }
    Ok(g)
}

/// A string/body connection as the retained modes see it through the omitted
/// modes' compliance at its station.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CondensedBridge {
    /// Omitted-mode compliance at the station, along the connection [m/N].
    pub residual_m_n: f64,
    /// k / (1 + k G) [N/m]: exact for the spring.
    pub stiffness_n_m: f64,
    /// c / (1 + k G)^2 [N s/m]: first order in c G omega.
    pub damping_n_s_m: f64,
}

/// A spring/dashpot of `stiffness_n_m` and `damping_n_s_m` in series with
/// `residual_m_n` of omitted-mode compliance along its own axis.
///
/// # Errors
/// A negative or unrepresentable residual.
pub fn condense(residual_m_n: f64, stiffness_n_m: f64, damping_n_s_m: f64)
    -> Result<CondensedBridge, ReciprocalStringError>
{
    let factor = 1.0 + stiffness_n_m*residual_m_n;
    if !(residual_m_n >= 0.0) || !factor.is_finite() {
        return Err(invalid("bridge residual compliance must be finite and nonnegative"));
    }
    Ok(CondensedBridge { residual_m_n, stiffness_n_m: stiffness_n_m/factor,
        damping_n_s_m: damping_n_s_m/(factor*factor) })
}

/// Condense one connection of `source` (its card `card`) for the residual.
///
/// # Errors
/// Those of `residual_compliance_m_n` and `condense`.
pub fn condense_bridge(source: &MaterialStringSource, card: &BowedStringCard,
    connection: &InteriorStringConnection) -> Result<CondensedBridge, ReciprocalStringError>
{
    let (_, axis) = directions(source);
    let station = connection.station_fraction;
    let g = dot2(axis, axis)*residual_compliance_m_n(card, station, station)?;
    condense(g, connection.stiffness_n_m, connection.damping_n_s_m)
}

// Bow tangent and bridge axis in the string's transverse frame. A one-plane
// source has only the first plane, as its projector does.
pub(super) fn directions(source: &MaterialStringSource) -> ([f64; 2], [f64; 2]) {
    match source {
        MaterialStringSource::WithPolarizations(p) => (p.directions().bow(), p.directions().bridge()),
        _ => ([1.0, 0.0], [1.0, 0.0]),
    }
}
fn dot2(a: [f64; 2], b: [f64; 2]) -> f64 { a[0]*b[0] + a[1]*b[1] }

/// One stop for the residual: its string, station and signed unit closure
/// direction in that string's transverse frame.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ResidualStop {
    /// Index in the original material string list.
    pub string: usize,
    /// Interior span fraction.
    pub station_fraction: f64,
    /// Closure direction, the stop side's sign included.
    pub direction: [f64; 2],
}

/// One string's single bridge connection as that string's stops see it.
#[derive(Clone, Copy, Debug)]
pub struct ResidualBridge<'a> {
    /// The string's own card: tension, EI, span and retained modes.
    pub card: &'a BowedStringCard,
    /// The bridge's interior span fraction.
    pub station_fraction: f64,
    /// The bridge axis in the string's transverse frame.
    pub axis: [f64; 2],
    /// The condensed connection (`condense`).
    pub condensed: CondensedBridge,
}

/// Every string's one bridge connection, condensed (`condense_bridge`).
///
/// # Errors
/// A string with other than one connection, mismatched cards, or a bad residual.
pub fn residual_bridges<'a>(strings: &[MaterialStringMember], cards: &'a [BowedStringCard])
    -> Result<Vec<ResidualBridge<'a>>, ReciprocalStringError>
{
    if strings.len() != cards.len() { return Err(invalid("residual bridges need one card per string")); }
    strings.iter().zip(cards).map(|(member, card)| {
        let [connection] = member.connections.as_slice() else {
            return Err(invalid("residual flexibility condenses exactly one bridge connection per string"));
        };
        Ok(ResidualBridge { card, station_fraction: connection.station_fraction,
            axis: directions(&member.source).1, condensed: condense_bridge(&member.source, card, connection)? })
    }).collect()
}

/// The omitted modes' coupling of point forces with one another and with their own
/// strings' bridges (`bridges`, one per string), condensed through the bridges.
/// Fixed bows reuse this force-law-independent projection with their tangents.
///
/// Returns, per stop, the bridge cross factor k_eff G_sb: the stop's closure
/// column gains minus this times its bridge link's column (string and body
/// parts). Also returns the stops' row-major compliance G_ss' - k_eff G_sb G_s'b
/// [m/N], zero across different strings, computed once per pair and mirrored.
///
/// # Errors
/// A stop naming an unknown string, or a bad residual.
pub fn stop_residual(bridges: &[ResidualBridge<'_>], stops: &[ResidualStop])
    -> Result<(Vec<f64>, Vec<f64>), ReciprocalStringError>
{
    let q = stops.len();
    let mut toward = Vec::with_capacity(q);
    for stop in stops {
        let bridge = bridges.get(stop.string).ok_or_else(|| invalid("residual stop names an unknown string"))?;
        toward.push(dot2(stop.direction, bridge.axis)
            *residual_compliance_m_n(bridge.card, stop.station_fraction, bridge.station_fraction)?);
    }
    let cross = stops.iter().zip(&toward).map(|(stop, g)| bridges[stop.string].condensed.stiffness_n_m*g).collect();
    let mut compliance = vec![0.0; q*q];
    for i in 0..q {
        for j in i..q {
            if stops[i].string != stops[j].string { continue; }
            let bridge = &bridges[stops[i].string];
            let own = dot2(stops[i].direction, stops[j].direction)
                *residual_compliance_m_n(bridge.card, stops[i].station_fraction, stops[j].station_fraction)?;
            let g = own - bridge.condensed.stiffness_n_m*(toward[i]*toward[j]);
            if !g.is_finite() { return Err(invalid("stop residual compliance is not representable")); }
            compliance[i*q+j] = g;
            compliance[j*q+i] = g;
        }
    }
    Ok((cross, compliance))
}

#[cfg(test)]
mod tests;
