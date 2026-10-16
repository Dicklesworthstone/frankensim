//! Independent speed derivatives of heterogeneous series/parallel fan banks.
//! This differentiates the nominal hydraulic equations, not tolerance corners.

use super::{compare_banks, effective_curve, interpolate, interpolate_inverse};
use crate::{FanArrangement, FanBank, LossResistance};
use fs_exec::Cx;
use fs_qty::VolumetricFlowRate;

/// Why no complete local hydraulic derivative can be returned.
#[derive(Debug, Clone, PartialEq)]
pub enum SpeedDerivativeError {
    /// Malformed controls, arithmetic, or work policy.
    Invalid(&'static str),
    /// The complete member-curve traversal exceeds the caller's allowance.
    PointBudget { required: usize, allowed: usize },
    /// A member is at a validity endpoint, knot with unequal slopes, or plateau
    /// whose inverse is not locally unique. Index is in the supplied bank order.
    NonSmooth { bank: usize },
    /// The supplied operating point does not close the nominal equations.
    Residual { relative: f64, tolerance: f64 },
    /// The work context was cancelled; no partial gradient is returned.
    Cancelled,
}
impl core::fmt::Display for SpeedDerivativeError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Invalid(why) => write!(f, "fan-speed derivative: {why}"),
            Self::PointBudget { required, allowed } => write!(f, "fan-speed derivative needs {required} curve points, limit {allowed}"),
            Self::NonSmooth { bank } => write!(f, "fan bank {bank} has no admitted two-sided speed derivative at this point"),
            Self::Residual { relative, tolerance } => write!(f, "fan-speed nominal residual {relative:e} exceeds {tolerance:e}"),
            Self::Cancelled => f.write_str("fan-speed derivative cancelled"),
        }
    }
}
impl core::error::Error for SpeedDerivativeError {}
type Result<T> = core::result::Result<T, SpeedDerivativeError>;
fn checked(x: f64) -> Result<f64> {
    if x.is_finite() { Ok(x) } else { Err(SpeedDerivativeError::Invalid("nonfinite hydraulic differential")) }
}
fn poll(cx: &Cx<'_>) -> Result<()> {
    cx.checkpoint().map_err(|_| SpeedDerivativeError::Cancelled)
}

/// Local derivatives of the complete nominal fan/loss intersection.
#[derive(Debug, Clone, PartialEq)]
pub struct FanSpeedFlowGradient {
    /// d ln(total flow) / d ln(member speed), in the supplied bank order.
    /// Every terminal flow of a fixed quadratic loss network has this same
    /// logarithmic derivative; member flows INSIDE a parallel fan do not.
    pub log_flow_per_log_speed: Vec<f64>,
    /// Recomputed relative nominal pressure residual (series), or summed
    /// member-flow residual (parallel). Not a derivative-error bound.
    pub relative_residual: f64,
}

/// Differentiate independently controlled fan banks at a retained operating flow.
///
/// `topology` describes the connection BETWEEN banks; each bank retains its own
/// series/parallel identical-fan arrangement, count, source curve and speed.
/// The passive network must have fixed equivalent resistance R and pressure
/// R*Q^2. For a member pressure slope b=dP/dQ, affinity supplies the fixed-flow
/// partial dP/dln(s)=2P-Q*b. Series pressure sums and parallel inverse-flow sums
/// are differentiated together with the loss pressure, not independently held
/// constant. No perturbed fan/network solve or matrix factorization is used.
///
/// The caller supplies the nominal retained flow, not a physical-uncertainty
/// corner. Its full nominal equations must close within `relative_tolerance`.
/// All member speed/flow domains must be interior. Unequal-slope knots and
/// nonunique parallel inverses refuse; redundant equal-slope knots are allowed.
/// Reductions use the producer's canonical bank order; returned controls use
/// the caller's order. Up to 64 banks and `max_curve_points` total points are
/// traversed, with cancellation between banks and bounded curve tiles.
///
/// This is an Estimated derivative at the supplied nominal point, NOT a
/// certified derivative interval, sensitivity of manufacturer-tolerance bands,
/// derivative of a changing topology/resistance, or finite-change guarantee.
/// The caller still composes air-capacity AND heat-transfer changes with a
/// complete thermal dual; this hydraulic factor alone is not dT/dspeed.
///
/// # Errors
/// Invalid policy/point, exceeded work, nonsmooth/out-of-domain members,
/// nominal residual failure, nonfinite arithmetic, or cancellation.
#[allow(clippy::too_many_arguments)]
pub fn operating_flow_speed_gradient(
    cx: &Cx<'_>, banks: &[FanBank], topology: FanArrangement,
    flow: VolumetricFlowRate, resistance: LossResistance,
    relative_tolerance: f64, max_curve_points: usize,
) -> Result<FanSpeedFlowGradient> {
    poll(cx)?;
    let q = flow.value();
    let r = resistance.value();
    if banks.is_empty() || banks.len() > 64 || !(q.is_finite() && q > 0.0
        && r.is_finite() && r > 0.0 && relative_tolerance.is_finite()
        && relative_tolerance > 0.0 && relative_tolerance < 1.0) {
        return Err(SpeedDerivativeError::Invalid("one to 64 banks, positive flow/resistance and tolerance in (0,1) required"));
    }
    let mut required = 0_usize;
    for bank in banks {
        poll(cx)?;
        required = required.checked_add(bank.curve().points().len())
            .ok_or(SpeedDerivativeError::Invalid("curve-point count overflow"))?;
    }
    if required > max_curve_points {
        return Err(SpeedDerivativeError::PointBudget { required, allowed: max_curve_points });
    }
    // One bank has no inverse-flow coupling, even when its fans are parallel.
    let topology = if banks.len() == 1 { FanArrangement::Series } else { topology };
    let pressure = checked(checked(r*q)?*q)?;
    if pressure <= 0.0 { return Err(SpeedDerivativeError::Invalid("unrepresentable loss pressure")); }
    let mut order: Vec<_> = (0..banks.len()).collect();
    order.sort_by(|&a, &b| compare_banks(&banks[a], &banks[b]));
    let mut numerators = vec![0.0; banks.len()];
    let mut sum = 0.0;
    let mut response = 0.0;
    for i in order {
        poll(cx)?;
        let bank = &banks[i];
        let speed_bounds = bank.curve().model_card().validity.bounds()["speed_ratio"];
        if bank.speed_ratio() <= speed_bounds.0 || bank.speed_ratio() >= speed_bounds.1 {
            return Err(SpeedDerivativeError::NonSmooth { bank: i });
        }
        let curve = effective_curve(bank);
        for (at, &(x,y)) in curve.iter().enumerate() {
            if at % 256 == 0 { poll(cx)?; }
            checked(x)?; checked(y)?;
            if x < 0.0 || y < 0.0 { return Err(SpeedDerivativeError::Invalid("invalid scaled member curve")); }
        }
        let lower = checked(bank.curve().admissible_min_flow().value()*bank.flow_factor())?;
        let upper = curve.last().expect("validated fan points").0;
        let (member_q, member_p) = match topology {
            FanArrangement::Series => (q, interpolate(&curve, q)),
            FanArrangement::Parallel => {
                let p_high = interpolate(&curve, lower);
                let p_low = curve.last().expect("validated fan points").1;
                if !(pressure > p_low && pressure < p_high) {
                    return Err(SpeedDerivativeError::NonSmooth { bank: i });
                }
                (interpolate_inverse(&curve, pressure), pressure)
            }
        };
        if !(member_q > lower && member_q < upper) {
            return Err(SpeedDerivativeError::NonSmooth { bank: i });
        }
        let b = slope(cx, &curve, member_q, i)?;
        match topology {
            FanArrangement::Series => {
                numerators[i] = checked(2.0*member_p - checked(q*b)?)?;
                response = checked(response - checked(q*b)?)?;
                sum = checked(sum + member_p)?;
            }
            FanArrangement::Parallel => {
                if b >= 0.0 { return Err(SpeedDerivativeError::NonSmooth { bank: i }); }
                let pressure_response = checked(-2.0*pressure/b)?;
                numerators[i] = checked(member_q + pressure_response)?;
                response = checked(response + pressure_response)?;
                sum = checked(sum + member_q)?;
            }
        }
    }
    let target = match topology { FanArrangement::Series => pressure, FanArrangement::Parallel => q };
    let relative = checked((sum/target - 1.0).abs())?;
    if relative >= relative_tolerance {
        return Err(SpeedDerivativeError::Residual { relative, tolerance: relative_tolerance });
    }
    let denominator = checked(match topology {
        FanArrangement::Series => 2.0*pressure + response,
        FanArrangement::Parallel => q + response,
    })?;
    if denominator <= 0.0 { return Err(SpeedDerivativeError::Invalid("singular hydraulic tangent")); }
    for value in &mut numerators { *value = checked(*value/denominator)?; }
    poll(cx)?;
    Ok(FanSpeedFlowGradient { log_flow_per_log_speed: numerators, relative_residual: relative })
}

fn slope(cx: &Cx<'_>, curve: &[(f64,f64)], q: f64, bank: usize) -> Result<f64> {
    let mut selected = None::<f64>;
    for (i, pair) in curve.windows(2).enumerate() {
        if i % 256 == 0 { poll(cx)?; }
        let [(q0,p0),(q1,p1)] = [pair[0],pair[1]];
        if q < q0 || q > q1 { continue; }
        let next = checked((p1-p0)/(q1-q0))?;
        if next > 0.0 || selected.is_some_and(|old| old != next) {
            return Err(SpeedDerivativeError::NonSmooth { bank });
        }
        selected = Some(next);
    }
    selected.ok_or(SpeedDerivativeError::NonSmooth { bank })
}
