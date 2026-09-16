//! Fan-speed chain rule for the admitted single-fan quadratic-loss network.
//!
//! Fan affinity gives p(q,s)=s^2 p(q/s,1); every passive edge has drop R q|q|.
//! Thus q_e(s)=s q_e(1) and p_v(s)=s^2 p_v(1), including reverse-oriented
//! edges. Mixing fractions remain fixed, capacity rates and Reynolds numbers
//! scale with s, and dln(h)/dln(s)=dln(Nu)/dln(Re) at fixed fluid/geometry.
//! This is NOT a derivative for arbitrary pressure controls, flow-dependent
//! resistance, fan heating, multiple independent fans or topology switches.

use super::*;
use fs_airflow::graph::thermal::coupled_transport::sensitivity::CoupledObjective;
use fs_convection::{CorrelationId, NusseltEvaluation};
use fs_math::det;

/// Retain old conditional controls together with the new total speed control.
/// Deref preserves existing internal consumers of the conditional gradient.
#[derive(Debug)]
pub(super) struct CoolingGradient {
    thermal: CoupledGradient,
    speed: Option<SpeedGradient>,
}

impl std::ops::Deref for CoolingGradient {
    type Target = CoupledGradient;
    fn deref(&self) -> &Self::Target { &self.thermal }
}

#[derive(Debug)]
enum SpeedGradient {
    Available { flow: f64, convection: f64, total: f64 },
    Unavailable,
}

impl CoolingGradient {
    pub(super) fn speed_json(&self) -> Result<String> {
        match &self.speed {
            None => Ok("null".into()),
            Some(SpeedGradient::Unavailable) => Ok(
                "{\"status\":\"unavailable\",\"reason\":\"a convection card has no admitted Reynolds derivative here; conditional thermal gradients are retained\"}".into(),
            ),
            Some(SpeedGradient::Available { flow, convection, total }) => Ok(format!(
                "{{\"status\":\"available\",\"method\":\"fan-affinity-coupled-adjoint\",\"dobjective_dlog_speed_ratio_k\":{},\"capacity_contribution_k\":{},\"convection_contribution_k\":{},\"scope\":\"steady selected temperature objective; single affinity-scaled fan bank and fixed quadratic losses; full solid/air feedback and smooth-card Reynolds response; fixed geometry, material laws, source loads and fluid properties; local derivative only, not a finite speed-change or safety bound\"}}",
                num(*total)?, num(*flow)?, num(*convection)?,
            )),
        }
    }
}

/// Use exactly one coupled adjoint. The common-flow extension needs one extra
/// air-only reverse sweep, not additional perturbed hydraulic or solid solves.
#[allow(clippy::too_many_arguments)]
pub(super) fn pullback(
    request: &Request,
    cx: &Cx<'_>,
    binding: &CoupledLinearization<'_, '_>,
    objective: &CoupledObjective,
    names: &[&str],
    derivations: &[convection::Derived],
) -> Result<CoolingGradient> {
    let config = InterfaceSolveConfig {
        max_iterations: request.limits.derivative,
        absolute_tolerance: request.limits.relative,
        relative_tolerance: request.limits.relative,
        relaxation: request.limits.relaxation,
    };
    if request.fan.is_none() {
        return Ok(CoolingGradient {
            thermal: binding.pullback_iqn(cx, objective, config, acceleration::POLICY).map_err(producer)?,
            speed: None,
        });
    }
    let mut slopes = BTreeMap::new();
    for derived in derivations {
        poll(cx)?;
        let Some(slope) = reynolds_elasticity(&derived.nu)? else {
            return Ok(CoolingGradient {
                thermal: binding.pullback_iqn(cx, objective, config, acceleration::POLICY).map_err(producer)?,
                speed: Some(SpeedGradient::Unavailable),
            });
        };
        slopes.insert(derived.surface.as_str(), slope);
    }
    let result = binding.pullback_flow_scale_iqn(cx, objective, config, acceleration::POLICY)
        .map_err(producer)?;
    if names.len() != result.thermal.log_htc.len() {
        return Err(bad("fan derivative requires the exact coupled region ordering"));
    }
    let mut convection = 0.0;
    for (name, gradient) in names.iter().zip(&result.thermal.log_htc) {
        poll(cx)?;
        // Declared scalar h is independent of speed, not an unknown slope.
        let slope = slopes.get(name).copied().unwrap_or(0.0);
        convection = checked(convection + checked(gradient * slope)?)?;
    }
    let total = checked(result.log_flow_scale + convection)?;
    poll(cx)?;
    Ok(CoolingGradient {
        thermal: result.thermal,
        speed: Some(SpeedGradient::Available { flow: result.log_flow_scale, convection, total }),
    })
}

/// Differentiate the admitted formula at the actual producer's group point.
/// Nu itself comes from that producer, not a second h evaluation. The C0
/// rectangular developing table is deliberately withheld: choosing a chord
/// at a knot would invent a two-sided derivative. Other unhandled cards also
/// remain unavailable, never silently assigned zero speed response.
fn reynolds_elasticity(evaluation: &NusseltEvaluation) -> Result<Option<f64>> {
    let group = |name: &str| evaluation.groups().get(name).copied()
        .ok_or_else(|| bad(format!("convection derivative lacks producer group {name}")));
    let nu = evaluation.evidence().value;
    let value = match evaluation.card().id {
        CorrelationId::CircularDuctLaminarCwt | CorrelationId::RectangularDuctLaminarCwt => 0.0,
        CorrelationId::CircularDuctHausen => {
            let gz = group("Gz")?;
            let b = 0.04 * det::pow(gz, 2.0 / 3.0);
            let enhancement = 0.0668 * gz / (1.0 + b);
            enhancement / nu * (1.0 + b / 3.0) / (1.0 + b)
        }
        CorrelationId::DittusBoelter => 0.8,
        CorrelationId::Gnielinski => {
            let re = group("Re")?;
            let pr = group("Pr")?;
            let z = 0.79 * det::ln(re) - 1.64;
            let root_f_over_8 = 1.0 / (det::sqrt(8.0) * z);
            let correction = 12.7 * root_f_over_8 * (det::pow(pr, 2.0 / 3.0) - 1.0);
            let log_f_slope = -1.58 / z;
            re / (re - 1000.0) + log_f_slope * (1.0 - 0.5 * correction / (1.0 + correction))
        }
        _ => return Ok(None),
    };
    Ok(Some(checked(value)?))
}

fn checked(value: f64) -> Result<f64> {
    if value.is_finite() { Ok(value) } else { Err(producer("nonfinite fan-speed sensitivity")) }
}

#[cfg(test)]
mod tests {
    use super::*;
    use fs_convection::{CorrelationInputs, ThermalDirection, evaluate};

    #[test]
    fn analytic_reynolds_slopes_match_the_real_card_evaluator() {
        let eps = 1e-5_f64;
        for (card, re, pr) in [
            (CorrelationId::CircularDuctLaminarCwt, 1000.0, 0.72),
            (CorrelationId::RectangularDuctLaminarCwt, 1000.0, 0.72),
            (CorrelationId::CircularDuctHausen, 100.0, 0.72),
            (CorrelationId::CircularDuctHausen, 1000.0, 0.72),
            (CorrelationId::DittusBoelter, 50_000.0, 0.72),
            (CorrelationId::Gnielinski, 20_000.0, 0.72),
            (CorrelationId::Gnielinski, 20_000.0, 7.0),
        ] {
            let run = |r| evaluate(card, CorrelationInputs::forced(r, pr)
                .with_length_ratio(1000.0).with_aspect_ratio(0.5)
                .with_direction(ThermalDirection::HeatingFluid)).unwrap();
            let nominal = run(re);
            let plus = run(re * eps.exp()).evidence().value;
            let minus = run(re * (-eps).exp()).evidence().value;
            let expected = (plus.ln() - minus.ln()) / (2.0 * eps);
            let actual = reynolds_elasticity(&nominal).unwrap().unwrap();
            assert!((actual - expected).abs() < 1e-8, "{card:?}: {actual} != {expected}");
        }
    }

    #[test]
    fn source_table_does_not_acquire_a_fabricated_smooth_derivative() {
        let table = evaluate(CorrelationId::RectangularDuctLaminarCwtDevelopingPr072,
            CorrelationInputs::forced(1000.0, 0.72).with_length_ratio(24.0)
                .with_aspect_ratio(0.5)).unwrap();
        assert_eq!(reynolds_elasticity(&table).unwrap(), None);
    }
}
