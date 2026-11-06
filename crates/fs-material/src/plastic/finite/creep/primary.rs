//! Primary resistance in the native multiplicative finite return.
//!
//! Reuse the original scalar h(alpha) and logarithmic creep root. Finite
//! kinematics, isochoric metric reconstruction and dP/dF stay in FiniteJ2.
//! The inverse kinetic stress is q(c) = s0[eta*c/(dt*s0*h(old+c))]^(1/n).
//! Its integral, not q(c)*c or the secondary-only n/(n+1) fraction, supplies
//! the reduced Newton potential. No recoverable energy is assigned to h.

use crate::creep::{CreepError, StrainHardeningCreep};
use crate::plastic::finite::rate::FiniteOverstressJ2;
use crate::plastic::finite::{FiniteJ2Error, FiniteJ2State, RadialReturn, ReturnMap};
use fs_math::det;

fn kinetic_error(error: CreepError) -> FiniteJ2Error {
    match error {
        CreepError::Invalid(reason) => FiniteJ2Error::Numerical(reason),
        _ => FiniteJ2Error::Numerical("primary creep kinetic evaluation failed"),
    }
}

pub(in crate::plastic::finite) fn endpoint_stress(
    law: &FiniteOverstressJ2,
    increment: f64,
    old_alpha: f64,
) -> Result<f64, FiniteJ2Error> {
    if !old_alpha.is_finite() || old_alpha < 0.0 || !increment.is_finite() || increment < 0.0 {
        return Err(FiniteJ2Error::Parameters(
            "finite nonnegative creep history and increment required",
        ));
    }
    let Some(primary) = law.active_primary() else {
        return law.endpoint_overstress_pa(increment);
    };
    if increment == 0.0 {
        return Ok(0.0);
    }
    let (log_q, _) = inverse_rate(law, primary, increment, old_alpha)?;
    let q = det::exp(log_q);
    if !q.is_finite() || q <= 0.0 {
        return Err(FiniteJ2Error::Numerical(
            "primary endpoint stress is not representable",
        ));
    }
    Ok(q)
}

// Log stress and beta=-c*d(log h)/dalpha; dq/dc=q*(1+beta)/(n*c).
fn inverse_rate(
    law: &FiniteOverstressJ2,
    primary: StrainHardeningCreep,
    increment: f64,
    old_alpha: f64,
) -> Result<(f64, f64), FiniteJ2Error> {
    let alpha = old_alpha + increment;
    if !alpha.is_finite() {
        return Err(FiniteJ2Error::Numerical(
            "primary creep accumulation is not representable",
        ));
    }
    let (log_h, weight) = primary.value(alpha).map_err(kinetic_error)?;
    let p = law.parameters();
    let log_s0 = det::ln(p.stress_scale_pa);
    let log_q = log_s0
        + (det::ln(p.viscosity_pa_s) - det::ln(law.duration_s()) + det::ln(increment)
            - log_s0
            - log_h)
            / p.exponent;
    let beta = primary.sensitivity(alpha, increment, weight);
    Ok((log_q, beta))
}

// Shared primary relaxation; only the stress measure and elastic stiffness
// differ from the additive owner. The scalar rate law is not duplicated.
fn relaxation(
    law: &FiniteOverstressJ2,
    primary: StrainHardeningCreep,
    trial: f64,
    old_alpha: f64,
) -> Result<(f64, f64, f64, f64), FiniteJ2Error> {
    let h = 3.0 * law.backbone().parameters().shear_modulus_pa;
    let log_h = det::ln(h);
    let p = law.parameters();
    let returned = primary
        .relaxation(
            trial,
            log_h + det::ln(law.duration_s()) - det::ln(p.viscosity_pa_s),
            log_h,
            old_alpha,
            p.exponent,
            det::ln(p.stress_scale_pa),
        )
        .map_err(kinetic_error)?;
    if trial > 0.0 && returned.delta_alpha <= 0.0 {
        return Err(FiniteJ2Error::Numerical(
            "positive finite creep increment is below representable range",
        ));
    }
    let residual = h.mul_add(returned.delta_alpha, returned.q) - trial;
    if !residual.is_finite() || residual.abs() > 512.0 * f64::EPSILON * trial {
        return Err(FiniteJ2Error::Numerical(
            "finite primary return violates stress balance",
        ));
    }
    Ok((
        returned.delta_alpha,
        returned.q,
        returned.ratio,
        returned.radial_correction,
    ))
}

pub(in crate::plastic::finite) fn standalone(
    law: &FiniteOverstressJ2,
    primary: StrainHardeningCreep,
    f: &[f64; 9],
    old: &FiniteJ2State,
) -> Result<ReturnMap, FiniteJ2Error> {
    let backbone = law.backbone();
    if backbone.maximum_creep_alpha.is_some() {
        return Err(FiniteJ2Error::Parameters(
            "combined creep history requires its simultaneous return",
        ));
    }
    let trial = backbone.elastic_endpoint(f, old)?;
    let stress = trial
        .deviator
        .map(|x| 2.0 * backbone.parameters.shear_modulus_pa * x);
    let q_trial = det::sqrt(1.5) * super::super::norm(&stress);
    if !q_trial.is_finite() {
        return Err(FiniteJ2Error::Numerical(
            "primary trial stress is not finite",
        ));
    }
    let (creep, q, scale, correction) = relaxation(law, primary, q_trial, old.alpha)?;
    let heat = q * creep;
    let potential = rate_potential(law, primary, old.alpha, creep, q)?;
    backbone.finish_return(
        f,
        old,
        trial,
        q_trial,
        RadialReturn {
            // The original standalone finite history names its sole flow alpha.
            plastic_increment: creep,
            creep_increment: 0.0,
            scale,
            tangent_correction: correction,
            plastic_dissipation: heat,
            creep_dissipation: 0.0,
            primary_potential: 0.0,
            rate_potential: potential,
        },
    )
}

pub(in crate::plastic::finite) fn combined(
    law: &FiniteOverstressJ2,
    primary: StrainHardeningCreep,
    trial: f64,
    plastic_alpha: f64,
    creep_alpha: f64,
) -> Result<RadialReturn, FiniteJ2Error> {
    let backbone = law.backbone();
    let h = 3.0 * backbone.parameters.shear_modulus_pa;
    let yield_old = backbone.yield_radius(plastic_alpha);
    let gap = trial - yield_old;
    let resolved_gap = gap - 64.0 * f64::EPSILON * trial.max(yield_old);
    let inactive = resolved_gap <= 0.0 || {
        // Test the actual creep history at the yield crossing, not old history.
        let log_hardening = primary
            .value(creep_alpha + gap / h)
            .map_err(kinetic_error)?
            .0;
        let log_flow = det::ln(h) + law.log_increment_at_stress(yield_old) + log_hardening;
        let log_gap = det::ln(resolved_gap);
        let tie = 8.0 * f64::EPSILON * (1.0 + log_flow.abs() + log_gap.abs());
        log_flow >= log_gap || log_flow.is_finite() && log_gap - log_flow <= tie
    };
    let (creep, q, plastic, scale, correction) = if inactive {
        let (creep, q, scale, correction) = relaxation(law, primary, trial, creep_alpha)?;
        (creep, q, 0.0, scale, correction)
    } else {
        let (creep, q, plastic, radial) =
            active_combined(law, primary, trial, plastic_alpha, creep_alpha)?;
        let scale = q / trial;
        (creep, q, plastic, scale, radial - scale)
    };
    let balance = h.mul_add(creep + plastic, q) - trial;
    let radius = backbone.yield_radius(plastic_alpha + plastic);
    let yield_error = if plastic > 0.0 {
        (q - radius).abs()
    } else {
        (q - radius).max(0.0)
    };
    let tolerance = 512.0 * f64::EPSILON * trial.max(yield_old);
    if !balance.is_finite() || balance.abs() > tolerance || yield_error > tolerance {
        return Err(FiniteJ2Error::Numerical(
            "primary joint return violates stress balance or plastic yield",
        ));
    }
    let plastic_heat = backbone.parameters.yield_stress_pa * plastic;
    Ok(RadialReturn {
        plastic_increment: plastic,
        creep_increment: creep,
        scale,
        tangent_correction: correction,
        plastic_dissipation: plastic_heat,
        creep_dissipation: q * creep,
        primary_potential: plastic_heat,
        rate_potential: rate_potential(law, primary, creep_alpha, creep, q)?,
    })
}

// Eliminate plastic flow: p=(trial-q(c))/h-c. The remaining yield residual is
// strictly increasing in log(c), including perfect, linear and Voce hardening.
// This uses one scalar solve; there is no nested implicit creep solve at yield.
fn active_combined(
    law: &FiniteOverstressJ2,
    primary: StrainHardeningCreep,
    trial: f64,
    plastic_alpha: f64,
    creep_alpha: f64,
) -> Result<(f64, f64, f64, f64), FiniteJ2Error> {
    let backbone = law.backbone();
    let h = 3.0 * backbone.parameters.shear_modulus_pa;
    let yield_old = backbone.yield_radius(plastic_alpha);
    let maximum = (trial - yield_old) / h;
    if !maximum.is_finite() || maximum <= 0.0 {
        return Err(FiniteJ2Error::Numerical(
            "active primary creep bound is not representable",
        ));
    }
    let mut lower = det::ln(f64::from_bits(1));
    let mut upper = det::ln(maximum);
    if upper < lower {
        return Err(FiniteJ2Error::Numerical(
            "active primary creep increment is below representable range",
        ));
    }
    let initial = law.log_increment_at_stress(yield_old)
        + primary
            .value(creep_alpha + maximum)
            .map_err(kinetic_error)?
            .0;
    let mut y = initial.clamp(lower, upper);
    for _ in 0..96 {
        let creep = det::exp(y);
        let (log_q, beta) = inverse_rate(law, primary, creep, creep_alpha)?;
        // Upper bracket samples may exceed trial stress; they need no power,
        // plastic history or hardening-energy construction before refinement.
        if log_q >= det::ln(trial) {
            upper = y;
            y = 0.5 * lower + 0.5 * upper;
            continue;
        }
        let q = det::exp(log_q);
        let plastic = (trial - h.mul_add(creep, q)) / h;
        if plastic < 0.0 {
            upper = y;
            y = 0.5 * lower + 0.5 * upper;
            continue;
        }
        let alpha = plastic_alpha + plastic;
        let residual = q - backbone.yield_radius(alpha);
        let hardening = backbone.hardening_tangent(alpha);
        if residual.is_finite() && residual.abs() <= 128.0 * f64::EPSILON * trial {
            if q <= 0.0 || creep <= 0.0 || plastic <= 0.0 {
                return Err(FiniteJ2Error::Numerical(
                    "active primary joint increments are not representable",
                ));
            }
            let radial = if hardening == 0.0 {
                0.0
            } else {
                // dq/dtrial=1/(1+h/H+h/v), v=dq/dc. Keep all stiffness
                // ratios logarithmic, including subnormal positive H.
                let log_v = log_q - y + det::ln(1.0 + beta) - det::ln(law.parameters().exponent);
                let a = det::ln(h) - det::ln(hardening);
                let b = det::ln(h) - log_v;
                let shift = 0.0_f64.max(a).max(b);
                det::exp(-shift) / (det::exp(-shift) + det::exp(a - shift) + det::exp(b - shift))
            };
            return Ok((creep, q, plastic, radial));
        }
        if residual > 0.0 {
            upper = y;
        } else {
            lower = y;
        }
        let dq = q * ((1.0 + beta) / law.parameters().exponent);
        let derivative = dq + hardening * (creep + dq / h);
        let candidate = y - residual / derivative;
        y = if candidate.is_finite() && candidate > lower && candidate < upper && candidate != y {
            candidate
        } else {
            0.5 * lower + 0.5 * upper
        };
    }
    Err(FiniteJ2Error::Numerical(
        "bounded primary joint return did not converge",
    ))
}

fn rate_potential(
    law: &FiniteOverstressJ2,
    primary: StrainHardeningCreep,
    old_alpha: f64,
    increment: f64,
    q: f64,
) -> Result<f64, FiniteJ2Error> {
    if increment == 0.0 {
        return Ok(0.0);
    }
    let n = law.parameters().exponent;
    let alpha = old_alpha + increment;
    let (endpoint_log_h, weight) = primary.value(alpha).map_err(kinetic_error)?;
    // Exact cubic primitive for the common shifted reciprocal primary law.
    let fraction = if n == 1.0 && primary.exponent == 1.0 && primary.residual_rate_fraction == 0.0 {
        0.5 - primary.sensitivity(alpha, increment, weight) / 6.0
    } else {
        // c=increment*u^4 removes the inverse-rate endpoint singularity.
        // The normalized integrand is finite, nonnegative and bounded by 4.
        let integrand = |u: f64| -> Result<f64, FiniteJ2Error> {
            if u == 0.0 {
                return Ok(0.0);
            }
            let u2 = u * u;
            let log_h = primary
                .value(old_alpha + increment * u2 * u2)
                .map_err(kinetic_error)?
                .0;
            Ok(4.0 * u2 * u * det::exp((4.0 * det::ln(u) + endpoint_log_h - log_h) / n))
        };
        integrate_potential(integrand)?
    };
    let potential = (q * increment) * fraction;
    if !potential.is_finite()
        || potential <= 0.0
        || fraction > (1.0 / (1.0 + 1.0 / n)) * (1.0 + 512.0 * f64::EPSILON)
    {
        return Err(FiniteJ2Error::Numerical(
            "primary creep potential is not representable",
        ));
    }
    Ok(potential)
}

// Embedded Gauss 7 / Kronrod 15 rule (QUADPACK dqk15 nodes and weights).
// The error difference is an estimate, not an interval certificate. Fixed
// stack, depth and work caps bound every evaluation; unresolved sharp endpoint
// layers refuse instead of returning a false zero from two missed rules.
fn integrate_potential(
    f: impl Fn(f64) -> Result<f64, FiniteJ2Error>,
) -> Result<f64, FiniteJ2Error> {
    const X: [f64; 8] = [
        0.9914553711208126,
        0.9491079123427585,
        0.8648644233597691,
        0.7415311855993945,
        0.5860872354676911,
        0.4058451513773972,
        0.2077849550078985,
        0.0,
    ];
    const K: [f64; 8] = [
        0.02293532201052922,
        0.06309209262997855,
        0.1047900103222502,
        0.1406532597155259,
        0.1690047266392679,
        0.1903505780647854,
        0.2044329400752989,
        0.2094821410847278,
    ];
    const G: [f64; 4] = [
        0.1294849661688697,
        0.2797053914892767,
        0.3818300505051189,
        0.4179591836734694,
    ];
    let mut stack = [(0.0, 0.0, 0_u8); 30];
    stack[0] = (0.0, 1.0, 0);
    let mut count = 1;
    let mut total = 0.0;
    let mut integral_scale = 0.0_f64;
    for _ in 0..512 {
        if count == 0 {
            return Ok(total);
        }
        count -= 1;
        let (left, right, depth) = stack[count];
        let half = 0.5 * (right - left);
        let middle = 0.5 * left + 0.5 * right;
        let center = f(middle)?;
        let mut kronrod = K[7] * center;
        let mut gauss = G[3] * center;
        for i in 0..7 {
            let pair = f(middle - half * X[i])? + f(middle + half * X[i])?;
            kronrod += K[i] * pair;
            if i % 2 == 1 {
                gauss += G[i / 2] * pair;
            }
        }
        kronrod *= half;
        gauss *= half;
        integral_scale = integral_scale.max(kronrod).max(total);
        let resolved_endpoint = kronrod >= (right - left) * f(right)? * 1e-5;
        // Allocate an absolute budget from the whole integral as well as a
        // relative panel budget. The transformed integrand can have a weak
        // noninteger derivative at zero: its local relative quadrature error
        // does not decrease under subdivision even as its contribution vanishes.
        let tolerance = 2e-13 * (kronrod + (right - left) * integral_scale);
        if resolved_endpoint && (kronrod - gauss).abs() <= tolerance {
            total += kronrod;
        } else {
            if depth >= 28 || middle == left || middle == right || count + 2 > stack.len() {
                return Err(FiniteJ2Error::Numerical(
                    "primary potential quadrature reached its resolution limit",
                ));
            }
            stack[count] = (left, middle, depth + 1);
            stack[count + 1] = (middle, right, depth + 1);
            count += 2;
        }
    }
    Err(FiniteJ2Error::Numerical(
        "primary potential quadrature exhausted its bounded work",
    ))
}
