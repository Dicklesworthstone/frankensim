//! Bounded scalar heat/work balance; the caller owns all caloric evaluations.

use fs_math::dd::Dd;

#[derive(Clone, Copy, Debug)]
pub(super) struct Response {
    pub work_j: Dd,
    pub work_derivative_kg: Dd,
}

/// A composing thermal owner supplies its original integrated energy residual
/// and analytic tangent; pressure work retains its independent native account.
pub(super) struct EnergyResidual {
    pub residual_j: Dd,
    pub derivative_kg: Dd,
    pub work_j: Dd,
}

#[derive(Clone, Copy, Debug)]
pub(super) struct Solution {
    pub increment_j_kg: f64,
    pub enthalpy_increment_j: f64,
    pub work_j: f64,
    pub residual_j: f64,
    pub evaluations: usize,
}

#[derive(Debug)]
pub(super) enum Failure<E> {
    Invalid(&'static str),
    Evaluation(E),
    Balance { evaluations: usize, residual_j: f64 },
}

/// Safeguarded Newton on an explicitly authored sign bracket. Only a joule
/// residual accepts a result; a small interval or stalled iterate never does.
/// The callback evaluates the original state/law and checks cancellation.
pub(super) fn solve<E>(
    mass_kg: Dd,
    heat_j: f64,
    bracket_j_kg: [f64; 2],
    tolerance_j: f64,
    max_evaluations: usize,
    mut response: impl FnMut(f64) -> Result<Response, E>,
) -> Result<Solution, Failure<E>> {
    if !heat_j.is_finite() {
        return Err(Failure::Invalid("invalid slab heat input"));
    }
    solve_residual(
        mass_kg,
        heat_j / number(mass_kg),
        bracket_j_kg,
        tolerance_j,
        max_evaluations,
        |increment| {
            let value = response(increment)?;
            Ok(EnergyResidual {
                residual_j: mass_kg * Dd::from_f64(increment) + value.work_j - Dd::from_f64(heat_j),
                derivative_kg: mass_kg + value.work_derivative_kg,
                work_j: value.work_j,
            })
        },
    )
}

/// The same bracket, acceptance and work budget for a composed source response.
/// A finite joule residual is the only convergence criterion.
pub(super) fn solve_residual<E>(
    mass_kg: Dd,
    guess_j_kg: f64,
    bracket_j_kg: [f64; 2],
    tolerance_j: f64,
    max_evaluations: usize,
    mut response: impl FnMut(f64) -> Result<EnergyResidual, E>,
) -> Result<Solution, Failure<E>> {
    let mass = number(mass_kg);
    let [mut lower, mut upper] = bracket_j_kg;
    if !mass.is_finite()
        || mass <= 0.0
        || !lower.is_finite()
        || !upper.is_finite()
        || lower >= upper
        || !tolerance_j.is_finite()
        || tolerance_j < 0.0
        || max_evaluations < 2
    {
        return Err(Failure::Invalid("invalid slab heat bracket or work budget"));
    }
    let mut evaluations = 0;
    let mut evaluate = |increment: f64| -> Result<(Solution, f64), Failure<E>> {
        let value = response(increment).map_err(Failure::Evaluation)?;
        evaluations += 1;
        let enthalpy = mass_kg * Dd::from_f64(increment);
        let residual = number(value.residual_j);
        let derivative = number(value.derivative_kg);
        let enthalpy = number(enthalpy);
        let work = number(value.work_j);
        if [residual, derivative, enthalpy, work]
            .iter()
            .any(|v| !v.is_finite())
        {
            return Err(Failure::Invalid(
                "slab heat/work arithmetic is unrepresentable",
            ));
        }
        Ok((
            Solution {
                increment_j_kg: increment,
                enthalpy_increment_j: enthalpy,
                work_j: work,
                residual_j: residual,
                evaluations,
            },
            derivative,
        ))
    };
    let (left, _) = evaluate(lower)?;
    if left.residual_j.abs() <= tolerance_j {
        return Ok(left);
    }
    let (right, _) = evaluate(upper)?;
    if right.residual_j.abs() <= tolerance_j {
        return Ok(right);
    }
    let mut best = if left.residual_j.abs() < right.residual_j.abs() {
        left
    } else {
        right
    };
    if left.residual_j.is_sign_positive() == right.residual_j.is_sign_positive() {
        return Err(Failure::Balance {
            evaluations: 2,
            residual_j: best.residual_j,
        });
    }
    let mut lower_positive = left.residual_j.is_sign_positive();
    let mut x = if guess_j_kg > lower && guess_j_kg < upper {
        guess_j_kg
    } else {
        0.5 * lower + 0.5 * upper
    };
    for count in 2..max_evaluations {
        if x <= lower || x >= upper {
            return Err(Failure::Balance {
                evaluations: count,
                residual_j: best.residual_j,
            });
        }
        let (trial, derivative) = evaluate(x)?;
        if trial.residual_j.abs() < best.residual_j.abs() {
            best = trial;
        }
        if trial.residual_j.abs() <= tolerance_j {
            return Ok(trial);
        }
        if trial.residual_j.is_sign_positive() == lower_positive {
            lower = x;
            lower_positive = trial.residual_j.is_sign_positive();
        } else {
            upper = x;
        }
        let newton = x - trial.residual_j / derivative;
        // Keep a definite bracket contraction even for a nearly horizontal
        // or nonmonotone source branch. No derivative sign is presumed.
        let inner_lower = 0.9 * lower + 0.1 * upper;
        let inner_upper = 0.1 * lower + 0.9 * upper;
        x = if newton.is_finite() && newton > inner_lower && newton < inner_upper {
            newton
        } else {
            0.5 * lower + 0.5 * upper
        };
    }
    Err(Failure::Balance {
        evaluations: max_evaluations,
        residual_j: best.residual_j,
    })
}

/// Existing backward-Euler pressure work for a uniform affine source:
/// B*u*(u-u0), with B the ORIGINAL native face inertia at unit surface speed.
/// This function is only algebra; the native projection measures B and checks
/// the accepted work independently.
pub(super) fn pressure_response(
    inertia_kg: f64,
    surface_speed_m_s: f64,
    incoming_surface_speed_m_s: f64,
    surface_speed_derivative: f64,
) -> Response {
    let b = Dd::from_f64(inertia_kg);
    let u = Dd::from_f64(surface_speed_m_s);
    let old = Dd::from_f64(incoming_surface_speed_m_s);
    Response {
        work_j: b * u * (u - old),
        work_derivative_kg: b * (u + u - old) * Dd::from_f64(surface_speed_derivative),
    }
}

fn number(value: Dd) -> f64 {
    value.hi + value.lo
}

#[cfg(test)]
mod tests {
    use super::*;

    fn run(heat: f64, incoming: f64, bracket: [f64; 2]) -> Solution {
        solve(Dd::from_f64(2.0), heat, bracket, 1e-12, 80, |h| {
            Ok::<_, ()>(pressure_response(3.0, 0.1 * h, incoming, 0.1))
        })
        .unwrap()
    }

    #[test]
    fn heating_pays_for_actual_pressure_work() {
        // 2*10 + 3*1^2 = 23.
        let solved = run(23.0, 0.0, [0.0, 20.0]);
        assert!((solved.increment_j_kg - 10.0).abs() < 1e-12);
        assert!((solved.enthalpy_increment_j - 20.0).abs() < 1e-12);
        assert!((solved.work_j - 3.0).abs() < 1e-12);
        assert!(solved.evaluations <= 80);
    }

    #[test]
    fn cooling_supplies_contraction_motion_without_invented_heat() {
        // 2*(-10) + 3*(-1)^2 = -17.
        let solved = run(-17.0, 0.0, [-20.0, 0.0]);
        assert!((solved.increment_j_kg + 10.0).abs() < 1e-12);
        assert!((solved.work_j - 3.0).abs() < 1e-12);
    }

    #[test]
    fn retained_motion_can_return_signed_pressure_work() {
        // u0=2, u=1: 2*10 + 3*1*(1-2)=17.
        let solved = run(17.0, 2.0, [0.0, 20.0]);
        assert!((solved.increment_j_kg - 10.0).abs() < 1e-12);
        assert!((solved.work_j + 3.0).abs() < 1e-12);
    }

    #[test]
    fn flat_density_reduces_to_original_mass_times_enthalpy() {
        let solved = solve(Dd::from_f64(4.0), 12.0, [-5.0, 10.0], 0.0, 8, |_| {
            Ok::<_, ()>(pressure_response(7.0, 0.0, 3.0, 0.0))
        })
        .unwrap();
        assert_eq!(solved.increment_j_kg, 3.0);
        assert_eq!(solved.work_j, 0.0);
        assert_eq!(solved.residual_j, 0.0);
    }

    #[test]
    fn endpoint_and_zero_source_roots_are_admitted() {
        assert_eq!(run(0.0, 0.0, [0.0, 10.0]).increment_j_kg, 0.0);
        let upper = run(23.0, 0.0, [0.0, 10.0]);
        assert_eq!(upper.increment_j_kg, 10.0);
        assert_eq!(upper.evaluations, 2);
    }

    #[test]
    fn explicit_bracket_selects_a_decreasing_response_branch() {
        let solved = solve(Dd::from_f64(1.0), 3.0, [-2.0, 2.0], 0.0, 80, |h| {
            Ok::<_, ()>(Response {
                work_j: Dd::from_f64(3.0 - 2.0 * h),
                work_derivative_kg: Dd::from_f64(-2.0),
            })
        })
        .unwrap();
        assert_eq!(solved.increment_j_kg, 0.0);
    }

    #[test]
    fn unbracketed_and_exhausted_solves_publish_no_candidate() {
        assert!(matches!(
            solve(Dd::from_f64(2.0), 100.0, [0.0, 1.0], 1e-12, 80, |h| Ok::<
                _,
                (),
            >(
                pressure_response(3.0, h, 0.0, 1.0)
            )),
            Err(Failure::Balance { evaluations: 2, .. })
        ));
        assert!(matches!(
            solve(Dd::from_f64(2.0), 23.0, [0.0, 20.0], 1e-12, 2, |h| Ok::<
                _,
                (),
            >(
                pressure_response(3.0, 0.1 * h, 0.0, 0.1)
            )),
            Err(Failure::Balance { evaluations: 2, .. })
        ));
    }

    #[test]
    fn source_or_cancellation_refusal_is_not_replaced_by_a_trial() {
        let mut count = 0;
        let error = solve(Dd::from_f64(2.0), 23.0, [0.0, 20.0], 1e-12, 80, |h| {
            count += 1;
            if count == 3 {
                Err("cancelled")
            } else {
                Ok(pressure_response(3.0, 0.1 * h, 0.0, 0.1))
            }
        });
        assert!(matches!(error, Err(Failure::Evaluation("cancelled"))));
        assert_eq!(count, 3);
    }

    #[test]
    fn product_residual_keeps_small_heat_after_large_cancellation() {
        let mass = Dd::from_f64(1e16) + Dd::from_f64(1.0);
        let solved = solve(mass, 1.0, [1.0, 2.0], 0.0, 2, |_| {
            Ok::<_, ()>(Response {
                work_j: Dd::from_f64(-1e16),
                work_derivative_kg: Dd::ZERO,
            })
        })
        .unwrap();
        assert_eq!(solved.increment_j_kg, 1.0);
        assert_eq!(solved.residual_j, 0.0);
    }

    #[test]
    fn unrepresentable_arithmetic_and_invalid_budgets_refuse() {
        assert!(matches!(
            solve(Dd::from_f64(0.0), 1.0, [0.0, 1.0], 0.0, 2, |_| Ok::<_, ()>(
                pressure_response(1.0, 0.0, 0.0, 0.0)
            )),
            Err(Failure::Invalid(_))
        ));
        assert!(matches!(
            solve(Dd::from_f64(1.0), 1.0, [0.0, 1.0], 0.0, 2, |_| Ok::<_, ()>(
                pressure_response(f64::MAX, f64::MAX, 0.0, 1.0)
            )),
            Err(Failure::Invalid(_))
        ));
    }
}
