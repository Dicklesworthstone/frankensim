//! Cancellation-safe scalar spectra for passive multi-arm Maxwell materials.
//!
//! The underdamped pair is recovered from positive physical loss sums, rather
//! than subtracting relaxation poles from a trace. Every material state remains:
//! distinct active rates have interlacing real poles, and repeated or uncoupled
//! rates retain their original dark-mode multiplicities. No Prony arm is fitted,
//! truncated, or replaced by a frequency-sampled damping ratio.

use super::{MaxwellDampingBranch, ModalError, positive_ratio};
use crate::{quadratic_fraction, quadratic_power_scale};
use fs_math::{c64::C64, dd::Dd, det};

const MAX_ITERATIONS: usize = 64;
const MAX_STIFFNESS_RATIO: f64 = 0.0625;

struct Arm {
    rate: f64,
    gamma: (f64, i32),
    ratio: f64,
}

// Source products stay factored until the final physical result is rounded.
// This also preserves a finite decay when gamma/K or a pole shift underflows.
fn ratio_parts(numerator: &[f64], denominator: &[f64]) -> (f64, i32) {
    let mut mantissa = 1.0;
    let mut exponent = 0;
    for &value in numerator {
        if value == 0.0 {
            return (0.0, 0);
        }
        let (fraction, power) = quadratic_fraction(value);
        mantissa *= fraction;
        exponent += power;
    }
    for &value in denominator {
        let (fraction, power) = quadratic_fraction(value);
        mantissa /= fraction;
        exponent -= power;
    }
    (mantissa, exponent)
}

#[derive(Default)]
struct ScaledSum {
    value: Dd,
    exponent: i32,
}
impl ScaledSum {
    fn add(&mut self, (value, exponent): (f64, i32)) {
        if value == 0.0 {
            return;
        }
        if self.value == Dd::ZERO {
            self.value = Dd::from_f64(value);
            self.exponent = exponent;
        } else if exponent > self.exponent {
            self.value = Dd::from_pair(
                quadratic_power_scale(self.value.hi, self.exponent - exponent),
                quadratic_power_scale(self.value.lo, self.exponent - exponent),
            ) + Dd::from_f64(value);
            self.exponent = exponent;
        } else {
            self.value =
                self.value + Dd::from_f64(quadratic_power_scale(value, exponent - self.exponent));
        }
    }

    fn physical(&self, power: i32) -> f64 {
        quadratic_power_scale(self.value.to_f64(), self.exponent + power)
    }

    fn sign(&self) -> f64 {
        self.value.to_f64()
    }
}

fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), ModalError> {
    if cancelled() {
        Err(ModalError::Cancelled)
    } else {
        Ok(())
    }
}

fn unresolved() -> ModalError {
    ModalError::QuadraticEig {
        detail: "scalar multi-arm Maxwell poles require representable physical rates and decay"
            .into(),
    }
}

// H(u) = u - C/m + K/(m*u) + sum(gamma_j/(u-r_j)).
// Its sign locates a real pole -u between distinct active source rates. A
// common binary scale and DD summation avoid overflowing stiffness/u or a
// near-pole quotient before the terms cancel. No source coefficient is dropped.
fn real_sign(
    u: f64,
    damping: (f64, i32),
    stiffness: f64,
    mass: f64,
    arms: &[Arm],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<f64, ModalError> {
    let mut sum = ScaledSum::default();
    sum.add(quadratic_fraction(u));
    sum.add((-damping.0, damping.1));
    sum.add(ratio_parts(&[stiffness], &[mass, u]));
    for (i, arm) in arms.iter().enumerate() {
        if i % 64 == 0 {
            poll(cancelled)?;
        }
        if arm.gamma.0 == 0.0 {
            continue;
        }
        let difference = u - arm.rate;
        // Only strict interior points of a source-rate bracket are evaluated.
        if difference == 0.0 {
            return Err(unresolved());
        }
        let (fraction, exponent) = quadratic_fraction(difference);
        sum.add((
            (arm.gamma.0 / fraction).copysign(difference),
            arm.gamma.1 - exponent,
        ));
    }
    Ok(sum.sign())
}

/// Already-admitted scalar sources; None selects the original dense solver.
///
/// Write omega²=K/m, a=C/(m*omega), g_j=gamma_j/omega², R_j=r_j/omega.
/// For a<=1 and sum(g_j)<=1/16, let the pair be -omega*d +/- i*omega*v,
/// and P=d²+v². The original rational characteristic equation gives
///
///   2*d = a + sum(g_j*R_j / (R_j² - 2*R_j*d + P)),
///   P = 1 / (1 - sum(g_j / (R_j² - 2*R_j*d + P))).
///
/// This map preserves 0<=d<=0.55, 1<=P<=1.1. Its infinity-norm derivative
/// is below 0.4: denominators are >=1-0.55² and R/denominator<=1/0.9.
/// Thus 64 bounded contractions recover the full pair without a loss floor.
///
/// The remaining poles are real. Between each successive distinct active
/// rate the rational equation changes sign, providing N real roots in
/// addition to the proven conjugate pair. Each root u_j lies above
/// max(r_(j-1), r_j/(1+4*sum(g)/3)): u²-a*omega*u+omega²>=3*omega²/4.
/// These short relative brackets retain arbitrarily separated source rates.
pub(super) fn weak_scalar_multi_maxwell(
    mass: f64,
    damping: f64,
    stiffness: f64,
    branches: &[MaxwellDampingBranch],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Option<Vec<C64>>, ModalError> {
    poll(cancelled)?;
    if stiffness <= 0.0 || damping < 0.0 || branches.len() < 2 {
        return Ok(None);
    }
    let (kf, ke) = quadratic_fraction(stiffness);
    let (mf, me) = quadratic_fraction(mass);
    let exponent = ke - me;
    let omega = quadratic_power_scale(
        det::sqrt(
            (kf / mf)
                * if exponent.rem_euclid(2) == 0 {
                    1.0
                } else {
                    2.0
                },
        ),
        exponent.div_euclid(2),
    );
    if !omega.is_finite() || omega <= 0.0 {
        return Err(unresolved());
    }
    let a = positive_ratio(&[damping], &[mass, omega]);
    if a > 1.0 {
        return Ok(None);
    }
    let mut arms = Vec::with_capacity(branches.len());
    let mut total = Dd::ZERO;
    for branch in branches {
        poll(cancelled)?;
        let strain = branch.strain[0].abs();
        let ratio = positive_ratio(&[branch.stiffness, strain, strain], &[stiffness]);
        if !ratio.is_finite() || ratio > MAX_STIFFNESS_RATIO {
            return Ok(None);
        }
        total = total + Dd::from_f64(ratio);
        if total.to_f64() > MAX_STIFFNESS_RATIO {
            return Ok(None);
        }
        let rate = 1.0 / branch.relaxation_time_s;
        if !rate.is_finite() || rate <= 0.0 {
            return Err(unresolved());
        }
        arms.push(Arm {
            rate,
            gamma: ratio_parts(&[branch.stiffness, strain, strain], &[mass]),
            ratio,
        });
    }
    let total = total.to_f64();
    let mut order: Vec<usize> = (0..arms.len())
        .filter(|&i| arms[i].gamma.0 != 0.0)
        .collect();
    order.sort_by(|&i, &j| arms[i].rate.total_cmp(&arms[j].rate).then(i.cmp(&j)));
    // Distinct input times whose reciprocals round to one rate cannot define
    // separate strict brackets here. Do not silently merge their source states.
    for pair in order.windows(2) {
        if arms[pair[0]].rate == arms[pair[1]].rate
            && branches[pair[0]].relaxation_time_s != branches[pair[1]].relaxation_time_s
        {
            return Ok(None);
        }
    }

    let denominator = |arm: &Arm, d: f64, p: f64| {
        if arm.rate >= omega {
            let ratio = omega / arm.rate;
            1.0 - 2.0 * d * ratio + p * ratio * ratio
        } else {
            let ratio = arm.rate / omega;
            ratio * (ratio - 2.0 * d) + p
        }
    };
    let mut d = 0.5 * a;
    let mut p = 1.0;
    let mut converged = false;
    for _ in 0..MAX_ITERATIONS {
        poll(cancelled)?;
        let mut loss = Dd::ZERO;
        let mut storage = Dd::ZERO;
        for (i, arm) in arms.iter().enumerate() {
            if i % 64 == 0 {
                poll(cancelled)?;
            }
            let den = denominator(arm, d, p);
            let ratio = if arm.rate >= omega {
                omega / arm.rate
            } else {
                arm.rate / omega
            };
            let weight = arm.ratio / den;
            loss = loss + Dd::from_f64(weight * ratio);
            storage = storage
                + Dd::from_f64(if arm.rate >= omega {
                    weight * ratio * ratio
                } else {
                    weight
                });
        }
        let next_d = 0.5 * (a + loss.to_f64());
        let next_p = 1.0 / (1.0 - storage.to_f64());
        if !next_d.is_finite()
            || !next_p.is_finite()
            || !(0.0..=0.55).contains(&next_d)
            || !(1.0..=1.1).contains(&next_p)
        {
            return Ok(None);
        }
        converged = (next_d - d).abs() <= 8.0 * f64::EPSILON * next_d
            && (next_p - p).abs() <= 8.0 * f64::EPSILON * next_p;
        (d, p) = (next_d, next_p);
        if converged {
            break;
        }
    }
    if !converged {
        return Ok(None);
    }

    // Recover PHYSICAL 2*decay from the original dimensional factors. Rounding
    // small g_j or d to zero above must not erase representable material heat.
    let viscous = ratio_parts(&[damping], &[mass]);
    let mut decay = ScaledSum::default();
    decay.add(viscous);
    for (i, (arm, branch)) in arms.iter().zip(branches).enumerate() {
        if i % 64 == 0 {
            poll(cancelled)?;
        }
        let strain = branch.strain[0].abs();
        let den = denominator(arm, d, p);
        decay.add(if arm.rate >= omega {
            ratio_parts(
                &[branch.stiffness, strain, strain, branch.relaxation_time_s],
                &[mass, den],
            )
        } else {
            ratio_parts(
                &[branch.stiffness, strain, strain, arm.rate],
                &[stiffness, den],
            )
        });
    }
    let decay = decay.physical(-1);
    let frequency = omega * det::sqrt(p - (decay / omega) * (decay / omega));
    if !decay.is_finite()
        || !frequency.is_finite()
        || frequency <= 0.0
        || ((damping != 0.0 || !order.is_empty()) && decay == 0.0)
    {
        return Err(unresolved());
    }
    let mut poles = Vec::with_capacity(branches.len() + 2);
    poles.push(C64::new(-decay, -frequency));
    poles.push(C64::new(-decay, frequency));
    for arm in &arms {
        if arm.gamma.0 == 0.0 {
            poles.push(C64::from_re(-arm.rate));
        }
    }

    let mut previous_rate: f64 = 0.0;
    let mut position = 0;
    while position < order.len() {
        poll(cancelled)?;
        let rate = arms[order[position]].rate;
        let end = position
            + order[position..]
                .iter()
                .take_while(|&&i| arms[i].rate == rate)
                .count();
        // Equal-time active arms contribute one coupled pole and independent
        // stress differences with the original exact relaxation rate.
        poles.extend((position + 1..end).map(|_| C64::from_re(-rate)));
        let mut lower = previous_rate.max(rate / (1.0 + 4.0 * total / 3.0));
        let mut upper = rate;
        if lower < upper
            && lower != previous_rate
            && real_sign(lower, viscous, stiffness, mass, &arms, cancelled)? < 0.0
        {
            return Ok(None);
        }
        for _ in 0..MAX_ITERATIONS {
            poll(cancelled)?;
            let middle = lower + 0.5 * (upper - lower);
            if middle == lower || middle == upper {
                break;
            }
            let sign = real_sign(middle, viscous, stiffness, mass, &arms, cancelled)?;
            if sign == 0.0 {
                lower = middle;
                upper = middle;
                break;
            } else if sign > 0.0 {
                lower = middle;
            } else {
                upper = middle;
            }
        }
        poles.push(C64::from_re(-(lower + 0.5 * (upper - lower))));
        previous_rate = rate;
        position = end;
    }
    poles.sort_by(|a, b| a.re.total_cmp(&b.re).then_with(|| a.im.total_cmp(&b.im)));
    poll(cancelled)?;
    Ok(Some(poles))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn branch(stiffness: f64, tau: f64, strain: f64) -> MaxwellDampingBranch {
        MaxwellDampingBranch {
            stiffness,
            relaxation_time_s: tau,
            strain: vec![strain],
        }
    }
    fn solve(m: f64, c: f64, k: f64, arms: &[MaxwellDampingBranch]) -> Vec<C64> {
        weak_scalar_multi_maxwell(m, c, k, arms, &mut || false)
            .unwrap()
            .unwrap()
    }
    fn close(actual: f64, expected: f64, tolerance: f64) {
        assert!(
            (actual / expected - 1.0).abs() <= tolerance,
            "{actual:e} != {expected:e}"
        );
    }

    #[test]
    fn g1_multi_arm_weak_decay_is_not_lost_beside_separated_relaxation_rates() {
        for stiffness in [1e-20, 1e-100, 1e-300] {
            let arms = [branch(stiffness, 0.5, 1.0), branch(stiffness, 2.0, 1.0)];
            let poles = solve(1.0, 0.0, 1.0, &arms);
            let expected = 0.4 * stiffness;
            assert_eq!(poles.len(), 4);
            for pole in poles.iter().filter(|s| s.im != 0.0) {
                close(pole.re, -expected, 8e-15);
                close(pole.im.abs(), 1.0, 8e-15);
            }
        }
        let arms = [
            branch(1e-20, 1e-100, 1.0),
            branch(1e-20, 1.0, 1.0),
            branch(1e-20, 1e100, 1.0),
        ];
        let poles = solve(1.0, 0.0, 1.0, &arms);
        assert_eq!(poles.len(), 5);
        for pole in poles.iter().filter(|s| s.im != 0.0) {
            close(pole.re, -2.5e-21, 8e-15);
        }
        let real: Vec<_> = poles.iter().filter(|s| s.im == 0.0).map(|s| s.re).collect();
        for (actual, expected) in real.iter().zip([-1e100, -1.0, -1e-100]) {
            close(*actual, expected, 8e-15);
        }
    }

    #[test]
    fn g1_complete_quartic_matches_manufactured_finite_loss_roots() {
        let (alpha, beta, r1, r2, u1, u2) = (0.1_f64, 4.0, 1.0, 4.0, 0.995, 3.99);
        let d = 0.5 * (alpha + r1 + r2 - u1 - u2);
        let product = beta * r1 * r2 / (u1 * u2);
        let total = u1 * u2 + 2.0 * d * (u1 + u2) + product - beta - alpha * (r1 + r2) - r1 * r2;
        let weighted = 2.0 * d * u1 * u2 + product * (u1 + u2) - alpha * r1 * r2 - beta * (r1 + r2);
        let g1 = (weighted - total * r1) / (r2 - r1);
        let g2 = total - g1;
        let arms = [branch(g1, 1.0 / r1, 1.0), branch(g2, 1.0 / r2, 1.0)];
        let poles = solve(1.0, alpha, beta, &arms);
        let v = (product - d * d).sqrt();
        let expected = [
            C64::from_re(-u2),
            C64::from_re(-u1),
            C64::new(-d, -v),
            C64::new(-d, v),
        ];
        for (&actual, expected) in poles.iter().zip(expected) {
            assert!(
                (actual - expected).abs() <= 2e-14 * expected.abs(),
                "{actual:?} != {expected:?}"
            );
        }
        assert_eq!(poles, solve(1.0, alpha, beta, &arms));
    }

    #[test]
    fn g1_all_real_and_complex_poles_satisfy_the_original_multi_arm_polynomial() {
        let arms = [
            branch(0.01, 0.2, 1.0),
            branch(0.02, 0.7, -1.0),
            branch(0.03, 2.0, 1.0),
            branch(0.015, 4.0, 1.0),
        ];
        let poles = solve(2.0, 0.2, 8.0, &arms);
        assert_eq!(poles.len(), 6);
        for pole in poles {
            let mut product = C64::ONE;
            for arm in &arms {
                product = product * (pole + C64::from_re(1.0 / arm.relaxation_time_s));
            }
            let mechanical = (pole * pole).scale(2.0) + pole.scale(0.2) + C64::from_re(8.0);
            let mut value = mechanical * product;
            let mut scale = value.abs();
            for (i, arm) in arms.iter().enumerate() {
                let mut term = pole.scale(arm.stiffness);
                for (j, other) in arms.iter().enumerate() {
                    if i != j {
                        term = term * (pole + C64::from_re(1.0 / other.relaxation_time_s));
                    }
                }
                value = value + term;
                scale += term.abs();
            }
            assert!(value.abs() <= 2e-11 * scale, "{pole:?}: {value:?}/{scale}");
        }
    }

    #[test]
    fn g0_repeated_and_uncoupled_arms_retain_every_material_state() {
        let arms = [
            branch(0.02, 0.5, 1.0),
            branch(0.03, 0.5, -1.0),
            branch(2.0, 0.5, 0.0),
            branch(0.01, 3.0, 1.0),
        ];
        let poles = solve(1.0, 0.0, 4.0, &arms);
        assert_eq!(poles.len(), 6);
        assert_eq!(
            poles.iter().filter(|&&s| s == C64::from_re(-2.0)).count(),
            2
        );
        let combined = solve(
            1.0,
            0.0,
            4.0,
            &[branch(0.05, 0.5, 1.0), branch(0.01, 3.0, 1.0)],
        );
        let active: Vec<_> = poles.iter().filter(|&&s| s != C64::from_re(-2.0)).collect();
        for (&actual, expected) in active.iter().zip(combined) {
            assert!((*actual - expected).abs() < 2e-14 * expected.abs());
        }
    }

    #[test]
    fn g3_source_products_survive_units_and_dimensionless_underflow() {
        let reference = solve(
            1.0,
            0.0,
            4.0,
            &[branch(1e-20, 0.5, 1.0), branch(2e-20, 2.0, 1.0)],
        );
        for exponent in [-800, 800] {
            let mass = 2.0_f64.powi(exponent);
            let inverse = 2.0_f64.powi(-exponent);
            let poles = solve(
                mass,
                0.0,
                4.0 * mass,
                &[
                    branch(1e-20 * inverse, 0.5, mass),
                    branch(2e-20 * inverse, 2.0, mass),
                ],
            );
            for (&actual, expected) in poles.iter().zip(&reference) {
                close(actual.re, expected.re, 8e-15);
                assert_eq!(actual.im, expected.im);
            }
        }
        let poles = solve(
            1.0,
            0.0,
            1e200,
            &[branch(1e-200, 1e-100, 1.0), branch(1e-200, 2e-100, 1.0)],
        );
        for pole in poles.iter().filter(|s| s.im != 0.0) {
            close(pole.re, -4.5e-301, 8e-15);
        }
    }

    #[test]
    fn g0_admission_boundaries_and_cancellation_publish_no_partial_spectrum() {
        let arms = [branch(0.03125, 1.0, 1.0), branch(0.03125, 2.0, 1.0)];
        assert!(
            weak_scalar_multi_maxwell(1.0, 1.0, 1.0, &arms, &mut || false)
                .unwrap()
                .is_some()
        );
        assert!(
            weak_scalar_multi_maxwell(1.0, 1.01, 1.0, &arms, &mut || false)
                .unwrap()
                .is_none()
        );
        let excessive = [branch(0.03125, 1.0, 1.0), branch(0.03126, 2.0, 1.0)];
        assert!(
            weak_scalar_multi_maxwell(1.0, 0.0, 1.0, &excessive, &mut || false)
                .unwrap()
                .is_none()
        );
        for stop in [0, 4, 12] {
            let mut checks = 0;
            assert_eq!(
                weak_scalar_multi_maxwell(1.0, 0.0, 1.0, &arms, &mut || {
                    checks += 1;
                    checks > stop
                }),
                Err(ModalError::Cancelled)
            );
        }
    }
}
