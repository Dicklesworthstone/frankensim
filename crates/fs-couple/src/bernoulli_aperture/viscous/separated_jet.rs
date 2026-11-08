//! The separated-jet entrance law of a reed channel: van Zon, Hirschberg, Gilbert and Wijnands 1990, "Flow through the reed
//! channel of a single reed music instrument", J. Phys. Colloques 51 (C2) 821-824 (HAL jpa-00230506), eqs. 1-8 and Fig. 2.
//! A 2-D clarinet-mouthpiece model: height `h`, length `L`, width `w`, volume flow `phi`, mean speed `v = phi / (w h)`,
//! `Re = phi / (nu w) = v h / nu`, Bernoulli flow `phi_B = h w sqrt(2 dp / rho)` (eq. 2), `S` the length of the separated
//! region at the entrance.
//!
//! * `L / h < 3` (eq. 3): the jet never reattaches, `phi / phi_B = alpha` with the measured `alpha = 0.45 +- 0.05`, so
//!   `dp = rho v |v| / (2 alpha^2)`.
//! * `L / h >= 3`: the jet reattaches at `S` (the paper's fit is `S = 2 h`, within 3 % of the 2-D measurements for
//!   `Re <= 4000`, laminar). Boundary layers then grow (von Karman integral, linear profile, eq. 4)
//!   `(x - S) / h = (Re / 6) f(d)`, `f(d) = 4 d + 9 ln(1 - d) + 5 d / (1 - d)`, `d = delta / h`, and the core is a Bernoulli flow
//!   (eq. 5): `phi / phi_B = 1 - d(L)`, that is `dp = rho v^2 / (2 (1 - d)^2)`, with `d(L)` from `f(d) = 6 (L - S) / (h Re)`
//!   (an increasing function of `d`, solved here by Newton from the right). This holds while `d(L) <= d_c`,
//!   `d_c = (4/9) [1 - sqrt(5/32)]`, where the boundary-layer flow has the volume and momentum flux of a Poiseuille flow.
//! * Past that length `l_c` (`l_c - S = c Re h`, eq. 6, `c = f(d_c) / 6`) the flow is Poiseuille (eq. 1 per length):
//!   `dp = rho v^2 / (2 (1 - d_c)^2) + 12 mu (L - l_c) v / h^2`, which is `dp = TAIL rho v |v| + 12 mu (L - S) v / h^2` with
//!   `TAIL = 1 / (2 (1 - d_c)^2) - 12 c`. Continuous with the entrance branch at `d = d_c` (value, not slope).
//!
//! Print checks (HAL copy, md5 cf5ea392f21c41cc4ddbfb776f180b8f), settled with MB 12740:
//!
//! 1. Eq. 4 prints its last term as `5 delta / (h - delta)`: the braces are positive and expand as `d^2 / 2 + 2 d^3`. The
//!    `5 b / (b - h)` of Chatziioannou and van Walstijn 2008 (their eq. 3) is the misprint.
//! 2. Eq. 8 as printed lacks the 9 before the logarithm. `c` is taken from eq. 4 at `d_c` (`c = f(d_c) / 6 = 0.0159368`), which is
//!    what eq. 6 needs by its own definition of `l_c`.
//! 3. Eq. 7 as printed is NOT used, because it is inconsistent with eq. 5 at `l_c`: at `l_c = L` it gives `phi / phi_B = 0.903` against
//!    eq. 5's `0.731` (a 24 % step), and at `h Re / [6 (L - S)] = 10` it gives 0.896 against the 0.727 read off the solid curve of
//!    Fig. 2. The tail above is a CONSTRUCTION FROM EQS. 1, 5 AND 6 (eq. 5 for the core up to `l_c`, eq. 1 per length beyond it, eq. 6
//!    for `l_c`): continuous at `d_c` in value, and it reproduces the Fig. 2 solid curve at the fifteen read points of the test
//!    `the_law_reproduces_the_printed_theory_curve_of_fig_2` to 0.011 (tolerance 0.015, the read error of the print). No coefficient
//!    is fitted or borrowed.
//!
//! Limits kept in the open: the law is stationary and laminar (no unsteady channel inertia, no turbulence above `Re` about 4000),
//! 2-D (the printed fit is for a 3 cm square pipe with `h` = 1 mm and `L` = 4 or 8 mm, not shown to be this mouthpiece),
//! and the branch at `L / h = 3` is a step: the paper gives the two data sets (`L / h <= 3` and `>= 4`) and no blend. A massive reed
//! whose opening dwells near `h = L / 3` chatters across the step: a 1 mm channel on a 0.4 mm reed (`L / h = 2.5` at rest) made its
//! implicit contact solve report that it did not resolve the momentum balance (the reed tests' fixture channel is 2.5 mm), so a channel
//! with `L / h` at rest near or below 3 is outside what this law has been run on.

use fs_math::det;

/// `delta_c / h = (4/9) (1 - sqrt(5/32))` (printed under eq. 4).
pub(super) const DELTA_C: f64 = 0.268_762_352_212_867_81;
/// `f(delta_c)`, `f(d) = 4 d + 9 ln(1 - d) + 5 d / (1 - d)` (the braces of eq. 4).
pub(super) const F_C: f64 = 0.095_620_796_239_128_49;
/// `c = f(delta_c) / 6` (eq. 6 and eq. 8 with the printed 9 restored).
pub(super) const C: f64 = 0.015_936_799_373_188_08;
/// `1 / (2 (1 - delta_c)^2) - 12 c`: the coefficient of `rho v |v|` once the flow is Poiseuille.
pub(super) const TAIL: f64 = 0.743_847_343_454_391_3;
/// `alpha` of eq. 3 (`L / h < 3`): the measured `0.45 +- 0.05`.
pub(super) const CONTRACTION: f64 = 0.45;
/// `L / h` below which the jet does not reattach (eq. 3).
pub(super) const NO_REATTACHMENT: f64 = 3.0;
/// `S / h` (the paper's fit, within 3 % for `Re <= 4000`).
pub(super) const SEPARATION: f64 = 2.0;

/// Pressure drop `dp` [Pa] at speed `s = |v| >= 0`, with `d dp / d s` [Pa s / m] and `d dp / d h` at fixed speed [Pa / m].
#[derive(Debug, Clone, Copy, PartialEq)]
pub(super) struct Drop {
    pub dp: f64,
    pub d_dv: f64,
    pub d_dh: f64,
}

/// Which line of the paper's law a geometry and speed are on.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(super) enum Branch {
    /// `L / h < 3`: no reattachment (eq. 3).
    Contracted,
    /// `L / h >= 3`, boundary layers still growing at `L` (eqs. 4, 5).
    Entrance,
    /// `L / h >= 3`, Poiseuille past `l_c` (eqs. 1, 6).
    Developed,
}

/// `f(d)` of eq. 4. The closed form cancels its first-order terms (`4 - 9 + 5 = 0`), so a small `d` uses the series
/// `sum_{k >= 2} (5 k - 9) / k d^k`, whose terms are all positive.
pub(super) fn braces(d: f64) -> f64 {
    if d < 0.1 {
        let mut power = d * d;
        let mut sum = 0.0;
        for k in 2..=26 {
            sum += ((5.0 * f64::from(k) - 9.0) / f64::from(k)) * power;
            power *= d;
        }
        sum
    } else {
        (4.0 * d + 9.0 * det::ln(1.0 - d)) + 5.0 * d / (1.0 - d)
    }
}

/// `f'(d) = d (1 + 4 d) / (1 - d)^2`.
fn braces_slope(d: f64) -> f64 {
    let one_minus = 1.0 - d;
    (d * (1.0 + 4.0 * d)) / (one_minus * one_minus)
}

/// The `d` in `[0, d_c]` with `f(d) = y` for `y` in `[0, f(d_c)]`. `f` is increasing and convex and `f(d) >= d^2 / 2`, so Newton
/// from `min(d_c, sqrt(2 y))` stays at or right of the root and converges monotonically.
pub(super) fn thickness(y: f64) -> f64 {
    if y <= 0.0 {
        return 0.0;
    }
    if y >= F_C {
        return DELTA_C;
    }
    let mut d = DELTA_C.min(det::sqrt(2.0 * y));
    for _ in 0..64 {
        let step = (braces(d) - y) / braces_slope(d);
        d -= step;
        if step.abs() <= 4.0 * f64::EPSILON * d {
            break;
        }
    }
    d
}

/// The branch at speed `s > 0`. `length` and `h` positive, `rho` and `mu` the gas density and dynamic viscosity.
pub(super) fn branch(s: f64, h: f64, length: f64, rho: f64, mu: f64) -> Branch {
    if length < NO_REATTACHMENT * h {
        return Branch::Contracted;
    }
    let reynolds = (s * h) * (rho / mu);
    let x = length / h - SEPARATION;
    if C * reynolds >= x { Branch::Entrance } else { Branch::Developed }
}

/// Speed at which a reattached geometry passes from the Poiseuille branch to the entrance branch (`l_c = L`), `X nu / (c h)`.
pub(super) fn boundary_speed(h: f64, length: f64, rho: f64, mu: f64) -> f64 {
    let x = length / h - SEPARATION;
    (x * (mu / rho)) / (C * h)
}

/// The law at speed `s >= 0` (zero returns the rest slope: Poiseuille, or none on the contracted branch).
pub(super) fn drop_at(s: f64, h: f64, length: f64, rho: f64, mu: f64) -> Drop {
    let viscous_slope = 12.0 * mu * (length - SEPARATION * h) / (h * h);
    if length < NO_REATTACHMENT * h {
        let k = rho / (2.0 * CONTRACTION * CONTRACTION);
        return Drop { dp: (k * s) * s, d_dv: 2.0 * k * s, d_dh: 0.0 };
    }
    if s == 0.0 {
        return Drop { dp: 0.0, d_dv: viscous_slope, d_dh: 0.0 };
    }
    let nu = mu / rho;
    // d/dh of (L - 2 h) / h^2 = -2 (L - h) / h^3.
    let slope_h = -2.0 * (length - h) / (h * h * h);
    match branch(s, h, length, rho, mu) {
        Branch::Developed => Drop {
            dp: (TAIL * rho * s) * s + (viscous_slope * s),
            d_dv: 2.0 * TAIL * rho * s + viscous_slope,
            d_dh: 12.0 * mu * s * slope_h,
        },
        Branch::Entrance => {
            let reynolds = s * h / nu;
            let x = length / h - SEPARATION;
            let y = 6.0 * x / reynolds;
            let d = thickness(y);
            let u = 1.0 / (1.0 - d);
            let slope = braces_slope(d);
            // y = 6 nu (L - 2 h) / (s h^2): d y / d s = -y / s, d y / d h = 6 nu / s * (-2 (L - h) / h^3).
            let dy_dh = 6.0 * nu / s * slope_h;
            let half = 0.5 * rho * s;
            Drop {
                dp: (half * s) * (u * u),
                d_dv: (rho * s * u * u) * (1.0 - u * y / slope),
                d_dh: (2.0 * half * s * u * u * u) * dy_dh / slope,
            }
        }
        Branch::Contracted => unreachable!("handled above"),
    }
}

/// Positive root of `a s^2 + b s = drive` (`a > 0`, `b >= 0`, `drive >= 0`), rationalized and scaled like the fully developed
/// law's own root so that `b` is never squared against `drive`.
fn quadratic_speed(a: f64, b: f64, drive: f64) -> Option<f64> {
    if drive == 0.0 {
        return Some(0.0);
    }
    let q = 2.0 * det::sqrt(a) * det::sqrt(drive);
    let scale = b.max(q);
    if !(scale.is_finite() && scale > 0.0) {
        return None;
    }
    let (bn, qn) = (b / scale, q / scale);
    let speed = (2.0 * (drive / scale)) / (bn + det::sqrt(bn * bn + qn * qn));
    speed.is_finite().then_some(speed)
}

/// Speed `s >= 0` solving `dp(s) + wave_load s - recovery s^2 = drive` (`drive >= 0`): the slit law in series with the bore's
/// wave load `Z A` and, for forward flow through a double reed's staple, the pressure recovery `recovery = C_P rho / 2 (A / S_c)^2`
/// (`recovery < rho / 2`, so `g` stays increasing). `None` only if the arithmetic leaves the finite set or does not converge.
pub(super) fn speed_for_drive(
    h: f64, length: f64, rho: f64, mu: f64, wave_load: f64, recovery: f64, drive: f64,
) -> Option<f64> {
    if drive == 0.0 {
        return Some(0.0);
    }
    if length < NO_REATTACHMENT * h {
        let a = rho / (2.0 * CONTRACTION * CONTRACTION) - recovery;
        return quadratic_speed(a, wave_load, drive);
    }
    let developed_a = TAIL * rho - recovery;
    let developed_b = wave_load + 12.0 * mu * (length - SEPARATION * h) / (h * h);
    let boundary = boundary_speed(h, length, rho, mu);
    let at_boundary = (developed_a * boundary) * boundary + developed_b * boundary;
    if drive <= at_boundary {
        return quadratic_speed(developed_a, developed_b, drive);
    }
    // Entrance branch: dp(s) >= rho s^2 / 2, so the root lies below the root of the bound.
    let mut lo = boundary;
    let mut hi = quadratic_speed(0.5 * rho - recovery, wave_load, drive)?.max(boundary);
    let mut s = hi;
    for _ in 0..200 {
        let law = drop_at(s, h, length, rho, mu);
        let g = ((law.dp - (recovery * s) * s) + wave_load * s) - drive;
        let slope = (law.d_dv - 2.0 * recovery * s) + wave_load;
        if !g.is_finite() || !(slope > 0.0) {
            return None;
        }
        if g.abs() <= 8.0 * f64::EPSILON * drive {
            return Some(s);
        }
        if g > 0.0 { hi = s; } else { lo = s; }
        let newton = s - g / slope;
        s = if newton > lo && newton < hi { newton } else { f64::midpoint(lo, hi) };
        if s == lo || s == hi {
            return Some(s);
        }
    }
    None
}

#[cfg(test)]
#[allow(clippy::float_cmp)]
mod tests {
    use super::*;

    const RHO: f64 = 1.2;
    const MU: f64 = 1.81e-5;

    fn close(actual: f64, expected: f64, relative: f64) {
        assert!((actual - expected).abs() <= relative * expected.abs().max(1e-250), "actual={actual:e}, expected={expected:e}");
    }

    /// The universal ratio `phi / phi_B` at `x = h Re / (6 (L - S))` for `L / h = ratio`, through the pressure law.
    fn flow_ratio(x: f64, h: f64, ratio: f64) -> f64 {
        let length = ratio * h;
        let reynolds = 6.0 * (ratio - SEPARATION) * x;
        let s = reynolds * (MU / RHO) / h;
        let law = drop_at(s, h, length, RHO, MU);
        s / det::sqrt(2.0 * law.dp / RHO)
    }

    #[test]
    fn the_constants_are_the_printed_closed_forms() {
        let d_c = (4.0 / 9.0) * (1.0 - det::sqrt(5.0 / 32.0));
        close(DELTA_C, d_c, 1e-15);
        close(F_C, 4.0 * d_c + 9.0 * det::ln(1.0 - d_c) + 5.0 * d_c / (1.0 - d_c), 1e-13);
        close(C, F_C / 6.0, 1e-15);
        close(TAIL, 1.0 / (2.0 * (1.0 - d_c) * (1.0 - d_c)) - 12.0 * C, 1e-14);
        // The two conditions that define delta_c: a linear boundary layer has the flux of a Poiseuille flow of the same volume
        // flow (volume flow U (h - d) = phi / w, momentum U^2 (h - 4 d / 3) = (6/5) (phi / w)^2 / h).
        close(1.0 - 4.0 * d_c / 3.0, 1.2 * (1.0 - d_c) * (1.0 - d_c), 1e-14);
        assert_eq!(braces(0.0), 0.0);
        assert!(braces(0.05) > 0.0, "the braces of eq. 4 are positive: the printed 5 d / (h - d), not the 5 b / (b - h) misprint");
        close(braces(0.01), 0.5e-4 + 2e-6 + 2.75e-8 + 3.2e-10, 2e-3); // d^2/2 + 2 d^3 + 11 d^4 / 4 + 16 d^5 / 5
    }

    #[test]
    fn the_series_and_the_closed_form_of_the_braces_agree_where_they_meet() {
        for d in [0.05, 0.0999999, 0.1000001, 0.15, DELTA_C] {
            let closed = (4.0 * d + 9.0 * det::ln(1.0 - d)) + 5.0 * d / (1.0 - d);
            let mut power = d * d;
            let mut series = 0.0;
            for k in 2..=200 {
                series += ((5.0 * f64::from(k) - 9.0) / f64::from(k)) * power;
                power *= d;
            }
            // The closed form is used from 0.1 up; below it only this comparison sees it, with the cancellation of its first-order terms
            // (about 3e-13 of the value at d = 0.05 with the software logarithm).
            close(closed, series, 5e-12);
            close(braces(d), series, 3e-13);
        }
    }

    #[test]
    fn the_thickness_inverts_the_braces_over_twelve_decades() {
        for k in 0..=48 {
            let y = F_C * 10f64.powf(-0.25 * f64::from(k));
            let d = thickness(y);
            assert!((0.0..=DELTA_C).contains(&d));
            close(braces(d), y, 1e-12);
        }
        assert_eq!(thickness(F_C), DELTA_C);
        assert_eq!(thickness(1.0), DELTA_C);
        assert_eq!(thickness(0.0), 0.0);
    }

    #[test]
    fn the_poiseuille_limit_is_eq_1_over_the_reattached_length() {
        // Long channel, low Re: dp -> 12 mu (L - S) v / h^2, and over L within 2 h / L.
        for (ratio, tolerance) in [(50.0, 0.04), (1000.0, 0.002)] {
            let (h, v) = (1e-4, 1e-3 * (MU / RHO) / 1e-4); // Re = 1e-3
            let length = ratio * h;
            let law = drop_at(v, h, length, RHO, MU);
            close(law.dp, 12.0 * MU * (length - 2.0 * h) * v / (h * h), 1e-5); // the quadratic term is 1.3e-6 of it at Re = 1e-3
            close(law.dp, 12.0 * MU * length * v / (h * h), tolerance + 1e-6);
            close(law.d_dv, 12.0 * MU * (length - 2.0 * h) / (h * h), 1e-5);
            assert_eq!(branch(v, h, length, RHO, MU), Branch::Developed);
        }
        let rest = drop_at(0.0, 1e-4, 5e-3, RHO, MU);
        assert_eq!(rest.dp, 0.0);
        close(rest.d_dv, 12.0 * MU * (5e-3 - 2e-4) / 1e-8, 1e-14);
    }

    #[test]
    fn below_l_over_h_3_the_flow_is_alpha_times_the_bernoulli_flow() {
        let h = 4e-4;
        for ratio in [0.5, 1.0, 2.0, 2.999] {
            for v in [0.01, 1.0, 30.0] {
                let law = drop_at(v, h, ratio * h, RHO, MU);
                let phi_ratio = v / det::sqrt(2.0 * law.dp / RHO);
                close(phi_ratio, CONTRACTION, 1e-14);
                assert_eq!(branch(v, h, ratio * h, RHO, MU), Branch::Contracted);
                // independent of the viscosity and of the length
                assert_eq!(law.dp, drop_at(v, h, 0.5 * ratio * h, RHO, 7.0 * MU).dp);
            }
        }
        assert_ne!(branch(30.0, h, 3.0 * h, RHO, MU), Branch::Contracted, "L / h = 3 itself is reattached");
    }

    #[test]
    fn the_two_reattached_branches_meet_at_delta_c_in_value() {
        for (ratio, h) in [(3.5, 4e-4), (6.0, 1e-4), (20.0, 5e-5), (350.0, 4e-5)] {
            let length = ratio * h;
            let boundary = boundary_speed(h, length, RHO, MU);
            let eps = 1e-9;
            let below = drop_at(boundary * (1.0 - eps), h, length, RHO, MU);
            let above = drop_at(boundary * (1.0 + eps), h, length, RHO, MU);
            assert_eq!(branch(boundary * (1.0 - eps), h, length, RHO, MU), Branch::Developed);
            assert_eq!(branch(boundary * (1.0 + eps), h, length, RHO, MU), Branch::Entrance);
            close(below.dp, above.dp, 5e-9);
            // the value there is the Bernoulli head of the core, rho v^2 / (2 (1 - d_c)^2), with eq. 5: phi / phi_B = 1 - d_c
            close(above.dp, 0.5 * RHO * boundary * boundary / ((1.0 - DELTA_C) * (1.0 - DELTA_C)), 5e-9);
            // the slope has a kink (2.6 % in the coefficient of rho v): the paper's two regimes are matched in value only
            let kink = below.d_dv / (RHO * boundary) - above.d_dv / (RHO * boundary);
            assert!(kink > 0.02 && kink < 0.06, "{kink}");
        }
    }

    // Fig. 2, solid curve ("theory, S = 2 h"): phi / phi_B against x = h Re / [6 (L - S)], read off the 600 dpi print, axes
    // fitted to the left-edge ticks (1.0 at pixel 118, 0.1 per 72.25 px). The curve is one universal function of x on both
    // branches. Read tolerance 0.015 (about 11 px; the largest difference found is 0.011).
    const FIG_2: [(f64, f64); 15] = [
        (1.1, 0.432), (1.5, 0.485), (2.0, 0.530), (3.0, 0.588), (5.0, 0.650), (7.0, 0.691), (10.0, 0.727), (15.0, 0.763),
        (20.0, 0.788), (30.0, 0.811), (100.0, 0.879), (200.0, 0.907), (300.0, 0.921), (500.0, 0.938), (900.0, 0.954),
    ];

    #[test]
    fn the_law_reproduces_the_printed_theory_curve_of_fig_2() {
        let mut worst: f64 = 0.0;
        for (x, read) in FIG_2 {
            for (h, ratio) in [(1e-3, 6.0), (1e-4, 12.0), (3e-5, 40.0), (1e-4, 300.0)] {
                let ours = flow_ratio(x, h, ratio);
                worst = worst.max((ours - read).abs());
                assert!((ours - read).abs() <= 0.015, "x={x} L/h={ratio}: {ours} vs read {read}");
                // universal in x: the same value for every geometry
                close(ours, flow_ratio(x, 1e-4, 12.0), 1e-9);
            }
        }
        assert!(worst > 0.0 && worst < 0.012, "{worst}");
        // the Poiseuille end of the curve: phi / phi_B -> sqrt(x) / 2 at small x, 1 - d -> 1 at large x
        close(flow_ratio(1e-4, 1e-4, 12.0), 0.5 * det::sqrt(1e-4), 2e-3);
        assert!(flow_ratio(1e8, 1e-4, 12.0) > 0.999);
        close(flow_ratio(1.0 / F_C, 1e-4, 12.0), 1.0 - DELTA_C, 1e-8);
        // On the Poiseuille branch the curve is the closed form (phi / phi_B)^2 = x / (4 + 2 TAIL x); at x = 1 / f(d_c) it is (1 - d_c)^2.
        for x in [0.01, 0.3, 1.0, 4.0, 10.0] {
            close(flow_ratio(x, 1e-4, 12.0), det::sqrt(x / (4.0 + 2.0 * TAIL * x)), 1e-9);
        }
        close(det::sqrt((1.0 / F_C) / (4.0 + 2.0 * TAIL / F_C)), 1.0 - DELTA_C, 1e-13);
    }

    #[test]
    fn the_printed_eq_7_is_not_the_theory_of_fig_2() {
        // eq. 7 as printed with c from eq. 4: the flow ratio it gives at x = 10 and at the joint.
        let eq7 = |x: f64| {
            let ell = 6.0 * C * x;
            let q = 12.0 * C * (1.0 - DELTA_C) / (24.0 * C - 1.0);
            let r = (24.0 * C - 1.0) / (72.0 * (1.0 - DELTA_C));
            let delta = (2.0 * ell / q - ell * ell / (q * q)) / r;
            (ell / C) / det::sqrt(2.0 * delta)
        };
        assert!((eq7(10.0) - 0.727).abs() > 0.15, "{}", eq7(10.0));
        assert!((eq7(1.0 / F_C) - (1.0 - DELTA_C)).abs() > 0.15, "{}", eq7(1.0 / F_C));
    }

    #[test]
    fn the_drop_is_odd_increasing_and_its_tangents_match_differences() {
        for (h, length) in [(4e-5, 1.4e-3), (1e-4, 1.4e-3), (3e-4, 1.4e-3), (4e-4, 1.4e-3), (1e-4, 25e-3), (3e-4, 25e-3), (1e-4, 2e-4)] {
            let mut previous = 0.0;
            for k in 1..=400 {
                let s = 0.25 * f64::from(k);
                let law = drop_at(s, h, length, RHO, MU);
                assert!(law.dp > previous, "h={h} L={length} s={s}");
                assert!(law.d_dv > 0.0);
                previous = law.dp;
                if k % 40 == 0 {
                    let (es, eh) = (s * 1e-6, h * 1e-6);
                    let on_branch = branch(s * (1.0 + 1e-6), h, length, RHO, MU) == branch(s * (1.0 - 1e-6), h, length, RHO, MU)
                        && branch(s, h * (1.0 + 1e-6), length, RHO, MU) == branch(s, h * (1.0 - 1e-6), length, RHO, MU);
                    if !on_branch { continue; }
                    let fs = (drop_at(s + es, h, length, RHO, MU).dp - drop_at(s - es, h, length, RHO, MU).dp) / (2.0 * es);
                    let fh = (drop_at(s, h + eh, length, RHO, MU).dp - drop_at(s, h - eh, length, RHO, MU).dp) / (2.0 * eh);
                    close(law.d_dv, fs, 2e-6);
                    close(law.d_dh, fh, 2e-5);
                }
            }
        }
    }

    #[test]
    fn the_speed_solver_inverts_the_law_with_a_wave_load_and_a_recovery_on_every_branch() {
        let mut entrance = 0;
        let mut developed = 0;
        let mut contracted = 0;
        for (h, length) in [(4e-5, 1.4e-3), (1e-4, 1.4e-3), (3e-4, 1.4e-3), (4e-4, 1.4e-3), (6e-4, 1.4e-3), (1e-4, 25e-3), (3e-4, 25e-3)] {
            for wave_load in [0.0, 3e3, 4e6] {
                for recovery in [0.0, 0.3 * RHO / 2.0, 0.9 * RHO / 2.0] {
                    for s in [1e-9, 1e-3, 0.3, 4.0, 25.0, 90.0, 250.0] {
                        let law = drop_at(s, h, length, RHO, MU);
                        let drive = ((law.dp - recovery * s * s) + wave_load * s).max(0.0);
                        let solved = speed_for_drive(h, length, RHO, MU, wave_load, recovery, drive).unwrap();
                        close(solved, s, 1e-10);
                        match branch(s, h, length, RHO, MU) {
                            Branch::Entrance => entrance += 1,
                            Branch::Developed => developed += 1,
                            Branch::Contracted => contracted += 1,
                        }
                    }
                }
            }
        }
        assert!(entrance > 50 && developed > 100 && contracted > 30, "{entrance} {developed} {contracted}");
        assert_eq!(speed_for_drive(1e-4, 1.4e-3, RHO, MU, 1e3, 0.0, 0.0), Some(0.0));
    }
}
