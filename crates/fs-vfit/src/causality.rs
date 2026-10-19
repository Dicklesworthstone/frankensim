//! Causality of sampled responses by generalized dispersion relations with
//! subtractions (Triverio & Grivet-Talocia; IRP M2.5-W3b).
//!
//! On this crate's `s = i w` axis a causal real impulse response has
//! `H(-w) = conj H(w)`, and `H` is analytic for `Im w < 0`. With subtraction
//! points `w_1..w_n`, `Z(w) = PROD (w - w_k)` and the Lagrange basis `L_k`
//! on those points:
//!
//! `H(w) = SUM_k H(w_k) L_k(w) + (i/pi) Z(w) PV INT H(x) / (Z(x) (x - w)) dx`.
//!
//! The integral is truncated to the sampled band `[-W, W]` by Hermitian
//! extension, and its poles at `w` and the `w_k` are integrated analytically.
//! A sample is a violation when the reconstruction misses it by more than the
//! sum of three declared terms:
//! - truncation: the out-of-band tail, assuming
//!   `|H| <= tail_bound_ratio * max |H|` outside the band;
//! - discretization: twice the leading trapezoid error of the regular
//!   integrand, estimated from its second differences;
//! - data: `relative_data_tolerance * max |H|` per sample, propagated through
//!   the reconstruction by the triangle inequality.
//!
//! Points above `0.9 W` are not tested, because the truncation term grows
//! without bound at the band edge. A passing response is consistent with
//! causality under these assumptions; it is not a proof.

use fs_math::{c64::C64, det};

/// Declared assumptions of a causality test.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CausalityOptions {
    /// Subtraction points (1 to 8), spread evenly over the band.
    pub subtractions: usize,
    /// Out-of-band bound `|H| <= ratio * max_in_band |H|`, at least 1.
    pub tail_bound_ratio: f64,
    /// Per-sample data error bound relative to `max_in_band |H|`.
    pub relative_data_tolerance: f64,
}

/// Measured outcome of a passing causality test.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CausalityReport {
    /// Samples tested; zero for an identically zero response.
    pub tested: usize,
    /// Largest reconstruction error over its bound (below 1 to pass).
    pub worst_ratio: f64,
    /// Frequency of `worst_ratio` [rad/s].
    pub worst_omega: f64,
    /// Largest reconstruction error over `max |H|`.
    pub max_relative_error: f64,
}

/// Causality-test refusals.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum CausalityError {
    /// The grid or data cannot be tested.
    BadSample {
        /// What was wrong.
        what: &'static str,
    },
    /// Options outside their documented ranges.
    BadOptions,
    /// The reconstruction misses a sample by more than its bound.
    NonCausal {
        /// Row-major matrix entry (0 for a scalar response).
        entry: usize,
        /// Frequency of the worst violation [rad/s].
        omega: f64,
        /// Reconstruction error there.
        error: f64,
        /// Its bound.
        bound: f64,
    },
}

impl core::fmt::Display for CausalityError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::BadSample { what } => write!(f, "causality test refused: {what}"),
            Self::BadOptions => write!(f, "causality options out of range"),
            Self::NonCausal {
                entry,
                omega,
                error,
                bound,
            } => write!(
                f,
                "entry {entry} is not causal at {omega} rad/s: error {error:.3e} > bound {bound:.3e}"
            ),
        }
    }
}

impl std::error::Error for CausalityError {}

/// Test one sampled response `data[j] = H(i omega[j])`. The grid must be
/// uniform and start at 0 (a real DC sample) or at half a step.
///
/// # Errors
/// [`CausalityError::NonCausal`] at the worst violation, or a grid, data or
/// options refusal.
pub fn check_causality(
    omega: &[f64],
    data: &[C64],
    opts: &CausalityOptions,
) -> Result<CausalityReport, CausalityError> {
    check_entry(omega, data, opts, 0)
}

/// Test every entry of a sampled `ports x ports` matrix after extracting its
/// declared delay, `Z_ab exp(+i w tau_ab)`: the delay-free part that a
/// delayed fit realizes. A delay longer than the entry's first arrival
/// leaves an advanced, non-causal remainder and refuses.
///
/// # Errors
/// As for [`check_causality`], with the failing entry, plus mis-shaped
/// data or delays.
pub fn check_mimo_causality(
    omega: &[f64],
    data: &[Vec<C64>],
    ports: usize,
    delays: &[f64],
    opts: &CausalityOptions,
) -> Result<Vec<CausalityReport>, CausalityError> {
    let n = ports * ports;
    if ports == 0
        || data.len() != omega.len()
        || data.iter().any(|z| z.len() != n)
        || delays.len() != n
        || delays.iter().any(|t| !(t.is_finite() && *t >= 0.0))
    {
        return Err(CausalityError::BadSample {
            what: "matrix samples and delays must be ports x ports, delays finite and nonnegative",
        });
    }
    (0..n)
        .map(|e| {
            let entry: Vec<C64> = omega
                .iter()
                .zip(data)
                .map(|(&w, z)| {
                    let (s, c) = det::sin_cos(w * delays[e]);
                    z[e] * C64::new(c, s)
                })
                .collect();
            check_entry(omega, &entry, opts, e)
        })
        .collect()
}

#[allow(clippy::too_many_lines)] // one reconstruction and its three bounds
fn check_entry(
    omega: &[f64],
    data: &[C64],
    opts: &CausalityOptions,
    entry: usize,
) -> Result<CausalityReport, CausalityError> {
    let n_sub = opts.subtractions;
    if !(1..=8).contains(&n_sub)
        || !(opts.tail_bound_ratio >= 1.0 && opts.tail_bound_ratio.is_finite())
        || !(opts.relative_data_tolerance >= 0.0 && opts.relative_data_tolerance.is_finite())
    {
        return Err(CausalityError::BadOptions);
    }
    let bad = |what| Err(CausalityError::BadSample { what });
    let len = omega.len();
    if data.len() != len || len < 8 * (n_sub + 1) {
        return bad("need one sample per frequency and at least 8 (subtractions + 1)");
    }
    if data.iter().any(|z| !(z.re.is_finite() && z.im.is_finite())) {
        return bad("samples must be finite");
    }
    let step = omega[1] - omega[0];
    let top = omega[len - 1];
    if !(step > 0.0 && top.is_finite())
        || omega
            .iter()
            .enumerate()
            .any(|(j, &w)| (w - (omega[0] + j as f64 * step)).abs() > 1.0e-9 * top)
    {
        return bad("the grid must be uniform and increasing");
    }
    let dc = omega[0] == 0.0;
    if !dc && (omega[0] - 0.5 * step).abs() > 1.0e-9 * step {
        return bad("the grid must start at 0 or at half a step");
    }
    let scale = data.iter().map(|z| z.abs()).fold(0.0, f64::max);
    if scale == 0.0 {
        return Ok(CausalityReport {
            tested: 0,
            worst_ratio: 0.0,
            worst_omega: 0.0,
            max_relative_error: 0.0,
        });
    }
    if dc && data[0].im.abs() > 1.0e-12 * scale {
        return bad("a real response is real at DC");
    }
    // Hermitian extension onto a symmetric grid with trapezoid (DC grid) or
    // midpoint (half-step grid) weights over [-half_width, half_width].
    let mirror = if dc { len - 1 } else { len };
    let mut x = Vec::with_capacity(2 * len);
    let mut h = Vec::with_capacity(2 * len);
    for j in (usize::from(dc)..len).rev() {
        x.push(-omega[j]);
        h.push(data[j].conj());
    }
    x.extend_from_slice(omega);
    h.extend_from_slice(data);
    let mut weights = vec![step; x.len()];
    if dc {
        weights[0] = 0.5 * step;
        *weights.last_mut().expect("nonempty grid") = 0.5 * step;
    }
    let half_width = if dc { top } else { top + 0.5 * step };
    let sub: Vec<usize> = (0..n_sub)
        .map(|k| (k + 1) * (len - 1) / (n_sub + 1))
        .collect();
    let nodes: Vec<f64> = sub.iter().map(|&j| omega[j]).collect();
    let tail = opts.tail_bound_ratio * scale;
    let noise = opts.relative_data_tolerance * scale;
    let pi = core::f64::consts::PI;
    let mut report = CausalityReport {
        tested: 0,
        worst_ratio: 0.0,
        worst_omega: 0.0,
        max_relative_error: 0.0,
    };
    let mut worst = (0.0, 0.0);
    let mut r = vec![C64::ZERO; x.len()];
    for m in 1..len {
        if omega[m] > 0.9 * half_width || sub.contains(&m) {
            continue;
        }
        // Poles of the integrand: the subtraction points and omega[m].
        let s_idx: Vec<usize> = sub.iter().copied().chain([m]).collect();
        let s: Vec<f64> = s_idx.iter().map(|&j| omega[j]).collect();
        let s_sym: Vec<usize> = s_idx.iter().map(|&j| mirror + j).collect();
        let wprime: Vec<f64> = (0..=n_sub)
            .map(|p| {
                (0..=n_sub)
                    .filter(|&q| q != p)
                    .map(|q| s[p] - s[q])
                    .product()
            })
            .collect();
        let c: Vec<C64> = s_idx
            .iter()
            .zip(&wprime)
            .map(|(&j, &wp)| data[j].scale(1.0 / wp))
            .collect();
        // Regular part r = H / W - SUM c_p / (x - s_p), its limit at the
        // poles, and the first-order data gain of every term.
        let mut gain = 0.0;
        for i in 0..x.len() {
            let mut value;
            let mut g;
            if let Some(p) = s_sym.iter().position(|&k| k == i) {
                let derivative = (h[i + 1] - h[i - 1]).scale(0.5 / step);
                let inv: f64 = (0..=n_sub)
                    .filter(|&q| q != p)
                    .map(|q| 1.0 / (s[p] - s[q]))
                    .sum();
                value = (derivative - h[i].scale(inv)).scale(1.0 / wprime[p]);
                g = (1.0 / step + inv.abs()) / wprime[p].abs();
                for q in (0..=n_sub).filter(|&q| q != p) {
                    value = value - c[q].scale(1.0 / (s[p] - s[q]));
                    g += 1.0 / (wprime[q] * (s[p] - s[q])).abs();
                }
            } else {
                let big_w: f64 = s.iter().map(|&sp| x[i] - sp).product();
                value = h[i].scale(1.0 / big_w);
                g = 1.0 / big_w.abs();
                for p in 0..=n_sub {
                    value = value - c[p].scale(1.0 / (x[i] - s[p]));
                    g += 1.0 / (wprime[p] * (x[i] - s[p])).abs();
                }
            }
            r[i] = value;
            gain += weights[i] * g;
        }
        let mut integral = r
            .iter()
            .zip(&weights)
            .fold(C64::ZERO, |acc, (v, &wi)| acc + v.scale(wi));
        for p in 0..=n_sub {
            let log = det::ln((half_width - s[p]) / (half_width + s[p]));
            integral = integral + c[p].scale(log);
            gain += log.abs() / wprime[p].abs();
        }
        let z_m: f64 = nodes.iter().map(|&wk| omega[m] - wk).product();
        let mut interpolant = C64::ZERO;
        let mut lagrange = 0.0;
        for k in 0..n_sub {
            let l: f64 = (0..n_sub)
                .filter(|&j| j != k)
                .map(|j| (omega[m] - nodes[j]) / (nodes[k] - nodes[j]))
                .product();
            interpolant = interpolant + data[sub[k]].scale(l);
            lagrange += l.abs();
        }
        let reconstructed = interpolant + C64::new(0.0, z_m / pi) * integral;
        let error = (reconstructed - data[m]).abs();
        let second: f64 = r
            .windows(3)
            .map(|v| (v[2] - v[1].scale(2.0) + v[0]).abs())
            .sum();
        let discretization = 2.0 * z_m.abs() / pi * step * second / 12.0;
        let truncation = tail * z_m.abs() / pi * tail_integral(half_width, &nodes, omega[m]);
        let data_term = noise * (1.0 + lagrange + z_m.abs() / pi * gain);
        let bound = truncation + discretization + data_term;
        let ratio = error / bound;
        report.tested += 1;
        report.max_relative_error = report.max_relative_error.max(error / scale);
        if ratio > report.worst_ratio {
            report.worst_ratio = ratio;
            report.worst_omega = omega[m];
            worst = (error, bound);
        }
    }
    if report.worst_ratio > 1.0 {
        return Err(CausalityError::NonCausal {
            entry,
            omega: report.worst_omega,
            error: worst.0,
            bound: worst.1,
        });
    }
    Ok(report)
}

/// `INT_{|x| > W} dx / (|Z(x)| |x - w|)` with `x = W / u`, by composite
/// Simpson on `u in [0, 1]`; smooth there because `w <= 0.9 W` and the
/// subtraction points lie below `W`.
fn tail_integral(half_width: f64, nodes: &[f64], w: f64) -> f64 {
    let n = nodes.len() as i32;
    let f = |u: f64| {
        let upper: f64 = nodes
            .iter()
            .map(|&k| (half_width - u * k).abs())
            .product::<f64>()
            * (half_width - u * w).abs();
        let lower: f64 =
            nodes.iter().map(|&k| half_width + u * k).product::<f64>() * (half_width + u * w);
        half_width * det::powi(u, n - 1) * (1.0 / upper + 1.0 / lower)
    };
    let panels = 1000;
    let h = 1.0 / f64::from(panels);
    let mut sum = f(0.0) + f(1.0);
    for i in 1..panels {
        sum += if i % 2 == 1 { 4.0 } else { 2.0 } * f(f64::from(i) * h);
    }
    sum * h / 3.0
}

#[cfg(test)]
mod tests {
    use super::*;

    const TAU: f64 = core::f64::consts::TAU;

    fn opts() -> CausalityOptions {
        CausalityOptions {
            subtractions: 4,
            tail_bound_ratio: 2.0,
            relative_data_tolerance: 1.0e-12,
        }
    }

    /// `0.5 + 1500/(s + 2000)` plus two resonant conjugate pairs: causal.
    fn rational(w: f64) -> C64 {
        let s = C64::new(0.0, w);
        let mut h = C64::from_re(0.5) + C64::from_re(1500.0) * (s + C64::from_re(2000.0)).recip();
        for (p, r) in [
            (C64::new(-300.0, TAU * 400.0), C64::new(200.0, 50.0)),
            (C64::new(-800.0, TAU * 1500.0), C64::new(-300.0, 400.0)),
        ] {
            h = h + r * (s - p).recip() + r.conj() * (s - p.conj()).recip();
        }
        h
    }

    fn delayed(w: f64, tau: f64) -> C64 {
        let (s, c) = det::sin_cos(w * tau);
        rational(w) * C64::new(c, -s)
    }

    fn grid(dc: bool) -> Vec<f64> {
        let step = TAU * 5000.0 / 1024.0;
        (0..1025)
            .map(|j| (f64::from(j) + if dc { 0.0 } else { 0.5 }) * step)
            .collect()
    }

    fn report(label: &str, result: &Result<CausalityReport, CausalityError>) {
        println!(
            "{{\"suite\":\"fs-vfit-causality\",\"case\":\"{label}\",\"result\":\"{result:?}\"}}"
        );
    }

    /// G0: causal responses pass on both grids, delayed ones included.
    #[test]
    fn causal_and_delayed_responses_pass() {
        for dc in [true, false] {
            let omega = grid(dc);
            for tau in [0.0, 5.0e-4] {
                let data: Vec<C64> = omega.iter().map(|&w| delayed(w, tau)).collect();
                let result = check_causality(&omega, &data, &opts());
                report(&format!("causal dc={dc} tau={tau}"), &result);
                let passed = result.unwrap();
                assert!(passed.tested > 700 && passed.worst_ratio < 1.0);
            }
        }
    }

    /// G0: advances on the band's time scale (1/f_max = 0.2 ms and longer)
    /// refuse, and so does data left on the `exp(-i w t)` axis (conjugated).
    /// The anti-causal part enters the error as twice its miss from the
    /// subtraction interpolant, so a much shorter advance, smooth across the
    /// band, is not resolved; the 50 us probe is printed, not asserted.
    #[test]
    fn advanced_and_conjugated_responses_refuse() {
        let omega = grid(false);
        for tau in [5.0e-4, 2.0e-4] {
            let data: Vec<C64> = omega.iter().map(|&w| delayed(w, -tau)).collect();
            let result = check_causality(&omega, &data, &opts());
            report(&format!("advanced tau={tau}"), &result);
            assert!(matches!(
                result,
                Err(CausalityError::NonCausal { entry: 0, .. })
            ));
        }
        let probe: Vec<C64> = omega.iter().map(|&w| delayed(w, -5.0e-5)).collect();
        report(
            "advanced tau=5e-5 (below band resolution)",
            &check_causality(&omega, &probe, &opts()),
        );
        let data: Vec<C64> = omega.iter().map(|&w| rational(w).conj()).collect();
        let result = check_causality(&omega, &data, &opts());
        report("conjugated", &result);
        assert!(matches!(result, Err(CausalityError::NonCausal { .. })));
    }

    /// G0: a pulsating sphere's radiation impedance `rho c (s a/c)/(1 + s a/c)`
    /// is causal; the same data on the `exp(-i w t)` axis is not.
    #[test]
    fn sphere_radiation_impedance_and_its_convention_error() {
        let omega = grid(true);
        let (rho_c, a, c) = (413.0, 0.05, 343.0);
        let z = |w: f64| {
            let sa = C64::new(0.0, w * a / c);
            (sa * (C64::ONE + sa).recip()).scale(rho_c)
        };
        let causal: Vec<C64> = omega.iter().map(|&w| z(w)).collect();
        let result = check_causality(&omega, &causal, &opts());
        report("sphere", &result);
        assert!(result.unwrap().worst_ratio < 1.0);
        let wrong: Vec<C64> = causal.iter().map(|v| v.conj()).collect();
        let result = check_causality(&omega, &wrong, &opts());
        report("sphere conjugated", &result);
        assert!(matches!(result, Err(CausalityError::NonCausal { .. })));
    }

    /// G0: matrix entries are tested after their declared delays; an
    /// over-extracted delay names its entry.
    #[test]
    fn mimo_entries_with_declared_delays() {
        let omega = grid(false);
        let tau = 2.0e-4;
        let data: Vec<Vec<C64>> = omega
            .iter()
            .map(|&w| {
                let m = delayed(w, tau).scale(0.3);
                vec![rational(w), m, m, rational(w).scale(2.0)]
            })
            .collect();
        let reports = check_mimo_causality(&omega, &data, 2, &[0.0, tau, tau, 0.0], &opts());
        println!(
            "{{\"suite\":\"fs-vfit-causality\",\"case\":\"mimo\",\"result\":\"{reports:?}\"}}"
        );
        assert_eq!(reports.unwrap().len(), 4);
        let over =
            check_mimo_causality(&omega, &data, 2, &[0.0, 2.0 * tau, 2.0 * tau, 0.0], &opts());
        assert!(matches!(
            over,
            Err(CausalityError::NonCausal { entry: 1, .. })
        ));
        assert!(check_mimo_causality(&omega, &data, 2, &[0.0, -tau, tau, 0.0], &opts()).is_err());
    }

    /// G0: grids, data and options outside the documented ranges refuse.
    #[test]
    fn untestable_inputs_refuse() {
        let omega = grid(false);
        let data: Vec<C64> = omega.iter().map(|&w| rational(w)).collect();
        let mut shifted = omega.clone();
        shifted.iter_mut().for_each(|w| *w += 1.0);
        assert!(check_causality(&shifted, &data, &opts()).is_err());
        let mut uneven = omega.clone();
        uneven[10] += 1.0;
        assert!(check_causality(&uneven, &data, &opts()).is_err());
        assert!(check_causality(&omega[..20], &data[..20], &opts()).is_err());
        let dc = grid(true);
        let mut complex_dc: Vec<C64> = dc.iter().map(|&w| rational(w)).collect();
        complex_dc[0] = C64::new(1.0, 0.5);
        assert!(check_causality(&dc, &complex_dc, &opts()).is_err());
        for bad in [
            CausalityOptions {
                subtractions: 0,
                ..opts()
            },
            CausalityOptions {
                subtractions: 9,
                ..opts()
            },
            CausalityOptions {
                tail_bound_ratio: 0.5,
                ..opts()
            },
            CausalityOptions {
                relative_data_tolerance: -1.0,
                ..opts()
            },
        ] {
            assert_eq!(
                check_causality(&omega, &data, &bad),
                Err(CausalityError::BadOptions)
            );
        }
        let zero = vec![C64::ZERO; omega.len()];
        assert_eq!(check_causality(&omega, &zero, &opts()).unwrap().tested, 0);
    }
}
