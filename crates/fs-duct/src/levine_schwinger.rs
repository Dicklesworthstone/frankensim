//! Exact thin-wall, unflanged circular-pipe Wiener–Hopf termination.
//!
//! Levine & Schwinger, Phys. Rev. 73 (1948), 383–406, p.383 and Eq.(V.16),
//! DOI:10.1103/PhysRev.73.383. This evaluates their integrals, not a fitted
//! reflection filter or the long-wave radiation-impedance asymptote.
//!
//! Writing K=ka, the two end-correction integrals have weights
//! `1/[x sqrt(K²-x²)]` on (0,K) and `1/[x sqrt(x²+K²)]` on (0,infinity).
//! The numerators are ln[pi J1 sqrt(J1²+Y1²)] and -ln(2 I1 K1),
//! respectively; their sum is divided by pi. The logarithm of |R| is
//! -2K/pi times the first weight integrated against atan2(J1,-Y1).
//! In particular the second integral does NOT start at K. The propagating
//! integral uses x=K sin(theta), theta=(pi/2)(3u²-2u³); the infinite tail
//! uses x=(u/(1-u))^4 to
//! flatten both logarithmic endpoints without truncating the infinite domain.
//! Every quadrature is fixed-node Gauss–Legendre and refinable through `nodes`.

use crate::{DuctError, modal};
use fs_material::gas::GasState;
use fs_math::{c64::C64, dd::Dd, det};
use std::{
    f64::consts::{FRAC_PI_2, LN_2, PI},
    sync::LazyLock,
};

/// Smallest admitted fixed Gauss–Legendre rule.
pub const LS_MIN_NODES: usize = 16;
/// Work cap: each of the bounded nested rules uses at most this many nodes.
pub const LS_MAX_NODES: usize = 1024;
/// First positive J1 root, computed from the shared modal root finder.
/// This is a lazy value, rather than a transcribed cutoff constant.
pub static LS_MAX_KA: LazyLock<f64> = LazyLock::new(|| {
    let root = modal::j1_roots(1).expect("one J1 root is inside the modal cap")[0];
    let mut low = root - 1e-12;
    let mut high = root + 1e-12;
    for _ in 0..64 {
        let middle = 0.5 * (low + high);
        if middle == low || middle == high {
            break;
        }
        if modal::bessel_j1(middle) > 0.0 {
            low = middle;
        } else {
            high = middle;
        }
    }
    low
});

/// Explicit quadrature resolution. There is no implicit accuracy default.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct LsQuadrature {
    /// Actual nodes per rule, inclusive range `LS_MIN_NODES..=LS_MAX_NODES`.
    pub nodes: usize,
}

/// Thin-wall dominant-axisymmetric-mode solution under exp(-i omega t).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct LsSolution {
    /// Dimensionless acoustic radius, omega*a/c.
    pub ka: f64,
    /// Magnitude of the pressure/velocity-potential reflection coefficient.
    pub r_abs: f64,
    /// Log magnitude, retained for cancellation-safe low-frequency resistance.
    pub log_r_abs: f64,
    /// Acoustic end correction divided by radius, referenced to the mouth plane.
    pub l_over_a: f64,
    /// Actual quadrature resolution used.
    pub nodes: usize,
}

impl LsSolution {
    /// Complex pressure reflection `R = -|R| exp(2 i ka l/a)`.
    #[must_use]
    pub fn reflection(self) -> C64 {
        let phase = 2.0 * self.ka * self.l_over_a;
        C64::new(-self.r_abs * det::cos(phase), -self.r_abs * det::sin(phase))
    }

    /// Radiation impedance divided by rho*c/(pi*a²).
    /// The real part uses expm1, so 1-|R|² is not lost near zero frequency.
    #[must_use]
    pub fn normalized_impedance(self) -> C64 {
        let phase = 2.0 * self.ka * self.l_over_a;
        let c = det::cos(phase);
        let s = det::sin(phase);
        let denominator =
            (1.0 + self.r_abs * c) * (1.0 + self.r_abs * c) + (self.r_abs * s) * (self.r_abs * s);
        C64::new(
            -det::expm1(2.0 * self.log_r_abs) / denominator,
            -2.0 * self.r_abs * s / denominator,
        )
    }
}

// Legendre roots and weights on (0,1), with symmetric construction and fixed
// summation order. No endpoints are sampled, so logarithmic limits stay limits.
// Shared with the flanged-pipe radiation quadrature.
pub(crate) fn gauss_legendre(n: usize) -> Vec<(f64, f64)> {
    let mut rule = vec![(0.0, 0.0); n];
    for i in 0..n.div_ceil(2) {
        let mut z = det::cos(PI * (i as f64 + 0.75) / (n as f64 + 0.5));
        let mut derivative = 0.0;
        for _ in 0..32 {
            let mut p0 = 1.0;
            let mut p1 = z;
            for j in 2..=n {
                let jf = j as f64;
                let p = ((2.0 * jf - 1.0) * z * p1 - (jf - 1.0) * p0) / jf;
                p0 = p1;
                p1 = p;
            }
            derivative = n as f64 * (z * p1 - p0) / (z * z - 1.0);
            let step = p1 / derivative;
            z -= step;
            if step.abs() < 2.0 * f64::EPSILON {
                break;
            }
        }
        let weight = 1.0 / ((1.0 - z * z) * derivative * derivative);
        rule[i] = ((1.0 - z) * 0.5, weight);
        rule[n - 1 - i] = ((1.0 + z) * 0.5, weight);
    }
    rule
}

// Local cancellation-safe log(1+d), not a low-frequency model fit.
fn log_one_plus(d: f64) -> f64 {
    if d.abs() > 0.125 {
        return det::ln(1.0 + d);
    }
    let r = d / (2.0 + d);
    let rr = r * r;
    let mut power = r;
    let mut sum = r;
    for k in 1..16 {
        power *= rr;
        sum += power / (2 * k + 1) as f64;
    }
    2.0 * sum
}

// DLMF 10.8.1, order one. The harmonic-number series supplies the logarithmic
// Y1 part. The ordinary J1 value uses the shared owner; its scaled near-zero
// series is retained only where subtracting 1 from 2 J1/x loses the logarithm.
fn propagating_kernel(x: f64) -> (f64, f64) {
    if x == 0.0 {
        return (0.0, 0.0);
    }
    let q = -0.25 * x * x;
    let mut term = Dd::from_f64(0.5);
    let mut remainder = Dd::ZERO;
    let mut harmonic = 0.0;
    const GAMMA: f64 = 0.5772156649015329;
    let mut psi_sum = term * Dd::from_f64(1.0 - 2.0 * GAMMA);
    for k in 1..64 {
        let kf = k as f64;
        term = term * Dd::from_f64(q / (kf * (kf + 1.0)));
        harmonic += 1.0 / kf;
        remainder = remainder + term;
        psi_sum = psi_sum + term * Dd::from_f64(2.0 * harmonic + 1.0 / (kf + 1.0) - 2.0 * GAMMA);
        if term.hi.abs() < 1e-32 {
            break;
        }
    }
    let j_over_x = if x < 1e-4 {
        0.5 + remainder.to_f64()
    } else {
        modal::bessel_j1(x) / x
    };
    let delta_y = x * x * (-(det::ln(x) - LN_2) * j_over_x + 0.5 * psi_sum.to_f64());
    let negative_y_scaled = 1.0 + delta_y;
    let numerator = 0.5 * PI * x * x * j_over_x;
    // atan2 unwraps the arctangent continuously through the Y1 zero. A plain
    // atan(-J1/Y1) would jump by pi and give gain above ka~2.197.
    let phase = det::atan2(numerator, negative_y_scaled);
    let norm = det::hypot(negative_y_scaled, numerator);
    let log_argument = if x < 0.25 {
        let norm_delta = (2.0 * delta_y + delta_y * delta_y + numerator * numerator) / (norm + 1.0);
        log_one_plus(2.0 * remainder.to_f64() + 2.0 * j_over_x * norm_delta)
    } else {
        det::ln(2.0 * j_over_x * norm)
    };
    (phase, log_argument)
}

// -ln(2 I1 K1). At small arguments, DLMF 10.31.1 is arranged around
// 2 I1 K1 = 1+delta before taking the logarithm. This retains x² ln(x).
// For x>=1, K1 is evaluated by its scaled DLMF 10.32.9 integral:
// K1(x)e^x = integral_0^infinity 2 e^(-u²)(1+u²/x)/sqrt(2x+u²) du,
// after u=sqrt(2x)sinh(t/2). The u>12 tail is <1e-60 for x>=1.
fn log_inverse_ik(x: f64, rule: &[(f64, f64)]) -> f64 {
    const GAMMA: f64 = 0.5772156649015329;
    if x < 1.0 {
        let q = 0.25 * x * x;
        let mut term = Dd::from_f64(0.5);
        let mut remainder = Dd::ZERO;
        let mut psi_sum = term * Dd::from_f64(1.0 - 2.0 * GAMMA);
        let mut harmonic = 0.0;
        for k in 1..64 {
            let kf = k as f64;
            term = term * Dd::from_f64(q / (kf * (kf + 1.0)));
            harmonic += 1.0 / kf;
            remainder = remainder + term;
            psi_sum =
                psi_sum + term * Dd::from_f64(2.0 * harmonic + 1.0 / (kf + 1.0) - 2.0 * GAMMA);
            if term.hi.abs() < 1e-32 {
                break;
            }
        }
        let scaled_i = 0.5 + remainder.to_f64();
        let delta_k = x * x * ((det::ln(x) - LN_2) * scaled_i - 0.5 * psi_sum.to_f64());
        let delta = 2.0 * remainder.to_f64() + 2.0 * scaled_i * delta_k;
        return -log_one_plus(delta);
    }
    if x >= 32.0 {
        // DLMF 10.40.1/2. The 24-term expansion at x>=32 is far past
        // binary64 precision before its least term; no exp(+x) is formed.
        let mut term = 1.0;
        let mut i_sum = 1.0;
        let mut k_sum = 1.0;
        for k in 1..=24 {
            let odd = (2 * k - 1) as f64;
            term *= (4.0 - odd * odd) / (8.0 * k as f64 * x);
            k_sum += term;
            i_sum += if k % 2 == 0 { term } else { -term };
        }
        return det::ln(x) - det::ln(i_sum * k_sum);
    }
    let mut term = 0.5 * x;
    let mut i_sum = term;
    for k in 1..128 {
        let kf = k as f64;
        term *= 0.25 * x * x / (kf * (kf + 1.0));
        i_sum += term;
        if term < 1e-17 * i_sum {
            break;
        }
    }
    let i_scaled = i_sum * det::exp(-x);
    let mut k_scaled = Dd::ZERO;
    for &(v, w) in rule {
        let u = 12.0 * v;
        let uu = u * u;
        k_scaled = k_scaled
            + Dd::from_f64(w * 24.0 * det::exp(-uu) * (1.0 + uu / x) / det::sqrt(2.0 * x + uu));
    }
    -det::ln(2.0 * i_scaled * k_scaled.to_f64())
}

/// Evaluate the original LS Wiener–Hopf integrals below the first higher
/// axisymmetric cutoff. `nodes` is a numerical resolution, not a certificate.
pub fn levine_schwinger(ka: f64, q: LsQuadrature) -> Result<LsSolution, DuctError> {
    if !ka.is_finite() || ka <= 0.0 {
        return Err(DuctError::BadParameter {
            what: "Levine-Schwinger ka must be finite and positive",
        });
    }
    if ka >= *LS_MAX_KA {
        return Err(DuctError::LevineSchwingerCutoff {
            ka,
            cutoff: *LS_MAX_KA,
        });
    }
    if !(LS_MIN_NODES..=LS_MAX_NODES).contains(&q.nodes) {
        return Err(DuctError::BadParameter {
            what: "Levine-Schwinger nodes outside [16, 1024]",
        });
    }
    let rule = gauss_legendre(q.nodes);
    let mut modulus = Dd::ZERO;
    let mut length = Dd::ZERO;
    for &(u, w) in &rule {
        // Smooth both angular endpoints. At zero the length kernel is
        // theta*ln(theta): this map changes it to u³*ln(u), avoiding the
        // slow N^-4 convergence of an ordinary linear theta rule.
        let theta = FRAC_PI_2 * u * u * (3.0 - 2.0 * u);
        let jacobian = 3.0 * PI * u * (1.0 - u);
        let s = det::sin(theta);
        let x = ka * s;
        let (phase, logarithm) = propagating_kernel(x);
        modulus = modulus + Dd::from_f64(w * 6.0 * u * (1.0 - u) * phase / s);
        if x > 0.0 {
            length = length + Dd::from_f64(w * jacobian * logarithm / x);
        }
        let ratio = u / (1.0 - u);
        let xx = (ratio * ratio) * (ratio * ratio);
        let tail = 4.0 * log_inverse_ik(xx, &rule) / (u * (1.0 - u) * det::hypot(xx, ka));
        length = length + Dd::from_f64(w * tail);
    }
    // The -2/pi factor is included in the modulus rule above.
    let log_r_abs = -modulus.to_f64();
    let solution = LsSolution {
        ka,
        r_abs: det::exp(log_r_abs),
        log_r_abs,
        l_over_a: length.to_f64() / PI,
        nodes: q.nodes,
    };
    if !solution.l_over_a.is_finite() || !log_r_abs.is_finite() || log_r_abs > 0.0 {
        return Err(DuctError::BadParameter {
            what: "Levine-Schwinger integral left its finite passive domain",
        });
    }
    Ok(solution)
}

/// Radiation load in Pa s/m³ from the exact thin-wall LS reflection.
/// Gas properties are caller-supplied physical data, with radius in metres
/// and positive angular frequency in rad/s. No negative-frequency extension
/// or higher-mode load is inferred.
pub fn levine_schwinger_impedance(
    state: &GasState,
    radius: f64,
    omega: f64,
    q: LsQuadrature,
) -> Result<C64, DuctError> {
    crate::validate_gas(state)?;
    if [state.density, state.sound_speed, radius, omega]
        .iter()
        .any(|v| !v.is_finite() || *v <= 0.0)
    {
        return Err(DuctError::BadParameter {
            what: "Levine-Schwinger gas, radius and omega must be finite positive SI values",
        });
    }
    let ka = omega / state.sound_speed * radius;
    let zc = state.density * state.sound_speed / (PI * radius * radius);
    let z = levine_schwinger(ka, q)?.normalized_impedance().scale(zc);
    if !z.re.is_finite() || !z.im.is_finite() {
        return Err(DuctError::BadParameter {
            what: "Levine-Schwinger physical impedance is not representable",
        });
    }
    Ok(z)
}

#[cfg(test)]
mod tests {
    use super::*;
    use fs_material::gas::GasSpec;
    use fsci_runtime::RuntimeMode;
    use fsci_special::SpecialTensor;
    use std::time::Instant;

    // Independent 35-decimal-digit mpmath 1.3.0 tanh-sinh integration of the
    // ORIGINAL p.383/V.16 kernels, with mpmath J/Y/I/K, not these series/rules.
    // Eight dimensionless (ka, |R|, l/a) rows. The static value is 0.612701...;
    // the historical 0.6133 numerical value is corrected in Rienstra's
    // Fundamentals of Duct Acoustics, section on the unflanged open end:
    // https://sjoerdr.win.tue.nl/papers/VKI_Rienstra.pdf (delta=0.6127).
    const REFERENCE: [(f64, f64, f64); 8] = [
        (0.01, 0.9999500105996408, 0.6126826892082684),
        (0.1, 0.9950666976448565, 0.6109874925521944),
        (0.5, 0.8964413976209895, 0.5818364010713508),
        (1.0, 0.6951021732460863, 0.5274307840437657),
        (2.0, 0.3461760483069402, 0.4169443402868210),
        (3.0, 0.1543283211386729, 0.3082754595106364),
        (3.5, 0.1008588556982320, 0.2421584389555464),
        (3.8, 0.07777298733986563, 0.1707439491744777),
    ];

    fn air() -> GasState {
        GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 101325.0).unwrap()
    }
    fn relative(a: f64, b: f64) -> f64 {
        (a / b - 1.0).abs()
    }
    fn scalar(f: fn(&SpecialTensor, RuntimeMode) -> fsci_special::SpecialResult, x: f64) -> f64 {
        match f(&SpecialTensor::RealScalar(x), RuntimeMode::Strict).unwrap() {
            SpecialTensor::RealScalar(v) => v,
            other => panic!("non-scalar oracle: {other:?}"),
        }
    }
    // Independent double-exponential quadrature + Cephes-heritage Bessels.
    // This does not call the production GL rule or any production kernel.
    fn independent(ka: f64) -> (f64, f64) {
        let mut log_r = 0.0;
        let mut length = 0.0;
        let h = 1.0 / 64.0;
        for index in -384..=384 {
            let t = f64::from(index) * h;
            let et = det::exp(t);
            let emt = det::exp(-t);
            let sh = 0.5 * (et - emt);
            let ch = 0.5 * (et + emt);
            let v = FRAC_PI_2 * sh;
            if v.abs() > 20.0 {
                continue;
            }
            let cv = 0.5 * (det::exp(v) + det::exp(-v));
            let u = 0.5 * (1.0 + det::tanh(v));
            if u <= 0.0 || u >= 1.0 {
                continue;
            }
            let weight = h * PI * 0.25 * ch / (cv * cv);
            let theta = FRAC_PI_2 * u;
            let s = det::sin(theta);
            let x = ka * s;
            let j = scalar(fsci_special::j1, x);
            let y = scalar(fsci_special::y1, x);
            log_r -= weight * det::atan2(j, -y) / s;
            length += weight * FRAC_PI_2 * det::ln(PI * j * det::hypot(j, y)) / x;
            let ratio = u / (1.0 - u);
            let xx = (ratio * ratio) * (ratio * ratio);
            let product = 2.0 * fsci_special::i1e_scalar(xx) * fsci_special::k1e_scalar(xx);
            length += weight * 4.0 * (-det::ln(product)) / (u * (1.0 - u) * det::hypot(xx, ka));
        }
        (det::exp(log_r), length / PI)
    }

    #[test]
    fn ls_matches_reference_points() {
        for (ka, ref_r, ref_l) in REFERENCE {
            let start = Instant::now();
            let got = levine_schwinger(ka, LsQuadrature { nodes: 256 }).unwrap();
            assert!(relative(got.r_abs, ref_r) < 1e-3, "{got:?} vs {ref_r}");
            assert!(relative(got.l_over_a, ref_l) < 1e-3, "{got:?} vs {ref_l}");
            let (other_r, other_l) = independent(ka);
            let err = relative(got.r_abs, other_r).max(relative(got.l_over_a, other_l));
            assert!(
                err < 1e-8,
                "{got:?} vs independent ({other_r},{other_l}) error {err}"
            );
            eprintln!(
                "{{\"suite\":\"fs-duct\",\"case\":\"ls-reference\",\"inputs\":{{\"ka\":{ka},\"provenance\":\"Levine-Schwinger 1948 p.383/V.16; independent mpmath35 and fsci-special/tanh-sinh\"}},\"nodes\":256,\"observed\":{{\"r_abs\":{},\"l_over_a\":{}}},\"oracle\":{{\"r_abs\":{ref_r},\"l_over_a\":{ref_l}}},\"err\":{err},\"tol\":1e-8,\"verdict\":\"pass\",\"elapsed_us\":{}}}",
                got.r_abs,
                got.l_over_a,
                start.elapsed().as_micros()
            );
        }
    }

    #[test]
    fn ls_low_ka_limit() {
        // The exact integral has static limit 0.612701..., not the historical
        // four-digit approximation 0.6133. At ka=.05 the leading modulus-loss
        // asymptote differs by 0.00394305 relative; it is not a 1e-3 oracle.
        let mut previous = f64::INFINITY;
        for ka in [0.01, 0.005, 0.0025] {
            let s = levine_schwinger(ka, LsQuadrature { nodes: 256 }).unwrap();
            assert!(relative(s.l_over_a, 0.6127010359298) < 1e-3, "{s:?}");
            let err = relative(-det::expm1(s.log_r_abs), 0.5 * ka * ka);
            assert!(err < 1e-3 && err < previous, "{ka}: asymptotic error {err}");
            previous = err;
        }
        let s = levine_schwinger(0.05, LsQuadrature { nodes: 256 }).unwrap();
        assert!(relative(s.r_abs, 0.9987549288169632) < 1e-10);
        // Bound the leading expansion's measured finite-ka remainder;
        // it exceeds 1e-3 here even though the exact integral is accurate.
        let remainder = relative(-det::expm1(s.log_r_abs), 0.5 * 0.05 * 0.05);
        assert!((0.0039..0.0040).contains(&remainder), "{remainder}");
    }

    #[test]
    fn ls_historical_fit_retains_its_distinct_end_correction() {
        // The exact static value is also delta=0.6127 in Rienstra &
        // Hirschberg, An Introduction to Acoustics, sec.7.9 p.214:
        // https://sjoerdr.win.tue.nl/papers/boek.pdf (6 Sep 2026).
        // See V.V. Tupov's 2019 recalculation, AIP Conf. Proc.2195,020034,
        // "The calculation of channel end hole reflection ratio without
        // a flange", DOI:10.1063/1.5140134 (abstract reports 0.6127).
        let ka = 0.01;
        let exact = levine_schwinger(ka, LsQuadrature { nodes: 256 }).unwrap();
        let state = air();
        let radius = 0.008;
        let zc = state.density * state.sound_speed / (PI * radius * radius);
        let historical = crate::compact_radiation_impedance(
            crate::Termination::UnflangedOpen,
            &state,
            radius,
            ka * state.sound_speed / radius,
        )
        .unwrap();
        let historical_coefficient = -historical.im / (zc * ka);
        assert!(relative(historical_coefficient, 0.6133) < 1e-15);
        let disagreement = relative(exact.l_over_a, historical_coefficient);
        assert!((disagreement - 0.00100654).abs() < 1e-8, "{disagreement}");
    }

    #[test]
    fn ls_quadrature_ladder_converges() {
        let ka = 3.83;
        let mut values = Vec::new();
        let mut elapsed = Vec::new();
        for nodes in [32, 64, 128, 256] {
            let start = Instant::now();
            let s = levine_schwinger(ka, LsQuadrature { nodes }).unwrap();
            assert_eq!(s.nodes, nodes);
            values.push(s);
            elapsed.push(start.elapsed().as_micros());
        }
        for (s, elapsed_us) in values.iter().zip(elapsed) {
            let nodes = s.nodes;
            let finest = values[3];
            let error = (s.r_abs - finest.r_abs)
                .abs()
                .max((s.l_over_a - finest.l_over_a).abs());
            let verdict = if error < 1e-10 { "pass" } else { "fail" };
            eprintln!(
                "{{\"suite\":\"fs-duct\",\"case\":\"ls-ladder\",\"inputs\":{{\"ka\":{ka},\"provenance\":\"original Wiener-Hopf integrals\"}},\"axis\":\"quadrature-nodes\",\"rung\":{nodes},\"observed\":{{\"r_abs\":{},\"l_over_a\":{}}},\"oracle\":\"256-node rung\",\"err\":{error},\"tol\":1e-10,\"verdict\":\"{verdict}\",\"elapsed_us\":{elapsed_us}}}",
                s.r_abs, s.l_over_a
            );
        }
        let differences: Vec<_> = values
            .windows(2)
            .map(|v| {
                (v[0].r_abs - v[1].r_abs)
                    .abs()
                    .max((v[0].l_over_a - v[1].l_over_a).abs())
            })
            .collect();
        assert!(
            differences[1] <= differences[0] + 64.0 * f64::EPSILON,
            "{differences:?}"
        );
        assert!(
            differences[2] <= differences[1] + 64.0 * f64::EPSILON,
            "{differences:?}"
        );
        assert!(differences[2] < 1e-10, "{differences:?}");
        assert!(
            (values[0].l_over_a - values[3].l_over_a).abs() > 1e-12,
            "dead quadrature axis"
        );
    }

    #[test]
    fn ls_refuses_at_first_higher_mode_cutoff() {
        let q = LsQuadrature { nodes: 32 };
        let cap = *LS_MAX_KA;
        assert!(modal::bessel_j1(cap).abs() < 1e-13);
        assert!(levine_schwinger(cap * (1.0 - 1e-9), q).is_ok());
        for ka in [cap, cap + 1.0] {
            let error = levine_schwinger(ka, q).unwrap_err();
            assert!(matches!(error, DuctError::LevineSchwingerCutoff { .. }));
            assert!(error.to_string().contains("FS-DUCT-LS-CUTOFF"));
        }
        for ka in [0.0, -1.0, f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            assert!(
                levine_schwinger(ka, q)
                    .unwrap_err()
                    .to_string()
                    .contains("FS-DUCT-BAD-PARAMETER")
            );
        }
        assert!(
            levine_schwinger(
                1.0,
                LsQuadrature {
                    nodes: LS_MIN_NODES
                }
            )
            .is_ok()
        );
        assert!(
            levine_schwinger(
                1.0,
                LsQuadrature {
                    nodes: LS_MIN_NODES - 1
                }
            )
            .is_err()
        );
        assert!(
            levine_schwinger(
                1.0,
                LsQuadrature {
                    nodes: LS_MAX_NODES
                }
            )
            .is_ok()
        );
        assert!(
            levine_schwinger(
                1.0,
                LsQuadrature {
                    nodes: LS_MAX_NODES + 1
                }
            )
            .is_err()
        );
        for (radius, omega) in [
            (0.0, 1.0),
            (1.0, 0.0),
            (f64::NAN, 1.0),
            (1.0, f64::INFINITY),
        ] {
            assert!(levine_schwinger_impedance(&air(), radius, omega, q).is_err());
        }
    }

    #[test]
    fn ls_passive_and_convention_pinned() {
        let state = air();
        let radius = 0.008;
        for i in 1..=200 {
            let ka = 3.8 * f64::from(i) / 200.0;
            let q = LsQuadrature { nodes: 64 };
            let s = levine_schwinger(ka, q).unwrap();
            let z = levine_schwinger_impedance(&state, radius, ka * state.sound_speed / radius, q)
                .unwrap();
            assert!(s.r_abs <= 1.0 && s.r_abs >= 0.0);
            assert!(z.re >= 0.0, "{ka}: {z:?}");
            if ka < 1.0 {
                assert!(z.im < 0.0, "mass-like sign at {ka}: {z:?}");
            }
            let r = s.reflection();
            let direct = C64::new(1.0 + r.re, r.im) / C64::new(1.0 - r.re, -r.im);
            assert!((direct - s.normalized_impedance()).abs() < 1e-12);
        }
    }

    #[test]
    fn ls_low_ka_fit_model_form_band() {
        let error = |ka| {
            let exact = levine_schwinger(ka, LsQuadrature { nodes: 256 })
                .unwrap()
                .normalized_impedance();
            (C64::new(0.25 * ka * ka, -0.6133 * ka) - exact).abs() / exact.abs()
        };
        let mut previous = 0.0;
        for ka in [0.1, 0.25, 0.5, 1.0] {
            let start = Instant::now();
            let e = error(ka);
            assert!(e > previous);
            previous = e;
            eprintln!(
                "{{\"suite\":\"fs-duct\",\"case\":\"ls-model-form\",\"inputs\":{{\"ka\":{ka},\"provenance\":\"legacy 0.6133 low-ka fit vs exact WH\"}},\"observed\":{e},\"err\":{e},\"tol\":1e-3,\"verdict\":\"fail\",\"elapsed_us\":{}}}",
                start.elapsed().as_micros()
            );
        }
        let start = Instant::now();
        let mut low = 0.001;
        let mut high = 0.1;
        assert!(error(low) < 1e-3 && error(high) > 1e-3);
        for _ in 0..36 {
            let middle = 0.5 * (low + high);
            if error(middle) > 1e-3 {
                high = middle;
            } else {
                low = middle;
            }
        }
        eprintln!(
            "{{\"suite\":\"fs-duct\",\"case\":\"ls-fit-band-crossing\",\"inputs\":{{\"provenance\":\"first continuous crossing from low-ka limit\"}},\"observed\":{{\"ka_lo\":{low},\"ka_hi\":{high}}},\"err\":{{\"lo\":{},\"hi\":{}}},\"tol\":1e-3,\"verdict\":\"pass\",\"elapsed_us\":{}}}",
            error(low),
            error(high),
            start.elapsed().as_micros()
        );
    }

    #[test]
    fn ls_bitwise_repeat() {
        for ka in [0.01, 1.0, 2.0, 3.8] {
            let a = levine_schwinger(ka, LsQuadrature { nodes: 128 }).unwrap();
            let b = levine_schwinger(ka, LsQuadrature { nodes: 128 }).unwrap();
            assert_eq!(a.r_abs.to_bits(), b.r_abs.to_bits());
            assert_eq!(a.l_over_a.to_bits(), b.l_over_a.to_bits());
            assert_eq!(a.reflection().re.to_bits(), b.reflection().re.to_bits());
        }
    }
}
