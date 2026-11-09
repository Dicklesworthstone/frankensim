//! Delayed free-field mutual radiation between a flue's compact outlets.
//!
//! Each outlet's exterior pressure gains the free-field pressure of every
//! OTHER outlet's volume flow, p_m(t) = sum_n rho/(4 pi d_mn) dq_n/dt(t - d_mn/c)
//! (the off-diagonal of Z_mn = j w rho exp(-j k d)/(4 pi d)), realized
//! explicitly from accepted flows:
//!
//! p_m[k] = sum_n rho/(4 pi d) [s_n(k-N) - s_n(k-N-2)]/(2 dt),
//! s_n(j) = (1-f) q_n[j] + f q_n[j-1], N + f = d/(c dt) - 1.
//!
//! A central difference (zero at Nyquist, one sample of group delay) on a
//! linearly interpolated delay. With N >= 1 it reads only accepted flows, so
//! it enters each outlet's own terminal solve as a KNOWN series pressure and
//! no simultaneous multiport solve exists. The exact continuous term is not
//! passive together with the realized compact self loads (above about 120 kHz
//! at 288 kHz its resistance outgrows theirs); the central difference's
//! sin(w dt)/dt roll-off is the band limit, and passivity is checked, not
//! assumed: before any sample, Re Z = diag(realized self) + Re M must have a
//! Cholesky factor at every 5 Hz from 5 Hz to Nyquist, and M must stay within
//! 5 % of the exact term over the receiver's band.
//!
//! The `pritchard` form weights that term by the openings' own radiation
//! resistance: Pritchard's (1960) mutual impedance of small radiators, as
//! printed by Butler (2022, Eq. 1), Z_12 = R_11 [sin(kd) + j cos(kd)]/(kd) for
//! (ka)^2 << 1 and a/d << 1, with R_11 -> sqrt(R_ii R_jj) for unequal openings.
//! It is the free-field term times w_ij = sqrt(R_ii R_jj)/(rho w^2/(4 pi c)),
//! so Re Z = (rho w^2/(4 pi c)) diag(g) S diag(g) with S_ij = sinc(k d_ij) and
//! g_i^2 = R_ii/(rho w^2/(4 pi c)): positive semidefinite by the Schur product,
//! where compact unflanged self loads (whose resistance falls under a free-field
//! monopole's) make the plain free-field term indefinite once three or more
//! openings sit in a line within a wavelength. The weight is realized as a
//! symmetric FIR of half-length T = min(N - 1, 16) taps, from the self loads'
//! bilinear-realized resistance at the warped frequency, normalized to its
//! exact DC value and convolved into the pair's kernel; every tap still reads
//! an accepted flow. The same cold admission applies, its band gate against
//! the weighted exact term.
//!
//! Not modeled: scattering by the instrument body (the outlets sit in the
//! card's flat plane at their declared centers, in free field), the image of
//! a baffle, the mouth (it is not an outlet), and any mean flow.
use crate::flue_loop::FlueError;
use fs_math::c64::C64;
use fs_vfit::relaxation::{RelaxationImpedanceSpec, RelaxationSampling};
use fs_vfit::waveguide::WaveguideError;
use fs_vfit::waveguide::network::WaveguideNetwork;

/// Admission grid step [Hz].
pub const OUTLET_MUTUAL_GRID_HZ: f64 = 5.0;
/// Largest admitted relative error of the realized mutual term in band.
pub const OUTLET_MUTUAL_BAND_ERROR: f64 = 0.05;
/// Largest number of coupled outlets.
pub const MAX_MUTUAL_OUTLETS: usize = 64;
/// Largest retained flow history per outlet [samples].
pub const MAX_MUTUAL_HISTORY: usize = 1 << 16;
// Largest admitted whole delay, N + 3 <= MAX_MUTUAL_HISTORY - 1.
const MAX_DELAY: f64 = 65_532.0;
/// Largest half-length of a `pritchard` weight filter [samples].
pub const PRITCHARD_MAX_HALF_TAPS: usize = 16;
/// Smallest admitted half-length of a `pritchard` weight filter [samples].
pub const PRITCHARD_MIN_HALF_TAPS: usize = 8;
// Midpoints on (0, pi) of the weight filter's inverse transform.
const PRITCHARD_QUADRATURE: usize = 4096;

fn invalid(what: &'static str) -> FlueError { FlueError::Invalid { what } }

/// The mutual term a card row asks for.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum OutletMutualForm {
    /// The free-field monopole term j w rho exp(-j k d)/(4 pi d).
    FreeField,
    /// Pritchard's term, the free-field one times sqrt(R_ii R_jj)/(rho w^2/(4 pi c)).
    Pritchard,
}

/// A pair's `pritchard` weight filter and its convolved kernel.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct Weighted {
    /// Half-length T [samples].
    half: usize,
    /// Start of w[0..=T] in the weight store.
    weights: usize,
    /// Start of the 2T + 4 kernel taps in the tap store.
    taps: usize,
    /// Lag of the first kernel tap, N - T >= 1.
    first_lag: usize,
}

/// One coupled pair, in the declared outlet order (`outlets[0] < outlets[1]`).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct OutletMutualPair {
    /// Positions in the declared outlet list.
    pub outlets: [usize; 2],
    /// Center distance in the card's plane [m].
    pub distance_m: f64,
    /// Whole samples N of the interpolated delay, N >= 1.
    pub delay_samples: usize,
    /// Fraction f in [0, 1) of the interpolated delay.
    pub delay_fraction: f64,
    /// rho / (4 pi d) / (2 dt) [Pa s/m^3 per sample].
    gain: f64,
}

/// Cold admission numbers, retained for the render's record.
#[derive(Debug, Clone, PartialEq)]
pub struct OutletMutualRecord {
    /// Exterior density of the self loads [kg/m^3].
    pub density_kg_m3: f64,
    /// Exterior sound speed of the self loads [m/s].
    pub sound_speed_m_s: f64,
    /// Mechanical step [s].
    pub time_step_s: f64,
    /// Receiver band the 5 % error gate covers [Hz].
    pub band_hz: f64,
    /// Number of admission frequencies (5 Hz to below Nyquist).
    pub grid_points: usize,
    /// Smallest whole delay over the pairs [samples].
    pub min_delay_samples: usize,
    /// Largest relative error of the realized term against the exact one in band.
    pub max_band_error: f64,
    /// Frequency of that error [Hz].
    pub max_band_error_hz: f64,
    /// Smallest Cholesky pivot over its own diagonal entry, over the grid.
    pub min_pivot_ratio: f64,
    /// Frequency of that pivot [Hz].
    pub min_pivot_hz: f64,
    /// The admitted form; `max_band_error` is against its own exact term.
    pub form: OutletMutualForm,
    /// Smallest weight-filter half-length over the pairs (0 for `free-field`).
    pub min_weight_half_taps: usize,
    /// Largest relative distance of the realized term from the plain
    /// free-field term in band (equal to `max_band_error` for `free-field`).
    pub max_free_field_deviation: f64,
}

/// Explicit delayed mutual radiation: one flow ring per outlet.
#[derive(Debug, Clone)]
pub struct OutletMutual {
    nodes: Vec<usize>,
    pairs: Vec<OutletMutualPair>,
    len: usize,
    head: usize,
    rings: Vec<f64>,
    scratch: Vec<f64>,
    record: OutletMutualRecord,
    dt: f64,
    // `pritchard` only: one entry per pair, its weights and convolved taps.
    weighted: Vec<Weighted>,
    weights: Vec<f64>,
    taps: Vec<f64>,
}

/// lim Re Z(w)/w^2 at DC of a load whose resistance vanishes there as w^2:
/// no base resistance and no compliance term, at least one relaxation term.
fn dc_curvature(load: &RelaxationImpedanceSpec) -> Option<f64> {
    if load.base().resistance_pa_s_m3.abs() > 0.0 || !load.compliance_terms().is_empty() || load.terms().is_empty() {
        return None;
    }
    let c: f64 = load.terms().iter().map(|t| t.resistance_pa_s_m3 / (t.rate_per_s * t.rate_per_s)).sum();
    (c.is_finite() && c > 0.0).then_some(c)
}

/// Sum of R/p^4 over the relaxation terms of a load [`dc_curvature`] admits:
/// Re Z(w)/w^2 = sum R/(w^2 + p^2) falls from its DC value as -w^2 times this.
fn dc_quartic(load: &RelaxationImpedanceSpec) -> f64 {
    load.terms().iter().map(|t| {
        let p2 = t.rate_per_s * t.rate_per_s;
        t.resistance_pa_s_m3 / (p2 * p2)
    }).sum()
}

/// sqrt(Re Z_a Re Z_b)/(mono w^2) at angular frequency `w`, mono = rho/(4 pi c).
fn pritchard_weight(a: &RelaxationImpedanceSpec, b: &RelaxationImpedanceSpec, mono: f64, w: f64) -> f64 {
    fs_math::det::sqrt(self_resistance(a, w) * self_resistance(b, w)) / (mono * w * w)
}

/// Symmetric weight filter w[0..=half] of the bilinear-realized weight
/// W(theta) = pritchard_weight at (2/dt) tan(theta/2): its cosine series by the
/// midpoint rule on (0, pi), then the least change of the response (in L2 on
/// (0, pi)) that gives it W's exact DC value `dc` and exact DC curvature `bend`
/// = W''(0): w[0] += l1, w[m] += l1 - m^2 l2. The truncated series alone misses
/// that curvature, and a realized weight above the one the self resistances
/// imply leaves Re Z without a Cholesky factor wherever the exact matrix is
/// nearly singular (openings in a line, low frequency).
#[allow(clippy::too_many_arguments)]
fn pritchard_filter(a: &RelaxationImpedanceSpec, b: &RelaxationImpedanceSpec, mono: f64, dc: f64, bend: f64,
    dt: f64, half: usize, out: &mut [f64],
) {
    out.fill(0.0);
    #[allow(clippy::cast_precision_loss)]
    let q = PRITCHARD_QUADRATURE as f64;
    for k in 0..PRITCHARD_QUADRATURE {
        #[allow(clippy::cast_precision_loss)]
        let theta = core::f64::consts::PI * (k as f64 + 0.5) / q;
        let value = pritchard_weight(a, b, mono, 2.0 * fs_math::det::tan(0.5 * theta) / dt);
        let c1 = fs_math::det::cos(theta);
        let (mut prev, mut cur) = (c1, 1.0);
        for tap in out.iter_mut().take(half + 1) {
            *tap += value * cur;
            let next = 2.0 * c1 * cur - prev;
            (prev, cur) = (cur, next);
        }
    }
    for tap in out.iter_mut().take(half + 1) { *tap /= q; }
    // Value V = w[0] + 2 sum w[m], curvature K = -2 sum m^2 w[m]; solve
    // [2T+1, -2 S2; -2 S2, 2 S4] [l1; l2] = [dc - V; bend - K], S_k = sum m^k.
    let (mut value, mut curvature, mut s2, mut s4) = (out[0], 0.0, 0.0, 0.0);
    for (m, &tap) in out.iter().enumerate().take(half + 1).skip(1) {
        #[allow(clippy::cast_precision_loss)]
        let m2 = (m * m) as f64;
        value += 2.0 * tap;
        curvature -= 2.0 * m2 * tap;
        s2 += m2;
        s4 += m2 * m2;
    }
    #[allow(clippy::cast_precision_loss)]
    let t = (2 * half + 1) as f64;
    let det = 2.0 * t * s4 - 4.0 * s2 * s2;
    let (r1, r2) = (dc - value, bend - curvature);
    let l1 = (2.0 * s4 * r1 + 2.0 * s2 * r2) / det;
    let l2 = (2.0 * s2 * r1 + t * r2) / det;
    out[0] += l1;
    for (m, tap) in out.iter_mut().enumerate().take(half + 1).skip(1) {
        #[allow(clippy::cast_precision_loss)]
        let m2 = (m * m) as f64;
        *tap += l1 - m2 * l2;
    }
}

/// W(theta) = w[0] + 2 sum_m w[m] cos(m theta) of a symmetric filter.
fn filter_response(w: &[f64], theta: f64) -> f64 {
    let c1 = fs_math::det::cos(theta);
    let (mut prev, mut cur) = (c1, 1.0);
    let mut sum = 0.0;
    for (m, &tap) in w.iter().enumerate() {
        sum += if m == 0 { tap } else { 2.0 * tap * cur };
        let next = 2.0 * c1 * cur - prev;
        (prev, cur) = (cur, next);
    }
    sum
}

/// Bilinear-realized self resistance Re Z(w) of an admitted relaxation load.
fn self_resistance(load: &RelaxationImpedanceSpec, warped: f64) -> f64 {
    let w2 = warped * warped;
    let mut r = load.base().resistance_pa_s_m3;
    for term in load.terms() {
        let p2 = term.rate_per_s * term.rate_per_s;
        r += term.resistance_pa_s_m3 * (w2 / (w2 + p2));
    }
    for term in load.compliance_terms() {
        let p2 = term.rate_per_s * term.rate_per_s;
        r += term.resistance_pa_s_m3 * (p2 / (w2 + p2));
    }
    r
}

impl OutletMutualPair {
    /// The realized kernel at angular frequency `omega`, under exp(+j w t):
    /// gain (1 - z^-2) z^-N ((1-f) + f z^-1), z^-1 = exp(-j w dt).
    #[must_use]
    pub fn realized(&self, omega: f64, dt: f64) -> C64 {
        let theta = omega * dt;
        let (s1, c1) = fs_math::det::sin_cos(theta);
        let (s2, c2) = fs_math::det::sin_cos(2.0 * theta);
        #[allow(clippy::cast_precision_loss)] // delays are bounded by MAX_MUTUAL_HISTORY
        let (sn, cn) = fs_math::det::sin_cos(self.delay_samples as f64 * theta);
        let difference = C64::new(1.0 - c2, s2);
        let delay = C64::new(cn, -sn);
        let interp = C64::new((1.0 - self.delay_fraction) + self.delay_fraction * c1, -self.delay_fraction * s1);
        (difference * delay * interp).scale(self.gain)
    }

    /// The exact free-field term j w rho exp(-j w d/c)/(4 pi d).
    #[must_use]
    pub fn exact(&self, omega: f64, density: f64, sound_speed: f64) -> C64 {
        let (s, c) = fs_math::det::sin_cos(omega * self.distance_m / sound_speed);
        let k = omega * density / (4.0 * core::f64::consts::PI * self.distance_m);
        C64::new(k * s, k * c)
    }
}

impl OutletMutual {
    /// Admit the coupling of these outlets (terminal `nodes`, plane `centers`,
    /// realized self `loads`, all in the declared outlet order) for exterior
    /// `density` and `sound_speed`, step `dt` and receiver band `band_hz`.
    ///
    /// # Errors
    /// Fewer than two or more than 64 outlets, mismatched lists, a repeated
    /// node, nonfinite or nonpositive inputs, coincident centers, a pair closer
    /// than two samples of travel (not explicit), a load not sampled
    /// bilinear-midpoint, a history beyond 65536 samples, a band at or above
    /// Nyquist, a realized term more than 5 % from the exact one in band, or
    /// a Re Z without a Cholesky factor at some grid frequency (not passive).
    pub fn admit(nodes: &[usize], centers: &[[f64; 2]], loads: &[RelaxationImpedanceSpec],
        density: f64, sound_speed: f64, dt: f64, band_hz: f64, gate: &fs_exec::CancelGate,
    ) -> Result<Self, FlueError> {
        Self::admit_form(OutletMutualForm::FreeField, nodes, centers, loads, density, sound_speed, dt, band_hz, gate)
    }

    /// [`Self::admit`] for either form. `pritchard` also refuses a load whose
    /// resistance does not vanish as w^2 at DC (a base resistance or a
    /// compliance term), a pair closer than ten samples of travel (its weight
    /// filter would not fit inside the delay), and a realized term more than
    /// 5 % from the weighted exact term in band.
    ///
    /// # Errors
    /// As [`Self::admit`], plus the `pritchard` refusals above.
    #[allow(clippy::too_many_arguments, clippy::too_many_lines)]
    pub fn admit_form(form: OutletMutualForm, nodes: &[usize], centers: &[[f64; 2]], loads: &[RelaxationImpedanceSpec],
        density: f64, sound_speed: f64, dt: f64, band_hz: f64, gate: &fs_exec::CancelGate,
    ) -> Result<Self, FlueError> {
        let n = nodes.len();
        if n < 2 || n > MAX_MUTUAL_OUTLETS {
            return Err(invalid("outlet_mutual couples two to 64 outlets"));
        }
        if centers.len() != n || loads.len() != n {
            return Err(invalid("outlet_mutual needs one center and one self load per outlet"));
        }
        if (0..n).any(|i| nodes[i + 1..].contains(&nodes[i])) {
            return Err(invalid("outlet_mutual outlets must be distinct terminals"));
        }
        if ![density, sound_speed, dt, band_hz].iter().all(|v| v.is_finite() && *v > 0.0)
            || centers.iter().flatten().any(|v| !v.is_finite())
        {
            return Err(invalid("outlet_mutual needs a finite positive gas, clock and band and finite centers"));
        }
        let nyquist = 0.5 / dt;
        if band_hz >= nyquist {
            return Err(invalid("outlet_mutual band must lie below Nyquist"));
        }
        if loads.iter().any(|l| l.sampling() != RelaxationSampling::BilinearMidpoint) {
            return Err(invalid("outlet_mutual admits bilinear-midpoint outlet loads only"));
        }
        let mut pairs = Vec::new();
        pairs.try_reserve_exact(n * (n - 1) / 2).map_err(|_| invalid("outlet_mutual allocation failed"))?;
        let mut longest = 0usize;
        for i in 0..n {
            for j in i + 1..n {
                let d = fs_math::det::hypot(centers[i][0] - centers[j][0], centers[i][1] - centers[j][1]);
                if !d.is_finite() || d <= 0.0 {
                    return Err(invalid("outlet_mutual outlets must have distinct centers"));
                }
                let rest = d / (sound_speed * dt) - 1.0;
                if !rest.is_finite() || rest < 1.0 {
                    return Err(invalid("outlet_mutual pair is closer than two samples of travel: the kernel would not be explicit"));
                }
                if rest >= MAX_DELAY {
                    return Err(invalid("outlet_mutual pair delay exceeds the retained history"));
                }
                #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)] // 1 <= rest < 65532
                let whole = rest.floor() as usize;
                #[allow(clippy::cast_precision_loss)]
                let fraction = rest - whole as f64;
                longest = longest.max(whole);
                pairs.push(OutletMutualPair { outlets: [i, j], distance_m: d, delay_samples: whole,
                    delay_fraction: fraction,
                    gain: density / (4.0 * core::f64::consts::PI * d) / (2.0 * dt) });
            }
        }
        let mut len = longest + 3;
        let mono = density / (4.0 * core::f64::consts::PI * sound_speed);
        let (mut weighted, mut weights, mut taps) = (Vec::new(), Vec::new(), Vec::new());
        if form == OutletMutualForm::Pritchard {
            let mut curvature = Vec::with_capacity(n);
            // d ln(Re Z/w^2)/d(w^2) at DC, per load.
            let mut slope = Vec::with_capacity(n);
            for load in loads {
                let Some(c) = dc_curvature(load) else {
                    return Err(invalid("outlet_mutual pritchard needs self loads whose resistance vanishes as w^2 at DC (no base resistance or compliance term)"));
                };
                curvature.push(c);
                slope.push(-dc_quartic(load) / c);
            }
            for pair in &pairs {
                if pair.delay_samples < PRITCHARD_MIN_HALF_TAPS + 1 {
                    return Err(invalid("outlet_mutual pritchard pair is closer than ten samples of travel: its weight filter would not fit inside the delay"));
                }
                let half = (pair.delay_samples - 1).min(PRITCHARD_MAX_HALF_TAPS);
                let first_lag = pair.delay_samples - half;
                if pair.delay_samples + half + 3 > MAX_MUTUAL_HISTORY - 1 {
                    return Err(invalid("outlet_mutual pair delay exceeds the retained history"));
                }
                len = len.max(pair.delay_samples + half + 3);
                let [a, b] = pair.outlets;
                let entry = Weighted { half, weights: weights.len(), taps: taps.len(), first_lag };
                weights.try_reserve(half + 1).map_err(|_| invalid("outlet_mutual allocation failed"))?;
                weights.resize(entry.weights + half + 1, 0.0);
                let dc = fs_math::det::sqrt(curvature[a] * curvature[b]) / mono;
                // W = sqrt(f_a f_b)/mono, f = Re Z/w^2, w^2 = theta^2/dt^2 + O(theta^4):
                // W''(0) = dc (f_a'/f_a + f_b'/f_b)/dt^2.
                let bend = dc * (slope[a] + slope[b]) / (dt * dt);
                pritchard_filter(&loads[a], &loads[b], mono, dc, bend, dt, half, &mut weights[entry.weights..]);
                // Kernel taps at lags N - T + i: gain [(1-f), f, -(1-f), -f] at
                // lags N..N+3, convolved with w[|s|] at shifts s = -T..T.
                let f = pair.delay_fraction;
                let h4 = [1.0 - f, f, -(1.0 - f), -f];
                taps.try_reserve(2 * half + 4).map_err(|_| invalid("outlet_mutual allocation failed"))?;
                taps.resize(entry.taps + 2 * half + 4, 0.0);
                for (i, &h) in h4.iter().enumerate() {
                    for m in 0..=2 * half {
                        taps[entry.taps + i + m] += pair.gain * h * weights[entry.weights + m.abs_diff(half)];
                    }
                }
                weighted.push(entry);
            }
        }
        // Grid: 5 Hz steps strictly below Nyquist.
        #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
        let points = ((nyquist / OUTLET_MUTUAL_GRID_HZ).ceil() as usize).saturating_sub(1);
        let mut matrix = vec![0.0; n * n];
        let mut factor = vec![0.0; n * n];
        let (mut worst_err, mut worst_err_hz) = (0.0f64, 0.0);
        let mut worst_dev = 0.0f64;
        let (mut worst_pivot, mut worst_pivot_hz) = (f64::INFINITY, 0.0);
        for k in 1..=points {
            if k % 4096 == 0 && gate.is_requested() {
                return Err(invalid("outlet_mutual admission cancelled"));
            }
            #[allow(clippy::cast_precision_loss)]
            let hz = OUTLET_MUTUAL_GRID_HZ * k as f64;
            let omega = core::f64::consts::TAU * hz;
            let warped = 2.0 * fs_math::det::tan(0.5 * omega * dt) / dt;
            for (i, load) in loads.iter().enumerate() {
                matrix[i * n + i] = self_resistance(load, warped);
            }
            for (p, pair) in pairs.iter().enumerate() {
                let mut m = pair.realized(omega, dt);
                let [i, j] = pair.outlets;
                if let Some(entry) = weighted.get(p) {
                    m = m.scale(filter_response(&weights[entry.weights..=entry.weights + entry.half], omega * dt));
                }
                matrix[i * n + j] = m.re;
                matrix[j * n + i] = m.re;
                if hz <= band_hz {
                    let exact = pair.exact(omega, density, sound_speed);
                    let e = (m - exact).abs() / exact.abs();
                    if weighted.is_empty() {
                        if !e.is_finite() || e > worst_err { (worst_err, worst_err_hz) = (e, hz); }
                    } else {
                        if !e.is_finite() || e > worst_dev { worst_dev = e; }
                        let target = exact.scale(pritchard_weight(&loads[i], &loads[j], mono, omega));
                        let e = (m - target).abs() / target.abs();
                        if !e.is_finite() || e > worst_err { (worst_err, worst_err_hz) = (e, hz); }
                    }
                }
            }
            // Cholesky of the symmetric Re Z; every pivot must be positive.
            for j in 0..n {
                let mut pivot = matrix[j * n + j];
                for q in 0..j { pivot -= factor[j * n + q] * factor[j * n + q]; }
                let ratio = pivot / matrix[j * n + j];
                if !(pivot > 0.0) || !ratio.is_finite() {
                    return Err(invalid("outlet_mutual Re Z has no Cholesky factor at some frequency: not passive with these self loads"));
                }
                if ratio < worst_pivot { (worst_pivot, worst_pivot_hz) = (ratio, hz); }
                let diagonal = fs_math::det::sqrt(pivot);
                factor[j * n + j] = diagonal;
                for i in j + 1..n {
                    let mut v = matrix[i * n + j];
                    for q in 0..j { v -= factor[i * n + q] * factor[j * n + q]; }
                    factor[i * n + j] = v / diagonal;
                }
            }
        }
        if !(worst_err <= OUTLET_MUTUAL_BAND_ERROR) {
            return Err(invalid(if weighted.is_empty() {
                "outlet_mutual realized term is more than 5 % from the exact free-field term in band"
            } else {
                "outlet_mutual pritchard realized term is more than 5 % from its exact term in band"
            }));
        }
        let mut rings = Vec::new();
        rings.try_reserve_exact(n * len).map_err(|_| invalid("outlet_mutual allocation failed"))?;
        rings.resize(n * len, 0.0);
        let min_delay = pairs.iter().map(|p| p.delay_samples).min().unwrap_or(0);
        let min_half = weighted.iter().map(|w| w.half).min().unwrap_or(0);
        let deviation = if weighted.is_empty() { worst_err } else { worst_dev };
        Ok(Self { nodes: nodes.to_vec(), pairs, len, head: 0, rings, scratch: vec![0.0; n],
            record: OutletMutualRecord { density_kg_m3: density, sound_speed_m_s: sound_speed, time_step_s: dt,
                band_hz, grid_points: points, min_delay_samples: min_delay, max_band_error: worst_err,
                max_band_error_hz: worst_err_hz, min_pivot_ratio: worst_pivot, min_pivot_hz: worst_pivot_hz,
                form, min_weight_half_taps: min_half, max_free_field_deviation: deviation },
            dt, weighted, weights, taps })
    }

    /// The realized term of pair `index` at angular frequency `omega`: the
    /// pair's kernel, times its weight filter's response under `pritchard`.
    ///
    /// # Panics
    /// An `index` outside [`Self::pairs`].
    #[must_use]
    pub fn realized_term(&self, index: usize, omega: f64) -> C64 {
        let m = self.pairs[index].realized(omega, self.dt);
        match self.weighted.get(index) {
            Some(entry) => m.scale(filter_response(&self.weights[entry.weights..=entry.weights + entry.half], omega * self.dt)),
            None => m,
        }
    }

    /// Weight-filter half-length T of pair `index` (0 under `free-field`).
    #[must_use]
    pub fn weight_half_taps(&self, index: usize) -> usize { self.weighted.get(index).map_or(0, |w| w.half) }

    /// Smallest lag any kernel tap reads [samples]: N under `free-field`,
    /// N - T under `pritchard`; at least 1, so every read is an accepted flow.
    #[must_use]
    pub fn min_read_lag(&self) -> usize {
        if self.weighted.is_empty() {
            self.pairs.iter().map(|p| p.delay_samples).min().unwrap_or(0)
        } else {
            self.weighted.iter().map(|w| w.first_lag).min().unwrap_or(0)
        }
    }

    /// Terminal nodes, in the declared outlet order.
    #[must_use]
    pub fn nodes(&self) -> &[usize] { &self.nodes }
    /// Coupled pairs with their realized delays.
    #[must_use]
    pub fn pairs(&self) -> &[OutletMutualPair] { &self.pairs }
    /// Cold admission numbers.
    #[must_use]
    pub const fn record(&self) -> &OutletMutualRecord { &self.record }

    /// Flow q[now - lag] of outlet `outlet`, lag in 1..=len; zero before the start.
    fn past(&self, outlet: usize, lag: usize) -> f64 {
        self.rings[outlet * self.len + (self.head + self.len - lag) % self.len]
    }

    /// Exterior pressure of every outlet for the NEXT step [Pa], from the
    /// accepted flows only. `out` holds one value per outlet.
    pub fn exterior_pressures(&self, out: &mut [f64]) {
        out.fill(0.0);
        if !self.weighted.is_empty() {
            for (pair, entry) in self.pairs.iter().zip(&self.weighted) {
                let [a, b] = pair.outlets;
                let kernel = &self.taps[entry.taps..entry.taps + 2 * entry.half + 4];
                let (mut pa, mut pb) = (0.0, 0.0);
                for (i, &c) in kernel.iter().enumerate() {
                    pa += c * self.past(b, entry.first_lag + i);
                    pb += c * self.past(a, entry.first_lag + i);
                }
                out[a] += pa;
                out[b] += pb;
            }
            return;
        }
        for pair in &self.pairs {
            let (n, f) = (pair.delay_samples, pair.delay_fraction);
            let [a, b] = pair.outlets;
            let slope = |k: usize| (1.0 - f) * (self.past(k, n) - self.past(k, n + 2))
                + f * (self.past(k, n + 1) - self.past(k, n + 3));
            out[a] += pair.gain * slope(b);
            out[b] += pair.gain * slope(a);
        }
    }

    /// Retain one accepted step's outlet flows [m^3/s], in outlet order.
    pub fn push_flows(&mut self, flows: &[f64]) {
        for (k, &q) in flows.iter().enumerate().take(self.nodes.len()) {
            self.rings[k * self.len + self.head] = q;
        }
        self.head = (self.head + 1) % self.len;
    }

    /// Set every outlet's exterior pressure on the network for its next step.
    ///
    /// # Errors
    /// A nonfinite pressure or a terminal the network did not admit.
    pub(crate) fn apply(&mut self, physics: &mut WaveguideNetwork) -> Result<(), WaveguideError> {
        let mut out = core::mem::take(&mut self.scratch);
        self.exterior_pressures(&mut out);
        let result = self.nodes.iter().zip(&out)
            .try_for_each(|(&node, &p)| physics.set_terminal_exterior_pressure(node, p));
        self.scratch = out;
        result
    }

    /// Retain the flows the network just accepted at the outlets.
    pub(crate) fn record_step(&mut self, physics: &WaveguideNetwork) {
        for (k, &node) in self.nodes.iter().enumerate() {
            let q = physics.node_frame(node).map_or(0.0, |f| f.net_flow_into_node_m3_s);
            self.rings[k * self.len + self.head] = q;
        }
        self.head = (self.head + 1) % self.len;
    }
}

#[cfg(test)]
mod tests;
