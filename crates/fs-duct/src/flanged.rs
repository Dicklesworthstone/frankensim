//! Flanged circular-pipe radiation: one multimodal model for every `ka`.
//!
//! W. E. Zorumski, "Generalized radiation impedances and reflection
//! coefficients of circular and annular ducts", J. Acoust. Soc. Am.
//! 54(6), 1667–1673 (1973). A pipe of radius `a` ends in an infinite
//! rigid baffle. The rigid-duct axisymmetric modes
//! `ψ_n(r) = J0(α_n r/a)/J0(α_n)` (`J1(α_n) = 0`, `α_0 = 0`) are
//! orthonormal under the aperture mean, and their Hankel transforms are
//! closed form: `∫_0^a ψ_n J0(κr) r dr = a² D_n(κa)`,
//! `D_n(x) = x J1(x)/(x² − α_n²)`. With `K = ka` the King spectral
//! integral gives the aperture radiation matrix (`e^{-iωt}`, so a
//! mass-like reactance has `Im Z < 0`)
//!
//! `Z_mn/(ρc) = 2K [∫_0^K D_m D_n(√(K²−s²)) ds − i ∫_0^∞ D_m D_n(√(K²+s²)) ds]`.
//!
//! Both integrands depend on `x²` only, so they are analytic in `s`, and
//! fixed nodes give a load that is smooth in `K`. A causal line needs
//! that: a jump in `Z(ω)` is a two-sided `1/n` tail in its impulse
//! response. Higher duct modes leave the aperture back into the pipe with
//! `p_n = −ζ_n v_n`, `ζ_n/(ρc) = K/√(K² − α_n²)` (`−iK/√(α_n² − K²)` below
//! cutoff). The plane-mode load is the Schur complement
//! `Z_eff = Z_00 − Z_0E (Z_EE + diag ζ)^{-1} Z_E0`.
//! Across a stationary gas interface, pressure and normal velocity are
//! continuous. The radiation matrix uses the exterior `ρ_e, c_e`, while
//! the eliminated duct modes retain the interior `ρ_i, c_i`:
//! `ζ_n/(ρ_e c_e) = (ρ_i c_i)/(ρ_e c_e) K_i/√(K_i² − α_n²)`.
//!
//! `Z_00` is the closed-form piston ([`fs_phs::baffled_piston_impedance`]),
//! so one mode IS that piston. The mode ladder moves the low-`ka` end
//! correction from the piston's `8/(3π) a` to the flanged-pipe
//! `0.8216 a` (Norris & Sheng, J. Sound Vib. 135(1), 85–93, 1989).
//!
//! The reactive integral splits at `x = TAIL_X`, which lies above every
//! admitted `α_n`. Below the split it uses composite Gauss–Legendre.
//! Above it, `J1²` is the order-one Hankel expansion
//! `(1/(πx))[(P² + Q²) + Re(W e^{2ix})]`, `W = −2PQ + i(P² − Q²)`. The
//! non-oscillatory part integrates in `u = TAIL_X/x`. The oscillatory
//! part is rotated onto `TAIL_X + iy`, where it decays like `e^{−2y}`;
//! no pole or branch point lies in `Re x > TAIL_X`.

use crate::DuctError;
use crate::levine_schwinger::gauss_legendre;
use fs_la::eigen_complex::lu_complex;
use fs_material::gas::GasState;
use fs_math::{c64::C64, det};
use std::sync::LazyLock;

/// Production mode count, chosen from the measured mode ladder
/// (`flanged_mode_ladder_converges_to_0p8216`).
pub const FLANGED_MODES: usize = 16;
/// Largest admitted mode count (the tail split must clear `α_{N−1}`).
pub const FLANGED_MAX_MODES: usize = 24;
/// Largest admitted `ka`. A limit of the fixed node design, not of the
/// model: Zorumski's multimodal load is exact for every `ka`. Measured by
/// `flanged_design_is_measured_across_and_beyond_the_cap`, the quadrature
/// `Z_00` against the closed-form piston is 3.5e-14 at `ka = 20`, at most
/// 6.5e-14 through 26, then 2.9e-13 at 28, 2.0e-11 at 30 (past the 1e-12
/// gate), 3.6e-8 at 34 and 1.5e-5 at 38. The cap sits well inside the
/// plateau, not on its edge; above it the load refuses by name.
pub const FLANGED_MAX_KA: f64 = 24.0;

/// Reactive split, above `α_23 ≈ 72.9` with room for the rotated tail.
const TAIL_X: f64 = 86.0;
// Node design: Gauss–Legendre resolves e^{iωx} once n > e·ω·h/4 with h
// the half-width (error ~(e ω h/(4n))^{2n}). The near panels span ~7.2 in
// `s` against J1²'s period π, the real interval at most K/π periods.
// The closed-form piston cross-check pins the design at 1e-12.
const REAL_NODES: usize = 32;
const NEAR_PANELS: usize = 12;
const PANEL_NODES: usize = 16;
const TAIL_AVG_NODES: usize = 16;
const TAIL_OSC_NODES: usize = 24;
/// `|x − α_n|` below which `D_n` uses its removable-singularity expansion.
const ROOT_GUARD: f64 = 1.0e-5;

struct Rules {
    real: Vec<(f64, f64)>,
    panel: Vec<(f64, f64)>,
    avg: Vec<(f64, f64)>,
    osc: Vec<(f64, f64)>,
    /// `(α_n, J0(α_n))` for `n = 1..FLANGED_MAX_MODES`.
    roots: Vec<(f64, f64)>,
}

static RULES: LazyLock<Result<Rules, DuctError>> = LazyLock::new(|| {
    Ok(Rules {
        real: gauss_legendre(REAL_NODES),
        panel: gauss_legendre(PANEL_NODES),
        avg: gauss_legendre(TAIL_AVG_NODES),
        osc: gauss_legendre(TAIL_OSC_NODES),
        roots: j1_roots_with_j0(FLANGED_MAX_MODES - 1)?,
    })
});

/// Order-zero Bessel `J0` for the root table: the shared double-double
/// series up to 18, the Hankel expansion (optimally truncated) above it.
fn bessel_j0(x: f64) -> f64 {
    if x <= 18.0 {
        return crate::modal::bessel_j0(x);
    }
    let (p, q) = hankel_pq_real(0.0, x);
    let (s, c) = det::sin_cos(x);
    det::sqrt(1.0 / (core::f64::consts::PI * x)) * ((p + q) * c + (p - q) * s)
}

/// Hankel `P, Q` of order `ν` at real `x > 0` (DLMF 10.17.3), optimally truncated.
fn hankel_pq_real(nu: f64, x: f64) -> (f64, f64) {
    let mu = 4.0 * nu * nu;
    let (mut p, mut q) = (0.0, 0.0);
    let mut term = 1.0_f64;
    let mut last = f64::INFINITY;
    for m in 0..128_u32 {
        if term.abs() > last {
            break;
        }
        last = term.abs();
        match m % 4 {
            0 => p += term,
            1 => q += term,
            2 => p -= term,
            _ => q -= term,
        }
        if term.abs() < 1e-18 {
            break;
        }
        let mf = f64::from(m);
        term *= (mu - (2.0 * mf + 1.0) * (2.0 * mf + 1.0)) / (8.0 * (mf + 1.0)) / x;
    }
    (p, q)
}

/// Order-one Hankel `P, Q` at complex `x` with `|x| ≥ TAIL_X`. Fixed term
/// count: the terms fall below 1e-18 well before the series turns.
fn hankel_pq_one(x: C64) -> (C64, C64) {
    let inv = x.recip();
    let (mut p, mut q) = (C64::ZERO, C64::ZERO);
    let mut term = C64::ONE;
    for m in 0..14_u32 {
        match m % 4 {
            0 => p = p + term,
            1 => q = q + term,
            2 => p = p - term,
            _ => q = q - term,
        }
        let mf = f64::from(m);
        term = term * inv.scale((4.0 - (2.0 * mf + 1.0) * (2.0 * mf + 1.0)) / (8.0 * (mf + 1.0)));
    }
    (p, q)
}

/// The first `count` positive `J1` roots with `J0` at each, by bisection
/// inside `((n + 1/4)π ± 1/2)`.
fn j1_roots_with_j0(count: usize) -> Result<Vec<(f64, f64)>, DuctError> {
    let mut roots = Vec::with_capacity(count);
    for n in 1..=count {
        let centre = (n as f64 + 0.25) * core::f64::consts::PI;
        let (mut a, mut b) = (centre - 0.5, centre + 0.5);
        let mut fa = fs_math::special::bessel_j1(a);
        if fa * fs_math::special::bessel_j1(b) > 0.0 {
            return Err(DuctError::BadParameter {
                what: "flanged mode root bracket lost its sign change",
            });
        }
        for _ in 0..200 {
            let mid = f64::midpoint(a, b);
            if mid <= a || mid >= b {
                break;
            }
            let fm = fs_math::special::bessel_j1(mid);
            if fa * fm <= 0.0 {
                b = mid;
            } else {
                a = mid;
                fa = fm;
            }
        }
        let alpha = f64::midpoint(a, b);
        roots.push((alpha, bessel_j0(alpha)));
    }
    Ok(roots)
}

/// `D_n(x) = x J1(x)/(x² − α_n²)` for every admitted mode at one real node.
fn fill_d(x: f64, j1: f64, roots: &[(f64, f64)], d: &mut [f64]) {
    d[0] = if x == 0.0 { 0.5 } else { j1 / x };
    for (slot, &(alpha, j0)) in d[1..].iter_mut().zip(roots) {
        let delta = x - alpha;
        *slot = if delta.abs() < ROOT_GUARD {
            // J1(α+δ) = J0(α)[δ − δ²/(2α) + O(δ³)].
            x / (x + alpha) * j0 * (1.0 - delta / (2.0 * alpha))
        } else {
            x * j1 / (delta * (x + alpha))
        };
    }
}

/// Accumulate `weight · D_m D_n` into the packed upper triangle.
fn accumulate(acc: &mut [f64], d: &[f64], weight: f64) {
    let n = d.len();
    let mut idx = 0;
    for m in 0..n {
        let wm = weight * d[m];
        for dn in &d[m..n] {
            acc[idx] += wm * dn;
            idx += 1;
        }
    }
}

/// Dimensionless aperture radiation matrix `Z_mn/(ρc)` (packed upper
/// triangle, `modes(modes+1)/2` entries) from the quadrature, including
/// `Z_00`. Production replaces `Z_00` by the closed form. This is the seam
/// the closed-form cross-check certifies.
pub(crate) fn radiation_matrix(ka: f64, modes: usize) -> Result<Vec<C64>, DuctError> {
    let rules = RULES.as_ref().map_err(Clone::clone)?;
    let roots = &rules.roots[..modes - 1];
    let pairs = modes * (modes + 1) / 2;
    let mut re = vec![0.0_f64; pairs];
    let mut im = vec![0.0_f64; pairs];
    let mut d = vec![0.0_f64; modes];
    // Real part: s in (0, K), x = sqrt((K - s)(K + s)).
    for &(t, w) in &rules.real {
        let s = ka * t;
        let x = det::sqrt((ka - s) * (ka + s));
        fill_d(x, fs_math::special::bessel_j1(x), roots, &mut d);
        accumulate(&mut re, &d, ka * w);
    }
    // Reactive near part: s in (0, s_c), x = sqrt(K² + s²) in (K, TAIL_X).
    let s_c = det::sqrt((TAIL_X - ka) * (TAIL_X + ka));
    let width = s_c / NEAR_PANELS as f64;
    for panel in 0..NEAR_PANELS {
        let base = panel as f64 * width;
        for &(t, w) in &rules.panel {
            let s = t.mul_add(width, base);
            let x = det::sqrt(ka.mul_add(ka, s * s));
            fill_d(x, fs_math::special::bessel_j1(x), roots, &mut d);
            accumulate(&mut im, &d, width * w);
        }
    }
    // Reactive tail, non-oscillatory part: x = TAIL_X/u.
    let mut h = vec![0.0_f64; modes];
    for &(u, w) in &rules.avg {
        let x = TAIL_X / u;
        let (p, q) = hankel_pq_real(1.0, x);
        h[0] = 1.0 / x;
        for (slot, &(alpha, _)) in h[1..].iter_mut().zip(roots) {
            *slot = x / ((x - alpha) * (x + alpha));
        }
        let root = det::sqrt((x - ka) * (x + ka));
        let weight = w * (TAIL_X / (u * u)) * p.mul_add(p, q * q) / (core::f64::consts::PI * root);
        accumulate(&mut im, &h, weight);
    }
    // Reactive tail, oscillatory part: Re ∫ G e^{2ix} dx over (TAIL_X, ∞)
    // = Re[i e^{2i TAIL_X} ∫_0^∞ G(TAIL_X + iy) e^{-2y} dy], y = v/(1 − v).
    let mut osc = vec![C64::ZERO; pairs];
    let mut hc = vec![C64::ZERO; modes];
    for &(v, w) in &rules.osc {
        let y = v / (1.0 - v);
        let jacobian = 1.0 / ((1.0 - v) * (1.0 - v));
        let x = C64::new(TAIL_X, y);
        let (p, q) = hankel_pq_one(x);
        let wave = C64::new(0.0, 1.0) * (p * p - q * q) - (p * q).scale(2.0);
        let root = (x * x - C64::from_re(ka * ka)).sqrt();
        let common = wave / root.scale(core::f64::consts::PI);
        let x2 = x * x;
        hc[0] = x.recip();
        for (slot, &(alpha, _)) in hc[1..].iter_mut().zip(roots) {
            *slot = x / (x2 - C64::from_re(alpha * alpha));
        }
        let weight = common.scale(w * jacobian * det::exp(-2.0 * y));
        let mut idx = 0;
        for m in 0..modes {
            let wm = weight * hc[m];
            for hn in &hc[m..modes] {
                osc[idx] = osc[idx] + wm * *hn;
                idx += 1;
            }
        }
    }
    let (sin_c, cos_c) = det::sin_cos(2.0 * TAIL_X);
    let rotate = C64::new(0.0, 1.0) * C64::new(cos_c, sin_c);
    let two_k = 2.0 * ka;
    Ok((0..pairs)
        .map(|i| {
            let x = im[i] + (rotate * osc[i]).re;
            C64::new(two_k * re[i], -two_k * x)
        })
        .collect())
}

fn packed(i: usize, j: usize, modes: usize) -> usize {
    let (m, n) = if i <= j { (i, j) } else { (j, i) };
    m * modes - m * (m + 1) / 2 + n
}

/// Complete flanged-aperture radiation impedance in the duct modal basis.
///
/// Returns the row-major `modes x modes` matrix `P = Z U` [Pa s/m³],
/// where `psi_n(r) = J0(alpha_n r/radius)/|J0(alpha_n)|`,
/// `P_n = area^-1 integral p psi_n dA`, and
/// `U_n = integral v_out psi_n dA`. Thus the peak-phasor radiated power is
/// `Re(sum_n conj(U_n) P_n)/2`. All mutual and higher-mode entries are retained;
/// no outgoing pipe modes are eliminated as in [`flanged_pipe_impedance`].
///
/// The existing spectral quadrature uses wall-positive modes `J0/J0(alpha_n)`.
/// Its signed congruence `diag(sign J0) Z diag(sign J0)` and the single
/// division by aperture area preserve both pressure/flow pairing and power.
/// Mode zero is the closed-form baffled piston. The exterior gas alone owns
/// this operator; a connected bore must separately admit a common static
/// pressure and propagate with its own gas. Recompute at every retained-mode
/// rung rather than reusing a condensed scalar termination.
///
/// # Errors
/// Invalid gas, nonpositive or nonfinite radius/frequency, a mode count outside
/// `1..=FLANGED_MAX_MODES`, exterior `ka > FLANGED_MAX_KA`, or an
/// unrepresentable impedance. This finite-basis infinite-baffle operator does
/// not establish arbitrary-geometry or out-of-band radiation accuracy.
pub fn flanged_modal_impedance(
    exterior: &GasState,
    radius: f64,
    omega: f64,
    modes: usize,
) -> Result<Vec<C64>, DuctError> {
    crate::validate_gas(exterior)?;
    if !(radius > 0.0 && radius.is_finite() && omega > 0.0 && omega.is_finite()) {
        return Err(DuctError::BadParameter {
            what: "flanged modal radiation needs positive finite radius and frequency",
        });
    }
    if !(1..=FLANGED_MAX_MODES).contains(&modes) {
        return Err(DuctError::BadParameter {
            what: "flanged radiation mode count outside 1..=FLANGED_MAX_MODES",
        });
    }
    let ka = omega / exterior.sound_speed * radius;
    if !(ka <= FLANGED_MAX_KA) {
        return Err(DuctError::RadiationKaTooLarge { ka });
    }
    let z0 = exterior.density * exterior.sound_speed
        / (core::f64::consts::PI * radius * radius);
    if !(z0.is_finite() && z0 > 0.0) {
        return Err(DuctError::Singular);
    }
    let (r00, x00) =
        fs_phs::baffled_piston_impedance(exterior.density, exterior.sound_speed, radius, omega)
            .map_err(|_| DuctError::BadParameter {
                what: "flanged radiation piston term",
            })?;
    let mut matrix = vec![C64::new(r00, x00); modes * modes];
    if modes > 1 {
        let packed_matrix = radiation_matrix(ka, modes)?;
        let rules = RULES.as_ref().map_err(Clone::clone)?;
        let sign = |n: usize| {
            if n == 0 || rules.roots[n - 1].1 > 0.0 { 1.0 } else { -1.0 }
        };
        for m in 0..modes {
            for n in m..modes {
                if m == 0 && n == 0 {
                    continue;
                }
                let value = packed_matrix[packed(m, n, modes)].scale(z0 * sign(m) * sign(n));
                matrix[m * modes + n] = value;
                matrix[n * modes + m] = value;
            }
        }
    }
    if matrix.iter().any(|z| !z.re.is_finite() || !z.im.is_finite()) {
        return Err(DuctError::Singular);
    }
    Ok(matrix)
}

/// Plane-mode flanged-pipe load `Z_L = p/U` [Pa s/m³] at the mouth of a
/// circular pipe of `radius` in an infinite rigid baffle, from `modes`
/// aperture modes. One analytic model for every admitted `ka`.
///
/// # Errors
/// [`DuctError::RadiationKaTooLarge`] above [`FLANGED_MAX_KA`];
/// [`DuctError::BadParameter`] for invalid inputs, a mode count outside
/// `1..=FLANGED_MAX_MODES`, a failed piston evaluation, a singular modal
/// solve, or a non-passive result.
pub fn flanged_pipe_impedance(
    state: &GasState,
    radius: f64,
    omega: f64,
    modes: usize,
) -> Result<C64, DuctError> {
    flanged_pipe_impedance_with_exterior(state, state, radius, omega, modes)
}

/// Flanged aperture between a local duct gas and a homogeneous exterior.
/// The exterior owns radiation; the duct gas owns the eliminated higher
/// modes. Both coherent states must have the same static pressure. This is
/// linear acoustics at a stationary interface, without mixing or mean flow.
/// Identical states reproduce [`flanged_pipe_impedance`] bit for bit.
///
/// # Errors
/// Invalid states, unequal static pressure, invalid geometry/frequency/mode
/// count, either gas's `ka` above [`FLANGED_MAX_KA`], or a singular/non-passive load.
pub fn flanged_pipe_impedance_with_exterior(
    state: &GasState,
    exterior: &GasState,
    radius: f64,
    omega: f64,
    modes: usize,
) -> Result<C64, DuctError> {
    crate::validate_gas(state)?;
    crate::validate_gas(exterior)?;
    if state.pressure != exterior.pressure {
        return Err(DuctError::BadParameter {
            what: "flanged radiation interface requires equal static pressure",
        });
    }
    if !(radius > 0.0 && radius.is_finite() && omega > 0.0 && omega.is_finite()) {
        return Err(DuctError::BadParameter {
            what: "flanged radiation needs positive finite radius and frequency",
        });
    }
    if !(1..=FLANGED_MAX_MODES).contains(&modes) {
        return Err(DuctError::BadParameter {
            what: "flanged radiation mode count outside 1..=FLANGED_MAX_MODES",
        });
    }
    let ka = omega / exterior.sound_speed * radius;
    if !(ka <= FLANGED_MAX_KA) {
        return Err(DuctError::RadiationKaTooLarge { ka });
    }
    let internal_ka = omega / state.sound_speed * radius;
    if !(internal_ka <= FLANGED_MAX_KA) {
        return Err(DuctError::RadiationKaTooLarge { ka: internal_ka });
    }
    let key = (
        radius.to_bits(),
        omega.to_bits(),
        state.density.to_bits(),
        state.sound_speed.to_bits(),
        exterior.density.to_bits(),
        exterior.sound_speed.to_bits(),
        modes,
    );
    if let Some(&hit) = memo().get(&key) {
        return Ok(hit);
    }
    let load = flanged_load_uncached(state, exterior, radius, omega, ka, modes)?;
    let mut cache = memo();
    if cache.len() >= MEMO_CAP {
        cache.clear();
    }
    cache.insert(key, load);
    Ok(load)
}

/// Exact-input memo. The load is a pure function of these bits, so a hit
/// returns the same bits a recomputation would. Banks and corpus sweeps
/// revisit one hole radius on one frequency grid across fingerings.
type MemoKey = (u64, u64, u64, u64, u64, u64, usize);
const MEMO_CAP: usize = 1 << 20;
static MEMO: LazyLock<std::sync::Mutex<std::collections::HashMap<MemoKey, C64>>> =
    LazyLock::new(|| std::sync::Mutex::new(std::collections::HashMap::new()));

fn memo() -> std::sync::MutexGuard<'static, std::collections::HashMap<MemoKey, C64>> {
    MEMO.lock().unwrap_or_else(std::sync::PoisonError::into_inner)
}

fn flanged_load_uncached(
    state: &GasState,
    exterior: &GasState,
    radius: f64,
    omega: f64,
    ka: f64,
    modes: usize,
) -> Result<C64, DuctError> {
    let z0 = exterior.density * exterior.sound_speed / (core::f64::consts::PI * radius * radius);
    let (r00, x00) =
        fs_phs::baffled_piston_impedance(exterior.density, exterior.sound_speed, radius, omega)
            .map_err(|_| DuctError::BadParameter {
                what: "flanged radiation piston term",
            })?;
    let piston = C64::new(r00, x00);
    let load = if modes == 1 {
        piston
    } else {
        let rules = RULES.as_ref().map_err(Clone::clone)?;
        let z = radiation_matrix(ka, modes)?;
        let internal_ka = omega / state.sound_speed * radius;
        let impedance_ratio = state.density * state.sound_speed
            / (exterior.density * exterior.sound_speed);
        // Modes exactly at cutoff carry ζ = ∞ and drop out (the continuous limit).
        let keep: Vec<usize> = (1..modes)
            .filter(|&n| {
                let alpha = rules.roots[n - 1].0;
                (internal_ka - alpha) * (internal_ka + alpha) != 0.0
            })
            .collect();
        let e = keep.len();
        let mut a = vec![C64::ZERO; e * e];
        let mut b = vec![C64::ZERO; e];
        for (row, &m) in keep.iter().enumerate() {
            for (col, &n) in keep.iter().enumerate() {
                a[row * e + col] = z[packed(m, n, modes)];
            }
            let alpha = rules.roots[m - 1].0;
            let gap = (internal_ka - alpha) * (internal_ka + alpha);
            let zeta = if gap > 0.0 {
                C64::from_re(internal_ka / det::sqrt(gap))
            } else {
                C64::new(0.0, -internal_ka / det::sqrt(-gap))
            }.scale(impedance_ratio);
            a[row * e + row] = a[row * e + row] + zeta;
            b[row] = z[packed(0, m, modes)];
        }
        let lu = lu_complex(&a, e).map_err(|_| DuctError::BadParameter {
            what: "flanged radiation modal solve is singular",
        })?;
        let mut v = b.clone();
        lu.solve(&mut v);
        let correction = b.iter().zip(&v).fold(C64::ZERO, |acc, (&bi, &vi)| acc + bi * vi);
        piston - correction.scale(z0)
    };
    if !(load.re.is_finite() && load.im.is_finite() && load.re >= -1.0e-12 * load.abs()) {
        return Err(DuctError::BadParameter {
            what: "flanged radiation load is not finite and passive",
        });
    }
    Ok(load)
}

#[cfg(test)]
mod tests;
