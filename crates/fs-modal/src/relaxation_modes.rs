//! Frequency-dependent material damping with retained Maxwell memory.
//!
//! Each supplied branch has strain `e = b^T q`, stiffness `k_b > 0`, and
//! relaxation time `tau > 0`. Its stress obeys `s' = k_b e' - s/tau` and
//! opposes motion through `b s`. Together with the RELAXED stiffness K,
//! this realizes `M q'' + C q' + K q + sum(b s) = 0` exactly: its Laplace
//! damping is `C + sum(k_b tau b b^T / (1 + lambda tau))`.
//!
//! All `2*n + branches.len()` state poles are retained, including material
//! relaxation poles. Replacing the branches by a frequency-sampled damping
//! ratio would lose both this memory and the relaxed/unrelaxed stiffness
//! change. This dense spectrum is not a Rouché certificate or a time stepper.
//! Dense QR does not guarantee relative accuracy for very slow poles or tiny
//! real parts of oscillatory poles, especially across widely separated rates.
//! Exact source zero-count guards catch some lost poles, not every loss of
//! decay-rate accuracy. A free scalar coordinate with one branch instead
//! factors out its exact zero and uses stable quadratic root recovery.
//! A passive scalar resonator with one weak Maxwell arm instead solves the
//! small relaxation-pole shift directly, retaining weak mechanical decay
//! without subtracting nearly equal eigenvalues (see `weak_scalar_maxwell`).
//! Weak passive scalar multi-arm materials use positive physical loss sums for
//! the conjugate pair and source-rate brackets for every real material pole;
//! repeated-rate and uncoupled material states keep their full multiplicity.
//!
//! For symmetric K and positive-semidefinite symmetric part of C, the branch
//! storage `s^2/(2 k_b)` gives physical heat rate `s^2/(k_b tau)`; gyroscopic
//! entries in C do no work. The solver preserves the full source C, including
//! active entries, and does not certify the sign of its symmetric part.

use crate::{
    ModalError, quadratic_eigenvalues_with_cancel, quadratic_fraction, quadratic_mass_whiten,
    quadratic_power_scale, quadratic_scalar_eigenvalues,
};
use fs_la::{
    eigen_complex::{EigError, eig_with_cancel},
    factor::cholesky,
};
use fs_math::{c64::C64, dd::Dd, det};

mod scalar;
mod shapes;
pub use shapes::{
    MaxwellBranchDirection, MaxwellForceResidue, MaxwellMode, MaxwellModeError, MaxwellModeResidue,
    MaxwellPoleDirection, PoleOptions, PreparedMaxwellModes,
};

/// One physical scalar strain channel of a generalized Maxwell material.
/// Multiple channels may share a relaxation time; their mutual modal coupling
/// is retained through their original strain directions, without projection.
#[derive(Clone, Debug)]
pub struct MaxwellDampingBranch {
    /// Original generalized-coordinate strain direction b; exactly n entries.
    pub strain: Vec<f64>,
    /// Positive branch stiffness, conjugate to the supplied scalar strain.
    /// Omit zero-stiffness branches: they carry no material storage or stress.
    pub stiffness: f64,
    /// Positive physical relaxation time [s].
    pub relaxation_time_s: f64,
}

/// Compute the complete damped spectrum of a mechanical system with Maxwell
/// memory, using the default per-window QR budget of 120 sweeps.
///
/// Source M, C, and RELAXED K are dense row-major n-by-n matrices. M must be
/// symmetric positive definite. C retains every entry (including gyroscopic
/// terms); K may retain unstable stiffness. The branches add their unrelaxed
/// stiffness dynamically, so it must not also be included in K.
///
/// Poles use the convention `exp(lambda*t)` in physical [1/s] and are sorted
/// by real then imaginary part. Every state is included. Dense work is cubic
/// in the full state count; callers must bound their source model size. Dense
/// QR has no relative decay-rate certificate for weak or widely separated
/// losses; source-proven excessive zero counts refuse instead of publishing.
///
/// # Errors
/// Invalid shapes, nonfinite coefficients, nonpositive branch data, a failed
/// mass factor, unrepresentable scaling/poles, a source-inconsistent zero-pole
/// count, or the dense QR sweep limit.
pub fn maxwell_eigenvalues(
    m: &[f64],
    c: &[f64],
    k: &[f64],
    n: usize,
    branches: &[MaxwellDampingBranch],
) -> Result<Vec<C64>, ModalError> {
    maxwell_eigenvalues_with_cancel(m, c, k, n, branches, 120, || false)
}

/// Cancellable, explicitly sweep-bounded Maxwell spectrum. Cancellation is
/// checked during preparation and inside the shared dense eigensolver; no
/// partial spectrum accompanies a refusal.
///
/// # Errors
/// The same refusals as [`maxwell_eigenvalues`], a zero sweep budget, or
/// [`ModalError::Cancelled`].
#[allow(clippy::too_many_lines, clippy::float_cmp)] // One atomic realization; exact source symmetry/zero identities.
pub fn maxwell_eigenvalues_with_cancel(
    m: &[f64],
    c: &[f64],
    k: &[f64],
    n: usize,
    branches: &[MaxwellDampingBranch],
    max_sweeps_per_window: usize,
    mut cancelled: impl FnMut() -> bool,
) -> Result<Vec<C64>, ModalError> {
    let invalid = |field, index| ModalError::DampingInput { field, index };
    let scaling = || ModalError::QuadraticEig {
        detail: "Maxwell state scaling requires representable coefficients and poles".into(),
    };
    if cancelled() {
        return Err(ModalError::Cancelled);
    }
    let nn = n
        .checked_mul(n)
        .ok_or(invalid("Maxwell matrix extent", n))?;
    let dim = n
        .checked_mul(2)
        .and_then(|x| x.checked_add(branches.len()))
        .ok_or(invalid("Maxwell state extent", n))?;
    let extent = dim
        .checked_mul(dim)
        .ok_or(invalid("Maxwell state extent", n))?;
    if m.len() != nn || c.len() != nn || k.len() != nn || (n == 0 && !branches.is_empty()) {
        return Err(invalid("Maxwell matrix shape", n));
    }
    if max_sweeps_per_window == 0 {
        return Err(invalid("Maxwell QR sweep budget", 0));
    }
    for i in 0..n {
        if !m[i * n + i].is_finite() || m[i * n + i] <= 0.0 {
            return Err(ModalError::MassNotSpd { negative: 0 });
        }
        for j in 0..n {
            if !m[i * n + j].is_finite() || m[i * n + j] != m[j * n + i] {
                return Err(ModalError::MassNotSpd { negative: 0 });
            }
        }
    }
    if let Some(i) = c.iter().chain(k).position(|v| !v.is_finite()) {
        return Err(invalid("finite Maxwell mechanical coefficients", i));
    }
    for (index, branch) in branches.iter().enumerate() {
        if cancelled() {
            return Err(ModalError::Cancelled);
        }
        if branch.strain.len() != n || branch.strain.iter().any(|x| !x.is_finite()) {
            return Err(invalid("Maxwell strain direction", index));
        }
        if !branch.stiffness.is_finite() || branch.stiffness <= 0.0 {
            return Err(invalid("positive Maxwell branch stiffness", index));
        }
        if !branch.relaxation_time_s.is_finite() || branch.relaxation_time_s <= 0.0 {
            return Err(invalid("positive Maxwell relaxation time", index));
        }
    }
    // Preserve the quadratic owner's analytic scalar and diagonal paths when
    // there is no memory, rather than sending them through an augmented QR.
    if branches.is_empty() {
        return quadratic_eigenvalues_with_cancel(m, c, k, n, max_sweeps_per_window, cancelled);
    }
    // The exact free-coordinate root factors out before any dense deflation.
    // QR can otherwise turn a small but physical relaxation pole into a second
    // zero, either through its threshold or cancellation in its 2x2 formula.
    if n == 1 && k[0] == 0.0 && branches.len() == 1 {
        let poles = free_scalar_maxwell(m[0], c[0], &branches[0])?;
        if cancelled() {
            return Err(ModalError::Cancelled);
        }
        return Ok(poles);
    }
    if n == 1 && k[0] > 0.0 && c[0] >= 0.0 && branches.len() == 1 {
        if let Some(poles) = weak_scalar_maxwell(m[0], c[0], k[0], &branches[0])? {
            if cancelled() {
                return Err(ModalError::Cancelled);
            }
            return Ok(poles);
        }
    }
    if n == 1 && branches.len() > 1 {
        if let Some(poles) =
            scalar::weak_scalar_multi_maxwell(m[0], c[0], k[0], branches, &mut cancelled)?
        {
            return Ok(poles);
        }
    }
    let exponents: Vec<_> = (0..n)
        .map(|i| quadratic_fraction(m[i * n + i]).1.div_euclid(2))
        .collect();
    let mut time_exponent = i32::MIN;
    for i in 0..n {
        if cancelled() {
            return Err(ModalError::Cancelled);
        }
        for j in 0..n {
            let coordinate = exponents[i] + exponents[j];
            if c[i * n + j] != 0.0 {
                time_exponent = time_exponent.max(quadratic_fraction(c[i * n + j]).1 - coordinate);
            }
            if k[i * n + j] != 0.0 {
                time_exponent = time_exponent
                    .max((quadratic_fraction(k[i * n + j]).1 - coordinate + 1).div_euclid(2));
            }
        }
    }
    for branch in branches {
        time_exponent = time_exponent.max(-quadratic_fraction(branch.relaxation_time_s).1);
        for (i, &strain) in branch.strain.iter().enumerate() {
            if strain != 0.0 {
                time_exponent =
                    time_exponent.max(coupling_parts(branch.stiffness, strain).1 - exponents[i]);
            }
        }
    }
    let checked_scale = |value: f64, exponent| {
        let result = quadratic_power_scale(value, exponent);
        if !result.is_finite() || (value != 0.0 && result == 0.0) {
            Err(scaling())
        } else {
            Ok(result)
        }
    };
    let transform = |source: &[f64], order| -> Result<Vec<f64>, ModalError> {
        source
            .iter()
            .enumerate()
            .map(|(index, &value)| {
                checked_scale(
                    value,
                    -exponents[index / n] - exponents[index % n] - order * time_exponent,
                )
            })
            .collect()
    };
    let mass = transform(m, 0)?;
    let chol = cholesky(&mass, n).map_err(|_| ModalError::MassNotSpd { negative: 0 })?;
    let mut damping = transform(c, 1)?;
    let mut stiffness = transform(k, 2)?;
    quadratic_mass_whiten(&mut damping, n, &chol);
    if cancelled() {
        return Err(ModalError::Cancelled);
    }
    quadratic_mass_whiten(&mut stiffness, n, &chol);
    if damping.iter().chain(&stiffness).any(|x| !x.is_finite()) {
        return Err(scaling());
    }
    // Energy coordinates z=s/sqrt(k_b), and one shared time scale, yield
    // antisymmetric velocity/memory blocks [-g; g^T]. Keeping both signs and
    // using the same rounded g preserves the physical power exchange.
    let mut state = vec![C64::ZERO; extent];
    for i in 0..n {
        state[i * dim + n + i] = C64::ONE;
        for j in 0..n {
            state[(n + i) * dim + j] = C64::from_re(-stiffness[i * n + j]);
            state[(n + i) * dim + n + j] = C64::from_re(-damping[i * n + j]);
        }
    }
    for (branch_index, branch) in branches.iter().enumerate() {
        if cancelled() {
            return Err(ModalError::Cancelled);
        }
        let memory = 2 * n + branch_index;
        let (tau_fraction, tau_exponent) = quadratic_fraction(branch.relaxation_time_s);
        state[memory * dim + memory] = C64::from_re(-checked_scale(
            1.0 / tau_fraction,
            -tau_exponent - time_exponent,
        )?);
        let mut coupling = vec![0.0; n];
        for i in 0..n {
            if branch.strain[i] != 0.0 {
                let (fraction, exponent) = coupling_parts(branch.stiffness, branch.strain[i]);
                coupling[i] = checked_scale(fraction, exponent - exponents[i] - time_exponent)?;
            }
            for j in 0..i {
                coupling[i] = (-chol.l(i, j)).mul_add(coupling[j], coupling[i]);
            }
            coupling[i] /= chol.l(i, i);
            if !coupling[i].is_finite() {
                return Err(scaling());
            }
            state[(n + i) * dim + memory] = C64::from_re(-coupling[i]);
            state[memory * dim + n + i] = C64::from_re(coupling[i]);
        }
    }
    let mut poles =
        eig_with_cancel(&state, dim, max_sweeps_per_window, &mut cancelled).map_err(|error| {
            match error {
                EigError::Cancelled => ModalError::Cancelled,
                EigError::NotConverged(failure) => ModalError::QuadraticEig {
                    detail: format!("Maxwell: {failure:?}"),
                },
                EigError::Invalid { what } => ModalError::QuadraticEig {
                    detail: what.into(),
                },
            }
        })?;
    for pole in &mut poles {
        pole.re = checked_scale(pole.re, time_exponent)?;
        pole.im = checked_scale(pole.im, time_exponent)?;
    }
    // These are exact source identities, not numerical rank estimates. At
    // lambda=0 every finite-tau memory block is invertible, so triangular K
    // with nonzero diagonal permits no zero pole. For a free scalar, the
    // first nonzero coefficient is C + sum(k_b*tau*b²); its sign can prove
    // one zero without forming products that could overflow or underflow.
    let zero_limit = if n == 1 && k[0] == 0.0 {
        let coupled = branches.iter().any(|branch| branch.strain[0] != 0.0);
        if c[0] > 0.0 || (c[0] == 0.0 && coupled) || (c[0] != 0.0 && !coupled) {
            Some(1)
        } else if c[0] == 0.0 && !coupled {
            Some(2)
        } else {
            None // Active C can cancel the material contribution exactly.
        }
    } else if (0..n).all(|i| k[i * n + i] != 0.0)
        && ((0..n).all(|i| (0..i).all(|j| k[i * n + j] == 0.0))
            || (0..n).all(|i| (0..i).all(|j| k[j * n + i] == 0.0)))
    {
        Some(0)
    } else {
        None
    };
    if let Some(limit) = zero_limit {
        let zeros = poles.iter().filter(|pole| **pole == C64::ZERO).count();
        if zeros > limit {
            return Err(ModalError::QuadraticEig {
                detail: format!(
                    "Maxwell QR returned {zeros} exact zero poles; source permits at most {limit}"
                ),
            });
        }
    }
    poles.sort_by(|a, b| a.re.total_cmp(&b.re).then_with(|| a.im.total_cmp(&b.im)));
    if cancelled() {
        return Err(ModalError::Cancelled);
    }
    Ok(poles)
}

// A product/quotient of nonnegative finite source factors. Carry the binary
// exponent separately so a change of mechanical units cannot first overflow
// b², m*omega or the dimensional branch stiffness. There are at most four
// numerator and three denominator factors here, so mantissa products are safe.
#[allow(clippy::float_cmp)] // Exact zero factors carry exact uncoupled/source identities.
fn positive_ratio(numerator: &[f64], denominator: &[f64]) -> f64 {
    let mut value = 1.0;
    let mut exponent = 0;
    for &factor in numerator {
        if factor == 0.0 {
            return 0.0;
        }
        let (fraction, power) = quadratic_fraction(factor);
        value *= fraction;
        exponent += power;
    }
    for &factor in denominator {
        let (fraction, power) = quadratic_fraction(factor);
        value /= fraction;
        exponent -= power;
    }
    quadratic_power_scale(value, exponent)
}

// Let omega²=K/m, a=C/(m*omega), g=k_b*b²/K and r=1/tau. The real
// relaxation pole is -u=-r*(1-t), where
//   t = g / ((r/omega)²*(1-t)² - a*(r/omega)*(1-t) + 1 + g).
// For 0<=a<=1 and 0<=g<=1/8, this maps [0,1/7] into itself and has
// derivative magnitude <=4/9. Indeed D=x²-a*x+1>=3/4 and
// |x*(2*x-a)|/D<=8/3. The other two roots are therefore an underdamped
// conjugate pair: their total decay is C/m + r*t, without trace subtraction.
// Their product is omega²/(1-t); r*t/omega<=g gives an imaginary-part floor.
//
// These restrictions select a numerically stable route, not a material-law
// approximation. All other parameter regimes keep the complete dense owner.
// Work is bounded by 64 scalar contractions, independent of frequency/rate
// separation. The physical shift r*t is evaluated from original factors even
// when the dimensionless t or g is smaller than the least positive f64.
#[allow(clippy::float_cmp, clippy::too_many_lines)] // Complete source admission and bounded scalar pole recovery.
fn weak_scalar_maxwell(
    mass: f64,
    damping: f64,
    stiffness: f64,
    branch: &MaxwellDampingBranch,
) -> Result<Option<Vec<C64>>, ModalError> {
    let invalid = || ModalError::QuadraticEig {
        detail: "scalar Maxwell decay requires representable poles".into(),
    };
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
        return Err(invalid());
    }
    let strain = branch.strain[0].abs();
    let a = positive_ratio(&[damping], &[mass, omega]);
    let g = positive_ratio(&[branch.stiffness, strain, strain], &[stiffness]);
    if a > 1.0 || g > 0.125 {
        return Ok(None);
    }
    let rate = 1.0 / branch.relaxation_time_s;
    if !rate.is_finite() {
        return Err(invalid());
    }
    let fast = rate >= omega;
    let ratio = if fast { omega / rate } else { rate / omega };
    let denominator = |t: f64| {
        let remaining = 1.0 - t;
        if fast {
            remaining * (remaining - a * ratio) + (1.0 + g) * ratio * ratio
        } else {
            let x = ratio * remaining;
            x * (x - a) + 1.0 + g
        }
    };
    let mut t = 0.0;
    let mut converged = false;
    for _ in 0..64 {
        let next = if fast { g * ratio * ratio } else { g } / denominator(t);
        converged = next == t || (next - t).abs() <= 4.0 * f64::EPSILON * next;
        t = next;
        if converged {
            break;
        }
    }
    if !converged {
        return Ok(None);
    }
    let shift = if fast {
        positive_ratio(
            &[branch.stiffness, strain, strain, branch.relaxation_time_s],
            &[mass, denominator(t)],
        )
    } else {
        positive_ratio(
            &[branch.stiffness, strain, strain, rate],
            &[stiffness, denominator(t)],
        )
    };
    let viscous = positive_ratio(&[damping], &[mass]);
    let total_decay = viscous + shift;
    let decay = if total_decay.is_finite() {
        0.5 * total_decay
    } else {
        0.5 * viscous + 0.5 * shift
    };
    let scaled_decay = decay / omega;
    let frequency = omega * det::sqrt(1.0 / (1.0 - t) - scaled_decay * scaled_decay);
    let real = -rate * (1.0 - t);
    if ![decay, frequency, real].into_iter().all(f64::is_finite)
        || frequency <= 0.0
        || (strain != 0.0 && shift == 0.0)
        || ((strain != 0.0 || damping != 0.0) && decay == 0.0)
    {
        return Err(invalid());
    }
    let mut poles = vec![
        C64::from_re(real),
        C64::new(-decay, -frequency),
        C64::new(-decay, frequency),
    ];
    poles.sort_by(|a, b| a.re.total_cmp(&b.re).then_with(|| a.im.total_cmp(&b.im)));
    Ok(Some(poles))
}

// lambda * [m*tau*lambda² + (m+c*tau)*lambda + c+k_b*b²*tau].
// Preserve the slow root with the original quadratic product-of-roots owner.
// Factored source products avoid overflowing b² or losing m*tau before the
// common coefficient scale is chosen. DD accumulation retains signed-C
// cancellation before each final coefficient is rounded to f64.
#[allow(clippy::too_many_lines)] // Factoring, coefficient admission and stable root recovery share one source equation.
fn free_scalar_maxwell(
    mass: f64,
    damping: f64,
    branch: &MaxwellDampingBranch,
) -> Result<Vec<C64>, ModalError> {
    fn product(factors: &[f64]) -> (Dd, i32) {
        let mut value = Dd::ONE;
        let mut exponent = 0;
        for &factor in factors {
            if factor == 0.0 {
                return (Dd::ZERO, 0);
            }
            let (fraction, power) = quadratic_fraction(factor);
            value = value * Dd::from_f64(fraction.copysign(factor));
            exponent += power;
        }
        (value, exponent)
    }
    fn sum(left: (Dd, i32), right: (Dd, i32)) -> (f64, i32) {
        let exponent = if left.0 == Dd::ZERO {
            right.1
        } else if right.0 == Dd::ZERO {
            left.1
        } else {
            left.1.max(right.1)
        };
        let scaled = |(value, power): (Dd, i32)| {
            Dd::from_pair(
                quadratic_power_scale(value.hi, power - exponent),
                quadratic_power_scale(value.lo, power - exponent),
            )
        };
        let value = (scaled(left) + scaled(right)).to_f64();
        if value == 0.0 {
            (0.0, 0)
        } else {
            let (fraction, power) = quadratic_fraction(value);
            (fraction.copysign(value), power + exponent)
        }
    }
    let tau = branch.relaxation_time_s;
    let strain = branch.strain[0];
    let coefficients = [
        sum(product(&[mass, tau]), (Dd::ZERO, 0)),
        sum(product(&[mass]), product(&[damping, tau])),
        sum(
            product(&[damping]),
            product(&[branch.stiffness, strain, strain, tau]),
        ),
    ];
    let minimum = coefficients
        .iter()
        .filter(|x| x.0 != 0.0)
        .map(|x| x.1)
        .min()
        .unwrap();
    let maximum = coefficients
        .iter()
        .filter(|x| x.0 != 0.0)
        .map(|x| x.1)
        .max()
        .unwrap();
    let invalid = || ModalError::QuadraticEig {
        detail: "free Maxwell quadratic requires representable coefficients and poles".into(),
    };
    // Mantissas have magnitude in [1,2); keep the largest safely finite and
    // never silently erase a nonzero coefficient at the subnormal boundary.
    let lower = maximum - 1022;
    let upper = minimum + 1074;
    if lower > upper {
        return Err(invalid());
    }
    let shift = (minimum + maximum).div_euclid(2).clamp(lower, upper);
    let mut scaled = [0.0; 3];
    for (index, &(fraction, power)) in coefficients.iter().enumerate() {
        scaled[index] = quadratic_power_scale(fraction, power - shift);
        if !scaled[index].is_finite() || (fraction != 0.0 && scaled[index] == 0.0) {
            return Err(invalid());
        }
    }
    let mut poles = quadratic_scalar_eigenvalues(scaled[0], scaled[1], scaled[2])?;
    poles.push(C64::ZERO);
    poles.sort_by(|a, b| a.re.total_cmp(&b.re).then_with(|| a.im.total_cmp(&b.im)));
    Ok(poles)
}

// sqrt(k)*b as a small mantissa plus exponent, without overflowing the
// dimensional product before mass/time equilibration. Both inputs nonzero.
fn coupling_parts(stiffness: f64, strain: f64) -> (f64, i32) {
    let (k_fraction, k_exponent) = quadratic_fraction(stiffness);
    let (b_fraction, b_exponent) = quadratic_fraction(strain);
    let root = det::sqrt(
        k_fraction
            * if k_exponent.rem_euclid(2) == 0 {
                1.0
            } else {
                2.0
            },
    );
    (
        (root * b_fraction).copysign(strain),
        k_exponent.div_euclid(2) + b_exponent,
    )
}
