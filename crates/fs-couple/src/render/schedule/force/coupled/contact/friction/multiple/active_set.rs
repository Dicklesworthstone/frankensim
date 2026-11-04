//! Active-set rescue for the joint multi-port friction solve (frankensim-o7qbp).
//!
//! The simultaneous friction owner normally advances with a Gauss-Seidel sweep:
//! each port re-chooses its stick/slide branch against the other ports' current
//! held tractions, and the step commits when one full sweep is jointly
//! consistent. The sweep can OSCILLATE between candidate partitions near a
//! stick<->slip transition and never land a clean sweep, so it exhausts
//! `max_sweeps` and refuses a step that does have a solution (o7qbp: the 5-cell
//! compliant bow ribbon at the first Helmholtz transition; enumerating all 3^p
//! partitions at that step found exactly ONE consistent joint solution).
//!
//! More than one consistent partition can exist. Sample 451 of the 160 N stop
//! bank in `normal_support_press` has three. The walk below never consults the
//! ports' history: it starts at all-stuck and flips the worst violator, so the
//! partition it keeps is a function of the step alone. At that step it keeps
//! SS+, which continues the previous motion; the walk itself makes nothing of
//! that general.
//!
//! When the walk exhausts on a plain line of at most five ports, the plain
//! rescue has a last resort (`solve_active_set_with_hysteresis`, MagentaBadger
//! 16900/17076). It solves every one of the 3^p partitions, with each sliding
//! row written in its branch speed so the zero-slip kink is a boundary rather
//! than a fold. It keeps only partitions that pass the walk's own residual gate
//! and `consistent_at`, then takes the one nearest the previous step's
//! branches. That choice is the hysteresis rule of McIntyre, Schumacher &
//! Woodhouse (1983), not a tie-break: the contact keeps its state while that
//! state exists, and when it ceases to exist the contact goes to the state that
//! remains. A step the walk resolves never reaches this. A step with no
//! consistent partition, or with a tie the rule cannot break, keeps the walk's
//! refusal unchanged.
//!
//! This module solves a consistent partition directly. It runs ONLY after the
//! Gauss-Seidel sweep has failed, so every step the sweep already converges is
//! bit-identical. It fixes a partition, joint-Newtons the effective load line
//! for it, checks the Karush-Kuhn-Tucker friction conditions, flips the one
//! worst-violated port, and repeats under a finite pivot budget; it refuses on
//! budget rather than looping. It is deterministic: one thread, fixed order, a
//! plain partial-pivot elimination, no reductions whose order depends on layout.
//!
//! Effective load line (per port i): `slip_i = f[i] - sum_j a[i*p+j] * t[j]`,
//! where `a` already carries each port's series (bow-hair) compliance on its
//! diagonal and `f` already carries each port's stored-stretch free-slip shift.
//! Every row is on the bank's one line: the endpoint line (end slip speed times
//! dt, MagentaBadger 8060), or the displacement line for a bank with a series
//! drive, a shear layer or a bound bow residual on any port.
//!
//! A stop bank coupled to the ports (MagentaBadger 12431, Proposal A) makes that
//! line move with the tractions: `slip(t) = f - a t + T R(t)`, with each
//! contact-loaded port's normal load `N_i(t) = R_c(i)(t)`. `R(t)` is the stops'
//! own reaction at held tractions `t`, so the stops are condensed out and are
//! not walk cells. `R(t)` is unique, because `I + D Nn` is a P-matrix
//! (`D = diag(phi') >= 0`, `Nn` symmetric positive semidefinite). Nothing at run
//! time relies on that argument: every accepted state passes the stops' own
//! check. The walk itself is unchanged. Only the partition solve differs. It
//! uses the condensed line's derivative, `-a + T (I + D Nn)^-1 D Nt` for the
//! slips, and the rows of `(I + D Nn)^-1 D Nt` for contact-loaded normals. A
//! partition whose stops can't be resolved is named `stops-unresolved` in the
//! refusal's payload.

use super::*;

/// The walk's name for a partition whose stops refused to resolve at its held tractions.
pub(super) const STOPS_UNRESOLVED: &str = "stops-unresolved";
const NEWTON_NONCONVERGENCE: &str = "newton-nonconvergence";
/// A partition Newton that failed with nothing to measure the miss by.
const DIVERGED: (&str, f64) = (NEWTON_NONCONVERGENCE, f64::INFINITY);

/// One evaluation of a condensed line (a stop bank coupled to the ports): every
/// port's slip and normal load at tractions `t`, with their derivatives in `t`.
pub(super) struct LinePoint {
    pub slip: Vec<f64>,
    pub normal: Vec<f64>,
    /// d slip_i / d t_j, row-major p x p.
    pub dslip: Vec<f64>,
    /// d normal_i / d t_j, row-major p x p. A prescribed load's row is zero.
    pub dnormal: Vec<f64>,
}

/// A load line whose slips and normal loads move with the tractions.
pub(super) trait CondensedLine {
    /// The line at tractions `t`, or the stops' refusal to resolve there.
    fn at(&self, t: &[f64]) -> Result<LinePoint, ModalFrictionError>;
}

#[cfg(test)]
thread_local! {
    /// The raw number of the last stops' refusal the walk met on this thread: an
    /// unconverged root's signed residual, or a broken budget's value. Only the
    /// tests read it, to pin a refusal's size to the root that refused.
    pub(super) static LAST_STOPS_REFUSAL: core::cell::Cell<Option<f64>> = const { core::cell::Cell::new(None) };
}

#[cfg(test)]
fn record_refusal(point: &Result<LinePoint, ModalFrictionError>) {
    if let Err(error) = point {
        let raw = match error {
            ModalFrictionError::Network(ModalCouplingError::ContactSolve { residual_n, .. }) => Some(*residual_n),
            ModalFrictionError::Network(ModalCouplingError::Budget { value, .. }) => Some(*value),
            _ => None,
        };
        LAST_STOPS_REFUSAL.with(|last| last.set(raw));
    }
}

#[cfg(not(test))]
fn record_refusal(_: &Result<LinePoint, ModalFrictionError>) {}

/// Every read of a condensed line by the walk. Outside the tests it is `line.at(t)`.
fn condensed_at(line: &dyn CondensedLine, t: &[f64]) -> Result<LinePoint, ModalFrictionError> {
    let point = line.at(t);
    record_refusal(&point);
    point
}

/// The line a walk runs on: the plain fixed line, or a condensed one.
enum Line<'a> {
    Fixed { a: &'a [f64], f: &'a [f64], normal: &'a [f64] },
    Condensed(&'a dyn CondensedLine),
}

impl Line<'_> {
    /// Every port's slip and normal load at tractions `t`.
    fn state(&self, t: &[f64]) -> Result<(Vec<f64>, Vec<f64>), ModalFrictionError> {
        match self {
            Line::Fixed { a, f, normal } => Ok((slips(a, f, t), normal.to_vec())),
            Line::Condensed(line) => condensed_at(*line, t).map(|point| (point.slip, point.normal)),
        }
    }
    /// The tractions of partition `cells`, or the name of why it did not resolve
    /// and by how much (unbounded when nothing measures the miss).
    fn solve(&self, law: &[&FrictionLaw], dt: f64, cells: &[Branch]) -> Result<Vec<f64>, (&'static str, f64)> {
        match self {
            Line::Fixed { a, f, normal } => solve_partition(a, f, normal, law, dt, cells, 100).ok_or(DIVERGED),
            Line::Condensed(line) => solve_condensed(*line, law, dt, cells, 100),
        }
    }
}

/// How far a stops' refusal missed: the residual of an unconverged root, the
/// value of a broken budget, otherwise unbounded.
fn refusal_size(error: &ModalFrictionError) -> f64 {
    match error {
        ModalFrictionError::Network(ModalCouplingError::ContactSolve { residual_n, .. }) => residual_n.abs(),
        ModalFrictionError::Network(ModalCouplingError::Budget { value, .. }) => *value,
        _ => f64::INFINITY,
    }
}

/// The worst Karush-Kuhn-Tucker violation of a candidate active set: which port,
/// which friction condition it broke, and by how much. Carried in a refusal so a
/// reader can see why the active-set rescue could not land a partition.
#[derive(Clone, Copy, Debug)]
pub(super) struct KktViolation {
    pub port: usize,
    pub condition: &'static str,
    pub size: f64,
}

/// Result of the active-set rescue: either a consistent partition, or an
/// exhaustion carrying the pivots spent, the last active set tried (as a compact
/// S/+/- string) and the worst KKT violation.
pub(super) enum Resolution {
    Solved(Vec<f64>, Vec<ModalFrictionRegime>),
    Exhausted {
        pivots: usize,
        last_set: String,
        worst: KktViolation,
    },
}

/// Compact S (stuck) / + (slide positive) / - (slide negative) encoding of a set.
fn encode(cells: &[Branch]) -> String {
    cells
        .iter()
        .map(|b| match b {
            Branch::Stuck => 'S',
            Branch::SlidePos => '+',
            Branch::SlideNeg => '-',
        })
        .collect()
}

/// One port's assigned branch in a candidate active set.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub(super) enum Branch {
    /// Zero relative slip; the held traction must stay inside the static cone.
    Stuck,
    /// Sliding with positive relative slip; traction is the positive kinetic value.
    SlidePos,
    /// Sliding with negative relative slip; traction is the negative kinetic value.
    SlideNeg,
}

/// The realized relative slip of every port under tractions `t`.
fn slips(a: &[f64], f: &[f64], t: &[f64]) -> Vec<f64> {
    let p = f.len();
    (0..p)
        .map(|i| {
            let mut s = f[i];
            for j in 0..p {
                s -= a[i * p + j] * t[j];
            }
            s
        })
        .collect()
}

/// Solve the dense `n x n` system `m x = b` by partial-pivot Gaussian
/// elimination. Deterministic; returns None on a singular pivot. The condensed
/// stops (`ObstacleSet::condensed`) use it for `(I + D Nn)^-1 D Nt`.
pub(super) fn dense_solve(m: &[f64], b: &[f64], n: usize) -> Option<Vec<f64>> {
    let mut aug = vec![0.0f64; n * (n + 1)];
    for i in 0..n {
        for j in 0..n {
            aug[i * (n + 1) + j] = m[i * n + j];
        }
        aug[i * (n + 1) + n] = b[i];
    }
    for c in 0..n {
        let mut piv = c;
        let mut best = aug[c * (n + 1) + c].abs();
        for r in (c + 1)..n {
            let v = aug[r * (n + 1) + c].abs();
            if v > best {
                best = v;
                piv = r;
            }
        }
        if best < 1e-300 {
            return None;
        }
        if piv != c {
            for k in c..=n {
                aug.swap(c * (n + 1) + k, piv * (n + 1) + k);
            }
        }
        let d = aug[c * (n + 1) + c];
        for r in 0..n {
            if r == c {
                continue;
            }
            let factor = aug[r * (n + 1) + c] / d;
            if factor == 0.0 {
                continue;
            }
            for k in c..=n {
                aug[r * (n + 1) + k] -= factor * aug[c * (n + 1) + k];
            }
        }
    }
    let mut x = vec![0.0f64; n];
    for i in 0..n {
        let d = aug[i * (n + 1) + i];
        if d == 0.0 {
            return None;
        }
        x[i] = aug[i * (n + 1) + n] / d;
    }
    Some(x)
}

/// The kinetic traction a sliding port carries at relative slip `slip`: magnitude
/// `mu_k(|slip|/dt) * normal`, signed by `sign`, the direction the PARTITION
/// assigns the port (+1 positive slip, -1 negative), so the friction force is
/// dissipative when the assignment is consistent.
///
/// The sign is the assignment, not the sign of the current iterate's slip. The
/// traction jumps by `2 * mu * normal` across `slip = 0`, and a consistent
/// partition can sit within ~1e-10 m of it (a bow cell barely sliding): signing
/// by the iterate made the full Newton step cycle between the two sides of the
/// jump forever, and the one consistent partition was refused.
fn kinetic_traction(law: &FrictionLaw, normal: f64, sign: f64, slip: f64, dt: f64) -> Result<f64, ModalFrictionError> {
    let speed = (slip / dt).abs();
    let mu = kinetic_mu(law, speed)?;
    Ok(sign * mu * normal)
}

/// The direction a partition assigns a sliding port: +1 for positive slip, -1 for negative.
fn branch_sign(b: Branch) -> f64 {
    if b == Branch::SlideNeg { -1.0 } else { 1.0 }
}

/// Joint-Newton the tractions for a FIXED partition `cells`. Stuck ports pin
/// `slip_i = 0`; sliding ports pin `t_i = kinetic_traction(branch sign, |slip_i|)`.
/// Returns the tractions, or None if the Newton does not resolve.
///
/// A stuck row `slip_i = f_i - sum_j a_ij t_j` is linear in `t`, so its Jacobian
/// row is exactly `-a_i`. It is not differenced: with step `h`, the forward
/// difference of a row of magnitude `f_i` carries an absolute error of about
/// `eps * f_i / h`, which for a near-rank-deficient `a` (many rigid-ribbon cells
/// on few string modes) reached 0.7 of the smallest eigenvalue of `a`, and Newton
/// stalled above the 1e-12 m gate. Sliding rows use the shared constitutive
/// derivative too. One value/slope evaluation per sliding row replaces the
/// traction-sized differences, which can cross a friction knee or round to zero.
fn solve_partition(
    a: &[f64],
    f: &[f64],
    normal: &[f64],
    law: &[&FrictionLaw],
    dt: f64,
    cells: &[Branch],
    max_newton: usize,
) -> Option<Vec<f64>> {
    let p = f.len();
    let mut t = vec![0.0f64; p];
    let evaluate = |t: &[f64]| -> Option<(Vec<f64>, Vec<f64>)> {
        let slip = slips(a, f, t);
        let mut r = vec![0.0f64; p];
        let mut jac = vec![0.0f64; p * p];
        for i in 0..p {
            match cells[i] {
                Branch::Stuck => {
                    r[i] = slip[i];
                    for j in 0..p {
                        jac[i * p + j] = -a[i * p + j];
                    }
                }
                Branch::SlidePos | Branch::SlideNeg => {
                    let sign = branch_sign(cells[i]);
                    let speed = (slip[i] / dt).abs();
                    let (mu, slope) = law[i].kinetic_coefficient_and_slope(speed).ok()?;
                    r[i] = t[i] - sign * mu * normal[i];
                    let along = slip_direction(slip[i]) / dt;
                    for j in 0..p {
                        let identity = if i == j { 1.0 } else { 0.0 };
                        jac[i * p + j] = identity
                            + sign * slope * along * a[i * p + j] * normal[i];
                    }
                }
            }
        }
        Some((r, jac))
    };
    for _ in 0..max_newton {
        let (r, jac) = evaluate(&t)?;
        if r.iter().all(|x| x.abs() <= 1e-18) {
            break;
        }
        let neg: Vec<f64> = r.iter().map(|x| -x).collect();
        let delta = dense_solve(&jac, &neg, p)?;
        for i in 0..p {
            t[i] += delta[i];
            if !t[i].is_finite() {
                return None;
            }
        }
    }
    let (r, _) = evaluate(&t)?;
    if r.iter().all(|x| x.abs() <= 1e-12) {
        Some(t)
    } else {
        None
    }
}

/// A law's kinetic coefficient at slip speed `speed`, exactly as `kinetic_traction` takes it.
fn kinetic_mu(law: &FrictionLaw, speed: f64) -> Result<f64, ModalFrictionError> {
    Ok(if speed == 0.0 { law.zero_speed_kinetic_coefficient()? } else { law.kinetic_coefficient(speed)? })
}

/// A sliding iterate at the cusp uses the symmetric generalized derivative of
/// `abs(slip)`. The partition still owns the traction sign and final mode check.
fn slip_direction(slip: f64) -> f64 {
    if slip > 0.0 { 1.0 } else if slip < 0.0 { -1.0 } else { 0.0 }
}

/// Joint-Newton the tractions of a FIXED partition on a condensed line. The rows
/// are those of `solve_partition`: a stuck port pins `slip_i(t) = 0`, a sliding
/// port pins `t_i = kinetic_traction(branch sign, N_i(t), |slip_i(t)|)`.
///
/// No row is differenced through the stops. Their reactions are converged only
/// to the stops' tolerance, so a difference over a small traction step would be
/// mostly that tolerance. Every row instead takes the line's own derivative. A
/// sliding row adds the original law's analytic slope in slip speed:
/// `d/dt_j (sign mu N_i) = sign (mu' (sgn slip_i / dt) dslip_ij N_i + mu dN_ij)`.
/// Returns the tractions, or why the partition did not resolve and by how much:
/// a stops' refusal is sized by `refusal_size` (the residual of an unconverged
/// root, the value of a broken budget), a Newton failure is unbounded.
fn solve_condensed(
    line: &dyn CondensedLine,
    law: &[&FrictionLaw],
    dt: f64,
    cells: &[Branch],
    max_newton: usize,
) -> Result<Vec<f64>, (&'static str, f64)> {
    let p = cells.len();
    let evaluate = |t: &[f64]| -> Result<(Vec<f64>, Vec<f64>), (&'static str, f64)> {
        let point = condensed_at(line, t).map_err(|error| (STOPS_UNRESOLVED, refusal_size(&error)))?;
        let mut r = vec![0.0f64; p];
        let mut jac = vec![0.0f64; p * p];
        for i in 0..p {
            if cells[i] == Branch::Stuck {
                r[i] = point.slip[i];
                jac[i * p..(i + 1) * p].copy_from_slice(&point.dslip[i * p..(i + 1) * p]);
                continue;
            }
            let sign = branch_sign(cells[i]);
            let speed = (point.slip[i] / dt).abs();
            let (mu, slope) = law[i].kinetic_coefficient_and_slope(speed).map_err(|_| DIVERGED)?;
            r[i] = t[i] - sign * mu * point.normal[i];
            // d|slip|/d slip; a port at exactly zero slip contributes no slope.
            let along = slip_direction(point.slip[i]) / dt;
            for j in 0..p {
                let d = slope * along * point.dslip[i * p + j] * point.normal[i] + mu * point.dnormal[i * p + j];
                let identity = if i == j { 1.0 } else { 0.0 };
                jac[i * p + j] = identity - sign * d;
            }
        }
        Ok((r, jac))
    };
    let mut t = vec![0.0f64; p];
    for _ in 0..max_newton {
        let (r, jac) = evaluate(&t)?;
        if r.iter().all(|x| x.abs() <= 1e-18) {
            break;
        }
        let neg: Vec<f64> = r.iter().map(|x| -x).collect();
        let delta = dense_solve(&jac, &neg, p).ok_or(DIVERGED)?;
        for i in 0..p {
            t[i] += delta[i];
            if !t[i].is_finite() {
                return Err(DIVERGED);
            }
        }
    }
    let (r, _) = evaluate(&t)?;
    if r.iter().all(|x| x.abs() <= 1e-12) {
        Ok(t)
    } else {
        Err(DIVERGED)
    }
}

/// Branch -> owner regime.
fn regime_of(b: Branch) -> ModalFrictionRegime {
    match b {
        Branch::Stuck => ModalFrictionRegime::Sticking,
        Branch::SlidePos | Branch::SlideNeg => ModalFrictionRegime::Sliding,
    }
}

/// Is partition `cells` with tractions `t` a consistent (KKT) joint solution?
/// Stuck ports stay inside the static cone; sliding ports actually slip in their
/// assigned direction and the traction opposes that slip (dissipative). The walk
/// uses `consistent_at` on its line's own evaluation; the tests use this form.
#[cfg(test)]
fn is_consistent(
    a: &[f64],
    f: &[f64],
    normal: &[f64],
    law: &[&FrictionLaw],
    dt: f64,
    cells: &[Branch],
    t: &[f64],
) -> Result<bool, ModalFrictionError> {
    consistent_at(&slips(a, f, t), normal, law, dt, cells, t)
}

/// `is_consistent` on slips and normal loads already evaluated at `t`.
fn consistent_at(
    slip: &[f64],
    normal: &[f64],
    law: &[&FrictionLaw],
    dt: f64,
    cells: &[Branch],
    t: &[f64],
) -> Result<bool, ModalFrictionError> {
    let p = slip.len();
    for i in 0..p {
        match cells[i] {
            Branch::Stuck => {
                let cap = law[i].static_coefficient()? * normal[i];
                if t[i].abs() > cap * (1.0 + 1e-9) {
                    return Ok(false);
                }
            }
            Branch::SlidePos => {
                if !(slip[i] > 0.0) || (slip[i] / dt).abs() < 1e-12 || t[i] <= 0.0 {
                    return Ok(false);
                }
            }
            Branch::SlideNeg => {
                if !(slip[i] < 0.0) || (slip[i] / dt).abs() < 1e-12 || t[i] >= 0.0 {
                    return Ok(false);
                }
            }
        }
    }
    Ok(true)
}

/// Active-set rescue. Returns `(tractions, regimes)` for the first consistent
/// partition it reaches, or None when the pivot budget is spent without consistency.
///
/// `seed` is the sweep's last (oscillating) regimes. It is unused: the walk
/// always starts at all-stuck.
#[allow(clippy::too_many_arguments)]
pub(super) fn solve_active_set(
    a: &[f64],
    f: &[f64],
    normal: &[f64],
    law: &[&FrictionLaw],
    dt: f64,
    seed: &[ModalFrictionRegime],
    budget: usize,
) -> Result<Resolution, ModalFrictionError> {
    let _ = seed;
    walk_line(&Line::Fixed { a, f, normal }, f.len(), law, dt, budget)
}

/// The same walk on a condensed line: a stop bank coupled to the ports.
pub(super) fn solve_active_set_condensed(
    line: &dyn CondensedLine,
    ports: usize,
    law: &[&FrictionLaw],
    dt: f64,
    budget: usize,
) -> Result<Resolution, ModalFrictionError> {
    walk_line(&Line::Condensed(line), ports, law, dt, budget)
}

/// Largest port count the last-resort enumeration runs on: 3^5 = 243 partitions.
pub(super) const ENUMERATION_MAX_PORTS: usize = 5;

/// The plain rescue: the walk, then its last resort. Only when the walk exhausts
/// on at most [`ENUMERATION_MAX_PORTS`] ports does [`nearest_consistent_partition`]
/// run. A step the walk resolves returns the walk's tractions bit for bit. A step
/// with no consistent partition, or with a tie the hysteresis rule cannot break,
/// returns the walk's own exhaustion, payload and all.
///
/// `previous` is each port's branch on the previous accepted step: a sliding
/// port in the direction of its traction, and any other port stuck. Before the
/// first accepted step, the bow is at rest and every port counts as stuck.
#[allow(clippy::too_many_arguments)]
pub(super) fn solve_active_set_with_hysteresis(
    a: &[f64],
    f: &[f64],
    normal: &[f64],
    law: &[&FrictionLaw],
    dt: f64,
    seed: &[ModalFrictionRegime],
    previous: &[Branch],
    budget: usize,
) -> Result<Resolution, ModalFrictionError> {
    let walked = solve_active_set(a, f, normal, law, dt, seed, budget)?;
    if matches!(walked, Resolution::Solved(..)) || f.len() > ENUMERATION_MAX_PORTS {
        return Ok(walked);
    }
    Ok(match nearest_consistent_partition(a, f, normal, law, dt, previous)? {
        Some((t, regimes)) => Resolution::Solved(t, regimes),
        None => walked,
    })
}

/// Every one of the 3^p partitions of a plain line, each solved once by
/// [`solve_partition_by_speed`] and put to the walk's own acceptance: the
/// residual on the actual slip within 1e-12, then [`consistent_at`]. Returns the
/// consistent partition that differs from `previous` in the fewest ports. It
/// returns `None` when no partition is consistent, or when two or more tie for
/// nearest. It never returns a partition that failed either test.
///
/// The selection is the hysteresis rule of McIntyre, Schumacher & Woodhouse
/// (1983, "On the oscillations of musical instruments", J. Acoust. Soc. Am.
/// 74:1325), discussed in Woodhouse (2014, "The acoustics of the violin: a
/// review", Rep. Prog. Phys. 77:115901) §2.2. Where the falling friction curve
/// cuts the load line more than once, the contact keeps its current state while
/// that state exists. When the state ceases to exist, the contact jumps to the
/// state that remains. Among several consistent joint states, the one nearest
/// the previous step's is that rule. The selection is physics, so a tie it
/// cannot break refuses instead of picking one.
///
/// Each partition is solved from one start, so a second root on the same
/// partition is not looked for.
pub(super) fn nearest_consistent_partition(
    a: &[f64],
    f: &[f64],
    normal: &[f64],
    law: &[&FrictionLaw],
    dt: f64,
    previous: &[Branch],
) -> Result<Option<(Vec<f64>, Vec<ModalFrictionRegime>)>, ModalFrictionError> {
    let p = f.len();
    if p > ENUMERATION_MAX_PORTS || previous.len() != p {
        return Ok(None);
    }
    let branches = [Branch::Stuck, Branch::SlidePos, Branch::SlideNeg];
    let count = (0..p).fold(1usize, |n, _| n * 3);
    let mut nearest: Option<(usize, Vec<f64>, Vec<Branch>)> = None;
    let mut tied = false;
    for code in 0..count {
        let mut rest = code;
        let cells: Vec<Branch> = (0..p)
            .map(|_| {
                let b = branches[rest % 3];
                rest /= 3;
                b
            })
            .collect();
        let Some(t) = solve_partition_by_speed(a, f, normal, law, dt, &cells, 100) else {
            continue;
        };
        if !consistent_at(&slips(a, f, &t), normal, law, dt, &cells, &t)? {
            continue;
        }
        let distance = cells.iter().zip(previous).filter(|(c, b)| c != b).count();
        let best = nearest.as_ref().map(|(d, ..)| *d);
        if best == Some(distance) {
            tied = true;
        } else if best.is_none_or(|d| distance < d) {
            nearest = Some((distance, t, cells));
            tied = false;
        }
    }
    Ok(match nearest {
        Some((_, t, cells)) if !tied => Some((t, cells.iter().map(|&b| regime_of(b)).collect())),
        _ => None,
    })
}

/// The slip speed past which a sliding row's own residual is monotone. It comes
/// from the law's own slope and the port's diagonal compliance `a_ii`, for any
/// `FrictionLaw`.
///
/// On a sliding row the residual is `t_i - sign mu(|slip_i| / dt) N_i`. With the
/// other tractions held, its derivative in `t_i` on the branch's side is
/// `1 + mu'(v) N_i a_ii / dt`. Where the law falls faster than `-dt / (N_i a_ii)`
/// the row is not monotone: it has a minimum, and a local maximum at the
/// zero-slip kink. On a step where the stuck state has just ceased to exist,
/// that maximum sits just below zero: it is the ghost of the near intersection.
/// A Newton started on the near side cycles on it. That happened on the violin
/// E5 onset step with the measured rosin row (WG2025 eq. 13, where
/// `-mu' N a_ii / dt` is 1.72 at rest; the hyperbolic 0.8/0.4/0.2 control gives
/// 0.078).
///
/// The existing single-port load-line owner locates the last stationary point
/// using the law's actual speed scale. Twice that speed seeds the rising side,
/// avoiding the zero Jacobian at the stationary point itself. A clamped linear
/// weakening law uses its exact zero-coefficient knee; flat and rising laws use
/// zero. There is no fixed speed grid that can miss a narrow weakening region.
///
/// No claim. `a_ii` is the port's own compliance, not the line condensed onto
/// the port by the partition's other rows. For a symmetric positive definite `a`,
/// condensing stuck rows out lowers the effective compliance, so the condensed
/// critical speed is at most this one, and the start lies past it. That is an
/// argument, not a proof, and it says nothing about sliding neighbours.
///
/// It is checked on five captured steps. On the violin E5 eq. 13 onset step it
/// finds the one consistent state. On this module's o7qbp, stop-bank and A4
/// three- and five-cell steps it finds the same consistent partitions as the
/// walk's own solve.
fn critical_speed(law: &FrictionLaw, normal: f64, a_ii: f64, dt: f64) -> Result<f64, ModalFrictionError> {
    law.validate()?;
    if !normal.is_finite() || normal < 0.0 || !a_ii.is_finite() || a_ii < 0.0
        || !dt.is_finite() || dt <= 0.0 {
        return Err(fs_tribo::TriboError::InvalidInput { field: "friction_load_line" }.into());
    }
    if normal == 0.0 || a_ii == 0.0 {
        return Ok(0.0);
    }
    let mobility = a_ii / dt;
    let gain = normal * mobility;
    if !mobility.is_finite() || !gain.is_finite() {
        return Err(fs_tribo::TriboError::InvalidInput { field: "friction_load_line_scale" }.into());
    }
    let turning = match *law {
        FrictionLaw::Coulomb { .. } => 0.0,
        FrictionLaw::VelocityDependent { mu_zero, slope_per_speed, .. } => {
            if slope_per_speed < 0.0 && -slope_per_speed * gain >= 1.0 {
                mu_zero / -slope_per_speed
            } else {
                0.0
            }
        }
        FrictionLaw::Stribeck { .. } | FrictionLaw::Hyperbolic { .. }
        | FrictionLaw::ShiftedPower { .. } => {
            let line = crate::stribeck_friction::StribeckLoadLine::new(law, normal, mobility)
                .map_err(|error| match error {
                    crate::stribeck_friction::ImplicitFrictionError::Constitutive(error) => error,
                    _ => fs_tribo::TriboError::InvalidInput { field: "friction_load_line_turning_speed" },
                })?;
            let (_, last) = line.critical_slips_m_s();
            if last == f64::INFINITY { 0.0 } else { last }
        }
    };
    let seed = 2.0 * turning;
    if !seed.is_finite() || seed < 0.0 {
        return Err(fs_tribo::TriboError::InvalidInput { field: "friction_load_line_turning_speed" }.into());
    }
    Ok(seed)
}

/// Joint-Newton a FIXED partition of a plain line, with every sliding row written
/// in its branch speed `s_i >= 0` instead of its traction. The row's traction
/// `t_i = sign mu(s_i) N_i` is explicit, and the row is
/// `slip_i(t) - sign s_i dt = 0`. No `|slip|` is taken, so the zero-slip kink is
/// the boundary `s_i = 0`, not a fold in the residual. Stuck rows are
/// `slip_i(t) = 0` in `t_i`, as in `solve_partition`.
///
/// Starts: stuck tractions at 0, and each sliding speed at the larger of its
/// free-slip speed `|f_i| / dt` and [`critical_speed`]. No step takes a speed
/// below a tenth of its current value, so a row stays on its branch. The
/// Jacobian uses the exact line (`-a_ij dt_j / dx_j`) and the original law's
/// analytic slope, without differencing nearly equal friction coefficients.
///
/// Iteration stops when no unknown moves, or after `max_newton` steps. The
/// result then faces `solve_partition`'s own gate: the residual on the actual
/// slip within 1e-12 on every row.
fn solve_partition_by_speed(
    a: &[f64],
    f: &[f64],
    normal: &[f64],
    law: &[&FrictionLaw],
    dt: f64,
    cells: &[Branch],
    max_newton: usize,
) -> Option<Vec<f64>> {
    let p = f.len();
    let sliding = |i: usize| cells[i] != Branch::Stuck;
    let mut x = vec![0.0f64; p];
    for i in 0..p {
        if sliding(i) {
            x[i] = (f[i] / dt).abs().max(critical_speed(law[i], normal[i], a[i * p + i], dt).ok()?);
        }
    }
    // The tractions and their derivatives in the unknowns.
    let tractions = |x: &[f64]| -> Option<(Vec<f64>, Vec<f64>)> {
        let mut t = vec![0.0f64; p];
        let mut dt_dx = vec![1.0f64; p];
        for i in 0..p {
            if sliding(i) {
                let sign = branch_sign(cells[i]);
                let (mu, slope) = law[i].kinetic_coefficient_and_slope(x[i]).ok()?;
                t[i] = sign * mu * normal[i];
                dt_dx[i] = sign * slope * normal[i];
            } else {
                t[i] = x[i];
            }
        }
        Some((t, dt_dx))
    };
    for _ in 0..max_newton {
        let (t, dt_dx) = tractions(&x)?;
        let slip = slips(a, f, &t);
        let mut neg = vec![0.0f64; p];
        let mut jac = vec![0.0f64; p * p];
        for i in 0..p {
            let row = if sliding(i) { slip[i] - branch_sign(cells[i]) * x[i] * dt } else { slip[i] };
            neg[i] = -row;
            for j in 0..p {
                jac[i * p + j] = -a[i * p + j] * dt_dx[j];
            }
            if sliding(i) {
                jac[i * p + i] -= branch_sign(cells[i]) * dt;
            }
        }
        let delta = dense_solve(&jac, &neg, p)?;
        let mut moved = false;
        for i in 0..p {
            let next = x[i] + delta[i];
            if !next.is_finite() {
                return None;
            }
            let next = if sliding(i) { next.max(0.1 * x[i]) } else { next };
            moved |= next != x[i];
            x[i] = next;
        }
        if !moved {
            break;
        }
    }
    let (t, _) = tractions(&x)?;
    let slip = slips(a, f, &t);
    for i in 0..p {
        let r = if sliding(i) {
            t[i] - kinetic_traction(law[i], normal[i], branch_sign(cells[i]), slip[i], dt).ok()?
        } else {
            slip[i]
        };
        if !(r.abs() <= 1e-12) {
            return None;
        }
    }
    Some(t)
}

fn walk_line(line: &Line<'_>, p: usize, law: &[&FrictionLaw], dt: f64, budget: usize)
    -> Result<Resolution, ModalFrictionError>
{
    // First partition: every port stuck. A stuck port that cannot hold will be
    // flipped to the sliding direction its over-cone traction indicates.
    let mut cells = vec![Branch::Stuck; p];
    let mut seen: Vec<Vec<Branch>> = Vec::new();
    let mut pivots = 0usize;
    // The worst KKT violation seen on the most recent partition, carried into an
    // exhaustion so the refusal can name port, condition and size.
    let mut worst = KktViolation { port: 0, condition: "none", size: 0.0 };
    for _ in 0..budget {
        pivots += 1;
        if seen.contains(&cells) {
            // Revisited a partition: the single-flip walk is cycling; stop.
            return Ok(Resolution::Exhausted { pivots, last_set: encode(&cells), worst });
        }
        seen.push(cells.clone());
        let t = match line.solve(law, dt, &cells) {
            Ok(t) => t,
            Err((condition, size)) => {
                // This partition did not resolve; perturb by flipping the port
                // with the largest |slip| that is pinned stuck, else give up.
                let slip = match line.state(&vec![0.0; p]) {
                    Ok((slip, _)) => slip,
                    Err(error) => {
                        worst = KktViolation { port: 0, condition: STOPS_UNRESOLVED, size: refusal_size(&error) };
                        return Ok(Resolution::Exhausted { pivots, last_set: encode(&cells), worst });
                    }
                };
                let mut flip = None;
                let mut w = 0.0;
                for i in 0..p {
                    if cells[i] == Branch::Stuck && slip[i].abs() > w {
                        w = slip[i].abs();
                        flip = Some(i);
                    }
                }
                match flip {
                    Some(i) => {
                        worst = KktViolation { port: i, condition, size: w };
                        cells[i] = if slip[i] >= 0.0 { Branch::SlidePos } else { Branch::SlideNeg };
                        continue;
                    }
                    None => {
                        worst = KktViolation { port: 0, condition, size };
                        return Ok(Resolution::Exhausted { pivots, last_set: encode(&cells), worst });
                    }
                }
            }
        };
        let (slip, normal) = match line.state(&t) {
            Ok(state) => state,
            Err(error) => {
                worst = KktViolation { port: 0, condition: STOPS_UNRESOLVED, size: refusal_size(&error) };
                return Ok(Resolution::Exhausted { pivots, last_set: encode(&cells), worst });
            }
        };
        if consistent_at(&slip, &normal, law, dt, &cells, &t)? {
            let regimes = cells.iter().map(|&b| regime_of(b)).collect();
            return Ok(Resolution::Solved(t, regimes));
        }
        // Flip the single worst KKT violator, recording port/condition/size.
        let mut flip = None;
        let mut wsize = 0.0f64;
        let mut wport = 0usize;
        let mut wcond = "none";
        for i in 0..p {
            let (viol, cond, want) = match cells[i] {
                Branch::Stuck => {
                    let cap = law[i].static_coefficient()? * normal[i];
                    (
                        t[i].abs() - cap,
                        "stuck-exceeds-static-cone",
                        if t[i] >= 0.0 { Branch::SlidePos } else { Branch::SlideNeg },
                    )
                }
                Branch::SlidePos => {
                    let bad = if slip[i] <= 0.0 || t[i] <= 0.0 { slip[i].abs().max(1.0) } else { -1.0 };
                    (bad, "slide-pos-wrong-direction-or-nondissipative", Branch::Stuck)
                }
                Branch::SlideNeg => {
                    let bad = if slip[i] >= 0.0 || t[i] >= 0.0 { slip[i].abs().max(1.0) } else { -1.0 };
                    (bad, "slide-neg-wrong-direction-or-nondissipative", Branch::Stuck)
                }
            };
            if viol > wsize {
                wsize = viol;
                wport = i;
                wcond = cond;
                flip = Some((i, want));
            }
        }
        worst = KktViolation { port: wport, condition: wcond, size: wsize };
        match flip {
            Some((i, want)) => cells[i] = want,
            None => return Ok(Resolution::Exhausted { pivots, last_set: encode(&cells), worst }),
        }
    }
    Ok(Resolution::Exhausted { pivots, last_set: encode(&cells), worst })
}

#[cfg(test)]
mod tests {
    use super::*;

    // Captured step-106 operating point of the 5-cell compliant bow ribbon
    // (frankensim-o7qbp): the effective load line the Gauss-Seidel sweep could
    // not converge. Enumerating all 3^5 partitions here must find exactly ONE
    // consistent joint solution, and the active-set solve must land it.
    const FREE_SLIP: [f64; 5] = [
        2.13421513728807844e-7,
        4.21264052786064505e-7,
        3.17417008614100996e-7,
        1.09274888317735950e-7,
        4.97456614456249488e-9,
    ];
    const COMPLIANCE: [f64; 25] = [
        5.64306391073245562e-6, 5.63755801463918332e-6, 5.64108199600910941e-6, 5.64350046616687824e-6, 5.64238912765379553e-6,
        5.63755801463918248e-6, 5.64153713689234851e-6, 5.64031646319536489e-6, 5.63326154555429124e-6, 5.62742755872163133e-6,
        5.64108199600910856e-6, 5.64031646319536235e-6, 5.64146937748549164e-6, 5.63915254790404680e-6, 5.63568001574791498e-6,
        5.64350046616687909e-6, 5.63326154555429039e-6, 5.63915254790404680e-6, 5.64630049046423115e-6, 5.64754857299451424e-6,
        5.64238912765379553e-6, 5.62742755872163133e-6, 5.63568001574791498e-6, 5.64754857299451424e-6, 5.65115279603956879e-6,
    ];
    const SERIES_COMPLIANCE: f64 = 7.25332907672419810e-5;
    const FREE_SLIP_SHIFT: [f64; 5] = [
        2.71280494184456747e-7,
        2.58005920731509035e-6,
        7.10037885325489961e-7,
        -2.46886722478804812e-8,
        -3.16847756698896144e-7,
    ];
    const NORMAL: [f64; 5] = [0.03, 0.03, 0.03, 0.03, 0.03];
    const DT: f64 = 2.08333333333333326e-5;

    fn law() -> FrictionLaw {
        FrictionLaw::Stribeck {
            static_mu: 1.2,
            kinetic_mu: 0.35,
            characteristic_speed: 0.05,
            viscous_per_speed: 0.0,
        }
    }

    /// G1: a fixed 1e-11 N finite difference vanishes next to these otherwise
    /// ordinary finite tractions. The constant law has an exact unit Jacobian.
    #[test]
    fn large_contact_loads_keep_the_exact_coulomb_partition() {
        let law = FrictionLaw::Coulomb { static_mu: 0.8, kinetic_mu: 0.4 };
        let a = [1e-8, 2e-9, 2e-9, 1e-8];
        let f = [0.02, -0.04];
        let n = [1e6, 2e6];
        let cells = [Branch::SlidePos, Branch::SlideNeg];
        let t = solve_partition(&a, &f, &n, &[&law; 2], 0.01, &cells, 2)
            .expect("finite Coulomb loads do not require resolvable traction perturbations");
        assert_eq!(t, [400_000.0, -800_000.0]);
        assert!(is_consistent(&a, &f, &n, &[&law; 2], 0.01, &cells, &t).unwrap());
    }

    /// G1: the derivative retains both actual slip response and changing load.
    #[test]
    fn rate_and_normal_load_jacobians_recover_independent_roots() {
        let law = FrictionLaw::VelocityDependent {
            static_mu: 0.8, mu_zero: 0.2, slope_per_speed: 0.4,
        };
        let cells = [Branch::SlidePos];
        let plain = solve_partition(&[0.2], &[0.5], &[2.0], &[&law], 0.1, &cells, 4)
            .expect("affine velocity-dependent friction");
        assert!((plain[0] - 4.4 / 2.6).abs() < 1e-13);

        struct ChangingLoad;
        impl CondensedLine for ChangingLoad {
            fn at(&self, t: &[f64]) -> Result<LinePoint, ModalFrictionError> {
                Ok(LinePoint {
                    slip: vec![0.5 - 0.2 * t[0]],
                    normal: vec![1.0 + 0.2 * t[0]],
                    dslip: vec![-0.2],
                    dnormal: vec![0.2],
                })
            }
        }
        // F = (1 + 0.2 F) * (0.2 + 0.4 * (0.5 - 0.2 F) / 0.1).
        // Thus 0.16 F² + 1.36 F - 2.2 = 0, independently of Newton.
        // det-ok: test-only independent quadratic-root oracle.
        let expected = 4.4 / (1.36 + (1.36_f64.powi(2) + 4.0 * 0.16 * 2.2).sqrt());
        let condensed = solve_condensed(&ChangingLoad, &[&law], 0.1, &cells, 32)
            .expect("joint normal-load and slip derivative");
        assert!((condensed[0] - expected).abs() < 1e-13);
    }

    /// G1/G3: weakening folds follow the material speed scale, including speeds
    /// far below the former 1e-6 m/s scan and above its fixed upper limit.
    #[test]
    fn sliding_seed_follows_the_actual_friction_speed_scale() {
        for c in [1e-200, 1e-12, 1e-6, 1.0, 1e6, 1e200] {
            let hyperbolic = FrictionLaw::Hyperbolic {
                static_mu: 0.8, kinetic_mu: 0.4, characteristic_speed: c,
            };
            let seed = critical_speed(&hyperbolic, 1.0, 9.0 * c / 0.4, 1.0).unwrap();
            assert!((seed / c - 4.0).abs() < 1e-13);

            let stribeck = FrictionLaw::Stribeck {
                static_mu: 0.8, kinetic_mu: 0.4, characteristic_speed: c,
                viscous_per_speed: 0.0,
            };
            let gain = 10.0 * c / 0.4;
            let seed = critical_speed(&stribeck, 1.0, gain, 1.0).unwrap();
            let (_, at_fold) = stribeck.kinetic_coefficient_and_slope(seed / 2.0).unwrap();
            let (_, after_fold) = stribeck.kinetic_coefficient_and_slope(seed).unwrap();
            assert!((1.0 + gain * at_fold).abs() < 2e-14);
            assert!(1.0 + gain * after_fold > 0.0);
        }
        let clamped = FrictionLaw::VelocityDependent {
            static_mu: 0.8, mu_zero: 0.4, slope_per_speed: -2.0,
        };
        assert_eq!(critical_speed(&clamped, 1.0, 1.0, 1.0).unwrap(), 0.4);
        assert_eq!(critical_speed(&clamped, 1.0, 0.1, 1.0).unwrap(), 0.0);
        assert_eq!(critical_speed(&clamped, 0.0, 1.0, 1.0).unwrap(), 0.0);
    }

    fn operating_point() -> ([f64; 25], [f64; 5]) {
        let mut a = COMPLIANCE;
        for i in 0..5 {
            a[i * 5 + i] += SERIES_COMPLIANCE;
        }
        let mut f = [0.0; 5];
        for i in 0..5 {
            f[i] = FREE_SLIP[i] + FREE_SLIP_SHIFT[i];
        }
        (a, f)
    }

    fn all_branches() -> [Branch; 3] {
        [Branch::Stuck, Branch::SlidePos, Branch::SlideNeg]
    }

    #[test]
    fn enumerating_the_243_partitions_finds_exactly_one_consistent_solution() {
        let (a, f) = operating_point();
        let l = law();
        let laws: Vec<&FrictionLaw> = (0..5).map(|_| &l).collect();
        let bs = all_branches();
        let mut consistent = Vec::new();
        for i0 in bs {
            for i1 in bs {
                for i2 in bs {
                    for i3 in bs {
                        for i4 in bs {
                            let cells = [i0, i1, i2, i3, i4];
                            if let Some(t) =
                                solve_partition(&a, &f, &NORMAL, &laws, DT, &cells, 200)
                            {
                                if is_consistent(&a, &f, &NORMAL, &laws, DT, &cells, &t).unwrap() {
                                    consistent.push((cells, t));
                                }
                            }
                        }
                    }
                }
            }
        }
        assert_eq!(
            consistent.len(),
            1,
            "expected a unique consistent joint solution at the o7qbp step-106 state, got {}: {:?}",
            consistent.len(),
            consistent.iter().map(|(c, _)| c).collect::<Vec<_>>()
        );
        // The unique partition is port 0 stuck, port 1 sliding +, ports 2-4 stuck.
        let (cells, _t) = &consistent[0];
        assert_eq!(
            *cells,
            [Branch::Stuck, Branch::SlidePos, Branch::Stuck, Branch::Stuck, Branch::Stuck]
        );
    }

    // Sample 451 of the 160 N stop bank in normal_support_press, dumped where the
    // sweep refused. The stops have no coupling to these ports, so these rows are
    // the whole friction step. The sweep alternates SSS and SS+ here.
    const STOP_DT: f64 = 1.0 / 48_000.0;
    const STOP_LINE: [f64; 9] = [
        6.5438179355632094e-6, 1.9929157425147624e-6, 1.9363011049747153e-6,
        1.9929157425147624e-6, 6.598357687970561e-6, 1.957313263922187e-6,
        1.9363011049747153e-6, 1.957313263922187e-6, 6.485641401497002e-6,
    ];
    const STOP_FREE_SLIP: [f64; 3] = [-1.6153907305051958e-8, -3.306641346093483e-8, 1.3310924282490778e-8];
    const STOP_NORMAL: [f64; 3] = [0.007177012064098477, 0.008213326310517676, 0.005221382740318865];
    const STOP_SS_PLUS: [f64; 3] = [-0.00151066, -0.00517458, 0.00208855];

    #[test]
    fn a_stop_bank_step_has_three_consistent_partitions_and_the_walk_lands_on_ss_plus() {
        use Branch::{SlideNeg, SlidePos, Stuck};
        let law = FrictionLaw::Coulomb { static_mu: 0.8, kinetic_mu: 0.4 };
        let laws = [&law, &law, &law];
        let (a, f, n) = (&STOP_LINE, &STOP_FREE_SLIP, &STOP_NORMAL);
        let cap = |i: usize| 0.8 * n[i];
        // The force a port needs to stick with the others held: the sweep's own test.
        let required = |t: &[f64], i: usize| t[i] + slips(a, f, t)[i] / a[i * 3 + i];
        let assert_near = |got: &[f64], want: &[f64]| {
            for (g, w) in got.iter().zip(want) {
                assert!((g - w).abs() < 1e-8, "force {g} against {w}");
            }
        };

        // Port forces [N], whether the partition is consistent, and the ports the
        // sweep can't leave on their branch.
        let rows: [([Branch; 3], [f64; 3], bool, &[usize]); 4] = [
            ([Stuck, Stuck, Stuck], [-0.00202881, -0.00569704, 0.00437739], false, &[2]),
            ([Stuck, Stuck, SlidePos], STOP_SS_PLUS, true, &[2]),
            ([Stuck, SlideNeg, Stuck], [-0.00259823, -0.00328533, 0.00381956], true, &[1]),
            ([Stuck, SlideNeg, SlidePos], [-0.00208603, -0.00328533, 0.00208855], true, &[1, 2]),
        ];
        for (cells, forces, consistent, moved) in &rows {
            let t = solve_partition(a, f, n, &laws, STOP_DT, cells, 100).expect("every row solves");
            assert_near(&t, forces);
            assert_eq!(is_consistent(a, f, n, &laws, STOP_DT, cells, &t).unwrap(), *consistent, "{cells:?}");
            // A stuck port past its cone slides; a sliding port inside it re-sticks.
            // No partition keeps every port where the sweep's rule puts it.
            for &i in *moved {
                assert_eq!(required(&t, i).abs() <= cap(i), cells[i] != Stuck, "{cells:?} port {i}");
            }
        }
        let all_stuck = solve_partition(a, f, n, &laws, STOP_DT, &[Stuck; 3], 100).unwrap();
        assert!((all_stuck[2] / cap(2) - 1.048).abs() < 5e-4, "SSS is 4.8 % over its cap");

        let bs = all_branches();
        let mut consistent = Vec::new();
        for code in 0..27usize {
            let cells: Vec<Branch> = (0..3u32).map(|i| bs[code / 3usize.pow(i) % 3]).collect();
            if let Some(t) = solve_partition(a, f, n, &laws, STOP_DT, &cells, 100) {
                if is_consistent(a, f, n, &laws, STOP_DT, &cells, &t).unwrap() {
                    consistent.push(cells);
                }
            }
        }
        assert_eq!(consistent.len(), 3, "{consistent:?}");
        for (cells, _, is_one, _) in &rows {
            assert_eq!(consistent.iter().any(|c| c == cells), *is_one, "{cells:?}");
        }

        let seed = [ModalFrictionRegime::Sticking; 3];
        match solve_active_set(a, f, n, &laws, STOP_DT, &seed, 4 * 3 + 16).unwrap() {
            Resolution::Solved(t, regimes) => {
                assert_near(&t, &STOP_SS_PLUS);
                assert_eq!(regimes, [ModalFrictionRegime::Sticking, ModalFrictionRegime::Sticking,
                    ModalFrictionRegime::Sliding]);
            }
            Resolution::Exhausted { last_set, .. } => panic!("the walk stopped at {last_set}"),
        }
    }

    #[test]
    fn active_set_solve_lands_the_unique_partition() {
        let (a, f) = operating_point();
        let l = law();
        let laws: Vec<&FrictionLaw> = (0..5).map(|_| &l).collect();
        let seed = vec![ModalFrictionRegime::Sticking; 5];
        let (t, regimes) = match solve_active_set(&a, &f, &NORMAL, &laws, DT, &seed, 64).unwrap() {
            Resolution::Solved(t, regimes) => (t, regimes),
            Resolution::Exhausted { pivots, last_set, worst } => panic!(
                "active-set rescue found no consistent solution: {pivots} pivots, last set {last_set}, worst {worst:?}"
            ),
        };
        assert_eq!(
            regimes,
            vec![
                ModalFrictionRegime::Sticking,
                ModalFrictionRegime::Sliding,
                ModalFrictionRegime::Sticking,
                ModalFrictionRegime::Sticking,
                ModalFrictionRegime::Sticking,
            ]
        );
        // Port 1 slides within the kinetic band; the stuck ports stay in the cone.
        let cap = l.static_coefficient().unwrap() * NORMAL[0];
        for (i, &ti) in t.iter().enumerate() {
            if regimes[i] == ModalFrictionRegime::Sticking {
                assert!(ti.abs() <= cap * (1.0 + 1e-9), "port {i} stuck traction {ti} exceeds the cone {cap}");
            }
        }
    }

    #[test]
    fn active_set_refuses_when_the_pivot_budget_is_spent() {
        let (a, f) = operating_point();
        let l = law();
        let laws: Vec<&FrictionLaw> = (0..5).map(|_| &l).collect();
        let seed = vec![ModalFrictionRegime::Sticking; 5];
        // From the all-stuck start the unique partition S+SSS needs more than one
        // pivot (flip port 1, then re-solve), so budget 1 cannot resolve it and
        // the rescue refuses with a payload; the real default budget resolves it.
        match solve_active_set(&a, &f, &NORMAL, &laws, DT, &seed, 1).unwrap() {
            Resolution::Exhausted { pivots, last_set, worst } => {
                assert!(pivots <= 1, "budget 1 must spend at most one pivot, spent {pivots}");
                assert!(!last_set.is_empty(), "the refusal must carry the last active set");
                assert!(worst.condition != "none", "the refusal must carry a KKT violation");
            }
            Resolution::Solved(..) => panic!("budget 1 must not resolve the partition"),
        }
        assert!(
            matches!(
                solve_active_set(&a, &f, &NORMAL, &laws, DT, &seed, 64).unwrap(),
                Resolution::Solved(..)
            ),
            "the default budget must resolve the partition"
        );
    }

    // Captured refusal states of the finite-width violin bow on the A4 card, with 3
    // and 5 rigid-ribbon cells on 25 bending modes (material_strings at 91048d5c8,
    // 2026-10-30): the compliance matrix, free slip and normal force at the step where
    // the Gauss-Seidel sweep stalled (at 5 cells its iteration matrix has spectral
    // radius 1 - 1.5e-7, and the matrix has condition number 6e7) and this rescue
    // then refused. Each state has exactly ONE consistent partition. The reference
    // tractions are that solution, found over all 3^p partitions at 50 digits with an
    // analytic-Jacobian Newton. These states belong to a ribbon that is under-resolved
    // (25 modes against 2 mm cells); they pin the per-partition solve, not the sound.
    const A4_3_COMPLIANCE: [f64; 9] = [
        3.242189682783081e-6, 3.0060271872024254e-6, 3.071974318659879e-6,
        3.0060271872024258e-6, 3.1408969890100562e-6, 2.536988762713617e-6,
        3.0719743186598804e-6, 2.536988762713616e-6, 3.2461790651167124e-6,
    ];
    const A4_3_FREE_SLIP: [f64; 3] = [7.718593353040716e-8, 6.022912021836881e-8, 8.809302586191596e-8];
    const A4_3_NORMAL: [f64; 3] = [8.88e-2, 8.88e-2, 8.88e-2];
    const A4_3_REFERENCE: [f64; 3] = [-7.1039995226896366e-2, 2.9784357894763439e-2, 7.1039998232043456e-2];
    const A4_5_COMPLIANCE: [f64; 25] = [
        3.242189682783081e-6, 2.9144763779873244e-6, 3.147149255586746e-6, 3.1885709159069045e-6, 2.9911340259655605e-6,
        2.914476377987324e-6, 3.121797588234734e-6, 3.082077947364439e-6, 2.6346157436085015e-6, 2.2690673081597915e-6,
        3.147149255586745e-6, 3.082077947364439e-6, 3.1845775820084684e-6, 2.972113669697493e-6, 2.6760888082238244e-6,
        3.1885709159069045e-6, 2.6346157436085015e-6, 2.972113669697493e-6, 3.2611797107557466e-6, 3.1825532487603138e-6,
        2.9911340259655605e-6, 2.269067308159791e-6, 2.6760888082238235e-6, 3.1825532487603138e-6, 3.23169229485642e-6,
    ];
    const A4_5_FREE_SLIP: [f64; 5] = [1.8601190476190475e-8; 5];
    const A4_5_NORMAL: [f64; 5] = [5.328000000000001e-2; 5];
    const A4_5_REFERENCE: [f64; 5] = [
        -1.0664531455004529e-2,
        4.1891461592371590e-2,
        -4.2623999997376573e-2,
        3.6760181756948487e-3,
        1.7889132963283344e-2,
    ];
    const A4_DT: f64 = 1.0 / 144_000.0;
    const A4_TOL_N: f64 = 1e-12;

    fn a4_law() -> FrictionLaw {
        FrictionLaw::Stribeck {
            static_mu: 0.8,
            kinetic_mu: 0.4,
            characteristic_speed: 0.1,
            viscous_per_speed: 0.0,
        }
    }

    /// Every partition that solves AND is consistent, with its tractions.
    fn consistent_partitions(a: &[f64], f: &[f64], normal: &[f64]) -> Vec<(Vec<Branch>, Vec<f64>)> {
        let p = f.len();
        let l = a4_law();
        let laws: Vec<&FrictionLaw> = (0..p).map(|_| &l).collect();
        let bs = all_branches();
        let mut found = Vec::new();
        for code in 0..3usize.pow(p as u32) {
            let cells: Vec<Branch> = (0..p)
                .map(|i| bs[code / 3usize.pow(i as u32) % 3])
                .collect();
            if let Some(t) = solve_partition(a, f, normal, &laws, A4_DT, &cells, 100) {
                if is_consistent(a, f, normal, &laws, A4_DT, &cells, &t).unwrap() {
                    found.push((cells, t));
                }
            }
        }
        found
    }

    fn assert_matches_reference(t: &[f64], reference: &[f64]) {
        for (i, (&got, &want)) in t.iter().zip(reference).enumerate() {
            assert!(
                (got - want).abs() <= A4_TOL_N,
                "port {i}: traction {got} differs from the 50-digit reference {want} by more than {A4_TOL_N}"
            );
        }
    }

    fn walk(a: &[f64], f: &[f64], normal: &[f64], budget: usize) -> Resolution {
        let l = a4_law();
        let laws: Vec<&FrictionLaw> = (0..f.len()).map(|_| &l).collect();
        let seed = vec![ModalFrictionRegime::Sticking; f.len()];
        solve_active_set(a, f, normal, &laws, A4_DT, &seed, budget).unwrap()
    }

    #[test]
    fn a_sliding_traction_takes_the_partitions_direction_not_the_iterates_slip() {
        let l = a4_law();
        // 1e-10 m of slip is 7e-5 m/s, far below the Stribeck speed, so mu is static
        // to 1e-6 and the traction jump across slip = 0 is ~2 mu N = 0.14 N.
        let positive = kinetic_traction(&l, 0.0888, 1.0, -1e-10, A4_DT).unwrap();
        let negative = kinetic_traction(&l, 0.0888, -1.0, 1e-10, A4_DT).unwrap();
        assert!(
            positive > 0.07 && negative < -0.07,
            "the assigned sign must win: {positive}, {negative}"
        );
        assert_eq!(positive, -negative, "the magnitude depends on |slip| only");
    }

    #[test]
    fn a4_three_cell_step_has_one_consistent_partition_and_it_solves() {
        // The sliding partition sits ~2.5e-10 m from slip = 0: a Newton signed by the
        // iterate's slip cycled across the jump and never solved it.
        let found = consistent_partitions(&A4_3_COMPLIANCE, &A4_3_FREE_SLIP, &A4_3_NORMAL);
        assert_eq!(
            found.len(),
            1,
            "expected a unique consistent partition, got {:?}",
            found.iter().map(|(c, _)| c).collect::<Vec<_>>()
        );
        let (cells, t) = &found[0];
        assert_eq!(*cells, [Branch::SlideNeg, Branch::Stuck, Branch::SlidePos]);
        assert_matches_reference(t, &A4_3_REFERENCE);
    }

    #[test]
    fn a4_five_cell_step_has_one_consistent_partition_and_the_all_stuck_start_solves() {
        let l = a4_law();
        let laws: Vec<&FrictionLaw> = (0..5).map(|_| &l).collect();
        // The all-stuck partition is the walk's first try. The matrix's smallest
        // eigenvalue is 2.4e-13 m/N, and a difference-quotient Jacobian (h = 1e-11 N)
        // is wrong by about 0.7 of that: Newton stalled near 5e-12 m against the
        // 1e-12 m gate, so the solve reported failure instead of the true answer,
        // that port 0 would need ~10.8 N against a 0.043 N cone.
        let stuck = [Branch::Stuck; 5];
        let (a, f, n) = (&A4_5_COMPLIANCE, &A4_5_FREE_SLIP, &A4_5_NORMAL);
        let t = solve_partition(a, f, n, &laws, A4_DT, &stuck, 100)
            .expect("the all-stuck partition is linear and must solve");
        let cap = l.static_coefficient().unwrap() * n[0];
        let reach = t[0].abs();
        assert!(reach > 100.0 * cap, "all-stuck needs far more than the cone holds: {reach} vs {cap}");
        let consistent = is_consistent(a, f, n, &laws, A4_DT, &stuck, &t).unwrap();
        assert!(!consistent, "the all-stuck partition is not consistent here");
        let found = consistent_partitions(&A4_5_COMPLIANCE, &A4_5_FREE_SLIP, &A4_5_NORMAL);
        assert_eq!(
            found.len(),
            1,
            "expected a unique consistent partition, got {:?}",
            found.iter().map(|(c, _)| c).collect::<Vec<_>>()
        );
        let (cells, t) = &found[0];
        assert_eq!(
            *cells,
            [Branch::Stuck, Branch::Stuck, Branch::SlideNeg, Branch::Stuck, Branch::Stuck]
        );
        assert_matches_reference(t, &A4_5_REFERENCE);
    }

    #[test]
    fn the_walk_lands_each_a4_partition_in_its_pivot_count_and_refuses_one_pivot_short() {
        let cases: [(&[f64], &[f64], &[f64], &[f64], usize, &[ModalFrictionRegime]); 2] = [
            (
                &A4_3_COMPLIANCE, &A4_3_FREE_SLIP, &A4_3_NORMAL, &A4_3_REFERENCE, 3,
                &[ModalFrictionRegime::Sliding, ModalFrictionRegime::Sticking, ModalFrictionRegime::Sliding],
            ),
            (
                &A4_5_COMPLIANCE, &A4_5_FREE_SLIP, &A4_5_NORMAL, &A4_5_REFERENCE, 4,
                &[
                    ModalFrictionRegime::Sticking,
                    ModalFrictionRegime::Sticking,
                    ModalFrictionRegime::Sliding,
                    ModalFrictionRegime::Sticking,
                    ModalFrictionRegime::Sticking,
                ],
            ),
        ];
        for (a, f, normal, reference, pivots, regimes) in cases {
            match walk(a, f, normal, pivots - 1) {
                Resolution::Exhausted { pivots: spent, last_set, worst } => {
                    assert!(spent < pivots, "a budget of {} must stop short of the {pivots} pivots needed", pivots - 1);
                    assert!(!last_set.is_empty() && worst.condition != "none");
                }
                Resolution::Solved(..) => panic!("{} pivots cannot reach a partition that needs {pivots}", pivots - 1),
            }
            for budget in [pivots, 4 * f.len() + 16] {
                match walk(a, f, normal, budget) {
                    Resolution::Solved(t, got) => {
                        assert_eq!(got, regimes, "budget {budget}");
                        assert_matches_reference(&t, reference);
                    }
                    Resolution::Exhausted { pivots: spent, last_set, worst } => panic!(
                        "budget {budget} must resolve the step: {spent} pivots, last set {last_set}, worst {worst:?}"
                    ),
                }
            }
        }
    }

    #[test]
    fn a_state_with_no_consistent_partition_is_refused_not_invented() {
        // A port with zero normal load: stuck needs a nonzero traction inside a zero
        // cone, and sliding needs a strictly positive traction from a zero load. No
        // partition is consistent, so the rescue must exhaust, not return one.
        let l = a4_law();
        let laws = [&l];
        let seed = vec![ModalFrictionRegime::Sticking];
        match solve_active_set(&[3.0e-6], &[1.0e-7], &[0.0], &laws, A4_DT, &seed, 8).unwrap() {
            Resolution::Exhausted { last_set, worst, .. } => {
                assert!(!last_set.is_empty() && worst.condition != "none");
            }
            Resolution::Solved(t, _) => panic!("a zero-load port has no consistent partition, got tractions {t:?}"),
        }
    }

    // Captured refusal state of the violin E5 sitebow card with the measured rosin
    // row `friction rosin-wg2025-eq13` (3-cell ribbon, 144 kHz). It was dumped by a
    // scratch capture build (MagentaBadger 16900, 2026-11-03) at the step after
    // accepted sample 10788, where the sweep failed and the walk refused: "5/28
    // pivots, last active set S-+, worst KKT violation at port 0
    // (slide-pos-wrong-direction-or-nondissipative) size 1e0".
    //
    // The previous step was SS+. Port 1 now needs 2.6e-5 N more than its cone
    // holds, so SS+ has ceased to exist. Solving every partition from a grid of
    // starts finds exactly one consistent state, S-+: port 1 breaks away at
    // -3.1 mm/s. The walk visits S-+, but its partition Newton cycles across
    // port 1's zero-slip kink.
    const E5_LINE: [f64; 9] = [
        5.1734155698535495e-6, 2.9809085288452823e-6, 3.0751603963720907e-6,
        2.9809085288452823e-6, 5.099163599374646e-6, 2.718315640216884e-6,
        3.075160396372091e-6, 2.7183156402168843e-6, 5.262846745918796e-6,
    ];
    const E5_FREE_SLIP: [f64; 3] = [2.747519612686609e-7, -1.5907963135027875e-8, 3.627451607228994e-7];
    const E5_NORMAL: [f64; 3] = [0.0529; 3];
    const E5_DT: f64 = 1.0 / 144_000.0;
    // The S-+ root found by the multi-start |slip| Newton, independent of the speed solve.
    const E5_S_MINUS_PLUS: [f64; 3] = [0.05895761563373555, -0.05818005488303967, 0.04655664028650187];

    /// Every partition the speed solve resolves to a consistent state, with its tractions.
    fn consistent_by_speed(a: &[f64], f: &[f64], n: &[f64], laws: &[&FrictionLaw], dt: f64) -> Vec<(Vec<Branch>, Vec<f64>)> {
        let p = f.len();
        let bs = all_branches();
        let mut found = Vec::new();
        for code in 0..3usize.pow(p as u32) {
            let cells: Vec<Branch> = (0..p).map(|i| bs[code / 3usize.pow(i as u32) % 3]).collect();
            if let Some(t) = solve_partition_by_speed(a, f, n, laws, dt, &cells, 100) {
                if is_consistent(a, f, n, laws, dt, &cells, &t).unwrap() {
                    found.push((cells, t));
                }
            }
        }
        found
    }

    #[test]
    fn the_eq13_onset_step_refuses_on_the_walk_and_its_one_consistent_state_is_rescued() {
        use Branch::{SlideNeg, SlidePos, Stuck};
        let law = FrictionLaw::rosin_woodhouse_galluzzo_2025_eq13();
        let laws = [&law, &law, &law];
        let (a, f, n) = (&E5_LINE, &E5_FREE_SLIP, &E5_NORMAL);
        let seed = [ModalFrictionRegime::Sticking, ModalFrictionRegime::Sticking, ModalFrictionRegime::Sliding];
        let previous = [Stuck, Stuck, SlidePos];
        let s_minus_plus = [Stuck, SlideNeg, SlidePos];

        // Before: the walk refuses with the render's own payload, because its
        // partition solve cannot resolve the one consistent state.
        match solve_active_set(a, f, n, &laws, E5_DT, &seed, 28).unwrap() {
            Resolution::Exhausted { pivots, last_set, worst } => assert_eq!(
                (pivots, last_set.as_str(), worst.port, worst.condition, worst.size),
                (5, "S-+", 0, "slide-pos-wrong-direction-or-nondissipative", 1.0)
            ),
            Resolution::Solved(t, _) => panic!("the walk is the captured refusal, got tractions {t:?}"),
        }
        assert!(
            solve_partition(a, f, n, &laws, E5_DT, &s_minus_plus, 100).is_none(),
            "the walk's partition solve cycles on S-+"
        );
        let cap = law.static_coefficient().unwrap() * n[1];
        let ss_plus = solve_partition(a, f, n, &laws, E5_DT, &[Stuck, Stuck, SlidePos], 100).unwrap();
        assert!(ss_plus[1].abs() > cap && ss_plus[1].abs() < cap * (1.0 + 1e-3), "SS+ is just past its cone");

        let found = consistent_by_speed(a, f, n, &laws, E5_DT);
        assert_eq!(
            found.iter().map(|(c, _)| c.clone()).collect::<Vec<_>>(),
            vec![s_minus_plus.to_vec()],
            "the speed solve finds exactly the one consistent state"
        );

        // After: the rescue lands it, and it passes the unchanged gate and KKT test.
        match solve_active_set_with_hysteresis(a, f, n, &laws, E5_DT, &seed, &previous, 28).unwrap() {
            Resolution::Solved(t, regimes) => {
                assert_eq!(
                    regimes,
                    [ModalFrictionRegime::Sticking, ModalFrictionRegime::Sliding, ModalFrictionRegime::Sliding]
                );
                assert!(is_consistent(a, f, n, &laws, E5_DT, &s_minus_plus, &t).unwrap());
                let slip = slips(a, f, &t);
                assert!(slip[0].abs() <= 1e-12, "port 0 stays stuck: slip {}", slip[0]);
                for i in 1..3 {
                    let kt = kinetic_traction(&law, n[i], branch_sign(s_minus_plus[i]), slip[i], E5_DT).unwrap();
                    assert!((t[i] - kt).abs() <= 1e-12, "port {i}: traction {} against kinetic {kt}", t[i]);
                }
                for (i, (&got, &want)) in t.iter().zip(&E5_S_MINUS_PLUS).enumerate() {
                    assert!((got - want).abs() <= 1e-12, "port {i}: {got} against the multi-start root {want}");
                }
                let breakaway = slip[1] / E5_DT;
                assert!((breakaway + 3.102e-3).abs() < 1e-5, "port 1 breaks away at {breakaway} m/s");
            }
            Resolution::Exhausted { pivots, last_set, worst } => panic!(
                "the one consistent state must be rescued: {pivots} pivots, last set {last_set}, worst {worst:?}"
            ),
        }
    }

    #[test]
    fn a_step_with_no_consistent_partition_keeps_the_walks_refusal_through_the_last_resort() {
        // The zero-load port above: no partition is consistent, so the last resort
        // finds nothing and the walk's exhaustion comes back unchanged.
        let l = a4_law();
        let laws = [&l];
        let seed = vec![ModalFrictionRegime::Sticking];
        let (a, f, n) = ([3.0e-6], [1.0e-7], [0.0]);
        assert!(nearest_consistent_partition(&a, &f, &n, &laws, A4_DT, &[Branch::Stuck]).unwrap().is_none());
        let walked = solve_active_set(&a, &f, &n, &laws, A4_DT, &seed, 8).unwrap();
        let rescued = solve_active_set_with_hysteresis(&a, &f, &n, &laws, A4_DT, &seed, &[Branch::Stuck], 8).unwrap();
        match (walked, rescued) {
            (
                Resolution::Exhausted { pivots: p0, last_set: s0, worst: w0 },
                Resolution::Exhausted { pivots: p1, last_set: s1, worst: w1 },
            ) => assert_eq!(
                (p0, s0, w0.port, w0.condition, w0.size.to_bits()),
                (p1, s1, w1.port, w1.condition, w1.size.to_bits())
            ),
            _ => panic!("a zero-load port has no consistent partition: both must refuse"),
        }
    }

    #[test]
    fn a_step_the_walk_resolves_returns_the_walks_tractions_bit_for_bit() {
        // The last resort runs only after an exhaustion. `previous` is set unlike
        // every walk's answer, so any leak of the enumeration would show.
        let check = |a: &[f64], f: &[f64], n: &[f64], laws: &[&FrictionLaw], dt: f64, budget: usize| {
            let p = f.len();
            let seed = vec![ModalFrictionRegime::Sticking; p];
            let previous = vec![Branch::SlideNeg; p];
            let walked = solve_active_set(a, f, n, laws, dt, &seed, budget).unwrap();
            let rescued = solve_active_set_with_hysteresis(a, f, n, laws, dt, &seed, &previous, budget).unwrap();
            match (walked, rescued) {
                (Resolution::Solved(t0, r0), Resolution::Solved(t1, r1)) => {
                    assert_eq!(r0, r1);
                    let bits = |t: &[f64]| t.iter().map(|x| x.to_bits()).collect::<Vec<_>>();
                    assert_eq!(bits(&t0), bits(&t1));
                }
                _ => panic!("the walk resolves this step, and so must the rescue"),
            }
        };
        let (a, f) = operating_point();
        let stribeck = law();
        check(&a, &f, &NORMAL, &[&stribeck; 5], DT, 64);
        let coulomb = FrictionLaw::Coulomb { static_mu: 0.8, kinetic_mu: 0.4 };
        check(&STOP_LINE, &STOP_FREE_SLIP, &STOP_NORMAL, &[&coulomb; 3], STOP_DT, 4 * 3 + 16);
        let a4 = a4_law();
        check(&A4_3_COMPLIANCE, &A4_3_FREE_SLIP, &A4_3_NORMAL, &[&a4; 3], A4_DT, 4 * 3 + 16);
        check(&A4_5_COMPLIANCE, &A4_5_FREE_SLIP, &A4_5_NORMAL, &[&a4; 5], A4_DT, 4 * 5 + 16);
    }

    #[test]
    fn the_speed_solve_finds_the_walks_consistent_partitions_on_the_captured_steps() {
        // The start behind `critical_speed` rests on an argument, not a proof.
        // These are the steps it was checked on: on each, the speed solve finds
        // the same consistent partitions as the walk's own solve, with the same
        // tractions.
        let same = |a: &[f64], f: &[f64], n: &[f64], laws: &[&FrictionLaw], dt: f64| -> usize {
            let p = f.len();
            let bs = all_branches();
            let mut walk_solve = Vec::new();
            for code in 0..3usize.pow(p as u32) {
                let cells: Vec<Branch> = (0..p).map(|i| bs[code / 3usize.pow(i as u32) % 3]).collect();
                if let Some(t) = solve_partition(a, f, n, laws, dt, &cells, 200) {
                    if is_consistent(a, f, n, laws, dt, &cells, &t).unwrap() {
                        walk_solve.push((cells, t));
                    }
                }
            }
            let by_speed = consistent_by_speed(a, f, n, laws, dt);
            assert_eq!(
                by_speed.iter().map(|(c, _)| c.clone()).collect::<Vec<_>>(),
                walk_solve.iter().map(|(c, _)| c.clone()).collect::<Vec<_>>()
            );
            for ((cells, t), (_, w)) in by_speed.iter().zip(&walk_solve) {
                for (i, (&got, &want)) in t.iter().zip(w).enumerate() {
                    assert!((got - want).abs() <= 1e-9, "{cells:?} port {i}: {got} against {want}");
                }
            }
            by_speed.len()
        };
        let (a, f) = operating_point();
        let stribeck = law();
        assert_eq!(same(&a, &f, &NORMAL, &[&stribeck; 5], DT), 1);
        let coulomb = FrictionLaw::Coulomb { static_mu: 0.8, kinetic_mu: 0.4 };
        assert_eq!(same(&STOP_LINE, &STOP_FREE_SLIP, &STOP_NORMAL, &[&coulomb; 3], STOP_DT), 3);
        let a4 = a4_law();
        assert_eq!(same(&A4_3_COMPLIANCE, &A4_3_FREE_SLIP, &A4_3_NORMAL, &[&a4; 3], A4_DT), 1);
        assert_eq!(same(&A4_5_COMPLIANCE, &A4_5_FREE_SLIP, &A4_5_NORMAL, &[&a4; 5], A4_DT), 1);

        // On the stop bank's three states, the previous step decides. SS+ continues
        // an SS+ step. From rest, SS+ and S-S each differ in one port: the rule
        // cannot break that tie, so it refuses rather than pick one.
        use Branch::{SlidePos, Stuck};
        let laws = [&coulomb; 3];
        let (a, f, n) = (&STOP_LINE, &STOP_FREE_SLIP, &STOP_NORMAL);
        let (_, regimes) = nearest_consistent_partition(a, f, n, &laws, STOP_DT, &[Stuck, Stuck, SlidePos])
            .unwrap()
            .expect("SS+ is nearest an SS+ step");
        assert_eq!(regimes, [ModalFrictionRegime::Sticking, ModalFrictionRegime::Sticking, ModalFrictionRegime::Sliding]);
        assert!(nearest_consistent_partition(a, f, n, &laws, STOP_DT, &[Stuck; 3]).unwrap().is_none());
    }
}
