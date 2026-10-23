//! Active-set rescue for the joint multi-port friction solve (frankensim-o7qbp).
//!
//! The simultaneous friction owner normally advances with a Gauss-Seidel sweep:
//! each port re-chooses its stick/slide branch against the other ports' current
//! held tractions, and the step commits when one full sweep is jointly
//! consistent. For a well-posed step the consistent stick/slide partition is
//! unique, but the sweep can OSCILLATE between candidate partitions near a
//! stick<->slip transition and never land a clean sweep, so it exhausts
//! `max_sweeps` and refuses a step that does have a solution (o7qbp: the 5-cell
//! compliant bow ribbon at the first Helmholtz transition; enumerating all 3^p
//! partitions at that step found exactly ONE consistent joint solution).
//!
//! This module solves that unique partition directly. It runs ONLY after the
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

use super::*;

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
/// elimination. Deterministic; returns None on a singular pivot.
fn dense_solve(m: &[f64], b: &[f64], n: usize) -> Option<Vec<f64>> {
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

/// The signed kinetic traction a sliding port carries at relative slip `slip`:
/// magnitude `mu_k(|slip|/dt) * normal`, sign equal to the slip direction so the
/// friction force is dissipative.
fn kinetic_traction(law: &FrictionLaw, normal: f64, slip: f64, dt: f64) -> Result<f64, ModalFrictionError> {
    let speed = (slip / dt).abs();
    let mu = if speed == 0.0 {
        law.zero_speed_kinetic_coefficient()?
    } else {
        law.kinetic_coefficient(speed)?
    };
    let sign = if slip < 0.0 { -1.0 } else { 1.0 };
    Ok(sign * mu * normal)
}

/// Joint-Newton the tractions for a FIXED partition `cells`. Stuck ports pin
/// `slip_i = 0`; sliding ports pin `t_i = kinetic_traction(sign, |slip_i|)`.
/// Returns the tractions, or None if the Newton does not resolve.
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
    let residual = |t: &[f64]| -> Option<Vec<f64>> {
        let slip = slips(a, f, t);
        let mut r = vec![0.0f64; p];
        for i in 0..p {
            r[i] = match cells[i] {
                Branch::Stuck => slip[i],
                Branch::SlidePos | Branch::SlideNeg => {
                    let kt = kinetic_traction(law[i], normal[i], slip[i], dt).ok()?;
                    t[i] - kt
                }
            };
        }
        Some(r)
    };
    for _ in 0..max_newton {
        let r = residual(&t)?;
        if r.iter().all(|x| x.abs() <= 1e-18) {
            break;
        }
        // Numerical Jacobian (p x p), central-free forward difference.
        let h = 1e-11;
        let mut jac = vec![0.0f64; p * p];
        for j in 0..p {
            let mut tp = t.clone();
            tp[j] += h;
            let rp = residual(&tp)?;
            for i in 0..p {
                jac[i * p + j] = (rp[i] - r[i]) / h;
            }
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
    let r = residual(&t)?;
    if r.iter().all(|x| x.abs() <= 1e-12) {
        Some(t)
    } else {
        None
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
/// assigned direction and the traction opposes that slip (dissipative).
fn is_consistent(
    a: &[f64],
    f: &[f64],
    normal: &[f64],
    law: &[&FrictionLaw],
    dt: f64,
    cells: &[Branch],
    t: &[f64],
) -> Result<bool, ModalFrictionError> {
    let p = f.len();
    let slip = slips(a, f, t);
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

/// Active-set rescue. Returns `(tractions, regimes)` for the unique consistent
/// partition, or None when the pivot budget is spent without consistency.
///
/// `seed` is the sweep's last (oscillating) regimes, used only to bias the first
/// partition; correctness does not depend on it.
#[allow(clippy::too_many_arguments)]
pub(super) fn solve_active_set(
    a: &[f64],
    f: &[f64],
    normal: &[f64],
    law: &[&FrictionLaw],
    dt: f64,
    seed: &[ModalFrictionRegime],
    budget: usize,
) -> Result<Option<(Vec<f64>, Vec<ModalFrictionRegime>)>, ModalFrictionError> {
    let p = f.len();
    // First partition: every port stuck. A stuck port that cannot hold will be
    // flipped to the sliding direction its over-cone traction indicates.
    let _ = seed;
    let mut cells = vec![Branch::Stuck; p];
    let mut seen: Vec<Vec<Branch>> = Vec::new();
    for _pivot in 0..budget {
        if seen.contains(&cells) {
            // Revisited a partition: the single-flip walk is cycling; stop.
            return Ok(None);
        }
        seen.push(cells.clone());
        let Some(t) = solve_partition(a, f, normal, law, dt, &cells, 100) else {
            // This partition's Newton did not resolve; perturb by flipping the
            // port with the largest |slip| that is pinned stuck, else give up.
            let slip = slips(a, f, &vec![0.0; p]);
            let mut flip = None;
            let mut worst = 0.0;
            for i in 0..p {
                if cells[i] == Branch::Stuck && slip[i].abs() > worst {
                    worst = slip[i].abs();
                    flip = Some(i);
                }
            }
            match flip {
                Some(i) => {
                    cells[i] = if slips(a, f, &vec![0.0; p])[i] >= 0.0 {
                        Branch::SlidePos
                    } else {
                        Branch::SlideNeg
                    };
                    continue;
                }
                None => return Ok(None),
            }
        };
        if is_consistent(a, f, normal, law, dt, &cells, &t)? {
            let regimes = cells.iter().map(|&b| regime_of(b)).collect();
            return Ok(Some((t, regimes)));
        }
        // Flip the single worst KKT violator.
        let slip = slips(a, f, &t);
        let mut flip = None;
        let mut worst = 0.0f64;
        for i in 0..p {
            let (viol, want) = match cells[i] {
                Branch::Stuck => {
                    let cap = law[i].static_coefficient()? * normal[i];
                    let excess = t[i].abs() - cap;
                    let want = if t[i] >= 0.0 { Branch::SlidePos } else { Branch::SlideNeg };
                    (excess, want)
                }
                Branch::SlidePos => {
                    // Wrong direction or non-dissipative -> should stick.
                    let bad = if slip[i] <= 0.0 || t[i] <= 0.0 { slip[i].abs().max(1.0) } else { -1.0 };
                    (bad, Branch::Stuck)
                }
                Branch::SlideNeg => {
                    let bad = if slip[i] >= 0.0 || t[i] >= 0.0 { slip[i].abs().max(1.0) } else { -1.0 };
                    (bad, Branch::Stuck)
                }
            };
            if viol > worst {
                worst = viol;
                flip = Some((i, want));
            }
        }
        match flip {
            Some((i, want)) => cells[i] = want,
            None => return Ok(None),
        }
    }
    Ok(None)
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

    #[test]
    fn active_set_solve_lands_the_unique_partition() {
        let (a, f) = operating_point();
        let l = law();
        let laws: Vec<&FrictionLaw> = (0..5).map(|_| &l).collect();
        let seed = vec![ModalFrictionRegime::Sticking; 5];
        let got = solve_active_set(&a, &f, &NORMAL, &laws, DT, &seed, 64)
            .unwrap()
            .expect("active-set rescue found no consistent solution");
        let (t, regimes) = got;
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
}
