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
    let mu = if speed == 0.0 {
        law.zero_speed_kinetic_coefficient()?
    } else {
        law.kinetic_coefficient(speed)?
    };
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
/// stalled above the 1e-12 m gate. Sliding rows keep the difference: their law
/// derivative would otherwise be needed for every `FrictionLaw` variant.
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
                    let sign = branch_sign(cells[i]);
                    let kt = kinetic_traction(law[i], normal[i], sign, slip[i], dt).ok()?;
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
        for i in 0..p {
            if cells[i] == Branch::Stuck {
                for j in 0..p {
                    jac[i * p + j] = -a[i * p + j];
                }
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
) -> Result<Resolution, ModalFrictionError> {
    let p = f.len();
    // First partition: every port stuck. A stuck port that cannot hold will be
    // flipped to the sliding direction its over-cone traction indicates.
    let _ = seed;
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
        let Some(t) = solve_partition(a, f, normal, law, dt, &cells, 100) else {
            // This partition's Newton did not resolve; perturb by flipping the
            // port with the largest |slip| that is pinned stuck, else give up.
            let slip = slips(a, f, &vec![0.0; p]);
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
                    worst = KktViolation { port: i, condition: "newton-nonconvergence", size: w };
                    cells[i] = if slip[i] >= 0.0 { Branch::SlidePos } else { Branch::SlideNeg };
                    continue;
                }
                None => {
                    worst = KktViolation { port: 0, condition: "newton-nonconvergence", size: f64::INFINITY };
                    return Ok(Resolution::Exhausted { pivots, last_set: encode(&cells), worst });
                }
            }
        };
        if is_consistent(a, f, normal, law, dt, &cells, &t)? {
            let regimes = cells.iter().map(|&b| regime_of(b)).collect();
            return Ok(Resolution::Solved(t, regimes));
        }
        // Flip the single worst KKT violator, recording port/condition/size.
        let slip = slips(a, f, &t);
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
}
