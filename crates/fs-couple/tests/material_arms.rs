//! The Maxwell arms' chains (`relaxation/arms.rs`), each held to the plain loop of its order bit for bit, and the stage
//! held to the loop it replaced, which added and checked every term: the same storage, the same forces and the same
//! refusals.
//!
//! fs-couple's lib test crate does not build under the memory guard, so this target compiles the kernel's own source
//! file a second time. That file names nothing of its crate.

#[path = "../src/render/schedule/force/coupled/relaxation/arms.rs"]
mod arms;

use arms::{arm_chains, staged_arms};

/// `len` values spread over (-1, 1) with mixed exponents and signs, the same for the same seed.
fn mixed(len: usize, seed: u64) -> Vec<f64> {
    let mut state = seed.wrapping_mul(0x9E37_79B9_7F4A_7C15) | 1;
    (0..len)
        .map(|_| {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            (state >> 11) as f64 / (1_u64 << 53) as f64 * 2.0 - 1.0
        })
        .collect()
}

#[test]
fn each_arms_chain_has_the_bits_of_its_own_fold_for_any_stride_and_length() {
    let bits = |v: &[f64]| v.iter().map(|v| v.to_bits()).collect::<Vec<_>>();
    for (seed, count) in (1..).zip([0, 1, 7, 8, 9, 17]) {
        for modes in [0, 1, 5, 64] {
            let q = mixed(modes, seed + 50);
            // Rows shorter than, as long as and longer than `q` reaches at each stride, so the arms of a block stop apart.
            let lengths = [modes / 2, modes, 2 * modes, 2 * modes + 3];
            let rows: Vec<_> = (0..count)
                .map(|a| mixed(lengths[a % 4], seed + a as u64))
                .collect();
            let (starts, scales) = (mixed(count, seed + 60), mixed(count, seed + 70));
            let arms: Vec<_> = rows
                .iter()
                .enumerate()
                .map(|(a, row)| (starts[a], scales[a], &row[..]))
                .collect();
            for (stride, chains) in [
                (1, arm_chains::<1>(&arms, &q)),
                (2, arm_chains::<2>(&arms, &q)),
            ] {
                let fold = |(start, scale, row): &(f64, f64, &[f64])| {
                    let terms = row.iter().step_by(stride).zip(&q);
                    terms.fold(*start, |sum, (b, x)| sum + scale * b * x)
                };
                let folds: Vec<_> = arms.iter().map(fold).collect();
                assert_eq!(bits(&chains), bits(&folds), "{count} {modes} {stride}");
            }
        }
    }
}

fn finite(value: f64) -> Result<f64, ()> {
    value.is_finite().then_some(value).ok_or(())
}

/// The stage as relaxation.rs wrote it at 8fa3491d19: each arm's strain over its whole projection, every term
/// checked, then its mismatch, the storage and the forces.
fn plain_stage(
    arms: &[(f64, &[f64])],
    memory: &[f64],
    x: &[f64],
    forces: &mut [f64],
) -> Result<f64, ()> {
    let mut energy = 0.0;
    for (&(root, projection), &z) in arms.iter().zip(memory) {
        let mut terms = projection.iter().zip(x);
        let strain = terms.try_fold(0.0, |sum, (b, x)| finite(sum + root * b * x))?;
        let mismatch = finite(strain - z)?;
        energy = finite(energy + 0.5 * mismatch * mismatch)?;
        for (force, &b) in forces.iter_mut().zip(projection.iter().step_by(2)) {
            *force = finite(*force - root * b * mismatch)?;
        }
    }
    Ok(energy)
}

#[test]
fn the_stage_has_the_bits_and_the_refusals_of_the_loop_that_checks_and_adds_every_term() {
    let (count, modes) = (19, 9);
    let roots: Vec<_> = mixed(count, 7).iter().map(|r| 1.25 + 0.75 * r).collect();
    // Mixed displacement coefficients; velocity coefficients +0 and -0, as the admission keeps them.
    let mut rows: Vec<_> = (0..count)
        .map(|a| mixed(2 * modes, 17 + a as u64))
        .collect();
    for row in &mut rows {
        for (i, b) in row.iter_mut().enumerate().skip(1).step_by(2) {
            *b = [0.0, -0.0][i / 2 % 2];
        }
    }
    let (memory, forces) = (mixed(count, 12), mixed(modes, 9));
    // Velocities of every finite size; modes 1 and 2 far out, mode 3 at rest. Then every displacement +-0 as well,
    // so that each strain is a sum of signed zeros, which a chain from +0.0 leaves +0.
    let mut x = mixed(2 * modes, 8);
    for (i, v) in x.iter_mut().enumerate().skip(1).step_by(2) {
        *v = [*v, -0.0, 0.0, 1e300, -f64::MAX][i / 2 % 5];
    }
    (x[2], x[4], x[6]) = (1e10, 1e10, 0.0);
    let mut zeros = x.clone();
    for (i, q) in zeros.iter_mut().enumerate().step_by(2) {
        *q = [0.0, -0.0][i / 2 % 2];
    }
    // The first case changes nothing. Each other changes one arm, in the first block, the second or the last: its
    // strain overflows to +inf; to +inf and then to -inf (NaN); its kept strain is +inf; its mismatch, about 1e150
    // with the storage finite, overflows its fourth force after its first three have moved.
    let cases: [fn(&mut [f64], &mut f64); 5] = [
        |_, _| {},
        |p, _| p[2] = 1e300,
        |p, _| (p[2], p[4]) = (1e300, -1e300),
        |_, z| *z = f64::INFINITY,
        |p, z| (p[6], *z) = (1e200, -1e150),
    ];
    let bits = |v: &[f64]| v.iter().map(|v| v.to_bits()).collect::<Vec<_>>();
    for (case, change) in cases.iter().enumerate() {
        let states = [(3, &x), (9, &x), (18, &x), (18, &zeros)];
        for (arm, x) in states.into_iter().take(if case == 0 { 4 } else { 3 }) {
            let (mut rows, mut memory) = (rows.clone(), memory.clone());
            change(&mut rows[arm], &mut memory[arm]);
            let arms: Vec<_> = roots.iter().zip(&rows).map(|(r, p)| (*r, &p[..])).collect();
            let (mut new, mut old) = (forces.clone(), forces.clone());
            let staged = staged_arms(&arms, &memory, x, &mut new, finite);
            let plain = plain_stage(&arms, &memory, x, &mut old);
            let outcome = |stage: Result<f64, ()>| stage.map(f64::to_bits);
            let at = format!("case {case}, arm {arm}");
            assert_eq!(
                (outcome(staged), bits(&new)),
                (outcome(plain), bits(&old)),
                "{at}"
            );
            // The loop accepts the first case and refuses the others, after the arms before it staged their forces.
            assert_eq!(plain.is_ok(), case == 0, "{at}");
            assert!(case == 0 || new != forces, "{at}");
        }
    }
}
