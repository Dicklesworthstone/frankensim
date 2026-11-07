//! The tridiagonal solve against the dense owner: every nonzero result must be bit-identical to
//! `lu_complex(dense with structural zeros).solve`, the singular step must be the same, and the cases must really
//! reach the swap, tie and singular paths (a test that only sees diagonally dominant matrices proves nothing).
use super::super::{EigFailure, lu_complex};
use super::*;

struct Rng(u64);
impl Rng {
    fn next(&mut self) -> u64 {
        // xorshift64*: deterministic, no dependency.
        self.0 ^= self.0 >> 12;
        self.0 ^= self.0 << 25;
        self.0 ^= self.0 >> 27;
        self.0.wrapping_mul(0x2545_F491_4F6C_DD1D)
    }
    fn unit(&mut self) -> f64 {
        (self.next() >> 11) as f64 / (1u64 << 53) as f64 * 2.0 - 1.0
    }
    fn complex(&mut self) -> C64 {
        C64::new(self.unit(), self.unit())
    }
    /// A small integer in -2..=2: ties, exact zeros and singular steps.
    fn small(&mut self) -> f64 {
        (self.next() % 5) as f64 - 2.0
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
enum Flavor {
    /// Diagonally dominant: no row swap.
    Dominant,
    /// Sub-diagonal far above the diagonal: a swap at nearly every step, so multipliers wander down.
    SubDominant,
    /// Row scales drawn from 1e-2, 1 and 1e2.
    Mixed,
    /// Some diagonal entries 1e-9 or exactly zero.
    TinyDiagonal,
    /// Entries from -2..=2 on both parts: exact magnitude ties, exact zeros, singular matrices.
    Integers,
}
const FLAVORS: [Flavor; 5] =
    [Flavor::Dominant, Flavor::SubDominant, Flavor::Mixed, Flavor::TinyDiagonal, Flavor::Integers];

struct Tri {
    lower: Vec<C64>,
    diag: Vec<C64>,
    upper: Vec<C64>,
}

fn tridiagonal(n: usize, flavor: Flavor, rng: &mut Rng) -> Tri {
    let mut entry = |scale: f64| -> C64 {
        if flavor == Flavor::Integers { C64::new(rng.small(), rng.small()) } else { rng.complex().scale(scale) }
    };
    let mut diag = Vec::with_capacity(n);
    let mut lower = Vec::new();
    let mut upper = Vec::new();
    for i in 0..n {
        let scale = match flavor {
            Flavor::Mixed => [1.0e-2, 1.0, 1.0e2][i % 3],
            _ => 1.0,
        };
        let d = match flavor {
            Flavor::Dominant => entry(1.0) + C64::new(4.0, 0.5),
            Flavor::SubDominant => entry(0.1),
            Flavor::TinyDiagonal => match i % 4 {
                1 => entry(1.0e-9),
                3 => C64::ZERO,
                _ => entry(1.0),
            },
            _ => entry(scale),
        };
        diag.push(d);
        if i + 1 < n {
            lower.push(match flavor {
                Flavor::SubDominant => entry(10.0),
                _ => entry(scale),
            });
            upper.push(entry(scale));
        }
    }
    Tri { lower, diag, upper }
}

fn dense(t: &Tri) -> Vec<C64> {
    let n = t.diag.len();
    let mut m = vec![C64::ZERO; n * n];
    for i in 0..n {
        m[i * n + i] = t.diag[i];
        if i + 1 < n {
            m[i * n + i + 1] = t.upper[i];
            m[(i + 1) * n + i] = t.lower[i];
        }
    }
    m
}

/// Bits with the sign of a zero removed: the one thing the sparse graph is allowed to differ in.
fn bits(v: &[C64]) -> Vec<(u64, u64)> {
    v.iter().map(|c| ((c.re + 0.0).to_bits(), (c.im + 0.0).to_bits())).collect()
}

#[derive(Default)]
struct Coverage {
    cases: usize,
    solved: usize,
    singular: usize,
    swapped: usize,
    wandering: usize,
}

/// One matrix, `columns` right-hand sides: dense against tridiagonal.
fn compare(t: &Tri, columns: usize, rng: &mut Rng, cover: &mut Coverage) {
    let n = t.diag.len();
    let rhs: Vec<C64> = (0..n * columns).map(|_| rng.complex()).collect();
    cover.cases += 1;
    let mut sparse = rhs.clone();
    let got = solve_tridiagonal_complex(&t.lower, &t.diag, &t.upper, &mut sparse, columns);
    match lu_complex(&dense(t), n) {
        Err(EigFailure { window_hi }) => {
            assert_eq!(got, Err(LuComplexError::Singular { index: window_hi }), "n = {n}");
            assert_eq!(bits(&sparse), bits(&rhs), "a singular matrix must leave the right-hand sides alone");
            cover.singular += 1;
        }
        Ok(lu) => {
            assert_eq!(got, Ok(()), "n = {n}: the dense factorization succeeded");
            let mut expected = rhs.clone();
            for column in expected.chunks_exact_mut(n) {
                lu.solve(column);
            }
            assert_eq!(bits(&sparse), bits(&expected), "n = {n}, {columns} right-hand sides");
            cover.solved += 1;
            if lu.perm.iter().enumerate().any(|(i, &p)| i != p) {
                cover.swapped += 1;
            }
            // A multiplier below the first sub-diagonal exists only after consecutive row swaps.
            if (0..n).any(|i| (0..i.saturating_sub(1)).any(|k| lu.data[i * n + k] != C64::ZERO)) {
                cover.wandering += 1;
            }
        }
    }
}

#[test]
fn tridiagonal_solve_is_bit_identical_to_the_dense_owner() {
    let mut cover = Coverage::default();
    let sizes: Vec<usize> = (1..=40).chain([57, 63, 64, 65, 66, 100, 113, 128, 140]).collect();
    for &n in &sizes {
        let seeds = if n > 64 { 2 } else { 6 };
        for (f, &flavor) in FLAVORS.iter().enumerate() {
            for seed in 0..seeds {
                let mut rng = Rng(0x9E37_79B9_7F4A_7C15 ^ ((n as u64) << 20) ^ ((f as u64) << 8) ^ seed as u64);
                let t = tridiagonal(n, flavor, &mut rng);
                for columns in [1, 3, 5] {
                    compare(&t, columns, &mut rng, &mut cover);
                }
            }
        }
    }
    // The cases must exercise what the argument is about.
    assert!(cover.solved > 1000, "{} solved", cover.solved);
    assert!(cover.singular >= 10, "{} singular", cover.singular);
    assert!(cover.swapped > 500, "{} with a row swap", cover.swapped);
    assert!(cover.wandering >= 100, "{} with a multiplier below the sub-diagonal", cover.wandering);
}

#[test]
fn nonsingular_solutions_have_a_small_residual() {
    // Guards against the dense owner and the sparse routine agreeing on a wrong answer.
    let mut rng = Rng(0xC0FF_EE00_1234_5678);
    for n in [1, 2, 3, 7, 31, 100] {
        let t = tridiagonal(n, Flavor::Mixed, &mut rng);
        let b: Vec<C64> = (0..n).map(|_| rng.complex()).collect();
        let mut x = b.clone();
        if solve_tridiagonal_complex(&t.lower, &t.diag, &t.upper, &mut x, 1).is_err() {
            continue;
        }
        let m = dense(&t);
        let mut worst = 0.0_f64;
        let mut scale = 0.0_f64;
        for i in 0..n {
            let mut s = C64::ZERO;
            for j in 0..n {
                s = s + m[i * n + j] * x[j];
            }
            worst = worst.max((s - b[i]).abs());
            scale = scale.max(b[i].abs() + m[i * n..(i + 1) * n].iter().zip(&x).map(|(a, y)| (*a * *y).abs()).sum::<f64>());
        }
        assert!(worst <= 1.0e-9 * scale, "n = {n}: residual {worst:e} against {scale:e}");
    }
}

#[test]
fn empty_and_single_systems() {
    let mut none: Vec<C64> = Vec::new();
    assert_eq!(solve_tridiagonal_complex(&[], &[], &[], &mut none, 3), Ok(()));
    let d = C64::new(3.0, -4.0);
    let mut x = vec![C64::new(1.0, 2.0), C64::new(-5.0, 0.5)];
    let expected: Vec<C64> = x.iter().map(|&v| v / d).collect();
    assert_eq!(solve_tridiagonal_complex(&[], &[d], &[], &mut x, 2), Ok(()));
    assert_eq!(bits(&x), bits(&expected));
    let mut y = vec![C64::new(1.0, 2.0)];
    assert_eq!(
        solve_tridiagonal_complex(&[], &[C64::ZERO], &[], &mut y, 1),
        Err(LuComplexError::Singular { index: 0 })
    );
}

#[test]
fn a_zero_pivot_after_elimination_names_its_step() {
    // [[1, 1, 0], [1, 1, 0], [0, 0, 1]]: step 0 leaves a zero at (1, 1) and the candidate (2, 1) is zero too, so the
    // matrix is singular at step 1, as the dense owner reports.
    let one = C64::ONE;
    let mut x = vec![one, one, one];
    assert_eq!(solve_tridiagonal_complex(&[one, C64::ZERO], &[one, one, one], &[one, C64::ZERO], &mut x, 1), Err(LuComplexError::Singular { index: 1 }));
    assert_eq!(bits(&x), bits(&[one, one, one]));
}

#[test]
#[should_panic(expected = "rhs must hold n * columns entries")]
fn a_wrong_right_hand_side_length_is_refused() {
    let mut x = vec![C64::ONE; 5];
    let _ = solve_tridiagonal_complex(&[C64::ONE; 2], &[C64::ONE; 3], &[C64::ONE; 2], &mut x, 2);
}
