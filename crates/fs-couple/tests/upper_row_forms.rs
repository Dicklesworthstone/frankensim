//! The prepared operator's two orders by rows (the loss forms over the upper triangle, the storage over the whole
//! square), each held to a plain loop of the same order.
//!
//! fs-couple's lib test crate does not build under the memory guard, so this target compiles the kernel's own
//! source file a second time. That file names nothing of its crate.

#[path = "../src/render/schedule/force/coupled/exact/prepared/upper_rows.rs"]
mod upper_rows;

use upper_rows::{bitwise_symmetric, identity_rows, row_storage, upper_row_form};

/// Sides around the blocks of eight rows, and one near a hundred.
const SIDES: [usize; 7] = [1, 7, 8, 9, 16, 17, 100];

/// Numbers in (-1, 1) of both signs.
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

/// `A^T A` of a mixed-sign `A` with `side + 3` rows, each entry computed once and written to both halves: positive
/// semidefinite, of both signs off the diagonal, and equal to its transpose bit for bit.
fn gramian(side: usize, seed: u64) -> Vec<f64> {
    let a = mixed((side + 3) * side, seed);
    let mut w = vec![0.0; side * side];
    for i in 0..side {
        for j in i..side {
            let value: f64 = a.chunks_exact(side).map(|row| row[i] * row[j]).sum();
            w[i * side + j] = value;
            w[j * side + i] = value;
        }
    }
    w
}

/// The same order as the kernel's, one row at a time: row i sums its columns j > i from left to right, and the rows
/// are added in row order.
fn plain(w: &[f64], z: &[f64]) -> (f64, f64) {
    let (mut sum, mut magnitude) = (0.0, 0.0);
    for (i, (row, left)) in w.chunks_exact(z.len()).zip(z).enumerate() {
        let (mut partial, mut size) = (0.0, 0.0);
        for (value, right) in row[i + 1..].iter().zip(&z[i + 1..]) {
            let term = value * right;
            partial += term;
            size += term.abs();
        }
        let diagonal = row[i] * left;
        sum += left * (diagonal + 2.0 * partial);
        magnitude += left.abs() * (diagonal.abs() + 2.0 * size);
    }
    (sum, magnitude)
}

/// Today's order: every term of the square, row-major, in one chain.
fn chain(w: &[f64], z: &[f64]) -> (f64, f64) {
    let (mut sum, mut magnitude) = (0.0, 0.0);
    for (row, left) in w.chunks_exact(z.len()).zip(z) {
        for (value, right) in row.iter().zip(z) {
            let term = left * value * right;
            sum += term;
            magnitude += term.abs();
        }
    }
    (sum, magnitude)
}

#[test]
fn the_upper_row_form_has_the_bits_of_a_plain_loop_of_its_order_and_agrees_with_the_chain() {
    for (seed, side) in (1..).zip(SIDES) {
        let (w, z) = (gramian(side, seed), mixed(side, seed + 100));
        assert!(
            w.iter().any(|v| *v < 0.0) || side == 1,
            "side {side}: no negative entry"
        );
        let (sum, magnitude) = upper_row_form(&w, &z);
        let (plain_sum, plain_magnitude) = plain(&w, &z);
        assert_eq!(sum.to_bits(), plain_sum.to_bits(), "side {side}: sum");
        assert_eq!(
            magnitude.to_bits(),
            plain_magnitude.to_bits(),
            "side {side}: magnitude"
        );
        let (chain_sum, chain_magnitude) = chain(&w, &z);
        assert!(sum > 0.0 && chain_sum > 0.0);
        assert!(
            (sum - chain_sum).abs() <= 1e-13 * chain_sum,
            "side {side}: {sum:e} vs {chain_sum:e}"
        );
        let off = (magnitude - chain_magnitude).abs();
        assert!(
            off <= 1e-13 * chain_magnitude,
            "side {side}: {magnitude:e} vs {chain_magnitude:e}"
        );
    }
}

#[test]
fn a_matrix_is_symmetric_only_if_every_pair_has_the_same_bits() {
    for (seed, side) in (1..).zip(SIDES) {
        let w = gramian(side, seed);
        assert!(bitwise_symmetric(&w, side), "side {side}");
        assert!(
            !bitwise_symmetric(&w, side + 1),
            "side {side}: a wrong side"
        );
        if side > 1 {
            // One entry of the lower triangle moved by one unit in its last place.
            let mut moved = w.clone();
            let at = (side - 1) * side;
            moved[at] = f64::from_bits(moved[at].to_bits() + 1);
            assert!(!bitwise_symmetric(&moved, side), "side {side}: one ulp");
            // The upper-row form does not read the lower triangle, which is why the check comes first.
            let z = mixed(side, seed + 100);
            assert_eq!(upper_row_form(&moved, &z), upper_row_form(&w, &z));
        }
    }
    // Zero and its negative are two bit patterns.
    assert!(!bitwise_symmetric(&[1.0, 0.0, -0.0, 1.0], 2));
}

#[test]
fn rows_are_the_identitys_only_if_every_entry_has_its_bits() {
    for (seed, side) in (1..).zip(SIDES) {
        for from in [0, side / 2, side - 1] {
            // Mixed rows, then the identity's rows from `from` on.
            let mut square = mixed(side * side, seed);
            for (at, entry) in square.iter_mut().enumerate().skip(from * side) {
                *entry = f64::from(u8::from(at / side == at % side));
            }
            assert!(identity_rows(&square, side, from), "{side} {from}");
            assert!(!identity_rows(&square, side + 1, from), "{side}");
            // In the last row: the one moved by a unit in its last place, then, off the diagonal, a zero made the
            // smallest positive number and the negative zero.
            let last = (side - 1) * side;
            let one = (last + side - 1, f64::from_bits(1.0_f64.to_bits() + 1));
            let wrong = [one, (last, f64::from_bits(1)), (last, -0.0)];
            for (at, value) in wrong.into_iter().take(if side > 1 { 3 } else { 1 }) {
                let mut moved = square.clone();
                moved[at] = value;
                assert!(!identity_rows(&moved, side, from), "{side} {from} {at}");
            }
        }
    }
    // A row before `from` is not read, and from `side` on there is no row to hold.
    assert!(identity_rows(&[f64::NAN, 2.0, 0.0, 1.0], 2, 1));
    assert!(identity_rows(&[f64::NAN; 4], 2, 2));
}

/// A stiffness as `mechanical_stiffness` leaves it: a positive diagonal, then three links added entry by entry as
/// `k * c[i] * c[j]`. That product is `(k c_i) c_j` at (i, j) and `(k c_j) c_i` at (j, i).
fn stiffness(side: usize, seed: u64) -> Vec<f64> {
    let mut k = vec![0.0; side * side];
    for (i, d) in mixed(side, seed).iter().enumerate() {
        k[i * side + i] = 2.0 + d;
    }
    for (link, column) in mixed(3 * side, seed + 300).chunks_exact(side).enumerate() {
        for i in 0..side {
            for j in 0..side {
                k[i * side + j] += (0.3 + link as f64) * column[i] * column[j];
            }
        }
    }
    k
}

/// The same order as the kernel's, one row at a time: row i sums every column from left to right, and the rows are
/// added in row order.
fn plain_storage(k: &[f64], modes: &[f64]) -> f64 {
    let side = modes.len() / 2;
    let mut energy = 0.0;
    for i in 0..side {
        let mut sum = 0.0;
        for j in 0..side {
            sum += k[i * side + j] * modes[2 * j];
        }
        energy += 0.5 * modes[2 * i + 1] * modes[2 * i + 1] + 0.5 * modes[2 * i] * sum;
    }
    energy
}

/// Today's order: each mode's kinetic term and then its row of the square, in one chain.
fn chain_storage(k: &[f64], modes: &[f64]) -> f64 {
    let side = modes.len() / 2;
    let mut energy = 0.0;
    for i in 0..side {
        energy += 0.5 * modes[2 * i + 1].powi(2);
        for j in 0..side {
            energy += 0.5 * modes[2 * i] * k[i * side + j] * modes[2 * j];
        }
    }
    energy
}

#[test]
fn the_row_storage_has_the_bits_of_a_plain_loop_of_its_order_and_agrees_with_the_chain() {
    for (seed, side) in (1..).zip(SIDES) {
        let (k, modes) = (stiffness(side, seed), mixed(2 * side, seed + 200));
        // The kernel reads the whole square, so it is held on a square that differs from its transpose.
        assert!(
            !bitwise_symmetric(&k, side) || side == 1,
            "side {side}: the stiffness equals its transpose"
        );
        let rows = row_storage(&k, &modes);
        assert_eq!(
            rows.to_bits(),
            plain_storage(&k, &modes).to_bits(),
            "side {side}"
        );
        let chain = chain_storage(&k, &modes);
        assert!(
            chain > 0.0 && (rows - chain).abs() <= 1e-13 * chain,
            "side {side}: {rows:e} vs {chain:e}"
        );
    }
    assert_eq!(row_storage(&[], &[]).to_bits(), 0.0_f64.to_bits());
}
