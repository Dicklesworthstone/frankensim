//! The prepared operator's two orders by rows (the loss forms over the upper triangle, the storage over the whole
//! square), each held to a plain loop of the same order. The forms and the products read copies by blocks of eight
//! rows (`UpperBlocks`, `RowBlocks`); the plain loops read the row-major square. The forms are held to it on one
//! thread and with their rows shared out to several. The held-force defect is held to the fold that always computes
//! its bounds.
//!
//! fs-couple's lib test crate does not build under the memory guard, so this target compiles the kernel's own
//! source file a second time. That file names nothing of its crate.

#[path = "../src/render/schedule/force/coupled/exact/prepared/upper_rows.rs"]
mod upper_rows;

use std::cell::Cell;

use upper_rows::{
    RowBlocks, StateForms, UpperBlocks, bitwise_symmetric, held_force_defect, identity_rows,
    product_from, row_from, row_storage, state_sums, upper_row_form, upper_row_form_from,
    upper_row_forms_from,
};

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
        let (sum, magnitude) = upper_row_form(&UpperBlocks::new(&w, side), &z);
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
            // The upper-row form's copy does not read the lower triangle, which is why the check comes first.
            let z = mixed(side, seed + 100);
            let form = |w: &[f64]| upper_row_form(&UpperBlocks::new(w, side), &z);
            assert_eq!(form(&moved), form(&w));
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

/// The defect of two forces past their bounds with the bounds always computed, the fold the prepared owner made
/// before `held_force_defect`.
fn bounded_defect(a: &[f64], b: &[f64], bounds: &[f64]) -> f64 {
    a.iter()
        .zip(b)
        .zip(bounds)
        .fold(0.0_f64, |m, ((x, y), r)| m.max((x - y).abs() - r))
}

#[test]
fn a_held_force_defect_computes_its_bounds_only_past_the_tolerance_and_refuses_as_with_them() {
    // Up to nine forces from 1e-300 to 1e12, the second equal, a few units in the last place off, or 1e-6 off; bounds
    // of zero or up to twice each difference.
    let scales = [1e-300, 1e-17, 1.0, 1e12].into_iter().cycle();
    for ((seed, n), scale) in (0..600_u64).zip((1..=9).cycle()).zip(scales) {
        let a: Vec<f64> = mixed(n, seed).iter().map(|x| x * scale).collect();
        let b: Vec<f64> = a
            .iter()
            .zip(mixed(n, seed + 1000))
            .map(|(&x, u)| match seed % 3 {
                0 => x,
                1 => x + 4.0 * u * x.abs() * f64::EPSILON,
                _ => x * (1.0 + 1e-6 * u),
            })
            .collect();
        let bounds: Vec<f64> = a
            .iter()
            .zip(&b)
            .zip(mixed(n, seed + 2000))
            .map(|((x, y), u)| {
                if seed % 5 == 0 {
                    0.0
                } else {
                    (x - y).abs() * 2.0 * u.abs()
                }
            })
            .collect();
        let plain = a
            .iter()
            .zip(&b)
            .fold(0.0_f64, |m, (x, y)| m.max((x - y).abs()));
        let below = if plain > 0.0 {
            f64::from_bits(plain.to_bits() - 1)
        } else {
            0.0
        };
        for tolerance in [0.0, below, plain, 0.5 * plain, 2.0 * plain] {
            let ran = Cell::new(false);
            let defect = held_force_defect::<()>(&a, &b, tolerance, || {
                ran.set(true);
                Ok(bounds.clone())
            })
            .unwrap();
            let full = bounded_defect(&a, &b, &bounds);
            let at = format!("seed {seed} tolerance {tolerance:e} plain {plain:e}");
            assert_eq!(defect > tolerance, full > tolerance, "{at}");
            assert_eq!(ran.get(), plain > tolerance, "{at}");
            let expected = if plain > tolerance { full } else { plain };
            assert_eq!(defect.to_bits(), expected.to_bits(), "{at}");
        }
    }
    // An error of the bounds is met only past the tolerance: within it the plain defect stands.
    let (a, b) = ([1.0, 2.0], [1.0, 2.0 + 1e-12]);
    let plain = (2.0 + 1e-12) - 2.0;
    assert_eq!(
        held_force_defect(&a, &b, 1e-9, || Err::<Vec<f64>, _>("overflow")),
        Ok(plain)
    );
    assert_eq!(
        held_force_defect(&a, &b, 1e-13, || Err::<Vec<f64>, _>("overflow")),
        Err("overflow")
    );
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

/// One left-to-right sum for each whole row, from 0.0: the bits a product must have.
fn plain_product(matrix: &[f64], z: &[f64]) -> Vec<u64> {
    let row_sum = |row: &[f64]| row.iter().zip(z).fold(0.0, |s, (a, b)| s + a * b);
    matrix
        .chunks_exact(z.len())
        .map(|row| row_sum(row).to_bits())
        .collect()
}

#[test]
fn a_product_continued_from_a_states_kept_sums_has_the_bits_of_a_plain_product() {
    let bits = |product: Vec<f64>| product.iter().map(|v| v.to_bits()).collect::<Vec<_>>();
    // Each product takes the next of 1, 2, 3 and 7 threads, so two products in a row run on two counts.
    let mut counts = [1, 2, 3, 7].into_iter().cycle();
    let mut next_product = |slot: &mut _, blocks: &_, z: &[f64], state| {
        product_from(slot, blocks, z, state, counts.next().unwrap_or(1))
    };
    for (seed, side) in (1..).zip(SIDES) {
        // All the rows of a square, then its first two thirds (67 rows at side 100: not whole blocks of eight).
        for rows in [side, (2 * side).div_ceil(3)] {
            let matrix = mixed(rows * side, seed);
            let blocks = RowBlocks::new(&matrix, side);
            // No state, a third of the columns, all but the last, all of them.
            for state in [0, side / 3, side - 1, side] {
                let x = mixed(state, seed + 10);
                let z = |x: &[f64], k: u64| [x, &mixed(side - state, seed + 20 + k)].concat();
                let mut slot = None;
                // The first product computes the state's sums; the next two continue from them under other forces.
                for k in 0..3 {
                    let z = z(&x, k);
                    let product = bits(next_product(&mut slot, &blocks, &z, state));
                    let plain = plain_product(&matrix, &z);
                    assert_eq!(product, plain, "{side} {rows} {state} {k}");
                }
                // A state that differs in the last bit of one entry, the first and then the last, is another state.
                for at in [0, state.saturating_sub(1)].into_iter().take(state.min(2)) {
                    let mut y = x.clone();
                    y[at] = f64::from_bits(y[at].to_bits() ^ 1);
                    let z = z(&y, 0);
                    let product = bits(next_product(&mut slot, &blocks, &z, state));
                    let plain = plain_product(&matrix, &z);
                    assert_eq!(product, plain, "{side} {rows} {state} at {at}");
                }
            }
        }
    }
}

#[test]
fn a_row_continued_from_a_states_kept_sums_has_the_bits_of_a_plain_product() {
    // A state's sums are kept on the next of 1, 2, 3 and 7 threads, by `state_sums` or by the product; a call on a
    // state that is held starts none.
    let mut counts = [1, 2, 3, 7].into_iter().cycle();
    for (seed, side) in (1..).zip(SIDES) {
        // All the rows of a square, then its first two thirds (67 rows at side 100: not whole blocks of eight).
        for rows in [side, (2 * side).div_ceil(3)] {
            let matrix = mixed(rows * side, seed);
            let blocks = RowBlocks::new(&matrix, side);
            let check = |sums: &[f64], z: &[f64], state, plain: &[u64], case: &str| {
                for (index, row) in matrix.chunks_exact(side).enumerate() {
                    let continued = row_from(sums, index, row, z, state).to_bits();
                    assert_eq!(continued, plain[index], "{case} row {index}");
                }
            };
            // No state, a third of the columns, all but the last, all of them.
            for state in [0, side / 3, side - 1, side] {
                let x = mixed(state, seed + 10);
                let force = |k: u64| mixed(side - state, seed + 20 + k);
                let mut slot = None;
                // The first z keeps the state's sums; the next two continue from them under other forces, and so does
                // a product on the same state.
                for k in 0..3 {
                    let z = [&x[..], &force(k)].concat();
                    let case = format!("{side} {rows} {state} {k}");
                    let plain = plain_product(&matrix, &z);
                    let threads = counts.next().unwrap_or(1);
                    let sums = state_sums(&mut slot, &blocks, &x, threads).to_vec();
                    check(&sums, &z, state, &plain, &case);
                    let product = product_from(&mut slot, &blocks, &z, state, 1);
                    let bits: Vec<_> = product.iter().map(|v| v.to_bits()).collect();
                    assert_eq!(bits, plain, "{case} product");
                }
                // A state that differs in the last bit of one entry, the first and then the last, is another state:
                // the product keeps its sums, and `state_sums` reads them.
                for at in [0, state.saturating_sub(1)].into_iter().take(state.min(2)) {
                    let mut y = x.clone();
                    y[at] = f64::from_bits(y[at].to_bits() ^ 1);
                    let z = [&y[..], &force(0)].concat();
                    let case = format!("{side} {rows} {state} at {at}");
                    let plain = plain_product(&matrix, &z);
                    let threads = counts.next().unwrap_or(1);
                    let product = product_from(&mut slot, &blocks, &z, state, threads);
                    let bits: Vec<_> = product.iter().map(|v| v.to_bits()).collect();
                    assert_eq!(bits, plain, "{case} product");
                    let sums = state_sums(&mut slot, &blocks, &y, 1).to_vec();
                    check(&sums, &z, state, &plain, &case);
                }
            }
        }
    }
}

#[test]
fn a_form_continued_from_a_states_kept_sums_has_the_bits_of_a_plain_loop() {
    let bits = |(sum, magnitude): (f64, f64)| (sum.to_bits(), magnitude.to_bits());
    for (seed, side) in (1..).zip(SIDES) {
        let w = gramian(side, seed);
        let upper = UpperBlocks::new(&w, side);
        // No state; lengths below, at and above a block of eight; two thirds of the side; all but the last; all.
        for state in [0, 7, 8, 9, 2 * side / 3, side - 1, side] {
            let state = state.min(side);
            let x = mixed(state, seed + 10);
            let z = |x: &[f64], k: u64| [x, &mixed(side - state, seed + 20 + k)].concat();
            let mut slot = None;
            // The first form computes the state's sums; the next two continue them under other forces.
            for k in 0..3 {
                let z = z(&x, k);
                let form = bits(upper_row_form_from(&mut slot, &upper, &z, state));
                assert_eq!(form, bits(plain(&w, &z)), "{side} {state} {k}");
            }
            // A state that differs in one entry, the first and then the last, is another state: by its last bit, and
            // by a half, which moves every row's kept sum (one bit of one term can vanish in a sum of many).
            for at in [0, state.saturating_sub(1)].into_iter().take(state.min(2)) {
                for other in [f64::from_bits(x[at].to_bits() ^ 1), x[at] + 0.5] {
                    let mut y = x.clone();
                    y[at] = other;
                    let z = z(&y, 0);
                    let form = bits(upper_row_form_from(&mut slot, &upper, &z, state));
                    assert_eq!(form, bits(plain(&w, &z)), "{side} {state} at {at}");
                }
            }
        }
    }
}

#[test]
fn forms_shared_out_to_threads_have_the_bits_of_a_plain_loop() {
    let bits = |(sum, magnitude): (f64, f64)| (sum.to_bits(), magnitude.to_bits());
    // Sides of one block, of a whole and a short block (a split inside the last block), and of many blocks.
    for (seed, side) in (1..).zip([1, 7, 9, 61, 100]) {
        let w: Vec<_> = (0..3).map(|k| gramian(side, 10 * seed + k)).collect();
        let upper: Vec<_> = w.iter().map(|w| UpperBlocks::new(w, side)).collect();
        for state in [0, side / 2, side] {
            let x = mixed(state, seed + 30);
            for threads in [1, 2, 3, 7] {
                let mut slots: Vec<Option<StateForms>> = (0..3).map(|_| None).collect();
                // The first sweep computes the state's sums; the second continues them under another force.
                for k in 0..2 {
                    let z = [&x[..], &mixed(side - state, seed + 40 + k)].concat();
                    let mut forms: Vec<_> = upper.iter().zip(&mut slots).collect();
                    let forms = upper_row_forms_from(&mut forms, &z, state, threads);
                    let plain: Vec<_> = w.iter().map(|w| bits(plain(w, &z))).collect();
                    let forms: Vec<_> = forms.into_iter().map(bits).collect();
                    assert_eq!(forms, plain, "{side} {state} {threads} threads, sweep {k}");
                }
            }
        }
    }
}
