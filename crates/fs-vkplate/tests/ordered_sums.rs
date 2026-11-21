//! `ordered_sums::sums_of_products` against the plain `.sum()` it stands in for, bit for bit: whole blocks
//! of eight rows, the rows after the last whole block, and blocks whose rows differ in length. The file
//! names nothing of the crate (the module is included by path), so it also compiles alone with `rustc --test`.

#[path = "../src/ordered_sums.rs"]
mod ordered_sums;

/// Values over seven decades with both signs and some signed zeros, so that a sum taken in another
/// order would round differently.
fn mixed(len: usize, seed: u64) -> Vec<f64> {
    let mut s = seed.wrapping_mul(0x9e37_79b9_7f4a_7c15) | 1;
    (0..len)
        .map(|_| {
            s ^= s << 13;
            s ^= s >> 7;
            s ^= s << 17;
            match s % 23 {
                0 => 0.0,
                1 => -0.0,
                _ => {
                    let unit = (s >> 11) as f64 / (1u64 << 53) as f64 - 0.5;
                    unit * 10f64.powi((s % 7) as i32 - 3)
                }
            }
        })
        .collect()
}

fn plain(row: &[f64], x: &[f64]) -> f64 {
    row.iter().zip(x).map(|(a, b)| a * b).sum()
}

fn check(rows: &[Vec<f64>], x: &[f64]) {
    let mut out = vec![f64::NAN; rows.len()];
    ordered_sums::sums_of_products(|r| rows[r].as_slice(), x, &mut out);
    for (r, row) in rows.iter().enumerate() {
        assert_eq!(out[r].to_bits(), plain(row, x).to_bits(), "row {r} of {}, length {}", rows.len(), row.len());
    }
}

#[test]
fn equal_rows_have_the_bits_of_their_plain_sums() {
    for count in [0, 1, 7, 8, 9, 15, 16, 17, 347] {
        for len in [0, 1, 5, 8, 13, 64, 4197] {
            let rows: Vec<_> = (0..count).map(|r| mixed(len, 1000 * count as u64 + r as u64)).collect();
            check(&rows, &mixed(len, 7 + len as u64));
        }
    }
}

#[test]
fn rows_longer_or_shorter_than_x_stop_where_zip_stops() {
    let x = mixed(100, 3);
    // a whole block cut at x's end, a block mixing rows shorter and longer than x, and a remainder that mixes them too
    let mut rows: Vec<_> = (0..8).map(|r| mixed(130, 50 + r)).collect();
    rows.extend((0..8).map(|r| mixed([40, 100, 130, 0][r as usize % 4], 60 + r)));
    rows.extend((0..5).map(|r| mixed([99, 101, 3][r as usize % 3], 70 + r)));
    check(&rows, &x);
}

#[test]
fn a_matrix_read_by_rows_has_the_bits_of_its_row_sums() {
    // the modal-space product's form: n rows of one row-major n x n matrix against a vector of length n
    for n in [1, 8, 11, 347] {
        let matrix = mixed(n * n, 11 + n as u64);
        let v = mixed(n, 13);
        let mut out = vec![f64::NAN; n];
        ordered_sums::sums_of_products(|i| &matrix[i * n..(i + 1) * n], &v, &mut out);
        for i in 0..n {
            let want: f64 = (0..n).map(|j| matrix[i * n + j] * v[j]).sum();
            assert_eq!(out[i].to_bits(), want.to_bits(), "n {n}, row {i}");
        }
    }
}
