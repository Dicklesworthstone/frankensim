//! Const-evaluated Cayley tables (plan §7.7): the multiplication tables
//! for Cl(3,0,1) and Cl(4,1) are computed AT COMPILE TIME from the metric
//! signatures — no runtime blade bookkeeping, no handwritten tables to
//! typo. Blades are basis-vector bitmasks in canonical ascending order;
//! every product below is a table lookup plus a fused sign.
//!
//! Audit: `cargo test -p fs-ga table::` prints/checks the generated
//! tables against first-principles identities (metric squares, anti-
//! commutation, associativity at the blade level).

/// One Cayley-table entry: `e_a * e_b = sign * e_blade` (sign 0 when the
/// degenerate metric annihilates the product).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Term {
    /// −1, 0, or +1.
    pub sign: i8,
    /// Result blade bitmask.
    pub blade: u8,
}

/// Sign of reordering the concatenation of blades `a`,`b` (each in
/// canonical ascending order) into canonical order: (−1)^(#transpositions).
#[must_use]
pub const fn reorder_sign(a: u32, b: u32) -> i8 {
    let mut rest = a >> 1;
    let mut swaps = 0u32;
    while rest != 0 {
        swaps += (rest & b).count_ones();
        rest >>= 1;
    }
    if swaps & 1 == 0 { 1 } else { -1 }
}

/// Blade grade (number of basis vectors).
#[must_use]
pub const fn grade(blade: u32) -> u32 {
    blade.count_ones()
}

/// Reverse sign for a blade of grade g: (−1)^(g(g−1)/2).
#[must_use]
pub const fn reverse_sign(blade: u32) -> i8 {
    let g = grade(blade);
    if (g * (g.wrapping_sub(1)) / 2) & 1 == 0 {
        1
    } else {
        -1
    }
}

/// Grade-involution sign: (−1)^g.
#[must_use]
pub const fn involute_sign(blade: u32) -> i8 {
    if grade(blade) & 1 == 0 { 1 } else { -1 }
}

/// Build the full geometric-product Cayley table for a diagonal metric.
/// `D` basis vectors, `B = 2^D` blades — both given explicitly because
/// stable const generics cannot derive one from the other.
#[must_use]
#[allow(clippy::cast_possible_truncation)] // B <= 256 keeps blades in u8
pub const fn build_table<const D: usize, const B: usize>(metric: [i8; D]) -> [[Term; B]; B] {
    let mut table = [[Term { sign: 0, blade: 0 }; B]; B];
    let mut i = 0;
    while i < B {
        let mut j = 0;
        while j < B {
            let a = i as u32;
            let b = j as u32;
            let mut sign = reorder_sign(a, b) as i32;
            let common = a & b;
            let mut k = 0;
            while k < D {
                if (common >> k) & 1 == 1 {
                    sign *= metric[k] as i32;
                }
                k += 1;
            }
            table[i][j] = Term {
                sign: sign as i8,
                blade: (a ^ b) as u8,
            };
            j += 1;
        }
        i += 1;
    }
    table
}

/// PGA Cl(3,0,1): basis (e0, e1, e2, e3) with e0² = 0 (the projective /
/// degenerate direction), e1² = e2² = e3² = +1. Bit k ↔ e_k.
pub const PGA_METRIC: [i8; 4] = [0, 1, 1, 1];
/// Number of PGA blades.
pub const PGA_BLADES: usize = 16;
/// The PGA Cayley table as a `const` (readable by downstream const-fn
/// codegen, e.g. the monomorphized motor kernels).
pub const PGA_TABLE_CONST: [[Term; PGA_BLADES]; PGA_BLADES] =
    build_table::<4, PGA_BLADES>(PGA_METRIC);
/// The compile-time PGA Cayley table (single runtime copy).
pub static PGA_TABLE: [[Term; PGA_BLADES]; PGA_BLADES] = PGA_TABLE_CONST;

/// CGA Cl(4,1): basis (e1, e2, e3, e+, e−) with e1²=e2²=e3²=e+²=+1,
/// e−²=−1. Bit 0..=2 ↔ e1..=e3, bit 3 ↔ e+, bit 4 ↔ e−.
pub const CGA_METRIC: [i8; 5] = [1, 1, 1, 1, -1];
/// Number of CGA blades.
pub const CGA_BLADES: usize = 32;
/// The compile-time CGA Cayley table.
pub static CGA_TABLE: [[Term; CGA_BLADES]; CGA_BLADES] = build_table::<5, CGA_BLADES>(CGA_METRIC);

/// Right complement: blade `rc` with `a ∧ rc = +pseudoscalar` (used for
/// the Poincaré duality that makes PGA's regressive product / join work
/// despite the degenerate metric).
#[must_use]
#[allow(clippy::cast_possible_truncation)] // full < 256
pub const fn right_complement(a: u32, full: u32) -> Term {
    let c = full ^ a;
    Term {
        sign: reorder_sign(a, c),
        blade: c as u8,
    }
}

/// Left complement: blade `lc` with `lc ∧ a = +pseudoscalar`.
#[must_use]
#[allow(clippy::cast_possible_truncation)] // full < 256
pub const fn left_complement(a: u32, full: u32) -> Term {
    let c = full ^ a;
    Term {
        sign: reorder_sign(c, a),
        blade: c as u8,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn metric_squares_come_out_of_the_generator() {
        // e0² = 0 in PGA; e1²..e3² = +1.
        assert_eq!(PGA_TABLE[1][1], Term { sign: 0, blade: 0 });
        for k in 1..4 {
            let b = 1usize << k;
            assert_eq!(PGA_TABLE[b][b], Term { sign: 1, blade: 0 });
        }
        // CGA: e−² = −1, others +1.
        let em = 1usize << 4;
        assert_eq!(CGA_TABLE[em][em], Term { sign: -1, blade: 0 });
        for k in 0..4 {
            let b = 1usize << k;
            assert_eq!(CGA_TABLE[b][b], Term { sign: 1, blade: 0 });
        }
    }

    #[test]
    fn distinct_basis_vectors_anticommute() {
        for i in 0..4u32 {
            for j in 0..4u32 {
                if i == j {
                    continue;
                }
                let (a, b) = (1usize << i, 1usize << j);
                let ab = PGA_TABLE[a][b];
                let ba = PGA_TABLE[b][a];
                assert_eq!(ab.blade, ba.blade);
                assert_eq!(ab.sign, -ba.sign, "e{i} e{j} must anticommute");
            }
        }
    }

    #[test]
    fn blade_level_associativity_both_algebras() {
        // (e_a e_b) e_c == e_a (e_b e_c) including the fused metric signs —
        // exhaustively over all blade triples of both algebras.
        fn check<const B: usize>(table: &[[Term; B]; B]) {
            for a in 0..B {
                for b in 0..B {
                    let ab = table[a][b];
                    for c in 0..B {
                        let left = {
                            let t = table[ab.blade as usize][c];
                            (i32::from(ab.sign) * i32::from(t.sign), t.blade)
                        };
                        let bc = table[b][c];
                        let right = {
                            let t = table[a][bc.blade as usize];
                            (i32::from(bc.sign) * i32::from(t.sign), t.blade)
                        };
                        if left.0 == 0 && right.0 == 0 {
                            continue; // both annihilated
                        }
                        assert_eq!(left, right, "associativity broke at ({a},{b},{c})");
                    }
                }
            }
        }
        check(&PGA_TABLE);
        check(&CGA_TABLE);
    }

    #[test]
    fn complements_wedge_to_the_pseudoscalar() {
        let full = (PGA_BLADES - 1) as u32;
        for a in 0..PGA_BLADES as u32 {
            let rc = right_complement(a, full);
            // a ∧ rc has no common bits, so its sign is the pure reorder
            // sign — which the complement definition sets to +1.
            assert_eq!(reorder_sign(a, u32::from(rc.blade)) * rc.sign, 1);
            let lc = left_complement(a, full);
            assert_eq!(reorder_sign(u32::from(lc.blade), a) * lc.sign, 1);
        }
    }

    #[test]
    fn reverse_and_involution_signs() {
        // Grade pattern for reverse: + + − − + (mod 4 in grade).
        assert_eq!(reverse_sign(0b0000), 1);
        assert_eq!(reverse_sign(0b0001), 1);
        assert_eq!(reverse_sign(0b0011), -1);
        assert_eq!(reverse_sign(0b0111), -1);
        assert_eq!(reverse_sign(0b1111), 1);
        assert_eq!(involute_sign(0b0001), -1);
        assert_eq!(involute_sign(0b0011), 1);
    }
}
