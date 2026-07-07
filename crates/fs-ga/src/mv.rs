//! Dense multivector types over the const-evaluated Cayley tables. One
//! macro instantiation per algebra — fully monomorphized, no runtime
//! blade bookkeeping (the tables and complements are compile-time data).
//!
//! Determinism (P2): every product accumulates in a fixed (row-major
//! blade-index) order, so results are bit-identical run to run and
//! platform to platform (pure f64 adds/muls).

use crate::table::{Term, grade, involute_sign, left_complement, reverse_sign, right_complement};

macro_rules! algebra {
    ($(#[$meta:meta])* $name:ident, $blades:expr, $table:path) => {
        $(#[$meta])*
        #[derive(Debug, Clone, Copy, PartialEq)]
        pub struct $name(pub [f64; $blades]);

        impl Default for $name {
            fn default() -> Self {
                Self::zero()
            }
        }

        impl $name {
            /// Number of basis blades.
            pub const BLADES: usize = $blades;
            /// Pseudoscalar blade index (all basis vectors).
            pub const PSEUDO: usize = $blades - 1;

            /// The zero multivector.
            #[must_use]
            pub fn zero() -> Self {
                Self([0.0; $blades])
            }

            /// A scalar.
            #[must_use]
            pub fn scalar(s: f64) -> Self {
                let mut v = [0.0; $blades];
                v[0] = s;
                Self(v)
            }

            /// One basis blade with a coefficient.
            ///
            /// # Panics
            /// If `blade >= Self::BLADES`.
            #[must_use]
            pub fn blade(blade: usize, coeff: f64) -> Self {
                assert!(blade < $blades, "blade index out of range");
                let mut v = [0.0; $blades];
                v[blade] = coeff;
                Self(v)
            }

            /// Scalar (grade-0) part.
            #[must_use]
            pub fn scalar_part(&self) -> f64 {
                self.0[0]
            }

            /// Geometric product, fixed accumulation order.
            #[must_use]
            pub fn gp(&self, rhs: &Self) -> Self {
                let mut out = [0.0f64; $blades];
                for (i, &a) in self.0.iter().enumerate() {
                    if a == 0.0 {
                        continue;
                    }
                    for (j, &b) in rhs.0.iter().enumerate() {
                        let t: Term = $table[i][j];
                        if t.sign != 0 && b != 0.0 {
                            out[t.blade as usize] += f64::from(t.sign) * a * b;
                        }
                    }
                }
                Self(out)
            }

            /// Outer (wedge) product: the grade-raising, metric-free part.
            #[must_use]
            pub fn wedge(&self, rhs: &Self) -> Self {
                let mut out = [0.0f64; $blades];
                for (i, &a) in self.0.iter().enumerate() {
                    if a == 0.0 {
                        continue;
                    }
                    for (j, &b) in rhs.0.iter().enumerate() {
                        if i & j != 0 || b == 0.0 {
                            continue; // shared basis vector ⇒ wedge dies
                        }
                        let t: Term = $table[i][j];
                        out[t.blade as usize] += f64::from(t.sign) * a * b;
                    }
                }
                Self(out)
            }

            /// Left contraction `self ⌋ rhs` (grade-lowering inner product).
            #[must_use]
            pub fn lcontract(&self, rhs: &Self) -> Self {
                let mut out = [0.0f64; $blades];
                for (i, &a) in self.0.iter().enumerate() {
                    if a == 0.0 {
                        continue;
                    }
                    let gi = grade(i as u32);
                    for (j, &b) in rhs.0.iter().enumerate() {
                        let gj = grade(j as u32);
                        if gj < gi || b == 0.0 {
                            continue;
                        }
                        let t: Term = $table[i][j];
                        if t.sign != 0 && grade(u32::from(t.blade)) == gj - gi {
                            out[t.blade as usize] += f64::from(t.sign) * a * b;
                        }
                    }
                }
                Self(out)
            }

            /// Regressive product (join/meet dual): `J⁻¹(J(a) ∧ J(b))`
            /// with `J` the right complement — Poincaré duality, valid in
            /// degenerate metrics.
            #[must_use]
            pub fn vee(&self, rhs: &Self) -> Self {
                let dual_lhs = self.right_complement();
                let dual_rhs = rhs.right_complement();
                dual_lhs.wedge(&dual_rhs).left_complement()
            }

            /// Right complement (`a ∧ rc(a) = +I` per blade).
            #[must_use]
            pub fn right_complement(&self) -> Self {
                let full = ($blades - 1) as u32;
                let mut out = [0.0f64; $blades];
                for (i, &a) in self.0.iter().enumerate() {
                    let t = right_complement(i as u32, full);
                    out[t.blade as usize] += f64::from(t.sign) * a;
                }
                Self(out)
            }

            /// Left complement (inverse of the right complement).
            #[must_use]
            pub fn left_complement(&self) -> Self {
                let full = ($blades - 1) as u32;
                let mut out = [0.0f64; $blades];
                for (i, &a) in self.0.iter().enumerate() {
                    let t = left_complement(i as u32, full);
                    out[t.blade as usize] += f64::from(t.sign) * a;
                }
                Self(out)
            }

            /// Reverse (anti-automorphism: reverses blade factor order).
            #[must_use]
            pub fn reverse(&self) -> Self {
                let mut out = self.0;
                for (i, v) in out.iter_mut().enumerate() {
                    *v *= f64::from(reverse_sign(i as u32));
                }
                Self(out)
            }

            /// Grade involution (negates odd grades).
            #[must_use]
            pub fn involute(&self) -> Self {
                let mut out = self.0;
                for (i, v) in out.iter_mut().enumerate() {
                    *v *= f64::from(involute_sign(i as u32));
                }
                Self(out)
            }

            /// Keep one grade, zero the rest.
            #[must_use]
            pub fn grade_part(&self, g: u32) -> Self {
                let mut out = [0.0f64; $blades];
                for (i, &v) in self.0.iter().enumerate() {
                    if grade(i as u32) == g {
                        out[i] = v;
                    }
                }
                Self(out)
            }

            /// Squared norm: `⟨a ã⟩₀` (can be negative in mixed signature).
            #[must_use]
            pub fn norm_sq(&self) -> f64 {
                self.gp(&self.reverse()).scalar_part()
            }

            /// Multivector sum.
            #[must_use]
            pub fn add(&self, rhs: &Self) -> Self {
                let mut out = self.0;
                for (o, r) in out.iter_mut().zip(rhs.0.iter()) {
                    *o += r;
                }
                Self(out)
            }

            /// Multivector difference.
            #[must_use]
            pub fn sub(&self, rhs: &Self) -> Self {
                let mut out = self.0;
                for (o, r) in out.iter_mut().zip(rhs.0.iter()) {
                    *o -= r;
                }
                Self(out)
            }

            /// Uniform scale.
            #[must_use]
            pub fn scale(&self, s: f64) -> Self {
                let mut out = self.0;
                for o in out.iter_mut() {
                    *o *= s;
                }
                Self(out)
            }

            /// Largest absolute component (∞-norm) — the test metric.
            #[must_use]
            pub fn max_abs(&self) -> f64 {
                let mut m = 0.0f64;
                for &v in &self.0 {
                    if v.abs() > m {
                        m = v.abs();
                    }
                }
                m
            }
        }
    };
}

algebra!(
    /// A dense Cl(3,0,1) (PGA) multivector: 16 blade coefficients indexed
    /// by basis bitmask over (e0, e1, e2, e3), canonical ascending order.
    Pga,
    16,
    crate::table::PGA_TABLE
);

algebra!(
    /// A dense Cl(4,1) (CGA) multivector: 32 blade coefficients indexed
    /// by basis bitmask over (e1, e2, e3, e+, e−), canonical ascending
    /// order.
    Cga,
    32,
    crate::table::CGA_TABLE
);

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn complement_round_trip_is_identity() {
        for i in 0..Pga::BLADES {
            let a = Pga::blade(i, 2.5);
            let back = a.right_complement().left_complement();
            assert_eq!(back, a, "complement round trip broke at blade {i}");
        }
    }

    #[test]
    fn wedge_matches_gp_on_disjoint_blades() {
        // e1 ∧ e2 == e1 e2 (no shared vectors, no metric contraction).
        let e1 = Pga::blade(0b0010, 1.0);
        let e2 = Pga::blade(0b0100, 1.0);
        assert_eq!(e1.wedge(&e2), e1.gp(&e2));
        // e1 ∧ e1 == 0 while e1 e1 == 1.
        assert_eq!(e1.wedge(&e1), Pga::zero());
        assert_eq!(e1.gp(&e1).scalar_part().to_bits(), 1.0f64.to_bits());
    }

    #[test]
    fn reverse_is_an_antiautomorphism() {
        // rev(ab) = rev(b) rev(a) on a couple of dense multivectors.
        let mut a = Pga::zero();
        let mut b = Pga::zero();
        for i in 0..Pga::BLADES {
            a.0[i] = (i as f64).sin_like();
            b.0[i] = (i as f64 + 3.0).sin_like();
        }
        let lhs = a.gp(&b).reverse();
        let rhs = b.reverse().gp(&a.reverse());
        assert!(lhs.sub(&rhs).max_abs() < 1e-12);
    }

    trait SinLike {
        fn sin_like(self) -> f64;
    }
    impl SinLike for f64 {
        fn sin_like(self) -> f64 {
            // Deterministic pseudo-values without pulling in trig.
            let x = (self * 12.9898 + 78.233) % 3.7;
            x - 1.85
        }
    }
}
