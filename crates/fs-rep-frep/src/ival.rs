//! The auto-derived INTERVAL evaluator: `Frep::interval(box)` returns a
//! range guaranteed to contain `f(p)` for every `p` in the box (the G0
//! containment law, frep-001). Every primitive and transform uses conservative
//! outward interval arithmetic. Rotations enclose Rodrigues directly; they do
//! not rely on a rounded corner AABB or a platform-libm ULP assumption.
//! Booleans use monotonicity: `min`/`smin` are nondecreasing in both
//! arguments, so endpoint evaluation is an inclusion. The minimal local
//! interval kit rounds every arithmetic endpoint outward; rotation trig reuses
//! fs-ivl's deterministic fs-math ULP budgets.

use crate::{BoolStyle, Frep, Node, NodeId, bool_signs};
use fs_geom::{Aabb, Point3, Vec3};
use fs_ivl::Interval as CertifiedInterval;

/// Closed interval `[lo, hi]`.
#[derive(Debug, Clone, Copy)]
pub(crate) struct Iv {
    pub lo: f64,
    pub hi: f64,
}

impl Iv {
    const WHOLE: Iv = Iv {
        lo: f64::NEG_INFINITY,
        hi: f64::INFINITY,
    };

    fn new(lo: f64, hi: f64) -> Iv {
        Iv { lo, hi }
    }

    fn point(value: f64) -> Iv {
        Iv::new(value, value)
    }

    fn down(value: f64) -> f64 {
        if value.is_finite() {
            value.next_down()
        } else {
            value
        }
    }

    fn up(value: f64) -> f64 {
        if value.is_finite() {
            value.next_up()
        } else {
            value
        }
    }

    fn outward(lo: f64, hi: f64) -> Iv {
        if lo.is_nan() || hi.is_nan() {
            Iv::WHOLE
        } else {
            Iv::new(Iv::down(lo), Iv::up(hi))
        }
    }

    fn add(self, other: Iv) -> Iv {
        Iv::outward(self.lo + other.lo, self.hi + other.hi)
    }

    fn sub(self, other: Iv) -> Iv {
        Iv::outward(self.lo - other.hi, self.hi - other.lo)
    }

    fn mul(self, other: Iv) -> Iv {
        let products = [
            self.lo * other.lo,
            self.lo * other.hi,
            self.hi * other.lo,
            self.hi * other.hi,
        ];
        if products.iter().any(|value| value.is_nan()) {
            return Iv::WHOLE;
        }
        let lo = products.iter().copied().fold(f64::INFINITY, f64::min);
        let hi = products.iter().copied().fold(f64::NEG_INFINITY, f64::max);
        Iv::outward(lo, hi)
    }

    fn add_c(self, c: f64) -> Iv {
        self.add(Iv::point(c))
    }

    fn neg(self) -> Iv {
        Iv::new(-self.hi, -self.lo)
    }

    fn scale_pos(self, s: f64) -> Iv {
        Iv::outward(self.lo * s, self.hi * s)
    }

    fn div_pos(self, divisor: f64) -> Iv {
        Iv::outward(self.lo / divisor, self.hi / divisor)
    }

    fn contains_zero(self) -> bool {
        self.lo <= 0.0 && 0.0 <= self.hi
    }

    fn div(self, other: Iv) -> Iv {
        if other.contains_zero() {
            return Iv::WHOLE;
        }
        let quotients = [
            self.lo / other.lo,
            self.lo / other.hi,
            self.hi / other.lo,
            self.hi / other.hi,
        ];
        if quotients.iter().any(|value| value.is_nan()) {
            return Iv::WHOLE;
        }
        let lo = quotients.iter().copied().fold(f64::INFINITY, f64::min);
        let hi = quotients.iter().copied().fold(f64::NEG_INFINITY, f64::max);
        Iv::outward(lo, hi)
    }

    fn is_finite(self) -> bool {
        self.lo.is_finite() && self.hi.is_finite()
    }

    fn abs(self) -> Iv {
        if self.lo >= 0.0 {
            self
        } else if self.hi <= 0.0 {
            self.neg()
        } else {
            Iv::new(0.0, self.hi.max(-self.lo))
        }
    }

    fn sq(self) -> Iv {
        let a = self.abs();
        let lo = if a.lo == 0.0 {
            0.0
        } else {
            Iv::down(a.lo * a.lo).max(0.0)
        };
        Iv::new(lo, Iv::up(a.hi * a.hi))
    }

    fn sqrt(self) -> Iv {
        let lo = if self.lo <= 0.0 {
            0.0
        } else {
            Iv::down(self.lo.sqrt()).max(0.0)
        };
        Iv::new(lo, Iv::up(self.hi.max(0.0).sqrt()))
    }

    fn max_c(self, c: f64) -> Iv {
        Iv::new(self.lo.max(c), self.hi.max(c))
    }

    fn min_c(self, c: f64) -> Iv {
        Iv::new(self.lo.min(c), self.hi.min(c))
    }

    fn max_iv(self, o: Iv) -> Iv {
        Iv::new(self.lo.max(o.lo), self.hi.max(o.hi))
    }

    fn min_iv(self, o: Iv) -> Iv {
        Iv::new(self.lo.min(o.lo), self.hi.min(o.hi))
    }

    /// `smin` is nondecreasing in both arguments (its partials are the
    /// convex weights), so endpoint evaluation is an inclusion.
    fn smin_iv(self, o: Iv, r: f64) -> Iv {
        fn point_smin(a: f64, b: f64, r: f64) -> Iv {
            let a = Iv::point(a);
            let b = Iv::point(b);
            let h = Iv::point(r).sub(a.sub(b).abs()).max_c(0.0).div_pos(r);
            let correction = h.sq().scale_pos(r).scale_pos(0.25);
            a.min_iv(b).sub(correction)
        }
        let lower = point_smin(self.lo, o.lo, r);
        let upper = point_smin(self.hi, o.hi, r);
        Iv::new(lower.lo, upper.hi)
    }

    /// `hypot`-style √(a² + b²) inclusion.
    fn hypot_iv(self, o: Iv) -> Iv {
        self.sq().add(o.sq()).sqrt()
    }
}

/// Component intervals of `p − c` for `p` in the box.
fn delta_iv(b: &Aabb, c: Point3) -> [Iv; 3] {
    [
        Iv::new(b.min.x, b.max.x).sub(Iv::point(c.x)),
        Iv::new(b.min.y, b.max.y).sub(Iv::point(c.y)),
        Iv::new(b.min.z, b.max.z).sub(Iv::point(c.z)),
    ]
}

/// Outward-rounded `|p - c|` range over a box.
fn dist_iv(b: &Aabb, c: Point3) -> Iv {
    let [x, y, z] = delta_iv(b, c);
    x.sq().add(y.sq()).add(z.sq()).sqrt()
}

fn rotation_trig(angle: f64) -> (Iv, Iv) {
    if angle.is_finite() {
        let angle = CertifiedInterval::point(angle);
        let sine = angle.sin();
        let cosine = angle.cos();
        (
            Iv::new(sine.lo(), sine.hi()),
            Iv::new(cosine.lo(), cosine.hi()),
        )
    } else {
        (Iv::new(-1.0, 1.0), Iv::new(-1.0, 1.0))
    }
}

/// Interval matrix for the real Rodrigues map represented by `axis, angle`.
/// The coefficient intervals share the same fs-ivl trig enclosure as field
/// evaluation; treating their dependencies independently only widens results.
fn rodrigues_matrix(axis: Vec3, angle: f64) -> [[Iv; 3]; 3] {
    let (sine, cosine) = rotation_trig(angle);
    let one_minus_cosine = Iv::point(1.0).sub(cosine);
    let zero = Iv::point(0.0);
    let a = [Iv::point(axis.x), Iv::point(axis.y), Iv::point(axis.z)];
    let cross = [
        [zero, a[2].neg(), a[1]],
        [a[2], zero, a[0].neg()],
        [a[1].neg(), a[0], zero],
    ];
    core::array::from_fn(|row| {
        core::array::from_fn(|column| {
            let diagonal = if row == column { cosine } else { zero };
            diagonal
                .add(cross[row][column].mul(sine))
                .add(a[row].mul(a[column]).mul(one_minus_cosine))
        })
    })
}

fn inverse_3x3(matrix: [[Iv; 3]; 3]) -> Option<[[Iv; 3]; 3]> {
    let cofactor_00 = matrix[1][1]
        .mul(matrix[2][2])
        .sub(matrix[1][2].mul(matrix[2][1]));
    let cofactor_01 = matrix[1][2]
        .mul(matrix[2][0])
        .sub(matrix[1][0].mul(matrix[2][2]));
    let cofactor_02 = matrix[1][0]
        .mul(matrix[2][1])
        .sub(matrix[1][1].mul(matrix[2][0]));
    let determinant = matrix[0][0]
        .mul(cofactor_00)
        .add(matrix[0][1].mul(cofactor_01))
        .add(matrix[0][2].mul(cofactor_02));
    if determinant.contains_zero() || !determinant.is_finite() {
        return None;
    }

    let adjugate = [
        [
            cofactor_00,
            matrix[0][2]
                .mul(matrix[2][1])
                .sub(matrix[0][1].mul(matrix[2][2])),
            matrix[0][1]
                .mul(matrix[1][2])
                .sub(matrix[0][2].mul(matrix[1][1])),
        ],
        [
            cofactor_01,
            matrix[0][0]
                .mul(matrix[2][2])
                .sub(matrix[0][2].mul(matrix[2][0])),
            matrix[0][2]
                .mul(matrix[1][0])
                .sub(matrix[0][0].mul(matrix[1][2])),
        ],
        [
            cofactor_02,
            matrix[0][1]
                .mul(matrix[2][0])
                .sub(matrix[0][0].mul(matrix[2][1])),
            matrix[0][0]
                .mul(matrix[1][1])
                .sub(matrix[0][1].mul(matrix[1][0])),
        ],
    ];
    let inverse = adjugate.map(|row| row.map(|entry| entry.div(determinant)));
    inverse
        .iter()
        .flatten()
        .all(|entry| entry.is_finite())
        .then_some(inverse)
}

fn infinite_aabb() -> Aabb {
    Aabb::WHOLE_SPACE
}

/// Enclose the preimage `{p | R(-angle) p in child}` used by `Node::Rotate`.
/// Inverting the certified inverse-map matrix avoids assuming that rounded
/// Rodrigues coefficients remain exactly orthogonal. If regularity cannot be
/// certified, the only honest support is the whole space.
#[allow(clippy::float_cmp)] // Exact zero is the identity map, including -0.0.
pub(crate) fn rotation_preimage_support(child: &Aabb, axis: Vec3, angle: f64) -> Aabb {
    if angle == 0.0 {
        return *child;
    }
    let inverse_map = rodrigues_matrix(axis, -angle);
    let Some(preimage_map) = inverse_3x3(inverse_map) else {
        return infinite_aabb();
    };
    let child_coordinates = [
        Iv::new(child.min.x, child.max.x),
        Iv::new(child.min.y, child.max.y),
        Iv::new(child.min.z, child.max.z),
    ];
    let preimage: [Iv; 3] = core::array::from_fn(|row| {
        preimage_map[row][0]
            .mul(child_coordinates[0])
            .add(preimage_map[row][1].mul(child_coordinates[1]))
            .add(preimage_map[row][2].mul(child_coordinates[2]))
    });
    if !preimage.iter().all(|coordinate| coordinate.is_finite()) {
        return infinite_aabb();
    }
    Aabb::new(
        Point3::new(preimage[0].lo, preimage[1].lo, preimage[2].lo),
        Point3::new(preimage[0].hi, preimage[1].hi, preimage[2].hi),
    )
}

impl Frep {
    /// Range guaranteed to contain `f(p)` for all `p ∈ region`.
    #[must_use]
    pub fn interval(&self, region: &Aabb) -> (f64, f64) {
        let iv = self.iv_at(self.root(), region);
        (iv.lo, iv.hi)
    }

    #[allow(clippy::too_many_lines)] // One exhaustive evaluator keeps every Node enclosure rule co-located.
    fn iv_at(&self, id: NodeId, b: &Aabb) -> Iv {
        match self.nodes()[id.0 as usize] {
            Node::Sphere { center, radius } => dist_iv(b, center).add_c(-radius),
            Node::HalfSpace { normal, offset } => {
                let mut value = Iv::point(-offset);
                for (n, bmin, bmax) in [
                    (normal.x, b.min.x, b.max.x),
                    (normal.y, b.min.y, b.max.y),
                    (normal.z, b.min.z, b.max.z),
                ] {
                    let coordinate = Iv::new(bmin, bmax);
                    let product = if n >= 0.0 {
                        coordinate.scale_pos(n)
                    } else {
                        coordinate.neg().scale_pos(-n)
                    };
                    value = value.add(product);
                }
                value
            }
            Node::BoxPrim { center, half } => {
                let d = delta_iv(b, center);
                let q = [
                    d[0].abs().add_c(-half.x),
                    d[1].abs().add_c(-half.y),
                    d[2].abs().add_c(-half.z),
                ];
                let out = [q[0].max_c(0.0), q[1].max_c(0.0), q[2].max_c(0.0)];
                let norm = out[0].sq().add(out[1].sq()).add(out[2].sq()).sqrt();
                let inner = q[0].max_iv(q[1]).max_iv(q[2]).min_c(0.0);
                norm.add(inner)
            }
            Node::Torus {
                center,
                major,
                minor,
            } => {
                let d = delta_iv(b, center);
                let ring = d[0].hypot_iv(d[1]).add_c(-major);
                ring.hypot_iv(d[2]).add_c(-minor)
            }
            Node::Cylinder { center, radius } => {
                let d = delta_iv(b, center);
                d[0].hypot_iv(d[1]).add_c(-radius)
            }
            Node::Translate { child, offset } => {
                let x = Iv::new(b.min.x, b.max.x).add_c(-offset.x);
                let y = Iv::new(b.min.y, b.max.y).add_c(-offset.y);
                let z = Iv::new(b.min.z, b.max.z).add_c(-offset.z);
                let shifted =
                    Aabb::new(Point3::new(x.lo, y.lo, z.lo), Point3::new(x.hi, y.hi, z.hi));
                self.iv_at(child, &shifted)
            }
            Node::Rotate { child, axis, angle } => {
                let v = [
                    Iv::new(b.min.x, b.max.x),
                    Iv::new(b.min.y, b.max.y),
                    Iv::new(b.min.z, b.max.z),
                ];
                let a = [Iv::point(axis.x), Iv::point(axis.y), Iv::point(axis.z)];
                // Evaluate the inverse rotation with the same trig enclosure
                // used by the support preimage proof.
                let (sine, cosine) = rotation_trig(-angle);
                let one_minus_cosine = Iv::point(1.0).sub(cosine);
                let cross = [
                    a[1].mul(v[2]).sub(a[2].mul(v[1])),
                    a[2].mul(v[0]).sub(a[0].mul(v[2])),
                    a[0].mul(v[1]).sub(a[1].mul(v[0])),
                ];
                let dot = a[0].mul(v[0]).add(a[1].mul(v[1])).add(a[2].mul(v[2]));
                let rotated: [Iv; 3] = core::array::from_fn(|component| {
                    v[component]
                        .mul(cosine)
                        .add(cross[component].mul(sine))
                        .add(a[component].mul(dot).mul(one_minus_cosine))
                });
                let mapped = Aabb::new(
                    Point3::new(rotated[0].lo, rotated[1].lo, rotated[2].lo),
                    Point3::new(rotated[0].hi, rotated[1].hi, rotated[2].hi),
                );
                self.iv_at(child, &mapped)
            }
            Node::Scale { child, factor } => {
                let x = Iv::new(b.min.x, b.max.x).div_pos(factor);
                let y = Iv::new(b.min.y, b.max.y).div_pos(factor);
                let z = Iv::new(b.min.z, b.max.z).div_pos(factor);
                let shrunk =
                    Aabb::new(Point3::new(x.lo, y.lo, z.lo), Point3::new(x.hi, y.hi, z.hi));
                self.iv_at(child, &shrunk).scale_pos(factor)
            }
            Node::Offset { child, distance } => self.iv_at(child, b).add_c(-distance),
            Node::Bool {
                op,
                style,
                a,
                b: rhs,
            } => {
                let (sa, sb, sr) = bool_signs(op);
                let ia = if sa < 0.0 {
                    self.iv_at(a, b).neg()
                } else {
                    self.iv_at(a, b)
                };
                let ib = if sb < 0.0 {
                    self.iv_at(rhs, b).neg()
                } else {
                    self.iv_at(rhs, b)
                };
                let m = match style {
                    BoolStyle::Hard => ia.min_iv(ib),
                    BoolStyle::Blend { radius } => ia.smin_iv(ib, radius),
                };
                if sr < 0.0 { m.neg() } else { m }
            }
        }
    }
}
