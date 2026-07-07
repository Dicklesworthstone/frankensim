//! Plane-based PGA Cl(3,0,1) (plan §7.7): points/lines/planes as blades,
//! MOTORS for rigid motion (rotation + translation in one even-grade
//! versor — no gimbal coordinates anywhere), join (∨) and meet (∧) for
//! incidence, and closed-form screw exp/log for interpolation.
//!
//! Conventions (ganja/PGA-course component layout on canonical blades):
//! plane `a·x+b·y+c·z+d = 0` ↔ `a e1 + b e2 + c e3 + d e0`; point
//! `(x,y,z)` ↔ `e123 + x e032 + y e013 + z e021`. Trig comes from
//! `fs_math::det` so motors are bit-identical across platforms.

use crate::GaError;
use crate::mv::Pga;
use crate::table::{PGA_TABLE_CONST, grade};
use fs_math::det;

/// Blade index of the PGA pseudoscalar e0123.
pub const E0123: usize = 15;
const E12: usize = 0b0110;
const E13: usize = 0b1010; // e31 = −e13
const E23: usize = 0b1100;
const E01: usize = 0b0011;
const E02: usize = 0b0101;
const E03: usize = 0b1001;
const E012: usize = 0b0111;
const E013: usize = 0b1011;
const E023: usize = 0b1101;
const E123: usize = 0b1110;

/// Angle threshold under which a screw is treated as a pure translation
/// (the exp/log series is exact to f64 below it).
const TINY_ANGLE: f64 = 1e-9;

/// A Euclidean point (homogeneous trivector).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Point {
    /// Cartesian x.
    pub x: f64,
    /// Cartesian y.
    pub y: f64,
    /// Cartesian z.
    pub z: f64,
}

impl Point {
    /// Construct from Cartesian coordinates.
    #[must_use]
    pub fn new(x: f64, y: f64, z: f64) -> Self {
        Point { x, y, z }
    }

    /// The unit-weight trivector encoding.
    #[must_use]
    pub fn to_mv(self) -> Pga {
        let mut v = Pga::zero();
        v.0[E123] = 1.0;
        v.0[E023] = -self.x;
        v.0[E013] = self.y;
        v.0[E012] = -self.z;
        v
    }

    /// Decode a trivector back to Cartesian coordinates (dividing out the
    /// homogeneous weight).
    ///
    /// # Errors
    /// [`GaError::IdealPoint`] when the weight (e123 coefficient) is zero
    /// — an ideal (at-infinity) point has no Cartesian form.
    pub fn from_mv(v: &Pga) -> Result<Self, GaError> {
        let w = v.0[E123];
        if w == 0.0 {
            return Err(GaError::IdealPoint);
        }
        Ok(Point {
            x: -v.0[E023] / w,
            y: v.0[E013] / w,
            z: -v.0[E012] / w,
        })
    }
}

/// An oriented plane `a·x + b·y + c·z + d = 0` (grade-1 vector).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Plane {
    /// Normal x component.
    pub a: f64,
    /// Normal y component.
    pub b: f64,
    /// Normal z component.
    pub c: f64,
    /// Offset.
    pub d: f64,
}

impl Plane {
    /// Construct from the implicit equation coefficients.
    #[must_use]
    pub fn new(a: f64, b: f64, c: f64, d: f64) -> Self {
        Plane { a, b, c, d }
    }

    /// The grade-1 encoding.
    #[must_use]
    pub fn to_mv(self) -> Pga {
        let mut v = Pga::zero();
        v.0[0b0010] = self.a;
        v.0[0b0100] = self.b;
        v.0[0b1000] = self.c;
        v.0[0b0001] = self.d;
        v
    }

    /// Decode a grade-1 multivector.
    #[must_use]
    pub fn from_mv(v: &Pga) -> Self {
        Plane {
            a: v.0[0b0010],
            b: v.0[0b0100],
            c: v.0[0b1000],
            d: v.0[0b0001],
        }
    }

    /// Signed incidence measure `a·x + b·y + c·z + d` for a point,
    /// computed as the (plane ∧ point) pseudoscalar coefficient.
    #[must_use]
    pub fn incidence(&self, p: Point) -> f64 {
        self.to_mv().wedge(&p.to_mv()).0[E0123]
    }
}

/// A line as a grade-2 blade: Euclidean part = direction bivector,
/// ideal part = moment. Built by joining points or meeting planes.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Line(pub Pga);

impl Line {
    /// Join: the line through two points (regressive product).
    #[must_use]
    pub fn through(p: Point, q: Point) -> Self {
        Line(p.to_mv().vee(&q.to_mv()).grade_part(2))
    }

    /// Meet: the intersection line of two planes (outer product).
    #[must_use]
    pub fn meet(p: Plane, q: Plane) -> Self {
        Line(p.to_mv().wedge(&q.to_mv()).grade_part(2))
    }

    /// ∞-norm of the join with a point: zero iff the point lies on the
    /// line.
    #[must_use]
    pub fn incidence(&self, p: Point) -> f64 {
        self.0.vee(&p.to_mv()).max_abs()
    }
}

/// Build the Euclidean rotation bivector for a (not necessarily unit)
/// axis direction through the origin: `nx e23 + ny e31 + nz e12`.
#[must_use]
pub fn axis_bivector(nx: f64, ny: f64, nz: f64) -> Pga {
    let mut v = Pga::zero();
    v.0[E23] = nx;
    v.0[E13] = -ny; // e31 stored on canonical e13 with flipped sign
    v.0[E12] = nz;
    v
}

/// Build the ideal (translation) bivector `tx e01 + ty e02 + tz e03`.
#[must_use]
pub fn ideal_bivector(tx: f64, ty: f64, tz: f64) -> Pga {
    let mut v = Pga::zero();
    v.0[E01] = tx;
    v.0[E02] = ty;
    v.0[E03] = tz;
    v
}

/// Even-subalgebra blade indices in ascending order — the eight motor
/// components (scalar, e01, e02, e12, e03, e13, e23, e0123).
pub const EVEN_BLADES: [usize; 8] = [0, 3, 5, 6, 9, 10, 12, 15];
/// Odd blade indices (grades 1 and 3) — motor × point intermediates.
const ODD_BLADES: [usize; 8] = [1, 2, 4, 7, 8, 11, 13, 14];
/// Grade-3 (point trivector) blade indices.
const POINT_BLADES: [usize; 4] = [7, 11, 13, 14];

/// A compact table entry: `sign * component[idx]` (sign 0 ⇒ dead term).
#[derive(Debug, Clone, Copy)]
struct KTerm {
    sign: i8,
    idx: u8,
}

const fn position(list: &[usize], blade: usize) -> u8 {
    let mut k = 0;
    while k < list.len() {
        if list[k] == blade {
            #[allow(clippy::cast_possible_truncation)] // lists are tiny
            return k as u8;
        }
        k += 1;
    }
    255
}

/// Monomorphized even ⊗ even kernel (motor composition): generated at
/// compile time from the Cayley table — 64 fused sign/index terms, no
/// runtime blade bookkeeping.
static EVEN_GP: [[KTerm; 8]; 8] = build_even_gp();

const fn build_even_gp() -> [[KTerm; 8]; 8] {
    let mut t = [[KTerm { sign: 0, idx: 0 }; 8]; 8];
    let mut i = 0;
    while i < 8 {
        let mut j = 0;
        while j < 8 {
            let term = PGA_TABLE_CONST[EVEN_BLADES[i]][EVEN_BLADES[j]];
            t[i][j] = KTerm {
                sign: term.sign,
                idx: position(&EVEN_BLADES, term.blade as usize),
            };
            j += 1;
        }
        i += 1;
    }
    t
}

/// Monomorphized even ⊗ point kernel (first half of the sandwich).
static EVEN_POINT_GP: [[KTerm; 4]; 8] = build_even_point_gp();

const fn build_even_point_gp() -> [[KTerm; 4]; 8] {
    let mut t = [[KTerm { sign: 0, idx: 0 }; 4]; 8];
    let mut i = 0;
    while i < 8 {
        let mut j = 0;
        while j < 4 {
            let term = PGA_TABLE_CONST[EVEN_BLADES[i]][POINT_BLADES[j]];
            t[i][j] = KTerm {
                sign: term.sign,
                idx: position(&ODD_BLADES, term.blade as usize),
            };
            j += 1;
        }
        i += 1;
    }
    t
}

/// Monomorphized odd ⊗ even kernel restricted to grade-3 output (second
/// half of the sandwich; non-point targets are dead terms).
static ODD_EVEN_TO_POINT_GP: [[KTerm; 8]; 8] = build_odd_even_gp();

const fn build_odd_even_gp() -> [[KTerm; 8]; 8] {
    let mut t = [[KTerm { sign: 0, idx: 0 }; 8]; 8];
    let mut i = 0;
    while i < 8 {
        let mut j = 0;
        while j < 8 {
            let term = PGA_TABLE_CONST[ODD_BLADES[i]][EVEN_BLADES[j]];
            if grade(term.blade as u32) == 3 {
                t[i][j] = KTerm {
                    sign: term.sign,
                    idx: position(&POINT_BLADES, term.blade as usize),
                };
            }
            j += 1;
        }
        i += 1;
    }
    t
}

/// A motor: an even-grade versor encoding a full rigid motion (screw).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Motor(pub Pga);

impl Motor {
    /// The identity motion.
    #[must_use]
    pub fn identity() -> Self {
        Motor(Pga::scalar(1.0))
    }

    /// Rotation by `angle` (radians) about the unit `axis` through the
    /// origin.
    #[must_use]
    pub fn rotor(axis: [f64; 3], angle: f64) -> Self {
        exp_bivector(&axis_bivector(axis[0], axis[1], axis[2]).scale(-angle / 2.0))
    }

    /// Translation by `(tx, ty, tz)`.
    #[must_use]
    pub fn translator(tx: f64, ty: f64, tz: f64) -> Self {
        exp_bivector(&ideal_bivector(tx, ty, tz).scale(-0.5))
    }

    /// The eight even-subalgebra components in [`EVEN_BLADES`] order.
    #[must_use]
    fn even(&self) -> [f64; 8] {
        let mut e = [0.0f64; 8];
        for (out, &blade) in e.iter_mut().zip(EVEN_BLADES.iter()) {
            *out = self.0.0[blade];
        }
        e
    }

    fn from_even(e: [f64; 8]) -> Motor {
        let mut v = Pga::zero();
        for (&val, &blade) in e.iter().zip(EVEN_BLADES.iter()) {
            v.0[blade] = val;
        }
        Motor(v)
    }

    /// Compose: `self` applied AFTER `rhs` (like matrix multiplication).
    /// Runs on the monomorphized even ⊗ even kernel (64 fused terms).
    #[must_use]
    pub fn compose(&self, rhs: &Motor) -> Self {
        let a = self.even();
        let b = rhs.even();
        let mut out = [0.0f64; 8];
        for (i, &ai) in a.iter().enumerate() {
            for (j, &bj) in b.iter().enumerate() {
                let t = EVEN_GP[i][j];
                out[t.idx as usize] += f64::from(t.sign) * ai * bj;
            }
        }
        Motor::from_even(out)
    }

    /// Reverse — the inverse for unit motors.
    #[must_use]
    pub fn reverse(&self) -> Self {
        Motor(self.0.reverse())
    }

    /// Apply the rigid motion to a point (sandwich `M P M̃`) on the
    /// monomorphized even⊗point and odd⊗even kernels.
    ///
    /// # Errors
    /// [`GaError::IdealPoint`] if the result has zero weight (cannot
    /// happen for unit motors on Euclidean points; guards drifted input).
    pub fn transform_point(&self, p: Point) -> Result<Point, GaError> {
        let m = self.even();
        // Reverse signs on the even subalgebra: grades (0,2,2,2,2,2,2,4)
        // in EVEN_BLADES order → (+,−,−,−,−,−,−,+).
        let mrev = [m[0], -m[1], -m[2], -m[3], -m[4], -m[5], -m[6], m[7]];
        // Point trivector components in POINT_BLADES (7, 11, 13, 14) order.
        let pv = [-p.z, p.y, -p.x, 1.0];
        let mut mid = [0.0f64; 8];
        for (i, &mi) in m.iter().enumerate() {
            for (j, &pj) in pv.iter().enumerate() {
                let t = EVEN_POINT_GP[i][j];
                mid[t.idx as usize] += f64::from(t.sign) * mi * pj;
            }
        }
        let mut out = [0.0f64; 4];
        for (i, &qi) in mid.iter().enumerate() {
            for (j, &mj) in mrev.iter().enumerate() {
                let t = ODD_EVEN_TO_POINT_GP[i][j];
                if t.sign != 0 {
                    out[t.idx as usize] += f64::from(t.sign) * qi * mj;
                }
            }
        }
        let w = out[3];
        if w == 0.0 {
            return Err(GaError::IdealPoint);
        }
        Ok(Point {
            x: -out[2] / w,
            y: out[1] / w,
            z: -out[0] / w,
        })
    }

    /// Reference sandwich on the dense 16-component path (the kernel
    /// cross-check; also handles non-motor even elements).
    ///
    /// # Errors
    /// [`GaError::IdealPoint`] on zero-weight results.
    pub fn transform_point_dense(&self, p: Point) -> Result<Point, GaError> {
        let sandwich = self.0.gp(&p.to_mv()).gp(&self.0.reverse());
        Point::from_mv(&sandwich.grade_part(3))
    }

    /// Apply the rigid motion to a plane.
    #[must_use]
    pub fn transform_plane(&self, p: Plane) -> Plane {
        let sandwich = self.0.gp(&p.to_mv()).gp(&self.0.reverse());
        Plane::from_mv(&sandwich.grade_part(1))
    }

    /// Versor-renormalization policy: divides out the full `M M̃ = a + b I`
    /// residue (both the scalar drift and the pseudoscalar drift), so long
    /// products stay unit to first order. Returns the drift magnitude
    /// `max(|a−1|, |b|)` for ledger statistics.
    #[must_use]
    pub fn renormalize(&mut self) -> f64 {
        let mm = self.0.gp(&self.0.reverse());
        let a = mm.scalar_part();
        let b = mm.0[E0123];
        let drift = (a - 1.0).abs().max(b.abs());
        let inv_sqrt_a = 1.0 / det::sqrt(a);
        // (a + bI)^(−1/2) = a^(−1/2) · (1 − (b / 2a) I)   [I² = 0]
        let corr = Pga::scalar(inv_sqrt_a).sub(&Pga::blade(E0123, inv_sqrt_a * b / (2.0 * a)));
        self.0 = self.0.gp(&corr);
        drift
    }

    /// Unit-motor defect `‖M M̃ − 1‖∞` (0 for exact motors).
    #[must_use]
    pub fn unit_defect(&self) -> f64 {
        self.0
            .gp(&self.0.reverse())
            .sub(&Pga::scalar(1.0))
            .max_abs()
    }

    /// Deterministic screw interpolation `M₀ · exp(t · log(M̃₀ M₁))`:
    /// constant-speed rigid interpolation with no gimbal pathologies.
    #[must_use]
    pub fn slerp(&self, other: &Motor, t: f64) -> Motor {
        let rel = self.reverse().compose(other);
        let b = motor_log(&rel);
        self.compose(&exp_bivector(&b.scale(t)))
    }
}

/// Exponential of a PGA bivector via exact screw decomposition
/// `B = θ ℓ + d ℓ*` (ℓ a unit line, ℓ* = ℓI its ideal companion):
/// `exp B = (cos θ + sin θ · ℓ)(1 + d ℓ*)`.
#[must_use]
pub fn exp_bivector(b: &Pga) -> Motor {
    let b = b.grade_part(2);
    let b_sq = b.gp(&b); // = −θ² − 2θd·I  (scalar + pseudoscalar only)
    let s = b_sq.scalar_part();
    let theta = det::sqrt((-s).max(0.0));
    if theta < TINY_ANGLE {
        // Ideal (or vanishing) bivector: (B)² = 0 ⇒ exp is exactly 1 + B.
        return Motor(Pga::scalar(1.0).add(&b));
    }
    let d = -b_sq.0[E0123] / (2.0 * theta);
    // ℓ = B (θ − d I) / θ²   [(θ + dI)⁻¹ = (θ − dI)/θ² since I² = 0]
    let ell = b
        .gp(&Pga::scalar(theta).sub(&Pga::blade(E0123, d)))
        .scale(1.0 / (theta * theta));
    let ell_star = ell.gp(&Pga::blade(E0123, 1.0));
    let rot = Pga::scalar(det::cos(theta)).add(&ell.scale(det::sin(theta)));
    let trans = Pga::scalar(1.0).add(&ell_star.scale(d));
    Motor(rot.gp(&trans))
}

/// Principal logarithm of a unit motor (inverse of [`exp_bivector`];
/// returns the bivector `θ ℓ + d ℓ*`).
#[must_use]
pub fn motor_log(m: &Motor) -> Pga {
    // Double cover: −M is the same motion; take the branch with s₀ ≥ 0.
    let mv = if m.0.scalar_part() < 0.0 {
        m.0.scale(-1.0)
    } else {
        m.0
    };
    let s0 = mv.scalar_part();
    let p = mv.0[E0123];
    let b2 = mv.grade_part(2);
    let ne = det::sqrt(b2.0[E12] * b2.0[E12] + b2.0[E13] * b2.0[E13] + b2.0[E23] * b2.0[E23]);
    if ne < TINY_ANGLE {
        // Pure translator (or identity): M = 1 + B with B ideal.
        return b2.scale(1.0 / s0);
    }
    let theta = det::atan2(ne, s0);
    let (sin_t, cos_t) = (det::sin(theta), det::cos(theta));
    let d = -p / sin_t;
    // Euclidean part of ℓ from the Euclidean part of ⟨M⟩₂ = sinθ·ℓ + d cosθ·ℓ*.
    let mut ell_e = Pga::zero();
    ell_e.0[E12] = b2.0[E12] / sin_t;
    ell_e.0[E13] = b2.0[E13] / sin_t;
    ell_e.0[E23] = b2.0[E23] / sin_t;
    let ell_star = ell_e.gp(&Pga::blade(E0123, 1.0)); // ℓ* = ℓI (ideal only)
    let mut ideal = Pga::zero();
    for idx in [E01, E02, E03] {
        ideal.0[idx] = (b2.0[idx] - d * cos_t * ell_star.0[idx]) / sin_t;
    }
    ell_e.add(&ideal).scale(theta).add(&ell_star.scale(d))
}
