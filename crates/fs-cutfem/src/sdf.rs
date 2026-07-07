//! The level-set interface contract: a [`CutSdf`] is a scalar field φ
//! (negative inside, positive outside) that can CERTIFY its range over
//! an axis-aligned box via fs-ivl outward-rounded intervals. The
//! enclosure is what makes cut classification sound: a cell is labeled
//! Inside/Outside only when the certified range excludes zero, so no
//! cell is ever misclassified — tangency and near-tangency land in the
//! conservative Cut class (the G0 battery's adversarial cases).

use fs_ivl::Interval;

/// A level-set function with certified box enclosures.
pub trait CutSdf {
    /// φ at a point (negative = inside the domain Ω).
    fn value(&self, p: [f64; 2]) -> f64;
    /// ∇φ at a point (used for interface normals; need not be unit —
    /// consumers normalize).
    fn gradient(&self, p: [f64; 2]) -> [f64; 2];
    /// A certified enclosure of φ over the box `[lo, hi]`: the returned
    /// interval CONTAINS φ(p) for every p in the box (the containment
    /// law — conformance-tested against dense sampling).
    fn enclose(&self, lo: [f64; 2], hi: [f64; 2]) -> Interval;
}

/// Euclidean circle: φ = |p − c| − r (an exact SDF).
#[derive(Debug, Clone, Copy)]
pub struct Circle {
    /// Center.
    pub center: [f64; 2],
    /// Radius (> 0).
    pub radius: f64,
}

impl CutSdf for Circle {
    fn value(&self, p: [f64; 2]) -> f64 {
        let dx = p[0] - self.center[0];
        let dy = p[1] - self.center[1];
        (dx * dx + dy * dy).sqrt() - self.radius
    }

    fn gradient(&self, p: [f64; 2]) -> [f64; 2] {
        let dx = p[0] - self.center[0];
        let dy = p[1] - self.center[1];
        let n = (dx * dx + dy * dy).sqrt();
        if n < 1e-300 {
            // The center is a gradient singularity; any unit vector is
            // a valid subgradient direction for our consumers.
            return [1.0, 0.0];
        }
        [dx / n, dy / n]
    }

    fn enclose(&self, lo: [f64; 2], hi: [f64; 2]) -> Interval {
        let dx = (Interval::new(lo[0], hi[0]) - Interval::point(self.center[0])).abs();
        let dy = (Interval::new(lo[1], hi[1]) - Interval::point(self.center[1])).abs();
        let d2 = dx * dx + dy * dy;
        d2.sqrt() - Interval::point(self.radius)
    }
}

/// Half-plane: φ = n·p − offset (exact SDF when `normal` is unit; a
/// valid level set for any nonzero `normal`).
#[derive(Debug, Clone, Copy)]
pub struct HalfPlane {
    /// Outward normal (need not be unit).
    pub normal: [f64; 2],
    /// Signed offset: the zero set is `{ p : normal·p = offset }`.
    pub offset: f64,
}

impl CutSdf for HalfPlane {
    fn value(&self, p: [f64; 2]) -> f64 {
        self.normal[0] * p[0] + self.normal[1] * p[1] - self.offset
    }

    fn gradient(&self, _p: [f64; 2]) -> [f64; 2] {
        self.normal
    }

    fn enclose(&self, lo: [f64; 2], hi: [f64; 2]) -> Interval {
        Interval::point(self.normal[0]) * Interval::new(lo[0], hi[0])
            + Interval::point(self.normal[1]) * Interval::new(lo[1], hi[1])
            - Interval::point(self.offset)
    }
}
