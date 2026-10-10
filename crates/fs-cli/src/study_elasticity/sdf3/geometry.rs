//! Physical SI implicit domains. The legacy unit-cube field keeps its arithmetic.
use super::*;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(super) enum FixedFace {
    Left,
    Right,
    Front,
    Back,
    Bottom,
}
impl FixedFace {
    pub(super) fn contains(self, p: [f64; 3], bounds: ([f64; 3], [f64; 3])) -> bool {
        let (axis, upper) = match self {
            Self::Left => (0, false),
            Self::Right => (0, true),
            Self::Front => (1, false),
            Self::Back => (1, true),
            Self::Bottom => (2, false),
        };
        let (lo, hi) = bounds;
        // Octree endpoint interpolation may round lo + (hi - lo) differently
        // from hi. Both name the same terminal lattice plane, not a finite band.
        if upper { p[axis] == hi[axis] || p[axis] == lo[axis] + (hi[axis] - lo[axis]) }
        else { p[axis] == lo[axis] }
    }
}

pub(super) struct PhysicalDomain {
    pub bounds: ([f64; 3], [f64; 3]),
    pub height: f64,
    pub curvature: f64,
}
impl CutSdf3 for PhysicalDomain {
    fn value(&self, p: [f64; 3]) -> f64 {
        let (lo, hi) = self.bounds;
        p[2] - lo[2] - self.height
            - self.curvature * (p[0] - lo[0]) * (hi[0] - p[0])
    }
    fn enclose(&self, lo: [f64; 3], hi: [f64; 3]) -> Interval {
        let point = |x| Interval::new(x, x);
        let x = Interval::new(lo[0], hi[0]);
        let a = point(self.bounds.0[0]);
        let b = point(self.bounds.1[0]);
        let at = |x: Interval| point(self.curvature) * (x - a) * (b - x);
        let curve = if lo[0] == hi[0] {
            at(x)
        } else {
            // Repeated x in the natural interval extension invents material in
            // cells above a curved face. A certified slope sign places extrema
            // at the endpoints; boxes crossing the apex keep the natural bound.
            let slope = point(self.curvature) * ((b - x) - (x - a));
            if slope.lo() >= 0.0 {
                Interval::new(at(point(lo[0])).lo(), at(point(hi[0])).hi())
            } else if slope.hi() <= 0.0 {
                Interval::new(at(point(hi[0])).lo(), at(point(lo[0])).hi())
            } else {
                at(x)
            }
        };
        Interval::new(lo[2], hi[2]) - point(self.bounds.0[2]) - point(self.height) - curve
    }
    fn derivative_enclose(&self, lo: [f64; 3], hi: [f64; 3], axis: HeightAxis) -> Interval {
        let point = |x| Interval::new(x, x);
        match axis {
            HeightAxis::X => point(-self.curvature)
                * ((point(self.bounds.1[0]) - Interval::new(lo[0], hi[0]))
                    - (Interval::new(lo[0], hi[0]) - point(self.bounds.0[0]))),
            HeightAxis::Y => point(0.0),
            HeightAxis::Z => point(1.0),
        }
    }
}

pub(super) fn validate(spec: &Spec) -> Result<()> {
    let invalid = |what| fail("cli-study-sdf3-input", what);
    if !spec.physical {
        if spec.bounds != ([0.0; 3], [1.0; 3]) || spec.fixed != FixedFace::Left
            || !(0.1..=0.8).contains(&spec.height)
            || !(0.0..=0.4).contains(&spec.curvature)
            || !(0.001..=1.0).contains(&spec.radius)
        {
            return Err(invalid("legacy curved-height-sdf requires its original unit cube, left clamp and parameter envelope"));
        }
        return Ok(());
    }
    let (lo, hi) = spec.bounds;
    let spans: [f64; 3] = std::array::from_fn(|a| hi[a] - lo[a]);
    for a in 0..3 {
        if !lo[a].is_finite() || !hi[a].is_finite()
            || !(1e-6..=1000.0).contains(&spans[a])
            || lo[a].abs().max(hi[a].abs()) > spans[a] * 1e8
        {
            return Err(invalid("physical bounds require ordered finite spans in [1e-6,1000] m and resolvable offsets"));
        }
    }
    let shortest = spans.iter().copied().fold(f64::INFINITY, f64::min);
    let longest = spans.iter().copied().fold(0.0_f64, f64::max);
    if longest / shortest > 64.0
        || !(0.001 * shortest..=longest).contains(&spec.radius)
    {
        return Err(invalid("physical domain requires aspect <=64 and a positive filter radius within the domain envelope"));
    }
    if let Some(domain) = &spec.constructive {
        let range = domain.enclose(lo, hi);
        if !range.lo().is_finite() || !range.hi().is_finite() {
            return Err(invalid("constructive field overflows on the declared physical bounds"));
        }
        if range.lo() >= 0.0 {
            return Err(invalid("constructive field proves no material interior in the declared box"));
        }
        // A sign range does not prove connectivity, support, or resolved cuts.
        // The unchanged bounded quadrature and elasticity owner admit those.
    } else if !(0.1..=0.8).contains(&(spec.height / spans[2]))
        || !(0.0..=0.4).contains(&(spec.curvature * spans[0] * spans[0] / spans[2]))
    {
        return Err(invalid("physical graph requires height/Lz in [0.1,0.8] and curvature*Lx^2/Lz in [0,0.4]"));
    }
    Ok(())
}

#[cfg(test)]
#[test]
fn physical_graph_certifies_the_empty_cell_above_the_pressure_fixture() {
    let domain = PhysicalDomain {
        bounds: ([0.0; 3], [0.1, 0.05, 0.04]), height: 0.028, curvature: 1.0,
    };
    let range = domain.enclose([0.0, 0.0, 0.03], [0.025, 0.0125, 0.04]);
    assert!(range.lo() > 0.0, "the true face stays below z=0.03 in this cell");
    let apex = domain.enclose([0.04, 0.0, 0.03], [0.06, 0.0125, 0.04]);
    assert!(apex.lo() < 0.0 && apex.hi() > 0.0, "do not remove the interior apex");
    for i in 0..=16 {
        let x = 0.1 * i as f64 / 16.0;
        let lo = [x, 0.0, 0.02];
        let hi = [(x + 0.01).min(0.1), 0.05, 0.04];
        let range = domain.enclose(lo, hi);
        for px in [lo[0], lo[0].midpoint(hi[0]), hi[0]] {
            for z in [lo[2], hi[2]] {
                let value = domain.value([px, 0.02, z]);
                assert!(range.lo() <= value && value <= range.hi());
            }
        }
    }
}

#[cfg(test)]
#[path = "geometry_tests.rs"]
mod tests;
