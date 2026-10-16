use super::*;
use fs_ivl::{Sign, orient2d};
use fs_rep_nurbs::{KnotVector, NurbsCurve, Rat, TrimLoop, TrimmedPatch};

/// Measured unclosed polygon loops: outer first, then real openings.
/// Exact rational degree-one trim curves represent these same segments.
#[derive(Debug)]
pub struct OutlineWithHoles {
    pub(super) label: String,
    loops: Vec<Vec<[MeasuredDimension; 2]>>,
    vertices: Vec<Vec<[f64; 2]>>,
    patch: TrimmedPatch,
}
fn signed_area(p: &[[f64; 2]]) -> f64 {
    // Translate before summing to reduce cancellation under rigid translation.
    let o = p[0];
    (1..p.len() - 1)
        .map(|i| (p[i][0] - o[0]) * (p[i + 1][1] - o[1]) - (p[i][1] - o[1]) * (p[i + 1][0] - o[0]))
        .sum::<f64>()
        / 2.0
}
fn on(a: [f64; 2], b: [f64; 2], p: [f64; 2]) -> bool {
    orient2d(a, b, p) == Sign::Zero
        && (0..2).all(|i| p[i] >= a[i].min(b[i]) && p[i] <= a[i].max(b[i]))
}
fn crosses(a: [f64; 2], b: [f64; 2], c: [f64; 2], d: [f64; 2]) -> bool {
    let ab = [orient2d(a, b, c), orient2d(a, b, d)];
    let cd = [orient2d(c, d, a), orient2d(c, d, b)];
    on(a, b, c)
        || on(a, b, d)
        || on(c, d, a)
        || on(c, d, b)
        || (ab[0] != Sign::Zero
            && ab[1] != Sign::Zero
            && cd[0] != Sign::Zero
            && cd[1] != Sign::Zero
            && ab[0] != ab[1]
            && cd[0] != cd[1])
}
fn inside(p: [f64; 2], polygon: &[[f64; 2]]) -> bool {
    let mut winding = 0;
    for i in 0..polygon.len() {
        let a = polygon[i];
        let b = polygon[(i + 1) % polygon.len()];
        if on(a, b, p) {
            return false;
        }
        if a[1] <= p[1] && b[1] > p[1] && orient2d(a, b, p) == Sign::Positive {
            winding += 1;
        }
        if a[1] > p[1] && b[1] <= p[1] && orient2d(a, b, p) == Sign::Negative {
            winding -= 1;
        }
    }
    winding != 0
}
/// Exact binary64 -> bounded i128 rational, with a structured refusal outside
/// the representation owner's exactness domain; no decimal quantization.
fn rational(x: f64) -> Result<Rat, GeneratorError> {
    if x == 0.0 {
        return Ok(Rat::int(0));
    }
    let bits = x.to_bits();
    let exponent = ((bits >> 52) & 2047) as i32;
    let mut mantissa = bits & ((1u64 << 52) - 1);
    if exponent != 0 {
        mantissa |= 1u64 << 52;
    }
    let mut power = if exponent == 0 {
        -1074
    } else {
        exponent - 1023 - 52
    };
    let zeros = mantissa.trailing_zeros();
    mantissa >>= zeros;
    power += zeros as i32;
    let magnitude = i128::from(mantissa);
    let sign = if x < 0.0 { -1 } else { 1 };
    if power >= 0 {
        let factor = 1i128
            .checked_shl(power as u32)
            .filter(|v| *v > 0)
            .ok_or_else(|| bad("trim coordinate exceeds rational domain"))?;
        Ok(Rat::new(
            magnitude
                .checked_mul(factor)
                .ok_or_else(|| bad("trim coordinate exceeds rational domain"))?
                * sign,
            1,
        ))
    } else {
        let den = 1i128
            .checked_shl((-power) as u32)
            .filter(|v| *v > 0)
            .ok_or_else(|| bad("trim coordinate exceeds rational domain"))?;
        Ok(Rat::new(magnitude * sign, den))
    }
}
impl OutlineWithHoles {
    pub fn new(
        label: String,
        outer: Vec<[MeasuredDimension; 2]>,
        holes: Vec<Vec<[MeasuredDimension; 2]>>,
        budget: GeneratorBudget,
    ) -> Result<Self, GeneratorError> {
        label_valid(&label)?;
        let count = holes
            .iter()
            .try_fold(outer.len(), |a, h| a.checked_add(h.len()))
            .ok_or_else(|| bad("outline count overflow"))?;
        cap(count.saturating_mul(2), budget.max_dimensions, "dimensions")?;
        cap(
            count
                .checked_mul(count)
                .ok_or_else(|| bad("outline work overflow"))?,
            budget.max_work,
            "loop intersection work",
        )?;
        let mut loops = vec![outer];
        loops.extend(holes);
        let mut vertices = Vec::with_capacity(loops.len());
        for points in &loops {
            if points.len() < 3 {
                return Err(bad("loop needs three distinct vertices"));
            }
            let mut p = Vec::with_capacity(points.len());
            for point in points {
                for d in point {
                    d.validate(&label, LENGTH)?;
                    if d.value.value.abs() > 1e12 {
                        return Err(bad(
                            "outline coordinate exceeds exact-predicate arithmetic envelope",
                        ));
                    }
                }
                p.push([point[0].value.value, point[1].value.value]);
            }
            vertices.push(p);
        }
        for (index, p) in vertices.iter().enumerate() {
            for i in 0..p.len() {
                if p[i] == p[(i + 1) % p.len()] {
                    return Err(bad("repeated adjacent contour vertex"));
                }
                let prev = p[(i + p.len() - 1) % p.len()];
                let next = p[(i + 1) % p.len()];
                if orient2d(prev, p[i], next) == Sign::Zero
                    && (on(prev, p[i], next) || on(p[i], next, prev))
                {
                    return Err(GeneratorError::SelfIntersection { loop_index: index });
                }
                for j in i + 1..p.len() {
                    if j == i + 1 || (i == 0 && j == p.len() - 1) {
                        continue;
                    }
                    if crosses(p[i], p[(i + 1) % p.len()], p[j], p[(j + 1) % p.len()]) {
                        return Err(GeneratorError::SelfIntersection { loop_index: index });
                    }
                }
            }
            if signed_area(p) == 0.0 {
                return Err(bad("zero-area contour"));
            }
        }
        for h in 1..vertices.len() {
            if !vertices[h].iter().all(|p| inside(*p, &vertices[0])) {
                return Err(GeneratorError::HoleOutsideOutline { hole: h - 1 });
            }
            for k in 0..h {
                let a = &vertices[h];
                let b = &vertices[k];
                let intersects = (0..a.len()).any(|i| {
                    (0..b.len())
                        .any(|j| crosses(a[i], a[(i + 1) % a.len()], b[j], b[(j + 1) % b.len()]))
                });
                if intersects {
                    return Err(if k == 0 {
                        GeneratorError::HoleOutsideOutline { hole: h - 1 }
                    } else {
                        GeneratorError::OverlappingHoles
                    });
                }
                if k > 0 && (inside(a[0], b) || inside(b[0], a)) {
                    return Err(GeneratorError::OverlappingHoles);
                }
            }
        }
        let p = &vertices[0];
        let rectangle = p.len() == 4
            && (0..4).all(|i| {
                let a = p[i];
                let b = p[(i + 1) % 4];
                let c = p[(i + 2) % 4];
                (b[0] - a[0]) * (c[0] - b[0]) + (b[1] - a[1]) * (c[1] - b[1]) == 0.0
            });
        if rectangle && !fixture(&label) {
            return Err(GeneratorError::FixtureRequired);
        }
        let mut trims = Vec::with_capacity(loops.len());
        for (index, p) in vertices.iter().enumerate() {
            let mut controls = p
                .iter()
                .map(|p| Ok([rational(p[0])?, rational(p[1])?]))
                .collect::<Result<Vec<_>, GeneratorError>>()?;
            if (signed_area(p) > 0.0) != (index == 0) {
                controls.reverse();
            }
            controls.push(controls[0]);
            let mut knots = vec![Rat::int(0)];
            knots.extend((0..controls.len()).map(|i| Rat::new(i as i128, 1)));
            knots.push(Rat::new((controls.len() - 1) as i128, 1));
            let curve = NurbsCurve::new(
                KnotVector::new(knots, 1).map_err(|e| bad(e.to_string()))?,
                &controls,
                &vec![Rat::int(1); controls.len()],
            )
            .map_err(|e| bad(e.to_string()))?;
            trims.push(TrimLoop::new(curve).map_err(|e| bad(e.to_string()))?);
        }
        Ok(Self {
            label,
            loops,
            vertices,
            patch: TrimmedPatch::new(trims),
        })
    }
    /// Raw observations remain in their original order and binary64 encoding.
    pub fn contours(&self) -> &[Vec<[MeasuredDimension; 2]>] {
        &self.loops
    }
    /// Direct Polyline input to fs-mesh::planar, without a second mesher.
    pub fn vertices_m(&self) -> &[Vec<[f64; 2]>] {
        &self.vertices
    }
    pub fn trimmed_patch(&self) -> &TrimmedPatch {
        &self.patch
    }
    pub fn area(&self) -> Result<DerivedDimension, GeneratorError> {
        let mut value = 0.0;
        let mut enclosure = Interval::point(0.0);
        let mut inputs = Vec::new();
        for (index, points) in self.loops.iter().enumerate() {
            let sign = if signed_area(&self.vertices[index]) > 0.0 {
                1.0
            } else {
                -1.0
            };
            let multiplier = if index == 0 { sign } else { -sign };
            value += if index == 0 {
                signed_area(&self.vertices[index]).abs()
            } else {
                -signed_area(&self.vertices[index]).abs()
            };
            let mut a = Interval::point(0.0);
            let o = &points[0];
            for i in 1..points.len() - 1 {
                let p = &points[i];
                let q = &points[i + 1];
                a = a + (p[0].interval() - o[0].interval()) * (q[1].interval() - o[1].interval())
                    - (p[1].interval() - o[1].interval()) * (q[0].interval() - o[0].interval());
            }
            enclosure = enclosure + a * Interval::point(multiplier / 2.0);
            inputs.extend(points.iter().flat_map(|p| p.iter().cloned()));
        }
        DerivedDimension::new(value, enclosure, AREA, inputs)
    }
    pub fn geometry_source(&self) -> GeometrySource {
        let mut w = TableWriter::new(&self.label);
        w.count("loops", self.loops.len());
        for (i, p) in self.loops.iter().enumerate() {
            w.count(&format!("loop{i}"), p.len());
            for point in p {
                for d in point {
                    w.dimension(d);
                }
            }
        }
        w.source("OutlineWithHoles.v1")
    }
    pub(super) fn read(
        label: String,
        r: &mut TableReader<'_>,
        b: GeneratorBudget,
    ) -> Result<Self, GeneratorError> {
        let n = r.count("loops")?;
        cap(n, b.max_dimensions / 6, "loops")?;
        if n == 0 {
            return Err(bad("missing outer contour"));
        }
        let mut loops = Vec::new();
        let mut total = 0usize;
        for i in 0..n {
            let count = r.count(&format!("loop{i}"))?;
            total = total
                .checked_add(count.saturating_mul(2))
                .ok_or_else(|| bad("dimension count overflow"))?;
            cap(total, b.max_dimensions, "dimensions")?;
            loops.push(
                (0..count)
                    .map(|_| Ok([r.dimension()?, r.dimension()?]))
                    .collect::<Result<Vec<_>, GeneratorError>>()?,
            );
        }
        let outer = loops.remove(0);
        Self::new(label, outer, loops, b)
    }
}
