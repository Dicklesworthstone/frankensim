use super::*;
use crate::acoustic::CylinderSegment;
use fs_rep_nurbs::{KnotVector, NurbsCurve};

/// Linear radius or cubic Hermite radius with explicitly measured dr/dx.
/// Cubics are represented and evaluated by fs-rep-nurbs, not another spline.
#[derive(Debug, Clone, PartialEq)]
pub enum BoreInterpolant {
    PiecewiseConical,
    CubicHermite { slopes: Vec<MeasuredDimension> },
}
impl BoreInterpolant {
    /// Interpolation order for smooth data and exact supplied derivatives.
    pub fn formal_order(&self) -> u32 {
        match self {
            Self::PiecewiseConical => 2,
            Self::CubicHermite { .. } => 4,
        }
    }
}
#[derive(Debug)]
pub struct BoreProfile {
    pub(super) label: String,
    stations: Vec<(MeasuredDimension, MeasuredDimension)>,
    interpolant: BoreInterpolant,
    curve: NurbsCurve<f64, 1>,
}
#[derive(Debug, Clone, PartialEq)]
pub struct GeneratedSegment {
    pub value: CylinderSegment,
    pub length: DerivedDimension,
    pub inlet_radius: DerivedDimension,
    pub outlet_radius: DerivedDimension,
}
#[derive(Debug, Clone, PartialEq)]
pub struct BoreSampling {
    /// Every original station remains, in order, with inserted samples between.
    pub stations: Vec<[DerivedDimension; 2]>,
    pub segments: Vec<GeneratedSegment>,
}
impl BoreProfile {
    pub fn new(
        label: String,
        stations: Vec<(MeasuredDimension, MeasuredDimension)>,
        interpolant: BoreInterpolant,
        budget: GeneratorBudget,
    ) -> Result<Self, GeneratorError> {
        label_valid(&label)?;
        let n = stations.len();
        let count = n
            .checked_mul(
                if matches!(&interpolant, BoreInterpolant::CubicHermite { .. }) {
                    3
                } else {
                    2
                },
            )
            .ok_or_else(|| bad("station count overflow"))?;
        cap(count, budget.max_dimensions, "dimensions")?;
        cap(n.saturating_mul(4), budget.max_work, "spline controls")?;
        let axis: Vec<_> = stations.iter().map(|s| s.0.clone()).collect();
        checked_axis(&axis, &label)?;
        for (_, r) in &stations {
            r.validate(&label, LENGTH)?;
            if r.bound.lower.value <= 0.0 {
                return Err(bad("positive bore radius bound required"));
            }
        }
        let x: Vec<_> = axis.iter().map(|d| d.value.value).collect();
        let r: Vec<_> = stations.iter().map(|s| s.1.value.value).collect();
        let (knots, controls) = match &interpolant {
            BoreInterpolant::PiecewiseConical => (
                linear_knots(&x)?,
                r.iter().map(|v| [*v]).collect::<Vec<_>>(),
            ),
            BoreInterpolant::CubicHermite { slopes } => {
                if slopes.len() != n {
                    return Err(bad("one measured slope per station required"));
                }
                for s in slopes {
                    s.validate(&label, Dims::NONE)?;
                }
                let mut knots = vec![x[0]; 4];
                for &v in &x[1..n - 1] {
                    knots.extend([v; 3]);
                }
                knots.extend([x[n - 1]; 4]);
                let mut controls = vec![[r[0]]];
                for i in 0..n - 1 {
                    let h = x[i + 1] - x[i];
                    let a = r[i] + h * slopes[i].value.value / 3.0;
                    let b = r[i + 1] - h * slopes[i + 1].value.value / 3.0;
                    // Convex-hull positivity also applies to the uncertain net.
                    let hi = axis[i + 1].interval() - axis[i].interval();
                    let ai =
                        stations[i].1.interval() + hi * slopes[i].interval() / Interval::point(3.0);
                    let bi = stations[i + 1].1.interval()
                        - hi * slopes[i + 1].interval() / Interval::point(3.0);
                    if ai.lo() <= 0.0 || bi.lo() <= 0.0 {
                        return Err(bad("cubic radius positivity not enclosed"));
                    }
                    controls.extend([[a], [b], [r[i + 1]]]);
                }
                (
                    KnotVector::new(knots, 3).map_err(|e| bad(e.to_string()))?,
                    controls,
                )
            }
        };
        let linear = (1..n - 1)
            .all(|i| (r[i] - r[0]) * (x[n - 1] - x[0]) == (r[n - 1] - r[0]) * (x[i] - x[0]));
        let ideal = match &interpolant {
            BoreInterpolant::PiecewiseConical => linear,
            BoreInterpolant::CubicHermite { slopes } => {
                linear
                    && slopes
                        .iter()
                        .all(|s| s.value.value == (r[n - 1] - r[0]) / (x[n - 1] - x[0]))
            }
        };
        if ideal && !fixture(&label) {
            return Err(GeneratorError::FixtureRequired);
        }
        let curve = NurbsCurve::new(knots, &controls, &vec![1.0; controls.len()])
            .map_err(|e| bad(e.to_string()))?;
        Ok(Self {
            label,
            stations,
            interpolant,
            curve,
        })
    }
    pub fn stations(&self) -> &[(MeasuredDimension, MeasuredDimension)] {
        &self.stations
    }
    pub fn interpolant(&self) -> &BoreInterpolant {
        &self.interpolant
    }
    pub fn curve(&self) -> &NurbsCurve<f64, 1> {
        &self.curve
    }
    fn inputs(&self) -> Vec<MeasuredDimension> {
        let mut v: Vec<_> = self
            .stations
            .iter()
            .flat_map(|s| [s.0.clone(), s.1.clone()])
            .collect();
        if let BoreInterpolant::CubicHermite { slopes } = &self.interpolant {
            v.extend(slopes.iter().cloned());
        }
        v
    }
    fn at(&self, i: usize, t: f64) -> Result<[DerivedDimension; 2], GeneratorError> {
        let (x0, r0) = &self.stations[i];
        let (x1, r1) = &self.stations[i + 1];
        let x = if t == 0.0 {
            x0.value.value
        } else if t == 1.0 {
            x1.value.value
        } else {
            (1.0 - t) * x0.value.value + t * x1.value.value
        };
        let radius = if t == 0.0 {
            r0.value.value
        } else if t == 1.0 {
            r1.value.value
        } else {
            self.curve.eval(x).map_err(|e| bad(e.to_string()))?[0]
        };
        let ti = Interval::point(t);
        let ui = Interval::point(1.0 - t);
        let xb = ui * x0.interval() + ti * x1.interval();
        let rb = match &self.interpolant {
            BoreInterpolant::PiecewiseConical => ui * r0.interval() + ti * r1.interval(),
            BoreInterpolant::CubicHermite { slopes } => {
                let h = x1.interval() - x0.interval();
                let a = r0.interval() + h * slopes[i].interval() / Interval::point(3.0);
                let b = r1.interval() - h * slopes[i + 1].interval() / Interval::point(3.0);
                // Bernstein weights are nonnegative throughout each span.
                ui * ui * ui * r0.interval()
                    + Interval::point(3.0) * ui * ui * ti * a
                    + Interval::point(3.0) * ui * ti * ti * b
                    + ti * ti * ti * r1.interval()
            }
        };
        Ok([
            DerivedDimension::new(x, xb, LENGTH, self.inputs())?,
            DerivedDimension::new(radius, rb, LENGTH, self.inputs())?,
        ])
    }
    /// Query nominal station coordinates; the interval follows uncertain
    /// material stations at the same normalized span coordinate.
    pub fn radius_at(&self, x: f64) -> Result<DerivedDimension, GeneratorError> {
        if let Some(i) = self
            .stations
            .iter()
            .position(|s| s.0.value.value.to_bits() == x.to_bits())
        {
            let r = &self.stations[i].1;
            return DerivedDimension::new(r.value.value, r.interval(), LENGTH, self.inputs());
        }
        let i = self
            .stations
            .windows(2)
            .position(|s| x > s[0].0.value.value && x < s[1].0.value.value)
            .ok_or_else(|| bad("bore query outside station domain"))?;
        Ok(self.at(
            i,
            (x - self.stations[i].0.value.value)
                / (self.stations[i + 1].0.value.value - self.stations[i].0.value.value),
        )?[1]
            .clone())
    }
    pub fn sample(
        &self,
        rung: SamplingRung,
        budget: GeneratorBudget,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<BoreSampling, GeneratorError> {
        let m = rung.subdivisions_per_span;
        if m == 0 {
            return Err(bad("positive sampling rung required"));
        }
        let n = (self.stations.len() - 1)
            .checked_mul(m)
            .and_then(|n| n.checked_add(1))
            .ok_or_else(|| bad("sampling count overflow"))?;
        cap(n, budget.max_samples, "samples")?;
        let input_count = self.stations.len()
            * if matches!(self.interpolant, BoreInterpolant::CubicHermite { .. }) {
                3
            } else {
                2
            };
        cap(
            n.checked_mul(input_count)
                .and_then(|work| work.checked_mul(5))
                .ok_or_else(|| bad("sampling work overflow"))?,
            budget.max_work,
            "sample provenance work",
        )?;
        if cancelled() {
            return Err(GeneratorError::Cancelled);
        }
        let mut stations = Vec::with_capacity(n);
        for i in 0..self.stations.len() - 1 {
            for j in 0..m {
                if cancelled() {
                    return Err(GeneratorError::Cancelled);
                }
                stations.push(self.at(i, j as f64 / m as f64)?);
            }
        }
        stations.push(self.at(self.stations.len() - 2, 1.0)?);
        let mut segments = Vec::with_capacity(n - 1);
        for pair in stations.windows(2) {
            if cancelled() {
                return Err(GeneratorError::Cancelled);
            }
            let a = &pair[0];
            let b = &pair[1];
            let length = b[0].value.value - a[0].value.value;
            let enclosure = Interval::new(b[0].lower.value, b[0].upper.value)
                - Interval::new(a[0].lower.value, a[0].upper.value);
            segments.push(GeneratedSegment {
                value: CylinderSegment::taper(a[1].value.value, b[1].value.value, length),
                length: DerivedDimension::new(length, enclosure, LENGTH, self.inputs())?,
                inlet_radius: a[1].clone(),
                outlet_radius: b[1].clone(),
            });
        }
        Ok(BoreSampling { stations, segments })
    }
    /// Exact frustum sum for piecewise-conical input. Cubics must select an
    /// explicit conical sampling rung instead; no hidden geometry reduction.
    pub fn volume(&self) -> Result<DerivedDimension, GeneratorError> {
        if !matches!(self.interpolant, BoreInterpolant::PiecewiseConical) {
            return Err(bad("cubic volume requires an explicit sampling rung"));
        }
        let mut value = 0.0;
        let mut bound = Interval::point(0.0);
        for w in self.stations.windows(2) {
            let h = w[1].0.value.value - w[0].0.value.value;
            let a = w[0].1.value.value;
            let b = w[1].1.value.value;
            value += std::f64::consts::PI * h * (a * a + a * b + b * b) / 3.0;
            let hi = w[1].0.interval() - w[0].0.interval();
            let ai = w[0].1.interval();
            let bi = w[1].1.interval();
            bound = bound
                + Interval::new(
                    std::f64::consts::PI.next_down(),
                    std::f64::consts::PI.next_up(),
                ) * hi
                    * (ai * ai + ai * bi + bi * bi)
                    / Interval::point(3.0);
        }
        DerivedDimension::new(value, bound, VOLUME, self.inputs())
    }
    pub fn geometry_source(&self) -> GeometrySource {
        let mut w = TableWriter::new(&self.label);
        w.count("stations", self.stations.len());
        w.count(
            "cubic",
            usize::from(matches!(
                self.interpolant,
                BoreInterpolant::CubicHermite { .. }
            )),
        );
        for (x, r) in &self.stations {
            w.dimension(x);
            w.dimension(r);
        }
        if let BoreInterpolant::CubicHermite { slopes } = &self.interpolant {
            for s in slopes {
                w.dimension(s);
            }
        }
        w.source("BoreProfile.v1")
    }
    pub(super) fn read(
        label: String,
        r: &mut TableReader<'_>,
        b: GeneratorBudget,
    ) -> Result<Self, GeneratorError> {
        let n = r.count("stations")?;
        let cubic = r.count("cubic")?;
        if cubic > 1 {
            return Err(bad("invalid bore interpolant"));
        }
        cap(n.saturating_mul(2 + cubic), b.max_dimensions, "dimensions")?;
        let stations = (0..n)
            .map(|_| Ok((r.dimension()?, r.dimension()?)))
            .collect::<Result<Vec<_>, GeneratorError>>()?;
        let interpolant = if cubic == 0 {
            BoreInterpolant::PiecewiseConical
        } else {
            BoreInterpolant::CubicHermite {
                slopes: (0..n)
                    .map(|_| r.dimension())
                    .collect::<Result<Vec<_>, _>>()?,
            }
        };
        Self::new(label, stations, interpolant, b)
    }
}
