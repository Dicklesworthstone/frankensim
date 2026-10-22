use super::*;
use fs_rep_nurbs::NurbsSurface;
use std::sync::Arc;

/// Measured rectangular contour table, represented by a tensor degree-one
/// NURBS graph. This explicitly declares second-order interpolation, with C0
/// cell seams. A denser measured table refines the geometry; no sample is
/// decimated. Thickness is independent geometry data, never inferred.
#[derive(Debug)]
pub struct ArchingSurface {
    pub(super) label: String,
    x: Vec<MeasuredDimension>,
    y: Vec<MeasuredDimension>,
    heights: Vec<MeasuredDimension>,
    thickness: Option<Vec<MeasuredDimension>>,
    reference_state: ReferenceState,
    surface: Arc<NurbsSurface<f64>>,
    thickness_surface: Option<Arc<NurbsSurface<f64>>>,
}
/// The degree-one tensor graph through every table value, with physical x/y
/// as its parameters. Heights and thickness share this one interpolant.
fn graph(
    xs: &[f64],
    ys: &[f64],
    values: &[MeasuredDimension],
) -> Result<Arc<NurbsSurface<f64>>, GeneratorError> {
    let points: Vec<Vec<_>> = (0..xs.len())
        .map(|i| {
            (0..ys.len())
                .map(|j| [xs[i], ys[j], values[j * xs.len() + i].value.value])
                .collect()
        })
        .collect();
    let weights = vec![vec![1.0; ys.len()]; xs.len()];
    Ok(Arc::new(
        NurbsSurface::new(linear_knots(xs)?, linear_knots(ys)?, &points, &weights)
            .map_err(|e| bad(e.to_string()))?,
    ))
}
impl ArchingSurface {
    pub const FORMAL_ORDER: u32 = 2;
    pub fn new(
        label: String,
        x: Vec<MeasuredDimension>,
        y: Vec<MeasuredDimension>,
        heights: Vec<MeasuredDimension>,
        thickness: Option<Vec<MeasuredDimension>>,
        reference_state: ReferenceState,
        budget: GeneratorBudget,
    ) -> Result<Self, GeneratorError> {
        label_valid(&label)?;
        let n = x
            .len()
            .checked_mul(y.len())
            .ok_or_else(|| bad("contour count overflow"))?;
        let count = n
            .checked_mul(1 + usize::from(thickness.is_some()))
            .and_then(|n| n.checked_add(x.len()))
            .and_then(|n| n.checked_add(y.len()))
            .ok_or_else(|| bad("dimension count overflow"))?;
        cap(count, budget.max_dimensions, "dimensions")?;
        cap(n, budget.max_work, "surface controls")?;
        checked_axis(&x, &label)?;
        checked_axis(&y, &label)?;
        if heights.len() != n || thickness.as_ref().is_some_and(|v| v.len() != n) {
            return Err(bad(
                "every grid height and optional thickness must be supplied",
            ));
        }
        for h in &heights {
            h.validate(&label, LENGTH)?;
        }
        let z0 = heights[0].value.value;
        let dx = x[1].value.value - x[0].value.value;
        let dy = y[1].value.value - y[0].value.value;
        let zx = (heights[1].value.value - z0) / dx;
        let zy = (heights[x.len()].value.value - z0) / dy;
        let plane = heights.iter().enumerate().all(|(k, h)| {
            let expected = z0
                + zx * (x[k % x.len()].value.value - x[0].value.value)
                + zy * (y[k / x.len()].value.value - y[0].value.value);
            (h.value.value - expected).abs()
                <= 64.0 * f64::EPSILON * h.value.value.abs().max(expected.abs())
        });
        if plane && !fixture(&label) {
            return Err(GeneratorError::FixtureRequired);
        }
        if let Some(values) = &thickness {
            for h in values {
                h.validate(&label, LENGTH)?;
                if h.bound.lower.value <= 0.0 {
                    return Err(bad("positive thickness bound required"));
                }
            }
        }
        let xs: Vec<_> = x.iter().map(|d| d.value.value).collect();
        let ys: Vec<_> = y.iter().map(|d| d.value.value).collect();
        let surface = graph(&xs, &ys, &heights)?;
        let thickness_surface = thickness
            .as_deref()
            .map(|values| graph(&xs, &ys, values))
            .transpose()?;
        Ok(Self {
            label,
            x,
            y,
            heights,
            thickness,
            reference_state,
            surface,
            thickness_surface,
        })
    }
    /// The source label: a ledger-backed description or a `test-fixture:`.
    pub fn label(&self) -> &str {
        &self.label
    }
    /// Shared immutable height graph, lowered by fs-plate's
    /// `ArchingSurface::from_measured`.
    pub fn surface(&self) -> &Arc<NurbsSurface<f64>> {
        &self.surface
    }
    /// The measured graduation on the same degree-one graph, lowered by
    /// fs-plate's `ThicknessField::from_measured`; `None` when unmeasured.
    pub fn thickness_surface(&self) -> Option<&Arc<NurbsSurface<f64>>> {
        self.thickness_surface.as_ref()
    }
    pub fn axes(&self) -> (&[MeasuredDimension], &[MeasuredDimension]) {
        (&self.x, &self.y)
    }
    /// Row-major y/x, retaining all source observations bit-exactly.
    pub fn contour_heights(&self) -> &[MeasuredDimension] {
        &self.heights
    }
    pub fn thickness(&self) -> Option<&[MeasuredDimension]> {
        self.thickness.as_deref()
    }
    pub fn reference_state(&self) -> ReferenceState {
        self.reference_state
    }
    fn inputs(&self, thickness: bool) -> Vec<MeasuredDimension> {
        self.x
            .iter()
            .chain(&self.y)
            .chain(if thickness {
                self.thickness.as_deref().unwrap_or(&[])
            } else {
                &self.heights
            })
            .cloned()
            .collect()
    }
    fn cell(axis: &[MeasuredDimension], v: f64) -> Result<(usize, f64), GeneratorError> {
        if !v.is_finite() || v < axis[0].value.value || v > axis[axis.len() - 1].value.value {
            return Err(bad("query outside measured contour domain"));
        }
        let i = axis
            .windows(2)
            .position(|w| v <= w[1].value.value)
            .ok_or_else(|| bad("missing contour cell"))?;
        Ok((
            i,
            (v - axis[i].value.value) / (axis[i + 1].value.value - axis[i].value.value),
        ))
    }
    fn field_at(
        &self,
        p: [f64; 2],
        values: &[MeasuredDimension],
        thickness: bool,
    ) -> Result<DerivedDimension, GeneratorError> {
        let (i, t) = Self::cell(&self.x, p[0])?;
        let (j, u) = Self::cell(&self.y, p[1])?;
        if let (Some(ix), Some(iy)) = (
            self.x
                .iter()
                .position(|v| v.value.value.to_bits() == p[0].to_bits()),
            self.y
                .iter()
                .position(|v| v.value.value.to_bits() == p[1].to_bits()),
        ) {
            let d = &values[iy * self.x.len() + ix];
            return DerivedDimension::new(
                d.value.value,
                d.interval(),
                LENGTH,
                self.inputs(thickness),
            );
        }
        let mut enclosure = Interval::point(0.0);
        for (dx, wx) in [(0, 1.0 - t), (1, t)] {
            for (dy, wy) in [(0, 1.0 - u), (1, u)] {
                let d = &values[(j + dy) * self.x.len() + i + dx];
                enclosure = enclosure + Interval::point(wx) * Interval::point(wy) * d.interval();
            }
        }
        // The value is the same graph fs-plate lowers, never a second fit.
        let chart = if thickness {
            self.thickness_surface
                .as_ref()
                .ok_or_else(|| bad("no measured thickness field"))?
        } else {
            &self.surface
        };
        let value = chart.eval(p[0], p[1]).map_err(|e| bad(e.to_string()))?[2];
        DerivedDimension::new(value, enclosure, LENGTH, self.inputs(thickness))
    }
    /// Bounds follow the same normalized cell coordinate as uncertain axes;
    /// they do not certify registration at a fixed unknown physical position.
    pub fn height_at(&self, p: [f64; 2]) -> Result<DerivedDimension, GeneratorError> {
        self.field_at(p, &self.heights, false)
    }
    pub fn thickness_at(&self, p: [f64; 2]) -> Result<DerivedDimension, GeneratorError> {
        self.field_at(
            p,
            self.thickness
                .as_deref()
                .ok_or_else(|| bad("no measured thickness field"))?,
            true,
        )
    }
    pub fn geometry_source(&self) -> GeometrySource {
        let mut w = TableWriter::new(&self.label);
        w.count("x", self.x.len());
        w.count("y", self.y.len());
        w.count("thickness", usize::from(self.thickness.is_some()));
        w.count(
            "reference",
            match self.reference_state {
                ReferenceState::Unloaded => 0,
                ReferenceState::Strung => 1,
                ReferenceState::Assembled => 2,
                ReferenceState::Played => 3,
            },
        );
        for d in self.x.iter().chain(&self.y).chain(&self.heights) {
            w.dimension(d);
        }
        if let Some(v) = &self.thickness {
            for d in v {
                w.dimension(d);
            }
        }
        w.source("ArchingSurface.v1")
    }
    pub(super) fn read(
        label: String,
        r: &mut TableReader<'_>,
        b: GeneratorBudget,
    ) -> Result<Self, GeneratorError> {
        let nx = r.count("x")?;
        let ny = r.count("y")?;
        let th = r.count("thickness")?;
        if th > 1 {
            return Err(bad("invalid thickness presence"));
        }
        let state = match r.count("reference")? {
            0 => ReferenceState::Unloaded,
            1 => ReferenceState::Strung,
            2 => ReferenceState::Assembled,
            3 => ReferenceState::Played,
            _ => return Err(bad("invalid reference state")),
        };
        let n = nx
            .checked_mul(ny)
            .ok_or_else(|| bad("grid count overflow"))?;
        cap(
            nx.saturating_add(ny)
                .saturating_add(n.saturating_mul(1 + th)),
            b.max_dimensions,
            "dimensions",
        )?;
        let x = (0..nx)
            .map(|_| r.dimension())
            .collect::<Result<Vec<_>, _>>()?;
        let y = (0..ny)
            .map(|_| r.dimension())
            .collect::<Result<Vec<_>, _>>()?;
        let heights = (0..n)
            .map(|_| r.dimension())
            .collect::<Result<Vec<_>, _>>()?;
        let thickness = if th == 0 {
            None
        } else {
            Some(
                (0..n)
                    .map(|_| r.dimension())
                    .collect::<Result<Vec<_>, _>>()?,
            )
        };
        Self::new(label, x, y, heights, thickness, state, b)
    }
}
