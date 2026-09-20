//! Numerical oriented surface rules for the same implicit field as `quad3`.
//!
//! Monotone height graphs give dS = |grad phi| / |d_height phi| dA and the
//! outward normal grad(phi)/|grad(phi)| for the negative-inside convention.
//! Box and derivative admission and root contraction reuse the bulk machinery.
//! Surface area, normals and crossing selection are NUMERICAL, not enclosures.
//! In particular, endpoint scalar samples resolve crossings only after their
//! interval uncertainty is smaller than the requested geometric resolution.
//! A sampled zero at a face is owned by its material-side box, never both.
use super::*;

/// Additional admission policy for surface normals. Geometric resolution and
/// all point/box/field work limits come from the shared QuadratureControl3.
#[derive(Debug, Clone, Copy)]
pub struct SurfaceOptions3 {
    /// Maximum point-gradient interval width divided by the admitted lower
    /// bound on the height derivative. A broad derivative refuses a normal;
    /// its midpoint is never silently accepted as accurate.
    pub relative_normal_tolerance: f64,
}
impl Default for SurfaceOptions3 {
    fn default() -> Self { Self { relative_normal_tolerance: 1e-6 } }
}
impl SurfaceOptions3 {
    pub(crate) fn validate(self) -> Result<(), QuadratureError3> {
        if !self.relative_normal_tolerance.is_finite()
            || self.relative_normal_tolerance <= 0.0 || self.relative_normal_tolerance >= 1.0 {
            return Err(QuadratureError3::Invalid("surface normal tolerance must lie in (0,1)"));
        }
        Ok(())
    }
}

/// A numerical interface point with positive area measure. Normal orientation
/// follows phi<0, not the coordinate direction chosen for height integration.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SurfacePoint3 {
    pub position: [f64; 3],
    pub normal: [f64; 3],
    pub weight: f64,
}
/// Completed interface rule. Empty means no sampled crossing, NOT a certificate
/// that an unresolved surface patch is absent. No rule escapes on a returned error.
#[derive(Debug, Clone)]
pub struct SurfaceRules3 {
    points: Vec<SurfacePoint3>,
}
impl SurfaceRules3 {
    #[must_use] pub fn points(&self) -> &[SurfacePoint3] { &self.points }
    /// Numerical area, not a certified enclosure or a box-boundary area.
    #[must_use] pub fn area(&self) -> f64 { self.points.iter().map(|p| p.weight).sum() }
}

/// Integrate only the zero-level interface inside a box. Artificial background
/// faces on which phi is strictly negative are NOT part of this rule. The
/// scalar `value` and interval/derivative producers must describe the same pure
/// field. Endpoint values outside their own enclosures are refused.
///
/// Material-side ownership counts a zero on an upper height face only when
/// phi increases, and on a lower height face only when phi decreases. This
/// avoids double-counting a plane coincident with a subdivision/grid face.
/// Endpoint decisions within the root tolerance remain numerical decisions,
/// not certified root-existence claims. Bulk volume enclosures are unchanged.
///
/// `work.points` and its allowance count BOTH bulk and surface emissions when
/// the control is shared. Failed candidates do not refund consumed work.
pub fn surface_cell_rules3(sdf: &dyn CutSdf3, cell: HexCell, options: SurfaceOptions3,
    control: &mut QuadratureControl3<'_>) -> Result<SurfaceRules3, QuadratureError3> {
    options.validate()?;
    control.poll()?;
    let mut result = SurfaceRules3 { points: Vec::new() };
    visit_surface(sdf, cell, control.options.depth, options, control, &mut result)?;
    if !result.area().is_finite() { return Err(QuadratureError3::Invalid("surface area overflow")); }
    control.poll()?;
    Ok(result)
}

// Interpret an endpoint only when a zero-containing enclosure is sufficiently
// narrow in geometric units. Use the same global endpoint coordinate in both
// incident boxes, so deterministic samples give complementary ownership.
fn endpoint(sdf: &dyn CutSdf3, p: [f64; 3], resolution: f64, slope: f64,
    control: &mut QuadratureControl3<'_>) -> Result<f64, QuadratureError3> {
    let enclosure = control.field(|| sdf.enclose(p, p))?;
    if enclosure.lo() <= 0.0 && enclosure.hi() >= 0.0 {
        let uncertainty = enclosure.hi() / slope - enclosure.lo() / slope;
        if !uncertainty.is_finite() || uncertainty > resolution {
            return Err(QuadratureError3::Invalid("surface endpoint sign unresolved at requested resolution"));
        }
    }
    control.poll()?;
    control.work.field_evaluations = control.work.field_evaluations.checked_add(1)
        .ok_or(QuadratureError3::Invalid("field counter overflow"))?;
    let value = sdf.value(p);
    control.poll()?;
    if !value.is_finite() || value < enclosure.lo() || value > enclosure.hi() {
        return Err(QuadratureError3::Invalid("surface scalar value disagrees with its enclosure"));
    }
    Ok(value)
}

fn visit_surface(sdf: &dyn CutSdf3, cell: HexCell, depth: u32, options: SurfaceOptions3,
    control: &mut QuadratureControl3<'_>, out: &mut SurfaceRules3) -> Result<(), QuadratureError3> {
    control.poll()?;
    if control.work.boxes >= control.options.max_boxes { return Err(QuadratureError3::BoxBudget); }
    control.work.boxes += 1;
    box_volume(cell)?;
    let (lo, hi) = (cell.lo(), cell.hi());
    let sign = control.field(|| sdf.enclose(lo, hi))?;
    if sign.lo() >= 0.0 || sign.hi() < 0.0 { return Ok(()); }
    if depth > 0 {
        let mid: [f64; 3] = std::array::from_fn(|a| f64::midpoint(lo[a], hi[a]));
        if (0..3).any(|a| mid[a] <= lo[a] || mid[a] >= hi[a]) {
            return Err(QuadratureError3::Invalid("surface subdivision cannot advance"));
        }
        for octant in 0..8 {
            let a = std::array::from_fn(|i| if octant & (1 << i) == 0 { lo[i] } else { mid[i] });
            let b = std::array::from_fn(|i| if octant & (1 << i) == 0 { mid[i] } else { hi[i] });
            let child = HexCell::try_new(a, b).map_err(|_| QuadratureError3::Invalid("invalid surface subdivision"))?;
            visit_surface(sdf, child, depth-1, options, control, out)?;
        }
        return Ok(());
    }
    let (axis, increasing, slope) = select_height(sdf, cell, control)?;
    let bases: Vec<_> = (0..3).filter(|&i| i != axis).collect();
    let (a, b) = (bases[0], bases[1]);
    let sa = 0.5*(hi[a]-lo[a]); let sb = 0.5*(hi[b]-lo[b]);
    let resolution = (hi[axis]-lo[axis])*control.options.root_relative_tolerance;
    if !resolution.is_finite() || resolution <= 0.0 {
        return Err(QuadratureError3::Invalid("surface root resolution under/overflow"));
    }
    for (ga, wa) in GAUSS { for (gb, wb) in GAUSS {
        let mut p = lo;
        p[a] = f64::midpoint(lo[a], hi[a]) + sa*ga;
        p[b] = f64::midpoint(lo[b], hi[b]) + sb*gb;
        p[axis] = lo[axis]; let left = endpoint(sdf, p, resolution, slope, control)?;
        p[axis] = hi[axis]; let right = endpoint(sdf, p, resolution, slope, control)?;
        if (increasing && left > right) || (!increasing && left < right) {
            return Err(QuadratureError3::Invalid("surface samples contradict monotone derivative"));
        }
        let crossing = if increasing { left < 0.0 && right >= 0.0 }
            else { left >= 0.0 && right < 0.0 };
        if !crossing { continue; }
        p[axis] = if increasing && right == 0.0 { hi[axis] }
            else if !increasing && left == 0.0 { lo[axis] }
            else { height_cut(sdf, cell, p, axis, increasing, slope, control)? };
        let mut gradient = [0.0; 3];
        for i in 0..3 {
            let d = control.field(|| sdf.derivative_enclose(p, p, AXES[i]))?;
            let width = d.hi()/slope - d.lo()/slope;
            if !width.is_finite() || width > options.relative_normal_tolerance {
                return Err(QuadratureError3::Invalid("surface point-gradient enclosure too broad"));
            }
            gradient[i] = f64::midpoint(d.lo(), d.hi());
        }
        if (increasing && gradient[axis] <= 0.0) || (!increasing && gradient[axis] >= 0.0) {
            return Err(QuadratureError3::Invalid("surface normal contradicts height direction"));
        }
        // Normalize relative to the height derivative first; this avoids
        // squaring potentially enormous/small dimensional gradient values.
        let oriented: [f64; 3] = gradient.map(|v| v/gradient[axis].abs());
        let jacobian = oriented[0].hypot(oriented[1]).hypot(oriented[2]);
        let normal = oriented.map(|v| v/jacobian);
        let weight = (wa*wb*sa*sb)*jacobian;
        if !p.iter().all(|v| v.is_finite()) || !normal.iter().all(|v| v.is_finite())
            || !weight.is_finite() || weight <= 0.0 {
            return Err(QuadratureError3::Invalid("unrepresentable oriented surface rule"));
        }
        if control.work.points >= control.options.max_points { return Err(QuadratureError3::PointBudget); }
        control.work.points += 1;
        out.points.push(SurfacePoint3 { position: p, normal, weight });
    } }
    Ok(())
}
