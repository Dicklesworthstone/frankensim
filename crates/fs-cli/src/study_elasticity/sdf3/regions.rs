//! SI-box authoring of whole-cell, non-design material regions.
//! Boundaries must coincide with the INITIAL octree lattice; no centroid mask
//! or partial-cell rasterization silently changes the requested region. The
//! core adaptive owner inherits labels after this initial binding.
use super::*;
use fs_ir::ast::{Node, NodeKind};
use fs_topopt::sdf3::{AdaptiveSdf3Elasticity, PhysicalRegion3};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(super) struct RegionBox {
    pub lower: [u32; 3],
    pub upper: [u32; 3],
    pub kind: PhysicalRegion3,
}
fn invalid(message: impl Into<String>) -> Failure {
    fail("cli-study-sdf3-regions", message)
}
fn vector(node: &Node) -> Result<[f64; 3]> {
    let values = super::super::list(node, "physical region coordinate")?;
    if values.len() != 3 { return Err(invalid("a region coordinate requires three SI lengths")); }
    Ok([
        super::super::number(&values[0], "region x")?,
        super::super::number(&values[1], "region y")?,
        super::super::number(&values[2], "region z")?,
    ])
}
fn lattice_point(value: f64, axis: usize, spec: &Spec) -> Result<u32> {
    let cells = 1_u32 << spec.level;
    let lo = spec.bounds.0[axis];
    let hi = spec.bounds.1[axis];
    if !value.is_finite() { return Err(invalid("physical region coordinates must be finite")); }
    (0..=cells).find(|&i| {
        let plane = if i == cells { hi } else { lo + (hi - lo) * (f64::from(i) / f64::from(cells)) };
        value == plane
    }).ok_or_else(|| invalid("region boundaries must lie on initial octree planes inside the declared SI bounds"))
}

pub(super) fn parse(node: &Node, spec: &Spec) -> Result<Vec<RegionBox>> {
    if !(1..=2).contains(&spec.level) { return Err(invalid("region authoring requires an admitted initial octree")); }
    let values = super::super::list(node, "design-regions")?;
    if values.len() != 5
        || !matches!(&values[0].kind, NodeKind::Symbol(s) if s == "design-regions")
        || !matches!(&values[1].kind, NodeKind::Keyword(s) if s == "solid")
        || !matches!(&values[3].kind, NodeKind::Keyword(s) if s == "void")
    {
        return Err(invalid("expected (design-regions :solid (BOX ...) :void (BOX ...)), with no extra fields"));
    }
    let mut regions = Vec::new();
    for (node, kind) in [(&values[2], PhysicalRegion3::Solid), (&values[4], PhysicalRegion3::Void)] {
        let boxes = super::super::list(node, "region box list")?;
        if regions.len() + boxes.len() > 32 { return Err(invalid("at most 32 physical region boxes are admitted")); }
        for node in boxes {
            let points = super::super::list(node, "region box")?;
            if points.len() != 2 { return Err(invalid("a region box requires exactly lower and upper SI coordinates")); }
            let lo = vector(&points[0])?;
            let hi = vector(&points[1])?;
            let mut lower = [0; 3]; let mut upper = [0; 3];
            for axis in 0..3 {
                lower[axis] = lattice_point(lo[axis], axis, spec)?;
                upper[axis] = lattice_point(hi[axis], axis, spec)?;
                if lower[axis] >= upper[axis] { return Err(invalid("physical region boxes must have positive ordered spans")); }
            }
            let region = RegionBox { lower, upper, kind };
            if regions.iter().any(|old: &RegionBox| old.kind != kind && (0..3).all(|axis|
                old.lower[axis] < upper[axis] && lower[axis] < old.upper[axis])) {
                return Err(invalid("solid and void region interiors overlap"));
            }
            regions.push(region);
        }
    }
    Ok(regions)
}

pub(super) fn bind<O: AdaptiveSdf3Elasticity>(
    study: CutDensityStudy3<O>, spec: &Spec,
) -> Result<CutDensityStudy3<O>> {
    if spec.regions.is_empty() { return Ok(study); }
    let leaves = study.operator().adaptive().leaves();
    let mut labels = vec![PhysicalRegion3::Design; leaves.len()];
    for region in &spec.regions {
        let mut hits = 0;
        for (index, leaf) in leaves.iter().enumerate() {
            let level = u32::from(leaf.level());
            if level < spec.level { return Err(invalid("region binding cannot coarsen the declared initial lattice")); }
            let ancestor = leaf.index().map(|i| i >> (level - spec.level));
            if (0..3).all(|axis| region.lower[axis] <= ancestor[axis] && ancestor[axis] < region.upper[axis]) {
                if labels[index] != PhysicalRegion3::Design && labels[index] != region.kind {
                    return Err(invalid("contradictory physical region assignments"));
                }
                labels[index] = region.kind;
                hits += 1;
            }
        }
        if hits == 0 { return Err(invalid("a declared region contains no active material cells; adjust its box or implicit domain")); }
    }
    let study = study.with_physical_regions(labels)
        .map_err(|e| invalid(format!("cannot bind physical regions: {e:?}")))?;
    if spec.stress.is_none() && study.prescribed_solid_fraction() > spec.volume + 1e-8 {
        return Err(invalid("prescribed solid cells alone exceed the material-volume allowance"));
    }
    Ok(study)
}

#[cfg(test)]
mod tests;
