//! Native declarations lower to the existing cell-covering fixed-node owner.
//! The canonical source, not a recovered candidate, defines prescribed values.
use super::*;
use std::fmt::Write as _;
use fs_topols::design_regions::{DesignPhase, PreparedDesignRegions, prepare_design_regions_controlled};

pub(super) fn parse_regions(fields: &[Node]) -> Result<Vec<DesignRegion>> {
    let Some(pair) = fields.windows(2).find(|pair|
        matches!(&pair[0].kind, NodeKind::Keyword(key) if key == "design-regions"))
    else { return Ok(Vec::new()) };
    let entries = list(&pair[1], "design-regions")?;
    if !(1..=64).contains(&entries.len()) {
        return Err(malformed("design-regions requires 1..=64 explicit rectangles"));
    }
    let mut regions = Vec::with_capacity(entries.len());
    for entry in entries {
        let entry = list(entry, "design region")?;
        if entry.len() != 9
            || !matches!(&entry[0].kind, NodeKind::Symbol(value) if value == "region")
            || ["phase", "lower", "upper", "phi-margin"].iter().enumerate()
                .any(|(i, key)| !matches!(&entry[1 + 2 * i].kind, NodeKind::Keyword(value) if value.as_str() == *key))
        {
            return Err(malformed("each region requires (region :phase material|void :lower (x y) :upper (x y) :phi-margin value), in that order"));
        }
        let phase = match &entry[2].kind {
            NodeKind::Symbol(value) if value == "material" => DesignPhase::Material,
            NodeKind::Symbol(value) if value == "void" => DesignPhase::Void,
            _ => return Err(malformed("design-region phase must be material or void")),
        };
        let lower = super::super::super::pair(&entry[4], "region lower")?;
        let upper = super::super::super::pair(&entry[6], "region upper")?;
        let margin = super::super::super::number(&entry[8], "region phi-margin")?;
        regions.push(DesignRegion::new(phase, lower, upper, margin)
            .map_err(|error| malformed(&error.to_string()))?);
    }
    Ok(regions)
}

fn phase(region: DesignRegion) -> &'static str {
    match region.phase() { DesignPhase::Material => "material", DesignPhase::Void => "void" }
}

// This optional final optimizer field closes both its list and the optimizer.
// Omitting it retains the original canonical source byte-for-byte.
pub(super) fn canonical(regions: &[DesignRegion], out: &mut String) {
    let _ = writeln!(out, "    :design-regions (");
    for &region in regions {
        let [x0, y0] = region.lower();
        let [x1, y1] = region.upper();
        let _ = writeln!(out, "      (region :phase {} :lower ({} {}) :upper ({} {}) :phi-margin {})",
            phase(region), canonical_float(x0), canonical_float(y0),
            canonical_float(x1), canonical_float(y1), canonical_float(region.margin()));
    }
    let _ = writeln!(out, "    ))");
}

fn json(regions: &[DesignRegion]) -> String {
    let rows: Vec<_> = regions.iter().map(|&region| {
        let [x0, y0] = region.lower();
        let [x1, y1] = region.upper();
        format!("{{\"phase\":\"{}\",\"lower\":[{x0:.17e},{y0:.17e}],\"upper\":[{x1:.17e},{y1:.17e}],\"phi_margin\":{:.17e}}}",
            phase(region), region.margin())
    }).collect();
    format!("[{}]", rows.join(","))
}

pub(super) fn json_field(regions: &[DesignRegion]) -> String {
    if regions.is_empty() { String::new() }
    else { format!(",\"design_regions\":{}", json(regions)) }
}

pub(super) fn check_retained(value: &JsonValue, regions: &[DesignRegion]) -> Result<()> {
    let retained = value.get("design_regions");
    if regions.is_empty() {
        if retained.is_none() { return Ok(()) }
    } else if retained == Some(&document(json(regions).as_bytes())?) {
        return Ok(());
    }
    Err(malformed("retained design regions differ from the canonical study"))
}

pub(super) fn prepare<B>(spec: &ElasticitySpec, regions: &[DesignRegion],
    mut control: impl FnMut(DesignRegionStage) -> ControlFlow<B>,
) -> Result<ControlFlow<B, PreparedDesignRegions>> {
    if let ControlFlow::Break(reason) = control(DesignRegionStage::Prepare) {
        return Ok(ControlFlow::Break(reason));
    }
    let initial = initial_phi(spec);
    let fixed: Vec<_> = initial.nodes().iter().copied().enumerate()
        .filter(|(i, _)| i % (initial.n() + 1) == 0 || i % (initial.n() + 1) == initial.n()).collect();
    // Strictly interior holes leave complete material boundary traces. Preserve
    // their original nodal bits; conflicting regions must never erase a load.
    if fixed.iter().any(|(_, value)| !value.is_finite() || *value >= 0.0) {
        return Err(malformed("declared plate must retain material on both fixed boundary traces"));
    }
    prepare_design_regions_controlled(&initial, &fixed, regions, control)
        .map_err(|error| malformed(&error.to_string()))
}
