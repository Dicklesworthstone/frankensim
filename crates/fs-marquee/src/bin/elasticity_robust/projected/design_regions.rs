//! Bounded authored non-design regions, baked into existing fixed-node state.
use super::{numbers, write_field, writer};
use fs_topols::design_regions::{DesignPhase, DesignRegion, PreparedDesignRegions, prepare_design_regions};
use fs_topols::GridSdf;
use std::error::Error;
use std::fs::File;
use std::io::{Read, Write};
use std::path::Path;

pub(super) fn options(args: &[String]) -> Result<(Vec<String>, Option<String>), Box<dyn Error>> {
    if args.len() > 32 { return Err("too many projected study arguments".into()); }
    let mut remaining = Vec::new();
    let mut path = None;
    let mut index = 0;
    while index < args.len() {
        if args[index] == "--design-regions" {
            if path.is_some() { return Err("duplicate --design-regions".into()); }
            let value = args.get(index + 1).ok_or("--design-regions requires a CSV path")?;
            if value.is_empty() || value.starts_with("--") {
                return Err("--design-regions requires a CSV path, not another option".into());
            }
            path = Some(value.clone());
            index += 2;
        } else {
            remaining.push(args[index].clone());
            index += 1;
        }
    }
    Ok((remaining, path))
}

fn parse(text: &str) -> Result<Vec<DesignRegion>, Box<dyn Error>> {
    let mut regions = Vec::new();
    for (line_number, line) in text.lines().enumerate() {
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') { continue; }
        if regions.len() == 64 { return Err("at most 64 design regions are admitted".into()); }
        let fields: Vec<_> = line.split(',').map(str::trim).collect();
        if fields.len() != 6 {
            return Err(format!("design-region line {} requires phase,x_min,y_min,x_max,y_max,phi_margin", line_number + 1).into());
        }
        let phase = match fields[0] {
            "material" => DesignPhase::Material,
            "void" => DesignPhase::Void,
            _ => return Err(format!("design-region line {} phase must be material or void", line_number + 1).into()),
        };
        regions.push(DesignRegion::new(phase,
            [fields[1].parse()?, fields[2].parse()?],
            [fields[3].parse()?, fields[4].parse()?], fields[5].parse()?)?);
    }
    if regions.is_empty() { return Err("design-region CSV contains no regions".into()); }
    Ok(regions)
}

pub(super) struct Authoring {
    records: Vec<DesignRegion>,
    pub prepared: PreparedDesignRegions,
}

pub(super) fn load(
    path: &Path, field: &GridSdf, fixed: &[(usize, f64)],
) -> Result<Authoring, Box<dyn Error>> {
    const MAX_BYTES: u64 = 1_048_576;
    let mut text = String::new();
    File::open(path)?.take(MAX_BYTES + 1).read_to_string(&mut text)?;
    if text.len() as u64 > MAX_BYTES { return Err("design-region CSV exceeds 1 MiB".into()); }
    let records = parse(&text)?;
    let prepared = prepare_design_regions(field, fixed, &records)?;
    Ok(Authoring { records, prepared })
}

impl Authoring {
    /// Called by the original runner only after the complete baseline is admitted.
    pub fn export(&self, output: &Path) -> Result<String, Box<dyn Error>> {
        let mut file = writer(&output.join("design-regions.csv"))?;
        writeln!(file, "# phase,x_min,y_min,x_max,y_max,phi_margin")?;
        let mut entries = Vec::with_capacity(self.records.len());
        for (record, cells) in self.records.iter().zip(&self.prepared.covered_cells) {
            let phase = match record.phase() { DesignPhase::Material => "material", DesignPhase::Void => "void" };
            let [x0, y0] = record.lower();
            let [x1, y1] = record.upper();
            let margin = record.margin();
            writeln!(file, "{phase},{x0:.17e},{y0:.17e},{x1:.17e},{y1:.17e},{margin:.17e}")?;
            entries.push(format!(
                "{{\"phase\":\"{phase}\",\"lower\":{},\"upper\":{},\"phi_margin\":{margin:.17e},\"covered_cells\":[{},{},{},{}]}}",
                numbers(&record.lower()), numbers(&record.upper()), cells[0], cells[1], cells[2], cells[3],
            ));
        }
        file.flush()?;
        write_field(&output.join("design-region-level-set.csv"), &self.prepared.geometry)?;
        Ok(format!(
            "{{\"schema\":\"fixed-design-regions-v1\",\"units\":\"normalized_geometry_and_phi\",\"coverage\":\"whole_intersected_cells\",\"material_nodes\":{},\"void_nodes\":{},\"total_fixed_nodes\":{},\"changed_input_nodes\":{},\"manufacturing_certificate\":false,\"regions\":[{}]}}",
            self.prepared.material_nodes, self.prepared.void_nodes,
            self.prepared.fixed_nodes.len(), self.prepared.changed_nodes, entries.join(","),
        ))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn options_do_not_consume_stress_or_pause_and_reject_duplicates_before_reads() {
        let args = ["out", "loads.csv", "--stress-limit", "100", "--design-regions", "regions.csv", "--pause-after", "0"].map(str::to_string);
        let (rest, path) = options(&args).unwrap();
        assert_eq!(path.as_deref(), Some("regions.csv"));
        assert_eq!(rest, ["out", "loads.csv", "--stress-limit", "100", "--pause-after", "0"]);
        for bad in [vec!["--design-regions"], vec!["--design-regions", "--checkpoint"],
            vec!["--design-regions", "a", "--design-regions", "b"]]
        { assert!(options(&bad.into_iter().map(str::to_string).collect::<Vec<_>>()).is_err()); }
    }

    #[test]
    fn explicit_phase_bounds_margin_and_record_limits_are_checked() {
        let row = "void,0.5,0.5,0.625,0.625,0.02\n";
        let records = parse(&format!("# clearance\n{row}")).unwrap();
        assert_eq!(records[0].phase(), DesignPhase::Void);
        assert!(parse(&row.repeat(64)).is_ok());
        assert!(parse(&row.repeat(65)).is_err());
        for bad in ["# empty", "void,0,0,1,1", "solid,0,0,1,1,0.1",
            "void,0,0,1,1,NaN", "void,0,0,1,1,-0.1", "void,0.5,0.5,0.4,0.6,0.1"]
        { assert!(parse(bad).is_err(), "{bad}"); }
    }
}
