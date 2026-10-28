//! Scientific fields projected from the sealed selected design, without a solve.
//! The cells are the active octree background, not a reconstructed cut surface.

use fs_viz::vtu::{CellType, DataArray, DataAssociation, DataValues, UnstructuredGrid, VtuWriter};

use super::{
    JsonValue as J, LOAD_SPACE3_DRIVER, MAX_ARTIFACT_BYTES, Result, SDF3_DRIVER, STRESS3_DRIVER,
    fail,
};

const MAX_SCALARS: usize = 1_048_576;
const MAX_VTU_BYTES: usize = 64 * 1024 * 1024;
pub(super) const SCOPE: &str = "Estimated fields on the retained selected design's active octree background in metres. Cells show their full background boxes; cut_volume_m3 is the retained integrated material-domain measure, not a clipped surface. Density thresholds and displaced views are visualization choices, not manufactured geometry or continuum certificates. Load-space displacement fields use the retained eroded material.";

pub(super) enum Projection {
    Available(String),
    Unavailable(&'static str),
}

fn bad(message: impl Into<String>) -> super::Failure {
    fail("cli-study-sdf3-field-export", message)
}

fn array<'a>(value: &'a J, name: &str) -> Result<&'a [J]> {
    value
        .get(name)
        .and_then(J::as_array)
        .ok_or_else(|| bad(format!("retained design lacks array {name}")))
}

fn finite(value: &J) -> Result<f64> {
    value
        .as_f64()
        .filter(|v| v.is_finite())
        .ok_or_else(|| bad("retained field contains a nonfinite or nonnumeric value"))
}

fn charge(remaining: &mut usize, count: usize) -> Result<()> {
    *remaining = remaining
        .checked_sub(count)
        .ok_or_else(|| bad("retained field exceeds the 1048576-scalar export envelope"))?;
    Ok(())
}

fn vectors(
    grid: &mut UnstructuredGrid,
    name: String,
    rows: &[J],
    remaining: &mut usize,
) -> Result<()> {
    let count = grid
        .points
        .len()
        .checked_mul(3)
        .ok_or_else(|| bad("displacement component count overflows"))?;
    if rows.len() != count {
        return Err(bad(
            "retained displacement does not cover every physical node",
        ));
    }
    charge(remaining, count)?;
    let values = rows.iter().map(finite).collect::<Result<Vec<_>>>()?;
    grid.add_array(DataArray {
        name,
        association: DataAssociation::PointData,
        components: 3,
        unit: Some("m".into()),
        values: DataValues::Float64(values),
    });
    Ok(())
}

fn selected<'a>(driver: &str, design: &'a J) -> Result<std::result::Result<&'a J, &'static str>> {
    match (driver, design.str_field("schema")) {
        (SDF3_DRIVER, Some("sdf3-study-design.v1")) => match design.str_field("state") {
            Some("solved") => Ok(Ok(design)),
            Some("no-solved-baseline") => Ok(Err("no-solved-baseline")),
            _ => Err(bad("unrecognized retained adaptive design state")),
        },
        (STRESS3_DRIVER, Some("sdf3-stress-design.v1"))
        | (LOAD_SPACE3_DRIVER, Some("sdf3-load-space-design.v1")) => {
            if design.str_field("state") == Some("no-solved-baseline") {
                return Ok(Err("no-solved-baseline"));
            }
            match design.get("selected_feasible") {
                Some(J::Bool(false)) => return Ok(Err("selected-design-not-feasible")),
                Some(J::Bool(true)) => {}
                _ => return Err(bad("retained selected design lacks its feasibility state")),
            }
            let selection = design.str_field("selected_design");
            let recognized = selection == Some("last-accepted")
                || selection
                    == Some(if driver == LOAD_SPACE3_DRIVER {
                        "least-material-feasible"
                    } else {
                        "best-feasible"
                    });
            if !recognized {
                return Err(bad("retained design has no recognized selected incumbent"));
            }
            design
                .get("selected")
                .filter(|v| v.as_object().is_some())
                .map(Ok)
                .ok_or_else(|| bad("retained selected design fields are missing"))
        }
        _ => Err(bad(
            "retained design schema does not match its study driver",
        )),
    }
}

/// No solver feature is needed to view an already retained physical field.
pub(super) fn render(driver: &str, bytes: &[u8]) -> Result<Projection> {
    if bytes.len() as u64 > MAX_ARTIFACT_BYTES {
        return Err(bad("retained design exceeds its artifact envelope"));
    }
    let design = J::parse(std::str::from_utf8(bytes).map_err(|e| bad(e.to_string()))?)
        .map_err(|e| bad(e.to_string()))?;
    let fields = match selected(driver, &design)? {
        Ok(fields) => fields,
        Err(reason) => return Ok(Projection::Unavailable(reason)),
    };
    let nodes = array(fields, "physical_nodes_m")?;
    let cells = array(fields, "cells")?;
    let connectivity = array(fields, "cell_nodes")?;
    if nodes.is_empty() || cells.is_empty() || connectivity.len() != cells.len() {
        return Err(bad(
            "retained physical nodes, cells and connectivity do not form a nonempty mesh",
        ));
    }
    let mut remaining = MAX_SCALARS;
    charge(
        &mut remaining,
        nodes
            .len()
            .checked_mul(3)
            .ok_or_else(|| bad("node count overflows"))?,
    )?;
    charge(
        &mut remaining,
        cells
            .len()
            .checked_mul(8)
            .ok_or_else(|| bad("cell count overflows"))?,
    )?;
    let mut grid = UnstructuredGrid::new();
    for node in nodes {
        let xyz = node
            .as_array()
            .filter(|p| p.len() == 3)
            .ok_or_else(|| bad("physical node must have exactly three coordinates"))?;
        grid.add_point(finite(&xyz[0])?, finite(&xyz[1])?, finite(&xyz[2])?);
    }
    for cell in connectivity {
        let corners = cell
            .as_array()
            .filter(|p| p.len() == 8)
            .ok_or_else(|| bad("retained octree cell must have eight corners"))?;
        let mut ids = [0usize; 8];
        for (id, value) in ids.iter_mut().zip(corners) {
            *id = value
                .number_raw()
                .and_then(|raw| raw.parse::<usize>().ok())
                .filter(|id| *id < grid.points.len())
                .ok_or_else(|| {
                    bad("retained cell corner is not an in-range physical node index")
                })?;
        }
        // fs-cutfem::Octree3::corners uses x/y/z endpoint bits. This is
        // VTK_VOXEL order; VTK_HEXAHEDRON would swap corners 2/3 and 6/7.
        let lo = grid.points[ids[0]];
        let hi = grid.points[ids[7]];
        for axis in 0..3 {
            if lo[axis] >= hi[axis] {
                return Err(bad("retained octree cell has collapsed or reversed bounds"));
            }
            for (corner, &id) in ids.iter().enumerate() {
                let expected = if corner & (1 << axis) == 0 {
                    lo[axis]
                } else {
                    hi[axis]
                };
                if grid.points[id][axis] != expected {
                    return Err(bad(
                        "retained cell is not an axis-aligned octree voxel in canonical corner order",
                    ));
                }
            }
        }
        grid.add_cell(CellType::Voxel, &ids);
    }
    let densities: &[&str] = if driver == LOAD_SPACE3_DRIVER {
        &[
            "raw_density",
            "eroded_density",
            "nominal_density",
            "dilated_density",
        ]
    } else {
        &["raw_density", "projected_density"]
    };
    for &name in densities.iter().chain(std::iter::once(&"cut_volume_m3")) {
        charge(&mut remaining, cells.len())?;
        let values = cells
            .iter()
            .map(|cell| {
                let value = finite(
                    cell.get(name)
                        .ok_or_else(|| bad(format!("retained cell lacks {name}")))?,
                )?;
                if value < 0.0 || (name != "cut_volume_m3" && value > 1.0) {
                    return Err(bad(format!(
                        "retained {name} is outside its physical range"
                    )));
                }
                Ok(value)
            })
            .collect::<Result<Vec<_>>>()?;
        grid.add_array(
            DataArray::new_cell_scalar(name, values).with_unit(if name == "cut_volume_m3" {
                "m^3"
            } else {
                "1"
            }),
        );
    }
    let (field_name, prefix, layout) = if driver == LOAD_SPACE3_DRIVER {
        (
            "basis_displacements_m",
            "basis_displacement",
            "basis column then node-major xyz, on eroded material",
        )
    } else {
        (
            "displacements_m",
            "displacement_case",
            "load-case then node-major xyz",
        )
    };
    if fields.str_field("displacement_layout") != Some(layout) {
        return Err(bad(
            "retained displacement layout differs from the declared study schema",
        ));
    }
    let cases = array(fields, field_name)?;
    if cases.is_empty() || cases.len() > 64 {
        return Err(bad("retained displacement case count must be in 1..=64"));
    }
    for (index, field) in cases.iter().enumerate() {
        let values = field
            .as_array()
            .ok_or_else(|| bad("retained displacement case is not an array"))?;
        vectors(
            &mut grid,
            format!("{prefix}_{index}"),
            values,
            &mut remaining,
        )?;
    }
    if driver == LOAD_SPACE3_DRIVER {
        vectors(
            &mut grid,
            "worst_displacement".into(),
            array(fields, "worst_displacement_m")?,
            &mut remaining,
        )?;
    } else if driver == STRESS3_DRIVER {
        let stress = array(fields, "cell_relaxed_max_pa")?;
        if stress.len() != cells.len() {
            return Err(bad("retained sampled stress does not cover every cell"));
        }
        charge(&mut remaining, stress.len())?;
        let values = stress
            .iter()
            .map(|v| {
                let value = finite(v)?;
                if value < 0.0 {
                    return Err(bad("retained sampled stress is negative"));
                }
                Ok(value)
            })
            .collect::<Result<Vec<_>>>()?;
        grid.add_array(
            DataArray::new_cell_scalar("sampled_relaxed_max_stress", values).with_unit("Pa"),
        );
    }
    let xml = VtuWriter::write_ascii(&grid).map_err(|e| bad(e.to_string()))?;
    if xml.len() > MAX_VTU_BYTES {
        return Err(bad("VTU exceeds the 64 MiB export envelope"));
    }
    Ok(Projection::Available(xml))
}

#[cfg(test)]
#[path = "field_export/tests.rs"]
mod tests;
