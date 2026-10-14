//! Ordered interface packs lowered through the shared electrical-interface owner.
//! Material data and query semantics stay owned by fs-matdb/fs-material.

use std::collections::BTreeSet;
use std::path::Path;

use fs_exec::CancelGate;
use fs_lbm::conjugate::{ChtError, VoxelDomain};
use fs_lbm::conjugate::material_coupling::interfaces::MaterialContact;
use fs_matdb::{NormalizedInterfacePack, QueryPoint};
use fs_material::electrical_interface::resolve_electrical_interface_state;
use fs_material::state_point::MaterialPropertySelection;
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_qty::{QuantitySpec, Temperature};

use super::source::{Problem, Source, invalid, poll};
use super::{J, MAX_CELLS, MAX_PACK_BYTES, MAX_TOTAL_PACK_BYTES, array, fields,
    integer, quantity, quote, read, scalar, state_point, string};

pub(super) struct InterfaceSource {
    pack: NormalizedInterfacePack,
    faces: Vec<[usize; 2]>,
    point: QueryPoint,
    temperature_axes: [String; 2],
    temperature_quantity: Option<QuantitySpec>,
    heat_fraction_to_a: f64,
}

/// Binding a pack explicitly assigns its ordered surface/finish descriptors
/// to these faces. Actual bulk states are checked here; finishes are declared,
/// not measured or inferred. Temperatures are adjacent cell-state samples,
/// not microscopic interface temperatures; the source must support that model.
#[allow(clippy::too_many_arguments, clippy::too_many_lines)]
pub(super) fn load(
    root: &J,
    scene_path: &Path,
    domain: &VoxelDomain,
    materials: &[Source],
    material_indices: &[usize],
    total_bytes: &mut u64,
) -> Result<Vec<InterfaceSource>, String> {
    if root.get("electrical_interfaces").is_none() { return Ok(Vec::new()); }
    let mut sources = Vec::new();
    let mut used_faces = BTreeSet::new();
    for definition in array(root, "electrical_interfaces", 64)? {
        fields(definition, &["pack", "faces", "temperature_axes",
            "temperature_kind", "state", "heat_fraction_to_a"])?;
        let path = scene_path.parent().unwrap_or(Path::new(".")).join(string(definition, "pack")?);
        let bytes = read(&path, MAX_PACK_BYTES.min(MAX_TOTAL_PACK_BYTES - *total_bytes))?;
        *total_bytes += bytes.len() as u64;
        let pack = NormalizedInterfacePack::from_bytes(&bytes).map_err(|e| e.to_string())?;
        let fraction = scalar(definition, "heat_fraction_to_a")?;
        if !(0.0..=1.0).contains(&fraction) {
            return Err("heat_fraction_to_a must be in [0,1]".into());
        }
        let axes = array(definition, "temperature_axes", 2)?;
        if axes.len() != 2 { return Err("temperature_axes needs distinct A and B cell-state axes".into()); }
        let mut names = [String::new(), String::new()];
        for (name, value) in names.iter_mut().zip(axes) {
            *name = value.as_str().filter(|s| !s.trim().is_empty())
                .ok_or("temperature axes must be nonempty strings")?.to_owned();
        }
        let point = state_point(definition)?;
        if names[0] == names[1] || names.iter().any(|axis| point.axes().contains_key(axis)) {
            return Err("two temperature axes must be distinct and absent from the fixed state".into());
        }
        let temperature_quantity = quantity(string(definition, "temperature_kind")?, Temperature::DIMS)?;
        let absolute = QuantitySpec::semantic(SemanticType::new(QuantityKind::AbsoluteTemperature, ValueForm::Static));
        if temperature_quantity.is_some_and(|q| q != absolute && q != QuantitySpec::dimensional(Temperature::DIMS)) {
            return Err("interface temperatures must use absolute kelvin, not temperature differences".into());
        }
        let input_faces = array(definition, "faces", MAX_CELLS)?;
        if input_faces.is_empty() { return Err("an interface must bind at least one face".into()); }
        let mut faces = Vec::with_capacity(input_faces.len());
        for input in input_faces {
            let pair = input.as_array().filter(|p| p.len() == 2).ok_or("an interface face needs two cell indices")?;
            let cells = [integer(&pair[0])?, integer(&pair[1])?];
            if cells.iter().any(|&c| c >= domain.cell_count()) {
                return Err("interface cell lies outside the domain".into());
            }
            let (a, b) = (domain.coords(cells[0]), domain.coords(cells[1]));
            if (0..3).map(|axis| a[axis].abs_diff(b[axis])).sum::<usize>() != 1 {
                return Err("interface cells must share exactly one voxel face".into());
            }
            if !used_faces.insert((cells[0].min(cells[1]), cells[0].max(cells[1]))) {
                return Err("an electrical interface face is bound more than once".into());
            }
            if used_faces.len() > MAX_CELLS { return Err("electrical interfaces exceed the face budget".into()); }
            for (surface, &cell) in [pack.card().surface_a(), pack.card().surface_b()].iter().zip(&cells) {
                let material = &materials[material_indices[cell]].pack;
                if &surface.material != material.card().id() {
                    return Err(format!("interface surface order/material state does not match cell {cell}"));
                }
            }
            faces.push(cells);
        }
        sources.push(InterfaceSource {
            pack, faces, point, temperature_axes: names, temperature_quantity,
            heat_fraction_to_a: fraction,
        });
    }
    Ok(sources)
}

pub(super) fn resolve(
    sources: &[InterfaceSource], temperatures: &[f64], gate: &CancelGate,
) -> Result<Vec<MaterialContact>, ChtError> {
    let mut resolved = Vec::new();
    for source in sources {
        for &cells in &source.faces {
            if resolved.len().is_multiple_of(256) { poll(gate)?; }
            let mut point = source.point.clone();
            for (axis, &cell) in source.temperature_axes.iter().zip(&cells) {
                let temperature = temperatures.get(cell).copied()
                    .filter(|t| t.is_finite() && *t > 0.0)
                    .ok_or_else(|| invalid("interface.temperature", "missing or invalid cell-state temperature"))?;
                point = match source.temperature_quantity {
                    Some(q) => point.with_quantity(axis, q, temperature),
                    None => point.with(axis, temperature),
                }.map_err(|e| invalid("interface.state", e))?;
            }
            // The material owner alone defines the canonical property, units,
            // source selection and admissibility. This layer binds geometry.
            let state = resolve_electrical_interface_state(
                source.pack.card(), &point, MaterialPropertySelection::SingleClaimOnly,
            ).map_err(|e| invalid("interface.electrical_resistance", format!("cells {cells:?}: {e}")))?;
            resolved.push(MaterialContact::try_new(state, cells,
                [source.pack.card().surface_a(), source.pack.card().surface_b()],
                source.heat_fraction_to_a)?);
        }
    }
    poll(gate)?;
    Ok(resolved)
}

pub(super) fn render(problem: &Problem, resolved: &[MaterialContact]) -> String {
    // The immutable bindings always emit the same contact order/count; the
    // coupled owner retains exactly its final successful electrical input.
    let count: usize = problem.electrical_interfaces.iter().map(|s| s.faces.len()).sum();
    assert_eq!(count, resolved.len(), "contact report must cover every source binding");
    let definitions = problem.electrical_interfaces.iter()
        .flat_map(|source| source.faces.iter().map(move |_| source));
    let rows = definitions.zip(resolved).map(|(source, entry)| {
        let face = entry.face();
        format!(concat!("{{\"cells\":{:?},\"resistance_ohm_m2\":{},\"heat_fraction_to_a\":{},",
            "\"pack_identity\":{},\"interface_card_identity\":{},\"state_identity\":{}}}"),
            face.cells, face.resistance_ohm_m2, face.heat_fraction_to_first,
            quote(&source.pack.content_hash().to_hex()), quote(&source.pack.card().content_hash().to_hex()),
            quote(&entry.source().resolved().identity().to_hex()))
    }).collect::<Vec<_>>().join(",");
    format!("[{rows}]")
}
