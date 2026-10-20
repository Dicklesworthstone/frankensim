//! Ordered interface packs shared by electrical and thermal response bindings.
//! Geometry, surface roles and state coordinates have one parser; distinct
//! constitutive responses retain their existing fs-material/fs-conduction owners.

use std::collections::BTreeSet;
use std::path::Path;

use fs_conduction::interface::{AREA_SPECIFIC_THERMAL_RESISTANCE_PROPERTY, InterfaceResistance};
use fs_exec::CancelGate;
use fs_lbm::conjugate::{ChtError, EnergySolution, FaceContactResistance, VoxelDomain};
use fs_lbm::conjugate::material_coupling::interfaces::MaterialContact;
use fs_matdb::{ClaimId, NormalizedInterfacePack, QueryPoint};
use fs_material::electrical_interface::resolve_electrical_interface_state;
use fs_material::state_point::MaterialPropertySelection;
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_qty::{QuantitySpec, Temperature};

use super::source::{Problem, ResultField, Source, invalid, poll};
use super::{J, MAX_CELLS, MAX_PACK_BYTES, MAX_TOTAL_PACK_BYTES, array, fields,
    integer, quantity, quote, read, scalar, state_point, string};

enum Response {
    Electrical { heat_fraction_to_a: f64 },
    Thermal { claim: ClaimId },
}

pub(super) struct InterfaceSource {
    pack: NormalizedInterfacePack,
    faces: Vec<[usize; 2]>,
    point: QueryPoint,
    temperature_axes: [String; 2],
    temperature_quantity: Option<QuantitySpec>,
    response: Response,
}

impl InterfaceSource {
    fn point_at(&self, cells: [usize; 2], temperatures: &[f64]) -> Result<QueryPoint, ChtError> {
        let mut point = self.point.clone();
        for (axis, cell) in self.temperature_axes.iter().zip(cells) {
            let temperature = temperatures.get(cell).copied()
                .filter(|t| t.is_finite() && *t > 0.0)
                .ok_or_else(|| invalid("interface.temperature", "missing or invalid cell-state temperature"))?;
            point = match self.temperature_quantity {
                Some(q) => point.with_quantity(axis, q, temperature),
                None => point.with(axis, temperature),
            }.map_err(|e| invalid("interface.state", e))?;
        }
        Ok(point)
    }
}

/// Binding assigns the pack's ordered material/finish descriptors to the faces.
/// Temperatures are adjacent cell-state samples, not inferred microscopic values.
pub(super) fn load(
    root: &J, scene_path: &Path, domain: &VoxelDomain, materials: &[Source],
    material_indices: &[usize], total_bytes: &mut u64,
) -> Result<Vec<InterfaceSource>, String> {
    load_kind(root, scene_path, domain, materials, material_indices, total_bytes, false)
}

pub(super) fn load_thermal(
    root: &J, scene_path: &Path, domain: &VoxelDomain, materials: &[Source],
    material_indices: &[usize], total_bytes: &mut u64,
) -> Result<Vec<InterfaceSource>, String> {
    load_kind(root, scene_path, domain, materials, material_indices, total_bytes, true)
}

#[allow(clippy::too_many_arguments, clippy::too_many_lines)]
fn load_kind(
    root: &J, scene_path: &Path, domain: &VoxelDomain, materials: &[Source],
    material_indices: &[usize], total_bytes: &mut u64, thermal: bool,
) -> Result<Vec<InterfaceSource>, String> {
    let key = if thermal { "thermal_interfaces" } else { "electrical_interfaces" };
    if root.get(key).is_none() { return Ok(Vec::new()); }
    let mut sources = Vec::new();
    let mut used_faces = BTreeSet::new();
    for definition in array(root, key, 64)? {
        let allowed: &[&str] = if thermal {
            &["pack", "faces", "temperature_axes", "temperature_kind", "state"]
        } else {
            &["pack", "faces", "temperature_axes", "temperature_kind", "state", "heat_fraction_to_a"]
        };
        fields(definition, allowed)?;
        let path = scene_path.parent().unwrap_or(Path::new(".")).join(string(definition, "pack")?);
        let remaining = MAX_TOTAL_PACK_BYTES.checked_sub(*total_bytes)
            .ok_or("interface packs exceed the aggregate byte budget")?;
        let bytes = read(&path, MAX_PACK_BYTES.min(remaining))?;
        *total_bytes += bytes.len() as u64;
        let pack = NormalizedInterfacePack::from_bytes(&bytes).map_err(|e| e.to_string())?;
        let response = if thermal {
            // Select once, independently of temperature. A trajectory may not
            // switch to a different claim when the first claim's support ends.
            let claims = pack.card().claims_for(AREA_SPECIFIC_THERMAL_RESISTANCE_PROPERTY);
            if claims.len() != 1 {
                return Err("thermal interface requires exactly one canonical thermal-resistance claim".into());
            }
            Response::Thermal { claim: claims[0].0 }
        } else {
            let fraction = scalar(definition, "heat_fraction_to_a")?;
            if !(0.0..=1.0).contains(&fraction) {
                return Err("heat_fraction_to_a must be in [0,1]".into());
            }
            Response::Electrical { heat_fraction_to_a: fraction }
        };
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
                return Err(format!("an interface face is bound more than once in {key}"));
            }
            if used_faces.len() > MAX_CELLS { return Err(format!("{key} exceeds the face budget")); }
            for (surface, &cell) in [pack.card().surface_a(), pack.card().surface_b()].iter().zip(&cells) {
                let material = &materials[material_indices[cell]].pack;
                if &surface.material != material.card().id() {
                    return Err(format!("interface surface order/material state does not match cell {cell}"));
                }
            }
            faces.push(cells);
        }
        sources.push(InterfaceSource {
            pack, faces, point, temperature_axes: names, temperature_quantity, response,
        });
    }
    Ok(sources)
}

pub(super) fn resolve(
    sources: &[InterfaceSource], temperatures: &[f64], gate: &CancelGate,
) -> Result<Vec<MaterialContact>, ChtError> {
    let mut resolved = Vec::new();
    for source in sources {
        let Response::Electrical { heat_fraction_to_a } = source.response else {
            return Err(invalid("interface.response", "thermal law cannot supply electrical resistance"));
        };
        for &cells in &source.faces {
            if resolved.len().is_multiple_of(256) { poll(gate)?; }
            let point = source.point_at(cells, temperatures)?;
            let state = resolve_electrical_interface_state(
                source.pack.card(), &point, MaterialPropertySelection::SingleClaimOnly,
            ).map_err(|e| invalid("interface.electrical_resistance", format!("cells {cells:?}: {e}")))?;
            resolved.push(MaterialContact::try_new(state, cells,
                [source.pack.card().surface_a(), source.pack.card().surface_b()], heat_fraction_to_a)?);
        }
    }
    poll(gate)?;
    Ok(resolved)
}

/// Final candidate source state, distinct from the coefficient used by the
/// preceding linear solve. The caller checks their difference before accepting.
pub(super) struct ThermalInterfaceState {
    pub cells: [usize; 2],
    pub resistance: InterfaceResistance,
}
impl ThermalInterfaceState {
    pub fn face(&self) -> FaceContactResistance {
        FaceContactResistance { cells: self.cells,
            resistance_m2_k_w: self.resistance.value_m2_k_per_w() }
    }
}

pub(super) fn resolve_thermal(
    sources: &[InterfaceSource], temperatures: &[f64], gate: &CancelGate,
) -> Result<Vec<ThermalInterfaceState>, ChtError> {
    let mut resolved = Vec::new();
    for source in sources {
        let Response::Thermal { claim } = source.response else {
            return Err(invalid("interface.response", "electrical law cannot supply thermal resistance"));
        };
        for &cells in &source.faces {
            if resolved.len().is_multiple_of(256) { poll(gate)?; }
            let point = source.point_at(cells, temperatures)?;
            let resistance = InterfaceResistance::from_card_pinned(
                source.pack.pack_id(), source.pack.card(), &point, claim,
            ).map_err(|e| invalid("interface.thermal_resistance", format!("cells {cells:?}: {e}")))?;
            resolved.push(ThermalInterfaceState { cells, resistance });
        }
    }
    poll(gate)?;
    Ok(resolved)
}

pub(super) fn thermal_defect(states: &[ThermalInterfaceState], solution: &EnergySolution) -> Result<f64, ChtError> {
    let mut defect = 0.0_f64;
    for state in states {
        let [a, b] = state.cells;
        let used = solution.face_contact_resistance_m2_k_w.get(&(a.min(b), a.max(b)))
            .ok_or_else(|| invalid("interface.thermal_consistency", "source joint is absent from the solved field"))?;
        let next = state.resistance.value_m2_k_per_w();
        let local = (next - used).abs() / next.max(*used);
        if !local.is_finite() { return Err(invalid("interface.thermal_consistency", "nonfinite contact defect")); }
        defect = defect.max(local);
    }
    Ok(defect)
}

pub(super) fn render(problem: &Problem, resolved: &[MaterialContact]) -> String {
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

pub(super) fn render_thermal(problem: &Problem, result: &ResultField) -> String {
    let count: usize = problem.thermal_interfaces.iter().map(|s| s.faces.len()).sum();
    assert_eq!(count, result.thermal_interfaces.len(), "thermal report must cover every binding");
    let definitions = problem.thermal_interfaces.iter()
        .flat_map(|source| source.faces.iter().map(move |_| source));
    let rows = definitions.zip(&result.thermal_interfaces).map(|(source, state)| {
        let [a, b] = state.cells;
        let used = result.coupled.thermal.face_contact_resistance_m2_k_w[&(a.min(b), a.max(b))];
        format!(concat!("{{\"cells\":{:?},\"resistance_m2_k_w\":{},\"solved_resistance_m2_k_w\":{},",
            "\"pack_identity\":{},\"interface_card_identity\":{},\"source_receipt_identity\":{}}}"),
            state.cells, state.resistance.value_m2_k_per_w(), used,
            quote(&source.pack.content_hash().to_hex()), quote(&state.resistance.card_identity().to_hex()),
            quote(&state.resistance.receipt().content_hash().to_hex()))
    }).collect::<Vec<_>>().join(",");
    format!("[{rows}]")
}
