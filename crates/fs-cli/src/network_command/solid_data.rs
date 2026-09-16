//! Declared heterogeneous materials and component heating for cooling-network.
//! Reuses fs-conduction's checked element assignment and power-preserving P1
//! source projection. Nodal component footprints are not sharp cellwise sources:
//! their support extends over incident tetrahedra, including material interfaces.
//! Constant isotropic/tensor/orthotropic materials and bounded scalar k(T)
//! curves share the same FEM and derivative owners. No scalar averaging or
//! frozen-temperature substitution occurs.

mod constitutive;

use constitutive::Conductivity;
use fs_conduction::{ComponentPower, ConductivityModel, ElementMaterials, MaterialId,
    MaterialTable, PowerAudit, PowerMap, PowerUncertainty};
use super::*;

#[derive(Debug)]
struct MaterialDeclaration { name: String, conductivity: Conductivity, source: String }

#[derive(Debug)]
pub(super) struct SolidData {
    pub element_materials: Option<ElementMaterials>,
    pub nodal_source: Option<ScalarField>,
    pub power: Option<PowerAudit>,
    pub component_map: Option<PowerMap>,
    materials: Vec<MaterialDeclaration>,
}

impl SolidData {
    /// Exactly one material mode and one heating mode. Legacy scalar inputs
    /// keep the same uniform operator; mixed spellings never silently override.
    pub fn parse(solid: &J, mesh: &ConductionMesh) -> Result<(Self, f64, f64)> {
        let mut data = Self { element_materials: None, nodal_source: None,
            power: None, component_map: None, materials: Vec::new() };
        let uniform_k = solid.get("conductivity_w_m_k");
        let table = solid.get("materials");
        let assignment = solid.get("element_materials");
        let fallback = match (uniform_k, table, assignment) {
            (Some(k), None, None) => positive(k, "conductivity_w_m_k")?,
            (None, Some(table), Some(assignment)) => {
                let mut declarations = BTreeMap::new();
                for row in array(table, "materials", 4096)? {
                    object(row, &["name", "conductivity_w_m_k", "conductivity_tensor_w_m_k", "orthotropic", "conductivity_curve", "source"], "material")?;
                    let name = string(get(row, "name")?, "material.name")?;
                    let conductivity = Conductivity::parse(row)?;
                    let source = string(get(row, "source")?, "material.source")?;
                    if declarations.insert(name.clone(), MaterialDeclaration { name, conductivity, source }).is_some() {
                        return Err(bad("duplicate material name"));
                    }
                }
                if declarations.is_empty() { return Err(bad("materials must not be empty")); }
                data.materials = declarations.into_values().collect();
                let ids: BTreeMap<_, _> = data.materials.iter().enumerate()
                    .map(|(i, m)| (m.name.as_str(), MaterialId(i as u32))).collect();
                let rows = array(assignment, "element_materials", mesh.element_count())?;
                if rows.len() != mesh.element_count() {
                    return Err(bad("element_materials requires exactly one material name per tetrahedron"));
                }
                let of_element = rows.iter().map(|row| {
                    let name = string(row, "element material")?;
                    ids.get(name.as_str()).copied().ok_or_else(|| bad(format!("unknown element material {name}")))
                }).collect::<Result<Vec<_>>>()?;
                let entries = data.materials.iter().enumerate().map(|(i, material)| {
                    Ok((MaterialId(i as u32), material.conductivity.model()?))
                }).collect::<Result<Vec<_>>>()?;
                let assigned = ElementMaterials::new(MaterialTable::new(entries).map_err(producer)?, of_element).map_err(producer)?;
                assigned.validate_for(mesh).map_err(producer)?;
                data.element_materials = Some(assigned);
                // The fallback slot is ignored whenever element assignment is
                // present. It is not a homogenized or averaged conductivity;
                // the actual operator always receives the full material law.
                data.materials[0].conductivity.inactive_scalar()
            }
            _ => return Err(bad("use either conductivity_w_m_k or both materials and element_materials, never a mixture")),
        };
        let source = match (solid.get("source_w_m3"), solid.get("component_power")) {
            (Some(value), None) => number(value, "source_w_m3")?,
            (None, Some(value)) => {
                object(value, &["total_w", "relative_tolerance", "components"], "component_power")?;
                let total = number(get(value, "total_w")?, "component_power.total_w")?;
                let tolerance = number(get(value, "relative_tolerance")?, "component_power.relative_tolerance")?;
                if !(0.0..1.0).contains(&tolerance) { return Err(bad("power relative_tolerance must be in [0,1)")); }
                let mut components = Vec::new();
                for row in array(get(value, "components")?, "components", 4096)? {
                    object(row, &["name", "watts", "vertices"], "component")?;
                    let name = string(get(row, "name")?, "component.name")?;
                    let watts = number(get(row, "watts")?, "component.watts")?;
                    let mut vertices = BTreeSet::new();
                    for vertex in array(get(row, "vertices")?, "component.vertices", mesh.vertex_count())? {
                        let vertex = integer(vertex, "component vertex", mesh.vertex_count() - 1)?;
                        if !vertices.insert(vertex) { return Err(bad(format!("component {name} repeats vertex {vertex}"))); }
                    }
                    components.push(ComponentPower::new(name, watts, PowerUncertainty::Unstated,
                        vertices.into_iter().collect()).map_err(producer)?);
                }
                let map = PowerMap::new(components, total).map_err(producer)?;
                let (source, audit) = map.volumetric_source(mesh, tolerance).map_err(producer)?;
                source.validate("component power source", mesh.vertex_count()).map_err(producer)?;
                data.nodal_source = Some(source);
                data.power = Some(audit);
                data.component_map = Some(map);
                // No additional uniform background source was declared.
                0.0
            }
            _ => return Err(bad("use exactly one of source_w_m3 and component_power")),
        };
        Ok((data, fallback, source))
    }

    /// Preserve declarations and delivered powers alongside the solved field.
    /// This is caller-supplied data, not material-card or uncertainty evidence.
    pub fn render(&self, uniform_k: f64, uniform_source: f64) -> Result<String> {
        let materials = if let Some(assignment) = &self.element_materials {
            let rows = self.materials.iter().map(|m| Ok(format!(
                "{{\"name\":{},{},\"source\":{}}}",
                quote(&m.name), m.conductivity.render_fields()?, quote(&m.source))))
                .collect::<Result<Vec<_>>>()?.join(",");
            let names = assignment.of_element().iter().map(|id| quote(&self.materials[id.0 as usize].name))
                .collect::<Vec<_>>().join(",");
            format!("{{\"mode\":\"element-materials\",\"authority\":\"caller-declared\",\"materials\":[{rows}],\"element_materials\":[{names}]}}")
        } else { format!("{{\"mode\":\"uniform\",\"authority\":\"caller-declared\",\"conductivity_w_m_k\":{}}}", num(uniform_k)?) };
        let heating = if let Some(audit) = &self.power {
            let rows = audit.rows().iter().map(|row| Ok(format!(
                "{{\"name\":{},\"declared_w\":{},\"delivered_w\":{},\"bound_volume_m3\":{},\"density_w_m3\":{}}}",
                quote(row.name()), num(row.declared_w())?, num(row.delivered_w())?,
                num(row.bound_volume_m3())?, num(row.density_w_per_m3())?)))
                .collect::<Result<Vec<_>>>()?.join(",");
            format!("{{\"mode\":\"component-power\",\"projection\":\"lumped-volume-normalized nodal P1 support\",\"declared_total_w\":{},\"delivered_total_w\":{},\"uncertainty_w\":null,\"components\":[{rows}]}}",
                num(audit.declared_total_w())?, num(audit.delivered_total_w())?)
        } else { format!("{{\"mode\":\"uniform\",\"source_w_m3\":{}}}", num(uniform_source)?) };
        Ok(format!("{{\"constitutive\":{materials},\"heating\":{heating}}}"))
    }
}

#[cfg(test)]
mod tests;
