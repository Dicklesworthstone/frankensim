//! Mesh-bound heterogeneous caloric, conductivity and initial-state inputs.
//!
//! Each material owns its phase chart and frozen reference density. Different
//! charts cannot share an enthalpy vertex: the native heterogeneous preparation
//! enforces this instead of inventing a mixture. Contact laws are separate native
//! inputs; this file does not turn touching surfaces into insulated gaps.
use std::collections::{BTreeMap, BTreeSet};
use crate::specimen::Result;
use fs_blake3::ContentHash;
use fs_conduction::{ConductionMesh, ConductivityModel};
use fs_conduction::material::{ElementMaterials, MaterialId, MaterialTable};
use fs_conduction::transient::enthalpy::EnthalpyBudget;
use fs_conduction::transient::enthalpy::heterogeneous::{
    HeterogeneousEnthalpyBackwardEuler, ReferenceEnthalpyMaterial,
};
use fs_exec::Cx;
use fs_material::phase::{EnthalpyPhaseKnot, EquilibriumEnthalpyPhaseCurve};

/// Total declaration bounds, not per-material allowances or a solver byte cap.
#[derive(Debug, Clone, Copy)]
pub struct ScanMaterialBudget {
    pub max_bytes: usize,
    pub max_materials: usize,
    pub max_phase_knots: usize,
}

/// Owned constitutive inputs bound to the exact original mesh. Getters expose
/// immutable data only, so conductivity, storage and initial-state assignments
/// cannot diverge between candidates or between a fit and its replay.
pub struct ScanMaterialMap<'m> {
    mesh: &'m ConductionMesh,
    conductivity: ElementMaterials,
    curves: Vec<EquilibriumEnthalpyPhaseCurve>,
    densities: Vec<f64>,
    element_ids: Vec<usize>,
    initial_h: Vec<f64>,
    reference_mass_kg: f64,
    identity: ContentHash,
}

impl<'m> ScanMaterialMap<'m> {
    /// FS_SCAN_MATERIALS_V1 rows:
    /// material ID K_W_m_K RHO_kg_m3 INITIAL_H_J_kg
    /// phase_knot ID H_J_kg T_K LIQUID_MASS_FRACTION
    /// element TET_INDEX MATERIAL_ID
    /// initial VERTEX_INDEX H_J_kg (optional explicit override).
    ///
    /// Every tetrahedron must be assigned exactly once. IDs need not be dense;
    /// native storage indices are constructed deterministically in ID order.
    /// All declared materials must be used. Knot order is authored order, never
    /// sorted to repair an invalid phase law. All units are SI and all properties
    /// are caller declarations, not calibrated material-card claims.
    pub fn from_text(cx: &Cx<'_>, mesh: &'m ConductionMesh, text: &str,
        budget: ScanMaterialBudget) -> Result<Self>
    {
        poll(cx)?;
        if text.len() > budget.max_bytes { return Err("material file exceeds its byte allowance".into()); }
        let mut lines = text.lines().enumerate().filter_map(|(i, line)| {
            let line = line.split('#').next().unwrap_or("").trim();
            (!line.is_empty()).then_some((i + 1, line))
        });
        if lines.next().map(|(_, line)| line) != Some("FS_SCAN_MATERIALS_V1") {
            return Err("material file must begin with FS_SCAN_MATERIALS_V1".into());
        }
        let mut materials = BTreeMap::<u32, [f64; 3]>::new();
        let mut knots = BTreeMap::<u32, Vec<[f64; 3]>>::new();
        let mut assignments = vec![None; mesh.element_count()];
        let mut overrides = BTreeMap::<usize, f64>::new();
        let mut knot_count = 0usize;
        for (line_number, line) in lines {
            poll(cx)?;
            let result = (|| -> Result<()> {
                let fields: Vec<_> = line.split_whitespace().take(7).collect();
                match fields.as_slice() {
                    ["material", id, k, rho, initial] => {
                        if materials.len() >= budget.max_materials { return Err("material count allowance exceeded".into()); }
                        let id: u32 = id.parse()?;
                        let values = [number(k)?, number(rho)?, number(initial)?];
                        if values[0] <= 0.0 || values[1] <= 0.0 { return Err("conductivity and reference density must be positive".into()); }
                        if materials.insert(id, values).is_some() { return Err("duplicate material ID".into()); }
                    }
                    ["phase_knot", id, h, t, fraction] => {
                        if knot_count >= budget.max_phase_knots { return Err("total phase-knot allowance exceeded".into()); }
                        let id: u32 = id.parse()?;
                        if !materials.contains_key(&id) { return Err("declare the material before its phase knots".into()); }
                        knots.entry(id).or_default().push([number(h)?, number(t)?, number(fraction)?]);
                        knot_count += 1;
                    }
                    ["element", element, id] => {
                        let element: usize = element.parse()?;
                        let id: u32 = id.parse()?;
                        let slot = assignments.get_mut(element).ok_or("material assignment names an absent tetrahedron")?;
                        if slot.replace(id).is_some() { return Err("duplicate tetrahedron material assignment".into()); }
                    }
                    ["initial", vertex, h] => {
                        let vertex: usize = vertex.parse()?;
                        if vertex >= mesh.vertex_count() { return Err("initial enthalpy names an absent vertex".into()); }
                        if overrides.insert(vertex, number(h)?).is_some() { return Err("duplicate initial-enthalpy override".into()); }
                    }
                    _ => return Err("unknown material row or wrong field count".into()),
                }
                Ok(())
            })();
            result.map_err(|e| format!("material line {line_number}: {e}"))?;
        }
        if materials.is_empty() { return Err("material file needs a nonempty material table".into()); }
        let assigned = assignments.into_iter().collect::<Option<Vec<_>>>()
            .ok_or("every tetrahedron needs one explicit material assignment")?;
        let used: BTreeSet<_> = assigned.iter().copied().collect();
        if used != materials.keys().copied().collect() {
            return Err("material assignments contain an unknown ID or a declared material is unused".into());
        }
        let identity = fs_blake3::hash_domain("scan-fit-material-map-v1", text.as_bytes());
        let mut curves = Vec::with_capacity(materials.len());
        let mut densities = Vec::with_capacity(materials.len());
        let mut initial = Vec::with_capacity(materials.len());
        let mut intervals = Vec::with_capacity(materials.len());
        let mut models = Vec::with_capacity(materials.len());
        let mut native_ids = BTreeMap::new();
        for (&id, &[k, rho, h0]) in &materials {
            poll(cx)?;
            let rows = knots.remove(&id).ok_or("each material requires its own phase knots")?;
            if rows.len() < 2 { return Err("each material requires at least two phase knots".into()); }
            let span = [rows[0][0], rows[rows.len() - 1][0]];
            if h0 < span[0] || h0 > span[1] { return Err("initial enthalpy is outside its material's caloric interval".into()); }
            let mut key = [0u8; 36];
            key[..32].copy_from_slice(&identity.0);
            key[32..].copy_from_slice(&id.to_le_bytes());
            let curve = EquilibriumEnthalpyPhaseCurve::try_new(
                fs_blake3::hash_domain("scan-fit-material-curve-v1", &key),
                rows.into_iter().map(|[h, t, f]| EnthalpyPhaseKnot {
                    specific_enthalpy_j_kg: h, temperature_k: t,
                    liquid_mass_fraction: f, bulk_density_kg_m3: rho,
                }).collect())?;
            native_ids.insert(id, curves.len());
            curves.push(curve);
            densities.push(rho);
            initial.push(h0);
            intervals.push(span);
            models.push((MaterialId(id), ConductivityModel::isotropic_declared(k)?));
        }
        let element_ids: Vec<_> = assigned.iter().map(|id| native_ids[id]).collect();
        let conductivity = ElementMaterials::new(MaterialTable::new(models)?,
            assigned.into_iter().map(MaterialId).collect())?;
        let references: Vec<_> = curves.iter().zip(&densities).map(|(curve, &rho)|
            ReferenceEnthalpyMaterial { curve, reference_density_kg_m3: rho }).collect();
        // This is native storage admission, not another thermal advance. It
        // checks chart ownership and supplies the same lumped reference masses.
        let prepared = HeterogeneousEnthalpyBackwardEuler::new(cx, mesh, &references,
            &element_ids, EnthalpyBudget {
                max_vertices: mesh.vertex_count(), max_elements: mesh.element_count(),
            })?;
        let mut initial_h = Vec::with_capacity(mesh.vertex_count());
        let mut reference_mass_kg = 0.0;
        for (vertex, (&id, &mass)) in prepared.vertex_material_ids().iter()
            .zip(prepared.reference_nodal_masses_kg()).enumerate()
        {
            if vertex % 256 == 0 { poll(cx)?; }
            let h = overrides.get(&vertex).copied().unwrap_or(initial[id]);
            if h < intervals[id][0] || h > intervals[id][1] {
                return Err(format!("initial enthalpy at vertex {vertex} is outside its assigned material chart").into());
            }
            initial_h.push(h);
            reference_mass_kg += mass;
        }
        if !reference_mass_kg.is_finite() || reference_mass_kg <= 0.0 {
            return Err("total reference mass is not representable".into());
        }
        drop(prepared);
        drop(references);
        poll(cx)?;
        Ok(Self { mesh, conductivity, curves, densities, element_ids,
            initial_h, reference_mass_kg, identity })
    }

    pub fn mesh(&self) -> &'m ConductionMesh { self.mesh }
    pub fn conductivity(&self) -> &ElementMaterials { &self.conductivity }
    pub fn element_material_ids(&self) -> &[usize] { &self.element_ids }
    pub fn initial_specific_enthalpy_j_kg(&self) -> &[f64] { &self.initial_h }
    pub fn reference_mass_kg(&self) -> f64 { self.reference_mass_kg }
    pub fn material_count(&self) -> usize { self.curves.len() }
    pub fn identity(&self) -> ContentHash { self.identity }

    /// Borrow each original chart; no material selection or curve averaging.
    pub fn references(&self) -> Vec<ReferenceEnthalpyMaterial<'_>> {
        self.curves.iter().zip(&self.densities).map(|(curve, &rho)|
            ReferenceEnthalpyMaterial { curve, reference_density_kg_m3: rho }).collect()
    }
}

fn number(text: &str) -> Result<f64> {
    let value: f64 = text.parse()?;
    if !value.is_finite() { return Err("material values must be finite SI numbers".into()); }
    Ok(value)
}
fn poll(cx: &Cx<'_>) -> Result<()> {
    cx.checkpoint().map_err(|_| "material binding cancelled".into())
}

#[cfg(test)]
mod tests;
