//! Original source coefficients change without resetting the accepted solid.
//!
//! The prescribed state change is an external parameter operation at fixed
//! geometry, mass, displacement and physical Maxwell spring strains. Its signed
//! storage change is external work, not dissipation or heat. This deliberately
//! does not claim thermal expansion, entropy closure or a solved heat feedback.

use super::*;
use fs_matdb::{MaterialCard, QueryPoint};
use fs_material::visco_tensor::source::ResolvedTensorMaterial;
use fs_solid::linear3::visco::{SolidProjectionBudget, TetViscoMaterial};

/// One immutable material source and its explicit reference orientation.
#[derive(Clone, Debug)]
pub struct SolidMaterialSource {
    pub material: ResolvedTensorMaterial,
    pub material_to_mesh: [[f64; 3]; 3],
}

/// Complete proposed continuation. The original accepted history is untouched
/// on success as well as on refusal; publishing this value is the caller's act.
pub struct SolidMaterialChange<'a, 's> {
    pub history: TensorSolidHistory<'a, 's>,
    pub parameter_work_j: f64,
}

impl<'s> TensorSolidBody<'s> {
    /// Bind exact source data to every original mesh element. Original density,
    /// law and source identity must match the material that actually assembled
    /// the reference mass and quadrature; a caller label is insufficient.
    pub fn from_material_sources(
        projection: ProjectedViscoSolid<'s>,
        sources: Vec<SolidMaterialSource>,
        budget: TensorMaterialBudget,
        gate: &CancelGate,
    ) -> SolidResult<Self> {
        poll(gate)?;
        if sources.len() != projection.element_count() {
            return Err(bad("one original material source per solid element required").into());
        }
        for (element, source) in sources.iter().enumerate() {
            poll(gate)?;
            let supplied =
                TetViscoMaterial::from_source(&source.material, source.material_to_mesh)?;
            let original = projection
                .source()
                .reference_material(element)
                .ok_or_else(|| bad("missing original solid element material"))?;
            if supplied.law != original.law
                || supplied.density_kg_m3.to_bits() != original.density_kg_m3.to_bits()
                || supplied.material_state_identity != original.material_state_identity
            {
                return Err(bad(
                    "solid source must be the exact original material and reference density",
                )
                .into());
            }
            for (point, &owner) in projection
                .material_points()
                .iter()
                .zip(projection.point_elements())
            {
                if owner == element && point.material() != &supplied.law {
                    return Err(bad("projected solid law differs from its original source").into());
                }
            }
        }
        let mut body = Self::new(projection, budget, gate)?;
        body.sources = Some(sources);
        Ok(body)
    }

    /// Current complete source queries in original element order.
    pub fn material_sources(&self) -> Option<&[SolidMaterialSource]> {
        self.sources.as_deref()
    }

    /// Prepare a complete replacement operator from the SAME original cards and
    /// selected property claims. Every physical query is explicit. The original
    /// basis and reference inertia survive; current density must remain equal
    /// because fixed geometry plus a changed density would change specimen mass.
    ///
    /// This evaluates endpoints only. The caller owns the prescribed parameter
    /// path and its applicability; no interpolation across unsupported source
    /// states, thermal strain, phase transition or basis-error claim is inferred.
    pub fn at_query_points(
        &self,
        cards: &[&MaterialCard],
        points: &[QueryPoint],
        projection_budget: SolidProjectionBudget,
        dynamics_budget: TensorMaterialBudget,
        gate: &CancelGate,
    ) -> SolidResult<Self> {
        self.at_query_points_with_thermal_strains(
            cards,
            points,
            None,
            projection_budget,
            dynamics_budget,
            gate,
        )
    }

    pub(crate) fn at_query_points_with_thermal_strains(
        &self,
        cards: &[&MaterialCard],
        points: &[QueryPoint],
        strains: Option<&[fs_solid::linear3::TetThermalStrainState]>,
        projection_budget: SolidProjectionBudget,
        dynamics_budget: TensorMaterialBudget,
        gate: &CancelGate,
    ) -> SolidResult<Self> {
        poll(gate)?;
        if self.thermal_expansion_identity.is_some() && strains.is_none() {
            return Err(bad("source-bound expansion must be refreshed by its original thermal owner").into());
        }
        let current = self
            .sources
            .as_ref()
            .ok_or_else(|| bad("material refresh requires original solid source bindings"))?;
        if cards.len() != current.len() || points.len() != current.len() {
            return Err(
                bad("complete original cards and physical queries per element required").into(),
            );
        }
        let mut sources = Vec::with_capacity(current.len());
        let mut laws = Vec::with_capacity(current.len());
        for ((source, &card), point) in current.iter().zip(cards).zip(points) {
            poll(gate)?;
            let material = source.material.at_query_point(card, point)?;
            if material.density_kg_m3().to_bits() != source.material.density_kg_m3().to_bits() {
                return Err(bad("density-changing source update requires a conservative geometry/mass evolution owner").into());
            }
            laws.push(material.rotated_law(source.material_to_mesh)?);
            sources.push(SolidMaterialSource {
                material,
                material_to_mesh: source.material_to_mesh,
            });
        }
        let mut projection = self
            .projection
            .with_element_laws(&laws, projection_budget, gate)?;
        if let Some(strains) = strains {
            if strains.len() != sources.len()
                || strains.iter().zip(&sources).any(|(strain, source)| {
                    strain.elastic_material_state_identity() != source.material.state().identity().0
                })
            {
                return Err(bad(
                    "thermal eigenstrain must retain the exact rebuilt material identity",
                )
                .into());
            }
            projection = projection.with_element_thermal_strains(strains, gate)?;
        }
        let mut body = Self::new(projection, dynamics_budget, gate)?;
        body.sources = Some(sources);
        body.thermal_expansion_identity = self.thermal_expansion_identity.clone();
        poll(gate)?;
        Ok(body)
    }
}

impl<'s> TensorSolidHistory<'_, 's> {
    /// Continue at the SAME physical instant with rebuilt source coefficients.
    /// Source branch selectors/order and orientation must be unchanged. Physical
    /// Maxwell spring strains survive, so a stiffness change changes their
    /// stress/storage without incorrectly keeping old energy-factor coordinates.
    /// The native exact dynamics then uses new stiffness and relaxation times.
    ///
    /// Displacement, momentum, time, step count and accumulated heat/force work
    /// are retained. Signed coefficient work has a distinct cumulative account.
    /// Changed singular factors at either endpoint, a foreign source/geometry,
    /// insufficient work allowance, material-domain exit or cancellation returns
    /// no continuation.
    pub fn rebind_materials<'b>(
        &self,
        next: &'b TensorSolidBody<'s>,
        budget: MaterialChangeBudget,
        gate: &CancelGate,
    ) -> SolidResult<SolidMaterialChange<'b, 's>> {
        poll(gate)?;
        let previous =
            self.solid.sources.as_ref().ok_or_else(|| {
                bad("live material continuation requires original source bindings")
            })?;
        let sources = next
            .sources
            .as_ref()
            .ok_or_else(|| bad("replacement solid lacks original source bindings"))?;
        if !core::ptr::eq(self.solid.projection.source(), next.projection.source())
            || self.solid.projection.basis() != next.projection.basis()
            || self.solid.projection.point_elements() != next.projection.point_elements()
            || previous.len() != sources.len()
            || self.solid.thermal_expansion_identity.is_some()
                && self.solid.thermal_expansion_identity != next.thermal_expansion_identity
        {
            return Err(bad(
                "material continuation cannot change original geometry, supports or basis",
            )
            .into());
        }
        for (old, new) in previous.iter().zip(sources) {
            poll(gate)?;
            let a = old.material.state();
            let b = new.material.state();
            if old.material_to_mesh != new.material_to_mesh
                || old.material.properties() != new.material.properties()
                || old.material.density_kg_m3().to_bits() != new.material.density_kg_m3().to_bits()
                || a.card_identity() != b.card_identity()
                || a.material() != b.material()
                || a.axis_quantities() != b.axis_quantities()
                || !a
                    .query_point()
                    .iter()
                    .map(|p| &p.0)
                    .eq(b.query_point().iter().map(|p| &p.0))
                || a.properties().len() != b.properties().len()
                || a.properties().iter().zip(b.properties()).any(|(x, y)| {
                    x.requirement() != y.requirement()
                        || x.answer().receipt.selected != y.answer().receipt.selected
                })
            {
                return Err(bad("continuation must preserve original tensor source claims, branch identities and reference frame").into());
            }
        }
        let (inner, parameter_work_j) = self.inner.rebind_body(&next.body, budget, gate)?;
        poll(gate)?;
        Ok(SolidMaterialChange {
            history: TensorSolidHistory { solid: next, inner },
            parameter_work_j,
        })
    }
}
