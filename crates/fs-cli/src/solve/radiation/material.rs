//! Explicit local-temperature emissivity on the original enthalpy trace.
//! Lowering retains the selected source; binding borrows the actual storage
//! charts and physical boundary instead of reconstructing caloric objects.

use fs_conduction::{MaterialAmbientRadiationPatch, MaterialPointwiseAmbientRadiation};
use fs_conduction::transient::enthalpy::heterogeneous::HeterogeneousEnthalpyBackwardEuler;
use fs_material::phase::EquilibriumEnthalpyPhaseCurve;
use fs_material::state_point::MaterialPropertySelection;
use fs_material::thermal::source::ThermalMaterialPoint;
use fs_matdb::{InterpolationPolicy, MaterialCard, PropertyClaim, PropertyValue, QueryPoint};
use fs_project::{RadiatingSurface, RadiationSpatialModel, TEMPERATURE_AXIS};

use super::*;

/// One immutable source per authored surface, never one cloned card per node.
pub(super) struct MaterialAmbientSource {
    material: MaterialCard,
    context: QueryPoint,
    selection: MaterialPropertySelection,
}

impl MaterialAmbientSource {
    fn point(&self) -> ThermalMaterialPoint<'_> {
        ThermalMaterialPoint {
            material: &self.material,
            context: &self.context,
            temperature_axis: TEMPERATURE_AXIS,
            liquid_fraction_axis: None,
        }
    }
}

fn bad(what: impl Into<String>) -> SolveRefusal {
    conduction_error(
        "cli-solve-conduction-material-emissivity",
        what,
        "declare pointwise radiation with primary enthalpy and a static scalar or linear T emissivity source on each surface's original storage-material card",
    )
}

fn physical(surface: &str, what: impl Into<String>) -> ConductionError {
    ConductionError::Radiation {
        surface: surface.into(),
        what: what.into(),
        fix: "each local-temperature surface must use one original storage material and its exact emissivity card; split traces at material boundaries".into(),
    }
}

fn poll(cx: &Cx<'_>, at: usize) -> Result<(), ConductionError> {
    cx.checkpoint().map_err(|_| ConductionError::Cancelled {
        stage: "native-material-ambient-binding", at,
    })
}

/// Refuse incompatible source intent before a fixed-emissivity path can run.
pub(super) fn admit(setup: &ConductionSetup) -> Result<(), SolveRefusal> {
    if let Some(radiation) = &setup.radiation
        && radiation.emissivity_response.is_some()
        && (radiation.spatial_model != Some(RadiationSpatialModel::Pointwise)
            || !setup.transient.as_ref().is_some_and(|time|
                time.enthalpy.is_some() && time.capacities.is_empty()))
    {
        return Err(bad("local-temperature emissivity requires pointwise radiation and primary total-enthalpy storage"));
    }
    Ok(())
}

pub(super) fn lower_source(
    surface: &RadiatingSurface,
    cards: &CardPackSet,
) -> Result<(SurfaceEmissivity, MaterialAmbientSource), SolveRefusal> {
    let library = cards.library();
    let card = library.material(&surface.card).ok_or_else(||
        bad(format!("surface `{}` references absent emissivity card {}", surface.name, surface.card)))?;
    let context = QueryPoint::new().with(TEMPERATURE_AXIS, surface.query_temperature.value)
        .map_err(|error| bad(format!("surface `{}` source query refused: {error}", surface.name)))?;
    let selection = match surface.claim.as_deref() {
        Some(pin) => MaterialPropertySelection::PinnedByProperty(vec![(
            fs_conduction::SURFACE_EMISSIVITY_PROPERTY.into(),
            parse_claim_id(pin, &format!("radiating surface {}", surface.name))?,
        )]),
        None => MaterialPropertySelection::SingleClaimOnly,
    };
    let (emissivity, _) = SurfaceEmissivity::from_condition(
        &surface.name,
        ThermalMaterialPoint {
            material: card, context: &context,
            temperature_axis: TEMPERATURE_AXIS, liquid_fraction_axis: None,
        },
        selection,
    ).map_err(|error| bad(format!("surface `{}` emissivity source refused: {error}", surface.name)))?;
    let selected = emissivity.receipt().selected;
    let claim = card.claims().claim(selected)
        .ok_or_else(|| bad("selected emissivity claim is absent from its original card"))?;
    match (&claim.value, &claim.interpolation) {
        (PropertyValue::Scalar { .. }, InterpolationPolicy::ConstantWithinValidity) => {}
        (PropertyValue::Curve { abscissa, .. }, InterpolationPolicy::LinearInside)
            if abscissa == TEMPERATURE_AXIS => {}
        _ => return Err(bad(format!(
            "surface `{}` local-temperature emissivity requires a static scalar or original linear T curve",
            surface.name,
        ))),
    }
    Ok((emissivity, MaterialAmbientSource {
        material: card.clone(), context,
        selection: MaterialPropertySelection::PinnedByProperty(vec![(
            fs_conduction::SURFACE_EMISSIVITY_PROPERTY.into(), selected,
        )]),
    }))
}

impl LoweredRadiation {
    /// Original source controls, including zero-step/initial-peak adjoints.
    /// Fixed snapshots and absent patch indices have no evolving source law.
    pub(in crate::solve) fn material_ambient_claim(&self, patch: usize) -> Option<&PropertyClaim> {
        if !self.is_material_ambient() { return None; }
        self.material_sources.get(patch)?.material.claims()
            .claim(self.surfaces.get(patch)?.emissivity.receipt().selected)
    }

    /// Bind one owner to the exact current boundary and storage chart objects.
    /// Reuse this entry point for forward, endpoint-audit and adjoint responses.
    pub(in crate::solve) fn bind_material_ambient<'a>(
        &'a self,
        cx: &Cx<'_>,
        problem: ConductionProblem<'a>,
        engine: &'a HeterogeneousEnthalpyBackwardEuler<'_, '_>,
    ) -> Result<Option<MaterialPointwiseAmbientRadiation<'a>>, ConductionError> {
        if !self.is_material_ambient() { return Ok(None); }
        poll(cx, 0)?;
        if !self.is_pointwise() || self.is_enclosure()
            || self.material_sources.len() != self.surfaces.len()
            || self.patches.len() != self.surfaces.len()
        {
            return Err(physical("material-ambient", "incomplete local-temperature ambient lowering"));
        }
        if !std::ptr::eq(problem.mesh, engine.mesh()) {
            return Err(physical("material-ambient", "emissivity and enthalpy require the same original mesh"));
        }
        let mut regions = vec![None; problem.boundary.region_names().len()];
        let mut identities = Vec::with_capacity(self.surfaces.len());
        let mut owners: Vec<Option<(usize, &EquilibriumEnthalpyPhaseCurve)>> =
            vec![None; self.surfaces.len()];
        for (index, (surface, source)) in self.surfaces.iter().zip(&self.material_sources).enumerate() {
            poll(cx, index)?;
            let region = problem.boundary.region_names().iter()
                .position(|name| name == &surface.target)
                .ok_or_else(|| physical(&surface.name, "local-temperature target has no physical boundary"))?;
            if regions[region].replace(index).is_some() {
                return Err(physical(&surface.name, "duplicate local-temperature radiation target"));
            }
            identities.push(source.material.content_hash());
        }
        for (slot, face) in problem.mesh.boundary().iter().enumerate() {
            if slot % 256 == 0 { poll(cx, slot)?; }
            let Some(index) = problem.boundary.region_for(slot).and_then(|region| regions[region])
                else { continue; };
            let surface = &self.surfaces[index];
            for vertex in face.vertices {
                let vertex = vertex as usize;
                let material_id = *engine.vertex_material_ids().get(vertex)
                    .ok_or_else(|| physical(&surface.name, "surface vertex has no storage-material owner"))?;
                let curve = engine.phase_curve_at(vertex)
                    .ok_or_else(|| physical(&surface.name, "surface vertex has no original enthalpy chart"))?;
                if curve.material_card_identity() != identities[index] {
                    return Err(physical(&surface.name,
                        "surface emissivity and enthalpy require the same original material card"));
                }
                if let Some((expected_id, expected_curve)) = owners[index] {
                    if material_id != expected_id || !std::ptr::eq(curve, expected_curve) {
                        return Err(physical(&surface.name,
                            "one local-temperature surface crosses distinct storage-material owners"));
                    }
                } else {
                    owners[index] = Some((material_id, curve));
                }
            }
        }
        let mut patches = Vec::with_capacity(self.surfaces.len());
        for (index, (surface, source)) in self.surfaces.iter().zip(&self.material_sources).enumerate() {
            poll(cx, index)?;
            let (_, curve) = owners[index].ok_or_else(||
                physical(&surface.name, "local-temperature target has no exterior area"))?;
            patches.push(MaterialAmbientRadiationPatch {
                region: &surface.target,
                ambient_temperature_k: self.patches[index].ambient_temperature_k(),
                curve,
                material: source.point(),
                selection: source.selection.clone(),
            });
        }
        MaterialPointwiseAmbientRadiation::new(cx, problem.mesh, problem.boundary, &patches).map(Some)
    }
}
