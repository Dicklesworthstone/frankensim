//! Native history controls of the actual local P1 radiation operator.
//! Geometry, quadrature, reflected exchange and parameter contractions are
//! owned by fs-conduction. This adapter preserves source identities and units.

use fs_conduction::{
    PointwiseAmbientRadiation, PointwiseEnclosureRadiation, PointwiseEnclosureTolerance,
};
use fs_conduction::radiation::{GrayDiffuseEnclosure, OpenDiffuseEnclosure};
use fs_conduction::transient::enthalpy::adjoint::pointwise as core_pointwise;
use fs_conduction::transient::enthalpy::heterogeneous::EnthalpyRadiation;

use super::{
    ConductionProblem, Cx, EnthalpyStepConfig, EnthalpyStepLinearization,
    EnthalpyStepSolution, HeterogeneousEnthalpyBackwardEuler,
    HeterogeneousEnthalpyStepLinearization, Instant, SolveRefusal, ThermalInterfaces,
    TemperatureContacts,
    bad, failed, finite, json_string, lower, number, poll, radiation, zeros,
};

mod material;
pub(in crate::solve::transient) use material::storage_bytes as material_storage_bytes;

pub(super) fn model_workspace_entries(
    nodes: usize, faces: usize, lowering: &radiation::LoweredRadiation,
) -> Result<usize, SolveRefusal> {
    if lowering.is_material() {
        material::workspace_entries(nodes, faces, lowering)
    } else {
        workspace_entries(nodes, faces, lowering.surface_count())
    }
}

pub(super) fn workspace_entries(
    nodes: usize, faces: usize, patches: usize,
) -> Result<usize, SolveRefusal> {
    core_pointwise::workspace_entries(nodes, faces, patches).map_err(failed)
}

/// One active local/reflected tangent, native control scratch and retained
/// pointwise reports. The common enthalpy allowance separately admits complete
/// checkpoints, conductivity/contact assembly and the thermal Krylov basis.
/// Face geometry is sparse; no node-by-node radiation matrix is admitted.
pub(in crate::solve::transient) fn storage_bytes(
    nodes: usize, faces: usize, patches: usize, fine_steps: usize, open_bath: bool,
) -> Result<usize, SolveRefusal> {
    if fine_steps == 0 {
        return Err(bad("pointwise history controls need a nonempty fine time grid"));
    }
    // Include the completed bath row, its retained factorization and the
    // single shared control in both forward admission and reverse execution.
    let open_bytes = if open_bath {
        patches.checked_add(1).and_then(|p| p.checked_mul(p))
            .and_then(|entries| entries.checked_mul(64))
            .and_then(|bytes| bytes.checked_add(512))
            .ok_or_else(|| bad("open-bath reverse storage allowance overflows"))?
    } else { 0 };
    workspace_entries(nodes, faces, patches)?
        .checked_mul(std::mem::size_of::<f64>())
        .and_then(|bytes| bytes.checked_add(open_bytes))
        .and_then(|bytes| nodes.checked_mul(64).and_then(|extra| bytes.checked_add(extra)))
        .and_then(|bytes| patches.checked_mul(4096).and_then(|extra| bytes.checked_add(extra)))
        .and_then(|bytes| fine_steps.checked_mul(patches)
            .and_then(|rows| rows.checked_mul(2048))
            .and_then(|extra| bytes.checked_add(extra)))
        .filter(|&bytes| isize::try_from(bytes).is_ok())
        .ok_or_else(|| bad("pointwise enthalpy reverse workspace allowance overflows"))
}

/// Bind the exact original convection boundary and unchanged enthalpy field.
/// The core rechecks the accepted endpoint under its original physical target,
/// retains the complete local/reflected tangent, and binds reference density
/// against the same immutable old enthalpy.
#[allow(clippy::too_many_arguments)]
pub(super) fn linearize<'m>(
    cx: &Cx<'_>, engine: &HeterogeneousEnthalpyBackwardEuler<'m, '_>,
    problem: ConductionProblem<'_>, interfaces: Option<&ThermalInterfaces>,
    contacts: Option<&TemperatureContacts<'_, '_>>, max_contact_pairs: usize,
    old_h: &[f64], dt_s: f64, config: EnthalpyStepConfig,
    lowering: &radiation::LoweredRadiation, enclosure: Option<&GrayDiffuseEnclosure>,
    open_enclosure: Option<&OpenDiffuseEnclosure>,
    accepted: EnthalpyStepSolution, physical_target_j: f64,
    max_radiation_entries: usize, deadline: Option<(Instant, f64)>,
) -> Result<HeterogeneousEnthalpyStepLinearization<'m>, SolveRefusal> {
    poll(cx, deadline)?;
    if contacts.is_some() && interfaces.is_some()
        || !lowering.is_pointwise()
        || lowering.is_enclosure() != (enclosure.is_some() || open_enclosure.is_some())
        || enclosure.is_some() && open_enclosure.is_some()
        || lowering.open_ambient_temperature_k().is_some() != open_enclosure.is_some()
    {
        return Err(bad("pointwise history linearization lost its declared spatial model"));
    }
    let tolerance = PointwiseEnclosureTolerance {
        absolute_heat_w: lowering.config.balance_tolerance_w,
        relative_heat: lowering.config.balance_relative_tolerance,
    };
    let result = if lowering.is_material_ambient() {
        if enclosure.is_some() || open_enclosure.is_some() {
            return Err(bad("local-temperature emissivity requires its original ambient owner"));
        }
        let owner = lowering.bind_material_ambient(cx, problem, engine).map_err(lower)?
            .ok_or_else(|| bad("material history linearization lost its original source owner"))?;
        if let Some(contacts) = contacts {
            engine.linearize_temperature_contact_endpoint_with_radiation(
                cx, problem, contacts, old_h, dt_s, config,
                EnthalpyRadiation::MaterialAmbient(&owner), accepted,
                physical_target_j, max_radiation_entries, max_contact_pairs,
            )
        } else {
            engine.linearize_accepted_with_material_ambient_radiation(
                cx, problem, interfaces, old_h, dt_s, config, &owner, accepted,
                physical_target_j, max_radiation_entries,
            )
        }
    } else if lowering.is_material_enclosure() {
        // The material owner borrows the exact prepared geometry and factors.
        // Keep that base alive until the complete endpoint has been retained;
        // the returned linearization borrows only the original thermal mesh.
        let base = if let Some(model) = open_enclosure {
            PointwiseEnclosureRadiation::new_open(
                cx, problem.mesh, problem.boundary, model,
                lowering.open_ambient_temperature_k()
                    .ok_or_else(|| bad("material enclosure history lost its black environment"))?,
                tolerance,
            )
        } else {
            PointwiseEnclosureRadiation::new(
                cx, problem.mesh, problem.boundary,
                enclosure.ok_or_else(|| bad("material enclosure history lost its closed view factors"))?,
                tolerance,
            )
        }.map_err(lower)?;
        let owner = lowering.bind_material_enclosure(cx, problem, engine, &base).map_err(lower)?
            .ok_or_else(|| bad("material enclosure history lost its original source owner"))?;
        if let Some(contacts) = contacts {
            engine.linearize_temperature_contact_endpoint_with_radiation(
                cx, problem, contacts, old_h, dt_s, config,
                EnthalpyRadiation::MaterialEnclosure(&owner), accepted,
                physical_target_j, max_radiation_entries, max_contact_pairs,
            )
        } else {
            engine.linearize_accepted_with_material_enclosure_radiation(
                cx, problem, interfaces, old_h, dt_s, config, &owner, accepted,
                physical_target_j, max_radiation_entries,
            )
        }
    } else if let Some(model) = open_enclosure {
        let owner = PointwiseEnclosureRadiation::new_open(
            cx, problem.mesh, problem.boundary, model,
            lowering.open_ambient_temperature_k()
                .ok_or_else(|| bad("open pointwise history lost its black environment"))?,
            tolerance,
        ).map_err(lower)?;
        if let Some(contacts) = contacts {
            engine.linearize_temperature_contact_endpoint_with_radiation(
                cx, problem, contacts, old_h, dt_s, config,
                EnthalpyRadiation::PointwiseEnclosure(&owner), accepted,
                physical_target_j, max_radiation_entries, max_contact_pairs,
            )
        } else {
            engine.linearize_accepted_with_pointwise_enclosure_radiation(
                cx, problem, interfaces, old_h, dt_s, config, &owner, accepted,
                physical_target_j, max_radiation_entries,
            )
        }
    } else if let Some(model) = enclosure {
        let owner = PointwiseEnclosureRadiation::new(
            cx, problem.mesh, problem.boundary, model, tolerance,
        ).map_err(lower)?;
        if let Some(contacts) = contacts {
            engine.linearize_temperature_contact_endpoint_with_radiation(
                cx, problem, contacts, old_h, dt_s, config,
                EnthalpyRadiation::PointwiseEnclosure(&owner), accepted,
                physical_target_j, max_radiation_entries, max_contact_pairs,
            )
        } else {
            engine.linearize_accepted_with_pointwise_enclosure_radiation(
                cx, problem, interfaces, old_h, dt_s, config, &owner, accepted,
                physical_target_j, max_radiation_entries,
            )
        }
    } else {
        let owner = PointwiseAmbientRadiation::new(
            cx, problem.mesh, problem.boundary, &lowering.patches,
        ).map_err(lower)?;
        if let Some(contacts) = contacts {
            engine.linearize_temperature_contact_endpoint_with_radiation(
                cx, problem, contacts, old_h, dt_s, config,
                EnthalpyRadiation::PointwiseAmbient(&owner), accepted,
                physical_target_j, max_radiation_entries, max_contact_pairs,
            )
        } else {
            engine.linearize_accepted_with_pointwise_ambient_radiation(
                cx, problem, interfaces, old_h, dt_s, config, &owner, accepted,
                physical_target_j, max_radiation_entries,
            )
        }
    }.map_err(failed)?;
    poll(cx, deadline)?;
    Ok(result)
}

struct Surface {
    target: String,
    source_fields: String,
    reservoir_k: Option<f64>,
    emissivity_bar: f64,
    reservoir_bar: f64,
}

struct FixedControls {
    surfaces: Vec<Surface>,
    view_factor_identity: Option<String>,
    /// The open bath is one shared physical control, not one per surface.
    black_environment: Option<BlackEnvironment>,
    steps: usize,
}

struct BlackEnvironment {
    temperature_k: f64,
    temperature_bar: f64,
}

enum ControlModel {
    Fixed(FixedControls),
    Material(material::Controls),
}

pub(super) struct Controls(ControlModel);

impl Controls {
    pub(super) fn new(
        cx: &Cx<'_>, lowering: &radiation::LoweredRadiation,
        enclosure: Option<&GrayDiffuseEnclosure>, open_enclosure: Option<&OpenDiffuseEnclosure>,
        deadline: Option<(Instant, f64)>,
    ) -> Result<Self, SolveRefusal> {
        poll(cx, deadline)?;
        if lowering.is_material() {
            Ok(Self(ControlModel::Material(material::Controls::new(
                cx, lowering, enclosure, open_enclosure, deadline,
            )?)))
        } else {
            Ok(Self(ControlModel::Fixed(FixedControls::new(lowering, enclosure, open_enclosure)?)))
        }
    }

    pub(super) fn add(
        &mut self, cx: &Cx<'_>, step: &EnthalpyStepLinearization<'_>,
        nodal_load_adjoint: &[f64], deadline: Option<(Instant, f64)>,
    ) -> Result<(), SolveRefusal> {
        match &mut self.0 {
            ControlModel::Fixed(controls) => controls.add(cx, step, nodal_load_adjoint, deadline),
            ControlModel::Material(controls) => controls.add(cx, step, nodal_load_adjoint, deadline),
        }
    }

    pub(super) fn fragment(
        &self, cx: &Cx<'_>, expected_steps: usize, max_radiation_entries: usize,
        workspace_bytes: usize, deadline: Option<(Instant, f64)>,
    ) -> Result<String, SolveRefusal> {
        poll(cx, deadline)?;
        match &self.0 {
            ControlModel::Fixed(controls) => controls.fragment(expected_steps, max_radiation_entries, workspace_bytes),
            ControlModel::Material(controls) => controls.fragment(cx, expected_steps, max_radiation_entries, workspace_bytes, deadline),
        }
    }
}

impl FixedControls {
    /// Keep independent authored surfaces even when they share a source card.
    /// An initial selected peak retains these identities with zero derivatives.
    pub(super) fn new(
        lowering: &radiation::LoweredRadiation, enclosure: Option<&GrayDiffuseEnclosure>,
        open_enclosure: Option<&OpenDiffuseEnclosure>,
    ) -> Result<Self, SolveRefusal> {
        let count = lowering.surface_count();
        let enclosed = enclosure.is_some() || open_enclosure.is_some();
        if !lowering.is_pointwise() || lowering.is_enclosure() != enclosed
            || enclosure.is_some() && open_enclosure.is_some()
            || lowering.open_ambient_temperature_k().is_some() != open_enclosure.is_some()
            || !(1..=64).contains(&count)
            || enclosure.is_some_and(|model| model.surfaces().len() != count)
            || open_enclosure.is_some_and(|model| model.surfaces().len() != count)
            || !enclosed && lowering.patches.len() != count
        {
            return Err(bad("pointwise controls lost their complete declared surface layout"));
        }
        let mut surfaces = Vec::with_capacity(count);
        for (index, (name, target)) in lowering.surface_labels().enumerate() {
            let physical_surfaces = enclosure.map(GrayDiffuseEnclosure::surfaces)
                .or_else(|| open_enclosure.map(OpenDiffuseEnclosure::surfaces));
            let (emissivity, reservoir_k) = if let Some(model_surfaces) = physical_surfaces {
                let surface = &model_surfaces[index];
                if target != surface.name() {
                    return Err(bad("pointwise enclosure controls differ from the view-factor order"));
                }
                (surface.emissivity(), None)
            } else {
                let patch = &lowering.patches[index];
                if target != patch.region() {
                    return Err(bad("pointwise ambient controls differ from their declared trace order"));
                }
                (patch.emissivity(), Some(patch.ambient_temperature_k()))
            };
            let directions = if emissivity.value() == 1.0 {
                "nonpositive-at-emissivity-one"
            } else {
                "two-sided-within-zero-exclusive-one-inclusive"
            };
            let source_fields = format!(
                "\"name\":{},\"target\":{},\"card\":{},\"material_state\":{},\"emissivity_receipt\":{},\"emissivity\":{},\"query_temperature_k\":{},\"emissivity_feasible_directions\":{}",
                json_string(name), json_string(target),
                json_string(&emissivity.card_identity().to_hex()),
                json_string(emissivity.material_state()),
                json_string(&emissivity.receipt().content_hash().to_hex()),
                number(emissivity.value())?, number(emissivity.temperature_k())?,
                json_string(directions),
            );
            surfaces.push(Surface {
                target: target.to_owned(), source_fields, reservoir_k,
                emissivity_bar: 0.0, reservoir_bar: 0.0,
            });
        }
        Ok(Self {
            surfaces,
            view_factor_identity: enclosure.map(|model| model.view_factors().identity().to_hex())
                .or_else(|| open_enclosure.map(|model| model.view_factors().identity().to_hex())),
            black_environment: lowering.open_ambient_temperature_k().map(|temperature_k|
                BlackEnvironment { temperature_k, temperature_bar: 0.0 }),
            steps: 0,
        })
    }

    pub(super) fn add(
        &mut self, cx: &Cx<'_>, step: &EnthalpyStepLinearization<'_>,
        nodal_load_adjoint: &[f64], deadline: Option<(Instant, f64)>,
    ) -> Result<(), SolveRefusal> {
        poll(cx, deadline)?;
        // The transport owner supplies P_free * dt * lambda. Outward heat
        // enters the physical residual positively, so its objective bar has
        // the opposite sign. Neither parameter owner introduces another dt.
        let mut outward_weights = zeros(nodal_load_adjoint.len())?;
        for (vertex, (&load, weight)) in nodal_load_adjoint.iter()
            .zip(&mut outward_weights).enumerate()
        {
            if vertex % 256 == 0 { poll(cx, deadline)?; }
            *weight = finite(-load)?;
        }
        let gradient = step.pointwise_parameter_pullback(cx, &outward_weights).map_err(failed)?;
        let count = self.surfaces.len();
        let enclosure = self.view_factor_identity.is_some();
        let reservoirs = if self.black_environment.is_some() { 1 }
            else if enclosure { 0 } else { count };
        if gradient.radiation_regions.len() != count || gradient.emissivities.len() != count
            || gradient.reservoir_temperatures.len() != reservoirs
        {
            return Err(bad("pointwise parameter pullback lost its complete physical control layout"));
        }
        for (index, surface) in self.surfaces.iter_mut().enumerate() {
            poll(cx, deadline)?;
            if gradient.radiation_regions[index] != surface.target {
                return Err(bad("pointwise parameter pullback changed its declared surface order"));
            }
            surface.emissivity_bar = finite(
                surface.emissivity_bar + finite(gradient.emissivities[index])?,
            )?;
            if !enclosure {
                surface.reservoir_bar = finite(
                    surface.reservoir_bar + finite(gradient.reservoir_temperatures[index])?,
                )?;
            }
        }
        if let Some(bath) = &mut self.black_environment {
            bath.temperature_bar = finite(
                bath.temperature_bar + finite(gradient.reservoir_temperatures[0])?,
            )?;
        }
        self.steps = self.steps.checked_add(1)
            .ok_or_else(|| bad("pointwise parameter history count overflows"))?;
        poll(cx, deadline)
    }

    pub(super) fn fragment(
        &self, expected_steps: usize, max_radiation_entries: usize, workspace_bytes: usize,
    ) -> Result<String, SolveRefusal> {
        if self.steps != expected_steps {
            return Err(bad("pointwise controls do not cover the selected reverse prefix"));
        }
        let rows = self.surfaces.iter().map(|surface| {
            let reservoir = surface.reservoir_k.map(|temperature| -> Result<String, SolveRefusal> {
                Ok(format!(
                    ",\"reservoir_temperature_k\":{},\"dtemperature_dreservoir_temperature_k_per_k\":{}",
                    number(temperature)?, number(surface.reservoir_bar)?,
                ))
            }).transpose()?.unwrap_or_default();
            Ok(format!(
                "{{{},\"dtemperature_demissivity_k\":{}{reservoir}}}",
                surface.source_fields, number(surface.emissivity_bar)?,
            ))
        }).collect::<Result<Vec<_>, SolveRefusal>>()?;
        let (method, parameters) = if self.black_environment.is_some() {
            ("complete-pointwise-p1-open-gray-diffuse-enthalpy-tangent",
             "absolute-per-surface-emissivity-and-shared-black-environment-temperature-throughout-selected-trajectory")
        } else if self.view_factor_identity.is_some() {
            ("complete-pointwise-p1-closed-gray-diffuse-enthalpy-tangent",
             "absolute-per-surface-emissivity-throughout-selected-trajectory")
        } else {
            ("complete-pointwise-p1-gray-enthalpy-tangent",
             "absolute-per-surface-emissivity-and-reservoir-temperature-throughout-selected-trajectory")
        };
        let view_factors = self.view_factor_identity.as_ref().map(|identity| format!(
            ",\"view_factor_identity\":{}", json_string(identity),
        )).unwrap_or_default();
        let black_environment = self.black_environment.as_ref().map(|bath| -> Result<String, SolveRefusal> {
            Ok(format!(
                ",\"black_environment\":{{\"temperature_k\":{},\"dtemperature_dambient_temperature_k_per_k\":{},\"feasible_directions\":\"two-sided-within-positive-kelvin\"}}",
                number(bath.temperature_k)?, number(bath.temperature_bar)?,
            ))
        }).transpose()?.unwrap_or_default();
        Ok(format!(
            ",\"radiation\":{{\"method\":\"{method}\",\"spatial_model\":\"pointwise\",\"patch_count\":{},\"control_steps\":{},\"max_radiation_workspace_entries\":{max_radiation_entries},\"additional_workspace_bytes\":{workspace_bytes},\"parameters\":\"{parameters}\",\"parameter_pullback\":\"prepared-pointwise-outward-load-contraction\",\"parameter_krylov_iterations\":0,\"parameter_linear_solves\":0,\"residual_target\":\"original-accepted-forward-physical-tolerance\"{view_factors}{black_environment},\"surfaces\":[{}],\"scope\":{}}}",
            self.surfaces.len(), self.steps, rows.join(","), json_string(SCOPE),
        ))
    }
}

/// Keep the common history/control contract while describing the selected
/// spatial owner accurately; historical patch-mean reports remain unchanged.
pub(super) fn history_scope(common: &str, local_material: bool, material_enclosure: bool) -> String {
    if local_material {
        // The common text describes frozen patch controls. Replace that entire
        // radiation paragraph so a source-ordinate report cannot inherit an
        // independent-per-patch or fixed-emissivity claim. These anchors belong
        // to the caller's internal SCOPE constant, never user-provided text.
        let start = common.find("Declared closed enclosures retain")
            .expect("common enthalpy scope retains its radiation paragraph");
        let end = start + common[start..].find("An isothermal plateau has zero direct temperature sensitivity")
            .expect("common enthalpy scope retains its energetic-history paragraph");
        let scope = if material_enclosure { material::ENCLOSURE_SCOPE } else { material::SCOPE };
        return format!("{}{} {}", &common[..start], scope, &common[end..]);
    }
    common.replace(
        "Their absolute per-surface emissivity controls use the shared radiosity transpose with the same masked time-scaled thermal multiplier; view factors and geometry stay fixed.",
        "Their absolute per-surface emissivity controls contract the prepared reflected response and the local P1 emission with the same masked time-scaled thermal multiplier; view factors and geometry stay fixed.",
    ).replace(
        "Declared ambient radiation uses its complete nonlinear area-mean state feedback in the same h/source/reference-density transpose; the report also differentiates independent absolute surface emissivities and reservoir temperatures held throughout the selected trajectory.",
        "Declared pointwise radiation uses the exact local P1 fourth-power emission tangent and, for closed or open enclosures, the complete nonsymmetric reflected irradiation tangent in the same h/source/reference-density transpose. Independent absolute surface emissivities and ambient reservoir temperatures are held throughout the selected trajectory. Open enclosures expose one shared black-environment temperature control whose response includes every direct and reflected escape path; it is separate from the air and prescribed-convection reservoirs.",
    )
}

const SCOPE: &str = "The complete local P1 fourth-power emission and reflected irradiation enter the same accepted enthalpy transpose as storage, conductivity, contact and physical convection feedback. Parameter partials contract the prepared pointwise outward load with minus the masked time-scaled thermal multiplier. This adds no parameter Krylov or linear solve; construction and physical audits still perform their declared quadrature and prepared radiosity work. Open exchange has one shared black-environment temperature control through the actual reflected bath response; its heat never becomes a solid degree of freedom or an air-convection control. Emissivity is frozen at its declared source query, every radiating node retains the claim's temperature support, and query temperature, geometry, view factors and time grid stay fixed. Shared material cards do not merge distinct surface controls. An initial selected peak has zero future radiation sensitivity. These are discrete branch derivatives, not temperature-dependent emissivity, continuum bounds or physical validation.";

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn pointwise_workspace_admits_independent_surface_layout_and_checkpoint_reports() {
        let bare = workspace_entries(12, 16, 2).unwrap() * std::mem::size_of::<f64>();
        let one = storage_bytes(12, 16, 2, 1, false).unwrap();
        let several = storage_bytes(12, 16, 2, 8, false).unwrap();
        assert!(one > bare && several > one);
        assert!(storage_bytes(12, 16, 2, 8, true).unwrap() > several);
        assert!(storage_bytes(12, 16, 0, 1, false).is_err());
        assert!(storage_bytes(12, 16, 65, 1, true).is_err());
        assert!(storage_bytes(12, 16, 2, 0, false).is_err());
        assert!(storage_bytes(12, 16, 2, usize::MAX, false).is_err());
    }

    #[test]
    fn pointwise_scope_preserves_history_controls_without_a_patch_mean_claim() {
        let scope = history_scope(super::super::SCOPE, false, false);
        assert!(scope.contains("local P1 fourth-power emission tangent"));
        assert!(scope.contains("complete nonsymmetric reflected irradiation tangent"));
        assert!(scope.contains("Initial controls are additive nodal or regional specific-enthalpy offsets"));
        assert!(!scope.contains("nonlinear area-mean state feedback"));
        assert!(!scope.contains("controls use the shared radiosity transpose"));
    }

    #[test]
    fn material_scope_retains_history_and_names_original_shared_source_controls() {
        let scope = history_scope(super::super::SCOPE, true, false);
        assert!(scope.contains("original pinned scalar or linear T-curve emissivity"));
        assert!(scope.contains("exact material-card and selected-claim identity"));
        assert!(scope.contains("Initial controls are additive nodal or regional specific-enthalpy offsets"));
        assert!(scope.contains("h history carry remains energetic"));
        assert!(!scope.contains("sharing a card does not merge controls"));
        assert!(!scope.contains("independent absolute surface emissivities"));
        assert!(!scope.contains("controls use the shared radiosity transpose"));
    }

    #[test]
    fn material_enclosure_scope_keeps_reflected_source_controls_and_one_open_bath() {
        let scope = history_scope(super::super::SCOPE, true, true);
        assert!(scope.contains("current reflected radiosity"));
        assert!(scope.contains("one shared open black-bath temperature"));
        assert!(scope.contains("retained radiosity transpose backsolve per consumed endpoint"));
        assert!(scope.contains("h history carry remains energetic"));
        assert!(!scope.contains("independent absolute surface emissivities"));
        assert!(!scope.contains("Per-surface black-bath temperatures remain independent"));
    }
}

