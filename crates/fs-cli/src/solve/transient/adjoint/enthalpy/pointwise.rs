//! Native history controls of the actual local P1 radiation operator.
//! Geometry, quadrature, reflected exchange and parameter contractions are
//! owned by fs-conduction. This adapter preserves source identities and units.

use fs_conduction::{
    PointwiseAmbientRadiation, PointwiseEnclosureRadiation, PointwiseEnclosureTolerance,
};
use fs_conduction::radiation::GrayDiffuseEnclosure;
use fs_conduction::transient::enthalpy::adjoint::pointwise as core_pointwise;

use super::{
    ConductionProblem, Cx, EnthalpyStepConfig, EnthalpyStepLinearization,
    EnthalpyStepSolution, HeterogeneousEnthalpyBackwardEuler,
    HeterogeneousEnthalpyStepLinearization, Instant, SolveRefusal, ThermalInterfaces,
    bad, failed, finite, json_string, lower, number, poll, radiation, zeros,
};

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
    nodes: usize, faces: usize, patches: usize, fine_steps: usize,
) -> Result<usize, SolveRefusal> {
    if fine_steps == 0 {
        return Err(bad("pointwise history controls need a nonempty fine time grid"));
    }
    workspace_entries(nodes, faces, patches)?
        .checked_mul(std::mem::size_of::<f64>())
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
    old_h: &[f64], dt_s: f64, config: EnthalpyStepConfig,
    lowering: &radiation::LoweredRadiation, enclosure: Option<&GrayDiffuseEnclosure>,
    accepted: EnthalpyStepSolution, physical_target_j: f64,
    max_radiation_entries: usize, deadline: Option<(Instant, f64)>,
) -> Result<HeterogeneousEnthalpyStepLinearization<'m>, SolveRefusal> {
    poll(cx, deadline)?;
    if !lowering.is_pointwise() || lowering.is_enclosure() != enclosure.is_some() {
        return Err(bad("pointwise history linearization lost its declared spatial model"));
    }
    let result = if let Some(model) = enclosure {
        let owner = PointwiseEnclosureRadiation::new(
            cx, problem.mesh, problem.boundary, model,
            PointwiseEnclosureTolerance {
                absolute_heat_w: lowering.config.balance_tolerance_w,
                relative_heat: lowering.config.balance_relative_tolerance,
            },
        ).map_err(lower)?;
        engine.linearize_accepted_with_pointwise_enclosure_radiation(
            cx, problem, interfaces, old_h, dt_s, config, &owner, accepted,
            physical_target_j, max_radiation_entries,
        )
    } else {
        let owner = PointwiseAmbientRadiation::new(
            cx, problem.mesh, problem.boundary, &lowering.patches,
        ).map_err(lower)?;
        engine.linearize_accepted_with_pointwise_ambient_radiation(
            cx, problem, interfaces, old_h, dt_s, config, &owner, accepted,
            physical_target_j, max_radiation_entries,
        )
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

pub(super) struct Controls {
    surfaces: Vec<Surface>,
    view_factor_identity: Option<String>,
    steps: usize,
}

impl Controls {
    /// Keep independent authored surfaces even when they share a source card.
    /// An initial selected peak retains these identities with zero derivatives.
    pub(super) fn new(
        lowering: &radiation::LoweredRadiation, enclosure: Option<&GrayDiffuseEnclosure>,
    ) -> Result<Self, SolveRefusal> {
        let count = lowering.surface_count();
        if !lowering.is_pointwise() || lowering.is_enclosure() != enclosure.is_some()
            || !(1..=64).contains(&count)
            || enclosure.is_some_and(|model| model.surfaces().len() != count)
            || enclosure.is_none() && lowering.patches.len() != count
        {
            return Err(bad("pointwise controls lost their complete declared surface layout"));
        }
        let mut surfaces = Vec::with_capacity(count);
        for (index, (name, target)) in lowering.surface_labels().enumerate() {
            let (emissivity, reservoir_k) = if let Some(model) = enclosure {
                let surface = &model.surfaces()[index];
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
            view_factor_identity: enclosure.map(|model| model.view_factors().identity().to_hex()),
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
        if gradient.radiation_regions.len() != count || gradient.emissivities.len() != count
            || gradient.reservoir_temperatures.len() != (if enclosure { 0 } else { count })
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
        let (method, parameters) = if self.view_factor_identity.is_some() {
            ("complete-pointwise-p1-closed-gray-diffuse-enthalpy-tangent",
             "absolute-per-surface-emissivity-throughout-selected-trajectory")
        } else {
            ("complete-pointwise-p1-gray-enthalpy-tangent",
             "absolute-per-surface-emissivity-and-reservoir-temperature-throughout-selected-trajectory")
        };
        let view_factors = self.view_factor_identity.as_ref().map(|identity| format!(
            ",\"view_factor_identity\":{}", json_string(identity),
        )).unwrap_or_default();
        Ok(format!(
            ",\"radiation\":{{\"method\":\"{method}\",\"spatial_model\":\"pointwise\",\"patch_count\":{},\"control_steps\":{},\"max_radiation_workspace_entries\":{max_radiation_entries},\"additional_workspace_bytes\":{workspace_bytes},\"parameters\":\"{parameters}\",\"parameter_pullback\":\"prepared-pointwise-outward-load-contraction\",\"parameter_krylov_iterations\":0,\"parameter_linear_solves\":0,\"residual_target\":\"original-accepted-forward-physical-tolerance\"{view_factors},\"surfaces\":[{}],\"scope\":{}}}",
            self.surfaces.len(), self.steps, rows.join(","), json_string(SCOPE),
        ))
    }
}

/// Keep the common history/control contract while describing the selected
/// spatial owner accurately; historical patch-mean reports remain unchanged.
pub(super) fn history_scope(common: &str) -> String {
    common.replace(
        "Their absolute per-surface emissivity controls use the shared radiosity transpose with the same masked time-scaled thermal multiplier; view factors and geometry stay fixed.",
        "Their absolute per-surface emissivity controls contract the prepared reflected response and the local P1 emission with the same masked time-scaled thermal multiplier; view factors and geometry stay fixed.",
    ).replace(
        "Declared ambient radiation uses its complete nonlinear area-mean state feedback in the same h/source/reference-density transpose; the report also differentiates independent absolute surface emissivities and reservoir temperatures held throughout the selected trajectory.",
        "Declared pointwise radiation uses the exact local P1 fourth-power emission tangent and, for closed enclosures, the complete nonsymmetric reflected irradiation tangent in the same h/source/reference-density transpose. Independent absolute surface emissivities and ambient reservoir temperatures are held throughout the selected trajectory.",
    )
}

const SCOPE: &str = "The complete local P1 fourth-power emission and reflected irradiation enter the same accepted enthalpy transpose as storage, conductivity, contact and physical convection feedback. Parameter partials contract the prepared pointwise outward load with minus the masked time-scaled thermal multiplier. This adds no parameter Krylov or linear solve; construction and physical audits still perform their declared quadrature and prepared radiosity work. Emissivity is frozen at its declared source query, every radiating node retains the claim's temperature support, and query temperature, geometry, view factors and time grid stay fixed. Shared material cards do not merge distinct surface controls. An initial selected peak has zero future radiation sensitivity. These are discrete branch derivatives, not temperature-dependent emissivity, continuum bounds or physical validation.";

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn pointwise_workspace_admits_independent_surface_layout_and_checkpoint_reports() {
        let bare = workspace_entries(12, 16, 2).unwrap() * std::mem::size_of::<f64>();
        let one = storage_bytes(12, 16, 2, 1).unwrap();
        let several = storage_bytes(12, 16, 2, 8).unwrap();
        assert!(one > bare && several > one);
        assert!(storage_bytes(12, 16, 0, 1).is_err());
        assert!(storage_bytes(12, 16, 65, 1).is_err());
        assert!(storage_bytes(12, 16, 2, 0).is_err());
        assert!(storage_bytes(12, 16, 2, usize::MAX).is_err());
    }

    #[test]
    fn pointwise_scope_preserves_history_controls_without_a_patch_mean_claim() {
        let scope = history_scope(super::super::SCOPE);
        assert!(scope.contains("local P1 fourth-power emission tangent"));
        assert!(scope.contains("complete nonsymmetric reflected irradiation tangent"));
        assert!(scope.contains("Initial controls are additive nodal or regional specific-enthalpy offsets"));
        assert!(!scope.contains("nonlinear area-mean state feedback"));
        assert!(!scope.contains("controls use the shared radiosity transpose"));
    }
}

