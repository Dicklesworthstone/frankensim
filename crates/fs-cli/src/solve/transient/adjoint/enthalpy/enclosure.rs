//! Emissivity controls of the complete native enclosure history transpose.
//! The core owns both the physical thermal tangent and the radiosity transpose.

use fs_conduction::radiation::{EnclosureRadiosityLinearization, GrayDiffuseEnclosure};

use super::{
    ConductionProblem, Cx, Instant, LinearConfig, SolveRefusal, bad, finite,
    json_string, lower, number, poll, radiation, zeros,
};

pub(super) struct Controls {
    emissivities: Vec<f64>,
    steps: usize,
    radiosity_iterations: usize,
    maximum_relative_residual: f64,
}

impl Controls {
    pub(super) fn new(model: &GrayDiffuseEnclosure) -> Result<Self, SolveRefusal> {
        if !(2..=64).contains(&model.surfaces().len()) {
            return Err(bad("native enclosure controls require 2..=64 ordered surfaces"));
        }
        Ok(Self {
            emissivities: zeros(model.surfaces().len())?, steps: 0,
            radiosity_iterations: 0, maximum_relative_residual: 0.0,
        })
    }

    /// The caller has solved the COMPLETE physical transpose. load_weights
    /// is dt*lambda with prescribed residual rows masked by that same owner.
    /// Outward flux enters the residual positively, so its control bar is -G^T
    /// load_weights. Reflections are differentiated by the shared radiosity owner.
    #[allow(clippy::too_many_arguments)]
    pub(super) fn add(
        &mut self,
        cx: &Cx<'_>,
        model: &GrayDiffuseEnclosure,
        problem: ConductionProblem<'_>,
        temperature: &[f64],
        load_weights: &[f64],
        heat_tolerance_w: f64,
        linear: LinearConfig,
        deadline: Option<(Instant, f64)>,
    ) -> Result<(), SolveRefusal> {
        poll(cx, deadline)?;
        if temperature.len() != problem.mesh.vertex_count()
            || load_weights.len() != temperature.len()
            || model.surfaces().len() != self.emissivities.len()
        {
            return Err(bad("enclosure controls lost their accepted endpoint or patch order"));
        }
        let mut means = Vec::with_capacity(self.emissivities.len());
        let mut flux_weights = zeros(self.emissivities.len())?;
        for (index, surface) in model.surfaces().iter().enumerate() {
            poll(cx, deadline)?;
            means.push(surface.mean_temperature(problem.mesh, temperature).map_err(lower)?);
            for (ordinal, &slot) in surface.face_slots().iter().enumerate() {
                if ordinal % 256 == 0 { poll(cx, deadline)?; }
                let face = &problem.mesh.boundary()[slot];
                for &vertex in &face.vertices {
                    flux_weights[index] = finite(
                        flux_weights[index] - finite(load_weights[vertex as usize])? * face.area / 3.0,
                    )?;
                }
            }
        }
        let radiosity = EnclosureRadiosityLinearization::new(
            cx, model, &means, heat_tolerance_w,
        ).map_err(lower)?;
        let gradient = radiosity.pullback(
            cx, &flux_weights, linear.tolerance, linear.max_iterations,
        ).map_err(lower)?;
        if gradient.log_emissivities.len() != self.emissivities.len() {
            return Err(bad("enclosure emissivity transpose lost its declared ordering"));
        }
        for ((sum, &bar), surface) in self.emissivities.iter_mut()
            .zip(&gradient.log_emissivities).zip(model.surfaces())
        {
            *sum = finite(*sum + finite(bar / surface.emissivity().value())?)?;
        }
        self.radiosity_iterations = self.radiosity_iterations.checked_add(gradient.iterations)
            .ok_or_else(|| bad("enclosure radiosity reverse iteration count overflows"))?;
        self.maximum_relative_residual = self.maximum_relative_residual.max(gradient.relative_residual);
        self.steps = self.steps.checked_add(1)
            .ok_or_else(|| bad("enclosure reverse step count overflows"))?;
        poll(cx, deadline)
    }

    pub(super) fn iterations(&self) -> usize { self.radiosity_iterations }

    pub(super) fn relative_residual(&self) -> f64 { self.maximum_relative_residual }

    pub(super) fn fragment(
        &self,
        lowering: &radiation::LoweredRadiation,
        model: &GrayDiffuseEnclosure,
        expected_steps: usize,
        max_feedback_entries: usize,
        workspace_bytes: usize,
    ) -> Result<String, SolveRefusal> {
        if self.steps != expected_steps || lowering.surface_count() != self.emissivities.len()
            || model.surfaces().len() != self.emissivities.len()
        {
            return Err(bad("enclosure controls do not cover the selected reverse prefix"));
        }
        let mut rows = Vec::with_capacity(self.emissivities.len());
        for (index, ((name, target), surface)) in lowering.surface_labels()
            .zip(model.surfaces()).enumerate()
        {
            if target != surface.name() {
                return Err(bad("enclosure control identity differs from the view-factor order"));
            }
            let emissivity = surface.emissivity();
            let directions = if emissivity.value() == 1.0 {
                "nonpositive-at-emissivity-one"
            } else {
                "two-sided-within-zero-exclusive-one-inclusive"
            };
            rows.push(format!(
                "{{\"name\":{},\"target\":{},\"card\":{},\"material_state\":{},\"emissivity_receipt\":{},\"emissivity\":{},\"query_temperature_k\":{},\"dtemperature_demissivity_k\":{},\"emissivity_feasible_directions\":{}}}",
                json_string(name), json_string(target),
                json_string(&emissivity.card_identity().to_hex()),
                json_string(emissivity.material_state()),
                json_string(&emissivity.receipt().content_hash().to_hex()),
                number(emissivity.value())?, number(emissivity.temperature_k())?,
                number(self.emissivities[index])?, json_string(directions),
            ));
        }
        Ok(format!(
            ",\"radiation\":{{\"method\":\"complete-closed-gray-diffuse-enthalpy-tangent\",\"patch_count\":{},\"control_steps\":{},\"max_feedback_entries\":{},\"additional_workspace_bytes\":{},\"radiosity_adjoint_krylov_iterations\":{},\"maximum_radiosity_dual_relative_residual\":{},\"parameters\":\"absolute-per-surface-emissivity-throughout-selected-trajectory\",\"view_factor_identity\":{},\"surfaces\":[{}]}}",
            self.emissivities.len(), self.steps, max_feedback_entries, workspace_bytes,
            self.radiosity_iterations, number(self.maximum_relative_residual)?,
            json_string(&model.view_factors().identity().to_hex()), rows.join(","),
        ))
    }
}
