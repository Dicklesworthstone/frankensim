//! Causal original shell loss heats its original through-thickness caloric cells.
//!
//! Material coefficients and reference geometry stay fixed. The thermal body
//! resolves its original nonlinear conduction and finite caloric storage; it
//! does not feed temperature back into strain, relaxation or material stiffness.
//! Each cell covers the entire original patch area: lateral temperature is
//! homogenized by that existing thermal discretization, never by averaging
//! mechanical layer laws or their independent heat densities.

use fs_exec::Cx;
use fs_matdb::{MaterialCard, QueryPoint};
use fs_material::laminate::{
    LaminaLayer, LaminateError, LaminateSection, LayerLaw, MAX_LAMINA_LAYERS, ReferenceSurface,
};
use fs_material::thermal::conduction::LayerConductionConfig;
use fs_material::thermal::transport::stack::transient::{
    BoundLayeredThermalBody,
    implicit::{
        LayeredThermalSolveControls, LayeredThermalSolveError, LayeredThermalSolveReport,
        LayeredThermalStep, LayeredThermalWorkspace,
    },
};
use fs_material::visco_tensor::source::{ResolvedTensorMaterial, TensorSourceError};
use fs_math::{c64::C64, dd::Dd};
use fs_plate::{
    PlateError,
    shell::{
        CstDktMass, ShellMesh, ShellSupport,
        profile::ProfileBudget,
        survey::validate_midsurface,
        viscosity::{
            ShellViscosityBudget,
            harmonic::{
                TensorPronyShell, TensorPronyShellSection, TensorShellHarmonicError,
                TensorShellMaterialCausalHistory, TensorShellViscoelasticCausalField,
                TensorShellViscoelasticCausalObserver,
            },
        },
    },
};
use fs_thermal::enthalpy::{DepositionError, DepositionSimulation};
use std::sync::Arc;

/// One physical thermal layer's original complete tensor/density selection.
/// Thickness and subdivisions come from the bound thermal body. Orientation is
/// measured from each facet's ex axis, as in the original laminate owner.
#[derive(Clone, Debug)]
pub struct HarmonicLayerSource {
    pub card: Arc<MaterialCard>,
    pub material: ResolvedTensorMaterial,
    pub orientation_rad: f64,
}

#[derive(Debug)]
pub enum LayeredShellHeatError {
    Invalid(&'static str),
    Source(TensorSourceError),
    Laminate(LaminateError),
    Geometry(PlateError),
    Shell(TensorShellHarmonicError),
    Thermal(DepositionError),
    HeatSolve(LayeredThermalSolveError),
    Cancelled,
}
impl core::fmt::Display for LayeredShellHeatError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "causal layered shell heat: {self:?}")
    }
}
impl std::error::Error for LayeredShellHeatError {}
type Result<T> = core::result::Result<T, LayeredShellHeatError>;

fn poll(cx: &Cx<'_>) -> Result<()> {
    cx.checkpoint()
        .map_err(|_| LayeredShellHeatError::Cancelled)
}
fn thermal_poll(cx: &Cx<'_>) -> core::result::Result<(), DepositionError> {
    cx.checkpoint().map_err(|_| DepositionError::Cancelled)
}
fn thermal_error(error: DepositionError) -> LayeredShellHeatError {
    match error {
        DepositionError::Cancelled => LayeredShellHeatError::Cancelled,
        other => LayeredShellHeatError::Thermal(other),
    }
}
fn shell_error(error: TensorShellHarmonicError) -> LayeredShellHeatError {
    match error {
        TensorShellHarmonicError::Plate(PlateError::Cancelled) => LayeredShellHeatError::Cancelled,
        other => LayeredShellHeatError::Shell(other),
    }
}
fn invalid(what: &'static str) -> LayeredShellHeatError {
    LayeredShellHeatError::Invalid(what)
}
fn same(a: f64, b: f64) -> bool {
    a.is_finite()
        && b.is_finite()
        && a > 0.0
        && b > 0.0
        && (a - b).abs() <= 64.0 * f64::EPSILON * a.max(b)
}

/// Original source-bound shell and the exact finite thermal stack it heats.
/// An authored shell with coincident numerical coefficients cannot acquire
/// source provenance through this interface: construction uses the original
/// resolved tensor sources, density and thermal subdivisions directly.
pub struct CausalLayeredShell<'b, 'a> {
    thermal: &'b BoundLayeredThermalBody<'a>,
    sources: Vec<HarmonicLayerSource>,
    cell_sources: Vec<usize>,
    cell_volumes_m3: Vec<f64>,
    shell: TensorPronyShell,
}

/// Independent original constitutive and caloric accounts for this interval.
/// The native thermal account includes generated heat as incoming thermal heat;
/// `other_external_heat_j` removes that internal transfer exactly once.
pub struct CausalLayeredShellHeatReport {
    pub material: TensorShellViscoelasticCausalField,
    /// Original thermal cell order, including each physical subdivision [J].
    pub cell_heat_j: Vec<f64>,
    pub thermal: LayeredThermalSolveReport,
    pub deposited_heat_j: f64,
    pub heat_handoff_residual_j: f64,
    pub other_external_heat_j: f64,
    /// Constitutive plus caloric storage minus stress work and external heat.
    /// Original compensated residuals preserve small changes on large references.
    pub constitutive_thermal_energy_residual_j: f64,
}

/// Both unpublished endpoints. Refusal leaves both accepted owners untouched;
/// choosing this complete candidate adopts each material joule exactly once.
pub struct CausalLayeredShellHeatCandidate<'s, 't> {
    history: TensorShellMaterialCausalHistory<'s>,
    thermal: DepositionSimulation<'t>,
    report: CausalLayeredShellHeatReport,
}
impl<'s, 't> CausalLayeredShellHeatCandidate<'s, 't> {
    pub const fn history(&self) -> &TensorShellMaterialCausalHistory<'s> {
        &self.history
    }
    pub const fn thermal(&self) -> &DepositionSimulation<'t> {
        &self.thermal
    }
    pub const fn report(&self) -> &CausalLayeredShellHeatReport {
        &self.report
    }
    pub fn into_parts(
        self,
    ) -> (
        TensorShellMaterialCausalHistory<'s>,
        DepositionSimulation<'t>,
        CausalLayeredShellHeatReport,
    ) {
        (self.history, self.thermal, self.report)
    }
}

impl<'b, 'a> CausalLayeredShell<'b, 'a> {
    /// Build an original CST/DKT laminate over a planar through-thickness body.
    /// Each original thermal subdivision becomes one full tensor sublamina with
    /// the same physical source and orientation. Facet area, every cell volume
    /// and reference mass must agree. All original material selections must be
    /// temperature independent over the body's complete admitted thermal span.
    /// Nonzero Kelvin terms are retained; their causal plane-stress admission
    /// belongs to the shell/material observer prepared from `shell()`.
    ///
    /// The original shell budgets admit mesh extent and native assembly. A
    /// conservative source-copy score 4096*facets*sum(cells*(branches+1)) must
    /// also fit max_work before layer expansion; no cell or branch is dropped.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        thermal: &'b BoundLayeredThermalBody<'a>,
        sources: Vec<HarmonicLayerSource>,
        mesh: ShellMesh,
        boundary_nodes: &[usize],
        support: ShellSupport,
        reference: ReferenceSurface,
        reference_omega_rad_s: f64,
        budget: ShellViscosityBudget,
        cx: &Cx<'_>,
    ) -> Result<Self> {
        poll(cx)?;
        let body = thermal.prepared_body();
        let conduction = body.conduction();
        if !conduction.geometry().is_planar()
            || body.cell_count() > MAX_LAMINA_LAYERS
            || sources.len() != conduction.layers().len()
            || mesh.tris.is_empty()
            || mesh.tris.len() > budget.max_elements
            || mesh
                .nodes
                .len()
                .checked_mul(6)
                .is_none_or(|n| n > budget.max_full_dofs)
        {
            return Err(invalid(
                "complete bounded source layers and original planar thermal geometry required",
            ));
        }
        let copy_work = sources.iter().zip(body.layer_cell_counts()).try_fold(
            0usize,
            |sum, (source, &count)| {
                count
                    .checked_mul(source.material.law().branches().len() + 1)
                    .and_then(|v| v.checked_mul(4096))
                    .and_then(|v| v.checked_mul(mesh.tris.len()))
                    .and_then(|v| sum.checked_add(v))
            },
        );
        if copy_work.is_none_or(|work| work > budget.max_work) {
            return Err(invalid(
                "original layered shell source expansion exceeds its work budget",
            ));
        }
        validate_midsurface(
            &mesh,
            ProfileBudget {
                max_nodes: budget.max_full_dofs / 6,
                max_triangles: budget.max_elements,
                max_feature_evaluations: 0,
            },
        )
        .map_err(LayeredShellHeatError::Geometry)?;
        let mut area = Dd::ZERO;
        for element in 0..mesh.tris.len() {
            poll(cx)?;
            area = area
                + Dd::from_f64(
                    mesh.facet(element)
                        .map_err(LayeredShellHeatError::Geometry)?
                        .area_m2,
                );
        }
        if !same(area.to_f64(), conduction.area_m2()) {
            return Err(invalid(
                "thermal patch area differs from original shell facet area",
            ));
        }
        let geometries: Vec<_> = body.cell_geometries().collect();
        let mut cell_sources = Vec::with_capacity(body.cell_count());
        let mut layers = Vec::with_capacity(body.cell_count());
        for (layer, ((source, original), &count)) in sources
            .iter()
            .zip(conduction.layers())
            .zip(body.layer_cell_counts())
            .enumerate()
        {
            poll(cx)?;
            let axis = original.temperature_axis();
            let state = source.material.state();
            let context = &original.query_endpoints()[0];
            if source.card.content_hash() != original.material_card_identity()
                || state.card_identity() != source.card.content_hash()
                || !same_context(&source.material, context, axis)
                || !source.orientation_rad.is_finite()
            {
                return Err(invalid(
                    "original tensor/thermal card, fixed query or layer orientation differs",
                ));
            }
            // The original body permits capacity and conductivity to name
            // extra fixed coordinates independently. This combined specimen
            // must establish the complete SAME physical context for both.
            for component in body.storages()[layer].components() {
                if component.material_card_identity() != state.card_identity()
                    || component.temperature_axis() != axis
                    || component
                        .query_endpoints()
                        .iter()
                        .any(|point| !same_context(&source.material, point, axis))
                {
                    return Err(invalid(
                        "caloric and mechanical layers differ in their complete original source context",
                    ));
                }
            }
            for bundle in
                core::iter::once(state).chain(source.material.kelvin_source().map(|k| k.state()))
            {
                for property in bundle.properties() {
                    poll(cx)?;
                    if !property.has_only_fixed_external_conditions(Some(axis), true)
                        || property.source_axes().iter().any(|(a, _)| a == axis)
                    {
                        return Err(invalid(
                            "causal self-heating requires validity-wide temperature-independent material coefficients and density",
                        ));
                    }
                }
            }
            for temperature in body.temperature_span_k() {
                check_source(source, axis, temperature, cx)?;
            }
            let thickness = original.thickness_m() / count as f64;
            let volume = area.to_f64() * thickness;
            let mass = volume * source.material.density_kg_m3();
            for _ in 0..count {
                poll(cx)?;
                let cell = layers.len();
                if !same(volume, geometries[cell].volume_m3)
                    || !same(mass, body.layer_cell_masses_kg()[layer])
                {
                    return Err(invalid(
                        "original shell and thermal cell reference mass or volume differs",
                    ));
                }
                layers.push(LaminaLayer {
                    thickness_m: thickness,
                    orientation_rad: source.orientation_rad,
                    density_kg_m3: source.material.density_kg_m3(),
                    material: Some(state.clone()),
                    law: LayerLaw::Prony(source.material.law().clone()),
                });
                cell_sources.push(layer);
            }
        }
        let section =
            LaminateSection::try_new(layers, reference).map_err(LayeredShellHeatError::Laminate)?;
        let mut sections = Vec::with_capacity(mesh.tris.len());
        for _ in &mesh.tris {
            poll(cx)?;
            sections.push(
                TensorPronyShellSection::from_laminate(section.clone())
                    .map_err(LayeredShellHeatError::Geometry)?,
            );
        }
        let shell = TensorPronyShell::assemble(
            mesh,
            sections,
            boundary_nodes,
            support,
            CstDktMass::Consistent,
            reference_omega_rad_s,
            budget,
            cx,
        )
        .map_err(shell_error)?;
        poll(cx)?;
        Ok(Self {
            thermal,
            sources,
            cell_sources,
            cell_volumes_m3: geometries.iter().map(|g| g.volume_m3).collect(),
            shell,
        })
    }

    pub const fn shell(&self) -> &TensorPronyShell {
        &self.shell
    }
    pub fn source_layers(&self) -> &[HarmonicLayerSource] {
        &self.sources
    }
    /// Original physical source-layer index for each retained thermal cell.
    pub fn cell_source_layers(&self) -> &[usize] {
        &self.cell_sources
    }

    fn admit_state(&self, state: &DepositionSimulation<'_>, cx: &Cx<'_>) -> Result<()> {
        let temperatures = self
            .thermal
            .cell_temperatures(state, || thermal_poll(cx))
            .map_err(thermal_error)?;
        for (&temperature, &layer) in temperatures.iter().zip(&self.cell_sources) {
            let axis = self.thermal.prepared_body().conduction().layers()[layer].temperature_axis();
            check_source(&self.sources[layer], axis, temperature, cx)?;
        }
        Ok(())
    }

    /// Evolve original memory once and deposit only its independently integrated
    /// physical material heat into the original thermal cells. Startup, coherent
    /// interference and held-displacement relaxation retain their true memories.
    /// Thermal iteration reuses that one fixed heat input; it never reruns or
    /// debits the mechanical trajectory. Equal source laws are not sufficient:
    /// the observer must borrow this exact owned shell.
    ///
    /// Banks/bias follow the causal shell's represented-coordinate convention.
    /// The original full nodal and material continuity checks remain active.
    /// The original thermal clock is f64: both endpoints must equal the complete
    /// DD material clock without rounding. Pending events retain native refusal.
    ///
    /// Actual joules become interval-mean powers on the original backward-Euler
    /// thermal residual; refine time intervals to resolve beat-driven transport.
    /// The body homogenizes lateral temperature over its admitted patch area.
    /// No source-temperature feedback, thermal strain, kinetic/boundary work,
    /// force equilibrium, or reciprocal thermoelastic entropy is inferred.
    #[allow(clippy::too_many_arguments)]
    pub fn heat_candidate<'s, 't>(
        &self,
        thermal: &DepositionSimulation<'t>,
        target_s: f64,
        observer: &TensorShellViscoelasticCausalObserver<'s>,
        accepted: &TensorShellMaterialCausalHistory<'s>,
        frequencies: &[f64],
        displacements: &[&[C64]],
        bias: &[f64],
        conduction: LayerConductionConfig,
        controls: LayeredThermalSolveControls,
        workspace: &mut LayeredThermalWorkspace,
        energy_tolerance_j: f64,
        cx: &Cx<'_>,
    ) -> Result<CausalLayeredShellHeatCandidate<'s, 't>> {
        poll(cx)?;
        let duration = target_s - thermal.time_s();
        if !core::ptr::eq(observer.source(), &self.shell)
            || !target_s.is_finite()
            || !duration.is_finite()
            || duration <= 0.0
            || accepted.time_s_dd() != Dd::from_f64(thermal.time_s())
            || accepted.time_s_dd() + Dd::from_f64(duration) != Dd::from_f64(target_s)
            || !energy_tolerance_j.is_finite()
            || energy_tolerance_j <= 0.0
            || !controls.total_tolerance_j.is_finite()
            || controls.total_tolerance_j > energy_tolerance_j
        {
            return Err(invalid(
                "same original shell, shared complete clocks and bounded thermal energy controls required",
            ));
        }
        self.admit_state(thermal, cx)?;
        let native = LayeredThermalStep::new(self.thermal, thermal, target_s, conduction, || {
            thermal_poll(cx)
        })
        .map_err(thermal_error)?;
        let step = observer
            .advance_window(accepted, frequencies, displacements, bias, duration, cx)
            .map_err(shell_error)?;
        let field = step.field;
        let n = self.cell_sources.len();
        let mut heat = vec![Dd::ZERO; n];
        let mut volumes = vec![Dd::ZERO; n];
        let mut counts = vec![0usize; n];
        for point in &field.points {
            poll(cx)?;
            if point.layer >= n
                || point.element >= self.shell.mesh().tris.len()
                || !point.volume_m3.is_finite()
                || point.volume_m3 <= 0.0
                || !point.material.material_heat_j_m3.is_finite()
                || point.material.material_heat_j_m3 < 0.0
            {
                return Err(invalid(
                    "causal shell heat lost its original positive layer quadrature",
                ));
            }
            let point_heat =
                Dd::from_f64(point.volume_m3) * Dd::from_f64(point.material.material_heat_j_m3);
            if point.material.material_heat_j_m3 > 0.0 && point_heat == Dd::ZERO {
                return Err(invalid(
                    "original positive layer heat is below representable cell energy",
                ));
            }
            heat[point.layer] = heat[point.layer] + point_heat;
            volumes[point.layer] = volumes[point.layer] + Dd::from_f64(point.volume_m3);
            counts[point.layer] += 1;
        }
        let mut cell_heat_j = Vec::with_capacity(n);
        let mut deposited = Dd::ZERO;
        for cell in 0..n {
            poll(cx)?;
            let energy = heat[cell].to_f64();
            let power = energy / duration;
            if counts[cell] != 6 * self.shell.mesh().tris.len()
                || !same(volumes[cell].to_f64(), self.cell_volumes_m3[cell])
                || !energy.is_finite()
                || energy < 0.0
                || !power.is_finite()
                || energy > 0.0 && power == 0.0
            {
                return Err(invalid(
                    "incomplete thermal cell volume or unrepresentable original material heat",
                ));
            }
            // Same interval-power conversion as the original thermal owner.
            deposited = deposited + Dd::from_f64(power) * Dd::from_f64(duration);
            cell_heat_j.push(energy);
        }
        let handoff = (deposited - Dd::from_f64(field.material_heat_j)).to_f64();
        if !handoff.is_finite()
            || !field.energy_balance_residual_j.is_finite()
            || handoff.abs() + field.energy_balance_residual_j.abs() > energy_tolerance_j
        {
            return Err(invalid(
                "original constitutive balance or cell heat handoff exceeds its allowance",
            ));
        }
        let native = native
            .with_cell_energy_j(&cell_heat_j, || thermal_poll(cx))
            .map_err(thermal_error)?;
        let (next_thermal, _, report) = native
            .solve(controls, workspace, || thermal_poll(cx))
            .map_err(|error| match error {
                LayeredThermalSolveError::Thermal(DepositionError::Cancelled) => {
                    LayeredShellHeatError::Cancelled
                }
                other => LayeredShellHeatError::HeatSolve(other),
            })?;
        self.admit_state(&next_thermal, cx)?;
        let external = (Dd::from_f64(report.energy.external_heat_j) - deposited).to_f64();
        let combined = (Dd::from_f64(report.energy.residual_j)
            - Dd::from_f64(field.energy_balance_residual_j)
            + Dd::from_f64(handoff))
        .to_f64();
        if !external.is_finite()
            || !combined.is_finite()
            || report.energy.incoming_mass_kg != 0.0
            || report.energy.incoming_energy_j != 0.0
            || report.energy.residual_j.abs()
                + handoff.abs()
                + field.energy_balance_residual_j.abs()
                > energy_tolerance_j
        {
            return Err(invalid(
                "causal shell/caloric energy account exceeds its physical allowance",
            ));
        }
        poll(cx)?;
        Ok(CausalLayeredShellHeatCandidate {
            history: step.history,
            thermal: next_thermal,
            report: CausalLayeredShellHeatReport {
                material: field,
                cell_heat_j,
                thermal: report,
                deposited_heat_j: deposited.to_f64(),
                heat_handoff_residual_j: handoff,
                other_external_heat_j: external,
                constitutive_thermal_energy_residual_j: combined,
            },
        })
    }
}

// Source axes establish temperature independence; equal samples alone cannot.
// Re-evaluation also preserves original selection, query and calibration domain.
fn same_context(source: &ResolvedTensorMaterial, context: &QueryPoint, axis: &str) -> bool {
    source.state().axis_quantities() == context.axis_quantities()
        && source
            .state()
            .query_point()
            .iter()
            .filter(|(a, _)| a != axis)
            .map(|(a, v)| (a, *v))
            .eq(context
                .axes()
                .iter()
                .filter(|(a, _)| a.as_str() != axis)
                .map(|(a, v)| (a, *v)))
}

fn check_source(
    source: &HarmonicLayerSource,
    axis: &str,
    temperature: f64,
    cx: &Cx<'_>,
) -> Result<()> {
    poll(cx)?;
    let query = source
        .material
        .state()
        .query()
        .with(axis, temperature)
        .map_err(|_| invalid("invalid original layer source temperature"))?;
    let refreshed = source
        .material
        .at_query_point(&source.card, &query)
        .map_err(LayeredShellHeatError::Source)?;
    if refreshed.law() != source.material.law()
        || refreshed.density_kg_m3().to_bits() != source.material.density_kg_m3().to_bits()
    {
        return Err(invalid(
            "thermal proposal changed an original frozen shell material or density",
        ));
    }
    Ok(())
}
