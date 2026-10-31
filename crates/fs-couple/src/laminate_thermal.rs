//! Spatial heat conduction driving the original layered MITC shell.
//!
//! Every original thermal control volume becomes one mechanical sublamina at
//! its actual temperature. Stiffness and integrated expansion retain their own
//! source claims; reference geometry, mass and material directions stay fixed.
//! A through-thickness gradient is never replaced by an average temperature.
//!
//! This is a one-way, quasistatic thermoelastic approximation: the original
//! thermal owner advances heat and the native shell solves compatible motion.
//! Its heat audit excludes mechanical feedback, thermal entropy of elasticity,
//! and stiffness-change work. The two accounts are exposed separately, never
//! advertised as a reciprocal energy balance. Mechanical interfaces are perfectly
//! bonded; the original thermal contact resistances remain independently owned.

use fs_material::laminate::thermal::source::{
    LaminaSourceInput, LaminateThermalSource, LaminateThermalSourceError, ResolvedLaminateThermal,
};
use fs_material::laminate::{MAX_LAMINA_LAYERS, ReferenceSurface};
use fs_material::thermal::conduction::LayerConductionConfig;
use fs_material::thermal::transport::stack::transient::{
    BoundLayeredThermalBody,
    implicit::{
        LayeredThermalSolveControls, LayeredThermalSolveError, LayeredThermalSolveReport,
        LayeredThermalStep, LayeredThermalWorkspace,
    },
};
use fs_math::det;
use fs_plate::{
    PlateError,
    shell::{
        ShellMesh,
        mitc::{
            MitcBudget, MitcShellModel, NodalDirectors, ShellConstraint, ShellMass,
            assemble_mitc3p,
            thermal::{LaminateThermalField, ThermalShellResponse},
        },
    },
};
use fs_thermal::enthalpy::{DepositionError, DepositionSimulation};

/// Explicit reference mesh, directors, supports and numerical shell choices.
#[derive(Clone, Debug)]
pub struct LaminateShellDiscretization {
    pub mesh: ShellMesh,
    pub directors: NodalDirectors,
    pub constraints: Vec<ShellConstraint>,
    /// Whole-section x direction relative to each facet ex/ey frame [rad].
    pub section_angles_rad: Vec<f64>,
    pub mass: ShellMass,
    pub budget: MitcBudget,
    pub drilling_alpha: f64,
}

/// Caller-selected validity envelope for this linear kinematic approximation.
/// These are numerical/model limits, not invented material failure strengths.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct LaminateShellLimits {
    pub maximum_translation_m: f64,
    pub maximum_rotation_rad: f64,
    /// Norm of represented elastic xx,yy,xy,xz,yz strain; thickness-normal
    /// strain is not inferred from a reduced plane-stress material source.
    pub maximum_reduced_elastic_strain: f64,
}

#[derive(Debug)]
pub enum LaminateShellError {
    Invalid(&'static str),
    Source(LaminateThermalSourceError),
    Thermal(DepositionError),
    HeatSolve(LayeredThermalSolveError),
    Shell(PlateError),
    Envelope {
        quantity: &'static str,
        actual: f64,
        limit: f64,
    },
    Cancelled,
}
impl core::fmt::Display for LaminateShellError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "FS-COUPLE-LAMINATE-THERMAL: {self:?}")
    }
}
impl std::error::Error for LaminateShellError {}

/// A single reference shell patch and its original axial thermal stack.
/// Temperature varies through thickness; lateral heat transport belongs to a
/// different thermal discretization. The stack area equals the mesh facet area.
pub struct LayeredThermalShell<'b, 'a> {
    thermal: &'b BoundLayeredThermalBody<'a>,
    source: LaminateThermalSource,
    discretization: LaminateShellDiscretization,
    limits: LaminateShellLimits,
}

/// One accepted source temperature field and its compatible native deformation.
/// The native pencil can serve existing vibration consumers; it does not include
/// geometric prestress stiffness or transfer a previous vibration history.
pub struct LaminateShellSnapshot {
    time_s: f64,
    temperatures_k: Vec<f64>,
    material: ResolvedLaminateThermal,
    model: MitcShellModel,
    response: ThermalShellResponse,
    strain_bound: f64,
}
impl LaminateShellSnapshot {
    pub const fn time_s(&self) -> f64 {
        self.time_s
    }
    pub fn temperatures_k(&self) -> &[f64] {
        &self.temperatures_k
    }
    pub const fn material(&self) -> &ResolvedLaminateThermal {
        &self.material
    }
    pub const fn model(&self) -> &MitcShellModel {
        &self.model
    }
    pub const fn response(&self) -> &ThermalShellResponse {
        &self.response
    }
    pub const fn reduced_elastic_strain_upper_bound(&self) -> f64 {
        self.strain_bound
    }
}

/// Atomic endpoint of a native heat step followed by compatible deformation.
/// `heat_report` is the thermal owner's heat account only, not combined energy.
pub struct LaminateShellStep<'m> {
    pub thermal: DepositionSimulation<'m>,
    pub shell: LaminateShellSnapshot,
    pub heat_report: LayeredThermalSolveReport,
}

impl<'b, 'a> LayeredThermalShell<'b, 'a> {
    /// Bind one mechanical source input per original thermal material layer.
    /// Each thermal subdivision is retained as a separate sublamina; the
    /// existing 128-lamina limit is checked before expansion or allocation.
    pub fn try_new(
        thermal: &'b BoundLayeredThermalBody<'a>,
        layers: Vec<LaminaSourceInput>,
        reference: ReferenceSurface,
        discretization: LaminateShellDiscretization,
        limits: LaminateShellLimits,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, LaminateShellError> {
        poll(&mut cancelled)?;
        let body = thermal.prepared_body();
        if body.cell_count() > MAX_LAMINA_LAYERS {
            return Err(LaminateShellError::Invalid(
                "thermal resolution exceeds native laminate layer budget",
            ));
        }
        if layers.len() != body.layer_cell_counts().len() {
            return Err(LaminateShellError::Invalid(
                "one mechanical source per original thermal layer required",
            ));
        }
        if ![
            limits.maximum_translation_m,
            limits.maximum_rotation_rad,
            limits.maximum_reduced_elastic_strain,
            discretization.drilling_alpha,
        ]
        .iter()
        .all(|v| v.is_finite() && *v > 0.0)
        {
            return Err(LaminateShellError::Invalid(
                "positive finite shell limits and drilling strength required",
            ));
        }
        let mesh = &discretization.mesh;
        if mesh.nodes.len() > discretization.budget.max_nodes
            || mesh.tris.len() > discretization.budget.max_triangles
            || discretization.section_angles_rad.len() != mesh.tris.len()
            || discretization
                .section_angles_rad
                .iter()
                .any(|v| !v.is_finite())
        {
            return Err(LaminateShellError::Invalid(
                "shell mesh/angle count exceeds declared budget",
            ));
        }
        fs_plate::shell::survey::validate_midsurface(
            mesh,
            fs_plate::shell::profile::ProfileBudget {
                max_nodes: discretization.budget.max_nodes,
                max_triangles: discretization.budget.max_triangles,
                max_feature_evaluations: 0,
            },
        )
        .map_err(LaminateShellError::Shell)?;
        let (mut area, mut correction) = (0.0, 0.0);
        for i in 0..mesh.tris.len() {
            poll(&mut cancelled)?;
            let value = mesh.facet(i).map_err(LaminateShellError::Shell)?.area_m2 - correction;
            let next = area + value;
            correction = (next - area) - value;
            area = next;
        }
        if !same(area, body.conduction().area_m2()) {
            return Err(LaminateShellError::Invalid(
                "thermal area must equal the original shell patch area",
            ));
        }
        let mut inputs = Vec::with_capacity(body.cell_count());
        for ((input, layer), &count) in layers
            .into_iter()
            .zip(body.conduction().layers())
            .zip(body.layer_cell_counts())
        {
            poll(&mut cancelled)?;
            let context = &layer.query_endpoints()[0];
            if input.card.content_hash() != layer.material_card_identity()
                || !same(input.thickness_m, layer.thickness_m())
                || layer.temperature_axis() != "T"
                || input.reference.axis_quantities() != context.axis_quantities()
                || !input.reference.axes().keys().eq(context.axes().keys())
                || input
                    .reference
                    .axes()
                    .iter()
                    .any(|(axis, value)| axis != "T" && context.axes().get(axis) != Some(value))
            {
                return Err(LaminateShellError::Invalid(
                    "mechanical and thermal layer card, geometry or fixed context differ",
                ));
            }
            for _ in 0..count {
                inputs.push(LaminaSourceInput {
                    thickness_m: input.thickness_m / count as f64,
                    ..input.clone()
                });
            }
        }
        let source = LaminateThermalSource::try_new(inputs, reference, &mut cancelled)
            .map_err(LaminateShellError::Source)?;
        let mut cell = 0;
        for (layer, &count) in body.layer_cell_counts().iter().enumerate() {
            for _ in 0..count {
                let original = &source.layers()[cell];
                let mass = area * original.input().thickness_m * original.reference_density_kg_m3();
                if !same(mass, body.layer_cell_masses_kg()[layer]) {
                    return Err(LaminateShellError::Invalid(
                        "thermal and mechanical reference cell mass differ",
                    ));
                }
                cell += 1;
            }
        }
        poll(&mut cancelled)?;
        Ok(Self {
            thermal,
            source,
            discretization,
            limits,
        })
    }

    pub const fn source(&self) -> &LaminateThermalSource {
        &self.source
    }
    pub const fn discretization(&self) -> &LaminateShellDiscretization {
        &self.discretization
    }
    pub const fn limits(&self) -> LaminateShellLimits {
        self.limits
    }

    /// Read any accepted original thermal state and solve the current native
    /// shell under explicit world nodal loads. Reassembly uses current elastic
    /// properties while retaining original geometry and every mass moment.
    /// Failure or cancellation leaves the input thermal timeline unchanged.
    pub fn resolve(
        &self,
        thermal: &DepositionSimulation<'_>,
        external_load_6n: &[f64],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<LaminateShellSnapshot, LaminateShellError> {
        poll(&mut cancelled)?;
        let temperatures = self
            .thermal
            .cell_temperatures(thermal, || thermal_poll(&mut cancelled))
            .map_err(LaminateShellError::Thermal)?;
        let material = self
            .source
            .at_temperatures(&temperatures, &mut cancelled)
            .map_err(LaminateShellError::Source)?;
        let prepared = material.prepared();
        let mut fields = Vec::with_capacity(self.discretization.mesh.tris.len());
        let mut sections = Vec::with_capacity(fields.capacity());
        for &angle in &self.discretization.section_angles_rad {
            poll(&mut cancelled)?;
            let field = LaminateThermalField::try_new(prepared.clone(), angle)
                .map_err(LaminateShellError::Shell)?;
            sections.push(*field.section());
            fields.push(field);
        }
        let d = &self.discretization;
        let model = assemble_mitc3p(
            &d.mesh,
            &d.directors,
            &sections,
            &d.constraints,
            d.mass,
            d.budget,
            d.drilling_alpha,
            &mut cancelled,
        )
        .map_err(LaminateShellError::Shell)?;
        let loading = model
            .prepare_laminate_thermal(&fields, &mut cancelled)
            .map_err(LaminateShellError::Shell)?;
        let response = loading
            .solve_static(external_load_6n, &mut cancelled)
            .map_err(LaminateShellError::Shell)?;
        for nodal in response.displacement_6n.chunks_exact(6) {
            poll(&mut cancelled)?;
            envelope(
                "translation",
                norm(&nodal[..3]),
                self.limits.maximum_translation_m,
            )?;
            envelope(
                "rotation",
                norm(&nodal[3..]),
                self.limits.maximum_rotation_rad,
            )?;
        }
        let mut strain_bound = 0.0_f64;
        for (element, tri) in d.mesh.tris.iter().enumerate() {
            poll(&mut cancelled)?;
            let nodal = core::array::from_fn(|i| response.displacement_6n[6 * tri[i / 6] + i % 6]);
            let strain = loading
                .strain_field(element, &nodal, &mut cancelled)
                .map_err(LaminateShellError::Shell)?;
            for layer in 0..self.source.layers().len() {
                let bound = fields[element]
                    .elastic_strain_bound(layer, &strain, &mut cancelled)
                    .map_err(LaminateShellError::Shell)?;
                envelope(
                    "reduced elastic strain",
                    bound,
                    self.limits.maximum_reduced_elastic_strain,
                )?;
                strain_bound = strain_bound.max(bound);
            }
        }
        poll(&mut cancelled)?;
        Ok(LaminateShellSnapshot {
            time_s: thermal.time_s(),
            temperatures_k: temperatures,
            material,
            model,
            response,
            strain_bound,
        })
    }

    /// Advance an isolated original stack with its sparse implicit heat owner,
    /// then solve its compatible shell. Both endpoints are returned only after
    /// source, mechanical and envelope admission succeeds. Existing heat events
    /// and timestep rules remain owned by `LayeredThermalStep`.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_to<'m>(
        &self,
        initial: &DepositionSimulation<'m>,
        target_s: f64,
        external_load_6n: &[f64],
        conduction: LayerConductionConfig,
        controls: LayeredThermalSolveControls,
        workspace: &mut LayeredThermalWorkspace,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<LaminateShellStep<'m>, LaminateShellError> {
        poll(&mut cancelled)?;
        let step = LayeredThermalStep::new(self.thermal, initial, target_s, conduction, || {
            thermal_poll(&mut cancelled)
        })
        .map_err(LaminateShellError::Thermal)?;
        let (thermal, _, heat_report) = step
            .solve(controls, workspace, || thermal_poll(&mut cancelled))
            .map_err(LaminateShellError::HeatSolve)?;
        let shell = self.resolve(&thermal, external_load_6n, &mut cancelled)?;
        poll(&mut cancelled)?;
        Ok(LaminateShellStep {
            thermal,
            shell,
            heat_report,
        })
    }
}

fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), LaminateShellError> {
    if cancelled() {
        Err(LaminateShellError::Cancelled)
    } else {
        Ok(())
    }
}
fn thermal_poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), DepositionError> {
    if cancelled() {
        Err(DepositionError::Cancelled)
    } else {
        Ok(())
    }
}
fn same(a: f64, b: f64) -> bool {
    a.is_finite()
        && b.is_finite()
        && a > 0.0
        && b > 0.0
        && (a - b).abs() <= 128.0 * f64::EPSILON * a.max(b)
}
fn norm(values: &[f64]) -> f64 {
    let scale = values.iter().fold(0.0_f64, |m, v| m.max(v.abs()));
    if scale == 0.0 {
        0.0
    } else {
        scale * det::sqrt(values.iter().map(|v| (v / scale).powi(2)).sum())
    }
}
fn envelope(quantity: &'static str, actual: f64, limit: f64) -> Result<(), LaminateShellError> {
    if !actual.is_finite() || actual > limit {
        Err(LaminateShellError::Envelope {
            quantity,
            actual,
            limit,
        })
    } else {
        Ok(())
    }
}
