//! Source-backed natural strain in the existing conservative TET10 owner.
//!
//! Material laws remain at their admitted held temperature and density remains
//! at the original source datum. A declared manufacturing map is composed before
//! the source thermal map: `G = G_thermal * G_before_thermal`. The same total G
//! enters the native objective material work, exact-volume constraint and
//! original finite-bulk storage. Original geometry and consistent mass do not
//! change. No constitutive parameters, frame or temperature path are inferred.
//!
//! This is a held-reference conservative elastic trajectory. Prescribed
//! temperature evolution, heat, viscous memory and reference-change work need
//! their own admitted events; replacing a field does not mutate an accepted body.

use super::{FiniteSolidMaterial, ThermalReference};
use fs_exec::Cx;
use fs_material::MaterialError;
use fs_material::hyper::Hyperelastic;
use fs_material::hyper::discrete_gradient::DiscreteGradientSettings;
use fs_material::hyper::reference::ReferenceDeformation;
use fs_solid::hyper3::mixed::conservative::{
    ConservativeFormulation, ConservativeReferenceField, ConservativeTetDefinition,
    ConservativeTetEnergy, ConservativeTetError, ConservativeTetSolid, InertialEnergy,
};
use fs_solid::hyper3::mixed::{MixedLawField, MixedUpSettings};
use fs_solid::linear3::TetAssemblyBudget;
use fs_solver::nonlinear::NewtonReport;

#[derive(Debug)]
pub enum ConservativeThermalError {
    Invalid(&'static str),
    Reference {
        element: usize,
        source: MaterialError,
    },
    Native(ConservativeTetError),
    Cancelled,
}

impl From<ConservativeTetError> for ConservativeThermalError {
    fn from(value: ConservativeTetError) -> Self {
        Self::Native(value)
    }
}

impl core::fmt::Display for ConservativeThermalError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(
            f,
            "source-backed conservative thermal solid refused: {self:?}"
        )
    }
}

impl core::error::Error for ConservativeThermalError {}

fn checkpoint(cx: &Cx<'_>) -> Result<(), ConservativeThermalError> {
    cx.checkpoint()
        .map_err(|_| ConservativeThermalError::Cancelled)
}

fn compose_reference(
    material: &FiniteSolidMaterial,
    before_thermal: ReferenceDeformation,
    element: usize,
) -> Result<ReferenceDeformation, ConservativeThermalError> {
    let thermal = match material.thermal_reference() {
        None => ReferenceDeformation::default(),
        Some(ThermalReference::Mean(source)) => {
            let s = source.stretch();
            ReferenceDeformation::try_new([s, 0.0, 0.0, 0.0, s, 0.0, 0.0, 0.0, s])
                .map_err(|source| ConservativeThermalError::Reference { element, source })?
        }
        Some(ThermalReference::Expansion(source)) => *source.reference_deformation(),
    };
    before_thermal
        .followed_by(thermal)
        .map_err(|source| ConservativeThermalError::Reference { element, source })
}

#[derive(Debug)]
enum SourceField {
    Uniform {
        material: FiniteSolidMaterial,
        before_thermal: ReferenceDeformation,
        reference: ReferenceDeformation,
    },
    PerElement {
        materials: Vec<FiniteSolidMaterial>,
        laws: Vec<Hyperelastic>,
        before_thermal: Vec<ReferenceDeformation>,
        references: Vec<ReferenceDeformation>,
    },
}

/// Immutable source owners, composed natural strains and original inventory.
/// Uniform fields keep one source/law/reference; heterogeneous fields retain
/// exact original element order. Native law copies are prepared only once.
#[derive(Debug)]
pub struct ConservativeThermalField {
    source: SourceField,
    reference_density_kg_m3: Vec<f64>,
}

/// The existing native controls, with no inferred formulation or quadrature.
#[derive(Clone, Copy, Debug)]
pub struct ConservativeThermalPreparation {
    pub formulation: ConservativeFormulation,
    pub budget: TetAssemblyBudget,
    pub quadrature_points_per_dir: usize,
    pub maximum_quadrature_points: usize,
    pub maximum_mass_assembly_bytes: usize,
    pub material_settings: DiscreteGradientSettings,
}

impl ConservativeThermalField {
    /// Repeat one admitted source over an explicit nonempty original mesh.
    /// The output frame of `before_thermal` must be the thermal source's input
    /// frame; the caller supplies that physical convention explicitly.
    pub fn uniform(
        material: FiniteSolidMaterial,
        before_thermal: ReferenceDeformation,
        element_count: usize,
        cx: &Cx<'_>,
    ) -> Result<Self, ConservativeThermalError> {
        checkpoint(cx)?;
        if element_count == 0 {
            return Err(ConservativeThermalError::Invalid(
                "at least one original element is required",
            ));
        }
        let reference = compose_reference(&material, before_thermal, 0)?;
        let density = material.reference_density_kg_m3();
        let mut reference_density_kg_m3 = Vec::with_capacity(element_count);
        for element in 0..element_count {
            if element % 256 == 0 {
                checkpoint(cx)?;
            }
            reference_density_kg_m3.push(density);
        }
        checkpoint(cx)?;
        Ok(Self {
            source: SourceField::Uniform {
                material,
                before_thermal,
                reference,
            },
            reference_density_kg_m3,
        })
    }

    /// Keep one admitted source and one declared prethermal map per element.
    /// No constitutive law or reference strain is averaged across an interface.
    pub fn per_element(
        materials: Vec<FiniteSolidMaterial>,
        before_thermal: &[ReferenceDeformation],
        cx: &Cx<'_>,
    ) -> Result<Self, ConservativeThermalError> {
        checkpoint(cx)?;
        if materials.is_empty() || before_thermal.len() != materials.len() {
            return Err(ConservativeThermalError::Invalid(
                "one source and prethermal reference per original element required",
            ));
        }
        let mut laws = Vec::with_capacity(materials.len());
        let mut references = Vec::with_capacity(materials.len());
        let mut reference_density_kg_m3 = Vec::with_capacity(materials.len());
        for (element, (material, &reference)) in materials.iter().zip(before_thermal).enumerate() {
            checkpoint(cx)?;
            references.push(compose_reference(material, reference, element)?);
            laws.push(material.law().law().clone());
            reference_density_kg_m3.push(material.reference_density_kg_m3());
        }
        checkpoint(cx)?;
        Ok(Self {
            source: SourceField::PerElement {
                materials,
                laws,
                before_thermal: before_thermal.to_vec(),
                references,
            },
            reference_density_kg_m3,
        })
    }

    pub fn element_count(&self) -> usize {
        self.reference_density_kg_m3.len()
    }

    pub fn reference_density_kg_m3(&self) -> &[f64] {
        &self.reference_density_kg_m3
    }

    pub fn material_at(&self, element: usize) -> Option<&FiniteSolidMaterial> {
        if element >= self.element_count() {
            return None;
        }
        Some(match &self.source {
            SourceField::Uniform { material, .. } => material,
            SourceField::PerElement { materials, .. } => &materials[element],
        })
    }

    pub fn reference_before_thermal_at(&self, element: usize) -> Option<&ReferenceDeformation> {
        if element >= self.element_count() {
            return None;
        }
        Some(match &self.source {
            SourceField::Uniform { before_thermal, .. } => before_thermal,
            SourceField::PerElement { before_thermal, .. } => &before_thermal[element],
        })
    }

    pub fn reference_at(&self, element: usize) -> Option<&ReferenceDeformation> {
        if element >= self.element_count() {
            return None;
        }
        Some(match &self.source {
            SourceField::Uniform { reference, .. } => reference,
            SourceField::PerElement { references, .. } => &references[element],
        })
    }

    /// Prepare the actual conservative owner with the source field's original
    /// density and complete G. Input displacement, velocity and pressure must
    /// already satisfy the chosen native formulation; none is projected.
    #[allow(clippy::too_many_arguments)]
    pub fn prepare<'a>(
        &'a self,
        nodes_m: &'a [[f64; 3]],
        elements: &'a [[usize; 10]],
        displacement_m: &[[f64; 3]],
        velocity_mps: &[[f64; 3]],
        initial_pressure_pa: &[f64],
        controls: ConservativeThermalPreparation,
        cx: &Cx<'_>,
    ) -> Result<ConservativeThermalSolid<'a>, ConservativeThermalError> {
        checkpoint(cx)?;
        if elements.len() != self.element_count() {
            return Err(ConservativeThermalError::Invalid(
                "source field must cover every original element exactly",
            ));
        }
        let (material, references) = match &self.source {
            SourceField::Uniform {
                material,
                reference,
                ..
            } => (
                MixedLawField::Uniform(material.law().law()),
                ConservativeReferenceField::Uniform(reference),
            ),
            SourceField::PerElement {
                laws, references, ..
            } => (
                MixedLawField::PerElement(laws),
                ConservativeReferenceField::PerElement(references),
            ),
        };
        let definition = ConservativeTetDefinition {
            nodes_m,
            elements,
            material,
            reference_density_kg_m3: &self.reference_density_kg_m3,
            formulation: controls.formulation,
            budget: controls.budget,
            quadrature_points_per_dir: controls.quadrature_points_per_dir,
            maximum_quadrature_points: controls.maximum_quadrature_points,
            maximum_mass_assembly_bytes: controls.maximum_mass_assembly_bytes,
            material_settings: controls.material_settings,
        };
        let native = definition.prepare_with_reference_deformations(
            references,
            displacement_m,
            velocity_mps,
            initial_pressure_pa,
            cx,
        )?;
        checkpoint(cx)?;
        Ok(ConservativeThermalSolid {
            field: self,
            native,
        })
    }
}

/// Source binding and native accepted conservative trajectory share a lifetime.
#[derive(Debug)]
pub struct ConservativeThermalSolid<'a> {
    field: &'a ConservativeThermalField,
    native: ConservativeTetSolid<'a>,
}

#[derive(Debug)]
pub struct ConservativeThermalStep<'a> {
    pub solid: ConservativeThermalSolid<'a>,
    pub newton: Option<NewtonReport>,
    pub reactions_n: Vec<(usize, f64)>,
    pub energy: ConservativeTetEnergy,
    pub inertia: InertialEnergy,
    pub min_det_f: f64,
    pub pressure_residual_m3: Vec<f64>,
}

impl<'a> ConservativeThermalSolid<'a> {
    pub fn native(&self) -> &ConservativeTetSolid<'a> {
        &self.native
    }
    pub fn source_field(&self) -> &'a ConservativeThermalField {
        self.field
    }

    pub fn reference_acceleration_load_n(
        &self,
        acceleration_m_s2: &[[f64; 3]],
        cx: &Cx<'_>,
    ) -> Result<Vec<(usize, f64)>, ConservativeThermalError> {
        self.native
            .reference_acceleration_load_n(acceleration_m_s2, cx)
            .map_err(Into::into)
    }

    /// Advance the same held source/reference under native objective work,
    /// pressure, momentum and accumulated-energy gates. A refusal leaves this
    /// accepted source binding and every native state component unchanged.
    #[allow(clippy::too_many_arguments)]
    pub fn advance(
        &self,
        dt_s: f64,
        prescribed_m: &[(usize, f64)],
        nodal_forces_n: &[(usize, f64)],
        body_force_n_m3: [f64; 3],
        settings: MixedUpSettings,
        maximum_work_error_j: f64,
        cx: &Cx<'_>,
    ) -> Result<ConservativeThermalStep<'a>, ConservativeThermalError> {
        let step = self.native.advance(
            dt_s,
            prescribed_m,
            nodal_forces_n,
            body_force_n_m3,
            settings,
            maximum_work_error_j,
            cx,
        )?;
        Ok(ConservativeThermalStep {
            solid: Self {
                field: self.field,
                native: step.solid,
            },
            newton: step.newton,
            reactions_n: step.reactions_n,
            energy: step.energy,
            inertia: step.inertia,
            min_det_f: step.min_det_f,
            pressure_residual_m3: step.pressure_residual_m3,
        })
    }
}
