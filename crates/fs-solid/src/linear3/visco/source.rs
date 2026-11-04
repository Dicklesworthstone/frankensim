//! Source coefficient changes on the original sparse finite-element histories.
//!
//! An endpoint refresh is an externally prescribed parameter operation at fixed
//! geometry, reference mass, displacement and physical Maxwell spring strains.
//! Its signed storage change is parameter work, never dashpot heat. Coefficients
//! are then held during each existing exact-ramp material/time step. Refining
//! that splitting is necessary for continuously varying material parameters;
//! neither reciprocal thermal feedback nor an unsupported source path is inferred.

use super::*;
use crate::linear3::QuadraticThermalStrainField;
use fs_matdb::{MaterialCard, QueryPoint};
use fs_material::visco_tensor::source::{ResolvedTensorMaterial, TensorSourceError};
use std::sync::Arc;

/// One original material source and its fixed material-to-mesh frame.
#[derive(Clone, Debug)]
pub struct ViscoMaterialSource {
    pub material: ResolvedTensorMaterial,
    pub material_to_mesh: [[f64; 3]; 3],
}

/// Preserve the original failing source or native finite-element owner.
#[derive(Debug)]
pub enum ViscoSourceError {
    Native(ViscoElementError),
    Source(TensorSourceError),
    Invalid(&'static str),
}
impl core::fmt::Display for ViscoSourceError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "native viscoelastic material continuation: {self:?}")
    }
}
impl std::error::Error for ViscoSourceError {}
impl From<ViscoElementError> for ViscoSourceError {
    fn from(error: ViscoElementError) -> Self {
        Self::Native(error)
    }
}
impl From<TensorSourceError> for ViscoSourceError {
    fn from(error: TensorSourceError) -> Self {
        Self::Source(error)
    }
}
impl From<ViscoError> for ViscoSourceError {
    fn from(error: ViscoError) -> Self {
        Self::Native(error.into())
    }
}

pub(super) struct SourceBinding {
    // Only at_query_points can share this origin. Its geometry, B/B-bar columns,
    // reference inventory, source selectors, branch order and frames cannot change.
    reference: Arc<[ViscoMaterialSource]>,
    current: Arc<[ViscoMaterialSource]>,
}

fn invalid(what: &'static str) -> ViscoSourceError {
    ViscoSourceError::Invalid(what)
}
fn finite(value: f64) -> Result<f64, ViscoSourceError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(ViscoElementError::NonFiniteResult.into())
    }
}

impl ViscoElements {
    /// Bind the exact sources that assembled this original solid, once.
    /// Every tensor, branch, rate, reference density and material identity must
    /// match; labels alone cannot establish a continuation of live memories.
    pub fn with_material_sources(
        mut self,
        sources: Vec<ViscoMaterialSource>,
        cx: &Cx<'_>,
    ) -> Result<Self, ViscoSourceError> {
        checkpoint(cx)?;
        if self.sources.is_some() || sources.len() != self.elements.len() {
            return Err(invalid(
                "one original source per element and an unbound solid required",
            ));
        }
        for (element, source) in self.elements.iter().zip(&sources) {
            checkpoint(cx)?;
            let supplied =
                TetViscoMaterial::from_source(&source.material, source.material_to_mesh)?;
            if supplied.law != element.law
                || supplied.density_kg_m3.to_bits() != element.density_kg_m3.to_bits()
                || supplied.material_state_identity != element.material_state_identity
            {
                return Err(invalid(
                    "source law, frame, identity and density must match the original solid",
                ));
            }
        }
        let reference: Arc<[ViscoMaterialSource]> = sources.into();
        self.sources = Some(SourceBinding {
            current: reference.clone(),
            reference,
        });
        checkpoint(cx)?;
        Ok(self)
    }

    /// Current pinned coefficient queries, in original element order.
    pub fn material_sources(&self) -> Option<&[ViscoMaterialSource]> {
        self.sources
            .as_ref()
            .map(|binding| binding.current.as_ref())
    }

    /// The original source states that supplied the retained reference mass.
    pub fn reference_material_sources(&self) -> Option<&[ViscoMaterialSource]> {
        self.sources
            .as_ref()
            .map(|binding| binding.reference.as_ref())
    }

    /// Resolve prescribed endpoints from the SAME original source selections.
    ///
    /// The existing native source owner retains branch order, schema and Kelvin
    /// selections. Material axes and the original B/B-bar quadrature are fixed.
    /// Rebuild only equilibrium stiffness; consistent mass, reference density,
    /// constraints and quadrature weights remain the original inventory. A
    /// current source density is retained as metadata and never becomes mass.
    /// Each element has one actual coefficient-query temperature. The immutable
    /// reference source remains the inventory datum, not the current thermal
    /// state; continuation checks the new query against all held thermal samples.
    ///
    /// Each endpoint must be supported. This API does not certify an intervening
    /// path, interpolate missing data, transfer histories or solve heat balance.
    pub fn at_query_points(
        &self,
        cards: &[&MaterialCard],
        points: &[QueryPoint],
        cx: &Cx<'_>,
    ) -> Result<Self, ViscoSourceError> {
        checkpoint(cx)?;
        let binding = self
            .sources
            .as_ref()
            .ok_or_else(|| invalid("material refresh requires original source bindings"))?;
        if cards.len() != self.elements.len() || points.len() != self.elements.len() {
            return Err(invalid(
                "complete original cards and physical queries per element required",
            ));
        }
        let mut sources = Vec::with_capacity(self.elements.len());
        let mut elements = Vec::with_capacity(self.elements.len());
        let mut tensors = Vec::with_capacity(self.elements.len());
        for (((old, source), &card), point) in self
            .elements
            .iter()
            .zip(binding.current.iter())
            .zip(cards)
            .zip(points)
        {
            checkpoint(cx)?;
            let material = source.material.at_query_point(card, point)?;
            let law = material.rotated_law(source.material_to_mesh)?;
            if law.branches().len() != old.law.branches().len()
                || material.properties() != source.material.properties()
            {
                return Err(invalid(
                    "source refresh must retain original branch identities and ordering",
                ));
            }
            tensors.push(*law.c_eq());
            let mut quadrature = Vec::with_capacity(old.points.len());
            for point in &old.points {
                checkpoint(cx)?;
                quadrature.push(IntegrationPoint {
                    b: point.b.clone(),
                    compatible_gradients: point.compatible_gradients.clone(),
                    volume_m3: point.volume_m3,
                });
            }
            elements.push(Element {
                dofs: old.dofs.clone(),
                law,
                density_kg_m3: old.density_kg_m3,
                material_state_identity: material.state().identity().0,
                points: quadrature,
            });
            sources.push(ViscoMaterialSource {
                material,
                material_to_mesh: source.material_to_mesh,
            });
        }
        let (_, stiffness) = super::quasistatic::assemble_tangents(self, &tensors, cx)?;
        let original = &self.equilibrium;
        let equilibrium = TetElasticAssembly {
            stiffness,
            mass: original.mass.clone(),
            free_dofs: original.free_dofs.clone(),
            total_mass_kg: original.total_mass_kg,
            minimum_scaled_jacobian: original.minimum_scaled_jacobian,
            element_volumes_m3: original.element_volumes_m3.clone(),
        };
        checkpoint(cx)?;
        Ok(Self {
            elements,
            equilibrium,
            full_dofs: self.full_dofs,
            internal_variables: self.internal_variables,
            static_constraint_error: self.static_constraint_error.clone(),
            sources: Some(SourceBinding {
                reference: binding.reference.clone(),
                current: sources.into(),
            }),
        })
    }
}

/// Complete continuation at the SAME physical instant. The input history is
/// untouched on both success and refusal; the caller publishes this candidate.
pub struct ViscoMaterialChange<'a> {
    pub history: ViscoElementState<'a>,
    /// Signed change of stored material energy at fixed physical strains [J].
    /// This is external parameter work; no time passes and no heat is generated.
    pub parameter_work_j: f64,
}

impl ViscoElementState<'_> {
    /// Cumulative signed coefficient work since the original relaxed state [J].
    /// Together with energy().work_j it supplies storage change and dashpot heat.
    pub const fn parameter_work_j(&self) -> f64 {
        self.parameter_work_j
    }

    /// Continue the original physical spring strains with new source coefficients.
    /// Displacement, branch memories, accumulated strain work and heat survive.
    /// An accepted thermal history requires explicit replacement receipts through
    /// rebind_materials_with_thermal_strain; it cannot be silently relabeled.
    pub fn rebind_materials<'b>(
        &self,
        next: &'b ViscoElements,
        cx: &Cx<'_>,
    ) -> Result<ViscoMaterialChange<'b>, ViscoSourceError> {
        self.rebind_materials_impl(next, None, cx)
    }

    /// Rebind coefficients while holding the entire accepted eigenstrain field.
    /// Every replacement receipt must retain the same raw free strain, current
    /// and reference temperatures, expansion-law identity and strain ceiling.
    /// Only its elastic-material identity changes to the new pinned source.
    /// Changes of free strain remain finite-duration inputs to the native history
    /// owner; this instantaneous parameter operation does not bypass Kelvin loss.
    /// The new element coefficient query must carry canonical T equal to EVERY
    /// held sample's current temperature. Within-element temperature variation
    /// requires a more resolved coefficient model; it is never averaged here.
    /// A thermal split first advances the actual free-strain field under held
    /// coefficients, then applies this fixed-strain coefficient work operation.
    pub fn rebind_materials_with_thermal_strain<'b, 't>(
        &self,
        next: &'b ViscoElements,
        thermal: impl Into<QuadraticThermalStrainField<'t>>,
        cx: &Cx<'_>,
    ) -> Result<ViscoMaterialChange<'b>, ViscoSourceError> {
        self.rebind_materials_impl(next, Some(thermal.into()), cx)
    }

    pub(super) fn rebind_materials_impl<'b>(
        &self,
        next: &'b ViscoElements,
        thermal: Option<QuadraticThermalStrainField<'_>>,
        cx: &Cx<'_>,
    ) -> Result<ViscoMaterialChange<'b>, ViscoSourceError> {
        checkpoint(cx)?;
        let old = self
            .system
            .sources
            .as_ref()
            .ok_or_else(|| invalid("live continuation requires original source bindings"))?;
        let new = next
            .sources
            .as_ref()
            .ok_or_else(|| invalid("replacement solid lacks original source bindings"))?;
        if !Arc::ptr_eq(&old.reference, &new.reference) {
            return Err(invalid(
                "continuation must retain the exact source solid, reference mass and frames",
            ));
        }
        let thermal_strains = match (&self.thermal_strains, thermal) {
            (None, None) => None,
            (Some(old), Some(field)) => Some(old.rebound(next, field, cx)?),
            _ => {
                return Err(invalid(
                    "accepted thermal history requires complete held-strain replacement receipts",
                ));
            }
        };
        next.validate_displacement(&self.displacement_m, || cx.checkpoint().is_err())?;
        let mut histories = Vec::with_capacity(self.histories.len());
        let mut stored_j = 0.0;
        let mut parameter_work_j = 0.0;
        for (element_index, ((previous, element), accepted)) in self
            .system
            .elements
            .iter()
            .zip(&next.elements)
            .zip(&self.histories)
            .enumerate()
        {
            checkpoint(cx)?;
            // Rates and Kelvin viscosity do not store energy. Keep exact zero
            // parameter work when the physical spring tensors are unchanged.
            let same_storage = previous.law.c_eq() == element.law.c_eq()
                && previous
                    .law
                    .branches()
                    .iter()
                    .zip(element.law.branches())
                    .all(|(a, b)| a.c_pa == b.c_pa);
            let source_temperature = if thermal_strains.is_some() {
                Some(new.current[element_index].material.state().state_coordinate("T")
                    .ok_or_else(|| invalid("thermal continuation requires the original canonical T coefficient query"))?)
            } else {
                None
            };
            let mut points = Vec::with_capacity(accepted.len());
            for (point_index, (point, old_state)) in element.points.iter().zip(accepted).enumerate()
            {
                checkpoint(cx)?;
                element.law.validate_mandel_strain(&old_state.strain)?;
                if let Some(field) = &thermal_strains {
                    if field
                        .at(element_index, point_index)
                        .temperature_k()
                        .to_bits()
                        != source_temperature
                            .expect("thermal query admitted")
                            .to_bits()
                    {
                        return Err(invalid(
                            "element coefficient temperature must match every held thermal quadrature state",
                        ));
                    }
                    let total = point.strain(&element.dofs, &self.displacement_m);
                    let raw = point.compatible_strain(&element.dofs, &self.displacement_m, total);
                    let free = field.at(element_index, point_index).free_strain_mandel();
                    element
                        .law
                        .validate_mandel_strain(&core::array::from_fn(|i| raw[i] - free[i]))?;
                }
                let mut state = old_state.clone();
                if previous.law != element.law {
                    // Native zero-duration evaluation refreshes the storage law
                    // without relaxing a branch. No force observation is emitted;
                    // subsequent finite-duration steps own their actual Kelvin rate.
                    element.law.step(&mut state, old_state.strain, 0.0)?;
                    state.q.clone_from(&old_state.q);
                    state.strain = old_state.strain;
                    state.work = old_state.work;
                    state.dissipated = old_state.dissipated;
                    if same_storage {
                        state.stored = old_state.stored;
                    }
                }
                let work = finite(point.volume_m3 * finite(state.stored - old_state.stored)?)?;
                parameter_work_j = finite(parameter_work_j + work)?;
                stored_j = finite(stored_j + finite(point.volume_m3 * state.stored)?)?;
                points.push(state);
            }
            histories.push(points);
        }
        let cumulative = finite(self.parameter_work_j + parameter_work_j)?;
        checkpoint(cx)?;
        Ok(ViscoMaterialChange {
            history: ViscoElementState {
                system: next,
                histories,
                energy: ViscoEnergy {
                    stored_j,
                    ..self.energy
                },
                displacement_m: self.displacement_m.clone(),
                thermal_strains,
                eigenstrain_work_j: self.eigenstrain_work_j,
                parameter_work_j: cumulative,
            },
            parameter_work_j,
        })
    }
}
