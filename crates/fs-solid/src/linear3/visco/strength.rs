//! Source-bound strength on the original tensor-viscoelastic stress histories.
//!
//! Every retained integration point supplies its actual endpoint stress,
//! including Maxwell memory, Kelvin stress from the last accepted strain ramp,
//! and imposed thermal eigenstrain. Recomputing equilibrium elasticity from
//! nodal displacement would omit those physical contributions. The same
//! material axes used by the native constitutive assembly define strength.
//!
//! This observer applies the shared directional tensile or 3-D Tsai-Wu
//! criterion. It does not advance time, relax a history, evolve damage, or
//! certify peaks between integration points or accepted temporal endpoints.

use super::{ViscoElementState, ViscoElements};
use crate::strength::orthotropic::{OrthotropicStrengthCriterion, OrthotropicStrengthObservation};
use crate::strength::{TensileAssessmentError, admit_same_state, checkpoint};
use fs_blake3::ContentHash;

struct Binding<C> {
    criterion: C,
    world_to_material: [[f64; 3]; 3],
    temperature_k: Option<f64>,
    material_state_identity: ContentHash,
}

/// Criteria prepared against one exact native geometry/material/frame owner.
/// The source owner remains borrowed, so its full coefficient receipts and
/// original integration-point ordering cannot be replaced behind this binding.
pub struct PreparedViscoStrength<'a, C: OrthotropicStrengthCriterion> {
    system: &'a ViscoElements,
    bindings: Vec<Binding<C>>,
    samples: usize,
}

/// One actual native endpoint stress in its original material axes.
#[derive(Clone, Debug, PartialEq)]
pub struct ViscoStrengthSample<A: OrthotropicStrengthObservation> {
    pub element: usize,
    pub quadrature_point: usize,
    /// Current native material coefficient query, before frame rotation.
    pub material_state_identity: ContentHash,
    /// Original directional strength claims and complete physical query.
    pub limit_source_identity: ContentHash,
    pub assessment: A,
}

/// Complete unaveraged observation, in native element/point order.
#[derive(Clone, Debug, PartialEq)]
pub struct ViscoStrengthAssessment<A: OrthotropicStrengthObservation> {
    samples: Vec<ViscoStrengthSample<A>>,
    governing: Option<usize>,
}

impl<A: OrthotropicStrengthObservation> ViscoStrengthAssessment<A> {
    pub fn samples(&self) -> &[ViscoStrengthSample<A>] {
        &self.samples
    }

    /// First native point with the largest utilization; ties are deterministic.
    pub fn governing_sample(&self) -> Option<&ViscoStrengthSample<A>> {
        self.governing.map(|point| &self.samples[point])
    }

    pub fn reaches_limit(&self) -> bool {
        self.governing_sample()
            .is_some_and(|point| point.assessment.reaches_limit())
    }
}

impl ViscoElements {
    /// Bind one explicit criterion per original source-bearing element.
    ///
    /// Call `with_material_sources` when constructing the mechanics. Every
    /// criterion must select the SAME immutable material card and complete
    /// typed physical query as that element's current coefficients. Rotation
    /// is the original material-to-mesh frame, never a principal-stress frame.
    ///
    /// A source refresh creates a new native owner and requires a new prepared
    /// strength binding with correspondingly refreshed criteria. Original
    /// reference mass is not compared to a current source density or changed.
    pub fn prepare_orthotropic_strength<C: OrthotropicStrengthCriterion>(
        &self,
        criteria: Vec<C>,
        cx: &fs_exec::Cx<'_>,
    ) -> Result<PreparedViscoStrength<'_, C>, TensileAssessmentError> {
        checkpoint(cx)?;
        let sources = self
            .material_sources()
            .ok_or(TensileAssessmentError::Invalid(
                "viscoelastic strength requires the original material source bindings",
            ))?;
        if criteria.len() != self.elements.len() {
            return Err(TensileAssessmentError::Invalid(
                "one strength criterion per original viscoelastic element required",
            ));
        }
        let mut bindings = Vec::new();
        bindings
            .try_reserve_exact(criteria.len())
            .map_err(|_| TensileAssessmentError::Budget)?;
        let mut samples = 0_usize;
        for ((criterion, source), element) in criteria.into_iter().zip(sources).zip(&self.elements)
        {
            checkpoint(cx)?;
            admit_same_state(source.material.state(), criterion.resolved())?;
            samples = samples
                .checked_add(element.points.len())
                .ok_or(TensileAssessmentError::Budget)?;
            bindings.push(Binding {
                temperature_k: source.material.state().state_coordinate("T"),
                material_state_identity: source.material.state().identity(),
                criterion,
                world_to_material: core::array::from_fn(|i| {
                    core::array::from_fn(|j| source.material_to_mesh[j][i])
                }),
            });
        }
        checkpoint(cx)?;
        Ok(PreparedViscoStrength {
            system: self,
            bindings,
            samples,
        })
    }
}

impl<'a, C: OrthotropicStrengthCriterion> PreparedViscoStrength<'a, C> {
    /// Native coefficient/geometry owner and its complete source receipts.
    pub const fn system(&self) -> &'a ViscoElements {
        self.system
    }

    /// Selected original source criterion for one element, if present.
    pub fn criterion(&self, element: usize) -> Option<&C> {
        self.bindings.get(element).map(|binding| &binding.criterion)
    }

    /// Observe every actual accepted constitutive stress without a new solve.
    ///
    /// The state must belong to this exact source/geometry owner; equal
    /// displacement or stiffness is insufficient to rebind live histories.
    /// The original P1, quadratic, hexahedral and B-bar quadrature semantics
    /// are preserved. `maximum_samples` bounds work and output allocation
    /// before observation, with cancellation per native point and before return.
    ///
    /// Only the result vector allocates. No source lookup, tangent, operator,
    /// material update, strain reconstruction or intermediate stress field is
    /// built. Kelvin stress is the existing accepted ramp-endpoint stress; no
    /// instantaneous rate or unsampled time evolution is inferred. If imposed
    /// thermal strains are present, each original sample temperature must match
    /// the current coefficient/strength query. A stale-temperature strength
    /// cannot silently assess a prescribed thermal path with frozen coefficients.
    /// After a law-changing Kelvin coefficient refresh, the next successful
    /// positive-duration native step must supply a complete stress observation.
    pub fn assess(
        &self,
        state: &ViscoElementState<'_>,
        maximum_samples: usize,
        cx: &fs_exec::Cx<'_>,
    ) -> Result<ViscoStrengthAssessment<C::Assessment>, TensileAssessmentError> {
        checkpoint(cx)?;
        if !core::ptr::eq(self.system, state.system()) {
            return Err(TensileAssessmentError::Invalid(
                "strength and constitutive histories require the same native viscoelastic owner",
            ));
        }
        if !state.endpoint_stress_observed {
            return Err(TensileAssessmentError::Invalid(
                "Kelvin coefficient refresh requires a new native endpoint stress observation",
            ));
        }
        if self.samples > maximum_samples {
            return Err(TensileAssessmentError::Budget);
        }
        let mut samples = Vec::new();
        samples
            .try_reserve_exact(self.samples)
            .map_err(|_| TensileAssessmentError::Budget)?;
        let mut governing = None;
        let mut largest = -1.0;
        for (element, (histories, binding)) in
            state.histories.iter().zip(&self.bindings).enumerate()
        {
            for (quadrature_point, history) in histories.iter().enumerate() {
                checkpoint(cx)?;
                if let Some(thermal) = &state.thermal_strains {
                    let temperature = thermal.at(element, quadrature_point).temperature_k();
                    if !binding
                        .temperature_k
                        .is_some_and(|t| t.to_bits() == temperature.to_bits())
                    {
                        return Err(TensileAssessmentError::Invalid(
                            "strength temperature must match each native thermal sample",
                        ));
                    }
                }
                // Native histories use work-conjugate Mandel coordinates.
                // Convert shear exactly once before the physical tensor
                // rotation and shared strength criterion.
                let material_stress = fs_material::tensor::rotate(
                    &fs_material::tensor::from_mandel(&history.stress),
                    &binding.world_to_material,
                );
                let assessment = binding
                    .criterion
                    .assess(&material_stress)
                    .map_err(TensileAssessmentError::Material)?;
                if assessment.utilization() > largest {
                    largest = assessment.utilization();
                    governing = Some(samples.len());
                }
                samples.push(ViscoStrengthSample {
                    element,
                    quadrature_point,
                    material_state_identity: binding.material_state_identity,
                    limit_source_identity: binding.criterion.resolved().identity(),
                    assessment,
                });
            }
        }
        checkpoint(cx)?;
        Ok(ViscoStrengthAssessment { samples, governing })
    }
}
