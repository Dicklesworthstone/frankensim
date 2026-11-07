//! Tensile strength of the original force-based material frame.
//!
//! Every prismatic member has constant signed axial force and affine bending
//! moment `M(xi) = (xi - 1) M_i + xi M_j`. For each original concentric layer,
//! its greatest tensile normal stress is therefore attained at a member end:
//! `E_i (N/EA + r_i |M|/EI)`. Checking both ends retains this complete spatial
//! envelope without sampling an arbitrary number of interior sections or
//! averaging opposite bending moments. Large rigid rotations change neither
//! the local section forces nor this scalar material-axis comparison.
//!
//! These are observations of accepted elastic material members, including
//! static and Newmark endpoints. They do not model plastic redistribution,
//! rupture, fatigue, shear failure, or buckling. Timoshenko members retain the
//! same axial/bending equilibrium envelope, but their shear stresses are not
//! included in this explicitly normal-stress-only criterion.

use fs_exec::Cx;
use fs_material::section::strength::{AxialTensileCapacity, ConcentricTensileAssessment};
use fs_material::state_point::MaterialStatePointError;

use super::Frame2d;

/// One physical member end with all of its original layer strengths.
#[derive(Debug, Clone, PartialEq)]
pub struct FrameTensileSample<'a> {
    pub element: usize,
    /// Zero is the original start node, one the original end node.
    pub end: usize,
    pub axial_force_n: f64,
    /// Signed section bending moment, distinct from the basic nodal moment.
    pub bending_moment_n_m: f64,
    pub capacity: &'a AxialTensileCapacity,
    pub assessment: ConcentricTensileAssessment,
}

/// Immutable observations in original element order, start before end.
#[derive(Debug, Clone, PartialEq)]
pub struct FrameTensileAssessment<'a> {
    samples: Vec<FrameTensileSample<'a>>,
    governing_sample: Option<usize>,
}

impl<'a> FrameTensileAssessment<'a> {
    #[must_use]
    pub fn samples(&self) -> &[FrameTensileSample<'a>] {
        &self.samples
    }

    /// Largest tensile utilization; the first sample wins exact ties.
    #[must_use]
    pub fn governing_sample(&self) -> Option<&FrameTensileSample<'a>> {
        self.governing_sample.map(|index| &self.samples[index])
    }

    #[must_use]
    pub fn reaches_limit(&self) -> bool {
        self.governing_sample()
            .is_some_and(|sample| sample.assessment.utilization() >= 1.0)
    }
}

/// Source, storage, numerical or cancellation refusal; no partial assessment.
#[derive(Debug, Clone, PartialEq)]
pub enum FrameTensileError {
    Invalid(&'static str),
    SourceMismatch {
        element: usize,
    },
    /// Material viscosity needs its own per-layer rate stresses, not E times strain.
    ViscousMember {
        element: usize,
    },
    Budget,
    Allocation,
    Material {
        element: usize,
        source: MaterialStatePointError,
    },
    Cancelled,
}

impl core::fmt::Display for FrameTensileError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(
            f,
            "native material-frame tensile assessment refused: {self:?}"
        )
    }
}
impl std::error::Error for FrameTensileError {}

impl Frame2d {
    /// Compare accepted local axial/bending stresses with each layer's limit.
    ///
    /// Supply exactly one capacity for every original material member. Each
    /// capacity must retain that member's complete elastic section, geometry,
    /// material cards and physical state. A raw-coefficient or plastic fiber
    /// member refuses: its stress cannot be reconstructed from this elastic
    /// section law. Material-viscous members also refuse because the retained
    /// elastic force omits their per-layer rate stresses. A temperature refresh
    /// needs a refreshed capacity, even when its scalar strength is unchanged.
    ///
    /// `maximum_layer_samples` bounds twice the total layer count before output
    /// allocation. No tangent, constitutive update or extra solve is performed.
    /// Cancellation and all refusals leave the accepted frame unchanged.
    pub fn assess_tensile_strength<'a>(
        &self,
        capacities: &'a [AxialTensileCapacity],
        maximum_layer_samples: usize,
        cx: &Cx<'_>,
    ) -> Result<FrameTensileAssessment<'a>, FrameTensileError> {
        checkpoint(cx)?;
        if capacities.len() != self.elements.len() {
            return Err(FrameTensileError::Invalid(
                "exactly one tensile capacity is required per original material member",
            ));
        }
        let mut layers = 0_usize;
        for (element, (member, capacity)) in self.elements.iter().zip(capacities).enumerate() {
            checkpoint(cx)?;
            if member.material_section.as_ref() != Some(capacity.section()) {
                return Err(FrameTensileError::SourceMismatch { element });
            }
            if member.viscous_section.is_some() {
                return Err(FrameTensileError::ViscousMember { element });
            }
            layers = capacity
                .section()
                .layers()
                .len()
                .checked_mul(2)
                .and_then(|count| layers.checked_add(count))
                .filter(|&count| count <= maximum_layer_samples)
                .ok_or(FrameTensileError::Budget)?;
        }
        let sample_count = capacities
            .len()
            .checked_mul(2)
            .ok_or(FrameTensileError::Budget)?;
        let mut samples = Vec::new();
        samples
            .try_reserve_exact(sample_count)
            .map_err(|_| FrameTensileError::Allocation)?;
        let mut governing_sample = None;
        let mut largest = -1.0;
        for (element, (member, capacity)) in self.elements.iter().zip(capacities).enumerate() {
            checkpoint(cx)?;
            let [axial_force_n, moment_i, moment_j] = member.commit.q;
            for (end, bending_moment_n_m) in [-moment_i, moment_j].into_iter().enumerate() {
                let assessment = capacity
                    .assess_axial_bending(axial_force_n, [bending_moment_n_m, 0.0])
                    .map_err(|source| FrameTensileError::Material { element, source })?;
                if assessment.utilization() > largest {
                    largest = assessment.utilization();
                    governing_sample = Some(samples.len());
                }
                samples.push(FrameTensileSample {
                    element,
                    end,
                    axial_force_n,
                    bending_moment_n_m,
                    capacity,
                    assessment,
                });
            }
        }
        checkpoint(cx)?;
        Ok(FrameTensileAssessment {
            samples,
            governing_sample,
        })
    }
}

fn checkpoint(cx: &Cx<'_>) -> Result<(), FrameTensileError> {
    cx.checkpoint().map_err(|_| FrameTensileError::Cancelled)
}
