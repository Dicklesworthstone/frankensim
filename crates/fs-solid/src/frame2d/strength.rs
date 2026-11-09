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
//! Native Kelvin--Voigt members additionally use the accepted nodal velocity
//! and original reference interpolation. A layer's normal stress is
//! `E_i epsilon + eta_ai epsilon_dot +/- r_i |E_i kappa + eta_bi kappa_dot|`.
//! Signed elastic and viscous bending contributions combine before taking an
//! absolute value. Their affine spatial variation still has an exact end
//! envelope. Distinct layer viscosities are never redistributed through E/EA.
//!
//! These are observations of accepted source material members, including
//! static and Newmark endpoints. They do not model plastic redistribution,
//! rupture, fatigue, shear failure, or buckling. Timoshenko members retain the
//! same axial/bending equilibrium envelope, but their shear stresses are not
//! included in this explicitly normal-stress-only criterion. Authored Rayleigh
//! damping is not a source material normal stress. Native source viscosity is
//! linear-reference Kelvin--Voigt, with no inferred thermal or history stress.

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
    /// Signed section axial force, including original material viscosity.
    pub axial_force_n: f64,
    /// Signed section bending moment, including original material viscosity;
    /// distinct from the basic elastic nodal moment.
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
    /// Original reference geometry or accepted velocity cannot produce a rate.
    Kinematics {
        element: usize,
        reason: String,
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
    /// section law. Native material-viscous members combine their original
    /// axial and bending viscosities with the actual accepted velocity field,
    /// retaining each layer's distinct signed rate stresses. Shear-deformable
    /// members use their original Timoshenko curvature interpolation. A
    /// temperature refresh needs a refreshed capacity, even when its scalar
    /// strength is unchanged; no thermal or hereditary stress is inferred.
    ///
    /// `maximum_layer_samples` bounds twice the total layer count before output
    /// allocation. No tangent, constitutive update or extra solve is performed.
    /// Cancellation and all refusals leave the accepted frame unchanged.
    #[allow(clippy::too_many_lines)]
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
            if member.sections.iter().any(|section| {
                matches!(section, super::SectionModel::ThermalElastic(_))
            }) {
                return Err(FrameTensileError::Invalid(
                    "thermal frame layer stresses require thermal section recovery; resultant-only strength omits thermal mismatch",
                ));
            }
            if member.material_section.as_ref() != Some(capacity.section()) {
                return Err(FrameTensileError::SourceMismatch { element });
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
            let [elastic_axial_force, moment_i, moment_j] = member.commit.q;
            let viscous = if let Some(source) = &member.viscous_section {
                // Reuse the exact work-conjugate columns used by native C v.
                // They include rigid translation/rotation cancellation and
                // refuse unsupported corotational viscous geometry.
                let kinematic_error = |error: super::SolidError| FrameTensileError::Kinematics {
                    element,
                    reason: error.to_string(),
                };
                let (dofs, columns, length) = self
                    .reference_member_coordinates(element)
                    .map_err(kinematic_error)?;
                let coordinate_rates = columns.map(|column| {
                    column
                        .into_iter()
                        .zip(dofs)
                        .map(|(weight, dof)| weight * self.vel[dof])
                        .sum::<f64>()
                });
                let kinematics = source
                    .reference_kinematics(length)
                    .map_err(kinematic_error)?;
                let measures = kinematics
                    .strain_measures(coordinate_rates)
                    .map_err(kinematic_error)?;
                let symmetric =
                    3.0 * kinematics.bending_fraction() * (coordinate_rates[1] / length / length);
                let antisymmetric = coordinate_rates[2] / length / length;
                let curvature_rates = [-(symmetric + antisymmetric), symmetric - antisymmetric];
                let prepared = capacity
                    .clone()
                    .with_kelvin_voigt(source.axial(), source.bending())
                    .map_err(|source| FrameTensileError::Material { element, source })?;
                Some((prepared, measures.axial_strain, curvature_rates))
            } else {
                None
            };
            for (end, elastic_moment) in [-moment_i, moment_j].into_iter().enumerate() {
                checkpoint(cx)?;
                let (assessment, axial_force_n, bending_moment_n_m) = match &viscous {
                    Some((prepared, axial_rate, curvature_rates)) => {
                        let section = capacity.section();
                        let axial_strain = elastic_axial_force / section.axial_stiffness_n();
                        let curvature = elastic_moment / section.bending_stiffness_n_m2();
                        if (elastic_axial_force != 0.0 && axial_strain == 0.0)
                            || (elastic_moment != 0.0 && curvature == 0.0)
                        {
                            return Err(FrameTensileError::Invalid(
                                "elastic material section strain is unrepresentable",
                            ));
                        }
                        let assessment = prepared.assess_strain_rates(
                            axial_strain,
                            *axial_rate,
                            [curvature, 0.0],
                            [curvature_rates[end], 0.0],
                        );
                        (
                            assessment,
                            prepared
                                .axial_viscous_rigidity_n_s()
                                .mul_add(*axial_rate, elastic_axial_force),
                            prepared
                                .bending_viscous_stiffness_n_m2_s()
                                .mul_add(curvature_rates[end], elastic_moment),
                        )
                    }
                    None => (
                        capacity.assess_axial_bending(elastic_axial_force, [elastic_moment, 0.0]),
                        elastic_axial_force,
                        elastic_moment,
                    ),
                };
                let assessment =
                    assessment.map_err(|source| FrameTensileError::Material { element, source })?;
                if !axial_force_n.is_finite() || !bending_moment_n_m.is_finite() {
                    return Err(FrameTensileError::Invalid(
                        "total material section force is unrepresentable",
                    ));
                }
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
