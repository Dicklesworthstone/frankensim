//! Principal tension from native axial force, bending and circular torsion.
//!
//! This opt-in observation combines the original segment's force and all three
//! moment resultants with its exact sourced shear section. The existing axial-
//! normal-stress assessment keeps its original meaning. Temperature mismatch
//! retains the original thermal section and strain bound. Pin-ended links stay
//! axial, at their actual force per strand; no shear law is invented for them.
//!
//! Saint-Venant torsion is admitted only for the existing bonded concentric
//! circular sections. Transverse force is NOT a local shear-stress field: the
//! actual omitted force is reported, and any nonzero component marks incomplete
//! transverse-shear assessment. Even roundoff-size forces retain that marker;
//! no physical insignificance tolerance is inferred. This is a selected local
//! maximum-principal-tension criterion, not a ductile yield or fracture law.

use fs_exec::CancelGate;
use fs_material::section::ConcentricShearSection;
use fs_material::section::strength::AxialTensileCapacity;
use fs_material::section::strength::combined::ConcentricPrincipalTensileAssessment;

use super::{
    RodTensileMember, RodTensileSample, SpatialMaterialRod, SpatialRodError, checked_layers,
    checkpoint,
};

/// One native rod segment assessed in axial force, bending and torsion.
#[derive(Clone, Debug, PartialEq)]
pub struct RodPrincipalTensileSample<'a> {
    segment: usize,
    capacity: &'a AxialTensileCapacity,
    shear_section: &'a ConcentricShearSection,
    total_axial_strain: f64,
    force_body_n: [f64; 3],
    moment_body_n_m: [f64; 3],
    assessment: ConcentricPrincipalTensileAssessment,
}

impl<'a> RodPrincipalTensileSample<'a> {
    /// Original material-interval index.
    #[must_use]
    pub const fn segment(&self) -> usize {
        self.segment
    }
    /// Original elastic section and explicitly selected per-layer strengths.
    #[must_use]
    pub const fn capacity(&self) -> &'a AxialTensileCapacity {
        self.capacity
    }
    /// Exact native shear-section source, including its selected shear law.
    #[must_use]
    pub const fn shear_section(&self) -> &'a ConcentricShearSection {
        self.shear_section
    }
    /// Native centerline strain before subtracting thermal eigenstrain.
    #[must_use]
    pub const fn total_axial_strain(&self) -> f64 {
        self.total_axial_strain
    }
    /// Signed native body-frame force [N]: axial, transverse-y, transverse-z.
    #[must_use]
    pub const fn force_body_n(&self) -> [f64; 3] {
        self.force_body_n
    }
    /// Signed native body-frame moment [N m]: torsion, bending-y, bending-z.
    #[must_use]
    pub const fn moment_body_n_m(&self) -> [f64; 3] {
        self.moment_body_n_m
    }
    /// Transverse forces [N] whose local shear stress is absent from this model.
    #[must_use]
    pub const fn unassessed_transverse_force_n(&self) -> [f64; 2] {
        [self.force_body_n[1], self.force_body_n[2]]
    }
    /// Whether any native transverse force is nonzero, without a hidden cutoff.
    #[must_use]
    pub fn has_unassessed_transverse_shear(&self) -> bool {
        self.force_body_n[1] != 0.0 || self.force_body_n[2] != 0.0
    }
    /// Per-layer principal tension from axial/bending stress and torsional shear.
    #[must_use]
    pub const fn assessment(&self) -> &ConcentricPrincipalTensileAssessment {
        &self.assessment
    }
}

/// Combined principal-tension observations plus original axial strand checks.
///
/// Original section and capacity sources are borrowed, without cloning material
/// property bundles. Observations describe only the configuration evaluated by
/// this call; they establish neither equilibrium nor a full failure envelope.
#[derive(Clone, Debug, PartialEq)]
pub struct RodCombinedTensileAssessment<'a> {
    segments: Vec<RodPrincipalTensileSample<'a>>,
    links: Vec<RodTensileSample<'a>>,
    governing: RodTensileMember,
    utilization: f64,
}

impl<'a> RodCombinedTensileAssessment<'a> {
    /// Native segments in their original spatial order.
    #[must_use]
    pub fn segment_samples(&self) -> &[RodPrincipalTensileSample<'a>] {
        &self.segments
    }
    /// Pin-ended links in attachment order; each sample describes one strand.
    #[must_use]
    pub fn link_samples(&self) -> &[RodTensileSample<'a>] {
        &self.links
    }
    /// Greatest utilization; segments precede links and the first exact tie wins.
    #[must_use]
    pub const fn governing_member(&self) -> RodTensileMember {
        self.governing
    }
    /// Governing member's first material layer with maximum utilization.
    #[must_use]
    pub fn governing_layer(&self) -> usize {
        match self.governing {
            RodTensileMember::Segment(index) => self.segments[index].assessment.governing_layer(),
            RodTensileMember::TensionLink(index) => self.links[index].assessment.governing_layer(),
        }
    }
    /// Maximum utilization of this explicitly selected partial stress criterion.
    #[must_use]
    pub const fn utilization(&self) -> f64 {
        self.utilization
    }
    /// Whether this selected criterion reaches a sourced limit, without tolerance.
    #[must_use]
    pub fn reaches_limit(&self) -> bool {
        self.utilization >= 1.0
    }
    /// Any nonzero transverse force whose local shear remains unassessed.
    /// A false value does not establish a complete material failure envelope.
    #[must_use]
    pub fn has_unassessed_transverse_shear(&self) -> bool {
        self.segments
            .iter()
            .any(RodPrincipalTensileSample::has_unassessed_transverse_shear)
    }
}

impl SpatialMaterialRod {
    /// Select the maximum-principal-tension criterion for native axial force,
    /// both bending moments and Saint-Venant torque on each circular section.
    ///
    /// Supply one capacity per original segment and attached tension bundle.
    /// Every capacity must retain that member's exact elastic section and state.
    /// The rod's original shear-section owner supplies G and GJ; the stiffness
    /// correction kGA is never used to invent a local transverse shear stress.
    /// Native transverse force is retained as an explicit unassessed component.
    /// Links have no torque in their pin-ended model and use their source-bound
    /// axial capacity at the native force per strand.
    ///
    /// Thermal layers retain their actual mismatch stresses and original strain
    /// limit, including the material criterion's outer-fiber torsional strain
    /// bound. The layer budget is checked before allocation; each member polls
    /// cancellation. No mechanical solve, tangent or state update is performed.
    /// The returned source references borrow this rod as well as the capacities.
    ///
    /// # Errors
    /// Count, budget, source, strain-domain and arithmetic refusals or cancellation
    /// publish no partial assessment and leave the accepted rod untouched.
    pub fn assess_axial_bending_torsion_strength<'a>(
        &'a self,
        segment_capacities: &'a [AxialTensileCapacity],
        link_capacities: &'a [AxialTensileCapacity],
        maximum_layer_samples: usize,
        gate: &CancelGate,
    ) -> Result<RodCombinedTensileAssessment<'a>, SpatialRodError> {
        checkpoint(gate)?;
        if segment_capacities.len() != self.sections.len()
            || link_capacities.len() != self.tension_links.len()
        {
            return Err(SpatialRodError::Invalid(
                "exactly one tensile capacity is required per rod segment and attached tension link",
            ));
        }
        let mut layers = 0;
        for (capacity, section) in segment_capacities.iter().zip(&self.sections) {
            checkpoint(gate)?;
            if capacity.section() != section.elastic() {
                return Err(SpatialRodError::Invalid(
                    "combined tensile capacity must retain the exact native section and material state",
                ));
            }
            layers = checked_layers(layers, capacity, maximum_layer_samples)?;
        }
        for (capacity, link) in link_capacities.iter().zip(&self.tension_links) {
            checkpoint(gate)?;
            if capacity.section() != link.section() {
                return Err(SpatialRodError::Invalid(
                    "link tensile capacity must retain the exact strand source bundle and geometry",
                ));
            }
            layers = checked_layers(layers, capacity, maximum_layer_samples)?;
        }

        let mut segments = Vec::with_capacity(segment_capacities.len());
        let mut links = Vec::with_capacity(link_capacities.len());
        let mut governing = RodTensileMember::Segment(0);
        let mut utilization = 0.0;
        for (segment, capacity) in segment_capacities.iter().enumerate() {
            checkpoint(gate)?;
            let (strain, curvature) = self.segment_strains(segment)?;
            let (force, moment) = self.resultants_from_strains(segment, strain, curvature)?;
            let shear = &self.sections[segment];
            let bending = [moment[1], moment[2]];
            let assessment = match &self.thermal {
                Some(field) => field.assess_tensile_with_torsion(
                    segment, strain[0], bending, moment[0], shear, capacity,
                )?,
                None => {
                    capacity.assess_axial_bending_torsion(shear, force[0], bending, moment[0])?
                }
            };
            if assessment.utilization() > utilization {
                utilization = assessment.utilization();
                governing = RodTensileMember::Segment(segment);
            }
            segments.push(RodPrincipalTensileSample {
                segment,
                capacity,
                shear_section: shear,
                total_axial_strain: strain[0],
                force_body_n: force,
                moment_body_n_m: moment,
                assessment,
            });
        }
        for (index, (link, capacity)) in self.tension_links.iter().zip(link_capacities).enumerate()
        {
            checkpoint(gate)?;
            let native = link.state_for(&self.rod)?;
            let axial_force_n = native.tension_n / f64::from(link.parallel_strands());
            if native.tension_n > 0.0 && axial_force_n == 0.0 {
                return Err(SpatialRodError::Invalid(
                    "tension-link force per strand is unrepresentable",
                ));
            }
            let assessment = capacity.assess_axial_bending(axial_force_n, [0.0; 2])?;
            if assessment.utilization() > utilization {
                utilization = assessment.utilization();
                governing = RodTensileMember::TensionLink(index);
            }
            links.push(RodTensileSample {
                member: RodTensileMember::TensionLink(index),
                capacity,
                total_axial_strain: native.tensile_strain,
                axial_force_n,
                bending_moment_n_m: [0.0; 2],
                assessment,
            });
        }
        checkpoint(gate)?;
        Ok(RodCombinedTensileAssessment {
            segments,
            links,
            governing,
            utilization,
        })
    }
}
