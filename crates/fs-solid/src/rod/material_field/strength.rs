//! Tensile normal-stress limits on the original heterogeneous rod mechanics.
//!
//! Signed axial force and both bending moments come from the same accepted
//! Cosserat strains used in equilibrium. Thermal sections retain their actual
//! per-layer stresses, including self-equilibrated mismatch at zero resultant.
//! Attached pin-ended bundles use the native force of ONE strand, not the whole
//! bundle's force applied to one strand's section. No mechanical solve or tangent
//! is assembled by this observation.
//!
//! These are local axial normal-stress comparisons on concentric, bonded,
//! small-strain sections. Transverse shear, torsional principal tension,
//! interface stresses, yield redistribution, buckling and rupture evolution are
//! not assessed. A limit observation does not sever a segment or a link.
//! [`combined`] adds a separate opt-in principal-tension criterion including
//! sourced circular torsion; it does not redefine this axial-only observation.

/// Opt-in principal tension from native axial force, bending and circular torsion.
pub mod combined;

use fs_exec::CancelGate;
use fs_material::section::strength::{AxialTensileCapacity, ConcentricTensileAssessment};

use super::{SpatialMaterialRod, SpatialRodError, checkpoint};

/// Original spatial member that supplies a tensile observation.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum RodTensileMember {
    /// A rod segment in its original material-interval order.
    Segment(usize),
    /// An attached tension bundle in its original attachment order.
    TensionLink(usize),
}

/// One physical section and its source-bound tensile stress extrema.
///
/// A tension-link sample describes one of its identical parallel strands. Its
/// force is the bundle force divided by the original strand count; its strain
/// and stress do not multiply with that count.
#[derive(Clone, Debug, PartialEq)]
pub struct RodTensileSample<'a> {
    member: RodTensileMember,
    capacity: &'a AxialTensileCapacity,
    total_axial_strain: f64,
    axial_force_n: f64,
    bending_moment_n_m: [f64; 2],
    assessment: ConcentricTensileAssessment,
}

impl<'a> RodTensileSample<'a> {
    /// Native member that supplied this section's stress state.
    #[must_use]
    pub const fn member(&self) -> RodTensileMember {
        self.member
    }
    /// Exact original section and every explicitly selected strength source.
    #[must_use]
    pub const fn capacity(&self) -> &'a AxialTensileCapacity {
        self.capacity
    }
    /// Native centerline strain, before subtracting any thermal eigenstrain.
    #[must_use]
    pub const fn total_axial_strain(&self) -> f64 {
        self.total_axial_strain
    }
    /// Signed native axial force [N]; force per strand for a tension bundle.
    #[must_use]
    pub const fn axial_force_n(&self) -> f64 {
        self.axial_force_n
    }
    /// Native body-frame bending resultants about the transverse axes [N m].
    /// Pin-ended links have zero bending moments in their declared model.
    #[must_use]
    pub const fn bending_moment_n_m(&self) -> [f64; 2] {
        self.bending_moment_n_m
    }
    /// Normal-stress extrema and utilization of each original material layer.
    #[must_use]
    pub const fn assessment(&self) -> &ConcentricTensileAssessment {
        &self.assessment
    }
}

/// Transactional tensile observation of every segment and attached member.
///
/// Capacities are borrowed so repeated observations retain full source receipts
/// without cloning material cards or their property bundles. The values are a
/// snapshot of the configuration at the time of assessment, not a future-state
/// certificate or a claim that the configuration satisfies equilibrium.
#[derive(Clone, Debug, PartialEq)]
pub struct RodTensileAssessment<'a> {
    samples: Vec<RodTensileSample<'a>>,
    governing: usize,
}

impl<'a> RodTensileAssessment<'a> {
    /// All segments first, then all attached links, in their original order.
    #[must_use]
    pub fn samples(&self) -> &[RodTensileSample<'a>] {
        &self.samples
    }
    /// Greatest utilization; the first member and layer win exact ties.
    #[must_use]
    pub fn governing_sample(&self) -> &RodTensileSample<'a> {
        &self.samples[self.governing]
    }
    /// Maximum positive axial normal stress divided by its own sourced limit.
    #[must_use]
    pub fn utilization(&self) -> f64 {
        self.governing_sample().assessment.utilization()
    }
    /// Whether some sampled layer reaches or exceeds its selected limit.
    #[must_use]
    pub fn reaches_limit(&self) -> bool {
        self.utilization() >= 1.0
    }
}

impl SpatialMaterialRod {
    /// Assess the current native section stresses against explicitly sourced
    /// tensile limits, without changing the accepted mechanical configuration.
    ///
    /// Supply exactly one capacity for every segment and every attached bundle
    /// (`&[]` for a rod without links). Each capacity must retain the EXACT
    /// elastic section, layer geometry and material states of that member.
    /// A temperature refresh therefore needs capacities for the refreshed
    /// sections; a previous state's numerically equal strength is not reused.
    ///
    /// `maximum_layer_samples` bounds the total number of material layers across
    /// these sections before allocation. Identical parallel strands share one
    /// layer assessment. Cancellation is checked for each member, and no partial
    /// report is returned. Thermal assessment also enforces the original thermal
    /// field's strain limit on the complete axial-plus-bending fiber extrema.
    ///
    /// This checks axial normal stress only. A low utilization is not a complete
    /// multiaxial strength assessment when transverse shear or torsion is present.
    ///
    /// # Errors
    /// Missing/extra capacities, foreign or stale sections, insufficient budget,
    /// source/domain/numerical refusal and cancellation leave this rod untouched.
    pub fn assess_tensile_strength<'a>(
        &self,
        segment_capacities: &'a [AxialTensileCapacity],
        link_capacities: &'a [AxialTensileCapacity],
        maximum_layer_samples: usize,
        gate: &CancelGate,
    ) -> Result<RodTensileAssessment<'a>, SpatialRodError> {
        checkpoint(gate)?;
        if segment_capacities.len() != self.sections.len()
            || link_capacities.len() != self.tension_links.len()
        {
            return Err(SpatialRodError::Invalid(
                "exactly one tensile capacity is required per rod segment and attached tension link",
            ));
        }

        let mut layers = 0usize;
        for (capacity, section) in segment_capacities.iter().zip(&self.sections) {
            checkpoint(gate)?;
            if capacity.section() != section.elastic() {
                return Err(SpatialRodError::Invalid(
                    "segment tensile capacity must retain the exact elastic source bundle and geometry",
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

        let mut samples = Vec::with_capacity(segment_capacities.len() + link_capacities.len());
        for (segment, capacity) in segment_capacities.iter().enumerate() {
            checkpoint(gate)?;
            let (strain, curvature) = self.segment_strains(segment)?;
            let (force, moment) = self.resultants_from_strains(segment, strain, curvature)?;
            let bending = [moment[1], moment[2]];
            let assessment = match &self.thermal {
                Some(field) => field.assess_tensile(segment, strain[0], bending, capacity)?,
                None => capacity.assess_axial_bending(force[0], bending)?,
            };
            samples.push(RodTensileSample {
                member: RodTensileMember::Segment(segment),
                capacity,
                total_axial_strain: strain[0],
                axial_force_n: force[0],
                bending_moment_n_m: bending,
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
            samples.push(RodTensileSample {
                member: RodTensileMember::TensionLink(index),
                capacity,
                total_axial_strain: native.tensile_strain,
                axial_force_n,
                bending_moment_n_m: [0.0; 2],
                assessment,
            });
        }
        let mut governing = 0;
        for (index, sample) in samples.iter().enumerate().skip(1) {
            if sample.assessment.utilization() > samples[governing].assessment.utilization() {
                governing = index;
            }
        }
        checkpoint(gate)?;
        Ok(RodTensileAssessment { samples, governing })
    }
}

fn checked_layers(
    current: usize,
    capacity: &AxialTensileCapacity,
    maximum: usize,
) -> Result<usize, SpatialRodError> {
    current
        .checked_add(capacity.section().layers().len())
        .filter(|&count| count <= maximum)
        .ok_or(SpatialRodError::Invalid(
            "rod tensile layer-sample budget exceeded",
        ))
}
