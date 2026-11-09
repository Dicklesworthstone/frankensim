//! Source-pinned strength of accepted thermal laminate frame stations.
//!
//! Maximum-stress checks use the shared material owner's actual free-sided
//! beam face stresses, including residual thermal stress at zero net N/M.
//! The critical sample is selected among the existing Lobatto stations; no
//! continuous maximum between independently prescribed station temperatures is
//! inferred. Exceedance is an indication and never changes native stiffness,
//! material, damage, equilibrium or heat.

use super::{Frame2d, SectionModel, checkpoint, invalid, source_error};
use crate::SolidError;
use fs_exec::Cx;
use fs_material::laminate::strength::{
    LaminaFace, LaminaFailureMode, LaminaStrengthAssessment, LaminaStrengthSelection,
    PreparedLaminateStrength,
};
use std::sync::Arc;

/// Governing accepted layer face in deterministic original traversal order.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FrameLaminateStrengthAssessment {
    /// Original frame element index.
    pub element: usize,
    /// Original Gauss–Lobatto station index, including both member ends.
    pub station: usize,
    /// Station coordinate x/L in the original member, from its start node.
    pub station_coordinate: f64,
    /// Original laminate layer index.
    pub layer: usize,
    /// This layer's own bottom or top face.
    pub face: LaminaFace,
    /// Actual material-axis stress and all five independent mode ratios.
    pub assessment: LaminaStrengthAssessment,
}

impl FrameLaminateStrengthAssessment {
    /// Largest actual face demand divided by its current selected strength.
    #[must_use]
    pub const fn utilization(&self) -> f64 {
        self.assessment.utilization()
    }

    /// Governing signed material-axis failure mode; None at exactly zero stress.
    #[must_use]
    pub fn critical_mode(&self) -> Option<LaminaFailureMode> {
        self.assessment.critical_mode()
    }

    /// Equality/exceedance flags a strength limit, without evolving damage.
    #[must_use]
    pub fn reaches_limit(&self) -> bool {
        self.assessment.reaches_limit()
    }
}

impl Frame2d {
    /// Bind maximum-stress property roles to every original layer of one
    /// existing thermal laminate member at each station's CURRENT temperature.
    /// Source cards, complete condition points and orientations come from the
    /// existing mechanical snapshots, not from independent strength inputs.
    ///
    /// Initial selection occurs once per station. Thereafter the exact selected
    /// claims refresh together with mechanics through
    /// [`Self::with_laminate_temperatures`]. A second binding refuses instead of
    /// replacing source identities. Pose, velocity, forces and inertia remain
    /// intact. Failure or cancellation returns no partially bound candidate.
    pub fn with_laminate_strength(
        &self,
        element: usize,
        selections: &[LaminaStrengthSelection],
        cx: &Cx<'_>,
    ) -> Result<Self, SolidError> {
        checkpoint(cx)?;
        let member = self
            .elements
            .get(element)
            .ok_or_else(|| invalid("laminate strength member index is out of range"))?;
        let original = member.laminate_thermal.as_ref().ok_or_else(|| {
            invalid("laminate strength requires an original source-bound thermal member")
        })?;
        if original.strength.is_some() {
            return Err(invalid(
                "laminate strength claims are already bound and cannot be reselected",
            ));
        }
        let mut selected = Vec::with_capacity(member.sections.len());
        for section in &member.sections {
            checkpoint(cx)?;
            let SectionModel::LaminateThermal(section) = section else {
                return Err(invalid(
                    "laminate strength must retain every original thermal station",
                ));
            };
            if !Arc::ptr_eq(&section.binding, original) || section.strength.is_some() {
                return Err(invalid(
                    "laminate strength cannot replace an existing station source",
                ));
            }
            let prepared =
                PreparedLaminateStrength::resolve((*section.resolved).clone(), selections, || {
                    cx.checkpoint().is_err()
                })
                .map_err(|error| source_error(error, cx))?;
            selected.push(Arc::new(prepared));
        }
        let mut binding = (**original).clone();
        binding.strength = Some(selected.clone());
        let binding = Arc::new(binding);
        let mut candidate = self.clone();
        for (section, prepared) in candidate.elements[element]
            .sections
            .iter_mut()
            .zip(selected)
        {
            let SectionModel::LaminateThermal(section) = section else {
                unreachable!("every original station was admitted above")
            };
            section.binding = Arc::clone(&binding);
            section.strength = Some(prepared);
        }
        candidate.elements[element].laminate_thermal = Some(binding);
        checkpoint(cx)?;
        Ok(candidate)
    }

    /// Assess all accepted original stations of one selected laminate member.
    /// Uses actual layer thermal self-stress even when net force/moment vanish.
    /// No source queries, trigonometry or allocation occur during assessment.
    /// Exact ties retain the earlier station, layer, bottom face and failure
    /// mode according to the shared material owner's deterministic ordering.
    pub fn assess_laminate_element_strength(
        &self,
        element: usize,
        cx: &Cx<'_>,
    ) -> Result<FrameLaminateStrengthAssessment, SolidError> {
        checkpoint(cx)?;
        let member = self
            .elements
            .get(element)
            .ok_or_else(|| invalid("laminate strength member index is out of range"))?;
        let original = member.laminate_thermal.as_ref().ok_or_else(|| {
            invalid("laminate assessment requires an original source-bound thermal member")
        })?;
        if original.strength.is_none() {
            return Err(invalid(
                "laminate assessment requires explicitly bound strength claims",
            ));
        }
        let mut worst: Option<FrameLaminateStrengthAssessment> = None;
        for (station, section) in member.sections.iter().enumerate() {
            checkpoint(cx)?;
            let SectionModel::LaminateThermal(section) = section else {
                return Err(invalid(
                    "laminate assessment encountered a foreign station model",
                ));
            };
            if !Arc::ptr_eq(&section.binding, original) {
                return Err(invalid(
                    "laminate assessment encountered a foreign station binding",
                ));
            }
            let strength = section.strength.as_ref().ok_or_else(|| {
                invalid("laminate assessment requires every original station strength")
            })?;
            let (eps0, kappa) = member.commit.d[station];
            section.admit_strain(eps0, kappa)?;
            let response = strength
                .assess_beam(&section.prepared, eps0, kappa, || cx.checkpoint().is_err())
                .map_err(|error| source_error(error, cx))?;
            let current = FrameLaminateStrengthAssessment {
                element,
                station,
                station_coordinate: member.xi[station],
                layer: response.layer,
                face: response.face,
                assessment: response.assessment,
            };
            if worst
                .as_ref()
                .is_none_or(|old| current.utilization() > old.utilization())
            {
                worst = Some(current);
            }
        }
        checkpoint(cx)?;
        worst.ok_or_else(|| invalid("laminate strength member has no original stations"))
    }

    /// Governing accepted face over EVERY member of a laminate frame.
    /// Refuses an empty frame or any member without its own bound source
    /// strengths, so an unchecked member cannot silently appear safe. Mixed
    /// frames can assess individual members through the element method.
    pub fn assess_laminate_frame_strength(
        &self,
        cx: &Cx<'_>,
    ) -> Result<FrameLaminateStrengthAssessment, SolidError> {
        checkpoint(cx)?;
        let mut worst: Option<FrameLaminateStrengthAssessment> = None;
        for element in 0..self.elements.len() {
            let current = self.assess_laminate_element_strength(element, cx)?;
            if worst
                .as_ref()
                .is_none_or(|old| current.utilization() > old.utilization())
            {
                worst = Some(current);
            }
        }
        checkpoint(cx)?;
        worst.ok_or_else(|| invalid("laminate frame strength requires at least one member"))
    }
}
