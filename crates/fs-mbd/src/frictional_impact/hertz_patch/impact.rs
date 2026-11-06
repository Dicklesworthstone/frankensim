//! One coupled normal and finite-patch friction impact on original bodies.
//!
//! The Hertz radius and pressure-impulse shape are declared effective impact
//! data. Normal restitution and sliding/spin share one four-coordinate body
//! mobility and one normal impulse. This instantaneous model supplies no contact
//! duration, evolving Hertz radius, elastic history or thermal allocation.

use super::{
    HertzImpulsePatch, MassProperties, PreparedRigidHertzImpulse, RigidBodyState,
    RigidHertzImpulseContact, RigidHertzImpulseControls, RigidHertzImpulseCounterpart,
    RigidHertzImpulseError, RigidHertzImpulseReport,
};
use fs_contact::tangential::impact::hertz::{
    HertzImpactControls, HertzImpactError, HertzImpactInput, HertzImpactReport, solve_hertz_impact,
};
use fs_math::dd::Dd;

/// Native-body gates and coupled contact controls. The two embedded patch
/// controls must agree, so preparation, contact projection and native admission
/// retain the same physical impulse bounds and constitutive tolerances.
#[derive(Clone, Copy, Debug)]
pub struct RigidHertzImpactControls {
    pub native: RigidHertzImpulseControls,
    pub impact: HertzImpactControls,
}

/// Complete coupled contact and original-body result. Both nested reports keep
/// their measured residuals; mechanical loss is not assigned to a temperature.
#[derive(Clone, Debug, PartialEq)]
#[must_use]
pub struct RigidHertzImpactReport {
    pub impact: HertzImpactReport,
    pub mechanics: RigidHertzImpulseReport,
    /// Native normal velocity minus the contact law's restitution target [m/s].
    /// An inactive separating contact may retain a strictly positive gap.
    pub normal_law_gap_m_s: f64,
    /// Actual normal impulse times that native gap [J].
    pub normal_complementarity_j: f64,
}

/// Failed preparation, constitutive resolution or native audit never publishes
/// any of the caller's momenta.
#[derive(Debug)]
pub enum RigidHertzImpactError {
    Mechanics(RigidHertzImpulseError),
    Contact(HertzImpactError),
    Cancelled,
    Audit(Box<RigidHertzImpactReport>),
}

impl From<RigidHertzImpulseError> for RigidHertzImpactError {
    fn from(error: RigidHertzImpulseError) -> Self {
        match error {
            RigidHertzImpulseError::Cancelled => Self::Cancelled,
            other => Self::Mechanics(other),
        }
    }
}
impl core::fmt::Display for RigidHertzImpactError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "native coupled Hertz impact: {self:?}")
    }
}
impl std::error::Error for RigidHertzImpactError {}

impl PreparedRigidHertzImpulse {
    /// Original response in work-conjugate coordinates
    /// `[normal, tangent_x, tangent_y, radius * spin]` [1/kg].
    /// All normal/tangent/spin cross terms come from the same full body inertia.
    #[must_use]
    pub fn impact_mobility_per_kg(&self) -> [[f64; 4]; 4] {
        let mut matrix = [[0.0; 4]; 4];
        matrix[0][0] = self.normal_mobility;
        for i in 0..3 {
            matrix[0][i + 1] = self.normal_cross[i];
            matrix[i + 1][0] = self.normal_cross[i];
            for j in 0..3 {
                matrix[i + 1][j + 1] = self.mobility[i][j];
            }
        }
        matrix
    }

    /// Resolve restitution and finite-patch friction on this immutable native
    /// snapshot, then apply their summed impulse once to each candidate body.
    /// Repeated coefficient/restitution trials reuse the original geometry;
    /// radius changes require a new preparation.
    pub fn solve_impact(
        &self,
        law: &HertzImpulsePatch,
        restitution_coefficient: f64,
        controls: HertzImpactControls,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<RigidHertzImpactReport, RigidHertzImpactError> {
        if cancelled() {
            return Err(RigidHertzImpactError::Cancelled);
        }
        if controls.patch != self.controls.contact {
            return Err(RigidHertzImpulseError::Invalid(
                "prepared and coupled Hertz patch controls must agree",
            )
            .into());
        }
        if law.radius_m().to_bits() != self.radius_m.to_bits() {
            return Err(RigidHertzImpulseError::Invalid(
                "prepared Hertz radius must remain unchanged",
            )
            .into());
        }
        let impact = solve_hertz_impact(
            law,
            &HertzImpactInput {
                velocity_before_m_s: [
                    self.before.normal,
                    self.before.tangent[0],
                    self.before.tangent[1],
                    self.radius_m * self.before.spin,
                ],
                scaled_mobility_per_kg: self.impact_mobility_per_kg(),
                restitution_coefficient,
            },
            controls,
            &mut cancelled,
        )
        .map_err(|error| match error {
            HertzImpactError::Cancelled => RigidHertzImpactError::Cancelled,
            other => RigidHertzImpactError::Contact(other),
        })?;
        let mechanics = self.apply_solution(
            impact.patch,
            impact.scaled_impulse_ns[0],
            impact.dissipated_j,
            self.controls
                .energy_tolerance_j
                .min(controls.energy_tolerance_j),
            &mut cancelled,
        )?;
        let gap = (Dd::from_f64(mechanics.normal_velocity_after_m_s)
            - Dd::from_f64(impact.normal_target_m_s))
        .to_f64();
        let complementarity =
            (Dd::from_f64(impact.scaled_impulse_ns[0]) * Dd::from_f64(gap)).to_f64();
        let report = RigidHertzImpactReport {
            impact,
            mechanics,
            normal_law_gap_m_s: gap,
            normal_complementarity_j: complementarity,
        };
        let tolerance = self
            .controls
            .velocity_tolerance_m_s
            .min(controls.normal_velocity_tolerance_m_s);
        if !gap.is_finite()
            || !complementarity.is_finite()
            || gap < -tolerance
            || (report.impact.scaled_impulse_ns[0] > 0.0 && gap.abs() > tolerance)
            || complementarity.abs()
                > self
                    .controls
                    .energy_tolerance_j
                    .min(controls.energy_tolerance_j)
        {
            return Err(RigidHertzImpactError::Audit(Box::new(report)));
        }
        if cancelled() {
            return Err(RigidHertzImpactError::Cancelled);
        }
        Ok(report)
    }
}

/// Publish one coupled Hertz impact atomically on the original one or two
/// bodies. The supplied restitution is resolved together with friction; no
/// separate caller-applied normal impulse belongs to this operation.
pub fn apply_rigid_hertz_impact(
    bodies: &mut [RigidBodyState],
    properties: &[MassProperties],
    contact: RigidHertzImpulseContact<'_>,
    restitution_coefficient: f64,
    controls: RigidHertzImpactControls,
    mut cancelled: impl FnMut() -> bool,
) -> Result<RigidHertzImpactReport, RigidHertzImpactError> {
    let prepared = PreparedRigidHertzImpulse::new(
        bodies,
        properties,
        contact,
        controls.native,
        &mut cancelled,
    )?;
    let report = prepared.solve_impact(
        contact.law,
        restitution_coefficient,
        controls.impact,
        &mut cancelled,
    )?;
    if cancelled() {
        return Err(RigidHertzImpactError::Cancelled);
    }
    bodies[contact.point.body] = report.mechanics.point_impulse.state_after;
    if let (RigidHertzImpulseCounterpart::Body(point), Some(receipt)) =
        (contact.other, report.mechanics.counterpart_impulse)
    {
        bodies[point.body] = receipt.state_after;
    }
    Ok(report)
}
