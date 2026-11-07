//! Sustained contact with original full-coordinate prescribed-support work.

use super::*;
use crate::dynamics::sustained::{self as contact, Native, State, Work};
use contact::{
    SolidSustainedContact, SolidSustainedControls, SolidSustainedError, SolidSustainedPath,
    SolidSustainedStep,
};

// Binding one authored support-acceleration endpoint never mutates the source
// owner. It supplies the same input to every native contact response and retry.
pub(super) struct Interval<'a> {
    pub(super) body: &'a mut SupportedSpectralSolid,
    pub(super) acceleration_end: &'a [f64],
}
pub(super) struct Workspace {
    pub(super) trial: Option<SupportedIntervalTrial>,
}

impl SupportedSpectralSolid {
    /// Resolve distinct static and kinetic friction ellipses while retaining
    /// the actual prescribed support history, full reactions and actuator work.
    /// Each kinetic pair weakens the corresponding static principal strengths
    /// on their unchanged material axes. Every coupled regime trial starts at
    /// the same accepted native state and shares the original work budgets.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_sustained_static_kinetic_elliptic_contacts<P, E>(
        &mut self,
        support_acceleration_end_mps2: &[f64],
        contacts: &[contact::elliptic::SolidSustainedEllipticContact<'_>],
        kinetic_coefficients_mu: &[[f64; 2]],
        initial_forces_n: &[[f64; 3]],
        external_force_n: &[f64],
        controls: SolidSustainedControls,
        cx: &Cx<'_>,
        admit_path: impl FnOnce(
            &contact::elliptic::SolidSustainedEllipticPath<'_, '_, Self>,
        ) -> core::result::Result<P, E>,
    ) -> core::result::Result<
        contact::elliptic::SolidSustainedEllipticStep<P, SupportedSolidStep>,
        SolidSustainedError<E, SupportedSolidError>,
    > {
        contact::elliptic::advance_static_kinetic(
            &mut Interval {
                body: self,
                acceleration_end: support_acceleration_end_mps2,
            },
            contacts,
            kinetic_coefficients_mu,
            initial_forces_n,
            external_force_n,
            controls,
            cx,
            admit_path,
        )
    }

    /// Advance an oriented Coulomb ellipse through the same original complete
    /// response, work gates and atomic acceptance as sustained disk contact.
    /// Physical contacts/forces reach the mandatory whole-path callback; the
    /// returned nonlinear diagnostics explicitly retain solver coordinates.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_sustained_elliptic_contacts<P, E>(
        &mut self,
        support_acceleration_end_mps2: &[f64],
        contacts: &[contact::elliptic::SolidSustainedEllipticContact<'_>],
        initial_forces_n: &[[f64; 3]],
        external_force_n: &[f64],
        controls: SolidSustainedControls,
        cx: &Cx<'_>,
        admit_path: impl FnOnce(
            &contact::elliptic::SolidSustainedEllipticPath<'_, '_, Self>,
        ) -> core::result::Result<P, E>,
    ) -> core::result::Result<
        contact::elliptic::SolidSustainedEllipticStep<P, SupportedSolidStep>,
        SolidSustainedError<E, SupportedSolidError>,
    > {
        contact::elliptic::advance(
            &mut Interval {
                body: self,
                acceleration_end: support_acceleration_end_mps2,
            },
            contacts,
            initial_forces_n,
            external_force_n,
            controls,
            cx,
            admit_path,
        )
    }

    /// Advance held external/contact forces while support acceleration follows
    /// the original linear interpolation to the supplied endpoint.
    ///
    /// Contact covectors and external forces use ALL `represented_dofs()`,
    /// including prescribed coordinates. Actual support motion participates in
    /// contact velocity/slip and work. Contact loads at a prescribed coordinate
    /// change its reaction, without inventing mobility or dropping its work.
    ///
    /// The original cubic free source retains K_fc and M_fc, every free mode,
    /// full-row reactions and independently integrated support-actuator work.
    /// The shared sustained-contact solve and whole-path admission add no
    /// penalty stiffness, impulse approximation or thermal loss partition.
    /// Native free/full work gates, contact-mode gates and final cancellation
    /// precede simultaneous publication of both original states and accounts.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_sustained_contacts<P, E>(
        &mut self,
        support_acceleration_end_mps2: &[f64],
        contacts: &[SolidSustainedContact<'_>],
        initial_forces_n: &[[f64; 3]],
        external_force_n: &[f64],
        controls: SolidSustainedControls,
        cx: &Cx<'_>,
        admit_path: impl FnOnce(&SolidSustainedPath<'_, '_, Self>) -> core::result::Result<P, E>,
    ) -> core::result::Result<
        SolidSustainedStep<P, SupportedSolidStep>,
        SolidSustainedError<E, SupportedSolidError>,
    > {
        contact::advance(
            &mut Interval {
                body: self,
                acceleration_end: support_acceleration_end_mps2,
            },
            contacts,
            initial_forces_n,
            external_force_n,
            controls,
            cx,
            admit_path,
        )
    }

    /// Resolve static capture and kinetic sliding with original support motion
    /// and actuator work. Contact coefficients are static, and kinetic bounds
    /// are supplied in the same contact order with `0 <= mu_k <= mu_s`.
    /// All coupled regime re-solves retain the same predecessor, support input
    /// and global response/work budgets. The path callback observes the selected
    /// bounds and must still exclude interior stopping or mode changes.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_sustained_static_kinetic_contacts<P, E>(
        &mut self,
        support_acceleration_end_mps2: &[f64],
        contacts: &[SolidSustainedContact<'_>],
        kinetic_coefficients: &[f64],
        initial_forces_n: &[[f64; 3]],
        external_force_n: &[f64],
        controls: SolidSustainedControls,
        cx: &Cx<'_>,
        admit_path: impl FnOnce(&SolidSustainedPath<'_, '_, Self>) -> core::result::Result<P, E>,
    ) -> core::result::Result<
        SolidSustainedStep<P, SupportedSolidStep>,
        SolidSustainedError<E, SupportedSolidError>,
    > {
        contact::advance_with_kinetic(
            &mut Interval {
                body: self,
                acceleration_end: support_acceleration_end_mps2,
            },
            contacts,
            Some(kinetic_coefficients),
            initial_forces_n,
            external_force_n,
            controls,
            cx,
            admit_path,
        )
    }
}

impl Native for Interval<'_> {
    type Owner = SupportedSpectralSolid;
    type Step = SupportedSolidStep;
    type Error = SupportedSolidError;
    type Workspace = Workspace;

    fn owner(&self) -> &SupportedSpectralSolid {
        self.body
    }
    fn state(&self) -> State<'_> {
        State {
            displacement: &self.body.displacement,
            velocity: &self.body.velocity,
            mass: &self.body.mass,
            stiffness: &self.body.stiffness,
            dt_s: self.body.dt_s,
        }
    }
    fn admit(&self, cx: &Cx<'_>) -> Result<()> {
        poll(cx)?;
        if self.acceleration_end.len() != self.body.support_rows.len()
            || self.acceleration_end.iter().any(|value| !value.is_finite())
        {
            return Err(invalid(
                "complete finite support acceleration endpoint required",
            ));
        }
        Ok(())
    }
    fn cancelled() -> Self::Error {
        SupportedSolidError::Cancelled
    }
    fn is_cancelled(error: &Self::Error) -> bool {
        matches!(error, SupportedSolidError::Cancelled)
    }
    fn invalid(what: &'static str) -> Self::Error {
        super::invalid(what)
    }
    fn workspace_bytes(&self) -> u128 {
        // Conservative simultaneous old/new native cubic trial, full-coordinate
        // reaction/work scratch and retained response. No prepared matrix copy.
        8 * (128 * self.body.coordinates.len() as u128 + 64)
    }
    fn response_multiply_terms(&self) -> u128 {
        let n = self.body.coordinates.len() as u128;
        // The original full supported interval bound; each response also keeps
        // the original free and full-row work/reaction admission active.
        64 * n * n + 8192 * n
    }
    fn workspace(&self) -> Result<Workspace> {
        Ok(Workspace { trial: None })
    }
    fn evaluate(&self, force: &[f64], scratch: &mut Workspace, cx: &Cx<'_>) -> Result<()> {
        let trial = self
            .body
            .trial_affine_supports(self.acceleration_end, force, force, cx)?;
        scratch.trial = Some(trial);
        Ok(())
    }
    fn displacement<'a>(&self, scratch: &'a Workspace) -> &'a [f64] {
        &scratch
            .trial
            .as_ref()
            .expect("successful native response precedes observation")
            .displacement
    }
    fn velocity<'a>(&self, scratch: &'a Workspace) -> &'a [f64] {
        &scratch
            .trial
            .as_ref()
            .expect("successful native response precedes observation")
            .velocity
    }
    fn work(&self, _: &[f64], scratch: &Workspace, _: &Cx<'_>) -> Result<Work> {
        let energy = scratch
            .trial
            .as_ref()
            .expect("successful native response precedes observation")
            .output
            .energy;
        Ok(Work {
            applied_j: energy.applied_work_j,
            viscosity_j: 0.0,
            support_j: energy.support_work_j,
        })
    }
    fn accept(&mut self, scratch: Workspace, _: &[f64], cx: &Cx<'_>) -> Result<SupportedSolidStep> {
        poll(cx)?;
        Ok(self.body.accept_affine_supports(
            scratch
                .trial
                .expect("successful final native response precedes acceptance"),
        ))
    }
}
