//! Sustained contact with original full-coordinate prescribed-support work.

use super::*;
use crate::dynamics::sustained::{self as contact, Native, State, Work};
use contact::{
    SolidSustainedContact, SolidSustainedControls, SolidSustainedError, SolidSustainedPath,
    SolidSustainedStep,
};

// Binding one authored support-acceleration endpoint never mutates the source
// owner. It supplies the same input to every native contact response and retry.
struct Interval<'a> {
    body: &'a mut SupportedSpectralSolid,
    acceleration_end: &'a [f64],
}
struct Workspace {
    trial: Option<SupportedIntervalTrial>,
}

impl SupportedSpectralSolid {
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
