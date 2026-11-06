//! Original Archard removal coupled to a finite-stiffness loading device.
//!
//! At a held device position, wear changes both contact force and approach.
//! The loading spring releases stored energy; the held actuator performs no
//! external work. The native continuation retains this spring-to-contact work
//! separately from contact storage, plastic work, and recession parameter work.

use fs_exec::Cx;

use crate::normal_patch::{
    FiniteGapPlasticDrive, FiniteGapPlasticDriveReceipt, FiniteGapPlasticLoadControls,
    FiniteGapPlasticReceipt,
};

use super::integration::{IntegratedWearAdvance, WearIntegrationControls, WearIntegrationError};
use super::{
    SourcePlasticArchardWear, SourcePlasticDriveError, SourcePlasticWearAdvance,
    SourcePlasticWearError, SourcePlasticWearState, checkpoint,
};

/// Device loading with separate original wear profiles and inventories retained.
#[derive(Clone, Debug, PartialEq)]
pub struct SourcePlasticWearNormalDriveAdvance {
    /// Complete successor contact and wear history.
    pub candidate: SourcePlasticWearState,
    /// Original native pressure, contact storage, and plastic work.
    pub receipt: FiniteGapPlasticReceipt,
    /// Finite loading-device balance and spring storage.
    pub drive: FiniteGapPlasticDriveReceipt,
}

/// A local wear step and its held-device equilibrium.
#[derive(Clone, Debug, PartialEq)]
pub struct SourcePlasticWearDriveAdvance {
    /// Both original body removals, irreversible state, and separate work accounts.
    pub wear: SourcePlasticWearAdvance,
    /// Contact/spring balance and remaining device storage.
    pub drive: FiniteGapPlasticDriveReceipt,
}

impl SourcePlasticArchardWear {
    /// Change the loading-device endpoint without resetting accepted wear.
    ///
    /// This does not infer external actuator work along an unspecified movement.
    /// The original remaining homogeneous material depth constrains the response.
    pub fn advance_normal_with_drive(
        &self,
        before: &SourcePlasticWearState,
        drive: FiniteGapPlasticDrive,
        controls: FiniteGapPlasticLoadControls,
        cx: &Cx<'_>,
    ) -> Result<SourcePlasticWearNormalDriveAdvance, SourcePlasticWearError> {
        self.admit(before)?;
        checkpoint(cx)?;
        let next = self
            .contact
            .advance_with_drive(&before.contact, drive, controls, cx)?;
        self.admit_receded_profile(&next.candidate, &before.recession_m, cx)?;
        Ok(SourcePlasticWearNormalDriveAdvance {
            candidate: SourcePlasticWearState {
                owner: self.identity,
                contact: next.candidate,
                volumes: before.volumes,
                recession_m: before.recession_m.clone(),
            },
            receipt: next.receipt,
            drive: next.drive,
        })
    }

    /// Remove material with the loading-device stiffness and position held fixed.
    ///
    /// The predecessor must already satisfy this same drive line. Every trial
    /// combines the new recession with the original plastic history directly.
    /// Each body's own original coefficient and hardness set its removal; the
    /// contact cap is not substituted for either wear denominator.
    pub fn advance_wear_with_drive(
        &self,
        before: &SourcePlasticWearState,
        local_unsigned_slip_m: &[f64],
        drive: FiniteGapPlasticDrive,
        controls: FiniteGapPlasticLoadControls,
        cx: &Cx<'_>,
    ) -> Result<SourcePlasticWearDriveAdvance, SourcePlasticWearError> {
        self.advance_wear_with_drive_over_interval(
            before,
            local_unsigned_slip_m,
            None,
            drive,
            controls,
            cx,
        )
    }

    /// Retain the loading spring while sampling wear at actual local mean speed.
    /// The positive duration belongs to the supplied irreversible sliding path.
    pub fn advance_wear_with_drive_during(
        &self,
        before: &SourcePlasticWearState,
        local_unsigned_slip_m: &[f64],
        duration_s: f64,
        drive: FiniteGapPlasticDrive,
        controls: FiniteGapPlasticLoadControls,
        cx: &Cx<'_>,
    ) -> Result<SourcePlasticWearDriveAdvance, SourcePlasticWearError> {
        self.advance_wear_with_drive_over_interval(
            before,
            local_unsigned_slip_m,
            Some(duration_s),
            drive,
            controls,
            cx,
        )
    }

    fn advance_wear_with_drive_over_interval(
        &self,
        before: &SourcePlasticWearState,
        local_unsigned_slip_m: &[f64],
        duration_s: Option<f64>,
        drive: FiniteGapPlasticDrive,
        controls: FiniteGapPlasticLoadControls,
        cx: &Cx<'_>,
    ) -> Result<SourcePlasticWearDriveAdvance, SourcePlasticWearError> {
        self.admit(before)?;
        checkpoint(cx)?;
        self.contact.admit_load_controls(controls)?;
        let (wear, drive) = self.advance_wear_during_with_recession(
            before,
            local_unsigned_slip_m,
            duration_s,
            cx,
            |increment| {
                let response = self
                    .contact
                    .model()
                    .recede_with_drive(&before.contact.contact, increment, drive, controls, cx)
                    .map_err(SourcePlasticDriveError::from)?;
                Ok((response.continuation, response.drive))
            },
        )?;
        Ok(SourcePlasticWearDriveAdvance { wear, drive })
    }

    /// Integrate a complete local sliding path with finite-stiffness load feedback.
    ///
    /// The same error-controlled driver used for fixed approach and force
    /// resolves pressure, indentation, recession, and rigid approach together.
    /// The drive is fixed over this path. Summed normal work is energy released
    /// by its spring, not work of a moving actuator. Only accepted fine histories
    /// contribute to the work and removed-volume accounts.
    pub fn integrate_wear_with_drive(
        &self,
        before: &SourcePlasticWearState,
        total_local_slip_m: &[f64],
        drive: FiniteGapPlasticDrive,
        drive_controls: FiniteGapPlasticLoadControls,
        controls: WearIntegrationControls,
        cx: &Cx<'_>,
    ) -> Result<IntegratedWearAdvance, WearIntegrationError> {
        self.contact
            .admit_load_controls(drive_controls)
            .map_err(SourcePlasticWearError::from)?;
        self.integrate_wear_with(before, total_local_slip_m, controls, cx, |state, slip| {
            self.advance_wear_with_drive(state, slip, drive, drive_controls, cx)
                .map(|advance| advance.wear)
        })
    }
}
