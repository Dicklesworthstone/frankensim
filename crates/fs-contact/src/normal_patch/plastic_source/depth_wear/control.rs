//! Original depth profiles under force or finite loading-spring control.
//!
//! Recession and exposure of the next hardness field both use the requested
//! native boundary. Each stage retains its own mechanical and control receipt;
//! no held-approach intermediate indentation is accepted.

use fs_exec::Cx;

use super::{
    DepthHardnessNormalAdvance, DepthHardnessWearAdvance, DepthHardnessWearError, RemovalRule,
    SourceDepthHardnessWear, checkpoint, within,
};
use crate::normal_patch::{
    FiniteGapPlasticDrive, FiniteGapPlasticDriveError, FiniteGapPlasticDriveReceipt,
    FiniteGapPlasticLoadControls, FiniteGapPlasticLoadError, FiniteGapPlasticLoadReceipt,
};

/// Normal force adjustment at unchanged exposed material depths.
#[derive(Clone, Debug, PartialEq)]
pub struct DepthHardnessNormalLoadAdvance {
    /// Complete original-source normal successor and native plastic response.
    pub normal: DepthHardnessNormalAdvance,
    /// Requested and achieved force, approach and native search work.
    pub load: FiniteGapPlasticLoadReceipt,
}

/// Normal loading-spring adjustment at unchanged exposed material depths.
#[derive(Clone, Debug, PartialEq)]
pub struct DepthHardnessNormalDriveAdvance {
    /// Complete original-source normal successor and native plastic response.
    pub normal: DepthHardnessNormalAdvance,
    /// Actual spring/contact balance and reversible spring storage.
    pub drive: FiniteGapPlasticDriveReceipt,
}

/// Removal and hardness exposure at one maintained compressive resultant.
#[derive(Clone, Debug, PartialEq)]
pub struct DepthHardnessWearLoadAdvance {
    /// Complete candidate, actual body inventories and both native work accounts.
    pub wear: DepthHardnessWearAdvance,
    /// Force equilibrium after removal at the original hardness field.
    pub recession_load: FiniteGapPlasticLoadReceipt,
    /// Force equilibrium after exposure of the successor hardness field.
    pub hardness_load: FiniteGapPlasticLoadReceipt,
}

/// Removal and hardness exposure at one held finite-stiffness drive position.
#[derive(Clone, Debug, PartialEq)]
pub struct DepthHardnessWearDriveAdvance {
    /// Complete candidate, actual body inventories and both native work accounts.
    pub wear: DepthHardnessWearAdvance,
    /// Spring balance and storage after removal at the original hardness field.
    pub recession_drive: FiniteGapPlasticDriveReceipt,
    /// Spring balance and storage after exposure of the successor hardness field.
    pub hardness_drive: FiniteGapPlasticDriveReceipt,
}

impl From<FiniteGapPlasticLoadError> for DepthHardnessWearError {
    fn from(error: FiniteGapPlasticLoadError) -> Self {
        Self::Load(error)
    }
}

impl From<FiniteGapPlasticDriveError> for DepthHardnessWearError {
    fn from(error: FiniteGapPlasticDriveError) -> Self {
        Self::Drive(error)
    }
}

impl SourceDepthHardnessWear {
    /// Adjust compressive force using the currently exposed hardness field.
    ///
    /// Original body inventories and irreversible indentation remain attached
    /// to the accepted history. The full search must fit the original approach
    /// envelope. An unresolved force path supplies no invented actuator work.
    pub fn advance_normal_at_load(
        &self,
        target_force_n: f64,
        controls: FiniteGapPlasticLoadControls,
        cx: &Cx<'_>,
    ) -> Result<DepthHardnessNormalLoadAdvance, DepthHardnessWearError> {
        self.admit_control_bounds(controls, cx)?;
        let native = self
            .model
            .advance_at_load(&self.contact, target_force_n, controls, cx)?;
        self.admit(
            &self.model,
            &native.contact.candidate,
            &self.recession_m,
            cx,
        )?;
        let mut candidate = self.clone();
        candidate.contact = native.contact.candidate;
        Ok(DepthHardnessNormalLoadAdvance {
            normal: DepthHardnessNormalAdvance {
                candidate,
                receipt: native.contact.receipt,
            },
            load: native.load,
        })
    }

    /// Adjust the finite loading spring using the currently exposed hardness.
    ///
    /// The native spring stiffness determines approach even at a fully plastic
    /// pressure plateau. An endpoint position change supplies no actuator path,
    /// inertia or damping; original body removal histories remain unchanged.
    pub fn advance_normal_with_drive(
        &self,
        drive: FiniteGapPlasticDrive,
        controls: FiniteGapPlasticLoadControls,
        cx: &Cx<'_>,
    ) -> Result<DepthHardnessNormalDriveAdvance, DepthHardnessWearError> {
        self.admit_control_bounds(controls, cx)?;
        let native = self
            .model
            .advance_with_drive(&self.contact, drive, controls, cx)?;
        self.admit(
            &self.model,
            &native.contact.candidate,
            &self.recession_m,
            cx,
        )?;
        let mut candidate = self.clone();
        candidate.contact = native.contact.candidate;
        Ok(DepthHardnessNormalDriveAdvance {
            normal: DepthHardnessNormalAdvance {
                candidate,
                receipt: native.contact.receipt,
            },
            drive: native.drive,
        })
    }

    /// Integrate original depth-hardness through actual removal while keeping
    /// the predecessor's compressive resultant fixed at both native stages.
    ///
    /// The original H(z) integral is exact for held beginning pressure and k.
    /// The complete contact still needs slip refinement when pressure migrates
    /// between cells. Exposing a layer unable to support the force refuses the
    /// whole step; neither a trial indentation nor body volume is committed.
    /// Both normal-work accounts use the maintained force times approach change.
    pub fn advance_wear_integrated_at_load(
        &self,
        actual_unsigned_slip_m: &[f64],
        target_force_n: f64,
        controls: FiniteGapPlasticLoadControls,
        cx: &Cx<'_>,
    ) -> Result<DepthHardnessWearLoadAdvance, DepthHardnessWearError> {
        self.admit_control_bounds(controls, cx)?;
        let mut prepared =
            self.prepare_wear(actual_unsigned_slip_m, RemovalRule::IntegratedProfile, cx)?;
        let removed = self.model.recede_at_load(
            &self.contact,
            &prepared.relative_increment_m,
            target_force_n,
            controls,
            cx,
        )?;
        let model = self.model_after_removal(&mut prepared, &removed.continuation, cx)?;
        let refreshed = model.refresh_from_at_load(
            &self.model,
            &removed.continuation.candidate,
            target_force_n,
            controls,
            cx,
        )?;
        Ok(DepthHardnessWearLoadAdvance {
            wear: self.finish_wear(
                prepared,
                model,
                removed.continuation,
                refreshed.continuation,
                cx,
            )?,
            recession_load: removed.load,
            hardness_load: refreshed.load,
        })
    }

    /// Integrate original depth-hardness while holding a finite loading spring
    /// and its position fixed through removal and exposure of the next layer.
    ///
    /// Beginning pressure and k supply each body's exact depth integral. Native
    /// spring equilibrium then supplies approach, pressure and irreversible
    /// indentation at both stages. The sum of their normal work equals the
    /// spring's energy release; it is not an inferred actuator or heat input.
    /// Refine slip to resolve pressure changes within a finite wear interval.
    pub fn advance_wear_integrated_with_drive(
        &self,
        actual_unsigned_slip_m: &[f64],
        drive: FiniteGapPlasticDrive,
        controls: FiniteGapPlasticLoadControls,
        cx: &Cx<'_>,
    ) -> Result<DepthHardnessWearDriveAdvance, DepthHardnessWearError> {
        self.admit_control_bounds(controls, cx)?;
        let mut prepared =
            self.prepare_wear(actual_unsigned_slip_m, RemovalRule::IntegratedProfile, cx)?;
        let removed = self.model.recede_with_drive(
            &self.contact,
            &prepared.relative_increment_m,
            drive,
            controls,
            cx,
        )?;
        let model = self.model_after_removal(&mut prepared, &removed.continuation, cx)?;
        let refreshed = model.refresh_from_with_drive(
            &self.model,
            &removed.continuation.candidate,
            drive,
            controls,
            cx,
        )?;
        Ok(DepthHardnessWearDriveAdvance {
            wear: self.finish_wear(
                prepared,
                model,
                removed.continuation,
                refreshed.continuation,
                cx,
            )?,
            recession_drive: removed.drive,
            hardness_drive: refreshed.drive,
        })
    }

    pub(super) fn admit_control_bounds(
        &self,
        controls: FiniteGapPlasticLoadControls,
        cx: &Cx<'_>,
    ) -> Result<(), DepthHardnessWearError> {
        checkpoint(cx)?;
        within(
            "load-control maximum approach [m]",
            controls.maximum_approach_m,
            self.sources.envelope.maximum_approach_m,
        )?;
        Ok(())
    }
}
