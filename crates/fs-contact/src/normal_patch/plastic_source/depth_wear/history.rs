//! Independently accepted native preconditioning on retained worn surfaces.
//!
//! Replacement profiles have already passed complete J2/cyclic admission and
//! continuous positive depth support. This adapter retains original source/test
//! pins, family/branch order, fixed environment and nondecreasing native alpha.
//! Each updated cap is evaluated at the body's actual retained exposed depth.
//! The existing normal continuation owns all ensuing indentation, external work
//! and signed parameter work. No bulk deformation, heat, residual stress or
//! plastic work is transferred, and no intermediate preconditioning path is
//! certified. The independently accepted alpha must represent the whole column.

use std::sync::Arc;

use fs_blake3::ContentHash;
use fs_exec::Cx;
use fs_material::hardness::depth::ResolvedDepthHardness;
use fs_material::hardness::plastic_history::ACCUMULATED_EQUIVALENT_PLASTIC_STRAIN_AXIS as ALPHA;

use super::{
    DepthHardnessWearError, SourceDepthHardnessWear, Sources, admit_depth_sources, checkpoint,
    depth_source_identity, invalid,
};
use crate::normal_patch::{
    FiniteGapPlasticContinuation, FiniteGapPlasticContinuationReceipt, FiniteGapPlasticDrive,
    FiniteGapPlasticDriveReceipt, FiniteGapPlasticLoadControls, FiniteGapPlasticLoadReceipt,
    FiniteGapPlasticModel,
};

/// Complete refreshed contact/wear owner and its original native work account.
#[derive(Clone, Debug, PartialEq)]
pub struct DepthHardnessHistoryRefresh {
    /// Accepted source/configuration identity before the profile replacement.
    pub original_source_identity: ContentHash,
    /// Complete successor, preserving removed volume and both exposed surfaces.
    pub candidate: SourceDepthHardnessWear,
    /// Native storage, irreversible indentation and signed parameter work.
    pub receipt: FiniteGapPlasticContinuationReceipt,
}

/// Profile replacement directly at the maintained compressive resultant.
#[derive(Clone, Debug, PartialEq)]
pub struct DepthHardnessHistoryLoadRefresh {
    /// Complete profile/contact/wear successor and native work account.
    pub history: DepthHardnessHistoryRefresh,
    /// Original native force equilibrium, approach and search receipt.
    pub load: FiniteGapPlasticLoadReceipt,
}

/// Profile replacement at the original finite loading spring and position.
#[derive(Clone, Debug, PartialEq)]
pub struct DepthHardnessHistoryDriveRefresh {
    /// Complete profile/contact/wear successor and native work account.
    pub history: DepthHardnessHistoryRefresh,
    /// Original native spring equilibrium and reversible spring storage.
    pub drive: FiniteGapPlasticDriveReceipt,
}

struct PreparedProfiles {
    sources: Arc<Sources>,
    model: FiniteGapPlasticModel,
}

impl SourceDepthHardnessWear {
    /// Adopt complete independently accepted H(z, alpha) profiles at the held
    /// approach. Native alpha may increase while its original hardness source
    /// softens or hardens. Rehardening cannot heal an already retained dent.
    ///
    /// Both ordered bodies must keep their original material, indentation test,
    /// selected claim, complete native source family, fixed conditions and depth
    /// ceiling. An ordinary fixed-environment side must be exactly unchanged.
    /// This endpoint-only property transfer preserves all removal inventories.
    pub fn refresh_hardness_profiles(
        &self,
        hardness: [ResolvedDepthHardness; 2],
        cx: &Cx<'_>,
    ) -> Result<DepthHardnessHistoryRefresh, DepthHardnessWearError> {
        let prepared = self.prepare_history_profiles(hardness, cx)?;
        let next = prepared
            .model
            .refresh_from(&self.model, &self.contact, cx)?;
        self.finish_history_profiles(prepared, next, cx)
    }

    /// Adopt successor profiles directly at the maintained force. The accepted
    /// predecessor must already support this same force, and a softened profile
    /// unable to support it refuses before any intermediate dent is published.
    pub fn refresh_hardness_profiles_at_load(
        &self,
        hardness: [ResolvedDepthHardness; 2],
        force_n: f64,
        controls: FiniteGapPlasticLoadControls,
        cx: &Cx<'_>,
    ) -> Result<DepthHardnessHistoryLoadRefresh, DepthHardnessWearError> {
        self.admit_control_bounds(controls, cx)?;
        let prepared = self.prepare_history_profiles(hardness, cx)?;
        let next = prepared.model.refresh_from_at_load(
            &self.model,
            &self.contact,
            force_n,
            controls,
            cx,
        )?;
        Ok(DepthHardnessHistoryLoadRefresh {
            history: self.finish_history_profiles(prepared, next.continuation, cx)?,
            load: next.load,
        })
    }

    /// Adopt successor profiles while holding the accepted finite spring and
    /// position. Its reversible storage supplies normal work; parameter work
    /// remains a separate signed account and is not interpreted as heat.
    pub fn refresh_hardness_profiles_with_drive(
        &self,
        hardness: [ResolvedDepthHardness; 2],
        drive: FiniteGapPlasticDrive,
        controls: FiniteGapPlasticLoadControls,
        cx: &Cx<'_>,
    ) -> Result<DepthHardnessHistoryDriveRefresh, DepthHardnessWearError> {
        self.admit_control_bounds(controls, cx)?;
        let prepared = self.prepare_history_profiles(hardness, cx)?;
        let next = prepared.model.refresh_from_with_drive(
            &self.model,
            &self.contact,
            drive,
            controls,
            cx,
        )?;
        Ok(DepthHardnessHistoryDriveRefresh {
            history: self.finish_history_profiles(prepared, next.continuation, cx)?,
            drive: next.drive,
        })
    }

    fn prepare_history_profiles(
        &self,
        hardness: [ResolvedDepthHardness; 2],
        cx: &Cx<'_>,
    ) -> Result<PreparedProfiles, DepthHardnessWearError> {
        checkpoint(cx)?;
        for side in 0..2 {
            admit_successor(&self.sources.hardness[side], &hardness[side])?;
        }
        let (k, selected) = admit_depth_sources(
            &self.sources.normal,
            &self.sources.coefficients,
            &self.sources.coefficient_keys,
            &hardness,
            self.sources.envelope,
            self.sources.maximum_recession_m,
        )?;
        if hardness == self.sources.hardness {
            return Ok(PreparedProfiles {
                sources: Arc::clone(&self.sources),
                model: self.model.clone(),
            });
        }
        let mut caps = Vec::with_capacity(self.contact.pressure_pa().len());
        for cell in 0..self.contact.pressure_pa().len() {
            if cell % 64 == 0 {
                checkpoint(cx)?;
            }
            caps.push(
                hardness[0]
                    .hardness_pa_at(self.recession_m[0][cell])?
                    .min(hardness[1].hardness_pa_at(self.recession_m[1][cell])?),
            );
        }
        let model = self.model.with_material_properties(
            self.model.config().reduced_modulus_pa,
            caps,
            cx,
        )?;
        let identity = depth_source_identity(
            &self.sources.normal,
            &self.sources.coefficients,
            &hardness,
            &selected,
            &model,
            self.sources.envelope,
            self.sources.maximum_recession_m,
        );
        let sources = Arc::new(Sources {
            normal: self.sources.normal.clone(),
            coefficients: self.sources.coefficients.clone(),
            coefficient_keys: self.sources.coefficient_keys.clone(),
            k,
            hardness,
            envelope: self.sources.envelope,
            maximum_recession_m: self.sources.maximum_recession_m,
            identity,
        });
        checkpoint(cx)?;
        Ok(PreparedProfiles { sources, model })
    }

    fn finish_history_profiles(
        &self,
        prepared: PreparedProfiles,
        next: FiniteGapPlasticContinuation,
        cx: &Cx<'_>,
    ) -> Result<DepthHardnessHistoryRefresh, DepthHardnessWearError> {
        self.admit(&prepared.model, &next.candidate, &self.recession_m, cx)?;
        checkpoint(cx)?;
        Ok(DepthHardnessHistoryRefresh {
            original_source_identity: self.source_identity(),
            candidate: Self {
                sources: prepared.sources,
                model: prepared.model,
                contact: next.candidate,
                volumes: self.volumes,
                recession_m: self.recession_m.clone(),
            },
            receipt: next.receipt,
        })
    }
}

fn admit_successor(
    before: &ResolvedDepthHardness,
    after: &ResolvedDepthHardness,
) -> Result<(), DepthHardnessWearError> {
    match (before.accepted_history(), after.accepted_history()) {
        (None, None) if before == after => Ok(()),
        (Some(old), Some(next)) => {
            let a = before.reference_source();
            let b = after.reference_source();
            let ar = a.resolved();
            let br = b.resolved();
            if old.native_source_identity() != next.native_source_identity()
                || old.mechanical_source() != next.mechanical_source()
                || next.alpha() < old.alpha()
                || before.maximum_depth_m() != after.maximum_depth_m()
                || a.key() != b.key()
                || a.property().answer().receipt.selected != b.property().answer().receipt.selected
                || ar.card_identity() != br.card_identity()
                || ar.material() != br.material()
                || ar.axis_quantities() != br.axis_quantities()
                || !ar
                    .query_point()
                    .iter()
                    .filter(|(name, _)| name.as_str() != ALPHA)
                    .eq(br
                        .query_point()
                        .iter()
                        .filter(|(name, _)| name.as_str() != ALPHA))
            {
                return Err(invalid(
                    "profile refresh requires original test, claim, native family, conditions and depth support with nondecreasing alpha",
                ));
            }
            Ok(())
        }
        _ => Err(invalid(
            "profile refresh must retain native-history ownership and leave each ordinary fixed side unchanged",
        )),
    }
}
