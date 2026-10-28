//! Original sourced circular Hertz geometry with exact gross sliding and spin.
//!
//! Admission reuses `SourceMindlinContact`; no normal physics, property resolver,
//! interface coefficient or source authority is recreated. This selects a
//! distinct zero-storage constitutive rung. It is not a transfer of retained
//! Mindlin shear: construction requires a source-matching unloaded checkpoint,
//! borrows it unchanged, and refuses any stored shear. Subsequent motion remains
//! owned by this sliding rung; no history migration back to Mindlin is provided.

use fs_blake3::ContentHash;
use fs_material::state_point::ResolvedInterfaceStatePoint;
use fs_tribo::TangentialLaw;
use fs_tribo::sliding_patch::{HertzSlidingPatch, HertzSlidingPatchParameters, HertzSlidingWrench};

use super::{SourceMindlinContact, SourceMindlinError, SourceMindlinState};
use crate::normal_patch::PointNormalPatchReceipt;
use crate::normal_patch::source::SourcePointContact;

/// Explicit held rigid patch motion and its original fixed source conditions.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SourceHertzSlidingInput {
    /// Supplied normal force [N], required to equal the original Hertz response.
    pub normal_force_n: f64,
    /// Centre velocity of side A relative to side B in the patch tangent axes [m/s].
    pub slip_velocity_m_s: [f64; 2],
    /// Relative angular velocity about the positive patch normal [rad/s].
    pub spin_rad_s: f64,
    /// Interval over which the supplied motion is held [s].
    pub dt_s: f64,
    /// Required original source temperature [K], without an inferred default.
    pub contact_temperature_k: Option<f64>,
    /// This frozen-query rung refuses a nominal-pressure refresh.
    pub nominal_pressure_pa: Option<f64>,
}

/// Exact sourced wrench and physical heat for one declared held-motion interval.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SourceHertzSlidingResponse {
    /// Identity of the original bound normal model and its ordered source states.
    pub normal_model_identity: ContentHash,
    /// Identity of the original scalar interface-friction query and receipts.
    pub friction_state_identity: ContentHash,
    /// Resultant on side A, opposing its relative sliding and spin.
    pub wrench: HertzSlidingWrench,
    /// Integral of physical local sliding power over the held interval [J].
    pub dissipated_j: f64,
    /// Original admitted held-motion interval [s].
    pub dt_s: f64,
}

/// Source-backed, fixed-load, constant-coefficient circular Hertz sliding patch.
#[derive(Clone, Debug, PartialEq)]
pub struct SourceHertzSlidingPatch {
    source: SourceMindlinContact,
    law: HertzSlidingPatch,
}

impl SourceHertzSlidingPatch {
    /// Select gross sliding only from a matching checkpoint with no stored shear.
    ///
    /// The checkpoint is not consumed or changed. This does not model a physical
    /// transition from elastic partial slip to sliding, release stored energy,
    /// or provide an inverse history transfer. A fresh source's initial state
    /// admits a fresh gross-sliding model; loaded elastic checkpoints refuse.
    pub fn from_unloaded_source(
        source: SourceMindlinContact,
        checkpoint: &SourceMindlinState,
    ) -> Result<Self, SourceMindlinError> {
        if !source.owns_state(checkpoint)
            || source.stored_energy_j(checkpoint) != 0.0
            || checkpoint.endpoint_force_n() != 0.0
        {
            return Err(SourceMindlinError::Invalid(
                "Hertz sliding requires its original source and an unloaded shear checkpoint",
            ));
        }
        let parameters = source.law().parameters();
        let law = HertzSlidingPatch::new(HertzSlidingPatchParameters {
            normal_force_n: source.normal_response().normal_force_n,
            contact_radius_m: source.normal_response().patch_radius_m,
            friction_coefficient: parameters.friction_coefficient,
        })?;
        Ok(Self { source, law })
    }

    /// Complete original source admission, including both bound bulk materials.
    #[must_use]
    pub const fn source(&self) -> &SourceMindlinContact {
        &self.source
    }

    /// Original normal model, geometry, material data and operating envelope.
    #[must_use]
    pub const fn normal_source(&self) -> &SourcePointContact {
        self.source.normal_source()
    }

    /// Original circular Hertz load, area, pressure and footprint receipt.
    #[must_use]
    pub const fn normal_response(&self) -> &PointNormalPatchReceipt {
        self.source.normal_response()
    }

    /// Original explicitly selected friction claim, query, uncertainty and authority.
    #[must_use]
    pub const fn friction_source(&self) -> &ResolvedInterfaceStatePoint {
        self.source.friction_source()
    }

    /// Native exact pressure-integrated constitutive law, without source promotion.
    #[must_use]
    pub const fn law(&self) -> &HertzSlidingPatch {
        &self.law
    }

    /// Evaluate prescribed two-axis slip/spin at the retained original source point.
    ///
    /// Both sides' original temperature and all external query conditions stay
    /// frozen. No replacement pressure, coefficient or normal geometry is accepted.
    /// Only the supplied motion is held over `dt_s`; this does not advance bodies,
    /// integrate a changing velocity path or allocate the exported heat.
    pub fn evaluate(
        &self,
        input: SourceHertzSlidingInput,
    ) -> Result<SourceHertzSlidingResponse, SourceMindlinError> {
        if input.normal_force_n != self.normal_response().normal_force_n
            || input.contact_temperature_k != Some(self.source.temperature_kelvin())
            || input.nominal_pressure_pa.is_some()
        {
            return Err(SourceMindlinError::Invalid(
                "Hertz sliding requires its frozen normal load, temperature and pressure query",
            ));
        }
        self.normal_source().validate_kinematics(
            self.normal_response().approach_m,
            0.0,
            input.dt_s,
        )?;
        let wrench = self
            .law
            .evaluate(input.slip_velocity_m_s, input.spin_rad_s)?;
        let dissipated_j = wrench.dissipated_power_w * input.dt_s;
        if !dissipated_j.is_finite() || dissipated_j < 0.0 {
            return Err(SourceMindlinError::Invalid(
                "Hertz sliding heat is not finite nonnegative work",
            ));
        }
        Ok(SourceHertzSlidingResponse {
            normal_model_identity: self.normal_source().model().identity(),
            friction_state_identity: self.friction_source().identity(),
            wrench,
            dissipated_j,
            dt_s: input.dt_s,
        })
    }
}
