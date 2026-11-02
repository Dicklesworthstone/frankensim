//! Source-bound fretting through simultaneous Hertz loading and scalar shear.
//!
//! The complete accepted pressure/sliding path supplies Archard's load-distance
//! integral. Endpoint load times macroscopic travel cannot represent unloading
//! microslip, first contact or lift-off. Original normal potential, tangential
//! storage, irreversible work and wear volumes remain separate quantities.

use super::*;
use crate::interface_binding::mindlin::loading::{
    SourceHertzMindlinContact, SourceHertzMindlinState, path::SourceMindlinPathInput,
};
use crate::normal_patch::PointNormalPatchReceipt;

/// Original fixed material/Archard claims with a changing Hertz contact owner.
#[derive(Clone, Debug, PartialEq)]
pub struct SourceHertzMindlinWear {
    wear: SourceMindlinWear,
    contact: SourceHertzMindlinContact,
}

/// Complete source-bound normal/shear history and both cumulative wear volumes.
#[derive(Clone, Debug, PartialEq)]
pub struct SourceHertzMindlinWearState {
    owner: WearOwner,
    contact: SourceHertzMindlinState,
    volumes: [WearState; 2],
}
impl SourceHertzMindlinWearState {
    pub const fn contact_state(&self) -> &SourceHertzMindlinState {
        &self.contact
    }
    /// Cumulative removed volumes in original A/B roles [m^3].
    pub fn volumes_m3(&self) -> [f64; 2] {
        self.volumes.map(|volume| volume.volume_m3())
    }
}

/// One candidate carries native normal/shear work and both source wear updates.
#[derive(Clone, Debug, PartialEq)]
pub struct SourceHertzMindlinWearAdvance {
    pub tangential: TangentialStep<SourceHertzMindlinWearState>,
    pub normal_response: PointNormalPatchReceipt,
    pub normal_work_j: f64,
    pub normal_stored_change_j: f64,
    pub wear: [WearAdvance; 2],
    /// Time and area integral of p times unsigned local sliding [N m].
    pub pressure_weighted_sliding_j: f64,
}

impl SourceHertzMindlinWear {
    /// Reuse the original k/H claims and admit their contact for changing load.
    ///
    /// The fixed constructor already pins the wear query to the original
    /// interface and each hardness query to its original bulk state. The
    /// changing-load owner additionally refuses pressure, footprint, rate and
    /// other load-dependent coordinates. Temperature, mu, k and H stay fixed.
    /// No accepted history is consumed or transferred by this constructor.
    pub fn from_fixed_source(wear: SourceMindlinWear) -> Result<Self, SourceFrettingError> {
        let contact = SourceHertzMindlinContact::from_fixed_source(wear.contact.clone())?;
        Ok(Self { wear, contact })
    }

    pub const fn contact(&self) -> &SourceHertzMindlinContact {
        &self.contact
    }
    /// The complete original coefficient/hardness sources and fixed geometry.
    pub const fn fixed_source(&self) -> &SourceMindlinWear {
        &self.wear
    }

    /// Original sourced normal preload, zero shear and zero removed volume.
    pub fn initial_state(&self) -> SourceHertzMindlinWearState {
        SourceHertzMindlinWearState {
            owner: self.wear.owner,
            contact: self.contact.initial_state(),
            volumes: [WearState::default(); 2],
        }
    }

    /// Attach declared existing inventories to their original changing-load state.
    pub fn state_from_contact(
        &self,
        contact: SourceHertzMindlinState,
        volumes: [WearState; 2],
    ) -> Result<SourceHertzMindlinWearState, SourceFrettingError> {
        if !self.contact.stored_energy_j(&contact).is_finite() {
            return Err(invalid(
                "changing-load wear requires original contact history",
            ));
        }
        Ok(SourceHertzMindlinWearState {
            owner: self.wear.owner,
            contact,
            volumes,
        })
    }

    /// Propose both original wear volumes from a mechanical owner's path integral.
    ///
    /// This volume-only leaf performs no contact solve, history transfer or
    /// acceptance. The caller owns `integral p*abs(ds) dA` [N m] for the same
    /// admitted affine Hertz/shear phase and must accept these inventories with
    /// that phase's complete native candidate. An endpoint load multiplied by
    /// offered motion is not this integral. Both original k/H roles and the
    /// accepted contact identity are checked; either surface's refusal leaves
    /// the entire borrowed checkpoint unchanged.
    pub fn advance_volume_candidates(
        &self,
        state: &SourceHertzMindlinWearState,
        pressure_weighted_sliding_j: f64,
    ) -> Result<([WearState; 2], [WearAdvance; 2]), SourceFrettingError> {
        if state.owner != self.wear.owner
            || !self.contact.stored_energy_j(&state.contact).is_finite()
        {
            return Err(invalid(
                "native fretting requires original contact history and ordered k/H claims",
            ));
        }
        self.wear
            .advance_volumes(state.volumes, pressure_weighted_sliding_j)
    }

    /// Advance one affine indentation/solved-shear path and both wear inventories.
    ///
    /// The original source path resolves normal pressure and compliant shear
    /// once. Positive-mu wear uses its irreversible work divided by mu; the
    /// explicitly frictionless path uses the exact sourced Hertz mean load and
    /// accepted shear travel. Normal potential work is never treated as heat or
    /// wear. Open motion has no wear; unloading can wear at zero imposed shear.
    ///
    /// The original scalar/isothermal, rate, force and history limits remain.
    /// Cancellation or either surface's failure preserves the entire accepted
    /// checkpoint. Volumes diagnose small removal without changing geometry,
    /// mass, normal compliance, material properties or the heat partition.
    pub fn advance_coupled_path(
        &self,
        state: &SourceHertzMindlinWearState,
        input: &SourceMindlinPathInput,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<SourceHertzMindlinWearAdvance, SourceFrettingError> {
        poll(cancelled)?;
        if state.owner != self.wear.owner {
            return Err(invalid(
                "changing-load wear requires its original ordered k/H claims",
            ));
        }
        let phase =
            self.contact
                .advance_coupled_path_with_sliding(&state.contact, input, cancelled)?;
        let measure = phase.pressure_weighted_sliding_j();
        poll(cancelled)?;
        let (volumes, wear) = self.wear.advance_volumes(state.volumes, measure)?;
        let path = phase.into_path_step();
        let candidate = path.tangential;
        poll(cancelled)?;
        Ok(SourceHertzMindlinWearAdvance {
            tangential: TangentialStep {
                traction_n: candidate.traction_n,
                slip_increment_m: candidate.slip_increment_m,
                state: SourceHertzMindlinWearState {
                    owner: self.wear.owner,
                    contact: candidate.state,
                    volumes,
                },
                regime: candidate.regime,
                stored_change_j: candidate.stored_change_j,
                dissipated_j: candidate.dissipated_j,
                closure_residual_j: candidate.closure_residual_j,
                compatibility_residual_m: candidate.compatibility_residual_m,
                iterations: candidate.iterations,
            },
            normal_response: path.normal_response,
            normal_work_j: path.normal_work_j,
            normal_stored_change_j: path.normal_stored_change_j,
            wear,
            pressure_weighted_sliding_j: measure,
        })
    }
}
