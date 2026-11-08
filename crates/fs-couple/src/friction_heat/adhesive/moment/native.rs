//! Native source/history dispatch; neither mode accepts caller-authored losses.

use super::{AdhesiveMomentHeatError, AdhesiveMomentMechanicalControls};
use fs_contact::normal_patch::adhesive::{JkrResponse, SourceJkrSphereContact};
use fs_exec::CancelGate;
use fs_mbd::tangential_friction::RigidContactBinding;
use fs_mbd::tangential_friction::adhesive::{rolling::*, torsion::*};
use fs_mbd::{MassProperties, RigidBodyState, Wrench};
use fs_tribo::partial_slip::TangentFrame;

/// An original physical contact and its sealed native memory, including preload.
/// The contact never changes when its history advances through the heat owner.
#[derive(Clone, Copy, Debug)]
pub enum NativeAdhesiveMoment<'s> {
    Rolling {
        contact: RigidAdhesiveRollingContact<'s>,
        state: RigidAdhesiveRollingState,
    },
    Torsion {
        contact: RigidAdhesiveTorsionContact<'s>,
        state: RigidAdhesiveTorsionState,
    },
}

#[derive(Clone, Debug, PartialEq)]
pub enum AdhesiveMomentMechanicalReport {
    Rolling(RigidAdhesiveRollingReport),
    Torsion(RigidAdhesiveTorsionReport),
}

pub(super) struct NativeStep<'s> {
    pub bodies_after: [RigidBodyState; 2],
    pub after: NativeAdhesiveMoment<'s>,
    pub report: AdhesiveMomentMechanicalReport,
}
pub(super) struct Balance {
    pub kinetic_change_j: f64,
    pub storage_change_j: f64,
    pub irreversible_loss_j: f64,
    pub external_work_j: f64,
    pub residual_j: f64,
}
impl AdhesiveMomentMechanicalReport {
    pub(super) fn balance(&self) -> Balance {
        match self {
            Self::Rolling(r) => Balance {
                kinetic_change_j: r.kinetic_energy_change_j,
                storage_change_j: r.contact.rolling().storage_change_j,
                irreversible_loss_j: r.contact.rolling().irreversible_loss_j,
                external_work_j: r.external_work_j,
                residual_j: r.total_energy_residual_j,
            },
            Self::Torsion(r) => Balance {
                kinetic_change_j: r.kinetic_energy_change_j,
                storage_change_j: r.contact.torsion().storage_change_j,
                irreversible_loss_j: r.contact.torsion().irreversible_loss_j,
                external_work_j: r.external_work_j,
                residual_j: r.total_energy_residual_j,
            },
        }
    }
}

impl<'s> NativeAdhesiveMoment<'s> {
    pub(super) fn normal_source(self) -> &'s SourceJkrSphereContact {
        match self {
            Self::Rolling { contact, .. } => contact.source.normal_source(),
            Self::Torsion { contact, .. } => contact.source.normal_source(),
        }
    }
    pub(super) fn normal_response(self) -> &'s JkrResponse {
        match self {
            Self::Rolling { contact, .. } => contact.source.normal_response(),
            Self::Torsion { contact, .. } => contact.source.normal_response(),
        }
    }
    pub(super) fn binding(self) -> RigidContactBinding {
        match self {
            Self::Rolling { contact, .. } => contact.binding,
            Self::Torsion { contact, .. } => contact.binding,
        }
    }
    pub(super) fn frame(self) -> TangentFrame {
        match self {
            Self::Rolling { contact, .. } => contact.frame,
            Self::Torsion { contact, .. } => contact.frame,
        }
    }
    pub(super) fn stored_energy_j<E>(self) -> Result<f64, AdhesiveMomentHeatError<E>> {
        match self {
            Self::Rolling { contact, state } => Ok(contact
                .source
                .stored_energy_j(&state.source_state())
                .map_err(RigidAdhesiveRollingError::from)?),
            Self::Torsion { contact, state } => Ok(contact
                .source
                .stored_energy_j(&state.source_state())
                .map_err(RigidAdhesiveTorsionError::from)?),
        }
    }
    #[allow(clippy::too_many_arguments)]
    pub(super) fn prepare<E>(
        self,
        bodies: &[RigidBodyState],
        properties: &[MassProperties],
        wrenches: [Wrench; 2],
        dt_s: f64,
        controls: AdhesiveMomentMechanicalControls,
        gate: &CancelGate,
    ) -> Result<NativeStep<'s>, AdhesiveMomentHeatError<E>> {
        Ok(match self {
            Self::Rolling { contact, state } => {
                let next = prepare_rigid_adhesive_rolling_with_wrenches(
                    bodies, properties, &state, contact, wrenches, dt_s, controls, gate,
                )?;
                NativeStep {
                    bodies_after: next.bodies_after,
                    after: Self::Rolling {
                        contact,
                        state: next.state_after,
                    },
                    report: AdhesiveMomentMechanicalReport::Rolling(next.report),
                }
            }
            Self::Torsion { contact, state } => {
                let next = prepare_rigid_adhesive_torsion_with_wrenches(
                    bodies, properties, &state, contact, wrenches, dt_s, controls, gate,
                )?;
                NativeStep {
                    bodies_after: next.bodies_after,
                    after: Self::Torsion {
                        contact,
                        state: next.state_after,
                    },
                    report: AdhesiveMomentMechanicalReport::Torsion(next.report),
                }
            }
        })
    }
}
