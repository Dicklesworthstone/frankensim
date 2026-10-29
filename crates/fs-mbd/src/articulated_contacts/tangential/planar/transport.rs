//! Objective retained shear on the original articulated contact proposals.
//!
//! Accepted history carries its physical frame. The shared tribology owner
//! transports it by the minimum normal rotation and the explicit integrated
//! common material spin, then the original single or simultaneous mechanical
//! proposal acts at the supplied frozen endpoint configuration. This is a
//! transport-then-shear split, not integration of the moving contact geometry.
//! Tangent gauge changes add no material rotation; common spin adds no slip.
//! Endpoint normals and one spin angle do not reconstruct a general path of
//! noncommuting normal tilt and spin; the caller must resolve that path in time.

use super::super::{ArticulatedTangentialDrive, admit};
use super::multiple::{
    ArticulatedMultiTangentialConfig, ArticulatedMultiTangentialResult, ArticulatedTangentialSite,
    advance_articulated_tangential_contacts,
};
use super::{
    ArticulatedPlanarTangentialConfig, ArticulatedPlanarTangentialContact,
    ArticulatedPlanarTangentialResult, admit_frame, advance_articulated_planar_tangential,
};
use crate::articulated_contacts::{
    ArticulatedContactError, ArticulatedContactSystem, finite, poll,
};
use fs_tribo::TangentialStep;
use fs_tribo::partial_slip::TangentFrame;
use fs_tribo::tangential::transport::{
    FramedTangentialState, ObjectiveTangentialLaw, TangentialTransportConfig,
};

/// Single native body/history candidate, with the accepted frame inseparable
/// from the new constitutive checkpoint at `mechanics.contact.state`.
#[derive(Clone, Debug)]
pub struct ArticulatedFramedPlanarTangentialResult<S> {
    /// Unchanged mechanical shear audit, relative to the transported checkpoint.
    pub mechanics: ArticulatedPlanarTangentialResult<FramedTangentialState<S>>,
    /// Transported minus originally accepted reversible energy [J]. This is
    /// roundoff, never heat, coefficient work or a physical energy source.
    pub transport_energy_residual_j: f64,
    /// Original accepted storage change + heat - actual contact work [J].
    pub accepted_contact_work_residual_j: f64,
    /// Actual kinetic + original accepted storage change + heat - support work [J].
    pub accepted_energy_balance_residual_j: f64,
}

/// One interface at its supplied endpoint configuration and its old framed
/// checkpoint. Different sites may tilt, change tangent gauge and spin differently.
#[derive(Debug)]
pub struct ArticulatedFramedTangentialSite<'a, L: ObjectiveTangentialLaw> {
    pub contact: ArticulatedPlanarTangentialContact,
    pub law: &'a L,
    pub history: &'a FramedTangentialState<L::State>,
    pub drive: ArticulatedTangentialDrive,
    /// Integrated common material-frame rotation about the endpoint normal [rad].
    /// This is not relative spin, slip or an inferred tangent-gauge rotation.
    pub common_spin_angle_rad: f64,
}

/// All body responses and all newly framed histories, in original site order.
#[derive(Clone, Debug)]
pub struct ArticulatedFramedMultiTangentialResult<S> {
    /// Original fully coupled mechanical candidate; each contact state includes
    /// its endpoint frame. No cross-contact mobility is removed by transport.
    pub mechanics: ArticulatedMultiTangentialResult<FramedTangentialState<S>>,
    /// Numerical energy changes from transport, one per contact [J].
    pub transport_energy_residuals_j: Vec<f64>,
    /// Original accepted storage changes + heat - actual per-contact work [J].
    pub accepted_contact_work_residuals_j: Vec<f64>,
    /// Actual kinetic + all original accepted storage changes + heat - support work [J].
    pub accepted_energy_balance_residual_j: f64,
}

/// Transport the same accepted history, then propose one original shear step.
///
/// `site.contact.frame_world` is the endpoint working frame. The shared transport is
/// an isometry using the minimum rotation between the accepted and endpoint
/// normals, followed by the supplied common spin. Antipodal/invalid transport
/// refuses. No angular slip, heat or coefficient change is fabricated. The
/// caller owns geometry/pose evolution, contact identity, the common-spin
/// integral and refinement of this split across finite rotations/reversals.
///
/// Accepted bodies and history are borrowed. All transport, native mechanics,
/// original-state energy checks and cancellation checks precede return of the
/// framed candidate. Transport's tolerance cannot widen the mechanical budget.
pub fn advance_articulated_framed_planar_tangential<L: ObjectiveTangentialLaw>(
    systems: &[ArticulatedContactSystem<'_>],
    site: ArticulatedFramedTangentialSite<'_, L>,
    config: ArticulatedPlanarTangentialConfig,
    transport_config: TangentialTransportConfig,
    mut cancelled: impl FnMut() -> bool,
) -> Result<ArticulatedFramedPlanarTangentialResult<L::State>, ArticulatedContactError> {
    poll(&mut cancelled)?;
    let ArticulatedFramedTangentialSite {
        contact,
        law,
        history,
        drive,
        common_spin_angle_rad,
    } = site;
    admit(systems.len(), drive, config.mechanical)?;
    config.coupling.validate()?;
    transport_config.validate()?;
    admit_frame(contact.frame_world)?;
    let transported = history.transport(
        law,
        contact.frame_world,
        common_spin_angle_rad,
        transport_config,
        &mut cancelled,
    )?;
    let candidate = advance_articulated_planar_tangential(
        systems,
        contact,
        law,
        transported.checkpoint.state(),
        drive,
        config,
        &mut cancelled,
    )?;
    let drift = transported.stored_energy_residual_j;
    let work = accepted_residual(
        candidate.contact_work_residual_j,
        drift,
        config.mechanical.energy_tolerance_j,
    )?;
    let balance = accepted_residual(
        candidate.energy_balance_residual_j,
        drift,
        config.mechanical.energy_tolerance_j,
    )?;
    let mechanics = ArticulatedPlanarTangentialResult {
        contact: frame_step(law, contact.frame_world, candidate.contact)?,
        responses: candidate.responses,
        mobility_kg_inverse: candidate.mobility_kg_inverse,
        compliance_m_per_n: candidate.compliance_m_per_n,
        relative_speed_before_m_s: candidate.relative_speed_before_m_s,
        relative_speed_after_m_s: candidate.relative_speed_after_m_s,
        slip_increment_m: candidate.slip_increment_m,
        compatibility_residual_m: candidate.compatibility_residual_m,
        coupling_iterations: candidate.coupling_iterations,
        law_evaluations: candidate.law_evaluations,
        contact_energy_before_j: candidate.contact_energy_before_j,
        contact_energy_after_j: candidate.contact_energy_after_j,
        kinetic_energy_change_j: candidate.kinetic_energy_change_j,
        support_work_j: candidate.support_work_j,
        contact_work_residual_j: candidate.contact_work_residual_j,
        energy_balance_residual_j: candidate.energy_balance_residual_j,
        linear_momentum_residual: candidate.linear_momentum_residual,
        angular_momentum_residual: candidate.angular_momentum_residual,
        endpoint_velocity_residual_m_s: candidate.endpoint_velocity_residual_m_s,
    };
    poll(&mut cancelled)?;
    Ok(ArticulatedFramedPlanarTangentialResult {
        mechanics,
        transport_energy_residual_j: drift,
        accepted_contact_work_residual_j: work,
        accepted_energy_balance_residual_j: balance,
    })
}

/// Transport every retained history and propose all interfaces simultaneously.
///
/// Count and full-matrix budgets are admitted before allocating transported
/// checkpoints. All sites use one duration; every transport and every original
/// constitutive response must succeed before any body/history candidate is
/// returned. The unchanged native owner assembles the complete contact mobility
/// and accumulates all impulses once. Geometry and time integration remain with
/// the caller under the same split boundary as the single-contact entry point.
pub fn advance_articulated_framed_tangential_contacts<L: ObjectiveTangentialLaw>(
    systems: &[ArticulatedContactSystem<'_>],
    sites: &[ArticulatedFramedTangentialSite<'_, L>],
    config: ArticulatedMultiTangentialConfig,
    transport_config: TangentialTransportConfig,
    mut cancelled: impl FnMut() -> bool,
) -> Result<ArticulatedFramedMultiTangentialResult<L::State>, ArticulatedContactError> {
    poll(&mut cancelled)?;
    config.coupling.validate()?;
    transport_config.validate()?;
    let count = sites.len();
    if count == 0 || count > config.coupling.maximum_contacts {
        return Err(ArticulatedContactError::Invalid(
            "framed shear contact count exceeds its explicit domain",
        ));
    }
    // The shared validated count cap is at most 128, so these products fit.
    if (2 * count) * (2 * count) > config.coupling.maximum_matrix_terms {
        return Err(ArticulatedContactError::Invalid(
            "framed shear contact-space matrix exceeds its allocation budget",
        ));
    }
    let dt = sites[0].drive.dt_s;
    for site in sites {
        poll(&mut cancelled)?;
        admit(systems.len(), site.drive, config.mechanical)?;
        if site.drive.dt_s != dt {
            return Err(ArticulatedContactError::Invalid(
                "framed shear contacts require one physical interval",
            ));
        }
        admit_frame(site.contact.frame_world)?;
    }
    let mut transported = Vec::with_capacity(count);
    let mut drifts = Vec::with_capacity(count);
    let mut total_drift = 0.0;
    for site in sites {
        poll(&mut cancelled)?;
        let step = site.history.transport(
            site.law,
            site.contact.frame_world,
            site.common_spin_angle_rad,
            transport_config,
            &mut cancelled,
        )?;
        total_drift = finite(total_drift + step.stored_energy_residual_j)?;
        drifts.push(step.stored_energy_residual_j);
        transported.push(step.checkpoint);
    }
    let working: Vec<_> = sites
        .iter()
        .zip(&transported)
        .map(|(site, history)| ArticulatedTangentialSite {
            contact: site.contact,
            law: site.law,
            state: history.state(),
            drive: site.drive,
        })
        .collect();
    let candidate =
        advance_articulated_tangential_contacts(systems, &working, config, &mut cancelled)?;
    let balance = accepted_residual(
        candidate.energy_balance_residual_j,
        total_drift,
        config.mechanical.energy_tolerance_j,
    )?;
    let mut work = Vec::with_capacity(count);
    let mut contacts = Vec::with_capacity(count);
    for (index, (site, step)) in sites.iter().zip(candidate.contacts).enumerate() {
        poll(&mut cancelled)?;
        work.push(accepted_residual(
            candidate.contact_work_residual_j[index],
            drifts[index],
            config.mechanical.energy_tolerance_j,
        )?);
        contacts.push(frame_step(site.law, site.contact.frame_world, step)?);
    }
    let mechanics = ArticulatedMultiTangentialResult {
        contacts,
        responses: candidate.responses,
        mobility_kg_inverse: candidate.mobility_kg_inverse,
        compliance_m_per_n: candidate.compliance_m_per_n,
        relative_speed_before_m_s: candidate.relative_speed_before_m_s,
        relative_speed_after_m_s: candidate.relative_speed_after_m_s,
        slip_increment_m: candidate.slip_increment_m,
        compatibility_residual_m: candidate.compatibility_residual_m,
        contact_energy_before_j: candidate.contact_energy_before_j,
        contact_energy_after_j: candidate.contact_energy_after_j,
        contact_work_residual_j: candidate.contact_work_residual_j,
        dissipated_energy_j: candidate.dissipated_energy_j,
        kinetic_energy_change_j: candidate.kinetic_energy_change_j,
        support_work_j: candidate.support_work_j,
        energy_balance_residual_j: candidate.energy_balance_residual_j,
        linear_momentum_residual: candidate.linear_momentum_residual,
        angular_momentum_residual: candidate.angular_momentum_residual,
        endpoint_velocity_residual_m_s: candidate.endpoint_velocity_residual_m_s,
        coupling_sweeps: candidate.coupling_sweeps,
        law_evaluations: candidate.law_evaluations,
    };
    poll(&mut cancelled)?;
    Ok(ArticulatedFramedMultiTangentialResult {
        mechanics,
        transport_energy_residuals_j: drifts,
        accepted_contact_work_residuals_j: work,
        accepted_energy_balance_residual_j: balance,
    })
}

fn accepted_residual(
    mechanical: f64,
    transport: f64,
    tolerance: f64,
) -> Result<f64, ArticulatedContactError> {
    let residual = finite(mechanical + transport)?;
    if residual.abs() > tolerance {
        return Err(ArticulatedContactError::Unresolved(
            "transported shear violates the original accepted-state work budget",
        ));
    }
    Ok(residual)
}

// Move the one original candidate into its new frame after all mechanics.
// This performs no constitutive evaluation and creates no second history copy.
fn frame_step<L: ObjectiveTangentialLaw>(
    law: &L,
    frame: TangentFrame,
    step: TangentialStep<L::State>,
) -> Result<TangentialStep<FramedTangentialState<L::State>>, ArticulatedContactError> {
    Ok(TangentialStep {
        state: FramedTangentialState::new(law, step.state, frame)?,
        traction_n: step.traction_n,
        slip_increment_m: step.slip_increment_m,
        regime: step.regime,
        stored_change_j: step.stored_change_j,
        dissipated_j: step.dissipated_j,
        closure_residual_j: step.closure_residual_j,
        compatibility_residual_m: step.compatibility_residual_m,
        iterations: step.iterations,
    })
}
