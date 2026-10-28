//! Simultaneous retained shear through the original articulated mass operator.
//!
//! All declared two-axis contacts share one frozen configuration and duration.
//! Their complete 2N-row J M^-1 J^T includes cross-contact and cross-axis motion.
//! The shared tribology owner resolves every original constitutive checkpoint
//! against that operator before one accumulated native impulse application.
//! Independent endpoint velocities then recover each physical work port.

use super::super::{ArticulatedTangentialConfig, ArticulatedTangentialDrive, admit};
use super::{ArticulatedPlanarTangentialContact, admit_frame, prepare_axis};
use crate::articulated::FreeFloatingImpulseResponse;
use crate::articulated_contacts::{
    ArticulatedContactError, ArticulatedContactSystem, EventResponse, EventWorkspace, finite, poll,
};
use fs_ga::Vec3;
use fs_math::det;
use fs_tribo::tangential::multiple::{
    MultipleTangentialConfig, PreparedTangentialContacts, TangentialContactInput,
};
use fs_tribo::{TangentialLaw, TangentialStep, TangentialStepInput};

/// One declared interface with its original law, accepted history and conditions.
#[derive(Debug)]
pub struct ArticulatedTangentialSite<'a, L: TangentialLaw> {
    /// Original material points and physical constitutive tangent frame.
    pub contact: ArticulatedPlanarTangentialContact,
    /// Original constitutive owner; coefficients may differ between contacts.
    pub law: &'a L,
    /// Accepted checkpoint borrowed unchanged by every coupled trial.
    pub state: &'a L::State,
    /// Per-interface load and environment; all sites must share exactly one dt.
    pub drive: ArticulatedTangentialDrive,
}

/// Explicit mechanical and collective constitutive work/storage budgets.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ArticulatedMultiTangentialConfig {
    /// Per-contact geometry/force/speed and total native energy allowances.
    pub mechanical: ArticulatedTangentialConfig,
    /// Contact count, dense contact-space size and collective solve controls.
    pub coupling: MultipleTangentialConfig,
}

/// Complete simultaneous body/history candidate in the original contact order.
#[derive(Clone, Debug)]
pub struct ArticulatedMultiTangentialResult<S> {
    /// Original native articulated updates, in supplied system order.
    pub responses: Vec<FreeFloatingImpulseResponse>,
    /// Every original law candidate; publish all histories and body states together.
    pub contacts: Vec<TangentialStep<S>>,
    /// Actual 2N-by-2N row-major J M^-1 J^T [1/kg]. Rows are contact 0 axes 0/1,
    /// then contact 1 axes 0/1, and so on; no cross block is discarded.
    pub mobility_kg_inverse: Vec<f64>,
    /// Original full row-major dt² W / 2 [m/N], moved from the shared owner.
    pub compliance_m_per_n: Vec<f64>,
    /// Original point-minus-counterpart tangent speeds [m/s] at each interface.
    pub relative_speed_before_m_s: Vec<[f64; 2]>,
    /// Speeds reconstructed from the actual native endpoint bodies [m/s].
    pub relative_speed_after_m_s: Vec<[f64; 2]>,
    /// Actual counterpart-minus-point midpoint displacements [m].
    pub slip_increment_m: Vec<[f64; 2]>,
    /// Actual displacement minus original constitutive displacement [m].
    pub compatibility_residual_m: Vec<[f64; 2]>,
    /// Original reversible contact storage before and after [J].
    pub contact_energy_before_j: Vec<f64>,
    pub contact_energy_after_j: Vec<f64>,
    /// Per-interface storage change + heat - F dot actual slip [J].
    pub contact_work_residual_j: Vec<f64>,
    /// Sum of physical exported contact heat [J]; no impact loss is relabeled.
    pub dissipated_energy_j: f64,
    /// Actual total articulated kinetic-energy change [J].
    pub kinetic_energy_change_j: f64,
    /// Signed work supplied by all prescribed supports [J].
    pub support_work_j: f64,
    /// Actual kinetic + all contact storage changes + heat - support work [J].
    pub energy_balance_residual_j: f64,
    /// Actual world momentum changes after external support impulses [N s, N m s].
    pub linear_momentum_residual: Vec3,
    pub angular_momentum_residual: Vec3,
    /// Actual versus full native mobility prediction discrepancy [m/s].
    pub endpoint_velocity_residual_m_s: f64,
    /// Collective solve work, including original-law trials that did not publish.
    pub coupling_sweeps: usize,
    pub law_evaluations: usize,
}

/// Advance all declared retained-shear contacts through one native body update.
///
/// The site count is 1..=the explicit contact cap. The complete (2N)² matrix
/// allocation is admitted before any native probe. Each site retains its own
/// load, thermal/pressure condition, constitutive owner and history. Every site
/// must use the same physical duration and an unchanged physical frame/history
/// association. Shared-body coupling is never replaced by independent contacts.
///
/// One failed source, history, geometry, solve or work check rejects the whole
/// candidate. Cancellation is observed during assembly, original-law evaluations,
/// native impulse application and final actual-state checks. Accepted mechanics
/// and all histories are borrowed throughout. The caller owns normal response,
/// pose integration, contact discovery, interface identity/ownership, frame
/// transport and time refinement across intrastep branch changes. This shear
/// operation introduces no normal constraint, spin/rolling law or thermal store.
pub fn advance_articulated_tangential_contacts<L: TangentialLaw>(
    systems: &[ArticulatedContactSystem<'_>],
    sites: &[ArticulatedTangentialSite<'_, L>],
    config: ArticulatedMultiTangentialConfig,
    mut cancelled: impl FnMut() -> bool,
) -> Result<ArticulatedMultiTangentialResult<L::State>, ArticulatedContactError> {
    poll(&mut cancelled)?;
    config.coupling.validate()?;
    let count = sites.len();
    if count == 0 || count > config.coupling.maximum_contacts {
        return Err(ArticulatedContactError::Invalid(
            "simultaneous shear contact count exceeds its explicit domain",
        ));
    }
    // The shared cap is at most 128 contacts. Check before constructing any
    // square matrix or invoking a native inverse-mass probe.
    let dimension = 2 * count;
    if dimension * dimension > config.coupling.maximum_matrix_terms {
        return Err(ArticulatedContactError::Invalid(
            "simultaneous shear contact-space matrix exceeds its allocation budget",
        ));
    }
    let mechanical = config.mechanical;
    let dt = sites[0].drive.dt_s;
    let mut energies_before = Vec::with_capacity(count);
    for site in sites {
        poll(&mut cancelled)?;
        admit(systems.len(), site.drive, mechanical)?;
        if site.drive.dt_s != dt {
            return Err(ArticulatedContactError::Invalid(
                "simultaneous shear contacts require one physical interval",
            ));
        }
        let energy = finite(site.law.stored_energy_j(site.state))?;
        if energy < 0.0 {
            return Err(ArticulatedContactError::Invalid("negative contact storage"));
        }
        energies_before.push(energy);
    }
    let mut workspace = EventWorkspace::new(systems, &mut cancelled)?;
    let mut rows = Vec::with_capacity(dimension);
    let mut before_speed = Vec::with_capacity(count);
    let mut inputs = Vec::with_capacity(count);
    for site in sites {
        poll(&mut cancelled)?;
        let axes = admit_frame(site.contact.frame_world)?;
        let mut speed = [0.0; 2];
        for i in 0..2 {
            let row = prepare_axis(
                site.contact,
                axes[i],
                &workspace,
                mechanical.geometry_tolerance_m,
            )?;
            speed[i] = finite(row.velocity(&workspace.before))?;
            rows.push(row);
        }
        before_speed.push(speed);
        inputs.push(TangentialContactInput {
            law: site.law,
            state: site.state,
            input: TangentialStepInput {
                normal_n: site.drive.normal_force_n,
                slip_increment_m: [finite(-dt * speed[0])?, finite(-dt * speed[1])?],
                dt_s: dt,
                contact_temperature_k: site.drive.contact_temperature_k,
                nominal_pressure_pa: site.drive.nominal_pressure_pa,
                compliance_m_per_n: 0.0,
                solve: mechanical.solve,
            },
        });
    }
    let mobility = workspace.delassus(systems, &rows, &mut cancelled)?;
    let compliance = mobility
        .iter()
        .map(|w| finite(0.5 * dt * dt * w))
        .collect::<Result<Vec<_>, _>>()?;
    let prepared =
        PreparedTangentialContacts::new(count, compliance, config.coupling, &mut cancelled)?;
    let coupled = prepared.advance(&inputs, &mut cancelled)?;
    poll(&mut cancelled)?;
    let mut energies_after = Vec::with_capacity(count);
    let mut impulses = Vec::with_capacity(dimension);
    let mut stored_change = 0.0;
    let mut heat = 0.0;
    for (index, (site, step)) in sites.iter().zip(&coupled.contacts).enumerate() {
        poll(&mut cancelled)?;
        let energy = finite(site.law.stored_energy_j(&step.state))?;
        let change = finite(energy - energies_before[index])?;
        let dissipated = finite(step.dissipated_j)?;
        if energy < 0.0 || dissipated < 0.0 {
            return Err(ArticulatedContactError::Invalid(
                "negative contact storage or heat",
            ));
        }
        let force = finite(det::hypot(step.traction_n[0], step.traction_n[1]))?;
        if force > mechanical.solve.maximum_force_n
            || step.iterations > mechanical.solve.max_iterations
        {
            return Err(ArticulatedContactError::Unresolved(
                "constitutive effort/work cap exceeded",
            ));
        }
        for i in 0..2 {
            if finite(step.compatibility_residual_m[i])?.abs() > mechanical.displacement_tolerance_m
            {
                return Err(ArticulatedContactError::Unresolved(
                    "collective constitutive compatibility",
                ));
            }
            impulses.push(finite(dt * step.traction_n[i])?);
        }
        for residual in [
            finite(change - step.stored_change_j)?,
            finite(step.closure_residual_j)?,
            finite(coupled.contact_work_residuals_j[index])?,
        ] {
            if residual.abs() > mechanical.energy_tolerance_j {
                return Err(ArticulatedContactError::Unresolved(
                    "collective constitutive contact energy",
                ));
            }
        }
        stored_change = finite(stored_change + change)?;
        heat = finite(heat + dissipated)?;
        energies_after.push(energy);
    }
    for residual in [
        finite(stored_change - coupled.stored_change_j)?,
        finite(heat - coupled.dissipated_j)?,
        finite(coupled.body_work_residual_j)?,
    ] {
        if residual.abs() > mechanical.energy_tolerance_j {
            return Err(ArticulatedContactError::Unresolved(
                "collective constitutive energy total",
            ));
        }
    }
    let mut predicted = Vec::with_capacity(dimension);
    for (row, values) in mobility.chunks_exact(dimension).enumerate() {
        let mut speed = before_speed[row / 2][row % 2];
        for (&w, &impulse) in values.iter().zip(&impulses) {
            speed = finite(speed + w * impulse)?;
        }
        predicted.push(speed);
    }
    let mechanics = workspace.apply(
        systems,
        &rows,
        EventResponse {
            impulses: &impulses,
            predicted_speeds: &predicted,
            stored_energy_change: stored_change,
            dissipated_energy: heat,
            velocity_tolerance: mechanical.velocity_tolerance_m_s,
            energy_tolerance: mechanical.energy_tolerance_j,
        },
        &mut cancelled,
    )?;
    let mut after_speed = Vec::with_capacity(count);
    let mut slips = Vec::with_capacity(count);
    let mut compatibility = Vec::with_capacity(count);
    let mut work_residuals = Vec::with_capacity(count);
    for (index, step) in coupled.contacts.iter().enumerate() {
        poll(&mut cancelled)?;
        let after = [
            mechanics.speeds_after[2 * index],
            mechanics.speeds_after[2 * index + 1],
        ];
        let mut slip = [0.0; 2];
        let mut residual = [0.0; 2];
        for i in 0..2 {
            slip[i] = finite(-dt * (0.5 * before_speed[index][i] + 0.5 * after[i]))?;
            residual[i] = finite(slip[i] - step.slip_increment_m[i])?;
            if residual[i].abs() > mechanical.displacement_tolerance_m {
                return Err(ArticulatedContactError::Unresolved(
                    "actual coupled shear displacement",
                ));
            }
        }
        let work = finite(
            energies_after[index] - energies_before[index] + step.dissipated_j
                - step.traction_n[0] * slip[0]
                - step.traction_n[1] * slip[1],
        )?;
        if finite(det::hypot(slip[0], slip[1]) / dt)? > mechanical.maximum_slip_speed_m_s
            || work.abs() > mechanical.energy_tolerance_j
        {
            return Err(ArticulatedContactError::Unresolved(
                "actual coupled shear speed/work",
            ));
        }
        after_speed.push(after);
        slips.push(slip);
        compatibility.push(residual);
        work_residuals.push(work);
    }
    poll(&mut cancelled)?;
    Ok(ArticulatedMultiTangentialResult {
        responses: mechanics.responses,
        contacts: coupled.contacts,
        mobility_kg_inverse: mobility,
        compliance_m_per_n: prepared.into_compliance_m_per_n(),
        relative_speed_before_m_s: before_speed,
        relative_speed_after_m_s: after_speed,
        slip_increment_m: slips,
        compatibility_residual_m: compatibility,
        contact_energy_before_j: energies_before,
        contact_energy_after_j: energies_after,
        contact_work_residual_j: work_residuals,
        dissipated_energy_j: heat,
        kinetic_energy_change_j: mechanics.kinetic_energy_change,
        support_work_j: mechanics.plane_work,
        energy_balance_residual_j: mechanics.energy_balance_residual,
        linear_momentum_residual: mechanics.linear_momentum_residual,
        angular_momentum_residual: mechanics.angular_momentum_residual,
        endpoint_velocity_residual_m_s: mechanics.endpoint_velocity_residual,
        coupling_sweeps: coupled.sweeps,
        law_evaluations: coupled.law_evaluations,
    })
}
