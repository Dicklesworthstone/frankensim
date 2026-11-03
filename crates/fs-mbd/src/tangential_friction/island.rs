//! Native shared-body mechanics, independent of the circular traction law.
//!
//! Sparse inverse-mass Gram rows retain complete inertias and physical arms.
//! All patch impulses sum before one update of each body, then actual native
//! endpoint velocities independently supply work, recoil and event checks.

use super::hertz::multiple::{RigidHertzIslandContactReport, RigidHertzIslandReport};
use super::hertz::{
    RigidHertzControls, RigidHertzError, finite, points, relative_point, relative_spin,
};
use super::{NativeImpulsePhase, ResultantImpulse, apply_resultants, components, magnitude, poll};
use crate::{MassProperties, PointKinematics, RigidBodyState, Vec3};
use fs_exec::CancelGate;
use fs_tribo::sliding_patch::contact::multiple::MultipleHertzStep;
use fs_tribo::sliding_patch::contact::{HertzWrenchMode, HertzWrenchStep};

#[derive(Clone, Debug)]
pub(super) struct PatchMechanics {
    pub(super) bodies: [usize; 2],
    pub(super) arms: [Vec3; 2],
    pub(super) axes: [Vec3; 3],
    pub(super) before: [PointKinematics; 2],
    pub(super) rates_before: [f64; 3],
}

pub(super) type PatchFactor = [[[f64; 6]; 3]; 2];

pub(super) fn factor(
    pair: [RigidBodyState; 2],
    properties: [MassProperties; 2],
    mechanics: &PatchMechanics,
    radius: f64,
) -> PatchFactor {
    let mut factor = [[[0.0; 6]; 3]; 2];
    for (side, sign) in [(0, 1.0), (1, -1.0)] {
        let inverse_mass_root = properties[side].mass().sqrt().recip();
        let orientation = pair[side].pose().orientation();
        for axis in 0..3 {
            let direction = orientation.rotate_world_to_body(mechanics.axes[axis]);
            let translation = if axis < 2 {
                components(mechanics.axes[axis].scale(sign * inverse_mass_root))
            } else {
                [0.0; 3]
            };
            let torque = if axis < 2 {
                mechanics.arms[side].cross(direction).scale(sign)
            } else {
                direction.scale(sign * radius)
            };
            let rotation = components(properties[side].angular_mobility_factor(torque));
            factor[side][axis][..3].copy_from_slice(&translation);
            factor[side][axis][3..].copy_from_slice(&rotation);
        }
    }
    factor
}

pub(super) fn mobility(
    mechanics: &[PatchMechanics],
    factors: &[PatchFactor],
    gate: &CancelGate,
) -> Result<Vec<f64>, RigidHertzError> {
    if factors.len() != mechanics.len() {
        return Err(RigidHertzError::Invalid(
            "one original factor per contact required",
        ));
    }
    let rows = 3 * mechanics.len();
    let mut mobility = vec![0.0; rows * rows];
    for i in 0..rows {
        poll(gate)?;
        for j in 0..=i {
            let mut value = 0.0;
            for side_i in 0..2 {
                for side_j in 0..2 {
                    if mechanics[i / 3].bodies[side_i] == mechanics[j / 3].bodies[side_j] {
                        value += factors[i / 3][side_i][i % 3]
                            .iter()
                            .zip(&factors[j / 3][side_j][j % 3])
                            .map(|(a, b)| a * b)
                            .sum::<f64>();
                    }
                }
            }
            let value = finite(value)?;
            mobility[i * rows + j] = value;
            mobility[j * rows + i] = value;
        }
    }
    Ok(mobility)
}

pub(super) fn apply_wrenches(
    before: &[RigidBodyState],
    properties: &[MassProperties],
    mechanics: &[PatchMechanics],
    contacts: &[HertzWrenchStep],
    dt_s: f64,
    gate: &CancelGate,
) -> Result<NativeImpulsePhase, RigidHertzError> {
    if mechanics.len() != contacts.len() {
        return Err(RigidHertzError::Invalid(
            "one original wrench per contact required",
        ));
    }
    let mut resultants = vec![ResultantImpulse::default(); before.len()];
    for (mechanics, contact) in mechanics.iter().zip(contacts) {
        poll(gate)?;
        let impulse = mechanics.axes[0]
            .scale(dt_s * contact.force_n[0])
            .add(mechanics.axes[1].scale(dt_s * contact.force_n[1]));
        let couple = mechanics.axes[2].scale(dt_s * contact.spin_moment_nm);
        for (side, sign) in [(0, 1.0), (1, -1.0)] {
            let body = mechanics.bodies[side];
            resultants[body].add(
                before[body],
                mechanics.arms[side],
                impulse.scale(sign),
                couple.scale(sign),
            );
        }
    }
    Ok(apply_resultants(before, properties, &resultants, gate)?)
}

/// Original law evaluated independently at actual native midpoint motion.
pub(super) struct NativeTraction {
    pub(super) force_n: [f64; 2],
    pub(super) spin_moment_nm: f64,
    pub(super) dissipated_power_w: f64,
}

#[allow(clippy::too_many_arguments)]
pub(super) fn contact_report(
    patch_id: &str,
    mechanics: &PatchMechanics,
    after: &[RigidBodyState],
    properties: &[MassProperties],
    radius: f64,
    dt_s: f64,
    contact: &HertzWrenchStep,
    native: impl FnOnce([f64; 2], f64) -> Result<NativeTraction, RigidHertzError>,
) -> Result<RigidHertzIslandContactReport, RigidHertzError> {
    let current = points(
        [after[mechanics.bodies[0]], after[mechanics.bodies[1]]],
        [
            properties[mechanics.bodies[0]],
            properties[mechanics.bodies[1]],
        ],
        mechanics.arms,
    )?;
    let relative_before = relative_point(mechanics.before);
    let relative_after = relative_point(current);
    let spin_before = relative_spin(mechanics.before);
    let spin_after = relative_spin(current);
    let midpoint = relative_before.add(relative_after).scale(0.5);
    let midpoint_spin = spin_before
        .add(spin_after)
        .scale(0.5)
        .dot(mechanics.axes[2]);
    let force = mechanics.axes[0]
        .scale(contact.force_n[0])
        .add(mechanics.axes[1].scale(contact.force_n[1]));
    let couple = mechanics.axes[2].scale(contact.spin_moment_nm);
    let impulse = force.scale(dt_s);
    let angular_impulse = couple.scale(dt_s);
    let midpoint_work = impulse.dot(midpoint) + dt_s * contact.spin_moment_nm * midpoint_spin;
    let slip_after = [
        relative_after.dot(mechanics.axes[0]),
        relative_after.dot(mechanics.axes[1]),
    ];
    let mut heat = 0.0;
    let mut force_residual = 0.0;
    let mut moment_residual = 0.0;
    if contact.mode != HertzWrenchMode::Sticking {
        let response = native(
            [
                midpoint.dot(mechanics.axes[0]),
                midpoint.dot(mechanics.axes[1]),
            ],
            midpoint_spin,
        )?;
        heat = finite(dt_s * response.dissipated_power_w)?;
        force_residual = (response.force_n[0] - contact.force_n[0])
            .hypot(response.force_n[1] - contact.force_n[1]);
        moment_residual = (response.spin_moment_nm - contact.spin_moment_nm).abs();
    }
    let endpoint_velocity_residual = (slip_after[0]
        - (2.0 * contact.slip_displacement_m[0] / dt_s - mechanics.rates_before[0]))
        .abs()
        .max(
            (slip_after[1]
                - (2.0 * contact.slip_displacement_m[1] / dt_s - mechanics.rates_before[1]))
                .abs(),
        );
    let endpoint_spin_residual = (spin_after.dot(mechanics.axes[2])
        - (2.0 * contact.spin_angle_rad / dt_s - mechanics.rates_before[2] / radius))
        .abs();
    let assisting_work = [
        force.dot(relative_before),
        force.dot(relative_after),
        couple.dot(spin_before),
        couple.dot(spin_after),
    ]
    .into_iter()
    .fold(0.0, f64::max)
        * dt_s;
    Ok(RigidHertzIslandContactReport {
        patch_id: patch_id.to_owned(),
        contact: *contact,
        slip_velocity_before_m_s: [mechanics.rates_before[0], mechanics.rates_before[1]],
        slip_velocity_after_m_s: slip_after,
        spin_before_rad_s: spin_before.dot(mechanics.axes[2]),
        spin_after_rad_s: spin_after.dot(mechanics.axes[2]),
        normal_velocity_before_m_s: relative_before.dot(mechanics.axes[2]),
        normal_velocity_after_m_s: relative_after.dot(mechanics.axes[2]),
        actual_midpoint_displacement_m: [
            dt_s * midpoint.dot(mechanics.axes[0]),
            dt_s * midpoint.dot(mechanics.axes[1]),
        ],
        actual_midpoint_spin_angle_rad: dt_s * midpoint_spin,
        force_impulse_world_ns: impulse,
        spin_angular_impulse_world_nms: angular_impulse,
        midpoint_impulse_work_j: midpoint_work,
        dissipated_energy_j: heat,
        contact_work_residual_j: midpoint_work + contact.contact_work_j,
        contact_heat_residual_j: heat - contact.dissipated_energy_j,
        native_force_residual_n: force_residual,
        native_spin_moment_residual_nm: moment_residual,
        endpoint_velocity_residual_m_s: endpoint_velocity_residual,
        endpoint_spin_residual_rad_s: endpoint_spin_residual,
        maximum_assisting_work_j: assisting_work,
    })
}

pub(super) fn contact_admitted(
    report: &RigidHertzIslandContactReport,
    controls: RigidHertzControls,
    held_normal: bool,
) -> bool {
    [
        report.midpoint_impulse_work_j,
        report.dissipated_energy_j,
        report.native_force_residual_n,
        report.native_spin_moment_residual_nm,
        report.endpoint_velocity_residual_m_s,
        report.endpoint_spin_residual_rad_s,
        report.maximum_assisting_work_j,
        report.normal_velocity_after_m_s,
    ]
    .into_iter()
    .all(f64::is_finite)
        && [
            report.contact_work_residual_j,
            report.contact_heat_residual_j,
        ]
        .into_iter()
        .all(|x| x.is_finite() && x.abs() <= controls.energy_tolerance_j)
        && report.native_force_residual_n <= controls.contact.force_tolerance_n
        && report.native_spin_moment_residual_nm <= controls.contact.moment_tolerance_nm
        && report.endpoint_velocity_residual_m_s <= controls.velocity_tolerance_m_s
        && report.endpoint_spin_residual_rad_s <= controls.angular_velocity_tolerance_rad_s
        && (!held_normal
            || report.normal_velocity_after_m_s.abs() <= controls.velocity_tolerance_m_s)
}

pub(super) fn report(
    contacts: Vec<RigidHertzIslandContactReport>,
    native: MultipleHertzStep,
    impulses: &NativeImpulsePhase,
    gate: &CancelGate,
) -> Result<RigidHertzIslandReport, RigidHertzError> {
    poll(gate)?;
    let midpoint_work = contacts
        .iter()
        .map(|r| r.midpoint_impulse_work_j)
        .sum::<f64>();
    let heat = contacts.iter().map(|r| r.dissipated_energy_j).sum::<f64>();
    Ok(RigidHertzIslandReport {
        contacts,
        native,
        kinetic_before_j: impulses.kinetic_before_j,
        kinetic_after_j: impulses.kinetic_after_j,
        kinetic_change_j: impulses.kinetic_change_j,
        midpoint_impulse_work_j: midpoint_work,
        dissipated_energy_j: heat,
        work_energy_residual_j: impulses.kinetic_change_j - midpoint_work,
        total_energy_residual_j: impulses.kinetic_change_j + heat,
        linear_momentum_change_ns: impulses.linear_momentum_change_ns,
        angular_momentum_change_nms: impulses.angular_momentum_change_nms,
    })
}

pub(super) fn report_admitted(
    report: &RigidHertzIslandReport,
    controls: RigidHertzControls,
) -> bool {
    [
        report.kinetic_before_j,
        report.kinetic_after_j,
        report.midpoint_impulse_work_j,
        report.dissipated_energy_j,
        magnitude(report.linear_momentum_change_ns),
        magnitude(report.angular_momentum_change_nms),
    ]
    .into_iter()
    .all(f64::is_finite)
        && [
            report.work_energy_residual_j,
            report.total_energy_residual_j,
            report.dissipated_energy_j - report.native.dissipated_energy_j,
        ]
        .into_iter()
        .all(|x| x.is_finite() && x.abs() <= controls.energy_tolerance_j)
        && magnitude(report.linear_momentum_change_ns) <= controls.linear_momentum_tolerance_ns
        && magnitude(report.angular_momentum_change_nms) <= controls.angular_momentum_tolerance_nms
}
