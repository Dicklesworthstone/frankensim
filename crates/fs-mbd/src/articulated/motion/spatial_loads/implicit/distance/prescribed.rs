//! Native mixed material/joint impulse events on a fixed or driven support.
//!
//! The original prescribed-base ABA supplies every inverse-mass column and the
//! final support reaction. Material geometry, row normalization, rank admission,
//! and factorization are shared with the free mechanism. Scalar root joints are
//! retained; no floating surrogate support or force-duration conversion is used.

use super::*;
use crate::articulated::motion::prescribed::PrescribedBaseMotionState;
use crate::articulated::{
    PrescribedBaseImpulseResponse, apply_prescribed_base_impulses, prescribed_base_impulse_delta,
};

/// One combined event with the complete native support impulse and signed work.
///
/// `source_work_j`, `constraint_work_j`, and `work_partition_residual_j` partition
/// only `impulse.applied_impulse_work`. The separate `impulse.support_work` is
/// paired with the unchanged prescribed base twist and remains a physical source.
pub type PrescribedBilateralImpulseResponse =
    BilateralImpulseResponse<PrescribedBaseImpulseResponse>;

/// Apply supplied body/joint impulses while existing material and joint rows react.
///
/// All rows participate in one coupled native mobility solve, including spherical
/// closure couples at finite admitted gaps and body/joint cross responses. The
/// source probe is an inertia action, with no endpoint speed admission. Only the
/// actual summed event receives native endpoint and representability checks.
/// Its six-component support impulse is exerted on the tree, including reactions
/// transmitted through a mobile root joint. The base boundary and positions are
/// unchanged; its declared smooth acceleration contributes no finite impulse.
///
/// Initial positions and rates must meet the declared tolerances. A zero source
/// still admits geometry and combined rank, but applies zero reactions and does
/// not repair finite admitted defects. Empty constraints delegate to the native
/// prescribed event. At most 32 independent scalar rows and 4096 links are
/// admitted. This advances no time and performs no contact activation or repair.
#[allow(clippy::too_many_arguments, clippy::too_many_lines)]
pub fn apply_impulses(
    model: &ArticulatedModel,
    state: &PrescribedBaseMotionState,
    joint_impulses: &[f64],
    body_impulses: &[Wrench],
    loops: &[BilateralLoop],
    joint_rows: &[JointConstraint],
    controls: BilateralImpulseConfig,
    mut cancelled: impl FnMut() -> bool,
) -> Result<PrescribedBilateralImpulseResponse, BilateralLoopError> {
    poll(&mut cancelled)?;
    let config = BilateralLoopConfig {
        position_tolerance_m: controls.position_tolerance_m,
        velocity_tolerance_mps: controls.velocity_tolerance_mps,
        minimum_scaled_pivot: controls.minimum_scaled_pivot,
        ..Default::default()
    };
    validate_bilateral_controls(config)?;
    if !state.time_s.is_finite() {
        return Err(BilateralLoopError::Invalid("finite event time required"));
    }
    if loops.len() > MAX_BILATERAL_ROWS
        || joint_rows.len() > MAX_BILATERAL_ROWS
        || loops.iter().map(|c| c.rows()).sum::<usize>() + joint_rows.len() > MAX_BILATERAL_ROWS
    {
        return Err(BilateralLoopError::Invalid(
            "at most 32 independent bilateral rows",
        ));
    }
    if model.link_count() > MAX_FREE_FLOATING_STEP_LINKS {
        return Err(FreeFloatingMotionError::TooManyLinks {
            got: model.link_count(),
            maximum: MAX_FREE_FLOATING_STEP_LINKS,
        }
        .into());
    }
    if loops.is_empty() && joint_rows.is_empty() {
        let impulse = apply_prescribed_base_impulses(
            model,
            state.base,
            &state.joint_position,
            &state.joint_velocity,
            joint_impulses,
            body_impulses,
        )?;
        poll(&mut cancelled)?;
        return Ok(BilateralImpulseResponse {
            source_work_j: impulse.applied_impulse_work,
            impulse,
            rows: Vec::new(),
            loop_multipliers_ns: Vec::new(),
            joint_multipliers: Vec::new(),
            constraint_body_impulses: vec![Wrench::default(); model.link_count()],
            constraint_joint_impulses: vec![0.; model.dof_count()],
            position_residual_m: Vec::new(),
            velocity_residual_before_mps: Vec::new(),
            velocity_residual_after_mps: Vec::new(),
            point_residual_before: Vec::new(),
            point_residual_after: Vec::new(),
            joint_position_residuals: Vec::new(),
            joint_velocity_residuals_before: Vec::new(),
            joint_velocity_residuals_after: Vec::new(),
            minimum_scaled_pivot: 1.,
            reciprocal_condition: 1.,
            constraint_work_j: 0.,
            row_constraint_work_j: 0.,
            work_partition_residual_j: 0.,
        });
    }
    model.validate_generalized("joint_impulses", joint_impulses)?;
    crate::articulated::require_len("body_impulses", body_impulses.len(), model.link_count())?;
    crate::articulated::validate_twist(state.base.twist_body, "base.twist_body")?;
    crate::articulated::validate_twist(state.base.acceleration_body, "base.acceleration_body")?;
    for (index, impulse) in body_impulses.iter().enumerate() {
        if index % 64 == 0 {
            poll(&mut cancelled)?;
        }
        crate::articulated::validate_wrench(*impulse, "body_impulses")?;
    }
    admit_loops(model, loops)?;
    let joint_before = admit_constraints(
        joint_rows,
        &state.joint_position,
        &state.joint_velocity,
        &mut cancelled,
    )?;
    // Relative shape kinematics retain mobile root joints and avoid subtracting
    // a common world translation or velocity when constructing material rows.
    let shape = forward_kinematics(
        model,
        BaseState::stationary(Se3::identity()),
        &state.joint_position,
        &state.joint_velocity,
    )?;
    let base_view = FreeFloatingBaseState::new(state.base.world_from_base, state.base.twist_body);
    let before = geometry(&shape, base_view, loops)?;
    check_position(&before, config, "initial")?;
    check_velocity(&before, config, "initial")?;
    let mobility = prescribed_mobility(
        model,
        &state.joint_position,
        &before,
        joint_rows,
        config,
        &mut cancelled,
    )?;
    let reciprocal_condition = mobility.reciprocal_condition(&mut cancelled)?;
    let nonzero_source = joint_impulses.iter().any(|j| *j != 0.)
        || body_impulses
            .iter()
            .any(|j| j.to_array().into_iter().any(|x| x != 0.));
    let loop_count = before.rates.len();
    let mut multipliers = vec![0.; loop_count + joint_rows.len()];
    if nonzero_source {
        poll(&mut cancelled)?;
        let source = prescribed_base_impulse_delta(
            model,
            &state.joint_position,
            joint_impulses,
            body_impulses,
        )?
        .velocity;
        let mut rates = before.rates.clone();
        for (rate, terms) in rates.iter_mut().zip(&before.terms) {
            poll(&mut cancelled)?;
            *rate = (Dd::from_f64(*rate)
                + sum_products(
                    terms.iter().flat_map(|(_, wrench)| wrench.to_array()),
                    terms
                        .iter()
                        .flat_map(|(link, _)| source.body[*link].to_array()),
                ))
            .to_f64();
        }
        let joint_change = apply_rows(joint_rows, &source.joint, &mut cancelled)?;
        rates.extend(
            joint_before
                .velocity
                .iter()
                .zip(joint_change)
                .map(|(v, dv)| (Dd::from_f64(*v) + Dd::from_f64(dv)).to_f64()),
        );
        validate_slice(&rates, "bilateral_impulse.source_rates")?;
        multipliers = mobility.cancel_rates(&rates)?;
    }
    let (constraint_body_impulses, constraint_joint_impulses) =
        assemble_impulses(model, &before, joint_rows, &multipliers, &mut cancelled)?;
    let total_joint: Vec<_> = joint_impulses
        .iter()
        .zip(&constraint_joint_impulses)
        .map(|(source, reaction)| source + reaction)
        .collect();
    let mut total_body = Vec::with_capacity(model.link_count());
    for (index, (source, reaction)) in body_impulses
        .iter()
        .zip(&constraint_body_impulses)
        .enumerate()
    {
        if index % 64 == 0 {
            poll(&mut cancelled)?;
        }
        total_body.push(Wrench::new(
            source.torque + reaction.torque,
            source.force + reaction.force,
        ));
    }
    poll(&mut cancelled)?;
    let impulse = apply_prescribed_base_impulses(
        model,
        state.base,
        &state.joint_position,
        &state.joint_velocity,
        &total_joint,
        &total_body,
    )?;
    poll(&mut cancelled)?;
    let after_shape = forward_kinematics(
        model,
        BaseState::stationary(Se3::identity()),
        &state.joint_position,
        &impulse.joint_velocity_after,
    )?;
    let after = geometry(&after_shape, base_view, loops)?;
    check_position(&after, config, "endpoint")?;
    check_velocity(&after, config, "endpoint")?;
    let joint_after = admit_constraints(
        joint_rows,
        &state.joint_position,
        &impulse.joint_velocity_after,
        &mut cancelled,
    )?;
    let before_kin = forward_kinematics(
        model,
        state.base,
        &state.joint_position,
        &state.joint_velocity,
    )?;
    let after_kin = forward_kinematics(
        model,
        impulse.base_after,
        &state.joint_position,
        &impulse.joint_velocity_after,
    )?;
    let source_work = event_work(
        joint_impulses,
        body_impulses,
        &state.joint_velocity,
        &impulse.joint_velocity_after,
        &before_kin,
        &after_kin,
        &mut cancelled,
    )?;
    let constraint_work = event_work(
        &constraint_joint_impulses,
        &constraint_body_impulses,
        &state.joint_velocity,
        &impulse.joint_velocity_after,
        &before_kin,
        &after_kin,
        &mut cancelled,
    )?;
    let row_before = before.rates.iter().chain(&joint_before.velocity).copied();
    let row_after = after.rates.iter().chain(&joint_after.velocity).copied();
    let row_constraint_work_j = multipliers
        .iter()
        .zip(row_before.zip(row_after))
        .fold(Dd::ZERO, |sum, (lambda, (old, new))| {
            sum + Dd::from_f64(*lambda) * (Dd::from_f64(0.5 * old) + Dd::from_f64(0.5 * new))
        })
        .to_f64();
    let work_partition_residual_j =
        (Dd::from_f64(impulse.applied_impulse_work) - source_work - constraint_work).to_f64();
    let source_work_j = source_work.to_f64();
    let constraint_work_j = constraint_work.to_f64();
    validate_slice(
        &[
            source_work_j,
            constraint_work_j,
            row_constraint_work_j,
            work_partition_residual_j,
        ],
        "bilateral_impulse.work",
    )?;
    let joint_multipliers = multipliers.split_off(loop_count);
    poll(&mut cancelled)?;
    Ok(BilateralImpulseResponse {
        impulse,
        rows: before.rows,
        loop_multipliers_ns: multipliers,
        joint_multipliers,
        constraint_body_impulses,
        constraint_joint_impulses,
        position_residual_m: before.gaps,
        velocity_residual_before_mps: before.rates,
        velocity_residual_after_mps: after.rates,
        point_residual_before: before.points.iter().map(|p| p.residual).collect(),
        point_residual_after: after.points.iter().map(|p| p.residual).collect(),
        joint_position_residuals: joint_before.position,
        joint_velocity_residuals_before: joint_before.velocity,
        joint_velocity_residuals_after: joint_after.velocity,
        minimum_scaled_pivot: mobility.minimum_pivot,
        reciprocal_condition,
        source_work_j,
        constraint_work_j,
        row_constraint_work_j,
        work_partition_residual_j,
    })
}

fn prescribed_mobility(
    model: &ArticulatedModel,
    joint_position: &[f64],
    geometry: &Geometry,
    joint_rows: &[JointConstraint],
    config: BilateralLoopConfig,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Mobility, BilateralLoopError> {
    let loop_count = geometry.terms.len();
    let k = loop_count + joint_rows.len();
    let mut matrix = vec![0.; k * k];
    let mut joints = vec![0.; model.dof_count()];
    let mut body = vec![Wrench::default(); model.link_count()];
    for column in 0..k {
        poll(cancelled)?;
        let terms = geometry.terms.get(column);
        if let Some(terms) = terms {
            for (link, wrench) in terms {
                body[*link] = *wrench;
            }
        } else {
            joints.copy_from_slice(&joint_rows[column - loop_count].coefficients);
        }
        let response =
            prescribed_base_impulse_delta(model, joint_position, &joints, &body)?.velocity;
        for (row, row_terms) in geometry.terms.iter().enumerate() {
            matrix[row * k + column] = row_terms
                .iter()
                .map(|(link, w)| w.pairing(response.body[*link]))
                .sum();
        }
        for (row, value) in apply_rows(joint_rows, &response.joint, cancelled)?
            .into_iter()
            .enumerate()
        {
            matrix[(loop_count + row) * k + column] = value;
        }
        if let Some(terms) = terms {
            for (link, _) in terms {
                body[*link] = Wrench::default();
            }
        } else {
            joints.fill(0.);
        }
    }
    Mobility::from_matrix(
        model,
        geometry,
        config,
        "initial",
        matrix,
        k,
        Vec::new(),
        cancelled,
    )
}
