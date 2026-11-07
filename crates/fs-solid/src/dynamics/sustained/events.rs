//! Located Coulomb stops on a rigid contact trace of the original solid.
//!
//! The complete solid may retain elastic and viscous spectator motion. Each
//! contact mass-response column must lie in the original K and C nullspaces;
//! the contact velocity is then affine under a held force. This permits an
//! analytic stopping time and exact force moments without a new exponential,
//! a velocity reset, or freezing the remaining native coordinates.

use super::*;

mod audit;
mod phase;

/// One fixed sliding axis in the contact's orthonormal tangent plane.
/// The axis is normalized internally; its sign does not change the physics.
#[derive(Clone, Copy, Debug)]
pub struct CollinearCoulombContact<'a> {
    /// `friction_coefficient` is the static coefficient. Initial gap must be zero.
    pub contact: SolidSustainedContact<'a>,
    pub kinetic_coefficient: f64,
    pub tangent_axis: [f64; 2],
}

/// Physical mode over a segment with constant force and affine contact velocity.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum CoulombPhaseMode {
    Sticking,
    Sliding,
}

/// An actual segment of the original interval; no state is published at its end.
#[derive(Debug)]
pub struct SolidCoulombPhase {
    pub mode: CoulombPhaseMode,
    pub start_offset_s: f64,
    /// Rounded duration; integration retains the full difference between the
    /// two representable segment endpoints, including its low component.
    pub duration_s: f64,
    pub friction_coefficient: f64,
    pub force_n: [f64; 3],
    /// Normal and two tangent components, relative to the prescribed surface.
    pub relative_velocity_before_mps: [f64; 3],
    pub relative_velocity_after_mps: [f64; 3],
    pub relative_displacement_m: [f64; 3],
    pub normal_gap_before_m: f64,
    pub normal_gap_after_m: f64,
    pub frictional_work_j: f64,
    /// Fresh original contact-law solve at this segment's actual duration.
    pub solve: SolveReport,
    pub(crate) duration: Dd,
}

/// Complete staged motion for the caller's geometry/feature admission.
/// The rigid-trace test proves the contact law's affine path, not completeness
/// of the contact set or validity of keeping a geometric feature/frame fixed.
pub struct SolidCoulombEventPath<'a, 'contact, B = ExactLinearSolid> {
    pub before: &'a B,
    pub contact: &'a CollinearCoulombContact<'contact>,
    pub phases: &'a [SolidCoulombPhase],
    pub displacement_after_m: &'a [f64],
    pub velocity_after_mps: &'a [f64],
    pub displacement_change_m: &'a [f64],
    pub external_force_n: &'a [f64],
    pub duration_s: f64,
}

/// One native clock advance, including a located stop and any static recapture
/// or reversed sliding. Relative loss includes both sides of a reversal even
/// when their net tangential displacement cancels.
#[derive(Debug)]
pub struct SolidCoulombEventStep<P, S = ExactLinearSolidStep> {
    pub native: S,
    pub contact_id: u64,
    pub phases: Vec<SolidCoulombPhase>,
    pub path_admission: P,
    /// Offset from the original interval's start, including an endpoint stop.
    pub stopping_time_offset_s: Option<f64>,
    /// Actual residual at the located stop; the native velocity is never clipped.
    pub stopping_speed_residual_mps: Option<f64>,
    /// Includes discarded static-capture trials as well as accepted phases.
    pub response_evaluations: usize,
    pub static_capture_trials: usize,
    pub displacement_change_m: Vec<f64>,
    pub momentum_change_kg_m_s: Vec<f64>,
    pub sliding_distance_m: f64,
    pub external_work_j: f64,
    pub prescribed_surface_work_j: f64,
    /// Independently integrated work supplied by native prescribed supports.
    /// Zero for owners without moving support coordinates.
    pub support_work_j: f64,
    pub normal_contact_work_j: f64,
    pub frictional_work_j: f64,
    pub frictional_dissipation_j: f64,
    pub energy_balance_error_j: f64,
}

impl ExactLinearSolid {
    /// Advance through a friction stop and static recapture or reversed sliding.
    ///
    /// A single initially closed contact must have relative normal speed within
    /// the declared velocity tolerance (actual residual motion is retained),
    /// collinear tangential motion, and independent rigid mass-response traces.
    /// Their products with the original stiffness and viscosity must be exactly
    /// zero. Other native elastic/viscous modes continue unchanged. Forces stay
    /// finite, and a changing load, detachment or rotating slip direction refuses.
    ///
    /// The original held external-force flow is evaluated once. Actual contact
    /// force impulses and first moments supply its rigid-coordinate correction.
    /// Every phase uses the shared force solver and original work/energy gates;
    /// all phases and static trials share the response and product budgets.
    pub fn advance_collinear_coulomb_events<P, E>(
        &mut self,
        contact: &CollinearCoulombContact<'_>,
        initial_force_n: [f64; 3],
        external_force_n: &[f64],
        controls: SolidSustainedControls,
        cx: &Cx<'_>,
        admit_path: impl FnOnce(&SolidCoulombEventPath<'_, '_>) -> Result<P, E>,
    ) -> Result<SolidCoulombEventStep<P>, SolidSustainedError<E>> {
        advance(
            self,
            contact,
            initial_force_n,
            external_force_n,
            controls,
            cx,
            admit_path,
        )
    }
}

pub(crate) struct TraceResponse {
    pub columns: [Vec<f64>; 3],
    pub modal_columns: Option<[Vec<f64>; 3]>,
    /// Additional original-row support reaction per contact force component.
    /// Free rows are zero; absent for unconstrained native coordinates.
    pub support_reaction_columns: Option<[Vec<f64>; 3]>,
    /// Acceleration of the trace due to prescribed support motion.
    pub acceleration_offset: [f64; 3],
}

pub(crate) struct ContactForcePath<'a, 'contact> {
    pub contact: &'a SolidSustainedContact<'contact>,
    pub phases: &'a [SolidCoulombPhase],
    pub external: &'a [f64],
    pub impulse: [Dd; 3],
    pub displacement_moment: [Dd; 3],
}

pub(crate) trait EventNative: Native {
    fn trace_response(
        &self,
        weights: [&[f64]; 3],
        cx: &Cx<'_>,
    ) -> Result<TraceResponse, Self::Error>;
    fn damping(&self) -> Option<&[f64]> {
        None
    }
    fn add_contact_moments(
        &self,
        scratch: &mut Self::Workspace,
        response: &TraceResponse,
        path: &ContactForcePath<'_, '_>,
        cx: &Cx<'_>,
    ) -> Result<f64, Self::Error>;
    fn accept_force_path(
        &mut self,
        scratch: Self::Workspace,
        supplied_j: f64,
        viscosity_j: f64,
        cx: &Cx<'_>,
    ) -> Result<Self::Step, Self::Error>;
}

pub(crate) fn advance<B: EventNative, P, E>(
    body: &mut B,
    law: &CollinearCoulombContact<'_>,
    initial_force: [f64; 3],
    external: &[f64],
    controls: SolidSustainedControls,
    cx: &Cx<'_>,
    admit_path: impl FnOnce(&SolidCoulombEventPath<'_, '_, B::Owner>) -> Result<P, E>,
) -> Result<SolidCoulombEventStep<P, B::Step>, SolidSustainedError<E, B::Error>> {
    poll(cx)?;
    body.admit(cx).map_err(native_error::<E, B>)?;
    let before = body.state();
    let n = before.displacement.len();
    let h = before.dt_s;
    let c = &law.contact;
    let axis_norm = law.tangent_axis[0].hypot(law.tangent_axis[1]);
    if external.len() != n
        || external.iter().any(|x| !x.is_finite())
        || initial_force.iter().any(|x| !x.is_finite())
        || c.gap_m != 0.0
        || c.velocity_weights
            .iter()
            .any(|row| row.len() != n || row.iter().any(|x| !x.is_finite()))
        || c.prescribed_surface_velocity_mps
            .iter()
            .any(|x| !x.is_finite())
        || !c.friction_coefficient.is_finite()
        || c.friction_coefficient < 0.0
        || !law.kinetic_coefficient.is_finite()
        || law.kinetic_coefficient < 0.0
        || law.kinetic_coefficient > c.friction_coefficient
        || !axis_norm.is_finite()
        || axis_norm <= 0.0
        || [
            controls.displacement_tolerance_m,
            controls.velocity_tolerance_mps,
            controls.force_tolerance_n,
            controls.complementarity_tolerance_j,
            controls.maximum_work_error_j,
        ]
        .iter()
        .any(|x| !x.is_finite() || *x < 0.0)
    {
        return Err(SolidSustainedError::Invalid(
            "finite complete rigid contact trace, closed gap, axis, friction pair and tolerances required",
        ));
    }
    let axis = law.tangent_axis.map(|x| x / axis_norm);
    let nn = n as u128;
    bound(
        "Coulomb event scratch bytes",
        body.workspace_bytes() + 8 * (32 * nn + 256),
        controls.maximum_workspace_bytes as u128,
    )?;
    let native_terms = body.response_multiply_terms() + 64 * nn * nn + 512 * nn;
    bound(
        "Coulomb event multiply terms",
        native_terms + 256,
        controls.maximum_response_multiply_terms as u128,
    )?;
    let response = body
        .trace_response(c.velocity_weights, cx)
        .map_err(native_error::<E, B>)?;
    audit::rigid_response::<E, B::Error>(
        &before,
        body.damping(),
        c.velocity_weights,
        &response,
        cx,
    )?;
    let mut mobility = [[0.0; 3]; 3];
    let mut acceleration = [0.0; 3];
    let mut velocity = [0.0; 3];
    for a in 0..3 {
        acceleration[a] = finite(
            dot(&response.columns[a], external) + Dd::from_f64(response.acceleration_offset[a]),
        )?;
        velocity[a] = finite(
            dot(c.velocity_weights[a], before.velocity)
                - Dd::from_f64(c.prescribed_surface_velocity_mps[a]),
        )?;
        for b in 0..3 {
            mobility[a][b] = finite(dot(c.velocity_weights[a], &response.columns[b]))?;
        }
    }
    if velocity[0] < -controls.velocity_tolerance_mps {
        return Err(SolidSustainedError::ImpactRequired {
            contact_id: c.id,
            normal_speed_mps: velocity[0],
        });
    }
    if velocity[0].abs() > controls.velocity_tolerance_mps {
        return Err(mode(
            c.id,
            "a persistently closed contact must start within the normal velocity allowance",
        ));
    }
    let matrix: Vec<_> = mobility.iter().flatten().copied().collect();
    fs_la::factor::cholesky(&matrix, 3)
        .map_err(|_| mode(c.id, "independent positive contact mobility required"))?;
    phase::collinear::<E, B::Error>(c.id, axis, velocity, controls.velocity_tolerance_mps)?;
    let mut scratch = body.workspace().map_err(native_error::<E, B>)?;
    body.evaluate(external, &mut scratch, cx)
        .map_err(native_error::<E, B>)?;
    audit::free_trace::<E, B::Error>(
        &before,
        body.displacement(&scratch),
        body.velocity(&scratch),
        c,
        velocity,
        acceleration,
        controls,
    )?;
    let baseline_velocity = body.velocity(&scratch).to_vec();
    let native_work = body
        .work(external, &scratch, cx)
        .map_err(native_error::<E, B>)?;
    let mut engine = phase::Engine {
        mobility,
        acceleration,
        axis,
        contact: c,
        kinetic: law.kinetic_coefficient,
        duration: h,
        controls,
        native_terms,
        evaluations: 0,
        static_trials: 0,
        cx,
    };
    let path = engine.resolve::<B, E>(velocity, initial_force)?;
    let mut impulse = [Dd::ZERO; 3];
    let mut moment = [Dd::ZERO; 3];
    for phase in &path.phases {
        let dt = phase.duration;
        let first =
            (Dd::from_f64(h) - Dd::from_f64(phase.start_offset_s) - Dd::from_f64(0.5) * dt) * dt;
        for a in 0..3 {
            impulse[a] = impulse[a] + Dd::from_f64(phase.force_n[a]) * dt;
            moment[a] = moment[a] + Dd::from_f64(phase.force_n[a]) * first;
        }
    }
    let support_work = body
        .add_contact_moments(
            &mut scratch,
            &response,
            &ContactForcePath {
                contact: c,
                phases: &path.phases,
                external,
                impulse,
                displacement_moment: moment,
            },
            cx,
        )
        .map_err(native_error::<E, B>)?;
    let du: Vec<_> = body
        .displacement(&scratch)
        .iter()
        .zip(before.displacement)
        .map(|(after, before)| after - before)
        .collect();
    let measured = audit::final_path(
        &before,
        body.displacement(&scratch),
        body.velocity(&scratch),
        &baseline_velocity,
        &du,
        c,
        &path.phases,
        &response,
        impulse,
        external,
        native_work.viscosity_j,
        support_work,
        controls,
        cx,
    )?;
    let path_admission = admit_path(&SolidCoulombEventPath {
        before: body.owner(),
        contact: law,
        phases: &path.phases,
        displacement_after_m: body.displacement(&scratch),
        velocity_after_mps: body.velocity(&scratch),
        displacement_change_m: &du,
        external_force_n: external,
        duration_s: h,
    })
    .map_err(SolidSustainedError::Path)?;
    poll(cx)?;
    let native = body
        .accept_force_path(scratch, measured.supplied, native_work.viscosity_j, cx)
        .map_err(native_error::<E, B>)?;
    Ok(SolidCoulombEventStep {
        native,
        contact_id: c.id,
        phases: path.phases,
        path_admission,
        stopping_time_offset_s: path.stop_time,
        stopping_speed_residual_mps: path.stop_residual,
        response_evaluations: engine.evaluations,
        static_capture_trials: engine.static_trials,
        displacement_change_m: du,
        momentum_change_kg_m_s: measured.momentum,
        sliding_distance_m: measured.distance,
        external_work_j: measured.external,
        prescribed_surface_work_j: measured.surface,
        support_work_j: support_work,
        normal_contact_work_j: measured.normal,
        frictional_work_j: measured.friction,
        frictional_dissipation_j: (-measured.friction).max(0.0),
        energy_balance_error_j: measured.energy_error,
    })
}
