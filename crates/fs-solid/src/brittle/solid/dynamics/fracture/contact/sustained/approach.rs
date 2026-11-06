//! Released-face arrival while existing contacts retain native sustained forces.
//!
//! Only initially resting, active and matched rows enter the force solve.
//! Other released pairs receive no contact load during a prefix. Every complete
//! constrained response is then checked against ALL released openings. Unsafe
//! unforced endpoints bound a native prefix search; no failed force/material
//! solve is reinterpreted as an impact or repaired by a free-flight substitute.

use fs_contact::sustained::solve::{MAX_CONTACTS, SolveError};
use fs_exec::Cx;

use super::{
    BrittleFaceSustainedControls, BrittleFaceSustainedError, BrittleFaceSustainedStep,
    SustainedScope, pair, resting, validate_controls,
};
use crate::brittle::solid::dynamics::BrittleTetMotion;
use crate::brittle::solid::dynamics::fracture::contact::BrittleFaceImpactRow;
use crate::brittle::solid::dynamics::fracture::contact::approach::{
    BrittleFaceApproachControls, BrittleFaceApproachError, BrittleFaceContactLocation, PrefixLoads,
    ReleasedFaceMotionView, scan_released_face_motion, validate,
};
use crate::brittle::solid::fracture::require_released_faces;
use crate::brittle::solid::{BrittleTetError, checkpoint};
use crate::plastic3::IncrementSettings;

/// One physically accepted constrained prefix, or an already approaching pair.
#[derive(Debug)]
pub struct BrittleFaceSustainedApproachResult<Step> {
    pub advance: Option<Step>,
    pub contact: Option<BrittleFaceContactLocation>,
    /// Completed native endpoint solves, including unsafe discarded endpoints.
    pub trials: usize,
    /// Every native mechanical response spent by those force solves, including
    /// Newton probes, feasibility corrections and final original replays.
    pub response_evaluations: usize,
    pub numerical_tangent_evaluations: usize,
}

pub type BrittleFaceSustainedApproachStep<'a> =
    BrittleFaceSustainedApproachResult<BrittleFaceSustainedStep<'a>>;

/// A complete original force response, including private source candidates when
/// present. The locator retains only a safe complete response, never its fields.
pub(crate) trait NativeSustainedApproachTrial {
    fn face_motion(&self) -> ReleasedFaceMotionView<'_>;
    fn rows(&self) -> &[BrittleFaceImpactRow];
    fn physical_forces_n(&self) -> &[[f64; 3]];
    fn response_evaluations(&self) -> usize;
    fn tangent_evaluations(&self) -> usize;
}

impl NativeSustainedApproachTrial for BrittleFaceSustainedStep<'_> {
    fn face_motion(&self) -> ReleasedFaceMotionView<'_> {
        self.motion.released_face_motion_view()
    }
    fn rows(&self) -> &[BrittleFaceImpactRow] {
        &self.rows
    }
    fn physical_forces_n(&self) -> &[[f64; 3]] {
        &self.solve.forces_n
    }
    fn response_evaluations(&self) -> usize {
        self.solve.response_evaluations
    }
    fn tangent_evaluations(&self) -> usize {
        self.numerical_tangent_evaluations
    }
}

#[derive(Debug)]
pub enum BrittleFaceSustainedApproachError {
    Invalid(&'static str),
    Approach(BrittleFaceApproachError),
    Sustained(BrittleFaceSustainedError),
    ResponseBudget { evaluations: usize, maximum: usize },
}
impl core::fmt::Display for BrittleFaceSustainedApproachError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "constrained released-face approach refused: {self:?}")
    }
}
impl std::error::Error for BrittleFaceSustainedApproachError {}
impl From<BrittleFaceApproachError> for BrittleFaceSustainedApproachError {
    fn from(error: BrittleFaceApproachError) -> Self {
        Self::Approach(error)
    }
}
impl From<BrittleFaceSustainedError> for BrittleFaceSustainedApproachError {
    fn from(error: BrittleFaceSustainedError) -> Self {
        Self::Sustained(error)
    }
}
impl From<BrittleTetError> for BrittleFaceSustainedApproachError {
    fn from(error: BrittleTetError) -> Self {
        Self::Approach(error.into())
    }
}

impl<'a> BrittleTetMotion<'a> {
    /// Locate arrival without releasing the constraints already at rest.
    ///
    /// Applied forces are affine endpoint pairs over the requested interval;
    /// every actual native prefix receives its own exact force mean and support
    /// endpoint. `sustained.work` gates constrained responses. Both work-control
    /// sets and matching activation/speed allowances are admitted before entry.
    /// The sustained response-evaluation budget applies to the ENTIRE search;
    /// `approach.maximum_trials` separately caps endpoint solves.
    ///
    /// An initially approaching active pair returns an immediate contact. In
    /// its absence, at least one matched active pair at normal rest is required.
    /// That fixed initial subset alone receives native sustained normal forces;
    /// it may unload and separate. Initially open pairs remain unforced, even
    /// when an unpublished trial crosses them. All accepted endpoints satisfy
    /// strict nonnegative gaps, loaded-row mode/work checks and the native guard.
    ///
    /// The accepted path is the native affine displacement trace. The returned
    /// safe/unsafe endpoint bracket has the original locator's numerical meaning,
    /// not a globally earliest continuous-collision claim. Source, constitutive,
    /// force-solve, mode and work refusals are terminal. Retrying uses this same
    /// immutable predecessor; only the returned prefix consumes a physical step.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_toward_released_face_contact_with_sustained(
        &self,
        duration_s: f64,
        prescribed_end_m: &[(usize, f64)],
        force_intervals_n: &[(usize, [f64; 2])],
        settings: IncrementSettings,
        approach: BrittleFaceApproachControls,
        sustained: BrittleFaceSustainedControls,
        cx: &Cx<'_>,
    ) -> Result<BrittleFaceSustainedApproachStep<'a>, BrittleFaceSustainedApproachError> {
        self.advance_toward_released_face_contact_with_sustained_from(
            duration_s,
            prescribed_end_m,
            force_intervals_n,
            None,
            settings,
            approach,
            sustained,
            cx,
        )
    }

    /// Seeds identify actual paired rows; changing subsets never shifts forces
    /// onto unrelated rows. They affect only the numerical initial iterate.
    #[allow(clippy::too_many_arguments)]
    pub(in crate::brittle::solid::dynamics::fracture::contact) fn advance_toward_released_face_contact_with_sustained_from(
        &self,
        duration_s: f64,
        prescribed_end_m: &[(usize, f64)],
        force_intervals_n: &[(usize, [f64; 2])],
        seed: Option<(&[BrittleFaceImpactRow], &[[f64; 3]])>,
        settings: IncrementSettings,
        approach: BrittleFaceApproachControls,
        sustained: BrittleFaceSustainedControls,
        cx: &Cx<'_>,
    ) -> Result<BrittleFaceSustainedApproachStep<'a>, BrittleFaceSustainedApproachError> {
        use BrittleFaceSustainedApproachError as A;
        if seed.is_some_and(|(_, forces)| forces.iter().any(|f| f[1] != 0.0 || f[2] != 0.0)) {
            return Err(A::Invalid(
                "normal force seeds must have zero tangent components",
            ));
        }
        approach_released_faces_with_forces(
            self.released_face_motion_view(),
            duration_s,
            prescribed_end_m,
            force_intervals_n,
            seed,
            settings,
            approach,
            sustained,
            cx,
            || {
                let mut rows = Vec::new();
                for face in self.solid.released_faces.iter() {
                    checkpoint(cx)?;
                    for vertex in 0..3 {
                        let nodes = [face.paired_nodes[0][vertex], face.paired_nodes[1][vertex]];
                        if nodes[0] == nodes[1] {
                            continue;
                        }
                        let row = BrittleFaceImpactRow {
                            face_id: face.id,
                            vertex,
                            nodes,
                        };
                        if !resting(pair(self, row, face.normal)?, sustained) {
                            continue;
                        }
                        if rows.len() == sustained.solve.maximum_contacts {
                            return Err(A::Invalid(
                                "all initially resting rows must fit the native force budget",
                            ));
                        }
                        rows.push(row);
                    }
                }
                Ok(rows)
            },
            |dt, prescribed, forces, _rows, physical_seeds, controls| {
                let normal_seeds: Vec<_> = physical_seeds.iter().map(|f| f[0]).collect();
                self.sustained_candidate(
                    dt,
                    prescribed,
                    forces,
                    Some(&normal_seeds),
                    SustainedScope::InitialResting,
                    settings,
                    controls,
                    cx,
                )
                .map_err(A::from)
            },
            |error| {
                if let A::Sustained(BrittleFaceSustainedError::Solve(source)) = error
                    && let SolveError::ResponseBudget { evaluations, .. } = source.as_ref()
                {
                    Some(*evaluations)
                } else {
                    None
                }
            },
        )
    }
}

/// One bounded constrained-prefix locator for original normal, directional and
/// sourced force responses. Each candidate and its private samples remain tied
/// together; only a safe actual endpoint can be returned to a physical owner.
#[allow(clippy::too_many_arguments, clippy::too_many_lines)]
pub(crate) fn approach_released_faces_with_forces<
    Step: NativeSustainedApproachTrial,
    E: From<BrittleFaceSustainedApproachError>,
>(
    initial: ReleasedFaceMotionView<'_>,
    duration_s: f64,
    prescribed_end_m: &[(usize, f64)],
    force_intervals_n: &[(usize, [f64; 2])],
    seed: Option<(&[BrittleFaceImpactRow], &[[f64; 3]])>,
    settings: IncrementSettings,
    approach: BrittleFaceApproachControls,
    sustained: BrittleFaceSustainedControls,
    cx: &Cx<'_>,
    select_rows: impl FnOnce() -> Result<Vec<BrittleFaceImpactRow>, E>,
    mut candidate: impl FnMut(
        f64,
        &[(usize, f64)],
        &[(usize, f64)],
        &[BrittleFaceImpactRow],
        &[[f64; 3]],
        BrittleFaceSustainedControls,
    ) -> Result<Step, E>,
    spent_budget: impl Fn(&E) -> Option<usize>,
) -> Result<BrittleFaceSustainedApproachResult<Step>, E> {
    use BrittleFaceSustainedApproachError as A;
    checkpoint(cx).map_err(A::from)?;
    let end_time_s = initial.time_s + duration_s;
    validate(duration_s, initial.time_s, end_time_s, settings, approach).map_err(A::from)?;
    validate(
        duration_s,
        initial.time_s,
        end_time_s,
        settings,
        BrittleFaceApproachControls {
            work: sustained.work,
            ..approach
        },
    )
    .map_err(A::from)?;
    validate_controls(sustained).map_err(A::from)?;
    if approach.activation_tolerance_m != sustained.activation_tolerance_m
        || approach.velocity_tolerance_mps != sustained.velocity_tolerance_mps
        || sustained.solve.maximum_response_evaluations == 0
        || !(1..=MAX_CONTACTS).contains(&sustained.solve.maximum_contacts)
        || !(1..=256).contains(&sustained.solve.maximum_iterations)
        || sustained.solve.newton.linear_restart > 64
        || sustained.solve.newton.max_linear_cycles > 128
        || !sustained.solve.difference_step_n.is_finite()
        || sustained.solve.difference_step_n <= 0.0
        || !sustained.solve.force_residual_tolerance_n.is_finite()
        || sustained.solve.force_residual_tolerance_n <= 0.0
    {
        return Err(A::Invalid(
            "compatible localization gates and bounded native force controls required",
        )
        .into());
    }
    sustained.solve.newton.validate().map_err(|error| {
        A::Sustained(BrittleFaceSustainedError::Solve(Box::new(
            SolveError::Newton(error),
        )))
    })?;
    if let Some((rows, forces)) = seed {
        if rows.len() != forces.len()
            || rows.len() > MAX_CONTACTS
            || forces
                .iter()
                .any(|f| f.iter().any(|v| !v.is_finite()) || f[0] < 0.0)
            || rows
                .iter()
                .enumerate()
                .any(|(i, row)| rows[..i].contains(row))
        {
            return Err(A::Invalid(
                    "distinct paired physical force seeds with finite nonnegative normal components required",
                ).into());
        }
    }
    let mut loads = PrefixLoads::new(
        initial.displacement_m.as_flattened(),
        prescribed_end_m,
        force_intervals_n,
        cx,
    )
    .map_err(A::from)?;
    let scan = |motion: ReleasedFaceMotionView<'_>| {
        scan_released_face_motion(
            motion,
            approach.activation_tolerance_m,
            approach.velocity_tolerance_mps,
            cx,
        )
    };
    if let Some(contact) = scan(initial).map_err(A::from)? {
        return Ok(BrittleFaceSustainedApproachResult {
            advance: None,
            contact: Some(contact),
            trials: 0,
            response_evaluations: 0,
            numerical_tangent_evaluations: 0,
        });
    }
    let rows = select_rows()?;
    let mut physical_seeds: Vec<_> = rows
        .iter()
        .map(|row| {
            seed.and_then(|(old_rows, forces)| {
                old_rows
                    .iter()
                    .position(|old| old == row)
                    .map(|i| forces[i])
            })
            .unwrap_or([0.0; 3])
        })
        .collect();
    if rows.is_empty() {
        return Err(
            A::Invalid("at least one initially matched active resting pair is required").into(),
        );
    }
    if rows.len() > sustained.solve.maximum_contacts || rows.len() > MAX_CONTACTS {
        return Err(A::Invalid("initial constrained rows exceed the force budget").into());
    }

    let maximum = sustained.solve.maximum_response_evaluations;
    let mut response_evaluations = 0usize;
    let mut trials = 0usize;
    let mut tangent_evaluations = 0usize;
    let mut lower_time_s = initial.time_s;
    let mut upper_time_s = end_time_s;
    let mut target = end_time_s;
    let mut bracketed = false;
    let mut safe = None;
    let mut witness = None;
    loop {
        checkpoint(cx).map_err(A::from)?;
        if trials == approach.maximum_trials {
            return Err(A::from(BrittleFaceApproachError::WorkLimit {
                trials,
                lower_time_s,
                upper_time_s,
            })
            .into());
        }
        if response_evaluations == maximum {
            return Err(A::ResponseBudget {
                evaluations: response_evaluations,
                maximum,
            }
            .into());
        }
        let fraction = if target == end_time_s {
            1.0
        } else {
            (target - initial.time_s) / (end_time_s - initial.time_s)
        };
        loads.at(fraction, cx).map_err(A::from)?;
        let mut controls = sustained;
        controls.solve.maximum_response_evaluations = maximum - response_evaluations;
        trials += 1;
        let step = candidate(
            target - initial.time_s,
            &loads.prescribed,
            &loads.forces,
            &rows,
            &physical_seeds,
            controls,
        )
        .map_err(|error| {
            if let Some(evaluations) = spent_budget(&error) {
                return A::ResponseBudget {
                    evaluations: response_evaluations + evaluations,
                    maximum,
                }
                .into();
            }
            error
        })?;
        response_evaluations = response_evaluations
            .checked_add(step.response_evaluations())
            .ok_or(A::Invalid("native response count overflow"))?;
        tangent_evaluations = tangent_evaluations
            .checked_add(step.tangent_evaluations())
            .ok_or(A::Invalid("completed tangent count overflow"))?;
        if step.face_motion().time_s != target
            || step.rows() != rows
            || step.physical_forces_n().len() != rows.len()
            || response_evaluations > maximum
        {
            return Err(A::Invalid(
                "native constrained prefix escaped its clock, rows or response allowance",
            )
            .into());
        }
        for (value, force) in physical_seeds.iter_mut().zip(step.physical_forces_n()) {
            *value = *force;
        }
        match scan(step.face_motion()) {
            Ok(contact) => {
                require_released_faces(
                    initial.displacement_m.len(),
                    initial.released_faces,
                    step.face_motion().displacement_m,
                    cx,
                )
                .map_err(A::from)?;
                if !bracketed {
                    return Ok(BrittleFaceSustainedApproachResult {
                        advance: Some(step),
                        contact,
                        trials,
                        response_evaluations,
                        numerical_tangent_evaluations: tangent_evaluations,
                    });
                }
                lower_time_s = target;
                safe = Some(step);
                witness = contact;
            }
            Err(error @ BrittleFaceApproachError::Overlap { face, vertex, .. }) => {
                if rows
                    .iter()
                    .any(|row| row.face_id == face && row.vertex == vertex)
                {
                    // A participating-force failure is never a new arrival.
                    return Err(A::from(error).into());
                }
                bracketed = true;
                upper_time_s = target;
            }
            Err(error) => return Err(A::from(error).into()),
        }
        if upper_time_s - lower_time_s <= approach.time_tolerance_s
            && let Some(mut contact) = witness
        {
            contact.bracket_end_time_s = upper_time_s;
            checkpoint(cx).map_err(A::from)?;
            return Ok(BrittleFaceSustainedApproachResult {
                advance: safe,
                contact: Some(contact),
                trials,
                response_evaluations,
                numerical_tangent_evaluations: tangent_evaluations,
            });
        }
        target = lower_time_s + 0.5 * (upper_time_s - lower_time_s);
        if target <= lower_time_s || target >= upper_time_s {
            return Err(A::from(BrittleFaceApproachError::Unresolved {
                lower_time_s,
                upper_time_s,
                reason: "no representable constrained approaching prefix within both allowances",
            })
            .into());
        }
    }
}
