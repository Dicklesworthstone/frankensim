//! Bounded native midpoint prefixes approaching a released-face normal event.
//!
//! Every trial starts from the same accepted motion and uses the original
//! material, inertia, Newton and opening guards. A successful prefix has
//! nonnegative paired openings along its affine displacement path. Bisection
//! bounds sampled endpoint admissibility; it does not certify the earliest
//! collision of a continuous or nonmonotone trajectory. Neither penetration
//! projection nor interpolation of material history supplies a candidate.

use fs_exec::Cx;

use crate::brittle::solid::dynamics::midpoint::{
    BrittleTetMidpointControls, BrittleTetMidpointStep,
};
use crate::brittle::solid::dynamics::{BrittleTetDynamicError, BrittleTetMotion};
use crate::brittle::solid::fracture::ReleasedBrittleFace;
use crate::brittle::solid::path::interpolate;
use crate::brittle::solid::{BrittleTetError, checkpoint};
use crate::plastic3::IncrementSettings;

/// Original midpoint work gates and explicit localization allowances.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BrittleFaceApproachControls {
    pub work: BrittleTetMidpointControls,
    /// Strictly positive distance at which an approaching pair may activate.
    pub activation_tolerance_m: f64,
    /// A contact witness must have normal speed below the negative allowance.
    pub velocity_tolerance_mps: f64,
    /// Strictly positive maximum safe/unsafe native endpoint time separation.
    pub time_tolerance_s: f64,
    /// Counts every native attempt, including rejected overlap endpoints.
    pub maximum_trials: usize,
}

/// One approaching witness; the later impact still gathers ALL active pairs.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BrittleFaceContactLocation {
    pub face_id: u64,
    pub vertex: usize,
    pub gap_m: f64,
    pub relative_normal_speed_mps: f64,
    pub time_s: f64,
    /// A sampled unsafe native endpoint, when bisection was necessary.
    /// Equals `time_s` for initial proximity or a successful full endpoint;
    /// that equality is not an exact physical collision-time certificate.
    pub bracket_end_time_s: f64,
}

/// Complete native prefix, or an approaching contact at the original clock.
#[derive(Debug)]
pub struct BrittleFaceApproachStep<'a> {
    /// `None` only for an initially active approaching pair; no step is spent.
    pub advance: Option<BrittleTetMidpointStep<'a>>,
    /// `None` means the complete requested native interval was accepted.
    pub contact: Option<BrittleFaceContactLocation>,
    pub trials: usize,
    /// All successfully completed native trials, including discarded prefixes.
    /// Native refusals do not expose their internal fallback counts.
    pub numerical_tangent_evaluations: usize,
}

#[derive(Clone, Debug, PartialEq)]
pub enum BrittleFaceApproachError {
    Invalid(&'static str),
    Endpoint(BrittleTetDynamicError),
    Overlap {
        face: u64,
        vertex: usize,
        gap_m: f64,
    },
    WorkLimit {
        trials: usize,
        lower_time_s: f64,
        upper_time_s: f64,
    },
    Unresolved {
        lower_time_s: f64,
        upper_time_s: f64,
        reason: &'static str,
    },
}
impl core::fmt::Display for BrittleFaceApproachError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "released brittle face approach refused: {self:?}")
    }
}
impl std::error::Error for BrittleFaceApproachError {}
impl From<BrittleTetError> for BrittleFaceApproachError {
    fn from(error: BrittleTetError) -> Self {
        Self::Endpoint(error.into())
    }
}

impl<'a> BrittleTetMotion<'a> {
    /// Advance a native midpoint prefix toward approaching released faces.
    ///
    /// Forces are affine endpoint pairs over the requested interval; each
    /// trial uses their exact prefix mean. Prescribed displacements follow the
    /// affine path from this accepted displacement to the supplied endpoint.
    /// Every shortened solve replays the original accepted material history.
    /// All failed or discarded trials leave this motion unchanged. Only the
    /// returned prefix consumes one physical step; initial contact consumes none.
    ///
    /// A nonnegative full endpoint is accepted directly. Otherwise bounded
    /// bisection retains a safe native prefix until an approaching witness is
    /// within the distance allowance AND the sampled endpoint bracket meets
    /// the time allowance. Arbitrarily small negative gaps always refuse the
    /// candidate, including gaps within the solid's opening roundoff allowance.
    /// Cancellation, material/solver/work refusals and exhausted budgets are
    /// terminal; only explicit endpoint overlap narrows the contact bracket.
    ///
    /// These are discrete affine-path and sampled endpoint claims, not global
    /// continuous collision detection. Intermediate closing and reopening of
    /// a different trajectory is not detected. A closed pair at rest driven
    /// inward may have no admissible prefix and refuses; this method does not
    /// supply persistent contact forces or reset velocity to manufacture one.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_toward_released_face_contact(
        &self,
        duration_s: f64,
        prescribed_end_m: &[(usize, f64)],
        force_intervals_n: &[(usize, [f64; 2])],
        settings: IncrementSettings,
        controls: BrittleFaceApproachControls,
        cx: &Cx<'_>,
    ) -> Result<BrittleFaceApproachStep<'a>, BrittleFaceApproachError> {
        let result = approach_released_faces(
            self.released_face_motion_view(),
            duration_s,
            prescribed_end_m,
            force_intervals_n,
            settings,
            controls,
            cx,
            |dt_s, prescribed, forces| {
                self.advance_midpoint(dt_s, prescribed, forces, settings, controls.work, cx)
            },
        )?;
        Ok(BrittleFaceApproachStep {
            advance: result.advance,
            contact: result.contact,
            trials: result.trials,
            numerical_tangent_evaluations: result.numerical_tangent_evaluations,
        })
    }

    pub(crate) fn released_face_motion_view(&self) -> ReleasedFaceMotionView<'_> {
        ReleasedFaceMotionView {
            displacement_m: &self.solid.solution.displacement_m,
            velocity_mps: &self.velocity_mps,
            released_faces: &self.solid.released_faces,
            time_s: self.time_s,
        }
    }
}

/// Physical fields needed by event localization, without a material conversion.
#[derive(Clone, Copy, Debug)]
pub(crate) struct ReleasedFaceMotionView<'a> {
    pub(crate) displacement_m: &'a [[f64; 3]],
    pub(crate) velocity_mps: &'a [[f64; 3]],
    pub(crate) released_faces: &'a [ReleasedBrittleFace],
    pub(crate) time_s: f64,
}

/// A complete native prefix keeps its original constitutive and work receipt.
pub(crate) struct NativeFaceApproach<Step> {
    pub(crate) advance: Option<Step>,
    pub(crate) contact: Option<BrittleFaceContactLocation>,
    pub(crate) trials: usize,
    pub(crate) numerical_tangent_evaluations: usize,
}

/// Read-only recovery from the actual native mechanical or thermal endpoint.
pub(crate) trait NativeFaceApproachTrial {
    fn face_motion(&self) -> ReleasedFaceMotionView<'_>;
    fn tangent_evaluations(&self) -> usize;
}

impl NativeFaceApproachTrial for BrittleTetMidpointStep<'_> {
    fn face_motion(&self) -> ReleasedFaceMotionView<'_> {
        self.motion.released_face_motion_view()
    }

    fn tangent_evaluations(&self) -> usize {
        self.numerical_tangent_evaluations
    }
}

/// One bounded localization algorithm for both original material owners.
/// Each callback replays its accepted native history; only a returned prefix
/// can consume physical time. No interpolated constitutive state enters here.
#[allow(clippy::too_many_arguments, clippy::too_many_lines)]
pub(crate) fn approach_released_faces<Step: NativeFaceApproachTrial>(
    initial: ReleasedFaceMotionView<'_>,
    duration_s: f64,
    prescribed_end_m: &[(usize, f64)],
    force_intervals_n: &[(usize, [f64; 2])],
    settings: IncrementSettings,
    controls: BrittleFaceApproachControls,
    cx: &Cx<'_>,
    mut advance: impl FnMut(
        f64,
        &[(usize, f64)],
        &[(usize, f64)],
    ) -> Result<Step, BrittleTetDynamicError>,
) -> Result<NativeFaceApproach<Step>, BrittleFaceApproachError> {
    use BrittleFaceApproachError as E;
    checkpoint(cx)?;
    let end_time_s = initial.time_s + duration_s;
    validate(duration_s, initial.time_s, end_time_s, settings, controls)?;
    let mut loads = PrefixLoads::new(
        initial.displacement_m.as_flattened(),
        prescribed_end_m,
        force_intervals_n,
        cx,
    )?;
    let scan = |motion: ReleasedFaceMotionView<'_>| {
        scan_released_face_motion(
            motion,
            controls.activation_tolerance_m,
            controls.velocity_tolerance_mps,
            cx,
        )
    };
    if let Some(contact) = scan(initial)? {
        return Ok(NativeFaceApproach {
            advance: None,
            contact: Some(contact),
            trials: 0,
            numerical_tangent_evaluations: 0,
        });
    }

    let mut lower_time_s = initial.time_s;
    let mut upper_time_s = end_time_s;
    let mut target = end_time_s;
    let mut bracketed = false;
    let mut safe = None;
    let mut witness = None;
    let mut trials = 0;
    let mut tangent_evaluations = 0_usize;
    loop {
        checkpoint(cx)?;
        if trials == controls.maximum_trials {
            return Err(E::WorkLimit {
                trials,
                lower_time_s,
                upper_time_s,
            });
        }
        let fraction = if target == end_time_s {
            1.0
        } else {
            (target - initial.time_s) / (end_time_s - initial.time_s)
        };
        loads.at(fraction, cx)?;
        trials += 1;
        let result = advance(target - initial.time_s, &loads.prescribed, &loads.forces);
        match result {
            Ok(step) => {
                tangent_evaluations = tangent_evaluations
                    .checked_add(step.tangent_evaluations())
                    .ok_or(E::Invalid("completed tangent count overflow"))?;
                if step.face_motion().time_s != target {
                    return Err(E::Invalid("native clock did not reach requested prefix"));
                }
                match scan(step.face_motion()) {
                    Ok(contact) => {
                        if !bracketed {
                            return Ok(NativeFaceApproach {
                                advance: Some(step),
                                contact,
                                trials,
                                numerical_tangent_evaluations: tangent_evaluations,
                            });
                        }
                        lower_time_s = target;
                        safe = Some(step);
                        witness = contact;
                    }
                    Err(E::Overlap { .. }) => {
                        bracketed = true;
                        upper_time_s = target;
                    }
                    Err(error) => return Err(error),
                }
            }
            Err(BrittleTetDynamicError::Solid(BrittleTetError::ReleasedFaceClosure { .. })) => {
                bracketed = true;
                upper_time_s = target;
            }
            Err(error) => return Err(E::Endpoint(error)),
        }
        if upper_time_s - lower_time_s <= controls.time_tolerance_s
            && let Some(mut contact) = witness
        {
            contact.bracket_end_time_s = upper_time_s;
            return Ok(NativeFaceApproach {
                advance: safe,
                contact: Some(contact),
                trials,
                numerical_tangent_evaluations: tangent_evaluations,
            });
        }
        target = lower_time_s + 0.5 * (upper_time_s - lower_time_s);
        if target <= lower_time_s || target >= upper_time_s {
            return Err(E::Unresolved {
                lower_time_s,
                upper_time_s,
                reason: "no representable safe approaching prefix within both allowances",
            });
        }
    }
}

/// Original native-motion entry, retained for contact-time mode composition.
pub(super) fn scan_released_faces(
    motion: &BrittleTetMotion<'_>,
    activation_tolerance_m: f64,
    velocity_tolerance_mps: f64,
    cx: &Cx<'_>,
) -> Result<Option<BrittleFaceContactLocation>, BrittleFaceApproachError> {
    scan_released_face_motion(
        motion.released_face_motion_view(),
        activation_tolerance_m,
        velocity_tolerance_mps,
        cx,
    )
}

/// Scan every non-tied pair before returning the first ordered witness.
/// A later negative pair cannot be hidden by an earlier active approach.
pub(crate) fn scan_released_face_motion(
    motion: ReleasedFaceMotionView<'_>,
    activation_tolerance_m: f64,
    velocity_tolerance_mps: f64,
    cx: &Cx<'_>,
) -> Result<Option<BrittleFaceContactLocation>, BrittleFaceApproachError> {
    let mut witness = None;
    for face in motion.released_faces {
        checkpoint(cx)?;
        for vertex in 0..3 {
            let [a, b] = [face.paired_nodes[0][vertex], face.paired_nodes[1][vertex]];
            if a == b {
                continue;
            }
            // Equal reference coordinates cancel before arithmetic, just as
            // in the original released-face impact owner.
            let separation = core::array::from_fn(|axis| {
                motion.displacement_m[b][axis] - motion.displacement_m[a][axis]
            });
            let relative_velocity = core::array::from_fn(|axis| {
                motion.velocity_mps[b][axis] - motion.velocity_mps[a][axis]
            });
            let gap_m = super::dot(face.normal, separation);
            let relative_normal_speed_mps = super::dot(face.normal, relative_velocity);
            if !gap_m.is_finite() || !relative_normal_speed_mps.is_finite() {
                return Err(BrittleFaceApproachError::Invalid(
                    "finite paired opening and normal speed required",
                ));
            }
            if gap_m < 0.0 {
                return Err(BrittleFaceApproachError::Overlap {
                    face: face.id,
                    vertex,
                    gap_m,
                });
            }
            if witness.is_none()
                && gap_m <= activation_tolerance_m
                && relative_normal_speed_mps < -velocity_tolerance_mps
            {
                witness = Some(BrittleFaceContactLocation {
                    face_id: face.id,
                    vertex,
                    gap_m,
                    relative_normal_speed_mps,
                    time_s: motion.time_s,
                    bracket_end_time_s: motion.time_s,
                });
            }
        }
    }
    Ok(witness)
}

pub(crate) fn validate(
    duration_s: f64,
    start: f64,
    end: f64,
    settings: IncrementSettings,
    controls: BrittleFaceApproachControls,
) -> Result<(), BrittleFaceApproachError> {
    if !duration_s.is_finite()
        || duration_s <= 0.0
        || !end.is_finite()
        || end <= start
        || controls.maximum_trials == 0
        || !controls.activation_tolerance_m.is_finite()
        || controls.activation_tolerance_m <= 0.0
        || !controls.time_tolerance_s.is_finite()
        || controls.time_tolerance_s <= 0.0
        || [
            controls.velocity_tolerance_mps,
            controls.work.material_energy_tolerance_j,
            controls.work.equilibrium_work_tolerance_j,
            controls.work.kinetic_work_tolerance_j,
            controls.work.total_energy_tolerance_j,
        ]
        .iter()
        .any(|value| !value.is_finite() || *value < 0.0)
        || controls
            .work
            .difference_step_strain
            .is_some_and(|h| !h.is_finite() || h <= 0.0)
        || settings.max_newton_iterations == 0
        || !settings.relative_tolerance.is_finite()
        || settings.relative_tolerance <= 0.0
        || !settings.linear_tolerance.is_finite()
        || settings.linear_tolerance <= 0.0
    {
        return Err(BrittleFaceApproachError::Invalid(
            "representable positive duration, localization allowances and native work/solver controls required",
        ));
    }
    Ok(())
}

struct PrefixLoads<'a> {
    initial_u: &'a [f64],
    prescribed_end: &'a [(usize, f64)],
    force_intervals: &'a [(usize, [f64; 2])],
    prescribed: Vec<(usize, f64)>,
    forces: Vec<(usize, f64)>,
}
impl<'a> PrefixLoads<'a> {
    fn new(
        initial_u: &'a [f64],
        prescribed_end: &'a [(usize, f64)],
        force_intervals: &'a [(usize, [f64; 2])],
        cx: &Cx<'_>,
    ) -> Result<Self, BrittleFaceApproachError> {
        use BrittleFaceApproachError as E;
        if prescribed_end.len() > initial_u.len() || force_intervals.len() > initial_u.len() {
            return Err(E::Invalid("load entries exceed original DOF count"));
        }
        let mut prescribed = Vec::new();
        let mut forces = Vec::new();
        prescribed
            .try_reserve_exact(prescribed_end.len())
            .and_then(|()| forces.try_reserve_exact(force_intervals.len()))
            .map_err(|_| E::Invalid("prefix load allocation refused"))?;
        for &(dof, end) in prescribed_end {
            checkpoint(cx)?;
            if dof >= initial_u.len() || !end.is_finite() {
                return Err(E::Invalid(
                    "finite prescribed endpoint on existing DOF required",
                ));
            }
            prescribed.push((dof, initial_u[dof]));
        }
        for &(dof, [start, end]) in force_intervals {
            checkpoint(cx)?;
            if dof >= initial_u.len() || !start.is_finite() || !end.is_finite() {
                return Err(E::Invalid(
                    "finite force endpoints on existing DOF required",
                ));
            }
            forces.push((dof, start));
        }
        Ok(Self {
            initial_u,
            prescribed_end,
            force_intervals,
            prescribed,
            forces,
        })
    }

    fn at(&mut self, fraction: f64, cx: &Cx<'_>) -> Result<(), BrittleFaceApproachError> {
        for (entry, &(dof, end)) in self.prescribed.iter_mut().zip(self.prescribed_end) {
            checkpoint(cx)?;
            entry.1 = interpolate(self.initial_u[dof], end, fraction);
        }
        for (entry, &(_, [start, end])) in self.forces.iter_mut().zip(self.force_intervals) {
            checkpoint(cx)?;
            // Each trial begins at the original interval start; affine force
            // means are sampled at half its actual endpoint fraction.
            entry.1 = interpolate(start, end, 0.5 * fraction);
        }
        Ok(())
    }
}
