//! Support-driven quasistatic approach of open tensile-fracture faces.
//!
//! Every candidate replays the original `advance_path` from the same accepted
//! history. The native path owns equilibrium, damage evolution, physical work
//! and adaptive subdivision. Its endpoint hook imposes nonnegative openings
//! at EVERY represented endpoint, even when the original face would permit a
//! negative roundoff allowance. Affine P1 face openings then remain nonnegative
//! on each represented displacement segment. No positions or histories are
//! interpolated to manufacture a contacting state.
//!
//! Bisection brackets sampled native path admissibility. Nonlinear equilibrium
//! can close and reopen between load samples; neither a globally earliest
//! continuous contact nor a monotone equilibrium branch is certified. Already
//! active pairs require the existing contact owner. Unsupported free fragments
//! need dynamics: solver, material and work refusals are not treated as contact.

use fs_exec::Cx;

use crate::brittle::solid::fracture::{ReleasedBrittleFace, require_released_faces};
use crate::brittle::solid::path::{
    BrittleTetPathControls, BrittleTetPathError, BrittleTetPathStep, interpolate,
};
use crate::brittle::solid::{BrittleTetError, BrittleTetSolid, checkpoint};
use crate::plastic3::IncrementSettings;

/// Native path controls and finite sampled-contact localization allowances.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BrittleFacePathApproachControls {
    pub path: BrittleTetPathControls,
    /// Positive maximum gap of a returned contact witness, metres. Use an
    /// allowance no larger than the following contact solve's activation limit.
    pub activation_tolerance_m: f64,
    /// Positive maximum accepted/unsafe fraction separation, at most one.
    pub fraction_tolerance: f64,
    /// Every attempted native prefix counts, including rejected closure paths.
    /// At most this many times `path.maximum_trials` native solves are possible.
    pub maximum_prefix_trials: usize,
}

/// Actual native endpoint witness; not an interpolated zero-gap configuration.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BrittleFacePathContactLocation {
    pub face_id: u64,
    pub vertex: usize,
    pub gap_m: f64,
    /// Accepted prefix of the original requested affine load interval.
    pub fraction: f64,
    /// Sampled unsafe prefix, or `fraction` when the full endpoint was safe.
    /// Equality is not a claim that the physical contact load is exact.
    pub bracket_end_fraction: f64,
}

/// One complete native prefix with its original work and irreversible history.
#[derive(Debug)]
pub struct BrittleFacePathApproachStep<'a> {
    pub advance: BrittleTetPathStep<'a>,
    /// A missing contact means the entire requested interval was accepted.
    pub contact: Option<BrittleFacePathContactLocation>,
    /// Multiply each native path record's fraction by this value to express it
    /// in the original requested interval. Native records remain unmodified.
    pub fraction: f64,
    pub prefix_trials: usize,
}

/// Refusal publishes neither a partially advanced path nor a projected contact.
#[derive(Clone, Debug, PartialEq)]
pub enum BrittleFacePathApproachError {
    Invalid(&'static str),
    Path(BrittleTetPathError),
    AlreadyActive {
        face: u64,
        vertex: usize,
        gap_m: f64,
    },
    WorkLimit {
        prefix_trials: usize,
        lower_fraction: f64,
        upper_fraction: f64,
    },
    Unresolved {
        lower_fraction: f64,
        upper_fraction: f64,
    },
}
impl core::fmt::Display for BrittleFacePathApproachError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "quasistatic released-face approach refused: {self:?}")
    }
}
impl std::error::Error for BrittleFacePathApproachError {}
impl From<BrittleTetPathError> for BrittleFacePathApproachError {
    fn from(error: BrittleTetPathError) -> Self {
        Self::Path(error)
    }
}
impl From<BrittleTetError> for BrittleFacePathApproachError {
    fn from(error: BrittleTetError) -> Self {
        Self::Path(error.into())
    }
}

impl<'a> BrittleTetSolid<'a> {
    /// Resolve a nonpenetrating equilibrium prefix toward an open crack's
    /// contact, retaining the same source-bound native path and work receipts.
    ///
    /// Supports interpolate from the accepted displacement to their requested
    /// targets. Forces interpolate between their explicit endpoint values.
    /// Every trial replays the SAME accepted history; no discarded prefix can
    /// change damage. The native path checks initial force equilibrium, applies
    /// its damage/work subdivision and enforces strict openings at each
    /// accepted endpoint. Its endpoint-work budget is scaled by the actual
    /// prefix fraction of the original requested interval.
    ///
    /// Only a finite negative released-face opening narrows the upper bracket.
    /// Other native refusals propagate unchanged. Success returns either the
    /// whole open interval or a nonnegative contact witness within both caller
    /// allowances. Pass `advance.solid` to `advance_released_face_equilibrium`
    /// for subsequent touching-face loads; that owner still checks its own
    /// geometry and contact-work tolerances. It is never bypassed here.
    ///
    /// All initially non-tied pairs must be farther apart than the activation
    /// allowance. Already active constraints need contact reactions, not this
    /// unconstrained path. Loading reversals require separate intervals.
    pub fn advance_toward_released_face_contact(
        &self,
        prescribed_end_m: &[(usize, f64)],
        force_intervals_n: &[(usize, [f64; 2])],
        settings: IncrementSettings,
        controls: BrittleFacePathApproachControls,
        cx: &Cx<'_>,
    ) -> Result<BrittleFacePathApproachStep<'a>, BrittleFacePathApproachError> {
        use BrittleFacePathApproachError as E;
        checkpoint(cx)?;
        if !controls.activation_tolerance_m.is_finite()
            || controls.activation_tolerance_m <= 0.0
            || !controls.fraction_tolerance.is_finite()
            || !(0.0 < controls.fraction_tolerance && controls.fraction_tolerance <= 1.0)
            || controls.maximum_prefix_trials == 0
            || controls
                .maximum_prefix_trials
                .checked_mul(controls.path.maximum_trials)
                .is_none()
        {
            return Err(E::Invalid(
                "positive finite contact allowances and bounded native work required",
            ));
        }
        // A stricter observation policy does not alter any retained face or
        // constitutive source. The original opening kernel checks every row.
        let mut strict_faces = buffer(self.released_faces.len())?;
        let mut has_pairs = false;
        for face in self.released_faces.iter() {
            checkpoint(cx)?;
            strict_faces.push(ReleasedBrittleFace {
                opening_tolerance_m: 0.0,
                ..*face
            });
            has_pairs |= (0..3).any(|i| face.paired_nodes[0][i] != face.paired_nodes[1][i]);
        }
        let scan = |solid: &BrittleTetSolid<'a>, fraction, upper| {
            contact_witness(solid, controls.activation_tolerance_m, fraction, upper, cx)
        };
        if let Some(contact) = scan(self, 0.0, 0.0)? {
            return Err(E::AlreadyActive {
                face: contact.face_id,
                vertex: contact.vertex,
                gap_m: contact.gap_m,
            });
        }
        if !has_pairs {
            return Err(E::Invalid("open non-tied released faces required"));
        }
        let old = self.solution.displacement_m.as_flattened();
        if prescribed_end_m.len() > old.len() || force_intervals_n.len() > old.len() {
            return Err(E::Invalid("load entries exceed the original DOF count"));
        }
        let mut prescribed = buffer(prescribed_end_m.len())?;
        for &(dof, target) in prescribed_end_m {
            checkpoint(cx)?;
            if dof >= old.len() || !target.is_finite() {
                return Err(E::Invalid(
                    "finite prescribed targets on original coordinates required",
                ));
            }
            prescribed.push((dof, old[dof]));
        }
        let mut forces = buffer(force_intervals_n.len())?;
        for &(dof, values) in force_intervals_n {
            checkpoint(cx)?;
            if dof >= old.len() || values.iter().any(|v| !v.is_finite()) {
                return Err(E::Invalid(
                    "finite force endpoints on original coordinates required",
                ));
            }
            forces.push((dof, [values[0]; 2]));
        }
        let mut lower = 0.0;
        let mut upper = 1.0;
        let mut target = 1.0;
        let mut bracketed = false;
        let mut safe = None;
        let mut witness = None;
        let mut trials = 0;
        loop {
            checkpoint(cx)?;
            if trials == controls.maximum_prefix_trials {
                return Err(E::WorkLimit {
                    prefix_trials: trials,
                    lower_fraction: lower,
                    upper_fraction: upper,
                });
            }
            for (index, (entry, &(dof, end))) in
                prescribed.iter_mut().zip(prescribed_end_m).enumerate()
            {
                if index % 256 == 0 {
                    checkpoint(cx)?;
                }
                entry.1 = interpolate(old[dof], end, target);
            }
            for (index, (entry, &(_, [start, end]))) in
                forces.iter_mut().zip(force_intervals_n).enumerate()
            {
                if index % 256 == 0 {
                    checkpoint(cx)?;
                }
                entry.1 = [start, interpolate(start, end, target)];
            }
            let mut path = controls.path;
            path.work.maximum_endpoint_work_defect_j *= target;
            trials += 1;
            let result = self.advance_path_admitted(
                &prescribed,
                &forces,
                settings,
                path,
                cx,
                &mut |candidate| {
                    require_released_faces(
                        candidate.solution.displacement_m.len(),
                        &strict_faces,
                        &candidate.solution.displacement_m,
                        cx,
                    )
                },
            );
            match result {
                Ok(step) => {
                    let contact = scan(&step.solid, target, target)?;
                    if !bracketed {
                        return Ok(BrittleFacePathApproachStep {
                            advance: step,
                            contact,
                            fraction: target,
                            prefix_trials: trials,
                        });
                    }
                    lower = target;
                    safe = Some(step);
                    witness = contact;
                }
                Err(BrittleTetPathError::Endpoint(BrittleTetError::ReleasedFaceClosure {
                    opening_m,
                    ..
                })) if opening_m.is_finite() && opening_m < 0.0 => {
                    bracketed = true;
                    upper = target;
                }
                Err(error) => return Err(error.into()),
            }
            if upper - lower <= controls.fraction_tolerance
                && let Some(mut contact) = witness
            {
                contact.bracket_end_fraction = upper;
                return Ok(BrittleFacePathApproachStep {
                    advance: safe.expect("contact witness belongs to a native accepted prefix"),
                    contact: Some(contact),
                    fraction: lower,
                    prefix_trials: trials,
                });
            }
            target = lower + 0.5 * (upper - lower);
            if target <= lower || target >= upper {
                return Err(E::Unresolved {
                    lower_fraction: lower,
                    upper_fraction: upper,
                });
            }
        }
    }
}

fn buffer<T>(capacity: usize) -> Result<Vec<T>, BrittleFacePathApproachError> {
    let mut buffer = Vec::new();
    buffer
        .try_reserve_exact(capacity)
        .map_err(|_| BrittleTetPathError::Allocation)?;
    Ok(buffer)
}

fn contact_witness(
    solid: &BrittleTetSolid<'_>,
    tolerance_m: f64,
    fraction: f64,
    bracket_end_fraction: f64,
    cx: &Cx<'_>,
) -> Result<Option<BrittleFacePathContactLocation>, BrittleTetError> {
    let mut first = None;
    for face in solid.released_faces.iter() {
        checkpoint(cx)?;
        for vertex in 0..3 {
            let [a, b] = [face.paired_nodes[0][vertex], face.paired_nodes[1][vertex]];
            if a == b {
                continue;
            }
            let mut gap_m = 0.0;
            for axis in 0..3 {
                gap_m = face.normal[axis].mul_add(
                    solid.solution.displacement_m[b][axis] - solid.solution.displacement_m[a][axis],
                    gap_m,
                );
            }
            if !gap_m.is_finite() || gap_m < 0.0 {
                return Err(BrittleTetError::ReleasedFaceClosure {
                    face: face.id,
                    vertex,
                    opening_m: gap_m,
                    tolerance_m: 0.0,
                });
            }
            if first.is_none() && gap_m <= tolerance_m {
                first = Some(BrittleFacePathContactLocation {
                    face_id: face.id,
                    vertex,
                    gap_m,
                    fraction,
                    bracket_end_fraction,
                });
            }
        }
    }
    Ok(first)
}
