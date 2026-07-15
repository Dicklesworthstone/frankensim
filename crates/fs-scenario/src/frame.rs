//! Reference frames, including MOVING frames (rotating vessels, tilt
//! schedules, accelerating platforms). Frame poses are fs-ga MOTORS —
//! rigid motions with no gimbal coordinates — and frame chains compose by
//! motor multiplication, checked for cycles and dangling parents at
//! validation time (an inconsistent frame chain is a scenario bug, not a
//! solver crash).

use crate::ScenarioError;
use crate::scenario::{ValidationError, Violation, sort_validation_index};
use crate::signal::TimeSignal;
use fs_ga::{Motor, Quat, Vec3};
use fs_qty::{Dims, QtyAny};

/// A frame identity (0 is the world frame).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct FrameId(pub u32);

/// The world (inertial root) frame.
pub const WORLD: FrameId = FrameId(0);

/// Angular-rate dimensions (rad/s ⇒ s⁻¹ in SI exponents).
pub const RATE_DIMS: Dims = Dims([0, 0, -1, 0, 0, 0]);

/// How a frame moves relative to its parent.
#[derive(Debug, Clone, PartialEq)]
pub enum FrameMotion {
    /// A constant rigid offset.
    Fixed {
        /// Orientation relative to the parent.
        orientation: Quat,
        /// Translation relative to the parent (m).
        translation: Vec3,
    },
    /// Steady rotation about an axis through a center point.
    Rotating {
        /// Rotation axis (unit, in the parent frame).
        axis: [f64; 3],
        /// A point on the axis (m, parent frame).
        center: Vec3,
        /// Angular rate (rad/s).
        rate: QtyAny,
    },
    /// A scheduled tilt: angle(t) about an axis through a center — the
    /// vessel pour `(ramp 0deg 65deg 3s)` lowers to this.
    Tilt {
        /// Rotation axis (unit, in the parent frame).
        axis: [f64; 3],
        /// A point on the axis (m, parent frame).
        center: Vec3,
        /// The tilt angle schedule (radians ⇒ dimensionless).
        angle: TimeSignal,
    },
}

/// One frame in the tree.
#[derive(Debug, Clone, PartialEq)]
pub struct Frame {
    /// Identity (must be nonzero; 0 is the world).
    pub id: FrameId,
    /// Human/IR name.
    pub name: String,
    /// Parent frame (WORLD terminates every valid chain).
    pub parent: FrameId,
    /// Motion relative to the parent.
    pub motion: FrameMotion,
}

/// The scenario's frame tree (the world frame is implicit).
#[derive(Debug, Clone, Default, PartialEq)]
pub struct FrameTree {
    /// All non-world frames.
    pub frames: Vec<Frame>,
}

fn rotation_about(axis: [f64; 3], center: Vec3, angle: f64) -> Motor {
    let rot = Motor::rotor(axis, angle);
    let to_center = Motor::translator(center.x, center.y, center.z);
    let back = Motor::translator(-center.x, -center.y, -center.z);
    to_center.compose(&rot).compose(&back)
}

fn vec3_is_finite(vector: Vec3) -> bool {
    vector.x.is_finite() && vector.y.is_finite() && vector.z.is_finite()
}

fn axis_squared_norm(axis: &[f64; 3]) -> f64 {
    axis[0] * axis[0] + axis[1] * axis[1] + axis[2] * axis[2]
}

fn axis_is_unit(axis: &[f64; 3]) -> bool {
    let squared_norm = axis_squared_norm(axis);
    axis.iter().all(|component| component.is_finite())
        && squared_norm.is_finite()
        && (squared_norm - 1.0).abs() <= 1e-9
}

fn quat_squared_norm(quaternion: Quat) -> f64 {
    quaternion.w * quaternion.w
        + quaternion.x * quaternion.x
        + quaternion.y * quaternion.y
        + quaternion.z * quaternion.z
}

fn quat_is_unit(quaternion: Quat) -> bool {
    let squared_norm = quat_squared_norm(quaternion);
    [quaternion.w, quaternion.x, quaternion.y, quaternion.z]
        .iter()
        .all(|component| component.is_finite())
        && squared_norm.is_finite()
        && (squared_norm - 1.0).abs() <= 1e-9
}

fn motor_is_finite(motor: &Motor) -> bool {
    motor.0.0.iter().all(|coefficient| coefficient.is_finite())
}

fn reserve_frame_validation<T>(
    values: &mut Vec<T>,
    requested: usize,
    resource: &'static str,
) -> Result<(), ValidationError> {
    values
        .try_reserve_exact(requested)
        .map_err(|_| ValidationError::AllocationRefused {
            resource,
            requested,
        })
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum FrameIdLookup {
    Missing,
    Unique(usize),
    Ambiguous(usize),
}

fn frame_id_lookup(index: &[(u32, usize)], id: u32) -> FrameIdLookup {
    let position = index.partition_point(|(candidate, _)| *candidate < id);
    let Some((candidate, row)) = index.get(position) else {
        return FrameIdLookup::Missing;
    };
    if *candidate != id {
        return FrameIdLookup::Missing;
    }
    if index.get(position + 1).is_some_and(|(next, _)| *next == id) {
        FrameIdLookup::Ambiguous(*row)
    } else {
        FrameIdLookup::Unique(*row)
    }
}

fn first_frame_id_row(index: &[(u32, usize)], id: u32) -> Option<usize> {
    match frame_id_lookup(index, id) {
        FrameIdLookup::Missing => None,
        FrameIdLookup::Unique(row) | FrameIdLookup::Ambiguous(row) => Some(row),
    }
}

fn first_frame_name_row(index: &[(&str, usize)], name: &str) -> Option<usize> {
    let position = index.partition_point(|(candidate, _)| *candidate < name);
    index
        .get(position)
        .and_then(|(candidate, row)| (*candidate == name).then_some(*row))
}

impl FrameTree {
    /// An empty tree (world only).
    #[must_use]
    pub fn new() -> Self {
        FrameTree { frames: Vec::new() }
    }

    /// Add a frame.
    pub fn add(&mut self, frame: Frame) {
        self.frames.push(frame);
    }

    fn find(&self, id: FrameId) -> Option<&Frame> {
        self.frames.iter().find(|f| f.id == id)
    }

    fn find_unique(&self, id: FrameId) -> Result<&Frame, ScenarioError> {
        let mut matches = self.frames.iter().filter(|frame| frame.id == id);
        let frame = matches.next().ok_or_else(|| ScenarioError::Frame {
            what: format!("frame id {} not found", id.0),
        })?;
        if matches.next().is_some() {
            return Err(ScenarioError::Frame {
                what: format!(
                    "frame id {} is ambiguous because it is defined more than once",
                    id.0
                ),
            });
        }
        Ok(frame)
    }

    /// The motor mapping this frame's coordinates into its PARENT's at
    /// time `t` (seconds).
    ///
    /// # Errors
    /// [`ScenarioError`] for bad schedules or dimension defects.
    pub fn local_pose(frame: &Frame, t: f64) -> Result<Motor, ScenarioError> {
        if !t.is_finite() {
            return Err(ScenarioError::Evaluate {
                what: format!("frame {:?}: non-finite evaluation time {t}", frame.name),
            });
        }
        let motor = match &frame.motion {
            FrameMotion::Fixed {
                orientation,
                translation,
            } => {
                if !quat_is_unit(*orientation) {
                    return Err(ScenarioError::Frame {
                        what: format!(
                            "frame {:?}: fixed orientation is non-finite or not unit length",
                            frame.name
                        ),
                    });
                }
                if !vec3_is_finite(*translation) {
                    return Err(ScenarioError::Frame {
                        what: format!("frame {:?}: fixed translation is non-finite", frame.name),
                    });
                }
                Motor::from_parts(*orientation, *translation)
            }
            FrameMotion::Rotating { axis, center, rate } => {
                if !axis_is_unit(axis) {
                    return Err(ScenarioError::Frame {
                        what: format!(
                            "frame {:?}: rotating axis is non-finite or not unit length",
                            frame.name
                        ),
                    });
                }
                if !vec3_is_finite(*center) {
                    return Err(ScenarioError::Frame {
                        what: format!("frame {:?}: rotating center is non-finite", frame.name),
                    });
                }
                if rate.dims != RATE_DIMS {
                    return Err(ScenarioError::Dimensions {
                        context: format!("frame {:?} angular rate", frame.name),
                        expected: RATE_DIMS.0,
                        got: rate.dims.0,
                    });
                }
                if !rate.value.is_finite() {
                    return Err(ScenarioError::Evaluate {
                        what: format!("frame {:?}: angular rate is non-finite", frame.name),
                    });
                }
                let angle = rate.value * t;
                if !angle.is_finite() {
                    return Err(ScenarioError::Evaluate {
                        what: format!(
                            "frame {:?}: angular-rate/time product overflowed",
                            frame.name
                        ),
                    });
                }
                rotation_about(*axis, *center, angle)
            }
            FrameMotion::Tilt {
                axis,
                center,
                angle,
            } => {
                if !axis_is_unit(axis) {
                    return Err(ScenarioError::Frame {
                        what: format!(
                            "frame {:?}: tilt axis is non-finite or not unit length",
                            frame.name
                        ),
                    });
                }
                if !vec3_is_finite(*center) {
                    return Err(ScenarioError::Frame {
                        what: format!("frame {:?}: tilt center is non-finite", frame.name),
                    });
                }
                let a = angle.eval(t)?;
                if !a.dims.is_none() {
                    return Err(ScenarioError::Dimensions {
                        context: format!("frame {:?} tilt angle", frame.name),
                        expected: Dims::NONE.0,
                        got: a.dims.0,
                    });
                }
                rotation_about(*axis, *center, a.value)
            }
        };
        if !motor_is_finite(&motor) {
            return Err(ScenarioError::Frame {
                what: format!(
                    "frame {:?}: finite inputs produced non-finite motor coefficients",
                    frame.name
                ),
            });
        }
        Ok(motor)
    }

    /// The motor mapping `id`'s coordinates into WORLD coordinates at
    /// time `t`, composed down the parent chain.
    ///
    /// # Errors
    /// [`ScenarioError::Frame`] on unresolvable or cyclic chains.
    pub fn world_pose(&self, id: FrameId, t: f64) -> Result<Motor, ScenarioError> {
        if !t.is_finite() {
            return Err(ScenarioError::Evaluate {
                what: format!("frame id {}: non-finite evaluation time {t}", id.0),
            });
        }
        let mut pose = Motor::identity();
        let mut current = id;
        let mut hops = 0usize;
        while current != WORLD {
            let frame = self.find_unique(current)?;
            pose = Self::local_pose(frame, t)?.compose(&pose);
            if !motor_is_finite(&pose) {
                return Err(ScenarioError::Frame {
                    what: format!(
                        "frame id {}: parent-chain composition produced non-finite motor coefficients",
                        id.0
                    ),
                });
            }
            current = frame.parent;
            hops += 1;
            if hops > self.frames.len() {
                return Err(ScenarioError::Frame {
                    what: format!("cyclic frame chain reached from id {}", id.0),
                });
            }
        }
        Ok(pose)
    }

    fn validation_cycle_flags(
        &self,
        first_by_id: &[(u32, usize)],
        checkpoint: &mut impl FnMut(&'static str) -> Result<(), ValidationError>,
    ) -> Result<Vec<bool>, ValidationError> {
        // Each frame has at most one parent, so a tri-color walk visits every
        // storage row at most once. Nodes leading into a cycle are marked too:
        // their parent chain also never reaches WORLD.
        let frame_count = self.frames.len();
        let mut state = Vec::new();
        reserve_frame_validation(&mut state, frame_count, "frame cycle states")?;
        for _ in 0..frame_count {
            checkpoint("frame cycle state initialization")?;
            state.push(0u8);
        }
        let mut reaches_cycle = Vec::new();
        reserve_frame_validation(&mut reaches_cycle, frame_count, "frame cycle result flags")?;
        for _ in 0..frame_count {
            checkpoint("frame cycle result initialization")?;
            reaches_cycle.push(false);
        }
        let mut path = Vec::new();
        reserve_frame_validation(&mut path, frame_count, "frame cycle path")?;
        for start in 0..frame_count {
            checkpoint("frame cycle start")?;
            if state[start] != 0 {
                continue;
            }
            path.clear();
            let mut current = Some(start);
            let cyclic = loop {
                checkpoint("frame cycle traversal")?;
                let Some(index) = current else {
                    break false;
                };
                match state[index] {
                    0 => {
                        state[index] = 1;
                        path.push(index);
                        let frame = &self.frames[index];
                        current = if frame.id == WORLD || frame.parent == WORLD {
                            None
                        } else {
                            match frame_id_lookup(first_by_id, frame.parent.0) {
                                FrameIdLookup::Unique(parent) => Some(parent),
                                FrameIdLookup::Missing | FrameIdLookup::Ambiguous(_) => None,
                            }
                        };
                    }
                    1 => break true,
                    2 => break reaches_cycle[index],
                    _ => unreachable!("frame validation color is internal"),
                }
            };
            for index in path.drain(..).rev() {
                checkpoint("frame cycle finalization")?;
                reaches_cycle[index] = cyclic;
                state[index] = 2;
            }
        }
        Ok(reaches_cycle)
    }

    /// Structural validation: nonempty unique names, unique nonzero ids,
    /// resolvable acyclic parents, unit axes, and angle-schedule dimensions.
    ///
    /// Identity lookup is indexed once and cycle detection is a tri-color
    /// parent walk, avoiding the former repeated linear scans per chain.
    pub fn check(&self, out: &mut Vec<Violation>) {
        let mut checkpoint = |_: &'static str| Ok(());
        if let Err(error) = self.check_with_checkpoint(out, &mut checkpoint) {
            out.push(error.into_violation());
        }
    }

    pub(crate) fn check_with_checkpoint(
        &self,
        out: &mut Vec<Violation>,
        checkpoint: &mut impl FnMut(&'static str) -> Result<(), ValidationError>,
    ) -> Result<(), ValidationError> {
        let frame_count = self.frames.len();
        let mut first_by_id = Vec::new();
        reserve_frame_validation(&mut first_by_id, frame_count, "frame id index")?;
        let mut first_by_name = Vec::new();
        reserve_frame_validation(&mut first_by_name, frame_count, "frame name index")?;
        for (index, frame) in self.frames.iter().enumerate() {
            checkpoint("frame index")?;
            first_by_id.push((frame.id.0, index));
            first_by_name.push((frame.name.as_str(), index));
        }
        sort_validation_index(&mut first_by_id, "frame id index sort", checkpoint)?;
        sort_validation_index(&mut first_by_name, "frame name index sort", checkpoint)?;
        let reaches_cycle = self.validation_cycle_flags(&first_by_id, checkpoint)?;

        for (i, f) in self.frames.iter().enumerate() {
            checkpoint("frame validation")?;
            let ctx = format!("frame {:?}", f.name);
            if f.id == WORLD {
                out.push(Violation {
                    code: "frame-id-zero",
                    what: format!("{ctx}: id 0 is reserved for the world frame"),
                    fix: "renumber the frame with a nonzero id".to_string(),
                });
            }
            if f.name.is_empty() {
                out.push(Violation {
                    code: "frame-name-empty",
                    what: format!("frame storage row {i} has an empty name"),
                    fix: "give every frame a nonempty exact UTF-8 name".to_string(),
                });
            }
            let first_id_row = first_frame_id_row(&first_by_id, f.id.0).unwrap_or(i);
            if first_id_row != i {
                out.push(Violation {
                    code: "frame-id-duplicate",
                    what: format!(
                        "{ctx}: id {} first appears at frame row {} and repeats at row {i}",
                        f.id.0, first_id_row
                    ),
                    fix: "give every frame a unique id".to_string(),
                });
            }
            let first_name_row = first_frame_name_row(&first_by_name, f.name.as_str()).unwrap_or(i);
            if first_name_row != i {
                out.push(Violation {
                    code: "frame-name-duplicate",
                    what: format!(
                        "{ctx}: name first appears at frame row {} and repeats at row {i}",
                        first_name_row
                    ),
                    fix: "give every frame a unique name".to_string(),
                });
            }
            if f.parent != WORLD {
                match frame_id_lookup(&first_by_id, f.parent.0) {
                    FrameIdLookup::Missing => out.push(Violation {
                        code: "frame-parent-missing",
                        what: format!("{ctx}: parent id {} does not exist", f.parent.0),
                        fix: "point the frame at an existing parent or at the world (0)"
                            .to_string(),
                    }),
                    FrameIdLookup::Ambiguous(first) => out.push(Violation {
                        code: "frame-parent-ambiguous",
                        what: format!(
                            "{ctx}: parent id {} is defined more than once, first at frame row {first}",
                            f.parent.0
                        ),
                        fix: "give every frame a unique id before resolving parent chains"
                            .to_string(),
                    }),
                    FrameIdLookup::Unique(_) => {}
                }
            }
            match &f.motion {
                FrameMotion::Fixed {
                    orientation,
                    translation,
                } => {
                    if !quat_is_unit(*orientation) {
                        out.push(Violation {
                            code: "frame-orientation-invalid",
                            what: format!(
                                "{ctx}: fixed orientation has squared norm {} or a non-finite component",
                                quat_squared_norm(*orientation)
                            ),
                            fix: "supply a finite unit quaternion".to_string(),
                        });
                    }
                    if !vec3_is_finite(*translation) {
                        out.push(Violation {
                            code: "frame-translation-nonfinite",
                            what: format!("{ctx}: fixed translation is non-finite"),
                            fix: "replace every translation component with a finite value"
                                .to_string(),
                        });
                    }
                }
                FrameMotion::Rotating { axis, center, rate } => {
                    check_axis(axis, &ctx, out);
                    check_center(*center, &ctx, out);
                    if rate.dims != RATE_DIMS {
                        out.push(dims_violation(&ctx, "angular rate", RATE_DIMS, rate.dims));
                    }
                    if !rate.value.is_finite() {
                        out.push(Violation {
                            code: "frame-rate-nonfinite",
                            what: format!("{ctx}: angular rate {} is non-finite", rate.value),
                            fix: "replace the angular rate with a finite value".to_string(),
                        });
                    }
                }
                FrameMotion::Tilt {
                    axis,
                    center,
                    angle,
                } => {
                    check_axis(axis, &ctx, out);
                    check_center(*center, &ctx, out);
                    angle.check_with_checkpoint(&ctx, out, checkpoint)?;
                    if !angle.dims().is_none() {
                        out.push(dims_violation(&ctx, "tilt angle", Dims::NONE, angle.dims()));
                    }
                }
            }
            if reaches_cycle[i] {
                out.push(Violation {
                    code: "frame-chain-cyclic",
                    what: format!("{ctx}: parent chain never reaches the world frame"),
                    fix: "break the cycle so every chain terminates at frame 0".to_string(),
                });
            }
            checkpoint("frame validation")?;
        }
        Ok(())
    }
}

fn check_axis(axis: &[f64; 3], ctx: &str, out: &mut Vec<Violation>) {
    let squared_norm = axis_squared_norm(axis);
    if !axis_is_unit(axis) {
        out.push(Violation {
            code: "frame-axis-not-unit",
            what: format!(
                "{ctx}: rotation axis has squared norm {squared_norm} or a non-finite component"
            ),
            fix: "supply a finite rotation axis normalized to unit length".to_string(),
        });
    }
}

fn check_center(center: Vec3, ctx: &str, out: &mut Vec<Violation>) {
    if !vec3_is_finite(center) {
        out.push(Violation {
            code: "frame-center-nonfinite",
            what: format!("{ctx}: rotation center is non-finite"),
            fix: "replace every center component with a finite value".to_string(),
        });
    }
}

fn dims_violation(ctx: &str, quantity: &str, expected: Dims, got: Dims) -> Violation {
    Violation {
        code: "frame-dims",
        what: format!(
            "{ctx}: {quantity} has dimensions {:?}, expected {:?}",
            got.0, expected.0
        ),
        fix: format!("express the {quantity} in coherent SI units"),
    }
}

#[cfg(test)]
mod validation_internal_tests {
    use super::{Frame, FrameId, FrameMotion, FrameTree, WORLD, reserve_frame_validation};
    use crate::scenario::ValidationError;
    use fs_ga::{Quat, Vec3};

    fn fixed_frame(id: u32, parent: FrameId) -> Frame {
        Frame {
            id: FrameId(id),
            name: format!("frame-{id}"),
            parent,
            motion: FrameMotion::Fixed {
                orientation: Quat::identity(),
                translation: Vec3::new(0.0, 0.0, 0.0),
            },
        }
    }

    #[test]
    fn injected_cancellation_reaches_the_frame_cycle_walk() {
        let mut tree = FrameTree::new();
        tree.add(fixed_frame(1, WORLD));
        tree.add(fixed_frame(2, FrameId(1)));
        let mut findings = Vec::new();
        let mut visited = Vec::new();

        let result = tree.check_with_checkpoint(&mut findings, &mut |phase| {
            visited.push(phase);
            if phase == "frame cycle traversal" {
                Err(ValidationError::Cancelled {
                    phase,
                    completed: 0,
                    planned: 2,
                })
            } else {
                Ok(())
            }
        });

        assert!(matches!(
            result,
            Err(ValidationError::Cancelled {
                phase: "frame cycle traversal",
                completed: 0,
                planned: 2,
            })
        ));
        assert!(findings.is_empty());
        assert_eq!(&visited[..2], ["frame index", "frame index"]);
        assert!(visited.contains(&"frame id index sort"));
        assert!(visited.contains(&"frame name index sort"));
        assert_eq!(
            &visited[visited.len() - 2..],
            ["frame cycle start", "frame cycle traversal"]
        );
    }

    #[test]
    fn frame_scratch_capacity_overflow_is_typed() {
        let mut scratch = Vec::<u8>::new();
        assert!(matches!(
            reserve_frame_validation(&mut scratch, usize::MAX, "frame test scratch"),
            Err(ValidationError::AllocationRefused {
                resource: "frame test scratch",
                requested: usize::MAX,
            })
        ));
        assert!(scratch.is_empty());
    }

    #[test]
    fn duplicate_parent_ids_are_ambiguous_under_every_declaration_order() {
        let root = Frame {
            id: FrameId(1),
            name: "root-definition".to_string(),
            parent: WORLD,
            motion: fixed_frame(1, WORLD).motion,
        };
        let back_edge = Frame {
            id: FrameId(1),
            name: "back-edge-definition".to_string(),
            parent: FrameId(2),
            motion: fixed_frame(1, WORLD).motion,
        };
        let child = Frame {
            id: FrameId(2),
            name: "child".to_string(),
            parent: FrameId(1),
            motion: fixed_frame(2, WORLD).motion,
        };

        for frames in [
            [root.clone(), back_edge.clone(), child.clone()],
            [back_edge.clone(), root.clone(), child.clone()],
        ] {
            let tree = FrameTree {
                frames: frames.into(),
            };
            let mut findings = Vec::new();
            tree.check(&mut findings);
            assert_eq!(
                findings
                    .iter()
                    .map(|finding| finding.code)
                    .collect::<Vec<_>>(),
                ["frame-id-duplicate", "frame-parent-ambiguous"]
            );
            assert!(
                findings
                    .iter()
                    .all(|finding| finding.code != "frame-chain-cyclic"),
                "an ambiguous parent must not be resolved through an arbitrary first row"
            );
        }
    }
}
