//! Reference frames, including MOVING frames (rotating vessels, tilt
//! schedules, accelerating platforms). Frame poses are fs-ga MOTORS —
//! rigid motions with no gimbal coordinates — and frame chains compose by
//! motor multiplication, checked for cycles and dangling parents at
//! validation time (an inconsistent frame chain is a scenario bug, not a
//! solver crash).

use crate::ScenarioError;
use crate::scenario::Violation;
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

    /// The motor mapping this frame's coordinates into its PARENT's at
    /// time `t` (seconds).
    ///
    /// # Errors
    /// [`ScenarioError`] for bad schedules or dimension defects.
    pub fn local_pose(frame: &Frame, t: f64) -> Result<Motor, ScenarioError> {
        match &frame.motion {
            FrameMotion::Fixed {
                orientation,
                translation,
            } => Ok(Motor::from_parts(*orientation, *translation)),
            FrameMotion::Rotating { axis, center, rate } => {
                if rate.dims != RATE_DIMS {
                    return Err(ScenarioError::Dimensions {
                        context: format!("frame {:?} angular rate", frame.name),
                        expected: RATE_DIMS.0,
                        got: rate.dims.0,
                    });
                }
                Ok(rotation_about(*axis, *center, rate.value * t))
            }
            FrameMotion::Tilt {
                axis,
                center,
                angle,
            } => {
                let a = angle.eval(t)?;
                if !a.dims.is_none() {
                    return Err(ScenarioError::Dimensions {
                        context: format!("frame {:?} tilt angle", frame.name),
                        expected: Dims::NONE.0,
                        got: a.dims.0,
                    });
                }
                Ok(rotation_about(*axis, *center, a.value))
            }
        }
    }

    /// The motor mapping `id`'s coordinates into WORLD coordinates at
    /// time `t`, composed down the parent chain.
    ///
    /// # Errors
    /// [`ScenarioError::Frame`] on unresolvable or cyclic chains.
    pub fn world_pose(&self, id: FrameId, t: f64) -> Result<Motor, ScenarioError> {
        let mut pose = Motor::identity();
        let mut current = id;
        let mut hops = 0usize;
        while current != WORLD {
            let frame = self.find(current).ok_or_else(|| ScenarioError::Frame {
                what: format!("frame id {} not found", current.0),
            })?;
            pose = Self::local_pose(frame, t)?.compose(&pose);
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

    /// Structural validation: unique nonzero ids, unique names, resolvable
    /// acyclic parents, unit axes, angle-schedule dimensions.
    pub fn check(&self, out: &mut Vec<Violation>) {
        for (i, f) in self.frames.iter().enumerate() {
            let ctx = format!("frame {:?}", f.name);
            if f.id == WORLD {
                out.push(Violation {
                    code: "frame-id-zero",
                    what: format!("{ctx}: id 0 is reserved for the world frame"),
                    fix: "renumber the frame with a nonzero id".to_string(),
                });
            }
            if self.frames[..i].iter().any(|g| g.id == f.id) {
                out.push(Violation {
                    code: "frame-id-duplicate",
                    what: format!("{ctx}: id {} is already taken", f.id.0),
                    fix: "give every frame a unique id".to_string(),
                });
            }
            if self.frames[..i].iter().any(|g| g.name == f.name) {
                out.push(Violation {
                    code: "frame-name-duplicate",
                    what: format!("{ctx}: name is already taken"),
                    fix: "give every frame a unique name".to_string(),
                });
            }
            if f.parent != WORLD && self.find(f.parent).is_none() {
                out.push(Violation {
                    code: "frame-parent-missing",
                    what: format!("{ctx}: parent id {} does not exist", f.parent.0),
                    fix: "point the frame at an existing parent or at the world (0)".to_string(),
                });
            }
            match &f.motion {
                FrameMotion::Fixed { .. } => {}
                FrameMotion::Rotating { axis, rate, .. } => {
                    check_axis(axis, &ctx, out);
                    if rate.dims != RATE_DIMS {
                        out.push(dims_violation(&ctx, "angular rate", RATE_DIMS, rate.dims));
                    }
                }
                FrameMotion::Tilt { axis, angle, .. } => {
                    check_axis(axis, &ctx, out);
                    angle.check(&ctx, out);
                    if !angle.dims().is_none() {
                        out.push(dims_violation(&ctx, "tilt angle", Dims::NONE, angle.dims()));
                    }
                }
            }
            // Cycle detection: walk the parent chain STRUCTURALLY — no pose
            // evaluation. Using `world_pose(...).is_err()` conflated a genuine
            // cycle with a `local_pose` failure (a bad motion dimension or an
            // unevaluable angle signal, e.g. an empty table), so a well-formed
            // acyclic frame with an independently-reported signal defect got a
            // spurious `frame-chain-cyclic` that misdirected the repair. A
            // dangling parent is reported separately (`frame-parent-missing`)
            // and is not a cycle here.
            if self.find(f.parent).is_some() {
                let mut current = f.id;
                let mut hops = 0usize;
                let cyclic = loop {
                    if current == WORLD {
                        break false;
                    }
                    match self.find(current) {
                        None => break false, // chain left the graph — not a cycle
                        Some(frame) => current = frame.parent,
                    }
                    hops += 1;
                    if hops > self.frames.len() {
                        break true;
                    }
                };
                if cyclic {
                    out.push(Violation {
                        code: "frame-chain-cyclic",
                        what: format!("{ctx}: parent chain never reaches the world frame"),
                        fix: "break the cycle so every chain terminates at frame 0".to_string(),
                    });
                }
            }
        }
    }
}

fn check_axis(axis: &[f64; 3], ctx: &str, out: &mut Vec<Violation>) {
    let n2 = axis[0] * axis[0] + axis[1] * axis[1] + axis[2] * axis[2];
    if (n2 - 1.0).abs() > 1e-9 {
        out.push(Violation {
            code: "frame-axis-not-unit",
            what: format!("{ctx}: rotation axis has squared norm {n2}"),
            fix: "normalize the rotation axis to unit length".to_string(),
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
