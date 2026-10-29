//! Changing rigid connectivity without losing spin or orbital momentum.
//!
//! Component COM poses and full native inertias define one composite body.
//! Release inherits its rigid velocity field; internal impulses then act at
//! shared physical points, with optional equal/opposite free couples. Latching
//! projects arbitrary component velocities onto their common rigid motion at
//! fixed poses, preserving spatial momentum and reporting the removed relative
//! kinetic energy. The caller owns the physical decision to connect/disconnect:
//! this module does not discover contact, move geometry, or partition heat.
//!
//! Inertia uses R I R^T + m (|r|² 1 - r r^T). Angular momentum is always spin
//! plus r cross p about the same point. See Drake's spatial-inertia and spatial-
//! momentum definitions: https://drake.mit.edu/doxygen_cxx/group__multibody__spatial__inertia.html.

use crate::{DynamicsError, MassProperties, Pose, RigidBodyState, UnitQuaternion, Vec3};

/// A component's original mass properties and COM pose in a common reference.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AssemblyPart {
    /// Full inertia about this component's COM, in its own body axes.
    pub properties: MassProperties,
    /// Body-to-reference orientation and COM position in reference axes.
    pub pose_in_reference: Pose,
}

/// Explicit work caps and absolute acceptance tolerances in SI units.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AssemblyControls {
    /// Maximum components; checked before allocating their staging vectors.
    pub maximum_parts: usize,
    /// Maximum internal impulse pairs in one release.
    pub maximum_impulses: usize,
    /// Maximum component of actual total linear-momentum error [N s].
    pub linear_momentum_tolerance_ns: f64,
    /// Maximum component of spin plus orbital momentum error [N m s].
    pub angular_momentum_tolerance_nms: f64,
    /// Maximum actual kinetic/work balance error [J].
    pub energy_tolerance_j: f64,
}

/// Equal/opposite impulses at one assembly-COM-relative physical point.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SeparationImpulse {
    /// Index receiving the stated force impulse and free couple.
    pub part_a: usize,
    /// Distinct index receiving their exact negatives.
    pub part_b: usize,
    /// Shared application point, relative to the assembly COM in its axes [m].
    pub point_assembly_m: Vec3,
    /// Force impulse on A, expressed in assembly axes [N s].
    pub impulse_on_a_assembly_ns: Vec3,
    /// Free angular impulse on A, excluding the point-force moment [N m s].
    pub couple_on_a_assembly_nms: Vec3,
}

/// Residuals reconstructed from actual returned native states.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AssemblyBalance {
    /// Final minus initial total world linear momentum [N s].
    pub linear_momentum_residual_world_ns: Vec3,
    /// Final minus initial spin plus orbital momentum about the event COM [N m s].
    pub angular_momentum_residual_world_nms: Vec3,
    /// Kinetic change minus supplied work, or plus removed relative energy [J].
    pub energy_residual_j: f64,
}

/// Every original component survives release as a native rigid-body state.
#[derive(Clone, Debug, PartialEq)]
pub struct AssemblyRelease {
    /// States in the original part order, each with its own body axes and COM.
    pub bodies: Vec<RigidBodyState>,
    /// Work supplied by all internal impulses together [J]. Shared-body cross
    /// terms are included; individual pair energies must not be added instead.
    pub supplied_work_j: f64,
    /// Kinetic energy before release in the parent-COM translating frame [J].
    pub kinetic_before_j: f64,
    /// Kinetic energy after release in the same inertial frame [J].
    pub kinetic_after_j: f64,
    /// Actual momentum and work accounting.
    pub balance: AssemblyBalance,
}

/// A perfectly inelastic rigid lock at unchanged component poses.
#[derive(Clone, Debug, PartialEq)]
pub struct AssemblyLatch {
    /// Composite mass distribution retaining all original component geometry.
    pub assembly: RigidAssembly,
    /// Momentum-preserving composite state at its COM.
    pub state: RigidBodyState,
    /// Nonnegative kinetic energy removed from relative motion [J], evaluated
    /// as a mass-metric squared residual. No thermal destination is assumed.
    pub removed_relative_energy_j: f64,
    /// Actual momentum and energy accounting in a common translating frame.
    pub balance: AssemblyBalance,
}

/// A cached composite inertia and its COM-relative component geometry.
#[derive(Clone, Debug, PartialEq)]
pub struct RigidAssembly {
    parts: Vec<AssemblyPart>,
    properties: MassProperties,
    center_in_reference: Vec3,
}

/// Admission, arithmetic, work-cap, cancellation, or physical-balance refusal.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum AssemblyError {
    /// Invalid input or an unrepresentable derived quantity.
    Invalid(&'static str),
    /// The original native rigid owner refused.
    Dynamics(DynamicsError),
    /// Requested components or impulses exceed their explicit cap.
    Budget {
        /// Component or impulse count.
        requested: usize,
        /// Declared limit.
        maximum: usize,
    },
    /// Actual output failed a declared balance tolerance; nothing is published.
    Balance(AssemblyBalance),
    /// Cancellation observed before final publication.
    Cancelled,
}
impl From<DynamicsError> for AssemblyError {
    fn from(error: DynamicsError) -> Self {
        Self::Dynamics(error)
    }
}
impl core::fmt::Display for AssemblyError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "rigid assembly: {self:?}")
    }
}
impl std::error::Error for AssemblyError {}

impl RigidAssembly {
    /// Assemble arbitrary oriented full inertias using the parallel-axis law.
    /// The supplied reference may have any origin; returned part positions are
    /// shifted to the computed COM, with reference axes unchanged. Assembly is
    /// O(parts), with one native tensor admission and no dense 6N matrix.
    pub fn new(
        parts: &[AssemblyPart],
        controls: AssemblyControls,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, AssemblyError> {
        poll(&mut cancelled)?;
        controls.validate(parts.len(), 0)?;
        let mut mass = Sum::default();
        for part in parts {
            poll(&mut cancelled)?;
            part.properties.validate()?;
            mass.add(part.properties.mass());
        }
        let mass = finite(mass.value())?;
        let anchor = parts[0].pose_in_reference.position_world();
        let mut offset = VectorSum::default();
        for part in parts {
            poll(&mut cancelled)?;
            offset.add(
                vector(part.pose_in_reference.position_world().sub(anchor))?
                    .scale(part.properties.mass() / mass),
            );
        }
        let offset = vector(offset.value())?;
        let center = vector(anchor.add(offset))?;
        let mut moments = [Sum::default(); 6];
        let mut centered = Vec::with_capacity(parts.len());
        for part in parts {
            poll(&mut cancelled)?;
            // Subtract the small reference-relative offset before adding a
            // potentially large absolute reference origin back anywhere.
            let r = vector(
                part.pose_in_reference
                    .position_world()
                    .sub(anchor)
                    .sub(offset),
            )?;
            let q = compose(
                part.pose_in_reference.orientation(),
                part.properties.principal_axes_to_body(),
            )?;
            let axes = [
                array(q.rotate_body_to_world(Vec3::new(1.0, 0.0, 0.0))),
                array(q.rotate_body_to_world(Vec3::new(0.0, 1.0, 0.0))),
                array(q.rotate_body_to_world(Vec3::new(0.0, 0.0, 1.0))),
            ];
            let principal = array(part.properties.principal_inertia_body());
            let u = array(vector(r.scale(part.properties.mass().sqrt()))?);
            for (k, (row, col)) in ENTRIES.into_iter().enumerate() {
                let rotated = principal[0] * axes[0][row] * axes[0][col]
                    + principal[1] * axes[1][row] * axes[1][col]
                    + principal[2] * axes[2][row] * axes[2][col];
                let shifted = if row == col {
                    u[(row + 1) % 3].powi(2) + u[(row + 2) % 3].powi(2)
                } else {
                    -u[row] * u[col]
                };
                moments[k].add(finite(rotated + shifted)?);
            }
            centered.push(AssemblyPart {
                properties: part.properties,
                pose_in_reference: Pose::new(r, part.pose_in_reference.orientation())?,
            });
        }
        let mut inertia = [0.0; 9];
        for (k, (row, col)) in ENTRIES.into_iter().enumerate() {
            let value = finite(moments[k].value())?;
            inertia[3 * row + col] = value;
            inertia[3 * col + row] = value;
        }
        let properties =
            MassProperties::from_inertia_tensor(mass, Vec3::ZERO, fs_ga::Mat3 { m: inertia })?;
        poll(&mut cancelled)?;
        Ok(Self {
            parts: centered,
            properties,
            center_in_reference: center,
        })
    }

    /// Composite mass and full inertia about the assembly COM in reference axes.
    pub const fn properties(&self) -> MassProperties {
        self.properties
    }

    /// Original reference-coordinate location of the computed composite COM.
    pub const fn center_of_mass_in_reference(&self) -> Vec3 {
        self.center_in_reference
    }

    /// Original part order with COM-relative positions and unchanged body axes.
    pub fn parts(&self) -> &[AssemblyPart] {
        &self.parts
    }

    /// Detach all parts, optionally applying prescribed internal impulse pairs.
    /// Parent pose is the assembly COM and its orientation maps assembly axes
    /// to world. Every part inherits v_com + omega cross r and the same world
    /// angular velocity before recoil. Inputs are borrowed and never mutated.
    pub fn release(
        &self,
        parent: RigidBodyState,
        impulses: &[SeparationImpulse],
        controls: AssemblyControls,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<AssemblyRelease, AssemblyError> {
        poll(&mut cancelled)?;
        controls.validate(self.parts.len(), impulses.len())?;
        let mut bodies = vec![parent; self.parts.len()];
        let mut deltas = vec![ImpulseTotals::default(); self.parts.len()];
        let result = release_components_into(
            &self.parts,
            self.properties,
            parent,
            impulses,
            &mut deltas,
            &mut bodies,
            &mut cancelled,
        )?;
        controls.admit(result.balance)?;
        poll(&mut cancelled)?;
        Ok(AssemblyRelease {
            bodies,
            supplied_work_j: result.supplied_work_j,
            kinetic_before_j: result.kinetic_before_j,
            kinetic_after_j: result.kinetic_after_j,
            balance: result.balance,
        })
    }

    /// Lock all relative motion instantaneously at unchanged world poses.
    /// The chosen assembly orientation affects coordinates only. This is the
    /// perfectly inelastic rigid-connection law, not a contact/restitution or
    /// finite-duration attachment model. Geometry/contact admissibility and
    /// the destination of removed mechanical energy remain with the caller.
    pub fn latch(
        bodies: &[RigidBodyState],
        properties: &[MassProperties],
        orientation: UnitQuaternion,
        controls: AssemblyControls,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<AssemblyLatch, AssemblyError> {
        poll(&mut cancelled)?;
        controls.validate(bodies.len(), 0)?;
        if bodies.len() != properties.len() {
            return Err(AssemblyError::Invalid(
                "matching states and mass properties required",
            ));
        }
        let anchor = bodies[0].pose().position_world();
        let inverse = orientation.conjugate();
        let mut parts = Vec::with_capacity(bodies.len());
        for (body, property) in bodies.iter().zip(properties) {
            poll(&mut cancelled)?;
            parts.push(AssemblyPart {
                properties: *property,
                pose_in_reference: Pose::new(
                    inverse.rotate_body_to_world(vector(body.pose().position_world().sub(anchor))?),
                    compose(inverse, body.pose().orientation())?,
                )?,
            });
        }
        let assembly = Self::new(&parts, controls, &mut cancelled)?;
        let center =
            vector(anchor.add(orientation.rotate_body_to_world(assembly.center_in_reference)))?;
        let boost = bodies[0].center_of_mass_velocity_world(properties[0])?;
        let before = totals(
            bodies,
            properties.iter().copied(),
            center,
            boost,
            &mut cancelled,
        )?;
        let state = RigidBodyState::new(
            Pose::new(center, orientation)?,
            before.linear_world,
            orientation.rotate_world_to_body(before.angular_world),
        )?;
        let after = totals(
            &[state],
            [assembly.properties].into_iter(),
            center,
            boost,
            &mut cancelled,
        )?;
        let velocity_relative = state
            .linear_momentum_world()
            .sub(boost.scale(assembly.properties.mass()))
            .scale(assembly.properties.mass().recip());
        let omega = orientation.rotate_body_to_world(
            assembly
                .properties
                .angular_velocity_body_checked(state.angular_momentum_body())?,
        );
        let mut removed = Sum::default();
        for (body, property) in bodies.iter().zip(properties) {
            poll(&mut cancelled)?;
            let arm = body.pose().position_world().sub(center);
            let relative = body
                .linear_momentum_world()
                .sub(boost.scale(property.mass()))
                .scale(property.mass().recip());
            let dv = relative.sub(velocity_relative.add(omega.cross(arm)));
            let dw = property
                .angular_velocity_body_checked(body.angular_momentum_body())?
                .sub(body.pose().orientation().rotate_world_to_body(omega));
            removed.add(finite(
                0.5 * (property.mass() * dv.dot(dv)
                    + dw.dot(property.angular_momentum_body_checked(dw)?)),
            )?);
        }
        let removed = finite(removed.value())?;
        if removed < 0.0 {
            return Err(AssemblyError::Invalid("negative relative kinetic energy"));
        }
        let balance = AssemblyBalance {
            linear_momentum_residual_world_ns: after.linear_world.sub(before.linear_world),
            angular_momentum_residual_world_nms: after.angular_world.sub(before.angular_world),
            energy_residual_j: finite(after.kinetic - before.kinetic + removed)?,
        };
        controls.admit(balance)?;
        poll(&mut cancelled)?;
        Ok(AssemblyLatch {
            assembly,
            state,
            removed_relative_energy_j: removed,
            balance,
        })
    }
}

// The rocket adapter has already admitted its dry/fuel partition and common-
// COM inertia. Sharing this owner keeps point impulses and spin/orbit accounting
// identical while retaining its specific inventory and inertia-law state.
#[derive(Clone, Copy, Default)]
pub(crate) struct ImpulseTotals([VectorSum; 2]);
pub(crate) struct ReleaseReceipt {
    pub supplied_work_j: f64,
    pub kinetic_before_j: f64,
    pub kinetic_after_j: f64,
    pub balance: AssemblyBalance,
}

pub(crate) fn release_components_into(
    parts: &[AssemblyPart],
    parent_properties: MassProperties,
    parent: RigidBodyState,
    impulses: &[SeparationImpulse],
    deltas: &mut [ImpulseTotals],
    bodies: &mut [RigidBodyState],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<ReleaseReceipt, AssemblyError> {
    if deltas.len() != parts.len() || bodies.len() != parts.len() {
        return Err(AssemblyError::Invalid(
            "matching native release workspace required",
        ));
    }
    for delta in deltas.iter_mut() {
        poll(cancelled)?;
        *delta = ImpulseTotals::default();
    }
    for impulse in impulses {
        poll(cancelled)?;
        if impulse.part_a >= parts.len()
            || impulse.part_b >= parts.len()
            || impulse.part_a == impulse.part_b
            || !impulse.point_assembly_m.is_finite()
            || !impulse.impulse_on_a_assembly_ns.is_finite()
            || !impulse.couple_on_a_assembly_nms.is_finite()
        {
            return Err(AssemblyError::Invalid(
                "finite impulse at one point on distinct existing parts required",
            ));
        }
        for (i, sign) in [(impulse.part_a, 1.0), (impulse.part_b, -1.0)] {
            let j = impulse.impulse_on_a_assembly_ns.scale(sign);
            let arm = impulse
                .point_assembly_m
                .sub(parts[i].pose_in_reference.position_world());
            let k = arm
                .cross(j)
                .add(impulse.couple_on_a_assembly_nms.scale(sign));
            deltas[i].0[0].add(vector(j)?);
            deltas[i].0[1].add(vector(k)?);
        }
    }
    let rotation = parent.pose().orientation();
    let velocity = parent.center_of_mass_velocity_world(parent_properties)?;
    let omega = parent_properties.angular_velocity_body_checked(parent.angular_momentum_body())?;
    let mut work = Sum::default();
    for ((part, delta), body) in parts.iter().zip(deltas).zip(bodies.iter_mut()) {
        poll(cancelled)?;
        let r = part.pose_in_reference.position_world();
        let axes = part.pose_in_reference.orientation();
        let j = vector(delta.0[0].value())?;
        let k = axes.rotate_world_to_body(vector(delta.0[1].value())?);
        let spin = part
            .properties
            .angular_momentum_body_checked(axes.rotate_world_to_body(omega))?
            .add(k);
        let relative = omega.cross(r).scale(part.properties.mass()).add(j);
        let linear = velocity
            .scale(part.properties.mass())
            .add(rotation.rotate_body_to_world(relative));
        *body = RigidBodyState::new(
            Pose::new(
                parent.pose().point_world_from_body(r)?,
                compose(rotation, axes)?,
            )?,
            linear,
            spin,
        )?;
        work.add(finite(
            0.5 * (j.dot(j.scale(part.properties.mass().recip()))
                + k.dot(part.properties.angular_velocity_body_checked(k)?)),
        )?);
    }
    let work = finite(work.value())?;
    release_accounting(parts, parent_properties, parent, bodies, work, cancelled)
}

// Also used after a domain adapter converts the native momenta to its own
// stored velocity fields, so receipts describe the actual rounded outputs.
pub(crate) fn release_accounting(
    parts: &[AssemblyPart],
    parent_properties: MassProperties,
    parent: RigidBodyState,
    bodies: &[RigidBodyState],
    work: f64,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<ReleaseReceipt, AssemblyError> {
    let velocity = parent.center_of_mass_velocity_world(parent_properties)?;
    let center = parent.pose().position_world();
    let before = totals(
        &[parent],
        [parent_properties].into_iter(),
        center,
        velocity,
        cancelled,
    )?;
    let after = totals(
        &bodies,
        parts.iter().map(|p| p.properties),
        center,
        velocity,
        cancelled,
    )?;
    let balance = AssemblyBalance {
        linear_momentum_residual_world_ns: after.linear_world.sub(before.linear_world),
        angular_momentum_residual_world_nms: after.angular_world.sub(before.angular_world),
        energy_residual_j: finite(after.kinetic - before.kinetic - work)?,
    };
    Ok(ReleaseReceipt {
        supplied_work_j: work,
        kinetic_before_j: before.kinetic,
        kinetic_after_j: after.kinetic,
        balance,
    })
}

impl AssemblyControls {
    fn validate(self, parts: usize, impulses: usize) -> Result<(), AssemblyError> {
        if parts == 0
            || self.maximum_parts == 0
            || [
                self.linear_momentum_tolerance_ns,
                self.angular_momentum_tolerance_nms,
                self.energy_tolerance_j,
            ]
            .into_iter()
            .any(|v| !v.is_finite() || v < 0.0)
        {
            return Err(AssemblyError::Invalid(
                "nonempty assembly, positive part cap, finite nonnegative balance gates required",
            ));
        }
        for (requested, maximum) in [
            (parts, self.maximum_parts),
            (impulses, self.maximum_impulses),
        ] {
            if requested > maximum {
                return Err(AssemblyError::Budget { requested, maximum });
            }
        }
        Ok(())
    }
    fn admit(self, balance: AssemblyBalance) -> Result<(), AssemblyError> {
        vector(balance.linear_momentum_residual_world_ns)?;
        vector(balance.angular_momentum_residual_world_nms)?;
        finite(balance.energy_residual_j)?;
        if maximum(balance.linear_momentum_residual_world_ns) > self.linear_momentum_tolerance_ns
            || maximum(balance.angular_momentum_residual_world_nms)
                > self.angular_momentum_tolerance_nms
            || balance.energy_residual_j.abs() > self.energy_tolerance_j
        {
            return Err(AssemblyError::Balance(balance));
        }
        Ok(())
    }
}

struct Totals {
    linear_world: Vec3,
    angular_world: Vec3,
    kinetic: f64,
}
fn totals(
    bodies: &[RigidBodyState],
    properties: impl Iterator<Item = MassProperties>,
    center: Vec3,
    boost: Vec3,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Totals, AssemblyError> {
    let mut linear = VectorSum::default();
    let mut angular = VectorSum::default();
    let mut first_moment = VectorSum::default();
    let mut kinetic = Sum::default();
    for (body, property) in bodies.iter().zip(properties) {
        poll(cancelled)?;
        let p = body.linear_momentum_world();
        let relative = vector(p.sub(boost.scale(property.mass())))?;
        let r = vector(body.pose().position_world().sub(center))?;
        let spin = body.angular_momentum_body();
        linear.add(p);
        angular.add(vector(
            body.pose()
                .orientation()
                .rotate_body_to_world(spin)
                .add(r.cross(relative)),
        )?);
        first_moment.add(vector(r.scale(property.mass()))?);
        kinetic.add(finite(
            0.5 * (relative.dot(relative.scale(property.mass().recip()))
                + spin.dot(property.angular_velocity_body_checked(spin)?)),
        )?);
    }
    Ok(Totals {
        linear_world: vector(linear.value())?,
        angular_world: vector(angular.value().add(first_moment.value().cross(boost)))?,
        kinetic: finite(kinetic.value())?,
    })
}

const ENTRIES: [(usize, usize); 6] = [(0, 0), (1, 1), (2, 2), (0, 1), (0, 2), (1, 2)];
fn compose(a: UnitQuaternion, b: UnitQuaternion) -> Result<UnitQuaternion, DynamicsError> {
    let q = a.multiply(b);
    UnitQuaternion::new(q.w, q.x, q.y, q.z)
}
fn array(v: Vec3) -> [f64; 3] {
    [v.x, v.y, v.z]
}
fn maximum(v: Vec3) -> f64 {
    v.x.abs().max(v.y.abs()).max(v.z.abs())
}
fn finite(v: f64) -> Result<f64, AssemblyError> {
    if v.is_finite() {
        Ok(v)
    } else {
        Err(AssemblyError::Invalid("unrepresentable assembly quantity"))
    }
}
fn vector(v: Vec3) -> Result<Vec3, AssemblyError> {
    if v.is_finite() {
        Ok(v)
    } else {
        Err(AssemblyError::Invalid("unrepresentable assembly vector"))
    }
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), AssemblyError> {
    if cancelled() {
        Err(AssemblyError::Cancelled)
    } else {
        Ok(())
    }
}

// Fixed input order with compensated sums: no N-by-N matrices or reductions
// over enormous common translational kinetic energy are needed.
#[derive(Clone, Copy, Debug, Default)]
struct Sum {
    value: f64,
    correction: f64,
}
impl Sum {
    fn add(&mut self, x: f64) {
        let next = self.value + x;
        self.correction += if self.value.abs() >= x.abs() {
            (self.value - next) + x
        } else {
            (x - next) + self.value
        };
        self.value = next;
    }
    fn value(self) -> f64 {
        self.value + self.correction
    }
}
#[derive(Clone, Copy, Debug, Default)]
struct VectorSum([Sum; 3]);
impl VectorSum {
    fn add(&mut self, v: Vec3) {
        for (sum, x) in self.0.iter_mut().zip(array(v)) {
            sum.add(x);
        }
    }
    fn value(self) -> Vec3 {
        Vec3::new(self.0[0].value(), self.0[1].value(), self.0[2].value())
    }
}
