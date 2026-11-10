//! Owner slice: the native fs-mbd rigid body on the composition kernel.
//!
//! `fs_mbd::component::RigidBody` is the first physics owner behind the
//! Component contract. It passes the conformance check, and the kernel joins
//! it to a spring defined in this test with no code written for the pair.
//! Run with `cargo test -p fs-system --features slices`.

// Exact float comparisons here assert exact zeros.
#![allow(clippy::float_cmp)]

use core::num::NonZeroUsize;

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, Cancelled, Cx, ExecMode, StreamKey};
use fs_mbd::component::RigidBody;
use fs_mbd::{
    Gravity, MassProperties, Pose, RigidBodyIntegrator, RigidBodyState, UnitQuaternion, Vec3,
    Wrench,
};
use fs_ports::component::conformance::{Case, Observed, check_component};
use fs_ports::component::{Class, Component, PortRole, PortSpec, Proposal, Receipt};
use fs_ports::{PortKind, StableId};
use fs_system::SystemBuilder;

fn with_cx<R>(gate: &CancelGate, f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let pool = ArenaPool::new(ArenaConfig::default());
    pool.scope(|arena| {
        let key = StreamKey {
            seed: 1,
            kernel_id: 2,
            tile: 3,
            iteration: 4,
        };
        let cx = Cx::new(gate, arena, key, Budget::INFINITE, ExecMode::Deterministic);
        f(&cx)
    })
}

const MASS_KG: f64 = 2.0;

fn properties() -> MassProperties {
    MassProperties::new(MASS_KG, Vec3::ZERO, Vec3::new(0.02, 0.05, 0.03)).expect("properties")
}

fn body(gravity: Gravity) -> RigidBody {
    RigidBody::new(RigidBodyIntegrator::new(gravity), properties()).expect("ports")
}

fn norm3(v: &[f64]) -> f64 {
    (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt()
}

fn state(position: Vec3, momentum: Vec3, spin: Vec3) -> RigidBodyState {
    let pose = Pose::new(position, UnitQuaternion::IDENTITY).expect("pose");
    RigidBodyState::new(pose, momentum, spin).expect("state")
}

#[test]
fn the_native_rigid_body_conforms_under_force_torque_and_gravity() {
    let gate = CancelGate::new_clock_free();
    let gravity = Gravity::new(Vec3::new(0.0, 0.0, -9.81)).expect("gravity");
    let component = body(gravity);
    // A tumbling asymmetric body: the rotational solve is nonlinear.
    let start = state(
        Vec3::new(0.3, -0.1, 1.2),
        Vec3::new(0.4, -0.2, 0.1),
        Vec3::new(0.010, -0.020, 0.015),
    );
    let stale = state(
        Vec3::new(9.0, 9.0, 9.0),
        Vec3::new(1.0, 2.0, 3.0),
        Vec3::ZERO,
    );
    let inputs = [1.5, -0.5, 2.0, 0.010, 0.020, -0.005];
    let case = Case {
        state: &start,
        stale: &stale,
        dt_s: 1.0e-3,
        inputs: &inputs,
    };
    let observed = with_cx(&gate, |cx| check_component(&component, &case, cx)).expect("conforms");
    let Observed::Stepped {
        receipt,
        audit,
        outputs,
    } = observed
    else {
        panic!("the rigid body takes the step");
    };
    // Independent of the component and of its tolerance: the integrator's own
    // account of the work of this force and torque. The kernel books the same
    // work from the two values on each wire.
    let wrench = Wrench {
        force_world: Vec3::new(inputs[0], inputs[1], inputs[2]),
        torque_body: Vec3::new(inputs[3], inputs[4], inputs[5]),
    };
    let native = RigidBodyIntegrator::new(gravity)
        .step(start, properties(), wrench, 1.0e-3)
        .expect("native step");
    let paired_j = 1.0e-3
        * (norm3(&inputs[..3]) * norm3(&outputs[..3]) + norm3(&inputs[3..]) * norm3(&outputs[3..]));
    assert!(native.work.applied_work_j.abs() > 0.1 * paired_j);
    assert!(
        (audit.port_energy_in_j - native.work.applied_work_j).abs() <= 1.0e-13 * paired_j,
        "{} vs {}",
        audit.port_energy_in_j,
        native.work.applied_work_j
    );
    assert_eq!(receipt.dissipated_j, 0.0);
    // What the midpoint step fails to conserve is reported, and it is small.
    let exchanged_j = audit.port_energy_in_j.abs() + receipt.stored_change_j.abs();
    assert!(receipt.numerical_loss_j.abs() <= 1.0e-9 * exchanged_j);
    // The force port returns the mean centre-of-mass velocity: the momentum
    // gains (force + weight) dt, and the mean is half way.
    let expected = 0.4 / MASS_KG + 0.5 * 1.0e-3 * 1.5 / MASS_KG;
    assert!(
        (outputs[0] - expected).abs() <= 1.0e-12,
        "{} vs {expected}",
        outputs[0]
    );
    let expected_z = 0.1 / MASS_KG + 0.5 * 1.0e-3 * (2.0 / MASS_KG - 9.81);
    assert!(
        (outputs[2] - expected_z).abs() <= 1.0e-12,
        "{} vs {expected_z}",
        outputs[2]
    );
}

/// Deterministic values in (-1, 1).
struct Sequence(u64);

impl Sequence {
    fn next(&mut self) -> f64 {
        self.0 ^= self.0 << 13;
        self.0 ^= self.0 >> 7;
        self.0 ^= self.0 << 17;
        (self.0 >> 11) as f64 / (1_u64 << 53) as f64 * 2.0 - 1.0
    }

    fn vector(&mut self, scale: f64) -> Vec3 {
        Vec3::new(
            scale * self.next(),
            scale * self.next(),
            scale * self.next(),
        )
    }
}

#[test]
fn the_rigid_body_conforms_across_heights_speeds_frames_and_loads() {
    // One operating point proves little. These states include the ones that
    // defeat a tolerance scaled by the energy exchanged in the step: free
    // fall far from the origin (kinetic energy and potential cancel), a
    // principal frame rotated against the body frame, and a torque nearly
    // perpendicular to the spin.
    let gate = CancelGate::new_clock_free();
    let gravity = Gravity::new(Vec3::new(0.0, 0.0, -9.81)).expect("gravity");
    let rotated = MassProperties::from_principal_frame(
        MASS_KG,
        Vec3::ZERO,
        Vec3::new(0.02, 0.05, 0.03),
        UnitQuaternion::new(0.9, 0.2, -0.3, 0.1).expect("orientation"),
    )
    .expect("rotated properties");
    let bodies = [
        body(gravity),
        RigidBody::new(RigidBodyIntegrator::new(gravity), rotated).expect("ports"),
    ];
    let stale = state(
        Vec3::new(9.0, 9.0, 9.0),
        Vec3::new(1.0, 2.0, 3.0),
        Vec3::ZERO,
    );
    let mut sequence = Sequence(0x9E37_79B9_7F4A_7C15);
    let mut refused = Vec::new();
    for case_index in 0..400_usize {
        let component = &bodies[case_index % 2];
        // Heights to 500 m and fall speeds to 100 m/s.
        let height = 500.0 * sequence.next();
        let fall = 100.0 * sequence.next();
        let start = state(
            Vec3::new(sequence.next(), sequence.next(), height),
            Vec3::new(sequence.next(), sequence.next(), MASS_KG * fall),
            sequence.vector(0.05),
        );
        let mut inputs = [0.0_f64; 6];
        match case_index % 4 {
            // Free fall: no load at all.
            0 => {}
            // A torque about one axis only, against a spin in the other two.
            1 => inputs[3] = 1.0e-3 * sequence.next(),
            // A general wrench.
            _ => {
                for value in &mut inputs {
                    *value = sequence.next();
                }
            }
        }
        let dt_s = if case_index % 3 == 0 { 1.0e-2 } else { 1.0e-3 };
        let case = Case {
            state: &start,
            stale: &stale,
            dt_s,
            inputs: &inputs,
        };
        if let Err(failure) = with_cx(&gate, |cx| check_component(component, &case, cx)) {
            refused.push(format!("case {case_index}: {failure}"));
        }
    }
    assert!(
        refused.is_empty(),
        "{} refused, first: {}",
        refused.len(),
        refused[0]
    );
}

#[test]
fn a_body_in_free_fall_follows_the_parabola_for_seconds() {
    // A constant force is integrated exactly by the midpoint rule, so the
    // trajectory is the parabola to roundoff. The body starts at rest at the
    // origin, where its stored energy is zero and stays the difference of two
    // growing terms.
    let gate = CancelGate::new_clock_free();
    let gravity = Gravity::new(Vec3::new(0.0, 0.0, -9.81)).expect("gravity");
    let mut builder = SystemBuilder::new();
    let rigid = builder
        .add(
            "body",
            body(gravity),
            state(Vec3::ZERO, Vec3::ZERO, Vec3::new(0.010, -0.020, 0.015)),
        )
        .expect("add body");
    for name in ["force", "torque"] {
        builder
            .boundary(
                builder.port(rigid, name).expect("port"),
                Box::new(|_, _, value| value.fill(0.0)),
            )
            .expect("boundary");
    }
    let mut system = builder.seal().expect("seal");
    let (dt, steps) = (1.0e-3, 5_000_u32);
    with_cx(&gate, |cx| {
        for step in 0..steps {
            system
                .step(dt, cx)
                .unwrap_or_else(|error| panic!("step {step}: {error}"));
        }
    });
    let time = f64::from(steps) * dt;
    let end = *system.state::<RigidBodyState>(rigid).expect("state");
    let (fallen, speed) = (0.5 * 9.81 * time * time, 9.81 * time);
    let z = end.pose().position_world().z;
    let p = end.linear_momentum_world().z;
    assert!((z + fallen).abs() <= 1.0e-10 * fallen, "{z} vs -{fallen}");
    assert!((p + MASS_KG * speed).abs() <= 1.0e-10 * MASS_KG * speed);
    // Nothing did work on the body: kinetic energy grew by what the
    // potential lost, to a part in 1e10 of the kinetic energy.
    let kinetic_j = 0.5 * MASS_KG * speed * speed;
    let ledger = system.ledger();
    assert_eq!(ledger.boundary_energy_in_j(), 0.0);
    assert!(ledger.stored_change_j().abs() <= 1.0e-10 * kinetic_j);
    assert!(ledger.numerical_loss_j().abs() <= 1.0e-10 * kinetic_j);
}

/// Isotropic three-axis spring to the origin: receives the flow of its end,
/// returns the effort of the mean displacement.
struct Spring3 {
    stiffness_n_m: f64,
    ports: [PortSpec; 1],
}

impl Spring3 {
    fn new(stiffness_n_m: f64) -> Self {
        let three = NonZeroUsize::MIN.saturating_add(2);
        let port = PortSpec::power(
            StableId::new("end").expect("id"),
            PortKind::MechanicalForceVelocity,
            PortRole::FlowIn,
            three,
        )
        .expect("port");
        Self {
            stiffness_n_m,
            ports: [port],
        }
    }

    fn energy(&self, q: &[f64; 3]) -> f64 {
        0.5 * self.stiffness_n_m * q.iter().map(|c| c * c).sum::<f64>()
    }
}

impl Component for Spring3 {
    type State = [f64; 3];
    type Error = Cancelled;

    fn ports(&self) -> &[PortSpec] {
        &self.ports
    }

    fn class(&self) -> Class {
        Class::ExactLinear
    }

    fn stored_energy_j(&self, q: &[f64; 3]) -> Result<f64, Cancelled> {
        Ok(self.energy(q))
    }

    fn propose(
        &self,
        old: &[f64; 3],
        dt_s: f64,
        inputs: &[f64],
        new: &mut [f64; 3],
        outputs: &mut [f64],
        cx: &Cx<'_>,
    ) -> Result<Proposal, Cancelled> {
        cx.checkpoint()?;
        let mut stored_change_j = 0.0;
        for axis in 0..3 {
            new[axis] = old[axis] + dt_s * inputs[axis];
            outputs[axis] = 0.5 * self.stiffness_n_m * (old[axis] + new[axis]);
            stored_change_j += outputs[axis] * (new[axis] - old[axis]);
        }
        let scale = self.energy(old) + self.energy(new) + stored_change_j.abs();
        Ok(Proposal::Step(Receipt::lossless(
            stored_change_j,
            64.0 * f64::EPSILON * scale,
        )))
    }
}

#[test]
fn a_tumbling_rigid_body_on_a_spring_is_composed_by_the_kernel() {
    let gate = CancelGate::new_clock_free();
    let stiffness = 800.0;
    let (x0, p0) = (Vec3::new(0.02, -0.01, 0.015), Vec3::new(0.4, -0.2, 0.1));
    let spin = Vec3::new(0.010, -0.020, 0.015);

    let mut builder = SystemBuilder::new();
    let rigid = builder
        .add("body", body(Gravity::ZERO), state(x0, p0, spin))
        .expect("add body");
    // The spring pulls the body towards the origin: its displacement is minus
    // the body position, and its flow is minus the body velocity.
    let spring = builder
        .add("spring", Spring3::new(stiffness), [-x0.x, -x0.y, -x0.z])
        .expect("add spring");
    let force = builder.port(rigid, "force").expect("force port");
    builder
        .bond(force, builder.port(spring, "end").expect("end port"))
        .expect("bond");
    builder
        .boundary(
            builder.port(rigid, "torque").expect("torque port"),
            Box::new(|_, _, torque| torque.fill(0.0)),
        )
        .expect("boundary");
    let mut system = builder.seal().expect("seal");
    assert_eq!(system.plan_shape(), (1, 1));

    let (dt, steps) = (5.0e-4, 4_000_u32);
    let iterations = with_cx(&gate, |cx| {
        let mut iterations = 0;
        for _ in 0..steps {
            iterations += system.step(dt, cx).expect("step").iterations;
        }
        iterations
    });
    // The body offers no tangent, so the kernel differences its port map; the
    // translational loop is linear and still settles in a few iterations.
    assert!(
        iterations <= 4 * steps,
        "{iterations} iterations for {steps} steps"
    );

    // Translation: the monolithic midpoint rule for m x'' = -k x, per axis.
    let half = 0.5 * dt / MASS_KG;
    let advance = |mut x: f64, mut v: f64| {
        for _ in 0..steps {
            let mean = (v - half * stiffness * x) / (1.0 + 0.5 * dt * half * stiffness);
            x += dt * mean;
            v = 2.0 * mean - v;
        }
        (x, v)
    };
    let end = *system.state::<RigidBodyState>(rigid).expect("body state");
    let position = end.pose().position_world();
    let momentum = end.linear_momentum_world();
    for (got_x, got_p, x, p) in [
        (position.x, momentum.x, x0.x, p0.x),
        (position.y, momentum.y, x0.y, p0.y),
        (position.z, momentum.z, x0.z, p0.z),
    ] {
        let (want_x, want_v) = advance(x, p / MASS_KG);
        assert!(
            (got_x - want_x).abs() <= 1.0e-8 * 0.03,
            "{got_x} vs {want_x}"
        );
        assert!(
            (got_p - MASS_KG * want_v).abs() <= 1.0e-8 * 0.5,
            "{got_p} vs {want_v}"
        );
    }

    // Rotation is torque-free: the size of the body angular momentum stays.
    let size = |v: Vec3| (v.x * v.x + v.y * v.y + v.z * v.z).sqrt();
    let spin_after = end.angular_momentum_body();
    assert!((size(spin_after) - size(spin)).abs() <= 1.0e-10 * size(spin));
    // The body did tumble: an asymmetric top does not keep its spin axis.
    assert!((spin_after.x - spin.x).abs() > 1.0e-4);

    // The ledger: the zero torque did no work, and what the two components
    // store changed only by the body's reported numerical loss.
    let ledger = system.ledger();
    let start_j = 0.5 * (p0.x * p0.x + p0.y * p0.y + p0.z * p0.z) / MASS_KG
        + 0.5 * stiffness * (x0.x * x0.x + x0.y * x0.y + x0.z * x0.z);
    assert_eq!(ledger.boundary_energy_in_j(), 0.0);
    assert!(ledger.energy_residual_j().abs() <= 1.0e-10 * start_j);
    assert!(ledger.connection_defect_j().abs() <= 1.0e-9 * start_j);
    assert!((ledger.stored_change_j() + ledger.numerical_loss_j()).abs() <= 1.0e-9 * start_j);
    assert!(ledger.numerical_loss_j().abs() <= 1.0e-8 * start_j);
    assert_eq!(ledger.nodes[rigid.index()].dissipated_j.value(), 0.0);
}

#[test]
fn a_body_far_out_on_a_level_surface_of_an_oblique_field_keeps_its_account() {
    // Second independent review, 2026-11-10: with the tolerance sized by the
    // magnitude of the potential's total, 278 of these 2,000 steps failed
    // the audit, because the three terms of `-m g . x` are large and cancel.
    let gate = CancelGate::new_clock_free();
    let gravity = Gravity::new(Vec3::new(6.9367, 0.0, -6.9367)).expect("gravity");
    let component = body(gravity);
    let mut current = state(
        Vec3::new(1000.0, 0.0, 1000.0),
        Vec3::new(0.4, -0.2, 0.1),
        Vec3::new(0.010, -0.020, 0.015),
    );
    let stale = current;
    let inputs = [0.0; 6];
    let mut outputs = [0.0; 6];
    with_cx(&gate, |cx| {
        for step in 0..2_000 {
            let case = Case {
                state: &current,
                stale: &stale,
                dt_s: 1.0e-3,
                inputs: &inputs,
            };
            check_component(&component, &case, cx)
                .unwrap_or_else(|failure| panic!("step {step}: {failure}"));
            let mut next = current;
            component
                .propose(&current, 1.0e-3, &inputs, &mut next, &mut outputs, cx)
                .expect("step");
            current = next;
        }
    });
}
