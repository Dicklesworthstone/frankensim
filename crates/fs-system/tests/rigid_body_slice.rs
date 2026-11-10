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

fn body(gravity: Gravity) -> RigidBody {
    let properties =
        MassProperties::new(MASS_KG, Vec3::ZERO, Vec3::new(0.02, 0.05, 0.03)).expect("properties");
    RigidBody::new(RigidBodyIntegrator::new(gravity), properties).expect("ports")
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
    // The wire values reproduce the integrator's own applied work: the first
    // law closes inside the owner's tolerance, with nothing booked as heat.
    assert!(audit.balance_defect_j.abs() <= receipt.energy_tolerance_j);
    assert_eq!(receipt.dissipated_j, 0.0);
    // What the midpoint step fails to conserve is reported, and it is small.
    let exchanged_j = audit.port_energy_in_j.abs() + receipt.stored_change_j.abs();
    assert!(receipt.numerical_loss_j.abs() <= 1.0e-9 * exchanged_j);
    // The force port returns the mean centre-of-mass velocity.
    let expected = 0.4 / MASS_KG + 0.5 * 1.0e-3 * 1.5 / MASS_KG;
    assert!(
        (outputs[0] - expected).abs() <= 1.0e-12,
        "{} vs {expected}",
        outputs[0]
    );
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
