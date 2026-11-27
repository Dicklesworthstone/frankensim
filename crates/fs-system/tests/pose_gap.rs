//! `PoseGap` (`frankensim-rc-1109-root-3gc12.13.1`, O1): the gap from a point on a moving body to a
//! fixed face, as the body's pose changes, under a contact.
//!
//! A rotor of inertia `I` on a torsional spring `k` carries an arm of length `r`, whose tip meets
//! a fixed wall at the turn `theta_c`. Released at rest from the turn `theta0`, it swings as
//! `theta0 cos(omega t)`, `omega = sqrt(k / I)`, until the tip reaches the wall at
//! `t_c = acos(theta_c / theta0) / omega`. The tip bounces off a lossless Hertz contact and the
//! swing returns to `theta0`. The contact's closure is the geometry's gap at every step, to a
//! rounding bound that grows with the step count; so is the energy with the contact's potential
//! taken at the geometry's penetration.

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_dcontact::Obstacle;
use fs_dcontact::component::{ContactInterface, ContactSolve, ContactState};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_ports::PortKind;
use fs_ports::component::conformance::{Case, check_component, check_port_tangent};
use fs_ports::component::{Component, LoadLine};
use fs_system::elements::{Mass, Motion, Pose, PoseGap, Spring};
use fs_system::{AdmissionError, NodeId, SystemBuilder};

const INERTIA_KG_M2: f64 = 1.0e-3;
const STIFFNESS_N_M_RAD: f64 = 1.0;
const ARM_M: f64 = 0.1;
/// The turn the rotor is released from [rad].
const THETA0: f64 = 0.3;
/// The turn at which the arm's tip reaches the wall [rad].
const THETA_C: f64 = 0.15;
const DT_S: f64 = 1.0e-5;
/// Past the bounce and the swing's return to its turning point.
const STEPS: usize = 8000;
/// Hertz stiffness [N/m^1.5]: `4/3 E* sqrt(R)` with `E* = 0.8 GPa` and `R = 3 mm`.
const HERTZ: f64 = 4.0 / 3.0 * 0.8e9 * 0.054_772_255_750_516_61;
const SOLVE: ContactSolve = ContactSolve {
    max_iterations: 100,
    maximum_force_n: 1.0e6,
    maximum_penetration_m: 1.0e-3,
    force_absolute_tolerance_n: 1.0e-9,
    force_relative_tolerance: 1.0e-12,
};

fn with_cx<T>(f: impl FnOnce(&Cx<'_>) -> T) -> T {
    let gate = CancelGate::new_clock_free();
    let pool = ArenaPool::new(ArenaConfig::default());
    pool.scope(|arena| {
        let key = StreamKey {
            seed: 7,
            kernel_id: 8,
            tile: 9,
            iteration: 10,
        };
        f(&Cx::new(
            &gate,
            arena,
            key,
            Budget::INFINITE,
            ExecMode::Deterministic,
        ))
    })
}

/// The tip at `r (cos theta, sin theta)` about a pivot held at the origin, and the wall's face on
/// the line `x = r cos(theta_c)`, solid beyond it: the gap is `r (cos(theta_c) - cos(theta))`.
fn reach() -> PoseGap {
    let face = [ARM_M * THETA_C.cos(), 0.0];
    PoseGap::new([0.0, 0.0], [0.0, 1.0], [ARM_M, 0.0], face, [-1.0, 0.0]).unwrap()
}

fn wall() -> ContactInterface {
    let law = Obstacle::new(
        vec![-1.0],
        1,
        1,
        vec![0.0],
        vec![1.0],
        HERTZ,
        1.5,
        "wall".into(),
    );
    ContactInterface::new(law.unwrap(), SOLVE).unwrap()
}

/// A step's record. The energies are changes since the start: the kernel's with the contact's
/// own store, and the geometry's with the contact's potential at the penetration `-gap`.
struct Step {
    theta: f64,
    closure_m: f64,
    gap_m: f64,
    reaction_n: f64,
    kernel_j: f64,
    geometric_j: f64,
}

fn swing(count: usize) -> Vec<Step> {
    let rotational = PortKind::RotationalTorqueAngularVelocity;
    let at_rest = Motion {
        position_m: 0.0,
        velocity_m_s: 0.0,
    };
    // The spring's end and the gap's turn both receive minus the rotor's rate, so the turn is
    // the spring's extension: theta0 at the start.
    let start = Pose {
        s_m: 0.0,
        theta_rad: THETA0,
    };
    let shape = reach();
    let gap0 = shape.gap_m(&start);
    let mut b = SystemBuilder::new();
    let inertia = Mass::in_domain(rotational, INERTIA_KG_M2).unwrap();
    let rotor = b.add_with_tangent("rotor", inertia, at_rest).unwrap();
    let torsion = Spring::in_domain(rotational, STIFFNESS_N_M_RAD, 0.0).unwrap();
    let spring = b.add_with_tangent("spring", torsion, THETA0).unwrap();
    let arm = b.add_with_tangent("arm", shape.clone(), start).unwrap();
    let contact = b
        .add_on_line("wall", wall(), ContactState::at(-gap0))
        .unwrap();
    let p = |b: &SystemBuilder, node, name| b.port(node, name).unwrap();
    let rest = || Box::new(|_: f64, _: f64, v: &mut [f64]| v[0] = 0.0);
    let members = [p(&b, spring, "a"), p(&b, arm, "theta")];
    b.common_flow(p(&b, rotor, "p"), &members).unwrap();
    b.boundary(p(&b, spring, "b"), rest()).unwrap();
    b.boundary(p(&b, arm, "s"), rest()).unwrap();
    b.bond(p(&b, contact, "a"), p(&b, arm, "gap")).unwrap();
    b.boundary(p(&b, contact, "b"), rest()).unwrap();
    let mut system = b.seal().unwrap();
    let probe = wall();
    let potential = |gap: f64| probe.stored_energy_j(&ContactState::at(-gap)).unwrap();
    let base_j = potential(gap0);
    with_cx(|cx| {
        (0..count)
            .map(|_| {
                system.step(DT_S, cx).unwrap();
                let pose = *system.state::<Pose>(arm).unwrap();
                let touch = system.state::<ContactState>(contact).unwrap();
                let ledger = system.ledger();
                let stored = |node: NodeId| ledger.nodes[node.index()].stored_change_j.value();
                let gap_m = shape.gap_m(&pose);
                let moving_j = stored(rotor) + stored(spring);
                Step {
                    theta: pose.theta_rad,
                    closure_m: touch.closure_m,
                    gap_m,
                    reaction_n: touch.reaction_n,
                    kernel_j: moving_j + stored(contact),
                    geometric_j: moving_j + potential(gap_m) - base_j,
                }
            })
            .collect()
    })
}

#[test]
fn the_tip_meets_the_wall_on_time_on_the_geometry_and_the_swing_returns() {
    let steps = swing(STEPS);
    let omega = (STIFFNESS_N_M_RAD / INERTIA_KG_M2).sqrt();
    let t_c = (THETA_C / THETA0).acos() / omega;
    let first = 1 + steps
        .iter()
        .position(|s| s.closure_m > 0.0)
        .expect("the tip reaches the wall");
    assert_eq!(
        first,
        (t_c / DT_S).ceil() as usize,
        "the first step the contact closes"
    );
    // The rounding each step adds to the closure and to the gap evaluated at the pose, both of
    // the size of the gap's terms: |g0| + |c| + |k|, with the pivot held.
    let scale_m = 2.0 * ARM_M;
    let total_j = 0.5 * STIFFNESS_N_M_RAD * THETA0 * THETA0;
    let (mut worst_closure, mut worst_kernel, mut worst_geometric) = (0.0_f64, 0.0_f64, 0.0_f64);
    let mut peak_n = 0.0_f64;
    for (k, step) in (1_u32..).zip(&steps) {
        let closure_allowed = 8.0 * f64::EPSILON * scale_m * f64::from(k + 1);
        let off = (step.closure_m + step.gap_m).abs() / closure_allowed;
        assert!(
            off <= 1.0,
            "step {k}: closure {} m, gap {} m",
            step.closure_m,
            step.gap_m
        );
        peak_n = peak_n.max(step.reaction_n);
        let kernel_allowed = 64.0 * f64::EPSILON * total_j * f64::from(k);
        let kernel = step.kernel_j.abs() / kernel_allowed;
        let geometric = step.geometric_j.abs() / (kernel_allowed + peak_n * closure_allowed);
        assert!(
            kernel <= 1.0,
            "step {k}: the kernel's energy moved {} J",
            step.kernel_j
        );
        assert!(
            geometric <= 1.0,
            "step {k}: the energy moved {} J",
            step.geometric_j
        );
        worst_closure = worst_closure.max(off);
        worst_kernel = worst_kernel.max(kernel);
        worst_geometric = worst_geometric.max(geometric);
    }
    let last = steps.iter().rposition(|s| s.reaction_n > 0.0).unwrap();
    let (turn, peak) = steps[last..]
        .iter()
        .enumerate()
        .fold((0, f64::MIN), |best, (i, s)| {
            if s.theta > best.1 { (i, s.theta) } else { best }
        });
    assert!(
        last + turn + 1 < steps.len(),
        "the swing turns inside the window"
    );
    // At the step nearest the turning point the rotor still moves by up to omega theta0 omega
    // dt / 2: the deficit of the sampled peak is at most omega^2 theta0 dt^2 / 8.
    let sampled = omega * omega * THETA0 * DT_S * DT_S / 8.0;
    let deficit = THETA0 - peak;
    assert!(
        (-1.0e-12..=sampled).contains(&deficit),
        "the swing returns to {peak} rad, {deficit} short of {THETA0}"
    );
    println!(
        "contact from step {first} (t_c / dt = {:.3}) to step {}, peak {peak_n:.3} N; worst closure \
         {worst_closure:.3}, kernel energy {worst_kernel:.3}, geometric energy {worst_geometric:.3} \
         of their bounds; return deficit {deficit:.3e} rad of {sampled:.3e}",
        t_c / DT_S,
        last + 1,
    );
}

/// The contact reads mean flows only, so the kernel joins it to the gap, whose end it does not
/// tell. A contact that reads an end input at `a` needs a neighbour that tells it, and is refused.
#[test]
fn a_contact_that_reads_its_end_is_refused_beside_the_gap() {
    let start = Pose {
        s_m: 0.0,
        theta_rad: THETA0,
    };
    let shape = reach();
    let closure = ContactState::at(-shape.gap_m(&start));
    let reading = wall().reading_ends([true, false]).unwrap();
    assert!(wall().reads_mean_flows_only() && !reading.reads_mean_flows_only());
    let mut b = SystemBuilder::new();
    let arm = b.add_with_tangent("arm", shape, start).unwrap();
    let contact = b.add_on_line("wall", reading, closure).unwrap();
    let p = |b: &SystemBuilder, node, name| b.port(node, name).unwrap();
    let rest = || Box::new(|_: f64, _: f64, v: &mut [f64]| v[0] = 0.0);
    b.boundary(p(&b, arm, "s"), rest()).unwrap();
    b.boundary(p(&b, arm, "theta"), rest()).unwrap();
    b.bond(p(&b, contact, "a"), p(&b, arm, "gap")).unwrap();
    b.boundary(p(&b, contact, "b"), rest()).unwrap();
    match b.seal() {
        Err(AdmissionError::LineEnd { port, reason }) => {
            assert_eq!(port, "wall.a");
            assert!(reason.contains("does not tell its end"), "{reason}");
        }
        Err(other) => panic!("{other}"),
        Ok(_) => panic!("sealed: an end input beside a neighbour that does not tell its end"),
    }
}

#[test]
fn pose_gap_geometry_tangents_and_refusals() {
    let (pivot, axis, point, face, normal) = (
        [0.02, -0.01],
        [0.6, 0.8],
        [0.07, 0.03],
        [0.1, 0.05],
        [-0.8, 0.6],
    );
    let shape = PoseGap::new(pivot, axis, point, face, normal).unwrap();
    // The gap from the geometry itself: P = p0 + s u + R(theta) d, g = n . (P - A).
    let geometry = |pose: &Pose| {
        let (sin, cos) = pose.theta_rad.sin_cos();
        let x = pivot[0] + pose.s_m * axis[0] + cos * point[0] - sin * point[1];
        let y = pivot[1] + pose.s_m * axis[1] + sin * point[0] + cos * point[1];
        normal[0] * (x - face[0]) + normal[1] * (y - face[1])
    };
    let old = Pose {
        s_m: 0.004,
        theta_rad: 0.4,
    };
    let stale = Pose {
        s_m: 1.0,
        theta_rad: -2.0,
    };
    assert!((shape.gap_m(&old) - geometry(&old)).abs() <= 8.0 * f64::EPSILON * 0.2);
    let dt_s = 1.0e-3;
    with_cx(|cx| {
        for inputs in [[0.3, -2.0, 15.0], [0.3, 0.0, 15.0], [0.0, 1.0e-9, 2.0]] {
            let case = Case {
                state: &old,
                stale: &stale,
                dt_s,
                inputs: &inputs,
            };
            check_component(&shape, &case, cx).unwrap();
            check_port_tangent(&shape, &case, 1.0e-6, 1.0e-12, cx).unwrap();
            let (mut new, mut out) = (stale, [0.0; 3]);
            shape
                .propose(&old, dt_s, &inputs, &mut new, &mut out, cx)
                .unwrap();
            // The gap's rate over the step is its change; the efforts are -R J, so the ports'
            // powers cancel.
            let change = geometry(&new) - geometry(&old);
            assert!(
                (dt_s * out[2] - change).abs() <= 16.0 * f64::EPSILON * 0.2,
                "{inputs:?}"
            );
            let power: f64 = out.iter().zip(&inputs).map(|(e, f)| e * f).sum();
            assert!(power.abs() <= 8.0 * f64::EPSILON * 15.0 * 2.0, "{power}");
            if inputs[1] == 0.0 {
                // A pose that turns by nothing: the slope at the pose, d g / d theta.
                let (sin, cos) = old.theta_rad.sin_cos();
                let tip = [
                    -sin * point[0] - cos * point[1],
                    cos * point[0] - sin * point[1],
                ];
                let slope = normal[0] * tip[0] + normal[1] * tip[1];
                assert!((out[1] + inputs[2] * slope).abs() <= 8.0 * f64::EPSILON * 15.0 * 0.1);
            }
        }
    });
    let mut words = Vec::new();
    assert!(shape.encode_state(&old, &mut words));
    assert_eq!(shape.decode_state(&words), Some(old));
    let bad = Pose {
        s_m: f64::NAN,
        theta_rad: 0.0,
    };
    assert!(!shape.encode_state(&bad, &mut Vec::new()));
    assert!(PoseGap::new(pivot, [0.0, 0.0], point, face, normal).is_err());
    assert!(PoseGap::new(pivot, axis, point, face, [0.0, 0.0]).is_err());
    assert!(PoseGap::new([f64::NAN, 0.0], axis, point, face, normal).is_err());
}
