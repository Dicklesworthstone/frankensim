//! Joint normal/traction consistency and centered energy checks.

use super::contact_island_tests::{dense_solve, loaded_body};
use super::friction_contact_tests::{point_with_prediction, turning_energy, turning_oscillator};
use super::*;

fn solve_joint(
    bodies: &mut [Body],
    contacts: &mut [Contact],
    frictions: &mut [Friction],
    k: f64,
) -> FrictionSolveReport {
    friction_response::solve_contacts_with_friction_checked(
        bodies,
        contacts,
        frictions,
        &mut ContactQueue::default(),
        &mut contact_response::Workspace::default(),
        &mut friction_response::Workspace::default(),
        k,
        ContactSolveControls::default(),
        FrictionSolveControls::default(),
    )
    .unwrap()
    .1
}

fn centered_slip(bodies: &[Body], terms: &[(usize, usize, f64)]) -> f64 {
    terms
        .iter()
        .map(|&(b, i, w)| w * (bodies[b].wn[i] - bodies[b].wp[i]))
        .sum()
}

fn normal_pose(bodies: &[Body], c: &Contact, previous: bool) -> f64 {
    let value = |b: usize, i: usize| {
        if previous {
            bodies[b].wp[i]
        } else {
            bodies[b].wn[i]
        }
    };
    c.a.iter().map(|&(b, i, w)| w * value(b, i)).sum::<f64>()
        - c.b.iter().map(|&(b, i, w)| w * value(b, i)).sum::<f64>()
        - c.gap
}

fn check_normal_laws(bodies: &[Body], contacts: &[Contact], k: f64) {
    for c in contacts.iter().filter(|c| c.active) {
        let (d, previous) = (normal_pose(bodies, c, false), normal_pose(bodies, c, true));
        let Law::Power {
            k: stiffness,
            alpha,
            mu,
        } = c.law;
        let expected = match &c.crush {
            Some(crush) => crush_phi(d, previous, stiffness, crush, mu, k).0,
            None => match c.retardation_s {
                Some(retardation) => viscoelastic_phi(d, previous, stiffness, retardation, k).0,
                None => contact_phi(d, previous, stiffness, alpha, mu, k).0,
            },
        };
        assert!((c.force - expected).abs() < 3e-9 * expected.abs().max(1.0));
        assert!((c.d_now - d).abs() < 2e-13);
    }
}

#[test]
fn joint_contact_recovers_analytic_sticking_and_sliding_after_normal_feedback() {
    let k = 0.001;
    for mu in [0.12, 0.2] {
        let (mut bodies, mut contacts, mut frictions) = turning_oscillator(k);
        frictions[0].sx.insert(0, (0, 0, 1.0));
        frictions[0].mu = mu;
        let initial = turning_energy(&bodies, k, false);
        let predicted_tangent = bodies[1].wn[0];
        let report = solve_joint(&mut bodies, &mut contacts, &mut frictions, k);
        let (normal, tangent, stuck) = if mu < 0.125 {
            let n = 1000.0 / (1.5 + 0.5 * mu);
            (n, -mu * n, false)
        } else {
            (640.0, -80.0, true)
        };
        // Eliminate the three scalar mass equations and the linear potential
        // secant independently of the joint nonlinear solver.
        let applied_tangent = (bodies[1].wn[0] - predicted_tangent) / (k * k);
        assert!((contacts[0].force - normal).abs() < 2e-7);
        assert!((applied_tangent - tangent).abs() < 2e-7);
        assert!((frictions[0].force - tangent.abs()).abs() < 2e-7);
        assert_eq!(frictions[0].stuck, stuck);
        check_normal_laws(&bodies, &contacts, k);
        let slip = centered_slip(&bodies, &frictions[0].sx);
        if stuck {
            assert!(slip.abs() < 2e-13);
        } else {
            assert!(slip > 5e-6);
        }
        let heat = -0.5 * applied_tangent * slip;
        assert!(heat >= -1e-12);
        assert!((report.dissipated_j - heat).abs() < 3e-11);
        assert!((turning_energy(&bodies, k, true) + heat - initial).abs() < 3e-10);
        assert!(report.normal_residual <= 1.0 && report.friction_residual <= 1.0);
    }
}

#[test]
fn joint_hertz_friction_closes_sustained_resonance_storage_and_heat() {
    let (k, preload, equilibrium, spring) = (1e-4, 1000.0, 0.001_f64, 5e5);
    let stiffness = preload / equilibrium.powf(1.5);
    let mut bodies = vec![
        point_with_prediction("preloaded Hertz", equilibrium, equilibrium, preload, k),
        point_with_prediction("tangent resonance", 0.0003, 0.00027, -spring * 0.0003, k),
        point_with_prediction("transverse tangent", 0.0, 0.0, 0.0, k),
    ];
    let mut contacts = vec![Contact::new(
        vec![(0, 0, 1.0)],
        vec![],
        0.0,
        Law::Power {
            k: stiffness,
            alpha: 1.5,
            mu: 0.0,
        },
    )];
    let mut frictions = vec![Friction {
        contact: 0,
        sx: vec![(0, 0, 1.0), (1, 0, 1.0)],
        sy: vec![(2, 0, 1.0)],
        mu: 0.18,
        force: 0.0,
        stuck: false,
    }];
    let potential = |q: f64| stiffness * q.max(0.0).powf(2.5) / 2.5 - preload * q;
    let energy = |bs: &[Body], after: bool| {
        let a = |i: usize| if after { bs[i].wn[0] } else { bs[i].w[0] };
        let b = |i: usize| if after { bs[i].w[0] } else { bs[i].wp[0] };
        let kinetic: f64 = (0..3).map(|i| 0.5 * ((a(i) - b(i)) / k).powi(2)).sum();
        kinetic + 0.5 * (potential(a(0)) + potential(b(0))) - potential(equilibrium)
            + 0.5 * spring * a(1) * b(1)
    };
    let initial = energy(&bodies, false);
    let (mut queue, mut normal_work, mut friction_work) = (
        ContactQueue::default(),
        contact_response::Workspace::default(),
        friction_response::Workspace::default(),
    );
    let (mut heat, mut worst, mut sliding) = (0.0, 0.0_f64, 0);
    for _ in 0..600 {
        bodies[0].f[0] = preload;
        bodies[1].f[0] = -spring * bodies[1].w[0];
        for body in &mut bodies {
            body.predict(k);
        }
        let predicted = [bodies[1].wn[0], bodies[2].wn[0]];
        let report = friction_response::solve_contacts_with_friction_checked(
            &mut bodies,
            &mut contacts,
            &mut frictions,
            &mut queue,
            &mut normal_work,
            &mut friction_work,
            k,
            ContactSolveControls::default(),
            FrictionSolveControls::default(),
        )
        .unwrap()
        .1;
        check_normal_laws(&bodies, &contacts, k);
        let traction = [
            (bodies[1].wn[0] - predicted[0]) / (k * k),
            (bodies[2].wn[0] - predicted[1]) / (k * k),
        ];
        let slip = [
            centered_slip(&bodies, &frictions[0].sx),
            centered_slip(&bodies, &frictions[0].sy),
        ];
        let measured_heat = -0.5 * (traction[0] * slip[0] + traction[1] * slip[1]);
        assert!(measured_heat >= -1e-12);
        assert!((report.dissipated_j - measured_heat).abs() < 2e-11);
        assert!(traction[0].hypot(traction[1]) <= frictions[0].mu * contacts[0].force + 2e-6);
        heat += measured_heat;
        worst = worst.max(((energy(&bodies, true) + heat - initial) / initial).abs());
        sliding += usize::from(!frictions[0].stuck && frictions[0].force > 0.0);
        for body in &mut bodies {
            body.rotate();
        }
    }
    assert!(sliding > 0 && heat > 0.0);
    assert!(
        worst < 3e-8,
        "joint Hertz/friction energy closure {worst:e}"
    );
    eprintln!(
        "joint Hertz/friction: 600-step relative energy closure {worst:.12e}, heat {heat:.12e} J"
    );
}

#[test]
fn joint_loaded_response_recovers_a_prescribed_dense_nodal_solution() {
    let k = 0.002;
    let mut body = loaded_body(k);
    let load = body.head_radiation.as_ref().unwrap();
    let mut operator = vec![0.0; 16];
    for i in 0..4 {
        for j in 0..4 {
            operator[i * 4 + j] = body.air_mass * body.h.powi(4) / (body.area.powi(2) * k * k)
                + body.g[i] * body.g[j] / 4.0;
            if i == j {
                operator[i * 4 + j] += 1.0 / body.cf[i][0];
            }
            for r in 0..load.cols.len() {
                for s in 0..load.cols.len() {
                    operator[i * 4 + j] += load.cols[r][i]
                        * load.cols[s][j]
                        * (load.mass[r * load.cols.len() + s] / (k * k)
                            + load.damping[r * load.cols.len() + s] / (2.0 * k));
                }
            }
        }
    }
    let n = [1.0, -0.2, 0.3, 0.0];
    let tx = [0.2, 1.0, -0.3, 0.4];
    let ty = [-0.1, 0.2, 0.5, 1.0];
    let final_pose: [f64; 4] = [0.0008, 0.0002, 0.0004, -0.0001];
    let previous: [f64; 4] = [0.0007, 0.00015, 0.0003, -0.00012];
    let dot = |a: &[f64; 4], b: &[f64; 4]| a.iter().zip(b).map(|(a, b)| a * b).sum::<f64>();
    let (d, dp) = (dot(&n, &final_pose), dot(&n, &previous));
    let stiffness = 2e7;
    let normal = stiffness * (d.powf(2.5) - dp.powf(2.5)) / (2.5 * (d - dp));
    let delta = std::array::from_fn(|i| final_pose[i] - previous[i]);
    let slip = [dot(&tx, &delta), dot(&ty, &delta)];
    let cap = 0.25 * normal;
    let tangent = [
        cap * slip[0] / slip[0].hypot(slip[1]),
        cap * slip[1] / slip[0].hypot(slip[1]),
    ];
    let reaction: Vec<_> = (0..4)
        .map(|i| n[i] * normal + tx[i] * tangent[0] + ty[i] * tangent[1])
        .collect();
    let correction = dense_solve(operator, reaction);
    body.wp = previous.to_vec();
    body.wn = (0..4).map(|i| final_pose[i] + correction[i]).collect();
    let mut contacts = vec![Contact::new(
        vec![(0, 0, 1.0), (0, 2, 0.3)],
        vec![(0, 1, 0.2)],
        0.0,
        Law::Power {
            k: stiffness,
            alpha: 1.5,
            mu: 0.0,
        },
    )];
    let terms = |column: &[f64; 4]| column.iter().enumerate().map(|(i, &w)| (0, i, w)).collect();
    let mut frictions = vec![Friction {
        contact: 0,
        sx: terms(&tx),
        sy: terms(&ty),
        mu: 0.25,
        force: 0.0,
        stuck: false,
    }];
    let report = solve_joint(
        std::slice::from_mut(&mut body),
        &mut contacts,
        &mut frictions,
        k,
    );
    for (actual, expected) in body.wn.iter().zip(final_pose) {
        assert!((actual - expected).abs() < 3e-12);
    }
    assert!((contacts[0].force / normal - 1.0).abs() < 3e-9);
    assert!((frictions[0].force / cap - 1.0).abs() < 3e-9);
    assert!(!frictions[0].stuck);
    check_normal_laws(std::slice::from_ref(&body), &contacts, k);
    let heat = 0.5 * (tangent[0] * slip[0] + tangent[1] * slip[1]);
    assert!((report.dissipated_j - heat).abs() < 3e-10);
}

fn coupled_pair(k: f64) -> (Vec<Body>, Vec<Contact>, Vec<Friction>) {
    let bodies = vec![
        point_with_prediction("normal a", 0.001, 0.001, 0.0, k),
        point_with_prediction("tangent a x", 0.001, 0.0001, -1000.0, k),
        point_with_prediction("tangent a y", 0.0002, 0.0, 0.0, k),
        point_with_prediction("normal b", 0.0008, 0.0008, 0.0, k),
        point_with_prediction("tangent b x", -0.0005, -0.0002, 250.0, k),
        point_with_prediction("tangent b y", -0.00015, 0.0, 0.0, k),
    ];
    let law = Law::Power {
        k: 1e6,
        alpha: 1.0,
        mu: 0.0,
    };
    let contacts = vec![
        Contact::new(vec![(0, 0, 1.0)], vec![], 0.0, law),
        Contact::new(vec![(3, 0, 1.0)], vec![], 0.0, law),
    ];
    let frictions = vec![
        Friction {
            contact: 0,
            sx: vec![(1, 0, 1.0), (3, 0, 0.2)],
            sy: vec![(2, 0, 1.0)],
            mu: 0.1,
            force: 0.0,
            stuck: false,
        },
        Friction {
            contact: 1,
            sx: vec![(4, 0, 1.0), (0, 0, 0.3)],
            sy: vec![(5, 0, 1.0), (2, 0, 0.2)],
            mu: 0.15,
            force: 0.0,
            stuck: false,
        },
    ];
    (bodies, contacts, frictions)
}

#[test]
fn interacting_friction_disks_are_order_independent_and_obey_both_final_laws() {
    let k = 0.001;
    let mut results = Vec::new();
    for reverse in [false, true] {
        let (mut bodies, mut contacts, mut frictions) = coupled_pair(k);
        let predicted: Vec<_> = bodies.iter().map(|b| b.wn[0]).collect();
        if reverse {
            contacts.reverse();
            frictions.reverse();
            for f in &mut frictions {
                f.contact = 1 - f.contact;
            }
        }
        let report = solve_joint(&mut bodies, &mut contacts, &mut frictions, k);
        if reverse {
            contacts.reverse();
            frictions.reverse();
            for f in &mut frictions {
                f.contact = 1 - f.contact;
            }
        }
        check_normal_laws(&bodies, &contacts, k);
        let response = |i: usize| (predicted[i] - bodies[i].wn[0]) / (k * k);
        let tangent = [
            [response(1), response(2) - 0.2 * response(5)],
            [response(4), response(5)],
        ];
        let mut heat = 0.0;
        for (f, q) in frictions.iter().zip(tangent) {
            let s = [centered_slip(&bodies, &f.sx), centered_slip(&bodies, &f.sy)];
            let magnitude = q[0].hypot(q[1]);
            let cap = f.mu * contacts[f.contact].force;
            assert!(magnitude <= cap + 2e-7);
            assert!((f.force - magnitude).abs() < 2e-7);
            if f.stuck {
                assert!(s[0].hypot(s[1]) < 2e-12);
            } else {
                assert!((magnitude - cap).abs() < 2e-7);
                assert!((q[0] * s[1] - q[1] * s[0]).abs() < 3e-10);
                assert!(q[0] * s[0] + q[1] * s[1] >= -1e-12);
            }
            heat += 0.5 * (q[0] * s[0] + q[1] * s[1]);
        }
        assert!((report.dissipated_j - heat).abs() < 3e-10);
        results.push((bodies, contacts, frictions));
    }
    for i in 0..6 {
        assert!((results[0].0[i].wn[0] - results[1].0[i].wn[0]).abs() < 3e-12);
    }
    for i in 0..2 {
        assert!((results[0].1[i].force - results[1].1[i].force).abs() < 3e-7);
        assert!((results[0].2[i].force - results[1].2[i].force).abs() < 3e-7);
    }
}

#[test]
fn refused_joint_step_preserves_every_contact_friction_and_crush_state() {
    let k = 0.001;
    let (mut bodies, mut contacts, mut frictions) = coupled_pair(k);
    contacts[0].law = Law::Power {
        k: 1e7,
        alpha: 1.5,
        mu: 0.0,
    };
    contacts[0].crush = Some(Crush::new(0.01, 1e5).unwrap());
    contacts[0].crush.as_mut().unwrap().deepest = 0.0002;
    for c in &mut contacts {
        c.force = 5.0;
        c.peak_force = 8.0;
        c.d_now = -0.001;
    }
    for f in &mut frictions {
        f.force = 4.0;
        f.stuck = true;
    }
    let old_contacts = contacts.clone();
    let old_friction = frictions.clone();
    let old_pose: Vec<_> = bodies.iter().map(|b| b.wn.clone()).collect();
    let mut controls = FrictionSolveControls::default();
    controls.max_newton_iterations = 0;
    let result = friction_response::solve_contacts_with_friction_checked(
        &mut bodies,
        &mut contacts,
        &mut frictions,
        &mut ContactQueue::default(),
        &mut contact_response::Workspace::default(),
        &mut friction_response::Workspace::default(),
        k,
        ContactSolveControls::default(),
        controls,
    );
    assert!(result.is_err());
    for (body, old) in bodies.iter().zip(old_pose) {
        assert_eq!(body.wn, old);
    }
    for (c, old) in contacts.iter().zip(old_contacts) {
        assert_eq!(
            (c.force, c.peak_force, c.d_now, c.crush),
            (old.force, old.peak_force, old.d_now, old.crush)
        );
    }
    for (f, old) in frictions.iter().zip(old_friction) {
        assert_eq!((f.force, f.stuck), (old.force, old.stuck));
    }
}
