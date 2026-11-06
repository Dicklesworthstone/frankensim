//! Physics checks for simultaneous normal contacts on the actual body adapter.

use super::*;
use contact_island::Controls;
use contact_response::{Workspace, solve_contacts_checked};

fn point_body(mass: f64, k: f64, position: f64, previous: f64) -> Body {
    let mut body = Body::point("contact verification", mass);
    body.w[0] = position;
    body.wp[0] = previous;
    body.predict(k);
    body
}

fn chain_contacts() -> Vec<Contact> {
    (0..2)
        .map(|i| {
            Contact::new(
                vec![(i, 0, 1.0)],
                vec![(i + 1, 0, 1.0)],
                0.0,
                Law::Power {
                    k: 1e6,
                    alpha: 1.5,
                    mu: 0.0,
                },
            )
        })
        .collect()
}

fn penetration(bodies: &[Body], contact: &Contact, frame: usize) -> f64 {
    let value = |b: usize, i: usize| match frame {
        0 => bodies[b].wn[i],
        1 => bodies[b].w[i],
        _ => bodies[b].wp[i],
    };
    contact
        .a
        .iter()
        .map(|&(b, i, w)| w * value(b, i))
        .sum::<f64>()
        - contact
            .b
            .iter()
            .map(|&(b, i, w)| w * value(b, i))
            .sum::<f64>()
        - contact.gap
}

fn elastic_energy(bodies: &[Body], contacts: &[Contact], k: f64, after: bool) -> f64 {
    let kinetic: f64 = bodies
        .iter()
        .map(|body| {
            let v = if after {
                (body.wn[0] - body.w[0]) / k
            } else {
                (body.w[0] - body.wp[0]) / k
            };
            0.5 * body.mass[0] * v * v
        })
        .sum();
    kinetic
        + contacts
            .iter()
            .map(|contact| {
                let Law::Power {
                    k: stiffness,
                    alpha,
                    ..
                } = contact.law;
                let potential = |d: f64| stiffness * d.max(0.0).powf(alpha + 1.0) / (alpha + 1.0);
                0.5 * (potential(penetration(bodies, contact, if after { 0 } else { 1 }))
                    + potential(penetration(bodies, contact, if after { 1 } else { 2 })))
            })
            .sum::<f64>()
}

fn check_final_laws(bodies: &[Body], contacts: &[Contact], k: f64) {
    for contact in contacts.iter().filter(|c| c.active) {
        let d = penetration(bodies, contact, 0);
        let previous = penetration(bodies, contact, 2);
        let Law::Power {
            k: stiffness,
            alpha,
            mu,
        } = contact.law;
        let expected = if let Some(crush) = &contact.crush {
            crush_phi(d, previous, stiffness, crush, mu, k).0
        } else if let Some(retardation) = contact.retardation_s {
            viscoelastic_phi(d, previous, stiffness, retardation, k).0
        } else {
            contact_phi(d, previous, stiffness, alpha, mu, k).0
        };
        assert!(contact.force >= 0.0);
        assert!(
            (contact.force - expected).abs() <= 5e-12 + 3e-10 * expected.abs(),
            "final force {} differs from its final-pose law {expected}",
            contact.force
        );
        assert!(
            (contact.d_now - d).abs() < 2e-14,
            "stale penetration diagnostic"
        );
    }
}

#[test]
fn simultaneous_three_mass_impact_preserves_momentum_and_staggered_energy() {
    let k = 0.001;
    let mut bodies = vec![
        point_body(1.0, k, 0.004, 0.004),
        point_body(0.7, k, 0.0, 0.0),
        point_body(1.5, k, -0.003, -0.003),
    ];
    let mut contacts = chain_contacts();
    let initial = elastic_energy(&bodies, &contacts, k, false);
    let mut queue = ContactQueue::default();
    let mut workspace = Workspace::default();
    let mut worst = 0.0_f64;
    for step in 0..600 {
        for body in &mut bodies {
            body.predict(k);
        }
        solve_contacts_checked(
            &mut bodies,
            &mut contacts,
            &mut queue,
            &mut workspace,
            k,
            Controls::default(),
        )
        .unwrap();
        check_final_laws(&bodies, &contacts, k);
        let energy = elastic_energy(&bodies, &contacts, k, true);
        worst = worst.max(((energy - initial) / initial).abs());
        let momentum: f64 = bodies
            .iter()
            .map(|b| b.mass[0] * (b.wn[0] - b.w[0]) / k)
            .sum();
        assert!(momentum.abs() < 2e-11, "step {step}: momentum {momentum:e}");
        for body in &mut bodies {
            std::mem::swap(&mut body.wp, &mut body.w);
            std::mem::swap(&mut body.w, &mut body.wn);
        }
    }
    eprintln!("simultaneous three-mass 600-step relative energy closure {worst:.12e}");
    assert!(contacts.iter().all(|c| c.peak_force > 0.0));
    assert!(worst < 2e-9, "lossless impact energy drift {worst:e}");
}

#[test]
fn one_contact_activates_an_open_contact_independently_of_contact_order() {
    let k = 0.001;
    let initial = vec![
        point_body(1.0, k, 0.01, 0.01),
        point_body(1.0, k, 0.0, 0.0),
        point_body(1.0, k, 0.0001, 0.0001),
    ];
    let mut results = Vec::new();
    for reverse in [false, true] {
        let mut bodies = initial.clone();
        let mut contacts = chain_contacts();
        assert!(penetration(&bodies, &contacts[1], 0) < 0.0);
        if reverse {
            contacts.reverse();
        }
        solve_contacts_checked(
            &mut bodies,
            &mut contacts,
            &mut ContactQueue::default(),
            &mut Workspace::default(),
            k,
            Controls::default(),
        )
        .unwrap();
        check_final_laws(&bodies, &contacts, k);
        if reverse {
            contacts.reverse();
        }
        assert!(
            contacts[1].force > 7.0,
            "cross-contact motion must close the initial gap"
        );
        assert!(contacts[1].d_now > 0.0);
        results.push((bodies, contacts));
    }
    for i in 0..3 {
        assert!((results[0].0[i].wn[0] - results[1].0[i].wn[0]).abs() < 2e-14);
    }
    for i in 0..2 {
        assert!((results[0].1[i].force - results[1].1[i].force).abs() < 2e-8);
    }
    eprintln!(
        "simultaneous activation forces {:?}",
        results[0].1.iter().map(|c| c.force).collect::<Vec<_>>()
    );
}

#[test]
fn simultaneous_hunt_crossley_and_retardation_close_dissipative_work() {
    let k = 0.001;
    let mut bodies = vec![
        point_body(1.0, k, 0.004, 0.004),
        point_body(0.7, k, 0.0, 0.0),
        point_body(1.5, k, -0.003, -0.003),
    ];
    let mut contacts = chain_contacts();
    contacts[0].law = Law::Power {
        k: 1e6,
        alpha: 1.5,
        mu: 0.4,
    };
    contacts[1].retardation_s = Some(0.0003);
    let initial = elastic_energy(&bodies, &contacts, k, false);
    let mut queue = ContactQueue::default();
    let mut workspace = Workspace::default();
    let (mut heat, mut worst) = (0.0, 0.0_f64);
    for _ in 0..160 {
        for body in &mut bodies {
            body.predict(k);
        }
        solve_contacts_checked(
            &mut bodies,
            &mut contacts,
            &mut queue,
            &mut workspace,
            k,
            Controls::default(),
        )
        .unwrap();
        check_final_laws(&bodies, &contacts, k);
        for contact in &contacts {
            let d = penetration(&bodies, contact, 0);
            let previous = penetration(&bodies, contact, 2);
            let elastic = 1e6 * power_secant(d, previous, 1.5).0;
            let loss = 0.5 * (contact.force - elastic) * (d - previous);
            assert!(loss > -2e-12, "contact generated work {loss:e}");
            heat += loss;
        }
        let energy = elastic_energy(&bodies, &contacts, k, true);
        worst = worst.max(((energy + heat - initial) / initial).abs());
        for body in &mut bodies {
            std::mem::swap(&mut body.wp, &mut body.w);
            std::mem::swap(&mut body.w, &mut body.wn);
        }
    }
    eprintln!(
        "simultaneous HC/retardation relative storage+heat closure {worst:.12e}; dissipated {} J",
        heat
    );
    assert!(heat > 0.0);
    assert!(worst < 2e-9);
}

pub(super) fn loaded_body(k: f64) -> Body {
    // Manufactured passive, nonorthogonal radiation columns; the operator,
    // cached inverse and contact laws remain the production implementations.
    let mut body = Body::point("loaded contact fixture", 1.0);
    body.n = 4;
    body.h = 0.6_f64.sqrt();
    body.area = 1.0;
    body.air_mass = 0.12;
    body.mass = vec![1.2, 0.7, 1.5, 0.9];
    body.s0 = vec![0.0; 4];
    body.w = vec![0.0; 4];
    body.wp = vec![0.0; 4];
    body.wn = vec![0.0; 4];
    body.tmp = vec![0.0; 4];
    body.g = vec![120.0, -80.0, 30.0, 50.0];
    body.sav = true;
    body.lap = Csr::from_rows(&vec![BTreeMap::new(); 4]);
    body.head_radiation = Some(Box::new(radiation::centered_tests::manufactured_load(
        vec![vec![1.0, 0.2, -0.1, 0.7], vec![0.4, 0.9, 0.3, -0.2]],
        vec![0.3, 0.15],
        vec![120.0, 75.0],
    )));
    body.coefs(k);
    body
}

/// Independent full nodal elimination, not a second Woodbury implementation.
pub(super) fn dense_solve(mut a: Vec<f64>, mut b: Vec<f64>) -> Vec<f64> {
    let n = b.len();
    for j in 0..n {
        let p = (j..n)
            .max_by(|&i, &k| a[i * n + j].abs().total_cmp(&a[k * n + j].abs()))
            .unwrap();
        for k in 0..n {
            a.swap(j * n + k, p * n + k);
        }
        b.swap(j, p);
        let diagonal = a[j * n + j];
        for k in j..n {
            a[j * n + k] /= diagonal;
        }
        b[j] /= diagonal;
        for i in 0..n {
            if i == j {
                continue;
            }
            let coefficient = a[i * n + j];
            for k in j..n {
                a[i * n + k] -= coefficient * a[j * n + k];
            }
            b[i] -= coefficient * b[j];
        }
    }
    b
}

#[test]
fn coupled_loaded_head_matches_full_dense_sav_air_and_radiation_response() {
    let k = 0.002;
    let mut body = loaded_body(k);
    let load = body.head_radiation.as_ref().unwrap();
    let modes = load.cols.len();
    let mut operator = vec![0.0; 16];
    for i in 0..4 {
        for j in 0..4 {
            let mut entry = body.air_mass * body.h.powi(4) / (k * k) + body.g[i] * body.g[j] / 4.0;
            if i == j {
                entry += 1.0 / body.cf[i][0];
            }
            for r in 0..modes {
                for s in 0..modes {
                    entry += load.cols[r][i]
                        * (load.mass[r * modes + s] / (k * k)
                            + load.damping[r * modes + s] / (2.0 * k))
                        * load.cols[s][j];
                }
            }
            operator[i * 4 + j] = entry;
        }
    }
    let ports = [[1.0, -0.2, 0.3, 0.0], [0.0, 0.3, 0.5, -0.8]];
    let response: Vec<_> = ports
        .iter()
        .map(|p| dense_solve(operator.clone(), p.to_vec()))
        .collect();
    let mut contacts = vec![
        Contact::new(
            vec![(0, 0, 0.4), (0, 0, 0.6), (0, 2, 0.3)],
            vec![(0, 1, 0.2)],
            0.0,
            Law::Power {
                k: 1e8,
                alpha: 1.5,
                mu: 0.0,
            },
        ),
        Contact::new(
            vec![(0, 1, 0.3), (0, 2, 0.5)],
            vec![(0, 3, 0.8)],
            0.0,
            Law::Power {
                k: 7e7,
                alpha: 1.5,
                mu: 0.0,
            },
        ),
    ];
    // Check every self and cross entry through the actual prepared inverse.
    for p in 0..2 {
        body.footprint_response(0, &contacts[p].a, 1.0, &contacts[p].b, -1.0, k);
        for q in 0..2 {
            let actual = dot(&ports[q], &body.tmp);
            let expected = dot(&ports[q], &response[p]);
            assert!((actual - expected).abs() < 2e-12 * expected.abs().max(1e-12));
        }
    }
    body.wp = vec![0.0002, 0.0003, 0.0001, -0.0001];
    body.wn = vec![0.0004, 0.0002, 0.0003, -0.0002];
    let predicted = body.wn.clone();
    let factor = body
        .head_radiation
        .as_ref()
        .unwrap()
        .prepared
        .as_ref()
        .unwrap()
        .chol
        .as_ptr();
    solve_contacts_checked(
        std::slice::from_mut(&mut body),
        &mut contacts,
        &mut ContactQueue::default(),
        &mut Workspace::default(),
        k,
        Controls::default(),
    )
    .unwrap();
    check_final_laws(std::slice::from_ref(&body), &contacts, k);
    for i in 0..4 {
        let expected = predicted[i]
            - (0..2)
                .map(|p| response[p][i] * contacts[p].force)
                .sum::<f64>();
        assert!((body.wn[i] - expected).abs() < 2e-14);
    }
    assert_eq!(
        factor,
        body.head_radiation
            .as_ref()
            .unwrap()
            .prepared
            .as_ref()
            .unwrap()
            .chol
            .as_ptr()
    );
    assert!(contacts.iter().all(|c| c.force > 0.0));
}

#[test]
fn contact_refusal_preserves_predicted_pose_forces_and_crush_memory() {
    let k = 0.001;
    let mut bodies = vec![
        point_body(1.0, k, 0.004, 0.004),
        point_body(0.7, k, 0.0, 0.0),
        point_body(1.5, k, -0.003, -0.003),
    ];
    let mut contacts = chain_contacts();
    contacts[0].crush = Some(Crush::new(0.01, 3e6).unwrap());
    contacts[0].peak_force = 9.0;
    contacts[0].d_now = -0.001;
    let old_contacts = contacts.clone();
    let old_pose: Vec<_> = bodies.iter().map(|b| b.wn.clone()).collect();
    let controls = Controls {
        max_newton_iterations: 0,
        ..Controls::default()
    };
    let refused = solve_contacts_checked(
        &mut bodies,
        &mut contacts,
        &mut ContactQueue::default(),
        &mut Workspace::default(),
        k,
        controls,
    );
    assert!(refused.is_err());
    for (body, previous) in bodies.iter().zip(old_pose) {
        assert_eq!(body.wn, previous);
    }
    for (contact, previous) in contacts.iter().zip(old_contacts) {
        assert_eq!(contact.force, previous.force);
        assert_eq!(contact.peak_force, previous.peak_force);
        assert_eq!(contact.d_now, previous.d_now);
        assert_eq!(contact.crush, previous.crush);
    }
}

#[test]
fn simultaneous_loaded_contacts_close_energy_and_radiation_work() {
    let k = 1.0 / 48_000.0;
    let mut body = loaded_body(k);
    // A frozen arbitrary SAV gradient is useful for an inverse check, but
    // not an energy trajectory. Here storage is the explicit physical mass
    // (dry + modal air + volume air), linear springs and Hertz potentials.
    body.sav = false;
    let load = body.head_radiation.as_ref().unwrap();
    let modes = load.cols.len();
    let mut mass = vec![0.0; 16];
    for i in 0..4 {
        for j in 0..4 {
            mass[4 * i + j] = body.air_mass * body.h.powi(4);
            if i == j {
                mass[4 * i + j] += body.mass[i];
            }
            for r in 0..modes {
                for s in 0..modes {
                    mass[4 * i + j] += load.cols[r][i] * load.mass[r * modes + s] * load.cols[s][j];
                }
            }
        }
    }
    let stiffness: Vec<_> = body
        .mass
        .iter()
        .map(|m| m * (2.0 * std::f64::consts::PI * 700.0).powi(2))
        .collect();
    let ports = [[1.0, -0.2, 0.3, 0.0], [0.0, 0.3, 0.5, -0.8]];
    let potential = |d: f64| 0.4 * 1e10 * d.max(0.0).powf(2.5);
    let energy = |w: &[f64], wp: &[f64]| {
        let v: Vec<_> = w.iter().zip(wp).map(|(a, b)| (a - b) / k).collect();
        0.5 * (0..4)
            .map(|i| (0..4).map(|j| v[i] * mass[4 * i + j] * v[j]).sum::<f64>())
            .sum::<f64>()
            + 0.5 * (0..4).map(|i| stiffness[i] * w[i] * wp[i]).sum::<f64>()
            + 0.5
                * ports
                    .iter()
                    .map(|p| potential(dot(p, w)) + potential(dot(p, wp)))
                    .sum::<f64>()
    };
    let law = Law::Power {
        k: 1e10,
        alpha: 1.5,
        mu: 0.0,
    };
    let mut contacts = vec![
        Contact::new(
            vec![(0, 0, 0.4), (0, 0, 0.6), (0, 2, 0.3)],
            vec![(0, 1, 0.2)],
            0.0,
            law,
        ),
        Contact::new(vec![(0, 1, 0.3), (0, 2, 0.5)], vec![(0, 3, 0.8)], 0.0, law),
    ];
    body.w = vec![2e-5, 0.0, 1e-5, -2e-5];
    body.wp = body.w.clone();
    let initial = energy(&body.w, &body.wp);
    let (mut heat, mut worst) = (0.0, 0.0_f64);
    let mut queue = ContactQueue::default();
    let mut workspace = Workspace::default();
    for _ in 0..1_200 {
        let h2 = body.h * body.h;
        let mu = body.air_mass / (body.area * body.area * k * k);
        let volume_previous = h2
            * body
                .w
                .iter()
                .zip(&body.wp)
                .map(|(a, b)| 2.0 * a - b)
                .sum::<f64>();
        let load = body.head_radiation.as_deref_mut().unwrap();
        load.load_rhs(&body.w, &body.wp, k);
        for i in 0..4 {
            body.wn[i] = body.mass[i] * (2.0 * body.w[i] - body.wp[i]) / (k * k)
                - stiffness[i] * body.w[i]
                + load.rhs[i]
                + mu * h2 * volume_previous;
        }
        solve_body_response(
            &body.cf,
            &body.g,
            false,
            h2,
            mu,
            Some(load),
            &mut body.wn,
            k,
        );
        solve_contacts_checked(
            std::slice::from_mut(&mut body),
            &mut contacts,
            &mut queue,
            &mut workspace,
            k,
            Controls::default(),
        )
        .unwrap();
        check_final_laws(std::slice::from_ref(&body), &contacts, k);
        heat += body
            .head_radiation
            .as_deref_mut()
            .unwrap()
            .radiated(&body.wn, k);
        let storage = energy(&body.wn, &body.w);
        worst = worst.max(((storage + heat - initial) / initial).abs());
        std::mem::swap(&mut body.wp, &mut body.w);
        std::mem::swap(&mut body.w, &mut body.wn);
    }
    eprintln!(
        "simultaneous loaded 1200-step relative storage+radiation closure {worst:.12e}; radiated {} J",
        heat
    );
    assert!(heat > 0.0);
    assert!(contacts.iter().all(|c| c.peak_force > 0.0));
    assert!(
        worst < 2e-9,
        "loaded simultaneous contact energy defect {worst:e}"
    );
}

#[test]
fn instrument_contact_refusal_preserves_clock_output_and_requires_restoration() {
    let k = 0.001;
    let mut instrument = Instrument::new(1.0 / k);
    instrument.add(point_body(1.0, k, 0.004, 0.004));
    instrument.add(point_body(0.7, k, 0.0, 0.0));
    instrument.add(point_body(1.5, k, -0.003, -0.003));
    instrument.contacts = chain_contacts();
    instrument.mics.push(Mic::new([0.0, 0.0, 1.0]));
    instrument.prepare();
    let mut restored = instrument.clone();
    let previous_frames: Vec<_> = instrument
        .bodies
        .iter()
        .map(|b| (b.w.clone(), b.wp.clone()))
        .collect();
    let previous_ring = instrument.mics[0].ring.clone();
    instrument.contact_solver_controls.max_newton_iterations = 0;
    let mut output = [17.25];
    let error = instrument.try_step(&mut output).unwrap_err();
    assert!(matches!(error, ContactSolveError::NoConvergence { .. }));
    assert_eq!(instrument.step_index, 0);
    assert_eq!(instrument.time(), 0.0);
    assert_eq!(output, [17.25]);
    assert_eq!(instrument.mics[0].ring, previous_ring);
    for (body, (w, wp)) in instrument.bodies.iter().zip(previous_frames) {
        assert_eq!(body.w, w);
        assert_eq!(body.wp, wp);
    }
    assert!(
        instrument
            .contacts
            .iter()
            .all(|c| c.force == 0.0 && c.peak_force == 0.0)
    );

    // Increasing the budget cannot resume histories advanced before refusal.
    // A second attempt returns before prediction or microphone publication.
    instrument.contact_solver_controls = ContactSolveControls::default();
    let refused_pose: Vec<_> = instrument.bodies.iter().map(|b| b.wn.clone()).collect();
    output[0] = -23.5;
    assert_eq!(instrument.try_step(&mut output), Err(error));
    assert_eq!(instrument.step_index, 0);
    assert_eq!(instrument.time(), 0.0);
    assert_eq!(output, [-23.5]);
    assert_eq!(instrument.mics[0].ring, previous_ring);
    for (body, pose) in instrument.bodies.iter().zip(refused_pose) {
        assert_eq!(body.wn, pose);
    }

    restored.try_step(&mut output).unwrap();
    assert_eq!(restored.step_index, 1);
    assert_eq!(restored.time(), k);
    assert!(output[0].is_finite());
    assert!(restored.contacts.iter().all(|c| c.force > 0.0));
}

#[test]
fn independent_scalar_candidates_match_coupled_physics_and_keep_sparse_publication() {
    let k = 0.001;
    let initial: Vec<_> = (0..8)
        .map(|i| {
            let position = 0.001 + 0.0001 * i as f64;
            point_body(0.7 + 0.1 * i as f64, k, position, position)
        })
        .collect();
    let contacts: Vec<_> = (0..7)
        .map(|i| {
            Contact::new(
                vec![(i, 0, 1.0), (7, 0, 0.0)],
                vec![],
                0.0,
                Law::Power {
                    k: 1e6 + 1e5 * i as f64,
                    alpha: 1.5,
                    mu: 0.0,
                },
            )
        })
        .collect();
    let mut results = Vec::new();
    for limit in [32, 31] {
        let mut bodies = initial.clone();
        bodies[7].tmp[0] = 91.0;
        let mut contacts = contacts.clone();
        let controls = Controls {
            max_newton_iterations: limit,
            ..Controls::default()
        };
        let report = solve_contacts_checked(
            &mut bodies,
            &mut contacts,
            &mut ContactQueue::default(),
            &mut Workspace::default(),
            k,
            controls,
        )
        .unwrap();
        check_final_laws(&bodies, &contacts, k);
        assert_eq!(bodies[7].wn, initial[7].wn);
        if limit == 32 {
            assert_eq!(report.cg_iterations, 0);
            assert_eq!(report.compliance_applications, 1);
            assert_eq!(
                bodies[7].tmp[0], 91.0,
                "sparse publication must skip zero weights"
            );
        } else {
            assert!(report.cg_iterations > 0);
        }
        results.push((bodies, contacts));
    }
    for i in 0..7 {
        assert!((results[0].0[i].wn[0] - results[1].0[i].wn[0]).abs() < 2e-14);
        assert!((results[0].1[i].force - results[1].1[i].force).abs() < 2e-8);
    }
}

#[test]
fn loaded_scalar_path_requires_true_contact_independence() {
    let k = 0.002;
    let mut initial = loaded_body(k);
    initial.wp = vec![0.0002, 0.0003, 0.0001, -0.0001];
    initial.wn = vec![0.0004, 0.0002, 0.0003, -0.0002];
    let law = Law::Power {
        k: 1e8,
        alpha: 1.5,
        mu: 0.0,
    };
    let contact = Contact::new(vec![(0, 0, 1.0), (0, 2, 0.3)], vec![(0, 1, 0.2)], 0.0, law);
    let mut results = Vec::new();
    for limit in [32, 31] {
        let mut body = initial.clone();
        let mut contacts = vec![contact.clone()];
        let controls = Controls {
            max_newton_iterations: limit,
            ..Controls::default()
        };
        let report = solve_contacts_checked(
            std::slice::from_mut(&mut body),
            &mut contacts,
            &mut ContactQueue::default(),
            &mut Workspace::default(),
            k,
            controls,
        )
        .unwrap();
        check_final_laws(std::slice::from_ref(&body), &contacts, k);
        if limit == 32 {
            assert_eq!(report.cg_iterations, 0);
            assert_eq!(report.compliance_applications, 2);
        } else {
            assert!(report.cg_iterations > 0);
        }
        results.push((body, contacts));
    }
    for i in 0..4 {
        assert!((results[0].0.wn[i] - results[1].0.wn[i]).abs() < 2e-14);
    }
    assert!((results[0].1[0].force - results[1].1[0].force).abs() < 2e-8);

    // Distinct nodes are coupled through the same nonlocal loaded body.
    let mut contacts = vec![
        Contact::new(vec![(0, 0, 1.0)], vec![], 0.0, law),
        Contact::new(vec![(0, 1, 1.0)], vec![], 0.0, law),
    ];
    let mut body = initial.clone();
    let report = solve_contacts_checked(
        std::slice::from_mut(&mut body),
        &mut contacts,
        &mut ContactQueue::default(),
        &mut Workspace::default(),
        k,
        Controls::default(),
    )
    .unwrap();
    check_final_laws(std::slice::from_ref(&body), &contacts, k);
    assert!(
        report.cg_iterations > 0,
        "nonlocal cross-compliance cannot use independent scalar laws"
    );

    // Repeated terms also use the aggregate operator, so a diagonal sum of
    // squared terms cannot silently replace the square of their signed sum.
    let mut body = point_body(1.0, k, 0.001, 0.001);
    let mut contacts = vec![Contact::new(
        vec![(0, 0, 0.4), (0, 0, 0.6)],
        vec![],
        0.0,
        law,
    )];
    let report = solve_contacts_checked(
        std::slice::from_mut(&mut body),
        &mut contacts,
        &mut ContactQueue::default(),
        &mut Workspace::default(),
        k,
        Controls::default(),
    )
    .unwrap();
    check_final_laws(std::slice::from_ref(&body), &contacts, k);
    assert!(report.cg_iterations > 0);

    let mut body = point_body(1.0, k, 0.0, 0.0);
    let mut contacts = vec![Contact::new(
        vec![(0, 0, 1.0)],
        vec![(0, 0, 1.0)],
        -0.001,
        law,
    )];
    let report = solve_contacts_checked(
        std::slice::from_mut(&mut body),
        &mut contacts,
        &mut ContactQueue::default(),
        &mut Workspace::default(),
        k,
        Controls::default(),
    )
    .unwrap();
    check_final_laws(std::slice::from_ref(&body), &contacts, k);
    assert_eq!(
        body.wn[0], 0.0,
        "opposite signed terms have zero net footprint"
    );
    // Zero C needs no Newton correction, but must still use the aggregate
    // operator: its nonzero-force action and final publication check. The
    // all-zero warm force needs no operator application.
    assert_eq!(report.cg_iterations, 0);
    assert_eq!(report.compliance_applications, 2);
    assert!((contacts[0].force - 1e8 * 0.001_f64.powf(1.5)).abs() < 1e-9);
}

fn mixed_contact_components(k: f64) -> (Vec<Body>, Vec<Contact>) {
    let mut dry = loaded_body(k);
    dry.sav = false;
    dry.air_mass = 0.0;
    dry.head_radiation = None;
    dry.wp = vec![0.004, 0.00013, 0.0007, 0.0008];
    dry.w = dry.wp.clone();
    dry.wn = dry.wp.clone();
    dry.tmp = vec![35.0, 91.0, -42.0, 17.0];
    let mut loaded_a = loaded_body(k);
    loaded_a.wp = vec![0.0001, 0.0002, 0.0001, 0.0];
    loaded_a.wn = vec![0.0002, 0.0001, 0.0001, -0.0001];
    let mut loaded_b = loaded_body(k);
    loaded_b.wp = vec![0.0002, 0.0001, 0.0001, -0.0001];
    loaded_b.wn = vec![0.0001, 0.0002, 0.0002, 0.0];
    let bodies = vec![
        dry,
        point_body(0.7, k, 0.0, 0.0),
        point_body(1.5, k, -0.003, -0.003),
        loaded_a,
        loaded_b,
    ];
    let mut contacts = chain_contacts();
    contacts.push(Contact::new(
        vec![(0, 2, 1.0)],
        vec![(3, 0, 1.0), (3, 2, 0.2)],
        0.0,
        Law::Power {
            k: 5e7,
            alpha: 1.5,
            mu: 0.0,
        },
    ));
    contacts.push(Contact::new(
        vec![(0, 3, 1.0)],
        vec![(4, 1, 0.7), (4, 3, -0.1)],
        0.0,
        Law::Power {
            k: 2e7,
            alpha: 1.5,
            mu: 0.0,
        },
    ));
    (bodies, contacts)
}

/// Assemble every mutual response, then solve one global nonlinear system.
/// This oracle has no contact graph, component detection or scalar shortcut.
fn global_contact_oracle(
    bodies: &[Body],
    contacts: &[Contact],
    k: f64,
) -> (Vec<Vec<f64>>, Vec<f64>) {
    let count = contacts.len();
    let mut probe = bodies.to_vec();
    let mut images = Vec::with_capacity(count);
    let mut compliance = vec![0.0; count * count];
    for (column, contact) in contacts.iter().enumerate() {
        for (bi, body) in probe.iter_mut().enumerate() {
            body.footprint_response(bi, &contact.a, 1.0, &contact.b, -1.0, k);
        }
        for (row, receiver) in contacts.iter().enumerate() {
            compliance[row * count + column] = receiver
                .a
                .iter()
                .map(|&(b, i, w)| w * probe[b].tmp[i])
                .sum::<f64>()
                - receiver
                    .b
                    .iter()
                    .map(|&(b, i, w)| w * probe[b].tmp[i])
                    .sum::<f64>();
        }
        images.push(probe.iter().map(|b| b.tmp.clone()).collect::<Vec<_>>());
    }
    let predicted: Vec<_> = contacts.iter().map(|c| penetration(bodies, c, 0)).collect();
    let previous: Vec<_> = contacts.iter().map(|c| penetration(bodies, c, 2)).collect();
    let warm: Vec<_> = contacts.iter().map(|c| c.force).collect();
    let diagonal: Vec<_> = (0..count).map(|i| compliance[i * count + i]).collect();
    let mut global = contact_island::Workspace::default();
    global
        .solve(
            &predicted,
            &warm,
            &diagonal,
            Controls::default(),
            1.0,
            |x, y| {
                for (i, value) in y.iter_mut().enumerate() {
                    *value = dot(&compliance[i * count..(i + 1) * count], x);
                }
            },
            |i, d| {
                let Law::Power {
                    k: stiffness,
                    alpha,
                    mu,
                } = contacts[i].law;
                contact_phi(d, previous[i], stiffness, alpha, mu, k)
            },
        )
        .unwrap();
    let forces = global.forces().to_vec();
    let pose = bodies
        .iter()
        .enumerate()
        .map(|(b, body)| {
            body.wn
                .iter()
                .enumerate()
                .map(|(i, predicted)| {
                    predicted - (0..count).map(|p| images[p][b][i] * forces[p]).sum::<f64>()
                })
                .collect()
        })
        .collect();
    (pose, forces)
}

#[test]
fn independent_islands_match_global_core_when_sharing_distinct_dry_body_nodes() {
    let k = 0.001;
    let (initial, contacts) = mixed_contact_components(k);
    let (expected_pose, expected_forces) = global_contact_oracle(&initial, &contacts, k);
    for reverse in [false, true] {
        let mut bodies = initial.clone();
        let mut contacts = contacts.clone();
        if reverse {
            contacts.reverse();
        }
        let report = solve_contacts_checked(
            &mut bodies,
            &mut contacts,
            &mut ContactQueue::default(),
            &mut Workspace::default(),
            k,
            Controls::default(),
        )
        .unwrap();
        check_final_laws(&bodies, &contacts, k);
        if reverse {
            contacts.reverse();
        }
        for (body, expected) in bodies.iter().zip(&expected_pose) {
            for (actual, expected) in body.wn.iter().zip(expected) {
                assert!(
                    (actual - expected).abs() < 3e-14,
                    "component correction contaminated another response: {actual:e} vs {expected:e}"
                );
            }
        }
        for (contact, expected) in contacts.iter().zip(&expected_forces) {
            assert!((contact.force - expected).abs() < 2e-8 + 2e-10 * expected.abs());
        }
        assert_eq!(
            bodies[0].wn[1], initial[0].wn[1],
            "an untouched diagonal node must not receive another island's scratch"
        );
        assert!(
            report.cg_iterations > 0,
            "the shared-node Hertz pair needs the coupled owner"
        );
    }
}

#[test]
fn later_island_refusal_publishes_no_earlier_pose_force_or_crush_history() {
    let k = 0.001;
    let controls = Controls {
        max_newton_iterations: 0,
        ..Controls::default()
    };
    let early_body = point_body(1.0, k, 0.004, 0.004);
    let law = Law::Power {
        k: 1e6,
        alpha: 1.5,
        mu: 0.0,
    };
    let crush = Crush::new(0.01, 3e6).unwrap();
    let frozen = Crush {
        deepest: 0.004,
        ..crush
    };
    let mut early = Contact::new(vec![(0, 0, 1.0)], vec![], 0.0, law)
        .crushing(crush)
        .unwrap();
    early.force = crush_force(0.004, k * k, 0.0, 0.004, law, &frozen, true, k);
    early.peak_force = 1.0;
    early.d_now = -0.001;
    // Prove this first island can accept the same zero-Newton budget and
    // would publish nonzero displacement, force diagnostics and crush memory.
    let mut isolated_body = early_body.clone();
    let mut isolated_contact = early.clone();
    solve_contacts_checked(
        std::slice::from_mut(&mut isolated_body),
        std::slice::from_mut(&mut isolated_contact),
        &mut ContactQueue::default(),
        &mut Workspace::default(),
        k,
        controls,
    )
    .unwrap();
    assert!(isolated_body.wn[0] < early_body.wn[0]);
    assert!(isolated_contact.peak_force > early.peak_force);
    assert!(isolated_contact.crush.unwrap().deepest > crush.deepest);

    let mut bodies = vec![
        early_body,
        point_body(1.0, k, 0.004, 0.004),
        point_body(0.7, k, 0.0, 0.0),
        point_body(1.5, k, -0.003, -0.003),
    ];
    let mut contacts = vec![early];
    contacts.extend(chain_contacts().into_iter().map(|mut c| {
        for (b, _, _) in c.a.iter_mut().chain(&mut c.b) {
            *b += 1;
        }
        c
    }));
    let previous_pose: Vec<_> = bodies.iter().map(|b| b.wn.clone()).collect();
    let previous_contacts = contacts.clone();
    let result = solve_contacts_checked(
        &mut bodies,
        &mut contacts,
        &mut ContactQueue::default(),
        &mut Workspace::default(),
        k,
        controls,
    );
    assert!(matches!(
        result,
        Err(contact_island::SolveError::NoConvergence { .. })
    ));
    for (body, previous) in bodies.iter().zip(previous_pose) {
        assert_eq!(body.wn, previous);
    }
    for (contact, previous) in contacts.iter().zip(previous_contacts) {
        assert_eq!(contact.force, previous.force);
        assert_eq!(contact.peak_force, previous.peak_force);
        assert_eq!(contact.d_now, previous.d_now);
        assert_eq!(contact.crush, previous.crush);
    }
}

/// The normal law's version of the same defect: a warm force at the edge of the
/// core's acceptance passes the core, and roundoff can refuse it at the rounded
/// pose. That refusal goes back to Newton until the audit passes. Without budget
/// it stays a named refusal that publishes nothing.
#[test]
fn rounded_pose_refusal_at_the_acceptance_edge_iterates_until_the_audit_passes() {
    let k = 1e-5;
    let step = |state: &mut (Vec<Body>, Vec<Contact>), controls| {
        solve_contacts_checked(
            &mut state.0,
            &mut state.1,
            &mut ContactQueue::default(),
            &mut Workspace::default(),
            k,
            controls,
        )
    };
    let mut retried = 0;
    for member in 0..32 {
        let offset = 1e-2 * (1.0 + f64::from(member) / 32.0);
        let penetration = 1.5e-8;
        let mut contacts = chain_contacts();
        for c in &mut contacts {
            c.law = Law::Power {
                k: 1e8,
                alpha: 1.0,
                mu: 0.0,
            };
        }
        let pristine = (
            vec![
                point_body(1.0, k, offset + penetration, offset + penetration),
                point_body(1.0, k, offset, offset),
                point_body(1.0, k, offset - penetration, offset - penetration),
            ],
            contacts,
        );
        let mut settled = pristine.clone();
        step(&mut settled, Controls::default()).unwrap();
        let warmed = |excess: f64| {
            let mut state = pristine.clone();
            for (c, old) in state.1.iter_mut().zip(&settled.1) {
                c.force = old.force * (1.0 + excess);
            }
            state
        };
        // Bisect for the largest warm error the core accepts without a step.
        let (mut inside, mut outside) = (0.0, 1e-6);
        for _ in 0..48 {
            let excess = 0.5 * (inside + outside);
            let report = step(&mut warmed(excess), Controls::default()).unwrap();
            if report.newton_iterations == 0 {
                inside = excess;
            } else {
                outside = excess;
            }
        }
        let mut state = warmed(outside);
        let report = step(&mut state, Controls::default()).unwrap();
        if report.rounded_pose_retries == 0 {
            continue;
        }
        retried += 1;
        assert!(report.newton_iterations > 0);
        assert!(report.force_residual <= 1.0 && report.displacement_residual <= 1.0);

        let mut refused = warmed(outside);
        let exhausted = Controls {
            max_newton_iterations: 0,
            ..Controls::default()
        };
        let error = step(&mut refused, exhausted).unwrap_err();
        assert!(matches!(
            error,
            ContactSolveError::NoConvergence {
                site: ContactRefusalSite::RoundedPose,
                newton_iterations: 0,
                ..
            }
        ));
        for (body, old) in refused.0.iter().zip(&pristine.0) {
            assert_eq!(body.wn, old.wn);
        }
        for (c, old) in refused.1.iter().zip(&warmed(outside).1) {
            assert_eq!(c.force, old.force);
        }
    }
    assert!(retried > 0, "no family member reached the rounded-pose refusal");
}
