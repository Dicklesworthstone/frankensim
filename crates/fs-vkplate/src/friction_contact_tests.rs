//! Contact/friction checks at the centered velocity used by the body energy balance.

use super::*;

pub(super) fn point_with_prediction(
    name: &str,
    position: f64,
    previous: f64,
    force: f64,
    k: f64,
) -> Body {
    let mut body = Body::point(name, 1.0);
    body.w[0] = position;
    body.wp[0] = previous;
    body.f[0] = force;
    body.predict(k);
    body
}

pub(super) fn turning_oscillator(k: f64) -> (Vec<Body>, Vec<Contact>, Vec<Friction>) {
    let stiffness = 1e6;
    let bodies = vec![
        point_with_prediction("normal", 0.001, 0.001, 0.0, k),
        point_with_prediction("tangent x", 0.001, 0.0001, -stiffness * 0.001, k),
        point_with_prediction("tangent y", 0.0, 0.0, 0.0, k),
    ];
    let contacts = vec![Contact::new(
        vec![(0, 0, 1.0)],
        vec![],
        0.0,
        Law::Power {
            k: stiffness,
            alpha: 1.0,
            mu: 0.0,
        },
    )];
    let frictions = vec![Friction {
        contact: 0,
        sx: vec![(1, 0, 1.0)],
        sy: vec![(2, 0, 1.0)],
        mu: 0.075,
        force: 0.0,
        stuck: false,
    }];
    (bodies, contacts, frictions)
}

pub(super) fn turning_energy(bodies: &[Body], k: f64, after: bool) -> f64 {
    let kinetic: f64 = bodies
        .iter()
        .map(|b| {
            let increment = if after {
                b.wn[0] - b.w[0]
            } else {
                b.w[0] - b.wp[0]
            };
            0.5 * b.mass[0] * (increment / k).powi(2)
        })
        .sum();
    let normal = &bodies[0];
    let (a, b) = if after {
        (normal.wn[0], normal.w[0])
    } else {
        (normal.w[0], normal.wp[0])
    };
    let contact = 0.25 * 1e6 * (a.max(0.0).powi(2) + b.max(0.0).powi(2));
    let tangent = &bodies[1];
    let spring = 0.5
        * 1e6
        * if after {
            tangent.wn[0] * tangent.w[0]
        } else {
            tangent.w[0] * tangent.wp[0]
        };
    kinetic + contact + spring
}

#[test]
fn centered_friction_cannot_pump_a_turning_resonant_oscillator() {
    // The last half-step and the centered velocity have opposite signs near
    // this turning point. Friction based on wn-w would apply +50 N and add
    // 0.02125 J. Centered work requires -50 N and removes 0.01875 J.
    let k = 0.001;
    let (mut bodies, mut contacts, mut frictions) = turning_oscillator(k);
    let initial = turning_energy(&bodies, k, false);
    assert!((initial - 0.955).abs() < 2e-14);
    contact_response::solve_contacts_checked(
        &mut bodies,
        &mut contacts,
        &mut ContactQueue::default(),
        &mut contact_response::Workspace::default(),
        k,
        ContactSolveControls::default(),
    )
    .unwrap();
    assert!((contacts[0].force - 2000.0 / 3.0).abs() < 2e-9);
    let predicted_tangent = bodies[1].wn[0];
    assert!(predicted_tangent - bodies[1].w[0] < 0.0);
    assert!(predicted_tangent - bodies[1].wp[0] > 0.0);

    solve_frictions(&mut bodies, &contacts, &mut frictions, k);

    // Recover the signed applied force from the body update, independently
    // of the friction diagnostic (which stores only its magnitude).
    let tangent_force = bodies[1].mass[0] * (bodies[1].wn[0] - predicted_tangent) / (k * k);
    let centered_increment = bodies[1].wn[0] - bodies[1].wp[0];
    let heat = -0.5 * tangent_force * centered_increment;
    let final_energy = turning_energy(&bodies, k, true);
    assert!(!frictions[0].stuck);
    assert!((tangent_force + 50.0).abs() < 2e-9);
    assert!((frictions[0].force - 50.0).abs() < 2e-9);
    assert!((bodies[1].wn[0] - 0.00085).abs() < 2e-15);
    assert!(heat > 0.0);
    assert!((heat - 0.01875).abs() < 2e-12);
    assert!((final_energy + heat - initial).abs() < 3e-12);
    eprintln!(
        "turning oscillator: energy {initial:.12} -> {final_energy:.12} J, centered friction heat {heat:.12} J"
    );
}
