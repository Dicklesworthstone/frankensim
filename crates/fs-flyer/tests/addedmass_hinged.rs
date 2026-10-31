//! G0/G1/G3: finite hinge kinematics, configuration inertia, and joint power.
//! Repro: cargo test -p fs-flyer --test addedmass_hinged

use fs_flyer::addedmass::hinged::{HingedStrip, assemble_hinged_strips};
use fs_flyer::addedmass::{
    AeroGeneralizedLoads, MAX_CONTROL_COORDS, MAX_STRIPS, Strip, assemble_analytic_strip,
    solve_joint,
};
use fs_flyer::spine::RigidBody;

const RHO: f64 = 1.3;

fn plate() -> Strip {
    Strip {
        name: "hinged plate",
        chord_m: 1.2,
        span_m: 2.6,
        position_m: [0.8, 0.0, 0.0],
        normal: [0.0, 0.0, 1.0],
        span_axis: [0.0, 1.0, 0.0],
        control_coord: None,
        control_gain: 0.0,
        control_angular_gain: [0.0; 3],
    }
}

fn close(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance * a.abs().max(b.abs()).max(1.0),
        "{a} != {b}"
    );
}

fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

fn matrix(loads: &AeroGeneralizedLoads) -> Vec<Vec<f64>> {
    let nc = loads.m_added_cc.len();
    let mut m = vec![vec![0.0; 6 + nc]; 6 + nc];
    for i in 0..6 {
        for j in 0..6 {
            m[i][j] = loads.m_added_rr[i][j];
        }
        for k in 0..nc {
            m[i][6 + k] = loads.m_added_rc[k][i];
            m[6 + k][i] = loads.m_added_rc[k][i];
        }
    }
    for k in 0..nc {
        for j in 0..nc {
            m[6 + k][6 + j] = loads.m_added_cc[k][j];
        }
    }
    m
}

fn apply(m: &[Vec<f64>], u: &[f64]) -> Vec<f64> {
    m.iter()
        .map(|row| row.iter().zip(u).map(|(x, y)| x * y).sum())
        .collect()
}

#[test]
fn finite_pitch_hinge_has_exact_geometry_and_centrifugal_reaction() {
    // Independent planar rigid rotation, not an infinitesimal gain model.
    let source = plate();
    let radius = source.position_m[0];
    let ma = RHO * std::f64::consts::PI * source.chord_m.powi(2) * source.span_m / 4.0;
    let pitch = ma * source.chord_m.powi(2) / 32.0;
    let hinge = HingedStrip::new(source, [0.0; 3], [0.0, 1.0, 0.0], 0).unwrap();
    for q in [-1.2_f64, 0.0, 0.47, 1.57] {
        let strip = hinge.at_angle(q).unwrap();
        close(strip.position_m[0], radius * q.cos(), 2e-13);
        close(strip.position_m[2], -radius * q.sin(), 2e-13);
        close(strip.normal[0], q.sin(), 2e-13);
        close(strip.normal[2], q.cos(), 2e-13);
        close(strip.control_gain, -radius, 2e-13);
        let rate = 1.7;
        let loads = assemble_hinged_strips(
            RHO,
            &[],
            std::slice::from_ref(&hinge),
            &[q],
            &[0.0; 6],
            &[rate],
        )
        .unwrap();
        close(loads.m_added_cc[0][0], ma * radius * radius + pitch, 2e-13);
        // The fluid momentum is tangent to the circular centroid path.
        // Its derivative requires a centripetal force; the body receives
        // the opposite, outward inertial reaction even with nu=0.
        close(
            loads.q_added_bias[0],
            ma * radius * rate * rate * q.cos(),
            3e-13,
        );
        close(
            loads.q_added_bias[2],
            -ma * radius * rate * rate * q.sin(),
            3e-13,
        );
        for i in [1, 3, 4, 5, 6] {
            close(loads.q_added_bias[i], 0.0, 3e-13);
        }
    }
    let at_zero = hinge.at_angle(0.0).unwrap();
    let fixed = assemble_analytic_strip(RHO, &[at_zero], 1, &[0.0; 6], &[0.0]).unwrap();
    let moving = assemble_hinged_strips(RHO, &[], &[hinge], &[0.0], &[0.0; 6], &[0.0]).unwrap();
    assert_eq!(fixed, moving);
}

#[test]
fn released_hinge_has_kinetic_gradient_even_when_its_rate_is_zero() {
    let source = plate();
    let ma = RHO * std::f64::consts::PI * source.chord_m.powi(2) * source.span_m / 4.0;
    let hinge = HingedStrip::new(source, [0.0; 3], [0.0, 1.0, 0.0], 0).unwrap();
    // Keep the flow away from the normal-aligned stationary point of T(q).
    let angle = 0.8_f64;
    let nu = [0.6, 0.0, 1.2, 0.0, 0.0, 0.0];
    let normal_speed = nu[0] * angle.sin() + nu[2] * angle.cos();
    let tangent_speed = nu[0] * angle.cos() - nu[2] * angle.sin();
    for rate in [0.0, -0.7, 1.2] {
        let loads = assemble_hinged_strips(
            RHO,
            &[],
            std::slice::from_ref(&hinge),
            &[angle],
            &nu,
            &[rate],
        )
        .unwrap();
        close(
            loads.q_added_bias[6],
            ma * normal_speed * tangent_speed,
            3e-13,
        );
        assert!(loads.q_added_bias[6].abs() > 0.1);
    }
}

#[test]
fn analytic_configuration_bias_matches_lagrange_equations_and_joint_energy() {
    // The independent oracle differentiates the entire mass matrix in tests
    // only, then uses Euler-Lagrange plus body momentum transport. It checks
    // every force component, not only the scalar power identity.
    let mut second = plate();
    second.position_m = [-0.5, 0.4, 0.3];
    second.chord_m = 0.7;
    let hinges = [
        HingedStrip::new(plate(), [0.1, 0.3, -0.2], [0.0, 1.0, 0.0], 0).unwrap(),
        HingedStrip::new(second, [-0.2, 0.1, 0.0], [0.36, 0.48, 0.8], 1).unwrap(),
    ];
    let angles = [0.35, -0.51];
    let rates = [0.7, -0.2];
    let nu = [2.0, -0.3, 0.5, 0.4, -0.2, 0.6];
    let u: Vec<_> = nu.into_iter().chain(rates).collect();
    let mut loads = assemble_hinged_strips(RHO, &[], &hinges, &angles, &nu, &rates).unwrap();
    let m = matrix(&loads);
    let mut derivative = vec![vec![vec![0.0; 8]; 8]; 2];
    for k in 0..2 {
        let h = 1e-5;
        let mut plus = angles;
        let mut minus = angles;
        plus[k] += h;
        minus[k] -= h;
        let p = matrix(&assemble_hinged_strips(RHO, &[], &hinges, &plus, &nu, &rates).unwrap());
        let n = matrix(&assemble_hinged_strips(RHO, &[], &hinges, &minus, &nu, &rates).unwrap());
        for i in 0..8 {
            for j in 0..8 {
                derivative[k][i][j] = (p[i][j] - n[i][j]) / (2.0 * h);
            }
        }
    }
    let mut mdot = vec![vec![0.0; 8]; 8];
    for i in 0..8 {
        for j in 0..8 {
            mdot[i][j] = (0..2).map(|k| rates[k] * derivative[k][i][j]).sum();
        }
    }
    let mdotu = apply(&mdot, &u);
    let momentum = apply(&m, &u);
    let v = [nu[0], nu[1], nu[2]];
    let w = [nu[3], nu[4], nu[5]];
    let linear = [momentum[0], momentum[1], momentum[2]];
    let angular = [momentum[3], momentum[4], momentum[5]];
    let force = cross(w, linear);
    let moment = cross(w, angular);
    let transport = cross(v, linear);
    let mut expected: Vec<_> = mdotu.iter().map(|x| -x).collect();
    for i in 0..3 {
        expected[i] -= force[i];
        expected[3 + i] -= moment[i] + transport[i];
    }
    for k in 0..2 {
        expected[6 + k] += 0.5
            * u.iter()
                .zip(apply(&derivative[k], &u))
                .map(|(x, y)| x * y)
                .sum::<f64>();
    }
    for i in 0..8 {
        close(loads.q_added_bias[i], expected[i], 2e-8);
    }
    let metric_power = 0.5 * u.iter().zip(&mdotu).map(|(x, y)| x * y).sum::<f64>();
    close(
        u.iter()
            .zip(&loads.q_added_bias)
            .map(|(x, y)| x * y)
            .sum::<f64>(),
        -metric_power,
        2e-8,
    );

    // The actual shared joint solve must close instantaneous mechanical
    // power with rigid Euler/Coriolis bias included exactly once.
    let body = RigidBody::new(3.0, [1.2, 1.7, 1.9]).unwrap();
    let rigid_force = cross(w, v.map(|x| 3.0 * x));
    let rigid_moment = cross(w, [1.2 * w[0], 1.7 * w[1], 1.9 * w[2]]);
    for i in 0..3 {
        loads.q_rigid_nonaccel[i] = -rigid_force[i];
        loads.q_rigid_nonaccel[3 + i] = -rigid_moment[i];
    }
    let (a, b) = solve_joint(&body, &[0.2, 0.3], &loads).unwrap();
    let acceleration: Vec<_> = a.into_iter().chain(b).collect();
    let mut total = m;
    for (i, inertia) in [3.0, 3.0, 3.0, 1.2, 1.7, 1.9, 0.2, 0.3]
        .into_iter()
        .enumerate()
    {
        total[i][i] += inertia;
    }
    let acceleration_power = u
        .iter()
        .zip(apply(&total, &acceleration))
        .map(|(x, y)| x * y)
        .sum::<f64>();
    close(acceleration_power + metric_power, 0.0, 3e-8);
}

#[test]
fn subdividing_a_moving_strip_preserves_full_inertia_and_configuration_bias() {
    // All children share the same physical hinge. Their centroid gains are
    // independently derived, including rotation about an oblique axis.
    let source = plate();
    let axis = [0.36, 0.48, 0.8];
    let pivot = [0.1, -0.2, 0.3];
    let parent = HingedStrip::new(source.clone(), pivot, axis, 0).unwrap();
    let nu = [0.4, -0.2, 1.3, 0.6, 0.2, -0.3];
    let expected = assemble_hinged_strips(RHO, &[], &[parent], &[0.8], &nu, &[1.4]).unwrap();
    for count in [2, 5, 13] {
        let children: Vec<_> = (0..count)
            .map(|i| {
                let mut child = source.clone();
                let offset = source.span_m * ((i as f64 + 0.5) / count as f64 - 0.5);
                child.position_m[1] += offset;
                child.span_m /= count as f64;
                HingedStrip::new(child, pivot, axis, 0).unwrap()
            })
            .collect();
        let actual = assemble_hinged_strips(RHO, &[], &children, &[0.8], &nu, &[1.4]).unwrap();
        for (row_a, row_b) in matrix(&actual).iter().zip(matrix(&expected)) {
            for (a, b) in row_a.iter().zip(row_b) {
                close(*a, b, 3e-12);
            }
        }
        for (a, b) in actual.q_added_bias.iter().zip(&expected.q_added_bias) {
            close(*a, *b, 3e-12);
        }
    }
}

#[test]
fn invalid_hinges_and_combined_resource_caps_refuse_before_assembly() {
    for axis in [[0.0; 3], [0.0, 2.0, 0.0], [f64::NAN, 0.0, 1.0]] {
        assert_eq!(
            HingedStrip::new(plate(), [0.0; 3], axis, 0)
                .unwrap_err()
                .code,
            "hinge-invalid"
        );
    }
    let last =
        HingedStrip::new(plate(), [0.0; 3], [0.0, 1.0, 0.0], MAX_CONTROL_COORDS - 1).unwrap();
    assert!(
        assemble_hinged_strips(
            RHO,
            &[],
            &[last],
            &[0.0; MAX_CONTROL_COORDS],
            &[0.0; 6],
            &[0.0; MAX_CONTROL_COORDS]
        )
        .is_ok()
    );
    assert_eq!(
        HingedStrip::new(plate(), [0.0; 3], [0.0, 1.0, 0.0], MAX_CONTROL_COORDS)
            .unwrap_err()
            .code,
        "control-coords-outside-domain"
    );
    let hinge = HingedStrip::new(plate(), [0.0; 3], [0.0, 1.0, 0.0], 0).unwrap();
    assert_eq!(
        hinge.at_angle(f64::NAN).unwrap_err().code,
        "non-finite-input"
    );
    let fixed = vec![plate(); MAX_STRIPS];
    assert!(
        assemble_hinged_strips(
            RHO,
            &fixed[..MAX_STRIPS - 1],
            std::slice::from_ref(&hinge),
            &[0.0],
            &[0.0; 6],
            &[0.0]
        )
        .is_ok()
    );
    assert_eq!(
        assemble_hinged_strips(RHO, &fixed, &[hinge], &[0.0], &[0.0; 6], &[0.0])
            .unwrap_err()
            .code,
        "strip-count-exceeded"
    );
    let mut preassigned = plate();
    preassigned.control_gain = 0.2;
    assert_eq!(
        HingedStrip::new(preassigned, [0.0; 3], [0.0, 1.0, 0.0], 0)
            .unwrap_err()
            .code,
        "hinge-invalid"
    );
    assert_eq!(
        assemble_hinged_strips(
            RHO,
            &[],
            &[],
            &[0.0; MAX_CONTROL_COORDS + 1],
            &[0.0; 6],
            &[0.0; MAX_CONTROL_COORDS + 1]
        )
        .unwrap_err()
        .code,
        "control-coords-outside-domain"
    );
}
