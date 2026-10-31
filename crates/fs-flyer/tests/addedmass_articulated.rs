//! G0/G1/G3: the native articulated geometry drives the existing fluid law.
use fs_flyer::addedmass::articulated::{LinkedStrip, assemble_articulated_strips};
use fs_flyer::addedmass::hinged::{HingedStrip, assemble_hinged_strips};
use fs_flyer::addedmass::{
    AeroGeneralizedLoads, MAX_CONTROL_COORDS, MAX_STRIPS, Strip, assemble_analytic_strip,
};
use fs_ga::{Mat3, Se3, So3, So3Tangent, Twist, Vec3};
use fs_mbd::articulated::generalized::{GeneralizedKinematics, free_floating_terms};
use fs_mbd::articulated::{
    ArticulatedModel, BaseState, FreeFloatingBaseState, JointModel, Link, SpatialInertia,
    forward_kinematics,
};

fn vector(a: [f64; 3]) -> Vec3 {
    Vec3::new(a[0], a[1], a[2])
}
fn array(v: Vec3) -> [f64; 3] {
    [v.x, v.y, v.z]
}
fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance * (1.0 + expected.abs()),
        "{actual:e} != {expected:e}"
    );
}
fn inertia() -> SpatialInertia {
    SpatialInertia::new(
        2.0,
        Vec3::new(0.1, -0.03, 0.07),
        Mat3 {
            m: [0.4, 0.02, -0.01, 0.02, 0.5, 0.03, -0.01, 0.03, 0.6],
        },
    )
    .unwrap()
}
fn transform(r: [f64; 3], p: [f64; 3]) -> Se3 {
    Se3::from_parts(So3::exp(So3Tangent::new(vector(r))).unwrap(), vector(p)).unwrap()
}
fn strip(position: [f64; 3]) -> Strip {
    Strip {
        name: "plate",
        chord_m: 0.75,
        span_m: 1.2,
        position_m: position,
        normal: [0.6, 0.0, 0.8],
        span_axis: [0.0, 1.0, 0.0],
        control_coord: None,
        control_gain: 0.0,
        control_angular_gain: [0.0; 3],
    }
}
fn base(nu: &[f64; 6]) -> FreeFloatingBaseState {
    FreeFloatingBaseState::new(
        transform([0.2, -0.3, 0.1], [2.0, -3.0, 1.0]),
        Twist::new(vector([nu[3], nu[4], nu[5]]), vector([nu[0], nu[1], nu[2]])),
    )
}
fn kinematics(
    model: &ArticulatedModel,
    q: &[f64],
    v: &[f64],
    nu: &[f64; 6],
) -> GeneralizedKinematics {
    free_floating_terms(model, base(nu), q, v, || false)
        .unwrap()
        .kinematics
}
fn full(loads: &AeroGeneralizedLoads) -> Vec<Vec<f64>> {
    let nc = loads.m_added_cc.len();
    let mut m = vec![vec![0.0; 6 + nc]; 6 + nc];
    for i in 0..6 {
        for j in 0..6 {
            m[i][j] = loads.m_added_rr[i][j];
        }
    }
    for i in 0..nc {
        for j in 0..6 {
            m[6 + i][j] = loads.m_added_rc[i][j];
            m[j][6 + i] = loads.m_added_rc[i][j];
        }
        for j in 0..nc {
            m[6 + i][6 + j] = loads.m_added_cc[i][j];
        }
    }
    m
}
fn compare(a: &AeroGeneralizedLoads, b: &AeroGeneralizedLoads, tolerance: f64) {
    for (row, expected) in full(a).iter().zip(full(b)) {
        for (x, y) in row.iter().zip(expected) {
            close(*x, y, tolerance);
        }
    }
    for (x, y) in a.q_added_bias.iter().zip(&b.q_added_bias) {
        close(*x, *y, tolerance);
    }
}
fn kinetic(m: &[Vec<f64>], u: &[f64]) -> f64 {
    let mut value = 0.0;
    for i in 0..u.len() {
        for j in 0..u.len() {
            value += 0.5 * u[i] * m[i][j] * u[j];
        }
    }
    value
}

#[test]
fn fixed_root_and_offset_link_recover_the_original_body_strip_model() {
    let offset = transform([0.13, -0.21, 0.31], [0.4, -0.2, 0.6]);
    let model = ArticulatedModel::new(vec![
        Link::new("root", None, Se3::identity(), JointModel::FIXED, inertia()),
        Link::new(
            "fixed surface",
            Some(0),
            offset,
            JointModel::FIXED,
            inertia(),
        ),
    ])
    .unwrap();
    let a = strip([0.1, 0.3, -0.2]);
    let b = strip([-0.2, 0.15, 0.08]);
    let mut b_body = b.clone();
    b_body.position_m = array(offset.transform_point(vector(b.position_m)).unwrap());
    b_body.normal = array(offset.rotation().rotate(vector(b.normal)).unwrap());
    b_body.span_axis = array(offset.rotation().rotate(vector(b.span_axis)).unwrap());
    let linked = [
        LinkedStrip::new(0, a.clone()).unwrap(),
        LinkedStrip::new(1, b).unwrap(),
    ];
    assert_eq!(linked[1].link(), 1);
    let nu = [1.3, -0.7, 0.5, 0.4, -0.2, 0.3];
    let actual =
        assemble_articulated_strips(1.18, &linked, &kinematics(&model, &[], &[], &nu)).unwrap();
    let expected = assemble_analytic_strip(1.18, &[a, b_body], 0, &nu, &[]).unwrap();
    compare(&actual, &expected, 2e-13);
    assert_eq!(actual.q_rigid_nonaccel, [0.0; 6]);
}

#[test]
fn finite_single_hinge_recovers_all_existing_mass_and_bias_components() {
    let pivot = [0.3, -0.2, 0.1];
    let axis = [0.36, 0.48, 0.8];
    let model = ArticulatedModel::new(vec![
        Link::new("root", None, Se3::identity(), JointModel::FIXED, inertia()),
        Link::new(
            "canard",
            Some(0),
            transform([0.0; 3], pivot),
            JointModel::revolute(vector(axis), None).unwrap(),
            inertia(),
        ),
    ])
    .unwrap();
    let reference = strip([0.8, -0.1, 0.15]);
    let mut local = reference.clone();
    local.position_m = core::array::from_fn(|i| reference.position_m[i] - pivot[i]);
    let linked = [LinkedStrip::new(1, local).unwrap()];
    let hinged = [HingedStrip::new(reference, pivot, axis, 0).unwrap()];
    let nu = [1.2, -0.4, 0.7, 0.28, -0.43, 0.19];
    for angle in [0.0, 0.63, -1.1] {
        for rate in [0.0, -0.78] {
            let map = kinematics(&model, &[angle], &[rate], &nu);
            let actual = assemble_articulated_strips(1.18, &linked, &map).unwrap();
            let expected =
                assemble_hinged_strips(1.18, &[], &hinged, &[angle], &nu, &[rate]).unwrap();
            compare(&actual, &expected, 3e-13);
        }
    }
}

fn chain(mixed: bool) -> ArticulatedModel {
    let first = if mixed {
        JointModel::prismatic(Vec3::new(0.3, 0.4, 0.5), None).unwrap()
    } else {
        JointModel::revolute(Vec3::new(0.0, 1.0, 0.0), None).unwrap()
    };
    let second = if mixed {
        JointModel::helical(Vec3::new(0.36, 0.48, 0.8), 0.12, None).unwrap()
    } else {
        JointModel::revolute(Vec3::new(0.36, 0.48, 0.8), None).unwrap()
    };
    ArticulatedModel::new(vec![
        Link::new("root", None, Se3::identity(), JointModel::FIXED, inertia()),
        Link::new(
            "upstream",
            Some(0),
            transform([0.1, -0.2, 0.3], [0.35, -0.27, 0.18]),
            first,
            inertia(),
        ),
        Link::new(
            "downstream",
            Some(1),
            transform([-0.2, 0.15, 0.07], [0.55, 0.08, -0.13]),
            second,
            inertia(),
        ),
    ])
    .unwrap()
}

#[test]
fn serial_and_mixed_joints_match_lagrange_bias_and_independent_strip_energy() {
    let q = [0.47, -0.36];
    let nu = [2.4, -0.8, 1.4, 0.56, -0.86, 0.38];
    let geometry = [
        strip([0.2, -0.3, 0.1]),
        strip([0.3, 0.08, -0.11]),
        strip([0.45, -0.1, 0.2]),
    ];
    let linked: Vec<_> = geometry
        .iter()
        .enumerate()
        .map(|(i, s)| LinkedStrip::new(i, s.clone()).unwrap())
        .collect();
    let rho = 1.18;
    for mixed in [false, true] {
        let model = chain(mixed);
        for rates in [[0.8, -0.6], [0.0, 0.0]] {
            let map = kinematics(&model, &q, &rates, &nu);
            let loads = assemble_articulated_strips(rho, &linked, &map).unwrap();
            let m = full(&loads);
            let mut u = nu.to_vec();
            u.extend_from_slice(&rates);
            assert!(
                m[6][7].abs() > 1e-3,
                "cross-control inertia is required: {:?}",
                loads.m_added_cc
            );
            for i in 0..8 {
                for j in 0..8 {
                    assert_eq!(m[i][j], m[j][i]);
                }
            }

            // Independent local energy: integrate actual normal velocity at
            // two span Gauss points and retain the exact section pitch inertia.
            let b = base(&nu);
            let poses = forward_kinematics(
                &model,
                BaseState::prescribed(b.world_from_base, b.twist_body, Twist::zero()),
                &q,
                &rates,
            )
            .unwrap();
            let mut physical_energy = 0.0;
            for (i, s) in geometry.iter().enumerate() {
                let v = poses.body_twist[i];
                let ma = rho * std::f64::consts::PI * s.chord_m.powi(2) * s.span_m / 4.0;
                for sign in [-1.0, 1.0] {
                    let point = vector(s.position_m)
                        + vector(s.span_axis).scale(sign * s.span_m / (2.0 * 3.0_f64.sqrt()));
                    let speed = (v.linear + v.angular.cross(point)).dot(vector(s.normal));
                    physical_energy += 0.25 * ma * speed * speed;
                }
                physical_energy += 0.5 * ma * s.chord_m.powi(2) / 32.0
                    * v.angular.dot(vector(s.span_axis)).powi(2);
            }
            close(kinetic(&m, &u), physical_energy, 3e-13);

            // Test-only configuration differentiation: independent Lagrange
            // equation oracle for every body and control bias component.
            let eps = 2e-6;
            let mut partials = Vec::new();
            for k in 0..2 {
                let mut plus = q;
                plus[k] += eps;
                let mut minus = q;
                minus[k] -= eps;
                let mp = full(
                    &assemble_articulated_strips(
                        rho,
                        &linked,
                        &kinematics(&model, &plus, &rates, &nu),
                    )
                    .unwrap(),
                );
                let mm = full(
                    &assemble_articulated_strips(
                        rho,
                        &linked,
                        &kinematics(&model, &minus, &rates, &nu),
                    )
                    .unwrap(),
                );
                let mut derivative = vec![vec![0.0; 8]; 8];
                for i in 0..8 {
                    for j in 0..8 {
                        derivative[i][j] = (mp[i][j] - mm[i][j]) / (2.0 * eps);
                    }
                }
                partials.push(derivative);
            }
            let mut mdot = vec![vec![0.0; 8]; 8];
            for i in 0..8 {
                for j in 0..8 {
                    mdot[i][j] = rates[0] * partials[0][i][j] + rates[1] * partials[1][i][j];
                }
            }
            let momentum: Vec<f64> = m
                .iter()
                .map(|row| row.iter().zip(&u).map(|(a, b)| a * b).sum())
                .collect();
            let omega = vector([nu[3], nu[4], nu[5]]);
            let velocity = vector([nu[0], nu[1], nu[2]]);
            let p = vector([momentum[0], momentum[1], momentum[2]]);
            let h = vector([momentum[3], momentum[4], momentum[5]]);
            let transport = [
                array(omega.cross(p)),
                array(omega.cross(h) + velocity.cross(p)),
            ]
            .concat();
            for i in 0..8 {
                let mut expected = -mdot[i].iter().zip(&u).map(|(a, b)| a * b).sum::<f64>();
                if i < 6 {
                    expected -= transport[i];
                } else {
                    expected += kinetic(&partials[i - 6], &u);
                }
                close(loads.q_added_bias[i], expected, 4e-9);
            }
            let power: f64 = loads.q_added_bias.iter().zip(&u).map(|(a, b)| a * b).sum();
            close(power, -kinetic(&mdot, &u), 3e-9);
            if rates == [0.0; 2] {
                assert!(
                    loads.q_added_bias[6..].iter().any(|x| x.abs() > 0.02),
                    "released joints retain kinetic-gradient torque: mixed={mixed}, bias={:?}",
                    &loads.q_added_bias[6..]
                );
            }
        }
    }
}

#[test]
fn downstream_strip_subdivision_preserves_all_serial_inertia_and_bias() {
    let model = chain(false);
    let map = kinematics(
        &model,
        &[0.6, -0.7],
        &[0.9, -0.4],
        &[0.8, 0.2, -0.1, 0.3, -0.5, 0.4],
    );
    let whole = strip([0.45, -0.1, 0.2]);
    let expected =
        assemble_articulated_strips(1.18, &[LinkedStrip::new(2, whole.clone()).unwrap()], &map)
            .unwrap();
    for count in [2, 5, 13] {
        let pieces: Vec<_> = (0..count)
            .map(|i| {
                let mut s = whole.clone();
                s.span_m /= count as f64;
                let offset = ((i as f64 + 0.5) / count as f64 - 0.5) * whole.span_m;
                s.position_m =
                    array(vector(whole.position_m) + vector(whole.span_axis).scale(offset));
                LinkedStrip::new(2, s).unwrap()
            })
            .collect();
        compare(
            &assemble_articulated_strips(1.18, &pieces, &map).unwrap(),
            &expected,
            4e-13,
        );
    }
}

#[test]
fn linked_geometry_and_combined_assembly_caps_refuse_invalid_bindings() {
    let mut preset = strip([0.0; 3]);
    preset.control_gain = 0.1;
    assert!(LinkedStrip::new(0, preset).is_err());
    let mut preset = strip([0.0; 3]);
    preset.control_coord = Some(0);
    assert!(LinkedStrip::new(0, preset).is_err());
    let model = chain(false);
    let map = kinematics(&model, &[0.0; 2], &[0.0; 2], &[0.0; 6]);
    assert!(
        assemble_articulated_strips(1.18, &[LinkedStrip::new(3, strip([0.0; 3])).unwrap()], &map)
            .is_err()
    );
    for rho in [0.0, -1.0, f64::NAN] {
        assert!(assemble_articulated_strips(rho, &[], &map).is_err());
    }
    let strips = vec![LinkedStrip::new(0, strip([0.0; 3])).unwrap(); MAX_STRIPS];
    assert!(assemble_articulated_strips(1.18, &strips, &map).is_ok());
    let mut too_many = strips;
    too_many.push(LinkedStrip::new(1, strip([0.0; 3])).unwrap());
    assert_eq!(
        assemble_articulated_strips(1.18, &too_many, &map)
            .unwrap_err()
            .code,
        "strip-count-exceeded"
    );
    for nc in [MAX_CONTROL_COORDS, MAX_CONTROL_COORDS + 1] {
        let mut links = vec![Link::new(
            "root",
            None,
            Se3::identity(),
            JointModel::FIXED,
            inertia(),
        )];
        for i in 0..nc {
            links.push(Link::new(
                format!("joint {i}"),
                Some(i),
                Se3::identity(),
                JointModel::revolute(Vec3::new(0.0, 1.0, 0.0), None).unwrap(),
                inertia(),
            ));
        }
        let model = ArticulatedModel::new(links).unwrap();
        let map = kinematics(&model, &vec![0.0; nc], &vec![0.0; nc], &[0.0; 6]);
        let result = assemble_articulated_strips(1.18, &[], &map);
        if nc == MAX_CONTROL_COORDS {
            assert!(result.is_ok());
        } else {
            assert_eq!(result.unwrap_err().code, "control-coords-outside-domain");
        }
    }
}
