//! G1/G3/G4: instantaneous elliptic loss, full mobility and retained history.
use super::*;
use fs_tribo::tangential::elliptic::impulse::*;
use fs_tribo::tangential::transport::{ObjectiveTangentialLaw, TangentialHistoryRotation};

fn controls() -> EllipticImpulseControls {
    EllipticImpulseControls {
        maximum_iterations: 128,
        maximum_impulse_ns: 100.0,
        impulse_absolute_tolerance_ns: 1e-11,
        impulse_relative_tolerance: 1e-12,
        velocity_tolerance_m_s: 1e-10,
        energy_tolerance_j: 1e-10,
    }
}
fn impact(mu: [f64; 2]) -> EllipticCoulomb {
    EllipticCoulomb::new(mu, mu, [1.0, 0.0]).unwrap()
}
fn request(velocity: [f64; 2], mobility: [[f64; 2]; 2]) -> EllipticImpulseInput {
    EllipticImpulseInput {
        normal_impulse_ns: 2.0,
        velocity_before_m_s: velocity,
        mobility_per_kg: mobility,
    }
}
fn matvec(m: [[f64; 2]; 2], v: [f64; 2]) -> [f64; 2] {
    [
        m[0][0] * v[0] + m[0][1] * v[1],
        m[1][0] * v[0] + m[1][1] * v[1],
    ]
}
fn dot(a: [f64; 2], b: [f64; 2]) -> f64 {
    a[0] * b[0] + a[1] * b[1]
}
fn kinetic(w: [[f64; 2]; 2], v: [f64; 2]) -> f64 {
    let determinant = w[0][0] * w[1][1] - w[0][1] * w[1][0];
    0.5 * (w[1][1] * v[0] * v[0] - 2.0 * w[0][1] * v[0] * v[1] + w[0][0] * v[1] * v[1])
        / determinant
}
fn audit(input: EllipticImpulseInput, report: EllipticImpulseReport, tolerance: f64) {
    let change = matvec(input.mobility_per_kg, report.impulse_ns);
    close2(
        report.velocity_after_m_s,
        [
            input.velocity_before_m_s[0] + change[0],
            input.velocity_before_m_s[1] + change[1],
        ],
        tolerance,
    );
    let delta = kinetic(input.mobility_per_kg, report.velocity_after_m_s)
        - kinetic(input.mobility_per_kg, input.velocity_before_m_s);
    close(report.kinetic_energy_change_j, delta, tolerance);
    close(delta + report.dissipated_j, 0.0, tolerance);
    assert!(report.dissipated_j >= 0.0);
    close(report.energy_residual_j, 0.0, tolerance);
}

#[test]
fn isotropic_and_principal_axis_impulses_match_independent_analytic_arrest_and_loss() {
    let w = [[1.75, 0.0], [0.0, 1.75]];
    for (mu, expected_impulse, expected_after, loss) in
        [(0.2, -4.0 / 7.0, 0.0, 2.0 / 7.0), (0.05, -0.2, 0.65, 0.165)]
    {
        let law = impact([mu; 2]);
        let input = EllipticImpulseInput {
            normal_impulse_ns: 4.0,
            ..request([1.0, 0.0], w)
        };
        let report = law
            .solve_impulse(&law.initial_state(), &input, controls(), &mut || false)
            .unwrap();
        close2(report.impulse_ns, [expected_impulse, 0.0], 2e-15);
        close2(report.velocity_after_m_s, [expected_after, 0.0], 2e-15);
        close(report.dissipated_j, loss, 2e-15);
        assert_eq!(report.iterations, 0);
        assert_eq!(
            report.regime,
            if expected_after == 0.0 {
                FrictionRegime::Sticking
            } else {
                FrictionRegime::Sliding
            }
        );
        audit(input, report, 3e-15);
    }
    let law = impact([0.4, 0.1]);
    for (v, p, end, loss) in [
        ([3.0, 0.0], [-0.8, 0.0], [1.4, 0.0], 1.76),
        ([0.0, 3.0], [0.0, -0.2], [0.0, 2.6], 0.56),
    ] {
        let input = request(v, [[2.0, 0.0], [0.0, 2.0]]);
        let report = law
            .solve_impulse(&law.initial_state(), &input, controls(), &mut || false)
            .unwrap();
        close2(report.impulse_ns, p, 3e-15);
        close2(report.velocity_after_m_s, end, 3e-15);
        close(report.dissipated_j, loss, 3e-15);
        audit(input, report, 4e-15);
    }
}

#[test]
fn off_axis_full_mobility_satisfies_independent_ellipse_kkt_and_energy_minimum() {
    let law = impact([0.6, 0.25]);
    let w = [[1.8, 0.35], [0.35, 0.9]];
    let endpoint = [0.7, -0.4];
    let caps: [f64; 2] = [1.2, 0.5];
    // Choose a final slip and construct the exact associated support point.
    let support = ((caps[0] * endpoint[0]).powi(2) + (caps[1] * endpoint[1]).powi(2)).sqrt();
    let reaction = [
        caps[0].powi(2) * endpoint[0] / support,
        caps[1].powi(2) * endpoint[1] / support,
    ];
    let recoil = matvec(w, reaction);
    let before = [endpoint[0] + recoil[0], endpoint[1] + recoil[1]];
    let input = request(before, w);
    let report = law
        .solve_impulse(&law.initial_state(), &input, controls(), &mut || false)
        .unwrap();
    close2(report.impulse_ns, reaction.map(|x| -x), 5e-11);
    close2(report.velocity_after_m_s, endpoint, 8e-11);
    assert!(report.iterations > 0 && report.iterations <= 128);
    assert!(report.impulse_ns[0] * endpoint[1] - report.impulse_ns[1] * endpoint[0] > 0.1);
    let objective = |p: [f64; 2]| dot(p, before) + 0.5 * dot(p, matvec(w, p));
    let accepted = objective(report.impulse_ns);
    for k in 0..128 {
        let angle = std::f64::consts::TAU * k as f64 / 128.0;
        let trial = [caps[0] * angle.cos(), caps[1] * angle.sin()];
        assert!(objective(trial) >= accepted - 2e-11);
    }
    assert!(report.capacity_residual_ns <= 2e-11);
    assert!(report.stationarity_residual_ns <= 2e-11);
    audit(input, report, 2e-10);

    let reaction = [0.2, -0.1];
    let input = request(matvec(w, reaction), w);
    let stopped = law
        .solve_impulse(&law.initial_state(), &input, controls(), &mut || false)
        .unwrap();
    assert_eq!(stopped.regime, FrictionRegime::Sticking);
    close2(stopped.impulse_ns, reaction.map(|x| -x), 3e-15);
    close2(stopped.velocity_after_m_s, [0.0; 2], 3e-15);
    close(
        stopped.dissipated_j,
        0.5 * dot(reaction, input.velocity_before_m_s),
        3e-15,
    );
    audit(input, stopped, 4e-15);
}

#[test]
fn retained_texture_rotation_and_mass_velocity_scaling_preserve_physical_impulse() {
    let law = impact([0.6, 0.25]);
    let state = law.initial_state();
    let w = [[1.8, 0.35], [0.35, 0.9]];
    let input = request([3.0, -1.0], w);
    let baseline = law
        .solve_impulse(&state, &input, controls(), &mut || false)
        .unwrap();
    let frame = TangentFrame::new([0.0, 0.0, 1.0], [1.0, 0.0, 0.0]).unwrap();
    for angle in [0.43_f64, 0.91] {
        let rotation = TangentialHistoryRotation::between(frame, frame, angle).unwrap();
        let rotated = law.transport_state(&state, rotation).unwrap();
        let [c, s] = rotation.cosine_sine();
        let [a, b, d] = [w[0][0], w[0][1], w[1][1]];
        let cross = c * s * (a - d) + (c * c - s * s) * b;
        let matrix = [
            [c * c * a - 2.0 * c * s * b + s * s * d, cross],
            [cross, s * s * a + 2.0 * c * s * b + c * c * d],
        ];
        let request = request(rotation.apply(input.velocity_before_m_s).unwrap(), matrix);
        let report = law
            .solve_impulse(&rotated, &request, controls(), &mut || false)
            .unwrap();
        close2(
            report.impulse_ns,
            rotation.apply(baseline.impulse_ns).unwrap(),
            6e-11,
        );
        close2(
            report.velocity_after_m_s,
            rotation.apply(baseline.velocity_after_m_s).unwrap(),
            1e-10,
        );
        close(report.dissipated_j, baseline.dissipated_j, 2e-10);
        audit(request, report, 3e-10);
    }
    for (mass, speed) in [(1e-6, 1e3), (1e6, 1e-3)] {
        let factor = mass * speed;
        let energy = mass * speed * speed;
        let scaled = EllipticImpulseInput {
            normal_impulse_ns: input.normal_impulse_ns * factor,
            velocity_before_m_s: input.velocity_before_m_s.map(|v| v * speed),
            mobility_per_kg: w.map(|row| row.map(|x| x / mass)),
        };
        let mut tolerance = controls();
        tolerance.maximum_impulse_ns *= factor;
        tolerance.impulse_absolute_tolerance_ns *= factor;
        tolerance.velocity_tolerance_m_s *= speed;
        tolerance.energy_tolerance_j *= energy;
        let report = law
            .solve_impulse(&state, &scaled, tolerance, &mut || false)
            .unwrap();
        close2(
            report.impulse_ns.map(|p| p / factor),
            baseline.impulse_ns,
            4e-11,
        );
        close2(
            report.velocity_after_m_s.map(|v| v / speed),
            baseline.velocity_after_m_s,
            8e-11,
        );
        close(report.dissipated_j / energy, baseline.dissipated_j, 2e-10);
        audit(scaled, report, 3e-10 * energy);
    }
}

#[test]
fn zero_axes_and_normal_impulse_preserve_free_motion_and_strict_mobility_limits() {
    let law = impact([0.0, 0.4]);
    let w = [[2.0, 0.5], [0.5, 1.0]];
    for (v, impulse, endpoint, loss) in [
        ([3.0, 0.2], [0.0, -0.2], [2.9, 0.0], 0.02),
        ([3.0, 2.0], [0.0, -0.8], [2.6, 1.2], 1.28),
    ] {
        let input = request(v, w);
        let report = law
            .solve_impulse(&law.initial_state(), &input, controls(), &mut || false)
            .unwrap();
        close2(report.impulse_ns, impulse, 3e-15);
        close2(report.velocity_after_m_s, endpoint, 3e-15);
        close(report.dissipated_j, loss, 3e-15);
        assert_eq!(report.iterations, 0);
        audit(input, report, 3e-15);
    }
    for (law, normal, velocity) in [
        (impact([0.0; 2]), 2.0, [3.0, -1.0]),
        (impact([0.6, 0.25]), 0.0, [3.0, -1.0]),
        (impact([0.6, 0.25]), 2.0, [0.0; 2]),
    ] {
        let input = EllipticImpulseInput {
            normal_impulse_ns: normal,
            ..request(velocity, w)
        };
        let report = law
            .solve_impulse(&law.initial_state(), &input, controls(), &mut || false)
            .unwrap();
        assert_eq!(report.impulse_ns, [0.0; 2]);
        assert_eq!(report.velocity_after_m_s, velocity);
        assert_eq!(report.dissipated_j, 0.0);
    }
    for bad in [
        [[1.0, 1.0], [1.0, 1.0]],
        [[1.0, 0.2], [0.3, 1.0]],
        [[-1.0, 0.0], [0.0, 1.0]],
    ] {
        assert!(
            law.solve_impulse(
                &law.initial_state(),
                &request([1.0; 2], bad),
                controls(),
                &mut || false
            )
            .is_err()
        );
    }
    let mut limited = controls();
    limited.maximum_impulse_ns = 0.1;
    assert!(
        law.solve_impulse(
            &law.initial_state(),
            &request([3.0, 2.0], w),
            limited,
            &mut || false
        )
        .is_err()
    );
    let overflow = impact([2.0; 2]);
    let huge = EllipticImpulseInput {
        normal_impulse_ns: f64::MAX,
        ..request([1.0; 2], w)
    };
    assert!(
        overflow
            .solve_impulse(&overflow.initial_state(), &huge, controls(), &mut || false)
            .is_err()
    );
}

#[test]
fn impact_never_advances_original_history_and_refusals_cancellation_retry_are_atomic() {
    let law = impact([0.6, 0.25]);
    let history = law
        .advance(&law.initial_state(), &super::input([0.04, 0.03], 0.0))
        .unwrap()
        .state;
    assert_ne!(history.slip_m(), [0.0; 2]);
    assert_eq!(history.contact_temperature_k(), Some(300.0));
    let original = history.clone();
    let input = request([3.0, -1.0], [[1.8, 0.35], [0.35, 0.9]]);
    let first = law
        .solve_impulse(&history, &input, controls(), &mut || false)
        .unwrap();
    assert_eq!(history, original);
    for cancellation_poll in [1, 2, 3] {
        let mut polls = 0;
        let result = law.solve_impulse(&history, &input, controls(), &mut || {
            polls += 1;
            polls == cancellation_poll
        });
        assert!(matches!(result, Err(EllipticImpulseError::Cancelled)));
        assert_eq!(history, original);
    }
    let mut limited = controls();
    limited.maximum_iterations = 1;
    assert!(
        law.solve_impulse(&history, &input, limited, &mut || false)
            .is_err()
    );
    let ambiguous = EllipticCoulomb::new([0.8, 0.4], [0.6, 0.25], [1.0, 0.0]).unwrap();
    assert!(matches!(
        ambiguous.solve_impulse(&ambiguous.initial_state(), &input, controls(), &mut || {
            false
        }),
        Err(EllipticImpulseError::UnequalCoefficients)
    ));
    assert!(matches!(
        impact([0.2; 2]).solve_impulse(&history, &input, controls(), &mut || false),
        Err(EllipticImpulseError::Constitutive(_))
    ));
    for invalid in [f64::NAN, -1.0, f64::INFINITY] {
        let bad = EllipticImpulseInput {
            normal_impulse_ns: invalid,
            ..input
        };
        assert!(
            law.solve_impulse(&history, &bad, controls(), &mut || false)
                .is_err()
        );
    }
    for _ in 0..2 {
        assert_eq!(
            law.solve_impulse(&history, &input, controls(), &mut || false)
                .unwrap(),
            first
        );
        assert_eq!(history, original);
    }
}
