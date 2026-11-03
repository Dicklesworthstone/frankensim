//! Full circulatory-force receipts: signed lift, lateral force, frame
//! covariance and analytic section loading. These checks do not add a drag
//! model or extend the nonlinear body-xz section equation to vertical surfaces.
//! Repro: cargo test -p fs-wing --test force_vector_battery

use fs_wing::images::{CertifiedGround, solve_weissinger_ground};
use fs_wing::nonlinear::{
    InfluenceOperator, StripRegime, StripSpec, solve_nonlinear_with_incident,
};
use fs_wing::{Panel, SolveReport, SurfaceId, flat_surface, solve_weissinger_linear};

const RHO: f64 = 1.225;

fn close(actual: f64, expected: f64, tolerance: f64) {
    let scale = actual.abs().max(expected.abs()).max(1.0);
    assert!(
        actual.is_finite()
            && expected.is_finite()
            && (actual - expected).abs() <= tolerance * scale,
        "actual {actual:e}, expected {expected:e}, tolerance {tolerance:e}"
    );
}

fn vector_close(actual: [f64; 3], expected: [f64; 3], tolerance: f64) {
    for axis in 0..3 {
        close(actual[axis], expected[axis], tolerance);
    }
}

fn perpendicular(force: [f64; 3], velocity: [f64; 3]) {
    let power = force[0] * velocity[0] + force[1] * velocity[1] + force[2] * velocity[2];
    let force_norm = force.iter().map(|x| x * x).sum::<f64>().sqrt();
    let speed = velocity.iter().map(|x| x * x).sum::<f64>().sqrt();
    assert!(force_norm > 1.0, "fixture must carry a nonzero force");
    assert!(
        power.abs() <= 1.0e-12 * force_norm * speed,
        "circulatory force does work along uniform flow: F dot V = {power:e}"
    );
}

fn rotate_panels(panels: &[Panel], rotate: impl Fn([f64; 3]) -> [f64; 3]) -> Vec<Panel> {
    panels
        .iter()
        .map(|panel| Panel {
            a: rotate(panel.a),
            b: rotate(panel.b),
            ctrl: rotate(panel.ctrl),
            normal: rotate(panel.normal),
            ..*panel
        })
        .collect()
}

fn scalar_projection(report: &SolveReport) {
    assert_eq!(
        report.total_force_n[2].to_bits(),
        (-report.total_lift_n).to_bits(),
        "the legacy lift receipt must remain the exact negative body-z component"
    );
    close(
        report.surface_lift_n.iter().map(|(_, lift)| lift).sum(),
        report.total_lift_n,
        1.0e-13,
    );
}

#[test]
fn linear_signed_lift_rotates_into_vertical_fin_and_canted_wing_force() {
    let panels = flat_surface(SurfaceId::WingLower, 6.0, 1.0, 0.0, 0.0, 6, 2).unwrap();
    for incidence in [-0.8, 0.8] {
        let velocity = [12.0, 0.0, incidence];
        let reference = solve_weissinger_linear(&panels, velocity, RHO).unwrap();
        scalar_projection(&reference);
        perpendicular(reference.total_force_n, velocity);
        assert!(reference.total_lift_n * incidence > 0.0);

        // Proper rotations about x: an exact quarter turn produces a
        // vertical fin; the 3-4-5 rotation produces a canted lifting surface.
        for (sin, cos) in [(1.0, 0.0), (0.6, 0.8)] {
            let roll = |v: [f64; 3]| [v[0], cos * v[1] - sin * v[2], sin * v[1] + cos * v[2]];
            let mut rotated_panels = rotate_panels(&panels, &roll);
            if cos == 0.0 {
                for panel in &mut rotated_panels {
                    panel.surface = SurfaceId::Vertical;
                }
            }
            let rotated_velocity = roll(velocity);
            let rotated = solve_weissinger_linear(&rotated_panels, rotated_velocity, RHO).unwrap();
            for (actual, expected) in rotated.gamma.iter().zip(&reference.gamma) {
                close(*actual, *expected, 2.0e-10);
            }
            vector_close(
                rotated.total_force_n,
                roll(reference.total_force_n),
                2.0e-10,
            );
            scalar_projection(&rotated);
            perpendicular(rotated.total_force_n, rotated_velocity);
            assert!(rotated.total_force_n[1] * incidence > 0.0);
            if cos == 0.0 {
                close(rotated.total_lift_n, 0.0, 1.0e-13);
                close(rotated.total_force_n[1], reference.total_lift_n, 2.0e-10);
            }
        }
    }
}

#[test]
fn ground_force_rotates_with_the_wing_while_preserving_the_plane() {
    let flat = flat_surface(SurfaceId::WingLower, 6.0, 1.0, 0.0, 0.0, 6, 2).unwrap();
    let roll = |v: [f64; 3]| [v[0], 0.8 * v[1] - 0.6 * v[2], 0.6 * v[1] + 0.8 * v[2]];
    let panels = rotate_panels(&flat, &roll);
    let velocity = roll([12.0, 0.0, 0.8]);
    // The canted span lies within z = +/-1.8 m, wholly above this plane.
    let ground = CertifiedGround {
        z_m: 3.0,
        certificate_slope: 0.0,
        certificate_rms_m: 0.0,
    };
    let reference = solve_weissinger_ground(&panels, velocity, RHO, &ground).unwrap();
    assert!(reference.total_force_n[1].abs() > 1.0);
    assert!(reference.total_force_n[2].abs() > 1.0);
    scalar_projection(&reference);
    perpendicular(reference.total_force_n, velocity);

    // A yaw rotation leaves z and the certified horizontal plane unchanged.
    // Rotating only the wing about x would change its relation to the ground.
    let yaw = |v: [f64; 3]| [-v[1], v[0], v[2]];
    let rotated_panels = rotate_panels(&panels, &yaw);
    let rotated_velocity = yaw(velocity);
    let rotated = solve_weissinger_ground(&rotated_panels, rotated_velocity, RHO, &ground).unwrap();
    for (actual, expected) in rotated.gamma.iter().zip(&reference.gamma) {
        close(*actual, *expected, 2.0e-10);
    }
    vector_close(rotated.total_force_n, yaw(reference.total_force_n), 2.0e-10);
    scalar_projection(&rotated);
    perpendicular(rotated.total_force_n, rotated_velocity);
}

#[test]
fn nonlinear_force_matches_signed_section_load_with_local_flow_and_wash() {
    let span = 2.0;
    let chord = 1.0;
    let panels = flat_surface(SurfaceId::WingLower, span, chord, 0.0, 0.0, 1, 1).unwrap();
    let strips = vec![StripSpec {
        panel_indices: vec![0],
        chord_m: chord,
        twist_rad: 0.0,
    }];
    let reference_velocity = [12.0, 0.0, 0.8];
    let operator = InfluenceOperator::build(&panels, reference_velocity, RHO).unwrap();
    let incident = [[14.0, 2.0, 1.0]];
    let wash = [3.0];
    let local: [f64; 3] = [17.0, 2.0, 1.0];
    let section_speed = (local[0] * local[0] + local[2] * local[2]).sqrt();
    let dynamic_pressure = 0.5 * RHO * section_speed * section_speed;

    // A constant section coefficient makes q*S*cl an independent analytic
    // oracle despite the solved trailing-sheet induction. Spanwise velocity
    // is present but cannot inflate the section speed or its lift.
    for cl in [-0.8, 0.8] {
        let closure = |_: usize, _: f64| (cl, StripRegime::Attached);
        let report = solve_nonlinear_with_incident(
            &operator,
            &panels,
            &strips,
            reference_velocity,
            RHO,
            &closure,
            &incident,
            Some(&wash),
            None,
        )
        .unwrap();
        let signed_lift = dynamic_pressure * span * chord * cl;
        let expected = [
            signed_lift * local[2] / section_speed,
            0.0,
            -signed_lift * local[0] / section_speed,
        ];
        vector_close(report.total_force_n, expected, 2.0e-9);
        perpendicular(report.total_force_n, local);
        assert_eq!(
            report.total_force_n[2].to_bits(),
            (-report.total_lift_n).to_bits()
        );
        assert!(report.total_lift_n * cl > 0.0);
    }
}
