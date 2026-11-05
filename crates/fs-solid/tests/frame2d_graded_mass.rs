//! Original graded beam density in the native mass, modal, impulse and time owners.

use fs_solid::frame2d::{FrameInertiaProfile, FrameMassFormulation};
use fs_solid::{Frame2d, Geometry, Rayleigh, SectionModel};

const LENGTH: f64 = 3.0;
const EA: f64 = 120.0;
const UNDAMPED: Rayleigh = Rayleigh { a0: 0.0, a1: 0.0 };

fn section() -> SectionModel {
    SectionModel::Elastic { ea: EA, ei: 40.0 }
}
fn profile() -> FrameInertiaProfile {
    FrameInertiaProfile::new([2.0, 6.0], [0.2, 0.8]).unwrap()
}
fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(a, b)| a * b).sum()
}
fn near(a: f64, b: f64, scale: f64) {
    assert!(
        (a - b).abs() < 3e-12 * scale,
        "{a:e} versus {b:e}, scale={scale:e}"
    );
}
fn build(form: FrameMassFormulation, inertia: FrameInertiaProfile) -> Frame2d {
    let mut frame = Frame2d::with_mass_formulation(Geometry::Linear, form);
    frame.add_node(0.0, 0.0);
    frame.add_node(LENGTH, 0.0);
    frame
        .add_element_with_inertia_profile(0, 1, 5, inertia, &section)
        .unwrap();
    frame
}
fn axial(inertia: FrameInertiaProfile) -> Frame2d {
    let mut frame = build(FrameMassFormulation::ConsistentReference, inertia);
    frame.fix(0, [true; 3]);
    frame.fix(1, [false, true, true]);
    frame
}

#[test]
fn full_rigid_mass_centroid_and_angular_inertia_follow_the_density_profile() {
    let start = [0.7, -0.4];
    // Independent integrals of mu(s)=2+4s and J(s)=0.2+0.6s, x=L*s.
    let mass = 12.0;
    let first = 21.0;
    let second = 45.0;
    let spin = 1.5;
    for angle in [0.0_f64, 0.7, 1.9] {
        for reversed in [false, true] {
            let tangent = [angle.cos(), angle.sin()];
            let end = [
                start[0] + LENGTH * tangent[0],
                start[1] + LENGTH * tangent[1],
            ];
            let mut frame = Frame2d::with_mass_formulation(
                Geometry::Linear,
                FrameMassFormulation::ConsistentReference,
            );
            frame.add_node(start[0], start[1]);
            frame.add_node(end[0], end[1]);
            let (i, j, inertia) = if reversed {
                (
                    1,
                    0,
                    FrameInertiaProfile::new([6.0, 2.0], [0.8, 0.2]).unwrap(),
                )
            } else {
                (0, 1, profile())
            };
            frame
                .add_element_with_inertia_profile(i, j, 5, inertia, &section)
                .unwrap();
            assert_eq!(frame.reference_member_inertia(0).unwrap(), inertia);
            let center_moment = [
                mass * start[0] + first * tangent[0],
                mass * start[1] + first * tangent[1],
            ];
            let angular = mass * (start[0] * start[0] + start[1] * start[1])
                + 2.0 * first * (start[0] * tangent[0] + start[1] * tangent[1])
                + second
                + spin;
            let expected = [
                [mass, 0.0, -center_moment[1]],
                [0.0, mass, center_moment[0]],
                [-center_moment[1], center_moment[0], angular],
            ];
            let rigid = [
                [1.0, 0.0, 0.0, 1.0, 0.0, 0.0],
                [0.0, 1.0, 0.0, 0.0, 1.0, 0.0],
                [-start[1], start[0], 1.0, -end[1], end[0], 1.0],
            ];
            for a in 0..3 {
                for b in 0..3 {
                    near(
                        dot(&rigid[a], &frame.inertial_forces(&rigid[b]).unwrap()),
                        expected[a][b],
                        (expected[a][a] * expected[b][b]).sqrt(),
                    );
                }
            }
        }
    }
}

fn polynomial_work(coefficients: &[f64], moments: [f64; 2]) -> f64 {
    // Direct monomial integration of (m0 + (m1-m0)*s) * a(s)^2.
    let mut total = 0.0;
    for (i, &a) in coefficients.iter().enumerate() {
        for (j, &b) in coefficients.iter().enumerate() {
            total += a
                * b
                * (moments[0] / (i + j + 1) as f64
                    + (moments[1] - moments[0]) / (i + j + 2) as f64);
        }
    }
    total
}

#[test]
fn exact_density_velocity_correlation_matches_independent_polynomial_energy() {
    let frame = build(FrameMassFormulation::ConsistentReference, profile());
    for v in [
        [0.3, -0.7, 1.1, -0.2, 0.8, -0.9],
        [0.0, 0.0, 1.0, 0.0, 0.0, -1.0],
        [0.7, -1.2, 0.0, 0.7, -1.2, 0.0],
        [1.0, 0.2, -0.4, 0.3, -0.8, 0.7],
    ] {
        let axial = [v[0], v[3] - v[0]];
        let transverse = [
            v[1],
            LENGTH * v[2],
            -3.0 * v[1] - 2.0 * LENGTH * v[2] + 3.0 * v[4] - LENGTH * v[5],
            2.0 * v[1] + LENGTH * v[2] - 2.0 * v[4] + LENGTH * v[5],
        ];
        let slope = [
            transverse[1] / LENGTH,
            2.0 * transverse[2] / LENGTH,
            3.0 * transverse[3] / LENGTH,
        ];
        let exact = LENGTH
            * (polynomial_work(&axial, [2.0, 6.0])
                + polynomial_work(&transverse, [2.0, 6.0])
                + polynomial_work(&slope, [0.2, 0.8]));
        near(dot(&v, &frame.inertial_forces(&v).unwrap()), exact, exact);
    }
}

#[test]
fn equal_total_mass_profiles_change_the_actual_mode_and_impulse_response() {
    for (mu, modal_mass) in [([2.0, 6.0], 5.0), ([6.0, 2.0], 3.0)] {
        let mut frame = axial(FrameInertiaProfile::new(mu, [0.0; 2]).unwrap());
        let frequencies = frame.natural_frequencies().unwrap();
        assert_eq!(frequencies.len(), 1);
        // N_tip=x/L: M_ff = L*(mu_i + 3*mu_j)/12; K_ff=EA/L.
        near(
            frequencies[0].powi(2),
            EA / LENGTH / modal_mass,
            EA / LENGTH / modal_mass,
        );
        let event = frame
            .apply_impulse(&[0.0, 0.0, 0.0, 10.0, 0.0, 0.0], &[0.0; 6], UNDAMPED, 0.01)
            .unwrap();
        near(
            event.velocity_change[3],
            10.0 / modal_mass,
            10.0 / modal_mass,
        );
        // The fixed-end reaction retains M_if=L*(mu_i+mu_j)/12=2.
        near(
            event.support_impulse[0],
            20.0 / modal_mass,
            20.0 / modal_mass,
        );
        near(
            event.kinetic_energy_after_j,
            50.0 / modal_mass,
            50.0 / modal_mass,
        );
        near(
            event.impulse_work_j,
            event.kinetic_energy_after_j,
            event.kinetic_energy_after_j,
        );
    }
}

#[test]
fn loaded_newmark_uses_the_same_profile_for_weight_base_forcing_and_motion() {
    let mut reference = axial(profile());
    let gravity = reference
        .inertial_forces(&[-0.2, 0.0, 0.0, -0.2, 0.0, 0.0])
        .unwrap();
    // integral mu*N_tip dx=7 kg, so gravity supplies -1.4 N to the free DOF.
    near(gravity[3], -1.4, 1.4);
    reference.static_load(&gravity, 1).unwrap();
    let stiffness = EA / LENGTH;
    let equilibrium = -1.4 / stiffness;
    near(reference.displacements()[3], equilibrium, equilibrium.abs());
    let omega = (stiffness / 5.0).sqrt();
    let duration = 0.73 * 2.0 * std::f64::consts::PI / omega;
    let expected = equilibrium - 7.0 / stiffness * (1.0 - (omega * duration).cos());
    let mut previous = f64::INFINITY;
    for steps in [80, 160, 320] {
        let mut frame = reference.clone();
        let history = frame
            .newmark_with_load(
                &vec![1.0; steps + 1],
                duration / steps as f64,
                (1.0, 0.0),
                UNDAMPED,
                &gravity,
            )
            .unwrap();
        let error = (history.u[steps][3] - expected).abs();
        assert!(
            error < 0.30 * previous,
            "steps={steps}, error={error:e}, previous={previous:e}"
        );
        assert!(history.energy_balance_error() < 1e-8);
        assert_eq!(frame.reference_member_inertia(0).unwrap(), profile());
        previous = error;
    }
    assert!(previous < 1e-4 * expected.abs(), "error={previous:e}");
}

#[test]
fn lumped_profile_preserves_total_mass_first_moment_and_section_spin() {
    let frame = build(FrameMassFormulation::Lumped, profile());
    let expected = [5.0, 5.0, 0.6, 7.0, 7.0, 0.9];
    for i in 0..6 {
        let mut q = [0.0; 6];
        q[i] = 1.0;
        let force = frame.inertial_forces(&q).unwrap();
        for j in 0..6 {
            near(
                force[j],
                if i == j { expected[i] } else { 0.0 },
                expected[i],
            );
        }
    }
    let gravity = frame
        .inertial_forces(&[0.0, -1.0, 0.0, 0.0, -1.0, 0.0])
        .unwrap();
    near(gravity[1] + gravity[4], -12.0, 12.0);
    near(LENGTH * gravity[4], -21.0, 21.0);
}

#[test]
fn uniform_profiles_keep_the_native_constant_paths_bitwise() {
    let uniform = FrameInertiaProfile::new([4.0; 2], [0.5; 2]).unwrap();
    for formulation in [
        FrameMassFormulation::Lumped,
        FrameMassFormulation::ConsistentReference,
    ] {
        let mut legacy = Frame2d::with_mass_formulation(Geometry::Linear, formulation);
        legacy.add_node(0.0, 0.0);
        legacy.add_node(LENGTH, 0.0);
        legacy
            .add_element_with_rotary_inertia(0, 1, 5, 4.0, 0.5, &section)
            .unwrap();
        let mut graded = build(formulation, uniform);
        for i in 0..6 {
            let mut q = [0.0; 6];
            q[i] = 1.0;
            assert_eq!(
                legacy
                    .inertial_forces(&q)
                    .unwrap()
                    .iter()
                    .map(|v| v.to_bits())
                    .collect::<Vec<_>>(),
                graded
                    .inertial_forces(&q)
                    .unwrap()
                    .iter()
                    .map(|v| v.to_bits())
                    .collect::<Vec<_>>()
            );
        }
        for frame in [&mut legacy, &mut graded] {
            frame.fix(0, [true; 3]);
            frame.fix(1, [false, true, true]);
        }
        assert_eq!(
            legacy.natural_frequencies().unwrap(),
            graded.natural_frequencies().unwrap()
        );
        assert_eq!(
            legacy
                .newmark(
                    &[0.0, 0.5, 1.0, -0.3, -0.7, 0.0],
                    0.02,
                    (1.0, 0.0),
                    UNDAMPED
                )
                .unwrap(),
            graded
                .newmark(
                    &[0.0, 0.5, 1.0, -0.3, -0.7, 0.0],
                    0.02,
                    (1.0, 0.0),
                    UNDAMPED
                )
                .unwrap()
        );
        assert_eq!(graded.reference_member_inertia(0).unwrap(), uniform);
    }
}

#[test]
fn bad_profile_and_unrepresentable_inventory_refuse_without_frame_changes() {
    for (mass, spin) in [
        ([-1.0, 2.0], [0.0; 2]),
        ([1.0, 2.0], [0.0, -1.0]),
        ([f64::NAN, 2.0], [0.0; 2]),
        ([1.0, 2.0], [0.0, f64::INFINITY]),
        ([0.0, 2.0], [1.0, 0.0]),
    ] {
        assert!(FrameInertiaProfile::new(mass, spin).is_err());
    }
    let mut frame = axial(profile());
    let probe = [0.3, -0.7, 1.1, -0.2, 0.8, -0.9];
    let original = frame.inertial_forces(&probe).unwrap();
    assert!(
        frame
            .add_element_with_inertia_profile(0, 1, 2, profile(), &section)
            .is_err()
    );
    assert!(
        frame
            .add_element_with_inertia_profile(0, 2, 5, profile(), &section)
            .is_err()
    );
    let huge = FrameInertiaProfile::new([1.0, f64::MAX], [0.0; 2]).unwrap();
    assert!(
        frame
            .add_element_with_inertia_profile(0, 1, 5, huge, &section)
            .is_err()
    );
    assert_eq!(frame.inertial_forces(&probe).unwrap(), original);
    assert_eq!(frame.reference_member_inertia(0).unwrap(), profile());
    assert!(frame.reference_member_inertia(1).is_err());
    assert_eq!(frame.displacements(), &[0.0; 6]);
    assert_eq!(frame.velocities(), &[0.0; 6]);
}
