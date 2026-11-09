//! Original cold/hot force endpoints, actual inertia and native continuation.

use super::*;
use fs_solid::Rayleigh;
use fs_solid::frame2d::laminate_thermal::dynamics::FrameLaminateNewmarkControls;

const LENGTH: f64 = 0.2;
const WIDTH: f64 = 0.002;
const DT: f64 = 0.0001;

fn controls() -> FrameLaminateNewmarkControls {
    FrameLaminateNewmarkControls {
        maximum_dofs: 6,
        maximum_members: 1,
        maximum_force_residual_n: 1e-7,
        maximum_moment_residual_n_m: 1e-8,
        maximum_displacement_defect_m: 1e-15,
        maximum_rotation_defect_rad: 1e-15,
        maximum_energy_defect_j: 1e-11,
    }
}
fn specimen(mass: FrameMassFormulation, cx: &Cx<'_>) -> (Frame2d, Ply) {
    let p = Ply {
        h: 0.003,
        rho: 2700.0,
        young: 2e9,
        poisson: 0.2,
        alpha: 1e-5,
        softening: 0.1,
    };
    let mut frame = Frame2d::with_mass_formulation(Geometry::Linear, mass);
    frame.add_node(0.0, 0.0);
    frame.add_node(LENGTH, 0.0);
    frame.fix(0, [true; 3]);
    frame.fix(1, [false, true, true]);
    frame
        .add_laminate_thermal_element(
            0,
            1,
            typed_source(&[p], 0.0, false),
            WIDTH,
            &vec![vec![300.0]; 3],
            0.02,
            cx,
        )
        .unwrap();
    (
        frame
            .with_laminate_thermoelasticity(0, &[[260.0, 400.0]], cx)
            .unwrap(),
        p,
    )
}
fn close(actual: f64, expected: f64, scale: f64) {
    assert!(
        (actual - expected).abs() <= 2e-9 * expected.abs().max(scale),
        "{actual:e} != {expected:e}"
    );
}

#[test]
fn changing_temperature_uses_original_force_and_retains_motion_with_both_native_masses() {
    with_cx(|cx| {
        for formulation in [
            FrameMassFormulation::Lumped,
            FrameMassFormulation::ConsistentReference,
        ] {
            let (mut frame, p) = specimen(formulation, cx);
            let mass = p.rho
                * WIDTH
                * p.h
                * LENGTH
                * match formulation {
                    FrameMassFormulation::Lumped => 0.5,
                    FrameMassFormulation::ConsistentReference => 1.0 / 3.0,
                };
            let (mut old_t, mut u, mut v) = (300.0, 0.0, 0.0);
            for t in [320.0, 340.0, 315.0] {
                // Independent one-free-DOF Newmark equations. In particular,
                // the initial acceleration belongs to old_t, never t.
                let k0 = young(p, old_t) * WIDTH * p.h / LENGTH;
                let k1 = young(p, t) * WIDTH * p.h / LENGTH;
                let f0 = k0 * (u - LENGTH * free(p, old_t));
                let a0 = -f0 / mass;
                let predictor = u + DT * v + 0.25 * DT * DT * a0;
                let a1 = -k1 * (predictor - LENGTH * free(p, t)) / (mass + 0.25 * DT * DT * k1);
                let u1 = predictor + 0.25 * DT * DT * a1;
                let v1 = v + 0.5 * DT * (a0 + a1);
                let f1 = k1 * (u1 - LENGTH * free(p, t));
                let before = frame.displacements().to_vec();
                let step = frame
                    .laminate_newmark_candidate(&[vec![vec![t]; 3]], &[0.0; 6], DT, controls(), cx)
                    .unwrap();
                assert_eq!(frame.displacements(), before);
                close(step.frame.displacements()[3], u1, 1e-8);
                close(step.frame.velocities()[3], v1, 1e-4);
                close(step.frame.accelerations()[3], a1, 1.0);
                close(step.internal_work_j, 0.5 * (f0 + f1) * (u1 - u), 1e-7);
                close(
                    step.kinetic_energy_change_j,
                    0.5 * mass * (v1 + v) * (v1 - v),
                    1e-7,
                );
                assert!(step.energy_balance_residual_j.abs() < 1e-11);
                assert_eq!(step.external_work_j, 0.0);
                assert_eq!(step.history.damping_energy, vec![0.0, 0.0]);
                assert_eq!(step.history.material_damping_energy, vec![0.0, 0.0]);
                let unit = [0.0, 0.0, 0.0, 1.0, 0.0, 0.0];
                close(step.frame.inertial_forces(&unit).unwrap()[3], mass, mass);
                frame = step.frame;
                (old_t, u, v) = (t, u1, v1);
            }
        }
    });
}

#[test]
fn isothermal_path_matches_existing_newmark_and_failed_hot_candidate_publishes_nothing() {
    with_cx(|cx| {
        let (mut frame, _) = specimen(FrameMassFormulation::ConsistentReference, cx);
        let damping = Rayleigh { a0: 0.0, a1: 0.0 };
        let velocity = [0.0, 0.0, 0.0, 0.1, 0.0, 0.0];
        let impulse = frame.inertial_forces(&velocity).unwrap();
        frame
            .apply_impulse(&impulse, &[0.0; 6], damping, 1e-10)
            .unwrap();
        let mut original = frame.clone();
        let history = original
            .newmark_with_load(&[0.0, 0.0], DT, (0.0, 0.0), damping, &[0.0; 6])
            .unwrap();
        let native = frame
            .laminate_newmark_candidate(&[vec![vec![300.0]; 3]], &[0.0; 6], DT, controls(), cx)
            .unwrap();
        for i in 0..6 {
            close(
                native.frame.displacements()[i],
                original.displacements()[i],
                1e-10,
            );
            close(native.frame.velocities()[i], original.velocities()[i], 1e-7);
        }
        close(native.internal_work_j, history.internal_work[1], 1e-7);
        let u = frame.displacements().to_vec();
        let v = frame.velocities().to_vec();
        assert!(
            frame
                .laminate_newmark_candidate(&[vec![vec![450.0]; 3]], &[0.0; 6], DT, controls(), cx)
                .is_err()
        );
        let mut small = controls();
        small.maximum_dofs = 5;
        assert!(
            frame
                .laminate_newmark_candidate(&[vec![vec![300.0]; 3]], &[0.0; 6], DT, small, cx)
                .is_err()
        );
        assert_eq!(frame.displacements(), u);
        assert_eq!(frame.velocities(), v);
    });
}
