//! Independent finite-force contact checks. Point fixtures have a complete
//! analytic constant-acceleration path; the elastic fixture admits its one
//! sliding mode by a conservative force/velocity bound over the whole interval.

use fs_contact::sustained::NaturalMapScales;
use fs_solid::frame2d::FrameMassFormulation;
use fs_solid::frame2d::sustained::{
    FrameSustainedContact, FrameSustainedControls, FrameSustainedError, FrameSustainedPath,
    FrameSustainedStep,
};
use fs_solid::{Frame2d, Geometry, Rayleigh, SectionModel};
use std::cell::Cell;

const ZERO: Rayleigh = Rayleigh { a0: 0.0, a1: 0.0 };
const NORMAL: [f64; 3] = [0.0, 1.0, 0.0];
const TANGENT: [f64; 3] = [1.0, 0.0, 0.0];

fn close(a: f64, b: f64) {
    assert!((a - b).abs() <= 2e-8 * (1.0 + b.abs()), "{a:e} != {b:e}");
}
fn point() -> Frame2d {
    let mut frame = Frame2d::new(Geometry::Linear);
    frame.add_node(0.0, 0.0);
    frame.add_mass(0, 2.0);
    frame.fix(0, [false, false, true]);
    frame
}
fn velocity(frame: &mut Frame2d, v: [f64; 2]) {
    let jump = [
        2.0 * (v[0] - frame.velocities()[0]),
        2.0 * (v[1] - frame.velocities()[1]),
        0.0,
    ];
    frame.apply_impulse(&jump, &[0.0; 3], ZERO, 0.125).unwrap();
    close(frame.velocities()[0], v[0]);
    close(frame.velocities()[1], v[1]);
}
fn control(h: f64) -> FrameSustainedControls {
    FrameSustainedControls {
        duration_s: h,
        ..Default::default()
    }
}
fn contact(h: f64, mu: f64, obstacle: [f64; 2]) -> FrameSustainedContact<'static> {
    FrameSustainedContact {
        id: 7,
        normal_covector: &NORMAL,
        tangent_covector: &TANGENT,
        gap_m: 0.0,
        obstacle_velocity: obstacle,
        friction_coefficient: mu,
        scales: NaturalMapScales {
            normal_n_per_m: 4.0 / (h * h),
            tangent_n_per_m: 4.0 / (h * h),
        },
    }
}
// Actual physical point, orthonormal plane traces, and one held load. Here the
// complete original trajectory is quadratic; checking the quadratic minimum
// and affine slip velocities is an analytic whole-path test, not sampling.
fn admit_point(path: &FrameSustainedPath<'_, '_>) -> Result<&'static str, &'static str> {
    let h = path.duration_s;
    let c = path.contacts[0];
    let force = path.forces_n[0];
    let mut acceleration = [0.0; 3];
    for (i, a) in acceleration.iter_mut().enumerate().take(2) {
        *a = (path.external_load[i]
            + c.normal_covector[i] * force[0]
            + c.tangent_covector[i] * force[1])
            / 2.0;
        close(path.after.accelerations()[i], *a);
        close(
            path.after.velocities()[i],
            path.before.velocities()[i] + h * *a,
        );
        close(
            path.actual_displacement[i],
            h * path.before.velocities()[i] + 0.5 * h * h * *a,
        );
    }
    let dot = |row: &[f64], x: &[f64]| row.iter().zip(x).map(|(a, b)| a * b).sum::<f64>();
    let vn = dot(c.normal_covector, path.before.velocities()) - c.obstacle_velocity[0];
    let an = dot(c.normal_covector, &acceleration);
    let gap = |t: f64| c.gap_m + vn * t + 0.5 * an * t * t;
    let mut least = gap(0.0).min(gap(h));
    if an > 0.0 {
        let t = -vn / an;
        if t > 0.0 && t < h {
            least = least.min(gap(t));
        }
    }
    if least < -1e-9 {
        return Err("quadratic path penetrates plane");
    }
    if force[0] > 1e-7 && (vn.abs() > 1e-7 || an.abs() > 1e-7 || c.gap_m > 1e-9) {
        return Err("loaded point does not dwell");
    }
    let vt = dot(c.tangent_covector, path.before.velocities()) - c.obstacle_velocity[1];
    let at = dot(c.tangent_covector, &acceleration);
    if force[0] * c.friction_coefficient > 1e-7 {
        if vt * (vt + h * at) < -1e-12 {
            return Err("affine slip crosses a stop");
        }
        let capacity = c.friction_coefficient * force[0];
        if force[1].abs() < capacity - 1e-7 && (vt.abs() > 1e-7 || at.abs() > 1e-7) {
            return Err("interior cone needs complete sticking");
        }
    }
    Ok("analytic point path")
}
fn step(
    frame: &mut Frame2d,
    c: FrameSustainedContact<'_>,
    warm: [f64; 2],
    load: [f64; 3],
    h: f64,
) -> FrameSustainedStep<&'static str> {
    frame
        .advance_sustained_contacts(&[c], &[warm], &load, control(h), || false, admit_point)
        .unwrap()
}
fn ledger(step: &FrameSustainedStep<&str>) {
    close(step.energy_balance_residual_j, 0.0);
    close(step.native_energy_balance_residual_j, 0.0);
    close(step.normal_contact_work_j, 0.0);
    close(step.passivity_residual_j, 0.0);
    assert_eq!(step.response_evaluations, step.solve.response_evaluations);
    assert!(step.response_evaluations >= 2);
    assert!(step.solve.forces_n.iter().all(|v| v[2] == 0.0));
    assert!(step.solve.paths_m.iter().all(|v| v[2] == 0.0));
}

#[test]
fn loaded_rest_sticking_and_threshold_have_no_artificial_impulse_or_heat() {
    for tangent_load in [0.0, 5.0, 10.0] {
        let mut frame = point();
        let s = step(
            &mut frame,
            contact(0.125, 0.5, [0.0; 2]),
            [20.0, -tangent_load],
            [tangent_load, -20.0, 0.0],
            0.125,
        );
        assert_eq!(frame.displacements(), &[0.0; 3]);
        assert_eq!(frame.velocities(), &[0.0; 3]);
        assert_eq!(frame.accelerations(), &[0.0; 3]);
        close(s.solve.forces_n[0][0], 20.0);
        close(s.solve.forces_n[0][1], -tangent_load);
        close(s.frictional_dissipation_j, 0.0);
        assert_eq!(s.generalized_momentum_change, vec![0.0; 3]);
        assert_eq!(s.path_admission, "analytic point path");
        ledger(&s);
    }
}

#[test]
fn native_force_solve_recovers_unknown_resting_reactions_and_counts_probes() {
    let mut frame = point();
    let s = step(
        &mut frame,
        contact(0.125, 0.5, [0.0; 2]),
        [0.0, 0.0],
        [5.0, -20.0, 0.0],
        0.125,
    );
    close(s.solve.forces_n[0][0], 20.0);
    close(s.solve.forces_n[0][1], -5.0);
    assert!(s.response_evaluations > 3);
    close(frame.velocities()[0], 0.0);
    close(frame.velocities()[1], 0.0);
    ledger(&s);
}

#[test]
fn sliding_and_breakaway_match_continuous_force_acceleration_and_work() {
    for v0 in [0.0, 1.0] {
        let h = 0.125;
        let mut frame = point();
        velocity(&mut frame, [v0, 0.0]);
        let s = step(
            &mut frame,
            contact(h, 0.5, [0.0; 2]),
            [20.0, -10.0],
            [12.0, -20.0, 0.0],
            h,
        );
        let dx = h * v0 + 0.5 * h * h;
        close(frame.accelerations()[0], 1.0);
        close(frame.velocities()[0], v0 + h);
        close(frame.displacements()[0], dx);
        close(s.generalized_momentum_change[0], 2.0 * h);
        close(s.kinetic_energy_change_j, (v0 + h).powi(2) - v0.powi(2));
        close(s.external_work_j, 12.0 * dx);
        close(s.frictional_dissipation_j, 10.0 * dx);
        close(s.obstacle_work_j, 0.0);
        ledger(&s);
    }
}

#[test]
fn separating_load_detaches_and_carries_no_tensile_or_friction_force() {
    let h = 0.125;
    let mut frame = point();
    let s = step(
        &mut frame,
        contact(h, 0.5, [0.0; 2]),
        [20.0, 0.0],
        [0.0, 6.0, 0.0],
        h,
    );
    close(s.solve.forces_n[0][0], 0.0);
    close(s.solve.forces_n[0][1], 0.0);
    close(frame.accelerations()[1], 3.0);
    close(frame.velocities()[1], 3.0 * h);
    close(s.solve.paths_m[0][0], 1.5 * h * h);
    close(s.external_work_j, 9.0 * h * h);
    close(s.frictional_dissipation_j, 0.0);
    ledger(&s);
}

#[test]
fn moving_belt_supplies_world_work_and_only_relative_slip_becomes_heat() {
    let h = 0.125;
    let mut frame = point();
    let s = step(
        &mut frame,
        contact(h, 0.5, [0.0, 2.0]),
        [20.0, 10.0],
        [0.0, -20.0, 0.0],
        h,
    );
    close(frame.accelerations()[0], 5.0);
    close(frame.velocities()[0], 0.625);
    close(frame.displacements()[0], 0.0390625);
    close(s.solve.paths_m[0][1], -0.2109375);
    close(s.obstacle_work_j, 2.5);
    close(s.frictional_dissipation_j, 2.109375);
    close(s.kinetic_energy_change_j, 0.390625);
    ledger(&s);
}

#[test]
fn lifting_obstacle_work_balances_gravity_without_normal_heat() {
    let h = 0.125;
    let mut frame = point();
    velocity(&mut frame, [0.0, 0.25]);
    let s = step(
        &mut frame,
        contact(h, 0.5, [0.25, 0.0]),
        [20.0, 0.0],
        [0.0, -20.0, 0.0],
        h,
    );
    close(frame.velocities()[1], 0.25);
    close(frame.accelerations()[1], 0.0);
    close(s.obstacle_work_j, 0.625);
    close(s.external_work_j, -0.625);
    close(s.kinetic_energy_change_j, 0.0);
    close(s.frictional_dissipation_j, 0.0);
    ledger(&s);
}

#[test]
fn an_interval_crossing_a_stop_refuses_even_when_net_slip_keeps_its_sign() {
    let h = 0.3;
    let mut frame = point();
    velocity(&mut frame, [1.0, 0.0]);
    // True stop is .2 s; the inadmissible held-friction continuation ends at
    // -.5 m/s yet still has +.075 m net displacement.
    close(1.0 * h - 0.5 * 5.0 * h * h, 0.075);
    let before = format!("{frame:?}");
    let called = Cell::new(false);
    let error = frame
        .advance_sustained_contacts(
            &[contact(h, 0.5, [0.0; 2])],
            &[[20.0, -10.0]],
            &[0.0, -20.0, 0.0],
            control(h),
            || false,
            |_| {
                called.set(true);
                Ok::<_, &str>(())
            },
        )
        .unwrap_err();
    assert!(
        matches!(error, FrameSustainedError::ContactMode { .. }),
        "{error:?}"
    );
    assert!(!called.get());
    assert_eq!(format!("{frame:?}"), before);
}

#[test]
fn stopped_point_can_continue_in_stick_with_recomputed_finite_acceleration() {
    let mut frame = point();
    velocity(&mut frame, [1.0, 0.0]);
    let stop = step(
        &mut frame,
        contact(0.2, 0.5, [0.0; 2]),
        [20.0, -10.0],
        [0.0, -20.0, 0.0],
        0.2,
    );
    close(frame.velocities()[0], 0.0);
    close(frame.displacements()[0], 0.1);
    close(stop.frictional_dissipation_j, 1.0);
    let hold = step(
        &mut frame,
        contact(0.125, 0.5, [0.0; 2]),
        [20.0, -3.0],
        [3.0, -20.0, 0.0],
        0.125,
    );
    close(frame.velocities()[0], 0.0);
    close(frame.accelerations()[0], 0.0);
    close(frame.displacements()[0], 0.1);
    close(hold.frictional_dissipation_j, 0.0);
    frame
        .newmark_with_load(&[0.0; 2], 0.125, (0.0, 0.0), ZERO, &[0.0; 3])
        .unwrap();
    close(frame.velocities()[0], 0.0);
    close(frame.accelerations()[0], 0.0);
}

fn elastic_frame() -> Frame2d {
    let mut frame =
        Frame2d::with_mass_formulation(Geometry::Linear, FrameMassFormulation::ConsistentReference);
    frame.add_node(0.0, 0.0);
    frame.add_node(1.0, 0.0);
    frame.fix(0, [true; 3]);
    frame.fix(1, [false, false, true]);
    frame.add_mass(1, 2.0);
    frame
        .add_element(0, 1, 3, 6.0, &|| SectionModel::Elastic {
            ea: 4000.0,
            ei: 1000.0,
        })
        .unwrap();
    frame
}
fn sdof(q: f64, v: f64, force: f64, h: f64) -> (f64, f64, f64) {
    let (m, k, c) = (4.0, 4000.0, 4.8);
    let a0 = (force - k * q - c * v) / m;
    let a1 = (force - k * (q + h * v + 0.25 * h * h * a0) - c * (v + 0.5 * h * a0))
        / (m + 0.5 * h * c + 0.25 * h * h * k);
    (
        q + h * v + 0.25 * h * h * (a0 + a1),
        v + 0.5 * h * (a0 + a1),
        a1,
    )
}

#[test]
fn original_consistent_mass_support_reaction_elastic_damping_energy_and_continuation() {
    let h = 0.001;
    let damping = Rayleigh { a0: 0.2, a1: 0.001 };
    let mut frame = elastic_frame();
    frame
        .static_load(&[0.0, 0.0, 0.0, 4.0, 0.0, 0.0], 1)
        .unwrap();
    frame
        .apply_impulse(&[0.0, 0.0, 0.0, 4.0, 0.0, 0.0], &[0.0; 6], damping, h)
        .unwrap();
    let q0 = frame.displacements()[3];
    let v0 = frame.velocities()[3];
    close(q0, 0.001);
    close(v0, 1.0);
    let normal = [0.0, 0.0, 0.0, 0.0, 1.0, 0.0];
    let tangent = [0.0, 0.0, 0.0, 1.0, 0.0, 0.0];
    let c = FrameSustainedContact {
        id: 31,
        normal_covector: &normal,
        tangent_covector: &tangent,
        gap_m: 0.0,
        obstacle_velocity: [0.0; 2],
        friction_coefficient: 0.25,
        scales: NaturalMapScales {
            normal_n_per_m: 8e6,
            tangent_n_per_m: 8e6,
        },
    };
    let ctrl = FrameSustainedControls {
        damping,
        ..control(h)
    };
    let s = frame
        .advance_sustained_contacts(
            &[c],
            &[[10.0, -2.5]],
            &[0.0, 0.0, 0.0, 8.0, -10.0, 0.0],
            ctrl,
            || false,
            |path| {
                // Normal/axial blocks decouple in this straight linear member.
                // Normal load cancels identically. For the physical axial ODE,
                // q in [0,.004], v in [0,2] gives |a| <= (5.5+16+9.6)/4.
                // Over .001 s this bound cannot reach v=0 or 2 or leave q range.
                // Thus constant sliding direction holds over the entire interval.
                let bound = (5.5 + 4000.0 * 0.004 + 4.8 * 2.0) / 4.0;
                assert!(v0 - h * bound > 0.0 && v0 + h * bound < 2.0);
                assert!(q0 > 0.0 && q0 + 2.0 * h < 0.004);
                close(path.after.displacements()[4], 0.0);
                Ok::<_, &str>("bounded elastic path")
            },
        )
        .unwrap();
    let (q1, v1, a1) = sdof(q0, v0, 5.5, h);
    close(frame.displacements()[3], q1);
    close(frame.velocities()[3], v1);
    close(frame.accelerations()[3], a1);
    let dx = q1 - q0;
    let vm = 0.5 * (v0 + v1);
    let qm = 0.5 * (q0 + q1);
    close(s.generalized_momentum_change[3], 4.0 * (v1 - v0));
    close(s.generalized_momentum_change[0], v1 - v0);
    // Original consistent M_cf=1, K_cf=-4000, C_cf=.2-4=-3.8.
    close(s.support_impulse[0], v1 - v0 - h * (4000.0 * qm + 3.8 * vm));
    close(s.support_impulse[3], 0.0);
    close(s.strain_energy_change_j, 2000.0 * (q1 * q1 - q0 * q0));
    close(s.kinetic_energy_change_j, 2.0 * (v1 * v1 - v0 * v0));
    close(s.rayleigh_damping_energy_j, 4.8 * h * vm * vm);
    close(s.material_damping_energy_j, 0.0);
    close(s.external_work_j, 8.0 * dx);
    close(s.frictional_dissipation_j, 2.5 * dx);
    ledger(&s);
    let continued = sdof(q1, v1, 5.5, h);
    frame
        .newmark_with_load(
            &[0.0; 2],
            h,
            (0.0, 0.0),
            damping,
            &[0.0, 0.0, 0.0, 5.5, 0.0, 0.0],
        )
        .unwrap();
    close(frame.displacements()[3], continued.0);
    close(frame.velocities()[3], continued.1);
    close(frame.accelerations()[3], continued.2);
}

#[test]
fn pretranslated_small_representable_motion_uses_physical_displacement_allowance() {
    let mut frame = point();
    velocity(&mut frame, [0.1, 0.0]);
    frame
        .newmark_with_load(&[0.0; 2], 1.0, (0.0, 0.0), ZERO, &[0.0; 3])
        .unwrap();
    velocity(&mut frame, [1e-4, 0.0]);
    let before = frame.displacements()[0];
    let s = step(
        &mut frame,
        contact(0.01, 0.0, [0.0; 2]),
        [20.0, 0.0],
        [0.0, -20.0, 0.0],
        0.01,
    );
    close(before, 0.1);
    close(s.displacement_change[0], 1e-6);
    assert!(s.displacement_change[0] > 0.0);
    close(frame.velocities()[0], 1e-4);
    ledger(&s);
}

#[test]
fn lost_motion_on_large_stored_displacement_cannot_publish_fictitious_contact_work() {
    let mut frame = point();
    velocity(&mut frame, [1e16, 0.0]);
    frame
        .newmark_with_load(&[0.0; 2], 1.0, (0.0, 0.0), ZERO, &[0.0; 3])
        .unwrap();
    velocity(&mut frame, [0.0, 0.0]);
    velocity(&mut frame, [1.0, 0.0]);
    let before = format!("{frame:?}");
    let called = Cell::new(false);
    let error = frame
        .advance_sustained_contacts(
            &[contact(0.125, 0.0, [0.0; 2])],
            &[[20.0, 0.0]],
            &[0.0, -20.0, 0.0],
            control(0.125),
            || false,
            |_| {
                called.set(true);
                Ok::<_, &str>(())
            },
        )
        .unwrap_err();
    assert!(
        matches!(
            error,
            FrameSustainedError::Balance {
                quantity: "actual Newmark displacement",
                ..
            }
        ),
        "{error:?}"
    );
    assert!(!called.get());
    assert_eq!(format!("{frame:?}"), before);
}

#[test]
fn path_refusal_and_late_cancellation_preserve_complete_native_state() {
    let mut frame = point();
    let before = format!("{frame:?}");
    let c = contact(0.125, 0.5, [0.0; 2]);
    let error = frame
        .advance_sustained_contacts(
            &[c],
            &[[20.0, 0.0]],
            &[0.0, -20.0, 0.0],
            control(0.125),
            || false,
            |_| Err::<(), _>("physical feature coverage incomplete"),
        )
        .unwrap_err();
    assert!(matches!(
        error,
        FrameSustainedError::Path("physical feature coverage incomplete")
    ));
    assert_eq!(format!("{frame:?}"), before);
    let admitted = Cell::new(false);
    let error = frame
        .advance_sustained_contacts(
            &[c],
            &[[20.0, 0.0]],
            &[0.0, -20.0, 0.0],
            control(0.125),
            || admitted.get(),
            |_| {
                admitted.set(true);
                Ok::<_, &str>(())
            },
        )
        .unwrap_err();
    assert!(matches!(error, FrameSustainedError::Cancelled));
    assert!(admitted.get());
    assert_eq!(format!("{frame:?}"), before);
}

#[test]
fn response_and_workspace_budgets_refuse_atomically() {
    for which in 0..3 {
        let mut frame = point();
        let before = format!("{frame:?}");
        let mut ctrl = control(0.125);
        match which {
            0 => ctrl.solve.maximum_response_evaluations = 1,
            1 => ctrl.maximum_workspace_entries = 1,
            _ => ctrl.maximum_dofs = 2,
        }
        let error = frame
            .advance_sustained_contacts(
                &[contact(0.125, 0.5, [0.0; 2])],
                &[[0.0; 2]],
                &[0.0, -20.0, 0.0],
                ctrl,
                || false,
                |_| Ok::<_, &str>(()),
            )
            .unwrap_err();
        assert!(
            matches!(
                error,
                FrameSustainedError::Solve(_) | FrameSustainedError::Budget { .. }
            ),
            "{error:?}"
        );
        assert_eq!(format!("{frame:?}"), before);
    }
}

#[test]
fn closing_impact_bad_domains_and_massless_motion_refuse_before_path_admission() {
    for which in 0..5 {
        let mut frame = point();
        let mut c = contact(0.125, 0.5, [0.0; 2]);
        let mut ctrl = control(0.125);
        match which {
            0 => velocity(&mut frame, [0.0, -1.0]),
            1 => c.gap_m = -1e-15,
            2 => c.friction_coefficient = -0.1,
            3 => ctrl.duration_s = 0.0,
            _ => {
                frame = Frame2d::new(Geometry::Linear);
                frame.add_node(0.0, 0.0);
                frame.fix(0, [false, false, true]);
            }
        }
        let before = format!("{frame:?}");
        let called = Cell::new(false);
        let result = frame.advance_sustained_contacts(
            &[c],
            &[[20.0, 0.0]],
            &[0.0, -20.0, 0.0],
            ctrl,
            || false,
            |_| {
                called.set(true);
                Ok::<_, &str>(())
            },
        );
        assert!(result.is_err());
        assert!(!called.get());
        assert_eq!(format!("{frame:?}"), before);
    }
}

#[test]
fn pre_rotated_motion_uses_radians_independently_of_contact_metres() {
    let h = 0.01;
    let mut original =
        Frame2d::with_mass_formulation(Geometry::Linear, FrameMassFormulation::ConsistentReference);
    original.add_node(0.0, 0.0);
    original.add_node(1.0, 0.0);
    original.fix(0, [true; 3]);
    original.fix(1, [true, true, false]);
    original
        .add_element(0, 1, 3, 6.0, &|| SectionModel::Elastic {
            ea: 4000.0,
            ei: 1000.0,
        })
        .unwrap();
    let load = [0.0, 0.0, 0.0, 0.0, 0.0, 400.0];
    original.static_load(&load, 1).unwrap();
    // Original Hermite consistent rotary entry is rho*A*L^3/105.
    let inertia = 6.0 / 105.0;
    original
        .apply_impulse(&[0.0, 0.0, 0.0, 0.0, 0.0, inertia * 1e-4], &load, ZERO, h)
        .unwrap();
    close(original.displacements()[5], 0.1);
    close(original.velocities()[5], 1e-4);
    let normal = [0.0, 0.0, 0.0, 0.0, 1.0, 0.0];
    let tangent = [0.0, 0.0, 0.0, 1.0, 0.0, 0.0];
    let c = FrameSustainedContact {
        id: 101,
        normal_covector: &normal,
        tangent_covector: &tangent,
        gap_m: 0.0,
        obstacle_velocity: [0.0; 2],
        friction_coefficient: 0.0,
        scales: NaturalMapScales {
            normal_n_per_m: 1.0,
            tangent_n_per_m: 1.0,
        },
    };
    // Contact is at the translationally fixed node: its physical trace stays
    // exactly fixed at all times, independently of the allowed section spin.
    let mut frame = original.clone();
    let ctrl = FrameSustainedControls {
        displacement_tolerance_m: 0.0,
        rotation_tolerance_rad: 1e-9,
        ..control(h)
    };
    let s = frame
        .advance_sustained_contacts(
            &[c],
            &[[0.0; 2]],
            &load,
            ctrl,
            || false,
            |path| {
                assert_eq!(&path.actual_displacement[3..5], &[0.0; 2]);
                Ok::<_, &str>("fixed node contact")
            },
        )
        .unwrap();
    assert!(s.displacement_change[5] > 0.0);
    assert_eq!(s.displacement_change[3], 0.0);
    ledger(&s);
    let before = format!("{original:?}");
    let error = original
        .advance_sustained_contacts(
            &[c],
            &[[0.0; 2]],
            &load,
            FrameSustainedControls {
                displacement_tolerance_m: 1e-9,
                rotation_tolerance_rad: 0.0,
                ..control(h)
            },
            || false,
            |_| Ok::<_, &str>(()),
        )
        .unwrap_err();
    assert!(
        matches!(
            error,
            FrameSustainedError::Balance {
                quantity: "actual Newmark displacement",
                ..
            }
        ),
        "{error:?}"
    );
    assert_eq!(format!("{original:?}"), before);
}
