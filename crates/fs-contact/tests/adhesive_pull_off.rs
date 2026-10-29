//! G1/G3/G4: original mixed-mode potential, apparatus work and irreversible release.
#[path = "adhesive_dmt_coupled/common.rs"]
mod common;

use core::f64::consts::PI;
use fs_contact::normal_patch::adhesive::fixture::JkrSpringFixture;
use fs_contact::normal_patch::adhesive::mixed_mode::loading::JkrShearBoundaryKind;
use fs_contact::normal_patch::adhesive::mixed_mode::pull_off::*;
use fs_contact::normal_patch::adhesive::mixed_mode::{JkrMixedModeEnvelope, SourceJkrMixedMode};
use fs_contact::normal_patch::adhesive::{JkrEnvelope, SourceJkrSphereContact};
use fs_exec::CancelGate;

const W: f64 = 0.05;
const R: f64 = 0.01;
const E: f64 = 1.0 / (0.9375 / 1e11 + 0.9375 / 2e11);

fn length() -> f64 {
    (PI * W * R * R / (8.0 * E)).cbrt()
}
fn depth() -> f64 {
    length() * length() / R
}
fn stiffness() -> f64 {
    E * length()
}
fn sources(
    strength: f64,
) -> (
    SourceJkrMixedMode,
    fs_contact::interface_binding::adhesive_shear::SourceAdhesiveShear,
) {
    // Reuse original manufactured cards/resolvers. Only the constitutive model,
    // sphere radius and attraction-range applicability are selected here.
    let strength_pa = strength * W * R / (length() * length());
    let (normal, shear) = common::sources(W, strength_pa, 0.5 * strength_pa);
    let normal = SourceJkrSphereContact::new(
        normal.normal_state().clone(),
        R,
        1e-12,
        JkrEnvelope {
            maximum_indentation_m: 100.0 * depth(),
            half_space_depth_m: 0.1,
            layer_thickness_m: 0.1,
            maximum_patch_to_radius: 0.1,
            maximum_patch_to_depth: 0.1,
            maximum_patch_to_layer: 0.1,
            maximum_nominal_strain: 0.1,
            maximum_repulsive_pressure_to_yield: 0.1,
            maximum_process_zone_to_patch: 0.01,
            minimum_tabor_parameter: 5.0,
        },
    )
    .unwrap();
    (
        SourceJkrMixedMode::new(
            normal,
            JkrMixedModeEnvelope {
                maximum_tangential_displacement_to_radius: 0.1,
                maximum_mean_shear_stress_to_yield: 0.1,
            },
        )
        .unwrap(),
        shear,
    )
}
fn limits() -> JkrPullOffLimits {
    JkrPullOffLimits {
        maximum_contact_radius_m: 4.7 * length(),
        maximum_absolute_normal_stage_m: 1000.0 * depth(),
        maximum_iterations: 96,
        maximum_evaluations: 256,
    }
}
fn fixture(strength: f64, kn: f64, kt: f64, stage: [f64; 2]) -> JkrPullOffFixture {
    let (contact, shear) = sources(strength);
    JkrPullOffFixture::new(
        contact,
        shear,
        kn * stiffness(),
        kt * stiffness(),
        stage.map(|u| u * depth()),
        limits(),
        &CancelGate::new(),
    )
    .unwrap()
}
fn close(a: f64, b: f64, scale: f64, rel: f64) {
    assert!(
        (a - b).abs() <= scale * rel,
        "actual={a:e}, expected={b:e}, scale={scale:e}, tolerance={rel:e}"
    );
}

#[test]
fn g1_zero_shear_recovers_original_finite_normal_fixture_and_release() {
    for kn in [0.2, 4.0, 200.0] {
        let law = fixture(100.0, kn, 4.0, [0.0; 2]);
        let normal =
            JkrSpringFixture::new(law.contact().normal_source().clone(), kn * stiffness()).unwrap();
        let event = law.first_boundary();
        assert_eq!(event.kind, JkrShearBoundaryKind::RadiusInstability);
        assert!(event.relaxed_separation_admissible);
        close(
            event.response.state().normal_stage_m(),
            normal.minimum_attached_stage_m(),
            depth(),
            2e-12,
        );
        close(
            event.response.contact.unwrap().contact_radius_m,
            normal.critical_contact().contact_radius_m,
            length(),
            5e-14,
        );
        close(
            event.response.stage_force_n[0],
            normal.critical_contact().normal_force_n,
            W * R,
            3e-11,
        );
        for fraction in [0.05, 0.3, 0.8] {
            let stage = event.response.state().normal_stage_m()
                + fraction
                    * (law.initial_response().state().normal_stage_m()
                        - event.response.state().normal_stage_m());
            let mixed = law.evaluate_attached(stage, &CancelGate::new()).unwrap();
            let pure = normal.evaluate_attached(stage).unwrap();
            close(
                mixed.tip_coordinates_m[0],
                pure.contact.indentation_m,
                depth(),
                2e-11,
            );
            close(
                mixed.stage_force_n[0],
                pure.contact.normal_force_n,
                W * R,
                2e-10,
            );
            close(
                mixed.total_stored_energy_j,
                pure.total_stored_energy_j,
                W * length() * length(),
                1e-10,
            );
            close(
                mixed.stage_stiffness_n_per_m.unwrap()[0][0],
                pure.stage_tangent_n_per_m.unwrap(),
                stiffness(),
                1e-10,
            );
        }
    }
}

#[test]
fn g1_two_spring_equilibrium_and_full_stage_tangent_come_from_original_potential() {
    let shift = [3.0, -2.0];
    let law = fixture(100.0, 4.0, 4.0, shift);
    let stage = 0.4 * law.first_boundary().response.state().normal_stage_m()
        + 0.6 * law.initial_response().state().normal_stage_m();
    let center = law.evaluate_attached(stage, &CancelGate::new()).unwrap();
    let contact = center.contact.unwrap();
    let a = contact.contact_radius_m;
    let delta = center.tip_coordinates_m[0];
    let u = contact.tangential_displacement_m;
    let k = 8.0 * law.contact().reduced_shear_modulus_pa();
    let kn = law.normal_stiffness_n_per_m();
    let kt = law.tangential_stiffness_n_per_m();
    let normal_elastic =
        E * a * ((delta - a * a / (3.0 * R)).powi(2) + 4.0 * a.powi(4) / (45.0 * R * R));
    let shear_elastic = 0.5 * k * a * (u[0] * u[0] + u[1] * u[1]);
    let springs = 0.5 * kn * (stage - delta).powi(2)
        + 0.5 * kt * ((shift[0] * depth() - u[0]).powi(2) + (shift[1] * depth() - u[1]).powi(2));
    close(
        center.total_stored_energy_j,
        normal_elastic + shear_elastic - PI * W * a * a + springs,
        W * length() * length(),
        1e-10,
    );
    let h = 4.0 * E * a * (a * a / R - delta) / R - 2.0 * PI * W;
    let h_machine = h
        - (2.0 * E * (a * a / R - delta)).powi(2) / (2.0 * E * a + kn)
        - k * k * (u[0] * u[0] + u[1] * u[1]) / (k * a + kt);
    close(
        center.apparatus_radius_curvature_n_per_m,
        h_machine,
        W,
        2e-12,
    );
    let tangent = center.stage_stiffness_n_per_m.unwrap();
    let step = 1e-4 * depth();
    for j in 0..3 {
        let (left, right) = if j == 0 {
            (
                law.evaluate_attached(stage - step, &CancelGate::new())
                    .unwrap(),
                law.evaluate_attached(stage + step, &CancelGate::new())
                    .unwrap(),
            )
        } else {
            let mut minus = shift;
            let mut plus = shift;
            minus[j - 1] -= step / depth();
            plus[j - 1] += step / depth();
            (
                fixture(100.0, 4.0, 4.0, minus)
                    .evaluate_attached(stage, &CancelGate::new())
                    .unwrap(),
                fixture(100.0, 4.0, 4.0, plus)
                    .evaluate_attached(stage, &CancelGate::new())
                    .unwrap(),
            )
        };
        close(
            (right.total_stored_energy_j - left.total_stored_energy_j) / (2.0 * step),
            center.stage_force_n[j],
            W * R,
            3e-6,
        );
        for i in 0..3 {
            close(
                (right.stage_force_n[i] - left.stage_force_n[i]) / (2.0 * step),
                tangent[i][j],
                stiffness(),
                2e-7,
            );
            close(tangent[i][j], tangent[j][i], stiffness(), 1e-14);
        }
    }
    for residual in center.spring_force_residual_n {
        close(residual, 0.0, W * R, 1e-10);
    }
}

#[test]
fn g1_withdrawal_localizes_mixed_mode_separation_and_transfers_energy_once() {
    let law = fixture(100.0, 4.0, 4.0, [3.0, 2.0]);
    let before = law.initial_response();
    let boundary = law.first_boundary();
    assert!(boundary.relaxed_separation_admissible);
    let target = boundary.response.state().normal_stage_m() - 5.0 * depth();
    let whole = law
        .advance(before.state(), target, &CancelGate::new())
        .unwrap();
    let event = whole.event.unwrap();
    assert_eq!(event.kind, JkrPullOffEventKind::Separation);
    assert_eq!(
        event.normal_stage_m,
        boundary.response.state().normal_stage_m()
    );
    assert_eq!(
        event.before_tip_coordinates_m,
        boundary.response.tip_coordinates_m
    );
    assert_eq!(event.after_tip_coordinates_m[0], event.normal_stage_m);
    assert_eq!(
        [
            event.after_tip_coordinates_m[1],
            event.after_tip_coordinates_m[2]
        ],
        law.tangential_stage_m()
    );
    assert_eq!(whole.response.state().branch(), JkrPullOffBranch::Detached);
    assert!(whole.response.contact.is_none());
    assert_eq!(whole.untraversed_normal_stage_m, 0.0);
    assert_eq!(whole.response.stage_force_n, [0.0; 3]);
    assert_eq!(whole.response.total_stored_energy_j, 0.0);
    assert_eq!(
        whole.released_energy_j,
        boundary.response.total_stored_energy_j
    );
    assert_eq!(
        whole.response.state().accumulated_released_energy_j(),
        whole.released_energy_j
    );
    close(
        whole.prescribed_work_j,
        whole.response.total_stored_energy_j - before.total_stored_energy_j
            + whole.released_energy_j,
        W * length() * length(),
        1e-14,
    );
    let repeat = law
        .advance(whole.response.state(), target - depth(), &CancelGate::new())
        .unwrap();
    assert_eq!(repeat.event, None);
    assert_eq!(repeat.released_energy_j, 0.0);
    assert_eq!(repeat.prescribed_work_j, 0.0);
    assert_eq!(
        repeat.response.state().accumulated_released_energy_j(),
        whole.released_energy_j
    );
    let returns = law
        .advance(repeat.response.state(), 0.0, &CancelGate::new())
        .unwrap();
    assert_eq!(
        returns.response.state().branch(),
        JkrPullOffBranch::Detached
    );
    assert!(matches!(
        law.advance(returns.response.state(), depth(), &CancelGate::new()),
        Err(JkrPullOffError::DetachedOverlap)
    ));

    // Independent quadrature of N dD along a, finite at the stage fold.
    let a0 = before.contact.unwrap().contact_radius_m;
    let a1 = boundary.response.contact.unwrap().contact_radius_m;
    let k = 8.0 * law.contact().reduced_shear_modulus_pa();
    let kn = law.normal_stiffness_n_per_m();
    let kt = law.tangential_stiffness_n_per_m();
    let dt = law.tangential_stage_m()[0].hypot(law.tangential_stage_m()[1]);
    let n = 40000;
    let da = (a1 - a0) / f64::from(n);
    let mut work = 0.0;
    for i in 0..n {
        let a = a0 + (f64::from(i) + 0.5) * da;
        let u = kt * dt / (kt + k * a);
        let opening = ((2.0 * PI * W * a - 0.5 * k * u * u) / E).sqrt();
        let delta = a * a / R - opening;
        let force = 2.0 * E * a * (delta - a * a / (3.0 * R));
        let opening_derivative =
            (2.0 * PI * W + k * k * u * u / (kt + k * a)) / (2.0 * E * opening);
        let delta_derivative = 2.0 * a / R - opening_derivative;
        let force_derivative = 2.0 * E * (delta - a * a / R) + 2.0 * E * a * delta_derivative;
        work += force * (delta_derivative + force_derivative / kn) * da;
    }
    close(whole.prescribed_work_j, work, W * length() * length(), 3e-7);
    let middle =
        0.5 * (before.state().normal_stage_m() + boundary.response.state().normal_stage_m());
    let first = law
        .advance(before.state(), middle, &CancelGate::new())
        .unwrap();
    let second = law
        .advance(first.response.state(), target, &CancelGate::new())
        .unwrap();
    assert_eq!(second.response, whole.response);
    assert_eq!(second.released_energy_j, whole.released_energy_j);
    close(
        first.prescribed_work_j + second.prescribed_work_j,
        whole.prescribed_work_j,
        W * length() * length(),
        3e-13,
    );
}

#[test]
fn g1_weak_junction_and_interpenetrating_relaxed_geometry_stop_at_actual_handoff() {
    let weak = fixture(0.1, 4.0, 4.0, [3.0, 2.0]);
    let boundary = weak.first_boundary();
    assert_eq!(boundary.kind, JkrShearBoundaryKind::GrossSlip);
    assert!(!boundary.relaxed_separation_admissible);
    close(boundary.response.static_margin_n, 0.0, W * R, 1e-12);
    let target = boundary.response.state().normal_stage_m() - 3.0 * depth();
    let step = weak
        .advance(weak.initial_response().state(), target, &CancelGate::new())
        .unwrap();
    assert_eq!(step.event.unwrap().kind, JkrPullOffEventKind::GrossSlip);
    assert_eq!(step.response, boundary.response);
    assert!(step.untraversed_normal_stage_m < 0.0);
    assert_eq!(step.released_energy_j, 0.0);
    assert!(step.response.contact.is_some());

    let held = fixture(100.0, 4.0, 4.0, [10.0, 0.0]);
    let fold = held.first_boundary();
    assert_eq!(fold.kind, JkrShearBoundaryKind::RadiusInstability);
    assert!(fold.response.state().normal_stage_m() > 0.0);
    assert!(!fold.relaxed_separation_admissible);
    let stop = held
        .advance(
            held.initial_response().state(),
            -depth(),
            &CancelGate::new(),
        )
        .unwrap();
    assert_eq!(
        stop.event.unwrap().kind,
        JkrPullOffEventKind::RadiusInstability
    );
    assert_eq!(stop.response, fold.response);
    assert_eq!(stop.released_energy_j, 0.0);
    assert!(stop.untraversed_normal_stage_m < 0.0);
}

#[test]
fn g3_tangent_rotation_and_apparatus_compliance_change_physics_consistently() {
    let a = fixture(100.0, 4.0, 4.0, [3.0, 4.0]);
    let b = fixture(100.0, 4.0, 4.0, [-4.0, 3.0]);
    let stage = 0.5
        * (a.initial_response().state().normal_stage_m()
            + a.first_boundary().response.state().normal_stage_m());
    let x = a.evaluate_attached(stage, &CancelGate::new()).unwrap();
    let y = b.evaluate_attached(stage, &CancelGate::new()).unwrap();
    assert_eq!(
        x.contact.unwrap().contact_radius_m,
        y.contact.unwrap().contact_radius_m
    );
    assert_eq!(x.total_stored_energy_j, y.total_stored_energy_j);
    assert_eq!(x.stage_force_n[0], y.stage_force_n[0]);
    assert_eq!(x.stage_force_n[1], y.stage_force_n[2]);
    assert_eq!(x.stage_force_n[2], -y.stage_force_n[1]);
    let mut old_stage = f64::NEG_INFINITY;
    for kn in [0.2, 1.0, 4.0, 200.0] {
        let law = fixture(100.0, kn, 4.0, [2.0, 0.0]);
        let fold = law.first_boundary().response;
        assert!(fold.state().normal_stage_m() > old_stage);
        old_stage = fold.state().normal_stage_m();
        close(fold.apparatus_radius_curvature_n_per_m, 0.0, W, 2e-12);
        assert!(fold.stage_stiffness_n_per_m.is_none());
    }
}

#[test]
fn g0_g4_source_budget_envelope_and_cancel_refusals_leave_checkpoint_for_retry() {
    let law = fixture(100.0, 4.0, 4.0, [3.0, 2.0]);
    let state = law.initial_response().state();
    let target =
        0.5 * (state.normal_stage_m() + law.first_boundary().response.state().normal_stage_m());
    let expected = law.advance(state, target, &CancelGate::new()).unwrap();
    let other = fixture(100.0, 5.0, 4.0, [3.0, 2.0]);
    assert_eq!(
        other.advance(state, target, &CancelGate::new()),
        Err(JkrPullOffError::StateMismatch)
    );
    let cancelled = CancelGate::new();
    cancelled.request();
    assert_eq!(
        law.advance(state, target, &cancelled),
        Err(JkrPullOffError::Cancelled)
    );
    assert_eq!(
        law.advance(state, f64::NAN, &CancelGate::new()),
        Err(JkrPullOffError::OutsideStageEnvelope)
    );
    assert_eq!(
        law.advance(
            state,
            2.0 * limits().maximum_absolute_normal_stage_m,
            &CancelGate::new()
        ),
        Err(JkrPullOffError::OutsideStageEnvelope)
    );
    assert_eq!(
        law.advance(state, target, &CancelGate::new()).unwrap(),
        expected
    );
    assert_eq!(state, law.initial_response().state());
    let (normal, shear) = sources(100.0);
    let (_, wrong) = sources(99.0);
    assert!(matches!(
        JkrPullOffFixture::new(
            normal.clone(),
            wrong,
            4.0 * stiffness(),
            4.0 * stiffness(),
            [0.0; 2],
            limits(),
            &CancelGate::new()
        ),
        Err(JkrPullOffError::StateMismatch)
    ));
    let mut small = limits();
    small.maximum_evaluations = 1;
    assert!(matches!(
        JkrPullOffFixture::new(
            normal.clone(),
            shear.clone(),
            4.0 * stiffness(),
            4.0 * stiffness(),
            [0.0; 2],
            small,
            &CancelGate::new()
        ),
        Err(JkrPullOffError::Budget)
    ));
    assert!(matches!(
        JkrPullOffFixture::new(
            normal,
            shear,
            4.0 * stiffness(),
            4.0 * stiffness(),
            [0.0; 2],
            limits(),
            &cancelled
        ),
        Err(JkrPullOffError::Cancelled)
    ));
}
