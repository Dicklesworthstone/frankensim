//! The retained material source drives both the spatial plate and contact law.
#[path = "../examples/material_plate_impact.rs"]
#[allow(dead_code)]
mod example;

use example::{
    GAP_M, STEP_S, THICKNESS_M, build_impact, build_plate, initial_state, method, synthetic_source,
};
use fs_plate_transient::NoPlateLoad;
use fs_time::galpha::{SecondOrderProblem, SecondOrderState, TimeSolveError};
use fs_tribo::InputAuthority;

fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 3.0e-10 * actual.abs().max(expected.abs()).max(1.0e-25),
        "{actual:.17e} != {expected:.17e}",
    );
}

#[test]
fn sourced_plate_contact_matches_original_force_potential_and_loss() {
    let source = synthetic_source("synthetic response", 1.0e7, 900.0, 0.4).unwrap();
    let (mesh, plate) = build_plate(&source).unwrap();
    let no_load = NoPlateLoad;
    let impact = build_impact(&mesh, &plate, source.clone(), &no_load).unwrap();
    assert_eq!(impact.source(), Some(&source));
    assert_eq!(source.envelope().half_space_depth_m, THICKNESS_M);
    assert_eq!(source.envelope().layer_thickness_m, THICKNESS_M);
    let mut q = vec![0.0; impact.dimension()];
    q[impact.plate_dofs()] = GAP_M + 3.0e-6;
    for rate in [-0.07, 0.0, 0.07] {
        let mut v = vec![0.0; impact.dimension()];
        v[impact.plate_dofs()] = rate;
        let expected = source
            .evaluate("plate sample", 3.0e-6, rate, STEP_S)
            .unwrap();
        let observed = impact.probe(&q, &v).unwrap();
        close(observed.normal_force_n, expected.normal_force_n);
        close(observed.contact_energy_j, expected.reversible_energy_j);
        close(
            observed.contact_dissipated_power_w,
            expected.dissipated_power_w,
        );
        assert!(observed.contact_dissipated_power_w >= 0.0);
        assert_eq!(expected.authority, InputAuthority::SyntheticFixture);
    }
}

#[test]
fn a_bulk_material_swap_changes_plate_stiffness_mass_and_contact_together() {
    let soft = synthetic_source("synthetic soft", 1.0e7, 900.0, 0.4).unwrap();
    let stiff = synthetic_source("synthetic stiff", 3.0e7, 1800.0, 0.4).unwrap();
    let (mesh_a, plate_a) = build_plate(&soft).unwrap();
    let (mesh_b, plate_b) = build_plate(&stiff).unwrap();
    assert_eq!(plate_a.free, plate_b.free);
    assert_eq!(plate_a.dof_map, plate_b.dof_map);
    for i in 0..plate_a.free {
        close(plate_b.k.get(i, i), 3.0 * plate_a.k.get(i, i));
        close(plate_b.m.get(i, i), 2.0 * plate_a.m.get(i, i));
    }
    assert_ne!(soft.model().identity(), stiff.model().identity());
    let ratio = stiff.model().normal_state().elastic().reduced_modulus_pa()
        / soft.model().normal_state().elastic().reduced_modulus_pa();
    close(
        stiff.stiffness_n_per_m_3_2(),
        ratio * soft.stiffness_n_per_m_3_2(),
    );
    assert!(ratio > 1.0);
    let no_load = NoPlateLoad;
    let a = build_impact(&mesh_a, &plate_a, soft, &no_load).unwrap();
    let b = build_impact(&mesh_b, &plate_b, stiff, &no_load).unwrap();
    let mut q = vec![0.0; a.dimension()];
    let mut v = vec![0.0; a.dimension()];
    q[a.plate_dofs()] = GAP_M + 3.0e-6;
    v[a.plate_dofs()] = 0.04;
    let pa = a.probe(&q, &v).unwrap();
    let pb = b.probe(&q, &v).unwrap();
    close(pb.normal_force_n, ratio * pa.normal_force_n);
    close(pb.contact_energy_j, ratio * pa.contact_energy_j);
    close(
        pb.contact_dissipated_power_w,
        ratio * pa.contact_dissipated_power_w,
    );
}

#[test]
fn source_step_penetration_and_rate_limits_refuse_without_state_publication() {
    let source = synthetic_source("synthetic limits", 1.0e7, 900.0, 0.4).unwrap();
    let envelope = source.envelope();
    let (mesh, plate) = build_plate(&source).unwrap();
    let no_load = NoPlateLoad;
    let impact = build_impact(&mesh, &plate, source, &no_load).unwrap();
    let forcing = vec![0.0; impact.dimension()];
    for case in 0..3 {
        let mut state = initial_state(&impact).unwrap();
        let (step, expected_reason) = match case {
            0 => (
                1.01 * envelope.maximum_step_s,
                "source point sample duration",
            ),
            1 => {
                state.q[impact.plate_dofs()] = GAP_M + 1.01 * envelope.maximum_indentation_m;
                (STEP_S, "source point indentation")
            }
            _ => {
                state.v[impact.plate_dofs()] = -1.01 * envelope.maximum_absolute_rate_m_per_s;
                (STEP_S, "source point absolute rate")
            }
        };
        let before = state.clone();
        let error = method(impact.dimension(), step)
            .step(&mut state, &impact, &forcing)
            .unwrap_err();
        assert!(error.to_string().contains(expected_reason), "{error}");
        assert_eq!(state, before);
    }
}

#[test]
fn a_valid_implicit_stage_cannot_publish_an_out_of_envelope_endpoint() {
    let source = synthetic_source("synthetic endpoint", 1.0e7, 900.0, 0.4).unwrap();
    let step = source.envelope().maximum_step_s;
    let (mesh, plate) = build_plate(&source).unwrap();
    let no_load = NoPlateLoad;
    let impact = build_impact(&mesh, &plate, source, &no_load).unwrap();
    let n = impact.dimension();
    let body = impact.plate_dofs();
    let zero = vec![0.0; n];
    let mut state = SecondOrderState::new(0.0, &zero, &zero, &zero);
    let before = state.clone();
    let mut forcing = vec![0.0; n];
    forcing[body] = 1.8;
    // With rho=1 the stage velocity is half the endpoint velocity. The
    // force produces about .09 m/s at the stage and .18 m/s at the endpoint;
    // the source rate ceiling is .1 m/s. Contact reduces these only slightly.
    let mut stage_q = vec![0.0; n];
    let mut stage_v = vec![0.0; n];
    stage_q[body] = 4.5e-6;
    stage_v[body] = 0.09;
    impact.validate_state(&stage_q, &stage_v).unwrap();
    let error = method(n, step)
        .step(&mut state, &impact, &forcing)
        .unwrap_err();
    assert!(
        matches!(
            &error,
            TimeSolveError::ModelRefused { phase: "accepted endpoint", reason }
                if reason.contains("source point absolute rate")
        ),
        "{error}"
    );
    assert_eq!(state, before);
}
