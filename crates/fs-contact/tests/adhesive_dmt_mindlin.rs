//! G1/G3/G4: sourced adhesion, actual-area microslip, reversible work and refusals.
#[path = "adhesive_dmt_mindlin/common.rs"]
mod common;
use common::*;
use core::f64::consts::PI;
use fs_contact::normal_patch::dmt::mindlin::*;
use fs_material::state_point::MaterialPropertySelection;
use fs_qty::Pressure;
use fs_tribo::{InputAuthority, TangentialLaw, TangentialSolveConfig, TangentialStepInput};

fn input(indentation: f64, slip: f64, compliance: f64) -> SourceMindlinPathInput {
    SourceMindlinPathInput {
        indentation_m: indentation,
        free_slip_increment_m: slip,
        dt_s: 1e-6,
        compliance_m_per_n: compliance,
        solve: TangentialSolveConfig {
            max_iterations: 96,
            maximum_force_n: 1e-3,
            force_absolute_tolerance_n: 1e-20,
            force_relative_tolerance: 1e-12,
        },
    }
}
fn close(actual: f64, expected: f64, absolute: f64, relative: f64) {
    assert!(
        (actual - expected).abs() <= absolute + relative * actual.abs().max(expected.abs()),
        "{actual:e} != {expected:e}"
    );
}
fn scales(source: &SourceDmtMindlinContact, state: &DmtMindlinState) -> (f64, f64, f64) {
    let p = source.law().parameters();
    let compliance: f64 = p
        .shear_moduli_pa
        .into_iter()
        .zip(p.poisson_ratios)
        .map(|(g, nu)| (2.0 - nu) / g)
        .sum();
    let k = 8.0 * state.contact_radius_m() / compliance;
    let capacity = p.friction_coefficient * state.normal_response().repulsive_normal_force_n;
    (k, capacity, 1.5 * capacity / k)
}

#[test]
fn zero_and_tensile_net_load_retain_independent_hertz_microslip_and_surface_traction() {
    let source = source(0.01, Some(0.6), "initial").unwrap();
    let normal = source.normal_source();
    let pull = -normal.force_controlled_pull_off_n();
    for load in [0.0, -0.7 * pull] {
        let response = normal.evaluate_contact_load(load).unwrap();
        let state = source.initialize(response.indentation_m).unwrap();
        let (k, capacity, d) = scales(&source, &state);
        close(capacity, 0.6 * (load + pull), 1e-22, 3e-14);
        let u = 0.3 * d;
        let expected_endpoint = capacity * (1.0 - 0.7_f64.powf(1.5));
        let expected_average = capacity * (1.0 - (0.4 / 0.3) * (1.0 - 0.7_f64.powf(2.5)));
        let trial = source
            .trial_coupled_path(&state, &input(response.indentation_m, u, 0.0), &mut || {
                false
            })
            .unwrap();
        let step = source
            .advance_coupled_path(&state, &input(response.indentation_m, u, 0.0), &mut || {
                false
            })
            .unwrap();
        close(
            step.tangential.state.endpoint_force_n(),
            expected_endpoint,
            1e-22,
            2e-13,
        );
        close(
            step.tangential.traction_n[0],
            expected_average,
            1e-22,
            2e-13,
        );
        close(
            trial.tangential.average_force_n,
            expected_average,
            1e-22,
            2e-13,
        );
        close(
            source.law().parameters().friction_coefficient,
            0.6,
            0.0,
            0.0,
        );
        close(step.mean_normal_force_n, load, 1e-20, 1e-13);
        assert!(k > 0.0 && step.tangential.dissipated_j > 0.0);
        assert_eq!(step.normal_work_j, 0.0);
        assert_eq!(
            step.normal_response.authority,
            InputAuthority::SyntheticFixture
        );
        assert_eq!(
            step.normal_response.normal_state_identity(),
            normal.normal_state().identity()
        );

        let a = state.contact_radius_m();
        let c = a * 0.7_f64.sqrt();
        for radius in [0.0, 0.5 * c, 0.5 * (a + c)] {
            let field = source
                .surface_response(&step.tangential.state, radius)
                .unwrap();
            let e = normal.normal_state().elastic().reduced_modulus_pa();
            let pressure = 2.0 * e / (PI * R) * (a * a - radius * radius).sqrt();
            let corrective = if radius < c {
                2.0 * e / (PI * R) * (c * c - radius * radius).sqrt()
            } else {
                0.0
            };
            close(field.normal_pressure_pa, pressure, 1e-6, 2e-13);
            close(
                field.tangential_traction_pa,
                0.6 * (pressure - corrective),
                1e-6,
                2e-12,
            );
        }
        close(
            step.tangential.stored_change_j + step.tangential.dissipated_j,
            expected_average * u,
            1e-31,
            2e-12,
        );
    }
}

#[test]
fn evolving_adhesive_area_uses_original_memory_and_pays_unloading_microslip_once() {
    let source = source(0.01, Some(0.6), "initial").unwrap();
    let initial = source.initialize(0.25 * Z).unwrap();
    let (k, _, d) = scales(&source, &initial);
    let mut state = source
        .advance_coupled_path(&initial, &input(0.25 * Z, 0.8 * d, 0.0), &mut || false)
        .unwrap()
        .tangential
        .state;
    let unloaded = source
        .advance_coupled_path(&state, &input(0.04 * Z, 0.0, 0.0), &mut || false)
        .unwrap();
    assert!(unloaded.tangential.dissipated_j > 0.0);
    close(
        unloaded.tangential.stored_change_j,
        -unloaded.tangential.dissipated_j,
        1e-30,
        2e-12,
    );
    assert_eq!(unloaded.tangential.slip_increment_m, [0.0; 2]);
    state = unloaded.tangential.state;
    for (delta, du) in [
        (0.8 * Z, 0.4 * d),
        (0.08 * Z, -0.5 * d),
        (0.4 * Z, -1.5 * d),
        (0.0, 0.2 * d),
        (0.3 * Z, 0.1 * d),
    ] {
        let saved = state.clone();
        let request = input(delta, du, 0.7 / k);
        let trial = source
            .trial_coupled_path(&state, &request, &mut || false)
            .unwrap();
        let step = source
            .advance_coupled_path(&state, &request, &mut || false)
            .unwrap();
        let original = source
            .law()
            .advance_coupled_path(
                state.constitutive_state(),
                &TangentialStepInput {
                    normal_n: step.normal_response.repulsive_normal_force_n,
                    slip_increment_m: [du, 0.0],
                    dt_s: request.dt_s,
                    contact_temperature_k: Some(source.temperature_kelvin()),
                    nominal_pressure_pa: None,
                    compliance_m_per_n: request.compliance_m_per_n,
                    solve: request.solve,
                },
                &mut || false,
            )
            .unwrap();
        assert_eq!(state, saved);
        assert_eq!(
            step.tangential.state.constitutive_state(),
            &original.tangential.state
        );
        assert_eq!(step.tangential.traction_n, original.tangential.traction_n);
        assert_eq!(
            step.tangential.dissipated_j,
            original.tangential.dissipated_j
        );
        assert_eq!(
            trial.tangential.average_force_n,
            step.tangential.traction_n[0]
        );
        assert_eq!(
            trial.tangential.stored_change_j,
            step.tangential.stored_change_j
        );
        close(
            step.mean_normal_force_n * (delta - state.indentation_m()),
            step.normal_work_j,
            1e-30,
            2e-14,
        );
        close(
            step.normal_work_j,
            step.normal_response.reversible_energy_j - state.normal_response().reversible_energy_j,
            1e-30,
            2e-13,
        );
        close(
            step.tangential.slip_increment_m[0]
                + request.compliance_m_per_n * step.tangential.traction_n[0],
            du,
            1e-20,
            2e-12,
        );
        state = step.tangential.state;
    }
}

#[test]
fn adhesion_changes_normal_work_and_force_controlled_area_without_inventing_a_shear_coefficient() {
    let first = source(0.01, Some(0.6), "initial").unwrap();
    let second = source(0.02, Some(0.6), "initial").unwrap();
    let delta0 = 0.03 * Z;
    let delta1 = 0.12 * Z;
    let a = first.initialize(delta0).unwrap();
    let b = second.initialize(delta0).unwrap();
    let request = input(delta1, 0.01 * Z, 0.0);
    let x = first
        .advance_coupled_path(&a, &request, &mut || false)
        .unwrap();
    let y = second
        .advance_coupled_path(&b, &request, &mut || false)
        .unwrap();
    assert_eq!(x.tangential.traction_n, y.tangential.traction_n);
    assert_eq!(x.tangential.dissipated_j, y.tangential.dissipated_j);
    let extra_pull = 2.0 * PI * R * 0.01;
    close(
        y.mean_normal_force_n - x.mean_normal_force_n,
        -extra_pull,
        1e-20,
        2e-13,
    );
    close(
        y.normal_work_j - x.normal_work_j,
        -extra_pull * (delta1 - delta0),
        1e-30,
        2e-13,
    );
    let r0 = first
        .normal_source()
        .evaluate_contact_load(0.0)
        .unwrap()
        .contact_radius_m;
    let r1 = second
        .normal_source()
        .evaluate_contact_load(0.0)
        .unwrap()
        .contact_radius_m;
    close(r1 / r0, 2.0_f64.cbrt(), 1e-14, 1e-14);
}

#[test]
fn normal_mean_force_matches_exact_work_across_gap_contact_and_adjacent_endpoints() {
    let normal = normal(0.01, Some(0.6), "initial");
    for (from, to) in [
        (-2.0 * Z, -1.5 * Z),
        (-1.8 * Z, -0.3 * Z),
        (-0.8 * Z, 0.4 * Z),
        (0.0, 0.7 * Z),
        (0.2 * Z, 0.8 * Z),
        (0.8 * Z, 0.2 * Z),
        (-Z, -Z),
        (0.0, 0.0),
    ] {
        let mean = normal.mean_force_between(from, to).unwrap();
        let work = normal.reversible_work_between(from, to).unwrap();
        close(mean * (to - from), work, 1e-30, 3e-14);
        close(normal.mean_force_between(to, from).unwrap(), mean, 0.0, 0.0);
        if from == to {
            assert_eq!(mean, normal.evaluate(from).unwrap().normal_force_n);
        }
    }
    for from in [-0.5 * Z, 0.5 * Z] {
        let to = f64::from_bits(from.to_bits() + 1);
        let mean = normal.mean_force_between(from, to).unwrap();
        close(
            mean,
            normal.evaluate(from).unwrap().normal_force_n,
            1e-20,
            5e-15,
        );
        close(
            mean * (to - from),
            normal.reversible_work_between(from, to).unwrap(),
            1e-45,
            5e-15,
        );
    }
}

#[test]
fn original_sources_gap_boundary_and_failed_proposals_preserve_accepted_memory() {
    let source = source(0.01, Some(0.6), "initial").unwrap();
    let state = source.initialize(0.1 * Z).unwrap();
    let saved = state.clone();
    let request = input(0.2 * Z, 0.01 * Z, 0.0);
    let expected = source
        .advance_coupled_path(&state, &request, &mut || false)
        .unwrap();
    for foreign in [
        common::source(0.02, Some(0.6), "initial").unwrap(),
        common::source(0.01, Some(0.6), "different specimen history").unwrap(),
    ] {
        assert_eq!(
            foreign.trial_coupled_path(&state, &request, &mut || false),
            Err(DmtMindlinError::StateMismatch)
        );
    }
    assert!(common::source(0.01, None, "initial").is_err());
    assert!(
        SourceDmtMindlinContact::resolve(
            normal_with_dimensions(0.01, Some(0.6), "initial", Pressure::DIMS),
            requirement(),
            MaterialPropertySelection::SingleClaimOnly,
            64
        )
        .is_err()
    );
    assert!(source.initialize(-0.1 * Z).is_err());
    for invalid in [
        SourceMindlinPathInput {
            indentation_m: -0.1 * Z,
            ..request
        },
        SourceMindlinPathInput {
            indentation_m: 2.0 * Z,
            ..request
        },
        SourceMindlinPathInput {
            dt_s: 0.0,
            ..request
        },
        SourceMindlinPathInput {
            compliance_m_per_n: -1.0,
            ..request
        },
    ] {
        assert!(
            source
                .advance_coupled_path(&state, &invalid, &mut || false)
                .is_err()
        );
        assert!(
            source
                .trial_coupled_path(&state, &invalid, &mut || false)
                .is_err()
        );
    }
    assert!(
        source
            .advance_coupled_path(&state, &request, &mut || true)
            .is_err()
    );
    assert!(
        source
            .trial_coupled_path(&state, &request, &mut || true)
            .is_err()
    );
    assert_eq!(state, saved);
    assert_eq!(
        source
            .advance_coupled_path(&state, &request, &mut || false)
            .unwrap(),
        expected
    );

    let zero = common::source(0.01, Some(0.0), "initial").unwrap();
    let state = zero.initialize(0.1 * Z).unwrap();
    let step = zero
        .advance_coupled_path(&state, &request, &mut || false)
        .unwrap();
    assert_eq!(step.tangential.traction_n, [0.0; 2]);
    assert_eq!(step.tangential.dissipated_j, 0.0);
    assert_eq!(step.tangential.state.tangential_energy_j(), 0.0);
    assert_ne!(step.normal_work_j, 0.0);
}
