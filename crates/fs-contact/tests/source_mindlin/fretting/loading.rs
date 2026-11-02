//! Complete load-path wear, source work separation and joint acceptance.
use super::*;
use fs_contact::interface_binding::mindlin::fretting::loading::SourceHertzMindlinWear;
use fs_contact::interface_binding::mindlin::loading::path::SourceMindlinPathInput;

fn request(law: &SourceHertzMindlinWear, indentation: f64, slip: f64) -> SourceMindlinPathInput {
    SourceMindlinPathInput {
        indentation_m: indentation,
        free_slip_increment_m: slip,
        dt_s: 1e-3,
        compliance_m_per_n: 0.0,
        solve: step(law.fixed_source().contact(), slip).solve,
    }
}

#[test]
fn unloading_and_lift_off_wear_follow_retained_local_sliding_not_endpoint_load() {
    let f = fixture(0.6, COEFFICIENTS, HARDNESS);
    let gross = f.wear.contact().law().gross_slip_displacement_m();
    let law = SourceHertzMindlinWear::from_fixed_source(f.wear).unwrap();
    let initial = law.initial_state();
    let mut state = initial.clone();
    let mut integrated = 0.0;
    for (indentation, slip) in [(5e-5, gross), (5e-6, 0.0), (0.0, 0.0), (0.0, -gross)] {
        let input = request(&law, indentation, slip);
        let original = law
            .contact()
            .advance_coupled_path(state.contact_state(), &input, &mut || false)
            .unwrap();
        let candidate = law
            .advance_coupled_path(&state, &input, &mut || false)
            .unwrap();
        let shear = &candidate.tangential;
        assert_eq!(shear.state.contact_state(), &original.tangential.state);
        assert_eq!(shear.traction_n, original.tangential.traction_n);
        assert_eq!(shear.dissipated_j, original.tangential.dissipated_j);
        assert_eq!(shear.stored_change_j, original.tangential.stored_change_j);
        assert_eq!(candidate.normal_response, original.normal_response);
        assert_eq!(candidate.normal_work_j, original.normal_work_j);
        assert_eq!(candidate.normal_work_j, candidate.normal_stored_change_j);
        let measure = candidate.pressure_weighted_sliding_j;
        near(measure, shear.dissipated_j / 0.6);
        integrated += measure;
        for i in 0..2 {
            near(
                candidate.wear[i].increment_m3,
                COEFFICIENTS[i] * measure / HARDNESS[i],
            );
            near(
                shear.state.volumes_m3()[i],
                COEFFICIENTS[i] * integrated / HARDNESS[i],
            );
        }
        if slip == 0.0 {
            assert_eq!(shear.slip_increment_m, [0.0; 2]);
            assert!(
                measure > 0.0,
                "unloading retained shear must produce fretting"
            );
            assert!(shear.stored_change_j < 0.0);
            assert!(candidate.normal_work_j < 0.0);
            assert!((shear.stored_change_j + shear.dissipated_j).abs() < 1e-13);
            if indentation == 0.0 {
                assert_eq!(candidate.normal_response.normal_force_n, 0.0);
                assert_eq!(shear.state.contact_state().history_intervals(), 0);
            }
        }
        if state.contact_state().normal_force_n() == 0.0 && indentation == 0.0 {
            assert_eq!(measure, 0.0);
            assert_eq!(shear.state.volumes_m3(), state.volumes_m3());
        }
        state = candidate.tangential.state;
    }
    // Recontact has no abandoned shear memory, while wear inventory survives.
    let recontact = law
        .advance_coupled_path(&state, &request(&law, 2e-5, 0.0), &mut || false)
        .unwrap();
    assert_eq!(recontact.pressure_weighted_sliding_j, 0.0);
    assert_eq!(recontact.tangential.state.volumes_m3(), state.volumes_m3());
    assert_eq!(recontact.tangential.traction_n, [0.0; 2]);
    assert_eq!(law.initial_state(), initial);
}

#[test]
fn frictionless_wear_integrates_nonlinear_pressure_through_grow_detach_and_onset() {
    let f = fixture(0.0, COEFFICIENTS, HARDNESS);
    let law = SourceHertzMindlinWear::from_fixed_source(f.wear).unwrap();
    let mut state = law.initial_state();
    let n0 = state.contact_state().normal_force_n();
    // d1=4*d0 gives N1=8*N0 and mean N=(2/5)*(32-1)/(4-1)*N0.
    for (indentation, travel, mean) in [
        (8e-5, -1e-5, (62.0 / 15.0) * n0),
        (0.0, 2e-5, 0.4 * 8.0 * n0),
        (0.0, -3e-5, 0.0),
        (2e-5, 1e-5, 0.4 * n0),
    ] {
        let result = law
            .advance_coupled_path(&state, &request(&law, indentation, travel), &mut || false)
            .unwrap();
        let integral = mean * travel.abs();
        near(result.pressure_weighted_sliding_j, integral);
        assert_eq!(result.tangential.dissipated_j, 0.0);
        assert_eq!(result.tangential.traction_n, [0.0; 2]);
        assert_eq!(result.tangential.slip_increment_m, [travel, 0.0]);
        for i in 0..2 {
            near(
                result.wear[i].increment_m3,
                COEFFICIENTS[i] * integral / HARDNESS[i],
            );
        }
        assert_eq!(result.normal_work_j, result.normal_stored_change_j);
        if indentation == 8e-5 {
            let endpoint_mean = 0.5 * (n0 + result.normal_response.normal_force_n);
            assert!((endpoint_mean - mean).abs() > 0.05 * mean);
        }
        state = result.tangential.state;
    }
}

#[test]
fn fixed_load_reduction_and_compliant_unloading_preserve_declared_wear_inventories() {
    let f = fixture(0.6, COEFFICIENTS, HARDNESS);
    let fixed = f.wear;
    let law = SourceHertzMindlinWear::from_fixed_source(fixed.clone()).unwrap();
    let initial = law.initial_state();
    let gross = fixed.contact().law().gross_slip_displacement_m();
    let mut input = request(&law, initial.contact_state().indentation_m(), 0.4 * gross);
    input.compliance_m_per_n = 0.5 / fixed.contact().law().initial_stiffness_n_per_m();
    let loaded = law
        .advance_coupled_path(&initial, &input, &mut || false)
        .unwrap();
    let fixed_answer = fixed
        .advance_candidate(
            &fixed.initial_state(),
            &TangentialStepInput {
                dt_s: input.dt_s,
                compliance_m_per_n: input.compliance_m_per_n,
                ..step(fixed.contact(), input.free_slip_increment_m)
            },
            &mut || false,
        )
        .unwrap();
    near(
        loaded.pressure_weighted_sliding_j,
        fixed_answer.pressure_weighted_sliding_j,
    );
    near(
        loaded.tangential.traction_n[0],
        fixed_answer.tangential.traction_n[0],
    );
    for i in 0..2 {
        near(
            loaded.wear[i].increment_m3,
            fixed_answer.wear[i].increment_m3,
        );
    }
    assert_eq!(loaded.normal_work_j, 0.0);
    let old_volumes = [2e-15, 5e-15];
    let state = law
        .state_from_contact(
            loaded.tangential.state.contact_state().clone(),
            old_volumes.map(|v| WearState::new(v).unwrap()),
        )
        .unwrap();
    input.indentation_m = 5e-6;
    input.free_slip_increment_m = 0.0;
    let recoil = law
        .advance_coupled_path(&state, &input, &mut || false)
        .unwrap();
    assert!(recoil.tangential.slip_increment_m[0] < 0.0);
    assert!(recoil.pressure_weighted_sliding_j > 0.0);
    assert!(recoil.tangential.stored_change_j < 0.0);
    for i in 0..2 {
        near(
            recoil.tangential.state.volumes_m3()[i],
            old_volumes[i] + recoil.wear[i].increment_m3,
        );
    }
    assert_eq!(state.volumes_m3(), old_volumes);
}

#[test]
fn late_cancellation_foreign_sources_and_second_surface_overflow_are_atomic() {
    let f = fixture(0.6, COEFFICIENTS, HARDNESS);
    let law = SourceHertzMindlinWear::from_fixed_source(f.wear).unwrap();
    let state = law.initial_state();
    let saved = state.clone();
    let gross = law
        .fixed_source()
        .contact()
        .law()
        .gross_slip_displacement_m();
    let input = request(&law, 5e-5, 0.6 * gross);
    let mut polls = 0;
    let expected = law
        .advance_coupled_path(&state, &input, &mut || {
            polls += 1;
            false
        })
        .unwrap();
    assert!(polls > 3);
    let mut reached = 0;
    let stopped = law.advance_coupled_path(&state, &input, &mut || {
        reached += 1;
        reached == polls
    });
    assert_eq!(stopped.unwrap_err(), SourceFrettingError::Cancelled);
    assert_eq!(reached, polls);
    let other = SourceHertzMindlinWear::from_fixed_source(
        fixture(0.6, [3e-6, COEFFICIENTS[1]], HARDNESS).wear,
    )
    .unwrap();
    assert!(
        other
            .advance_coupled_path(&state, &input, &mut || false)
            .is_err()
    );
    assert!(
        other
            .state_from_contact(state.contact_state().clone(), [WearState::default(); 2])
            .is_err()
    );
    assert_eq!(state, saved);
    assert_eq!(
        law.advance_coupled_path(&state, &input, &mut || false)
            .unwrap(),
        expected
    );
    // A is finite before B overflows: no shear, normal or wear checkpoint changes.
    let huge = SourceHertzMindlinWear::from_fixed_source(
        fixture(0.6, [COEFFICIENTS[0], 1e308], [HARDNESS[0], 1e-300]).wear,
    )
    .unwrap();
    let old = huge.initial_state();
    let before = old.clone();
    let bad = request(&huge, 5e-6, gross);
    assert!(
        huge.advance_coupled_path(&old, &bad, &mut || false)
            .is_err()
    );
    assert_eq!(old, before);
    let rest = huge
        .advance_coupled_path(&old, &request(&huge, 2e-5, 0.0), &mut || false)
        .unwrap();
    assert_eq!(rest.wear.map(|w| w.increment_m3), [0.0; 2]);
    assert_eq!(rest.tangential.state.volumes_m3(), [0.0; 2]);
}
