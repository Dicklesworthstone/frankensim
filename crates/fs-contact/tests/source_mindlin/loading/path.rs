//! G1/G3/G4: original source admission reaches simultaneous body-coupled shear.

#[path = "path/sliding.rs"]
mod sliding;
use super::*;
use fs_contact::interface_binding::mindlin::loading::path::{
    SourceMindlinPathError, SourceMindlinPathInput,
};
use fs_tribo::mindlin::loading::path::MindlinPathError;

fn input(law: &SourceHertzMindlinContact, indentation: f64, slip: f64) -> SourceMindlinPathInput {
    SourceMindlinPathInput {
        indentation_m: indentation,
        free_slip_increment_m: slip,
        dt_s: 1e-3,
        compliance_m_per_n: 1e-4,
        solve: shear(law, &law.initial_state(), slip).solve,
    }
}

#[test]
fn original_sources_drive_simultaneous_recoil_and_separate_work() {
    let law = source(1e7, 0.6);
    let initial = law.initial_state();
    let slip = 0.8 * law.fixed_source().law().gross_slip_displacement_m();
    let mut state = initial.clone();
    for (indentation, released) in [(8e-5, slip), (5e-6, -slip), (0.0, slip)] {
        let request = input(&law, indentation, released);
        let result = law
            .advance_coupled_path(&state, &request, &mut || false)
            .unwrap();
        let step = &result.tangential;
        let actual = step.slip_increment_m[0];
        let work = step.traction_n[0] * actual;
        assert!(
            (actual + request.compliance_m_per_n * step.traction_n[0] - released).abs() < 1e-13
        );
        assert!(step.dissipated_j >= 0.0);
        assert!((step.stored_change_j + step.dissipated_j - work).abs() < 1e-13);
        common::close(
            result.normal_work_j,
            step.state.normal_energy_j() - state.normal_energy_j(),
        );
        assert_eq!(result.normal_work_j, result.normal_stored_change_j);
        assert_eq!(
            step.state.normal_force_n(),
            result.normal_response.normal_force_n
        );
        assert_eq!(
            result.normal_response.authority,
            law.fixed_source().normal_response().authority
        );
        assert_eq!(
            result.normal_response.uncertainty,
            law.fixed_source().normal_response().uncertainty
        );
        if indentation == 8e-5 {
            common::close(step.state.normal_force_n(), 8.0 * initial.normal_force_n());
            assert!(actual > 0.0 && actual < released);
        }
        state = step.state.clone();
    }
    assert_eq!(state.history_intervals(), 0);
    assert_eq!(state.normal_energy_j(), 0.0);
    assert_eq!(initial, law.initial_state());
}

#[test]
fn source_path_delegates_prescribed_motion_and_keeps_refusals_atomic() {
    let law = source(1e7, 0.6);
    let state = law.initial_state();
    let slip = 0.4 * law.fixed_source().law().gross_slip_displacement_m();
    let mut request = input(&law, 5e-5, slip);
    request.compliance_m_per_n = 0.0;
    let accepted = law
        .advance_coupled_path(&state, &request, &mut || false)
        .unwrap();
    let direct = law
        .law()
        .advance_prescribed_path(
            state.constitutive_state(),
            &TangentialStepInput {
                normal_n: accepted.normal_response.normal_force_n,
                slip_increment_m: [slip, 0.0],
                dt_s: request.dt_s,
                contact_temperature_k: Some(law.temperature_kelvin()),
                nominal_pressure_pa: None,
                compliance_m_per_n: 0.0,
                solve: request.solve,
            },
            &mut || false,
        )
        .unwrap();
    assert_eq!(
        accepted.tangential.state.constitutive_state(),
        &direct.tangential.state
    );
    assert_eq!(accepted.tangential.traction_n, direct.tangential.traction_n);
    assert_eq!(
        accepted.tangential.dissipated_j,
        direct.tangential.dissipated_j
    );

    for (indentation, dt) in [
        (2e-4, 1e-3),
        (8e-5, 1e-6),
        (-1e-5, 1e-3),
        (5e-5, 0.0),
        (5e-5, 2e-3),
        (f64::NAN, 1e-3),
    ] {
        let bad = SourceMindlinPathInput {
            indentation_m: indentation,
            dt_s: dt,
            ..request
        };
        assert!(
            law.advance_coupled_path(&state, &bad, &mut || false)
                .is_err()
        );
    }
    for compliance in [-1.0, f64::NAN, f64::INFINITY] {
        let bad = SourceMindlinPathInput {
            compliance_m_per_n: compliance,
            ..request
        };
        assert!(
            law.advance_coupled_path(&state, &bad, &mut || false)
                .is_err()
        );
    }
    assert!(matches!(
        law.advance_coupled_path(&state, &request, &mut || true),
        Err(SourceMindlinPathError::Path(MindlinPathError::Cancelled))
    ));
    let foreign = SourceHertzMindlinContact::from_fixed_source(
        resolve(normal_with_history(
            1e7,
            Some(0.6),
            Dims::NONE,
            false,
            "foreign history",
        ))
        .unwrap(),
    )
    .unwrap();
    assert!(
        foreign
            .advance_coupled_path(&state, &request, &mut || false)
            .is_err()
    );
    assert_eq!(state, law.initial_state());
    assert_eq!(
        accepted,
        law.advance_coupled_path(&state, &request, &mut || false)
            .unwrap()
    );
}
