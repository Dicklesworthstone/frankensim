//! Original source receipts reach the native radial fields without a new law.

use super::*;
use fs_contact::interface_binding::mindlin::loading::SourceHertzMindlinContact;
use fs_tribo::mindlin::surface::MindlinSurfaceResponse;

#[test]
fn source_hertz_pressure_material_stiffness_and_friction_reach_local_fields() {
    let mut peaks = Vec::new();
    for (young, mu) in [(1e7, 0.6), (3e7, 0.6), (1e7, 0.2)] {
        let law = resolve(normal(young, Some(mu), Dims::NONE, false)).unwrap();
        let state = law
            .advance(
                &law.initial_state(),
                &step(&law, 1.2 * law.law().gross_slip_displacement_m()),
            )
            .unwrap()
            .state;
        let saved_state = state.clone();
        let saved_source = law.clone();
        let a = law.normal_response().patch_radius_m;
        let peak = law.normal_response().pressure.peak_pressure_pa;
        peaks.push(peak);
        for fraction in [0.0_f64, 0.4, 0.8] {
            let response = law.surface_response(&state, fraction * a).unwrap();
            common::close(
                response.normal_pressure_pa,
                peak * (1.0 - fraction * fraction).sqrt(),
            );
            common::close(
                response.tangential_traction_pa,
                mu * response.normal_pressure_pa,
            );
            assert_eq!(
                response,
                law.law()
                    .surface_response(state.constitutive_state(), fraction * a)
                    .unwrap()
            );
        }
        let exterior = law.surface_response(&state, 2.0 * a).unwrap();
        assert_eq!(exterior.normal_pressure_pa, 0.0);
        assert_eq!(exterior.tangential_traction_pa, 0.0);
        assert!(exterior.elastic_displacement_m > 0.0);
        assert_eq!(state, saved_state);
        assert_eq!(law, saved_source);
    }
    assert!(peaks[1] > peaks[0]);
    assert_eq!(peaks[2], peaks[0]);
}

#[test]
fn changing_patch_open_and_recontact_fields_retain_original_source_and_memory() {
    let fixed = resolve(normal(1e7, Some(0.6), Dims::NONE, false)).unwrap();
    let law = SourceHertzMindlinContact::from_fixed_source(fixed.clone()).unwrap();
    let initial = law.initial_state();
    let input = step(&fixed, 0.6 * fixed.law().gross_slip_displacement_m());
    let loaded = law.advance(&initial, &input).unwrap().state;
    let grown = law.normal_step(&loaded, 8e-5, 1e-3).unwrap();
    let radius = 0.4 * grown.state.contact_radius_m();
    let response = law.surface_response(&grown.state, radius).unwrap();
    common::close(
        response.normal_pressure_pa,
        grown.normal_response.pressure.peak_pressure_pa * (1.0_f64 - 0.16).sqrt(),
    );
    assert_eq!(
        response,
        law.law()
            .surface_response(grown.state.constitutive_state(), radius)
            .unwrap()
    );
    assert!(response.tangential_traction_pa > 0.0);
    assert_eq!(
        grown.normal_response.input_source_id,
        fixed.normal_response().input_source_id
    );
    let open = law.normal_step(&grown.state, 0.0, 1e-3).unwrap().state;
    for radius in [
        0.0,
        initial.contact_radius_m(),
        2.0 * grown.state.contact_radius_m(),
    ] {
        assert_eq!(
            law.surface_response(&open, radius).unwrap(),
            MindlinSurfaceResponse::default()
        );
    }
    let recontact = law
        .normal_step(&open, initial.indentation_m(), 1e-3)
        .unwrap()
        .state;
    let restored = law.surface_response(&recontact, 0.0).unwrap();
    common::close(
        restored.normal_pressure_pa,
        fixed.normal_response().pressure.peak_pressure_pa,
    );
    assert_eq!(restored.tangential_traction_pa, 0.0);
    assert_eq!(restored.elastic_displacement_m, 0.0);
    let repeated = law.advance(&recontact, &input).unwrap().state;
    for fraction in [0.0, 0.5, 1.5] {
        let radius = fraction * initial.contact_radius_m();
        assert_eq!(
            law.surface_response(&repeated, radius).unwrap(),
            law.surface_response(&loaded, radius).unwrap()
        );
    }
    assert_eq!(law.fixed_source(), &fixed);
}

#[test]
fn equal_physics_different_source_history_and_invalid_radii_do_not_change_state() {
    let law = resolve(normal(1e7, Some(0.6), Dims::NONE, false)).unwrap();
    let foreign = resolve(normal_with_history(
        1e7,
        Some(0.6),
        Dims::NONE,
        false,
        "different original surface history",
    ))
    .unwrap();
    assert_eq!(law.law(), foreign.law());
    let input = step(&law, 0.3 * law.law().gross_slip_displacement_m());
    let fixed_state = law.advance(&law.initial_state(), &input).unwrap().state;
    let changing = SourceHertzMindlinContact::from_fixed_source(law.clone()).unwrap();
    let changing_foreign = SourceHertzMindlinContact::from_fixed_source(foreign.clone()).unwrap();
    let state = changing
        .advance(&changing.initial_state(), &input)
        .unwrap()
        .state;
    let saved_fixed = fixed_state.clone();
    let saved_changing = state.clone();
    for radius in [0.0, f64::NAN] {
        assert!(matches!(
            foreign.surface_response(&fixed_state, radius),
            Err(SourceMindlinError::Invalid(_))
        ));
        assert!(matches!(
            changing_foreign.surface_response(&state, radius),
            Err(SourceMindlinError::Invalid(_))
        ));
    }
    for radius in [-1.0, f64::INFINITY, f64::NAN] {
        assert!(law.surface_response(&fixed_state, radius).is_err());
        assert!(changing.surface_response(&state, radius).is_err());
    }
    assert_eq!(fixed_state, saved_fixed);
    assert_eq!(state, saved_changing);
    let fixed = law.surface_response(&fixed_state, 0.0).unwrap();
    let evolved = changing.surface_response(&state, 0.0).unwrap();
    // The changing-load owner reconstructs its radius from the original Hertz
    // load; that independent material computation agrees to floating roundoff.
    common::close(fixed.normal_pressure_pa, evolved.normal_pressure_pa);
    common::close(fixed.tangential_traction_pa, evolved.tangential_traction_pa);
    common::close(fixed.elastic_displacement_m, evolved.elastic_displacement_m);
}
