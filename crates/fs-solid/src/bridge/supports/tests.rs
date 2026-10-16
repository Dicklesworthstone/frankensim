use super::*;
use crate::bridge::tests::{birth, config as beam_config, material, temperatures};

fn config() -> Config {
    Config { contact: GramOrthantConfig { max_variables: 8, max_columns: 128,
        max_setup_products: 8192, max_sweeps: 4096, gradient_tolerance: 1e-13 },
        gap_tolerance_m: 1e-11, complementarity_tolerance_j: 1e-13 }
}
fn stop(x: f64, lower: f64) -> Stop {
    Stop { id: 1, reference_x_m: x, minimum_deflection_m: lower }
}
fn near(a: f64, b: f64, tolerance: f64) {
    assert!((a-b).abs() < tolerance, "{a:e} != {b:e}");
}

// G1: continuum elastic references; point-contact error is independently refined.
#[test]
fn propped_cantilever_and_midspan_support_reactions_refine_to_analytic_limits() {
    let law = material(1e8, vec![]);
    let weight = 1200.0*1e-6*9.81*0.02;
    for (support, x, fraction) in [(Support::CantileverLeft, 0.02, 3.0/8.0),
        (Support::SimplySupported, 0.01, 5.0/8.0)] {
        let mut errors = Vec::new();
        for n in [8, 16, 32] {
            let state = birth(n, support, &law);
            let (next, report) = state.advance_supported_controlled(1.0, &temperatures(&state, 300.0),
                9.81, &vec![0.0; n as usize], &[stop(x, 0.0)], beam_config(), config(), &mut || false).unwrap();
            let force = report.reactions[0].upward_force_n;
            errors.push((force/weight-fraction).abs());
            near(next.deformed_position_m(x, 0.0).unwrap()[1], 0.0, 1e-12);
            near(report.beam.vertical_reactions_n.iter().sum::<f64>()+force, weight, 1e-12);
            near(report.beam.load_work_residual_j, 0.0, 1e-12);
            near(report.beam.energy_j.residual, 0.0, 1e-12);
        }
        assert!(errors[0]/errors[1] > 3.5 && errors[1]/errors[2] > 3.5, "{errors:?}");
    }
}

// G0: supported admission must not require an inadmissible unsupported material trial.
#[test]
fn supports_reduce_real_sag_and_removal_restores_the_free_solution() {
    let law = material(1e8, vec![]);
    let state = birth(16, Support::CantileverLeft, &law);
    let t = temperatures(&state, 300.0);
    let cfg = BeamConfig { maximum_deflection_m: 1e-5, ..beam_config() };
    assert!(state.advance_controlled(1.0, &t, 9.81, &[0.0;16], cfg, &mut || false).is_err());
    let (supported, r) = state.advance_supported_controlled(1.0, &t, 9.81, &[0.0;16],
        &[stop(0.02, 0.0)], cfg, config(), &mut || false).unwrap();
    assert!(r.beam.maximum_sag_m < cfg.maximum_deflection_m);
    assert!(r.reactions[0].upward_force_n > 0.0);
    let (released, rr) = supported.advance_supported_controlled(2.0, &t, 9.81, &[0.0;16],
        &[], beam_config(), config(), &mut || false).unwrap();
    let free = state.advance_controlled(2.0, &t, 9.81, &[0.0;16], beam_config(), &mut || false).unwrap();
    assert_eq!(released, free.0);
    near(rr.beam.maximum_sag_m, free.1.maximum_sag_m, 1e-14);
}

#[test]
fn unloading_can_lift_off_without_an_adhesive_or_negative_support_force() {
    let law = material(1e8, vec![]);
    let state = birth(8, Support::SimplySupported, &law);
    let t = temperatures(&state, 300.0);
    let supports = [stop(0.01, 0.0)];
    let (held, _) = state.advance_supported_controlled(1.0, &t, 9.81, &[0.0;8],
        &supports, beam_config(), config(), &mut || false).unwrap();
    let (lifted, r) = held.advance_supported_controlled(2.0, &t, 0.0, &[-0.1;8],
        &supports, beam_config(), config(), &mut || false).unwrap();
    assert_eq!(r.reactions[0].upward_force_n, 0.0);
    assert!(r.reactions[0].gap_m > 0.0);
    assert!(lifted.deformed_position_m(0.01, 0.0).unwrap()[1] > 0.0);
    assert_eq!(r.reaction_work_j, 0.0);
}

#[test]
fn multiple_stops_retain_force_balance_moving_support_work_and_order_invariance() {
    let law = material(1e8, vec![]);
    let state = birth(12, Support::SimplySupported, &law);
    let t = temperatures(&state, 300.0);
    let a = stop(0.006, -1e-6);
    let b = Stop { id: 2, ..stop(0.014, -0.5e-6) };
    let far = Stop { id: 3, ..stop(0.01, -0.001) };
    let solve = |stops: &[Stop]| state.advance_supported_controlled(1.0, &t, 9.81,
        &[0.0;12], stops, beam_config(), config(), &mut || false).unwrap();
    let (held, r) = solve(&[b, far, a]);
    let (ordered, ro) = solve(&[a, b, far]);
    assert_eq!(held, ordered); assert_eq!(r.reactions, ro.reactions);
    assert_eq!(r.reactions[2].upward_force_n, 0.0);
    for row in &r.reactions { assert!(row.gap_m >= -config().gap_tolerance_m); }
    let force: f64 = r.reactions.iter().map(|row| row.upward_force_n).sum();
    near(force+r.beam.vertical_reactions_n.iter().sum::<f64>(), state.mass_kg()*9.81, 1e-12);
    let (_, moved) = held.advance_supported_controlled(2.0, &t, 9.81, &[0.0;12],
        &[Stop { minimum_deflection_m: 0.0, ..a }, Stop { minimum_deflection_m: 0.0, ..b }],
        beam_config(), config(), &mut || false).unwrap();
    assert!(moved.reaction_work_j > 0.0);
    near(moved.beam.load_work_residual_j, 0.0, 1e-12);
}

// G1/G5: retained Maxwell history can bring a previously clear stop into contact.
#[test]
fn creep_engages_a_stop_without_resetting_material_memory_and_replays() {
    let law = material(1e7, vec![(9e7, 0.2)]);
    let mut state = birth(8, Support::SimplySupported, &law);
    let t = temperatures(&state, 300.0);
    let supports = [stop(0.01, -5e-6)];
    let mut last_force = 0.0;
    for i in 1..=10 {
        let (next, r) = state.advance_supported_controlled(f64::from(i)*0.1, &t, 9.81,
            &[0.0;8], &supports, beam_config(), config(), &mut || false).unwrap();
        if i == 1 { assert_eq!(r.reactions[0].upward_force_n, 0.0); }
        last_force = r.reactions[0].upward_force_n;
        state = next;
    }
    assert!(last_force > 0.0);
    assert!(state.cells()[0].section().layers()[&7].material_points()[0].branch_strains()[0] != 0.0);
    let checkpoint = state.clone();
    let expected = state.advance_supported_controlled(1.1, &t, 9.81, &[0.0;8],
        &supports, beam_config(), config(), &mut || false).unwrap();
    let actual = checkpoint.advance_supported_controlled(1.1, &t, 9.81, &[0.0;8],
        &supports, beam_config(), config(), &mut || false).unwrap();
    assert_eq!(expected.0, actual.0);
    assert_eq!(expected.1.reactions, actual.1.reactions);
}

// G4: no contact or section probe publishes part of a failed physical increment.
#[test]
fn admission_exhaustion_and_late_cancellation_leave_all_histories_unchanged() {
    let law = material(1e8, vec![]);
    let state = birth(8, Support::SimplySupported, &law);
    let before = state.clone();
    let t = temperatures(&state, 300.0);
    let solve = |stops: &[Stop], cfg| state.advance_supported_controlled(1.0, &t, 9.81,
        &[0.0;8], stops, beam_config(), cfg, &mut || false);
    assert!(solve(&[stop(0.01, 0.0), stop(0.012, 0.0)], config()).is_err());
    assert!(solve(&[stop(0.0, 0.0)], config()).is_err());
    assert!(solve(&[stop(f64::NAN, 0.0)], config()).is_err());
    let cfg = Config { contact: GramOrthantConfig { max_setup_products: 7, ..config().contact }, ..config() };
    assert!(solve(&[stop(0.01, 0.0)], cfg).is_err());
    let cfg = Config { contact: GramOrthantConfig { max_sweeps: 1, ..config().contact }, ..config() };
    assert!(matches!(solve(&[stop(0.009, 0.0), Stop { id: 2, ..stop(0.011, 0.0) }], cfg),
        Err(Error::Contact(GramOrthantError::NotConverged { .. }))));
    let mut total_polls = 0;
    state.advance_supported_controlled(1.0, &t, 9.81, &[0.0;8],
        &[stop(0.01, 0.0)], beam_config(), config(), &mut || { total_polls += 1; false }).unwrap();
    for limit in [0, total_polls/3, 2*total_polls/3, total_polls-1] {
        let mut calls = 0;
        assert!(state.advance_supported_controlled(1.0, &t, 9.81, &[0.0;8],
            &[stop(0.01, 0.0)], beam_config(), config(), &mut || { calls += 1; calls > limit }).is_err());
    }
    assert_eq!(state, before);
}

#[test]
fn condensed_section_response_matches_the_actual_thermal_history_update() {
    let law = material(1e8, vec![(2e8, 0.5)]);
    let state = birth(8, Support::SimplySupported, &law)
        .deposit_bonded_layer(8, crate::deposition::LayerGeometry { area_m2: 0.8e-6,
            centroid_z_m: 0.001, thickness_m: 0.001 }, 1200.0, &law,
            &(0..8).map(|id| (id, 450.0)).collect(), &mut || false).unwrap();
    let old = state.advance_controlled(0.5, &temperatures(&state, 350.0), 9.81,
        &[0.0;8], beam_config(), &mut || false).unwrap().0;
    let section = old.cells()[3].section();
    let before = section.clone();
    let t = BTreeMap::from([(7, 320.0), (8, 390.0)]);
    let response = section.bending_response(1.0, &t, || Ok(())).unwrap();
    for m in [-1e-6, 0.0, 1e-6] {
        let control = crate::deposition::Control { axial: crate::deposition::AxialControl::ForceN(0.0),
            bending: crate::deposition::BendingControl::MomentNm(m) };
        let (next, _) = section.advance_to(1.0, &t, control, beam_config().section, || Ok(())).unwrap();
        near(next.curvature_per_m(), response.zero_moment_curvature_per_m+response.compliance_per_n_m2*m, 1e-12);
    }
    assert_eq!(*section, before);
}
