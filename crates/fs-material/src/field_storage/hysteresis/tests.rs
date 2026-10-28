use super::*;

fn branch(chi: f64, hc: f64, saturation: f64) -> HysteresisBranch {
    HysteresisBranch {
        susceptibility_h_m: chi,
        pinning_field_a_m: hc,
        saturation_polarization_t: saturation,
    }
}

fn law() -> ScalarHysteresis {
    ScalarHysteresis::new(0.01, 30.0, vec![branch(0.2, 2.0, 1.0)]).unwrap()
}

fn close(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance,
        "{a:.17e} != {b:.17e}, tol={tolerance:e}"
    );
}

#[test]
fn analytic_pinning_slide_and_saturation_have_distinct_energy_and_heat() {
    let law = law();
    let initial = law.initial_state();
    let pinned = law.candidate_field(&initial, 1.0).unwrap();
    close(pinned.state.flux_density_t(), 0.01, 1e-15);
    assert_eq!(pinned.state.polarizations_t(), &[0.0]);
    assert_eq!(pinned.irreversible_heat_j_m3, 0.0);
    close(pinned.magnetic_work_j_m3, 0.005, 1e-15);

    let slide = law.candidate_field(&initial, 3.0).unwrap();
    close(slide.state.polarizations_t()[0], 0.2, 1e-15);
    close(slide.state.flux_density_t(), 0.23, 1e-15);
    close(slide.state.stored_energy_j_m3(), 0.145, 1e-15);
    close(slide.irreversible_heat_j_m3, 0.4, 1e-15);
    // Independent piecewise integral: mu*integral_0^2 H dH plus
    // (mu+chi)*integral_2^3 H dH. Endpoint trapezoids miss pinning.
    close(slide.magnetic_work_j_m3, 0.01 * 2.0 + 0.21 * 2.5, 1e-15);
    close(slide.incremental_permeability_h_m, 0.21, 1e-15);

    let saturated = law.candidate_field(&initial, 10.0).unwrap();
    assert_eq!(saturated.state.polarizations_t(), &[1.0]);
    close(saturated.state.flux_density_t(), 1.1, 1e-15);
    close(saturated.state.stored_energy_j_m3(), 3.0, 1e-15);
    close(saturated.irreversible_heat_j_m3, 2.0, 1e-15);
    close(saturated.incremental_permeability_h_m, 0.01, 1e-15);
}

#[test]
fn closed_major_loop_area_is_exact_pinning_heat() {
    let law = law();
    let top = law
        .candidate_field(&law.initial_state(), 10.0)
        .unwrap()
        .state;
    let down = law.candidate_field(&top, -10.0).unwrap();
    let up = law.candidate_field(&down.state, 10.0).unwrap();
    assert_eq!(up.state, top);
    close(
        down.irreversible_heat_j_m3 + up.irreversible_heat_j_m3,
        8.0,
        1e-14,
    );
    close(down.magnetic_work_j_m3 + up.magnetic_work_j_m3, 8.0, 1e-14);

    // Independent polygon at the analytical pinning and saturation corners.
    let corners = [
        (10.0, 1.1),
        (3.0, 1.03),
        (-7.0, -1.07),
        (-10.0, -1.1),
        (-3.0, -1.03),
        (7.0, 1.07),
        (10.0, 1.1),
    ];
    let area: f64 = corners
        .windows(2)
        .map(|p| 0.5 * (p[0].0 + p[1].0) * (p[1].1 - p[0].1))
        .sum();
    close(area, 8.0, 1e-14);
    close(area, down.magnetic_work_j_m3 + up.magnetic_work_j_m3, 1e-14);
}

#[test]
fn branch_memory_retains_remanence_and_closes_nested_minor_loops() {
    let law = ScalarHysteresis::new(
        0.01,
        30.0,
        vec![
            branch(0.2, 2.0, 1.0),
            branch(0.1, 4.0, 0.7),
            branch(0.05, 0.0, 0.4),
        ],
    )
    .unwrap();
    let high = law
        .candidate_field(&law.initial_state(), 15.0)
        .unwrap()
        .state;
    let remanent = law.candidate_field(&high, 0.0).unwrap();
    assert!(remanent.state.flux_density_t() > 0.0);
    assert!(remanent.irreversible_heat_j_m3 > 0.0);
    let low = law.candidate_field(&high, 1.0).unwrap().state;
    let middle = law.candidate_field(&low, 6.0).unwrap().state;
    let returned = law.candidate_field(&middle, 1.0).unwrap().state;
    assert_eq!(returned, low);
    let restored = law.candidate_field(&returned, 15.0).unwrap().state;
    assert_eq!(restored, high);
    assert_ne!(
        law.candidate_field(&law.initial_state(), 1.0)
            .unwrap()
            .state,
        low
    );
}

#[test]
fn monotone_partition_changes_neither_history_nor_heat_or_work() {
    let law = ScalarHysteresis::new(
        0.01,
        30.0,
        vec![branch(0.2, 2.0, 1.0), branch(0.1, 5.0, 0.7)],
    )
    .unwrap();
    let coarse = law.candidate_field(&law.initial_state(), 12.0).unwrap();
    let mut state = law.initial_state();
    let (mut heat, mut work) = (0.0, 0.0);
    for i in 1..=120 {
        let step = law.candidate_field(&state, f64::from(i) * 0.1).unwrap();
        heat += step.irreversible_heat_j_m3;
        work += step.magnetic_work_j_m3;
        state = step.state;
    }
    assert_eq!(state, coarse.state);
    close(heat, coarse.irreversible_heat_j_m3, 5e-14);
    close(work, coarse.magnetic_work_j_m3, 5e-14);
    for _ in 0..10 {
        let hold = law.candidate_field(&state, state.field_a_m()).unwrap();
        assert_eq!(hold.state, state);
        assert_eq!(hold.irreversible_heat_j_m3, 0.0);
        assert_eq!(hold.magnetic_work_j_m3, 0.0);
    }
}

#[test]
fn flux_inverse_recovers_the_same_reversal_history_and_heat() {
    let law = ScalarHysteresis::new(
        0.001,
        40.0,
        vec![
            branch(0.2, 2.0, 1.0),
            branch(0.04, 3.0, 0.3),
            branch(0.1, 0.0, 0.6),
        ],
    )
    .unwrap();
    let mut state = law.initial_state();
    for h in [7.0, 3.0, -4.0, 2.0, -10.0, 10.0, 0.0] {
        let direct = law.candidate_field(&state, h).unwrap();
        let inverse = law
            .candidate_flux_density(&state, direct.state.flux_density_t())
            .unwrap();
        close(inverse.state.field_a_m(), h, 1e-10);
        for (&a, &b) in inverse
            .state
            .polarizations_t()
            .iter()
            .zip(direct.state.polarizations_t())
        {
            close(a, b, 2e-12);
        }
        close(inverse.magnetic_work_j_m3, direct.magnetic_work_j_m3, 2e-11);
        close(
            inverse.irreversible_heat_j_m3,
            direct.irreversible_heat_j_m3,
            2e-11,
        );
        state = direct.state;
    }
}

#[test]
fn frozen_history_tangent_matches_independent_perturbations() {
    let law = law();
    let initial = law.initial_state();
    for h in [1.0, 4.0, 10.0, -1.0, -4.0, -10.0] {
        let point = law.candidate_field(&initial, h).unwrap();
        let plus = law
            .candidate_field(&initial, h + 1e-5)
            .unwrap()
            .state
            .flux_density_t();
        let minus = law
            .candidate_field(&initial, h - 1e-5)
            .unwrap()
            .state
            .flux_density_t();
        close(
            (plus - minus) / 2e-5,
            point.incremental_permeability_h_m,
            2e-10,
        );
    }
}

#[test]
fn zero_pinning_is_reversible_and_empty_branches_recover_linear_storage() {
    for law in [
        ScalarHysteresis::new(0.01, 30.0, vec![]).unwrap(),
        ScalarHysteresis::new(0.01, 30.0, vec![branch(0.2, 0.0, 1.0)]).unwrap(),
    ] {
        let mut state = law.initial_state();
        let mut work = 0.0;
        for h in [1.0, 10.0, -10.0, 0.0] {
            let step = law.candidate_field(&state, h).unwrap();
            assert_eq!(step.irreversible_heat_j_m3, 0.0);
            work += step.magnetic_work_j_m3;
            state = step.state;
        }
        close(work, 0.0, 1e-14);
        assert_eq!(state, law.initial_state());
    }
}

#[test]
fn inverse_does_not_erase_tiny_below_pinning_flux() {
    let law = ScalarHysteresis::new(1e-6, 1000.0, vec![branch(0.01, 30.0, 1.0)]).unwrap();
    for b in [1e-20, -1e-20] {
        let step = law.candidate_flux_density(&law.initial_state(), b).unwrap();
        close(step.state.flux_density_t(), b, 1e-34);
        close(step.state.field_a_m(), b / 1e-6, 1e-28);
        assert_eq!(step.irreversible_heat_j_m3, 0.0);
    }
}

#[test]
fn refusals_preserve_history_and_do_not_rebind_parameters() {
    let law = law();
    let state = law
        .candidate_field(&law.initial_state(), 10.0)
        .unwrap()
        .state;
    let before = state.clone();
    for h in [f64::NAN, f64::INFINITY, 30.01, -30.01] {
        assert!(law.candidate_field(&state, h).is_err());
    }
    for b in [f64::NAN, f64::INFINITY, 100.0, -100.0] {
        assert!(law.candidate_flux_density(&state, b).is_err());
    }
    let different = ScalarHysteresis::new(0.01, 30.0, vec![branch(0.2, 3.0, 1.0)]).unwrap();
    assert!(different.candidate_field(&state, 10.0).is_err());
    assert_eq!(state, before);
    for bad in [
        branch(0.0, 2.0, 1.0),
        branch(0.2, -1.0, 1.0),
        branch(0.2, 2.0, f64::NAN),
    ] {
        assert!(ScalarHysteresis::new(0.01, 30.0, vec![bad]).is_err());
    }
    assert!(ScalarHysteresis::new(0.01, 30.0, vec![branch(0.2, 2.0, 1.0); 1025]).is_err());
    assert!(ScalarHysteresis::new(0.0, 30.0, vec![]).is_err());
}

#[test]
fn uniform_core_links_current_flux_and_heat_with_signed_winding_orientation() {
    let positive = UniformHystereticCore::new(law(), 100.0, 1e-4, 0.2).unwrap();
    let negative = UniformHystereticCore::new(law(), -100.0, 1e-4, 0.2).unwrap();
    let driven = positive
        .candidate_current(&positive.initial_state(), 0.02)
        .unwrap();
    close(driven.flux_linkage_wb, 0.011, 1e-15);
    close(driven.current_a, 0.02, 1e-15);
    close(driven.core_heat_j, 2.0 * 2e-5, 1e-18);
    close(driven.stored_energy_j, 3.0 * 2e-5, 1e-18);
    close(
        driven.mean_current_a * driven.flux_linkage_wb,
        driven.magnetic_work_j,
        1e-18,
    );
    close(
        driven.magnetic_work_j,
        driven.stored_energy_j + driven.core_heat_j,
        1e-18,
    );
    let inverted = positive
        .candidate_linkage(&positive.initial_state(), driven.flux_linkage_wb)
        .unwrap();
    close(inverted.current_a, driven.current_a, 1e-12);
    close(inverted.core_heat_j, driven.core_heat_j, 1e-16);
    let reversed = negative
        .candidate_current(&negative.initial_state(), -0.02)
        .unwrap();
    close(reversed.flux_linkage_wb, -driven.flux_linkage_wb, 1e-15);
    close(reversed.current_a, -driven.current_a, 1e-15);
    close(reversed.core_heat_j, driven.core_heat_j, 1e-18);
    close(
        reversed.incremental_inductance_h,
        driven.incremental_inductance_h,
        1e-15,
    );
}

#[test]
fn core_loop_heat_is_the_analytic_saturation_area_times_actual_volume() {
    // Explicit synthetic SI coefficients; no empirical material validation claim.
    let law =
        ScalarHysteresis::new(1.25663706212e-6, 5000.0, vec![branch(0.005, 40.0, 1.6)]).unwrap();
    let core = UniformHystereticCore::new(law, 200.0, 2e-4, 0.3).unwrap();
    let top = core.candidate_current(&core.initial_state(), 1.5).unwrap();
    let down = core.candidate_current(&top.material.state, -1.5).unwrap();
    let up = core.candidate_current(&down.material.state, 1.5).unwrap();
    let analytic_heat = 4.0 * 40.0 * 1.6 * (2e-4 * 0.3);
    close(down.core_heat_j + up.core_heat_j, analytic_heat, 1e-16);
    close(
        down.magnetic_work_j + up.magnetic_work_j,
        analytic_heat,
        1e-16,
    );
    close(up.stored_energy_j, top.stored_energy_j, 1e-16);
    assert_eq!(up.material.state, top.material.state);
}
