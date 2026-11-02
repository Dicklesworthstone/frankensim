//! G1 exact internally heated slabs, G3 conservation/reversal, G4 refusal.
use super::*;

fn symmetric_stack() -> LayeredWall {
    LayeredWall::new(vec![constant("heated slab", 2.0, 4.0)], vec![], [0.0; 2]).unwrap()
}

fn asymmetric_stack() -> LayeredWall {
    LayeredWall::new(
        vec![
            constant("heated", 1.0, 1.0),
            constant("conductor", 1.0, 2.0),
        ],
        vec![0.5],
        [0.5, 1.0],
    )
    .unwrap()
}

#[test]
fn g1_generated_slab_resolves_the_interior_peak_and_two_outward_heat_rates() {
    let stack = symmetric_stack();
    let answer = with_cx(|cx| stack.solve_with_generation(cx, [300.0; 2], 3.0, &[8.0])).unwrap();
    // Independent closed form: T(x)=300+2x-x^2, q(x)=-8+8x.
    let layer = &answer.layer_temperatures[0];
    assert_eq!(layer.face_temperatures_k, [300.0; 2]);
    assert_eq!(layer.face_fluxes_w_per_m2, [-8.0, 8.0]);
    close(layer.maximum_temperature_k, 301.0);
    close(layer.maximum_position_m, 1.0);
    close(layer.minimum_temperature_k, 300.0);
    close(layer.mean_temperature_k, 300.0 + 2.0 / 3.0);
    assert_eq!(answer.reservoir_heat_rates_w, [24.0; 2]);
    close(answer.generated_heat_rate_w, 48.0);
    close(answer.energy_residual_w, 0.0);
    close(answer.maximum_flux_residual_w_per_m2, 0.0);
    close(answer.right_boundary_residual_k, 0.0);
    assert_eq!(answer.wall, stack);
}

#[test]
fn g1_heterogeneous_generation_uses_the_local_flux_at_each_contact_and_film() {
    let answer =
        with_cx(|cx| asymmetric_stack().solve_with_generation(cx, [300.0; 2], 7.0, &[2.0, 0.0]))
            .unwrap();
    // Independently solved series balance: R=7/2, generation drop=5 K.
    // q_left=-10/7, q_right=4/7. The first layer's q vanishes at x=5/7.
    close(answer.boundary_fluxes_w_per_m2[0], -10.0 / 7.0);
    close(answer.boundary_fluxes_w_per_m2[1], 4.0 / 7.0);
    close(answer.reservoir_heat_rates_w[0], 10.0);
    close(answer.reservoir_heat_rates_w[1], 4.0);
    close(answer.generated_heat_rate_w, 14.0);
    let first = &answer.layer_temperatures[0];
    let second = &answer.layer_temperatures[1];
    close(first.face_temperatures_k[0], 300.0 + 5.0 / 7.0);
    close(first.face_temperatures_k[1], 300.0 + 8.0 / 7.0);
    close(second.face_temperatures_k[0], 300.0 + 6.0 / 7.0);
    close(second.face_temperatures_k[1], 300.0 + 4.0 / 7.0);
    close(first.maximum_position_m, 5.0 / 7.0);
    close(first.maximum_temperature_k, 300.0 + 60.0 / 49.0);
    close(
        (first.face_temperatures_k[1] - second.face_temperatures_k[0]) / 0.5,
        4.0 / 7.0,
    );
    close(answer.energy_residual_w, 0.0);
    close(answer.maximum_flux_residual_w_per_m2, 0.0);
    close(answer.right_boundary_residual_k, 0.0);
}

#[test]
fn g1_signed_sources_resolve_internal_minima_and_zero_net_generation() {
    let removal =
        with_cx(|cx| symmetric_stack().solve_with_generation(cx, [300.0; 2], 3.0, &[-8.0]))
            .unwrap();
    let layer = &removal.layer_temperatures[0];
    close(layer.minimum_temperature_k, 299.0);
    close(layer.minimum_position_m, 1.0);
    close(layer.maximum_temperature_k, 300.0);
    close(layer.mean_temperature_k, 300.0 - 2.0 / 3.0);
    assert_eq!(removal.reservoir_heat_rates_w, [-24.0; 2]);
    close(removal.generated_heat_rate_w, -48.0);

    let stack = LayeredWall::new(
        vec![constant("source", 1.0, 1.0), constant("sink", 1.0, 1.0)],
        vec![0.0],
        [0.0; 2],
    )
    .unwrap();
    let answer =
        with_cx(|cx| stack.solve_with_generation(cx, [300.0; 2], 1.0, &[2.0, -2.0])).unwrap();
    // Zero net heat does not imply constant temperatures or zero local heat.
    assert_eq!(answer.generated_heat_rate_w, 0.0);
    assert_eq!(answer.reservoir_heat_rates_w, [1.0, -1.0]);
    close(answer.layer_temperatures[0].maximum_temperature_k, 300.25);
    close(answer.layer_temperatures[1].minimum_temperature_k, 299.75);
    close(answer.energy_residual_w, 0.0);
}

#[test]
fn g3_generation_free_path_matches_the_original_constant_conduction_owner() {
    let stack = wall();
    for temperatures in [[340.0, 280.0], [280.0, 340.0], [310.0, 310.0]] {
        let original = with_cx(|cx| stack.solve(cx, temperatures, 2.5)).unwrap();
        let answer =
            with_cx(|cx| stack.solve_with_generation(cx, temperatures, 2.5, &[0.0; 2])).unwrap();
        for flux in answer.boundary_fluxes_w_per_m2 {
            close(flux, original.heat_flux_w_per_m2);
        }
        close(answer.reservoir_heat_rates_w[0], -original.heat_rate_w);
        close(answer.reservoir_heat_rates_w[1], original.heat_rate_w);
        for (new, old) in answer
            .layer_temperatures
            .iter()
            .zip(&original.layer_temperatures)
        {
            close(new.face_temperatures_k[0], old.left_k);
            close(new.face_temperatures_k[1], old.right_k);
            close(new.mean_temperature_k, 0.5 * (old.left_k + old.right_k));
        }
    }
}

#[test]
fn g3_generated_stack_reversal_area_scaling_and_subdivision_preserve_physics() {
    let stack = asymmetric_stack();
    let original =
        with_cx(|cx| stack.solve_with_generation(cx, [310.0, 300.0], 7.0, &[2.0, 0.0])).unwrap();
    let reversed = LayeredWall::new(
        stack.layers().iter().rev().cloned().collect(),
        vec![0.5],
        [1.0, 0.5],
    )
    .unwrap();
    let backward =
        with_cx(|cx| reversed.solve_with_generation(cx, [300.0, 310.0], 7.0, &[0.0, 2.0])).unwrap();
    close(
        backward.reservoir_heat_rates_w[0],
        original.reservoir_heat_rates_w[1],
    );
    close(
        backward.reservoir_heat_rates_w[1],
        original.reservoir_heat_rates_w[0],
    );
    for (a, b) in original
        .layer_temperatures
        .iter()
        .zip(backward.layer_temperatures.iter().rev())
    {
        close(a.face_temperatures_k[0], b.face_temperatures_k[1]);
        close(a.face_temperatures_k[1], b.face_temperatures_k[0]);
        close(a.maximum_temperature_k, b.maximum_temperature_k);
        close(a.minimum_temperature_k, b.minimum_temperature_k);
        close(a.mean_temperature_k, b.mean_temperature_k);
    }
    let larger =
        with_cx(|cx| stack.solve_with_generation(cx, [310.0, 300.0], 28.0, &[2.0, 0.0])).unwrap();
    assert_eq!(larger.layer_temperatures, original.layer_temperatures);
    for side in 0..2 {
        close(
            larger.reservoir_heat_rates_w[side],
            4.0 * original.reservoir_heat_rates_w[side],
        );
    }

    let split = LayeredWall::new(
        vec![constant("left", 0.5, 4.0), constant("right", 1.5, 4.0)],
        vec![0.0],
        [0.0; 2],
    )
    .unwrap();
    let answer = with_cx(|cx| split.solve_with_generation(cx, [300.0; 2], 3.0, &[8.0; 2])).unwrap();
    assert_eq!(answer.reservoir_heat_rates_w, [24.0; 2]);
    close(answer.layer_temperatures[0].face_temperatures_k[1], 300.75);
    close(answer.layer_temperatures[1].maximum_temperature_k, 301.0);
    close(answer.layer_temperatures[1].maximum_position_m, 0.5);
    close(
        (0.5 * answer.layer_temperatures[0].mean_temperature_k
            + 1.5 * answer.layer_temperatures[1].mean_temperature_k)
            / 2.0,
        300.0 + 2.0 / 3.0,
    );
}

#[test]
fn g0_generation_retains_original_card_and_checks_interior_source_support() {
    let card = source_card(0.2);
    let stack = LayeredWall::new(vec![sourced(&card)], vec![], [0.0; 2]).unwrap();
    let result =
        with_cx(|cx| stack.solve_with_generation(cx, [300.0, 320.0], 1.5, &[2000.0])).unwrap();
    assert_eq!(result.wall, stack);
    close(result.generated_heat_rate_w, 300.0);
    close(result.layer_temperatures[0].maximum_temperature_k, 324.5);
    for receipt in result.wall.layers()[0].conductivity().receipts() {
        card.claims().verify_receipt(receipt).unwrap();
    }
    let bounded = PlanarThermalLayer::new(
        "bounded",
        2.0,
        ConductivityTable::declared_curve(vec![(299.5, 4.0), (300.5, 4.0)]).unwrap(),
    )
    .unwrap();
    let wall = LayeredWall::new(vec![bounded], vec![], [0.0; 2]).unwrap();
    // Both prescribed faces are supported, but the interior peak/trough is not.
    for (source, expected) in [(8.0, 301.0), (-8.0, 299.0)] {
        let error =
            with_cx(|cx| wall.solve_with_generation(cx, [300.0; 2], 1.0, &[source])).unwrap_err();
        assert!(
            matches!(error, ConductionError::OutsideTemperatureSpan { temperature, .. } if temperature == expected)
        );
    }
}

#[test]
fn g0_g4_generation_refuses_invalid_sources_nonlinear_coefficients_and_cancellation() {
    let stack = symmetric_stack();
    for sources in [vec![], vec![1.0, 2.0], vec![f64::NAN], vec![f64::INFINITY]] {
        assert!(with_cx(|cx| stack.solve_with_generation(cx, [300.0; 2], 1.0, &sources)).is_err());
    }
    assert!(with_cx(|cx| stack.solve_with_generation(cx, [300.0; 2], 0.0, &[8.0])).is_err());
    assert!(with_cx(|cx| stack.solve_with_generation(cx, [300.0; 2], 1.0, &[-3000.0])).is_err());
    let nonlinear = PlanarThermalLayer::new(
        "nonlinear",
        1.0,
        ConductivityTable::declared_curve(vec![(280.0, 1.0), (400.0, 2.0)]).unwrap(),
    )
    .unwrap();
    let nonlinear = LayeredWall::new(vec![nonlinear], vec![], [0.0; 2]).unwrap();
    assert!(matches!(
        with_cx(|cx| nonlinear.solve_with_generation(cx, [300.0; 2], 1.0, &[1.0])),
        Err(ConductionError::Conductivity { .. })
    ));
    assert!(matches!(
        with_cancelled_cx(|cx| stack.solve_with_generation(cx, [300.0; 2], 1.0, &[8.0])),
        Err(ConductionError::Cancelled { .. })
    ));
    assert_eq!(stack, symmetric_stack());
}
