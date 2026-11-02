//! G0/G1/G3/G4: internal generation with the original nonlinear material law.
use super::*;
use fs_conduction::layered::NonlinearGenerationOptions;

fn generated_options(interval: [f64; 2]) -> NonlinearGenerationOptions {
    NonlinearGenerationOptions {
        temperature_interval_k: interval,
        solver: options(),
    }
}

fn heated_affine() -> LayeredWall {
    LayeredWall::new(
        vec![curve(
            "k=T/300",
            2.0,
            vec![(280.0, 280.0 / 300.0), (800.0, 800.0 / 300.0)],
        )],
        vec![],
        [0.0; 2],
    )
    .unwrap()
}

#[test]
fn g1_heated_affine_law_recovers_the_nonlinear_interior_profile() {
    let wall = heated_affine();
    let result = with_cx(|cx| {
        wall.solve_temperature_dependent_with_generation(
            cx,
            [300.0; 2],
            2.0,
            &[600.0],
            generated_options([280.0, 800.0]),
        )
    })
    .unwrap();
    // Independently: T²=300²+600*(600*x-300*x²), from k=T/300.
    let row = &result.layer_temperatures[0];
    close(result.boundary_fluxes_w_per_m2[0], -600.0);
    close(result.boundary_fluxes_w_per_m2[1], 600.0);
    close(result.reservoir_heat_rates_w[0], 1200.0);
    close(result.reservoir_heat_rates_w[1], 1200.0);
    close(result.generated_heat_rate_w, 2400.0);
    close(result.energy_residual_w, 0.0);
    close(row.maximum_temperature_k, 270_000.0_f64.sqrt());
    close(row.maximum_position_m, 1.0);
    close(row.minimum_temperature_k, 300.0);
    for x in [0.0, 0.25, 0.5, 1.0, 1.75, 2.0] {
        close(
            with_cx(|cx| result.temperature_at(cx, 0, x)).unwrap(),
            (90_000.0 + 360_000.0 * x - 180_000.0 * x * x).sqrt(),
        );
    }
    assert!(
        row.maximum_temperature_k < 600.0,
        "freezing k at its boundary value would predict 600 K"
    );
    assert_eq!(result.wall, wall);
    assert_eq!(result.iterations, 1);
}

#[test]
fn g1_interior_source_knot_and_negative_generation_are_resolved() {
    let wall = LayeredWall::new(
        vec![curve(
            "kink",
            2.0,
            vec![(290.0, 1.0), (300.0, 1.0), (302.0, 3.0), (310.0, 3.0)],
        )],
        vec![],
        [0.0; 2],
    )
    .unwrap();
    let answer = with_cx(|cx| {
        wall.solve_temperature_dependent_with_generation(
            cx,
            [300.0; 2],
            1.0,
            &[20.0],
            generated_options([290.0, 310.0]),
        )
    })
    .unwrap();
    // Potential rise at the centre is s*L²/8=10: 4 below302 K and 6 above.
    close(answer.layer_temperatures[0].maximum_temperature_k, 304.0);
    close(
        with_cx(|cx| answer.temperature_at(cx, 0, 0.5)).unwrap(),
        302.0 + 3.5 / 3.0,
    );
    let cold = with_cx(|cx| {
        heated_affine().solve_temperature_dependent_with_generation(
            cx,
            [600.0; 2],
            1.0,
            &[-600.0],
            generated_options([280.0, 800.0]),
        )
    })
    .unwrap();
    close(
        cold.layer_temperatures[0].minimum_temperature_k,
        180_000.0_f64.sqrt(),
    );
    close(cold.layer_temperatures[0].minimum_position_m, 1.0);
    close(cold.reservoir_heat_rates_w[0], -600.0);
    close(cold.reservoir_heat_rates_w[1], -600.0);
    close(cold.energy_residual_w, 0.0);
}

fn mixed_generated() -> (LayeredWall, [f64; 2]) {
    let wall = LayeredWall::new(
        vec![
            curve(
                "nonlinear heater",
                1.0,
                vec![(280.0, 280.0 / 300.0), (400.0, 400.0 / 300.0)],
            ),
            PlanarThermalLayer::new(
                "constant sink",
                0.5,
                ConductivityTable::declared(2.0).unwrap(),
            )
            .unwrap(),
        ],
        vec![0.25],
        [0.5, 0.25],
    )
    .unwrap();
    // Independently construct from inlet q=-1: T_left=300.5, mean q=1,
    // K drop=1, q after layer=3; contact drop=.75; sink mean q=2.5
    // and drop=.625; right q=2, film drop=.5.
    (wall, [300.0, (300.5_f64.powi(2) - 600.0).sqrt() - 1.875])
}

#[test]
fn g1_g3_generated_stack_closes_local_contacts_films_and_signed_energy() {
    let (wall, reservoirs) = mixed_generated();
    let opts = generated_options([280.0, 400.0]);
    let answer = with_cx(|cx| {
        wall.solve_temperature_dependent_with_generation(cx, reservoirs, 4.0, &[4.0, -2.0], opts)
    })
    .unwrap();
    close(answer.boundary_fluxes_w_per_m2[0], -1.0);
    close(answer.boundary_fluxes_w_per_m2[1], 2.0);
    close(answer.reservoir_heat_rates_w[0], 4.0);
    close(answer.reservoir_heat_rates_w[1], 8.0);
    close(answer.generated_heat_rate_w, 12.0);
    close(answer.energy_residual_w, 0.0);
    let a = &answer.layer_temperatures[0];
    let b = &answer.layer_temperatures[1];
    close(a.face_temperatures_k[0], 300.5);
    close(a.face_temperatures_k[1], (300.5_f64.powi(2) - 600.0).sqrt());
    close(a.maximum_position_m, 0.25);
    close(a.maximum_temperature_k, (300.5_f64.powi(2) + 75.0).sqrt());
    close(a.face_temperatures_k[1] - b.face_temperatures_k[0], 0.75);
    close(b.face_temperatures_k[0] - b.face_temperatures_k[1], 0.625);
    assert!(answer.maximum_flux_residual_w_per_m2 <= 1.01e-8);
    assert!(answer.right_boundary_residual_k.abs() <= opts.solver.temperature_tolerance_k);
    assert!(answer.left_flux_bracket_w_per_m2[0] <= -1.0);
    assert!(answer.left_flux_bracket_w_per_m2[1] >= -1.0);
    let reflected = LayeredWall::new(
        wall.layers().iter().cloned().rev().collect(),
        vec![0.25],
        [0.25, 0.5],
    )
    .unwrap();
    let reverse = with_cx(|cx| {
        reflected.solve_temperature_dependent_with_generation(
            cx,
            [reservoirs[1], reservoirs[0]],
            4.0,
            &[-2.0, 4.0],
            opts,
        )
    })
    .unwrap();
    close(reverse.boundary_fluxes_w_per_m2[0], -2.0);
    close(reverse.boundary_fluxes_w_per_m2[1], 1.0);
    for (forward, backward) in answer
        .layer_temperatures
        .iter()
        .zip(reverse.layer_temperatures.iter().rev())
    {
        close(
            forward.face_temperatures_k[0],
            backward.face_temperatures_k[1],
        );
        close(
            forward.face_temperatures_k[1],
            backward.face_temperatures_k[0],
        );
        close(
            forward.maximum_temperature_k,
            backward.maximum_temperature_k,
        );
    }
    let large = with_cx(|cx| {
        wall.solve_temperature_dependent_with_generation(cx, reservoirs, 12.0, &[4.0, -2.0], opts)
    })
    .unwrap();
    assert_eq!(large.layer_temperatures, answer.layer_temperatures);
    close(
        large.generated_heat_rate_w,
        3.0 * answer.generated_heat_rate_w,
    );
    assert_eq!(
        answer,
        with_cx(|cx| wall.solve_temperature_dependent_with_generation(
            cx,
            reservoirs,
            4.0,
            &[4.0, -2.0],
            opts,
        ))
        .unwrap()
    );
}

#[test]
fn g3_generation_subdivision_and_zero_source_agree_with_existing_owners() {
    let wall = heated_affine();
    let options = generated_options([280.0, 800.0]);
    let split = LayeredWall::new(
        vec![
            curve(
                "left",
                0.5,
                wall.layers()[0].conductivity().knots().to_vec(),
            ),
            curve(
                "right",
                1.5,
                wall.layers()[0].conductivity().knots().to_vec(),
            ),
        ],
        vec![0.0],
        [0.0; 2],
    )
    .unwrap();
    let answer = with_cx(|cx| {
        split.solve_temperature_dependent_with_generation(
            cx,
            [300.0; 2],
            2.0,
            &[600.0, 600.0],
            options,
        )
    })
    .unwrap();
    close(answer.boundary_fluxes_w_per_m2[0], -600.0);
    close(
        answer.layer_temperatures[1].maximum_temperature_k,
        270_000.0_f64.sqrt(),
    );
    close(answer.layer_temperatures[1].maximum_position_m, 0.5);
    let heated_zero = with_cx(|cx| {
        wall.solve_temperature_dependent_with_generation(cx, [600.0, 300.0], 2.0, &[0.0], options)
    })
    .unwrap();
    let existing =
        with_cx(|cx| wall.solve_temperature_dependent(cx, [600.0, 300.0], 2.0, options.solver))
            .unwrap();
    close(
        heated_zero.boundary_fluxes_w_per_m2[0],
        existing.heat_flux_w_per_m2,
    );
    close(
        heated_zero.boundary_fluxes_w_per_m2[1],
        existing.heat_flux_w_per_m2,
    );
    let constant = LayeredWall::new(
        vec![
            PlanarThermalLayer::new("constant", 2.0, ConductivityTable::declared(4.0).unwrap())
                .unwrap(),
        ],
        vec![],
        [0.5, 0.25],
    )
    .unwrap();
    let analytic =
        with_cx(|cx| constant.solve_with_generation(cx, [300.0, 310.0], 3.0, &[8.0])).unwrap();
    let nonlinear = with_cx(|cx| {
        constant.solve_temperature_dependent_with_generation(
            cx,
            [300.0, 310.0],
            3.0,
            &[8.0],
            generated_options([280.0, 400.0]),
        )
    })
    .unwrap();
    for i in 0..2 {
        close(
            nonlinear.boundary_fluxes_w_per_m2[i],
            analytic.boundary_fluxes_w_per_m2[i],
        );
    }
    close(
        nonlinear.layer_temperatures[0].maximum_temperature_k,
        analytic.layer_temperatures[0].maximum_temperature_k,
    );
}

#[test]
fn g0_g4_generation_refuses_hidden_extrema_unsupported_domains_and_unfunded_work() {
    let wall = LayeredWall::new(
        vec![curve(
            "kink",
            2.0,
            vec![(290.0, 1.0), (300.0, 1.0), (302.0, 3.0), (310.0, 3.0)],
        )],
        vec![],
        [0.0; 2],
    )
    .unwrap();
    // Both prescribed faces are supported, but the actual peak is304 K.
    assert!(
        with_cx(|cx| wall.solve_temperature_dependent_with_generation(
            cx,
            [300.0; 2],
            1.0,
            &[20.0],
            generated_options([290.0, 303.0]),
        ))
        .is_err()
    );
    let (mixed, reservoirs) = mixed_generated();
    let mut short = generated_options([280.0, 400.0]);
    short.solver.max_iterations = 1;
    assert!(
        with_cx(|cx| mixed.solve_temperature_dependent_with_generation(
            cx,
            reservoirs,
            1.0,
            &[4.0, -2.0],
            short,
        ))
        .is_err()
    );
    for interval in [
        [280.0, 310.0],
        [300.0, 300.0],
        [310.0, 300.0],
        [0.0, 310.0],
        [290.0, f64::INFINITY],
    ] {
        assert!(
            with_cx(|cx| wall.solve_temperature_dependent_with_generation(
                cx,
                [300.0; 2],
                1.0,
                &[20.0],
                generated_options(interval),
            ))
            .is_err()
        );
    }
    for source in [vec![], vec![f64::NAN], vec![f64::INFINITY], vec![1.0, 2.0]] {
        assert!(
            with_cx(|cx| wall.solve_temperature_dependent_with_generation(
                cx,
                [300.0; 2],
                1.0,
                &source,
                generated_options([290.0, 310.0]),
            ))
            .is_err()
        );
    }
    assert!(matches!(
        with_cancelled_cx(|cx| wall.solve_temperature_dependent_with_generation(
            cx,
            [300.0; 2],
            1.0,
            &[20.0],
            generated_options([290.0, 310.0]),
        )),
        Err(ConductionError::Cancelled { .. })
    ));
    let solved = with_cx(|cx| {
        wall.solve_temperature_dependent_with_generation(
            cx,
            [300.0; 2],
            1.0,
            &[20.0],
            generated_options([290.0, 310.0]),
        )
    })
    .unwrap();
    assert!(with_cx(|cx| solved.temperature_at(cx, 1, 0.0)).is_err());
    assert!(with_cx(|cx| solved.temperature_at(cx, 0, -0.1)).is_err());
    assert!(with_cx(|cx| solved.temperature_at(cx, 0, f64::NAN)).is_err());
    assert!(matches!(
        with_cancelled_cx(|cx| solved.temperature_at(cx, 0, 1.0)),
        Err(ConductionError::Cancelled { .. })
    ));
}
