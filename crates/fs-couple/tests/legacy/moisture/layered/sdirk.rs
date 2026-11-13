//! Independent physical exchange and storage checks for layered SDIRK2.
use super::*;

const GAMMA: f64 = 1.0 - core::f64::consts::FRAC_1_SQRT_2;

fn stage_moisture(report: &MoistureLayeredReport, cell: usize) -> f64 {
    report.transport_end[cell].query().axes()["moisture_content"]
}

#[test]
fn layered_sdirk_converges_to_exact_partition_exchange_with_finite_contact() {
    with_cx(|cx| {
        let rho = [400.0, 900.0];
        let slope = [0.2, 0.5];
        let length: [f64; 2] = [0.08, 0.18];
        let diffusion = [0.0002, 0.00005];
        let initial = [0.8, 0.15];
        let area = 0.001;
        let h = 0.04;
        let mass: [f64; 2] = std::array::from_fn(|i| rho[i] * area * length[i]);
        let capacity: [f64; 2] = std::array::from_fn(|i| mass[i] * slope[i]);
        let half: [f64; 2] =
            std::array::from_fn(|i| 2.0 * rho[i] * area * diffusion[i] * slope[i] / length[i]);
        let conductance = 1.0 / (1.0 / half[0] + 1.0 / (h * area) + 1.0 / half[1]);
        let rate = conductance * (1.0 / capacity[0] + 1.0 / capacity[1]);
        let duration = 1.0 / rate;
        let mean =
            (capacity[0] * initial[0] + capacity[1] * initial[1]) / (capacity[0] + capacity[1]);
        let exact_difference = (initial[0] - initial[1]) * (-rate * duration).exp();
        let exact = [
            mean + capacity[1] / (capacity[0] + capacity[1]) * exact_difference,
            mean - capacity[0] / (capacity[0] + capacity[1]) * exact_difference,
        ];
        let original = MoistureLayeredSlab::new(
            (0..2)
                .map(|i| {
                    layer(
                        material(rho[i], vec![(0.0, 0.0), (1.0, slope[i])], [diffusion[i]; 2]),
                        length[i],
                        vec![slope[i] * initial[i]],
                        cx,
                    )
                })
                .collect(),
            2,
            cx,
        )
        .unwrap()
        .with_transfer(
            [MoistureActivityTransfer::Sealed; 2],
            vec![MoistureActivityTransfer::ConductanceKgM2S(h)],
            cx,
        )
        .unwrap();
        let mut errors = Vec::new();
        for count in [4, 8, 16] {
            let dt = duration / f64::from(count);
            let mut state = original.clone();
            for _ in 0..count {
                let before = state.clone();
                let (next, report) = state
                    .advance_sdirk2(dt, &sealed(), config(), nonlinear(), cx)
                    .unwrap();
                let first: [f64; 2] =
                    std::array::from_fn(|i| stage_moisture(&report.first_stage, i) / slope[i]);
                let first_rate = conductance * (first[0] - first[1]);
                let last_rate = conductance * (next.activity()[0] - next.activity()[1]);
                let transfer = dt * ((1.0 - GAMMA) * first_rate + GAMMA * last_rate);
                near(
                    report.first_stage.interfaces[0].water_left_to_right_kg,
                    GAMMA * dt * first_rate,
                    2e-12,
                );
                near(
                    report.endpoint.interfaces[0].water_left_to_right_kg,
                    transfer,
                    2e-12,
                );
                near(
                    next.layers()[0].slab().water_mass_kg()
                        - before.layers()[0].slab().water_mass_kg(),
                    -transfer,
                    2e-12,
                );
                near(
                    next.layers()[1].slab().water_mass_kg()
                        - before.layers()[1].slab().water_mass_kg(),
                    transfer,
                    2e-12,
                );
                near(
                    report.endpoint.interfaces[0].relative_activity,
                    next.activity()[0] - last_rate / half[0],
                    2e-12,
                );
                near(
                    report.endpoint.interfaces[0].right_relative_activity,
                    next.activity()[1] + last_rate / half[1],
                    2e-12,
                );
                near(next.water_mass_kg(), original.water_mass_kg(), 2e-12);
                assert_eq!(report.endpoint.boundary_water_kg, [0.0; 2]);
                near(report.endpoint.interval_s[0], before.time_s(), 0.0);
                near(report.endpoint.interval_s[1], before.time_s() + dt, 0.0);
                near(
                    report.first_stage.interval_s[1],
                    before.time_s() + GAMMA * dt,
                    0.0,
                );
                state = next;
            }
            errors.push(
                (state.activity()[0] - exact[0])
                    .abs()
                    .max((state.activity()[1] - exact[1]).abs()),
            );
        }
        for pair in errors.windows(2) {
            assert!(
                pair[0] / pair[1] > 3.8 && pair[0] / pair[1] < 4.3,
                "{errors:?}"
            );
        }
        assert_eq!(original.time_s(), 0.0);
    });
}

#[test]
fn layered_sdirk_advances_actual_water_across_a_sorption_capacity_kink() {
    with_cx(|cx| {
        // Constant D makes this single actual cell obey M'=k*(Meq-M),
        // even though M(a) changes slope sixfold between the two stages.
        let old_m = 0.02;
        let equilibrium_m = 0.28;
        let length: f64 = 0.5;
        let diffusion = 0.02;
        let dt = 0.2 / (2.0 * diffusion / length.powi(2));
        let original = MoistureLayeredSlab::new(
            vec![layer(
                material(
                    500.0,
                    vec![(0.0, 0.0), (0.4, 0.04), (1.0, 0.4)],
                    [diffusion; 2],
                ),
                length,
                vec![old_m],
                cx,
            )],
            1,
            cx,
        )
        .unwrap();
        let faces = [
            MoistureLayeredFace::Sealed,
            MoistureLayeredFace::PrescribedRelativeHumidity(0.8),
        ];
        let controls = MoistureNonlinearConfig {
            relaxation: 0.5,
            maximum_iterations: 256,
            ..nonlinear()
        };
        let (next, report) = original
            .advance_sdirk2(dt, &faces, config(), controls, cx)
            .unwrap();
        let z = 0.2;
        let first_m = (old_m + GAMMA * z * equilibrium_m) / (1.0 + GAMMA * z);
        let amplification = (1.0 - (1.0 - 2.0 * GAMMA) * z) / (1.0 + GAMMA * z).powi(2);
        let final_m = equilibrium_m + amplification * (old_m - equilibrium_m);
        assert!(first_m < 0.04 && final_m > 0.04);
        near(stage_moisture(&report.first_stage, 0), first_m, 3e-12);
        near(next.layers()[0].slab().moisture()[0], final_m, 3e-12);
        near(next.activity()[0], 0.4 + (final_m - 0.04) / 0.6, 6e-12);
        near(
            report.endpoint.boundary_water_kg[1],
            original.dry_mass_kg() * (final_m - old_m),
            3e-12,
        );
        near(
            report.endpoint.stored_water_kg[0],
            original.water_mass_kg(),
            0.0,
        );
        near(
            report.endpoint.stored_water_kg[1],
            next.water_mass_kg(),
            0.0,
        );
        assert_eq!(report.endpoint.boundary_water_kg[0], 0.0);
        assert!(report.endpoint.nonlinear.relative_residual < controls.relative_tolerance);
    });
}

#[test]
fn nonlinear_layered_sdirk_weights_each_original_interface_and_surface_flux_once() {
    with_cx(|cx| {
        let rho = [500.0, 1200.0];
        let length = [0.1, 0.25];
        let area = 0.001;
        let h_surface = 0.02;
        let h_interface = 0.01;
        let sorption = [
            vec![(0.0, 0.0), (0.3, 0.03), (1.0, 0.31)],
            vec![(0.0, 0.0), (0.5, 0.2), (1.0, 0.5)],
        ];
        // Two different affine D(M) laws, integrated independently below.
        let intercept = [0.0, 0.0001];
        let gradient = [0.005, 0.002];
        let initial = [0.15, 0.04];
        let original = MoistureLayeredSlab::new(
            (0..2)
                .map(|i| {
                    layer(
                        material(
                            rho[i],
                            sorption[i].clone(),
                            [intercept[i], intercept[i] + 0.6 * gradient[i]],
                        ),
                        length[i],
                        vec![initial[i]],
                        cx,
                    )
                })
                .collect(),
            2,
            cx,
        )
        .unwrap()
        .with_transfer(
            [
                MoistureActivityTransfer::ConductanceKgM2S(h_surface),
                MoistureActivityTransfer::Sealed,
            ],
            vec![MoistureActivityTransfer::ConductanceKgM2S(h_interface)],
            cx,
        )
        .unwrap();
        let dt = 2.0;
        let faces = [
            MoistureLayeredFace::PrescribedRelativeHumidity(0.85),
            MoistureLayeredFace::Sealed,
        ];
        let controls = MoistureNonlinearConfig {
            relaxation: 0.5,
            ..nonlinear()
        };
        let (next, report) = original
            .advance_sdirk2(dt, &faces, config(), controls, cx)
            .unwrap();
        let rates = |stage: &MoistureLayeredReport| {
            let interface = &stage.interfaces[0];
            let m: [f64; 2] = std::array::from_fn(|i| stage_moisture(stage, i));
            let potential = |i: usize, x: f64| intercept[i] * x + 0.5 * gradient[i] * x * x;
            let left = 2.0 * rho[0] * area / length[0]
                * (potential(0, m[0]) - potential(0, interface.moisture[0]));
            let right = 2.0 * rho[1] * area / length[1]
                * (potential(1, interface.moisture[1]) - potential(1, m[1]));
            let film = h_interface
                * area
                * (interface.relative_activity - interface.right_relative_activity);
            near(left, right, 2e-13);
            near(left, film, 2e-13);
            let surface = stage.surface_activities[0].unwrap();
            let surface_m = original.layers()[0].sorption().moisture(surface).unwrap();
            let incoming = h_surface * area * (0.85 - surface);
            near(
                incoming,
                2.0 * rho[0] * area / length[0] * (potential(0, surface_m) - potential(0, m[0])),
                2e-13,
            );
            [incoming, film]
        };
        let first = rates(&report.first_stage);
        let last = rates(&report.endpoint);
        let weighted: [f64; 2] =
            std::array::from_fn(|i| dt * ((1.0 - GAMMA) * first[i] + GAMMA * last[i]));
        near(
            report.first_stage.boundary_water_kg[0],
            GAMMA * dt * first[0],
            2e-12,
        );
        near(
            report.first_stage.interfaces[0].water_left_to_right_kg,
            GAMMA * dt * first[1],
            2e-12,
        );
        near(report.endpoint.boundary_water_kg[0], weighted[0], 2e-12);
        near(
            report.endpoint.interfaces[0].water_left_to_right_kg,
            weighted[1],
            2e-12,
        );
        near(
            next.layers()[0].slab().water_mass_kg() - original.layers()[0].slab().water_mass_kg(),
            weighted[0] - weighted[1],
            2e-12,
        );
        near(
            next.layers()[1].slab().water_mass_kg() - original.layers()[1].slab().water_mass_kg(),
            weighted[1],
            2e-12,
        );
        near(
            next.water_mass_kg() - original.water_mass_kg(),
            weighted[0],
            2e-12,
        );
        assert_eq!(report.endpoint.boundary_water_kg[1], 0.0);
    });
}

#[test]
fn layered_sdirk_rejects_a_nonmonotone_stage_and_preserves_stiff_partition_equilibrium() {
    with_cx(|cx| {
        let original = MoistureLayeredSlab::new(
            vec![layer(
                material(500.0, vec![(0.0, 0.0), (1.0, 0.2)], [0.01; 2]),
                1.0,
                vec![0.1],
                cx,
            )],
            1,
            cx,
        )
        .unwrap();
        let faces = [
            MoistureLayeredFace::Sealed,
            MoistureLayeredFace::PrescribedRelativeHumidity(0.0),
        ];
        assert!(
            original
                .advance(GAMMA * 500.0, &faces, config(), nonlinear(), cx)
                .is_ok()
        );
        assert!(
            matches!(original.advance_sdirk2(500.0, &faces, config(), nonlinear(), cx), Err(MoistureError::ActivityStageRange { cell: 0, activity, bounds }) if activity < 0.0 && bounds == [0.0, 0.5])
        );
        assert_eq!(original.layers()[0].slab().moisture(), &[0.1]);
        assert_eq!(original.activity(), &[0.5]);
        assert_eq!(original.time_s(), 0.0);

        let equilibrium = MoistureLayeredSlab::new(
            vec![
                layer(
                    material(400.0, vec![(0.0, 0.0), (1.0, 0.25)], [1e100; 2]),
                    0.2,
                    vec![0.125; 2],
                    cx,
                ),
                layer(
                    material(900.0, vec![(0.0, 0.0), (1.0, 0.5)], [1e100; 2]),
                    0.3,
                    vec![0.25; 3],
                    cx,
                ),
            ],
            5,
            cx,
        )
        .unwrap();
        let (next, report) = equilibrium
            .advance_sdirk2(1e100, &sealed(), config(), nonlinear(), cx)
            .unwrap();
        assert_eq!(next.activity(), equilibrium.activity());
        assert_eq!(next.layers()[0].slab().moisture(), &[0.125; 2]);
        assert_eq!(next.layers()[1].slab().moisture(), &[0.25; 3]);
        assert_eq!(report.endpoint.boundary_water_kg, [0.0; 2]);
        assert_eq!(report.endpoint.interfaces[0].water_left_to_right_kg, 0.0);
        assert_eq!(report.first_stage.interfaces[0].water_left_to_right_kg, 0.0);
        near(next.water_mass_kg(), equilibrium.water_mass_kg(), 0.0);
    });
}
