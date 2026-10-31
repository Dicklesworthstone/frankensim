//! G1 nonlinear endpoint equation and G3 transport through existing consumers.
use super::*;

fn nonlinear() -> MoistureNonlinearConfig {
    MoistureNonlinearConfig {
        relative_tolerance: 1e-12,
        maximum_iterations: 128,
        relaxation: 1.0,
    }
}

#[test]
fn g1_nonlinear_one_cell_uses_endpoint_diffusivity_and_endpoint_water_flux() {
    with_cx(|cx| {
        let lagged = slab(card([0.01, 0.04], 0.3), 1.0, vec![0.1], cx);
        let implicit = lagged
            .clone()
            .with_implicit_diffusivity(nonlinear())
            .unwrap();
        let (next, report) = implicit.advance(10.0, &faces(0.2), config(), cx).unwrap();
        // M-.1=40*(.01+.1*M)*(.2-M), hence 4*M²+.6*M-.18=0.
        // The physical positive root is .15; it is not the lagged .14444… .
        near(next.moisture()[0], 0.15, 1e-12);
        near(report.boundary_water_kg[0], 0.0125, 1e-12);
        near(report.boundary_water_kg[1], 0.0125, 1e-12);
        near(report.water_balance_residual_kg, 0.0, 1e-12);
        let nonlinear = report.nonlinear.unwrap();
        assert!(nonlinear.iterations > 1);
        assert!(nonlinear.relative_residual < 1e-12);
        assert!(report.linear.converged_euclidean());
        let original = lagged.advance(10.0, &faces(0.2), config(), cx).unwrap();
        assert!(original.0.moisture()[0] < next.moisture()[0] - 0.005);
        assert!(original.1.nonlinear.is_none());
        assert_eq!(lagged.definition_identity(), implicit.definition_identity());
        assert_eq!(implicit.moisture(), &[0.1]);
        assert_eq!(implicit.time_s(), 0.0);
        let replay = implicit.advance(10.0, &faces(0.2), config(), cx).unwrap();
        assert_eq!(replay.0.moisture(), next.moisture());
    });
}

#[test]
fn g1_fine_stiff_grid_uses_one_linear_sweep_and_preserves_sealed_inventory() {
    with_cx(|cx| {
        let n = 256;
        let length = 1.0;
        let dt = 0.2;
        let diffusion = 0.01;
        let wave = |i: usize| (std::f64::consts::PI * (i as f64 + 0.5) / n as f64).cos();
        let initial = (0..n).map(|i| 0.15 + 0.05 * wave(i)).collect();
        let state = slab(card([diffusion; 2], 0.3), length, initial, cx)
            .with_implicit_diffusivity(MoistureNonlinearConfig {
                maximum_iterations: 1,
                ..nonlinear()
            })
            .unwrap();
        let mut controls = config();
        controls.linear.max_iterations = 1;
        controls.linear.max_corrections = 0;
        let (next, report) = state
            .advance(
                dt,
                &[MoistureFace::Sealed, MoistureFace::Sealed],
                controls,
                cx,
            )
            .unwrap();
        // Independent Neumann cosine eigenmode of the cell-centered Laplacian.
        let dx = length / n as f64;
        let decay = 1.0
            / (1.0
                + 4.0 * dt * diffusion / (dx * dx)
                    * (std::f64::consts::PI / (2.0 * n as f64)).sin().powi(2));
        for (i, value) in next.moisture().iter().enumerate() {
            near(*value, 0.15 + 0.05 * decay * wave(i), 2e-13);
        }
        near(next.water_mass_kg(), state.water_mass_kg(), 1e-12);
        assert_eq!(report.nonlinear.unwrap().tridiagonal_solves, 1);
        assert_eq!(report.linear.iters, 1);
    });
}

#[test]
fn g3_nonlinear_diffusion_closes_with_the_same_finite_vapor_inventory() {
    with_cx(|cx| {
        let initial_rh = 0.8;
        let slope = 0.25;
        let dt = 10.0;
        let air = atmosphere::atmosphere(initial_rh, 1.0);
        let state = slab(card([0.01, 0.04], 0.3), 1.0, vec![0.1], cx)
            .with_implicit_diffusivity(nonlinear())
            .unwrap();
        let boundaries = [
            Some(atmosphere::sorption(slope)),
            Some(atmosphere::sorption(slope)),
        ];
        let candidate = air
            .advance(
                &state,
                dt,
                &boundaries,
                config(),
                atmosphere::controls(),
                cx,
            )
            .unwrap();
        let report = candidate.step();
        let vapor_capacity = report.vapor_mass_kg[0] / initial_rh;
        let dry = state.dry_mass_kg();
        // Eliminate RH with C*RH + dry*M = C*RH0 + dry*M0,
        // then solve the independent scalar quadratic of the joint BE law.
        let a = slope * initial_rh + slope * dry * 0.1 / vapor_capacity;
        let b = 1.0 + slope * dry / vapor_capacity;
        let qa = 4.0 * dt * 0.1 * b;
        let qb = 1.0 - 4.0 * dt * (0.1 * a - 0.01 * b);
        let qc = -(0.1 + 4.0 * dt * 0.01 * a);
        let expected = -2.0 * qc / (qb + (qb * qb - 4.0 * qa * qc).sqrt());
        near(candidate.moisture().moisture()[0], expected, 2e-11);
        near(
            candidate.atmosphere().relative_humidity().unwrap(),
            initial_rh - dry * (expected - 0.1) / vapor_capacity,
            2e-10,
        );
        near(report.water_balance_residual_kg, 0.0, 1e-13);
        assert!(report.transport.nonlinear.is_some());
        near(
            candidate.atmosphere().gas().mass_kg() + candidate.moisture().wet_mass_kg(),
            air.gas().mass_kg() + state.wet_mass_kg(),
            2e-13,
        );
    });
}

#[test]
fn nonlinear_work_refusal_does_not_publish_partial_physical_state() {
    with_cx(|cx| {
        let state = slab(card([0.01, 0.04], 0.3), 1.0, vec![0.1], cx)
            .with_implicit_diffusivity(MoistureNonlinearConfig {
                maximum_iterations: 1,
                ..nonlinear()
            })
            .unwrap();
        assert!(matches!(
            state.advance(10.0, &faces(0.2), config(), cx),
            Err(MoistureError::NonlinearNotConverged(_))
        ));
        assert_eq!(state.moisture(), &[0.1]);
        assert_eq!(state.time_s(), 0.0);
        assert!(
            state
                .clone()
                .with_implicit_diffusivity(MoistureNonlinearConfig {
                    relaxation: 0.0,
                    ..nonlinear()
                })
                .is_err()
        );
    });
}
