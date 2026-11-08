use super::*;
use fs_flux::ns::newton::{FluxSdirk2Error, FluxSdirk2Stage};
use std::cell::Cell;

#[test]
fn sdirk2_viscous_mode_is_second_order_and_damps_unresolved_stiff_decay() {
    let mesh = mesh(1);
    let system = FluxSystem::new(&mesh);
    let params = FluxParams {
        sigma: 15.0,
        tol: 1e-12,
        ..FluxParams::default()
    };
    let moment = 2 * mesh
        .edges
        .iter()
        .position(|edge| edge.tris.1 != usize::MAX)
        .unwrap()
        + 1;
    let nu = 0.03;
    // Independent integrals of the sole free divergence-free BDM1 mode on
    // two unit-square triangles: M=24, viscous K/nu=sigma*384-288.
    // Turning off convection isolates physical linear viscous decay.
    let k = 6.0 * core::f64::consts::SQRT_2;
    let mass = k * k * 2.0 / 6.0;
    let penalty = k * k * (4.0 / 3.0 + 4.0);
    let strain = 4.0 * k * k;
    let decay = nu * (params.sigma * penalty - strain) / mass;
    let gamma = 1.0 - 1.0 / 2.0_f64.sqrt();
    let amplification = |z: f64| (1.0 - (1.0 - 2.0 * gamma) * z) / (1.0 + gamma * z).powi(2);
    let law = |_, _, _, _| FluxViscosityTangent {
        nu,
        derivative: 0.0,
    };
    let controls = FluxNewtonControls {
        residual_tolerance: 1e-11,
        ..FluxNewtonControls::default()
    };
    let mut errors = Vec::new();
    for count in [4, 8, 16] {
        let timestep = 1.0 / (decay * f64::from(count));
        let mut previous = vec![0.0; system.n_u];
        previous[moment] = 0.1;
        for _ in 0..count {
            let step = system
                .sdirk2_rheology_newton_step(
                    params,
                    &law,
                    &|_, _| [0.0; 2],
                    &|_, _| [0.0; 2],
                    &previous,
                    timestep,
                    false,
                    controls,
                )
                .unwrap();
            let expected_first = previous[moment] / (1.0 + gamma * decay * timestep);
            let expected = previous[moment] * amplification(decay * timestep);
            assert!((step.first_stage.solution.x[moment] - expected_first).abs() < 5e-11);
            assert!((step.endpoint.solution.x[moment] - expected).abs() < 5e-11);
            let (_, divergence) = system.velocity_error(&step.endpoint.solution.x, &|_| [0.0; 2]);
            assert!(divergence < 1e-10);
            previous.copy_from_slice(&step.endpoint.solution.x[..system.n_u]);
        }
        errors.push((previous[moment] - 0.1 * (-1.0_f64).exp()).abs());
    }
    for pair in errors.windows(2) {
        let ratio = pair[0] / pair[1];
        assert!(
            (3.8..4.3).contains(&ratio),
            "second-order decay errors: {errors:?}"
        );
    }

    let mut previous = vec![0.0; system.n_u];
    previous[moment] = 0.1;
    let original = previous.clone();
    let step = system
        .sdirk2_rheology_newton_step(
            params,
            &law,
            &|_, _| [0.0; 2],
            &|_, _| [0.0; 2],
            &previous,
            1e6 / decay,
            false,
            controls,
        )
        .unwrap();
    assert!((step.endpoint.solution.x[moment] - 0.1 * amplification(1e6)).abs() < 5e-11);
    assert!(step.endpoint.solution.x[moment].abs() < 1e-5 * previous[moment].abs());
    assert_eq!(previous, original);
}

#[test]
fn sdirk2_uses_each_stages_rheology_force_and_boundary_time() {
    let mesh = mesh(3);
    let system = FluxSystem::new(&mesh);
    let params = FluxParams {
        tol: 1e-12,
        sigma: 15.0,
        ..FluxParams::default()
    };
    let time = 0.2;
    let timestep = 0.4;
    let slope = 0.7;
    let amplitude = |t: f64| 1.0 + slope * t;
    let previous = couette_moments(&mesh, amplitude(time));
    let original = previous.clone();
    let reads = [Cell::new(0_usize), Cell::new(0_usize)];
    let index = |stage| match stage {
        FluxSdirk2Stage::First => 0,
        FluxSdirk2Stage::Endpoint => 1,
    };
    let gradient = |stage| [1.5, -0.7][index(stage)];
    let law = |stage, _, _, rate: f64| {
        let i = index(stage);
        reads[i].set(reads[i].get() + 1);
        let base = [0.2, 0.6][i];
        FluxViscosityTangent {
            nu: base * (1.0 + rate * rate),
            derivative: 2.0 * base * rate,
        }
    };
    let force = |stage, p: [f64; 2]| [slope * p[1] + gradient(stage), 0.0];
    let boundary = |stage: FluxSdirk2Stage, p: [f64; 2]| {
        [
            amplitude(time + stage.time_fraction() * timestep) * p[1],
            0.0,
        ]
    };
    let result = system
        .sdirk2_rheology_newton_step(
            params,
            &law,
            &force,
            &boundary,
            &previous,
            timestep,
            true,
            FluxNewtonControls {
                residual_tolerance: 1e-10,
                ..FluxNewtonControls::default()
            },
        )
        .unwrap();
    for (stage, step) in [
        (FluxSdirk2Stage::First, &result.first_stage),
        (FluxSdirk2Stage::Endpoint, &result.endpoint),
    ] {
        let expected = |p| boundary(stage, p);
        let (error, divergence) = system.velocity_error(&step.solution.x, &expected);
        assert!(error < 3e-9 && divergence < 3e-9);
        let shear = amplitude(time + stage.time_fraction() * timestep);
        assert!(
            system
                .shear_rates(&step.solution.x)
                .iter()
                .all(|rate| (rate - shear).abs() < 3e-8)
        );
        for (pressure, point) in system
            .pressure(&step.solution.x)
            .iter()
            .zip(&mesh.centroids)
        {
            assert!((pressure - gradient(stage) * (point[0] - 0.5)).abs() < 2e-7);
        }
        assert!(reads[index(stage)].get() > 0);
    }
    assert_eq!(previous, original);
    let invalid = system.sdirk2_rheology_newton_step(
        params,
        &|_, _, _, _| panic!("invalid timestep queried viscosity"),
        &|_, _| panic!("invalid timestep queried force"),
        &|_, _| panic!("invalid timestep queried boundary"),
        &previous,
        0.0,
        true,
        FluxNewtonControls::default(),
    );
    assert!(matches!(
        invalid,
        Err(FluxSdirk2Error {
            stage: None,
            error: FluxNewtonError::Invalid(_)
        })
    ));
}
