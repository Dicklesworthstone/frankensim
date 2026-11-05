use super::*;
use crate::resistive_film::cavitation::CompressibleCavitatingFilm;

fn model(closed: bool, bulk: f64, alpha: f64) -> CompressibleCavitatingFilm {
    CompressibleCavitatingFilm::new(
        CavitatingFilm::new(
            geometry(closed).with_barus_viscosity(alpha).unwrap(),
            liquid(),
            config(),
        )
        .unwrap(),
        bulk,
    )
    .unwrap()
}

fn approx(actual: f64, expected: f64, relative: f64) {
    assert!(
        (actual - expected).abs() <= relative * expected.abs().max(1e-12),
        "actual={actual:.16e}, expected={expected:.16e}"
    );
}

// Independent dimensional Simpson integration of rho(p)/mu(p). The native
// pressure-flow code never uses this quadrature or the standard-library exp.
fn pressure_mass_flow(p: f64, q: f64, bulk: f64, alpha: f64) -> f64 {
    let count = 1024;
    let dx = (p - q) / count as f64;
    let medium = liquid();
    let density_fluidity = |pressure: f64| {
        let gauge = pressure - medium.ambient_pressure_pa;
        let rho = medium.density_kg_m3 * (gauge / bulk).exp();
        rho / (alpha * gauge).exp()
    };
    let mut integral = density_fluidity(q) + density_fluidity(p);
    for i in 1..count {
        integral += if i % 2 == 0 { 2.0 } else { 4.0 } * density_fluidity(q + i as f64 * dx);
    }
    1e-10 * dx * integral / 3.0
}

#[test]
fn g1_sealed_compression_release_and_reformation_preserve_original_mass() {
    let bulk = 2e9;
    let film = model(true, bulk, 0.0);
    let original = film.initial_state(&[0.0]).unwrap();
    let saved = original.clone();
    let compressed = film
        .advance(&original, &[0.0001], &[0.0], 0.01, || false)
        .unwrap();
    let pressure = liquid().ambient_pressure_pa + bulk * (1.0_f64 / 0.9).ln();
    approx(compressed.state.pressure_pa()[0], pressure, 2e-12);
    approx(
        compressed.generalized_force[0],
        0.001 * (pressure - 100_000.0),
        2e-12,
    );
    assert_eq!(compressed.state.liquid_fraction(), &[1.0]);
    assert_eq!(compressed.channel_mass_flow_kg_s, vec![0.0]);
    approx(compressed.report.liquid_mass_after_kg, 0.001, 2e-12);
    assert!(compressed.report.mass_residual_kg.abs() < 1e-16);

    let released = film
        .advance(&compressed.state, &[-0.001], &[0.0], 0.03, || false)
        .unwrap();
    let rho_c = 1000.0 * ((2000.0 - 100_000.0) / bulk).exp();
    approx(
        released.state.liquid_fraction()[0],
        0.001 / (2e-6 * rho_c),
        2e-12,
    );
    assert_eq!(released.state.pressure_pa(), &[2000.0]);
    approx(released.report.liquid_mass_after_kg, 0.001, 2e-12);
    let reformed = film
        .advance(&released.state, &[0.0], &[0.0], 0.1, || false)
        .unwrap();
    approx(reformed.state.pressure_pa()[0], 100_000.0, 2e-9);
    assert_eq!(reformed.state.liquid_fraction(), &[1.0]);
    assert_eq!(original, saved);
}

#[test]
fn g1_connected_compression_refills_a_cavity_with_actual_barus_mass_transport() {
    let bulk = 1e8;
    let alpha = 2e-8;
    let graph = ResistiveFilm::new(
        vec![
            FilmCell {
                area_m2: 0.001,
                gap: gap(0.001, 1.0),
            },
            FilmCell {
                area_m2: 0.001,
                gap: gap(0.001, 0.0),
            },
        ],
        vec![FilmChannel {
            from: 0,
            to: Some(1),
            width_m: 0.0012,
            length_m: 0.1,
            gap: gap(0.001, 0.0),
        }],
        1,
        0.01,
        FilmLimits {
            minimum_cell_gap_m: 1e-8,
            maximum_gap_m: 0.1,
            maximum_pressure_pa: 1e9,
        },
    )
    .unwrap()
    .with_barus_viscosity(alpha)
    .unwrap();
    let film = CompressibleCavitatingFilm::new(
        CavitatingFilm::new(graph, liquid(), config()).unwrap(),
        bulk,
    )
    .unwrap();
    let old = film.initialize(&[0.0], &[1.0, 0.3]).unwrap();
    let dt = 0.001;
    let next = film.advance(&old, &[0.0002], &[0.0], dt, || false).unwrap();
    let mut low = 2000.0;
    let mut high = 100_000_000.0;
    for _ in 0..100 {
        let p = 0.5 * (low + high);
        let stored = 1000.0 * ((p - 100_000.0) / bulk).exp() * 8e-7;
        if stored + dt * pressure_mass_flow(p, 2000.0, bulk, alpha) > 0.001 {
            high = p;
        } else {
            low = p;
        }
    }
    let expected_pressure = 0.5 * (low + high);
    let transferred = dt * pressure_mass_flow(expected_pressure, 2000.0, bulk, alpha);
    approx(next.state.pressure_pa()[0], expected_pressure, 2e-9);
    assert_eq!(next.state.pressure_pa()[1], 2000.0);
    approx(next.channel_mass_flow_kg_s[0], transferred / dt, 2e-9);
    let rho_c = 1000.0 * (-98_000.0 / bulk).exp();
    approx(
        next.state.liquid_fraction()[1],
        0.3 + transferred / (rho_c * 1e-6),
        2e-9,
    );
    approx(
        next.report.liquid_mass_after_kg,
        0.001 + rho_c * 3e-7,
        2e-11,
    );
    assert_eq!(next.report.boundary_mass_out_kg, 0.0);
    for i in 0..2 {
        let old_mass = if i == 0 { 0.001 } else { rho_c * 3e-7 };
        let volume = if i == 0 { 8e-7 } else { 1e-6 };
        let rho = 1000.0 * ((next.state.pressure_pa()[i] - 100_000.0) / bulk).exp();
        let out = if i == 0 {
            next.channel_mass_flow_kg_s[0]
        } else {
            -next.channel_mass_flow_kg_s[0]
        };
        let residual = rho * volume * next.state.liquid_fraction()[i] - old_mass + dt * out;
        assert!(residual.abs() < 2e-14, "cell {i}: {residual}");
    }
}

#[test]
fn g1_supplied_compressible_throughflow_and_viscosity_refresh_match_independent_integral() {
    let bulk = 2e7;
    // Include alpha*B=1, the logarithmic removable limit of the primitive.
    for alpha in [0.0, 5e-8, 1e-7] {
        let native = super::reservoirs::supplied_film(alpha);
        let film = CompressibleCavitatingFilm::new(native, bulk).unwrap();
        for (left, right) in [(2e6, 2000.0), (2000.0, 2e6)] {
            let supplied = film
                .with_reservoirs(&[
                    Some(CavitationReservoir {
                        pressure_pa: left,
                        liquid_fraction: 1.0,
                    }),
                    Some(CavitationReservoir {
                        pressure_pa: right,
                        liquid_fraction: 1.0,
                    }),
                ])
                .unwrap();
            let mut low = left.min(right);
            let mut high = left.max(right);
            for _ in 0..100 {
                let p = 0.5 * (low + high);
                let balance = pressure_mass_flow(p, left, bulk, alpha)
                    + pressure_mass_flow(p, right, bulk, alpha);
                if balance > 0.0 {
                    high = p;
                } else {
                    low = p;
                }
            }
            let pressure = 0.5 * (low + high);
            let state = supplied
                .initialize_with_pressure(&[0.0], &[1.0], &[pressure])
                .unwrap();
            let first = supplied
                .advance(&state, &[0.0], &[0.0; 2], 0.01, || false)
                .unwrap();
            approx(first.state.pressure_pa()[0], pressure, 2e-11);
            let flow = pressure_mass_flow(pressure, left, bulk, alpha);
            approx(first.channel_mass_flow_kg_s[0], flow, 2e-11);
            approx(first.channel_mass_flow_kg_s[1], -flow, 2e-11);
            assert!(first.report.mass_residual_kg.abs() < 1e-14);
            let thicker = supplied.with_channel_viscosities(&[0.02, 0.02]).unwrap();
            let second = thicker
                .advance(&first.state, &[0.0], &[0.0; 2], 0.01, || false)
                .unwrap();
            approx(second.state.pressure_pa()[0], pressure, 2e-11);
            approx(second.channel_mass_flow_kg_s[0], flow / 2.0, 2e-11);
        }
    }
}

#[test]
fn g3_native_starvation_advection_and_incompressible_limit_are_recovered() {
    let native = two_cells();
    let film = CompressibleCavitatingFilm::new(native.clone(), 2e9).unwrap();
    let state = film.initialize(&[0.0], &[0.2, 0.1]).unwrap();
    let old = native.initialize(&[0.0], &[0.2, 0.1]).unwrap();
    let flux = [-1e-5, 1e-5, 1e-5];
    let result = film.advance(&state, &[0.0], &flux, 0.1, || false).unwrap();
    let original = native.advance(&old, &[0.0], &flux, 0.1, || false).unwrap();
    for i in 0..2 {
        approx(
            result.state.liquid_fraction()[i],
            original.state.liquid_fraction()[i],
            2e-12,
        );
    }
    for i in 0..3 {
        approx(
            result.channel_mass_flow_kg_s[i],
            1000.0 * original.channel_volume_flow_m3_s[i],
            2e-12,
        );
    }
    let mut previous_error = f64::INFINITY;
    for bulk in [1e8, 1e10, 1e12] {
        let film = model(false, bulk, 0.0);
        let old = film.initial_state(&[0.0]).unwrap();
        let result = film
            .advance(&old, &[0.0001], &[0.0], 0.01, || false)
            .unwrap();
        // Exact one-cell barotropic BE balance, independent of the native
        // coordinate or nonlinear solve, and well conditioned as B increases.
        let exact = 100_000.0 + bulk * (1e-7 / (9e-7 + 1e-12 * bulk)).ln_1p();
        approx(result.state.pressure_pa()[0], exact, 2e-10);
        let error = (result.state.pressure_pa()[0] - 200_000.0).abs();
        assert!(error < previous_error / 50.0);
        previous_error = error;
    }
    assert!(previous_error < 0.1);
}

#[test]
fn g0_g4_domain_budget_identity_and_late_cancellation_are_atomic() {
    let native = CavitatingFilm::new(geometry(false), liquid(), config()).unwrap();
    for bulk in [0.0, -1.0, f64::NAN, f64::INFINITY, 1e-10] {
        assert!(CompressibleCavitatingFilm::new(native.clone(), bulk).is_err());
    }
    let film = model(false, 2e9, 1e-7);
    let state = film.initial_state(&[0.0]).unwrap();
    let saved = state.clone();
    let native_state = native.initial_state(&[0.0]).unwrap();
    assert!(
        film.advance(&native_state, &[0.0], &[0.0], 0.01, || false)
            .is_err()
    );
    assert!(
        native
            .advance(&state, &[0.0], &[0.0], 0.01, || false)
            .is_err()
    );
    assert!(
        film.initialize_with_pressure(&[0.0], &[0.5], &[100_000.0])
            .is_err()
    );
    assert!(film.initialize_with_pressure(&[0.0], &[1.0], &[]).is_err());
    assert!(film.with_channel_viscosities(&[]).is_err());
    assert!(
        film.with_reservoirs(&[Some(CavitationReservoir {
            pressure_pa: 1000.0,
            liquid_fraction: 1.0,
        })])
        .is_err()
    );
    let mut polls = 0;
    film.advance(&state, &[0.0001], &[0.0], 0.01, || {
        polls += 1;
        false
    })
    .unwrap();
    let mut actual = 0;
    assert!(
        film.advance(&state, &[0.0001], &[0.0], 0.01, || {
            actual += 1;
            actual == polls
        })
        .is_err()
    );
    assert!(
        film.advance(&state, &[0.0001], &[0.0], 0.01, || true)
            .is_err()
    );
    assert_eq!(state, saved);

    let sealed = model(true, 2e9, 0.0);
    let old = sealed.initial_state(&[0.0]).unwrap();
    assert!(
        sealed
            .advance(&old, &[0.0009], &[0.0], 0.01, || false)
            .is_err()
    );
    assert!(
        sealed
            .advance(&old, &[0.0], &[1e-6], 0.01, || false)
            .is_err()
    );
    let limited = CompressibleCavitatingFilm::new(
        CavitatingFilm::new(
            geometry(false).with_barus_viscosity(1e-7).unwrap(),
            liquid(),
            CavitationSolveConfig {
                maximum_iterations: 1,
                ..config()
            },
        )
        .unwrap(),
        2e9,
    )
    .unwrap();
    let old = limited.initial_state(&[0.0]).unwrap();
    assert!(
        limited
            .advance(&old, &[0.0002], &[0.0], 0.01, || false)
            .is_err()
    );
}
