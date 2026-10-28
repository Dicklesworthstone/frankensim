//! G1/G3 finite-body vacuum exchange with independent scalar energy balances.
use super::*;

fn boundary(fixture: &Fixture, laws: [ThermalBc; 2]) -> ThermalBoundary {
    let [left, right] = laws;
    ThermalBoundaryBuilder::new(&fixture.mesh)
        .region("left", |face| face.vertices.iter().all(|v| *v < 4), left)
        .unwrap()
        .region("right", |face| face.vertices.iter().all(|v| *v >= 4), right)
        .unwrap()
        .finish()
        .unwrap()
}

fn problem<'a>(
    fixture: &'a Fixture,
    boundary: &'a ThermalBoundary,
    source: &'a ScalarField,
) -> ConductionProblem<'a> {
    ConductionProblem {
        boundary,
        ..fixture.problem(source)
    }
}

fn sensible_history(temperature: [f64; 2]) -> Vec<f64> {
    temperature
        .into_iter()
        .flat_map(|t| [CP * (t - 200.0); 4])
        .collect()
}

// Eliminate each body's temperature from its own finite-capacity balance.
// The remaining scalar unknown Q is gray-body exchange, positive left-to-right.
// This never invokes the FEM assembler or the production radiosity solver.
fn sensible_reference(old: [f64; 2], htc: [f64; 2], prescribed_flux: [f64; 2]) -> ([f64; 2], f64) {
    let area = 3.0_f64.sqrt();
    let capacity = RHO * CP * 2.0_f64.sqrt() / 12.0;
    let denominator = htc.map(|h| capacity + DT * area * h);
    let free: [f64; 2] = std::array::from_fn(|i| {
        (capacity * old[i] + DT * area * (300.0 * htc[i] - prescribed_flux[i])) / denominator[i]
    });
    let temperatures = |heat: f64| {
        [
            free[0] - DT * heat / denominator[0],
            free[1] + DT * heat / denominator[1],
        ]
    };
    let equilibrium_heat = (free[0] - free[1]) / (DT / denominator[0] + DT / denominator[1]);
    let (mut lo, mut hi) = (equilibrium_heat.min(0.0), equilibrium_heat.max(0.0));
    let factor = area * STEFAN_BOLTZMANN_W_M2_K4 / (1.0 / EPS[0] + 1.0 / EPS[1] - 1.0);
    for _ in 0..90 {
        let heat = (lo + hi) / 2.0;
        let t = temperatures(heat);
        if heat - factor * (t[0].powi(4) - t[1].powi(4)) > 0.0 {
            hi = heat;
        } else {
            lo = heat;
        }
    }
    let heat = (lo + hi) / 2.0;
    (temperatures(heat), heat)
}

#[test]
fn insulated_finite_bodies_exchange_radiation_without_air_or_external_heat() {
    let fixture = Fixture::new(false);
    let boundary = boundary(&fixture, [ThermalBc::adiabatic(), ThermalBc::adiabatic()]);
    let source = ScalarField::Uniform(0.0);
    with_cx(|cx| {
        let engine = fixture.engine(cx, false);
        for initial in [[450.0, 300.0], [300.0, 450.0], [375.0, 375.0]] {
            let old = sensible_history(initial);
            let history = old.clone();
            let got = engine
                .advance_with_enclosure_radiation(
                    cx,
                    problem(&fixture, &boundary, &source),
                    None,
                    &old,
                    DT,
                    step_config(),
                    &fixture.enclosure,
                    radiation_config(),
                )
                .unwrap();
            let (expected, heat) = sensible_reference(initial, [0.0; 2], [0.0; 2]);
            for (body, &temperature) in expected.iter().enumerate() {
                for &actual in &got.conduction.temperature[4 * body..4 * body + 4] {
                    close(actual, temperature, 2e-10);
                }
            }
            close(got.radiation.radiosity.net_outward_heat_w[0], heat, 2e-9);
            close(got.radiation.radiosity.net_outward_heat_w[1], -heat, 2e-9);
            close(got.conduction.stored_energy_change_j, 0.0, 1e-8);
            close(got.physical_energy_residual_j, 0.0, 1e-8);
            assert!(got.physical_residual_norm_j <= got.physical_residual_tolerance_j);
            assert_eq!(got.prescribed_neumann_out_w, 0.0);
            assert_eq!(got.convective_out_w, 0.0);
            assert!(got.convective_robin_fluxes.is_empty());
            assert!(got.conduction.robin_fluxes.is_empty());
            assert!(
                got.combined_boundary
                    .conditions()
                    .iter()
                    .all(|law| matches!(law, ThermalBc::Neumann { .. }))
            );
            assert_eq!(old, history);
        }
    });
}

#[test]
fn vacuum_receiver_melts_at_fixed_temperature_using_the_donors_finite_enthalpy() {
    let fixture = Fixture::new(true);
    let boundary = boundary(&fixture, [ThermalBc::adiabatic(), ThermalBc::adiabatic()]);
    let source = ScalarField::Uniform(0.0);
    let old = [
        3500.0, 3500.0, 3500.0, 3500.0, 1500.0, 1500.0, 1500.0, 1500.0,
    ];
    let mass = RHO * 2.0_f64.sqrt() / 12.0;
    let factor = 3.0_f64.sqrt() * STEFAN_BOLTZMANN_W_M2_K4 / (1.0 / EPS[0] + 1.0 / EPS[1] - 1.0);
    let (mut lo, mut hi) = (350.0_f64, 400.0_f64);
    for _ in 0..90 {
        let t = (lo + hi) / 2.0;
        let residual = mass * CP * (t - 400.0) + DT * factor * (t.powi(4) - 350.0_f64.powi(4));
        if residual > 0.0 {
            hi = t;
        } else {
            lo = t;
        }
    }
    let donor_t = (lo + hi) / 2.0;
    let heat = factor * (donor_t.powi(4) - 350.0_f64.powi(4));
    with_cx(|cx| {
        let got = fixture
            .engine(cx, true)
            .advance_with_enclosure_radiation(
                cx,
                problem(&fixture, &boundary, &source),
                None,
                &old,
                DT,
                step_config(),
                &fixture.enclosure,
                radiation_config(),
            )
            .unwrap();
        for node in 0..4 {
            close(got.conduction.temperature[node], donor_t, 2e-10);
            close(
                got.conduction.specific_enthalpy_j_kg[node],
                3500.0 - DT * heat / mass,
                2e-10,
            );
            assert_eq!(got.conduction.temperature[node + 4], 350.0);
            close(
                got.conduction.specific_enthalpy_j_kg[node + 4],
                1500.0 + DT * heat / mass,
                2e-10,
            );
            close(
                got.conduction.liquid_mass_fraction[node + 4],
                (500.0 + DT * heat / mass) / 2000.0,
                2e-11,
            );
        }
        assert!(got.conduction.liquid_mass_fraction[4] > 0.25);
        close(got.conduction.stored_energy_change_j, 0.0, 1e-8);
        close(got.physical_energy_residual_j, 0.0, 1e-8);
    });
}

#[test]
fn authored_flux_survives_exchange_and_mixed_boundary_totals_count_heat_once() {
    let fixture = Fixture::new(false);
    let source = ScalarField::Uniform(0.0);
    let initial = [440.0, 300.0];
    let old = sensible_history(initial);
    let area = 3.0_f64.sqrt();
    with_cx(|cx| {
        for (htc, imposed) in [([0.0, HTC], [-25.0, 0.0]), ([0.0; 2], [-25.0, 7.0])] {
            let laws = std::array::from_fn(|i| {
                if htc[i] > 0.0 {
                    ThermalBc::robin(htc[i], 300.0).unwrap()
                } else {
                    ThermalBc::neumann(imposed[i]).unwrap()
                }
            });
            let boundary = boundary(&fixture, laws);
            let got = fixture
                .engine(cx, false)
                .advance_with_enclosure_radiation(
                    cx,
                    problem(&fixture, &boundary, &source),
                    None,
                    &old,
                    DT,
                    step_config(),
                    &fixture.enclosure,
                    radiation_config(),
                )
                .unwrap();
            let (expected, heat) = sensible_reference(initial, htc, imposed);
            for (body, &temperature) in expected.iter().enumerate() {
                for &actual in &got.conduction.temperature[4 * body..4 * body + 4] {
                    close(actual, temperature, 2e-10);
                }
            }
            let prescribed = area * (imposed[0] + imposed[1]);
            let convection =
                area * (htc[0] * (expected[0] - 300.0) + htc[1] * (expected[1] - 300.0));
            close(got.prescribed_neumann_out_w, prescribed, 1e-12);
            close(got.convective_out_w, convection, 2e-9);
            close(
                got.conduction.stored_energy_change_j,
                -DT * (prescribed + convection),
                1e-8,
            );
            close(got.physical_energy_residual_j, 0.0, 1e-8);
            close(got.radiation.radiosity.net_outward_heat_w[0], heat, 2e-9);
            if htc[1] > 0.0 {
                assert_eq!(got.convective_robin_fluxes.len(), 1);
                assert_eq!(got.convective_robin_fluxes[0].region, "right");
                // Internal heat enters the combined Neumann total but must
                // never enter the independently integrated external account.
                close(got.conduction.neumann_out_w, prescribed + heat, 2e-9);
                assert!((got.conduction.neumann_out_w - prescribed).abs() > 1.0);
            }
        }
    });
}

#[test]
fn relabeling_and_permuting_radiating_surfaces_preserves_the_finite_body_solution() {
    let fixture = Fixture::new(false);
    let source = ScalarField::Uniform(0.0);
    let original = boundary(&fixture, [ThermalBc::adiabatic(), ThermalBc::adiabatic()]);
    let renamed = ThermalBoundaryBuilder::new(&fixture.mesh)
        .region(
            "receiver",
            |f| f.vertices.iter().all(|v| *v >= 4),
            ThermalBc::adiabatic(),
        )
        .unwrap()
        .region(
            "emitter",
            |f| f.vertices.iter().all(|v| *v < 4),
            ThermalBc::adiabatic(),
        )
        .unwrap()
        .finish()
        .unwrap();
    let receiver = RadiationSurface::new(
        &fixture.mesh,
        "receiver",
        |f| f.vertices.iter().all(|v| *v >= 4),
        emissivity("receiver", EPS[1]),
    )
    .unwrap();
    let emitter = RadiationSurface::new(
        &fixture.mesh,
        "emitter",
        |f| f.vertices.iter().all(|v| *v < 4),
        emissivity("emitter", EPS[0]),
    )
    .unwrap();
    let factors = ViewFactorMatrix::infinite_parallel_plates(receiver.area_m2()).unwrap();
    let enclosure = GrayDiffuseEnclosure::new(vec![receiver, emitter], factors).unwrap();
    let old = sensible_history([450.0, 300.0]);
    with_cx(|cx| {
        let engine = fixture.engine(cx, false);
        let before = engine
            .advance_with_enclosure_radiation(
                cx,
                problem(&fixture, &original, &source),
                None,
                &old,
                DT,
                step_config(),
                &fixture.enclosure,
                radiation_config(),
            )
            .unwrap();
        let after = engine
            .advance_with_enclosure_radiation(
                cx,
                problem(&fixture, &renamed, &source),
                None,
                &old,
                DT,
                step_config(),
                &enclosure,
                radiation_config(),
            )
            .unwrap();
        for (&a, &b) in before
            .conduction
            .specific_enthalpy_j_kg
            .iter()
            .zip(&after.conduction.specific_enthalpy_j_kg)
        {
            close(a, b, 2e-10);
        }
        assert!(after.radiation.radiosity.net_outward_heat_w[0] < 0.0);
        assert!(after.radiation.radiosity.net_outward_heat_w[1] > 0.0);
        close(after.conduction.stored_energy_change_j, 0.0, 1e-8);
    });
}

#[test]
fn vacuum_source_derivative_retains_reflected_exchange_and_refuses_fake_air_ports() {
    let fixture = Fixture::new(false);
    let boundary = boundary(&fixture, [ThermalBc::adiabatic(), ThermalBc::adiabatic()]);
    let source = ScalarField::Uniform(0.0);
    let old = sensible_history([450.0, 300.0]);
    with_cx(|cx| {
        let engine = fixture.engine(cx, false);
        let accepted = engine
            .advance_with_enclosure_radiation(
                cx,
                problem(&fixture, &boundary, &source),
                None,
                &old,
                DT,
                step_config(),
                &fixture.enclosure,
                radiation_config(),
            )
            .unwrap();
        let linear = engine
            .linearize_accepted_with_enclosure_radiation(
                cx,
                problem(&fixture, &boundary, &source),
                None,
                &old,
                DT,
                step_config(),
                &fixture.enclosure,
                radiation_config(),
                accepted.conduction,
                32,
            )
            .unwrap();
        let goal = [0.0, 0.0, 0.0, 0.0, 0.25, 0.25, 0.25, 0.25];
        let seed = linear.temperature_pullback(cx, &goal).unwrap();
        let policy = LinearConfig {
            tolerance: 1e-11,
            max_iterations: 128,
            restart: 8,
        };
        let gradient = linear.pullback(cx, &seed, policy).unwrap();
        let objective = |power| {
            let source = ScalarField::Nodal(vec![power, power, power, power, 0.0, 0.0, 0.0, 0.0]);
            let got = engine
                .advance_with_enclosure_radiation(
                    cx,
                    problem(&fixture, &boundary, &source),
                    None,
                    &old,
                    DT,
                    step_config(),
                    &fixture.enclosure,
                    radiation_config(),
                )
                .unwrap();
            dot(&got.conduction.temperature, &goal)
        };
        let finite_difference = (objective(0.2) - objective(-0.2)) / 0.4;
        let derivative: f64 = gradient.source_density[..4].iter().sum();
        assert!(derivative > 0.0);
        close(derivative, finite_difference, 2e-8);
        assert!(linear.robin_response(cx, &["left"], policy).is_err());
        let history = old.clone();
        let mut exhausted = radiation_config();
        exhausted.max_iterations = 1;
        assert!(
            engine
                .advance_with_enclosure_radiation(
                    cx,
                    problem(&fixture, &boundary, &source),
                    None,
                    &old,
                    DT,
                    step_config(),
                    &fixture.enclosure,
                    exhausted,
                )
                .is_err()
        );
        with_cancelled_cx(|cancelled| {
            assert!(
                engine
                    .advance_with_enclosure_radiation(
                        cancelled,
                        problem(&fixture, &boundary, &source),
                        None,
                        &old,
                        DT,
                        step_config(),
                        &fixture.enclosure,
                        radiation_config(),
                    )
                    .is_err()
            )
        });
        assert_eq!(old, history);
    });
}
