//! Depleted pump photons, finite native inventories and accepted thermal history.

use super::*;
use fs_render::volumes::heating::feedback::IncidentBeam;
use fs_render::volumes::heating::luminescence::saturation::beam::{
    BeamEmitterCrossSections, PumpColumnLimits,
};

fn owner(
    medium: &AbsorbingGrid<'_>,
    populations: &[[f64; 3]],
    lifetimes_s: [f64; 2],
    crossing_s: f64,
) -> (
    DepositionSimulation<'static>,
    SaturableLuminescentHeating<2>,
) {
    let thermal = thermal(
        (0..medium.cell_count())
            .map(|i| cell(i as u64 + 1, 1.0, 1.0, 0.0))
            .collect(),
        &SENSIBLE,
        vec![],
    );
    let pump =
        PrescribedPhotonPump::try_new(ContentHash([160; 32]), [400.0, 450.0], [0.0; 2], [0.0; 2])
            .unwrap();
    let bindings: Vec<_> = populations
        .iter()
        .copied()
        .enumerate()
        .map(|(grid_cell, initial_populations)| SaturableEmitterCell {
            grid_cell,
            thermal_cell: CellId(grid_cell as u64 + 1),
            initial_populations,
            pump,
        })
        .collect();
    let triplet = SingletTripletSource::try_new(
        source(600.0, lifetimes_s[1], [0.0, 1.0], 161),
        crossing_s,
        ArrheniusQuenching::try_new(0.0, 0.0).unwrap(),
        SingletTripletLimits::default(),
    )
    .unwrap();
    let heating = SaturableLuminescentHeating::try_new(
        medium,
        source(500.0, lifetimes_s[0], [1.0, 0.0], 162),
        triplet,
        &bindings,
        &thermal,
        medium.cell_count(),
        1e-7,
        None,
        || Ok(()),
    )
    .unwrap();
    (thermal, heating)
}

#[allow(clippy::too_many_arguments)]
fn connect(
    cx: &Cx<'_>,
    medium: &AbsorbingGrid<'_>,
    thermal: &DepositionSimulation<'_>,
    heating: SaturableLuminescentHeating<2>,
    origin: Point3,
    direction: Vec3,
    powers_w: [f64; 2],
    cross_sections_m2: &[[[f64; 2]; 2]],
) -> SaturableLuminescentHeating<2> {
    let sections: Vec<_> = cross_sections_m2
        .iter()
        .copied()
        .enumerate()
        .map(|(cell, values)| {
            BeamEmitterCrossSections::authored(cell, ContentHash([163; 32]), values).unwrap()
        })
        .collect();
    heating
        .with_depleting_beam(
            cx,
            medium,
            ContentHash([164; 32]),
            IncidentBeam {
                origin,
                unit_direction: direction,
                max_distance_m: f64::INFINITY,
                wavelengths_nm: [400.0, 450.0],
                incident_power_w: powers_w,
            },
            &sections,
            thermal,
            medium.cell_count(),
            PumpColumnLimits {
                relative_tolerance: 2e-5,
                maximum_cell_evaluations: 1_000_000,
                ..PumpColumnLimits::default()
            },
            || Ok(()),
        )
        .unwrap()
}

fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance * actual.abs().max(expected.abs()).max(1e-30),
        "{actual:.17e} != {expected:.17e} within {tolerance:e}",
    );
}

fn beam_balance(heating: &SaturableLuminescentHeating<2>) {
    let pump = heating.depleting_pump_energy().unwrap();
    for lane in 0..2 {
        near(
            pump.incident_j[lane],
            pump.transmitted_j[lane] + pump.absorbed_j[lane] + pump.passive_absorbed_j[lane],
            3e-10,
        );
        assert!(pump.residual_j[lane].abs() < 1e-7);
    }
    relative(
        pump.absorbed_j.iter().sum(),
        heating.absorbed_pump_energy_j().iter().sum(),
    );
    relative(
        pump.passive_absorbed_j.iter().sum(),
        heating.energy().passive_pump_absorption_j,
    );
    let energy = heating.energy();
    relative(
        energy.initial_excitation_j
            + pump.absorbed_j.iter().sum::<f64>()
            + pump.passive_absorbed_j.iter().sum::<f64>(),
        energy.stored_excitation_j
            + energy.escaped_j.iter().sum::<f64>()
            + energy.unresolved_radiation_j.iter().sum::<f64>()
            + heating.thermal_input_j(),
    );
    assert!(energy.residual_j.abs() < 1e-7);
    assert!(energy.absolute_roundoff_j < 1e-7);
}

#[test]
fn g1_depleting_beam_bleaches_dense_native_cells_in_actual_ray_order() {
    with_cx(|cx| {
        // The first inventory deliberately has zero host density: molecular
        // absorption depends on its actual emitter count, not host mass.
        let density = [0.0, 1.0];
        let medium = grid(cx, &density, transparent_host());
        let (thermal, heating) = owner(&medium, &[[COUNT, 0.0, 0.0]; 2], [1e12; 2], 0.0);
        let cross_sections = [[[1.0 / COUNT, 0.0], [0.0; 2]]; 2];
        let powers = [5.0 * COUNT * HC / 400.0, 0.0];
        let forward = connect(
            cx,
            &medium,
            &thermal,
            heating.clone(),
            Point3::new(-1.0, 0.5, 0.5),
            Vec3::new(1.0, 0.0, 0.0),
            powers,
            &cross_sections,
        );
        let reverse = connect(
            cx,
            &medium,
            &thermal,
            heating,
            Point3::new(3.0, 0.5, 0.5),
            Vec3::new(-1.0, 0.0, 0.0),
            powers,
            &cross_sections,
        );
        let (warm, a, audit) = thermal
            .advance_observed(&forward, 1.0, config(1.0), || Ok(()))
            .unwrap();
        let (_, b, _) = thermal
            .advance_observed(&reverse, 1.0, config(1.0), || Ok(()))
            .unwrap();

        // Independent one-shot fluence solution: for optical thickness d,
        // F_out = ln[1 + exp(-d) (exp(F_in)-1)], with F measured in 1/sigma.
        // A 1e12 s lifetime makes decay negligible at this test's tolerance.
        let middle = (1.0 + (-1.0_f64).exp() * 5.0_f64.exp_m1()).ln();
        let outgoing = (1.0 + (-2.0_f64).exp() * 5.0_f64.exp_m1()).ln();
        let expected_ground = [COUNT * (middle - 4.0), COUNT * (1.0 + outgoing - middle)];
        for (index, expected) in expected_ground.into_iter().enumerate() {
            let populations = a.populations(index).unwrap();
            near(populations[0], expected, 2e-4);
            near(populations[1], COUNT - expected, 2e-4);
            relative(populations.iter().sum(), COUNT);
            assert_eq!(populations[2], 0.0);
            near(b.populations(1 - index).unwrap()[0], populations[0], 2e-4);
        }
        assert!(a.populations(0).unwrap()[0] < a.populations(1).unwrap()[0]);
        near(
            a.depleting_pump_energy().unwrap().transmitted_j[0],
            outgoing * COUNT * HC / 400.0,
            6e-5,
        );
        assert!(
            a.depleting_pump_energy().unwrap().transmitted_j[0]
                > 3.0 * powers[0] * (-2.0_f64).exp()
        );
        assert!(a.depleting_pump_energy().unwrap().accepted_kinetic_substeps > 2);
        relative(warm.energy_j(), a.thermal_input_j());
        relative(audit.external_heat_j, a.thermal_input_j());
        beam_balance(&a);
        beam_balance(&b);
    });
}

#[test]
fn g1_depleting_beam_competing_targets_and_passive_host_heat_keep_distinct_accounts() {
    with_cx(|cx| {
        let density = [1.0; 2];
        let absorption = MassAbsorptionSpectrum::source_table(
            &[(380.0, 0.5), (500.0, 0.5), (650.0, 0.0), (900.0, 0.0)],
            ContentHash([165; 32]),
        )
        .unwrap();
        let medium = grid(cx, &density, absorption);
        let (thermal, heating) = owner(&medium, &[[0.0; 3], [COUNT, 0.0, 0.0]], [1e12; 2], 0.0);
        let powers = [2.0 * COUNT * HC / 400.0, 3.0 * COUNT * HC / 450.0];
        let heating = connect(
            cx,
            &medium,
            &thermal,
            heating,
            Point3::new(-1.0, 0.5, 0.5),
            Vec3::new(1.0, 0.0, 0.0),
            powers,
            &[
                [[0.0; 2]; 2],
                [[0.6 / COUNT, 0.4 / COUNT], [0.1 / COUNT, 0.9 / COUNT]],
            ],
        );
        let (warm, accepted, audit) = thermal
            .advance_observed(&heating, 1.0, config(1.0), || Ok(()))
            .unwrap();
        let pump = accepted.depleting_pump_energy().unwrap();
        // The upstream zero-population host has exactly tau=.5. Its heat is
        // independent of the nonlinear downstream inventory.
        let upstream_heat = powers.iter().sum::<f64>() * -(-0.5_f64).exp_m1();
        relative(
            4.0 * (warm.temperature_k(CellId(1)).unwrap() - 300.0),
            upstream_heat,
        );
        assert_eq!(accepted.populations(0), Some([0.0; 3]));
        assert!(accepted.energy().passive_pump_absorption_j > upstream_heat);

        // Equal summed cross sections and passive depths preserve a 2:3
        // photon-lane ratio. The target hazards are (2*.6+3*.1):(2*.4+3*.9)
        // = 3:7, even as the shared ground inventory bleaches.
        let populations = accepted.populations(1).unwrap();
        relative(populations[1] / populations[2], 3.0 / 7.0);
        relative(populations.iter().sum(), COUNT);
        relative(
            (pump.absorbed_j[0] / (HC / 400.0)) / (pump.absorbed_j[1] / (HC / 450.0)),
            2.0 / 3.0,
        );
        relative(
            (pump.passive_absorbed_j[0] / (HC / 400.0))
                / (pump.passive_absorbed_j[1] / (HC / 450.0)),
            2.0 / 3.0,
        );
        relative(
            (pump.transmitted_j[0] / (HC / 400.0)) / (pump.transmitted_j[1] / (HC / 450.0)),
            2.0 / 3.0,
        );
        relative(
            accepted.energy().stored_excitation_j,
            populations[1] * HC / 500.0 + populations[2] * HC / 600.0,
        );
        relative(warm.energy_j(), accepted.thermal_input_j());
        relative(audit.external_heat_j, accepted.thermal_input_j());
        beam_balance(&accepted);
    });
}

#[test]
fn g3_depleting_beam_pulse_dark_afterglow_and_retries_share_transactional_accounts() {
    with_cx(|cx| {
        let density = [1.0];
        let medium = grid(cx, &density, transparent_host());
        let (thermal, heating) = owner(&medium, &[[COUNT, 0.0, 0.0]], [1.0, 5.0], 3.0);
        let heating = connect(
            cx,
            &medium,
            &thermal,
            heating,
            Point3::new(-1.0, 0.5, 0.5),
            Vec3::new(1.0, 0.0, 0.0),
            [5.0 * COUNT * HC / 400.0, 0.0],
            &[[[1.0 / COUNT, 0.0], [0.0; 2]]],
        );
        let original_beam = heating.depleting_pump_energy().unwrap().clone();
        let mut candidate = heating.clone();
        candidate
            .step_sources(&thermal, 0.0, 0.2, 0.2, &mut || Ok(()))
            .unwrap();
        let retry = candidate
            .step_sources(&thermal, 0.0, 0.1, 0.1, &mut || Ok(()))
            .unwrap();
        let fresh = heating
            .clone()
            .step_sources(&thermal, 0.0, 0.1, 0.1, &mut || Ok(()))
            .unwrap();
        assert_eq!(retry, fresh);
        assert_eq!(candidate.depleting_pump_energy().unwrap(), &original_beam);
        assert_eq!(candidate.energy(), heating.energy());
        assert_eq!(candidate.populations(0), heating.populations(0));
        let mut checkpoints = 0;
        assert!(matches!(
            thermal.advance_observed(&heating, 0.5, config(0.5), || {
                checkpoints += 1;
                if checkpoints > 20 {
                    Err(DepositionError::Cancelled)
                } else {
                    Ok(())
                }
            }),
            Err(OpticalHeatingError::Thermal(DepositionError::Cancelled))
        ));
        assert!(checkpoints > 20);
        assert_eq!(heating.depleting_pump_energy().unwrap(), &original_beam);
        assert_eq!(heating.time_s(), 0.0);
        assert_eq!(thermal.energy_j(), 0.0);

        let heating = heating.with_pump_window(0.25, 1.0).unwrap();
        let (_, whole, audit) = thermal
            .advance_observed(&heating, 3.0, config(3.0), || Ok(()))
            .unwrap();
        let (middle, charged, _) = thermal
            .advance_observed(&heating, 1.0, config(3.0), || Ok(()))
            .unwrap();
        let (_, split, _) = middle
            .advance_observed(&charged, 3.0, config(3.0), || Ok(()))
            .unwrap();
        assert_eq!(audit.substeps, 3);
        assert_eq!(whole.populations(0), split.populations(0));
        assert_eq!(whole.energy(), split.energy());
        assert_eq!(whole.depleting_pump_energy(), split.depleting_pump_energy());
        let pump = whole.depleting_pump_energy().unwrap();
        let charged_pump = charged.depleting_pump_energy().unwrap();
        relative(pump.incident_j[0], 0.75 * 5.0 * COUNT * HC / 400.0);
        assert_eq!(pump.incident_j, charged_pump.incident_j);
        assert_eq!(pump.transmitted_j, charged_pump.transmitted_j);
        assert_eq!(pump.absorbed_j, charged_pump.absorbed_j);
        assert_eq!(
            whole.absorbed_pump_energy_j(),
            charged.absorbed_pump_energy_j()
        );
        assert!(whole.energy().stored_excitation_j < charged.energy().stored_excitation_j);
        assert!(
            whole.energy().escaped_j.iter().sum::<f64>()
                > charged.energy().escaped_j.iter().sum::<f64>()
        );
        let view = LuminescenceGrid::try_new(
            cx,
            &medium,
            whole.emission_source_identity(),
            [650.0, 800.0],
            whole.cell_emission_w(),
            1,
        )
        .unwrap();
        let seen = view
            .transport(
                cx,
                Point3::new(-1.0, 0.5, 0.5),
                Vec3::new(1.0, 0.0, 0.0),
                f64::INFINITY,
            )
            .unwrap();
        for lane in 0..2 {
            relative(
                seen.emitted_radiance_w_m2_sr[lane],
                whole.cell_emission_w()[0].1[lane] / (4.0 * core::f64::consts::PI),
            );
            assert!(seen.emitted_radiance_w_m2_sr[lane] > 0.0);
        }
        beam_balance(&whole);
    });
}

#[test]
fn g3_depleting_beam_empty_columns_transmit_without_spurious_pump() {
    with_cx(|cx| {
        let density = [0.0];
        let medium = grid(cx, &density, transparent_host());
        for misses in [false, true] {
            let initial = if misses {
                vec![[0.0, COUNT, 0.0]]
            } else {
                vec![]
            };
            let (thermal, heating) = owner(&medium, &initial, [1.0, 5.0], 0.0);
            let heating = connect(
                cx,
                &medium,
                &thermal,
                heating,
                Point3::new(-1.0, if misses { 2.0 } else { 0.5 }, 0.5),
                Vec3::new(1.0, 0.0, 0.0),
                [2.0, 3.0],
                &[],
            );
            let (warm, accepted, _) = thermal
                .advance_observed(&heating, 1.0, config(1.0), || Ok(()))
                .unwrap();
            let pump = accepted.depleting_pump_energy().unwrap();
            relative(pump.incident_j[0], 2.0);
            relative(pump.incident_j[1], 3.0);
            assert_eq!(pump.incident_j, pump.transmitted_j);
            assert_eq!(pump.absorbed_j, [0.0; 2]);
            assert_eq!(pump.passive_absorbed_j, [0.0; 2]);
            assert_eq!(accepted.absorbed_pump_energy_j(), [0.0; 2]);
            assert_eq!(accepted.energy().prompt_heat_j, 0.0);
            if misses {
                relative(
                    accepted.populations(0).unwrap()[1],
                    COUNT * (-1.0_f64).exp(),
                );
                assert!(accepted.energy().escaped_j[0] > 0.0);
                assert!(warm.energy_j() > 0.0);
            } else {
                assert!(accepted.cell_emission_w().is_empty());
                assert_eq!(accepted.energy().stored_excitation_j, 0.0);
                assert_eq!(warm.energy_j(), 0.0);
            }
            beam_balance(&accepted);
        }
    });
}
