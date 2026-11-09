//! Complete elastic walks share the same evolving G/S/T and heat accounts.

use super::*;

fn medium<'a>(
    cx: &Cx<'_>,
    density: &'a [f64],
    passive: f64,
    scattering: f64,
) -> ScatteringGrid<'a> {
    ScatteringGrid::try_new(
        cx,
        grid(cx, density, passive_spectrum(passive)),
        MassScatteringSpectrum::source_table(
            &[(380.0, scattering), (800.0, scattering)],
            ContentHash([82; 32]),
        )
        .unwrap(),
        PhaseFunction::Isotropic,
    )
    .unwrap()
}

#[allow(clippy::too_many_arguments)]
fn owner(
    cx: &Cx<'_>,
    medium: &ScatteringGrid<'_>,
    populations: &[[f64; 3]],
    thermal: &DepositionSimulation<'_>,
    sections: [[f64; 2]; 2],
    geometry: PhotonRecyclingGeometry,
    events: u32,
) -> Result<SaturableLuminescentHeating<2>, OpticalHeatingError> {
    let (singlet, triplet) = sources();
    SaturableLuminescentHeating::try_new_with_scattered_photon_recycling(
        cx,
        medium,
        singlet,
        triplet,
        &bindings(populations),
        thermal,
        populations.len(),
        1e-7,
        &molecule(populations.len(), sections),
        geometry,
        events,
        controls(),
        || Ok(()),
    )
}

#[test]
fn g1_zero_scattering_preserves_ballistic_identity_populations_and_heat() {
    with_cx(|cx| {
        let density = [1.0; 2];
        let medium = medium(cx, &density, 0.2, 0.0);
        let initial = [[0.5 * COUNT, 0.5 * COUNT, 0.0], [COUNT, 0.0, 0.0]];
        let sections = [[0.8 / COUNT, 0.0], [0.0; 2]];
        let (thermal, ballistic) = super::owner(cx, medium.absorption_grid(), &initial, sections);
        let elastic = owner(cx, &medium, &initial, &thermal, sections, geometry(), 0).unwrap();
        assert_eq!(elastic.emission_source_identity(), ballistic.emission_source_identity());
        let (warm_a, accepted_a, _) = thermal
            .advance_observed(&ballistic, 0.4, config(0.4), || Ok(()))
            .unwrap();
        let (warm_b, accepted_b, _) = thermal
            .advance_observed(&elastic, 0.4, config(0.4), || Ok(()))
            .unwrap();
        assert_eq!(warm_a.energy_j(), warm_b.energy_j());
        assert_eq!(accepted_a.energy(), accepted_b.energy());
        for cell in 0..2 {
            assert_eq!(accepted_a.populations(cell), accepted_b.populations(cell));
        }
    });
}

#[test]
fn g3_elastic_return_reexcites_original_populations_and_deposits_passive_heat() {
    with_cx(|cx| {
        let density = [1.0; 2];
        let medium = medium(cx, &density, 0.2, 1.5);
        let initial = [[0.5 * COUNT, 0.5 * COUNT, 0.0], [COUNT, 0.0, 0.0]];
        let thermal = thermal(
            vec![cell(1, 1.0, 1.0, 0.0), cell(2, 1.0, 1.0, 0.0)],
            &SENSIBLE,
            vec![],
        );
        let geometry = PhotonRecyclingGeometry {
            maximum_total_segments: 16_384,
            ..geometry()
        };
        let sections = [[0.8 / COUNT, 0.0], [0.0; 2]];
        let heating = owner(cx, &medium, &initial, &thermal, sections, geometry, 128).unwrap();
        let replay = owner(cx, &medium, &initial, &thermal, sections, geometry, 128).unwrap();
        assert_eq!(heating.emission_source_identity(), replay.emission_source_identity());
        let (warm, accepted, audit) = thermal
            .advance_observed(&heating, 0.4, config(0.4), || Ok(()))
            .unwrap();
        let (_, repeated, _) = thermal
            .advance_observed(&replay, 0.4, config(0.4), || Ok(()))
            .unwrap();
        assert_eq!(accepted.energy(), repeated.energy());
        assert!(accepted.populations(1).unwrap()[1] > 0.0);
        assert!(accepted.energy().recycled_excitation_j > 0.0);
        assert!(accepted.energy().secondary_absorption_heat_j > 0.0);
        assert_eq!(accepted.absorbed_pump_energy_j(), [0.0; 2]);
        for cell in 0..2 {
            let populations = accepted.populations(cell).unwrap();
            assert!(populations.iter().all(|value| value.is_finite() && *value >= 0.0));
            relative(populations.iter().sum(), COUNT);
            assert_eq!(accepted.populations(cell), repeated.populations(cell));
        }
        relative(warm.energy_j(), accepted.thermal_input_j());
        relative(audit.external_heat_j, accepted.thermal_input_j());
        balance(&accepted);

        let transparent = owner(cx, &medium, &initial, &thermal, [[0.0; 2]; 2], geometry, 128).unwrap();
        let (_, reference, _) = thermal
            .advance_observed(&transparent, 0.4, config(0.4), || Ok(()))
            .unwrap();
        relative(reference.populations(0).unwrap()[1], initial[0][1] * (-0.4_f64).exp());
        assert_eq!(reference.populations(1).unwrap()[1], 0.0);
        assert_eq!(reference.energy().recycled_excitation_j, 0.0);
        assert!(reference.energy().secondary_absorption_heat_j > 0.0);
        balance(&reference);
    });
}

#[test]
fn g0_elastic_binding_refuses_unfinished_histories_and_aggregate_budgets() {
    with_cx(|cx| {
        let density = [1.0; 2];
        let medium = medium(cx, &density, 0.2, 1e6);
        let initial = [[0.5 * COUNT, 0.5 * COUNT, 0.0], [COUNT, 0.0, 0.0]];
        let thermal = thermal(
            vec![cell(1, 1.0, 1.0, 0.0), cell(2, 1.0, 1.0, 0.0)],
            &SENSIBLE,
            vec![],
        );
        for (geometry, events) in [
            (geometry(), 0),
            (PhotonRecyclingGeometry { maximum_paths: 127, ..geometry() }, 128),
            (PhotonRecyclingGeometry { maximum_total_segments: 1, ..geometry() }, 128),
        ] {
            assert!(owner(cx, &medium, &initial, &thermal, [[0.0; 2]; 2], geometry, events).is_err());
            assert_eq!(thermal.energy_j(), 0.0);
            assert_eq!(density, [1.0; 2]);
        }
    });
}
