//! G1/G3/G4: finite photon return, native heat and accepted active opacity.

use super::beams::{accounts, bind, ray, sections};
use super::*;
use fs_render::volumes::heating::luminescence::saturation::annihilation::{
    AnnihilationLimits, PhosphorescenceAnnihilation, TripletAnnihilation,
};
use fs_render::volumes::heating::luminescence::saturation::beam::DepletingPumpBeam;
use fs_render::volumes::heating::luminescence::saturation::recycling::{
    PhotonRecyclingGeometry, RecyclingEmitterCrossSections, SaturableRecyclingLimits,
};

mod scattering;

// This fixture emits at 800 nm, beyond the legacy visible constant spectrum.
// Its passive host is explicitly authored across both emission nodes.
fn passive_spectrum(coefficient: f64) -> MassAbsorptionSpectrum {
    MassAbsorptionSpectrum::source_table(
        &[(380.0, coefficient), (800.0, coefficient)],
        ContentHash([81; 32]),
    )
    .unwrap()
}

fn geometry() -> PhotonRecyclingGeometry {
    PhotonRecyclingGeometry {
        histories_per_emitter: 16,
        maximum_paths: 128,
        maximum_total_segments: 512,
        seed: 37,
        source_id: 38,
        minimum_temperature_k: 290.0,
        maximum_temperature_k: 800.0,
        fraction_roundoff: 1e-12,
    }
}

fn controls() -> SaturableRecyclingLimits {
    SaturableRecyclingLimits {
        pump: super::beams::limits(),
        annihilation: AnnihilationLimits {
            relative_tolerance: 2e-5,
            absolute_energy_tolerance_j: 1e-10,
            maximum_substeps: 2048,
        },
        maximum_feedback_iterations: 96,
        maximum_secant_history: 8,
    }
}

fn sources() -> (LuminescenceSource<2>, SingletTripletSource<2>) {
    let make = |level, lifetime, weights, id| {
        LuminescenceSource::try_new(
            FluorescenceSource::try_new(
                ContentHash([id; 32]),
                (380.0, 650.0),
                (290.0, 800.0),
                1.0,
                [650.0, 800.0],
                weights,
            )
            .unwrap(),
            level,
            lifetime,
        )
        .unwrap()
    };
    (
        make(650.0, 1.0, [1.0, 0.0], 73),
        SingletTripletSource::try_new(
            make(800.0, 5.0, [0.0, 1.0], 74),
            0.0,
            ArrheniusQuenching::try_new(0.0, 0.0).unwrap(),
            SingletTripletLimits::default(),
        )
        .unwrap(),
    )
}

fn bindings(populations: &[[f64; 3]]) -> Vec<SaturableEmitterCell> {
    let pump =
        PrescribedPhotonPump::try_new(ContentHash([75; 32]), [400.0, 450.0], [0.0; 2], [0.0; 2])
            .unwrap();
    populations
        .iter()
        .enumerate()
        .map(|(grid_cell, &initial_populations)| SaturableEmitterCell {
            grid_cell,
            thermal_cell: CellId(grid_cell as u64 + 1),
            initial_populations,
            pump,
        })
        .collect()
}

fn molecule(cells: usize, values: [[f64; 2]; 2]) -> Vec<RecyclingEmitterCrossSections<2>> {
    (0..cells)
        .map(|grid_cell| {
            RecyclingEmitterCrossSections::authored(grid_cell, ContentHash([76; 32]), values)
                .unwrap()
        })
        .collect()
}

fn owner(
    cx: &Cx<'_>,
    medium: &AbsorbingGrid<'_>,
    populations: &[[f64; 3]],
    values: [[f64; 2]; 2],
) -> (
    DepositionSimulation<'static>,
    SaturableLuminescentHeating<2>,
) {
    let thermal = thermal(
        (0..populations.len())
            .map(|index| cell(index as u64 + 1, 1.0, 1.0, 0.0))
            .collect(),
        &SENSIBLE,
        vec![],
    );
    let (singlet, triplet) = sources();
    let heating = SaturableLuminescentHeating::try_new_with_photon_recycling(
        cx,
        medium,
        singlet,
        triplet,
        &bindings(populations),
        &thermal,
        populations.len(),
        1e-7,
        &molecule(populations.len(), values),
        geometry(),
        controls(),
        || Ok(()),
    )
    .unwrap();
    (thermal, heating)
}

#[test]
fn g1_zero_active_return_retains_analytic_decay_and_routes_passive_host_heat_once() {
    with_cx(|cx| {
        for coefficient in [0.0, 0.4] {
            let density = [1.0];
            let medium = grid(
                cx,
                &density,
                passive_spectrum(coefficient),
            );
            let initial = [0.5 * COUNT, 0.3 * COUNT, 0.2 * COUNT];
            let (thermal, heating) = owner(cx, &medium, &[initial], [[0.0; 2]; 2]);
            let (warm, accepted, audit) = thermal
                .advance_observed(&heating, 0.4, config(0.4), || Ok(()))
                .unwrap();
            relative(
                accepted.populations(0).unwrap()[1],
                initial[1] * (-0.4_f64).exp(),
            );
            relative(
                accepted.populations(0).unwrap()[2],
                initial[2] * (-0.08_f64).exp(),
            );
            assert_eq!(accepted.absorbed_pump_energy_j(), [0.0; 2]);
            assert_eq!(accepted.energy().recycled_excitation_j, 0.0);
            assert_eq!(accepted.energy().passive_pump_absorption_j, 0.0);
            if coefficient == 0.0 {
                relative(
                    accepted.energy().escaped_j[0],
                    initial[1] * -(-0.4_f64).exp_m1() * HC / 650.0,
                );
                relative(
                    accepted.energy().escaped_j[1],
                    initial[2] * -(-0.08_f64).exp_m1() * HC / 800.0,
                );
                assert_eq!(accepted.thermal_input_j(), 0.0);
            } else {
                assert!(accepted.energy().secondary_absorption_heat_j > 0.0);
                relative(
                    accepted.thermal_input_j(),
                    accepted.energy().secondary_absorption_heat_j,
                );
                relative(warm.energy_j(), accepted.thermal_input_j());
                relative(audit.external_heat_j, accepted.thermal_input_j());
            }
            balance(&accepted);
        }
    });
}

#[test]
fn g3_recycled_light_excites_an_unpumped_cell_and_endpoint_camera_uses_its_ground_opacity() {
    with_cx(|cx| {
        let density = [1.0; 2];
        let medium = grid(cx, &density, passive_spectrum(0.0));
        let initial = [[0.5 * COUNT, 0.5 * COUNT, 0.0], [COUNT, 0.0, 0.0]];
        let sigma = 4.0 / COUNT;
        let (thermal, heating) = owner(cx, &medium, &initial, [[sigma, 0.0], [0.0; 2]]);
        let (_, transparent) = owner(cx, &medium, &initial, [[0.0; 2]; 2]);
        let (_, accepted, _) = thermal
            .advance_observed(&heating, 0.4, config(0.4), || Ok(()))
            .unwrap();
        let (_, reference, _) = thermal
            .advance_observed(&transparent, 0.4, config(0.4), || Ok(()))
            .unwrap();
        assert!(accepted.populations(1).unwrap()[1] > 0.0);
        assert_eq!(reference.populations(1).unwrap()[1], 0.0);
        assert!(accepted.energy().recycled_excitation_j > 0.0);
        assert!(accepted.energy().stored_excitation_j > reference.energy().stored_excitation_j);
        assert!(accepted.energy().escaped_j[0] < reference.energy().escaped_j[0]);
        assert_eq!(accepted.absorbed_pump_energy_j(), [0.0; 2]);
        let view = accepted.recycling_emission_grid(cx, &medium, 2).unwrap();
        let observed = view
            .transport(
                cx,
                Point3::new(-1.0, 0.5, 0.5),
                Vec3::new(1.0, 0.0, 0.0),
                f64::INFINITY,
            )
            .unwrap();
        let ground: f64 = (0..2)
            .map(|cell| accepted.populations(cell).unwrap()[0])
            .sum();
        relative(observed.transmittance[0], (-sigma * ground).exp());
        assert_eq!(observed.transmittance[1], 1.0);
        let passive = LuminescenceGrid::try_new(
            cx,
            &medium,
            accepted.emission_source_identity(),
            [650.0, 800.0],
            accepted.cell_emission_w(),
            2,
        )
        .unwrap()
        .transport(
            cx,
            Point3::new(-1.0, 0.5, 0.5),
            Vec3::new(1.0, 0.0, 0.0),
            f64::INFINITY,
        )
        .unwrap();
        assert!(observed.emitted_radiance_w_m2_sr[0] < passive.emitted_radiance_w_m2_sr[0]);
        for cell in 0..2 {
            relative(accepted.populations(cell).unwrap().iter().sum(), COUNT);
        }
        balance(&accepted);
    });
}

#[test]
fn g3_shared_pumps_pairs_recycling_and_retry_keep_one_native_inventory_and_external_ledger() {
    with_cx(|cx| {
        let density = [1.0];
        let medium = grid(cx, &density, passive_spectrum(0.0));
        let (thermal, base) = owner(
            cx,
            &medium,
            &[[COUNT, 0.0, 0.0]],
            [[0.0, 2.0 / COUNT], [0.0; 2]],
        );
        let cross_sections = sections(1, [[1.0 / COUNT, 0.0], [0.0; 2]]);
        let beams = [false, true].map(|reverse| DepletingPumpBeam {
            source_identity: ContentHash([77; 32]),
            beam: ray(1, reverse, [400.0, 450.0], [COUNT * HC / 400.0, 0.0]),
            cross_sections: &cross_sections,
            pump_window_s: None,
        });
        // Adding external beams after recycling must rebuild the shared model,
        // while pair laws remain attached to the same original inventory.
        let heating = bind(cx, &medium, &thermal, base, &beams)
            .with_triplet_annihilation(
                &[PhosphorescenceAnnihilation {
                    grid_cell: 0,
                    law: TripletAnnihilation::try_new(1.0 / COUNT, 1.0, 0.5).unwrap(),
                }],
                controls().annihilation,
            )
            .unwrap()
            .with_pump_window(0.0, 0.2)
            .unwrap();
        let mut candidate = heating.clone();
        candidate
            .step_sources(&thermal, 0.0, 0.2, 0.2, &mut || Ok(()))
            .unwrap();
        let retried = candidate
            .step_sources(&thermal, 0.0, 0.1, 0.1, &mut || Ok(()))
            .unwrap();
        let fresh = heating
            .clone()
            .step_sources(&thermal, 0.0, 0.1, 0.1, &mut || Ok(()))
            .unwrap();
        assert_eq!(retried, fresh);
        assert_eq!(candidate.populations(0), heating.populations(0));
        assert_eq!(candidate.energy(), heating.energy());
        assert_eq!(
            candidate.depleting_pump_report(),
            heating.depleting_pump_report()
        );
        let mut polls = 0;
        assert!(matches!(
            thermal.advance_observed(&heating, 0.2, config(0.2), || {
                polls += 1;
                if polls > 40 {
                    Err(DepositionError::Cancelled)
                } else {
                    Ok(())
                }
            }),
            Err(OpticalHeatingError::Thermal(DepositionError::Cancelled))
        ));
        let (warm, accepted, audit) = thermal
            .advance_observed(&heating, 0.2, config(0.2), || Ok(()))
            .unwrap();
        assert!(accepted.energy().recycled_excitation_j > 0.0);
        assert!(accepted.energy().secondary_absorption_heat_j > 0.0);
        assert!(accepted.energy().annihilation_heat_j > 0.0);
        assert!(accepted.populations(0).unwrap()[2] > 0.0);
        assert_eq!(accepted.absorbed_pump_energy_j()[1], 0.0);
        relative(audit.external_heat_j, accepted.thermal_input_j());
        accounts(&accepted);
        let (_, dark, _) = warm
            .advance_observed(&accepted, 0.4, config(0.2), || Ok(()))
            .unwrap();
        assert_eq!(
            dark.absorbed_pump_energy_j(),
            accepted.absorbed_pump_energy_j()
        );
        assert!(dark.energy().recycled_excitation_j > accepted.energy().recycled_excitation_j);
        accounts(&dark);

        let (singlet, triplet) = sources();
        let mut exhausted = geometry();
        exhausted.maximum_paths = 1;
        assert!(
            SaturableLuminescentHeating::try_new_with_photon_recycling(
                cx,
                &medium,
                singlet,
                triplet,
                &bindings(&[[COUNT, 0.0, 0.0]]),
                &thermal,
                1,
                1e-7,
                &molecule(1, [[0.0; 2]; 2]),
                exhausted,
                controls(),
                || Ok(()),
            )
            .is_err()
        );
    });
}

#[test]
fn g1_prescribed_excitation_with_transparent_recycling_keeps_analytic_local_accounts() {
    with_cx(|cx| {
        let density = [1.0];
        let medium = grid(cx, &density, passive_spectrum(0.0));
        let thermal = thermal(vec![cell(1, 1.0, 1.0, 0.0)], &SENSIBLE, vec![]);
        let mut binding = bindings(&[[COUNT, 0.0, 0.0]]);
        binding[0].pump = PrescribedPhotonPump::try_new(
            ContentHash([78; 32]),
            [400.0, 450.0],
            [1e-20; 2],
            [0.0, 1.5e20],
        )
        .unwrap();
        let (singlet, triplet) = sources();
        let heating = SaturableLuminescentHeating::try_new_with_photon_recycling(
            cx,
            &medium,
            singlet,
            triplet,
            &binding,
            &thermal,
            1,
            1e-7,
            &molecule(1, [[0.0; 2]; 2]),
            geometry(),
            controls(),
            || Ok(()),
        )
        .unwrap()
        .with_pump_window(0.0, 0.2)
        .unwrap();
        let (warm, charged, audit) = thermal
            .advance_observed(&heating, 0.2, config(0.2), || Ok(()))
            .unwrap();
        let duration = 0.2_f64;
        let rate = 1.5;
        let decay = 0.2;
        let a = rate + decay;
        let fraction = -(-a * duration).exp_m1();
        let expected_triplet = COUNT * rate / a * fraction;
        let integrated_triplets = COUNT * rate / a * (duration - fraction / a);
        let local_events = rate * (COUNT * duration - integrated_triplets);
        relative(charged.populations(0).unwrap()[2], expected_triplet);
        relative(
            charged.absorbed_pump_energy_j()[1],
            local_events * HC / 450.0,
        );
        assert_eq!(charged.absorbed_pump_energy_j()[0], 0.0);
        assert_eq!(charged.energy().recycled_excitation_j, 0.0);
        assert!(charged.depleting_pump_report().is_none());
        relative(audit.external_heat_j, charged.thermal_input_j());
        balance(&charged);
        let (_, dark, _) = warm
            .advance_observed(&charged, 0.4, config(0.2), || Ok(()))
            .unwrap();
        relative(
            dark.populations(0).unwrap()[2],
            expected_triplet * (-0.04_f64).exp(),
        );
        assert_eq!(
            dark.absorbed_pump_energy_j(),
            charged.absorbed_pump_energy_j()
        );
        balance(&dark);
    });
}

#[test]
fn g3_local_triplet_excitation_competes_with_a_singlet_beam_pairs_and_photon_return() {
    with_cx(|cx| {
        for recycling in [false, true] {
            let density = [1.0];
            let medium = grid(cx, &density, passive_spectrum(0.0));
            let thermal = thermal(vec![cell(1, 1.0, 1.0, 0.0)], &SENSIBLE, vec![]);
            let mut binding = bindings(&[[COUNT, 0.0, 0.0]]);
            binding[0].pump = PrescribedPhotonPump::try_new(
                ContentHash([79; 32]),
                [400.0, 450.0],
                [1e-20; 2],
                [0.0, 1.5e20],
            )
            .unwrap();
            let (singlet, triplet) = sources();
            let mut base = SaturableLuminescentHeating::try_new(
                &medium,
                singlet,
                triplet,
                &binding,
                &thermal,
                1,
                1e-7,
                None,
                || Ok(()),
            )
            .unwrap();
            if recycling {
                base = base
                    .with_photon_recycling(
                        cx,
                        &medium,
                        &molecule(1, [[0.0, 2.0 / COUNT], [0.0; 2]]),
                        &thermal,
                        geometry(),
                        controls(),
                        || Ok(()),
                    )
                    .unwrap();
            }
            let cross_sections = sections(1, [[1.0 / COUNT, 0.0], [0.0; 2]]);
            let power = COUNT * HC / 400.0;
            let heating = bind(
                cx,
                &medium,
                &thermal,
                base,
                &[DepletingPumpBeam {
                    source_identity: ContentHash([80; 32]),
                    beam: ray(1, false, [400.0, 450.0], [power, 0.0]),
                    cross_sections: &cross_sections,
                    pump_window_s: Some((0.0, 0.1)),
                }],
            )
            .with_triplet_annihilation(
                &[PhosphorescenceAnnihilation {
                    grid_cell: 0,
                    law: TripletAnnihilation::try_new(1.0 / COUNT, 1.0, 0.5).unwrap(),
                }],
                controls().annihilation,
            )
            .unwrap()
            .with_pump_window(0.0, 0.2)
            .unwrap();
            let mut candidate = heating.clone();
            candidate
                .step_sources(&thermal, 0.0, 0.08, 0.08, &mut || Ok(()))
                .unwrap();
            let retry = candidate
                .step_sources(&thermal, 0.0, 0.04, 0.04, &mut || Ok(()))
                .unwrap();
            let fresh = heating
                .clone()
                .step_sources(&thermal, 0.0, 0.04, 0.04, &mut || Ok(()))
                .unwrap();
            assert_eq!(retry, fresh);
            assert_eq!(candidate.populations(0), heating.populations(0));
            assert_eq!(candidate.energy(), heating.energy());
            assert_eq!(
                candidate.depleting_pump_report(),
                heating.depleting_pump_report()
            );
            let (warm, charged, audit) = thermal
                .advance_observed(&heating, 0.2, config(0.2), || Ok(()))
                .unwrap();
            let beam = &charged.depleting_pump_report().unwrap().beams[0];
            relative(beam.incident_j[0], 0.1 * power);
            relative(
                beam.incident_j[0],
                beam.transmitted_j[0] + beam.absorbed_j[0],
            );
            // Only the declared beam supplies external G->S photons. Only
            // the maintained local field supplies external G->T photons;
            // recycled G->T photons must not enter either external account.
            relative(charged.absorbed_pump_energy_j()[0], beam.absorbed_j[0]);
            assert!(charged.absorbed_pump_energy_j()[1] > 0.0);
            assert!(charged.energy().annihilation_heat_j > 0.0);
            assert!(charged.populations(0).unwrap()[2] > 0.0);
            if recycling {
                assert!(charged.energy().recycled_excitation_j > 0.0);
                assert!(charged.energy().secondary_absorption_heat_j > 0.0);
            } else {
                assert_eq!(charged.energy().recycled_excitation_j, 0.0);
            }
            relative(audit.external_heat_j, charged.thermal_input_j());
            balance(&charged);
            let (_, dark, _) = warm
                .advance_observed(&charged, 0.4, config(0.2), || Ok(()))
                .unwrap();
            assert_eq!(
                dark.absorbed_pump_energy_j(),
                charged.absorbed_pump_energy_j()
            );
            assert_eq!(
                dark.depleting_pump_report().unwrap().beams[0].incident_j,
                beam.incident_j
            );
            balance(&dark);
        }
    });
}
