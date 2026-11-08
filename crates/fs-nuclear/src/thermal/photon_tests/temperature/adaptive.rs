//! G1/G3/G4: native time control from actual absorbed photon energy.

use super::*;

fn adaptive_controls(energy_j: f64) -> NuclearAdaptiveControls {
    NuclearAdaptiveControls {
        step: controls(1.0),
        initial_dt_s: 1.0,
        minimum_dt_s: 1e-6,
        maximum_trials: 1000,
        maximum_substeps: 2000,
        tolerances: NuclearTimeTolerances {
            // Deliberately loose: only the transported photon joules can drive
            // refinement of this otherwise exactly integrated decay source.
            specific_energy_j_kg: 1e6,
            temperature_k: 1e6,
            nuclei: 1e5,
            nuclear_energy_j: energy_j,
            relative: 0.0,
            ..NuclearTimeTolerances::default()
        },
    }
}

fn owner(thermal: &DepositionSimulation<'_>, growing: bool) -> NuclearHeating<DecayCellSource> {
    let mut layers = temperature_layers(growing, 350.0);
    let rear = material(3, [4.0, 4.0], false);
    // Rear optical depth 16: the total absorbed energy barely changes when
    // front shielding changes, but the energy partition between cells does.
    layers[2].material = Arc::new(
        TemperaturePhotonMassMaterial::new(
            "manufactured opaque rear absorber".into(),
            [300.0, 350.0]
                .map(|temperature_k| PhotonMassTemperatureKnot {
                    temperature_k,
                    material: rear.as_ref().clone(),
                })
                .to_vec(),
            PhotonSlabLimits::default(),
            || Ok(()),
        )
        .unwrap(),
    );
    let partition = partition();
    plain(thermal, &partition)
        .with_temperature_photon_slab(
            thermal,
            "manufactured adaptive radioactive shielding".into(),
            geometry(),
            layers,
            vec![DecayRadiationEmitter {
                cell_id: CellId(1),
                partition,
            }],
            vec![PhotonOrdinate {
                mu: 1.0,
                weight: 1.0,
            }],
            PhotonSlabLimits::default(),
            None,
            PhotonFeedbackBudget::default(),
            || Ok(()),
        )
        .unwrap()
        .0
}

#[test]
fn g1_adaptive_photon_absorption_controls_receiver_heat_with_loose_thermal_scales() {
    let original = thermal_grid(cells(1.0));
    let source = owner(&original, true);
    let photons = emitted_j(1.0);
    let exact_theta = ((-0.5_f64).exp() - (-0.5_f64).exp_m1() * (photons / 8.0).exp()).ln() / 0.5;
    let (coarse_thermal, coarse, _) = source
        .advance_to(&original, 1.0, controls(1.0), || Ok(()))
        .unwrap();
    let (fine_thermal, fine, _) = source
        .advance_to(&original, 1.0, controls(0.5), || Ok(()))
        .unwrap();
    assert!(
        (coarse.audit().radiation_deposited_j - fine.audit().radiation_deposited_j).abs() < 1e-6
    );
    assert!(
        2.0 * (coarse_thermal.temperature_k(CellId(2)).unwrap()
            - fine_thermal.temperature_k(CellId(2)).unwrap())
        .abs()
            > 0.2
    );
    let mut errors = Vec::new();
    for tolerance in [1e-3, 2.5e-4] {
        let (heated, evolved, report) = source
            .advance_adaptive_to(&original, 1.0, adaptive_controls(tolerance), || Ok(()))
            .unwrap();
        let theta = heated.temperature_k(CellId(2)).unwrap() - 300.0;
        errors.push((theta - exact_theta).abs());
        assert!(report.accuracy_rejections > 0);
        assert!(report.accepted_comparisons > 1);
        assert!(report.maximum_accepted_normalized_error <= 1.0);
        assert_eq!(report.thermal.substeps, 2 * report.accepted_comparisons);
        assert_eq!(evolved.audit().substeps, report.thermal.substeps);
        let shield_heat = 2.0 * theta;
        let rear_heat = (0.5 * photons - shield_heat) * -(-16.0_f64).exp_m1();
        close(
            heated.temperature_k(CellId(3)).unwrap(),
            300.0 + rear_heat / 16.0,
            2e-11,
        );
        close(evolved.audit().released.radiation_j, photons, 2e-10);
        close(
            evolved.audit().radiation_deposited_j,
            shield_heat + rear_heat,
            2e-10,
        );
        close(
            evolved.audit().radiation_escaped_j,
            photons - shield_heat - rear_heat,
            2e-10,
        );
        close(
            report.thermal.external_heat_j,
            shield_heat + rear_heat,
            2e-10,
        );
        let populations = evolved.source(CellId(1)).unwrap().inventory.populations();
        close(populations[0], INITIAL * 0.5, 100.0);
        close(
            populations[1],
            INITIAL * 0.5 * std::f64::consts::LN_2,
            100.0,
        );
        let last = evolved.last_radiation_deposition().unwrap();
        let release =
            <DecayCellSource as NuclearHeatSource>::released(evolved.last_step(CellId(1)).unwrap());
        close(last.emitted_j, release.radiation_j, 2e-12);
        assert!(last.emitted_j < photons);
    }
    assert!(
        errors[1] < 0.008 && errors[1] < 0.65 * errors[0],
        "{errors:?}"
    );
}

#[test]
fn g3_adaptive_photon_absorption_keeps_fixed_opacity_event_cuts_and_interval_energy() {
    let original = thermal(
        cells(1.0),
        &[
            (MaterialId(1), &HOT),
            (MaterialId(2), &COOL),
            (MaterialId(3), &HOT),
        ],
        vec![TimedEvent {
            time_s: 0.375,
            event: Event::Boundary(
                CellId(2),
                Boundary {
                    ambient_k: 300.0,
                    htc_w_m2_k: 0.0,
                    reservoir: None,
                    source_w: 2.0,
                },
            ),
        }],
    );
    let source = owner(&original, false);
    let (heated, evolved, report) = source
        .advance_adaptive_to(&original, 1.0, adaptive_controls(1e-7), || Ok(()))
        .unwrap();
    assert!(report.native_interval_cuts > 0);
    assert_eq!(report.accuracy_rejections, 0);
    assert_eq!(report.accepted_comparisons, 2);
    assert_eq!(report.thermal.substeps, 4);
    assert_eq!(report.thermal.events, 1);
    let photons = emitted_j(1.0);
    let shield = 0.5 * photons * -(-0.5_f64).exp_m1();
    let rear = 0.5 * photons * (-0.5_f64).exp() * -(-16.0_f64).exp_m1();
    let boundary_heat = 2.0 * (1.0 - 0.375);
    close(
        heated.temperature_k(CellId(2)).unwrap(),
        300.0 + (shield + boundary_heat) / 2.0,
        2e-12,
    );
    close(
        heated.temperature_k(CellId(3)).unwrap(),
        300.0 + rear / 16.0,
        2e-12,
    );
    close(evolved.audit().radiation_deposited_j, shield + rear, 2e-11);
    close(evolved.audit().thermal_credit_j, shield + rear, 2e-11);
    close(
        report.thermal.external_heat_j,
        shield + rear + boundary_heat,
        2e-11,
    );
    // The final native interval is [0.6875, 1]. The accepted advance account
    // includes all four physical intervals, while this observation includes
    // only the final one and excludes the independent boundary heat.
    let last_photons = emitted_j(1.0) - emitted_j(0.6875);
    let last = evolved.last_radiation_deposition().unwrap();
    close(last.emitted_j, last_photons, 2e-11);
    close(
        last.receiver_energy_j[&2],
        0.5 * last_photons * -(-0.5_f64).exp_m1(),
        2e-11,
    );
    close(
        last.receiver_energy_j[&3],
        0.5 * last_photons * (-0.5_f64).exp() * -(-16.0_f64).exp_m1(),
        2e-11,
    );
    assert!(last.deposited_j < evolved.audit().radiation_deposited_j);
}

#[test]
fn g4_adaptive_photon_absorption_refusal_preserves_accepted_transport_and_nuclei() {
    let original = thermal_grid(cells(1.0));
    let source = owner(&original, true);
    let (middle, accepted, _) = source
        .advance_to(&original, 0.25, controls(0.25), || Ok(()))
        .unwrap();
    let radiation = accepted.last_radiation_deposition().unwrap().clone();
    let response = Arc::clone(accepted.grouped_radiation_response().unwrap());
    let temperatures = accepted.prepared_photon_temperatures_k().unwrap().to_vec();
    let work = accepted.photon_feedback_work();
    let populations = accepted
        .source(CellId(1))
        .unwrap()
        .inventory
        .populations()
        .to_vec();
    let audit = accepted.audit();
    let mut controls = adaptive_controls(1e-6);
    controls.minimum_dt_s = 1.0;
    let failure = accepted
        .advance_adaptive_to(&middle, 1.25, controls, || Ok(()))
        .err()
        .unwrap();
    assert!(
        matches!(failure, NuclearAdaptiveError::TimeAccuracy { error, .. }
        if error.nuclear_energy_difference_j > 1e-3 && error.normalized > 1.0)
    );
    assert_eq!(accepted.time_s(), 0.25);
    assert_eq!(middle.time_s(), 0.25);
    assert_eq!(
        accepted.source(CellId(1)).unwrap().inventory.populations(),
        populations.as_slice()
    );
    assert_eq!(accepted.audit(), audit);
    assert_eq!(accepted.last_radiation_deposition(), Some(&radiation));
    assert_eq!(accepted.photon_feedback_work(), work);
    assert_eq!(
        accepted.prepared_photon_temperatures_k().unwrap(),
        temperatures.as_slice()
    );
    assert!(Arc::ptr_eq(
        &response,
        accepted.grouped_radiation_response().unwrap()
    ));
}
