use super::*;
use fs_blake3::ContentHash;
use fs_render::volumes::heating::fluorescence::{FluorescenceSource, FluorescentOpticalHeating};
use fs_render::volumes::heating::illumination::schedule::{
    SpectralPowerSchedule, SpectralPowerSegment,
};
use fs_render::volumes::heating::luminescence::{LuminescenceSource, LuminescentOpticalHeating};
use fs_render::volumes::transport::luminescence::LuminescenceGrid;
use fs_thermal::enthalpy::DepositionObserver;
use fs_thermal::luminescence::ArrheniusQuenching;

fn source(lifetime: f64, yield_: f64) -> LuminescenceSource<2> {
    LuminescenceSource::try_new(
        FluorescenceSource::try_new(
            ContentHash([111; 32]),
            (380.0, 500.0),
            (290.0, 800.0),
            yield_,
            [600.0, 800.0],
            [1.0, 3.0],
        )
        .unwrap(),
        500.0,
        lifetime,
    )
    .unwrap()
}

fn spectrum(absorption: f64) -> MassAbsorptionSpectrum {
    MassAbsorptionSpectrum::source_table(
        &[
            (380.0, absorption),
            (500.0, absorption),
            (600.0, 0.0),
            (800.0, 0.0),
        ],
        ContentHash([112; 32]),
    )
    .unwrap()
}

fn setup(
    cx: &Cx<'_>,
    source: LuminescenceSource<2>,
) -> (
    DepositionSimulation<'static>,
    LuminescentOpticalHeating<1, 2>,
) {
    let density = [1.0];
    let medium = grid(cx, &density, spectrum(1e4));
    let deposition = medium
        .deposit_beam_power(
            cx,
            Point3::new(0.0, 0.5, 0.5),
            Vec3::new(1.0, 0.0, 0.0),
            1.0,
            &[400.0],
            &[10.0],
            1,
        )
        .unwrap();
    let thermal = thermal(vec![cell(1, 1.0, 1.0, 0.0)], &SENSIBLE, vec![]);
    let heating = LuminescentOpticalHeating::try_new(
        &medium,
        source,
        deposition,
        &[(0, CellId(1))],
        &thermal,
        1,
        1e-12,
        1e-9,
        || Ok(()),
    )
    .unwrap();
    (thermal, heating)
}

fn pulse(end_s: f64) -> SpectralPowerSchedule<1> {
    SpectralPowerSchedule::try_new(
        &[SpectralPowerSegment {
            start_s: 0.0,
            end_s,
            start_power_w: [10.0],
            end_power_w: [10.0],
        }],
        1,
    )
    .unwrap()
}

fn assert_partition(heating: &LuminescentOpticalHeating<1, 2>) {
    let e = heating.energy();
    let heat = e.prompt_heat_j + e.nonradiative_heat_j + e.relaxation_heat_j;
    close(heat, heating.pump_energy().thermal_input_j);
    close(
        heating.pump_energy().absorbed_j[0],
        e.escaped_j.iter().sum::<f64>() + heat + e.stored_excitation_j,
    );
    assert!(e.residual_j.abs() < 1e-10);
    assert!(e.absolute_roundoff_j < 1e-10);
}

#[test]
fn g1_continuous_pump_has_an_exact_excited_store_and_separate_heat_channels() {
    with_cx(|cx| {
        let (thermal, heating) = setup(cx, source(2.0, 0.75));
        let (next, accepted, audit) = thermal
            .advance_observed(&heating, 1.0, config(1.0), || Ok(()))
            .unwrap();
        let stored = 16.0 * (1.0 - (-0.5_f64).exp());
        let released = 8.0 - stored;
        let emitted = [
            released * 0.75 * 0.25 * 500.0 / 600.0,
            released * 0.75 * 0.75 * 500.0 / 800.0,
        ];
        close(accepted.energy().stored_excitation_j, stored);
        close(accepted.stored_energy_j(0).unwrap(), stored);
        assert_eq!(accepted.stored_energy_j(1), None);
        close(accepted.energy().prompt_heat_j, 2.0);
        close(accepted.energy().nonradiative_heat_j, 0.25 * released);
        for (actual, expected) in accepted.energy().escaped_j.iter().zip(emitted) {
            close(*actual, expected);
        }
        close(
            accepted.energy().relaxation_heat_j,
            0.75 * released - emitted.iter().sum::<f64>(),
        );
        close(
            audit.external_heat_j,
            accepted.pump_energy().thermal_input_j,
        );
        close(
            next.temperature_k(CellId(1)).unwrap(),
            300.0 + audit.external_heat_j / 4.0,
        );
        // The thermal state never receives the excitation store as heat.
        assert!(audit.external_heat_j < 10.0 - emitted.iter().sum::<f64>());
        assert_partition(&accepted);
        assert_eq!(heating.energy().stored_excitation_j, 0.0);
        assert_eq!(thermal.temperature_k(CellId(1)).unwrap(), 300.0);
    });
}

#[test]
fn g1_pulse_edge_is_exact_and_dark_afterglow_releases_stored_excitation() {
    with_cx(|cx| {
        let (thermal, heating) = setup(cx, source(2.0, 0.75));
        let heating = heating.with_pump_schedule(pulse(1.0)).unwrap();
        let (next, accepted, audit) = thermal
            .advance_observed(&heating, 10.0, config(10.0), || Ok(()))
            .unwrap();
        assert_eq!(audit.substeps, 2); // exactly [0,1] and [1,10], no lifetime CFL
        close(accepted.pump_energy().incident_j[0], 10.0);
        let stored = 16.0 * (1.0 - (-0.5_f64).exp()) * (-4.5_f64).exp();
        close(accepted.energy().stored_excitation_j, stored);
        close(accepted.energy().prompt_heat_j, 2.0);
        let released = 8.0 - stored;
        close(accepted.energy().nonradiative_heat_j, released * 0.25);
        close(
            accepted.cell_emission_w()[0].1[0],
            stored / 2.0 * 0.75 * 0.25 * 500.0 / 600.0,
        );
        assert!(accepted.cell_emission_w()[0].1[0] > 0.0);
        assert!(accepted.last_mean_cell_emission_w()[0].1[0] > accepted.cell_emission_w()[0].1[0]);
        let density = [1.0];
        let medium = grid(cx, &density, spectrum(1e4));
        let view = LuminescenceGrid::try_new(
            cx,
            &medium,
            accepted.source().fluorescence().source_identity(),
            *accepted.source().fluorescence().wavelengths_nm(),
            accepted.cell_emission_w(),
            1,
        )
        .unwrap();
        let ray = view
            .transport(
                cx,
                Point3::new(-1.0, 0.5, 0.5),
                Vec3::new(1.0, 0.0, 0.0),
                3.0,
            )
            .unwrap();
        close(
            ray.emitted_radiance_w_m2_sr[0],
            accepted.cell_emission_w()[0].1[0] / (4.0 * core::f64::consts::PI),
        );
        close(next.energy_j(), accepted.pump_energy().thermal_input_j);
        assert_partition(&accepted);
    });
}

#[test]
fn g3_partitioning_and_restart_preserve_pulse_fluence_heat_light_and_store() {
    with_cx(|cx| {
        let (thermal, heating) = setup(cx, source(0.3, 0.9));
        let heating = heating.with_pump_schedule(pulse(0.75)).unwrap();
        let (whole, one, _) = thermal
            .advance_observed(&heating, 2.0, config(2.0), || Ok(()))
            .unwrap();
        let (fine, many, _) = thermal
            .advance_observed(&heating, 2.0, config(0.03125), || Ok(()))
            .unwrap();
        let (middle, saved, _) = thermal
            .advance_observed(&heating, 0.5, config(0.125), || Ok(()))
            .unwrap();
        let (resumed, replay, _) = middle
            .advance_observed(&saved, 2.0, config(0.125), || Ok(()))
            .unwrap();
        for other in [&many, &replay] {
            close(
                other.energy().stored_excitation_j,
                one.energy().stored_excitation_j,
            );
            close(
                other.pump_energy().thermal_input_j,
                one.pump_energy().thermal_input_j,
            );
            for (a, b) in other.energy().escaped_j.iter().zip(one.energy().escaped_j) {
                close(*a, b);
            }
            close(other.pump_energy().absorbed_j[0], 7.5);
            assert_partition(other);
        }
        close(whole.energy_j(), fine.energy_j());
        close(whole.energy_j(), resumed.energy_j());
        let (_, duplicate, _) = middle
            .advance_observed(&saved, 2.0, config(0.125), || Ok(()))
            .unwrap();
        assert_eq!(duplicate.energy(), replay.energy());
        assert_eq!(duplicate.cell_emission_w(), replay.cell_emission_w());
    });
}

#[test]
fn g1_short_lifetime_recovers_existing_prompt_fluorescence_without_short_steps() {
    with_cx(|cx| {
        let density = [1.0];
        let medium = grid(cx, &density, spectrum(1e4));
        let deposition = medium
            .deposit_beam_power(
                cx,
                Point3::new(0.0, 0.5, 0.5),
                Vec3::new(1.0, 0.0, 0.0),
                1.0,
                &[400.0],
                &[10.0],
                1,
            )
            .unwrap();
        let (thermal, delayed) = setup(cx, source(1e-10, 0.8));
        let prompt = FluorescentOpticalHeating::try_new(
            &medium,
            *delayed.source().fluorescence(),
            deposition,
            &[(0, CellId(1))],
            &thermal,
            1,
            1e-12,
            || Ok(()),
        )
        .unwrap();
        let (_, delayed, audit) = thermal
            .advance_observed(&delayed, 1.0, config(1.0), || Ok(()))
            .unwrap();
        assert_eq!(audit.substeps, 1);
        close(delayed.energy().stored_excitation_j, 8e-10);
        for (a, b) in delayed
            .energy()
            .escaped_j
            .iter()
            .zip(prompt.power().escaped_w)
        {
            assert!((*a - b).abs() < 1e-9);
        }
        assert!(
            (delayed.pump_energy().thermal_input_j - prompt.pump_power().thermal_input_w).abs()
                < 1e-9
        );
        assert_partition(&delayed);
    });
}

#[test]
fn g1_quenching_changes_lifetime_and_yield_together_and_heats_the_native_cell() {
    with_cx(|cx| {
        let base = source(2.0, 0.75);
        let quenched = base.with_quenching(ArrheniusQuenching::try_new(1.5, 0.0).unwrap());
        let decay = quenched.decay_at_temperature(300.0).unwrap();
        close(decay.lifetime_s(), 0.5);
        close(decay.quantum_yield(), 0.1875);
        close(decay.quantum_yield() / decay.lifetime_s(), 0.75 / 2.0);
        let (thermal, a) = setup(cx, base);
        let (_, b) = setup(cx, quenched);
        let (_, a, _) = thermal
            .advance_observed(&a, 4.0, config(4.0), || Ok(()))
            .unwrap();
        let (hotter, b, audit) = thermal
            .advance_observed(&b, 4.0, config(4.0), || Ok(()))
            .unwrap();
        let stored = 4.0 * (1.0 - (-8.0_f64).exp());
        close(b.energy().stored_excitation_j, stored);
        close(
            b.energy().nonradiative_heat_j,
            (32.0 - stored) * (1.0 - 0.1875),
        );
        assert!(b.energy().escaped_j[0] < a.energy().escaped_j[0]);
        assert!(b.energy().nonradiative_heat_j > a.energy().nonradiative_heat_j);
        close(hotter.energy_j(), audit.external_heat_j);
        assert_partition(&b);
    });
}

#[test]
fn g0_zero_and_unit_yields_keep_dark_heat_and_photon_counts_distinct() {
    with_cx(|cx| {
        for yield_ in [0.0, 1.0] {
            let (thermal, heating) = setup(cx, source(0.1, yield_));
            let heating = heating.with_pump_schedule(pulse(0.2)).unwrap();
            let (_, accepted, _) = thermal
                .advance_observed(&heating, 1.0, config(1.0), || Ok(()))
                .unwrap();
            let e = accepted.energy();
            let decayed_j = 1.6 - e.stored_excitation_j;
            // Multiplication by wavelength converts relative energies back to
            // photon counts; h*c cancels from this independent number audit.
            let emitted_photon_energy_nm = e.escaped_j[0] * 600.0 + e.escaped_j[1] * 800.0;
            close(emitted_photon_energy_nm, decayed_j * 500.0 * yield_);
            close(e.nonradiative_heat_j, decayed_j * (1.0 - yield_));
            if yield_ == 0.0 {
                assert_eq!(accepted.cell_emission_w()[0].1, [0.0; 2]);
            }
            assert_partition(&accepted);
        }
    });
}

#[test]
fn g3_distinct_optical_stores_can_share_one_thermal_receiver() {
    with_cx(|cx| {
        let density = [1.0, 1.0];
        let medium = grid(cx, &density, spectrum(core::f64::consts::LN_2));
        let deposition = medium
            .deposit_beam_power(
                cx,
                Point3::new(0.0, 0.5, 0.5),
                Vec3::new(1.0, 0.0, 0.0),
                2.0,
                &[400.0],
                &[10.0],
                2,
            )
            .unwrap();
        let thermal = thermal(vec![cell(1, 2.0, 2.0, 0.0)], &SENSIBLE, vec![]);
        let heating = LuminescentOpticalHeating::try_new(
            &medium,
            source(2.0, 0.75),
            deposition,
            &[(0, CellId(1)), (1, CellId(1))],
            &thermal,
            2,
            1e-12,
            1e-9,
            || Ok(()),
        )
        .unwrap();
        let (next, accepted, audit) = thermal
            .advance_observed(&heating, 1.0, config(1.0), || Ok(()))
            .unwrap();
        close(
            accepted.stored_energy_j(0).unwrap(),
            2.0 * accepted.stored_energy_j(1).unwrap(),
        );
        assert_eq!(accepted.cell_emission_w().len(), 2);
        close(
            accepted.cell_emission_w()[0].1[0],
            2.0 * accepted.cell_emission_w()[1].1[0],
        );
        close(accepted.pump_energy().absorbed_j[0], 7.5);
        close(
            next.temperature_k(CellId(1)).unwrap(),
            300.0 + audit.external_heat_j / 8.0,
        );
        assert_partition(&accepted);
    });
}

#[test]
fn g0_authored_level_and_pulse_domains_refuse_inconsistent_sources() {
    let base = *source(1.0, 0.8).fluorescence();
    for level in [499.0, 601.0, f64::NAN, f64::INFINITY] {
        assert!(LuminescenceSource::try_new(base, level, 1.0).is_err());
    }
    for lifetime in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(LuminescenceSource::try_new(base, 500.0, lifetime).is_err());
    }
    with_cx(|cx| {
        let (_, heating) = setup(cx, source(1.0, 0.8));
        let ramp = SpectralPowerSchedule::try_new(
            &[SpectralPowerSegment {
                start_s: 0.0,
                end_s: 1.0,
                start_power_w: [0.0],
                end_power_w: [10.0],
            }],
            1,
        )
        .unwrap();
        assert!(heating.clone().with_pump_schedule(ramp).is_err());
        let over = SpectralPowerSchedule::try_new(
            &[SpectralPowerSegment {
                start_s: 0.0,
                end_s: 1.0,
                start_power_w: [11.0],
                end_power_w: [11.0],
            }],
            1,
        )
        .unwrap();
        assert!(heating.with_pump_schedule(over).is_err());
    });
}

#[test]
fn g4_late_cancellation_and_temperature_refusal_leave_all_stores_unchanged() {
    with_cx(|cx| {
        let (thermal, heating) = setup(cx, source(1.0, 0.8));
        let mut calls = 0;
        assert!(
            thermal
                .advance_observed(&heating, 1.0, config(0.01), || {
                    calls += 1;
                    if calls == 150 {
                        Err(DepositionError::Cancelled)
                    } else {
                        Ok(())
                    }
                })
                .is_err()
        );
        assert_eq!(heating.energy().stored_excitation_j, 0.0);
        assert_eq!(heating.energy().escaped_j, [0.0; 2]);
        assert_eq!(heating.time_s(), 0.0);
        assert_eq!(thermal.temperature_k(CellId(1)).unwrap(), 300.0);
        let narrow = LuminescenceSource::try_new(
            FluorescenceSource::try_new(
                ContentHash([113; 32]),
                (380.0, 500.0),
                (299.0, 300.05),
                0.8,
                [600.0, 800.0],
                [1.0, 3.0],
            )
            .unwrap(),
            500.0,
            1.0,
        )
        .unwrap();
        let (thermal, heating) = setup(cx, narrow);
        assert!(
            thermal
                .advance_observed(&heating, 1.0, config(0.01), || Ok(()))
                .is_err()
        );
        assert_eq!(thermal.time_s(), 0.0);
        assert_eq!(heating.energy().stored_excitation_j, 0.0);
        assert_eq!(heating.pump_energy().thermal_input_j, 0.0);
    });
}

#[test]
fn g0_interval_dependent_heat_cannot_be_reported_as_a_complete_instantaneous_balance() {
    with_cx(|cx| {
        let (thermal, heating) = setup(cx, source(1.0, 0.8));
        // A native interval is needed for decay-integrated heat, even when the
        // initial population is zero. The observer advertises that requirement.
        assert!(thermal.power_balance_observed(&heating, || Ok(())).is_err());
        assert_eq!(thermal.time_s(), heating.time_s());
    });
}

#[test]
fn g4_implicit_retry_replaces_tentative_energy_without_advancing_the_store() {
    with_cx(|cx| {
        let (thermal, heating) = setup(cx, source(1.0, 0.8));
        let mut candidate = heating.clone();
        let rejected = candidate
            .step_sources(&thermal, 0.0, 1.0, 1.0, &mut || Ok(()))
            .unwrap();
        let retried = candidate
            .step_sources(&thermal, 0.0, 0.25, 0.25, &mut || Ok(()))
            .unwrap();
        let mut direct = heating.clone();
        let short = direct
            .step_sources(&thermal, 0.0, 0.25, 0.25, &mut || Ok(()))
            .unwrap();
        assert_eq!(retried, short);
        assert_ne!(rejected, retried);
        assert_eq!(candidate.energy().stored_excitation_j, 0.0);
        let (next, accepted, _) = thermal
            .advance_observed(&heating, 0.25, config(0.25), || Ok(()))
            .unwrap();
        candidate
            .step(
                0.0,
                0.25,
                &std::collections::BTreeMap::from([(CellId(1), 300.0)]),
                &next,
                &mut || Ok(()),
            )
            .unwrap();
        assert_eq!(candidate.energy(), accepted.energy());
        assert_eq!(candidate.pump_energy(), accepted.pump_energy());
        assert_eq!(candidate.cell_emission_w(), accepted.cell_emission_w());
    });
}

#[test]
fn g0_unrepresentable_pump_scaling_cannot_silently_lose_transmitted_energy() {
    with_cx(|cx| {
        let density = [1.0];
        let medium = grid(cx, &density, spectrum(0.0));
        let deposition = medium
            .deposit_beam_power(
                cx,
                Point3::new(0.0, 0.5, 0.5),
                Vec3::new(1.0, 0.0, 0.0),
                1.0,
                &[400.0],
                &[1e308],
                1,
            )
            .unwrap();
        let thermal = thermal(vec![cell(1, 1.0, 1.0, 0.0)], &SENSIBLE, vec![]);
        let schedule = SpectralPowerSchedule::try_new(
            &[SpectralPowerSegment {
                start_s: 0.0,
                end_s: 1.0,
                start_power_w: [1e-20],
                end_power_w: [1e-20],
            }],
            1,
        )
        .unwrap();
        let heating = LuminescentOpticalHeating::try_new(
            &medium,
            source(1.0, 0.8),
            deposition,
            &[],
            &thermal,
            1,
            1e-12,
            1e-9,
            || Ok(()),
        )
        .unwrap()
        .with_pump_schedule(schedule)
        .unwrap();
        assert!(
            thermal
                .advance_observed(&heating, 1.0, config(1.0), || Ok(()))
                .is_err()
        );
        assert_eq!(heating.pump_energy().incident_j, [0.0]);
    });
}
