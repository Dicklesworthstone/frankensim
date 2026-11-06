use super::*;
use fs_blake3::ContentHash;
use fs_render::volumes::heating::exchange::RadiativeThermalCell;
use fs_render::volumes::heating::fluorescence::FluorescenceSource;
use fs_render::volumes::heating::illumination::schedule::{
    SpectralPowerSchedule, SpectralPowerSegment,
};
use fs_render::volumes::heating::luminescence::{
    AbsorbingLuminescenceEnvironment, LuminescenceAbsorptionTolerance, LuminescenceSource,
    LuminescentOpticalHeating,
};
use fs_render::volumes::transport::scattering::luminescence::{
    FractionEstimate, LuminescenceAbsorption, LuminescenceTransferOperator,
    LuminescenceTransferRow, LuminescenceTransportConfig,
};
use fs_thermal::enthalpy::DepositionObserver;
use fs_thermal::luminescence::recycling::RecyclingLimits;

const PUMP_W: f64 = 12.0;
const STORE_PUMP_W: f64 = 8.0; // 12 W at 400 nm relaxed to the 600 nm state.
const PROMPT_HEAT_W: f64 = 4.0;

fn source_with_domain(
    tau: f64,
    phi: f64,
    excitation_max_nm: f64,
    state_nm: f64,
) -> LuminescenceSource<1> {
    LuminescenceSource::try_new(
        FluorescenceSource::try_new(
            ContentHash([181; 32]),
            (380.0, excitation_max_nm),
            (290.0, 800.0),
            phi,
            [600.0],
            [1.0],
        )
        .unwrap(),
        state_nm,
        tau,
    )
    .unwrap()
}

fn source(tau: f64, phi: f64) -> LuminescenceSource<1> {
    source_with_domain(tau, phi, 600.0, 600.0)
}

fn fraction(value: f64) -> FractionEstimate<1> {
    FractionEstimate {
        mean: [value],
        standard_error: [0.0],
    }
}

/// Independent directed-probability fixtures isolate the kinetics/thermal
/// coupling. Grid-derived transport probabilities have their own ray tests.
/// Rows are emitters; columns are original receiving optical cells. Supplying
/// fewer rows than columns deliberately leaves unpumped emitters unprepared.
fn passive(
    cx: &Cx<'_>,
    source: LuminescenceSource<1>,
    absorption: &[&[f64]],
    unresolved: &[f64],
) -> (
    DepositionSimulation<'static>,
    LuminescentOpticalHeating<1, 1>,
) {
    let count = absorption[0].len();
    let density = vec![1.0; count];
    let medium = ScatteringGrid::try_new(
        cx,
        grid(
            cx,
            &density,
            MassAbsorptionSpectrum::source_table(
                &[(380.0, 1e4), (900.0, 1e4)],
                ContentHash([182; 32]),
            )
            .unwrap(),
        ),
        MassScatteringSpectrum::source_table(&[(380.0, 0.0), (900.0, 0.0)], ContentHash([183; 32]))
            .unwrap(),
        PhaseFunction::Isotropic,
    )
    .unwrap();
    let thermal = thermal(
        (0..count)
            .map(|i| cell(i as u64 + 1, 1.0, 1.0, 0.0))
            .collect(),
        &SENSIBLE,
        vec![],
    );
    let operator = LuminescenceTransferOperator {
        medium_identity: medium.content_identity(),
        absorption_identity: medium.absorption_grid().content_identity(),
        wavelength_bounds_nm: medium.wavelength_bounds_nm(),
        wavelengths_nm: [600.0],
        grid_cell_count: count,
        material_cells: (0..count).collect(),
        rows: absorption
            .iter()
            .enumerate()
            .map(|(emitter, row)| {
                let absorbed = row.iter().sum::<f64>();
                LuminescenceTransferRow {
                    emitter_cell: emitter,
                    absorption: row
                        .iter()
                        .enumerate()
                        .filter(|(_, value)| **value > 0.0)
                        .map(|(receiver, &value)| LuminescenceAbsorption {
                            cell_index: receiver,
                            fraction: fraction(value),
                        })
                        .collect(),
                    absorbed: fraction(absorbed),
                    escape: fraction(1.0 - absorbed - unresolved[emitter]),
                    unresolved: fraction(unresolved[emitter]),
                    precision_limited_histories: [0],
                }
            })
            .collect(),
        traced_segments: 2 * absorption.len() as u64,
        config: LuminescenceTransportConfig {
            histories_per_emitter_lane: 2,
            max_collisions_per_history: 2,
            max_total_segments: 4 * absorption.len() as u64,
            max_material_cells: count,
            max_emitter_cells: count,
            max_absorption_links: count * count,
            seed: 37,
            source_id: 19,
        },
    };
    let bindings: Vec<_> = (0..count)
        .map(|i| RadiativeThermalCell {
            grid_cell: i,
            thermal_cell: CellId(i as u64 + 1),
            minimum_temperature_k: 290.0,
            maximum_temperature_k: 800.0,
        })
        .collect();
    let environment = AbsorbingLuminescenceEnvironment::try_new(
        &medium,
        operator,
        &bindings,
        &thermal,
        count,
        LuminescenceAbsorptionTolerance {
            fraction_roundoff: 1e-13,
            max_fraction_standard_error: 0.0,
            max_unresolved_energy_j: 100.0,
        },
        || Ok(()),
    )
    .unwrap();
    // The pump stops at the first cell's far face. All other original cells
    // begin completely unpumped, even though they can absorb emitted photons.
    let deposition = medium
        .absorption_grid()
        .deposit_beam_power(
            cx,
            Point3::new(0.0, 0.5, 0.5),
            Vec3::new(1.0, 0.0, 0.0),
            1.0,
            &[400.0],
            &[PUMP_W],
            count,
        )
        .unwrap();
    let owner = LuminescentOpticalHeating::try_new_with_emission_transport(
        &medium,
        source,
        deposition,
        &[(0, CellId(1))],
        &thermal,
        count,
        1e-12,
        1e-9,
        environment,
        || Ok(()),
    )
    .unwrap();
    (thermal, owner)
}

fn pulse(end_s: f64) -> SpectralPowerSchedule<1> {
    SpectralPowerSchedule::try_new(
        &[SpectralPowerSegment {
            start_s: 0.0,
            end_s,
            start_power_w: [PUMP_W],
            end_power_w: [PUMP_W],
        }],
        1,
    )
    .unwrap()
}

fn emitted(powers: &[(usize, [f64; 1])], cell: usize) -> f64 {
    powers.iter().find(|(index, _)| *index == cell).unwrap().1[0]
}

fn partition(owner: &LuminescentOpticalHeating<1, 1>) {
    let e = owner.energy();
    assert_eq!(e.secondary_absorption_heat_j, 0.0);
    assert_eq!(e.relaxation_heat_j, 0.0);
    let heat = e.prompt_heat_j + e.nonradiative_heat_j;
    close(heat, owner.pump_energy().thermal_input_j);
    close(
        owner.pump_energy().absorbed_j[0],
        heat + e.stored_excitation_j + e.escaped_j[0] + e.unresolved_radiation_j[0],
    );
    // Recycling is an internal transfer; adding it to the final pump balance
    // would count the same excitation repeatedly as it emits and returns.
    close(
        e.emitted_j[0],
        e.recycled_excitation_j + e.escaped_j[0] + e.unresolved_radiation_j[0],
    );
    assert!(e.residual_j.abs() < 1e-9);
}

/// Independent scalar solution U'=P-kU with U(0)=0 and its residence integral.
fn pumped_mode(pump: f64, rate: f64, duration: f64) -> (f64, f64) {
    let response = -(-rate * duration).exp_m1() / rate;
    (pump * response, pump * (duration - response) / rate)
}

#[test]
fn g1_self_recycling_uses_effective_lifetime_and_keeps_all_energy_channels_distinct() {
    with_cx(|cx| {
        let (tau, phi, returned, tail, time) = (2.0, 0.8, 0.75, 0.125, 3.0);
        let (thermal, passive) = passive(cx, source(tau, phi), &[&[returned]], &[tail]);
        let owner = passive
            .with_resonant_recycling(&thermal, RecyclingLimits::default(), || Ok(()))
            .unwrap();
        let (next, accepted, audit) = thermal
            .advance_observed(&owner, time, config(time), || Ok(()))
            .unwrap();
        let rate = (1.0 - phi * returned) / tau;
        let (stored, residence) = pumped_mode(STORE_PUMP_W, rate, time);
        let released = residence / tau;
        let produced = phi * released;
        close(accepted.stored_energy_j(0).unwrap(), stored);
        close(accepted.energy().stored_excitation_j, stored);
        close(accepted.energy().emitted_j[0], produced);
        close(accepted.energy().recycled_excitation_j, returned * produced);
        close(
            accepted.energy().escaped_j[0],
            (1.0 - returned - tail) * produced,
        );
        close(accepted.energy().unresolved_radiation_j[0], tail * produced);
        close(
            accepted.energy().nonradiative_heat_j,
            (1.0 - phi) * released,
        );
        close(accepted.energy().prompt_heat_j, PROMPT_HEAT_W * time);
        close(emitted(accepted.cell_emission_w(), 0), phi * stored / tau);
        close(
            emitted(accepted.last_mean_cell_emission_w(), 0),
            produced / time,
        );
        assert!(emitted(accepted.cell_emission_w(), 0) > produced / time);
        close(
            audit.external_heat_j,
            PROMPT_HEAT_W * time + (1.0 - phi) * released,
        );
        close(
            next.temperature_k(CellId(1)).unwrap(),
            300.0 + audit.external_heat_j / 4.0,
        );
        assert_eq!(audit.substeps, 1);
        partition(&accepted);
        assert_eq!(owner.energy().stored_excitation_j, 0.0);
        assert_eq!(thermal.time_s(), 0.0);
    });
}

#[test]
fn g1_cross_recycling_excites_an_unpumped_cell_and_matches_pulse_dark_eigenmodes() {
    with_cx(|cx| {
        let (tau, phi, returned, pulse_end, dark) = (1.2, 0.75, 0.6, 1.25, 3.0);
        let (thermal, passive) = passive(
            cx,
            source(tau, phi),
            &[&[0.0, returned], &[returned, 0.0]],
            &[0.0, 0.0],
        );
        assert_eq!(passive.stored_energy_j(1), None);
        let owner = passive
            .with_resonant_recycling(&thermal, RecyclingLimits::default(), || Ok(()))
            .unwrap()
            .with_pump_schedule(pulse(pulse_end))
            .unwrap();
        assert_eq!(owner.stored_energy_j(1), Some(0.0));
        let (middle, first, _) = thermal
            .advance_observed(&owner, pulse_end, config(pulse_end), || Ok(()))
            .unwrap();
        // Sum and difference diagonalize this symmetric two-cell problem.
        let slow = (1.0 - phi * returned) / tau;
        let fast = (1.0 + phi * returned) / tau;
        let (sum, sum_integral) = pumped_mode(STORE_PUMP_W, slow, pulse_end);
        let (difference, difference_integral) = pumped_mode(STORE_PUMP_W, fast, pulse_end);
        let pulse_store = [(sum + difference) / 2.0, (sum - difference) / 2.0];
        let pulse_integral = [
            (sum_integral + difference_integral) / 2.0,
            (sum_integral - difference_integral) / 2.0,
        ];
        for i in 0..2 {
            close(first.stored_energy_j(i).unwrap(), pulse_store[i]);
            close(
                emitted(first.cell_emission_w(), i),
                phi * pulse_store[i] / tau,
            );
            close(
                emitted(first.last_mean_cell_emission_w(), i),
                phi * pulse_integral[i] / (tau * pulse_end),
            );
            let prompt = if i == 0 {
                PROMPT_HEAT_W * pulse_end
            } else {
                0.0
            };
            close(
                middle.temperature_k(CellId(i as u64 + 1)).unwrap(),
                300.0 + (prompt + (1.0 - phi) * pulse_integral[i] / tau) / 4.0,
            );
        }
        assert!(pulse_store[1] > 0.0);
        assert!(emitted(first.cell_emission_w(), 1) > 0.0);

        let end = pulse_end + dark;
        // The target is the absolute native end time, not the dark duration.
        let (split, after, _) = middle
            .advance_observed(&first, end, config(dark), || Ok(()))
            .unwrap();
        let sum_end = sum * (-slow * dark).exp();
        let difference_end = difference * (-fast * dark).exp();
        let dark_sum_integral = sum * -(-slow * dark).exp_m1() / slow;
        let dark_difference_integral = difference * -(-fast * dark).exp_m1() / fast;
        let dark_integral = [
            (dark_sum_integral + dark_difference_integral) / 2.0,
            (dark_sum_integral - dark_difference_integral) / 2.0,
        ];
        let final_store = [
            (sum_end + difference_end) / 2.0,
            (sum_end - difference_end) / 2.0,
        ];
        for i in 0..2 {
            close(after.stored_energy_j(i).unwrap(), final_store[i]);
            close(
                emitted(after.cell_emission_w(), i),
                phi * final_store[i] / tau,
            );
            close(
                emitted(after.last_mean_cell_emission_w(), i),
                phi * dark_integral[i] / (tau * dark),
            );
            let prompt = if i == 0 {
                PROMPT_HEAT_W * pulse_end
            } else {
                0.0
            };
            close(
                split.temperature_k(CellId(i as u64 + 1)).unwrap(),
                300.0 + (prompt + (1.0 - phi) * (pulse_integral[i] + dark_integral[i]) / tau) / 4.0,
            );
        }
        let produced = phi * (sum_integral + dark_sum_integral) / tau;
        close(after.energy().emitted_j[0], produced);
        close(after.energy().recycled_excitation_j, returned * produced);
        close(after.energy().escaped_j[0], (1.0 - returned) * produced);
        close(after.pump_energy().incident_j[0], PUMP_W * pulse_end);
        let (whole, one, audit) = thermal
            .advance_observed(&owner, end, config(end), || Ok(()))
            .unwrap();
        assert_eq!(audit.substeps, 2);
        assert_eq!(one.time_s(), end);
        close(
            one.energy().stored_excitation_j,
            after.energy().stored_excitation_j,
        );
        close(one.energy().emitted_j[0], after.energy().emitted_j[0]);
        close(
            one.energy().recycled_excitation_j,
            after.energy().recycled_excitation_j,
        );
        close(whole.energy_j(), split.energy_j());
        partition(&first);
        partition(&after);
        partition(&one);
    });
}

#[test]
fn g1_asymmetric_recycling_preserves_source_receiver_orientation_in_each_cell() {
    with_cx(|cx| {
        let (tau, phi, time): (f64, f64, f64) = (1.0, 0.8, 2.0);
        let (thermal, passive) = passive(
            cx,
            source(tau, phi),
            &[&[0.0, 0.5], &[0.25, 0.0]],
            &[0.0, 0.0],
        );
        let owner = passive
            .with_resonant_recycling(&thermal, RecyclingLimits::default(), || Ok(()))
            .unwrap();
        let (next, accepted, audit) = thermal
            .advance_observed(&owner, time, config(time), || Ok(()))
            .unwrap();

        // The original directed rows give U0'=P-kU0+bU1 and U1'=aU0-kU1.
        // The eigenvalues depend on a*b, but the receiving cell's amplitude
        // depends on a, so transposing the transfer matrix cannot pass this case.
        let k = 1.0 / tau;
        let a = phi * 0.5 / tau;
        let b = phi * 0.25 / tau;
        let g = (a * b).sqrt();
        let (slow, slow_integral) = pumped_mode(STORE_PUMP_W, k - g, time);
        let (fast, fast_integral) = pumped_mode(STORE_PUMP_W, k + g, time);
        let stores = [(slow + fast) / 2.0, a * (slow - fast) / (2.0 * g)];
        let residence = [
            (slow_integral + fast_integral) / 2.0,
            a * (slow_integral - fast_integral) / (2.0 * g),
        ];
        let produced = residence.map(|integral| phi * integral / tau);
        let dark_heat = residence.map(|integral| (1.0 - phi) * integral / tau);
        for i in 0..2 {
            close(accepted.stored_energy_j(i).unwrap(), stores[i]);
            close(
                emitted(accepted.cell_emission_w(), i),
                phi * stores[i] / tau,
            );
            close(
                emitted(accepted.last_mean_cell_emission_w(), i),
                produced[i] / time,
            );
            let prompt = if i == 0 { PROMPT_HEAT_W * time } else { 0.0 };
            close(
                next.temperature_k(CellId(i as u64 + 1)).unwrap(),
                300.0 + (prompt + dark_heat[i]) / 4.0,
            );
        }
        close(accepted.energy().emitted_j[0], produced[0] + produced[1]);
        close(
            accepted.energy().recycled_excitation_j,
            0.5 * produced[0] + 0.25 * produced[1],
        );
        close(
            accepted.energy().escaped_j[0],
            0.5 * produced[0] + 0.75 * produced[1],
        );
        close(
            accepted.energy().nonradiative_heat_j,
            dark_heat[0] + dark_heat[1],
        );
        close(
            audit.external_heat_j,
            PROMPT_HEAT_W * time + dark_heat[0] + dark_heat[1],
        );

        // Reconstruct incoming excitation from each native cell's store and
        // integrated emission, independently of the analytic eigenmode values.
        // U_i = pump_i + incoming_i - emitted_i/phi for this resonant source.
        let actual_produced =
            [0, 1].map(|i| emitted(accepted.last_mean_cell_emission_w(), i) * time);
        let incoming_0 =
            accepted.stored_energy_j(0).unwrap() - STORE_PUMP_W * time + actual_produced[0] / phi;
        let incoming_1 = accepted.stored_energy_j(1).unwrap() + actual_produced[1] / phi;
        close(incoming_0, 0.25 * actual_produced[1]);
        close(incoming_1, 0.5 * actual_produced[0]);
        assert_eq!(audit.substeps, 1);
        partition(&accepted);
    });
}

#[test]
fn g1_closed_unit_return_stores_only_the_pump_while_counting_repeated_emission() {
    with_cx(|cx| {
        let tau = 0.5;
        let (thermal, passive) = passive(cx, source(tau, 1.0), &[&[1.0]], &[0.0]);
        let owner = passive
            .with_resonant_recycling(&thermal, RecyclingLimits::default(), || Ok(()))
            .unwrap()
            .with_pump_schedule(pulse(4.0))
            .unwrap();
        let (next, accepted, audit) = thermal
            .advance_observed(&owner, 6.0, config(6.0), || Ok(()))
            .unwrap();
        let stored = STORE_PUMP_W * 4.0;
        // Linear growth during pumping, then a constant excitation store while
        // each emitted photon returns: no artificial decay or heat per cycle.
        let emitted_j = (0.5 * STORE_PUMP_W * 4.0 * 4.0 + stored * 2.0) / tau;
        close(accepted.energy().stored_excitation_j, stored);
        close(accepted.energy().emitted_j[0], emitted_j);
        close(accepted.energy().recycled_excitation_j, emitted_j);
        assert!(accepted.energy().recycled_excitation_j > accepted.pump_energy().absorbed_j[0]);
        assert_eq!(accepted.energy().escaped_j, [0.0]);
        assert_eq!(accepted.energy().unresolved_radiation_j, [0.0]);
        assert_eq!(accepted.energy().nonradiative_heat_j, 0.0);
        close(accepted.energy().prompt_heat_j, 4.0 * PROMPT_HEAT_W);
        close(emitted(accepted.cell_emission_w(), 0), stored / tau);
        close(
            emitted(accepted.last_mean_cell_emission_w(), 0),
            stored / tau,
        );
        close(audit.external_heat_j, 4.0 * PROMPT_HEAT_W);
        close(next.temperature_k(CellId(1)).unwrap(), 304.0);
        assert_eq!(audit.substeps, 2);
        partition(&accepted);
    });
}

#[test]
fn g0_recycling_requires_complete_states_supported_resonance_and_initial_admission() {
    with_cx(|cx| {
        let (thermal, missing) = passive(cx, source(1.0, 0.8), &[&[0.0, 0.5]], &[0.0]);
        assert!(
            missing
                .with_resonant_recycling(&thermal, RecyclingLimits::default(), || Ok(()))
                .is_err()
        );
        for incompatible in [
            source_with_domain(1.0, 0.8, 500.0, 600.0),
            source_with_domain(1.0, 0.8, 500.0, 550.0),
        ] {
            let (thermal, owner) = passive(cx, incompatible, &[&[0.5]], &[0.0]);
            assert!(
                owner
                    .with_resonant_recycling(&thermal, RecyclingLimits::default(), || Ok(()))
                    .is_err()
            );
        }
        let (thermal, owner) = passive(
            cx,
            source(1.0, 0.8),
            &[&[0.0, 0.5], &[0.5, 0.0]],
            &[0.0, 0.0],
        );
        assert!(
            owner
                .clone()
                .with_resonant_recycling(
                    &thermal,
                    RecyclingLimits {
                        maximum_nodes: 1,
                        ..RecyclingLimits::default()
                    },
                    || Ok(()),
                )
                .is_err()
        );
        assert!(matches!(
            owner
                .clone()
                .with_resonant_recycling(&thermal, RecyclingLimits::default(), || Err(
                    DepositionError::Cancelled
                ),),
            Err(OpticalHeatingError::Thermal(DepositionError::Cancelled))
        ));
        let mut tentative = owner.clone();
        tentative
            .step_sources(&thermal, 0.0, 0.25, 0.25, &mut || Ok(()))
            .unwrap();
        assert!(
            tentative
                .with_resonant_recycling(&thermal, RecyclingLimits::default(), || Ok(()))
                .is_err()
        );
        let (next, advanced, _) = thermal
            .advance_observed(&owner, 0.25, config(0.25), || Ok(()))
            .unwrap();
        assert!(
            advanced
                .with_resonant_recycling(&next, RecyclingLimits::default(), || Ok(()))
                .is_err()
        );
        assert_eq!(owner.time_s(), 0.0);
        assert_eq!(owner.energy().emitted_j, [0.0]);
        assert_eq!(thermal.time_s(), 0.0);
    });
}

#[test]
fn g4_recycling_candidate_retry_clone_and_late_cancellation_preserve_accepted_state() {
    with_cx(|cx| {
        let (thermal, passive) = passive(
            cx,
            source(1.0, 0.8),
            &[&[0.1, 0.5], &[0.5, 0.1]],
            &[0.0, 0.0],
        );
        let owner = passive
            .with_resonant_recycling(&thermal, RecyclingLimits::default(), || Ok(()))
            .unwrap();
        let mut candidate = owner.clone();
        let rejected = candidate
            .step_sources(&thermal, 0.0, 1.0, 1.0, &mut || Ok(()))
            .unwrap();
        let retried = candidate
            .step_sources(&thermal, 0.0, 0.25, 0.25, &mut || Ok(()))
            .unwrap();
        let mut direct = owner.clone();
        let short = direct
            .step_sources(&thermal, 0.0, 0.25, 0.25, &mut || Ok(()))
            .unwrap();
        assert_eq!(retried, short);
        assert_ne!(rejected, retried);
        assert_eq!(candidate.energy(), owner.energy());
        assert_eq!(candidate.stored_energy_j(1), Some(0.0));
        let (next, accepted, _) = thermal
            .advance_observed(&owner, 0.25, config(0.25), || Ok(()))
            .unwrap();
        candidate
            .step(
                0.0,
                0.25,
                &std::collections::BTreeMap::from([(CellId(1), 300.0), (CellId(2), 300.0)]),
                &next,
                &mut || Ok(()),
            )
            .unwrap();
        assert_eq!(candidate.energy(), accepted.energy());
        assert_eq!(candidate.pump_energy(), accepted.pump_energy());
        assert_eq!(candidate.cell_emission_w(), accepted.cell_emission_w());
        assert_eq!(
            candidate.last_mean_cell_emission_w(),
            accepted.last_mean_cell_emission_w()
        );

        // Count a successful native run, then cancel immediately before its
        // final checkpoints. This exercises late rollback without hard-coding
        // how many checkpoints the reusable numerical owner happens to take.
        let snapshot = accepted.clone();
        let before_thermal_energy = next.energy_j();
        let mut total_checks = 0;
        let (complete, baseline, _) = next
            .advance_observed(&accepted, 1.0, config(0.125), || {
                total_checks += 1;
                Ok(())
            })
            .unwrap();
        assert!(total_checks > 2);
        let cancel_at = total_checks - 1;
        let mut checks = 0;
        let refused = next.advance_observed(&accepted, 1.0, config(0.125), || {
            checks += 1;
            if checks == cancel_at {
                Err(DepositionError::Cancelled)
            } else {
                Ok(())
            }
        });
        assert!(matches!(
            refused,
            Err(OpticalHeatingError::Thermal(DepositionError::Cancelled))
        ));
        assert_eq!(checks, cancel_at);
        assert_eq!(accepted.time_s(), 0.25);
        assert_eq!(accepted.energy(), snapshot.energy());
        assert_eq!(accepted.pump_energy(), snapshot.pump_energy());
        assert_eq!(accepted.cell_emission_w(), snapshot.cell_emission_w());
        assert_eq!(next.time_s(), 0.25);
        assert_eq!(next.energy_j(), before_thermal_energy);
        let (replayed_thermal, replayed, _) = next
            .advance_observed(&accepted, 1.0, config(0.125), || Ok(()))
            .unwrap();
        assert_eq!(baseline.energy(), replayed.energy());
        assert_eq!(baseline.cell_emission_w(), replayed.cell_emission_w());
        close(complete.energy_j(), replayed_thermal.energy_j());
        partition(&accepted);
        partition(&replayed);
    });
}
