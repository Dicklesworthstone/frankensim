//! Coupled phosphorescent photon return, receiver heat and dark persistence.

#[path = "phosphorescence_recycling/annihilation.rs"]
mod annihilation;

use super::*;
use fs_blake3::ContentHash;
use fs_render::volumes::heating::exchange::RadiativeThermalCell;
use fs_render::volumes::heating::fluorescence::FluorescenceSource;
use fs_render::volumes::heating::illumination::schedule::{
    SpectralPowerSchedule, SpectralPowerSegment,
};
use fs_render::volumes::heating::luminescence::{
    AbsorbingLuminescenceEnvironment, LuminescenceAbsorptionTolerance, LuminescenceSource,
    LuminescentOpticalHeating, PhosphorescenceReceiver, SingletTripletLimits, SingletTripletSource,
};
use fs_render::volumes::transport::scattering::luminescence::{
    FractionEstimate, LuminescenceAbsorption, LuminescenceTransferOperator,
    LuminescenceTransferRow, LuminescenceTransportConfig,
};
use fs_thermal::enthalpy::DepositionObserver;
use fs_thermal::luminescence::ArrheniusQuenching;

type Heating = LuminescentOpticalHeating<1, 2>;

fn fraction(mean: [f64; 2]) -> FractionEstimate<2> {
    FractionEstimate {
        mean,
        standard_error: [0.0; 2],
    }
}

fn source(
    wavelengths: [f64; 2],
    level: f64,
    tau: f64,
    phi: f64,
    weights: [f64; 2],
    id: u8,
) -> LuminescenceSource<2> {
    LuminescenceSource::try_new(
        FluorescenceSource::try_new(
            ContentHash([id; 32]),
            (380.0, 500.0),
            (290.0, 800.0),
            phi,
            wavelengths,
            weights,
        )
        .unwrap(),
        level,
        tau,
    )
    .unwrap()
}

/// Manufactured directed fractions isolate kinetics. The transport owner has
/// independent ray tests; these rows are not measurements of this geometry.
fn setup(
    cx: &Cx<'_>,
    wavelengths: [f64; 2],
    rows: &[Vec<[f64; 2]>],
    tails: &[[f64; 2]],
    crossing: f64,
    triplet_yield: f64,
) -> (DepositionSimulation<'static>, Heating) {
    let count = rows[0].len();
    let density = vec![1.0; count];
    let medium = ScatteringGrid::try_new(
        cx,
        grid(
            cx,
            &density,
            MassAbsorptionSpectrum::source_table(
                &[(380.0, 1e4), (900.0, 1e4)],
                ContentHash([229; 32]),
            )
            .unwrap(),
        ),
        MassScatteringSpectrum::source_table(&[(380.0, 0.0), (900.0, 0.0)], ContentHash([230; 32]))
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
        wavelengths_nm: wavelengths,
        grid_cell_count: count,
        material_cells: (0..count).collect(),
        rows: rows
            .iter()
            .enumerate()
            .map(|(emitter, row)| {
                let absorbed = [
                    row.iter().map(|v| v[0]).sum::<f64>(),
                    row.iter().map(|v| v[1]).sum::<f64>(),
                ];
                LuminescenceTransferRow {
                    emitter_cell: emitter,
                    absorption: row
                        .iter()
                        .enumerate()
                        .filter(|(_, value)| value.iter().any(|x| *x > 0.0))
                        .map(|(receiver, &value)| LuminescenceAbsorption {
                            cell_index: receiver,
                            fraction: fraction(value),
                        })
                        .collect(),
                    absorbed: fraction(absorbed),
                    escape: fraction([
                        1.0 - absorbed[0] - tails[emitter][0],
                        1.0 - absorbed[1] - tails[emitter][1],
                    ]),
                    unresolved: fraction(tails[emitter]),
                    precision_limited_histories: [0; 2],
                }
            })
            .collect(),
        traced_segments: 4 * rows.len() as u64,
        config: LuminescenceTransportConfig {
            histories_per_emitter_lane: 2,
            max_collisions_per_history: 2,
            max_total_segments: 8 * rows.len() as u64,
            max_material_cells: count,
            max_emitter_cells: count,
            max_absorption_links: count * count,
            seed: 73,
            source_id: 29,
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
    let deposition = medium
        .absorption_grid()
        .deposit_beam_power(
            cx,
            Point3::new(0.0, 0.5, 0.5),
            Vec3::new(1.0, 0.0, 0.0),
            1.0,
            &[400.0],
            &[10.0],
            count,
        )
        .unwrap();
    let singlet = source(wavelengths, 500.0, 1.0, 1.0, [1.0, 0.0], 227);
    let triplet = source(wavelengths, 600.0, 5.0, triplet_yield, [0.0, 1.0], 228);
    let heating = Heating::try_new_with_emission_transport(
        &medium,
        singlet,
        deposition,
        &[(0, CellId(1))],
        &thermal,
        count,
        1e-12,
        1e-8,
        environment,
        || Ok(()),
    )
    .unwrap()
    .with_singlet_triplet(
        SingletTripletSource::try_new(
            triplet,
            crossing,
            ArrheniusQuenching::try_new(0.0, 0.0).unwrap(),
            SingletTripletLimits::default(),
        )
        .unwrap(),
        &thermal,
        || Ok(()),
    )
    .unwrap();
    (thermal, heating)
}

fn receiver(grid_cell: usize, lane: usize, triplet_fraction: f64) -> PhosphorescenceReceiver<2> {
    let mut fractions = [[0.0; 2]; 2];
    fractions[lane][1] = triplet_fraction;
    PhosphorescenceReceiver {
        grid_cell,
        excitation_fractions: fractions,
    }
}

fn pulse(heating: Heating, end_s: f64) -> Heating {
    heating
        .with_pump_schedule(
            SpectralPowerSchedule::try_new(
                &[SpectralPowerSegment {
                    start_s: 0.0,
                    end_s,
                    start_power_w: [10.0],
                    end_power_w: [10.0],
                }],
                1,
            )
            .unwrap(),
        )
        .unwrap()
}

fn partition(heating: &Heating) {
    let e = heating.energy();
    let heat = e.prompt_heat_j
        + e.nonradiative_heat_j
        + e.annihilation_heat_j
        + e.relaxation_heat_j
        + e.crossing_heat_j
        - e.reverse_crossing_heat_absorbed_j
        + e.secondary_absorption_heat_j;
    close(heat, heating.pump_energy().thermal_input_j);
    close(
        heating.pump_energy().absorbed_j[0],
        e.stored_excitation_j
            + heat
            + e.escaped_j.iter().sum::<f64>()
            + e.unresolved_radiation_j.iter().sum::<f64>(),
    );
    // Internal absorption partitions into stored excitation and receiver heat.
    close(
        e.emitted_j.iter().sum::<f64>(),
        e.recycled_excitation_j
            + e.secondary_absorption_heat_j
            + e.escaped_j.iter().sum::<f64>()
            + e.unresolved_radiation_j.iter().sum::<f64>(),
    );
    assert!(e.absolute_roundoff_j < 1e-8);
}

#[test]
fn g1_triplet_self_recycling_changes_afterglow_decay_inside_each_interval() {
    with_cx(|cx| {
        let (thermal, owner) = setup(
            cx,
            [500.0, 600.0],
            &[vec![[0.0, 0.75]]],
            &[[0.0; 2]],
            3.0,
            0.5,
        );
        let owner = owner
            .with_singlet_triplet_recycling(&[receiver(0, 1, 1.0)], &thermal, || Ok(()))
            .unwrap();
        let owner = pulse(owner, 1.0);
        let (charged_thermal, charged, _) = thermal
            .advance_observed(&owner, 1.0, config(1.0), || Ok(()))
            .unwrap();
        let beta: f64 = 0.2 - 0.75 * 0.1;
        let ns = 2.0 * (1.0 - (-4.0_f64).exp());
        let nt = 5.0
            * ((1.0 - (-beta).exp()) / beta - ((-beta).exp() - (-4.0_f64).exp()) / (4.0 - beta));
        let stores = charged.state_stored_energy_j(0).unwrap();
        close(stores[0], ns);
        close(stores[1], nt);
        let (_, dark, _) = charged_thermal
            .advance_observed(&charged, 5.0, config(4.0), || Ok(()))
            .unwrap();
        let fast = (-16.0_f64).exp();
        let slow = (-4.0 * beta).exp();
        let expected_t = nt * slow + 2.5 * ns * (slow - fast) / (4.0 - beta);
        close(dark.state_stored_energy_j(0).unwrap()[1], expected_t);
        close(dark.cell_emission_w()[0].1[1], 0.1 * expected_t);
        close(
            dark.energy().recycled_excitation_j,
            0.75 * dark.energy().emitted_j[1],
        );
        close(dark.energy().secondary_absorption_heat_j, 0.0);
        partition(&charged);
        partition(&dark);
    });
}

#[test]
fn g1_singlet_photons_excite_unpumped_triplets_and_heat_passive_receivers() {
    with_cx(|cx| {
        let rows = [
            vec![[0.0, 0.0], [0.5, 0.0], [0.25, 0.0]],
            vec![[0.0; 2]; 3],
            vec![[0.0; 2]; 3],
        ];
        let (thermal, owner) = setup(cx, [550.0, 700.0], &rows, &[[0.0; 2]; 3], 0.0, 0.5);
        let owner = owner
            .with_singlet_triplet_recycling(&[receiver(1, 0, 0.8)], &thermal, || Ok(()))
            .unwrap();
        let (heated, accepted, _) = thermal
            .advance_observed(&owner, 1.0, config(1.0), || Ok(()))
            .unwrap();
        let singlet = 8.0 * (1.0 - (-1.0_f64).exp());
        let release = 8.0 - singlet;
        let triplet = (8.0 / 3.0)
            * ((1.0 - (-0.2_f64).exp()) / 0.2 - ((-0.2_f64).exp() - (-1.0_f64).exp()) / 0.8);
        close(accepted.state_stored_energy_j(0).unwrap()[0], singlet);
        close(accepted.state_stored_energy_j(0).unwrap()[1], 0.0);
        close(accepted.state_stored_energy_j(1).unwrap()[0], 0.0);
        close(accepted.state_stored_energy_j(1).unwrap()[1], triplet);
        assert!(accepted.state_stored_energy_j(2).is_none());
        let emitted_s = release * 500.0 / 550.0;
        let incoming = release / 3.0;
        close(accepted.energy().recycled_excitation_j, incoming);
        close(
            accepted.energy().secondary_absorption_heat_j,
            0.75 * emitted_s - incoming,
        );
        close(
            heated.temperature_k(CellId(3)).unwrap(),
            300.0 + 0.25 * emitted_s / 4.0,
        );
        let triplet_decay = incoming - triplet;
        let receiver_heat = 0.5 * emitted_s - incoming
            + 0.5 * triplet_decay
            + 0.5 * triplet_decay * (1.0 - 600.0 / 700.0);
        close(
            heated.temperature_k(CellId(2)).unwrap(),
            300.0 + receiver_heat / 4.0,
        );
        partition(&accepted);
    });
}

#[test]
fn g1_closed_triplet_photon_cycles_keep_external_energy_and_stable_afterglow() {
    with_cx(|cx| {
        let rows = [vec![[0.0; 2], [0.0, 1.0]], vec![[0.0, 1.0], [0.0; 2]]];
        let (thermal, owner) = setup(cx, [500.0, 600.0], &rows, &[[0.0; 2]; 2], 3.0, 1.0);
        let owner = owner
            .with_singlet_triplet_recycling(
                &[receiver(0, 1, 1.0), receiver(1, 1, 1.0)],
                &thermal,
                || Ok(()),
            )
            .unwrap();
        let owner = pulse(owner, 1.0);
        let (_, accepted, _) = thermal
            .advance_observed(&owner, 20.0, config(20.0), || Ok(()))
            .unwrap();
        close(accepted.energy().stored_excitation_j, 5.0);
        close(accepted.energy().escaped_j[1], 0.0);
        close(accepted.energy().secondary_absorption_heat_j, 0.0);
        assert!(accepted.energy().emitted_j[1] > 10.0);
        assert!(accepted.state_stored_energy_j(1).unwrap()[1] > 0.0);
        close(
            accepted.energy().emitted_j[1],
            accepted.energy().recycled_excitation_j,
        );
        partition(&accepted);
    });
}

#[test]
fn g0_recycling_refuses_uphill_photons_and_excess_excitation_probability() {
    with_cx(|cx| {
        let (thermal, owner) = setup(
            cx,
            [500.0, 600.0],
            &[vec![[0.5, 0.5]]],
            &[[0.0; 2]],
            3.0,
            0.5,
        );
        let uphill = PhosphorescenceReceiver {
            grid_cell: 0,
            excitation_fractions: [[0.0, 0.0], [1.0, 0.0]],
        };
        assert!(
            owner
                .clone()
                .with_singlet_triplet_recycling(&[uphill], &thermal, || Ok(()),)
                .is_err()
        );
        let excessive = PhosphorescenceReceiver {
            grid_cell: 0,
            excitation_fractions: [[0.8, 0.8], [0.0, 0.0]],
        };
        assert!(
            owner
                .with_singlet_triplet_recycling(&[excessive], &thermal, || Ok(()),)
                .is_err()
        );
    });
}

#[test]
fn g4_recycling_retry_and_cancellation_do_not_publish_tentative_populations() {
    with_cx(|cx| {
        let (thermal, owner) = setup(
            cx,
            [500.0, 600.0],
            &[vec![[0.0, 0.75]]],
            &[[0.0; 2]],
            3.0,
            0.5,
        );
        let owner = owner
            .with_singlet_triplet_recycling(&[receiver(0, 1, 1.0)], &thermal, || Ok(()))
            .unwrap();
        let mut candidate = owner.clone();
        candidate
            .step_sources(&thermal, 0.0, 1.0, 1.0, &mut || Ok(()))
            .unwrap();
        let retry = candidate
            .step_sources(&thermal, 0.0, 0.5, 0.5, &mut || Ok(()))
            .unwrap();
        let fresh = owner
            .clone()
            .step_sources(&thermal, 0.0, 0.5, 0.5, &mut || Ok(()))
            .unwrap();
        assert_eq!(retry, fresh);
        assert_eq!(
            candidate.state_stored_energy_j(0),
            owner.state_stored_energy_j(0)
        );
        assert_eq!(candidate.energy(), owner.energy());
        let mut cancelled = owner.clone();
        let mut calls = 0;
        let result = cancelled.step_sources(&thermal, 0.0, 1.0, 1.0, &mut || {
            calls += 1;
            if calls == 12 {
                Err(DepositionError::Cancelled)
            } else {
                Ok(())
            }
        });
        assert!(result.is_err());
        assert_eq!(
            cancelled.state_stored_energy_j(0),
            owner.state_stored_energy_j(0)
        );
        assert_eq!(cancelled.energy(), owner.energy());
    });
}
