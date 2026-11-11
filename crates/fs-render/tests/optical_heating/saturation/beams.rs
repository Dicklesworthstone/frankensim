//! Shared inventories under simultaneous incoherent and counterpropagating pumps.

use super::*;
use fs_render::volumes::heating::feedback::IncidentBeam;
use fs_render::volumes::heating::luminescence::saturation::beam::{
    BeamEmitterCrossSections, DepletingPumpBeam, PumpBeamsLimits, PumpColumnLimits,
};

pub(super) fn inventory(
    medium: &AbsorbingGrid<'_>,
    populations: &[[f64; 3]],
    lifetimes: [f64; 2],
    crossing: f64,
) -> (
    DepositionSimulation<'static>,
    SaturableLuminescentHeating<2>,
) {
    let thermal = thermal(
        (0..medium.cell_count())
            .map(|index| cell(index as u64 + 1, 1.0, 1.0, 0.0))
            .collect(),
        &SENSIBLE,
        vec![],
    );
    let pump =
        PrescribedPhotonPump::try_new(ContentHash([174; 32]), [400.0, 450.0], [0.0; 2], [0.0; 2])
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
        source(600.0, lifetimes[1], [0.0, 1.0], 175),
        crossing,
        ArrheniusQuenching::try_new(0.0, 0.0).unwrap(),
        SingletTripletLimits::default(),
    )
    .unwrap();
    let heating = SaturableLuminescentHeating::try_new(
        medium,
        source(500.0, lifetimes[0], [1.0, 0.0], 176),
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

pub(super) fn sections(cells: usize, values: [[f64; 2]; 2]) -> Vec<BeamEmitterCrossSections> {
    (0..cells)
        .map(|index| {
            BeamEmitterCrossSections::authored(index, ContentHash([177; 32]), values).unwrap()
        })
        .collect()
}

pub(super) fn ray(
    cells: usize,
    reverse: bool,
    wavelengths_nm: [f64; 2],
    incident_power_w: [f64; 2],
) -> IncidentBeam<2> {
    IncidentBeam {
        origin: Point3::new(if reverse { cells as f64 + 1.0 } else { -1.0 }, 0.5, 0.5),
        unit_direction: Vec3::new(if reverse { -1.0 } else { 1.0 }, 0.0, 0.0),
        max_distance_m: f64::INFINITY,
        wavelengths_nm,
        incident_power_w,
    }
}

// A pumped pair history takes the smaller of the column's and the pair law's
// substep budgets (fs-thermal CONTRACT). The pair fixtures that share these
// limits author 4096 and 2048 pair substeps, so the column states 4096 rather
// than its default 1024, which would silently cap both.
pub(super) fn limits() -> PumpBeamsLimits {
    PumpBeamsLimits {
        column: PumpColumnLimits {
            relative_tolerance: 2e-5,
            maximum_substeps: 4096,
            maximum_cell_evaluations: 2_000_000,
            ..PumpColumnLimits::default()
        },
        maximum_coupling_iterations: 128,
    }
}

pub(super) fn bind(
    cx: &Cx<'_>,
    medium: &AbsorbingGrid<'_>,
    thermal: &DepositionSimulation<'_>,
    heating: SaturableLuminescentHeating<2>,
    beams: &[DepletingPumpBeam<'_>],
) -> SaturableLuminescentHeating<2> {
    heating
        .with_depleting_beams(
            cx,
            medium,
            beams,
            thermal,
            beams.len(),
            beams.len() * medium.cell_count(),
            limits(),
            || Ok(()),
        )
        .unwrap()
}

fn near(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 2e-4 * actual.abs().max(expected.abs()).max(1e-30),
        "{actual:.17e} != {expected:.17e}"
    );
}

pub(super) fn accounts(heating: &SaturableLuminescentHeating<2>) {
    let report = heating.depleting_pump_report().unwrap();
    let mut active = 0.0;
    let mut passive = 0.0;
    for beam in &report.beams {
        for lane in 0..2 {
            relative(
                beam.incident_j[lane],
                beam.transmitted_j[lane] + beam.absorbed_j[lane] + beam.passive_absorbed_j[lane],
            );
            assert!(beam.residual_j[lane].abs() < 1e-7);
            active += beam.absorbed_j[lane];
            passive += beam.passive_absorbed_j[lane];
        }
    }
    relative(active, heating.absorbed_pump_energy_j().iter().sum());
    relative(passive, heating.energy().passive_pump_absorption_j);
    let energy = heating.energy();
    relative(
        energy.initial_excitation_j + active + passive,
        energy.stored_excitation_j
            + energy.escaped_j.iter().sum::<f64>()
            + energy.unresolved_radiation_j.iter().sum::<f64>()
            + heating.thermal_input_j(),
    );
    assert!(energy.residual_j.abs() < 1e-7);
    assert!(energy.absolute_roundoff_j < 1e-7);
}

#[test]
fn g1_simultaneous_pumps_split_one_field_without_duplicating_emitters_or_source_ids() {
    with_cx(|cx| {
        let density = [1.0];
        let medium = grid(cx, &density, transparent_host());
        let (thermal, heating) = inventory(&medium, &[[COUNT, 0.0, 0.0]], [1e12; 2], 0.0);
        let cross_sections = sections(1, [[1.0 / COUNT, 0.0], [0.0; 2]]);
        let watt = COUNT * HC / 400.0;
        let beams = [3.0, 2.0].map(|weight| DepletingPumpBeam {
            // Two rays may originate from one source. The identity must not
            // deduplicate their separately authored physical powers.
            source_identity: ContentHash([178; 32]),
            beam: ray(1, false, [400.0, 450.0], [weight * watt, 0.0]),
            cross_sections: &cross_sections,
            pump_window_s: None,
        });
        let simultaneous = bind(cx, &medium, &thermal, heating.clone(), &beams);
        let combined = heating
            .with_depleting_beam(
                cx,
                &medium,
                ContentHash([178; 32]),
                ray(1, false, [400.0, 450.0], [5.0 * watt, 0.0]),
                &cross_sections,
                &thermal,
                1,
                limits().column,
                || Ok(()),
            )
            .unwrap();
        let (warm, accepted, audit) = thermal
            .advance_observed(&simultaneous, 0.4, config(0.4), || Ok(()))
            .unwrap();
        let (_, reference, _) = thermal
            .advance_observed(&combined, 0.4, config(0.4), || Ok(()))
            .unwrap();
        for (actual, expected) in accepted
            .populations(0)
            .unwrap()
            .into_iter()
            .zip(reference.populations(0).unwrap())
        {
            near(actual, expected);
        }
        // Independent one-shot fluence solution with initial optical depth1
        // and dimensionless total incident fluence2. Lifetime loss is <1e-12.
        let outgoing = (1.0 + (-1.0_f64).exp() * 2.0_f64.exp_m1()).ln();
        near(
            accepted.populations(0).unwrap()[0],
            COUNT * (outgoing - 1.0),
        );
        relative(accepted.populations(0).unwrap().iter().sum(), COUNT);
        let report = accepted.depleting_pump_report().unwrap();
        assert_eq!(report.beams.len(), 2);
        assert_eq!(
            report.beams[0].source_identity,
            report.beams[1].source_identity
        );
        assert!(accepted.depleting_pump_energy().is_none());
        assert!(reference.depleting_pump_energy().is_some());
        relative(
            report.beams[0].absorbed_j[0] / report.beams[1].absorbed_j[0],
            1.5,
        );
        near(
            report.beams.iter().map(|beam| beam.transmitted_j[0]).sum(),
            outgoing * watt,
        );
        relative(warm.energy_j(), accepted.thermal_input_j());
        relative(audit.external_heat_j, accepted.thermal_input_j());
        accounts(&accepted);
        accounts(&reference);
    });
}

#[test]
fn g1_counterpropagating_colors_share_ground_inventory_and_preserve_beam_order() {
    with_cx(|cx| {
        let density = [0.0, 1.0];
        let medium = grid(cx, &density, transparent_host());
        let (thermal, heating) = inventory(&medium, &[[COUNT, 0.0, 0.0]; 2], [1e12; 2], 0.0);
        let singlet = sections(2, [[1.0 / COUNT, 0.0], [0.0; 2]]);
        let triplet = sections(2, [[0.0; 2], [0.0, 1.0 / COUNT]]);
        let beams = [
            DepletingPumpBeam {
                source_identity: ContentHash([179; 32]),
                beam: ray(2, false, [400.0, 410.0], [2.0 * COUNT * HC / 400.0, 0.0]),
                cross_sections: &singlet,
                pump_window_s: None,
            },
            DepletingPumpBeam {
                source_identity: ContentHash([180; 32]),
                beam: ray(2, true, [390.0, 450.0], [0.0, 2.0 * COUNT * HC / 450.0]),
                cross_sections: &triplet,
                pump_window_s: None,
            },
        ];
        let original = bind(cx, &medium, &thermal, heating.clone(), &beams);
        let reordered = bind(
            cx,
            &medium,
            &thermal,
            heating,
            &[beams[1].clone(), beams[0].clone()],
        );
        let (_, a, _) = thermal
            .advance_observed(&original, 0.4, config(0.4), || Ok(()))
            .unwrap();
        let (_, b, _) = thermal
            .advance_observed(&reordered, 0.4, config(0.4), || Ok(()))
            .unwrap();
        let left = a.populations(0).unwrap();
        let right = a.populations(1).unwrap();
        near(left[0], right[0]);
        near(left[1], right[2]);
        near(left[2], right[1]);
        assert!(left[1] > left[2]);
        assert!(right[2] > right[1]);
        for index in 0..2 {
            relative(a.populations(index).unwrap().iter().sum(), COUNT);
            for (x, y) in a
                .populations(index)
                .unwrap()
                .into_iter()
                .zip(b.populations(index).unwrap())
            {
                near(x, y);
            }
        }
        let report = a.depleting_pump_report().unwrap();
        let reverse = b.depleting_pump_report().unwrap();
        assert_eq!(report.beams[0].wavelengths_nm, [400.0, 410.0]);
        assert_eq!(report.beams[1].wavelengths_nm, [390.0, 450.0]);
        assert_eq!(
            reverse.beams[0].source_identity,
            report.beams[1].source_identity
        );
        assert_eq!(
            reverse.beams[1].source_identity,
            report.beams[0].source_identity
        );
        near(
            report.beams[0].absorbed_j[0] / (HC / 400.0),
            report.beams[1].absorbed_j[1] / (HC / 450.0),
        );
        for beam in 0..2 {
            for lane in 0..2 {
                near(
                    report.beams[beam].absorbed_j[lane],
                    reverse.beams[1 - beam].absorbed_j[lane],
                );
                near(
                    report.beams[beam].transmitted_j[lane],
                    reverse.beams[1 - beam].transmitted_j[lane],
                );
            }
        }
        accounts(&a);
        accounts(&b);
    });
}

#[test]
fn g3_simultaneous_pump_windows_intersect_the_global_gate_and_retry_atomically() {
    with_cx(|cx| {
        let density = [1.0];
        let medium = grid(cx, &density, transparent_host());
        let (thermal, heating) = inventory(&medium, &[[COUNT, 0.0, 0.0]], [1.0, 5.0], 3.0);
        let cross_sections = sections(1, [[1.0 / COUNT, 0.0], [0.0; 2]]);
        let powers = [2.0 * COUNT * HC / 400.0, 3.0 * COUNT * HC / 400.0];
        let beams = [
            DepletingPumpBeam {
                source_identity: ContentHash([181; 32]),
                beam: ray(1, false, [400.0, 450.0], [powers[0], 0.0]),
                cross_sections: &cross_sections,
                pump_window_s: Some((0.0, 0.7)),
            },
            DepletingPumpBeam {
                source_identity: ContentHash([182; 32]),
                beam: ray(1, true, [400.0, 450.0], [powers[1], 0.0]),
                cross_sections: &cross_sections,
                pump_window_s: Some((0.5, 1.2)),
            },
        ];
        let heating = bind(cx, &medium, &thermal, heating, &beams);
        let initial_report = heating.depleting_pump_report().unwrap().clone();
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
        assert_eq!(candidate.depleting_pump_report().unwrap(), &initial_report);
        assert_eq!(candidate.populations(0), heating.populations(0));
        let mut checkpoints = 0;
        assert!(matches!(
            thermal.advance_observed(&heating, 0.2, config(0.2), || {
                checkpoints += 1;
                if checkpoints > 25 {
                    Err(DepositionError::Cancelled)
                } else {
                    Ok(())
                }
            }),
            Err(OpticalHeatingError::Thermal(DepositionError::Cancelled))
        ));
        assert!(checkpoints > 25);
        assert_eq!(heating.depleting_pump_report().unwrap(), &initial_report);

        let heating = heating.with_pump_window(0.2, 1.0).unwrap();
        let (_, whole, audit) = thermal
            .advance_observed(&heating, 1.5, config(1.5), || Ok(()))
            .unwrap();
        let (middle, first, _) = thermal
            .advance_observed(&heating, 0.7, config(1.5), || Ok(()))
            .unwrap();
        let (_, split, _) = middle
            .advance_observed(&first, 1.5, config(1.5), || Ok(()))
            .unwrap();
        assert!(audit.substeps >= 5);
        assert_eq!(whole.populations(0), split.populations(0));
        assert_eq!(whole.energy(), split.energy());
        assert_eq!(whole.depleting_pump_report(), split.depleting_pump_report());
        let report = whole.depleting_pump_report().unwrap();
        relative(report.duration_s, 1.5);
        for (beam, power) in report.beams.iter().zip(powers) {
            relative(beam.incident_j[0], 0.5 * power);
        }
        let (dark_start, charged, _) = thermal
            .advance_observed(&heating, 1.2, config(1.5), || Ok(()))
            .unwrap();
        let (_, dark, _) = dark_start
            .advance_observed(&charged, 1.5, config(1.5), || Ok(()))
            .unwrap();
        for (before, after) in charged
            .depleting_pump_report()
            .unwrap()
            .beams
            .iter()
            .zip(&dark.depleting_pump_report().unwrap().beams)
        {
            assert_eq!(before.incident_j, after.incident_j);
            assert_eq!(before.transmitted_j, after.transmitted_j);
            assert_eq!(before.absorbed_j, after.absorbed_j);
        }
        assert!(dark.energy().stored_excitation_j < charged.energy().stored_excitation_j);
        assert!(
            dark.energy().escaped_j.iter().sum::<f64>()
                > charged.energy().escaped_j.iter().sum::<f64>()
        );
        assert!(dark.cell_emission_w()[0].1.iter().all(|power| *power > 0.0));
        accounts(&whole);
        accounts(&dark);
    });
}

#[test]
fn g1_simultaneous_passive_pump_absorption_is_deposited_once_per_native_cell() {
    with_cx(|cx| {
        let density = [1.0];
        let medium = grid(
            cx,
            &density,
            MassAbsorptionSpectrum::source_table(
                &[(380.0, 0.3), (500.0, 0.3), (650.0, 0.0), (900.0, 0.0)],
                ContentHash([183; 32]),
            )
            .unwrap(),
        );
        let (thermal, heating) = inventory(&medium, &[[0.0; 3]], [1.0, 5.0], 0.0);
        let cross_sections = sections(1, [[0.0; 2]; 2]);
        let powers = [[1.0, 2.0], [3.0, 4.0]];
        let beams = [
            DepletingPumpBeam {
                source_identity: ContentHash([184; 32]),
                beam: ray(1, false, [400.0, 450.0], powers[0]),
                cross_sections: &cross_sections,
                pump_window_s: None,
            },
            DepletingPumpBeam {
                source_identity: ContentHash([185; 32]),
                beam: ray(1, true, [410.0, 420.0], powers[1]),
                cross_sections: &cross_sections,
                pump_window_s: None,
            },
        ];
        let heating = bind(cx, &medium, &thermal, heating, &beams);
        let (warm, accepted, audit) = thermal
            .advance_observed(&heating, 0.4, config(0.4), || Ok(()))
            .unwrap();
        let fraction = -(-0.3_f64).exp_m1();
        for (beam, powers) in accepted
            .depleting_pump_report()
            .unwrap()
            .beams
            .iter()
            .zip(powers)
        {
            for (lane, power) in powers.into_iter().enumerate() {
                relative(beam.passive_absorbed_j[lane], power * 0.4 * fraction);
                relative(beam.transmitted_j[lane], power * 0.4 * (-0.3_f64).exp());
            }
            assert_eq!(beam.absorbed_j, [0.0; 2]);
        }
        let heat = 4.0 * fraction;
        relative(warm.energy_j(), heat);
        relative(audit.external_heat_j, heat);
        relative(accepted.energy().passive_pump_absorption_j, heat);
        assert_eq!(accepted.populations(0), Some([0.0; 3]));
        assert_eq!(accepted.absorbed_pump_energy_j(), [0.0; 2]);
        assert_eq!(accepted.energy().emitted_j, [0.0; 2]);
        accounts(&accepted);
    });
}

#[test]
fn g3_local_excitation_survives_outside_beams_and_uses_the_original_global_pulse() {
    with_cx(|cx| {
        let density = [1.0, 1.0];
        let medium = grid(cx, &density, transparent_host());
        let thermal = thermal(
            vec![cell(1, 1.0, 1.0, 0.0), cell(2, 1.0, 1.0, 0.0)],
            &SENSIBLE,
            vec![],
        );
        let pump = PrescribedPhotonPump::try_new(
            ContentHash([186; 32]),
            [400.0, 450.0],
            [1e-20; 2],
            [2e20, 0.0],
        )
        .unwrap();
        let bindings: Vec<_> = (0..2)
            .map(|grid_cell| SaturableEmitterCell {
                grid_cell,
                thermal_cell: CellId(grid_cell as u64 + 1),
                initial_populations: [COUNT, 0.0, 0.0],
                pump,
            })
            .collect();
        let triplet = SingletTripletSource::try_new(
            source(600.0, 5.0, [0.0, 1.0], 187),
            0.0,
            ArrheniusQuenching::try_new(0.0, 0.0).unwrap(),
            SingletTripletLimits::default(),
        )
        .unwrap();
        let heating = SaturableLuminescentHeating::try_new(
            &medium,
            source(500.0, 1.0, [1.0, 0.0], 188),
            triplet,
            &bindings,
            &thermal,
            2,
            1e-7,
            None,
            || Ok(()),
        )
        .unwrap();
        let cross_sections = sections(1, [[1.0 / COUNT, 0.0], [0.0; 2]]);
        // Travel in y through cell 0, leaving cell 1 outside every beam.
        let power = COUNT * HC / 400.0;
        let mut incident = ray(2, false, [400.0, 450.0], [power, 0.0]);
        incident.origin = Point3::new(0.5, -1.0, 0.5);
        incident.unit_direction = Vec3::new(0.0, 1.0, 0.0);
        let heating = bind(
            cx,
            &medium,
            &thermal,
            heating,
            &[DepletingPumpBeam {
                source_identity: ContentHash([189; 32]),
                beam: incident,
                cross_sections: &cross_sections,
                pump_window_s: Some((0.4, 0.6)),
            }],
        )
        .with_pump_window(0.2, 0.8)
        .unwrap();
        let (warm, whole, audit) = thermal
            .advance_observed(&heating, 1.0, config(1.0), || Ok(()))
            .unwrap();
        let (middle, first, _) = thermal
            .advance_observed(&heating, 0.6, config(1.0), || Ok(()))
            .unwrap();
        let (_, split, _) = middle
            .advance_observed(&first, 1.0, config(1.0), || Ok(()))
            .unwrap();
        assert_eq!(whole.populations(0), split.populations(0));
        assert_eq!(whole.populations(1), split.populations(1));
        assert_eq!(whole.energy(), split.energy());
        assert_eq!(whole.depleting_pump_report(), split.depleting_pump_report());
        assert!(audit.substeps >= 5);
        // The unvisited cell has W=2/s, radiative decay=1/s. It is pumped for
        // .6 s under the GLOBAL window, then decays for .2 s. The individual
        // beam's narrower window must not suppress this local excitation.
        let expected = COUNT * (2.0 / 3.0) * -(-1.8_f64).exp_m1() * (-0.2_f64).exp();
        relative(whole.populations(1).unwrap()[1], expected);
        assert!(whole.populations(0).unwrap()[1] > whole.populations(1).unwrap()[1]);
        for index in 0..2 {
            relative(whole.populations(index).unwrap().iter().sum(), COUNT);
        }
        let beam = &whole.depleting_pump_report().unwrap().beams[0];
        relative(beam.incident_j[0], power * (0.6 - 0.4));
        relative(
            beam.incident_j[0],
            beam.transmitted_j[0] + beam.absorbed_j[0],
        );
        assert!(whole.absorbed_pump_energy_j()[0] > beam.absorbed_j[0]);
        relative(warm.energy_j(), whole.thermal_input_j());
        balance(&whole);
    });
}
