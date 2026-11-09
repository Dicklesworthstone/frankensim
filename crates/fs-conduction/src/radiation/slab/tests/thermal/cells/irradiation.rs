use super::*;
use fs_thermal::enthalpy::conduction::implicit::source::{
    HeatSourceSchedule, HeatSourceScheduleKnot,
};

fn schedule(points: &[(f64, f64)]) -> HeatSourceSchedule {
    HeatSourceSchedule::new(
        points
            .iter()
            .map(|&(time_s, power_multiplier)| HeatSourceScheduleKnot {
                time_s,
                power_multiplier,
            })
            .collect(),
        points.len(),
        || Ok(()),
    )
    .unwrap()
}

fn gray<'a>(
    response: &SlabCellResponse,
    law: &'a Caloric,
    old: &DepositionSimulation<'_>,
    area: f64,
) -> GraySlabThermal<'a> {
    GraySlabThermal::new_cells_full_spectrum(
        response,
        area,
        response
            .slab()
            .exterior_refractive_indices()
            .map(|n| n * n * exitance(500.0)),
        &bindings(law, response.cell_count()),
        old,
        || Ok(()),
    )
    .unwrap()
}

#[test]
fn g1_oblique_pulse_heats_original_depth_cells_with_exact_native_joules() {
    let slab = slab(&[0.8]);
    let response = partition(&slab, 0.25);
    let law = Caloric {
        capacity: 1000.0,
        datum: 0.0,
        latent: 20_000.0,
    };
    let old = cell_state(&law, &response, &[510_000.0; 2], 0.01);
    let beam = slab.beam_response(0.5, 0.35, || false).unwrap();
    let source = gray(&response, &law, &old, 0.01)
        .with_beams(
            vec![SlabThermalBeam {
                id: 7,
                band: 0,
                response: beam,
                irradiation_w_m2: [1000.0, 0.0],
                schedule: schedule(&[(0.0, 0.0), (0.1, 2.0), (0.4, 0.0), (0.6, 0.0)]),
            }],
            4,
            &old,
            || Ok(()),
        )
        .unwrap();
    let report = source
        .beam_illumination()
        .unwrap()
        .interval(0.0, 0.6, 0.6, || Ok(()))
        .unwrap();
    // Independent matched-index Beer integral: pulse area is0.4s and normal
    // face power is10W. No second cosine or diffuse angular average belongs here.
    let incident = 4.0;
    let expected = [
        incident * -(-0.8_f64 * 0.25 / 0.5).exp_m1(),
        incident * (-0.8_f64 * 0.25 / 0.5).exp() * -(-0.8_f64 * 0.75 / 0.5).exp_m1(),
    ];
    close(report.incident_j[0], incident, 1e-14);
    close(report.absorbed_j, expected.iter().sum(), 2e-13);
    assert_eq!(report.reflected_j, [0.0; 2]);
    close(report.transmitted_j[1], incident * (-1.6_f64).exp(), 1e-13);
    close(report.balance_residual_j, 0.0, 1e-13);
    let mut coarse_config = config();
    coarse_config.max_dt_s = 0.6;
    let (next, _, audit) = old
        .advance_observed(&source, 0.6, coarse_config, || Ok(()))
        .unwrap();
    close(audit.external_heat_j, report.absorbed_j, 2e-12);
    for i in 0..2 {
        let id = CellId(i as u64);
        let mass = old.cells()[&id].mass_kg;
        close(
            next.cells()[&id].specific_h_j_kg,
            old.cells()[&id].specific_h_j_kg + expected[i] / mass,
            2e-9,
        );
        assert_eq!(next.temperature_k(id).unwrap(), 500.0);
        assert_eq!(next.cells()[&id].geometry, old.cells()[&id].geometry);
    }
    assert!(next.links().is_empty());
    assert_eq!(source.slab().layer_count(), 1);
    assert!(source.requires_step_interval());
    assert!(source.endpoint_transport(&old, || Ok(())).is_err());
}

#[test]
fn g1_multiple_polarized_angles_keep_asymmetric_layer_and_directional_accounts() {
    let slab = asymmetric::between(&[0.7, 2.0], &[1.5, 1.8], [1.3, 1.1]);
    let response = partition(&slab, 0.3);
    let law = Caloric {
        capacity: 1000.0,
        datum: 0.0,
        latent: 20_000.0,
    };
    let old = cell_state(&law, &response, &[510_000.0; 4], 0.02);
    let front = SlabThermalBeam {
        id: 20,
        band: 0,
        response: slab.beam_response_from_side(0, 0.7, 1.0, || false).unwrap(),
        irradiation_w_m2: [20.0, 0.0],
        schedule: schedule(&[(0.0, 0.0), (0.2, 1.0), (1.0, 3.0)]),
    };
    let back = SlabThermalBeam {
        id: 10,
        band: 0,
        response: slab
            .beam_response_from_side(1, 0.85, 0.0, || false)
            .unwrap(),
        irradiation_w_m2: [0.0, 70.0],
        schedule: schedule(&[(0.0, 2.0), (1.0, 2.0)]),
    };
    let source = gray(&response, &law, &old, 0.02)
        .with_beams(vec![front.clone(), back.clone()], 5, &old, || Ok(()))
        .unwrap();
    let prepared = source.beam_illumination().unwrap();
    assert_eq!(prepared.beams().map(|b| b.id).collect::<Vec<_>>(), [10, 20]);
    let report = prepared.interval(0.0, 1.0, 1.0, || Ok(())).unwrap();
    close(report.incident_j[0], 0.02 * 20.0 * 1.7, 1e-14);
    close(report.incident_j[1], 0.02 * 70.0 * 2.0, 1e-14);
    let mut expected_layers = [0.0; 2];
    for (beam, integral) in [(&front, 1.7), (&back, 2.0)] {
        let mut powers = [0.0; 2];
        beam.response
            .exchange(beam.irradiation_w_m2, &mut powers)
            .unwrap();
        for i in 0..2 {
            expected_layers[i] += 0.02 * integral * powers[i];
        }
    }
    for (layer, expected) in expected_layers.iter().enumerate() {
        close(
            report.mean_sources_w[&CellId(2 * layer as u64)]
                + report.mean_sources_w[&CellId(2 * layer as u64 + 1)],
            *expected,
            2e-12,
        );
    }
    close(report.balance_residual_j, 0.0, 2e-12);
    let (_, _, audit) = old
        .advance_observed(&source, 1.0, config(), || Ok(()))
        .unwrap();
    close(audit.external_heat_j, report.absorbed_j, 2e-11);
    assert_eq!(
        prepared.beams().next().unwrap().response.layers().as_ptr(),
        slab.layers().as_ptr()
    );
}

#[test]
fn g1_spectral_be_uses_pulse_integral_and_actual_endpoint_planck_exchange() {
    let slabs = [
        asymmetric::between(&[0.4], &[1.5], [1.2, 1.7]),
        asymmetric::between(&[0.9], &[1.8], [1.4, 1.1]),
    ];
    let responses = slabs.each_ref().map(|s| partition(s, 0.4));
    let windows = [
        BlackbodyBand::new(1e-6, 5e-6).unwrap(),
        BlackbodyBand::new(8e-6, 20e-6).unwrap(),
    ];
    let law = Caloric {
        capacity: 1000.0,
        datum: 0.0,
        latent: 0.0,
    };
    let known = cell_state(&law, &responses[0], &[550_000.0, 350_000.0], 0.01);
    let base = SpectralSlabThermal::new_cells(
        (0..2)
            .map(|b| SlabThermalCellBand {
                window: windows[b],
                response: responses[b].clone(),
                irradiation_w_m2: [0.0; 2],
            })
            .collect(),
        0.01,
        &bindings(&law, 2),
        &known,
        || Ok(()),
    )
    .unwrap();
    let thermal = base.endpoint_transport(&known, || Ok(())).unwrap();
    let source = base
        .with_beams(
            (0..2)
                .map(|b| SlabThermalBeam {
                    id: b as u64,
                    band: b,
                    response: slabs[b]
                        .beam_response_from_side(b, 0.8, b as f64, || false)
                        .unwrap(),
                    irradiation_w_m2: if b == 0 { [5000.0, 0.0] } else { [0.0, 7000.0] },
                    // Both endpoints are dark; an endpoint-amplitude approximation loses
                    // the ENTIRE supplied pulse. Its integral is2s at unit base power.
                    schedule: schedule(&[(0.0, 0.0), (0.25, 4.0), (1.0, 0.0)]),
                })
                .collect(),
            6,
            &known,
            || Ok(()),
        )
        .unwrap();
    let beam = source
        .beam_illumination()
        .unwrap()
        .interval(0.0, 1.0, 1.0, || Ok(()))
        .unwrap();
    let mut increments = [0.0; 2];
    for i in 0..2 {
        let id = CellId(i as u64);
        increments[i] =
            thermal.heat.bulk.sources_w.get(&id).copied().unwrap_or(0.0) + beam.mean_sources_w[&id];
    }
    assert_eq!(thermal.transfers.len(), 1);
    increments[0] -= thermal.transfers[0].power_w;
    increments[1] += thermal.transfers[0].power_w;
    let old_energies: Vec<_> = (0..2)
        .map(|i| {
            let cell = &known.cells()[&CellId(i as u64)];
            cell.specific_h_j_kg - increments[i] / cell.mass_kg
        })
        .collect();
    let old = cell_state(&law, &responses[0], &old_energies, 0.01);
    let step = BackwardEulerStep::new(&old, 1.0, || Ok(())).unwrap();
    let trial = step
        .evaluate_with_transport(
            &increments,
            |candidate| {
                source.endpoint_transport_over_interval(candidate, 0.0, 1.0, 1.0, || Ok(()))
            },
            || Ok(()),
        )
        .unwrap();
    assert!(trial.maximum_cell_residual_j() < 1e-8);
    let (next, (), audit) = trial.accept(1e-8, 1e-8, || Ok(())).unwrap();
    close(next.temperature_k(CellId(0)).unwrap(), 550.0, 1e-12);
    close(next.temperature_k(CellId(1)).unwrap(), 350.0, 1e-12);
    close(
        audit.external_heat_j,
        thermal.heat.bulk.sources_w.values().sum::<f64>() + beam.absorbed_j,
        2e-11,
    );
    assert!(source.endpoint_transport(&known, || Ok(())).is_err());
    assert!(
        source
            .endpoint_transport_over_interval(&known, 0.0, 1.0, 1.0, || Ok(()))
            .is_err()
    );
    assert_eq!(source.bands()[0].optics.layer_count(), 1);
}

#[test]
fn g1_closed_snell_channel_reflects_energy_without_native_heat() {
    let slab = asymmetric::between(&[0.9], &[1.0], [2.0, 1.0]);
    let response = partition(&slab, 0.5);
    let law = Caloric {
        capacity: 1000.0,
        datum: 0.0,
        latent: 0.0,
    };
    let old = cell_state(&law, &response, &[500_000.0; 2], 0.01);
    let beam = slab.beam_response(0.5, 0.5, || false).unwrap();
    assert_eq!(beam.exterior_direction_cosines()[1], 0.0);
    let input = SlabThermalBeam {
        id: 0,
        band: 0,
        response: beam,
        irradiation_w_m2: [100.0, 0.0],
        schedule: schedule(&[(0.0, 1.0), (1.0, 1.0)]),
    };
    let source = gray(&response, &law, &old, 0.01)
        .with_beams(vec![input.clone()], 2, &old, || Ok(()))
        .unwrap();
    let report = source
        .beam_illumination()
        .unwrap()
        .interval(0.0, 1.0, 1.0, || Ok(()))
        .unwrap();
    assert!(report.mean_sources_w.is_empty());
    assert_eq!(report.absorbed_j, 0.0);
    assert_eq!(report.incident_j, report.reflected_j);
    assert_eq!(report.transmitted_j, [0.0; 2]);
    let (next, _, audit) = old
        .advance_observed(&source, 1.0, config(), || Ok(()))
        .unwrap();
    assert_eq!(audit.external_heat_j, 0.0);
    assert_eq!(next.cells()[&CellId(0)].specific_h_j_kg, 500_000.0);
    let mut invalid = input;
    invalid.irradiation_w_m2[1] = 1.0;
    assert!(
        gray(&response, &law, &old, 0.01)
            .with_beams(vec![invalid], 2, &old, || Ok(()))
            .is_err()
    );
}

#[test]
fn g3_weak_narrow_beam_cells_retain_positive_heat_and_original_time_domain() {
    let slab = slab(&[1e-20]);
    let response = partition(&slab, 1e-9);
    let law = Caloric {
        capacity: 1000.0,
        datum: 0.0,
        latent: 0.0,
    };
    let old = cell_state(&law, &response, &[500_000.0; 2], 0.01);
    let source = gray(&response, &law, &old, 0.01)
        .with_beams(
            vec![SlabThermalBeam {
                id: 0,
                band: 0,
                response: slab.beam_response(0.5, 0.0, || false).unwrap(),
                irradiation_w_m2: [1e20, 0.0],
                schedule: schedule(&[(0.0, 1.0), (1.0, 1.0)]),
            }],
            2,
            &old,
            || Ok(()),
        )
        .unwrap();
    let illumination = source.beam_illumination().unwrap();
    let receipt = illumination.interval(0.0, 1.0, 1.0, || Ok(())).unwrap();
    let value = receipt.mean_sources_w[&CellId(0)];
    assert!(value > 0.0);
    close(value / 2e-11, 1.0, 1e-12);
    assert!(illumination.interval(0.0, 1.1, 1.1, || Ok(())).is_err());
    assert!(illumination.interval(0.0, 0.5, 1.0, || Ok(())).is_err());
    let mut polls = 0;
    assert_eq!(
        illumination.interval(0.0, 1.0, 1.0, || {
            polls += 1;
            if polls == 3 {
                Err(DepositionError::Cancelled)
            } else {
                Ok(())
            }
        }),
        Err(DepositionError::Cancelled)
    );
    assert_eq!(
        illumination.interval(0.0, 1.0, 1.0, || Ok(())).unwrap(),
        receipt
    );
}

#[test]
fn g3_beams_admit_original_sources_and_bounded_work_before_attachment() {
    let slab = slab(&[0.5]);
    let response = partition(&slab, 0.5);
    let law = Caloric {
        capacity: 1000.0,
        datum: 0.0,
        latent: 20_000.0,
    };
    let old = cell_state(&law, &response, &[510_000.0; 2], 0.01);
    let source = gray(&response, &law, &old, 0.01);
    let beam = SlabThermalBeam {
        id: 1,
        band: 0,
        response: slab.beam_response(1.0, 0.5, || false).unwrap(),
        irradiation_w_m2: [0.0; 2],
        schedule: schedule(&[(0.0, 0.0), (1.0, 0.0)]),
    };
    assert!(
        source
            .clone()
            .with_beams(vec![beam.clone()], 1, &old, || Ok(()))
            .is_err()
    );
    assert!(
        source
            .clone()
            .with_beams(vec![beam.clone(); 2], 4, &old, || Ok(()))
            .is_err()
    );
    let mut wrong = beam.clone();
    wrong.band = 1;
    assert!(
        source
            .clone()
            .with_beams(vec![wrong], 2, &old, || Ok(()))
            .is_err()
    );
    let other = AbsorbingSlab::new(
        vec![
            AbsorbingLayer::new(
                1.0,
                0.5,
                [100.0, 3000.0],
                "different original optical source",
            )
            .unwrap(),
        ],
        SlabQuadrature::default(),
        || false,
    )
    .unwrap();
    let mut wrong = beam.clone();
    wrong.response = other.beam_response(1.0, 0.5, || false).unwrap();
    assert!(
        source
            .clone()
            .with_beams(vec![wrong], 2, &old, || Ok(()))
            .is_err()
    );
    let too_many: Vec<_> = (0..=MAX_SLAB_THERMAL_BEAMS)
        .map(|id| {
            let mut b = beam.clone();
            b.id = id as u64;
            b
        })
        .collect();
    assert!(
        source
            .clone()
            .with_beams(too_many, 1000, &old, || Ok(()))
            .is_err()
    );
    // An underflowing positive schedule mean must never be reported as a
    // physically dark beam, even when a later dimensional scale is large.
    let mut tiny = beam.clone();
    tiny.irradiation_w_m2 = [1e20, 0.0];
    tiny.schedule = schedule(&[(0.0, f64::from_bits(1)), (1.0, 0.0)]);
    let tiny = source
        .clone()
        .with_beams(vec![tiny], 2, &old, || Ok(()))
        .unwrap();
    assert_eq!(
        tiny.beam_illumination()
            .unwrap()
            .interval(0.0, 1.0, 1.0, || Ok(())),
        Err(DepositionError::Arithmetic)
    );
    let attached = source
        .clone()
        .with_beams(vec![beam], 2, &old, || Ok(()))
        .unwrap();
    assert!(
        old.advance_observed(&attached, 1.1, config(), || Ok(()))
            .is_err()
    );
    assert!(source.beam_illumination().is_none());
    assert_eq!(old.time_s(), 0.0);
}
