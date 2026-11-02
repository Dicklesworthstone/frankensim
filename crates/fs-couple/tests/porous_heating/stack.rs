use super::*;
use fs_duct::{compose::ComposeOptions, porous::PorousLayerStack};

const FIRST: f64 = 0.03;
const SECOND: f64 = LENGTH - FIRST;

fn layers<'m>(materials: &'m [PorousFluidLaw; 2]) -> [PorousLayer<'m>; 2] {
    let gas = gas();
    std::array::from_fn(|i| {
        PorousLayer::try_new(
            &materials[i],
            &gas,
            core::f64::consts::TAU * 1100.0,
            [FIRST, SECOND][i],
            AREA,
            gas.characteristic_impedance / AREA,
        )
        .unwrap()
    })
}
fn materials() -> [PorousFluidLaw; 2] {
    [
        material(),
        PorousFluidLaw::try_new(PorousCard::basotect_foam02(), 1.81e-5).unwrap(),
    ]
}
fn partition(n: usize) -> (Vec<f64>, [Vec<PorousHeatCell>; 2]) {
    let mut ends = Vec::new();
    let parts = std::array::from_fn(|layer| {
        (1..=n)
            .map(|i| {
                let length = [FIRST, SECOND][layer];
                let end_depth_m = if i == n {
                    length
                } else {
                    length * i as f64 / n as f64
                };
                ends.push(if layer == 0 {
                    end_depth_m
                } else {
                    FIRST + end_depth_m
                });
                PorousHeatCell {
                    cell: id(layer * n + i - 1),
                    end_depth_m,
                }
            })
            .collect()
    });
    (ends, parts)
}
fn multi_initial(
    mut schedule: Vec<TimedEvent>,
    second_start: usize,
) -> DepositionSimulation<'static> {
    for event in &mut schedule[second_start..] {
        if let Event::Birth(cell) = &mut event.event {
            cell.material_id = MaterialId(9);
            cell.mass_kg *= 2.0;
        }
    }
    DepositionSimulation::new(
        schedule,
        &[(MaterialId(5), &SENSIBLE), (MaterialId(9), &LATENT)],
        DepositionBudget {
            max_cells: 32,
            max_links: 32,
            max_events: 64,
            max_substeps: 10000,
        },
    )
    .unwrap()
    .advance_to(0.0, config(1.0), || Ok(()))
    .unwrap()
    .0
}

#[test]
fn g1_heterogeneous_stack_heats_each_native_material_with_its_internal_coherent_field() {
    let materials = materials();
    let layers = layers(&materials);
    let stack = PorousLayerStack::try_new(&layers, ComposeOptions::default(), || false).unwrap();
    let drive = [C64::from_re(100.0), C64::new(25.0, -50.0)];
    let excitation = stack.prepare_excitation(drive, || false).unwrap();
    let mut deposited = Vec::new();
    for n in [1, 3, 7] {
        let (ends, parts) = partition(n);
        let mut h = vec![0.0; n];
        h.extend(vec![100.0; n]);
        let thermal = multi_initial(events(&ends, &h, &vec![0.0; 2 * n], 0.0), n);
        let heating = PorousStackHeating::try_new(
            &excitation,
            &thermal,
            &[&parts[0], &parts[1]],
            2 * n,
            || Ok(()),
        )
        .unwrap();
        assert!(heating.layer_cells(2).is_none());
        for (i, field) in excitation.layer_excitations().iter().enumerate() {
            for cell in heating.layer_cells(i).unwrap() {
                let expected = quadrature(field, cell.start_depth_m, cell.end_depth_m);
                close(
                    cell.dissipation.viscous_power_w,
                    expected[0],
                    3e-11 * expected[0],
                );
                close(
                    cell.dissipation.thermal_power_w,
                    expected[1],
                    3e-11 * expected[1],
                );
            }
        }
        let step = heating
            .advance_to(&thermal, 0.25, config(0.037), || Ok(()))
            .unwrap();
        for (i, cell) in heating.cells().iter().enumerate() {
            let before = &thermal.cells()[&cell.cell];
            let after = &step.thermal.cells()[&cell.cell];
            let expected = 0.25 * cell.dissipation.absorbed_power_w;
            close(
                (after.specific_h_j_kg - before.specific_h_j_kg) * before.mass_kg,
                expected,
                2e-14,
            );
            close(step.cell_energy_j[i].1, expected, 2e-14);
            assert_eq!(after.material_id, before.material_id);
            assert_eq!(after.boundary.source_w, 0.0);
            if i >= n {
                assert!(after.specific_h_j_kg > 100.0);
                assert_eq!(step.thermal.temperature_k(cell.cell).unwrap(), 300.5);
            }
        }
        close(
            step.acoustic.incident_j,
            0.25 * excitation.boundary_power().incident_power_w,
            2e-14,
        );
        close(
            step.acoustic.outgoing_j,
            0.25 * excitation.boundary_power().outgoing_power_w,
            2e-14,
        );
        close(step.acoustic.residual_j, 0.0, 3e-14);
        close(
            step.acoustic.viscous_j + step.acoustic.thermal_j,
            step.acoustic.deposited_j,
            2e-14,
        );
        close(
            step.heat_report.external_heat_j,
            step.acoustic.deposited_j,
            2e-14,
        );
        deposited.push(step.acoustic.deposited_j);
    }
    for total in deposited {
        close(
            total,
            0.25 * excitation.dissipation().absorbed_power_w,
            3e-14,
        );
    }
    // The second slab's drive includes interface reflection and transmission;
    // applying the external drive to each isolated slab cannot produce its heat.
    let isolated = layers[1]
        .prepare_excitation(drive)
        .unwrap()
        .dissipation_between(0.0, SECOND)
        .unwrap();
    let actual = excitation.layer_excitations()[1]
        .dissipation_between(0.0, SECOND)
        .unwrap();
    assert!(
        (isolated.absorbed_power_w - actual.absorbed_power_w).abs()
            > 0.01 * actual.absorbed_power_w
    );
}

#[test]
fn g1_stack_contacts_authored_heat_and_continuation_share_one_native_balance() {
    let materials = materials();
    let layers = layers(&materials);
    let stack = PorousLayerStack::try_new(&layers, ComposeOptions::default(), || false).unwrap();
    let excitation = stack
        .prepare_excitation([C64::from_re(100.0), C64::ZERO], || false)
        .unwrap();
    let (ends, parts) = partition(1);
    let mut schedule = events(&ends, &[0.0, 100.0], &[0.3, -0.2], 0.1);
    let Event::Birth(cell) = &schedule[0].event else {
        panic!("birth")
    };
    let changed = Boundary {
        source_w: 0.6,
        ..cell.boundary
    };
    schedule.push(TimedEvent {
        time_s: 0.1,
        event: Event::Boundary(id(0), changed),
    });
    let thermal = multi_initial(schedule, 1);
    let heating =
        PorousStackHeating::try_new(&excitation, &thermal, &[&parts[0], &parts[1]], 2, || Ok(()))
            .unwrap();
    let first = heating
        .advance_to(&thermal, 0.1, config(0.1), || Ok(()))
        .unwrap();
    assert_eq!(first.heat_report.substeps, 1);
    // The thermal contact carries 0.1*(300.5-300) W across the material face.
    for (i, external) in [0.35, -0.25].into_iter().enumerate() {
        let before = &thermal.cells()[&id(i)];
        let expected = 0.1 * (heating.cells()[i].dissipation.absorbed_power_w + external);
        close(
            (first.thermal.cells()[&id(i)].specific_h_j_kg - before.specific_h_j_kg)
                * before.mass_kg,
            expected,
            2e-14,
        );
    }
    assert_eq!(first.thermal.links(), thermal.links());
    assert_eq!(first.heat_report.events, 1);
    assert_eq!(first.thermal.cells()[&id(0)].boundary.source_w, 0.6);
    let second = heating
        .advance_to(&first.thermal, 0.3, config(0.1), || Ok(()))
        .unwrap();
    let together = heating
        .advance_to(&thermal, 0.3, config(0.1), || Ok(()))
        .unwrap();
    for cell in [id(0), id(1)] {
        close(
            second.thermal.cells()[&cell].specific_h_j_kg,
            together.thermal.cells()[&cell].specific_h_j_kg,
            2e-12,
        );
    }
    close(
        together.heat_report.external_heat_j,
        together.acoustic.deposited_j + 0.1 * 0.1 + 0.4 * 0.2,
        2e-14,
    );
    close(
        first.acoustic.deposited_j + second.acoustic.deposited_j,
        together.acoustic.deposited_j,
        2e-14,
    );
    assert_eq!(thermal.time_s(), 0.0);
}

#[test]
fn g0_one_layer_stack_retains_original_heat_bits_and_cancelled_stack_retains_state() {
    let material = material();
    let layers = [layer(&material)];
    let stack = PorousLayerStack::try_new(&layers, ComposeOptions::default(), || false).unwrap();
    let excitation = stack
        .prepare_excitation([C64::from_re(100.0), C64::new(0.0, 20.0)], || false)
        .unwrap();
    let ends = [LENGTH / 2.0, LENGTH];
    let thermal = initial(events(&ends, &[20.0, 0.0], &[0.3, -0.2], 1.0), &SENSIBLE);
    let routes = routes(&ends);
    let single = PorousHeating::try_new(
        &excitation.layer_excitations()[0],
        &thermal,
        &routes,
        2,
        || Ok(()),
    )
    .unwrap();
    let stack_heat =
        PorousStackHeating::try_new(&excitation, &thermal, &[&routes], 2, || Ok(())).unwrap();
    let a = single
        .advance_to(&thermal, 0.3, config(0.03), || Ok(()))
        .unwrap();
    let b = stack_heat
        .advance_to(&thermal, 0.3, config(0.03), || Ok(()))
        .unwrap();
    assert_eq!(a.thermal.cells(), b.thermal.cells());
    assert_eq!(a.heat_report, b.heat_report);
    assert_eq!(a.acoustic, b.acoustic);
    assert_eq!(a.cell_energy_j, b.cell_energy_j);
    let mut polls = 0;
    assert!(matches!(
        stack_heat.advance_to(&thermal, 0.3, config(0.03), || {
            polls += 1;
            if polls > 40 {
                Err(DepositionError::Cancelled)
            } else {
                Ok(())
            }
        }),
        Err(PorousHeatingError::Thermal(DepositionError::Cancelled))
    ));
    let retry = stack_heat
        .advance_to(&thermal, 0.3, config(0.03), || Ok(()))
        .unwrap();
    assert_eq!(retry.thermal.cells(), b.thermal.cells());
    assert_eq!(thermal.time_s(), 0.0);
}

#[test]
fn g0_stack_mapping_requires_complete_distinct_adjoining_material_volumes() {
    let materials = materials();
    let layers = layers(&materials);
    let stack = PorousLayerStack::try_new(&layers, ComposeOptions::default(), || false).unwrap();
    let excitation = stack.prepare_excitation([C64::ZERO; 2], || false).unwrap();
    let (ends, parts) = partition(1);
    let schedule = events(&ends, &[0.0, 100.0], &[0.0, 0.0], 0.0);
    let thermal = multi_initial(schedule.clone(), 1);
    for bad in [&[][..], &[&parts[0][..]][..], &[&parts[0][..], &[][..]][..]] {
        assert!(PorousStackHeating::try_new(&excitation, &thermal, bad, 2, || Ok(())).is_err());
    }
    assert!(
        PorousStackHeating::try_new(&excitation, &thermal, &[&parts[0], &parts[1]], 1, || Ok(()))
            .is_err()
    );
    for mutation in 0..3 {
        let mut changed = schedule.clone();
        let Event::Birth(cell) = &mut changed[1].event else {
            panic!("birth")
        };
        match mutation {
            0 => {
                cell.geometry.start_m[2] += 0.01;
                cell.geometry.end_m[2] += 0.01;
            }
            1 => cell.body_id += 1,
            _ => {
                cell.geometry.end_m = [
                    cell.geometry.start_m[0] - 0.8 * SECOND,
                    cell.geometry.start_m[1] + 0.6 * SECOND,
                    cell.geometry.start_m[2],
                ]
            }
        }
        let changed = multi_initial(changed, 1);
        assert!(
            PorousStackHeating::try_new(&excitation, &changed, &[&parts[0], &parts[1]], 2, || Ok(
                ()
            ))
            .is_err()
        );
    }
    let heating =
        PorousStackHeating::try_new(&excitation, &thermal, &[&parts[0], &parts[1]], 2, || Ok(()))
            .unwrap();
    let zero = heating
        .advance_to(&thermal, 0.3, config(0.03), || Ok(()))
        .unwrap();
    let native = thermal.advance_to(0.3, config(0.03), || Ok(())).unwrap();
    assert_eq!(zero.thermal.cells(), native.0.cells());
    assert_eq!(zero.heat_report, native.1);
    assert_eq!(zero.acoustic.deposited_j, 0.0);
    let mut changed = schedule;
    let Event::Birth(cell) = &mut changed[1].event else {
        panic!("birth")
    };
    cell.mass_kg *= 1.1;
    let changed = multi_initial(changed, 1);
    assert!(
        heating
            .advance_to(&changed, 0.3, config(0.03), || Ok(()))
            .is_err()
    );
}
