//! Original 3D vector memory participates in native circuit terminal work.
#[allow(dead_code)]
#[path = "hysteretic_common.rs"]
mod data;

use data::*;
use fs_em::{
    spatial::vector_potential::induction::{
        hysteretic::winding::{SpatialHystereticWindingCircuit, SpatialHystereticWindingError},
        nonlinear::winding::SpatialNonlinearWindings,
    },
    winding::{
        end::EndWindingMode,
        nonlinear::{NonlinearWindingPort, resistors::GraphResistorMaterial},
    },
};
use fs_phs::circuit::{Branch, CircuitGraph};

#[test]
fn voltage_driven_play_loop_separates_pinning_eddy_copper_graph_and_end_energy() {
    let card = material::card(500.0);
    let source = material::conductor(&card, 300.0);
    let problem = cube(&card);
    let original = circuit(&problem, &source);
    let mut accepted = original.clone();
    let mut old = Mode::default();
    let mut minor_return = None;
    near(
        accepted.resistance_ohm()[0],
        8.0 * material::resistivity(300.0),
        1e-12,
    );
    near(accepted.end_inductance_matrix_h()[0], END_L, 1e-15);
    for (index, h) in [12.0, 0.0, 6.0, 1.0, 6.0].into_iter().enumerate() {
        let expected = at_field(old, 300.0, h);
        let input = voltage(old, expected, 300.0, GRAPH_R);
        let candidate = accepted
            .step_to_candidate((index + 1) as f64 * DT, &[input], controls(), &mut || false)
            .unwrap();
        let next = candidate.circuit();
        let step = candidate.step();
        assert_field(next.field(), expected);
        let core = core_heat(old, expected, 300.0);
        let pinning = pinning_heat(old, expected);
        let copper = DT * 8.0 * material::resistivity(300.0) * expected.i.powi(2);
        let graph = DT * GRAPH_R * expected.i.powi(2);
        let change = stored(expected) - stored(old);
        let field_numerical = expected.h * (expected.b - old.b) / 8.0 - change - pinning;
        near(next.currents_a()[0], expected.i, 3e-7);
        near(
            next.flux_linkages_wb()[0],
            expected.b / 96.0 + END_L * expected.i,
            3e-7,
        );
        near(
            step.terminal_voltage_v[0],
            input - GRAPH_R * expected.i,
            3e-7,
        );
        near(step.field.joule_heat_j, core, 3e-8);
        near(step.field.hysteresis_heat_j, pinning, 3e-8);
        near(step.field.heat_j, core + pinning, 3e-8);
        near(step.field.numerical_dissipation_j, field_numerical, 3e-8);
        near(
            step.cells.iter().map(|cell| cell.heat_j).sum(),
            copper,
            4e-7,
        );
        near(step.graph_resistors[0].heat_j, graph, 4e-7);
        near(step.heat_j, core + pinning + copper + graph, 6e-7);
        near(
            step.stored_before_j,
            stored(old) + 0.5 * END_L * old.i.powi(2),
            4e-7,
        );
        near(
            step.stored_after_j,
            stored(expected) + 0.5 * END_L * expected.i.powi(2),
            4e-7,
        );
        near(
            step.numerical_dissipation_j,
            field_numerical + 0.5 * END_L * (expected.i - old.i).powi(2),
            4e-7,
        );
        near(step.supplied_work_j, DT * input * expected.i, 7e-7);
        near(step.energy_residual_j, 0.0, controls().energy_tolerance_j);
        assert_eq!(step.field.cells.len(), 4);
        assert!(step.field.maximum_charge_residual_a < 2e-8);
        for cell in &step.field.cells {
            near(cell.heat_j, core / 4.0, 3e-8);
        }
        if index == 0 {
            assert!(core > 0.0 && pinning > 0.0 && copper > 0.0 && graph > 0.0);
            assert!(field_numerical > 0.0);
        } else {
            assert_eq!(step.graph_materials.len(), 1);
        }
        if index == 2 {
            minor_return = Some(next.field().cell_states().to_vec());
        } else if index == 4 {
            for (returned, previous) in next
                .field()
                .cell_states()
                .iter()
                .zip(minor_return.as_ref().unwrap())
            {
                for axis in 0..3 {
                    near(
                        returned.polarizations_t()[0][axis],
                        previous.polarizations_t()[0][axis],
                        3e-8,
                    );
                }
            }
        }
        assert_eq!(accepted.time_s(), index as f64 * DT);
        old = expected;
        accepted = candidate.into_parts().0;
        if index == 0 {
            let memory = accepted.field().cell_states().to_vec();
            accepted = accepted
                .with_graph_resistors(
                    vec![GraphResistorMaterial {
                        branch: 1,
                        conductor: source
                            .clone()
                            .into_uniform_conductor(GRAPH_R / material::resistivity(300.0), 1.0)
                            .unwrap(),
                    }],
                    &mut || false,
                )
                .unwrap();
            assert_eq!(accepted.field().cell_states(), memory);
            assert_eq!(accepted.end_inductance(), original.end_inductance());
        }
    }
    assert_eq!(original.time_s(), 0.0);
    assert_eq!(original.currents_a(), [0.0]);
    assert_field(original.field(), Mode::default());
}

#[test]
fn current_driven_open_secondary_retains_remanent_flux_and_mutual_end_emf() {
    let card = material::card(500.0);
    let source = material::conductor(&card, 300.0);
    let problem = cube(&card);
    let first = profile(&problem);
    let profiles = vec![first.clone(), first.iter().map(|j| j.map(|v| -v)).collect()];
    let mut accepted = SpatialHystereticWindingCircuit::new_with_end_modes(
        field(&problem, &source, false),
        CircuitGraph {
            node_count: 3,
            branches: vec![
                (1, 0, Branch::CurrentSource { port: 0 }),
                (1, 0, Branch::Inductor { henries: f64::NAN }),
                (2, 0, Branch::CurrentSource { port: 1 }),
                (2, 0, Branch::Inductor { henries: f64::NAN }),
            ],
            transformers: vec![],
        },
        SpatialNonlinearWindings {
            branch_indices: vec![1, 3],
            profiles_per_m2: profiles.clone(),
            end_inductance_h: vec![0.0; 2],
        },
        vec![
            NonlinearWindingPort::VoltageInput {
                initial_voltage_v: 0.0,
                residual_scale_ohm: 1.0,
                current_tolerance_a: 2e-7,
            };
            2
        ],
        vec![EndWindingMode {
            inductance_h: 0.03,
            turns: vec![1.0, -2.0],
        }],
        strands(&profiles, &source, 0.5),
        initial(2),
        config(),
        &mut || false,
    )
    .unwrap();
    let mut old = Mode::default();
    let mut old_linkage = [0.0; 2];
    let mut old_end = 0.0;
    for (index, imposed) in [[144.0, 0.0], [0.0, 0.0], [60.0, -12.0]]
        .into_iter()
        .enumerate()
    {
        let h = (imposed[0] - imposed[1]) / 12.0;
        let z = material::play(old.z, h);
        let next = Mode {
            h,
            z,
            b: material::MU * h + z,
            i: imposed[0] - imposed[1],
        };
        let end = imposed[0] - 2.0 * imposed[1];
        let linkage = [next.b / 96.0 + 0.03 * end, -next.b / 96.0 - 0.06 * end];
        let candidate = accepted
            .step_with_conductors_candidate(
                (index + 1) as f64 * DT,
                vec![None; 48],
                accepted.conductors().to_vec(),
                &imposed,
                controls(),
                &mut || false,
            )
            .unwrap();
        let step = candidate.step();
        assert_field(candidate.circuit().field(), next);
        for winding in 0..2 {
            near(step.currents_a[winding], imposed[winding], 2e-7);
            near(step.flux_linkages_wb[1][winding], linkage[winding], 3e-7);
            near(
                step.terminal_voltage_v[winding],
                16.0 * material::resistivity(300.0) * imposed[winding]
                    + (linkage[winding] - old_linkage[winding]) / DT,
                3e-7,
            );
        }
        let pinning = pinning_heat(old, next);
        let copper =
            DT * 16.0 * material::resistivity(300.0) * (imposed[0].powi(2) + imposed[1].powi(2));
        near(step.field.joule_heat_j, 0.0, 1e-14);
        near(step.field.hysteresis_heat_j, pinning, 3e-8);
        near(step.heat_j, copper + pinning, 4e-7);
        near(step.stored_after_j, stored(next) + 0.015 * end * end, 4e-7);
        near(
            step.numerical_dissipation_j,
            next.h * (next.b - old.b) / 8.0 - stored(next) + stored(old) - pinning
                + 0.015 * (end - old_end).powi(2),
            4e-7,
        );
        near(step.energy_residual_j, 0.0, controls().energy_tolerance_j);
        if index == 0 {
            assert!(step.terminal_voltage_v[1].abs() > 1.0);
            near(
                step.cells
                    .iter()
                    .filter(|cell| cell.winding == 1)
                    .map(|cell| cell.heat_j)
                    .sum(),
                0.0,
                1e-14,
            );
        } else if index == 1 {
            near(next.b, material::CHI * material::PIN, 1e-14);
            assert!(step.field.hysteresis_heat_j > 0.0);
        }
        old = next;
        old_linkage = linkage;
        old_end = end;
        accepted = candidate.into_parts().0;
    }
}

#[test]
fn cancelled_or_foreign_source_trials_leave_magnetized_native_history_unchanged() {
    let card = material::card(500.0);
    let source = material::conductor(&card, 300.0);
    let problem = cube(&card);
    let accepted = circuit(&problem, &source)
        .step_to_candidate(DT, &[180.0], controls(), &mut || false)
        .unwrap()
        .into_parts()
        .0;
    let before = (
        accepted.state().to_vec(),
        accepted.field().cell_states().to_vec(),
        accepted.flux_linkages_wb().to_vec(),
    );
    let mut polls = 0;
    let good = accepted
        .step_to_candidate(2.0 * DT, &[-90.0], controls(), &mut || {
            polls += 1;
            false
        })
        .unwrap();
    assert!(polls > 20);
    let mut seen = 0;
    assert!(matches!(
        accepted.step_to_candidate(2.0 * DT, &[-90.0], controls(), &mut || {
            seen += 1;
            seen >= polls - 1
        }),
        Err(SpatialHystereticWindingError::Cancelled)
    ));
    let foreign = material::card_with_heat_capacity(500.0, CAPACITY + 1.0);
    let foreign_source = material::conductor(&foreign, 300.0);
    let core = (0..48)
        .map(|cell| (26..30).contains(&cell).then(|| foreign_source.clone()))
        .collect();
    assert!(
        accepted
            .step_with_conductors_candidate(
                2.0 * DT,
                core,
                accepted.conductors().to_vec(),
                &[-90.0],
                controls(),
                &mut || false,
            )
            .is_err()
    );
    assert_eq!(accepted.state(), before.0);
    assert_eq!(accepted.field().cell_states(), before.1);
    assert_eq!(accepted.flux_linkages_wb(), before.2);
    assert_eq!(accepted.time_s(), DT);
    let retry = accepted
        .step_to_candidate(2.0 * DT, &[-90.0], controls(), &mut || false)
        .unwrap();
    assert_eq!(retry.circuit().state(), good.circuit().state());
    assert_eq!(
        retry.circuit().field().cell_states(),
        good.circuit().field().cell_states()
    );
}
