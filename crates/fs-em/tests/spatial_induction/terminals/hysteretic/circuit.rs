//! Native capacitor storage drives actual distributed current and vector memory.
#[path = "../circuit/common.rs"]
mod common;
use super::*;
use common::{bracket, config, controls, initial};
use fs_em::spatial::vector_potential::induction::hysteretic::terminals::circuit::{
    SpatialHystereticTerminalCircuit, SpatialHystereticTerminalCircuitError,
};
use fs_phs::circuit::{Branch, CircuitGraph};

const CAPACITANCE: f64 = 2.0;
const INITIAL_CHARGE: f64 = 56.0;

fn capacitor<'m>(
    problem: &'m fs_em::spatial::vector_potential::induction::hysteretic::SpatialHystereticProblem,
    source: fs_material::conductor::ResolvedConductorState,
    reference: f64,
) -> SpatialHystereticTerminalCircuit<'m> {
    let field = SpatialHystereticTerminalInduction::new(
        data::initial(problem, Some(source)),
        &electrodes(problem),
        preparation(),
        &mut || false,
    )
    .unwrap();
    SpatialHystereticTerminalCircuit::new(
        field,
        CircuitGraph {
            node_count: 2,
            branches: vec![(
                1,
                0,
                Branch::Capacitor {
                    farads: CAPACITANCE,
                },
            )],
            transformers: vec![],
        },
        vec![1, 0],
        vec![],
        initial(vec![], vec![INITIAL_CHARGE], vec![]),
        config(reference),
        &mut || false,
    )
    .unwrap()
}

#[test]
fn capacitor_discharge_retains_vector_memory_and_cancels_reference_loss() {
    let card = data::card(500.0);
    let (positions, tets) = grid(1);
    let problem = data::problem(&card, positions, tets);
    let edge = problem.mesh().edges.binary_search(&[0, 7]).unwrap();
    let sigma = resistivity(300.0).recip();
    for reference in [0.01, 1.0, 100.0] {
        let mut circuit = capacitor(&problem, conductor(&card, 300.0), reference);
        let mut old = Endpoint::default();
        let mut charge = INITIAL_CHARGE;
        let (mut dissipated, mut numerical_loss, mut pinning_total) = (0.0, 0.0, 0.0);
        for _ in 0..14 {
            // Native capacitor: q_next=q_old-dt*I and U=(q_old+q_next)/(2C).
            // The field retains backward-Euler H and the original scalar play.
            let h = bracket(-29.0, 29.0, |h| {
                let next = at_field(old, sigma, DT, h);
                next.voltage + DT * next.incoming / (2.0 * CAPACITANCE) - charge / CAPACITANCE
            });
            let next = at_field(old, sigma, DT, h);
            let next_charge = charge - DT * next.incoming;
            let heat = joule(old, next, sigma, DT);
            let pinning = PIN * (next.z - old.z).abs();
            let change = stored(next.h, next.z) - stored(old.h, old.z);
            let numerical = next.h * (next.b - old.b) - change - pinning;
            let proposal = circuit
                .step_candidate(
                    &[],
                    &[[0.0; 3]; 6],
                    DT,
                    controls(field_controls()),
                    &mut || false,
                )
                .unwrap();
            let step = proposal.step();
            assert_memory(proposal.circuit().field().cell_states(), next);
            near(
                proposal.circuit().field().potential_wb()[edge],
                0.5 * next.b,
                5e-8,
            );
            near(step.terminal_voltage_v[0], next.voltage, 5e-8);
            near(step.terminal_current_a[0], next.incoming, 5e-7);
            near(
                proposal.circuit().native_state()[proposal.circuit().dae().charge_index[0]],
                next_charge,
                5e-8,
            );
            near(
                step.circuit_stored_change_j,
                next.voltage * (next_charge - charge),
                5e-7,
            );
            near(
                step.circuit_stored_after_j,
                next_charge.powi(2) / (2.0 * CAPACITANCE),
                5e-7,
            );
            near(step.field.joule_heat_j, heat, 5e-7);
            near(step.field.hysteresis_heat_j, pinning, 5e-8);
            near(step.field.stored_change_j, change, 5e-8);
            near(step.numerical_dissipation_j, numerical, 5e-8);
            near(step.heat_j, heat + pinning, 5e-7);
            near(step.resistor_heat_j, 0.0, 1e-14);
            near(step.external_work_j, 0.0, 1e-14);
            near(step.energy_residual_j, 0.0, 2e-8);
            assert!(step.graph_resistors.is_empty());
            assert!(step.cancelled_reference_heat_j > 0.0);
            dissipated += heat + pinning;
            numerical_loss += numerical;
            pinning_total += pinning;
            near(
                step.stored_after_j + dissipated + numerical_loss,
                INITIAL_CHARGE.powi(2) / (2.0 * CAPACITANCE),
                8e-7,
            );
            old = next;
            charge = next_charge;
            circuit = proposal.into_parts().0;
        }
        // The discharged circuit leaves a real remanent core; resetting the
        // return-map memory on each circuit probe cannot reproduce this path.
        assert!(old.z > 0.3 && old.z < 0.5);
        assert!(old.h.abs() < 0.05 && charge.abs() < 0.02);
        assert!(pinning_total > 3.0 && numerical_loss > 1.0);
    }
}

#[test]
fn late_cancellation_and_trial_exhaustion_preserve_native_storage_and_vector_history() {
    let card = data::card(500.0);
    let (positions, tets) = grid(1);
    let problem = data::problem(&card, positions, tets);
    let mut circuit = capacitor(&problem, conductor(&card, 300.0), 0.7);
    for _ in 0..2 {
        circuit = circuit
            .step_candidate(
                &[],
                &[[0.0; 3]; 6],
                DT,
                controls(field_controls()),
                &mut || false,
            )
            .unwrap()
            .into_parts()
            .0;
    }
    let before = (
        circuit.native_state().to_vec(),
        circuit.field().potential_wb().to_vec(),
        circuit.field().cell_states().to_vec(),
        circuit.time_s(),
    );
    let mut polls = 0;
    let expected = circuit
        .step_candidate(
            &[],
            &[[0.0; 3]; 6],
            DT,
            controls(field_controls()),
            &mut || {
                polls += 1;
                false
            },
        )
        .unwrap();
    assert!(polls > 10);
    let mut reached = 0;
    assert!(matches!(
        circuit.step_candidate(
            &[],
            &[[0.0; 3]; 6],
            DT,
            controls(field_controls()),
            &mut || {
                reached += 1;
                reached >= polls - 1
            },
        ),
        Err(SpatialHystereticTerminalCircuitError::Cancelled)
    ));
    let mut limited = controls(field_controls());
    limited.maximum_field_trials = 1;
    assert!(matches!(
        circuit.step_candidate(&[], &[[0.0; 3]; 6], DT, limited, &mut || false,),
        Err(SpatialHystereticTerminalCircuitError::FieldTrialBudget)
    ));
    assert_eq!(circuit.native_state(), before.0);
    assert_eq!(circuit.field().potential_wb(), before.1);
    assert_eq!(circuit.field().cell_states(), before.2);
    assert_eq!(circuit.time_s(), before.3);
    let retry = circuit
        .step_candidate(
            &[],
            &[[0.0; 3]; 6],
            DT,
            controls(field_controls()),
            &mut || false,
        )
        .unwrap();
    assert_eq!(
        retry.circuit().native_state(),
        expected.circuit().native_state()
    );
    assert_eq!(
        retry.circuit().field().potential_wb(),
        expected.circuit().field().potential_wb()
    );
    assert_eq!(
        retry.circuit().field().cell_states(),
        expected.circuit().field().cell_states()
    );
    assert_eq!(retry.step().heat_j, expected.step().heat_j);
}
