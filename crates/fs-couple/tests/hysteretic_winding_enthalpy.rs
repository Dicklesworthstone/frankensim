//! Real voltage circuit, conducting core, vector memory and finite heat storage.
//! Scalar P1/play balances are independent of the production Newton/field owners.
#![cfg(feature = "magnetic-hysteresis")]

#[path = "hysteretic_induction/common.rs"]
mod common;
use common::{
    CELL, CP, DIRECTIONS, MU, PIN, VOLUME, bisect, boundary, caloric, card, cells, conductor,
    field_controls, near, paired_squares, play, point, problem, resistivity, routes, storage,
    thermal,
};
use fs_couple::circuit_enthalpy::{
    CircuitEnthalpyError, CircuitEnthalpySolve, ResistorCell, magnetic::MagneticRegion,
    material::MaterialCell, winding::WindingThermalNetwork,
};
use fs_em::{
    planar::{
        eddy::{EddyConductor, hysteretic::PlanarHystereticEddy},
        nonlinear::hysteretic::PlanarHystereticMagnetostatics,
    },
    winding::{
        hysteretic::{HystereticWindingCircuit, HystereticWindingControls},
        nonlinear::{NonlinearWindingConfig, NonlinearWindingInitial, NonlinearWindings},
        ohmic::WindingConductor,
    },
};
use fs_matdb::MaterialCard;
use fs_phs::circuit::{Branch, CircuitGraph};
use fs_solver::NewtonKrylovConfig;

const DT: f64 = 0.02;
const DEPTH: f64 = VOLUME;
const CAPACITY: f64 = CP * 2.0 * DEPTH;
const END: f64 = 0.032;
const GRAPH_R: f64 = 0.3;
const GRAPH: [ResistorCell; 1] = [ResistorCell {
    branch: 2,
    cell: CELL,
}];

fn controls() -> CircuitEnthalpySolve {
    CircuitEnthalpySolve {
        maximum_circuit_trials: 4000,
        ..common::controls()
    }
}
fn winding_controls() -> HystereticWindingControls {
    HystereticWindingControls {
        field: field_controls(),
        maximum_field_trials: 1000,
        maximum_newton_iterations: 30,
        difference_step_a: 1e-5,
        terminal_voltage_tolerance_v: 1e-8,
        power_tolerance_w: 1e-8,
        energy_tolerance_j: 2e-8,
        newton: NewtonKrylovConfig {
            absolute_tolerance: 1e-10,
            relative_tolerance: 1e-12,
            forcing_minimum: 1e-4,
            linear_restart: 8,
            ..NewtonKrylovConfig::default()
        },
    }
}
fn circuit<'a>(
    card: &MaterialCard,
    owner: &'a PlanarHystereticMagnetostatics,
) -> HystereticWindingCircuit<'a> {
    let source = conductor(card, 300.0);
    let core = [1, 7].map(|triangle| EddyConductor {
        triangle,
        region: 7,
        conductor: source.clone(),
    });
    let copper: Vec<_> = (0..8)
        .filter(|triangle| !matches!(*triangle, 1 | 7))
        .map(|triangle| WindingConductor {
            triangle,
            conductor: source.clone(),
        })
        .collect();
    let insulating: Vec<_> = copper.iter().map(|cell| cell.triangle).collect();
    let mut profile = vec![0.0; 8];
    for cell in &copper {
        profile[cell.triangle] = if cell.triangle < 4 { 8.0 } else { -8.0 };
    }
    let zero = owner.demagnetized(0.0, &mut || false).unwrap();
    let field = PlanarHystereticEddy::new(
        &zero,
        core.to_vec(),
        &insulating,
        &boundary(6),
        0.0,
        field_controls(),
        &mut || false,
    )
    .unwrap();
    HystereticWindingCircuit::new(
        field,
        CircuitGraph {
            node_count: 3,
            branches: vec![
                (1, 0, Branch::VoltageSource { port: 0 }),
                (1, 2, Branch::Inductor { henries: 999.0 }),
                (2, 0, Branch::Resistor { ohms: GRAPH_R }),
            ],
            transformers: vec![],
        },
        NonlinearWindings {
            branch_indices: vec![1],
            profiles_per_m2: vec![profile],
            end_inductance_h: vec![END],
        },
        copper,
        NonlinearWindingInitial {
            currents_a: vec![0.0],
            remaining_fluxes_wb: vec![],
            charges_c: vec![],
            external_inputs: vec![0.0],
        },
        NonlinearWindingConfig {
            max_windings: 4,
            max_profile_values: 32,
            current_balance_tolerance_turns: 1e-10,
        },
        &mut || false,
    )
    .unwrap()
}
fn layout<'a>(
    cells: &'a [MaterialCell<'a>],
    regions: &'a [MagneticRegion],
) -> WindingThermalNetwork<'a> {
    WindingThermalNetwork {
        cells,
        links: &[],
        magnetic_regions: regions,
        graph_resistors: &GRAPH,
        volume_tolerance_m3: 1e-12,
    }
}
#[derive(Clone, Copy)]
struct History {
    t: f64,
    i: f64,
    h: f64,
    z: f64,
    b: f64,
}
struct Reference {
    state: History,
    copper: f64,
    core: f64,
    pinning: f64,
    graph: f64,
    stored_before: f64,
    stored_after: f64,
    numerical: f64,
    supplied: f64,
}
impl Reference {
    fn heat(&self) -> f64 {
        self.copper + self.core + self.pinning + self.graph
    }
}
fn winding_resistance(t: f64) -> f64 {
    96.0 * DEPTH * resistivity(t)
}
// Eliminate current from the two-centre diffusion equation, then solve the
// physical voltage equation. Exact P1 moments give sigma*depth/(48*dt)*dB^2.
fn electrical(old: History, t: f64, drive: f64) -> History {
    let sigma = 1.0 / resistivity(t);
    let current = |h: f64| h + sigma / (96.0 * DT) * (MU * h + play(old.z, h) - old.b);
    let h = bisect(-29.0, 29.0, |h| {
        let i = current(h);
        (winding_resistance(t) + GRAPH_R) * i
            + (2.0 * DEPTH * (MU * h + play(old.z, h) - old.b) + END * (i - old.i)) / DT
            - drive
    });
    let z = play(old.z, h);
    History {
        t,
        i: current(h),
        h,
        z,
        b: MU * h + z,
    }
}
fn reference(old: History, drive: f64) -> Reference {
    let at = |t: f64| {
        let state = electrical(old, t, drive);
        let copper = DT * winding_resistance(t) * state.i * state.i;
        let core = DEPTH / (48.0 * DT * resistivity(t)) * (state.b - old.b).powi(2);
        let pinning = 2.0 * DEPTH * PIN * (state.z - old.z).abs();
        let graph = DT * GRAPH_R * state.i * state.i;
        (state, copper, core, pinning, graph)
    };
    let t = bisect(old.t, 500.0, |t| {
        let (_, copper, core, pinning, graph) = at(t);
        CAPACITY * (t - old.t) - copper - core - pinning - graph
    });
    let (state, copper, core, pinning, graph) = at(t);
    let stored_before = 2.0 * DEPTH * storage(old.h, old.z) + 0.5 * END * old.i * old.i;
    let stored_after = 2.0 * DEPTH * storage(state.h, state.z) + 0.5 * END * state.i * state.i;
    let numerical = 2.0
        * DEPTH
        * (state.h * (state.b - old.b) - storage(state.h, state.z) + storage(old.h, old.z))
        - pinning
        + 0.5 * END * (state.i - old.i).powi(2);
    Reference {
        state,
        copper,
        core,
        pinning,
        graph,
        stored_before,
        stored_after,
        numerical,
        supplied: DT * drive * state.i,
    }
}

#[test]
fn voltage_reversal_solves_current_core_diffusion_memory_and_temperature_together() {
    let card = card(Some(500.0));
    let curve = caloric(&card);
    let context = point(300.0);
    let owner = problem(&card, paired_squares());
    let mut state = circuit(&card, &owner);
    let mut heat = thermal(&curve, 2.0 * DEPTH);
    let cells = cells(&card, &curve, &context);
    let regions = routes(8);
    let network = layout(&cells, &regions);
    let mut old = History {
        t: 300.0,
        i: 0.0,
        h: 0.0,
        z: 0.0,
        b: 0.0,
    };
    let (mut supplied, mut numerical) = (0.0, 0.0);
    for (index, drive) in [20.0, 0.0, -30.0, 0.0].into_iter().enumerate() {
        let expected = reference(old, drive);
        let target = (index + 1) as f64 * DT;
        let next = network
            .step_hysteretic_candidate(
                &heat,
                target,
                &state,
                &[drive],
                controls(),
                winding_controls(),
                &[0.0],
                || false,
            )
            .unwrap();
        let report = next.report();
        let actual = &report.electrical;
        let e = expected.state;
        near(next.thermal().temperature_k(CELL).unwrap(), e.t);
        near(next.circuit().currents_a()[0], e.i);
        near(next.circuit().resistance_ohm()[0], winding_resistance(e.t));
        near(next.circuit().field().potential_wb_m()[6], 0.5 * e.b);
        near(next.circuit().field().potential_wb_m()[7], -0.5 * e.b);
        near(
            next.circuit().flux_linkages_wb()[0],
            2.0 * DEPTH * e.b + END * e.i,
        );
        near(actual.terminal_voltage_v[0], drive - GRAPH_R * e.i);
        near(actual.field.joule_heat_j, expected.core);
        near(actual.field.hysteresis_heat_j, expected.pinning);
        near(
            actual.cells.iter().map(|loss| loss.heat_j).sum(),
            expected.copper,
        );
        assert_eq!(actual.graph_resistors.len(), 1);
        near(actual.graph_resistors[0].heat_j, expected.graph);
        near(actual.heat_j, expected.heat());
        near(report.deposited_joule_j, expected.heat());
        near(report.thermal.stored_energy_change_j, expected.heat());
        near(actual.stored_before_j, expected.stored_before);
        near(actual.stored_after_j, expected.stored_after);
        near(actual.numerical_dissipation_j, expected.numerical);
        near(actual.supplied_work_j, expected.supplied);
        near(report.total_energy_residual_j, 0.0);
        near(report.non_electrical_heat_j, 0.0);
        assert_eq!(actual.field.cells.len(), 2);
        assert_eq!(actual.cells.len(), 6);
        for loss in &actual.field.cells {
            assert!(matches!(loss.source.triangle, 1 | 7));
            near(loss.source.conductor.resistivity_ohm_m(), resistivity(e.t));
        }
        for loss in &actual.cells {
            assert!(!matches!(loss.source.triangle, 1 | 7));
            near(loss.source.conductor.resistivity_ohm_m(), resistivity(e.t));
        }
        for (cell, memory) in next
            .circuit()
            .field()
            .magnetic_state()
            .cell_states()
            .iter()
            .enumerate()
        {
            let sign = if cell < 4 { 1.0 } else { -1.0 };
            for axis in 0..2 {
                near(
                    memory.field_a_m()[axis],
                    sign * e.h * DIRECTIONS[cell % 4][axis],
                );
                near(
                    memory.polarizations_t()[0][axis],
                    sign * e.z * DIRECTIONS[cell % 4][axis],
                );
            }
        }
        if index == 0 {
            assert!(expected.core > 0.01 && expected.copper > 0.01 && expected.pinning > 0.01);
            assert!(expected.numerical > 0.01);
            assert!(e.t > old.t + 1.0);
            assert!(electrical(old, old.t, drive).i > e.i + 1e-3);
        }
        supplied += expected.supplied;
        numerical += expected.numerical;
        near(
            actual.stored_after_j + CAPACITY * (e.t - 300.0) + numerical,
            supplied,
        );
        assert_eq!(state.time_s(), heat.time_s());
        assert_eq!(state.field().time_s(), heat.time_s());
        near(heat.temperature_k(CELL).unwrap(), old.t);
        let (next_heat, next_state, _) = next.into_parts();
        heat = next_heat;
        state = next_state;
        old = e;
    }
    // Reverse excitation leaves the actually retained vector polarization reversed.
    assert!(state.field().magnetic_state().cell_states()[0].polarizations_t()[0][0] < -0.05);
}

#[test]
fn late_cancel_and_nested_budget_preserve_circuit_field_and_finite_caloric_histories() {
    let card = card(Some(500.0));
    let curve = caloric(&card);
    let context = point(300.0);
    let owner = problem(&card, paired_squares());
    let initial = circuit(&card, &owner);
    let heat = thermal(&curve, 2.0 * DEPTH);
    let cells = cells(&card, &curve, &context);
    let regions = routes(8);
    let network = layout(&cells, &regions);
    let (heat, state, _) = network
        .step_hysteretic_candidate(
            &heat,
            DT,
            &initial,
            &[20.0],
            controls(),
            winding_controls(),
            &[0.0],
            || false,
        )
        .unwrap()
        .into_parts();
    let before = (
        heat.cells().clone(),
        state.state().to_vec(),
        state.currents_a().to_vec(),
        state.flux_linkages_wb().to_vec(),
        state.field().potential_wb_m().to_vec(),
        state.field().magnetic_state().cell_states().to_vec(),
    );
    let mut polls = 0;
    let expected = network
        .step_hysteretic_candidate(
            &heat,
            2.0 * DT,
            &state,
            &[-20.0],
            controls(),
            winding_controls(),
            &[0.0],
            || {
                polls += 1;
                false
            },
        )
        .unwrap();
    assert!(polls > 20);
    let mut limited = controls();
    limited.maximum_circuit_trials = 1;
    assert!(matches!(
        network.step_hysteretic_candidate(
            &heat,
            2.0 * DT,
            &state,
            &[-20.0],
            limited,
            winding_controls(),
            &[0.0],
            || false
        ),
        Err(CircuitEnthalpyError::CircuitBudget)
    ));
    let mut seen = 0;
    assert!(matches!(
        network.step_hysteretic_candidate(
            &heat,
            2.0 * DT,
            &state,
            &[-20.0],
            controls(),
            winding_controls(),
            &[0.0],
            || {
                seen += 1;
                seen >= polls
            }
        ),
        Err(CircuitEnthalpyError::Cancelled)
    ));
    assert_eq!(heat.cells(), &before.0);
    assert_eq!(state.state(), before.1);
    assert_eq!(state.currents_a(), before.2);
    assert_eq!(state.flux_linkages_wb(), before.3);
    assert_eq!(state.field().potential_wb_m(), before.4);
    assert_eq!(state.field().magnetic_state().cell_states(), before.5);
    assert_eq!(heat.time_s(), DT);
    assert_eq!(state.time_s(), DT);
    assert_eq!(state.field().time_s(), DT);
    let retry = network
        .step_hysteretic_candidate(
            &heat,
            2.0 * DT,
            &state,
            &[-20.0],
            controls(),
            winding_controls(),
            &[0.0],
            || false,
        )
        .unwrap();
    assert_eq!(retry.thermal().cells(), expected.thermal().cells());
    assert_eq!(retry.circuit().state(), expected.circuit().state());
    assert_eq!(
        retry.circuit().currents_a(),
        expected.circuit().currents_a()
    );
    assert_eq!(
        retry.circuit().flux_linkages_wb(),
        expected.circuit().flux_linkages_wb()
    );
    assert_eq!(
        retry.circuit().field().potential_wb_m(),
        expected.circuit().field().potential_wb_m()
    );
    assert_eq!(
        retry.circuit().field().magnetic_state().cell_states(),
        expected.circuit().field().magnetic_state().cell_states()
    );
}
