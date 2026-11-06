//! Total-current solid electrodes deposit physical Joule and pinning losses.
use super::*;
#[path = "../../../../fs-em/tests/spatial_induction/terminals/hysteretic_common.rs"]
mod reference;

use fs_em::spatial::vector_potential::induction::{
    SpatialInductionLimits, hysteretic::terminals::SpatialHystereticTerminalInduction,
};
use reference::{Endpoint, assert_memory, current_endpoint, electrodes, joule};

#[test]
fn total_current_electrodes_and_retained_play_follow_finite_temperature_and_open_circuit_relaxation()
 {
    const CAPACITY: f64 = 10.0;
    let card = data::card_with_heat_capacity(500.0, CAPACITY);
    let curve = caloric(&card);
    let problem = cube(&card);
    let query = point(300.0);
    let cells = cells(&card, &curve, &query);
    let regions = regions();
    let network = SpatialMagneticThermalNetwork {
        cells: &cells,
        links: &[],
        magnetic_regions: &regions,
        volume_tolerance_m3: 1e-12,
    };
    let limits = SpatialInductionLimits {
        maximum_unknowns: 2000,
        maximum_conducting_cells: 100,
        maximum_assembly_terms: 100_000,
    };
    let mut state = SpatialHystereticTerminalInduction::new(
        data::initial(&problem, Some(conductor(&card, 300.0))),
        &electrodes(&problem),
        limits,
        &mut || false,
    )
    .unwrap()
    .into_current_driven(limits, &mut || false)
    .unwrap();
    let mut heat = thermal_with_capacity(&curve, 1.0, CAPACITY);
    let mut old = Endpoint::default();
    let mut old_t = 300.0;
    for (index, imposed) in [248.16, 0.0].into_iter().enumerate() {
        let physical = |t: f64| {
            let sigma = resistivity(t).recip();
            let next = current_endpoint(old, sigma, DT, imposed);
            joule(old, next, sigma, DT) + PIN * (next.z - old.z).abs()
        };
        let t = bracket(old_t, 480.0, |t| CAPACITY * (t - old_t) - physical(t));
        let sigma = resistivity(t).recip();
        let expected = current_endpoint(old, sigma, DT, imposed);
        let heat_joule = joule(old, expected, sigma, DT);
        let pinning = PIN * (expected.z - old.z).abs();
        let change = stored(expected.h, expected.z) - stored(old.h, old.z);
        let numerical = expected.h * (expected.b - old.b) - change - pinning;
        let before = (
            heat.cells().clone(),
            state.cell_states().to_vec(),
            state.time_s(),
        );
        let candidate = network
            .step_hysteretic_current_candidate(
                &heat,
                (index + 1) as f64 * DT,
                &state,
                &[imposed, -imposed],
                &[[0.0; 3]; 6],
                thermal_controls(),
                field_controls(),
                &[0.0],
                || false,
            )
            .unwrap();
        let report = candidate.report();
        let step = &report.coupled.electrical;
        near(
            candidate.thermal().temperature_k(CellId(7)).unwrap(),
            t,
            4e-7,
        );
        assert_memory(candidate.field().cell_states(), expected);
        near(step.terminals[0].incoming_current_a, imposed, 4e-7);
        near(step.terminals[1].incoming_current_a, -imposed, 4e-7);
        near(
            step.terminals[1].voltage_to_reference_v,
            -expected.voltage,
            4e-7,
        );
        near(step.terminal_work_j, DT * expected.voltage * imposed, 5e-7);
        near(step.joule_heat_j, heat_joule, 5e-7);
        near(step.hysteresis_heat_j, pinning, 3e-8);
        near(step.heat_j, heat_joule + pinning, 5e-7);
        near(step.stored_change_j, change, 3e-8);
        near(step.numerical_dissipation_j, numerical, 3e-8);
        near(candidate.deposited_heat_j(), heat_joule + pinning, 5e-7);
        near(
            report.coupled.thermal.stored_energy_change_j,
            heat_joule + pinning,
            5e-7,
        );
        near(
            report.total_energy_residual_j,
            0.0,
            thermal_controls().total_energy_tolerance_j,
        );
        for source in candidate.field().conductors().iter().flatten() {
            near(source.resistivity_ohm_m(), resistivity(t), 1e-12);
            near(source.material().query().axes()["T"], t, 4e-7);
        }
        if index == 0 {
            assert!(pinning > 1.0 && heat_joule > 100.0 && numerical > 1.0);
            assert!(t > 330.0 && t < 400.0);
            let cold = current_endpoint(old, resistivity(300.0).recip(), DT, imposed);
            assert!(expected.voltage > cold.voltage);
        } else {
            assert!(pinning > 0.0 && heat_joule > 0.0 && change < 0.0);
            assert!(expected.voltage.abs() > 1e-3);
            near(step.terminal_work_j, 0.0, 1e-12);
        }
        assert_eq!(heat.cells(), &before.0);
        assert_eq!(state.cell_states(), before.1);
        assert_eq!(state.time_s(), before.2);
        old = expected;
        old_t = t;
        let (next_heat, next_state, _) = candidate.into_parts();
        heat = next_heat;
        state = next_state;
    }
}
