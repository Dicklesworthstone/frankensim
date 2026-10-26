//! Analytic two-species recirculation, using the real thermochemical/thermal owners.
use super::*;
use crate::thermochemical_gas::tests::model;
use crate::thermochemical_gas::energy::GasEnergySolve;
use crate::fluid_transport::FluidTransport;
use fs_thermochem::UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R;

const FLOW: f64 = 5e-5;

fn controls() -> SteadyGasSolve {
    SteadyGasSolve {
        caloric: GasEnergySolve { minimum_temperature_k: 275.0,
            maximum_temperature_k: 350.0, energy_tolerance_j_kg: 1e-8, maximum_iterations: 64 },
        maximum_dt_s: 0.05, maximum_relaxations: 2000, maximum_substeps: 100000,
        steady_power_tolerance_w: 1e-8, energy_tolerance_j: 1e-8, maximum_mean_mach: 0.05,
    }
}

fn inlet(cell: usize, temperature_k: f64, fractions: Vec<f64>) -> GasInjection {
    GasInjection { cell, mass_flow_kg_s: FLOW, state: PrimitiveGasState {
        temperature_k, composition: Composition::new(CompositionBasis::MoleFraction, fractions).unwrap(),
    } }
}

fn fixture() -> SteadyGasNetwork {
    SteadyGasNetwork {
        source: model(vec![0.75, 0.25], 0.0), source_conditions: QueryPoint::new(),
        static_pressure_pa: 101325.0, phase: ThermochemicalPhase::FrozenGas,
        cells: (0..2).map(|index| GasNetworkCell {
            geometry: CellGeometry { start_m: [index as f64 * 0.1, 0.0, 0.0],
                end_m: [(index + 1) as f64 * 0.1, 0.0, 0.0],
                volume_m3: 2e-5, surface_area_m2: 0.005 },
            flow_area_m2: 2e-4, wall_conductance_w_k: 0.0,
            wall_temperature_k: 290.0, initial_temperature_k: 310.0,
        }).collect(),
        flows: vec![
            MassFlow { id: 7, a: CellId(0), b: CellId(1), kg_s: 3.0 * FLOW },
            MassFlow { id: 11, a: CellId(1), b: CellId(0), kg_s: 2.0 * FLOW },
        ],
        inlets: vec![inlet(0, 330.0, vec![1.0, 0.0]), inlet(1, 290.0, vec![0.0, 1.0])],
        outlets: vec![GasOutlet { cell: 1, mass_flow_kg_s: 2.0 * FLOW }],
    }
}

#[test]
fn recirculation_preserves_species_and_enthalpy_and_publishes_one_physical_medium() {
    let problem = fixture();
    let solved = problem.solve(controls(), || false).unwrap();
    let cp_a = R * 3.5 / 0.028;
    let cp_b = R * 2.5 / 0.004;
    // Species flux: 3*y0 = 1+2*y1 and 4*y1 = 3*y0.
    // Energy uses each species' own enthalpy, not a temperature average.
    for (cell, a_fraction) in [(0, 2.0 / 3.0), (1, 0.5)] {
        let fractions = &solved.species().cells()[&CellId(cell as u64)];
        assert!((fractions[0] - a_fraction).abs() < 2e-13);
        let expected_t = (a_fraction * cp_a * 330.0 + (1.0 - a_fraction) * cp_b * 290.0)
            / (a_fraction * cp_a + (1.0 - a_fraction) * cp_b);
        let state = &solved.states()[cell];
        assert!((state.state().temperature - expected_t).abs() < 2e-5);
        state.state().validate_coherence().unwrap();
        let transport = state.transport_state().unwrap();
        assert_eq!(transport.density_kg_m3(), state.state().density);
        assert_eq!(transport.dynamic_viscosity_pa_s(), state.state().dynamic_viscosity);
        assert_eq!(state.state().pressure, problem.static_pressure_pa);
        assert!(solved.cell_power_residuals_w()[cell].abs() <= controls().steady_power_tolerance_w);
    }
    assert_eq!(solved.mass_flow_kg_s(), &[3.0 * FLOW, 4.0 * FLOW]);
    assert!(solved.maximum_energy_residual_j() <= controls().energy_tolerance_j);
    assert!(solved.maximum_power_residual_bound_w() <= controls().steady_power_tolerance_w);
    assert!(solved.substeps() > 0);
    // Neither branch forwards a temperature guessed from a different mixture.
    assert_ne!(solved.states()[0].state().sound_speed, solved.states()[1].state().sound_speed);
}

#[test]
fn data_only_gas_and_pressure_substitutions_reach_heat_flow_and_acoustics() {
    let mut problem = fixture();
    for cell in &mut problem.cells { cell.wall_conductance_w_k = 0.03; }
    for supply in &mut problem.inlets { supply.state.temperature_k = 330.0; }
    for supply in &mut problem.inlets {
        supply.state.composition = Composition::new(CompositionBasis::MoleFraction, vec![1.0, 0.0]).unwrap();
    }
    let a = problem.solve(controls(), || false).unwrap();
    for supply in &mut problem.inlets {
        supply.state.composition = Composition::new(CompositionBasis::MoleFraction, vec![0.0, 1.0]).unwrap();
    }
    let b = problem.solve(controls(), || false).unwrap();
    for cell in 0..2 {
        let sa = a.states()[cell].state();
        let sb = b.states()[cell].state();
        assert!(sb.temperature > sa.temperature + 3.0);
        assert!(sb.sound_speed > 2.0 * sa.sound_speed);
        assert!(sb.density < sa.density / 4.0);
        assert_ne!(sb.dynamic_viscosity, sa.dynamic_viscosity);
        assert_ne!(sb.thermal_conductivity, sa.thermal_conductivity);
    }
    problem.static_pressure_pa *= 2.0;
    let compressed = problem.solve(controls(), || false).unwrap();
    for cell in 0..2 {
        assert!((compressed.states()[cell].state().temperature - b.states()[cell].state().temperature).abs() < 2e-4);
        assert!((compressed.states()[cell].state().density / b.states()[cell].state().density - 2.0).abs() < 2e-6);
    }
    assert!(compressed.maximum_mean_mach() < 0.51 * b.maximum_mean_mach());
}

#[test]
fn connection_order_and_signed_orientation_preserve_the_same_material_states() {
    let problem = fixture();
    let original = problem.solve(controls(), || false).unwrap();
    let mut reordered = problem.clone();
    reordered.flows.reverse();
    for flow in &mut reordered.flows {
        core::mem::swap(&mut flow.a, &mut flow.b);
        flow.kg_s = -flow.kg_s;
    }
    assert_eq!(original, reordered.solve(controls(), || false).unwrap());
}

#[test]
fn network_refuses_unanchored_or_unbalanced_flow_and_unpublished_states() {
    let problem = fixture();
    let before = problem.clone();
    assert!(matches!(problem.solve(controls(), || true),
        Err(SteadyGasError::Thermal(fs_thermal::enthalpy::DepositionError::Cancelled))));
    assert!(matches!(problem.solve(SteadyGasSolve { maximum_relaxations: 1, ..controls() }, || false),
        Err(SteadyGasError::NotConverged { .. })));
    assert_eq!(problem, before);
    let mut bad_flow = problem.clone();
    bad_flow.outlets[0].mass_flow_kg_s *= 0.5;
    assert!(bad_flow.solve(controls(), || false).is_err());
    let mut sealed = problem.clone();
    sealed.inlets.clear();
    sealed.outlets.clear();
    sealed.flows[0].kg_s = 2.0 * FLOW;
    assert!(sealed.solve(controls(), || false).is_err());
    let mut source_gap = problem.clone();
    source_gap.inlets[0].state.temperature_k = 360.0;
    assert!(source_gap.solve(controls(), || false).is_err());
    let mut too_fast = problem.clone();
    too_fast.cells[0].flow_area_m2 = 1e-9;
    assert!(matches!(too_fast.solve(controls(), || false), Err(SteadyGasError::MeanFlow { .. })));
}
