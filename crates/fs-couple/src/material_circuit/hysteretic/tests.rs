use super::thermal::ThermalHystereticCircuit;
use super::*;
use crate::cavity_phs::piston::thermal::source::ThermalMaterialPoint;
use crate::circuit_enthalpy::{
    CircuitEnthalpySolve, ResistorCell,
    material::{CellResistor, MaterialCell, MaterialNetwork},
};
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, SelectionPolicy, UncertaintyModel,
};
use fs_material::{
    conductor::{ElectricalPropertyBasis, resolve_uniform_conductor},
    field_storage::hysteresis::{HysteresisBranch, ScalarHysteresis},
    phase::{HeatCapacitySourceQuery, ResolvedHeatCapacityCurve, SolidLiquidPhase},
    state_point::{MaterialPropertySelection, SPECIFIC_HEAT_CAPACITY_DIMS},
};
use fs_qty::{Density, QuantitySpec, Temperature};
use fs_thermal::enthalpy::*;

mod heat;

mod source;
mod shared;

fn near(a: f64, b: f64, tol: f64) {
    assert!((a - b).abs() <= tol, "{a:.15e} != {b:.15e}, tol {tol:e}");
}
fn point(t: f64) -> QueryPoint {
    QueryPoint::new()
        .with_quantity("T", QuantitySpec::dimensional(Temperature::DIMS), t)
        .unwrap()
}
// Deliberately scaled analytic cards, not measured magnetic/copper data.
fn material(resistance: f64, alpha: f64) -> (MaterialCard, ResolvedHeatCapacityCurve) {
    let heat = PropertyKey::new("specific_heat_capacity", SPECIFIC_HEAT_CAPACITY_DIMS);
    let density = PropertyKey::new("density", Density::DIMS);
    let electrical = PropertyKey::new(
        "electrical_resistivity",
        ElectricalPropertyBasis::Resistivity.dimensions(),
    );
    let mut claims = ClaimSet::new();
    for (key, a, b) in [
        (heat.clone(), 1.0, 1.0),
        (density.clone(), 1.0, 1.0),
        (
            electrical,
            resistance * (1.0 - 100.0 * alpha),
            resistance * (1.0 + 1700.0 * alpha),
        ),
    ] {
        let dims = key.quantity().dims();
        claims
            .insert_claim(PropertyClaim {
                key,
                value: PropertyValue::Curve {
                    abscissa: "T".into(),
                    abscissa_dims: Temperature::DIMS,
                    knots: vec![(200.0, a), (2000.0, b)],
                    dims,
                },
                validity: ValidityDomain::unconstrained().with_quantity(
                    "T",
                    QuantitySpec::dimensional(Temperature::DIMS),
                    200.0,
                    2000.0,
                ),
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::LinearInside,
                observations: vec![],
                provenance: Provenance {
                    source: "synthetic hysteretic circuit fixture".into(),
                    license: "CC0-1.0".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic".into(),
            phase: "solid".into(),
            process: "authored test".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    let curve = ResolvedHeatCapacityCurve::try_new(
        &card,
        HeatCapacitySourceQuery {
            heat_capacity: &heat,
            density: &density,
            temperature_axis: "T",
            lower: &point(200.0),
            upper: &point(2000.0),
            selection: SelectionPolicy::SingleClaimOnly,
        },
        SolidLiquidPhase::Solid,
        0.0,
        1e-6,
    )
    .unwrap();
    (card, curve)
}
fn resistor(branch: usize, card: &MaterialCard) -> MaterialResistor {
    MaterialResistor {
        branch,
        conductor: resolve_uniform_conductor(
            card,
            &point(300.0),
            MaterialPropertySelection::SingleClaimOnly,
            1.0,
            1.0,
        )
        .unwrap(),
    }
}
fn core(hysteretic: bool) -> UniformHystereticCore {
    UniformHystereticCore::new(
        ScalarHysteresis::new(
            0.01,
            30.0,
            if hysteretic {
                vec![HysteresisBranch {
                    susceptibility_h_m: 0.2,
                    pinning_field_a_m: 2.0,
                    saturation_polarization_t: 1.0,
                }]
            } else {
                vec![]
            },
        )
        .unwrap(),
        1.0,
        1.0,
        1.0,
    )
    .unwrap()
}
fn rl_graph() -> CircuitGraph {
    CircuitGraph {
        node_count: 3,
        branches: vec![
            (1, 0, Branch::VoltageSource { port: 0 }),
            (1, 2, Branch::Resistor { ohms: 999.0 }),
            (2, 0, Branch::Inductor { henries: 999.0 }),
        ],
        transformers: vec![],
    }
}
fn initial(
    core: &UniformHystereticCore,
    current: f64,
    inputs: Vec<f64>,
    charges: Vec<f64>,
) -> HystereticCircuitInitial {
    HystereticCircuitInitial {
        core: core
            .candidate_current(&core.initial_state(), current)
            .unwrap()
            .material
            .state,
        remaining_fluxes_wb: vec![],
        charges_c: charges,
        external_inputs: inputs,
        time_s: 0.0,
    }
}
fn controls() -> HystereticCircuitSolve {
    HystereticCircuitSolve {
        maximum_circuit_trials: 256,
        maximum_newton_iterations: 32,
        difference_step_wb: 1e-6,
        terminal_voltage_tolerance_v: 1e-8,
        energy_tolerance_j: 1e-8,
        circuit: MaterialCircuitSolve {
            energy_tolerance_j: 1e-9,
            power_tolerance_w: 1e-9,
        },
        newton: NewtonKrylovConfig {
            absolute_tolerance: 1e-12,
            relative_tolerance: 1e-13,
            linear_restart: 2,
            forcing_minimum: 1e-4,
            ..NewtonKrylovConfig::default()
        },
    }
}

#[test]
fn zero_pinning_linear_core_recovers_original_native_rlc_storage_and_work() {
    let (card, _) = material(0.8, 0.0);
    let core = UniformHystereticCore::new(
        ScalarHysteresis::new(0.5, 100.0, vec![]).unwrap(),
        1.0,
        1.0,
        1.0,
    )
    .unwrap();
    // Deliberately interleave capacitor and resistor before the original L.
    let graph = CircuitGraph {
        node_count: 3,
        branches: vec![
            (2, 0, Branch::Capacitor { farads: 0.2 }),
            (1, 0, Branch::VoltageSource { port: 0 }),
            (1, 2, Branch::Resistor { ohms: 0.8 }),
            (2, 0, Branch::Inductor { henries: 0.5 }),
        ],
        transformers: vec![],
    };
    let sources = vec![resistor(2, &card)];
    let mut native =
        MaterialCircuit::new(graph.clone(), sources.clone(), vec![0.0], vec![0.0], &[0.0]).unwrap();
    // Both models are initially unpowered; the held 1-V source starts at t=0.
    let start = initial(&core, 0.0, vec![0.0], vec![0.0]);
    let mut history = HystereticCircuit::new(graph, 3, core, sources.clone(), start).unwrap();
    for tick in 1..=20 {
        let a = native
            .step_candidate(sources.clone(), &[1.0], 0.005, controls().circuit, || false)
            .unwrap();
        let b = history
            .step_candidate(
                f64::from(tick) * 0.005,
                sources.clone(),
                &[1.0],
                controls(),
                || false,
            )
            .unwrap();
        near(
            b.circuit().flux_linkage_wb(),
            a.circuit().fluxes_wb()[0],
            1e-10,
        );
        near(
            b.circuit().native_circuit().charges_c()[0],
            a.circuit().charges_c()[0],
            1e-10,
        );
        near(b.step().supplied_work_j, a.step().native.supplied, 1e-10);
        near(b.step().circuit.heat_j, a.step().heat_j, 1e-10);
        near(b.step().core.core_heat_j, 0.0, 0.0);
        near(b.step().energy_residual_j, 0.0, 1e-9);
        native = a.into_parts().0;
        history = b.into_parts().0;
    }
}

#[test]
fn voltage_driven_major_loop_returns_memory_and_deposits_analytic_pinning_area() {
    let core = core(true);
    let start = initial(&core, 10.0, vec![0.0], vec![]);
    let graph = CircuitGraph {
        node_count: 2,
        branches: vec![
            (1, 0, Branch::VoltageSource { port: 0 }),
            (1, 0, Branch::Inductor { henries: 999.0 }),
        ],
        transformers: vec![],
    };
    let mut circuit = HystereticCircuit::new(graph, 1, core, vec![], start).unwrap();
    let top = circuit.magnetic_state().clone();
    let mut heat = 0.0;
    let mut supplied = 0.0;
    for (target, voltage) in [(1.0, -2.2), (2.0, 2.2)] {
        let next = circuit
            .step_candidate(target, vec![], &[voltage], controls(), || false)
            .unwrap();
        near(next.step().terminal_voltage_v, voltage, 1e-12);
        heat += next.step().core.core_heat_j;
        supplied += next.step().supplied_work_j;
        near(next.step().circuit.heat_j, 0.0, 0.0);
        near(next.step().numerical_dissipation_j, 0.0, 0.0);
        circuit = next.into_parts().0;
    }
    near(heat, 4.0 * 2.0 * 1.0 * 1.0, 1e-9);
    near(supplied, heat, 1e-9);
    near(circuit.flux_linkage_wb(), 1.1, 1e-10);
    near(circuit.current_a(), 10.0, 1e-9);
    assert_eq!(
        circuit.magnetic_state().polarizations_t(),
        top.polarizations_t()
    );
    // Demagnetizing current to zero preserves remanent flux, then zero voltage
    // holds the actual memory rather than resetting the core to a virgin state.
    let remanent = circuit
        .step_candidate(3.0, vec![], &[-0.7], controls(), || false)
        .unwrap()
        .into_parts()
        .0;
    near(remanent.current_a(), 0.0, 1e-9);
    near(remanent.flux_linkage_wb(), 0.4, 1e-10);
    let held = remanent
        .step_candidate(4.0, vec![], &[0.0], controls(), || false)
        .unwrap();
    near(held.circuit().flux_linkage_wb(), 0.4, 1e-10);
    near(held.step().core.core_heat_j, 0.0, 1e-14);
}

#[test]
fn active_hysteretic_rl_segment_converges_at_second_order_to_independent_solution() {
    let (card, _) = material(1.0, 0.0);
    // Starting at H=3 keeps the branch sliding, below its H=7 saturation.
    // The independent ODE is (mu+chi) dI/dt + R I = V.
    let expected = 4.0 - (-0.2_f64 / 0.21).exp();
    let mut prior = None;
    for count in [10, 20, 40] {
        let core = core(true);
        let start = initial(&core, 3.0, vec![4.0], vec![]);
        let sources = vec![resistor(1, &card)];
        let mut circuit =
            HystereticCircuit::new(rl_graph(), 2, core, sources.clone(), start).unwrap();
        let mut core_heat = 0.0;
        for tick in 1..=count {
            let step = circuit
                .step_candidate(
                    0.2 * f64::from(tick) / f64::from(count),
                    sources.clone(),
                    &[4.0],
                    controls(),
                    || false,
                )
                .unwrap();
            core_heat += step.step().core.core_heat_j;
            circuit = step.into_parts().0;
        }
        let error = (circuit.current_a() - expected).abs();
        if let Some(previous) = prior {
            let ratio: f64 = previous / error;
            assert!(ratio > 3.8 && ratio < 4.2, "second-order ratio {ratio}");
        }
        near(core_heat, 2.0 * 0.2 * (circuit.current_a() - 3.0), 2e-10);
        prior = Some(error);
    }
    assert!(prior.unwrap() < 2e-5);
}

#[test]
fn detectable_reversals_and_failed_trials_never_advance_memory() {
    let (card, _) = material(1.0, 0.0);
    let linear = UniformHystereticCore::new(
        ScalarHysteresis::new(1.0, 100.0, vec![]).unwrap(),
        1.0,
        1.0,
        1.0,
    )
    .unwrap();
    let start = initial(&linear, 0.0, vec![1.0], vec![]);
    let sources = vec![resistor(1, &card)];
    let circuit = HystereticCircuit::new(rl_graph(), 2, linear, sources.clone(), start).unwrap();
    let saved = circuit.clone();
    // A step longer than 2L/R overshoots the RL equilibrium in midpoint
    // stepping; endpoint V then reverses, invalidating a monotone core segment.
    assert!(matches!(
        circuit.step_candidate(4.0, sources.clone(), &[1.0], controls(), || false),
        Err(HystereticCircuitError::MagneticReversal { .. })
    ));
    assert!(matches!(
        circuit.step_candidate(0.1, sources.clone(), &[1.0], controls(), || true),
        Err(HystereticCircuitError::Cancelled)
    ));
    let mut bounded = controls();
    bounded.maximum_circuit_trials = 4;
    assert!(matches!(
        circuit.step_candidate(0.1, sources.clone(), &[1.0], bounded, || false),
        Err(HystereticCircuitError::TrialBudget)
    ));
    assert_eq!(circuit, saved);
    let mut refined = circuit;
    for tick in 1..=4 {
        refined = refined
            .step_candidate(f64::from(tick), sources.clone(), &[1.0], controls(), || {
                false
            })
            .unwrap()
            .into_parts()
            .0;
    }
    assert!(refined.current_a() > 0.98 && refined.current_a() < 1.0);
    let core = core(true);
    let start = initial(&core, 0.0, vec![0.0], vec![]);
    let graph = CircuitGraph {
        node_count: 2,
        branches: vec![
            (1, 0, Branch::VoltageSource { port: 0 }),
            (1, 0, Branch::Inductor { henries: 1.0 }),
        ],
        transformers: vec![],
    };
    let limited = HystereticCircuit::new(graph, 1, core, vec![], start).unwrap();
    let saved = limited.clone();
    assert!(matches!(
        limited.step_candidate(1.0, vec![], &[200.0], controls(), || false),
        Err(HystereticCircuitError::Material(_))
    ));
    assert_eq!(limited, saved);
}

fn thermal_initial(curves: &[ResolvedHeatCapacityCurve; 2]) -> DepositionSimulation<'_> {
    let events = (0..2)
        .map(|i| TimedEvent {
            time_s: 0.0,
            event: Event::Birth(DepositedCell {
                id: CellId(10 + i),
                segment_id: i,
                body_id: i,
                material_id: MaterialId(i),
                mass_kg: 1.0,
                specific_h_j_kg: 100.0,
                geometry: CellGeometry {
                    start_m: [i as f64, 0.0, 0.0],
                    end_m: [i as f64 + 1.0, 0.0, 0.0],
                    volume_m3: 1.0,
                    surface_area_m2: 6.0,
                },
                boundary: Boundary {
                    ambient_k: 300.0,
                    htc_w_m2_k: 0.0,
                    reservoir: None,
                    source_w: 0.0,
                },
            }),
        })
        .collect();
    DepositionSimulation::new(
        events,
        &[
            (MaterialId(0), curves[0].curve() as &dyn ThermalLaw),
            (MaterialId(1), curves[1].curve() as &dyn ThermalLaw),
        ],
        DepositionBudget {
            max_cells: 2,
            max_links: 0,
            max_events: 2,
            max_substeps: 1000,
        },
    )
    .unwrap()
    .advance_to(
        0.0,
        StepConfig {
            max_dt_s: 0.01,
            energy_tolerance_j: 1e-9,
        },
        || Ok(()),
    )
    .unwrap()
    .0
}
fn heat_controls() -> CircuitEnthalpySolve {
    CircuitEnthalpySolve {
        maximum_cells: 2,
        maximum_circuit_trials: 4000,
        maximum_newton_iterations: 16,
        difference_step_j: 1e-4,
        cell_energy_tolerance_j: 1e-8,
        thermal_energy_tolerance_j: 1e-8,
        total_energy_tolerance_j: 2e-8,
        circuit: controls().circuit,
        newton: NewtonKrylovConfig {
            absolute_tolerance: 1e-10,
            relative_tolerance: 1e-13,
            forcing_minimum: 1e-4,
            linear_restart: 4,
            ..NewtonKrylovConfig::default()
        },
    }
}
fn cells<'a>(
    cards: &'a [MaterialCard; 2],
    curves: &'a [ResolvedHeatCapacityCurve; 2],
    context: &'a QueryPoint,
) -> [MaterialCell<'a>; 2] {
    [0, 1].map(|i| MaterialCell {
        cell: CellId(10 + i as u64),
        material_id: MaterialId(i as u64),
        curve: curves[i].curve(),
        point: ThermalMaterialPoint {
            material: &cards[i],
            context,
            temperature_axis: "T",
            liquid_fraction_axis: None,
        },
        conductivity: None,
    })
}

#[test]
fn finite_core_and_coil_heat_advance_once_with_voltage_resistance_feedback() {
    let (a, ca) = material(1.0, 0.0);
    let (b, cb) = material(1.0, 0.005);
    let cards = [a, b];
    let curves = [ca, cb];
    let context = point(300.0);
    let cells = cells(&cards, &curves, &context);
    let resistors = [CellResistor {
        destination: ResistorCell {
            branch: 1,
            cell: CellId(11),
        },
        selection: MaterialPropertySelection::SingleClaimOnly,
        length_m: 1.0,
        area_m2: 1.0,
    }];
    let network = MaterialNetwork {
        cells: &cells,
        resistors: &resistors,
        links: &[],
    };
    let bases = [(1, ElectricalPropertyBasis::Resistivity)];
    let mut thermal = thermal_initial(&curves);
    let sources = network
        .prepare_with_bases(&thermal, &bases, || false)
        .unwrap()
        .resolve(&thermal, || false)
        .unwrap()
        .resistors;
    let core = core(true);
    let start = initial(&core, 0.0, vec![5.0], vec![]);
    let circuit = HystereticCircuit::new(rl_graph(), 2, core, sources, start).unwrap();
    let mut coupled = ThermalHystereticCircuit::bind(
        circuit,
        &thermal,
        &network,
        &bases,
        CellId(10),
        [290.0, 310.0],
        || false,
    )
    .unwrap();
    let initial_energy = thermal.energy_j();
    let (mut copper, mut pinning, mut supplied) = (0.0, 0.0, 0.0);
    for tick in 1..=6 {
        let old = coupled.circuit().clone();
        let saved = thermal.cells().clone();
        let next = coupled
            .step_candidate(
                &thermal,
                0.05 * f64::from(tick),
                &[5.0],
                heat_controls(),
                controls(),
                &[0.0; 2],
                || false,
            )
            .unwrap();
        let report = next.report();
        let step = &report.electrical;
        let coil_t = next.thermal().temperature_k(CellId(11)).unwrap();
        near(
            step.circuit.resistors[0]
                .material
                .conductor
                .resistance_ohm(),
            1.0 + 0.005 * (coil_t - 300.0),
            1e-12,
        );
        near(
            report.deposited_joule_j,
            step.circuit.heat_j + step.core.core_heat_j,
            1e-12,
        );
        near(report.total_energy_residual_j, 0.0, 2e-8);
        near(step.numerical_dissipation_j, 0.0, 0.0);
        copper += step.circuit.heat_j;
        pinning += step.core.core_heat_j;
        supplied += step.supplied_work_j;
        assert_eq!(coupled.circuit(), &old);
        assert_eq!(thermal.cells(), &saved);
        let (t, c, _) = next.into_parts();
        thermal = t;
        coupled = c;
    }
    assert!(copper > 0.1 && pinning > 0.1);
    near(
        thermal.temperature_k(CellId(10)).unwrap() - 300.0,
        pinning,
        2e-8,
    );
    near(
        thermal.temperature_k(CellId(11)).unwrap() - 300.0,
        copper,
        2e-8,
    );
    near(thermal.energy_j() - initial_energy, copper + pinning, 3e-8);
    near(
        coupled.circuit().core_stored_energy_j() + copper + pinning,
        supplied,
        4e-8,
    );
    assert!(
        coupled.circuit().native_circuit().resistors()[0]
            .conductor
            .resistance_ohm()
            > 1.0
    );
    assert_eq!(coupled.circuit().time_s(), thermal.time_s());
}

#[test]
fn thermal_domain_failure_preserves_both_histories_and_initial_source_binding() {
    let (a, ca) = material(1.0, 0.0);
    let (b, cb) = material(1.0, 0.005);
    let cards = [a, b];
    let curves = [ca, cb];
    let context = point(300.0);
    let cells = cells(&cards, &curves, &context);
    let resistors = [CellResistor {
        destination: ResistorCell {
            branch: 1,
            cell: CellId(11),
        },
        selection: MaterialPropertySelection::SingleClaimOnly,
        length_m: 1.0,
        area_m2: 1.0,
    }];
    let network = MaterialNetwork {
        cells: &cells,
        resistors: &resistors,
        links: &[],
    };
    let bases = [(1, ElectricalPropertyBasis::Resistivity)];
    let thermal = thermal_initial(&curves);
    let sources = network
        .prepare_with_bases(&thermal, &bases, || false)
        .unwrap()
        .resolve(&thermal, || false)
        .unwrap()
        .resistors;
    let core = core(true);
    let start = initial(&core, 0.0, vec![5.0], vec![]);
    let circuit = HystereticCircuit::new(rl_graph(), 2, core, sources, start).unwrap();
    let coupled = ThermalHystereticCircuit::bind(
        circuit,
        &thermal,
        &network,
        &bases,
        CellId(10),
        [299.0, 300.05],
        || false,
    )
    .unwrap();
    let before = coupled.circuit().clone();
    let saved = thermal.cells().clone();
    assert!(
        coupled
            .step_candidate(
                &thermal,
                0.05,
                &[5.0],
                heat_controls(),
                controls(),
                &[0.0; 2],
                || false
            )
            .is_err()
    );
    assert_eq!(coupled.circuit(), &before);
    assert_eq!(thermal.cells(), &saved);
    assert_eq!(thermal.time_s(), 0.0);
}
