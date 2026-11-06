//! Scaled analytical fixtures, not measured ferroelectric material fits.
#![cfg(feature = "experimental-dielectric-hysteresis")]

use fs_couple::{
    circuit_enthalpy::{
        CircuitEnthalpyError, CircuitEnthalpySolve, ResistorCell,
        contacts::{heat::MaterialHeatNetwork, material::ElectricalContactNetwork},
        material::{CellResistor, ConductionLink, MaterialCell, MaterialNetwork},
    },
    material_circuit::{
        MaterialCircuitSolve,
        dielectric_hysteretic::{
            DielectricCircuitInitial, DielectricCircuitSolve, DielectricHystereticCircuit,
            thermal::ThermalDielectricHystereticCircuit,
        },
    },
};
use fs_evidence::ValidityDomain;
use fs_exec::{CancelGate, Cx};
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, SelectionPolicy, UncertaintyModel,
};
use fs_material::{
    conductor::ElectricalPropertyBasis,
    field_storage::{
        ELECTRIC_PERMITTIVITY_DIMS,
        dielectric::{
            hysteresis::{
                UniformHystereticCapacitor,
                source::{
                    DielectricHysteresisBranchSources, DielectricHysteresisSources,
                    ResolvedDielectricHysteresis,
                },
            },
            reversible::{ELECTRIC_DISPLACEMENT_DIMS, ELECTRIC_FIELD_DIMS},
        },
    },
    phase::{HeatCapacitySourceQuery, ResolvedHeatCapacityCurve, SolidLiquidPhase},
    state_point::{
        MaterialPropertySelection, SPECIFIC_HEAT_CAPACITY_DIMS, THERMAL_CONDUCTIVITY_DIMS,
    },
    thermal::source::ThermalMaterialPoint,
};
use fs_phs::circuit::{Branch, CircuitGraph};
use fs_qty::{Density, QuantitySpec, Temperature};
use fs_solver::NewtonKrylovConfig;
use fs_thermal::enthalpy::*;
use std::sync::Arc;

fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:e} != {expected:e}, tolerance {tolerance:e}"
    );
}
fn point(temperature: f64) -> QueryPoint {
    QueryPoint::new()
        .with_quantity(
            "T",
            QuantitySpec::dimensional(Temperature::DIMS),
            temperature,
        )
        .unwrap()
}
fn roles() -> DielectricHysteresisSources {
    DielectricHysteresisSources {
        background_permittivity: PropertyKey::new(
            "switching_background",
            ELECTRIC_PERMITTIVITY_DIMS,
        ),
        field_limit: PropertyKey::new("switching_limit", ELECTRIC_FIELD_DIMS),
        branches: vec![DielectricHysteresisBranchSources {
            susceptibility: PropertyKey::new("switching_slope", ELECTRIC_PERMITTIVITY_DIMS),
            pinning_field: PropertyKey::new("switching_pinning", ELECTRIC_FIELD_DIMS),
            saturation_polarization: PropertyKey::new(
                "switching_saturation",
                ELECTRIC_DISPLACEMENT_DIMS,
            ),
        }],
    }
}
fn conductivity() -> PropertyKey {
    PropertyKey::new("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS)
}
fn material(
    alpha: f64,
    explicit_fit_temperature: bool,
) -> (MaterialCard, ResolvedHeatCapacityCurve) {
    let heat = PropertyKey::new("specific_heat_capacity", SPECIFIC_HEAT_CAPACITY_DIMS);
    let density = PropertyKey::new("density", Density::DIMS);
    let electrical = PropertyKey::new(
        "electrical_resistivity",
        ElectricalPropertyBasis::Resistivity.dimensions(),
    );
    let provenance = Provenance {
        source: "analytical dielectric/enthalpy fixture".into(),
        license: "CC0-1.0".into(),
        artifact: None,
    };
    let mut claims = ClaimSet::new();
    for (key, a, b) in [
        (heat.clone(), 1.0, 1.0),
        (density.clone(), 1.0, 1.0),
        (conductivity(), 1.0, 1.0),
        (electrical, 1.0 - 10.0 * alpha, 1.0 + 100.0 * alpha),
    ] {
        let dims = key.dims();
        claims
            .insert_claim(PropertyClaim {
                key,
                value: PropertyValue::Curve {
                    abscissa: "T".into(),
                    abscissa_dims: Temperature::DIMS,
                    knots: vec![(290.0, a), (400.0, b)],
                    dims,
                },
                validity: ValidityDomain::unconstrained().with_quantity(
                    "T",
                    QuantitySpec::dimensional(Temperature::DIMS),
                    290.0,
                    400.0,
                ),
                interpolation: InterpolationPolicy::LinearInside,
                uncertainty: UncertaintyModel::Unstated,
                observations: vec![],
                provenance: provenance.clone(),
            })
            .unwrap();
    }
    let fit = roles();
    for (key, value) in [
        (fit.background_permittivity, 0.01),
        (fit.field_limit, 30.0),
        (fit.branches[0].susceptibility.clone(), 0.2),
        (fit.branches[0].pinning_field.clone(), 2.0),
        (fit.branches[0].saturation_polarization.clone(), 1.0),
    ] {
        let dims = key.dims();
        claims
            .insert_claim(PropertyClaim {
                key,
                value: PropertyValue::Scalar { value, dims },
                validity: if explicit_fit_temperature {
                    ValidityDomain::unconstrained().with_quantity(
                        "T",
                        QuantitySpec::dimensional(Temperature::DIMS),
                        290.0,
                        330.0,
                    )
                } else {
                    ValidityDomain::unconstrained()
                },
                interpolation: InterpolationPolicy::ConstantWithinValidity,
                uncertainty: UncertaintyModel::Unstated,
                observations: vec![],
                provenance: provenance.clone(),
            })
            .unwrap();
    }
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic dielectric".into(),
            phase: "solid".into(),
            process: "scaled analytical fixture".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    let caloric = ResolvedHeatCapacityCurve::try_new(
        &card,
        HeatCapacitySourceQuery {
            heat_capacity: &heat,
            density: &density,
            temperature_axis: "T",
            lower: &point(290.0),
            upper: &point(400.0),
            selection: SelectionPolicy::SingleClaimOnly,
        },
        SolidLiquidPhase::Solid,
        0.0,
        1e-6,
    )
    .unwrap();
    (card, caloric)
}
fn capacitor(card: &MaterialCard) -> UniformHystereticCapacitor {
    UniformHystereticCapacitor::from_source(
        ResolvedDielectricHysteresis::resolve(
            Arc::new(card.clone()),
            &point(300.0),
            roles(),
            MaterialPropertySelection::SingleClaimOnly,
        )
        .unwrap(),
        1.0,
        1.0,
    )
    .unwrap()
}
fn thermal_initial(
    curve: &ResolvedHeatCapacityCurve,
    count: u64,
    linked: bool,
) -> DepositionSimulation<'_> {
    let mut events: Vec<_> = (0..count)
        .map(|i| TimedEvent {
            time_s: 0.0,
            event: Event::Birth(DepositedCell {
                id: CellId(i),
                segment_id: i,
                body_id: i,
                material_id: MaterialId(0),
                mass_kg: 1.0,
                specific_h_j_kg: 10.0,
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
    if linked {
        events.push(TimedEvent {
            time_s: 0.0,
            event: Event::Contact(ThermalLink {
                id: 7,
                a: CellId(0),
                b: CellId(1),
                area_m2: 1.0,
                conductance_w_k: 999.0,
            }),
        });
    }
    DepositionSimulation::new(
        events,
        &[(MaterialId(0), curve.curve() as &dyn ThermalLaw)],
        DepositionBudget {
            max_cells: 2,
            max_links: 1,
            max_events: 3,
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
fn electrical_controls() -> DielectricCircuitSolve {
    DielectricCircuitSolve {
        maximum_circuit_trials: 256,
        maximum_newton_iterations: 32,
        difference_step_c: 1e-6,
        terminal_current_tolerance_a: 1e-8,
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
fn heat_controls() -> CircuitEnthalpySolve {
    CircuitEnthalpySolve {
        maximum_cells: 2,
        maximum_circuit_trials: 4000,
        maximum_newton_iterations: 16,
        difference_step_j: 1e-4,
        cell_energy_tolerance_j: 1e-8,
        thermal_energy_tolerance_j: 1e-8,
        total_energy_tolerance_j: 2e-8,
        circuit: electrical_controls().circuit,
        newton: NewtonKrylovConfig {
            absolute_tolerance: 1e-10,
            relative_tolerance: 1e-13,
            linear_restart: 4,
            forcing_minimum: 1e-4,
            ..NewtonKrylovConfig::default()
        },
    }
}
fn initial(
    capacitor: &UniformHystereticCapacitor,
    voltage: f64,
    current: f64,
) -> DielectricCircuitInitial {
    DielectricCircuitInitial {
        dielectric: capacitor
            .candidate_voltage(&capacitor.initial_state(), voltage)
            .unwrap()
            .material
            .state,
        fluxes_wb: vec![],
        remaining_charges_c: vec![],
        external_inputs: vec![current],
        time_s: 0.0,
    }
}

#[test]
fn finite_two_cell_leakage_and_switching_heat_match_independent_feedback_balance() {
    let mut voltages = Vec::new();
    for alpha in [0.0, 0.05] {
        let (card, caloric) = material(alpha, true);
        let context = point(300.0);
        let conductivity = conductivity();
        let cells = [0, 1].map(|i| MaterialCell {
            cell: CellId(i),
            material_id: MaterialId(0),
            curve: caloric.curve(),
            point: ThermalMaterialPoint {
                material: &card,
                context: &context,
                temperature_axis: "T",
                liquid_fraction_axis: None,
            },
            conductivity: Some((&conductivity, MaterialPropertySelection::SingleClaimOnly)),
        });
        let resistors = [CellResistor {
            destination: ResistorCell {
                branch: 1,
                cell: CellId(0),
            },
            selection: MaterialPropertySelection::SingleClaimOnly,
            length_m: 1.0,
            area_m2: 1.0,
        }];
        let links = [ConductionLink {
            link: 7,
            distance_a_m: 0.5,
            distance_b_m: 0.5,
        }];
        let network = MaterialNetwork {
            cells: &cells,
            resistors: &resistors,
            links: &links,
        };
        let thermal = thermal_initial(&caloric, 2, true);
        let capacitor = capacitor(&card);
        let start = initial(&capacitor, 3.0, 4.0);
        let graph = CircuitGraph {
            node_count: 2,
            branches: vec![
                (1, 0, Branch::CurrentSource { port: 0 }),
                (1, 0, Branch::Resistor { ohms: 999.0 }),
                (1, 0, Branch::Capacitor { farads: 999.0 }),
            ],
            transformers: vec![],
        };
        let native = DielectricHystereticCircuit::new(
            graph,
            2,
            capacitor,
            network.resolve(&thermal, || false).unwrap().resistors,
            start,
        )
        .unwrap();
        let coupled = ThermalDielectricHystereticCircuit::bind(
            native,
            &thermal,
            &network,
            &[(1, ElectricalPropertyBasis::Resistivity)],
            CellId(0),
            [290.0, 330.0],
            || false,
        )
        .unwrap();
        let original_circuit = coupled.circuit().clone();
        let original_cells = thermal.cells().clone();
        let dt = 0.02;
        let bases = [(1, ElectricalPropertyBasis::Resistivity)];
        let electrical = ElectricalContactNetwork {
            bulk: &network,
            contacts: &[],
            electrical_bases: &bases,
        };
        let heat = MaterialHeatNetwork {
            electrical: &electrical,
            thermal_contacts: &[],
            radiation: None,
        };
        let gate = CancelGate::new();
        let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
        let next = pool.scope(|arena| {
            let cx = Cx::new(
                &gate,
                arena,
                fs_exec::StreamKey {
                    seed: 17,
                    kernel_id: 1,
                    tile: 0,
                    iteration: 0,
                },
                fs_exec::Budget::INFINITE,
                fs_exec::ExecMode::Deterministic,
            );
            let prepared = heat.prepare(&cx, &thermal).unwrap();
            coupled
                .step_candidate_with_heat_assembly(
                    &cx,
                    &prepared,
                    &thermal,
                    dt,
                    &[4.0],
                    heat_controls(),
                    electrical_controls(),
                    &[0.0; 2],
                )
                .unwrap()
        });

        // Independent fixed-active-set equations: Ceff*dV/dt = I - Vmean/R(T).
        // Hswitch=Ec*chi*dV; Hohmic=dt*Vmean²/R. Two unit-capacity cells
        // exchange dt*G*(T0-T1), with G=1 from two half-unit conductive paths.
        let oracle = |theta: f64| {
            let resistance = 1.0 + alpha * theta;
            let dv = dt * (4.0 - 3.0 / resistance) / (0.21 + dt / (2.0 * resistance));
            let mean = 3.0 + 0.5 * dv;
            let switching = 0.4 * dv;
            let joule = dt * mean * mean / resistance;
            (
                dv,
                switching,
                joule,
                (switching + joule) * (1.0 + dt) / (1.0 + 2.0 * dt),
            )
        };
        let (mut low, mut high) = (0.0, 1.0);
        for _ in 0..60 {
            let theta = 0.5 * (low + high);
            if theta < oracle(theta).3 {
                low = theta;
            } else {
                high = theta;
            }
        }
        let theta = 0.5 * (low + high);
        let (dv, switching, joule, _) = oracle(theta);
        let step = &next.report().electrical;
        near(next.circuit().circuit().voltage_v(), 3.0 + dv, 2e-8);
        near(step.dielectric.switching_heat_j, switching, 2e-9);
        near(step.circuit.heat_j, joule, 2e-9);
        near(next.deposited_heat_j(), switching + joule, 2e-9);
        near(
            next.thermal().temperature_k(CellId(0)).unwrap(),
            300.0 + theta,
            2e-8,
        );
        near(
            next.thermal().temperature_k(CellId(1)).unwrap(),
            300.0 + theta * dt / (1.0 + dt),
            2e-8,
        );
        near(next.report().total_energy_residual_j, 0.0, 2e-8);
        near(step.numerical_dissipation_j, 0.0, 0.0);
        let source = next.circuit().circuit().capacitor().source().unwrap();
        assert_eq!(source.material().card_identity(), card.content_hash());
        assert_eq!(
            source.material().query(),
            point(next.thermal().temperature_k(CellId(0)).unwrap())
        );
        assert_eq!(source.law(), coupled.circuit().capacitor().law());
        assert_eq!(coupled.circuit(), &original_circuit);
        assert_eq!(thermal.cells(), &original_cells);
        near(thermal.time_s(), 0.0, 0.0);
        voltages.push(next.circuit().circuit().voltage_v());
        let (warm, accepted, _) = next.into_parts();
        let continued = accepted
            .step_candidate(
                &warm,
                2.0 * dt,
                &[4.0],
                heat_controls(),
                electrical_controls(),
                &[0.0; 2],
                || false,
            )
            .unwrap();
        assert!(continued.circuit().circuit().charge_c() > accepted.circuit().charge_c());
        assert!(
            continued.thermal().temperature_k(CellId(1)).unwrap()
                > warm.temperature_k(CellId(1)).unwrap()
        );
        near(continued.report().total_energy_residual_j, 0.0, 2e-8);
    }
    assert!(
        voltages[1] > voltages[0] + 1e-3,
        "R(T) must change the real charging trajectory"
    );
}

#[test]
fn complete_major_loop_heats_finite_dielectric_by_exact_loop_area_and_keeps_memory() {
    let (card, caloric) = material(0.0, true);
    let context = point(300.0);
    let cells = [MaterialCell {
        cell: CellId(0),
        material_id: MaterialId(0),
        curve: caloric.curve(),
        point: ThermalMaterialPoint {
            material: &card,
            context: &context,
            temperature_axis: "T",
            liquid_fraction_axis: None,
        },
        conductivity: None,
    }];
    let network = MaterialNetwork {
        cells: &cells,
        resistors: &[],
        links: &[],
    };
    let mut thermal = thermal_initial(&caloric, 1, false);
    let capacitor = capacitor(&card);
    let start = initial(&capacitor, 10.0, 0.0);
    let graph = CircuitGraph {
        node_count: 2,
        branches: vec![
            (1, 0, Branch::CurrentSource { port: 0 }),
            (1, 0, Branch::Capacitor { farads: 999.0 }),
        ],
        transformers: vec![],
    };
    let circuit = DielectricHystereticCircuit::new(graph, 1, capacitor, vec![], start).unwrap();
    let memory = circuit.dielectric_state().clone();
    let mut coupled = ThermalDielectricHystereticCircuit::bind(
        circuit,
        &thermal,
        &network,
        &[],
        CellId(0),
        [290.0, 330.0],
        || false,
    )
    .unwrap();
    let (mut deposited, mut supplied) = (0.0, 0.0);
    for (target, current) in [(1.0, -2.2), (2.0, 2.2)] {
        let next = coupled
            .step_candidate(
                &thermal,
                target,
                &[current],
                heat_controls(),
                electrical_controls(),
                &[0.0],
                || false,
            )
            .unwrap();
        near(next.report().electrical.circuit.heat_j, 0.0, 0.0);
        near(next.report().total_energy_residual_j, 0.0, 2e-8);
        deposited += next.deposited_heat_j();
        supplied += next.report().electrical.supplied_work_j;
        (thermal, coupled, _) = next.into_parts();
    }
    near(deposited, 4.0 * 2.0 * 1.0 * 1.0, 2e-8);
    near(supplied, deposited, 2e-8);
    near(thermal.temperature_k(CellId(0)).unwrap(), 308.0, 2e-8);
    near(coupled.circuit().voltage_v(), 10.0, 2e-8);
    assert_eq!(
        coupled.circuit().dielectric_state().polarizations_c_m2(),
        memory.polarizations_c_m2()
    );
    let held = coupled
        .step_candidate(
            &thermal,
            3.0,
            &[0.0],
            heat_controls(),
            electrical_controls(),
            &[0.0],
            || false,
        )
        .unwrap();
    near(held.deposited_heat_j(), 0.0, 1e-12);
    near(
        held.thermal().temperature_k(CellId(0)).unwrap(),
        308.0,
        2e-8,
    );
    assert_eq!(
        held.circuit().circuit().dielectric_state(),
        coupled.circuit().dielectric_state()
    );
}

#[test]
fn cancellation_budgets_and_temperature_refusal_publish_neither_heat_nor_memory() {
    let (card, caloric) = material(0.0, true);
    let context = point(300.0);
    let cells = [MaterialCell {
        cell: CellId(0),
        material_id: MaterialId(0),
        curve: caloric.curve(),
        point: ThermalMaterialPoint {
            material: &card,
            context: &context,
            temperature_axis: "T",
            liquid_fraction_axis: None,
        },
        conductivity: None,
    }];
    let network = MaterialNetwork {
        cells: &cells,
        resistors: &[],
        links: &[],
    };
    let thermal = thermal_initial(&caloric, 1, false);
    let cap = capacitor(&card);
    let start = initial(&cap, 3.0, 0.0);
    let graph = CircuitGraph {
        node_count: 2,
        branches: vec![
            (1, 0, Branch::CurrentSource { port: 0 }),
            (1, 0, Branch::Capacitor { farads: 999.0 }),
        ],
        transformers: vec![],
    };
    let circuit = DielectricHystereticCircuit::new(graph, 1, cap, vec![], start).unwrap();
    let coupled = ThermalDielectricHystereticCircuit::bind(
        circuit,
        &thermal,
        &network,
        &[],
        CellId(0),
        [299.0, 300.01],
        || false,
    )
    .unwrap();
    let saved = coupled.circuit().clone();
    let cells_before = thermal.cells().clone();
    let mut checks = 0;
    assert!(matches!(
        coupled.step_candidate(
            &thermal,
            0.001,
            &[0.21],
            heat_controls(),
            electrical_controls(),
            &[0.0],
            || {
                checks += 1;
                checks > 8
            },
        ),
        Err(CircuitEnthalpyError::Cancelled)
    ));
    let mut bounded = heat_controls();
    bounded.maximum_circuit_trials = 4;
    assert!(matches!(
        coupled.step_candidate(
            &thermal,
            0.001,
            &[0.21],
            bounded,
            electrical_controls(),
            &[0.0],
            || false,
        ),
        Err(CircuitEnthalpyError::CircuitBudget)
    ));
    assert!(
        coupled
            .step_candidate(
                &thermal,
                0.1,
                &[0.21],
                heat_controls(),
                electrical_controls(),
                &[0.0],
                || false,
            )
            .is_err()
    );
    assert_eq!(coupled.circuit(), &saved);
    assert_eq!(thermal.cells(), &cells_before);
    let retry = coupled
        .step_candidate(
            &thermal,
            0.001,
            &[0.21],
            heat_controls(),
            electrical_controls(),
            &[0.0],
            || false,
        )
        .unwrap();
    let repeat = coupled
        .step_candidate(
            &thermal,
            0.001,
            &[0.21],
            heat_controls(),
            electrical_controls(),
            &[0.0],
            || false,
        )
        .unwrap();
    assert_eq!(retry.circuit().circuit(), repeat.circuit().circuit());
    assert_eq!(retry.thermal().cells(), repeat.thermal().cells());
    near(retry.deposited_heat_j(), 0.0004, 2e-9);
}

#[test]
fn binding_requires_explicit_fit_temperature_and_exact_original_caloric_card() {
    for explicit in [false, true] {
        let (card, caloric) = material(0.0, explicit);
        let context = point(300.0);
        let cells = [MaterialCell {
            cell: CellId(0),
            material_id: MaterialId(0),
            curve: caloric.curve(),
            point: ThermalMaterialPoint {
                material: &card,
                context: &context,
                temperature_axis: "T",
                liquid_fraction_axis: None,
            },
            conductivity: None,
        }];
        let network = MaterialNetwork {
            cells: &cells,
            resistors: &[],
            links: &[],
        };
        let thermal = thermal_initial(&caloric, 1, false);
        // A separate alpha changes original card identity even though the
        // selected switching constants and caloric numbers are identical.
        let other = material(0.01, true).0;
        let cap = capacitor(if explicit { &other } else { &card });
        let start = initial(&cap, 0.0, 0.0);
        let graph = CircuitGraph {
            node_count: 2,
            branches: vec![
                (1, 0, Branch::CurrentSource { port: 0 }),
                (1, 0, Branch::Capacitor { farads: 1.0 }),
            ],
            transformers: vec![],
        };
        let circuit = DielectricHystereticCircuit::new(graph, 1, cap, vec![], start).unwrap();
        assert!(
            ThermalDielectricHystereticCircuit::bind(
                circuit,
                &thermal,
                &network,
                &[],
                CellId(0),
                [290.0, 330.0],
                || false,
            )
            .is_err()
        );
    }
}
