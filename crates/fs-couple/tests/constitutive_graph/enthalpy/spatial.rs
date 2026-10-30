//! Native spatial balances, source-data substitutions and transactional histories.
//! The unit-tetra reference below uses explicit P1 integrals, not the assembler.
use super::*;
use fs_conduction::transient::enthalpy::heterogeneous::{
    EnthalpyRadiation, HeterogeneousEnthalpyBackwardEuler, ReferenceEnthalpyMaterial,
};
use fs_conduction::transient::enthalpy::{EnthalpyBudget, EnthalpyStepConfig};
use fs_conduction::{
    ConductionMesh, ConductionProblem, ConductivityModel, ConductivityTable, ScalarField,
    ThermalBc, ThermalBoundary, ThermalBoundaryBuilder,
};
use fs_couple::constitutive_graph::spatial::{
    self as spatial_coupling, GraphSpatialError, GraphSpatialHeatPort, GraphSpatialInput,
    GraphSpatialSolve, GraphSpatialTemperature,
};
use fs_rep_mesh::TetComplex;

fn mesh(two_materials: bool) -> ConductionMesh {
    let mut positions = vec![
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 0.0, 1.0],
    ];
    let mut tets = vec![[0, 1, 2, 3]];
    if two_materials {
        positions.extend([
            [2.0, 0.0, 0.0],
            [3.0, 0.0, 0.0],
            [2.0, 1.0, 0.0],
            [2.0, 0.0, 1.0],
        ]);
        tets.push([4, 5, 6, 7]);
    }
    ConductionMesh::new(TetComplex::from_tets(positions.len(), tets), positions).unwrap()
}

fn spatial_controls() -> GraphSpatialSolve {
    GraphSpatialSolve {
        maximum_vertices: 8,
        maximum_binding_terms: 64,
        maximum_endpoint_trials: 500,
        maximum_graph_evaluations: 5000,
        difference_step_j: 1e-5,
        node_tolerance_j: 1e-8,
        thermal: EnthalpyStepConfig {
            newton: NewtonKrylovConfig {
                absolute_tolerance: 1e-11,
                relative_tolerance: 1e-12,
                linear_restart: 8,
                max_linear_cycles: 8,
                forcing_minimum: 1e-10,
                forcing_maximum: 1e-4,
                ..NewtonKrylovConfig::default()
            },
            max_newton_iterations: 32,
            energy_tolerance_j: 1e-8,
        },
        graph: controls(),
    }
}

fn prepare<'m, 'c>(
    cx: &Cx<'_>,
    mesh: &'m ConductionMesh,
    materials: &[ReferenceEnthalpyMaterial<'c>],
    ids: &[usize],
) -> HeterogeneousEnthalpyBackwardEuler<'m, 'c> {
    HeterogeneousEnthalpyBackwardEuler::new(
        cx,
        mesh,
        materials,
        ids,
        EnthalpyBudget {
            max_vertices: 8,
            max_elements: 2,
        },
    )
    .unwrap()
}

fn conductivity(card: &MaterialCard) -> ConductivityModel {
    ConductivityModel::isotropic(
        ConductivityTable::from_claims_at_query_points(
            card.claims(),
            &PropertyKey::new("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS),
            "T",
            &[point(200.0), point(800.0)],
            SelectionPolicy::SingleClaimOnly,
        )
        .unwrap(),
    )
}

fn problem<'a>(
    mesh: &'a ConductionMesh,
    boundary: &'a ThermalBoundary,
    conductivity: &'a ConductivityModel,
    source: &'a ScalarField,
) -> ConductionProblem<'a> {
    ConductionProblem {
        mesh,
        boundary,
        material: conductivity,
        element_materials: None,
        source,
    }
}

fn input<'a>(
    graph: &'a ConstitutiveGraph,
    history: &'a [f64],
    feedback: &'a Feedback,
    temperatures: &'a [GraphSpatialTemperature],
    heat: &'a [GraphSpatialHeatPort],
) -> GraphSpatialInput<'a> {
    GraphSpatialInput {
        graph,
        version: graph.state_schema().version(),
        state: history,
        fixed: &feedback.fixed,
        unknowns: &feedback.unknowns,
        equations: &feedback.equations,
        initial_coordinates: if feedback.unknowns.is_empty() {
            &[]
        } else {
            &[0.0]
        },
        temperatures,
        durations: &feedback.durations,
        heat,
    }
}

fn resistor_sample(vertices: Vec<(usize, f64)>) -> [GraphSpatialTemperature; 1] {
    [GraphSpatialTemperature {
        input: key("resistor", "temperature_change"),
        vertices,
        coordinate: ThermalCoordinate::ChangeFrom(300.0),
    }]
}

// V=1/6, mass=1/24 per vertex, grad N0=(-1,-1,-1), grad Nj=e_j.
// The normalized P1 shape N0 supplies loads P*(2,1,1,1)/5, NOT P*(1,0,0,0).
// Conservation fixes mean T at each scalar power trial, hence k(mean T).
// A star-matrix elimination then solves all four temperatures; bisection closes
// P=V²/R(T_sample). Neither native transport assembly nor Newton is replayed.
fn reference(old: [f64; 4], dt: f64, r0: f64, cp: f64, k_scale: f64) -> ([f64; 4], f64) {
    let capacity = cp / 24.0;
    let at_power = |power: f64| {
        let mean = old.iter().sum::<f64>() / 4.0 + dt * power / (4.0 * capacity);
        let a = dt * k_scale * (1.0 + 0.002 * (mean - 300.0)) / 6.0;
        let rhs: [f64; 4] = std::array::from_fn(|i| {
            capacity * old[i] + dt * power * if i == 0 { 0.4 } else { 0.2 }
        });
        let first = (rhs[0] + a * (rhs[1] + rhs[2] + rhs[3]) / (capacity + a))
            / (capacity + 3.0 * a - 3.0 * a * a / (capacity + a));
        [
            first,
            (rhs[1] + a * first) / (capacity + a),
            (rhs[2] + a * first) / (capacity + a),
            (rhs[3] + a * first) / (capacity + a),
        ]
    };
    let mut lower = 0.0;
    let mut upper = VOLTAGE * VOLTAGE / (r0 - 100.0 * RESISTANCE_SLOPE);
    for _ in 0..100 {
        let power = (lower + upper) * 0.5;
        let t = at_power(power);
        let resistance = r0 + RESISTANCE_SLOPE * (0.75 * t[0] + 0.25 * t[1] - 300.0);
        if power * resistance < VOLTAGE * VOLTAGE {
            lower = power;
        } else {
            upper = power;
        }
    }
    let power = (lower + upper) * 0.5;
    (at_power(power), power)
}

fn sourced_trial(cx: &Cx<'_>, r0: f64, cp: f64, k_scale: f64) -> (Vec<f64>, f64) {
    let mesh = mesh(false);
    let (card, caloric) = material(r0, cp, k_scale);
    let storage = prepare(
        cx,
        &mesh,
        &[ReferenceEnthalpyMaterial {
            curve: caloric.curve(),
            reference_density_kg_m3: 1.0,
        }],
        &[0],
    );
    let boundary = ThermalBoundaryBuilder::new(&mesh)
        .adiabatic_remainder()
        .finish()
        .unwrap();
    let conductivity = conductivity(&card);
    assert!(!conductivity.receipts().is_empty());
    let source = ScalarField::Uniform(0.0);
    let graph = ohmic_graph(&card);
    let feedback = voltage_feedback();
    let samples = resistor_sample(vec![(0, 0.75), (1, 0.25)]);
    let heat = [GraphSpatialHeatPort {
        output: key("resistor", "heat"),
        distribution: ScalarField::Nodal(vec![1.0, 0.0, 0.0, 0.0]),
    }];
    let old_t = [340.0, 300.0, 310.0, 320.0];
    let old_h = old_t.map(|t| cp * (t - 200.0));
    let dt = 0.03;
    let candidate = spatial_coupling::step_candidate(
        cx,
        &storage,
        problem(&mesh, &boundary, &conductivity, &source),
        None,
        &old_h,
        dt,
        EnthalpyRadiation::None,
        input(&graph, &[], &feedback, &samples, &heat),
        spatial_controls(),
        &[0.0; 4],
    )
    .unwrap();
    let report = candidate.report();
    let t = &report.thermal.temperature;
    let (expected, power) = reference(old_t, dt, r0, cp, k_scale);
    for i in 0..4 {
        close(t[i], expected[i]);
        close(
            candidate.specific_enthalpy_j_kg()[i],
            cp * (expected[i] - 200.0),
        );
        close(report.thermal.residual_j[i], 0.0);
        close(storage.reference_nodal_masses_kg()[i], 1.0 / 24.0);
    }
    let sampled = 0.75 * t[0] + 0.25 * t[1];
    let current = candidate.graph().external[&key("resistor", "current")];
    close(
        candidate.graph().external[&key("resistor", "temperature_change")],
        sampled - 300.0,
    );
    close(current, power / VOLTAGE);
    close(report.graph_heat_j, dt * power);
    close(report.thermal.stored_energy_change_j, dt * power);
    close(report.other_external_heat_j, 0.0);
    assert!(report.endpoint_trials > 1 && report.maximum_node_residual_j < 1e-8);
    assert_eq!(old_h, old_t.map(|t| cp * (t - 200.0)));
    (t.clone(), current)
}

#[test]
fn spatial_material_feedback_matches_independent_p1_storage_and_source_integration() {
    with_cx(|cx, _| {
        let ordinary = sourced_trial(cx, 2.0, 2.0, 0.2);
        let resistive = sourced_trial(cx, 4.0, 2.0, 0.2);
        let capacious = sourced_trial(cx, 2.0, 4.0, 0.2);
        assert!(
            ordinary
                .0
                .windows(2)
                .any(|pair| (pair[0] - pair[1]).abs() > 1.0)
        );
        assert!(resistive.0.iter().zip(&ordinary.0).all(|(a, b)| a < b));
        assert!(resistive.1 < ordinary.1);
        assert!(capacious.0[0] < ordinary.0[0] && capacious.1 > ordinary.1);
    });
}

#[test]
fn spatial_conductivity_data_changes_native_temperature_field_and_electric_current() {
    with_cx(|cx, _| {
        let low = sourced_trial(cx, 2.0, 2.0, 0.2);
        let high = sourced_trial(cx, 2.0, 2.0, 3.0);
        assert!(high.0[0] < low.0[0]);
        assert!(high.0[1..].iter().zip(&low.0[1..]).all(|(a, b)| a > b));
        assert!(high.1 > low.1);
    });
}

#[test]
fn spatial_heterogeneous_latent_progress_preserves_charts_and_commits_graph_memory_once() {
    with_cx(|cx, _| {
        let mesh = mesh(true);
        let latent = latent_curve();
        let (card, sensible) = material(2.0, 2.0, 1.0);
        let storage = prepare(
            cx,
            &mesh,
            &[
                ReferenceEnthalpyMaterial {
                    curve: &latent,
                    reference_density_kg_m3: 2.0,
                },
                ReferenceEnthalpyMaterial {
                    curve: sensible.curve(),
                    reference_density_kg_m3: 4.0,
                },
            ],
            &[0, 1],
        );
        let boundary = ThermalBoundaryBuilder::new(&mesh)
            .adiabatic_remainder()
            .finish()
            .unwrap();
        let conductivity = conductivity(&card);
        let source = ScalarField::Uniform(0.0);
        let mut graph = ohmic_graph(&card);
        let observed = Rc::new(RefCell::new(Vec::new()));
        graph
            .add_node(
                "dose",
                Box::new(HeatDose {
                    declaration: declaration(
                        "spatial-heat-dose",
                        &[("power", W), ("dt", Dims([0, 0, 1, 0, 0, 0]))],
                        &[("heat", W)],
                        &["dose_j"],
                    ),
                    observed_history: observed.clone(),
                }),
            )
            .unwrap();
        graph
            .connect(Edge {
                from: "resistor".into(),
                from_port: "heat".into(),
                to: "dose".into(),
                to_port: "power".into(),
            })
            .unwrap();
        observed.borrow_mut().clear();
        let mut feedback = voltage_feedback();
        feedback.durations = vec![key("dose", "dt")];
        let samples = resistor_sample(vec![(0, 0.4), (4, 0.6)]);
        let heat = [GraphSpatialHeatPort {
            output: key("dose", "heat"),
            distribution: ScalarField::Nodal(vec![1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0]),
        }];
        let old = [10.0, 10.0, 10.0, 10.0, 400.0, 400.0, 400.0, 400.0];
        let history = [5.0];
        let first = spatial_coupling::step_candidate(
            cx,
            &storage,
            problem(&mesh, &boundary, &conductivity, &source),
            None,
            &old,
            0.01,
            EnthalpyRadiation::None,
            input(&graph, &history, &feedback, &samples, &heat),
            spatial_controls(),
            &[0.0; 8],
        )
        .unwrap();
        for i in 0..8 {
            close(first.report().thermal.temperature[i], 400.0);
            close(
                first.specific_enthalpy_j_kg()[i],
                if i < 4 { 11.0 } else { 400.0 },
            );
            close(
                first.report().liquid_mass_fraction[i],
                if i < 4 { 0.11 } else { 0.0 },
            );
        }
        close(
            first.graph().external[&key("resistor", "temperature_change")],
            100.0,
        );
        close(first.graph().graph.next_state[0], 5.0 + 1.0 / 3.0);
        assert!(observed.borrow().len() > 1 && observed.borrow().iter().all(|&h| h == history[0]));
        assert!(std::ptr::eq(storage.phase_curve_at(0).unwrap(), &latent));
        assert!(std::ptr::eq(
            storage.phase_curve_at(4).unwrap(),
            sensible.curve()
        ));
        let (accepted_h, accepted_graph, _) = first.into_parts();
        let saved_h = accepted_h.clone();
        let saved_memory = accepted_graph.graph.next_state.clone();
        observed.borrow_mut().clear();
        let second = spatial_coupling::step_candidate(
            cx,
            &storage,
            problem(&mesh, &boundary, &conductivity, &source),
            None,
            &accepted_h,
            0.02,
            EnthalpyRadiation::None,
            input(&graph, &saved_memory, &feedback, &samples, &heat),
            spatial_controls(),
            &[0.0; 8],
        )
        .unwrap();
        for i in 0..4 {
            close(second.specific_enthalpy_j_kg()[i], 13.0);
            close(second.report().thermal.temperature[i], 400.0);
            close(second.report().liquid_mass_fraction[i], 0.13);
        }
        close(second.graph().graph.next_state[0], 6.0);
        close(second.graph().external[&key("dose", "dt")], 0.02);
        close(second.report().thermal.stored_energy_change_j, 2.0 / 3.0);
        assert!(
            observed.borrow().len() > 1 && observed.borrow().iter().all(|&h| h == saved_memory[0])
        );
        assert_eq!(accepted_h, saved_h);
        assert_eq!(accepted_graph.graph.next_state, saved_memory);
        assert_eq!(history, [5.0]);
        assert_eq!(old, [10.0, 10.0, 10.0, 10.0, 400.0, 400.0, 400.0, 400.0]);
    });
}

#[test]
fn spatial_prescribed_field_keeps_bath_reactions_and_admits_original_graph_sources() {
    with_cx(|cx, _| {
        let mesh = mesh(false);
        let (card, caloric) = material(2.0, 2.0, 1.0);
        let storage = prepare(
            cx,
            &mesh,
            &[ReferenceEnthalpyMaterial {
                curve: caloric.curve(),
                reference_density_kg_m3: 1.0,
            }],
            &[0],
        );
        let bath = [350.0, 340.0, 330.0, 320.0];
        let boundary = ThermalBoundaryBuilder::new(&mesh)
            .remainder(
                "finite source / prescribed bath",
                ThermalBc::Dirichlet {
                    temperature: ScalarField::Nodal(bath.to_vec()),
                },
            )
            .unwrap()
            .finish()
            .unwrap();
        let conductivity = conductivity(&card);
        let source = ScalarField::Uniform(6.0); // One watt on the unit tetrahedron.
        let graph = ohmic_graph(&card);
        let feedback = voltage_feedback();
        let samples = resistor_sample(vec![(0, 0.75), (1, 0.25)]);
        let heat = [GraphSpatialHeatPort {
            output: key("resistor", "heat"),
            distribution: ScalarField::Uniform(1.0),
        }];
        let old = [200.0, 220.0, 240.0, 260.0];
        let candidate = spatial_coupling::step_candidate(
            cx,
            &storage,
            problem(&mesh, &boundary, &conductivity, &source),
            None,
            &old,
            0.1,
            EnthalpyRadiation::None,
            input(&graph, &[], &feedback, &samples, &heat),
            spatial_controls(),
            &[],
        )
        .unwrap();
        let report = candidate.report();
        let power = VOLTAGE * VOLTAGE / (2.0 + RESISTANCE_SLOPE * 47.5);
        let stored = 80.0 / 12.0;
        for (&temperature, prescribed) in report.thermal.temperature.iter().zip(bath) {
            close(temperature, prescribed);
        }
        assert_eq!(report.thermal.free_dofs, 0);
        assert_eq!(report.endpoint_trials, 1);
        close(report.thermal.stored_energy_change_j, stored);
        close(report.graph_heat_j, 0.1 * power);
        close(report.thermal.source_w, 1.0 + power);
        close(report.thermal.dirichlet_in_w, stored / 0.1 - power - 1.0);
        close(report.other_external_heat_j, stored - 0.1 * power);
        close(report.thermal.energy_residual_j, 0.0);
        let mut wrong_history = input(&graph, &[], &feedback, &samples, &heat);
        wrong_history.version += 1;
        assert!(matches!(
            spatial_coupling::step_candidate(
                cx,
                &storage,
                problem(&mesh, &boundary, &conductivity, &source),
                None,
                &old,
                0.1,
                EnthalpyRadiation::None,
                wrong_history,
                spatial_controls(),
                &[],
            ),
            Err(GraphSpatialError::Graph(_))
        ));
        let unsupported = ConductivityModel::isotropic(
            ConductivityTable::declared_curve(vec![(200.0, 1.0), (310.0, 1.0)]).unwrap(),
        );
        assert!(
            spatial_coupling::step_candidate(
                cx,
                &storage,
                problem(&mesh, &boundary, &unsupported, &source),
                None,
                &old,
                0.1,
                EnthalpyRadiation::None,
                input(&graph, &[], &feedback, &samples, &heat),
                spatial_controls(),
                &[],
            )
            .is_err()
        );
        assert_eq!(old, [200.0, 220.0, 240.0, 260.0]);
    });
}

#[test]
fn spatial_invalid_bindings_work_bounds_and_cancellation_leave_both_histories_unchanged() {
    with_cx(|cx, gate| {
        let mesh = mesh(false);
        let (card, caloric) = material(2.0, 2.0, 1.0);
        let storage = prepare(
            cx,
            &mesh,
            &[ReferenceEnthalpyMaterial {
                curve: caloric.curve(),
                reference_density_kg_m3: 1.0,
            }],
            &[0],
        );
        let boundary = ThermalBoundaryBuilder::new(&mesh)
            .adiabatic_remainder()
            .finish()
            .unwrap();
        let conductivity = conductivity(&card);
        let source = ScalarField::Uniform(0.0);
        let (graph, observed) = dose_graph();
        let feedback = dose_feedback();
        let heat = [GraphSpatialHeatPort {
            output: key("dose", "heat"),
            distribution: ScalarField::Uniform(1.0),
        }];
        let old = [200.0, 220.0, 240.0, 260.0];
        let history = [5.0];
        let run = |heat: &[GraphSpatialHeatPort], controls: GraphSpatialSolve| {
            spatial_coupling::step_candidate(
                cx,
                &storage,
                problem(&mesh, &boundary, &conductivity, &source),
                None,
                &old,
                0.001,
                EnthalpyRadiation::None,
                input(&graph, &history, &feedback, &[], heat),
                controls,
                &[0.0; 4],
            )
        };
        for shape in [
            ScalarField::Uniform(0.0),
            ScalarField::Nodal(vec![-1.0, 2.0, 2.0, 2.0]),
        ] {
            assert!(matches!(
                run(
                    &[GraphSpatialHeatPort {
                        output: key("dose", "heat"),
                        distribution: shape
                    }],
                    spatial_controls()
                ),
                Err(GraphSpatialError::Invalid(_))
            ));
        }
        assert!(matches!(
            run(&[heat[0].clone(), heat[0].clone()], spatial_controls()),
            Err(GraphSpatialError::Invalid(_))
        ));
        let mut bound = spatial_controls();
        bound.maximum_binding_terms = 0;
        assert!(matches!(
            run(&heat, bound),
            Err(GraphSpatialError::Invalid(_))
        ));
        bound = spatial_controls();
        bound.maximum_endpoint_trials = 1;
        assert!(matches!(
            run(&heat, bound),
            Err(GraphSpatialError::EndpointBudget)
        ));
        bound = spatial_controls();
        bound.maximum_graph_evaluations = 1;
        assert!(matches!(
            run(&heat, bound),
            Err(GraphSpatialError::GraphBudget
                | GraphSpatialError::Graph(MaterialGraphError::EvaluationBudget))
        ));
        assert!(observed.borrow().iter().all(|&state| state == 5.0));
        gate.request();
        assert!(matches!(
            run(&heat, spatial_controls()),
            Err(GraphSpatialError::Cancelled)
        ));
        assert_eq!(history, [5.0]);
        assert_eq!(old, [200.0, 220.0, 240.0, 260.0]);
    });
}
