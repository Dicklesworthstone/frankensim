//! Coupled Newton must solve the physical source law, including dR/dT.
use super::*;
use super::temperature::source_point;
use fs_conduction::interface::temperature::{
    TemperatureContact, TemperatureContactCoordinates, TemperatureContacts,
};
use fs_conduction::interface::temperature::steady::SteadyContactConfig;
use fs_conduction::interface::temperature::steady::newton::{
    ContactSteadyNewtonSolution, SteadyContactNewtonConfig,
};
use fs_conduction::solve::{
    MeanTemperatureContactSolution, MeanTemperatureNewtonConfig,
    solve_mean_temperature_newton_with_contacts,
};

// Two unit slabs, k=10 W/(m K), with 300/330 K reservoirs. Each bulk
// resistance is 0.1 K/W. This positive source law has deliberately strong
// feedback, but is smooth throughout its admitted thermal neighbourhood.
fn steep_card() -> InterfaceSystemCard {
    card_with_knots(vec![(309.0, 0.04), (312.0, 0.22)])
}

fn card_with_knots(knots: Vec<(f64, f64)>) -> InterfaceSystemCard {
    use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
    let dims = AREA_SPECIFIC_THERMAL_RESISTANCE_DIMS;
    let kelvin = fs_qty::QuantitySpec::semantic(SemanticType::new(
        QuantityKind::AbsoluteTemperature, ValueForm::Static,
    ));
    let mut claims = ClaimSet::new();
    claims.insert_claim(PropertyClaim {
        key: PropertyKey::new(AREA_SPECIFIC_THERMAL_RESISTANCE_PROPERTY, dims),
        value: PropertyValue::Curve {
            abscissa: "TA".into(),
            abscissa_dims: fs_qty::Temperature::DIMS,
            knots,
            dims,
        },
        validity: ValidityDomain::unconstrained()
            .with_quantity("TA", kelvin, 309.0, 312.0)
            .with_quantity("TB", kelvin, 300.0, 340.0)
            .with_quantity("normal_pressure",
                fs_qty::QuantitySpec::dimensional(fs_qty::Pressure::DIMS), 1e6, 1e6),
        uncertainty: UncertaintyModel::Unstated,
        interpolation: InterpolationPolicy::LinearInside,
        observations: vec![],
        provenance: Provenance {
            source: "manufactured steep contact law, not measured data".into(),
            license: "internal-test-use".into(),
            artifact: None,
        },
    }).unwrap();
    let side = |name: &str| SurfaceSpec {
        material: MaterialStateId {
            chemistry: name.into(), phase: "solid".into(),
            process: "fixture".into(), revision: 0,
        },
        texture_frame: "fixed reference normal".into(),
    };
    InterfaceSystemCard::assemble(side("A"), side("B"), SystemContext {
        medium: "dry".into(), third_body: None,
        environment: "test".into(), history: "unaged".into(),
    }, claims, vec![]).unwrap()
}

fn reference_point() -> QueryPoint {
    source_point().with("TA", 311.0).unwrap().with("TB", 319.0).unwrap()
}

fn contacts<'m, 'c>(
    cx: &fs_exec::Cx<'_>, mesh: &'m ConductionMesh, bc: &ThermalBoundary,
    card: &'c InterfaceSystemCard, split: usize, matching: bool,
) -> TemperatureContacts<'m, 'c> {
    let (a, b) = sides(mesh, split, 0.0);
    if matching {
        let pairs = ThermalInterfaces::coincident_face_pairs(mesh).unwrap().into_iter()
            .map(|pair| if a.contains(&pair.side_a) { pair } else {
                InterfaceFacePair { side_a: pair.side_b, side_b: pair.side_a }
            }).collect();
        let law = TemperatureContact::new_two_sided(
            "bond", pairs, card, &reference_point(), ["TA", "TB"],
            SelectionPolicy::SingleClaimOnly,
        ).unwrap();
        TemperatureContacts::new(cx, mesh, bc, vec![], vec![law], 10_000).unwrap()
    } else {
        let law = TemperatureContact::new_nonmatching(
            "bond", a, b, card, &reference_point(),
            TemperatureContactCoordinates::TwoSided(["TA".into(), "TB".into()]),
            SelectionPolicy::SingleClaimOnly, options(),
        ).unwrap();
        TemperatureContacts::with_nonmatching(
            cx, mesh, bc, vec![], vec![], vec![law], 10_000, 10_000,
        ).unwrap()
    }
}

fn config() -> SteadyContactNewtonConfig {
    let mut result = SteadyContactNewtonConfig::default();
    // A small physical problem needs accurate directions to distinguish its
    // exact Newton response from a damped, frozen-resistance fixed point.
    result.newton.forcing_minimum = 1e-12;
    result.newton.forcing_maximum = 1e-10;
    result.newton.linear_restart = 64;
    result.newton.absolute_tolerance = 1e-9;
    result.max_iterations = 8;
    result.max_residual_evaluations = 64;
    result.max_face_pairs = 10_000;
    result.residual_tolerance_w = 1e-8;
    result.energy_tolerance_w = 1e-8;
    result
}

fn picard_solver() -> SolveConfig {
    let mut result = SolveConfig::default();
    result.nonlinearity = Nonlinearity::FixedPoint { relaxation: 1.0, max_backtracks: 16 };
    result.linear.tolerance = 1e-12;
    result.linear.max_iterations = 20_000;
    result.stop.residual_rtol = 1e-13;
    result.stop.step_atol = 0.0;
    result
}

// Independent one-dimensional bulk law for heat Q flowing from hot B to
// cold A. Inserting these interface traces into the original source gives
// Q[0.2 + R(300+0.1Q)] = 30, or .006Q^2 - .3Q - 30 = 0.
// Q=100 W is the unique root with TA in [309,312]; the other root is -50 W.
fn field(mesh: &ConductionMesh, split: usize, heat: f64) -> Vec<f64> {
    mesh.positions().iter().enumerate().map(|(i, p)| {
        if i < split { 300.0 + 0.1 * heat * p[0] }
        else { 330.0 - 0.1 * heat * (2.0 - p[0]) }
    }).collect()
}

fn check_solution(
    cx: &fs_exec::Cx<'_>, problem: ConductionProblem<'_>, split: usize,
    card: &InterfaceSystemCard, set: &TemperatureContacts<'_, '_>,
    solved: &ContactSteadyNewtonSolution,
) {
    let mesh = problem.mesh;
    let expected = field(mesh, split, 100.0);
    for (&actual, &expected) in solved.temperature.iter().zip(&expected) {
        near(actual, expected, 2e-7);
    }
    assert_eq!(solved.temperature.len(), expected.len());
    assert!(solved.newton.converged);
    assert!(solved.newton.diagnosis.is_none());
    assert!(solved.newton.iterations > 0);
    assert!(solved.newton.iterations <= config().max_iterations);
    assert_eq!(solved.newton.history.len(), solved.newton.iterations);
    assert!(solved.jacobian_assemblies > 0);
    assert!(solved.jacobian_assemblies <= solved.newton.iterations);
    assert!(solved.residual_evaluations >= solved.newton.iterations + 1);
    assert!(solved.residual_evaluations <= config().max_residual_evaluations);

    // Independently re-query the original source at the returned field and
    // assemble its physical residual. Neither the Newton stopping flag nor
    // global cancellation is an oracle for local heat balance.
    let physical = set.evaluate(cx, problem.boundary, &solved.temperature).unwrap();
    let system = assemble_operator_with_interfaces(
        cx, mesh, problem.boundary, problem.material, problem.source,
        &solved.temperature, physical.frozen_interfaces(),
    ).unwrap();
    let mut residual = vec![0.0; mesh.vertex_count()];
    system.operator.spmv(&solved.temperature, &mut residual);
    for (value, load) in residual.iter_mut().zip(&system.load) { *value -= load; }
    let dofs = DofMap::new(problem.boundary, mesh.vertex_count()).unwrap();
    assert_eq!(solved.free_residual_w.len(), dofs.free().len());
    for (&i, &(reported_i, value)) in dofs.free().iter().zip(&solved.free_residual_w) {
        assert_eq!(reported_i, i);
        near(value, residual[i], 2e-10);
        assert!(residual[i].abs() <= config().residual_tolerance_w);
    }
    let norm = dofs.free().iter().map(|&i| residual[i] * residual[i]).sum::<f64>().sqrt();
    near(solved.residual_norm_w, norm, 2e-10);
    assert!(norm <= config().residual_tolerance_w);
    near(solved.energy.closure_w, 0.0, config().energy_tolerance_w);
    near(solved.energy.source_w, 0.0, 1e-12);
    near(solved.energy.neumann_out_w, 0.0, 1e-12);
    assert!(solved.robin_fluxes.is_empty());
    for (x, expected) in [(0.0, -100.0), (2.0, 100.0)] {
        let heat = solved.dirichlet_reactions_w.iter()
            .filter(|(i, _)| on_box_face(mesh.positions()[*i][0], x))
            .map(|(_, heat)| heat).sum();
        near(heat, expected, 2e-6);
    }
    let flux = solved.physical_contacts.frozen_interfaces()
        .fluxes(&solved.temperature).unwrap();
    assert_eq!(flux.len(), 1);
    near(flux[0].heat_rate_a_to_b_w, -100.0, 2e-6);
    near(flux[0].area_m2, 1.0, 2e-12);

    // Internal heat is equal and opposite on the original independent traces,
    // down to each P1 vertex's independently integrated area/3 load.
    let linear = set.linearize_temperature(
        cx, problem.boundary, &solved.temperature, config().max_face_pairs,
    ).unwrap();
    let mut expected_load = vec![0.0; mesh.vertex_count()];
    let (a, b) = sides(mesh, split, 0.0);
    for (faces, flux) in [(a, -100.0), (b, 100.0)] {
        for slot in faces {
            let face = &mesh.boundary()[slot];
            for vertex in face.vertices {
                expected_load[vertex as usize] += flux * face.area / 3.0;
            }
        }
    }
    for (&actual, &expected) in linear.nodal_outflow_w().iter().zip(&expected_load) {
        near(actual, expected, 2e-6);
    }
    near(linear.nodal_outflow_w()[..split].iter().sum(), -100.0, 2e-6);
    near(linear.nodal_outflow_w()[split..].iter().sum(), 100.0, 2e-6);
    near(linear.nodal_outflow_w().iter().sum(), 0.0, 2e-10);

    // The accepted source evidence retains the actual local temperatures,
    // original ordered card, selected claim, and fixed nonthermal condition.
    let pinned = InterfaceResistance::from_card(
        "bond", card, &reference_point(), SelectionPolicy::SingleClaimOnly,
    ).unwrap();
    let states = solved.physical_contacts.states().iter()
        .map(|state| (state.side_a_temperature_k, state.side_b_temperature_k, &state.resistance))
        .chain(solved.physical_contacts.quadrature_states().iter()
            .map(|state| (state.side_a_temperature_k, state.side_b_temperature_k, &state.resistance)));
    let mut count = 0;
    for (ta, tb, resistance) in states {
        count += 1;
        near(ta, 310.0, 2e-7);
        near(tb, 320.0, 2e-7);
        near(resistance.value_m2_k_per_w(), 0.1, 2e-8);
        assert_eq!(resistance.card_identity(), card.content_hash());
        assert_eq!(resistance.receipt().selected, pinned.receipt().selected);
        for (axis, value) in [("TA", ta), ("TB", tb), ("normal_pressure", 1e6)] {
            assert!(resistance.receipt().query_point.contains(&(axis.into(), value)));
        }
    }
    assert!(count > 0);
}

#[test]
fn a_live_contact_cannot_hide_an_unrelated_floating_body_from_newton() {
    with_cx(|cx| {
        for matching in [true, false] {
            let (base, split) = mesh(if matching { [1, 1, 1] } else { [1, 2, 3] }, 0.0);
            let mut positions = base.positions().to_vec();
            let offset = positions.len();
            let mut tets = base.complex().tets.clone();
            let (extra, points) = box_grid([1, 1, 1], [1.0; 3]);
            tets.extend(extra.tets.into_iter().map(|tet| tet.map(|v| v + offset as u32)));
            positions.extend(points.into_iter().map(|p| [p[0] + 4.0, p[1], p[2]]));
            let mesh = ConductionMesh::new(TetComplex::from_tets(positions.len(), tets), positions).unwrap();
            let bc = boundary(&mesh);
            let card = steep_card();
            let set = contacts(cx, &mesh, &bc, &card, split, matching);
            let material = ConductivityModel::isotropic_declared(10.0).unwrap();
            let source = ScalarField::Uniform(0.0);
            let mut initial = field(&mesh, split, 100.0);
            initial[offset..].fill(315.0);
            // The first pair is already at its exact nonlinear-contact
            // equilibrium, and the floating third body's residual is zero.
            let error = set.solve_steady_newton(cx, ConductionProblem {
                mesh: &mesh, boundary: &bc, material: &material, source: &source, element_materials: None,
            }, &initial, config()).unwrap_err();
            assert!(matches!(error, ConductionError::SingularPureNeumann));
        }
    });
}

#[test]
fn g1_coupled_newton_resolves_stiff_contact_feedback_on_both_trace_meshes() {
    with_cx(|cx| {
        for matching in [true, false] {
            let (mesh, split) = mesh(if matching { [1, 1, 1] } else { [1, 2, 3] }, 0.0);
            let bc = boundary(&mesh);
            let card = steep_card();
            let set = contacts(cx, &mesh, &bc, &card, split, matching);
            let material = ConductivityModel::isotropic_declared(10.0).unwrap();
            let source = ScalarField::Uniform(0.0);
            let problem = ConductionProblem {
                mesh: &mesh, boundary: &bc, material: &material,
                element_materials: None, source: &source,
            };
            let initial = field(&mesh, split, 110.0);
            let before = initial.clone();

            // Picard freezes R=.16, giving Q=30/.36 and TA=308+1/3 K,
            // outside this source's support. Its declared budget is ample;
            // silently extrapolating the source cannot rescue that response.
            let frozen = set.evaluate(cx, &bc, &initial).unwrap();
            let response = solve_with_interfaces(
                cx, problem, frozen.frozen_interfaces(), picard_solver(),
            ).unwrap();
            for (i, &value) in response.temperature.iter().enumerate() {
                near(value, field(&mesh, split, 30.0 / 0.36)[i], 2e-7);
            }
            assert!(set.evaluate(cx, &bc, &response.temperature).is_err());
            assert!(set.solve_steady(cx, problem, &initial, picard_solver(), SteadyContactConfig {
                max_iterations: 64, relaxation: 1.0,
                residual_tolerance_w: 1e-8, energy_tolerance_w: 1e-8,
            }).is_err());

            let solved = set.solve_steady_newton(cx, problem, &initial, config()).unwrap();
            check_solution(cx, problem, split, &card, &set, &solved);
            assert_eq!(initial, before);

            // At Q=110, the independent reduced residual is
            // F=Q-(30-.2Q)/(.006Q-.5)=60 and F'=33/8.
            // Hence the FIRST actual full Newton direction has dQ=-160/11.
            // A frozen contact Jacobian cannot pass this assertion even if
            // damping eventually finds the same root.
            let next = field(&mesh, split, 110.0 - 160.0 / 11.0);
            let dofs = DofMap::new(&bc, mesh.vertex_count()).unwrap();
            let expected_norm = dofs.free().iter()
                .map(|&i| (next[i] - initial[i]).powi(2)).sum::<f64>().sqrt();
            near(solved.newton.history[0].newton_step_norm, expected_norm, 2e-7);
            near(solved.newton.history[0].step_length, 1.0, 1e-14);
            if matching {
                assert!(!solved.physical_contacts.states().is_empty());
                assert!(solved.physical_contacts.quadrature_states().is_empty());
            } else {
                assert!(solved.physical_contacts.states().is_empty());
                assert!(!solved.physical_contacts.quadrature_states().is_empty());
            }
        }
    });
}

#[test]
fn g1_newton_stalled_at_field_resolution_is_judged_by_the_physical_gates() {
    // A numerical target far below the residual's rounding floor cannot be
    // met. Newton reaches the floor, and its last line search is rejected with
    // a smallest trial that moves no free temperature. That field passes the
    // unchanged physical gates and is published with its stall diagnosis.
    // With the physical gates at the same unreachable target, the same field
    // is refused by those gates.
    let floor = || {
        let mut result = config();
        result.newton.absolute_tolerance = 1e-30;
        result.max_iterations = 64;
        result.max_residual_evaluations = 4096;
        result
    };
    with_cx(|cx| {
        for matching in [true, false] {
            let (mesh, split) = mesh(if matching { [1, 1, 1] } else { [1, 2, 3] }, 0.0);
            let bc = boundary(&mesh);
            let card = steep_card();
            let set = contacts(cx, &mesh, &bc, &card, split, matching);
            let material = ConductivityModel::isotropic_declared(10.0).unwrap();
            let source = ScalarField::Uniform(0.0);
            let problem = ConductionProblem {
                mesh: &mesh, boundary: &bc, material: &material,
                element_materials: None, source: &source,
            };
            let initial = field(&mesh, split, 110.0);

            let stalled = set.solve_steady_newton(cx, problem, &initial, floor()).unwrap();
            assert!(!stalled.newton.converged);
            assert_eq!(
                stalled.newton.diagnosis,
                Some(fs_solver::NewtonStallDiagnosis::LineSearchRejected),
            );
            assert!(stalled.residual_norm_w <= floor().residual_tolerance_w);
            assert!(stalled.energy.closure_w.abs() <= floor().energy_tolerance_w);
            let converged = set.solve_steady_newton(cx, problem, &initial, config()).unwrap();
            for (&a, &b) in stalled.temperature.iter().zip(&converged.temperature) {
                near(a, b, 1e-9);
            }

            let mut unreachable = floor();
            unreachable.residual_tolerance_w = 1e-30;
            unreachable.energy_tolerance_w = 1e-30;
            let refused = set.solve_steady_newton(cx, problem, &initial, unreachable).unwrap_err();
            assert!(format!("{refused:?}").contains("fails physical gates"), "{refused:?}");
        }
    });
}

#[test]
fn g1_newton_backtracks_source_domain_trials_without_losing_the_pinned_law() {
    with_cx(|cx| {
        for matching in [true, false] {
            let (mesh, split) = mesh(if matching { [1, 1, 1] } else { [1, 2, 3] }, 0.0);
            let bc = boundary(&mesh);
            let card = steep_card();
            let set = contacts(cx, &mesh, &bc, &card, split, matching);
            let material = ConductivityModel::isotropic_declared(10.0).unwrap();
            let source = ScalarField::Uniform(0.0);
            let problem = ConductionProblem {
                mesh: &mesh, boundary: &bc, material: &material,
                element_materials: None, source: &source,
            };
            let initial = field(&mesh, split, 119.0);
            let before = initial.clone();
            // Exact scalar full Newton trial: Q=2711700/31449, whose
            // TA=308.6225317... K is unsupported. Half the step is supported.
            let full_trial = field(&mesh, split, 2_711_700.0 / 31_449.0);
            assert!(set.evaluate(cx, &bc, &full_trial).is_err());
            let half_trial = field(&mesh, split, (119.0 + 2_711_700.0 / 31_449.0) / 2.0);
            assert!(set.evaluate(cx, &bc, &half_trial).is_ok());

            let solved = set.solve_steady_newton(cx, problem, &initial, config()).unwrap();
            check_solution(cx, problem, split, &card, &set, &solved);
            near(solved.newton.history[0].step_length, 0.5, 1e-14);
            assert!(solved.residual_evaluations >= solved.newton.iterations + 2);
            assert_eq!(initial, before);
        }
    });
}

#[test]
fn g1_steady_newton_reuses_primal_tangents_at_source_endpoints_and_corners() {
    with_cx(|cx| {
        let (mesh, split) = mesh([1, 1, 1], 0.0);
        let bc = boundary(&mesh);
        let material = ConductivityModel::isotropic_declared(10.0).unwrap();
        let source = ScalarField::Uniform(0.0);
        let problem = ConductionProblem {
            mesh: &mesh, boundary: &bc, material: &material,
            element_materials: None, source: &source,
        };
        for (card, initial_heat) in [
            (steep_card(), 90.0),
            (steep_card(), 120.0),
            // The root remains Q=100 W on the original lower segment. Only
            // the initial TA=311 K is a genuine source-slope corner.
            (card_with_knots(vec![(309.0, 0.04), (311.0, 0.16), (312.0, 0.26)]), 110.0),
        ] {
            let set = contacts(cx, &mesh, &bc, &card, split, true);
            let initial = field(&mesh, split, initial_heat);
            let before = initial.clone();
            assert!(set.evaluate(cx, &bc, &initial).is_ok());
            // A forward tangent must not silently relax the independent
            // two-sided sensitivity contract at these very same states.
            assert!(set.linearize_temperature(cx, &bc, &initial, 10_000).is_err());
            let solved = set.solve_steady_newton(cx, problem, &initial, config()).unwrap();
            check_solution(cx, problem, split, &card, &set, &solved);
            assert_eq!(initial, before);
        }
    });
}

#[test]
fn g0_newton_work_source_and_cancellation_refusals_are_transactional() {
    with_cx(|cx| {
        let (mesh, split) = mesh([1, 2, 3], 0.0);
        let bc = boundary(&mesh);
        let card = steep_card();
        let set = contacts(cx, &mesh, &bc, &card, split, false);
        let material = ConductivityModel::isotropic_declared(10.0).unwrap();
        let source = ScalarField::Uniform(0.0);
        let problem = ConductionProblem {
            mesh: &mesh, boundary: &bc, material: &material,
            element_materials: None, source: &source,
        };
        let initial = field(&mesh, split, 110.0);
        let before = initial.clone();
        assert!(set.solve_steady_newton(cx, problem, &initial, SteadyContactNewtonConfig {
            max_residual_evaluations: 1, ..config()
        }).is_err());
        assert!(set.solve_steady_newton(cx, problem, &initial, SteadyContactNewtonConfig {
            max_iterations: 1, ..config()
        }).is_err());
        assert!(set.solve_steady_newton(cx, problem, &initial, SteadyContactNewtonConfig {
            max_face_pairs: 1, ..config()
        }).is_err());
        let unsupported = field(&mesh, split, 80.0);
        let unsupported_before = unsupported.clone();
        assert!(set.solve_steady_newton(cx, problem, &unsupported, config()).is_err());
        assert_eq!(unsupported, unsupported_before);
        support::with_cancelled_cx(|cancelled| {
            assert!(matches!(
                set.solve_steady_newton(cancelled, problem, &initial, config()),
                Err(ConductionError::Cancelled { .. })
            ));
        });
        assert_eq!(initial, before);
        let solved = set.solve_steady_newton(cx, problem, &initial, config()).unwrap();
        check_solution(cx, problem, split, &card, &set, &solved);
    });
}

fn mean_boundary(mesh: &ConductionMesh) -> ThermalBoundary {
    ThermalBoundaryBuilder::new(mesh)
        .region(
            "out",
            |f| on_box_face(f.centroid[0], 0.0),
            ThermalBc::neumann(100.0).unwrap(),
        )
        .unwrap()
        .region(
            "in",
            |f| on_box_face(f.centroid[0], 2.0),
            ThermalBc::neumann(-100.0).unwrap(),
        )
        .unwrap()
        .adiabatic_remainder()
        .finish()
        .unwrap()
}

fn mean_config(initial: &[f64]) -> MeanTemperatureNewtonConfig {
    let mut result = MeanTemperatureNewtonConfig::default();
    result.initial_temperature = Some(initial.to_vec());
    result.linear.tolerance = 1e-12;
    result.linear.max_iterations = 2048;
    result.linear.restart = 64;
    result.stop.residual_rtol = 1e-12;
    result.stop.residual_atol = 1e-9;
    result.stop.step_atol = 0.0;
    result.stop.max_iterations = 12;
    result
}

// Independent affine states with mean 315 K. The Neumann boundaries do not
// fix the far-end temperatures: trace temperature and bulk heat are separate.
fn mean_field(mesh: &ConductionMesh, split: usize, ta: f64, heat: f64) -> Vec<f64> {
    mesh.positions()
        .iter()
        .enumerate()
        .map(|(i, p)| {
            let trace = if i < split { ta } else { 630.0 - ta };
            trace + 0.1 * heat * (p[0] - 1.0)
        })
        .collect()
}

fn mean_residual(
    cx: &fs_exec::Cx<'_>,
    problem: ConductionProblem<'_>,
    set: &TemperatureContacts<'_, '_>,
    temperature: &[f64],
) -> Vec<f64> {
    let physical = set.evaluate(cx, problem.boundary, temperature).unwrap();
    let system = assemble_operator_with_interfaces(
        cx,
        problem.mesh,
        problem.boundary,
        problem.material,
        problem.source,
        temperature,
        physical.frozen_interfaces(),
    )
    .unwrap();
    let mut residual = vec![0.0; temperature.len()];
    system.operator.spmv(temperature, &mut residual);
    for (value, load) in residual.iter_mut().zip(&system.load) {
        *value -= load;
    }
    residual
}

fn check_mean_solution(
    cx: &fs_exec::Cx<'_>,
    problem: ConductionProblem<'_>,
    split: usize,
    card: &InterfaceSystemCard,
    set: &TemperatureContacts<'_, '_>,
    solved: &MeanTemperatureContactSolution,
) {
    let mesh = problem.mesh;
    let temperature = &solved.conduction.temperature;
    let report = &solved.conduction.report;
    // Q=100 is imposed by the flux boundaries. Equal slab volumes give
    // TA+TB=630; TB-TA=100*(.04+.06*(TA-309)) then fixes TA=310, TB=320.
    assert_eq!(temperature.len(), mesh.vertex_count());
    for (&actual, expected) in temperature.iter().zip(field(mesh, split, 100.0)) {
        near(actual, expected, 2e-7);
    }
    let volume_mean = mesh
        .complex()
        .tets
        .iter()
        .enumerate()
        .map(|(e, tet)| {
            mesh.element_volume(e) * tet.iter().map(|&v| temperature[v as usize]).sum::<f64>() / 4.0
        })
        .sum::<f64>()
        / mesh.total_volume();
    near(volume_mean, 315.0, 2e-10);

    // Check EVERY original row, including the algebraically eliminated row 0,
    // against a fresh source query at the unchanged returned field.
    let residual = mean_residual(cx, problem, set, temperature);
    for value in &residual {
        assert!(value.abs() <= 1e-8, "unbalanced original row: {value} W");
    }
    let norm = residual
        .iter()
        .map(|value| value * value)
        .sum::<f64>()
        .sqrt();
    near(report.final_residual, norm, 2e-10);
    assert!(norm <= 1e-8);
    assert!(report.final_residual <= report.residual_threshold);
    assert_eq!(report.stop_reason, StopReason::ResidualTolerance);
    assert!(report.iterations > 0 && report.iterations <= 12);
    assert_eq!(report.free_dofs, mesh.vertex_count());
    assert!(
        DofMap::new(problem.boundary, mesh.vertex_count())
            .unwrap()
            .fixed()
            .is_empty()
    );
    assert_eq!(report.energy.dirichlet_in_w, 0.0);
    near(report.energy.source_w, 0.0, 1e-12);
    near(report.energy.neumann_out_w, 0.0, 1e-10);
    near(report.energy.closure_w, 0.0, 1e-10);
    assert!(report.robin_fluxes.is_empty());

    let flux = solved
        .physical_contacts
        .frozen_interfaces()
        .fluxes(temperature)
        .unwrap();
    assert_eq!(flux.len(), 1);
    assert_eq!(report.interface_fluxes.len(), 1);
    near(flux[0].heat_rate_a_to_b_w, -100.0, 2e-6);
    near(report.interface_fluxes[0].heat_rate_a_to_b_w, -100.0, 2e-6);
    near(flux[0].area_m2, 1.0, 2e-12);
    let contact = set
        .linearize_temperature(cx, problem.boundary, temperature, 10_000)
        .unwrap();
    near(
        contact.nodal_outflow_w()[..split].iter().sum(),
        -100.0,
        2e-6,
    );
    near(contact.nodal_outflow_w()[split..].iter().sum(), 100.0, 2e-6);
    near(contact.nodal_outflow_w().iter().sum(), 0.0, 2e-10);

    let pinned = InterfaceResistance::from_card(
        "bond",
        card,
        &reference_point(),
        SelectionPolicy::SingleClaimOnly,
    )
    .unwrap();
    let states = solved
        .physical_contacts
        .states()
        .iter()
        .map(|state| {
            (
                state.side_a_temperature_k,
                state.side_b_temperature_k,
                &state.resistance,
            )
        })
        .chain(
            solved
                .physical_contacts
                .quadrature_states()
                .iter()
                .map(|state| {
                    (
                        state.side_a_temperature_k,
                        state.side_b_temperature_k,
                        &state.resistance,
                    )
                }),
        );
    let mut count = 0;
    for (ta, tb, resistance) in states {
        count += 1;
        near(ta, 310.0, 2e-7);
        near(tb, 320.0, 2e-7);
        near(resistance.value_m2_k_per_w(), 0.1, 2e-8);
        assert_eq!(resistance.card_identity(), card.content_hash());
        assert_eq!(resistance.receipt().selected, pinned.receipt().selected);
        for (axis, value) in [("TA", ta), ("TB", tb), ("normal_pressure", 1e6)] {
            assert!(
                resistance
                    .receipt()
                    .query_point
                    .contains(&(axis.into(), value))
            );
        }
    }
    assert!(count > 0);
}

#[test]
fn g1_mean_temperature_newton_resolves_live_contact_without_a_thermostat() {
    with_cx(|cx| {
        for matching in [true, false] {
            let (mesh, split) = mesh(if matching { [1, 1, 1] } else { [1, 2, 3] }, 0.0);
            let bc = mean_boundary(&mesh);
            let card = steep_card();
            let set = contacts(cx, &mesh, &bc, &card, split, matching);
            let material = ConductivityModel::isotropic_declared(10.0).unwrap();
            let source = ScalarField::Uniform(0.0);
            let problem = ConductionProblem {
                mesh: &mesh,
                boundary: &bc,
                material: &material,
                element_materials: None,
                source: &source,
            };
            let initial = field(&mesh, split, 110.0);
            let before = initial.clone();

            // The complete Newton full step has TA=309.4 and Q=100. Its
            // contact heat is 175 W, so Armijo accepts the half step instead.
            // Its independent residual pins dR/dT, not merely the final root.
            let half_trial = mean_field(&mesh, split, 310.2, 105.0);
            let expected = mean_residual(cx, problem, &set, &half_trial);
            let expected_norm = expected.iter().map(|r| r * r).sum::<f64>().sqrt();
            let mut one_step = mean_config(&initial);
            one_step.stop.max_iterations = 1;
            match solve_mean_temperature_newton_with_contacts(
                cx, problem, &set, 315.0, 10_000, one_step,
            ) {
                Err(ConductionError::NotConverged {
                    iterations: 1,
                    residual,
                    ..
                }) => {
                    near(residual, expected_norm, 2e-7);
                }
                other => panic!("expected the first complete contact Newton step, got {other:?}"),
            }
            let solved = solve_mean_temperature_newton_with_contacts(
                cx,
                problem,
                &set,
                315.0,
                10_000,
                mean_config(&initial),
            )
            .unwrap();
            check_mean_solution(cx, problem, split, &card, &set, &solved);
            if !matching {
                // The finer B trace has more nodes, but does not own more volume.
                let nodal_mean =
                    solved.conduction.temperature.iter().sum::<f64>() / mesh.vertex_count() as f64;
                assert!((nodal_mean - 315.0).abs() > 1.0);
            }
            assert_eq!(initial, before);
        }
    });
}

#[test]
fn g1_mean_temperature_contact_backtracks_original_source_domain_refusals() {
    with_cx(|cx| {
        for matching in [true, false] {
            let (mesh, split) = mesh(if matching { [1, 1, 1] } else { [1, 2, 3] }, 0.0);
            let bc = mean_boundary(&mesh);
            let card = steep_card();
            let set = contacts(cx, &mesh, &bc, &card, split, matching);
            let material = ConductivityModel::isotropic_declared(10.0).unwrap();
            let source = ScalarField::Uniform(0.0);
            let problem = ConductionProblem {
                mesh: &mesh,
                boundary: &bc,
                material: &material,
                element_materials: None,
                source: &source,
            };
            let initial = field(&mesh, split, 119.0);
            let before = initial.clone();
            // At TA=311.9, TB=318.1, R=.214: the complete constrained
            // Newton correction is dTA=(6.2-100*.214)*.214/(2*.214+.06*6.2)
            // = -4.066 K. The bulk flux reaches 100 W on the full step.
            let full_trial = mean_field(&mesh, split, 307.834, 100.0);
            let half_trial = mean_field(&mesh, split, 309.867, 109.5);
            assert!(set.evaluate(cx, &bc, &full_trial).is_err());
            assert!(set.evaluate(cx, &bc, &half_trial).is_ok());
            let mut no_backtracking = mean_config(&initial);
            no_backtracking.line_search.max_backtracks = 0;
            assert!(matches!(
                solve_mean_temperature_newton_with_contacts(
                    cx, problem, &set, 315.0, 10_000, no_backtracking,
                ),
                Err(ConductionError::Interface { interface, what, .. })
                    if interface == "bond" && what.contains("fs-matdb refused")
            ));
            assert_eq!(initial, before);

            let expected = mean_residual(cx, problem, &set, &half_trial);
            let expected_norm = expected.iter().map(|r| r * r).sum::<f64>().sqrt();
            let mut one_step = mean_config(&initial);
            one_step.stop.max_iterations = 1;
            match solve_mean_temperature_newton_with_contacts(
                cx, problem, &set, 315.0, 10_000, one_step,
            ) {
                Err(ConductionError::NotConverged {
                    iterations: 1,
                    residual,
                    ..
                }) => {
                    near(residual, expected_norm, 2e-7);
                }
                other => panic!("expected the source-admitted half step, got {other:?}"),
            }
            let solved = solve_mean_temperature_newton_with_contacts(
                cx,
                problem,
                &set,
                315.0,
                10_000,
                mean_config(&initial),
            )
            .unwrap();
            check_mean_solution(cx, problem, split, &card, &set, &solved);
            assert_eq!(initial, before);
        }
    });
}

#[test]
fn g0_mean_temperature_contact_work_balance_and_cancellation_limits_hold() {
    with_cx(|cx| {
        let (mesh, split) = mesh([1, 2, 3], 0.0);
        let bc = mean_boundary(&mesh);
        let card = steep_card();
        let set = contacts(cx, &mesh, &bc, &card, split, false);
        let material = ConductivityModel::isotropic_declared(10.0).unwrap();
        let source = ScalarField::Uniform(0.0);
        let problem = ConductionProblem {
            mesh: &mesh,
            boundary: &bc,
            material: &material,
            element_materials: None,
            source: &source,
        };
        let initial = field(&mesh, split, 110.0);
        let before = initial.clone();
        let error = solve_mean_temperature_newton_with_contacts(
            cx,
            problem,
            &set,
            315.0,
            1,
            mean_config(&initial),
        )
        .unwrap_err();
        assert!(error.to_string().contains("contact tangent exceeds"));
        let unbalanced = ScalarField::Uniform(1.0);
        let error = solve_mean_temperature_newton_with_contacts(
            cx,
            ConductionProblem {
                source: &unbalanced,
                ..problem
            },
            &set,
            315.0,
            10_000,
            mean_config(&initial),
        )
        .unwrap_err();
        assert!(error.to_string().contains("no steady heat balance"));
        support::with_cancelled_cx(|cancelled| {
            assert!(matches!(
                solve_mean_temperature_newton_with_contacts(
                    cancelled,
                    problem,
                    &set,
                    315.0,
                    10_000,
                    mean_config(&initial),
                ),
                Err(ConductionError::Cancelled { .. })
            ));
        });
        assert_eq!(initial, before);
    });
}
