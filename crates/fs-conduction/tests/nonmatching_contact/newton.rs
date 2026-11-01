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

// Two unit slabs, k=10 W/(m K), with 300/330 K reservoirs. Each bulk
// resistance is 0.1 K/W. This positive source law has deliberately strong
// feedback, but is smooth throughout its admitted thermal neighbourhood.
fn steep_card() -> InterfaceSystemCard {
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
            knots: vec![(309.0, 0.04), (312.0, 0.22)],
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
