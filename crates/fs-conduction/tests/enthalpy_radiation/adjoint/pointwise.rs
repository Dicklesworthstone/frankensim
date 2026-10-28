//! Exact local-emission reverse against independently perturbed physical solves.
use super::*;
use fs_conduction::PointwiseAmbientRadiation;
use fs_conduction::transient::enthalpy::EnthalpyStepSolution;
use fs_conduction::transient::enthalpy::adjoint::pointwise::workspace_entries;

fn controls() -> EnthalpyStepConfig {
    let mut config = config();
    config.newton.absolute_tolerance = 1e-10;
    config.newton.relative_tolerance = 0.0;
    config
}

fn assigned<'a>(
    fixture: &'a Fixture,
    cx: &fs_exec::Cx<'_>,
    density: f64,
) -> HeterogeneousEnthalpyBackwardEuler<'a, 'a> {
    HeterogeneousEnthalpyBackwardEuler::new(
        cx,
        &fixture.mesh,
        &[
            ReferenceEnthalpyMaterial {
                curve: &fixture.curve,
                reference_density_kg_m3: 3.0 * RHO,
            },
            ReferenceEnthalpyMaterial {
                curve: &fixture.curve,
                reference_density_kg_m3: density,
            },
        ],
        &[1],
        EnthalpyBudget { max_elements: 2, ..budget() },
    )
    .unwrap()
}

#[allow(clippy::too_many_arguments)]
fn endpoint(
    fixture: &Fixture,
    cx: &fs_exec::Cx<'_>,
    boundary: &ThermalBoundary,
    old: &[f64],
    source: &[f64],
    ambient: f64,
    epsilon: f64,
    density: f64,
) -> EnthalpyStepSolution {
    let source = ScalarField::nodal("heat density", 4, source.to_vec()).unwrap();
    let patches = [patch_with_emissivity(ambient, 520.0, epsilon)];
    let radiation = PointwiseAmbientRadiation::new(cx, &fixture.mesh, boundary, &patches).unwrap();
    assigned(fixture, cx, density)
        .advance_with_pointwise_ambient_radiation(
            cx,
            ConductionProblem { boundary, ..fixture.problem(&source) },
            None,
            old,
            DT,
            controls(),
            &radiation,
        )
        .unwrap()
}

fn objective(solution: &EnthalpyStepSolution) -> f64 {
    solution.temperature.iter().zip(GOAL).map(|(t, w)| w * (t - 350.0)).sum()
}

#[test]
fn local_emission_reverse_matches_nonuniform_physics_and_storage_controls() {
    let fixture = Fixture::new();
    with_cx(|cx| {
        for ambient in [300.0, 500.0] {
            let source = ScalarField::nodal("heat density", 4, SOURCE.to_vec()).unwrap();
            let patches = [patch_with_emissivity(ambient, 520.0, EPSILON)];
            let radiation = PointwiseAmbientRadiation::new(
                cx, &fixture.mesh, &fixture.boundary, &patches,
            ).unwrap();
            let engine = assigned(&fixture, cx, RHO);
            let accepted = endpoint(
                &fixture, cx, &fixture.boundary, &OLD, &SOURCE, ambient, EPSILON, RHO,
            );
            let entries = workspace_entries(4, fixture.mesh.boundary().len(), 1).unwrap();
            let lin = engine.linearize_accepted_with_pointwise_ambient_radiation(
                cx, fixture.problem(&source), None, &OLD, DT, controls(), &radiation,
                accepted.clone(), controls().newton.absolute_tolerance, entries,
            ).unwrap();
            assert_eq!(lin.primal().specific_enthalpy_j_kg, accepted.specific_enthalpy_j_kg);
            assert_eq!(lin.primal().temperature, accepted.temperature);
            assert_eq!(lin.primal().liquid_mass_fraction, accepted.liquid_mass_fraction);
            assert_eq!(lin.primal().pointwise_radiation, accepted.pointwise_radiation);
            let range = accepted.temperature.iter().copied().fold(f64::NEG_INFINITY, f64::max)
                - accepted.temperature.iter().copied().fold(f64::INFINITY, f64::min);
            assert!(range > 100.0);
            let transport = lin.transport();
            let u = [0.2, -0.5, 0.7, 0.3];
            let v = [-1.0, 0.1, 0.6, 0.8];
            let jv = transport.apply_jacobian(cx, &v).unwrap();
            let jtu = transport.apply_jacobian_transpose(cx, &u).unwrap();
            assert!((dot(&u, &jv) - dot(&jtu, &v)).abs() < 2e-12);
            let assembly = radiation.assemble(cx, &accepted.temperature).unwrap();
            let mut rv = [0.0; 4];
            let mut rtu = [0.0; 4];
            assembly.tangent_apply(cx, &v, &mut rv).unwrap();
            assembly.tangent_apply_transpose(cx, &u, &mut rtu).unwrap();
            assert!((dot(&u, &rv) - dot(&rtu, &v)).abs() < 2e-11);
            let seed = transport.temperature_pullback(cx, &GOAL).unwrap();
            let gradient = lin.pullback(cx, &seed, linear()).unwrap();
            let outward_weights: Vec<_> = transport
                .nodal_load_pullback(cx, &gradient.transport.adjoint).unwrap()
                .into_iter().map(|x| -x).collect();
            let parameters = transport.pointwise_parameter_pullback(cx, &outward_weights).unwrap();
            assert_eq!(parameters.radiation_regions, ["exterior"]);
            assert!(gradient.transport.relative_residual < linear().tolerance);
            assert_eq!(gradient.material_reference_density[0], 0.0);
            assert_eq!(gradient.material_reference_density[1], gradient.element_reference_density[0]);
            let history_direction = [0.3, -0.7, 0.2, 0.5];
            agree(
                dot(&gradient.transport.previous_specific_enthalpy, &history_direction),
                five_point(0.03, |delta| {
                    let old = std::array::from_fn::<_, 4, _>(|i| OLD[i] + delta * history_direction[i]);
                    objective(&endpoint(&fixture, cx, &fixture.boundary, &old, &SOURCE, ambient, EPSILON, RHO))
                }),
                "pointwise history",
            );
            let source_direction = [0.5, 0.2, -0.4, 0.8];
            agree(
                dot(&gradient.transport.source_density, &source_direction),
                five_point(0.1, |delta| {
                    let q = std::array::from_fn::<_, 4, _>(|i| SOURCE[i] + delta * source_direction[i]);
                    objective(&endpoint(&fixture, cx, &fixture.boundary, &OLD, &q, ambient, EPSILON, RHO))
                }),
                "pointwise source",
            );
            agree(
                gradient.material_reference_density[1],
                five_point(0.01, |delta| objective(&endpoint(
                    &fixture, cx, &fixture.boundary, &OLD, &SOURCE, ambient, EPSILON, RHO + delta,
                ))),
                "pointwise reference density",
            );
            agree(
                parameters.emissivities[0],
                five_point(1e-4, |delta| objective(&endpoint(
                    &fixture, cx, &fixture.boundary, &OLD, &SOURCE, ambient, EPSILON + delta, RHO,
                ))),
                "pointwise emissivity",
            );
            agree(
                parameters.reservoir_temperatures[0],
                five_point(0.05, |delta| objective(&endpoint(
                    &fixture, cx, &fixture.boundary, &OLD, &SOURCE, ambient + delta, EPSILON, RHO,
                ))),
                "pointwise reservoir",
            );
            let density_action = lin.apply_reference_density_jacobian(cx, &[0.3]).unwrap();
            agree(
                -dot(&gradient.transport.adjoint, &density_action),
                0.3 * gradient.element_reference_density[0],
                "pointwise density residual transpose",
            );
            assert_eq!(lin.pullback(cx, &seed, linear()).unwrap(), gradient);
        }
    });
}

fn bath_boundary(fixture: &Fixture, temperature: f64) -> ThermalBoundary {
    ThermalBoundaryBuilder::new(&fixture.mesh)
        .region("bath", |face| face.centroid[2] == 0.0, ThermalBc::dirichlet(temperature).unwrap())
        .unwrap()
        .remainder("exterior", ThermalBc::robin(HTC, 300.0).unwrap())
        .unwrap().finish().unwrap()
}

#[test]
fn pointwise_bath_lift_keeps_fixed_columns_and_masks_parameter_load_rows() {
    let fixture = Fixture::new();
    let boundary = bath_boundary(&fixture, 400.0);
    with_cx(|cx| {
        let source = ScalarField::nodal("heat density", 4, SOURCE.to_vec()).unwrap();
        let problem = ConductionProblem { boundary: &boundary, ..fixture.problem(&source) };
        let patches = [patch_with_emissivity(500.0, 520.0, EPSILON)];
        let radiation = PointwiseAmbientRadiation::new(cx, &fixture.mesh, &boundary, &patches).unwrap();
        let accepted = endpoint(&fixture, cx, &boundary, &OLD, &SOURCE, 500.0, EPSILON, RHO);
        assert!(accepted.dirichlet_in_w.abs() > 1.0);
        let engine = assigned(&fixture, cx, RHO);
        let lin = engine.linearize_accepted_with_pointwise_ambient_radiation(
            cx, problem, None, &OLD, DT, controls(), &radiation, accepted,
            controls().newton.absolute_tolerance,
            workspace_entries(4, fixture.mesh.boundary().len(), 1).unwrap(),
        ).unwrap();
        let transport = lin.transport();
        let seed = transport.temperature_pullback(cx, &GOAL).unwrap();
        let gradient = lin.pullback(cx, &seed, linear()).unwrap();
        let nodal_load = transport.nodal_load_pullback(cx, &gradient.transport.adjoint).unwrap();
        for &(vertex, _) in boundary.dirichlet() {
            let vertex = vertex as usize;
            assert_eq!(nodal_load[vertex], 0.0);
            assert_eq!(gradient.transport.previous_specific_enthalpy[vertex], 0.0);
        }
        agree(
            gradient.transport.prescribed_temperature.iter().sum(),
            five_point(0.02, |delta| {
                let shifted = bath_boundary(&fixture, 400.0 + delta);
                objective(&endpoint(&fixture, cx, &shifted, &OLD, &SOURCE, 500.0, EPSILON, RHO))
            }),
            "pointwise prescribed bath lift",
        );
        let outward_weights: Vec<_> = nodal_load.into_iter().map(|x| -x).collect();
        let parameters = transport.pointwise_parameter_pullback(cx, &outward_weights).unwrap();
        agree(
            parameters.emissivities[0],
            five_point(1e-4, |delta| objective(&endpoint(
                &fixture, cx, &boundary, &OLD, &SOURCE, 500.0, EPSILON + delta, RHO,
            ))),
            "pointwise emissivity with prescribed rows",
        );
    });
}

#[test]
fn latent_pointwise_controls_change_stored_enthalpy_without_changing_temperature() {
    let fixture = Fixture::new();
    with_cx(|cx| {
        let old = [2000.0; 4];
        let source = ScalarField::Uniform(0.0);
        let patches = [patch_with_emissivity(500.0, 520.0, EPSILON)];
        let radiation = PointwiseAmbientRadiation::new(
            cx, &fixture.mesh, &fixture.boundary, &patches,
        ).unwrap();
        let accepted = endpoint(&fixture, cx, &fixture.boundary, &old, &[0.0; 4], 500.0, EPSILON, RHO);
        assert_eq!(accepted.temperature, [350.0; 4]);
        assert!(accepted.specific_enthalpy_j_kg.iter().all(|h| *h > 2000.0 && *h < 3000.0));
        let engine = assigned(&fixture, cx, RHO);
        let lin = engine.linearize_accepted_with_pointwise_ambient_radiation(
            cx, fixture.problem(&source), None, &old, DT, controls(), &radiation, accepted,
            controls().newton.absolute_tolerance,
            workspace_entries(4, fixture.mesh.boundary().len(), 1).unwrap(),
        ).unwrap();
        let transport = lin.transport();
        let temperature_seed = transport.temperature_pullback(cx, &[1.0; 4]).unwrap();
        assert_eq!(temperature_seed, [0.0; 4]);
        assert_eq!(lin.pullback(cx, &temperature_seed, linear()).unwrap().transport.adjoint, [0.0; 4]);
        let carry = lin.pullback(cx, &[1.0; 4], linear()).unwrap();
        let outward_weights: Vec<_> = transport.nodal_load_pullback(cx, &carry.transport.adjoint)
            .unwrap().into_iter().map(|x| -x).collect();
        let parameters = transport.pointwise_parameter_pullback(cx, &outward_weights).unwrap();
        assert!(parameters.emissivities[0] > 1.0);
        agree(
            parameters.emissivities[0],
            five_point(1e-4, |delta| endpoint(
                &fixture, cx, &fixture.boundary, &old, &[0.0; 4], 500.0, EPSILON + delta, RHO,
            ).specific_enthalpy_j_kg.iter().sum()),
            "pointwise latent enthalpy carry",
        );
    });
}

#[test]
fn pointwise_reverse_rechecks_the_endpoint_and_preserves_it_after_refusal() {
    let fixture = Fixture::new();
    with_cx(|cx| {
        let source = ScalarField::nodal("heat density", 4, SOURCE.to_vec()).unwrap();
        let patches = [patch_with_emissivity(300.0, 520.0, EPSILON)];
        let radiation = PointwiseAmbientRadiation::new(
            cx, &fixture.mesh, &fixture.boundary, &patches,
        ).unwrap();
        let accepted = endpoint(&fixture, cx, &fixture.boundary, &OLD, &SOURCE, 300.0, EPSILON, RHO);
        let engine = assigned(&fixture, cx, RHO);
        let entries = workspace_entries(4, fixture.mesh.boundary().len(), 1).unwrap();
        let target = controls().newton.absolute_tolerance;
        let bind = |cx: &fs_exec::Cx<'_>, state: EnthalpyStepSolution, target, entries| {
            engine.linearize_accepted_with_pointwise_ambient_radiation(
                cx, fixture.problem(&source), None, &OLD, DT, controls(), &radiation,
                state, target, entries,
            )
        };
        assert!(bind(cx, accepted.clone(), target, entries - 1).is_err());
        assert!(bind(cx, accepted.clone(), f64::NAN, entries).is_err());
        let mut invalid = controls();
        invalid.newton.relative_tolerance = f64::NAN;
        assert!(engine.linearize_accepted_with_pointwise_ambient_radiation(
            cx, fixture.problem(&source), None, &OLD, DT, invalid, &radiation,
            accepted.clone(), target, entries,
        ).is_err(), "an explicit physical target must not bypass Newton configuration admission");
        let mut stale = accepted.clone();
        stale.specific_enthalpy_j_kg[0] += 10.0;
        assert!(bind(cx, stale, target, entries).is_err());
        assert!(engine.linearize_accepted(
            cx, fixture.problem(&source), None, &OLD, DT, controls(), accepted.clone(),
        ).is_err());
        let lin = bind(cx, accepted.clone(), target, entries).unwrap();
        let seed = lin.transport().temperature_pullback(cx, &GOAL).unwrap();
        let gradient = lin.pullback(cx, &seed, linear()).unwrap();
        assert!(lin.pullback(cx, &seed, LinearConfig { max_iterations: 0, ..linear() }).is_err());
        with_cancelled_cx(|cancelled| {
            assert!(bind(cancelled, accepted.clone(), target, entries).is_err());
            assert!(lin.pullback(cancelled, &seed, linear()).is_err());
            assert!(lin.transport().pointwise_parameter_pullback(cancelled, &[1.0; 4]).is_err());
        });
        assert_eq!(lin.primal().specific_enthalpy_j_kg, accepted.specific_enthalpy_j_kg);
        assert_eq!(lin.primal().temperature, accepted.temperature);
        assert_eq!(lin.pullback(cx, &seed, linear()).unwrap(), gradient);
        assert!(workspace_entries(usize::MAX, 4, 1).is_err());
        assert!(workspace_entries(4, 4, 0).is_err());
        assert!(workspace_entries(4, 4, 65).is_err());
    });
}
