//! Reflected local-emission transpose and physical enthalpy-history controls.
use super::*;
use fs_conduction::transient::enthalpy::EnthalpyStepSolution;
use fs_conduction::transient::enthalpy::adjoint::pointwise::workspace_entries;
use fs_conduction::transient::enthalpy::heterogeneous::{
    HeterogeneousEnthalpyBackwardEuler, ReferenceEnthalpyMaterial,
};

const OLD: [f64; 8] = [600.0, 1800.0, 3400.0, 4000.0, 900.0, 1700.0, 3100.0, 3900.0];
const SOURCE: [f64; 8] = [12.0, -3.0, 8.0, 5.0, -4.0, 7.0, 2.0, 9.0];
const GOAL: [f64; 8] = [0.7, -0.4, 0.2, 1.1, -0.3, 0.6, 0.9, -0.2];
const DENSITIES: [f64; 2] = [RHO, 1.7 * RHO];

fn controls() -> EnthalpyStepConfig {
    let mut config = step_config();
    config.newton.absolute_tolerance = 1e-10;
    config.newton.relative_tolerance = 0.0;
    config
}

fn linear() -> LinearConfig {
    LinearConfig { tolerance: 1e-11, max_iterations: 128, restart: 8 }
}

fn assigned<'a>(
    fixture: &'a Fixture,
    cx: &fs_exec::Cx<'_>,
    densities: [f64; 2],
    phase_k: bool,
) -> HeterogeneousEnthalpyBackwardEuler<'a, 'a> {
    let engine = HeterogeneousEnthalpyBackwardEuler::new(
        cx,
        &fixture.mesh,
        &[
            ReferenceEnthalpyMaterial {
                curve: &fixture.curve,
                reference_density_kg_m3: densities[0],
            },
            ReferenceEnthalpyMaterial {
                curve: &fixture.curve,
                reference_density_kg_m3: densities[1],
            },
            ReferenceEnthalpyMaterial {
                curve: &fixture.curve,
                reference_density_kg_m3: 3.0 * RHO,
            },
        ],
        &[0, 1],
        EnthalpyBudget { max_vertices: 8, max_elements: 3 },
    ).unwrap();
    if phase_k {
        engine.with_phase_conductivity(cx, vec![
            LiquidMassFractionConductivity::declared(0.7, 2.2).unwrap(),
            LiquidMassFractionConductivity::declared(1.1, 0.6).unwrap(),
        ]).unwrap()
    } else {
        engine
    }
}

fn model(fixture: &Fixture, epsilon: [f64; 2]) -> GrayDiffuseEnclosure {
    let left = RadiationSurface::new(
        &fixture.mesh, "left", |face| face.vertices.iter().all(|v| *v < 4),
        emissivity("left", epsilon[0]),
    ).unwrap();
    let right = RadiationSurface::new(
        &fixture.mesh, "right", |face| face.vertices.iter().all(|v| *v >= 4),
        emissivity("right", epsilon[1]),
    ).unwrap();
    let factors = ViewFactorMatrix::infinite_parallel_plates(left.area_m2()).unwrap();
    GrayDiffuseEnclosure::new(vec![left, right], factors).unwrap()
}

#[allow(clippy::too_many_arguments)]
fn endpoint(
    fixture: &Fixture,
    cx: &fs_exec::Cx<'_>,
    boundary: &ThermalBoundary,
    old: &[f64],
    source: &[f64],
    epsilon: [f64; 2],
    densities: [f64; 2],
    phase_k: bool,
) -> EnthalpyStepSolution {
    let source = ScalarField::nodal("heat density", 8, source.to_vec()).unwrap();
    let enclosure = model(fixture, epsilon);
    let radiation = PointwiseEnclosureRadiation::new(
        cx, &fixture.mesh, boundary, &enclosure, PointwiseEnclosureTolerance::default(),
    ).unwrap();
    assigned(fixture, cx, densities, phase_k)
        .advance_with_pointwise_enclosure_radiation(
            cx, problem(fixture, boundary, &source), None, old, DT, controls(), &radiation,
        ).unwrap()
}

fn objective(solution: &EnthalpyStepSolution) -> f64 {
    solution.temperature.iter().zip(GOAL).map(|(t, w)| w * (t - 350.0)).sum()
}

fn five_point(step: f64, f: impl Fn(f64) -> f64) -> f64 {
    (-f(2.0 * step) + 8.0 * f(step) - 8.0 * f(-step) + f(-2.0 * step)) / (12.0 * step)
}

fn agree(actual: f64, expected: f64, label: &str) {
    assert!(
        (actual - expected).abs() <= 3e-7 * (1.0 + expected.abs()),
        "{label}: implicit {actual:.12e}, independent perturbation {expected:.12e}"
    );
}

#[test]
fn reflected_local_emission_exposes_the_genuine_temperature_transpose() {
    let fixture = Fixture::new(false);
    let boundary = boundary(&fixture, false);
    let temperature = [300.0, 310.0, 500.0, 550.0, 300.0, 320.0, 340.0, 360.0];
    let u = [0.2, -0.4, 0.6, 0.3, -0.1, 0.7, -0.2, 0.5];
    let v = [-0.3, 0.8, 0.1, -0.5, 0.6, 0.2, 0.4, -0.7];
    with_cx(|cx| {
        let radiation = PointwiseEnclosureRadiation::new(
            cx, &fixture.mesh, &boundary, &fixture.enclosure,
            PointwiseEnclosureTolerance::default(),
        ).unwrap();
        let assembly = radiation.assemble(cx, &temperature).unwrap();
        let mut jv = [0.0; 8];
        let mut jtu = [0.0; 8];
        let mut wrong = [0.0; 8];
        assembly.tangent_apply(cx, &v, &mut jv).unwrap();
        assembly.tangent_apply_transpose(cx, &u, &mut jtu).unwrap();
        assembly.tangent_apply(cx, &u, &mut wrong).unwrap();
        close(dot(&u, &jv), dot(&jtu, &v), 2e-13);
        assert!(wrong.iter().zip(jtu).any(|(a, b)| (a - b).abs() > 0.1));
        let delta = 1e-3;
        let plus = std::array::from_fn::<_, 8, _>(|i| temperature[i] + delta * v[i]);
        let minus = std::array::from_fn::<_, 8, _>(|i| temperature[i] - delta * v[i]);
        // The existing oracle uses exact simplex moments and an explicit
        // two-plate radiosity solution, independent of production quadrature.
        let plus = oracle(&fixture, &plus).0;
        let minus = oracle(&fixture, &minus).0;
        for i in 0..8 {
            close(jv[i], (plus[i] - minus[i]) / (2.0 * delta), 2e-8);
        }
        with_cancelled_cx(|cancelled| {
            assert!(assembly.tangent_apply_transpose(cancelled, &u, &mut jtu).is_err());
        });
    });
}

#[test]
fn accepted_pointwise_enclosure_reverse_keeps_reflection_phase_and_density_feedback() {
    let fixture = Fixture::new(true);
    let boundary = boundary(&fixture, true);
    with_cx(|cx| {
        let source = ScalarField::nodal("heat density", 8, SOURCE.to_vec()).unwrap();
        let radiation = PointwiseEnclosureRadiation::new(
            cx, &fixture.mesh, &boundary, &fixture.enclosure,
            PointwiseEnclosureTolerance::default(),
        ).unwrap();
        let accepted = endpoint(&fixture, cx, &boundary, &OLD, &SOURCE, EPS, DENSITIES, true);
        let engine = assigned(&fixture, cx, DENSITIES, true);
        let lin = engine.linearize_accepted_with_pointwise_enclosure_radiation(
            cx, problem(&fixture, &boundary, &source), None, &OLD, DT, controls(), &radiation,
            accepted.clone(), controls().newton.absolute_tolerance,
            workspace_entries(8, fixture.mesh.boundary().len(), 2).unwrap(),
        ).unwrap();
        assert_eq!(lin.primal().specific_enthalpy_j_kg, accepted.specific_enthalpy_j_kg);
        assert_eq!(lin.primal().temperature, accepted.temperature);
        assert_eq!(lin.primal().liquid_mass_fraction, accepted.liquid_mass_fraction);
        assert_eq!(lin.primal().pointwise_enclosure_radiation, accepted.pointwise_enclosure_radiation);
        assert_eq!(lin.primal().neumann_out_w, accepted.neumann_out_w);
        assert_eq!(lin.primal().robin_out_w, accepted.robin_out_w);
        assert!(accepted.pointwise_enclosure_radiation.as_ref().unwrap().patches.iter()
            .map(|patch| patch.outward_heat_w.abs()).sum::<f64>() > 1.0);
        let range = accepted.temperature.iter().copied().fold(f64::NEG_INFINITY, f64::max)
            - accepted.temperature.iter().copied().fold(f64::INFINITY, f64::min);
        assert!(range > 80.0);
        let transport = lin.transport();
        let u = [0.2, -0.5, 0.7, 0.3, -0.4, 0.8, -0.1, 0.6];
        let v = [-1.0, 0.1, 0.6, 0.8, 0.4, -0.2, 0.5, -0.3];
        let jv = transport.apply_jacobian(cx, &v).unwrap();
        let jtu = transport.apply_jacobian_transpose(cx, &u).unwrap();
        close(dot(&u, &jv), dot(&jtu, &v), 2e-13);
        let seed = transport.temperature_pullback(cx, &GOAL).unwrap();
        let gradient = lin.pullback(cx, &seed, linear()).unwrap();
        let outward_weights: Vec<_> = transport.nodal_load_pullback(cx, &gradient.transport.adjoint)
            .unwrap().into_iter().map(|x| -x).collect();
        let parameters = transport.pointwise_parameter_pullback(cx, &outward_weights).unwrap();
        assert_eq!(parameters.radiation_regions, ["left", "right"]);
        assert!(parameters.reservoir_temperatures.is_empty());
        assert_eq!(gradient.material_reference_density.len(), 3);
        assert_eq!(gradient.material_reference_density[2], 0.0);
        assert_eq!(&gradient.material_reference_density[..2], gradient.element_reference_density.as_slice());
        assert!(gradient.transport.relative_residual < linear().tolerance);
        let actual_dual = transport.apply_jacobian_transpose(cx, &gradient.transport.adjoint).unwrap();
        assert!(actual_dual.iter().zip(&seed).all(|(a, b)| (a - b).abs() < 2e-12));
        let history_direction = [0.3, -0.7, 0.2, 0.5, -0.1, 0.8, -0.4, 0.6];
        agree(
            dot(&gradient.transport.previous_specific_enthalpy, &history_direction),
            five_point(0.03, |delta| {
                let old = std::array::from_fn::<_, 8, _>(|i| OLD[i] + delta * history_direction[i]);
                objective(&endpoint(&fixture, cx, &boundary, &old, &SOURCE, EPS, DENSITIES, true))
            }),
            "pointwise enclosure history",
        );
        let source_direction = [0.5, 0.2, -0.4, 0.8, 0.6, -0.3, 0.9, 0.1];
        agree(
            dot(&gradient.transport.source_density, &source_direction),
            five_point(0.1, |delta| {
                let q = std::array::from_fn::<_, 8, _>(|i| SOURCE[i] + delta * source_direction[i]);
                objective(&endpoint(&fixture, cx, &boundary, &OLD, &q, EPS, DENSITIES, true))
            }),
            "pointwise enclosure source",
        );
        for patch in 0..2 {
            agree(
                parameters.emissivities[patch],
                five_point(1e-4, |delta| {
                    let mut epsilon = EPS;
                    epsilon[patch] += delta;
                    objective(&endpoint(&fixture, cx, &boundary, &OLD, &SOURCE, epsilon, DENSITIES, true))
                }),
                "pointwise reflected emissivity",
            );
            agree(
                gradient.material_reference_density[patch],
                five_point(0.01, |delta| {
                    let mut densities = DENSITIES;
                    densities[patch] += delta;
                    objective(&endpoint(&fixture, cx, &boundary, &OLD, &SOURCE, EPS, densities, true))
                }),
                "pointwise heterogeneous density",
            );
        }
        let density_direction = [0.3, -0.4];
        let density_action = lin.apply_reference_density_jacobian(cx, &density_direction).unwrap();
        agree(
            -dot(&gradient.transport.adjoint, &density_action),
            dot(&gradient.element_reference_density, &density_direction),
            "pointwise enclosure density residual transpose",
        );
        assert_eq!(lin.pullback(cx, &seed, linear()).unwrap(), gradient);
    });
}

#[test]
fn reflected_radiation_controls_enter_latent_carry_even_when_temperature_seed_is_zero() {
    let fixture = Fixture::new(true);
    let boundary = boundary(&fixture, false);
    let old = [4000.0, 4000.0, 4000.0, 4000.0, 2000.0, 2000.0, 2000.0, 2000.0];
    with_cx(|cx| {
        let source = ScalarField::Uniform(0.0);
        let radiation = PointwiseEnclosureRadiation::new(
            cx, &fixture.mesh, &boundary, &fixture.enclosure,
            PointwiseEnclosureTolerance::default(),
        ).unwrap();
        let accepted = endpoint(&fixture, cx, &boundary, &old, &[0.0; 8], EPS, DENSITIES, false);
        assert_eq!(&accepted.temperature[4..], &[350.0; 4]);
        assert!(accepted.specific_enthalpy_j_kg[4..].iter().all(|h| *h > 2000.0 && *h < 3000.0));
        assert!(accepted.energy_residual_j.abs() <= controls().energy_tolerance_j);
        let engine = assigned(&fixture, cx, DENSITIES, false);
        let lin = engine.linearize_accepted_with_pointwise_enclosure_radiation(
            cx, problem(&fixture, &boundary, &source), None, &old, DT, controls(), &radiation,
            accepted, controls().newton.absolute_tolerance,
            workspace_entries(8, fixture.mesh.boundary().len(), 2).unwrap(),
        ).unwrap();
        let receiver = [0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0];
        let transport = lin.transport();
        let temperature_seed = transport.temperature_pullback(cx, &receiver).unwrap();
        assert_eq!(temperature_seed, [0.0; 8]);
        assert_eq!(lin.pullback(cx, &temperature_seed, linear()).unwrap().transport.adjoint, [0.0; 8]);
        let carry = lin.pullback(cx, &receiver, linear()).unwrap();
        let outward_weights: Vec<_> = transport.nodal_load_pullback(cx, &carry.transport.adjoint)
            .unwrap().into_iter().map(|x| -x).collect();
        let parameters = transport.pointwise_parameter_pullback(cx, &outward_weights).unwrap();
        for patch in 0..2 {
            assert!(parameters.emissivities[patch] > 1.0);
            agree(
                parameters.emissivities[patch],
                five_point(1e-4, |delta| {
                    let mut epsilon = EPS;
                    epsilon[patch] += delta;
                    endpoint(&fixture, cx, &boundary, &old, &[0.0; 8], epsilon, DENSITIES, false)
                        .specific_enthalpy_j_kg[4..].iter().sum()
                }),
                "reflected latent enthalpy carry",
            );
        }
    });
}

#[test]
fn pointwise_enclosure_binding_checks_the_real_endpoint_and_bounded_work() {
    let fixture = Fixture::new(true);
    let boundary = boundary(&fixture, true);
    with_cx(|cx| {
        let source = ScalarField::nodal("heat density", 8, SOURCE.to_vec()).unwrap();
        let radiation = PointwiseEnclosureRadiation::new(
            cx, &fixture.mesh, &boundary, &fixture.enclosure,
            PointwiseEnclosureTolerance::default(),
        ).unwrap();
        let accepted = endpoint(&fixture, cx, &boundary, &OLD, &SOURCE, EPS, DENSITIES, true);
        let engine = assigned(&fixture, cx, DENSITIES, true);
        let entries = workspace_entries(8, fixture.mesh.boundary().len(), 2).unwrap();
        let target = controls().newton.absolute_tolerance;
        let bind = |cx: &fs_exec::Cx<'_>, state: EnthalpyStepSolution, target, entries| {
            engine.linearize_accepted_with_pointwise_enclosure_radiation(
                cx, problem(&fixture, &boundary, &source), None, &OLD, DT, controls(), &radiation,
                state, target, entries,
            )
        };
        assert!(bind(cx, accepted.clone(), target, entries - 1).is_err());
        assert!(bind(cx, accepted.clone(), -1.0, entries).is_err());
        let mut stale = accepted.clone();
        stale.specific_enthalpy_j_kg[0] += 10.0;
        assert!(bind(cx, stale, target, entries).is_err());
        assert!(engine.linearize_accepted(
            cx, problem(&fixture, &boundary, &source), None, &OLD, DT, controls(), accepted.clone(),
        ).is_err());
        let lin = bind(cx, accepted.clone(), target, entries).unwrap();
        let seed = lin.transport().temperature_pullback(cx, &GOAL).unwrap();
        let gradient = lin.pullback(cx, &seed, linear()).unwrap();
        assert!(lin.transport().pointwise_parameter_pullback(cx, &[1.0; 7]).is_err());
        with_cancelled_cx(|cancelled| {
            assert!(bind(cancelled, accepted.clone(), target, entries).is_err());
            assert!(lin.pullback(cancelled, &seed, linear()).is_err());
            assert!(lin.transport().pointwise_parameter_pullback(cancelled, &[1.0; 8]).is_err());
        });
        assert_eq!(lin.primal().specific_enthalpy_j_kg, accepted.specific_enthalpy_j_kg);
        assert_eq!(lin.primal().pointwise_enclosure_radiation, accepted.pointwise_enclosure_radiation);
        assert_eq!(lin.pullback(cx, &seed, linear()).unwrap(), gradient);
    });
}
