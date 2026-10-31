//! G1/G3 open-cavity emission, reflection, thermal storage and derivatives.
//! View factors are authored algebra fixtures, not inferred from these bodies.
use super::*;
use fs_conduction::radiation::OpenDiffuseEnclosure;

const VIEW: f64 = 0.25;

fn model(fixture: &Fixture, cx: &fs_exec::Cx<'_>, epsilon: [f64; 2]) -> OpenDiffuseEnclosure {
    let surfaces = ["left", "right"]
        .into_iter()
        .enumerate()
        .map(|(body, name)| {
            RadiationSurface::new(
                &fixture.mesh,
                name,
                |face| face.vertices.iter().all(|&v| v as usize / 4 == body),
                emissivity(name, epsilon[body]),
            )
            .unwrap()
        })
        .collect();
    OpenDiffuseEnclosure::new(
        cx,
        surfaces,
        vec![vec![0.0, VIEW], vec![VIEW, 0.0]],
        ViewFactorEvidence::Analytic {
            geometry: "declared-partial-view-algebra-fixture".into(),
        },
        ViewFactorTolerance::default(),
    )
    .unwrap()
}

// The local open model goes out of scope here. The prepared spatial owner must
// retain its factors/evidence without borrowing that temporary declaration.
fn spatial<'a>(
    fixture: &'a Fixture,
    boundary: &'a ThermalBoundary,
    cx: &fs_exec::Cx<'_>,
    epsilon: [f64; 2],
    ambient: f64,
) -> PointwiseEnclosureRadiation<'a> {
    let open = model(fixture, cx, epsilon);
    PointwiseEnclosureRadiation::new_open(
        cx,
        &fixture.mesh,
        boundary,
        &open,
        ambient,
        PointwiseEnclosureTolerance::default(),
    )
    .unwrap()
}

// Independent simplex moments and explicit 2x2 elimination, not production
// quadrature, shifted Planck differences, response factors or heat diagnostics.
fn reference(
    fixture: &Fixture,
    temperature: &[f64],
    epsilon: [f64; 2],
    ambient: f64,
) -> (Vec<f64>, [f64; 2], f64) {
    let mut emission = [0.0; 8];
    let mut weights = [0.0; 8];
    let mut area = [0.0; 2];
    for face in fixture.mesh.boundary() {
        let v = face.vertices.map(|v| v as usize);
        let body = v[0] / 4;
        area[body] += face.area;
        let values = moments(face.area, v.map(|vertex| temperature[vertex]));
        for i in 0..3 {
            emission[v[i]] += STEFAN_BOLTZMANN_W_M2_K4 * values[i];
            weights[v[i]] += face.area / 3.0;
        }
    }
    let mean: [f64; 2] = std::array::from_fn(|body| {
        emission[4 * body..4 * body + 4].iter().sum::<f64>() / area[body]
    });
    let eb = STEFAN_BOLTZMANN_W_M2_K4 * ambient.powi(4);
    let rhs: [f64; 2] = std::array::from_fn(|body| {
        epsilon[body] * mean[body] + (1.0 - epsilon[body]) * (1.0 - VIEW) * eb
    });
    let j0 = (rhs[0] + (1.0 - epsilon[0]) * VIEW * rhs[1])
        / (1.0 - (1.0 - epsilon[0]) * (1.0 - epsilon[1]) * VIEW * VIEW);
    let j1 = rhs[1] + (1.0 - epsilon[1]) * VIEW * j0;
    let irradiation = [VIEW * j1 + (1.0 - VIEW) * eb, VIEW * j0 + (1.0 - VIEW) * eb];
    let load: Vec<_> = (0..8)
        .map(|v| epsilon[v / 4] * (emission[v] - weights[v] * irradiation[v / 4]))
        .collect();
    let heat = std::array::from_fn(|body| load[4 * body..4 * body + 4].iter().sum());
    let ambient_heat = (1.0 - VIEW) * (area[0] * (eb - j0) + area[1] * (eb - j1));
    (load, heat, ambient_heat)
}

#[test]
fn local_hot_spots_keep_open_reflection_with_exact_tangent_transpose_and_bath_derivative() {
    let fixture = Fixture::new(false);
    let boundary = boundary(&fixture, false);
    let temperature = [300.0, 310.0, 500.0, 550.0, 300.0, 320.0, 340.0, 360.0];
    let u = [0.2, -0.4, 0.6, 0.3, -0.1, 0.7, -0.2, 0.5];
    let v = [-0.3, 0.8, 0.1, -0.5, 0.6, 0.2, 0.4, -0.7];
    with_cx(|cx| {
        let owner = spatial(&fixture, &boundary, cx, EPS, 280.0);
        let got = owner.assemble(cx, &temperature).unwrap();
        let (load, heat, bath) = reference(&fixture, &temperature, EPS, 280.0);
        for i in 0..8 {
            close(got.outward_load_w[i], load[i], 3e-13);
        }
        for i in 0..2 {
            close(got.report.patches[i].outward_heat_w, heat[i], 3e-13);
        }
        close(got.report.outward_heat_w, heat.iter().sum(), 3e-13);
        let ambient = got.report.ambient.as_ref().unwrap();
        close(ambient.outward_heat_w, bath, 3e-13);
        assert_eq!(ambient.temperature_k, 280.0);
        close(got.report.system_energy_closure_w, 0.0, 2e-10);
        assert!(got.report.outward_heat_w > 0.0 && ambient.outward_heat_w < 0.0);
        let means = [
            temperature[..4].iter().sum::<f64>() / 4.0,
            temperature[4..].iter().sum::<f64>() / 4.0,
        ];
        let mean_only = model(&fixture, cx, EPS).solve(cx, &means, 280.0).unwrap();
        assert!(got.report.patches[0].outward_heat_w > 1.1 * mean_only.surface_heat_w()[0]);

        let mut ju = [0.0; 8];
        let mut jtv = [0.0; 8];
        got.tangent_apply(cx, &u, &mut ju).unwrap();
        got.tangent_apply_transpose(cx, &v, &mut jtv).unwrap();
        close(dot(&v, &ju), dot(&jtv, &u), 3e-13);
        let delta = 1e-3;
        let plus: [f64; 8] = std::array::from_fn(|i| temperature[i] + delta * u[i]);
        let minus: [f64; 8] = std::array::from_fn(|i| temperature[i] - delta * u[i]);
        let lp = reference(&fixture, &plus, EPS, 280.0).0;
        let lm = reference(&fixture, &minus, EPS, 280.0).0;
        for i in 0..8 {
            close(ju[i], (lp[i] - lm[i]) / (2.0 * delta), 3e-8);
        }
        for i in 0..8 {
            let mut direction = [0.0; 8];
            direction[i] = 1.0;
            let mut action = [0.0; 8];
            got.tangent_apply(cx, &direction, &mut action).unwrap();
            close(got.temperature_diagonal_w_k[i], action[i], 2e-13);
        }
        let parameters = got.parameter_pullback(cx, &v).unwrap();
        assert_eq!(parameters.reservoir_temperatures.len(), 1);
        let objective =
            |epsilon, bath| dot(&v, &reference(&fixture, &temperature, epsilon, bath).0);
        close(
            parameters.reservoir_temperatures[0],
            (objective(EPS, 280.0 + delta) - objective(EPS, 280.0 - delta)) / (2.0 * delta),
            3e-8,
        );
        for i in 0..2 {
            let mut ep = EPS;
            ep[i] += 1e-5;
            let mut em = EPS;
            em[i] -= 1e-5;
            close(
                parameters.emissivities[i],
                (objective(ep, 280.0) - objective(em, 280.0)) / 2e-5,
                3e-8,
            );
        }
        assert_eq!(got.report, owner.assemble(cx, &temperature).unwrap().report);
    });
}

#[test]
fn open_cavity_cools_finite_unanchored_storage_and_preserves_each_nodal_energy_equation() {
    let fixture = Fixture::new(false);
    let source = ScalarField::Uniform(0.0);
    let initial = [300.0, 310.0, 500.0, 550.0, 300.0, 320.0, 340.0, 360.0];
    let old: Vec<_> = initial.iter().map(|t| CP * (t - 200.0)).collect();
    let history = old.clone();
    let volume = 2.0_f64.sqrt() / 12.0;
    let area = 3.0_f64.sqrt();
    with_cx(|cx| {
        for mixed in [false, true] {
            let boundary = boundary(&fixture, mixed);
            let owner = spatial(&fixture, &boundary, cx, EPS, 280.0);
            let got = fixture
                .engine(cx, false)
                .advance_with_pointwise_enclosure_radiation(
                    cx,
                    problem(&fixture, &boundary, &source),
                    None,
                    &old,
                    DT,
                    step_config(),
                    &owner,
                )
                .unwrap();
            let (radiative, _, bath) = reference(&fixture, &got.temperature, EPS, 280.0);
            for body in 0..2 {
                let sum = got.temperature[4 * body..4 * body + 4].iter().sum::<f64>();
                let conductivity = 2.0 + 6.0 * (sum / 4.0 - 240.0) / 320.0;
                for v in 4 * body..4 * body + 4 {
                    let conduction = conductivity * volume * (2.0 * got.temperature[v] - 0.5 * sum);
                    let external = if !mixed {
                        0.0
                    } else if body == 0 {
                        -25.0 * area / 4.0
                    } else {
                        HTC * area / 24.0 * (2.0 * got.temperature[v] + sum - 1800.0)
                    };
                    close(
                        RHO * volume / 4.0 * (got.specific_enthalpy_j_kg[v] - old[v])
                            + DT * (conduction + radiative[v] + external),
                        0.0,
                        3e-8,
                    );
                }
            }
            let prescribed = if mixed { -25.0 * area } else { 0.0 };
            let convection = if mixed {
                HTC * area * (got.temperature[4..].iter().sum::<f64>() / 4.0 - 300.0)
            } else {
                0.0
            };
            close(got.neumann_out_w, prescribed, 1e-13);
            close(got.robin_out_w, convection, 1e-12);
            close(
                got.stored_energy_change_j,
                -DT * (prescribed + convection - bath),
                1e-8,
            );
            close(got.energy_residual_j, 0.0, 1e-8);
            assert!(got.stored_energy_change_j < 0.0);
            let report = got.pointwise_enclosure_radiation.as_ref().unwrap();
            close(report.outward_heat_w, -bath, 3e-13);
            close(report.ambient.as_ref().unwrap().outward_heat_w, bath, 3e-13);
            assert_eq!(old, history);
        }
    });
}

#[test]
fn changing_black_bath_supplies_or_removes_real_latent_heat_without_an_artificial_anchor() {
    let fixture = Fixture::new(true);
    let boundary = boundary(&fixture, false);
    let source = ScalarField::Uniform(0.0);
    let old = [2000.0; 8];
    let mass = RHO * 2.0_f64.sqrt() / 12.0;
    with_cx(|cx| {
        let cool = spatial(&fixture, &boundary, cx, EPS, 300.0);
        let warm = cool.with_ambient_temperature_k(400.0).unwrap();
        for (bath, owner) in [(300.0, &cool), (400.0, &warm)] {
            let (_, heat, ambient) = reference(&fixture, &[350.0; 8], EPS, bath);
            let got = fixture
                .engine(cx, true)
                .advance_with_pointwise_enclosure_radiation(
                    cx,
                    problem(&fixture, &boundary, &source),
                    None,
                    &old,
                    DT,
                    step_config(),
                    owner,
                )
                .unwrap();
            for body in 0..2 {
                let h = 2000.0 - DT * heat[body] / mass;
                for v in 4 * body..4 * body + 4 {
                    close(got.specific_enthalpy_j_kg[v], h, 2e-12);
                    assert_eq!(got.temperature[v], 350.0);
                }
            }
            close(got.stored_energy_change_j, DT * ambient, 1e-8);
            assert_eq!(got.robin_out_w, 0.0);
            assert_eq!(got.neumann_out_w, 0.0);
            close(got.energy_residual_j, 0.0, 1e-8);
        }
        let equilibrated = cool
            .with_ambient_temperature_k(350.0)
            .unwrap()
            .assemble(cx, &[350.0; 8])
            .unwrap();
        assert!(equilibrated.outward_load_w.iter().all(|&q| q == 0.0));
        assert_eq!(equilibrated.report.ambient.unwrap().outward_heat_w, 0.0);
        assert!(cool.with_ambient_temperature_k(f64::NAN).is_err());
        let closed = PointwiseEnclosureRadiation::new(
            cx,
            &fixture.mesh,
            &boundary,
            &fixture.enclosure,
            PointwiseEnclosureTolerance::default(),
        )
        .unwrap();
        assert!(closed.with_ambient_temperature_k(350.0).is_err());
    });
    assert_eq!(old, [2000.0; 8]);
}

#[test]
fn accepted_open_enthalpy_adjoint_differentiates_the_actual_black_bath_control() {
    let fixture = Fixture::new(false);
    let boundary = boundary(&fixture, false);
    let source = ScalarField::Uniform(0.0);
    let old = [
        1000.0, 1200.0, 2800.0, 3400.0, 1100.0, 1200.0, 1300.0, 1400.0,
    ];
    let goal = [0.2, -0.3, 0.4, 0.1, 0.5, -0.2, 0.1, 0.3];
    with_cx(|cx| {
        let owner = spatial(&fixture, &boundary, cx, EPS, 280.0);
        let engine = fixture.engine(cx, false);
        let problem = problem(&fixture, &boundary, &source);
        let accepted = engine
            .advance_with_pointwise_enclosure_radiation(
                cx,
                problem,
                None,
                &old,
                DT,
                step_config(),
                &owner,
            )
            .unwrap();
        let linearization = engine
            .linearize_accepted_with_pointwise_enclosure_radiation(
                cx,
                problem,
                None,
                &old,
                DT,
                step_config(),
                &owner,
                accepted,
                step_config().newton.absolute_tolerance,
                100_000,
            )
            .unwrap();
        let transport = &linearization;
        let seed = transport.temperature_pullback(cx, &goal).unwrap();
        let dual = linearization
            .pullback(
                cx,
                &seed,
                LinearConfig {
                    tolerance: 1e-12,
                    max_iterations: 10000,
                    restart: 8,
                },
            )
            .unwrap();
        let weights: Vec<_> = transport
            .nodal_load_pullback(cx, &dual.adjoint)
            .unwrap()
            .into_iter()
            .map(|value| -value)
            .collect();
        let derivative = transport
            .pointwise_parameter_pullback(cx, &weights)
            .unwrap();
        let objective = |bath| {
            let radiation = owner.with_ambient_temperature_k(bath).unwrap();
            let got = engine
                .advance_with_pointwise_enclosure_radiation(
                    cx,
                    problem,
                    None,
                    &old,
                    DT,
                    step_config(),
                    &radiation,
                )
                .unwrap();
            dot(&goal, &got.temperature)
        };
        close(
            derivative.reservoir_temperatures[0],
            (objective(280.01) - objective(279.99)) / 0.02,
            3e-7,
        );
    });
}

#[test]
fn one_fully_exposed_patch_matches_local_ambient_radiation_and_refuses_invalid_states() {
    use fs_conduction::{AmbientRadiationPatch, PointwiseAmbientRadiation};
    let fixture = Fixture::new(false);
    let boundary = boundary(&fixture, false);
    let temperature = [300.0, 310.0, 500.0, 550.0, 300.0, 320.0, 340.0, 360.0];
    with_cx(|cx| {
        let open = OpenDiffuseEnclosure::new(
            cx,
            vec![fixture.enclosure.surfaces()[0].clone()],
            vec![vec![0.0]],
            ViewFactorEvidence::Analytic {
                geometry: "fully-exposed-patch".into(),
            },
            ViewFactorTolerance::default(),
        )
        .unwrap();
        let owner = PointwiseEnclosureRadiation::new_open(
            cx,
            &fixture.mesh,
            &boundary,
            &open,
            3.0,
            PointwiseEnclosureTolerance::default(),
        )
        .unwrap();
        let patches =
            [AmbientRadiationPatch::new("left", emissivity("left", EPS[0]), 3.0).unwrap()];
        let ambient =
            PointwiseAmbientRadiation::new(cx, &fixture.mesh, &boundary, &patches).unwrap();
        let a = owner.assemble(cx, &temperature).unwrap();
        let b = ambient.assemble(cx, &temperature).unwrap();
        for i in 0..8 {
            close(a.outward_load_w[i], b.outward_load_w[i], 3e-13);
        }
        let mut invalid = temperature;
        invalid[0] = 600.0;
        assert!(owner.assemble(cx, &invalid).is_err());
        with_cancelled_cx(|cancelled| {
            assert!(matches!(
                owner.assemble(cancelled, &temperature),
                Err(ConductionError::Cancelled { .. })
            ));
            assert!(matches!(
                PointwiseEnclosureRadiation::new_open(
                    cancelled,
                    &fixture.mesh,
                    &boundary,
                    &open,
                    280.0,
                    PointwiseEnclosureTolerance::default()
                ),
                Err(ConductionError::Cancelled { .. })
            ));
        });
    });
}
