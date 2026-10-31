//! G1/G3 selective surfaces, independent spectral/space integration, finite
//! latent storage and the complete accepted endpoint bath sensitivity.
use super::*;
use fs_conduction::radiation::OpenDiffuseEnclosure;
use fs_conduction::{
    BlackbodyBand, PointwisePatchEmissivity, SpectralPointwiseControls, SpectralSurfaceEmissivity,
    SurfaceSpectralBand, SurfaceSpectralBands,
};

const WINDOWS: [[f64; 2]; 2] = [[0.05e-6, 8.0e-6], [8.0e-6, 2.0e-3]];
const SELECTIVE: [[f64; 2]; 2] = [[0.12, 0.82], [0.91, 0.18]];
const VIEW: f64 = 0.25;

fn spectrum(values: [[f64; 2]; 2]) -> SurfaceSpectralBands {
    SurfaceSpectralBands::new(
        WINDOWS
            .into_iter()
            .zip(values)
            .map(|(wavelengths_m, eps)| SurfaceSpectralBand {
                wavelengths_m,
                surfaces: eps
                    .into_iter()
                    .enumerate()
                    .map(|(i, value)| {
                        SpectralSurfaceEmissivity::new(
                            value,
                            [200.0, 1600.0],
                            format!("authored opaque diffuse band fixture, surface {i}"),
                        )
                        .unwrap()
                    })
                    .collect(),
            })
            .collect(),
    )
    .unwrap()
}

fn controls() -> SpectralPointwiseControls {
    SpectralPointwiseControls {
        absolute_spatial_power_w: 1e-9,
        absolute_spatial_tangent_w_k: 1e-11,
        relative_spatial_tolerance: 2e-10,
        ..SpectralPointwiseControls::default()
    }
}

fn open(fixture: &Fixture, cx: &fs_exec::Cx<'_>) -> OpenDiffuseEnclosure {
    OpenDiffuseEnclosure::new(
        cx,
        fixture.enclosure.surfaces().to_vec(),
        vec![vec![0.0, VIEW], vec![VIEW, 0.0]],
        ViewFactorEvidence::Analytic {
            geometry: "authored partial-view spectral fixture".into(),
        },
        ViewFactorTolerance::default(),
    )
    .unwrap()
}

fn owner<'a>(
    fixture: &'a Fixture,
    boundary: &'a ThermalBoundary,
    cx: &fs_exec::Cx<'_>,
    values: [[f64; 2]; 2],
    bath: f64,
) -> PointwiseEnclosureRadiation<'a> {
    PointwiseEnclosureRadiation::new_open_spectral(
        cx,
        &fixture.mesh,
        boundary,
        &open(fixture, cx),
        bath,
        &spectrum(values),
        controls(),
        PointwiseEnclosureTolerance::default(),
    )
    .unwrap()
}

// Independent SI Planck quadrature in log wavelength, not a production
// moment/tail formula or its cancellation-stable band difference.
fn spectral_integral(window: [f64; 2], temperature: f64) -> f64 {
    const H: f64 = 6.626_070_15e-34;
    const C: f64 = 299_792_458.0;
    const K: f64 = 1.380_649e-23;
    let (a, b) = (window[0].ln(), window[1].ln());
    let panels = 4096;
    let step = (b - a) / panels as f64;
    let mut sum = 0.0;
    for i in 0..=panels {
        let wavelength = (a + step * i as f64).exp();
        let x = H * C / (wavelength * K * temperature);
        let value = if x > 700.0 {
            0.0
        } else {
            2.0 * std::f64::consts::PI * H * C * C / wavelength.powi(4) / x.exp_m1()
        };
        sum += value
            * if i == 0 || i == panels {
                1.0
            } else if i % 2 == 0 {
                2.0
            } else {
                4.0
            };
    }
    step * sum / 3.0
}

#[test]
fn selective_closed_plates_match_independent_spectral_resistances_without_using_gray_cards() {
    let fixture = Fixture::new(false);
    let boundary = boundary(&fixture, false);
    with_cx(|cx| {
        let data = spectrum(SELECTIVE);
        let owner = PointwiseEnclosureRadiation::new_spectral(
            cx,
            &fixture.mesh,
            &boundary,
            &fixture.enclosure,
            &data,
            controls(),
            PointwiseEnclosureTolerance::default(),
        )
        .unwrap();
        // 1000 K is outside the geometry carrier's gray-card support. The
        // independently supplied spectral source is the actual optical owner.
        let got = owner
            .assemble(
                cx,
                &[1000.0, 1000.0, 1000.0, 1000.0, 300.0, 300.0, 300.0, 300.0],
            )
            .unwrap();
        let area = fixture.enclosure.surfaces()[0].area_m2();
        let mut total = 0.0;
        for band in 0..2 {
            let [left, right] = SELECTIVE[band];
            let heat = area
                * (spectral_integral(WINDOWS[band], 1000.0)
                    - spectral_integral(WINDOWS[band], 300.0))
                / (1.0 / left + 1.0 / right - 1.0);
            let report = &got.report.spectral.as_ref().unwrap().bands[band];
            close(report.patches[0].outward_heat_w, heat, 2e-9);
            close(report.patches[1].outward_heat_w, -heat, 2e-9);
            close(report.system_energy_closure_w, 0.0, 1e-8);
            total += heat;
        }
        close(got.report.patches[0].outward_heat_w, total, 2e-9);
        for v in 0..4 {
            close(got.outward_load_w[v], total / 4.0, 2e-9);
        }
        close(got.report.outward_heat_w, 0.0, 1e-8);
        assert!(
            got.report
                .patches
                .iter()
                .all(|p| p.emissivity.gray().is_none())
        );
        assert!(
            matches!(&got.report.patches[0].emissivity, PointwisePatchEmissivity::Spectral(v) if v.len() == 2)
        );
        assert!(
            got.report
                .spectral
                .as_ref()
                .unwrap()
                .maximum_omitted_blackbody_fraction
                < 1e-4
        );
        assert_eq!(
            got.report,
            owner
                .assemble(
                    cx,
                    &[1000.0, 1000.0, 1000.0, 1000.0, 300.0, 300.0, 300.0, 300.0]
                )
                .unwrap()
                .report
        );
    });
}

// Independent product Gauss-Legendre/Duffy surface integration. The Planck
// constitutive primitive is separately checked above against direct SI light.
fn reference(
    fixture: &Fixture,
    temperature: &[f64],
    values: [[f64; 2]; 2],
    bath: f64,
) -> (Vec<f64>, f64) {
    const GAUSS: [(f64, f64); 8] = [
        (-0.960_289_856_497_536_3, 0.101_228_536_290_376_3),
        (-0.796_666_477_413_626_7, 0.222_381_034_453_374_5),
        (-0.525_532_409_916_329, 0.313_706_645_877_887_3),
        (-0.183_434_642_495_649_8, 0.362_683_783_378_362),
        (0.183_434_642_495_649_8, 0.362_683_783_378_362),
        (0.525_532_409_916_329, 0.313_706_645_877_887_3),
        (0.796_666_477_413_626_7, 0.222_381_034_453_374_5),
        (0.960_289_856_497_536_3, 0.101_228_536_290_376_3),
    ];
    let mut load = vec![0.0; 8];
    let mut bath_heat = 0.0;
    for band in 0..2 {
        let law = BlackbodyBand::new(WINDOWS[band][0], WINDOWS[band][1]).unwrap();
        let mut emission = [0.0; 8];
        let mut weights = [0.0; 8];
        let mut area = [0.0; 2];
        for face in fixture.mesh.boundary() {
            let v = face.vertices.map(|v| v as usize);
            area[v[0] / 4] += face.area;
            for i in 0..3 {
                weights[v[i]] += face.area / 3.0;
            }
            for (x, wx) in GAUSS {
                for (y, wy) in GAUSS {
                    let (u, w) = ((1.0 + x) / 2.0, (1.0 + y) / 2.0);
                    let bary = [1.0 - u, u * (1.0 - w), u * w];
                    let t = (0..3).map(|i| bary[i] * temperature[v[i]]).sum::<f64>();
                    let e = law.evaluate(t).unwrap().power_w_m2 * face.area * u * wx * wy / 2.0;
                    for i in 0..3 {
                        emission[v[i]] += bary[i] * e;
                    }
                }
            }
        }
        let mean: [f64; 2] =
            std::array::from_fn(|i| emission[4 * i..4 * i + 4].iter().sum::<f64>() / area[i]);
        let epsilon = values[band];
        let eb = law.evaluate(bath).unwrap().power_w_m2;
        let rhs: [f64; 2] =
            std::array::from_fn(|i| epsilon[i] * mean[i] + (1.0 - epsilon[i]) * (1.0 - VIEW) * eb);
        let j0 = (rhs[0] + (1.0 - epsilon[0]) * VIEW * rhs[1])
            / (1.0 - (1.0 - epsilon[0]) * (1.0 - epsilon[1]) * VIEW * VIEW);
        let j1 = rhs[1] + (1.0 - epsilon[1]) * VIEW * j0;
        let irradiation = [VIEW * j1 + (1.0 - VIEW) * eb, VIEW * j0 + (1.0 - VIEW) * eb];
        for v in 0..8 {
            load[v] += epsilon[v / 4] * (emission[v] - weights[v] * irradiation[v / 4]);
        }
        bath_heat += (1.0 - VIEW) * (area[0] * (eb - j0) + area[1] * (eb - j1));
    }
    (load, bath_heat)
}

#[test]
fn nonuniform_spectral_open_exchange_has_consistent_tangent_transpose_and_band_controls() {
    let fixture = Fixture::new(false);
    let boundary = boundary(&fixture, false);
    let t = [300.0, 320.0, 480.0, 540.0, 310.0, 340.0, 370.0, 420.0];
    let u = [0.2, -0.4, 0.6, 0.3, -0.1, 0.7, -0.2, 0.5];
    let v = [-0.3, 0.8, 0.1, -0.5, 0.6, 0.2, 0.4, -0.7];
    with_cx(|cx| {
        let radiation = owner(&fixture, &boundary, cx, SELECTIVE, 290.0);
        let got = radiation.assemble(cx, &t).unwrap();
        let (expected, bath) = reference(&fixture, &t, SELECTIVE, 290.0);
        for i in 0..8 {
            close(got.outward_load_w[i], expected[i], 4e-8);
        }
        close(
            got.report.ambient.as_ref().unwrap().outward_heat_w,
            bath,
            4e-8,
        );
        close(got.report.system_energy_closure_w, 0.0, 1e-9);
        let mut ju = [0.0; 8];
        let mut jtv = [0.0; 8];
        got.tangent_apply(cx, &u, &mut ju).unwrap();
        got.tangent_apply_transpose(cx, &v, &mut jtv).unwrap();
        close(dot(&v, &ju), dot(&jtv, &u), 2e-12);
        let delta = 1e-3;
        let plus = radiation
            .assemble(cx, &std::array::from_fn::<_, 8, _>(|i| t[i] + delta * u[i]))
            .unwrap();
        let minus = radiation
            .assemble(cx, &std::array::from_fn::<_, 8, _>(|i| t[i] - delta * u[i]))
            .unwrap();
        for i in 0..8 {
            close(
                ju[i],
                (plus.outward_load_w[i] - minus.outward_load_w[i]) / (2.0 * delta),
                2e-7,
            );
        }
        let gradients = got.parameter_pullback(cx, &v).unwrap();
        assert!(gradients.emissivities.is_empty());
        assert_eq!(gradients.spectral_emissivities.len(), 2);
        assert_eq!(gradients.reservoir_temperatures.len(), 1);
        for band in 0..2 {
            for patch in 0..2 {
                let objective = |epsilon| {
                    let mut values = SELECTIVE;
                    values[band][patch] = epsilon;
                    dot(
                        &v,
                        &owner(&fixture, &boundary, cx, values, 290.0)
                            .assemble(cx, &t)
                            .unwrap()
                            .outward_load_w,
                    )
                };
                let d = 1e-5;
                close(
                    gradients.spectral_emissivities[band].emissivities[patch],
                    (objective(SELECTIVE[band][patch] + d) - objective(SELECTIVE[band][patch] - d))
                        / (2.0 * d),
                    2e-7,
                );
            }
        }
        let plus = radiation
            .with_ambient_temperature_k(290.001)
            .unwrap()
            .assemble(cx, &t)
            .unwrap();
        let minus = radiation
            .with_ambient_temperature_k(289.999)
            .unwrap()
            .assemble(cx, &t)
            .unwrap();
        close(
            gradients.reservoir_temperatures[0],
            (dot(&v, &plus.outward_load_w) - dot(&v, &minus.outward_load_w)) / 0.002,
            2e-7,
        );
        let spatial = got.report.spectral.as_ref().unwrap().spatial;
        assert!(spatial.planck_evaluations >= 35 * 2 * fixture.mesh.boundary().len());
        assert!(spatial.estimated_emission_error_w.is_finite());
    });
}

#[test]
fn spectral_load_enters_each_original_enthalpy_equation_and_the_accepted_bath_adjoint() {
    let fixture = Fixture::new(false);
    let boundary = boundary(&fixture, true);
    let source = ScalarField::Uniform(0.0);
    let initial = [300.0, 320.0, 480.0, 540.0, 310.0, 340.0, 370.0, 420.0];
    let old: Vec<_> = initial.iter().map(|t| CP * (t - 200.0)).collect();
    let history = old.clone();
    let volume = 2.0_f64.sqrt() / 12.0;
    let area = 3.0_f64.sqrt();
    let goal = [0.2, -0.3, 0.4, 0.1, 0.5, -0.2, 0.1, 0.3];
    with_cx(|cx| {
        let radiation = owner(&fixture, &boundary, cx, SELECTIVE, 290.0);
        let engine = fixture.engine(cx, false);
        let problem = problem(&fixture, &boundary, &source);
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
        let (heat, bath) = reference(&fixture, &got.temperature, SELECTIVE, 290.0);
        for body in 0..2 {
            let sum = got.temperature[4 * body..4 * body + 4].iter().sum::<f64>();
            let conductivity = 2.0 + 6.0 * (sum / 4.0 - 240.0) / 320.0;
            for v in 4 * body..4 * body + 4 {
                let conduction = conductivity * volume * (2.0 * got.temperature[v] - 0.5 * sum);
                let external = if body == 0 {
                    -25.0 * area / 4.0
                } else {
                    HTC * area / 24.0 * (2.0 * got.temperature[v] + sum - 1800.0)
                };
                close(
                    RHO * volume / 4.0 * (got.specific_enthalpy_j_kg[v] - old[v])
                        + DT * (conduction + heat[v] + external),
                    0.0,
                    2e-6,
                );
            }
        }
        close(got.neumann_out_w, -25.0 * area, 1e-13);
        close(
            got.robin_out_w,
            HTC * area * (got.temperature[4..].iter().sum::<f64>() / 4.0 - 300.0),
            1e-12,
        );
        close(
            got.stored_energy_change_j,
            -DT * (got.neumann_out_w + got.robin_out_w - bath),
            2e-6,
        );
        close(got.energy_residual_j, 0.0, 1e-8);
        let linear = engine
            .linearize_accepted_with_pointwise_enclosure_radiation(
                cx,
                problem,
                None,
                &old,
                DT,
                step_config(),
                &radiation,
                got,
                step_config().newton.absolute_tolerance,
                100_000,
            )
            .unwrap();
        let seed = linear.temperature_pullback(cx, &goal).unwrap();
        let dual = linear
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
        let weights: Vec<_> = linear
            .nodal_load_pullback(cx, &dual.adjoint)
            .unwrap()
            .into_iter()
            .map(|x| -x)
            .collect();
        let derivative = linear.pointwise_parameter_pullback(cx, &weights).unwrap();
        let objective = |bath| {
            let owner = radiation.with_ambient_temperature_k(bath).unwrap();
            let next = engine
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
            dot(&goal, &next.temperature)
        };
        close(
            derivative.reservoir_temperatures[0],
            (objective(290.01) - objective(289.99)) / 0.02,
            4e-7,
        );
        assert_eq!(old, history);
    });
}

#[test]
fn spectral_bath_transfers_real_latent_heat_and_respects_source_coverage_and_work_limits() {
    let fixture = Fixture::new(true);
    let boundary = boundary(&fixture, false);
    let source = ScalarField::Uniform(0.0);
    let old = [2000.0; 8];
    let mass = RHO * 2.0_f64.sqrt() / 12.0;
    with_cx(|cx| {
        let radiation = owner(&fixture, &boundary, cx, SELECTIVE, 300.0);
        for bath in [300.0, 400.0] {
            let owner = radiation.with_ambient_temperature_k(bath).unwrap();
            let (load, ambient) = reference(&fixture, &[350.0; 8], SELECTIVE, bath);
            let got = fixture
                .engine(cx, true)
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
            for body in 0..2 {
                let h = 2000.0 - DT * load[4 * body..4 * body + 4].iter().sum::<f64>() / mass;
                for v in 4 * body..4 * body + 4 {
                    close(got.specific_enthalpy_j_kg[v], h, 2e-11);
                    assert_eq!(got.temperature[v], 350.0);
                }
            }
            close(got.stored_energy_change_j, DT * ambient, 1e-8);
        }
        let equilibrium = radiation
            .with_ambient_temperature_k(350.0)
            .unwrap()
            .assemble(cx, &[350.0; 8])
            .unwrap();
        assert!(equilibrium.outward_load_w.iter().all(|&q| q == 0.0));
        assert!(radiation.assemble(cx, &[1700.0; 8]).is_err());
        let small = PointwiseEnclosureRadiation::new_open_spectral(
            cx,
            &fixture.mesh,
            &boundary,
            &open(&fixture, cx),
            300.0,
            &spectrum(SELECTIVE),
            SpectralPointwiseControls {
                maximum_spatial_evaluations: 35,
                ..controls()
            },
            PointwiseEnclosureTolerance::default(),
        )
        .unwrap();
        assert!(small.assemble(cx, &[350.0; 8]).is_err());
        let visible = SurfaceSpectralBands::new(vec![SurfaceSpectralBand {
            wavelengths_m: [400e-9, 700e-9],
            surfaces: spectrum(SELECTIVE).bands()[0].surfaces.clone(),
        }])
        .unwrap();
        let uncovered = PointwiseEnclosureRadiation::new_open_spectral(
            cx,
            &fixture.mesh,
            &boundary,
            &open(&fixture, cx),
            300.0,
            &visible,
            controls(),
            PointwiseEnclosureTolerance::default(),
        )
        .unwrap();
        assert!(uncovered.assemble(cx, &[350.0; 8]).is_err());
        let reflected = owner(&fixture, &boundary, cx, [[0.0; 2]; 2], 300.0)
            .assemble(cx, &[400.0; 8])
            .unwrap();
        assert!(reflected.outward_load_w.iter().all(|&q| q == 0.0));
        assert!(
            PointwiseEnclosureRadiation::new_spectral(
                cx,
                &fixture.mesh,
                &boundary,
                &fixture.enclosure,
                &spectrum([[0.0; 2]; 2]),
                controls(),
                PointwiseEnclosureTolerance::default()
            )
            .is_err()
        );
    });
    with_cancelled_cx(|cx| {
        assert!(matches!(
            PointwiseEnclosureRadiation::new_spectral(
                cx,
                &fixture.mesh,
                &boundary,
                &fixture.enclosure,
                &spectrum(SELECTIVE),
                controls(),
                PointwiseEnclosureTolerance::default()
            ),
            Err(ConductionError::Cancelled { .. })
        ));
    });
    assert_eq!(old, [2000.0; 8]);
}
