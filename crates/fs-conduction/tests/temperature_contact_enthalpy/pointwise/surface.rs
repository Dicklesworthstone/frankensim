//! Source-surface fixed points retain original pointwise radiation and history.
use super::*;
use fs_conduction::interface::temperature::enthalpy::ContactEnthalpyStep;
use fs_conduction::radiation::{
    ViewFactorEvidence, ViewFactorTolerance, open::OpenDiffuseEnclosure,
};

fn surface_card(lower: f64, upper: f64) -> InterfaceSystemCard {
    let dims = fs_conduction::AREA_SPECIFIC_THERMAL_RESISTANCE_DIMS;
    let mut claims = ClaimSet::new();
    claims
        .insert_claim(PropertyClaim {
            key: PropertyKey::new(
                fs_conduction::AREA_SPECIFIC_THERMAL_RESISTANCE_PROPERTY,
                dims,
            ),
            value: PropertyValue::Surface2D {
                axes: ["T".into(), "normal_pressure".into()],
                axis_dims: [
                    fs_conduction::TEMPERATURE_DIMS,
                    fs_qty::Dims([-1, 1, -2, 0, 0, 0]),
                ],
                samples: [
                    [lower, 0.5e6],
                    [upper, 0.5e6],
                    [lower, 1.5e6],
                    [upper, 1.5e6],
                ]
                .into_iter()
                .map(|p| (p, 0.1 + 0.006 * (p[0] - 280.) + 1e-8 * (p[1] - 1e6)))
                .collect(),
                triangles: vec![[0, 1, 2], [1, 3, 2]],
                dims,
            },
            validity: ValidityDomain::unconstrained()
                .with("T", lower, upper)
                .with("normal_pressure", 0.5e6, 1.5e6),
            uncertainty: UncertaintyModel::Unstated,
            interpolation: InterpolationPolicy::LinearInside,
            observations: vec![],
            provenance: Provenance {
                source: "manufactured surface/contact/pointwise-radiation reference".into(),
                license: "internal-test-use".into(),
                artifact: None,
            },
        })
        .unwrap();
    let side = |name: &str| SurfaceSpec {
        material: MaterialStateId {
            chemistry: name.into(),
            phase: "solid".into(),
            process: "fixture".into(),
            revision: 0,
        },
        texture_frame: "fixed normal".into(),
    };
    InterfaceSystemCard::assemble(
        side("A"),
        side("B"),
        SystemContext {
            medium: "dry".into(),
            third_body: None,
            environment: "vacuum".into(),
            history: "unaged".into(),
        },
        claims,
        vec![],
    )
    .unwrap()
}

fn engine<'m, 'c>(
    cx: &Cx<'_>,
    mesh: &'m ConductionMesh,
    cold: &'c EquilibriumEnthalpyPhaseCurve,
    hot: &'c EquilibriumEnthalpyPhaseCurve,
) -> HeterogeneousEnthalpyBackwardEuler<'m, 'c> {
    HeterogeneousEnthalpyBackwardEuler::new(
        cx,
        mesh,
        &[
            ReferenceEnthalpyMaterial {
                curve: cold,
                reference_density_kg_m3: 1.,
            },
            ReferenceEnthalpyMaterial {
                curve: hot,
                reference_density_kg_m3: 2.,
            },
        ],
        &[0, 1],
        EnthalpyBudget {
            max_vertices: 8,
            max_elements: 2,
        },
    )
    .unwrap()
}

fn audit(accepted: &ContactEnthalpyStep, h: &[f64; 8], t: &[f64; 8]) {
    for v in 0..8 {
        close(accepted.conduction.specific_enthalpy_j_kg[v], h[v], 3e-6);
        close(accepted.conduction.temperature[v], t[v], 3e-7);
    }
    assert!(accepted.iterations > 1 && accepted.iterations <= coupling().max_iterations);
    assert!(accepted.nonlinear_iterations > accepted.conduction.newton.iterations);
    assert!(
        accepted.nonlinear_iterations <= coupling().max_iterations * config().max_newton_iterations
    );
    assert!(
        accepted.krylov_iterations
            >= accepted
                .conduction
                .newton
                .history
                .iter()
                .map(|row| row.linear_iterations)
                .sum::<usize>()
    );
    assert!(
        accepted.krylov_iterations
            <= coupling().max_iterations
                * config().max_newton_iterations
                * config().newton.linear_restart
                * config().newton.max_linear_cycles
    );
    assert!(accepted.physical.residual_norm_j <= coupling().residual_tolerance_j);
    close(
        accepted.physical.energy_residual_j,
        0.,
        config().energy_tolerance_j,
    );
    assert!(accepted.radiation.is_none() && accepted.combined_boundary.is_none());
    assert_eq!(accepted.radiation_iterations, 0);
    let state = &accepted.physical_contacts.states()[0];
    close(
        state.resistance.value_m2_k_per_w(),
        0.1 + 0.006 * (state.interface_temperature_k - 280.),
        2e-14,
    );
}

#[test]
fn g1_surface_pointwise_ambient_recovers_simplex_moments_and_latent_heat() {
    context(&CancelGate::new_clock_free(), |cx| {
        let mesh = mesh();
        let boundary = ambient_boundary(&mesh);
        let card = surface_card(280., 380.);
        let contacts =
            TemperatureContacts::new(cx, &mesh, &boundary, vec![], vec![law(&mesh, &card)], 1)
                .unwrap();
        let material = ConductivityModel::isotropic_declared(2.).unwrap();
        let source = ScalarField::Uniform(0.);
        let problem = ConductionProblem {
            mesh: &mesh,
            boundary: &boundary,
            material: &material,
            element_materials: None,
            source: &source,
        };
        for ambient in [300., 340.] {
            let load = ambient_load(ambient);
            let heat: f64 = load.iter().sum();
            let convective = convection_load();
            let original = manufactured_history(std::array::from_fn(|v| load[v] + convective[v]));
            let patches =
                [
                    AmbientRadiationPatch::new("exterior", emissivity("exterior", 0.2), ambient)
                        .unwrap(),
                ];
            let owner = PointwiseAmbientRadiation::new(cx, &mesh, &boundary, &patches).unwrap();
            for offsets in [[0., 0.], [1200., -800.]] {
                let cold = shifted_chart(false, offsets[0]);
                let hot = shifted_chart(true, offsets[1]);
                let engine = engine(cx, &mesh, &cold, &hot);
                let old: Vec<_> = original
                    .iter()
                    .enumerate()
                    .map(|(v, h)| h + offsets[v / 4])
                    .collect();
                let before = old.clone();
                let accepted = contacts
                    .advance_enthalpy_with_pointwise_ambient_radiation(
                        cx,
                        &engine,
                        problem,
                        &old,
                        DT,
                        config(),
                        coupling(),
                        &owner,
                    )
                    .unwrap();
                assert_eq!(old, before);
                audit(
                    &accepted,
                    &std::array::from_fn(|v| TARGET_H[v] + offsets[v / 4]),
                    &TARGET_T,
                );
                close(accepted.conduction.liquid_mass_fraction[0], 0.2, 3e-8);
                close(accepted.conduction.liquid_mass_fraction[5], 0.5, 3e-8);
                close(accepted.physical.convective_out_w, 10., 1e-7);
                close(accepted.physical.radiation_out_w, heat, 1e-7);
                close(
                    accepted.physical.stored_energy_change_j,
                    -DT * (10. + heat),
                    1e-8,
                );
                close(
                    accepted
                        .conduction
                        .pointwise_radiation
                        .as_ref()
                        .unwrap()
                        .outward_heat_w,
                    heat,
                    1e-7,
                );
                assert_eq!(
                    accepted.physical_contacts.states()[0]
                        .resistance
                        .card_identity(),
                    card.content_hash()
                );
            }
        }
    });
}

// Exact local emission moments and an independent two-by-two open radiosity
// solve. The bath is an external reservoir, never a caloric degree of freedom.
fn open_load(epsilon: [f64; 2], view: f64, ambient: f64) -> [f64; 8] {
    let moments = [COLD_FACE, HOT_FACE].map(|face| {
        emission_moments(face.map(|v| TARGET_T[v])).map(|q| q * STEFAN_BOLTZMANN_W_M2_K4)
    });
    let emission = moments.map(|m| m.iter().sum::<f64>() / 0.5);
    let bath = STEFAN_BOLTZMANN_W_M2_K4 * ambient.powi(4);
    let escape = 1. - view;
    let a = std::array::from_fn::<_, 2, _>(|i| {
        epsilon[i] * emission[i] + (1. - epsilon[i]) * escape * bath
    });
    let j0 = (a[0] + (1. - epsilon[0]) * view * a[1])
        / (1. - (1. - epsilon[0]) * (1. - epsilon[1]) * view * view);
    let j1 = a[1] + (1. - epsilon[1]) * view * j0;
    let incoming = [view * j1 + escape * bath, view * j0 + escape * bath];
    let mut load = [0.; 8];
    for (patch, face) in [COLD_FACE, HOT_FACE].into_iter().enumerate() {
        for i in 0..3 {
            load[face[i]] = epsilon[patch] * (moments[patch][i] - 0.5 / 3. * incoming[patch]);
        }
    }
    load
}

#[test]
fn g1_g3_surface_pointwise_enclosures_keep_local_balance_and_one_external_bath() {
    context(&CancelGate::new_clock_free(), |cx| {
        let mesh = mesh();
        let boundary = enclosure_boundary(&mesh);
        let card = surface_card(280., 380.);
        let contacts =
            TemperatureContacts::new(cx, &mesh, &boundary, vec![], vec![law(&mesh, &card)], 1)
                .unwrap();
        let cold = chart(false);
        let hot = chart(true);
        let engine = engine(cx, &mesh, &cold, &hot);
        let material = ConductivityModel::isotropic_declared(2.).unwrap();
        let source = ScalarField::Uniform(0.);
        let problem = ConductionProblem {
            mesh: &mesh,
            boundary: &boundary,
            material: &material,
            element_materials: None,
            source: &source,
        };
        let epsilon = [0.2, 0.7];
        let closed = enclosure(&mesh, epsilon);
        let open = OpenDiffuseEnclosure::new(
            cx,
            closed.surfaces().to_vec(),
            vec![vec![0., 0.35], vec![0.35, 0.]],
            ViewFactorEvidence::Analytic {
                geometry: "prescribed equal-area numerical exchange".into(),
            },
            ViewFactorTolerance::default(),
        )
        .unwrap();
        for is_open in [false, true] {
            let owner = if is_open {
                PointwiseEnclosureRadiation::new_open(
                    cx,
                    &mesh,
                    &boundary,
                    &open,
                    300.,
                    PointwiseEnclosureTolerance::default(),
                )
                .unwrap()
            } else {
                PointwiseEnclosureRadiation::new(
                    cx,
                    &mesh,
                    &boundary,
                    &closed,
                    PointwiseEnclosureTolerance::default(),
                )
                .unwrap()
            };
            let load = if is_open {
                open_load(epsilon, 0.35, 300.)
            } else {
                enclosure_load(epsilon)
            };
            let old = manufactured_history(load);
            let before = old.clone();
            if !is_open {
                let initial: Vec<_> = old
                    .iter()
                    .enumerate()
                    .map(|(v, &h)| {
                        engine
                            .phase_curve_at(v)
                            .unwrap()
                            .state_at_specific_enthalpy(h)
                            .unwrap()
                            .temperature_k()
                    })
                    .collect();
                let frozen = contacts.evaluate(cx, &boundary, &initial).unwrap();
                let trial = engine
                    .advance_with_pointwise_enclosure_radiation(
                        cx,
                        problem,
                        Some(frozen.frozen_interfaces()),
                        &old,
                        DT,
                        config(),
                        &owner,
                    )
                    .unwrap();
                let actual = contacts
                    .evaluate(cx, &boundary, &trial.temperature)
                    .unwrap();
                let physical = engine
                    .evaluate_endpoint_with_pointwise_enclosure_radiation(
                        cx,
                        problem,
                        Some(actual.frozen_interfaces()),
                        &old,
                        &trial.specific_enthalpy_j_kg,
                        DT,
                        config(),
                        &owner,
                    )
                    .unwrap();
                close(physical.energy_residual_j, 0., 1e-8);
                close(physical.radiation_out_w, 0., 1e-8);
                assert!(physical.residual_norm_j > coupling().residual_tolerance_j * 100.);
                assert!(
                    contacts
                        .advance_enthalpy_with_pointwise_enclosure_radiation(
                            cx,
                            &engine,
                            problem,
                            &old,
                            DT,
                            config(),
                            ContactIterationConfig {
                                max_iterations: 1,
                                ..coupling()
                            },
                            &owner
                        )
                        .is_err(),
                    "zero net enclosure heat must not admit an unbalanced contact field"
                );
            }
            let accepted = contacts
                .advance_enthalpy_with_pointwise_enclosure_radiation(
                    cx,
                    &engine,
                    problem,
                    &old,
                    DT,
                    config(),
                    coupling(),
                    &owner,
                )
                .unwrap();
            assert_eq!(old, before);
            audit(&accepted, &TARGET_H, &TARGET_T);
            let report = accepted
                .conduction
                .pointwise_enclosure_radiation
                .as_ref()
                .unwrap();
            let heat: f64 = load.iter().sum();
            close(accepted.physical.radiation_out_w, heat, 2e-7);
            close(accepted.physical.stored_energy_change_j, -DT * heat, 3e-8);
            close(report.outward_heat_w, heat, 2e-7);
            for (patch, face) in [COLD_FACE, HOT_FACE].into_iter().enumerate() {
                close(
                    report.patches[patch].outward_heat_w,
                    face.iter().map(|&v| load[v]).sum(),
                    2e-7,
                );
            }
            if is_open {
                assert!(heat.abs() > 1.);
                close(report.ambient.as_ref().unwrap().outward_heat_w, -heat, 2e-7);
            } else {
                assert!(report.ambient.is_none());
                assert!(report.patches[0].outward_heat_w < -10.);
                close(heat, 0., 1e-12);
            }
        }
    });
}

fn transport_load(t: [f64; 8], mut load: [f64; 8]) -> [f64; 8] {
    // Independent tetrahedral stiffness and consistent triangular interface.
    for start in [0, 4] {
        load[start] += (3. * t[start] - t[start + 1] - t[start + 2] - t[start + 3]) / 3.;
        for local in 1..4 {
            load[start + local] += (t[start + local] - t[start]) / 3.;
        }
    }
    let mean = (0..3).map(|i| t[i] + t[i + 4]).sum::<f64>() / 6.;
    let resistance = 0.1 + 0.006 * (mean - 280.);
    let delta: [f64; 3] = std::array::from_fn(|i| t[i] - t[i + 4]);
    for i in 0..3 {
        let heat = (delta.iter().sum::<f64>() + delta[i]) / (24. * resistance);
        load[i] += heat;
        load[i + 4] -= heat;
    }
    load
}

#[test]
fn g1_surface_pointwise_radiation_keeps_prescribed_node_storage_in_both_jump_directions() {
    context(&CancelGate::new_clock_free(), |cx| {
        let mesh = mesh();
        let boundary = ThermalBoundaryBuilder::new(&mesh)
            .region(
                "bath",
                |f| {
                    f.vertices
                        .iter()
                        .all(|&v| COLD_FACE.contains(&(v as usize)))
                },
                ThermalBc::dirichlet(315.).unwrap(),
            )
            .unwrap()
            .region(
                "radiator",
                |f| f.vertices.iter().all(|&v| HOT_FACE.contains(&(v as usize))),
                ThermalBc::adiabatic(),
            )
            .unwrap()
            .adiabatic_remainder()
            .finish()
            .unwrap();
        let card = surface_card(280., 380.);
        let contacts =
            TemperatureContacts::new(cx, &mesh, &boundary, vec![], vec![law(&mesh, &card)], 1)
                .unwrap();
        let cold = chart(false);
        let hot = chart(true);
        let engine = engine(cx, &mesh, &cold, &hot);
        let material = ConductivityModel::isotropic_declared(2.).unwrap();
        let source = ScalarField::Uniform(0.);
        let problem = ConductionProblem {
            mesh: &mesh,
            boundary: &boundary,
            material: &material,
            element_materials: None,
            source: &source,
        };
        let patches =
            [AmbientRadiationPatch::new("radiator", emissivity("radiator", 0.2), 300.).unwrap()];
        let owner = PointwiseAmbientRadiation::new(cx, &mesh, &boundary, &patches).unwrap();
        let target_t = [315., 315., 310., 315., 340., 345., 360., 350.];
        let target_h = [250., 250., 150., 250., 400., 450., 700., 550.];
        let mut radiation = [0.; 8];
        let moments = emission_moments(HOT_FACE.map(|v| target_t[v]));
        for i in 0..3 {
            radiation[HOT_FACE[i]] =
                0.2 * STEFAN_BOLTZMANN_W_M2_K4 * (moments[i] - 0.5 / 3. * 300_f64.powi(4));
        }
        let heat: f64 = radiation.iter().sum();
        let load = transport_load(target_t, radiation);
        for bath_old_h in [0., 400.] {
            let old: Vec<_> = (0..8)
                .map(|v| {
                    if COLD_FACE.contains(&v) {
                        bath_old_h
                    } else {
                        target_h[v] + DT * load[v] / if v < 4 { 1. / 24. } else { 1. / 12. }
                    }
                })
                .collect();
            let before = old.clone();
            let prescribed_storage = COLD_FACE
                .iter()
                .map(|&v| (target_h[v] - old[v]) / (24. * DT))
                .sum::<f64>();
            assert!(prescribed_storage.abs() > 90.);
            let reaction = prescribed_storage + COLD_FACE.iter().map(|&v| load[v]).sum::<f64>();
            let accepted = contacts
                .advance_enthalpy_with_pointwise_ambient_radiation(
                    cx,
                    &engine,
                    problem,
                    &old,
                    DT,
                    config(),
                    coupling(),
                    &owner,
                )
                .unwrap();
            assert_eq!(old, before);
            audit(&accepted, &target_h, &target_t);
            close(accepted.physical.dirichlet_in_w, reaction, 2e-7);
            close(accepted.physical.radiation_out_w, heat, 2e-7);
            close(
                accepted.physical.stored_energy_change_j,
                DT * (reaction - heat),
                3e-8,
            );
        }
    });
}

#[test]
fn g4_surface_pointwise_late_source_exit_work_and_owner_refusals_preserve_history() {
    context(&CancelGate::new_clock_free(), |cx| {
        let mesh = mesh();
        let boundary = ambient_boundary(&mesh);
        let card = surface_card(280., 380.);
        let contacts =
            TemperatureContacts::new(cx, &mesh, &boundary, vec![], vec![law(&mesh, &card)], 1)
                .unwrap();
        let cold = chart(false);
        let hot = chart(true);
        let engine = engine(cx, &mesh, &cold, &hot);
        let material = ConductivityModel::isotropic_declared(2.).unwrap();
        let source = ScalarField::Uniform(0.);
        let problem = ConductionProblem {
            mesh: &mesh,
            boundary: &boundary,
            material: &material,
            element_materials: None,
            source: &source,
        };
        let radiation = ambient_load(300.);
        let convection = convection_load();
        let old = manufactured_history(std::array::from_fn(|v| radiation[v] + convection[v]));
        let before = old.clone();
        let patches =
            [AmbientRadiationPatch::new("exterior", emissivity("exterior", 0.2), 300.).unwrap()];
        let owner = PointwiseAmbientRadiation::new(cx, &mesh, &boundary, &patches).unwrap();
        let initial: Vec<_> = old
            .iter()
            .enumerate()
            .map(|(v, &h)| {
                engine
                    .phase_curve_at(v)
                    .unwrap()
                    .state_at_specific_enthalpy(h)
                    .unwrap()
                    .temperature_k()
            })
            .collect();
        let frozen = contacts.evaluate(cx, &boundary, &initial).unwrap();
        let trial = engine
            .advance_with_pointwise_ambient_radiation(
                cx,
                problem,
                Some(frozen.frozen_interfaces()),
                &old,
                DT,
                config(),
                &owner,
            )
            .unwrap();
        let actual = contacts
            .evaluate(cx, &boundary, &trial.temperature)
            .unwrap();
        let initial_query = frozen.states()[0].interface_temperature_k;
        let actual_query = actual.states()[0].interface_temperature_k;
        let half_span = (actual_query - initial_query).abs() / 2.;
        assert!(half_span > 1e-4);
        let narrow = surface_card(initial_query - half_span, initial_query + half_span);
        let point = QueryPoint::new()
            .with("T", initial_query)
            .unwrap()
            .with("normal_pressure", 1e6)
            .unwrap();
        let law = TemperatureContact::new(
            "joint",
            ThermalInterfaces::coincident_face_pairs(&mesh).unwrap(),
            &narrow,
            &point,
            "T",
            SelectionPolicy::SingleClaimOnly,
        )
        .unwrap();
        let bounded = TemperatureContacts::new(cx, &mesh, &boundary, vec![], vec![law], 1).unwrap();
        assert!(bounded.evaluate(cx, &boundary, &initial).is_ok());
        assert!(bounded.evaluate(cx, &boundary, &trial.temperature).is_err());
        let refused = bounded.advance_enthalpy_with_pointwise_ambient_radiation(
            cx,
            &engine,
            problem,
            &old,
            DT,
            config(),
            coupling(),
            &owner,
        );
        assert!(
            matches!(refused, Err(EnthalpyError::Conduction(_))),
            "late contact query: {refused:?}"
        );
        assert_eq!(old, before);
        assert!(
            contacts
                .advance_enthalpy_with_pointwise_ambient_radiation(
                    cx,
                    &engine,
                    problem,
                    &old,
                    DT,
                    config(),
                    ContactIterationConfig {
                        max_iterations: usize::MAX,
                        ..coupling()
                    },
                    &owner
                )
                .is_err()
        );
        let other_boundary = ambient_boundary(&mesh);
        let wrong = PointwiseAmbientRadiation::new(cx, &mesh, &other_boundary, &patches).unwrap();
        assert!(matches!(
            contacts.advance_enthalpy_with_pointwise_ambient_radiation(
                cx,
                &engine,
                problem,
                &old,
                DT,
                config(),
                coupling(),
                &wrong
            ),
            Err(EnthalpyError::InvalidInput(
                "pointwise radiation must retain the exact transport mesh and boundary"
            ))
        ));
        let gate = CancelGate::new_clock_free();
        gate.request();
        context(&gate, |cancelled| {
            assert!(matches!(
                contacts.advance_enthalpy_with_pointwise_ambient_radiation(
                    cancelled,
                    &engine,
                    problem,
                    &old,
                    DT,
                    config(),
                    coupling(),
                    &owner
                ),
                Err(EnthalpyError::Conduction(ConductionError::Cancelled { .. }))
            ))
        });
        assert_eq!(old, before);
        let retried = contacts
            .advance_enthalpy_with_pointwise_ambient_radiation(
                cx,
                &engine,
                problem,
                &old,
                DT,
                config(),
                coupling(),
                &owner,
            )
            .unwrap();
        audit(&retried, &TARGET_H, &TARGET_T);
    });
}
