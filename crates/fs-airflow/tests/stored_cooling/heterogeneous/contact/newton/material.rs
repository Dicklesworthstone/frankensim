//! Original emissivity feedback participates in the finite coolant transaction.
use super::*;
use fs_conduction::{MaterialAmbientRadiationPatch, MaterialEnclosureRadiationPatch,
    SurfaceEmissivity};
use fs_material::thermal::source::ThermalMaterialPoint;

fn emissive_card(scale: f64) -> MaterialCard {
    let base = solid_card(scale);
    let mut claims = ClaimSet::new();
    for (_, claim) in base.claims().claims_ordered() {
        claims.insert_claim(claim.clone()).unwrap();
    }
    claims
        .insert_claim(PropertyClaim {
            key: PropertyKey::new(SURFACE_EMISSIVITY_PROPERTY, Dims::NONE),
            value: PropertyValue::Curve {
                abscissa: "T".into(),
                abscissa_dims: Temperature::DIMS,
                knots: vec![(250.0, 0.2), (450.0, 0.8)],
                dims: Dims::NONE,
            },
            validity: ValidityDomain::unconstrained()
                .with("T", 250.0, 450.0)
                .with("pressure", 100000.0, 100000.0),
            interpolation: InterpolationPolicy::LinearInside,
            uncertainty: UncertaintyModel::Unstated,
            observations: Vec::new(),
            provenance: Provenance {
                source: "synthetic temperature-dependent finish".into(),
                license: "test-only".into(),
                artifact: None,
            },
        })
        .unwrap();
    MaterialCard::assemble(base.id().clone(), claims, Vec::new()).unwrap()
}

// Exact P1 mean(T^p) from the same simplex monomial identity as emission().
// This is independent of the production degree-six face quadrature.
fn mean_power(t: [f64; 3], degree: i32) -> f64 {
    let mut sum = 0.0;
    for a in 0..=degree {
        for b in 0..=degree - a {
            sum += t[0].powi(a) * t[1].powi(b) * t[2].powi(degree - a - b);
        }
    }
    2.0 * sum / f64::from((degree + 1) * (degree + 2))
}

#[test]
fn material_emissivity_refreshes_on_the_coupled_endpoint_without_heating_coolant_twice() {
    fixture(1.0, |cx, _, network, pump, uniform, _, material| {
        let (mesh, ids) = adjacent_bodies();
        let cards = [emissive_card(1.0), emissive_card(2.0)];
        let charts = cards
            .each_ref()
            .map(|card| caloric(card, SolidLiquidPhase::Solid));
        let storage = prepare(cx, &mesh, &ids, &charts);
        let source = smooth_contact();
        let contacts = contacts(cx, &mesh, &source);
        let old = unequal_history(&storage, &charts);
        let tank_h = uniform.reservoir().state().specific_enthalpy_j_kg();
        let mut state = HeterogeneousStoredCoolingState::try_new(
            cx,
            uniform.reservoir().clone(),
            &storage,
            old.clone(),
        )
        .unwrap();
        let context = point(300.0);
        let patches: Vec<_> = NAMES
            .into_iter()
            .enumerate()
            .map(|(index, region)| MaterialAmbientRadiationPatch {
                region,
                ambient_temperature_k: 250.0,
                curve: charts[index].curve(),
                material: ThermalMaterialPoint {
                    material: &cards[index],
                    context: &context,
                    temperature_axis: "T",
                    liquid_fraction_axis: None,
                },
                selection: MaterialPropertySelection::SingleClaimOnly,
            })
            .collect();
        let heat = left_heat(&mesh, 20.0);
        let result = network
            .advance_stored_solid_with_pump_and_contact_newton(
                cx,
                &mut state,
                DT,
                0.0,
                &[320.0],
                pump,
                StoredSolidProblem {
                    mesh: &mesh,
                    material,
                    source: &heat,
                    element_materials: None,
                    interfaces: None,
                },
                &contacts,
                StoredContactRadiation::MaterialAmbient(&patches),
                newton_controls(),
                |_, rows| radiating_boundary(&mesh, rows),
            )
            .unwrap();
        let response = &result.step.liquid.circuit.state.solid.state;
        let solid = &response.thermal.conduction;
        let report = solid.pointwise_radiation.as_ref().unwrap();
        let expected: f64 = mesh
            .boundary()
            .iter()
            .filter(|face| face.outward_normal[2] < -0.99)
            .map(|face| {
                let t = face
                    .vertices
                    .map(|vertex| solid.temperature[vertex as usize]);
                let (t1, t4, t5) = (mean_power(t, 1), mean_power(t, 4), mean_power(t, 5));
                // epsilon(T)=0.2+0.003(T-250), evaluated on the full P1 trace.
                face.area
                    * STEFAN_BOLTZMANN_W_M2_K4
                    * (0.2 * (t4 - 250.0_f64.powi(4))
                        + 0.003 * (t5 - 250.0 * t4 - 250.0_f64.powi(4) * (t1 - 250.0)))
            })
            .sum();
        let frozen: f64 = (0..2)
            .map(|body| {
                let (area, emitted) = emission(&mesh, &solid.temperature, body);
                area * 0.35 * (emitted - STEFAN_BOLTZMANN_W_M2_K4 * 250.0_f64.powi(4))
            })
            .sum();
        assert!((response.thermal.physical.radiation_out_w - expected).abs() < 1e-10);
        assert!((report.outward_heat_w - expected).abs() < 1e-10);
        assert!(
            (expected - frozen).abs() > 1e-3,
            "material feedback collapsed to its 300 K snapshot"
        );
        for (index, patch) in report.patches.iter().enumerate() {
            assert_eq!(patch.patch.emissivity().temperature_k(), 300.0);
            assert!(!patch.material_states.is_empty());
            for sample in &patch.material_states {
                let t = solid.temperature[sample.vertex];
                assert_eq!(sample.phase.temperature_k(), t);
                assert_eq!(
                    sample.emissivity.card_identity(),
                    cards[index].content_hash()
                );
                assert!((sample.emissivity.value() - (0.2 + 0.003 * (t - 250.0))).abs() < 1e-13);
            }
        }
        assert_eq!(response.external_convection_out_w, 0.0);
        assert_eq!(response.coolant_robin_fluxes.len(), 1);
        assert!((DT * response.coolant_out_w - result.step.liquid.returned_heat_j).abs() < 2e-6);
        assert!((stored_change(&state, &old, tank_h) - DT * (20.0 - expected)).abs() < 2e-6);
        assert_eq!(
            state.solid_specific_enthalpy_j_kg(),
            solid.specific_enthalpy_j_kg.as_slice()
        );
        assert_eq!(
            result.total_enthalpy_solves,
            result.step.liquid.total_solid_solves
        );
    });
}

#[test]
fn changing_material_reflection_retains_internal_exchange_and_signed_bath_heat() {
    for open in [false, true] {
        fixture(1.0, |cx, _, network, pump, uniform, _, material| {
            let (mesh, ids) = adjacent_bodies();
            let cards = [emissive_card(1.0), emissive_card(2.0)];
            let charts = cards.each_ref().map(|card| caloric(card, SolidLiquidPhase::Solid));
            let storage = prepare(cx, &mesh, &ids, &charts);
            let source = smooth_contact();
            let contacts = contacts(cx, &mesh, &source);
            let old = unequal_history(&storage, &charts);
            let tank_h = uniform.reservoir().state().specific_enthalpy_j_kg();
            let mut state = HeterogeneousStoredCoolingState::try_new(
                cx, uniform.reservoir().clone(), &storage, old.clone(),
            ).unwrap();
            let context = point(300.0);
            let patches: Vec<_> = NAMES.into_iter().enumerate().map(|(index, region)| {
                MaterialEnclosureRadiationPatch {
                    region, curve: charts[index].curve(),
                    material: ThermalMaterialPoint {
                        material: &cards[index], context: &context,
                        temperature_axis: "T", liquid_fraction_axis: None,
                    },
                    selection: MaterialPropertySelection::SingleClaimOnly,
                }
            }).collect();
            let surfaces = || patches.iter().enumerate().map(|(body, patch)| {
                RadiationSurface::new(
                    &mesh, patch.region,
                    |face| face.outward_normal[2] < -0.99
                        && face.vertices.iter().all(|&v| v as usize / 8 == body),
                    SurfaceEmissivity::from_condition(patch.region, patch.material,
                        patch.selection.clone()).unwrap().0,
                ).unwrap()
            }).collect();
            let closed = GrayDiffuseEnclosure::new(surfaces(),
                ViewFactorMatrix::infinite_parallel_plates(0.003).unwrap()).unwrap();
            let partial = OpenDiffuseEnclosure::new(
                cx, surfaces(), vec![vec![0.0, 0.25], vec![0.25, 0.0]],
                ViewFactorEvidence::Analytic { geometry: "declared partial-view material fixture".into() },
                ViewFactorTolerance::default(),
            ).unwrap();
            let radiation = if open {
                StoredContactRadiation::MaterialOpenEnclosure {
                    enclosure: &partial, patches: &patches, ambient_temperature_k: 420.0,
                    tolerance: PointwiseEnclosureTolerance::default(),
                }
            } else {
                StoredContactRadiation::MaterialClosedEnclosure {
                    enclosure: &closed, patches: &patches,
                    tolerance: PointwiseEnclosureTolerance::default(),
                }
            };
            let heat = left_heat(&mesh, 20.0);
            let result = network.advance_stored_solid_with_pump_and_contact_newton(
                cx, &mut state, DT, 0.0, &[320.0], pump,
                StoredSolidProblem { mesh: &mesh, material, source: &heat,
                    element_materials: None, interfaces: None },
                &contacts, radiation, newton_controls(),
                |_, rows| radiating_boundary(&mesh, rows),
            ).unwrap();
            let response = &result.step.liquid.circuit.state.solid.state;
            let solid = &response.thermal.conduction;
            let report = solid.pointwise_enclosure_radiation.as_ref().unwrap();
            // Independent exact P1 moments retain the correlation epsilon(T)*T^4.
            let moments: [(f64, f64, f64); 2] = std::array::from_fn(|body| {
                let (mut area, mut epsilon, mut emission) = (0.0, 0.0, 0.0);
                for face in mesh.boundary().iter().filter(|face| {
                    face.outward_normal[2] < -0.99
                        && face.vertices.iter().all(|&v| v as usize / 8 == body)
                }) {
                    let t = face.vertices.map(|v| solid.temperature[v as usize]);
                    area += face.area;
                    epsilon += face.area * (0.2 + 0.003 * (mean_power(t, 1) - 250.0));
                    emission += face.area * STEFAN_BOLTZMANN_W_M2_K4
                        * ((0.2 - 0.003 * 250.0) * mean_power(t, 4)
                            + 0.003 * mean_power(t, 5));
                }
                (area, epsilon / area, emission / area)
            });
            let view = if open { 0.25 } else { 1.0 };
            let bath = STEFAN_BOLTZMANN_W_M2_K4 * 420.0_f64.powi(4);
            let reflection = moments.map(|(_, e, _)| (1.0 - e) * view);
            let rhs = moments.map(|(_, e, b)| b + (1.0 - e) * (1.0 - view) * bath);
            let denominator = 1.0 - reflection[0] * reflection[1];
            let j = [(rhs[0] + reflection[0] * rhs[1]) / denominator,
                (rhs[1] + reflection[1] * rhs[0]) / denominator];
            let q: [f64; 2] = std::array::from_fn(|i| moments[i].0
                * (j[i] - view * j[1 - i] - (1.0 - view) * bath));
            for (i, patch) in report.patches.iter().enumerate() {
                assert!((patch.outward_heat_w - q[i]).abs() < 2e-10);
                assert!((patch.radiosity_w_m2 - j[i]).abs() < 1e-8);
                assert!((patch.material_mean_emissivity.unwrap() - moments[i].1).abs() < 1e-13);
                assert!(!patch.material_states.is_empty());
                for sample in &patch.material_states {
                    assert_eq!(sample.emissivity.card_identity(), cards[i].content_hash());
                    assert!((sample.emissivity.value()
                        - (0.2 + 0.003 * (solid.temperature[sample.vertex] - 250.0))).abs() < 1e-13);
                }
            }
            assert!((moments[0].1 - 0.35).abs() > 0.01,
                "reflection retained the preparation emissivity");
            let external = q[0] + q[1];
            if open {
                assert!(external < 0.0);
                assert!((report.ambient.as_ref().unwrap().outward_heat_w + external).abs() < 2e-10);
            } else {
                assert!(q[0] > 0.0 && q[1] < 0.0);
                assert!(external.abs() < 2e-10);
                assert!(report.ambient.is_none());
            }
            assert!((response.thermal.physical.radiation_out_w - external).abs() < 2e-10);
            assert!((stored_change(&state, &old, tank_h) - DT * (20.0 - external)).abs() < 2e-6);
            assert!((result.step.solid_external_heat_j - DT * (20.0 - external)).abs() < 1e-8);
            assert!((DT * response.coolant_out_w - result.step.liquid.returned_heat_j).abs() < 2e-6);
            assert!(result.step.interface_energy_residual_j.abs() < 2e-6);
        });
    }
}
