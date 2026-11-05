//! Native correlation/source films close the original finite gas/solid residual.
use super::*;
use fs_couple::air_path::input::convection::{
    ForcedGasBoundary, ForcedGasGeometry, GasConvectionBoundary, NaturalConvectionGeometry,
    NaturalGasBoundary,
};
use fs_material::gas_species::input::temperature::IsobaricGasLaw;

fn bound<'a>(
    thermal: &DepositionSimulation<'_>,
    card: &'a MaterialCard,
    caloric: &'a ResolvedHeatCapacityCurve,
    context: &'a QueryPoint,
    cx: &Cx<'_>,
) -> SolidEnthalpy<'a> {
    SolidEnthalpy::try_new(
        thermal,
        &PAIR_NODES,
        &PAIR_TETS,
        vec![
            super::super::source(30, card, caloric.curve(), context),
            super::super::source(10, card, caloric.curve(), context),
        ],
        TetAssemblyBudget::standard(),
        0.1,
        config().mesh_frame,
        cx,
    )
    .unwrap()
}

fn sphere(area: f64) -> GasConvectionBoundary {
    GasConvectionBoundary::Natural(NaturalGasBoundary {
        geometry: NaturalConvectionGeometry::Sphere,
        characteristic_length_m: 0.4,
        exposed_area_m2: area,
        gravity_m_s2: 0.0,
    })
}
fn cylinder(area: f64) -> GasConvectionBoundary {
    GasConvectionBoundary::Forced(ForcedGasBoundary {
        geometry: ForcedGasGeometry::CylinderCrossflow,
        characteristic_length_m: 0.03,
        exposed_area_m2: area,
        speed_m_s: 3.0,
    })
}

fn sphere_gas(t: f64, volume: f64) -> IsochoricGasParcel {
    // The parent fixture's M=0.028 gives Pr=0.6929, below the sphere card's
    // strict 0.7 floor. Reuse its native manufactured source at M=0.027,
    // giving Pr=0.7185; retain the same k, source conditions and NASA9 form.
    let (species, transport) = gas_fixture::component("sphere-film-fixture", 0.027, 0.0);
    let model = ThermochemicalGasModel::new(
        FrozenIdealGasMixtureModelV1::new(
            vec![species],
            Composition::new(CompositionBasis::MoleFraction, vec![1.0]).unwrap(),
        )
        .unwrap(),
        vec![transport],
        None,
    )
    .unwrap();
    IsochoricGasParcel::new(model, gas(t, volume).point().clone(), volume).unwrap()
}
fn films(boundary: fn(f64) -> GasConvectionBoundary, distance: f64) -> [SolidGasContact; 2] {
    contacts([0.1, 0.2]).map(|mut contact| {
        contact.geometry.gas_distance_m = 0.0;
        contact.geometry.cell_distance_m = distance;
        contact.convection = Some(SolidGasConvection {
            boundary: boundary(contact.geometry.area_m2),
            solid_profile: SolidGasConductivityProfile::CellPhase,
            maximum_surface_iterations: 96,
            surface_power_tolerance_w: 1e-10,
        });
        contact
    })
}

#[test]
fn natural_sphere_limit_matches_series_solution_and_reverses_finite_heat() {
    with_cx(|cx| {
        // The existing sphere card has a sourced Nu=2 limit at Ra=0. The
        // original synthetic gas has constant k=0.03, so h=2*k/D=0.15 exactly.
        // No constant-h callback or substitute gas is inserted into the solver.
        let k = 0.01;
        let (card, caloric) = transport_material(k);
        let context = point(340.0);
        let dt = 2.0;
        for (temperature, bulk) in [([340.0, 360.0], 450.0), ([440.0, 460.0], 300.0)] {
            for distance in [0.0, 0.02] {
                let thermal = initial_pair(caloric.curve(), temperature, [0.0; 2]);
                let gas = sphere_gas(bulk, 0.05);
                let saved = gas.clone();
                let source = bound(&thermal, &card, &caloric, &context, cx);
                let patches = films(sphere, distance);
                let next = source
                    .step_with_isochoric_gas(
                        &thermal,
                        &gas,
                        dt,
                        0.0,
                        &patches,
                        solve(),
                        &[0.0; 3],
                        cx,
                    )
                    .unwrap();
                let conductance = [0.1, 0.2].map(|area| area / (distance / k + 1.0 / 0.15));
                let expected = reference(
                    [temperature[0], temperature[1], bulk],
                    [10.0, 10.0, 2.5 * 100000.0 * 0.05 / bulk],
                    conductance,
                    dt,
                    [0.0; 3],
                );
                near(
                    next.thermal().temperature_k(CellId(30)).unwrap(),
                    expected[0],
                    2e-7,
                );
                near(
                    next.thermal().temperature_k(CellId(10)).unwrap(),
                    expected[1],
                    2e-7,
                );
                near(next.gas().state().state().temperature, expected[2], 2e-7);
                for sample in &next.report().contacts {
                    let i = if sample.cell == CellId(30) { 0 } else { 1 };
                    let patch = &patches[i];
                    let film = sample.convection.as_ref().unwrap();
                    let q = conductance[i] * (expected[i] - expected[2]);
                    near(sample.solid_to_gas_heat_j, dt * q, 2e-7);
                    near(sample.conductance_w_k, conductance[i], 1e-12);
                    near(
                        film.film.nusselt.as_ref().unwrap().evidence().value,
                        2.0,
                        0.0,
                    );
                    near(
                        film.film.wall_temperature_k,
                        expected[2] + q / (0.15 * patch.geometry.area_m2),
                        2e-7,
                    );
                    assert!(film.surface_balance_residual_w.abs() <= 1e-10);
                    assert!(film.film.entropy_generation_w_k >= 0.0);
                    assert!(film.solid_layer_entropy_generation_w_k >= 0.0);
                    if distance == 0.0 {
                        near(
                            film.film.wall_temperature_k,
                            sample.solid_temperature_k,
                            0.0,
                        );
                        assert_eq!(film.surface_iterations, 0);
                    } else {
                        assert!(
                            (film.film.wall_temperature_k - sample.solid_temperature_k).abs() > 1.0
                        );
                    }
                }
                near(next.report().external_thermal_heat_j, 0.0, 0.0);
                near(next.report().external_gas_heat_j, 0.0, 0.0);
                near(
                    next.report().thermal.stored_energy_change_j
                        + next.report().gas.stored_energy_change_j,
                    0.0,
                    2e-7,
                );
                near(
                    next.report().gas.supplied_heat_j,
                    next.report().solid_to_gas_heat_j,
                    2e-7,
                );
                assert_eq!(gas, saved);
                assert_eq!(thermal.time_s(), 0.0);
            }
        }
    });
}

fn nonlinear_material() -> (MaterialCard, ResolvedHeatCapacityCurve) {
    let (original, _) = transport_material(1.0);
    let mut claims = ClaimSet::new();
    for (_, claim) in original.claims().claims_ordered() {
        let mut claim = claim.clone();
        if claim.key.name() == "thermal_conductivity" {
            claim.value = PropertyValue::Curve {
                abscissa: "T".into(),
                abscissa_dims: Temperature::DIMS,
                knots: vec![(280.0, 0.28), (400.0, 0.4), (600.0, 0.6)],
                dims: THERMAL_CONDUCTIVITY_DIMS,
            };
            claim.interpolation = InterpolationPolicy::LinearInside;
        }
        claims.insert_claim(claim).unwrap();
    }
    let card = MaterialCard::assemble(original.id().clone(), claims, vec![]).unwrap();
    let caloric = caloric_for(&card);
    (card, caloric)
}

#[test]
fn forced_film_uses_actual_surface_and_source_gas_endpoint_pressure() {
    with_cx(|cx| {
        let (card, caloric) = nonlinear_material();
        let context = point(340.0);
        let thermal = initial_pair(caloric.curve(), [340.0, 360.0], [0.0; 2]);
        let gas = gas(450.0, 0.05);
        let source = bound(&thermal, &card, &caloric, &context, cx);
        let patches = films(cylinder, 0.02);
        let next = source
            .step_with_isochoric_gas(&thermal, &gas, 1.0, 0.0, &patches, solve(), &[0.0; 3], cx)
            .unwrap();
        let endpoint_law =
            IsobaricGasLaw::try_from_model(next.gas().model().clone(), next.gas().point().clone())
                .unwrap();
        let initial_law =
            IsobaricGasLaw::try_from_model(gas.model().clone(), gas.point().clone()).unwrap();
        assert!((next.gas().state().state().pressure - gas.state().state().pressure).abs() > 1.0);
        for sample in &next.report().contacts {
            let patch = patches.iter().find(|patch| patch.id == sample.id).unwrap();
            let film = sample.convection.as_ref().unwrap();
            let surface = film.film.wall_temperature_k;
            let tc = sample.solid_temperature_k;
            let tg = next.gas().state().state().temperature;
            let k = sample
                .source
                .property("thermal_conductivity")
                .unwrap()
                .value_si();
            near(k, 0.001 * tc, 1e-14);
            assert_eq!(sample.source.card_identity(), card.content_hash());
            near(sample.source.query().axes()["T"], tc, 0.0);
            assert!(tc < surface && surface < tg);
            let fresh = patch
                .convection
                .unwrap()
                .boundary
                .evaluate(&endpoint_law, surface)
                .unwrap();
            assert_eq!(film.film.property_gas, fresh.property_gas);
            assert_eq!(film.film.nusselt, fresh.nusselt);
            near(film.film.heat_into_wall_w, fresh.heat_into_wall_w, 0.0);
            near(sample.solid_to_gas_heat_j, -fresh.heat_into_wall_w, 1e-10);
            near(
                film.solid_layer_power_w,
                patch.geometry.area_m2 * k / patch.geometry.cell_distance_m * (tc - surface),
                1e-10,
            );
            assert!(film.surface_balance_residual_w.abs() <= 1e-10);
            // The safeguarded source secant should close this smooth nonlinear
            // film without spending the full bisection-scale query budget.
            assert!(film.surface_iterations <= 24);
            near(
                film.film.property_gas.state().pressure,
                next.gas().state().state().pressure,
                1e-9,
            );
            near(
                film.film.property_gas.state().temperature,
                0.5 * (surface + tg),
                1e-12,
            );
            let stale = patch
                .convection
                .unwrap()
                .boundary
                .evaluate(&initial_law, surface)
                .unwrap();
            assert!((stale.conductance_w_k - film.film.conductance_w_k).abs() > 1e-4);
        }
        assert!(next.report().combined_energy_residual_j.abs() <= solve().total_energy_tolerance_j);
        near(
            next.report().thermal.stored_energy_change_j + next.report().gas.stored_energy_change_j,
            0.0,
            2e-7,
        );
        assert_eq!(next.gas().mass_kg(), gas.mass_kg());
        assert_eq!(next.gas().volume_m3(), gas.volume_m3());
        assert_eq!(next.gas().point().axes().get("material-lot"), Some(&7.0));
    });
}

fn integrated_films(distance: f64, span: [f64; 2]) -> [SolidGasContact; 2] {
    films(sphere, distance).map(|mut patch| {
        patch.convection.as_mut().unwrap().solid_profile =
            SolidGasConductivityProfile::IntegratedPurePhase {
                phase: SolidLiquidPhase::Solid,
                temperature_range_k: span,
            };
        patch
    })
}

#[test]
fn integrated_source_profile_closes_exact_affine_heat_in_both_directions() {
    with_cx(|cx| {
        let (card, caloric) = nonlinear_material();
        let context = point(340.0);
        let distance = 0.5;
        let dt = 1.0;
        for (temperature, bulk) in [([340.0, 360.0], 500.0), ([440.0, 460.0], 300.0)] {
            let thermal = initial_pair(caloric.curve(), temperature, [0.0; 2]);
            let gas = sphere_gas(bulk, 0.05);
            let saved = gas.clone();
            let source = bound(&thermal, &card, &caloric, &context, cx);
            // Only the actual material profile must lie in this source span.
            // The finite gas is deliberately outside it in both directions.
            let patches = integrated_films(distance, [320.0, 480.0]);
            let next = source.step_with_isochoric_gas(
                &thermal, &gas, dt, 0.0, &patches, solve(), &[0.0; 3], cx,
            ).unwrap();
            let held = source.step_with_isochoric_gas(
                &thermal, &gas, dt, 0.0, &films(sphere, distance), solve(), &[0.0; 3], cx,
            ).unwrap();
            let tg = next.gas().state().state().temperature;
            assert!(tg < 320.0 || tg > 480.0);
            for sample in &next.report().contacts {
                let patch = patches.iter().find(|p| p.id == sample.id).unwrap();
                let contact = sample.convection.as_ref().unwrap();
                let tc = sample.solid_temperature_k;
                // Independent analytic root for k(T)=0.001*T and the original
                // native sphere limit h=0.15. No material integral/solver helper
                // participates in this expected answer.
                let a = 0.0005 / distance;
                let rhs = a * tc * tc + 0.15 * tg;
                let ts = 2.0 * rhs / (0.15 + (0.15_f64.powi(2) + 4.0 * a * rhs).sqrt());
                let q = 0.15 * patch.geometry.area_m2 * (ts - tg);
                near(contact.film.wall_temperature_k, ts, 2e-8);
                near(contact.solid_layer_power_w,
                    patch.geometry.area_m2 / distance * 0.0005 * (tc - ts) * (tc + ts), 1e-9);
                near(sample.solid_to_gas_heat_j, dt * q, 1e-9);
                near(sample.conductance_w_k, -contact.film.heat_into_wall_w / (tc - tg), 1e-14);
                assert!(contact.surface_balance_residual_w.abs() <= 1e-10);
                assert!(contact.solid_layer_entropy_generation_w_k > 0.0);
                assert!(contact.film.entropy_generation_w_k > 0.0);
                let line = contact.solid_conductivity_profile.as_ref().unwrap();
                assert_eq!(line.axis, "T");
                assert_eq!(line.knots.iter().map(|(t, _)| *t).collect::<Vec<_>>(),
                    vec![320.0, 400.0, 480.0]);
                let receipt = sample.source.property("thermal_conductivity").unwrap().answer();
                assert!(line.knots.iter().all(|(_, answer)| answer.receipt.selected == receipt.receipt.selected));
                assert_eq!(sample.source.card_identity(), card.content_hash());
                assert_eq!(sample.source.query().axes().get("lot"), Some(&7.0));
                let old = held.report().contacts.iter().find(|old| old.id == sample.id).unwrap();
                assert!((old.solid_to_gas_heat_j - sample.solid_to_gas_heat_j).abs() > 1e-4);
                assert!(old.convection.as_ref().unwrap().solid_conductivity_profile.is_none());
            }
            near(next.report().thermal.stored_energy_change_j + next.report().gas.stored_energy_change_j,
                0.0, solve().total_energy_tolerance_j);
            near(next.report().gas.supplied_heat_j, next.report().solid_to_gas_heat_j, 2e-7);
            assert_eq!(thermal.time_s(), 0.0);
            assert_eq!(gas, saved);
        }
    });
}

#[test]
fn integrated_profile_matches_constant_series_and_refuses_unsupported_surfaces_atomically() {
    with_cx(|cx| {
        let (card, caloric) = transport_material(0.01);
        let context = point(340.0);
        let thermal = initial_pair(caloric.curve(), [340.0, 360.0], [0.0; 2]);
        let before = thermal.cells().clone();
        let gas = sphere_gas(500.0, 0.05);
        let saved = gas.clone();
        let source = bound(&thermal, &card, &caloric, &context, cx);
        for span in [[350.0, 480.0], [320.0, 361.0]] {
            // The first excludes a cell; the second admits the cells but the
            // original film requires a surface outside material support.
            assert!(source.step_with_isochoric_gas(
                &thermal, &gas, 1.0, 0.0, &integrated_films(0.02, span),
                solve(), &[0.0; 3], cx,
            ).is_err());
        }
        let mut mixed = integrated_films(0.02, [280.0, 600.0]);
        mixed[0].convection.as_mut().unwrap().solid_profile =
            SolidGasConductivityProfile::IntegratedPurePhase {
                phase: SolidLiquidPhase::SolidLiquid, temperature_range_k: [280.0, 600.0],
            };
        assert!(matches!(source.step_with_isochoric_gas(
            &thermal, &gas, 1.0, 0.0, &mixed, solve(), &[0.0; 3], cx,
        ), Err(SolidGasError::Invalid(_))));
        let next = source.step_with_isochoric_gas(
            &thermal, &gas, 1.0, 0.0, &integrated_films(0.02, [280.0, 600.0]),
            solve(), &[0.0; 3], cx,
        ).unwrap();
        let conductance = [0.1, 0.2].map(|area| area / (0.02 / 0.01 + 1.0 / 0.15));
        let expected = reference([340.0, 360.0, 500.0], [10.0, 10.0, 25.0],
            conductance, 1.0, [0.0; 3]);
        near(next.thermal().temperature_k(CellId(30)).unwrap(), expected[0], 2e-7);
        near(next.thermal().temperature_k(CellId(10)).unwrap(), expected[1], 2e-7);
        near(next.gas().state().state().temperature, expected[2], 2e-7);
        assert_eq!(thermal.cells(), &before);
        assert_eq!(gas, saved);
    });
}

#[test]
fn area_distance_and_unclosed_surface_refuse_without_publishing_heat() {
    with_cx(|cx| {
        let (card, caloric) = nonlinear_material();
        let context = point(340.0);
        let thermal = initial_pair(caloric.curve(), [340.0, 360.0], [0.0; 2]);
        let before = thermal.cells().clone();
        let gas = gas(450.0, 0.05);
        let saved = gas.clone();
        let source = bound(&thermal, &card, &caloric, &context, cx);
        let patches = films(cylinder, 0.02);
        for fault in 0..4 {
            let mut bad = patches.clone();
            match fault {
                0 => bad[0].convection.as_mut().unwrap().boundary = cylinder(0.2),
                1 => bad[0].geometry.gas_distance_m = 0.001,
                2 => bad[0].geometry.cell_distance_m = -0.01,
                _ => {
                    bad[0]
                        .convection
                        .as_mut()
                        .unwrap()
                        .surface_power_tolerance_w = f64::NAN
                }
            }
            assert!(matches!(
                source.step_with_isochoric_gas(
                    &thermal,
                    &gas,
                    1.0,
                    0.0,
                    &bad,
                    solve(),
                    &[0.0; 3],
                    cx,
                ),
                Err(SolidGasError::Invalid(_))
            ));
        }
        let mut exhausted = patches.clone();
        for patch in &mut exhausted {
            let film = patch.convection.as_mut().unwrap();
            film.maximum_surface_iterations = 1;
            film.surface_power_tolerance_w = 0.0;
        }
        match source.step_with_isochoric_gas(
            &thermal,
            &gas,
            1.0,
            0.0,
            &exhausted,
            solve(),
            &[0.0; 3],
            cx,
        ) {
            Err(SolidGasError::SurfaceBalance {
                residual_w,
                iterations,
            }) => {
                assert!(residual_w.is_finite() && residual_w != 0.0);
                assert_eq!(iterations, 1);
            }
            _ => panic!("one nonlinear surface trial must not pass an exact flux gate"),
        }
        let recovered = source
            .step_with_isochoric_gas(&thermal, &gas, 1.0, 0.0, &patches, solve(), &[0.0; 3], cx)
            .unwrap();
        assert!(recovered.report().solid_to_gas_heat_j < 0.0);
        assert_eq!(thermal.cells(), &before);
        assert_eq!(thermal.time_s(), 0.0);
        assert_eq!(gas, saved);
    });
}
