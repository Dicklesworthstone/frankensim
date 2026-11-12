//! Source profiles and reciprocal storage on the distributed fin extension.
use super::*;
use fs_airflow::graph::thermal::medium::convection::finite_gas::fin::nonlinear::{
    NonlinearSourceFinError, NonlinearSourceFinSample,
};
use fs_convection::fin::nonlinear::{NonlinearFinControls, NonlinearFinError};

fn controls() -> NonlinearFinControls {
    NonlinearFinControls {
        segments: 32,
        maximum_segments: 256,
        maximum_iterations: 256,
        maximum_source_evaluations: 100_000,
        absolute_heat_tolerance_w: 1e-8,
        relative_heat_tolerance: 1e-11,
        absolute_temperature_tolerance_k: 1e-8,
        relaxation: 0.7,
        maximum_transverse_biot: 0.1,
    }
}

fn check_local_balance(
    sample: &NonlinearSourceFinSample<
        fs_airflow::graph::thermal::medium::convection::finite_gas::film::FilmGasConvectionSample<
            fs_material::thermochemical_gas::ResolvedThermochemicalGasState,
            Temperature,
        >,
    >,
    material: &ResolvedConductionLayer<'_>,
    gas: &fs_material::thermochemical_gas::ResolvedThermochemicalGasState,
) {
    let geometry = sample.geometry;
    let nodes = &sample.fin.nodes;
    let dx = geometry.fin().length_m() / (nodes.len() - 1) as f64;
    let perimeter = geometry.fin().exposed_perimeter_m();
    let bulk = sample.bulk_temperature_k;
    let mut side = 0.0;
    for (i, node) in nodes.iter().enumerate() {
        let local = &node.properties.source;
        let film_state = local.property_gas.as_ref().unwrap();
        assert_eq!(&local.interface.gas, gas);
        assert_eq!(film_state.source_packs(), gas.source_packs());
        near(
            film_state.state().temperature,
            0.5 * (node.point.temperature_k + bulk),
            1e-11,
        );
        near(
            local.interface.body.value(),
            node.point.temperature_k,
            1e-12,
        );
        let original_h = plate_reference(film_state, gas.state().density);
        near(node.properties.coefficient_w_m2_k, original_h, 2e-12);
        near(
            node.properties.conductivity_w_m_k,
            material
                .conductivity_curve()
                .conductivity_w_m_k(node.point.temperature_k)
                .unwrap(),
            1e-13,
        );
        let half = if i == 0 || i + 1 == nodes.len() {
            0.5
        } else {
            1.0
        };
        side += original_h * perimeter * dx * half * (node.point.temperature_k - bulk);
    }
    near(side, sample.fin.side_heat_w, 1e-10);
    let tip_heat = sample.fin.tip.as_ref().map_or(0.0, |tip| {
        plate_reference(
            tip.properties.source.property_gas.as_ref().unwrap(),
            gas.state().density,
        ) * geometry.fin().cross_section_m2()
            * (tip.point.temperature_k - bulk)
    });
    near(tip_heat, sample.fin.tip_heat_w, 1e-10);
    // Exact source Kirchhoff difference and the base half control volume,
    // independently reconstructed from the accepted original temperatures.
    let root_edge = geometry.fin().cross_section_m2() / dx
        * material
            .conductivity_curve()
            .integral_w_m(nodes[1].point.temperature_k, nodes[0].point.temperature_k)
            .unwrap();
    let half_base = nodes[0].properties.coefficient_w_m2_k
        * perimeter
        * 0.5
        * dx
        * (sample.base_temperature_k - bulk);
    near(sample.fin.base_heat_w, root_edge + half_base, 2e-10);
    near(
        sample.fin.base_heat_w - side - tip_heat,
        sample.fin.energy_residual_w,
        2e-10,
    );
    let bare = nodes[0].properties.coefficient_w_m2_k
        * geometry.bare_base_area_m2()
        * (sample.base_temperature_k - bulk);
    near(
        sample.applied_heat_into_gas_w,
        bare + f64::from(geometry.count()) * sample.fin.base_heat_w,
        2e-10,
    );
    assert!(Arc::ptr_eq(
        &sample.conductivity_source,
        material.source_line()
    ));
}

#[test]
fn original_local_sources_and_conservative_spatial_heat_survive_both_directions_and_equilibrium() {
    // Affine transport claims, value*(T+100)/500. With claims exactly
    // proportional to T, the prescribed-mass-flux flat-plate h is the same at
    // every film temperature (rho/mu ~ T^-2, so Nu ~ 1/T and h = Nu k/L is
    // constant), and root and tip could differ only by rounding. With the
    // intercept h ~ sqrt((T+100)/T): about 1 % from root to tip here.
    let model = source::model_transport_affine(1.0);
    let gas = model.evaluate(&source::point(350.0, 100000.0)).unwrap();
    let card = solid(true, 5.0);
    let curve = caloric(&card);
    let material = layer(&card, &curve, [250.0, 440.0]);
    for tip in [FinTipExposure::Adiabatic, FinTipExposure::Exposed] {
        let geometry = array(tip);
        let owner = film(&model, geometry)
            .with_prescribed_mass_flow(RATE, FLOW_AREA)
            .unwrap()
            .with_nonlinear_straight_fins(&material, geometry, controls())
            .unwrap();
        assert!(core::ptr::eq(owner.material(), &material));
        for base in [420.0, 290.0, 350.0] {
            let result = owner
                .resolve(&gas, Temperature::new(base), || Ok(()))
                .unwrap();
            check_local_balance(&result.source, &material, &gas);
            assert!(result.source.fin.source_evaluations <= controls().maximum_source_evaluations);
            assert!(result.source.array_absolute_heat_defect_w < 2e-8);
            if base == 350.0 {
                assert_eq!(result.source.applied_heat_into_gas_w, 0.0);
                assert!(
                    result
                        .source
                        .fin
                        .nodes
                        .iter()
                        .all(|n| n.point.temperature_k == base)
                );
                assert!(result.conductance_w_k > 0.0);
            } else {
                assert_eq!(
                    result.source.applied_heat_into_gas_w.is_sign_positive(),
                    base > 350.0
                );
                assert!(result.source.fin.nodes.windows(2).all(|n| {
                    if base > 350.0 {
                        n[0].point.temperature_k >= n[1].point.temperature_k
                    } else {
                        n[0].point.temperature_k <= n[1].point.temperature_k
                    }
                }));
                assert_ne!(
                    result.source.fin.nodes[0].properties.coefficient_w_m2_k,
                    result
                        .source
                        .fin
                        .nodes
                        .last()
                        .unwrap()
                        .properties
                        .coefficient_w_m2_k
                );
            }
        }
    }
}

#[test]
fn nonlinear_source_fin_uses_actual_k_knots_and_preserves_support_work_and_cancel_refusals() {
    let model = source::model_transport(1.0, true);
    let gas = model.evaluate(&source::point(350.0, 100000.0)).unwrap();
    let geometry = StraightFinArray::try_new(
        StraightFinGeometry::try_new(0.08, 0.000002, 0.042).unwrap(),
        10,
        AREA - 10.0 * 0.000002,
        FinTipExposure::Adiabatic,
    )
    .unwrap();
    let selection = film(&model, geometry)
        .with_prescribed_mass_flow(RATE, FLOW_AREA)
        .unwrap();
    let card = conductivity_with_interior_value(4.0);
    let curve = caloric(&card);
    let material = layer(&card, &curve, [250.0, 440.0]);
    let owner = selection
        .with_nonlinear_straight_fins(&material, geometry, controls())
        .unwrap();
    let result = owner
        .resolve(&gas, Temperature::new(400.0), || Ok(()))
        .unwrap();
    check_local_balance(&result.source, &material, &gas);
    let fixed = selection
        .with_straight_fins(
            &material,
            geometry,
            SourceFinLimits {
                maximum_relative_conductivity_variation: 4.0,
                ..limits()
            },
        )
        .unwrap()
        .resolve(&gas, Temperature::new(400.0), || Ok(()))
        .unwrap();
    assert!((result.conductance_w_k - fixed.conductance_w_k).abs() > 0.01 * fixed.conductance_w_k);
    // Both endpoint values can be k=1 while an original interior knot changes
    // the axial Kirchhoff integral; the nonlinear result retains that response.
    assert!(result.source.fin.nodes.last().unwrap().point.temperature_k < 375.0);

    let strict = selection
        .with_nonlinear_straight_fins(
            &material,
            geometry,
            NonlinearFinControls {
                maximum_transverse_biot: 1e-12,
                ..controls()
            },
        )
        .unwrap();
    assert!(
        strict
            .resolve(&gas, Temperature::new(400.0), || Ok(()))
            .is_err()
    );
    let limited = selection.with_nonlinear_straight_fins(
        &material,
        geometry,
        NonlinearFinControls {
            maximum_source_evaluations: 1,
            ..controls()
        },
    );
    // Capacity admission may reject immediately or the first original sweep
    // refuses; neither may fabricate a converged physical sample.
    match limited {
        Ok(owner) => assert!(
            owner
                .resolve(&gas, Temperature::new(400.0), || Ok(()))
                .is_err()
        ),
        Err(_) => {}
    }
    let narrow = layer(&card, &curve, [390.0, 440.0]);
    let unsupported = selection
        .with_nonlinear_straight_fins(&narrow, geometry, controls())
        .unwrap();
    assert!(matches!(
        unsupported.resolve(&gas, Temperature::new(400.0), || Ok(())),
        Err(NonlinearSourceFinError::Solve(NonlinearFinError::Source {
            error: SourceFinError::Thermal(fs_thermal::enthalpy::DepositionError::MaterialDomain),
            ..
        }))
    ));
    let mut polls = 0;
    assert!(matches!(
        owner.resolve(&gas, Temperature::new(400.0), || {
            polls += 1;
            if polls == 4 {
                Err(fs_thermal::enthalpy::DepositionError::Cancelled)
            } else {
                Ok(())
            }
        }),
        Err(NonlinearSourceFinError::Source(SourceFinError::Thermal(
            fs_thermal::enthalpy::DepositionError::Cancelled,
        )))
    ));
    let retry = owner
        .resolve(&gas, Temperature::new(400.0), || Ok(()))
        .unwrap();
    assert_eq!(retry.conductance_w_k, result.conductance_w_k);
}

#[test]
fn nonlinear_fin_profile_drives_the_original_finite_body_and_gas_energy_once() {
    scope(|cx, gate| {
        let old = IsochoricGasParcel::new(
            source::model_transport(1.0, true),
            source::point(350.0, 100000.0),
            0.001,
        )
        .unwrap();
        let saved = old.clone();
        let boundary = uniform_boundary();
        let mut environment = uniform_environment(&old, &boundary);
        environment.exchange.heat_bracket_j = [-20.0, 20.0];
        let wall = WallConvectionSolve {
            temperature_increment_bracket_k: [-2.0, 2.0],
            temperature_tolerance_k: 1e-8,
            maximum_iterations: 128,
        };
        let geometry = array(FinTipExposure::Exposed);
        let card = solid(true, 5.0);
        let curve = caloric(&card);
        let material = layer(&card, &curve, [250.0, 440.0]);
        let owner = film(old.model(), geometry)
            .with_prescribed_mass_flow(RATE, FLOW_AREA)
            .unwrap()
            .with_nonlinear_straight_fins(&material, geometry, controls())
            .unwrap();
        let body = UniformEnthalpyStepInput {
            surface_area_m2: geometry.base_area_m2(),
            duration_s: 0.1,
            ..thermal(&curve, 420.0)
        };
        let next = owner
            .isochoric_body_candidate(cx, &environment, body, wall)
            .unwrap();
        let report = &next.report;
        let endpoint = &report.convection.source.endpoint_source;
        assert_eq!(endpoint.base_temperature_k, next.state.temperature_k());
        check_local_balance(endpoint, &material, report.next_gas.state());
        assert!(
            endpoint.array_absolute_heat_defect_w * body.duration_s
                <= environment.exchange.source_tolerance_j
        );
        near(
            report.solid_to_gas_heat_j,
            body.duration_s * endpoint.applied_heat_into_gas_w,
            2.0 * environment.exchange.source_tolerance_j,
        );
        let body_change = body.mass_kg
            * (next.state.specific_enthalpy_j_kg() - body.initial.specific_enthalpy_j_kg());
        let gas_change = old.mass_kg()
            * (report.next_gas.state().specific_internal_energy_j_kg()
                - old.state().specific_internal_energy_j_kg());
        near(
            body_change + gas_change,
            0.0,
            environment.exchange.energy_tolerance_j,
        );
        near(
            gas_change,
            report.solid_to_gas_heat_j,
            environment.exchange.energy_tolerance_j,
        );
        assert!(report.entropy.convection_production_j_k >= 0.0);
        assert_eq!(old, saved);
        assert!(
            owner
                .resolve_isochoric_coefficient(
                    cx,
                    old.state(),
                    geometry.wetted_area_m2(),
                    Temperature::new(420.0)
                )
                .is_err()
        );
        let radiation = LumpedThermalEnvironment {
            transport: LumpedThermalTransport::try_declared(20.0, 0.2).unwrap(),
            ..boundary.clone()
        };
        assert!(
            owner
                .isochoric_body_candidate(cx, &uniform_environment(&old, &radiation), body, wall)
                .is_err()
        );
        gate.request();
        assert!(matches!(
            owner.isochoric_body_candidate(cx, &environment, body, wall),
            Err(ConvectionError::Thermal(ThermalPistonError::Cancelled))
        ));
        assert_eq!(old, saved);
    });
}
