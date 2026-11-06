//! Original material surfaces and independently closed finite-wall heat balances.

use super::*;
use fs_conduction::enthalpy::conduction::implicit::boundary::{
    GrayAmbientBoundary, GrayAmbientControls, GrayTwoBathBoundary, ImplicitBoundaryGeometry,
};
use fs_material::thermal::conduction::{
    LayerConductionConfig, network::ImplicitLayerBoundary, source::ResolvedConductionLayer,
};

fn pure_solid_curve(material: &MaterialCard) -> EquilibriumEnthalpyPhaseCurve {
    EquilibriumEnthalpyPhaseCurve::try_single_phase(
        material.content_hash(),
        SolidLiquidPhase::Solid,
        [(0.0, 300.0), (200.0, 500.0)]
            .into_iter()
            .map(|(h, temperature)| EnthalpyPhaseKnot {
                specific_enthalpy_j_kg: h,
                temperature_k: temperature,
                liquid_mass_fraction: 0.0,
                bulk_density_kg_m3: 2.0,
            })
            .collect(),
    )
    .unwrap()
}

fn layer<'a>(
    material: &'a MaterialCard,
    caloric: &'a EquilibriumEnthalpyPhaseCurve,
    context: &'a QueryPoint,
    range: [f64; 2],
) -> ResolvedConductionLayer<'a> {
    ResolvedConductionLayer::new(
        caloric,
        ThermalMaterialPoint {
            material,
            context,
            temperature_axis: "T",
            liquid_fraction_axis: Some("f"),
        },
        &PropertyKey::new("k", KD),
        ClaimSelection::Policy(SelectionPolicy::SingleClaimOnly),
        SolidLiquidPhase::Solid,
        range,
        || Ok(()),
    )
    .unwrap()
}

#[test]
fn original_layer_integral_and_mixed_faces_close_graded_wall_calorimetry() {
    let fluid_card = card("T", &[(300.0, 2.0), (500.0, 2.0)]);
    let fluid_curve = curve(
        &fluid_card,
        &[(0.0, 300.0, 1.0, 2.0), (200.0, 500.0, 1.0, 2.0)],
    );
    let fluid = source(&fluid_card, &fluid_curve, [0.0, 200.0]);
    let wall_card = card("T", &[(300.0, 1.0), (500.0, 3.0)]);
    let wall_curve = pure_solid_curve(&wall_card);
    let wall = source(&wall_card, &wall_curve, [0.0, 200.0]);
    let context = QueryPoint::new()
        .with_quantity("T", QuantitySpec::dimensional(Temperature::DIMS), 300.0)
        .unwrap()
        .with_quantity("f", QuantitySpec::dimensional(Dims::NONE), 0.0)
        .unwrap();
    let material_layer = layer(&wall_card, &wall_curve, &context, [300.0, 500.0]);
    let geometry = ImplicitBoundaryGeometry {
        area_m2: 1.0,
        half_distance_m: 1.0,
        surface_resistance_m2_k_w: 4.0 / 7.0,
    };
    let config = LayerConductionConfig {
        max_iterations: 96,
        absolute_heat_tolerance_w: 1e-11,
        relative_heat_tolerance: 1e-13,
    };
    let physical = ImplicitLayerBoundary::new(&material_layer, geometry, config).unwrap();
    let face = bath(1, Face3::XMax, 300.0, 4.0 / 7.0);
    let selected = [ThermalVolumeBoundaryFace::Law {
        face,
        law: &physical,
    }];
    let plain = [bath(1, Face3::YMax, 410.0, 0.75)];
    let boundaries = ThermalVolumeBoundaries {
        baths: &plain,
        sources: &selected,
    };
    let mut domain = VoxelDomain::graded([vec![1.0, 2.0], vec![1.0], vec![1.0]]).unwrap();
    domain.set(1, 0, 0, Voxel::Solid(0));
    let field = FlowField::quiescent(&domain);
    let phases = [
        Some(fluid_curve.state_at_specific_enthalpy(268.0 / 3.0).unwrap()),
        None,
    ];
    let walls = [ThermalVolumeSolidCell {
        cell: 1,
        conductivity: &wall,
        mass_kg: 4.0,
        phase: wall_curve
            .state_at_specific_enthalpy(100.0 + 97.0 / 48.0)
            .unwrap(),
        net_enthalpy_increment_j: 0.0,
    }];
    let run = |boundaries: ThermalVolumeBoundaries<'_>, workspace: &mut _, gate: &_| {
        field.advance_isotropic_implicit_conjugate_thermal_volume_with_baths(
            &domain,
            [FlowFace::Wall; 6],
            &fluid,
            boundaries,
            &walls,
            &[],
            input(&[2.0, 0.0], &phases, &[0.0; 2]),
            controls(0.1),
            solve_controls(),
            workspace,
            gate,
        )
    };
    let gate = CancelGate::new_clock_free();
    let mut workspace = ImplicitThermalVolumeWorkspace::new(2).unwrap();
    let result = run(boundaries, &mut workspace, &gate).unwrap();

    // Independent endpoint: Twall=400, Tfluid=390 and k(T)=(T-200)/100.
    // The physical X surface is 350 K: (A/d)*integral_350^400 k(T)dT=87.5 W
    // and (350-300)/(4/7)=87.5 W through its film. Bulk-frozen k(400)
    // would instead remove 280/3 W, so that approximation cannot pass.
    // The plain Y surface has A=2,d=.5,R=.75: G=2, heat inward=20 W.
    // The shared interface has G=1/(.5/2+1/2)=4/3, heat to fluid=40/3 W.
    // dt=.1 therefore gives fluid +4/3 J, wall -97/12 J, baths -27/4 J.
    close(
        result.thermal.volume.transport.specific_enthalpy_j_kg[0],
        90.0,
    );
    close(
        result.solid_endpoints[0].phase.specific_enthalpy_j_kg(),
        100.0,
    );
    close(result.conductive_enthalpy_increment_j[0], 4.0 / 3.0);
    close(
        result.solid_endpoints[0].conductive_enthalpy_increment_j,
        -97.0 / 12.0,
    );
    close(result.report.solid_to_fluid_heat_j, 4.0 / 3.0);
    close(result.report.solid_enthalpy_change_j, -97.0 / 12.0);
    close(result.report.solid_boundary_heat_j, -27.0 / 4.0);
    close(result.report.absolute_solid_boundary_heat_j, 43.0 / 4.0);
    close(
        result.report.conduction.absolute_boundary_heat_j,
        43.0 / 4.0,
    );
    close(result.report.conduction.boundary_heat_j[1], -35.0 / 4.0);
    close(result.report.conduction.boundary_heat_j[3], 2.0);
    close(
        result.report.conduction.solve.energy.external_heat_j,
        -27.0 / 4.0,
    );
    close(
        result.thermal.volume.report.material_enthalpy_change_j
            + result.report.solid_enthalpy_change_j,
        -27.0 / 4.0,
    );
    assert_eq!(result.thermal.volume.transport.mass_kg, [2.0, 0.0]);
    assert!(result.thermal.volume.phase_states[1].is_none());
    assert_eq!(result.solid_endpoints[0].mass_kg, 4.0);
    assert_eq!(result.thermal.volume.transport.report.mass_routes, 0);
    assert_eq!(result.report.conduction.solve.energy.substeps, 1);
    let exchanges = &result.report.conduction.bath_exchanges;
    assert_eq!(exchanges.len(), 2);
    for (face, temperature, heat) in [(Face3::XMax, 300.0, -35.0 / 4.0), (Face3::YMax, 410.0, 2.0)]
    {
        let exchange = exchanges
            .iter()
            .find(|entry| entry.cell == 1 && entry.face == face)
            .unwrap();
        assert_eq!(exchange.channel, 0);
        close(exchange.temperature_k, temperature);
        close(exchange.heat_j, heat);
    }

    // Equal-content caloric objects remain distinct owners. A narrower source
    // admits the bulk endpoint but cannot support the actual 350 K surface.
    let foreign_curve = pure_solid_curve(&wall_card);
    assert_eq!(foreign_curve.identity(), wall_curve.identity());
    let foreign_layer = layer(&wall_card, &foreign_curve, &context, [300.0, 500.0]);
    let limited_layer = layer(&wall_card, &wall_curve, &context, [390.0, 500.0]);
    for unsupported in [&foreign_layer, &limited_layer] {
        let law = ImplicitLayerBoundary::new(unsupported, geometry, config).unwrap();
        let sources = [ThermalVolumeBoundaryFace::Law { face, law: &law }];
        assert!(
            run(
                ThermalVolumeBoundaries {
                    baths: &plain,
                    sources: &sources
                },
                &mut workspace,
                &gate
            )
            .is_err()
        );
    }
    let occupied = [ThermalVolumeBoundaryFace::Law {
        face: bath(1, Face3::XMin, 300.0, 4.0 / 7.0),
        law: &physical,
    }];
    assert!(
        run(
            ThermalVolumeBoundaries {
                baths: &plain,
                sources: &occupied
            },
            &mut workspace,
            &gate
        )
        .is_err()
    );
    assert!(
        run(
            ThermalVolumeBoundaries {
                baths: &[face],
                sources: &selected
            },
            &mut workspace,
            &gate
        )
        .is_err()
    );
    let cancelled = CancelGate::new_clock_free();
    cancelled.request();
    assert!(matches!(
        run(boundaries, &mut workspace, &cancelled),
        Err(ChtError::Cancelled)
    ));
    assert_eq!(result, run(boundaries, &mut workspace, &gate).unwrap());
    assert_eq!(field, FlowField::quiescent(&domain));
    close(walls[0].phase.specific_enthalpy_j_kg(), 100.0 + 97.0 / 48.0);
}

const SIGMA: f64 = 5.670_374_419e-8;

fn gray_emissivity() -> f64 {
    100.0 / (SIGMA * (350.0_f64.powi(4) - 300.0_f64.powi(4)))
}

#[test]
fn native_gray_surface_cools_finite_wall_at_its_physical_surface_temperature() {
    let material = card("T", &[(300.0, 2.0), (500.0, 2.0)]);
    let fluid_curve = curve(
        &material,
        &[(0.0, 300.0, 1.0, 2.0), (200.0, 500.0, 1.0, 2.0)],
    );
    let fluid = source(&material, &fluid_curve, [0.0, 200.0]);
    let wall_curve = pure_solid_curve(&material);
    let wall = source(&material, &wall_curve, [0.0, 200.0]);
    let gray = GrayAmbientBoundary::new(
        &wall,
        gray_emissivity(),
        2.0,
        GrayAmbientControls::default(),
    )
    .unwrap();
    let sources = [ThermalVolumeBoundaryFace::Law {
        face: bath(1, Face3::XMax, 300.0, 0.0),
        law: &gray,
    }];
    let mut domain = VoxelDomain::new(2, 1, 1, 1.0).unwrap();
    domain.set(1, 0, 0, Voxel::Solid(0));
    let phases = [
        Some(fluid_curve.state_at_specific_enthalpy(100.0).unwrap()),
        None,
    ];
    let walls = [ThermalVolumeSolidCell {
        cell: 1,
        conductivity: &wall,
        mass_kg: 2.0,
        phase: wall_curve.state_at_specific_enthalpy(110.0).unwrap(),
        net_enthalpy_increment_j: 0.0,
    }];
    let result = FlowField::quiescent(&domain)
        .advance_isotropic_implicit_conjugate_thermal_volume_with_baths(
            &domain,
            [FlowFace::Wall; 6],
            &fluid,
            ThermalVolumeBoundaries {
                baths: &[],
                sources: &sources,
            },
            &walls,
            &[],
            input(&[2.0, 0.0], &phases, &[0.0; 2]),
            controls(0.1),
            solve_controls(),
            &mut ImplicitThermalVolumeWorkspace::new(2).unwrap(),
            &CancelGate::new_clock_free(),
        )
        .unwrap();
    // At Tcell=400 the HALF-cell G=2/.5=4 conducts 200 W to Ts=350.
    // Convection and gray radiation each remove exactly 100 W to 300 K.
    // M=2,Cp=1,dt=.1 takes Told=410 to400; fluid ends400, hence interface Q=0.
    close(result.solid_endpoints[0].phase.temperature_k(), 400.0);
    close(
        result.thermal.volume.phase_states[0]
            .unwrap()
            .temperature_k(),
        400.0,
    );
    close(result.report.solid_to_fluid_heat_j, 0.0);
    close(result.report.solid_enthalpy_change_j, -20.0);
    close(result.report.solid_boundary_heat_j, -20.0);
    close(result.report.absolute_solid_boundary_heat_j, 20.0);
    close(result.report.conduction.solve.energy.external_heat_j, -20.0);
    assert_eq!(result.thermal.volume.transport.mass_kg, [2.0, 0.0]);
    let exchanges = &result.report.conduction.bath_exchanges;
    assert_eq!(exchanges.len(), 1);
    assert_eq!(
        (exchanges[0].cell, exchanges[0].face, exchanges[0].channel),
        (1, Face3::XMax, 0)
    );
    close(exchanges[0].temperature_k, 300.0);
    close(exchanges[0].heat_j, -20.0);
}

#[test]
fn opposing_native_gray_baths_retain_both_heats_when_wall_receives_zero_net_heat() {
    let material = card("T", &[(300.0, 2.0), (500.0, 2.0)]);
    let fluid_curve = curve(
        &material,
        &[(0.0, 300.0, 1.0, 2.0), (200.0, 500.0, 1.0, 2.0)],
    );
    let fluid = source(&material, &fluid_curve, [0.0, 200.0]);
    let wall_curve = pure_solid_curve(&material);
    let wall = source(&material, &wall_curve, [0.0, 200.0]);
    let gray = GrayTwoBathBoundary::new(
        &wall,
        gray_emissivity(),
        2.0,
        300.0,
        GrayAmbientControls::default(),
    )
    .unwrap();
    let sources = [ThermalVolumeBoundaryFace::Exchange {
        face: bath(1, Face3::XMax, 400.0, 0.0),
        law: &gray,
    }];
    let mut domain = VoxelDomain::new(2, 1, 1, 1.0).unwrap();
    domain.set(1, 0, 0, Voxel::Solid(0));
    let phases = [
        Some(fluid_curve.state_at_specific_enthalpy(50.0).unwrap()),
        None,
    ];
    let walls = [ThermalVolumeSolidCell {
        cell: 1,
        conductivity: &wall,
        mass_kg: 2.0,
        phase: wall_curve.state_at_specific_enthalpy(50.0).unwrap(),
        net_enthalpy_increment_j: 0.0,
    }];
    let result = FlowField::quiescent(&domain)
        .advance_isotropic_implicit_conjugate_thermal_volume_with_baths(
            &domain,
            [FlowFace::Wall; 6],
            &fluid,
            ThermalVolumeBoundaries {
                baths: &[],
                sources: &sources,
            },
            &walls,
            &[],
            input(&[2.0, 0.0], &phases, &[0.0; 2]),
            controls(0.1),
            solve_controls(),
            &mut ImplicitThermalVolumeWorkspace::new(2).unwrap(),
            &CancelGate::new_clock_free(),
        )
        .unwrap();
    // Tcell=Ts=350 K. The400 K gas supplies2*(400-350)=100 W;
    // the300 K radiative bath receives epsilon*sigma*(350^4-300^4)=100 W.
    // Both bodies stay unchanged, but dt=.1 exchanges +10/-10 J through
    // ONE half-cell/face. Gross heat must not collapse to its zero net value.
    close(result.solid_endpoints[0].phase.temperature_k(), 350.0);
    close(
        result.thermal.volume.phase_states[0]
            .unwrap()
            .temperature_k(),
        350.0,
    );
    close(result.report.solid_to_fluid_heat_j, 0.0);
    close(result.report.solid_enthalpy_change_j, 0.0);
    close(result.report.solid_boundary_heat_j, 0.0);
    close(result.report.conduction.boundary_heat_j[1], 0.0);
    close(result.report.absolute_solid_boundary_heat_j, 20.0);
    close(result.report.conduction.absolute_boundary_heat_j, 20.0);
    close(result.report.conduction.solve.energy.external_heat_j, 0.0);
    assert_eq!(result.thermal.volume.transport.mass_kg, [2.0, 0.0]);
    let exchanges = &result.report.conduction.bath_exchanges;
    assert_eq!(exchanges.len(), 2);
    for (channel, temperature, heat) in [(0, 400.0, 10.0), (1, 300.0, -10.0)] {
        let exchange = exchanges
            .iter()
            .find(|entry| entry.channel == channel)
            .unwrap();
        assert_eq!((exchange.cell, exchange.face), (1, Face3::XMax));
        close(exchange.temperature_k, temperature);
        close(exchange.heat_j, heat);
    }
}
