//! Direct original-source volume demand and internal gas-face pressure checks.
use super::*;
use fs_lbm::conjugate::DensityProjectionControls;
use fs_lbm::material::free_surface::{
    MaterialFreeSurface3, MaterialFreeSurfaceVolumeControls3, MaterialFreeSurfaceVolumeInput3,
};
use fs_material::phase::{HeatCapacitySourceQuery, ResolvedHeatCapacityCurve};

#[derive(Clone, Copy)]
struct Case {
    hot_density: f64,
    hot: bool,
    datum: f64,
    fill: f64,
    full_domain: bool,
    numerical_density: f64,
}

impl Default for Case {
    fn default() -> Self {
        Self {
            hot_density: 800.0,
            hot: false,
            datum: 0.0,
            fill: 1.0,
            full_domain: false,
            numerical_density: 1.0,
        }
    }
}

fn controls() -> MaterialFreeSurfaceVolumeControls3 {
    MaterialFreeSurfaceVolumeControls3 {
        projection: DensityProjectionControls {
            duration_s: 0.001,
            relative_tolerance: 1e-12,
            maximum_volume_defect_m3_s: 1e-18,
            maximum_energy_defect_j: 1e-18,
        },
        max_cells: 512,
        maximum_initial_volume_defect_m3: 1e-20,
        maximum_source_density_defect_kg_m3: 1e-10,
        maximum_enthalpy_defect_j: 1e-12,
    }
}

fn with_case(
    case: Case,
    run: impl FnOnce(
        &MaterialFreeSurface3<'_>,
        &LiquidTemperatureLaw<'_>,
        &EquilibriumEnthalpyPhaseCurve,
        &CancelGate,
    ),
) {
    let card = card(0.1, case.hot_density);
    let temperature = if case.hot { 500.0 } else { 300.0 };
    let law = LiquidTemperatureLaw::try_new(
        &card,
        &point(temperature),
        "T",
        LiquidViscosityBasis::Dynamic,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    let source = resolve_liquid_surface_state(
        &card,
        &point(temperature),
        LiquidViscosityBasis::Dynamic,
        LiquidCapillarity::Neglected,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    let chart = ResolvedHeatCapacityCurve::try_new_exact(
        &card,
        HeatCapacitySourceQuery {
            heat_capacity: &PropertyKey::with_quantity(
                "specific_heat_capacity",
                QuantitySpec::dimensional(Dims([2, 0, -2, -1, 0, 0])),
            ),
            density: &PropertyKey::with_quantity(
                "density",
                QuantitySpec::dimensional(Dims([-3, 1, 0, 0, 0, 0])),
            ),
            temperature_axis: "T",
            lower: &point(300.0),
            upper: &point(500.0),
            selection: fs_matdb::SelectionPolicy::SingleClaimOnly,
        },
        SolidLiquidPhase::Liquid,
        case.datum,
    )
    .unwrap();
    let mut cells = if case.full_domain {
        initial(true)
    } else {
        vec![MaterialFreeSurfaceInitialCell3 {
            coordinates: [3; 3],
            fill_fraction: case.fill,
            velocity_m_s: [0.0; 3],
        }]
    };
    for cell in &mut cells {
        cell.velocity_m_s = [0.0; 3];
    }
    let gate = CancelGate::new_clock_free();
    with_cx(&gate, |cx| {
        let scaling = plan_free_surface_from_material(&source, 0.001, 0.001, 0.1).unwrap();
        let flow = if case.numerical_density == 1.0 {
            scaling
                .free_surface3_from_cells(domain(), &cells, cx)
                .unwrap()
        } else {
            let pressure =
                100_000.0 + (case.numerical_density - 1.0) * source.density_kg_m3() / 3.0;
            scaling
                .with_pressure_reference(100_000.0)
                .unwrap()
                .free_surface3_from_cells_with_pressure(domain(), &cells, 0.1, |_| pressure, cx)
                .unwrap()
        };
        run(&flow, &law, chart.curve(), &gate);
    });
}

#[test]
fn g3_actual_source_density_drives_signed_flux_at_fixed_mass_and_enthalpy_datum() {
    for hot_density in [800.0, 1250.0] {
        for hot in [false, true] {
            for datum in [0.0, 1e12] {
                let case = Case {
                    hot_density,
                    hot,
                    datum,
                    ..Case::default()
                };
                with_case(case, |flow, law, chart, gate| {
                    let rho0 = if hot { hot_density } else { 1000.0 };
                    let rho1 = if hot { 1000.0 } else { hot_density };
                    let h0 = datum + if hot { 400_000.0 } else { 0.0 };
                    let mass = flow.region_mass_kg(|x, y, z| [x, y, z] == [3; 3]).unwrap();
                    let requested = mass * if hot { -400_000.0 } else { 400_000.0 };
                    let input = [MaterialFreeSurfaceVolumeInput3 {
                        coordinates: [3; 3],
                        state: chart.state_at_specific_enthalpy(h0).unwrap(),
                        net_enthalpy_increment_j: requested,
                    }];
                    let bits = flow.lattice().grid().state_bits();
                    let ledger = flow.ledger_mass_kg().unwrap().to_bits();
                    let result = flow
                        .project_enthalpy_volume_demand(
                            law,
                            chart,
                            &input,
                            |_, _| 0.0,
                            controls(),
                            gate,
                        )
                        .unwrap();
                    let cell = &result.cells[&[3; 3]];
                    let volume = 1e-9;
                    let delta_volume = mass * (1.0 / rho1 - 1.0 / rho0);
                    let target = delta_volume / controls().projection.duration_s;
                    near(cell.mass_kg, mass, 1e-14);
                    near(cell.represented_enthalpy_increment_j, requested, 1e-14);
                    near(cell.initial_volume_m3, volume, 1e-14);
                    near(cell.proposed_volume_m3, mass / rho1, 1e-14);
                    near(cell.volume_change_m3, delta_volume, 1e-13);
                    near(cell.initial_material.density_kg_m3(), rho0, 1e-14);
                    near(cell.proposed_material.density_kg_m3(), rho1, 1e-14);
                    assert_eq!(result.domain.fluid_count(), 1);
                    assert_eq!(result.gas_pressure_faces.len(), 6);
                    let c = result.domain.index(3, 3, 3);
                    near(result.volume_source_m3_s[c], target, 1e-13);
                    near(
                        result.projection.field.net_outflow(&result.domain, c),
                        target,
                        1e-12,
                    );
                    // Six half-cell outlets: G=2 dt dx/rho0 each, p=S/(6G).
                    let pressure = rho0 * target / (12.0 * 0.001 * 0.001);
                    near(result.projection.pressure_correction_pa[c], pressure, 1e-12);
                    for opening in &result.gas_pressure_faces {
                        near(
                            result.projection.field.outward(
                                &result.domain,
                                opening.cell,
                                opening.face as usize,
                            ),
                            target / 6.0,
                            1e-12,
                        );
                    }
                    assert!(result.projection.report.energy_defect_j.abs() < 1e-18);
                    assert_eq!(flow.lattice().grid().state_bits(), bits);
                    assert_eq!(flow.ledger_mass_kg().unwrap().to_bits(), ledger);
                    assert_eq!(
                        flow.sample_at_si(3, 3, 3).unwrap().unwrap().fill_fraction,
                        1.0
                    );
                });
            }
        }
    }
}

#[test]
fn g0_partial_geometry_and_numerical_pressure_mass_are_not_thermodynamic_expansion() {
    for case in [
        Case {
            fill: 0.5,
            ..Case::default()
        },
        Case {
            numerical_density: 1.04,
            ..Case::default()
        },
    ] {
        with_case(case, |flow, law, chart, gate| {
            let before = flow.lattice().grid().state_bits();
            let input = [MaterialFreeSurfaceVolumeInput3 {
                coordinates: [3; 3],
                state: chart.state_at_specific_enthalpy(0.0).unwrap(),
                net_enthalpy_increment_j: 0.0,
            }];
            assert!(
                flow.project_enthalpy_volume_demand(
                    law,
                    chart,
                    &input,
                    |_, _| 0.0,
                    controls(),
                    gate,
                )
                .is_err()
            );
            assert_eq!(flow.lattice().grid().state_bits(), before);
        });
    }
}

#[test]
fn g0_sealed_expansion_source_mismatch_and_cancellation_publish_nothing() {
    with_case(
        Case {
            full_domain: true,
            ..Case::default()
        },
        |flow, law, chart, gate| {
            let mass = flow.scaling().mass_unit_kg();
            let input: Vec<_> = initial(true)
                .into_iter()
                .map(|cell| MaterialFreeSurfaceVolumeInput3 {
                    coordinates: cell.coordinates,
                    state: chart.state_at_specific_enthalpy(0.0).unwrap(),
                    net_enthalpy_increment_j: mass * 400_000.0,
                })
                .collect();
            let before = flow.lattice().grid().state_bits();
            assert!(
                flow.project_enthalpy_volume_demand(
                    law,
                    chart,
                    &input,
                    |_, _| 0.0,
                    controls(),
                    gate,
                )
                .is_err()
            );
            assert_eq!(flow.lattice().grid().state_bits(), before);
        },
    );
    with_case(Case::default(), |flow, law, chart, gate| {
        let before = flow.lattice().grid().state_bits();
        let mass = flow.scaling().mass_unit_kg();
        let input = [MaterialFreeSurfaceVolumeInput3 {
            coordinates: [3; 3],
            state: chart.state_at_specific_enthalpy(0.0).unwrap(),
            net_enthalpy_increment_j: mass * 400_000.0,
        }];
        assert!(
            flow.project_enthalpy_volume_demand(law, chart, &[], |_, _| 0.0, controls(), gate,)
                .is_err()
        );
        // A second chart with a different hot density cannot replace the source.
        let wrong = caloric(law.card());
        let mismatched = [MaterialFreeSurfaceVolumeInput3 {
            state: wrong.state_at_specific_enthalpy(0.0).unwrap(),
            ..input[0]
        }];
        assert!(
            flow.project_enthalpy_volume_demand(
                law,
                &wrong,
                &mismatched,
                |_, _| 0.0,
                controls(),
                gate,
            )
            .is_err()
        );
        gate.request();
        assert!(
            flow.project_enthalpy_volume_demand(law, chart, &input, |_, _| 0.0, controls(), gate,)
                .is_err()
        );
        assert_eq!(flow.lattice().grid().state_bits(), before);
    });
}

#[test]
fn g3_retained_cp_density_source_admits_interior_steps_with_explicit_density_budget() {
    for hot in [false, true] {
        for magnitude in [0.1, 12345.6789] {
            with_case(
                Case {
                    hot,
                    ..Case::default()
                },
                |flow, law, chart, gate| {
                    assert!(!chart.source_receipts().is_empty());
                    let mass = flow.scaling().mass_unit_kg();
                    let h0 = if hot { 400_000.0 } else { 0.0 };
                    let dh = if hot { -magnitude } else { magnitude };
                    let input = [MaterialFreeSurfaceVolumeInput3 {
                        coordinates: [3; 3],
                        state: chart.state_at_specific_enthalpy(h0).unwrap(),
                        net_enthalpy_increment_j: mass * dh,
                    }];
                    let result = flow
                        .project_enthalpy_volume_demand(
                            law,
                            chart,
                            &input,
                            |_, _| 0.0,
                            controls(),
                            gate,
                        )
                        .unwrap();
                    let cell = &result.cells[&[3; 3]];
                    let represented_dh = cell.proposed.specific_enthalpy_j_kg() - h0;
                    let t0 = if hot { 500.0 } else { 300.0 };
                    near(
                        cell.proposed.temperature_k(),
                        t0 + represented_dh / 2000.0,
                        1e-14,
                    );
                    let rho0 = if hot { 800.0 } else { 1000.0 };
                    let delta_rho = -represented_dh / 2000.0;
                    let rho1 = rho0 + delta_rho;
                    // Independent affine rho(T), Cp=2000 reference. Very small
                    // density differences expose their own floating-point resolution.
                    let expected_volume_change = (mass / rho0) * (-delta_rho / rho1);
                    near(cell.volume_change_m3, expected_volume_change, 1e-8);
                    let c = result.domain.index(3, 3, 3);
                    near(
                        result.projection.field.net_outflow(&result.domain, c),
                        cell.volume_change_m3 / controls().projection.duration_s,
                        1e-12,
                    );
                    assert!(
                        result.maximum_source_density_defect_kg_m3
                            <= controls().maximum_source_density_defect_kg_m3
                    );
                    if result.maximum_source_density_defect_kg_m3 > 0.0 {
                        // T(h) and source rho(T) take independent rounding paths.
                        // Any measured discrepancy requires an explicit positive
                        // allowance; exact agreement may also round identically.
                        let strict = MaterialFreeSurfaceVolumeControls3 {
                            maximum_source_density_defect_kg_m3: 0.0,
                            ..controls()
                        };
                        assert!(matches!(
                            flow.project_enthalpy_volume_demand(
                                law,
                                chart,
                                &input,
                                |_, _| 0.0,
                                strict,
                                gate,
                            ),
                            Err(FreeSurfaceEnthalpyError::Invalid(
                                "caloric/source density discrepancy exceeds its physical accuracy budget"
                            ))
                        ));
                    }
                },
            );
        }
    }
}
