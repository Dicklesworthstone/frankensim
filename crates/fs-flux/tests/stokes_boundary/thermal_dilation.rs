//! G1 thermal expansion, source compatibility and exact native divergence.

use super::*;
use fs_blake3::ContentHash;
use fs_flux::ns::boundary::thermal_dilation::{CaloricDilationCell, CaloricDilationError};
use fs_material::phase::{EnthalpyPhaseKnot, EquilibriumEnthalpyPhaseCurve, SolidLiquidPhase};

fn caloric(final_density: f64, phase: SolidLiquidPhase) -> EquilibriumEnthalpyPhaseCurve {
    let fraction = if phase == SolidLiquidPhase::Liquid {
        1.0
    } else {
        0.0
    };
    EquilibriumEnthalpyPhaseCurve::try_single_phase(
        ContentHash([0x94; 32]),
        phase,
        vec![
            EnthalpyPhaseKnot {
                temperature_k: 300.0,
                specific_enthalpy_j_kg: 0.0,
                liquid_mass_fraction: fraction,
                bulk_density_kg_m3: 1000.0,
            },
            EnthalpyPhaseKnot {
                temperature_k: 400.0,
                specific_enthalpy_j_kg: 100.0,
                liquid_mass_fraction: fraction,
                bulk_density_kg_m3: final_density,
            },
        ],
    )
    .unwrap()
}

fn cell(curve: &EquilibriumEnthalpyPhaseCurve, rate: f64) -> CaloricDilationCell<'_> {
    CaloricDilationCell {
        caloric: curve,
        specific_enthalpy_j_kg: 0.0,
        material_specific_enthalpy_rate_w_kg: rate,
    }
}

#[test]
fn heated_open_liquid_expands_and_negative_expansivity_or_cooling_reverses_flow() {
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 4, 3));
    let system = FluxSystem::new(&mesh);
    let zero = |_| [0.0; 2];
    let traction = |_, _| [0.0; 2];
    let boundaries: Vec<_> = mesh
        .edges
        .iter()
        .enumerate()
        .filter(|(e, _)| mesh.is_boundary(*e))
        .map(|(edge, data)| {
            if data.mid[0] == 1.0 {
                StokesBoundary::TotalTraction {
                    edge,
                    value: &traction,
                }
            } else if data.mid[0] == 0.0 {
                StokesBoundary::Velocity { edge, value: &zero }
            } else {
                StokesBoundary::NavierSlip {
                    edge,
                    value: &zero,
                    tangential_resistance_m_s: 0.0,
                }
            }
        })
        .collect();
    for (density, enthalpy, rate, expected) in [
        (800.0, 0.0, 40.0, 0.1),
        (800.0, 50.0, -45.0, -0.1),
        (4000.0 / 3.0, 0.0, 40.0, -0.1),
    ] {
        let curve = caloric(density, SolidLiquidPhase::Liquid);
        let source = vec![
            CaloricDilationCell {
                specific_enthalpy_j_kg: enthalpy,
                ..cell(&curve, rate)
            };
            mesh.tris.len()
        ];
        let solution = system
            .assemble_stokes_caloric_dilation(
                FluxParams::default(),
                &|_, _| 2.0,
                &zero,
                &boundaries,
                &source,
            )
            .unwrap()
            .solve(2e-11)
            .unwrap();
        for point in [[0.1, 0.2], [0.5, 0.5], [0.9, 0.8]] {
            let velocity = system.velocity_at(&solution.flow.field.x, point);
            assert!(
                (velocity[0] - expected * point[0]).abs() < 2e-11,
                "{velocity:?}"
            );
            assert!(velocity[1].abs() < 2e-11);
        }
        for pressure in solution.flow.pressure {
            assert!((pressure - 4.0 * expected).abs() < 2e-10, "{pressure}");
        }
        assert!(solution.maximum_dilation_defect_per_s < 2e-11);
        for response in solution.source {
            assert!((response.dilation_per_s - expected).abs() < 1e-15);
            assert_eq!(
                response.state.material_card_identity(),
                ContentHash([0x94; 32])
            );
        }
    }
}

#[test]
fn opposed_nonuniform_heating_and_cooling_flow_inside_a_sealed_domain() {
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 4, 2));
    let system = FluxSystem::new(&mesh);
    let zero = |_| [0.0; 2];
    let boundaries: Vec<_> = mesh
        .edges
        .iter()
        .enumerate()
        .filter(|(e, _)| mesh.is_boundary(*e))
        .map(|(edge, _)| StokesBoundary::NavierSlip {
            edge,
            value: &zero,
            tangential_resistance_m_s: 0.0,
        })
        .collect();
    let curve = caloric(800.0, SolidLiquidPhase::Liquid);
    let source: Vec<_> = mesh
        .tris
        .iter()
        .map(|tri| {
            let x = tri.iter().map(|&n| mesh.verts[n][0]).sum::<f64>() / 3.0;
            if x < 0.5 {
                cell(&curve, 40.0)
            } else {
                CaloricDilationCell {
                    specific_enthalpy_j_kg: 100.0,
                    ..cell(&curve, -50.0)
                }
            }
        })
        .collect();
    let solution = system
        .assemble_stokes_caloric_dilation(
            FluxParams::default(),
            &|_, _| 1.0,
            &zero,
            &boundaries,
            &source,
        )
        .unwrap()
        .solve(2e-11)
        .unwrap();
    for point in [[0.2, 0.2], [0.4, 0.8], [0.7, 0.4], [0.9, 0.3]] {
        let u = system.velocity_at(&solution.flow.field.x, point);
        assert!((u[0] - 0.1 * point[0].min(1.0 - point[0])).abs() < 2e-11);
        assert!(u[1].abs() < 2e-11);
    }
    assert_eq!(
        solution.flow.pressure_reference,
        StokesPressureReference::MeanZero
    );
}

#[test]
fn rigid_cavities_refuse_net_expansion_per_component_even_when_global_sources_cancel() {
    let mesh = separate_domains();
    let system = FluxSystem::new(&mesh);
    let zero = |_| [0.0; 2];
    let boundaries = velocity_boundaries(&mesh, &zero);
    let curve = caloric(800.0, SolidLiquidPhase::Liquid);
    for reverse_second in [false, true] {
        let source: Vec<_> = system
            .pressure_components()
            .iter()
            .map(|&component| {
                if reverse_second && component == 1 {
                    CaloricDilationCell {
                        specific_enthalpy_j_kg: 100.0,
                        ..cell(&curve, -50.0)
                    }
                } else {
                    cell(&curve, 40.0)
                }
            })
            .collect();
        let result = system.assemble_stokes_caloric_dilation(
            FluxParams::default(),
            &|_, _| 1.0,
            &zero,
            &boundaries,
            &source,
        );
        assert!(matches!(
            result,
            Err(CaloricDilationError::Flow(
                StokesBoundaryError::IncompatibleFlux { component: 0, .. }
            ))
        ));
    }
}

#[test]
fn compatible_moving_boundary_and_zero_source_reuse_the_original_stokes_operator() {
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 3, 2));
    let system = FluxSystem::new(&mesh);
    let curve = caloric(800.0, SolidLiquidPhase::Liquid);
    let force = |_| [0.3, -0.2];
    for rate in [0.0, 40.0] {
        let expected = 0.0025 * rate;
        let velocity = |p: [f64; 2]| [expected * p[0], 0.0];
        let boundaries = velocity_boundaries(&mesh, &velocity);
        let source = vec![cell(&curve, rate); mesh.tris.len()];
        let thermal = system
            .assemble_stokes_caloric_dilation(
                FluxParams::default(),
                &|_, _| 1.0,
                &force,
                &boundaries,
                &source,
            )
            .unwrap();
        let solution = thermal.solve(2e-11).unwrap();
        let u = system.velocity_at(&solution.flow.field.x, [0.6, 0.7]);
        assert!((u[0] - expected * 0.6).abs() < 2e-11);
        assert!(u[1].abs() < 2e-11);
        if rate == 0.0 {
            let ordinary = system
                .assemble_stokes_mixed_variable_viscosity(
                    FluxParams::default(),
                    &|_, _| 1.0,
                    &force,
                    &boundaries,
                )
                .unwrap();
            assert_eq!(thermal.native().rhs(), ordinary.rhs());
            assert_eq!(solution.flow.field.x, ordinary.solve().unwrap().field.x);
        }
    }
}

#[test]
fn unsupported_phases_and_nonfinite_or_out_of_domain_caloric_rates_refuse() {
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 1, 1));
    let system = FluxSystem::new(&mesh);
    let zero = |_| [0.0; 2];
    let boundaries = velocity_boundaries(&mesh, &zero);
    let solid = caloric(800.0, SolidLiquidPhase::Solid);
    let liquid = caloric(800.0, SolidLiquidPhase::Liquid);
    let mut knots = liquid.knots().to_vec();
    knots.insert(1, EnthalpyPhaseKnot {
        temperature_k: 350.0,
        specific_enthalpy_j_kg: 50.0,
        liquid_mass_fraction: 1.0,
        bulk_density_kg_m3: 950.0,
    });
    let kink = EquilibriumEnthalpyPhaseCurve::try_single_phase(
        ContentHash([0x94; 32]), SolidLiquidPhase::Liquid, knots,
    ).unwrap();
    for input in [
        cell(&solid, 0.0),
        cell(&liquid, f64::NAN),
        cell(&liquid, -1.0),
        CaloricDilationCell {
            specific_enthalpy_j_kg: 100.0,
            ..cell(&liquid, 1.0)
        },
        CaloricDilationCell {
            specific_enthalpy_j_kg: 50.0,
            ..cell(&kink, -1.0)
        },
        CaloricDilationCell {
            specific_enthalpy_j_kg: 101.0,
            ..cell(&liquid, 0.0)
        },
    ] {
        assert!(
            system
                .assemble_stokes_caloric_dilation(
                    FluxParams::default(),
                    &|_, _| 1.0,
                    &zero,
                    &boundaries,
                    &[input; 2],
                )
                .is_err()
        );
    }
}
