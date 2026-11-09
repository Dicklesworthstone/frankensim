//! Original-source buoyancy and native caloric inventory regressions.

use super::super::boussinesq::BoussinesqConfig;
use super::*;

fn curve() -> ThermalCarreauYasuda {
    ThermalCarreauYasuda::new(
        CarreauYasuda::new(2.0, 0.2, 1.2, 0.45, 2.0, [0.0, 1e4]).unwrap(),
        TemperatureShift::new(
            TemperatureShiftModel::Arrhenius {
                activation_energy_j_mol: 12_000.0,
            },
            300.0,
            [280.0, 360.0],
        )
        .unwrap(),
        ViscosityShiftMode::VerticalAndRate,
    )
}

fn config(gravity_m_s2: [f64; 2]) -> BoussinesqConfig<2> {
    BoussinesqConfig {
        gravity_m_s2,
        max_relative_density_departure: 0.05,
        max_relative_specific_heat_departure: 0.05,
    }
}

fn declared_source(
    temperature: f64,
    density: f64,
    specific_heat: f64,
) -> Result<ThermalFluidTransportState, MaterialError> {
    ThermalFluidTransportState::try_new(
        FluidTransportState::try_new(density, specific_heat, 0.5, 2.0)?,
        temperature,
        None,
    )
}

#[test]
fn zero_density_departure_agrees_with_original_thermal_rheology_owner() {
    with_cx(|cx, _| {
        let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 1, 1));
        let system = FluxSystem::new(&mesh);
        let caloric = Caloric { datum: 1e5 };
        let (old, cells) = inventory(&mesh, &caloric, &[298.0, 308.0], false);
        let material = reference(2.0);
        let owner = LaminarEnthalpy::new(
            cx,
            &system,
            &material,
            MaterialFluxNumerics {
                linear_tolerance: 1e-12,
                ..MaterialFluxNumerics::default()
            },
            &old,
            map(&cells),
            DEPTH,
            2 * mesh.edges.len(),
        )
        .unwrap();
        let velocity = |_| [0.2, 0.0];
        let previous = moments(&mesh, &velocity);
        let boundaries: Vec<_> = (0..mesh.edges.len())
            .filter(|&edge| mesh.is_boundary(edge))
            .map(|edge| MaterialNavierStokesBoundary::Velocity {
                edge,
                value: &velocity,
            })
            .collect();
        let input = LaminarEnthalpyInput {
            history: LaminarEnthalpyHistory {
                velocity_m_s: &previous,
                time_s: 0.0,
            },
            force_density: &zero,
            boundaries: &boundaries,
            slip_heat_routes: &[],
        };
        let rheology = curve();
        let resolve = |point: MaterialFluxRheologyPoint, t: f64, _: &FluidTransportState| {
            rheology
                .evaluate(t, point.shear_rate_per_s)?
                .bind_thermal(source(t, 2.0))
        };
        let plain = owner
            .advance(
                cx,
                &old,
                0.04,
                input,
                controls(),
                &[0.0; 2],
                &resolve,
                |_, _, _| Ok::<_, ()>(exchange(&mesh, caloric.h(320.0))),
            )
            .unwrap();
        let buoyant = owner
            .advance_boussinesq(
                cx,
                &old,
                0.04,
                input,
                BoussinesqConfig {
                    max_relative_density_departure: 0.0,
                    max_relative_specific_heat_departure: 0.0,
                    ..config([0.0, -9.0])
                },
                controls(),
                &[0.0; 2],
                &resolve,
                |_, _, _| Ok::<_, ()>(exchange(&mesh, caloric.h(320.0))),
            )
            .unwrap();
        for (&a, &b) in plain
            .flow()
            .step
            .solution
            .velocity_dofs()
            .iter()
            .zip(buoyant.flow().step.solution.velocity_dofs())
        {
            near(a, b, 2e-9);
        }
        for &cell in &cells {
            near(
                plain.thermal().temperature_k(cell).unwrap(),
                buoyant.thermal().temperature_k(cell).unwrap(),
                2e-8,
            );
        }
        let report = buoyant.report().buoyancy.as_ref().unwrap();
        assert_eq!(report.cell_force_density_n_m3, vec![[0.0; 2]; 2]);
        assert_eq!(report.endpoint_power_w, 0.0);
        assert_eq!(report.endpoint_work_j, 0.0);
        assert_eq!(report.maximum_relative_density_departure, 0.0);
        assert_eq!(report.maximum_relative_specific_heat_departure, 0.0);
        near(
            buoyant.report().thermal.stored_energy_change_j,
            plain.report().thermal.stored_energy_change_j,
            2e-8,
        );
        assert_eq!(
            buoyant
                .thermal()
                .material_caloric_identity(MATERIAL)
                .unwrap(),
            old.material_caloric_identity(MATERIAL).unwrap()
        );
    });
}

#[test]
fn heated_original_density_drives_endpoint_force_and_signed_work_without_extra_heat() {
    with_cx(|cx, _| {
        for gravity in [[0.0, -9.0], [0.0, 9.0]] {
            let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 1, 1));
            let system = FluxSystem::new(&mesh);
            let caloric = Caloric { datum: 1e5 };
            let (old, cells) = inventory(&mesh, &caloric, &[300.0; 2], false);
            let original_cells = old.cells().clone();
            let material = reference(2.0);
            let owner = LaminarEnthalpy::new(
                cx,
                &system,
                &material,
                MaterialFluxNumerics {
                    linear_tolerance: 1e-12,
                    ..MaterialFluxNumerics::default()
                },
                &old,
                map(&cells),
                DEPTH,
                2 * mesh.edges.len(),
            )
            .unwrap();
            let dt = 0.2;
            let endpoint_temperature = 310.0;
            let density = RHO * 0.99;
            let buoyancy = gravity.map(|g| (density - RHO) * g);
            let held_force = [0.12, 0.30];
            let expected_velocity = std::array::from_fn::<_, 2, _>(|axis| {
                dt * (held_force[axis] + buoyancy[axis]) / RHO
            });
            let velocity = |_| expected_velocity;
            let traction = |_: [f64; 2], _: [f64; 2]| [0.0; 2];
            let force = |_| held_force;
            let previous = vec![0.0; system.n_u];
            let boundaries: Vec<_> = (0..mesh.edges.len())
                .filter(|&edge| mesh.is_boundary(edge))
                .map(|edge| MaterialNavierStokesBoundary::OpenReservoir {
                    edge,
                    traction: &traction,
                    velocity: &velocity,
                })
                .collect();
            let input = LaminarEnthalpyInput {
                history: LaminarEnthalpyHistory {
                    velocity_m_s: &previous,
                    time_s: 0.0,
                },
                force_density: &force,
                boundaries: &boundaries,
                slip_heat_routes: &[],
            };
            let rheology = curve();
            let solved = owner
                .advance_boussinesq(
                    cx,
                    &old,
                    dt,
                    input,
                    config(gravity),
                    controls(),
                    &[0.0; 2],
                    |point, t, _| {
                        // Authored small-departure fixture, carried by the ORIGINAL
                        // thermal Carreau binding. Native calorics retain CP below.
                        let original = declared_source(
                            t,
                            RHO * (1.0 - 0.001 * (t - 300.0)),
                            CP * (1.0 + 0.002 * (t - 300.0)),
                        )?;
                        rheology
                            .evaluate(t, point.shear_rate_per_s)?
                            .bind_thermal(original)
                    },
                    |_, _, _| {
                        let mut response = exchange(&mesh, caloric.h(endpoint_temperature));
                        for &cell in &cells {
                            response
                                .heat
                                .bulk
                                .sources_w
                                .insert(cell, old.cells()[&cell].mass_kg * CP * 10.0 / dt);
                        }
                        Ok::<_, ()>(response)
                    },
                )
                .unwrap();
            let (error, divergence) =
                system.velocity_error(&solved.flow().step.solution.discrete().x, &velocity);
            assert!(error < 2e-8 && divergence < 2e-8, "{error}, {divergence}");
            let report = solved.report().buoyancy.as_ref().unwrap();
            let power =
                DEPTH * (buoyancy[0] * expected_velocity[0] + buoyancy[1] * expected_velocity[1]);
            near(report.endpoint_power_w, power, 2e-8);
            near(report.endpoint_work_j, dt * power, 2e-8);
            assert!(report.endpoint_work_j.abs() > 1e-5);
            assert_eq!(
                report.endpoint_work_j.is_sign_positive(),
                power.is_sign_positive()
            );
            near(report.maximum_relative_density_departure, 0.01, 1e-9);
            near(report.maximum_relative_specific_heat_departure, 0.02, 1e-9);
            for (cell, sample) in solved.flow().samples().iter().enumerate() {
                near(
                    solved.thermal().temperature_k(cells[cell]).unwrap(),
                    310.0,
                    2e-8,
                );
                near(report.cell_force_density_n_m3[cell][0], buoyancy[0], 1e-9);
                near(report.cell_force_density_n_m3[cell][1], buoyancy[1], 1e-9);
                let projection = sample.material().source();
                let original = projection.material().source();
                near(
                    original.rheology().temperature().temperature_k(),
                    310.0,
                    2e-8,
                );
                near(original.source().temperature_k(), 310.0, 2e-8);
                near(
                    projection
                        .material()
                        .thermal_state()
                        .transport()
                        .density_kg_m3(),
                    density,
                    1e-9,
                );
                near(
                    projection.resolved_transport().specific_heat_j_kg_k(),
                    CP * 1.02,
                    1e-9,
                );
                assert_eq!(sample.transport().density_kg_m3(), RHO);
                assert_eq!(sample.transport().specific_heat_j_kg_k(), CP);
                near(
                    sample.transport().dynamic_viscosity_pa_s(),
                    original.transport_state().unwrap().dynamic_viscosity_pa_s(),
                    1e-12,
                );
                assert_eq!(
                    sample
                        .material()
                        .thermal_state()
                        .volumetric_expansion_per_k(),
                    None
                );
            }
            // Uniform endpoint enthalpy has zero net stream heat. Uniform
            // translation has zero strain heat. Buoyancy work is NOT added to
            // the separately authored 24 J of native caloric heating.
            let expected_heat = RHO * DEPTH * CP * 10.0;
            near(
                solved.report().thermal.stored_energy_change_j,
                expected_heat,
                2e-8,
            );
            near(solved.report().other_external_heat_j, expected_heat, 2e-8);
            near(solved.report().bulk_viscous_heat_j, 0.0, 1e-12);
            near(solved.report().slip_heat_j, 0.0, 1e-12);
            near(solved.report().external_stream_enthalpy_j, 0.0, 2e-8);
            near(solved.report().thermal_balance_residual_j, 0.0, 1e-8);
            assert_eq!(
                solved
                    .thermal()
                    .material_caloric_identity(MATERIAL)
                    .unwrap(),
                old.material_caloric_identity(MATERIAL).unwrap()
            );
            assert_eq!(old.cells(), &original_cells);
            assert_eq!(old.time_s(), 0.0);
            assert!(previous.iter().all(|&value| value == 0.0));
        }
    });
}

#[test]
fn departure_source_and_late_cancellation_refusals_preserve_both_original_histories() {
    with_cx(|cx, gate| {
        let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 1, 1));
        let system = FluxSystem::new(&mesh);
        let caloric = Caloric { datum: 0.0 };
        let (old, cells) = inventory(&mesh, &caloric, &[300.0; 2], false);
        let material = reference(2.0);
        let owner = LaminarEnthalpy::new(
            cx,
            &system,
            &material,
            MaterialFluxNumerics {
                linear_tolerance: 1e-12,
                ..MaterialFluxNumerics::default()
            },
            &old,
            map(&cells),
            DEPTH,
            2 * mesh.edges.len(),
        )
        .unwrap();
        let velocity = |_| [0.1, 0.0];
        let previous = moments(&mesh, &velocity);
        let original_velocity = previous.clone();
        let original_cells = old.cells().clone();
        let boundaries: Vec<_> = (0..mesh.edges.len())
            .filter(|&edge| mesh.is_boundary(edge))
            .map(|edge| MaterialNavierStokesBoundary::Velocity {
                edge,
                value: &velocity,
            })
            .collect();
        let input = LaminarEnthalpyInput {
            history: LaminarEnthalpyHistory {
                velocity_m_s: &previous,
                time_s: 0.0,
            },
            force_density: &zero,
            boundaries: &boundaries,
            slip_heat_routes: &[],
        };
        let rheology = curve();
        for mode in 0..3 {
            let refused = owner.advance_boussinesq(
                cx,
                &old,
                0.04,
                input,
                config([0.0, -9.0]),
                controls(),
                &[0.0; 2],
                |point, t, _| {
                    if mode == 2 {
                        return Err(material_error("original source outside support"));
                    }
                    let original = declared_source(
                        t,
                        if mode == 0 { 1.1 * RHO } else { RHO },
                        if mode == 1 { 1.1 * CP } else { CP },
                    )?;
                    rheology
                        .evaluate(t, point.shear_rate_per_s)?
                        .bind_thermal(original)
                },
                |_, _, _| Ok::<_, ()>(exchange(&mesh, caloric.h(300.0))),
            );
            assert!(matches!(refused, Err(LaminarEnthalpyError::Momentum(_))));
            assert_eq!(previous, original_velocity);
            assert_eq!(old.cells(), &original_cells);
            assert_eq!(old.time_s(), 0.0);
        }
        let resolve = |point: MaterialFluxRheologyPoint, t: f64, _: &FluidTransportState| {
            rheology
                .evaluate(t, point.shear_rate_per_s)?
                .bind_thermal(source(t, 2.0))
        };
        let retry = owner
            .advance_boussinesq(
                cx,
                &old,
                0.04,
                input,
                config([0.0, -9.0]),
                controls(),
                &[0.0; 2],
                &resolve,
                |_, _, _| Ok::<_, ()>(exchange(&mesh, caloric.h(300.0))),
            )
            .unwrap();
        near(retry.report().thermal_balance_residual_j, 0.0, 1e-8);
        let reached_flow = std::cell::Cell::new(false);
        let cancelled = owner.advance_boussinesq(
            cx,
            &old,
            0.04,
            input,
            config([0.0, -9.0]),
            controls(),
            &[0.0; 2],
            &resolve,
            |_, _, flow| {
                assert_eq!(flow.samples().len(), cells.len());
                reached_flow.set(true);
                gate.request();
                Ok::<_, ()>(exchange(&mesh, caloric.h(300.0)))
            },
        );
        assert!(reached_flow.get());
        assert!(matches!(cancelled, Err(LaminarEnthalpyError::Cancelled)));
        assert_eq!(previous, original_velocity);
        assert_eq!(old.cells(), &original_cells);
        assert_eq!(old.time_s(), 0.0);
    });
}
