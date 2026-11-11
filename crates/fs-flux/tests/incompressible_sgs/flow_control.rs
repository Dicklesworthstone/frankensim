//! Prescribed integral flows, physical pressure and native endpoint turbulence.
use super::*;
use fs_flux::material::{
    MaterialFluxError, MaterialFluxNumerics, MaterialFluxRheologyPoint, MaterialFluxSystem,
    MaterialNavierStokesBoundary as Boundary, MaterialSgsConfig, MaterialSgsError,
};
use fs_material::MaterialError;
use fs_material::fluid_transport::{FluidTransport, FluidTransportState};
use std::cell::Cell;
use std::rc::Rc;

fn controls() -> NavierStokesControls {
    NavierStokesControls {
        maximum_iterations: 120,
        relative_tolerance: 1e-10,
        relaxation: 0.8,
    }
}
fn numerics() -> MaterialFluxNumerics {
    MaterialFluxNumerics {
        linear_tolerance: 1e-12,
        ..MaterialFluxNumerics::default()
    }
}
fn molecular(rate: f64) -> f64 {
    0.1 + 0.5 / (1.0 + 0.6 * rate)
}
fn transport(rho: f64, nu: f64) -> FluidTransportState {
    FluidTransportState::try_new(rho, 1000.0, 0.3, rho * nu).unwrap()
}
struct Sample {
    point: MaterialFluxRheologyPoint,
    state: FluidTransportState,
    reads: Rc<Cell<usize>>,
}
impl FluidTransport for Sample {
    fn transport_state(&self) -> Result<FluidTransportState, MaterialError> {
        self.reads.set(self.reads.get() + 1);
        Ok(self.state)
    }
}
fn inlet(mesh: &TriMesh) -> Vec<usize> {
    mesh.edges
        .iter()
        .enumerate()
        .filter(|(edge, face)| mesh.is_boundary(*edge) && face.normal[0] < -0.5)
        .map(|(edge, _)| edge)
        .collect()
}

#[test]
fn controlled_sgs_pressure_closes_steady_and_bdf_endpoints_with_original_material_and_reactions() {
    let length = 2.0;
    let height = 0.7;
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(length, height, 3, 2));
    let system = FluxSystem::new(&mesh);
    let inlet = inlet(&mesh);
    let base = |p: [f64; 2]| [0.4 + 0.2 * p[0] + 0.3 * p[1], 0.1 - 0.2 * p[1]];
    let width = 0.4;
    let widths = vec![width; mesh.tris.len()];
    let time = 0.3;
    let dt = 0.1;
    let rho = 7.0;
    let reference = transport(rho, molecular(0.0));
    let owner = MaterialFluxSystem::new(&system, &reference, numerics()).unwrap();
    for model in [
        SgsModel::Smagorinsky { cs: 0.5 },
        SgsModel::Wale { cw: 0.5 },
        SgsModel::Vreman { cv: 0.07 },
    ] {
        // True steady, BE, equal-step BDF2, and unequal-interval BDF2.
        for scheme in 0..4 {
            let previous_dt = if scheme == 3 { 0.15 } else { dt };
            let amplitude = |t: f64| {
                if scheme == 0 {
                    1.0
                } else {
                    1.0 + 0.2 * t + if scheme >= 2 { 0.3 * t * t } else { 0.0 }
                }
            };
            let scale = amplitude(time);
            let derivative = if scheme == 0 {
                0.0
            } else {
                0.2 + if scheme >= 2 { 0.6 * time } else { 0.0 }
            };
            let old_scale = if scheme == 0 {
                0.6
            } else {
                amplitude(time - dt)
            };
            let previous = moments(&mesh, &|p| base(p).map(|u| old_scale * u));
            let older = moments(&mesh, &|p| {
                base(p).map(|u| amplitude(time - dt - previous_dt) * u)
            });
            let saved = [previous.clone(), older.clone()];
            let exact = |p| base(p).map(|u| scale * u);
            let rate = 0.5 * scale;
            let gradient = [[0.2 * scale, 0.3 * scale], [0.0, -0.2 * scale]];
            let nu_t = trace_free_planar_reference(model, width, gradient);
            assert!(nu_t > 0.0);
            let total_mu = rho * (molecular(rate) + nu_t);
            // The missing uniform acceleration in force must be supplied by
            // the UNKNOWN inlet pressure. The length factor is physical.
            let missing_acceleration = 0.25 + derivative;
            let pressure_drop = rho * missing_acceleration * length;
            let force = |p: [f64; 2]| {
                let u = base(p);
                [
                    rho * (derivative * u[0] + scale * scale * (0.04 * p[0] + 0.11)
                        - missing_acceleration),
                    rho * (derivative * u[1] + scale * scale * (0.04 * p[1] - 0.02)) - 3.0,
                ]
            };
            for policy in [
                OpenBoundaryPolicy::NaturalTraction,
                OpenBoundaryPolicy::DissipativeBackflow,
            ] {
                let traction = |p: [f64; 2], n: [f64; 2]| {
                    let pressure = 7.0 - 3.0 * p[1];
                    let u = exact(p);
                    let backflow = if policy == OpenBoundaryPolicy::DissipativeBackflow {
                        rho * (u[0] * n[0] + u[1] * n[1]).min(0.0)
                    } else {
                        0.0
                    };
                    [
                        (0.4 * scale * total_mu - pressure) * n[0] + 0.3 * scale * total_mu * n[1]
                            - backflow * u[0],
                        0.3 * scale * total_mu * n[0] + (-0.4 * scale * total_mu - pressure) * n[1]
                            - backflow * u[1],
                    ]
                };
                let boundaries: Vec<_> = mesh
                    .edges
                    .iter()
                    .enumerate()
                    .filter(|(edge, _)| mesh.is_boundary(*edge))
                    .map(|(edge, face)| {
                        if face.normal[0].abs() > 0.5 {
                            Boundary::Open {
                                edge,
                                traction: &traction,
                                backflow: policy,
                            }
                        } else {
                            Boundary::Velocity {
                                edge,
                                value: &exact,
                            }
                        }
                    })
                    .collect();
                let target_flow = -scale * height * (0.4 + 0.15 * height);
                let reads = Rc::new(Cell::new(0));
                let mut first_rates = vec![None; mesh.tris.len()];
                let mut queries = 0;
                let resolve = |point: MaterialFluxRheologyPoint, retained: &FluidTransportState| {
                    assert_eq!(retained, &reference);
                    first_rates[point.cell].get_or_insert(point.shear_rate_per_s);
                    queries += 1;
                    Ok(Sample {
                        point,
                        state: transport(rho, molecular(point.shear_rate_per_s)),
                        reads: Rc::clone(&reads),
                    })
                };
                let config = MaterialSgsConfig {
                    model,
                    filter_widths_m: &widths,
                };
                let solved = match scheme {
                    0 => owner.solve_navier_stokes_mixed_rheology_sgs_flow_control(
                        config,
                        &force,
                        &boundaries,
                        &[&inlet],
                        &[target_flow],
                        &previous,
                        controls(),
                        1e-12,
                        resolve,
                    ),
                    1 => owner.bdf1_navier_stokes_mixed_rheology_sgs_flow_control(
                        config,
                        &force,
                        &boundaries,
                        &[&inlet],
                        &[target_flow],
                        &previous,
                        dt,
                        controls(),
                        1e-12,
                        resolve,
                    ),
                    _ => owner.bdf2_navier_stokes_mixed_rheology_sgs_variable_step_flow_control(
                        config,
                        &force,
                        &boundaries,
                        &[&inlet],
                        &[target_flow],
                        &previous,
                        Some((&older, previous_dt)),
                        dt,
                        controls(),
                        1e-12,
                        resolve,
                    ),
                }
                .unwrap();
                let flow = &solved.flow.flow;
                let solution = &flow.step.solution;
                let (error, divergence) = system.velocity_error(&solution.discrete().x, &exact);
                assert!(
                    error < 2e-8 && divergence < 2e-8,
                    "{model:?}, {scheme}: {error:e}, {divergence:e}"
                );
                close(solved.pressure_offsets_pa[0], pressure_drop, 3e-7);
                close(solved.outward_volume_rates_m2_s[0], target_flow, 2e-11);
                assert!(solved.maximum_flow_residual_m2_s < 2e-11);
                assert!(flow.step.nonlinear_residual <= controls().relative_tolerance);
                assert!(flow.step.relative_update <= controls().relative_tolerance);
                assert!(
                    flow.step
                        .open_boundary_policies
                        .iter()
                        .all(|(_, actual)| *actual == policy)
                );
                for (&pressure, p) in solution.pressure_pa().iter().zip(&mesh.centroids) {
                    close(
                        pressure,
                        7.0 - 3.0 * p[1] + pressure_drop * (1.0 - p[0] / length),
                        3e-7,
                    );
                }
                let rates = system.shear_rates(&solution.discrete().x);
                for (cell, sample) in flow.samples().iter().enumerate() {
                    assert_eq!(sample.point(), sample.material().point);
                    assert_eq!(
                        sample.point().shear_rate_per_s.to_bits(),
                        rates[cell].to_bits()
                    );
                    close(first_rates[cell].unwrap(), 0.5 * old_scale, 2e-12);
                    close(
                        sample.transport().dynamic_viscosity_pa_s(),
                        rho * molecular(rate),
                        2e-7,
                    );
                    close(solved.flow.sgs.cells[cell].eddy_viscosity_m2_s, nu_t, 2e-8);
                }
                close(
                    solution.viscous_dissipation().total_power_w_per_m(),
                    length * height * rho * molecular(rate) * rate * rate,
                    3e-7,
                );
                close(
                    solved.flow.modeled_transfer_w_per_m,
                    length * height * rho * nu_t * rate * rate,
                    3e-7,
                );
                // Reactions must include the solved port pressure and TOTAL
                // molecular+SGS stress, while heating remains molecular-only.
                let load = solution.boundary_loads(&inlet, [0.0; 2], &[None]).unwrap();
                close(
                    load.stress.force_n_per_m[0],
                    height * (0.4 * scale * total_mu - 7.0 - pressure_drop) + 1.5 * height * height,
                    5e-7,
                );
                close(
                    load.stress.force_n_per_m[1],
                    0.3 * scale * total_mu * height,
                    5e-7,
                );
                close(
                    load.stress.moment_n,
                    0.5 * height * height * (7.0 + pressure_drop - 0.4 * scale * total_mu)
                        - height.powi(3),
                    5e-7,
                );
                assert_eq!(
                    reads.get(),
                    queries,
                    "reporting must retain endpoint sources"
                );
                assert_eq!(previous, saved[0]);
                assert_eq!(older, saved[1]);
            }
        }
    }
}

#[test]
fn fixed_flow_turbulent_channel_requires_more_pressure_and_retains_si_pressure_work() {
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(3.0, 1.0, 6, 3));
    let system = FluxSystem::new(&mesh);
    let inlet = inlet(&mesh);
    let zero = |_: [f64; 2]| [0.0; 2];
    let traction = |_: [f64; 2], _: [f64; 2]| [0.0; 2];
    let boundaries: Vec<_> = mesh
        .edges
        .iter()
        .enumerate()
        .filter(|(edge, _)| mesh.is_boundary(*edge))
        .map(|(edge, face)| {
            if face.normal[0].abs() > 0.5 {
                Boundary::Open {
                    edge,
                    traction: &traction,
                    backflow: OpenBoundaryPolicy::NaturalTraction,
                }
            } else {
                Boundary::Velocity { edge, value: &zero }
            }
        })
        .collect();
    let initial = vec![0.0; system.n_u];
    let widths = vec![0.4; mesh.tris.len()];
    let requested = -0.3;
    let mut normalized: Option<Vec<f64>> = None;
    for rho in [2.0, 7.0] {
        let reference = transport(rho, 0.2);
        let owner = MaterialFluxSystem::new(&system, &reference, numerics()).unwrap();
        let mut required_pressure = Vec::new();
        for coefficient in [0.0, 0.6] {
            let solved = owner
                .solve_navier_stokes_mixed_rheology_sgs_flow_control(
                    MaterialSgsConfig {
                        model: SgsModel::Smagorinsky { cs: coefficient },
                        filter_widths_m: &widths,
                    },
                    &zero,
                    &boundaries,
                    &[&inlet],
                    &[requested],
                    &initial,
                    controls(),
                    1e-12,
                    |_, _| Ok(reference),
                )
                .unwrap();
            let pressure = solved.pressure_offsets_pa[0];
            assert!(pressure > 0.0 && -pressure * requested > 0.0);
            close(solved.outward_volume_rates_m2_s[0], requested, 2e-11);
            let solution = &solved.flow.flow.step.solution;
            let integrated: f64 = inlet
                .iter()
                .map(|&edge| mesh.edges[edge].len * solution.discrete().x[2 * edge])
                .sum();
            close(integrated, requested, 2e-11);
            let low = inlet
                .iter()
                .map(|&e| solution.discrete().x[2 * e])
                .fold(f64::INFINITY, f64::min);
            let high = inlet
                .iter()
                .map(|&e| solution.discrete().x[2 * e])
                .fold(f64::NEG_INFINITY, f64::max);
            assert!(high - low > 0.01, "the port profile must remain resolved");
            let load = solution.boundary_loads(&inlet, [0.0; 2], &[None]).unwrap();
            close(load.stress.force_n_per_m[0], -pressure, 2e-7);
            close(load.stress.force_n_per_m[1], 0.0, 2e-7);
            if coefficient == 0.0 {
                assert_eq!(solved.flow.modeled_transfer_w_per_m, 0.0);
            } else {
                assert!(solved.flow.modeled_transfer_w_per_m > 0.0);
            }
            required_pressure.push(pressure / rho);
        }
        assert!(
            required_pressure[1] > 1.05 * required_pressure[0],
            "SGS must change the required pressure at fixed flow: {required_pressure:?}"
        );
        if let Some(expected) = &normalized {
            for (actual, expected) in required_pressure.iter().zip(expected) {
                close(*actual, *expected, 2e-8);
            }
        } else {
            normalized = Some(required_pressure);
        }
    }
    assert_eq!(initial, vec![0.0; system.n_u]);
}

#[test]
fn zero_model_matches_original_bdf2_flow_control_and_refusals_keep_histories() {
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 0.7, 2, 2));
    let system = FluxSystem::new(&mesh);
    let inlet = inlet(&mesh);
    let outlet: Vec<_> = mesh
        .edges
        .iter()
        .enumerate()
        .filter(|(edge, face)| mesh.is_boundary(*edge) && face.normal[0] > 0.5)
        .map(|(edge, _)| edge)
        .collect();
    let zero = |_: [f64; 2]| [0.0; 2];
    let traction = |_: [f64; 2], _: [f64; 2]| [0.0; 2];
    let boundaries: Vec<_> = mesh
        .edges
        .iter()
        .enumerate()
        .filter(|(edge, _)| mesh.is_boundary(*edge))
        .map(|(edge, face)| {
            if face.normal[0].abs() > 0.5 {
                Boundary::Open {
                    edge,
                    traction: &traction,
                    backflow: OpenBoundaryPolicy::NaturalTraction,
                }
            } else {
                Boundary::Velocity { edge, value: &zero }
            }
        })
        .collect();
    let previous = vec![0.0; system.n_u];
    let older = vec![0.0; system.n_u];
    let saved = [previous.clone(), older.clone()];
    let reference = transport(3.0, 0.2);
    let owner = MaterialFluxSystem::new(&system, &reference, numerics()).unwrap();
    let widths = vec![0.4; mesh.tris.len()];
    let baseline = owner
        .bdf2_navier_stokes_mixed_rheology_variable_step_flow_control(
            &zero,
            &boundaries,
            &[&inlet],
            &[-0.15],
            &previous,
            Some((&older, 0.08)),
            0.05,
            controls(),
            1e-12,
            |_, _| Ok(reference),
        )
        .unwrap();
    for model in [
        SgsModel::Smagorinsky { cs: 0.0 },
        SgsModel::Wale { cw: 0.0 },
        SgsModel::Vreman { cv: 0.0 },
    ] {
        let solved = owner
            .bdf2_navier_stokes_mixed_rheology_sgs_variable_step_flow_control(
                MaterialSgsConfig {
                    model,
                    filter_widths_m: &widths,
                },
                &zero,
                &boundaries,
                &[&inlet],
                &[-0.15],
                &previous,
                Some((&older, 0.08)),
                0.05,
                controls(),
                1e-12,
                |_, _| Ok(reference),
            )
            .unwrap();
        assert_eq!(
            solved.flow.flow.step.solution.discrete().x,
            baseline.flow.step.solution.discrete().x
        );
        assert_eq!(solved.pressure_offsets_pa, baseline.pressure_offsets_pa);
        assert_eq!(
            solved.outward_volume_rates_m2_s,
            baseline.outward_volume_rates_m2_s
        );
        assert_eq!(solved.flow.modeled_transfer_w_per_m, 0.0);
    }
    let calls = Cell::new(0);
    let config = MaterialSgsConfig {
        model: SgsModel::Smagorinsky { cs: 0.5 },
        filter_widths_m: &widths,
    };
    let refused = owner.bdf2_navier_stokes_mixed_rheology_sgs_variable_step_flow_control(
        config,
        &zero,
        &boundaries,
        &[&inlet],
        &[-0.15],
        &previous,
        Some((&older, 0.01)),
        0.05,
        controls(),
        1e-12,
        |_, _| {
            calls.set(calls.get() + 1);
            Ok(reference)
        },
    );
    assert!(matches!(
        refused,
        Err(MaterialSgsError::Material(MaterialFluxError::Boundary(_)))
    ));
    assert_eq!(
        calls.get(),
        0,
        "temporal admission must precede source queries"
    );
    let refused = owner.bdf1_navier_stokes_mixed_rheology_sgs_flow_control(
        config,
        &zero,
        &boundaries,
        &[&inlet, &outlet],
        &[-0.15, 0.15],
        &previous,
        0.05,
        controls(),
        1e-12,
        |_, _| Ok(reference),
    );
    assert!(matches!(
        refused,
        Err(MaterialSgsError::Material(MaterialFluxError::Boundary(_)))
    ));
    let refused = owner.bdf1_navier_stokes_mixed_rheology_sgs_flow_control(
        config,
        &zero,
        &boundaries,
        &[&inlet],
        &[-0.15],
        &previous,
        0.05,
        controls(),
        1e-12,
        |point, _| Ok(transport(if point.cell == 1 { 4.0 } else { 3.0 }, 0.2)),
    );
    assert!(matches!(refused, Err(MaterialSgsError::Material(_))));
    let refused = owner.bdf1_navier_stokes_mixed_rheology_sgs_flow_control(
        config,
        &zero,
        &boundaries,
        &[&inlet],
        &[-0.15],
        &previous,
        0.05,
        NavierStokesControls {
            maximum_iterations: 1,
            ..controls()
        },
        1e-12,
        |_, _| Ok(reference),
    );
    assert!(matches!(
        refused,
        Err(MaterialSgsError::Material(MaterialFluxError::Boundary(_)))
    ));
    assert_eq!(previous, saved[0]);
    assert_eq!(older, saved[1]);
}
