#![cfg(feature = "continuum")]
//! G1/G3: shear-dependent flow through open traction boundaries.

use fs_flux::material::MaterialStokesBoundary;
use fs_flux::ns::boundary::{
    StokesBoundary, StokesBoundaryError, StokesPressureReference, StokesRheologyControls,
};
use fs_flux::{
    FluxParams, FluxRheologyError, FluxSystem, MaterialFluxError, MaterialFluxNumerics,
    MaterialFluxRheologyPoint, MaterialFluxSystem, TriMesh,
};
use fs_material::MaterialError;
use fs_material::fluid_transport::{FluidTransport, FluidTransportState};
use fs_solid::Mesh2;

fn transport(rho: f64, mu: f64) -> FluidTransportState {
    FluidTransportState::try_new(rho, 1000.0, 0.5, mu).unwrap()
}

#[test]
fn g1_open_shear_thinning_channel_matches_independent_profile() {
    // nu=nu0/(1+lambda*gamma) has strictly increasing shear stress.
    // Independently invert tau=G*z: gamma=G*z/(nu0-lambda*G*z),
    // then integrate from z to the no-slip wall. No imposed inlet velocity.
    let nu0: f64 = 0.4;
    let lambda: f64 = 0.5;
    let gradient: f64 = 0.4;
    let half: f64 = 0.5;
    let b = lambda * gradient;
    let exact = |point: [f64; 2]| {
        let z = (point[1] - half).abs();
        [
            -(half - z) / lambda
                - nu0 / (lambda * lambda * gradient) * ((nu0 - b * half) / (nu0 - b * z)).ln(),
            0.0,
        ]
    };
    let exact_flux = 2.0
        * gradient
        * (-half * half / (2.0 * b)
            - nu0 * half / (b * b)
            - nu0 * nu0 / (b * b * b) * ((nu0 - b * half) / nu0).ln());
    let mut errors = Vec::new();
    for n in [4, 8] {
        let mesh = TriMesh::from_mesh2(&Mesh2::triangles(2.0, 1.0, 2 * n, n));
        let sys = FluxSystem::new(&mesh);
        let params = FluxParams {
            nu: f64::NAN,
            tol: 1e-12,
            ..FluxParams::default()
        };
        let controls = StokesRheologyControls {
            maximum_iterations: 60,
            relative_tolerance: 1e-9,
            relaxation: 1.0,
        };
        let wall = |_: [f64; 2]| [0.0; 2];
        let traction = |point: [f64; 2], normal: [f64; 2]| {
            let p = 5.0 - gradient * point[0];
            let shear = gradient * (half - point[1]);
            [
                -p * normal[0] + shear * normal[1],
                shear * normal[0] - p * normal[1],
            ]
        };
        let boundaries: Vec<_> = mesh
            .edges
            .iter()
            .enumerate()
            .filter(|(e, _)| mesh.is_boundary(*e))
            .map(|(edge, g)| {
                if g.normal[0].abs() > 0.5 {
                    StokesBoundary::TotalTraction {
                        edge,
                        value: &traction,
                    }
                } else {
                    StokesBoundary::Velocity { edge, value: &wall }
                }
            })
            .collect();
        let initial = vec![0.0; sys.n_u];
        let solved = sys
            .solve_stokes_mixed_rheology(
                params,
                &|_, _, gamma| nu0 / (1.0 + lambda * gamma),
                &|_| [0.0; 2],
                &boundaries,
                &initial,
                controls,
            )
            .unwrap();
        assert!(solved.solution.field.converged(controls.relative_tolerance));
        assert!(solved.relative_update <= controls.relative_tolerance);
        assert!(solved.nonlinear_residual <= controls.relative_tolerance);
        assert!(solved.iterations <= controls.maximum_iterations);
        assert_eq!(
            solved.solution.pressure_reference,
            StokesPressureReference::TractionDatum
        );
        let (error, div) = sys.velocity_error(&solved.solution.field.x, &exact);
        let mut inlet = 0.0;
        let mut outlet = 0.0;
        for (edge, g) in mesh.edges.iter().enumerate() {
            if !mesh.is_boundary(edge) {
                continue;
            }
            let q = g.len * solved.solution.field.x[2 * edge];
            if g.normal[0] < -0.5 {
                inlet -= q;
            } else if g.normal[0] > 0.5 {
                outlet += q;
            }
        }
        let mean = solved
            .solution
            .pressure
            .iter()
            .zip(&mesh.areas)
            .map(|(p, a)| p * a)
            .sum::<f64>()
            / mesh.areas.iter().sum::<f64>();
        eprintln!(
            "nonlinear open channel n={n}: iterations={}, L2={error:e}, div={div:e}, Q={outlet:e}, exactQ={exact_flux:e}, meanp={mean:e}, residual={:e}",
            solved.iterations, solved.nonlinear_residual
        );
        assert!(div < 1e-8 && (inlet - outlet).abs() < 1e-9);
        assert!((mean - 4.6).abs() < 0.01);
        if n == 8 {
            assert!(error < 0.004, "fine velocity error {error:e}");
            assert!(
                (outlet / exact_flux - 1.0).abs() < 0.03,
                "fine flow rate {outlet:e}"
            );
        } else {
            let frozen = sys
                .assemble_stokes_mixed_variable_viscosity(
                    params,
                    &|_, _| nu0,
                    &|_| [0.0; 2],
                    &boundaries,
                )
                .unwrap()
                .solve()
                .unwrap();
            assert!(sys.velocity_error(&frozen.field.x, &exact).0 > 2.0 * error);
        }
        assert_eq!(initial, vec![0.0; sys.n_u]);
        errors.push(error);
    }
    let ratio = errors[0] / errors[1];
    assert!(
        (3.48..4.60).contains(&ratio),
        "second-order errors={errors:?}, ratio={ratio}"
    );
}

struct Response {
    state: FluidTransportState,
    point: MaterialFluxRheologyPoint,
    temperature_k: f64,
    query: usize,
}
impl FluidTransport for Response {
    fn transport_state(&self) -> Result<FluidTransportState, MaterialError> {
        Ok(self.state)
    }
}

#[test]
fn g1_material_open_thickening_flow_retains_sources_stress_and_pressure() {
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 2, 2));
    let sys = FluxSystem::new(&mesh);
    // mu=K(T)*(1+gamma^2). Exact shear=1Pa requires gamma=3 below
    // and gamma=1 above y=.5. The undamped local fixed point is unstable
    // in the lower layer; bounded relaxation must close the actual stress.
    let velocity = |p: [f64; 2]| [0.2 + 3.0 * p[1].min(0.5) + (p[1] - 0.5).max(0.0), 0.0];
    let traction = |_: [f64; 2], n: [f64; 2]| [-7.0 * n[0] + n[1], n[0] - 7.0 * n[1]];
    let boundaries: Vec<_> = mesh
        .edges
        .iter()
        .enumerate()
        .filter(|(e, _)| mesh.is_boundary(*e))
        .map(|(edge, g)| {
            if g.normal[0] < -0.5 {
                MaterialStokesBoundary::Velocity {
                    edge,
                    value: &velocity,
                }
            } else {
                MaterialStokesBoundary::TotalTraction {
                    edge,
                    value: &traction,
                }
            }
        })
        .collect();
    for rho in [2.0, 7.0] {
        let reference = transport(rho, 0.3);
        let flow = MaterialFluxSystem::new(
            &sys,
            &reference,
            MaterialFluxNumerics {
                linear_tolerance: 1e-12,
                ..MaterialFluxNumerics::default()
            },
        )
        .unwrap();
        let initial = vec![0.0; sys.n_u];
        let mut queries = 0;
        let controls = StokesRheologyControls {
            maximum_iterations: 80,
            relative_tolerance: 1e-9,
            relaxation: 0.5,
        };
        let solved = flow
            .solve_stokes_mixed_rheology(
                &|_| [0.0; 2],
                &boundaries,
                &initial,
                controls,
                |point, reference| {
                    queries += 1;
                    let lower = point.position_m[1] < 0.5;
                    let k = if lower { 1.0 / 30.0 } else { 0.5 };
                    Ok(Response {
                        state: transport(
                            reference.density_kg_m3(),
                            k * (1.0 + point.shear_rate_per_s * point.shear_rate_per_s),
                        ),
                        point,
                        temperature_k: if lower { 330.0 } else { 300.0 },
                        query: queries,
                    })
                },
            )
            .unwrap();
        let (error, div) = sys.velocity_error(&solved.discrete().x, &velocity);
        let rates = sys.shear_rates(&solved.discrete().x);
        let mut stress_error = 0.0_f64;
        for (cell, sample) in solved.samples().iter().enumerate() {
            assert_eq!(sample.point(), sample.material().point);
            assert_eq!(
                sample.point().shear_rate_per_s.to_bits(),
                rates[cell].to_bits()
            );
            assert_eq!(sample.point().position_m, mesh.centroids[cell]);
            assert_eq!(sample.transport().density_kg_m3(), rho);
            assert_eq!(
                sample.material().temperature_k,
                if mesh.centroids[cell][1] < 0.5 {
                    330.0
                } else {
                    300.0
                }
            );
            assert!(sample.material().query <= queries);
            stress_error = stress_error
                .max((sample.transport().dynamic_viscosity_pa_s() * rates[cell] - 1.0).abs());
            let expected_heat = if mesh.centroids[cell][1] < 0.5 {
                3.0
            } else {
                1.0
            };
            assert!(
                (solved.viscous_dissipation().cell_heating_w_m3()[cell] - expected_heat).abs()
                    < 4e-7
            );
        }
        let pressure_error = solved
            .pressure_pa()
            .iter()
            .fold(0.0_f64, |e, p| e.max((p - 7.0).abs()));
        // Independent integral of the two affine velocity pieces.
        let exact_energy = 0.5
            * rho
            * ((1.7_f64.powi(3) - 0.2_f64.powi(3)) / 9.0
                + (2.2_f64.powi(3) - 1.7_f64.powi(3)) / 3.0);
        eprintln!(
            "material open thickening rho={rho}: iterations={}, L2={error:e}, div={div:e}, stress={stress_error:e}, pressure={pressure_error:e}, residual={:e}, queries={queries}",
            solved.iterations(),
            solved.nonlinear_residual()
        );
        assert!(error < 1e-7 && div < 1e-9 && stress_error < 1e-7 && pressure_error < 1e-7);
        assert!((solved.kinetic_energy_j_per_m() / exact_energy - 1.0).abs() < 1e-7);
        // tau=1 Pa over area one: half at gamma=3/s, half at gamma=1/s.
        assert!((solved.viscous_dissipation().total_power_w_per_m() - 2.0).abs() < 4e-7);
        assert_eq!(
            solved.pressure_reference(),
            StokesPressureReference::TractionDatum
        );
        assert!(queries <= mesh.tris.len() * (solved.iterations() + 1));
        assert_eq!(initial, vec![0.0; sys.n_u]);
    }
}

#[test]
fn material_mixed_refusals_preserve_source_error_and_input_state() {
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 2, 2));
    let sys = FluxSystem::new(&mesh);
    let reference = transport(2.0, 1.0);
    let flow = MaterialFluxSystem::new(&sys, &reference, MaterialFluxNumerics::default()).unwrap();
    let initial = vec![0.0; sys.n_u];
    let velocity = |p: [f64; 2]| [p[1], 0.0];
    let traction = |_: [f64; 2], n: [f64; 2]| [-7.0 * n[0] + n[1], n[0] - 7.0 * n[1]];
    let boundaries: Vec<_> = mesh
        .edges
        .iter()
        .enumerate()
        .filter(|(e, _)| mesh.is_boundary(*e))
        .map(|(edge, g)| {
            if g.normal[0] < -0.5 {
                MaterialStokesBoundary::Velocity {
                    edge,
                    value: &velocity,
                }
            } else {
                MaterialStokesBoundary::TotalTraction {
                    edge,
                    value: &traction,
                }
            }
        })
        .collect();
    let refusal = MaterialError::State {
        what: "outside original material support".into(),
    };
    let controls = StokesRheologyControls::default();
    let mut queries = 0;
    let failed =
        flow.solve_stokes_mixed_rheology(&|_| [0.0; 2], &boundaries, &initial, controls, |_, _| {
            queries += 1;
            Err::<FluidTransportState, _>(refusal.clone())
        });
    assert!(matches!(failed,Err(MaterialFluxError::Material(e)) if e==refusal));
    assert_eq!(queries, 1);
    queries = 0;
    let failed = flow.solve_stokes_mixed_rheology(
        &|_| [0.0; 2],
        &boundaries,
        &initial,
        controls,
        |point, _| {
            queries += 1;
            if point.shear_rate_per_s > 1e-6 {
                Err(refusal.clone())
            } else {
                Ok(reference)
            }
        },
    );
    assert!(matches!(failed,Err(MaterialFluxError::Material(e)) if e==refusal));
    assert!(queries > mesh.tris.len() && queries <= 2 * mesh.tris.len());
    let failed =
        flow.solve_stokes_mixed_rheology(&|_| [0.0; 2], &boundaries, &initial, controls, |_, _| {
            Ok(transport(3.0, 1.0))
        });
    assert!(matches!(failed, Err(MaterialFluxError::Invalid(_))));
    let failed = flow.solve_stokes_mixed_rheology(
        &|_| [0.0; 2],
        &boundaries,
        &initial,
        StokesRheologyControls {
            maximum_iterations: 1,
            ..controls
        },
        |_, _| Ok(reference),
    );
    assert!(matches!(
        failed,
        Err(MaterialFluxError::Boundary(StokesBoundaryError::Rheology(
            FluxRheologyError::NonlinearConvergence { iterations: 1, .. }
        )))
    ));
    assert_eq!(initial, vec![0.0; sys.n_u]);
}

#[test]
fn constant_law_mixed_rheology_reduces_to_the_checked_linear_flow() {
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 2, 2));
    let sys = FluxSystem::new(&mesh);
    let params = FluxParams {
        nu: 0.7,
        ..FluxParams::default()
    };
    let velocity = |p: [f64; 2]| [0.2 + 1.3 * p[1], 0.0];
    let traction = |_: [f64; 2], n: [f64; 2]| [-7.0 * n[0] + 0.91 * n[1], 0.91 * n[0] - 7.0 * n[1]];
    let boundaries: Vec<_> = mesh
        .edges
        .iter()
        .enumerate()
        .filter(|(e, _)| mesh.is_boundary(*e))
        .map(|(edge, g)| {
            if g.normal[0] < -0.5 {
                StokesBoundary::Velocity {
                    edge,
                    value: &velocity,
                }
            } else {
                StokesBoundary::TotalTraction {
                    edge,
                    value: &traction,
                }
            }
        })
        .collect();
    let linear = sys
        .assemble_stokes_mixed(params, &|_| [0.0; 2], &boundaries)
        .unwrap()
        .solve()
        .unwrap();
    let solved = sys
        .solve_stokes_mixed_rheology(
            params,
            &|_, _, _| params.nu,
            &|_| [0.0; 2],
            &boundaries,
            &vec![0.0; sys.n_u],
            StokesRheologyControls::default(),
        )
        .unwrap();
    assert_eq!(solved.iterations, 2);
    assert_eq!(
        solved.solution.pressure_reference,
        linear.pressure_reference
    );
    assert!(
        linear
            .field
            .x
            .iter()
            .zip(&solved.solution.field.x)
            .all(|(a, b)| a.to_bits() == b.to_bits())
    );
}

#[test]
fn prescribed_shear_thinning_flow_recovers_independent_required_pressure() {
    // Invert tau=nu0*gamma/(1+lambda*gamma) for tau=G*|y-h/2|
    // and integrate the resulting velocity and volume rate analytically.
    // The port fixes only that integral; its profile and pressure are unknown.
    let nu0: f64 = 0.4;
    let lambda: f64 = 0.5;
    let gradient: f64 = 0.4;
    let half: f64 = 0.5;
    let b = lambda * gradient;
    let exact = |point: [f64; 2]| {
        let z = (point[1] - half).abs();
        [
            -(half - z) / lambda
                - nu0 / (lambda * lambda * gradient) * ((nu0 - b * half) / (nu0 - b * z)).ln(),
            0.0,
        ]
    };
    let exact_flux = 2.0
        * gradient
        * (-half * half / (2.0 * b)
            - nu0 * half / (b * b)
            - nu0 * nu0 / (b * b * b) * ((nu0 - b * half) / nu0).ln());
    let mut pressure_errors = Vec::new();
    for n in [4, 8] {
        let mesh = TriMesh::from_mesh2(&Mesh2::triangles(2.0, 1.0, 2 * n, n));
        let sys = FluxSystem::new(&mesh);
        let params = FluxParams {
            nu: f64::NAN,
            tol: 1e-12,
            ..FluxParams::default()
        };
        let controls = StokesRheologyControls {
            maximum_iterations: 80,
            relative_tolerance: 1e-9,
            relaxation: 0.65,
        };
        let wall = |_: [f64; 2]| [0.0; 2];
        let traction = |point: [f64; 2], normal: [f64; 2]| {
            let shear = gradient * (half - point[1]);
            [
                -4.2 * normal[0] + shear * normal[1],
                shear * normal[0] - 4.2 * normal[1],
            ]
        };
        let inlet: Vec<_> = mesh
            .edges
            .iter()
            .enumerate()
            .filter(|(edge, g)| mesh.is_boundary(*edge) && g.normal[0] < -0.5)
            .map(|(edge, _)| edge)
            .collect();
        let boundaries: Vec<_> = mesh
            .edges
            .iter()
            .enumerate()
            .filter(|(edge, _)| mesh.is_boundary(*edge))
            .map(|(edge, g)| {
                if g.normal[0].abs() > 0.5 {
                    StokesBoundary::TotalTraction {
                        edge,
                        value: &traction,
                    }
                } else {
                    StokesBoundary::Velocity { edge, value: &wall }
                }
            })
            .collect();
        let initial = vec![0.0; sys.n_u];
        let controlled = sys
            .solve_stokes_mixed_rheology_flow_control(
                params,
                &|_, _, gamma| nu0 / (1.0 + lambda * gamma),
                &|_| [0.0; 2],
                &boundaries,
                &[&inlet],
                &[-exact_flux],
                &initial,
                controls,
                1e-12,
            )
            .unwrap();
        let solved = &controlled.flow;
        assert!(solved.solution.field.converged(controls.relative_tolerance));
        assert!(solved.relative_update <= controls.relative_tolerance);
        assert!(solved.nonlinear_residual <= controls.relative_tolerance);
        assert!(controlled.maximum_flow_residual_m2_s < 2e-12);
        assert!((controlled.outward_volume_rates_m2_s[0] + exact_flux).abs() < 2e-12);
        let (error, divergence) = sys.velocity_error(&solved.solution.field.x, &exact);
        assert!(divergence < 1e-8);
        let pressure_error = (controlled.pressure_offsets_m2_s2[0] - 0.8).abs();
        pressure_errors.push(pressure_error);
        if n == 8 {
            assert!(error < 0.004, "resolved velocity error={error:e}");
            assert!(
                pressure_error < 0.04,
                "required pressure error={pressure_error:e}"
            );
        } else {
            // Freezing the zero-shear viscosity gives a different required
            // pressure even though it also enforces the same total flow.
            let frozen = sys
                .assemble_stokes_mixed_variable_viscosity(
                    params,
                    &|_, _| nu0,
                    &|_| [0.0; 2],
                    &boundaries,
                )
                .unwrap()
                .prepare_flow_control(&[&inlet])
                .unwrap()
                .solve(&[-exact_flux], 1e-12)
                .unwrap();
            assert!(
                (frozen.pressure_offsets_m2_s2[0] - 0.8).abs() > 2.0 * pressure_error,
                "a frozen material must not satisfy the nonlinear pressure reference"
            );
        }
        assert_eq!(initial, vec![0.0; sys.n_u]);
    }
    assert!(
        pressure_errors[1] < 0.7 * pressure_errors[0],
        "required pressure must approach the independent channel law: {pressure_errors:?}"
    );
}

#[test]
fn material_flow_control_closes_layered_shear_and_required_pressure_together() {
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 2, 2));
    let sys = FluxSystem::new(&mesh);
    // mu=K(T)*(1+gamma²), with gamma=3 below and gamma=1 above:
    // the continuous shear stress is 1 Pa and the integral flow is 1.45 m²/s.
    let velocity = |p: [f64; 2]| [0.2 + 3.0 * p[1].min(0.5) + (p[1] - 0.5).max(0.0), 0.0];
    let traction = |_: [f64; 2], n: [f64; 2]| [-7.0 * n[0] + n[1], n[0] - 7.0 * n[1]];
    let inlet: Vec<_> = mesh
        .edges
        .iter()
        .enumerate()
        .filter(|(edge, g)| mesh.is_boundary(*edge) && g.normal[0] < -0.5)
        .map(|(edge, _)| edge)
        .collect();
    let boundaries: Vec<_> = mesh
        .edges
        .iter()
        .enumerate()
        .filter(|(edge, _)| mesh.is_boundary(*edge))
        .map(|(edge, g)| {
            if g.normal[1].abs() > 0.5 {
                MaterialStokesBoundary::Velocity {
                    edge,
                    value: &velocity,
                }
            } else {
                MaterialStokesBoundary::TotalTraction {
                    edge,
                    value: &traction,
                }
            }
        })
        .collect();
    let initial = vec![0.0; sys.n_u];
    let controls = StokesRheologyControls {
        maximum_iterations: 100,
        relative_tolerance: 1e-9,
        relaxation: 0.5,
    };
    for density in [2.0, 7.0] {
        let reference = transport(density, 0.3);
        let flow = MaterialFluxSystem::new(
            &sys,
            &reference,
            MaterialFluxNumerics {
                linear_tolerance: 1e-12,
                ..MaterialFluxNumerics::default()
            },
        )
        .unwrap();
        let mut queries = 0;
        let controlled = flow
            .solve_stokes_mixed_rheology_flow_control(
                &|_| [-1.7, 0.0],
                &boundaries,
                &[&inlet],
                &[-1.45],
                &initial,
                controls,
                1e-12,
                |point, reference| {
                    queries += 1;
                    let lower = point.position_m[1] < 0.5;
                    let k = if lower { 1.0 / 30.0 } else { 0.5 };
                    Ok(Response {
                        state: transport(
                            reference.density_kg_m3(),
                            k * (1.0 + point.shear_rate_per_s.powi(2)),
                        ),
                        point,
                        temperature_k: if lower { 330.0 } else { 300.0 },
                        query: queries,
                    })
                },
            )
            .unwrap();
        let solved = &controlled.flow;
        let (error, divergence) = sys.velocity_error(&solved.discrete().x, &velocity);
        assert!(error < 1e-7 && divergence < 1e-9, "{error}, {divergence}");
        assert!((controlled.pressure_offsets_pa[0] - 1.7).abs() < 1e-7);
        assert!(controlled.maximum_flow_residual_m2_s < 3e-12);
        assert!(solved.discrete().converged(controls.relative_tolerance));
        assert!(solved.relative_update() <= controls.relative_tolerance);
        assert!(solved.nonlinear_residual() <= controls.relative_tolerance);
        let rates = sys.shear_rates(&solved.discrete().x);
        for (cell, sample) in solved.samples().iter().enumerate() {
            assert_eq!(sample.point(), sample.material().point);
            assert_eq!(
                sample.point().shear_rate_per_s.to_bits(),
                rates[cell].to_bits()
            );
            assert_eq!(sample.transport().density_kg_m3(), density);
            assert_eq!(sample.point().position_m, mesh.centroids[cell]);
            assert!(sample.material().query <= queries);
            assert_eq!(
                sample.material().temperature_k,
                if mesh.centroids[cell][1] < 0.5 {
                    330.0
                } else {
                    300.0
                }
            );
            assert!((sample.transport().dynamic_viscosity_pa_s() * rates[cell] - 1.0).abs() < 1e-7);
            assert!(
                (solved.pressure_pa()[cell] - 7.0 - 1.7 * (1.0 - mesh.centroids[cell][0])).abs()
                    < 1e-7
            );
        }
        let exact_energy = 0.5
            * density
            * ((1.7_f64.powi(3) - 0.2_f64.powi(3)) / 9.0
                + (2.2_f64.powi(3) - 1.7_f64.powi(3)) / 3.0);
        assert!((solved.kinetic_energy_j_per_m() / exact_energy - 1.0).abs() < 1e-7);
        assert!((solved.viscous_dissipation().total_power_w_per_m() - 2.0).abs() < 4e-7);
        let load = solved.boundary_loads(&inlet, [0.0; 2], &[None]).unwrap();
        assert!((load.stress.force_n_per_m[0] + 8.7).abs() < 1e-7);
        assert!((load.stress.force_n_per_m[1] - 1.0).abs() < 1e-7);
        assert!((load.stress.moment_n - 4.35).abs() < 1e-7);
        assert!(queries <= mesh.tris.len() * (solved.iterations() + 1));
        assert_eq!(initial, vec![0.0; sys.n_u]);
    }
}

#[test]
fn flow_control_refuses_source_and_nonlinear_failures_without_partial_state() {
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 2, 2));
    let sys = FluxSystem::new(&mesh);
    let reference = transport(2.0, 1.0);
    let flow = MaterialFluxSystem::new(&sys, &reference, MaterialFluxNumerics::default()).unwrap();
    let initial = vec![0.0; sys.n_u];
    let velocity = |p: [f64; 2]| [p[1], 0.0];
    let traction = |_: [f64; 2], n: [f64; 2]| [-7.0 * n[0] + n[1], n[0] - 7.0 * n[1]];
    let inlet: Vec<_> = mesh
        .edges
        .iter()
        .enumerate()
        .filter(|(edge, g)| mesh.is_boundary(*edge) && g.normal[0] < -0.5)
        .map(|(edge, _)| edge)
        .collect();
    let boundaries: Vec<_> = mesh
        .edges
        .iter()
        .enumerate()
        .filter(|(edge, _)| mesh.is_boundary(*edge))
        .map(|(edge, g)| {
            if g.normal[1].abs() > 0.5 {
                MaterialStokesBoundary::Velocity {
                    edge,
                    value: &velocity,
                }
            } else {
                MaterialStokesBoundary::TotalTraction {
                    edge,
                    value: &traction,
                }
            }
        })
        .collect();
    let refusal = MaterialError::State {
        what: "outside original material support".into(),
    };
    let controls = StokesRheologyControls::default();
    let mut queries = 0;
    let failed = flow.solve_stokes_mixed_rheology_flow_control(
        &|_| [0.0; 2],
        &boundaries,
        &[&inlet],
        &[-0.5],
        &initial,
        controls,
        1e-12,
        |_, _| {
            queries += 1;
            Err::<FluidTransportState, _>(refusal.clone())
        },
    );
    assert!(matches!(failed, Err(MaterialFluxError::Material(error)) if error == refusal));
    assert_eq!(queries, 1);
    queries = 0;
    let failed = flow.solve_stokes_mixed_rheology_flow_control(
        &|_| [0.0; 2],
        &boundaries,
        &[&inlet],
        &[-0.5],
        &initial,
        controls,
        1e-12,
        |point, _| {
            queries += 1;
            if point.shear_rate_per_s > 1e-6 {
                Err(refusal.clone())
            } else {
                Ok(reference)
            }
        },
    );
    assert!(matches!(failed, Err(MaterialFluxError::Material(error)) if error == refusal));
    assert!(queries > mesh.tris.len() && queries <= 2 * mesh.tris.len());
    let failed = flow.solve_stokes_mixed_rheology_flow_control(
        &|_| [0.0; 2],
        &boundaries,
        &[&inlet],
        &[-0.5],
        &initial,
        controls,
        1e-12,
        |_, _| Ok(transport(3.0, 1.0)),
    );
    assert!(matches!(failed, Err(MaterialFluxError::Invalid(_))));
    let failed = flow.solve_stokes_mixed_rheology_flow_control(
        &|_| [0.0; 2],
        &boundaries,
        &[&inlet],
        &[-0.5],
        &initial,
        StokesRheologyControls {
            maximum_iterations: 1,
            ..controls
        },
        1e-12,
        |_, _| Ok(reference),
    );
    assert!(matches!(
        failed,
        Err(MaterialFluxError::Boundary(StokesBoundaryError::Rheology(
            FluxRheologyError::NonlinearConvergence { iterations: 1, .. }
        )))
    ));
    queries = 0;
    let failed = flow.solve_stokes_mixed_rheology_flow_control(
        &|_| [0.0; 2],
        &boundaries,
        &[&inlet],
        &[f64::NAN],
        &initial,
        controls,
        1e-12,
        |_, _| {
            queries += 1;
            Ok(reference)
        },
    );
    assert!(matches!(
        failed,
        Err(MaterialFluxError::Boundary(
            StokesBoundaryError::InvalidInput(_)
        ))
    ));
    assert_eq!(queries, 0);
    assert_eq!(initial, vec![0.0; sys.n_u]);
}

#[test]
#[allow(clippy::too_many_lines)]
fn g1_mixed_newton_resolves_pressure_driven_thickening_and_datum_shifts() {
    use fs_flux::{FluxNewtonControls, FluxViscosityTangent};

    // Exact affine extension between two pressure-loaded faces. The other
    // faces are stationary symmetry walls, not imposed tangential velocities.
    // u=(x,-y), gamma=2, nu=0.2*(1+gamma^2)=1, p=P.
    // Right/top reservoir pressures are P-2 and P+2. Their net work is 4,
    // exactly the bulk strain dissipation nu*gamma^2 on the unit square.
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 2, 2));
    let sys = FluxSystem::new(&mesh);
    let initial = vec![0.0; sys.n_u];
    let params = FluxParams {
        nu: f64::NAN,
        tol: 1e-12,
        ..FluxParams::default()
    };
    let controls = FluxNewtonControls {
        residual_tolerance: 2e-8,
        ..FluxNewtonControls::default()
    };
    let law = |_: usize, _: [f64; 2], rate: f64| FluxViscosityTangent {
        nu: 0.2 * (1.0 + rate * rate),
        derivative: 0.4 * rate,
    };
    let wall = |_: [f64; 2]| [0.0; 2];
    for pressure in [7.0, 100_007.0] {
        let traction = |_: [f64; 2], normal: [f64; 2]| {
            let reservoir = pressure + if normal[0] > 0.5 { -2.0 } else { 2.0 };
            [-reservoir * normal[0], -reservoir * normal[1]]
        };
        let boundaries: Vec<_> = mesh
            .edges
            .iter()
            .enumerate()
            .filter(|(edge, _)| mesh.is_boundary(*edge))
            .map(|(edge, geometry)| {
                if geometry.normal[0] > 0.5 || geometry.normal[1] > 0.5 {
                    StokesBoundary::TotalTraction {
                        edge,
                        value: &traction,
                    }
                } else {
                    StokesBoundary::NavierSlip {
                        edge,
                        value: &wall,
                        tangential_resistance_m_s: 0.0,
                    }
                }
            })
            .collect();
        if pressure == 7.0 {
            // With uniform strain, undamped Picard has the scalar map
            // a_next=5/(1+4*a^2), whose derivative at a=1 is -8/5.
            let picard = sys.solve_stokes_mixed_rheology(
                params,
                &|cell, point, rate| law(cell, point, rate).nu,
                &|_| [0.0; 2],
                &boundaries,
                &initial,
                StokesRheologyControls {
                    maximum_iterations: 12,
                    relative_tolerance: controls.residual_tolerance,
                    relaxation: 1.0,
                },
            );
            assert!(matches!(
                picard,
                Err(StokesBoundaryError::Rheology(
                    FluxRheologyError::NonlinearConvergence { .. }
                ))
            ));
        }
        let solved = sys
            .solve_stokes_mixed_rheology_newton(
                params,
                &law,
                &|_| [0.0; 2],
                &boundaries,
                &initial,
                controls,
            )
            .unwrap();
        assert!(solved.iterations > 1 && solved.iterations <= controls.maximum_iterations);
        assert!(solved.nonlinear_residual <= controls.residual_tolerance);
        assert!(solved.solution.field.converged(controls.residual_tolerance));
        assert_eq!(
            solved.solution.pressure_reference,
            StokesPressureReference::TractionDatum
        );
        let (error, divergence) =
            sys.velocity_error(&solved.solution.field.x, &|x| [x[0], -x[1]]);
        assert!(
            error < 2e-7 && divergence < 2e-7,
            "{error:e}, {divergence:e}"
        );
        for p in &solved.solution.pressure {
            assert!(
                (p - pressure).abs() < 2e-6,
                "pressure datum {p} versus {pressure}"
            );
        }
        let mut right_flux = 0.0;
        let mut top_flux = 0.0;
        for (edge, geometry) in mesh.edges.iter().enumerate() {
            if mesh.is_boundary(edge) {
                let flux = geometry.len * solved.solution.field.x[2 * edge];
                if geometry.normal[0] > 0.5 {
                    right_flux += flux;
                } else if geometry.normal[1] > 0.5 {
                    top_flux += flux;
                } else {
                    assert!(flux.abs() < 1e-8);
                }
            }
        }
        assert!((right_flux - 1.0).abs() < 2e-7);
        assert!((top_flux + 1.0).abs() < 2e-7);
        assert!((right_flux + top_flux).abs() < 1e-9);
        let rates = sys.shear_rates(&solved.solution.field.x);
        let bulk: f64 = rates
            .iter()
            .zip(&mesh.areas)
            .map(|(rate, area)| 0.2 * (1.0 + rate * rate) * rate * rate * area)
            .sum();
        let pressure_work = -(pressure - 2.0) * right_flux - (pressure + 2.0) * top_flux;
        assert!((bulk - 4.0).abs() < 3e-6);
        assert!((pressure_work - bulk).abs() < 3e-6);
        assert_eq!(initial, vec![0.0; sys.n_u]);
    }
}

#[test]
fn g1_mixed_newton_preserves_layer_stress_and_physical_navier_slip() {
    use fs_flux::{FluxNewtonControls, FluxViscosityTangent};

    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 2, 2));
    let sys = FluxSystem::new(&mesh);
    // nu=K*(1+gamma^2): gamma=3 below and gamma=1 above y=1/2
    // both give shear stress 3. A bottom resistance beta=12 gives slip=1/4.
    // The top shear work 3*(1/4+3/2+1/2)=6.75 splits into bulk=6 and slip=.75.
    let coefficients: Vec<f64> = mesh
        .tris
        .iter()
        .map(|tri| {
            let y = tri.iter().map(|&node| mesh.verts[node][1]).sum::<f64>() / 3.0;
            if y < 0.5 {
                0.1
            } else {
                1.5
            }
        })
        .collect();
    let law = |cell: usize, _: [f64; 2], rate: f64| FluxViscosityTangent {
        nu: coefficients[cell] * (1.0 + rate * rate),
        derivative: 2.0 * coefficients[cell] * rate,
    };
    let wall = |_: [f64; 2]| [0.0; 2];
    let traction =
        |_: [f64; 2], n: [f64; 2]| [-7.0 * n[0] + 3.0 * n[1], 3.0 * n[0] - 7.0 * n[1]];
    let boundaries: Vec<_> = mesh
        .edges
        .iter()
        .enumerate()
        .filter(|(edge, _)| mesh.is_boundary(*edge))
        .map(|(edge, geometry)| {
            if geometry.normal[1] < -0.5 {
                StokesBoundary::NavierSlip {
                    edge,
                    value: &wall,
                    tangential_resistance_m_s: 12.0,
                }
            } else {
                StokesBoundary::TotalTraction {
                    edge,
                    value: &traction,
                }
            }
        })
        .collect();
    let solved = sys
        .solve_stokes_mixed_rheology_newton(
            FluxParams {
                nu: f64::NAN,
                tol: 1e-12,
                ..FluxParams::default()
            },
            &law,
            &|_| [0.0; 2],
            &boundaries,
            &vec![0.0; sys.n_u],
            FluxNewtonControls::default(),
        )
        .unwrap();
    let exact = |point: [f64; 2]| {
        [
            0.25 + 3.0 * point[1].min(0.5) + (point[1] - 0.5).max(0.0),
            0.0,
        ]
    };
    let (error, divergence) = sys.velocity_error(&solved.solution.field.x, &exact);
    assert!(error < 2e-7 && divergence < 2e-7, "{error:e}, {divergence:e}");
    let rates = sys.shear_rates(&solved.solution.field.x);
    let mut bulk = 0.0;
    for (cell, &rate) in rates.iter().enumerate() {
        let expected = if coefficients[cell] == 0.1 { 3.0 } else { 1.0 };
        assert!((rate - expected).abs() < 3e-7);
        let nu = law(cell, [0.0; 2], rate).nu;
        assert!((nu * rate - 3.0).abs() < 2e-6);
        bulk += mesh.areas[cell] * nu * rate * rate;
    }
    let slip: f64 = solved
        .solution
        .reactions
        .edges()
        .iter()
        .map(|edge| edge.load([0.0; 2], None).unwrap().slip_dissipation)
        .sum();
    assert!((bulk - 6.0).abs() < 4e-6);
    assert!((slip - 0.75).abs() < 4e-7);
    assert!((bulk + slip - 6.75).abs() < 4e-6);
}

#[test]
fn mixed_newton_refuses_invalid_tangents_and_preserves_the_initial_state() {
    use fs_flux::ns::linearization::FluxLinearizationError;
    use fs_flux::{FluxNewtonControls, FluxNewtonError, FluxViscosityTangent};

    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 1, 1));
    let sys = FluxSystem::new(&mesh);
    let wall = |_: [f64; 2]| [0.0; 2];
    let traction = |_: [f64; 2], n: [f64; 2]| [-n[0], -n[1]];
    let boundaries: Vec<_> = mesh
        .edges
        .iter()
        .enumerate()
        .filter(|(edge, _)| mesh.is_boundary(*edge))
        .map(|(edge, geometry)| {
            if geometry.normal[1] < -0.5 {
                StokesBoundary::Velocity { edge, value: &wall }
            } else {
                StokesBoundary::TotalTraction {
                    edge,
                    value: &traction,
                }
            }
        })
        .collect();
    let initial = vec![0.0; sys.n_u];
    let refused = sys.solve_stokes_mixed_rheology_newton(
        FluxParams::default(),
        &|_, _, _| FluxViscosityTangent {
            nu: 1.0,
            derivative: f64::NAN,
        },
        &|_| [0.0; 2],
        &boundaries,
        &initial,
        FluxNewtonControls::default(),
    );
    assert!(matches!(
        refused,
        Err(FluxNewtonError::Linearization(
            FluxLinearizationError::Constitutive { .. }
        ))
    ));
    assert_eq!(initial, vec![0.0; sys.n_u]);
}
