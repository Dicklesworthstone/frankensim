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
