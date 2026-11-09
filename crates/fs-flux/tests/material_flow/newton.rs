//! G1/G3/G4: original constitutive tangents reach the actual nonlinear solver.

use super::*;
use fs_flux::MaterialFluxViscosityTangent;
use fs_flux::ns::linearization::FluxViscosityTangent;
use fs_flux::ns::newton::{FluxNewtonControls, FluxNewtonError};
use fs_material::fluid_transport::rheology::{
    CarreauYasuda, ThermalCarreauYasuda, ViscosityShiftMode,
};
use fs_material::visco::temperature::{TemperatureShift, TemperatureShiftModel};

#[path = "newton/mixed.rs"]
mod mixed;
#[path = "newton/mixed_transient.rs"]
mod mixed_transient;
#[path = "newton/sdirk2.rs"]
mod sdirk2;

fn bdf2_couette_moments(mesh: &TriMesh, shear: f64) -> Vec<f64> {
    mesh.edges
        .iter()
        .flat_map(|edge| {
            let a = mesh.verts[edge.verts.0];
            let b = mesh.verts[edge.verts.1];
            // Exact mean and centered first moment of the affine normal trace.
            [
                shear * edge.normal[0] * (a[1] + b[1]) / 2.0,
                shear * edge.normal[0] * (b[1] - a[1]) / 12.0,
            ]
        })
        .collect()
}

#[test]
fn bdf2_material_quadratic_shear_retains_endpoint_source_heat_and_tangents() {
    let mesh = mesh();
    let system = FluxSystem::new(&mesh);
    let amplitude = |t: f64| 1.0 + t + t * t;
    let older = bdf2_couette_moments(&mesh, amplitude(0.1));
    let previous = bdf2_couette_moments(&mesh, amplitude(0.3));
    let saved_histories = (older.clone(), previous.clone());
    let endpoint_time = 0.6;
    let shear = amplitude(endpoint_time);
    let velocity = |p: [f64; 2]| [shear * p[1], 0.0];
    for density in [2.0, 7.0] {
        let reference = transport(density, 0.2 * density);
        let owner = MaterialFluxSystem::new(
            &system,
            &reference,
            MaterialFluxNumerics {
                penalty: 15.0,
                linear_tolerance: 1e-12,
                ..MaterialFluxNumerics::default()
            },
        )
        .unwrap();
        // Uniform shear makes div(tau)=0 and u.grad(u)=0. The physical
        // force is exactly rho*u_t at the new endpoint, not at either history.
        let force = |p: [f64; 2]| [density * (1.0 + 2.0 * endpoint_time) * p[1], 0.0];
        let reads = Rc::new(Cell::new(0));
        let mut queries = 0;
        let mut first_rates = vec![None; mesh.tris.len()];
        let result = owner
            .bdf2_rheology_newton_variable_step(
                &force,
                &velocity,
                &previous,
                Some((&older, 0.2)),
                0.3,
                FluxNewtonControls {
                    residual_tolerance: 1e-10,
                    ..FluxNewtonControls::default()
                },
                |point, retained| {
                    assert_eq!(retained, &reference);
                    first_rates[point.cell].get_or_insert(point.shear_rate_per_s);
                    queries += 1;
                    let rate = point.shear_rate_per_s;
                    Ok(MaterialFluxViscosityTangent {
                        material: RheologyResponse {
                            transport: transport(density, density * 0.2 * (1.0 + rate * rate)),
                            point,
                            temperature_k: 315.0,
                            reads: Rc::clone(&reads),
                        },
                        derivative_pa_s2: density * 0.4 * rate,
                    })
                },
            )
            .unwrap();
        // The material is first sampled at the actual previous velocity,
        // never at BDF2's extrapolated inertial history.
        for rate in first_rates {
            near(rate.unwrap(), amplitude(0.3), 1e-10);
        }
        assert_eq!(reads.get(), queries);
        let step = &result.step;
        let (error, divergence) = system.velocity_error(&step.discrete().x, &velocity);
        assert!(error < 2e-9 && divergence < 2e-9);
        let mu = density * 0.2 * (1.0 + shear * shear);
        let heating = mu * shear * shear;
        near(
            step.kinetic_energy_j_per_m(),
            density * shear * shear / 6.0,
            1e-8,
        );
        near(
            step.viscous_dissipation().total_power_w_per_m(),
            heating,
            1e-6,
        );
        for (cell, sample) in step.samples().iter().enumerate() {
            assert_eq!(sample.point(), sample.material().point);
            assert_eq!(sample.material().temperature_k, 315.0);
            near(sample.point().shear_rate_per_s, shear, 2e-8);
            near(sample.transport().dynamic_viscosity_pa_s(), mu, 1e-7);
            near(
                result.viscosity_derivatives_pa_s2()[cell],
                density * 0.4 * shear,
                1e-7,
            );
            near(
                step.viscous_dissipation().cell_heating_w_m3()[cell],
                heating,
                1e-6,
            );
            near(step.pressure_pa()[cell], 0.0, density * 1e-7);
        }
    }
    assert_eq!((older, previous), saved_histories);
}

#[test]
fn bdf2_material_invalid_history_precedes_callbacks_and_source_errors_stay_terminal() {
    let mesh = mesh();
    let system = FluxSystem::new(&mesh);
    let reference = transport(2.0, 0.4);
    let owner =
        MaterialFluxSystem::new(&system, &reference, MaterialFluxNumerics::default()).unwrap();
    let previous = vec![0.0; system.n_u];
    let older = previous.clone();
    let controls = FluxNewtonControls::default();
    for (history, previous_dt) in [(&older[..], 0.1), (&older[..system.n_u - 1], 0.3)] {
        let result = owner.bdf2_rheology_newton_variable_step(
            &|_| panic!("invalid BDF2 history queried force"),
            &|_| panic!("invalid BDF2 history queried boundary"),
            &previous,
            Some((history, previous_dt)),
            0.3,
            controls,
            |_, _| -> Result<MaterialFluxViscosityTangent<FluidTransportState>, MaterialError> {
                panic!("invalid BDF2 history queried material")
            },
        );
        assert!(matches!(
            result,
            Err(MaterialFluxError::Newton(FluxNewtonError::History(_)))
        ));
    }

    let refusal = MaterialError::State {
        what: "source refused the BDF2 trial shear".into(),
    };
    let mut refused = false;
    let mut queries = 0;
    let force = |p: [f64; 2]| [4.0 * (p[1] - 0.5), -4.0 * (p[0] - 0.5)];
    let result = owner.bdf2_rheology_newton_variable_step(
        &force,
        &no_slip,
        &previous,
        Some((&older, 0.1)),
        0.15,
        controls,
        |point, _| {
            assert!(!refused, "a refused trial must end material queries");
            queries += 1;
            if point.shear_rate_per_s > 0.0 {
                refused = true;
                return Err(refusal.clone());
            }
            Ok(MaterialFluxViscosityTangent {
                material: reference,
                derivative_pa_s2: 0.0,
            })
        },
    );
    assert!(matches!(result, Err(MaterialFluxError::Material(error)) if error == refusal));
    assert!(
        refused && queries > mesh.tris.len(),
        "exercise a genuine trial refusal"
    );
    assert!(previous.iter().chain(&older).all(|value| *value == 0.0));
}

#[test]
fn steady_material_couette_retains_physical_viscosity_heat_and_tangents() {
    let mesh = mesh();
    let system = FluxSystem::new(&mesh);
    let initial = vec![0.0; system.n_u];
    let velocity = |x: [f64; 2]| [2.0 * x[1], 0.0];
    for density in [2.0, 7.0] {
        let reference = transport(density, 0.2 * density);
        let owner = MaterialFluxSystem::new(
            &system,
            &reference,
            MaterialFluxNumerics {
                linear_tolerance: 1e-12,
                ..MaterialFluxNumerics::default()
            },
        )
        .unwrap();
        let result = owner
            .steady_rheology_newton(
                &|_| [0.0; 2],
                &velocity,
                &initial,
                false,
                FluxNewtonControls {
                    residual_tolerance: 1e-10,
                    ..FluxNewtonControls::default()
                },
                |point, retained| {
                    assert_eq!(retained, &reference);
                    let rate = point.shear_rate_per_s;
                    Ok(MaterialFluxViscosityTangent {
                        material: transport(density, density * 0.2 * (1.0 + rate * rate)),
                        derivative_pa_s2: density * 0.4 * rate,
                    })
                },
            )
            .unwrap();
        let step = &result.step;
        let (error, divergence) = system.velocity_error(&step.discrete().x, &velocity);
        assert!(error < 1e-9 && divergence < 1e-9);
        near(step.kinetic_energy_j_per_m(), density * 2.0 / 3.0, 1e-8);
        for (cell, sample) in step.samples().iter().enumerate() {
            near(sample.point().shear_rate_per_s, 2.0, 1e-8);
            near(sample.transport().dynamic_viscosity_pa_s(), density, 1e-8);
            near(
                result.viscosity_derivatives_pa_s2()[cell],
                0.8 * density,
                1e-8,
            );
            near(
                step.viscous_dissipation().cell_heating_w_m3()[cell],
                4.0 * density,
                1e-7,
            );
        }
        near(
            step.viscous_dissipation().total_power_w_per_m(),
            4.0 * density,
            1e-7,
        );
    }
    assert!(initial.iter().all(|value| *value == 0.0));
}

#[test]
fn source_newton_solves_stiff_viscosity_and_retains_exact_endpoint_heat_and_tangents() {
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 4, 4));
    let system = FluxSystem::new(&mesh);
    let previous = vec![0.0; system.n_u];
    let acceleration = |p: [f64; 2]| [100.0 * (p[1] - 0.5), -100.0 * (p[0] - 0.5)];
    let law = |_: usize, _: [f64; 2], rate: f64| {
        let nu = 0.1 * (1.0 + rate * rate).powi(2);
        FluxViscosityTangent {
            nu,
            derivative: 4.0 * rate * nu / (1.0 + rate * rate),
        }
    };
    let controls = FluxNewtonControls {
        maximum_iterations: 30,
        maximum_backtracks: 16,
        residual_tolerance: 1e-9,
    };
    let params = fs_flux::FluxParams {
        tol: 1e-11,
        sigma: 15.0,
        ..fs_flux::FluxParams::default()
    };
    let native = system
        .bdf1_rheology_newton_step(
            params,
            &law,
            &acceleration,
            &no_slip,
            &previous,
            0.3,
            controls,
        )
        .unwrap();
    for density in [2.0, 7.0] {
        let reference = transport(density, density * 0.1);
        let owner = MaterialFluxSystem::new(
            &system,
            &reference,
            MaterialFluxNumerics {
                penalty: params.sigma,
                linear_tolerance: params.tol,
                max_linear_iterations: params.max_iters,
            },
        )
        .unwrap();
        let force = |p| {
            let a = acceleration(p);
            [density * a[0], density * a[1]]
        };
        let reads = Rc::new(Cell::new(0));
        let mut queries = 0;
        let result = owner
            .bdf1_rheology_newton_step(
                &force,
                &no_slip,
                &previous,
                0.3,
                controls,
                |point, retained| {
                    assert_eq!(retained, &reference);
                    queries += 1;
                    let tangent = law(point.cell, point.position_m, point.shear_rate_per_s);
                    Ok(MaterialFluxViscosityTangent {
                        material: RheologyResponse {
                            transport: transport(density, density * tangent.nu),
                            point,
                            temperature_k: 310.0,
                            reads: Rc::clone(&reads),
                        },
                        derivative_pa_s2: density * tangent.derivative,
                    })
                },
            )
            .unwrap();
        let step = &result.step;
        assert!(step.iterations < 20 && result.backtracks > 0);
        assert!(step.nonlinear_residual <= controls.residual_tolerance);
        assert_eq!(
            reads.get(),
            queries,
            "one transport read per original owner"
        );
        for (actual, expected) in step.discrete().x.iter().zip(&native.solution.x) {
            near(*actual, *expected, 2e-9);
        }
        near(
            step.kinetic_energy_j_per_m(),
            density * system.kinetic_energy(&native.solution.x),
            1e-9,
        );
        let rates = system.shear_rates(&step.discrete().x);
        let native_pressure = system.pressure(&native.solution.x);
        let mut power = 0.0;
        for (cell, sample) in step.samples().iter().enumerate() {
            let rate = rates[cell];
            let tangent = law(cell, mesh.centroids[cell], rate);
            assert_eq!(sample.point(), sample.material().point);
            assert_eq!(sample.point().shear_rate_per_s, rate);
            assert_eq!(sample.material().temperature_k, 310.0);
            near(
                result.viscosity_derivatives_pa_s2()[cell],
                density * tangent.derivative,
                0.0,
            );
            let heat = density * tangent.nu * rate * rate;
            near(
                step.viscous_dissipation().cell_heating_w_m3()[cell],
                heat,
                1e-9,
            );
            near(
                step.pressure_pa()[cell],
                density * native_pressure[cell],
                2e-8,
            );
            power += heat * mesh.areas[cell];
        }
        near(
            step.viscous_dissipation().total_power_w_per_m(),
            power,
            1e-9,
        );
        // This is the established stiff physical fixture, not an artificially
        // restricted domain: the same Picard budget cannot solve its equations.
        assert!(
            owner
                .bdf1_rheology_step(&force, &no_slip, &previous, 0.3, 30, 1e-9, |point, _| Ok(
                    transport(
                        density,
                        density * law(point.cell, point.position_m, point.shear_rate_per_s,).nu
                    )
                ),)
                .is_err()
        );
    }
    assert!(previous.iter().all(|value| *value == 0.0));
}

#[test]
fn thermal_carreau_sources_and_derivatives_are_retained_at_the_newton_endpoint() {
    let mesh = mesh();
    let system = FluxSystem::new(&mesh);
    let reference = transport(2.0, 1.6);
    let owner = MaterialFluxSystem::new(
        &system,
        &reference,
        MaterialFluxNumerics {
            linear_tolerance: 1e-12,
            ..MaterialFluxNumerics::default()
        },
    )
    .unwrap();
    let law = ThermalCarreauYasuda::new(
        CarreauYasuda::new(1.6, 0.2, 3.0, 0.4, 2.0, [0.0, 1e4]).unwrap(),
        TemperatureShift::new(
            TemperatureShiftModel::Arrhenius {
                activation_energy_j_mol: 12_000.0,
            },
            300.0,
            [280.0, 340.0],
        )
        .unwrap(),
        ViscosityShiftMode::VerticalAndRate,
    );
    let previous = vec![0.0; system.n_u];
    let force = |p: [f64; 2]| [4.0 * (p[1] - 0.5), -4.0 * (p[0] - 0.5)];
    let result = owner
        .bdf1_rheology_newton_step(
            &force,
            &no_slip,
            &previous,
            0.1,
            FluxNewtonControls {
                residual_tolerance: 1e-10,
                ..FluxNewtonControls::default()
            },
            |point, retained| {
                let state = law.evaluate(310.0, point.shear_rate_per_s)?;
                Ok(MaterialFluxViscosityTangent {
                    material: state.bind(*retained)?,
                    derivative_pa_s2: state.shear_derivative_pa_s2().unwrap(),
                })
            },
        )
        .unwrap();
    assert!(result.step.kinetic_energy_j_per_m() > 1e-5);
    let picard = owner
        .bdf1_rheology_step(
            &force,
            &no_slip,
            &previous,
            0.1,
            100,
            1e-10,
            |point, retained| law.evaluate(310.0, point.shear_rate_per_s)?.bind(*retained),
        )
        .unwrap();
    for (actual, expected) in result.step.discrete().x.iter().zip(&picard.discrete().x) {
        near(*actual, *expected, 2e-9);
    }
    for (cell, sample) in result.step.samples().iter().enumerate() {
        let state = sample.material().rheology();
        assert_eq!(sample.material().source(), &reference);
        assert_eq!(state.shear_rate_per_s(), sample.point().shear_rate_per_s);
        assert_eq!(state.temperature().temperature_k(), 310.0);
        assert_eq!(
            Some(result.viscosity_derivatives_pa_s2()[cell]),
            state.shear_derivative_pa_s2()
        );
    }
}

#[test]
fn source_and_tangent_refusals_stop_newton_without_altering_history() {
    let mesh = mesh();
    let system = FluxSystem::new(&mesh);
    let reference = transport(2.0, 0.4);
    let owner =
        MaterialFluxSystem::new(&system, &reference, MaterialFluxNumerics::default()).unwrap();
    let previous = vec![0.0; system.n_u];
    let force = |p: [f64; 2]| [4.0 * (p[1] - 0.5), -4.0 * (p[0] - 0.5)];
    let controls = FluxNewtonControls::default();
    let refusal = MaterialError::State {
        what: "source is unavailable at trial temperature".into(),
    };
    let mut queries = 0;
    let mut refused = false;
    let result =
        owner.bdf1_rheology_newton_step(&force, &no_slip, &previous, 0.1, controls, |point, _| {
            assert!(!refused, "no more source queries after a refused trial");
            queries += 1;
            if point.shear_rate_per_s > 0.0 {
                refused = true;
                return Err(refusal.clone());
            }
            Ok(MaterialFluxViscosityTangent {
                material: reference,
                derivative_pa_s2: 0.0,
            })
        });
    assert!(matches!(result, Err(MaterialFluxError::Material(error)) if error == refusal));
    assert!(
        refused && queries > mesh.tris.len(),
        "exercise actual line-search source failure"
    );
    for (density, derivative) in [(3.0, 0.0), (2.0, f64::NAN), (2.0, 1.0)] {
        let result =
            owner.bdf1_rheology_newton_step(&force, &no_slip, &previous, 0.1, controls, |_, _| {
                Ok(MaterialFluxViscosityTangent {
                    material: transport(density, 0.4),
                    derivative_pa_s2: derivative,
                })
            });
        assert!(result.is_err());
    }
    let result = owner.bdf1_rheology_newton_step(
        &force,
        &no_slip,
        &previous,
        0.0,
        controls,
        |_, _| -> Result<MaterialFluxViscosityTangent<FluidTransportState>, MaterialError> {
            panic!("invalid timestep must refuse before any constitutive query")
        },
    );
    assert!(matches!(
        result,
        Err(MaterialFluxError::Newton(FluxNewtonError::Invalid(_)))
    ));
    assert!(previous.iter().all(|value| *value == 0.0));
}

#[test]
fn newton_endpoint_tolerance_is_independent_of_the_tight_linear_correction() {
    let mesh = mesh();
    let system = FluxSystem::new(&mesh);
    let params = fs_flux::FluxParams {
        tol: 1e-12,
        sigma: 15.0,
        ..fs_flux::FluxParams::default()
    };
    let previous = vec![0.0; system.n_u];
    let initial = vec![0.0; system.n];
    let timestep = 0.1;
    let acceleration = |p: [f64; 2]| [p[1] - 0.5, -(p[0] - 0.5)];
    let law = |_: usize, _: [f64; 2], rate: f64| FluxViscosityTangent {
        nu: 0.2 * (1.0 + rate * rate),
        derivative: 0.4 * rate,
    };
    // Independently solve the original first Newton correction at the tight
    // linear tolerance. Its actual nonlinear endpoint has a larger defect.
    let first = system
        .linearize_rheology(
            params,
            &law,
            &acceleration,
            &no_slip,
            &initial,
            true,
            1.0 / timestep,
            Some(&previous),
        )
        .unwrap();
    let rhs: Vec<_> = first.residual.iter().map(|r| -r).collect();
    let correction = system.solve_linear(&first.jacobian, &rhs, params);
    assert!(correction.converged(params.tol));
    let endpoint = system
        .linearize_rheology(
            params,
            &law,
            &acceleration,
            &no_slip,
            &correction.x,
            true,
            1.0 / timestep,
            Some(&previous),
        )
        .unwrap();
    let defect = endpoint
        .relative_residual
        .max(endpoint.scaled_relative_residual);
    assert!(defect > 100.0 * params.tol);
    let controls = FluxNewtonControls {
        maximum_iterations: 1,
        maximum_backtracks: 0,
        residual_tolerance: 2.0 * defect,
    };
    assert!(
        controls.residual_tolerance < first.relative_residual.max(first.scaled_relative_residual),
        "the initial state must require a real Newton correction"
    );
    let reference = transport(2.0, 0.4);
    let owner = MaterialFluxSystem::new(
        &system,
        &reference,
        MaterialFluxNumerics {
            penalty: params.sigma,
            linear_tolerance: params.tol,
            max_linear_iterations: params.max_iters,
        },
    )
    .unwrap();
    let force = |point| {
        let a = acceleration(point);
        [2.0 * a[0], 2.0 * a[1]]
    };
    let result = owner
        .bdf1_rheology_newton_step(
            &force,
            &no_slip,
            &previous,
            timestep,
            controls,
            |point, _| {
                let value = law(point.cell, point.position_m, point.shear_rate_per_s);
                Ok(MaterialFluxViscosityTangent {
                    material: transport(2.0, 2.0 * value.nu),
                    derivative_pa_s2: 2.0 * value.derivative,
                })
            },
        )
        .unwrap();
    assert_eq!(result.step.iterations, 1);
    assert_eq!(result.backtracks, 0);
    let accepted = result.step.discrete();
    assert!(accepted.converged(controls.residual_tolerance));
    assert!(!accepted.converged(params.tol));
    for (actual, expected) in accepted.x.iter().zip(&correction.x) {
        near(*actual, *expected, 2e-13);
    }
    // The accepted state still solves that first linear correction to the
    // original tight tolerance; nonlinear acceptance did not relax the solve.
    let mut product = vec![0.0; system.n];
    first.jacobian.spmv(&accepted.x, &mut product);
    let linear_defect = product
        .iter()
        .zip(&rhs)
        .fold(0.0_f64, |norm, (value, load)| norm.hypot(value - load));
    let load_norm = rhs.iter().fold(0.0_f64, |norm, value| norm.hypot(*value));
    assert!(linear_defect / load_norm <= params.tol);
    // Preserve both measured nonlinear residuals instead of disguising them
    // as linear defects or reducing them to pass physical output conversion.
    near(accepted.rel_residual, endpoint.relative_residual, 2e-13);
    near(
        accepted.scaled_rel_residual,
        endpoint.scaled_relative_residual,
        2e-13,
    );
    assert_eq!(
        result.step.nonlinear_residual,
        accepted.rel_residual.max(accepted.scaled_rel_residual)
    );
    assert!(previous.iter().all(|value| *value == 0.0));
}
