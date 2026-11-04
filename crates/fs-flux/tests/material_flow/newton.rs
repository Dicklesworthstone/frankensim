//! G1/G3/G4: original constitutive tangents reach the actual nonlinear solver.

use super::*;
use fs_flux::MaterialFluxViscosityTangent;
use fs_flux::ns::linearization::FluxViscosityTangent;
use fs_flux::ns::newton::{FluxNewtonControls, FluxNewtonError};
use fs_material::fluid_transport::rheology::{
    CarreauYasuda, ThermalCarreauYasuda, ViscosityShiftMode,
};
use fs_material::visco::temperature::{TemperatureShift, TemperatureShiftModel};

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
