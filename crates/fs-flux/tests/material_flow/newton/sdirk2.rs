use super::*;
use fs_flux::ns::newton::{FluxSdirk2Error, FluxSdirk2Stage};

#[test]
fn sdirk2_material_stages_retain_distinct_sources_and_physical_work() {
    let mesh = mesh();
    let system = FluxSystem::new(&mesh);
    let time = 0.2;
    let timestep = 0.4;
    let slope = 0.7;
    let amplitude = |t: f64| 1.0 + slope * t;
    let previous = bdf2_couette_moments(&mesh, amplitude(time));
    let original = previous.clone();
    let index = |stage| match stage {
        FluxSdirk2Stage::First => 0,
        FluxSdirk2Stage::Endpoint => 1,
    };
    let gradient_pa_m = |stage| [3.0, 7.0][index(stage)];
    let boundary = |stage: FluxSdirk2Stage, p: [f64; 2]| {
        [
            amplitude(time + stage.time_fraction() * timestep) * p[1],
            0.0,
        ]
    };
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
        let force = |stage, p: [f64; 2]| [density * slope * p[1] + gradient_pa_m(stage), 0.0];
        let mut first_rates = [vec![None; mesh.tris.len()], vec![None; mesh.tris.len()]];
        let reads = Rc::new(Cell::new(0));
        let mut queries = 0;
        let result = owner
            .sdirk2_rheology_newton_step(
                &force,
                &boundary,
                &previous,
                timestep,
                true,
                FluxNewtonControls {
                    residual_tolerance: 1e-10,
                    ..FluxNewtonControls::default()
                },
                |stage, point, retained| {
                    assert_eq!(retained, &reference);
                    let i = index(stage);
                    first_rates[i][point.cell].get_or_insert(point.shear_rate_per_s);
                    queries += 1;
                    let rate = point.shear_rate_per_s;
                    let base = density * [0.2, 0.6][i];
                    Ok(MaterialFluxViscosityTangent {
                        material: RheologyResponse {
                            transport: transport(density, base * (1.0 + rate * rate)),
                            point,
                            temperature_k: [300.0, 330.0][i],
                            reads: Rc::clone(&reads),
                        },
                        derivative_pa_s2: 2.0 * base * rate,
                    })
                },
            )
            .unwrap();
        assert_eq!(
            reads.get(),
            queries,
            "one transport read per original stage owner"
        );
        let first_time = time + FluxSdirk2Stage::First.time_fraction() * timestep;
        for cell in 0..mesh.tris.len() {
            near(first_rates[0][cell].unwrap(), amplitude(time), 1e-10);
            near(first_rates[1][cell].unwrap(), amplitude(first_time), 2e-8);
        }
        for (stage, result) in [
            (FluxSdirk2Stage::First, &result.first_stage),
            (FluxSdirk2Stage::Endpoint, &result.endpoint),
        ] {
            let i = index(stage);
            let shear = amplitude(time + stage.time_fraction() * timestep);
            let mu = density * [0.2, 0.6][i] * (1.0 + shear * shear);
            let heat = mu * shear * shear;
            let step = &result.step;
            let (error, divergence) =
                system.velocity_error(&step.discrete().x, &|p| boundary(stage, p));
            assert!(error < 3e-9 && divergence < 3e-9);
            near(
                step.kinetic_energy_j_per_m(),
                density * shear * shear / 6.0,
                2e-8,
            );
            near(step.viscous_dissipation().total_power_w_per_m(), heat, 2e-6);
            for (cell, sample) in step.samples().iter().enumerate() {
                assert_eq!(sample.point(), sample.material().point);
                assert_eq!(sample.material().temperature_k, [300.0, 330.0][i]);
                near(sample.point().shear_rate_per_s, shear, 3e-8);
                near(sample.transport().dynamic_viscosity_pa_s(), mu, 2e-7);
                near(
                    result.viscosity_derivatives_pa_s2()[cell],
                    2.0 * density * [0.2, 0.6][i] * shear,
                    2e-7,
                );
                near(
                    step.viscous_dissipation().cell_heating_w_m3()[cell],
                    heat,
                    2e-6,
                );
                near(
                    step.pressure_pa()[cell],
                    gradient_pa_m(stage) * (mesh.centroids[cell][0] - 0.5),
                    2e-7 * density,
                );
            }
        }
    }
    assert_eq!(previous, original);
}

#[test]
fn sdirk2_material_second_stage_refusal_preserves_error_and_old_history() {
    let mesh = mesh();
    let system = FluxSystem::new(&mesh);
    let reference = transport(2.0, 0.4);
    let owner =
        MaterialFluxSystem::new(&system, &reference, MaterialFluxNumerics::default()).unwrap();
    let previous = vec![0.0; system.n_u];
    let refusal = MaterialError::State {
        what: "endpoint temperature is outside the material source".into(),
    };
    let mut first_queries = 0;
    let mut endpoint_queries = 0;
    let force = |_, p: [f64; 2]| [4.0 * (p[1] - 0.5), -4.0 * (p[0] - 0.5)];
    let result = owner.sdirk2_rheology_newton_step(
        &force,
        &|_, p| no_slip(p),
        &previous,
        0.1,
        true,
        FluxNewtonControls::default(),
        |stage, _, _| {
            if stage == FluxSdirk2Stage::Endpoint {
                endpoint_queries += 1;
                assert_eq!(endpoint_queries, 1, "a source refusal is terminal");
                return Err(refusal.clone());
            }
            first_queries += 1;
            Ok(MaterialFluxViscosityTangent {
                material: reference,
                derivative_pa_s2: 0.0,
            })
        },
    );
    assert!(first_queries >= mesh.tris.len() && endpoint_queries == 1);
    assert!(matches!(result, Err(FluxSdirk2Error {
        stage: Some(FluxSdirk2Stage::Endpoint),
        error: MaterialFluxError::Material(error),
    }) if error == refusal));
    assert!(previous.iter().all(|value| *value == 0.0));
    let invalid = owner.sdirk2_rheology_newton_step(
        &|_, _| panic!("invalid step queried force"),
        &|_, _| panic!("invalid step queried boundary"),
        &previous,
        0.0,
        true,
        FluxNewtonControls::default(),
        |_, _, _| -> Result<MaterialFluxViscosityTangent<FluidTransportState>, MaterialError> {
            panic!("invalid step queried the material source")
        },
    );
    assert!(matches!(
        invalid,
        Err(FluxSdirk2Error {
            stage: None,
            error: MaterialFluxError::Newton(FluxNewtonError::Invalid(_)),
        })
    ));
}
