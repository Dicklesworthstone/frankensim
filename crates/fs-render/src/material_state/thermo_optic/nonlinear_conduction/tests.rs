use super::*;
use crate::charts::Ray;
use crate::dielectric::{BeerLambertAbsorption, DielectricGlass, GlassProvenance, SellmeierIor};
use crate::gradient_index::bundle::{
    BundleOptions, ProjectionOptions, SecondOrderIndexField, propagate_bundle,
};
use crate::gradient_index::{GradientIndexField, GradientIndexOptions};
use crate::material_state::thermo_optic::MaterialThermoOpticField;
use fs_blake3::ContentHash;
use fs_conduction::ConductivityTable;
use fs_conduction::layered::{
    LayeredWall, NonlinearGenerationOptions, NonlinearWallOptions, PlanarThermalLayer,
};
use fs_exec::{Budget, CancelGate, ExecMode, StreamKey};
use fs_material::thermo_optic::{SchottThermoOpticLaw, ThermoOpticCoefficients};

fn with_cx<R>(f: impl FnOnce(&CancelGate, &Cx<'_>) -> R) -> R {
    let gate = CancelGate::new();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 73,
                kernel_id: 1,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        f(&gate, &cx)
    })
}
fn options(interval: [f64; 2]) -> NonlinearGenerationOptions {
    NonlinearGenerationOptions {
        temperature_interval_k: interval,
        solver: NonlinearWallOptions {
            max_iterations: 128,
            temperature_tolerance_k: 1.0e-10,
            flux_absolute_tolerance_w_per_m2: 1.0e-8,
            flux_relative_tolerance: 1.0e-11,
        },
    }
}
fn wall(knots: Vec<(f64, f64)>) -> LayeredWall {
    LayeredWall::new(
        vec![
            PlanarThermalLayer::new(
                "original source",
                2.0,
                ConductivityTable::declared_curve(knots).unwrap(),
            )
            .unwrap(),
        ],
        vec![],
        [0.0; 2],
    )
    .unwrap()
}
fn slab(cx: &Cx<'_>, source: f64) -> NonlinearGeneratedLayeredWallSolution {
    wall(vec![(280.0, 280.0 / 300.0), (800.0, 800.0 / 300.0)])
        .solve_temperature_dependent_with_generation(
            cx,
            [300.0; 2],
            1.0,
            &[source],
            options([280.0, 800.0]),
        )
        .unwrap()
}
fn reference() -> DielectricGlass {
    DielectricGlass::new(
        SellmeierIor::try_new_in_band([1.1, 0.2, 0.8], [0.01, 0.04, 100.0], (1000.0, 2000.0))
            .unwrap(),
        BeerLambertAbsorption::try_constant(0.2).unwrap(),
        GlassProvenance::Custom,
    )
}
fn law(reference_temperature: f64) -> SchottThermoOpticLaw {
    SchottThermoOpticLaw::try_new(
        ThermoOpticCoefficients {
            d: [0.012, 0.0, 0.0],
            e: [0.003, 0.0],
            resonance_um: 0.0,
        },
        reference_temperature,
        [290.0, 400.0],
        [1000.0, 2000.0],
        ContentHash([40; 32]),
    )
    .unwrap()
}
fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:.17e} != {expected:.17e}"
    );
}

#[test]
fn nonlinear_temperature_gradient_and_curvature_rotate_from_the_original_source() {
    with_cx(|_, cx| {
        let solution = slab(cx, 80.0);
        let field = PlanarNonlinearConductionTemperatureField::try_new(
            cx,
            &solution,
            0,
            Point3::new(2.0, -3.0, 4.0),
            Vec3::new(0.6, 0.8, 0.0),
        )
        .unwrap();
        let p = Point3::new(1.74, -2.18, 4.0); // Original offset0.5 plus transverse0.7.
        let temperature = 108_000.0_f64.sqrt();
        let gradient = 12_000.0 / temperature;
        let curvature = -24_000.0 / temperature - 12_000.0_f64.powi(2) / temperature.powi(3);
        let jet = field.temperature_jet(p).unwrap();
        near(jet.sample.temperature_k, temperature, 1.0e-12);
        near(jet.sample.gradient_k_per_m.x, 0.6 * gradient, 1.0e-12);
        near(jet.sample.gradient_k_per_m.y, 0.8 * gradient, 1.0e-12);
        assert_eq!(jet.sample.gradient_k_per_m.z, 0.0);
        for (actual, factor) in jet.diagonal_k_per_m2.into_iter().zip([0.36, 0.64, 0.0]) {
            near(actual, factor * curvature, 1.0e-12);
        }
        for (actual, factor) in jet.off_diagonal_k_per_m2.into_iter().zip([0.48, 0.0, 0.0]) {
            near(actual, factor * curvature, 1.0e-12);
        }
        assert_eq!(field.temperature_gradient(p).unwrap(), jet.sample);
        assert_eq!(field.temperature_k(p).unwrap(), jet.sample.temperature_k);
        assert!(core::ptr::eq(field.solution(), &solution));
        assert!(!field.is_spatially_uniform());
    });
}

#[test]
fn nonlinear_conduction_drives_chromatic_native_infrared_bundle_focusing() {
    with_cx(|_, cx| {
        let solution = slab(cx, 80.0);
        let temperature = PlanarNonlinearConductionTemperatureField::try_new(
            cx,
            &solution,
            0,
            Point3::new(-1.0, 0.0, 0.0),
            Vec3::new(1.0, 0.0, 0.0),
        )
        .unwrap();
        let centre_temperature = 114_000.0_f64.sqrt();
        let reference = reference();
        let law = law(centre_temperature);
        let field = MaterialThermoOpticField::try_new(
            &reference,
            &law,
            &temperature,
            ContentHash([41; 32]),
        )
        .unwrap();
        // Independent high-precision Sellmeier phase/group reference values.
        for (wavelength, phase, group) in [
            (1200.0, 1.5171414837510317, 1.5339380276614663),
            (1550.0, 1.5127097597193700, 1.5313771849709665),
            (1800.0, 1.5096746814659043, 1.5319525185732263),
        ] {
            let prepared = field.prepare(wavelength).unwrap();
            let length = 0.8;
            let path = propagate_bundle(
                &prepared,
                Ray {
                    origin: Point3::new(0.0, 0.0, 0.0),
                    dir: Vec3::new(0.0, 0.0, 1.0),
                },
                Vec3::new(1.0, 0.0, 0.0),
                wavelength,
                length,
                BundleOptions {
                    integration: GradientIndexOptions {
                        max_step_m: 0.04,
                        relative_tolerance: 2.0e-12,
                        position_tolerance_m: 2.0e-12,
                        direction_tolerance: 2.0e-12,
                        optical_path_tolerance_m: 2.0e-12,
                        optical_depth_tolerance: 2.0e-12,
                        ..GradientIndexOptions::default()
                    },
                    position_derivative_tolerance_m: 2.0e-12,
                    direction_derivative_tolerance: 2.0e-12,
                },
            )
            .unwrap();
            // From T²=114000-24000x²: T'=0 and T''=-24000/Tc on axis.
            // At the authored reference T: n_T=(n-1/n)/2*(D0+E0/lambda²).
            // The native Jacobi equations give X''=Hxx/n X and Y''=0.
            let lambda_um = wavelength * 1.0e-3;
            let dn_dt = 0.5 * (phase - 1.0 / phase) * (0.012 + 0.003 / (lambda_um * lambda_um));
            let frequency = (dn_dt * 24_000.0 / (centre_temperature * phase)).sqrt();
            let expected_area = length * (frequency * length).sin() / frequency;
            let projection = path
                .project_at_endpoint(Vec3::new(0.0, 0.0, 1.0), ProjectionOptions::default())
                .unwrap();
            near(projection.area_per_solid_angle_m2(), expected_area, 1.0e-10);
            assert!((expected_area - length * length).abs() > 0.01);
            near(
                path.propagation().optical_path_length_m,
                phase * length,
                2.0e-12,
            );
            near(
                path.propagation().group_optical_path_length_m.unwrap(),
                group * length,
                2.0e-12,
            );
            near(
                path.propagation().power_transmittance(),
                (-0.2 * length).exp(),
                2.0e-12,
            );
        }
    });
}

#[test]
fn nonlinear_source_knot_keeps_native_bending_but_refuses_undefined_focusing() {
    with_cx(|_, cx| {
        let wall = wall(vec![(290.0, 1.0), (302.0, 1.0), (304.0, 3.0), (310.0, 3.0)]);
        let solution = wall
            .solve_temperature_dependent_with_generation(
                cx,
                [304.0, 300.0],
                1.0,
                &[0.0],
                options([290.0, 304.0]),
            )
            .unwrap();
        let temperature = PlanarNonlinearConductionTemperatureField::try_new(
            cx,
            &solution,
            0,
            Point3::new(0.0, 0.0, 0.0),
            Vec3::new(1.0, 0.0, 0.0),
        )
        .unwrap();
        let knot = Point3::new(4.0 / 3.0, 0.0, 0.0);
        assert_eq!(temperature.temperature_k(knot), Ok(302.0));
        near(
            temperature
                .temperature_gradient(knot)
                .unwrap()
                .gradient_k_per_m
                .x,
            -3.0,
            1.0e-12,
        );
        assert_eq!(
            temperature.temperature_jet(knot),
            Err(GradientIndexError::InvalidSample)
        );
        let reference = reference();
        let law = law(302.0);
        let field = MaterialThermoOpticField::try_new(
            &reference,
            &law,
            &temperature,
            ContentHash([42; 32]),
        )
        .unwrap();
        let prepared = field.prepare(1550.0).unwrap();
        assert!(prepared.sample(knot, 1550.0).is_ok());
        assert_eq!(
            prepared.sample_second_order(knot, 1550.0),
            Err(GradientIndexError::InvalidSample)
        );
        let zero_flux = wall
            .solve_temperature_dependent_with_generation(
                cx,
                [300.0; 2],
                1.0,
                &[4.0],
                options([290.0, 310.0]),
            )
            .unwrap();
        let field = PlanarNonlinearConductionTemperatureField::try_new(
            cx,
            &zero_flux,
            0,
            Point3::new(0.0, 0.0, 0.0),
            Vec3::new(1.0, 0.0, 0.0),
        )
        .unwrap();
        let peak = field.temperature_jet(Point3::new(1.0, 0.0, 0.0)).unwrap();
        assert_eq!(peak.sample.temperature_k, 302.0);
        assert_eq!(peak.diagonal_k_per_m2, [-4.0, 0.0, 0.0]);
    });
}

#[test]
fn nonlinear_field_retains_layer_support_and_runtime_cancellation() {
    with_cx(|gate, cx| {
        let solution = slab(cx, 0.0);
        let origin = Point3::new(0.0, 0.0, 0.0);
        let axis = Vec3::new(1.0, 0.0, 0.0);
        for invalid in [
            Vec3::new(2.0, 0.0, 0.0),
            Vec3::new(0.0, 0.0, 0.0),
            Vec3::new(f64::NAN, 0.0, 0.0),
        ] {
            assert!(
                PlanarNonlinearConductionTemperatureField::try_new(
                    cx, &solution, 0, origin, invalid
                )
                .is_err()
            );
        }
        let field =
            PlanarNonlinearConductionTemperatureField::try_new(cx, &solution, 0, origin, axis)
                .unwrap();
        assert_eq!(field.temperature_k(origin), Ok(300.0));
        assert_eq!(field.temperature_k(Point3::new(2.0, 0.0, 0.0)), Ok(300.0));
        assert!(!field.is_spatially_uniform());
        for p in [
            Point3::new(-0.01, 0.0, 0.0),
            Point3::new(2.01, 0.0, 0.0),
            Point3::new(1.0, f64::NAN, 0.0),
        ] {
            assert_eq!(
                field.temperature_k(p),
                Err(GradientIndexError::OutsideFieldDomain)
            );
        }
        gate.request();
        assert_eq!(
            field.temperature_k(origin),
            Err(GradientIndexError::Cancelled)
        );
        assert_eq!(
            field.temperature_gradient(origin),
            Err(GradientIndexError::Cancelled)
        );
        assert_eq!(
            field.temperature_jet(origin),
            Err(GradientIndexError::Cancelled)
        );
    });
}
