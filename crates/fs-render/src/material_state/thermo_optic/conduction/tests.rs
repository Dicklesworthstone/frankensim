use super::*;
use crate::charts::Ray;
use crate::dielectric::{BeerLambertAbsorption, DielectricGlass, GlassProvenance, SellmeierIor};
use crate::gradient_index::GradientIndexOptions;
use crate::gradient_index::bundle::{BundleOptions, ProjectionOptions, propagate_bundle};
use crate::material_state::thermo_optic::MaterialThermoOpticField;
use fs_blake3::ContentHash;
use fs_conduction::ConductivityTable;
use fs_conduction::layered::{LayeredWall, PlanarThermalLayer};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_material::thermo_optic::{SchottThermoOpticLaw, ThermoOpticCoefficients};

fn with_cx<R>(f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate = CancelGate::new();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 71,
                kernel_id: 1,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        f(&cx)
    })
}

fn layer(name: &str, thickness: f64, k: f64) -> PlanarThermalLayer {
    PlanarThermalLayer::new(name, thickness, ConductivityTable::declared(k).unwrap()).unwrap()
}
fn slab(source: f64) -> GeneratedLayeredWallSolution {
    let wall = LayeredWall::new(vec![layer("heated", 2.0, 4.0)], vec![], [0.0; 2]).unwrap();
    with_cx(|cx| wall.solve_with_generation(cx, [300.0; 2], 1.0, &[source])).unwrap()
}
fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:.17e} != {expected:.17e}"
    );
}

#[test]
fn original_conduction_temperature_and_derivatives_rotate_into_world_coordinates() {
    let solution = slab(8.0);
    let field = PlanarConductionTemperatureField::try_new(
        &solution,
        0,
        Point3::new(2.0, -3.0, 4.0),
        Vec3::new(0.6, 0.8, 0.0),
    )
    .unwrap();
    // Offset0.5m along the heat-flow axis plus0.7m in the plane.
    let p = Point3::new(1.74, -2.18, 4.0);
    let jet = field.temperature_jet(p).unwrap();
    // Independent solved profile T(x)=300+2x-x², T'=2-2x, T''=-2.
    near(jet.sample.temperature_k, 300.75, 1.0e-13);
    near(jet.sample.gradient_k_per_m.x, 0.6, 1.0e-14);
    near(jet.sample.gradient_k_per_m.y, 0.8, 1.0e-14);
    assert_eq!(jet.sample.gradient_k_per_m.z, 0.0);
    for (actual, expected) in jet.diagonal_k_per_m2.into_iter().zip([-0.72, -1.28, 0.0]) {
        near(actual, expected, 1.0e-14);
    }
    for (actual, expected) in jet.off_diagonal_k_per_m2.into_iter().zip([-0.96, 0.0, 0.0]) {
        near(actual, expected, 1.0e-14);
    }
    assert_eq!(field.temperature_k(p).unwrap(), jet.sample.temperature_k);
    assert_eq!(field.temperature_gradient(p).unwrap(), jet.sample);
    assert!(core::ptr::eq(field.solution(), &solution));
    assert!(!field.is_spatially_uniform());
}

#[test]
fn conduction_layer_domains_preserve_contact_jumps_without_cross_layer_smoothing() {
    let wall = LayeredWall::new(
        vec![layer("left", 1.0, 1.0), layer("right", 1.0, 2.0)],
        vec![0.5],
        [0.5, 1.0],
    )
    .unwrap();
    let solution =
        with_cx(|cx| wall.solve_with_generation(cx, [300.0; 2], 7.0, &[2.0, 0.0])).unwrap();
    let axis = Vec3::new(1.0, 0.0, 0.0);
    let left =
        PlanarConductionTemperatureField::try_new(&solution, 0, Point3::new(0.0, 0.0, 0.0), axis)
            .unwrap();
    let right =
        PlanarConductionTemperatureField::try_new(&solution, 1, Point3::new(1.0, 0.0, 0.0), axis)
            .unwrap();
    let contact = Point3::new(1.0, 0.0, 0.0);
    let a = left.temperature_gradient(contact).unwrap();
    let b = right.temperature_gradient(contact).unwrap();
    near(a.temperature_k - b.temperature_k, 2.0 / 7.0, 1.0e-13);
    near(a.gradient_k_per_m.x, -4.0 / 7.0, 1.0e-14);
    near(b.gradient_k_per_m.x, -2.0 / 7.0, 1.0e-14);
    assert_eq!(
        left.temperature_k(Point3::new(1.01, 0.0, 0.0)),
        Err(GradientIndexError::OutsideFieldDomain)
    );
    assert_eq!(
        right.temperature_k(Point3::new(0.99, 0.0, 0.0)),
        Err(GradientIndexError::OutsideFieldDomain)
    );
}

#[test]
fn solved_conduction_curvature_drives_native_infrared_thermal_lens_focusing() {
    let solution = slab(8.0);
    let temperature = PlanarConductionTemperatureField::try_new(
        &solution,
        0,
        Point3::new(-1.0, 0.0, 0.0),
        Vec3::new(1.0, 0.0, 0.0),
    )
    .unwrap();
    let reference = DielectricGlass::new(
        SellmeierIor::try_new_in_band([1.25, 0.0, 0.0], [0.0; 3], (1000.0, 2000.0)).unwrap(),
        BeerLambertAbsorption::try_constant(0.2).unwrap(),
        GlassProvenance::Custom,
    );
    // Synthetic reference law: A=5/12 and D0=.024/K give dn/dT=.01/K.
    let law = SchottThermoOpticLaw::try_new(
        ThermoOpticCoefficients {
            d: [0.024, 0.0, 0.0],
            e: [0.0; 2],
            resonance_um: 0.0,
        },
        301.0,
        [299.0, 302.0],
        [1000.0, 2000.0],
        ContentHash([30; 32]),
    )
    .unwrap();
    let field =
        MaterialThermoOpticField::try_new(&reference, &law, &temperature, ContentHash([31; 32]))
            .unwrap();
    let prepared = field.prepare(1550.0).unwrap();
    let length = 0.8;
    let path = propagate_bundle(
        &prepared,
        Ray {
            origin: Point3::new(0.0, 0.0, 0.0),
            dir: Vec3::new(0.0, 0.0, 1.0),
        },
        Vec3::new(1.0, 0.0, 0.0),
        1550.0,
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
    // Original solved T=301-x² gives n=1.5-.01x². On the axial ray,
    // X''=(-.02/1.5)X and Y''=0, so dA/dOmega=L*sin(omega L)/omega.
    let frequency = (0.02_f64 / 1.5).sqrt();
    let expected_area = length * (frequency * length).sin() / frequency;
    let projection = path
        .project_at_endpoint(Vec3::new(0.0, 0.0, 1.0), ProjectionOptions::default())
        .unwrap();
    near(projection.area_per_solid_angle_m2(), expected_area, 1.0e-10);
    assert!((expected_area - length * length).abs() > 1.0e-4);
    near(
        path.propagation().optical_path_length_m,
        1.5 * length,
        2.0e-12,
    );
    near(
        path.propagation().group_optical_path_length_m.unwrap(),
        1.5 * length,
        2.0e-12,
    );
    near(
        path.propagation().power_transmittance(),
        (-0.2 * length).exp(),
        2.0e-12,
    );
    assert_eq!(
        prepared
            .temperature_k(Point3::new(0.0, 0.0, length))
            .unwrap(),
        301.0
    );
}

#[test]
fn invalid_placement_and_even_constant_temperature_keep_the_original_layer_domain() {
    let solution = slab(0.0);
    let origin = Point3::new(0.0, 0.0, 0.0);
    for axis in [
        Vec3::new(0.0, 0.0, 0.0),
        Vec3::new(2.0, 0.0, 0.0),
        Vec3::new(f64::NAN, 0.0, 0.0),
    ] {
        assert!(PlanarConductionTemperatureField::try_new(&solution, 0, origin, axis).is_err());
    }
    let axis = Vec3::new(1.0, 0.0, 0.0);
    assert!(PlanarConductionTemperatureField::try_new(&solution, 1, origin, axis).is_err());
    let field = PlanarConductionTemperatureField::try_new(&solution, 0, origin, axis).unwrap();
    assert!(!field.is_spatially_uniform());
    assert_eq!(field.temperature_k(origin), Ok(300.0));
    assert_eq!(field.temperature_k(Point3::new(2.0, 0.0, 0.0)), Ok(300.0));
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
}
