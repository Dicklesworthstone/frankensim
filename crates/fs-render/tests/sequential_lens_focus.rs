//! Prescribed-object focusing retains the native lens and exact refraction path.
#![cfg(feature = "tracer")]

use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_geom::{Point3, Vec3};
use fs_render::camera::{
    Aperture, CameraError, CameraProjection, LensCoating, LensSample, LensSurface,
    PhysicalCamera, SequentialLens,
};
use fs_render::charts::Ray;
use fs_render::conductor::{
    CauchyCoatingLayer, CoherentDielectricCoating, CoherentPeriodicCoating,
    ConductorDataStatus, ConductorSource,
};
use fs_render::dielectric::{
    BeerLambertAbsorption, CauchyIor, DielectricCoating, DielectricGlass, GlassProvenance,
};

fn with_cx(f: impl FnOnce(&Cx<'_>)) {
    let gate = CancelGate::new_clock_free();
    fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default()).scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 29,
                kernel_id: 14,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        f(&cx);
    });
}

fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:e} != {expected:e}, tolerance {tolerance:e}"
    );
}

fn glass(a: f64, b: f64) -> DielectricGlass {
    DielectricGlass::new(
        CauchyIor::try_new(a, b, 0.0).unwrap(),
        BeerLambertAbsorption::CLEAR,
        GlassProvenance::Custom,
    )
}

fn prescription(
    scale: f64,
    sensor: Option<DielectricGlass>,
    interior: DielectricGlass,
    scene: Option<DielectricGlass>,
) -> SequentialLens {
    SequentialLens::try_new_with_sensor_medium(
        550.0,
        sensor,
        &[
            LensSurface::try_new(
                0.08 * scale, 20.0 / scale, 0.012 * scale, Some(interior),
            )
            .unwrap(),
            LensSurface::try_new(0.09 * scale, -20.0 / scale, 0.012 * scale, scene)
                .unwrap(),
        ],
    )
    .unwrap()
}

fn exact_axis_intercept(cx: &Cx<'_>, lens: &SequentialLens, wavelength_nm: f64) -> f64 {
    let ray = lens
        .trace_ray(
            cx,
            Ray {
                origin: Point3::new(0.0, 0.0, 0.0),
                dir: Vec3::new(1e-7, 0.0, 1.0),
            },
            wavelength_nm,
        )
        .unwrap()
        .unwrap()
        .ray;
    ray.origin.z - ray.origin.x * ray.dir.z / ray.dir.x
}

#[test]
fn spectral_and_immersed_refocusing_reaches_the_original_camera_and_exact_rays() {
    with_cx(|cx| {
        for immersed in [false, true] {
            let original = if immersed {
                prescription(
                    1.0, Some(glass(1.1, 0.003)), glass(1.8, 0.01),
                    Some(glass(1.3, 0.002)),
                )
            } else {
                prescription(1.0, None, glass(1.48, 0.006), None)
            };
            let unchanged = original.clone();
            for wavelength in [450.0, 650.0] {
                let focused = original
                    .refocused_for_object(wavelength, 0.8, [-0.06, 0.02], 1e-11)
                    .unwrap();
                let lens = focused.lens();
                assert_eq!(focused.target_distance_m(), 0.8);
                near(focused.achieved_distance_m(), 0.8, 1e-11);
                assert_eq!(lens.reference_wavelength_nm(), wavelength);
                for (old, new) in original.surfaces().iter().zip(lens.surfaces()) {
                    near(
                        new.vertex_z_m() - old.vertex_z_m(),
                        focused.axial_shift_m(),
                        2e-17,
                    );
                    assert_eq!(new.curvature_per_m(), old.curvature_per_m());
                    assert_eq!(new.medium_after(), old.medium_after());
                }
                near(exact_axis_intercept(cx, lens, wavelength), 0.8, 3e-8);
                let expected = lens.clone();
                let camera = PhysicalCamera::try_look_at(
                    Point3::new(0.0, 0.0, 0.0),
                    Point3::new(0.0, 0.0, 1.0),
                    Vec3::new(0.0, 1.0, 0.0),
                    CameraProjection::try_focal_sensor(0.05, 0.024).unwrap(),
                    1.0,
                    Aperture::try_circular(0.001).unwrap(),
                )
                .unwrap()
                .with_sequential_lens(focused.into_lens())
                .unwrap();
                near(camera.focus_distance_m(), 0.8, 1e-11);
                assert_eq!(camera.sequential_lens(), Some(&expected));
                let ray = camera
                    .generate_sequential_spectral_rays_from_tangent_offsets(
                        cx, 0.0, 0.0, LensSample::CENTER, &[wavelength],
                    )
                    .unwrap()[0]
                    .unwrap();
                assert!(ray.weight > 0.0);
                near(ray.ray.dir.z, 1.0, 1e-14);
            }
            assert_eq!(original, unchanged);
        }
    });
}

#[test]
fn mechanical_travel_selects_both_physical_conjugates_and_prefers_the_nearest() {
    let original = prescription(1.0, None, glass(1.5, 0.0), None);
    let near_focus = original
        .refocused_for_object(550.0, 0.5, [-0.07, 0.5], 1e-11)
        .unwrap();
    let far_focus = original
        .refocused_for_object(550.0, 0.5, [0.1, 0.5], 1e-11)
        .unwrap();
    // Independent thick-lens conjugates for n=1.5, R=50 mm, thickness=10 mm.
    near(near_focus.axial_shift_m(), -0.024802528961023074, 3e-16);
    near(far_focus.axial_shift_m(), 0.35480252896102316, 3e-16);
    assert!(near_focus.axial_shift_m().abs() < far_focus.axial_shift_m().abs());
    with_cx(|cx| {
        for focused in [&near_focus, &far_focus] {
            near(focused.achieved_distance_m(), 0.5, 1e-11);
            near(exact_axis_intercept(cx, focused.lens(), 550.0), 0.5, 3e-8);
        }
    });
    assert_eq!(
        original.refocused_for_object(550.0, 0.5, [0.0, 0.1], 1e-11),
        Err(CameraError::SequentialLensOutsideModel)
    );
}

#[test]
fn refocusing_preserves_aspheres_stops_original_films_and_ghost_proposals() {
    let source = ConductorSource::try_new(
        fs_blake3::ContentHash([17; 32]),
        ConductorDataStatus::Representative,
    )
    .unwrap();
    let coating_law = CauchyIor::try_constant(1.3).unwrap();
    let compact = DielectricCoating::try_new(
        CoherentPeriodicCoating::try_new(
            &[CauchyCoatingLayer::try_new(coating_law, 0.01, 100e-9).unwrap()],
            1,
            source,
        )
        .unwrap(),
    )
    .unwrap();
    let film = LensCoating::try_new(&[
        CoherentDielectricCoating::try_new(coating_law, 120e-9, source)
            .unwrap()
            .into(),
    ])
    .unwrap();
    let interior = glass(1.5, 0.0);
    let original = SequentialLens::try_new(
        550.0,
        &[
            LensSurface::try_new_aspheric(
                0.08, 20.0, -0.5, [50.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                0.012, Some(interior),
            )
            .unwrap()
            .with_coating(compact),
            LensSurface::try_new(0.085, 0.0, 0.008, Some(interior)).unwrap(),
            LensSurface::try_new(0.09, -20.0, 0.012, None).unwrap(),
        ],
    )
    .unwrap()
    .with_spectral_coating(2, film)
    .unwrap()
    .with_two_reflection_ghosts(0.7)
    .unwrap();
    let focused = original
        .refocused_for_object(650.0, 0.8, [-0.06, 0.02], 1e-11)
        .unwrap();
    let lens = focused.lens();
    assert_eq!(lens.ghost_direct_probability(), Some(0.7));
    assert_eq!(lens.wavelength_bounds_nm(), original.wavelength_bounds_nm());
    assert_eq!(lens.sensor_medium(), original.sensor_medium());
    assert_eq!(lens.scene_medium(), original.scene_medium());
    for (index, (old, new)) in original.surfaces().iter().zip(lens.surfaces()).enumerate() {
        near(
            new.vertex_z_m() - old.vertex_z_m(), focused.axial_shift_m(), 2e-17,
        );
        assert_eq!(new.conic_constant(), old.conic_constant());
        assert_eq!(new.even_aspheric_coefficients(), old.even_aspheric_coefficients());
        assert_eq!(new.curvature_per_m(), old.curvature_per_m());
        assert_eq!(new.clear_radius_m(), old.clear_radius_m());
        assert_eq!(new.medium_after(), old.medium_after());
        assert_eq!(new.coating(), old.coating());
        assert_eq!(lens.spectral_coating(index), original.spectral_coating(index));
    }
    assert_eq!(lens.surfaces()[1].curvature_per_m(), 0.0);
    assert_eq!(original.reference_wavelength_nm(), 550.0);
}

#[test]
fn explicit_domain_travel_and_residual_limits_refuse_without_fabricating_a_focus() {
    let lens = prescription(1.0, None, glass(1.5, 0.0), None);
    let focus = lens.reference_focus().unwrap().1;
    let unchanged = lens.refocused_for_object(550.0, focus, [0.0, 0.0], 0.0).unwrap();
    assert_eq!(unchanged.axial_shift_m(), 0.0);
    assert_eq!(unchanged.into_lens(), lens);
    // A currently virtual conjugate is not a barrier to repositioning a
    // positive-power prescription into a real, physically accessible focus.
    let virtual_lens = SequentialLens::try_new(
        550.0,
        &[
            LensSurface::try_new(0.03, 20.0, 0.012, Some(glass(1.5, 0.0))).unwrap(),
            LensSurface::try_new(0.04, -20.0, 0.012, None).unwrap(),
        ],
    )
    .unwrap();
    assert_eq!(
        virtual_lens.reference_focus(),
        Err(CameraError::SequentialLensOutsideModel)
    );
    let real_focus = virtual_lens
        .refocused_for_object(550.0, 0.8, [0.0, 0.1], 1e-11)
        .unwrap();
    near(real_focus.achieved_distance_m(), 0.8, 1e-11);
    for distance in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert_eq!(
            lens.refocused_for_object(550.0, distance, [-0.1, 0.1], 1e-10),
            Err(CameraError::InvalidFocusDistance)
        );
    }
    for (wavelength, distance, range, tolerance) in [
        (550.0, 0.8, [0.0, 0.0], 1e-10),
        (550.0, 0.05, [-0.08, 0.1], 1e-10),
        (550.0, 0.8, [0.1, -0.1], 1e-10),
        (550.0, 0.8, [f64::NEG_INFINITY, 0.1], 1e-10),
        (550.0, 0.8, [-0.1, 0.1], -1.0),
        (550.0, 0.8, [-0.1, 0.1], f64::NAN),
        (379.0, 0.8, [-0.1, 0.1], 1e-10),
        (781.0, 0.8, [-0.1, 0.1], 1e-10),
    ] {
        assert_eq!(
            lens.refocused_for_object(wavelength, distance, range, tolerance),
            Err(CameraError::SequentialLensOutsideModel)
        );
    }
}

#[test]
fn physical_length_scaling_preserves_the_conjugate_solution() {
    let original = prescription(1.0, None, glass(1.5, 0.0), None);
    let shift = original
        .refocused_for_object(550.0, 0.8, [-0.06, 0.02], 1e-11)
        .unwrap()
        .axial_shift_m();
    with_cx(|cx| {
        for scale in [1e-3, 1.0, 1e3] {
            let lens = prescription(scale, None, glass(1.5, 0.0), None);
            let focused = lens
                .refocused_for_object(
                    550.0, 0.8 * scale, [-0.06 * scale, 0.02 * scale], 1e-11 * scale,
                )
                .unwrap();
            near(focused.axial_shift_m(), shift * scale, 3e-15 * scale);
            near(focused.achieved_distance_m(), 0.8 * scale, 1e-11 * scale);
            near(
                exact_axis_intercept(cx, focused.lens(), 550.0),
                0.8 * scale,
                3e-8 * scale,
            );
        }
    });
}
