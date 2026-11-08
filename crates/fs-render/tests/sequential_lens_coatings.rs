//! Original spectral multilayers on real lens paths, including nonvisible bands.
#![cfg(feature = "tracer")]

use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_geom::{Point3, Vec3};
use fs_render::camera::{
    Aperture, CameraError, CameraProjection, LensCoating, LensCoatingLayer, LensSample,
    LensSurface, PhysicalCamera, SequentialLens,
};
use fs_render::charts::Ray;
use fs_render::conductor::{
    CauchyCoatingLayer, CoatingIorSample, CoherentAbsorbingCoating, CoherentCoatingStack,
    CoherentDielectricCoating, CoherentPeriodicCoating, ConductorDataStatus, ConductorSource,
};
use fs_render::dielectric::{
    BeerLambertAbsorption, CauchyIor, DielectricCoating, DielectricCoatingStack, DielectricGlass,
    DielectricIorSample, GlassProvenance, SellmeierIor, TabulatedIor,
};
use fs_render::polarization::Stokes;

fn with_cx(f: impl FnOnce(&Cx<'_>)) {
    let gate = CancelGate::new_clock_free();
    fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default()).scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 41,
                kernel_id: 11,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        f(&cx);
    });
}

fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() < 8e-13,
        "{actual:e} != {expected:e}"
    );
}

fn source(tag: u8) -> ConductorSource {
    ConductorSource::try_new(
        fs_blake3::ContentHash([tag; 32]),
        ConductorDataStatus::Representative,
    )
    .unwrap()
}

fn constant_law(index: f64, band: (f64, f64)) -> SellmeierIor {
    SellmeierIor::try_new_in_band([index * index - 1.0, 0.0, 0.0], [0.0; 3], band).unwrap()
}

fn broad_coating(index: f64, thickness: f64, band: (f64, f64)) -> LensCoating {
    LensCoating::try_new(&[CoherentDielectricCoating::try_new(
        constant_law(index, band),
        thickness,
        source(1),
    )
    .unwrap()
    .into()])
    .unwrap()
}

fn compact(thickness: f64) -> DielectricCoating {
    DielectricCoating::try_new(
        CoherentPeriodicCoating::try_new(
            &[
                CauchyCoatingLayer::try_new(CauchyIor::try_constant(1.3).unwrap(), 0.02, thickness)
                    .unwrap(),
            ],
            1,
            source(2),
        )
        .unwrap(),
    )
    .unwrap()
}

fn material(index: f64, sigma: f64) -> DielectricGlass {
    DielectricGlass::new(
        constant_law(index, (400.0, 1800.0)),
        BeerLambertAbsorption::try_constant(sigma).unwrap(),
        GlassProvenance::Custom,
    )
}

fn biconvex() -> SequentialLens {
    SequentialLens::try_new(
        550.0,
        &[
            LensSurface::try_new(0.08, 20.0, 0.012, Some(material(1.5, 7.0))).unwrap(),
            LensSurface::try_new(0.09, -20.0, 0.012, None).unwrap(),
        ],
    )
    .unwrap()
}

fn camera(lens: SequentialLens) -> PhysicalCamera {
    PhysicalCamera::try_look_at(
        Point3::new(0.0, 0.0, 0.0),
        Point3::new(0.0, 0.0, 1.0),
        Vec3::new(0.0, 1.0, 0.0),
        CameraProjection::try_focal_sensor(0.05, 0.024).unwrap(),
        1.0,
        Aperture::try_circular(0.001).unwrap(),
    )
    .unwrap()
    .with_sequential_lens(lens)
    .unwrap()
}

fn mixed_layers() -> [LensCoatingLayer; 3] {
    let transparent = CoherentDielectricCoating::try_new(
        SellmeierIor::try_new_in_band([0.8, 0.0, 0.0], [0.02, 0.0, 0.0], (500.0, 1700.0)).unwrap(),
        160e-9,
        source(3),
    )
    .unwrap();
    let indices = [1.8, 1.85, 1.78, 1.9, 1.82, 1.91, 1.84, 1.95, 1.86];
    let table = TabulatedIor::try_new_in_band(
        std::array::from_fn(|i| {
            DielectricIorSample::try_new_in_band(500.0 + 150.0 * i as f64, indices[i]).unwrap()
        }),
        150.0,
    )
    .unwrap();
    let tabulated = CoherentDielectricCoating::try_new(table, 90e-9, source(4)).unwrap();
    let absorbing = CoherentAbsorbingCoating::try_new_in_band(
        std::array::from_fn(|i| {
            CoatingIorSample::try_new_in_band(
                600.0 + 112.5 * i as f64,
                2.2 - 0.04 * i as f64,
                0.04 + 0.008 * i as f64,
            )
            .unwrap()
        }),
        120e-9,
        112.5,
        source(5),
    )
    .unwrap();
    [transparent.into(), absorbing.into(), tabulated.into()]
}

#[test]
fn original_sellmeier_table_and_absorbing_films_reach_nonvisible_lens_paths() {
    let layers = mixed_layers();
    let coating = LensCoating::try_new(&layers).unwrap();
    let borrowed: Vec<_> = layers.iter().map(LensCoatingLayer::as_borrowed).collect();
    let native = DielectricCoatingStack::try_mixed_in_band(&borrowed, (600.0, 1500.0)).unwrap();
    assert_eq!(coating.layers(), &layers);
    assert_eq!(coating.wavelength_bounds_nm(), (600.0, 1500.0));
    assert_eq!(
        coating.content_identity(),
        CoherentCoatingStack::try_new(&borrowed)
            .unwrap()
            .content_identity()
    );
    let bulk = DielectricGlass::new(
        SellmeierIor::try_new_in_band([1.1, 0.12, 0.0], [0.01, 0.04, 0.0], (400.0, 1800.0))
            .unwrap(),
        BeerLambertAbsorption::try_constant(0.7).unwrap(),
        GlassProvenance::Custom,
    );
    let bare = SequentialLens::try_new(
        1064.0,
        &[
            LensSurface::try_new(0.02, 0.0, 0.1, Some(bulk)).unwrap(),
            LensSurface::try_new(0.03, 0.0, 0.1, None).unwrap(),
        ],
    )
    .unwrap();
    let lens = bare
        .clone()
        .with_spectral_coating(0, coating.clone())
        .unwrap()
        .with_spectral_coating(1, coating)
        .unwrap();
    assert_eq!(lens.wavelength_bounds_nm(), (600.0, 1500.0));
    let reversed = LensCoating::try_new(&[layers[2], layers[1], layers[0]]).unwrap();
    let reversed = lens.clone().with_spectral_coating(1, reversed).unwrap();
    let si = 0.65_f64;
    let ci = (1.0 - si * si).sqrt();
    let ray = Ray {
        origin: Point3::new(0.0, 0.0, 0.0),
        dir: Vec3::new(si, 0.0, ci),
    };
    with_cx(|cx| {
        for wavelength in [600.0, 825.0, 1064.0, 1500.0] {
            let n = bulk.ior().eval(wavelength).unwrap();
            let ct = (1.0 - (si / n).powi(2)).sqrt();
            let first = native
                .fresnel_polarized(wavelength, false, ci, 1.0, n)
                .unwrap();
            let second = native
                .fresnel_polarized(wavelength, false, ct, n, 1.0)
                .unwrap();
            let transmitted = second
                .transmission
                .apply(first.transmission.apply(Stokes::unpolarized(1.0)));
            let beer = (-0.7 * 0.01 / ct).exp();
            let actual = lens.trace_ray(cx, ray, wavelength).unwrap().unwrap();
            let uncoated = bare.trace_ray(cx, ray, wavelength).unwrap().unwrap();
            close(actual.transmittance, transmitted.i * beer);
            assert_eq!(actual.ray, uncoated.ray);
            assert_eq!(actual.distance_m, uncoated.distance_m);
            assert_eq!(actual.optical_path_m, uncoated.optical_path_m);
            assert!(actual.transmittance < uncoated.transmittance);
            if wavelength == 1064.0 {
                let reversed = reversed.trace_ray(cx, ray, wavelength).unwrap().unwrap();
                assert!((actual.transmittance - reversed.transmittance).abs() > 1e-4);
                let independent_means =
                    first.transmission.0[0][0] * second.transmission.0[0][0] * beer;
                assert!((actual.transmittance - independent_means).abs() > 1e-6);
            }
        }
        for wavelength in [599.0, 1501.0] {
            assert_eq!(
                lens.trace_ray(cx, ray, wavelength),
                Err(CameraError::SequentialLensOutsideModel)
            );
        }
    });
}

#[test]
fn replacing_visible_coatings_restores_bulk_band_and_quarter_wave_ar() {
    let bare = biconvex();
    let surfaces: Vec<_> = bare
        .surfaces()
        .iter()
        .map(|surface| surface.with_coating(compact(100e-9)))
        .collect();
    let compact_lens = SequentialLens::try_new(550.0, &surfaces).unwrap();
    assert_eq!(compact_lens.wavelength_bounds_nm(), (400.0, 780.0));
    let film_n = 1.5_f64.sqrt();
    let ar = broad_coating(film_n, 1200e-9 / (4.0 * film_n), (400.0, 2000.0));
    let one = compact_lens
        .clone()
        .with_spectral_coating(0, ar.clone())
        .unwrap();
    assert_eq!(one.wavelength_bounds_nm(), (400.0, 780.0));
    let both = one.with_spectral_coating(1, ar).unwrap();
    assert_eq!(both.wavelength_bounds_nm(), (400.0, 1800.0));
    assert_eq!(both.reference_focus(), bare.reference_focus());
    assert_eq!(compact_lens.wavelength_bounds_nm(), (400.0, 780.0));
    with_cx(|cx| {
        let ray = Ray {
            origin: Point3::new(0.0, 0.0, 0.0),
            dir: Vec3::new(0.0, 0.0, 1.0),
        };
        let actual = both.trace_ray(cx, ray, 1200.0).unwrap().unwrap();
        let uncoated = bare.trace_ray(cx, ray, 1200.0).unwrap().unwrap();
        let beer = (-0.07_f64).exp();
        close(actual.transmittance, beer);
        close(uncoated.transmittance, 0.96_f64.powi(2) * beer);
        assert_eq!(actual.ray, uncoated.ray);
        assert_eq!(actual.optical_path_m, uncoated.optical_path_m);
    });
}

#[test]
fn only_effective_films_control_timing_and_domain_before_clipping() {
    let bare = biconvex();
    let surfaces: Vec<_> = bare
        .surfaces()
        .iter()
        .map(|surface| surface.with_coating(compact(100e-9)))
        .collect();
    let compact_lens = SequentialLens::try_new(550.0, &surfaces).unwrap();
    let zero = broad_coating(1.3, 0.0, (450.0, 1600.0));
    let zero = compact_lens
        .with_spectral_coating(0, zero.clone())
        .unwrap()
        .with_spectral_coating(1, zero)
        .unwrap();
    assert_eq!(zero.wavelength_bounds_nm(), (450.0, 1600.0));
    let active = zero
        .clone()
        .with_spectral_coating(0, broad_coating(1.3, 100e-9, (450.0, 1600.0)))
        .unwrap();
    let (bare, zero, active) = (camera(bare), camera(zero), camera(active));
    with_cx(|cx| {
        let expected = bare
            .generate_sequential_spectral_timed_beams_from_tangent_offsets(
                cx,
                0.0,
                0.0,
                LensSample::CENTER,
                &[1200.0],
            )
            .unwrap()[0]
            .unwrap();
        let actual = zero
            .generate_sequential_spectral_timed_beams_from_tangent_offsets(
                cx,
                0.0,
                0.0,
                LensSample::CENTER,
                &[1200.0],
            )
            .unwrap()[0]
            .unwrap();
        assert_eq!(actual.ray(), expected.ray());
        assert_eq!(actual.power(), expected.power());
        assert_eq!(actual.group_delay_s(), expected.group_delay_s());
        for tangent in [0.0, 20.0] {
            assert_eq!(
                active.generate_sequential_spectral_timed_beams_from_tangent_offsets(
                    cx,
                    tangent,
                    0.0,
                    LensSample::CENTER,
                    &[1200.0],
                ),
                Err(CameraError::SequentialLensOutsideModel)
            );
            assert_eq!(
                zero.generate_sequential_spectral_beams_from_tangent_offsets(
                    cx,
                    tangent,
                    0.0,
                    LensSample::CENTER,
                    &[1700.0],
                ),
                Err(CameraError::SequentialLensOutsideModel)
            );
        }
    });
    assert_eq!(
        LensCoating::try_new(&[]),
        Err(CameraError::InvalidSequentialLens)
    );
    let nir_only = broad_coating(1.3, 0.0, (1000.0, 1600.0));
    assert_eq!(
        biconvex().with_spectral_coating(0, nir_only),
        Err(CameraError::InvalidSequentialLens)
    );
}
