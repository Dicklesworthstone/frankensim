//! Full Stokes coupling through the real lens, not independent scalar losses.

use super::*;
use crate::camera::{LensSurface, SequentialLens};
use crate::charts::Ray;
use crate::conductor::{CoatingIorSample, ConductorDataStatus, ConductorSource};
use crate::dielectric::{BeerLambertAbsorption, DielectricGlass, GlassProvenance, SellmeierIor};
use crate::polarization::Stokes;
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_geom::{Point3, Vec3};

#[test]
fn mixed_nonvisible_films_retain_signed_phase_for_all_stokes_components() {
    let source =
        ConductorSource::try_new(ContentHash([7; 32]), ConductorDataStatus::Representative)
            .unwrap();
    let transparent = CoherentDielectricCoating::try_new(
        SellmeierIor::try_new_in_band([0.9, 0.05, 0.0], [0.02, 0.09, 0.0], (800.0, 1600.0))
            .unwrap(),
        170e-9,
        source,
    )
    .unwrap();
    let absorbing = CoherentAbsorbingCoating::try_new_in_band(
        std::array::from_fn(|i| {
            CoatingIorSample::try_new_in_band(
                800.0 + 100.0 * i as f64,
                2.1 - 0.03 * i as f64,
                0.04 + 0.01 * i as f64,
            )
            .unwrap()
        }),
        130e-9,
        100.0,
        source,
    )
    .unwrap();
    let layers: [LensCoatingLayer; 2] = [transparent.into(), absorbing.into()];
    let coating = LensCoating::try_new(&layers).unwrap();
    let borrowed = [
        CoatingLayer::Dielectric(&transparent),
        CoatingLayer::Absorbing(&absorbing),
    ];
    let native = DielectricCoatingStack::try_mixed_in_band(&borrowed, (800.0, 1600.0)).unwrap();
    let glass = DielectricGlass::new(
        SellmeierIor::try_new_in_band([1.25, 0.0, 0.0], [0.0; 3], (800.0, 1600.0)).unwrap(),
        BeerLambertAbsorption::try_constant(0.4).unwrap(),
        GlassProvenance::Custom,
    );
    let bare = SequentialLens::try_new(
        1064.0,
        &[
            LensSurface::try_new(0.02, 0.0, 0.1, Some(glass)).unwrap(),
            LensSurface::try_new(0.05, 0.0, 0.1, None).unwrap(),
        ],
    )
    .unwrap();
    let seed = Stokes {
        i: 0.6,
        q: 0.1,
        u: 0.2,
        v: 0.25,
    };
    let ci = 0.8;
    let ct = (1.0_f64 - (0.6_f64 / 1.5).powi(2)).sqrt();
    let first = native
        .fresnel_polarized(1064.0, false, ci, 1.0, 1.5)
        .unwrap();
    assert!(first.transmission.0[2][3].abs() > 1e-4);
    let attenuation = (-0.4 * 0.03 / ct).exp();
    let gate = CancelGate::new_clock_free();
    fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default()).scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 47,
                kernel_id: 11,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        for reverse_front in [false, true] {
            let front = if reverse_front {
                LensCoating::try_new(&[layers[1], layers[0]]).unwrap()
            } else {
                coating.clone()
            };
            let lens = bare
                .clone()
                .with_spectral_coating(0, coating.clone())
                .unwrap()
                .with_spectral_coating(1, front)
                .unwrap();
            let second = native
                .fresnel_polarized(1064.0, reverse_front, ct, 1.5, 1.0)
                .unwrap();
            let expected = second.transmission.apply(first.transmission.apply(seed));
            let actual = lens
                .trace_seeded_ray(
                    &cx,
                    Ray {
                        origin: Point3::new(0.0, 0.0, 0.0),
                        dir: Vec3::new(0.6, 0.0, ci),
                    },
                    1064.0,
                    Some(Vec3::new(0.0, 1.0, 0.0)),
                    seed,
                    false,
                )
                .unwrap()
                .unwrap();
            for (actual, expected) in [
                actual.power.i,
                actual.power.q,
                actual.power.u,
                actual.power.v,
            ]
            .into_iter()
            .zip([expected.i, expected.q, expected.u, expected.v])
            {
                assert!(
                    (actual - expected * attenuation).abs() < 8e-13,
                    "{actual:e} != {expected:e}"
                );
            }
            assert!((actual.distance_m - (0.02 / ci + 0.03 / ct)).abs() < 8e-15);
        }
    });
}
