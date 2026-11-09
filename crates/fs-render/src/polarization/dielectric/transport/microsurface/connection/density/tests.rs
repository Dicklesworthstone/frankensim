use super::*;
use crate::conductor::{ConductorDataStatus, ConductorIorSample, ConductorOptics, ConductorSource};
use crate::dielectric::{BeerLambertAbsorption, CauchyIor, DielectricGlass, GlassProvenance};
use crate::polarization::Stokes;
use crate::polarization::dielectric::transport::microsurface::SmithWalkSettings;
use crate::surface_detail::BrushingPattern;
use core::f64::consts::{PI, TAU};
use fs_exec::{CancelGate, Cx, ExecMode, StreamKey};
use fs_rand::philox::philox4x32_10;

const Z: Vec3 = Vec3::new(0.0, 0.0, 1.0);
fn input() -> PolarizedBeam {
    PolarizedBeam::try_new(
        Point3::new(0.0, 0.0, 0.0),
        PolarizationFrame::try_new(Z.scale(-1.0), Vec3::new(1.0, 0.0, 0.0)).unwrap(),
        550.0,
        DielectricGlass::new(
            CauchyIor::try_constant(1.0).unwrap(),
            BeerLambertAbsorption::CLEAR,
            GlassProvenance::Custom,
        ),
        Stokes::try_new(1.0, 0.3, -0.4, 0.5).unwrap(),
    )
    .unwrap()
}
fn metal() -> ConductorOptics {
    ConductorOptics::try_new(
        core::array::from_fn(|i| {
            ConductorIorSample::try_new(380.0 + 50.0 * i as f64, 0.2, 3.0).unwrap()
        }),
        ConductorSource::try_new(
            fs_blake3::ContentHash([47; 32]),
            ConductorDataStatus::Representative,
        )
        .unwrap(),
    )
    .unwrap()
}
fn with_cx(f: impl FnOnce(&Cx<'_>)) {
    let gate = CancelGate::new();
    fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default()).scope(|arena| {
        f(&Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 41,
                kernel_id: 7,
                tile: 0,
                iteration: 0,
            },
            asupersync::types::Budget::INFINITE,
            ExecMode::Deterministic,
        ));
    });
}
fn near(a: f64, b: f64) {
    assert!((a - b).abs() < 5e-12 * b.abs().max(1.0), "{a} != {b}");
}
fn draw(path: u32, facet: u32, domain: u32) -> [f64; 3] {
    let b = philox4x32_10([path, facet, domain, 0], [17, 43]);
    core::array::from_fn(|i| (f64::from(b[i]) + 0.5) / 4_294_967_296.0)
}

#[test]
fn g1_first_escape_density_integrates_to_the_analytic_single_reflection_probability() {
    let beam = input();
    let surface = ConductorSurface::try_rough(1.0).unwrap();
    let count = 16_384;
    let mut integral = 0.0;
    for i in 0..count {
        let cosine = (i as f64 + 0.5) / count as f64;
        let out = Vec3::new((1.0 - cosine * cosine).sqrt(), 0.0, cosine);
        let q = beam
            .smith_conductor_first_reflection_density(Z, surface, out)
            .unwrap();
        // Independent closed form for alpha=1 and normal incidence.
        near(q, cosine / (2.0 * PI * (1.0 + cosine)));
        integral += q * TAU / count as f64;
    }
    assert!((integral - (1.0 - fs_math::det::ln(2.0))).abs() < 2e-9);
    assert!(
        integral < 0.31,
        "this subdensity must never be normalized to one"
    );
    assert_eq!(
        beam.smith_conductor_first_reflection_density(Z, surface, Z.scale(-1.0))
            .unwrap(),
        0.0
    );
    assert_eq!(
        beam.smith_conductor_first_reflection_density(Z, ConductorSurface::SMOOTH, Z)
            .unwrap(),
        0.0
    );
    assert!(
        beam.smith_conductor_first_reflection_density(Z, surface, Vec3::new(f64::NAN, 0.0, 1.0))
            .is_err()
    );
}

#[test]
fn g3_response_and_endpoint_share_a_deterministic_density_in_original_brushed_frames() {
    with_cx(|cx| {
        let beam = input();
        let optics = metal();
        let out = Vec3::new(0.3, -0.4, (0.75_f64).sqrt());
        let source = PolarizationFrame::try_new(out.scale(-1.0), Vec3::new(1.0, 0.0, 0.0)).unwrap();
        let point = Point3::new(0.2, 0.1, 0.0);
        let settings = SmithWalkSettings::try_new(64).unwrap();
        let ordinary = ConductorSurface::try_rough(0.25).unwrap();
        let expected = beam
            .smith_conductor_first_reflection_density(Z, ordinary, out)
            .unwrap();
        for detail in [
            SurfaceDetail::none(0.25),
            SurfaceDetail::try_new(0.25, 4.0, BrushingPattern::Radial, 0.0, None).unwrap(),
            SurfaceDetail::try_new(0.25, 4.0, BrushingPattern::Radial, 0.2, None).unwrap(),
        ] {
            let q = beam
                .instanced_smith_brushed_first_reflection_density(
                    &detail,
                    RigidTransform::identity(),
                    point,
                    Z,
                    out,
                )
                .unwrap();
            for history in 0..16 {
                let response = beam
                    .evaluate_instanced_smith_brushed_conductor_connection(
                        cx,
                        &optics,
                        &detail,
                        RigidTransform::identity(),
                        point,
                        Z,
                        source,
                        Stokes::unpolarized(1.0),
                        settings,
                        |facet| draw(history, facet, 1),
                    )
                    .unwrap();
                near(response.first_reflection_density(), q);
            }
        }
        let fallback = beam
            .instanced_smith_brushed_first_reflection_density(
                &SurfaceDetail::none(0.25),
                RigidTransform::identity(),
                point,
                Z,
                out,
            )
            .unwrap();
        assert_eq!(fallback, expected);
        let direct = beam
            .evaluate_smith_conductor_connection(
                cx,
                Z,
                &optics,
                ordinary,
                source,
                Stokes::unpolarized(1.0),
                settings,
                |_| [0.5, 0.0, 0.0],
            )
            .unwrap();
        assert_eq!(direct.first_reflection_density(), expected);
    });
}

#[test]
fn g1_full_native_polarized_reflection_and_direction_only_partition_have_same_expectation() {
    with_cx(|cx| {
        let beam = input();
        let optics = metal();
        let surface = ConductorSurface::try_rough(0.7).unwrap();
        let settings = SmithWalkSettings::try_new(128).unwrap();
        let n = 32_768;
        let p = 1.0 / (2.0 * PI);
        let mut difference = 0.0;
        let mut square = 0.0;
        let mut higher = 0;
        for path in 0..n {
            let continuation = beam
                .sample_smith_conductor(cx, Z, &optics, surface, settings, |facet| {
                    draw(path, facet, 2)
                })
                .unwrap();
            higher += u32::from(continuation.facet_events() > 1);
            let escaped = continuation.beam().ray().dir;
            let q = beam
                .smith_conductor_first_reflection_density(Z, surface, escaped)
                .unwrap();
            let w_camera = q * q / (p * p + q * q);
            // Independent unit-hemisphere source proposal, with a nonconstant
            // radiance field. Both strategies integrate this same physical sky.
            let values = draw(path, 0, 3);
            let c = values[0];
            let (sine, cosine) = fs_math::det::sin_cos(TAU * values[1]);
            let radial = (1.0 - c * c).sqrt();
            let out = Vec3::new(radial * cosine, radial * sine, c);
            let frame =
                PolarizationFrame::try_new(out.scale(-1.0), Vec3::new(0.0, 1.0, 0.0)).unwrap();
            let response = beam
                .evaluate_smith_conductor_connection(
                    cx,
                    Z,
                    &optics,
                    surface,
                    frame,
                    Stokes::unpolarized(1.0),
                    settings,
                    |facet| draw(path, facet, 4),
                )
                .unwrap_or_else(|error| panic!("path {path}, source cosine {c}: {error:?}"));
            let q = response.first_reflection_density();
            let w_light = p * p / (p * p + q * q);
            let radiance = |direction: Vec3| 1.0 + 3.0 * direction.x * direction.x;
            let path_only = continuation.beam().power().i * radiance(escaped);
            let combined =
                path_only * w_camera + response.factor_per_sr() * c * radiance(out) * w_light / p;
            let d = combined - path_only;
            difference += d;
            square += d * d;
        }
        let mean = difference / n as f64;
        let variance = (square / n as f64 - mean * mean).max(0.0);
        let error = (variance / n as f64).sqrt();
        assert!(
            higher > n / 4,
            "exercise multiple reflections, not just the analytic first order"
        );
        assert!(
            mean.abs() < 6.0 * error + 5e-4,
            "mean delta {mean}, standard error {error}"
        );
    });
}

#[test]
fn g0_grazing_connection_terms_survive_until_the_complete_positive_sum_is_rounded() {
    use super::super::{PositiveTerm, Sum};
    let model = SmithGgx::try_new(1.0, 1.0).unwrap();
    let incoming = Z.scale(-1.0);
    let out = Vec3::new((1.0 - 1e-6_f64).sqrt(), 0.0, 0.001);
    assert!(matches!(
        model.connection_at_depth(incoming, 3.0, out),
        Err(crate::conductor::microsurface::SmithError::NumericalRange)
    ));
    let connection = model
        .log_connection_at_depth(incoming, 3.0, out)
        .unwrap()
        .unwrap();
    // For alpha=1, A(-incoming)=1 and Lambda(out)=(1-c)/(2c).
    let expected = -fs_math::det::ln(4.0 * PI * out.z) - 3.0 * (1.0 - out.z) / (2.0 * out.z);
    near(connection.log_factor_per_sr, expected);
    let mut sum = Sum::default();
    sum.add(PositiveTerm::from_log_factor(connection.log_factor_per_sr, 0.5).unwrap())
        .unwrap();
    assert!(
        sum.value().is_err(),
        "a wholly unrepresentable positive response must refuse"
    );
    sum.add(PositiveTerm::Linear(0.25)).unwrap();
    assert_eq!(sum.value().unwrap(), 0.25);
    let mut subnormals = Sum::default();
    for _ in 0..16 {
        subnormals.add(PositiveTerm::Log(-745.5)).unwrap();
    }
    let expected = fs_math::det::exp(-745.5 + fs_math::det::ln(16.0));
    assert!(expected > 0.0);
    assert_eq!(subnormals.value().unwrap(), expected);
    let normal = model
        .log_connection_at_depth(incoming, 0.3, Z)
        .unwrap()
        .unwrap();
    near(
        fs_math::det::exp(normal.log_factor_per_sr),
        model
            .connection_at_depth(incoming, 0.3, Z)
            .unwrap()
            .unwrap()
            .factor_per_sr,
    );
}
