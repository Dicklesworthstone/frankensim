//! Conformance battery for as-built FIELD extraction (bead
//! `frankensim-extreal-program-f85xj.12.3`).
//!
//! The tests are organized around what could silently go wrong rather than
//! around the API surface:
//!
//! - G0 algebra: deviations equal injected offsets; the composed half-width is
//!   exactly the declared quadrature; roughness matches closed forms.
//! - INDEPENDENT CHECK: the registration sensitivity is verified against a
//!   finite difference of the deviation under an explicitly reconstructed
//!   perturbed pose, so a wrong covariance pivot cannot pass.
//! - ADVERSARIAL: near-edge-on surfaces are refused and counted, not silently
//!   measured; a rank-deficient station set refuses instead of returning an
//!   arbitrary member of the solution family; unfiltered roughness is shown to
//!   actually count waviness, proving the documented boundary is real.
//! - G5 determinism: identical inputs give identical values and identity, and
//!   content changes move the identity.

use fs_asbuilt::field::{
    DeviationField, FieldError, ProbeModel, ProfileForm, SurfaceSample, ThicknessField, fit_form,
    profile_statistics,
};
use fs_asbuilt::propagate::CoveragePolicy;
use fs_asbuilt::rigid3::{
    CalibratedRigid3Registration, Covariance3, CrossFiducialModel3, Fiducial3, MetrologyModel3,
    Point3, estimate_calibrated_rigid3,
};
use fs_asbuilt::uncertainty::HuberPolicy;
use fs_evidence::uncertainty::TermValue;
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};

// ---------------------------------------------------------------------------
// Harness.
// ---------------------------------------------------------------------------

fn with_cx<R>(cancelled: bool, f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate = CancelGate::new_clock_free();
    if cancelled {
        gate.request();
    }
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    let clock = fs_exec::VirtualClock::new();
    let result = pool.scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 0x0C12_0300,
                kernel_id: 12,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        )
        .with_time_source(&clock);
        f(&cx)
    });
    assert!(
        pool.stats().quiescent(),
        "arena must be quiescent after scope"
    );
    result
}

fn with_default_cx<R>(f: impl FnOnce(&Cx<'_>) -> R) -> R {
    with_cx(false, f)
}

fn p3(x: f64, y: f64, z: f64) -> Point3 {
    Point3::new(x, y, z).expect("finite test point")
}

/// Rodrigues rotation about a unit axis.
fn rodrigues(axis: [f64; 3], angle: f64) -> [[f64; 3]; 3] {
    let norm = (axis[0] * axis[0] + axis[1] * axis[1] + axis[2] * axis[2]).sqrt();
    let a = [axis[0] / norm, axis[1] / norm, axis[2] / norm];
    let (s, c) = angle.sin_cos();
    let v = 1.0 - c;
    [
        [
            c + v * a[0] * a[0],
            v * a[0] * a[1] - s * a[2],
            v * a[0] * a[2] + s * a[1],
        ],
        [
            v * a[0] * a[1] + s * a[2],
            c + v * a[1] * a[1],
            v * a[1] * a[2] - s * a[0],
        ],
        [
            v * a[0] * a[2] - s * a[1],
            v * a[1] * a[2] + s * a[0],
            c + v * a[2] * a[2],
        ],
    ]
}

fn apply(rotation: &[[f64; 3]; 3], translation: [f64; 3], point: [f64; 3]) -> [f64; 3] {
    let mut out = [0.0f64; 3];
    for row in 0..3 {
        out[row] = rotation[row][0] * point[0]
            + rotation[row][1] * point[1]
            + rotation[row][2] * point[2]
            + translation[row];
    }
    out
}

/// Eight design fiducials on a box plus an interior point: a well-conditioned,
/// non-coplanar set.
fn fiducial_design() -> [[f64; 3]; 8] {
    [
        [0.0, 0.0, 0.0],
        [2.0, 0.0, 0.0],
        [0.0, 3.0, 0.0],
        [2.0, 3.0, 0.0],
        [0.0, 0.0, 1.0],
        [2.0, 0.0, 1.0],
        [0.0, 3.0, 1.0],
        [0.5, 1.0, 0.6],
    ]
}

/// Build a calibrated registration for a declared true pose, with EXACT
/// fiducial measurements. The covariance is then driven purely by the declared
/// metrology model, which keeps the sensitivity tests analytic.
fn registration_for(
    rotation: &[[f64; 3]; 3],
    translation: [f64; 3],
    sigma: f64,
    cx: &Cx<'_>,
) -> CalibratedRigid3Registration {
    let mut fiducials = Vec::new();
    for design in fiducial_design() {
        let measured = apply(rotation, translation, design);
        fiducials.push(Fiducial3::new(
            p3(design[0], design[1], design[2]),
            p3(measured[0], measured[1], measured[2]),
        ));
    }
    let variance = sigma * sigma;
    let covariances = vec![
        Covariance3::new(variance, 0.0, 0.0, variance, 0.0, variance)
            .expect("isotropic covariance");
        fiducials.len()
    ];
    let model = MetrologyModel3::new(
        covariances,
        CrossFiducialModel3::Independent,
        HuberPolicy::Disabled,
        "field-test-calibration/rev1",
    )
    .expect("metrology model");
    estimate_calibrated_rigid3(&fiducials, &model, cx).expect("calibrated registration")
}

fn identity_registration(cx: &Cx<'_>) -> CalibratedRigid3Registration {
    registration_for(
        &rodrigues([0.0, 0.0, 1.0], 0.0),
        [0.0, 0.0, 0.0],
        1.0e-3,
        cx,
    )
}

fn coverage() -> CoveragePolicy {
    CoveragePolicy::new(0.95, 2.0).expect("coverage policy")
}

/// A probe looking straight down -z, so an upward-facing plate is at normal
/// incidence.
fn down_probe(lateral_sigma: f64, range_sigma: f64) -> ProbeModel {
    ProbeModel::new([0.0, 0.0, -1.0], 0.1, lateral_sigma, range_sigma).expect("probe model")
}

/// Upward-facing plate stations on a grid, each displaced outward by
/// `deviation(x, y)`.
fn plate_samples(
    stations: &[[f64; 2]],
    deviation: impl Fn(f64, f64) -> f64,
) -> (Vec<SurfaceSample>, Vec<[f64; 2]>) {
    let normal = [0.0, 0.0, 1.0];
    let mut samples = Vec::new();
    for station in stations {
        let (x, y) = (station[0], station[1]);
        let d = deviation(x, y);
        samples.push(SurfaceSample::new(p3(x, y, 0.0), normal, p3(x, y, d)).expect("plate sample"));
    }
    (samples, stations.to_vec())
}

fn grid_stations() -> Vec<[f64; 2]> {
    let mut stations = Vec::new();
    for i in 0..5 {
        for j in 0..5 {
            stations.push([f64::from(i) * 0.5, f64::from(j) * 0.75]);
        }
    }
    stations
}

// ---------------------------------------------------------------------------
// G0: deviation algebra.
// ---------------------------------------------------------------------------

#[test]
fn deviation_recovers_an_injected_analytic_offset() {
    with_default_cx(|cx| {
        let registration = identity_registration(cx);
        let probe = down_probe(0.0, 0.0);
        let stations = grid_stations();
        // A deviation that varies over the plate so a constant-offset bug
        // cannot pass.
        let (samples, _) = plate_samples(&stations, |x, y| 0.01 * x - 0.004 * y + 0.002);

        let field = DeviationField::extract(&samples, &registration, &probe, coverage(), cx)
            .expect("field extracts");

        assert_eq!(field.points().len(), stations.len());
        assert!(field.grazing().is_empty());
        for (point, station) in field.points().iter().zip(stations.iter()) {
            let expected = 0.01 * station[0] - 0.004 * station[1] + 0.002;
            assert!(
                (point.deviation() - expected).abs() < 1e-12,
                "deviation {} != injected {expected}",
                point.deviation()
            );
            // Normal incidence: the incidence cosine is exactly one.
            assert!((point.incidence_cosine() - 1.0).abs() < 1e-12);
        }
    });
}

#[test]
fn deviation_is_measured_along_the_registered_normal_under_a_rotated_pose() {
    with_default_cx(|cx| {
        // A general pose: the deviation must still be the outward normal
        // component, which only works if the normal is rotated with the part.
        let rotation = rodrigues([0.3, -0.5, 0.81], 0.4);
        let translation = [1.25, -0.75, 2.5];
        let registration = registration_for(&rotation, translation, 1.0e-3, cx);

        let design = [0.5, 1.0, 0.0];
        let normal = [0.0, 0.0, 1.0];
        let injected = 0.037;
        // Displace the measured point along the ROTATED normal, plus a purely
        // tangential slide that must not register as deviation.
        let rotated_normal = apply(&rotation, [0.0, 0.0, 0.0], normal);
        let tangent = apply(&rotation, [0.0, 0.0, 0.0], [1.0, 0.0, 0.0]);
        let base = apply(&rotation, translation, design);
        let measured = [
            base[0] + injected * rotated_normal[0] + 0.21 * tangent[0],
            base[1] + injected * rotated_normal[1] + 0.21 * tangent[1],
            base[2] + injected * rotated_normal[2] + 0.21 * tangent[2],
        ];

        let sample = SurfaceSample::new(
            p3(design[0], design[1], design[2]),
            normal,
            p3(measured[0], measured[1], measured[2]),
        )
        .expect("sample");
        // A probe aligned with the rotated normal keeps this at normal
        // incidence.
        let probe = ProbeModel::new(
            [-rotated_normal[0], -rotated_normal[1], -rotated_normal[2]],
            0.1,
            0.0,
            0.0,
        )
        .expect("probe");

        let field = DeviationField::extract(&[sample], &registration, &probe, coverage(), cx)
            .expect("field extracts");
        let deviation = field.points()[0].deviation();
        assert!(
            (deviation - injected).abs() < 1e-9,
            "rotated-pose deviation {deviation} != injected {injected}; \
             a tangential slide must not appear as deviation"
        );
    });
}

// ---------------------------------------------------------------------------
// INDEPENDENT CHECK: the pose sensitivity and its pivot.
// ---------------------------------------------------------------------------

#[test]
fn pose_sensitivity_matches_a_finite_difference_about_the_covariance_pivot() {
    // The crate documents that a sensitivity computed about a different pivot
    // is SILENTLY WRONG. This test is the independent oracle: it reconstructs
    // the perturbed pose explicitly from the documented parameterization (a
    // left rotation-vector perturbation about the image of the weighted design
    // centroid), recomputes the deviation exactly, and compares against the
    // analytic gradient. A wrong pivot `c'` produces a FIRST-ORDER error
    // `dr . ((c' - c) x n)`, which this tolerance rejects.
    with_default_cx(|cx| {
        let rotation = rodrigues([0.2, 0.7, -0.68], 0.35);
        let translation = [0.4, 1.1, -0.6];
        let registration = registration_for(&rotation, translation, 1.0e-3, cx);
        let pivot = registration.rotation_pivot();
        let fitted_rotation = *registration.registration().rotation();
        let fitted_translation = registration.registration().translation();

        // A station well away from the pivot, so the moment arm is large and a
        // pivot error cannot hide.
        let design = [2.0, 3.0, 1.0];
        let normal = {
            let n = [0.6f64, -0.8, 0.0];
            let norm = (n[0] * n[0] + n[1] * n[1] + n[2] * n[2]).sqrt();
            [n[0] / norm, n[1] / norm, n[2] / norm]
        };
        let gradient = registration
            .normal_deviation_sensitivity(p3(design[0], design[1], design[2]), normal)
            .expect("sensitivity");

        let delta_t = [7.0e-7f64, -3.0e-7, 5.0e-7];
        let delta_r = [2.0e-7f64, 6.0e-7, -4.0e-7];

        // Nominal image and the pivot image.
        let nominal_image = apply(&fitted_rotation, fitted_translation, design);
        let pivot_image = apply(&fitted_rotation, fitted_translation, pivot);
        // Perturbed image: rotate about the pivot image by delta_r, then
        // translate by delta_t.
        let magnitude =
            (delta_r[0] * delta_r[0] + delta_r[1] * delta_r[1] + delta_r[2] * delta_r[2]).sqrt();
        let perturbation = rodrigues(delta_r, magnitude);
        let relative = [
            nominal_image[0] - pivot_image[0],
            nominal_image[1] - pivot_image[1],
            nominal_image[2] - pivot_image[2],
        ];
        let rotated_relative = apply(&perturbation, [0.0, 0.0, 0.0], relative);
        let perturbed_image = [
            pivot_image[0] + rotated_relative[0] + delta_t[0],
            pivot_image[1] + rotated_relative[1] + delta_t[1],
            pivot_image[2] + rotated_relative[2] + delta_t[2],
        ];

        // The deviation of a FIXED measured point moves by minus the image
        // displacement projected on the registered normal.
        let rotated_normal = apply(&fitted_rotation, [0.0, 0.0, 0.0], normal);
        let exact_change = -((perturbed_image[0] - nominal_image[0]) * rotated_normal[0]
            + (perturbed_image[1] - nominal_image[1]) * rotated_normal[1]
            + (perturbed_image[2] - nominal_image[2]) * rotated_normal[2]);
        // The gradient is d(deviation)/d(pose); the image moving by +d moves
        // the deviation by -d, so compare against the negated linear form.
        let predicted = -(gradient[0] * delta_t[0]
            + gradient[1] * delta_t[1]
            + gradient[2] * delta_t[2]
            + gradient[3] * delta_r[0]
            + gradient[4] * delta_r[1]
            + gradient[5] * delta_r[2]);

        let residual = (exact_change - predicted).abs();
        // Second-order truncation is O(|delta|^2) ~ 1e-13; a pivot error would
        // be first order, ~1e-7 * arm.
        assert!(
            residual < 1.0e-12,
            "sensitivity residual {residual:e} exceeds second-order truncation; \
             exact {exact_change:e} vs predicted {predicted:e}"
        );
    });
}

#[test]
fn a_pure_translation_moves_every_deviation_by_the_same_normal_component() {
    // Sanity companion to the finite-difference test: with no rotation the
    // pivot is irrelevant, so any translation-block error shows up alone.
    with_default_cx(|cx| {
        let registration = identity_registration(cx);
        for design in [[0.0, 0.0, 0.0], [2.0, 3.0, 1.0], [0.5, 1.0, 0.6]] {
            let gradient = registration
                .normal_deviation_sensitivity(p3(design[0], design[1], design[2]), [0.0, 0.0, 1.0])
                .expect("sensitivity");
            assert!((gradient[0]).abs() < 1e-12);
            assert!((gradient[1]).abs() < 1e-12);
            assert!(
                (gradient[2] - 1.0).abs() < 1e-12,
                "a +z normal must have unit tz sensitivity"
            );
        }
    });
}

// ---------------------------------------------------------------------------
// Uncertainty composition.
// ---------------------------------------------------------------------------

#[test]
fn half_width_is_exactly_the_declared_coverage_scaled_quadrature() {
    with_default_cx(|cx| {
        let registration = identity_registration(cx);
        let probe = down_probe(2.0e-4, 5.0e-4);
        let stations = grid_stations();
        let (samples, _) = plate_samples(&stations, |x, _| 0.001 * x);

        let field = DeviationField::extract(&samples, &registration, &probe, coverage(), cx)
            .expect("field extracts");

        for point in field.points() {
            let (scan, registration_sigma, ambiguity) = (
                point.scan_sigma(),
                point.registration_sigma(),
                point.ambiguity_sigma(),
            );
            let expected = coverage().factor()
                * (scan * scan + registration_sigma * registration_sigma + ambiguity * ambiguity)
                    .sqrt();
            assert!(
                (point.half_width() - expected).abs() <= 1e-15 * expected.max(1.0),
                "half width {} is not the declared quadrature {expected}",
                point.half_width()
            );
            assert!(
                point.registration_sigma() > 0.0,
                "a calibrated pose must contribute a positive sigma"
            );
        }
    });
}

#[test]
fn oblique_incidence_shrinks_range_noise_and_grows_correspondence_ambiguity() {
    // The two terms move in OPPOSITE directions with incidence. A model that
    // simply scaled everything by 1/c would pass a one-sided test.
    with_default_cx(|cx| {
        let registration = identity_registration(cx);
        let probe = down_probe(1.0e-4, 1.0e-3);

        let mut previous_scan = f64::INFINITY;
        let mut previous_ambiguity = 0.0f64;
        // Tilt the surface normal progressively away from the probe axis.
        for step in 0..5 {
            let angle = 0.2 * f64::from(step);
            let normal = [angle.sin(), 0.0, angle.cos()];
            let sample =
                SurfaceSample::new(p3(0.5, 0.5, 0.0), normal, p3(0.5, 0.5, 0.0)).expect("sample");
            let field = DeviationField::extract(&[sample], &registration, &probe, coverage(), cx)
                .expect("field extracts");
            let point = field.points()[0];

            assert!(
                point.scan_sigma() < previous_scan,
                "range-noise contribution must shrink as incidence becomes oblique"
            );
            assert!(
                point.ambiguity_sigma() >= previous_ambiguity,
                "correspondence ambiguity must grow as incidence becomes oblique"
            );
            previous_scan = point.scan_sigma();
            previous_ambiguity = point.ambiguity_sigma();
        }
        assert!(
            previous_ambiguity > 0.0,
            "an oblique sample must carry a positive ambiguity term"
        );
    });
}

#[test]
fn a_zero_noise_probe_leaves_only_the_registration_term() {
    with_default_cx(|cx| {
        let registration = identity_registration(cx);
        let probe = down_probe(0.0, 0.0);
        let sample = SurfaceSample::new(p3(1.0, 1.0, 0.0), [0.0, 0.0, 1.0], p3(1.0, 1.0, 0.01))
            .expect("sample");
        let field = DeviationField::extract(&[sample], &registration, &probe, coverage(), cx)
            .expect("field extracts");
        let point = field.points()[0];
        assert_eq!(point.scan_sigma(), 0.0);
        assert_eq!(point.ambiguity_sigma(), 0.0);
        assert!(point.registration_sigma() > 0.0);
        assert!(
            (point.half_width() - coverage().factor() * point.registration_sigma()).abs() < 1e-18
        );
    });
}

#[test]
fn the_measurement_term_brackets_the_half_width_not_the_value() {
    with_default_cx(|cx| {
        let registration = identity_registration(cx);
        let probe = down_probe(1.0e-4, 5.0e-4);
        let stations = grid_stations();
        // Tilted stations give a spread of incidence, hence a real bracket.
        let (mut samples, _) = plate_samples(&stations, |x, _| 0.001 * x);
        samples.push(
            SurfaceSample::new(
                p3(9.0, 9.0, 0.0),
                [0.5, 0.0, 0.75_f64.sqrt()],
                p3(9.0, 9.0, 0.0),
            )
            .expect("tilted sample"),
        );

        let field = DeviationField::extract(&samples, &registration, &probe, coverage(), cx)
            .expect("field extracts");
        let term = field
            .measurement_term("as-built-field-uncertainty")
            .expect("term admits");

        // The bracket is [min, max] of the per-point half-widths. It brackets
        // the HALF-WIDTH itself, not the deviation value, so aggregation reads
        // the UPPER endpoint; a midpoint reading would report less uncertainty
        // than the field carries.
        let TermValue::IntervalBound { lower, upper } = term.value() else {
            panic!("a composed field term must be an interval bound");
        };
        assert!(
            (*upper - field.max_half_width()).abs() < 1e-18,
            "the upper endpoint must be the largest per-point half-width"
        );
        assert!(
            (*lower - field.min_half_width()).abs() < 1e-18,
            "the lower endpoint must be the smallest per-point half-width"
        );
        assert!(
            field.min_half_width() < field.max_half_width(),
            "a spread of incidence must produce a real bracket, not a point"
        );
    });
}

// ---------------------------------------------------------------------------
// ADVERSARIAL: grazing incidence.
// ---------------------------------------------------------------------------

#[test]
fn near_edge_on_walls_are_flagged_and_counted_not_dropped() {
    with_default_cx(|cx| {
        let registration = identity_registration(cx);
        let probe = down_probe(1.0e-4, 1.0e-4);

        // Three good upward stations and two near-vertical wall stations whose
        // normals are almost perpendicular to the probe axis.
        let mut samples = Vec::new();
        for i in 0..3 {
            samples.push(
                SurfaceSample::new(
                    p3(f64::from(i), 0.0, 0.0),
                    [0.0, 0.0, 1.0],
                    p3(f64::from(i), 0.0, 0.005),
                )
                .expect("flat sample"),
            );
        }
        for i in 0..2 {
            // |n . v| = 0.05, below the 0.1 floor.
            let normal = [(1.0f64 - 0.05 * 0.05).sqrt(), 0.0, 0.05];
            samples.push(
                SurfaceSample::new(
                    p3(f64::from(i), 5.0, 0.0),
                    normal,
                    p3(f64::from(i), 5.0, 0.0),
                )
                .expect("wall sample"),
            );
        }

        let field = DeviationField::extract(&samples, &registration, &probe, coverage(), cx)
            .expect("field extracts");

        assert_eq!(
            field.points().len(),
            3,
            "only the flat stations are measured"
        );
        assert_eq!(
            field.grazing(),
            &[3, 4],
            "refused stations are retained by SUPPLIED index, not discarded"
        );
        assert_eq!(
            field.sample_count(),
            5,
            "admitted plus refused must account for every supplied sample"
        );
        assert!(field.disposition(3).is_some());
        assert!(field.disposition(0).is_none());
    });
}

#[test]
fn an_entirely_edge_on_scan_refuses_rather_than_returning_an_empty_field() {
    with_default_cx(|cx| {
        let registration = identity_registration(cx);
        let probe = down_probe(1.0e-4, 1.0e-4);
        let normal = [(1.0f64 - 0.01 * 0.01).sqrt(), 0.0, 0.01];
        let samples: Vec<_> = (0..4)
            .map(|i| {
                SurfaceSample::new(
                    p3(f64::from(i), 0.0, 0.0),
                    normal,
                    p3(f64::from(i), 0.0, 0.0),
                )
                .expect("wall sample")
            })
            .collect();

        let error = DeviationField::extract(&samples, &registration, &probe, coverage(), cx)
            .expect_err("an all-grazing scan has no field");
        assert_eq!(error, FieldError::NoAdmittedSamples { grazing: 4 });
    });
}

// ---------------------------------------------------------------------------
// Thickness and gap.
// ---------------------------------------------------------------------------

#[test]
fn thickness_recovers_an_injected_bond_line_change() {
    with_default_cx(|cx| {
        let registration = identity_registration(cx);
        let probe = down_probe(0.0, 0.0);
        let nominal_thickness = 0.25;

        // Top face outward +z at z = 0.25; bottom face outward -z at z = 0.
        // Grow the bond line by 0.02 at the top and 0.01 at the bottom.
        let top_growth = 0.02;
        let bottom_growth = 0.01;
        let side_a = vec![
            SurfaceSample::new(
                p3(0.5, 0.5, nominal_thickness),
                [0.0, 0.0, 1.0],
                p3(0.5, 0.5, nominal_thickness + top_growth),
            )
            .expect("top sample"),
        ];
        let side_b = vec![
            SurfaceSample::new(
                p3(0.5, 0.5, 0.0),
                [0.0, 0.0, -1.0],
                p3(0.5, 0.5, -bottom_growth),
            )
            .expect("bottom sample"),
        ];

        let field = ThicknessField::extract(
            &side_a,
            &side_b,
            &[nominal_thickness],
            &registration,
            &probe,
            coverage(),
            cx,
        )
        .expect("thickness extracts");

        let point = field.points()[0];
        let expected = nominal_thickness + top_growth + bottom_growth;
        assert!(
            (point.thickness() - expected).abs() < 1e-12,
            "as-built thickness {} != {expected}",
            point.thickness()
        );
        assert!((point.change() - (top_growth + bottom_growth)).abs() < 1e-12);
    });
}

#[test]
fn thickness_is_first_order_insensitive_to_rigid_pose_error() {
    // The load-bearing correctness claim of the paired path: both faces ride
    // ONE pose, so a rigid motion cannot change the distance between them. The
    // summed-sensitivity composition must see that cancellation; adding the
    // two faces' sigmas in quadrature would manufacture uncertainty the rigid
    // motion cannot produce.
    with_default_cx(|cx| {
        let rotation = rodrigues([0.1, 0.9, -0.42], 0.25);
        let registration = registration_for(&rotation, [0.3, -0.2, 0.8], 2.0e-3, cx);
        let probe = ProbeModel::new([0.0, 0.0, -1.0], 0.05, 0.0, 0.0).expect("probe");

        let nominal_thickness = 0.4;
        let station = [1.0, 1.5];
        let side_a = vec![
            SurfaceSample::new(
                p3(station[0], station[1], nominal_thickness),
                [0.0, 0.0, 1.0],
                p3(station[0], station[1], nominal_thickness),
            )
            .expect("top"),
        ];
        let side_b = vec![
            SurfaceSample::new(
                p3(station[0], station[1], 0.0),
                [0.0, 0.0, -1.0],
                p3(station[0], station[1], 0.0),
            )
            .expect("bottom"),
        ];

        let thickness = ThicknessField::extract(
            &side_a,
            &side_b,
            &[nominal_thickness],
            &registration,
            &probe,
            coverage(),
            cx,
        )
        .expect("thickness extracts");

        // Each face on its own carries a real pose sigma.
        let single = DeviationField::extract(&side_a, &registration, &probe, coverage(), cx)
            .expect("single-face field");
        let face_sigma = single.points()[0].registration_sigma();
        assert!(face_sigma > 0.0, "one face must have a positive pose sigma");

        let paired = thickness.points()[0].registration_sigma();
        assert!(
            paired < face_sigma * 1e-3,
            "paired pose sigma {paired:e} must nearly cancel against the \
             single-face sigma {face_sigma:e}; quadrature would have given \
             about {:e}",
            face_sigma * std::f64::consts::SQRT_2
        );
    });
}

#[test]
fn a_thickness_station_is_refused_when_either_face_is_edge_on() {
    with_default_cx(|cx| {
        let registration = identity_registration(cx);
        let probe = down_probe(0.0, 0.0);
        let wall_normal = [(1.0f64 - 0.02 * 0.02).sqrt(), 0.0, -0.02];

        let side_a = vec![
            SurfaceSample::new(p3(0.0, 0.0, 0.2), [0.0, 0.0, 1.0], p3(0.0, 0.0, 0.2))
                .expect("good top"),
            SurfaceSample::new(p3(1.0, 0.0, 0.2), [0.0, 0.0, 1.0], p3(1.0, 0.0, 0.2))
                .expect("good top"),
        ];
        let side_b = vec![
            SurfaceSample::new(p3(0.0, 0.0, 0.0), [0.0, 0.0, -1.0], p3(0.0, 0.0, 0.0))
                .expect("good bottom"),
            // Second station's lower face is edge-on.
            SurfaceSample::new(p3(1.0, 0.0, 0.0), wall_normal, p3(1.0, 0.0, 0.0))
                .expect("edge-on bottom"),
        ];

        let field = ThicknessField::extract(
            &side_a,
            &side_b,
            &[0.2, 0.2],
            &registration,
            &probe,
            coverage(),
            cx,
        )
        .expect("thickness extracts");

        assert_eq!(field.points().len(), 1);
        assert_eq!(field.grazing(), &[1]);
    });
}

#[test]
fn thickness_refuses_mismatched_pairings() {
    with_default_cx(|cx| {
        let registration = identity_registration(cx);
        let probe = down_probe(0.0, 0.0);
        let a = vec![
            SurfaceSample::new(p3(0.0, 0.0, 0.2), [0.0, 0.0, 1.0], p3(0.0, 0.0, 0.2)).expect("top"),
        ];
        let b = Vec::new();
        let error = ThicknessField::extract(&a, &b, &[0.2], &registration, &probe, coverage(), cx)
            .expect_err("unpaired faces refuse");
        assert_eq!(
            error,
            FieldError::LengthMismatch {
                field: "side_b",
                expected: 1,
                found: 0,
            }
        );
    });
}

// ---------------------------------------------------------------------------
// Warpage form fit.
// ---------------------------------------------------------------------------

#[test]
fn form_fit_recovers_an_injected_quadratic_bow() {
    with_default_cx(|cx| {
        let registration = identity_registration(cx);
        let probe = down_probe(0.0, 0.0);
        let stations = grid_stations();
        // A saddle-free bow plus tilt and piston.
        let bow = |x: f64, y: f64| 0.003 + 0.001 * x - 0.0005 * y + 0.0008 * x * x + 0.0004 * y * y;
        let (samples, coordinates) = plate_samples(&stations, bow);

        let field = DeviationField::extract(&samples, &registration, &probe, coverage(), cx)
            .expect("field extracts");
        let fit = fit_form(&field, &coordinates, 2).expect("order-2 fit");

        assert_eq!(fit.order(), 2);
        assert_eq!(fit.point_count(), stations.len());
        assert!(
            fit.rms_residual() < 1e-12,
            "an exactly quadratic field must be fit exactly, residual {}",
            fit.rms_residual()
        );
        // The fitted surface must reproduce the injected field pointwise.
        for station in &stations {
            let expected = bow(station[0], station[1]);
            assert!(
                (fit.evaluate(*station) - expected).abs() < 1e-12,
                "fit does not reproduce the injected bow at {station:?}"
            );
        }
        // Warpage is the fitted surface's span over the stations.
        let mut low = f64::INFINITY;
        let mut high = f64::NEG_INFINITY;
        for station in &stations {
            let value = bow(station[0], station[1]);
            low = low.min(value);
            high = high.max(value);
        }
        assert!((fit.peak_to_valley() - (high - low)).abs() < 1e-12);
    });
}

#[test]
fn a_too_low_order_reports_what_it_missed_instead_of_hiding_it() {
    with_default_cx(|cx| {
        let registration = identity_registration(cx);
        let probe = down_probe(0.0, 0.0);
        let stations = grid_stations();
        let (samples, coordinates) =
            plate_samples(&stations, |x, y| 0.001 * x * x + 0.0006 * y * y);

        let field = DeviationField::extract(&samples, &registration, &probe, coverage(), cx)
            .expect("field extracts");
        let plane = fit_form(&field, &coordinates, 1).expect("order-1 fit");
        let quadratic = fit_form(&field, &coordinates, 2).expect("order-2 fit");

        assert!(
            plane.rms_residual() > 1e-4,
            "a plane cannot absorb a quadratic bow, yet residual is {}",
            plane.rms_residual()
        );
        assert!(quadratic.rms_residual() < 1e-12);
        assert!(plane.max_abs_residual() > plane.rms_residual());
    });
}

#[test]
fn form_fit_refuses_a_station_set_that_cannot_determine_the_order() {
    with_default_cx(|cx| {
        let registration = identity_registration(cx);
        let probe = down_probe(0.0, 0.0);
        // Every station on one line: y is constant, so no v-direction term is
        // determined. An order-2 fit is rank deficient.
        let stations: Vec<[f64; 2]> = (0..8).map(|i| [f64::from(i) * 0.25, 1.0]).collect();
        let (samples, coordinates) = plate_samples(&stations, |x, _| 0.001 * x);

        let field = DeviationField::extract(&samples, &registration, &probe, coverage(), cx)
            .expect("field extracts");
        let error =
            fit_form(&field, &coordinates, 2).expect_err("a collinear set cannot fix curvature");
        assert_eq!(error, FieldError::RankDeficientFit { order: 2 });

        // Order 1 must refuse too, and that is the STRICTER, correct answer:
        // the total-degree-1 basis carries a v-direction tilt term, and
        // stations on a single y = const line determine nothing about it.
        // Returning a fit here would mean returning an arbitrary member of the
        // solution family with a fabricated cross-line tilt.
        assert_eq!(
            fit_form(&field, &coordinates, 1)
                .expect_err("a cross-line tilt is undetermined on collinear stations"),
            FieldError::RankDeficientFit { order: 1 }
        );

        // Order 0 asks only for a piston, which one line of stations does
        // determine, so the refusal is specific rather than blanket.
        let piston = fit_form(&field, &coordinates, 0).expect("a piston is determined");
        assert_eq!(piston.coefficients().len(), 1);
    });
}

#[test]
fn form_fit_refuses_an_underdetermined_station_set() {
    with_default_cx(|cx| {
        let registration = identity_registration(cx);
        let probe = down_probe(0.0, 0.0);
        let stations = vec![[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]];
        let (samples, coordinates) = plate_samples(&stations, |x, y| 0.001 * (x + y));

        let field = DeviationField::extract(&samples, &registration, &probe, coverage(), cx)
            .expect("field extracts");
        let error = fit_form(&field, &coordinates, 2)
            .expect_err("three points cannot determine six coefficients");
        assert_eq!(
            error,
            FieldError::UnderdeterminedFit {
                points: 3,
                coefficients: 6,
            }
        );
    });
}

#[test]
fn form_fit_refuses_an_order_above_the_admitted_maximum() {
    with_default_cx(|cx| {
        let registration = identity_registration(cx);
        let probe = down_probe(0.0, 0.0);
        let stations = grid_stations();
        let (samples, coordinates) = plate_samples(&stations, |x, _| 0.001 * x);
        let field = DeviationField::extract(&samples, &registration, &probe, coverage(), cx)
            .expect("field extracts");
        let error = fit_form(&field, &coordinates, 9).expect_err("order is capped");
        assert_eq!(error, FieldError::FormOrderTooHigh { order: 9, max: 3 });
    });
}

// ---------------------------------------------------------------------------
// Roughness.
// ---------------------------------------------------------------------------

#[test]
fn square_wave_roughness_matches_its_closed_form() {
    // A zero-mean square wave of amplitude a has Ra = a, Rq = a, Rt = 2a, and
    // every segment spans 2a so Rz = 2a.
    let amplitude = 0.4;
    let heights: Vec<f64> = (0..64)
        .map(|i| if i % 2 == 0 { amplitude } else { -amplitude })
        .collect();
    let stats = profile_statistics(&heights, ProfileForm::Mean, 4).expect("statistics");

    assert!(
        (stats.ra() - amplitude).abs() < 1e-12,
        "Ra = {}",
        stats.ra()
    );
    assert!(
        (stats.rq() - amplitude).abs() < 1e-12,
        "Rq = {}",
        stats.rq()
    );
    assert!((stats.rt() - 2.0 * amplitude).abs() < 1e-12);
    assert!((stats.rz() - 2.0 * amplitude).abs() < 1e-12);
    assert_eq!(stats.segments(), 4);
    assert_eq!(stats.point_count(), 64);
}

#[test]
fn declared_form_removal_cancels_a_tilt() {
    // The same square wave, ridden onto a steep ramp. Line form removal must
    // return the untilted statistics; mean-only removal must not.
    let amplitude = 0.4;
    let flat: Vec<f64> = (0..64)
        .map(|i| if i % 2 == 0 { amplitude } else { -amplitude })
        .collect();
    let tilted: Vec<f64> = flat
        .iter()
        .enumerate()
        .map(|(i, value)| value + 0.05 * i as f64)
        .collect();

    let removed = profile_statistics(&tilted, ProfileForm::Line, 4).expect("line removal");
    let baseline = profile_statistics(&flat, ProfileForm::Line, 4).expect("baseline");
    assert!(
        (removed.ra() - baseline.ra()).abs() < 1e-12,
        "line removal must cancel the ramp: {} vs {}",
        removed.ra(),
        baseline.ra()
    );

    // Mean-only removal leaves the ramp behind. A centred ramp of half-span
    // `a` has mean |deviation| of `a/2`, here about 0.79 against the 0.4
    // square-wave baseline, so the inflation is roughly 2x — real, and the
    // point of declaring the form.
    let unremoved = profile_statistics(&tilted, ProfileForm::Mean, 4).expect("mean removal");
    assert!(
        unremoved.ra() > baseline.ra() * 1.5,
        "mean-only removal must NOT absorb a ramp: {} vs baseline {}",
        unremoved.ra(),
        baseline.ra()
    );
}

#[test]
fn a_parabolic_form_is_removed_only_at_the_declared_order() {
    let heights: Vec<f64> = (0..64)
        .map(|i| {
            let x = i as f64 / 63.0;
            0.5 * x * x + if i % 2 == 0 { 0.1 } else { -0.1 }
        })
        .collect();
    let parabola = profile_statistics(&heights, ProfileForm::Parabola, 4).expect("parabola");
    let line = profile_statistics(&heights, ProfileForm::Line, 4).expect("line");
    // The parabolic form comes out and the +/-0.1 ripple survives. Not to
    // machine precision: on a finite grid the alternating ripple is not
    // exactly orthogonal to the quadratic basis, so the least-squares form
    // absorbs a little of it (about 7e-4 relative here). Asserting 1e-12
    // would be asserting a false orthogonality.
    assert!(
        (parabola.ra() - 0.1).abs() < 1e-3,
        "the parabolic form must be removed leaving the ripple, Ra = {}",
        parabola.ra()
    );
    // Ra is the WRONG statistic to detect the leftover curvature, and that is
    // worth pinning rather than working around: when a smooth residual `f`
    // stays inside the ripple amplitude, the alternating pair contributes
    // `(|r + f| + |-r + f|) / 2 = r` exactly, so the mean-absolute measure
    // cancels the curvature and both orders report about 0.1. The total height
    // does not cancel, so Rt is what shows the under-removed form.
    assert!(
        (line.ra() - parabola.ra()).abs() < 1e-2,
        "Ra is expected to be blind here: {} vs {}",
        line.ra(),
        parabola.ra()
    );
    assert!(
        line.rt() > parabola.rt() * 1.25,
        "a straight mean line leaves curvature in the total height: {} vs {}",
        line.rt(),
        parabola.rt()
    );
}

#[test]
fn unfiltered_statistics_really_do_count_waviness() {
    // The documented no-claim is that these are NOT ISO 4287 values because no
    // lambda_c filter separates roughness from waviness. Prove the boundary is
    // real rather than merely asserted: adding a long-wavelength component
    // that a roughness filter would remove must inflate Ra here.
    let ripple: Vec<f64> = (0..256)
        .map(|i| if i % 2 == 0 { 0.05 } else { -0.05 })
        .collect();
    let with_waviness: Vec<f64> = ripple
        .iter()
        .enumerate()
        .map(|(i, value)| {
            let phase = std::f64::consts::TAU * i as f64 / 256.0;
            value + 0.5 * phase.sin()
        })
        .collect();

    let clean = profile_statistics(&ripple, ProfileForm::Line, 8).expect("clean");
    let wavy = profile_statistics(&with_waviness, ProfileForm::Line, 8).expect("wavy");
    assert!(
        wavy.ra() > clean.ra() * 3.0,
        "unfiltered Ra must absorb waviness ({} vs {}); if this ever stops \
         being true the crate is silently filtering and the no-claim is stale",
        wavy.ra(),
        clean.ra()
    );
}

#[test]
fn profile_statistics_refuse_degenerate_requests() {
    let heights = [0.1, -0.1, 0.2, -0.2];
    assert_eq!(
        profile_statistics(&heights, ProfileForm::Mean, 0).expect_err("zero segments"),
        FieldError::InvalidScalar {
            field: "segments",
            requirement: "at least one",
        }
    );
    assert_eq!(
        profile_statistics(&heights, ProfileForm::Mean, 9).expect_err("too few points"),
        FieldError::ProfileTooShort { points: 4, need: 9 }
    );
    let non_finite = [0.1, f64::NAN, 0.2, -0.2];
    assert!(matches!(
        profile_statistics(&non_finite, ProfileForm::Mean, 2),
        Err(FieldError::InvalidScalar { .. })
    ));
}

#[test]
fn rz_folds_a_trailing_partial_segment_into_the_last_full_one() {
    // 10 points into 3 segments: every point must be counted exactly once.
    let heights: Vec<f64> = (0..10).map(|i| if i == 9 { 5.0 } else { 0.0 }).collect();
    let stats = profile_statistics(&heights, ProfileForm::Mean, 3).expect("statistics");
    // The tall final point sits in the folded remainder, so it must influence
    // Rz; if the remainder were dropped, the last segment would be flat.
    assert!(
        stats.rz() > 0.0,
        "the trailing remainder must be counted, not dropped"
    );
    assert_eq!(stats.point_count(), 10);
}

// ---------------------------------------------------------------------------
// G5 determinism, identity, cancellation, admission.
// ---------------------------------------------------------------------------

#[test]
fn field_extraction_is_deterministic_across_independent_runs() {
    let build = || {
        with_default_cx(|cx| {
            let registration = identity_registration(cx);
            let probe = down_probe(1.0e-4, 5.0e-4);
            let stations = grid_stations();
            let (samples, _) = plate_samples(&stations, |x, y| 0.002 * x - 0.001 * y);
            DeviationField::extract(&samples, &registration, &probe, coverage(), cx)
                .expect("field extracts")
        })
    };
    let first = build();
    let second = build();
    assert_eq!(first.identity(), second.identity());
    assert_eq!(first, second);
}

#[test]
fn the_identity_moves_when_the_measured_content_moves() {
    with_default_cx(|cx| {
        let registration = identity_registration(cx);
        let probe = down_probe(1.0e-4, 5.0e-4);
        let stations = grid_stations();

        let (baseline_samples, _) = plate_samples(&stations, |x, _| 0.002 * x);
        let (perturbed_samples, _) = plate_samples(&stations, |x, y| {
            0.002 * x + if y == 0.0 { 1e-9 } else { 0.0 }
        });

        let baseline =
            DeviationField::extract(&baseline_samples, &registration, &probe, coverage(), cx)
                .expect("baseline");
        let perturbed =
            DeviationField::extract(&perturbed_samples, &registration, &probe, coverage(), cx)
                .expect("perturbed");

        assert_ne!(
            baseline.identity(),
            perturbed.identity(),
            "a one-nanometre content change must move the identity"
        );
    });
}

#[test]
fn the_identity_binds_the_refusal_set() {
    // Two scans with identical ADMITTED points but different refusals must not
    // share an identity: what was refused is part of what the field claims.
    with_default_cx(|cx| {
        let registration = identity_registration(cx);
        let probe = down_probe(0.0, 0.0);
        let good = SurfaceSample::new(p3(0.0, 0.0, 0.0), [0.0, 0.0, 1.0], p3(0.0, 0.0, 0.01))
            .expect("good");
        let wall_normal = [(1.0f64 - 0.01 * 0.01).sqrt(), 0.0, 0.01];
        let wall =
            SurfaceSample::new(p3(1.0, 0.0, 0.0), wall_normal, p3(1.0, 0.0, 0.0)).expect("wall");

        let without = DeviationField::extract(&[good], &registration, &probe, coverage(), cx)
            .expect("without refusal");
        let with = DeviationField::extract(&[good, wall], &registration, &probe, coverage(), cx)
            .expect("with refusal");

        assert_eq!(without.points(), with.points());
        assert_ne!(
            without.identity(),
            with.identity(),
            "the refusal set must be bound into the identity"
        );
    });
}

#[test]
fn cancellation_publishes_no_partial_field() {
    let registration = with_default_cx(identity_registration);
    with_cx(true, |cx| {
        let probe = down_probe(0.0, 0.0);
        let stations = grid_stations();
        let (samples, _) = plate_samples(&stations, |x, _| 0.001 * x);
        let error = DeviationField::extract(&samples, &registration, &probe, coverage(), cx)
            .expect_err("a cancelled scan publishes nothing");
        assert!(matches!(error, FieldError::Cancelled { .. }));
    });
}

#[test]
fn an_empty_sample_set_refuses() {
    with_default_cx(|cx| {
        let registration = identity_registration(cx);
        let probe = down_probe(0.0, 0.0);
        assert_eq!(
            DeviationField::extract(&[], &registration, &probe, coverage(), cx)
                .expect_err("no samples"),
            FieldError::NoSamples
        );
    });
}

#[test]
fn a_non_unit_normal_is_refused_at_declaration() {
    assert!(
        SurfaceSample::new(p3(0.0, 0.0, 0.0), [0.0, 0.0, 0.5], p3(0.0, 0.0, 0.0)).is_err(),
        "a half-length normal would silently halve every deviation"
    );
    assert!(
        SurfaceSample::new(p3(0.0, 0.0, 0.0), [0.0, 0.0, f64::NAN], p3(0.0, 0.0, 0.0)).is_err()
    );
    assert!(SurfaceSample::new(p3(0.0, 0.0, 0.0), [0.0, 0.0, 1.0], p3(0.0, 0.0, 0.0)).is_ok());
}

#[test]
fn probe_admission_rejects_degenerate_declarations() {
    assert!(
        ProbeModel::new([0.0, 0.0, 2.0], 0.1, 0.0, 0.0).is_err(),
        "non-unit direction"
    );
    assert!(
        ProbeModel::new([0.0, 0.0, 1.0], 0.0, 0.0, 0.0).is_err(),
        "zero floor"
    );
    assert!(
        ProbeModel::new([0.0, 0.0, 1.0], 1.5, 0.0, 0.0).is_err(),
        "floor above one"
    );
    assert!(
        ProbeModel::new([0.0, 0.0, 1.0], 0.1, -1.0, 0.0).is_err(),
        "negative sigma"
    );
    assert!(ProbeModel::new([0.0, 0.0, 1.0], 0.1, 0.0, 0.0).is_ok());
}

#[test]
fn a_field_cites_the_registration_it_was_extracted_against() {
    with_default_cx(|cx| {
        let registration = identity_registration(cx);
        let probe = down_probe(0.0, 0.0);
        let sample = SurfaceSample::new(p3(0.0, 0.0, 0.0), [0.0, 0.0, 1.0], p3(0.0, 0.0, 0.01))
            .expect("sample");
        let field = DeviationField::extract(&[sample], &registration, &probe, coverage(), cx)
            .expect("field extracts");
        assert_eq!(field.registration_ref(), registration.model_identity());
        assert_ne!(
            field.identity(),
            registration.model_identity(),
            "the field identity is its own, not a copy of the citation"
        );
    });
}

#[test]
fn resolution_screens_compare_the_signal_against_its_own_half_width() {
    with_default_cx(|cx| {
        let registration = identity_registration(cx);
        // A deliberately noisy probe so the small deviation is unresolved.
        let probe = down_probe(0.0, 0.05);
        let small = SurfaceSample::new(p3(0.0, 0.0, 0.0), [0.0, 0.0, 1.0], p3(0.0, 0.0, 1e-6))
            .expect("small");
        let large = SurfaceSample::new(p3(1.0, 0.0, 0.0), [0.0, 0.0, 1.0], p3(1.0, 0.0, 5.0))
            .expect("large");
        let field = DeviationField::extract(&[small, large], &registration, &probe, coverage(), cx)
            .expect("field extracts");
        assert!(
            !field.points()[0].resolved(),
            "1 um under a 50 mm sigma is not resolved"
        );
        assert!(
            field.points()[1].resolved(),
            "5 m over a 50 mm sigma is resolved"
        );
    });
}
