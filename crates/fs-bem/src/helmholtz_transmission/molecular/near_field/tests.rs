//! G0/G1/G3 fixed-trace integration, real-owner parity and refusal coverage.
use super::super::tests::{plane, state};
use super::*;
use crate::helmholtz::bc::{IncidentSource, Monopole, PlaneWave};
use crate::helmholtz::transmission::FluidInterface;

fn options() -> MolecularNearFieldOptions {
    MolecularNearFieldOptions {
        quadrature: crate::near_field::Options {
            relative_tolerance: 1e-8,
            ..crate::near_field::Options::default()
        },
        maximum_row_bytes: 4_000_000,
    }
}

fn close(a: C64, b: C64, tolerance: f64) {
    assert!(
        (a - b).abs() < tolerance * (1.0 + b.abs()),
        "{a:?} != {b:?}"
    );
}

#[test]
fn empty_reservoir_adaptive_fields_recover_original_real_owner_on_both_sides() {
    let gate = CancelGate::new();
    let surface = SpherePanels::icosphere(1.0, 0).unwrap();
    let outside = state(293.15, 101325.0, None);
    let inside = state(383.15, 101325.0, None);
    let omega = 200.0;
    let molecular = MolecularFluidInterface::new(
        &surface,
        &outside,
        &inside,
        omega,
        Options::default(),
        &gate,
    )
    .unwrap();
    let original = FluidInterface::new(
        &surface,
        outside.gas(),
        inside.gas(),
        omega,
        Options::default(),
        &gate,
    )
    .unwrap();
    let actual = molecular.scatter_field(plane(), 0.01, 1000, &gate).unwrap();
    let expected = original.scatter_field(plane(), 0.01, &gate).unwrap();
    let exterior = [[0.2, 0.1, 1.08], [1.4, -0.3, 0.2]];
    let adaptive = actual
        .scattered_pressure_velocity_adaptive_at(&exterior, 0.01, 1000, options(), &gate)
        .unwrap();
    let reference = Geometry::new(&surface, &exterior, 0.01)
        .unwrap()
        .prepare_velocity(
            omega / outside.gas().sound_speed,
            Medium::from(outside.gas()),
            options().quadrature,
        )
        .unwrap()
        .evaluate_velocity(expected.solution().scattered())
        .unwrap();
    for i in 0..exterior.len() {
        close(
            adaptive.scalar.pressure[i],
            reference.scalar.pressure[i],
            3e-11,
        );
        for axis in 0..3 {
            close(
                adaptive.particle_velocity_m_s[i][axis],
                reference.particle_velocity_m_s[i][axis],
                3e-11,
            );
        }
    }
    let total = actual
        .total_pressure_velocity_adaptive_at(&exterior, 0.01, 1000, options(), &gate)
        .unwrap();
    let total_reference = expected
        .total_pressure_at(&exterior, 0.01, options().quadrature)
        .unwrap();
    for (i, (&a, &b)) in total
        .scalar
        .pressure
        .iter()
        .zip(&total_reference.pressure)
        .enumerate()
    {
        // The original scalar-only tree need not refine as far as the paired
        // pressure/gradient tree. Compare using both declared discrepancies.
        let quadrature = total.scalar.quadrature_error_estimate_pa[i]
            + total_reference.quadrature_error_estimate_pa[i];
        assert!((a - b).abs() < 4.0 * quadrature + 3e-11 * (1.0 + b.abs()));
    }
    let interior = [[0.2, 0.1, 0.3], [-0.25, 0.1, -0.4]];
    let transmitted = actual
        .interior_pressure_velocity_adaptive_at(&interior, 0.01, 1000, options(), &gate)
        .unwrap();
    let transmitted_reference = expected
        .solution()
        .interior_pressure_at(&interior, 0.01, options().quadrature)
        .unwrap();
    for (i, (&a, &b)) in transmitted
        .scalar
        .pressure
        .iter()
        .zip(&transmitted_reference.pressure)
        .enumerate()
    {
        let quadrature = transmitted.scalar.quadrature_error_estimate_pa[i]
            + transmitted_reference.quadrature_error_estimate_pa[i];
        assert!((a - b).abs() < 4.0 * quadrature + 3e-11 * (1.0 + b.abs()));
    }
}

// Independent composite midpoint integration of the ORIGINAL flat triangles.
// No production Green, adaptive tree, or prepared-row evaluator is called.
// Uniform 4/5/6 subdivision ladders resolve the fixed physical p/v traces.
#[allow(clippy::too_many_arguments)]
fn reference_integral(
    surface: &SpherePanels,
    pressure: &[C64],
    velocity: &[C64],
    k: C64,
    rho_omega: f64,
    point: [f64; 3],
    depth: usize,
    sign: f64,
) -> [C64; 4] {
    let mut result = [C64::ZERO; 4];
    for (panel, &triangle) in surface.triangles().unwrap().iter().enumerate() {
        integrate_triangle(
            triangle,
            surface.normals()[panel],
            pressure[panel],
            velocity[panel],
            k,
            rho_omega,
            point,
            depth,
            &mut result,
        );
    }
    result[0] = result[0].scale(sign);
    for u in &mut result[1..] {
        *u = *u * C64::new(0.0, -sign / rho_omega);
    }
    result
}

#[allow(clippy::too_many_arguments)]
fn integrate_triangle(
    triangle: [[f64; 3]; 3],
    normal: [f64; 3],
    p: C64,
    v: C64,
    k: C64,
    rho_omega: f64,
    point: [f64; 3],
    depth: usize,
    sum: &mut [C64; 4],
) {
    let [a, b, c] = triangle;
    if depth != 0 {
        let midpoint = |a: [f64; 3], b: [f64; 3]| std::array::from_fn(|i| 0.5 * (a[i] + b[i]));
        let (ab, ac, bc) = (midpoint(a, b), midpoint(a, c), midpoint(b, c));
        for subtriangle in [[a, ab, ac], [ab, b, bc], [ac, bc, c], [ab, bc, ac]] {
            integrate_triangle(
                subtriangle,
                normal,
                p,
                v,
                k,
                rho_omega,
                point,
                depth - 1,
                sum,
            );
        }
        return;
    }
    let subtract = |a: [f64; 3], b: [f64; 3]| std::array::from_fn::<_, 3, _>(|i| a[i] - b[i]);
    let (ab, ac) = (subtract(b, a), subtract(c, a));
    let cross = [
        ab[1] * ac[2] - ab[2] * ac[1],
        ab[2] * ac[0] - ab[0] * ac[2],
        ab[0] * ac[1] - ab[1] * ac[0],
    ];
    let area = 0.5 * cross.iter().map(|x| x * x).sum::<f64>().sqrt();
    let centre = std::array::from_fn(|i| (a[i] + b[i] + c[i]) / 3.0);
    let delta = subtract(point, centre);
    let r = delta.iter().map(|x| x * x).sum::<f64>().sqrt();
    let unit = delta.map(|x| x / r);
    let projection: f64 = normal.iter().zip(unit).map(|(a, b)| a * b).sum();
    let phase = k.re * r;
    let amplitude = fs_math::det::exp(-k.im * r) / (4.0 * core::f64::consts::PI * r);
    let g = C64::new(fs_math::det::cos(phase), fs_math::det::sin(phase)).scale(amplitude);
    let first = C64::new(-k.im - 1.0 / r, k.re);
    let radial = first * first + C64::from_re(1.0 / (r * r)) - first.scale(1.0 / r);
    let q = v * C64::new(0.0, rho_omega);
    sum[0] = sum[0] + (-g * first.scale(projection) * p - g * q).scale(area);
    for axis in 0..3 {
        let ds = g * first.scale(unit[axis]);
        let dd = -g * (radial.scale(projection * unit[axis]) + first.scale(normal[axis] / r));
        sum[axis + 1] = sum[axis + 1] + (dd * p - ds * q).scale(area);
    }
}

fn relative(reference: &[C64], actual: &[C64]) -> f64 {
    let numerator: f64 = reference
        .iter()
        .zip(actual)
        .map(|(&a, &b)| (a - b).norm_sq())
        .sum();
    let denominator: f64 = reference.iter().map(|a| a.norm_sq()).sum();
    (numerator / denominator.max(1e-40)).sqrt()
}

#[test]
fn adaptive_pair_resolves_nearby_flat_triangle_integrals_with_physical_interior_momentum() {
    let gate = CancelGate::new();
    let surface = SpherePanels::icosphere(1.0, 0).unwrap();
    let outside = state(293.15, 101325.0, Some(0.5));
    let inside = state(443.15, 101325.0, Some(5.0));
    let omega = 200.0;
    let owner = MolecularFluidInterface::new(
        &surface,
        &outside,
        &inside,
        omega,
        Options::default(),
        &gate,
    )
    .unwrap();
    let sources = [
        IncidentSource::PlaneWave(PlaneWave {
            direction: [0.0, 0.0, 1.0],
            amplitude_pa: C64::new(1.0, 0.3),
        }),
        IncidentSource::Monopole(Monopole {
            position_m: [2.0, 1.5, 0.7],
            volume_velocity_m3_s: C64::new(0.003, -0.001),
        }),
    ];
    let field = owner
        .scatter_field(IncidentField::Sources(&sources), 0.01, 1000, &gate)
        .unwrap();
    for interior in [false, true] {
        let sign = if interior { -1.0 } else { 1.0 };
        let point = std::array::from_fn(|i| {
            surface.centroids()[0][i] + sign * 0.04 * surface.normals()[0][i]
        });
        let actual = if interior {
            field.interior_pressure_velocity_adaptive_at(&[point], 0.01, 1000, options(), &gate)
        } else {
            field.scattered_pressure_velocity_adaptive_at(&[point], 0.01, 1000, options(), &gate)
        }
        .unwrap();
        let values = [
            actual.scalar.pressure[0],
            actual.particle_velocity_m_s[0][0],
            actual.particle_velocity_m_s[0][1],
            actual.particle_velocity_m_s[0][2],
        ];
        let (source, p, v) = if interior {
            (
                &inside,
                field.total_pressure_pa(),
                field.total_normal_velocity_m_s(),
            )
        } else {
            (
                &outside,
                field.scattered_pressure_pa(),
                field.scattered_normal_velocity_m_s(),
            )
        };
        let mut last = [f64::INFINITY; 2];
        for depth in 4..=6 {
            let expected = reference_integral(
                &surface,
                p,
                v,
                source.relaxing_wavenumber(omega).unwrap(),
                source.gas().density * omega,
                point,
                depth,
                sign,
            );
            let error = [
                relative(&expected[..1], &values[..1]),
                relative(&expected[1..], &values[1..]),
            ];
            for axis in 0..2 {
                assert!(
                    error[axis] < 0.7 * last[axis],
                    "interior={interior} depth={depth} {error:?} after {last:?}"
                );
            }
            last = error;
        }
        assert!(last[0] < 0.004 && last[1] < 0.004, "{last:?}");
        assert!(actual.scalar.quadrature_error_estimate_pa[0] > 0.0);
        assert!(
            actual.quadrature_error_estimate_m_s[0]
                .iter()
                .all(|e| e.is_finite() && *e >= 0.0)
        );
        if !interior {
            let centroid = field
                .scattered_pressure_velocity_at(&[point], 0.01, 1000, &gate)
                .unwrap()[0];
            assert!(
                relative(&values[1..], &centroid.velocity_m_s) > 0.01,
                "the near-panel fixture must distinguish centroid observation"
            );
        }
    }
}

#[test]
fn adaptive_admission_and_exact_transparency_preserve_source_and_cache() {
    let gate = CancelGate::new();
    let surface = SpherePanels::icosphere(1.0, 0).unwrap();
    let outside = state(293.15, 101325.0, Some(0.5));
    let inside = state(383.15, 101325.0, Some(5.0));
    let owner = MolecularFluidInterface::new(
        &surface,
        &outside,
        &inside,
        200.0,
        Options::default(),
        &gate,
    )
    .unwrap();
    let field = owner.scatter_field(plane(), 0.01, 1000, &gate).unwrap();
    let point = std::array::from_fn(|i| surface.centroids()[0][i] + 0.03 * surface.normals()[0][i]);
    for refused in [
        MolecularNearFieldOptions {
            maximum_row_bytes: 1,
            ..options()
        },
        MolecularNearFieldOptions {
            quadrature: crate::near_field::Options {
                maximum_kernel_evaluations: 80,
                ..options().quadrature
            },
            ..options()
        },
        MolecularNearFieldOptions {
            quadrature: crate::near_field::Options {
                maximum_depth: 0,
                ..options().quadrature
            },
            ..options()
        },
    ] {
        assert!(
            field
                .total_pressure_velocity_adaptive_at(&[point], 0.01, 1000, refused, &gate)
                .is_err()
        );
    }
    let first = field
        .total_pressure_velocity_adaptive_at(&[point], 0.01, 1000, options(), &gate)
        .unwrap();
    let replay = field
        .total_pressure_velocity_adaptive_at(&[point], 0.01, 1000, options(), &gate)
        .unwrap();
    assert_eq!(first.scalar.pressure, replay.scalar.pressure);
    assert_eq!(first.particle_velocity_m_s, replay.particle_velocity_m_s);
    let cancelled = CancelGate::new();
    cancelled.request();
    assert!(matches!(
        field.total_pressure_velocity_adaptive_at(&[point], 0.01, 1000, options(), &cancelled),
        Err(HelmholtzError::Cancelled)
    ));

    let transparent = MolecularFluidInterface::new(
        &surface,
        &outside,
        &outside,
        200.0,
        Options::default(),
        &gate,
    )
    .unwrap();
    let incident = plane();
    let exact = transparent
        .scatter_field(incident, 0.01, 20, &gate)
        .unwrap();
    for (point, interior) in [([0.1, 0.2, -0.1], true), ([0.1, 0.2, 1.1], false)] {
        let actual = if interior {
            exact.interior_pressure_velocity_adaptive_at(&[point], 0.01, 1000, options(), &gate)
        } else {
            exact.total_pressure_velocity_adaptive_at(&[point], 0.01, 1000, options(), &gate)
        }
        .unwrap();
        let expected = incident
            .sample_wave(
                transparent.k_o,
                outside.gas().density * 200.0,
                None,
                point,
                &mut || false,
            )
            .unwrap();
        assert_eq!(actual.scalar.pressure, [expected.pressure_pa]);
        assert_eq!(actual.particle_velocity_m_s, [expected.velocity_m_s]);
        assert_eq!(actual.scalar.quadrature_error_estimate_pa, [0.0]);
        assert_eq!(actual.quadrature_error_estimate_m_s, [[0.0; 3]]);
    }
    assert!(
        exact
            .total_pressure_velocity_adaptive_at(
                &[[0.1, 0.2, 1.1]; 21],
                0.01,
                1000,
                options(),
                &gate
            )
            .is_err()
    );
    assert!(
        exact
            .interior_pressure_velocity_adaptive_at(
                &[[0.0, 0.0, 1.1]],
                0.01,
                1000,
                options(),
                &gate
            )
            .is_err()
    );
}
