use super::*;
use crate::dielectric::{BeerLambertAbsorption, CauchyIor, GlassProvenance};

fn glass(n: f64) -> DielectricGlass {
    DielectricGlass::new(
        CauchyIor::try_constant(n).unwrap(),
        BeerLambertAbsorption::CLEAR,
        GlassProvenance::Custom,
    )
}
fn plane() -> NormalDifferential {
    NormalDifferential {
        normal: Vec3::new(0.0, 0.0, -1.0),
        derivative_per_m: [[0.0; 3]; 3],
    }
}
fn close(actual: Vec3, expected: Vec3, tolerance: f64) {
    assert!(
        vec_norm(sub(actual, expected)) < tolerance,
        "{actual:?} != {expected:?}"
    );
}

#[test]
fn g1_equal_local_indices_still_refract_neighboring_rays_when_spatial_gradients_differ() {
    let material = glass(1.5);
    let incident = HomogeneousIndexField::new(&material);
    let ray = Ray {
        origin: Point3::new(0.0, 0.0, 0.0),
        dir: Vec3::new(0.6, 0.0, 0.8),
    };
    let bundle = propagate_bundle(
        &incident,
        ray,
        Vec3::new(0.8, 0.0, -0.6),
        550.0,
        1.0,
        BundleOptions::default(),
    )
    .unwrap();
    let target =
        GradientIndexSample::try_new(1.5, Vec3::new(0.2, 0.1, 0.0), 0.0, Some(1.5)).unwrap();
    let seed = bundle
        .smooth_interface(
            plane(),
            target,
            DielectricEvent::Transmission,
            InterfaceOptions::default(),
        )
        .unwrap();
    close(seed.ray().dir, ray.dir, 1e-14);
    let columns = seed.angular_derivatives();
    // At z=.8 the surface-arrival derivatives are (1.25,0,0) and
    // (0,1,0). Differentiating n_i*t_parallel=n_t*t'_parallel gives
    // dt'_x=.8-(.6/1.5)*.2*1.25=.7 and dt'_x(second)=-.04.
    // Unit-direction tangency then fixes the z components independently.
    close(
        columns[0].position_m_per_rad,
        Vec3::new(1.25, 0.0, 0.0),
        1e-14,
    );
    close(
        columns[0].direction_per_rad,
        Vec3::new(0.7, 0.0, -0.525),
        1e-14,
    );
    close(
        columns[1].direction_per_rad,
        Vec3::new(-0.04, 1.0, 0.03),
        1e-14,
    );
}

/// Exact affine-index trajectory to z=.8, using conserved transverse optical
/// momenta and hyperbolic functions. No production propagation or differential
/// is used to determine this reference.
fn affine_reference(direction: Vec3) -> (Point3, Vec3, f64, f64) {
    let n0 = 1.5;
    let gradient = 0.2;
    let z = 0.8;
    let py = n0 * direction.y;
    let pz = n0 * direction.z;
    let transverse = py.hypot(pz);
    let u0 = (n0 * direction.x / transverse).asinh();
    let u = u0 + gradient * z / pz;
    let n = transverse * u.cosh();
    let point = Point3::new((n - n0) / gradient, py * z / pz, z);
    let tangent = Vec3::new(u.tanh(), py / n, pz / n);
    let distance = transverse * (u.sinh() - u0.sinh()) / gradient;
    (point, tangent, n, distance)
}

#[test]
fn g2_both_index_gradients_and_curved_arrival_match_independent_affine_snell_rays() {
    let material = glass(1.5);
    let field = AffineIndexField::try_new(
        &material,
        Point3::new(0.0, 0.0, 0.0),
        Vec3::new(0.2, 0.0, 0.0),
        Aabb::new(Point3::new(-2.0, -2.0, -1.0), Point3::new(3.0, 3.0, 3.0)),
    )
    .unwrap();
    let direction = Vec3::new(0.6, 0.0, 0.8);
    let (point, _, _, distance) = affine_reference(direction);
    let mut options = BundleOptions::default();
    options.integration.position_tolerance_m = 1e-13;
    options.integration.direction_tolerance = 1e-13;
    options.integration.relative_tolerance = 1e-13;
    options.position_derivative_tolerance_m = 1e-13;
    options.direction_derivative_tolerance = 1e-13;
    let bundle = propagate_bundle(
        &field,
        Ray {
            origin: Point3::new(0.0, 0.0, 0.0),
            dir: direction,
        },
        Vec3::new(0.8, 0.0, -0.6),
        550.0,
        distance,
        options,
    )
    .unwrap();
    close(
        bundle.propagation().ray.origin.delta_from(point),
        Vec3::new(0.0, 0.0, 0.0),
        2e-12,
    );
    let target_gradient = Vec3::new(0.1, 0.08, 0.0);
    let target = GradientIndexSample::try_new(1.6, target_gradient, 0.0, Some(1.6)).unwrap();
    let mapped = bundle
        .smooth_interface(
            plane(),
            target,
            DielectricEvent::Transmission,
            InterfaceOptions::default(),
        )
        .unwrap();
    let reference = |direction| {
        let (position, tangent, incident, _) = affine_reference(direction);
        let target = 1.6 + target_gradient.dot(position.delta_from(point));
        let x = incident * tangent.x / target;
        let y = incident * tangent.y / target;
        (position, Vec3::new(x, y, (1.0 - x * x - y * y).sqrt()))
    };
    let step = 1e-5;
    for (axis, actual) in bundle
        .launch_basis()
        .into_iter()
        .zip(mapped.angular_derivatives())
    {
        let plus = reference(unit(sub(direction, axis.scale(-step))).unwrap());
        let minus = reference(unit(sub(direction, axis.scale(step))).unwrap());
        close(
            actual.position_m_per_rad,
            plus.0.delta_from(minus.0).scale(0.5 / step),
            2e-8,
        );
        close(
            actual.direction_per_rad,
            sub(plus.1, minus.1).scale(0.5 / step),
            2e-8,
        );
    }
}
