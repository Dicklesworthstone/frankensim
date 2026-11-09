use super::*;
use crate::conductor::microsurface::SmithEvent;
use core::f64::consts::PI;
use fs_rand::philox::philox4x32_10;

fn draw(path: u32, event: u32, domain: u32) -> [f64; 3] {
    let b = philox4x32_10([path, event, domain, 0], [19, 41]);
    core::array::from_fn(|i| (f64::from(b[i]) + 0.5) / 4_294_967_296.0)
}
fn near(a: f64, b: f64, tolerance: f64) {
    assert!((a - b).abs() <= tolerance, "{a:.17e} != {b:.17e}");
}
fn lambda(v: Vec3, ax: f64, ay: f64) -> f64 {
    0.5 * ((1.0 + (ax * ax * v.x * v.x + ay * ay * v.y * v.y) / (v.z * v.z)).sqrt() - 1.0)
}
fn reference_first(wi: Vec3, wo: Vec3, ax: f64, ay: f64) -> f64 {
    let h = unit(Vec3::new(wi.x + wo.x, wi.y + wo.y, wi.z + wo.z)).unwrap();
    let d = 1.0 / (PI * ax * ay * ((h.x / ax).powi(2) + (h.y / ay).powi(2) + h.z * h.z).powi(2));
    d / (4.0 * wi.z * wo.z * (1.0 + lambda(wi, ax, ay) + lambda(wo, ax, ay)))
}
#[test]
fn g1_first_connection_is_analytic_correlated_smith_and_reciprocal() {
    for (ax, ay) in [(1.0, 1.0), (0.07, 0.8), (1e-4, 1e-4)] {
        let model = SmithGgx::try_new(ax, ay).unwrap();
        for wi in [
            Vec3::new(0.0, 0.0, 1.0),
            unit(Vec3::new(0.8, 0.3, 0.2)).unwrap(),
        ] {
            for wo in [
                Vec3::new(0.0, 0.0, 1.0),
                unit(Vec3::new(-0.1, 0.6, 0.4)).unwrap(),
            ] {
                let a = model
                    .single_scattering_connection(wi.scale(-1.0), wo)
                    .unwrap()
                    .unwrap();
                let b = model
                    .single_scattering_connection(wo.scale(-1.0), wi)
                    .unwrap()
                    .unwrap();
                let r = reference_first(wi, wo, ax, ay);
                near(a.factor_per_sr, r, 3e-12 * r.max(1.0));
                near(
                    a.factor_per_sr,
                    b.factor_per_sr,
                    3e-13 * a.factor_per_sr.max(1.0),
                );
            }
        }
    }
}
#[test]
fn g1_depth_connection_integrates_to_the_first_order_law() {
    let model = SmithGgx::try_new(0.3, 0.9).unwrap();
    let incoming = unit(Vec3::new(0.7, 0.2, -0.4)).unwrap();
    let outgoing = unit(Vec3::new(-0.4, 0.3, 0.6)).unwrap();
    let expected = model
        .single_scattering_connection(incoming, outgoing)
        .unwrap()
        .unwrap()
        .factor_per_sr;
    // Deterministic midpoint integration in the first free-flight uniform.
    let n = 65_536;
    let mut sum = 0.0;
    for i in 0..n {
        let step = model
            .sample(incoming, 0.0, [(i as f64 + 0.5) / n as f64, 0.3, 0.7])
            .unwrap();
        let SmithEvent::Reflection { depth, .. } = step.event else {
            panic!("incident flight escaped")
        };
        sum += model
            .connection_at_depth(incoming, depth, outgoing)
            .unwrap()
            .unwrap()
            .factor_per_sr;
    }
    near(sum / n as f64, expected, 4e-6 * expected);
}

// Scalar constant facet reflectance is an independent statistical fixture for
// the geometric connection, including all internal returns. No optical owner
// or fitted furnace compensation enters this test.
fn estimate(model: SmithGgx, wi: Vec3, wo: Vec3, rho: f64, path: u32, domain: u32) -> f64 {
    let mut sum = rho
        * model
            .single_scattering_connection(wi.scale(-1.0), wo)
            .unwrap()
            .unwrap()
            .factor_per_sr;
    let mut direction = wi.scale(-1.0);
    let mut depth = 0.0;
    let mut power = 1.0;
    for event in 0..128 {
        match model
            .sample(direction, depth, draw(path, event, domain))
            .unwrap()
            .event
        {
            SmithEvent::Escape => return sum,
            SmithEvent::Reflection {
                depth: next,
                direction: next_direction,
                ..
            } => {
                if event > 0 {
                    if let Some(c) = model.connection_at_depth(direction, next, wo).unwrap() {
                        sum += power * rho * c.factor_per_sr;
                    }
                }
                direction = next_direction;
                depth = next;
                power *= rho;
            }
        }
    }
    panic!("unresolved path");
}
fn moments(values: impl Iterator<Item = f64>, n: usize) -> (f64, f64) {
    let mut sum = 0.0;
    let mut squares = 0.0;
    for x in values {
        sum += x;
        squares += x * x;
    }
    let mean = sum / n as f64;
    (
        mean,
        ((squares / n as f64 - mean * mean).max(0.0) / n as f64).sqrt(),
    )
}
#[test]
fn g1_full_connection_white_furnace_and_absorption_match_sampled_escape() {
    let n = 65_536;
    let model = SmithGgx::try_new(1.0, 1.0).unwrap();
    let wi = Vec3::new(0.0, 0.0, 1.0);
    for rho in [1.0, 0.35] {
        let (response, error) = moments(
            (0..n).map(|p| {
                let u = draw(p as u32, 0, 91);
                let (s, c) = det::sin_cos(2.0 * PI * u[1]);
                let wo = Vec3::new(u[0].sqrt() * c, u[0].sqrt() * s, (1.0 - u[0]).sqrt());
                PI * estimate(model, wi, wo, rho, p as u32, 7)
            }),
            n,
        );
        let (escaped, escape_error) = moments(
            (0..n).map(|p| {
                let mut w = wi.scale(-1.0);
                let mut d = 0.0;
                let mut power = 1.0;
                for e in 0..128 {
                    match model.sample(w, d, draw(p as u32, e, 33)).unwrap().event {
                        SmithEvent::Escape => return power,
                        SmithEvent::Reflection {
                            depth, direction, ..
                        } => {
                            d = depth;
                            w = direction;
                            power *= rho;
                        }
                    }
                }
                panic!("unresolved path")
            }),
            n,
        );
        let tolerance = 7.0 * (error * error + escape_error * escape_error).sqrt() + 0.002;
        near(response, escaped, tolerance);
        if rho == 1.0 {
            near(escaped, 1.0, 0.0);
        }
    }
}
#[test]
fn g3_full_connection_reciprocity_holds_with_anisotropic_internal_returns() {
    let model = SmithGgx::try_new(0.25, 0.85).unwrap();
    let wi = unit(Vec3::new(0.8, 0.2, 0.4)).unwrap();
    let wo = unit(Vec3::new(-0.3, 0.6, 0.5)).unwrap();
    let n = 65_536;
    let (a, sa) = moments(
        (0..n).map(|p| estimate(model, wi, wo, 0.8, p as u32, 12)),
        n,
    );
    let (b, sb) = moments(
        (0..n).map(|p| estimate(model, wo, wi, 0.8, p as u32, 24)),
        n,
    );
    near(a, b, 7.0 * (sa * sa + sb * sb).sqrt() + 0.0005);
}
#[test]
fn g0_connection_horizons_and_invalid_geometry_are_not_fake_densities() {
    let model = SmithGgx::try_new(0.3, 0.7).unwrap();
    let down = Vec3::new(0.0, 0.0, -1.0);
    let up = down.scale(-1.0);
    assert_eq!(
        model.single_scattering_connection(down, down).unwrap(),
        None
    );
    assert_eq!(
        model
            .connection_at_depth(down, 0.3, Vec3::new(1.0, 0.0, 0.0))
            .unwrap(),
        None
    );
    assert_eq!(model.connection_at_depth(up, 0.3, up).unwrap(), None);
    assert_eq!(
        model.connection_at_depth(down, -1.0, up),
        Err(SmithError::InvalidDepth)
    );
    assert_eq!(
        model.connection_at_depth(Vec3::new(f64::NAN, 0.0, -1.0), 0.0, up),
        Err(SmithError::InvalidDirection)
    );
}
