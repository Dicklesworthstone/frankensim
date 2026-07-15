//! Gauntlet G3 relations for production rendering kernels.
//!
//! The declared relation here supplements the chart backend's existing
//! translated-scene frame-invariance pin rather than replacing it.

use fs_propcheck::metamorphic::{RelationCase, Tolerance, check_relation, unit_rescaling};
use fs_render::Lambertian;

#[test]
fn g3_lambertian_furnace_tracks_radiance_unit_rescaling() {
    let operator =
        |&(albedo, incident): &(f64, f64)| Lambertian { albedo }.furnace_radiance(incident, 16);
    let relation = unit_rescaling(
        "lambertian-radiance-scale-equivariance",
        Tolerance::AbsoluteRelative {
            max_abs: 2.0e-12,
            max_relative: 2.0e-12,
        },
        |&(albedo, incident): &(f64, f64), &exponent: &i64| {
            let scale = 2.0f64.powi(exponent as i32);
            (albedo, incident * scale)
        },
        |&base: &f64, &transformed: &f64, &exponent: &i64, tolerance| {
            let scale = 2.0f64.powi(exponent as i32);
            tolerance.evaluate_scalar(base * scale, transformed)
        },
    );

    check_relation(
        "fs-render::Lambertian::furnace_radiance",
        0x2ACE_0006,
        384,
        |stream| {
            RelationCase::new(
                (stream.f64_in(0.05, 0.95), stream.f64_in(0.125, 32.0)),
                stream.int_in(-3, 3),
            )
        },
        &operator,
        &relation,
    );
}
