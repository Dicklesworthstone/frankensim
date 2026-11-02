//! Independent centre-regular Fourier and volume identities.
use super::{RadialKind, SolidCore};

fn close(actual: f64, expected: f64) {
    if expected == 0.0 {
        assert_eq!(actual, 0.0);
        return;
    }
    assert!(
        (actual - expected).abs() <= 64.0 * f64::EPSILON * expected.abs(),
        "{actual} != {expected}"
    );
}

#[test]
fn solid_core_matches_regular_flux_volume_and_kirchhoff_integrals() {
    for (kind, dimension, area_factor) in [
        (
            RadialKind::Cylindrical,
            2.0,
            2.0 * core::f64::consts::PI * 0.7,
        ),
        (RadialKind::Spherical, 3.0, 4.0 * core::f64::consts::PI),
    ] {
        for radius in [0.001_f64, 0.2, 2.0, 1000.0] {
            let core = SolidCore::new(kind, radius).unwrap();
            let outer_area = area_factor * radius.powi(dimension as i32 - 1);
            for source in [-17.0, 0.0, 21.0] {
                let generated = source * core.volume();
                close(
                    generated * outer_area,
                    source * area_factor * radius.powi(dimension as i32) / dimension,
                );
                for fraction in [0.0, 1e-12, 0.125, 0.5, 1.0] {
                    let r = radius * fraction;
                    close(core.local_flux(generated, r), source * r / dimension);
                    close(
                        core.potential_drop(generated, r),
                        source * r * r / (2.0 * dimension),
                    );
                    close(
                        core.local_flux(generated, r) * core.area_ratio(r) * outer_area,
                        source * area_factor * r.powi(dimension as i32) / dimension,
                    );
                }
            }
        }
    }
}

#[test]
fn solid_core_rejects_invalid_geometry_and_preserves_physical_scaling() {
    for kind in [RadialKind::Cylindrical, RadialKind::Spherical] {
        for radius in [
            0.0,
            -1.0,
            f64::NAN,
            f64::INFINITY,
            f64::NEG_INFINITY,
            f64::from_bits(1),
        ] {
            assert!(SolidCore::new(kind, radius).is_none());
        }
        let a = SolidCore::new(kind, 0.3).unwrap();
        let b = SolidCore::new(kind, 0.6).unwrap();
        for fraction in [0.0, 0.1, 0.4, 1.0] {
            let ga = 12.0 * a.volume();
            let gb = 12.0 * b.volume();
            close(
                b.local_flux(gb, fraction * b.radius()),
                2.0 * a.local_flux(ga, fraction * a.radius()),
            );
            close(
                b.potential_drop(gb, fraction * b.radius()),
                4.0 * a.potential_drop(ga, fraction * a.radius()),
            );
            close(
                b.area_ratio(fraction * b.radius()),
                a.area_ratio(fraction * a.radius()),
            );
        }
    }
}
