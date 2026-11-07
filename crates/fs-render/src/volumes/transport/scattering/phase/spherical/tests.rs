use super::*;

#[test]
fn intermediate_underflow_does_not_erase_a_representable_polarized_response() {
    let peak = SphericalPhaseValue::try_new([1e300, 1e-30], C64::ZERO).unwrap();
    let table = SphericalPhaseTable::try_new(&[
        (-1.0, SphericalPhaseValue::ZERO),
        (1.0 - 4.0 * f64::EPSILON, SphericalPhaseValue::ZERO),
        (1.0, peak),
    ])
    .unwrap();
    let normalized = table.scalar_phase().eval_cos(1.0).unwrap();
    let expected = (1e-30 * normalized) / 5e299;
    let weak = table.evaluate(1.0).unwrap().powers_per_sr()[1];
    assert!(weak > 0.0 && expected > 0.0);
    assert!((weak / expected - 1.0).abs() < 1e-7);
    let tiny = f64::from_bits(1);
    let value = SphericalPhaseValue::try_new([2.0, 2.0], C64::ZERO).unwrap();
    assert_eq!(
        value.apply(Stokes::unpolarized(tiny)).unwrap(),
        Stokes::unpolarized(f64::from_bits(2))
    );
    assert_eq!(product_ratio(tiny, 2.0, 1.0).unwrap(), f64::from_bits(2));
    assert_eq!(
        product_ratio(1.0, 2.0, tiny),
        Err(VolumeTransportError::Unrepresentable)
    );
    assert_eq!(
        product_ratio(tiny, 0.5, 1.0),
        Err(VolumeTransportError::Unrepresentable)
    );
    for (a, b, c, expected) in [
        (1e300, 1e-300, 1e-300, 1e300),
        (1e-300, 1e300, 1e300, 1e-300),
        (-1e-300, 1e300, 1e300, -1e-300),
    ] {
        assert!((product_ratio(a, b, c).unwrap() / expected - 1.0).abs() < 5e-15);
    }
}

#[test]
fn weak_spherical_channel_survives_rounded_mueller_difference() {
    let value = SphericalPhaseValue::try_new([2.0, 2e-240], C64::new(0.0, 2e-120)).unwrap();
    let table = SphericalPhaseTable::try_new(&[(-1.0, value), (1.0, value)]).unwrap();
    let response = table.evaluate(0.3).unwrap();
    let p = Stokes {
        i: 1.0,
        q: -1.0,
        u: 0.0,
        v: 0.0,
    };
    assert_eq!(response.mueller().apply(p).i, 0.0);
    let out = response.apply(p).unwrap();
    assert!(out.i > 0.0);
    assert_eq!(out.i, response.powers_per_sr()[1]);
    assert_eq!(out.q, -out.i);
}

#[test]
fn convex_coherency_interpolation_depolarizes_without_amplitude_interference() {
    let a = SphericalPhaseValue::try_new([1.0, 1.0], C64::ONE).unwrap();
    let b = SphericalPhaseValue::try_new([1.0, 1.0], -C64::ONE).unwrap();
    let table = SphericalPhaseTable::try_new(&[(-1.0, a), (1.0, b)]).unwrap();
    let mid = table.evaluate(0.0).unwrap();
    assert_eq!(mid.cross_per_sr(), C64::ZERO);
    let out = mid
        .apply(Stokes {
            i: 1.0,
            q: 0.0,
            u: 1.0,
            v: 0.0,
        })
        .unwrap();
    assert!(out.i > 0.0);
    assert_eq!([out.q, out.u, out.v], [0.0; 3]);
    let area = table
        .scalar_phase()
        .knots_per_sr()
        .windows(2)
        .map(|p| core::f64::consts::PI * (p[1].0 - p[0].0) * (p[0].1 + p[1].1))
        .sum::<f64>();
    assert!((area - 1.0).abs() < 4e-15);
}

#[test]
fn circular_and_linear_channels_follow_the_declared_jones_phase() {
    let value = SphericalPhaseValue::try_new([1.0, 1.0], C64::new(0.0, 1.0)).unwrap();
    let linear = value
        .apply(Stokes {
            i: 1.0,
            q: 0.0,
            u: 1.0,
            v: 0.0,
        })
        .unwrap();
    assert_eq!(
        linear,
        Stokes {
            i: 1.0,
            q: 0.0,
            u: 0.0,
            v: 1.0
        }
    );
    let circular = value
        .apply(Stokes {
            i: 1.0,
            q: 0.0,
            u: 0.0,
            v: 1.0,
        })
        .unwrap();
    assert_eq!(
        circular,
        Stokes {
            i: 1.0,
            q: 0.0,
            u: -1.0,
            v: 0.0
        }
    );
    let dark = SphericalPhaseValue::try_new([1.0, 0.0], C64::ZERO)
        .unwrap()
        .apply(Stokes {
            i: 1.0,
            q: -1.0,
            u: 0.0,
            v: 0.0,
        })
        .unwrap();
    assert_eq!(dark, Stokes::unpolarized(0.0));
}

#[test]
fn coherency_admission_refuses_unphysical_and_unsupported_tables() {
    assert!(SphericalPhaseValue::try_new([1.0, 0.0], C64::ONE).is_err());
    assert!(SphericalPhaseValue::try_new([1.0, 1.0], C64::new(1.01, 0.0)).is_err());
    assert!(SphericalPhaseValue::try_new([-1.0, 1.0], C64::ZERO).is_err());
    assert!(
        SphericalPhaseTable::try_new(&[
            (-1.0, SphericalPhaseValue::ZERO),
            (1.0, SphericalPhaseValue::ZERO)
        ])
        .is_err()
    );
    let value = SphericalPhaseValue::try_new([1.0, 1.0], C64::ONE).unwrap();
    assert!(SphericalPhaseTable::try_new(&[(0.0, value), (1.0, value)]).is_err());
    assert!(
        SphericalPhaseTable::try_new(&[(-1.0, value), (1.0, value)])
            .unwrap()
            .evaluate(f64::NAN)
            .is_err()
    );
}
