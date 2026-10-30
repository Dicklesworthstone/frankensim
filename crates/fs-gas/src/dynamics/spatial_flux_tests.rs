use super::*;

fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 3e-10 * (1.0 + expected.abs()),
        "actual={actual:e}, expected={expected:e}",
    );
}

fn components(f: SpatialFaceFlux) -> [f64; 6] {
    [
        f.mass,
        f.normal_momentum,
        f.transverse_momentum[0],
        f.transverse_momentum[1],
        f.energy,
        f.signal_speed,
    ]
}

fn trace(q: [f64; 5]) -> SpatialPrimitive {
    SpatialPrimitive {
        density: q[0],
        normal_velocity: q[1],
        transverse_velocity: [q[2], q[3]],
        pressure: q[4],
    }
}

fn conserved(g: &FrozenGas, p: SpatialPrimitive) -> SpatialConserved {
    g.spatial_conserved(
        p.density,
        p.normal_velocity,
        p.transverse_velocity,
        p.pressure,
    )
    .unwrap()
}

#[test]
fn g0_spatial_flux_reduces_bitwise_to_the_original_planar_kernel() {
    let g = gas(0.4, 1e18);
    for (left, right) in [
        ([1.0, 0.0, 3.0, 1.0], [4.0, 0.0, -4.0, 1.0]),
        ([1.0, 0.2, 3.0, 1.0], [0.125, 0.2, -4.0, 0.1]),
        ([1.0, -2.0, 0.5, 0.4], [1.0, 2.0, -0.7, 0.4]),
        ([1.2, -8.0, 0.3, 1.7], [0.8, -10.0, -0.1, 1.3]),
    ] {
        let planar = |q: [f64; 4]| g.normal_conserved(q[0], q[1], q[2], q[3]).unwrap();
        let spatial = |q: [f64; 4]| g.spatial_conserved(q[0], q[1], [q[2], 0.0], q[3]).unwrap();
        let p = g
            .hllc_normal_flux(planar(left), planar(right), || false)
            .unwrap();
        let s = g
            .hllc_spatial_flux(spatial(left), spatial(right), || false)
            .unwrap();
        assert_eq!(
            components(s),
            [
                p.mass,
                p.normal_momentum,
                p.transverse_momentum,
                0.0,
                p.energy,
                p.signal_speed
            ],
        );
        assert_eq!(s.used_hll_fallback, p.used_hll_fallback);
    }
}

#[test]
fn g0_spatial_shear_contacts_advect_both_tangents_and_full_energy() {
    let g = gas(0.4, -1e18);
    for u in [-8.0, -0.2, 0.0, 0.7, 10.0] {
        let left = g.spatial_conserved(1.0, u, [3.0, -2.0], 1.0).unwrap();
        let right = g.spatial_conserved(0.3, u, [-4.0, 1.5], 1.0).unwrap();
        let flux = g.hllc_spatial_flux(left, right, || false).unwrap();
        let upwind = if u >= 0.0 { left } else { right };
        close(flux.mass, upwind.normal_momentum);
        close(flux.normal_momentum, upwind.normal_momentum * u + 1.0);
        for i in 0..2 {
            close(
                flux.transverse_momentum[i],
                upwind.transverse_momentum[i] * u,
            );
        }
        close(flux.energy, (upwind.energy + 1.0) * u);
        assert!(!flux.used_hll_fallback);
        close(g.spatial_primitive(upwind, || false).unwrap().pressure, 1.0);
    }
}

#[test]
fn g3_spatial_hllc_and_hll_are_covariant_under_transverse_rotations() {
    let g = gas(0.0, 0.0);
    let rotate = |v: [f64; 2]| [0.6 * v[0] - 0.8 * v[1], 0.8 * v[0] + 0.6 * v[1]];
    for (left, right, fallback) in [
        (
            [0.9, -0.1, 1.5, -0.3, 2.0],
            [0.35, 0.2, -2.0, 0.7, 0.2],
            false,
        ),
        (
            [1.0, -2.0, 0.5, -0.1, 0.4],
            [1.0, 2.0, -0.7, 0.9, 0.4],
            true,
        ),
    ] {
        let (left, right) = (conserved(&g, trace(left)), conserved(&g, trace(right)));
        let original = g.hllc_spatial_flux(left, right, || false).unwrap();
        let rotated = g
            .hllc_spatial_flux(
                SpatialConserved {
                    transverse_momentum: rotate(left.transverse_momentum),
                    ..left
                },
                SpatialConserved {
                    transverse_momentum: rotate(right.transverse_momentum),
                    ..right
                },
                || false,
            )
            .unwrap();
        assert_eq!(original.used_hll_fallback, fallback);
        assert_eq!(rotated.used_hll_fallback, fallback);
        close(rotated.mass, original.mass);
        close(rotated.normal_momentum, original.normal_momentum);
        close(rotated.energy, original.energy);
        close(rotated.signal_speed, original.signal_speed);
        for (actual, expected) in rotated
            .transverse_momentum
            .into_iter()
            .zip(rotate(original.transverse_momentum))
        {
            close(actual, expected);
        }
    }
}

#[test]
fn g3_spatial_hllc_and_hll_are_covariant_under_transverse_boosts() {
    let g = gas(0.0, 0.0);
    for (left, right, fallback) in [
        (
            [0.9, -0.1, 1.5, -0.3, 2.0],
            [0.35, 0.2, -2.0, 0.7, 0.2],
            false,
        ),
        (
            [1.0, -2.0, 0.5, -0.1, 0.4],
            [1.0, 2.0, -0.7, 0.9, 0.4],
            true,
        ),
    ] {
        let (left, right) = (conserved(&g, trace(left)), conserved(&g, trace(right)));
        let original = g.hllc_spatial_flux(left, right, || false).unwrap();
        for boost in [[-3.0, 2.0], [2.0, -1.0]] {
            let half_speed_squared = 0.5 * (boost[0] * boost[0] + boost[1] * boost[1]);
            let shift = |q: SpatialConserved| SpatialConserved {
                transverse_momentum: core::array::from_fn(|i| {
                    q.transverse_momentum[i] + boost[i] * q.density
                }),
                energy: q.energy
                    + boost[0] * q.transverse_momentum[0]
                    + boost[1] * q.transverse_momentum[1]
                    + half_speed_squared * q.density,
                ..q
            };
            let shifted = g
                .hllc_spatial_flux(shift(left), shift(right), || false)
                .unwrap();
            assert_eq!(original.used_hll_fallback, fallback);
            assert_eq!(shifted.used_hll_fallback, fallback);
            close(shifted.mass, original.mass);
            close(shifted.normal_momentum, original.normal_momentum);
            close(shifted.signal_speed, original.signal_speed);
            for i in 0..2 {
                close(
                    shifted.transverse_momentum[i],
                    original.transverse_momentum[i] + boost[i] * original.mass,
                );
            }
            close(
                shifted.energy,
                original.energy
                    + boost[0] * original.transverse_momentum[0]
                    + boost[1] * original.transverse_momentum[1]
                    + half_speed_squared * original.mass,
            );
        }
    }
}

#[test]
fn g0_spatial_hll_source_admission_accounts_for_both_shear_energies() {
    let g = gas(0.0, 0.0);
    let left = g.spatial_conserved(1.0, -2.0, [12.0, 12.0], 0.4).unwrap();
    let right = g.spatial_conserved(1.0, 2.0, [-12.0, -12.0], 0.4).unwrap();
    assert_eq!(
        g.hllc_spatial_flux(left, right, || false),
        Err(Error::Invalid(
            "HLL intermediate state outside admitted source domain",
        ))
    );
    // Either shear component alone lies inside the caloric domain; admitting
    // each independently would miss their combined conversion to heat.
    for tangent in [[12.0, 0.0], [0.0, 12.0]] {
        let left = g.spatial_conserved(1.0, -2.0, tangent, 0.4).unwrap();
        let right = g
            .spatial_conserved(1.0, 2.0, tangent.map(|x| -x), 0.4)
            .unwrap();
        assert!(
            g.hllc_spatial_flux(left, right, || false)
                .unwrap()
                .used_hll_fallback
        );
    }
}

#[test]
fn g0_spatial_stationary_wall_pressure_is_invariant_to_both_shear_components() {
    let g = gas(0.4, 1e18);
    for normal in [-1.0, 1.0] {
        for velocity in [-0.4, 0.0, 2.0] {
            let p = trace([1.2, normal * velocity, 5.0, -3.0, 1.0]);
            let spatial = g
                .stationary_wall_spatial_flux(conserved(&g, p), normal, || false)
                .unwrap();
            let fast = g
                .stationary_wall_spatial_flux_from_primitive(p, normal, || false)
                .unwrap();
            let reference = g
                .stationary_wall_normal_flux(
                    g.normal_conserved(1.2, normal * velocity, 0.0, 1.0)
                        .unwrap(),
                    normal,
                    || false,
                )
                .unwrap();
            close(spatial.normal_momentum, reference.normal_momentum);
            close(spatial.signal_speed, reference.signal_speed);
            assert_eq!(spatial.mass, 0.0);
            assert_eq!(spatial.transverse_momentum, [0.0; 2]);
            assert_eq!(spatial.energy, 0.0);
            for (actual, expected) in components(fast).into_iter().zip(components(spatial)) {
                close(actual, expected);
            }
        }
    }
}

#[test]
fn g3_spatial_primitive_faces_reuse_original_thermodynamics_on_both_flux_branches() {
    for cp in [0.0, 0.4] {
        let g = gas(cp, -1e18);
        for (left, right) in [
            ([1.0, 0.0, 3.0, -2.0, 1.0], [4.0, 0.0, -4.0, 3.0, 1.0]),
            ([1.0, 0.2, 3.0, -2.0, 1.0], [0.125, 0.2, -4.0, 3.0, 0.1]),
            ([1.0, -2.0, 0.5, -0.1, 0.4], [1.0, 2.0, -0.7, 0.9, 0.4]),
            ([1.2, -8.0, 0.3, -0.1, 1.7], [0.8, -10.0, -0.1, 0.9, 1.3]),
        ] {
            let (left, right) = (trace(left), trace(right));
            let original = g
                .hllc_spatial_flux(conserved(&g, left), conserved(&g, right), || false)
                .unwrap();
            let fast = g
                .hllc_spatial_flux_from_primitive(left, right, || false)
                .unwrap();
            assert_eq!(fast.used_hll_fallback, original.used_hll_fallback);
            for (actual, expected) in components(fast).into_iter().zip(components(original)) {
                close(actual, expected);
            }
        }
    }
}

#[test]
fn g4_spatial_traces_refuse_lost_heat_invalid_sources_and_cancellation() {
    let g = gas(0.4, 0.0);
    let p = trace([1.0, 0.2, 0.3, -0.1, 1.0]);
    let q = conserved(&g, p);
    for invalid in [
        SpatialPrimitive {
            pressure: 100.0,
            ..p
        },
        SpatialPrimitive { density: 0.0, ..p },
        SpatialPrimitive {
            transverse_velocity: [0.0, f64::NAN],
            ..p
        },
        SpatialPrimitive {
            transverse_velocity: [0.0, 1e12],
            ..p
        },
    ] {
        assert!(
            g.spatial_conserved(
                invalid.density,
                invalid.normal_velocity,
                invalid.transverse_velocity,
                invalid.pressure
            )
            .is_err()
        );
        assert!(
            g.hllc_spatial_flux_from_primitive(invalid, p, || false)
                .is_err()
        );
        assert!(
            g.stationary_wall_spatial_flux_from_primitive(invalid, 1.0, || false)
                .is_err()
        );
    }
    for invalid in [
        SpatialConserved { density: 0.0, ..q },
        SpatialConserved { energy: -1.0, ..q },
        SpatialConserved {
            transverse_momentum: [0.0, f64::INFINITY],
            ..q
        },
    ] {
        assert!(g.spatial_primitive(invalid, || false).is_err());
    }
    assert!(g.stationary_wall_spatial_flux(q, 0.0, || false).is_err());
    assert_eq!(g.spatial_primitive(q, || true), Err(Error::Cancelled));
    assert_eq!(g.hllc_spatial_flux(q, q, || true), Err(Error::Cancelled));
    assert_eq!(
        g.hllc_spatial_flux_from_primitive(p, p, || true),
        Err(Error::Cancelled)
    );
    assert_eq!(
        g.stationary_wall_spatial_flux(q, 1.0, || true),
        Err(Error::Cancelled)
    );
    assert_eq!(
        g.stationary_wall_spatial_flux_from_primitive(p, 1.0, || true),
        Err(Error::Cancelled)
    );
    let mut polls = 0;
    assert_eq!(
        g.hllc_spatial_flux(q, q, || {
            polls += 1;
            polls == 4
        }),
        Err(Error::Cancelled)
    );
}
