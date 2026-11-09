use super::*;

fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 3e-10 * (1.0 + expected.abs()),
        "actual={actual:e}, expected={expected:e}",
    );
}

#[test]
fn g1_compact_spatial_reservoir_matches_stagnation_choking_and_transports_full_shear() {
    let gas = gas(0.0, -1e18);
    let throat_temperature = 2.0_f64 / 2.4;
    let throat_density = throat_temperature.powf(2.5);
    let throat_pressure = throat_temperature.powf(3.5);
    let throat_speed = (1.4 * throat_temperature).sqrt();
    let reservoir = gas.conserved(1.0, 0.0, 1.0).unwrap();
    let enthalpy = reservoir.energy + 1.0;
    for normal in [-1.0, 1.0] {
        for shear in [[2.0, -3.0], [-3.0, 2.0]] {
            let inlet = SpatialPrimitive {
                density: 0.01,
                normal_velocity: 0.0,
                transverse_velocity: shear,
                pressure: 0.01,
            };
            let mass = -normal * throat_density * throat_speed;
            for flux in [
                gas.reservoir_spatial_flux(conserved(&gas, inlet), reservoir, normal, || false)
                    .unwrap(),
                gas.reservoir_spatial_flux_from_primitive(inlet, reservoir, normal, || false)
                    .unwrap(),
            ] {
                close(flux.mass, mass);
                close(
                    flux.normal_momentum,
                    throat_density * throat_speed.powi(2) + throat_pressure,
                );
                assert_eq!(flux.transverse_momentum, [0.0; 2]);
                close(flux.energy, mass * enthalpy);
            }
            let outlet = SpatialPrimitive {
                density: 1.0,
                normal_velocity: normal * 0.4,
                transverse_velocity: shear,
                pressure: 1.0,
            };
            let bath = gas.conserved(0.8, 0.0, 0.8).unwrap();
            let axial = gas
                .reservoir_flux(
                    gas.conserved(1.0, normal * 0.4, 1.0).unwrap(),
                    bath,
                    normal,
                    || false,
                )
                .unwrap();
            for flux in [
                gas.reservoir_spatial_flux(conserved(&gas, outlet), bath, normal, || false)
                    .unwrap(),
                gas.reservoir_spatial_flux_from_primitive(outlet, bath, normal, || false)
                    .unwrap(),
            ] {
                close(flux.mass, axial.mass);
                close(flux.normal_momentum, axial.momentum);
                for (component, velocity) in flux.transverse_momentum.into_iter().zip(shear) {
                    close(component, axial.mass * velocity);
                }
                close(
                    flux.energy,
                    axial.energy + 0.5 * axial.mass * (shear[0].powi(2) + shear[1].powi(2)),
                );
            }
        }
    }
}

#[test]
fn g4_compact_spatial_reservoir_refuses_nonstationary_and_unresolved_source_states() {
    let gas = gas(0.4, 0.0);
    let interior = trace([1.0, -0.4, 2.0, -3.0, 1.0]);
    let reservoir = gas.conserved(2.0, 0.0, 2.0).unwrap();
    // Inward interior motion is now physically resolved by a complete reservoir.
    assert!(
        gas.reservoir_spatial_flux_from_primitive(interior, reservoir, 1.0, || false)
            .is_ok()
    );
    for invalid in [
        Conserved {
            momentum: 0.1,
            ..reservoir
        },
        Conserved {
            density: -1.0,
            ..reservoir
        },
        Conserved {
            energy: f64::NAN,
            ..reservoir
        },
    ] {
        assert!(
            gas.reservoir_spatial_flux_from_primitive(interior, invalid, 1.0, || false)
                .is_err()
        );
    }
    let unresolved = SpatialPrimitive {
        transverse_velocity: [1e20, 0.0],
        ..interior
    };
    assert!(
        gas.reservoir_spatial_flux_from_primitive(unresolved, reservoir, 1.0, || false)
            .is_err()
    );
    assert!(
        gas.reservoir_spatial_flux_from_primitive(interior, reservoir, 0.0, || false)
            .is_err()
    );
    assert_eq!(
        gas.reservoir_spatial_flux_from_primitive(interior, reservoir, 1.0, || true),
        Err(Error::Cancelled),
    );
}

#[test]
fn g1_spatial_characteristic_outlet_passes_simple_acoustic_waves_and_shear_contacts() {
    let gas = gas(0.0, -1e18);
    let c0 = 1.4_f64.sqrt();
    for normal in [-1.0, 1.0] {
        for temperature in [0.95_f64, 1.0, 1.05] {
            // Exact gamma=1.4 outgoing simple wave: p=T^3.5, rho=T^2.5,
            // u-u_ref=5(c-c_ref). Shear velocities are arbitrary outgoing
            // contact data and must not alter the pressure characteristic.
            let density = temperature.powf(2.5);
            let pressure = temperature.powf(3.5);
            let velocity = normal * (0.4 + 5.0 * c0 * (temperature.sqrt() - 1.0));
            for tangents in [[0.0, 0.0], [2.0, -3.0], [-3.0, 2.0]] {
                let primitive = SpatialPrimitive {
                    density,
                    normal_velocity: velocity,
                    transverse_velocity: tangents,
                    pressure,
                };
                let state = conserved(&gas, primitive);
                let expected = [
                    state.normal_momentum,
                    state.normal_momentum * velocity + pressure,
                    state.transverse_momentum[0] * velocity,
                    state.transverse_momentum[1] * velocity,
                    (state.energy + pressure) * velocity,
                ];
                for flux in [
                    gas.characteristic_outlet_spatial_flux(
                        state,
                        normal,
                        1.0,
                        normal * 0.4,
                        || false,
                    )
                    .unwrap(),
                    gas.characteristic_outlet_spatial_flux_from_primitive(
                        primitive,
                        normal,
                        1.0,
                        normal * 0.4,
                        || false,
                    )
                    .unwrap(),
                ] {
                    for (actual, reference) in components(flux).into_iter().zip(expected) {
                        close(actual, reference);
                    }
                }
            }
        }
    }
}

#[test]
fn g4_spatial_characteristic_outlet_refuses_backflow_and_lost_caloric_precision() {
    let gas = gas(0.0, 0.0);
    let primitive = trace([1.0, 0.4, 2.0, -3.0, 1.0]);
    assert_eq!(
        gas.characteristic_outlet_spatial_flux_from_primitive(primitive, 1.0, 1.0, 0.4, || true),
        Err(Error::Cancelled)
    );
    for invalid in [
        SpatialPrimitive {
            normal_velocity: -0.4,
            ..primitive
        },
        SpatialPrimitive {
            transverse_velocity: [f64::NAN, 0.0],
            ..primitive
        },
        SpatialPrimitive {
            transverse_velocity: [1e20, 0.0],
            ..primitive
        },
    ] {
        assert!(
            gas.characteristic_outlet_spatial_flux_from_primitive(
                invalid, 1.0, 1.0, 0.4, || false,
            )
            .is_err()
        );
    }
    for (normal, pressure, velocity) in [(0.0, 1.0, 0.4), (1.0, 0.0, 0.4), (1.0, 1.0, -0.4)] {
        assert!(
            gas.characteristic_outlet_spatial_flux_from_primitive(
                primitive,
                normal,
                pressure,
                velocity,
                || false,
            )
            .is_err()
        );
    }
}

#[test]
fn g3_pressure_outlet_retains_both_shear_components_and_full_spatial_work() {
    let gas = gas(0.4, -1e18);
    for normal in [-1.0, 1.0] {
        let axial = gas.conserved(1.0, normal * 0.4, 1.0).unwrap();
        let reference = gas
            .pressure_outlet_flux(axial, normal, 0.8, || false)
            .unwrap();
        for tangent in [[0.0, 0.0], [2.0, -3.0], [-3.0, 2.0]] {
            let primitive = SpatialPrimitive {
                density: 1.0,
                normal_velocity: normal * 0.4,
                transverse_velocity: tangent,
                pressure: 1.0,
            };
            let state = conserved(&gas, primitive);
            let flux = gas
                .pressure_outlet_spatial_flux(state, normal, 0.8, || false)
                .unwrap();
            close(flux.mass, reference.mass);
            close(flux.normal_momentum, reference.momentum);
            for (index, velocity) in tangent.into_iter().enumerate() {
                close(flux.transverse_momentum[index], reference.mass * velocity);
            }
            close(
                flux.energy,
                reference.energy + 0.5 * reference.mass * (tangent[0].powi(2) + tangent[1].powi(2)),
            );
            let reconstructed = gas
                .pressure_outlet_spatial_flux_from_primitive(primitive, normal, 0.8, || false)
                .unwrap();
            for (actual, expected) in components(reconstructed).into_iter().zip(components(flux)) {
                close(actual, expected);
            }
            assert!(!flux.used_hll_fallback);
        }
    }
}

#[test]
fn g4_spatial_pressure_outlet_preserves_original_trace_admission() {
    let gas = gas(0.0, 0.0);
    let primitive = trace([1.0, 0.4, 2.0, -3.0, 1.0]);
    assert_eq!(
        gas.pressure_outlet_spatial_flux_from_primitive(primitive, 1.0, 0.8, || true),
        Err(Error::Cancelled)
    );
    for invalid in [
        SpatialPrimitive {
            transverse_velocity: [f64::NAN, 0.0],
            ..primitive
        },
        SpatialPrimitive {
            transverse_velocity: [1e20, 0.0],
            ..primitive
        },
        SpatialPrimitive {
            normal_velocity: -0.4,
            ..primitive
        },
    ] {
        assert!(
            gas.pressure_outlet_spatial_flux_from_primitive(invalid, 1.0, 1.0, || false)
                .is_err()
        );
    }
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
fn g3_ale_spatial_momentum_and_work_are_covariant_on_contact_shock_and_hll() {
    for cp in [0.0, 0.4] {
        let g = gas(cp, 1e18);
        for (l, r) in [
            ([1.0, 0.2, 0.3, -0.4, 1.0], [0.3, 0.2, -0.7, 0.8, 1.0]),
            ([1.0, 0.0, 0.3, -0.4, 1.0], [0.125, 0.0, -0.7, 0.8, 0.1]),
            ([1.0, -2.0, 0.3, -0.4, 0.4], [1.0, 2.0, -0.7, 0.8, 0.4]),
        ] {
            let (left, right) = (trace(l), trace(r));
            let speed = 0.6;
            let f = g
                .moving_hllc_spatial_flux_from_primitive(left, right, speed, || false)
                .unwrap();
            // A uniform boost of both fluids AND the face preserves the
            // relative Riemann problem, including the chosen fallback.
            let boost = [1.25, -0.5, 0.75];
            let shift = |p: SpatialPrimitive| SpatialPrimitive {
                normal_velocity: p.normal_velocity + boost[0],
                transverse_velocity: [
                    p.transverse_velocity[0] + boost[1],
                    p.transverse_velocity[1] + boost[2],
                ],
                ..p
            };
            let transformed = g
                .moving_hllc_spatial_flux_from_primitive(
                    shift(left),
                    shift(right),
                    speed + boost[0],
                    || false,
                )
                .unwrap();
            close(transformed.mass, f.mass);
            close(
                transformed.normal_momentum,
                f.normal_momentum + boost[0] * f.mass,
            );
            for axis in 0..2 {
                close(
                    transformed.transverse_momentum[axis],
                    f.transverse_momentum[axis] + boost[axis + 1] * f.mass,
                );
            }
            close(
                transformed.energy,
                f.energy
                    + boost[0] * f.normal_momentum
                    + boost[1] * f.transverse_momentum[0]
                    + boost[2] * f.transverse_momentum[1]
                    + 0.5 * boost.iter().map(|b| b * b).sum::<f64>() * f.mass,
            );
            close(transformed.signal_speed, f.signal_speed);
            assert_eq!(transformed.used_hll_fallback, f.used_hll_fallback);
            if l[1] == -2.0 {
                assert!(f.used_hll_fallback);
            }
            let ql = conserved(&g, left);
            let qr = conserved(&g, right);
            let decoded = g.moving_hllc_spatial_flux(ql, qr, speed, || false).unwrap();
            for (actual, expected) in components(decoded).into_iter().zip(components(f)) {
                close(actual, expected);
            }
            assert_eq!(
                g.moving_hllc_spatial_flux(ql, qr, 0.0, || false).unwrap(),
                g.hllc_spatial_flux(ql, qr, || false).unwrap()
            );
        }
    }
}

#[test]
fn g0_ale_spatial_reduces_to_planar_and_has_the_comoving_wall_limit() {
    let g = gas(0.4, -1e18);
    let left = trace([1.1, 0.4, -0.3, 0.0, 1.0]);
    let right = trace([0.7, -0.2, 0.8, 0.0, 0.5]);
    let planar = |p: SpatialPrimitive| NormalPrimitive {
        density: p.density,
        normal_velocity: p.normal_velocity,
        transverse_velocity: p.transverse_velocity[0],
        pressure: p.pressure,
    };
    for speed in [-1.0, 0.0, 0.7] {
        let p = g
            .moving_hllc_normal_flux_from_primitive(planar(left), planar(right), speed, || false)
            .unwrap();
        let s = g
            .moving_hllc_spatial_flux_from_primitive(left, right, speed, || false)
            .unwrap();
        assert_eq!(
            components(s),
            [
                p.mass,
                p.normal_momentum,
                p.transverse_momentum,
                0.0,
                p.energy,
                p.signal_speed,
            ]
        );
        assert_eq!(s.used_hll_fallback, p.used_hll_fallback);
    }
    let left = SpatialPrimitive {
        transverse_velocity: [-0.3, 0.9],
        ..left
    };
    let right = SpatialPrimitive {
        density: 0.3,
        transverse_velocity: [0.8, -0.5],
        ..left
    };
    let flux = g
        .moving_hllc_spatial_flux_from_primitive(left, right, left.normal_velocity, || false)
        .unwrap();
    assert_eq!(flux.mass, 0.0);
    assert_eq!(flux.transverse_momentum, [0.0; 2]);
    close(flux.normal_momentum, left.pressure);
    close(flux.energy, left.pressure * left.normal_velocity);
    let wall = g
        .moving_wall_spatial_flux_from_primitive(left, 1.0, left.normal_velocity, || false)
        .unwrap();
    close(wall.normal_momentum, flux.normal_momentum);
    close(wall.energy, flux.energy);
}

#[test]
fn g4_ale_faces_refuse_nonfinite_speeds_caloric_loss_and_cancelled_work() {
    let g = gas(0.4, 0.0);
    let p = trace([1.0, 0.3, -0.2, 0.4, 1.0]);
    let q = conserved(&g, p);
    for velocity in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY, 1e100] {
        assert!(
            g.moving_hllc_spatial_flux(q, q, velocity, || false)
                .is_err()
        );
        assert!(
            g.moving_hllc_spatial_flux_from_primitive(p, p, velocity, || false)
                .is_err()
        );
    }
    assert_eq!(
        g.moving_hllc_spatial_flux(q, q, 0.1, || true),
        Err(Error::Cancelled)
    );
    let mut polls = 0;
    assert_eq!(
        g.moving_hllc_spatial_flux_from_primitive(p, p, 0.1, || {
            polls += 1;
            polls == 5
        }),
        Err(Error::Cancelled)
    );
    // Failed and cancelled evaluations have no retained state or scratch;
    // the original traces still supply the identical next accepted flux.
    assert_eq!(
        g.moving_hllc_spatial_flux(q, q, 0.1, || false).unwrap(),
        g.moving_hllc_spatial_flux(q, q, 0.1, || false).unwrap()
    );
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
