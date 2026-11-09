use super::*;

fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 2e-10 * (1.0 + expected.abs()),
        "actual={actual:e}, expected={expected:e}",
    );
}

#[test]
fn g3_compact_planar_reservoir_selects_the_physical_shear_donor_across_reversal() {
    for cp in [0.0, 0.4] {
        let gas = gas(cp, 1e18);
        for normal in [-1.0, 1.0] {
            for (density, pressure, speed, reservoir_pressure) in [
                (1.0, 1.0, 0.4, 0.8),
                (1.0, 1.0, 0.4, 2.0),
                (0.01, 0.01, 0.0, 1.0),
                (1.0, 1.0, 0.4, 0.001),
                (1.0, 1.0, 3.0, 1.0),
            ] {
                let velocity = normal * speed;
                let reservoir = gas
                    .conserved(reservoir_pressure, 0.0, reservoir_pressure)
                    .unwrap();
                let axial = gas.conserved(density, velocity, pressure).unwrap();
                let reference = gas
                    .reservoir_flux(axial, reservoir, normal, || false)
                    .unwrap();
                for tangent in [-3.0, 0.0, 2.0] {
                    let primitive = NormalPrimitive {
                        density,
                        normal_velocity: velocity,
                        transverse_velocity: tangent,
                        pressure,
                    };
                    let state = gas
                        .normal_conserved(density, velocity, tangent, pressure)
                        .unwrap();
                    let donor_tangent = if normal * reference.mass < 0.0 {
                        0.0
                    } else {
                        tangent
                    };
                    for flux in [
                        gas.reservoir_normal_flux(state, reservoir, normal, || false)
                            .unwrap(),
                        gas.reservoir_normal_flux_from_primitive(
                            primitive,
                            reservoir,
                            normal,
                            || false,
                        )
                        .unwrap(),
                    ] {
                        close(flux.mass, reference.mass);
                        close(flux.normal_momentum, reference.momentum);
                        close(flux.transverse_momentum, reference.mass * donor_tangent);
                        close(
                            flux.energy,
                            reference.energy + 0.5 * reference.mass * donor_tangent.powi(2),
                        );
                        close(flux.signal_speed, reference.signal_speed);
                        assert!(!flux.used_hll_fallback);
                    }
                }
            }
        }
    }
}

#[test]
fn g0_ale_planar_flux_retains_stationary_and_original_axial_transport() {
    let g = gas(0.4, 1e18);
    for (left, right) in [
        ([1.0, 0.2, 1.0], [0.3, 0.2, 1.0]),
        ([1.0, 0.0, 1.0], [0.125, 0.0, 0.1]),
        ([1.0, -2.0, 0.4], [1.0, 2.0, 0.4]),
    ] {
        let l = g.conserved(left[0], left[1], left[2]).unwrap();
        let r = g.conserved(right[0], right[1], right[2]).unwrap();
        let pl = g.primitive(l, || false).unwrap();
        let pr = g.primitive(r, || false).unwrap();
        for speed in [-0.7, 0.0, 0.9] {
            let planar = g
                .moving_hllc_normal_flux(l.into(), r.into(), speed, || false)
                .unwrap();
            let axial = g
                .moving_interface(l, pl, r, pr, speed, &mut || false)
                .unwrap();
            assert_eq!(FaceFlux::from(planar), axial);
            assert_eq!(planar.transverse_momentum, 0.0);
            if speed == 0.0 {
                assert_eq!(
                    planar,
                    g.hllc_normal_flux(l.into(), r.into(), || false).unwrap()
                );
            }
        }
    }
}

#[test]
fn g1_ale_planar_transport_uses_relative_advection_and_laboratory_pressure_work() {
    let g = gas(0.4, -1e18);
    let p = NormalPrimitive {
        density: 1.2,
        normal_velocity: 0.3,
        transverse_velocity: -0.7,
        pressure: 1.7,
    };
    let q = g
        .normal_conserved(
            p.density,
            p.normal_velocity,
            p.transverse_velocity,
            p.pressure,
        )
        .unwrap();
    let fixed = g.hllc_normal_flux_from_primitive(p, p, || false).unwrap();
    for speed in [-0.7, 0.0, p.normal_velocity, 1.3] {
        let f = g
            .moving_hllc_normal_flux_from_primitive(p, p, speed, || false)
            .unwrap();
        close(f.mass, fixed.mass - speed * q.density);
        close(
            f.normal_momentum,
            fixed.normal_momentum - speed * q.normal_momentum,
        );
        close(
            f.transverse_momentum,
            fixed.transverse_momentum - speed * q.transverse_momentum,
        );
        close(f.energy, fixed.energy - speed * q.energy);
        let from_conserved = g.moving_hllc_normal_flux(q, q, speed, || false).unwrap();
        for (actual, expected) in [
            (from_conserved.mass, f.mass),
            (from_conserved.normal_momentum, f.normal_momentum),
            (from_conserved.transverse_momentum, f.transverse_momentum),
            (from_conserved.energy, f.energy),
        ] {
            close(actual, expected);
        }
    }
    // A face moving with a pressure/shear contact carries no material mass
    // or shear momentum, but it still transfers mechanical pressure work.
    let other = NormalPrimitive {
        density: 0.4,
        transverse_velocity: 1.1,
        ..p
    };
    let f = g
        .moving_hllc_normal_flux_from_primitive(p, other, p.normal_velocity, || false)
        .unwrap();
    assert_eq!(f.mass, 0.0);
    assert_eq!(f.transverse_momentum, 0.0);
    close(f.normal_momentum, p.pressure);
    close(f.energy, p.pressure * p.normal_velocity);
    assert!(!f.used_hll_fallback);
    for normal in [-1.0, 1.0] {
        let wall = g
            .moving_wall_normal_flux_from_primitive(p, normal, p.normal_velocity, || false)
            .unwrap();
        close(wall.normal_momentum, f.normal_momentum);
        close(wall.energy, f.energy);
    }
}

#[test]
fn g3_planar_characteristic_outlet_transports_shear_through_each_acoustic_regime() {
    for cp in [0.0, 0.4] {
        let gas = gas(cp, 1e18);
        for normal in [-1.0, 1.0] {
            // Smooth subsonic discharge, sonic expansion and supersonic
            // causality share the existing independently verified axial law.
            for (speed, reference_pressure) in [(0.4, 0.9), (0.4, 0.01), (3.0, 1e40)] {
                let velocity = normal * speed;
                let axial = gas.conserved(1.0, velocity, 1.0).unwrap();
                let reference = gas
                    .characteristic_outlet_flux(
                        axial,
                        normal,
                        reference_pressure,
                        normal * 0.2,
                        || false,
                    )
                    .unwrap();
                for tangent in [-3.0, 0.0, 2.0] {
                    let primitive = NormalPrimitive {
                        density: 1.0,
                        normal_velocity: velocity,
                        transverse_velocity: tangent,
                        pressure: 1.0,
                    };
                    let state = gas.normal_conserved(1.0, velocity, tangent, 1.0).unwrap();
                    for flux in [
                        gas.characteristic_outlet_normal_flux(
                            state,
                            normal,
                            reference_pressure,
                            normal * 0.2,
                            || false,
                        )
                        .unwrap(),
                        gas.characteristic_outlet_normal_flux_from_primitive(
                            primitive,
                            normal,
                            reference_pressure,
                            normal * 0.2,
                            || false,
                        )
                        .unwrap(),
                    ] {
                        close(flux.mass, reference.mass);
                        close(flux.normal_momentum, reference.momentum);
                        close(flux.transverse_momentum, reference.mass * tangent);
                        close(
                            flux.energy,
                            reference.energy + 0.5 * reference.mass * tangent * tangent,
                        );
                        close(flux.signal_speed, reference.signal_speed);
                        assert!(!flux.used_hll_fallback);
                    }
                }
            }
        }
    }
}

#[test]
fn g3_pressure_outlet_transports_planar_shear_energy_and_reconstruction_consistently() {
    let gas = gas(0.4, 1e18);
    for normal in [-1.0, 1.0] {
        let axial = gas.conserved(1.0, normal * 0.4, 1.0).unwrap();
        let reference = gas
            .pressure_outlet_flux(axial, normal, 0.8, || false)
            .unwrap();
        for tangent in [-3.0, 0.0, 2.0] {
            let state = gas
                .normal_conserved(1.0, normal * 0.4, tangent, 1.0)
                .unwrap();
            let flux = gas
                .pressure_outlet_normal_flux(state, normal, 0.8, || false)
                .unwrap();
            close(flux.mass, reference.mass);
            close(flux.normal_momentum, reference.momentum);
            close(flux.transverse_momentum, reference.mass * tangent);
            close(
                flux.energy,
                reference.energy + 0.5 * reference.mass * tangent * tangent,
            );
            let reconstructed = gas
                .pressure_outlet_normal_flux_from_primitive(
                    NormalPrimitive {
                        density: 1.0,
                        normal_velocity: normal * 0.4,
                        transverse_velocity: tangent,
                        pressure: 1.0,
                    },
                    normal,
                    0.8,
                    || false,
                )
                .unwrap();
            close(reconstructed.mass, flux.mass);
            close(reconstructed.normal_momentum, flux.normal_momentum);
            close(reconstructed.transverse_momentum, flux.transverse_momentum);
            close(reconstructed.energy, flux.energy);
            close(reconstructed.signal_speed, flux.signal_speed);
        }
    }
}

#[test]
fn g0_planar_flux_retains_original_one_dimensional_results() {
    let g = gas(0.0, 1e16);
    for (left, right) in [
        ([1.0, 0.0, 1.0], [4.0, 0.0, 1.0]),
        ([1.0, 0.0, 1.0], [0.125, 0.0, 0.1]),
        ([1.0, -2.0, 0.4], [1.0, 2.0, 0.4]),
        ([1.2, 10.0, 1.7], [0.8, 8.0, 1.3]),
    ] {
        let l = g.conserved(left[0], left[1], left[2]).unwrap();
        let r = g.conserved(right[0], right[1], right[2]).unwrap();
        let axial = g.hllc_flux(l, r, || false).unwrap();
        let planar = g.hllc_normal_flux(l.into(), r.into(), || false).unwrap();
        assert_eq!(FaceFlux::from(planar), axial);
        assert_eq!(planar.transverse_momentum, 0.0);
    }
}

#[test]
fn g0_planar_shear_contact_advects_full_energy_without_pressure_work_duplication() {
    let g = gas(0.4, -1e16);
    for u in [-8.0, -0.2, 0.0, 0.7, 10.0] {
        let l = g.normal_conserved(1.0, u, 3.0, 1.0).unwrap();
        let r = g.normal_conserved(0.3, u, -4.0, 1.0).unwrap();
        let f = g.hllc_normal_flux(l, r, || false).unwrap();
        let upwind = if u >= 0.0 { l } else { r };
        close(f.mass, upwind.normal_momentum);
        close(f.normal_momentum, upwind.normal_momentum * u + 1.0);
        close(f.transverse_momentum, upwind.transverse_momentum * u);
        close(f.energy, (upwind.energy + 1.0) * u);
        assert!(!f.used_hll_fallback);
        let p = g.normal_primitive(upwind, || false).unwrap();
        close(p.velocity, u);
        close(p.pressure, 1.0);
    }
}

#[test]
fn g0_exact_moving_pressure_and_shear_contacts_use_the_physical_upwind_flux() {
    for cp in [0.0, 0.4] {
        let g = gas(cp, -1e16);
        for pressure in [0.7, 1.3, 2.9, 7.1] {
            for velocity in [-0.31, 0.17] {
                let left = NormalPrimitive {
                    density: 1.0,
                    normal_velocity: velocity,
                    transverse_velocity: 0.43,
                    pressure,
                };
                let right = NormalPrimitive {
                    density: 2.0,
                    transverse_velocity: -0.91,
                    ..left
                };
                let upwind = if velocity > 0.0 { left } else { right };
                let q = g
                    .normal_conserved(
                        upwind.density,
                        upwind.normal_velocity,
                        upwind.transverse_velocity,
                        upwind.pressure,
                    )
                    .unwrap();
                let pl = g.from_rho_u_p(left.density, velocity, pressure).unwrap();
                let pr = g.from_rho_u_p(right.density, velocity, pressure).unwrap();
                assert_eq!(pl.pressure, pr.pressure);
                let flux = g
                    .hllc_normal_flux_from_primitive(left, right, || false)
                    .unwrap();
                assert_eq!(flux.mass, q.normal_momentum);
                assert_eq!(
                    flux.normal_momentum,
                    q.normal_momentum * velocity + pl.pressure
                );
                assert_eq!(flux.transverse_momentum, q.transverse_momentum * velocity);
                assert_eq!(flux.energy, (q.energy + pl.pressure) * velocity);
                assert_eq!(
                    flux.signal_speed,
                    velocity.abs() + pl.sound_speed.max(pr.sound_speed)
                );
                assert!(!flux.used_hll_fallback);
            }
        }
    }
}

#[test]
fn g0_exact_moving_contact_at_the_cold_source_endpoint_needs_no_hll_fallback() {
    let g = gas(0.0, 0.0);
    // Both traces are admitted and the denser trace is exactly at Tmin.
    // Reconstructing its unchanged HLLC star rounds zero sensible heat
    // negative; the spurious HLL fallback reverses this positive mass flux.
    let left = NormalPrimitive {
        density: 500.0,
        normal_velocity: 0.001,
        transverse_velocity: 0.0,
        pressure: 1.0,
    };
    let right = NormalPrimitive {
        density: 1000.0,
        ..left
    };
    let pl = g
        .from_rho_u_p(left.density, left.normal_velocity, left.pressure)
        .unwrap();
    let pr = g
        .from_rho_u_p(right.density, right.normal_velocity, right.pressure)
        .unwrap();
    assert_eq!(pl.pressure, pr.pressure);
    assert_eq!(pr.temperature, g.domain().minimum_temperature);
    let q = g
        .normal_conserved(left.density, left.normal_velocity, 0.0, left.pressure)
        .unwrap();
    let flux = g
        .hllc_normal_flux_from_primitive(left, right, || false)
        .unwrap();
    assert!(!flux.used_hll_fallback, "{flux:?}");
    assert_eq!(flux.mass, q.normal_momentum);
    assert_eq!(flux.energy, (q.energy + pl.pressure) * pl.velocity);
}

#[test]
fn g3_planar_hllc_and_hll_are_covariant_under_transverse_boosts() {
    let g = gas(0.0, 0.0);
    for (left, right, fallback) in [
        ([0.9, -0.1, 1.5, 2.0], [0.35, 0.2, -2.0, 0.2], false),
        ([1.0, -2.0, 0.5, 0.4], [1.0, 2.0, -0.7, 0.4], true),
    ] {
        let l = g
            .normal_conserved(left[0], left[1], left[2], left[3])
            .unwrap();
        let r = g
            .normal_conserved(right[0], right[1], right[2], right[3])
            .unwrap();
        let f = g.hllc_normal_flux(l, r, || false).unwrap();
        assert_eq!(f.used_hll_fallback, fallback);
        for boost in [-3.0, 2.0] {
            let shift = |q: NormalConserved| NormalConserved {
                transverse_momentum: q.transverse_momentum + boost * q.density,
                energy: q.energy + boost * q.transverse_momentum + 0.5 * boost * boost * q.density,
                ..q
            };
            let shifted = g.hllc_normal_flux(shift(l), shift(r), || false).unwrap();
            assert_eq!(shifted.used_hll_fallback, fallback);
            close(shifted.signal_speed, f.signal_speed);
            close(shifted.mass, f.mass);
            close(shifted.normal_momentum, f.normal_momentum);
            close(
                shifted.transverse_momentum,
                f.transverse_momentum + boost * f.mass,
            );
            close(
                shifted.energy,
                f.energy + boost * f.transverse_momentum + 0.5 * boost * boost * f.mass,
            );
        }
    }
}

#[test]
fn g0_planar_hll_accounts_for_shear_heating_in_source_domain_admission() {
    let g = gas(0.0, 0.0);
    // Both traces are cold and admissible. The near-vacuum expansion forces
    // HLL, whose mixed state converts opposing transverse KE to heat. Its
    // temperature exceeds this source's upper limit, even as waves expand.
    let left = g.normal_conserved(1.0, -2.0, 30.0, 0.4).unwrap();
    let right = g.normal_conserved(1.0, 2.0, -30.0, 0.4).unwrap();
    assert_eq!(
        g.hllc_normal_flux(left, right, || false),
        Err(Error::Invalid(
            "HLL intermediate state outside admitted source domain"
        )),
    );
    // Discarding transverse energy would incorrectly admit that same face.
    assert!(
        g.hllc_flux(
            g.conserved(1.0, -2.0, 0.4).unwrap(),
            g.conserved(1.0, 2.0, 0.4).unwrap(),
            || false,
        )
        .unwrap()
        .used_hll_fallback
    );
}

#[test]
fn g0_planar_slip_wall_retains_caloric_shock_pressure_and_zero_work() {
    let g = gas(0.4, 1e18);
    for normal in [-1.0, 1.0] {
        for velocity in [-0.4, 0.0, 2.0] {
            let axial = g
                .normal_conserved(1.2, normal * velocity, 0.0, 1.0)
                .unwrap();
            let shear = g
                .normal_conserved(1.2, normal * velocity, 5.0, 1.0)
                .unwrap();
            let reference = g
                .stationary_wall_normal_flux(axial, normal, || false)
                .unwrap();
            let flux = g
                .stationary_wall_normal_flux(shear, normal, || false)
                .unwrap();
            close(flux.normal_momentum, reference.normal_momentum);
            close(flux.signal_speed, reference.signal_speed);
            assert_eq!(flux.mass, 0.0);
            assert_eq!(flux.transverse_momentum, 0.0);
            assert_eq!(flux.energy, 0.0);
        }
    }
}

#[test]
fn g4_planar_flux_refuses_erased_heat_invalid_states_and_cancellation() {
    let g = gas(0.4, 0.0);
    assert!(g.normal_conserved(1.0, 0.0, 1e12, 1.0).is_err());
    assert!(g.normal_conserved(1.0, 0.0, f64::NAN, 1.0).is_err());
    let q = g.normal_conserved(1.0, 0.2, 0.3, 1.0).unwrap();
    assert!(
        g.normal_primitive(NormalConserved { density: 0.0, ..q }, || false)
            .is_err()
    );
    assert!(
        g.normal_primitive(NormalConserved { energy: -1.0, ..q }, || false)
            .is_err()
    );
    assert!(g.stationary_wall_normal_flux(q, 0.0, || false).is_err());
    assert_eq!(g.normal_primitive(q, || true), Err(Error::Cancelled));
    assert_eq!(g.hllc_normal_flux(q, q, || true), Err(Error::Cancelled));
    assert_eq!(
        g.stationary_wall_normal_flux(q, 1.0, || true),
        Err(Error::Cancelled)
    );
    let mut polls = 0;
    assert_eq!(
        g.hllc_normal_flux(q, q, || {
            polls += 1;
            polls == 4
        }),
        Err(Error::Cancelled)
    );
}

#[test]
fn g3_source_owned_primitive_faces_match_conservative_hllc_hll_and_wall_fluxes() {
    let trace = |a: [f64; 4]| NormalPrimitive {
        density: a[0],
        normal_velocity: a[1],
        transverse_velocity: a[2],
        pressure: a[3],
    };
    let components = |f: NormalFaceFlux| {
        [
            f.mass,
            f.normal_momentum,
            f.transverse_momentum,
            f.energy,
            f.signal_speed,
        ]
    };
    for cp in [0.0, 0.4] {
        let g = gas(cp, 1e18);
        let conservative = |p: NormalPrimitive| {
            g.normal_conserved(
                p.density,
                p.normal_velocity,
                p.transverse_velocity,
                p.pressure,
            )
            .unwrap()
        };
        for (left, right) in [
            ([1.0, 0.0, 3.0, 1.0], [4.0, 0.0, -4.0, 1.0]),
            ([1.0, 0.2, 3.0, 1.0], [0.125, 0.2, -4.0, 0.1]),
            ([1.0, -2.0, 0.5, 0.4], [1.0, 2.0, -0.7, 0.4]),
            ([1.2, -8.0, 0.3, 1.7], [0.8, -10.0, -0.1, 1.3]),
        ] {
            let (left, right) = (trace(left), trace(right));
            let original = g
                .hllc_normal_flux(conservative(left), conservative(right), || false)
                .unwrap();
            let reused = g
                .hllc_normal_flux_from_primitive(left, right, || false)
                .unwrap();
            assert_eq!(original.used_hll_fallback, reused.used_hll_fallback);
            for (actual, expected) in components(reused).into_iter().zip(components(original)) {
                assert!((actual - expected).abs() <= 3e-11 * (1.0 + expected.abs()));
            }
        }
        for normal in [-1.0, 1.0] {
            for velocity in [-0.4, 0.0, 2.0] {
                let p = trace([1.2, normal * velocity, 5.0, 1.0]);
                let original = g
                    .stationary_wall_normal_flux(conservative(p), normal, || false)
                    .unwrap();
                let reused = g
                    .stationary_wall_normal_flux_from_primitive(p, normal, || false)
                    .unwrap();
                for (actual, expected) in components(reused).into_iter().zip(components(original)) {
                    assert!((actual - expected).abs() <= 3e-11 * (1.0 + expected.abs()));
                }
            }
        }
    }
}

#[test]
fn g4_primitive_faces_keep_source_precision_domain_and_cancellation_admission() {
    let g = gas(0.4, 0.0);
    let p = NormalPrimitive {
        density: 1.0,
        normal_velocity: 0.2,
        transverse_velocity: 0.3,
        pressure: 1.0,
    };
    for invalid in [
        NormalPrimitive {
            pressure: 100.0,
            ..p
        },
        NormalPrimitive { pressure: 0.0, ..p },
        NormalPrimitive { density: -1.0, ..p },
        NormalPrimitive {
            normal_velocity: f64::NAN,
            ..p
        },
        NormalPrimitive {
            transverse_velocity: 1e12,
            ..p
        },
    ] {
        assert!(
            g.hllc_normal_flux_from_primitive(invalid, p, || false)
                .is_err()
        );
        assert!(
            g.hllc_normal_flux_from_primitive(p, invalid, || false)
                .is_err()
        );
        assert!(
            g.stationary_wall_normal_flux_from_primitive(invalid, 1.0, || false)
                .is_err()
        );
    }
    for cutoff in [1, 2, 3] {
        let mut polls = 0;
        assert_eq!(
            g.hllc_normal_flux_from_primitive(p, p, || {
                polls += 1;
                polls == cutoff
            }),
            Err(Error::Cancelled)
        );
    }
    for cutoff in [1, 2] {
        let mut polls = 0;
        assert_eq!(
            g.stationary_wall_normal_flux_from_primitive(p, 1.0, || {
                polls += 1;
                polls == cutoff
            }),
            Err(Error::Cancelled)
        );
    }
    let g = gas(0.0, 0.0);
    let left = NormalPrimitive {
        density: 1.0,
        normal_velocity: -2.0,
        transverse_velocity: 30.0,
        pressure: 0.4,
    };
    let right = NormalPrimitive {
        normal_velocity: 2.0,
        transverse_velocity: -30.0,
        ..left
    };
    assert_eq!(
        g.hllc_normal_flux_from_primitive(left, right, || false),
        Err(Error::Invalid(
            "HLL intermediate state outside admitted source domain"
        ))
    );
}
