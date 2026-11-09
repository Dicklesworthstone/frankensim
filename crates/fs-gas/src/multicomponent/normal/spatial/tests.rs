use super::*;
use crate::multicomponent::{
    NormalPrimitive,
    tests::{gas, source},
};
use fs_thermochem::UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R;

fn mixture(formation: f64) -> MixtureGas {
    gas(vec![
        source("A", R, 3.5, 0.0, formation),
        source("B", 2.0 * R, 3.5, 0.0, -2.0 * formation),
    ])
}
fn state(g: &MixtureGas, d: [f64; 2], u: f64, v: [f64; 2], p: f64) -> [f64; 6] {
    let mut q = [0.0; 6];
    g.spatial_conserved_into(
        SpatialPrimitive {
            partial_densities: &d,
            normal_velocity: u,
            transverse_velocity: v,
            pressure: p,
        },
        &mut q,
    )
    .unwrap();
    q
}
fn face(g: &MixtureGas, l: &[f64; 6], r: &[f64; 6]) -> ([f64; 6], FaceFlux) {
    let mut out = [0.0; 6];
    let report = g
        .hllc_spatial_flux(from_slice(l, 2), from_slice(r, 2), &mut out, || false)
        .unwrap();
    (out, report)
}
fn close(a: f64, b: f64) {
    assert!((a - b).abs() <= 3e-10 * (1.0 + b.abs()), "{a:e} != {b:e}");
}

#[test]
fn g0_spatial_mixture_reduces_bitwise_to_original_planar_flux() {
    let g = mixture(1e17);
    for (ld, lu, lv, lp, rd, ru, rv, rp) in [
        ([1.0, 0.0], 0.0, 0.0, 1.0, [0.0, 4.0], 0.0, 0.0, 1.0),
        ([0.8, 0.2], 0.0, 0.2, 1.0, [0.03, 0.095], 0.0, -0.7, 0.1),
        ([1.0, 0.0], -2.0, 0.5, 0.4, [0.0, 1.0], 2.0, -0.7, 0.4),
        ([0.2, 1.0], 10.0, -2.0, 1.7, [0.4, 0.4], 8.0, 3.0, 1.3),
    ] {
        let mut l = [0.0; 5];
        let mut r = [0.0; 5];
        for (d, u, v, p, out) in [(&ld, lu, lv, lp, &mut l), (&rd, ru, rv, rp, &mut r)] {
            g.normal_conserved_into(
                NormalPrimitive {
                    partial_densities: d,
                    normal_velocity: u,
                    transverse_velocity: v,
                    pressure: p,
                },
                out,
            )
            .unwrap();
        }
        let sl = state(&g, ld, lu, [lv, 0.0], lp);
        let sr = state(&g, rd, ru, [rv, 0.0], rp);
        assert_eq!(&sl[..4], &l[..4]);
        assert_eq!(sl[5], l[4]);
        assert_eq!(&sr[..4], &r[..4]);
        assert_eq!(sr[5], r[4]);
        let mut planar = [0.0; 5];
        let receipt = g
            .hllc_normal_flux(
                super::super::from_slice(&l, 2),
                super::super::from_slice(&r, 2),
                &mut planar,
                || false,
            )
            .unwrap();
        let (actual, report) = face(&g, &sl, &sr);
        assert_eq!(report, receipt);
        assert_eq!(&actual[..4], &planar[..4]);
        assert_eq!(actual[4], 0.0);
        assert_eq!(actual[5], planar[4]);
    }
}

#[test]
fn g1_three_momenta_and_each_species_share_the_material_contact_flux() {
    let g = mixture(-1e17);
    for u in [-8.0, -0.2, 0.0, 0.7, 10.0] {
        let l = state(&g, [1.0, 0.0], u, [3.0, -1.0], 1.0);
        let r = state(&g, [0.05, 0.25], u, [-4.0, 2.0], 1.0);
        let (f, report) = face(&g, &l, &r);
        assert!(!report.used_hll_fallback);
        let q = if u >= 0.0 { l } else { r };
        for k in [0, 1, 3, 4] {
            close(f[k], q[k] * u);
        }
        close(f[2], q[2] * u + 1.0);
        close(f[5], (q[5] + 1.0) * u);
    }
}

#[test]
fn g3_tangent_rotation_normal_reversal_and_two_component_boost_are_covariant() {
    let g = mixture(0.0);
    for (l, r, fallback) in [
        (
            state(&g, [0.8, 0.1], -0.1, [1.5, -0.7], 2.0),
            state(&g, [0.1, 0.25], 0.2, [-2.0, 1.3], 0.2),
            false,
        ),
        (
            state(&g, [1.0, 0.0], -2.0, [0.5, -0.1], 0.4),
            state(&g, [0.0, 1.0], 2.0, [-0.7, 0.2], 0.4),
            true,
        ),
    ] {
        let (f, report) = face(&g, &l, &r);
        assert_eq!(report.used_hll_fallback, fallback);
        let rotate = |mut q: [f64; 6]| {
            [q[3], q[4]] = [0.6 * q[3] - 0.8 * q[4], 0.8 * q[3] + 0.6 * q[4]];
            q
        };
        let (rotated, actual) = face(&g, &rotate(l), &rotate(r));
        assert_eq!(actual.used_hll_fallback, fallback);
        close(actual.signal_speed, report.signal_speed);
        for (a, b) in rotated.into_iter().zip(rotate(f)) {
            close(a, b);
        }
        let reverse = |mut q: [f64; 6]| {
            q[2] = -q[2];
            q
        };
        let (reversed, actual) = face(&g, &reverse(r), &reverse(l));
        assert_eq!(actual.used_hll_fallback, fallback);
        close(actual.signal_speed, report.signal_speed);
        for k in 0..6 {
            close(reversed[k], if k == 2 { f[k] } else { -f[k] });
        }
        for b in [[-3.0, 2.0], [0.7, -1.2]] {
            let boost = |mut q: [f64; 6]| {
                let rho = q[0] + q[1];
                q[5] += b[0] * q[3] + b[1] * q[4] + 0.5 * (b[0] * b[0] + b[1] * b[1]) * rho;
                q[3] += b[0] * rho;
                q[4] += b[1] * rho;
                q
            };
            let (shifted, actual) = face(&g, &boost(l), &boost(r));
            assert_eq!(actual.used_hll_fallback, fallback);
            close(actual.signal_speed, report.signal_speed);
            for k in 0..3 {
                close(shifted[k], f[k]);
            }
            for k in 0..2 {
                close(shifted[3 + k], f[3 + k] + b[k] * (f[0] + f[1]));
            }
            close(
                shifted[5],
                f[5] + b[0] * f[3]
                    + b[1] * f[4]
                    + 0.5 * (b[0] * b[0] + b[1] * b[1]) * (f[0] + f[1]),
            );
        }
    }
}

#[test]
fn g0_hll_mixing_admits_the_heat_of_both_transverse_components_together() {
    let g = mixture(0.0);
    // Each shear alone fits the source temperature range. Both together do not.
    for v in [[12.0, 0.0], [0.0, 12.0]] {
        let l = state(&g, [1.0, 0.0], -2.0, v, 0.4);
        let r = state(&g, [0.0, 1.0], 2.0, [-v[0], -v[1]], 0.4);
        let (_, report) = face(&g, &l, &r);
        assert!(report.used_hll_fallback);
    }
    let l = state(&g, [1.0, 0.0], -2.0, [12.0, 12.0], 0.4);
    let r = state(&g, [0.0, 1.0], 2.0, [-12.0, -12.0], 0.4);
    let mut out = [0.0; 6];
    assert_eq!(
        g.hllc_spatial_flux(from_slice(&l, 2), from_slice(&r, 2), &mut out, || false),
        Err(Error::Invalid(
            "mixture HLL intermediate state outside source domain"
        ))
    );
}

#[test]
fn g0_spatial_variable_cp_and_original_species_energy_gauges_are_retained() {
    let variable = gas(vec![
        source("A", R, 3.5, 0.01, 1e17),
        source("B", 2.0 * R, 4.0, 0.05, -2e17),
    ]);
    for d in [[0.2, 0.8], [1.0, 0.0], [0.0, 1.0]] {
        let q = state(&variable, d, 2.0, [-3.0, 4.0], 20.0);
        let p = variable
            .spatial_primitive(from_slice(&q, 2), || false)
            .unwrap();
        let c = d[0] + 0.5 * d[1];
        let t = 20.0 / c;
        let cv = d[0] * (2.5 + 0.01 * t) + 0.5 * d[1] * (3.0 + 0.05 * t);
        close(p.temperature, t);
        close(p.pressure, 20.0);
        close(p.sound_speed.powi(2), (1.0 + c / cv) * 20.0);
    }
    let g = mixture(0.0);
    let shifted = mixture(1e17);
    let l = state(&g, [0.6, 0.4], 0.3, [-0.4, 1.1], 2.0);
    let r = state(&g, [0.0, 0.5], -0.1, [0.7, -0.3], 0.7);
    assert_eq!(l, state(&shifted, [0.6, 0.4], 0.3, [-0.4, 1.1], 2.0));
    assert_eq!(face(&g, &l, &r), face(&shifted, &l, &r));
    let (f, _) = face(&g, &l, &r);
    let physical = |owner: &MixtureGas| {
        f[5] + f[0] * owner.energy_references()[0] + f[1] * owner.energy_references()[1]
    };
    let expected = f[0] * (shifted.energy_references()[0] - g.energy_references()[0])
        + f[1] * (shifted.energy_references()[1] - g.energy_references()[1]);
    close(physical(&shifted) - physical(&g), expected);
}

#[test]
fn g0_g4_spatial_wall_uses_original_pressure_and_all_zero_transfer_axes() {
    let g = mixture(0.0);
    let d = [0.2, 1.0];
    for normal in [-1.0, 1.0] {
        for u in [-0.4, 0.0, 2.0] {
            let l = state(&g, d, normal * u, [5.0, -3.0], 1.0);
            let q = from_slice(&l, 2);
            let primitive = SpatialPrimitive {
                partial_densities: &d,
                normal_velocity: normal * u,
                transverse_velocity: [5.0, -3.0],
                pressure: 1.0,
            };
            let mut out = [0.0; 6];
            let report = g
                .stationary_wall_spatial_flux(q, normal, &mut out, || false)
                .unwrap();
            let mut traced = [0.0; 6];
            let trace_report = g
                .stationary_wall_spatial_flux_from_primitive(primitive, normal, &mut traced, || {
                    false
                })
                .unwrap();
            for k in 0..6 {
                close(out[k], traced[k]);
            }
            close(report.signal_speed, trace_report.signal_speed);
            let axial = g.conserved(&d, normal * u, 1.0).unwrap();
            let p = g.primitive(axial.as_cell(), || false).unwrap();
            let mut reference = [0.0; 4];
            g.wall_flux(axial.as_cell(), p, normal, &mut reference, &mut || false)
                .unwrap();
            close(out[2], reference[2]);
            for k in [0, 1, 3, 4, 5] {
                assert_eq!(out[k], 0.0);
            }
            g.hllc_spatial_flux_from_primitive(primitive, primitive, &mut traced, || false)
                .unwrap();
            let (expected, _) = face(&g, &l, &l);
            for k in 0..6 {
                close(traced[k], expected[k]);
            }
        }
    }
}

#[test]
fn g0_g4_spatial_precision_axes_source_and_cancellation_refuse_without_input_change() {
    let g = mixture(0.0);
    let d = [0.2, 1.0];
    let primitive = SpatialPrimitive {
        partial_densities: &d,
        normal_velocity: 0.2,
        transverse_velocity: [0.3, 0.4],
        pressure: 1.0,
    };
    let mut out = [0.0; 6];
    for bad in [[f64::NAN, 0.0], [0.0, f64::INFINITY], [0.0, 1e12]] {
        assert!(
            g.spatial_conserved_into(
                SpatialPrimitive {
                    transverse_velocity: bad,
                    ..primitive
                },
                &mut out
            )
            .is_err()
        );
    }
    let l = state(&g, d, 0.2, [0.3, 0.4], 1.0);
    let q = from_slice(&l, 2);
    assert_eq!(g.spatial_primitive(q, || true), Err(Error::Cancelled));
    assert_eq!(
        g.hllc_spatial_flux(q, q, &mut out, || true),
        Err(Error::Cancelled)
    );
    assert!(
        g.stationary_wall_spatial_flux(q, 0.0, &mut out, || false)
            .is_err()
    );
    assert!(g.hllc_spatial_flux(q, q, &mut out[..5], || false).is_err());
    assert!(g.spatial_conserved_into(primitive, &mut out[..5]).is_err());
    let mut polls = 0;
    assert_eq!(
        g.hllc_spatial_flux(q, q, &mut out, || {
            polls += 1;
            polls == 3
        }),
        Err(Error::Cancelled)
    );
    assert!(
        g.spatial_primitive(
            SpatialCell {
                partial_densities: &[-0.1, 1.0],
                ..q
            },
            || false
        )
        .is_err()
    );
    assert!(
        g.spatial_primitive(
            SpatialCell {
                partial_densities: &[0.0, 0.0],
                ..q
            },
            || false
        )
        .is_err()
    );
    assert_eq!(l, state(&g, d, 0.2, [0.3, 0.4], 1.0));
}

#[test]
fn g0_spatial_face_retains_the_complete_bounded_species_axis() {
    let n = crate::multicomponent::MAX_SPECIES;
    let g = gas((0..n)
        .map(|i| source(&format!("S{i}"), R, 3.5, 0.0, 0.0))
        .collect());
    let densities = vec![1.0 / n as f64; n];
    let input = SpatialPrimitive {
        partial_densities: &densities,
        normal_velocity: 0.5,
        transverse_velocity: [-0.3, 0.7],
        pressure: 1.0,
    };
    let mut q = vec![0.0; n + 4];
    g.spatial_conserved_into(input, &mut q).unwrap();
    let mut out = vec![0.0; n + 4];
    g.hllc_spatial_flux(from_slice(&q, n), from_slice(&q, n), &mut out, || false)
        .unwrap();
    for i in 0..n {
        close(out[i], densities[i] * 0.5);
    }
    close(out[n], 1.25);
    close(out[n + 1], -0.15);
    close(out[n + 2], 0.35);
    close(out[n + 3], (q[n + 3] + 1.0) * 0.5);
    assert!(g.spatial_conserved_into(input, &mut out[..n + 3]).is_err());
    out.push(0.0);
    assert!(g.spatial_conserved_into(input, &mut out).is_err());
}

fn moving_face(g: &MixtureGas, l: &[f64; 6], r: &[f64; 6], velocity: f64) -> ([f64; 6], FaceFlux) {
    let mut out = [0.0; 6];
    let report = g
        .moving_hllc_spatial_flux(
            from_slice(l, 2),
            from_slice(r, 2),
            velocity,
            &mut out,
            || false,
        )
        .unwrap();
    (out, report)
}

#[test]
fn g0_g1_moving_spatial_reduces_to_planar_and_retains_uniform_pressure_work() {
    let g = mixture(1e17);
    for (ld, lu, lv, lp, rd, ru, rv, rp) in [
        ([1.0, 0.0], 0.0, 3.0, 1.0, [0.0, 4.0], 0.0, -4.0, 1.0),
        ([0.8, 0.2], 0.0, 0.2, 1.0, [0.03, 0.095], 0.0, -0.7, 0.1),
        ([1.0, 0.0], -2.0, 0.5, 0.4, [0.0, 1.0], 2.0, -0.7, 0.4),
    ] {
        let mut l = [0.0; 5];
        let mut r = [0.0; 5];
        for (d, u, v, p, out) in [(&ld, lu, lv, lp, &mut l), (&rd, ru, rv, rp, &mut r)] {
            g.normal_conserved_into(
                NormalPrimitive {
                    partial_densities: d,
                    normal_velocity: u,
                    transverse_velocity: v,
                    pressure: p,
                },
                out,
            )
            .unwrap();
        }
        let sl = [l[0], l[1], l[2], l[3], 0.0, l[4]];
        let sr = [r[0], r[1], r[2], r[3], 0.0, r[4]];
        for velocity in [-0.5, 0.0, 0.75] {
            let mut planar = [0.0; 5];
            let expected = g
                .moving_hllc_normal_flux(
                    super::super::from_slice(&l, 2),
                    super::super::from_slice(&r, 2),
                    velocity,
                    &mut planar,
                    || false,
                )
                .unwrap();
            let (actual, report) = moving_face(&g, &sl, &sr, velocity);
            assert_eq!(report, expected);
            assert_eq!(&actual[..4], &planar[..4]);
            assert_eq!(actual[4], 0.0);
            assert_eq!(actual[5], planar[4]);
            if velocity == 0.0 {
                assert_eq!((actual, report), face(&g, &sl, &sr));
            }
        }
    }
    let d = [0.25, 0.75];
    let primitive = SpatialPrimitive {
        partial_densities: &d,
        normal_velocity: 0.5,
        transverse_velocity: [1.25, -0.75],
        pressure: 2.0,
    };
    let q = state(&g, d, 0.5, [1.25, -0.75], 2.0);
    let thermo = g.spatial_primitive(from_slice(&q, 2), || false).unwrap();
    for velocity in [-1.0, 0.0, 0.5, 2.0] {
        let (f, report) = moving_face(&g, &q, &q, velocity);
        let relative = 0.5 - velocity;
        for k in [0, 1, 3, 4] {
            close(f[k], q[k] * relative);
        }
        close(f[2], q[2] * relative + 2.0);
        close(f[5], q[5] * relative + 2.0 * 0.5);
        close(report.signal_speed, relative.abs() + thermo.sound_speed);
        let mut traced = [0.0; 6];
        let trace_report = g
            .moving_hllc_spatial_flux_from_primitive(
                primitive,
                primitive,
                velocity,
                &mut traced,
                || false,
            )
            .unwrap();
        for k in 0..6 {
            close(traced[k], f[k]);
        }
        close(trace_report.signal_speed, report.signal_speed);
        if velocity == 0.0 {
            let mut stationary = [0.0; 6];
            let original = g
                .hllc_spatial_flux_from_primitive(primitive, primitive, &mut stationary, || false)
                .unwrap();
            assert_eq!((traced, trace_report), (stationary, original));
        }
    }
}

#[test]
fn g1_g3_moving_spatial_shear_contact_and_vector_boost_keep_one_species_wave_solution() {
    let g = mixture(0.0);
    for velocity in [-2.0, 0.75] {
        let l = state(&g, [1.0, 0.0], velocity, [3.0, -1.0], 1.0);
        let r = state(&g, [0.125, 0.25], velocity, [-4.0, 2.0], 1.0);
        let (f, report) = moving_face(&g, &l, &r, velocity);
        assert!(!report.used_hll_fallback);
        for k in [0, 1, 3, 4] {
            close(f[k], 0.0);
        }
        close(f[2], 1.0);
        close(f[5], velocity);
    }
    let velocity = 0.5;
    for (l, r, fallback) in [
        (
            state(&g, [0.8, 0.1], 0.4, [1.5, -0.7], 2.0),
            state(&g, [0.1, 0.25], 0.7, [-2.0, 1.3], 0.2),
            false,
        ),
        (
            state(&g, [1.0, 0.0], -1.5, [0.5, -0.1], 0.4),
            state(&g, [0.0, 1.0], 2.5, [-0.7, 0.2], 0.4),
            true,
        ),
    ] {
        let (f, report) = moving_face(&g, &l, &r, velocity);
        assert_eq!(report.used_hll_fallback, fallback);
        assert_eq!((f, report), moving_face(&mixture(1e17), &l, &r, velocity));
        for b in [[-3.0, 2.0, -1.0], [0.75, -1.25, 0.5]] {
            let boost = |mut q: [f64; 6]| {
                let rho = q[0] + q[1];
                q[5] += b[0] * q[2]
                    + b[1] * q[3]
                    + b[2] * q[4]
                    + 0.5 * b.iter().map(|v| v * v).sum::<f64>() * rho;
                for k in 0..3 {
                    q[2 + k] += b[k] * rho;
                }
                q
            };
            let (shifted, actual) = moving_face(&g, &boost(l), &boost(r), velocity + b[0]);
            assert_eq!(actual.used_hll_fallback, fallback);
            close(actual.signal_speed, report.signal_speed);
            close(shifted[0], f[0]);
            close(shifted[1], f[1]);
            let mass = f[0] + f[1];
            for k in 0..3 {
                close(shifted[2 + k], f[2 + k] + b[k] * mass);
            }
            close(
                shifted[5],
                f[5] + b[0] * f[2]
                    + b[1] * f[3]
                    + b[2] * f[4]
                    + 0.5 * b.iter().map(|v| v * v).sum::<f64>() * mass,
            );
        }
        let reverse = |mut q: [f64; 6]| {
            q[2] = -q[2];
            q
        };
        let (reversed, actual) = moving_face(&g, &reverse(r), &reverse(l), -velocity);
        assert_eq!(actual.used_hll_fallback, fallback);
        close(actual.signal_speed, report.signal_speed);
        for k in 0..6 {
            close(reversed[k], if k == 2 { f[k] } else { -f[k] });
        }
    }
}

#[test]
fn g0_g4_moving_faces_preserve_source_domain_dimension_and_cancellation_refusals() {
    let g = mixture(0.0);
    let d = [0.2, 1.0];
    let planar = NormalPrimitive {
        partial_densities: &d,
        normal_velocity: 0.2,
        transverse_velocity: 0.3,
        pressure: 1.0,
    };
    let spatial = SpatialPrimitive {
        partial_densities: &d,
        normal_velocity: 0.2,
        transverse_velocity: [0.3, 0.4],
        pressure: 1.0,
    };
    let mut l = [0.0; 5];
    g.normal_conserved_into(planar, &mut l).unwrap();
    let q = super::super::from_slice(&l, 2);
    let s = state(&g, d, 0.2, [0.3, 0.4], 1.0);
    let qs = from_slice(&s, 2);
    let mut out = [0.0; 5];
    let mut sout = [0.0; 6];
    for bad in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY, 1e12] {
        assert!(
            g.moving_hllc_normal_flux(q, q, bad, &mut out, || false)
                .is_err()
        );
        assert!(
            g.moving_hllc_spatial_flux(qs, qs, bad, &mut sout, || false)
                .is_err()
        );
        assert!(
            g.moving_hllc_normal_flux_from_primitive(planar, planar, bad, &mut out, || false,)
                .is_err()
        );
        assert!(
            g.moving_hllc_spatial_flux_from_primitive(spatial, spatial, bad, &mut sout, || false,)
                .is_err()
        );
    }
    assert!(
        g.moving_hllc_normal_flux(q, q, 0.5, &mut out[..4], || false)
            .is_err()
    );
    assert!(
        g.moving_hllc_spatial_flux(qs, qs, 0.5, &mut sout[..5], || false)
            .is_err()
    );
    for bad in [[-0.1, 1.0], [0.0, 0.0]] {
        assert!(
            g.moving_hllc_normal_flux_from_primitive(
                NormalPrimitive {
                    partial_densities: &bad,
                    ..planar
                },
                planar,
                0.5,
                &mut out,
                || false,
            )
            .is_err()
        );
        assert!(
            g.moving_hllc_spatial_flux_from_primitive(
                SpatialPrimitive {
                    partial_densities: &bad,
                    ..spatial
                },
                spatial,
                0.5,
                &mut sout,
                || false,
            )
            .is_err()
        );
    }
    assert_eq!(
        g.moving_hllc_normal_flux(q, q, 0.5, &mut out, || true),
        Err(Error::Cancelled)
    );
    assert_eq!(
        g.moving_hllc_spatial_flux(qs, qs, 0.5, &mut sout, || true),
        Err(Error::Cancelled)
    );
    assert_eq!(
        g.moving_hllc_normal_flux_from_primitive(planar, planar, 0.5, &mut out, || true,),
        Err(Error::Cancelled)
    );
    assert_eq!(
        g.moving_hllc_spatial_flux_from_primitive(spatial, spatial, 0.5, &mut sout, || true,),
        Err(Error::Cancelled)
    );
    let mut polls = 0;
    assert_eq!(
        g.moving_hllc_spatial_flux_from_primitive(spatial, spatial, 0.5, &mut sout, || {
            polls += 1;
            polls == 3
        },),
        Err(Error::Cancelled)
    );

    // Relative-frame conversion must not bypass source-limited HLL shear heating.
    let left = state(&g, [1.0, 0.0], -1.5, [12.0, 12.0], 0.4);
    let right = state(&g, [0.0, 1.0], 2.5, [-12.0, -12.0], 0.4);
    assert_eq!(
        g.moving_hllc_spatial_flux(
            from_slice(&left, 2),
            from_slice(&right, 2),
            0.5,
            &mut sout,
            || false,
        ),
        Err(Error::Invalid(
            "mixture HLL intermediate state outside source domain"
        ))
    );
    for v in [[12.0, 0.0], [0.0, 12.0]] {
        let left = state(&g, [1.0, 0.0], -1.5, v, 0.4);
        let right = state(&g, [0.0, 1.0], 2.5, [-v[0], -v[1]], 0.4);
        assert!(moving_face(&g, &left, &right, 0.5).1.used_hll_fallback);
    }
}
