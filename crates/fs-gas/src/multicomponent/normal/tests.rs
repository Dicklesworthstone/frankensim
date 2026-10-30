use super::*;
use crate::multicomponent::tests::{gas, source};
use fs_thermochem::UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R;

fn mixture(formation: f64) -> MixtureGas {
    gas(vec![
        source("A", R, 3.5, 0.0, formation),
        source("B", 2.0 * R, 3.5, 0.0, -2.0 * formation),
    ])
}
fn state(g: &MixtureGas, d: [f64; 2], u: f64, v: f64, p: f64) -> [f64; 5] {
    let mut q = [0.0; 5];
    g.normal_conserved_into(
        NormalPrimitive {
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
fn face(g: &MixtureGas, l: &[f64; 5], r: &[f64; 5]) -> ([f64; 5], FaceFlux) {
    let mut out = [0.0; 5];
    let report = g
        .hllc_normal_flux(from_slice(l, 2), from_slice(r, 2), &mut out, || false)
        .unwrap();
    (out, report)
}
fn close(a: f64, b: f64) {
    assert!((a - b).abs() <= 2e-10 * (1.0 + b.abs()), "{a:e} != {b:e}");
}

#[test]
fn g0_planar_mixture_retains_exact_original_axial_fluxes() {
    let g = mixture(1e17);
    for (ld, lu, lp, rd, ru, rp) in [
        ([1.0, 0.0], 0.0, 1.0, [0.0, 4.0], 0.0, 1.0),
        ([0.8, 0.2], 0.0, 1.0, [0.03, 0.095], 0.0, 0.1),
        ([1.0, 0.0], -2.0, 0.4, [0.0, 1.0], 2.0, 0.4),
        ([0.2, 1.0], 10.0, 1.7, [0.4, 0.4], 8.0, 1.3),
    ] {
        let l = g.conserved(&ld, lu, lp).unwrap();
        let r = g.conserved(&rd, ru, rp).unwrap();
        let mut axial = [0.0; 4];
        let receipt = g
            .hllc_flux(l.as_cell(), r.as_cell(), &mut axial, || false)
            .unwrap();
        let mut planar = [0.0; 5];
        let actual = g
            .hllc_normal_flux(l.as_cell().into(), r.as_cell().into(), &mut planar, || {
                false
            })
            .unwrap();
        assert_eq!(actual, receipt);
        assert_eq!(&planar[..3], &axial[..3]);
        assert_eq!(planar[3], 0.0);
        assert_eq!(planar[4], axial[3]);
    }
}

#[test]
fn g1_planar_material_contact_transports_each_species_and_full_energy() {
    let g = mixture(-1e17);
    for u in [-8.0, -0.2, 0.0, 0.7, 10.0] {
        let l = state(&g, [1.0, 0.0], u, 3.0, 1.0);
        let r = state(&g, [0.05, 0.25], u, -4.0, 1.0);
        let (f, report) = face(&g, &l, &r);
        assert!(!report.used_hll_fallback);
        let q = if u >= 0.0 { l } else { r };
        close(f[0], q[0] * u);
        close(f[1], q[1] * u);
        close(f[2], q[2] * u + 1.0);
        close(f[3], q[3] * u);
        close(f[4], (q[4] + 1.0) * u);
        let p = g.normal_primitive(from_slice(&q, 2), || false).unwrap();
        close(p.pressure, 1.0);
        close(p.temperature, 1.0 / (q[0] + 0.5 * q[1]));
        close(p.sound_speed.powi(2), 1.4 / (q[0] + q[1]));
    }
}

#[test]
fn g3_transverse_boost_covariance_includes_hll_shear_energy() {
    let g = mixture(0.0);
    for (l, r, fallback) in [
        (
            state(&g, [0.8, 0.1], -0.1, 1.5, 2.0),
            state(&g, [0.1, 0.25], 0.2, -2.0, 0.2),
            false,
        ),
        (
            state(&g, [1.0, 0.0], -2.0, 0.5, 0.4),
            state(&g, [0.0, 1.0], 2.0, -0.7, 0.4),
            true,
        ),
    ] {
        let (f, report) = face(&g, &l, &r);
        assert_eq!(report.used_hll_fallback, fallback);
        for b in [-3.0, 2.0] {
            let boost = |mut q: [f64; 5]| {
                let rho = q[0] + q[1];
                q[4] += b * q[3] + 0.5 * b * b * rho;
                q[3] += b * rho;
                q
            };
            let (shifted, actual) = face(&g, &boost(l), &boost(r));
            assert_eq!(actual.used_hll_fallback, fallback);
            close(actual.signal_speed, report.signal_speed);
            for k in 0..3 {
                close(shifted[k], f[k]);
            }
            close(shifted[3], f[3] + b * (f[0] + f[1]));
            close(shifted[4], f[4] + b * f[3] + 0.5 * b * b * (f[0] + f[1]));
        }
    }
    // Original formation references are transported per species, while the
    // conservative sensible-energy flux itself is unchanged by their magnitude.
    let l = state(&g, [0.6, 0.4], 0.3, -0.4, 2.0);
    let r = state(&g, [0.0, 0.5], -0.1, 0.7, 0.7);
    let (f, _) = face(&g, &l, &r);
    let shifted = mixture(1e17);
    let (f2, _) = face(&shifted, &l, &r);
    assert_eq!(f, f2);
    let physical = |owner: &MixtureGas| {
        f[4] + f[0] * owner.energy_references()[0] + f[1] * owner.energy_references()[1]
    };
    let expected = f[0] * (shifted.energy_references()[0] - g.energy_references()[0])
        + f[1] * (shifted.energy_references()[1] - g.energy_references()[1]);
    close(physical(&shifted) - physical(&g), expected);
}

#[test]
fn g0_full_shear_heating_must_fit_the_original_caloric_domain() {
    let g = mixture(0.0);
    let l = state(&g, [1.0, 0.0], -2.0, 30.0, 0.4);
    let r = state(&g, [0.0, 1.0], 2.0, -30.0, 0.4);
    let mut out = [0.0; 5];
    assert!(matches!(
        g.hllc_normal_flux(from_slice(&l, 2), from_slice(&r, 2), &mut out, || false),
        Err(Error::Invalid(
            "mixture HLL intermediate state outside source domain"
        ))
    ));
    let (_, axial) = face(
        &g,
        &state(&g, [1.0, 0.0], -2.0, 0.0, 0.4),
        &state(&g, [0.0, 1.0], 2.0, 0.0, 0.4),
    );
    assert!(axial.used_hll_fallback);
}

#[test]
fn g0_g4_original_wall_pressure_primitive_parity_and_refusals() {
    let g = mixture(0.0);
    let d = [0.2, 1.0];
    for normal in [-1.0, 1.0] {
        for u in [-0.4, 0.0, 2.0] {
            let l = state(&g, d, normal * u, 5.0, 1.0);
            let q = from_slice(&l, 2);
            let primitive = NormalPrimitive {
                partial_densities: &d,
                normal_velocity: normal * u,
                transverse_velocity: 5.0,
                pressure: 1.0,
            };
            let mut out = [0.0; 5];
            let report = g
                .stationary_wall_normal_flux(q, normal, &mut out, || false)
                .unwrap();
            let mut traced = [0.0; 5];
            let trace_report = g
                .stationary_wall_normal_flux_from_primitive(primitive, normal, &mut traced, || {
                    false
                })
                .unwrap();
            for k in 0..5 {
                close(out[k], traced[k]);
            }
            close(report.signal_speed, trace_report.signal_speed);
            let axial = g.conserved(&d, normal * u, 1.0).unwrap();
            let p = g.primitive(axial.as_cell(), || false).unwrap();
            let mut reference = [0.0; 4];
            g.wall_flux(axial.as_cell(), p, normal, &mut reference, &mut || false)
                .unwrap();
            close(out[2], reference[2]);
            for k in [0, 1, 3, 4] {
                assert_eq!(out[k], 0.0);
            }
            g.hllc_normal_flux_from_primitive(primitive, primitive, &mut traced, || false)
                .unwrap();
            let (expected, _) = face(&g, &l, &l);
            for k in 0..5 {
                close(traced[k], expected[k]);
            }
        }
    }
    let input = NormalPrimitive {
        partial_densities: &d,
        normal_velocity: 0.0,
        transverse_velocity: 1e12,
        pressure: 1.0,
    };
    let mut out = [0.0; 5];
    assert!(g.normal_conserved_into(input, &mut out).is_err());
    assert!(
        g.normal_conserved_into(
            NormalPrimitive {
                transverse_velocity: f64::NAN,
                ..input
            },
            &mut out
        )
        .is_err()
    );
    let l = state(&g, d, 0.2, 0.3, 1.0);
    let q = from_slice(&l, 2);
    assert_eq!(g.normal_primitive(q, || true), Err(Error::Cancelled));
    assert_eq!(
        g.hllc_normal_flux(q, q, &mut out, || true),
        Err(Error::Cancelled)
    );
    assert!(
        g.stationary_wall_normal_flux(q, 0.0, &mut out, || false)
            .is_err()
    );
    assert!(g.hllc_normal_flux(q, q, &mut out[..4], || false).is_err());
    let mut polls = 0;
    assert_eq!(
        g.hllc_normal_flux(q, q, &mut out, || {
            polls += 1;
            polls == 3
        }),
        Err(Error::Cancelled)
    );
    // An error does not mutate either borrowed source trace.
    assert_eq!(l, state(&g, d, 0.2, 0.3, 1.0));
}
