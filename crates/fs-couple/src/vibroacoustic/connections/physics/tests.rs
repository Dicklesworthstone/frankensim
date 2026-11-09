use super::*;
use crate::vibroacoustic::openings::{
    CavityOpeningResponse, CompactOpeningFrequency, CoupledOpeningFrequency,
    DistributedOpeningFrequency,
};
use fs_la::eigen_complex::lu_complex;

fn close(a: C64, b: C64) {
    assert!((a - b).abs() < 4e-12 * b.abs().max(1e-25), "{a:?} != {b:?}");
}
fn cylinder(omega: f64, length: f64, z: f64, attenuation: f64, reference: f64) -> TwoPort {
    TwoPort::from_uniform_wave(
        &fs_duct::SegmentWave {
            wavenumber: C64::new(omega, attenuation),
            characteristic_impedance: C64::from_re(z),
            specific_impedance: C64::from_re(z),
            shear_number: f64::INFINITY,
        },
        1.0,
        length,
        reference,
    )
    .unwrap()
}
fn link(left: usize, right: usize) -> CavityDuctConnection {
    CavityDuctConnection {
        left_compartment: left,
        right_compartment: right,
        left_pressure_projection: vec![1.0],
        right_pressure_projection: vec![1.0],
        radius_m: 1.0,
        length_m: 1.0,
        end_corrections: [None, None],
    }
}
#[allow(clippy::too_many_arguments)]
fn solve(
    internal: &InternalDuctFrequency<'_>,
    exterior: Option<&CoupledOpeningFrequency<'_>>,
    omega: f64,
    ws: &[f64],
    wa: &[f64],
    damping: &[f64],
    bulk: f64,
    lambdas: &[f64],
    coupling: &[f64],
    force: &[C64],
    source: &[C64],
) -> (
    Vec<C64>,
    Vec<C64>,
    Option<CavityOpeningResponse>,
    InternalCavityDuctResponse,
) {
    let (s, a) = (ws.len(), wa.len());
    let offset = s + a + exterior.map_or(0, |x| x.unknown_count());
    let n = offset + internal.unknown_count();
    let mut matrix = vec![C64::ZERO; n * n];
    let mut rhs = vec![C64::ZERO; n];
    for r in 0..s {
        matrix[r * n + r] = C64::new(ws[r] * ws[r] - omega * omega, -omega * damping[r]);
        rhs[r] = force[r];
        for q in 0..a {
            matrix[r * n + s + q] = C64::from_re(-coupling[r * a + q]);
        }
    }
    for q in 0..a {
        matrix[(s + q) * n + s + q] = C64::from_re(wa[q] * wa[q] - omega * omega);
        let scale = bulk / lambdas[q];
        for r in 0..s {
            matrix[(s + q) * n + r] = C64::from_re(-scale * omega * omega * coupling[r * a + q]);
        }
        rhs[s + q] = C64::new(0.0, -omega * scale) * source[q];
    }
    if let Some(law) = exterior {
        law.add_operator(omega, bulk, lambdas, &mut matrix, s, a, n)
            .unwrap();
    }
    internal
        .add_operator(omega, bulk, lambdas, &mut matrix, s, a, offset, n)
        .unwrap();
    lu_complex(&matrix, n).unwrap().solve(&mut rhs);
    let b = rhs[..s].to_vec();
    let p = rhs[s..s + a].to_vec();
    let external = exterior.map(|law| law.response(omega, &b, &p, &rhs[s + a..offset]).unwrap());
    let inside = internal.response(&p, &rhs[offset..]).unwrap();
    (b, p, external, inside)
}

#[test]
fn finite_duct_between_sealed_compartments_matches_independent_continuum_transfer() {
    let left = vec![vec![1.0], vec![0.0]];
    let right = vec![vec![0.0], vec![1.0]];
    let links = [link(0, 1)];
    let (bulk, cl, cr, z) = (2.0, 0.4, 0.7, 3.1);
    let sources = [C64::new(0.2, -0.1), C64::new(-0.05, 0.03)];
    for omega in [
        0.9,
        core::f64::consts::FRAC_PI_2,
        core::f64::consts::PI,
        core::f64::consts::TAU,
    ] {
        for alpha in [0.0, 0.13] {
            let ducts = [cylinder(omega, 1.0, z, alpha, 0.8)];
            let law = InternalDuctFrequency {
                left_projection: &left,
                right_projection: &right,
                connections: &links,
                compartments: 2,
                ducts: &ducts,
                tolerance: 3e-12,
            };
            let (_, p, _, out) = solve(
                &law,
                None,
                omega,
                &[8.0],
                &[0.0, 0.0],
                &[0.0],
                bulk,
                &[cl * bulk, cr * bulk],
                &[0.0, 0.0],
                &[C64::ZERO],
                &sources,
            );
            // Original axial p/U continuum, evaluated independently of native S.
            let (sin, cos) = (fs_math::det::sin(omega), fs_math::det::cos(omega));
            let e = fs_math::det::exp(alpha);
            let (ch, sh) = (0.5 * (e + 1.0 / e), 0.5 * (e - 1.0 / e));
            let c = C64::new(cos * ch, -sin * sh);
            let s = C64::new(sin * ch, cos * sh);
            let (a, b, cc, d) = (c, C64::new(0.0, -z) * s, C64::new(0.0, -1.0 / z) * s, c);
            let iw = C64::new(0.0, omega);
            let aa = a - iw.scale(cr) * b;
            let bb = cc - iw.scale(cr) * d;
            let pr = (sources[0] + (d - iw.scale(cl) * b) * sources[1]) / (bb - iw.scale(cl) * aa);
            let qr = -sources[1] - iw.scale(cr) * pr;
            let pl = a * pr + b * qr;
            let ql = cc * pr + d * qr;
            close(p[0], pl);
            close(p[1], pr);
            close(out.left_volume_flow_m3_s[0], ql);
            close(out.right_volume_flow_m3_s[0], qr);
            let input = 0.5 * (p[0] * sources[0].conj() + p[1] * sources[1].conj()).re;
            assert!((input - out.dissipation_w).abs() < 3e-13);
            assert!((out.net_cavity_work_w - out.dissipation_w).abs() < 3e-13);
            if alpha > 0.0 {
                assert!(out.dissipation_w > 0.0);
            }
            if omega == core::f64::consts::FRAC_PI_2 {
                assert!((ql - qr).abs() > 0.1 * ql.abs().max(qr.abs()));
            }
        }
    }
}

fn exterior_matrix(omega: f64) -> Vec<C64> {
    let a = [0.3, 1.0, 0.4];
    let b = [1.0, -0.2, 0.8];
    let h = [0.4, -0.5, 0.9];
    (0..9)
        .map(|at| {
            let (i, j) = (at / 3, at % 3);
            C64::new(
                a[i] * a[j] + b[i] * b[j] + if i == j { 0.1 } else { 0.0 },
                -omega * (0.1 * h[i] * h[j] + if i == j { 0.2 } else { 0.0 }),
            )
        })
        .collect()
}
#[test]
fn internal_graph_and_joint_exterior_close_work_and_preserve_orientation_and_basis() {
    let omega = 1.7;
    let left = vec![vec![1.0], vec![0.0]];
    let right = vec![vec![0.0], vec![1.0]];
    let links = [link(0, 1)];
    let ducts = [cylinder(omega, 1.2, 1.8, 0.1, 2.7)];
    let internal = InternalDuctFrequency {
        left_projection: &left,
        right_projection: &right,
        connections: &links,
        compartments: 2,
        ducts: &ducts,
        tolerance: 3e-12,
    };
    let h = vec![vec![1.0, 0.0], vec![0.0, 1.0]];
    let extducts = [
        cylinder(omega, 0.6, 2.1, 0.07, 0.8),
        cylinder(omega, 0.9, 1.7, 0.11, 2.3),
    ];
    let z = exterior_matrix(omega);
    let exterior = CoupledOpeningFrequency::Distributed(DistributedOpeningFrequency {
        structural_modes: 1,
        pressure_projection: &h,
        exterior: &z,
        ducts: &extducts,
        tolerance: 3e-12,
    });
    let force = [C64::new(0.7, 0.3)];
    let source = [C64::new(0.2, -0.1), C64::new(-0.03, 0.04)];
    let coupling = [0.4, -0.3];
    let (b, p, ext, out) = solve(
        &internal,
        Some(&exterior),
        omega,
        &[3.2],
        &[0.0, 0.0],
        &[0.11],
        2.4,
        &[1.2, 0.8],
        &coupling,
        &force,
        &source,
    );
    let ext = ext.unwrap();
    let v = C64::new(0.0, -omega) * b[0];
    let pressure_force = p[0].scale(coupling[0]) + p[1].scale(coupling[1]);
    let effort = (0..3).fold(C64::ZERO, |sum, j| {
        sum + z[j] * ext.radiation_port_velocities[j]
    });
    close(
        C64::new(3.2 * 3.2 - omega * omega, -omega * 0.11) * b[0] - pressure_force + effort,
        force[0],
    );
    let acoustic = 0.5 * (p[0] * source[0].conj() + p[1] * source[1].conj()).re;
    let input = 0.5 * (force[0] * v.conj()).re + acoustic;
    let structural = 0.5 * 0.11 * v.norm_sq();
    assert!(
        (input - structural - ext.neck_dissipation_w - ext.exterior_power_w - out.dissipation_w)
            .abs()
            < 3e-14
    );
    let interface = -0.5 * (pressure_force * v.conj()).re;
    assert!(
        (interface + acoustic - ext.cavity_outward_power_w - out.net_cavity_work_w).abs() < 3e-14
    );
    assert!(
        (out.compartment_outward_power_w.iter().sum::<f64>() - out.dissipation_w).abs() < 3e-14
    );
    assert!(out.dissipation_w > 0.0 && ext.neck_dissipation_w > 0.0 && ext.exterior_power_w > 0.0);

    let compact_necks = [C64::new(0.2, -omega * 0.4), C64::new(0.3, -omega * 0.6)];
    let compact = CoupledOpeningFrequency::Compact(CompactOpeningFrequency {
        structural_modes: 1,
        pressure_projection: &h,
        exterior: &z,
        neck: &compact_necks,
    });
    let (cb, cp, ce, ci) = solve(
        &internal,
        Some(&compact),
        omega,
        &[3.2],
        &[0.0, 0.0],
        &[0.11],
        2.4,
        &[1.2, 0.8],
        &coupling,
        &force,
        &source,
    );
    let ce = ce.unwrap();
    let cv = C64::new(0.0, -omega) * cb[0];
    let cin = 0.5 * (force[0] * cv.conj() + cp[0] * source[0].conj() + cp[1] * source[1].conj()).re;
    assert!(
        (cin - 0.5 * 0.11 * cv.norm_sq()
            - ce.neck_dissipation_w
            - ce.exterior_power_w
            - ci.dissipation_w)
            .abs()
            < 3e-14
    );
    assert_eq!(ce.cavity_volume_flow_m3_s, ce.volume_flow_m3_s);

    // Reverse the physical duct, not either pressure basis or source phase.
    let s = ducts[0].scattering;
    let reversed = [TwoPort::new([s[3], s[2], s[1], s[0]], ducts[0].reference_impedance).unwrap()];
    let reverse_links = [link(1, 0)];
    let reverse = InternalDuctFrequency {
        left_projection: &right,
        right_projection: &left,
        connections: &reverse_links,
        compartments: 2,
        ducts: &reversed,
        tolerance: 3e-12,
    };
    let (rb, rp, re, ri) = solve(
        &reverse,
        Some(&exterior),
        omega,
        &[3.2],
        &[0.0, 0.0],
        &[0.11],
        2.4,
        &[1.2, 0.8],
        &coupling,
        &force,
        &source,
    );
    close(rb[0], b[0]);
    close(rp[0], p[0]);
    close(rp[1], p[1]);
    close(ri.left_volume_flow_m3_s[0], -out.right_volume_flow_m3_s[0]);
    close(ri.right_volume_flow_m3_s[0], -out.left_volume_flow_m3_s[0]);
    close(re.unwrap().volume_flow_m3_s[0], ext.volume_flow_m3_s[0]);

    let sl = vec![vec![-3.0], vec![0.0]];
    let sr = vec![vec![0.0], vec![2.0]];
    let sh = vec![vec![-3.0, 0.0], vec![0.0, 2.0]];
    let scaled = InternalDuctFrequency {
        left_projection: &sl,
        right_projection: &sr,
        connections: &links,
        compartments: 2,
        ducts: &ducts,
        tolerance: 3e-12,
    };
    let scaled_ext = CoupledOpeningFrequency::Distributed(DistributedOpeningFrequency {
        structural_modes: 1,
        pressure_projection: &sh,
        exterior: &z,
        ducts: &extducts,
        tolerance: 3e-12,
    });
    let (sb, sp, se, si) = solve(
        &scaled,
        Some(&scaled_ext),
        omega,
        &[3.2],
        &[0.0, 0.0],
        &[0.11],
        2.4,
        &[1.2 * 9.0, 0.8 * 4.0],
        &[-1.2, -0.6],
        &force,
        &[source[0].scale(-3.0), source[1].scale(2.0)],
    );
    close(sb[0], b[0]);
    close(sp[0].scale(-3.0), p[0]);
    close(sp[1].scale(2.0), p[1]);
    close(si.left_volume_flow_m3_s[0], out.left_volume_flow_m3_s[0]);
    close(si.right_volume_flow_m3_s[0], out.right_volume_flow_m3_s[0]);
    close(se.unwrap().volume_flow_m3_s[0], ext.volume_flow_m3_s[0]);
    let mut changed = out.clone();
    changed.right_volume_flow_m3_s[0] = -changed.right_volume_flow_m3_s[0];
    assert!((flow_delta(Some(&out), Some(&changed)).unwrap() - 2.0).abs() < 1e-14);
}
#[test]
fn opaque_internal_duct_absorbs_from_one_compartment_without_forcing_the_other() {
    let omega = 1.3;
    let left = vec![vec![1.0], vec![0.0]];
    let right = vec![vec![0.0], vec![1.0]];
    let links = [link(0, 1)];
    let z = 2.7;
    let ducts = [cylinder(omega, 1.0, z, 800.0, 1.2)];
    assert!(ducts[0].abcd().is_err());
    let law = InternalDuctFrequency {
        left_projection: &left,
        right_projection: &right,
        connections: &links,
        compartments: 2,
        ducts: &ducts,
        tolerance: 3e-12,
    };
    let source = C64::new(0.2, 0.1);
    let (_, p, _, out) = solve(
        &law,
        None,
        omega,
        &[8.0],
        &[0.0, 0.0],
        &[0.0],
        2.0,
        &[0.8, 1.4],
        &[0.0, 0.0],
        &[C64::ZERO],
        &[source, C64::ZERO],
    );
    let q = source / (C64::ONE - C64::new(0.0, omega * 0.4 * z));
    close(out.left_volume_flow_m3_s[0], q);
    close(p[0], q.scale(z));
    assert_eq!(p[1], C64::ZERO);
    assert_eq!(out.right_volume_flow_m3_s[0], C64::ZERO);
    assert_eq!(out.compartment_outward_power_w[1], 0.0);
    assert!((0.5 * (p[0] * source.conj()).re - out.dissipation_w).abs() < 2e-14);
}

#[test]
fn combined_basis_preserves_local_pressure_support_and_all_compression_coordinates() {
    use crate::vibroacoustic::CavityModes;
    use crate::vibroacoustic::connections::basis::CavityCompartmentBasis;
    let first = CavityModes {
        omegas: vec![0.0, 4.0],
        lambdas: vec![1.2, 0.6],
        interface: vec![vec![-2.0, -2.0], vec![0.7, -0.3]],
        loss_factor: 0.0,
        rho0: 1.2,
        c0: 343.0,
    };
    let second = CavityModes {
        omegas: vec![0.0, 3.0],
        lambdas: vec![0.8, 0.4],
        interface: vec![vec![1.0, 1.0, 1.0], vec![0.2, -0.5, 0.7]],
        loss_factor: 0.0,
        rho0: 1.2,
        c0: 343.0,
    };
    let original = [first, second];
    let combined = CavityCompartmentBasis::try_new(&original, 20).unwrap();
    assert!(CavityCompartmentBasis::try_new(&original, 19).is_err());
    assert!(matches!(
        CavityCompartmentBasis::try_new_controlled(&original, 20, &mut || true),
        Err(VibroError::Cancelled)
    ));
    assert_eq!(&combined.modes().omegas[..2], &[0.0, 0.0]);
    assert_eq!(combined.modes().omegas, vec![0.0, 0.0, 3.0, 4.0]);
    let local = [
        [C64::new(0.3, 0.2), C64::new(-0.1, 0.5)],
        [C64::new(-0.4, 0.3), C64::new(0.2, -0.1)],
    ];
    let mut global = vec![C64::ZERO; 4];
    for c in 0..2 {
        for q in 0..2 {
            global[combined.mode_indices(c).unwrap()[q]] = local[c][q];
        }
    }
    let mut reference_norm = 0.0;
    for c in 0..2 {
        let start = combined.interface_ranges()[c].start;
        for i in 0..original[c].interface[0].len() {
            let actual = (0..4).fold(C64::ZERO, |sum, q| {
                sum + global[q].scale(combined.modes().interface[q][start + i])
            });
            let expected = (0..2).fold(C64::ZERO, |sum, q| {
                sum + local[c][q].scale(original[c].interface[q][i])
            });
            close(actual, expected);
        }
        for q in 0..2 {
            reference_norm += original[c].lambdas[q] * local[c][q].norm_sq();
        }
    }
    let combined_norm = (0..4)
        .map(|q| combined.modes().lambdas[q] * global[q].norm_sq())
        .sum::<f64>();
    assert!((combined_norm - reference_norm).abs() < 1e-15);
}
