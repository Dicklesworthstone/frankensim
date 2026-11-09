use super::*;
use fs_la::eigen_complex::lu_complex;

fn close(a: C64, b: C64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance * b.abs().max(1e-25),
        "{a:?} != {b:?}"
    );
}

#[allow(clippy::too_many_arguments)]
fn solve(
    law: &DistributedOpeningFrequency<'_>,
    omega: f64,
    ws: &[f64],
    wa: &[f64],
    damping: &[f64],
    eta: f64,
    bulk: f64,
    lambdas: &[f64],
    coupling: &[f64],
    force: &[C64],
    source: &[C64],
) -> (Vec<C64>, Vec<C64>, CavityOpeningResponse) {
    let (s, a) = (ws.len(), wa.len());
    let n = s + a + law.unknown_count();
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
        matrix[(s + q) * n + s + q] = C64::new(wa[q] * wa[q] - omega * omega, -eta * wa[q] * wa[q]);
        let scale = bulk / lambdas[q];
        for r in 0..s {
            matrix[(s + q) * n + r] = C64::from_re(-scale * omega * omega * coupling[r * a + q]);
        }
        rhs[s + q] = C64::new(0.0, -omega * scale) * source[q];
    }
    law.add_operator(omega, bulk, lambdas, &mut matrix, s, a)
        .unwrap();
    lu_complex(&matrix, n).unwrap().solve(&mut rhs);
    let b = rhs[..s].to_vec();
    let p = rhs[s..s + a].to_vec();
    let openings = law.response(omega, &b, &p, &rhs[s + a..]).unwrap();
    (b, p, openings)
}

fn cylinder(omega: f64, length: f64, impedance: f64, attenuation: f64, reference: f64) -> TwoPort {
    let wave = fs_duct::SegmentWave {
        wavenumber: C64::new(omega, attenuation),
        characteristic_impedance: C64::from_re(impedance),
        specific_impedance: C64::from_re(impedance),
        shear_number: f64::INFINITY,
    };
    TwoPort::from_uniform_wave(&wave, 1.0, length, reference).unwrap()
}

#[test]
fn finite_compressibility_retains_quarter_half_and_full_wave_transfer() {
    let projection = vec![vec![1.0]];
    let (bulk, volume, z, load) = (2.0, 0.4, 3.2, 1.3);
    let source = C64::new(0.03, 0.01);
    let exterior = [C64::ZERO, C64::ZERO, C64::ZERO, C64::from_re(load)];
    for phase in [
        0.05,
        core::f64::consts::FRAC_PI_2,
        core::f64::consts::PI,
        core::f64::consts::TAU,
    ] {
        let omega = phase;
        // Independent continuum pressure/flow solution, not an S->Y inverse.
        let (sine, cosine) = (fs_math::det::sin(phase), fs_math::det::cos(phase));
        let axial_a = C64::new(load * cosine, -z * sine);
        let b = C64::new(cosine, -load / z * sine);
        for entrance_mass in [0.0, 0.35] {
            let entrance = C64::new(0.0, -omega * entrance_mass);
            let a = axial_a + entrance * b;
            let qo = source / (b - C64::new(0.0, omega * volume / bulk) * a);
            for reference in [0.3, 2.1, 7.0] {
                let cylinder = cylinder(omega, 1.0, z, 0.0, reference);
                let series =
                    TwoPort::from_abcd([C64::ONE, entrance, C64::ZERO, C64::ONE], reference)
                        .unwrap();
                let ducts = [series.star(cylinder, 1e8).unwrap().0];
                let law = DistributedOpeningFrequency {
                    structural_modes: 1,
                    pressure_projection: &projection,
                    exterior: &exterior,
                    ducts: &ducts,
                    tolerance: 2e-12,
                };
                let (_, p, out) = solve(
                    &law,
                    omega,
                    &[13.0],
                    &[0.0],
                    &[0.0],
                    0.0,
                    bulk,
                    &[volume],
                    &[0.0],
                    &[C64::ZERO],
                    &[source],
                );
                close(p[0], a * qo, 3e-13);
                close(out.volume_flow_m3_s[0], qo, 3e-13);
                close(out.cavity_volume_flow_m3_s[0], b * qo, 3e-13);
                let input = 0.5 * (p[0] * source.conj()).re;
                assert!(
                    (input - out.exterior_power_w - out.neck_dissipation_w).abs() < 1e-12 * input
                );
                assert!(out.neck_dissipation_w.abs() < 1e-12 * input);
                assert_eq!(out.distributed_duct_power.len(), 1);
                close(
                    out.distributed_duct_power[0].inward_volume_flow_m3_s[1],
                    -qo,
                    3e-13,
                );
                if phase == core::f64::consts::FRAC_PI_2 {
                    assert!((out.cavity_volume_flow_m3_s[0] - qo).abs() > 0.5 * qo.abs());
                    let compact_z = C64::new(load, -omega * (z + entrance_mass));
                    let compact_q =
                        source / (C64::ONE - C64::new(0.0, omega * volume / bulk) * compact_z);
                    assert!((compact_q - qo).abs() > 0.25 * qo.abs());
                }
            }
        }
    }
}

#[test]
fn opaque_native_duct_keeps_finite_absorption_without_transmission_division() {
    let omega = 1.7;
    let (bulk, volume, z) = (2.0, 0.4, 3.2);
    let projection = vec![vec![1.0]];
    let exterior = [C64::ZERO, C64::ZERO, C64::ZERO, C64::from_re(1.3)];
    let ducts = [cylinder(omega, 1.0, z, 800.0, 2.1)];
    assert_eq!(ducts[0].scattering[1], C64::ZERO);
    assert!(ducts[0].abcd().is_err());
    let law = DistributedOpeningFrequency {
        structural_modes: 1,
        pressure_projection: &projection,
        exterior: &exterior,
        ducts: &ducts,
        tolerance: 2e-12,
    };
    let source = C64::new(0.03, 0.01);
    let (_, p, out) = solve(
        &law,
        omega,
        &[13.0],
        &[0.0],
        &[0.0],
        0.0,
        bulk,
        &[volume],
        &[0.0],
        &[C64::ZERO],
        &[source],
    );
    let qi = source / (C64::ONE - C64::new(0.0, omega * volume / bulk * z));
    close(p[0], qi.scale(z), 3e-13);
    close(out.cavity_volume_flow_m3_s[0], qi, 3e-13);
    assert!(out.volume_flow_m3_s[0].abs() < 1e-25);
    assert!(out.neck_dissipation_w > 0.0);
    let input = 0.5 * (p[0] * source.conj()).re;
    assert!((input - out.neck_dissipation_w - out.exterior_power_w).abs() < 1e-13 * input);
}

fn joint_impedance(omega: f64) -> Vec<C64> {
    let a = [0.2, 0.5, 1.0, -0.3];
    let b = [1.0, 0.2, -0.4, 0.8];
    let h = [1.0, -0.5, 0.8, 0.3];
    let masses = [0.2, 0.3, 0.6, 0.4];
    (0..16)
        .map(|at| {
            let (i, j) = (at / 4, at % 4);
            C64::new(
                a[i] * a[j] + b[i] * b[j] + if i == j { 0.05 } else { 0.0 },
                -omega * (0.2 * h[i] * h[j] + if i == j { masses[i] } else { 0.0 }),
            )
        })
        .collect()
}

#[test]
fn two_ended_mutual_loading_closes_physical_work_and_preserves_signed_cavity_basis() {
    let omega = 1.9;
    let ws = [2.3, 4.1];
    let wa = [0.0, 3.7];
    let damping = [0.13, 0.07];
    let (eta, bulk) = (0.04, 2.4);
    let lambdas = [1.2, 0.8];
    let coupling = [0.5, -0.2, 0.3, 0.6];
    let projection = vec![vec![1.0, 1.0], vec![0.4, -0.8]];
    let ducts = [
        cylinder(omega, 0.8, 2.1, 0.12, 0.8),
        cylinder(omega, 1.3, 1.7, 0.19, 3.0),
    ];
    let exterior = joint_impedance(omega);
    let law = DistributedOpeningFrequency {
        structural_modes: 2,
        pressure_projection: &projection,
        exterior: &exterior,
        ducts: &ducts,
        tolerance: 2e-12,
    };
    let force = [C64::new(0.8, 0.3), C64::new(-0.4, 0.5)];
    let source = [C64::new(0.2, -0.1), C64::new(0.1, 0.08)];
    let (b, p, out) = solve(
        &law, omega, &ws, &wa, &damping, eta, bulk, &lambdas, &coupling, &force, &source,
    );
    let velocity: Vec<_> = b.iter().map(|&x| C64::new(0.0, -omega) * x).collect();
    let effort: Vec<_> = (0..4)
        .map(|i| {
            (0..4).fold(C64::ZERO, |sum, j| {
                sum + exterior[i * 4 + j] * out.radiation_port_velocities[j]
            })
        })
        .collect();
    let (mut input, mut heat, mut cavity_heat, mut interface) = (0.0, 0.0, 0.0, 0.0);
    for r in 0..2 {
        let cavity_force = (0..2).fold(C64::ZERO, |sum, q| sum + p[q].scale(coupling[r * 2 + q]));
        let dynamic = C64::new(ws[r] * ws[r] - omega * omega, -omega * damping[r]) * b[r];
        close(dynamic - cavity_force + effort[r], force[r], 3e-13);
        input += 0.5 * (force[r] * velocity[r].conj()).re;
        heat += 0.5 * damping[r] * velocity[r].norm_sq();
        interface -= 0.5 * (cavity_force * velocity[r].conj()).re;
    }
    let mut acoustic = 0.0;
    let mut duct_work = 0.0;
    for q in 0..2 {
        acoustic += 0.5 * (p[q] * source[q].conj()).re;
        cavity_heat += 0.5 * eta * lambdas[q] / bulk * wa[q] * wa[q] / omega * p[q].norm_sq();
        duct_work += 0.5
            * (out.cavity_pressure_pa[q] * out.cavity_volume_flow_m3_s[q].conj()
                - out.exterior_pressure_pa[q] * out.volume_flow_m3_s[q].conj())
            .re;
        close(
            out.distributed_duct_power[q].inward_volume_flow_m3_s[0],
            out.cavity_volume_flow_m3_s[q],
            4e-13,
        );
        close(
            out.distributed_duct_power[q].inward_volume_flow_m3_s[1],
            -out.volume_flow_m3_s[q],
            4e-13,
        );
        assert!(out.distributed_duct_power[q].absorbed_power_w > 0.0);
    }
    let external = 0.5
        * (0..4)
            .map(|j| (effort[j] * out.radiation_port_velocities[j].conj()).re)
            .sum::<f64>();
    assert!((external - out.exterior_power_w).abs() < 2e-14);
    assert!((duct_work - out.neck_dissipation_w).abs() < 2e-14);
    assert!((input + acoustic - heat - cavity_heat - duct_work - external).abs() < 3e-14);
    assert!((interface + acoustic - cavity_heat - out.cavity_outward_power_w).abs() < 3e-14);

    let factors = [-3.0, 2.0];
    let scaled_projection: Vec<_> = (0..2)
        .map(|q| projection[q].iter().map(|x| x * factors[q]).collect())
        .collect();
    let scaled_coupling: Vec<_> = (0..4).map(|i| coupling[i] * factors[i % 2]).collect();
    let scaled_norms = [lambdas[0] * 9.0, lambdas[1] * 4.0];
    let scaled_source = [source[0].scale(-3.0), source[1].scale(2.0)];
    let scaled = DistributedOpeningFrequency {
        structural_modes: 2,
        pressure_projection: &scaled_projection,
        exterior: &exterior,
        ducts: &ducts,
        tolerance: 2e-12,
    };
    let (new_b, new_p, new_out) = solve(
        &scaled,
        omega,
        &ws,
        &wa,
        &damping,
        eta,
        bulk,
        &scaled_norms,
        &scaled_coupling,
        &force,
        &scaled_source,
    );
    for r in 0..2 {
        close(new_b[r], b[r], 3e-13);
        close(new_p[r].scale(factors[r]), p[r], 3e-13);
        close(new_out.volume_flow_m3_s[r], out.volume_flow_m3_s[r], 3e-13);
        close(
            new_out.cavity_volume_flow_m3_s[r],
            out.cavity_volume_flow_m3_s[r],
            3e-13,
        );
    }
    let mut diagonal = exterior.clone();
    for i in 0..4 {
        for j in 0..4 {
            if i != j {
                diagonal[i * 4 + j] = C64::ZERO;
            }
        }
    }
    let wrong = DistributedOpeningFrequency {
        structural_modes: 2,
        pressure_projection: &projection,
        exterior: &diagonal,
        ducts: &ducts,
        tolerance: 2e-12,
    };
    let (_, _, wrong_out) = solve(
        &wrong, omega, &ws, &wa, &damping, eta, bulk, &lambdas, &coupling, &force, &source,
    );
    assert!(
        (wrong_out.volume_flow_m3_s[0] - out.volume_flow_m3_s[0]).abs()
            > 0.02 * out.volume_flow_m3_s[0].abs()
    );
}

#[test]
fn both_flows_survive_truncation_at_the_original_exterior_port_indices() {
    let omega = 1.7;
    let exterior = joint_impedance(omega);
    let projection = vec![vec![1.0, -0.6], vec![-0.2, 0.8]];
    let ducts = [
        cylinder(omega, 0.8, 2.1, 0.12, 0.8),
        cylinder(omega, 1.3, 1.7, 0.19, 3.0),
    ];
    let law = DistributedOpeningFrequency {
        structural_modes: 2,
        pressure_projection: &projection,
        exterior: &exterior,
        ducts: &ducts,
        tolerance: 2e-12,
    };
    let (b, p, out) = solve(
        &law,
        omega,
        &[3.0],
        &[0.0],
        &[0.1],
        0.0,
        2.4,
        &[1.2],
        &[0.3],
        &[C64::ONE],
        &[C64::new(0.1, 0.2)],
    );
    assert_eq!(out.radiation_port_velocities.len(), 4);
    assert_eq!(out.radiation_port_velocities[1], C64::ZERO);
    assert_eq!(&out.radiation_port_velocities[2..], out.volume_flow_m3_s);
    let indices = [0, 2, 3];
    let principal: Vec<_> = (0..9)
        .map(|at| exterior[indices[at / 3] * 4 + indices[at % 3]])
        .collect();
    let map = vec![projection[0].clone()];
    let reduced = DistributedOpeningFrequency {
        structural_modes: 1,
        pressure_projection: &map,
        exterior: &principal,
        ducts: &ducts,
        tolerance: 2e-12,
    };
    let (b2, p2, out2) = solve(
        &reduced,
        omega,
        &[3.0],
        &[0.0],
        &[0.1],
        0.0,
        2.4,
        &[1.2],
        &[0.3],
        &[C64::ONE],
        &[C64::new(0.1, 0.2)],
    );
    close(b[0], b2[0], 4e-14);
    close(p[0], p2[0], 4e-14);
    for j in 0..2 {
        close(out.volume_flow_m3_s[j], out2.volume_flow_m3_s[j], 4e-14);
        close(
            out.cavity_volume_flow_m3_s[j],
            out2.cavity_volume_flow_m3_s[j],
            4e-14,
        );
    }
    let mut inner_opposite = out.clone();
    for q in &mut inner_opposite.cavity_volume_flow_m3_s {
        *q = -*q;
    }
    let delta =
        super::super::super::physics::flow_delta(Some(&out), Some(&inner_opposite)).unwrap();
    assert!((delta - 2.0).abs() < 1e-14);
}
