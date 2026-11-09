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
    law: &OpeningFrequency<'_>,
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
    let n = s + a + law.count();
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

#[test]
fn sealed_compliance_and_explicit_neck_recover_the_helmholtz_transfer() {
    let omega = 1.7;
    let (bulk, volume) = (2.0, 3.0);
    let neck = [C64::new(0.4, -omega * 0.7)];
    let exterior = vec![C64::ZERO, C64::ZERO, C64::ZERO, C64::new(0.6, -omega * 0.2)];
    let projection = vec![vec![1.0]];
    let law = OpeningFrequency {
        structural_modes: 1,
        pressure_projection: &projection,
        exterior: &exterior,
        neck: &neck,
    };
    let injected = C64::new(0.032, 0.01);
    let (_, p, out) = solve(
        &law,
        omega,
        &[3.0],
        &[0.0],
        &[0.0],
        0.0,
        bulk,
        &[volume],
        &[0.0],
        &[C64::ZERO],
        &[injected],
    );
    let z = neck[0] + exterior[3];
    let exact_q = injected / (C64::ONE - C64::new(0.0, omega * volume / bulk) * z);
    close(out.volume_flow_m3_s[0], exact_q, 4e-14);
    close(p[0], z * exact_q, 4e-14);
    close(
        out.cavity_pressure_pa[0] - out.exterior_pressure_pa[0],
        neck[0] * exact_q,
        4e-14,
    );
    let input = 0.5 * (p[0] * injected.conj()).re;
    assert!((input - out.neck_dissipation_w - out.exterior_power_w).abs() < 1e-15);
    assert!(out.neck_dissipation_w > 0.0 && out.exterior_power_w > 0.0);
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
fn complete_mutual_loading_closes_original_structural_cavity_and_opening_work() {
    let omega = 1.9;
    let ws = [2.3, 4.1];
    let wa = [0.0, 3.7];
    let damping = [0.13, 0.07];
    let (eta, bulk) = (0.04, 2.4);
    let lambdas = [1.2, 0.8];
    let coupling = [0.5, -0.2, 0.3, 0.6];
    let projection = vec![vec![1.0, 1.0], vec![0.4, -0.8]];
    let neck = [C64::new(0.3, -omega * 0.8), C64::new(0.4, -omega * 1.1)];
    let exterior = joint_impedance(omega);
    let law = OpeningFrequency {
        structural_modes: 2,
        pressure_projection: &projection,
        exterior: &exterior,
        neck: &neck,
    };
    let force = [C64::new(0.8, 0.3), C64::new(-0.4, 0.5)];
    let source = [C64::new(0.2, -0.1), C64::new(0.1, 0.08)];
    let (b, p, out) = solve(
        &law, omega, &ws, &wa, &damping, eta, bulk, &lambdas, &coupling, &force, &source,
    );
    let velocity: Vec<_> = b.iter().map(|&x| C64::new(0.0, -omega) * x).collect();
    let effort: Vec<C64> = (0..4)
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
    for q in 0..2 {
        acoustic += 0.5 * (p[q] * source[q].conj()).re;
        cavity_heat += 0.5 * eta * lambdas[q] / bulk * wa[q] * wa[q] / omega * p[q].norm_sq();
    }
    for j in 0..2 {
        close(
            out.cavity_pressure_pa[j],
            effort[2 + j] + neck[j] * out.volume_flow_m3_s[j],
            3e-13,
        );
    }
    let external = 0.5
        * (0..4)
            .map(|j| (effort[j] * out.radiation_port_velocities[j].conj()).re)
            .sum::<f64>();
    assert!((external - out.exterior_power_w).abs() < 1e-14);
    assert!(
        (input + acoustic - heat - cavity_heat - out.neck_dissipation_w - external).abs() < 2e-14
    );
    assert!((interface + acoustic - cavity_heat - out.cavity_outward_power_w).abs() < 2e-14);

    // SAME physical pressure basis, rescaled and sign-reversed independently.
    let factors = [-3.0, 2.0];
    let scaled_projection: Vec<_> = (0..2)
        .map(|q| projection[q].iter().map(|x| x * factors[q]).collect())
        .collect();
    let scaled_coupling: Vec<_> = (0..4).map(|i| coupling[i] * factors[i % 2]).collect();
    let scaled_norms = [lambdas[0] * 9.0, lambdas[1] * 4.0];
    let scaled_source = [source[0].scale(-3.0), source[1].scale(2.0)];
    let scaled = OpeningFrequency {
        structural_modes: 2,
        pressure_projection: &scaled_projection,
        exterior: &exterior,
        neck: &neck,
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
        close(new_b[r], b[r], 2e-13);
        close(new_p[r].scale(factors[r]), p[r], 2e-13);
        close(new_out.volume_flow_m3_s[r], out.volume_flow_m3_s[r], 2e-13);
    }
    let mut diagonal = exterior.clone();
    for i in 0..4 {
        for j in 0..4 {
            if i != j {
                diagonal[i * 4 + j] = C64::ZERO;
            }
        }
    }
    let wrong = OpeningFrequency {
        structural_modes: 2,
        pressure_projection: &projection,
        exterior: &diagonal,
        neck: &neck,
    };
    let (_, _, wrong_out) = solve(
        &wrong, omega, &ws, &wa, &damping, eta, bulk, &lambdas, &coupling, &force, &source,
    );
    assert!(
        (wrong_out.volume_flow_m3_s[0] - out.volume_flow_m3_s[0]).abs()
            > 0.05 * out.volume_flow_m3_s[0].abs()
    );
}

#[test]
fn outward_wall_and_neck_flow_cancel_in_the_shared_compact_radiation_field() {
    let omega = 1e-3;
    let area = 0.5;
    // Exact rank-one radiation: P=R|area*v+Q|^2/2. Independent
    // positive source powers would destroy the sealed-volume cancellation.
    let exterior = [
        C64::from_re(1.0),
        C64::from_re(2.0),
        C64::from_re(2.0),
        C64::from_re(4.0),
    ];
    let neck = [C64::new(0.2, -omega * 0.6)];
    let projection = vec![vec![1.0]];
    let law = OpeningFrequency {
        structural_modes: 1,
        pressure_projection: &projection,
        exterior: &exterior,
        neck: &neck,
    };
    let (_, _, out) = solve(
        &law,
        omega,
        &[4.0],
        &[0.0],
        &[0.0],
        0.0,
        2.0,
        &[1.0],
        &[area],
        &[C64::ONE],
        &[C64::ZERO],
    );
    let v = out.radiation_port_velocities[0];
    let q = out.volume_flow_m3_s[0];
    let coherent = 2.0 * (v.scale(area) + q).norm_sq();
    let independent = 0.5 * v.norm_sq() + 2.0 * q.norm_sq();
    assert!((out.exterior_power_w / coherent - 1.0).abs() < 2e-11);
    assert!(out.exterior_power_w / independent < 1e-7);
    assert!(out.neck_dissipation_w > out.exterior_power_w);
}

#[test]
fn truncated_modes_keep_original_exterior_indices_and_all_opening_flows() {
    let omega = 1.7;
    let exterior = joint_impedance(omega);
    let projection = vec![vec![1.0, -0.6], vec![-0.2, 0.8]];
    let neck = [C64::new(0.3, -0.7), C64::new(0.2, -0.9)];
    let law = OpeningFrequency {
        structural_modes: 2,
        pressure_projection: &projection,
        exterior: &exterior,
        neck: &neck,
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
    let reduced = OpeningFrequency {
        structural_modes: 1,
        pressure_projection: &map,
        exterior: &principal,
        neck: &neck,
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
    close(b[0], b2[0], 2e-14);
    close(p[0], p2[0], 2e-14);
    for j in 0..2 {
        close(out.volume_flow_m3_s[j], out2.volume_flow_m3_s[j], 2e-14);
    }
    assert_eq!(flow_delta(Some(&out), Some(&out)).unwrap(), 0.0);
    let mut opposite = out.clone();
    for q in &mut opposite.volume_flow_m3_s {
        *q = -*q;
    }
    assert!((flow_delta(Some(&out), Some(&opposite)).unwrap() - 2.0).abs() < 1e-14);
}
