use super::{Boundary, dot};
use fs_math::c64::C64;

fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() < 4e-13 * (1.0 + expected.abs()),
        "{actual:e} != {expected:e}"
    );
}

#[test]
fn signed_mutual_debye_exchange_matches_original_storage_identity() {
    // Authored passive realization: H_load = z^T R z/2, A=-2I, B=I,
    // C=R=[[3,-1],[-1,2]], D=0. No storage law is inferred in production.
    let boundary = Boundary {
        host: 3,
        states: 2,
        ports: 2,
        columns: vec![1.0, -0.4, 0.3, 0.7, -0.6, 1.2],
        a: vec![-2.0, 0.0, 0.0, -2.0],
        b: vec![1.0, 0.0, 0.0, 1.0],
        c: vec![3.0, -1.0, -1.0, 2.0],
        d: vec![0.0; 4],
    };
    for (row, expected) in [[-3.4, 1.8], [-0.2, -1.1], [3.0, -3.0]].iter().enumerate() {
        let mut reaction = [0.0; 2];
        assert!(boundary.state_reaction_row(row, &mut reaction));
        close(reaction[0], expected[0]);
        close(reaction[1], expected[1]);
    }
    for (h, z) in [
        ([0.2, -0.7, 1.1], [0.5, -0.2]),
        ([-1.0, 0.4, 0.1], [-0.3, 0.8]),
    ] {
        let mut flow = [0.0; 2];
        let mut effort = [0.0; 2];
        let mut rate = [0.0; 2];
        let mut host_rate = [0.0; 3];
        assert!(boundary.flow_into(&h, &mut flow));
        close(flow[0], h[0] + 0.3 * h[1] - 0.6 * h[2]);
        close(flow[1], -0.4 * h[0] + 0.7 * h[1] + 1.2 * h[2]);
        assert!(boundary.effort_into(&z, &flow, &mut effort));
        assert!(boundary.rate_into(&z, &flow, &mut rate));
        let mut direct = [1.0; 2];
        let mut driven = [0.0; 2];
        assert!(boundary.direct_effort_into(&flow, &mut direct));
        assert!(boundary.driven_rate_into(&flow, &mut driven));
        assert_eq!(direct, [0.0; 2]);
        assert_eq!(driven, flow);
        for (row, value) in host_rate.iter_mut().enumerate() {
            let mut reaction = [0.0; 2];
            assert!(boundary.state_reaction_row(row, &mut reaction));
            *value = dot(&reaction, &z);
        }
        close(effort[0], 3.0 * z[0] - z[1]);
        close(effort[1], -z[0] + 2.0 * z[1]);
        close(rate[0], -2.0 * z[0] + flow[0]);
        close(rate[1], -2.0 * z[1] + flow[1]);
        let exchanged = dot(&flow, &effort);
        close(boundary.host_power(&h, &effort), exchanged);
        close(dot(&h, &host_rate), -exchanged);
        let storage_rate = dot(&effort, &rate);
        let heat = 2.0 * (3.0 * z[0] * z[0] - 2.0 * z[0] * z[1] + 2.0 * z[1] * z[1]);
        close(exchanged - storage_rate, heat);
        assert!(heat > 0.0);
    }
}

#[test]
fn original_pair_chart_matches_independent_complex_residue_impedance() {
    let pole = C64::new(-0.4, 1.8);
    let residues = [
        C64::new(0.4, 0.1),
        C64::new(-0.1, 0.03),
        C64::new(-0.1, 0.03),
        C64::new(0.3, -0.05),
    ];
    let direct = [0.9, 0.15, 0.15, 0.7];
    let boundary = Boundary {
        host: 2,
        states: 4,
        ports: 2,
        columns: vec![1.0, -0.4, 0.3, 0.8],
        a: vec![
            -0.4, 0.0, 1.8, 0.0, 0.0, -0.4, 0.0, 1.8, -1.8, 0.0, -0.4, 0.0, 0.0, -1.8, 0.0, -0.4,
        ],
        b: vec![1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0],
        c: vec![0.8, -0.2, 0.2, 0.06, -0.2, 0.6, 0.06, -0.1],
        d: direct.to_vec(),
    };
    for omega in [0.0, 0.2, 1.8, 6.0] {
        let s = C64::new(0.0, omega);
        let h = [C64::new(0.4, -0.2), C64::new(-0.3, 0.1)];
        let y = [h[0] + h[1].scale(0.3), h[0].scale(-0.4) + h[1].scale(0.8)];
        let den = (s - C64::from_re(pole.re)) * (s - C64::from_re(pole.re))
            + C64::from_re(pole.im * pole.im);
        let z = [
            y[0] * (s - C64::from_re(pole.re)) / den,
            y[1] * (s - C64::from_re(pole.re)) / den,
            y[0].scale(-pole.im) / den,
            y[1].scale(-pole.im) / den,
        ];
        for imaginary in [false, true] {
            let part = |z: C64| if imaginary { z.im } else { z.re };
            let hs = h.map(part);
            let zs = z.map(part);
            let mut flow = [0.0; 2];
            let mut effort = [0.0; 2];
            let mut rate = [0.0; 4];
            assert!(boundary.flow_into(&hs, &mut flow));
            assert!(boundary.effort_into(&zs, &flow, &mut effort));
            assert!(boundary.rate_into(&zs, &flow, &mut rate));
            for row in 0..2 {
                let expected = (0..2).fold(C64::ZERO, |sum, col| {
                    let r = residues[row * 2 + col];
                    sum + (r / (s - pole)
                        + r.conj() / (s - pole.conj())
                        + C64::from_re(direct[row * 2 + col]))
                        * y[col]
                });
                close(flow[row], part(y[row]));
                close(effort[row], part(expected));
            }
            for i in 0..4 {
                close(rate[i], part(s * z[i]));
            }
        }
    }
}

#[test]
fn cancellation_keeps_small_real_flow_and_work_without_an_absolute_floor() {
    let boundary = Boundary {
        host: 3,
        states: 0,
        ports: 1,
        columns: vec![1e16, 1.0, -1e16],
        a: vec![],
        b: vec![],
        c: vec![],
        d: vec![2.0],
    };
    let h = [1.0; 3];
    let mut flow = [0.0];
    let mut effort = [0.0];
    assert!(boundary.flow_into(&h, &mut flow));
    assert_eq!(flow, [1.0]);
    assert!(boundary.effort_into(&[], &flow, &mut effort));
    assert!(boundary.rate_into(&[], &flow, &mut []));
    assert_eq!(effort, [2.0]);
    assert_eq!(boundary.host_power(&h, &effort), 2.0);
}

#[test]
fn original_dark_state_evolves_and_stored_history_can_return_work() {
    let boundary = Boundary {
        host: 1,
        states: 2,
        ports: 1,
        columns: vec![1.0],
        a: vec![-2.0, 0.0, 0.0, -3.0],
        b: vec![1.0, 1.0],
        c: vec![2.0, 0.0],
        d: vec![0.0],
    };
    let mut effort = [0.0];
    let mut rate = [0.0; 2];
    assert!(boundary.effort_into(&[-1.0, 5.0], &[0.4], &mut effort));
    assert!(boundary.rate_into(&[-1.0, 5.0], &[0.4], &mut rate));
    assert_eq!(effort, [-2.0]);
    close(rate[0], 2.4);
    close(rate[1], -14.6);
    assert_eq!(dot(&[0.4], &effort), -0.8);
    assert!(!boundary.flow_into(&[], &mut [0.0]));
    assert!(!boundary.rate_into(&[0.0], &[0.0], &mut [0.0; 2]));
}
