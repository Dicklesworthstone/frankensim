//! G1/G3 physical load shifts, mutual impedance, port work and unit changes.
use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_math::c64::C64;
use fs_modal::harmonic::{
    HarmonicPortControls, HarmonicResponseControls, HarmonicResponseError, PreparedHarmonicResponse,
};

fn with_cx<R>(gate: &CancelGate, f: impl FnOnce(&Cx<'_>) -> R) -> R {
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        f(&Cx::new(
            gate,
            arena,
            StreamKey {
                seed: 91,
                kernel_id: 22,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        ))
    })
}
fn controls() -> HarmonicResponseControls {
    HarmonicResponseControls {
        maximum_coordinates: 8,
        maximum_payload_bytes: 1 << 20,
        maximum_relative_residual: 1e-12,
    }
}
fn ports() -> HarmonicPortControls {
    HarmonicPortControls {
        maximum_ports: 8,
        maximum_projection_terms: 10_000,
        maximum_payload_bytes: 1 << 20,
    }
}
fn close(a: C64, b: C64, tol: f64) {
    assert!((a - b).abs() <= tol * b.abs().max(1e-250), "{a:?} != {b:?}");
}
fn scalar(a: f64, b: f64, tol: f64) {
    close(C64::from_re(a), C64::from_re(b), tol);
}

#[test]
fn g1_physical_load_regularizes_bare_pole_and_shifts_the_resonance() {
    with_cx(&CancelGate::new(), |cx| {
        assert!(matches!(
            PreparedHarmonicResponse::try_new(&[2.0], &[0.0], &[8.0], 1, 2.0, controls(), cx),
            Err(HarmonicResponseError::Singular { .. })
        ));
        let force = C64::new(3.0, -0.7);
        // A port transformer B=3 exposes an original dashpot, inertance and
        // spring. Their physical contributions are B^2 times their own data.
        let (load_c, load_m, load_k) = (0.4, 0.7, 0.2);
        let loaded_natural = ((8.0_f64 + 9.0 * load_k) / (2.0 + 9.0 * load_m)).sqrt();
        for w in [0.01, 0.7, loaded_natural, 2.0, 4.0, 100.0] {
            let z = C64::new(load_c, load_k / w - w * load_m);
            let response = PreparedHarmonicResponse::try_with_port_impedance(
                &[2.0],
                &[0.0],
                &[8.0],
                1,
                &[3.0],
                1,
                &[z],
                w,
                controls(),
                ports(),
                cx,
            )
            .unwrap()
            .solve(&[force], cx)
            .unwrap();
            let expected = force
                / C64::new(
                    8.0 + 9.0 * load_k - (2.0 + 9.0 * load_m) * w * w,
                    -9.0 * load_c * w,
                );
            close(response.displacement[0], expected, 1e-13);
            let port = response.port_load.as_ref().unwrap();
            close(port.flow[0], response.velocity[0].scale(3.0), 1e-14);
            close(port.effort[0], z * port.flow[0], 1e-14);
            scalar(
                port.average_absorbed_power_w,
                0.5 * load_c * port.flow[0].norm_sq(),
                2e-13,
            );
            assert_eq!(response.average_viscous_power_w, 0.0);
            assert!(
                response.power_balance_residual_w.abs() < 2e-12 * port.average_absorbed_power_w
            );
            if w == loaded_natural {
                assert!((response.velocity[0] / force).im.abs() < 2e-15);
            }
        }
        // An arbitrarily large but finite velocity impedance does no DC work.
        let dc = PreparedHarmonicResponse::try_with_port_impedance(
            &[2.0],
            &[0.0],
            &[8.0],
            1,
            &[3.0],
            1,
            &[C64::from_re(1e250)],
            0.0,
            controls(),
            ports(),
            cx,
        )
        .unwrap()
        .solve(&[force], cx)
        .unwrap();
        close(dc.displacement[0], force.scale(1.0 / 8.0), 1e-15);
        assert_eq!(dc.port_load.unwrap().average_absorbed_power_w, 0.0);
    });
}

const MASS: [f64; 4] = [2.0, 0.15, 0.15, 1.1];
const DAMPING: [f64; 4] = [0.3, 0.04, 0.04, 0.2];
const STIFFNESS: [f64; 4] = [9.0, -0.8, -0.8, 5.0];
const MAP: [f64; 4] = [1.0, 0.2, -0.3, 1.5];
const LOAD: [C64; 4] = [
    C64::new(2.0, 0.4),
    C64::new(0.3, 0.2),
    C64::new(-0.1, 0.1),
    C64::new(1.0, -0.2),
];
const FORCE: [C64; 2] = [C64::new(0.8, -0.1), C64::new(-0.3, 0.7)];

#[test]
fn g1_complete_mutual_nonreciprocal_impedance_matches_independent_physical_rows() {
    with_cx(&CancelGate::new(), |cx| {
        for w in [0.1, 1.5, 2.1, 5.0] {
            let prepared = PreparedHarmonicResponse::try_with_port_impedance(
                &MASS,
                &DAMPING,
                &STIFFNESS,
                2,
                &MAP,
                2,
                &LOAD,
                w,
                controls(),
                ports(),
                cx,
            )
            .unwrap();
            let response = prepared.solve(&FORCE, cx).unwrap();
            let mut a = [C64::ZERO; 4];
            for i in 0..2 {
                for j in 0..2 {
                    a[i * 2 + j] = C64::new(
                        STIFFNESS[i * 2 + j] - w * w * MASS[i * 2 + j],
                        -w * DAMPING[i * 2 + j],
                    );
                    // Independent four-index physical contraction and Cramer's rule.
                    for p in 0..2 {
                        for q in 0..2 {
                            a[i * 2 + j] = a[i * 2 + j]
                                + C64::new(0.0, -w)
                                    * LOAD[p * 2 + q].scale(MAP[i * 2 + p] * MAP[j * 2 + q]);
                        }
                    }
                }
            }
            let determinant = a[0] * a[3] - a[1] * a[2];
            close(
                response.displacement[0],
                (FORCE[0] * a[3] - a[1] * FORCE[1]) / determinant,
                2e-13,
            );
            close(
                response.displacement[1],
                (a[0] * FORCE[1] - FORCE[0] * a[2]) / determinant,
                2e-13,
            );
            let load = response.port_load.as_ref().unwrap();
            for i in 0..2 {
                let mut force = C64::ZERO;
                for j in 0..2 {
                    force = force
                        + response.displacement[j]
                            * C64::new(
                                STIFFNESS[i * 2 + j] - w * w * MASS[i * 2 + j],
                                -w * DAMPING[i * 2 + j],
                            );
                    force = force + load.effort[j].scale(MAP[i * 2 + j]);
                }
                close(force, FORCE[i], 2e-13);
            }
            assert!(response.average_viscous_power_w >= 0.0);
            assert!(load.average_absorbed_power_w > 0.0);
            assert!(
                response.power_balance_residual_w.abs() < 1e-12 * response.average_input_power_w
            );
            let repeated = prepared.solve(&FORCE, cx).unwrap();
            assert_eq!(response.displacement, repeated.displacement);
        }
    });
}

#[test]
fn g3_coordinate_and_power_conjugate_port_unit_changes_preserve_response() {
    with_cx(&CancelGate::new(), |cx| {
        let w = 1.7;
        let original = PreparedHarmonicResponse::try_with_port_impedance(
            &MASS,
            &DAMPING,
            &STIFFNESS,
            2,
            &MAP,
            2,
            &LOAD,
            w,
            controls(),
            ports(),
            cx,
        )
        .unwrap()
        .solve(&FORCE, cx)
        .unwrap();
        let coordinate = [1e-40, 1e40];
        let port_unit = [1e30, 1e-30];
        let congruence = |source: &[f64; 4]| {
            let mut out = [0.0; 4];
            for i in 0..2 {
                for j in 0..=i {
                    out[i * 2 + j] = (source[i * 2 + j] * coordinate[i]) * coordinate[j];
                    out[j * 2 + i] = out[i * 2 + j];
                }
            }
            out
        };
        let mut b = [0.0; 4];
        let mut z = [C64::ZERO; 4];
        for i in 0..2 {
            for j in 0..2 {
                b[i * 2 + j] = (MAP[i * 2 + j] * coordinate[i]) * port_unit[j];
                z[i * 2 + j] = LOAD[i * 2 + j]
                    .scale(1.0 / port_unit[i])
                    .scale(1.0 / port_unit[j]);
            }
        }
        let transformed = PreparedHarmonicResponse::try_with_port_impedance(
            &congruence(&MASS),
            &congruence(&DAMPING),
            &congruence(&STIFFNESS),
            2,
            &b,
            2,
            &z,
            w,
            controls(),
            ports(),
            cx,
        )
        .unwrap()
        .solve(
            &[FORCE[0].scale(coordinate[0]), FORCE[1].scale(coordinate[1])],
            cx,
        )
        .unwrap();
        for i in 0..2 {
            close(
                transformed.displacement[i].scale(coordinate[i]),
                original.displacement[i],
                5e-13,
            );
            close(
                transformed.velocity[i].scale(coordinate[i]),
                original.velocity[i],
                5e-13,
            );
            let (t, o) = (
                transformed.port_load.as_ref().unwrap(),
                original.port_load.as_ref().unwrap(),
            );
            close(t.flow[i].scale(1.0 / port_unit[i]), o.flow[i], 5e-13);
            close(t.effort[i].scale(port_unit[i]), o.effort[i], 5e-13);
        }
        scalar(
            transformed.average_input_power_w,
            original.average_input_power_w,
            5e-13,
        );
        scalar(
            transformed.port_load.unwrap().average_absorbed_power_w,
            original.port_load.unwrap().average_absorbed_power_w,
            5e-13,
        );
    });
}

#[test]
fn g1_lossless_nonreciprocal_and_rank_deficient_load_channels_are_retained() {
    with_cx(&CancelGate::new(), |cx| {
        // A lossless gyrator has skew real Z and a zero Hermitian part. It
        // couples the two modes and regularizes the first bare pole without heat.
        let z = [C64::ZERO, C64::from_re(0.5), C64::from_re(-0.5), C64::ZERO];
        let r = PreparedHarmonicResponse::try_modal_with_port_impedance(
            &[2.0, 3.0],
            &[0.0; 4],
            &[1.0, 0.0, 0.0, 1.0],
            2,
            &z,
            2.0,
            controls(),
            ports(),
            cx,
        )
        .unwrap()
        .solve(&[C64::ONE, C64::ZERO], cx)
        .unwrap();
        let port = r.port_load.unwrap();
        assert!(port.average_absorbed_power_w.abs() < 1e-14);
        assert!(r.power_balance_residual_w.abs() < 1e-14);

        // A rank-one reactive impedance with a zero force-map column needs no
        // inverse load admittance and retains the exact unused port flow.
        let r = PreparedHarmonicResponse::try_modal_with_port_impedance(
            &[2.0],
            &[0.0],
            &[1.0, 0.0],
            2,
            &[C64::new(0.0, 1.0); 4],
            2.0,
            controls(),
            ports(),
            cx,
        )
        .unwrap()
        .solve(&[C64::ONE], cx)
        .unwrap();
        close(r.displacement[0], C64::from_re(0.5), 1e-15);
        assert_eq!(r.port_load.unwrap().flow[1], C64::ZERO);
    });
}

#[test]
fn g1_tiny_physical_load_resolves_high_q_without_a_damping_floor() {
    with_cx(&CancelGate::new(), |cx| {
        let w = 0.1;
        for loss in [1e-6, 1e-14, 1e-80] {
            let r = PreparedHarmonicResponse::try_modal_with_port_impedance(
                &[w],
                &[0.0],
                &[1.0],
                1,
                &[C64::from_re(loss)],
                w,
                controls(),
                ports(),
                cx,
            )
            .unwrap()
            .solve(&[C64::ONE], cx)
            .unwrap();
            close(r.velocity[0], C64::from_re(1.0 / loss), 2e-15);
            assert!(r.power_balance_residual_w.abs() < 1e-14 * r.average_input_power_w);
        }
        assert!(matches!(
            PreparedHarmonicResponse::try_modal_with_port_impedance(
                &[w],
                &[0.0],
                &[1.0],
                1,
                &[C64::ZERO],
                w,
                controls(),
                ports(),
                cx
            ),
            Err(HarmonicResponseError::Singular { .. })
        ));
    });
}

#[test]
fn g0_full_hermitian_passivity_shapes_budgets_and_cancellation_are_enforced() {
    let gate = CancelGate::new();
    with_cx(&gate, |cx| {
        let make = |z: &[C64], b: &[f64], p: usize, limit| {
            PreparedHarmonicResponse::try_with_port_impedance(
                &MASS,
                &DAMPING,
                &STIFFNESS,
                2,
                b,
                p,
                z,
                1.0,
                controls(),
                limit,
                cx,
            )
        };
        // Real(Z) alone is positive definite, but the complex Hermitian part
        // has eigenvalue -1. A real-only passivity check would miss this source.
        let active = [C64::ONE, C64::new(0.0, 2.0), C64::new(0.0, -2.0), C64::ONE];
        assert!(matches!(
            make(&active, &MAP, 2, ports()),
            Err(HarmonicResponseError::NotPositiveSemidefinite { .. })
        ));
        let semidefinite = [C64::ONE, C64::new(0.0, 1.0), C64::new(0.0, -1.0), C64::ONE];
        assert!(make(&semidefinite, &MAP, 2, ports()).is_ok());
        assert!(make(&LOAD[..3], &MAP, 2, ports()).is_err());
        assert!(make(&LOAD, &MAP[..3], 2, ports()).is_err());
        let mut bad = LOAD;
        bad[1].im = f64::NAN;
        assert!(make(&bad, &MAP, 2, ports()).is_err());
        assert!(make(&[], &[], 0, ports()).is_err());
        assert!(
            make(
                &LOAD,
                &MAP,
                2,
                HarmonicPortControls {
                    maximum_ports: 1,
                    ..ports()
                }
            )
            .is_err()
        );
        assert!(
            make(
                &LOAD,
                &MAP,
                2,
                HarmonicPortControls {
                    maximum_projection_terms: 15,
                    ..ports()
                }
            )
            .is_err()
        );
        assert!(
            make(
                &LOAD,
                &MAP,
                2,
                HarmonicPortControls {
                    maximum_payload_bytes: 1,
                    ..ports()
                }
            )
            .is_err()
        );
        let prepared = make(&LOAD, &MAP, 2, ports()).unwrap();
        let before = prepared.solve(&FORCE, cx).unwrap();
        assert!(prepared.solve(&[C64::ONE], cx).is_err());
        assert_eq!(
            before.displacement,
            prepared.solve(&FORCE, cx).unwrap().displacement
        );
        gate.request();
        assert!(matches!(
            prepared.solve(&FORCE, cx),
            Err(HarmonicResponseError::Cancelled)
        ));
        assert!(matches!(
            make(&LOAD, &MAP, 2, ports()),
            Err(HarmonicResponseError::Cancelled)
        ));
    });
}
