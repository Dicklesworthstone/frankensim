use super::WakeVortex;
use super::convection::advance;

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "actual={actual:e}, expected={expected:e}, tolerance={tolerance:e}"
    );
}

fn moments(wake: &[WakeVortex]) -> [f64; 3] {
    let mut out = [0.0; 3];
    for w in wake {
        out[0] += w.gamma * w.pos[0];
        out[1] += w.gamma * w.pos[1];
        out[2] += w.gamma * (w.pos[0] * w.pos[0] + w.pos[1] * w.pos[1]);
    }
    out
}

#[test]
fn co_rotating_wake_has_clockwise_phase_without_spiral_drift() {
    let mut wake = [
        WakeVortex {
            pos: [-0.5, 0.0],
            gamma: 1.0,
        },
        WakeVortex {
            pos: [0.5, 0.0],
            gamma: 1.0,
        },
    ];
    let dt = 0.1;
    let core = 0.2;
    for _ in 0..1000 {
        advance(&mut wake, [0.0; 2], 0.0, [0.25, 0.0], dt, core).unwrap();
    }
    let angle = -2.0 * 1000.0 * dt / (std::f64::consts::TAU * (1.0 + core * core));
    close(wake[1].pos[0], 0.5 * angle.cos(), 2e-12);
    close(wake[1].pos[1], 0.5 * angle.sin(), 2e-12);
    close(moments(&wake)[2], 0.5, 3e-13);
}

#[test]
fn translating_dipole_retains_the_bem_sign_and_freestream_boost() {
    let mut wake = [
        WakeVortex {
            pos: [-0.5, 0.0],
            gamma: 1.0,
        },
        WakeVortex {
            pos: [0.5, 0.0],
            gamma: -1.0,
        },
    ];
    let initial = wake;
    let flow = [0.6, -0.4];
    let core = 0.1;
    let dt = 0.04;
    let duration = 100.0 * dt;
    for _ in 0..100 {
        advance(&mut wake, flow, 0.0, [0.25, 0.0], dt, core).unwrap();
    }
    let induced_y = -1.0 / (std::f64::consts::TAU * (1.0 + core * core));
    for (w, start) in wake.iter().zip(initial) {
        close(w.pos[0], start.pos[0] + duration * flow[0], 2e-13);
        close(
            w.pos[1],
            start.pos[1] + duration * (flow[1] + induced_y),
            2e-13,
        );
        assert_eq!(w.gamma.to_bits(), start.gamma.to_bits());
    }
}

#[test]
fn mutual_induction_preserves_linear_and_angular_impulse_moments() {
    let mut wake = [
        WakeVortex {
            pos: [-0.9, 0.1],
            gamma: 0.8,
        },
        WakeVortex {
            pos: [0.5, 0.6],
            gamma: -1.3,
        },
        WakeVortex {
            pos: [0.9, -0.2],
            gamma: 0.5,
        },
        WakeVortex {
            pos: [0.2, -0.8],
            gamma: 0.0,
        },
    ];
    let initial = moments(&wake);
    for _ in 0..500 {
        advance(&mut wake, [0.0; 2], 0.0, [0.25, 0.0], 0.03, 0.2).unwrap();
    }
    for (actual, expected) in moments(&wake).into_iter().zip(initial) {
        close(actual, expected, 2e-12);
    }
}

#[test]
fn fixed_bound_source_follows_the_exact_regularized_circular_streamline() {
    let quarter = [0.25, -0.1];
    let mut wake = [WakeVortex {
        pos: [1.25, -0.1],
        gamma: -2.0,
    }];
    let bound = 0.7;
    let core = 0.1;
    for _ in 0..300 {
        advance(&mut wake, [0.0; 2], bound, quarter, 0.05, core).unwrap();
    }
    let angle = -bound * 15.0 / (std::f64::consts::TAU * (1.0 + core * core));
    close(wake[0].pos[0], quarter[0] + angle.cos(), 2e-13);
    close(wake[0].pos[1], quarter[1] + angle.sin(), 2e-13);
    assert_eq!(wake[0].gamma, -2.0);
}

fn velocity(wake: &[WakeVortex; 3], flow: [f64; 2], bound: f64, core: f64) -> [[f64; 2]; 3] {
    std::array::from_fn(|i| {
        let mut v = flow;
        for source in wake.iter().chain([&WakeVortex {
            pos: [0.25, 0.0],
            gamma: bound,
        }]) {
            let dx = wake[i].pos[0] - source.pos[0];
            let dy = wake[i].pos[1] - source.pos[1];
            let factor = source.gamma / (std::f64::consts::TAU * (dx * dx + dy * dy + core * core));
            v[0] += factor * dy;
            v[1] -= factor * dx;
        }
        v
    })
}

fn shifted(wake: [WakeVortex; 3], k: [[f64; 2]; 3], h: f64) -> [WakeVortex; 3] {
    std::array::from_fn(|i| WakeVortex {
        pos: std::array::from_fn(|axis| wake[i].pos[axis] + h * k[i][axis]),
        gamma: wake[i].gamma,
    })
}

#[test]
fn frozen_background_convection_is_reversible_and_second_order() {
    let initial = [
        WakeVortex {
            pos: [0.9, 0.2],
            gamma: -0.4,
        },
        WakeVortex {
            pos: [1.2, -0.3],
            gamma: -0.2,
        },
        WakeVortex {
            pos: [1.7, 0.1],
            gamma: 0.1,
        },
    ];
    let flow = [0.8, 0.1];
    let bound = 0.5;
    let core = 0.15;
    let mut reversed = initial;
    advance(&mut reversed, flow, bound, [0.25, 0.0], 0.2, core).unwrap();
    advance(&mut reversed, flow, bound, [0.25, 0.0], -0.2, core).unwrap();
    for (actual, expected) in reversed.iter().zip(initial) {
        for axis in 0..2 {
            close(actual.pos[axis], expected.pos[axis], 2e-14);
        }
    }

    // Independent small-step RK4 of the original clockwise Biot–Savart ODE.
    let mut reference = initial;
    let h = 0.5 / 2000.0;
    for _ in 0..2000 {
        let k1 = velocity(&reference, flow, bound, core);
        let k2 = velocity(&shifted(reference, k1, 0.5 * h), flow, bound, core);
        let k3 = velocity(&shifted(reference, k2, 0.5 * h), flow, bound, core);
        let k4 = velocity(&shifted(reference, k3, h), flow, bound, core);
        for i in 0..3 {
            for axis in 0..2 {
                reference[i].pos[axis] +=
                    h / 6.0 * (k1[i][axis] + 2.0 * k2[i][axis] + 2.0 * k3[i][axis] + k4[i][axis]);
            }
        }
    }
    let errors: [f64; 3] = [10, 20, 40].map(|steps| {
        let mut wake = initial;
        for _ in 0..steps {
            advance(
                &mut wake,
                flow,
                bound,
                [0.25, 0.0],
                0.5 / f64::from(steps),
                core,
            )
            .unwrap();
        }
        wake.iter()
            .zip(reference)
            .flat_map(|(a, b)| [a.pos[0] - b.pos[0], a.pos[1] - b.pos[1]])
            .map(|e| e * e)
            .sum::<f64>()
            .sqrt()
    });
    assert!(errors[2] < 2e-6, "errors={errors:?}");
    for pair in errors.windows(2) {
        assert!(
            (3.8..4.2).contains(&(pair[0] / pair[1])),
            "errors={errors:?}"
        );
    }
}
