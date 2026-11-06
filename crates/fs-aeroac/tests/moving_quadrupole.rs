//! Independent Green-Hessian and static-limit checks of moving physical stress.

use fs_aeroac::curle3d::moving::{MovingCompactLimits, MovingCompactObserver, UniformSourceMotion};
use fs_aeroac::curle3d::quadrupole::transient::{
    PreparedTransientQuadrupoleObserver, TransientQuadrupoleSource,
};
use fs_material::fluid_acoustics::FluidAcousticState;

const TENSOR: [f64; 6] = [0.7, -0.2, 0.4, 0.3, -0.1, 0.23];
const ENTRIES: [(usize, usize); 6] = [(0, 0), (1, 1), (2, 2), (0, 1), (0, 2), (1, 2)];

fn observer(x: [f64; 3], c: f64) -> MovingCompactObserver {
    MovingCompactObserver::new(
        FluidAcousticState::try_new(1.2, c).unwrap(),
        x,
        MovingCompactLimits {
            maximum_mach: 0.8,
            minimum_distance_m: 1e-6,
        },
    )
    .unwrap()
}

fn strength(t: f64, omega: f64, tensor: [f64; 6]) -> TransientQuadrupoleSource {
    let (s, c) = (omega * t).sin_cos();
    TransientQuadrupoleSource {
        stress_moment_n_m: tensor.map(|a| a * c),
        stress_rate_n_m_s: tensor.map(|a| -omega * a * s),
        stress_acceleration_n_m_s2: tensor.map(|a| -omega * omega * a * c),
        ..Default::default()
    }
}

fn close(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() < tolerance,
        "actual={a:e}, expected={b:e}, tolerance={tolerance:e}"
    );
}

#[test]
fn stationary_limit_retains_trace_off_diagonal_and_all_radial_orders() {
    for c in [317.0, 1490.0] {
        for range in [0.003, 0.2, 12.0] {
            for n in [[1.0, 0.0, 0.0], [0.6, 0.0, 0.8], [-0.48, 0.64, 0.6]] {
                let x = n.map(|a| range * a);
                let motion = UniformSourceMotion::new(0.0, [0.0; 3], [0.0; 3]).unwrap();
                let moving = observer(x, c)
                    .prepare_uniform_quadrupole(0.137, &motion, 1e-12)
                    .unwrap();
                let fixed = PreparedTransientQuadrupoleObserver::new(
                    FluidAcousticState::try_new(1.2, c).unwrap(),
                    x,
                    &[[0.0; 3]],
                    1,
                )
                .unwrap();
                // Isolate each derivative order, so near-field terms cannot
                // hide underneath a dominant radiation contribution.
                for order in 0..3 {
                    let mut q = TransientQuadrupoleSource::default();
                    match order {
                        0 => q.stress_moment_n_m = TENSOR,
                        1 => q.stress_rate_n_m_s = TENSOR,
                        _ => q.stress_acceleration_n_m_s2 = TENSOR,
                    }
                    let expected = fixed.evaluate_at_retarded(&[q]).unwrap().pressure_pa;
                    close(
                        moving.pressure_pa(q).unwrap(),
                        expected,
                        5e-13 * expected.abs().max(1e-20),
                    );
                }
            }
        }
    }
}

// Independent bisection of the light cone, followed by the scalar retarded
// Green potential. No analytic retarded derivative or production root formula.
fn potential(x: [f64; 3], t: f64, y0: [f64; 3], v: [f64; 3], omega: f64) -> f64 {
    let c = 343.0;
    let mut lower = t - 2.0;
    let mut upper = t;
    for _ in 0..100 {
        let tau = 0.5 * (lower + upper);
        let r: [f64; 3] = core::array::from_fn(|i| x[i] - y0[i] - v[i] * tau);
        let range = (r.iter().map(|r| r * r).sum::<f64>()).sqrt();
        if tau + range / c > t {
            upper = tau;
        } else {
            lower = tau;
        }
    }
    let tau = 0.5 * (lower + upper);
    let r: [f64; 3] = core::array::from_fn(|i| x[i] - y0[i] - v[i] * tau);
    let range = (r.iter().map(|r| r * r).sum::<f64>()).sqrt();
    let radial_velocity = (0..3).map(|i| r[i] * v[i]).sum::<f64>() / range;
    (omega * tau).cos() / (4.0 * std::f64::consts::PI * range * (1.0 - radial_velocity / c))
}

fn hessian_pressure(x: [f64; 3], t: f64, y0: [f64; 3], v: [f64; 3], omega: f64, h: f64) -> f64 {
    let mut pressure = 0.0;
    for (entry, (i, j)) in ENTRIES.into_iter().enumerate() {
        let derivative = if i == j {
            let mut plus = x;
            plus[i] += h;
            let mut minus = x;
            minus[i] -= h;
            (potential(plus, t, y0, v, omega) - 2.0 * potential(x, t, y0, v, omega)
                + potential(minus, t, y0, v, omega))
                / (h * h)
        } else {
            let mut pp = x;
            pp[i] += h;
            pp[j] += h;
            let mut pm = x;
            pm[i] += h;
            pm[j] -= h;
            let mut mp = x;
            mp[i] -= h;
            mp[j] += h;
            let mut mm = x;
            mm[i] -= h;
            mm[j] -= h;
            2.0 * (potential(pp, t, y0, v, omega)
                - potential(pm, t, y0, v, omega)
                - potential(mp, t, y0, v, omega)
                + potential(mm, t, y0, v, omega))
                / (4.0 * h * h)
        };
        pressure += TENSOR[entry] * derivative;
    }
    pressure
}

#[test]
fn moving_pressure_converges_to_independent_spatial_green_hessian() {
    let x = [4.0, -2.0, 3.0];
    let y0 = [-0.2, 0.1, 0.3];
    let t = 0.03;
    for v in [[80.0, 0.0, 0.0], [-80.0, 20.0, 0.0], [0.0, 100.0, -30.0]] {
        let motion = UniformSourceMotion::new(0.0, y0, v).unwrap();
        let prepared = observer(x, 343.0)
            .prepare_uniform_quadrupole(t, &motion, 1e-12)
            .unwrap();
        for omega in [0.0, 2.0 * std::f64::consts::PI * 40.0] {
            let exact = prepared
                .pressure_pa(strength(prepared.source_time_s(), omega, TENSOR))
                .unwrap();
            let errors = [0.02, 0.01, 0.005]
                .map(|h| (hessian_pressure(x, t, y0, v, omega, h) - exact).abs());
            assert!(
                errors[1] < 0.35 * errors[0] && errors[2] < 0.35 * errors[1],
                "v={v:?}, omega={omega}, errors={errors:?}"
            );
            assert!(
                errors[2] < 1e-5 * exact.abs().max(1e-4),
                "exact={exact:e}, errors={errors:?}"
            );
        }
    }
}

#[test]
fn radial_doppler_and_moving_microphone_use_two_medium_frame_velocities() {
    let c = 343.0;
    let source_speed = 0.2 * c;
    let motion = UniformSourceMotion::new(0.0, [0.0; 3], [source_speed, 0.0, 0.0]).unwrap();
    let p = observer([10.0, 0.0, 0.0], c)
        .prepare_uniform_quadrupole(0.02, &motion, 1e-12)
        .unwrap();
    close(p.propagation().source_time_rate(), 1.25, 1e-14);
    close(
        p.observed_source_time_rate([source_speed, 0.0, 0.0])
            .unwrap(),
        1.0,
        1e-14,
    );
    close(
        p.observed_source_time_rate([-0.1 * c, 0.0, 0.0]).unwrap(),
        1.1 / 0.8,
        1e-14,
    );

    // Numerical derivative along an independently moving microphone worldline.
    let vr = [-0.1 * c, 20.0, 5.0];
    let t = 0.04;
    let h = 1e-6;
    let x0 = [12.0, 3.0, 1.0];
    let prepare = |time: f64| {
        observer(core::array::from_fn(|i| x0[i] + vr[i] * time), c)
            .prepare_uniform_quadrupole(time, &motion, 1e-12)
            .unwrap()
    };
    let rate = (prepare(t + h).source_time_s() - prepare(t - h).source_time_s()) / (2.0 * h);
    close(
        prepare(t).observed_source_time_rate(vr).unwrap(),
        rate,
        2e-10,
    );
    assert!(p.observed_source_time_rate([c, 0.0, 0.0]).is_err());
    assert!(p.observed_source_time_rate([f64::NAN, 0.0, 0.0]).is_err());
}

#[test]
fn tensor_and_geometry_rotate_together_without_altering_pressure() {
    // A cyclic orthogonal rotation takes (x,y,z) -> (z,x,y), including the
    // off-diagonal tensor entries. Nonzero trace is intentionally retained.
    let rotate = |x: [f64; 3]| [x[2], x[0], x[1]];
    let x = [8.0, 3.0, -1.0];
    let v = [40.0, -20.0, 30.0];
    let y = [0.4, 0.1, -0.2];
    let tensor = [
        TENSOR[2], TENSOR[0], TENSOR[1], TENSOR[4], TENSOR[5], TENSOR[3],
    ];
    let a = observer(x, 343.0)
        .prepare_uniform_quadrupole(0.1, &UniformSourceMotion::new(0.0, y, v).unwrap(), 1e-12)
        .unwrap();
    let b = observer(rotate(x), 343.0)
        .prepare_uniform_quadrupole(
            0.1,
            &UniformSourceMotion::new(0.0, rotate(y), rotate(v)).unwrap(),
            1e-12,
        )
        .unwrap();
    close(
        a.pressure_pa(strength(a.source_time_s(), 250.0, TENSOR))
            .unwrap(),
        b.pressure_pa(strength(b.source_time_s(), 250.0, tensor))
            .unwrap(),
        1e-14,
    );
}

#[test]
fn pressure_rejects_invalid_kinematics_and_stress_without_inventing_an_integral() {
    let obs = observer([1.0, 0.0, 0.0], 343.0);
    assert!(
        obs.prepare_uniform_quadrupole(
            0.0,
            &UniformSourceMotion::new(0.0, [0.0; 3], [300.0, 0.0, 0.0]).unwrap(),
            1e-12
        )
        .is_err()
    );
    assert!(
        obs.prepare_uniform_quadrupole(
            0.0,
            &UniformSourceMotion::new(0.0, [1.0, 0.0, 0.0], [0.0; 3]).unwrap(),
            1e-12
        )
        .is_err()
    );
    let prepared = obs
        .prepare_uniform_quadrupole(
            0.1,
            &UniformSourceMotion::new(0.0, [0.0; 3], [0.0; 3]).unwrap(),
            1e-12,
        )
        .unwrap();
    let mut q = strength(0.0, 300.0, TENSOR);
    q.stress_integral_n_m_s = [f64::NAN; 6];
    assert!(prepared.pressure_pa(q).is_ok());
    q.stress_rate_n_m_s[2] = f64::INFINITY;
    assert!(prepared.pressure_pa(q).is_err());
}
