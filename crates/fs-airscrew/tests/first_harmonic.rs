//! G0/G1/G3: physical modal inertia, actual section work and objective frames.
use fs_airscrew::dynamic_inflow::UniformAxialInflow;
use fs_airscrew::dynamic_inflow::first_harmonic::{FirstHarmonicInflow, HarmonicInduction};
use fs_airscrew::oblique::{AzimuthQuadrature, DiskFrame, SpinDirection};
use fs_airscrew::{BladeStation, Rotor};

const RHO: f64 = 1.225;
const R: f64 = 1.3;
const AXIS: [f64; 3] = [0.0, 0.0, 1.0];
const RADIAL: [f64; 3] = [1.0, 0.0, 0.0];
const ZERO: [f64; 3] = [0.0; 3];

fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a.into_iter().zip(b).map(|(a, b)| a * b).sum()
}
fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
fn near(a: f64, b: f64, tol: f64) {
    assert!((a - b).abs() <= tol, "{a:e} != {b:e}, error={:e}", a - b);
}
fn near_vec(a: [f64; 3], b: [f64; 3], tol: f64) {
    for k in 0..3 {
        near(a[k], b[k], tol);
    }
}
fn field(mean: f64, x: f64, y: f64) -> HarmonicInduction {
    HarmonicInduction {
        mean_mps: mean,
        first_mps: [x, y],
    }
}
fn frame() -> DiskFrame {
    DiskFrame::new(AXIS, RADIAL, SpinDirection::AlongThrust).unwrap()
}
fn rotor() -> Rotor {
    Rotor {
        radius_m: R,
        n_blades: 3,
        camber_ratio: 0.02,
        stations: [0.25, 0.45, 0.65, 0.85]
            .map(|r| BladeStation {
                r_over_r: r,
                chord_m: 0.12 - 0.04 * r,
                beta_rad: 9.0_f64.atan2(60.0 * R * r) + 0.06,
            })
            .to_vec(),
    }
}

#[test]
fn normal_harmonic_gain_and_air_mass_match_independent_exponential_response() {
    let model = FirstHarmonicInflow::new(R, RHO).unwrap();
    let mass = 16.0 / 45.0 * RHO * R.powi(3);
    let damping = 0.5 * RHO * std::f64::consts::PI * R * R * (5.0 + 4.0);
    near(model.modal_mass_kg()[0], 8.0 / 3.0 * RHO * R.powi(3), 1e-15);
    near(model.modal_mass_kg()[1], mass, 1e-15);
    let mean_force = 2.0 * RHO * std::f64::consts::PI * R * R * 2.0 * 7.0;
    let forces = [mean_force, 0.2 * damping, -0.1 * damping];
    let start = field(2.0, -0.3, 0.2);
    let horizon = 0.08;
    let expected = [
        0.2 + (-0.3 - 0.2) * (-damping * horizon / mass).exp(),
        -0.1 + (0.2 + 0.1) * (-damping * horizon / mass).exp(),
    ];
    let run = |n: usize| {
        let mut state = start;
        let mut impulse = [0.0; 3];
        let mut work = 0.0;
        let mut loss = 0.0;
        for _ in 0..n {
            let step = model
                .step_held(state, 5.0, forces, 0.0, horizon / n as f64)
                .unwrap();
            for k in 0..3 {
                impulse[k] += step.driving_impulse_ns[k] - step.transport_impulse_ns[k];
                assert!(step.momentum_residual_ns[k].abs() < 5e-14);
            }
            work += step.induced_work_j;
            loss += step.transport_work_j;
            assert!(step.transport_work_j >= 0.0);
            assert!(step.energy_residual_j.abs() < 3e-14);
            state = step.induction;
        }
        let before = start.modes();
        let after = state.modes();
        let masses = model.modal_mass_kg();
        for k in 0..3 {
            near(masses[k] * (after[k] - before[k]), impulse[k], 2e-12);
        }
        let delta_energy: f64 = (0..3)
            .map(|k| 0.5 * masses[k] * (after[k] * after[k] - before[k] * before[k]))
            .sum();
        near(delta_energy, work - loss, 1e-12);
        (state.first_mps[0] - expected[0]).hypot(state.first_mps[1] - expected[1])
    };
    let errors = [run(16), run(32), run(64)];
    println!("normal first-harmonic analytic errors: {errors:?}");
    for pair in errors.windows(2) {
        assert!((3.8..4.2).contains(&(pair[0] / pair[1])));
    }
}

#[test]
fn rotating_disk_reference_transports_asymmetry_without_artificial_work() {
    let model = FirstHarmonicInflow::new(R, RHO).unwrap();
    let mass = 16.0 / 45.0 * RHO * R.powi(3);
    let damping = 0.5 * RHO * std::f64::consts::PI * R * R * (3.0 + 2.0);
    let thrust = 2.0 * RHO * std::f64::consts::PI * R * R * 1.0 * 4.0;
    let omega = 12.0;
    let horizon = 0.06;
    let mut state = field(1.0, 0.3, -0.2);
    for _ in 0..128 {
        let step = model
            .step_held(state, 3.0, [thrust, 0.0, 0.0], omega, horizon / 128.0)
            .unwrap();
        near(step.frame_work_j, 0.0, 5e-18);
        assert!(step.energy_residual_j.abs() < 2e-14);
        state = step.induction;
    }
    let (s, c) = (-omega * horizon).sin_cos();
    let decay = (-damping * horizon / mass).exp();
    near(state.first_mps[0], decay * (0.3 * c + 0.2 * s), 2e-6);
    near(state.first_mps[1], decay * (0.3 * s - 0.2 * c), 2e-6);
}

#[test]
fn zero_modes_retain_the_original_uniform_wake_and_oblique_section_arithmetic() {
    let uniform = UniformAxialInflow::new(R, RHO).unwrap();
    let modes = FirstHarmonicInflow::new(R, RHO).unwrap();
    let old = uniform.step_held(2.0, 7.0, 150.0, 0.007).unwrap();
    let new = modes
        .step_held(field(2.0, 0.0, 0.0), 7.0, [150.0, 0.0, 0.0], 4.0, 0.007)
        .unwrap();
    assert_eq!(old.induced_mps, new.induction.mean_mps);
    assert_eq!(new.induction.first_mps, [0.0, 0.0]);
    let rule = AzimuthQuadrature::new(32).unwrap();
    let old = rule
        .evaluate(
            &rotor(),
            RHO,
            frame(),
            [1.0, 2.0, 7.0],
            [0.1, -0.2, 3.0],
            2.0,
            60.0,
        )
        .unwrap();
    let new = rule
        .evaluate_first_harmonic(
            &rotor(),
            RHO,
            frame(),
            [1.0, 2.0, 7.0],
            [0.1, -0.2, 3.0],
            field(2.0, 0.0, 0.0),
            60.0,
        )
        .unwrap();
    assert_eq!(old, new);
}

// Independently integrate attached lift/drag with 8192 midpoint azimuths.
// This does not call the production section owner or its deterministic trig.
fn reference(induction: HarmonicInduction) -> ([f64; 3], [f64; 3], f64, f64) {
    let rotor = rotor();
    let mut force = ZERO;
    let mut moment = ZERO;
    let mut torque = 0.0;
    let mut work = 0.0;
    for (station, dr) in rotor
        .stations
        .iter()
        .zip([0.15, 0.2, 0.2, 0.25].map(|x| x * R))
    {
        let r = R * station.r_over_r;
        for i in 0..8192 {
            let angle = std::f64::consts::TAU * (i as f64 + 0.5) / 8192.0;
            let (s, c) = angle.sin_cos();
            let radial = [c, s, 0.0];
            let tangent = [-s, c, 0.0];
            let arm = radial.map(|x| x * r);
            let tilt = cross([0.2, -0.1, 0.0], arm);
            let incident = [4.0 + tilt[0], -1.0 + tilt[1], 7.0 + tilt[2]];
            let vi = induction.mean_mps
                + station.r_over_r * (induction.first_mps[0] * c + induction.first_mps[1] * s);
            let ua = incident[2] + vi;
            let ut = 60.0 * r + dot(incident, tangent);
            let speed = ua.hypot(ut);
            let alpha = station.beta_rad - ua.atan2(ut);
            assert!(alpha.abs() < 0.25);
            let cl = std::f64::consts::TAU * (alpha + 0.04);
            let cd = 0.012 + 0.8 * alpha * alpha;
            let q = 0.5 * RHO * speed * speed * station.chord_m * 3.0 * dr / 8192.0;
            let thrust = q * (cl * ut - cd * ua) / speed;
            let drag = q * (cl * ua + cd * ut) / speed;
            let f = [-drag * tangent[0], -drag * tangent[1], thrust];
            let m = cross(arm, f);
            for k in 0..3 {
                force[k] += f[k];
                moment[k] += m[k];
            }
            torque += r * drag;
            work += thrust * vi;
        }
    }
    (force, moment, torque, work)
}

#[test]
fn spatial_induction_enters_real_sections_and_conjugate_moment_work() {
    let u = field(2.0, 0.8, -0.5);
    let expected = reference(u);
    let load = AzimuthQuadrature::new(64)
        .unwrap()
        .evaluate_first_harmonic(
            &rotor(),
            RHO,
            frame(),
            [4.0, -1.0, 7.0],
            [0.2, -0.1, 5.0],
            u,
            60.0,
        )
        .unwrap();
    near_vec(load.force_n, expected.0, 3e-10);
    near_vec(
        load.hub_moment_nm,
        [expected.1[0], expected.1[1], 0.0],
        3e-10,
    );
    near(load.torque_nm, expected.2, 3e-10);
    near(load.induced_work_w, expected.3, 3e-10);
    let q = frame().wake_force_n(&load, R).unwrap();
    near(q[1], -expected.1[1] / R, 3e-10);
    near(q[2], expected.1[0] / R, 3e-10);
    near(dot(q, u.modes()), load.induced_work_w, 3e-11);
    assert!((load.induced_work_w - load.thrust_n * u.mean_mps).abs() > 0.1);
    near(load.power_residual_w, 0.0, 2e-9);
    near(
        load.kinematic_power_w + load.induced_work_w + load.profile_dissipation_w,
        0.0,
        2e-9,
    );
}

fn rotate(a: [f64; 3]) -> [f64; 3] {
    let (s, c) = 0.54_f64.sin_cos();
    let b = [c * a[0] - s * a[1], s * a[0] + c * a[1], a[2]];
    let (s, c) = 0.31_f64.sin_cos();
    [b[0], c * b[1] - s * b[2], s * b[1] + c * b[2]]
}

#[test]
fn common_rotation_and_disk_reference_rebasing_preserve_physical_loads() {
    let rule = AzimuthQuadrature::new(64).unwrap();
    let u = field(2.0, 0.8, -0.5);
    let velocity = [4.0, -1.0, 7.0];
    let omega = [0.2, -0.1, 5.0];
    let a = rule
        .evaluate_first_harmonic(&rotor(), RHO, frame(), velocity, omega, u, 60.0)
        .unwrap();
    let rotated = DiskFrame::new(rotate(AXIS), rotate(RADIAL), SpinDirection::AlongThrust).unwrap();
    let b = rule
        .evaluate_first_harmonic(
            &rotor(),
            RHO,
            rotated,
            rotate(velocity),
            rotate(omega),
            u,
            60.0,
        )
        .unwrap();
    near_vec(b.force_n, rotate(a.force_n), 4e-11);
    near_vec(b.hub_moment_nm, rotate(a.hub_moment_nm), 4e-11);
    near(a.induced_work_w, b.induced_work_w, 4e-11);
    let (s, c) = 0.31_f64.sin_cos();
    let other = DiskFrame::new(AXIS, [c, s, 0.0], SpinDirection::AlongThrust).unwrap();
    let other_u = frame().rebase_induction(other, u).unwrap();
    let b = rule
        .evaluate_first_harmonic(&rotor(), RHO, other, velocity, omega, other_u, 60.0)
        .unwrap();
    near_vec(b.force_n, a.force_n, 4e-11);
    near_vec(b.hub_moment_nm, a.hub_moment_nm, 4e-11);
    near(
        dot(frame().wake_force_n(&a, R).unwrap(), u.modes()),
        dot(other.wake_force_n(&b, R).unwrap(), other_u.modes()),
        4e-11,
    );
    let wf = FirstHarmonicInflow::new(R, RHO).unwrap();
    let sa = wf
        .step_held(u, 7.0, frame().wake_force_n(&a, R).unwrap(), 5.0, 0.004)
        .unwrap();
    let sb = wf
        .step_held(other_u, 7.0, other.wake_force_n(&b, R).unwrap(), 5.0, 0.004)
        .unwrap();
    near_vec(
        frame()
            .rebase_induction(other, sa.induction)
            .unwrap()
            .modes(),
        sb.induction.modes(),
        2e-12,
    );
    near(sa.kinetic_after_j, sb.kinetic_after_j, 2e-11);
}

#[test]
fn whole_disk_flow_and_late_endpoint_are_admitted_without_clipping() {
    let model = FirstHarmonicInflow::new(R, RHO).unwrap();
    let initial = field(1.0, 0.1, 0.2);
    assert!(model.rate(field(1.0, 3.1, 0.0), 2.0, ZERO, 0.0).is_err());
    assert!(model.rate(initial, 2.0, ZERO, f64::NAN).is_err());
    assert!(
        model
            .step_held(initial, 2.0, [100.0, 1e7, 0.0], 0.0, 0.01)
            .is_err()
    );
    assert!(
        model
            .step_held(initial, 2.0, [100.0, 0.0, 0.0], 1e300, 1e300)
            .is_err()
    );
    let expected = model
        .step_held(initial, 2.0, [100.0, 1.0, 0.0], 0.0, 0.001)
        .unwrap();
    assert_eq!(
        expected,
        model
            .step_held(initial, 2.0, [100.0, 1.0, 0.0], 0.0, 0.001)
            .unwrap()
    );
    assert!(FirstHarmonicInflow::new(0.0, RHO).is_err());
    assert!(AzimuthQuadrature::new(516).is_err());
    let loads = AzimuthQuadrature::new(8).unwrap().evaluate_first_harmonic(
        &rotor(),
        RHO,
        frame(),
        [0.0, 0.0, 7.0],
        ZERO,
        field(2.0, 10.0, 0.0),
        60.0,
    );
    assert!(loads.is_err());
}
