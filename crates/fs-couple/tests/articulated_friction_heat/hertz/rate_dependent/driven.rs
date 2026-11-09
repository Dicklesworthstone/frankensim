//! G1/G3/G4: sourced motor-driven rate curves, actual work/heat and thermal time.
#[path = "driven/events.rs"]
mod events;

use super::*;
use fs_couple::friction_heat::hertz::rate_dependent::driven::*;
use fs_mbd::Wrench;
use fs_tribo::sliding_patch::contact::driven::HertzDrivenMode;
use fs_tribo::sliding_patch::rate_dependent::driven::HertzRateDrivenControls;

fn driven_controls() -> RigidHertzRateDrivenHeatControls {
    let original = rate_controls();
    let mut heat = original.heat;
    // The long thermal phase needs a resolvable dimensioned Jacobian stencil;
    // the actual Newton, body, physical heat and work thresholds stay unchanged.
    heat.difference_step_j = 1e-5;
    RigidHertzRateDrivenHeatControls {
        heat,
        quadrature: original.quadrature,
        driven: HertzRateDrivenControls {
            integration: original.transient,
            rest_scaled_rate_tolerance_m_s: 1e-10,
        },
    }
}
fn motor(moment: f64) -> [Wrench; 2] {
    [1.0, -1.0].map(|sign| Wrench {
        force_world: Vec3::ZERO,
        torque_body: Vec3::new(0.0, 0.0, sign * moment),
    })
}
#[allow(clippy::too_many_arguments)]
fn phase<'m>(
    pair: &Pair,
    owner: &SourceThermalFriction,
    thermal: &DepositionSimulation<'m>,
    moment: f64,
    mode: HertzDrivenMode,
    dt: f64,
    resolve: impl Fn(&DepositionSimulation<'m>, f64) -> Result<RigidHertzHeatSource<f64>, &'static str>,
    gate: &CancelGate,
) -> Result<HeatedRigidHertzRateDrivenCandidate<'m, f64>, RigidHertzRateDrivenHeatError<&'static str>>
{
    propose_rigid_hertz_rate_driven_heat(
        &pair.bodies,
        &pair.properties,
        contact(pair, owner),
        motor(moment),
        mode,
        thermal,
        &routes(thermal)[0],
        thermal.time_s(),
        dt,
        driven_controls(),
        &[0.0; 3],
        300.0,
        resolve,
        gate,
    )
}

// The affine Hertz spin torque is A+B*omega, with the true pressure second
// moment2a²/5. Integrate omega, omega², heat and supplied work independently;
// heat is not reconstructed from the production kinetic-energy residual.
fn affine_driven(theta: f64, moment: f64, dt: f64) -> [f64; 5] {
    let a = effort_scale(true) * (MU[0] - SLOPE[0] * theta);
    let b = 10.0 * 0.4 * radius().powi(2) * (VISC + VISC_T * theta);
    let decay = mobility(true) * b;
    let terminal = (moment - a) / b;
    let delta = 3.0 - terminal;
    let e = (-decay * dt).exp();
    let first = -(-decay * dt).exp_m1() / decay;
    let second = -(-2.0 * decay * dt).exp_m1() / (2.0 * decay);
    let angle = terminal * dt + delta * first;
    let square = terminal * terminal * dt + 2.0 * terminal * delta * first + delta * delta * second;
    [
        terminal + delta * e,
        angle,
        a * angle + b * square,
        moment * angle,
        terminal,
    ]
}

#[test]
fn sourced_affine_motor_spin_matches_exponential_motion_and_independent_caloric_root() {
    let pair = Pair::new(0.0, 3.0);
    let owner = rate_owner(false, 320.0);
    let thermal = thermal([300.0; 3], So3::identity());
    let moment = 0.08;
    for dt in [0.2, 2.0] {
        let theta = temperature_root(dt, |t| affine_driven(t, moment, dt)[2]);
        let expected = affine_driven(theta, moment, dt);
        let out = phase(
            &pair,
            &owner,
            &thermal,
            moment,
            HertzDrivenMode::Sliding,
            dt,
            source,
            &CancelGate::new(),
        )
        .unwrap();
        let native = &out.mechanical.mechanical;
        close(out.source, 300.0 + theta);
        close(native.driven.spin_rate_after_rad_s, expected[0]);
        close(native.driven.spin_increment_rad, expected[1]);
        close(out.heat.dissipated_j, expected[2]);
        close(native.external_work_j, expected[3]);
        close(
            native.kinetic_change_j,
            0.5 * (expected[0].powi(2) - 9.0) / mobility(true),
        );
        assert!(3.0 < expected[0] && expected[0] < expected[4]);
        assert!(native.external_work_j > out.heat.dissipated_j);
        assert_eq!(native.driven.mode, HertzDrivenMode::Sliding);
        assert!(
            out.mechanical.adaptive_error_estimates_available || out.mechanical.exact_invariant
        );
        close(out.combined_energy_residual_j, 0.0);
        for (i, t) in thermal_response(dt, theta).0.into_iter().enumerate() {
            close(
                out.thermal.temperature_k(CellId(i as u64)).unwrap(),
                300.0 + t,
            );
        }
        assert_eq!(out.friction_source.owner_identity(), owner.identity());
        assert_eq!(
            out.friction_source.friction_source().state_coordinate("T"),
            Some(out.source)
        );
        let FrictionLaw::VelocityDependent {
            mu_zero,
            slope_per_speed,
            ..
        } = out.friction_source.law().friction_law()
        else {
            panic!("original curve changed")
        };
        close(*mu_zero, MU[0] - SLOPE[0] * theta);
        close(*slope_per_speed, VISC + VISC_T * theta);
    }
}

// Original analytic Gaussian/beta pressure moment, with a separate fixed-grid
// RK4 of moving rate, angle and traction work. No production pressure rule,
// native driven integrator, or native work residual enters this oracle.
fn stribeck_driven(theta: f64, moment: f64, dt: f64) -> [f64; 3] {
    let derivative = |y: [f64; 3]| {
        let torque = stribeck_spin_torque(theta, y[0]);
        [mobility(true) * (moment - torque), y[0], torque * y[0]]
    };
    let mut y = [3.0, 0.0, 0.0];
    let h = dt / 800.0;
    for _ in 0..800 {
        let a = derivative(y);
        let b = derivative(std::array::from_fn(|i| y[i] + 0.5 * h * a[i]));
        let c = derivative(std::array::from_fn(|i| y[i] + 0.5 * h * b[i]));
        let d = derivative(std::array::from_fn(|i| y[i] + h * c[i]));
        y = std::array::from_fn(|i| y[i] + h * (a[i] + 2.0 * b[i] + 2.0 * c[i] + d[i]) / 6.0);
    }
    y
}

#[test]
fn driven_weakening_curve_retains_spatial_speed_and_separate_motor_work() {
    let pair = Pair::new(0.0, 3.0);
    let owner = rate_owner(true, 320.0);
    let thermal = thermal([300.0; 3], So3::identity());
    let (moment, dt) = (0.09, 0.05);
    let theta = temperature_root(dt, |t| stribeck_driven(t, moment, dt)[2]);
    let expected = stribeck_driven(theta, moment, dt);
    let out = phase(
        &pair,
        &owner,
        &thermal,
        moment,
        HertzDrivenMode::Sliding,
        dt,
        source,
        &CancelGate::new(),
    )
    .unwrap();
    let native = &out.mechanical.mechanical;
    close(out.source, 300.0 + theta);
    close(native.driven.spin_rate_after_rad_s, expected[0]);
    close(native.driven.spin_increment_rad, expected[1]);
    close(out.heat.dissipated_j, expected[2]);
    close(native.external_work_j, moment * expected[1]);
    close(out.combined_energy_residual_j, 0.0);
    close(energy(&out.thermal), expected[2]);
    let centre = 3.0 + mobility(true) * (moment - effort_scale(true) * (0.8 - 0.005 * theta)) * dt;
    assert!((native.driven.spin_rate_after_rad_s - centre).abs() > 5e-4);
    assert!(matches!(
        out.friction_source.law().friction_law(),
        FrictionLaw::Stribeck { .. }
    ));
    assert_eq!(out.friction_source.friction_source().properties().len(), 4);
}

// Continuous original three-cell temperatures and affine local-speed spin.
// Face elimination retains BOTH temperature-dependent coefficients. Each term
// in the heat and motor accounts is integrated as its own physical power.
fn continuous_driven(total: f64, steps: usize) -> [f64; 7] {
    let derivative = |y: [f64; 7]| {
        let (g0, g1, g02) = (
            K[0] * AREA / A,
            K[1] * AREA / A,
            AREA / (A / K[0] + A / K[2]),
        );
        let dry_power = effort_scale(true) * y[0];
        let viscous_power = 10.0 * 0.4 * radius().powi(2) * y[0] * y[0];
        let theta = (g0 * y[1] + g1 * y[2] + MU[0] * dry_power + VISC * viscous_power)
            / (g0 + g1 + SLOPE[0] * dry_power - VISC_T * viscous_power);
        let torque = effort_scale(true) * (MU[0] - SLOPE[0] * theta)
            + 10.0 * 0.4 * radius().powi(2) * (VISC + VISC_T * theta) * y[0];
        [
            mobility(true) * (0.08 - torque),
            (g0 * (theta - y[1]) + g02 * (y[3] - y[1])) / CAP[0],
            g1 * (theta - y[2]) / CAP[1],
            g02 * (y[1] - y[3]) / CAP[2],
            y[0],
            torque * y[0],
            0.08 * y[0],
        ]
    };
    let mut y = [3.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0];
    let h = total / steps as f64;
    for _ in 0..steps {
        let a = derivative(y);
        let b = derivative(std::array::from_fn(|i| y[i] + 0.5 * h * a[i]));
        let c = derivative(std::array::from_fn(|i| y[i] + 0.5 * h * b[i]));
        let d = derivative(std::array::from_fn(|i| y[i] + h * c[i]));
        y = std::array::from_fn(|i| y[i] + h * (a[i] + 2.0 * b[i] + 2.0 * c[i] + d[i]) / 6.0);
    }
    y
}

#[test]
fn full_rate_driven_thermal_march_refines_continuous_coefficients_heat_and_work() {
    let pair = Pair::new(0.0, 3.0);
    let owner = rate_owner(false, 320.0);
    let thermal = thermal([300.0; 3], So3::identity());
    let total = 0.4;
    let reference = continuous_driven(total, 8_000);
    let check = continuous_driven(total, 4_000);
    for i in 0..7 {
        assert!((reference[i] - check[i]).abs() < 1e-9);
    }
    let mut errors = Vec::new();
    for scale in [1.0, 0.0625] {
        let mut control = march_controls();
        control.tolerance.linear_momentum_ns *= scale;
        control.tolerance.angular_momentum_nms *= scale;
        control.tolerance.cell_energy_j *= scale;
        control.tolerance.scaled_travel_m *= scale;
        control.tolerance.dissipated_j *= scale;
        let out = propose_rigid_hertz_rate_driven_heat_march(
            &pair.bodies,
            &pair.properties,
            contact(&pair, &owner),
            motor(0.08),
            HertzDrivenMode::Sliding,
            &thermal,
            &routes(&thermal)[0],
            total,
            driven_controls(),
            control,
            300.0,
            source,
            &CancelGate::new(),
        )
        .unwrap();
        assert_eq!(out.mode_after, HertzDrivenMode::Sliding);
        let m = out.march;
        let HertzHeatMarchMotion::RateDriven(last) = &m.phases.last().unwrap().mechanical else {
            panic!("wrong native owner")
        };
        let rate_error = (last.mechanical.driven.spin_rate_after_rad_s - reference[0]).abs();
        let temperature_error = (0..3)
            .map(|i| {
                (m.thermal.temperature_k(CellId(i as u64)).unwrap() - 300.0 - reference[i + 1])
                    .abs()
            })
            .fold(0.0_f64, f64::max);
        errors.push([
            rate_error,
            temperature_error,
            (m.dissipated_j - reference[5]).abs(),
            (m.external_mechanical_work_j - reference[6]).abs(),
        ]);
        close(m.external_mechanical_work_j, 0.08 * m.spin_increment_rad);
        close(m.thermal_energy_change_j, m.dissipated_j);
        for residual in m.energy_residual_j {
            assert!(residual.abs() < control.energy_tolerance_j);
        }
        for p in &m.phases {
            assert_eq!(p.friction_source.owner_identity(), owner.identity());
            close(
                p.friction_source
                    .friction_source()
                    .property("slope")
                    .unwrap()
                    .value_si(),
                VISC + VISC_T * (p.source - 300.0),
            );
            close(
                p.friction_source
                    .friction_source()
                    .property("kinetic")
                    .unwrap()
                    .value_si(),
                MU[0] - SLOPE[0] * (p.source - 300.0),
            );
        }
        eprintln!(
            "rate-driven thermal: scale={scale}, phases={}, rejected={}, trials={}, errors={:?}",
            m.phases.len(),
            m.rejected_attempts,
            m.thermal_trials,
            errors.last().unwrap()
        );
    }
    for i in 0..4 {
        assert!(errors[1][i] < 0.55 * errors[0][i], "{errors:?}");
    }
    assert!(
        errors[1][0] < 2e-5 && errors[1][1] < 2e-4 && errors[1][2] < 2e-6 && errors[1][3] < 2e-7,
        "{errors:?}"
    );
    assert_eq!(thermal.time_s(), 0.0);
}

#[test]
fn rate_driven_static_conduction_and_late_source_refusals_preserve_joint_state() {
    let owner = rate_owner(false, 320.0);
    let resting = Pair::new(0.0, 0.0);
    let cooling = thermal([310.0, 290.0, 300.0], So3::identity());
    let held = phase(
        &resting,
        &owner,
        &cooling,
        0.07,
        HertzDrivenMode::Resting,
        0.02,
        source,
        &CancelGate::new(),
    )
    .unwrap();
    assert_eq!(
        held.mechanical.mechanical.driven.mode,
        HertzDrivenMode::Resting
    );
    assert_eq!(held.heat.dissipated_j, 0.0);
    close(held.mechanical.mechanical.external_work_j, 0.0);
    assert!(held.heat.interfaces[0].body_heat_j[0] < 0.0);
    assert!(held.heat.interfaces[0].body_heat_j[1] > 0.0);

    let pair = Pair::new(0.0, 3.0);
    let thermal = thermal([300.0; 3], So3::identity());
    let before = thermal.cells().clone();
    let narrow = rate_owner(false, 300.01);
    assert!(matches!(
        phase(
            &pair,
            &narrow,
            &thermal,
            0.08,
            HertzDrivenMode::Sliding,
            0.02,
            source,
            &CancelGate::new()
        ),
        Err(RigidHertzHeatError::FrictionSource(_))
    ));
    let run = |duration,
               resolve: &dyn Fn(
        &DepositionSimulation<'_>,
        f64,
    ) -> Result<RigidHertzHeatSource<f64>, &'static str>,
               gate: &CancelGate| {
        let mut control = march_controls();
        control.initial_step_s = 0.02;
        control.maximum_step_s = 0.02;
        propose_rigid_hertz_rate_driven_heat_march(
            &pair.bodies,
            &pair.properties,
            contact(&pair, &owner),
            motor(0.08),
            HertzDrivenMode::Sliding,
            &thermal,
            &routes(&thermal)[0],
            duration,
            driven_controls(),
            control,
            300.0,
            resolve,
            gate,
        )
    };
    let first = run(0.02, &source, &CancelGate::new()).unwrap();
    let reference = run(0.06, &source, &CancelGate::new()).unwrap();
    let calls = Cell::new(0);
    let refused = run(
        0.06,
        &|state, t| {
            calls.set(calls.get() + 1);
            if calls.get() > first.march.thermal_trials {
                Err("late rate-driven source")
            } else {
                source(state, t)
            }
        },
        &CancelGate::new(),
    );
    assert!(matches!(
        refused,
        Err(HertzHeatMarchError::RateDriven(
            RigidHertzHeatError::Source("late rate-driven source")
        ))
    ));
    assert_eq!(calls.get(), first.march.thermal_trials + 1);
    let gate = CancelGate::new();
    calls.set(0);
    assert!(
        run(
            0.06,
            &|state, t| {
                calls.set(calls.get() + 1);
                if calls.get() > first.march.thermal_trials {
                    gate.request();
                }
                source(state, t)
            },
            &gate
        )
        .is_err()
    );
    assert_eq!(thermal.cells(), &before);
    assert_eq!(thermal.time_s(), 0.0);
    assert_eq!(pair.bodies, Pair::new(0.0, 3.0).bodies);
    let retry = run(0.06, &source, &CancelGate::new()).unwrap();
    assert_eq!(retry.mode_after, reference.mode_after);
    assert_eq!(retry.march.bodies_after, reference.march.bodies_after);
    assert_eq!(retry.march.thermal.cells(), reference.march.thermal.cells());
    assert_eq!(
        retry.march.external_mechanical_work_j,
        reference.march.external_mechanical_work_j
    );
}
