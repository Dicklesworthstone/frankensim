//! G1/G3/G4: original temperature-resolved rate curves, real heat, and retry.
#[path = "rate_dependent/driven.rs"]
mod driven;

use super::super::{bind_dry_interface_system_card, curve};
use super::*;
use fs_contact::interface_binding::thermal_friction::ThermalFrictionAxes;
use fs_contact::interface_binding::thermal_friction::rate_dependent::{
    SourceThermalFriction, ThermalFrictionCurve,
};
use fs_couple::friction_heat::hertz::march::*;
use fs_couple::friction_heat::hertz::rate_dependent::*;
use fs_matdb::{InterfaceSystemCard, PropertyKey, PropertyValue, QueryPoint};
use fs_material::state_point::{
    MaterialPropertySelection, ScalarAdmissibility, ScalarPropertyRequirement,
};
use fs_mbd::tangential_friction::hertz::transient::HertzTransientControls;
use fs_qty::Dims;
use fs_tribo::sliding_patch::rate_dependent::HertzFrictionQuadrature;
use fs_tribo::sliding_patch::{HertzSlidingPatch, HertzSlidingPatchParameters};
use fs_tribo::{FrictionLaw, InputAuthority};

const VISC: f64 = 3.0;
const VISC_T: f64 = 0.05;
const INV_SPEED: Dims = Dims([-1, 0, 1, 0, 0, 0]);
const SPEED_DIM: Dims = Dims([1, 0, -1, 0, 0, 0]);

fn rate_owner(stribeck: bool, maximum_temperature: f64) -> SourceThermalFriction {
    let original = law(0);
    let card = original.interface().card();
    let mut claims = card.claims().clone();
    for (name, dims, values) in [
        (
            "slope",
            INV_SPEED,
            if stribeck {
                [0.2, 0.4]
            } else {
                [VISC - 20.0 * VISC_T, VISC + 20.0 * VISC_T]
            },
        ),
        ("speed", SPEED_DIM, [0.08, 0.12]),
    ] {
        let mut claim = curve(name, 0, values);
        claim.key = PropertyKey::new(name, dims);
        if let PropertyValue::Curve { dims: unit, .. } = &mut claim.value {
            *unit = dims;
        }
        claim.validity =
            fs_evidence::ValidityDomain::unconstrained().with("T", 280.0, maximum_temperature);
        claims.insert_claim(claim).unwrap();
    }
    let card = InterfaceSystemCard::assemble(
        card.surface_a().clone(),
        card.surface_b().clone(),
        card.context().clone(),
        claims,
        vec![],
    )
    .unwrap();
    let requirement = |name, dims| {
        ScalarPropertyRequirement::try_new(name, dims, ScalarAdmissibility::NonNegative).unwrap()
    };
    let selected = if stribeck {
        ThermalFrictionCurve::Stribeck {
            static_mu: requirement("static", Dims::NONE),
            kinetic_mu: requirement("kinetic", Dims::NONE),
            characteristic_speed: requirement("speed", SPEED_DIM),
            viscous_per_speed: requirement("slope", INV_SPEED),
        }
    } else {
        ThermalFrictionCurve::VelocityDependent {
            static_mu: requirement("static", Dims::NONE),
            mu_zero: requirement("kinetic", Dims::NONE),
            slope_per_speed: requirement("slope", INV_SPEED),
        }
    };
    SourceThermalFriction::resolve(
        bind_dry_interface_system_card(&card, InputAuthority::SyntheticFixture).unwrap(),
        QueryPoint::new().with("T", 300.0).unwrap(),
        selected,
        MaterialPropertySelection::SingleClaimOnly,
        ThermalFrictionAxes {
            temperature: "T".into(),
            nominal_pressure: None,
        },
    )
    .unwrap()
}
fn contact<'a>(pair: &Pair, owner: &'a SourceThermalFriction) -> RigidHertzRateHeatContact<'a> {
    RigidHertzRateHeatContact {
        binding: pair.binding,
        frame: TangentFrame::new([0.0, 0.0, 1.0], [1.0, 0.0, 0.0]).unwrap(),
        source: owner,
        normal_force_n: 10.0,
        contact_radius_m: radius(),
        nominal_pressure_pa: None,
    }
}
fn rate_controls() -> RigidHertzRateHeatControls {
    RigidHertzRateHeatControls {
        heat: controls(),
        quadrature: HertzFrictionQuadrature {
            radial_panels: 32,
            angular_panels: 8,
        },
        transient: HertzTransientControls {
            maximum_step_s: 0.02,
            minimum_step_s: 1e-15,
            maximum_attempts: 100_000,
            maximum_evaluations: 400_001,
            scaled_rate_tolerance_m_s: 1e-11,
            scaled_impulse_tolerance_ns: 1e-10,
            scaled_travel_tolerance_m: 1e-10,
            heat_tolerance_j: 1e-12,
            stopping_time_tolerance_s: 1e-7,
        },
    }
}
fn phase<'m>(
    pair: &Pair,
    owner: &SourceThermalFriction,
    thermal: &DepositionSimulation<'m>,
    route: &FrictionHeatRoute,
    dt: f64,
    resolve: impl Fn(&DepositionSimulation<'m>, f64) -> Result<RigidHertzHeatSource<f64>, &'static str>,
    gate: &CancelGate,
) -> Result<HeatedRigidHertzRateCandidate<'m, f64>, RigidHertzRateHeatError<&'static str>> {
    propose_rigid_hertz_rate_heat(
        &pair.bodies,
        &pair.properties,
        contact(pair, owner),
        thermal,
        route,
        thermal.time_s(),
        dt,
        rate_controls(),
        &[0.0; 3],
        300.0,
        resolve,
        gate,
    )
}

// Independent affine-rate scalar motion. The pure-spin viscous moment uses
// the original Hertz pressure second moment 2a²/5, not centre slip speed.
fn affine_motion(spin: bool, rate: f64, theta: f64, dt: f64) -> (f64, f64, f64) {
    let mu = MU[0] - SLOPE[0] * theta;
    let k = VISC + VISC_T * theta;
    let a = mobility(spin) * effort_scale(spin) * mu;
    let b = mobility(spin) * 10.0 * k * if spin { 0.4 * radius().powi(2) } else { 1.0 };
    let stop = (b * rate / a).ln_1p() / b;
    let moving = dt.min(stop);
    let fraction = -(-b * moving).exp_m1();
    let final_rate = if dt >= stop {
        0.0
    } else {
        rate - (rate + a / b) * fraction
    };
    let travel = (rate + a / b) * fraction / b - a * moving / b;
    (final_rate, travel, stop)
}
fn thermal_response(dt: f64, theta: f64) -> ([f64; 3], f64) {
    let d0 = dt * K[0] * AREA / A;
    let d1 = dt * K[1] * AREA / A;
    let d02 = dt * AREA / (A / K[0] + A / K[2]);
    let c0 = CAP[0] + d02 * CAP[2] / (CAP[2] + d02);
    let t0 = d0 * theta / (c0 + d0);
    (
        [t0, d1 * theta / (CAP[1] + d1), d02 * t0 / (CAP[2] + d02)],
        (d0 * c0 / (c0 + d0) + d1 * CAP[1] / (CAP[1] + d1)) * theta,
    )
}
fn temperature_root(dt: f64, heat: impl Fn(f64) -> f64) -> f64 {
    let (mut low, mut high) = (0.0, 20.0);
    assert!(thermal_response(dt, high).1 > heat(high));
    for _ in 0..60 {
        let mid = 0.5 * (low + high);
        if thermal_response(dt, mid).1 < heat(mid) {
            low = mid;
        } else {
            high = mid;
        }
    }
    0.5 * (low + high)
}

#[test]
fn sourced_affine_slide_and_spin_heat_match_independent_motion_and_thermal_roots() {
    for (spin, initial, durations) in [(false, 0.5, [0.01, 0.3]), (true, 3.0, [0.2, 12.0])] {
        let pair = Pair::new(
            if spin { 0.0 } else { initial },
            if spin { initial } else { 0.0 },
        );
        let owner = rate_owner(false, 320.0);
        let thermal = thermal([300.0; 3], So3::identity());
        let route = routes(&thermal)[0].clone();
        for dt in durations {
            let heat = |theta| {
                let end = affine_motion(spin, initial, theta, dt).0;
                0.5 * (initial * initial - end * end) / mobility(spin)
            };
            let theta = temperature_root(dt, heat);
            let expected = affine_motion(spin, initial, theta, dt);
            let out = phase(
                &pair,
                &owner,
                &thermal,
                &route,
                dt,
                source,
                &CancelGate::new(),
            )
            .unwrap_or_else(|error| panic!("affine Hertz thermal phase: spin={spin}, dt={dt}, initial temperature=300 K: {error:?}"));
            let motion = &out.mechanical.mechanical;
            close(out.source, 300.0 + theta);
            close(out.heat.dissipated_j, heat(theta));
            close(
                if spin {
                    motion.spin_rad_s[1]
                } else {
                    motion.slip_velocity_m_s[1][0]
                },
                expected.0,
            );
            let travel = if spin {
                motion.transient.spin_increment_rad
            } else {
                motion.transient.slip_increment_m[0]
            };
            assert!(
                (travel - expected.1).abs() < 2e-8,
                "{travel} vs {}",
                expected.1
            );
            for (i, t) in thermal_response(dt, theta).0.into_iter().enumerate() {
                close(
                    out.thermal.temperature_k(CellId(i as u64)).unwrap(),
                    300.0 + t,
                );
            }
            close(energy(&out.thermal) - energy(&thermal), heat(theta));
            close(out.combined_energy_residual_j, 0.0);
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
                panic!("curve changed");
            };
            close(*mu_zero, MU[0] - SLOPE[0] * theta);
            close(*slope_per_speed, VISC + VISC_T * theta);
            if dt > expected.2 {
                let stop = motion.transient.stopping_time_bracket_s.unwrap();
                assert!(stop[0] <= expected.2 + 1e-7 && expected.2 <= stop[1] + 1e-7);
                assert_eq!(motion.transient.velocity_after_m_s, [0.0; 2]);
                assert_eq!(motion.transient.spin_rate_after_rad_s, 0.0);
            }
        }
    }
}

// For pure spin, the Gaussian Hertz moment is
// (3π/16) * 1F1(3/2;3;-(ωa/v0)^2). Its rapidly convergent beta-integral
// series is independent of the production radial Gauss rule.
fn stribeck_spin(theta: f64, dt: f64) -> f64 {
    let torque = |w: f64| -mobility(true) * stribeck_spin_torque(theta, w);
    let h = dt / 800.0;
    let mut w = 3.0;
    for _ in 0..800 {
        let a = torque(w);
        let b = torque(w + 0.5 * h * a);
        let c = torque(w + 0.5 * h * b);
        let d = torque(w + h * c);
        w += h * (a + 2.0 * b + 2.0 * c + d) / 6.0;
    }
    w
}

fn stribeck_spin_torque(theta: f64, w: f64) -> f64 {
    let z = w * radius() / (0.1 + 0.001 * theta);
    let (mut term, mut sum) = (1.0, 1.0);
    for n in 0..64 {
        term *= -z * z * (n as f64 + 1.5) / ((n as f64 + 3.0) * (n as f64 + 1.0));
        sum += term;
        if term.abs() < 1e-18 {
            break;
        }
    }
    let mu_s = 0.8 - 0.005 * theta;
    let mu_k = MU[0] - SLOPE[0] * theta;
    let k = 0.3 + 0.005 * theta;
    10.0 * (radius() * 3.0 * core::f64::consts::PI / 16.0 * (mu_k + (mu_s - mu_k) * sum)
        + 0.4 * radius().powi(2) * k * w)
}

#[test]
fn temperature_resolved_stribeck_spin_uses_the_entire_pressure_footprint() {
    let pair = Pair::new(0.0, 3.0);
    let owner = rate_owner(true, 320.0);
    let thermal = thermal([300.0; 3], So3::identity());
    let route = routes(&thermal)[0].clone();
    let dt = 0.02;
    let heat = |theta| {
        let end = stribeck_spin(theta, dt);
        0.5 * (9.0 - end * end) / mobility(true)
    };
    let theta = temperature_root(dt, heat);
    let out = phase(
        &pair,
        &owner,
        &thermal,
        &route,
        dt,
        source,
        &CancelGate::new(),
    )
    .unwrap();
    close(out.source, 300.0 + theta);
    close(
        out.mechanical.mechanical.spin_rad_s[1],
        stribeck_spin(theta, dt),
    );
    close(out.heat.dissipated_j, heat(theta));
    let centre_coefficient = 0.8 - 0.005 * theta;
    let centre_speed_surrogate =
        3.0 - mobility(true) * effort_scale(true) * centre_coefficient * dt;
    assert!((out.mechanical.mechanical.spin_rad_s[1] - centre_speed_surrogate).abs() > 5e-4);
    assert!(matches!(
        out.friction_source.law().friction_law(),
        FrictionLaw::Stribeck { .. }
    ));
    assert_eq!(out.friction_source.friction_source().properties().len(), 4);
    close(energy(&out.thermal), heat(theta));
}

fn continuous_rate_reference(total: f64) -> [f64; 8] {
    let unit = HertzSlidingPatch::new(HertzSlidingPatchParameters {
        normal_force_n: 10.0,
        contact_radius_m: radius(),
        friction_coefficient: 1.0,
    })
    .unwrap();
    let derivative = |y: [f64; 8]| {
        let wrench = unit.evaluate([y[0], 0.0], y[1]).unwrap();
        let j = 0.4 * radius().powi(2);
        let quadratic = y[0] * y[0] + j * y[1] * y[1];
        let p = wrench.dissipated_power_w;
        let (g0, g1, g02) = (
            K[0] * AREA / A,
            K[1] * AREA / A,
            AREA / (A / K[0] + A / K[2]),
        );
        let theta = (g0 * y[2] + g1 * y[3] + MU[0] * p + 10.0 * VISC * quadratic)
            / (g0 + g1 + SLOPE[0] * p - 10.0 * VISC_T * quadratic);
        let (mu, k) = (MU[0] - SLOPE[0] * theta, VISC + VISC_T * theta);
        [
            mobility(false) * (mu * wrench.force_n[0] - 10.0 * k * y[0]),
            mobility(true) * (mu * wrench.spin_moment_nm - 10.0 * k * j * y[1]),
            (g0 * (theta - y[2]) + g02 * (y[4] - y[2])) / CAP[0],
            g1 * (theta - y[3]) / CAP[1],
            g02 * (y[2] - y[4]) / CAP[2],
            y[0],
            y[1],
            mu * p + 10.0 * k * quadratic,
        ]
    };
    let mut y = [0.1, 3.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0];
    let dt = total / 20_000.0;
    for _ in 0..20_000 {
        let a = derivative(y);
        let b = derivative(std::array::from_fn(|i| y[i] + 0.5 * dt * a[i]));
        let c = derivative(std::array::from_fn(|i| y[i] + 0.5 * dt * b[i]));
        let d = derivative(std::array::from_fn(|i| y[i] + dt * c[i]));
        y = std::array::from_fn(|i| y[i] + dt * (a[i] + 2.0 * b[i] + 2.0 * c[i] + d[i]) / 6.0);
    }
    y
}
fn march_controls() -> HertzHeatMarchControls {
    HertzHeatMarchControls {
        initial_step_s: 0.08,
        minimum_step_s: 1e-7,
        maximum_step_s: 0.08,
        maximum_attempts: 2_000,
        maximum_accepted_phases: 4_000,
        maximum_thermal_trials: 200_000,
        tolerance: HertzHeatMarchTolerance {
            linear_momentum_ns: 5e-6,
            angular_momentum_nms: 5e-7,
            cell_energy_j: 5e-7,
            scaled_travel_m: 5e-7,
            dissipated_j: 5e-7,
            relative: 0.0,
        },
        energy_tolerance_j: 5e-8,
    }
}

#[test]
fn original_full_rate_curve_march_tracks_continuous_thermal_feedback() {
    let pair = Pair::new(0.1, 3.0);
    let owner = rate_owner(false, 320.0);
    let thermal = thermal([300.0; 3], So3::identity());
    let route = routes(&thermal)[0].clone();
    let reference = continuous_rate_reference(0.08);
    let out = propose_rigid_hertz_rate_heat_march(
        &pair.bodies,
        &pair.properties,
        contact(&pair, &owner),
        &thermal,
        &route,
        0.08,
        rate_controls(),
        march_controls(),
        300.0,
        source,
        &CancelGate::new(),
    )
    .unwrap();
    let HertzHeatMarchMotion::RateDependent(m) = &out.phases.last().unwrap().mechanical else {
        panic!("wrong owner");
    };
    assert!((m.mechanical.slip_velocity_m_s[1][0] - reference[0]).abs() < 1e-5);
    assert!((m.mechanical.spin_rad_s[1] - reference[1]).abs() < 1e-4);
    assert!((out.slip_increment_m[0] - reference[5]).abs() < 2e-6);
    assert!((out.spin_increment_rad - reference[6]).abs() < 1e-5);
    assert!((out.dissipated_j - reference[7]).abs() < 1e-5);
    for i in 0..3 {
        let error =
            (out.thermal.temperature_k(CellId(i as u64)).unwrap() - 300.0 - reference[i + 2]).abs();
        assert!(error < 2e-4, "cell {i}: {error}");
    }
    close(out.thermal_energy_change_j, out.dissipated_j);
    close(out.kinetic_energy_change_j, -out.dissipated_j);
    assert!(out.rejected_attempts > 0);
    for p in &out.phases {
        assert_eq!(p.friction_source.owner_identity(), owner.identity());
        assert_eq!(
            p.friction_source.friction_source().state_coordinate("T"),
            Some(p.source)
        );
        close(
            p.friction_source
                .friction_source()
                .property("slope")
                .unwrap()
                .value_si(),
            VISC + VISC_T * (p.source - 300.0),
        );
    }
    eprintln!(
        "rate/thermal march: phases={}, rejected={}, trials={}",
        out.phases.len(),
        out.rejected_attempts,
        out.thermal_trials
    );
}

#[test]
fn rate_source_domains_late_failure_and_cancel_leave_original_owners_retryable() {
    let pair = Pair::new(0.1, 3.0);
    let owner = rate_owner(false, 320.0);
    let thermal = thermal([300.0; 3], So3::identity());
    let before = thermal.cells().clone();
    let route = routes(&thermal)[0].clone();
    let reference = phase(
        &pair,
        &owner,
        &thermal,
        &route,
        0.02,
        source,
        &CancelGate::new(),
    )
    .unwrap();
    let narrow = rate_owner(false, 300.01);
    assert!(matches!(
        phase(
            &pair,
            &narrow,
            &thermal,
            &route,
            0.02,
            source,
            &CancelGate::new()
        ),
        Err(RigidHertzHeatError::FrictionSource(_))
    ));
    let calls = Cell::new(0usize);
    let failed = phase(
        &pair,
        &owner,
        &thermal,
        &route,
        0.02,
        |state, t| {
            calls.set(calls.get() + 1);
            if calls.get() == reference.trials {
                Err("late rate/heat source")
            } else {
                source(state, t)
            }
        },
        &CancelGate::new(),
    );
    assert!(matches!(
        failed,
        Err(RigidHertzHeatError::Source("late rate/heat source"))
    ));
    let cancel = CancelGate::new();
    calls.set(0);
    assert!(
        phase(
            &pair,
            &owner,
            &thermal,
            &route,
            0.02,
            |state, t| {
                calls.set(calls.get() + 1);
                if calls.get() == reference.trials {
                    cancel.request();
                }
                source(state, t)
            },
            &cancel
        )
        .is_err()
    );
    assert_eq!(thermal.time_s(), 0.0);
    assert_eq!(thermal.cells(), &before);
    assert_eq!(pair.bodies, Pair::new(0.1, 3.0).bodies);
    let retry = phase(
        &pair,
        &owner,
        &thermal,
        &route,
        0.02,
        source,
        &CancelGate::new(),
    )
    .unwrap();
    assert_eq!(retry.bodies_after, reference.bodies_after);
    assert_eq!(retry.thermal.cells(), reference.thermal.cells());
    assert_eq!(retry.friction_source, reference.friction_source);
}
