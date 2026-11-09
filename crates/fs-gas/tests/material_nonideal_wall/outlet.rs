//! Original-source backpressure waves and conservative field integration.
use super::*;
use fs_gas::dynamics::{FaceFlux, StepReport};
use fs_qty::Pressure;

fn physical(gas: &RealGas, q: Conserved) -> FaceFlux {
    let p = gas.primitive(q, || false).unwrap();
    FaceFlux {
        mass: q.momentum,
        momentum: q.momentum * p.velocity + p.pressure,
        energy: (q.energy + p.pressure) * p.velocity,
        signal_speed: p.velocity.abs() + p.sound_speed,
        used_hll_fallback: false,
    }
}

fn flux_close(actual: FaceFlux, expected: FaceFlux, relative: f64) {
    for (a, b) in [
        (actual.mass, expected.mass),
        (actual.momentum, expected.momentum),
        (actual.energy, expected.energy),
    ] {
        near(a, b, relative * (1.0 + b.abs()));
    }
    assert!(!actual.used_hll_fallback);
    assert!(actual.signal_speed > 0.0 && actual.signal_speed.is_finite());
}

#[test]
fn dense_backpressure_drives_incoming_shock_against_supersonic_discharge() {
    for datum in [0.0, 1e12] {
        let gas = gas([datum; 2]);
        let velocity = 600.0;
        let downstream_velocity = velocity - (869.094_143_566_673 - 319.396_301_986_079_2);
        let expected = gas
            .conserved(
                Density::new(27.210_526_175_865_127),
                Velocity::new(downstream_velocity),
                Temperature::new(632.774_712_247_970_2),
                || false,
            )
            .unwrap();
        let right = gas
            .pressure_outlet_flux(
                state(&gas, velocity, 400.0),
                1.0,
                Pressure::new(6_182_880.252_416_264),
                controls(),
                || false,
            )
            .unwrap();
        flux_close(right, physical(&gas, expected), 2e-8);
        let left = gas
            .pressure_outlet_flux(
                state(&gas, -velocity, 400.0),
                -1.0,
                Pressure::new(6_182_880.252_416_264),
                controls(),
                || false,
            )
            .unwrap();
        near(left.mass, -right.mass, 0.0);
        near(left.momentum, right.momentum, 0.0);
        near(left.energy, -right.energy, 0.0);
    }
}

#[test]
fn fully_outgoing_supersonic_waves_leave_the_interior_flux_unchanged() {
    let gas = gas([0.0; 2]);
    let q = state(&gas, 600.0, 400.0);
    let p = gas.primitive(q, || false).unwrap();
    for pressure in [0.0, 0.01 * p.pressure, p.pressure, 1.01 * p.pressure] {
        flux_close(
            gas.pressure_outlet_flux(q, 1.0, Pressure::new(pressure), controls(), || false)
                .unwrap(),
            physical(&gas, q),
            0.0,
        );
    }
}

fn dilute() -> (RealGas, RealGasWallControl, Conserved) {
    let rho = 1e-5;
    let gas = RealGas::new(
        common::model([0.0; 2], [0.0; 2], [200.0, 2000.0]),
        point(800.0),
        Density::new(rho),
        RealGasDomain {
            density_bounds_kg_m3: [rho * 0.01, rho * 5.0],
            energy: energy(),
        },
        || false,
    )
    .unwrap();
    let mut control = controls();
    control.minimum_rarefaction_density_kg_m3 = rho * 0.1;
    control.compression.momentum_tolerance_pa = 1e-10;
    let q = gas
        .conserved(
            Density::new(rho),
            Velocity::new(100.0),
            Temperature::new(800.0),
            || false,
        )
        .unwrap();
    (gas, control, q)
}

#[test]
fn dilute_unchoked_expansion_converges_to_the_independent_ideal_riemann_invariant() {
    let (gas, control, q) = dilute();
    let p = gas.primitive(q, || false).unwrap();
    let ratio: f64 = 0.8;
    // Authored dilute source has Cp/R=4, hence gamma=4/3. This oracle
    // integrates du=-a d(rho)/rho analytically, independently of PR quadrature.
    let temperature = p.temperature * ratio.powf(0.25);
    let density = p.density * ratio.powf(0.75);
    let velocity = p.velocity + 6.0 * p.sound_speed * (1.0 - ratio.powf(0.125));
    let expected = gas
        .conserved(
            Density::new(density),
            Velocity::new(velocity),
            Temperature::new(temperature),
            || false,
        )
        .unwrap();
    let actual = gas
        .pressure_outlet_flux(q, 1.0, Pressure::new(ratio * p.pressure), control, || false)
        .unwrap();
    flux_close(actual, physical(&gas, expected), 3e-7);
}

#[test]
fn sonic_discharge_saturates_for_lower_backpressure_and_matches_the_dilute_limit() {
    let (gas, control, q) = dilute();
    let p = gas.primitive(q, || false).unwrap();
    let sonic = (6.0 * p.sound_speed + p.velocity) / 7.0;
    let temperature_ratio = (sonic / p.sound_speed).powi(2);
    let density = p.density * temperature_ratio.powi(3);
    let expected = gas
        .conserved(
            Density::new(density),
            Velocity::new(sonic),
            Temperature::new(p.temperature * temperature_ratio),
            || false,
        )
        .unwrap();
    let a = gas
        .pressure_outlet_flux(q, 1.0, Pressure::new(0.1 * p.pressure), control, || false)
        .unwrap();
    let b = gas
        .pressure_outlet_flux(q, 1.0, Pressure::new(0.01 * p.pressure), control, || false)
        .unwrap();
    flux_close(a, physical(&gas, expected), 4e-7);
    flux_close(b, a, 2e-8);
    let vacuum = gas
        .pressure_outlet_flux(q, 1.0, Pressure::new(0.0), control, || false)
        .unwrap();
    flux_close(vacuum, a, 2e-8);
}

#[test]
fn weak_real_gas_backpressure_is_continuous_with_native_acoustic_impedance() {
    let gas = gas([0.0; 2]);
    let q = state(&gas, 100.0, 800.0);
    let p = gas.primitive(q, || false).unwrap();
    for sign in [-1.0, 1.0] {
        let increment = sign * 0.01;
        let flux = gas
            .pressure_outlet_flux(
                q,
                1.0,
                Pressure::new(p.pressure + increment),
                controls(),
                || false,
            )
            .unwrap();
        let velocity = p.velocity - increment / (p.density * p.sound_speed);
        let density = p.density + increment / p.sound_speed.powi(2);
        near(flux.mass, density * velocity, 2e-6);
        // Finite tiny pressure increments must produce a nonidentity response.
        assert_ne!(flux.mass, q.momentum);
    }
}

fn check_balance(report: StepReport) {
    for (before, after, residual) in [
        (
            report.before.mass,
            report.after.mass,
            report.conservation_residual.mass,
        ),
        (
            report.before.momentum,
            report.after.momentum,
            report.conservation_residual.momentum,
        ),
        (
            report.before.energy,
            report.after.energy,
            report.conservation_residual.energy,
        ),
    ] {
        near(
            residual,
            0.0,
            1024.0 * f64::EPSILON * (1.0 + before.abs() + after.abs()),
        );
    }
}

#[test]
fn pressure_outlet_drives_two_stage_discharge_with_complete_conservation_and_atomic_cancel() {
    let gas = gas([0.0; 2]);
    let q = state(&gas, 100.0, 800.0);
    let p = gas.primitive(q, || false).unwrap();
    let mut flow = RealGasEuler1d::new(
        gas,
        Length::new(0.1),
        0.02,
        vec![q; 4],
        RealGasBoundaries::Ends {
            left: RealGasBoundary::Exterior(q),
            right: RealGasBoundary::PressureOutlet {
                pressure_pa: 0.9 * p.pressure,
                control: controls(),
            },
        },
        || false,
    )
    .unwrap();
    let step = StepControl {
        cfl: 0.4,
        reconstruction: Reconstruction::Muscl,
        maximum_halvings: 6,
    };
    for _ in 0..3 {
        let report = flow.step(Time::new(1e-6), step, || false).unwrap();
        check_balance(report);
        assert!(report.after.mass < report.before.mass);
        assert!(-report.right_boundary_change.mass > report.left_boundary_change.mass);
    }
    let saved = flow.states().to_vec();
    let time = flow.time();
    let mut calls = 0;
    assert!(matches!(
        flow.step(Time::new(1e-6), step, || {
            calls += 1;
            calls == 30
        }),
        Err(RealGasFlowError::Cancelled)
    ));
    assert_eq!(flow.states(), saved);
    assert_eq!(flow.time(), time);
}

#[test]
fn pressure_only_boundary_refuses_backflow_and_unsupported_source_paths() {
    let gas = gas([0.0; 2]);
    let q = state(&gas, 50.0, 400.0);
    assert!(
        gas.pressure_outlet_flux(
            q,
            1.0,
            Pressure::new(6_182_880.252_416_264),
            controls(),
            || false
        )
        .is_err()
    );
    let p = gas.primitive(q, || false).unwrap();
    assert!(
        gas.pressure_outlet_flux(
            state(&gas, -1.0, 400.0),
            1.0,
            Pressure::new(p.pressure),
            controls(),
            || false
        )
        .is_err()
    );
    let mut short = controls();
    short.minimum_rarefaction_density_kg_m3 = 9.99;
    assert!(
        gas.pressure_outlet_flux(q, 1.0, Pressure::new(0.1 * p.pressure), short, || false)
            .is_err()
    );
    assert!(
        gas.pressure_outlet_flux(q, 1.0, Pressure::new(0.0), short, || false)
            .is_err()
    );
}

#[test]
fn sonic_face_remains_supported_when_the_remote_fan_is_below_the_source_temperature() {
    let mut energy = energy();
    energy.minimum_temperature_k = 700.0;
    let gas = RealGas::new(
        common::model([0.0; 2], [0.0; 2], [200.0, 2000.0]),
        point(800.0),
        Density::new(10.0),
        RealGasDomain {
            density_bounds_kg_m3: [0.005, 50.0],
            energy,
        },
        || false,
    )
    .unwrap();
    let resting = state(&gas, 0.0, 800.0);
    let p = gas.primitive(resting, || false).unwrap();
    let q = state(&gas, 0.9 * p.sound_speed, 800.0);
    let mut control = controls();
    control.compression.energy = energy;
    let first = gas
        .pressure_outlet_flux(q, 1.0, Pressure::new(0.1 * p.pressure), control, || false)
        .unwrap();
    let lower = gas
        .pressure_outlet_flux(q, 1.0, Pressure::new(0.01 * p.pressure), control, || false)
        .unwrap();
    flux_close(first, lower, 2e-8);
    let vacuum = gas
        .pressure_outlet_flux(q, 1.0, Pressure::new(0.0), control, || false)
        .unwrap();
    flux_close(vacuum, first, 2e-8);
    assert!(first.mass > q.momentum);
}

#[test]
fn weak_backpressure_shock_enters_sonic_and_slightly_supersonic_outlets() {
    let (gas, control, reference) = dilute();
    let p = gas.primitive(reference, || false).unwrap();
    let acoustic_velocity = 0.75 * 32.0 * f64::EPSILON.sqrt() * p.sound_speed;
    for offset in [0.0, 0.25, 2.0] {
        let q = gas
            .conserved(
                Density::new(p.density),
                Velocity::new(p.sound_speed + offset * acoustic_velocity),
                Temperature::new(p.temperature),
                || false,
            )
            .unwrap();
        let inside = gas.primitive(q, || false).unwrap();
        let increment = inside.density * inside.sound_speed * acoustic_velocity;
        let flux = gas
            .pressure_outlet_flux(
                q,
                1.0,
                Pressure::new(inside.pressure + increment),
                control,
                || false,
            )
            .unwrap();
        // Independent gamma=4/3 dilute Rankine--Hugoniot wave speed:
        // u1/a1 = sqrt(1 + (gamma+1)*Delta p/(2*rho1*a1²)).
        let shock_speed = inside.velocity
            - inside.sound_speed
                * (1.0 + (7.0 / 6.0) * acoustic_velocity / inside.sound_speed).sqrt();
        if offset < 1.0 {
            assert!(shock_speed < 0.0);
            // Delta(mass flux) = shock_speed * Delta(rho) < 0. The
            // perturbation is resolved but small enough to use the weak path.
            assert!(flux.mass < q.momentum, "incoming weak shock was ignored");
        } else {
            assert!(shock_speed > 0.0);
            flux_close(flux, physical(&gas, q), 0.0);
        }
    }
}
