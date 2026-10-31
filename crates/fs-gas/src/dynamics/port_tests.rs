//! End ports: caller-owned stage fluxes, the public end-wave and nozzle
//! relations, and transactional local cell sources.

use super::super::tests::gas;
use super::*;
use fs_qty::Length;

const AREA: f64 = 0.7;

fn tube(gas: FrozenGas, n: usize, boundaries: Boundaries) -> Euler1d {
    let q = (0..n)
        .map(|i| {
            let x = (i as f64 + 0.5) / n as f64;
            let bump = (-((x - 0.3) / 0.05).powi(2)).exp();
            gas.conserved(1.0 + 0.2 * bump, 0.3, 10.0 + 3.0 * bump)
                .unwrap()
        })
        .collect();
    Euler1d::new(
        gas,
        Length::new(1.0 / n as f64),
        AREA,
        q,
        boundaries,
        || false,
    )
    .unwrap()
}

fn ends(left: Boundary, right: Boundary) -> Boundaries {
    Boundaries::Ends { left, right }
}

fn at_rest(gas: &FrozenGas, density: f64, velocity: f64, pressure: f64) -> Primitive {
    gas.primitive(gas.conserved(density, velocity, pressure).unwrap(), || {
        false
    })
    .unwrap()
}

#[test]
fn wall_ports_reproduce_reflecting_ends_bit_for_bit_with_ordered_stages() {
    let wall = gas(0.02, 0.0);
    let mut declared = tube(
        gas(0.02, 0.0),
        64,
        ends(Boundary::Reflecting, Boundary::Reflecting),
    );
    let mut ported = tube(
        gas(0.02, 0.0),
        64,
        ends(Boundary::Transmissive, Boundary::Transmissive),
    );
    for _ in 0..60 {
        let start = ported.time().value();
        let mut offered = Vec::new();
        let a = declared
            .step(Time::new(1.0), StepControl::default(), || false)
            .unwrap();
        let b = ported
            .step_with_ports(
                Time::new(1.0),
                StepControl::default(),
                |s| {
                    offered.push((s.stage, s.time, s.dt));
                    Ok([
                        Some(wall.wall_flux(s.left.0, s.left.1, -1.0, &mut || false)?),
                        Some(wall.wall_flux(s.right.0, s.right.1, 1.0, &mut || false)?),
                    ])
                },
                || false,
            )
            .unwrap();
        assert_eq!(a.elapsed, b.elapsed);
        assert_eq!(declared.states(), ported.states());
        let dt = b.elapsed.value();
        assert_eq!(offered, vec![(0, start, dt), (1, start + dt, dt)]);
        assert_eq!(b.left_boundary_change.mass, 0.0);
        assert_eq!(b.right_boundary_change.energy, 0.0);
    }
}

#[test]
fn an_injecting_port_is_the_whole_end_exchange_and_conserves() {
    let source = gas(0.02, 0.0);
    let mut flow = tube(
        gas(0.02, 0.0),
        48,
        ends(Boundary::Transmissive, Boundary::Reflecting),
    );
    let face = at_rest(&source, 1.1, 0.05, 10.4);
    let flux = source.inflow_flux(face, 10.6, 0.0).unwrap();
    assert_eq!(
        flux.energy,
        flux.mass * source.stagnation_enthalpy(10.6).unwrap()
    );
    for _ in 0..20 {
        let before = flow.totals();
        let report = flow
            .step_with_ports(
                Time::new(1.0),
                StepControl::default(),
                |_| Ok([Some(flux), None]),
                || false,
            )
            .unwrap();
        let dt = report.elapsed.value();
        let injected = dt * AREA * flux.mass;
        assert!((report.left_boundary_change.mass - injected).abs() <= 1e-15 * injected);
        assert!(
            (report.left_boundary_change.energy - dt * AREA * flux.energy).abs()
                <= 1e-14 * (dt * AREA * flux.energy).abs()
        );
        assert_eq!(report.right_boundary_change.mass, 0.0);
        let after = flow.totals();
        assert!((after.mass - before.mass - injected).abs() < 3e-13);
        assert!(report.conservation_residual.energy.abs() < 3e-11);
    }
}

#[test]
fn a_refusing_port_leaves_field_and_time_and_periodic_intervals_refuse_ports() {
    let mut flow = tube(
        gas(0.02, 0.0),
        16,
        ends(Boundary::Reflecting, Boundary::Reflecting),
    );
    let states = flow.states().to_vec();
    let mut calls = 0;
    let control = StepControl {
        maximum_halvings: 2,
        ..StepControl::default()
    };
    let refused = flow.step_with_ports(
        Time::new(1.0),
        control,
        |_| {
            calls += 1;
            Err(Error::Invalid("port refused"))
        },
        || false,
    );
    assert!(matches!(
        refused,
        Err(Error::StepRetries { attempts: 3, .. })
    ));
    assert_eq!(calls, 3);
    assert_eq!(flow.states(), &states[..]);
    assert_eq!(flow.time().value(), 0.0);
    let mut periodic = tube(gas(0.02, 0.0), 16, Boundaries::Periodic);
    assert!(matches!(
        periodic.step_with_ports(
            Time::new(1.0),
            StepControl::default(),
            |_| Ok([None; 2]),
            || false
        ),
        Err(Error::Invalid(_))
    ));
}

#[test]
fn end_wave_has_the_acoustic_limit_and_stops_gas_at_the_wall_pressure() {
    let source = gas(0.02, 0.0);
    let inside = at_rest(&source, 1.0, 0.4, 10.0);
    let impedance = inside.density * inside.sound_speed;
    for normal in [-1.0, 1.0] {
        let same = source
            .end_wave(inside, normal, inside.pressure, || false)
            .unwrap();
        assert_eq!(same.primitive.velocity, inside.velocity);
        assert!((same.velocity_per_pressure + normal / impedance).abs() <= 1e-12 / impedance);
        // A weak compression at the end changes velocity by dp/(rho c).
        let dp = 1e-5 * inside.pressure;
        let weak = source
            .end_wave(inside, normal, inside.pressure + dp, || false)
            .unwrap();
        let expected = inside.velocity - normal * dp / impedance;
        assert!((weak.primitive.velocity - expected).abs() <= 1e-4 * dp / impedance);
    }
    // Gas moving at +0.4 toward a right wall stops at the wall solver's pressure.
    let q = source.conserved(1.0, 0.4, 10.0).unwrap();
    let wall = source.wall_flux(q, inside, 1.0, &mut || false).unwrap();
    let mut p = inside.pressure;
    for _ in 0..40 {
        let wave = source.end_wave(inside, 1.0, p, || false).unwrap();
        let step = wave.primitive.velocity / wave.velocity_per_pressure;
        p -= step;
        if step.abs() <= 1e-14 * p {
            break;
        }
    }
    assert!(
        (p - wall.momentum).abs() <= 1e-9 * wall.momentum,
        "{p} {}",
        wall.momentum
    );
}

#[test]
fn nozzle_exit_has_the_incompressible_bernoulli_limit_and_face_flux_is_consistent() {
    let source = gas(0.02, 0.0);
    let donor = at_rest(&source, 1.0, 0.0, 10.0);
    let drop = 1e-4 * donor.pressure;
    let (exit, choked) = source
        .nozzle_exit(donor, donor.pressure - drop, || false)
        .unwrap();
    assert!(!choked);
    let bernoulli = (2.0 * drop / donor.density).sqrt();
    assert!((exit.velocity - bernoulli).abs() <= 1e-4 * bernoulli);
    let back = donor.pressure - drop;
    assert!((exit.pressure - back).abs() <= 1e-14 * back);
    // The exit's own stagnation enthalpy is the donor's.
    let h_exit = source.stagnation_enthalpy(exit.temperature).unwrap();
    let h_donor = source.stagnation_enthalpy(donor.temperature).unwrap();
    let kinetic = 0.5 * exit.velocity * exit.velocity;
    assert!((h_exit + kinetic - h_donor).abs() <= 1e-9 * kinetic);
    // A face carrying one state has the interface flux of that state on both sides.
    let face = at_rest(&source, 1.2, -0.3, 11.0);
    let q = source
        .conserved(face.density, face.velocity, face.pressure)
        .unwrap();
    let euler = source.face_flux(face, 0.0).unwrap();
    let riemann = source.interface(q, face, q, face, &mut || false).unwrap();
    for (a, b) in [
        (euler.mass, riemann.mass),
        (euler.momentum, riemann.momentum),
        (euler.energy, riemann.energy),
    ] {
        assert!((a - b).abs() <= 1e-12 * (1.0 + b.abs()), "{a} {b}");
    }
}

#[test]
fn cell_sources_are_transactional_and_exact() {
    let mut flow = tube(
        gas(0.02, 0.0),
        24,
        ends(Boundary::Reflecting, Boundary::Reflecting),
    );
    let states = flow.states().to_vec();
    let refused = flow.apply_cell_sources(
        |i, q, _| {
            Ok(Conserved {
                density: if i == 17 { -q.density } else { q.density },
                ..q
            })
        },
        || false,
    );
    assert!(refused.is_err());
    assert_eq!(flow.states(), &states[..]);
    let before = flow.totals();
    flow.apply_cell_sources(
        |_, q, p| {
            // Stop each cell and keep its energy: kinetic becomes internal.
            assert!(p.temperature > 0.0);
            Ok(Conserved { momentum: 0.0, ..q })
        },
        || false,
    )
    .unwrap();
    let after = flow.totals();
    assert_eq!(after.mass, before.mass);
    assert_eq!(after.energy, before.energy);
    assert_eq!(after.momentum, 0.0);
    for q in flow.states() {
        let p = flow.gas().primitive(*q, || false).unwrap();
        assert_eq!(p.velocity, 0.0);
    }
}
