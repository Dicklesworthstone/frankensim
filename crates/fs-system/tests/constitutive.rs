//! Nonlinear electrical storage follows its source curve through arbitrary
//! kernel connections. References below use independent piecewise physics.

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_ports::PortKind;
use fs_ports::component::conformance::{Case, Observed, check_component, check_port_tangent};
use fs_ports::component::{Component, PortTangent};
use fs_system::constitutive::{
    CapacitorChargeCurve, ConstitutiveError, InductorFluxCurve, NonlinearCapacitor,
    NonlinearInductor,
};
use fs_system::elements::{Admittance, HeatCapacity};
use fs_system::{System, SystemBuilder, SystemError};

fn with_cx<R>(f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate = CancelGate::new_clock_free();
    let pool = ArenaPool::new(ArenaConfig::default());
    pool.scope(|arena| {
        let key = StreamKey { seed: 41, kernel_id: 37, tile: 0, iteration: 0 };
        f(&Cx::new(&gate, arena, key, Budget::INFINITE, ExecMode::Deterministic))
    })
}

fn knots() -> Vec<(f64, f64)> {
    vec![(-2.0, -6.0), (-1.0, -2.0), (0.0, 0.0), (1.0, 2.0), (2.0, 6.0)]
}

// Constitutive law for this fixture: differential storage is 1/2 up to
// unit charge/flux and 1/4 outside. Integrating it gives this Hamiltonian.
fn effort(x: f64) -> f64 {
    if x.abs() <= 1.0 { 2.0 * x } else { 4.0 * x - 2.0 * x.signum() }
}

fn energy(x: f64) -> f64 {
    if x.abs() <= 1.0 { x * x } else { 2.0 * x * x - 2.0 * x.abs() + 1.0 }
}

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!((actual - expected).abs() <= tolerance, "{actual} vs {expected}, tolerance {tolerance}");
}

fn conforms<C: PortTangent<State = f64>>(component: &C, old: f64, dt: f64, input: f64) {
    with_cx(|cx| {
        let case = Case { state: &old, stale: &17.0, dt_s: dt, inputs: &[input] };
        let observed = check_component(component, &case, cx).unwrap();
        let Observed::Stepped { outputs, .. } = observed else { panic!("unexpected shortening") };
        let next = old + dt * input;
        close(dt * input * outputs[0], energy(next) - energy(old), 2.0e-14);
        close(component.stored_energy_j(&next).unwrap(), energy(next), 2.0e-14);
        check_port_tangent(component, &case, 1.0e-6, 1.0e-9, cx).unwrap();
    });
}

#[test]
fn both_adapters_conform_across_knots_and_keep_the_original_endpoint_law() {
    let capacitor = NonlinearCapacitor::new(CapacitorChargeCurve::try_new(knots()).unwrap());
    let inductor = NonlinearInductor::new(InductorFluxCurve::try_new(knots()).unwrap());
    for (old, dt, input) in [
        (0.2, 0.1, 2.0),
        (1.2, 0.1, -4.0),
        (-1.2, 0.1, 4.0),
        (-1.5, 0.5, 6.0),
        (0.3, 0.1, 0.0),
    ] {
        conforms(&capacitor, old, dt, input);
        conforms(&inductor, old, dt, input);
        let next = old + dt * input;
        close(capacitor.curve().voltage_v(next).unwrap(), effort(next), 1.0e-14);
        close(inductor.curve().current_a(next).unwrap(), effort(next), 1.0e-14);
    }
}

fn rc(dt: f64) -> (f64, f64) {
    // V(q)=4q-2 until q=1 at t=0.4, then V(q)=2q. With G=1/2,
    // q(t)=1/2+(q0-1/2)e^(-2t), followed by e^(-(t-0.4)).
    let initial = 0.5 + 0.5 * 0.8_f64.exp();
    let mut builder = SystemBuilder::new();
    let capacitor = builder.add_with_tangent(
        "capacitor", NonlinearCapacitor::new(CapacitorChargeCurve::try_new(knots()).unwrap()), initial,
    ).unwrap();
    let resistor = builder.add_with_tangent(
        "resistor", Admittance::in_domain(PortKind::ElectricalVoltageCurrent, 0.5).unwrap(), (),
    ).unwrap();
    let thermal = builder.add_with_tangent("body", HeatCapacity::new(10.0).unwrap(), 300.0).unwrap();
    builder.common_effort(builder.port(capacitor, "p").unwrap(), &[builder.port(resistor, "a").unwrap()]).unwrap();
    builder.boundary(builder.port(resistor, "b").unwrap(), Box::new(|_, _, values| values[0] = 0.0)).unwrap();
    let heat = builder.port(resistor, "heat").unwrap();
    builder.bond(builder.port(thermal, "heat").unwrap(), heat).unwrap();
    let mut system = builder.seal().unwrap();
    assert!(system.guess(heat, &[300.0]));
    with_cx(|cx| {
        for _ in 0..(1.0 / dt).round() as usize { system.step(dt, cx).unwrap(); }
    });
    let charge = *system.state::<f64>(capacitor).unwrap();
    let loss = energy(initial) - energy(charge);
    close(system.ledger().nodes[resistor.index()].dissipated_j.value(), loss, 2.0e-11);
    close(system.ledger().nodes[thermal.index()].stored_change_j.value(), loss, 2.0e-11);
    close(10.0 * (system.state::<f64>(thermal).unwrap() - 300.0), loss, 1.0e-9);
    close(system.ledger().energy_residual_j(), 0.0, 2.0e-11);
    (charge, (charge - (-0.6_f64).exp()).abs())
}

#[test]
fn nonlinear_rc_matches_piecewise_exponential_physics_and_returns_its_energy_as_heat() {
    let (_, coarse) = rc(0.01);
    let (charge, fine) = rc(0.005);
    close(charge, (-0.6_f64).exp(), 8.0e-6);
    assert!(fine < 0.3 * coarse, "coarse {coarse}, fine {fine}");
}

#[test]
fn nonlinear_lc_keeps_both_storage_laws_and_matches_the_differential_equations() {
    let (mut q, mut flux) = (1.6, -0.4);
    let initial_energy = energy(q) + energy(flux);
    let mut builder = SystemBuilder::new();
    let capacitor = builder.add_with_tangent(
        "capacitor", NonlinearCapacitor::new(CapacitorChargeCurve::try_new(knots()).unwrap()), q,
    ).unwrap();
    let inductor = builder.add_with_tangent(
        "inductor", NonlinearInductor::new(InductorFluxCurve::try_new(knots()).unwrap()), flux,
    ).unwrap();
    builder.bond(builder.port(capacitor, "p").unwrap(), builder.port(inductor, "p").unwrap()).unwrap();
    let mut system = builder.seal().unwrap();
    with_cx(|cx| { for _ in 0..1000 { system.step(0.002, cx).unwrap(); } });
    // Independent RK4 of dq/dt=-I(flux), d(flux)/dt=V(q), at 200x
    // finer resolution. No native curve or discrete-gradient call is used.
    let h = 1.0e-5;
    for _ in 0..200_000 {
        let a = [-effort(flux), effort(q)];
        let b = [-effort(flux + 0.5*h*a[1]), effort(q + 0.5*h*a[0])];
        let c = [-effort(flux + 0.5*h*b[1]), effort(q + 0.5*h*b[0])];
        let d = [-effort(flux + h*c[1]), effort(q + h*c[0])];
        q += h * (a[0] + 2.0*b[0] + 2.0*c[0] + d[0]) / 6.0;
        flux += h * (a[1] + 2.0*b[1] + 2.0*c[1] + d[1]) / 6.0;
    }
    let got_q = *system.state::<f64>(capacitor).unwrap();
    let got_flux = *system.state::<f64>(inductor).unwrap();
    close(got_q, q, 3.0e-5);
    close(got_flux, flux, 3.0e-5);
    close(energy(got_q) + energy(got_flux), initial_energy, 2.0e-10);
    close(system.ledger().energy_residual_j(), 0.0, 2.0e-10);
}

fn limited_capacitor() -> (System, fs_system::NodeId) {
    let mut builder = SystemBuilder::new();
    let node = builder.add_with_tangent(
        "limited", NonlinearCapacitor::new(CapacitorChargeCurve::try_new(knots()).unwrap()), 0.0,
    ).unwrap();
    builder.boundary(builder.port(node, "p").unwrap(), Box::new(|_, _, values| values[0] = 10.0)).unwrap();
    (builder.seal().unwrap(), node)
}

#[test]
fn source_support_and_cancellation_preserve_the_accepted_system() {
    let (mut system, node) = limited_capacitor();
    with_cx(|cx| {
        assert!(matches!(system.step(1.0, cx), Err(SystemError::Component { .. })));
        assert_eq!(system.time_s().to_bits(), 0.0_f64.to_bits());
        assert_eq!(system.state::<f64>(node).unwrap().to_bits(), 0.0_f64.to_bits());
        close(system.ledger().energy_residual_j(), 0.0, 0.0);
        system.step(0.01, cx).unwrap();
        close(*system.state::<f64>(node).unwrap(), 0.1, 1.0e-16);
        cx.cancel_gate().request();
        let element = NonlinearInductor::new(InductorFluxCurve::try_new(knots()).unwrap());
        let mut next = 7.0;
        let mut output = [8.0];
        assert!(matches!(element.propose(&0.0, 0.01, &[1.0], &mut next, &mut output, cx), Err(ConstitutiveError::Cancelled)));
        close(next, 7.0, 0.0);
        close(output[0], 8.0, 0.0);
    });
}
