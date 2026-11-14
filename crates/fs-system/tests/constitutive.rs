//! Nonlinear electrical storage follows its source curve through arbitrary
//! kernel connections. References below use independent piecewise physics.

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_ports::PortKind;
use fs_ports::component::conformance::{
    Case, Failure, Observed, Reversal, ReversalAllowance, check_component, check_port_tangent,
    check_time_reversal,
};
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

/// Unit roundoff of f64.
const U: f64 = f64::EPSILON / 2.0;

/// `n u / (1 - n u)`: a product of `n` factors `(1 + d)`, `|d| <= u`, lies
/// within it of 1.
fn gamma(n: u32) -> f64 {
    let nu = f64::from(n) * U;
    nu / (1.0 - nu)
}

/// An adapter there and back under `reverse`, on a curve whose knots, segment
/// slopes and knot energies are exact in binary (this file's are small
/// integers). The owner then evaluates the piecewise-linear `e` itself, with
/// largest slope `L` over `k` segments and `|e(x)| <= L |x|`, and its exact
/// integral `H`. The conjugate mean `D(a, b) = (H(b) - H(a)) / (b - a)` is
/// symmetric in its ends and moves by at most `L / 2` per unit of either
/// (`|e(b) - D(a, b)| <= L |b - a| / 2`).
///
/// Rounding, counted from fs-phs's scalar curve, every operation exact times
/// `(1 + d)`, `|d| <= u`, with `m` the rounded midpoint, which lies between
/// `a` and `b`:
/// - `effort(m)` is within `gamma(2) L |m|`, its offset from the anchor knot
///   being at most `|m|`;
/// - in one segment the correction `(s / 2) ((a - m) + (b - m))` is within
///   `gamma(3) L |b - a| / 2`;
/// - across segments each `relative_energy` walk from `m`, of at most `k`
///   sub-steps, is within `gamma(2 k + 6) L` times its length squared, so the
///   correction, their difference over `b - a`, is within
///   `gamma(2 k + 8) L |b - a|`;
/// - so `fl D(a, b)` is within
///   `rho(a, b) = L (gamma(3) max(|a|, |b|) + gamma(2 k + 9) |b - a|)`;
/// - `fl H(q)`, a knot energy plus a trapezoid, is within `gamma(5) L q^2`.
///
/// The increment `t = fl(dt u)` is the same both ways up to its sign, so the
/// coordinate that should return (the capacitor's `q2`, the inductor's
/// `-flux2`) is `fl(fl(x0 + t) - t)`, within `delta = gamma(2) (|x0| + |x1|)`
/// of `x0`. Hence the bounds: observed `delta`; output
/// `rho(x0, x1) + rho(R(x1), x2) + L delta / 2`; energy
/// `L max(|x0|, |x2|) delta + gamma(5) L (x0^2 + x2^2)`. For the inductor the
/// output relation also needs an odd curve, `D(-a, -b) = -D(a, b)`.
fn adapter_reverses<C>(
    adapter: &C,
    knots: &[(f64, f64)],
    old: f64,
    dt: f64,
    input: f64,
    reverse: &dyn Fn(&f64) -> f64,
) -> Result<Reversal, Failure<ConstitutiveError>>
where
    C: Component<State = f64, Error = ConstitutiveError>,
{
    let largest = knots
        .windows(2)
        .map(|pair| (pair[1].1 - pair[0].1) / (pair[1].0 - pair[0].0))
        .fold(0.0_f64, f64::max);
    let segments = u32::try_from(knots.len() - 1).unwrap();
    let rho = |a: f64, b: f64| {
        largest * (gamma(3) * a.abs().max(b.abs()) + gamma(2 * segments + 9) * (b - a).abs())
    };
    let observe = |state: &f64| vec![*state];
    let allowance = |start: &f64, forward_end: &f64, back_end: &f64| {
        let (x0, x1, x2) = (*start, *forward_end, *back_end);
        let delta = gamma(2) * (x0.abs() + x1.abs());
        // The backward step runs from R(x1) to x2.
        let output = rho(x0, x1) + rho(reverse(&x1), x2) + 0.5 * largest * delta;
        ReversalAllowance {
            outputs: vec![output],
            observed: vec![delta],
            energy_j: largest * x0.abs().max(x2.abs()) * delta
                + gamma(5) * largest * (x0 * x0 + x2 * x2),
        }
    };
    let inputs = [input];
    let case = Case {
        state: &old,
        stale: &17.0,
        dt_s: dt,
        inputs: &inputs,
    };
    with_cx(|cx| check_time_reversal(adapter, &case, reverse, &observe, &allowance, cx))
}

#[test]
fn both_adapters_conform_across_knots_and_keep_the_original_endpoint_law() {
    let capacitor = NonlinearCapacitor::new(CapacitorChargeCurve::try_new(knots()).unwrap());
    let inductor = NonlinearInductor::new(InductorFluxCurve::try_new(knots()).unwrap());
    // Charge is even under time reversal and flux linkage odd.
    let charge = |q: &f64| *q;
    let flux = |linkage: &f64| -*linkage;
    // The largest gap / allowance of the reversals, per law, and their count.
    let (mut worst, mut reversals) = ([0.0_f64; 3], 0);
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
        for reversal in [
            adapter_reverses(&capacitor, &knots(), old, dt, input, &charge),
            adapter_reverses(&inductor, &knots(), old, dt, input, &flux),
        ] {
            let reversal = reversal
                .unwrap_or_else(|failure| panic!("reversal at {old}, {dt} s, {input}: {failure}"));
            // A gap of zero is no part of its allowance, also of an allowance of zero.
            let ratio = |gap: f64, bound: f64| if gap > 0.0 { gap / bound } else { 0.0 };
            let allowance = &reversal.allowance;
            worst[0] = worst[0].max(ratio(reversal.output_gaps[0], allowance.outputs[0]));
            worst[1] = worst[1].max(ratio(reversal.observed_gaps[0], allowance.observed[0]));
            worst[2] = worst[2].max(ratio(reversal.energy_gap_j, allowance.energy_j));
            reversals += 1;
        }
    }
    println!(
        "adapters: {reversals} reversals; largest gap / allowance: outputs {:.3e}, observed {:.3e}, energy {:.3e}",
        worst[0], worst[1], worst[2]
    );
    assert_eq!(reversals, 10);
}

#[test]
fn a_capacitor_reversal_that_negates_the_charge_is_refused_at_the_outputs() {
    // Taking charge for odd sends the backward step from -q1 under -i, so its
    // mean voltage has the wrong sign.
    let capacitor = NonlinearCapacitor::new(CapacitorChargeCurve::try_new(knots()).unwrap());
    let failure = adapter_reverses(&capacitor, &knots(), 0.2, 0.1, 2.0, &|q: &f64| -*q)
        .expect_err("charge is even");
    let Failure::Law { law, detail } = failure else {
        panic!("expected a law failure, got {failure}")
    };
    println!("refused: {law} ({detail})");
    assert_eq!(
        law,
        "the backward outputs are the forward outputs times their parity"
    );
    assert!(detail.starts_with("output 0:"), "{detail}");
}

#[test]
fn an_inductor_whose_curve_is_not_odd_is_refused_at_the_outputs() {
    // Slope 2 below unit flux on both sides, 4 above it only for positive
    // flux. Reversed, the motion follows -I(-flux), which is another law.
    let skewed = vec![
        (-2.0, -4.0),
        (-1.0, -2.0),
        (0.0, 0.0),
        (1.0, 2.0),
        (2.0, 6.0),
    ];
    let inductor = NonlinearInductor::new(InductorFluxCurve::try_new(skewed.clone()).unwrap());
    let flux = |linkage: &f64| -*linkage;
    let failure = adapter_reverses(&inductor, &skewed, 1.2, 0.1, -4.0, &flux)
        .expect_err("the curve is not odd");
    let Failure::Law { law, detail } = failure else {
        panic!("expected a law failure, got {failure}")
    };
    println!("refused: {law} ({detail})");
    assert_eq!(
        law,
        "the backward outputs are the forward outputs times their parity"
    );
    assert!(detail.starts_with("output 0:"), "{detail}");
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
