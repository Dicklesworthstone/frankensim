//! Public moving-aperture consumer tests; all material values are synthetic.
use fs_couple::bernoulli_aperture::BernoulliAperture;
use fs_couple::bernoulli_aperture::dynamic::{
    ApertureDrive, ApertureFrame, ApertureState, ApertureTerminal,
    DynamicAperture, DynamicApertureSpec,
};
use fs_couple::bernoulli_aperture::paired::BladePropagator;
use fs_dcontact::Obstacle;
use fs_exec::CancelGate;

fn spec(max_steps: u64) -> DynamicApertureSpec {
    DynamicApertureSpec {
        aperture: BernoulliAperture {
            rest_opening_m: 4e-4,
            width_m: 0.013,
            closing_pressure_pa: 6000.0,
        },
        mass_kg: 1e-5,
        stiffness_n_m: 500.0,
        damping_ratio: 0.35,
        density_kg_m3: 1.2,
        impedance_pa_s_m3: 1e6,
        time_step_s: 1e-6,
        max_steps,
    }
}

fn lay(stiffness: f64) -> Obstacle {
    Obstacle::new(vec![-1.0], 1, 1, vec![0.0], vec![1.0], stiffness, 2.0,
        "synthetic caller-owned contact fixture".into()).unwrap()
}

fn model(max_steps: u64) -> DynamicAperture {
    DynamicAperture::new(spec(max_steps), ApertureState {
        opening_m: 4e-4, opening_velocity_m_s: 0.0,
    }, lay(1e8)).unwrap()
}

fn bits(f: ApertureFrame) -> Vec<u64> {
    let mut values = vec![f.step];
    values.extend([
        f.time_s, f.state.opening_m, f.state.opening_velocity_m_s,
        f.midpoint_opening_m, f.outgoing_pressure_pa, f.bore_pressure_pa,
        f.jet_flow_m3_s, f.swept_flow_m3_s, f.bore_flow_m3_s,
        f.flow_residual_m3_s, f.stored_energy_j, f.storage_change_j,
        f.dissipated_energy_j, f.pressure_work_j,
    ].map(f64::to_bits));
    values
}

#[test]
fn public_step_uses_the_supplied_contact_law_and_preserves_its_identity() {
    let initial = ApertureState { opening_m: -1e-4, opening_velocity_m_s: -0.2 };
    // The pole-rate bound is about 12100 rad/s: midpoint is admitted only
    // below 3.3e-6 s, while the exact block takes the original 1e-5 s clock.
    for (propagator, dt) in [(BladePropagator::Midpoint, 3e-6), (BladePropagator::ExactZoh, 1e-5)] {
        let mut s = spec(1);
        s.time_step_s = dt;
        let mut soft = DynamicAperture::new_with_propagator(s, initial, lay(1e6), propagator).unwrap();
        let mut stiff = DynamicAperture::new_with_propagator(s, initial, lay(1e9), propagator).unwrap();
        assert_eq!(soft.contact_law().stiffness(), 1e6);
        assert_eq!(stiff.contact_law().stiffness(), 1e9);
        assert_eq!(soft.contact_law().provenance(), "synthetic caller-owned contact fixture");
        let a = soft.step(ApertureDrive::default()).unwrap();
        let b = stiff.step(ApertureDrive::default()).unwrap();
        assert!((a.state.opening_m - b.state.opening_m).abs() > 1e-9, "{propagator:?}");
        assert!((a.outgoing_pressure_pa - b.outgoing_pressure_pa).abs() > 1.0, "{propagator:?}");
        for f in [a, b] {
            let scale = f.stored_energy_j.abs() + f.storage_change_j.abs()
                + f.dissipated_energy_j.abs() + f.pressure_work_j.abs();
            assert!(f.balance_residual_j().abs() <= 2e-10 * scale, "{propagator:?}");
            assert!(f.dissipated_energy_j >= 0.0, "{propagator:?}");
        }
    }
    // The original clock is outside the midpoint band and refuses there.
    let mut s = spec(1);
    s.time_step_s = 1e-5;
    assert!(DynamicAperture::new_with_propagator(s, initial, lay(1e6), BladePropagator::Midpoint)
        .and_then(|mut m| m.step(ApertureDrive::default())).is_err());
}

/// Checks one free-opening step against an independently manufactured frame:
/// the caller supplies the end velocity and the pressure-face force that the
/// propagator must have held to reach `next_y`. The jet uses the
/// endpoint-average opening and the swept flow uses the step-mean velocity in
/// both propagators.
fn check_manufactured_frame(propagator: BladePropagator, dt: f64, (y, next_y, v): (f64, f64, f64),
    next_v: f64, held_force_n: f64) {
    let mut s = spec(1);
    s.time_step_s = dt;
    let vm = (next_y - y) / dt;
    let face = s.stiffness_n_m * s.aperture.rest_opening_m / s.aperture.closing_pressure_pa;
    let dp = -held_force_n / face;
    let jet = s.aperture.width_m * (0.5 * (y + next_y)) * dp.signum()
        * (2.0 * dp.abs() / s.density_kg_m3).sqrt();
    let incoming = 75.0;
    let body = 2e-7;
    let outgoing = incoming + s.impedance_pa_s_m3 * (jet + body - face * vm);
    let mut model = DynamicAperture::new_with_propagator(s, ApertureState {
        opening_m: y, opening_velocity_m_s: v,
    }, lay(1e8), propagator).unwrap();
    let f = model.step(ApertureDrive {
        upstream_pressure_pa: dp + outgoing + incoming,
        incoming_pressure_pa: incoming,
        body_flow_m3_s: body,
    }).unwrap();
    assert!((f.state.opening_m - next_y).abs() < 1e-13, "{propagator:?}");
    assert!((f.state.opening_velocity_m_s - next_v).abs() < 1e-9, "{propagator:?}");
    assert!((f.outgoing_pressure_pa - outgoing).abs() < 1e-8, "{propagator:?}");
    assert!((f.midpoint_opening_m - 0.5 * (y + next_y)).abs() < 1e-13, "{propagator:?}");
    assert!((f.jet_flow_m3_s - jet).abs() < 1e-14, "{propagator:?}");
    assert!(f.flow_residual_m3_s.abs() < 1e-14, "{propagator:?}");
    assert!(f.balance_residual_j().abs() < 1e-14, "{propagator:?}");
    assert_eq!(f.step, 1);
    assert_eq!(f.time_s.to_bits(), s.time_step_s.to_bits());
}

#[test]
fn physical_frame_matches_an_independent_manufactured_moving_aperture() {
    let s = spec(1);
    let (m, k) = (s.mass_kg, s.stiffness_n_m);
    let rest = s.aperture.rest_opening_m;
    let damping = 2.0 * s.damping_ratio * (k * m).sqrt();

    // Midpoint, inside its band: next_v = 2 vm - v and the midpoint momentum
    // balance m (next_v - v)/dt + k (y_mid - rest) + c vm = F.
    let (dt, y, v) = (2e-6, 4e-4, 0.2);
    let next_y = y + 0.19 * dt;
    let vm = (next_y - y) / dt;
    let next_v = 2.0 * vm - v;
    let force = m * (next_v - v) / dt + k * (0.5 * (y + next_y) - rest) + damping * vm;
    check_manufactured_frame(BladePropagator::Midpoint, dt, (y, next_y, v), next_v, force);

    // Exact, at omega*dt = 0.71.
    let (dt, y, next_y, v) = (1e-4, 4e-4, 3e-4, 0.2);
    let (force, next_v) = exact_held_force(&s, dt, (y, next_y, v));
    // At this clock the midpoint image (2 vm - v = -2.2 m/s) misses the exact
    // end velocity (about -1.94 m/s) by far more than the frame tolerance.
    assert!((next_v - (2.0 * (next_y - y) / dt - v)).abs() > 0.1);
    check_manufactured_frame(BladePropagator::ExactZoh, dt, (y, next_y, v), next_v, force);
}

/// The force an exact step must hold to move the opening from `y` to `next_y`,
/// and the end velocity it reaches. Closed-form underdamped response about the
/// shifted equilibrium u = x - rest - F/k:
/// u(t) = exp(-a t)[u0 (cos + a/b sin) + v0/b sin], with a = zeta*omega and
/// b = omega*sqrt(1 - zeta^2).
fn exact_held_force(s: &DynamicApertureSpec, dt: f64, (y, next_y, v): (f64, f64, f64)) -> (f64, f64) {
    let (k, rest) = (s.stiffness_n_m, s.aperture.rest_opening_m);
    let omega = (k / s.mass_kg).sqrt();
    let a = s.damping_ratio * omega;
    let b = omega * (1.0 - s.damping_ratio * s.damping_ratio).sqrt();
    let (e, c, sn) = ((-a * dt).exp(), (b * dt).cos(), (b * dt).sin());
    let free = |u0: f64, v0: f64| e * (u0 * (c + a / b * sn) + v0 * sn / b);
    let free_v = |u0: f64, v0: f64| e * (v0 * (c - a / b * sn) - u0 * omega * omega / b * sn);
    let force = k * (next_y - rest - free(y - rest, v)) / (1.0 - e * (c + a / b * sn));
    (force, free_v(y - rest - force / k, v))
}

#[test]
fn cancellation_and_budget_resume_preserve_the_exact_accepted_trajectory() {
    let inputs: Vec<_> = (0..32).map(|i| ApertureDrive {
        upstream_pressure_pa: 500.0 + f64::from(i) * 10.0,
        incoming_pressure_pa: 20.0,
        body_flow_m3_s: 2e-7,
    }).collect();
    let gate = CancelGate::new_clock_free();
    let mut uninterrupted = model(32);
    let mut expected = vec![ApertureFrame::default(); 32];
    assert_eq!(uninterrupted.advance_block(&inputs, &mut expected, &gate).unwrap().terminal,
        ApertureTerminal::Complete);
    let sentinel = ApertureFrame { step: u64::MAX, ..ApertureFrame::default() };
    let mut actual = vec![sentinel; 32];
    let mut resumed = model(12);
    resumed.advance_block(&inputs[..8], &mut actual[..8], &gate).unwrap();
    let cancel = CancelGate::new_clock_free();
    cancel.request();
    let paused = resumed.advance_block(&inputs[8..], &mut actual[8..], &cancel).unwrap();
    assert_eq!(paused.completed, 0);
    assert_eq!(paused.terminal, ApertureTerminal::Cancelled);
    assert_eq!(resumed.accepted_steps(), 8);
    assert!(actual[8..].iter().all(|f| *f == sentinel));
    let exhausted = resumed.advance_block(&inputs[8..], &mut actual[8..], &gate).unwrap();
    assert_eq!(exhausted.completed, 4);
    assert_eq!(exhausted.terminal, ApertureTerminal::BudgetExhausted);
    assert!(actual[12..].iter().all(|f| *f == sentinel));
    assert_eq!(resumed.state(), expected[11].state);
    resumed.extend_step_budget(32).unwrap();
    resumed.advance_block(&inputs[12..], &mut actual[12..], &gate).unwrap();
    assert_eq!(resumed.accepted_steps(), 32);
    for (a, b) in actual.into_iter().zip(expected) { assert_eq!(bits(a), bits(b)); }
}

#[test]
fn refusal_is_transactional_and_does_not_poison_a_valid_next_step() {
    let mut a = model(4);
    let mut b = model(4);
    let valid = ApertureDrive { upstream_pressure_pa: 500.0, ..ApertureDrive::default() };
    a.step(valid).unwrap();
    b.step(valid).unwrap();
    let before = a.state();
    for bad in [f64::NAN, f64::INFINITY, f64::MAX] {
        assert!(a.step(ApertureDrive { body_flow_m3_s: bad, ..valid }).is_err());
        assert_eq!(a.state(), before);
        assert_eq!(a.accepted_steps(), 1);
    }
    let mut empty = [];
    assert!(a.advance_block(&[valid], &mut empty, &CancelGate::new_clock_free()).is_err());
    assert_eq!(bits(a.step(valid).unwrap()), bits(b.step(valid).unwrap()));
    assert!(a.extend_step_budget(3).is_err());
    assert_eq!(a.spec().max_steps, 4);
}

#[test]
fn explicit_admission_rejects_hidden_defaults_and_malformed_contact() {
    let state = ApertureState { opening_m: 4e-4, opening_velocity_m_s: 0.0 };
    for field in 0..4 {
        let mut s = spec(4);
        match field {
            0 => s.stiffness_n_m = 0.0,
            1 => s.mass_kg = 0.0,
            2 => s.damping_ratio = f64::NAN,
            _ => s.max_steps = 0,
        }
        assert!(DynamicAperture::new(s, state, lay(1e8)).is_err());
    }
    let raw = Obstacle::from_raw_parts(vec![-1.0], 1, vec![], vec![], 1e8, 2.0,
        "malformed raw escape".into());
    assert!(DynamicAperture::new(spec(4), state, raw).is_err());
    let nonunit = Obstacle::new(vec![-2.0], 1, 1, vec![0.0], vec![1.0], 1e8, 2.0,
        "different coordinate".into()).unwrap();
    let admitted = DynamicAperture::new(spec(4), state, nonunit).unwrap();
    assert_eq!(admitted.contact_law().collocation(), [-2.0]);
    let too_many = Obstacle::new(vec![-1.0;4097],4097,1,vec![0.0;4097],vec![1.0;4097],
        1e8,2.0,"over-budget quadrature".into()).unwrap();
    assert!(DynamicAperture::new(spec(4),state,too_many).is_err());
}

#[test]
fn accepted_outgoing_waves_drive_a_causal_delayed_characteristic_consumer() {
    // A prescribed passive reflection and propagation delay, not a microphone
    // or a fitted physical instrument. Feed actual accepted output back in.
    let render = |delay: usize| {
        let mut model = model(128);
        let mut line = vec![0.0; delay];
        let mut pressure = Vec::new();
        for i in 0..128 {
            let incoming = -0.8 * line[i % delay];
            let frame = model.step(ApertureDrive {
                upstream_pressure_pa: 800.0,
                incoming_pressure_pa: incoming,
                body_flow_m3_s: 0.0,
            }).unwrap();
            line[i % delay] = frame.outgoing_pressure_pa;
            pressure.push(frame.bore_pressure_pa);
        }
        pressure
    };
    let short = render(16);
    let long = render(24);
    assert!(short.iter().any(|p| p.abs() > 1.0));
    for i in 0..16 { assert_eq!(short[i].to_bits(), long[i].to_bits()); }
    assert!(short[16..].iter().zip(&long[16..]).any(|(a, b)| (a-b).abs() > 1.0));
    assert_eq!(short.iter().map(|p| p.to_bits()).collect::<Vec<_>>(),
        render(16).iter().map(|p| p.to_bits()).collect::<Vec<_>>());
}

#[path = "dynamic_aperture/distributed.rs"]
mod distributed;
