use super::*;
use crate::dynamics::tests::gas;
use crate::dynamics::{Boundaries, Boundary, Euler1d, Mesh1d, StepControl};
use fs_qty::Time;

fn closed() -> Boundaries {
    Boundaries::Ends {
        left: Boundary::Reflecting,
        right: Boundary::Reflecting,
    }
}
fn open() -> Boundaries {
    Boundaries::Ends {
        left: Boundary::Transmissive,
        right: Boundary::Transmissive,
    }
}

// `n_left` cells of width `h` seeing area `a1`, then `n_right` cells seeing `a2`; the step is face `n_left`.
fn stepped(n_left: usize, n_right: usize, a1: f64, a2: f64, h: f64) -> (Mesh1d, Vec<(usize, f64)>) {
    let n = n_left + n_right;
    let positions = (0..=n).map(|i| h * i as f64).collect();
    let areas = (0..=n).map(|i| if i <= n_left { a1 } else { a2 }).collect();
    let volumes = (0..n)
        .map(|i| h * if i < n_left { a1 } else { a2 })
        .collect();
    (
        Mesh1d::new(positions, areas, volumes).unwrap(),
        vec![(n_left, a2)],
    )
}

#[test]
fn a_uniform_resting_gas_stays_at_rest_across_a_step_and_the_annulus_pushes_with_the_pressure() {
    let gas = gas(0.4, 0.0);
    let (mesh, steps) = stepped(30, 50, 1.0, 1.7, 0.02);
    let states: Vec<_> = (0..80)
        .map(|_| gas.conserved(0.8, 0.0, 1.2).unwrap())
        .collect();
    let original = states.clone();
    let mut flow =
        Euler1d::new_on_stepped_mesh(gas, mesh, &steps, states, closed(), || false).unwrap();
    assert_eq!(flow.radius_steps(), vec![(30, 1.7)]);
    for _ in 0..40 {
        let r = flow
            .step(Time::new(0.01), StepControl::default(), || false)
            .unwrap();
        // The closed tube's wall force is the pressure times the area change, the step annulus included.
        let expected = 1.2 * (1.7 - 1.0) * r.elapsed.value();
        assert!((r.area_momentum_impulse - expected).abs() < 2e-14);
        assert!(r.conservation_residual.mass.abs() < 2e-14);
        assert!(r.conservation_residual.momentum.abs() < 2e-14);
        assert!(r.conservation_residual.energy.abs() < 2e-13);
    }
    for (a, b) in flow.states().iter().zip(&original) {
        assert!((a.density - b.density).abs() < 3e-13);
        assert!(a.momentum.abs() < 3e-13);
        assert!((a.energy - b.energy).abs() < 3e-12);
    }
}

// A right-going half-cosine pressure pulse (a simple wave) centred on cell `centre`, `half_width` cells to a side.
fn right_going_pulse(
    gas: &crate::dynamics::FrozenGas,
    n: usize,
    centre: usize,
    half_width: usize,
    relative_amplitude: f64,
) -> Vec<crate::dynamics::Conserved> {
    let (rho0, p0) = (0.8, 1.2);
    let rest = gas.conserved(rho0, 0.0, p0).unwrap();
    let c0 = gas.primitive(rest, || false).unwrap().sound_speed;
    (0..n)
        .map(|i| {
            let x = (i as f64 - centre as f64) / half_width as f64;
            let dp = if x.abs() < 1.0 {
                relative_amplitude * p0 * 0.5 * (1.0 + (core::f64::consts::PI * x).cos())
            } else {
                0.0
            };
            gas.conserved(rho0 + dp / (c0 * c0), dp / (rho0 * c0), p0 + dp)
                .unwrap()
        })
        .collect()
}

#[test]
fn a_finite_pulse_through_a_step_closes_mass_and_energy_to_rounding_and_the_annulus_momentum() {
    // Bounds fixed before the first run (T6): a pulse of 5 percent of the pressure through a 1 -> 1.7 area step, open ends
    // until it has crossed: mass and energy residuals of the extensive balance within 1e-12 of the totals at every step,
    // the momentum residual (the annulus impulse is part of it) within 1e-9 of the pressure force scale, and every cell finite and positive.
    let gas = gas(0.4, 0.0);
    let (mesh, steps) = stepped(120, 180, 1.0, 1.7, 1.0 / 120.0);
    let states = right_going_pulse(&gas, 300, 50, 24, 0.05);
    let mut flow =
        Euler1d::new_on_stepped_mesh(gas, mesh, &steps, states, open(), || false).unwrap();
    let mut worst = [0.0_f64; 3];
    let mut count = 0;
    while flow.time().value() < 0.85 {
        let before = flow.totals();
        let r = flow
            .step(
                Time::new(0.85 - flow.time().value()),
                StepControl::default(),
                || false,
            )
            .unwrap();
        let scale = [
            before.mass.abs(),
            before.energy.abs(),
            1.2 * 1.7 * r.elapsed.value(),
        ];
        let residual = r.conservation_residual;
        worst[0] = worst[0].max(residual.mass.abs() / scale[0]);
        worst[1] = worst[1].max(residual.energy.abs() / scale[1]);
        worst[2] = worst[2].max(residual.momentum.abs() / scale[2]);
        count += 1;
        assert!(count < 10_000);
    }
    println!(
        "T6 worst relative residuals: mass {:e}, energy {:e}, momentum {:e}",
        worst[0], worst[1], worst[2]
    );
    assert!(
        worst[0] <= 1e-12 && worst[1] <= 1e-12 && worst[2] <= 1e-9,
        "{worst:?}"
    );
    assert!(
        flow.states()
            .iter()
            .all(|q| q.density.is_finite() && q.density > 0.0 && q.energy.is_finite())
    );
}

#[test]
fn a_small_pulse_reflects_and_transmits_with_the_plane_wave_coefficients_of_the_step() {
    // Bound fixed before the first run: the spatial integral of the pressure disturbance sent back, and the one sent on, over the
    // incident one, are (A1 - A2)/(A1 + A2) and 2 A1/(A1 + A2) within 3e-3 (the gas is the same on both sides, so the wave speed, and with
    // it the ratio of spatial to temporal integrals, is the same).
    let gas = gas(0.4, 0.0);
    let (a1, a2) = (1.0, 1.7);
    let (mesh, steps) = stepped(150, 250, a1, a2, 1.0 / 150.0);
    let states = right_going_pulse(&gas, 400, 50, 24, 1e-5);
    let rest = gas.conserved(0.8, 0.0, 1.2).unwrap();
    let c0 = gas.primitive(rest, || false).unwrap().sound_speed;
    let p0 = 1.2;
    let mut flow =
        Euler1d::new_on_stepped_mesh(gas.clone(), mesh, &steps, states, open(), || false).unwrap();
    let disturbance = |flow: &Euler1d, from: usize, to: usize| -> f64 {
        flow.states()[from..to]
            .iter()
            .map(|q| gas.primitive(*q, || false).unwrap().pressure - p0)
            .sum::<f64>()
    };
    let sent = disturbance(&flow, 0, 150);
    // The pulse reaches the step after (150 - 50) cells / (c0 h); run until it has left it on both sides.
    let end = 2.0 * 0.35 / c0 + 0.25;
    while flow.time().value() < end {
        flow.step(
            Time::new(end - flow.time().value()),
            StepControl::default(),
            || false,
        )
        .unwrap();
    }
    let returned = disturbance(&flow, 0, 150) / sent;
    let passed = disturbance(&flow, 150, 400) / sent;
    let gamma = (a1 - a2) / (a1 + a2);
    let tau = 2.0 * a1 / (a1 + a2);
    println!(
        "small pulse: reflected {returned:.6} (expected {gamma:.6}), transmitted {passed:.6} (expected {tau:.6})"
    );
    assert!(
        (returned - gamma).abs() <= 3e-3,
        "{returned} against {gamma}"
    );
    assert!((passed - tau).abs() <= 3e-3, "{passed} against {tau}");
}

#[test]
fn without_steps_the_stepped_constructor_is_the_plain_one_bit_for_bit_and_bad_step_lists_refuse() {
    let gas = gas(0.4, 0.0);
    let mk = || {
        let positions: Vec<_> = (0..=60).map(|i| (i as f64 / 60.0).powf(1.2)).collect();
        let areas: Vec<_> = positions.iter().map(|x| (0.5 + x).powi(2)).collect();
        let volumes = positions
            .windows(2)
            .map(|x| ((0.5 + x[1]).powi(3) - (0.5 + x[0]).powi(3)) / 3.0)
            .collect();
        Mesh1d::new(positions, areas, volumes).unwrap()
    };
    let states: Vec<_> = (0..60)
        .map(|i| {
            gas.conserved(
                0.8 + 0.01 * (i % 7) as f64,
                0.0,
                1.2 + 0.05 * (i % 5) as f64,
            )
            .unwrap()
        })
        .collect();
    let mut plain =
        Euler1d::new_on_mesh(gas.clone(), mk(), states.clone(), closed(), || false).unwrap();
    let mut empty =
        Euler1d::new_on_stepped_mesh(gas.clone(), mk(), &[], states.clone(), closed(), || false)
            .unwrap();
    for _ in 0..100 {
        let a = plain
            .step(Time::new(0.01), StepControl::default(), || false)
            .unwrap();
        let b = empty
            .step(Time::new(0.01), StepControl::default(), || false)
            .unwrap();
        assert_eq!(a.elapsed.value().to_bits(), b.elapsed.value().to_bits());
        assert_eq!(
            a.area_momentum_impulse.to_bits(),
            b.area_momentum_impulse.to_bits()
        );
    }
    for (a, b) in plain.states().iter().zip(empty.states()) {
        assert_eq!(a.density.to_bits(), b.density.to_bits());
        assert_eq!(a.momentum.to_bits(), b.momentum.to_bits());
        assert_eq!(a.energy.to_bits(), b.energy.to_bits());
    }
    // A step needs an interior face, increasing faces, and a positive area that differs from the area on its left.
    let same_area = mk().face_areas()[10];
    for bad in [
        vec![(0, 2.0)],
        vec![(60, 2.0)],
        vec![(10, 2.0), (10, 3.0)],
        vec![(20, 2.0), (10, 3.0)],
        vec![(10, f64::NAN)],
        vec![(10, 0.0)],
        vec![(10, -1.0)],
        vec![(10, same_area)],
    ] {
        let refused =
            Euler1d::new_on_stepped_mesh(gas.clone(), mk(), &bad, states.clone(), closed(), || {
                false
            });
        assert!(refused.is_err(), "{bad:?}");
    }
}

// Two rest states, the left at 1.2 and the right at `ratio` times that, across a contraction to 0.3 of the area.
fn contraction_flow(ratio: f64) -> (Euler1d, Vec<crate::dynamics::Conserved>) {
    let gas = gas(0.4, 0.0);
    let (mesh, steps) = stepped(30, 30, 1.0, 0.3, 0.02);
    let states: Vec<_> = (0..60)
        .map(|i| {
            gas.conserved(0.8, 0.0, if i < 30 { 1.2 } else { 1.2 * ratio })
                .unwrap()
        })
        .collect();
    let original = states.clone();
    (
        Euler1d::new_on_stepped_mesh(gas, mesh, &steps, states, open(), || false).unwrap(),
        original,
    )
}

#[test]
fn strong_pressure_drops_across_a_contraction_are_solved_and_close_mass_and_energy() {
    // Two earlier expectations of this test were wrong and are scored as misses in
    // brasstap/predictions-step-junction.txt: that a tenfold and then a twentyfold drop across a contraction
    // to 0.3 of the area leaves the junction no root and is refused. A subsonic root exists for every drop
    // (the throat Mach number only approaches 1 as the drop grows; a tenfold drop needs over a hundred chord
    // iterations in the constant-gamma model), so the steps are accepted. What is checked is what holds:
    // mass and energy close to rounding, densities stay finite and positive, and the clock advances.
    // Expectations fixed before the first run of this version.
    for ratio in [0.3, 0.1, 0.05] {
        let (mut flow, _) = contraction_flow(ratio);
        let before = flow.totals();
        let report = flow
            .step(Time::new(0.01), StepControl::default(), || false)
            .unwrap();
        assert!(report.conservation_residual.mass.abs() <= 1e-12 * before.mass.abs());
        assert!(report.conservation_residual.energy.abs() <= 1e-12 * before.energy.abs());
        assert!(
            flow.states()
                .iter()
                .all(|q| q.density.is_finite() && q.density > 0.0)
        );
        assert!(flow.time().value() > 0.0);
    }
}
