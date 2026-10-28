//! Physical geometry regressions over the same source-owned Euler solver.
use super::*;

fn cone_mesh(n: usize) -> Mesh1d {
    let positions: Vec<_> = (0..=n).map(|i| (i as f64 / n as f64).powf(1.3)).collect();
    let areas = positions.iter().map(|x| (0.5 + x).powi(2)).collect();
    let volumes = positions
        .windows(2)
        .map(|x| ((0.5 + x[1]).powi(3) - (0.5 + x[0]).powi(3)) / 3.0)
        .collect();
    Mesh1d::new(positions, areas, volumes).unwrap()
}
fn closed() -> Boundaries {
    Boundaries::Ends {
        left: Boundary::Reflecting,
        right: Boundary::Reflecting,
    }
}

#[test]
fn nonuniform_taper_preserves_rest_and_accounts_for_wall_pressure_impulse() {
    let gas = gas(0.4, 0.0);
    let mesh = cone_mesh(53);
    assert!(mesh.uniform_cell_width().is_none());
    assert!(mesh.uniform_cross_section().is_none());
    let states: Vec<_> = (0..mesh.cell_count())
        .map(|i| {
            // A stationary contact field, not just identical initial cells.
            gas.conserved(0.8 + 0.2 * (i as f64 / 53.0), 0.0, 1.2)
                .unwrap()
        })
        .collect();
    let original = states.clone();
    let mut flow = Euler1d::new_on_mesh(gas, mesh, states, closed(), || false).unwrap();
    for _ in 0..24 {
        let r = flow
            .step(Time::new(0.01), StepControl::default(), || false)
            .unwrap();
        let expected_impulse = 1.2 * (2.25 - 0.25) * r.elapsed.value();
        assert!((r.area_momentum_impulse - expected_impulse).abs() < 2e-14);
        assert!((r.boundary_change.momentum + r.area_momentum_impulse).abs() < 2e-14);
        assert!(r.conservation_residual.momentum.abs() < 2e-14);
    }
    for (&a, &b) in flow.states().iter().zip(&original) {
        assert!((a.density - b.density).abs() < 3e-13);
        assert!(a.momentum.abs() < 3e-13);
        assert!((a.energy - b.energy).abs() < 3e-12);
    }
}

#[test]
fn finite_amplitude_tapered_pulse_closes_mass_energy_and_geometric_momentum() {
    let gas = gas(0.35, 1e18);
    let mesh = cone_mesh(100);
    let states = mesh
        .positions()
        .windows(2)
        .map(|faces| {
            let x = 0.5 * (faces[0] + faces[1]);
            let p = 1.0 + 0.8 * (-((x - 0.4) / 0.08).powi(2)).exp();
            gas.conserved(1.0, 0.0, p).unwrap()
        })
        .collect();
    let mut flow = Euler1d::new_on_mesh(gas, mesh, states, closed(), || false).unwrap();
    let original = flow.totals();
    let mut impulse = 0.0;
    let mut steps = 0;
    while flow.time().value() < 0.12 {
        let r = flow
            .step(
                Time::new(0.12 - flow.time().value()),
                StepControl::default(),
                || false,
            )
            .unwrap();
        impulse += r.boundary_change.momentum + r.area_momentum_impulse;
        assert!(r.maximum_cfl <= 0.40000000000001);
        assert!(r.conservation_residual.mass.abs() < 5e-14);
        assert!(r.conservation_residual.momentum.abs() < 5e-14);
        assert!(r.conservation_residual.energy.abs() < 5e-13);
        steps += 1;
        assert!(steps < 5000);
    }
    let final_state = flow.totals();
    assert!((final_state.mass - original.mass).abs() < 2e-13);
    assert!((final_state.energy - original.energy).abs() < 2e-12);
    assert!((final_state.momentum - original.momentum - impulse).abs() < 2e-13);
    assert!(flow.states().iter().any(|q| q.momentum.abs() > 0.01));
}

#[test]
fn stretched_periodic_contact_advection_retains_second_order_accuracy() {
    fn error(n: usize) -> f64 {
        let gas = gas(0.0, 0.0);
        let wave = 2.0 * core::f64::consts::PI;
        let positions: Vec<_> = (0..=n)
            .map(|i| {
                let s = i as f64 / n as f64;
                s + 0.25 / wave * (wave * s).sin()
            })
            .collect();
        let volumes = positions.windows(2).map(|x| x[1] - x[0]).collect();
        let mesh = Mesh1d::new(positions, vec![1.0; n + 1], volumes).unwrap();
        let average_rho = |l: f64, r: f64, t: f64| {
            1.0 + 0.1 * ((wave * (l - 0.3 * t)).cos() - (wave * (r - 0.3 * t)).cos())
                / (wave * (r - l))
        };
        let states = mesh
            .positions()
            .windows(2)
            .map(|x| {
                gas.conserved(average_rho(x[0], x[1], 0.0), 0.3, 1.0)
                    .unwrap()
            })
            .collect();
        let mut flow =
            Euler1d::new_on_mesh(gas, mesh, states, Boundaries::Periodic, || false).unwrap();
        while flow.time().value() < 0.2 {
            flow.step(
                Time::new(0.2 - flow.time().value()),
                StepControl::default(),
                || false,
            )
            .unwrap();
        }
        flow.states()
            .iter()
            .zip(flow.mesh().positions().windows(2))
            .map(|(q, x)| (q.density - average_rho(x[0], x[1], 0.2)).abs() * (x[1] - x[0]))
            .sum()
    }
    let coarse = error(40);
    let fine = error(80);
    eprintln!("nonuniform contact density L1 40/80 cells: {coarse:.9e}/{fine:.9e}");
    assert!(
        fine < 0.35 * coarse,
        "loss of second order: {coarse} -> {fine}"
    );
}

#[test]
fn geometry_refuses_invalid_faces_and_mismatched_periodic_areas_transactionally() {
    assert!(Mesh1d::new(vec![0.0, 0.5, 0.5], vec![1.0; 3], vec![0.5; 2]).is_err());
    assert!(Mesh1d::new(vec![0.0, 0.5, 1.0], vec![1.0; 3], vec![0.0, 0.5]).is_err());
    let gas = gas(0.0, 0.0);
    let q = gas.conserved(1.0, 0.0, 1.0).unwrap();
    assert!(
        Euler1d::new_on_mesh(
            gas.clone(),
            cone_mesh(2),
            vec![q; 2],
            Boundaries::Periodic,
            || false
        )
        .is_err()
    );
    let mut flow = Euler1d::new_on_mesh(gas, cone_mesh(8), vec![q; 8], closed(), || false).unwrap();
    let before = flow.states().to_vec();
    let mut polls = 0;
    assert_eq!(
        flow.step(Time::new(0.01), StepControl::default(), || {
            polls += 1;
            polls > 25
        }),
        Err(Error::Cancelled)
    );
    assert_eq!(flow.states(), before);
    assert_eq!(flow.time().value(), 0.0);
}
