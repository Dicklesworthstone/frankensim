use super::*;
use crate::thermal::boundary::FiniteThermalWall;

fn limits() -> ThermalOpenLimits3 {
    ThermalOpenLimits3 {
        maximum_cells: 4096,
        maximum_storage_bytes: 4 * 1024 * 1024,
        temperature_range: [-10.0, 10.0],
        balance_tolerance: 2e-11,
    }
}
fn faces(axis: usize, inlet: f64, outlet: f64) -> ThermalOpenBoundaries3 {
    let mut faces = [ThermalOpenFace3::Periodic; 6];
    faces[2 * axis] = ThermalOpenFace3::Reservoir { temperature: inlet };
    faces[2 * axis + 1] = ThermalOpenFace3::Reservoir {
        temperature: outlet,
    };
    ThermalOpenBoundaries3 {
        faces,
        internal_wall: ThermalWall::Adiabatic,
    }
}
fn field(
    dims: [usize; 3],
    u: [f64; 3],
    temperature: f64,
    tau: f64,
    boundaries: ThermalOpenBoundaries3,
) -> ThermalOpenField3 {
    let cells = dims.into_iter().product();
    ThermalOpenField3::try_new::<ThermalOpenError3>(
        dims,
        &vec![false; cells],
        &vec![temperature; cells],
        tau,
        boundaries,
        limits(),
        0,
        |_| Ok(u),
        || Ok(()),
    )
    .unwrap()
}
fn step(field: &mut ThermalOpenField3, heat: &[f64], tau: Option<&[f64]>) -> ThermalOpenBalance3 {
    let trial = field
        .begin_step::<ThermalOpenError3>(heat, tau, || Ok(()))
        .unwrap();
    let balance = field
        .finish_step::<ThermalOpenError3>(heat, trial, || Ok(()))
        .unwrap();
    field.publish(field.completed_steps + 1, balance);
    balance
}
fn near(a: f64, b: f64, tolerance: f64) {
    assert!((a - b).abs() <= tolerance, "{a:.17e} != {b:.17e}");
}
fn bits(field: &ThermalOpenField3) -> Vec<u64> {
    field
        .populations
        .iter()
        .flatten()
        .map(|v| v.to_bits())
        .collect()
}

#[test]
fn danckwerts_link_closes_flux_and_preserves_signed_equilibrium() {
    for temperature in [-1.7, 0.0, 0.2, 1.3] {
        for u in [-0.249, -0.07, -0.001, 0.0, 0.002, 0.08, 0.249] {
            let outgoing = (0.125 - 0.5 * u) * temperature;
            let incoming = reservoir_incoming(outgoing, u, temperature);
            near(incoming, (0.125 + 0.5 * u) * temperature, 2e-16);
            near(incoming - outgoing, u * temperature, 2e-16);
        }
    }
    let outgoing = 0.183;
    assert_eq!(
        reservoir_incoming(outgoing, 0.0, 8.0).to_bits(),
        outgoing.to_bits()
    );
    assert_eq!(
        reservoir_incoming(outgoing, -0.04, -9.0).to_bits(),
        reservoir_incoming(outgoing, -0.04, 9.0).to_bits()
    );
    let incoming = reservoir_incoming(outgoing, -0.04, 2.0);
    near(
        incoming - outgoing,
        -0.04 * (incoming + outgoing) / 0.25,
        5e-17,
    );
    near(
        reservoir_incoming(outgoing, 0.04, 1.7) - outgoing,
        0.04 * 1.7,
        5e-17,
    );
}

#[test]
fn uniform_temperature_and_opposite_reservoir_heat_survive_all_axes_and_rates() {
    let dims = [5, 4, 3];
    let cells = dims.into_iter().product();
    for axis in 0..3 {
        for sign in [-1.0, 1.0] {
            for tau in [0.61, 0.9, 1.8] {
                let mut u = [0.013, -0.017, 0.011];
                u[axis] = sign * 0.043;
                let mut field = field(dims, u, 1.2, tau, faces(axis, 1.2, 1.2));
                for _ in 0..20 {
                    let receipt = step(&mut field, &vec![0.0; cells], None);
                    let area = cells / dims[axis];
                    near(
                        receipt.reservoirs[2 * axis].heat_in,
                        area as f64 * u[axis] * 1.2,
                        3e-14,
                    );
                    near(
                        receipt.reservoirs[2 * axis + 1].heat_in,
                        -(area as f64) * u[axis] * 1.2,
                        3e-14,
                    );
                    near(receipt.defect, 0.0, 3e-14);
                    for cell in 0..cells {
                        near(field.temperature(cell).unwrap(), 1.2, 3e-14);
                    }
                }
            }
        }
    }
}

#[test]
fn backflow_uses_each_faces_own_reservoir_per_cell() {
    let dims = [5, 4, 2];
    let cells = dims.into_iter().product();
    let velocity = |cell: usize| [if cell / dims[0] % 2 == 0 { 0.04 } else { -0.03 }, 0.0, 0.0];
    let mut field = ThermalOpenField3::try_new::<ThermalOpenError3>(
        dims,
        &vec![false; cells],
        &vec![1.0; cells],
        0.8,
        faces(0, 0.3, 1.7),
        limits(),
        0,
        |cell| Ok(velocity(cell)),
        || Ok(()),
    )
    .unwrap();
    let local_tau: Vec<_> = (0..cells)
        .map(|cell| 0.7 + 0.01 * (cell % 7) as f64)
        .collect();
    let receipt = step(&mut field, &vec![0.0; cells], Some(&local_tau));
    near(receipt.reservoirs[0].inflow_heat, 4.0 * 0.04 * 0.3, 2e-16);
    near(receipt.reservoirs[1].inflow_heat, 4.0 * 0.03 * 1.7, 2e-16);
    near(receipt.reservoirs[0].outflow_heat, -4.0 * 0.03, 2e-16);
    near(receipt.reservoirs[1].outflow_heat, -4.0 * 0.04, 2e-16);
    near(receipt.reservoirs[0].inward_volume, 4.0 * 0.04, 1e-16);
    near(receipt.reservoirs[1].inward_volume, 4.0 * 0.03, 1e-16);
    near(
        receipt.stored_change,
        receipt.reservoirs.iter().map(|f| f.heat_in).sum(),
        2e-15,
    );
}

fn heated_fixture() -> ThermalOpenField3 {
    let mut boundaries = faces(0, 0.4, 0.7);
    for face in [2, 3] {
        boundaries.faces[face] =
            ThermalOpenFace3::Wall(ThermalWall::FiniteSolid(FiniteThermalWall {
                temperature: if face == 2 { 1.2 } else { 0.7 },
                heat_capacity: 2.7,
                resistance: 0.4,
                heat_per_step: 0.0002,
            }));
    }
    field([6, 3, 2], [0.025, 0.0, 0.0], 0.4, 0.83, boundaries)
}

#[test]
fn heated_finite_walls_and_open_reservoirs_close_independent_storage() {
    let mut field = heated_fixture();
    let cells = field.populations.len();
    let heat: Vec<_> = (0..cells)
        .map(|cell| 0.001 * (1 + cell % 3) as f64)
        .collect();
    for _ in 0..30 {
        let before: f64 = field.populations.iter().flatten().sum();
        let receipt = step(&mut field, &heat, None);
        let after: f64 = field.populations.iter().flatten().sum();
        let mut solid_storage = 0.0;
        let mut solid_source = 0.0;
        for face in [2, 3] {
            let exchange = receipt.walls[face].unwrap().solid.unwrap();
            solid_storage += exchange.stored_change;
            solid_source += exchange.source;
            near(exchange.defect(), 0.0, 2e-14);
        }
        near(
            after - before + solid_storage,
            heat.iter().sum::<f64>()
                + solid_source
                + receipt.reservoirs.iter().map(|r| r.heat_in).sum::<f64>(),
            5e-14,
        );
        near(receipt.combined_defect, 0.0, 2e-14);
    }
    assert_ne!(
        field.walls[2][0], field.walls[2][1],
        "independent solid cells must see the developing fluid field"
    );
}

#[test]
fn late_refusal_and_cancellation_preserve_populations_wall_state_history_and_time() {
    let mut actual = heated_fixture();
    let mut expected = heated_fixture();
    let heat = vec![0.001; actual.populations.len()];
    step(&mut actual, &heat, None);
    step(&mut expected, &heat, None);
    let before = bits(&actual);
    let before_walls = actual.walls.clone();
    let before_history = format!("{:?}", actual.history);
    let before_receipt = actual.last_balance;
    let trial = actual
        .begin_step::<ThermalOpenError3>(&heat, None, || Ok(()))
        .unwrap();
    actual.velocities[0][0] = 0.25;
    assert!(matches!(
        actual.finish_step::<ThermalOpenError3>(&heat, trial, || Ok(())),
        Err(ThermalOpenError3::Velocity { .. })
    ));
    actual.velocities[0][0] = 0.025;
    assert_eq!(bits(&actual), before);
    assert_eq!(actual.walls, before_walls);
    assert_eq!(format!("{:?}", actual.history), before_history);
    assert_eq!(actual.last_balance, before_receipt);
    assert_eq!(actual.completed_steps, 1);
    let trial = actual
        .begin_step::<ThermalOpenError3>(&heat, None, || Ok(()))
        .unwrap();
    let mut polls = 0;
    let cancelled = actual.finish_step::<ThermalOpenError3>(&heat, trial, || {
        polls += 1;
        if polls == 4 {
            Err(ThermalOpenError3::Cancelled)
        } else {
            Ok(())
        }
    });
    assert!(matches!(cancelled, Err(ThermalOpenError3::Cancelled)));
    assert_eq!(bits(&actual), before);
    assert_eq!(actual.walls, before_walls);
    assert_eq!(format!("{:?}", actual.history), before_history);
    let receipt = step(&mut actual, &heat, None);
    let control = step(&mut expected, &heat, None);
    assert_eq!(bits(&actual), bits(&expected));
    assert_eq!(actual.walls, expected.walls);
    assert_eq!(receipt, control);
}

#[test]
fn obstacles_exchange_heat_on_owned_links_without_thermal_storage() {
    let dims = [6, 4, 3];
    let cells = dims.into_iter().product();
    let mut solid = vec![false; cells];
    solid[2 + 6 * (1 + 4)] = true;
    let mut boundaries = faces(0, 0.4, 0.9);
    boundaries.internal_wall = ThermalWall::HeatFlux(0.001);
    let mut field = ThermalOpenField3::try_new::<ThermalOpenError3>(
        dims,
        &solid,
        &vec![0.4; cells],
        0.9,
        boundaries,
        limits(),
        0,
        |_| Ok([0.0; 3]),
        || Ok(()),
    )
    .unwrap();
    let heat = vec![0.0; cells];
    let receipt = step(&mut field, &heat, None);
    near(receipt.internal_wall.unwrap().heat_in, 0.006, 1e-16);
    near(receipt.stored_change, 0.006, 1e-14);
    assert_eq!(field.temperature(2 + 6 * 5), None);
    let mut forbidden = heat;
    forbidden[2 + 6 * 5] = 1.0;
    assert!(
        field
            .begin_step::<ThermalOpenError3>(&forbidden, None, || Ok(()))
            .is_err()
    );
}

fn channel_error(n: usize) -> (f64, f64) {
    // Fixed physical Peclet number and total source rise on refinement.
    // Fluid nodes lie at half-link coordinates x=i+1/2 in [0,n].
    let alpha = 0.1;
    let u = 0.4 / n as f64;
    let source = 0.2 * u / n as f64;
    let inlet = 0.7;
    let mut field = field(
        [n, 1, 1],
        [u, 0.0, 0.0],
        inlet,
        0.5 + alpha / D3Q7.cs2,
        faces(0, inlet, 8.0),
    );
    let heat = vec![source; n];
    for _ in 0..(30 * n * n) {
        step(&mut field, &heat, None);
    }
    let mut squared = 0.0;
    for cell in 0..n {
        let x = cell as f64 + 0.5;
        let exact = inlet
            + source * x / u
            + source * alpha / (u * u) * (1.0 - (u * (x - n as f64) / alpha).exp());
        squared += (field.temperature(cell).unwrap() - exact).powi(2);
    }
    let receipt = field.last_balance.unwrap();
    let outflow_temperature = -receipt.reservoirs[1].heat_in / u;
    (
        (squared / n as f64).sqrt(),
        (outflow_temperature - (inlet + 0.2)).abs(),
    )
}

#[test]
fn heated_danckwerts_channel_converges_to_independent_continuum_solution() {
    let (coarse, coarse_outlet) = channel_error(12);
    let (fine, fine_outlet) = channel_error(24);
    assert!(
        fine < 0.6 * coarse,
        "RMS error did not refine: {coarse:e} -> {fine:e}"
    );
    assert!(fine < 0.004, "fine RMS temperature error {fine:e}");
    // The steady outlet enthalpy follows directly from the integrated source,
    // independently of the finite-grid interior spatial error.
    assert!(
        coarse_outlet < 2e-8 && fine_outlet < 2e-8,
        "outlet errors {coarse_outlet:e}, {fine_outlet:e}"
    );
}

mod inputs;

mod profiles;
