//! G1 rectangular Poiseuille flow; G3 shape and frame invariance; G4 rollback.
use super::*;
use crate::multicomponent::network::wall_exchange::{DevelopedWallModel, RectangularLaminarWall};

fn rectangle(width: f64, height: f64, nu: f64) -> RectangularLaminarWall {
    RectangularLaminarWall::try_new(
        Length::new(width),
        Length::new(height),
        nu,
        [0.0, 2300.0],
        0.3,
    )
    .unwrap()
}

fn passage(
    g: &MixtureGas,
    model: RectangularLaminarWall,
    density: f64,
    speed: f64,
) -> MixtureNetwork1d {
    let flow = MixtureEuler1d::new(
        g.clone(),
        Length::new(0.5),
        model.width().value() * model.height().value(),
        vec![g.conserved(&[density], speed, 2.0 * density).unwrap(); 2],
        Boundaries::Periodic,
        || false,
    )
    .unwrap();
    MixtureNetwork1d::new(vec![flow], Vec::new(), Vec::new(), || false).unwrap()
}

fn adiabatic() -> WallContact {
    WallContact {
        drag: WallDrag::None,
        thermal: WallThermal::Adiabatic,
    }
}

// Direct long Fourier sum from Mortensen et al. (2005), Eq. (11). Production
// instead sums six exponential differences from zeta(5); no shared helper.
fn fourier_factor(width: f64, height: f64) -> f64 {
    let aspect = width.min(height) / width.max(height);
    let sum: f64 = (0..4096)
        .map(|i| {
            let n = (2 * i + 1) as f64;
            (n * core::f64::consts::PI / (2.0 * aspect)).tanh() / n.powi(5)
        })
        .sum();
    1.0 - 192.0 * aspect * sum / core::f64::consts::PI.powi(5)
}

#[test]
fn square_poiseuille_resistance_converts_signed_kinetic_energy_without_a_second_heat_source() {
    let g = gas(1, false, 0.0);
    let model = rectangle(1.0, 1.0, 0.0);
    // Exact square-duct Darcy f*Re = 56.9083075391245..., instead of 64.
    let square_rate = 28.454_153_769_562_24;
    for (density, mu) in [(1.0, 0.004), (2.0, 0.004), (1.0, 0.012)] {
        for speed in [-0.2_f64, 0.0, 0.2] {
            let mut n = passage(&g, model, density, speed);
            let before = fields(&n);
            let mut wall = uniform_operator(&n, adiabatic())
                .with_rectangular_laminar_transport(
                    law(&g, mu, 0.0, 0.0, 0.0),
                    vec![Some(model)],
                    || false,
                )
                .unwrap();
            assert_eq!(wall.rectangular_walls(), Some([Some(model)].as_slice()));
            assert_eq!(wall.circular_walls(), None);
            assert_eq!(wall.transport_law().unwrap().source(), g.source());
            let dt = 0.2;
            let report = wall
                .advance_for(&mut n, Time::new(dt), controls(), || false)
                .unwrap();
            let rate = square_rate * mu / density;
            let velocity = speed * (-rate * dt).exp();
            for i in 0..2 {
                let state = n.branches()[0].cell(i).unwrap();
                let actual = g.primitive(state, || false).unwrap();
                close(actual.velocity, velocity, 2e-13);
                close(
                    actual.temperature,
                    2.0 + (speed * speed - velocity * velocity) / 5.0,
                    2e-12,
                );
                assert_eq!(state.partial_densities, before.0[0][i].partial_densities);
                assert_eq!(state.energy, before.0[0][i].energy);
                let WallDrag::Linear { rate_per_s } = wall.branch_contacts(0).unwrap()[i].drag
                else {
                    panic!(
                        "rectangular material resistance must use the original linear drag owner"
                    )
                };
                close(rate_per_s, rate, 2e-13);
            }
            close(
                report.drag_heating_j,
                0.5 * density * (speed * speed - velocity * velocity),
                2e-13,
            );
            assert_eq!(report.heat_to_gas_j, 0.0);
            assert_eq!(report.mechanical_work_to_gas_j, 0.0);
            close(report.energy_residual_j, 0.0, 1e-13);
        }
    }
}

#[test]
fn rectangular_side_walls_and_authored_heat_use_actual_geometry_and_survive_width_height_swap() {
    let g = gas(1, false, 0.0);
    let (width, height, nu) = (1.0, 0.4, 6.0);
    let (mu, conductivity, dt, speed) = (0.004, 0.02, 0.2, 0.2_f64);
    let rate = 12.0 * mu / (height * height * fourier_factor(width, height));
    let velocity = speed * (-rate * dt).exp();
    let heated_temperature = 2.0 + (speed * speed - velocity * velocity) / 5.0;
    let hydraulic_diameter = 2.0 * width * height / (width + height);
    let cell_volume = 0.5 * width * height;
    let wetted_area = 2.0 * (width + height) * 0.5;
    let conductance = nu * conductivity * wetted_area / hydraulic_diameter;
    let gas_capacity = cell_volume * 2.5;
    let q = dt * conductance * (5.0 - heated_temperature)
        / (1.0 + dt * conductance * (1.0 / gas_capacity + 1.0));
    let mut original = None;
    for (w, h) in [(width, height), (height, width)] {
        let model = rectangle(w, h, nu);
        let mut n = passage(&g, model, 1.0, speed);
        let mut wall = uniform_operator(&n, lumped())
            .with_rectangular_laminar_transport(
                law(&g, mu, 0.0, conductivity, 0.0),
                vec![Some(model)],
                || false,
            )
            .unwrap();
        let report = wall
            .advance_for(&mut n, Time::new(dt), controls(), || false)
            .unwrap();
        close(
            model.hydraulic_diameter().value(),
            hydraulic_diameter,
            1e-14,
        );
        for i in 0..2 {
            let actual = g
                .primitive(n.branches()[0].cell(i).unwrap(), || false)
                .unwrap();
            close(actual.velocity, velocity, 2e-13);
            close(
                actual.temperature,
                heated_temperature + q / gas_capacity,
                1e-10,
            );
            let WallThermal::Lumped {
                temperature_k,
                conductance_w_per_k,
                ..
            } = wall.branch_contacts(0).unwrap()[i].thermal
            else {
                panic!("lost finite wall")
            };
            close(temperature_k, 5.0 - q, 1e-10);
            close(conductance_w_per_k, conductance, 1e-13);
        }
        close(report.heat_to_gas_j, 2.0 * q, 1e-10);
        close(report.wall_energy_change_j, -2.0 * q, 1e-10);
        close(report.energy_residual_j, 0.0, 1e-10);
        if let Some((previous_fields, previous_report)) = &original {
            assert_eq!(fields(&n), *previous_fields);
            assert_eq!(report, *previous_report);
        } else {
            original = Some((fields(&n), report));
        }
    }
}

#[test]
fn rectangular_narrow_gap_converges_to_parallel_plate_resistance() {
    let g = gas(1, false, 0.0);
    let (h, mu) = (0.1, 0.004);
    let plate_rate = 12.0 * mu / (h * h);
    let mut previous_error = 1.0;
    for aspect in [0.1, 0.01, 0.001, 0.0001] {
        let model = rectangle(h / aspect, h, 0.0);
        let mut n = passage(&g, model, 1.0, 0.2);
        let mut wall = uniform_operator(&n, adiabatic())
            .with_rectangular_laminar_transport(
                law(&g, mu, 0.0, 0.0, 0.0),
                vec![Some(model)],
                || false,
            )
            .unwrap();
        wall.advance_for(&mut n, Time::new(0.01), controls(), || false)
            .unwrap();
        let WallDrag::Linear { rate_per_s } = wall.branch_contacts(0).unwrap()[0].drag else {
            panic!("expected rectangular Poiseuille rate")
        };
        close(
            rate_per_s,
            plate_rate / fourier_factor(h / aspect, h),
            1e-12,
        );
        let error = rate_per_s / plate_rate - 1.0;
        assert!(error > 0.0 && error < previous_error / 5.0);
        previous_error = error;
    }
    assert!(previous_error < 7e-5);
}

#[test]
fn rectangular_moving_wall_preserves_relative_resistance_and_work_under_a_boost() {
    let g = gas(1, false, 0.0);
    let model = rectangle(1.0, 0.4, 0.0);
    let (mu, dt) = (0.004, 0.2);
    let rate = 12.0 * mu / (0.16 * fourier_factor(1.0, 0.4));
    for boost in [0.0, 20.0] {
        for slip in [-0.2_f64, 0.2] {
            let wall_velocity = 0.2 + boost;
            let mut n = passage(&g, model, 1.0, wall_velocity + slip);
            let mut wall = uniform_operator(
                &n,
                WallContact {
                    drag: WallDrag::MovingLinear {
                        rate_per_s: 999.0,
                        wall_velocity_m_s: wall_velocity,
                    },
                    thermal: WallThermal::Adiabatic,
                },
            )
            .with_rectangular_laminar_transport(
                law(&g, mu, 0.0, 0.0, 0.0),
                vec![Some(model)],
                || false,
            )
            .unwrap();
            let report = wall
                .advance_for(&mut n, Time::new(dt), controls(), || false)
                .unwrap();
            let final_slip = slip * (-rate * dt).exp();
            for i in 0..2 {
                let actual = g
                    .primitive(n.branches()[0].cell(i).unwrap(), || false)
                    .unwrap();
                close(actual.velocity, wall_velocity + final_slip, 2e-13);
                close(
                    actual.temperature,
                    2.0 + (slip * slip - final_slip * final_slip) / 5.0,
                    2e-11,
                );
            }
            close(
                report.mechanical_work_to_gas_j,
                0.4 * wall_velocity * (final_slip - slip),
                2e-12,
            );
            close(
                report.drag_heating_j,
                0.2 * (slip * slip - final_slip * final_slip),
                2e-13,
            );
            close(report.energy_residual_j, 0.0, 2e-12);
        }
    }
}

#[test]
fn mixed_circular_rectangular_network_uses_each_shape_with_one_transport_owner() {
    let g = gas(1, false, 0.0);
    let rectangle = rectangle(1.0, 0.4, 0.0);
    let circle = CircularLaminarWall {
        nusselt_number: 0.0,
        ..geometry()
    };
    let models = vec![
        Some(DevelopedWallModel::Circular(CircularWallModel::Laminar(
            circle,
        ))),
        Some(DevelopedWallModel::Rectangular(rectangle)),
    ];
    let speeds = [0.2_f64, -0.2];
    let (mu, dt) = (0.004, 0.2);
    let rates = [
        32.0 * mu / circle.diameter.value().powi(2),
        12.0 * mu / (0.16 * fourier_factor(1.0, 0.4)),
    ];
    let mut n = network(&g, &speeds, false);
    let mut wall = uniform_operator(&n, adiabatic())
        .with_developed_transport(law(&g, mu, 0.0, 0.0, 0.0), models.clone(), || false)
        .unwrap();
    assert_eq!(wall.developed_walls(), Some(models.as_slice()));
    assert_eq!(wall.circular_walls(), None);
    assert_eq!(wall.rectangular_walls(), None);
    let report = wall
        .advance_for(&mut n, Time::new(dt), controls(), || false)
        .unwrap();
    let mut heat = 0.0;
    for branch in 0..2 {
        let velocity = speeds[branch] * (-rates[branch] * dt).exp();
        for cell in 0..2 {
            let actual = g
                .primitive(n.branches()[branch].cell(cell).unwrap(), || false)
                .unwrap();
            close(actual.velocity, velocity, 2e-13);
        }
        heat += 0.2 * (speeds[branch] * speeds[branch] - velocity * velocity);
    }
    close(report.drag_heating_j, heat, 2e-13);
    close(report.energy_residual_j, 0.0, 1e-13);

    // The second shape's regime refusal cannot publish the first shape's drag.
    let limited =
        RectangularLaminarWall::try_new(Length::new(1.0), Length::new(0.4), 0.0, [0.0, 1.0], 0.3)
            .unwrap();
    let mut n = network(&g, &speeds, false);
    let before = fields(&n);
    let mut wall = uniform_operator(&n, adiabatic())
        .with_developed_transport(
            law(&g, mu, 0.0, 0.0, 0.0),
            vec![models[0], Some(DevelopedWallModel::Rectangular(limited))],
            || false,
        )
        .unwrap();
    assert!(
        wall.advance_for(&mut n, Time::new(dt), controls(), || false)
            .is_err()
    );
    assert_eq!(fields(&n), before);
    assert_eq!(wall.last_report(), None);
}

#[test]
fn rectangular_geometry_source_and_reynolds_refusals_preserve_accepted_state() {
    for value in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(
            RectangularLaminarWall::try_new(
                Length::new(value),
                Length::new(0.4),
                0.0,
                [0.0, 2300.0],
                0.3,
            )
            .is_err()
        );
    }
    for (nu, band, mach) in [
        (-1.0, [0.0, 2300.0], 0.3),
        (1.0, [0.0, 2301.0], 0.3),
        (1.0, [-1.0, 2300.0], 0.3),
        (1.0, [1.0, 0.0], 0.3),
        (1.0, [0.0, 2300.0], 0.31),
    ] {
        assert!(
            RectangularLaminarWall::try_new(Length::new(1.0), Length::new(0.4), nu, band, mach)
                .is_err()
        );
    }
    let g = gas(1, false, 0.0);
    let model = rectangle(1.0, 0.4, 0.0);
    let mut n = passage(&g, model, 1.0, 0.2);
    for wrong in [rectangle(1.0, 0.5, 0.0), rectangle(0.5, 0.5, 0.0)] {
        assert!(
            uniform_operator(&n, adiabatic())
                .with_rectangular_laminar_transport(
                    law(&g, 0.004, 0.0, 0.0, 0.0),
                    vec![Some(wrong)],
                    || false,
                )
                .is_err()
        );
    }
    let other = gas(1, false, 1.0);
    assert!(
        uniform_operator(&n, adiabatic())
            .with_rectangular_laminar_transport(
                law(&other, 0.004, 0.0, 0.0, 0.0),
                vec![Some(model)],
                || false,
            )
            .is_err()
    );
    let limited =
        RectangularLaminarWall::try_new(Length::new(1.0), Length::new(0.4), 0.0, [0.0, 1.0], 0.3)
            .unwrap();
    let mut wall = uniform_operator(&n, lumped())
        .with_rectangular_laminar_transport(
            law(&g, 0.004, 0.0, 0.02, 0.0),
            vec![Some(limited)],
            || false,
        )
        .unwrap();
    let before = fields(&n);
    let contacts = wall.branch_contacts(0).unwrap().to_vec();
    assert!(matches!(
        wall.advance_for(&mut n, Time::new(0.2), controls(), || false),
        Err(WallExchangeError::TransportDomain {
            quantity: "Reynolds number",
            ..
        })
    ));
    assert_eq!(fields(&n), before);
    assert_eq!(wall.branch_contacts(0).unwrap(), contacts);
    assert_eq!(wall.last_report(), None);
}

#[test]
fn rectangular_late_cancellation_preserves_fields_contacts_and_report_then_reuses_owner() {
    let g = gas(1, false, 0.0);
    let model = rectangle(1.0, 0.4, 6.0);
    let prepare = || {
        let n = passage(&g, model, 1.0, 0.2);
        let wall = uniform_operator(&n, lumped())
            .with_rectangular_laminar_transport(
                law(&g, 0.004, 0.0, 0.02, 0.0),
                vec![Some(model)],
                || false,
            )
            .unwrap();
        (n, wall)
    };
    let (mut reference, mut reference_wall) = prepare();
    let mut calls = 0;
    let expected = reference_wall
        .advance_for(&mut reference, Time::new(0.2), controls(), || {
            calls += 1;
            false
        })
        .unwrap();
    assert!(calls > 8);
    let (mut n, mut wall) = prepare();
    let before = fields(&n);
    let contacts = wall.branch_contacts(0).unwrap().to_vec();
    let mut cancelled_calls = 0;
    assert!(
        wall.advance_for(&mut n, Time::new(0.2), controls(), || {
            cancelled_calls += 1;
            cancelled_calls >= 3 * calls / 4
        })
        .is_err()
    );
    assert_eq!(fields(&n), before);
    assert_eq!(wall.branch_contacts(0).unwrap(), contacts);
    assert_eq!(wall.last_report(), None);
    let accepted = wall
        .advance_for(&mut n, Time::new(0.2), controls(), || false)
        .unwrap();
    assert_eq!(accepted, expected);
    assert_eq!(fields(&n), fields(&reference));
}
