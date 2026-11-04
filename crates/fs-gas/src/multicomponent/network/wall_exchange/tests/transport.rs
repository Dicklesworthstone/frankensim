//! G1 analytic wall exchange and G3 source/domain substitution in the real owner.
use super::*;
use crate::multicomponent::dissipation::{ConstantTransport, TransportDomain};
use crate::multicomponent::network::dissipation::tests::law;
use crate::multicomponent::network::wall_exchange::{
    CircularLaminarWall, CircularTurbulentWall, CircularWallModel,
};

fn geometry() -> CircularLaminarWall {
    CircularLaminarWall {
        diameter: Length::new((0.4 / core::f64::consts::FRAC_PI_4).sqrt()),
        nusselt_number: 3.66,
        reynolds_number: [0.0, 2300.0],
        maximum_mach_number: 0.3,
    }
}

#[test]
fn moving_wall_material_refresh_uses_relative_flow_and_preserves_motion() {
    let g = gas(1, false, 0.0);
    let (mu, dt) = (0.004, 0.2);
    let diameter = geometry().diameter.value();
    let rate = 32.0 * mu / (diameter * diameter);
    for boost in [0.0, 20.0] {
        let speed = 0.2 + boost;
        let initial_velocity = 0.4 + boost;
        let mut n = network(&g, &[initial_velocity], false);
        let mut wall = uniform_operator(
            &n,
            WallContact {
                drag: WallDrag::MovingLinear {
                    rate_per_s: 999.0,
                    wall_velocity_m_s: speed,
                },
                thermal: WallThermal::Adiabatic,
            },
        )
        .with_circular_laminar_transport(law(&g, mu, 0.0, 0.0, 0.0), vec![Some(geometry())], || {
            false
        })
        .unwrap();
        let report = wall
            .advance_for(&mut n, Time::new(dt), controls(), || false)
            .unwrap();
        let slip = 0.2 * (-rate * dt).exp();
        let expected_temperature = 2.0 + (0.04 - slip * slip) / 5.0;
        for i in 0..2 {
            let state = n.branches()[0].cell(i).unwrap();
            let actual = g.primitive(state, || false).unwrap();
            close(actual.velocity, speed + slip, 3e-13);
            close(actual.temperature, expected_temperature, 3e-11);
            match wall.branch_contacts(0).unwrap()[i].drag {
                WallDrag::MovingLinear {
                    rate_per_s,
                    wall_velocity_m_s,
                } => {
                    close(rate_per_s, rate, 3e-13);
                    assert_eq!(wall_velocity_m_s, speed);
                }
                _ => panic!("material coefficient refresh discarded prescribed wall motion"),
            }
        }
        close(
            report.mechanical_work_to_gas_j,
            0.4 * speed * (slip - 0.2),
            3e-12,
        );
        close(report.drag_heating_j, 0.2 * (0.04 - slip * slip), 3e-13);
    }
}

#[test]
fn current_density_and_source_coefficients_drive_drag_heat_and_wall_storage() {
    let g = gas(1, false, 0.0);
    let diameter = geometry().diameter.value();
    let dt = 0.2;
    let speed = 0.2;
    // Density changes at fixed temperature, and viscosity/conductivity changes
    // independently, all use exactly the same existing conservative solver.
    for (density, mu, conductivity) in [(1.0, 0.004, 0.02), (2.0, 0.004, 0.02), (1.0, 0.012, 0.1)] {
        let flow = MixtureEuler1d::new(
            g.clone(),
            Length::new(0.5),
            0.4,
            vec![g.conserved(&[density], speed, 2.0 * density).unwrap(); 2],
            Boundaries::Periodic,
            || false,
        )
        .unwrap();
        let mut n = MixtureNetwork1d::new(vec![flow], Vec::new(), Vec::new(), || false).unwrap();
        let before = fields(&n);
        let mut wall = uniform_operator(&n, lumped())
            .with_circular_laminar_transport(
                law(&g, mu, 0.0, conductivity, 0.0),
                vec![Some(geometry())],
                || false,
            )
            .unwrap();
        assert_eq!(wall.transport_law().unwrap().source(), g.source());
        assert_eq!(
            wall.circular_walls().unwrap(),
            &[Some(CircularWallModel::Laminar(geometry()))]
        );
        let report = wall
            .advance_for(&mut n, Time::new(dt), controls(), || false)
            .unwrap();

        let rate = 32.0 * mu / (density * diameter * diameter);
        let final_speed = speed * (-rate * dt).exp();
        let heated_temperature = 2.0 + (speed * speed - final_speed * final_speed) / 5.0;
        let gas_capacity = density * 0.2 * 2.5;
        let conductance = 3.66 * conductivity * (4.0 * 0.2 / diameter) / diameter;
        let q = dt * conductance * (5.0 - heated_temperature)
            / (1.0 + dt * conductance * (1.0 / gas_capacity + 1.0));
        for i in 0..2 {
            let state = n.branches()[0].cell(i).unwrap();
            let primitive = g.primitive(state, || false).unwrap();
            close(primitive.velocity, final_speed, 1e-12);
            close(
                primitive.temperature,
                heated_temperature + q / gas_capacity,
                1e-10,
            );
            assert_eq!(state.partial_densities, before.0[0][i].partial_densities);
            close((state.energy - before.0[0][i].energy) * 0.2, q, 1e-10);
            let contact = wall.branch_contacts(0).unwrap()[i];
            match contact.drag {
                WallDrag::Linear { rate_per_s } => close(rate_per_s, rate, 1e-13),
                _ => panic!("expected current material-derived wall drag"),
            }
            match contact.thermal {
                WallThermal::Lumped {
                    temperature_k,
                    conductance_w_per_k,
                    ..
                } => {
                    close(temperature_k, 5.0 - q, 1e-10);
                    close(conductance_w_per_k, conductance, 1e-13);
                }
                _ => panic!("expected original finite wall storage"),
            }
        }
        close(report.heat_to_gas_j, 2.0 * q, 1e-10);
        close(report.wall_energy_change_j, -2.0 * q, 1e-10);
        close(report.energy_residual_j, 0.0, 1e-10);
        assert_eq!(report.external_heat_to_gas_j, 0.0);
        assert_eq!(n.time(), Time::new(0.0));
    }
}

#[test]
fn source_identity_geometry_and_current_flow_applicability_refuse_without_publication() {
    let g = gas(1, false, 0.0);
    let mut n = network(&g, &[0.2], false);
    // Same label, different complete NASA source: no name-only binding.
    let other = gas(1, false, 1.0);
    assert!(matches!(
        uniform_operator(&n, lumped()).with_circular_laminar_transport(
            law(&other, 0.004, 0.0, 0.02, 0.0),
            vec![Some(geometry())],
            || false,
        ),
        Err(WallExchangeError::Invalid(_))
    ));
    let mut wrong = geometry();
    wrong.diameter = Length::new(1.1 * wrong.diameter.value());
    assert!(matches!(
        uniform_operator(&n, lumped()).with_circular_laminar_transport(
            law(&g, 0.004, 0.0, 0.02, 0.0),
            vec![Some(wrong)],
            || false,
        ),
        Err(WallExchangeError::Invalid(_))
    ));
    let mut narrow = geometry();
    narrow.reynolds_number = [0.0, 1.0];
    let mut wall = uniform_operator(&n, lumped())
        .with_circular_laminar_transport(law(&g, 0.004, 0.0, 0.02, 0.0), vec![Some(narrow)], || {
            false
        })
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
fn source_domain_at_heated_endpoint_rolls_back_earlier_drag_and_wall_candidates() {
    let g = gas(1, false, 0.0);
    let mut n = network(&g, &[0.2], false);
    let coefficient = law(&g, 0.004, 0.0, 1.0, 0.0);
    let limited = ConstantTransport::new(
        g.source().clone(),
        "bounded-wall-properties",
        1,
        TransportDomain {
            minimum_temperature_k: 0.001,
            maximum_temperature_k: 2.2,
            minimum_pressure_pa: 1e-8,
            maximum_pressure_pa: 1e8,
        },
        coefficient.coefficients(),
        vec![0.0],
    )
    .unwrap();
    let adiabatic = WallContact {
        drag: WallDrag::None,
        thermal: WallThermal::Adiabatic,
    };
    let mut wall =
        NetworkWallOperator::new(&n, vec![vec![adiabatic, lumped()]], Vec::new(), || false)
            .unwrap()
            .with_circular_laminar_transport(limited, vec![Some(geometry())], || false)
            .unwrap();
    // Publish an admitted zero-duration source report first; the later refusal
    // must preserve that report as well as both accepted cells and wall states.
    wall.advance_for(&mut n, Time::new(0.0), controls(), || false)
        .unwrap();
    let before = fields(&n);
    let contacts = wall.branch_contacts(0).unwrap().to_vec();
    let report = wall.last_report();
    assert!(matches!(
        wall.advance_for(&mut n, Time::new(0.2), controls(), || false),
        Err(WallExchangeError::Transport(_))
    ));
    assert_eq!(fields(&n), before);
    assert_eq!(wall.branch_contacts(0).unwrap(), contacts);
    assert_eq!(wall.last_report(), report);
    assert_eq!(n.time(), Time::new(0.0));
}

fn turbulent_geometry(roughness: f64) -> CircularTurbulentWall {
    CircularTurbulentWall {
        diameter: geometry().diameter,
        equivalent_sand_roughness: Length::new(roughness * geometry().diameter.value()),
        nusselt_number: 3.66,
        reynolds_number: [4000.0, 1e8],
        maximum_mach_number: 0.3,
    }
}

// Independent bisection in the Darcy factor itself, with platform logarithms.
fn reference_colebrook(reynolds: f64, roughness: f64) -> f64 {
    let (mut lo, mut hi): (f64, f64) = (0.001, 0.2);
    for _ in 0..100 {
        let f = 0.5 * (lo + hi);
        let residual =
            1.0 / f.sqrt() + 2.0 * (roughness / 3.7 + 2.51 / (reynolds * f.sqrt())).log10();
        if residual > 0.0 {
            lo = f;
        } else {
            hi = f;
        }
    }
    0.5 * (lo + hi)
}

#[test]
fn turbulent_wall_g1_current_material_and_roughness_drive_native_work_and_heat() {
    let g = gas(1, false, 0.0);
    let dt = 0.4;
    let diameter = geometry().diameter.value();
    for (mu, roughness, wall_speed, slip) in [
        (1e-6, 0.0, 0.0, 0.2),
        (2e-6, 0.0, 0.0, 0.2),
        (1e-6, 0.01, 0.0, 0.2),
        (1e-6, 0.01, 20.0, 0.2),
        (1e-6, 0.01, 20.0, -0.2),
    ] {
        let mut n = network(&g, &[wall_speed + slip], false);
        let before = fields(&n);
        let model = turbulent_geometry(roughness);
        let mut wall = uniform_operator(
            &n,
            WallContact {
                drag: WallDrag::MovingDarcy {
                    friction_factor: 999.0,
                    hydraulic_diameter_m: diameter,
                    wall_velocity_m_s: wall_speed,
                },
                ..lumped()
            },
        )
        .with_circular_transport(
            law(&g, mu, 0.0, 0.02, 0.0),
            vec![Some(CircularWallModel::Turbulent(model))],
            || false,
        )
        .unwrap();
        let report = wall
            .advance_for(&mut n, Time::new(dt), controls(), || false)
            .unwrap();
        let factor = reference_colebrook(slip.abs() * diameter / mu, roughness);
        let final_slip = slip / (1.0 + 0.5 * factor * dt * slip.abs() / diameter);
        let heated_temperature = 2.0 + (slip * slip - final_slip * final_slip) / 5.0;
        let conductance = model.nusselt_number * 0.02 * 4.0 * 0.2 / diameter.powi(2);
        let q =
            dt * conductance * (5.0 - heated_temperature) / (1.0 + dt * conductance * (2.0 + 1.0));
        for i in 0..2 {
            let cell = n.branches()[0].cell(i).unwrap();
            let primitive = g.primitive(cell, || false).unwrap();
            close(primitive.velocity, wall_speed + final_slip, 3e-12);
            close(primitive.temperature, heated_temperature + 2.0 * q, 3e-10);
            assert_eq!(cell.partial_densities, before.0[0][i].partial_densities);
            let contact = wall.branch_contacts(0).unwrap()[i];
            let WallDrag::MovingDarcy {
                friction_factor,
                wall_velocity_m_s,
                ..
            } = contact.drag
            else {
                panic!("turbulent source discarded the original moving-wall law");
            };
            close(friction_factor, factor, 2e-13);
            assert_eq!(wall_velocity_m_s, wall_speed);
            close(wall_temperature(contact), 5.0 - q, 3e-10);
        }
        close(
            report.mechanical_work_to_gas_j,
            0.4 * wall_speed * (final_slip - slip),
            3e-12,
        );
        close(
            report.drag_heating_j,
            0.2 * (slip * slip - final_slip * final_slip),
            3e-12,
        );
        close(report.heat_to_gas_j, 2.0 * q, 3e-10);
        close(report.wall_energy_change_j, -2.0 * q, 3e-10);
        close(report.energy_residual_j, 0.0, 3e-10);
        assert_eq!(n.time(), Time::new(0.0));
    }
}

#[test]
fn turbulent_wall_g3_mixed_regimes_share_the_original_source_and_refresh_each_step() {
    let g = gas(1, false, 0.0);
    let mu = 1e-6;
    let diameter = geometry().diameter.value();
    let model = turbulent_geometry(0.001);
    let mut n = network(&g, &[0.001, 0.2], false);
    let mut wall = uniform_operator(
        &n,
        WallContact {
            drag: WallDrag::None,
            thermal: WallThermal::Adiabatic,
        },
    )
    .with_circular_transport(
        law(&g, mu, 0.0, 0.0, 0.0),
        vec![
            Some(CircularWallModel::Laminar(geometry())),
            Some(CircularWallModel::Turbulent(model)),
        ],
        || false,
    )
    .unwrap();
    let dt = 10.0;
    let mut expected_speed: f64 = 0.2;
    for _ in 0..3 {
        let factor = reference_colebrook(expected_speed * diameter / mu, 0.001);
        expected_speed /= 1.0 + 0.5 * factor * dt * expected_speed / diameter;
        wall.advance_for(&mut n, Time::new(dt), controls(), || false)
            .unwrap();
        let actual = g
            .primitive(n.branches()[1].cell(0).unwrap(), || false)
            .unwrap();
        close(actual.velocity, expected_speed, 2e-13);
        let WallDrag::Darcy {
            friction_factor, ..
        } = wall.branch_contacts(1).unwrap()[0].drag
        else {
            panic!("missing resolved Darcy source");
        };
        close(friction_factor, factor, 2e-13);
    }
    let laminar_speed = g
        .primitive(n.branches()[0].cell(0).unwrap(), || false)
        .unwrap()
        .velocity;
    close(
        laminar_speed,
        0.001 * (-32.0 * mu * 3.0 * dt / diameter.powi(2)).exp(),
        2e-14,
    );
}

#[test]
fn turbulent_wall_g4_domain_exit_preserves_all_accepted_cells_and_reports() {
    let g = gas(1, false, 0.0);
    let mu = 1e-6;
    let mut n = network(&g, &[0.2], false);
    let mut model = turbulent_geometry(0.01);
    let initial_re = 0.2 * model.diameter.value() / mu;
    model.reynolds_number = [0.9 * initial_re, 2.0 * initial_re];
    let mut wall = uniform_operator(&n, lumped())
        .with_circular_transport(
            law(&g, mu, 0.0, 0.02, 0.0),
            vec![Some(CircularWallModel::Turbulent(model))],
            || false,
        )
        .unwrap();
    wall.advance_for(&mut n, Time::new(0.0), controls(), || false)
        .unwrap();
    let before = fields(&n);
    let contacts = wall.branch_contacts(0).unwrap().to_vec();
    let report = wall.last_report();
    assert!(matches!(
        wall.advance_for(&mut n, Time::new(100.0), controls(), || false),
        Err(WallExchangeError::TransportDomain {
            quantity: "Reynolds number",
            ..
        })
    ));
    assert_eq!(fields(&n), before);
    assert_eq!(wall.branch_contacts(0).unwrap(), contacts);
    assert_eq!(wall.last_report(), report);
}

#[test]
fn turbulent_wall_g1_refinement_resolves_reynolds_dependent_braking() {
    let g = gas(1, false, 0.0);
    let mu = 1e-6;
    let diameter = geometry().diameter.value();
    let duration = 100.0;
    // Independent RK4 on du/dt=-f(Re(u))*u*|u|/(2D). This reference refreshes
    // f inside every stage; production freezes it over an exact Darcy interval.
    let derivative = |v: f64| {
        -0.5 * reference_colebrook(v.abs() * diameter / mu, 0.001) * v * v.abs() / diameter
    };
    let mut reference = 0.2;
    let h = duration / 512.0;
    for _ in 0..512 {
        let k1 = derivative(reference);
        let k2 = derivative(reference + 0.5 * h * k1);
        let k3 = derivative(reference + 0.5 * h * k2);
        let k4 = derivative(reference + h * k3);
        reference += h * (k1 + 2.0 * k2 + 2.0 * k3 + k4) / 6.0;
    }
    let mut errors = Vec::new();
    for steps in [8, 16, 32] {
        let mut n = network(&g, &[0.2], false);
        let mut wall = uniform_operator(
            &n,
            WallContact {
                drag: WallDrag::None,
                thermal: WallThermal::Adiabatic,
            },
        )
        .with_circular_transport(
            law(&g, mu, 0.0, 0.0, 0.0),
            vec![Some(CircularWallModel::Turbulent(turbulent_geometry(
                0.001,
            )))],
            || false,
        )
        .unwrap();
        for _ in 0..steps {
            wall.advance_for(
                &mut n,
                Time::new(duration / f64::from(steps)),
                controls(),
                || false,
            )
            .unwrap();
        }
        let actual = g
            .primitive(n.branches()[0].cell(0).unwrap(), || false)
            .unwrap()
            .velocity;
        errors.push((actual - reference).abs());
    }
    assert!(errors[2] < 1.3e-5);
    for pair in errors.windows(2) {
        assert!((1.9..2.1).contains(&(pair[0] / pair[1])), "{errors:?}");
    }
}
