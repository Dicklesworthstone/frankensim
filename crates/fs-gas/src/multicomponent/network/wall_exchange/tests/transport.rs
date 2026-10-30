//! G1 analytic wall exchange and G3 source/domain substitution in the real owner.
use super::*;
use crate::multicomponent::dissipation::{ConstantTransport, TransportDomain};
use crate::multicomponent::network::dissipation::tests::law;
use crate::multicomponent::network::wall_exchange::CircularLaminarWall;

fn geometry() -> CircularLaminarWall {
    CircularLaminarWall {
        diameter: Length::new((0.4 / core::f64::consts::FRAC_PI_4).sqrt()),
        nusselt_number: 3.66,
        reynolds_number: [0.0, 2300.0],
        maximum_mach_number: 0.3,
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
        assert_eq!(wall.circular_laminar_walls().unwrap(), &[Some(geometry())]);
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
