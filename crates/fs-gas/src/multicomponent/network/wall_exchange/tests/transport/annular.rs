//! Original gas/network owners consuming annular resistance and both wall areas.
use super::*;
use crate::multicomponent::network::wall_exchange::{
    AnnularLaminarWall, DevelopedWallModel, RectangularLaminarWall,
};

fn annulus(inner: f64, outer: f64, nu: f64) -> AnnularLaminarWall {
    AnnularLaminarWall::try_new(
        Length::new(inner),
        Length::new(outer),
        nu,
        [0.0, 2300.0],
        0.3,
    )
    .unwrap()
}

fn rate(inner: f64, outer: f64, mu: f64, density: f64) -> f64 {
    let difference = outer * outer - inner * inner;
    let flux_factor =
        outer.powi(4) - inner.powi(4) - difference * difference / (outer / inner).ln();
    8.0 * mu * difference / (density * flux_factor)
}

#[test]
fn annular_heat_includes_both_surfaces_and_preserves_original_gas_and_wall_storage() {
    let g = gas(1, false, 0.0);
    let (inner, outer, nu, mu, conductivity, dt) = (0.5, 1.0, 5.0, 0.004, 0.02, 0.2);
    let model = annulus(inner, outer, nu);
    let speed = -0.2_f64;
    let flow = MixtureEuler1d::new(
        g.clone(),
        Length::new(0.5),
        model.flow_area_m2(),
        vec![g.conserved(&[1.0], speed, 2.0).unwrap(); 2],
        Boundaries::Periodic,
        || false,
    )
    .unwrap();
    let mut n = MixtureNetwork1d::new(vec![flow], Vec::new(), Vec::new(), || false).unwrap();
    let before = fields(&n);
    let mut wall = uniform_operator(&n, lumped())
        .with_annular_laminar_transport(
            law(&g, mu, 0.0, conductivity, 0.0),
            vec![Some(model)],
            || false,
        )
        .unwrap();
    assert_eq!(wall.annular_walls(), Some([Some(model)].as_slice()));
    let report = wall
        .advance_for(&mut n, Time::new(dt), controls(), || false)
        .unwrap();
    let velocity = speed * (-rate(inner, outer, mu, 1.0) * dt).exp();
    let heated = 2.0 + (speed * speed - velocity * velocity) / 5.0;
    let volume = core::f64::consts::PI * (outer * outer - inner * inner) * 0.5;
    let gas_capacity = volume * 2.5;
    let wetted_area = core::f64::consts::TAU * (outer + inner) * 0.5;
    let conductance = nu * conductivity * wetted_area / (2.0 * (outer - inner));
    let q =
        dt * conductance * (5.0 - heated) / (1.0 + dt * conductance * (1.0 / gas_capacity + 1.0));
    for cell in 0..2 {
        let state = n.branches()[0].cell(cell).unwrap();
        let actual = g.primitive(state, || false).unwrap();
        close(actual.velocity, velocity, 2e-13);
        close(actual.temperature, heated + q / gas_capacity, 1e-10);
        assert_eq!(state.partial_densities, before.0[0][cell].partial_densities);
        close((state.energy - before.0[0][cell].energy) * volume, q, 1e-10);
        let WallThermal::Lumped {
            temperature_k,
            conductance_w_per_k,
            ..
        } = wall.branch_contacts(0).unwrap()[cell].thermal
        else {
            panic!("lost original finite wall")
        };
        close(temperature_k, 5.0 - q, 1e-10);
        close(conductance_w_per_k, conductance, 2e-13);
    }
    close(report.heat_to_gas_j, 2.0 * q, 1e-10);
    close(report.wall_energy_change_j, -2.0 * q, 1e-10);
    close(
        report.drag_heating_j,
        volume * (speed * speed - velocity * velocity),
        2e-13,
    );
    close(report.energy_residual_j, 0.0, 1e-10);
}

#[test]
fn mixed_annular_circular_rectangular_resistance_retains_relative_motion_and_atomic_refusal() {
    let g = gas(1, false, 0.0);
    let inner = 0.2_f64;
    let outer = (0.4 / core::f64::consts::PI + inner * inner).sqrt();
    let annular = annulus(inner, outer, 0.0);
    let rectangle = RectangularLaminarWall::try_new(
        Length::new(1.0),
        Length::new(0.4),
        0.0,
        [0.0, 2300.0],
        0.3,
    )
    .unwrap();
    let circle = CircularLaminarWall {
        nusselt_number: 0.0,
        ..geometry()
    };
    let models = vec![
        Some(DevelopedWallModel::Circular(CircularWallModel::Laminar(
            circle,
        ))),
        Some(DevelopedWallModel::Rectangular(rectangle)),
        Some(DevelopedWallModel::Annular(annular)),
    ];
    let (mu, dt, slip) = (0.004, 0.2, -0.2_f64);
    let final_slip = slip * (-rate(inner, outer, mu, 1.0) * dt).exp();
    for boost in [0.0, 20.0] {
        let wall_velocity = 0.2 + boost;
        let mut n = network(
            &g,
            &[
                wall_velocity + 0.1,
                wall_velocity + 0.15,
                wall_velocity + slip,
            ],
            false,
        );
        let contact = WallContact {
            drag: WallDrag::MovingLinear {
                rate_per_s: 999.0,
                wall_velocity_m_s: wall_velocity,
            },
            thermal: WallThermal::Adiabatic,
        };
        let mut wall = uniform_operator(&n, contact)
            .with_developed_transport(law(&g, mu, 0.0, 0.0, 0.0), models.clone(), || false)
            .unwrap();
        assert_eq!(wall.developed_walls(), Some(models.as_slice()));
        assert_eq!(wall.annular_walls(), None);
        let before = fields(&n);
        let report = wall
            .advance_for(&mut n, Time::new(dt), controls(), || false)
            .unwrap();
        for cell in 0..2 {
            let state = n.branches()[2].cell(cell).unwrap();
            let actual = g.primitive(state, || false).unwrap();
            close(actual.velocity, wall_velocity + final_slip, 2e-13);
            close(
                actual.temperature,
                2.0 + (slip * slip - final_slip * final_slip) / 5.0,
                2e-11,
            );
            close(
                state.energy - before.0[2][cell].energy,
                wall_velocity * (state.momentum - before.0[2][cell].momentum),
                2e-12,
            );
        }
        assert!(report.drag_heating_j >= 0.0);
        close(report.energy_residual_j, 0.0, 2e-11);
    }

    let limited =
        AnnularLaminarWall::try_new(Length::new(inner), Length::new(outer), 0.0, [0.0, 1.0], 0.3)
            .unwrap();
    let mut n = network(&g, &[0.2, 0.2, 0.2], false);
    let before = fields(&n);
    let mut wall = uniform_operator(&n, lumped())
        .with_developed_transport(
            law(&g, mu, 0.0, 0.02, 0.0),
            vec![
                models[0],
                models[1],
                Some(DevelopedWallModel::Annular(limited)),
            ],
            || false,
        )
        .unwrap();
    let contacts = (0..3)
        .map(|b| wall.branch_contacts(b).unwrap().to_vec())
        .collect::<Vec<_>>();
    assert!(matches!(
        wall.advance_for(&mut n, Time::new(dt), controls(), || false),
        Err(WallExchangeError::TransportDomain {
            quantity: "Reynolds number",
            ..
        })
    ));
    assert_eq!(fields(&n), before);
    for (branch, original) in contacts.iter().enumerate() {
        assert_eq!(wall.branch_contacts(branch).unwrap(), original);
    }
    assert_eq!(wall.last_report(), None);

    // Geometry belongs to the original mesh, not only to a plausible diameter.
    assert!(
        uniform_operator(&n, lumped())
            .with_annular_laminar_transport(
                law(&g, mu, 0.0, 0.02, 0.0),
                vec![Some(annulus(0.5, 1.0, 0.0)); 3],
                || false,
            )
            .is_err()
    );
}
