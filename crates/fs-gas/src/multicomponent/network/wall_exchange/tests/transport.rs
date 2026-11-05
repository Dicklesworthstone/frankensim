//! G1 analytic wall exchange and G3 source/domain substitution in the real owner.
use super::*;
use crate::multicomponent::dissipation::{ConstantTransport, TransportDomain};
use crate::multicomponent::network::dissipation::tests::law;
use crate::multicomponent::network::wall_exchange::{
    CircularChurchillWall, CircularGnielinskiWall, CircularLaminarWall, CircularTurbulentWall,
    CircularWallModel,
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

fn gnielinski_geometry() -> CircularGnielinskiWall {
    CircularGnielinskiWall {
        diameter: Length::new(0.02),
        reynolds_number: [1e4, 1e6],
        prandtl_number: [0.5, 200.0],
        maximum_mach_number: 0.3,
    }
}

fn thermal_pipe(g: &MixtureGas, densities: &[f64], velocities: [f64; 2]) -> MixtureNetwork1d {
    use fs_thermochem::UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R;
    let pressure = 2.0
        * densities
            .iter()
            .zip(g.source().models())
            .map(|(rho, model)| rho * R / model.molar_mass().value())
            .sum::<f64>();
    let flow = MixtureEuler1d::new(
        g.clone(),
        Length::new(0.2),
        core::f64::consts::FRAC_PI_4 * 0.02_f64.powi(2),
        velocities
            .into_iter()
            .map(|u| g.conserved(densities, u, pressure).unwrap())
            .collect(),
        Boundaries::Periodic,
        || false,
    )
    .unwrap();
    // The original mesh is twenty diameters long, with two equal wall areas.
    MixtureNetwork1d::new(vec![flow], Vec::new(), Vec::new(), || false).unwrap()
}

fn thermal_contact() -> WallContact {
    WallContact {
        drag: WallDrag::None,
        thermal: WallThermal::Lumped {
            temperature_k: 2.2,
            heat_capacity_j_per_k: 5e-4,
            conductance_w_per_k: 999.0,
        },
    }
}

// Independent platform-math correlation; the physical source tests below also
// solve the frozen Darcy and two finite heat-capacity equations independently.
fn reference_gnielinski(re: f64, pr: f64) -> (f64, f64) {
    let f = (0.79 * re.ln() - 1.64).powi(-2);
    let nu = (f / 8.0) * (re - 1000.0) * pr
        / (1.0 + 12.7 * (f / 8.0).sqrt() * (pr.cbrt().powi(2) - 1.0));
    (f, nu)
}

#[test]
fn gnielinski_g1_live_slip_drives_heat_and_the_same_darcy_factor_without_double_heating() {
    let g = gas(1, false, 0.0);
    let (mu, conductivity, dt, diameter) = (2e-7, 7e-7, 0.2, 0.02);
    let mut heating = Vec::new();
    for (wall_speed, slip) in [(0.0, 0.2), (0.0, 0.4), (20.0, 0.2), (20.0, -0.2)] {
        let mut n = thermal_pipe(&g, &[1.0], [wall_speed + slip; 2]);
        let before = fields(&n);
        let volume = n.branches()[0].mesh().volumes()[0];
        let mut wall = uniform_operator(
            &n,
            WallContact {
                drag: WallDrag::MovingDarcy {
                    friction_factor: 999.0,
                    hydraulic_diameter_m: diameter,
                    wall_velocity_m_s: wall_speed,
                },
                ..thermal_contact()
            },
        )
        .with_circular_transport(
            law(&g, mu, 0.0, conductivity, 0.0),
            vec![Some(CircularWallModel::Gnielinski(gnielinski_geometry()))],
            || false,
        )
        .unwrap();
        let report = wall
            .advance_for(&mut n, Time::new(dt), controls(), || false)
            .unwrap();
        // Cp=3.5 enters Pr, whereas Cv=2.5 determines the gas energy storage.
        let (factor, nu) =
            reference_gnielinski(slip.abs() * diameter / mu, mu * 3.5 / conductivity);
        let final_slip = slip / (1.0 + factor * dt * slip.abs() / (2.0 * diameter));
        let heated_temperature = 2.0 + (slip * slip - final_slip * final_slip) / 5.0;
        let conductance = nu * conductivity * 4.0 * volume / diameter.powi(2);
        let q = dt * conductance * (2.2 - heated_temperature)
            / (1.0 + dt * conductance * (1.0 / (2.5 * volume) + 1.0 / 5e-4));
        for i in 0..2 {
            let cell = n.branches()[0].cell(i).unwrap();
            let actual = g.primitive(cell, || false).unwrap();
            close(actual.velocity, wall_speed + final_slip, 3e-12);
            close(
                actual.temperature,
                heated_temperature + q / (2.5 * volume),
                2e-8,
            );
            assert_eq!(cell.partial_densities, before.0[0][i].partial_densities);
            let contact = wall.branch_contacts(0).unwrap()[i];
            let WallDrag::MovingDarcy {
                friction_factor,
                wall_velocity_m_s,
                ..
            } = contact.drag
            else {
                panic!("coupled heat closure discarded the moving Darcy source");
            };
            close(friction_factor, factor, 2e-13);
            assert_eq!(wall_velocity_m_s, wall_speed);
            let WallThermal::Lumped {
                temperature_k,
                conductance_w_per_k,
                ..
            } = contact.thermal
            else {
                panic!("coupled heat closure discarded finite wall storage");
            };
            close(conductance_w_per_k / conductance, 1.0, 2e-12);
            close(temperature_k, 2.2 - q / 5e-4, 2e-8);
            // Total gas energy gains only wall heat and prescribed wall work;
            // conversion of its existing slip kinetic energy is not a source.
            close(
                (cell.energy - before.0[0][i].energy) * volume,
                q + volume * wall_speed * (final_slip - slip),
                3e-12,
            );
        }
        close(
            report.drag_heating_j,
            volume * (slip * slip - final_slip * final_slip),
            3e-12,
        );
        close(
            report.mechanical_work_to_gas_j,
            2.0 * volume * wall_speed * (final_slip - slip),
            3e-12,
        );
        close(report.heat_to_gas_j, 2.0 * q, 3e-12);
        close(report.wall_energy_change_j, -2.0 * q, 3e-12);
        close(report.energy_residual_j, 0.0, 3e-12);
        assert_eq!(report.external_heat_to_gas_j, 0.0);
        heating.push(report.heat_to_gas_j);
    }
    assert!(
        heating[1] > 1.5 * heating[0],
        "turbulent heat did not respond to flow speed"
    );
    close(heating[2], heating[0], 3e-12);
    close(heating[3], heating[0], 3e-12);
}

#[test]
fn gnielinski_g3_original_mixture_cp_viscosity_and_conductivity_set_heat_transfer() {
    let g = gas(2, true, 0.0);
    let mut conductances = Vec::new();
    for (y, mu, conductivity) in [
        ([1.0, 0.0], 2e-7, 7e-7),
        ([0.0, 1.0], 2e-7, 7e-7),
        ([0.3, 0.7], 2e-7, 7e-7),
        ([0.3, 0.7], 3e-7, 7e-7),
        ([0.3, 0.7], 2e-7, 9e-7),
    ] {
        let mut n = thermal_pipe(&g, &y, [0.2; 2]);
        let volume = n.branches()[0].mesh().volumes()[0];
        let mut wall = uniform_operator(&n, thermal_contact())
            .with_circular_transport(
                law(&g, mu, 0.0, conductivity, 0.01),
                vec![Some(CircularWallModel::Gnielinski(gnielinski_geometry()))],
                || false,
            )
            .unwrap();
        wall.advance_for(&mut n, Time::new(0.1), controls(), || false)
            .unwrap();
        // Original NASA polynomials at T=2: Cp_A=3.5+0.1*T and
        // Cp_B=(4.3+0.05*T)/2, using the actual species molar masses.
        let cp = y[0] * 3.7 + y[1] * 2.2;
        let (_, nu) = reference_gnielinski(0.2 * 0.02 / mu, mu * cp / conductivity);
        let expected = nu * conductivity * 4.0 * volume / 0.02_f64.powi(2);
        let WallThermal::Lumped {
            conductance_w_per_k,
            ..
        } = wall.branch_contacts(0).unwrap()[0].thermal
        else {
            panic!("missing original thermal wall");
        };
        close(conductance_w_per_k / expected, 1.0, 2e-12);
        conductances.push(conductance_w_per_k);
    }
    assert!(conductances[0] > conductances[1]);
    for i in [3, 4] {
        assert!((conductances[i] / conductances[2] - 1.0).abs() > 0.01);
    }
}

#[test]
fn gnielinski_g4_reynolds_or_prandtl_endpoint_exit_rolls_back_every_cell_and_report() {
    for quantity in ["Reynolds number", "Prandtl number"] {
        let variable = quantity == "Prandtl number";
        let g = gas(1, variable, 0.0);
        let velocities = if variable { [0.2; 2] } else { [0.4, 0.2] };
        let mut n = thermal_pipe(&g, &[1.0], velocities);
        let mut model = gnielinski_geometry();
        if variable {
            model.prandtl_number[1] = 2e-7 * 3.7 / 7e-7 + 0.0001;
        } else {
            model.reynolds_number[0] = 19800.0;
        }
        // The first cell completes an admitted drag candidate. Only the second
        // cell's evolved endpoint leaves the declared correlation domain.
        let mut wall = NetworkWallOperator::new(
            &n,
            vec![vec![
                WallContact {
                    drag: WallDrag::None,
                    thermal: WallThermal::Adiabatic,
                },
                thermal_contact(),
            ]],
            Vec::new(),
            || false,
        )
        .unwrap()
        .with_circular_transport(
            law(&g, 2e-7, 0.0, 7e-7, 0.0),
            vec![Some(CircularWallModel::Gnielinski(model))],
            || false,
        )
        .unwrap();
        wall.advance_for(&mut n, Time::new(0.0), controls(), || false)
            .unwrap();
        let before = fields(&n);
        let contacts = wall.branch_contacts(0).unwrap().to_vec();
        let report = wall.last_report();
        let impulses = wall.accepted_branch_impulses().unwrap().to_vec();
        assert!(
            matches!(wall.advance_for(&mut n, Time::new(0.4), controls(), || false),
            Err(WallExchangeError::TransportDomain { quantity: found, .. }) if found == quantity)
        );
        assert_eq!(fields(&n), before);
        assert_eq!(wall.branch_contacts(0).unwrap(), contacts);
        assert_eq!(wall.last_report(), report);
        assert_eq!(wall.accepted_branch_impulses().unwrap(), impulses);
        assert_eq!(n.time(), Time::new(0.0));
    }
}

fn churchill_geometry(roughness: f64) -> CircularChurchillWall {
    CircularChurchillWall {
        diameter: Length::new(0.02),
        equivalent_sand_roughness: Length::new(roughness * 0.02),
        nusselt_number: 3.66,
        reynolds_number: [0.0, 1e8],
        maximum_mach_number: 0.3,
    }
}

// Independent original Darcy equation; only used at nonzero, moderate Re.
fn reference_churchill(re: f64, roughness: f64) -> f64 {
    let a = (-2.457 * ((7.0 / re).powf(0.9) + 0.27 * roughness).ln()).powi(16);
    let b = (37_530.0 / re).powi(16);
    8.0 * ((8.0 / re).powi(12) + (a + b).powf(-1.5)).powf(1.0 / 12.0)
}

#[test]
fn churchill_g1_native_braking_crosses_transition_and_refines_to_continuous_law() {
    let g = gas(1, false, 0.0);
    let (mu, duration, diameter) = (1e-6, 40.0, 0.02);
    let rhs = |u: f64| -reference_churchill(u * diameter / mu, 0.001) * u * u / (2.0 * diameter);
    let mut reference = 0.4;
    let h = duration / 8192.0;
    for _ in 0..8192 {
        let k1 = rhs(reference);
        let k2 = rhs(reference + h * k1 / 2.0);
        let k3 = rhs(reference + h * k2 / 2.0);
        let k4 = rhs(reference + h * k3);
        reference += h * (k1 + 2.0 * k2 + 2.0 * k3 + k4) / 6.0;
    }
    let mut errors = Vec::new();
    for steps in [32, 64, 128] {
        let mut n = thermal_pipe(&g, &[1.0], [0.4; 2]);
        let before = fields(&n);
        let volume = n.branches()[0].mesh().volumes()[0];
        let mut wall = uniform_operator(
            &n,
            WallContact {
                drag: WallDrag::None,
                thermal: WallThermal::Adiabatic,
            },
        )
        .with_circular_transport(
            law(&g, mu, 0.0, 0.0, 0.0),
            vec![Some(CircularWallModel::Churchill(churchill_geometry(
                0.001,
            )))],
            || false,
        )
        .unwrap();
        let (mut heat, mut impulse) = (0.0, 0.0);
        let mut regimes = [false; 3];
        for _ in 0..steps {
            let u = g
                .primitive(n.branches()[0].cell(0).unwrap(), || false)
                .unwrap()
                .velocity;
            let re = u * diameter / mu;
            regimes[if re >= 4000.0 {
                0
            } else if re >= 2100.0 {
                1
            } else {
                2
            }] = true;
            let report = wall
                .advance_for(
                    &mut n,
                    Time::new(duration / f64::from(steps)),
                    controls(),
                    || false,
                )
                .unwrap();
            heat += report.drag_heating_j;
            impulse += wall.accepted_branch_impulses().unwrap()[0];
            assert_eq!(report.heat_to_gas_j, 0.0);
            assert_eq!(report.mechanical_work_to_gas_j, 0.0);
        }
        assert!(regimes.into_iter().all(|visited| visited));
        let cell = n.branches()[0].cell(0).unwrap();
        let actual = g.primitive(cell, || false).unwrap().velocity;
        assert!(actual * diameter / mu < 1000.0);
        assert_eq!(cell.energy, before.0[0][0].energy);
        assert_eq!(cell.partial_densities, before.0[0][0].partial_densities);
        close(heat, volume * (0.4_f64.powi(2) - actual * actual), 2e-14);
        close(impulse, 2.0 * volume * (0.4 - actual), 2e-14);
        errors.push((actual - reference).abs());
    }
    assert!(errors[2] < 7.1e-5, "{errors:?}");
    for pair in errors.windows(2) {
        assert!((1.9..2.1).contains(&(pair[0] / pair[1])), "{errors:?}");
    }
}

#[test]
fn churchill_g3_live_material_roughness_and_slip_keep_one_energy_transaction() {
    let g = gas(1, false, 0.0);
    let (dt, conductivity, diameter) = (0.4, 7e-7, 0.02);
    for (density, mu, roughness, boost, slip) in [
        (1.0, 1e-6, 0.0, 0.0, 0.4),
        (1.0, 2e-6, 0.0, 0.0, 0.4),
        (2.0, 1e-6, 0.0, 0.0, 0.4),
        (1.0, 1e-6, 0.01, 0.0, 0.4),
        (1.0, 1e-6, 0.01, 20.0, 0.4),
        (1.0, 1e-6, 0.01, 20.0, -0.4),
    ] {
        let mut n = thermal_pipe(&g, &[density], [boost + slip; 2]);
        let before = fields(&n);
        let volume = n.branches()[0].mesh().volumes()[0];
        let model = churchill_geometry(roughness);
        let mut wall = uniform_operator(
            &n,
            WallContact {
                drag: WallDrag::MovingLinear {
                    rate_per_s: 999.0,
                    wall_velocity_m_s: boost,
                },
                ..thermal_contact()
            },
        )
        .with_circular_transport(
            law(&g, mu, 0.0, conductivity, 0.0),
            vec![Some(CircularWallModel::Churchill(model))],
            || false,
        )
        .unwrap();
        let re = density * slip.abs() * diameter / mu;
        let a = 32.0 * mu / (density * diameter * diameter);
        let b = (reference_churchill(re, roughness) - 64.0 / re) / (2.0 * diameter);
        let end_slip = slip * (-a * dt).exp() / (1.0 + b * slip.abs() * -(-a * dt).exp_m1() / a);
        let heated = 2.0 + (slip * slip - end_slip * end_slip) / 5.0;
        let gas_capacity = 2.5 * density * volume;
        let conductance = model.nusselt_number * conductivity * 4.0 * volume / diameter.powi(2);
        let q = dt * conductance * (2.2 - heated)
            / (1.0 + dt * conductance * (1.0 / gas_capacity + 1.0 / 5e-4));
        let report = wall
            .advance_for(&mut n, Time::new(dt), controls(), || false)
            .unwrap();
        for i in 0..2 {
            let cell = n.branches()[0].cell(i).unwrap();
            let primitive = g.primitive(cell, || false).unwrap();
            close(primitive.velocity, boost + end_slip, 3e-12);
            close(primitive.temperature, heated + q / gas_capacity, 2e-8);
            assert_eq!(cell.partial_densities, before.0[0][i].partial_densities);
            let contact = wall.branch_contacts(0).unwrap()[i];
            let WallDrag::LinearQuadratic {
                linear_rate_per_s,
                quadratic_per_m,
                wall_velocity_m_s,
            } = contact.drag
            else {
                panic!("Churchill discarded linear-quadratic wall mechanics")
            };
            close(linear_rate_per_s, a, 2e-13);
            close(quadratic_per_m, b, 2e-12);
            assert_eq!(wall_velocity_m_s, boost);
            close(wall_temperature(contact), 2.2 - q / 5e-4, 2e-8);
            close(
                (cell.energy - before.0[0][i].energy) * volume,
                q + density * volume * boost * (end_slip - slip),
                3e-12,
            );
        }
        close(report.heat_to_gas_j, 2.0 * q, 3e-12);
        close(report.wall_energy_change_j, -2.0 * q, 3e-12);
        close(
            report.drag_heating_j,
            density * volume * (slip * slip - end_slip * end_slip),
            3e-12,
        );
        close(
            report.mechanical_work_to_gas_j,
            2.0 * density * volume * boost * (end_slip - slip),
            3e-12,
        );
        close(report.energy_residual_j, 0.0, 3e-12);
        assert_eq!(report.external_heat_to_gas_j, 0.0);
    }
}

#[test]
fn churchill_g0_native_rest_laminar_limit_and_explicit_applicability() {
    let model = churchill_geometry(0.001);
    for re in [-1.0, 1e8 + 1.0, f64::NAN, f64::INFINITY] {
        assert!(model.darcy_reynolds(re).is_err());
    }
    for invalid in [
        CircularChurchillWall {
            diameter: Length::new(0.0),
            ..model
        },
        CircularChurchillWall {
            equivalent_sand_roughness: Length::new(0.051 * 0.02),
            ..model
        },
        CircularChurchillWall {
            nusselt_number: -1.0,
            ..model
        },
        CircularChurchillWall {
            reynolds_number: [4000.0, 1e8],
            ..model
        },
        CircularChurchillWall {
            maximum_mach_number: 0.31,
            ..model
        },
    ] {
        assert!(invalid.darcy_reynolds(0.0).is_err());
    }
    assert_eq!(model.darcy_reynolds(0.0).unwrap(), 64.0);
    let g = gas(1, false, 0.0);
    let mut n = thermal_pipe(&g, &[1.0], [0.0, 1e-7]);
    let before = fields(&n);
    let mut wall = uniform_operator(
        &n,
        WallContact {
            drag: WallDrag::None,
            thermal: WallThermal::Adiabatic,
        },
    )
    .with_circular_transport(
        law(&g, 1e-6, 0.0, 0.0, 0.0),
        vec![Some(CircularWallModel::Churchill(model))],
        || false,
    )
    .unwrap();
    wall.advance_for(&mut n, Time::new(40.0), controls(), || false)
        .unwrap();
    assert_eq!(fields(&n).0[0][0], before.0[0][0]);
    let velocity = g
        .primitive(n.branches()[0].cell(1).unwrap(), || false)
        .unwrap()
        .velocity;
    close(velocity, 1e-7 * (-3.2_f64).exp(), 1e-20);
}

#[test]
fn churchill_g4_late_material_refusal_preserves_all_cells_walls_and_reports() {
    let g = gas(1, false, 0.0);
    let mut n = thermal_pipe(&g, &[1.0], [0.2; 2]);
    let original = law(&g, 1e-6, 0.0, 1.0, 0.0);
    let limited = ConstantTransport::new(
        g.source().clone(),
        "bounded-Churchill-properties",
        1,
        TransportDomain {
            minimum_temperature_k: 0.001,
            maximum_temperature_k: 2.05,
            minimum_pressure_pa: 1e-8,
            maximum_pressure_pa: 1e8,
        },
        original.coefficients(),
        vec![0.0],
    )
    .unwrap();
    let mut wall = NetworkWallOperator::new(
        &n,
        vec![vec![
            WallContact {
                drag: WallDrag::None,
                thermal: WallThermal::Adiabatic,
            },
            thermal_contact(),
        ]],
        Vec::new(),
        || false,
    )
    .unwrap()
    .with_circular_transport(
        limited,
        vec![Some(CircularWallModel::Churchill(churchill_geometry(
            0.001,
        )))],
        || false,
    )
    .unwrap();
    wall.advance_for(&mut n, Time::new(0.0), controls(), || false)
        .unwrap();
    let before = fields(&n);
    let contacts = wall.branch_contacts(0).unwrap().to_vec();
    let report = wall.last_report();
    let impulses = wall.accepted_branch_impulses().unwrap().to_vec();
    assert!(matches!(
        wall.advance_for(&mut n, Time::new(0.4), controls(), || false),
        Err(WallExchangeError::Transport(_))
    ));
    assert_eq!(fields(&n), before);
    assert_eq!(wall.branch_contacts(0).unwrap(), contacts);
    assert_eq!(wall.last_report(), report);
    assert_eq!(wall.accepted_branch_impulses().unwrap(), impulses);
    assert_eq!(n.time(), Time::new(0.0));
}
