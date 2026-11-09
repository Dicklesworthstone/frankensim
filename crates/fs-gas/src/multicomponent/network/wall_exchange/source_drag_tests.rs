//! G1 original-source thermoviscous relaxation and G4 native publication.

mod recoil;

use super::*;
use crate::multicomponent::dissipation::{
    ConstantTransport, DissipationError, TransportCoefficients, TransportLaw,
};
use crate::multicomponent::network::dissipation::tests::{close, fields, gas, law};
use crate::multicomponent::wall_exchange::{SharedRecoilControl, SlidingSupport};
use crate::multicomponent::{Boundaries, MixtureEuler1d};
use fs_qty::{DynViscosity, Length, Pressure, Temperature};
use fs_thermochem::equilibrium::composition::IdealGasEquilibriumV1;
use fs_thermochem::transport::kinetic::ThermalConductivityQuantityV1;

// The retained original fixture has R=1 and Cv=2.5 J/(kg K).
struct TemperatureLaw {
    original: ConstantTransport,
    viscosity_slope: f64,
    conductivity_slope: f64,
    maximum_temperature: f64,
}

impl TransportLaw for TemperatureLaw {
    fn source(&self) -> &IdealGasEquilibriumV1 {
        self.original.source()
    }

    fn evaluate(
        &mut self,
        temperature: Temperature,
        pressure: Pressure,
        mole_fractions: &[f64],
        binary_diffusivities: &mut [f64],
        cancelled: &mut dyn FnMut() -> bool,
    ) -> Result<TransportCoefficients, DissipationError> {
        if temperature.value() > self.maximum_temperature {
            return Err(DissipationError::Invalid(
                "manufactured wall temperature domain",
            ));
        }
        let mut coefficients = self.original.evaluate(
            temperature,
            pressure,
            mole_fractions,
            binary_diffusivities,
            cancelled,
        )?;
        let rise = temperature.value() - 2.0;
        coefficients.dynamic_viscosity = DynViscosity::new(
            coefficients.dynamic_viscosity.value() * (1.0 + self.viscosity_slope * rise),
        );
        coefficients.thermal_conductivity = ThermalConductivityQuantityV1::new(
            coefficients.thermal_conductivity.value() * (1.0 + self.conductivity_slope * rise),
        );
        Ok(coefficients)
    }
}

fn source_controls() -> SourceCoupledDragControl {
    SourceCoupledDragControl {
        absolute_slip_tolerance_m_s: 1e-10,
        relative_slip_tolerance: 1e-8,
        absolute_heating_tolerance_j: 1e-12,
        relative_heating_tolerance: 1e-8,
        maximum_attempted_steps: 16384,
        maximum_source_evaluations: 131072,
    }
}

fn geometry() -> CircularLaminarWall {
    CircularLaminarWall {
        diameter: Length::new((0.4 / core::f64::consts::FRAC_PI_4).sqrt()),
        nusselt_number: 3.66,
        reynolds_number: [0.0, 2300.0],
        maximum_mach_number: 0.3,
    }
}

fn network(gas: &MixtureGas, velocities: &[f64]) -> MixtureNetwork1d {
    let branches = velocities
        .iter()
        .map(|&velocity| {
            MixtureEuler1d::new(
                gas.clone(),
                Length::new(0.5),
                0.4,
                vec![gas.conserved(&[1.0], velocity, 2.0).unwrap(); 2],
                Boundaries::Periodic,
                || false,
            )
            .unwrap()
        })
        .collect();
    MixtureNetwork1d::new(branches, Vec::new(), Vec::new(), || false).unwrap()
}

fn original_operator(
    network: &MixtureNetwork1d,
    contacts: Vec<Vec<WallContact>>,
    viscosity_slope: f64,
    conductivity_slope: f64,
    maximum_temperature: f64,
    branches: Vec<Option<CircularLaminarWall>>,
) -> NetworkWallOperator<'static> {
    NetworkWallOperator::new(network, contacts, Vec::new(), || false)
        .unwrap()
        .with_circular_laminar_transport(
            TemperatureLaw {
                original: law(network.gas(), 0.004, 0.0, 0.02, 0.0),
                viscosity_slope,
                conductivity_slope,
                maximum_temperature,
            },
            branches,
            || false,
        )
        .unwrap()
}

fn adiabatic() -> WallContact {
    WallContact {
        drag: WallDrag::None,
        thermal: WallThermal::Adiabatic,
    }
}

fn lumped() -> WallContact {
    WallContact {
        drag: WallDrag::None,
        thermal: WallThermal::Lumped {
            temperature_k: 5.0,
            heat_capacity_j_per_k: 1.0,
            conductance_w_per_k: 0.0,
        },
    }
}

#[test]
fn source_coupled_native_drag_follows_temperature_viscosity_and_exact_heat() {
    let gas = gas(1, false, 0.0);
    let mut network = network(&gas, &[0.4]);
    let before = fields(&network);
    let mut wall = original_operator(
        &network,
        vec![vec![adiabatic(); 2]],
        100.0,
        0.0,
        50.0,
        vec![Some(geometry())],
    )
    .with_source_coupled_drag(source_controls(), || false)
    .unwrap();
    let duration = 1.0;
    let report = wall
        .advance_for(
            &mut network,
            Time::new(duration),
            WallControl::default(),
            || false,
        )
        .unwrap();
    // With T=2+(u0²-u²)/5 and mu=mu0[1+100(T-2)], y=u² solves
    // y'=-2a(C-20y)y. This is an independent closed-form continuum solution.
    let a = 32.0 * 0.004 / geometry().diameter.value().powi(2);
    let initial_square = 0.16;
    let c = 1.0 + 20.0 * initial_square;
    let decay = (-2.0 * a * c * duration).exp();
    let expected_square = c * initial_square * decay / (c - 20.0 * initial_square * (1.0 - decay));
    let expected_speed = expected_square.sqrt();
    let frozen_speed = 0.4 * (-a * duration).exp();
    let mut total_loss = 0.0;
    for cell in 0..2 {
        let state = network.branches()[0].cell(cell).unwrap();
        let primitive = gas.primitive(state, || false).unwrap();
        close(primitive.velocity, expected_speed, 2e-7);
        assert!((primitive.velocity - frozen_speed).abs() > 0.05);
        close(
            primitive.temperature,
            2.0 + (initial_square - expected_square) / 5.0,
            1e-7,
        );
        assert_eq!(state.partial_densities, before.0[0][cell].partial_densities);
        assert_eq!(state.energy, before.0[0][cell].energy);
        total_loss += 0.1 * (initial_square - primitive.velocity.powi(2));
    }
    close(report.drag_heating_j, total_loss, 1e-12);
    assert_eq!(report.mechanical_work_to_gas_j, 0.0);
    assert_eq!(report.heat_to_gas_j, 0.0);
    assert_eq!(report.energy_residual_j, 0.0);
    assert_eq!(report.maximum_iterations, 0);
    assert_eq!(network.time(), Time::new(0.0));
}

#[test]
fn source_coupled_native_moving_wall_retains_work_slip_heat_and_boost_covariance() {
    let gas = gas(1, false, 0.0);
    let duration = 0.7;
    let rate = 32.0 * 0.004 / geometry().diameter.value().powi(2);
    let slip = 0.2 * (-rate * duration).exp();
    for boost in [0.0, 20.0] {
        let wall_speed = 0.2 + boost;
        let mut network = network(&gas, &[0.4 + boost]);
        let contact = WallContact {
            drag: WallDrag::MovingLinear {
                rate_per_s: 99.0,
                wall_velocity_m_s: wall_speed,
            },
            thermal: WallThermal::Adiabatic,
        };
        let mut wall = original_operator(
            &network,
            vec![vec![contact; 2]],
            0.0,
            0.0,
            50.0,
            vec![Some(geometry())],
        )
        .with_source_coupled_drag(source_controls(), || false)
        .unwrap();
        let report = wall
            .advance_for(
                &mut network,
                Time::new(duration),
                WallControl::default(),
                || false,
            )
            .unwrap();
        for cell in 0..2 {
            let primitive = gas
                .primitive(network.branches()[0].cell(cell).unwrap(), || false)
                .unwrap();
            close(primitive.velocity, wall_speed + slip, 2e-12);
            close(
                primitive.temperature,
                2.0 + (0.04 - slip * slip) / 5.0,
                2e-10,
            );
            assert_eq!(
                wall.branch_contacts(0).unwrap()[cell]
                    .drag
                    .moving_velocity(),
                Some(wall_speed)
            );
        }
        close(
            report.mechanical_work_to_gas_j,
            0.4 * wall_speed * (slip - 0.2),
            1e-11,
        );
        close(report.drag_heating_j, 0.2 * (0.04 - slip * slip), 1e-12);
        close(report.energy_residual_j, 0.0, 1e-11);
    }
}

#[test]
fn mechanical_endpoint_refreshes_original_local_film_before_finite_wall_heat() {
    let gas = gas(1, false, 0.0);
    let mut network = network(&gas, &[0.4]);
    let mut wall = original_operator(
        &network,
        vec![vec![lumped(); 2]],
        0.0,
        100.0,
        50.0,
        vec![Some(geometry())],
    )
    .with_source_coupled_drag(source_controls(), || false)
    .unwrap();
    let duration = 1.0;
    let diameter = geometry().diameter.value();
    let speed = 0.4 * (-32.0 * 0.004 * duration / diameter.powi(2)).exp();
    let heated_temperature = 2.0 + (0.16 - speed * speed) / 5.0;
    let conductivity = 0.02 * (1.0 + 100.0 * (heated_temperature - 2.0));
    let film = 3.66 * conductivity * (4.0 * 0.2 / diameter) / diameter;
    let heat = duration * film * (5.0 - heated_temperature) / (1.0 + 3.0 * duration * film);
    let report = wall
        .advance_for(
            &mut network,
            Time::new(duration),
            WallControl::default(),
            || false,
        )
        .unwrap();
    for cell in 0..2 {
        let primitive = gas
            .primitive(network.branches()[0].cell(cell).unwrap(), || false)
            .unwrap();
        close(primitive.velocity, speed, 2e-12);
        close(primitive.temperature, heated_temperature + heat / 0.5, 2e-9);
        let WallThermal::Lumped {
            temperature_k,
            conductance_w_per_k,
            ..
        } = wall.branch_contacts(0).unwrap()[cell].thermal
        else {
            panic!("lost original finite wall");
        };
        close(temperature_k, 5.0 - heat, 2e-9);
        close(conductance_w_per_k, film, 2e-11);
    }
    close(report.heat_to_gas_j, 2.0 * heat, 2e-9);
    close(report.wall_energy_change_j, -2.0 * heat, 2e-9);
    close(report.energy_residual_j, 0.0, 2e-10);
}

#[test]
fn source_refusal_and_final_cancellation_preserve_complete_native_transaction() {
    let gas = gas(1, false, 0.0);
    let prepare = |maximum_temperature| {
        let mut network = network(&gas, &[0.4]);
        let mut wall = original_operator(
            &network,
            vec![vec![adiabatic(), lumped()]],
            100.0,
            0.0,
            maximum_temperature,
            vec![Some(geometry())],
        )
        .with_source_coupled_drag(source_controls(), || false)
        .unwrap();
        wall.advance_for(&mut network, Time::new(0.0), WallControl::default(), || {
            false
        })
        .unwrap();
        (network, wall)
    };
    let (mut network, mut wall) = prepare(2.2);
    let before = fields(&network);
    let contacts = wall.branch_contacts(0).unwrap().to_vec();
    let report = wall.last_report();
    assert!(matches!(
        wall.advance_for(&mut network, Time::new(1.0), WallControl::default(), || {
            false
        }),
        Err(WallExchangeError::Transport(_))
    ));
    assert_eq!(fields(&network), before);
    assert_eq!(wall.branch_contacts(0).unwrap(), contacts);
    assert_eq!(wall.last_report(), report);

    let (mut probe_network, mut probe_wall) = prepare(50.0);
    let mut checkpoints = 0;
    probe_wall
        .advance_for(
            &mut probe_network,
            Time::new(1.0),
            WallControl::default(),
            || {
                checkpoints += 1;
                false
            },
        )
        .unwrap();
    let (mut network, mut wall) = prepare(50.0);
    let before = fields(&network);
    let contacts = wall.branch_contacts(0).unwrap().to_vec();
    let report = wall.last_report();
    let mut reached = 0;
    assert!(matches!(
        wall.advance_for(&mut network, Time::new(1.0), WallControl::default(), || {
            reached += 1;
            reached == checkpoints
        }),
        Err(WallExchangeError::Cancelled)
    ));
    assert_eq!(fields(&network), before);
    assert_eq!(wall.branch_contacts(0).unwrap(), contacts);
    assert_eq!(wall.last_report(), report);
    wall.advance_for(&mut network, Time::new(1.0), WallControl::default(), || {
        false
    })
    .unwrap();
    assert_eq!(fields(&network), fields(&probe_network));
    assert_eq!(wall.last_report(), probe_wall.last_report());
}

#[test]
fn independent_supports_bind_both_orders_and_shared_supports_refuse_selected_branches() {
    let gas = gas(1, false, 0.0);
    let network = network(&gas, &[0.4, 0.4]);
    let make = || {
        original_operator(
            &network,
            vec![vec![adiabatic(); 2]; 2],
            0.0,
            0.0,
            50.0,
            vec![Some(geometry()), None],
        )
    };
    let support = SlidingSupport {
        mass_kg: 1.0,
        velocity_m_s: 0.0,
    };
    let selected = || vec![vec![Some(support); 2], vec![None; 2]];
    assert!(
        make()
            .with_finite_supports(selected(), 1e-10, || false)
            .unwrap()
            .with_source_coupled_drag(source_controls(), || false)
            .is_ok()
    );
    assert!(
        make()
            .with_source_coupled_drag(source_controls(), || false)
            .unwrap()
            .with_finite_supports(selected(), 1e-10, || false)
            .is_ok()
    );
    assert!(
        make()
            .with_shared_branch_supports(
                vec![Some(support), None],
                SharedRecoilControl::default(),
                1e-10,
                || false,
            )
            .unwrap()
            .with_source_coupled_drag(source_controls(), || false)
            .is_err()
    );
    assert!(
        make()
            .with_source_coupled_drag(source_controls(), || false)
            .unwrap()
            .with_shared_branch_supports(
                vec![Some(support), None],
                SharedRecoilControl::default(),
                1e-10,
                || false,
            )
            .is_err()
    );
    assert!(
        make()
            .with_source_coupled_drag(source_controls(), || false)
            .unwrap()
            .with_finite_supports(vec![vec![None; 2], vec![Some(support); 2]], 1e-10, || false)
            .is_ok()
    );
    assert!(
        make()
            .with_shared_branch_supports(
                vec![None, Some(support)],
                SharedRecoilControl::default(),
                1e-10,
                || false,
            )
            .unwrap()
            .with_source_coupled_drag(source_controls(), || false)
            .is_ok()
    );
}
