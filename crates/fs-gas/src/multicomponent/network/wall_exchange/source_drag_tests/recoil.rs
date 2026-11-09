//! Independent finite-support source relaxation through the native transaction.

use super::*;
use crate::multicomponent::wall_exchange::MaterialWallCaloric;
use fs_thermal::enthalpy::{
    Boundary, CellGeometry, CellId, DepositedCell, DepositionError, MaterialId, ThermalLaw,
    ThermalStorageBasis,
};

const GAS_MASS: f64 = 0.2;
const SUPPORT_MASS: f64 = 0.3;

fn support(boost: f64) -> SlidingSupport {
    SlidingSupport {
        mass_kg: SUPPORT_MASS,
        velocity_m_s: boost,
    }
}

fn no_heat_geometry() -> CircularLaminarWall {
    CircularLaminarWall {
        nusselt_number: 0.0,
        ..geometry()
    }
}

fn bind_supports(
    wall: NetworkWallOperator<'static>,
    boost: f64,
    fractions: [f64; 2],
) -> NetworkWallOperator<'static> {
    wall.with_finite_supports(vec![vec![Some(support(boost)); 2]], 1e-10, || false)
        .unwrap()
        .with_finite_support_heat_partition(vec![fractions.to_vec()], || false)
        .unwrap()
        .with_source_coupled_drag(source_controls(), || false)
        .unwrap()
}

// Return gas speed, support speed and true relative-slip energy loss. The
// two-body center-of-mass solution is independent of production mechanics.
fn exact_constant(duration: f64, boost: f64) -> (f64, f64, f64) {
    let rate = 32.0 * 0.004 / geometry().diameter.value().powi(2);
    let kappa = 1.0 + GAS_MASS / SUPPORT_MASS;
    let slip = 0.4 * (-kappa * rate * duration).exp();
    let center = boost + GAS_MASS * 0.4 / (GAS_MASS + SUPPORT_MASS);
    let gas_speed = center + SUPPORT_MASS * slip / (GAS_MASS + SUPPORT_MASS);
    let support_speed = center - GAS_MASS * slip / (GAS_MASS + SUPPORT_MASS);
    let reduced = GAS_MASS * SUPPORT_MASS / (GAS_MASS + SUPPORT_MASS);
    let loss = 0.5 * reduced * (0.16 - slip * slip);
    (gas_speed, support_speed, loss)
}

#[test]
fn constant_source_recoil_retains_exact_two_body_partition_and_boost_covariance() {
    let gas = gas(1, false, 0.0);
    let duration = 0.7;
    for boost in [0.0, 20.0] {
        for fraction in [0.0, 0.35, 1.0] {
            let mut network = network(&gas, &[boost + 0.4]);
            let before = fields(&network);
            let mut wall = bind_supports(
                original_operator(
                    &network,
                    vec![vec![lumped(); 2]],
                    0.0,
                    0.0,
                    50.0,
                    vec![Some(no_heat_geometry())],
                ),
                boost,
                [fraction; 2],
            );
            let report = wall
                .advance_for(
                    &mut network,
                    Time::new(duration),
                    WallControl::default(),
                    || false,
                )
                .unwrap();
            let (gas_speed, support_speed, loss) = exact_constant(duration, boost);
            let wall_heat = (1.0 - fraction) * loss;
            let kinetic = 0.5 * SUPPORT_MASS * (support_speed - boost) * (support_speed + boost);
            let mut gas_energy_change = 0.0;
            for cell in 0..2 {
                let state = network.branches()[0].cell(cell).unwrap();
                let primitive = gas.primitive(state, || false).unwrap();
                close(primitive.velocity, gas_speed, 2e-11);
                close(primitive.temperature, 2.0 + fraction * loss / 0.5, 2e-9);
                let actual_support = wall.branch_supports(0).unwrap()[cell].unwrap();
                close(actual_support.velocity_m_s, support_speed, 2e-11);
                assert_eq!(actual_support.mass_kg, SUPPORT_MASS);
                let WallThermal::Lumped { temperature_k, .. } =
                    wall.branch_contacts(0).unwrap()[cell].thermal
                else {
                    panic!("lost original finite wall");
                };
                close(temperature_k, 5.0 + wall_heat, 2e-10);
                let exchange = wall.accepted_branch_support_exchanges(0).unwrap()[cell].unwrap();
                close(exchange.kinetic_energy_change_j, kinetic, 2e-10);
                close(exchange.drag_heating_j, loss, 2e-12);
                close(exchange.wall_slip_heating_j, wall_heat, 2e-12);
                close(exchange.momentum_residual_ns, 0.0, 2e-12);
                close(exchange.energy_residual_j, 0.0, 2e-10);
                assert_eq!(state.partial_densities, before.0[0][cell].partial_densities);
                gas_energy_change += 0.2 * (state.energy - before.0[0][cell].energy);
            }
            close(report.drag_heating_j, 2.0 * loss, 2e-12);
            close(report.mechanical_work_to_gas_j, -2.0 * kinetic, 2e-10);
            close(report.heat_to_gas_j, -2.0 * wall_heat, 2e-12);
            close(report.wall_energy_change_j, 2.0 * wall_heat, 2e-12);
            close(gas_energy_change + 2.0 * (kinetic + wall_heat), 0.0, 2e-10);
            close(report.energy_residual_j, 0.0, 2e-10);
        }
    }
}

// Original nonlinear wall u(T)=0.2*T² J/kg. Its 0.1 kg amount belongs
// to the authored 0.3 kg support; the binding adds no mechanical mass.
struct QuadraticWall;
impl ThermalLaw for QuadraticWall {
    fn storage_basis(&self) -> ThermalStorageBasis {
        ThermalStorageBasis::InternalEnergy
    }
    fn temperature_k(&self, energy: f64) -> Result<f64, DepositionError> {
        if !energy.is_finite() || energy <= 0.0 {
            return Err(DepositionError::MaterialDomain);
        }
        Ok((5.0 * energy).sqrt())
    }
    fn maximum_temperature_slope(&self) -> f64 {
        100.0
    }
}

fn original_wall() -> MaterialWallCaloric<'static> {
    MaterialWallCaloric::from_thermal_law(
        QuadraticWall,
        DepositedCell {
            id: CellId(17),
            segment_id: 3,
            body_id: 8,
            material_id: MaterialId(29),
            mass_kg: 0.1,
            specific_h_j_kg: 5.0,
            geometry: CellGeometry {
                start_m: [0.0; 3],
                end_m: [1.0, 0.0, 0.0],
                volume_m3: 0.01,
                surface_area_m2: 1.0,
            },
            boundary: Boundary {
                ambient_k: 2.0,
                htc_w_m2_k: 0.0,
                reservoir: None,
                source_w: 0.0,
            },
        },
        [0.01, 50.0],
    )
    .unwrap()
}

#[test]
fn current_gas_viscosity_recoils_support_and_deposits_original_material_slip_share() {
    let gas = gas(1, false, 0.0);
    let mut network = network(&gas, &[0.4]);
    let material = original_wall();
    let contact = material.initial_contact(0.0, WallDrag::None).unwrap();
    let fraction = 0.35;
    let mut wall = bind_supports(
        original_operator(
            &network,
            vec![vec![contact; 2]],
            100.0,
            0.0,
            50.0,
            vec![Some(no_heat_geometry())],
        )
        .with_material_heat_storage(
            vec![vec![Some(material.clone()), Some(material)]],
            Vec::new(),
            || false,
        )
        .unwrap(),
        0.0,
        [fraction; 2],
    );
    let duration = 1.0;
    let report = wall
        .advance_for(
            &mut network,
            Time::new(duration),
            WallControl::default(),
            || false,
        )
        .unwrap();
    let reduced = GAS_MASS * SUPPORT_MASS / (GAS_MASS + SUPPORT_MASS);
    let kappa = GAS_MASS / reduced;
    let rate = 32.0 * 0.004 / geometry().diameter.value().powi(2);
    // Tgas=2+f*mu_reduced*(w0²-w²)/(2*m*Cv). Consequently
    // (w²)'=-2*kappa*a*(C-B*w²)*w², again a closed-form solution.
    let b = 100.0 * fraction * reduced / (5.0 * GAS_MASS);
    let c = 1.0 + b * 0.16;
    let decay = (-2.0 * kappa * rate * c * duration).exp();
    let square = c * 0.16 * decay / (c - b * 0.16 * (1.0 - decay));
    let slip = square.sqrt();
    let gas_speed = 0.16 + 0.6 * slip;
    let support_speed = 0.16 - 0.4 * slip;
    let loss = 0.5 * reduced * (0.16 - square);
    let wall_heat = (1.0 - fraction) * loss;
    for cell in 0..2 {
        let primitive = gas
            .primitive(network.branches()[0].cell(cell).unwrap(), || false)
            .unwrap();
        close(primitive.velocity, gas_speed, 2e-7);
        close(primitive.temperature, 2.0 + fraction * loss / 0.5, 2e-8);
        close(
            wall.branch_supports(0).unwrap()[cell].unwrap().velocity_m_s,
            support_speed,
            2e-7,
        );
        let actual = wall.branch_contacts(0).unwrap()[cell];
        let WallThermal::Caloric {
            mass_kg,
            specific_energy_j_kg,
            ..
        } = actual.thermal
        else {
            panic!("lost original wall energy coordinate");
        };
        assert_eq!(mass_kg, 0.1);
        close(mass_kg * (specific_energy_j_kg - 5.0), wall_heat, 2e-9);
        let source = wall.branch_heat_storage(0).unwrap()[cell].as_ref().unwrap();
        close(
            source.temperature_k(actual.thermal).unwrap(),
            (25.0 + 50.0 * wall_heat).sqrt(),
            2e-8,
        );
        let exchange = wall.accepted_branch_support_exchanges(0).unwrap()[cell].unwrap();
        close(exchange.wall_slip_heating_j, wall_heat, 2e-9);
        close(exchange.energy_residual_j, 0.0, 2e-10);
    }
    close(report.drag_heating_j, 2.0 * loss, 2e-9);
    close(report.wall_energy_change_j, 2.0 * wall_heat, 2e-9);
    close(report.heat_to_gas_j, -2.0 * wall_heat, 2e-9);
    close(report.energy_residual_j, 0.0, 2e-10);
}

#[test]
fn finite_recoil_refreshes_local_film_after_gas_and_wall_friction_partition() {
    let gas = gas(1, false, 0.0);
    let mut network = network(&gas, &[0.4]);
    let fraction = 0.35;
    let duration = 1.0;
    let mut wall = bind_supports(
        original_operator(
            &network,
            vec![vec![lumped(); 2]],
            0.0,
            100.0,
            50.0,
            vec![Some(geometry())],
        ),
        0.0,
        [fraction; 2],
    );
    let (gas_speed, support_speed, loss) = exact_constant(duration, 0.0);
    let wall_heat = (1.0 - fraction) * loss;
    let gas_temperature = 2.0 + fraction * loss / 0.5;
    let wall_temperature = 5.0 + wall_heat;
    let conductivity = 0.02 * (1.0 + 100.0 * (gas_temperature - 2.0));
    let diameter = geometry().diameter.value();
    let film = 3.66 * conductivity * (4.0 * 0.2 / diameter) / diameter;
    let heat =
        duration * film * (wall_temperature - gas_temperature) / (1.0 + 3.0 * duration * film);
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
        close(primitive.velocity, gas_speed, 2e-12);
        close(primitive.temperature, gas_temperature + heat / 0.5, 2e-9);
        close(
            wall.branch_supports(0).unwrap()[cell].unwrap().velocity_m_s,
            support_speed,
            2e-12,
        );
        let WallThermal::Lumped {
            temperature_k,
            conductance_w_per_k,
            ..
        } = wall.branch_contacts(0).unwrap()[cell].thermal
        else {
            panic!("lost original finite wall");
        };
        close(temperature_k, wall_temperature - heat, 2e-9);
        close(conductance_w_per_k, film, 2e-11);
    }
    close(report.heat_to_gas_j, 2.0 * (heat - wall_heat), 2e-9);
    close(report.wall_energy_change_j, 2.0 * (wall_heat - heat), 2e-9);
    close(report.drag_heating_j, 2.0 * loss, 2e-12);
    close(report.energy_residual_j, 0.0, 2e-10);
}

#[test]
fn source_failure_and_final_cancellation_roll_back_support_heat_and_complete_receipts() {
    let gas = gas(1, false, 0.0);
    let prepare = |maximum_temperature| {
        let mut network = network(&gas, &[0.4]);
        let mut wall = bind_supports(
            original_operator(
                &network,
                vec![vec![adiabatic(), lumped()]],
                20.0,
                0.0,
                maximum_temperature,
                vec![Some(geometry())],
            ),
            0.0,
            [1.0, 0.35],
        );
        wall.advance_for(&mut network, Time::new(0.0), WallControl::default(), || {
            false
        })
        .unwrap();
        (network, wall)
    };
    let (mut network, mut wall) = prepare(2.2);
    let before = fields(&network);
    let contacts = wall.branch_contacts(0).unwrap().to_vec();
    let supports = wall.branch_supports(0).unwrap().to_vec();
    let exchanges = wall.accepted_branch_support_exchanges(0).unwrap().to_vec();
    let report = wall.last_report();
    assert!(matches!(
        wall.advance_for(&mut network, Time::new(1.0), WallControl::default(), || {
            false
        }),
        Err(WallExchangeError::Transport(_))
    ));
    assert_eq!(fields(&network), before);
    assert_eq!(wall.branch_contacts(0).unwrap(), contacts);
    assert_eq!(wall.branch_supports(0).unwrap(), supports);
    assert_eq!(
        wall.accepted_branch_support_exchanges(0).unwrap(),
        exchanges
    );
    assert_eq!(wall.last_report(), report);

    let (mut reference_network, mut reference_wall) = prepare(50.0);
    let mut checkpoints = 0;
    reference_wall
        .advance_for(
            &mut reference_network,
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
    let supports = wall.branch_supports(0).unwrap().to_vec();
    let exchanges = wall.accepted_branch_support_exchanges(0).unwrap().to_vec();
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
    assert_eq!(wall.branch_supports(0).unwrap(), supports);
    assert_eq!(
        wall.accepted_branch_support_exchanges(0).unwrap(),
        exchanges
    );
    assert_eq!(wall.last_report(), report);
    wall.advance_for(&mut network, Time::new(1.0), WallControl::default(), || {
        false
    })
    .unwrap();
    assert_eq!(fields(&network), fields(&reference_network));
    assert_eq!(
        wall.branch_supports(0).unwrap(),
        reference_wall.branch_supports(0).unwrap()
    );
    assert_eq!(
        wall.accepted_branch_support_exchanges(0).unwrap(),
        reference_wall.accepted_branch_support_exchanges(0).unwrap()
    );
    assert_eq!(wall.last_report(), reference_wall.last_report());
}
