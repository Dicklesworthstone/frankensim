//! Independent wall-work and thermal-equilibrium checks on resolved gas fields.

use super::*;
use crate::dynamics::hydrostatic::{
    AtmosphereQuadrature, HydrostaticReference, IsentropicAtmosphere, IsothermalAtmosphere,
};
use crate::dynamics::network::{BranchEnd, NetworkGravity, Plenum, PlenumConnection};
use crate::dynamics::{Boundaries, Boundary, Reconstruction, StepControl};
use fs_qty::{Acceleration, Length};
use fs_thermal::enthalpy::{
    Boundary as ThermalBoundary, CellGeometry, CellId, DepositedCell, DepositionError, MaterialId,
    ThermalLaw, ThermalStorageBasis,
};

const WALLS: Boundaries = Boundaries::Ends {
    left: Boundary::Reflecting,
    right: Boundary::Reflecting,
};

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "actual={actual:.16e}, expected={expected:.16e}, tolerance={tolerance:.3e}"
    );
}

fn gas(formation: f64) -> FrozenGas {
    let original = crate::dynamics::tests::gas(0.4, formation);
    let mut domain = original.domain();
    domain.absolute_energy_tolerance = 1e-13;
    domain.relative_energy_tolerance = 1e-14;
    FrozenGas::new(original.source().clone(), domain, || false).unwrap()
}

fn controls() -> WallControl {
    WallControl {
        absolute_energy_tolerance_j: 1e-11,
        relative_energy_tolerance: 1e-12,
        ..WallControl::default()
    }
}

fn contact(drag: WallDrag, thermal: WallThermal) -> BranchWallContact<'static> {
    BranchWallContact {
        contact: WallContact { drag, thermal },
        heat_storage: None,
    }
}

fn plenum_contact(thermal: WallThermal) -> PlenumThermalContact<'static> {
    PlenumThermalContact {
        thermal,
        heat_storage: None,
    }
}

fn periodic(velocities: &[f64]) -> GasNetwork1d {
    let gas = gas(0.0);
    let branches = velocities
        .iter()
        .map(|&velocity| {
            let state = gas.conserved(1.0, velocity, 2.0).unwrap();
            Euler1d::new(
                gas.clone(),
                Length::new(0.25),
                0.7,
                vec![state; 4],
                Boundaries::Periodic,
                || false,
            )
            .unwrap()
        })
        .collect();
    GasNetwork1d::new(branches, vec![], vec![], || false).unwrap()
}

fn atmosphere(
    gas: &FrozenGas,
    isentropic: bool,
    acceleration: f64,
    datum: f64,
    pressure: f64,
) -> Euler1d {
    let mesh = Mesh1d::new(
        vec![0.0, 0.1, 0.3, 0.6, 1.0],
        vec![0.7; 5],
        vec![0.07, 0.14, 0.21, 0.28],
    )
    .unwrap();
    if isentropic {
        Euler1d::new_isentropic_atmosphere(
            gas.clone(),
            mesh,
            IsentropicAtmosphere {
                acceleration: Acceleration::new(acceleration),
                reference_temperature: Temperature::new(2.0),
                reference_pressure: Pressure::new(pressure),
                reference_position: Length::new(datum),
            },
            AtmosphereQuadrature::default(),
            WALLS,
            || false,
        )
        .unwrap()
    } else {
        Euler1d::new_isothermal_atmosphere(
            gas.clone(),
            mesh,
            IsothermalAtmosphere {
                acceleration: Acceleration::new(acceleration),
                temperature: Temperature::new(2.0),
                reference_pressure: Pressure::new(pressure),
                reference_position: Length::new(datum),
            },
            WALLS,
            || false,
        )
        .unwrap()
    }
}

fn chain(datum: f64, formation: f64) -> GasNetwork1d {
    let gas = gas(formation);
    let lower = atmosphere(&gas, false, -0.5, 0.0, 3.0);
    let bottom = lower.hydrostatic_reference_face(0).unwrap();
    let middle = lower.hydrostatic_reference_face(4).unwrap();
    let upper = atmosphere(&gas, false, -0.5, 0.0, middle.1.pressure);
    let top = upper.hydrostatic_reference_face(4).unwrap();
    GasNetwork1d::new_with_gravity(
        vec![lower, upper],
        vec![
            Plenum::new(0.3, bottom.0).unwrap(),
            Plenum::new(0.3, middle.0).unwrap(),
            Plenum::new(0.3, top.0).unwrap(),
        ],
        vec![
            PlenumConnection {
                branch: 0,
                end: BranchEnd::Left,
                plenum: 0,
            },
            PlenumConnection {
                branch: 0,
                end: BranchEnd::Right,
                plenum: 1,
            },
            PlenumConnection {
                branch: 1,
                end: BranchEnd::Left,
                plenum: 1,
            },
            PlenumConnection {
                branch: 1,
                end: BranchEnd::Right,
                plenum: 2,
            },
        ],
        NetworkGravity {
            branch_offsets_j_kg: vec![datum, datum + 0.5],
            plenum_potentials_j_kg: vec![datum, datum + 0.5, datum + 1.0],
        },
        || false,
    )
    .unwrap()
}

fn driven(network: &GasNetwork1d) -> NetworkWallOperator<'static> {
    let branches = network
        .branches()
        .iter()
        .enumerate()
        .map(|(b, branch)| {
            Some(
                branch
                    .mesh()
                    .volumes()
                    .iter()
                    .map(|&volume| {
                        contact(
                            WallDrag::Linear { rate_per_s: 0.5 },
                            WallThermal::Lumped {
                                temperature_k: if b == 0 { 4.0 } else { 1.0 },
                                heat_capacity_j_per_k: 0.8 * volume,
                                conductance_w_per_k: 0.3 * volume,
                            },
                        )
                    })
                    .collect(),
            )
        })
        .collect();
    NetworkWallOperator::new(
        network,
        vec![plenum_contact(WallThermal::Adiabatic); network.plenums().len()],
        || false,
    )
    .unwrap()
    .with_branch_contacts(network, branches, || false)
    .unwrap()
}

fn wall_energy(walls: &NetworkWallOperator<'_>, count: usize) -> f64 {
    let energy = |contact: &WallContact| match contact.thermal {
        WallThermal::Lumped {
            temperature_k,
            heat_capacity_j_per_k,
            ..
        } => temperature_k * heat_capacity_j_per_k,
        WallThermal::Caloric {
            mass_kg,
            specific_energy_j_kg,
            ..
        } => mass_kg * specific_energy_j_kg,
        _ => 0.0,
    };
    walls.contacts().iter().map(energy).sum::<f64>()
        + (0..count)
            .flat_map(|b| walls.branch_contacts(b).unwrap_or(&[]))
            .map(energy)
            .sum::<f64>()
}

fn gas_and_potential(network: &GasNetwork1d) -> f64 {
    let frame = network.gravity().unwrap();
    let datum = frame.branch_offsets_j_kg[0];
    let mut value = 0.0;
    for (b, branch) in network.branches().iter().enumerate() {
        let (g, origin) = match branch.hydrostatic_reference().unwrap() {
            HydrostaticReference::Isothermal(r) => {
                (r.acceleration.value(), r.reference_position.value())
            }
            HydrostaticReference::Isentropic(r) => {
                (r.acceleration.value(), r.reference_position.value())
            }
        };
        for (i, q) in branch.states().iter().enumerate() {
            let x = 0.5 * (branch.mesh().positions()[i] + branch.mesh().positions()[i + 1]);
            let phi = (frame.branch_offsets_j_kg[b] - datum) - g * (x - origin);
            value += branch.mesh().volumes()[i] * (q.energy + q.density * phi);
        }
    }
    for (i, plenum) in network.plenums().iter().enumerate() {
        let q = plenum.state();
        value +=
            plenum.volume_m3() * (q.energy + q.density * (frame.plenum_potentials_j_kg[i] - datum));
    }
    value
}

// Independent inversion of the fixture's Cv=2.5+0.4*T integral. Its zero at
// Tmin cancels before this arithmetic; no production gas inversion is called.
fn temperature_from_caloric(value: f64) -> f64 {
    2.0 * value / (2.5 + (6.25 + 0.8 * value).sqrt())
}

fn thermal_root(caloric: f64, mass: f64, wall_t: f64, effective: f64) -> f64 {
    let a = 0.2 * mass;
    let b = 2.5 * mass + effective;
    let c = -mass * caloric - effective * wall_t;
    -2.0 * c / (b + (b * b - 4.0 * a * c).sqrt())
}

#[test]
fn distributed_wall_contacts_preserve_exact_thermal_and_hydrostatic_equilibrium() {
    for isentropic in [false, true] {
        for reconstruction in [Reconstruction::FirstOrder, Reconstruction::Muscl] {
            let gas = gas(1e18);
            let lower = atmosphere(&gas, isentropic, -0.5, 1.0, 3.0);
            let reversed = atmosphere(&gas, isentropic, 0.5, 1.0, 3.0);
            let node = lower.hydrostatic_reference_face(4).unwrap().0;
            let mut network = GasNetwork1d::new_with_gravity(
                vec![lower, reversed],
                vec![Plenum::new(0.3, node).unwrap()],
                vec![
                    PlenumConnection {
                        branch: 0,
                        end: BranchEnd::Right,
                        plenum: 0,
                    },
                    PlenumConnection {
                        branch: 1,
                        end: BranchEnd::Right,
                        plenum: 0,
                    },
                ],
                NetworkGravity {
                    branch_offsets_j_kg: vec![2.0_f64.powi(40); 2],
                    plenum_potentials_j_kg: vec![2.0_f64.powi(40)],
                },
                || false,
            )
            .unwrap();
            let thermal = if isentropic {
                WallThermal::Adiabatic
            } else {
                WallThermal::Lumped {
                    temperature_k: 2.0,
                    heat_capacity_j_per_k: 0.7,
                    conductance_w_per_k: 1e9,
                }
            };
            let branches = vec![
                Some(vec![
                    contact(WallDrag::Linear { rate_per_s: 1e4 }, thermal);
                    4
                ]);
                2
            ];
            let mut walls = NetworkWallOperator::new(
                &network,
                vec![plenum_contact(WallThermal::Lumped {
                    temperature_k: 2.0,
                    heat_capacity_j_per_k: 1.0,
                    conductance_w_per_k: 1e9,
                })],
                || false,
            )
            .unwrap()
            .with_branch_contacts(&network, branches, || false)
            .unwrap();
            let fields: Vec<_> = network
                .branches()
                .iter()
                .map(|b| b.states().to_vec())
                .collect();
            let contacts: Vec<_> = (0..2)
                .map(|b| walls.branch_contacts(b).unwrap().to_vec())
                .collect();
            for _ in 0..20 {
                let report = network
                    .step_with_walls(
                        Time::new(0.003),
                        StepControl {
                            reconstruction,
                            ..StepControl::default()
                        },
                        &mut walls,
                        controls(),
                        || false,
                    )
                    .unwrap();
                for (b, branch) in network.branches().iter().enumerate() {
                    assert_eq!(branch.states(), fields[b]);
                    assert_eq!(walls.branch_contacts(b).unwrap(), contacts[b]);
                    let balance = network.accepted_branch_reports().unwrap()[b];
                    assert_eq!(balance.wall_heat, 0.0);
                    assert_eq!(balance.wall_momentum_impulse, 0.0);
                    close(balance.conservation_residual.momentum, 0.0, 2e-13);
                }
                assert_eq!(network.plenums()[0].state(), node);
                let heat = report.thermal.unwrap();
                assert_eq!(heat.source.heat_to_gas_j, 0.0);
                assert_eq!(heat.source.drag_heating_j, 0.0);
                assert_eq!(heat.gas_wall_potential_conservation_residual_j, Some(0.0));
            }
        }
    }
}

#[test]
fn linear_and_darcy_loss_thermalize_kinetic_energy_without_extra_total_energy() {
    for drag in [
        WallDrag::Linear { rate_per_s: 2.0 },
        WallDrag::Darcy {
            friction_factor: 0.6,
            hydraulic_diameter_m: 0.3,
        },
    ] {
        let velocities = [0.7, -0.7];
        let mut network = periodic(&velocities);
        let initial: Vec<_> = network.branches().iter().map(|b| b.states()[0]).collect();
        let mut walls = NetworkWallOperator::new(&network, vec![], || false)
            .unwrap()
            .with_branch_contacts(
                &network,
                vec![Some(vec![contact(drag, WallThermal::Adiabatic); 4]); 2],
                || false,
            )
            .unwrap();
        let report = network
            .step_with_walls(
                Time::new(0.08),
                StepControl::default(),
                &mut walls,
                controls(),
                || false,
            )
            .unwrap();
        let dt = report.elapsed.value();
        let mut converted = 0.0;
        for (b, &velocity) in velocities.iter().enumerate() {
            let final_velocity = match drag {
                WallDrag::Linear { rate_per_s } => velocity * (-rate_per_s * dt).exp(),
                WallDrag::Darcy {
                    friction_factor,
                    hydraulic_diameter_m,
                } => {
                    velocity
                        / (1.0 + 0.5 * friction_factor * dt * velocity.abs() / hydraulic_diameter_m)
                }
                _ => unreachable!(),
            };
            let caloric = 2.5 * 2.0
                + 0.2 * 4.0
                + 0.5 * (velocity * velocity - final_velocity * final_velocity);
            let temperature = temperature_from_caloric(caloric);
            for &q in network.branches()[b].states() {
                assert_eq!(q.energy, initial[b].energy);
                assert_eq!(q.density, initial[b].density);
                let actual = network.gas().primitive(q, || false).unwrap();
                close(actual.velocity, final_velocity, 2e-14);
                close(actual.temperature, temperature, 2e-12);
                assert!(actual.temperature > 2.0);
            }
            let mass = 0.7;
            converted += 0.5 * mass * (velocity * velocity - final_velocity * final_velocity);
            let branch = network.accepted_branch_reports().unwrap()[b];
            close(
                branch.wall_momentum_impulse,
                mass * (final_velocity - velocity),
                2e-14,
            );
            assert!(branch.wall_momentum_impulse * velocity < 0.0);
            assert_eq!(branch.wall_heat, 0.0);
            close(branch.conservation_residual.momentum, 0.0, 2e-14);
            close(branch.conservation_residual.energy, 0.0, 2e-14);
        }
        let heat = report.thermal.unwrap();
        close(heat.source.drag_heating_j, converted, 3e-14);
        assert_eq!(heat.source.heat_to_gas_j, 0.0);
        assert_eq!(report.after.energy, report.before.energy);
        close(heat.gas_wall_conservation_residual_j, 0.0, 2e-14);
    }
}

#[test]
fn finite_wall_exchange_uses_drag_heated_original_gas_calorics_once() {
    let mut network = periodic(&[1.0]);
    let mut walls = NetworkWallOperator::new(&network, vec![], || false)
        .unwrap()
        .with_branch_contacts(
            &network,
            vec![Some(vec![
                contact(
                    WallDrag::Linear { rate_per_s: 3.0 },
                    WallThermal::Lumped {
                        temperature_k: 1.0,
                        heat_capacity_j_per_k: 0.5,
                        conductance_w_per_k: 2.0,
                    },
                );
                4
            ])],
            || false,
        )
        .unwrap();
    let initial_energy = network.totals().energy + wall_energy(&walls, 1);
    let report = network
        .step_with_walls(
            Time::new(0.04),
            StepControl::default(),
            &mut walls,
            controls(),
            || false,
        )
        .unwrap();
    let dt = report.elapsed.value();
    let velocity = (-3.0 * dt).exp();
    let caloric = 2.5 * 2.0 + 0.2 * 4.0 + 0.5 * (1.0 - velocity * velocity);
    let mass = 0.7 * 0.25;
    let integrated = 2.0 * dt;
    let effective = 0.5 * integrated / (0.5 + integrated);
    let temperature = thermal_root(caloric, mass, 1.0, effective);
    let heat = mass * (2.5 * temperature + 0.2 * temperature * temperature - caloric);
    assert!(heat < 0.0);
    for (q, wall) in network.branches()[0]
        .states()
        .iter()
        .zip(walls.branch_contacts(0).unwrap())
    {
        let state = network.gas().primitive(*q, || false).unwrap();
        close(state.velocity, velocity, 2e-14);
        close(state.temperature, temperature, 2e-10);
        let WallThermal::Lumped { temperature_k, .. } = wall.thermal else {
            panic!("finite wall");
        };
        close(temperature_k, 1.0 - heat / 0.5, 2e-10);
    }
    let branch = network.accepted_branch_reports().unwrap()[0];
    close(branch.wall_heat, 4.0 * heat, 2e-10);
    close(
        branch.after.energy - branch.before.energy,
        4.0 * heat,
        2e-10,
    );
    close(branch.conservation_residual.energy, 0.0, 2e-13);
    close(branch.conservation_residual.momentum, 0.0, 2e-13);
    close(
        network.totals().energy + wall_energy(&walls, 1),
        initial_energy,
        2e-10,
    );
    close(
        report.thermal.unwrap().source.drag_heating_j,
        0.5 * 0.7 * (1.0 - velocity * velocity),
        3e-14,
    );
}

#[test]
fn distributed_heat_drives_gravity_flow_with_closed_wall_energy_and_gauge_invariance() {
    let mut base = chain(0.0, 0.0);
    let shift = 2.0_f64.powi(40);
    let mut shifted = chain(shift, 1e18);
    let mut walls = driven(&base);
    let mut shifted_walls = driven(&shifted);
    let initial = gas_and_potential(&base) + wall_energy(&walls, 2);
    let mut transported = 0.0;
    let mut converted = 0.0;
    for _ in 0..50 {
        let a = base
            .step_with_walls(
                Time::new(0.002),
                StepControl::default(),
                &mut walls,
                controls(),
                || false,
            )
            .unwrap();
        let b = shifted
            .step_with_walls(
                Time::new(0.002),
                StepControl::default(),
                &mut shifted_walls,
                controls(),
                || false,
            )
            .unwrap();
        assert_eq!(a.after, b.after);
        assert_eq!(a.conservation_residual, b.conservation_residual);
        for (index, (a, b)) in base.branches().iter().zip(shifted.branches()).enumerate() {
            assert_eq!(a.states(), b.states());
            assert_eq!(
                walls.branch_contacts(index),
                shifted_walls.branch_contacts(index)
            );
            let report = base.accepted_branch_reports().unwrap()[index];
            close(report.conservation_residual.energy, 0.0, 3e-12);
            close(report.conservation_residual.momentum, 0.0, 3e-13);
        }
        assert_eq!(base.plenums(), shifted.plenums());
        let (a_heat, b_heat) = (a.thermal.unwrap(), b.thermal.unwrap());
        assert_eq!(a_heat.source, b_heat.source);
        assert_eq!(
            a_heat.gas_wall_potential_conservation_residual_j,
            b_heat.gas_wall_potential_conservation_residual_j
        );
        assert_eq!(
            b_heat.absolute_gas_wall_potential_conservation_residual_j,
            Some(
                a_heat.gas_wall_potential_conservation_residual_j.unwrap()
                    + shift * a.conservation_residual.mass
            )
        );
        assert_eq!(a_heat.source.external_heat_to_gas_j, 0.0);
        close(
            a_heat.gas_wall_potential_conservation_residual_j.unwrap(),
            0.0,
            3e-11,
        );
        converted += a_heat.source.drag_heating_j;
        transported += base
            .accepted_port_exchanges()
            .unwrap()
            .iter()
            .map(|p| p.mass_to_plenum.abs())
            .sum::<f64>();
    }
    assert!(
        transported > 1e-6,
        "heat must generate actual port exchange"
    );
    assert!(
        converted > 0.0,
        "the generated flow must experience real wall drag"
    );
    assert!(
        base.branches()
            .iter()
            .flat_map(|b| b.states())
            .any(|q| q.momentum.abs() > 1e-6)
    );
    close(
        gas_and_potential(&base) + wall_energy(&walls, 2),
        initial,
        3e-10,
    );
}

// Exact constitutive specimen retained by the original thermal-law owner:
// T(h)=3 on h in[2,5], with sensible branches on either side. No gas EOS stub.
struct LatentWall;
impl ThermalLaw for LatentWall {
    fn temperature_k(&self, h: f64) -> Result<f64, DepositionError> {
        if !(0.0..=10.0).contains(&h) {
            return Err(DepositionError::MaterialDomain);
        }
        Ok(if h < 2.0 {
            1.0 + h
        } else if h <= 5.0 {
            3.0
        } else {
            3.0 + 0.5 * (h - 5.0)
        })
    }
    fn maximum_temperature_slope(&self) -> f64 {
        1.0
    }
}

fn latent(id: u64, range: [f64; 2], conductance: f64) -> BranchWallContact<'static> {
    let source = MaterialWallCaloric::from_thermal_law(
        LatentWall,
        DepositedCell {
            id: CellId(id),
            segment_id: 1,
            body_id: 4,
            material_id: MaterialId(29),
            mass_kg: 0.1,
            specific_h_j_kg: 4.0,
            geometry: CellGeometry {
                start_m: [0.0; 3],
                end_m: [1.0, 0.0, 0.0],
                volume_m3: 0.001,
                surface_area_m2: 1.0,
            },
            boundary: ThermalBoundary {
                ambient_k: 3.0,
                htc_w_m2_k: 0.0,
                reservoir: None,
                source_w: 0.0,
            },
        },
        range,
    )
    .unwrap();
    BranchWallContact {
        contact: source
            .initial_contact(conductance, WallDrag::Linear { rate_per_s: 0.7 })
            .unwrap(),
        heat_storage: Some(source),
    }
}

#[test]
fn moving_branch_retains_original_latent_inventory_and_energy_convention() {
    let mut network = periodic(&[0.4]);
    let mut walls = NetworkWallOperator::new(&network, vec![], || false)
        .unwrap()
        .with_branch_contacts(
            &network,
            vec![Some((0..4).map(|i| latent(i, [0.0, 10.0], 0.2)).collect())],
            || false,
        )
        .unwrap();
    let initial = network.totals().energy + wall_energy(&walls, 1);
    let report = network
        .step_with_walls(
            Time::new(0.02),
            StepControl::default(),
            &mut walls,
            controls(),
            || false,
        )
        .unwrap();
    let dt = report.elapsed.value();
    let velocity = 0.4 * (-0.7 * dt).exp();
    let caloric = 2.5 * 2.0 + 0.2 * 4.0 + 0.5 * (0.16 - velocity * velocity);
    let mass = 0.7 * 0.25;
    let temperature = thermal_root(caloric, mass, 3.0, 0.2 * dt);
    let heat = mass * (2.5 * temperature + 0.2 * temperature * temperature - caloric);
    for (i, q) in network.branches()[0].states().iter().enumerate() {
        close(
            network.gas().primitive(*q, || false).unwrap().temperature,
            temperature,
            3e-9,
        );
        let original = walls.branch_heat_storage(0).unwrap()[i].as_ref().unwrap();
        let state = walls.branch_contacts(0).unwrap()[i].thermal;
        assert_eq!(original.temperature_k(state).unwrap(), 3.0);
        assert_eq!(original.storage_basis(), ThermalStorageBasis::Enthalpy);
        assert_eq!(original.initial_cell().unwrap().id, CellId(i as u64));
        let WallThermal::Caloric {
            specific_energy_j_kg,
            ..
        } = state
        else {
            panic!("latent inventory");
        };
        close(specific_energy_j_kg, 4.0 - heat / 0.1, 3e-8);
    }
    close(
        network.totals().energy + wall_energy(&walls, 1),
        initial,
        2e-10,
    );
    close(
        network.accepted_branch_reports().unwrap()[0].wall_heat,
        4.0 * heat,
        3e-9,
    );
}

#[test]
fn final_branch_source_poll_cancellation_preserves_network_walls_and_reports() {
    let mut probe = chain(0.0, 0.0);
    let mut network = chain(0.0, 0.0);
    let mut probe_walls = driven(&probe);
    let mut walls = driven(&network);
    for (network, walls) in [(&mut probe, &mut probe_walls), (&mut network, &mut walls)] {
        network
            .step_with_walls(
                Time::new(0.002),
                StepControl::default(),
                walls,
                controls(),
                || false,
            )
            .unwrap();
    }
    let fields: Vec<_> = network
        .branches()
        .iter()
        .map(|b| b.states().to_vec())
        .collect();
    let contacts: Vec<_> = (0..2)
        .map(|b| walls.branch_contacts(b).unwrap().to_vec())
        .collect();
    let plenums = network.plenums().to_vec();
    let time = network.time();
    let report = network.last_report();
    let branch_reports = network.accepted_branch_reports().unwrap().to_vec();
    let ports = network.accepted_port_exchanges().unwrap().to_vec();
    let source = walls.last_report();
    let mut polls = 0;
    probe
        .step_with_walls(
            Time::new(0.002),
            StepControl::default(),
            &mut probe_walls,
            controls(),
            || {
                polls += 1;
                false
            },
        )
        .unwrap();
    let mut observed = 0;
    assert!(matches!(
        network.step_with_walls(
            Time::new(0.002),
            StepControl::default(),
            &mut walls,
            controls(),
            || {
                observed += 1;
                observed == polls
            }
        ),
        Err(Error::Cancelled)
    ));
    assert_eq!(observed, polls);
    for (b, branch) in network.branches().iter().enumerate() {
        assert_eq!(branch.states(), fields[b]);
        assert_eq!(walls.branch_contacts(b).unwrap(), contacts[b]);
    }
    assert_eq!(network.plenums(), plenums);
    assert_eq!(network.time(), time);
    assert_eq!(network.last_report(), report);
    assert_eq!(network.accepted_branch_reports().unwrap(), branch_reports);
    assert_eq!(network.accepted_port_exchanges().unwrap(), ports);
    assert_eq!(walls.last_report(), source);
    network
        .step_with_walls(
            Time::new(0.002),
            StepControl::default(),
            &mut walls,
            controls(),
            || false,
        )
        .unwrap();
    close(network.totals().energy, probe.totals().energy, 2e-11);
    close(wall_energy(&walls, 2), wall_energy(&probe_walls, 2), 2e-11);
}

#[test]
fn later_branch_source_refusal_and_wrong_geometry_publish_nothing() {
    let mut network = periodic(&[0.4, -0.4]);
    let mut inputs = vec![
        Some(vec![
            contact(
                WallDrag::Linear { rate_per_s: 1.0 },
                WallThermal::Reservoir {
                    temperature_k: 4.0,
                    conductance_w_per_k: 1.0
                }
            );
            4
        ]);
        2
    ];
    inputs[1].as_mut().unwrap()[3] = latent(9, [3.999, 4.0], 1.0);
    let mut walls = NetworkWallOperator::new(&network, vec![], || false)
        .unwrap()
        .with_branch_contacts(&network, inputs, || false)
        .unwrap();
    let fields: Vec<_> = network
        .branches()
        .iter()
        .map(|b| b.states().to_vec())
        .collect();
    let contacts: Vec<_> = (0..2)
        .map(|b| walls.branch_contacts(b).unwrap().to_vec())
        .collect();
    assert!(
        network
            .step_with_walls(
                Time::new(0.04),
                StepControl {
                    maximum_halvings: 0,
                    ..StepControl::default()
                },
                &mut walls,
                controls(),
                || false
            )
            .is_err()
    );
    for (b, branch) in network.branches().iter().enumerate() {
        assert_eq!(branch.states(), fields[b]);
        assert_eq!(walls.branch_contacts(b).unwrap(), contacts[b]);
    }
    assert_eq!(network.time().value(), 0.0);
    assert_eq!(network.last_report(), None);
    assert_eq!(walls.last_report(), None);
    let gas = network.gas().clone();
    let mut different = GasNetwork1d::new(
        (0..2)
            .map(|_| {
                Euler1d::new(
                    gas.clone(),
                    Length::new(0.3),
                    0.7,
                    vec![gas.conserved(1.0, 0.4, 2.0).unwrap(); 4],
                    Boundaries::Periodic,
                    || false,
                )
                .unwrap()
            })
            .collect(),
        vec![],
        vec![],
        || false,
    )
    .unwrap();
    assert!(
        different
            .step_with_walls(
                Time::new(0.001),
                StepControl::default(),
                &mut walls,
                controls(),
                || false
            )
            .is_err()
    );
    assert_eq!(different.time().value(), 0.0);
    assert!(
        NetworkWallOperator::new(&network, vec![], || false)
            .unwrap()
            .with_branch_contacts(&network, vec![None], || false)
            .is_err()
    );
    assert!(
        NetworkWallOperator::new(&network, vec![], || false)
            .unwrap()
            .with_branch_contacts(&network, vec![Some(vec![]), None], || false)
            .is_err()
    );
}
