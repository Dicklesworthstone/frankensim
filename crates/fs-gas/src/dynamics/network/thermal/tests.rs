//! Independent thermal storage and equilibrium checks on original source owners.

use super::*;
use crate::dynamics::hydrostatic::IsothermalAtmosphere;
use crate::dynamics::network::{BranchEnd, NetworkGravity, Plenum, PlenumConnection};
use crate::dynamics::{Boundaries, Boundary, Euler1d, Mesh1d, StepControl};
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
        absolute_energy_tolerance_j: 1e-9,
        ..WallControl::default()
    }
}

fn isolated(temperature: f64, plenums: usize) -> GasNetwork1d {
    let gas = gas(0.0);
    let state = gas.conserved(1.0, 0.0, temperature).unwrap();
    let branch = Euler1d::new(gas, Length::new(1.0), 0.1, vec![state; 2], WALLS, || false).unwrap();
    GasNetwork1d::new(
        vec![branch],
        (0..plenums)
            .map(|_| Plenum::new(0.4, state).unwrap())
            .collect(),
        vec![],
        || false,
    )
    .unwrap()
}

fn contact(thermal: WallThermal) -> PlenumThermalContact<'static> {
    PlenumThermalContact {
        thermal,
        heat_storage: None,
    }
}

fn lumped(temperature: f64, capacity: f64, conductance: f64) -> PlenumThermalContact<'static> {
    contact(WallThermal::Lumped {
        temperature_k: temperature,
        heat_capacity_j_per_k: capacity,
        conductance_w_per_k: conductance,
    })
}

fn gas_temperature(network: &GasNetwork1d, plenum: usize) -> f64 {
    network
        .gas()
        .primitive(network.plenums()[plenum].state(), || false)
        .unwrap()
        .temperature
}

fn wall_temperature(operator: &PlenumThermalOperator<'_>, index: usize) -> f64 {
    match operator.contacts()[index].thermal {
        WallThermal::Lumped { temperature_k, .. }
        | WallThermal::Reservoir { temperature_k, .. } => temperature_k,
        state => operator.heat_storage()[index]
            .as_ref()
            .unwrap()
            .temperature_k(state)
            .unwrap(),
    }
}

// Independent original fixture integrals: R=1, Cv(T)=2.5+0.4*T.
fn sensible(temperature: f64) -> f64 {
    2.5 * (temperature - 0.001) + 0.2 * (temperature * temperature - 0.000_001)
}

fn thermal_root(initial_t: f64, mass: f64, wall_t: f64, effective: f64) -> f64 {
    let a = 0.2 * mass;
    let b = 2.5 * mass + effective;
    let c = -mass * (2.5 * initial_t + 0.2 * initial_t * initial_t) - effective * wall_t;
    -2.0 * c / (b + (b * b - 4.0 * a * c).sqrt())
}

fn column(gas: &FrozenGas, pressure: f64) -> Euler1d {
    Euler1d::new_isothermal_atmosphere(
        gas.clone(),
        Mesh1d::new(
            vec![0.0, 0.1, 0.3, 0.6, 1.0],
            vec![0.7; 5],
            vec![0.07, 0.14, 0.21, 0.28],
        )
        .unwrap(),
        IsothermalAtmosphere {
            acceleration: Acceleration::new(-0.5),
            temperature: Temperature::new(2.0),
            reference_pressure: Pressure::new(pressure),
            reference_position: Length::new(0.0),
        },
        WALLS,
        || false,
    )
    .unwrap()
}

fn chain(datum: f64, formation: f64) -> GasNetwork1d {
    let gas = gas(formation);
    let lower = column(&gas, 3.0);
    let bottom = lower.hydrostatic_reference_face(0).unwrap();
    let middle = lower.hydrostatic_reference_face(4).unwrap();
    let upper = column(&gas, middle.1.pressure);
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

fn driven(network: &GasNetwork1d) -> PlenumThermalOperator<'static> {
    PlenumThermalOperator::new(
        network,
        vec![
            lumped(4.0, 2.0, 0.8),
            contact(WallThermal::Adiabatic),
            contact(WallThermal::Reservoir {
                temperature_k: 1.0,
                conductance_w_per_k: 0.4,
            }),
        ],
        || false,
    )
    .unwrap()
}

#[test]
fn exact_hydrostatic_and_thermal_equilibrium_survives_stiff_contacts() {
    let mut network = chain(2.0_f64.powi(40), 1e18);
    let original: Vec<_> = network
        .branches()
        .iter()
        .map(|b| b.states().to_vec())
        .collect();
    let plenums = network.plenums().to_vec();
    let mut thermal = PlenumThermalOperator::new(
        &network,
        vec![
            lumped(2.0, 3.0, 1e12),
            contact(WallThermal::Reservoir {
                temperature_k: 2.0,
                conductance_w_per_k: 1e12,
            }),
            lumped(2.0, 4.0, 1e12),
        ],
        || false,
    )
    .unwrap();
    let contacts = thermal.contacts().to_vec();
    for _ in 0..25 {
        let report = network
            .step_with_thermal(
                Time::new(0.003),
                StepControl::default(),
                &mut thermal,
                controls(),
                || false,
            )
            .unwrap();
        for (branch, initial) in network.branches().iter().zip(&original) {
            assert_eq!(branch.states(), initial);
        }
        assert_eq!(network.plenums(), plenums);
        assert_eq!(thermal.contacts(), contacts);
        let heat = report.thermal.unwrap();
        assert_eq!(heat.source.heat_to_gas_j, 0.0);
        assert_eq!(heat.source.wall_energy_change_j, 0.0);
        assert_eq!(heat.source.maximum_iterations, 0);
        assert_eq!(heat.gas_wall_potential_conservation_residual_j, Some(0.0));
        assert_eq!(thermal.last_report(), Some(heat.source));
    }
}

#[test]
fn variable_cv_finite_walls_use_actual_accepted_duration_without_crossing() {
    for (initial_t, wall_t) in [(1.5, 5.0), (5.0, 1.5)] {
        for conductance in [0.3, 1e4] {
            let mut network = isolated(initial_t, 1);
            let mut thermal = PlenumThermalOperator::new(
                &network,
                vec![lumped(wall_t, 3.0, conductance)],
                || false,
            )
            .unwrap();
            let report = network
                .step_with_thermal(
                    Time::new(1000.0),
                    StepControl::default(),
                    &mut thermal,
                    controls(),
                    || false,
                )
                .unwrap();
            assert!(report.elapsed.value() < 1000.0);
            let heat = report.thermal.unwrap();
            assert_eq!(heat.source.elapsed, report.elapsed);
            let integrated = conductance * report.elapsed.value();
            let effective = 3.0 * integrated / (3.0 + integrated);
            let expected = thermal_root(initial_t, 0.4, wall_t, effective);
            let q = 0.4 * (sensible(expected) - sensible(initial_t));
            close(gas_temperature(&network, 0), expected, 3e-11);
            close(wall_temperature(&thermal, 0), wall_t - q / 3.0, 3e-11);
            close(heat.source.heat_to_gas_j, q, 1e-10);
            close(heat.source.wall_energy_change_j, -q, 1e-10);
            close(heat.gas_wall_conservation_residual_j, 0.0, 1e-10);
            let difference = wall_temperature(&thermal, 0) - gas_temperature(&network, 0);
            assert!(difference * (wall_t - initial_t) >= 0.0);
            assert_eq!(heat.source.external_heat_to_gas_j, 0.0);
        }
    }
    let mut network = isolated(1.5, 1);
    let mut thermal =
        PlenumThermalOperator::new(&network, vec![lumped(5.0, 3.0, 20.0)], || false).unwrap();
    for _ in 0..80 {
        network
            .step_with_thermal(
                Time::new(0.1),
                StepControl::default(),
                &mut thermal,
                controls(),
                || false,
            )
            .unwrap();
    }
    let equilibrium = thermal_root(1.5, 0.4, 5.0, 3.0);
    close(gas_temperature(&network, 0), equilibrium, 3e-9);
    close(wall_temperature(&thermal, 0), equilibrium, 3e-9);
}

// Constitutive test specimens passed through the ORIGINAL immutable thermal
// law/storage owner. Neither specimen replaces the gas NASA caloric source.
struct QuadraticWall;
impl ThermalLaw for QuadraticWall {
    fn storage_basis(&self) -> ThermalStorageBasis {
        ThermalStorageBasis::InternalEnergy
    }
    fn temperature_k(&self, energy: f64) -> Result<f64, DepositionError> {
        if !energy.is_finite() || energy <= 0.0 {
            return Err(DepositionError::MaterialDomain);
        }
        Ok((2.0 * energy).sqrt())
    }
    fn maximum_temperature_slope(&self) -> f64 {
        100.0
    }
}

struct LatentWall;
impl ThermalLaw for LatentWall {
    fn temperature_k(&self, energy: f64) -> Result<f64, DepositionError> {
        if !(0.0..=10.0).contains(&energy) {
            return Err(DepositionError::MaterialDomain);
        }
        Ok(if energy < 2.0 {
            1.0 + energy
        } else if energy <= 5.0 {
            3.0
        } else {
            3.0 + 0.5 * (energy - 5.0)
        })
    }
    fn maximum_temperature_slope(&self) -> f64 {
        1.0
    }
}

fn wall_cell(energy: f64) -> DepositedCell {
    DepositedCell {
        id: CellId(17),
        segment_id: 3,
        body_id: 8,
        material_id: MaterialId(29),
        mass_kg: 1.0,
        specific_h_j_kg: energy,
        geometry: CellGeometry {
            start_m: [0.0; 3],
            end_m: [1.0, 0.0, 0.0],
            volume_m3: 0.1,
            surface_area_m2: 1.0,
        },
        boundary: ThermalBoundary {
            ambient_k: 2.0,
            htc_w_m2_k: 0.0,
            reservoir: None,
            source_w: 0.0,
        },
    }
}

fn material_contact(
    source: MaterialWallCaloric<'static>,
    conductance: f64,
) -> PlenumThermalContact<'static> {
    PlenumThermalContact {
        thermal: source
            .initial_contact(conductance, WallDrag::None)
            .unwrap()
            .thermal,
        heat_storage: Some(source),
    }
}

fn wall_energy(operator: &PlenumThermalOperator<'_>, index: usize) -> f64 {
    match operator.contacts()[index].thermal {
        WallThermal::Caloric {
            mass_kg,
            specific_energy_j_kg,
            ..
        } => mass_kg * specific_energy_j_kg,
        WallThermal::Lumped {
            temperature_k,
            heat_capacity_j_per_k,
            ..
        } => temperature_k * heat_capacity_j_per_k,
        _ => 0.0,
    }
}

#[test]
fn original_nonlinear_wall_integral_controls_equilibration() {
    let mut network = isolated(2.0, 1);
    let source =
        MaterialWallCaloric::from_thermal_law(QuadraticWall, wall_cell(12.5), [0.01, 50.0])
            .unwrap();
    let mut thermal =
        PlenumThermalOperator::new(&network, vec![material_contact(source, 2.0)], || false)
            .unwrap();
    let report = network
        .step_with_thermal(
            Time::new(0.04),
            StepControl::default(),
            &mut thermal,
            controls(),
            || false,
        )
        .unwrap();
    let integrated = 2.0 * report.elapsed.value();
    // Independent elimination: gas's exact quadratic energy and the original
    // quadratic wall energy express BOTH temperatures as functions of Q.
    let mut lo = 0.0;
    let mut hi = integrated * 3.0;
    for _ in 0..160 {
        let q = 0.5 * (lo + hi);
        let d = 2.5 * 2.0 + 0.2 * 4.0 + q / 0.4;
        let tg = 2.0 * d / (2.5 + (6.25 + 0.8 * d).sqrt());
        let tw = (25.0 - 2.0 * q).sqrt();
        if q - integrated * (tw - tg) < 0.0 {
            lo = q;
        } else {
            hi = q;
        }
    }
    let q = 0.5 * (lo + hi);
    close(report.thermal.unwrap().source.heat_to_gas_j, q, 1e-9);
    close(wall_energy(&thermal, 0), 12.5 - q, 1e-9);
    close(wall_temperature(&thermal, 0), (25.0 - 2.0 * q).sqrt(), 1e-9);
    let source = thermal.heat_storage()[0].as_ref().unwrap();
    assert_eq!(source.storage_basis(), ThermalStorageBasis::InternalEnergy);
    assert_eq!(source.initial_cell().unwrap().material_id, MaterialId(29));
}

#[test]
fn original_latent_inventory_changes_without_an_apparent_heat_capacity() {
    let mut network = isolated(1.0, 1);
    let source =
        MaterialWallCaloric::from_thermal_law(LatentWall, wall_cell(4.0), [0.0, 10.0]).unwrap();
    let mut thermal =
        PlenumThermalOperator::new(&network, vec![material_contact(source, 3.0)], || false)
            .unwrap();
    let report = network
        .step_with_thermal(
            Time::new(0.04),
            StepControl::default(),
            &mut thermal,
            controls(),
            || false,
        )
        .unwrap();
    let expected_t = thermal_root(1.0, 0.4, 3.0, 3.0 * report.elapsed.value());
    let q = 0.4 * (sensible(expected_t) - sensible(1.0));
    close(gas_temperature(&network, 0), expected_t, 2e-9);
    close(wall_energy(&thermal, 0), 4.0 - q, 2e-9);
    assert_eq!(wall_temperature(&thermal, 0), 3.0);
    assert!(wall_energy(&thermal, 0) > 2.0 && wall_energy(&thermal, 0) < 4.0);
    assert_eq!(
        thermal.heat_storage()[0].as_ref().unwrap().storage_basis(),
        ThermalStorageBasis::Enthalpy
    );
    close(
        report.thermal.unwrap().gas_wall_conservation_residual_j,
        0.0,
        1e-10,
    );
}

// Independent centered gas+potential inventory from actual retained fields.
fn physical_storage(network: &GasNetwork1d) -> f64 {
    let frame = network.gravity().unwrap();
    let datum = frame.branch_offsets_j_kg[0];
    let mut energy = 0.0;
    for (b, branch) in network.branches().iter().enumerate() {
        for (i, state) in branch.states().iter().enumerate() {
            let x = 0.5 * (branch.mesh().positions()[i] + branch.mesh().positions()[i + 1]);
            let phi = (frame.branch_offsets_j_kg[b] - datum) + 0.5 * x;
            energy += branch.mesh().volumes()[i] * (state.energy + state.density * phi);
        }
    }
    for (i, plenum) in network.plenums().iter().enumerate() {
        let phi = frame.plenum_potentials_j_kg[i] - datum;
        energy += plenum.volume_m3() * (plenum.state().energy + plenum.state().density * phi);
    }
    energy
}

#[test]
fn heating_drives_gravitating_flow_with_conservative_wall_storage_and_datum_invariance() {
    let mut base = chain(0.0, 0.0);
    let shift = 2.0_f64.powi(40);
    let mut shifted = chain(shift, 1e18);
    let mut wall = driven(&base);
    let mut shifted_wall = driven(&shifted);
    let initial = physical_storage(&base) + wall_energy(&wall, 0);
    let mut reservoir_heat = 0.0;
    let mut transported = 0.0;
    for _ in 0..50 {
        let a = base
            .step_with_thermal(
                Time::new(0.002),
                StepControl::default(),
                &mut wall,
                controls(),
                || false,
            )
            .unwrap();
        let b = shifted
            .step_with_thermal(
                Time::new(0.002),
                StepControl::default(),
                &mut shifted_wall,
                controls(),
                || false,
            )
            .unwrap();
        assert_eq!(a.elapsed, b.elapsed);
        assert_eq!(a.after, b.after);
        assert_eq!(a.conservation_residual, b.conservation_residual);
        assert_eq!(
            a.gravity.unwrap().gravity_work,
            b.gravity.unwrap().gravity_work
        );
        for (a, b) in base.branches().iter().zip(shifted.branches()) {
            assert_eq!(a.states(), b.states());
        }
        assert_eq!(base.plenums(), shifted.plenums());
        assert_eq!(wall.contacts(), shifted_wall.contacts());
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
        reservoir_heat += a_heat.source.external_heat_to_gas_j;
        transported += base
            .accepted_port_exchanges()
            .unwrap()
            .iter()
            .map(|p| p.mass_to_plenum.abs())
            .sum::<f64>();
        close(
            a_heat.gas_wall_potential_conservation_residual_j.unwrap(),
            0.0,
            2e-10,
        );
        close(a.conservation_residual.energy, 0.0, 2e-10);
    }
    assert!(transported > 1e-6);
    assert!(reservoir_heat < -1e-3);
    close(
        physical_storage(&base) + wall_energy(&wall, 0) - initial,
        reservoir_heat,
        2e-9,
    );
}

#[test]
fn final_thermal_publication_cancellation_preserves_both_original_owners() {
    let mut probe = chain(0.0, 0.0);
    let mut network = chain(0.0, 0.0);
    let mut probe_wall = driven(&probe);
    let mut wall = driven(&network);
    for (network, thermal) in [(&mut probe, &mut probe_wall), (&mut network, &mut wall)] {
        network
            .step_with_thermal(
                Time::new(0.002),
                StepControl::default(),
                thermal,
                controls(),
                || false,
            )
            .unwrap();
    }
    let states: Vec<_> = network
        .branches()
        .iter()
        .map(|b| b.states().to_vec())
        .collect();
    let plenums = network.plenums().to_vec();
    let contacts = wall.contacts().to_vec();
    let report = network.last_report();
    let heat = wall.last_report();
    let time = network.time();
    let ports = network.accepted_port_exchanges().unwrap().to_vec();
    let mut polls = 0;
    probe
        .step_with_thermal(
            Time::new(0.002),
            StepControl::default(),
            &mut probe_wall,
            controls(),
            || {
                polls += 1;
                false
            },
        )
        .unwrap();
    let mut observed = 0;
    assert!(matches!(
        network.step_with_thermal(
            Time::new(0.002),
            StepControl::default(),
            &mut wall,
            controls(),
            || {
                observed += 1;
                observed == polls
            }
        ),
        Err(Error::Cancelled)
    ));
    assert_eq!(observed, polls);
    for (branch, before) in network.branches().iter().zip(states) {
        assert_eq!(branch.states(), before);
    }
    assert_eq!(network.plenums(), plenums);
    assert_eq!(wall.contacts(), contacts);
    assert_eq!(network.last_report(), report);
    assert_eq!(wall.last_report(), heat);
    assert_eq!(network.time(), time);
    assert_eq!(network.accepted_port_exchanges().unwrap(), ports);
    network
        .step_with_thermal(
            Time::new(0.002),
            StepControl::default(),
            &mut wall,
            controls(),
            || false,
        )
        .unwrap();
    assert_eq!(network.plenums(), probe.plenums());
    assert_eq!(wall.contacts(), probe_wall.contacts());
}

#[test]
fn later_original_source_refusal_rolls_back_earlier_plenum_heat() {
    let mut network = isolated(2.0, 2);
    let source =
        MaterialWallCaloric::from_thermal_law(LatentWall, wall_cell(4.0), [3.999, 4.0]).unwrap();
    let mut thermal = PlenumThermalOperator::new(
        &network,
        vec![
            contact(WallThermal::Reservoir {
                temperature_k: 4.0,
                conductance_w_per_k: 1.0,
            }),
            material_contact(source, 1.0),
        ],
        || false,
    )
    .unwrap();
    let before = network.plenums().to_vec();
    let contacts = thermal.contacts().to_vec();
    assert!(
        network
            .step_with_thermal(
                Time::new(0.04),
                StepControl {
                    maximum_halvings: 0,
                    ..StepControl::default()
                },
                &mut thermal,
                controls(),
                || false,
            )
            .is_err()
    );
    assert_eq!(network.plenums(), before);
    assert_eq!(thermal.contacts(), contacts);
    assert_eq!(network.time().value(), 0.0);
    assert_eq!(network.last_report(), None);
    assert_eq!(thermal.last_report(), None);
    assert!(PlenumThermalOperator::new(&network, vec![], || false).is_err());
    let mut other = isolated(2.0, 1);
    assert!(
        other
            .step_with_thermal(
                Time::new(0.04),
                StepControl::default(),
                &mut thermal,
                controls(),
                || false,
            )
            .is_err()
    );
}
