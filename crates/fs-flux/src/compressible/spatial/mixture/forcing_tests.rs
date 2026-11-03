use super::*;
use crate::compressible::mixture::tests::{close, gas, source};
use crate::compressible::spatial::AxisBoundaries;
use fs_thermochem::UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R;

const PERIODIC: Boundaries3d = Boundaries3d {
    x: AxisBoundaries::Periodic,
    y: AxisBoundaries::Periodic,
    z: AxisBoundaries::Periodic,
};

fn material(reference: f64) -> MixtureGas {
    gas(vec![
        source("A", R, 3.5, 0.0, reference),
        source("B", 2.0 * R, 3.5, 0.0, -2.0 * reference),
        source("absent", 0.5 * R, 3.5, 0.0, 3.0 * reference),
    ])
}

fn field(reference: f64, variable_density: bool) -> MixtureEuler3d {
    let gas = material(reference);
    let grid = CartesianGrid::new(3, 4, 3, 0.3, 0.4, 0.5).unwrap();
    let mut states = vec![0.0; grid.len() * 7];
    for (i, row) in states.chunks_exact_mut(7).enumerate() {
        let rho = if variable_density {
            1.0 + 0.25 * grid.coordinates(i)[1] as f64
        } else {
            2.0
        };
        state_into(
            &gas,
            &[0.4 * rho, 0.6 * rho, 0.0],
            if variable_density {
                [0.0; 3]
            } else {
                [0.03, -0.02, 0.01]
            },
            1.0,
            row,
        )
        .unwrap();
    }
    let bytes = MixtureEuler3d::required_workspace_bytes(&gas, &grid).unwrap();
    MixtureEuler3d::new(gas, grid, states, PERIODIC, bytes, || false).unwrap()
}

fn control() -> ForceControl3d {
    ForceControl3d { maximum_mach: 0.5 }
}

fn work(flow: &MixtureEuler3d) -> ForcingWorkspace3d {
    let bytes = flow.required_forcing_workspace_bytes().unwrap();
    let result = flow.forcing_workspace(bytes).unwrap();
    assert_eq!(result.workspace_bytes(), bytes);
    result
}

fn balance(report: &StepReport) {
    assert!(report.maximum_cfl <= StepControl::default().cfl * (1.0 + 1e-13));
    for residual in report.conservation_residual.species_masses {
        close(residual, 0.0, 3e-12);
    }
    for residual in report.conservation_residual.momentum {
        close(residual, 0.0, 3e-12);
    }
    close(report.conservation_residual.energy, 0.0, 3e-12);
    close(report.absolute_energy_residual_j, 0.0, 3e-10);
}

#[test]
fn g1_uniform_force_uses_exact_kinetic_work_and_original_caloric_gauge() {
    let acceleration = [0.2, -0.1, 0.05];
    let mut reference = None;
    for (formation, mode) in [(0.0, 0), (0.0, 1), (0.0, 2), (0.0, 3), (1e12, 0)] {
        let mut flow = field(formation, false);
        let mut workspace = work(&flow);
        let initial = flow
            .cell(0)
            .unwrap()
            .primitive(flow.gas(), || false)
            .unwrap();
        let initial_species = flow.cell(0).unwrap().partial_densities.to_vec();
        let acceleration_field = vec![acceleration; flow.cell_count()];
        let force_density = acceleration.map(|a| 2.0 * a);
        let force_field = vec![force_density; flow.cell_count()];
        let force = match mode {
            0 => BodyForce3d::UniformAcceleration(acceleration),
            1 => BodyForce3d::Acceleration(&acceleration_field),
            2 => BodyForce3d::UniformForceDensity(force_density),
            _ => BodyForce3d::ForceDensity(&force_field),
        };
        let report = flow
            .step_with_forcing(
                Time::new(0.002),
                StepControl::default(),
                force,
                control(),
                &mut workspace,
                || false,
            )
            .unwrap();
        balance(&report);
        let dt = report.elapsed.value();
        let mass = report.before.mass();
        let expected_work: f64 = (0..3)
            .map(|a| {
                mass * acceleration[a] * dt * (initial.velocity[a] + 0.5 * acceleration[a] * dt)
            })
            .sum();
        close(report.source_change.work, expected_work, 3e-13);
        for axis in 0..3 {
            close(
                report.source_change.impulse[axis],
                mass * acceleration[axis] * dt,
                3e-13,
            );
        }
        for i in 0..flow.cell_count() {
            let q = flow.cell(i).unwrap();
            let p = q.primitive(flow.gas(), || false).unwrap();
            for species in 0..3 {
                close(
                    q.partial_densities[species],
                    initial_species[species],
                    3e-13,
                );
            }
            assert_eq!(q.partial_densities[2], 0.0);
            close(p.temperature, initial.temperature, 3e-11);
            close(p.pressure, initial.pressure, 3e-11);
            for axis in 0..3 {
                close(
                    p.velocity[axis],
                    initial.velocity[axis] + acceleration[axis] * dt,
                    3e-13,
                );
            }
        }
        if let Some(expected) = &reference {
            assert_eq!(flow.states(), expected);
        } else {
            reference = Some(flow.states().to_vec());
        }
    }
}

#[test]
fn g1_acceleration_and_force_density_use_distinct_local_mass() {
    for density_source in [false, true] {
        let mut flow = field(0.0, true);
        let mut workspace = work(&flow);
        let original = flow.states().to_vec();
        let force = if density_source {
            BodyForce3d::UniformForceDensity([0.4, 0.0, 0.0])
        } else {
            BodyForce3d::UniformAcceleration([0.4, 0.0, 0.0])
        };
        let report = flow
            .step_with_forcing(
                Time::new(0.001),
                StepControl::default(),
                force,
                control(),
                &mut workspace,
                || false,
            )
            .unwrap();
        balance(&report);
        for (i, before) in original.chunks_exact(7).enumerate() {
            let p0 = cell(before, 3).primitive(flow.gas(), || false).unwrap();
            let q = flow.cell(i).unwrap();
            let p = q.primitive(flow.gas(), || false).unwrap();
            let acceleration = if density_source {
                0.4 / p0.density
            } else {
                0.4
            };
            close(p.velocity[0], acceleration * report.elapsed.value(), 3e-12);
            close(p.temperature, p0.temperature, 3e-11);
            close(p.pressure, p0.pressure, 3e-11);
            for k in 0..3 {
                close(q.partial_densities[k], before[k], 3e-12);
            }
        }
    }
}

#[test]
fn g1_force_aware_retry_restarts_both_kicks_from_accepted_composition() {
    let mut flow = field(0.0, false);
    let mut workspace = work(&flow);
    let initial = flow
        .cell(0)
        .unwrap()
        .primitive(flow.gas(), || false)
        .unwrap();
    let report = flow
        .step_with_forcing(
            Time::new(0.01),
            StepControl::default(),
            BodyForce3d::UniformAcceleration([100.0, 0.0, 0.0]),
            control(),
            &mut workspace,
            || false,
        )
        .unwrap();
    assert!(report.halvings > 0);
    assert!(report.elapsed.value() < 0.01);
    balance(&report);
    close(
        report.source_change.impulse[0],
        report.before.mass() * 100.0 * report.elapsed.value(),
        3e-12,
    );
    for i in 0..flow.cell_count() {
        let q = flow.cell(i).unwrap();
        let p = q.primitive(flow.gas(), || false).unwrap();
        close(
            p.velocity[0],
            initial.velocity[0] + 100.0 * report.elapsed.value(),
            3e-12,
        );
        close(p.temperature, initial.temperature, 3e-11);
        assert_eq!(q.partial_densities[2], 0.0);
        control()
            .admit_velocity::<Error>(p.velocity, p.sound_speed)
            .unwrap();
    }
}

#[test]
fn g0_zero_force_preserves_original_step_and_workspace_capacity() {
    let mut actual = field(0.0, true);
    let mut expected = field(0.0, true);
    let mut workspace = work(&actual);
    let bytes = workspace.workspace_bytes();
    let reference = expected
        .step(Time::new(0.002), StepControl::default(), || false)
        .unwrap();
    let report = actual
        .step_with_forcing(
            Time::new(0.002),
            StepControl::default(),
            BodyForce3d::UniformAcceleration([0.0; 3]),
            control(),
            &mut workspace,
            || false,
        )
        .unwrap();
    assert_eq!(actual.states(), expected.states());
    assert_eq!(actual.boundary_impulses, expected.boundary_impulses);
    assert_eq!(report, reference);
    assert_eq!(report.source_change, ForceChange3d::default());
    assert_eq!(workspace.workspace_bytes(), bytes);
}

#[test]
fn g4_force_budget_endpoint_mach_and_last_checkpoint_refuse_atomically() {
    let mut flow = field(0.0, false);
    let bytes = flow.required_forcing_workspace_bytes().unwrap();
    assert!(flow.forcing_workspace(bytes - 1).is_err());
    let mut workspace = work(&flow);
    let before = flow.states().to_vec();
    let before_time = flow.time();
    let before_loads = flow.boundary_impulses.clone();
    let no_retries = StepControl {
        maximum_halvings: 0,
        ..StepControl::default()
    };
    assert!(
        flow.step_with_forcing(
            Time::new(0.01),
            no_retries,
            BodyForce3d::UniformAcceleration([1e5, 0.0, 0.0]),
            control(),
            &mut workspace,
            || false,
        )
        .is_err()
    );
    let malformed = vec![[0.0; 3]; flow.cell_count() - 1];
    assert!(
        flow.step_with_forcing(
            Time::new(0.001),
            StepControl::default(),
            BodyForce3d::ForceDensity(&malformed),
            control(),
            &mut workspace,
            || false,
        )
        .is_err()
    );
    assert_eq!(flow.states(), before);
    assert_eq!(flow.time(), before_time);
    assert_eq!(flow.boundary_impulses, before_loads);
    assert_eq!(flow.last_step(), None);

    let force = BodyForce3d::UniformAcceleration([0.2, -0.1, 0.05]);
    let mut reference = field(0.0, false);
    let mut other_work = work(&reference);
    let mut polls = 0;
    let expected = reference
        .step_with_forcing(
            Time::new(0.002),
            StepControl::default(),
            force,
            control(),
            &mut other_work,
            || {
                polls += 1;
                false
            },
        )
        .unwrap();
    let mut seen = 0;
    assert!(matches!(
        flow.step_with_forcing(
            Time::new(0.002),
            StepControl::default(),
            force,
            control(),
            &mut workspace,
            || {
                seen += 1;
                seen == polls
            },
        ),
        Err(Error::Cancelled)
    ));
    assert_eq!(flow.states(), before);
    assert_eq!(flow.time(), before_time);
    assert_eq!(flow.boundary_impulses, before_loads);
    assert_eq!(flow.last_step(), None);
    let retry = flow
        .step_with_forcing(
            Time::new(0.002),
            StepControl::default(),
            force,
            control(),
            &mut workspace,
            || false,
        )
        .unwrap();
    assert_eq!(retry, expected);
    assert_eq!(flow.states(), reference.states());
}

#[cfg(feature = "experimental-sgs")]
#[test]
fn g3_forced_reacting_channel_keeps_molecular_sgs_and_chemical_energy_owners() {
    use super::super::molecular::sgs::{
        ChannelDynamicComponentsConfig, ChannelDynamicScalarConfig, ChannelDynamicSgsConfig,
        ChannelDynamicThermalSubgridConfig, MixtureChannelDynamicSpeciesSubgridConfig, SgsUse,
    };
    use super::super::molecular::{DissipationControl, MolecularTransport3d, MolecularWall3d};
    use super::super::reactive::{
        CellChemistry, IdealGasKineticsV1, IsochoricKineticsSolve, ReactingEuler3d,
        ReactiveStepControl,
    };
    use fs_gas::multicomponent::dissipation::{
        ConstantTransport, SpatialTransport, TransportCoefficients, TransportDomain,
        WallThermalCondition,
    };
    use fs_qty::{DynViscosity, Temperature};
    use fs_thermochem::kinetics::{
        ArrheniusRateV1, ElementaryReactionV1, ReactionProgressRateQuantityV1,
    };
    use fs_thermochem::transport::kinetic::ThermalConductivityQuantityV1;
    use fs_thermochem::{MolarEnergyQuantityV1, ReactionId, SpeciesId};
    let gas = gas(vec![
        source("A", R, 3.5, 0.0, 0.0),
        source("B", R, 3.5, 0.0, -4.0 * R),
    ]);
    let reaction = ElementaryReactionV1::new(
        ReactionId::new("forced-channel-isomerization").unwrap(),
        vec![(SpeciesId::new("A").unwrap(), 1)],
        vec![(SpeciesId::new("B").unwrap(), 1)],
        ArrheniusRateV1::new(
            Temperature::new(1.0),
            ReactionProgressRateQuantityV1::new(1.0),
            0.0,
            MolarEnergyQuantityV1::new(0.0),
        )
        .unwrap(),
    )
    .unwrap();
    let mechanism = IdealGasKineticsV1::new(
        gas.source().clone(),
        "forced-mixture-source",
        1,
        vec![reaction],
        1e-12,
    )
    .unwrap();
    let grid = CartesianGrid::new(5, 4, 5, 0.2, 0.2, 0.2).unwrap();
    let mut states = vec![0.0; grid.len() * 6];
    for (i, row) in states.chunks_exact_mut(6).enumerate() {
        let [x, y, z] = grid.cell_center(i).unwrap();
        let a = 0.5 + 0.1 * (6.0 * x + 3.0 * y).sin();
        let envelope = 4.0 * z * (1.0 - z);
        state_into(
            &gas,
            &[a, 1.0 - a],
            [
                0.02 * envelope * (7.0 * y).sin(),
                0.01 * envelope * (5.0 * x).cos(),
                0.0,
            ],
            1.0 + 0.05 * (4.0 * x + z).cos(),
            row,
        )
        .unwrap();
    }
    let boundaries = Boundaries3d {
        z: AxisBoundaries::Ends {
            lower: BoundaryCondition::Reflecting,
            upper: BoundaryCondition::Reflecting,
        },
        ..PERIODIC
    };
    let bytes = MixtureEuler3d::required_workspace_bytes(&gas, &grid).unwrap();
    let flow = MixtureEuler3d::new(gas.clone(), grid, states, boundaries, bytes, || false).unwrap();
    let molecular = ConstantTransport::new(
        gas.source().clone(),
        "forced-channel-molecular",
        1,
        TransportDomain {
            minimum_temperature_k: 0.001,
            maximum_temperature_k: 50.0,
            minimum_pressure_pa: 0.001,
            maximum_pressure_pa: 100.0,
        },
        TransportCoefficients {
            dynamic_viscosity: DynViscosity::new(0.02),
            bulk_viscosity: DynViscosity::new(0.003),
            thermal_conductivity: ThermalConductivityQuantityV1::new(0.03),
        },
        vec![0.0, 0.01, 0.01, 0.0],
    )
    .unwrap();
    let transport = SpatialTransport::new(gas.clone(), molecular, || false).unwrap();
    let config = MixtureChannelDynamicSpeciesSubgridConfig {
        transport: ChannelDynamicThermalSubgridConfig {
            momentum: ChannelDynamicSgsConfig {
                maximum_coefficient_squared: 1000.0,
                maximum_eddy_to_molecular: 1000.0,
                usage: SgsUse::Exploration,
            },
            heat: ChannelDynamicScalarConfig {
                maximum_coefficient: 1000.0,
                maximum_diffusivity_to_molecular: 1000.0,
                usage: SgsUse::Exploration,
            },
        },
        species: ChannelDynamicComponentsConfig {
            maximum_coefficient: 1000.0,
            usage: SgsUse::Exploration,
        },
        maximum_species_diffusivity_m2_s: 10.0,
    };
    let wall = MolecularWall3d::NoSlip {
        velocity: [0.0; 3],
        thermal: WallThermalCondition::Adiabatic,
    };
    let bytes = MolecularTransport3d::required_workspace_bytes_with_channel_dynamic_species_sgs(
        &flow, &transport,
    )
    .unwrap();
    let mut operator = MolecularTransport3d::with_channel_dynamic_species_sgs(
        &flow,
        transport,
        config,
        [wall; 2],
        bytes,
        || false,
    )
    .unwrap();
    let chemistry = IsochoricKineticsSolve {
        minimum_temperature_k: 0.001,
        maximum_temperature_k: 50.0,
        energy_tolerance_j: 1e-11,
        relative_tolerance: 2e-6,
        absolute_amount_tolerance_mol: 1e-13,
        entropy_tolerance_j_k: 1e-11,
        initial_step_s: 1e-5,
        maximum_step_s: 0.01,
        maximum_steps: 10000,
        ..IsochoricKineticsSolve::default()
    };
    let reactor = CellChemistry::new(gas, mechanism, chemistry, || false).unwrap();
    let bytes = ReactingEuler3d::required_workspace_bytes(&flow, &reactor).unwrap();
    let mut reacting = ReactingEuler3d::with_chemistry(flow, reactor, bytes, || false).unwrap();
    let mut workspace = work(reacting.flow());
    let complete_bytes =
        reacting.workspace_bytes() + operator.workspace_bytes() + workspace.workspace_bytes();
    let report = reacting
        .step_with_dissipation_and_forcing(
            &mut operator,
            Time::new(1e-4),
            ReactiveStepControl::default(),
            DissipationControl::default(),
            complete_bytes,
            BodyForce3d::UniformForceDensity([0.05, 0.0, 0.0]),
            control(),
            &mut workspace,
            || false,
        )
        .unwrap();
    assert!(report.dissipation.is_some());
    assert!(report.chemistry.species_change_kg[0] < 0.0);
    assert_eq!(report.source_change, report.transport.source_change);
    close(
        report.source_change.impulse[0],
        0.05 * grid.cell_volume() * grid.len() as f64 * report.elapsed.value(),
        3e-12,
    );
    for residual in report.conservation_residual.momentum {
        close(residual, 0.0, 3e-11);
    }
    for residual in report.conservation_residual.species_masses {
        close(residual, 0.0, 3e-11);
    }
    close(report.conservation_residual.energy, 0.0, 3e-10);
    close(report.absolute_energy_residual_j, 0.0, 3e-10);
    assert_eq!(reacting.flow().time(), report.elapsed);
    assert_eq!(
        operator.last_step(),
        report.dissipation.as_ref().map(|d| &d.second)
    );
}
