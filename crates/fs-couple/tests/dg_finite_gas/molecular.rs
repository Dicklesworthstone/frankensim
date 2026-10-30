//! Finite backgrounds pay the original spatial diffusion and molecular work.
use super::*;
use fs_couple::dg_port::DgPortGasUpdate;
use fs_couple::dg_port::finite_gas::transport::{
    FiniteGasDgTransportError, FiniteGasDgTransportInput, FiniteGasDgTransportSolve,
    advance_finite_gas_port_transport,
};
use fs_duct::dg::{DgWallResolution, PortInput};
use fs_material::relaxation::{
    InternalReservoir, RelaxationAbsorption, RelaxationProcess, RelaxationRateLaw,
    SpeciesRelaxation, VibrationalMode,
};
use fs_material::thermochemical_gas::energy::parcel::transfer::{
    nozzle::network::snapshot::SnapshotNozzleLink,
    species::network::{
        SpeciesNetworkSolve,
        transport::{self as material_transport, GasTransportPaths, GasTransportSolve},
    },
};
use fs_thermochem::SpeciesId;
use std::collections::BTreeMap;

const ERROR: f64 = 0.02;

fn solve(molecular: bool) -> FiniteGasDgSolve {
    FiniteGasDgSolve {
        loss: if molecular {
            FiniteGasDgLoss::UniformMolecularDiffusion {
                maximum_relative_error: ERROR,
            }
        } else {
            FiniteGasDgLoss::SpatialViscothermal {
                maximum_relative_error: ERROR,
            }
        },
        ..controls()
    }
}

fn warm(molecular: bool, fraction: f64) -> (DgPortLoad, ClosedGasParcel) {
    let (template, gas) = fixture_with_composition(false, composition(fraction));
    // Retain the authored spatial quadrature, even for this uniform cylinder.
    let mut duct = template
        .duct()
        .clone()
        .with_spatial_viscothermal(DgWallResolution {
            poles_per_field: 8,
            quadrature_points: 3,
            band_rad_s: [50.0, 2000.0].map(|f| core::f64::consts::TAU * f),
            maximum_relative_error: ERROR,
        })
        .unwrap();
    if molecular {
        let law = RelaxationAbsorption::new(
            ["A", "B"]
                .into_iter()
                .enumerate()
                .map(|(i, name)| SpeciesRelaxation {
                    species: SpeciesId::new(name).unwrap(),
                    processes: vec![RelaxationProcess {
                        reservoir: InternalReservoir::Vibrational(vec![VibrationalMode {
                            theta_k: 900.0 + 600.0 * i as f64,
                            degeneracy: 1,
                        }]),
                        rate: RelaxationRateLaw::PartnerIndependent {
                            p_tau_ref_pa_s: 3.0 + 7.0 * i as f64,
                            t_ref_k: 300.0,
                            exponent: 1.2,
                            t_min_k: 200.0,
                            t_max_k: 500.0,
                        },
                    }],
                })
                .collect(),
        )
        .unwrap();
        duct = duct.with_molecular_relaxation_field(law).unwrap();
    }
    let mut load = DgPortLoad::new_with_terminal(
        duct,
        DT,
        1e-11,
        MAX_DG_PROPAGATOR_DOFS,
        template.terminal().unwrap().clone(),
    )
    .unwrap()
    .with_terminal_channels(vec![*template.terminal_channels()[0].input()], CELLS)
    .unwrap();
    for _ in 0..24 {
        load.preview_step(FLOW).unwrap().accept();
    }
    assert!(load.cumulative_dissipated_j() > 0.0);
    if molecular {
        assert!(
            load.duct()
                .molecular_dissipated_power_w(load.state())
                .unwrap()
                > 0.0
        );
    }
    (load, gas)
}

fn native_reference(
    live: &DgPortLoad,
    marched: &DgPortLoad,
    report: &DgPortGasUpdate,
    molecular: bool,
) {
    let mut duct = marched.duct().clone();
    let mut state = marched.state().clone();
    marched.propagator().sync_nodal(&mut state).unwrap();
    let field = live.duct().gas_field().unwrap();
    let native = if molecular {
        duct.rebind_uniform_molecular_gas_field(&mut state, field, ERROR, 1e-5, &mut || false)
            .unwrap()
    } else {
        duct.rebind_spatial_viscothermal_gas_field(&mut state, field, ERROR, 1e-5, &mut || false)
            .unwrap()
    };
    assert_eq!(report.duct, native);
    assert_eq!(live.state(), &state);
    assert_eq!(live.duct().operator_identity(), duct.operator_identity());
    assert_eq!(live.duct().total_dofs(), marched.duct().total_dofs());
    assert_eq!(
        live.duct().molecular_relaxation_law(),
        marched.duct().molecular_relaxation_law()
    );
    assert_eq!(
        live.duct()
            .spatial_viscothermal_admission()
            .unwrap()
            .resolution,
        marched
            .duct()
            .spatial_viscothermal_admission()
            .unwrap()
            .resolution
    );
    assert_eq!(
        live.terminal().unwrap().pole_states(),
        marched.terminal().unwrap().pole_states()
    );
    assert_eq!(live.last_terminal_frames(), marched.last_terminal_frames());
    assert_eq!(
        live.last_terminal_channel_frames(),
        marched.last_terminal_channel_frames()
    );
    let old = marched.terminal_channels()[0].channel().unwrap();
    let new = live.terminal_channels()[0].channel().unwrap();
    assert_eq!(new.last_step(), old.last_step());
    let cell = channel_volume(false) / CELLS as f64;
    let mut kinetic_work = 0.0;
    for (&a, &b) in old.velocities_m_s().iter().zip(new.velocities_m_s()) {
        near(
            new.gas().density * b,
            old.gas().density * a,
            (old.gas().density * a).abs(),
        );
        kinetic_work += 0.5 * cell * (new.gas().density * b * b - old.gas().density * a * a);
    }
    near(
        report.parameter_work_j,
        native.parameter_work_j + kinetic_work,
        report.energy_before_j + report.energy_after_j,
    );
    assert_eq!(new.gas(), field.uniform_state().unwrap().state());
}

#[test]
fn g1_finite_heat_pays_native_spatial_and_molecular_work_without_erasing_memory() {
    for molecular in [false, true] {
        let (mut live, mut gas) = warm(molecular, 0.8);
        let (original, old_gas) = (live.clone(), gas.clone());
        let mut marched = live.clone();
        let frozen = marched.preview_step(FLOW).unwrap().accept();
        // Unchanged coefficients still pass native memory admission. Their
        // continuing state and next sample equal the original acoustic owner.
        let (mut identity, mut identity_gas) = (original.clone(), old_gas.clone());
        let identity_step = identity
            .step_with_finite_gas_controlled(
                &mut identity_gas,
                FLOW,
                FiniteGasDgSolve {
                    work_bracket_j: [0.0, 0.0],
                    ..solve(molecular)
                },
                |old, _, work, _| {
                    assert_eq!(work, 0.0);
                    Ok::<_, GasEnergyError>(FiniteGasDgBackground {
                        gas: old.clone(),
                        supplied_heat_j: 0.0,
                        transported_enthalpy_j: 0.0,
                        source: (),
                    })
                },
                &mut || false,
            )
            .unwrap();
        native_reference(&identity, &marched, &identity_step.refresh, molecular);
        assert_eq!(identity_gas, old_gas);
        assert_eq!(identity_step.refresh.parameter_work_j, 0.0);
        assert_eq!(
            identity.duct().operator_identity(),
            original.duct().operator_identity()
        );
        let mut same_sample = marched.clone();
        assert_eq!(
            identity.preview_step(0.0).unwrap().accept(),
            same_sample.preview_step(0.0).unwrap().accept()
        );
        let invalid = live.step_with_finite_gas_heat_controlled(
            &mut gas,
            FLOW,
            0.0,
            gas_solve(),
            controls(),
            &mut || false,
        );
        assert!(matches!(invalid, Err(FiniteGasDgError::Acoustic(_))));
        assert_eq!(gas, old_gas);
        unchanged(&live, &original);
        let report = live
            .step_with_finite_gas_heat_controlled(
                &mut gas,
                FLOW,
                HEAT,
                gas_solve(),
                solve(molecular),
                &mut || false,
            )
            .unwrap();
        native_reference(&live, &marched, &report.refresh, molecular);
        assert_eq!(live.state().step_count(), original.state().step_count() + 1);
        assert_eq!(gas.mass_kg(), old_gas.mass_kg());
        assert_eq!(gas.volume_m3(), volume(false));
        assert_eq!(
            live.duct().gas_field().unwrap().uniform_state(),
            Some(gas.state())
        );
        assert_eq!(
            report.acoustic.ledger.port_work_step_j,
            frozen.ledger.port_work_step_j
        );
        // Original source polynomial Cp/R = .8*3.5 + .2*4.5 + .001*T.
        let t = gas.state().state().temperature;
        let caloric = gas.mass_kg()
            * gas.state().state().specific_gas_constant
            * (2.7 * (t - 300.0) + 0.0005 * (t * t - 300.0 * 300.0));
        assert!(
            (caloric + report.refresh.parameter_work_j - HEAT).abs()
                < controls().energy_tolerance_j
        );
        assert!(report.refresh.parameter_work_j.abs() > 1e-13);
        let inlet = FLOW
            * frozen
                .port_outputs
                .iter()
                .find(|p| p.id == live.inlet())
                .unwrap()
                .delta_x;
        let exterior: f64 = live
            .last_terminal_frames()
            .iter()
            .map(|f| f.load_work_j)
            .sum();
        let residual = caloric + live.stored_energy_j() - original.stored_energy_j()
            + live.cumulative_dissipated_j()
            - original.cumulative_dissipated_j()
            + exterior
            - inlet
            - HEAT;
        assert!(
            residual.abs() < controls().energy_tolerance_j,
            "{residual:e}"
        );
        assert!(report.work_residual_j.abs() <= controls().work_tolerance_j);
        if molecular {
            // A pressure-only compression scaling misses the live entropy correction.
            let bulk = |load: &DgPortLoad| {
                load.duct().gas().density
                    * load
                        .duct()
                        .molecular_relaxation()
                        .unwrap()
                        .frozen_sound_speed_m_s()
                        .powi(2)
            };
            let ratio = bulk(&live) / bulk(&marched);
            let correction: f64 = live
                .state()
                .pressure()
                .iter()
                .zip(marched.state().pressure())
                .map(|(p, old)| (p - old * ratio).abs())
                .sum();
            assert!(correction > 1e-8);
        }
        let mut difference = 0.0;
        for _ in 0..8 {
            let updated = live.preview_step(0.0).unwrap().accept();
            let reference = marched.preview_step(0.0).unwrap().accept();
            difference += (updated.port_outputs[0].conjugate_rate
                - reference.port_outputs[0].conjugate_rate)
                .abs();
            assert_eq!(
                updated.ledger.cumulative_parameter_work_j,
                report.refresh.parameter_work_j
            );
        }
        assert!(difference > 1e-8);
    }
}

type Parcels = BTreeMap<u64, ClosedGasParcel>;
type Loads = BTreeMap<u64, DgPortLoad>;
type Inputs = BTreeMap<u64, FiniteGasDgTransportInput>;

fn network() -> (
    Parcels,
    Loads,
    Inputs,
    [SnapshotNozzleLink; 1],
    FiniteGasDgTransportSolve,
) {
    // Exact pure-species fractions distinguish admitted transport from changing
    // the identity of an already-live molecular reservoir.
    let (load, gas) = warm(true, 1.0);
    let reservoir = ClosedGasParcel::new(
        gas.model().clone(),
        point(300.0, 180_000.0),
        6.0 * volume(false),
    )
    .unwrap();
    let parcels = BTreeMap::from([(10, gas), (90, reservoir)]);
    let nodes = BTreeMap::from([(10, load)]);
    let inputs = parcels
        .keys()
        .map(|&id| {
            (
                id,
                FiniteGasDgTransportInput {
                    gas: gas_solve(),
                    supplied_heat_j: 0.0,
                    acoustic: (id == 10).then_some(solve(true)),
                    ports: nodes.get(&id).map_or_else(Vec::new, |load| {
                        vec![PortInput {
                            id: load.inlet(),
                            effort: 0.0,
                        }]
                    }),
                },
            )
        })
        .collect();
    let links = [SnapshotNozzleLink {
        id: 7,
        parcels: [90, 10],
        nozzle: GasNozzle {
            area_m2: 1e-8,
            discharge_coefficient: 0.8,
        },
        nozzle_solve: gas_solve(),
    }];
    let controls = FiniteGasDgTransportSolve {
        transport: GasTransportSolve {
            exchange: SpeciesNetworkSolve {
                maximum_species_fraction: 0.05,
                mass_tolerance_kg: 1e-16,
                energy_tolerance_j: 1e-11,
            },
            entropy_tolerance_j_k: 1e-12,
        },
        energy_tolerance_j: 5e-11,
    };
    (parcels, nodes, inputs, links, controls)
}

#[test]
fn g3_equal_composition_transport_pays_molecular_and_channel_work_once() {
    let (mut parcels, mut nodes, inputs, links, controls) = network();
    let (old_gas, old_load) = (parcels.clone(), nodes[&10].clone());
    let endpoints = inputs.iter().map(|(&id, input)| (id, input.gas)).collect();
    let native = material_transport::step_candidate_with_paths(
        &parcels,
        GasTransportPaths {
            nozzles: &links,
            ..Default::default()
        },
        DT,
        &endpoints,
        controls.transport,
        || false,
    )
    .unwrap();
    let mut marched = old_load.clone();
    marched.preview_step(0.0).unwrap().accept();
    let report = advance_finite_gas_port_transport(
        &mut parcels,
        &mut nodes,
        &inputs,
        &links,
        &[],
        DT,
        controls,
        &mut || false,
    )
    .unwrap();
    assert_eq!(report.transport, *native.step());
    let transferred = report.transport.nozzles[0].transferred_mass_kg;
    assert!(transferred > 0.0);
    assert!((parcels[&10].mass_kg() - old_gas[&10].mass_kg() - transferred).abs() < 1e-16);
    assert!((old_gas[&90].mass_kg() - parcels[&90].mass_kg() - transferred).abs() < 1e-16);
    for gas in parcels.values() {
        assert_eq!(
            gas.model().thermodynamics().mole_composition(),
            &composition(1.0)
        );
    }
    assert!(
        report
            .transport
            .exchange
            .fluxes
            .iter()
            .all(|flux| flux.transfer.species == SpeciesId::new("A").unwrap())
    );
    let live = &nodes[&10];
    native_reference(live, &marched, &report.acoustics[&10].refresh, true);
    let gas_change: f64 = parcels
        .iter()
        .map(|(id, gas)| energy_change(&old_gas[id], gas))
        .sum();
    let exterior: f64 = live
        .last_terminal_frames()
        .iter()
        .map(|frame| frame.load_work_j)
        .sum();
    let residual = gas_change + live.stored_energy_j() - old_load.stored_energy_j()
        + live.cumulative_dissipated_j()
        - old_load.cumulative_dissipated_j()
        + exterior;
    assert!(
        residual.abs() <= controls.energy_tolerance_j,
        "{residual:e}"
    );
    assert!(report.combined_energy_residual_j.abs() <= controls.energy_tolerance_j);
    assert!(
        report.acoustics[&10].refresh.terminal_channels[0]
            .transfer
            .unwrap()
            .parameter_work_j
            < 0.0
    );
}

#[test]
fn g4_composition_refusal_invalid_accuracy_and_final_cancel_preserve_every_owner() {
    let (mut parcels, mut nodes, mut inputs, links, controls) = network();
    let (original, original_load) = (parcels.clone(), nodes[&10].clone());
    for maximum_relative_error in [0.0, f64::NAN, 0.51] {
        inputs.get_mut(&10).unwrap().acoustic.as_mut().unwrap().loss =
            FiniteGasDgLoss::UniformMolecularDiffusion {
                maximum_relative_error,
            };
        let result = advance_finite_gas_port_transport(
            &mut parcels,
            &mut nodes,
            &inputs,
            &links,
            &[],
            DT,
            controls,
            &mut || false,
        );
        assert!(matches!(
            result,
            Err(FiniteGasDgTransportError::Acoustic { parcel: 10, .. })
        ));
        assert_eq!(parcels, original);
        unchanged(&nodes[&10], &original_load);
    }
    inputs.get_mut(&10).unwrap().acoustic = Some(solve(true));
    let foreign = ClosedGasParcel::new(
        parcels[&90]
            .model()
            .with_composition(composition(0.1))
            .unwrap(),
        parcels[&90].point().clone(),
        parcels[&90].volume_m3(),
    )
    .unwrap();
    let mut incompatible = parcels.clone();
    incompatible.insert(90, foreign);
    let checkpoint = incompatible.clone();
    let result = advance_finite_gas_port_transport(
        &mut incompatible,
        &mut nodes,
        &inputs,
        &links,
        &[],
        DT,
        controls,
        &mut || false,
    );
    assert!(matches!(
        result,
        Err(FiniteGasDgTransportError::Acoustic {
            parcel: 10,
            error: FiniteGasDgError::Acoustic(_),
        })
    ));
    assert_eq!(incompatible, checkpoint);
    unchanged(&nodes[&10], &original_load);

    let (mut expected_gas, mut expected_nodes) = (parcels.clone(), nodes.clone());
    let mut polls = 0;
    let expected = advance_finite_gas_port_transport(
        &mut expected_gas,
        &mut expected_nodes,
        &inputs,
        &links,
        &[],
        DT,
        controls,
        &mut || {
            polls += 1;
            false
        },
    )
    .unwrap();
    let mut seen = 0;
    let result = advance_finite_gas_port_transport(
        &mut parcels,
        &mut nodes,
        &inputs,
        &links,
        &[],
        DT,
        controls,
        &mut || {
            seen += 1;
            seen == polls
        },
    );
    assert!(matches!(result, Err(FiniteGasDgTransportError::Cancelled)));
    assert_eq!(seen, polls);
    assert_eq!(parcels, original);
    unchanged(&nodes[&10], &original_load);
    let retried = advance_finite_gas_port_transport(
        &mut parcels,
        &mut nodes,
        &inputs,
        &links,
        &[],
        DT,
        controls,
        &mut || false,
    )
    .unwrap();
    assert_eq!(parcels, expected_gas);
    unchanged(&nodes[&10], &expected_nodes[&10]);
    assert_eq!(retried.transport, expected.transport);
    assert_eq!(
        retried.acoustics[&10].acoustic,
        expected.acoustics[&10].acoustic
    );
    assert_eq!(
        retried.acoustics[&10].refresh,
        expected.acoustics[&10].refresh
    );
    for _ in 0..8 {
        let actual = nodes
            .get_mut(&10)
            .unwrap()
            .preview_step(0.0)
            .unwrap()
            .accept();
        let clean = expected_nodes
            .get_mut(&10)
            .unwrap()
            .preview_step(0.0)
            .unwrap()
            .accept();
        assert_eq!(actual, clean);
    }
}
