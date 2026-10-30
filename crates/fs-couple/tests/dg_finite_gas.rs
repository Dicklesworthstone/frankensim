//! Finite thermochemical storage drives the complete continuing DG load.
//! Synthetic source laws test implementation and conservation, not measured gases.

use fs_couple::dg_port::finite_gas::{
    FiniteGasDgBackground, FiniteGasDgError, FiniteGasDgLoss, FiniteGasDgPortStep, FiniteGasDgSolve,
};
use fs_couple::dg_port::{DgPortLoad, DgTerminalChannelInput, DgTerminalChannelKind};
use fs_duct::dg::{DgBoundary, DgDuct, DgMesh, MAX_DG_PROPAGATOR_DOFS, PortId};
use fs_duct::{Duct, Segment};
use fs_matdb::QueryPoint;
use fs_material::gas_field::{AxialGasField, PrimitiveGasState};
use fs_material::gas_species::input::ThermochemicalPhase;
use fs_material::thermochemical_gas::energy::parcel::ClosedGasParcel;
use fs_material::thermochemical_gas::energy::parcel::transfer::GasMassTransferSolve;
use fs_material::thermochemical_gas::energy::parcel::transfer::nozzle::{
    GasNozzle, GasNozzleStepSolve,
};
use fs_material::thermochemical_gas::energy::{GasEnergyError, GasEnergySolve};
use fs_material::thermochemical_gas::packs::from_packs;
use fs_math::c64::C64;
use fs_qty::QuantitySpec;
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_thermochem::{Composition, CompositionBasis};
use fs_vfit::mimo::MimoModel;
use fs_vfit::mimo::load::{MimoLoadLimits, PreparedMimoLoad};
use fs_vfit::model::PoleTerm;

#[path = "../src/render/schedule/reed/field_input/tests/support.rs"]
mod gas_support;

#[path = "dg_finite_gas/transport.rs"]
mod transport;

const DT: f64 = 1.0 / 96_000.0;
const RADIUS: f64 = 0.001;
const LENGTH: f64 = 0.05;
const CHANNEL_LENGTH: f64 = 0.004;
const CHANNEL_RADIUS: f64 = 0.0002;
const CHANNEL_WIDTH: f64 = 0.004;
const CHANNEL_HEIGHT: f64 = 0.0001;
const CELLS: usize = 12;
const FLOW: f64 = 1e-6;
const HEAT: f64 = 0.002;

fn composition(fraction: f64) -> Composition {
    Composition::new(
        CompositionBasis::MoleFraction,
        vec![fraction, 1.0 - fraction],
    )
    .unwrap()
}

fn point(temperature: f64, pressure: f64) -> QueryPoint {
    QueryPoint::new()
        .with_quantity(
            "temperature",
            QuantitySpec::semantic(SemanticType::new(
                QuantityKind::AbsoluteTemperature,
                ValueForm::Static,
            )),
            temperature,
        )
        .unwrap()
        .with_quantity(
            "pressure",
            QuantitySpec::semantic(SemanticType::new(QuantityKind::Pressure, ValueForm::Static)),
            pressure,
        )
        .unwrap()
}

fn channel_volume(circular: bool) -> f64 {
    if circular {
        core::f64::consts::PI * CHANNEL_RADIUS * CHANNEL_RADIUS * CHANNEL_LENGTH
    } else {
        CHANNEL_WIDTH * CHANNEL_HEIGHT * CHANNEL_LENGTH
    }
}

fn volume(circular: bool) -> f64 {
    core::f64::consts::PI * RADIUS * RADIUS * LENGTH + channel_volume(circular)
}

fn fixture(circular: bool) -> (DgPortLoad, ClosedGasParcel) {
    let a = gas_support::packs("A", 0.028, 3.5, 200.0);
    let b = gas_support::packs("B", 0.044, 4.5, 200.0);
    let model = from_packs(&[a.binding(), b.binding()], composition(0.8), None).unwrap();
    let field = AxialGasField::new(
        vec![0.0, LENGTH],
        101325.0,
        vec![
            PrimitiveGasState {
                temperature_k: 300.0,
                composition: composition(0.8)
            };
            2
        ],
        model.clone(),
        QueryPoint::new(),
        ThermochemicalPhase::FrozenGas,
    )
    .unwrap();
    let gas = ClosedGasParcel::new(model, point(300.0, 101325.0), volume(circular)).unwrap();
    assert_eq!(field.uniform_state(), Some(gas.state()));
    let duct = DgDuct::new_with_gas_field(
        &Duct {
            segments: vec![Segment::Cylinder {
                radius: RADIUS,
                length: LENGTH,
            }],
        },
        &field,
        DgMesh {
            elements: 2,
            order: 1,
        },
        DgBoundary::Port(PortId(0)),
        DgBoundary::Port(PortId(1)),
    )
    .unwrap();
    let terminal = PreparedMimoLoad::new(
        &MimoModel {
            ports: 1,
            poles: vec![PoleTerm::Real {
                pole: -2000.0,
                residue: 0.0,
            }],
            residues: vec![vec![C64::from_re(2e11)]],
            d: vec![1e8],
            e: vec![0.0],
            reciprocal: true,
        },
        DT,
        MimoLoadLimits {
            band_rad_s: [1.0, core::f64::consts::PI / DT],
            sampled_passivity_points: 257,
            max_states: 4,
            max_memory_bytes: 1 << 20,
        },
    )
    .unwrap();
    let kind = if circular {
        DgTerminalChannelKind::Circular {
            radius_m: CHANNEL_RADIUS,
            length_m: CHANNEL_LENGTH,
            radial_cells: CELLS,
            maximum_mach: 0.3,
        }
    } else {
        DgTerminalChannelKind::Slit {
            width_m: CHANNEL_WIDTH,
            height_m: CHANNEL_HEIGHT,
            length_m: CHANNEL_LENGTH,
            transverse_cells: CELLS,
            maximum_mach: 0.3,
        }
    };
    let load = DgPortLoad::new_with_terminal(duct, DT, 1e-11, MAX_DG_PROPAGATOR_DOFS, terminal)
        .unwrap()
        .with_terminal_channels(
            vec![DgTerminalChannelInput {
                port: PortId(1),
                kind,
            }],
            CELLS,
        )
        .unwrap();
    (load, gas)
}

fn controls() -> FiniteGasDgSolve {
    FiniteGasDgSolve {
        work_bracket_j: [-1e-5, 1e-5],
        work_tolerance_j: 1e-14,
        energy_tolerance_j: 2e-11,
        maximum_iterations: 80,
        maximum_temperature_increment_k: 80.0,
        loss: FiniteGasDgLoss::Lossless,
    }
}

fn gas_solve() -> GasEnergySolve {
    GasEnergySolve {
        minimum_temperature_k: 200.0,
        maximum_temperature_k: 900.0,
        energy_tolerance_j_kg: 1e-8,
        maximum_iterations: 128,
    }
}

fn frozen_step(
    load: &mut DgPortLoad,
    gas: &mut ClosedGasParcel,
    flow: f64,
) -> FiniteGasDgPortStep<()> {
    load.step_with_finite_gas_controlled(
        gas,
        flow,
        FiniteGasDgSolve {
            work_bracket_j: [0.0, 0.0],
            ..controls()
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
    .unwrap()
}

fn primed(circular: bool) -> (DgPortLoad, ClosedGasParcel) {
    let (mut load, mut gas) = fixture(circular);
    for _ in 0..24 {
        frozen_step(&mut load, &mut gas, FLOW);
    }
    assert!(load.stored_energy_j() > 0.0);
    assert!(load.terminal_channels()[0].stored_energy_j() > 0.0);
    assert!(
        load.terminal()
            .unwrap()
            .pole_states()
            .iter()
            .any(|p| p.re != 0.0 || p.im != 0.0)
    );
    (load, gas)
}

fn energy_change(old: &ClosedGasParcel, next: &ClosedGasParcel) -> f64 {
    let u0 = old.state().specific_internal_energy_j_kg();
    next.mass_kg() * (next.state().specific_internal_energy_j_kg() - u0)
        + (next.mass_kg() - old.mass_kg()) * u0
}

fn near(actual: f64, expected: f64, scale: f64) {
    assert!(
        (actual - expected).abs() <= 5e-11 * scale.max(1e-30),
        "actual {actual:e}, expected {expected:e}, scale {scale:e}"
    );
}

fn unchanged(load: &DgPortLoad, original: &DgPortLoad) {
    assert_eq!(load.state(), original.state());
    assert_eq!(load.duct().gas_field(), original.duct().gas_field());
    assert_eq!(
        load.duct().operator_identity(),
        original.duct().operator_identity()
    );
    assert_eq!(load.impedance_pa_s_m3(), original.impedance_pa_s_m3());
    assert_eq!(load.propagator().dt(), original.propagator().dt());
    assert_eq!(
        load.port_response().unwrap().history_pressure_pa,
        original.port_response().unwrap().history_pressure_pa
    );
    assert_eq!(
        load.terminal().unwrap().pole_states(),
        original.terminal().unwrap().pole_states()
    );
    assert_eq!(load.last_terminal_frames(), original.last_terminal_frames());
    assert_eq!(
        load.last_terminal_channel_frames(),
        original.last_terminal_channel_frames()
    );
    assert_eq!(load.stored_energy_j(), original.stored_energy_j());
    assert_eq!(
        load.cumulative_dissipated_j(),
        original.cumulative_dissipated_j()
    );
    assert_eq!(
        load.cumulative_parameter_work_j(),
        original.cumulative_parameter_work_j()
    );
    for (current, old) in load
        .terminal_channels()
        .iter()
        .zip(original.terminal_channels())
    {
        assert_eq!(current.input(), old.input());
        assert_eq!(current.gas(), old.gas());
        assert_eq!(
            current.channel().unwrap().velocities_m_s(),
            old.channel().unwrap().velocities_m_s()
        );
        assert_eq!(
            current.channel().unwrap().last_step(),
            old.channel().unwrap().last_step()
        );
        assert_eq!(
            current.cumulative_pressure_work_j(),
            old.cumulative_pressure_work_j()
        );
        assert_eq!(
            current.cumulative_parameter_work_j(),
            old.cumulative_parameter_work_j()
        );
        assert_eq!(
            current.cumulative_viscous_dissipation_j(),
            old.cumulative_viscous_dissipation_j()
        );
        assert_eq!(
            current.cumulative_numerical_dissipation_j(),
            old.cumulative_numerical_dissipation_j()
        );
    }
}

#[test]
fn g1_heat_pays_actual_wave_work_and_retains_terminal_motion_and_memory() {
    for circular in [false, true] {
        let (mut load, mut gas) = primed(circular);
        let (before, old_gas) = (load.clone(), gas.clone());
        let (mut marched, mut frozen_gas) = (load.clone(), gas.clone());
        let frozen = frozen_step(&mut marched, &mut frozen_gas, FLOW);
        assert_eq!(frozen_gas, old_gas);
        let step = load
            .step_with_finite_gas_heat_controlled(
                &mut gas,
                FLOW,
                HEAT,
                gas_solve(),
                controls(),
                &mut || false,
            )
            .unwrap();

        // Independent caloric polynomial and exact pressure-polynomial integral.
        // Both source species have Cp/R = constant + .001*T; their mole-weighted
        // constant is .8*3.5 + .2*4.5. Fixed mass and full physical volume retain rho.
        let rho = old_gas.state().state().density;
        let r = old_gas.state().state().specific_gas_constant;
        let cp0 = 0.8 * 3.5 + 0.2 * 4.5;
        let bulk = |t: f64| rho * r * t * (cp0 + 0.001 * t) / (cp0 - 1.0 + 0.001 * t);
        let area = core::f64::consts::PI * RADIUS * RADIUS;
        let pressure_energy: f64 = marched
            .state()
            .pressure()
            .chunks_exact(2)
            .map(|p| {
                area * (LENGTH / 2.0) / (6.0 * bulk(300.0))
                    * (p[0] * p[0] + p[0] * p[1] + p[1] * p[1])
            })
            .sum();
        let balance = |t: f64| {
            old_gas.mass_kg() * r * ((cp0 - 1.0) * (t - 300.0) + 0.0005 * (t * t - 300.0 * 300.0))
                + pressure_energy * (bulk(t) / bulk(300.0) - 1.0)
                - HEAT
        };
        let (mut lo, mut hi) = (300.0, 380.0);
        for _ in 0..64 {
            let mid = (lo + hi) / 2.0;
            if balance(mid) > 0.0 {
                hi = mid;
            } else {
                lo = mid;
            }
        }
        let temperature = gas.state().state().temperature;
        assert!((temperature - (lo + hi) / 2.0).abs() < 2e-7);
        let expected_work = pressure_energy * (bulk(temperature) / bulk(300.0) - 1.0);
        assert!((step.refresh.parameter_work_j - expected_work).abs() < 3e-13);
        assert!(step.refresh.parameter_work_j > 1e-13);
        let mut heat_only = old_gas.clone();
        heat_only.add_heat(HEAT, gas_solve(), || false).unwrap();
        assert!(heat_only.state().state().temperature > temperature);
        assert_eq!(gas.mass_kg(), old_gas.mass_kg());
        assert_eq!(gas.volume_m3(), volume(circular));
        assert_eq!(
            load.duct().gas_field().unwrap().uniform_state(),
            Some(gas.state())
        );
        assert_eq!(load.terminal_channels()[0].gas(), gas.state().state());
        assert_eq!(load.state().step_count(), before.state().step_count() + 1);
        assert_eq!(
            load.terminal().unwrap().pole_states(),
            marched.terminal().unwrap().pole_states()
        );
        assert_eq!(load.last_terminal_frames(), marched.last_terminal_frames());
        assert_eq!(
            load.last_terminal_channel_frames(),
            marched.last_terminal_channel_frames()
        );
        assert_eq!(
            step.acoustic.ledger.port_work_step_j,
            frozen.acoustic.ledger.port_work_step_j
        );
        for (actual, previous) in step
            .acoustic
            .port_outputs
            .iter()
            .zip(&frozen.acoustic.port_outputs)
        {
            assert_eq!(actual.delta_x, previous.delta_x);
            assert_eq!(
                actual.end_conjugate,
                load.duct()
                    .port_pressures(load.state())
                    .unwrap()
                    .iter()
                    .find(|(id, _)| *id == actual.id)
                    .unwrap()
                    .1
            );
        }
        assert_eq!(step.acoustic.ledger.stored_energy_j, load.stored_energy_j());
        assert_eq!(
            step.acoustic.ledger.cumulative_parameter_work_j,
            load.cumulative_parameter_work_j()
        );
        let inlet_work = FLOW
            * step
                .acoustic
                .port_outputs
                .iter()
                .find(|p| p.id == PortId(0))
                .unwrap()
                .delta_x;
        let exterior_work: f64 = load
            .last_terminal_frames()
            .iter()
            .map(|f| f.load_work_j)
            .sum();
        near(step.terminal_work_j, exterior_work, exterior_work.abs());
        let total = energy_change(&old_gas, &gas) + load.stored_energy_j()
            - before.stored_energy_j()
            + (load.cumulative_dissipated_j() - before.cumulative_dissipated_j())
            + exterior_work
            - inlet_work
            - HEAT;
        assert!(total.abs() <= controls().energy_tolerance_j, "{total:e}");
        assert!(step.work_residual_j.abs() <= controls().work_tolerance_j);
        assert!(step.combined_energy_residual_j.abs() <= controls().energy_tolerance_j);
    }
}

#[test]
fn g1_native_mass_and_species_exchange_pays_bore_and_resolved_channel_work_once() {
    for circular in [false, true] {
        let (mut load, mut gas) = primed(circular);
        let (before, old_gas) = (load.clone(), gas.clone());
        let (mut marched, mut frozen_gas) = (load.clone(), gas.clone());
        frozen_step(&mut marched, &mut frozen_gas, 0.0);
        let reservoir = ClosedGasParcel::new(
            gas.model().with_composition(composition(0.1)).unwrap(),
            point(400.0, 180000.0),
            4.0 * volume(circular),
        )
        .unwrap();
        let flow_solve = GasNozzleStepSolve {
            nozzle: gas_solve(),
            transfer: GasMassTransferSolve {
                donor: gas_solve(),
                receiver: gas_solve(),
                mass_tolerance_kg: 1e-16,
                energy_tolerance_j: 1e-11,
            },
            maximum_donor_fraction: 0.05,
        };
        let step = load
            .step_with_finite_gas_controlled(
                &mut gas,
                0.0,
                controls(),
                |old, dt, work, cancel| {
                    let (proposal, flow) = reservoir.through_nozzle_to(
                        old,
                        GasNozzle {
                            area_m2: 1e-6,
                            discharge_coefficient: 0.8,
                        },
                        dt,
                        flow_solve,
                        &mut *cancel,
                    )?;
                    let (donor, mut receiver, transfer) = proposal.into_parts();
                    receiver.replace_model_with_heat_work_volume(
                        receiver.model().clone(),
                        receiver.volume_m3(),
                        0.0,
                        work,
                        gas_solve(),
                        &mut *cancel,
                    )?;
                    Ok::<_, GasEnergyError>(FiniteGasDgBackground {
                        gas: receiver,
                        supplied_heat_j: 0.0,
                        transported_enthalpy_j: transfer.transferred_enthalpy_j,
                        source: (donor, transfer, flow),
                    })
                },
                &mut || false,
            )
            .unwrap();
        let (donor, transfer, flow) = &step.background;
        assert!(flow.mass_flow_kg_s > 0.0);
        assert_eq!(transfer.transferred_mass_kg, flow.mass_flow_kg_s * DT);
        assert!(gas.mass_kg() > old_gas.mass_kg());
        assert_ne!(
            gas.model().thermodynamics().mole_composition(),
            old_gas.model().thermodynamics().mole_composition()
        );
        let old_channel = marched.terminal_channels()[0].channel().unwrap();
        let new_channel = load.terminal_channels()[0].channel().unwrap();
        let (rho0, rho1) = (old_channel.gas().density, new_channel.gas().density);
        assert!(rho1 > rho0);
        let mut e0 = 0.0;
        let mut e1 = 0.0;
        for (i, (&v0, &v1)) in old_channel
            .velocities_m_s()
            .iter()
            .zip(new_channel.velocities_m_s())
            .enumerate()
        {
            let cell_volume = if circular {
                let inner = CHANNEL_RADIUS * i as f64 / CELLS as f64;
                let outer = CHANNEL_RADIUS * (i + 1) as f64 / CELLS as f64;
                core::f64::consts::PI * CHANNEL_LENGTH * (outer * outer - inner * inner)
            } else {
                channel_volume(false) / CELLS as f64
            };
            let momentum = rho0 * cell_volume * v0;
            near(rho1 * cell_volume * v1, momentum, momentum.abs());
            e0 += 0.5 * rho0 * cell_volume * v0 * v0;
            e1 += 0.5 * rho1 * cell_volume * v1 * v1;
        }
        assert!(e0 > 0.0);
        near(e1, e0 * rho0 / rho1, e0);
        let channel = step.refresh.terminal_channels[0].transfer.unwrap();
        near(channel.kinetic_energy_before_j, e0, e0);
        near(channel.kinetic_energy_after_j, e1, e1);
        near(channel.parameter_work_j, e1 - e0, e0 + e1);
        assert!(channel.parameter_work_j < 0.0);
        assert_eq!(new_channel.last_step(), old_channel.last_step());
        assert_eq!(
            load.terminal().unwrap().pole_states(),
            marched.terminal().unwrap().pole_states()
        );
        near(
            step.refresh.parameter_work_j,
            step.refresh.duct.parameter_work_j + channel.parameter_work_j,
            step.refresh.absolute_parameter_work_j,
        );
        assert_eq!(step.supplied_heat_j, 0.0);
        assert_eq!(step.transported_enthalpy_j, transfer.transferred_enthalpy_j);
        let total = energy_change(&old_gas, &gas)
            + energy_change(&reservoir, donor)
            + load.stored_energy_j()
            - before.stored_energy_j()
            + step.acoustic_dissipation_j
            + step.terminal_work_j;
        assert!(total.abs() <= controls().energy_tolerance_j, "{total:e}");
        assert!(step.work_residual_j.abs() <= controls().work_tolerance_j);
        assert_eq!(load.state().step_count(), before.state().step_count() + 1);
    }
}

#[test]
fn g4_stale_inventory_late_source_refusal_and_final_cancellation_preserve_exact_retry() {
    let (mut load, mut gas) = primed(true);
    let (original, original_gas) = (load.clone(), gas.clone());
    let mut wrong_volume = ClosedGasParcel::new(
        gas.model().clone(),
        gas.point().clone(),
        core::f64::consts::PI * RADIUS * RADIUS * LENGTH,
    )
    .unwrap();
    let mut stale =
        ClosedGasParcel::new(gas.model().clone(), point(301.0, 101325.0), volume(true)).unwrap();
    for bad in [&mut wrong_volume, &mut stale] {
        let before = bad.clone();
        let mut called = false;
        let result = load.step_with_finite_gas_controlled(
            bad,
            FLOW,
            controls(),
            |_, _, _, _| {
                called = true;
                Err::<FiniteGasDgBackground<()>, _>(GasEnergyError::Invalid("must not run"))
            },
            &mut || false,
        );
        assert!(matches!(result, Err(FiniteGasDgError::Invalid(_))));
        assert!(!called);
        assert_eq!(*bad, before);
        unchanged(&load, &original);
    }
    let mut trials = 0;
    let result = load.step_with_finite_gas_controlled(
        &mut gas,
        FLOW,
        controls(),
        |old, _, work, cancel| {
            trials += 1;
            if trials > 1 {
                return Err(GasEnergyError::Invalid("late physical source refusal"));
            }
            let mut proposal = old.clone();
            let source = proposal.replace_model_with_heat_work_volume(
                old.model().clone(),
                old.volume_m3(),
                HEAT,
                work,
                gas_solve(),
                cancel,
            )?;
            Ok(FiniteGasDgBackground {
                gas: proposal,
                supplied_heat_j: HEAT,
                transported_enthalpy_j: 0.0,
                source,
            })
        },
        &mut || false,
    );
    assert!(matches!(result, Err(FiniteGasDgError::Background(_))));
    assert_eq!(trials, 2);
    assert_eq!(gas, original_gas);
    unchanged(&load, &original);
    assert!(
        load.step_with_finite_gas_heat_controlled(
            &mut gas,
            FLOW,
            HEAT,
            gas_solve(),
            FiniteGasDgSolve {
                work_bracket_j: [0.0, 0.0],
                ..controls()
            },
            &mut || false,
        )
        .is_err()
    );
    assert_eq!(gas, original_gas);
    unchanged(&load, &original);
    let (mut expected, mut expected_gas) = (original.clone(), original_gas.clone());
    let mut polls = 0;
    let expected_step = expected
        .step_with_finite_gas_heat_controlled(
            &mut expected_gas,
            FLOW,
            HEAT,
            gas_solve(),
            controls(),
            &mut || {
                polls += 1;
                false
            },
        )
        .unwrap();
    let mut count = 0;
    assert!(matches!(
        load.step_with_finite_gas_heat_controlled(
            &mut gas,
            FLOW,
            HEAT,
            gas_solve(),
            controls(),
            &mut || {
                count += 1;
                count == polls
            },
        ),
        Err(FiniteGasDgError::Cancelled)
    ));
    assert_eq!(count, polls);
    assert_eq!(gas, original_gas);
    unchanged(&load, &original);
    let retry = load
        .step_with_finite_gas_heat_controlled(
            &mut gas,
            FLOW,
            HEAT,
            gas_solve(),
            controls(),
            &mut || false,
        )
        .unwrap();
    assert_eq!(gas, expected_gas);
    assert_eq!(retry.acoustic, expected_step.acoustic);
    assert_eq!(retry.refresh, expected_step.refresh);
    assert_eq!(retry.work_by_gas_j, expected_step.work_by_gas_j);
    unchanged(&load, &expected);
}
