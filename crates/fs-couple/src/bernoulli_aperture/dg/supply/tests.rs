//! Physical component connections; manufactured fixtures do not qualify a reed.
use super::*;
use crate::bernoulli_aperture::BernoulliAperture;
use crate::bernoulli_aperture::dg::LoadedApertureDg;
use crate::bernoulli_aperture::dynamic::{ApertureState, DynamicAperture, DynamicApertureSpec};
use crate::bernoulli_aperture::slit_channel::circular::{CircularChannel, CircularChannelSpec};
use crate::bernoulli_aperture::slit_channel::{SlitChannel, SlitChannelSpec, SlitExit};
use crate::bernoulli_aperture::tube::TubeDrive;
use fs_dcontact::Obstacle;
use fs_duct::dg::{DgBoundary, DgDuct, DgMesh, MAX_DG_PROPAGATOR_DOFS, PortId};
use fs_duct::{Duct, Segment};
use fs_material::gas::GasSpec;

fn gas(pressure: f64) -> GasState {
    GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, pressure).unwrap()
}

fn load(gas: &GasState, dt: f64, length: f64) -> ApertureDgLoad {
    let duct = DgDuct::new(
        &Duct {
            segments: vec![Segment::Cylinder {
                radius: 0.001,
                length,
            }],
        },
        gas,
        DgMesh {
            elements: 1,
            order: 2,
        },
        DgBoundary::Port(PortId(0)),
        DgBoundary::IdealOpen,
    )
    .unwrap()
    .with_viscothermal(8, [50.0, 2000.0], 0.01)
    .unwrap();
    ApertureDgLoad::new(duct, dt, 1e-11, MAX_DG_PROPAGATOR_DOFS).unwrap()
}

fn slit(gas: &GasState, dt: f64, cells: usize) -> SlitChannel {
    SlitChannel::from_gas(
        gas,
        SlitChannelSpec {
            length_m: 0.004,
            width_m: 0.002,
            height_m: 0.0002,
            transverse_cells: cells,
            sample_period_s: dt,
            exit: SlitExit::StaticPressure,
        },
    )
    .unwrap()
}

fn aperture(load: &ApertureDgLoad, impedance: f64) -> DynamicAperture {
    let contact = Obstacle::new(
        vec![-1.0],
        1,
        1,
        vec![0.0],
        vec![1.0],
        0.0,
        2.0,
        "DG supply manufactured no-contact fixture".into(),
    )
    .unwrap();
    DynamicAperture::new_with_propagator(
        DynamicApertureSpec {
            aperture: BernoulliAperture {
                rest_opening_m: 0.0002,
                width_m: 0.001,
                closing_pressure_pa: 100.0,
            },
            mass_kg: 1e-5,
            stiffness_n_m: 5.0,
            damping_ratio: 0.02,
            density_kg_m3: load.duct().inlet_gas().density,
            impedance_pa_s_m3: impedance,
            time_step_s: load.propagator().dt(),
            max_steps: 2048,
        },
        ApertureState {
            opening_m: 0.0002,
            opening_velocity_m_s: 0.0,
        },
        contact,
        BladePropagator::ExactZoh,
    )
    .unwrap()
}

fn system(channel: FixedChannel, tract: bool) -> LoadedApertureDg {
    let dt = channel.sample_period_s();
    let down = load(channel.gas(), dt, 0.05);
    let up = tract.then(|| load(channel.gas(), dt, 0.035));
    let (_, impedance) = inlet_state_with_supply(up.as_ref(), &down, &channel).unwrap();
    LoadedApertureDg::new_loaded_with_supply(aperture(&down, impedance), up, down, channel).unwrap()
}

fn single(model: &LoadedApertureDg) -> &DynamicAperture {
    let LoadedAperture::Single(source) = model.aperture() else {
        unreachable!()
    };
    source
}

#[test]
fn both_channel_geometries_share_current_flow_and_work_and_release_passively() {
    let gas = gas(101325.0);
    let dt = 1.0 / 96000.0;
    let circular: FixedChannel = CircularChannel::from_gas(
        &gas,
        CircularChannelSpec {
            length_m: 0.004,
            radius_m: 0.0002,
            radial_cells: 16,
            sample_period_s: dt,
            exit: SlitExit::StaticPressure,
        },
    )
    .unwrap()
    .into();
    for channel in [slit(&gas, dt, 16).into(), circular] {
        for tract in [false, true] {
            let mut model = system(channel.clone(), tract);
            let mut reversed = false;
            let mut heat = 0.0;
            for step in 0..256 {
                let before = model.stored_energy_j();
                let drive = TubeDrive {
                    upstream_pressure_pa: if step < 48 {
                        0.5
                    } else if step < 96 {
                        -0.5
                    } else {
                        0.0
                    },
                    body_flow_m3_s: if step < 16 { 2e-10 } else { 0.0 },
                };
                let frame = model
                    .step_with_forces(drive, [if step < 24 { 1e-6 } else { 0.0 }, 0.0])
                    .unwrap();
                let supply = frame.supply.unwrap();
                let through = frame.aperture.bore_flow_m3_s - drive.body_flow_m3_s;
                let scale = frame.stored_energy_j.abs()
                    + frame.storage_change_j.abs()
                    + frame.inferred_total_dissipation_j.abs()
                    + frame.upstream_work_j.abs()
                    + frame.body_work_j.abs()
                    + frame.aperture.mechanical_work_j.abs();
                assert!(frame.balance_residual_j().abs() <= 1e-8 * scale.max(1e-24));
                assert!(frame.port_work_residual_j.abs() <= 1e-10 * scale.max(1e-24));
                assert!(
                    (frame.stored_energy_j - model.stored_energy_j()).abs()
                        <= 1e-13 * scale.max(1e-24)
                );
                assert!(
                    (frame.storage_change_j - (frame.stored_energy_j - before)).abs()
                        <= 1e-8 * scale.max(1e-24)
                );
                assert!(
                    (supply.channel.flow_m3_s - through).abs() <= 1e-10 * through.abs().max(1e-24)
                );
                let pu = frame
                    .upstream
                    .as_ref()
                    .map_or(0.0, |up| up.inlet_pressure_pa);
                let face = drive.upstream_pressure_pa + pu - supply.pressure_drop_pa;
                let pressure_work = (face - frame.aperture.bore_pressure_pa) * through * dt;
                assert!(
                    (frame.aperture.pressure_work_j - pressure_work).abs()
                        <= 1e-9 * scale.max(1e-24)
                );
                if let Some(up) = &frame.upstream {
                    assert_eq!(up.upstream_face_pressure_pa, face);
                    assert_eq!(up.inlet_flow_m3_s, -through);
                }
                assert!(supply.channel.viscous_dissipation_j >= 0.0);
                assert!(supply.channel.numerical_dissipation_j >= 0.0);
                assert_eq!(supply.channel.exit_work_j, 0.0);
                if step >= 96 {
                    assert!(frame.stored_energy_j <= before + 1e-8 * scale.max(1e-24));
                }
                reversed |= through < 0.0;
                heat += supply.channel.viscous_dissipation_j;
            }
            assert!(reversed && heat > 0.0);
            assert!(model.supply_channel().unwrap().stored_energy_j() > 0.0);
        }
    }
}

// Independent continuum oracle: integrate the odd sine eigenfunctions of
// rho dv/dt = dp(t)/L + mu d²v/dy², v(0)=v(h)=0. Every measured same-step
// pressure is held over its interval. This checks the coupled startup/reversal
// against the analytic unsteady Poiseuille response, not another FV stepper.
fn analytic_flow_error(rate: usize, cells: usize) -> f64 {
    let gas = gas(101325.0);
    let dt = 1.0 / rate as f64;
    let passage = slit(&gas, dt, cells);
    let spec = *passage.spec();
    let mut model = system(passage.into(), true);
    let mut modes = [0.0; 128];
    let mut squared_error = 0.0;
    let mut squared_flow = 0.0;
    for step in 0..rate / 250 {
        let time = step as f64 * dt;
        let frame = model
            .step(TubeDrive {
                upstream_pressure_pa: if time < 0.0015 {
                    0.5
                } else if time < 0.003 {
                    -0.5
                } else {
                    0.0
                },
                body_flow_m3_s: 0.0,
            })
            .unwrap();
        let supply = frame.supply.unwrap();
        for (index, mode) in modes.iter_mut().enumerate() {
            let npi = (2 * index + 1) as f64 * core::f64::consts::PI;
            let decay = gas.dynamic_viscosity / gas.density * (npi / spec.height_m).powi(2);
            let keep = (-decay * dt).exp();
            *mode = keep * *mode
                + (1.0 - keep) / decay * (8.0 / npi.powi(2)) * supply.pressure_drop_pa
                    / (gas.density * spec.length_m);
        }
        let exact = modes.iter().sum::<f64>() * spec.width_m * spec.height_m;
        squared_error += (supply.channel.flow_m3_s - exact).powi(2);
        squared_flow += exact.powi(2);
    }
    (squared_error / squared_flow).sqrt()
}

#[test]
fn coupled_channel_converges_to_analytic_unsteady_poiseuille_flow() {
    let coarse = analytic_flow_error(96000, 16);
    let refined = analytic_flow_error(192000, 32);
    assert!(refined < 0.8 * coarse, "coarse={coarse}, refined={refined}");
    assert!(refined < 0.02, "refined relative RMS flow error={refined}");
}

#[test]
fn rejected_sample_or_incomplete_gas_transfer_preserves_every_accepted_owner() {
    let mut model = system(slit(&gas(101325.0), 1.0 / 96000.0, 16).into(), true);
    let mut replay = model.clone();
    let drive = TubeDrive {
        upstream_pressure_pa: 0.5,
        body_flow_m3_s: 2e-10,
    };
    for _ in 0..47 {
        assert_eq!(model.step(drive).unwrap(), replay.step(drive).unwrap());
    }
    let down = model.load().state().clone();
    let up = model.upstream_load().unwrap().state().clone();
    let fluid = model.supply_channel().unwrap().velocities_m_s().to_vec();
    let exchange = model.supply_channel().unwrap().last_step();
    let mechanics = single(&model).state();
    for (bad, forces) in [
        (
            TubeDrive {
                upstream_pressure_pa: f64::NAN,
                ..drive
            },
            [0.0; 2],
        ),
        (
            TubeDrive {
                body_flow_m3_s: f64::INFINITY,
                ..drive
            },
            [0.0; 2],
        ),
        (drive, [0.0, 1.0]),
    ] {
        assert!(model.step_with_forces(bad, forces).is_err());
        assert_eq!(model.load().state(), &down);
        assert_eq!(model.upstream_load().unwrap().state(), &up);
        assert_eq!(model.supply_channel().unwrap().velocities_m_s(), fluid);
        assert_eq!(model.supply_channel().unwrap().last_step(), exchange);
        assert_eq!(single(&model).state(), mechanics);
    }
    let error = model
        .rebind_regional_gas_fields_controlled(None, &[], 1.0, &mut || false)
        .unwrap_err();
    assert!(error.to_string().contains("supply passage"));
    assert_eq!(model.step(drive).unwrap(), replay.step(drive).unwrap());
}

#[test]
fn binding_requires_complete_gas_clock_affine_channel_and_combined_impedance() {
    let material = gas(101325.0);
    let dt = 1.0 / 96000.0;
    let down = load(&material, dt, 0.05);
    let supply: FixedChannel = slit(&material, dt, 16).into();
    let (_, impedance) = inlet_state_with_supply(None, &down, &supply).unwrap();
    assert!(
        LoadedApertureDg::new_loaded_with_supply(
            aperture(&down, down.impedance_pa_s_m3()),
            None,
            down.clone(),
            &supply,
        )
        .is_err()
    );
    let mut nonlinear = *supply.slit().unwrap().spec();
    nonlinear.exit = SlitExit::UniformKineticHead;
    for other in [
        slit(&gas(81060.0), dt, 16).into(),
        slit(&material, dt / 2.0, 16).into(),
        SlitChannel::from_gas(&material, nonlinear).unwrap().into(),
    ] {
        assert!(inlet_state_with_supply(None, &down, &other).is_err());
    }
    let mut advanced = down.clone();
    advanced.preview_step(1e-10).unwrap().accept();
    assert!(
        LoadedApertureDg::new_loaded_with_supply(
            aperture(&down, impedance),
            None,
            advanced,
            &supply,
        )
        .is_err()
    );
    let mut moving = supply;
    moving.step(1.0).unwrap();
    let initial_energy = moving.stored_energy_j();
    let mut model = system(moving, false);
    assert_eq!(model.stored_energy_j(), initial_energy);
    let frame = model.step(TubeDrive::default()).unwrap();
    assert!(frame.stored_energy_j < initial_energy);
    assert!(
        frame.aperture.bore_flow_m3_s > 0.0,
        "original inertial history survives binding"
    );
}
