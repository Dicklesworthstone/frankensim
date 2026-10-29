//! Continuing source gas, nonlinear flow and retained terminal memory share one update.
//! Synthetic constitutive packs exercise implementation, not measured instruments.
use fs_couple::bernoulli_aperture::BernoulliAperture;
use fs_couple::bernoulli_aperture::dg::{ApertureDg, LoadedApertureDg, series_inlet_state};
use fs_couple::bernoulli_aperture::dynamic::{ApertureState, DynamicAperture, DynamicApertureSpec};
use fs_couple::bernoulli_aperture::paired::{
    ApertureBlade, BladePropagator, BladeState, PairedAperture, PairedApertureLimits,
    PairedApertureSpec,
};
use fs_couple::bernoulli_aperture::tube::TubeDrive;
use fs_couple::bernoulli_aperture::two_sided::LoadedAperture;
use fs_couple::bernoulli_aperture::viscous::ViscousSlit;
use fs_couple::dg_port::{DgGasRefresh, DgPortLoad, DgTerminalChannelInput, DgTerminalChannelKind};
use fs_dcontact::Obstacle;
use fs_duct::dg::{DgBoundary, DgDuct, DgMesh, DgState, MAX_DG_PROPAGATOR_DOFS, PortId};
use fs_duct::{Duct, Segment};
use fs_matdb::QueryPoint;
use fs_material::gas::GasState;
use fs_material::gas_field::{AxialGasField, PrimitiveGasState};
use fs_material::gas_species::input::ThermochemicalPhase;
use fs_material::thermochemical_gas::packs::from_packs;
use fs_math::c64::C64;
use fs_thermochem::{Composition, CompositionBasis};
use fs_vfit::mimo::MimoModel;
use fs_vfit::mimo::load::{MimoLoadLimits, PreparedMimoLoad};
use fs_vfit::model::PoleTerm;

#[path = "../src/render/schedule/reed/field_input/tests/support.rs"]
mod gas_support;

const DT: f64 = 1.0 / 96_000.0;
const RADIUS: f64 = 0.001;
const LENGTH: f64 = 0.05;
const WIDTH: f64 = 0.001;
const SLIT_LENGTH: f64 = 0.001;

fn composition(fraction: f64) -> Composition {
    Composition::new(
        CompositionBasis::MoleFraction,
        vec![fraction, 1.0 - fraction],
    )
    .unwrap()
}

fn field(temperature: f64, fraction: f64, pressure: f64, length: f64) -> AxialGasField {
    let a = gas_support::packs("A", 0.028, 3.5, 200.0);
    let b = gas_support::packs("B", 0.044, 4.5, 200.0);
    let model = from_packs(&[a.binding(), b.binding()], composition(0.5), None).unwrap();
    let primitive = PrimitiveGasState {
        temperature_k: temperature,
        composition: composition(fraction),
    };
    AxialGasField::new(
        vec![0.0, length],
        pressure,
        vec![primitive; 2],
        model,
        QueryPoint::new(),
        ThermochemicalPhase::FrozenGas,
    )
    .unwrap()
}

fn refresh(diffusion: bool) -> DgGasRefresh {
    if diffusion {
        DgGasRefresh::UniformDiffusion {
            maximum_relative_error: 0.01,
        }
    } else {
        DgGasRefresh::Lossless
    }
}

fn load(field: &AxialGasField, diffusion: bool) -> DgPortLoad {
    let duct = DgDuct::new_with_gas_field(
        &Duct {
            segments: vec![Segment::Cylinder {
                radius: RADIUS,
                length: field.length_m(),
            }],
        },
        field,
        DgMesh {
            elements: 1,
            order: 2,
        },
        DgBoundary::Port(PortId(0)),
        DgBoundary::Port(PortId(1)),
    )
    .unwrap();
    let duct = if diffusion {
        duct.with_viscothermal(8, [50.0, 2000.0], 0.01).unwrap()
    } else {
        duct
    };
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
    DgPortLoad::new_with_terminal(duct, DT, 1e-11, MAX_DG_PROPAGATOR_DOFS, terminal).unwrap()
}

fn aperture(load: &DgPortLoad, impedance: f64) -> DynamicAperture {
    let contact = Obstacle::new(
        vec![-1.0],
        1,
        1,
        vec![0.0],
        vec![1.0],
        0.0,
        2.0,
        "gas intervention regression: explicit zero contact".into(),
    )
    .unwrap();
    let mut aperture = DynamicAperture::new_with_propagator(
        DynamicApertureSpec {
            aperture: BernoulliAperture {
                rest_opening_m: 0.0002,
                width_m: WIDTH,
                closing_pressure_pa: 10_000.0,
            },
            mass_kg: 1e-5,
            stiffness_n_m: 500.0,
            damping_ratio: 0.0,
            density_kg_m3: load.duct().inlet_gas().density,
            impedance_pa_s_m3: impedance,
            time_step_s: DT,
            max_steps: 256,
        },
        ApertureState {
            opening_m: 0.0002,
            opening_velocity_m_s: 0.0,
        },
        contact,
        BladePropagator::ExactZoh,
    )
    .unwrap();
    aperture
        .set_viscous_slit(ViscousSlit::from_gas(load.duct().inlet_gas(), SLIT_LENGTH).unwrap())
        .unwrap();
    aperture
}

fn model(field: &AxialGasField, diffusion: bool) -> ApertureDg {
    let load = load(field, diffusion);
    ApertureDg::new(aperture(&load, load.impedance_pa_s_m3()), load).unwrap()
}

fn paired_source(load: &DgPortLoad, impedance: f64) -> PairedAperture {
    let contact = Obstacle::new(
        vec![-1.0],
        1,
        1,
        vec![0.0],
        vec![1.0],
        0.0,
        2.0,
        "unequal gas-driven blades: explicit zero contact".into(),
    )
    .unwrap();
    let mut pair = PairedAperture::new_with_propagator(
        PairedApertureSpec {
            rest_opening_m: 0.0002,
            width_m: WIDTH,
            blades: [
                ApertureBlade {
                    mass_kg: 1e-5,
                    stiffness_n_m: 500.0,
                    damping_n_s_m: 0.0,
                    pressure_area_m2: 1e-5,
                },
                ApertureBlade {
                    mass_kg: 2e-5,
                    stiffness_n_m: 700.0,
                    damping_n_s_m: 0.0,
                    pressure_area_m2: 1.3e-5,
                },
            ],
            density_kg_m3: load.duct().inlet_gas().density,
            impedance_pa_s_m3: impedance,
            time_step_s: DT,
            max_steps: 256,
            limits: PairedApertureLimits {
                maximum_abs_displacement_m: 0.001,
                maximum_abs_velocity_m_s: 10.0,
                maximum_energy_j: 1.0,
                max_evaluations: 512,
                velocity_tolerance_m_s: 1e-15,
                relative_tolerance: 1e-13,
            },
        },
        [
            BladeState {
                displacement_m: -1e-6,
                velocity_m_s: 0.005,
            },
            BladeState {
                displacement_m: 2e-6,
                velocity_m_s: -0.004,
            },
        ],
        contact,
        BladePropagator::ExactZoh,
    )
    .unwrap();
    pair.set_viscous_slit(ViscousSlit::from_gas(load.duct().inlet_gas(), SLIT_LENGTH).unwrap())
        .unwrap();
    pair
}

fn pair(model: &LoadedApertureDg) -> &PairedAperture {
    let LoadedAperture::Paired(source) = model.aperture() else {
        unreachable!()
    };
    source
}

fn drive(sample: usize) -> TubeDrive {
    TubeDrive {
        upstream_pressure_pa: if sample < 40 { 30.0 } else { -10.0 },
        body_flow_m3_s: if sample < 8 { 2e-9 } else { 0.0 },
    }
}

fn close(actual: f64, expected: f64, scale: f64) {
    assert!(
        (actual - expected).abs() <= 2e-9 * scale.max(1e-28),
        "actual={actual:e}, expected={expected:e}, scale={scale:e}"
    );
}

fn bulk(gas: &GasState) -> f64 {
    gas.density * gas.sound_speed.powi(2)
}

// Independent exact three-point Gauss integration of the quadratic pressure
// and volume-flow interpolants. No production mass matrix/energy map is used.
fn wave_energy(state: &DgState, gas: &GasState, length: f64) -> f64 {
    assert_eq!(state.pressure().len(), 3);
    let node = (3.0_f64 / 5.0).sqrt();
    let area = core::f64::consts::PI * RADIUS.powi(2);
    let interpolate = |values: &[f64], x: f64| {
        0.5 * x * (x - 1.0) * values[0]
            + (1.0 - x * x) * values[1]
            + 0.5 * x * (x + 1.0) * values[2]
    };
    [(-node, 5.0 / 9.0), (0.0, 8.0 / 9.0), (node, 5.0 / 9.0)]
        .iter()
        .map(|&(x, weight)| {
            let p = interpolate(state.pressure(), x);
            let q = interpolate(state.volume_flow(), x);
            length * weight * (area / bulk(gas) * p * p + gas.density / area * q * q) / 4.0
        })
        .sum()
}

#[test]
fn g1_source_gas_changes_live_jet_pressure_and_impedance_without_resetting_histories() {
    let cold = field(300.0, 1.0, 101_325.0, LENGTH);
    let hot = field(420.0, 0.2, 150_000.0, LENGTH);
    let old_gas = *cold.uniform_state().unwrap().state();
    let new_gas = *hot.uniform_state().unwrap().state();
    assert_ne!(old_gas.dynamic_viscosity, new_gas.dynamic_viscosity);
    assert_ne!(old_gas.thermal_conductivity, new_gas.thermal_conductivity);
    for diffusion in [false, true] {
        let mut changed = model(&cold, diffusion);
        let mut control = model(&cold, diffusion);
        for sample in 0..24 {
            assert_eq!(
                changed.step(drive(sample)).unwrap(),
                control.step(drive(sample)).unwrap()
            );
        }
        let before = changed.duct_state().clone();
        let opening = changed.aperture().state();
        let mechanics_energy = changed.aperture().stored_energy_j();
        let poles = changed.load().terminal().unwrap().pole_states().to_vec();
        let terminal_frames = changed.load().last_terminal_frames().to_vec();
        assert!(opening.opening_velocity_m_s != 0.0);
        assert!(poles.iter().any(|x| x.norm_sq() > 0.0));
        let old_wave_energy = wave_energy(&before, &old_gas, LENGTH);
        if diffusion {
            assert!(before.stored_energy_j() - old_wave_energy > 1e-12 * before.stored_energy_j());
        } else {
            close(before.stored_energy_j(), old_wave_energy, old_wave_energy);
        }
        let update = changed
            .rebind_gas_field_controlled(&hot, refresh(diffusion), 1.0, &mut || false)
            .unwrap();
        assert_eq!(changed.load().duct().gas_field(), Some(&hot));
        assert_eq!(
            changed
                .load()
                .duct()
                .gas_field()
                .unwrap()
                .model()
                .source_packs(),
            cold.model().source_packs()
        );
        assert_eq!(changed.aperture().spec().density_kg_m3, new_gas.density);
        assert_eq!(
            changed.aperture().spec().impedance_pa_s_m3,
            changed.load().impedance_pa_s_m3()
        );
        assert_eq!(
            changed.aperture().viscous_slit().unwrap().gas(),
            Some(&new_gas)
        );
        assert_eq!(
            changed.aperture().viscous_slit().unwrap().length_m(),
            SLIT_LENGTH
        );
        assert_eq!(changed.aperture().state(), opening);
        assert_eq!(changed.aperture().stored_energy_j(), mechanics_energy);
        assert_eq!(changed.aperture().accepted_steps(), 24);
        assert_eq!(changed.duct_state().step_count(), 24);
        assert_eq!(
            changed.duct_state().initial_energy_j(),
            before.initial_energy_j()
        );
        assert_eq!(
            changed.duct_state().cumulative_dissipated_j(),
            before.cumulative_dissipated_j()
        );
        assert_eq!(changed.load().terminal().unwrap().pole_states(), poles);
        assert_eq!(changed.load().last_terminal_frames(), terminal_frames);
        assert_ne!(update.parameter_work_j, 0.0);
        close(
            update.parameter_work_j,
            changed.duct_state().stored_energy_j() - before.stored_energy_j(),
            before.stored_energy_j() + changed.duct_state().stored_energy_j(),
        );
        let compression_scale = before
            .pressure()
            .iter()
            .map(|p| p.abs() / bulk(&old_gas))
            .fold(0.0_f64, f64::max);
        let momentum_scale = before
            .volume_flow()
            .iter()
            .map(|q| old_gas.density * q.abs())
            .fold(0.0_f64, f64::max);
        for i in 0..3 {
            close(
                changed.duct_state().pressure()[i] / bulk(&new_gas),
                before.pressure()[i] / bulk(&old_gas),
                compression_scale,
            );
            close(
                new_gas.density * changed.duct_state().volume_flow()[i],
                old_gas.density * before.volume_flow()[i],
                momentum_scale,
            );
        }
        if !diffusion {
            close(
                update.parameter_work_j,
                wave_energy(changed.duct_state(), &new_gas, LENGTH) - old_wave_energy,
                old_wave_energy + update.energy_after_j,
            );
        }
        let mut difference = 0.0;
        let mut signal = 0.0;
        for sample in 24..64 {
            let frame = changed.step(drive(sample)).unwrap();
            let reference = control.step(drive(sample)).unwrap();
            let aperture = frame.aperture;
            let gap = aperture.midpoint_opening_m;
            assert!(gap > 0.0);
            let velocity = aperture.jet_flow_m3_s / (WIDTH * gap);
            let constitutive_drop = 12.0 * new_gas.dynamic_viscosity * SLIT_LENGTH * velocity
                / gap.powi(2)
                + 0.5 * new_gas.density * velocity * velocity.abs();
            let drop = drive(sample).upstream_pressure_pa - aperture.bore_pressure_pa;
            close(
                drop,
                constitutive_drop,
                drive(sample).upstream_pressure_pa.abs() + aperture.bore_pressure_pa.abs(),
            );
            // This fixture has no structural damping, relaxation or contact.
            close(
                aperture.dissipated_energy_j,
                drop * aperture.jet_flow_m3_s * DT,
                aperture.dissipated_energy_j,
            );
            close(
                frame.duct.port_outputs[0].conjugate_rate,
                aperture.bore_pressure_pa,
                aperture.bore_pressure_pa.abs(),
            );
            assert_eq!(
                frame.duct.ledger.cumulative_parameter_work_j,
                update.parameter_work_j
            );
            close(
                frame.balance_residual_j(),
                0.0,
                frame.stored_energy_j
                    + frame.upstream_work_j.abs()
                    + frame.inferred_total_dissipation_j,
            );
            difference += (aperture.bore_pressure_pa - reference.aperture.bore_pressure_pa).powi(2);
            signal += reference.aperture.bore_pressure_pa.powi(2);
        }
        assert!(difference > 1e-6 * signal);
        assert_eq!(changed.duct_state().step_count(), 64);
    }
}

#[test]
fn g4_late_refresh_refusal_preserves_sources_memory_and_exact_retry() {
    let cold = field(300.0, 1.0, 101_325.0, LENGTH);
    let hot = field(420.0, 0.2, 150_000.0, LENGTH);
    let mut attempted = model(&cold, true);
    let mut clean = model(&cold, true);
    for sample in 0..24 {
        assert_eq!(
            attempted.step(drive(sample)).unwrap(),
            clean.step(drive(sample)).unwrap()
        );
    }
    let before = attempted.duct_state().clone();
    let opening = attempted.aperture().state();
    let spec = *attempted.aperture().spec();
    let slit = attempted.aperture().viscous_slit();
    let poles = attempted.load().terminal().unwrap().pole_states().to_vec();
    let terminal_frames = attempted.load().last_terminal_frames().to_vec();
    let response = attempted.load().port_response().unwrap();
    let mut polls = 0;
    let accepted = clean
        .rebind_gas_field_controlled(&hot, refresh(true), 1.0, &mut || {
            polls += 1;
            false
        })
        .unwrap();
    assert!(polls > 1 && accepted.parameter_work_j != 0.0);
    let assert_unchanged = |model: &ApertureDg| {
        assert_eq!(model.duct_state(), &before);
        assert_eq!(model.load().duct().gas_field(), Some(&cold));
        assert_eq!(model.aperture().state(), opening);
        assert_eq!(model.aperture().spec(), &spec);
        assert_eq!(model.aperture().viscous_slit(), slit);
        assert_eq!(model.aperture().accepted_steps(), 24);
        assert_eq!(model.load().terminal().unwrap().pole_states(), poles);
        assert_eq!(model.load().last_terminal_frames(), terminal_frames);
        assert_eq!(model.load().port_response().unwrap(), response);
    };
    assert!(
        attempted
            .rebind_gas_field_controlled(&hot, refresh(true), 0.0, &mut || false)
            .is_err()
    );
    assert_unchanged(&attempted);
    let mut calls = 0;
    assert!(
        attempted
            .rebind_gas_field_controlled(&hot, refresh(true), 1.0, &mut || {
                calls += 1;
                calls == polls
            })
            .is_err()
    );
    assert_eq!(calls, polls);
    assert_unchanged(&attempted);
    assert_eq!(
        attempted
            .rebind_gas_field_controlled(&hot, refresh(true), 1.0, &mut || false)
            .unwrap(),
        accepted
    );
    for sample in 24..48 {
        assert_eq!(
            attempted.step(drive(sample)).unwrap(),
            clean.step(drive(sample)).unwrap()
        );
    }
    assert_eq!(attempted.duct_state(), clean.duct_state());
    assert_eq!(
        attempted.duct_state().cumulative_parameter_work_j(),
        accepted.parameter_work_j
    );
}

#[test]
fn g3_series_refresh_keeps_two_loaded_tracts_and_source_on_one_gas_epoch() {
    let old = [
        field(300.0, 1.0, 101_325.0, 0.035),
        field(300.0, 1.0, 101_325.0, LENGTH),
    ];
    let new = [
        field(420.0, 0.2, 150_000.0, 0.035),
        field(420.0, 0.2, 150_000.0, LENGTH),
    ];
    let build = || {
        let up = load(&old[0], true);
        let down = load(&old[1], true);
        let (_, impedance) = series_inlet_state(&up, &down).unwrap();
        LoadedApertureDg::new_loaded_with_upstream(paired_source(&down, impedance), up, down)
            .unwrap()
    };
    let mut changed = build();
    let mut clean = build();
    for sample in 0..24 {
        assert_eq!(
            changed.step(drive(sample)).unwrap(),
            clean.step(drive(sample)).unwrap()
        );
    }
    let old_states = [
        changed.upstream_load().unwrap().state().clone(),
        changed.duct_state().clone(),
    ];
    let poles = [
        changed
            .upstream_load()
            .unwrap()
            .terminal()
            .unwrap()
            .pole_states()
            .to_vec(),
        changed.load().terminal().unwrap().pole_states().to_vec(),
    ];
    let motion = *pair(&changed).blade_states();
    let blades = pair(&changed).spec().blades;
    let mechanical_energy = pair(&changed).stored_energy_j();
    assert_ne!(motion[0], motion[1]);
    assert!(
        changed
            .rebind_gas_field_controlled(&new[1], refresh(true), 1.0, &mut || false)
            .is_err()
    );
    // The upstream proposal succeeds first; the incompatible downstream gas
    // must still prevent publication of either trajectory or source law.
    assert!(
        changed
            .rebind_series_gas_fields_controlled(
                [&new[0], &old[1]],
                [refresh(true); 2],
                1.0,
                &mut || false,
            )
            .is_err()
    );
    assert_eq!(changed.upstream_load().unwrap().state(), &old_states[0]);
    assert_eq!(changed.duct_state(), &old_states[1]);
    assert_eq!(
        changed
            .upstream_load()
            .unwrap()
            .terminal()
            .unwrap()
            .pole_states(),
        poles[0]
    );
    assert_eq!(changed.load().terminal().unwrap().pole_states(), poles[1]);
    assert_eq!(pair(&changed).blade_states(), &motion);
    assert_eq!(
        pair(&changed).gas(),
        old[1].uniform_state().map(|state| state.state())
    );
    let report = changed
        .rebind_series_gas_fields_controlled(
            [&new[0], &new[1]],
            [refresh(true); 2],
            1.0,
            &mut || false,
        )
        .unwrap();
    assert_eq!(
        clean
            .rebind_series_gas_fields_controlled(
                [&new[0], &new[1]],
                [refresh(true); 2],
                1.0,
                &mut || false,
            )
            .unwrap(),
        report
    );
    let source = pair(&changed);
    let gas = *changed.load().duct().inlet_gas();
    let impedance =
        changed.upstream_load().unwrap().impedance_pa_s_m3() + changed.load().impedance_pa_s_m3();
    assert_eq!(changed.upstream_load().unwrap().duct().inlet_gas(), &gas);
    assert_eq!(source.blade_states(), &motion);
    assert_eq!(source.spec().blades, blades);
    assert_eq!(source.stored_energy_j(), mechanical_energy);
    assert_eq!(source.gas(), Some(&gas));
    assert_eq!(source.spec().density_kg_m3, gas.density);
    assert_eq!(source.spec().impedance_pa_s_m3, impedance);
    assert_eq!(source.viscous_slit().unwrap().gas(), Some(&gas));
    for sample in 24..48 {
        let previous = *pair(&changed).blade_states();
        let frame = changed.step(drive(sample)).unwrap();
        assert_eq!(frame, clean.step(drive(sample)).unwrap());
        let pressure_drop = frame.upstream.as_ref().unwrap().pressure_drop_pa;
        for i in 0..2 {
            // Independent closed-form forced oscillator for each actual blade;
            // neither the combined gap nor a second source solve is an oracle.
            let b = blades[i];
            let omega = (b.stiffness_n_m / b.mass_kg).sqrt();
            let equilibrium = -b.pressure_area_m2 * pressure_drop / b.stiffness_n_m;
            let (sin, cos) = (omega * DT).sin_cos();
            let displacement = equilibrium
                + (previous[i].displacement_m - equilibrium) * cos
                + previous[i].velocity_m_s / omega * sin;
            let velocity = -(previous[i].displacement_m - equilibrium) * omega * sin
                + previous[i].velocity_m_s * cos;
            let amplitude = previous[i].displacement_m.abs()
                + equilibrium.abs()
                + previous[i].velocity_m_s.abs() / omega;
            close(
                pair(&changed).blade_states()[i].displacement_m,
                displacement,
                amplitude,
            );
            close(
                pair(&changed).blade_states()[i].velocity_m_s,
                velocity,
                omega * amplitude,
            );
        }
        assert_eq!(
            frame.duct.ledger.cumulative_parameter_work_j,
            report[1].parameter_work_j
        );
        assert_eq!(
            frame
                .upstream
                .as_ref()
                .unwrap()
                .duct
                .ledger
                .cumulative_parameter_work_j,
            report[0].parameter_work_j
        );
        close(
            frame.balance_residual_j(),
            0.0,
            frame.stored_energy_j
                + frame.upstream_work_j.abs()
                + frame.inferred_total_dissipation_j,
        );
    }
}

#[test]
fn g1_g4_live_gas_reaches_moving_aperture_and_finite_terminal_momentum_together() {
    let cold = field(300.0, 1.0, 101_325.0, LENGTH);
    let hot = field(420.0, 0.2, 150_000.0, LENGTH);
    let old_gas = *cold.uniform_state().unwrap().state();
    let new_gas = *hot.uniform_state().unwrap().state();
    for kind in [
        DgTerminalChannelKind::Slit {
            width_m: 0.002,
            height_m: 4e-5,
            length_m: 0.001,
            transverse_cells: 8,
            maximum_mach: 0.2,
        },
        DgTerminalChannelKind::Circular {
            radius_m: 4e-5,
            length_m: 0.001,
            radial_cells: 8,
            maximum_mach: 0.2,
        },
    ] {
        let build = || {
            let load = load(&cold, false)
                .with_terminal_channels(
                    vec![DgTerminalChannelInput {
                        port: PortId(1),
                        kind,
                    }],
                    8,
                )
                .unwrap();
            ApertureDg::new(aperture(&load, load.impedance_pa_s_m3()), load).unwrap()
        };
        let mut changed = build();
        let mut clean = build();
        for sample in 0..24 {
            assert_eq!(
                changed.step(drive(sample)).unwrap(),
                clean.step(drive(sample)).unwrap()
            );
        }
        let old_state = changed.duct_state().clone();
        let old_motion = changed.aperture().state();
        let channel = &changed.load().terminal_channels()[0];
        let old_velocity = channel.channel().unwrap().velocities_m_s().to_vec();
        let old_channel_step = channel.channel().unwrap().last_step();
        let old_channel_energy = channel.stored_energy_j();
        let old_energy = changed.load().stored_energy_j();
        let old_poles = changed.load().terminal().unwrap().pole_states().to_vec();
        assert!(old_channel_energy > 0.0);
        let mut polls = 0;
        let report = clean
            .rebind_gas_field_controlled(&hot, DgGasRefresh::Lossless, 1.0, &mut || {
                polls += 1;
                false
            })
            .unwrap();
        let mut calls = 0;
        assert!(
            changed
                .rebind_gas_field_controlled(&hot, DgGasRefresh::Lossless, 1.0, &mut || {
                    calls += 1;
                    calls == polls
                },)
                .is_err()
        );
        assert_eq!(calls, polls);
        assert_eq!(changed.duct_state(), &old_state);
        assert_eq!(changed.aperture().state(), old_motion);
        assert_eq!(changed.load().terminal_channels()[0].gas(), &old_gas);
        assert_eq!(
            changed.load().terminal_channels()[0]
                .channel()
                .unwrap()
                .velocities_m_s(),
            old_velocity
        );
        assert_eq!(
            changed
                .rebind_gas_field_controlled(&hot, DgGasRefresh::Lossless, 1.0, &mut || false,)
                .unwrap(),
            report
        );
        let channel = &changed.load().terminal_channels()[0];
        assert_eq!(channel.gas(), &new_gas);
        assert_eq!(channel.channel().unwrap().gas(), &new_gas);
        assert_eq!(channel.channel().unwrap().last_step(), old_channel_step);
        assert_eq!(changed.aperture().state(), old_motion);
        assert_eq!(changed.aperture().spec().density_kg_m3, new_gas.density);
        assert_eq!(changed.load().terminal().unwrap().pole_states(), old_poles);
        assert_eq!(report.terminal_channels.len(), 1);
        let momentum_scale = old_velocity
            .iter()
            .map(|v| (old_gas.density * v).abs())
            .fold(0.0_f64, f64::max);
        for (&before, &after) in old_velocity
            .iter()
            .zip(channel.channel().unwrap().velocities_m_s())
        {
            close(
                new_gas.density * after,
                old_gas.density * before,
                momentum_scale,
            );
        }
        let expected_channel_energy = old_channel_energy * old_gas.density / new_gas.density;
        close(
            channel.stored_energy_j(),
            expected_channel_energy,
            expected_channel_energy,
        );
        close(
            report.parameter_work_j,
            changed.load().stored_energy_j() - old_energy,
            old_energy + changed.load().stored_energy_j(),
        );
        close(
            report.parameter_work_j - report.duct.parameter_work_j,
            expected_channel_energy - old_channel_energy,
            old_channel_energy + expected_channel_energy,
        );
        for sample in 24..48 {
            let frame = changed.step(drive(sample)).unwrap();
            assert_eq!(frame, clean.step(drive(sample)).unwrap());
            close(
                frame.duct.ledger.cumulative_parameter_work_j,
                report.parameter_work_j,
                report.absolute_parameter_work_j,
            );
            close(
                frame.balance_residual_j(),
                0.0,
                frame.stored_energy_j
                    + frame.upstream_work_j.abs()
                    + frame.inferred_total_dissipation_j,
            );
            close(
                frame.duct.ledger.energy_balance_residual_j,
                0.0,
                frame.duct.ledger.stored_energy_j
                    + frame.duct.ledger.cumulative_port_work_j.abs()
                    + frame.duct.ledger.cumulative_dissipated_j
                    + report.absolute_parameter_work_j,
            );
        }
    }
}
