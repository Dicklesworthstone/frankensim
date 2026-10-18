//! Focused constitutive/runtime regressions for the uniform viscous slit.

#![allow(clippy::float_cmp)] // Exact state preservation and equilibrium pins.

use fs_couple::bernoulli_aperture::{
    BernoulliAperture,
    dynamic::{ApertureDrive, ApertureState, DynamicAperture, DynamicApertureSpec},
    paired::{ApertureBlade, BladePropagator, BladeState, PairedAperture, PairedApertureLimits, PairedApertureSpec},
    modal::{ApertureMode, ModalAperture, ModalApertureSpec},
    viscous::ViscousSlit,
};
use fs_couple::unilateral_contact::slit_lay;

fn model(propagator: BladePropagator) -> DynamicAperture {
    let spec = DynamicApertureSpec {
        aperture: BernoulliAperture {
            rest_opening_m: 4e-5,
            width_m: 0.013,
            closing_pressure_pa: 6000.0,
        },
        mass_kg: 1e-5,
        stiffness_n_m: 500.0,
        damping_ratio: 0.35,
        density_kg_m3: 1.2,
        impedance_pa_s_m3: 1e6,
        time_step_s: 1e-6,
        max_steps: 4096,
    };
    DynamicAperture::new_with_propagator(
        spec,
        ApertureState { opening_m: spec.aperture.rest_opening_m, opening_velocity_m_s: 0.0 },
        slit_lay(1e8, 1.5).unwrap(),
        propagator,
    ).unwrap()
}

fn paired(propagator: BladePropagator) -> PairedAperture {
    PairedAperture::new_with_propagator(PairedApertureSpec {
        rest_opening_m: 4e-5, width_m: 0.013,
        blades: [ApertureBlade { mass_kg: 1e-5, stiffness_n_m: 500.0,
            damping_n_s_m: 0.03, pressure_area_m2: 8e-6 },
            ApertureBlade { mass_kg: 1.5e-5, stiffness_n_m: 850.0,
                damping_n_s_m: 0.025, pressure_area_m2: 6e-6 }],
        density_kg_m3: 1.2, impedance_pa_s_m3: 1e6, time_step_s: 1e-6, max_steps: 4096,
        limits: PairedApertureLimits { maximum_abs_displacement_m: 0.01,
            maximum_abs_velocity_m_s: 100.0, maximum_energy_j: 1.0,
            max_evaluations: 256, velocity_tolerance_m_s: 1e-13, relative_tolerance: 1e-11 },
    }, [BladeState::default(); 2], slit_lay(1e8, 1.5).unwrap(), propagator).unwrap()
}

fn assert_account(frame: &fs_couple::bernoulli_aperture::dynamic::ApertureFrame) {
    let flow_scale = frame.jet_flow_m3_s.abs()+frame.swept_flow_m3_s.abs()+frame.bore_flow_m3_s.abs()+1e-15;
    assert!(frame.flow_residual_m3_s.abs() < 2e-9*flow_scale);
    let work_scale = frame.storage_change_j.abs()+frame.dissipated_energy_j.abs()
        +frame.pressure_work_j.abs()+frame.mechanical_work_j.abs()+1e-15;
    assert!(frame.balance_residual_j().abs() < 2e-8*work_scale);
    assert!(frame.dissipated_energy_j >= 0.0);
}

#[test]
fn g1_independent_blades_use_one_viscous_channel_in_both_propagators() {
    for propagator in [BladePropagator::Midpoint, BladePropagator::ExactZoh] {
        let mut inviscid = paired(propagator);
        let mut viscous = paired(propagator);
        let law = ViscousSlit::new(1.81e-5, 1e-3).unwrap();
        viscous.set_viscous_slit(law).unwrap();
        let drive = ApertureDrive { upstream_pressure_pa: 2000.0, ..ApertureDrive::default() };
        let old = inviscid.step(drive).unwrap();
        let first = viscous.step(drive).unwrap();
        assert!(first.aperture.jet_flow_m3_s < old.aperture.jet_flow_m3_s);
        assert_ne!(first.blades, old.blades);
        for sample in 0..512 {
            let drive = ApertureDrive { upstream_pressure_pa: if sample < 384 { 20000.0 } else { -2000.0 },
                ..ApertureDrive::default() };
            let frame = viscous.step_with_forces(drive, [0.001, -0.0003]).unwrap();
            let expected = law.evaluate(viscous.spec().width_m, frame.aperture.midpoint_opening_m,
                viscous.spec().density_kg_m3, frame.pressure_drop_pa).unwrap();
            assert_eq!(frame.aperture.jet_flow_m3_s, expected.flow_m3_s);
            assert_account(&frame.aperture);
        }
        let state = *viscous.blade_states();
        let steps = viscous.accepted_steps();
        assert!(viscous.set_viscous_slit(law).is_err());
        assert_eq!(*viscous.blade_states(), state);
        assert_eq!(viscous.accepted_steps(), steps);
    }
}

#[test]
fn g3_modal_and_paired_representations_preserve_shared_viscous_flow_and_work() {
    let mut pair = paired(BladePropagator::Midpoint);
    let p = *pair.spec();
    let modes = p.blades.iter().enumerate().map(|(i,b)| ApertureMode {
        mass_kg:b.mass_kg, stiffness_n_m:b.stiffness_n_m, damping_n_s_m:b.damping_n_s_m,
        pressure_area_m2:b.pressure_area_m2, gap_projection:1.0,
        force_projection: if i==0 {[1.0,0.0]} else {[0.0,1.0]},
    }).collect();
    let mut modal = ModalAperture::new(ModalApertureSpec {
        rest_opening_m:p.rest_opening_m, width_m:p.width_m, modes, density_kg_m3:p.density_kg_m3,
        impedance_pa_s_m3:p.impedance_pa_s_m3, time_step_s:p.time_step_s,
        max_steps:p.max_steps, limits:p.limits,
    }, vec![BladeState::default();2], slit_lay(1e8,1.5).unwrap()).unwrap();
    let law = ViscousSlit::new(1.81e-5,1e-3).unwrap();
    pair.set_viscous_slit(law).unwrap(); modal.set_viscous_slit(law).unwrap();
    for sample in 0..512 {
        let input = ApertureDrive { upstream_pressure_pa:if sample<384 {20000.0} else {-2000.0},
            ..ApertureDrive::default() };
        let a = pair.step_with_forces(input,[0.001,-0.0003]).unwrap();
        let b = modal.step(input,[0.001,-0.0003]).unwrap();
        for (a,b) in a.blades.iter().zip(&b.modes) {
            assert!((a.displacement_m-b.displacement_m).abs()<1e-11);
            assert!((a.velocity_m_s-b.velocity_m_s).abs()<1e-6);
        }
        assert!((a.aperture.jet_flow_m3_s-b.aperture.jet_flow_m3_s).abs()<1e-12);
        assert_account(&b.aperture);
    }
    let state = modal.modal_states().to_vec();
    assert!(modal.set_viscous_slit(law).is_err());
    assert_eq!(modal.modal_states(),state);
}

#[test]
fn g1_viscosity_changes_the_solved_jet_and_preserves_local_work_accounting() {
    for propagator in [BladePropagator::Midpoint, BladePropagator::ExactZoh] {
        for upstream in [-2000.0, 2000.0] {
            let mut inviscid = model(propagator);
            let mut viscous = model(propagator);
            let law = ViscousSlit::new(1.81e-5, 1e-3).unwrap();
            viscous.set_viscous_slit(law).unwrap();
            assert_eq!(viscous.viscous_slit(), Some(law));
            assert_eq!(inviscid.viscous_slit(), None);
            let drive = ApertureDrive {
                upstream_pressure_pa: upstream,
                incoming_pressure_pa: 0.0,
                body_flow_m3_s: 0.0,
            };
            let original = inviscid.step(drive).unwrap();
            let first = viscous.step(drive).unwrap();
            assert!(first.jet_flow_m3_s.abs() < original.jet_flow_m3_s.abs());
            assert_ne!(first.state.opening_m, original.state.opening_m);
            for _ in 0..128 {
                let frame = viscous.step(drive).unwrap();
                let dp = upstream - frame.bore_pressure_pa;
                let flow = law.evaluate(0.013, frame.midpoint_opening_m, 1.2, dp).unwrap();
                assert_eq!(frame.jet_flow_m3_s, flow.flow_m3_s);
                let flow_scale = frame.jet_flow_m3_s.abs() + frame.swept_flow_m3_s.abs()
                    + frame.bore_flow_m3_s.abs();
                assert!(frame.flow_residual_m3_s.abs() <= 2e-10 * flow_scale.max(1e-15));
                let work_scale = frame.storage_change_j.abs() + frame.dissipated_energy_j.abs()
                    + frame.pressure_work_j.abs();
                assert!(frame.balance_residual_j().abs() <= 2e-9 * work_scale.max(1e-15));
                assert!(frame.dissipated_energy_j >= 0.0);
            }
        }
    }
}

#[test]
fn g0_refused_input_and_late_law_change_preserve_the_accepted_state() {
    let mut aperture = model(BladePropagator::Midpoint);
    let law = ViscousSlit::new(1.81e-5, 1e-3).unwrap();
    aperture.set_viscous_slit(law).unwrap();
    let initial = aperture.state();
    let unrepresentable = ViscousSlit::new(1e150, 1e150).unwrap();
    assert!(aperture.set_viscous_slit(unrepresentable).is_err());
    assert_eq!(aperture.state(), initial);
    assert_eq!(aperture.accepted_steps(), 0);
    assert_eq!(aperture.viscous_slit(), Some(law));
    aperture.step(ApertureDrive { upstream_pressure_pa: 1000.0, ..ApertureDrive::default() }).unwrap();
    let state = aperture.state();
    let count = aperture.accepted_steps();
    let energy = aperture.stored_energy_j();
    let other = ViscousSlit::new(2e-5, 2e-3).unwrap();
    assert!(aperture.set_viscous_slit(other).is_err());
    assert!(aperture.step(ApertureDrive {
        upstream_pressure_pa: f64::NAN,
        ..ApertureDrive::default()
    }).is_err());
    assert_eq!(aperture.state(), state);
    assert_eq!(aperture.accepted_steps(), count);
    assert_eq!(aperture.stored_energy_j(), energy);
    assert_eq!(aperture.viscous_slit(), Some(law));
}

#[test]
fn g1_zero_drive_does_not_invent_jet_flow_or_mechanical_energy() {
    for propagator in [BladePropagator::Midpoint, BladePropagator::ExactZoh] {
        let mut aperture = model(propagator);
        aperture.set_viscous_slit(ViscousSlit::new(1.81e-5, 1e-3).unwrap()).unwrap();
        for _ in 0..16 {
            let frame = aperture.step(ApertureDrive::default()).unwrap();
            assert_eq!(frame.jet_flow_m3_s, 0.0);
            assert_eq!(frame.swept_flow_m3_s, 0.0);
            assert_eq!(frame.stored_energy_j, 0.0);
            assert_eq!(frame.dissipated_energy_j, 0.0);
        }
    }
}

#[test]
fn g1_forced_viscous_contact_cycle_preserves_flow_and_work_in_both_propagators() {
    for propagator in [BladePropagator::Midpoint, BladePropagator::ExactZoh] {
        let mut aperture = model(propagator);
        let law = ViscousSlit::new(1.81e-5, 1e-3).unwrap();
        aperture.set_viscous_slit(law).unwrap();
        let mut transitions = [false; 4];
        for sample in 0..512 {
            let was_open = aperture.state().opening_m >= 0.0;
            let drive = ApertureDrive {
                upstream_pressure_pa: if sample < 256 { 20000.0 } else { -2000.0 },
                incoming_pressure_pa: 75.0,
                body_flow_m3_s: 2e-7,
            };
            let force_n = if sample % 64 < 32 { 0.001 } else { -0.001 };
            let frame = aperture.step_with_force(drive, force_n).unwrap();
            let is_open = frame.state.opening_m >= 0.0;
            transitions[usize::from(was_open) * 2 + usize::from(is_open)] = true;
            let dp = drive.upstream_pressure_pa - frame.bore_pressure_pa;
            let flow = law.evaluate(0.013, frame.midpoint_opening_m, 1.2, dp).unwrap();
            assert_eq!(frame.jet_flow_m3_s, flow.flow_m3_s);
            let flow_scale = frame.jet_flow_m3_s.abs() + frame.swept_flow_m3_s.abs()
                + frame.bore_flow_m3_s.abs() + drive.body_flow_m3_s.abs();
            assert!(frame.flow_residual_m3_s.abs() <= 2e-10 * flow_scale.max(1e-15));
            let work_scale = frame.storage_change_j.abs() + frame.dissipated_energy_j.abs()
                + frame.pressure_work_j.abs() + frame.mechanical_work_j.abs();
            assert!(frame.balance_residual_j().abs() <= 2e-9 * work_scale.max(1e-15));
            assert!(frame.dissipated_energy_j >= 0.0);
        }
        assert!(transitions.iter().all(|seen| *seen), "{propagator:?}: {transitions:?}");
    }
}

#[test]
fn source_bound_channels_retain_full_gas_and_refuse_acoustically_equal_substitutions() {
    use fs_material::gas::{GasSpec, GasState};
    let source = GasSpec::dry_air_ussa1976();
    let gas = GasState::try_new(&source, 293.15, 101325.0).unwrap();
    // A synthetic coherent gas isolates the full-state binding: rho and c
    // match, while temperature, pressure and viscosity are different.
    let ratio = 1.1;
    let other = GasState::try_new(&GasSpec { gamma: source.gamma/ratio, ..source },
        293.15*ratio, 101325.0*ratio).unwrap();
    assert!((other.density-gas.density).abs() < 64.0*f64::EPSILON*gas.density);
    assert!((other.sound_speed-gas.sound_speed).abs() < 64.0*f64::EPSILON*gas.sound_speed);
    assert_ne!(other.dynamic_viscosity, gas.dynamic_viscosity);
    let mut p = *paired(BladePropagator::Midpoint).spec();
    p.density_kg_m3 = gas.density;
    let modes = p.blades.iter().enumerate().map(|(i,b)| ApertureMode {
        mass_kg:b.mass_kg, stiffness_n_m:b.stiffness_n_m, damping_n_s_m:b.damping_n_s_m,
        pressure_area_m2:b.pressure_area_m2, gap_projection:1.0,
        force_projection:if i==0 {[1.0,0.0]} else {[0.0,1.0]},
    }).collect();
    let mut pair = PairedAperture::new(p, [BladeState::default();2], slit_lay(1e8,1.5).unwrap()).unwrap();
    let mut modal = ModalAperture::new(ModalApertureSpec {
        rest_opening_m:p.rest_opening_m, width_m:p.width_m, modes,
        density_kg_m3:p.density_kg_m3, impedance_pa_s_m3:p.impedance_pa_s_m3,
        time_step_s:p.time_step_s, max_steps:p.max_steps, limits:p.limits,
    }, vec![BladeState::default();2], slit_lay(1e8,1.5).unwrap()).unwrap();
    let law = ViscousSlit::from_gas(&gas, 1e-3).unwrap();
    let mismatched = ViscousSlit::from_gas(&other, 1e-3).unwrap();
    // This would pass the old density-only setter admission.
    assert!(mismatched.evaluate(p.width_m, p.rest_opening_m, gas.density, 0.0).is_ok());
    pair.set_viscous_slit(law).unwrap();
    modal.set_viscous_slit(law).unwrap();
    assert_eq!(pair.gas(), Some(&gas));
    assert_eq!(modal.gas(), Some(&gas));
    assert!(pair.set_viscous_slit(mismatched).is_err());
    assert!(modal.set_viscous_slit(mismatched).is_err());
    assert_eq!(pair.gas(), Some(&gas));
    assert_eq!(modal.gas(), Some(&gas));
    assert_eq!(pair.viscous_slit(), Some(law));
    assert_eq!(modal.viscous_slit(), Some(law));
    assert_eq!(pair.accepted_steps(), 0);
    assert_eq!(modal.accepted_steps(), 0);
    assert_eq!(pair.blade_states(), &[BladeState::default();2]);
    assert_eq!(modal.modal_states(), &[BladeState::default();2]);
    pair.set_material_viscous_slit(2e-3).unwrap();
    modal.set_material_viscous_slit(2e-3).unwrap();
    assert_eq!(pair.viscous_slit().unwrap().gas(), Some(&gas));
    assert_eq!(modal.viscous_slit().unwrap().gas(), Some(&gas));
}
