//! Focused constitutive/runtime regressions for the uniform viscous slit.

#![allow(clippy::float_cmp)] // Exact state preservation and equilibrium pins.

use fs_couple::bernoulli_aperture::{
    BernoulliAperture,
    dynamic::{ApertureDrive, ApertureState, DynamicAperture, DynamicApertureSpec},
    paired::BladePropagator,
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
