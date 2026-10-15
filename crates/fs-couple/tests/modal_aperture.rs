use fs_couple::bernoulli_aperture::dynamic::ApertureDrive;
use fs_couple::bernoulli_aperture::modal::{ApertureMode, ModalAperture, ModalApertureSpec};
use fs_couple::bernoulli_aperture::paired::{ApertureBlade, BladeState, PairedAperture,
    PairedApertureLimits, PairedApertureSpec};
use fs_dcontact::Obstacle;

fn limits() -> PairedApertureLimits {
    PairedApertureLimits { maximum_abs_displacement_m: 0.01, maximum_abs_velocity_m_s: 100.,
        maximum_energy_j: 1., max_evaluations: 256,
        velocity_tolerance_m_s: 1e-13, relative_tolerance: 1e-11 }
}
fn contact() -> Obstacle {
    Obstacle::new(vec![-1.], 1, 1, vec![0.], vec![1.], 1e8, 2.,
        "explicit synthetic modal slit contact".into()).unwrap()
}
fn spec() -> ModalApertureSpec {
    ModalApertureSpec { rest_opening_m: 4e-4, width_m: 0.013,
        modes: vec![
            ApertureMode { mass_kg: 1e-5, stiffness_n_m: 500., damping_n_s_m: 0.03,
                pressure_area_m2: 1.25e-4, gap_projection: 1., force_projection: [1., 0.] },
            ApertureMode { mass_kg: 8e-6, stiffness_n_m: 2500., damping_n_s_m: 0.02,
                pressure_area_m2: -4e-5, gap_projection: -0.2, force_projection: [-0.3, 0.7] },
            ApertureMode { mass_kg: 6e-6, stiffness_n_m: 6000., damping_n_s_m: 0.01,
                pressure_area_m2: 7e-5, gap_projection: 0., force_projection: [0.2, -0.4] },
        ], density_kg_m3: 1.2, impedance_pa_s_m3: 1e6, time_step_s: 1./96000.,
        max_steps: 4096, limits: limits() }
}
fn model(s: ModalApertureSpec) -> ModalAperture {
    let state = vec![BladeState::default(); s.modes.len()];
    ModalAperture::new(s, state, contact()).unwrap()
}
fn drive(p: f64) -> ApertureDrive {
    ApertureDrive { upstream_pressure_pa: p, incoming_pressure_pa: 0., body_flow_m3_s: 0. }
}

#[test]
fn two_uniform_modes_reproduce_original_independent_blade_mechanics() {
    let mut s = spec(); s.modes.truncate(2);
    s.modes[1].pressure_area_m2 = 8e-5; s.modes[1].gap_projection = 1.;
    s.modes[1].force_projection = [0., 1.];
    let b = |m: ApertureMode| ApertureBlade { mass_kg: m.mass_kg,
        stiffness_n_m: m.stiffness_n_m, damping_n_s_m: m.damping_n_s_m,
        pressure_area_m2: m.pressure_area_m2 };
    let mut old = PairedAperture::new(PairedApertureSpec {
        rest_opening_m: s.rest_opening_m, width_m: s.width_m,
        blades: [b(s.modes[0]), b(s.modes[1])], density_kg_m3: s.density_kg_m3,
        impedance_pa_s_m3: s.impedance_pa_s_m3, time_step_s: s.time_step_s,
        max_steps: s.max_steps, limits: s.limits,
    }, [BladeState::default(); 2], contact()).unwrap();
    let mut new = model(s);
    for n in 0..1024 {
        let input = drive(if n < 256 { 50. } else { 0. });
        let forces = if n < 123 { [0.001, -0.0003] } else { [0.; 2] };
        let a = old.step_with_forces(input, forces).unwrap();
        let b = new.step(input, forces).unwrap();
        for (x,y) in a.blades.iter().zip(&b.modes) {
            assert!((x.displacement_m-y.displacement_m).abs() < 1e-12);
            assert!((x.velocity_m_s-y.velocity_m_s).abs() < 1e-8);
        }
        assert!((a.aperture.outgoing_pressure_pa-b.aperture.outgoing_pressure_pa).abs() < 1e-6);
    }
}

#[test]
fn signed_basis_changes_leave_physical_gap_flow_and_energy_unchanged() {
    let s = spec(); let mut transformed = s.clone(); let scales = [-2., 0.5, -1.5];
    for (m,scale) in transformed.modes.iter_mut().zip(scales) {
        m.mass_kg *= scale*scale; m.stiffness_n_m *= scale*scale;
        m.damping_n_s_m *= scale*scale; m.pressure_area_m2 *= scale;
        m.gap_projection *= scale;
        for f in &mut m.force_projection { *f *= scale; }
    }
    let (mut a, mut b) = (model(s), model(transformed));
    for n in 0..512 {
        let input = drive(if n < 128 { 80. } else { 0. });
        let force = [if n < 200 { 0.001 } else { 0. }, -0.0002];
        let x = a.step(input,force).unwrap(); let y = b.step(input,force).unwrap();
        for ((u,v),scale) in x.modes.iter().zip(&y.modes).zip(scales) {
            assert!((u.displacement_m-v.displacement_m*scale).abs() < 1e-12);
            assert!((u.velocity_m_s-v.velocity_m_s*scale).abs() < 1e-8);
        }
        assert!((x.aperture.outgoing_pressure_pa-y.aperture.outgoing_pressure_pa).abs() < 1e-6);
        assert!((x.aperture.stored_energy_j-y.aperture.stored_energy_j).abs() < 1e-12);
    }
}

#[test]
fn zero_gap_projection_mode_still_sweeps_fluid_and_survives_source_release() {
    let mut a = model(spec());
    let mut reduced = spec(); reduced.modes.pop();
    let mut b = model(reduced);
    let mut changed = false;
    for n in 0..1024 {
        let input = drive(if n < 256 { 50. } else { 0. });
        let x = a.step(input,[0.;2]).unwrap(); let y = b.step(input,[0.;2]).unwrap();
        changed |= (x.aperture.outgoing_pressure_pa-y.aperture.outgoing_pressure_pa).abs() > 1e-6;
        assert!(x.aperture.balance_residual_j().abs() < 1e-12);
    }
    assert!(changed);
    assert!(a.modal_states()[2].displacement_m.abs() > 1e-12);
    assert!(a.stored_energy_j() > 0.);
}

#[test]
fn failed_modal_steps_publish_no_state_and_missing_ports_are_not_silently_ignored() {
    let mut a = model(spec()); a.step(drive(50.), [0.;2]).unwrap();
    let before = a.modal_states().to_vec(); let energy = a.stored_energy_j();
    for bad in [f64::NAN, f64::INFINITY, f64::MAX] {
        assert!(a.step(drive(bad), [0.;2]).is_err());
        assert_eq!(a.modal_states(), before); assert_eq!(a.stored_energy_j(),energy);
        assert_eq!(a.accepted_steps(),1);
    }
    let mut no_port = spec(); for m in &mut no_port.modes { m.force_projection[1] = 0.; }
    assert!(model(no_port).step(drive(50.), [0.,1.]).is_err());
    let s = spec(); assert!(ModalAperture::new(s, vec![], contact()).is_err());
    let mut too_many = spec(); too_many.modes.resize(33,too_many.modes[0]);
    assert!(ModalAperture::new(too_many,vec![BladeState::default();33],contact()).is_err());
    let mut short = spec(); short.max_steps = 1; let mut short = model(short);
    short.step(drive(50.),[0.;2]).unwrap();
    assert!(short.step(drive(50.),[0.;2]).is_err());
    short.extend_step_budget(2).unwrap(); short.step(drive(0.),[0.;2]).unwrap();
    assert_eq!(short.accepted_steps(),2);
}

#[path = "modal_aperture/network.rs"]
mod network;
