//! Thermal loss on actual shared string/shell/air mechanics; the coarse FE
//! fixture exercises integration and does not resolve audio boundary layers.
use super::*;
use fs_couple::bowed_string::runtime::reciprocal::bank::cavity::ShellCavityThermalOptions;

#[path = "thermal/molecular.rs"]
mod molecular;

fn pressure_resolution()->ShellCavityOptions {
    let mut o=resolution();o.maximum_hz=1.0;o
}
fn heat_resolution()->ShellCavityThermalOptions {
    ShellCavityThermalOptions { maximum_decay_rate_per_s:100.0,
        maximum_thermal_modes:32,maximum_thermal_dofs:64,
        maximum_projection_terms:10_000_000,maximum_unresolved_fraction:0.5 }
}

#[test]
fn shared_bank_thermal_history_changes_pressure_and_survives_gas_rebuilds() {
    with_cx(|cx,gate|with_volume(2,gas(101325.0),|fe| {
        let atmosphere=gas(101325.0);
        let mut adiabatic=bank(cx,gate,fe,atmosphere,700.0)
            .with_material_shell_cavity(fe,pressure_resolution(),gate).unwrap();
        let mut thermal=bank(cx,gate,fe,atmosphere,700.0)
            .with_material_shell_cavity_thermal(fe,pressure_resolution(),Vec::new(),heat_resolution(),gate).unwrap();
        let image=thermal.source().spatial_cavity().unwrap();
        let memory=image.thermal_loss().unwrap().memory_states();
        assert!(memory>0);
        assert_eq!(thermal.system().material_memory().len(),memory);
        assert_eq!(thermal.system().mode_count(),adiabatic.system().mode_count());
        let mut scratch=vec![0.0;image.observation_state_len()];
        let mut a=[0.0];let mut b=[0.0];
        assert!(image.coupling().pressures_into(&scratch,&mut a).is_err());
        let mut force=vec![0.0;thermal.system().mode_count()];force[4]=0.02;
        let released=FrictionDrive { speed_m_s:0.0,normal_force_n:0.0 };
        let mut difference=0.0;
        for _ in 0..128 {
            adiabatic.step(&force,released).unwrap();
            let frame=*thermal.step(&force,released).unwrap();
            assert!(frame.energy_residual_j.abs()<=frame.energy_tolerance_j);
            adiabatic.spatial_cavity_pressures_into(&mut scratch,&mut a).unwrap();
            thermal.spatial_cavity_pressures_into(&mut scratch,&mut b).unwrap();
            difference+=(a[0]-b[0]).abs();
        }
        assert!(thermal.system().material_memory().iter().any(|z|*z!=0.0));
        assert!(difference>0.0);

        let flow=CylinderFlowResolution { minimum_wavenumber_radius:0.01,
            maximum_wavenumber_radius:4.0,intervals:2 };
        let mut first=bank(cx,gate,fe,atmosphere,700.0)
            .with_material_shell_cavity_thermal(fe,pressure_resolution(),Vec::new(),heat_resolution(),gate).unwrap()
            .with_cylinder_gas_loading(flow,gate).unwrap();
        let mut second=bank(cx,gate,fe,atmosphere,700.0)
            .with_cylinder_gas_loading(flow,gate).unwrap()
            .with_material_shell_cavity_thermal(fe,pressure_resolution(),Vec::new(),heat_resolution(),gate).unwrap();
        assert_eq!(first.system().material_memory().len(),memory);
        assert_eq!(first.source().body_connections(),second.source().body_connections());
        let zero=vec![0.0;first.system().mode_count()];
        for _ in 0..16 {
            assert_eq!(first.step(&zero,drive()).unwrap(),second.step(&zero,drive()).unwrap());
        }
        assert_eq!(first.system().material_memory(),second.system().material_memory());
        assert!(first.with_material_shell_cavity_thermal(fe,pressure_resolution(),Vec::new(),heat_resolution(),gate).is_err());
    }));
}
