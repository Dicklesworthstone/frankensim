//! Load/contact design derivatives of an actual nonlinear preload.
//! This example emits optimizer-ready local derivatives; it does not claim to
//! identify real material parameters from synthetic data or to optimize audio.
use fs_couple::modal_acoustic_time::{ModalAcousticMode,ModalAcousticTimeBudget,ModalAcousticTimeModel};
use fs_couple::render::schedule::force::coupled::{CoupledModalSystem,ModalAttachment,ModalConnection,ModalCouplingConfig};
use fs_couple::render::schedule::force::coupled::contact::{ModalContact,ModalContactConfig};
use fs_couple::render::schedule::force::coupled::contact::multiple::MultiContactConfig;
use fs_couple::render::schedule::force::coupled::equilibrium::sensitivity::{EquilibriumLinearization,SensitivityBudget};
use fs_couple::render::schedule::force::coupled::equilibrium::sensitivity::objective::DisplacementTarget;
use fs_dcontact::Obstacle;
use fs_exec::CancelGate;
use fs_math::c64::C64;

fn run()->Result<(),Box<dyn std::error::Error>> {
    let gate=CancelGate::new();
    // No acoustic transfer or artificial startup: the question is static design.
    let mass=ModalAcousticTimeModel::try_free_mass(48000,0.04,0.0,0.0,ModalAcousticTimeBudget::audible_reference())?;
    let receiver=ModalAcousticTimeModel::try_new(48000,vec![ModalAcousticMode {
        angular_frequency_rad_s:100.0,damping_ratio:0.02,pressure_per_modal_velocity:C64::ZERO,
    }],ModalAcousticTimeBudget::audible_reference())?;
    let actuator=ModalAttachment {component:0,shapes:vec![5.0]}; // 1/sqrt(0.04 kg)
    let receiver_port=ModalAttachment {component:1,shapes:vec![1.0]};
    let support=ModalConnection {left:actuator.clone(),right:ModalAttachment {component:0,shapes:vec![0.0]},
        stiffness_n_m:400.0,damping_n_s_m:0.0,rest_extension_m:0.0};
    let contact=ModalContact {left:actuator.clone(),right:receiver_port.clone(),
        law:Obstacle::new(vec![-1.0],1,1,vec![0.0001],vec![1.0],3e6,1.5,"authored-design-example".into())?};
    let contacts=vec![(contact,ModalContactConfig {max_iterations:128,maximum_force_n:10000.0,maximum_penetration_m:0.1,
        force_absolute_tolerance_n:1e-12,force_relative_tolerance:1e-13})];
    let mut network=CoupledModalSystem::new(vec![mass,receiver],vec![support],ModalCouplingConfig {
        max_modes:8,max_connections:4,max_setup_terms:8192,nyquist_guard_fraction:0.9,
        maximum_total_energy_j:1000.0,maximum_abs_pressure_pa:1e6,maximum_abs_connection_force_n:10000.0,
        solve_relative_tolerance:1e-10,energy_absolute_tolerance_j:1e-11,energy_relative_tolerance:1e-9,
    },&gate)?;
    let load=[10.0,0.0]; // 2 physical newtons times the mass's shape 5.
    network.initialize_contact_equilibrium(&load,&contacts,MultiContactConfig {
        max_contacts:4,max_sweeps:128,max_setup_terms:8192,
    },&gate)?;
    let linear=EquilibriumLinearization::new(&network,&load,&contacts,SensitivityBudget {
        max_contacts:4,max_setup_terms:8192,max_query_terms:8192,minimum_contact_margin_m:1e-8,
    },&gate)?;
    let target=DisplacementTarget {attachment:receiver_port,target_m:0.00015,scale_m:0.0001,weight:1.0};
    let result=linear.displacement_objective(&[target],&[actuator],4,&gate)?;
    // No explicit dependence on these design parameters: total derivatives are
    // the negative residual pullback. Observation-map partials remain explicit.
    println!("{{\"scope\":\"local-static-reduced-model\",\"objective\":{:.17e},\"receiver_m\":{:.17e},\"dJ_dload_N\":{:.17e},\"dJ_dgap_m\":{:.17e},\"dJ_dcontact_stiffness\":{:.17e},\"dJ_dsupport_N_per_m\":{:.17e},\"adjoint_residual\":{:.17e},\"active_contacts\":{}}}",
        result.value,result.observations_m[0],result.physical_force_gradient[0],-result.residual_pullback.contacts[0].gap,
        -result.residual_pullback.contacts[0].stiffness,-result.residual_pullback.springs[0].stiffness,
        result.adjoint.relative_residual,linear.report().active_contacts);
    Ok(())
}
fn main() {
    if let Err(error)=run() {eprintln!("contact design sensitivity refused: {error}");std::process::exit(1);}
}
