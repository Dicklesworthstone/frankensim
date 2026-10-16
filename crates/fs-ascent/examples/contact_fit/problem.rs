//! Declared two-body test rig and displacement records; no optimizer here.
use fs_couple::modal_acoustic_time::{ModalAcousticMode, ModalAcousticTimeBudget, ModalAcousticTimeModel};
use fs_couple::render::schedule::force::coupled::{ModalAttachment, ModalConnection, ModalCouplingConfig};
use fs_couple::render::schedule::force::coupled::contact::{ModalContact, ModalContactConfig};
use fs_couple::render::schedule::force::coupled::contact::multiple::MultiContactConfig;
use fs_couple::render::schedule::force::coupled::equilibrium::sensitivity::SensitivityBudget;
use fs_couple::render::schedule::force::coupled::equilibrium::sensitivity::objective::DisplacementTarget;
use fs_couple::render::schedule::force::coupled::equilibrium::sensitivity::objective::design::*;
use fs_dcontact::Obstacle;
use fs_exec::CancelGate;
use fs_math::c64::C64;

// Planted parameters: support=600 N/m, K=1.8e8 N/m^2, gap=0.2 mm.
// These displacements came from the independent quadratic physical equations,
// not from a source specimen or from the optimizer's own forward callback.
pub(super) fn synthetic_data() -> Vec<[f64;3]> {
    vec![[0.8,0.0003190067261188236,0.00006085959643287059],
         [1.6,0.0004212400757625386,0.00013472559545424769],
         [2.5,0.0005284197010700035,0.0002182948179357998]]
}

pub(super) fn parse_data(bytes: &[u8]) -> Result<Vec<[f64;3]>, String> {
    if bytes.len()>8192 {return Err("observations exceed 8192 bytes".into());}
    let text=std::str::from_utf8(bytes).map_err(|_|"observations must be UTF-8")?;
    let mut records=Vec::new();
    for (line,row) in text.lines().enumerate() {
        if records.len()==32 {return Err("at most 32 independent load cases are admitted".into());}
        let mut fields=row.split_ascii_whitespace();let mut values=[0.0;3];
        for value in &mut values {
            *value=fields.next().and_then(|s|s.parse::<f64>().ok()).filter(|v|v.is_finite())
                .ok_or_else(||format!("line {} needs finite load_N mass_displacement_m receiver_displacement_m",line+1))?;
        }
        if fields.next().is_some() {return Err(format!("extra observation field on line {}",line+1));}
        records.push(values);
    }
    if records.is_empty() {return Err("at least one complete load case is required".into());}
    Ok(records)
}

pub(super) fn build(data:&[[f64;3]], scale_m:f64, gate:&CancelGate)
    -> Result<EquilibriumDesign,Box<dyn std::error::Error>>
{
    if data.is_empty() || data.len()>32 || data.iter().flatten().any(|x|!x.is_finite()) {
        return Err("bounded complete finite observation records are required".into());
    }
    let map=|component|ModalAttachment {component,shapes:vec![if component==0 {5.0}else{1.0}]};
    let models=vec![
        ModalAcousticTimeModel::try_free_mass(48000,0.04,0.0,0.0,ModalAcousticTimeBudget::audible_reference())?,
        ModalAcousticTimeModel::try_new(48000,vec![ModalAcousticMode {angular_frequency_rad_s:100.0,
            damping_ratio:0.02,pressure_per_modal_velocity:C64::ZERO}],ModalAcousticTimeBudget::audible_reference())?,
    ];
    let spring=ModalConnection {left:map(0),right:ModalAttachment {component:0,shapes:vec![0.0]},
        stiffness_n_m:400.0,damping_n_s_m:0.0,rest_extension_m:0.0};
    let contact=ModalContact {left:map(0),right:map(1),law:Obstacle::new(vec![-1.0],1,1,vec![0.0001],vec![1.0],
        1e8,2.0,"authored-fit-template-not-calibrated".into())?};
    let cases=data.iter().enumerate().map(|(i,row)|DesignLoadCase {name:format!("experiment-{i}"),
        loads:vec![DesignLoad {attachment:map(0),force_n:row[0]}],
        targets:(0..2).map(|component|DisplacementTarget {attachment:map(component),target_m:row[component+1],
            scale_m,weight:1.0}).collect()}).collect();
    let variables=vec![
        DesignVariable {name:"support_n_m".into(),reference:400.0,scale:400.0,minimum:200.0,maximum:1000.0,
            fields:vec![DesignField::SpringStiffness(0)]},
        DesignVariable {name:"contact_n_m2".into(),reference:1e8,scale:1e8,minimum:2e7,maximum:4e8,
            fields:vec![DesignField::ContactStiffness(0)]},
        DesignVariable {name:"gap_m".into(),reference:0.0001,scale:0.0002,minimum:0.00001,maximum:0.0005,
            fields:vec![DesignField::ContactGap(0)]},
    ];
    Ok(EquilibriumDesign::new(models,vec![spring],vec![(contact,ModalContactConfig {
        max_iterations:128,maximum_force_n:10000.0,maximum_penetration_m:0.1,
        force_absolute_tolerance_n:1e-12,force_relative_tolerance:1e-13})],cases,variables,DesignBudget {
        coupling:ModalCouplingConfig {max_modes:8,max_connections:4,max_setup_terms:16384,nyquist_guard_fraction:0.9,
            maximum_total_energy_j:1000.0,maximum_abs_pressure_pa:1e6,maximum_abs_connection_force_n:10000.0,
            solve_relative_tolerance:1e-10,energy_absolute_tolerance_j:1e-11,energy_relative_tolerance:1e-9},
        contact:MultiContactConfig {max_contacts:4,max_sweeps:128,max_setup_terms:16384},
        sensitivity:SensitivityBudget {max_contacts:4,max_setup_terms:16384,max_query_terms:16384,minimum_contact_margin_m:1e-9},
        max_cases:32,max_variables:8,max_bindings:16,max_ports_per_case:8,
    },gate)?)
}
