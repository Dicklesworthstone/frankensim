//! Actual cooling -> thermoviscoelastic sag -> later-pass thermal contact.
//! All inputs are synthetic. The later pass is externally supported/straight;
//! contact capture does not assert mechanical adhesion, support force or healing.
use std::collections::BTreeMap;
use std::error::Error;
use fs_conduction::transient::deposition::{Boundary, CellGeometry, CellId, ContactUpdate,
    DepositedCell, DepositionBudget, DepositionError, DepositionSimulation, Event,
    MaterialId, StepConfig, ThermalLaw, ThermalLink, TimedEvent};
use fs_material::thermovisco::{Limits, ThermalKnot, ThermoMaxwell};
use fs_solid::bridge::{Bridge, Budget, Config, Support};
use fs_solid::bridge::contact::{ContactFootprint, Pass};
use fs_solid::bridge::contact::nozzle::{PathConfig, PathPoint, Profile};
use fs_solid::bridge::supports::{self, Stop};
use fs_solver::orthant::GramOrthantConfig;
use fs_solid::deposition::{self, LayerGeometry};

struct Sensible;
impl ThermalLaw for Sensible {
    fn temperature_k(&self,h: f64) -> Result<f64,DepositionError> {
        if !(0.0..=600_000.0).contains(&h) { return Err(DepositionError::MaterialDomain); }
        Ok(300.0+h/2000.0)
    }
    fn maximum_temperature_slope(&self) -> f64 { 0.0005 }
}

#[derive(Debug)]
struct Summary { sag_m: f64, flat_area_m2: f64, contact_area_m2: f64,
    sagged_old_mean_k: f64, flat_old_mean_k: f64, thermal_residual_j: f64,
    supported_sag_m: f64, support_force_n: f64, supported_contact_m2: f64, supported_old_mean_k: f64,
    nozzle_free_gap_m: f64, nozzle_supported_gap_m: f64 }

fn run() -> Result<Summary,Box<dyn Error>> {
    const N: usize = 16;
    let dx = 0.02/N as f64;
    let heat = Sensible;
    let law = ThermoMaxwell::new(1e7,vec![(9e7,0.2)],vec![
        ThermalKnot { temperature_k: 300.0,modulus_scale: 1.0,relaxation_scale: 1.0,free_strain: 0.0 },
        ThermalKnot { temperature_k: 500.0,modulus_scale: 0.1,relaxation_scale: 0.1,free_strain: 0.002 },
    ],Limits { max_strain: 0.03,max_strain_rate_s: 10.0 })?;
    let mut bridge = Bridge::new(&(0..N).map(|i| (i as u64,dx)).collect::<Vec<_>>(),
        Support::SimplySupported,Budget { max_cells: N,max_layers_per_cell: 2,max_branch_states: 128 })?
        .deposit_bonded_layer(1,LayerGeometry { area_m2: 1e-6,centroid_z_m: 0.0,thickness_m: 0.001 },
            1200.0,&law,&(0..N).map(|i| (i as u64,500.0)).collect(),&mut || false)?;
    let flat = bridge.clone();
    let mut supported = bridge.clone();
    // An ideal thermally insulating roller at a MATERIAL reference station.
    // It has no finite thermal area and supplies no adhesion or tensile force.
    let stops = [Stop { id: 1, reference_x_m: bridge.length_m()/2.0, minimum_deflection_m: 0.0 }];
    let support_config = supports::Config {
        contact: GramOrthantConfig { max_variables: 1, max_columns: N, max_setup_products: N,
            max_sweeps: 64, gradient_tolerance: 1e-13 },
        gap_tolerance_m: 1e-10, complementarity_tolerance_j: 1e-10,
    };
    let mut schedule = Vec::new();
    for (road,time_s) in [(0u64,0.0),(1,5.0)] {
        for i in 0..N {
            let id = road*100+i as u64;
            schedule.push(TimedEvent { time_s,event: Event::Birth(DepositedCell {
                id: CellId(id),segment_id: road,body_id: road,material_id: MaterialId(1),
                mass_kg: 1200.0*1e-6*dx,specific_h_j_kg: 400_000.0,
                geometry: CellGeometry { start_m: [i as f64*dx,0.0,road as f64*0.001],
                    end_m: [(i+1) as f64*dx,0.0,road as f64*0.001],
                    volume_m3: 1e-6*dx,surface_area_m2: 0.004*dx+2e-6 },
                boundary: Boundary { ambient_k: 300.0,htc_w_m2_k: 20.0,reservoir: None,source_w: 0.0 },
            }) });
        }
        for i in 1..N {
            schedule.push(TimedEvent { time_s,event: Event::Contact(ThermalLink {
                id: road*100+i as u64,a: CellId(road*100+(i-1) as u64),b: CellId(road*100+i as u64),
                area_m2: 1e-6,conductance_w_k: 0.2*1e-6/dx,
            }) });
        }
    }
    let cfg = StepConfig { max_dt_s: 0.01,energy_tolerance_j: 1e-8 };
    let mut thermal = DepositionSimulation::new(schedule,&[(MaterialId(1),&heat)],DepositionBudget {
        max_cells: 2*N,max_links: 4*N,max_events: 8*N,max_substeps: 10000,
    })?.advance_to(0.0,cfg,|| Ok(()))?.0;
    let mechanical = Config { section: deposition::Config { force_tolerance_n: 1e-7,
        moment_tolerance_nm: 1e-9,energy_tolerance_j: 1e-8,max_strain: 0.03,max_rotation_rad: 0.1 },
        maximum_deflection_m: 0.001 };
    let mut sag = 0.0;
    let mut supported_sag_m = 0.0;
    let mut support_force_n = 0.0;
    for step in 1..=50 {
        let time_s = f64::from(step)/10.0;
        let (candidate,_) = thermal.advance_to(time_s,cfg,|| Ok(()))?;
        let temperatures = (0..N).map(|i| Ok(((i as u64,1),candidate.temperature_k(CellId(i as u64))?)))
            .collect::<Result<BTreeMap<_,_>,DepositionError>>()?;
        let (shape,report) = bridge.advance_controlled(time_s,&temperatures,9.81,&[0.0;N],mechanical,&mut || false)?;
        let (propped, supports) = supported.advance_supported_controlled(time_s, &temperatures,
            9.81, &[0.0;N], &stops, mechanical, support_config, &mut || false)?;
        // All three candidates succeed before publication. The independent
        // supported/free cases share the same heat inputs, not their stresses.
        thermal = candidate; bridge = shape; sag = report.maximum_sag_m;
        supported = propped; supported_sag_m = supports.beam.maximum_sag_m;
        support_force_n = supports.reactions[0].upward_force_n;
    }
    let pass = Pass { start_x_m: 0.0,end_x_m: 0.02,start_lower_z_m: 0.0005,end_lower_z_m: 0.0005,
        width_m: 0.001,capture_gap_m: 5e-6,penetration_tolerance_m: 1e-12 };
    let actual = bridge.next_pass_contact(pass,N,&mut || false)?;
    let counterfactual = flat.next_pass_contact(pass,N,&mut || false)?;
    let propped_contact = supported.next_pass_contact(pass,N,&mut || false)?;
    // A deliberately low planned pass: supported material can obstruct a path
    // that clears the sagged road. Check the whole finite tool, not only its tip.
    // Geometry is frozen at 5 s; these times label the planned motion only.
    let path = [PathPoint { time_s: 5.0, tip_m: [0.008,0.00045] },
        PathPoint { time_s: 5.02, tip_m: [0.012,0.00045] }];
    let nozzle = Profile { flat_radius_m: 0.0002, outer_radius_m: 0.0008, flank_slope: 0.5 };
    let limits = PathConfig { max_segments: 1, max_cell_visits: N, required_clearance_m: 1e-6 };
    let free_path = bridge.nozzle_path_clearance(&path,nozzle,limits,&mut || false)?;
    let supported_path = supported.nozzle_path_clearance(&path,nozzle,limits,&mut || false)?;
    assert_eq!(free_path.first_interfering_segment,None);
    assert_eq!(supported_path.first_interfering_segment,Some(0));
    let nozzle_free_gap_m = free_path.minimum.ok_or("no free-road footprint")?.clearance.gap_m;
    let nozzle_supported_gap_m = supported_path.minimum.ok_or("no supported-road footprint")?.clearance.gap_m;
    let area = |footprint: &ContactFootprint| footprint.rows.iter().map(|r| r.reference_area_m2).sum::<f64>();
    let updates = |footprint: &ContactFootprint| footprint.rows.iter().filter(|r| r.reference_area_m2>0.0)
        .map(|row| ContactUpdate { id: 1000+row.cell_id,expected: None,replacement: Some(ThermalLink {
            id: 1000+row.cell_id,a: CellId(row.cell_id),b: CellId(100+row.cell_id),
            area_m2: row.reference_area_m2,
            // Independently declared W/(m2 K) per ORIGINAL captured area.
            conductance_w_k: 2000.0*row.reference_area_m2,
        }) }).collect::<Vec<_>>();
    let simulate = |footprint: &ContactFootprint| -> Result<(f64,f64),DepositionError> {
        let contacts = thermal.replace_contacts(&updates(footprint),|| Ok(()))?;
        assert_eq!(contacts.cells(),thermal.cells());
        let (next,audit) = contacts.advance_to(6.0,cfg,|| Ok(()))?;
        let mean = (0..N).map(|i| next.temperature_k(CellId(i as u64)))
            .collect::<Result<Vec<_>,_>>()?.iter().sum::<f64>()/N as f64;
        Ok((mean,audit.residual_j))
    };
    let (sagged_old_mean_k,residual) = simulate(&actual)?;
    let (flat_old_mean_k,flat_residual) = simulate(&counterfactual)?;
    let (supported_old_mean_k,supported_residual) = simulate(&propped_contact)?;
    assert!(sag>pass.capture_gap_m);
    assert!(area(&actual)<0.75*area(&counterfactual));
    assert!(sagged_old_mean_k<flat_old_mean_k);
    assert!(support_force_n > 0.0 && supported_sag_m < sag);
    assert!(area(&propped_contact) > area(&actual));
    assert!(supported_old_mean_k > sagged_old_mean_k);
    Ok(Summary { sag_m: sag,flat_area_m2: area(&counterfactual),contact_area_m2: area(&actual),
        sagged_old_mean_k,flat_old_mean_k,
        thermal_residual_j: residual.abs().max(flat_residual.abs()).max(supported_residual.abs()),
        supported_sag_m, support_force_n, supported_contact_m2: area(&propped_contact), supported_old_mean_k,
        nozzle_free_gap_m, nozzle_supported_gap_m })
}
fn main() -> Result<(),Box<dyn Error>> {
    let s = run()?;
    println!("sag_m={:.8e}; reference_flat_contact_m2={:.8e}; deformed_contact_m2={:.8e}; old_mean_with_sag_k={:.8}; old_mean_flat_k={:.8}; thermal_defect_j={:.3e}",
        s.sag_m,s.flat_area_m2,s.contact_area_m2,s.sagged_old_mean_k,s.flat_old_mean_k,s.thermal_residual_j);
    println!("supported_sag_m={:.8e}; roller_force_n={:.8e}; supported_contact_m2={:.8e}; supported_old_mean_k={:.8}",
        s.supported_sag_m, s.support_force_n, s.supported_contact_m2, s.supported_old_mean_k);
    println!("synthetic bounded beam/thermal example; no adhesion or contact-force qualification");
    println!("planned_nozzle_minimum_gap_m: free={:.8e}; supported={:.8e}; frozen top-surface query, not impact dynamics",
        s.nozzle_free_gap_m, s.nozzle_supported_gap_m);
    Ok(())
}
#[cfg(test)]
mod tests {
    #[test]
    fn actual_sag_changes_later_contact_and_heat_not_just_display_geometry() {
        let s = super::run().unwrap();
        assert!(s.contact_area_m2>0.0 && s.contact_area_m2<s.flat_area_m2);
        assert!(s.flat_old_mean_k-s.sagged_old_mean_k>0.1);
        assert!(s.thermal_residual_j<1e-8);
        assert!(s.support_force_n > 0.0 && s.supported_sag_m < s.sag_m);
        assert!(s.supported_contact_m2 > s.contact_area_m2);
        assert!(s.supported_old_mean_k > s.sagged_old_mean_k);
        assert!(s.nozzle_free_gap_m > 1e-6);
        assert!(s.nozzle_supported_gap_m < 0.0);
    }
}
