//! Prescribed opening/reclosure history; synthetic material and heat-transfer data.
//! The path is an input, not a prediction of an equilibrium or actuator trajectory.
use std::error::Error;
use fs_conduction::transient::deposition::{Boundary, CellGeometry, CellId, DepositedCell,
    DepositionBudget, DepositionError, DepositionSimulation, Event, MaterialId, StepConfig,
    ThermalLaw, ThermalLink, TimedEvent};
use fs_conduction::transient::enthalpy::deposition::fracture::{ContactPolicy, FractureCooling, Interface};
use fs_conduction::transient::enthalpy::deposition::fracture::prescribed::{Config, Endpoint, Frame};
use fs_material::cohesive::BilinearCohesive;

struct Sensible;
impl ThermalLaw for Sensible {
    fn temperature_k(&self, h: f64) -> Result<f64, DepositionError> {
        if !(0.0..=600_000.0).contains(&h) { return Err(DepositionError::MaterialDomain); }
        Ok(300.0+h/2000.0)
    }
    fn maximum_temperature_slope(&self) -> f64 { 0.0005 }
}

fn main() -> Result<(),Box<dyn Error>> {
    let heat = Sensible;
    let cohesive = BilinearCohesive::new(1e9,2e5,100.0,2e9,0.01)?;
    let mut events = Vec::new();
    for (id,t) in [(1,500.0),(2,300.0)] {
        events.push(TimedEvent { time_s: 0.0,event: Event::Birth(DepositedCell {
            id: CellId(id),segment_id: id,body_id: id,material_id: MaterialId(1),mass_kg: 0.001,
            specific_h_j_kg: 2000.0*(t-300.0),
            geometry: CellGeometry { start_m: [0.0,0.0,id as f64*0.001],end_m: [0.01,0.0,id as f64*0.001],
                volume_m3: 1e-6,surface_area_m2: 0.001 },
            boundary: Boundary { ambient_k: 300.0,htc_w_m2_k: 0.0,reservoir: None,source_w: 0.0 },
        }) });
    }
    let link = ThermalLink { id: 7,a: CellId(1),b: CellId(2),area_m2: 0.0002,conductance_w_k: 0.2 };
    events.push(TimedEvent { time_s: 0.0,event: Event::Contact(link) });
    let step = StepConfig { max_dt_s: 0.001,energy_tolerance_j: 1e-7 };
    let thermal = DepositionSimulation::new(events,&[(MaterialId(1),&heat)],DepositionBudget {
        max_cells: 2,max_links: 1,max_events: 3,max_substeps: 50_000,
    })?.advance_to(0.0,step,|| Ok(()))?.0;
    let initial = FractureCooling::new(thermal,vec![Interface { reference: link,point: cohesive.bonded_point(),
        closed: true,policy: ContactPolicy { open_gap_m: 0.0001,reclose_gap_m: 0.0,
            reclosed_conductance_w_k: 0.05,maximum_slip_m: 0.005 } }],1,|| Ok(()))?;
    let open = [Endpoint { interface: 7,jump_m: [0.002,0.0,0.0] }];
    let closed = [Endpoint { interface: 7,jump_m: [-0.002,0.0,0.0] }];
    let frames = [Frame { time_s: 10.0,endpoints: &open },Frame { time_s: 20.0,endpoints: &closed },
        Frame { time_s: 30.0,endpoints: &open }];
    // Contact is bonded on [0,5], open on (5,15), reclosed on [15,25.25],
    // then open again. C1=C2=2 J/K; hot=400+100 exp(-integral G dt).
    let exact_hot = 400.0+100.0*(-1.5125_f64).exp();
    println!("synthetic prescribed path; no calibration, mechanical trajectory or continuum accuracy claim");
    for dt in [0.04,0.02,0.01,0.001] {
        let config = Config { thermal: StepConfig { max_dt_s: dt,..step },event_time_tolerance_s: 1e-10,
            local_energy_tolerance_j_m2: 1e-8,max_event_trials: 2048,max_crossings: 12,max_thermal_substeps: 50_000 };
        let (final_state,report) = initial.advance_prescribed_history(&frames,config,3,|| Ok(()))?;
        let hot = final_state.thermal().temperature_k(CellId(1))?;
        println!("dt_s={dt}; hot_k={hot:.12}; analytic_hot_k={exact_hot:.12}; error_k={:.6e}; thermal_defect_j={:.6e}; steps={}",
            (hot-exact_hot).abs(),report.thermal.residual_j,report.thermal.substeps);
        assert_eq!(report.contacts.opened,2); assert_eq!(report.contacts.reclosed,1);
        assert!(final_state.interfaces()[&7].point.failed());
        assert_eq!(final_state.interfaces()[&7].point.traction_pa(),[0.0;3]);
        assert!((report.contacts.fracture_dissipation_j-0.02).abs() < 1e-10);
        assert!((final_state.thermal().energy_j()-400.0).abs() < step.energy_tolerance_j);
        if dt == 0.001 {
            assert!((hot-exact_hot).abs() < 0.01);
            for crossing in &report.crossings {
                println!("interface={}; threshold={:?}; before_s={:.12}; after_s={:.12}",
                    crossing.interface,crossing.kind,crossing.before_s,crossing.after_s);
            }
            let (prefix,_) = initial.advance_prescribed_history(&frames[..1],config,1,|| Ok(()))?;
            let (restarted,_) = prefix.clone().advance_prescribed_history(&frames[1..],config,2,|| Ok(()))?;
            assert_eq!(final_state.thermal().cells(),restarted.thermal().cells());
            assert_eq!(final_state.interfaces()[&7].point,restarted.interfaces()[&7].point);
        }
    }
    Ok(())
}
