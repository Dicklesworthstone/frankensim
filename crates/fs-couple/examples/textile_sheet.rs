//! One imported material pack drives original DKT modes and through-thickness
//! transient heat. This example supplies geometry/boundaries, not fabric laws.
use std::{collections::BTreeMap,error::Error,io::Read};
use fs_conduction::transient::deposition::*;
use fs_couple::sheet_plate::{SheetPlate,SheetPlateBudget};
use fs_material::{sheet::{ResolvedSheet,BENDING_PROPERTIES,BENDING_DIMS,MEMBRANE_PROPERTIES,MEMBRANE_DIMS,SheetLaw},
    sheet::thermal::{SheetCaloric,conductivity_at_temperature},
    state_point::{ScalarPropertyRequirement as Requirement,ScalarAdmissibility as Admit,
        MaterialPropertySelection as Pick,ResolvedMaterialStatePoint,resolve_material_state_point,
        SPECIFIC_HEAT_CAPACITY_PROPERTY as CP,SPECIFIC_HEAT_CAPACITY_DIMS,
        THERMAL_CONDUCTIVITY_PROPERTY as K,THERMAL_CONDUCTIVITY_DIMS}};
use fs_matdb::{MaterialCard,MaterialStateId,NormalizedPack,QueryPoint,ClaimId};
use fs_qty::Dims;

#[path="textile_sheet/membrane.rs"]
mod membrane;

type Result<T> = std::result::Result<T,Box<dyn Error>>;
fn query(t:f64)->Result<QueryPoint> {
    Ok(QueryPoint::new().with("ftx_characterization",0.0)?.with("T",t)?)
}
fn resolve(card:&MaterialCard,t:f64)->Result<ResolvedSheet> {
    let mut requirements=vec![Requirement::try_new("areal_density",Dims([-2,1,0,0,0,0]),Admit::StrictlyPositive)?,
        Requirement::try_new("thickness",Dims([1,0,0,0,0,0]),Admit::StrictlyPositive)?,
        Requirement::try_new(CP,SPECIFIC_HEAT_CAPACITY_DIMS,Admit::StrictlyPositive)?,
        Requirement::try_new(K,THERMAL_CONDUCTIVITY_DIMS,Admit::StrictlyPositive)?];
    for key in BENDING_PROPERTIES {requirements.push(Requirement::try_new(key,BENDING_DIMS,Admit::Finite)?);}
    for (i,key) in MEMBRANE_PROPERTIES.iter().enumerate() {
        requirements.push(Requirement::try_new(*key,if i==2 {Dims::NONE}else{MEMBRANE_DIMS},Admit::Finite)?);
    }
    Ok(ResolvedSheet::try_from_material(resolve_material_state_point(card,&query(t)?,&requirements,Pick::SingleClaimOnly)?)?)
}
fn selected(sheet:&ResolvedSheet,key:&str)->Result<ClaimId> {
    Ok(sheet.material().property(key).ok_or("missing source property")?.answer().receipt.selected)
}

#[derive(Clone)]
struct ThroughThickness<'a> {
    sheet:&'a ResolvedSheet,card:&'a MaterialCard,conductivity:ClaimId,area:f64,half_thickness:f64,
    // Keep the ACTUAL source evaluations used by the last common thermal step.
    last:Vec<ResolvedMaterialStatePoint>,
}
impl DepositionObserver for ThroughThickness<'_> {
    type Error=Box<dyn Error>;
    fn event(&mut self,_:&Event,_:&DepositionSimulation<'_>)->Result<()> {Ok(())}
    fn step(&mut self,_:f64,_:f64,_:&BTreeMap<CellId,f64>,_:&DepositionSimulation<'_>,
        check:&mut dyn FnMut()->std::result::Result<(),DepositionError>)->Result<()> {check()?;Ok(())}
    fn conductances(&mut self,thermal:&DepositionSimulation<'_>,
        check:&mut dyn FnMut()->std::result::Result<(),DepositionError>)->Result<BTreeMap<u64,f64>> {
        check()?;let mut samples=Vec::with_capacity(2);let mut resistances=[0.0;2];
        for (i,r) in resistances.iter_mut().enumerate() {
            let sample=conductivity_at_temperature(self.sheet,self.card,self.conductivity,
                thermal.temperature_k(CellId(i as u64))?)?;
            *r=0.5*self.half_thickness/sample.property(K).ok_or("missing conductivity")?.value_si();
            samples.push(sample);check()?;
        }
        let conductance=self.area/(resistances[0]+resistances[1]);
        if !conductance.is_finite() || conductance<=0.0 {return Err("unrepresentable through-thickness conductance".into());}
        self.last=samples;Ok(BTreeMap::from([(0,conductance)]))
    }
}

fn run(card:&MaterialCard,lo:f64,hi:f64,duration:f64,grain:f64)->Result<()> {
    if ![lo,hi,duration,grain].iter().all(|v|v.is_finite()) || lo<=0.0 || hi<=lo || duration<=0.0 || duration>60.0 {
        return Err("positive increasing Kelvin bounds, duration in (0,60] seconds and finite grain angle required".into());
    }
    let sheet=resolve(card,lo)?;
    let caloric=SheetCaloric::try_new(sheet.clone(),card,selected(&sheet,CP)?,[lo,hi],0.0,0.01)?;
    let mesh=fs_plate::PlateMesh::rectangle(0.1,0.1,2,2);let n=mesh.tris.len();
    let plate=SheetPlate::assemble(mesh,vec![sheet.clone();n],vec![grain;n],
        &fs_plate::PlateMesh::rectangle_boundary(2,2),
        fs_plate::AssemblyOptions {pretension:0.0,support:fs_plate::EdgeSupport::Clamped},
        SheetPlateBudget {max_nodes:9,max_triangles:8},||false)?;
    let modes=fs_plate::modes(plate.model(),(0.0,1e10),
        &fs_modal::SliceOptions {mass_diagonal_equilibration:true,..Default::default()})?;
    if modes.modes.is_empty() {return Err("no admitted sheet mode in the declared angular-frequency-squared window".into());}
    println!("# reference_mass_kg={:.12e}",plate.reference_mass_kg());
    for (i,mode) in modes.modes.iter().enumerate() {
        println!("# mode_{}_Hz={:.12e}",i+1,mode.lambda.sqrt()/std::f64::consts::TAU);
    }
    // Independent tensile/bending capabilities coexist on the SAME source.
    let elastic=sheet.membrane()?.response([0.01,0.0,0.0])?;
    println!("# membrane_energy_at_E11_0.01_J_m2={:.12e}",elastic.energy_j_m2);
    let area=0.01;let depth=0.5*sheet.thickness_m();let mass=0.5*sheet.areal_density_kg_m2()*area;
    let h0=caloric.state_at_temperature(lo)?.specific_enthalpy_j_kg();let mut events=Vec::new();
    for i in 0..2 {
        events.push(TimedEvent {time_s:0.0,event:Event::Birth(DepositedCell {
            id:CellId(i),segment_id:i,body_id:0,material_id:MaterialId(0),mass_kg:mass,specific_h_j_kg:h0,
            geometry:CellGeometry {start_m:[0.05,0.05,i as f64*depth],end_m:[0.05,0.05,(i+1) as f64*depth],
                volume_m3:area*depth,surface_area_m2:2.0*area},
            boundary:Boundary {ambient_k:lo,htc_w_m2_k:0.0,source_w:0.0,
                reservoir:Some(Reservoir {area_m2:area,conductance_w_k:10.0*area,temperature_k:if i==0 {hi}else{lo}})},
        })});
    }
    events.push(TimedEvent {time_s:0.0,event:Event::Contact(ThermalLink {id:0,a:CellId(0),b:CellId(1),
        area_m2:area,conductance_w_k:0.0})});
    let config=StepConfig {max_dt_s:0.01,energy_tolerance_j:1e-9};
    let laws=[(MaterialId(0),caloric.curve() as &dyn ThermalLaw)];
    let initial=DepositionSimulation::new(events,&laws,
        DepositionBudget {max_cells:2,max_links:1,max_events:3,max_substeps:100_000})?
        .advance_to(0.0,config,||Ok(()))?.0;
    let mut heat=initial.clone();let mut external=0.0;
    let mut source=ThroughThickness {sheet:&sheet,card,conductivity:selected(&sheet,K)?,area,half_thickness:depth,last:vec![]};
    println!("time_s,hot_half_K,cold_half_K,stored_J,boundary_heat_J,residual_J");
    for i in 1..=20 {
        let target=duration*i as f64/20.0;
        let (next,owner,audit)=heat.advance_observed(&source,target,config,||Ok(()))?;
        external+=audit.external_heat_j;
        let stored=next.energy_j()-initial.energy_j();let residual=stored-external;
        if !residual.is_finite() || residual.abs()>2e-8 {return Err("sheet energy balance failed".into());}
        println!("{:.9},{:.9},{:.9},{:.12e},{:.12e},{:.12e}",target,next.temperature_k(CellId(0))?,
            next.temperature_k(CellId(1))?,stored,external,residual);
        heat=next;source=owner;
    }
    Ok(())
}
fn main()->Result<()> {
    let args:Vec<_>=std::env::args().collect();
    if args.len()!=6 && !(args.len()==4 && args[2]=="--conductivity-at")
        && !(args.len()==9 && args[2]=="--membrane") {
        return Err("usage: textile_sheet PACK T_LOW_K T_HIGH_K SECONDS GRAIN_RAD; or textile_sheet PACK --conductivity-at T_K; or textile_sheet PACK --membrane T_K STRETCH POINT_FORCE_N GRAVITY_M_S2 GRAIN_RAD MAX_ABS_STRAIN".into());
    }
    let mut bytes=Vec::new();std::fs::File::open(&args[1])?.take(64*1024*1024+1).read_to_end(&mut bytes)?;
    if bytes.len()>64*1024*1024 {return Err("material pack exceeds 64 MiB".into());}
    let pack=NormalizedPack::from_bytes(&bytes)?;
    // Keep the whole original pack (including statistics, normalization and
    // rights) alive beside the card. No statistical blocks become constants.
    let card=MaterialCard::assemble(MaterialStateId {chemistry:pack.pack_id().into(),phase:"characterized sheet".into(),
        process:"imported fixed specimen/conditioning".into(),revision:0},pack.claims().clone(),vec![])?;
    println!("# source_pack={:?}; source_artifact={:?}",pack.pack_id(),pack.source_artifact());
    if args.len()==9 {
        return membrane::run(&card,args[3].parse()?,args[4].parse()?,args[5].parse()?,
            args[6].parse()?,args[7].parse()?,args[8].parse()?);
    }
    if args.len()==4 {
        let t:f64=args[3].parse()?;
        let answer=resolve_material_state_point(&card,&query(t)?,
            &[Requirement::try_new(K,THERMAL_CONDUCTIVITY_DIMS,Admit::StrictlyPositive)?],Pick::SingleClaimOnly)?;
        println!("thermal_conductivity_W_m_K={:.12e}",answer.property(K).ok_or("missing conductivity")?.value_si());
        return Ok(());
    }
    run(&card,args[2].parse()?,args[3].parse()?,args[4].parse()?,args[5].parse()?)
}
