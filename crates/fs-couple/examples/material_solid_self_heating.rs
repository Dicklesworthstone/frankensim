//! Same-card inertia, elasticity, Kelvin loss and finite caloric storage.
//! A deliberately small, insulated tetrahedron with one free axial component.
//! --demo is manufactured data, not a measured material. --pack never falls back.
use std::collections::BTreeMap;
use std::io::{Read, Write};
use fs_alloc::{ArenaConfig, ArenaPool};
use fs_blake3::ContentHash;
use fs_couple::tensor_material::{InitialMemory, TensorMaterialBudget, solid::TensorSolidBody};
use fs_couple::tensor_material::solid::thermal::{ElementThermalCell, SolidThermalConfig};
use fs_evidence::ValidityDomain;
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_matdb::{ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId,
    NormalizedMaterialCardPack, PropertyClaim, PropertyKey, PropertyValue, Provenance,
    QueryPoint, SelectionPolicy, UncertaintyModel};
use fs_material::phase::{HeatCapacitySourceQuery, ResolvedHeatCapacityCurve, SolidLiquidPhase};
use fs_material::state_point::{MaterialPropertySelection, MaterialStatePointError,
    ResolvedMaterialStatePoint, ScalarAdmissibility, ScalarPropertyRequirement,
    SPECIFIC_HEAT_CAPACITY_DIMS, resolve_material_state_point};
use fs_material::thermal::source::{ThermalMaterialError, ThermalMaterialPoint};
use fs_material::visco::tensor::TensorProny;
use fs_qty::{Density, Dims, Pressure, Temperature};
use fs_solid::linear3::{TetAssemblyBudget, TetElasticMaterial,
    visco::{SolidMeshRef, SolidProjectionBudget, TetViscoField, TetViscoMaterial,
        ViscoElements, ViscoStateBudget}};
use fs_thermal::enthalpy::*;

type Result<T> = std::result::Result<T, Box<dyn std::error::Error>>;
const NODES: [[f64;3];4] = [[0.,0.,0.],[1.,0.,0.],[0.,1.,0.],[0.,0.,1.]];
const TETS: [[usize;4];1] = [[0,1,2,3]];
const FIXED: [usize;11] = [0,1,2,3,4,5,6,7,8,9,10];
const TIME_DIMS: Dims = Dims([0,0,1,0,0,0]);
const DEMO_TIME: &str = "demonstration_kelvin_time";
const INITIAL_Q_M: f64 = 0.001;
const INITIAL_V_M_S: f64 = 0.1;

fn point(t: f64) -> Result<QueryPoint> { Ok(QueryPoint::new().with("T", t)?) }
fn requirements(time_property: &str) -> Result<Vec<ScalarPropertyRequirement>> {
    Ok(vec![
        ScalarPropertyRequirement::try_new("density",Density::DIMS,ScalarAdmissibility::StrictlyPositive)?,
        ScalarPropertyRequirement::try_new("young_modulus",Pressure::DIMS,ScalarAdmissibility::StrictlyPositive)?,
        ScalarPropertyRequirement::try_new("poisson_ratio",Dims::NONE,
            ScalarAdmissibility::OpenInterval {lower:-1.,upper:0.5})?,
        ScalarPropertyRequirement::try_new(time_property,TIME_DIMS,ScalarAdmissibility::NonNegative)?,
    ])
}
fn value(source: &ResolvedMaterialStatePoint, name: &str) -> Result<f64> {
    Ok(source.property(name).ok_or("incomplete mechanical material bundle")?.value_si())
}

#[derive(Debug)]
enum SourceError {
    Thermal(DepositionError), Material(MaterialStatePointError), Static(ThermalMaterialError),
    Invalid(&'static str),
}
impl From<DepositionError> for SourceError {
    fn from(e:DepositionError)->Self {Self::Thermal(e)}
}
impl std::fmt::Display for SourceError {
    fn fmt(&self,f:&mut std::fmt::Formatter<'_>)->std::fmt::Result {
        match self {
            Self::Thermal(e)=>write!(f,"thermal source: {e}"),
            Self::Material(e)=>write!(f,"mechanical source: {e}"),
            Self::Static(e)=>write!(f,"static mechanical source: {e}"),
            Self::Invalid(s)=>f.write_str(s),
        }
    }
}
impl std::error::Error for SourceError {}

/// Recheck pinned original claims at EVERY actual local endpoint. This model
/// admits no temperature dependence in mechanical coefficients; even a small
/// changed value refuses instead of silently retaining a stale coefficient.
#[derive(Clone)]
struct FrozenSources<'a> {
    card: &'a MaterialCard,
    requirements: Vec<ScalarPropertyRequirement>,
    selection: MaterialPropertySelection,
    reference: ResolvedMaterialStatePoint,
    last_identity: ContentHash,
}
impl FrozenSources<'_> {
    fn admit(&mut self, thermal: &DepositionSimulation<'_>) -> std::result::Result<(),SourceError> {
        let t=thermal.temperature_k(CellId(30))?;
        let context=QueryPoint::new().with("T",t)
            .map_err(|_|SourceError::Invalid("invalid caloric temperature query"))?;
        let current=resolve_material_state_point(self.card,&context,&self.requirements,self.selection.clone())
            .map_err(SourceError::Material)?;
        let thermal_point=ThermalMaterialPoint {material:self.card,context:&context,
            temperature_axis:"T",liquid_fraction_axis:None};
        for original in self.reference.properties() {
            let property=current.property(original.requirement().name())
                .ok_or(SourceError::Invalid("missing pinned mechanical property"))?;
            thermal_point.admit_static(property).map_err(SourceError::Static)?;
            if property.value_si().to_bits()!=original.value_si().to_bits() {
                return Err(SourceError::Invalid("temperature-dependent mechanical source requires a coupled coefficient/history update; frozen model refused"));
            }
        }
        self.last_identity=current.identity();
        Ok(())
    }
}
impl DepositionObserver for FrozenSources<'_> {
    type Error=SourceError;
    fn sources(&mut self,thermal:&DepositionSimulation<'_>,checkpoint:&mut dyn FnMut()->std::result::Result<(),DepositionError>)
        ->std::result::Result<BTreeMap<CellId,f64>,Self::Error> {
        checkpoint()?;self.admit(thermal)?;Ok(BTreeMap::new())
    }
    fn event(&mut self,_:&Event,thermal:&DepositionSimulation<'_>)->std::result::Result<(),Self::Error> {self.admit(thermal)}
    fn step(&mut self,_:f64,_:f64,_:&BTreeMap<CellId,f64>,thermal:&DepositionSimulation<'_>,
        checkpoint:&mut dyn FnMut()->std::result::Result<(),DepositionError>)->std::result::Result<(),Self::Error> {
        checkpoint()?;self.admit(thermal)
    }
}

#[derive(Clone,Debug,PartialEq)]
struct Sample {
    time_s:f64, displacement_m:f64, velocity_m_s:f64, temperature_k:f64,
    mechanical_energy_j:f64, caloric_energy_change_j:f64, material_heat_j:f64,
    combined_residual_j:f64,
}
struct Run {
    card:ContentHash, reference:ContentHash, final_material:ContentHash,
    reference_mass_kg:f64, samples:Vec<Sample>,
}

fn simulate(card:&MaterialCard,time_property:&str,steps:usize,dt:f64,gate:&CancelGate)->Result<Run> {
    if !(1..=4096).contains(&steps) || !dt.is_finite() || dt<=0. || !(dt*steps as f64).is_finite() {
        return Err("1..=4096 steps and a positive finite representable interval required".into());
    }
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        let cx=Cx::new(gate,arena,StreamKey {seed:87,kernel_id:88,tile:0,iteration:0},Budget::INFINITE,ExecMode::Deterministic);
        cx.checkpoint().map_err(|_|"self-heating cancelled")?;
        let requirements=requirements(time_property)?;
        let reference=resolve_material_state_point(card,&point(300.)?,&requirements,MaterialPropertySelection::SingleClaimOnly)?;
        let rho=value(&reference,"density")?;
        let elastic=TetElasticMaterial::try_new(rho,value(&reference,"young_modulus")?,
            value(&reference,"poisson_ratio")?,reference.identity().0)?;
        let time=value(&reference,time_property)?;
        let stiffness=*elastic.stiffness_mandel_pa();
        // Explicit proportional isotropic Kelvin model, not a modal loss fit.
        let viscosity=if time==0. {None} else {Some(stiffness.map(|row|row.map(|v|time*v)))};
        let material=TetViscoMaterial {law:TensorProny::new(stiffness,viscosity,vec![],(0.,1e5))?,
            density_kg_m3:rho,material_state_identity:reference.identity().0};
        let source=ViscoElements::new(SolidMeshRef::Tet4 {nodes_m:&NODES,elements:&TETS},
            TetViscoField::Uniform(&material),&FIXED,TetAssemblyBudget::standard(),
            ViscoStateBudget {max_internal_variables:1},&cx)?;
        let volume=source.equilibrium().element_volumes_m3[0];
        let mass=source.equilibrium().total_mass_kg;
        if !mass.is_finite() || mass<=0. || (mass-rho*volume).abs()>512.*f64::EPSILON*mass {
            return Err("native reference mass and source-density volume disagree".into());
        }
        let projection=source.project_basis(&[vec![1.]],SolidProjectionBudget {max_coordinates:1,
            max_basis_scalars:1,max_material_points:1,max_scalar_terms:4096},&cx)?;
        let body=TensorSolidBody::new(projection,TensorMaterialBudget {max_state_coordinates:2,max_material_points:1},gate)?;
        let caloric=ResolvedHeatCapacityCurve::try_new(card,HeatCapacitySourceQuery {
            heat_capacity:&PropertyKey::new("specific_heat_capacity",SPECIFIC_HEAT_CAPACITY_DIMS),
            density:&PropertyKey::new("density",Density::DIMS),temperature_axis:"T",
            lower:&point(300.)?,upper:&point(350.)?,selection:SelectionPolicy::SingleClaimOnly,
        },SolidLiquidPhase::Solid,0.,1e-10)?;
        let thermal_config=StepConfig {max_dt_s:dt,energy_tolerance_j:1e-10};
        let mut thermal=DepositionSimulation::new(vec![TimedEvent {time_s:0.,event:Event::Birth(DepositedCell {
            id:CellId(30),segment_id:30,body_id:1,material_id:MaterialId(0),mass_kg:mass,
            specific_h_j_kg:caloric.state_at_temperature(300.)?.specific_enthalpy_j_kg(),
            geometry:CellGeometry {start_m:NODES[0],end_m:NODES[3],volume_m3:volume,
                surface_area_m2:1.5+0.5*3.0_f64.sqrt()},
            boundary:Boundary {ambient_k:300.,htc_w_m2_k:0.,reservoir:None,source_w:0.},
        })}],&[(MaterialId(0),caloric.curve() as &dyn ThermalLaw)],
            DepositionBudget {max_cells:1,max_links:0,max_events:1,max_substeps:4})?
            .advance_to(0.,thermal_config,||if gate.is_requested(){Err(DepositionError::Cancelled)}else{Ok(())})?.0;
        let mut observer=FrozenSources {card,requirements,selection:MaterialPropertySelection::PinnedByProperty(
            reference.properties().iter().map(|p|(p.requirement().name().to_owned(),p.answer().receipt.selected)).collect()),
            last_identity:reference.identity(),reference:reference.clone()};
        observer.admit(&thermal)?;
        let cells=[ElementThermalCell {cell:CellId(30),curve:caloric.curve(),temperature_interval_k:[300.,350.]}];
        let mut history=body.history(&[INITIAL_Q_M],&[INITIAL_V_M_S],InitialMemory::Relaxed,gate)?;
        let initial_energy=history.material_history().stored_energy_j();
        let initial_caloric=thermal.energy_j();
        let mut samples=vec![Sample {time_s:0.,displacement_m:INITIAL_Q_M,velocity_m_s:INITIAL_V_M_S,
            temperature_k:300.,mechanical_energy_j:initial_energy,caloric_energy_change_j:0.,material_heat_j:0.,combined_residual_j:0.}];
        for _ in 0..steps {
            let target=thermal.time_s()+dt;
            let candidate=history.advance_thermal_candidate(&thermal,&observer,target,&[[0.;3];4],&cells,
                SolidThermalConfig {thermal:thermal_config,energy_tolerance_j:1e-9},gate)?;
            let (next,next_thermal,next_observer,_)=candidate.into_parts();
            let mechanical=next.material_history().stored_energy_j();
            let caloric_change=next_thermal.energy_j()-initial_caloric;
            let residual=mechanical+caloric_change-initial_energy;
            if !residual.is_finite() || residual.abs()>(samples.len() as f64)*1e-9 {
                return Err("whole self-heating trajectory energy did not close".into());
            }
            samples.push(Sample {time_s:next_thermal.time_s(),displacement_m:next.material_history().displacement_m()[0],
                velocity_m_s:next.material_history().velocity_m_s()[0],temperature_k:next_thermal.temperature_k(CellId(30))?,
                mechanical_energy_j:mechanical,caloric_energy_change_j:caloric_change,
                material_heat_j:next.material_history().dissipated_heat_j(),combined_residual_j:residual});
            history=next;thermal=next_thermal;observer=next_observer;
        }
        cx.checkpoint().map_err(|_|"self-heating cancelled")?;
        Ok(Run {card:card.content_hash(),reference:reference.identity(),final_material:observer.last_identity,
            reference_mass_kg:mass,samples})
    })
}

fn demo(time_scale:f64,cp_scale:f64)->Result<MaterialCard> {
    if !time_scale.is_finite() || time_scale<0. || !cp_scale.is_finite() || cp_scale<=0. {
        return Err("finite nonnegative viscosity scale and positive heat-capacity scale required".into());
    }
    demo_data(2.,120.,120.,Some(0.01*time_scale),cp_scale)
}
fn demo_data(rho:f64,young:f64,young_end:f64,time:Option<f64>,cp:f64)->Result<MaterialCard> {
    let mut claims=ClaimSet::new();
    let mut rows=vec![("density",Density::DIMS,rho,rho),("young_modulus",Pressure::DIMS,young,young_end),
        ("poisson_ratio",Dims::NONE,0.,0.),("specific_heat_capacity",SPECIFIC_HEAT_CAPACITY_DIMS,cp,cp)];
    if let Some(time)=time {rows.push((DEMO_TIME,TIME_DIMS,time,time));}
    for (name,dims,a,b) in rows {
        claims.insert_claim(PropertyClaim {key:PropertyKey::new(name,dims),
            value:PropertyValue::Curve {abscissa:"T".into(),abscissa_dims:Temperature::DIMS,
                knots:vec![(300.,a),(350.,b)],dims},
            validity:ValidityDomain::unconstrained().with("T",300.,350.),
            uncertainty:UncertaintyModel::Unstated,interpolation:InterpolationPolicy::LinearInside,
            observations:vec![],provenance:Provenance {source:"manufactured self-heating fixture, not measurements".into(),
                license:"CC0-1.0".into(),artifact:None}})?;
    }
    Ok(MaterialCard::assemble(MaterialStateId {chemistry:"synthetic".into(),phase:"solid".into(),
        process:"manufactured integration fixture".into(),revision:0},claims,vec![])? )
}
fn main()->Result<()> {
    let args=std::env::args().skip(1).collect::<Vec<_>>();
    if args.len()!=6 {
        return Err("usage: material_solid_self_heating --demo VISCOSITY_SCALE CP_SCALE STEPS DT_S OUTPUT.csv | --pack MATERIAL.matpack KELVIN_TIME_PROPERTY STEPS DT_S OUTPUT.csv".into());
    }
    let (card,key)=match args[0].as_str() {
        "--demo"=>(demo(args[1].parse()?,args[2].parse()?)?,DEMO_TIME),
        "--pack"=>{
            let mut bytes=Vec::new();
            std::fs::File::open(&args[1])?.take(16*1024*1024+1).read_to_end(&mut bytes)?;
            if bytes.len()>16*1024*1024 {return Err("material pack exceeds 16 MiB".into());}
            (NormalizedMaterialCardPack::from_bytes(&bytes)?.card().clone(),args[2].as_str())
        }
        _=>return Err("select --demo or --pack explicitly".into()),
    };
    let result=simulate(&card,key,args[3].parse()?,args[4].parse()?,&CancelGate::new_clock_free())?;
    // Compute/admit the entire trajectory before creating any output. Existing
    // output paths are never overwritten. File I/O itself is not transactional.
    let mut out=std::io::BufWriter::new(std::fs::OpenOptions::new().write(true).create_new(true).open(&args[5])?);
    writeln!(out,"# card={:?}; reference={:?}; final_material={:?}; reference_mass_kg={:.17e}; explicit isotropic proportional Kelvin model; fixed reference tetrahedron/supports; no thermal expansion; no Maxwell branches; 300..350 K fully-solid assumption; no empirical qualification",
        result.card,result.reference,result.final_material,result.reference_mass_kg)?;
    writeln!(out,"time_s,displacement_m,velocity_m_s,temperature_K,mechanical_energy_J,caloric_energy_change_J,material_heat_J,combined_residual_J")?;
    for s in result.samples {
        writeln!(out,"{:.17e},{:.17e},{:.17e},{:.17e},{:.17e},{:.17e},{:.17e},{:.17e}",s.time_s,s.displacement_m,
            s.velocity_m_s,s.temperature_k,s.mechanical_energy_j,s.caloric_energy_change_j,s.material_heat_j,s.combined_residual_j)?;
    }
    out.flush()?;
    Ok(())
}

#[cfg(test)]
#[path="material_solid_self_heating/tests.rs"]
mod tests;
