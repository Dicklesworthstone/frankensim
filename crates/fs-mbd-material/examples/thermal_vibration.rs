//! Table-driven heat -> compatible stress -> signed vibration spectrum.
//! A fixed one-tetrahedron specimen is an integration demonstration, not a
//! spatially converged continuum, measured material or instrument model.
use std::io::{Read, Write};
use fs_alloc::{ArenaConfig,ArenaPool};
use fs_blake3::ContentHash;
use fs_exec::{Budget,CancelGate,Cx,ExecMode,StreamKey};
use fs_ga::Vec3;
use fs_material::state_point::MaterialPropertySelection;
use fs_mbd_material::{RegionBinding,ResolvedBody,VolumeMoments,resolve_body};
use fs_mbd_material::deformation::TetrahedralMass;
use fs_mbd_material::deformation::thermoelastic::ThermalSolveConfig;
use fs_mbd_material::thermal::heating::ThermalRegion;
use fs_mbd_material::thermal::heating::constrained::{ConstrainedHeatedBody,ConstrainedRegion};
use fs_modal::eigh_gen_dense;
use fs_solid::linear3::{TetAssemblyBudget,TetStaticSolveConfig};
use fs_thermal::enthalpy::{Boundary,CellGeometry,CellId,DepositionBudget,Event,MaterialId,StepConfig,TimedEvent};

#[path="thermal_vibration/source.rs"]
mod source;
use source::Source;
type Error=Box<dyn std::error::Error>;
fn invalid(message:&str)->Error { std::io::Error::new(std::io::ErrorKind::InvalidInput,message).into() }
const NODES:[[f64;3];4]=[[0.,0.,0.],[1.,0.,0.],[0.,1.,0.],[0.,0.,1.]];
const TETS:[[usize;4];1]=[[0,1,2,3]];
const FIXED:[usize;7]=[0,1,2,3,6,8,9];
const MAX_SHAPES:usize=20_000;
const MAX_OUTPUTS:usize=2_000;
fn reference(source:&Source)->Result<ResolvedBody,Error> {
    let point=source::point(source.reference_k)?;
    let geometry=VolumeMoments::tetrahedron(NODES.map(|p|Vec3::new(p[0],p[1],p[2])))?;
    Ok(resolve_body(&[RegionBinding {key:"specimen",geometry,material:&source.card,point:&point,
        selection:MaterialPropertySelection::SingleClaimOnly}],1,||false)?)
}
fn body<'a>(reference:&'a ResolvedBody,source:&'a Source,watts:f64,off:f64,dt:f64,cx:&Cx<'_>)
    ->Result<ConstrainedHeatedBody<'a>,Error> {
    let boundary=|source_w|Boundary {ambient_k:source.reference_k,htc_w_m2_k:1.,reservoir:None,source_w};
    let storage=ThermalRegion {key:"specimen",card:&source.card,caloric:&source.caloric,
        cell_id:CellId(0),material_id:MaterialId(0),segment_id:0,body_id:0,
        geometry:CellGeometry {start_m:[0.,0.,0.],end_m:[0.5,0.5,0.5],volume_m3:1./6.,
            surface_area_m2:1.5+0.5*3.0_f64.sqrt()},
        initial_specific_enthalpy_j_kg:0.,boundary:boundary(watts),
        expansion_selection:MaterialPropertySelection::SingleClaimOnly,fixed_point_m:Vec3::default()};
    let mesh=TetrahedralMass::bind(reference,&NODES,&TETS,&["specimen"],4,||false)?;
    Ok(ConstrainedHeatedBody::new(mesh,vec![ConstrainedRegion {storage,
        stiffness_model:fs_mbd_material::deformation::thermoelastic::ThermalStiffnessModel::Isotropic,
        expansion_model:fs_mbd_material::deformation::thermoelastic::ThermalExpansionModel::Isotropic,
        stiffness_selection:MaterialPropertySelection::SingleClaimOnly}],FIXED.to_vec(),
        ThermalSolveConfig {assembly:TetAssemblyBudget {maximum_nodes:4,maximum_tetrahedra:1,
            maximum_free_dofs:12,minimum_scaled_jacobian:1e-10},solve:TetStaticSolveConfig::default(),
            mesh_frame:ContentHash([37;32]),maximum_absolute_mandel_strain:0.05,maximum_nodal_displacement_m:0.1},
        vec![TimedEvent {time_s:off,event:Event::Boundary(CellId(0),boundary(0.))}],
        DepositionBudget {max_cells:1,max_links:0,max_events:2,max_substeps:MAX_SHAPES},
        StepConfig {max_dt_s:dt,energy_tolerance_j:1e-7},cx)?)
}

#[derive(Debug)]
struct Sample { time:f64, temperature:f64, energy:f64, mass:f64, volume:f64, stress:f64,
    lambda:f64, frequency:Option<f64>, operator:ContentHash }
fn sample(body:&ConstrainedHeatedBody<'_>,cx:&Cx<'_>)->Result<Sample,Error> {
    let operator=body.configuration().initial_stress_pencil(cx)?;
    let p=operator.pencil();
    let modes=eigh_gen_dense(&p.stiffness.to_dense(),&p.mass.to_dense(),p.free_dofs.len())?;
    let mode=modes.first().ok_or_else(||invalid("no free vibration mode"))?;
    if ![mode.lambda,mode.residual,mode.interval.0,mode.interval.1].iter().all(|v|v.is_finite()) {
        return Err(invalid("modal result is unrepresentable"));
    }
    let frequency=positive_frequency(mode.lambda,mode.interval.0);
    let configuration=body.configuration();
    Ok(Sample {time:body.thermal().time_s(),temperature:body.thermal().temperature_k(CellId(0))?,
        energy:body.thermal().energy_j(),mass:configuration.body().inertia().mass(),
        volume:configuration.body().regions()[0].geometry().volume_m3(),
        stress:configuration.stress().elements()[0].stress_mandel_pa()[0],lambda:mode.lambda,
        frequency,operator:operator.identity()})
}
fn positive_frequency(lambda:f64,lower:f64)->Option<f64> {
    (lambda.is_finite() && lower.is_finite() && lambda>0. && lower>0.)
        .then(||lambda.sqrt()/core::f64::consts::TAU)
}
fn row(out:&mut impl Write,s:&Sample,net_heat:f64)->Result<(),Error> {
    let frequency=s.frequency.map(|v|format!("{v:.12e}")).unwrap_or_default();
    writeln!(out,"{:.12e},{:.12e},{:.12e},{:.12e},{:.12e},{:.12e},{:.12e},{:.12e},{},{}",
        s.time,s.temperature,s.energy,net_heat,s.mass,s.volume,s.stress,s.lambda,frequency,s.operator)?;
    Ok(())
}
fn run(args:&[String],out:&mut impl Write)->Result<(),Error> {
    if args.len()!=5 {return Err(invalid("usage: thermal_vibration MATERIAL_FILE HEAT_W OFF_S END_S MAX_STEP_S"));}
    let [watts,off,end,dt]:[f64;4]=[args[1].parse()?,args[2].parse()?,args[3].parse()?,args[4].parse()?];
    if ![watts,off,end,dt].iter().all(|v|v.is_finite()) || watts<0. || off<=0. || end<=off || dt<=0.
        || end/dt>MAX_OUTPUTS as f64 {
        return Err(invalid("finite nonnegative heat, 0<off<end, positive step and at most 2000 output intervals required"));
    }
    let mut text=String::new();
    std::fs::File::open(&args[0])?.take(source::MAX_BYTES as u64+1).read_to_string(&mut text)?;
    let source=Source::parse(&text,&args[0])?;
    let reference=reference(&source)?;
    let gate=CancelGate::new();
    ArenaPool::new(ArenaConfig::default()).scope(|arena|->Result<(),Error> {
        let cx=Cx::new(&gate,arena,StreamKey {seed:0,kernel_id:37,tile:0,iteration:0},Budget::INFINITE,ExecMode::Deterministic);
        let mut body=body(&reference,&source,watts,off,dt,&cx)?;
        writeln!(out,"# material_card={}",source.card.content_hash())?;
        writeln!(out,"time_s,temperature_k,thermal_energy_j,net_external_heat_j,mass_kg,volume_m3,sigma_xx_pa,lowest_lambda_per_s2,lowest_positive_frequency_hz,operator_identity")?;
        row(out,&sample(&body,&cx)?,0.)?;
        let (mut net_heat,mut remaining)=(0.,MAX_SHAPES);
        for i in 1..=(end/dt).ceil().max(1.) as usize {
            let target=(i as f64*dt).min(end);
            let (next,report)=body.advance_to(target,StepConfig {max_dt_s:dt,energy_tolerance_j:1e-7},remaining,&cx)?;
            remaining=remaining.checked_sub(report.shape_solves).ok_or_else(||invalid("whole-run shape allowance exhausted"))?;
            net_heat+=report.thermal.external_heat_j;
            if !net_heat.is_finite() {return Err(invalid("net heat is unrepresentable"));}
            let observation=sample(&next,&cx)?;
            row(out,&observation,net_heat)?;
            body=next;
        }
        Ok(())
    })
}
fn main() {
    if let Err(error)=run(&std::env::args().skip(1).collect::<Vec<_>>(),&mut std::io::stdout().lock()) {
        eprintln!("thermal_vibration: {error}");std::process::exit(1);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    fn with_cx<R>(f:impl FnOnce(&Cx<'_>)->R)->R {
        let gate=CancelGate::new();ArenaPool::new(ArenaConfig::default()).scope(|arena|
            f(&Cx::new(&gate,arena,StreamKey {seed:0,kernel_id:37,tile:0,iteration:0},Budget::INFINITE,ExecMode::Deterministic)))
    }
    #[test]
    fn changing_only_the_material_file_cp_changes_the_solved_spectrum() {
        let a=Source::parse(include_str!("thermal_vibration/cp10.txt"),"same source label").unwrap();
        let b=Source::parse(include_str!("thermal_vibration/cp20.txt"),"same source label").unwrap();
        let mut observed=Vec::new();
        for source in [&a,&b] {
            let reference=reference(source).unwrap();let old=with_cx(|cx|body(&reference,source,100.,1.,0.05,cx)).unwrap();
            let (hot,_)=with_cx(|cx|old.advance_to(1.,StepConfig {max_dt_s:0.05,energy_tolerance_j:1e-7},100,cx)).unwrap();
            let hot_sample=with_cx(|cx|sample(&hot,cx)).unwrap();
            let (cool,_)=with_cx(|cx|hot.advance_to(2.,StepConfig {max_dt_s:0.05,energy_tolerance_j:1e-7},100,cx)).unwrap();
            let cooled=with_cx(|cx|sample(&cool,cx)).unwrap();
            assert!(cooled.temperature<hot_sample.temperature && cooled.lambda>hot_sample.lambda);
            assert_eq!(hot_sample.mass,1.);observed.push(hot_sample);
        }
        assert!(observed[0].temperature>observed[1].temperature);
        assert!(observed[0].lambda<observed[1].lambda);
    }
    #[test]
    fn malformed_tables_and_unresolved_frequency_are_not_silently_repaired() {
        let good=include_str!("thermal_vibration/cp10.txt");
        for bad in [good.replace("300 6","300 NaN"),good.replace("0.25","0.5"),
            good.replace("400 6","299 6"),good.replace("v1","v2")] {
            assert!(Source::parse(&bad,"invalid").is_err());
        }
        assert!(positive_frequency(-1.,-1.).is_none());
        assert!(positive_frequency(0.,0.).is_none());
        assert!(positive_frequency(1.,-0.1).is_none());
        assert_eq!(positive_frequency(4.,3.),Some(1./core::f64::consts::PI));
    }
    #[test]
    fn command_reads_the_table_and_emits_actual_heating_and_cooling() {
        let path=std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("examples/thermal_vibration/cp10.txt");
        let args=vec![path.to_str().unwrap().to_owned(),"100".into(),"1".into(),"2".into(),"0.1".into()];
        let mut out=Vec::new();run(&args,&mut out).unwrap();
        let text=String::from_utf8(out).unwrap();
        let rows:Vec<Vec<f64>>=text.lines().filter(|s|!s.starts_with('#')).skip(1)
            .map(|s|s.split(',').take(9).map(|x|x.parse().unwrap()).collect()).collect();
        assert_eq!(rows.len(),21);assert_eq!(rows[0][0],0.);assert_eq!(rows[20][0],2.);
        assert!(rows[10][1]>rows[20][1]);assert!(rows[10][7]<rows[20][7]);
        for r in rows {assert!((r[2]-r[3]).abs()<1e-7);assert_eq!(r[4],1.);}
        let mut bad=args.clone();bad[2]="3".into();let mut out=Vec::new();
        assert!(run(&bad,&mut out).is_err());assert!(out.is_empty());
    }

}
