//! Geometry-derived fluid reaction, not a pressure-only output effect.
//! One existing modal BEM batch gives both Z = G^T A P and every receiver's
//! pressure/velocity transfer. exp(-i omega t): the force opposing motion
//! adds -i omega Z to the same string/board dynamic stiffness.
use super::{bridge_response::{BridgeResponse,Response}, exterior_geometry::{Boundary,Specification,MAX_PANELS,MAX_HARMONIC_RECEIVERS,ReceiverSet,SourceSweep}};
#[cfg(test)]
use fs_bem::{helmholtz::{self,Formulation},radiation_policy::GeometryPolicy};
use fs_math::c64::C64;
use std::{f64::consts::TAU,fmt::Write};

fn finite(v:C64)->bool {v.re.is_finite() && v.im.is_finite()}

pub struct LoadingSample {
    /// Row-major generalized force / generalized velocity, in Bank's basis.
    pub impedance:Vec<C64>,
    /// Receiver-major complex Pa per unit generalized velocity.
    pub receiver_transfer:Vec<Vec<C64>>,
    pub minimum_ppw:f64,
    pub condition_lower_bound:f64,
}
impl LoadingSample {
    pub fn pressure(&self,response:&Response,w:f64)->Result<Vec<C64>,String> {
        if !w.is_finite() || w<=0. || self.receiver_transfer.iter()
            .any(|r|r.len()!=response.board_displacement.len()) {
            return Err("pressure projection has a different board basis or frequency".into());
        }
        let out:Vec<C64>=self.receiver_transfer.iter().map(|row|row.iter()
            .zip(&response.board_displacement).fold(C64::ZERO,|s,(h,q)|s+*h * *q*C64::new(0.,-w))).collect();
        if out.iter().any(|p|!finite(*p)) {return Err("loaded receiver pressure overflow".into());}
        Ok(out)
    }
}

/// Only the current frequency's r modal fields are solved, not an n-panel
/// impedance matrix. Quadrature appears on the FORCE side exactly once.
/// The supplied boundary must already be transformed by `loaded(model.bank())`.
/// No symmetrization, passive projection, loss insertion or resonance fitting.
pub fn sample(boundary:&Boundary,spec:&Specification,w:f64)->Result<LoadingSample,String> {
    LoadingSweep::new(boundary,spec,&[w])?.sample(w)
}

/// Immutable geometry, modal velocity fields and receiver admission shared by
/// a harmonic sweep or a loaded-playback fit. Every frequency still gets its
/// own operator, factorization and bounded receiver quadrature.
pub struct LoadingSweep<'a> {
    boundary:&'a Boundary,
    spec:&'a Specification,
    receivers:ReceiverSet<'a>,
    source:SourceSweep<'a>,
    fields:Vec<Vec<C64>>,
}
impl<'a> LoadingSweep<'a> {
    pub fn new(boundary:&'a Boundary,spec:&'a Specification,omega:&[f64])->Result<Self,String> {
        let count=boundary.weights.len();let panels=boundary.surface.areas().len();
        if !(1..=super::linear::MAX_BOARD_MODES).contains(&count) || !(1..=MAX_PANELS).contains(&panels)
            || boundary.weights.iter().any(|r|r.len()!=panels || r.iter().any(|v|!v.is_finite()))
            || omega.iter().any(|w|!w.is_finite() || *w<=0.
                || *w/TAU<spec.band_hz.0-1e-10 || *w/TAU>spec.band_hz.1+1e-10)
            || !spec.medium.density.is_finite() || spec.medium.density<=0.
            || !spec.medium.sound_speed.is_finite() || spec.medium.sound_speed<=0.
            || !spec.min_ppw.is_finite() || spec.min_ppw<6.
            || !(1..=MAX_HARMONIC_RECEIVERS).contains(&spec.receivers.len()) {
            return Err("invalid complete radiation-loading basis, medium or requested frequency".into());
        }
        let receivers=ReceiverSet::for_spec(boundary,spec)?;
        let source=SourceSweep::new(boundary,omega,spec.medium,1)?;
        let fields=boundary.weights.iter().map(|r|r.iter().map(|v|C64::new(*v,0.)).collect()).collect();
        Ok(Self {boundary,spec,receivers,source,fields})
    }
    pub fn delays_s(&self)->&[f64] {self.receivers.delays_s()}
    pub fn sample(&self,w:f64)->Result<LoadingSample,String> {
    let boundary=self.boundary;let spec=self.spec;
    let count=boundary.weights.len();let panels=boundary.surface.areas().len();
    if !w.is_finite() || w<=0. || w/TAU<spec.band_hz.0-1e-10 || w/TAU>spec.band_hz.1+1e-10 {
        return Err("radiation-loading frequency is outside the declared band".into());
    }
    let k=w/spec.medium.sound_speed;
    let evaluation=self.receivers.prepare(k)?;
    let refs:Vec<&[C64]>=self.fields.iter().map(Vec::as_slice).collect();
    let formulation=self.source.formulation(k).map_err(|e|e.to_string())?;
    let solutions=self.source.solve_batch(k,spec.medium,&refs)
        .map_err(|e|format!("radiation-load solve at {} Hz ({formulation:?}): {e}",w/TAU))?;
    let mut impedance=vec![C64::ZERO;count*count];
    let mut receiver_transfer=vec![vec![C64::ZERO;count];spec.receivers.len()];
    let mut minimum_ppw=f64::INFINITY;let mut condition_lower_bound=0.0_f64;
    for (j,solution) in solutions.iter().enumerate() {
        if solution.pressure.len()!=panels || !solution.panels_per_wavelength.is_finite()
            || solution.panels_per_wavelength<spec.min_ppw || !solution.condition_lower_bound.is_finite()
            || !solution.radiated_power_roundoff_interval.1.is_finite()
            || solution.radiated_power_roundoff_interval.1<0. {
            return Err(format!("radiation load at {} Hz ({formulation:?}) has unresolved power, conditioning or wavelength resolution: ppw={}, condition lower bound={}, power interval={:?}",
                w/TAU,solution.panels_per_wavelength,solution.condition_lower_bound,
                solution.radiated_power_roundoff_interval));
        }
        minimum_ppw=minimum_ppw.min(solution.panels_per_wavelength);
        condition_lower_bound=condition_lower_bound.max(solution.condition_lower_bound);
        for (i,row) in boundary.weights.iter().enumerate() {
            impedance[i*count+j]=row.iter().zip(boundary.surface.areas()).zip(&solution.pressure)
                .fold(C64::ZERO,|sum,((shape,area),p)|sum+p.scale(shape*area));
        }
        let pressure=evaluation.pressure(solution)?;
        for (row,p) in receiver_transfer.iter_mut().zip(pressure) {row[j]=p;}
    }
    if impedance.iter().chain(receiver_transfer.iter().flatten()).any(|v|!finite(*v)) {
        return Err("projected radiation impedance or receiver transfer is nonfinite".into());
    }
    Ok(LoadingSample {impedance,receiver_transfer,minimum_ppw,condition_lower_bound})
    }
}

/// Unit peak bridge force; every other course is a passive, undamped-by-key
/// string termination, never removed from the mechanical system. The one-way
/// comparison retains the same air transfer but omits its mechanical reaction.
/// Input/wood/string/radiation/cavity powers apply to the COUPLED columns only.
/// Completion is atomic at the string level: no partial CSV on a failed sweep.
pub fn sweep(boundary:&Boundary,model:&BridgeResponse,spec:&Specification,drive:u8)->Result<String,String> {
    model.bridge_row(drive)?;
    if boundary.weights.len()!=model.bank().board_count {return Err("radiation and string bank bases differ".into());}
    let omega=spec.omega();let prepared=LoadingSweep::new(boundary,spec,&omega)?;
    let cavity_column=if model.has_cavity() {",cavity_w"} else {""};
    let mut out=format!("# radiation-loaded bridge admittance; exp(-i omega t); unit peak force 1 N at key {drive}\n# acoustic source: {}\n# all {} scale courses and {} retained string coordinates; no hammer or key-damper contact\n# bridge values: m/s/N; receiver values: Pa/N; power columns: cycle-average W for 1 N peak\n# one_way columns use the SAME air transfer on mechanics without exterior radiation reaction; a selected internal cavity remains in both solutions\n# no fitted transfer, time-domain radiation feedback, full-band convergence or measured Steinway claim\nfrequency_hz,observable,index,real,imag,one_way_real,one_way_imag,input_w,wood_w,string_w,radiation_w,power_defect_w,backward_error,panels_per_wavelength,condition_lower_bound{cavity_column}\n",
        spec.source,model.keys().len(),model.bank().modes.len());
    for w in omega {
        let hz=w/TAU;let field=prepared.sample(w)?;
        let loaded=model.solve(hz,drive,C64::ONE,Some(&field.impedance))?;
        let unreacted=model.solve(hz,drive,C64::ONE,None)?;
        let pressure=field.pressure(&loaded,w)?;let reference=field.pressure(&unreacted,w)?;
        let mut row=|kind:&str,index:usize,value:C64,one_way:C64| {
            write!(out,"{hz:.17e},{kind},{index},{:.17e},{:.17e},{:.17e},{:.17e},{:.17e},{:.17e},{:.17e},{:.17e},{:.17e},{:.17e},{:.17e},{:.17e}",
                value.re,value.im,one_way.re,one_way.im,loaded.input_w,loaded.board_loss_w,loaded.string_loss_w,
                loaded.radiation_w,loaded.power_defect_w,loaded.backward_error,field.minimum_ppw,field.condition_lower_bound).unwrap();
            if model.has_cavity() {write!(out,",{:.17e}",loaded.cavity_loss_w).unwrap();}
            writeln!(out).unwrap();
        };
        for ((&key,&v),&base) in model.keys().iter().zip(&loaded.bridge_velocity).zip(&unreacted.bridge_velocity) {
            row("bridge",usize::from(key),v,base);
        }
        for (i,(&p,&base)) in pressure.iter().zip(&reference).enumerate() {row("receiver",i,p,base);}
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;
    use super::super::{board_geometry::motion::MotionSurface,exterior_geometry::tests as geometry,geometry::Course,linear::BoardMode};
    fn fixture()->(Boundary,BridgeResponse,Specification) {
        let spec=Specification::read(&format!("{}receiver-m,0.05,0.05,-1\n",geometry::specification())).unwrap();
        let mesh=fs_plate::ShellMesh::new(vec![[0.,0.,0.],[0.1,0.,0.],[0.1,0.1,0.],[0.,0.1,0.]],vec![[0,1,2],[0,2,3]]).unwrap();
        // A constant kinematic field is explicit here, not a plate eigenpair.
        let motion=MotionSurface::new(mesh,vec![vec![[0.,0.,1.,0.,0.,0.];4]]).unwrap();
        let mut board=BoardMode {frequency_hz:170.,damping_ratio:0.01,bridge:[0.;88],volume:0.01};board.bridge[48]=1.;
        let c=super::super::geometry::demonstration_scale().unwrap()[48];
        let model=BridgeResponse::new(&[Course {unison:1,duplex_length_m:0.,..c}],&[board],192_000,21_600.,4,true).unwrap();
        let obj=geometry::box_obj("skin",[0.,0.,-0.01],[0.1,0.1,0.02]);
        let boundary=Boundary::from_obj(&obj,&spec,&motion).unwrap().loaded(model.bank()).unwrap();
        (boundary,model,spec)
    }
    #[test]
    fn actual_bem_reaction_changes_mobility_and_closes_mechanical_radiation_power() {
        let (boundary,model,mut spec)=fixture();let w=TAU*200.;
        let a=sample(&boundary,&spec,w).unwrap();
        let loaded=model.solve(200.,69,C64::ONE,Some(&a.impedance)).unwrap();
        let vacuum=model.solve(200.,69,C64::ONE,None).unwrap();
        assert!((loaded.bridge_velocity[0]-vacuum.bridge_velocity[0]).abs()>1e-8*vacuum.bridge_velocity[0].abs());
        assert!(loaded.radiation_w>0.);assert!(loaded.input_w>0.);
        assert!(loaded.power_defect_w.abs()<1e-7*loaded.input_w);
        let p=a.pressure(&loaded,w).unwrap();assert!(p[0].abs()>0.);
        assert!((p[0]+p[1]).abs()<0.02*p[0].abs());
        // Fluid density changes the actual load, not an output-only gain.
        spec.medium.density*=2.;let b=sample(&boundary,&spec,w).unwrap();
        assert!((b.impedance[0]-a.impedance[0].scale(2.)).abs()<1e-10*a.impedance[0].abs());
        let twice=model.solve(200.,69,C64::ONE,Some(&b.impedance)).unwrap();
        assert!((twice.bridge_velocity[0]-loaded.bridge_velocity[0]).abs()>1e-8*loaded.bridge_velocity[0].abs());
        assert!(model.bank().q.iter().chain(&model.bank().v).all(|v|*v==0.));
    }
    #[test]
    fn modal_pressure_projection_preserves_panel_work_and_receiver_linearity() {
        let (boundary,model,spec)=fixture();let w=TAU*100.;let load=sample(&boundary,&spec,w).unwrap();
        let field:Vec<_>=boundary.weights[0].iter().map(|v|C64::new(*v,0.)).collect();
        let direct=helmholtz::solve_radiation_batch(&boundary.surface,w/spec.medium.sound_speed,
            spec.medium,&[field.as_slice()],Formulation::PlainCbie).unwrap();
        let flux=direct[0].pressure.iter().zip(&field).zip(boundary.surface.areas())
            .fold(C64::ZERO,|sum,((p,v),area)|sum+v.conj()* *p *C64::new(*area,0.));
        assert!((load.impedance[0]-flux).abs()<1e-13*(1.+flux.abs()));
        let a=model.solve(100.,69,C64::ONE,Some(&load.impedance)).unwrap();
        let b=model.solve(100.,69,C64::new(2.,0.),Some(&load.impedance)).unwrap();
        let pa=load.pressure(&a,w).unwrap();let pb=load.pressure(&b,w).unwrap();
        for (a,b) in pa.iter().zip(pb) {assert!((a.scale(2.)-b).abs()<1e-12*(1.+a.abs()));}
        assert!((b.input_w-4.*a.input_w).abs()<1e-12*(1.+a.input_w));
    }
    #[test]
    fn harmonic_source_sweep_preserves_impedance_receiver_rows_and_admission() {
        let (boundary,_,spec)=fixture();let omega=spec.omega();
        let prepared=LoadingSweep::new(&boundary,&spec,&omega).unwrap();
        assert_eq!(prepared.delays_s(),ReceiverSet::for_spec(&boundary,&spec).unwrap().delays_s());
        for w in omega {
            let cached=prepared.sample(w).unwrap();let direct=sample(&boundary,&spec,w).unwrap();
            assert_eq!(cached.impedance,direct.impedance);
            assert_eq!(cached.receiver_transfer,direct.receiver_transfer);
            assert_eq!(cached.minimum_ppw,direct.minimum_ppw);
            assert_eq!(cached.condition_lower_bound,direct.condition_lower_bound);
        }
        for w in [0.,f64::NAN,TAU*500.] {assert!(prepared.sample(w).is_err());}
        assert!(LoadingSweep::new(&boundary,&spec,&[]).is_err());
    }
    #[test]
    fn harmonic_array_preserves_source_impedance_and_individual_pressure() {
        let (boundary,model,mut spec)=fixture();let w=TAU*100.;
        let original=sample(&boundary,&spec,w).unwrap();
        spec.receivers=(0..MAX_HARMONIC_RECEIVERS).map(|i|
            [0.02+0.001*i as f64,0.04,1.]).collect();
        spec.near_field_receivers=true;
        let points=spec.receivers.clone();let array=sample(&boundary,&spec,w).unwrap();
        assert_eq!(array.receiver_transfer.len(),MAX_HARMONIC_RECEIVERS);
        assert_eq!(array.impedance,original.impedance);
        let response=model.solve(100.,69,C64::ONE,Some(&array.impedance)).unwrap();
        let pressures=array.pressure(&response,w).unwrap();
        for (i,point) in points.into_iter().enumerate() {
            spec.receivers=vec![point];let single=sample(&boundary,&spec,w).unwrap();
            assert_eq!(array.impedance,single.impedance);
            assert_eq!(array.receiver_transfer[i],single.receiver_transfer[0]);
            assert_eq!(pressures[i],single.pressure(&response,w).unwrap()[0]);
        }
        spec.receivers=vec![[0.05,0.05,1.];MAX_HARMONIC_RECEIVERS+1];
        assert!(sample(&boundary,&spec,w).is_err());
        spec.receivers.clear();assert!(sample(&boundary,&spec,w).is_err());
    }
    #[test]
    fn failed_geometry_frequency_or_drive_has_no_partial_response() {
        let (boundary,model,mut spec)=fixture();
        assert!(sweep(&boundary,&model,&spec,60).is_err());
        for w in [f64::NAN,0.,TAU*500.] {assert!(sample(&boundary,&spec,w).is_err());}
        spec.receivers[0]=boundary.center;
        assert!(sample(&boundary,&spec,TAU*100.).is_err());
        assert!(model.bank().q.iter().chain(&model.bank().v).all(|v|*v==0.));
    }
}

#[cfg(test)]
#[path="radiation_policy_consumer_tests.rs"]
mod policy_tests;
