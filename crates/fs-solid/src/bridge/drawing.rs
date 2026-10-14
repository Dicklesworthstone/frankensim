//! Declared drawing tension/slack for the history-dependent deposited beam.
//!
//! Positive drawing force is tension. Slack delegates exactly to the ordinary
//! zero-axial-resultant path. Tension enters the same implicit material section
//! update and the small-deflection beam-column equilibrium, so it changes both
//! stress history and transverse sag. Compression is refused rather than
//! silently entering an unverified buckling lane.

use std::collections::BTreeMap;

use super::{Bridge, Config as BeamConfig, Error, Report as BeamReport, Support,
    add_energy, finite, poll, positive};
use crate::deposition::{self, AxialControl, BendingControl, Control};
use fs_material::thermovisco::Energy;

#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Drawing { Slack, TensionN(f64) }

impl Drawing {
    fn force_n(self) -> Result<f64, Error> {
        match self {
            Self::Slack => Ok(0.0),
            Self::TensionN(value) => { positive(value)?; Ok(value) }
        }
    }
}

#[derive(Debug, Clone)]
pub struct Report {
    pub drawing: Drawing,
    pub drawing_force_n: f64,
    pub vertical_reactions_n: [f64; 2],
    pub clamp_reaction_nm: f64,
    pub moments_nm: Vec<f64>,
    pub maximum_sag_m: f64,
    pub maximum_absolute_deflection_m: f64,
    pub maximum_absolute_rotation_rad: f64,
    pub material_point_updates: usize,
    pub axial_extension_m: f64,
    pub energy_j: Energy,
    pub transverse_load_work_j: f64,
    pub axial_material_work_j: f64,
    pub geometric_bending_work_j: f64,
    pub minimum_cell_condensation: f64,
}

#[derive(Clone, Copy, Debug, Default)]
struct State { w: f64, theta: f64, moment: f64, shear: f64 }

#[derive(Clone, Copy, Debug)]
struct ProbeCell { h: f64, q: f64, free_kappa: f64, compliance: f64 }

#[derive(Clone, Debug)]
struct Propagation {
    end: State,
    starts: Vec<State>,
    average_moments: Vec<f64>,
    curvatures: Vec<f64>,
    minimum_denominator: f64,
}

impl Bridge<'_> {
    pub fn advance_drawn_controlled(&self, target_s: f64,
        temperatures_k: &BTreeMap<(u64,u64), f64>, gravity_m_s2: f64,
        additional_downward_n_m: &[f64], drawing: Drawing, config: BeamConfig,
        cancelled: &mut impl FnMut() -> bool) -> Result<(Self, Report), Error> {
        poll(cancelled)?;
        let tension=drawing.force_n()?;
        if tension==0.0 {
            let (next,beam)=self.advance_controlled(target_s,temperatures_k,
                gravity_m_s2,additional_downward_n_m,config,cancelled)?;
            return Ok((next,slack_report(drawing,beam)));
        }
        positive(finite(target_s-self.time_s)?)?;
        if !gravity_m_s2.is_finite() || gravity_m_s2<0.0 {
            return Err(Error::Invalid("gravity must be finite and nonnegative"));
        }
        positive(config.maximum_deflection_m)?;
        let count=self.cells.len();
        if additional_downward_n_m.len()!=count {
            return Err(Error::Invalid("one line load per axial cell required"));
        }
        let layers=self.cells[0].section.layers().len();
        if layers==0 || temperatures_k.len()!=count.checked_mul(layers).ok_or(Error::Budget)? {
            return Err(Error::Invalid("born span and complete endpoint temperatures required"));
        }

        let mut probes=Vec::with_capacity(count);
        for (cell,&added) in self.cells.iter().zip(additional_downward_n_m) {
            poll(cancelled)?;
            let q=finite(gravity_m_s2*cell.mass_per_length_kg_m+finite(added)?)?;
            let mut local=BTreeMap::new();
            for &layer in cell.section.layers().keys() {
                local.insert(layer,*temperatures_k.get(&(cell.id,layer))
                    .ok_or(Error::Invalid("missing section temperature"))?);
            }
            let response=cell.section.bending_response_at_force(target_s,&local,tension,
                || if cancelled(){Err(deposition::Error::Cancelled)}else{Ok(())})
                .map_err(|error|if error==deposition::Error::Cancelled{Error::Cancelled}
                    else{Error::Section{cell:cell.id,error}})?;
            probes.push(ProbeCell{h:cell.length_m(),q,
                free_kappa:response.zero_moment_curvature_per_m,
                compliance:response.compliance_per_n_m2});
        }

        let initial=solve_initial(self.support,&probes,tension)?;
        let predicted=propagate(initial,&probes,tension,true)?;
        let mut next=self.clone();
        let mut energy=Energy::default();
        let mut updates=0usize;
        let mut material_work=0.0;
        let mut axial_work=0.0;
        let mut section_config=config.section;
        section_config.energy_tolerance_j/=count as f64;
        for (index,cell) in next.cells.iter_mut().enumerate() {
            poll(cancelled)?;
            let old=&self.cells[index];
            let mut local=BTreeMap::new();
            for &layer in cell.section.layers().keys() {
                local.insert(layer,*temperatures_k.get(&(cell.id,layer))
                    .ok_or(Error::Invalid("missing section temperature"))?);
            }
            let moment=predicted.average_moments[index];
            let control=Control{axial:AxialControl::ForceN(tension),
                bending:BendingControl::MomentNm(moment)};
            let (section,report)=cell.section.advance_to(target_s,&local,control,section_config,
                ||if cancelled(){Err(deposition::Error::Cancelled)}else{Ok(())})
                .map_err(|error|if error==deposition::Error::Cancelled{Error::Cancelled}
                    else{Error::Section{cell:cell.id,error}})?;
            let h=cell.length_m();
            let axial=finite(tension*h*(section.axial_strain()-old.section.axial_strain()))?;
            let bending=finite(moment*h*(section.curvature_per_m()-old.section.curvature_per_m()))?;
            axial_work=finite(axial_work+axial)?;
            material_work=finite(material_work+axial+bending)?;
            let curvature_tolerance=1e-10_f64.max(
                config.section.moment_tolerance_nm*probes[index].compliance*16.0);
            if (section.curvature_per_m()-predicted.curvatures[index]).abs()>curvature_tolerance {
                return Err(Error::Invalid("drawing prediction disagrees with accepted section curvature"));
            }
            cell.section=section;
            add_energy(&mut energy,report.energy_j)?;
            updates=updates.checked_add(report.material_point_updates).ok_or(Error::Budget)?;
        }
        if finite(material_work-energy.mechanical_work)?.abs()>config.section.energy_tolerance_j {
            return Err(Error::Energy(material_work-energy.mechanical_work));
        }

        let mut max_abs=0.0_f64;
        let mut sag=0.0_f64;
        let mut rotation=0.0_f64;
        let mut u=0.0;
        let mut transverse_work=0.0;
        for (index,cell) in next.cells.iter_mut().enumerate() {
            poll(cancelled)?;
            let start=predicted.starts[index];
            cell.w_m=start.w; cell.slope=start.theta; cell.u_m=u;
            let h=cell.length_m();
            let end=cell.height(h);
            for value in [start.w,end] {
                max_abs=max_abs.max(value.abs()); sag=sag.max(-value);
            }
            let k=cell.section.curvature_per_m();
            if k!=0.0 {
                let s=-start.theta/k;
                if s>0.0 && s<h {
                    let value=cell.height(s);
                    max_abs=max_abs.max(value.abs()); sag=sag.max(-value);
                }
            }
            rotation=rotation.max(start.theta.abs()).max((start.theta+k*h).abs());
            u=finite(u+cell.section.axial_strain()*h)?;
            transverse_work=finite(transverse_work-probes[index].q
                *(cell.integral_height()-self.cells[index].integral_height()))?;
        }
        if max_abs>config.maximum_deflection_m || rotation>config.section.max_rotation_rad {
            return Err(Error::Regime("drawn beam leaves declared deflection/rotation range"));
        }
        next.time_s=target_s;
        let storage=finite(next.stored_energy_j()-self.stored_energy_j())?;
        for defect in [finite(storage-energy.stored_change)?,energy.residual] {
            if defect.abs()>config.section.energy_tolerance_j { return Err(Error::Energy(defect)); }
        }

        let q0=finite(initial.shear-tension*initial.theta)?;
        let ql=finite(predicted.end.shear-tension*predicted.end.theta)?;
        let (vertical,clamp)=match self.support {
            Support::SimplySupported=>([q0,-ql],0.0),
            Support::CantileverLeft=>([q0,0.0],-initial.moment),
        };
        let total_downward=probes.iter().map(|c|c.q*c.h).sum::<f64>();
        if finite(vertical[0]+vertical[1]-total_downward)?.abs()>config.section.force_tolerance_n {
            return Err(Error::Invalid("drawing beam vertical reaction imbalance"));
        }
        poll(cancelled)?;
        Ok((next,Report{drawing,drawing_force_n:tension,vertical_reactions_n:vertical,
            clamp_reaction_nm:clamp,moments_nm:predicted.average_moments,
            maximum_sag_m:sag,maximum_absolute_deflection_m:max_abs,
            maximum_absolute_rotation_rad:rotation,material_point_updates:updates,
            axial_extension_m:u,energy_j:energy,transverse_load_work_j:transverse_work,
            axial_material_work_j:axial_work,
            geometric_bending_work_j:finite(energy.mechanical_work-axial_work-transverse_work)?,
            minimum_cell_condensation:predicted.minimum_denominator}))
    }

    pub fn advance_drawn_to(&self,cx:&fs_exec::Cx<'_>,target_s:f64,
        temperatures_k:&BTreeMap<(u64,u64),f64>,gravity_m_s2:f64,
        additional_downward_n_m:&[f64],drawing:Drawing,config:BeamConfig)
        -> Result<(Self,Report),Error> {
        self.advance_drawn_controlled(target_s,temperatures_k,gravity_m_s2,
            additional_downward_n_m,drawing,config,&mut||cx.checkpoint().is_err())
    }
}

fn slack_report(drawing:Drawing,beam:BeamReport)->Report {
    Report{drawing,drawing_force_n:0.0,vertical_reactions_n:beam.vertical_reactions_n,
        clamp_reaction_nm:beam.clamp_reaction_nm,moments_nm:beam.moments_nm,
        maximum_sag_m:beam.maximum_sag_m,
        maximum_absolute_deflection_m:beam.maximum_absolute_deflection_m,
        maximum_absolute_rotation_rad:beam.maximum_absolute_rotation_rad,
        material_point_updates:beam.material_point_updates,axial_extension_m:0.0,
        energy_j:beam.energy_j,transverse_load_work_j:beam.load_work_j,
        axial_material_work_j:0.0,geometric_bending_work_j:0.0,
        minimum_cell_condensation:1.0}
}

fn solve_initial(support:Support,cells:&[ProbeCell],tension:f64)->Result<State,Error> {
    let zero=State::default();
    let base=propagate(zero,cells,tension,false)?.end;
    let (a,rhs,kind)=match support {
        Support::SimplySupported=>{
            let pa=propagate(State{theta:1.0,..zero},cells,tension,false)?.end;
            let pb=propagate(State{shear:1.0,..zero},cells,tension,false)?.end;
            ([[pa.w-base.w,pb.w-base.w],[pa.moment-base.moment,pb.moment-base.moment]],
             [-base.w,-base.moment],0_u8)
        }
        Support::CantileverLeft=>{
            let pa=propagate(State{moment:1.0,..zero},cells,tension,false)?.end;
            let pb=propagate(State{shear:1.0,..zero},cells,tension,false)?.end;
            let rb=base.shear-tension*base.theta;
            ([[pa.moment-base.moment,pb.moment-base.moment],
              [(pa.shear-tension*pa.theta)-rb,(pb.shear-tension*pb.theta)-rb]],
             [-base.moment,-rb],1_u8)
        }
    };
    let det=finite(a[0][0]*a[1][1]-a[0][1]*a[1][0])?;
    let scale=a.iter().flatten().fold(0.0_f64,|m,v|m.max(v.abs())).max(1.0);
    if det.abs()<=64.0*f64::EPSILON*scale*scale {
        return Err(Error::Invalid("drawing boundary solve is singular; refine or change support/tension"));
    }
    let x0=finite((rhs[0]*a[1][1]-a[0][1]*rhs[1])/det)?;
    let x1=finite((a[0][0]*rhs[1]-rhs[0]*a[1][0])/det)?;
    Ok(if kind==0{State{theta:x0,shear:x1,..zero}}
       else{State{moment:x0,shear:x1,..zero}})
}

fn propagate(initial:State,cells:&[ProbeCell],tension:f64,record:bool)
    -> Result<Propagation,Error> {
    let mut state=initial;
    let mut starts=if record{Vec::with_capacity(cells.len())}else{Vec::new()};
    let mut moments=if record{Vec::with_capacity(cells.len())}else{Vec::new()};
    let mut curvatures=if record{Vec::with_capacity(cells.len())}else{Vec::new()};
    let mut minimum=1.0_f64;
    for cell in cells {
        if record{starts.push(state);}
        let tc=finite(tension*cell.compliance)?;
        let denominator=finite(1.0-tc*cell.h*cell.h/6.0)?;
        if denominator<=64.0*f64::EPSILON {
            return Err(Error::Regime("drawing tension requires finer axial cells"));
        }
        minimum=minimum.min(denominator);
        let acceleration=finite((tension*cell.free_kappa
            +tc*(state.moment+0.5*state.shear*cell.h)-cell.q)/denominator)?;
        let average_moment=finite(state.moment+0.5*state.shear*cell.h
            +acceleration*cell.h*cell.h/6.0)?;
        let curvature=finite(cell.free_kappa+cell.compliance*average_moment)?;
        if record{moments.push(average_moment);curvatures.push(curvature);}
        state.w=finite(state.w+state.theta*cell.h+0.5*curvature*cell.h*cell.h)?;
        state.theta=finite(state.theta+curvature*cell.h)?;
        state.moment=finite(state.moment+state.shear*cell.h
            +0.5*acceleration*cell.h*cell.h)?;
        state.shear=finite(state.shear+acceleration*cell.h)?;
    }
    Ok(Propagation{end:state,starts,average_moments:moments,curvatures,
        minimum_denominator:minimum})
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::bridge::tests::{birth,config,material,temperatures};

    fn axial_force(cell:&super::super::Cell<'_>)->f64 {
        cell.section().layers().values().map(|layer|
            0.5*layer.geometry.area_m2
                *layer.material_points().iter().map(|p|p.stress_pa()).sum::<f64>()
        ).sum()
    }
    fn near(a:f64,b:f64,t:f64){assert!((a-b).abs()<t,"{a:e} != {b:e}");}

    #[test]
    fn slack_is_exactly_the_existing_beam_path() {
        let law=material(1e8,vec![(2e8,0.5)]);
        let old=birth(8,Support::SimplySupported,&law);
        let t=temperatures(&old,340.0);
        let plain=old.advance_controlled(0.25,&t,9.81,&[0.0;8],config(),&mut||false).unwrap();
        let drawn=old.advance_drawn_controlled(0.25,&t,9.81,&[0.0;8],
            Drawing::Slack,config(),&mut||false).unwrap();
        assert_eq!(plain.0,drawn.0);
        assert_eq!(plain.1.moments_nm,drawn.1.moments_nm);
        assert_eq!(plain.1.maximum_sag_m.to_bits(),drawn.1.maximum_sag_m.to_bits());
    }

    #[test]
    fn elastic_tension_reduces_sag_and_converges_to_beam_column_solution() {
        let law=material(1e8,vec![]);
        let tension=0.01;
        let ei=1e8*0.001*0.001_f64.powi(3)/12.0;
        let q=1200.0*1e-6*9.81;
        let l=0.02;
        let k=fs_math::det::sqrt(tension/ei);
        let exact=q*l*l/(8.0*tension)-q*ei/(tension*tension)
            *(1.0-1.0/fs_math::det::cosh(0.5*k*l));
        let mut errors=Vec::new();
        for n in [4_u32,8,16] {
            let old=birth(n,Support::SimplySupported,&law);
            let (next,r)=old.advance_drawn_controlled(1.0,&temperatures(&old,300.0),
                9.81,&vec![0.0;n as usize],Drawing::TensionN(tension),config(),&mut||false).unwrap();
            let sag=-next.deformed_position_m(l/2.0,0.0).unwrap()[1];
            errors.push((sag-exact).abs());
            near(r.vertical_reactions_n.iter().sum(),q*l,1e-10);
            for cell in next.cells(){near(axial_force(cell),tension,1e-8);}
        }
        assert!(errors[0]/errors[1]>3.9&&errors[1]/errors[2]>3.9);
        let old=birth(16,Support::SimplySupported,&law);
        let slack=old.advance_drawn_controlled(1.0,&temperatures(&old,300.0),
            9.81,&[0.0;16],Drawing::Slack,config(),&mut||false).unwrap().1.maximum_sag_m;
        assert!(exact<slack);
    }

    #[test]
    fn thermal_history_and_tension_share_the_same_material_state() {
        let law=material(2e7,vec![(8e7,0.3)]);
        let mut state=birth(8,Support::CantileverLeft,&law);
        for i in 1..=20 {
            let temperature=500.0-5.0*f64::from(i);
            state=state.advance_drawn_controlled(f64::from(i)/20.0,&temperatures(&state,temperature),
                9.81,&[0.0;8],Drawing::TensionN(0.002),config(),&mut||false).unwrap().0;
        }
        assert!(state.cells()[0].section().layers()[&7].material_points()[0]
            .branch_strains()[0]!=0.0);
        for cell in state.cells(){near(axial_force(cell),0.002,1e-7);}
    }

    #[test]
    fn compression_coarse_tension_and_cancellation_refuse_atomically() {
        let law=material(1e8,vec![]);
        let old=birth(4,Support::SimplySupported,&law);
        let before=old.clone();
        let t=temperatures(&old,300.0);
        assert!(old.advance_drawn_controlled(1.0,&t,9.81,&[0.0;4],
            Drawing::TensionN(-1.0),config(),&mut||false).is_err());
        assert!(old.advance_drawn_controlled(1.0,&t,9.81,&[0.0;4],
            Drawing::TensionN(1e9),config(),&mut||false).is_err());
        let mut calls=0;
        assert!(old.advance_drawn_controlled(1.0,&t,9.81,&[0.0;4],
            Drawing::TensionN(0.01),config(),&mut||{calls+=1;calls>5}).is_err());
        assert_eq!(old,before);
    }
}
