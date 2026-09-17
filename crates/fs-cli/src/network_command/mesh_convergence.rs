//! Observed mesh ladders for the existing steady coupled cooling model.
//! Optional goal-recovery marking uses the actual coupled adjoint to prioritize
//! local refinement. Scores only choose cells; measured objective differences
//! and a global confirmation govern stopping, never a claimed error bound.
use super::*;
use fs_mesh::{TetRefinement,TetRefinementError,TetRefinementLimits};
mod mark;
mod split;
mod physics;

#[derive(Debug)]
pub(super) struct Study {
    input: J,
    max_refinements: usize,
    consecutive: usize,
    tolerance_k: f64,
    limits: TetRefinementLimits,
    marking: Option<f64>,
}

fn exhausted(message: impl Into<String>)->Failure {
    Failure{code:"cooling-network-mesh-budget",message:message.into()}
}
fn refinement_error(error:TetRefinementError)->Failure {
    match error {
        TetRefinementError::OutputLimit=>exhausted("the next complete mesh exceeds the declared vertex/tetrahedron limit; no convergence result published"),
        TetRefinementError::Cancelled=>Failure{code:"cooling-network-cancelled",message:error.to_string()},
        _=>producer(error),
    }
}

impl Study {
    pub(super) fn parse(value:&J,root:&J)->Result<Self> {
        object(value,&["max_refinements","consecutive_passes","temperature_tolerance_k","max_vertices","max_tetrahedra","strategy","marking_fraction"],"mesh_convergence")?;
        if ["transient","design","fan_speed_design"].iter().any(|key|root.get(key).is_some())
            || get(get(root,"objective")?,"gradient")?!=&J::Bool(false) {
            return Err(bad("mesh_convergence requires a steady request with gradient=false and no design search"));
        }
        let marking=match value.get("strategy") {
            None | Some(J::Str(_)) if value.get("strategy").is_none()
                || value.str_field("strategy")==Some("uniform") => {
                if value.get("marking_fraction").is_some() { return Err(bad("marking_fraction requires strategy=goal-recovery")); }
                None
            }
            Some(J::Str(strategy)) if strategy=="goal-recovery" => {
                let fraction=positive(get(value,"marking_fraction")?,"mesh marking_fraction")?;
                if fraction>1.0 { return Err(bad("mesh marking_fraction must be in (0,1]")); }
                Some(fraction)
            }
            _=>return Err(bad("mesh strategy must be uniform or goal-recovery")),
        };
        let cap=if marking.is_some(){32}else{6};
        let max_refinements=count(get(value,"max_refinements")?,"mesh max_refinements",cap)?;
        let consecutive=count(get(value,"consecutive_passes")?,"mesh consecutive_passes",cap)?;
        if consecutive<2 || consecutive>max_refinements {
            return Err(bad("mesh convergence requires 2 <= consecutive_passes <= max_refinements"));
        }
        let tolerance_k=positive(get(value,"temperature_tolerance_k")?,"mesh temperature_tolerance_k")?;
        let limits=TetRefinementLimits {
            max_vertices:count(get(value,"max_vertices")?,"mesh max_vertices",20_000)?,
            max_tetrahedra:count(get(value,"max_tetrahedra")?,"mesh max_tetrahedra",100_000)?,
        };
        let solid=get(root,"solid")?;
        if array(get(solid,"vertices_m")?,"vertices_m",20_000)?.len()>limits.max_vertices
            || array(get(solid,"tetrahedra")?,"tetrahedra",100_000)?.len()>limits.max_tetrahedra {
            return Err(exhausted("base mesh already exceeds mesh_convergence output limits"));
        }
        let mut input=root.clone();
        members(&mut input)?.retain(|(key,_)|key!="mesh_convergence");
        Ok(Self{input,max_refinements,consecutive,tolerance_k,limits,marking})
    }

    pub(super) fn solve(&self,gate:&CancelGate)->Result<String> {
        let mut input=self.input.clone();
        let mut request=Request::parse(&encode(&input)?)?;
        let mut previous:Option<f64>=None;
        let mut original_power:Option<f64>=None;
        let mut streak=0;
        let mut history=Vec::new();
        let mut last_change=None;
        let mut total_solves=0_usize;
        let mut total_adjoint_sweeps=0_usize;
        let mut arrived_by="base";
        for level in 0..=self.max_refinements {
            let physics::Solved { output, objective_k:value, source_w:power,
                solid_solves:solves, marking, adjoint_sweeps } =
                with_context(&request,gate,|cx|physics::solve(&request,cx,self.marking))?;
            let base_power=*original_power.get_or_insert(power);
            if !(power-base_power).is_finite() || (power-base_power).abs()>request.limits.heat {
                return Err(producer("refinement changed integrated source power beyond the original watt tolerance"));
            }
            let change=previous.map(|prior|(value-prior).abs());
            if change.is_some_and(|d|!d.is_finite()) {return Err(producer("nonfinite mesh objective change"));}
            streak=if change.is_some_and(|d|d<=self.tolerance_k){streak+1}else{0};
            total_solves=total_solves.checked_add(solves).ok_or_else(||exhausted("mesh study work overflow"))?;
            total_adjoint_sweeps=total_adjoint_sweeps.checked_add(adjoint_sweeps).ok_or_else(||exhausted("mesh study adjoint work overflow"))?;
            let marker=match &marking {
                None=>"null".into(),
                Some(m)=>format!("{{\"method\":\"goal-weighted-gradient-recovery\",\"marked_cells\":{},\"normalized_score_sum\":{},\"captured_fraction\":{},\"adjoint_sweeps\":{},\"score_is_error_bound\":false}}",
                    m.cells.len(),num(m.total)?,num(m.captured_fraction)?,adjoint_sweeps),
            };
            history.push(format!("{{\"level\":{level},\"vertices\":{},\"tetrahedra\":{},\"objective_k\":{},\"successive_change_k\":{},\"source_w\":{},\"source_change_w\":{},\"solid_solves\":{solves},\"consecutive_passes\":{streak},\"arrived_by\":{},\"marking\":{marker}}}",
                request.mesh.vertex_count(),request.mesh.element_count(),num(value)?,optional(change)?,num(power)?,num(power-base_power)?,quote(arrived_by)));
            // Local agreement can miss an unmarked region. Before accepting an
            // adaptive result require one complete uniform refinement from the
            // locally agreed mesh. This is still an observed check, not a bound.
            let confirmed=self.marking.is_none() || arrived_by=="uniform";
            if streak>=self.consecutive && confirmed {
                let prefix=output.strip_suffix("}\n").ok_or_else(||bad("internal mesh-study result framing"))?;
                let resolved=encode(&input)?;
                let method=if self.marking.is_some(){"goal-recovery-edge-bisection"}else{"uniform-red-tet-refinement"};
                return Ok(format!("{prefix},\"mesh_convergence\":{{\"status\":\"successive-mesh-tolerance-met\",\"method\":{},\"meshes_solved\":{},\"refinements\":{level},\"temperature_tolerance_k\":{},\"required_consecutive_passes\":{},\"achieved_change_k\":{},\"total_solid_solves\":{total_solves},\"total_adjoint_sweeps\":{total_adjoint_sweeps},\"global_confirmation\":{},\"history\":[{}],\"resolved_request\":{},\"scope\":\"observed same-model successive-mesh agreement only; goal-recovery scores prioritize cells and are not a DWR or continuum error bound, maximum-norm certificate, or physical validation; local refinement is conforming but carries no shape-regularity theorem; an adaptive success includes a complete uniform-refinement comparison; base P1 source is prolonged without renormalization, material laws inherit by parent cell, matching contact traces remain separate; radiation patch partition and constitutive law remain unchanged and its total adjoint drives marking when present; solid work includes nested radiation and reconstruction solves; point-set objectives retain original vertices; resolved_request owns the published field's mesh and may be solved independently\"}}}}\n",
                    quote(method),history.len(),num(self.tolerance_k)?,self.consecutive,optional(change)?,
                    if self.marking.is_some(){"true"}else{"null"},history.join(","),resolved.trim_end()));
            }
            last_change=change;
            previous=Some(value);
            if level==self.max_refinements {break;}
            // No recovery signal requests a global probe, never early success.
            // A pending global confirmation also overrides local marking.
            let marks=marking.as_ref().filter(|m|!m.cells.is_empty() && streak<self.consecutive)
                .map(|m|m.cells.as_slice());
            input=with_context(&request,gate,|cx|refine_request_selected(cx,&input,&request,self.limits,marks))?;
            arrived_by=if marks.is_some(){"marked-edge-stars"}else{"uniform"};
            request=Request::parse(&encode(&input)?)?;
        }
        Err(exhausted(format!("mesh refinement budget exhausted after {} solved meshes: last change {:?} K, tolerance {} K, consecutive passes {streak}/{}; any required global confirmation is still part of the budget; no convergence result published",
            history.len(),last_change,self.tolerance_k,self.consecutive)))
    }
}

fn with_context<T>(request:&Request,gate:&CancelGate,run:impl FnOnce(&Cx<'_>)->Result<T>)->Result<T> {
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        let cx=Cx::new(gate,arena,StreamKey{seed:request.seed,kernel_id:717,tile:0,iteration:0},
            Budget::INFINITE,ExecMode::Deterministic);
        poll(&cx)?;run(&cx)
    })
}

/// Uniform compatibility entry point used by the existing transfer tests.
#[cfg(test)]
fn refine_request(cx:&Cx<'_>,root:&J,request:&Request,limits:TetRefinementLimits)->Result<J> {
    refine_request_selected(cx,root,request,limits,None)
}

/// Transfer the DECLARED problem, not just node coordinates. Recreating the
/// old component vertex set on a finer mesh would shrink its source support.
fn refine_request_selected(cx:&Cx<'_>,root:&J,request:&Request,limits:TetRefinementLimits,
    marked:Option<&[usize]>)->Result<J> {
    poll(cx)?;
    let solid=get(root,"solid")?;
    let old_vertices=request.mesh.vertex_count();
    let tets=array(get(solid,"tetrahedra")?,"tetrahedra",100_000)?.iter()
        .map(|t|indices::<4>(t,"tetrahedron",old_vertices)).collect::<Result<Vec<_>>>()?;
    let split=split::Split::build(cx,root,request,&tets,limits,marked)?;
    let mut next=root.clone();
    let solid=member_mut(&mut next,"solid")?;
    replace(solid,"vertices_m",J::Array(split.positions().iter()
        .map(|p|J::Array(p.iter().map(|&v|jnum(v)).collect())).collect()))?;
    replace(solid,"tetrahedra",J::Array(split.tetrahedra().iter().map(|t|jindices(t)).collect()))?;
    if let Some(assignment)=solid.get("element_materials") {
        let rows=array(assignment,"element_materials",tets.len())?;
        if rows.len()!=tets.len(){return Err(bad("material assignment changed during refinement"));}
        let refined=(0..split.tetrahedra().len()).map(|i|rows[split.parent(i)].clone()).collect();
        replace(solid,"element_materials",J::Array(refined))?;
    }
    if let Some(source)=&request.solid_data.nodal_source {
        let old:Vec<f64>=(0..old_vertices).map(|i|source.at(i)).collect();
        let values=split.prolongate(cx,&old)?;
        let fields=members(solid)?;
        fields.retain(|(key,_)|key!="component_power" && key!="nodal_source_w_m3");
        fields.push(("nodal_source_w_m3".into(),J::Array(values.into_iter().map(jnum).collect())));
    }
    let J::Array(surfaces)=member_mut(solid,"surfaces")? else {return Err(bad("surfaces must be an array"));};
    for surface in surfaces {
        poll(cx)?;
        let mut children=Vec::new();
        for face in array(get(surface,"faces")?,"surface faces",200_000)? {
            let face=indices::<3>(face,"surface face",old_vertices)?;
            children.extend(split.face_children(face)?.iter().map(|f|jindices(f)));
        }
        replace(surface,"faces",J::Array(children))?;
    }
    if solid.get("contacts").is_some() {
        let J::Array(contacts)=member_mut(solid,"contacts")? else {return Err(bad("contacts must be an array"));};
        for contact in contacts {
            let mut children=Vec::new();
            for pair in array(get(contact,"face_pairs")?,"contact face pairs",200_000)? {
                poll(cx)?;
                let a=indices::<3>(get(pair,"side_a")?,"side_a",old_vertices)?;
                let b=indices::<3>(get(pair,"side_b")?,"side_b",old_vertices)?;
                for (a,b) in split.contact_children(a,b)? {
                    children.push(J::Object(vec![("side_a".into(),jindices(&a)),("side_b".into(),jindices(&b))]));
                }
            }
            replace(contact,"face_pairs",J::Array(children))?;
        }
    }
    poll(cx)?;
    Ok(next)
}

fn members(value:&mut J)->Result<&mut Vec<(String,J)>> {
    if let J::Object(fields)=value{Ok(fields)}else{Err(bad("expected refinement object"))}
}
fn member_mut<'a>(value:&'a mut J,key:&str)->Result<&'a mut J> {
    members(value)?.iter_mut().find_map(|(name,value)|(name==key).then_some(value))
        .ok_or_else(||bad(format!("missing refinement field {key}")))
}
fn replace(value:&mut J,key:&str,new:J)->Result<()> {*member_mut(value,key)?=new;Ok(())}
fn jnum(value:f64)->J {J::Number{value,raw:value.to_string()}}
fn jindices<const N:usize>(values:&[u32;N])->J {
    J::Array(values.iter().map(|&v|J::Number{value:f64::from(v),raw:v.to_string()}).collect())
}
fn encode(value:&J)->Result<String> {
    fn write(value:&J,out:&mut String)->Result<()> {
        match value {
            J::Null=>out.push_str("null"),J::Bool(v)=>out.push_str(if *v{"true"}else{"false"}),
            J::Number{value,raw}=>{if !value.is_finite(){return Err(bad("nonfinite refined input"));}out.push_str(raw);},
            J::Str(s)=>out.push_str(&quote(s)),
            J::Array(values)=>{
                out.push('[');for(i,v)in values.iter().enumerate(){if i>0{out.push(',');}write(v,out)?;}out.push(']');
            },
            J::Object(fields)=>{
                out.push('{');for(i,(k,v))in fields.iter().enumerate(){if i>0{out.push(',');}out.push_str(&quote(k));out.push(':');write(v,out)?;}out.push('}');
            },
        }
        if out.len() as u64>MAX_INPUT_BYTES{return Err(exhausted("refined request exceeds the existing input byte limit"));}Ok(())
    }
    let mut result=String::new();write(value,&mut result)?;Ok(result)
}

#[cfg(test)]
mod tests;
