//! An observed uniform h-ladder for the existing steady cooling model.
//! Each rung executes the real coupled producer. Agreement between meshes is
//! an empirical stopping condition, not a continuum bound or DWR adaptation.
use super::*;
use fs_mesh::{TetRefinement,TetRefinementError,TetRefinementLimits};

#[derive(Debug)]
pub(super) struct Study {
    input: J,
    max_refinements: usize,
    consecutive: usize,
    tolerance_k: f64,
    limits: TetRefinementLimits,
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
        object(value,&["max_refinements","consecutive_passes","temperature_tolerance_k","max_vertices","max_tetrahedra"],"mesh_convergence")?;
        if ["transient","design","fan_speed_design"].iter().any(|key|root.get(key).is_some())
            || get(get(root,"objective")?,"gradient")?!=&J::Bool(false) {
            return Err(bad("mesh_convergence requires a steady request with gradient=false and no design search"));
        }
        let max_refinements=count(get(value,"max_refinements")?,"mesh max_refinements",6)?;
        let consecutive=count(get(value,"consecutive_passes")?,"mesh consecutive_passes",6)?;
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
        Ok(Self{input,max_refinements,consecutive,tolerance_k,limits})
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
        for level in 0..=self.max_refinements {
            let (output,value,power,solves)=with_context(&request,gate,|cx| {
                let flow=request.flow(cx)?;
                let coefficients=request.surfaces.iter().map(|s|(s.name.clone(),s.h)).collect();
                let evaluated=request.evaluate(cx,&flow,&coefficients,false)?;
                let result=render(&request,&flow,&evaluated)?;
                let result=match &request.fan {
                    Some(fan)=>fan.attach(result,&flow,fan.speed_ratio)?,None=>result,
                };
                poll(cx)?;
                Ok((result,evaluated.objective,evaluated.source_total_w,evaluated.coupled.iterations))
            })?;
            let base_power=*original_power.get_or_insert(power);
            if !(power-base_power).is_finite() || (power-base_power).abs()>request.limits.heat {
                return Err(producer("refinement changed integrated source power beyond the original watt tolerance"));
            }
            let change=previous.map(|prior|(value-prior).abs());
            if change.is_some_and(|d|!d.is_finite()) {return Err(producer("nonfinite mesh objective change"));}
            streak=if change.is_some_and(|d|d<=self.tolerance_k){streak+1}else{0};
            total_solves=total_solves.checked_add(solves).ok_or_else(||exhausted("mesh study work overflow"))?;
            history.push(format!("{{\"level\":{level},\"vertices\":{},\"tetrahedra\":{},\"objective_k\":{},\"successive_change_k\":{},\"source_w\":{},\"source_change_w\":{},\"solid_solves\":{solves},\"consecutive_passes\":{streak}}}",
                request.mesh.vertex_count(),request.mesh.element_count(),num(value)?,optional(change)?,num(power)?,num(power-base_power)?));
            if streak>=self.consecutive {
                let prefix=output.strip_suffix("}\n").ok_or_else(||bad("internal mesh-study result framing"))?;
                let resolved=encode(&input)?;
                return Ok(format!("{prefix},\"mesh_convergence\":{{\"status\":\"successive-mesh-tolerance-met\",\"method\":\"uniform-red-tet-refinement\",\"meshes_solved\":{},\"refinements\":{level},\"temperature_tolerance_k\":{},\"required_consecutive_passes\":{},\"achieved_change_k\":{},\"total_solid_solves\":{total_solves},\"history\":[{}],\"resolved_request\":{},\"scope\":\"observed same-model successive-mesh agreement only; not a continuum error bound, maximum-norm certificate, goal-oriented adaptation or physical validation; base P1 source is prolonged without renormalization, material laws inherit by parent cell, matching contact traces remain separate; point-set objectives retain original vertices; resolved_request owns the published field's mesh and may be solved independently\"}}}}\n",
                    history.len(),num(self.tolerance_k)?,self.consecutive,optional(change)?,history.join(","),resolved.trim_end()));
            }
            last_change=change;
            previous=Some(value);
            if level==self.max_refinements {break;}
            input=with_context(&request,gate,|cx|refine_request(cx,&input,&request,self.limits))?;
            request=Request::parse(&encode(&input)?)?;
        }
        Err(exhausted(format!("mesh refinement budget exhausted after {} solved meshes: last change {:?} K, tolerance {} K, consecutive passes {streak}/{}; no convergence result published",
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

/// Transfer the DECLARED problem, not just node coordinates. In particular,
/// recreating a component's original vertex set on a finer mesh would shrink
/// its support and silently turn a convergence study into a changed-load study.
fn refine_request(cx:&Cx<'_>,root:&J,request:&Request,limits:TetRefinementLimits)->Result<J> {
    poll(cx)?;
    let solid=get(root,"solid")?;
    let old_vertices=request.mesh.vertex_count();
    let tets=array(get(solid,"tetrahedra")?,"tetrahedra",100_000)?.iter()
        .map(|t|indices::<4>(t,"tetrahedron",old_vertices)).collect::<Result<Vec<_>>>()?;
    let split=TetRefinement::build(cx,request.mesh.positions(),&tets,limits).map_err(refinement_error)?;
    let mut next=root.clone();
    let solid=member_mut(&mut next,"solid")?;
    replace(solid,"vertices_m",J::Array(split.positions().iter()
        .map(|p|J::Array(p.iter().map(|&v|jnum(v)).collect())).collect()))?;
    replace(solid,"tetrahedra",J::Array(split.tetrahedra().iter().map(|t|jindices(t)).collect()))?;
    if let Some(assignment)=solid.get("element_materials") {
        let rows=array(assignment,"element_materials",tets.len())?;
        if rows.len()!=tets.len(){return Err(bad("material assignment changed during refinement"));}
        let refined=rows.iter().flat_map(|row|std::iter::repeat_n(row.clone(),8)).collect();
        replace(solid,"element_materials",J::Array(refined))?;
    }
    if let Some(source)=&request.solid_data.nodal_source {
        let old:Vec<f64>=(0..old_vertices).map(|i|source.at(i)).collect();
        let values=split.prolongate(cx,&old).map_err(refinement_error)?;
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
            children.extend(split.face_children(face).map_err(refinement_error)?.iter().map(|f|jindices(f)));
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
                let mut matched=[0;3];
                for (i,&vertex) in a.iter().enumerate() {
                    matched[i]=b.iter().copied().find(|&other|
                        request.mesh.positions()[vertex as usize]==request.mesh.positions()[other as usize])
                        .ok_or_else(||bad("mesh refinement requires exactly coincident matching contact vertices"))?;
                }
                let aa=split.face_children(a).map_err(refinement_error)?;
                let bb=split.face_children(matched).map_err(refinement_error)?;
                for (a,b) in aa.iter().zip(bb.iter()) {
                    children.push(J::Object(vec![("side_a".into(),jindices(a)),("side_b".into(),jindices(b))]));
                }
            }
            replace(contact,"face_pairs",J::Array(children))?;
        }
    }
    // Refined source/material/contact inputs go through the original parser.
    // Existing max_vertices IDs intentionally remain fixed observation points.
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
