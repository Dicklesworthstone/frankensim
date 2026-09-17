//! Actual-command regressions, not a standalone replacement contact solver.
#![cfg(unix)]
#[allow(dead_code)]
#[path="../src/json_read.rs"]
mod json;
mod mesh;
use json::JsonValue as J;
use std::io::Write;
use std::process::{Command,Output,Stdio};
use std::path::PathBuf;
use std::sync::atomic::{AtomicU64,Ordering};
const BASE:&str=include_str!(concat!(env!("CARGO_MANIFEST_DIR"),
    "/../../examples/cooling-network/nonmatching-contact-hotspot.json"));
fn num(value:f64)->J {J::Number{value,raw:value.to_string()}}
fn member<'a>(root:&'a mut J,key:&str)->&'a mut J {
    let J::Object(fields)=root else {panic!("object required")};
    &mut fields.iter_mut().find(|(k,_)|k==key).unwrap().1
}
fn put(root:&mut J,key:&str,value:J) {
    let J::Object(fields)=root else {panic!("object required")};
    if let Some((_,v))=fields.iter_mut().find(|(k,_)|k==key){*v=value;}
    else {fields.push((key.into(),value));}
}
fn remove(root:&mut J,key:&str) {
    let J::Object(fields)=root else {panic!("object required")};fields.retain(|(k,_)|k!=key);
}
fn rows(root:&mut J)->&mut Vec<J> {let J::Array(rows)=root else {panic!("array required")};rows}
fn contact(root:&mut J)->&mut J {&mut rows(member(member(root,"solid"),"contacts"))[0]}
fn text(root:&J)->String {
    match root {
        J::Null=>"null".into(),J::Bool(b)=>b.to_string(),J::Number{raw,..}=>raw.clone(),
        J::Str(s)=>format!("\"{}\"",s.replace('\\',"\\\\").replace('"',"\\\"").replace('\n',"\\n")),
        J::Array(a)=>format!("[{}]",a.iter().map(text).collect::<Vec<_>>().join(",")),
        J::Object(a)=>format!("{{{}}}",a.iter().map(|(k,v)|format!("{}:{}",text(&J::Str(k.clone())),text(v))).collect::<Vec<_>>().join(",")),
    }
}
fn output(input:&J)->Output {
    let mut child=Command::new(env!("CARGO_BIN_EXE_frankensim"))
        .args(["--json","cooling-network","/dev/stdin"])
        .stdin(Stdio::piped()).stdout(Stdio::piped()).stderr(Stdio::piped()).spawn().unwrap();
    child.stdin.take().unwrap().write_all(text(input).as_bytes()).unwrap();child.wait_with_output().unwrap()
}
fn success(output:&Output)->J {
    assert!(output.status.success(),"{}",String::from_utf8_lossy(&output.stderr));
    J::parse(std::str::from_utf8(&output.stdout).unwrap()).unwrap()
}
fn run(input:&J)->J {success(&output(input))}
fn n(root:&J,key:&str)->f64 {root.f64_field(key).unwrap()}
fn peak(root:&J)->f64 {n(root.get("objective").unwrap(),"value_k")}
fn near(a:f64,b:f64,tol:f64){assert!((a-b).abs()<tol,"{a} versus {b}");}
fn gradient(root:&J)->f64 {
    n(&root.path(&["contact_sensitivities","rows"]).unwrap().as_array().unwrap()[0],"dobjective_dlog_resistance_k")
}
fn nonlinear(root:&mut J) {
    let materials=rows(member(member(root,"solid"),"materials"));
    for (row,curve) in materials.iter_mut().zip([
        r#"{"temperature_k":[280,400],"conductivity_w_m_k":[8,20]}"#,
        r#"{"temperature_k":[280,400],"conductivity_w_m_k":[1,7]}"#,
    ]) {remove(row,"conductivity_w_m_k");put(row,"conductivity_curve",J::parse(curve).unwrap());}
}
fn radiation(root:&mut J) {
    put(root,"radiation",J::parse(r#"{"max_iterations":128,"temperature_tolerance_k":1e-9,"relaxation":0.5,"surfaces":[{"surface":"first-face","emissivity":0.8,"ambient_temperature_k":300,"source":"test reservoir"},{"surface":"last-face","emissivity":0.8,"ambient_temperature_k":300,"source":"test reservoir"}]}"#).unwrap());
}

#[test]
fn independent_mesh_hotspot_preserves_temperature_nodes_and_matches_the_reference() {
    let input=J::parse(BASE).unwrap();let result=run(&input);
    assert_eq!(result.get("solid_temperatures_k").unwrap().as_array().unwrap().len(),32);
    let flux=&result.get("contacts").unwrap().as_array().unwrap()[0];
    assert_eq!(flux.str_field("discretization"),Some("planar-common-refinement-P1"));
    assert_eq!(flux.get("face_pairs"),Some(&J::Null));
    near(n(flux,"area_m2"),0.01,1e-12);
    // Independent exact-rational overlap moments + direct P1 FEM and analytic
    // air elimination. These tolerances are test comparisons, not error bounds.
    near(peak(&result),301.810057315888,2e-6);
    near(n(flux,"heat_a_to_b_w"),0.359935562787,2e-6);
    near(gradient(&result),0.168944130362,2e-6);
    let external:f64=result.get("walls").unwrap().as_array().unwrap().iter()
        .map(|w|n(w,"outward_heat_w")).sum();
    near(external,1.0,1e-7);near(n(&result,"robin_out_w"),1.0,1e-7);
}

#[test]
fn contact_gradients_match_complete_perturbed_cooling_with_material_and_radiation_feedback() {
    for mode in 0..3 {
        let mut input=J::parse(BASE).unwrap();
        if mode>0{nonlinear(&mut input);}if mode>1{radiation(&mut input);}
        let result=run(&input);let h=1e-4_f64;
        let mut plus=input.clone();let mut minus=input.clone();
        put(contact(&mut plus),"resistance_m2_k_w",num(0.01*h.exp()));
        put(contact(&mut minus),"resistance_m2_k_w",num(0.01*(-h).exp()));
        for trial in [&mut plus,&mut minus]{put(member(trial,"objective"),"gradient",J::Bool(false));}
        near(gradient(&result),(peak(&run(&plus))-peak(&run(&minus)))/(2.0*h),2e-5);
        if mode>1 {
            let rad=result.get("radiation").unwrap();
            near(n(rad,"radiative_out_w")+n(rad,"convective_out_w"),1.0,1e-7);
            assert!(n(rad,"radiative_out_w")>0.01);
        }
    }
}

#[test]
fn reordered_faces_replay_and_swapped_sides_reverse_heat_not_physics() {
    let input=J::parse(BASE).unwrap();let original=output(&input);let parsed=success(&original);
    let mut reordered=input.clone();
    for key in ["side_a_faces","side_b_faces"] {
        rows(member(member(contact(&mut reordered),"nonmatching"),key)).reverse();
    }
    let replay=output(&reordered);success(&replay);assert_eq!(original.stdout,replay.stdout);
    let mut reverse=input;
    let sides=member(contact(&mut reverse),"nonmatching");
    let a=sides.get("side_a_faces").unwrap().clone();let b=sides.get("side_b_faces").unwrap().clone();
    put(sides,"side_a_faces",b);put(sides,"side_b_faces",a);
    put(contact(&mut reverse),"side_a_material",J::Str("substrate".into()));
    put(contact(&mut reverse),"side_b_material",J::Str("spreader".into()));
    let swapped=run(&reverse);near(peak(&parsed),peak(&swapped),2e-7);
    near(gradient(&parsed),gradient(&swapped),2e-7);
    let a=&parsed.get("contacts").unwrap().as_array().unwrap()[0];
    let b=&swapped.get("contacts").unwrap().as_array().unwrap()[0];
    near(n(a,"heat_a_to_b_w"),-n(b,"heat_a_to_b_w"),2e-7);
}

#[test]
fn radiating_repeated_trajectory_keeps_contact_in_forward_and_storage_adjoint() {
    let mut input=J::parse(BASE).unwrap();nonlinear(&mut input);radiation(&mut input);
    put(member(&mut input,"objective"),"gradient",J::Bool(false));
    put(&mut input,"transient",J::parse(r#"{"initial_temperature_k":300,"volumetric_heat_capacity_j_m3_k":2000000,"max_step_s":2,"max_steps":100,"intervals":[{"duration_s":2,"power_scale":1,"fan_speed_ratio":1},{"duration_s":2,"power_scale":0,"fan_speed_ratio":1.5}],"repeat":{"cycles":2,"max_total_steps":100},"nonlinear":{"max_iterations":32,"residual_rtol":1e-10,"residual_atol_j":1e-10,"armijo_c":1e-4,"shrink":0.5,"max_backtracks":24},"adjoint":{"qoi":"final","max_checkpoint_bytes":1048576}}"#).unwrap());
    let result=run(&input);let repeated=result.get("repeated_cycles").unwrap();
    near(n(repeated,"stored_energy_change_j")+n(repeated,"air_energy_gain_j")
        +n(repeated,"radiative_energy_loss_j"),n(repeated,"input_energy_j"),1e-5);
    let derivative=n(&repeated.path(&["adjoint","intervals"]).unwrap().as_array().unwrap()[0],
        "dtemperature_dpower_multiplier_k");
    let mut primal=input.clone();remove(member(&mut primal,"transient"),"adjoint");
    let plain=run(&primal);assert_eq!(plain.get("solid_temperatures_k"),result.get("solid_temperatures_k"));
    let mut plus=primal.clone();let mut minus=primal;
    let h=0.001;
    put(&mut rows(member(member(&mut plus,"transient"),"intervals"))[0],"power_scale",num(1.0+h));
    put(&mut rows(member(member(&mut minus,"transient"),"intervals"))[0],"power_scale",num(1.0-h));
    near(derivative,(peak(&run(&plus))-peak(&run(&minus)))/(2.0*h),2e-5);
}

#[test]
fn missing_coverage_real_gaps_and_exhausted_geometry_budgets_never_publish() {
    let input=J::parse(BASE).unwrap();let mut variants=Vec::new();
    let mut missing=input.clone();
    rows(member(member(contact(&mut missing),"nonmatching"),"side_b_faces")).pop();variants.push(missing);
    let mut gap=input.clone();
    for p in rows(member(member(&mut gap,"solid"),"vertices_m")).iter_mut().skip(8) {
        let x=rows(p)[0].as_f64().unwrap();rows(p)[0]=num(x+1e-4);
    }
    variants.push(gap);
    for key in ["max_pair_tests","max_overlap_triangles"] {
        let mut exhausted=input.clone();put(member(contact(&mut exhausted),"nonmatching"),key,num(1.0));variants.push(exhausted);
    }
    for variant in variants {
        let result=output(&variant);assert!(!result.status.success());assert!(result.stdout.is_empty());
    }
    let mut mesh=input.clone();put(member(&mut mesh,"objective"),"gradient",J::Bool(false));
    put(&mut mesh,"mesh_convergence",J::parse(r#"{"max_refinements":2,"consecutive_passes":2,"temperature_tolerance_k":1,"max_vertices":20000,"max_tetrahedra":100000}"#).unwrap());
    let result=output(&mesh);assert_eq!(result.status.code(),Some(6));assert!(result.stdout.is_empty());
    assert!(String::from_utf8_lossy(&result.stderr).contains("max_pair_tests=10000"));
}

static NEXT:AtomicU64=AtomicU64::new(0);
struct Scratch(PathBuf);
impl Scratch {
    fn new()->Self {
        let path=std::env::temp_dir().join(format!("fs-nonmatching-uq-{}-{}",std::process::id(),NEXT.fetch_add(1,Ordering::SeqCst)));
        std::fs::create_dir(&path).unwrap();Self(path)
    }
    fn uq(&self,plan:&str,extra:&[&str])->Output {
        std::fs::write(self.0.join("uq.json"),plan).unwrap();
        Command::new(env!("CARGO_BIN_EXE_frankensim")).current_dir(&self.0)
            .args(["--json","cooling-network-uq","base.json","uq.json"]).args(extra).output().unwrap()
    }
}
impl Drop for Scratch {fn drop(&mut self){let _=std::fs::remove_dir_all(&self.0);}}
#[test]
fn uncertain_nonmatching_resistance_changes_real_samples_and_resumes_exactly() {
    let dir=Scratch::new();let mut input=J::parse(BASE).unwrap();
    put(member(&mut input,"objective"),"gradient",J::Bool(false));
    std::fs::write(dir.0.join("base.json"),text(&input)).unwrap();
    let plan=r#"{"schema":"frankensim.cooling-network-uq.v1","seed":"73","samples":2,"wall_seconds":600,"temperature_limit_k":305,"correlation":{"kind":"independent"},"parameters":[{"target":{"kind":"contact-resistance","contact":"bondline"},"distribution":{"kind":"uniform","lo":0.005,"hi":0.025}}]}"#;
    let zero=plan.replace("\"lo\":0.005,\"hi\":0.025","\"lo\":0.01,\"hi\":0.01");
    let result=success(&dir.uq(&zero,&[]));near(n(&result,"mean_k"),peak(&run(&input)),1e-10);
    assert_eq!(n(&result,"std_dev_k"),0.0);
    let full=dir.uq(plan,&["--checkpoint","full.uqcp"]);let parsed=success(&full);
    assert!(n(&parsed,"std_dev_k")>1e-5);
    let part=dir.uq(plan,&["--checkpoint","part.uqcp","--max-new-samples","1"]);
    assert_eq!(part.status.code(),Some(6));
    let saved=std::fs::read(dir.0.join("part.uqcp")).unwrap();
    let done=dir.uq(plan,&["--resume","part.uqcp","--checkpoint","done.uqcp"]);success(&done);
    assert_eq!(full.stdout,done.stdout);
    assert_eq!(std::fs::read(dir.0.join("full.uqcp")).unwrap(),std::fs::read(dir.0.join("done.uqcp")).unwrap());
    assert_eq!(saved,std::fs::read(dir.0.join("part.uqcp")).unwrap());
}
