//! Run the actual cooling binary. Perturbed calls change complete physical
//! inputs and rerun the whole trajectory, not the adjoint's local residual.
#![cfg(unix)]
#[allow(dead_code)]
#[path="../src/json_read.rs"]
mod json;
use json::JsonValue as J;
use std::io::Write;
use std::process::{Command,Output,Stdio};
const BASE: &str = include_str!(concat!(env!("CARGO_MANIFEST_DIR"),
    "/../../examples/cooling-network/adjoint-component-contact-pulse.json"));
const NONMATCHING: &str = include_str!(concat!(env!("CARGO_MANIFEST_DIR"),
    "/../../examples/cooling-network/nonmatching-contact-hotspot.json"));
const COMPONENTS: [&str;3] = ["chip","memory","standby"];
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
fn run(input:&J)->J {
    let result=output(input);
    assert!(result.status.success(),"{}",String::from_utf8_lossy(&result.stderr));
    J::parse(std::str::from_utf8(&result.stdout).unwrap()).unwrap()
}
fn n(root:&J,key:&str)->f64 {root.f64_field(key).unwrap()}
fn near(a:f64,b:f64,tol:f64) {assert!((a-b).abs()<tol,"{a:e} versus {b:e}");}
fn section(doc:&J)->&J {doc.get("repeated_cycles").unwrap_or_else(||doc.get("transient").unwrap())}
fn adjoint(doc:&J)->&J {section(doc).get("adjoint").unwrap()}
fn value(doc:&J,qoi:&str)->f64 {
    if qoi=="final" {n(doc.get("objective").unwrap(),"value_k")}
    else {n(section(doc),"sampled_peak_objective_k")}
}
fn component<'a>(doc:&'a J,interval:usize,name:&str)->&'a J {
    adjoint(doc).path(&["component_power_sensitivities","intervals"]).unwrap().as_array().unwrap()[interval]
        .get("rows").unwrap().as_array().unwrap().iter()
        .find(|row|row.str_field("component")==Some(name)).unwrap()
}
fn contact_gradient(doc:&J,name:&str)->f64 {
    n(adjoint(doc).path(&["contact_resistance_sensitivities","rows"]).unwrap().as_array().unwrap().iter()
        .find(|row|row.str_field("contact")==Some(name)).unwrap(),"dtemperature_dlog_resistance_k")
}
fn primal(input:&J)->J {
    let mut input=input.clone();remove(member(&mut input,"transient"),"adjoint");input
}
fn case(qoi:&str)->J {
    let mut input=J::parse(BASE).unwrap();
    put(member(member(&mut input,"transient"),"adjoint"),"qoi",J::Str(qoi.into()));input
}
fn applied(input:&J,interval:usize)->J {
    let row=&input.path(&["transient","intervals"]).unwrap().as_array().unwrap()[interval];
    if let Some(powers)=row.get("component_powers_w") {return powers.clone();}
    let scale=n(row,"power_scale");
    J::Object(input.path(&["solid","component_power","components"]).unwrap().as_array().unwrap().iter()
        .map(|c|(c.str_field("name").unwrap().into(),num(scale*n(c,"watts")))).collect())
}
fn change_interval(input:&mut J,interval:usize,name:&str,watts:f64) {
    let mut powers=applied(input,interval);put(&mut powers,name,num(watts));
    let row=&mut rows(member(member(input,"transient"),"intervals"))[interval];
    remove(row,"power_scale");put(row,"component_powers_w",powers);
}
fn change_base(input:&mut J,name:&str,watts:f64) {
    let power=member(member(input,"solid"),"component_power");
    let components=rows(member(power,"components"));
    put(components.iter_mut().find(|c|c.str_field("name")==Some(name)).unwrap(),"watts",num(watts));
    let total=components.iter().map(|c|n(c,"watts")).sum();put(power,"total_w",num(total));
}
fn resistance(input:&mut J,value:f64) {
    put(&mut rows(member(member(input,"solid"),"contacts"))[0],"resistance_m2_k_w",num(value));
}
fn interval_fd(input:&J,qoi:&str,interval:usize,name:&str)->f64 {
    let baseline=applied(input,interval).get(name).unwrap().as_f64().unwrap();
    let h=1e-3;
    let mut plus=primal(input);change_interval(&mut plus,interval,name,baseline+h);
    let mut minus=primal(input);
    let denominator=if baseline>=h {
        change_interval(&mut minus,interval,name,baseline-h);2.0*h
    } else {h};
    (value(&run(&plus),qoi)-value(&run(&minus),qoi))/denominator
}

#[test]
fn every_component_has_its_own_absolute_watt_derivative_including_dormant_chips() {
    for qoi in ["final","sampled-peak"] {
        let input=case(qoi);let result=run(&input);
        for interval in 0..2 {
            let mut relative_sum=0.0;
            for name in COMPONENTS {
                let row=component(&result,interval,name);
                let actual=n(row,"dtemperature_dpower_w_k_per_w");
                near(n(row,"applied_power_w"),applied(&input,interval).get(name).unwrap().as_f64().unwrap(),1e-12);
                near(actual,interval_fd(&input,qoi,interval,name),5e-5);
                near(n(row,"dtemperature_dpower_multiplier_k"),actual*n(row,"applied_power_w"),1e-12);
                relative_sum+=n(row,"dtemperature_dpower_multiplier_k");
            }
            let old=&adjoint(&result).get("intervals").unwrap().as_array().unwrap()[interval];
            near(relative_sum,n(old,"dtemperature_dpower_multiplier_k"),1e-7);
        }
        let dormant=component(&result,0,"standby");
        assert_eq!(n(dormant,"applied_power_w"),0.0);
        assert_eq!(n(dormant,"dtemperature_dpower_multiplier_k"),0.0);
        assert!(n(dormant,"dtemperature_dpower_w_k_per_w").abs()>1e-4);
        if qoi=="sampled-peak" {
            // Independent P1 volume/face integration and direct implicit
            // transpose with analytic air elimination, not a second CLI call.
            near(value(&result,qoi),302.286668764270,2e-5);
            near(n(component(&result,0,"chip"),"dtemperature_dpower_w_k_per_w"),0.150444865984,2e-5);
            near(n(dormant,"dtemperature_dpower_w_k_per_w"),0.071983013651,2e-5);
        }
    }
}

#[test]
fn contacts_use_total_trajectory_adjoint_on_matching_and_nonmatching_radiating_solids() {
    for nonmatching in [false,true] {
        let mut input=case("final");
        if nonmatching {
            let mut independent=J::parse(NONMATCHING).unwrap();
            put(member(&mut independent,"objective"),"gradient",J::Bool(false));
            for key in ["component_power","materials"] {
                put(member(&mut independent,"solid"),key,input.get("solid").unwrap().get(key).unwrap().clone());
            }
            put(&mut independent,"radiation",input.get("radiation").unwrap().clone());
            let mut schedule=input.get("transient").unwrap().clone();
            remove(&mut schedule,"element_heat_capacities_j_m3_k");
            put(&mut schedule,"volumetric_heat_capacity_j_m3_k",num(2e6));
            put(&mut independent,"transient",schedule);input=independent;
        }
        let result=run(&input);let h=1e-3_f64;
        let mut plus=primal(&input);let mut minus=primal(&input);
        resistance(&mut plus,0.01*h.exp());resistance(&mut minus,0.01*(-h).exp());
        near(contact_gradient(&result,"bondline"),(value(&run(&plus),"final")-value(&run(&minus),"final"))/(2.0*h),3e-5);
        if nonmatching {
            near(n(component(&result,0,"standby"),"dtemperature_dpower_w_k_per_w"),
                interval_fd(&input,"final",0,"standby"),5e-5);
        }
        let mut ordinary=input.clone();
        let options=member(member(&mut ordinary,"transient"),"adjoint");
        remove(options,"component_power");remove(options,"contact_resistance");
        let existing=run(&ordinary);
        assert_eq!(existing.get("solid_temperatures_k"),result.get("solid_temperatures_k"));
        assert_eq!(existing.path(&["transient","history"]),result.path(&["transient","history"]));
        assert_eq!(section(&existing).get("cycles"),section(&result).get("cycles"));
        for key in ["reconstruction_solid_solves","reconstructed_solid_endpoints","adjoint_sweeps"] {
            assert_eq!(adjoint(&existing).get(key),adjoint(&result).get(key),"no solve per control");
        }
        assert!(adjoint(&existing).get("component_power_sensitivities").is_none());
        assert!(adjoint(&result).get("checkpoint_bytes").unwrap().as_f64().unwrap()
            >adjoint(&existing).get("checkpoint_bytes").unwrap().as_f64().unwrap());
    }
}

#[test]
fn repeated_shared_controls_equal_unrolled_phase_sums_and_contacts_can_be_independent() {
    let input=case("sampled-peak");let repeated=run(&input);
    let mut unrolled=input.clone();let schedule=member(&mut unrolled,"transient");
    remove(schedule,"repeat");
    let intervals=rows(member(schedule,"intervals"));let one=intervals.clone();intervals.extend(one);
    let explicit=run(&unrolled);
    assert_eq!(explicit.get("solid_temperatures_k"),repeated.get("solid_temperatures_k"));
    for interval in 0..2 {for name in COMPONENTS {
        near(n(component(&repeated,interval,name),"dtemperature_dpower_w_k_per_w"),
            n(component(&explicit,interval,name),"dtemperature_dpower_w_k_per_w")
                +n(component(&explicit,interval+2,name),"dtemperature_dpower_w_k_per_w"),1e-7);
    }}
    near(contact_gradient(&repeated,"bondline"),contact_gradient(&explicit,"bondline"),1e-7);
    let mut split=input.clone();let contacts=rows(member(member(&mut split,"solid"),"contacts"));
    let original=contacts[0].clone();let pairs=original.get("face_pairs").unwrap().as_array().unwrap();
    contacts.clear();
    for (i,pair) in pairs.iter().enumerate() {
        let mut row=original.clone();put(&mut row,"name",J::Str(format!("bond-{i}")));
        put(&mut row,"face_pairs",J::Array(vec![pair.clone()]));contacts.push(row);
    }
    let separated=run(&split);
    near(value(&separated,"sampled-peak"),value(&repeated,"sampled-peak"),2e-7);
    near(contact_gradient(&separated,"bond-0")+contact_gradient(&separated,"bond-1"),
        contact_gradient(&repeated,"bondline"),1e-7);
}

#[test]
fn base_watt_gradients_respect_absolute_workload_overrides() {
    let input=case("final");let result=run(&input);
    let base=adjoint(&result).path(&["component_power_sensitivities","base_rows"]).unwrap().as_array().unwrap();
    for row in base {
        let name=row.str_field("component").unwrap();let watts=n(row,"base_power_w");
        let expected=1.25*n(component(&result,0,name),"dtemperature_dpower_w_k_per_w");
        near(n(row,"dtemperature_dbase_power_w_k_per_w"),expected,1e-12);
        // The second interval explicitly overrides every base component.
        let h=1e-3;let mut plus=primal(&input);let mut minus=primal(&input);
        change_base(&mut plus,name,watts+h);
        let denominator=if watts>=h {change_base(&mut minus,name,watts-h);2.0*h}else{h};
        near(expected,(value(&run(&plus),"final")-value(&run(&minus),"final"))/denominator,5e-5);
    }
}

#[test]
fn initial_peaks_zero_future_controls_and_control_memory_is_admitted_before_publication() {
    let mut input=case("sampled-peak");
    put(&mut input,"objective",J::parse(r#"{"mean_wall_region":"first-face","gradient":false}"#).unwrap());
    let schedule=member(&mut input,"transient");put(schedule,"initial_temperature_k",num(350.0));
    for interval in rows(member(schedule,"intervals")) {
        remove(interval,"component_powers_w");put(interval,"power_scale",num(0.0));
    }
    let result=run(&input);assert_eq!(n(adjoint(&result),"state_index"),0.0);
    assert_eq!(n(adjoint(&result),"reconstructed_solid_endpoints"),0.0);
    assert_eq!(contact_gradient(&result,"bondline"),0.0);
    for interval in 0..2 {for name in COMPONENTS {
        assert_eq!(n(component(&result,interval,name),"dtemperature_dpower_w_k_per_w"),0.0);
    }}
    let mut plain=case("final");
    let options=member(member(&mut plain,"transient"),"adjoint");
    remove(options,"component_power");remove(options,"contact_resistance");
    let old=run(&plain);let bytes=n(adjoint(&old),"checkpoint_bytes");
    let mut exhausted=case("final");
    put(member(member(&mut exhausted,"transient"),"adjoint"),"max_checkpoint_bytes",num(bytes));
    let failure=output(&exhausted);assert_eq!(failure.status.code(),Some(6));assert!(failure.stdout.is_empty());
    assert!(String::from_utf8_lossy(&failure.stderr).contains("max_checkpoint_bytes"));
    let mut missing=case("final");
    let solid=member(&mut missing,"solid");remove(solid,"component_power");put(solid,"source_w_m3",num(0.0));
    for interval in rows(member(member(&mut missing,"transient"),"intervals")) {
        remove(interval,"component_powers_w");put(interval,"power_scale",num(0.0));
    }
    let rejected=output(&missing);assert!(!rejected.status.success());assert!(rejected.stdout.is_empty());
    assert!(String::from_utf8_lossy(&rejected.stderr).contains("component_power"));
    let mut malformed=case("final");
    put(member(member(&mut malformed,"transient"),"adjoint"),"component_power",num(1.0));
    assert!(!output(&malformed).status.success());
}
