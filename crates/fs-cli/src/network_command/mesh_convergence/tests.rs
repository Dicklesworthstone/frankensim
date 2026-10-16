use super::*;
const HOTSPOT:&str=include_str!(concat!(env!("CARGO_MANIFEST_DIR"),"/../../examples/cooling-network/fan-hotspot.json"));
const CONTACT:&str=include_str!(concat!(env!("CARGO_MANIFEST_DIR"),"/../../examples/cooling-network/adjoint-contact-pulse.json"));
fn limits()->TetRefinementLimits{TetRefinementLimits{max_vertices:20000,max_tetrahedra:100000}}
fn policy()->J{J::parse(r#"{"max_refinements":2,"consecutive_passes":2,"temperature_tolerance_k":100,"max_vertices":20000,"max_tetrahedra":100000}"#).unwrap()}
fn source_integral(root:&J,request:&Request)->f64 {
    let p=request.mesh.positions();
    array(get(get(root,"solid").unwrap(),"tetrahedra").unwrap(),"tets",100000).unwrap().iter().map(|t| {
        let t=indices::<4>(t,"tet",p.len()).unwrap();
        let a:[f64;3]=std::array::from_fn(|i|p[t[1] as usize][i]-p[t[0] as usize][i]);
        let b:[f64;3]=std::array::from_fn(|i|p[t[2] as usize][i]-p[t[0] as usize][i]);
        let c:[f64;3]=std::array::from_fn(|i|p[t[3] as usize][i]-p[t[0] as usize][i]);
        let volume=(a[0]*(b[1]*c[2]-b[2]*c[1])-a[1]*(b[0]*c[2]-b[2]*c[0])+a[2]*(b[0]*c[1]-b[1]*c[0])).abs()/6.;
        volume*t.iter().map(|&v|request.solid_data.nodal_source.as_ref().map_or(request.source,|s|s.at(v as usize))).sum::<f64>()/4.
    }).sum()
}
#[test]
fn refinement_prolongs_the_original_power_density_instead_of_shrinking_its_footprint() {
    let mut root=J::parse(HOTSPOT).unwrap();
    let mut request=Request::parse(HOTSPOT).unwrap();
    let original_density=request.solid_data.nodal_source.as_ref().unwrap().at(4);
    let original_component=get(get(&root,"solid").unwrap(),"component_power").unwrap().clone();
    let gate=CancelGate::new_clock_free();
    for _ in 0..2 {
        root=with_context(&request,&gate,|cx|refine_request(cx,&root,&request,limits())).unwrap();
        request=Request::parse(&encode(&root).unwrap()).unwrap();
        assert_eq!(request.solid_data.nodal_source.as_ref().unwrap().at(4),original_density);
        assert!((source_integral(&root,&request)-1.).abs()<1e-12);
        assert!(get(&root,"solid").unwrap().get("component_power").is_none());
    }
    let mut wrong=root.clone();let solid=member_mut(&mut wrong,"solid").unwrap();
    members(solid).unwrap().retain(|(key,_)|key!="nodal_source_w_m3");
    members(solid).unwrap().push(("component_power".into(),original_component));
    let wrong=Request::parse(&encode(&wrong).unwrap()).unwrap();
    assert!(wrong.solid_data.nodal_source.as_ref().unwrap().at(4)>original_density*10.,
        "the old vertex footprint must be a detectable changed-load control");
}
#[test]
fn nonlinear_material_assignments_and_reordered_contact_pairs_survive_refinement() {
    let mut root=J::parse(CONTACT).unwrap();members(&mut root).unwrap().retain(|(key,_)|key!="transient");
    let J::Array(contacts)=member_mut(member_mut(&mut root,"solid").unwrap(),"contacts").unwrap() else{panic!()};
    let J::Array(pairs)=member_mut(&mut contacts[0],"face_pairs").unwrap() else{panic!()};
    for pair in pairs {let J::Array(vertices)=member_mut(pair,"side_b").unwrap() else{panic!()};vertices.reverse();}
    let mut request=Request::parse(&encode(&root).unwrap()).unwrap();
    let original_materials=root.path(&["solid","materials"]).unwrap().clone();
    let gate=CancelGate::new_clock_free();
    for generation in 1..=2 {
        root=with_context(&request,&gate,|cx|refine_request(cx,&root,&request,limits())).unwrap();
        request=Request::parse(&encode(&root).unwrap()).unwrap();
        assert_eq!(root.path(&["solid","materials"]),Some(&original_materials));
        assert_eq!(request.mesh.element_count(),12*8_usize.pow(generation));
        assert!(request.contacts.is_some());
        assert!((source_integral(&root,&request)-20.).abs()<1e-10);
        let contacts=root.path(&["solid","contacts"]).unwrap().as_array().unwrap();
        assert_eq!(contacts[0].get("face_pairs").unwrap().as_array().unwrap().len(),2*4_usize.pow(generation));
    }
}
#[test]
fn actual_rungs_replay_the_resolved_request_and_need_two_measured_changes() {
    let mut root=J::parse(HOTSPOT).unwrap();members(&mut root).unwrap().push(("mesh_convergence".into(),policy()));
    let request=Request::parse(&encode(&root).unwrap()).unwrap();let gate=CancelGate::new_clock_free();
    let result=J::parse(&execute(&request,&gate).unwrap()).unwrap();
    let study=result.get("mesh_convergence").unwrap();
    assert_eq!(study.f64_field("meshes_solved"),Some(3.));
    let rows=study.get("history").unwrap().as_array().unwrap();
    for(i,row)in rows.iter().enumerate(){
        assert_eq!(row.f64_field("tetrahedra"),Some((12*8_usize.pow(i as u32)) as f64));
        assert!((row.f64_field("source_w").unwrap()-1.).abs()<1e-7);
    }
    let resolved=study.get("resolved_request").unwrap();assert!(resolved.get("mesh_convergence").is_none());
    let replay=Request::parse(&encode(resolved).unwrap()).unwrap();
    let replay=J::parse(&execute(&replay,&gate).unwrap()).unwrap();
    assert_eq!(result.get("solid_temperatures_k"),replay.get("solid_temperatures_k"));
    assert_eq!(result.get("objective"),replay.get("objective"));
}
#[test]
fn policy_and_cancellation_cannot_publish_an_unmeasured_convergence() {
    let mut root=J::parse(HOTSPOT).unwrap();
    let invalid=J::parse(r#"{"max_refinements":1,"consecutive_passes":1,"temperature_tolerance_k":100,"max_vertices":20000,"max_tetrahedra":100000}"#).unwrap();
    assert!(Study::parse(&invalid,&root).is_err());
    members(&mut root).unwrap().push(("mesh_convergence".into(),policy()));
    let request=Request::parse(&encode(&root).unwrap()).unwrap();let gate=CancelGate::new_clock_free();gate.request();
    assert_eq!(execute(&request,&gate).unwrap_err().code,"cooling-network-cancelled");
    replace(member_mut(&mut root,"objective").unwrap(),"gradient",J::Bool(true)).unwrap();
    assert!(Request::parse(&encode(&root).unwrap()).is_err());
}
#[test]
fn nodal_source_mode_is_explicit_finite_and_mutually_exclusive() {
    let mut root=J::parse(HOTSPOT).unwrap();let solid=member_mut(&mut root,"solid").unwrap();
    members(solid).unwrap().push(("nodal_source_w_m3".into(),J::Array(vec![jnum(1000.);12])));
    assert!(Request::parse(&encode(&root).unwrap()).is_err());
    members(member_mut(&mut root,"solid").unwrap()).unwrap().retain(|(key,_)|key!="component_power");
    let request=Request::parse(&encode(&root).unwrap()).unwrap();
    assert!((source_integral(&root,&request)-0.5).abs()<1e-12);
    assert!(request.solid_data.render(request.conductivity,request.source).unwrap().contains("nodal-p1-density"));
    replace(member_mut(&mut root,"solid").unwrap(),"nodal_source_w_m3",J::Array(vec![jnum(1.)])).unwrap();
    assert!(Request::parse(&encode(&root).unwrap()).is_err());
}
