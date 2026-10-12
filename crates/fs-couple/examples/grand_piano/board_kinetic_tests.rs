use super::*;

fn fixture() -> String {
    let mesh=PlateMesh::rectangle(0.7,0.45,2,2);
    let mut text=format!("{HEADER}\nsource,estimated,kinetic regression NOT measured piano\nsupport,clamped\npretension,0\ndamping,0.015\n");
    for (i,(x,y)) in mesh.nodes.iter().enumerate() {
        text.push_str(&format!("node,{i},{x},{y}\n"));
    }
    for (i,t) in mesh.tris.iter().enumerate() {
        text.push_str(&format!("triangle,{i},{},{},{},0.009,450,11000000000,800000000,0.3,600000000,0\n",t[0],t[1],t[2]));
    }
    // Clamp one edge only: a small complete pencil with nontrivial rotations.
    for (i,(_,y)) in mesh.nodes.iter().enumerate() {
        if *y==0.0 {text.push_str(&format!("fixed,{i}\n"));}
    }
    text.push_str("bridge,69,4,0.2,0.3,0.5\n");text
}
fn kinetic() -> BoardGeometry {
    BoardGeometry::read(&format!("{}panel-rotary-mass,consistent-cubic\n",fixture())).unwrap()
}
fn model(board: &BoardGeometry) -> PlateModel {
    board.chart.assemble_edge_cubic_kinetic_mass(&board.stiffeners,&board.supports,board.stiffener_mass).unwrap()
}
fn vector(model: &PlateModel, shapes: &[[f64;6]]) -> Vec<f64> {
    let mut q=vec![0.0;model.free];
    for (node,shape) in shapes.iter().enumerate() {
        for (c,value) in [shape[2],-shape[4],shape[3]].into_iter().enumerate() {
            if let Some(index)=model.dof_map[3*node+c] {q[index]=value;}
        }
    }
    q
}
fn dot(a:&[f64],b:&[f64])->f64 {a.iter().zip(b).map(|(a,b)|a*b).sum()}
fn close(a:f64,b:f64,scale:f64) {
    assert!((a-b).abs()<2e-8*scale.max(a.abs()).max(b.abs()),"{a} != {b}");
}
fn shared_projection(board:&BoardGeometry,prepared:&PreparedBoard) {
    let motion=prepared.motion.as_ref().unwrap();assert!(motion.is_edge_cubic());
    let site=&board.bridge_sites[0];let tri=board.chart.mesh.tris[site.triangle];
    let point=[(0..3).map(|i|site.weights[i]*board.chart.mesh.nodes[tri[i]].0).sum(),
        (0..3).map(|i|site.weights[i]*board.chart.mesh.nodes[tri[i]].1).sum(),0.0];
    let bridge=motion.normal_weights(point,[0.,0.,1.],0.).unwrap();
    for (i,mode) in prepared.modes.iter().enumerate() {
        close(bridge[i],mode.bridge[48],1.0);
        let volume: f64=prepared.surface.iter().map(|p|p.area_m2*p.mode_shape[i]).sum();
        close(volume,mode.volume,1.0);
    }
    for point in &prepared.surface {
        let weights=motion.normal_weights(point.position_m,[0.,0.,1.],0.).unwrap();
        for (&a,&b) in point.mode_shape.iter().zip(&weights) {close(a,b,1.0);}
    }
}

#[test]
fn rotary_selection_is_explicit_and_preserves_the_original_default_image() {
    let text=fixture();
    let old=BoardGeometry::read(&text).unwrap();
    let explicit=BoardGeometry::read(&format!("{text}panel-rotary-mass,lumped\n")).unwrap();
    assert!(!old.consistent_rotary_mass && !explicit.consistent_rotary_mass);
    for cubic in [false,true] {
        let a=old.prepare_inner(&[69],20_000.,true,true,false,cubic).unwrap();
        let b=explicit.prepare_inner(&[69],20_000.,true,true,false,cubic).unwrap();
        assert_eq!(a.frequency_intervals_hz,b.frequency_intervals_hz);
        assert_eq!(a.provenance,b.provenance);
        assert_eq!(a.motion.unwrap().shapes,b.motion.unwrap().shapes);
    }
    for row in ["panel-rotary-mass\n","panel-rotary-mass,other\n",
        "panel-rotary-mass,consistent-cubic,extra\n",
        "panel-rotary-mass,lumped\npanel-rotary-mass,consistent-cubic\n"] {
        assert!(BoardGeometry::read(&format!("{text}{row}")).is_err(),"{row}");
    }
    let board=kinetic();
    for retain_motion in [false,true] {for equilibrated in [false,true] {for linear in [false,true] {
        let error=board.prepare_inner(&[69],20_000.,retain_motion,equilibrated,linear,false).unwrap_err();
        assert!(error.contains("requires the edge-cubic"),"{error}");
    }}}
}

#[test]
fn full_kinetic_pencil_reaches_mass_normalized_modes_bridges_and_spatial_motion() {
    let board=kinetic();let model=model(&board);
    let old=BoardGeometry::read(&fixture()).unwrap()
        .prepare_with_motion_edge_cubic_transverse_mass(&[69],20_000.,true).unwrap();
    let prepared=board.prepare_with_motion_edge_cubic_transverse_mass(&[69],20_000.,true).unwrap();
    assert_eq!(prepared.modes.len(),model.free); // complete small physical pencil
    assert_eq!(prepared.mass_kg,old.mass_kg);assert_eq!(prepared.area_m2,old.area_m2);
    assert!(prepared.provenance.contains("analytic-gradient rotary inertia"));
    assert!(!prepared.provenance.contains("lumped rotary"));
    assert!(prepared.modes.iter().zip(&old.modes).any(|(a,b)|
        (a.frequency_hz/b.frequency_hz-1.0).abs()>1e-7));
    let motion=prepared.motion.as_ref().unwrap();
    let vectors:Vec<_>=motion.shapes.iter().map(|s|vector(&model,s)).collect();
    for (i,q) in vectors.iter().enumerate() {
        let mut mq=vec![0.;model.free];let mut kq=mq.clone();
        model.m.spmv(q,&mut mq);model.k.spmv(q,&mut kq);
        let lambda=(TAU*prepared.modes[i].frequency_hz).powi(2);
        for (j,p) in vectors.iter().enumerate() {close(dot(p,&mq),if i==j {1.}else{0.},1.0);}
        let residual=kq.iter().zip(&mq).map(|(k,m)|(k-lambda*m).powi(2)).sum::<f64>().sqrt();
        let scale=dot(&kq,&kq).sqrt()+lambda*dot(&mq,&mq).sqrt();
        assert!(residual/scale<1e-7,"mode {i}: {}",residual/scale);
    }
    shared_projection(&board,&prepared);
    let ordinary=board.prepare_edge_cubic_transverse_mass(&[69],20_000.,true).unwrap();
    assert!(ordinary.motion.is_none());
    assert_eq!(ordinary.frequency_intervals_hz,prepared.frequency_intervals_hz);
    for (a,b) in ordinary.modes.iter().zip(&prepared.modes) {
        assert_eq!(a.bridge,b.bridge);assert_eq!(a.volume,b.volume);
    }
}

#[test]
fn kinetic_basis_and_material_damping_survive_reduction_and_acoustic_refinement() {
    let board=kinetic();let model=model(&board);
    let source=board.prepare_with_motion_edge_cubic_transverse_mass(&[69],20_000.,true).unwrap();
    let source_q:Vec<_>=source.motion.as_ref().unwrap().shapes.iter().map(|s|vector(&model,s)).collect();
    let options=ritz::RitzOptions {max_modes:5,keep_low_modes:1,
        sample_hz:vec![100.,700.],transfer_tolerance:None};
    let reduced=board.prepare_reduced(&[69],20_000.,true,true,false,true,&options).unwrap();
    let r=reduced.modes.len();assert!(r>1 && r<source.modes.len());
    let q:Vec<_>=reduced.motion.as_ref().unwrap().shapes.iter().map(|s|vector(&model,s)).collect();
    let mq:Vec<_>=q.iter().map(|v| {let mut out=vec![0.;model.free];model.m.spmv(v,&mut out);out}).collect();
    let c=reduced.physical_damping.as_ref().unwrap();
    for i in 0..r {for j in 0..r {
        close(dot(&q[i],&mq[j]),if i==j {1.}else{0.},1.0);
        let expected:f64=source_q.iter().zip(&source.modes).map(|(v,m)|
            2.0*board.damping_ratio*TAU*m.frequency_hz*dot(v,&mq[i])*dot(v,&mq[j])).sum();
        close(c[r*i+j],expected,c[r*i+i].abs().max(c[r*j+j].abs()));
    }}
    shared_projection(&board,&reduced);
    let refined=kinetic().with_acoustic_refinement(1).unwrap()
        .prepare_reduced(&[69],20_000.,true,true,false,true,&options).unwrap();
    assert_eq!(refined.frequency_intervals_hz,reduced.frequency_intervals_hz);
    assert_eq!(refined.physical_damping,reduced.physical_damping);
    assert_eq!(refined.motion.as_ref().unwrap().shapes,reduced.motion.as_ref().unwrap().shapes);
    assert_eq!(refined.surface.len(),4*reduced.surface.len());
    shared_projection(&board,&refined);
}

#[test]
fn crowned_models_do_not_silently_discard_the_flat_kinetic_law() {
    use super::super::crowned_board;
    let source=fixture();let heights=vec![0.0;9];
    let crowned=crowned_board::elevate(&source,&heights,"test reference plane").unwrap();
    for choice in ["lumped","consistent-cubic"] {
        let input=format!("{crowned}panel-rotary-mass,{choice}\n");
        let error=crowned_board::CrownedBoard::read(&input).err().expect("flat inertia must refuse");
        assert!(error.contains("six-DOF kinetic field"),"{error}");
    }
}
