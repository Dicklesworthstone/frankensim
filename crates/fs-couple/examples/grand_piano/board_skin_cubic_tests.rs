use super::*;

// Independent exact mean of a Cartesian monomial over an arbitrary triangle.
// Expand each coordinate in barycentric form and integrate every product by
// the simplex factorial identity, without evaluating any motion shape/rule.
fn monomial_mean(vertices: &[[f64;3];3], axes: &[usize]) -> f64 {
    let factorial=[1.,1.,2.,6.,24.,120.];
    let mut value=0.;
    for mut assignment in 0..3_usize.pow(axes.len() as u32) {
        let mut counts=[0;3];let mut coefficient=1.;
        for &axis in axes {
            let vertex=assignment%3;assignment/=3;
            counts[vertex]+=1;coefficient*=vertices[vertex][axis];
        }
        value+=coefficient*2.*counts.iter().map(|&n|factorial[n]).product::<f64>()
            /factorial[axes.len()+2];
    }
    value
}

#[test]
fn cubic_skin_integrates_physical_motion_on_walls_subtriangles_and_varying_offsets() {
    let mesh=ShellMesh::new(vec![[0.,0.,0.],[1.,0.,0.],[0.,1.,0.]],vec![[0,1,2]]).unwrap();
    // Actual reconstruction: w=x^3-(1-x-y)*x*y; theta=(w_y,-w_x,0).
    let motion=MotionSurface::new_edge_cubic(mesh.clone(),
        vec![vec![[0.;6],[0.,0.,1.,0.,-3.,0.],[0.;6]]]).unwrap();
    let closed=Skin::build(&mesh,&[0.012],0.02,100).unwrap();
    let bary=[[0.7,0.1,0.2],[0.3,0.6,0.1],[0.25,0.2,0.55]];
    let z=[0.004,0.009,0.006];
    // Isolate a valid panel embedding inside one source facet, with a varying
    // physical arm. This projection test makes no closed-boundary claim for
    // the cropped panel; the preceding case exercises the generated closure.
    let cropped=Skin {
        vertices:(0..3).map(|i|[bary[i][1],bary[i][2],z[i]]).collect(),
        triangles:vec![[0,1,2]],embeddings:vec![Embedding {element:0,bary,z}],
        source_triangles:mesh.tris.clone(),source_nodes:mesh.nodes.clone(),
        section_volume_m3:0.,volume_m3:0.,maximum_thickness_change_m:0.,maximum_offset_m:0.009,
    };
    let mut old_rule_error=0.0_f64;
    for skin in [closed,cropped] {
        let triangles=skin.panel_triangles();
        let normals:Vec<_>=triangles.iter().map(|&[a,b,c]| {
            let n=cross(sub(b,a),sub(c,a));let length=dot(n,n).sqrt();n.map(|v|v/length)
        }).collect();
        let actual=skin.normal_weights(&motion,&normals).unwrap();
        for (f,vertices) in triangles.iter().enumerate() {
            let m=|axes:&[usize]|monomial_mean(vertices,axes);
            let mean=[
                -3.*m(&[0,0,2])-2.*m(&[0,1,2])-m(&[1,1,2])+m(&[1,2]),
                -m(&[0,0,2])-2.*m(&[0,1,2])+m(&[0,2]),
                m(&[0,0,0])+m(&[0,0,1])+m(&[0,1,1])-m(&[0,1]),
            ];
            let expected=dot(normals[f],mean);
            assert!((actual[0][f]-expected).abs()<2e-14,
                "panel {f} mean {} != exact polynomial {expected}",actual[0][f]);
            let site=&skin.embeddings[f];
            let old: f64=QUADRATURE.iter().map(|q| {
                let bary=std::array::from_fn(|i|(0..3).map(|j|q[j]*site.bary[j][i]).sum());
                let z=(0..3).map(|j|q[j]*site.z[j]).sum::<f64>();
                motion.project_at(site.element,bary,[0.,0.,z],normals[f]).unwrap().0[0]/3.
            }).sum();
            old_rule_error=old_rule_error.max((old-expected).abs());
        }
    }
    assert!(old_rule_error>1e-7,"the case must detect degree-two quadrature used for cubic skin motion");
}
