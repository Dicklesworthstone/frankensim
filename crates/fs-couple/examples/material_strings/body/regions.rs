//! Regions and spatial contacts use the SAME pack, plate and performer owners.
use super::*;

fn regional() -> String { include_str!("../regional.body").into() }
fn regional_ready(cx:&Cx<'_>,gate:&CancelGate,s:input::Spec,text:&str,right_loss:f64)->Ready {
    let profile=Profile::parse(text.as_bytes()).unwrap();
    s.compile_with_body(cx,gate,Some(profile),|name| {
        // Two file paths can name distinct original body packs with the SAME
        // elastic coefficients. Only the selected right-hand viscosity changes.
        Ok(if name=="right.fsmatpk" {pack("body.fsmatpk",right_loss,None)}
            else {pack(name,1.0,None)})
    }).unwrap()
}
fn close(a:&[f64],b:&[f64]) {
    assert_eq!(a.len(),b.len());
    for (&a,&b) in a.iter().zip(b) {
        assert!((a-b).abs()<=1e-11*(a.abs()+b.abs()).max(1e-12),"{a} != {b}");
    }
}

#[test]
fn subdividing_one_material_retains_actual_body_operators_and_pressure() {
    scoped(|cx,gate| {
        let mut original=compile(cx,gate,Some(profile()),1.0);
        let text=regional().replace("1 -0.23 2 4","1 0.17 2 4");
        let mut divided=regional_ready(cx,gate,spec(),&text,1.0);
        let source=divided.performance.bank().source().bank();
        let chart=source.material_chart().unwrap();
        let mut triangles:Vec<_>=chart.specimen.regions().iter().flat_map(|r|
            r.input.region.triangle_indices.iter().copied()).collect();
        triangles.sort_unstable();assert_eq!(triangles,(0..32).collect::<Vec<_>>());
        assert_eq!(divided.catalog.len(),2,"regions must reuse the original decoded pack");
        let ca=original.performance.bank().source().bank().material_chart_viscosity().unwrap().damping().matrix_per_s();
        let cb=source.material_chart_viscosity().unwrap().damping().matrix_per_s();close(ca,cb);
        let ma=original.performance.bank().system().components()[1].modes();
        let mb=divided.performance.bank().system().components()[1].modes();assert_eq!(ma,mb);
        close(&pressure(&mut original,128),&pressure(&mut divided,37));
    });
}

#[test]
fn changing_only_one_regions_viscosity_changes_sound_without_replacing_other_sources() {
    scoped(|cx,gate| {
        let text=regional().replace("region right @body","region right right.fsmatpk");
        let mut a=regional_ready(cx,gate,spec(),&text,1.0);
        let mut b=regional_ready(cx,gate,spec(),&text,9.0);
        assert_eq!(a.spec.bytes,b.spec.bytes);
        assert_eq!(a.body_profile.as_ref().unwrap().bytes,b.body_profile.as_ref().unwrap().bytes);
        assert_eq!(a.catalog.len(),3);assert_eq!(b.catalog.len(),3);
        assert_eq!(a.catalog["wire.fsmatpk"].pack.content_hash(),b.catalog["wire.fsmatpk"].pack.content_hash());
        assert_eq!(a.catalog["body.fsmatpk"].pack.content_hash(),b.catalog["body.fsmatpk"].pack.content_hash());
        assert_ne!(a.catalog["right.fsmatpk"].pack.content_hash(),b.catalog["right.fsmatpk"].pack.content_hash());
        let sa=a.performance.bank().source().bank();let sb=b.performance.bank().source().bank();
        let ra=sa.material_chart().unwrap().specimen.regions();
        let rb=sb.material_chart().unwrap().specimen.regions();
        assert_eq!(ra[0].input.material.identity(),rb[0].input.material.identity());
        assert_ne!(ra[1].input.material.identity(),rb[1].input.material.identity());
        assert_eq!(a.performance.bank().system().components()[1].modes(),
            b.performance.bank().system().components()[1].modes());
        assert_ne!(sa.material_chart_viscosity().unwrap().damping().matrix_per_s(),
            sb.material_chart_viscosity().unwrap().damping().matrix_per_s());
        assert!(pressure(&mut a,128).iter().zip(pressure(&mut b,37)).any(|(a,b)|(a-b).abs()>1e-12));
    });
}

#[test]
fn each_string_uses_its_own_body_location_on_one_retained_modal_basis() {
    scoped(|cx,gate| {
        let mut s=spec();s.strings.push(s.strings[0].clone());
        // The second string is not bowed or directly forced. Its source and
        // body connection remain in the single reciprocal mechanical state.
        let text=regional().replace("bridges 1\nbridge 2 2","bridges 2\nbridge 1 2\nbridge 3 2");
        let mut a=regional_ready(cx,gate,s.clone(),&text,1.0);
        let mut b=regional_ready(cx,gate,s.clone(),&text,1.0);
        let together=text.replace("bridge 3 2","bridge 1 2");
        let mut c=regional_ready(cx,gate,s,&together,1.0);
        let src=a.performance.bank().source().bank();
        let maps=src.material_chart_bridges().unwrap();
        assert_eq!(maps.weights.len(),2);assert_eq!(maps.weights[0][11],1.0);assert_eq!(maps.weights[1][13],1.0);
        assert_ne!(src.strings()[0].connections[0].body.shapes,src.strings()[1].connections[0].body.shapes);
        // An attachment-location change must not create an independent body or
        // retune its free eigenmodes. It changes the coupled reaction instead.
        assert_eq!(a.performance.bank().system().components().len(),4);
        assert_eq!(a.performance.bank().system().components()[2].modes(),c.performance.bank().system().components()[2].modes());
        let pa=pressure(&mut a,1);assert_eq!(pa,pressure(&mut b,128));
        assert!(pa.iter().zip(pressure(&mut c,37)).any(|(a,b)|(a-b).abs()>1e-12));
        assert!(a.performance.bank().system().components()[1].states().iter().any(|s|
            s.displacement_m_sqrt_kg!=0.0 || s.velocity_m_sqrt_kg_per_s!=0.0));
        let frame=a.performance.bank().system().last_frame().unwrap();
        assert!(frame.energy_residual_j.abs()<=frame.energy_tolerance_j);
    });
}

#[test]
fn regional_coverage_thickness_paths_and_performer_binding_cannot_be_silently_inferred() {
    let text=regional();
    for (from,to) in [("regions 2","regions 999999999"),("region right","region left"),
        ("2 4 0 4","1 4 0 4"),("2 4 0 4","3 4 0 4"),("2 4 0 4","2 5 0 4"),
        ("orthotropic12 1 -0.23","orthotropic12 0 -0.23"),
        ("region right @body","region right ../body.fsmatpk"),
        ("@body","body.fsmatpk"),("bridges 1","bridges 9"),("bridge 2 2","bridge 4 2")] {
        assert!(Profile::parse(text.replace(from,to).as_bytes()).is_err(),"accepted {to}");
    }
    scoped(|cx,gate| {
        let wrong=text.replace("bridges 1\nbridge 2 2","bridges 2\nbridge 1 2\nbridge 3 2");
        assert!(spec().compile_with_body(cx,gate,Some(Profile::parse(wrong.as_bytes()).unwrap()),
            |n|Ok(pack(n,1.0,None))).is_err());
        let separate=text.replace("region right @body","region right right.fsmatpk");
        for key in ["young_modulus_2",VISCOSITY[4]] {
            assert!(spec().compile_with_body(cx,gate,Some(Profile::parse(separate.as_bytes()).unwrap()),|n|
                Ok(if n=="right.fsmatpk" {pack("body.fsmatpk",1.0,Some(key))} else {pack(n,1.0,None)})).is_err());
        }
        let overflow=text.replace("orthotropic12 1 -0.23","orthotropic12 1e308 -0.23");
        let mut large=spec();large.body.thickness=10.0;
        assert!(large.compile_with_body(cx,gate,Some(Profile::parse(overflow.as_bytes()).unwrap()),
            |n|Ok(pack(n,1.0,None))).is_err());
    });
}

#[test]
fn regional_node_and_cell_maps_are_the_original_rectangle_geometry() {
    let BodyModel::Plate(p)=Profile::parse(regional().as_bytes()).unwrap().model else {
        panic!("expected regional plate profile");
    };
    let mesh=PlateMesh::rectangle(0.4,0.3,p.nx,p.ny);
    for r in p.region_inputs() {
        for e in p.triangles(&r) {
            let t=mesh.tris[e];
            let x=t.iter().map(|i|mesh.nodes[*i].0).sum::<f64>()/3.0;
            let y=t.iter().map(|i|mesh.nodes[*i].1).sum::<f64>()/3.0;
            assert!(x>0.4*r.x.0 as f64/p.nx as f64 && x<0.4*r.x.1 as f64/p.nx as f64);
            assert!(y>0.3*r.y.0 as f64/p.ny as f64 && y<0.3*r.y.1 as f64/p.ny as f64);
        }
    }
    let weights=p.footprint((1,2));let i=weights.iter().position(|w|*w==1.0).unwrap();
    assert_eq!(mesh.nodes[i],(0.1,0.15));assert_eq!(weights.iter().sum::<f64>(),1.0);
    // V1 keeps the first profile's ordinary ASCII-whitespace acceptance.
    let padded=include_str!("../orthotropic.body").replace("frankensim-string-body-v1","  frankensim-string-body-v1  ");
    assert!(Profile::parse(padded.as_bytes()).is_ok());
}
