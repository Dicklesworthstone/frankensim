use super::*;
use super::super::acceleration_fields;
use fs_bem::{panel3d::SpherePanels,helmholtz::{Medium,Formulation,solve_radiation_batch,exterior_pressure_at_points}};
use fs_couple::{vibroacoustic::CavityModes,render::plate::impact::cavity::{CavityCoupling,neck::CavityNeck}};
use fs_exec::CancelGate;
use fs_math::c64::C64;
use fs_plate::{ModePair,shell::{head::{TensionedDisk,TensionedDiskSpec},profile::ProfileBudget}};

fn boundary()->Boundary {
    let film=||TensionedDisk::new(TensionedDiskSpec {radius_m:0.1,thickness_m:0.0002,
        young_pa:4e9,poisson:0.38,density_kg_m3:1390.0,tension_n_m:1000.0,
        radial_intervals:2,azimuths:8},ProfileBudget {max_nodes:100,max_triangles:200,max_feature_evaluations:0}).unwrap();
    let films=[film(),film()];
    let modes=films.iter().map(|f| {
        let mut phi=vec![0.0;f.model.free];
        for (i,&(x,y)) in f.mesh.nodes.iter().enumerate() {if let Some(k)=f.model.dof_map[3*i] {
            phi[k]=1.0-(x*x+y*y)/0.01;
        }}
        vec![ModePair {lambda:1.0,phi,residual:0.0,interval:(1.0,1.0)}]
    }).collect::<Vec<_>>();
    // Kinematic surface fixture, not an eigenfrequency/calibration claim.
    Boundary::drum(&films,&modes,0.12,0.11).unwrap()
}
fn cavity()->CavityCoupling {
    let air=CavityModes {omegas:vec![0.0],lambdas:vec![0.004],interface:vec![vec![1.0]],
        loss_factor:0.0,rho0:1.2,c0:343.0};
    CavityCoupling::new(&air,3,&[0.0,-0.02,0.02],&[0.0]).unwrap().with_necks(vec![CavityNeck {
        area_m2:PI*0.003_f64.powi(2),effective_length_m:0.008,resistance_pa_s_m3:1000.0,
        pressure_shape_averages:vec![1.0],initial_volume_m3:0.0,initial_flow_m3_s:0.0,
    }],&CancelGate::new_clock_free()).unwrap()
}
fn port()->NeckRadiationPort {cavity().neck_radiation_port(0).unwrap()}
fn vented()->Boundary {boundary().with_sidewall_aperture(port(),PI/8.0,0.06,0.12).unwrap()}
fn flux(b:&Boundary,input:usize)->f64 {
    b.weights[input].iter().zip(&b.triangles).map(|(w,t)|w*area(*t)).sum()
}

#[test]
fn watertight_aperture_preserves_heads_and_projects_actual_neck_volume_flow_once() {
    let original=boundary();let added=vented();let port=port();let cavity=cavity();
    assert_eq!(added.state_modes,[1,2,3]);
    for i in 0..2 {assert!((flux(&added,i)-flux(&original,i)).abs()<1e-15);}
    assert!((flux(&added,2)-port.volume_weight_m2_per_sqrt_kg).abs()<1e-14);
    for (old,t) in original.triangles.iter().enumerate() {
        if original.weights.iter().any(|r|r[old]!=0.0) {
            let new=added.triangles.iter().position(|n|n==t).unwrap();
            for i in 0..2 {assert_eq!(added.weights[i][new],original.weights[i][old]);}
            assert_eq!(added.weights[2][new],0.0);
        }
    }
    let a=SpherePanels::from_triangles(original.triangles).unwrap();
    let b=SpherePanels::from_triangles(added.triangles.clone()).unwrap();
    let volume=|s:&SpherePanels|s.centroids().iter().zip(s.normals()).zip(s.areas())
        .map(|((p,n),a)|dot(*p,*n)*a/3.0).sum::<f64>();
    assert!((volume(&a)-volume(&b)).abs()<1e-14);
    for velocity in [-0.01,0.0,0.02] {
        let mut state=vec![0.0;2*cavity.total_modes()];state[2*port.coordinate+1]=velocity;
        let neck=cavity.neck_observation(&state,0).unwrap();
        assert!((velocity*flux(&added,2)-neck.volume_flow_m3_s).abs()<1e-15);
    }
    assert!(cavity.neck_radiation_port(1).is_err());
}

#[test]
fn circular_patch_area_refines_without_changing_the_outer_edge_vertices() {
    let outer=[[-0.04,-0.06,0.0],[0.04,-0.06,0.0],[0.04,0.0,0.0],
        [0.04,0.06,0.0],[-0.04,0.06,0.0],[-0.04,0.0,0.0]];
    let mut previous=f64::INFINITY;
    for segments in [16,32,64] {
        let p=partition(&outer,[0.0;3],[1.0,0.0,0.0],[0.0,1.0,0.0],0.01,segments).unwrap();
        let error=PI*0.01_f64.powi(2)-p.area_m2;
        assert!(error>0.0 && error<0.3*previous);previous=error;
        assert!((p.triangles.iter().map(|t|area(*t)).sum::<f64>()-0.08*0.12).abs()<1e-14);
        assert_eq!(perimeter(&p.triangles).unwrap().len(),outer.len());
    }
}

#[test]
fn aperture_is_a_real_bem_source_with_linear_superposition_and_compact_monopole_limit() {
    let b=vented();let medium=Medium::air();let w=TAU*40.0;
    let surface=SpherePanels::from_triangles(b.triangles.clone()).unwrap();
    // Unit generalized acceleration produces volume acceleration equal to the
    // admitted neck weight. No independent pressure or full-scale gain enters.
    let mut fields=acceleration_fields(&b.weights,w);
    let combined:Vec<_>=(0..surface.areas().len()).map(|i|
        fields[0][i].scale(0.3)+fields[1][i].scale(-0.2)+fields[2][i]).collect();
    fields.push(combined);
    let refs:Vec<_>=fields.iter().map(Vec::as_slice).collect();
    let solutions=solve_radiation_batch(&surface,w/medium.sound_speed,medium,&refs,Formulation::PlainCbie).unwrap();
    let pressures=solutions.iter().map(|s|exterior_pressure_at_points(&surface,s,medium,&[[5.0,0.0,0.0]]).unwrap()[0]).collect::<Vec<C64>>();
    let sum=pressures[0].scale(0.3)+pressures[1].scale(-0.2)+pressures[2];
    assert!((sum-pressures[3]).abs()<1e-8*sum.abs().max(1e-10));
    let compact=medium.density*port().volume_weight_m2_per_sqrt_kg/(4.0*PI*5.0);
    assert!(pressures[2].abs()>0.0);
    // Finite ka, finite distance and coarse constant-panel error are explicit;
    // the flux/coordinate test above is independent of this asymptotic screen.
    assert!((pressures[2].abs()-compact).abs()<0.2*compact);
}

#[test]
fn invalid_aperture_chart_or_source_address_refuses_without_silent_substitution() {
    let p=port();
    for (angle,z) in [(0.0,0.06),(PI/8.0,0.0),(PI/8.0,0.12),(f64::NAN,0.06)] {
        assert!(boundary().with_sidewall_aperture(p,angle,z,0.12).is_err());
    }
    for bad in [NeckRadiationPort {coordinate:1,..p},
        NeckRadiationPort {effective_length_m:0.012,..p},
        NeckRadiationPort {volume_weight_m2_per_sqrt_kg:0.0,..p},
        NeckRadiationPort {area_m2:f64::NAN,..p}] {
        assert!(boundary().with_sidewall_aperture(bad,PI/8.0,0.06,0.12).is_err());
    }
    let mut args=["snare-mic","64","--prescribed-vent-radiation","--cavity-modes"].map(String::from).to_vec();
    assert!(option(&mut args).unwrap());assert_eq!(args,["snare-mic","64","--cavity-modes"]);
    let mut duplicate=vec!["--prescribed-vent-radiation".into();2];let before=duplicate.clone();
    assert!(option(&mut duplicate).is_err());assert_eq!(duplicate,before);
}
