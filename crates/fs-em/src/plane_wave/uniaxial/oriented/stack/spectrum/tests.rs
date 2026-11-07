use super::*;
use super::super::{ResolvedLayer,C0};

const MU:f64=1.256_637_062_12e-6;
fn fixture()->Frequency {
    Frequency {omega:C0,ni:1.7,nt:1.3,mui:MU,mut_:MU,layers:vec![
        ResolvedLayer {no:1.4,ne:2.1,mu:MU,axis:[0.3,0.7,0.4],thickness:0.6,entrance:0.0,exit:0.6},
        ResolvedLayer {no:1.6,ne:2.2,mu:MU,axis:[0.6,-0.2,0.7],thickness:0.4,entrance:0.6,exit:1.0},
    ]}
}
fn rotate(v:[f64;2],c:f64,s:f64)->[f64;2] {[c*v[0]-s*v[1],s*v[0]+c*v[1]]}
fn near(a:C64,b:C64,t:f64) {assert!((a-b).abs()<t,"{a:?} != {b:?}");}

#[test]
fn multiple_layer_axes_and_all_mode_frames_rotate_together() {
    let mut data=fixture();
    let input=[
        UniaxialIncidentMode {tangential_phase_index:0.5,incidence_plane_axis:[1.0,0.0],electric_v_m:[C64::ONE,C64::new(0.0,0.3)]},
        UniaxialIncidentMode {tangential_phase_index:-0.4,incidence_plane_axis:[0.6,0.8],electric_v_m:[C64::new(0.2,-0.1),C64::ONE]},
    ];
    let modes=prepare_modes(&data,&input,&mut ||false).unwrap();
    let (s,c)=0.73_f64.sin_cos();
    for layer in &mut data.layers {
        let a=rotate([layer.axis[0],layer.axis[1]],c,s);
        layer.axis=[a[0],a[1],layer.axis[2]];
    }
    let rotated=input.map(|mode| UniaxialIncidentMode {
        incidence_plane_axis:rotate(mode.incidence_plane_axis,c,s),
        electric_v_m:[mode.electric_v_m[0].scale(c)-mode.electric_v_m[1].scale(s),
            mode.electric_v_m[0].scale(s)+mode.electric_v_m[1].scale(c)],..mode
    });
    let rotated=prepare_modes(&data,&rotated,&mut ||false).unwrap();
    for (region,z) in [(StackRegion::TotalIncidentSide,-0.2),(StackRegion::Layer(0),0.3),
        (StackRegion::Layer(1),0.8),(StackRegion::Transmitted,1.2)] {
        let point=[0.31,-0.27,z];let xy=rotate([point[0],point[1]],c,s);
        let a=field_at(&modes,region,point,&mut ||false).unwrap();
        let b=field_at(&rotated,region,[xy[0],xy[1],z],&mut ||false).unwrap();
        for (old,new) in [(a.electric_v_m,b.electric_v_m),(a.magnetic_a_m,b.magnetic_a_m)] {
            near(new[0],old[0].scale(c)-old[1].scale(s),4e-14);
            near(new[1],old[0].scale(s)+old[1].scale(c),4e-14);
            near(new[2],old[2],4e-14);
        }
    }
}

#[test]
fn matched_multilayer_source_retains_bright_dark_interference_and_longitudinal_h() {
    let mut data=fixture();data.ni=1.0;data.nt=1.0;
    for layer in &mut data.layers {layer.no=1.0;layer.ne=1.0;}
    let input=[-0.3,0.3].map(|p| UniaxialIncidentMode {
        tangential_phase_index:p,incidence_plane_axis:[1.0,0.0],electric_v_m:[C64::ZERO,C64::ONE],
    });
    let modes=prepare_modes(&data,&input,&mut ||false).unwrap();
    for (region,z) in [(StackRegion::Layer(0),0.2),(StackRegion::Layer(1),0.8),(StackRegion::Transmitted,1.5)] {
        let bright=field_at(&modes,region,[0.0,0.0,z],&mut ||false).unwrap();
        let dark=field_at(&modes,region,[core::f64::consts::PI/0.6,0.0,z],&mut ||false).unwrap();
        assert!((bright.electric_v_m[1].abs()-2.0).abs()<2e-14);
        assert!(dark.electric_v_m[1].abs()<2e-14);
        assert!(dark.magnetic_a_m[2].abs()>1e-4);
        assert!(dark.poynting_w_m2[2].abs()<1e-28);
        assert!((bright.poynting_w_m2[2]/(2.0*0.91_f64.sqrt()/(MU*C0))-1.0).abs()<2e-14);
    }
}

#[test]
fn preparation_caps_and_mid_sweep_cancellation_leave_no_partial_spectrum() {
    let data=fixture();
    let input=[UniaxialIncidentMode {tangential_phase_index:0.4,incidence_plane_axis:[1.0,0.0],
        electric_v_m:[C64::ONE,C64::ZERO]}];
    assert!(admit(2,1,UniaxialSpectrumLimits {maximum_modes:0,maximum_workspace_bytes:usize::MAX}).is_err());
    assert!(admit(2,1,UniaxialSpectrumLimits {maximum_modes:1,maximum_workspace_bytes:1}).is_err());
    let mut calls=0;
    assert!(matches!(prepare_modes(&data,&input,&mut ||{calls+=1;calls==8}),
        Err(crate::plane_wave::angular_spectrum::AngularSpectrumError::Cancelled)));
    let first=prepare_modes(&data,&input,&mut ||false).unwrap();
    let again=prepare_modes(&data,&input,&mut ||false).unwrap();
    assert_eq!(first[0].response,again[0].response);
    assert!(field_at(&first,StackRegion::Layer(0),[0.0,0.0,0.3],&mut ||true).is_err());
}
