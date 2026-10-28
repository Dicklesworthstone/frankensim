use super::*;
use crate::{Mechanics, Stroke, drum_spec, mallets, shaft_playing};
use fs_exec::CancelGate;

fn drum() -> drum_spec::Spec {
    let mut d = drum_spec::Spec { radial_intervals:1, azimuths:8, band_hz:[0.0,2000.0],
        ..drum_spec::Spec::reference() };
    for h in &mut d.heads { h.damping_ratio = 0.0; }
    d
}
fn material() -> barrel::Spec {
    barrel::Spec {material:barrel::Material::Isotropic {young_pa:2e7,poisson:0.3},density_kg_m3:1000.0,damping_ratio:0.0,
        axial_intervals:2,band_hz:[0.0,100000.0]}
}

#[test]
fn sidewall_selection_requires_one_finite_station_on_an_admitted_drum() {
    let mut args = vec!["drum".into(),"--barrel-strike".into(),"0.3".into(),"0.08".into()];
    let selected = option(&mut args).unwrap().unwrap();
    assert_eq!(args, ["drum"]);
    assert_eq!(selected, Spec {azimuth_rad:0.3,height_below_batter_m:0.08});
    assert!(option(&mut args).unwrap().is_none());
    for family in ["drum","drum-modal","drum-stretch","snare","snare-off"] {
        for suffix in ["","-wav","-mic"] { admit_command(Some(&selected),&format!("{family}{suffix}")).unwrap(); }
    }
    for command in ["hihat","splash-wav","unknown"] { assert!(admit_command(Some(&selected),command).is_err()); }
    for raw in [vec!["--barrel-strike"], vec!["--barrel-strike","0.3"],
        vec!["--barrel-strike","NaN","0.08"], vec!["--barrel-strike","7","0.08"],
        vec!["--barrel-strike","0","-0.08"], vec!["--barrel-strike","0","inf"],
        vec!["--barrel-strike","0","0.08","--barrel-strike","1","0.09"]] {
        let mut args:Vec<String> = raw.into_iter().map(str::to_owned).collect();
        let before = args.clone(); assert!(option(&mut args).is_err()); assert_eq!(args,before);
    }
    let d = drum(); let b = material().prepare(&d,1e-6,false).unwrap();
    for height in [0.0,-0.01,d.depth_m,d.depth_m+0.01,f64::NAN] {
        assert!(Spec {azimuth_rad:0.3,height_below_batter_m:height}.port(&b).is_err());
    }
    let material=material(); let empty=mallets::Selection::default(); let shafts=shaft_playing::Selection::default();
    let build=|barrel,position,mallets:&mallets::Selection| crate::drum_with_barrel(1,1e-6,false,false,
        None,false,Stroke {speed_m_s:0.0,position_m:position},false,None,Some(d.clone()),None,&[],0.0,None,false,None,
        mallets,&shafts,barrel,Some(selected));
    assert!(build(None,None,&empty).is_err(),"a rigid sidewall has no moving contact port");
    assert!(build(Some(&material),Some([0.06,0.0]),&empty).is_err(),"one hand cannot silently choose between two targets");
    let felt=mallets::Selection {first:Some(mallets::Spec::parse(include_str!("estimated-felt-mallet.fsmallet")).unwrap()),second:None};
    assert!(build(Some(&material),None,&felt).is_err(),"a flat-head felt footprint cannot be reused on an undeclared curved sidewall");
}

#[test]
fn located_contact_matches_actual_radiating_skin_and_keeps_the_offset_moment() {
    let d = drum(); let b = material().prepare(&d,1e-6,false).unwrap();
    let mut offset_effect = 0.0_f64;
    for facet in [0,3,16,21] {
        let points = b.triangles[facet].map(|i| b.outer_positions[i]);
        let centre:[f64;3] = core::array::from_fn(|i| points.iter().map(|p|p[i]/3.0).sum());
        let spec = Spec {azimuth_rad:centre[1].atan2(centre[0]),height_below_batter_m:0.5*d.depth_m-centre[2]};
        let hit = spec.locate(&b).unwrap(); let port = spec.port(&b).unwrap();
        assert_eq!(hit.facet,facet);
        for (actual,expected) in port.position_m.iter().zip(centre) { assert!((actual-expected).abs()<1e-13); }
        let (normal,_) = barrel::normal_area(points).unwrap();
        assert!(dot(normal,[centre[0],centre[1],0.0])>0.0);
        let midsurface = b.reduction.point_port(facet,[1.0/3.0;3],normal).unwrap();
        for (k,w) in port.weights.iter().enumerate() {
            assert!((w-b.outer_weights[k][facet]).abs()<1e-10*(1.0+w.abs()));
            offset_effect = offset_effect.max((w-midsurface[k]).abs());
        }
        // An inward force does negative work against outward skin velocity.
        // The contact transpose must retain that sign for every shell mode.
        let velocities:Vec<_> = (0..b.mode_count()).map(|k|0.001*(k+1) as f64).collect();
        let outward:f64 = port.weights.iter().zip(&velocities).map(|(w,v)|w*v).sum();
        let generalized:f64 = port.weights.iter().zip(&velocities).map(|(w,v)|-2.3*w*v).sum();
        assert!((generalized+2.3*outward).abs()<1e-12);
    }
    assert!(offset_effect>1e-10,"the finite skin force must retain shell rotations, not just midsurface translation");
    // Exact interior ring vertex: one physical contact, deterministic facet.
    let node = b.azimuths; let p = b.outer_positions[node];
    let spec = Spec {azimuth_rad:p[1].atan2(p[0]),height_below_batter_m:0.5*d.depth_m-p[2]};
    let expected = b.triangles.iter().position(|t|t.contains(&node)).unwrap();
    assert_eq!(spec.locate(&b).unwrap().facet,expected);
    assert_eq!(spec.port(&b).unwrap().weights,spec.port(&b).unwrap().weights);
}

#[test]
fn an_exterior_strike_drives_wall_pressure_and_head_motion_on_the_shared_clock() {
    use fs_couple::render::plate::impact::{ImpactSubstepConfig, striker::{RadiusStation, flexible::FlexibleStriker}};
    use fs_plate::shell::stiffened::beam::RoundBeamSpec;
    let d = drum(); let material = material(); let dt = 1e-6; let steps = 512;
    let b = material.prepare(&d,dt,false).unwrap();
    let selected = Spec {azimuth_rad:0.3,height_below_batter_m:0.4*d.depth_m};
    let port = selected.port(&b).unwrap();
    let (_,heads) = d.prepare(dt,false).unwrap(); let first = 1+heads.iter().map(Vec::len).sum::<usize>();
    let mut shafts = shaft_playing::Selection::default();
    for distributed in [false,true] {
        if distributed {
            shafts.first = Some(FlexibleStriker::new(&[(0.0,0.005),(0.4,0.005)].map(|(position_m,radius_m)|
                RadiusStation {position_m,radius_m}),RoundBeamSpec {young_pa:12e9,density_kg_m3:800.0,
                pivot_m:0.1,contact_m:0.39,hand_m:0.16,subdivisions:8,maximum_hz:3000.0,maximum_modes:17},0.001).unwrap());
        }
        let mut e = crate::drum_with_barrel(steps,dt,false,false,None,false,
            Stroke {speed_m_s:2.0,position_m:None},distributed,None,Some(d.clone()),None,&[],0.0,None,false,None,
            &mallets::Selection::default(),&shafts,Some(&material),Some(selected)).unwrap();
        let initial = match &e.system { Mechanics::Reference(s)=>s.stored_energy_j(), _=>panic!("joint contact owner") };
        let (q,v) = shaft_playing::tip_motion(&e,0).unwrap();
        assert!((q+0.0002).abs()<1e-14 && (v-2.0).abs()<1e-14,"selection must preserve the physical launch");
        e.system = e.system.into_analytic_nonlinear().unwrap().with_impact_substeps(
            ImpactSubstepConfig {max_depth:8,max_attempts:511}).unwrap();
        let gate = CancelGate::new_clock_free();
        let (mut loss,mut inward,mut pressure,mut head_speed,mut source_speed) = (0.0,0.0_f64,0.0_f64,0.0_f64,0.0_f64);
        for tick in 0..steps {
            if tick==128 {
                let before = e.system.state().to_vec(); let mut invalid = e.force.clone(); invalid[0]=f64::NAN;
                assert!(e.system.step(&invalid,&gate).is_err()); assert_eq!(e.system.state(),before);
                let stop=CancelGate::new_clock_free(); stop.request();
                assert!(e.system.step(&e.force,&stop).is_err()); assert_eq!(e.system.state(),before);
            }
            let frame = e.system.step(&e.force,&gate).unwrap(); loss += frame.dissipated_energy_j;
            assert!(frame.balance_residual_j.abs()<1e-7);
            assert!((frame.stored_energy_j+loss-initial).abs()<1e-6);
            let state=e.system.state();
            let outward:f64 = port.weights.iter().enumerate().map(|(k,w)|w*state[2*(first+k)+1]).sum();
            inward=inward.max(-outward);
            let p = if let Some(air)=&e.air {air.uniform_pressure(state).unwrap()}
                else {crate::cavity_pressure(e.pressure.as_ref().unwrap(),state)};
            pressure=pressure.max(p.abs());
            head_speed=head_speed.max((1..first).map(|i|state[2*i+1].abs()).fold(0.0,f64::max));
            for facet in 0..b.triangles.len() {
                let v:f64=(0..b.mode_count()).map(|k|b.outer_weights[k][facet]*state[2*(first+k)+1]).sum();
                source_speed=source_speed.max(v.abs());
            }
            if tick<64 {assert!(state[2..2*(first+b.mode_count())].iter().all(|v|v.abs()<1e-18),
                "there is no head contact or moving barrel before the exterior gap closes");}
        }
        assert!(inward>1e-6 && pressure>1e-5 && head_speed>1e-12 && source_speed>1e-6,
            "an inward sidewall strike must drive the actual shell, enclosed pressure, heads and radiation rows");
        assert!((shaft_playing::tip_motion(&e,0).unwrap().1-2.0).abs()>1e-4,
            "the stick must receive the wall's reciprocal contact reaction");
        if let Some(p)=&e.flexible_sticks[0] {assert!(p.observe(e.system.state()).unwrap().flexural_energy_j>1e-12);}
    }
}
