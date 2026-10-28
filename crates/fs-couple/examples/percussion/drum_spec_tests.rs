use super::*;
use crate::{Stroke, drum_with_spec, mechanics::Mechanics};
use fs_exec::CancelGate;

/// The shipped default instrument environment (ruling 6205), which the
/// built-in reference also reads.
const AMBIENT: &str = "293.15 101325 0.65";

fn text(s: Spec) -> String {
    let mut out=format!("{HEADER}\ngeometry,{:.17e},{:.17e},{:.17e}\n",
        s.radius_m,s.depth_m,s.outer_radius_m);
    for ((name,h),loss) in ["batter","resonant"].into_iter().zip(s.heads).zip(s.head_loss) {
        let loss=if loss==HeadLoss::Relaxation {RELAXATION_HEAD_LOSS} else {STATED_UNKNOWN_HEAD_LOSS};
        out.push_str(&format!("head,{name},{:.17e},{:.17e},{:.17e},{:.17e},{:.17e},{loss}\n",
            h.thickness_m,h.young_pa,h.poisson,h.density_kg_m3,h.tension_n_m));
    }
    out.push_str(&format!("mesh,{},{}\nband_hz,{:.17e},{:.17e}\n",
        s.radial_intervals,s.azimuths,s.band_hz[0],s.band_hz[1]));
    out.push_str(&format!("ambient,{AMBIENT}\n"));
    out
}
/// What `text` declares: the spec with its air read from the ambient row, and
/// each head's loss as its file token. A file cannot carry a fixture's ratio.
fn declared(s: Spec) -> Spec {
    let file_loss=|loss| if loss==HeadLoss::Relaxation {loss} else {HeadLoss::StatedUnknown};
    Spec {ambient:crate::ambient::Ambient::described(AMBIENT).unwrap(),
        heads:s.heads.map(|h| Head {damping_ratio:0.0,..h}),head_loss:s.head_loss.map(file_loss),..s}
}
fn small() -> Spec { Spec {radial_intervals:2,azimuths:16,..Spec::reference()} }
fn stroke() -> Stroke {Stroke {speed_m_s:4.0,position_m:Some([0.06,0.01])}}

#[test]
fn complete_materials_and_mesh_round_trip_without_defaults() {
    let mut spec=small();
    spec.heads[0].tension_n_m=2456.7;spec.heads[1].young_pa=3.2e9;
    spec.heads[1].damping_ratio=0.0;spec.depth_m=0.20;
    let file=text(spec.clone());
    assert_eq!(Spec::read(&file).unwrap(),declared(spec.clone()));
    let mut records:Vec<_>=file.lines().skip(1).collect();records.reverse();
    assert_eq!(Spec::read(&format!("# inputs, not certified measurements\n{HEADER}\n{}\n",records.join("\n"))).unwrap(),declared(spec.clone()));
    let film=spec.head(1).unwrap();
    assert_eq!(film.spec.young_pa,spec.heads[1].young_pa);
    assert_eq!(film.spec.thickness_m,spec.heads[1].thickness_m);
    assert_eq!(film.spec.density_kg_m3,spec.heads[1].density_kg_m3);
    assert_eq!(film.spec.poisson,spec.heads[1].poisson);
    assert_eq!(film.spec.tension_n_m,spec.heads[1].tension_n_m);
    assert_eq!(film.spec.radius_m,spec.radius_m);
    assert_eq!(film.mesh.nodes.len(),1+spec.radial_intervals*spec.azimuths);
}

#[test]
fn incomplete_duplicate_nonfinite_and_nonphysical_inputs_refuse() {
    let valid=text(small());
    for (i,_) in valid.lines().enumerate() {
        let missing=valid.lines().enumerate().filter(|(j,_)|*j!=i)
            .map(|(_,s)|s).collect::<Vec<_>>().join("\n");
        assert!(Spec::read(&missing).is_err(),"accepted missing record {i}");
    }
    for row in valid.lines().skip(1) {
        assert!(Spec::read(&format!("{valid}{row}\n")).is_err());
    }
    for bad in [String::new(),valid.replace(HEADER,"wrong-version"),
        format!("{valid}pitch,440\n"),format!("{valid}{}","#".repeat(MAX_BYTES)),
        valid.replace("mesh,2,16","mesh,2.5,16"),valid.replace("mesh,2,16","mesh,0,16"),
        valid.replace("mesh,2,16","mesh,33,16"),valid.replace("mesh,2,16","mesh,2,129"),
        valid.replace("head,batter,","head,unknown,"),
        valid.replacen("geometry,","geometry,NaN,",1)] {
        assert!(Spec::read(&bad).is_err(),"accepted malformed declaration");
    }
    for bad in [f64::NAN,f64::INFINITY,-1.0,0.0] {
        let mut s=small();s.heads[0].tension_n_m=bad;
        assert!(Spec::read(&text(s)).is_err());
        let mut s=small();s.radius_m=bad;
        assert!(Spec::read(&text(s)).is_err());
    }
    let mut s=small();s.outer_radius_m=s.radius_m;assert!(s.validate().is_err());
    let mut s=small();s.heads[0].damping_ratio=-0.01;assert!(s.validate().is_err());
    let mut s=small();s.heads[1].young_pa=0.0;assert!(s.validate().is_err());
    let mut s=small();s.heads[1].poisson=1.0;assert!(s.validate().is_err());
}

#[test]
fn a_head_states_where_its_loss_comes_from_and_an_authored_ratio_refuses() {
    let valid=text(small());
    let spec=Spec::read(&valid).unwrap();
    assert_eq!(spec.head_loss,[HeadLoss::StatedUnknown;2]);
    assert_eq!(spec.heads.map(|h|h.damping_ratio),[0.0;2]);
    assert!(spec.admit_head_loss(false).is_ok());
    // Every number refuses, zero included: zero drag is stated, not authored.
    for ratio in ["0.001","0","1e-3"] {
        let bad=valid.replacen(",stated-unknown",&format!(",{ratio}"),1);
        let refusal=Spec::read(&bad).unwrap_err().to_string();
        assert!(refusal.contains("has no source"),"{ratio}: {refusal}");
    }
    // No file reaches the library fixture: neither its name nor its ratio is a token.
    for token in ["fitted","fixture","Fixture","0.001"] {
        assert!(Spec::read(&valid.replacen(",stated-unknown",&format!(",{token}"),1)).is_err(),"{token}");
    }
    let relaxing=Spec::read(&valid.replace(",stated-unknown",",relaxation")).unwrap();
    assert_eq!(relaxing.head_loss,[HeadLoss::Relaxation;2]);
    assert!(relaxing.admit_head_loss(false).is_err());
    assert!(relaxing.admit_head_loss(true).is_ok());
    // Only the library fixture carries a ratio.
    let mut s=small();s.head_loss[0]=HeadLoss::StatedUnknown;assert!(s.validate().is_err());
    s.heads[0].damping_ratio=0.0;assert!(s.validate().is_ok());
    assert_eq!(HeadLoss::StatedUnknown.receipt("batter",0.0),
        "head batter film and edge loss: STATED UNKNOWN, entered as 0; its decays are upper bounds");
    // The shipped specimens state their losses.
    for (file,loss) in [(include_str!("estimated_drum.fsd"),HeadLoss::StatedUnknown),
        (include_str!("estimated-varied-tension.fsd"),HeadLoss::StatedUnknown),
        (include_str!("illustrative-ten-lug.fsd"),HeadLoss::StatedUnknown),
        (include_str!("estimated-relaxing-drum.fsd"),HeadLoss::Relaxation)] {
        assert_eq!(Spec::read(file).unwrap().head_loss,[loss;2]);
    }
}

#[test]
fn option_composes_with_playing_and_existing_physical_images() {
    let mut args:Vec<_>=["drum-stretch-mic","128","--drum-spec","instrument.fsd",
        "--cavity-modes","--prepared-nonlinear","--strike-speed-m-s","4"]
        .into_iter().map(String::from).collect();
    assert_eq!(option(&mut args).unwrap().as_deref(),Some("instrument.fsd"));
    assert!(crate::cavity::option(&mut args).unwrap());
    assert!(crate::mechanics::prepared_option(&mut args).unwrap());
    let (args,stroke)=crate::playing::parse(args).unwrap();
    assert_eq!(args,["drum-stretch-mic","128"]);assert_eq!(stroke.speed_m_s,4.0);
    for command in ["drum","drum-mic","drum-stretch-wav","drum-modal","snare","snare-off-mic"] {
        assert!(admit_command(Some("instrument.fsd"),command).is_ok());
    }
    for command in ["splash","splash-mic","unknown"] {
        assert!(admit_command(Some("instrument.fsd"),command).is_err());
    }
    for raw in ["--drum-spec", "--drum-spec --cavity-modes", "--drum-spec a --drum-spec b"] {
        let mut args:Vec<_>=raw.split_whitespace().map(String::from).collect();
        let saved=args.clone();assert!(option(&mut args).is_err());assert_eq!(args,saved);
    }
}

#[test]
fn clocks_acoustic_band_and_fixed_snare_geometry_are_not_silently_changed() {
    let spec=small();
    for dt in [0.0,-1.0,f64::NAN,0.001] {assert!(spec.admit_clock(dt,false).is_err());}
    let wide=Spec {band_hz:[20.0,2000.0],..spec.clone()};
    assert!(wide.admit_clock(2e-6,false).is_ok());
    assert!(wide.admit_clock(crate::acoustics::MECHANICAL_DT,true).is_err());
    let wires=crate::snare::SnareSet::reference(false);
    assert!(spec.admit_snare(wires).is_ok());
    assert!(Spec {radius_m:0.08,outer_radius_m:0.09,..spec}.admit_snare(wires).is_err());
}

/// A file cannot carry the library fixture's 0.001 ratio, so the original here
/// is the library reference with the loss its file declares (stated unknown).
#[test]
fn supplied_reference_keeps_the_original_contact_air_and_accepted_trajectory() {
    let supplied=Spec::read(&text(Spec::reference())).unwrap();
    assert_eq!(supplied,declared(Spec::reference()));
    let mut original=drum_with_spec(64,2e-6,false,false,None,false,stroke(),false,None,
        Some(declared(Spec::reference()))).unwrap();
    let mut imported=drum_with_spec(64,2e-6,false,false,None,false,stroke(),false,None,Some(supplied)).unwrap();
    assert_eq!(original.observer_a,imported.observer_a);
    assert_eq!(original.observer_b,imported.observer_b);
    assert_eq!(original.pressure.as_ref().unwrap().volume_m3,imported.pressure.as_ref().unwrap().volume_m3);
    assert_eq!(original.pressure.as_ref().unwrap().areas,imported.pressure.as_ref().unwrap().areas);
    assert_eq!(original.system.state(),imported.system.state());
    let gate=CancelGate::new_clock_free();
    for _ in 0..64 {
        let a=original.system.step(&original.force,&gate).unwrap();
        let b=imported.system.step(&imported.force,&gate).unwrap();
        assert_eq!(a.stored_energy_j.to_bits(),b.stored_energy_j.to_bits());
        assert_eq!(original.system.state(),imported.system.state());
    }
}

#[test]
fn tension_and_elasticity_change_the_actual_pencil_and_resonances() {
    let a=small();let mut b=a.clone();
    for h in &mut b.heads {h.young_pa*=4.0;h.tension_n_m*=4.0;}
    b.band_hz=a.band_hz.map(|f|2.0*f);
    let (fa,ma)=a.prepare(2e-6,false).unwrap();
    let (fb,mb)=b.prepare(2e-6,false).unwrap();
    for i in 0..2 {
        assert_eq!(fa[i].mass_kg,fb[i].mass_kg);
        assert_eq!(fa[i].mesh.nodes,fb[i].mesh.nodes);
        assert_eq!(ma[i].len(),mb[i].len());
        let n=fa[i].model.free;
        for r in 0..n {for c in 0..n {
            assert_eq!(fa[i].model.m.get(r,c),fb[i].model.m.get(r,c));
            let expected=4.0*fa[i].model.k.get(r,c);
            assert!((fb[i].model.k.get(r,c)-expected).abs()<1e-10*(1.0+expected.abs()));
        }}
        for (x,y) in ma[i].iter().zip(&mb[i]) {
            assert!((y.lambda/(4.0*x.lambda)-1.0).abs()<1e-5,"physical stiffness scaling must double frequency");
        }
    }
    let mut tuned=a;for h in &mut tuned.heads {h.tension_n_m*=2.0;}
    let (_,mt)=tuned.prepare(2e-6,false).unwrap();
    let lowest=|m:&[ModePair]|m.iter().map(|p|p.lambda).fold(f64::INFINITY,f64::min);
    for i in 0..2 {assert!(lowest(&mt[i])>lowest(&ma[i]));}
}

#[test]
fn depth_changes_real_air_feedback_without_retuning_the_heads() {
    let a=small();let b=Spec {depth_m:1.5*a.depth_m,..a.clone()};
    let build=|s|drum_with_spec(128,2e-6,false,false,None,false,stroke(),false,None,Some(s)).unwrap();
    let mut shallow=build(a);let mut deep=build(b);
    assert_eq!(shallow.observer_a,deep.observer_a);
    assert_eq!(shallow.system.state(),deep.system.state());
    let va=shallow.pressure.as_ref().unwrap();let vb=deep.pressure.as_ref().unwrap();
    assert_eq!(va.areas,vb.areas);
    assert!((vb.volume_m3/va.volume_m3-1.5).abs()<1e-14);
    let gate=CancelGate::new_clock_free();let mut changed=0.0_f64;
    for _ in 0..128 {
        for e in [&mut shallow,&mut deep] {
            let f=e.system.step(&e.force,&gate).unwrap();
            assert!(f.stored_energy_j.is_finite());assert!(f.balance_residual_j.abs()<2e-7);
        }
        for (x,y) in shallow.system.state().iter().zip(deep.system.state()) {changed=changed.max((x-y).abs());}
    }
    assert!(changed>1e-14,"changed cavity volume must feed back into the accepted head/contact motion");
}

#[test]
fn larger_supplied_heads_reach_their_matching_exterior_not_the_stock_shell() {
    let s=Spec {radius_m:0.20,outer_radius_m:0.21,depth_m:0.14,..small()};
    let e=drum_with_spec(2,crate::acoustics::MECHANICAL_DT,true,false,None,false,
        stroke(),false,None,Some(s.clone())).unwrap();
    assert!(e.acoustics.is_some());
    assert_eq!(e.pressure.as_ref().unwrap().volume_m3,s.volume_m3());
    // The old 0.1778m outer surface cannot enclose these 0.20m heads.
    let (films,modes)=s.prepare(2e-6,false).unwrap();
    assert!(crate::acoustics::Boundary::drum(&films,&modes,s.depth_m,0.1778).is_err());
}

#[test]
fn supplied_drum_retains_stretching_distributed_air_and_prepared_nonlinearity() {
    let s=Spec {depth_m:0.18,..small()};
    let mut e=drum_with_spec(64,2e-6,false,false,None,true,stroke(),true,None,Some(s)).unwrap();
    assert!(e.air.is_some());
    e.system=e.system.into_prepared_nonlinear().unwrap();
    assert!(matches!(&e.system,Mechanics::Nonlinear(_)));
    let gate=CancelGate::new_clock_free();let mut stretching=0.0_f64;
    for _ in 0..64 {
        let f=e.system.step(&e.force,&gate).unwrap();assert!(f.balance_residual_j.abs()<2e-7);
        stretching=stretching.max(e.system.membrane_observation(1).unwrap().stretching_energy_j);
        let (a,b)=e.air.as_ref().unwrap().points(e.system.state()).unwrap();
        assert!(a.is_finite() && b.is_finite());
    }
    assert!(stretching>0.0);
}


#[test]
fn multi_lug_harmonics_parse_into_equilibrated_prestress_and_change_modes_not_mass() {
    let base=small();
    let declared=format!("{}tension_harmonic,batter,8,120,-35\ntension_harmonic,resonant,10,60,20\n",text(base.clone()));
    let tuned=Spec::read(&declared).unwrap();
    assert_eq!(tuned.tension_variation[0].harmonics_n_m[6],[120.0,-35.0]);
    assert_eq!(tuned.tension_variation[1].harmonics_n_m[8],[60.0,20.0]);
    let (plain_f,plain_m)=base.prepare(2e-6,false).unwrap();
    let (tuned_f,tuned_m)=tuned.prepare(2e-6,false).unwrap();
    for i in 0..2 {
        assert_eq!(plain_f[i].model.m,tuned_f[i].model.m);
        assert_eq!(plain_f[i].mesh.nodes,tuned_f[i].mesh.nodes);
        assert_eq!(plain_f[i].model.dof_map,tuned_f[i].model.dof_map);
        assert_eq!(plain_m[i].len(),tuned_m[i].len());
        assert!(plain_m[i].iter().zip(&tuned_m[i]).any(|(a,b)|(a.lambda-b.lambda).abs()>1e-8*a.lambda.max(1.0)));
    }
}

#[test]
fn malformed_duplicate_or_overstrong_lug_harmonics_refuse_without_clipping() {
    let valid=text(small());
    for row in [
        "tension_harmonic,batter,1,10,0",
        "tension_harmonic,batter,13,10,0",
        "tension_harmonic,batter,8,NaN,0",
        "tension_harmonic,unknown,8,10,0",
    ] { assert!(Spec::read(&format!("{valid}{row}\n")).is_err(),"{row}"); }
    assert!(Spec::read(&format!("{valid}tension_harmonic,batter,8,10,0\ntension_harmonic,batter,8,20,0\n")).is_err());
    // Harmonic bound exceeds the 1500 N/m resonant base, so physical head
    // construction refuses rather than scaling or clipping the requested map.
    let s=Spec::read(&format!("{valid}tension_harmonic,resonant,8,1600,0\n")).unwrap();
    assert!(s.head(1).is_err());
}


#[test]
fn affine_and_harmonic_records_compose_in_either_order() {
    let base = small();
    let file = format!(
        "{}tension_harmonic,batter,2,45,-12\n\
         tension_variation,resonant,20,40,-2,-1,2,3,-4\n\
         tension_harmonic,batter,8,120,-35\n\
         tension_variation,batter,120,70,10,1,-2,3,4\n\
         tension_harmonic,resonant,12,10,-4\n",
        text(base.clone())
    );
    let parsed = Spec::read(&file).unwrap();
    let mut records: Vec<_> = file.lines().skip(1).collect();
    records.reverse();
    let reversed = Spec::read(&format!("{HEADER}\n{}\n", records.join("\n"))).unwrap();
    let mut expected = declared(base);
    expected.tension_variation[0].constant_n_m = [120.0, 70.0, 10.0];
    expected.tension_variation[0].gradient_n_m2 = [1.0, -2.0, 3.0, 4.0];
    expected.tension_variation[0].harmonics_n_m[0] = [45.0, -12.0];
    expected.tension_variation[0].harmonics_n_m[6] = [120.0, -35.0];
    expected.tension_variation[1].constant_n_m = [20.0, 40.0, -2.0];
    expected.tension_variation[1].gradient_n_m2 = [-1.0, 2.0, 3.0, -4.0];
    expected.tension_variation[1].harmonics_n_m[10] = [10.0, -4.0];
    assert_eq!(parsed, expected);
    assert_eq!(reversed, expected);
    // Check the actual FE pencils, not only the parser's stored coefficients.
    for head in 0..2 {
        let declared = parsed.head(head).unwrap();
        let physical = expected.head(head).unwrap();
        assert_eq!(declared.mesh.nodes, physical.mesh.nodes);
        assert_eq!(declared.model.dof_map, physical.model.dof_map);
        assert_eq!(declared.model.m, physical.model.m);
        assert_eq!(declared.model.k, physical.model.k);
    }
}

#[test]
fn tension_duplicates_are_detected_by_record_kind_even_for_zero_values() {
    let valid = text(small());
    for head in ["batter", "resonant"] {
        let affine = format!("tension_variation,{head},0,0,0,0,0,0,0\n");
        let harmonic = format!("tension_harmonic,{head},8,0,0\n");
        for rows in [
            format!("{affine}{affine}"),
            format!("{harmonic}{affine}{affine}"),
            format!("{affine}{harmonic}{affine}"),
        ] {
            let error = Spec::read(&format!("{valid}{rows}")).unwrap_err();
            assert!(error.to_string().contains("duplicate head tension variation"));
        }
        for rows in [
            format!("{harmonic}{affine}{harmonic}"),
            format!("{affine}{harmonic}{harmonic}"),
        ] {
            let error = Spec::read(&format!("{valid}{rows}")).unwrap_err();
            assert!(error.to_string().contains("duplicate head tension harmonic"));
        }
    }
}

#[test]
fn legacy_affine_and_explicit_zero_harmonics_preserve_the_head_pencil() {
    let file = format!(
        "{}tension_variation,batter,120,70,10,1,-2,3,4\n",
        text(small())
    );
    let affine = Spec::read(&file).unwrap();
    let mut expected = declared(small());
    expected.tension_variation[0].constant_n_m = [120.0, 70.0, 10.0];
    expected.tension_variation[0].gradient_n_m2 = [1.0, -2.0, 3.0, 4.0];
    assert_eq!(affine, expected);
    let explicit_zero = Spec::read(&format!(
        "{file}tension_harmonic,batter,2,-0,0\n\
         tension_harmonic,resonant,12,0,-0\n"
    )).unwrap();
    assert_eq!(explicit_zero, expected);
    for head in 0..2 {
        let before = expected.head(head).unwrap();
        let after = explicit_zero.head(head).unwrap();
        assert_eq!(before.model.free, after.model.free);
        for row in 0..before.model.free {
            for col in 0..before.model.free {
                assert_eq!(before.model.k.get(row, col).to_bits(), after.model.k.get(row, col).to_bits());
                assert_eq!(before.model.m.get(row, col).to_bits(), after.model.m.get(row, col).to_bits());
            }
        }
    }
}
