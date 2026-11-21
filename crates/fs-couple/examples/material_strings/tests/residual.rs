//! `residual-flexibility exact-static-green`: the omitted bending modes in series with the bridge and every fixed
//! stop (pexfp.299, MagentaBadger 7592). The exact-string oracle is the library's; this is the input-to-PCM path.
use super::*;
use crate::stops::Profile;
use fs_couple::bowed_string::runtime::reciprocal::bank::residual;

const STOP: &str = "frankensim-string-stops-v1\nstop support 0 0.43 positive fixed 0 5000 1 0.2 10 0.001\n";

fn flexible(text: &str) -> input::Spec {
    let text = text.replace("cavity neck", "residual-flexibility exact-static-green\ncavity neck");
    input::Spec::parse(text.as_bytes()).unwrap()
}
fn pluck() -> String { include_str!("../pluck.material-strings").replace("48000 4096 0.01 128", "48000 256 0.01 128") }
// From rest, the 0.002 N port load held for 64 samples and released: no static preload.
fn struck() -> String { pluck().replace("forces preload", "forces retain").replace("force 0 0 0", "force 64 0 0") }
fn stopped(cx: &Cx<'_>, gate: &CancelGate, spec: input::Spec) -> Result<Ready> {
    spec.compile_with_profiles(cx, gate, None, Some(Profile::parse(STOP.as_bytes()).unwrap()),
        |name| Ok(pack(name, 7800.0, 200_000.0, None)))
}

#[test]
fn the_record_softens_the_bridge_through_the_omitted_modes_and_keeps_the_basis() {
    scoped(|cx, gate| {
        let flexible = flexible(&pluck());
        assert!(flexible.residual_flexibility && !spec().residual_flexibility);
        let mut stiff = compile(cx, gate, spec(), 7800.0);
        let mut soft = compile(cx, gate, flexible, 7800.0);
        let bank = soft.performance.bank();
        assert!(bank.source().bank().residual_flexibility());
        for (a, b) in stiff.performance.bank().system().components().iter().zip(bank.system().components()) {
            assert_eq!(a.modes(), b.modes(), "the residual must not retune the retained basis");
        }
        let (member, card) = (&bank.source().bank().strings()[0], &bank.cards()[0]);
        let connection = &member.connections[0];
        let condensed = residual::condense_bridge(&member.source, card, connection).unwrap();
        let g = residual::residual_compliance_m_n(card, connection.station_fraction, connection.station_fraction).unwrap();
        assert_eq!(condensed.residual_m_n, g);
        assert!(g > 0.0 && condensed.stiffness_n_m < connection.stiffness_n_m);
        // The bridge term of the stiffness bound shrinks with it.
        assert!(bank.system().linear_frequency_upper_bound_rad_s()
            < stiff.performance.bank().system().linear_frequency_upper_bound_rad_s());
        let (x, y) = (pressure(&mut stiff, 128), pressure(&mut soft, 128));
        assert!(x.iter().zip(&y).any(|(x, y)| (x-y).abs() > 1e-12));
    });
}

#[test]
fn a_fixed_stop_gains_the_bridge_cross_term_and_its_own_omitted_compliance() {
    scoped(|cx, gate| {
        let mut plain = stopped(cx, gate, input::Spec::parse(struck().as_bytes()).unwrap()).unwrap();
        let mut soft = stopped(cx, gate, flexible(&struck())).unwrap();
        assert!(plain.performance.bank().system().obstacle_residual_compliance().is_none());
        let bank = soft.performance.bank();
        let member = &bank.source().bank().strings()[0];
        let bridges = residual::residual_bridges(bank.source().bank().strings(), bank.cards()).unwrap();
        let (cross, compliance) = residual::stop_residual(&bridges, &[residual::ResidualStop {
            string: 0, station_fraction: 0.43, direction: [1.0, 0.0] }]).unwrap();
        assert_eq!(bank.system().obstacle_residual_compliance(), Some(compliance.as_slice()));
        assert!(compliance[0] > 0.0 && cross[0] != 0.0);
        assert!(bank.system().obstacle_pairs().is_empty());
        // Closure minus the cross term times the bridge link. The stop keeps its own receiver
        // (motionless for a fixed stop); the bridge's body part is a third participant, the
        // signed correction +cross*B inside the complete map left - right + correction.
        let normal = bank.system().obstacle(0).unwrap();
        let bridge = &member.connections[0].body;
        let fixed = plain.performance.bank().system().obstacle(0).unwrap();
        assert!(fixed.right.shapes.iter().all(|x| *x == 0.0));
        assert_eq!(normal.right, fixed.right);
        let offset: usize = bank.system().components()[..bank.cards().len()+bridge.component].iter()
            .map(|c| c.modes().len()).sum();
        let column = bank.system().obstacle_normal_column(0).unwrap();
        assert_eq!(column[offset..offset+bridge.shapes.len()].to_vec(),
            bridge.shapes.iter().map(|b| cross[0]*b).collect::<Vec<_>>());
        assert_eq!(normal.left.shapes.len(), fixed.left.shapes.len());
        assert!(normal.left.shapes.iter().zip(&fixed.left.shapes).any(|(a, b)| a != b));
        let (x, y) = (pressure(&mut plain, 128), pressure(&mut soft, 128));
        assert!(y.iter().any(|p| p.abs() > 1e-12));
        assert!(x.iter().zip(&y).any(|(x, y)| (x-y).abs() > 1e-12));
    });
}

#[test]
fn a_released_preload_reaches_the_stopped_card_with_the_record() {
    scoped(|cx, gate| {
        let mut plain = stopped(cx, gate, input::Spec::parse(pluck().as_bytes()).unwrap()).unwrap();
        let mut soft = stopped(cx, gate, flexible(&pluck())).unwrap();
        assert!(soft.performance.bank().system().obstacle_residual_compliance().is_some());
        assert!(soft.performance.preload_energy_j().unwrap() > 0.0);
        let (x, y) = (pressure(&mut plain, 128), pressure(&mut soft, 128));
        assert!(y.iter().any(|p| p.abs() > 1e-12));
        assert!(x.iter().zip(&y).any(|(x, y)| (x-y).abs() > 1e-12));
    });
}

// The stopped force-loaded card with longitudinal waves for its stretching: the stop stays
// condensed, the nonlinear static preload settles it in series with the omitted modes, and
// the nonlinear steps play with their storage booked. The card releases its 0.5 N load at
// sample zero; `held` keeps that same load on its port instead, so the stop stays loaded.
// `stretched_stopped` is the same card on its own geometric stretching.
fn longitudinal_stopped(residual: bool, held: bool) -> input::Spec {
    nonlinear_stopped(residual, held, "longitudinal-waves 0.2 0.01 128 100000 1e-11 1e-10 2 1 2")
}
fn stretched_stopped(residual: bool, held: bool) -> input::Spec {
    nonlinear_stopped(residual, held, "geometric-stretching 0.2 128 100000 1e-11 1e-10")
}
fn nonlinear_stopped(residual: bool, held: bool, law: &str) -> input::Spec {
    let record = if residual { "residual-flexibility exact-static-green\n" } else { "" };
    let mut text = include_str!("../nonlinear-stopped-preload.material-strings").replace("48000 4096", "48000 129")
        .replace("cylinder-gas-loading 0.001 30 4\n", "")
        .replace("geometric-stretching 0.2 128 100000 1e-11 1e-10", &format!("{record}{law}"));
    if held {
        assert!(text.contains("port 0 0.2 centreline 0.5\nforce 0 0 0\n"), "the card's release moved");
        text = text.replace("port 0 0.2 centreline 0.5\nforce 0 0 0\n", "port 0 0.2 centreline 0.5\nforce 0 0 0.5\n");
    }
    input::Spec::parse(text.as_bytes()).unwrap()
}

#[test]
fn a_stopped_preload_with_longitudinal_waves_keeps_the_condensed_stop() {
    scoped(|cx, gate| {
        let mut soft = stopped(cx, gate, longitudinal_stopped(true, false)).unwrap();
        let mut plain = stopped(cx, gate, longitudinal_stopped(false, false)).unwrap();
        let bank = soft.performance.bank();
        assert!(bank.source().bank().longitudinal_waves().is_some());
        assert!(bank.source().bank().cylinder_gas_loading().is_none());
        assert!(bank.system().obstacle_residual_compliance().unwrap()[0] > 0.0);
        assert!(plain.performance.bank().system().obstacle_residual_compliance().is_none());
        let report = *bank.system().nonlinear_preload_diagnostics().unwrap();
        assert!(report.force_residual <= report.force_tolerance);
        let total = bank.system().total_energy_j().unwrap();
        let preload = soft.performance.preload_energy_j().unwrap();
        assert!(preload > 0.0 && (preload-total).abs() <= 1e-12*total, "preload {preload} against the ledger's {total}");
        let expected = pressure(&mut soft, 37);
        assert!(expected.iter().any(|p| p.abs() > 1e-12));
        assert_eq!(expected, pressure(&mut stopped(cx, gate, longitudinal_stopped(true, false)).unwrap(), 128));
        let frame = soft.performance.bank().system().last_frame().unwrap();
        assert!(frame.nonlinear.is_some(), "the steps left the nonlinear owner");
        assert!(frame.energy_residual_j.abs() <= frame.energy_tolerance_j);
        assert!(expected.iter().zip(&pressure(&mut plain, 128)).any(|(x, y)| (x-y).abs() > 1e-12));
    });
}

// The same card with its load held, on both strain laws, with and without the record: from
// the exact static equilibrium every nonlinear step plays, and the loaded stop still stores
// R G R / 2 in the omitted modes at the last frame.
#[test]
fn a_held_stopped_preload_plays_on_both_strain_laws() {
    scoped(|cx, gate| {
        for (law, card) in [("longitudinal-waves", longitudinal_stopped as fn(bool, bool) -> input::Spec),
            ("geometric-stretching", stretched_stopped)] {
            for residual in [true, false] {
                let case = format!("{law}, residual {}", if residual { "on" } else { "off" });
                let mut held = stopped(cx, gate, card(residual, true))
                    .unwrap_or_else(|error| panic!("{case}: the held card refused: {error:?}"));
                assert_eq!(held.performance.bank().system().obstacle_residual_compliance().is_some(), residual, "{case}");
                let report = *held.performance.bank().system().nonlinear_preload_diagnostics()
                    .unwrap_or_else(|| panic!("{case}: no preload diagnostics"));
                assert!(report.force_residual <= report.force_tolerance,
                    "{case}: preload force residual {} against {}", report.force_residual, report.force_tolerance);
                let mut out = vec![0.0; held.spec.samples as usize];
                for chunk in out.chunks_mut(128) {
                    held.performance.block(chunk).unwrap_or_else(|error| panic!("{case}: the held render refused: {error:?}"));
                }
                let frame = held.performance.bank().system().last_frame().unwrap_or_else(|| panic!("{case}: no last frame"));
                assert!(frame.nonlinear.is_some(), "{case}: the held steps left the nonlinear owner");
                assert!(frame.energy_residual_j.abs() <= frame.energy_tolerance_j,
                    "{case}: energy residual {} against {}", frame.energy_residual_j, frame.energy_tolerance_j);
                if residual {
                    assert!(frame.obstacle_residual_energy_j > 0.0, "{case}: the held stop stores nothing in the omitted modes");
                }
            }
        }
    });
}

#[test]
fn what_the_residual_does_not_own_refuses_by_name() {
    scoped(|cx, gate| {
        // A stop on the bridge's own receiver keeps that receiver as its opposite participant;
        // the bridge correction is added beside it, not substituted for it.
        let receiver = flexible(&struck()).compile_with_profiles(cx, gate, None,
            Some(Profile::parse(STOP.replace("fixed", "bridge:0").as_bytes()).unwrap()),
            |name| Ok(pack(name, 7800.0, 200_000.0, None))).unwrap();
        let bank = receiver.performance.bank();
        let bridge = &bank.source().bank().strings()[0].connections[0].body;
        assert_eq!(bank.system().obstacle(0).unwrap().right.component, bank.cards().len()+bridge.component);
        assert!(bank.system().obstacle_residual_compliance().is_some());
    });
    for record in ["residual-flexibility", "residual-flexibility exact-static-green extra", "residual-flexibility modal"] {
        let text = pluck().replace("cavity neck", &format!("{record}\ncavity neck"));
        assert!(input::Spec::parse(text.as_bytes()).is_err(), "accepted {record:?}");
    }
}

fn bowed() -> String {
    include_str!("../pluck-and-bow.material-strings")
        .replace("48000 4096 0.01 128","48000 256 0.01 128")
        .replace("forces preload 1 1\nport 0 0.2 surface 0.002\nforce 0 0 0","forces retain 0 0")
        .replace("event\t1e-2\t5e-3","event\t1e-3\t5e-4")
        .replace("event\t2e-2\t1e-2","event\t2e-3\t1e-3")
        .replace("event\t4e-2","event\t4e-3")
}
fn bowed_residual(text:&str)->input::Spec {
    let text=text.replace("cavity neck","residual-flexibility exact-static-green\nbow-residual-flexibility exact-static-green\ncavity neck");
    input::Spec::parse(text.as_bytes()).unwrap()
}

#[test]
fn bow_residual_uses_material_geometry_and_reaches_the_original_pcm_writer() {
    scoped(|cx,gate| {
        let mut plain=compile(cx,gate,flexible(&bowed()),7800.0);
        let mut corrected=compile(cx,gate,bowed_residual(&bowed()),7800.0);
        let bank=corrected.performance.bank();
        let bridges=residual::residual_bridges(bank.source().bank().strings(),bank.cards()).unwrap();
        let (cross,g)=residual::stop_residual(&bridges,&[residual::ResidualStop {
            string:0,station_fraction:0.2,direction:[1.0,0.0],
        }]).unwrap();
        assert_eq!(bank.system().friction_residual_compliance(),Some(g.as_slice()));
        assert!(g[0]>0.0 && cross[0]!=0.0);
        // The port keeps the plain card's receiving surface (none at a fixed bow); the bridge's
        // body part is the separate correction +cross*B inside the port's complete map.
        assert_eq!(bank.system().port(0).unwrap().right,plain.performance.bank().system().port(0).unwrap().right);
        let bridge=&bank.source().bank().strings()[0].connections[0].body;
        let offset:usize=bank.system().components()[..bank.cards().len()+bridge.component].iter()
            .map(|c|c.modes().len()).sum();
        let column=bank.system().friction_port_column(0).unwrap();
        assert_eq!(column[offset..offset+bridge.shapes.len()].to_vec(),bridge.shapes.iter().map(|b|cross[0]*b).collect::<Vec<_>>());
        assert_eq!(bank.system().components()[0].modes(),plain.performance.bank().system().components()[0].modes());
        let a=pressure(&mut plain,37);
        let b=pressure(&mut corrected,1);
        assert!(b.iter().any(|p|p.abs()>1e-12));
        assert!(a.iter().zip(&b).any(|(a,b)|(a-b).abs()>1e-12));
        let mut encoded=compile(cx,gate,bowed_residual(&bowed()),7800.0);
        let (wav,_)=encoded.write(Cursor::new(Vec::<u8>::new()),gate).unwrap();
        assert_eq!(wav.into_inner(),fs_couple::pcm_wav::encode_pcm16_wav(&b,48000,0.01).unwrap().0);
        let mut fine=bowed_residual(&bowed());fine.strings[0].bending_modes=8;
        let fine=compile(cx,gate,fine,7800.0);
        assert!(fine.performance.bank().system().friction_residual_compliance().unwrap()[0]<g[0]);
    });
}

// MagentaBadger 12160: beside a FIXED stop the bow's omitted modes and the stop's are one joint
// compliance. Its stop/bow block enters the joint solve, and the ledger books the cross part of
// the joint form f^T G f / 2, f = [t; -R], at every sample of a stick-slip run with the stop
// engaging and releasing. The engine refuses any sample whose balance does not close.
#[test]
fn a_fixed_stop_and_the_bow_share_one_omitted_compliance_and_close_the_ledger() {
    scoped(|cx,gate| {
        let mut ready=stopped(cx,gate,bowed_residual(&bowed())).unwrap();
        let (g_ss,g_sb,g_bb)={
            let bank=ready.performance.bank();
            let bridges=residual::residual_bridges(bank.source().bank().strings(),bank.cards()).unwrap();
            let stop=residual::ResidualStop {string:0,station_fraction:0.43,direction:[1.0,0.0]};
            let bow=residual::ResidualStop {string:0,station_fraction:0.2,direction:[1.0,0.0]};
            let (_,joint)=residual::stop_residual(&bridges,&[stop,bow]).unwrap();
            let system=bank.system();
            // Every block of the one joint matrix, bit for bit: the stop's own block is still
            // the stop-only binding's, and the stop/bow block is bound beside it.
            assert_eq!(system.obstacle_residual_compliance(),Some(&joint[..1]));
            assert_eq!(system.obstacle_port_residual_compliance(),Some(&joint[1..2]));
            assert_eq!(system.friction_residual_compliance(),Some(&joint[3..4]));
            assert!(joint[1]!=0.0 && joint[1].to_bits()==joint[2].to_bits(),"stop/bow block {joint:?}");
            (joint[0],joint[1],joint[3])
        };
        let mut out=vec![0.0;ready.spec.samples as usize];
        let (mut sticking,mut sliding,mut engaged,mut released)=(0,0,0,0);
        let mut accepted:Option<(f64,f64)>=None;
        for sample in 0..out.len() {
            ready.performance.block(&mut out[sample..sample+1]).unwrap_or_else(|e| panic!("sample {sample}: {e:?}"));
            let system=ready.performance.bank().system();
            let frame=system.last_frame().unwrap();
            let stop=system.last_obstacle_frames().unwrap()[0];
            let (t,r)=(frame.ports[0].traction_n,stop.normal_force_n);
            assert!(frame.energy_residual_j.abs()<=frame.energy_tolerance_j,"sample {sample}: {frame:?}");
            let close=|got:f64,want:f64,scale:f64,what:&str| assert!((got-want).abs()<=1e-12*scale+1e-300,
                "sample {sample} {what}: {got:e} against {want:e} (t {t:e}, R {r:e})");
            let (own,cross)=(0.5*r*(g_ss*r),r*(g_sb*t));
            close(stop.residual_energy_j,own-cross,own.abs()+cross.abs(),"stop storage with the cross part");
            close(frame.friction_residual_energy_j,0.5*t*g_bb*t,0.5*t*t*g_bb,"bow storage");
            if let Some((t0,r0))=accepted {
                let (dt,dr)=(t-t0,r-r0);
                let (own,cross)=(0.5*dr*(g_ss*dr),dr*(g_sb*dt));
                close(stop.residual_loss_j,own-cross,own.abs()+cross.abs(),"stop numerical loss with the cross part");
                // Only the joint numerical loss is a nonnegative form.
                let joint=frame.friction_residual_loss_j+stop.residual_loss_j;
                assert!(joint>=-1e-12*(0.5*dt*dt*g_bb+own.abs()+cross.abs()),"sample {sample}: joint loss {joint:e}");
            }
            accepted=Some((t,r));
            match format!("{:?}",frame.ports[0].regime).as_str() { "Sticking"=>sticking+=1, "Sliding"=>sliding+=1, _=>{} }
            if r>0.0 && t!=0.0 { engaged+=1; } else if engaged>0 && r==0.0 { released+=1; }
        }
        eprintln!("G_ss {g_ss:e} G_sb {g_sb:e} G_bb {g_bb:e}: sticking {sticking}, sliding {sliding}, \
            engaged under traction {engaged}, released after {released}");
        assert!(sticking>0 && sliding>0 && engaged>0,"stick-slip with the stop engaged: {sticking} {sliding} {engaged}");
        assert!(out.iter().any(|p|p.abs()>1e-12));
    });
}

// Bead zrvrx (MagentaBadger 36437, 36451): the bows' omitted modes beside bound longitudinal waves,
// alone, beside a fixed stop held loaded from a static preload (the route of all--violin-B6, whose
// finger stop is pressed by a held load at its station), and on all-B6's quasistatic
// `bow_compliance` ribbon, whose three hair drives are series springs. The records go in the parser's order,
// the bow tail right after the support tail and before the strain law; the other order does not
// parse. Every nonlinear step books the bow tail's F G F / 2 and dF G dF / 2, and the stop's share
// of the joint form, in its own balance. On some steps the tail's loss, the change of its storage
// and the stop's cross part each exceed the step's energy tolerance: a balance missing the loss,
// the storage at both ends, or the cross part at the start refuses that step.
const LONGITUDINAL: &str = "longitudinal-waves 0.2 0.01 128 100000 1e-11 1e-10 2 1 2\n";
fn bowed_longitudinal(bow_tail: bool, strain: bool, held: bool) -> String {
    let records = format!("residual-flexibility exact-static-green\n{}{}cavity neck",
        if bow_tail { "bow-residual-flexibility exact-static-green\n" } else { "" }, if strain { LONGITUDINAL } else { "" });
    let text = bowed().replace("cavity neck", &records);
    // 0.5 N held at the stop's station from the static preload: no later force rows.
    if held { text.replace("forces retain 0 0", "forces preload 1 0\nport 0 0.43 centreline 0.5") } else { text }
}
fn parsed(text: &str) -> input::Spec { input::Spec::parse(text.as_bytes()).unwrap() }
// Each sample: the nonlinear balance closes and the bow storage is the bound tail's t G t / 2 (and
// with one bow port, the stop's R G_ss R / 2 - R G_sb t). Counts the frames whose tail loss, tail
// storage change and stop/bow cross part exceed the frame's energy tolerance, the sticking and
// sliding frames of port 0, and the frames where a series hair dashpot dissipates.
fn tail_ledger(ready: &mut Ready, g: &[f64], stop: Option<(f64, f64)>) -> ([usize; 6], Vec<f64>) {
    let (mut counts, mut stored) = ([0; 6], 0.0);
    let mut out = vec![0.0; ready.spec.samples as usize];
    for sample in 0..out.len() {
        ready.performance.block(&mut out[sample..sample+1]).unwrap_or_else(|e| panic!("sample {sample}: {e:?}"));
        let system = ready.performance.bank().system();
        let frame = system.last_frame().unwrap();
        assert!(frame.nonlinear.is_some(), "sample {sample}: the step left the nonlinear owner");
        assert!(frame.energy_residual_j.abs() <= frame.energy_tolerance_j, "sample {sample}: {frame:?}");
        let (t, tolerance): (Vec<f64>, f64) = (frame.ports.iter().map(|p| p.traction_n).collect(), frame.energy_tolerance_j);
        let terms = || (0..t.len()).flat_map(|i| (0..t.len()).map(move |j| (i, j))).map(|(i, j)| 0.5*t[i]*g[i*t.len()+j]*t[j]);
        let (bow, scale) = (terms().sum::<f64>(), terms().map(f64::abs).sum::<f64>());
        assert!((frame.friction_residual_energy_j-bow).abs() <= 1e-12*scale+1e-300,
            "sample {sample}: bow storage {:e} against {bow:e}", frame.friction_residual_energy_j);
        if frame.friction_residual_loss_j > tolerance { counts[0] += 1; }
        if (frame.friction_residual_energy_j-stored).abs() > tolerance { counts[1] += 1; }
        stored = frame.friction_residual_energy_j;
        if let Some((g_ss, g_sb)) = stop {
            let stop = system.last_obstacle_frames().unwrap()[0];
            let (r, got) = (stop.normal_force_n, stop.residual_energy_j);
            let (own, cross) = (0.5*r*(g_ss*r), r*(g_sb*t[0]));
            assert!((got-(own-cross)).abs() <= 1e-12*(own.abs()+cross.abs())+1e-300,
                "sample {sample}: stop storage {got:e} against {:e}", own-cross);
            if cross.abs() > tolerance && got != 0.0 && bow > 0.0 { counts[2] += 1; }
        }
        match format!("{:?}", frame.ports[0].regime).as_str() { "Sticking" => counts[3] += 1, "Sliding" => counts[4] += 1, _ => {} }
        if frame.ports.iter().any(|port| port.series_drive_dissipation_j > 0.0) { counts[5] += 1; }
    }
    assert!(out.iter().any(|p| p.abs() > 1e-12));
    (counts, out)
}

#[test]
fn the_bow_tail_beside_longitudinal_waves_closes_every_nonlinear_step_and_changes_the_sound() {
    scoped(|cx, gate| {
        let swapped = bowed().replace("cavity neck", &format!(
            "residual-flexibility exact-static-green\n{LONGITUDINAL}bow-residual-flexibility exact-static-green\ncavity neck"));
        assert!(input::Spec::parse(swapped.as_bytes()).is_err(), "the bow tail after the strain law parsed");
        let mut ready = compile(cx, gate, parsed(&bowed_longitudinal(true, true, false)), 7800.0);
        assert!(ready.performance.bank().source().bank().longitudinal_waves().is_some());
        let g = ready.performance.bank().system().friction_residual_compliance().unwrap().to_vec();
        let (g_bb, (counts, sound)) = (g[0], tail_ledger(&mut ready, &g, None));
        eprintln!("G_bb {g_bb:e}: loss past the tolerance {}, storage moved past it {}, sticking {}, sliding {}",
            counts[0], counts[1], counts[3], counts[4]);
        assert!(g_bb > 0.0 && counts[0] > 0 && counts[1] > 0 && counts[3] > 0 && counts[4] > 0, "{counts:?}");
        // Removing either record changes the sound.
        for (bow_tail, strain) in [(false, true), (true, false)] {
            let other = pressure(&mut compile(cx, gate, parsed(&bowed_longitudinal(bow_tail, strain, false)), 7800.0), 1);
            assert!(sound.iter().zip(&other).any(|(a, b)| (a-b).abs() > 1e-12), "bow tail {bow_tail}, strain {strain}");
        }
    });
}

#[test]
fn the_bow_tail_beside_a_held_stop_and_longitudinal_waves_closes_the_joint_form_every_step() {
    scoped(|cx, gate| {
        let mut ready = stopped(cx, gate, parsed(&bowed_longitudinal(true, true, true))).unwrap();
        let (g_ss, g_sb, g_bb) = {
            let system = ready.performance.bank().system();
            let report = *system.nonlinear_preload_diagnostics().unwrap();
            assert!(report.force_residual <= report.force_tolerance);
            (system.obstacle_residual_compliance().unwrap()[0], system.obstacle_port_residual_compliance().unwrap()[0],
                system.friction_residual_compliance().unwrap()[0])
        };
        let (counts, sound) = tail_ledger(&mut ready, &[g_bb], Some((g_ss, g_sb)));
        eprintln!("G_ss {g_ss:e} G_sb {g_sb:e} G_bb {g_bb:e}: loss past the tolerance {}, storage moved past it {}, \
            cross part past it {}, sticking {}, sliding {}", counts[0], counts[1], counts[2], counts[3], counts[4]);
        assert!(g_ss > 0.0 && g_sb != 0.0 && g_bb > 0.0 && counts[0] > 0 && counts[1] > 0 && counts[2] > 0, "{counts:?}");
        let plain = pressure(&mut stopped(cx, gate, parsed(&bowed_longitudinal(false, true, true))).unwrap(), 1);
        assert!(sound.iter().zip(&plain).any(|(a, b)| (a-b).abs() > 1e-12), "the bow tail leaves the stopped sound");
    });
}

#[test]
fn the_bow_tail_beside_a_compliant_ribbon_and_longitudinal_waves_closes_every_nonlinear_step() {
    scoped(|cx, gate| {
        // The string (16 modes: four per contact cell), ribbon and hair rows of compliant-ribbon.material-strings,
        // with fibers.rs's synthetic hair pack.
        let card = include_str!("../compliant-ribbon.material-strings");
        let row = |head: &str| card.lines().find(|line| line.starts_with(head)).unwrap();
        let make = |bow_tail| parsed(&bowed_longitudinal(bow_tail, true, false)
            .replace("string wire.fsmatpk 0.5 0.0003 20 4 3 0.8 3000 0.1", row("string "))
            .replace("forces retain 0 0", &format!("{}\n{}\nforces retain 0 0", row("bow_ribbon"), row("bow_compliance"))))
            .compile(cx, gate, |name| Ok(if name.starts_with("bow-hair-") { super::fibers::quasistatic_pack(7e9, Some(1e6)) }
                else { pack(name, 7800.0, 200_000.0, None) })).unwrap();
        let mut ready = make(true);
        assert!(ready.performance.bank().source().compliant_ribbon().is_some());
        let g = ready.performance.bank().system().friction_residual_compliance().unwrap().to_vec();
        let (counts, sound) = tail_ledger(&mut ready, &g, None);
        eprintln!("G {g:?}: loss past the tolerance {}, storage moved past it {}, hair dashpot dissipating {}",
            counts[0], counts[1], counts[5]);
        assert!(g.len() == 9 && counts[0] > 0 && counts[1] > 0 && counts[5] > 0, "{counts:?}");
        let plain = pressure(&mut make(false), 1);
        assert!(sound.iter().zip(&plain).any(|(a, b)| (a-b).abs() > 1e-12), "the bow tail leaves the hair's sound");
    });
}

#[test]
fn bow_residual_requires_the_matching_support_tail_and_refuses_moving_stops() {
    scoped(|cx,gate| {
        // A stop that translates along its normal is refused from the declarations, before any
        // network is prepared; the engine also refuses that binding after the joint residual.
        let moving=format!("{STOP}motion 0.001 0.1 1\nknot 12 0\n");
        let error=bowed_residual(&bowed()).compile_with_profiles(cx,gate,None,Some(Profile::parse(moving.as_bytes()).unwrap()),
            |name| Ok(pack(name,7800.0,200_000.0,None))).err().unwrap().to_string();
        assert!(error.contains("requires fixed stops"),"{error}");
    });
    for record in ["bow-residual-flexibility exact-static-green","residual-flexibility exact-static-green\nbow-residual-flexibility guessed",
        "residual-flexibility exact-static-green\nbow-residual-flexibility exact-static-green extra"] {
        assert!(input::Spec::parse(bowed().replace("cavity neck",&format!("{record}\ncavity neck")).as_bytes()).is_err());
    }
}
