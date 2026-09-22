//! The actual main CLI preparation seams, not a separate piano implementation.
use super::*;
use std::fmt::Write;

fn options(args:&[&str])->Result<Options,String> {
    Options::parse(&args.iter().map(|s|s.to_string()).collect::<Vec<_>>())
}
fn shell()->String {
    let mut text=String::from("frankensim-board-geometry-si-v1\nsource,estimated,main CLI crown regression\nsupport,clamped\npretension,0\ndamping,0.01\n");
    for (i,p) in [[0.,0.],[1.,0.],[1.,1.],[0.,1.],[0.5,0.5]].iter().enumerate() {
        writeln!(text,"node,{i},{},{}",p[0],p[1]).unwrap();
    }
    for (i,t) in [[0,1,4],[1,2,4],[2,3,4],[3,0,4]].iter().enumerate() {
        writeln!(text,"triangle,{i},{},{},{},0.008,450,1e10,8e8,0.3,6e8,0.27",t[0],t[1],t[2]).unwrap();
    }
    text.push_str("fixed,0\nfixed,1\nfixed,2\nfixed,3\nbridge,69,0,0,0,1\n");
    crowned_board::elevate(&text,&[0.,0.,0.,0.,0.015],"authored test crown").unwrap()
}

#[test]
fn supplied_geometry_overrides_only_the_preset_board_not_the_performance_or_materials() {
    let o=options(&["--preset","steinway-d","--board-geometry","measured.fss",
        "--render","piano.wav","--midi","performance.mid","--midi-channel","2",
        "--midi-half-pedal","--hammers","felt.fsh","--dampers","pads.fspd",
        "--concert-pitch","442","--board-band-hz","800","--modes","128",
        "--microphone","0.3,0.7,1.2"]).unwrap();
    assert!(!o.uses_preset_board());
    assert_eq!(o.preset.as_deref(),Some("steinway-d"));
    assert_eq!(o.board_geometry.as_deref(),Some("measured.fss"));
    assert_eq!(o.hammers.as_deref(),Some("felt.fsh"));
    assert_eq!(o.dampers.as_deref(),Some("pads.fspd"));
    assert_eq!(o.midi_mapping.channel,1);assert!(o.midi_mapping.continuous_sustain);
    assert_eq!(o.microphone,Some([0.3,0.7,1.2]));
    assert_eq!(o.board_band_hz,800.);assert_eq!(o.modes,128);
    let scale=selected_scale(None,&o).unwrap();
    assert_eq!(scale.len(),88);
    assert!((scale[48].partial_hz(1,scale[48].tension_n)-442.).abs()<1e-8);
    assert!(options(&["--preset","steinway-d"]).unwrap().uses_preset_board());
}

#[test]
fn board_override_never_allows_conflicting_export_or_silent_preset_fallback() {
    for extra in [vec!["--dump-geometry","generated.fsb"],vec!["--dump-obj","generated.obj"],
        vec!["--mesh-divisions","8"],vec!["--board","modal.csv"],
        vec!["--render","supplied.fss"]] {
        let mut args=vec!["--preset","steinway-d","--board-geometry","supplied.fss"];
        args.extend(extra);assert!(options(&args).is_err(),"accepted {args:?}");
    }
    assert!(prepare_geometric_board("not a supplied board",&[69],400.).is_err());
    assert!(prepare_geometric_board(&shell(),&[60,69],400.).is_err());
    assert!(prepare_geometric_board(&shell().replace("pretension,0","pretension,10"),&[69],400.).is_err());
}

#[test]
fn crowned_geometry_custom_felt_and_spatial_release_share_the_existing_main_audio_path() {
    let o=options(&["--preset","steinway-d","--board-geometry","supplied.fss",
        "--render","piano.wav","--hammers","felt.fsh","--dampers","estimated",
        "--concert-pitch","442","--modes","12","--microphone","0.3,0.7,1.2"]).unwrap();
    let course=selected_scale(None,&o).unwrap()[48];
    let board=prepare_geometric_board(&shell(),&[69],o.board_band_hz).unwrap();
    let card="frankensim-hammer-materials-v1\nfelt,69,400000,0.2,2.5,3.2,0.25,0.8,2500000\n";
    let run=|split:bool| {
        let mut piano=prepare_instrument_with_hammers(vec![course],&board.modes,&o,Some(card)).unwrap();
        let damper=linear::dampers::Specification::estimated(&[course]).unwrap();
        piano.configure_dampers(&damper).unwrap();
        let score=performance::Performance::read("sample,event,key,value\n0,note_on,69,2\n600,note_off,69,0\n",&[69],1800).unwrap();
        let mut stream=audio::AudioStream::new(piano,score,Some(&board.surface),o.microphone.unwrap(),
            fs_bem::helmholtz::Medium::air(),o.observer_gain).unwrap();
        let mut output=vec![0.;1800];
        if split {for block in output.chunks_mut(127) {stream.render_block(block).unwrap();}}
        else {stream.render_block(&mut output).unwrap();}
        let piano=stream.instrument();
        assert!(piano.damper_resolution().is_some());
        assert!(piano.accounting.felt_loss_j>0.);
        assert_eq!(piano.accounting.felt_relaxation_loss_j,0.);
        assert!(piano.accounting.shank_loss_j>0.);
        assert!(piano.accounting.damper_loss_j>0.);
        assert!((piano.accounting.input_work_j-piano.energy_j()-piano.accounting.dissipated_j()).abs()<1e-7);
        assert!(output.iter().all(|p|p.is_finite()));
        assert!(output.iter().any(|p|p.abs()>1e-10));
        output
    };
    assert_eq!(run(false),run(true));
}
