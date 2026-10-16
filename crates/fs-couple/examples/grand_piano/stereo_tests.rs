//! One actual hammer/string/board engine, two independent spatial observers.
use super::*;
const POSITIONS:[[f64;3];2]=[[-0.4,0.1,0.8],[0.65,0.15,1.2]];
const SCORE:&str="sample,event,key,value\n0,note_on,69,2\n255,sustain,0,1\n511,note_off,69,0\n768,sustain,0,0\n";
fn parts()->(Instrument,Vec<SurfaceSample>) {
    let course=super::super::geometry::demonstration_scale().unwrap()[48];
    let board=super::super::board::demonstration();
    // Explicit kinematic surface fixture; not a measured/eigensolved panel.
    let surface=[-0.25,0.25].into_iter().enumerate().map(|(side,x)|SurfaceSample {
        position_m:[x,0.0,0.0],area_m2:0.1,
        mode_shape:board.iter().enumerate().map(|(i,b)|b.volume/0.2
            *if (i+side)%2==0 {1.0}else{-0.5}).collect(),
    }).collect();
    (Instrument::new(vec![course],&board,48_000,4,8,true).unwrap(),surface)
}
fn score(text:&str)->Performance {Performance::read(text,&[69],100_000).unwrap()}
fn stereo(text:&str,positions:[[f64;3];2])->AudioStream {
    let (instrument,surface)=parts();
    AudioStream::new_stereo(instrument,score(text),&surface,positions,Medium::air()).unwrap()
}
fn mono(position:[f64;3])->AudioStream {
    let (instrument,surface)=parts();
    AudioStream::new(instrument,score(SCORE),Some(&surface),position,Medium::air(),1000.0).unwrap()
}

#[test]
fn stereo_equals_two_independent_mics_but_advances_mechanics_only_once() {
    let mut actual=stereo(SCORE,POSITIONS);let mut left=mono(POSITIONS[0]);let mut right=mono(POSITIONS[1]);
    let frames=1400;let mut out=vec![0.0;2*frames];let mut a=vec![0.0;frames];let mut b=a.clone();
    actual.render_interleaved_block(&mut out).unwrap();
    left.render_block(&mut a).unwrap();right.render_block(&mut b).unwrap();
    for (i,frame) in out.chunks_exact(2).enumerate() {assert_eq!(frame,[a[i],b[i]]);}
    assert!(a.iter().zip(&b).any(|(l,r)|(l-r).abs()>1e-10),"not dual mono or static panning");
    assert!(a.iter().any(|p|p.abs()>1e-10) && b.iter().any(|p|p.abs()>1e-10));
    assert_eq!(actual.sample_position(),frames as u64);
    assert_eq!(actual.instrument.bank.q,left.instrument.bank.q);
    assert_eq!(actual.instrument.bank.v,right.instrument.bank.v);
    assert_eq!(actual.instrument.accounting.input_work_j,left.instrument.accounting.input_work_j);
    assert_eq!(actual.instrument.energy_j(),left.instrument.energy_j());
    let (wav,_)=fs_couple::pcm_wav::encode_pcm16_wav_interleaved(&out,48_000,2,2.0).unwrap();
    assert_eq!(wav.len(),44+4*frames);assert_eq!(&wav[22..24],&2_u16.to_le_bytes());
}

#[test]
fn stereo_repartitioning_and_coincident_receivers_keep_frame_and_filter_history() {
    let positions=[POSITIONS[0];2];let mut whole=stereo(SCORE,positions);let mut split=stereo(SCORE,positions);
    let mut a=vec![0.0;2400];let mut b=a.clone();let trace=split.board_trace.as_ptr();
    whole.render_interleaved_block(&mut a).unwrap();
    let mut start=0;
    for frames in [1,7,0,64,127,13,256].into_iter().cycle() {
        let end=(start+frames*2).min(b.len());
        split.render_interleaved_block(&mut b[start..end]).unwrap();start=end;
        if start==b.len(){break;}
    }
    assert_eq!(a,b);assert!(a.chunks_exact(2).all(|frame|frame[0]==frame[1]));
    assert_eq!(split.sample_position(),1200);assert_eq!(trace,split.board_trace.as_ptr());
    assert_eq!(whole.instrument.bank.q,split.instrument.bank.q);
}

#[test]
fn incomplete_stereo_buffers_do_not_consume_events_and_failures_silence_whole_frames() {
    let mut stream=stereo("sample,event,key,value\n0,note_on,69,2\n5,note_on,69,2\n",POSITIONS);
    let mut odd=[f64::NAN;3];assert!(stream.render_interleaved_block(&mut odd).is_err());
    assert_eq!(odd,[0.0;3]);assert_eq!(stream.sample_position(),0);
    assert!(stream.render_block(&mut [0.0;2]).is_err()); // No implicit downmix.
    assert_eq!(stream.instrument.accounting.input_work_j,0.0);
    let mut block=[f64::NAN;32];let error=stream.render_interleaved_block(&mut block).unwrap_err();
    assert_eq!((error.completed_frames,error.sample),(5,5));
    assert!(block.iter().all(|p|p.is_finite()));assert_eq!(&block[10..],&[0.0;22]);
    let q=stream.instrument.bank.q.clone();let v=stream.instrument.bank.v.clone();
    assert!(stream.render_interleaved_block(&mut block).is_err());
    assert_eq!(stream.sample_position(),5);assert_eq!(stream.instrument.bank.q,q);assert_eq!(stream.instrument.bank.v,v);
    assert_eq!(block,[0.0;32]);
}

#[test]
fn invalid_second_receiver_is_rejected_during_cold_preparation() {
    for position in [[0.0,0.0,0.0],[0.0,0.0,f64::NAN],[1e100,0.0,1.0]] {
        let (instrument,surface)=parts();
        assert!(AudioStream::new_stereo(instrument,score(SCORE),&surface,
            [POSITIONS[0],position],Medium::air()).is_err());
    }
}
