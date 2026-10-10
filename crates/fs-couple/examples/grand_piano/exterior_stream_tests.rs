use super::*;
use super::super::{DiscreteStateSpace,encode_pcm16_wav_interleaved,render};
use super::super::super::{board,geometry,engine::radiation::{Model,Pole}};

const EMPTY:&str="sample,event,key,value\n";
const STRIKE:&str="sample,event,key,value\n0,note_on,69,2\n";
fn score(text:&str)->Performance {Performance::read(text,&[69],100_000).unwrap()}
fn piano(substeps:usize,loaded:bool)->Instrument {
    let course=geometry::demonstration_scale().unwrap()[48];
    let board=board::demonstration();
    let mut piano=Instrument::new(vec![course],&board,RATE,substeps,8,true).unwrap();
    if loaded {
        piano.configure_radiation(&Model {ports:board.len(),poles:vec![
            Pole {omega:1200.,zeta:0.2,coupling:vec![180.,-90.,50.,20.]},
        ]}).unwrap();
    }
    piano
}
fn baked(channels:usize)->Baked {
    // Analytic causal receiver fixtures exercise the actual filter and delay
    // owners. BEM/fit accuracy is separately admitted by Baked::from_samples.
    Baked {filters:(0..channels).map(|channel|(0..4).map(|mode| {
        let sign=if channel==1 && mode%2==0 {-1.}else{1.};
        DiscreteStateSpace {n:1,a:vec![0.6+0.05*mode as f64],b:vec![0.1],
            c:vec![sign*(mode+1) as f64*0.03],d:sign*0.001,
            e_leftover:0.,t_s:1./f64::from(RATE)}
    }).collect()).collect(),delays_s:(0..channels).map(|c|0.002+0.00101*c as f64).collect(),
        maximum_error:0.,worst_rms_error:0.}
}

#[test]
fn arbitrary_blocks_preserve_pressure_controls_and_passive_feedback_history() {
    let events="sample,event,key,value\n0,note_on,69,2\n255,sustain,0,1\n511,note_off,69,0\n768,sustain,0,0\n";
    // Both one-way mono and passive-feedback stereo use the admitted 4x
    // mechanical fixture and retain its shared causal decimator history.
    for (channels,loaded) in [(1,false),(2,true)] {
        let substeps=4;
        let baked=baked(channels);let mut pa=piano(substeps,loaded);let mut pb=piano(substeps,loaded);
        let mut a=ExteriorStream::new(&mut pa,score(events),&baked).unwrap();
        let mut b=ExteriorStream::new(&mut pb,score(events),&baked).unwrap();
        let mut whole=vec![0.;1600*channels];a.render_interleaved_block(&mut whole).unwrap();
        let mut split=vec![0.;whole.len()];let mut start=0;
        for size in [1,7,64,0,127,13,256].into_iter().cycle() {
            let end=(start+size).min(1600);
            b.render_interleaved_block(&mut split[start*channels..end*channels]).unwrap();
            start=end;if start==1600 {break;}
        }
        assert_eq!(whole,split);assert!(whole.iter().any(|x|x.abs()>1e-10));
        assert_eq!(a.sample_position(),1600);assert_eq!(b.sample_position(),1600);
        assert_eq!(a.instrument().bank.q,b.instrument().bank.q);
        assert_eq!(a.instrument().bank.v,b.instrument().bank.v);
        assert_eq!(a.instrument().accounting.input_work_j,b.instrument().accounting.input_work_j);
        assert_eq!(a.instrument().accounting.dissipated_j(),b.instrument().accounting.dissipated_j());
        assert_eq!(a.instrument().energy_j(),b.instrument().energy_j());
        assert_eq!(a.instrument().radiation_energy_j(),b.instrument().radiation_energy_j());
        if loaded {assert!(a.instrument().radiation_energy_j()>0.);}
    }
}

#[test]
fn stereo_observes_one_trajectory_and_bad_buffers_leave_the_clock_retriable() {
    let mono=baked(1);let stereo=baked(2);let mut pa=piano(4,false);let mut pb=piano(4,false);
    let mut a=ExteriorStream::new(&mut pa,score(STRIKE),&mono).unwrap();
    let mut b=ExteriorStream::new(&mut pb,score(STRIKE),&stereo).unwrap();
    let mut partial=[f64::NAN;3];let e=b.render_interleaved_block(&mut partial).unwrap_err();
    assert_eq!((e.completed_frames,e.sample),(0,0));assert_eq!(partial,[0.;3]);
    let mut downmix=[f64::NAN;8];assert!(b.render_block(&mut downmix).is_err());
    assert_eq!(downmix,[0.;8]);assert_eq!(b.instrument().accounting.input_work_j,0.);
    assert!(!b.failed);assert_eq!(b.sample_position(),0);
    let mut left=[0.;1200];let mut pair=[0.;2400];
    a.render_block(&mut left).unwrap();b.render_interleaved_block(&mut pair).unwrap();
    for (value,frame) in left.iter().zip(pair.chunks_exact(2)) {assert_eq!(*value,frame[0]);}
    assert!(pair.chunks_exact(2).any(|f|f[0]!=f[1]));
    assert_eq!(a.instrument().bank.q,b.instrument().bank.q);
    assert_eq!(a.instrument().bank.v,b.instrument().bank.v);
    assert_eq!(a.instrument().accounting.input_work_j,b.instrument().accounting.input_work_j);
}

#[test]
fn live_gestures_match_scheduled_events_and_prepared_buffers_survive_empty_blocks() {
    let baked=baked(2);let mut pa=piano(4,false);let mut pb=piano(4,false);
    let mut a=ExteriorStream::new(&mut pa,score(EMPTY),&baked).unwrap();
    let mut b=ExteriorStream::new(&mut pb,score("sample,event,key,value\n64,note_on,69,2\n"),&baked).unwrap();
    let buffers=(a.trace.as_ptr(),a.acceleration.as_ptr(),a.previous.as_ptr(),a.receivers.as_ptr());
    a.render_interleaved_block(&mut[]).unwrap();b.render_block(&mut[]).unwrap();
    assert_eq!(a.sample_position(),0);assert_eq!(b.sample_position(),0);
    let mut silence=[f64::NAN;128];
    a.render_interleaved_block(&mut silence).unwrap();assert_eq!(silence,[0.;128]);
    b.render_interleaved_block(&mut silence).unwrap();assert_eq!(silence,[0.;128]);
    assert_eq!(b.instrument().accounting.input_work_j,0.);
    a.instrument_mut().note_on(69,2.).unwrap();
    let mut live=[0.;2400];let mut scheduled=[0.;2400];
    a.render_interleaved_block(&mut live).unwrap();b.render_interleaved_block(&mut scheduled).unwrap();
    assert_eq!(live,scheduled);assert!(live.iter().any(|x|x.abs()>1e-10));
    assert_eq!(a.sample_position(),1264);assert_eq!(b.sample_position(),1264);
    assert_eq!(a.instrument().bank.q,b.instrument().bank.q);
    assert_eq!(buffers,(a.trace.as_ptr(),a.acceleration.as_ptr(),a.previous.as_ptr(),a.receivers.as_ptr()));
}

#[test]
fn execution_failure_reports_frame_prefix_silences_suffix_and_cannot_resume() {
    let baked=baked(2);let mut piano=piano(4,false);
    let events="sample,event,key,value\n0,note_on,69,2\n500,note_on,69,2\n";
    let mut stream=ExteriorStream::new(&mut piano,score(events),&baked).unwrap();
    stream.render_interleaved_block(&mut[0.;400]).unwrap();
    let mut block=[f64::NAN;1000];let e=stream.render_interleaved_block(&mut block).unwrap_err();
    assert_eq!((e.completed_frames,e.sample),(300,500));assert!(stream.failed);
    assert!(block[..600].iter().all(|v|v.is_finite()));
    assert!(block[..600].iter().any(|v|v.abs()>1e-10));assert_eq!(block[600..],[0.;400]);
    let q=stream.instrument().bank.q.clone();let v=stream.instrument().bank.v.clone();
    let e=stream.render_interleaved_block(&mut block).unwrap_err();
    assert_eq!((e.completed_frames,e.sample),(0,500));assert_eq!(block,[0.;1000]);
    assert_eq!(stream.instrument().bank.q,q);assert_eq!(stream.instrument().bank.v,v);
}

#[test]
fn late_receiver_failure_never_publishes_half_a_frame_or_restarts_advanced_history() {
    let baked=baked(2);let mut piano=piano(4,false);
    let mut stream=ExteriorStream::new(&mut piano,score(STRIKE),&baked).unwrap();
    stream.render_interleaved_block(&mut[0.;1024]).unwrap();
    // Deliberately break only the second prepared observer. Mechanics and the
    // first observer will advance before its complete-basis check refuses.
    stream.receivers[1].filters.pop();
    let previous=stream.receivers[0].filters[0].state().to_vec();
    let q=stream.instrument().bank.q.clone();let v=stream.instrument().bank.v.clone();
    let mut output=[f64::NAN;16];let e=stream.render_interleaved_block(&mut output).unwrap_err();
    assert_eq!((e.completed_frames,e.sample),(0,512));assert_eq!(output,[0.;16]);
    assert!(stream.instrument().bank.q!=q || stream.instrument().bank.v!=v);
    assert_ne!(stream.receivers[0].filters[0].state(),previous.as_slice());
    let q=stream.instrument().bank.q.clone();let v=stream.instrument().bank.v.clone();
    assert!(stream.render_interleaved_block(&mut output).is_err());
    assert_eq!(stream.instrument().bank.q,q);assert_eq!(stream.instrument().bank.v,v);
    assert_eq!(stream.sample_position(),512);
}

#[test]
fn offline_wav_uses_the_same_physical_stream_and_refuses_actual_pcm_clipping() {
    let baked=baked(2);let mut pa=piano(4,false);let mut pb=piano(4,false);let mut pc=piano(4,false);
    let mut pressure=vec![0.;4800];
    {
        let mut stream=ExteriorStream::new(&mut pa,score(STRIKE),&baked).unwrap();
        for block in pressure.chunks_mut(128) {stream.render_interleaved_block(block).unwrap();}
    }
    let peak=pressure.iter().fold(0.0_f64,|m,v|m.max(v.abs()));assert!(peak>1e-10);
    let full_scale=2.*peak+1.;
    let (reference,clips)=encode_pcm16_wav_interleaved(&pressure,RATE,2,full_scale).unwrap();
    assert_eq!(clips,0);
    let rendered=render(&mut pb,score(STRIKE),2400,&baked,full_scale).unwrap();
    assert_eq!(rendered.wav,reference);assert_eq!(rendered.peak_pa,peak);
    assert_eq!(pa.bank.q,pb.bank.q);assert_eq!(pa.bank.v,pb.bank.v);
    assert!(rendered.report.contains("0 clips"));
    let error=render(&mut pc,score(STRIKE),2400,&baked,peak/4.).err().expect("clipping must refuse WAV output");
    assert!(error.contains("clipped samples") && error.contains("no WAV published"));
    // Output conversion changes no mechanics and never hides clipping by gain.
    assert_eq!(pa.bank.q,pc.bank.q);assert_eq!(pa.bank.v,pc.bank.v);
}

#[test]
fn incomplete_receivers_and_clock_overflow_refuse_without_advancing_mechanics() {
    let mut invalid=baked(2);invalid.delays_s.pop();let mut piano=piano(4,false);
    assert!(ExteriorStream::new(&mut piano,score(STRIKE),&invalid).is_err());
    assert_eq!(piano.accounting.input_work_j,0.);
    let baked=baked(1);let mut stream=ExteriorStream::new(&mut piano,score(EMPTY),&baked).unwrap();
    stream.sample=u64::MAX;let mut output=[f64::NAN;1];
    let e=stream.render_block(&mut output).unwrap_err();
    assert_eq!((e.completed_frames,e.sample),(0,u64::MAX));assert_eq!(output,[0.]);
    assert!(stream.failed);assert_eq!(stream.instrument().accounting.input_work_j,0.);
}
