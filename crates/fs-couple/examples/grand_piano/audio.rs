//! Prepared block-output composition. Geometry, eigensolves, file parsing and
//! allocation happen before construction, never in the successful block loop.
//! The same sample path serves offline WAV rendering and a host audio callback.
//! No audio-device backend or measured real-time deadline is implied here.
use super::{board_geometry::SurfaceSample, engine::Instrument, microphone::Microphone,
    performance::Performance};
use fs_bem::helmholtz::Medium;

#[derive(Debug)]
pub struct BlockError {
    /// Successfully emitted frames in THIS block; the remainder is silenced.
    pub completed_frames: usize,
    /// Absolute audio sample at which the failure occurred.
    pub sample: u64,
    pub message: String,
}
impl std::fmt::Display for BlockError {
    fn fmt(&self, f:&mut std::fmt::Formatter<'_>)->std::fmt::Result {
        write!(f,"audio sample {} after {} frames in block: {}",self.sample,self.completed_frames,self.message)
    }
}
impl std::error::Error for BlockError {}

pub struct AudioStream {
    instrument: Instrument, microphone: Option<Microphone>, right_microphone: Option<Microphone>, score: Performance,
    board_trace: Vec<f64>, diagnostic_gain: f64, sample: u64, failed: bool,
}
impl AudioStream {
    /// A surface selects physical multirate pressure; None selects the existing
    /// explicitly diagnostic volume observer. Gain affects ONLY the latter.
    pub fn new(instrument:Instrument, score:Performance, surface:Option<&[SurfaceSample]>,
        position_m:[f64;3], medium:Medium, diagnostic_gain:f64)->Result<Self,String> {
        if !diagnostic_gain.is_finite() || diagnostic_gain<=0.0 {
            return Err("diagnostic observer gain must be finite and positive".into());
        }
        let microphone=surface.map(|surface|Microphone::new_multirate(surface,&instrument.bank,
            instrument.sample_rate(),position_m,medium)).transpose()?;
        let board_trace=vec![0.0;instrument.board_trace_len()];
        Ok(Self{instrument,microphone,right_microphone:None,score,board_trace,diagnostic_gain,sample:0,failed:false})
    }
    /// Two physical Rayleigh receivers of ONE instrument and ONE control clock.
    /// The complete substep trace is computed once, then observed independently
    /// at each position. Construction is cold; no synthetic stereo gain is used.
    pub fn new_stereo(instrument:Instrument, score:Performance, surface:&[SurfaceSample],
        positions_m:[[f64;3];2], medium:Medium)->Result<Self,String> {
        let mut stream=Self::new(instrument,score,Some(surface),positions_m[0],medium,1.0)?;
        stream.right_microphone=Some(Microphone::new_multirate(surface,&stream.instrument.bank,
            stream.instrument.sample_rate(),positions_m[1],medium)?);
        Ok(stream)
    }
    pub fn channels(&self)->usize {if self.right_microphone.is_some() {2}else{1}}
    pub fn stereo_microphones(&self)->Option<[&Microphone;2]> {
        Some([self.microphone.as_ref()?,self.right_microphone.as_ref()?])
    }
    pub fn instrument(&self)->&Instrument {&self.instrument}
    /// Host gesture updates may be applied between blocks. Scheduled gestures
    /// remain sample-accurate inside blocks and use the same physical controls.
    pub fn instrument_mut(&mut self)->&mut Instrument {&mut self.instrument}
    pub fn microphone(&self)->Option<&Microphone> {self.microphone.as_ref()}
    pub fn sample_position(&self)->u64 {self.sample}

    fn next_frame(&mut self)->Result<[f64;2],String> {
        self.score.dispatch(self.sample,&mut self.instrument)?;
        let volume=if self.microphone.is_some() {
            self.instrument.step_with_board_trace(&mut self.board_trace)
        } else {self.instrument.step()}.map_err(|e|e.to_string())?;
        let pressure=match &mut self.microphone {
            Some(mic)=>mic.step_trace(&self.board_trace)?,
            None=>self.diagnostic_gain*volume,
        };
        let right=match &mut self.right_microphone {
            Some(mic)=>mic.step_trace(&self.board_trace)?,
            None=>0.0,
        };
        if !pressure.is_finite() || !right.is_finite() {return Err("audio observer overflow".into());}
        self.sample+=1;
        Ok([pressure,right])
    }

    /// Arbitrary block sizes, including zero. Controls and histories never reset
    /// at a block boundary. No allocation, locks, files or logging on success.
    /// On failure, the successful prefix remains valid, the rest becomes zero,
    /// and the stream latches failed. Reconstruct to recover: an observer failure
    /// can follow an accepted mechanics step, so resuming would misalign clocks.
    pub fn render_block(&mut self, output:&mut[f64])->Result<(),BlockError> {
        if output.is_empty() {return Ok(());}
        if self.channels()!=1 {
            output.fill(0.0);
            return Err(BlockError{completed_frames:0,sample:self.sample,
                message:"stereo stream needs render_interleaved_block; implicit downmix is refused".into()});
        }
        self.render_interleaved_block(output)
    }

    /// Frame-interleaved output: L,R,L,R for stereo, unchanged samples for mono.
    /// A partial frame refuses BEFORE events/mechanics, without faulting the
    /// stream. An execution failure publishes no part of the failed frame and
    /// latches the stream, including when the second receiver failed after the
    /// first advanced. completed_frames and sample are always FRAME counts.
    /// Successful stepping allocates nothing in this host.
    pub fn render_interleaved_block(&mut self, output:&mut[f64])->Result<(),BlockError> {
        if output.is_empty() {return Ok(());}
        let channels=self.channels();
        if output.len()%channels!=0 {
            output.fill(0.0);
            return Err(BlockError{completed_frames:0,sample:self.sample,
                message:"audio buffer ends inside a channel frame".into()});
        }
        let frames=output.len()/channels;
        if self.failed || self.sample.checked_add(frames as u64).is_none() {
            output.fill(0.0);self.failed=true;
            return Err(BlockError{completed_frames:0,sample:self.sample,
                message:"stream is faulted or its sample clock would overflow".into()});
        }
        for i in 0..frames {
            match self.next_frame() {
                Ok(values)=>output[i*channels..(i+1)*channels].copy_from_slice(&values[..channels]),
                Err(message)=>{
                    output[i*channels..].fill(0.0);self.failed=true;
                    return Err(BlockError{completed_frames:i,sample:self.sample,message});
                }
            }
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    fn stream(score:&str, physical:bool)->AudioStream {
        let course=super::super::geometry::demonstration_scale().unwrap()[48];
        let board=super::super::board::demonstration();
        let surface=[SurfaceSample{position_m:[0.0;3],area_m2:0.2,
            mode_shape:board.iter().map(|b|b.volume/0.2).collect()}];
        let piano=Instrument::new(vec![course],&board,48_000,4,8,true).unwrap();
        let performance=Performance::read(score,&[69],100_000).unwrap();
        AudioStream::new(piano,performance,physical.then_some(&surface[..]),[0.,0.,1.],Medium::air(),1000.).unwrap()
    }
    #[test]
    fn arbitrary_blocks_reproduce_the_same_physical_pressure_and_control_timing() {
        let score="sample,event,key,value\n0,note_on,69,2\n255,sustain,0,1\n511,note_off,69,0\n768,sustain,0,0\n";
        for physical in [false,true] {
            let mut a=stream(score,physical);let mut b=stream(score,physical);
            let mut whole=vec![0.0;1600];a.render_block(&mut whole).unwrap();
            let mut split=vec![0.0;1600];let mut start=0;
            for size in [1,7,64,0,127,13,256].into_iter().cycle() {
                let end=(start+size).min(split.len());b.render_block(&mut split[start..end]).unwrap();
                start=end;if start==split.len(){break;}
            }
            assert_eq!(whole,split);assert!(whole.iter().any(|x|x.abs()>1e-10));
            assert_eq!(a.sample_position(),1600);assert_eq!(a.instrument.bank.q,b.instrument.bank.q);
            assert_eq!(a.instrument.bank.v,b.instrument.bank.v);
        }
    }
    #[test]
    fn failure_keeps_the_valid_prefix_and_silences_the_rest_without_resuming_time() {
        let mut a=stream("sample,event,key,value\n0,note_on,69,2\n5,note_on,69,2\n",true);
        let mut block=[f64::NAN;32];let error=a.render_block(&mut block).unwrap_err();
        assert_eq!(error.completed_frames,5);assert_eq!(error.sample,5);
        assert!(block.iter().all(|v|v.is_finite()));assert!(block[5..].iter().all(|v|*v==0.0));
        let q=a.instrument.bank.q.clone();let v=a.instrument.bank.v.clone();
        assert!(a.render_block(&mut block).is_err());assert_eq!(a.sample_position(),5);
        assert_eq!(a.instrument.bank.q,q);assert_eq!(a.instrument.bank.v,v);
        assert!(block.iter().all(|v|*v==0.0));
    }
    #[test]
    fn empty_callbacks_consume_no_events_and_preserve_prepared_buffers() {
        let mut a=stream("sample,event,key,value\n64,note_on,69,2\n",true);
        let trace=a.board_trace.as_ptr();a.render_block(&mut[]).unwrap();assert_eq!(a.sample_position(),0);
        a.render_block(&mut[0.0;64]).unwrap();assert_eq!(a.instrument.accounting.input_work_j,0.0);
        a.render_block(&mut[0.0;1]).unwrap();assert!(a.instrument.accounting.input_work_j>0.0);
        assert_eq!(trace,a.board_trace.as_ptr());
    }
}

#[cfg(test)]
#[path = "stereo_tests.rs"]
mod stereo_tests;
