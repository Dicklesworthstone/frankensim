//! Prepared finite-body pressure blocks sharing one mechanical/control clock.
//! BEM solves, fitting and allocation precede construction; successful blocks
//! retain the complete modal/receiver/propagation history without allocating.
use super::{Baked,Decimator,Instrument,Performance,RATE,Receiver};
use super::super::audio::BlockError;

pub struct ExteriorStream<'piano,'baked> {
    piano:&'piano mut Instrument,
    score:Performance,
    receivers:Vec<Receiver<'baked>>,
    decimator:Decimator,
    trace:Vec<f64>,
    acceleration:Vec<f64>,
    previous:Vec<f64>,
    mechanical_rate:u32,
    sample:u64,
    failed:bool,
}
impl<'piano,'baked> ExteriorStream<'piano,'baked> {
    /// Start a prepared, fresh-at-rest instrument and its fitted receivers on
    /// the same 48 kHz clock. The baked model must describe this complete loaded
    /// board basis. Receiver filters borrow the immutable fitted coefficients.
    pub fn new(piano:&'piano mut Instrument,score:Performance,baked:&'baked Baked)
        ->Result<Self,String> {
        let modes=piano.bank.board_count;let trace_len=piano.board_trace_len();
        if piano.sample_rate()!=RATE || modes==0 || trace_len%modes!=0
            || !(1..=2).contains(&baked.filters.len())
            || baked.delays_s.len()!=baked.filters.len()
            || baked.filters.iter().any(|r|r.len()!=modes)
            || piano.bank.q.iter().chain(&piano.bank.v).any(|v|*v!=0.) {
            return Err("exterior stream requires a fresh resting piano and matching complete clocks/basis".into());
        }
        let substeps=trace_len/modes;
        if !(1..=16).contains(&substeps) || piano.bank.rate!=RATE*substeps as u32 {
            return Err("exterior stream requires matching mechanical and output clocks".into());
        }
        let receivers=baked.runtime()?;
        let decimator=Decimator::new(substeps,modes).map_err(|e|e.to_string())?;
        let mechanical_rate=piano.bank.rate;
        Ok(Self {piano,score,receivers,decimator,trace:vec![0.;trace_len],
            acceleration:vec![0.;trace_len],previous:vec![0.;modes],mechanical_rate,
            sample:0,failed:false})
    }
    pub fn channels(&self)->usize {self.receivers.len()}
    pub fn sample_position(&self)->u64 {self.sample}
    pub fn instrument(&self)->&Instrument {self.piano}
    /// Apply physical gestures between blocks. Replacing the prepared bank,
    /// its loaded basis or clocks requires constructing a new stream and bake.
    pub fn instrument_mut(&mut self)->&mut Instrument {self.piano}
    /// Causal observation latency, in addition to each receiver's flight delay.
    pub fn decimator_delay_frames(&self)->f64 {self.decimator.delay_output_frames()}

    fn next_frame(&mut self)->Result<[f64;2],String> {
        self.score.dispatch(self.sample,self.piano)?;
        self.piano.step_with_board_trace(&mut self.trace)
            .map_err(|e|format!("mechanics frame {}: {e}",self.sample))?;
        let modes=self.previous.len();
        for (frame,out) in self.trace.chunks_exact(modes)
            .zip(self.acceleration.chunks_exact_mut(modes)) {
            for i in 0..modes {
                out[i]=(frame[i]-self.previous[i])*f64::from(self.mechanical_rate);
                self.previous[i]=frame[i];
            }
        }
        let filtered=self.decimator.preview(&self.acceleration).map_err(|e|e.to_string())?;
        let mut next=[0.;2];
        for (value,receiver) in next.iter_mut().zip(&mut self.receivers) {
            *value=receiver.step(filtered)?;
        }
        self.decimator.commit();self.sample+=1;
        Ok(next)
    }

    /// Mono pressure in Pa. Empty blocks consume neither controls nor time.
    /// Stereo requires the interleaved API; implicit downmixing is refused.
    pub fn render_block(&mut self,output:&mut[f64])->Result<(),BlockError> {
        if output.is_empty() {return Ok(());}
        if self.channels()!=1 {
            output.fill(0.);
            return Err(BlockError {completed_frames:0,sample:self.sample,
                message:"stereo exterior stream needs render_interleaved_block; implicit downmix is refused".into()});
        }
        self.render_interleaved_block(output)
    }

    /// Complete frames in receiver order: L,R,L,R for stereo, mono otherwise.
    /// Arbitrary block lengths preserve all histories; success allocates
    /// nothing in this host. Output remains physical Pa without gain/clipping.
    ///
    /// A partial frame is refused before stepping and can be retried. Execution
    /// failures preserve the successful prefix, silence the failed frame and
    /// suffix, and latch the stream. Mechanics or an earlier receiver may have
    /// advanced before an observer fails; no cross-owner rollback is claimed.
    pub fn render_interleaved_block(&mut self,output:&mut[f64])->Result<(),BlockError> {
        if output.is_empty() {return Ok(());}
        let channels=self.channels();
        if output.len()%channels!=0 {
            output.fill(0.);
            return Err(BlockError {completed_frames:0,sample:self.sample,
                message:"exterior audio buffer ends inside a channel frame".into()});
        }
        let frames=output.len()/channels;
        if self.failed || self.sample.checked_add(frames as u64).is_none() {
            output.fill(0.);self.failed=true;
            return Err(BlockError {completed_frames:0,sample:self.sample,
                message:"exterior stream is faulted or its sample clock would overflow".into()});
        }
        if self.piano.sample_rate()!=RATE || self.piano.bank.rate!=self.mechanical_rate
            || self.piano.bank.board_count!=self.previous.len()
            || self.piano.board_trace_len()!=self.trace.len() {
            output.fill(0.);self.failed=true;
            return Err(BlockError {completed_frames:0,sample:self.sample,
                message:"exterior stream's prepared mechanical clocks or basis changed".into()});
        }
        for i in 0..frames {
            match self.next_frame() {
                Ok(values)=>output[i*channels..(i+1)*channels].copy_from_slice(&values[..channels]),
                Err(message)=>{
                    output[i*channels..].fill(0.);self.failed=true;
                    return Err(BlockError {completed_frames:i,sample:self.sample,message});
                }
            }
        }
        Ok(())
    }
}

#[cfg(test)]
#[path="exterior_stream_tests.rs"]
mod tests;
