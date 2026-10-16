//! Existing BEM transfers -> fs-vfit causal receivers -> piano substep motion -> Pa.
//! The offline fit is not a new acoustic solver. Fitting sees even frequencies;
//! odd frequencies are held out. No physics or validation bound changes on failure.
use super::{engine::Instrument,performance::Performance,exterior_geometry::{Samples,RATE}};
use fs_math::{c64::C64,det};
use fs_couple::pcm_wav::decimate::Decimator;
#[cfg(test)]
use fs_couple::pcm_wav::encode_pcm16_wav_interleaved;
use fs_vfit::{FitOptions,WeightPreset,vector_fit};
use fs_vfit::discretize::{bilinear_state_space,DiscreteStateSpace,DiscreteStateSpaceRuntime,DelayedFilter,DigitalFilter};

#[path="exterior_stream.rs"]
mod stream;
pub use stream::ExteriorStream;

const MAX_ERROR:f64=0.15;
const MAX_RMS:f64=0.05;
pub struct Baked {
    filters:Vec<Vec<DiscreteStateSpace>>,
    delays_s:Vec<f64>,
    pub maximum_error:f64,
    pub worst_rms_error:f64,
}
fn integer_flight_samples(delay_s:f64)->f64 {(delay_s*f64::from(RATE)).floor()}
fn fit(omega:&[f64],values:&[C64],delay:f64,order:usize)->Result<(DiscreteStateSpace,f64,f64),String> {
    let dt=1./f64::from(RATE);
    if values.len()!=omega.len() || omega.len()<17 || omega.len()>257 || omega.len()%2==0
        || omega.iter().enumerate().any(|(i,w)|!w.is_finite() || *w<=0. || *w*dt>=std::f64::consts::PI
            || (i>0 && *w<=omega[i-1]))
        || values.iter().any(|v|!v.re.is_finite() || !v.im.is_finite())
        || !delay.is_finite() || delay<2.*dt || delay>0.5
        || !(2..=32).contains(&order) || 2*order>(omega.len()+1)/2 {
        return Err("invalid complete exterior transfer, proper-fit or propagation budget".into());
    }
    // exp(-i omega t): peel only complete samples of guaranteed flight. The
    // residual fractional phase stays in the fitted/held-out transfer. Restoring
    // a fractional delay with DelayedFilter's linear interpolation would add an
    // unchecked low-pass (5.1 dB at 15 kHz for a half sample at 48 kHz).
    // Actual 1/r spreading and reactive near-field terms remain in the transfer.
    let integer_delay_s=integer_flight_samples(delay)*dt;
    let shifted:Vec<_>=omega.iter().zip(values).map(|(w,v)| {
        let phase = -w*integer_delay_s;*v*C64::new(det::cos(phase),det::sin(phase))
    }).collect();
    let scale=shifted.iter().map(|v|v.abs()).fold(0.0_f64,f64::max);
    if !scale.is_finite() {return Err("exterior transfer scale overflow".into());}
    if scale==0. {return Ok((DiscreteStateSpace {n:0,a:vec![],b:vec![],c:vec![],d:0.,e_leftover:0.,t_s:dt},0.,0.));}
    let warped:Vec<_>=omega.iter().step_by(2).map(|w|2./dt*det::tan(w*dt/2.)).collect();
    // fs-vfit uses exp(+i omega t), so conjugate AND Tustin-warp training data.
    let training:Vec<_>=shifted.iter().step_by(2).map(|v|C64::new(v.re/scale,-v.im/scale)).collect();
    let fitted=vector_fit(&warped,&training,&FitOptions {order,iterations:12,
        weights:WeightPreset::Uniform,fit_d:true,fit_e:false}).map_err(|e|e.to_string())?;
    if !fitted.model.is_stable() || fitted.model.e!=0. || !fitted.report.weighted_rms.is_finite() {
        return Err("exterior receiver fit is not finite, stable and proper".into());
    }
    let mut filter=bilinear_state_space(&fitted.model,dt,0.).map_err(|e|e.to_string())?;
    for c in &mut filter.c {*c*=scale;}filter.d*=scale;
    filter.try_runtime().map_err(|e|e.to_string())?;
    let mut maximum=0.0_f64;let mut square=0.;let mut held=0;
    for i in (1..omega.len()).step_by(2) {
        let error=(filter.eval(omega[i]).map_err(|e|e.to_string())?.conj()-shifted[i]).abs()/scale;
        if !error.is_finite() {return Err("nonfinite held-out exterior error".into());}
        maximum=maximum.max(error);square+=error*error;held+=1;
    }
    let rms=(square/held as f64).sqrt();
    if maximum>MAX_ERROR || rms>MAX_RMS {
        return Err(format!("exterior fit refuses held-out error max={maximum}, rms={rms}; limits={MAX_ERROR}/{MAX_RMS}; refine the supplied band/grid/order"));
    }
    Ok((filter,maximum,rms))
}
impl Baked {
    pub fn from_samples(samples:&Samples,order:usize)->Result<Self,String> {
        if !(1..=2).contains(&samples.values.len()) || samples.delays_s.len()!=samples.values.len()
            || samples.values[0].is_empty() || samples.values[0].len()>super::linear::MAX_BOARD_MODES
            || samples.values.iter().any(|r|r.len()!=samples.values[0].len()) {
            return Err("exterior transfer must retain every input at every receiver".into());
        }
        let mut filters=Vec::new();let mut maximum_error=0.0_f64;let mut worst_rms_error=0.0_f64;
        for (channel,rows) in samples.values.iter().enumerate() {
            let mut receiver=Vec::with_capacity(rows.len());
            for (input,row) in rows.iter().enumerate() {
                let (f,m,r)=fit(&samples.omega,row,samples.delays_s[channel],order)
                    .map_err(|e|format!("receiver {channel}, modal input {input}: {e}"))?;
                receiver.push(f);maximum_error=maximum_error.max(m);worst_rms_error=worst_rms_error.max(r);
            }
            filters.push(receiver);
        }
        Ok(Self {filters,delays_s:samples.delays_s.clone(),maximum_error,worst_rms_error})
    }
    fn runtime(&self)->Result<Vec<Receiver<'_>>,String> {
        self.filters.iter().zip(&self.delays_s).map(|(filters,&delay)| {
            let dt=1./f64::from(RATE);
            Ok(Receiver {filters:filters.iter().map(|f|f.try_runtime().map_err(|e|e.to_string())).collect::<Result<_,_>>()?,
                delay:DelayedFilter::new(integer_flight_samples(delay),DigitalFilter {sections:vec![],direct:1.,t_s:dt,prewarp:0.})
                    .map_err(|e|e.to_string())?})
        }).collect()
    }
}
struct Receiver<'a> {filters:Vec<DiscreteStateSpaceRuntime<'a>>,delay:DelayedFilter}
impl Receiver<'_> {
    fn step(&mut self,acceleration:&[f64])->Result<f64,String> {
        if acceleration.len()!=self.filters.len() || acceleration.iter().any(|a|!a.is_finite()) {
            return Err("finite exterior receiver needs every modal acceleration".into());
        }
        let mut pressure=0.;
        for (filter,&a) in self.filters.iter_mut().zip(acceleration) {
            pressure+=filter.step(a).map_err(|e|e.to_string())?;
        }
        if !pressure.is_finite() {return Err("exterior pressure sum overflow".into());}
        self.delay.push(pressure).map_err(|e|e.to_string())
    }
}
pub struct Rendered {pub wav:Vec<u8>,pub report:String,pub peak_pa:f64}

/// Offline, fresh-at-rest instrument, one event/decimation clock for all
/// receivers. Failures publish no WAV; accepted mechanics may have advanced,
/// so discard this candidate rather than claiming cross-owner rollback.
pub fn render(piano:&mut Instrument,score:Performance,frames:usize,baked:&Baked,full_scale_pa:f64,
    format:super::wav_output::Format)
    ->Result<Rendered,String> {
    if !(2400..=2_880_000).contains(&frames) || !full_scale_pa.is_finite() || full_scale_pa<=0. {
        return Err("exterior render requires a fresh resting piano and matching complete clocks/basis".into());
    }
    let (pressure,channels,decimator_delay)={
        let mut stream=ExteriorStream::new(piano,score,baked)?;
        let channels=stream.channels();let mut pressure=vec![0.;frames*channels];
        stream.render_interleaved_block(&mut pressure).map_err(|e|e.to_string())?;
        (pressure,channels,stream.decimator_delay_frames())
    };
    let wav=format.encode(&pressure,RATE,channels as u16,full_scale_pa)?;
    let peak_pa=pressure.iter().fold(0.0_f64,|m,v|m.max(v.abs()));
    let coupling=if piano.has_radiation() {
        "Passive acoustic feedback, second-order substep splitting"
    } else {"One-way acoustics, no radiation backreaction"};
    let encoding=format.label();
    let mut report=format!("{frames} frames, {channels} receivers, {RATE} Hz, peak {peak_pa:e} Pa; {encoding}, 0 clips at {full_scale_pa} Pa full scale, no normalization.\nHeld-out transfer max={}, worst modal RMS={}; flight lower bounds {:?} s, decimator delay {} frames.\nInput {} J; combined stored {} J; total loss {} J; combined closure {} J. Acoustic storage {} J and acoustic loss {} J (exterior) are included once in those totals. Cavity storage {} J and cavity momentum loss {} J are also included once. {coupling}. No flexible lid/cabinet, room or accuracy outside the sampled band.",
        baked.maximum_error,baked.worst_rms_error,baked.delays_s,decimator_delay,
        piano.accounting.input_work_j,piano.energy_j(),piano.accounting.dissipated_j(),
        piano.accounting.input_work_j-piano.energy_j()-piano.accounting.dissipated_j(),
        piano.radiation_energy_j(),piano.accounting.radiation_loss_j,
        piano.cavity_energy_j(),piano.accounting.cavity_loss_j);
    if format==super::wav_output::Format::Float32 && peak_pa>full_scale_pa {
        report.push_str("\nFloat WAV preserves pressure above nominal digital full scale; attenuate before fixed-point playback.");
    }
    Ok(Rendered {wav,report,peak_pa})
}

#[cfg(test)]
mod tests {
    use super::*;
    fn omega()->Vec<f64> {(0..41).map(|i|std::f64::consts::TAU*(40.+i as f64*15.)).collect()}
    #[test]
    fn phase_convention_and_held_out_gate_are_not_training_residuals() {
        let w=omega();let delay=144.5/f64::from(RATE);
        // Independent two-real-pole passive lowpass with known propagation.
        let values:Vec<_>=w.iter().map(|w| {
            let h=C64::new(300.,0.)/C64::new(300.,-*w)+C64::new(700.,0.)/C64::new(700.,-*w);
            let phase=w*delay;h*C64::new(phase.cos(),phase.sin())
        }).collect();
        let (filter,m,r)=fit(&w,&values,delay,4).unwrap();
        assert!(m<0.01 && r<0.005);
        for i in (1..w.len()).step_by(2) {
            let phase=w[i]*integer_flight_samples(delay)/f64::from(RATE);
            let actual=filter.eval(w[i]).unwrap().conj()*C64::new(phase.cos(),phase.sin());
            assert!((actual-values[i]).abs()<0.01);
        }
        let mut corrupt=values.clone();
        for v in corrupt.iter_mut().skip(1).step_by(2) {*v = -*v;}
        assert!(fit(&w,&corrupt,delay,4).is_err());
        assert!(fit(&w,&values,0.,4).is_err());
    }
    #[test]
    fn fractional_flight_preserves_high_partials_in_the_actual_receiver_runtime() {
        let rate=f64::from(RATE);
        let w:Vec<_>=(0..129).map(|i|std::f64::consts::TAU*(100.+i as f64*160.)).collect();
        // Analytic one-way pressure transfer with constant SI gain and known
        // travel. No fit coefficients or discrete-delay implementation form
        // the reference. Both fractional phase and level must reach playback.
        for fraction in [0.,0.25,0.5,0.875] {
            let delay=(144.+fraction)/rate;
            let values:Vec<_>=w.iter().map(|omega| {
                let phase=omega*delay;C64::new(0.8*phase.cos(),0.8*phase.sin())
            }).collect();
            let baked=Baked::from_samples(&Samples {omega:w.clone(),values:vec![vec![values]],
                delays_s:vec![delay],minimum_ppw:8.,maximum_condition_lower_bound:1.},16)
                .unwrap_or_else(|error|panic!("fraction={fraction}: {error}"));
            assert!(baked.maximum_error<0.05 && baked.worst_rms_error<0.025);
            for frequency in [4800.,15_000.,20_000.] {
                let mut receivers=baked.runtime().unwrap();
                let omega=std::f64::consts::TAU*frequency;
                let (mut error,mut signal)=(0.,0.);
                for n in 0..4096 {
                    let time=n as f64/rate;
                    let pressure=receivers[0].step(&[(omega*time).sin()]).unwrap();
                    if n>1024 {
                        let reference=0.8*(omega*(time-delay)).sin();
                        error+=(pressure-reference).powi(2);signal+=reference.powi(2);
                    }
                }
                let relative_rms=(error/signal).sqrt();
                assert!(relative_rms<0.05,"fraction={fraction}, frequency={frequency}, error={relative_rms}");
            }
        }
    }
    #[test]
    fn zero_transfer_remains_silent_and_nonfinite_samples_never_create_a_filter() {
        let w=omega();let zero=vec![C64::new(0.,0.);w.len()];
        let (filter,m,r)=fit(&w,&zero,0.003,4).unwrap();
        assert_eq!((m,r),(0.,0.));let mut runtime=filter.try_runtime().unwrap();
        for n in 0..100 {assert_eq!(runtime.step(n as f64).unwrap(),0.);}
        let mut bad=zero;bad[3]=C64::new(f64::NAN,0.);assert!(fit(&w,&bad,0.003,4).is_err());
    }
}
