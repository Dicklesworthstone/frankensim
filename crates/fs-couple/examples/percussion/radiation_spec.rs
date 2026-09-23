//! Explicit offline acoustic bandwidth and cold-work controls. These never
//! retune modes, alter a contact, or change a mechanical or PCM sample clock.
use super::*;
use std::io::Read;

pub const HEADER: &str = "frankensim-radiation-preparation-v1";
const MAX_BYTES: usize = 4096;
const MAX_REFINED_PANELS: usize = 4096;
const MAX_DENSE_WORK: u64 = 100_000_000_000_000;

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Spec {
    pub band_hz: [f64;2],
    pub training_intervals: usize,
    pub max_order: usize,
    pub subdivisions: u32,
    pub max_panels: usize,
    pub max_dense_work: u64,
}
impl Default for Spec {
    fn default() -> Self {
        Self { band_hz: [40.0,1640.0], training_intervals: 20, max_order: 8,
            subdivisions: 0, max_panels: MAX_PANELS, max_dense_work: 1_000_000_000_000 }
    }
}
impl Spec {
    pub fn validate(self) -> Result<(), Error> {
        if self.band_hz.iter().any(|v| !v.is_finite()) || self.band_hz[0] <= 0.0
            || self.band_hz[1] <= self.band_hz[0] || self.band_hz[1] >= 0.45*f64::from(OUTPUT_RATE)
            || !(2..=observer_fit::MAX_ORDER).contains(&self.max_order) || self.max_order%2 != 0
            || self.training_intervals > observer_fit::MAX_INTERVALS
            || self.training_intervals < 2*self.max_order || self.subdivisions > 4
            || self.max_panels == 0 || self.max_panels > MAX_REFINED_PANELS
            || self.max_dense_work == 0 || self.max_dense_work > MAX_DENSE_WORK
        {
            return Err("radiation preparation requires 0 < low < high < 21600 Hz, even order 2..=32, 2*order..=128 training intervals, 0..=4 subdivisions, 1..=4096 panels and positive bounded dense work".into());
        }
        Ok(())
    }
    pub fn parse(text: &str) -> Result<Self, Error> {
        if text.len() > MAX_BYTES { return Err("radiation specification exceeds 4096 bytes".into()); }
        let mut lines=text.lines().map(|l|l.split('#').next().unwrap_or("").trim()).filter(|l| !l.is_empty());
        if lines.next() != Some(HEADER) { return Err("radiation specification requires frankensim-radiation-preparation-v1".into()); }
        let mut band=None; let mut intervals=None; let mut order=None;
        let mut subdivisions=None; let mut panels=None; let mut work=None;
        for line in lines {
            let f: Vec<_>=line.split(',').map(str::trim).collect();
            match f.as_slice() {
                ["band_hz",a,b] if band.is_none() => band=Some([a.parse()?,b.parse()?]),
                ["training_intervals",v] if intervals.is_none() => intervals=Some(v.parse()?),
                ["max_order",v] if order.is_none() => order=Some(v.parse()?),
                ["subdivisions",v] if subdivisions.is_none() => subdivisions=Some(v.parse()?),
                ["max_panels",v] if panels.is_none() => panels=Some(v.parse()?),
                ["max_dense_work",v] if work.is_none() => work=Some(v.parse()?),
                _ => return Err("radiation specification has an unknown, duplicate or malformed record".into()),
            }
        }
        let s=Self { band_hz:band.ok_or("missing radiation band_hz")?,
            training_intervals:intervals.ok_or("missing radiation training_intervals")?,
            max_order:order.ok_or("missing radiation max_order")?,
            subdivisions:subdivisions.ok_or("missing radiation subdivisions")?,
            max_panels:panels.ok_or("missing radiation max_panels")?,
            max_dense_work:work.ok_or("missing radiation max_dense_work")? };
        s.validate()?; Ok(s)
    }
    pub fn load(path: &str) -> Result<Self, Error> {
        let mut text=String::new();
        std::fs::File::open(path)?.take((MAX_BYTES+1) as u64).read_to_string(&mut text)?;
        Self::parse(&text)
    }
    pub fn admit_command(self, command: &str) -> Result<(), Error> {
        self.validate()?;
        if !matches!(command,"splash-wav"|"splash-mic"|"drum-wav"|"drum-mic"|
            "drum-stretch-wav"|"drum-stretch-mic"|"drum-modal-wav"|"drum-modal-mic"|
            "snare-wav"|"snare-mic"|"snare-off-wav"|"snare-off-mic") {
            return Err("--radiation-spec requires an existing pressure-audio command".into());
        }
        Ok(())
    }
    pub(super) fn frequencies(self) -> Result<Vec<f64>, Error> {
        self.validate()?;
        let steps=4*self.training_intervals;
        let width=self.band_hz[1]-self.band_hz[0];
        let omega: Vec<_>=(0..=steps).map(|i|core::f64::consts::TAU*(self.band_hz[0]+width*i as f64/steps as f64)).collect();
        if omega.windows(2).any(|w|w[1]<=w[0]) {
            return Err("radiation frequency lattice is not representable at the requested spacing".into());
        }
        Ok(omega)
    }
    /// A conservative algebraic work counter, NOT measured flops or elapsed time.
    /// Bound all global refinement and dense frequency solves before allocation.
    pub(super) fn work(self, base_panels: usize, inputs: usize, receivers: usize) -> Result<(usize,u64), Error> {
        self.validate()?;
        if base_panels == 0 || inputs == 0 || inputs > MAX_INPUTS || !(1..=2).contains(&receivers) {
            return Err("radiation work requires nonempty bounded panels, source rows and receivers".into());
        }
        let panels=base_panels.checked_mul(4_usize.pow(self.subdivisions)).ok_or("radiation panel count overflow")?;
        if panels > self.max_panels { return Err("radiation refinement exceeds the declared panel budget".into()); }
        let n=panels as u128; let sources=inputs as u128; let listeners=receivers as u128;
        let units=(4*self.training_intervals+1) as u128*(n*n*n+2*n*n*sources+n*sources*listeners);
        if units > self.max_dense_work as u128 { return Err("radiation preparation exceeds the declared dense-work budget".into()); }
        Ok((panels,units as u64))
    }
    /// A wider observer must not silently exceed the supplied compact neck law.
    pub(super) fn admit_necks(self, experiment: &Experiment) -> Result<(), Error> {
        self.validate()?;
        if let Some(air)=&experiment.air {
            let k=core::f64::consts::TAU*self.band_hz[1]/Medium::air().sound_speed;
            for i in 0..air.coupling.neck_count() {
                let port=air.coupling.neck_radiation_port(i)?;
                let radius=(port.area_m2/core::f64::consts::PI).sqrt();
                if !radius.is_finite() || !port.effective_length_m.is_finite()
                    || k*radius.max(port.effective_length_m)>0.3 {
                    return Err("selected radiation band exceeds the existing compact-neck ka/kL limit".into());
                }
            }
        }
        Ok(())
    }
}

/// Transactional option removal: malformed files leave all arguments intact.
pub fn option(args: &mut Vec<String>) -> Result<Option<Spec>, Error> {
    let mut occurrences=args.iter().enumerate().filter(|(_,a)|a.as_str()=="--radiation-spec");
    let Some((index,_))=occurrences.next() else { return Ok(None); };
    if occurrences.next().is_some() { return Err("--radiation-spec may be supplied only once".into()); }
    let path=args.get(index+1).ok_or("--radiation-spec requires a file path")?;
    let spec=Spec::load(path)?;
    args.drain(index..index+2); Ok(Some(spec))
}

#[cfg(test)]
#[path = "radiation_spec/tests.rs"]
mod tests;
