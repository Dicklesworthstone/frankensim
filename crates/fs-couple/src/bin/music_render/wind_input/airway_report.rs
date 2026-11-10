//! The `upstream_airway` sidecar record of a flue card that selects one: the record and the realized airway's receipt,
//! read before the render, and the samples the windway passed no flow (`flue_loop::airway`, a one-way windway), read
//! from the voice after every callback. A card without the row gets no record and its renderer unwrapped.
use std::cell::Cell;
use std::rc::Rc;
use fs_couple::flue_loop::airway::FlueAirwayClosures;
use fs_couple::flue_loop::network::performance::FlueOutletPerformance;
use fs_couple::flue_loop::network::performance::file::FlueFilePerformance;
use fs_couple::pcm_wav::observation::PressureRenderer;
use fs_couple::render::{ControlDelta, RenderError};
use super::json_string;

/// The record's fields known before the render, and the closures the tap keeps.
pub(super) struct AirwayReport {
    head: String,
    seen: Rc<Cell<Option<FlueAirwayClosures>>>,
}

/// The flue renderer, unchanged, with the voice's closures read after every callback.
struct AirwayClosureTap {
    inner: FlueOutletPerformance,
    seen: Rc<Cell<Option<FlueAirwayClosures>>>,
}

/// The record and receipt of `p`'s airway, before it is turned into its renderer; `None` without the row. A line
/// admitted by the passivity gain (MB 23083) names it after `clock_hz`; a line admitted any other way writes the same
/// bytes as before.
pub(super) fn head(p: &FlueFilePerformance) -> Option<String> {
    let record = p.upstream_airway()?;
    let airway = p.renderer().voice().upstream_airway()?;
    let r = airway.receipt();
    let gain = r.reflectance.passivity_gain.map_or_else(String::new, |gain| {
        format!(",\"passivity_gain\":{{\"passes\":{},\"top_band_db\":{:e}}}", gain.passes, gain.top_band_db)
    });
    Some(format!(
        "\"source\":{},\"generations\":{},\"duct_sections\":{},\"sections\":{},\"wall\":{},\"two_way\":{},\"port_impedance_pa_s_m3\":{:e},\"dc_resistance_pa_s_m3\":{:e},\"same_sample_reflection\":{:e},\"reflectance_taps\":{},\"source_taps\":{},\"clock_hz\":{}{gain}",
        json_string(&record.source), record.generations, record.duct_sections, record.chain.sections().len(),
        record.wall.is_some(), record.two_way,
        r.port_impedance, r.dc_resistance, r.same_sample_reflection, r.reflectance.n_fft, r.source_taps,
        airway.sample_rate_hz()))
}

/// `renderer` behind the tap and its report when `head` is a record; `renderer` alone otherwise.
pub(super) fn tap(renderer: FlueOutletPerformance, head: Option<String>)
    -> (Box<dyn PressureRenderer>, Option<AirwayReport>)
{
    let Some(head) = head else { return (Box::new(renderer), None) };
    let seen = Rc::new(Cell::new(renderer.voice().upstream_airway_closures()));
    (Box::new(AirwayClosureTap { inner: renderer, seen: Rc::clone(&seen) }), Some(AirwayReport { head, seen }))
}

impl AirwayReport {
    /// The sidecar member, with its trailing comma: indices count the voice's own clock, `clock_hz`, from its first
    /// accepted sample.
    pub(super) fn json(&self) -> String {
        let index = |i: Option<u64>| i.map_or_else(|| "null".to_owned(), |i| i.to_string());
        let no_flow = self.seen.get().map_or_else(|| "null".to_owned(), |c| format!(
            "{{\"samples\":{},\"runs\":{},\"first_sample\":{},\"last_sample\":{},\"least_open_circuit_pa\":{:e},\"arrival_samples\":{:e},\"lung_before_arrival_pa\":{:e},\"least_before_arrival_pa\":{:e},\"samples_from_arrival\":{},\"first_sample_from_arrival\":{},\"least_from_arrival_pa\":{:e},\"least_flow_m3_s\":{:e},\"scope\":\"accepted samples whose open-circuit pressure at the windway was not positive: no jet; one way (two_way false) the windway passed no flow and the airway saw a closed end, two way it drew air back (least_flow_m3_s); kept apart before the source FIR's arrival (arrival_samples, where its band-limited pre-ringing of a sub-sample arrival lies) and from it\"}}",
            c.samples, c.runs, index(c.first_sample), index(c.last_sample), c.least_open_circuit_pa, c.arrival_samples,
            c.lung_before_arrival_pa, c.least_before_arrival_pa, c.samples_from_arrival, index(c.first_sample_from_arrival),
            c.least_from_arrival_pa, c.least_flow_m3_s));
        format!("\"upstream_airway\":{{{},\"no_flow\":{no_flow}}},", self.head)
    }
}

impl PressureRenderer for AirwayClosureTap {
    fn samples_rendered(&self) -> u64 { self.inner.samples_rendered() }
    fn max_block_len(&self) -> usize { self.inner.max_block_len() }
    fn validate_sample_rate(&self, sample_rate_hz: u32) -> Result<(), RenderError> {
        self.inner.validate_sample_rate(sample_rate_hz)
    }
    fn validate_sample_count(&self, samples: u64) -> Result<(), RenderError> {
        self.inner.validate_sample_count(samples)
    }
    fn validate_controls(&self, deltas: &[ControlDelta]) -> Result<(), RenderError> {
        self.inner.validate_controls(deltas)
    }
    fn apply_controls(&mut self, deltas: &[ControlDelta]) -> Result<(), RenderError> {
        self.inner.apply_controls(deltas)
    }
    fn block(&mut self, output: &mut [f64]) -> Result<(), RenderError> {
        self.inner.block(output)?;
        self.seen.set(self.inner.voice().upstream_airway_closures());
        Ok(())
    }
}
