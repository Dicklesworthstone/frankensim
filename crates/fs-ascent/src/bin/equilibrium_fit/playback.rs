//! Explicit forward audition of the final, independently rechecked fitted model.
//! This is not an acoustic fitting objective and does not synthesize transfers.
use super::{DesignControl, EquilibriumDesignFile, Limits, json_string, output, solve};
use fs_couple::pcm_wav::stream::{Pcm16WavStream, ScheduledWavProgress, render_scheduled_pcm16};
use fs_couple::render::schedule::force::ForceRenderConfig;
use fs_couple::render::schedule::force::coupled::equilibrium::sensitivity::objective::design::forward::playback::{
    CaseForceEvent, DesignPlaybackConfig,
};
use fs_exec::CancelGate;
use std::fmt::Write as _;
use std::path::Path;

#[derive(Default)]
pub(super) struct Builder {
    wav: Option<String>, case: Option<String>, samples: Option<u64>,
    release: Option<u64>, scale: Option<f64>, block: Option<usize>,
}
pub(super) struct PlaybackOptions {
    wav: String, case: String, samples: u64, release: u64, scale: f64, block: usize,
}
impl Builder {
    pub(super) fn set(&mut self, flag: &str, value: &str) -> Result<(), String> {
        match flag {
            "--playback-wav" => self.wav = Some(value.to_owned()),
            "--playback-case" => self.case = Some(value.to_owned()),
            "--playback-samples" => self.samples = Some(value.parse::<u64>().ok()
                .filter(|x| (1..=28_800_000).contains(x)).ok_or("playback samples must be in 1..=28800000")?),
            "--playback-release" => self.release = Some(value.parse().map_err(|_| "playback release needs a sample ordinal")?),
            "--playback-full-scale-pa" => self.scale = Some(value.parse::<f64>().ok()
                .filter(|x| x.is_finite() && *x > 0.0).ok_or("playback full-scale must be finite positive pascals")?),
            "--playback-block" => self.block = Some(value.parse::<usize>().ok()
                .filter(|x| (1..=65_536).contains(x)).ok_or("playback block must be in 1..=65536")?),
            _ => return Err(format!("unknown playback option {flag}")),
        }
        Ok(())
    }
    pub(super) fn finish(self) -> Result<Option<PlaybackOptions>, String> {
        if self.wav.is_none() && self.case.is_none() && self.samples.is_none()
            && self.release.is_none() && self.scale.is_none() && self.block.is_none() { return Ok(None); }
        let options = PlaybackOptions {
            wav: self.wav.ok_or("playback needs --playback-wav")?,
            case: self.case.ok_or("playback needs --playback-case")?,
            samples: self.samples.ok_or("playback needs --playback-samples")?,
            release: self.release.ok_or("playback needs --playback-release")?,
            scale: self.scale.ok_or("playback needs --playback-full-scale-pa")?,
            block: self.block.unwrap_or(512),
        };
        if options.release >= options.samples { return Err("playback release must lie inside the requested horizon".into()); }
        Ok(Some(options))
    }
}

pub(super) fn fit_and_render(loaded: &EquilibriumDesignFile, limits: Limits,
    options: &PlaybackOptions, gate: &CancelGate) -> Result<String, Box<dyn std::error::Error>>
{
    if limits.evaluations < 3 { return Err("playback needs evaluations for initialization, final audit and its preload".into()); }
    if Path::new(&options.wav).exists() { return Err("playback output already exists; refusing to replace it".into()); }
    let case = loaded.problem().load_cases().iter().position(|case| case.name == options.case)
        .ok_or("unknown playback case name")?;
    let loads = loaded.problem().load_cases()[case].loads.len();
    if loads == 0 { return Err("release playback requires an explicitly loaded experiment".into()); }
    // Reserve ONE physical candidate/preload within the original total limit.
    // The existing fit retains its own independent final family re-solve.
    let mut result = solve(loaded, Limits { evaluations: limits.evaluations-1, ..limits }, gate)?;
    let mut work = DesignControl::new(1, loaded.problem().load_cases().len());
    let events = (0..loads).map(|load| CaseForceEvent { sample: options.release, load, force_n: 0.0 }).collect();
    let playback = loaded.problem().playback_case(&result.report.solution.x, case, events,
        DesignPlaybackConfig { samples: options.samples, force: ForceRenderConfig {
            sample_rate_hz: loaded.model_info().sample_rate_hz, max_block: options.block,
            max_events: loads, max_controls: 262_144, max_projection_terms: 16_777_216,
        } }, &mut work, gate)?;
    let origin = playback.info().clone();
    if origin.physical_parameters != result.audited.physical_parameters {
        return Err("playback parameters differ from the audited accepted design".into());
    }
    result.work.evaluations += work.work().evaluations;
    result.work.case_solves += work.work().case_solves;
    let mut renderer = playback.into_renderer();
    let mut scratch = vec![0.0; options.block];
    if gate.is_requested() { return Err("playback cancelled before output creation".into()); }
    // Exclusive creation also closes the race after the early existence check.
    // A failure after this point can leave an incomplete WAV, never a success JSON.
    let file = std::fs::OpenOptions::new().write(true).create_new(true).open(&options.wav)?;
    let mut stream = Pcm16WavStream::new(file, origin.sample_rate_hz, options.scale, options.block)?;
    let progress = render_scheduled_pcm16(&mut renderer, &mut stream, gate, &mut scratch, options.samples)
        .map_err(|e| format!("playback failed; WAV may be incomplete: {e}"))?;
    if progress != (ScheduledWavProgress::Completed { samples: options.samples }) {
        return Err("playback cancelled; WAV is not finalized".into());
    }
    let (_, summary) = stream.finish()?;
    let mut text = output(loaded, &result);
    text.pop(); // existing complete object's final brace
    write!(&mut text, ",\"playback\":{{\"wav\":{},\"case\":{},\"scope\":\"forward-loaded-release-not-transient-fit\",\"sample_rate_hz\":{},\"samples\":{},\"release_before_sample\":{},\"block\":{},\"full_scale_pa\":{:.17e},\"clipped_samples\":{},\"initial_energy_j\":{:.17e},\"preload_evaluations\":{},\"preload_case_solves\":{}}}}}",
        json_string(&options.wav), json_string(&options.case), summary.sample_rate_hz, summary.samples,
        options.release, options.block, summary.full_scale_pa, summary.clipped_samples,
        origin.initial_energy_j, work.work().evaluations, work.work().case_solves).expect("String write");
    Ok(text)
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn playback_settings_are_explicit_complete_and_in_window() {
        assert!(Builder::default().finish().unwrap().is_none());
        let mut b=Builder::default();b.set("--playback-wav","out.wav").unwrap();assert!(b.finish().is_err());
        for (key,value) in [("--playback-samples","0"),("--playback-full-scale-pa","NaN"),
            ("--playback-block","65537"),("--playback-release","-1"),("--playback-typo","1")] {
            assert!(Builder::default().set(key,value).is_err());
        }
        let mut b=Builder::default();
        for (key,value) in [("--playback-wav","out.wav"),("--playback-case","experiment"),
            ("--playback-samples","37"),("--playback-release","37"),("--playback-full-scale-pa","1")] {
            b.set(key,value).unwrap();
        }
        assert!(b.finish().is_err());
    }
}
