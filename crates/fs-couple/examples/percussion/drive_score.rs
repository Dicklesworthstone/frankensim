//! Cold compilation of a rhythmic score into the original SI force program.
//! A stroke starts an authored HAND force, not a collision or a velocity reset.
//! Contact time, rebound and acoustic output remain results of the mechanics.
use super::{Error, Knot, Program, MAX_BYTES, MAX_KNOTS};
use std::collections::BTreeMap;

pub(super) const HEADER: &str = "frankensim-stick-score-v1";
const MAX_SHAPES: usize = 128;
const MAX_TEMPOS: usize = 1024;
const MAX_STROKES: usize = MAX_KNOTS / 2;
const MAX_OVERLAP_TERMS: usize = 2_000_000;

fn rows(text: &str) -> impl Iterator<Item = (usize, &str)> {
    text.lines().enumerate().filter_map(|(i, raw)| {
        let row = raw.split('#').next().unwrap_or("").trim();
        (!row.is_empty()).then_some((i + 1, row))
    })
}
pub(super) fn selected(text: &str) -> bool {
    rows(text).next().is_some_and(|(_, row)| row == HEADER)
}
struct Tempo { beat: f64, time_s: f64, seconds_per_beat: f64 }
struct Stroke<'a> { beat: f64, shape: &'a str, scale: f64 }

pub(super) fn parse(text: &str) -> Result<Program, Error> {
    if text.len() as u64 > MAX_BYTES { return Err("stick score exceeds 4 MiB".into()); }
    let mut rows = rows(text);
    if rows.next().map(|(_, row)| row) != Some(HEADER) {
        return Err("missing stick-score header".into());
    }
    let mut tempos: Vec<Tempo> = Vec::new();
    let mut shapes: BTreeMap<&str, Vec<Knot>> = BTreeMap::new();
    let mut strokes: Vec<Stroke<'_>> = Vec::new();
    let mut shape_knots = 0;
    for (line, row) in rows {
        let fields: Vec<_> = row.split(',').map(str::trim).take(7).collect();
        let bad = || format!("stick-score line {line}: invalid tempo, shape, stroke or roll record");
        let number = |i: usize| fields[i].parse::<f64>().map_err(|_| bad());
        match (fields[0], fields.len()) {
            ("tempo", 3) => {
                let beat = number(1)?; let bpm = number(2)?;
                let seconds_per_beat = 60.0 / bpm;
                if tempos.len() == MAX_TEMPOS || !beat.is_finite() || beat < 0.0
                    || !bpm.is_finite() || bpm <= 0.0
                    || !seconds_per_beat.is_finite() || seconds_per_beat <= 0.0 {
                    return Err(bad().into());
                }
                let time_s = if let Some(previous) = tempos.last() {
                    let time = previous.time_s + (beat - previous.beat) * previous.seconds_per_beat;
                    if beat <= previous.beat || !time.is_finite() || time <= previous.time_s {
                        return Err("tempo boundaries must advance in both beats and representable seconds".into());
                    }
                    time
                } else {
                    if beat != 0.0 { return Err("the first explicit tempo must be at beat zero".into()); }
                    0.0
                };
                tempos.push(Tempo { beat, time_s, seconds_per_beat });
            }
            ("shape", 4) => {
                let name = fields[1];
                if name.is_empty() || name.len() > 64
                    || !name.bytes().all(|c| c.is_ascii_alphanumeric() || c == b'_' || c == b'-')
                    || (!shapes.contains_key(name) && shapes.len() == MAX_SHAPES)
                    || shape_knots == MAX_KNOTS {
                    return Err(bad().into());
                }
                let time_s = number(2)?; let force_n = number(3)?;
                let shape = shapes.entry(name).or_default();
                if !time_s.is_finite() || time_s < 0.0 || !force_n.is_finite()
                    || shape.last().is_some_and(|k| time_s <= k.time_s) {
                    return Err(bad().into());
                }
                shape.push(Knot { time_s, force_n }); shape_knots += 1;
            }
            ("stroke", 4) | ("roll", 6) => {
                let beat = number(1)?;
                let (spacing, count, name, scale) = if fields[0] == "roll" {
                    (number(2)?, fields[3].parse::<usize>().map_err(|_| bad())?, fields[4], number(5)?)
                } else { (0.0, 1, fields[2], number(3)?) };
                if !beat.is_finite() || beat < 0.0 || !scale.is_finite() || scale < 0.0
                    || !spacing.is_finite() || spacing < 0.0
                    || (fields[0] == "roll" && spacing == 0.0) || count == 0
                    || count > MAX_STROKES - strokes.len() {
                    return Err(bad().into());
                }
                let mut previous = None;
                for i in 0..count {
                    let at = beat + i as f64 * spacing;
                    if !at.is_finite() || previous.is_some_and(|p| at <= p) {
                        return Err("roll loses representable beat spacing".into());
                    }
                    strokes.push(Stroke { beat: at, shape: name, scale }); previous = Some(at);
                }
            }
            _ => return Err(bad().into()),
        }
    }
    if tempos.is_empty() || shapes.is_empty() || strokes.is_empty() {
        return Err("stick score needs explicit tempo, force shapes and at least one stroke".into());
    }
    for shape in shapes.values() {
        if shape.len() < 2 || shape[0].time_s != 0.0 || shape[0].force_n != 0.0
            || shape.last().unwrap().force_n != 0.0 {
            return Err("each force shape starts at zero seconds/force and ends at zero force".into());
        }
    }

    // Expand only admitted gestures. They stay in physical SECONDS across a
    // tempo change; tempo changes onset spacing, never the hand's force law.
    strokes.sort_by(|a, b| a.beat.total_cmp(&b.beat));
    let mut previous_onset: Option<(f64, f64)> = None;
    let mut expanded: Vec<Knot> = Vec::new();
    let mut spans = Vec::new();
    for stroke in strokes {
        let shape = shapes.get(stroke.shape)
            .ok_or_else(|| format!("undefined force shape {}", stroke.shape))?;
        if shape.len() > MAX_KNOTS - expanded.len() {
            return Err("expanded stick score exceeds 65536 force knots".into());
        }
        let tempo = &tempos[tempos.partition_point(|t| t.beat <= stroke.beat) - 1];
        let start = tempo.time_s + (stroke.beat - tempo.beat) * tempo.seconds_per_beat;
        if !start.is_finite() || (stroke.beat > tempo.beat && start <= tempo.time_s) {
            return Err("stroke onset loses representable physical time".into());
        }
        if previous_onset.is_some_and(|(beat, time)| stroke.beat > beat && start <= time) {
            return Err("distinct score beats collapse to the same physical onset".into());
        }
        previous_onset = Some((stroke.beat, start));
        let begin = expanded.len();
        for knot in shape {
            let time_s = start + knot.time_s;
            let force_n = stroke.scale * knot.force_n;
            if !time_s.is_finite() || !force_n.is_finite()
                || (expanded.len() > begin && time_s <= expanded.last().unwrap().time_s) {
                return Err("shifted force shape overflows or loses temporal resolution".into());
            }
            expanded.push(Knot { time_s, force_n });
        }
        spans.push(begin..expanded.len());
    }

    // Superpose at the UNION of all breakpoints. Concatenating gestures would
    // lose overlaps; slope-prefix sums would accumulate a spurious force tail.
    // Local interpolation also gives exact zero on gaps and at outer endpoints.
    let mut times: Vec<_> = expanded.iter().map(|k| k.time_s).collect();
    times.sort_by(f64::total_cmp); times.dedup();
    let mut knots: Vec<_> = times.into_iter().map(|time_s| Knot { time_s, force_n: 0.0 }).collect();
    let mut work = 0;
    for span in spans {
        let shape = &expanded[span];
        let begin = knots.partition_point(|k| k.time_s < shape[0].time_s);
        let end = knots.partition_point(|k| k.time_s <= shape.last().unwrap().time_s);
        if end - begin > MAX_OVERLAP_TERMS - work {
            return Err("stick score exceeds bounded overlap-compilation work".into());
        }
        work += end - begin;
        let mut segment = 0;
        for knot in &mut knots[begin..end] {
            while segment + 1 < shape.len() - 1 && shape[segment + 1].time_s < knot.time_s {
                segment += 1;
            }
            let a = shape[segment]; let b = shape[segment + 1];
            let x = (knot.time_s - a.time_s) / (b.time_s - a.time_s);
            knot.force_n += (1.0 - x) * a.force_n + x * b.force_n;
            if !knot.force_n.is_finite() { return Err("overlapping hand forces overflow".into()); }
        }
    }
    // The existing Program admits the COMPLETE resulting force and duration
    // against the actual mass/clock before any step. No envelope is weakened.
    Ok(Program { knots })
}
