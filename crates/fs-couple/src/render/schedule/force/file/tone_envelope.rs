//! Sample-clock gain knots for physical tone forces, never output amplitude.
use super::Reader;

const MAX_ENVELOPE_POINTS: usize = 65_536;

struct GainKnot {
    sample: u64,
    gain: f64,
}

/// One immutable knot sequence and its accepted interval cursor.
pub(super) struct GainEnvelope {
    points: Vec<GainKnot>,
    interval: usize,
}

/// Original tone-row indexing, before any component grouping or equal-frequency
/// reduction. Undeclared tones retain unit gain and carrier phase.
pub(super) struct ToneEnvelopePlan {
    curves: Vec<Option<GainEnvelope>>,
    count: usize,
}

impl ToneEnvelopePlan {
    /// Optional `envelopes COUNT`, then unique `envelope TONE POINTS` and
    /// `gain SAMPLE GAIN` rows. First knot is zero, samples strictly increase
    /// through the declared horizon, and the final gain is held. Signed gains
    /// are finite physical force multipliers. Zero envelopes retain the old path.
    /// Declared payload is charged before allocation; caller bytes and allocator
    /// metadata keep their existing separate bounds.
    pub(super) fn read_optional(
        reader: &mut Reader<'_>,
        tone_count: usize,
        samples: u64,
        remaining_storage: &mut usize,
    ) -> Result<Option<Self>, String> {
        if reader
            .lines
            .clone()
            .next()
            .and_then(|line| line.split_ascii_whitespace().next())
            != Some("envelopes")
        {
            return Ok(None);
        }
        let count: usize = reader.one("envelopes").map_err(|e| e.to_string())?;
        if count > tone_count || samples == 0 {
            return Err(
                "envelopes require existing tone rows and a positive performance horizon".into(),
            );
        }
        if count == 0 {
            return Ok(None);
        }
        let metadata = tone_count
            .checked_mul(core::mem::size_of::<Option<GainEnvelope>>())
            .and_then(|bytes| bytes.checked_add(core::mem::size_of::<Self>()))
            .ok_or("envelope metadata size overflow")?;
        reserve(remaining_storage, metadata)?;
        let mut curves: Vec<Option<GainEnvelope>> = (0..tone_count).map(|_| None).collect();
        let mut total_points = 0usize;
        for _ in 0..count {
            let mut row = reader.row("envelope").map_err(|e| e.to_string())?;
            let tone: usize = row.parse().map_err(|e| e.to_string())?;
            let points = row
                .count(MAX_ENVELOPE_POINTS - total_points)
                .map_err(|e| e.to_string())?;
            row.finish().map_err(|e| e.to_string())?;
            if tone >= tone_count || curves[tone].is_some() || points == 0 {
                return Err(
                    "each envelope requires one unique original tone and nonempty gain knots"
                        .into(),
                );
            }
            reserve(
                remaining_storage,
                points
                    .checked_mul(core::mem::size_of::<GainKnot>())
                    .ok_or("envelope knot size overflow")?,
            )?;
            total_points += points;
            let mut knots: Vec<GainKnot> = Vec::with_capacity(points);
            for point in 0..points {
                let mut row = reader.row("gain").map_err(|e| e.to_string())?;
                let sample: u64 = row.parse().map_err(|e| e.to_string())?;
                let gain = row.scalar().map_err(|e| e.to_string())?;
                row.finish().map_err(|e| e.to_string())?;
                if sample > samples
                    || (point == 0 && sample != 0)
                    || knots
                        .last()
                        .is_some_and(|previous| previous.sample >= sample)
                {
                    return Err("gain knots must start at sample zero and strictly increase within the performance horizon".into());
                }
                knots.push(GainKnot { sample, gain });
            }
            curves[tone] = Some(GainEnvelope {
                points: knots,
                interval: 0,
            });
        }
        Ok(Some(Self { curves, count }))
    }

    pub(super) const fn count(&self) -> usize {
        self.count
    }

    pub(super) fn contains(&self, tone: usize) -> bool {
        self.curves[tone].is_some()
    }

    /// Whole-curve amplitude and change-per-sample bounds, independent of the
    /// accepted cursor. Used to admit frequency reconstruction before rendering.
    pub(super) fn bounds(&self, tone: usize) -> Result<[f64; 2], String> {
        match self
            .curves
            .get(tone)
            .ok_or("envelope tone index is outside its original source")?
        {
            Some(curve) => curve.bounds(),
            None => Ok([1.0, 0.0]),
        }
    }

    /// Gains at this sample's start and end. Both multiply coefficients in the
    /// SAME interval-start carrier basis; the end gain does not rotate phase.
    pub(super) fn endpoints(&self, tone: usize, sample: u64) -> Result<[f64; 2], String> {
        match self
            .curves
            .get(tone)
            .ok_or("envelope tone index is outside its original source")?
        {
            Some(curve) => curve.endpoints(sample),
            None => Ok([1.0, 1.0]),
        }
    }

    /// Publish cursor changes only with an accepted physical sample.
    pub(super) fn commit(&mut self, next_sample: u64) {
        for curve in self.curves.iter_mut().flatten() {
            curve.commit(next_sample);
        }
    }

    /// Move a curve into a component driver without cloning its knot payload.
    pub(super) fn take(&mut self, tone: usize) -> Option<GainEnvelope> {
        self.curves[tone].take()
    }
}

impl GainEnvelope {
    pub(super) fn bounds(&self) -> Result<[f64; 2], String> {
        let maximum_gain = self
            .points
            .iter()
            .map(|point| point.gain.abs())
            .fold(0.0, f64::max);
        let mut maximum_change = 0.0_f64;
        for pair in self.points.windows(2) {
            let samples = (pair[1].sample - pair[0].sample) as f64;
            // Divide before subtracting so widely separated opposite finite
            // gains can have a representable per-sample change.
            let change = (pair[1].gain / samples - pair[0].gain / samples).abs();
            if !change.is_finite() {
                return Err("physical gain change per sample is not representable".into());
            }
            maximum_change = maximum_change.max(change);
        }
        Ok([maximum_gain, maximum_change])
    }

    pub(super) fn endpoints(&self, sample: u64) -> Result<[f64; 2], String> {
        let end = sample
            .checked_add(1)
            .ok_or("envelope sample clock overflow")?;
        let left = &self.points[self.interval];
        let Some(right) = self.points.get(self.interval + 1) else {
            return Ok([left.gain, left.gain]);
        };
        if sample < left.sample || end > right.sample {
            return Err("envelope query does not match its accepted sample cursor".into());
        }
        let value = |at| {
            if at == left.sample {
                return left.gain;
            }
            if at == right.sample {
                return right.gain;
            }
            let t = (at - left.sample) as f64 / (right.sample - left.sample) as f64;
            // Convex form avoids overflow of right-left for opposite large
            // finite gains and preserves exact knot/zero endpoints.
            left.gain * (1.0 - t) + right.gain * t
        };
        let gains = [value(sample), value(end)];
        if gains.iter().any(|g| !g.is_finite()) {
            return Err("interpolated physical gain is not representable".into());
        }
        Ok(gains)
    }

    pub(super) fn commit(&mut self, next_sample: u64) {
        while self
            .points
            .get(self.interval + 1)
            .is_some_and(|point| point.sample <= next_sample)
        {
            self.interval += 1;
        }
    }
}

fn reserve(remaining: &mut usize, bytes: usize) -> Result<(), String> {
    *remaining = remaining
        .checked_sub(bytes)
        .ok_or("physical gain envelopes exceed source storage")?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn read(text: &str, storage: usize) -> Result<Option<ToneEnvelopePlan>, String> {
        ToneEnvelopePlan::read_optional(
            &mut Reader {
                lines: text.lines(),
                line: 0,
            },
            2,
            16,
            &mut { storage },
        )
    }

    #[test]
    fn knot_endpoints_default_gain_and_accepted_cursor_preserve_the_physical_clock() {
        let mut plan = read(
            "envelopes 1\nenvelope 1 4\ngain 0 0\ngain 3 0\ngain 7 2\ngain 11 -2\n",
            4096,
        )
        .unwrap()
        .unwrap();
        assert_eq!(plan.count(), 1);
        for sample in 0..16 {
            assert_eq!(plan.endpoints(0, sample).unwrap(), [1.0, 1.0]);
            let value = |at| {
                if at <= 3 {
                    0.0
                } else if at <= 7 {
                    (at - 3) as f64 * 0.5
                } else if at <= 11 {
                    2.0 - (at - 7) as f64
                } else {
                    -2.0
                }
            };
            let expected = [value(sample), value(sample + 1)];
            assert_eq!(plan.endpoints(1, sample).unwrap(), expected);
            // Repeated candidate queries publish no interval progress.
            assert_eq!(plan.endpoints(1, sample).unwrap(), expected);
            plan.commit(sample + 1);
        }
        assert!(read("envelopes 0\n", 0).unwrap().is_none());
        let mut reader = Reader {
            lines: "broadband 1\n".lines(),
            line: 0,
        };
        assert!(
            ToneEnvelopePlan::read_optional(&mut reader, 2, 16, &mut 0)
                .unwrap()
                .is_none()
        );
        assert_eq!(reader.lines.next(), Some("broadband 1"));
    }

    #[test]
    fn malformed_or_unbudgeted_envelopes_refuse_before_count_driven_allocation() {
        let valid = "envelopes 1\nenvelope 1 2\ngain 0 0\ngain 16 1\n";
        assert!(read(valid, 1).is_err());
        for bad in [
            valid.replace("envelopes 1", "envelopes 3"),
            valid.replace("envelope 1 2", "envelope 2 2"),
            valid.replace("envelope 1 2", "envelope 1 0"),
            valid.replace("envelope 1 2", "envelope 1 65537"),
            valid.replace("gain 0 0", "gain 1 0"),
            valid.replace("gain 16 1", "gain 0 1"),
            valid.replace("gain 16 1", "gain 17 1"),
            valid.replace("gain 16 1", "gain 16 NaN"),
            valid.replace("gain 16 1", "gain 16 1 extra"),
            "envelopes 2\nenvelope 1 1\ngain 0 1\nenvelope 1 1\ngain 0 1\n".into(),
        ] {
            assert!(read(&bad, 4096).is_err(), "admitted {bad}");
        }
    }
}
