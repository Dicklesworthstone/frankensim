//! Source geometry to the original bow-position program. This adapter expands
//! ribbon centres with the existing spatial owner; it never interpolates time,
//! alters modal coefficients, or advances a separate mechanical state.
use super::*;
use fs_couple::bowed_string::runtime::reciprocal::bank::multiple::{
    performance::bow_stations::{BowStationKnot, BowStationMotionLimits, BowStationProgram},
    ribbon::UniformBowRibbon,
};

const MAX_KNOTS: usize = 128;
const MAX_TRACKS: usize = 16;
const MAX_COORDINATES: usize = MAX_KNOTS * 63;
const MAX_UPDATE_TERMS: usize = 1_048_576;

/// Source endpoints in authored track order. A ribbon uses one centre per
/// endpoint; its complete footprint is projected before material-pack I/O.
#[derive(Clone, Debug)]
pub struct Input {
    limits: BowStationMotionLimits,
    knots: Vec<BowStationKnot>,
}
impl Input {
    pub(super) fn parse(reader: &mut Reader<'_>) -> Result<Option<Self>> {
        if reader
            .lines
            .clone()
            .next()
            .and_then(|line| line.split_ascii_whitespace().next())
            != Some("bow_station_motion")
        {
            return Ok(None);
        }
        let mut row = reader.row("bow_station_motion")?;
        let limits = BowStationMotionLimits {
            maximum_speed_m_s: row.positive()?,
            maximum_step_m: row.positive()?,
            max_update_terms: row.count(MAX_UPDATE_TERMS)?,
        };
        let count = row.count(MAX_KNOTS)?;
        row.end()?;
        limits.validate()?;
        if count == 0 {
            return Err("bow_station_motion requires at least one endpoint".into());
        }
        let mut knots = Vec::with_capacity(count);
        for _ in 0..count {
            let mut row = reader.row("bow_station_knot")?;
            let sample = row.parse()?;
            let mut station_fractions = Vec::new();
            while row.words.clone().next().is_some() {
                if station_fractions.len() == MAX_TRACKS {
                    return Err("bow station endpoints exceed the 16 authored-track limit".into());
                }
                station_fractions.push(row.scalar()?);
            }
            row.end()?;
            if station_fractions.is_empty() {
                return Err("bow_station_knot needs one station per original bow track".into());
            }
            knots.push(BowStationKnot {
                sample,
                station_fractions,
            });
        }
        Ok(Some(Self { limits, knots }))
    }

    fn program(&self, spec: &Spec, gate: &CancelGate) -> Result<BowStationProgram> {
        poll(gate)?;
        if (spec.contact_thermal && !spec.dry_contact_thermal) || spec.friction.shear_layers().is_some() {
            return Err("bow migration requires transported material shear/thermal patches; those combinations are not supported".into());
        }
        let mut locations = Vec::with_capacity(spec.schedule.tracks().len());
        let mut lengths = Vec::with_capacity(spec.schedule.tracks().len());
        for track in spec.schedule.tracks() {
            let string = fs_scenario::gesture::path::legacy_member(&track.target, "bow")
                .ok_or("bow station motion requires physical bow tracks")?;
            let GestureValue::Bow { station, .. } = track.initial else {
                return Err("bow station motion requires complete initial bow coordinates".into());
            };
            let wire = spec
                .strings
                .get(string)
                .ok_or("bow motion names an absent string")?;
            locations.push(StringBowLocation {
                string,
                station_fraction: station,
            });
            lengths.push(wire.length);
        }
        let original = locations
            .iter()
            .map(|location| location.station_fraction)
            .collect::<Vec<_>>();
        let mut program = BowStationProgram {
            knots: self.knots.clone(),
            limits: self.limits,
            max_knots: MAX_KNOTS,
            max_coordinate_values: MAX_COORDINATES,
        };
        // Validate the authored path first, so missing track coordinates cannot
        // disappear while a one-track ribbon is expanded to its actual ports.
        program.validate(&original, &lengths, spec.samples, spec.rate, gate)?;
        if let Some((width, cells)) = spec.bow_ribbon {
            if locations.len() != 1 {
                return Err("moving ribbon needs exactly one original centre track".into());
            }
            let initial = UniformBowRibbon::try_new(locations[0], lengths[0], width, cells)?;
            for knot in &mut program.knots {
                poll(gate)?;
                let centre = StringBowLocation {
                    string: locations[0].string,
                    station_fraction: knot.station_fractions[0],
                };
                knot.station_fractions =
                    UniformBowRibbon::try_new(centre, lengths[0], width, cells)?
                        .stations()
                        .to_vec();
            }
            program.validate(
                initial.stations(),
                &vec![lengths[0]; cells],
                spec.samples,
                spec.rate,
                gate,
            )?;
        }
        Ok(program)
    }
}

impl Spec {
    /// Recheck source geometry against the current description before pack I/O.
    /// Runtime binding repeats admission against the actual material bank.
    pub fn bow_station_program(&self, gate: &CancelGate) -> Result<Option<BowStationProgram>> {
        self.bow_station_motion
            .as_ref()
            .map(|motion| motion.program(self, gate))
            .transpose()
    }
}
