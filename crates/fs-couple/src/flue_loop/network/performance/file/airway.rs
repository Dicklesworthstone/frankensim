//! The optional `upstream_airway` record: the player's airway ahead of the windway (`flue_loop::airway`).
//!
//! `upstream_airway weibel-1963 saab-2021-table-1-B GENERATIONS SECTIONS SOURCE_TOKEN` directly after `jet_limits` and
//! the optional `jet_noise`, then SECTIONS records `airway_section LENGTH_M AREA_M2` from the windway toward the lungs,
//! then an optional `airway_wall SURFACE_DENSITY STIFFNESS_PER_AREA RESISTANCE`. Each section is the cylinder of its
//! area ([`AirwaySection::from_area`]). After the last one come the first GENERATIONS of Weibel's tree as Saab et al.
//! (2021) print it ([`AirwayChain::weibel_saab_2021_b`]), and the tree's far end holds the gesture's blowing pressure
//! as the lung pressure. Every tube takes Bessel losses in the source gas and carries the wall pin, rigid without the
//! row. The source token names where the duct's areas come from. The airway is realized on the mechanics clock, on the
//! characteristic line's ladder at its default tolerances with a ceiling of [`UPSTREAM_AIRWAY_MAX_TAPS`] taps. Omitting
//! the record leaves every result bit-identical.

use super::{FlueInputError, Rows, bad};
use crate::driving_point::LineLadder;
use crate::upstream_airway::{AirwayChain, AirwaySection};
use crate::voice_performance::{SUBGLOTTAL_AIRWAY_GEOMETRY, SUBGLOTTAL_AIRWAY_TREE};
use fs_duct::LossModel;
use fs_phs::WallPin;

/// Most duct sections one record admits.
pub const MAX_UPSTREAM_AIRWAY_SECTIONS: usize = 256;
/// The ladder's ceiling: the airway's own largest DFT length.
pub const UPSTREAM_AIRWAY_MAX_TAPS: usize = 131_072;

/// A parsed `upstream_airway` record.
#[derive(Debug, Clone, PartialEq)]
pub struct FlueAirwayRecord {
    /// The airway from the windway's inlet to the lungs: the duct's sections, then the tree's generations.
    pub chain: AirwayChain,
    /// How many of the tree's generations follow the duct.
    pub generations: usize,
    /// How many duct sections precede the tree.
    pub duct_sections: usize,
    /// Where the duct's areas come from, as the card names it.
    pub source: String,
    /// The pin every tube carries; `None` is rigid.
    pub wall: Option<WallPin>,
}

/// The ladder every flue airway is realized on.
pub(super) fn ladder() -> LineLadder {
    LineLadder { max_taps: UPSTREAM_AIRWAY_MAX_TAPS, ..LineLadder::default() }
}

/// Read the record if the next row starts it.
pub(super) fn read(rows: &mut Rows<'_>) -> Result<Option<FlueAirwayRecord>, FlueInputError> {
    if !rows.next_is("upstream_airway") {
        return Ok(None);
    }
    let mut row = rows.row("upstream_airway")?;
    let line = row.line;
    if row.fields.next() != Some(SUBGLOTTAL_AIRWAY_TREE) || row.fields.next() != Some(SUBGLOTTAL_AIRWAY_GEOMETRY) {
        return Err(bad(line, "upstream_airway takes `weibel-1963 saab-2021-table-1-B GENERATIONS SECTIONS SOURCE_TOKEN`"));
    }
    let generations: usize = row.number()?;
    let duct_sections: usize = row.number()?;
    let source = row.fields.next()
        .ok_or_else(|| bad(line, "upstream_airway needs an explicit source token"))?
        .to_string();
    row.end()?;
    if !(1..=24).contains(&generations) {
        return Err(bad(line, "upstream_airway generations must be in 1..=24"));
    }
    if !(1..=MAX_UPSTREAM_AIRWAY_SECTIONS).contains(&duct_sections) {
        return Err(bad(line, "upstream_airway takes 1 to 256 duct sections"));
    }
    let mut duct = Vec::with_capacity(duct_sections);
    for _ in 0..duct_sections {
        let mut row = rows.row("airway_section")?;
        let section = (row.scalar()?, row.scalar()?, row.line);
        row.end()?;
        duct.push(section);
    }
    let wall = if rows.next_is("airway_wall") {
        let mut row = rows.row("airway_wall")?;
        let pin = WallPin { surface_density: row.scalar()?, stiffness_per_area: row.scalar()?, resistance: row.scalar()? };
        let at = row.line;
        row.end()?;
        // The wall pin's own admission: a positive mass and stiffness, a resistance that is not negative.
        if !(pin.surface_density > 0.0 && pin.stiffness_per_area > 0.0 && pin.resistance >= 0.0) {
            return Err(bad(at, "airway_wall takes a positive surface density and stiffness and a non-negative resistance"));
        }
        Some(pin)
    } else {
        None
    };
    let sections = duct.into_iter()
        .map(|(length_m, area_m2, at)| AirwaySection::from_area(length_m, area_m2, 1, wall)
            .map_err(|_| bad(at, "airway_section needs a positive finite length and area")))
        .collect::<Result<Vec<_>, _>>()?;
    let tree = AirwayChain::weibel_saab_2021_b(generations, wall)
        .map_err(|_| bad(line, "upstream_airway generations must be in 1..=24"))?;
    let chain = AirwayChain::new(sections, LossModel::Bessel)
        .and_then(|duct| duct.ahead_of(&tree))
        .map_err(|_| bad(line, "upstream_airway's duct and tree do not form one airway"))?;
    Ok(Some(FlueAirwayRecord { chain, generations, duct_sections, source, wall }))
}
