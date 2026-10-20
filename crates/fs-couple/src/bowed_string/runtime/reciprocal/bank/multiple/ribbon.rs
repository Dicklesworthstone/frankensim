//! Spatially resolved friction across a prescribed rigid ribbon.
//!
//! The width is a physical interval along the string, discretized by midpoint
//! quadrature. Each cell has its OWN slip and traction in the existing joint
//! solve; the string velocity is never averaged before applying friction.
//! Every contact uses the original material-derived transverse/torsional map.
//! This is the rigid, uniformly loaded limit, not compliant bow hair, a thermal
//! rosin law or a validated violin. Refine cells, structural modes and time.
//! Finite-width and point contacts need not produce the same motion: see
//! Pitteroff & Woodhouse, Acta Acustica 84 (1998), 744–757.

use super::{MaterialBowContact, MultiBowStringBank, MultiBowStringBankSource};
use super::super::{ReciprocalStringBank, StringBowLocation};
use crate::bowed_string::runtime::reciprocal::{ReciprocalStringError, invalid, poll};
use crate::render::schedule::force::coupled::contact::friction::{
    FrictionDrive, multiple::MultiFrictionConfig,
};
use fs_exec::CancelGate;

/// Frozen uniform line load and its actual spatial contact quadrature.
/// Contact zero is the original centre; other cells follow from left to right.
#[derive(Clone, Debug, PartialEq)]
pub struct UniformBowRibbon {
    centre: StringBowLocation,
    width_m: f64,
    stations: Vec<f64>,
}

impl UniformBowRibbon {
    fn new(
        centre: StringBowLocation,
        length_m: f64,
        width_m: f64,
        cells: usize,
    ) -> Result<Self, ReciprocalStringError> {
        // Odd quadrature retains the already bound physical centre as port zero.
        if !width_m.is_finite() || width_m <= 0.0 || !length_m.is_finite()
            || length_m <= 0.0 || !(1..=63).contains(&cells) || cells % 2 == 0
        {
            return Err(invalid("ribbon needs positive finite width and an odd cell count in 1..=63"));
        }
        let half = 0.5 * width_m / length_m;
        let left = centre.station_fraction - half;
        let right = centre.station_fraction + half;
        if !left.is_finite() || !right.is_finite() || left <= 0.0 || right >= 1.0
            || !(left < centre.station_fraction && centre.station_fraction < right)
        {
            return Err(invalid("the complete finite ribbon must lie strictly inside its actual string"));
        }
        let step = width_m / length_m / cells as f64;
        let mid = cells / 2;
        let mut stations = Vec::with_capacity(cells);
        stations.push(centre.station_fraction);
        let mut previous = left;
        for cell in 0..cells {
            let station = centre.station_fraction + (cell as f64 - mid as f64) * step;
            if !(previous < station && station < right) {
                return Err(invalid("ribbon cells must have distinct representable interior stations"));
            }
            previous = station;
            if cell != mid { stations.push(station); }
        }
        Ok(Self { centre, width_m, stations })
    }

    /// Actual material string and fixed ribbon centre.
    #[must_use]
    pub const fn centre(&self) -> StringBowLocation { self.centre }
    /// Physical contact width along that string [m], not a tonal parameter.
    #[must_use]
    pub const fn width_m(&self) -> f64 { self.width_m }
    /// Actual contact stations, in the joint solver's port order.
    #[must_use]
    pub fn stations(&self) -> &[f64] { &self.stations }
    /// Equal midpoint-quadrature share of the total normal force.
    #[must_use]
    pub fn load_fraction(&self) -> f64 { 1.0 / self.stations.len() as f64 }

    /// Convert ONE total ribbon load to the original solver's per-cell drives.
    /// The prescribed speed is common; each cell still solves its own slip.
    /// Output storage is caller-owned and must cover every spatial cell.
    pub fn project_drive(
        &self, drive: FrictionDrive, output: &mut [FrictionDrive],
    ) -> Result<(), ReciprocalStringError> {
        if output.len() != self.stations.len() || !drive.speed_m_s.is_finite()
            || !drive.normal_force_n.is_finite() || drive.normal_force_n < 0.0
        {
            return Err(invalid("ribbon drive requires finite speed, nonnegative total load and one output per cell"));
        }
        output.fill(FrictionDrive { speed_m_s: drive.speed_m_s,
            normal_force_n: drive.normal_force_n * self.load_fraction() });
        Ok(())
    }
}

impl MultiBowStringBankSource {
    /// Finite-width source when this bank was constructed from a uniform ribbon.
    #[must_use]
    pub const fn ribbon(&self) -> Option<&UniformBowRibbon> { self.ribbon.as_ref() }
}

impl ReciprocalStringBank {
    /// Resolve a uniform finite-width ribbon on the existing bow's string.
    ///
    /// The original material/body/cavity state, clock and observer survive.
    /// Every new point port participates in the SAME simultaneous implicit
    /// friction solve. The width does not filter sound or prescribe oscillation.
    /// The supplied odd cell count is a numerical resolution, not a hair count.
    ///
    /// An already bound interface refuses: its nominal contact area belongs to
    /// the whole original contact. Bind explicit per-cell interface areas after
    /// this construction with `with_bow_interface`, retaining the original law.
    /// Hair compliance and a mechanically resolved normal-load distribution are
    /// separate physics; this method prescribes uniform loading and rigid speed.
    pub fn with_uniform_bow_ribbon(
        self, width_m: f64, cells: usize, joint: MultiFrictionConfig, gate: &CancelGate,
    ) -> Result<MultiBowStringBank, ReciprocalStringError> {
        poll(gate)?;
        if self.bow_interface().is_some() {
            return Err(invalid("bind ribbon interfaces per cell after resolving the contact area"));
        }
        if cells > joint.max_ports {
            return Err(invalid("ribbon cells exceed the declared joint contact budget"));
        }
        let centre = self.source.options.bow;
        let ribbon = UniformBowRibbon::new(centre, self.cards[centre.string].length_m, width_m, cells)?;
        let law = self.system.port().law.clone();
        let friction = self.source.options.friction;
        let additional = ribbon.stations[1..].iter().map(|&station_fraction| MaterialBowContact {
            location: StringBowLocation { string: centre.string, station_fraction },
            law: law.clone(), friction,
        }).collect();
        let mut bank = self.with_additional_bows(additional, joint, gate)?;
        bank.source.ribbon = Some(ribbon);
        Ok(bank)
    }
}
