//! A closed-open air column's plane-wave modes, and the column as a composition-kernel Component.
//!
//! The law is a copy of the resonator-tube law of fs-couple's `bar_strike` example
//! (`examples/bar_strike/tube.rs`: `open_end`, `isolated_modes`, `cylinders`, the drawn-bore walk and
//! `acoustic_mass`), operation for operation, so that both give the same bits on the same column; a
//! test in that example compares the two on every tube card. The example keeps its own copy.
//!
//! The column's frequency and loss are Doutaut, Matignon and Chaigne, JASA 104 (1998): Eq. 7 for the
//! viscothermal wall, Eq. 12 (Levine and Schwinger, unflanged) for the open end. Each mode's frequency
//! is a fixed point, since the wall slows the wave and the open end's correction both depend on it.
//!
//! With a drawn bore each mode is the same fixed point with the column's length replaced by its
//! plane-wave length in the mouth's units: a walk up the bore, each cylinder's phase at its own Eq. 7
//! speed, and at each change of radius the pressure and the volume velocity continuous. The open
//! end's correction is Eq. 12's at the mouth. Each mode's mass per unit volume velocity, its wall loss
//! (Eq. 7's viscous term on the mode's kinetic energy, its thermal term on the potential energy, each
//! cylinder at its own radius) and its open-end loss come from the mode's energy along the bore. Drawn
//! as one cylinder, the bore gives the uniform law bit for bit.
//!
//! [`TubeResonator`] holds a column's modes with their acoustic masses, the data a joint modal basis
//! takes, and [`TubeResonator::structure`] makes the column alone a
//! [`fs_nlmodal::component::ModalStructure`]: one mass-normalized mode per plane-wave mode, damped by
//! its wall and open-end rates, with a `mouth` port (pressure at the mouth in, volume flux into the
//! column out) and the structure's `heat` port.
//!
//! No claim, all as in the law copied:
//! - at a change of radius only the plane wave is carried, so the inertance of the step's evanescent
//!   field and the loss on the step's annular face are left out; sloped pieces are cylinders
//!   [`BORE_CUT_M`] long;
//! - Eq. 12 holds below `ka` [`LARGEST_KA`]; a mode above it is refused;
//! - the uniform law stops its fixed point at a relative change of 1e-12 or after 64 iterations, and
//!   after 64 it keeps the last iterate (the drawn-bore walk refuses instead);
//! - the column alone books its open end's radiated power as heat through the `heat` port: it has no
//!   radiation port, and the field it radiates is not computed here.

use core::num::NonZeroUsize;

use fs_material::gas::GasState;
use fs_math::det;
use fs_nlmodal::component::{
    ModalAcousticTimeBudget, ModalAttachment, ModalStructure, ModalStructureRefusal,
};
use fs_ports::{PortKind, StableId};

const TAU: f64 = core::f64::consts::TAU;

/// Eq. 7's characteristic lengths of viscous and heat-propagation effects [m], as printed.
pub const VISCOUS_LENGTH_M: f64 = 4e-8;
/// Eq. 7's characteristic length of heat-propagation effects [m], as printed.
pub const THERMAL_LENGTH_M: f64 = 5.6e-8;
/// Eq. 12 holds below this, and the plane-wave column well below its cut-off (Eq. 5: ka = 1.84).
pub const LARGEST_KA: f64 = 1.5;
/// A sloped piece of a drawn bore is cut into cylinders no longer than this [m], each at its mid
/// radius.
pub const BORE_CUT_M: f64 = 1e-5;

/// Why a column or its modes were refused.
#[derive(Debug, Clone, PartialEq)]
pub enum ColumnRefusal {
    /// The geometry is not a column: radius not positive, column not longer than its diameter, no
    /// mode asked for, or a value that is not finite.
    Geometry,
    /// The drawn bore does not run from the floor face (height 0) to the mouth (the column's length,
    /// at its radius), or has fewer than two points.
    BoreEnds,
    /// The drawn bore has a height that falls or a radius that is not positive.
    BoreShape,
    /// The gas has a sound speed, density or heat-capacity ratio that cannot carry the law.
    Gas,
    /// Mode `mode` (from 1) reaches `ka` at or above [`LARGEST_KA`].
    OpenEndRange {
        /// Mode number, from 1.
        mode: usize,
        /// The `ka` reached.
        ka: f64,
    },
    /// Mode `mode` (from 1) of a drawn bore found no fixed point in 64 iterations.
    NoFixedPoint {
        /// Mode number, from 1.
        mode: usize,
    },
}

impl core::fmt::Display for ColumnRefusal {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Geometry => write!(
                f,
                "tube needs a positive radius, a column longer than its diameter, one mode or more, and finite values"
            ),
            Self::BoreEnds => write!(
                f,
                "tube_bore_m runs from the floor face (height 0) to the mouth (tube_air_length_m, at tube_radius_m)"
            ),
            Self::BoreShape => write!(
                f,
                "tube_bore_m needs heights that never fall and positive radii"
            ),
            Self::Gas => write!(
                f,
                "the gas needs a finite positive sound speed and density and a heat-capacity ratio of 1 or more"
            ),
            // The same words as the law copied, so that a comparison can read them.
            Self::OpenEndRange { mode, ka } => write!(
                f,
                "tube mode {mode} reaches ka {ka:.3}: the open-end law holds below {LARGEST_KA}; lower tube_modes"
            ),
            Self::NoFixedPoint { mode } => {
                write!(f, "tube mode {mode} has no fixed point over the drawn bore")
            }
        }
    }
}

impl std::error::Error for ColumnRefusal {}

/// A closed-open column: closed at the floor face, open at the mouth.
#[derive(Debug, Clone, PartialEq)]
pub struct ClosedOpenColumn {
    radius_m: f64,
    air_length_m: f64,
    modes: usize,
    bore: Option<Vec<[f64; 2]>>,
}

impl ClosedOpenColumn {
    /// A column of mouth radius `radius_m` and air length `air_length_m` from the floor face to the
    /// mouth, carrying its first `modes` odd plane-wave modes. `bore`, if given, is the bore as drawn:
    /// points (height above the floor face, radius) [m] joined by straight lines, from the floor face
    /// (height 0) to the mouth (`air_length_m`, at `radius_m`); two points at one height are a step.
    /// Without it the bore is one cylinder of `radius_m`.
    ///
    /// # Errors
    /// [`ColumnRefusal::Geometry`], [`ColumnRefusal::BoreEnds`] or [`ColumnRefusal::BoreShape`].
    // The ends and the merge are exact comparisons of finite values, as in the card parser.
    #[allow(clippy::float_cmp)]
    pub fn new(
        radius_m: f64,
        air_length_m: f64,
        modes: usize,
        bore: Option<Vec<[f64; 2]>>,
    ) -> Result<Self, ColumnRefusal> {
        if !(radius_m.is_finite()
            && air_length_m.is_finite()
            && radius_m > 0.0
            && air_length_m > 2.0 * radius_m
            && modes >= 1)
        {
            return Err(ColumnRefusal::Geometry);
        }
        if let Some(points) = &bore {
            if points.len() < 2 || points.iter().flatten().any(|v| !v.is_finite()) {
                return Err(ColumnRefusal::BoreEnds);
            }
            let (first, last) = (points[0], points[points.len() - 1]);
            if !(first[0] == 0.0 && last[0] == air_length_m && last[1] == radius_m) {
                return Err(ColumnRefusal::BoreEnds);
            }
            // Every value is finite here, so the plain comparisons decide.
            if points.windows(2).any(|p| p[1][0] < p[0][0]) || points.iter().any(|p| p[1] <= 0.0) {
                return Err(ColumnRefusal::BoreShape);
            }
        }
        Ok(Self {
            radius_m,
            air_length_m,
            modes,
            bore,
        })
    }

    /// Mouth radius [m].
    #[must_use]
    pub fn radius_m(&self) -> f64 {
        self.radius_m
    }

    /// Air length from the floor face to the mouth [m].
    #[must_use]
    pub fn air_length_m(&self) -> f64 {
        self.air_length_m
    }

    /// Number of odd plane-wave modes carried.
    #[must_use]
    pub fn modes(&self) -> usize {
        self.modes
    }

    /// The drawn bore, if any.
    #[must_use]
    pub fn bore(&self) -> Option<&[[f64; 2]]> {
        self.bore.as_deref()
    }
}

/// One plane-wave mode of the isolated column.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct IsolatedMode {
    /// Angular frequency [rad/s].
    pub omega_rad_s: f64,
    /// Amplitude decay rate through the wall [1/s], Eq. 7: alpha c.
    pub wall_per_s: f64,
    /// Amplitude decay rate through the open end [1/s]: c Re(Z)/L, with Z normalized by Z0.
    pub radiation_per_s: f64,
    /// `ka` at the mode's frequency, with the mouth's radius.
    pub ka: f64,
    /// Column length with the open end's correction [m]. Over a drawn bore, the length of a column
    /// of the mouth's radius with this mode's mass per unit volume velocity at the open end's plane.
    pub acoustic_length_m: f64,
}

/// Eq. 12, normalized radiation impedance of an unflanged open end, real and imaginary parts, for
/// ka < 1.5.
#[must_use]
pub fn open_end(ka: f64) -> (f64, f64) {
    let (x2, x3, x4) = (ka * ka, ka * ka * ka, ka * ka * ka * ka);
    let log = det::ln(ka);
    (
        x2 / 4.0 + 0.0127 * x4 + 0.082 * x4 * log - 0.023 * x3 * x3,
        0.6133 * ka - 0.036 * x3 + 0.034 * x3 * log - 0.0187 * x4 * ka,
    )
}

fn check_gas(gas: &GasState) -> Result<(), ColumnRefusal> {
    if gas.sound_speed.is_finite()
        && gas.sound_speed > 0.0
        && gas.density.is_finite()
        && gas.density > 0.0
        && gas.gamma.is_finite()
        && gas.gamma >= 1.0
    {
        Ok(())
    } else {
        Err(ColumnRefusal::Gas)
    }
}

/// The column's first `modes` odd modes: each frequency is a fixed point, since the wall slows the
/// wave and the open end's correction both depend on it. With a drawn bore, the walk up the bore.
///
/// # Errors
/// [`ColumnRefusal::Gas`], [`ColumnRefusal::OpenEndRange`], or (drawn bore only)
/// [`ColumnRefusal::NoFixedPoint`].
pub fn isolated_modes(
    column: &ClosedOpenColumn,
    gas: &GasState,
) -> Result<Vec<IsolatedMode>, ColumnRefusal> {
    check_gas(gas)?;
    if column.bore.is_some() {
        return drawn_modes(column, gas);
    }
    let c0 = gas.sound_speed;
    let a = column.radius_m;
    let eps = VISCOUS_LENGTH_M.sqrt() + (gas.gamma - 1.0) * THERMAL_LENGTH_M.sqrt();
    let wall = eps / (a * core::f64::consts::SQRT_2);
    (0..column.modes)
        .map(|m| {
            let quarter = (2 * m + 1) as f64 * 0.25 * TAU;
            let mut w = quarter * c0 / column.air_length_m;
            let mut out = None;
            for _ in 0..64 {
                let ka = w * a / c0;
                if ka >= LARGEST_KA {
                    return Err(ColumnRefusal::OpenEndRange { mode: m + 1, ka });
                }
                let (re, im) = open_end(ka);
                let length = column.air_length_m + a * im / ka;
                let c = c0 * (1.0 - wall * (c0 / w).sqrt());
                let next = quarter * c / length;
                out = Some(IsolatedMode {
                    omega_rad_s: next,
                    wall_per_s: wall * (next / c0).sqrt() * c,
                    radiation_per_s: c * re / length,
                    ka,
                    acoustic_length_m: length,
                });
                if (next - w).abs() <= 1e-12 * w {
                    break;
                }
                w = next;
            }
            // At least one pass ran, so `out` is set; after 64 passes the last iterate is kept, as in
            // the law copied.
            out.ok_or(ColumnRefusal::NoFixedPoint { mode: m + 1 })
        })
        .collect()
}

/// The drawn bore as cylinders from the floor face up, [length, radius] [m]: the points joined by
/// straight lines, each sloped piece cut into cylinders [`BORE_CUT_M`] long or less at their mid
/// radii, neighbours of one radius made one. Without a drawn bore, the one cylinder of the column's
/// radius.
#[must_use]
// The flat-piece test and the merge of equal radii are exact, as in the law copied.
#[allow(clippy::float_cmp)]
pub fn cylinders(column: &ClosedOpenColumn) -> Vec<[f64; 2]> {
    let Some(points) = &column.bore else {
        return vec![[column.air_length_m, column.radius_m]];
    };
    let mut out: Vec<[f64; 2]> = Vec::new();
    for pair in points.windows(2) {
        let ([za, ra], [zb, rb]) = (pair[0], pair[1]);
        if zb == za {
            continue;
        }
        let pieces = if ra == rb {
            1
        } else {
            ((zb - za) / BORE_CUT_M).ceil() as usize
        };
        for j in 0..pieces {
            let piece = if ra == rb {
                [zb - za, ra]
            } else {
                let n = pieces as f64;
                [(zb - za) / n, ra + (rb - ra) * (j as f64 + 0.5) / n]
            };
            match out.last_mut() {
                Some(last) if last[1] == piece[1] => last[0] += piece[0],
                _ => out.push(piece),
            }
        }
    }
    out
}

/// One cylinder of the walk up a drawn bore at one frequency.
struct Run {
    /// [m]; the last carries the open end's correction.
    length: f64,
    radius: f64,
    /// Eq. 7's slowed speed at this radius [m/s].
    speed: f64,
    /// Eq. 7's wall factor at this radius.
    wall: f64,
    /// Pressure amplitude, the floor's being 1: p = A cos(phase), Z0 U = A sin(phase).
    amplitude: f64,
    /// Phase at the cylinder's lower end [rad].
    phase: f64,
}

/// [`isolated_modes`] over the drawn bore. Each mode is the law's fixed point with the column's
/// length replaced by the walk's length in the mouth's units. The mass, the wall loss and the
/// open-end loss come from the mode's energy over the walk, with the open end's correction carried
/// at the mouth's radius to the plane where the pressure is zero; the mass is per unit volume
/// velocity there, as the uniform law's is.
fn drawn_modes(
    column: &ClosedOpenColumn,
    gas: &GasState,
) -> Result<Vec<IsolatedMode>, ColumnRefusal> {
    let c0 = gas.sound_speed;
    let a = column.radius_m;
    let viscous = VISCOUS_LENGTH_M.sqrt();
    let thermal = (gas.gamma - 1.0) * THERMAL_LENGTH_M.sqrt();
    let eps = viscous + thermal;
    // Eq. 7's viscous share less its thermal share.
    let split = (viscous - thermal) / eps;
    let bore = cylinders(column);
    let walls: Vec<f64> = bore
        .iter()
        .map(|c| eps / (c[1] * core::f64::consts::SQRT_2))
        .collect();
    let pi = core::f64::consts::PI;
    (0..column.modes)
        .map(|m| {
            let quarter = (2 * m + 1) as f64 * 0.25 * TAU;
            let mut w = quarter * c0 / column.air_length_m;
            let mut runs: Vec<Run> = Vec::with_capacity(bore.len());
            for _ in 0..64 {
                let ka = w * a / c0;
                if ka >= LARGEST_KA {
                    return Err(ColumnRefusal::OpenEndRange { mode: m + 1, ka });
                }
                let (re, im) = open_end(ka);
                let root = (c0 / w).sqrt();
                let mouth = c0 * (1.0 - walls[walls.len() - 1] * root);
                // The walk, its length in the mouth's units.
                runs.clear();
                let (mut length, mut amplitude) = (0.0, 1.0);
                for (i, (c, wall)) in bore.iter().zip(&walls).enumerate() {
                    let speed = c0 * (1.0 - wall * root);
                    if i > 0 {
                        // A change of radius: p and U continuous.
                        let phase = w * length / mouth;
                        let turns = (phase / pi).floor();
                        let (s, co) = det::sin_cos(phase - turns * pi);
                        let below = bore[i - 1][1];
                        let ratio = (below * below) / (c[1] * c[1]);
                        amplitude *= det::hypot(co, ratio * s);
                        length = (turns * pi + det::atan2(ratio * s, co)) * mouth / w;
                    }
                    runs.push(Run {
                        length: c[0],
                        radius: c[1],
                        speed,
                        wall: *wall,
                        amplitude,
                        phase: w * length / mouth,
                    });
                    length += c[0] * (mouth / speed);
                }
                let correction = a * im / ka;
                let column_length = length + correction;
                let next = quarter * mouth / column_length;
                if (next - w).abs() <= 1e-12 * w {
                    let last = runs.len() - 1;
                    runs[last].length += correction;
                    let (area, peak) = (a * a, runs[last].amplitude);
                    let (mut kinetic, mut potential) = (0.0, 0.0);
                    let mut parts = Vec::with_capacity(runs.len());
                    for (i, r) in runs.iter().enumerate() {
                        let k = w / r.speed;
                        // The last run ends where the pressure is zero.
                        let top = if i == last {
                            0.0
                        } else {
                            det::sin(2.0 * (r.phase + k * r.length))
                        };
                        let swing = (top - det::sin(2.0 * r.phase)) / (4.0 * k);
                        let weight = (r.radius * r.radius / area)
                            * (r.amplitude / peak)
                            * (r.amplitude / peak);
                        let (ke, pe) = (
                            weight * (r.length / 2.0 - swing),
                            weight * (r.length / 2.0 + swing),
                        );
                        kinetic += ke;
                        potential += pe;
                        parts.push((ke, pe));
                    }
                    let total = kinetic + potential;
                    let mass_length = 2.0 * kinetic;
                    let rise = (next / c0).sqrt();
                    let wall_per_s = runs.iter().zip(&parts).fold(0.0, |sum, (r, (ke, pe))| {
                        sum + r.wall
                            * rise
                            * r.speed
                            * ((ke + pe) / total + split * (ke - pe) / total)
                    });
                    return Ok(IsolatedMode {
                        omega_rad_s: next,
                        wall_per_s,
                        radiation_per_s: mouth * re / mass_length,
                        ka,
                        acoustic_length_m: mass_length,
                    });
                }
                w = next;
            }
            Err(ColumnRefusal::NoFixedPoint { mode: m + 1 })
        })
        .collect()
}

/// A mode's acoustic mass per unit mouth volume velocity [kg/m^4]: rho L / (2 S), with the mode's own
/// acoustic length.
#[must_use]
pub fn acoustic_mass(gas: &GasState, column: &ClosedOpenColumn, mode: &IsolatedMode) -> f64 {
    let area = 0.5 * TAU * column.radius_m * column.radius_m;
    gas.density * mode.acoustic_length_m / (2.0 * area)
}

/// Why a [`TubeResonator`] or its structure was refused.
#[derive(Debug, Clone, PartialEq)]
pub enum TubeResonatorRefusal {
    /// The column or its modes were refused.
    Column(ColumnRefusal),
    /// A mode's total decay rate is not positive and below its angular frequency, or its acoustic
    /// mass is not positive and finite.
    Mode {
        /// Mode number, from 1.
        mode: usize,
    },
    /// The modal structure was refused.
    Structure(ModalStructureRefusal),
}

impl core::fmt::Display for TubeResonatorRefusal {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Column(refusal) => write!(f, "{refusal}"),
            Self::Mode { mode } => write!(
                f,
                "tube mode {mode} needs a positive underdamped decay and a positive finite acoustic mass"
            ),
            Self::Structure(refusal) => write!(f, "{refusal}"),
        }
    }
}

impl std::error::Error for TubeResonatorRefusal {}

impl From<ColumnRefusal> for TubeResonatorRefusal {
    fn from(refusal: ColumnRefusal) -> Self {
        Self::Column(refusal)
    }
}

/// A closed-open column's modes with their acoustic masses.
#[derive(Debug, Clone, PartialEq)]
pub struct TubeResonator {
    modes: Vec<IsolatedMode>,
    acoustic_masses_kg_m4: Vec<f64>,
}

impl TubeResonator {
    /// The column's modes by [`isolated_modes`] and their masses by [`acoustic_mass`].
    ///
    /// # Errors
    /// The law's refusals, or [`TubeResonatorRefusal::Mode`] for a mode whose total decay rate is not
    /// positive and below its angular frequency or whose mass is not positive and finite.
    pub fn new(column: &ClosedOpenColumn, gas: &GasState) -> Result<Self, TubeResonatorRefusal> {
        let modes = isolated_modes(column, gas)?;
        let mut masses = Vec::with_capacity(modes.len());
        for (index, mode) in modes.iter().enumerate() {
            let decay = mode.wall_per_s + mode.radiation_per_s;
            let mass = acoustic_mass(gas, column, mode);
            if !(decay > 0.0 && decay < mode.omega_rad_s && mass.is_finite() && mass > 0.0) {
                return Err(TubeResonatorRefusal::Mode { mode: index + 1 });
            }
            masses.push(mass);
        }
        Ok(Self {
            modes,
            acoustic_masses_kg_m4: masses,
        })
    }

    /// The column's modes.
    #[must_use]
    pub fn modes(&self) -> &[IsolatedMode] {
        &self.modes
    }

    /// Each mode's acoustic mass per unit mouth volume velocity [kg/m^4].
    #[must_use]
    pub fn acoustic_masses_kg_m4(&self) -> &[f64] {
        &self.acoustic_masses_kg_m4
    }

    /// Each mode as `(angular frequency [rad/s], viscous damping ratio)`: the ratio is the mode's
    /// wall and open-end amplitude decay rates over its angular frequency.
    #[must_use]
    pub fn modal_rows(&self) -> Vec<(f64, f64)> {
        self.modes
            .iter()
            .map(|m| {
                (
                    m.omega_rad_s,
                    (m.wall_per_s + m.radiation_per_s) / m.omega_rad_s,
                )
            })
            .collect()
    }

    /// The `mouth` row: mode `k`'s coordinate is its volume displacement into the column times the
    /// root of its acoustic mass, so the pressure at the mouth acts on it through `1 / sqrt(M_k)`
    /// [m^2/sqrt(kg)] and its velocity gives the volume flux into the column through the same
    /// number.
    #[must_use]
    pub fn mouth_weights(&self) -> Vec<f64> {
        self.acoustic_masses_kg_m4
            .iter()
            .map(|m| 1.0 / m.sqrt())
            .collect()
    }

    /// The column alone as a modal structure at the owner's nominal clock `sample_rate_hz`: ports
    /// `mouth` ([`PortKind::FluidPressureFlux`], the pressure at the mouth pushing into the column in,
    /// the mean volume flux into the column out) and `heat`.
    ///
    /// # Errors
    /// The structure's refusals (a mode above the clock's Nyquist guard among them).
    pub fn structure(
        &self,
        sample_rate_hz: u32,
        budget: ModalAcousticTimeBudget,
    ) -> Result<ModalStructure, TubeResonatorRefusal> {
        let mouth = ModalAttachment {
            id: StableId::new("mouth").map_err(|e| {
                TubeResonatorRefusal::Structure(ModalStructureRefusal::Identifier(e))
            })?,
            kind: PortKind::FluidPressureFlux,
            width: NonZeroUsize::MIN,
            weights: self.mouth_weights(),
        };
        ModalStructure::from_modes(sample_rate_hz, &self.modal_rows(), budget, vec![mouth])
            .map_err(TubeResonatorRefusal::Structure)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The record bar_strike's cards name (`dry-air-20c.gas`), read the same way.
    fn air() -> GasState {
        fs_material::gas_species::input::read_gas_state(
            &mut "293.15 101325 0".split_ascii_whitespace(),
            None,
        )
        .expect("dry air at 20 C")
    }

    /// La Favre's A3 resonator tube.
    fn a3() -> ClosedOpenColumn {
        let inch = 0.0254;
        ClosedOpenColumn::new((1.25 - 0.0625) * inch, 14.375 * inch, 5, None).expect("A3 column")
    }

    /// The printed laws at La Favre's A3, against the sums by hand: Eq. 7 wall, Eq. 12 open end, odd
    /// modes (the law copied carries the same test).
    #[test]
    fn isolated_column_follows_the_printed_laws() {
        let g = air();
        let t = a3();
        let modes = isolated_modes(&t, &g).expect("modes");
        let f: Vec<f64> = modes.iter().map(|m| m.omega_rad_s / TAU).collect();
        assert!((f[0] - 223.5).abs() < 1.0, "{f:?}");
        assert!((modes[0].wall_per_s + modes[0].radiation_per_s - 8.2).abs() < 0.3);
        // Eq. 12 at low ka reduces to Eq. 13.
        let (re, im) = open_end(0.01);
        assert!(
            (re / (1e-4 / 4.0) - 1.0).abs() < 1e-3 && (im / (0.6133 * 0.01) - 1.0).abs() < 1e-3
        );
        // Odd modes, each above its odd multiple of the first by the two printed effects and nothing
        // else: the open end's correction shrinks as ka rises, and the wall slows the low modes most.
        let wall = (4e-8_f64.sqrt() + (g.gamma - 1.0) * 5.6e-8_f64.sqrt())
            / (t.radius_m() * core::f64::consts::SQRT_2);
        let slowed = |w: f64| 1.0 - wall * (g.sound_speed / w).sqrt();
        for (m, mode) in modes.iter().enumerate() {
            let odd = (2 * m + 1) as f64;
            let by_laws = odd
                * (modes[0].acoustic_length_m / mode.acoustic_length_m)
                * (slowed(mode.omega_rad_s) / slowed(modes[0].omega_rad_s));
            let ratio = f[m] / f[0];
            assert!(
                (ratio / by_laws - 1.0).abs() < 1e-6,
                "mode {m} ratio {ratio} against {by_laws}"
            );
        }
        // Mode 3's ratio by hand from Eq. 7 and Eq. 12, worked apart from this code.
        assert!((f[3] / f[0] - 7.0369).abs() < 0.0005, "{f:?}");
        // The wall rises as the root of frequency, the open end much faster.
        assert!((modes[2].wall_per_s / modes[0].wall_per_s - 5.0_f64.sqrt()).abs() < 0.05);
        assert!(modes[2].radiation_per_s > 15.0 * modes[0].radiation_per_s);
    }

    /// The walk over a two-cylinder bore (125 mm at radius 20 mm under 250 mm at 30 mm) against the
    /// bore's own condition, solved apart from the walk: with each cylinder's Eq. 7 speed and Eq. 12's
    /// correction at the mouth, a mode is a root of
    /// `S2 cos(k1 l1) cos(k2 (l2 + delta)) - S1 sin(k1 l1) sin(k2 (l2 + delta))`. The mass length is
    /// checked against a midpoint sum of the kinetic energy over that mode's shape (the law copied
    /// carries the same test).
    #[test]
    fn a_two_cylinder_bore_meets_its_own_condition() {
        let g = air();
        let (l1, r1, l2, r2) = (0.125, 0.020, 0.25, 0.030);
        let t = ClosedOpenColumn::new(
            r2,
            l1 + l2,
            3,
            Some(vec![[0.0, r1], [l1, r1], [l1, r2], [l1 + l2, r2]]),
        )
        .expect("two cylinders");
        let c0 = g.sound_speed;
        let eps = 4e-8_f64.sqrt() + (g.gamma - 1.0) * 5.6e-8_f64.sqrt();
        let speed =
            |w: f64, r: f64| c0 * (1.0 - eps / (r * core::f64::consts::SQRT_2) * (c0 / w).sqrt());
        let shape = |w: f64| {
            let ka = w * r2 / c0;
            let (_, im) = open_end(ka);
            (w / speed(w, r1), w / speed(w, r2), l2 + r2 * im / ka)
        };
        let condition = |w: f64| {
            let (k1, k2, top) = shape(w);
            r2 * r2 * (k1 * l1).cos() * (k2 * top).cos()
                - r1 * r1 * (k1 * l1).sin() * (k2 * top).sin()
        };
        let modes = isolated_modes(&t, &g).expect("the walk");
        assert_eq!(modes.len(), 3);
        for (m, mode) in modes.iter().enumerate() {
            let w = mode.omega_rad_s;
            let (mut lo, mut hi) = (0.97 * w, 1.03 * w);
            assert!(
                condition(lo) * condition(hi) < 0.0,
                "mode {m}: no root bracketed"
            );
            for _ in 0..200 {
                let mid = 0.5 * (lo + hi);
                if condition(lo) * condition(mid) <= 0.0 {
                    hi = mid;
                } else {
                    lo = mid;
                }
            }
            let root = 0.5 * (lo + hi);
            assert!(
                (w / root - 1.0).abs() < 1e-10,
                "mode {m}: walk {w} rad/s against the condition's {root}"
            );
            let (k1, k2, top) = shape(root);
            let ratio = (r1 * r1) / (r2 * r2);
            let phase = (ratio * (k1 * l1).sin()).atan2((k1 * l1).cos());
            let upper = (k1 * l1).cos().hypot(ratio * (k1 * l1).sin());
            let n: i32 = 200_000;
            let sum = |length: f64, k: f64, start: f64| {
                (0..n)
                    .map(|j| {
                        (k * (f64::from(j) + 0.5) * length / f64::from(n) + start)
                            .sin()
                            .powi(2)
                    })
                    .sum::<f64>()
                    * length
                    / f64::from(n)
            };
            let mass_length =
                2.0 * (ratio / (upper * upper) * sum(l1, k1, 0.0) + sum(top, k2, phase));
            assert!(
                (mode.acoustic_length_m / mass_length - 1.0).abs() < 1e-9,
                "mode {m}: mass length {} against the sum's {mass_length}",
                mode.acoustic_length_m
            );
        }
    }

    #[test]
    fn one_cylinder_drawn_is_the_uniform_law_bit_for_bit() {
        let gas = air();
        let inch = 0.0254;
        for (radius, length, modes) in [
            (0.0265, 0.239, 4),
            (0.030, 0.304, 5),
            ((1.25 - 0.0625) * inch, 14.375 * inch, 5),
        ] {
            let uniform = ClosedOpenColumn::new(radius, length, modes, None).expect("column");
            let drawn = ClosedOpenColumn::new(
                radius,
                length,
                modes,
                Some(vec![[0.0, radius], [length, radius]]),
            )
            .expect("drawn column");
            let (a, b) = (
                isolated_modes(&uniform, &gas).expect("uniform modes"),
                isolated_modes(&drawn, &gas).expect("drawn modes"),
            );
            assert_eq!((a.len(), b.len()), (modes, modes));
            assert_eq!(cylinders(&drawn), [[length, radius]]);
            for (x, y) in a.iter().zip(&b) {
                for (p, q) in [
                    (x.omega_rad_s, y.omega_rad_s),
                    (x.wall_per_s, y.wall_per_s),
                    (x.radiation_per_s, y.radiation_per_s),
                    (x.ka, y.ka),
                    (x.acoustic_length_m, y.acoustic_length_m),
                ] {
                    assert_eq!(p.to_bits(), q.to_bits(), "radius {radius} length {length}");
                }
            }
        }
    }

    #[test]
    fn modes_ascend_and_refuse_past_the_open_end_law() {
        let gas = air();
        let column = ClosedOpenColumn::new(0.0265, 0.239, 4, None).expect("column");
        let modes = isolated_modes(&column, &gas).expect("modes");
        for pair in modes.windows(2) {
            assert!(pair[1].omega_rad_s > pair[0].omega_rad_s);
        }
        // ka reaches 1.5 near 3.1 kHz at this radius; the twelfth odd mode of a 0.239 m column is
        // far above that, so the law must refuse on the way.
        let many = ClosedOpenColumn::new(0.0265, 0.239, 12, None).expect("column");
        let Err(ColumnRefusal::OpenEndRange { mode, ka }) = isolated_modes(&many, &gas) else {
            panic!("expected the open-end range refusal");
        };
        assert!(ka >= LARGEST_KA, "mode {mode} ka {ka}");
        // The refused mode is numbered from 1. Asked for exactly that many modes, the law refuses the
        // same mode at the same ka; asked for the modes below it, it answers with all of them.
        let exactly = ClosedOpenColumn::new(0.0265, 0.239, mode, None).expect("column");
        assert!(matches!(
            isolated_modes(&exactly, &gas),
            Err(ColumnRefusal::OpenEndRange { mode: m, ka: k }) if m == mode && k.to_bits() == ka.to_bits()
        ));
        let below = ClosedOpenColumn::new(0.0265, 0.239, mode - 1, None).expect("column");
        let within = isolated_modes(&below, &gas).expect("the modes below the refused one");
        assert_eq!(within.len(), mode - 1);
        assert!(within.iter().all(|m| m.ka < LARGEST_KA));
    }

    #[test]
    fn geometry_and_gas_refusals() {
        assert_eq!(
            ClosedOpenColumn::new(0.0, 0.2, 1, None),
            Err(ColumnRefusal::Geometry)
        );
        assert_eq!(
            ClosedOpenColumn::new(0.1, 0.2, 1, None),
            Err(ColumnRefusal::Geometry)
        );
        assert_eq!(
            ClosedOpenColumn::new(0.02, 0.2, 0, None),
            Err(ColumnRefusal::Geometry)
        );
        assert_eq!(
            ClosedOpenColumn::new(0.02, f64::INFINITY, 1, None),
            Err(ColumnRefusal::Geometry)
        );
        assert_eq!(
            ClosedOpenColumn::new(0.02, 0.2, 1, Some(vec![[0.0, 0.02], [0.19, 0.02]])),
            Err(ColumnRefusal::BoreEnds)
        );
        assert_eq!(
            ClosedOpenColumn::new(
                0.02,
                0.2,
                1,
                Some(vec![[0.0, 0.02], [0.1, 0.02], [0.05, 0.02], [0.2, 0.02]])
            ),
            Err(ColumnRefusal::BoreShape)
        );
        let column = ClosedOpenColumn::new(0.02, 0.2, 1, None).expect("column");
        let mut gas = air();
        gas.sound_speed = 0.0;
        assert_eq!(isolated_modes(&column, &gas), Err(ColumnRefusal::Gas));
    }

    #[test]
    fn the_structure_carries_the_law_unchanged() {
        let gas = air();
        let column = ClosedOpenColumn::new(0.0265, 0.239, 4, None).expect("column");
        let tube = TubeResonator::new(&column, &gas).expect("tube");
        let structure = tube
            .structure(44_100, ModalAcousticTimeBudget::audible_reference())
            .expect("structure");
        assert_eq!(structure.modes(), 4);
        for ((mode, row), (law, mass)) in structure
            .model()
            .modes()
            .iter()
            .zip(tube.modal_rows())
            .zip(tube.modes().iter().zip(tube.acoustic_masses_kg_m4()))
        {
            assert_eq!(
                mode.angular_frequency_rad_s.to_bits(),
                law.omega_rad_s.to_bits()
            );
            assert_eq!(mode.angular_frequency_rad_s.to_bits(), row.0.to_bits());
            assert_eq!(mode.damping_ratio.to_bits(), row.1.to_bits());
            // rho L / (2 S) by hand, with S = pi a^2.
            let by_hand = gas.density * law.acoustic_length_m
                / (2.0 * core::f64::consts::PI * 0.0265 * 0.0265);
            assert!(
                (mass / by_hand - 1.0).abs() < 8.0 * f64::EPSILON,
                "{mass} {by_hand}"
            );
        }
        let ids: Vec<&str> = fs_ports::component::Component::ports(&structure)
            .iter()
            .map(|p| p.id().as_str())
            .collect();
        assert_eq!(ids, ["mouth", "heat"]);
        let mouth = &fs_ports::component::Component::ports(&structure)[0];
        assert_eq!(
            (mouth.width(), mouth.role()),
            (1, fs_ports::component::PortRole::EffortIn)
        );
    }

    /// One step from rest under a held mouth pressure: the mean volume flux the structure returns
    /// is the sum over the modes of `(1 / sqrt(M_k)) q_k / dt`, with `q_k` the damped oscillator's
    /// closed-form response to the held force `p / sqrt(M_k)`. The closed form loses about
    /// `eps / (omega dt)^2` to cancellation, 2e-13 for the first mode at this step; the bound is
    /// that with room for the owner's own rounding of the same cancellation.
    #[test]
    fn a_held_mouth_pressure_drives_each_mode_through_its_root_mass() {
        use fs_alloc::{ArenaConfig, ArenaPool};
        use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
        use fs_ports::component::{Component, Proposal};
        let gas = air();
        let column = ClosedOpenColumn::new(0.0265, 0.239, 4, None).expect("column");
        let tube = TubeResonator::new(&column, &gas).expect("tube");
        let structure = tube
            .structure(44_100, ModalAcousticTimeBudget::audible_reference())
            .expect("structure");
        let (pressure_pa, dt_s) = (0.2, 1.0 / 44_100.0);
        let start = structure.at_rest();
        let mut next = start.clone();
        let mut outputs = [0.0; 2];
        let gate = CancelGate::new_clock_free();
        let pool = ArenaPool::new(ArenaConfig::default());
        let proposal = pool.scope(|arena| {
            let key = StreamKey {
                seed: 1,
                kernel_id: 2,
                tile: 3,
                iteration: 4,
            };
            let cx = Cx::new(&gate, arena, key, Budget::INFINITE, ExecMode::Deterministic);
            structure.propose(
                &start,
                dt_s,
                &[pressure_pa, 293.15],
                &mut next,
                &mut outputs,
                &cx,
            )
        });
        assert!(matches!(proposal, Ok(Proposal::Step(_))), "{proposal:?}");
        let mut flux = 0.0;
        for (mode, mass) in tube.modes().iter().zip(tube.acoustic_masses_kg_m4()) {
            let (w, s) = (mode.omega_rad_s, mode.wall_per_s + mode.radiation_per_s);
            let wd = (w * w - s * s).sqrt();
            let force = pressure_pa / mass.sqrt();
            let q = force / (w * w)
                * (1.0 - (-s * dt_s).exp() * ((wd * dt_s).cos() + s / wd * (wd * dt_s).sin()));
            flux += q / mass.sqrt() / dt_s;
        }
        assert!(flux > 0.0);
        assert!(
            (outputs[0] / flux - 1.0).abs() < 1e-9,
            "structure {} against the closed form {flux}",
            outputs[0]
        );
    }
}
