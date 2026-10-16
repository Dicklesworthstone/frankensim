//! Free-surface LBM (bead tfz.19): Körner-lineage mass-tracking VOF.
//! Interface cells carry a partial mass m (fill fraction κ = m/ρ);
//! mass moves along streaming links with PAIRWISE-ANTISYMMETRIC
//! exchange terms (so the ledger balances exactly by construction),
//! missing populations from gas neighbors are reconstructed with the
//! atmospheric-pressure anti-population closure (optionally shifted
//! by a Laplace term from fill-field curvature), and cell conversions
//! (fluid ↔ interface ↔ gas) redistribute their excess/deficit mass
//! conservatively — with a carry accumulator so NOTHING is ever
//! silently dropped. The battery gates the ledger at 1e-10 relative
//! EVERY step; contact-line physics is deliberately MODEL-BRACKETED
//! (neutral vs wetting wall ghosts), not pretended-certain.

use crate::core2::{Cell, Grid};
use crate::{CS2, E, OPP, Q, equilibrium};
use fs_matdb::{
    InterfaceSystemCard, MatDbError, MaterialAnswer, PropertyValue, QueryPoint, SelectionPolicy,
};
use fs_qty::{Dims, QtyAny};

/// MatDB property name consumed for the lower wetting-hysteresis endpoint.
pub const RECEDING_CONTACT_ANGLE_PROPERTY: &str = "receding-contact-angle";
/// MatDB property name consumed for the upper wetting-hysteresis endpoint.
pub const ADVANCING_CONTACT_ANGLE_PROPERTY: &str = "advancing-contact-angle";
/// MatDB query/curve axis carrying absolute contact-line speed.
pub const DYNAMIC_WETTING_SPEED_AXIS: &str = "contact_line_speed";
/// Required SI dimensions of the dynamic-wetting speed axis (`m s^-1`).
pub const DYNAMIC_WETTING_SPEED_DIMS: Dims = Dims([1, 0, -1, 0, 0, 0]);

/// Wall wetting model for the curvature fill-field ghost — the
/// contact-line bracket per the plan's honesty clause.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ContactModel {
    /// Neutral (≈90°): wall ghost copies the adjacent cell's fill.
    Neutral,
    /// Wetting: wall ghost is full (φ = 1) — the fluid is drawn along
    /// the wall.
    Wetting,
}

/// Directional state selected from a dynamic wetting hysteresis bracket.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ContactLineRegime2 {
    /// The signed contact-line speed is below the negative pinning deadband.
    Receding,
    /// The signed contact-line speed lies inside the closed pinning deadband.
    Pinned,
    /// The signed contact-line speed is above the positive pinning deadband.
    Advancing,
}

/// System- and endpoint-provenance-bound wetting selection for one interface.
///
/// Both material answers are retained even when only one endpoint is selected.
/// A pinned answer deliberately has no selected angle: its admissible output is
/// the full receding-to-advancing bracket, not an invented static value.
#[derive(Debug, Clone, PartialEq)]
#[must_use]
pub struct DynamicWettingAnswer2 {
    /// Canonical identity of the ordered, history-bearing interface card.
    pub interface_system_hash: [u8; 32],
    /// Directional regime selected by the signed speed and deadband.
    pub regime: ContactLineRegime2,
    /// Caller-supplied signed contact-line speed in SI metres per second.
    pub signed_contact_line_speed: f64,
    /// Symmetric pinning deadband in SI metres per second.
    pub pinning_speed: f64,
    /// Receipt-bearing receding-angle answer at the absolute speed.
    pub receding: MaterialAnswer,
    /// Receipt-bearing advancing-angle answer at the absolute speed.
    pub advancing: MaterialAnswer,
    /// Selected endpoint angle in radians, or `None` while pinned.
    pub selected_angle_radians: Option<f64>,
}

/// Fail-closed diagnostics for dynamic wetting-card admission and selection.
#[derive(Debug, Clone, PartialEq)]
pub enum DynamicWettingError2 {
    /// The underlying material query refused.
    MatDb(MatDbError),
    /// A caller-supplied speed does not carry SI velocity dimensions.
    SpeedDimsMismatch {
        /// Name of the rejected input.
        input: &'static str,
        /// Dimensions carried by the input.
        dims: Dims,
    },
    /// The signed contact-line speed is not finite.
    NonFiniteSpeed {
        /// Exact rejected floating-point representation.
        bits: u64,
    },
    /// The pinning speed is negative or non-finite.
    InvalidPinningSpeed {
        /// Exact rejected floating-point representation.
        bits: u64,
    },
    /// A contact-angle claim is not dimensionless (radians).
    NonDimensionlessAngle {
        /// Name of the rejected property.
        property: &'static str,
        /// Dimensions carried by the material answer.
        dims: Dims,
    },
    /// A speed-dependent angle curve labels its abscissa with non-velocity dimensions.
    ContactLineSpeedDimsMismatch {
        /// Name of the rejected angle property.
        property: &'static str,
        /// Dimensions declared by the selected curve's speed abscissa.
        dims: Dims,
    },
    /// A contact-angle answer is not finite and strictly inside `(0, pi)`.
    AngleOutOfRange {
        /// Name of the rejected property.
        property: &'static str,
        /// Rejected SI value in radians.
        value: f64,
    },
    /// The card's lower endpoint exceeds its upper endpoint.
    ReversedHysteresis {
        /// Receding angle in radians.
        receding: f64,
        /// Advancing angle in radians.
        advancing: f64,
    },
}

impl core::fmt::Display for DynamicWettingError2 {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::MatDb(error) => write!(f, "dynamic wetting material query refused: {error}"),
            Self::SpeedDimsMismatch { input, dims } => write!(
                f,
                "dynamic wetting input '{input}' carries {dims:?}, expected SI velocity dimensions {DYNAMIC_WETTING_SPEED_DIMS:?}"
            ),
            Self::NonFiniteSpeed { bits } => write!(
                f,
                "signed contact-line speed is non-finite (bits {bits:#018x})"
            ),
            Self::InvalidPinningSpeed { bits } => write!(
                f,
                "pinning speed must be finite and nonnegative (bits {bits:#018x})"
            ),
            Self::NonDimensionlessAngle { property, dims } => write!(
                f,
                "dynamic wetting property '{property}' must be dimensionless radians, got {dims:?}"
            ),
            Self::ContactLineSpeedDimsMismatch { property, dims } => write!(
                f,
                "dynamic wetting property '{property}' declares '{DYNAMIC_WETTING_SPEED_AXIS}' with {dims:?}, expected SI velocity dimensions {DYNAMIC_WETTING_SPEED_DIMS:?}"
            ),
            Self::AngleOutOfRange { property, value } => write!(
                f,
                "dynamic wetting property '{property}' must be finite and strictly inside (0, pi) radians, got {value}"
            ),
            Self::ReversedHysteresis {
                receding,
                advancing,
            } => write!(
                f,
                "dynamic wetting hysteresis is reversed: receding angle {receding} exceeds advancing angle {advancing} radians"
            ),
        }
    }
}

impl core::error::Error for DynamicWettingError2 {
    fn source(&self) -> Option<&(dyn core::error::Error + 'static)> {
        match self {
            Self::MatDb(error) => Some(error),
            _ => None,
        }
    }
}

impl From<MatDbError> for DynamicWettingError2 {
    fn from(error: MatDbError) -> Self {
        Self::MatDb(error)
    }
}

/// Query and select a system- and endpoint-provenance-bound contact-angle bracket.
///
/// Both speed inputs are runtime-dimensioned [`QtyAny`] values in coherent SI
/// and must carry velocity dimensions. The query point is cloned and its
/// `contact_line_speed` axis is replaced by the absolute supplied speed. Both
/// receding and advancing claims must answer under the same explicit selection
/// policy. Positive speed selects the advancing endpoint, negative speed
/// selects the receding endpoint, and the closed deadband
/// `[-pinning_speed, pinning_speed]` remains unselected.
///
/// This adapter performs admission and selection only. It does not impose the
/// returned angle in `FreeSurface::phi_at` or claim a contact-line law. The
/// caller-supplied pinning deadband is retained but carries no independent
/// calibration evidence or MatDB usage receipt.
///
/// # Errors
/// [`DynamicWettingError2`] for invalid speeds, material-query refusal,
/// dimensional or angular invalidity, or a reversed hysteresis bracket.
pub fn query_dynamic_wetting2(
    card: &InterfaceSystemCard,
    point: &QueryPoint,
    signed_contact_line_speed: QtyAny,
    pinning_speed: QtyAny,
    policy: SelectionPolicy,
) -> Result<DynamicWettingAnswer2, DynamicWettingError2> {
    if signed_contact_line_speed.dims != DYNAMIC_WETTING_SPEED_DIMS {
        return Err(DynamicWettingError2::SpeedDimsMismatch {
            input: "signed_contact_line_speed",
            dims: signed_contact_line_speed.dims,
        });
    }
    if pinning_speed.dims != DYNAMIC_WETTING_SPEED_DIMS {
        return Err(DynamicWettingError2::SpeedDimsMismatch {
            input: "pinning_speed",
            dims: pinning_speed.dims,
        });
    }
    let signed_contact_line_speed = signed_contact_line_speed.value;
    let pinning_speed = pinning_speed.value;
    if !signed_contact_line_speed.is_finite() {
        return Err(DynamicWettingError2::NonFiniteSpeed {
            bits: signed_contact_line_speed.to_bits(),
        });
    }
    if !pinning_speed.is_finite() || pinning_speed < 0.0 {
        return Err(DynamicWettingError2::InvalidPinningSpeed {
            bits: pinning_speed.to_bits(),
        });
    }

    let query_point = point
        .clone()
        .with(DYNAMIC_WETTING_SPEED_AXIS, signed_contact_line_speed.abs())?;
    let receding = card
        .claims()
        .query(RECEDING_CONTACT_ANGLE_PROPERTY, &query_point, policy)?;
    let advancing = card
        .claims()
        .query(ADVANCING_CONTACT_ANGLE_PROPERTY, &query_point, policy)?;

    let receding_angle = checked_contact_angle(card, RECEDING_CONTACT_ANGLE_PROPERTY, &receding)?;
    let advancing_angle =
        checked_contact_angle(card, ADVANCING_CONTACT_ANGLE_PROPERTY, &advancing)?;
    if receding_angle > advancing_angle {
        return Err(DynamicWettingError2::ReversedHysteresis {
            receding: receding_angle,
            advancing: advancing_angle,
        });
    }

    let (regime, selected_angle_radians) = if signed_contact_line_speed > pinning_speed {
        (ContactLineRegime2::Advancing, Some(advancing_angle))
    } else if signed_contact_line_speed < -pinning_speed {
        (ContactLineRegime2::Receding, Some(receding_angle))
    } else {
        (ContactLineRegime2::Pinned, None)
    };

    Ok(DynamicWettingAnswer2 {
        interface_system_hash: card.content_hash().0,
        regime,
        signed_contact_line_speed,
        pinning_speed,
        receding,
        advancing,
        selected_angle_radians,
    })
}

fn checked_contact_angle(
    card: &InterfaceSystemCard,
    property: &'static str,
    answer: &MaterialAnswer,
) -> Result<f64, DynamicWettingError2> {
    let selected_claim =
        card.claims()
            .claim(answer.receipt.selected)
            .ok_or(DynamicWettingError2::MatDb(MatDbError::ReceiptMismatch {
                field: "selected",
            }))?;
    if let PropertyValue::Curve {
        abscissa,
        abscissa_dims,
        ..
    } = &selected_claim.value
        && abscissa == DYNAMIC_WETTING_SPEED_AXIS
        && *abscissa_dims != DYNAMIC_WETTING_SPEED_DIMS
    {
        return Err(DynamicWettingError2::ContactLineSpeedDimsMismatch {
            property,
            dims: *abscissa_dims,
        });
    }
    let sample = &answer.evidence.value;
    if sample.dims != Dims::NONE {
        return Err(DynamicWettingError2::NonDimensionlessAngle {
            property,
            dims: sample.dims,
        });
    }
    let value = sample.value;
    if !value.is_finite() || value <= 0.0 || value >= core::f64::consts::PI {
        return Err(DynamicWettingError2::AngleOutOfRange { property, value });
    }
    Ok(value)
}

/// Free-surface simulation state.
pub struct FreeSurface {
    /// The lattice (flags: Fluid/Interface/Gas/Wall).
    pub grid: Grid,
    /// Per-cell tracked mass (meaningful for Interface; Fluid uses
    /// Σf, Gas is 0).
    pub mass: Vec<f64>,
    /// Surface tension coefficient (0 = off).
    pub sigma: f64,
    /// Contact-line model (curvature ghost at walls).
    pub contact: ContactModel,
    /// Conversion-conservation carry (mass awaiting redistribution).
    pub carry: f64,
    /// Cell-conversion statistics (cumulative).
    pub conversions: ConversionStats,
    scratch: Vec<[f64; Q]>,
    fill_smooth: Vec<f64>,
}

/// Cumulative conversion ledger.
#[derive(Debug, Clone, Copy, Default)]
pub struct ConversionStats {
    /// Interface → fluid events.
    pub to_fluid: u64,
    /// Interface → gas events.
    pub to_gas: u64,
    /// Gas → interface events (closure repair).
    pub gas_to_interface: u64,
    /// Fluid → interface events (closure repair).
    pub fluid_to_interface: u64,
}

const EPS: f64 = 1e-3;

impl FreeSurface {
    /// Build from a grid whose flags are already set (Fluid regions,
    /// Gas elsewhere, Wall boundaries). Interface cells are inserted
    /// automatically between fluid and gas; masses initialized
    /// (fluid = ρ, interface = ρ/2).
    ///
    /// # Panics
    /// If a fluid cell touches a gas cell after interface insertion
    /// (impossible by construction).
    #[must_use]
    pub fn new(mut grid: Grid, sigma: f64, contact: ContactModel) -> FreeSurface {
        let (nx, ny) = (grid.nx, grid.ny);
        // Insert interface cells: any fluid cell with a gas neighbor.
        let mut promote = Vec::new();
        for y in 0..ny {
            for x in 0..nx {
                let i = grid.idx(x, y);
                if grid.flags[i] != Cell::Fluid {
                    continue;
                }
                for q in 1..Q {
                    if let Some(nb) = neighbor(&grid, x, y, q)
                        && grid.flags[nb] == Cell::Gas
                    {
                        promote.push(i);
                        break;
                    }
                }
            }
        }
        for i in promote {
            grid.flags[i] = Cell::Interface;
        }
        let mut mass = vec![0.0f64; nx * ny];
        for (i, m) in mass.iter_mut().enumerate().take(nx * ny) {
            match grid.flags[i] {
                Cell::Fluid => *m = grid.f[i].iter().sum(),
                Cell::Interface => *m = 0.5 * grid.f[i].iter().sum::<f64>(),
                _ => {}
            }
        }
        let fs = FreeSurface {
            grid,
            mass,
            sigma,
            contact,
            carry: 0.0,
            conversions: ConversionStats::default(),
            scratch: Vec::new(),
            fill_smooth: vec![0.0; nx * ny],
        };
        fs.assert_closure();
        fs
    }

    fn assert_closure(&self) {
        for y in 0..self.grid.ny {
            for x in 0..self.grid.nx {
                let i = self.grid.idx(x, y);
                if self.grid.flags[i] != Cell::Fluid {
                    continue;
                }
                for q in 1..Q {
                    if let Some(nb) = neighbor(&self.grid, x, y, q) {
                        assert!(
                            self.grid.flags[nb] != Cell::Gas,
                            "closure violated: fluid touches gas at ({x},{y})"
                        );
                    }
                }
            }
        }
    }

    /// Fill fraction of a cell (fluid 1, gas 0, interface m/ρ).
    #[must_use]
    pub fn fill(&self, i: usize) -> f64 {
        match self.grid.flags[i] {
            Cell::Fluid => 1.0,
            Cell::Interface => {
                let rho: f64 = self.grid.f[i].iter().sum();
                (self.mass[i] / rho.max(1e-12)).clamp(0.0, 1.0)
            }
            _ => 0.0,
        }
    }

    /// The strict ledger: Σ_fluid Σf + Σ_interface m + carry.
    #[must_use]
    pub fn ledger_mass(&self) -> f64 {
        let mut total = self.carry;
        for i in 0..self.grid.nx * self.grid.ny {
            match self.grid.flags[i] {
                Cell::Fluid => total += self.grid.f[i].iter().sum::<f64>(),
                Cell::Interface => total += self.mass[i],
                _ => {}
            }
        }
        total
    }

    /// Count of connected fluid+interface components (4-connectivity)
    /// — the breaking-jet fragment counter.
    #[must_use]
    pub fn fragment_count(&self) -> usize {
        let (nx, ny) = (self.grid.nx, self.grid.ny);
        let wet = |i: usize| matches!(self.grid.flags[i], Cell::Fluid | Cell::Interface);
        let mut seen = vec![false; nx * ny];
        let mut count = 0;
        for start in 0..nx * ny {
            if !wet(start) || seen[start] {
                continue;
            }
            count += 1;
            let mut stack = vec![start];
            seen[start] = true;
            while let Some(i) = stack.pop() {
                let (x, y) = (i % nx, i / nx);
                for q in [1usize, 2, 3, 4] {
                    if let Some(nb) = neighbor(&self.grid, x, y, q)
                        && wet(nb)
                        && !seen[nb]
                    {
                        seen[nb] = true;
                        stack.push(nb);
                    }
                }
            }
        }
        count
    }

    /// Smoothed fill field + curvature-adjusted reference density for
    /// interface reconstruction at cell i.
    fn reference_density(&self, x: usize, y: usize) -> f64 {
        if self.sigma == 0.0 {
            return 1.0;
        }
        let kappa = self.curvature(x, y);
        self.sigma.mul_add(kappa / CS2, 1.0)
    }

    /// Fill-field value with the wall ghost per contact model.
    fn phi_at(&self, i: usize, from: usize) -> f64 {
        match self.grid.flags[i] {
            Cell::Wall => match self.contact {
                ContactModel::Neutral => self.fill(from),
                ContactModel::Wetting => 1.0,
            },
            _ => self.fill_smooth[i],
        }
    }

    /// Curvature of the smoothed fill field at (x, y): div(n̂) with
    /// n̂ = −∇φ/|∇φ| (outward from fluid), central differences.
    fn curvature(&self, x: usize, y: usize) -> f64 {
        let g = &self.grid;
        let i = g.idx(x, y);
        // Normal components at the four face neighbors via one-sided
        // gradients of the smoothed field; divergence by differencing.
        let phi = |dx: i32, dy: i32| -> f64 {
            let xx = offset_coord(x, dx, g.nx, g.periodic_x);
            let yy = offset_coord(y, dy, g.ny, g.periodic_y);
            match (xx, yy) {
                (Some(a), Some(b)) => self.phi_at(g.idx(a, b), i),
                _ => self.fill_smooth[i],
            }
        };
        let nhat = |dx: i32, dy: i32| -> [f64; 2] {
            let gx = (phi(dx + 1, dy) - phi(dx - 1, dy)) / 2.0;
            let gy = (phi(dx, dy + 1) - phi(dx, dy - 1)) / 2.0;
            let m = gx.hypot(gy).max(1e-9);
            [-gx / m, -gy / m]
        };
        let div = (nhat(1, 0)[0] - nhat(-1, 0)[0]) / 2.0 + (nhat(0, 1)[1] - nhat(0, -1)[1]) / 2.0;
        div.clamp(-1.0, 1.0)
    }

    fn refresh_fill(&mut self) {
        let (nx, ny) = (self.grid.nx, self.grid.ny);
        let raw: Vec<f64> = (0..nx * ny).map(|i| self.fill(i)).collect();
        for y in 0..ny {
            for x in 0..nx {
                let i = self.grid.idx(x, y);
                if self.grid.flags[i] == Cell::Wall {
                    self.fill_smooth[i] = raw[i];
                    continue;
                }
                let mut acc = raw[i];
                let mut count = 1.0;
                for q in 1..Q {
                    if let Some(nb) = neighbor(&self.grid, x, y, q)
                        && self.grid.flags[nb] != Cell::Wall
                    {
                        acc += raw[nb];
                        count += 1.0;
                    }
                }
                self.fill_smooth[i] = acc / count;
            }
        }
    }

    /// One free-surface step: collide, mass exchange + reconstruction
    /// + stream, conversion cascade with conservative redistribution.
    pub fn step(&mut self) {
        let (nx, ny) = (self.grid.nx, self.grid.ny);
        self.refresh_fill();
        self.grid.collide_into(&mut self.scratch);
        let post = std::mem::take(&mut self.scratch);
        // Mass exchange (pairwise antisymmetric, from POST populations)
        // and streaming with gas reconstruction.
        let fills: Vec<f64> = (0..nx * ny).map(|i| self.fill(i)).collect();
        let mut new_f = self.grid.f.clone();
        for y in 0..ny {
            for x in 0..nx {
                let i = self.grid.idx(x, y);
                let flag = self.grid.flags[i];
                if !matches!(flag, Cell::Fluid | Cell::Interface) {
                    continue;
                }
                // Mass exchange along all links (interface only; fluid
                // cells' Σf tracks their mass through plain streaming).
                if flag == Cell::Interface {
                    let mut dm = 0.0f64;
                    for q in 1..Q {
                        if let Some(nb) = neighbor(&self.grid, x, y, q) {
                            let w = match self.grid.flags[nb] {
                                Cell::Fluid => 1.0,
                                Cell::Interface => f64::midpoint(fills[i], fills[nb]),
                                _ => 0.0,
                            };
                            if w > 0.0 {
                                dm += w * (post[nb][OPP[q]] - post[i][q]);
                            }
                        }
                    }
                    self.mass[i] += dm;
                }
                // Stream (pull) with reconstruction from gas sources.
                let mm = self.grid.moments(i);
                for q in 0..Q {
                    let src = self.grid.source(x, y, q);
                    new_f[i][q] = match src {
                        Some(s) if self.grid.flags[s] == Cell::Gas => {
                            let rho_ref = self.reference_density(x, y);
                            let eq = equilibrium(rho_ref, mm.u[0], mm.u[1]);
                            eq[q] + eq[OPP[q]] - post[i][OPP[q]]
                        }
                        Some(s) if self.grid.flags[s] == Cell::Wall => post[i][OPP[q]],
                        Some(s) => post[s][q],
                        None => post[i][OPP[q]],
                    };
                }
            }
        }
        self.grid.f = new_f;
        self.scratch = post;
        self.apply_conversions();
    }

    fn apply_conversions(&mut self) {
        let (nx, ny) = (self.grid.nx, self.grid.ny);
        // Fluid mass is Σf; interface mass tracked. Conversions:
        let mut excess_pool = std::mem::take(&mut self.carry);
        let mut to_fluid = Vec::new();
        let mut to_gas = Vec::new();
        for i in 0..nx * ny {
            if self.grid.flags[i] != Cell::Interface {
                continue;
            }
            let rho: f64 = self.grid.f[i].iter().sum();
            if self.mass[i] > (1.0 + EPS) * rho {
                to_fluid.push(i);
            } else if self.mass[i] < -EPS * rho {
                to_gas.push(i);
            }
        }
        // Interface → fluid: excess to the pool; gas neighbors become
        // interface (closure).
        for &i in &to_fluid {
            let rho: f64 = self.grid.f[i].iter().sum();
            excess_pool += self.mass[i] - rho;
            self.grid.flags[i] = Cell::Fluid;
            self.mass[i] = rho;
            self.conversions.to_fluid += 1;
            let (x, y) = (i % nx, i / nx);
            for q in 1..Q {
                if let Some(nb) = neighbor(&self.grid, x, y, q)
                    && self.grid.flags[nb] == Cell::Gas
                {
                    // Initialize from the average of wet neighbors.
                    let (nbx, nby) = (nb % nx, nb / nx);
                    let mut rho_avg = 0.0;
                    let mut u_avg = [0.0f64; 2];
                    let mut cnt = 0.0;
                    for q2 in 1..Q {
                        if let Some(nn) = neighbor(&self.grid, nbx, nby, q2)
                            && matches!(self.grid.flags[nn], Cell::Fluid | Cell::Interface)
                        {
                            let m2 = self.grid.moments(nn);
                            rho_avg += m2.rho;
                            u_avg[0] += m2.u[0];
                            u_avg[1] += m2.u[1];
                            cnt += 1.0;
                        }
                    }
                    if cnt > 0.0 {
                        rho_avg /= cnt;
                        u_avg[0] /= cnt;
                        u_avg[1] /= cnt;
                    } else {
                        rho_avg = 1.0;
                    }
                    self.grid.f[nb] = equilibrium(rho_avg, u_avg[0], u_avg[1]);
                    self.grid.flags[nb] = Cell::Interface;
                    self.mass[nb] = 0.0;
                    self.conversions.gas_to_interface += 1;
                }
            }
        }
        // Interface → gas: deficit to the pool; fluid neighbors become
        // interface (their Σf IS their mass — ledger unchanged).
        for &i in &to_gas {
            if self.grid.flags[i] != Cell::Interface {
                continue; // may have been re-flagged by the cascade
            }
            excess_pool += self.mass[i];
            self.grid.flags[i] = Cell::Gas;
            self.mass[i] = 0.0;
            self.conversions.to_gas += 1;
            let (x, y) = (i % nx, i / nx);
            for q in 1..Q {
                if let Some(nb) = neighbor(&self.grid, x, y, q)
                    && self.grid.flags[nb] == Cell::Fluid
                {
                    self.grid.flags[nb] = Cell::Interface;
                    self.mass[nb] = self.grid.f[nb].iter().sum();
                    self.conversions.fluid_to_interface += 1;
                }
            }
        }
        // Conservative redistribution of the pool over interface cells.
        let interfaces: Vec<usize> = (0..nx * ny)
            .filter(|&i| self.grid.flags[i] == Cell::Interface)
            .collect();
        if interfaces.is_empty() {
            self.carry = excess_pool;
        } else {
            let share = excess_pool / interfaces.len() as f64;
            for &i in &interfaces {
                self.mass[i] += share;
            }
        }
    }
}

fn offset_coord(coord: usize, delta: i32, n: usize, periodic: bool) -> Option<usize> {
    assert!(n > 0, "grid dimension must be positive");
    let mut out = coord;
    if delta >= 0 {
        for _ in 0..delta {
            if out + 1 == n {
                if periodic {
                    out = 0;
                } else {
                    return None;
                }
            } else {
                out += 1;
            }
        }
    } else {
        for _ in 0..(-delta) {
            if out == 0 {
                if periodic {
                    out = n - 1;
                } else {
                    return None;
                }
            } else {
                out -= 1;
            }
        }
    }
    Some(out)
}

/// Neighbor cell index in direction q (None across non-periodic
/// boundaries).
fn neighbor(grid: &Grid, x: usize, y: usize, q: usize) -> Option<usize> {
    let (ex, ey) = E[q];
    let xx = offset_coord(x, ex, grid.nx, grid.periodic_x)?;
    let yy = offset_coord(y, ey, grid.ny, grid.periodic_y)?;
    Some(grid.idx(xx, yy))
}

/// A closed-box dam-break fixture: walls on all four sides, a fluid
/// column of `a × 2a` cells in the lower-left corner, gas elsewhere,
/// gravity `g` pointing down.
#[must_use]
pub fn dam_break(
    nx: usize,
    ny: usize,
    a: usize,
    g: f64,
    sigma: f64,
    contact: ContactModel,
) -> FreeSurface {
    let mut grid = Grid::uniform(nx, ny, 0.55);
    grid.periodic_x = false;
    grid.periodic_y = false;
    grid.g = [0.0, -g];
    for i in 0..nx * ny {
        grid.flags[i] = Cell::Gas;
    }
    for x in 0..nx {
        let b = grid.idx(x, 0);
        grid.flags[b] = Cell::Wall;
        let t = grid.idx(x, ny - 1);
        grid.flags[t] = Cell::Wall;
    }
    for y in 0..ny {
        let l = grid.idx(0, y);
        grid.flags[l] = Cell::Wall;
        let r = grid.idx(nx - 1, y);
        grid.flags[r] = Cell::Wall;
    }
    for y in 1..=(2 * a).min(ny - 2) {
        for x in 1..=a.min(nx - 2) {
            let i = grid.idx(x, y);
            grid.flags[i] = Cell::Fluid;
        }
    }
    FreeSurface::new(grid, sigma, contact)
}

/// Surge-front x position (rightmost wet cell in the bottom fluid
/// row), in cells from the left wall.
#[must_use]
pub fn surge_front(fs: &FreeSurface) -> usize {
    let mut front = 0;
    for x in 1..fs.grid.nx - 1 {
        let i = fs.grid.idx(x, 1);
        if matches!(fs.grid.flags[i], Cell::Fluid | Cell::Interface) {
            front = x;
        }
    }
    front
}
