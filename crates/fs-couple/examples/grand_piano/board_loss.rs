//! Soundboard loss from material data (bead pexfp.373), declared in the board
//! geometry and applied on the complete certified source slice.
//!
//! A `damping,material` board names a loss factor for every share of its
//! elastic stiffness: the panel's four orthotropic bending parts (`panel-l`,
//! `panel-r`, `panel-lr`, `panel-12`) and each stiffener group's bending and
//! torsion (`<group>-bending`, `<group>-torsion`), every stiffener in exactly
//! one group. Rows (after the header, anywhere among the others):
//!
//! - `stiffener-group,name,index,...` stiffener rows by their 0-based order
//! - `loss-share,share,eta,f_lo_hz,f_hi_hz,strict|extrapolate-frequency,receipt`
//! - `loss-unknown,share,reason` (the share's loss is omitted: a lower bound)
//! - `loss-sensitivity,label,from_hz,to_hz,ratio,source` (printed only)
//! - `loss-proxy,share,eta_lo,eta_hi,f_lo_hz,f_hi_hz,composition-proxy,receipt`
//!   (printed only): a bracket of a `loss-unknown` share's loss from other
//!   species, never entering the damping
//!
//! Each loss factor is constant (hysteretic) over the frequencies its source
//! measured, `[f_lo, f_hi]`. A source mode outside that domain refuses by name
//! unless the share declares `extrapolate-frequency`; then the same constant is
//! taken there and every receipt says so, with the size of the possible error:
//! the loss each mode would have if extrapolated shares rose with frequency as
//! a `loss-sensitivity` law measured (`(f/from)^p`, `p = ln ratio / ln(to/from)`,
//! clamped to `[from, to]`). Nothing falls back to an authored damping ratio.
//!
//! The shares' modal strain energies are projected once by the plate loss
//! owner (`fs_couple::plate_viscosity::hysteretic`); the reduction consumes the
//! diagonal source damping `c_a = eta_a omega_a`, and the receipt prints each
//! mode's energy split (including the shares whose loss is unknown) and the
//! measured off-diagonal coupling the diagonal neglects. Where that coupling is
//! not small, each coupled cluster's exact eigenvalues are solved and every
//! mode's diagonal loss factor is printed against its exact one, under the
//! declared [`DIAGONAL_TOLERANCE`] (ruling 6431).
use fs_couple::plate_viscosity::hysteretic::{self, NeighbourCoupling, OrthotropicBendingParts, StiffnessShare};
use fs_plate::{PlateMesh, PlateModel, Stiffener};
use std::collections::BTreeSet;
use std::fmt::Write as _;

/// The panel's four orthotropic bending shares, in projection order.
pub const PANEL_SHARES: [&str; 4] = ["panel-l", "panel-r", "panel-lr", "panel-12"];
/// Upper bound on stiffener groups (each adds two shares).
pub const MAX_GROUPS: usize = 16;
/// Upper bound on printed sensitivity laws.
pub const MAX_SENSITIVITY: usize = 4;
/// Upper bound on composition-proxy brackets.
pub const MAX_PROXIES: usize = 8;
/// Two modes whose `|H_ab| / |omega_a^2 - omega_b^2|` reaches this join one
/// cluster for the exact check of the diagonal (ruling 6431).
pub const CLUSTER_PAIR_RATIO: f64 = 0.03;
/// The declared bound on `|eta_exact / eta_diagonal - 1|` (ruling 6431): below
/// every uncertainty the receipt already prints. A mode above it means the full
/// source loss matrix must enter the reduction before the board is compared
/// with any reference.
pub const DIAGONAL_TOLERANCE: f64 = 0.05;

/// One share's declared loss.
#[derive(Clone, Debug, PartialEq)]
pub enum ShareLoss {
    /// A sourced loss factor holding over `domain_hz`; `extrapolate` takes the
    /// same value at modes outside it (declared, and printed in every receipt).
    Claim { eta: f64, domain_hz: (f64, f64), extrapolate: bool, receipt: String },
    /// No source: the share's loss is omitted and every receipt says so.
    Unknown { reason: String },
}

/// A measured rise of loss with frequency, used only for the printed
/// sensitivity of extrapolated shares.
#[derive(Clone, Debug, PartialEq)]
pub struct RiseLaw {
    pub label: String,
    pub from_hz: f64,
    pub to_hz: f64,
    pub ratio: f64,
    pub source: String,
}

impl RiseLaw {
    /// Loss multiplier at `hz`: one at or below `from_hz`, the measured power
    /// law up to `to_hz`, and held at `ratio` above it.
    #[must_use]
    pub fn factor(&self, hz: f64) -> f64 {
        if hz <= self.from_hz { return 1.0; }
        let p = self.ratio.ln() / (self.to_hz / self.from_hz).ln();
        (hz.min(self.to_hz) / self.from_hz).powf(p)
    }
}

/// A composition-proxy bracket of one `loss-unknown` share's loss factor,
/// measured on other species over `domain_hz`: printed beside the lower bound,
/// never part of the board's damping.
#[derive(Clone, Debug, PartialEq)]
pub struct ProxyBracket {
    pub share: String,
    pub eta: (f64, f64),
    pub domain_hz: (f64, f64),
    pub receipt: String,
}

/// Every loss row of one `damping,material` board.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct LossDeclaration {
    shares: Vec<(String, ShareLoss)>,
    groups: Vec<(String, Vec<usize>)>,
    sensitivity: Vec<RiseLaw>,
    proxies: Vec<ProxyBracket>,
}

fn number(fields: &[&str], i: usize, what: &str) -> Result<f64, String> {
    let x = fields.get(i).ok_or_else(|| format!("{what} is missing"))?.parse::<f64>()
        .map_err(|_| format!("invalid {what}"))?;
    if x.is_finite() { Ok(x) } else { Err(format!("nonfinite {what}")) }
}
fn text(fields: &[&str], from: usize, what: &str) -> Result<String, String> {
    // The geometry reader trims every comma-separated field; rejoin with ", ".
    let joined = fields.get(from..).map(|rest| rest.join(", ")).unwrap_or_default();
    if joined.trim().is_empty() { Err(format!("{what} needs its source text")) } else { Ok(joined) }
}
fn slug(name: &str) -> bool {
    !name.is_empty() && name.len() <= 64 && name.bytes().all(|b| b.is_ascii_lowercase() || b.is_ascii_digit() || b == b'-' || b == b'_')
}

impl LossDeclaration {
    /// Read one loss row; `Ok(false)` when the row is not a loss row.
    ///
    /// # Errors
    /// Malformed fields, negative or nonfinite loss factors, empty or reversed
    /// domains, unknown flags and repeated names refuse.
    pub fn read_row(&mut self, f: &[&str]) -> Result<bool, String> {
        match f[0] {
            "stiffener-group" => {
                let name = f.get(1).copied().unwrap_or("");
                if !slug(name) || name.starts_with("panel") || self.groups.iter().any(|(n, _)| n == name) {
                    return Err("stiffener-group needs a unique lowercase name that is not a panel share".into());
                }
                if self.groups.len() >= MAX_GROUPS { return Err("too many stiffener groups".into()); }
                let members = f[2..].iter().map(|v| v.parse::<usize>().map_err(|_| "invalid stiffener index".to_string()))
                    .collect::<Result<Vec<_>, _>>()?;
                if members.is_empty() { return Err("stiffener-group needs at least one stiffener".into()); }
                self.groups.push((name.to_string(), members));
            }
            "loss-share" => {
                let share = f.get(1).copied().unwrap_or("");
                let eta = number(f, 2, "loss factor")?;
                let (lo, hi) = (number(f, 3, "frequency domain lower edge")?, number(f, 4, "frequency domain upper edge")?);
                let extrapolate = match f.get(5).copied() {
                    Some("strict") => false,
                    Some("extrapolate-frequency") => true,
                    _ => return Err("loss-share flag must be strict or extrapolate-frequency".into()),
                };
                if !(eta >= 0.0 && eta < 1.0) || !(lo > 0.0 && lo <= hi) {
                    return Err(format!("loss-share {share}: loss factor in [0,1) and a positive frequency domain lo <= hi"));
                }
                self.push(share, ShareLoss::Claim { eta, domain_hz: (lo, hi), extrapolate, receipt: text(f, 6, "loss-share")? })?;
            }
            "loss-unknown" => {
                let share = f.get(1).copied().unwrap_or("");
                self.push(share, ShareLoss::Unknown { reason: text(f, 2, "loss-unknown")? })?;
            }
            "loss-sensitivity" => {
                if self.sensitivity.len() >= MAX_SENSITIVITY { return Err("too many loss-sensitivity laws".into()); }
                let label = f.get(1).copied().unwrap_or("");
                let law = RiseLaw { label: label.to_string(), from_hz: number(f, 2, "sensitivity start")?,
                    to_hz: number(f, 3, "sensitivity end")?, ratio: number(f, 4, "sensitivity ratio")?,
                    source: text(f, 5, "loss-sensitivity")? };
                if !slug(label) || self.sensitivity.iter().any(|l| l.label == label)
                    || !(law.from_hz > 0.0 && law.to_hz > law.from_hz && law.ratio > 0.0) {
                    return Err("loss-sensitivity needs a unique label, 0 < from < to and a positive ratio".into());
                }
                self.sensitivity.push(law);
            }
            "loss-proxy" => {
                let share = f.get(1).copied().unwrap_or("");
                let eta = (number(f, 2, "proxy loss factor lower edge")?, number(f, 3, "proxy loss factor upper edge")?);
                let domain_hz = (number(f, 4, "proxy frequency domain lower edge")?,
                    number(f, 5, "proxy frequency domain upper edge")?);
                if f.get(6).copied() != Some("composition-proxy") {
                    return Err(format!("loss-proxy {share}: the flag must be composition-proxy"));
                }
                if !(eta.0 >= 0.0 && eta.0 <= eta.1 && eta.1 < 1.0) || !(domain_hz.0 > 0.0 && domain_hz.0 <= domain_hz.1) {
                    return Err(format!("loss-proxy {share}: a bracket 0 <= lo <= hi < 1 and a positive frequency domain lo <= hi"));
                }
                if !slug(share) || self.proxies.iter().any(|p| p.share == share) || self.proxies.len() >= MAX_PROXIES {
                    return Err(format!("loss-proxy '{share}' must name one share once, at most {MAX_PROXIES} brackets"));
                }
                self.proxies.push(ProxyBracket { share: share.to_string(), eta, domain_hz, receipt: text(f, 7, "loss-proxy")? });
            }
            _ => return Ok(false),
        }
        Ok(true)
    }

    fn push(&mut self, share: &str, loss: ShareLoss) -> Result<(), String> {
        if !slug(share) || self.shares.iter().any(|(n, _)| n == share) {
            return Err(format!("loss share '{share}' must be a unique lowercase name"));
        }
        self.shares.push((share.to_string(), loss));
        Ok(())
    }

    /// True when the declaration holds any loss row.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.shares.is_empty() && self.groups.is_empty() && self.sensitivity.is_empty() && self.proxies.is_empty()
    }

    /// The shares in projection order: the four panel parts, then each group's
    /// bending and torsion in declared group order.
    #[must_use]
    pub fn share_names(&self) -> Vec<String> {
        PANEL_SHARES.iter().map(|s| (*s).to_string())
            .chain(self.groups.iter().flat_map(|(g, _)| [format!("{g}-bending"), format!("{g}-torsion")]))
            .collect()
    }

    /// Every stiffener in exactly one group and every share declared exactly
    /// once, with nothing undeclared or extra.
    ///
    /// # Errors
    /// An ungrouped, repeated or out-of-range stiffener, a missing or extra
    /// share declaration refuse.
    pub fn validate(&self, stiffeners: usize) -> Result<(), String> {
        let mut seen = BTreeSet::new();
        for (group, members) in &self.groups {
            for &m in members {
                if m >= stiffeners || !seen.insert(m) {
                    return Err(format!("stiffener-group {group}: stiffener {m} is out of range or already grouped"));
                }
            }
        }
        if seen.len() != stiffeners {
            let missing = (0..stiffeners).find(|m| !seen.contains(m)).unwrap_or(0);
            return Err(format!("stiffener {missing} belongs to no stiffener-group, so its loss is undeclared"));
        }
        let names = self.share_names();
        for name in &names {
            if !self.shares.iter().any(|(n, _)| n == name) {
                return Err(format!("loss share '{name}' is undeclared: give it a loss-share or loss-unknown row"));
            }
        }
        if let Some((extra, _)) = self.shares.iter().find(|(n, _)| !names.contains(n)) {
            return Err(format!("loss share '{extra}' names no panel part or stiffener group"));
        }
        for proxy in &self.proxies {
            if !matches!(self.shares.iter().find(|(n, _)| *n == proxy.share), Some((_, ShareLoss::Unknown { .. }))) {
                return Err(format!("loss-proxy '{}' brackets only a loss-unknown share; a claimed or undeclared share refuses",
                    proxy.share));
            }
        }
        Ok(())
    }

    fn loss(&self, share: &str) -> &ShareLoss {
        &self.shares.iter().find(|(n, _)| n == share).expect("validated share").1
    }
}

/// One source mode's material loss and its receipt.
#[derive(Clone, Debug, PartialEq)]
pub struct ModeLoss {
    pub frequency_hz: f64,
    /// `eta_a = H_aa / omega_a^2`.
    pub loss_factor: f64,
    /// `eta_a / 2`, the ratio the authored board stated as one number.
    pub damping_ratio: f64,
    /// Each share's fraction of the mode's strain energy, in share order.
    pub fractions: Vec<f64>,
    /// The sum of the fractions of shares whose loss is unknown: the missing
    /// loss is at most their unknown loss factor times this.
    pub unknown_fraction: f64,
    /// Shares taken outside their measured frequencies at this mode.
    pub extrapolated: Vec<usize>,
    pub max_coupling: f64,
    pub neighbours: [Option<NeighbourCoupling>; 2],
    /// The loss factor if the extrapolated shares rose as each sensitivity law.
    pub sensitivity: Vec<f64>,
    /// The exact loss factor of the coupled cluster's eigenvalue matched to
    /// this mode (`loss_factor` itself for a mode in no cluster).
    pub exact_loss_factor: f64,
    /// The mode's coupled cluster, if any, and its eigenvector overlap.
    pub cluster: Option<usize>,
    pub overlap: f64,
    /// Second-order estimate of the relative change the pairs left between
    /// clusters would make to the exact loss factor.
    pub neglected: f64,
    /// The loss factor if every composition-proxy share took the lower, then
    /// the upper edge of its bracket (printed only; `None` without brackets).
    pub proxy_bracket: Option<(f64, f64)>,
}

impl ModeLoss {
    /// `eta_exact / eta_diagonal - 1`: how far the reduction's diagonal is off.
    #[must_use]
    pub fn diagonal_error(&self) -> f64 { self.exact_loss_factor / self.loss_factor - 1.0 }
}

/// The applied declaration on every source mode.
#[derive(Clone, Debug, PartialEq)]
pub struct MaterialLossReceipt {
    pub shares: Vec<(String, ShareLoss)>,
    pub sensitivity: Vec<RiseLaw>,
    pub proxies: Vec<ProxyBracket>,
    pub modes: Vec<ModeLoss>,
    /// The coupled clusters of the exact check, and the largest pair ratio left
    /// between them.
    pub clusters: Vec<Vec<usize>>,
    pub max_outside_ratio: f64,
    pub caveats: Vec<String>,
}

/// Per-band summary against a reference damping ratio.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BandSummary {
    pub band_hz: (f64, f64),
    pub modes: usize,
    pub mean_ratio: f64,
    pub min_ratio: f64,
    pub max_ratio: f64,
    /// The mean damping ratio with every composition-proxy share at the lower,
    /// then the upper edge of its bracket.
    pub proxy_mean_ratio: Option<(f64, f64)>,
}

impl MaterialLossReceipt {
    /// Damping-ratio statistics per band `[edges[i], edges[i+1])`; empty
    /// bands are reported with zero modes and NaN statistics.
    #[must_use]
    pub fn band_summary(&self, edges_hz: &[f64]) -> Vec<BandSummary> {
        edges_hz.windows(2).map(|w| {
            let band: Vec<&ModeLoss> = self.modes.iter().filter(|m| m.frequency_hz >= w[0] && m.frequency_hz < w[1]).collect();
            let ratios: Vec<f64> = band.iter().map(|m| m.damping_ratio).collect();
            let n = ratios.len();
            let proxy: Vec<(f64, f64)> = band.iter().filter_map(|m| m.proxy_bracket).collect();
            let mean = |edge: fn(&(f64, f64)) -> f64| proxy.iter().map(edge).sum::<f64>() / (2.0 * n as f64);
            BandSummary { band_hz: (w[0], w[1]), modes: n,
                mean_ratio: if n == 0 { f64::NAN } else { ratios.iter().sum::<f64>() / n as f64 },
                min_ratio: ratios.iter().copied().fold(f64::NAN, f64::min),
                max_ratio: ratios.iter().copied().fold(f64::NAN, f64::max),
                proxy_mean_ratio: (n > 0 && proxy.len() == n).then(|| (mean(|b| b.0), mean(|b| b.1))) }
        }).collect()
    }

    /// The modes whose diagonal loss factor is off its cluster-exact one by more
    /// than [`DIAGONAL_TOLERANCE`], and the largest `|diagonal_error|` with its mode.
    #[must_use]
    pub fn diagonal_check(&self) -> (Vec<usize>, f64, usize) {
        let above = (0..self.modes.len()).filter(|&a| self.modes[a].diagonal_error().abs() > DIAGONAL_TOLERANCE).collect();
        let (worst, mode) = self.modes.iter().enumerate().map(|(a, m)| (m.diagonal_error().abs(), a))
            .fold((0.0_f64, 0), |best, x| if x.0 > best.0 { x } else { best });
        (above, worst, mode)
    }

    /// One line per mode, tab-separated, with a header naming every column.
    #[must_use]
    pub fn tsv(&self) -> String {
        let mut out = String::from("mode\tfrequency_hz\tloss_factor\tdamping_ratio");
        for (name, _) in &self.shares { let _ = write!(out, "\tf_{name}"); }
        out.push_str("\tunknown_fraction\textrapolated\tmax_coupling\tlower_neighbour\tlower_ratio\tupper_neighbour\tupper_ratio");
        for law in &self.sensitivity { let _ = write!(out, "\tloss_factor_if_{}", law.label); }
        out.push_str("\texact_loss_factor\tdiagonal_error\tcluster\toverlap\tneglected\tproxy_loss_factor_lo\tproxy_loss_factor_hi\n");
        for (a, m) in self.modes.iter().enumerate() {
            let _ = write!(out, "{a}\t{:.6}\t{:.6e}\t{:.6e}", m.frequency_hz, m.loss_factor, m.damping_ratio);
            for f in &m.fractions { let _ = write!(out, "\t{f:.6e}"); }
            let ext = m.extrapolated.iter().map(|&c| self.shares[c].0.as_str()).collect::<Vec<_>>().join(";");
            let _ = write!(out, "\t{:.6e}\t{}\t{:.3e}", m.unknown_fraction, if ext.is_empty() { "-" } else { &ext }, m.max_coupling);
            for side in &m.neighbours {
                match side {
                    Some(c) => { let _ = write!(out, "\t{}\t{:.3e}", c.mode, c.ratio); }
                    None => out.push_str("\t-\t-"),
                }
            }
            for s in &m.sensitivity { let _ = write!(out, "\t{s:.6e}"); }
            let cluster = m.cluster.map_or_else(|| "-".to_string(), |k| k.to_string());
            let _ = write!(out, "\t{:.6e}\t{:+.3e}\t{cluster}\t{:.3}\t{:.1e}", m.exact_loss_factor, m.diagonal_error(),
                m.overlap, m.neglected);
            let _ = match m.proxy_bracket {
                Some((lo, hi)) => writeln!(out, "\t{lo:.6e}\t{hi:.6e}"),
                None => writeln!(out, "\t-\t-"),
            };
        }
        out
    }
}

/// What a `damping,material` board supplies to resolve its loss.
pub(super) struct BoardParts<'a> {
    pub mesh: &'a PlateMesh,
    /// Per triangle: E_L, E_R, nu_LR, G_LR [Pa], thickness [m], grain [rad].
    pub panel: &'a [[f64; 6]],
    pub stiffeners: &'a [Stiffener],
}

/// Resolve the declaration on the source modes: the shares' strain energies
/// projected once by the owner, each share's loss applied by its domain rule,
/// and the diagonal source damping `c_a` [1/s] the reduction consumes.
///
/// # Errors
/// A mode outside a strict share's measured frequencies refuses by name, as do
/// assembly, projection and non-dissipative results; nothing falls back.
pub(super) fn resolve(declaration: &LossDeclaration, parts: &BoardParts<'_>, model: &PlateModel,
    phi: &[Vec<f64>], omega: &[f64], assemble: &dyn Fn(&[Stiffener]) -> Result<PlateModel, String>)
    -> Result<(Vec<f64>, MaterialLossReceipt), String> {
    declaration.validate(parts.stiffeners.len())?;
    let names = declaration.share_names();
    // Panel bending parts per triangle, rotated as fs-plate rotates the section.
    let panel_parts = parts.panel.iter().map(|p| OrthotropicBendingParts::new(p[0], p[1], p[2], p[3], p[4], p[5])
        .map_err(|e| e.to_string())).collect::<Result<Vec<_>, _>>()?;
    let pick = |f: fn(&OrthotropicBendingParts) -> [f64; 9]| panel_parts.iter().map(f).collect::<Vec<_>>();
    let tensors = [pick(|p| p.longitudinal), pick(|p| p.radial), pick(|p| p.shear), pick(|p| p.poisson)];
    // Each group's bending (E doubled) and torsion (G doubled) assemblies.
    let mut doubled = Vec::with_capacity(2 * declaration.groups.len());
    for (_, members) in &declaration.groups {
        for torsion in [false, true] {
            let edited: Vec<Stiffener> = parts.stiffeners.iter().enumerate().map(|(i, s)| {
                let mut s = s.clone();
                if members.contains(&i) { if torsion { s.g *= 2.0 } else { s.e *= 2.0 } }
                s
            }).collect();
            doubled.push(assemble(&edited)?);
        }
    }
    let mut shares: Vec<(&str, StiffnessShare<'_>)> = Vec::with_capacity(names.len());
    for (name, t) in PANEL_SHARES.iter().zip(&tensors) { shares.push((name, StiffnessShare::Panel(t))); }
    for (name, plus) in names[4..].iter().zip(&doubled) {
        shares.push((name, StiffnessShare::Difference { plus, minus: model }));
    }
    let energies = hysteretic::project_energy_shares_offline(parts.mesh, model, &shares, phi, omega)
        .map_err(|e| format!("board material loss projection: {e}"))?;
    // Per-mode loss table by each share's domain rule.
    let n = omega.len();
    let mut eta = vec![vec![0.0; n]; names.len()];
    let mut extrapolated = vec![Vec::new(); n];
    for (c, name) in names.iter().enumerate() {
        if let ShareLoss::Claim { eta: value, domain_hz: (lo, hi), extrapolate, receipt } = declaration.loss(name) {
            for a in 0..n {
                let hz = omega[a] / std::f64::consts::TAU;
                if hz < *lo || hz > *hi {
                    if !extrapolate {
                        return Err(format!("board mode {a} ({hz:.2} Hz) lies outside the {lo:.2}-{hi:.2} Hz the '{name}' \
                            loss claim was measured over ({receipt}); declare extrapolate-frequency for that share or the \
                            board refuses; no authored damping replaces it"));
                    }
                    extrapolated[a].push(c);
                }
                eta[c][a] = *value;
            }
        }
    }
    let loss = energies.hysteretic(&eta).map_err(|e| format!("board material loss: {e}"))?;
    let exact = energies.cluster_exact(&eta, CLUSTER_PAIR_RATIO).map_err(|e| format!("board material loss check: {e}"))?;
    let unknown: Vec<usize> = (0..names.len()).filter(|&c| matches!(declaration.loss(&names[c]), ShareLoss::Unknown { .. })).collect();
    let proxied: Vec<(usize, &ProxyBracket)> = declaration.proxies.iter()
        .map(|p| (names.iter().position(|n| *n == p.share).expect("validated proxy share"), p)).collect();
    let modes = (0..n).map(|a| {
        let hz = omega[a] / std::f64::consts::TAU;
        let fractions: Vec<f64> = (0..names.len()).map(|c| energies.fraction(c, a)).collect();
        let sensitivity = declaration.sensitivity.iter().map(|law| {
            (0..names.len()).map(|c| {
                let scale = if extrapolated[a].contains(&c) { law.factor(hz) } else { 1.0 };
                eta[c][a] * scale * fractions[c]
            }).sum()
        }).collect();
        // Printed only: the declared loss plus each bracketed share at one edge.
        let proxy_bracket = (!proxied.is_empty()).then(|| {
            let add = |edge: fn(&ProxyBracket) -> f64| proxied.iter().map(|(c, p)| edge(p) * fractions[*c]).sum::<f64>();
            (loss.loss_factor[a] + add(|p| p.eta.0), loss.loss_factor[a] + add(|p| p.eta.1))
        });
        ModeLoss { frequency_hz: hz, loss_factor: loss.loss_factor[a], damping_ratio: loss.loss_factor[a] / 2.0,
            unknown_fraction: unknown.iter().map(|&c| fractions[c]).sum(), fractions,
            extrapolated: extrapolated[a].clone(), max_coupling: loss.max_coupling[a],
            neighbours: loss.neighbours[a], sensitivity, exact_loss_factor: exact.loss_factor[a],
            cluster: exact.cluster_of[a], overlap: exact.overlap[a], neglected: exact.neglected[a], proxy_bracket }
    }).collect::<Vec<ModeLoss>>();
    let mut caveats = Vec::new();
    for &c in &unknown {
        if let ShareLoss::Unknown { reason } = declaration.loss(&names[c]) {
            let most = modes.iter().map(|m| m.fractions[c]).fold(f64::NEG_INFINITY, f64::max);
            caveats.push(format!("'{}' loss unknown ({reason}); board damping is a lower bound: the missing loss is its \
                unknown loss factor times f_{} (largest {most:.3e})", names[c], names[c]));
        }
    }
    for (c, name) in names.iter().enumerate() {
        let outside: Vec<f64> = modes.iter().filter(|m| m.extrapolated.contains(&c)).map(|m| m.frequency_hz).collect();
        if let (Some(lo), ShareLoss::Claim { eta, domain_hz, .. }) = (outside.iter().copied().reduce(f64::min), declaration.loss(name)) {
            let hi = outside.iter().copied().fold(lo, f64::max);
            caveats.push(format!("'{name}' loss factor {eta:.5} taken outside its measured {:.2}-{:.2} Hz at {} modes \
                ({lo:.1}-{hi:.1} Hz): constant loss is assumed there; see the sensitivity columns", domain_hz.0,
                domain_hz.1, outside.len()));
        }
    }
    for (c, p) in &proxied {
        let most = modes.iter().map(|m| m.fractions[*c]).fold(f64::NEG_INFINITY, f64::max);
        let outside = modes.iter().filter(|m| m.frequency_hz < p.domain_hz.0 || m.frequency_hz > p.domain_hz.1).count();
        caveats.push(format!("composition-proxy, printed only: '{}' bracketed by loss factor {:.5}-{:.5} ({}), measured \
            at {:.0}-{:.0} Hz and taken constant at the {outside} of {} modes outside it (largest f_{} {most:.3e}); the \
            board's damping keeps it at the lower bound", p.share, p.eta.0, p.eta.1, p.receipt, p.domain_hz.0, p.domain_hz.1,
            modes.len(), p.share));
    }
    // Off-diagonal loss: sqrt(eta_c,a eta_c,b) is exact only where a share's loss is the same at both modes.
    let varying = (0..names.len()).filter(|&c| eta[c].iter().any(|v| *v != eta[c][0])).count();
    caveats.push(if varying == 0 {
        "off-diagonal loss H_ab takes sqrt(eta_c,a eta_c,b) as share c's loss between two modes: exact here, since \
            every share's loss factor is the same at every mode".to_string()
    } else {
        format!("off-diagonal loss H_ab takes sqrt(eta_c,a eta_c,b) as share c's loss between two modes: a stated \
            choice for the {varying} shares whose loss factor differs between modes")
    });
    // Shares in projection order, matching every mode's fraction columns.
    let shares = names.iter().map(|name| (name.clone(), declaration.loss(name).clone())).collect();
    let mut receipt = MaterialLossReceipt { shares, sensitivity: declaration.sensitivity.clone(),
        proxies: declaration.proxies.clone(), modes,
        clusters: exact.clusters, max_outside_ratio: exact.max_outside_ratio, caveats };
    let (above, worst, mode) = receipt.diagonal_check();
    let largest = receipt.clusters.iter().map(Vec::len).max().unwrap_or(0);
    let left_out = receipt.modes.iter().map(|m| m.neglected).fold(0.0_f64, f64::max);
    let mut check = format!("diagonal loss factor against the exact eigenvalues of each coupled cluster (modes joined \
        where |H_ab|/|w_a^2-w_b^2| >= {CLUSTER_PAIR_RATIO}: {} clusters, largest {largest} modes; between clusters the \
        largest pair ratio is {:.3e} and the largest second-order estimate of what is left out {left_out:.1e}): \
        largest |change| {:.2}% at mode {mode} ({:.1} Hz); {} modes above the declared {:.0}% tolerance",
        receipt.clusters.len(), receipt.max_outside_ratio, 100.0 * worst, receipt.modes[mode].frequency_hz, above.len(),
        100.0 * DIAGONAL_TOLERANCE);
    if !above.is_empty() {
        let list = above.iter().map(|&a| format!("{a} ({:.1} Hz, {:+.1}%)", receipt.modes[a].frequency_hz,
            100.0 * receipt.modes[a].diagonal_error())).collect::<Vec<_>>().join(", ");
        let _ = write!(check, " ({list}): the full source loss matrix must enter the reduction before this board is \
            compared with any reference");
    }
    receipt.caveats.push(check);
    Ok((loss.source_damping, receipt))
}
