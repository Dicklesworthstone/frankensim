//! Declared-band scalar radiation bakes using the existing solver and SH owner.
//!
//! Complete means two successive refinement changes meet the declared relative
//! tolerance at EVERY declared frequency, after at least three solved rungs.
//! It is finite-grid discretization evidence, not a uniform-in-frequency error
//! enclosure, a validated instrument, or a passive rational realization.
use super::{RadiationLoadRow, bake_radiation_load, bake_radiation_load_iterative};
use crate::helmholtz::iterative::{IterativeError, IterativeOptions};
use crate::helmholtz::{
    HelmholtzError, MAX_DENSE_PANELS, MAX_SH_DEGREE, MIN_PANELS_PER_WAVELENGTH, Medium,
};
use crate::panel3d::SpherePanels;
use fs_blake3::{ContentHash, hash_domain};
use fs_math::{c64::C64, det};
use std::time::Instant;

/// The physical reference surface; a body-velocity mode is explicitly not a duct port.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum ReferencePlane {
    DuctMouth {
        origin_m: [f64; 3],
        outward_normal: [f64; 3],
    },
    WholeSurface,
}

/// Scalar plane mode is implemented here. Higher-mode projection waits on M2.1-W3c.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum ModeContent {
    Plane,
    Azimuthal { modes: Vec<(i32, u32)> },
}

/// Immutable validated declaration, whose identity is computed before any solve.
#[derive(Clone, Debug)]
pub struct BakeDeclaration {
    geometry_hash: ContentHash,
    band_hz: (f64, f64),
    frequencies_hz: Vec<f64>,
    operator_tol: f64,
    reference_plane: ReferencePlane,
    modes: ModeContent,
    sh_degree: usize,
    medium: Medium,
    hash: ContentHash,
}

#[derive(Clone, Debug, PartialEq)]
pub enum BakeError {
    Invalid(&'static str),
    Solver(HelmholtzError),
    DeclarationMismatch,
    Dependency {
        route: &'static str,
    },
    NonPassive {
        at_hz: f64,
        margin: f64,
    },
    OutOfBand {
        at_hz: f64,
        proven_hz: Option<(f64, f64)>,
    },
    /// A refusal of the M2.5-W5 iterative route, with its diagnostic.
    IterativeRoute {
        at_hz: f64,
        what: String,
    },
    Cancelled,
}
impl From<HelmholtzError> for BakeError {
    fn from(e: HelmholtzError) -> Self {
        Self::Solver(e)
    }
}
impl core::fmt::Display for BakeError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "radiation bake: {self:?}")
    }
}
impl std::error::Error for BakeError {}

impl BakeDeclaration {
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        geometry_hash: ContentHash,
        band_hz: (f64, f64),
        frequencies_hz: Vec<f64>,
        operator_tol: f64,
        reference_plane: ReferencePlane,
        modes: ModeContent,
        sh_degree: usize,
        medium: Medium,
    ) -> Result<Self, BakeError> {
        if geometry_hash.as_bytes().iter().all(|b| *b == 0) {
            return Err(BakeError::Invalid(
                "source geometry hash is required, separately from mesh identity",
            ));
        }
        if !band_hz.0.is_finite()
            || !band_hz.1.is_finite()
            || band_hz.0 <= 0.0
            || band_hz.1 <= band_hz.0
            || !(2..=4096).contains(&frequencies_hz.len())
            || frequencies_hz.first() != Some(&band_hz.0)
            || frequencies_hz.last() != Some(&band_hz.1)
            || frequencies_hz.iter().any(|f| !f.is_finite())
            || !frequencies_hz.windows(2).all(|w| w[0] < w[1])
        {
            return Err(BakeError::Invalid(
                "finite positive band and 2..=4096 increasing frequencies including both endpoints are required",
            ));
        }
        if !operator_tol.is_finite()
            || operator_tol <= 0.0
            || operator_tol >= 1.0
            || sh_degree > MAX_SH_DEGREE
            || !medium.density.is_finite()
            || medium.density <= 0.0
            || !medium.sound_speed.is_finite()
            || medium.sound_speed <= 0.0
        {
            return Err(BakeError::Invalid(
                "positive relative tolerance below one, bounded SH degree and physical medium are required",
            ));
        }
        if let ReferencePlane::DuctMouth {
            origin_m,
            outward_normal,
        } = reference_plane
        {
            let norm = outward_normal.iter().map(|v| v * v).sum::<f64>();
            if origin_m
                .iter()
                .chain(&outward_normal)
                .any(|v| !v.is_finite())
                || (norm - 1.0).abs() > 1e-12
            {
                return Err(BakeError::Invalid(
                    "reference plane requires finite SI coordinates and a unit outward normal",
                ));
            }
        }
        if let ModeContent::Azimuthal { modes } = &modes {
            if modes.is_empty()
                || modes.len() > 256
                || modes[0] != (0, 0)
                || modes
                    .iter()
                    .enumerate()
                    .any(|(i, m)| modes[..i].contains(m))
            {
                return Err(BakeError::Invalid(
                    "modal content needs the plane mode first and at most 256 unique modes",
                ));
            }
        }
        let mut bytes = geometry_hash.as_bytes().to_vec();
        for value in [
            band_hz.0,
            band_hz.1,
            operator_tol,
            medium.density,
            medium.sound_speed,
        ] {
            bytes.extend_from_slice(&value.to_bits().to_le_bytes());
        }
        bytes.extend_from_slice(&(sh_degree as u64).to_le_bytes());
        bytes.extend_from_slice(&(frequencies_hz.len() as u64).to_le_bytes());
        for f in &frequencies_hz {
            bytes.extend_from_slice(&f.to_bits().to_le_bytes());
        }
        match reference_plane {
            ReferencePlane::WholeSurface => bytes.push(0),
            ReferencePlane::DuctMouth {
                origin_m,
                outward_normal,
            } => {
                bytes.push(1);
                for value in origin_m.into_iter().chain(outward_normal) {
                    bytes.extend_from_slice(&value.to_bits().to_le_bytes());
                }
            }
        }
        match &modes {
            ModeContent::Plane => bytes.push(0),
            ModeContent::Azimuthal { modes } => {
                bytes.push(1);
                bytes.extend_from_slice(&(modes.len() as u64).to_le_bytes());
                for (m, n) in modes {
                    bytes.extend_from_slice(&m.to_le_bytes());
                    bytes.extend_from_slice(&n.to_le_bytes());
                }
            }
        }
        let hash = hash_domain("org.frankensim.fs-bem.bake-declaration.v1", &bytes);
        Ok(Self {
            geometry_hash,
            band_hz,
            frequencies_hz,
            operator_tol,
            reference_plane,
            modes,
            sh_degree,
            medium,
            hash,
        })
    }
    pub const fn hash(&self) -> ContentHash {
        self.hash
    }
    pub const fn geometry_hash(&self) -> ContentHash {
        self.geometry_hash
    }
    pub const fn band_hz(&self) -> (f64, f64) {
        self.band_hz
    }
    pub fn frequencies_hz(&self) -> &[f64] {
        &self.frequencies_hz
    }
    pub const fn operator_tolerance(&self) -> f64 {
        self.operator_tol
    }
    pub const fn reference_plane(&self) -> ReferencePlane {
        self.reference_plane
    }
    pub const fn medium(&self) -> Medium {
        self.medium
    }
}

/// A rung regenerated directly from the same source geometry, never subdivided
/// from the previous rung's faceted surface. The source owns that invariant.
pub struct BakeMesh {
    pub surface: SpherePanels,
    pub driven: Vec<bool>,
}

/// Geometry-to-mesh seam. Counts must be exact BEFORE allocation, and every mesh
/// must describe this source identity and physical reference surface. This is a
/// trusted geometry-producer contract, not verification of an arbitrary producer.
pub trait PanelSource {
    fn geometry_hash(&self) -> ContentHash;
    fn reference_plane(&self) -> ReferencePlane;
    fn surface_area_m2(&self) -> f64;
    fn panel_count(&self, resolution: usize) -> Result<usize, BakeError>;
    fn mesh(&self, resolution: usize) -> Result<BakeMesh, BakeError>;
}

pub struct PanelLadder<'a> {
    source: &'a dyn PanelSource,
    resolutions: Vec<usize>,
}
impl<'a> PanelLadder<'a> {
    pub fn new(source: &'a dyn PanelSource, resolutions: Vec<usize>) -> Result<Self, BakeError> {
        if !(3..=16).contains(&resolutions.len()) || !resolutions.windows(2).all(|w| w[0] < w[1]) {
            return Err(BakeError::Invalid(
                "panel refinement needs 3..=16 increasing source resolutions",
            ));
        }
        let area = source.surface_area_m2();
        if !area.is_finite() || area <= 0.0 {
            return Err(BakeError::Invalid(
                "source surface area must be positive and finite",
            ));
        }
        let mut previous = 0;
        for &resolution in &resolutions {
            let panels = source.panel_count(resolution)?;
            if panels <= previous {
                return Err(BakeError::Invalid(
                    "panel counts must be positive and strictly increase",
                ));
            }
            previous = panels;
        }
        Ok(Self {
            source,
            resolutions,
        })
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub enum BakeStatus {
    Complete,
    NotConverged {
        at_hz: f64,
    },
    Blocked {
        at_hz: f64,
        needed_panels: usize,
        route: &'static str,
    },
    PartialBand {
        proven_hz: (f64, f64),
    },
}

#[derive(Clone, Debug)]
pub struct RefinementRow {
    pub frequency_hz: f64,
    pub rung: usize,
    pub panels: usize,
    pub source_area_m2: f64,
    /// Necessary area/phase estimate; not a sufficient accuracy or mesh count bound.
    pub minimum_panel_estimate: f64,
    pub panels_per_wavelength: Option<f64>,
    pub z_rel_change: Option<f64>,
    pub sh_rel_change: Option<f64>,
    pub status: &'static str,
    /// `"dense"` LU, or `"iterative"` (M2.5-W5) above [`MAX_DENSE_PANELS`].
    pub route: &'static str,
    pub elapsed_s: f64,
}

/// Only converged rows are retained for consumption. Failed/capped rungs remain
/// diagnostics and cannot enter the load table. The declaration is immutable.
#[derive(Clone, Debug)]
pub struct DeclaredBake {
    declaration: BakeDeclaration,
    rows: Vec<RadiationLoadRow>,
    refinement: Vec<RefinementRow>,
    status: BakeStatus,
}
impl DeclaredBake {
    pub fn declaration(&self) -> &BakeDeclaration {
        &self.declaration
    }
    pub fn rows(&self) -> &[RadiationLoadRow] {
        &self.rows
    }
    pub fn refinement(&self) -> &[RefinementRow] {
        &self.refinement
    }
    pub const fn status(&self) -> BakeStatus {
        self.status
    }
    pub fn require_declaration(&self, expected: ContentHash) -> Result<(), BakeError> {
        if expected != self.declaration.hash {
            return Err(BakeError::DeclarationMismatch);
        }
        Ok(())
    }
    /// Explicit partial-band export. A blocked bake is never silently delivered
    /// as its requested full band. At least two contiguous proven rows are needed.
    pub fn partial_band(&self, expected: ContentHash) -> Result<PartialBand<'_>, BakeError> {
        self.require_declaration(expected)?;
        if self.rows.len() < 2 {
            return Err(BakeError::OutOfBand {
                at_hz: self.declaration.band_hz.0,
                proven_hz: None,
            });
        }
        let band = (
            self.rows[0].omega / std::f64::consts::TAU,
            self.rows[self.rows.len() - 1].omega / std::f64::consts::TAU,
        );
        Ok(PartialBand {
            bake: self,
            proven_hz: band,
        })
    }
    /// Complete scalar data at the declared reference surface. Whole-surface body
    /// modes are not duct mouth loads and refuse this export.
    pub fn load_rows(&self, expected: ContentHash) -> Result<(Vec<f64>, Vec<C64>), BakeError> {
        self.require_declaration(expected)?;
        if self.status != BakeStatus::Complete {
            return Err(BakeError::Invalid("blocked bake is not a full-band load"));
        }
        if self.declaration.reference_plane == ReferencePlane::WholeSurface {
            return Err(BakeError::Invalid(
                "whole-surface radiation is not a duct reference plane",
            ));
        }
        Ok(self.scalar_rows())
    }
    fn scalar_rows(&self) -> (Vec<f64>, Vec<C64>) {
        (
            self.rows.iter().map(|r| r.omega).collect(),
            self.rows.iter().map(|r| r.z_acoustic).collect(),
        )
    }
    /// Host/command are supplied by the actual runner; library code never invents
    /// execution provenance or writes to stdout. Every row includes its declaration.
    pub fn receipt_lines(&self, host: &str, command: &str) -> Vec<String> {
        let hash = self.declaration.hash.to_hex();
        self.refinement.iter().map(|r|format!(
            "{{\"case\":\"bake-band-refinement\",\"inputs\":\"{}\",\"declaration_hash\":\"{hash}\",\"axis\":\"panels\",\"rung\":{},\"at_hz\":{:e},\"panels\":{},\"surface_area_m2\":{:e},\"minimum_panels\":{},\"ppw\":{},\"z_rel_change\":{},\"sh_rel_change\":{},\"observed\":{},\"err\":{},\"tol\":{:e},\"status\":\"{}\",\"route\":\"{}\",\"verdict\":\"{}\",\"host\":{},\"command\":{},\"elapsed_s\":{:e}}}",
            self.declaration.geometry_hash.to_hex(),r.rung,r.frequency_hz,r.panels,r.source_area_m2,
            json_number(Some(r.minimum_panel_estimate)),json_number(r.panels_per_wavelength),json_number(r.z_rel_change),
            json_number(r.sh_rel_change),json_number(r.z_rel_change),json_number(r.z_rel_change.zip(r.sh_rel_change).map(|(z,sh)|z.max(sh))),
            self.declaration.operator_tol,r.status,r.route,if r.status=="converged" {"pass"} else if r.z_rel_change.is_some() {"fail"} else {"refused"},
            json_string(host),json_string(command),r.elapsed_s)).collect()
    }
}

pub struct PartialBand<'a> {
    bake: &'a DeclaredBake,
    proven_hz: (f64, f64),
}
impl PartialBand<'_> {
    pub const fn status(&self) -> BakeStatus {
        BakeStatus::PartialBand {
            proven_hz: self.proven_hz,
        }
    }
    pub const fn band_refusal(&self) -> (f64, f64) {
        self.proven_hz
    }
    pub fn load_rows(&self) -> Result<(Vec<f64>, Vec<C64>), BakeError> {
        if self.bake.declaration.reference_plane == ReferencePlane::WholeSurface {
            return Err(BakeError::Invalid(
                "whole-surface radiation is not a duct reference plane",
            ));
        }
        Ok(self.bake.scalar_rows())
    }
}

fn json_number(value: Option<f64>) -> String {
    value
        .filter(|v| v.is_finite())
        .map_or_else(|| "null".into(), |v| format!("{v:e}"))
}
fn json_string(value: &str) -> String {
    let mut s = String::from("\"");
    for c in value.chars() {
        match c {
            '"' => s.push_str("\\\""),
            '\\' => s.push_str("\\\\"),
            '\n' => s.push_str("\\n"),
            '\r' => s.push_str("\\r"),
            '\t' => s.push_str("\\t"),
            c if c.is_control() => {
                use std::fmt::Write;
                let _ = write!(s, "\\u{:04x}", c as u32);
            }
            c => s.push(c),
        }
    }
    s.push('"');
    s
}

fn relative_change(a: &[C64], b: &[C64]) -> f64 {
    let scale = a.iter().chain(b).map(|z| z.abs()).fold(0.0_f64, f64::max);
    if scale == 0.0 {
        return 0.0;
    }
    let error = a
        .iter()
        .zip(b)
        .map(|(&x, &y)| ((x - y).abs() / scale).powi(2))
        .sum::<f64>();
    let norm = a
        .iter()
        .map(|z| (z.abs() / scale).powi(2))
        .sum::<f64>()
        .max(b.iter().map(|z| (z.abs() / scale).powi(2)).sum::<f64>());
    det::sqrt(error / norm)
}

fn dense_status(at_hz: f64, panels: usize) -> Option<BakeStatus> {
    (panels > MAX_DENSE_PANELS).then_some(BakeStatus::Blocked {
        at_hz,
        needed_panels: panels,
        route: "M2.5-W5",
    })
}

fn admit_reference_plane(mesh: &BakeMesh, plane: ReferencePlane) -> Result<(), BakeError> {
    if mesh.driven.len() != mesh.surface.areas().len() {
        return Err(BakeError::Invalid(
            "driven mask must match the generated panel surface",
        ));
    }
    if let ReferencePlane::DuctMouth {
        origin_m,
        outward_normal,
    } = plane
    {
        for (i, &driven) in mesh.driven.iter().enumerate() {
            if !driven {
                continue;
            }
            let point = mesh.surface.centroids()[i];
            let normal = mesh.surface.normals()[i];
            let offset = (0..3)
                .map(|j| (point[j] - origin_m[j]) * outward_normal[j])
                .sum::<f64>();
            let alignment = (0..3).map(|j| normal[j] * outward_normal[j]).sum::<f64>();
            let scale = point
                .iter()
                .chain(&origin_m)
                .map(|v| v.abs())
                .fold(1.0_f64, f64::max);
            if offset.abs() > 1e-12 * scale || alignment < 1.0 - 1e-12 {
                return Err(BakeError::Invalid(
                    "driven panels are not on the declared outward mouth reference plane; an explicit transport owner is required",
                ));
            }
        }
    }
    Ok(())
}

/// Execute the predeclared frequency grid, regenerating every rung from source.
/// The unchanged owner selects CBIE/Burton–Miller and computes all SH coefficients.
/// Polls cancellation between bounded existing solves; no mid-LU latency claim.
/// Dense only: a rung above [`MAX_DENSE_PANELS`] is `Blocked` and routed to
/// M2.5-W5 ([`bake_declared_routed`]).
pub fn bake_declared(
    declaration: &BakeDeclaration,
    ladder: &PanelLadder<'_>,
    cancelled: impl FnMut() -> bool,
) -> Result<DeclaredBake, BakeError> {
    bake_declared_routed(declaration, ladder, None, cancelled)
}

/// [`bake_declared`] with the M2.5-W5 route. With `iterative`, a rung above
/// [`MAX_DENSE_PANELS`] is solved by the matrix-free iterative owner under
/// those explicit tolerance, FMM-order and storage budgets instead of being
/// `Blocked`. Rungs within the dense cap stay dense. The same ladder,
/// operator tolerance and three-rung rule decide completeness, and every
/// refinement row records its route. `cancelled` is polled between rungs and
/// frequencies; one iterative rung runs to completion once started.
pub fn bake_declared_routed(
    declaration: &BakeDeclaration,
    ladder: &PanelLadder<'_>,
    iterative: Option<&IterativeOptions>,
    mut cancelled: impl FnMut() -> bool,
) -> Result<DeclaredBake, BakeError> {
    if declaration.geometry_hash != ladder.source.geometry_hash()
        || declaration.reference_plane != ladder.source.reference_plane()
    {
        return Err(BakeError::DeclarationMismatch);
    }
    if !matches!(declaration.modes, ModeContent::Plane) {
        return Err(BakeError::Dependency { route: "M2.1-W3c" });
    }
    let mut bake = DeclaredBake {
        declaration: declaration.clone(),
        rows: Vec::new(),
        refinement: Vec::new(),
        status: BakeStatus::Complete,
    };
    if cancelled() {
        return Err(BakeError::Cancelled);
    }
    // Three solved rungs are mandatory. If even the third source mesh exceeds
    // the dense budget, the declared ladder cannot qualify at ANY frequency;
    // do not spend two dense factorizations rediscovering that fact.
    let preflight_start = Instant::now();
    for (rung, &resolution) in ladder.resolutions.iter().take(3).enumerate() {
        let panels = ladder.source.panel_count(resolution)?;
        if iterative.is_none()
            && let Some(status) = dense_status(declaration.band_hz.0, panels)
        {
            let frequency = declaration.band_hz.0;
            let area = ladder.source.surface_area_m2();
            let h = declaration.medium.sound_speed / (frequency * MIN_PANELS_PER_WAVELENGTH);
            bake.status = status;
            bake.refinement.push(RefinementRow {
                frequency_hz: frequency,
                rung,
                panels,
                source_area_m2: area,
                minimum_panel_estimate: (area / h / h).ceil(),
                panels_per_wavelength: None,
                z_rel_change: None,
                sh_rel_change: None,
                status: "blocked",
                route: "dense",
                elapsed_s: preflight_start.elapsed().as_secs_f64(),
            });
            return Ok(bake);
        }
    }
    for &frequency in &declaration.frequencies_hz {
        let area = ladder.source.surface_area_m2();
        let h = declaration.medium.sound_speed / (frequency * MIN_PANELS_PER_WAVELENGTH);
        let estimate = (area / h / h).ceil();
        let mut previous: Option<RadiationLoadRow> = None;
        let mut last_good = false;
        let mut solved = 0;
        let mut last_z = None;
        let mut last_sh = None;
        let mut converged = false;
        for (rung, &resolution) in ladder.resolutions.iter().enumerate() {
            if cancelled() {
                return Err(BakeError::Cancelled);
            }
            let start = Instant::now();
            let panels = ladder.source.panel_count(resolution)?;
            let routed = iterative.filter(|_| panels > MAX_DENSE_PANELS);
            let mut record = RefinementRow {
                frequency_hz: frequency,
                rung,
                panels,
                source_area_m2: area,
                minimum_panel_estimate: estimate,
                panels_per_wavelength: None,
                z_rel_change: None,
                sh_rel_change: None,
                status: "refining",
                route: if routed.is_some() {
                    "iterative"
                } else {
                    "dense"
                },
                elapsed_s: 0.0,
            };
            if iterative.is_none()
                && let Some(status) = dense_status(frequency, panels)
            {
                record.status = "blocked";
                record.elapsed_s = start.elapsed().as_secs_f64();
                bake.refinement.push(record);
                bake.status = status;
                return Ok(bake);
            }
            let mesh = ladder.source.mesh(resolution)?;
            if mesh.surface.areas().len() != panels {
                return Err(BakeError::Invalid(
                    "source panel-count preflight disagrees with generated geometry",
                ));
            }
            admit_reference_plane(&mesh, declaration.reference_plane)?;
            let omega = [std::f64::consts::TAU * frequency];
            let source_id = declaration.hash.to_hex();
            let (medium, degree) = (declaration.medium, declaration.sh_degree);
            let result = match routed {
                Some(options) => bake_radiation_load_iterative(
                    &mesh.surface,
                    &mesh.driven,
                    &omega,
                    medium,
                    degree,
                    &source_id,
                    options,
                )
                .map_err(|e| match e {
                    IterativeError::Helmholtz(e) => BakeError::Solver(e),
                    other => BakeError::IterativeRoute {
                        at_hz: frequency,
                        what: other.to_string(),
                    },
                }),
                None => bake_radiation_load(
                    &mesh.surface,
                    &mesh.driven,
                    &omega,
                    medium,
                    degree,
                    &source_id,
                )
                .map_err(BakeError::Solver),
            };
            let mut row = match result {
                Ok(value) => value.rows,
                Err(BakeError::Solver(HelmholtzError::TooCoarse {
                    panels_per_wavelength,
                })) => {
                    record.status = "under-resolved";
                    record.panels_per_wavelength = Some(panels_per_wavelength);
                    record.elapsed_s = start.elapsed().as_secs_f64();
                    bake.refinement.push(record);
                    continue;
                }
                Err(e) => return Err(e),
            };
            let row = row.remove(0);
            if !row.z_acoustic.re.is_finite()
                || !row.z_acoustic.im.is_finite()
                || row
                    .directivity
                    .coefficients
                    .iter()
                    .any(|z| !z.re.is_finite() || !z.im.is_finite())
            {
                return Err(BakeError::Invalid("nonfinite computed radiation operator"));
            }
            if row.passivity_margin < -declaration.operator_tol * row.z_specific.abs() {
                return Err(BakeError::NonPassive {
                    at_hz: frequency,
                    margin: row.passivity_margin,
                });
            }
            record.panels_per_wavelength = Some(row.panels_per_wavelength);
            if let Some(old) = &previous {
                last_z = Some(relative_change(&[row.z_acoustic], &[old.z_acoustic]));
                last_sh = Some(relative_change(
                    &row.directivity.coefficients,
                    &old.directivity.coefficients,
                ));
                record.z_rel_change = last_z;
                record.sh_rel_change = last_sh;
            }
            solved += 1;
            let good = last_z.is_some_and(|e| e <= declaration.operator_tol)
                && last_sh.is_some_and(|e| e <= declaration.operator_tol);
            converged = solved >= 3 && good && last_good;
            record.status = if converged { "converged" } else { "refining" };
            record.elapsed_s = start.elapsed().as_secs_f64();
            bake.refinement.push(record);
            if converged {
                bake.rows.push(row);
                break;
            }
            previous = Some(row);
            last_good = good;
        }
        if !converged {
            bake.status = BakeStatus::NotConverged { at_hz: frequency };
            return Ok(bake);
        }
    }
    Ok(bake)
}

#[cfg(test)]
mod tests;
