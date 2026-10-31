//! Multi-aperture radiation impedance (M2.5-W3a): the self and mutual loading
//! of declared apertures, such as tone holes and a bell mouth, on one body.
//!
//! An aperture is a real velocity shape over a closed body's panels: outward
//! normal velocity per unit generalized aperture velocity. Under `exp(-i w t)`,
//! `Z[a][b] = SUM_i A_i shape_a(i) p_b(i)` is the generalized force on aperture
//! `a` per unit generalized velocity of aperture `b` [N s/m]. The body force is
//! `-Z v`, and a mass load has a negative imaginary part. Every row is
//! projected from the same surface-pressure owners as single radiation solves.
//! Nothing is symmetrized, clipped or diagonally loaded. Each frequency row
//! reports its reciprocity defect and the smallest eigenvalue of its Hermitian
//! part, so a consumer sees how far the discretization is from the exact
//! reciprocal, passive operator.
//!
//! Routes:
//! - [`ApertureRoute::Collocation`]: centroid Burton-Miller. Dense LU applies
//!   up to [`MAX_DENSE_PANELS`]; above it, the matrix-free FMM/GMRES owner
//!   applies when options are supplied, and otherwise the dense cap refuses.
//!   Reciprocity holds only to discretization error.
//! - [`ApertureRoute::Galerkin`]: the symmetric P1/P0 velocity-port image of
//!   `reciprocal_radiation`. It is reciprocal and passive by construction,
//!   within its restricted band.
//! - [`baffled_piston_array_impedance`]: coplanar circular pistons in an
//!   infinite rigid baffle. Self terms come from the closed-form owner; mutual
//!   terms are the Rayleigh double integral on a spectral disc-quadrature
//!   ladder.
use super::iterative::{IterativeError, IterativeOptions, IterativeRadiation};
use super::{
    Formulation, HelmholtzError, MAX_DENSE_PANELS, MAX_RADIATION_FIELDS_PER_BATCH, Medium,
    baffled_piston_face_impedance, gauss_legendre, green, solve_radiation_batch,
    surface_fingerprint,
};
use crate::panel3d::SpherePanels;
use crate::reciprocal_radiation::{Budget, Prepared};
use fs_exec::CancelGate;
use fs_math::{c64::C64, det};

/// One declared aperture on a closed radiating body.
#[derive(Clone, Debug)]
pub struct Aperture {
    /// Outward normal velocity per unit generalized aperture velocity, one
    /// real coefficient per panel, in the surface's panel order.
    pub shape: Vec<f64>,
    /// A point on the aperture's declared reference plane [m].
    pub reference_point: [f64; 3],
    /// Outward unit normal of that reference plane.
    pub normal: [f64; 3],
}

/// Which surface-pressure owner produces the aperture pressures.
#[derive(Clone, Debug)]
pub enum ApertureRoute {
    /// Centroid Burton-Miller collocation. `None` keeps the dense cap as a
    /// refusal; `Some` routes bodies above it to the iterative owner.
    Collocation { iterative: Option<IterativeOptions> },
    /// Symmetric P1/P0 Galerkin velocity ports, admitted for the highest
    /// requested wavenumber.
    Galerkin { budget: Budget },
}

/// One frequency row of the aperture matrix.
#[derive(Clone, Debug)]
pub struct ApertureRow {
    /// Angular frequency [rad/s].
    pub omega: f64,
    /// Row-major `apertures x apertures` impedance [N s/m].
    pub impedance: Vec<C64>,
    /// `"dense"`, `"iterative"`, `"galerkin"` or `"molecular-dense"`.
    pub route: &'static str,
    /// `||Z - Z^T||_F / ||Z||_F`, measured and never repaired.
    pub reciprocity_defect: f64,
    /// Smallest eigenvalue of `(Z + Z^H)/2` over `||Z||_F`. Mean radiated
    /// power is `v^H (Z + Z^H) v / 4`, so a negative value is a measured
    /// discrete power deficit, not something to clip.
    pub passivity_margin: f64,
}

/// A bound aperture matrix artifact over a frequency band.
#[derive(Clone, Debug)]
pub struct ApertureImpedance {
    /// Hash of the surface, aperture shapes and reference planes, medium,
    /// band, route, and the route's accuracy parameters (FMM order, GMRES
    /// tolerance, angular refinement): a cross-wiring guard, not a certificate.
    pub declaration: u64,
    /// Number of apertures.
    pub apertures: usize,
    /// Outward volume flow per unit generalized velocity of each aperture,
    /// `q_a = SUM_i A_i shape_a(i)` [m^3/s per m/s].
    pub volume_gains: Vec<f64>,
    /// One row per requested frequency, in input order.
    pub rows: Vec<ApertureRow>,
}

impl ApertureImpedance {
    /// Row `index` as an acoustic matrix `Z[a][b] / (q_a q_b)` [Pa s/m^3]:
    /// mean aperture pressure per outward volume flow, the plane-port
    /// impedance a duct network terminates its openings with
    /// (`fs_duct::network::CoupledTermination`).
    ///
    /// # Errors
    /// An out-of-range row, or a zero volume gain (a shape with no net flow
    /// has no plane-port meaning).
    pub fn acoustic(&self, index: usize) -> Result<Vec<C64>, HelmholtzError> {
        let row = self
            .rows
            .get(index)
            .ok_or(bad("aperture row index out of range"))?;
        if self.volume_gains.iter().any(|q| !(q.abs() > 0.0)) {
            return Err(bad(
                "an aperture with zero net volume flow has no plane-port impedance",
            ));
        }
        let n = self.apertures;
        Ok((0..n * n)
            .map(|i| {
                row.impedance[i].scale(1.0 / (self.volume_gains[i / n] * self.volume_gains[i % n]))
            })
            .collect())
    }
}

/// Refusals from either pressure owner.
#[derive(Clone, Debug)]
pub enum ApertureError {
    /// Validation, dense or Galerkin refusal.
    Helmholtz(HelmholtzError),
    /// Iterative-owner refusal (including non-convergence).
    Iterative(IterativeError),
}
impl core::fmt::Display for ApertureError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Helmholtz(e) => e.fmt(f),
            Self::Iterative(e) => e.fmt(f),
        }
    }
}
impl std::error::Error for ApertureError {}
impl From<HelmholtzError> for ApertureError {
    fn from(e: HelmholtzError) -> Self {
        Self::Helmholtz(e)
    }
}
impl From<IterativeError> for ApertureError {
    fn from(e: IterativeError) -> Self {
        match e {
            IterativeError::Helmholtz(h) => Self::Helmholtz(h),
            other => Self::Iterative(other),
        }
    }
}

fn bad(what: &'static str) -> HelmholtzError {
    HelmholtzError::BadParameter { what }
}

fn finite(z: C64) -> bool {
    z.re.is_finite() && z.im.is_finite()
}

pub(super) fn mix(hash: u64, word: u64) -> u64 {
    (hash ^ word).wrapping_mul(0x0000_0100_0000_01b3)
}

/// Build the aperture matrix of `surface` at every `omegas` entry.
///
/// # Errors
/// Invalid apertures (wrong length, nonfinite or all-zero shapes, non-unit
/// reference normals), an empty or non-increasing band, any refusal of the
/// selected pressure owner, or cooperative cancellation.
pub fn aperture_impedance(
    surface: &SpherePanels,
    apertures: &[Aperture],
    omegas: &[f64],
    medium: Medium,
    route: &ApertureRoute,
    gate: &CancelGate,
) -> Result<ApertureImpedance, ApertureError> {
    let panels = surface.areas().len();
    let n = apertures.len();
    admit_apertures(surface, apertures)?;
    if omegas.is_empty()
        || omegas.iter().any(|w| !(w.is_finite() && *w > 0.0))
        || omegas.windows(2).any(|w| w[1] <= w[0])
    {
        return Err(bad("aperture band must be finite, positive and strictly increasing").into());
    }
    if !(medium.density.is_finite()
        && medium.density > 0.0
        && medium.sound_speed.is_finite()
        && medium.sound_speed > 0.0)
    {
        return Err(bad("medium density and sound speed must be positive and finite").into());
    }
    let mut declaration = aperture_declaration(surface, apertures);
    for &v in omegas.iter().chain([&medium.density, &medium.sound_speed]) {
        declaration = mix(declaration, v.to_bits());
    }
    let fields: Vec<Vec<C64>> = apertures
        .iter()
        .map(|a| a.shape.iter().map(|&v| C64::from_re(v)).collect())
        .collect();
    let mut rows = Vec::with_capacity(omegas.len());
    match route {
        ApertureRoute::Collocation { iterative } => {
            declaration = mix(declaration, 1);
            if let Some(o) = iterative.as_ref().filter(|_| panels > MAX_DENSE_PANELS) {
                for word in [
                    o.fmm_order as u64,
                    o.tol.to_bits(),
                    u64::from(o.directional.is_some()),
                ] {
                    declaration = mix(declaration, word);
                }
            }
            for &omega in omegas {
                if gate.is_requested() {
                    return Err(HelmholtzError::Cancelled.into());
                }
                let k = omega / medium.sound_speed;
                let (pressures, label) = match iterative {
                    Some(options) if panels > MAX_DENSE_PANELS => {
                        let operator = IterativeRadiation::prepare(
                            surface,
                            k,
                            medium,
                            Formulation::BurtonMiller,
                            options,
                            gate,
                        )?;
                        let mut out = Vec::with_capacity(n);
                        for field in &fields {
                            let mut state = operator.start(field, gate)?;
                            operator.advance(
                                &mut state,
                                options.max_iters.saturating_add(1),
                                gate,
                            )?;
                            out.push(operator.finish(&state)?.0.pressure);
                        }
                        (out, "iterative")
                    }
                    _ => {
                        let refs: Vec<&[C64]> = fields.iter().map(Vec::as_slice).collect();
                        let solutions = solve_radiation_batch(
                            surface,
                            k,
                            medium,
                            &refs,
                            Formulation::BurtonMiller,
                        )?;
                        (solutions.into_iter().map(|s| s.pressure).collect(), "dense")
                    }
                };
                let pressure_refs: Vec<_> = pressures.iter().map(Vec::as_slice).collect();
                let impedance = project_pressures(surface, apertures, &pressure_refs)?;
                rows.push(row(omega, impedance, n, label));
            }
        }
        ApertureRoute::Galerkin { budget } => {
            declaration = mix(declaration, 2);
            declaration = mix(declaration, u64::from(budget.angular_refinement));
            let highest_k = omegas[omegas.len() - 1] / medium.sound_speed;
            let prepared = Prepared::new(surface, highest_k, *budget, || gate.is_requested())?;
            let shapes: Vec<&[f64]> = apertures.iter().map(|a| a.shape.as_slice()).collect();
            for &omega in omegas {
                let ports = prepared.solve_velocity_ports(
                    omega / medium.sound_speed,
                    medium,
                    &shapes,
                    || gate.is_requested(),
                )?;
                rows.push(row(omega, ports.impedance().to_vec(), n, "galerkin"));
            }
        }
    }
    let volume_gains = apertures
        .iter()
        .map(|a| {
            a.shape
                .iter()
                .zip(surface.areas())
                .map(|(s, area)| s * area)
                .sum()
        })
        .collect();
    Ok(ApertureImpedance {
        declaration,
        apertures: n,
        volume_gains,
        rows,
    })
}

pub(super) fn admit_apertures(
    surface: &SpherePanels,
    apertures: &[Aperture],
) -> Result<(), HelmholtzError> {
    let panels = surface.areas().len();
    if apertures.is_empty() || apertures.len() > MAX_RADIATION_FIELDS_PER_BATCH {
        return Err(bad(
            "aperture set must be nonempty and within MAX_RADIATION_FIELDS_PER_BATCH",
        ));
    }
    for aperture in apertures {
        if aperture.shape.len() != panels {
            return Err(HelmholtzError::ShapeMismatch {
                what: "aperture shape length must equal the panel count",
            });
        }
        if aperture.shape.iter().any(|v| !v.is_finite()) || aperture.shape.iter().all(|&v| v == 0.0)
        {
            return Err(bad(
                "aperture shapes must be finite and not identically zero",
            ));
        }
        let m = aperture.normal;
        let norm2 = m[0] * m[0] + m[1] * m[1] + m[2] * m[2];
        if aperture.reference_point.iter().any(|v| !v.is_finite())
            || !norm2.is_finite()
            || (norm2 - 1.0).abs() > 1.0e-12
        {
            return Err(bad(
                "aperture reference planes need a finite point and a unit normal",
            ));
        }
    }
    Ok(())
}

pub(super) fn aperture_declaration(surface: &SpherePanels, apertures: &[Aperture]) -> u64 {
    let mut declaration = surface_fingerprint(surface);
    for aperture in apertures {
        for &v in aperture
            .shape
            .iter()
            .chain(&aperture.reference_point)
            .chain(&aperture.normal)
        {
            declaration = mix(declaration, v.to_bits());
        }
    }
    declaration
}

pub(super) fn project_pressures(
    surface: &SpherePanels,
    apertures: &[Aperture],
    pressures: &[&[C64]],
) -> Result<Vec<C64>, HelmholtzError> {
    project_real_velocity_rows(
        surface,
        apertures.len(),
        |row, panel| apertures[row].shape[panel],
        pressures,
        None,
    )
}

/// Shared work-conjugate projection for original real velocity shapes. The
/// accessor avoids allocating another copy of a retained complex-real basis.
/// Callers own matrix/projection budgets; this owner retains the exact sum and
/// polls the optional gate between bounded panel tiles.
pub(crate) fn project_real_velocity_rows(
    surface: &SpherePanels,
    n: usize,
    shape: impl Fn(usize, usize) -> f64,
    pressures: &[&[C64]],
    gate: Option<&CancelGate>,
) -> Result<Vec<C64>, HelmholtzError> {
    if pressures.len() != n || pressures.iter().any(|p| p.len() != surface.areas().len()) {
        return Err(bad("source-pressure projection dimensions disagree"));
    }
    let mut impedance = vec![C64::ZERO; n * n];
    for (b, pressure) in pressures.iter().enumerate() {
        for a in 0..n {
            let mut force = C64::ZERO;
            for (panel, (&p, &area)) in pressure.iter().zip(surface.areas()).enumerate() {
                if panel % 256 == 0 && gate.is_some_and(CancelGate::is_requested) {
                    return Err(HelmholtzError::Cancelled);
                }
                force = force + p.scale(area * shape(a, panel));
            }
            if !finite(force) {
                return Err(bad("nonfinite aperture force projection"));
            }
            impedance[a * n + b] = force;
        }
    }
    Ok(impedance)
}

pub(super) fn row(omega: f64, impedance: Vec<C64>, n: usize, route: &'static str) -> ApertureRow {
    let (reciprocity_defect, passivity_margin) = reciprocity_and_passivity(&impedance, n);
    ApertureRow {
        omega,
        impedance,
        route,
        reciprocity_defect,
        passivity_margin,
    }
}

/// Measured reciprocity defect and normalized smallest Hermitian-part
/// eigenvalue. The Hermitian part is diagonalized through its real symmetric
/// embedding `[[Re H, -Im H], [Im H, Re H]]`, whose spectrum is that of `H`
/// with every eigenvalue doubled.
#[must_use]
pub fn reciprocity_and_passivity(impedance: &[C64], n: usize) -> (f64, f64) {
    let norm = det::sqrt(impedance.iter().map(|z| z.norm_sq()).sum::<f64>());
    let mut embedding = vec![0.0; 4 * n * n];
    for a in 0..n {
        for b in 0..n {
            let h = (impedance[a * n + b] + impedance[b * n + a].conj()).scale(0.5);
            embedding[a * 2 * n + b] = h.re;
            embedding[(a + n) * 2 * n + b + n] = h.re;
            embedding[a * 2 * n + b + n] = -h.im;
            embedding[(a + n) * 2 * n + b] = h.im;
        }
    }
    let (values, _) = fs_la::eigen::jacobi_eigh(&embedding, 2 * n);
    (reciprocity_defect(impedance, n), values[0] / norm)
}

/// The same original transpose defect without allocating an eigenproblem.
pub(crate) fn reciprocity_defect(impedance: &[C64], n: usize) -> f64 {
    let norm = det::sqrt(impedance.iter().map(|z| z.norm_sq()).sum::<f64>());
    let mut asymmetry = 0.0;
    for a in 0..n {
        for b in 0..n {
            asymmetry += (impedance[a * n + b] - impedance[b * n + a]).norm_sq();
        }
    }
    det::sqrt(asymmetry) / norm
}

/// One coplanar circular piston in an infinite rigid baffle.
#[derive(Clone, Copy, Debug)]
pub struct BaffledPiston {
    /// Centre in the baffle plane [m].
    pub center_m: [f64; 2],
    /// Radius [m].
    pub radius_m: f64,
}

/// Disc-quadrature orders `(radial, azimuthal)` of the mutual-term ladder.
const PISTON_LADDER: [(usize, usize); 4] = [(8, 16), (16, 32), (32, 64), (32, 128)];
/// Relative change between successive rungs that admits a mutual term.
/// Convergence is exponential, so the kept finer value is far more accurate.
pub const PISTON_MUTUAL_TOLERANCE: f64 = 1.0e-9;

fn disc_rule(piston: &BaffledPiston, radial: usize, azimuthal: usize) -> Vec<([f64; 2], f64)> {
    let (nodes, weights) = gauss_legendre(radial);
    let a = piston.radius_m;
    let mut rule = Vec::with_capacity(radial * azimuthal);
    for (&x, &w) in nodes.iter().zip(&weights) {
        let r = 0.5 * a * (x + 1.0);
        for j in 0..azimuthal {
            let theta = core::f64::consts::TAU * (j as f64 + 0.5) / azimuthal as f64;
            let (s, c) = det::sin_cos(theta);
            rule.push((
                [piston.center_m[0] + r * c, piston.center_m[1] + r * s],
                0.5 * a * w * r * core::f64::consts::TAU / azimuthal as f64,
            ));
        }
    }
    rule
}

/// Force/velocity impedance matrix [N s/m] of coplanar baffled pistons.
/// Diagonal entries are the closed-form face impedance times the piston area.
/// Off-diagonal entries integrate `-2 i w rho G(|x - y|)` over both discs with
/// tensor Gauss-Legendre (radius) times midpoint (azimuth) rules. The ladder
/// stops at the first two successive rungs that agree within
/// [`PISTON_MUTUAL_TOLERANCE`] and keeps the finer value. Both triangles of
/// the matrix are evaluated separately. The azimuthal rule converges like
/// `exp(-m acosh((a^2 + s^2)/(2 a s)))` for a disc of radius `a` seen from
/// distance `s`, so small gaps need the finest rungs.
///
/// # Errors
/// Invalid pistons, overlapping discs, a medium or wavenumber outside the
/// closed-form owner's domain, or a mutual quadrature that did not converge
/// (nearly touching discs).
pub fn baffled_piston_array_impedance(
    pistons: &[BaffledPiston],
    k: f64,
    medium: Medium,
) -> Result<Vec<C64>, HelmholtzError> {
    let n = pistons.len();
    if n == 0 || n > MAX_RADIATION_FIELDS_PER_BATCH {
        return Err(bad(
            "piston array must be nonempty and within MAX_RADIATION_FIELDS_PER_BATCH",
        ));
    }
    if pistons.iter().any(|p| {
        !(p.radius_m.is_finite() && p.radius_m > 0.0) || p.center_m.iter().any(|v| !v.is_finite())
    }) {
        return Err(bad("pistons need finite centres and positive finite radii"));
    }
    let omega_rho = k * medium.sound_speed * medium.density;
    let mut z = vec![C64::ZERO; n * n];
    for a in 0..n {
        let pa = &pistons[a];
        z[a * n + a] = baffled_piston_face_impedance(pa.radius_m, k, medium)?
            .scale(core::f64::consts::PI * pa.radius_m * pa.radius_m);
        for b in 0..n {
            if a == b {
                continue;
            }
            let pb = &pistons[b];
            let d = det::hypot(
                pa.center_m[0] - pb.center_m[0],
                pa.center_m[1] - pb.center_m[1],
            );
            if !(d > pa.radius_m + pb.radius_m) {
                return Err(bad("baffled pistons must not overlap or touch"));
            }
            let mut previous: Option<C64> = None;
            let mut accepted = None;
            for (radial, azimuthal) in PISTON_LADDER {
                if accepted.is_some() {
                    break;
                }
                let (ra, rb) = (
                    disc_rule(pa, radial, azimuthal),
                    disc_rule(pb, radial, azimuthal),
                );
                let mut sum = C64::ZERO;
                for &(x, wx) in &ra {
                    for &(y, wy) in &rb {
                        let r = det::hypot(x[0] - y[0], x[1] - y[1]);
                        sum = sum + green(k, r).scale(wx * wy);
                    }
                }
                let value = sum * C64::new(0.0, -2.0 * omega_rho);
                if !finite(value) {
                    return Err(bad("nonfinite baffled-piston mutual impedance"));
                }
                if let Some(p) = previous {
                    accepted = ((value - p).abs() <= PISTON_MUTUAL_TOLERANCE * value.abs())
                        .then_some(value);
                }
                previous = Some(value);
            }
            z[a * n + b] = accepted.ok_or(bad(
                "baffled-piston mutual quadrature did not converge (discs nearly touching)",
            ))?;
        }
    }
    Ok(z)
}

#[cfg(test)]
#[path = "aperture_impedance_tests.rs"]
mod tests;
