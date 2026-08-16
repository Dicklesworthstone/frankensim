//! Offline radiation-load bake: a frequency sweep over the Helmholtz
//! radiation solver producing an area-averaged `Z_L(omega)` table plus
//! per-frequency directivity for a driven mouth on an exterior mesh
//! (music bead `frankensim-music-t-brass-zl-mm-zolja`).
//!
//! The bake is OFFLINE plumbing, not new physics: each row is one
//! admitted [`crate::helmholtz::solve_radiation`] call (Burton–Miller
//! above `ka = 0.5` at the mouth, plain CBIE below it — the recorded
//! low-`ka` resistance-artifact boundary), area-averaged over the driven
//! panels exactly like the pulsating-sphere oracle's arbiter. Every row
//! retains the guardrail diagnostics (panels-per-wavelength, condition
//! lower bound, radiated power, passivity margin) so a consumer can see
//! WHY a table is trustworthy — and the sweep REFUSES past the mesh's
//! panels-per-wavelength bound instead of extrapolating.
//!
//! Units: `z_specific` is specific acoustic impedance [Pa s/m]
//! (rho*c scale); `z_acoustic` divides by the driven area
//! [Pa s/m^3] — the quantity a duct termination consumes. Same
//! `e^{-i omega t}` convention as fs-duct: mass-like reactance is
//! NEGATIVE imaginary, and no conjugation happens on the hop.

use fs_math::c64::C64;

use crate::helmholtz::{
    DirectivityTable, Formulation, HelmholtzError, Medium, directivity_sh_table, solve_radiation,
};
use crate::panel3d::SpherePanels;

/// One frequency row of a radiation-load bake.
#[derive(Debug, Clone)]
pub struct RadiationLoadRow {
    /// Angular frequency [rad/s].
    pub omega: f64,
    /// Wavenumber [1/m].
    pub k: f64,
    /// `ka` at the equivalent mouth radius `sqrt(S/pi)`.
    pub mouth_ka: f64,
    /// Area-averaged specific impedance over the driven panels [Pa s/m].
    pub z_specific: C64,
    /// `z_specific / S_mouth` [Pa s/m^3] — the duct-side load.
    pub z_acoustic: C64,
    /// Formulation actually used for this row.
    pub formulation: Formulation,
    /// Mesh resolution at this frequency (guardrail diagnostic).
    pub panels_per_wavelength: f64,
    /// Probe-based condition LOWER bound (a large value is a warning; a
    /// small value certifies nothing).
    pub condition_lower_bound: f64,
    /// Radiated power for the unit driven field [W].
    pub radiated_power: f64,
    /// Passivity margin `Re z_specific` [Pa s/m] (>= 0 for a passive row).
    pub passivity_margin: f64,
    /// Spherical-harmonic directivity of the radiated field.
    pub directivity: DirectivityTable,
}

/// A complete baked radiation load: the `Z_L(omega)` table with its
/// receipts. The table is the source of truth; consumers interpolate
/// and refuse out-of-table queries on their own side.
#[derive(Debug, Clone)]
pub struct RadiationLoadBake {
    /// Caller-supplied immutable source identity (mesh + profile).
    pub source_id: String,
    /// Driven (mouth) area [m^2].
    pub mouth_area_m2: f64,
    /// Total panel count of the exterior mesh.
    pub panel_count: usize,
    /// Number of driven panels.
    pub driven_count: usize,
    /// Deterministic surface fingerprint from the solver.
    pub surface_fingerprint: u64,
    /// Medium the bake ran in.
    pub medium: Medium,
    /// Frequency rows, strictly increasing in omega.
    pub rows: Vec<RadiationLoadRow>,
}

impl RadiationLoadBake {
    /// Smallest passivity margin across the table [Pa s/m].
    #[must_use]
    pub fn worst_passivity_margin(&self) -> f64 {
        self.rows
            .iter()
            .map(|r| r.passivity_margin)
            .fold(f64::INFINITY, f64::min)
    }

    /// Whether every row is passive (`Re z >= -tol`, `tol` relative to
    /// the row magnitude).
    #[must_use]
    pub fn is_passive(&self, relative_tolerance: f64) -> bool {
        self.rows
            .iter()
            .all(|r| r.passivity_margin >= -relative_tolerance * r.z_specific.abs())
    }

    /// JSON-lines receipt rows (frequency, panel counts, wavelength
    /// check, passivity margin) — the archived bake evidence.
    #[must_use]
    pub fn receipt_lines(&self) -> Vec<String> {
        self.rows
            .iter()
            .map(|r| {
                format!(
                    "{{\"suite\":\"fs-bem\",\"case\":\"radiation-load-bake\",\
                     \"source_id\":\"{}\",\"omega\":{:.6e},\"mouth_ka\":{:.4},\
                     \"panels\":{},\"panels_per_wavelength\":{:.2},\
                     \"formulation\":\"{}\",\"z_acoustic_re\":{:.6e},\
                     \"z_acoustic_im\":{:.6e},\"passivity_margin\":{:.3e},\
                     \"condition_lower_bound\":{:.3e},\"captured_fraction\":{:.4}}}",
                    self.source_id,
                    r.omega,
                    r.mouth_ka,
                    self.panel_count,
                    r.panels_per_wavelength,
                    match r.formulation {
                        Formulation::PlainCbie => "plain-cbie",
                        Formulation::BurtonMiller => "burton-miller",
                        Formulation::BurtonMillerWrongAlphaSign => "wrong-alpha",
                    },
                    r.z_acoustic.re,
                    r.z_acoustic.im,
                    r.passivity_margin,
                    r.condition_lower_bound,
                    r.directivity.captured_fraction,
                )
            })
            .collect()
    }
}

/// Bake `Z_L(omega)` for the driven mouth of an exterior mesh.
///
/// `driven[i]` marks the panels of the radiating mouth (unit outward
/// normal velocity there, rigid elsewhere). The formulation switches at
/// `ka = 0.5` on the equivalent mouth radius (plain CBIE below — the
/// Burton–Miller low-`ka` resistance artifact is a recorded boundary —
/// Burton–Miller above, which owns the fictitious-frequency band).
///
/// # Errors
/// [`HelmholtzError`] on refusal: empty/non-increasing omegas, mask
/// length mismatch, no driven panels, or any per-row solver refusal —
/// including [`HelmholtzError::TooCoarse`] past the mesh's
/// panels-per-wavelength bound (the sweep refuses; it never
/// extrapolates).
pub fn bake_radiation_load(
    surface: &SpherePanels,
    driven: &[bool],
    omegas: &[f64],
    medium: Medium,
    l_max: usize,
    source_id: &str,
) -> Result<RadiationLoadBake, HelmholtzError> {
    let n = surface.areas().len();
    if driven.len() != n {
        return Err(HelmholtzError::ShapeMismatch {
            what: "driven mask length must equal the panel count",
        });
    }
    if omegas.is_empty() {
        return Err(HelmholtzError::BadParameter {
            what: "at least one bake frequency",
        });
    }
    if !omegas
        .windows(2)
        .all(|w| w[0].is_finite() && w[1] > w[0])
        || !(omegas[0].is_finite() && omegas[0] > 0.0)
    {
        return Err(HelmholtzError::BadParameter {
            what: "bake frequencies must be finite, positive, strictly increasing",
        });
    }
    let mouth_area: f64 = surface
        .areas()
        .iter()
        .zip(driven)
        .filter(|&(_, &d)| d)
        .map(|(a, _)| *a)
        .sum();
    let driven_count = driven.iter().filter(|&&d| d).count();
    if driven_count == 0 {
        return Err(HelmholtzError::BadParameter {
            what: "at least one driven mouth panel",
        });
    }
    let mouth_radius = (mouth_area / core::f64::consts::PI).sqrt();
    let velocity: Vec<C64> = driven
        .iter()
        .map(|&d| if d { C64::ONE } else { C64::ZERO })
        .collect();
    let mut rows = Vec::with_capacity(omegas.len());
    let mut fingerprint = 0u64;
    for &omega in omegas {
        let k = omega / medium.sound_speed;
        let mouth_ka = k * mouth_radius;
        let formulation = if mouth_ka < 0.5 {
            Formulation::PlainCbie
        } else {
            Formulation::BurtonMiller
        };
        let solution = solve_radiation(surface, k, medium, &velocity, formulation)?;
        fingerprint = solution.surface_fingerprint;
        // Area-averaged p over the driven panels (v = 1 there): the same
        // arbiter the pulsating-sphere oracle uses.
        let mut num = C64::ZERO;
        for i in 0..n {
            if driven[i] {
                num = num + solution.pressure[i].scale(surface.areas()[i]);
            }
        }
        let z_specific = num.scale(1.0 / mouth_area);
        let z_acoustic = z_specific.scale(1.0 / mouth_area);
        let directivity = directivity_sh_table(surface, &solution, medium, l_max)?;
        rows.push(RadiationLoadRow {
            omega,
            k,
            mouth_ka,
            z_specific,
            z_acoustic,
            formulation,
            panels_per_wavelength: solution.panels_per_wavelength,
            condition_lower_bound: solution.condition_lower_bound,
            radiated_power: solution.radiated_power,
            passivity_margin: z_specific.re,
            directivity,
        });
    }
    Ok(RadiationLoadBake {
        source_id: source_id.to_string(),
        mouth_area_m2: mouth_area,
        panel_count: n,
        driven_count,
        surface_fingerprint: fingerprint,
        medium,
        rows,
    })
}

/// A lathed closed bell fixture: flare wall + throat cap (rigid) and a
/// MOUTH cap whose panels are the driven piston. Panel order follows
/// triangle order: walls, throat cap, then the mouth cap LAST.
#[derive(Debug, Clone)]
pub struct LathedBell {
    /// Oriented outward triangle soup for `SpherePanels::from_triangles`.
    pub triangles: Vec<[[f64; 3]; 3]>,
    /// Index of the first mouth-cap triangle (driven panels are
    /// `mouth_start..triangles.len()`).
    pub mouth_start: usize,
}

impl LathedBell {
    /// Driven mask aligned with the triangle/panel order.
    #[must_use]
    pub fn driven_mask(&self) -> Vec<bool> {
        (0..self.triangles.len())
            .map(|i| i >= self.mouth_start)
            .collect()
    }
}

/// Lathe an axisymmetric `(x, r(x))` flare profile into a CLOSED
/// oriented exterior soup: the flare wall and throat cap are rigid; the
/// mouth disc is the driven piston (the classic bell-mouth radiation
/// fixture). The exterior Helmholtz problem needs a closed surface, so
/// the open mouth is capped BY the piston itself. Windings are outward
/// (verified by the enclosed-volume sign in tests).
///
/// # Errors
/// [`HelmholtzError::BadParameter`] on a degenerate profile.
pub fn lathe_profile(
    profile: &[(f64, f64)],
    circumferential: usize,
) -> Result<LathedBell, HelmholtzError> {
    if profile.len() < 2 || circumferential < 8 {
        return Err(HelmholtzError::BadParameter {
            what: "lathe needs >= 2 profile stations and >= 8 around",
        });
    }
    if !profile
        .windows(2)
        .all(|w| w[1].0 > w[0].0 && w[0].1 > 0.0 && w[1].1 > 0.0)
    {
        return Err(HelmholtzError::BadParameter {
            what: "profile x must increase and radii stay positive",
        });
    }
    let ring = |x: f64, r: f64| -> Vec<[f64; 3]> {
        (0..circumferential)
            .map(|j| {
                let th = core::f64::consts::TAU * j as f64 / circumferential as f64;
                [x, r * th.cos(), r * th.sin()]
            })
            .collect()
    };
    let rings: Vec<Vec<[f64; 3]>> = profile.iter().map(|&(x, r)| ring(x, r)).collect();
    let mut triangles = Vec::new();
    // Wall quads, outward = away from the axis: (lo_j, lo_jn, hi_j) and
    // (lo_jn, hi_jn, hi_j) give normals along (-dr/dx, cos th, sin th).
    for w in rings.windows(2) {
        let (lo, hi) = (&w[0], &w[1]);
        for j in 0..circumferential {
            let jn = (j + 1) % circumferential;
            triangles.push([lo[j], lo[jn], hi[j]]);
            triangles.push([lo[jn], hi[jn], hi[j]]);
        }
    }
    // Throat cap at min x, outward normal -x: fan (center, t_jn, t_j).
    let throat_center = [profile[0].0, 0.0, 0.0];
    let first = &rings[0];
    for j in 0..circumferential {
        let jn = (j + 1) % circumferential;
        triangles.push([throat_center, first[jn], first[j]]);
    }
    // Mouth cap at max x, outward normal +x: fan (center, m_j, m_jn) —
    // the DRIVEN piston, kept LAST so the mask is a suffix.
    let mouth_start = triangles.len();
    let mouth_center = [profile[profile.len() - 1].0, 0.0, 0.0];
    let last = &rings[rings.len() - 1];
    for j in 0..circumferential {
        let jn = (j + 1) % circumferential;
        triangles.push([mouth_center, last[j], last[jn]]);
    }
    Ok(LathedBell {
        triangles,
        mouth_start,
    })
}
