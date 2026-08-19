//! PrelaunchPhase schema + terrain binding + FlatnessCertificate issuance
//! (bead wf-root-guzez.4.8.3, E3.4-iii). The E1.3 terrain grids feed
//! bilinear height/slope queries; the launch region's plane fit issues
//! (or refuses) the FlatnessCertificate that FlatPlaneVortexImageExact
//! and the rail dynamics require. The INTEGRATED equilibration
//! (HeldOnRailEquilibrated with real aero) is E4.6d — this is the schema
//! + constrained-state scaffold the plan asks for early.

use crate::Refusal;
use crate::rail::RailSpec;
use fs_blake3::hash_domain;

/// Certificate bands (declared; consumers cite the certificate, not raw
/// terrain). Slope as rise/run; residuals in metres.
pub const CERT_MAX_SLOPE: f64 = 0.005;
/// RMS-residual band [m].
pub const CERT_MAX_RMS_M: f64 = 1.2;
/// Grid-size cap per axis.
pub const MAX_GRID_N: usize = 512;
/// Identity domain for prelaunch digests.
pub const PRELAUNCH_DIGEST_DOMAIN: &str = "org.frankensim.fs-flyer.prelaunch.v1";

fn refuse(code: &'static str, message: String, repair: &str) -> Refusal {
    Refusal {
        code,
        message,
        ranked_repairs: vec![repair.into()],
    }
}

/// A row-major terrain heightfield (rows south→north, cols west→east),
/// uniform spacing — the E1.3 grid shape.
#[derive(Clone, Debug, PartialEq)]
pub struct TerrainGrid {
    spacing_m: f64,
    rows: Vec<Vec<f64>>,
}

impl TerrainGrid {
    /// Build from rows; validates shape and finiteness.
    ///
    /// # Errors
    /// `terrain-grid-invalid` (ragged/empty/too large/non-finite; caps at
    /// cap AND cap+1 by construction of the check).
    pub fn new(spacing_m: f64, rows: Vec<Vec<f64>>) -> Result<TerrainGrid, Refusal> {
        if !spacing_m.is_finite() || spacing_m <= 0.0 {
            return Err(refuse(
                "terrain-grid-invalid",
                format!("spacing {spacing_m}"),
                "positive",
            ));
        }
        let n = rows.len();
        if n < 2 || n > MAX_GRID_N {
            return Err(refuse(
                "terrain-grid-invalid",
                format!("{n} rows outside [2, {MAX_GRID_N}]"),
                "the E1.3 grids are 17x17",
            ));
        }
        let m = rows[0].len();
        if m < 2 || m > MAX_GRID_N || rows.iter().any(|r| r.len() != m) {
            return Err(refuse(
                "terrain-grid-invalid",
                "ragged or oversized columns".into(),
                "rectangular grid",
            ));
        }
        if rows.iter().flatten().any(|z| !z.is_finite()) {
            return Err(refuse(
                "terrain-grid-invalid",
                "non-finite elevation".into(),
                "clean the grid",
            ));
        }
        Ok(TerrainGrid { spacing_m, rows })
    }

    /// Bilinear height at (x east, y north) metres from the SW corner.
    ///
    /// # Errors
    /// `terrain-query-outside-domain` (the grid never extrapolates).
    pub fn height_m(&self, x_m: f64, y_m: f64) -> Result<f64, Refusal> {
        let (nr, nc) = (self.rows.len(), self.rows[0].len());
        let (fx, fy) = (x_m / self.spacing_m, y_m / self.spacing_m);
        if !(fx.is_finite() && fy.is_finite())
            || fx < 0.0
            || fy < 0.0
            || fx > (nc - 1) as f64
            || fy > (nr - 1) as f64
        {
            return Err(refuse(
                "terrain-query-outside-domain",
                format!("({x_m}, {y_m}) outside the tile"),
                "the tile never extrapolates; clamp the camera, refuse the physics",
            ));
        }
        let (c0, r0) = (fx.floor() as usize, fy.floor() as usize);
        let (c1, r1) = ((c0 + 1).min(nc - 1), (r0 + 1).min(nr - 1));
        let (tx, ty) = (fx - c0 as f64, fy - r0 as f64);
        let z00 = self.rows[r0][c0];
        let z01 = self.rows[r0][c1];
        let z10 = self.rows[r1][c0];
        let z11 = self.rows[r1][c1];
        Ok(z00 * (1.0 - tx) * (1.0 - ty)
            + z01 * tx * (1.0 - ty)
            + z10 * (1.0 - tx) * ty
            + z11 * tx * ty)
    }
}

/// The issued FlatnessCertificate over a rectangular region.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FlatnessCertificate {
    /// Fitted plane slope components (dz/dx, dz/dy).
    pub slope: (f64, f64),
    /// RMS residual [m].
    pub rms_residual_m: f64,
    /// Max |residual| [m].
    pub max_abs_residual_m: f64,
    /// Plane height at the region's SW corner [m].
    pub z0_m: f64,
}

/// Fit a plane over grid cells [r0..=r1]×[c0..=c1] and ISSUE the
/// certificate iff the declared bands hold.
///
/// # Errors
/// `flatness-region-invalid`;
/// `flatness-uncertifiable` naming the violated band (slope or residual)
/// — the ground-effect image model may then not be used there.
pub fn issue_flatness_certificate(
    grid: &TerrainGrid,
    r0: usize,
    r1: usize,
    c0: usize,
    c1: usize,
) -> Result<FlatnessCertificate, Refusal> {
    let (nr, nc) = (grid.rows.len(), grid.rows[0].len());
    if r0 >= r1 || c0 >= c1 || r1 >= nr || c1 >= nc {
        return Err(refuse(
            "flatness-region-invalid",
            format!("region ({r0}..{r1}, {c0}..{c1}) in {nr}x{nc}"),
            "a non-degenerate in-bounds rectangle",
        ));
    }
    let s = grid.spacing_m;
    // Least-squares plane z = a·x + b·y + c over the region's nodes.
    let (mut sx, mut sy, mut sz) = (0.0f64, 0.0f64, 0.0f64);
    let (mut sxx, mut syy, mut sxy, mut sxz, mut syz) = (0.0f64, 0.0, 0.0, 0.0, 0.0);
    let mut n = 0.0f64;
    for r in r0..=r1 {
        for c in c0..=c1 {
            let (x, y, z) = ((c - c0) as f64 * s, (r - r0) as f64 * s, grid.rows[r][c]);
            sx += x;
            sy += y;
            sz += z;
            sxx += x * x;
            syy += y * y;
            sxy += x * y;
            sxz += x * z;
            syz += y * z;
            n += 1.0;
        }
    }
    // Cramer solve of the 3x3 normal system.
    let det3 = |m: [[f64; 3]; 3]| -> f64 {
        m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
            - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
            + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0])
    };
    let a_m = [[sxx, sxy, sx], [sxy, syy, sy], [sx, sy, n]];
    let rhs = [sxz, syz, sz];
    let d = det3(a_m);
    let col = |i: usize| -> [[f64; 3]; 3] {
        let mut m = a_m;
        for r in 0..3 {
            m[r][i] = rhs[r];
        }
        m
    };
    let (a, b, c_off) = (det3(col(0)) / d, det3(col(1)) / d, det3(col(2)) / d);
    let (mut ss, mut mx) = (0.0f64, 0.0f64);
    for r in r0..=r1 {
        for c in c0..=c1 {
            let (x, y, z) = ((c - c0) as f64 * s, (r - r0) as f64 * s, grid.rows[r][c]);
            let resid = z - (a * x + b * y + c_off);
            ss += resid * resid;
            mx = mx.max(resid.abs());
        }
    }
    let rms = (ss / n).sqrt();
    let slope_mag = (a * a + b * b).sqrt();
    if slope_mag > CERT_MAX_SLOPE {
        return Err(refuse(
            "flatness-uncertifiable",
            format!("slope {slope_mag:.5} exceeds the {CERT_MAX_SLOPE} band"),
            "the flat-plane image model may not be used here; pick a flatter region",
        ));
    }
    if rms > CERT_MAX_RMS_M {
        return Err(refuse(
            "flatness-uncertifiable",
            format!("RMS residual {rms:.3} m exceeds the {CERT_MAX_RMS_M} m band"),
            "dune texture too strong for a certified flat plane",
        ));
    }
    Ok(FlatnessCertificate {
        slope: (a, b),
        rms_residual_m: rms,
        max_abs_residual_m: mx,
        z0_m: c_off,
    })
}

/// The PrelaunchPhase schema (scaffold; integrated equilibration = E4.6d).
#[derive(Clone, Debug, PartialEq)]
pub struct PrelaunchPhase {
    /// The rail specification.
    pub rail: RailSpec,
    /// Mean headwind sampled for the hold [m/s].
    pub headwind_mps: f64,
    /// Neutral control commands at the hold (canard, warp) [rad].
    pub controls_rad: [f64; 2],
    /// The launch region's certificate (issuance is admission).
    pub flatness: FlatnessCertificate,
}

impl PrelaunchPhase {
    /// Admit the phase (rail spec + finite fields; the certificate was
    /// already issued under its own bands).
    ///
    /// # Errors
    /// Rail refusals; `prelaunch-invalid` (non-finite/negative headwind).
    pub fn admit(&self) -> Result<(), Refusal> {
        self.rail.admit()?;
        if !self.headwind_mps.is_finite()
            || self.headwind_mps < 0.0
            || !self.controls_rad.iter().all(|v| v.is_finite())
        {
            return Err(refuse(
                "prelaunch-invalid",
                format!(
                    "headwind {}, controls {:?}",
                    self.headwind_mps, self.controls_rad
                ),
                "finite non-negative headwind; finite neutral controls",
            ));
        }
        Ok(())
    }

    /// Content digest (identity ingredient: the mode-complete tick-0 state
    /// freezes BEFORE RunIntentId mints; this digest covers the SCHEMA
    /// fields — the full state digest joins in E4.6d).
    #[must_use]
    pub fn digest(&self) -> String {
        let mut payload = Vec::new();
        for v in [
            self.rail.z_rail_m,
            self.rail.length_m,
            f64::from(self.rail.hysteresis_ticks),
            self.headwind_mps,
            self.controls_rad[0],
            self.controls_rad[1],
            self.flatness.slope.0,
            self.flatness.slope.1,
            self.flatness.rms_residual_m,
            self.flatness.z0_m,
        ] {
            payload.extend_from_slice(&v.to_bits().to_le_bytes());
        }
        hash_domain(PRELAUNCH_DIGEST_DOMAIN, &payload).to_hex()
    }
}
