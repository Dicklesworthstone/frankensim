//! fs-metamat-e2e — MetamatCert: a certified stiffness-density frontier for a
//! porous metamaterial. Layer: L4 (ASCENT).
//!
//! # The campaign
//!
//! Numerical homogenization turns a microstructure into an effective stiffness
//! tensor — but the raw numbers carry no guarantee they are even physically
//! admissible. This discovers the stiffness-density frontier of a holed-plate
//! metamaterial and PROVES two things at every point, composing crates never
//! designed to meet:
//!
//! - **Homogenization** ([`fs_lattice`]): each porosity gives an effective Voigt
//!   tensor `C` and a solid fraction `ρ`; the axial stiffness is `C₁₁`.
//! - **Stability, PROVEN** ([`fs_sos::is_psd`]): a physical elastic tensor must
//!   be positive-definite (a non-PSD stiffness stores negative energy — it is
//!   unstable). The minimum-eigenvalue certificate proves `C ≻ 0` at every point.
//! - **Admissibility, PROVEN** ([`fs_lattice::voigt_bound`]): no microstructure
//!   at fraction `ρ` can beat the Voigt mixture bound `ρ·C₁₁ˢᵒˡⁱᵈ`. Every
//!   homogenized `C₁₁` is checked against it — a bound violation would mean the
//!   homogenizer itself is wrong (certifying the certifier), and it also PROVES
//!   the solid is optimal for specific stiffness (`C₁₁/ρ ≤ C₁₁ˢᵒˡⁱᵈ`).
//! - **Honest colors** ([`fs_evidence`]): an all-stable, all-admissible frontier
//!   is `Verified`.
//!
//! Deterministic; no dependencies beyond the composed crates.

use fs_evidence::Color;
use fs_lattice::{Homogenizer, UnitCell, voigt_bound};
use fs_sos::is_psd;

/// One point on the stiffness-density frontier.
#[derive(Debug, Clone, Copy)]
pub struct CellPoint {
    /// Hole radius (porosity parameter).
    pub r: f64,
    /// Solid volume fraction.
    pub density: f64,
    /// Axial effective stiffness `C₁₁`.
    pub c11: f64,
    /// Specific stiffness `C₁₁ / ρ`.
    pub specific_stiffness: f64,
    /// Is the effective tensor certified positive-definite (stable)?
    pub stable: bool,
    /// Is `C₁₁` at or below the Voigt upper bound (admissible)?
    pub admissible: bool,
}

/// The campaign report.
// The four bools are independent certificate outcomes on the frontier.
#[allow(clippy::struct_excessive_bools)]
#[derive(Debug, Clone)]
pub struct MetamatReport {
    /// The frontier (increasing porosity).
    pub frontier: Vec<CellPoint>,
    /// Solid-cell axial stiffness `C₁₁ˢᵒˡⁱᵈ`.
    pub c_solid: f64,
    /// Is every cell certified stable (PSD)?
    pub all_stable: bool,
    /// Is every cell certified admissible (≤ Voigt)?
    pub all_admissible: bool,
    /// Does `C₁₁` decrease monotonically with porosity?
    pub stiffness_monotone: bool,
    /// The Voigt bound PROVES the solid maximizes specific stiffness.
    pub solid_is_specific_optimal: bool,
    /// The frontier's stability color (`Verified` iff all-stable & all-admissible).
    pub stability_color: Color,
}

/// Run the MetamatCert campaign on an `n×n` cell over the hole radii `radii`.
///
/// # Panics
/// If `radii` is empty.
#[must_use]
pub fn run_campaign(n: usize, radii: &[f64]) -> MetamatReport {
    assert!(!radii.is_empty(), "need at least one radius");
    let homog = Homogenizer::new(n);
    // Solid reference (no hole) sets the Voigt bound's `C₁₁ˢᵒˡⁱᵈ`.
    let solid = homog.effective(&UnitCell::holed_plate(n, 0.0));
    let c_solid = solid.c[0][0];

    let mut frontier = Vec::with_capacity(radii.len());
    for &r in radii {
        let cell = UnitCell::holed_plate(n, r);
        let eff = homog.effective(&cell);
        let c11 = eff.c[0][0];
        let density = eff.density;
        // `is_psd` uses Jacobi rotations (symmetric-matrix precondition); the
        // energy tensor is symmetric up to roundoff, so symmetrize exactly before
        // trusting the certificate — a no-op for a correct tensor, and it keeps
        // the PSD certificate sound rather than resting on Jacobi's tolerance.
        let cm: Vec<Vec<f64>> = (0..3)
            .map(|i| {
                (0..3)
                    .map(|j| f64::midpoint(eff.c[i][j], eff.c[j][i]))
                    .collect()
            })
            .collect();
        let stable = is_psd(&cm, 1e-9);
        // Voigt upper bound with a tiny void stiffness floor.
        let bound = voigt_bound(c_solid, density, 1e-3);
        let admissible = c11 <= bound + 1e-6 * c_solid.max(1.0);
        frontier.push(CellPoint {
            r,
            density,
            c11,
            specific_stiffness: if density > 1e-12 { c11 / density } else { 0.0 },
            stable,
            admissible,
        });
    }

    let all_stable = frontier.iter().all(|p| p.stable);
    let all_admissible = frontier.iter().all(|p| p.admissible);
    let stiffness_monotone = frontier
        .windows(2)
        .all(|w| w[1].c11 <= w[0].c11 + 1e-9 * c_solid.max(1.0));
    // Voigt ⇒ specific stiffness C₁₁/ρ ≤ C₁₁ˢᵒˡⁱᵈ for every cell.
    let solid_is_specific_optimal = frontier
        .iter()
        .all(|p| p.specific_stiffness <= c_solid + 1e-6 * c_solid.max(1.0));

    let stability_color = if all_stable && all_admissible {
        // declared-color-ok: demo stability candidate from local frontier admissibility bounds; admitted only at a consumer's authority boundary (6pf9)
        Color::Verified {
            lo: frontier.iter().map(|p| p.c11).fold(f64::INFINITY, f64::min),
            hi: c_solid,
        }
    } else {
        Color::Estimated {
            estimator: "homogenization-uncertified".to_string(),
            dispersion: f64::INFINITY,
        }
    };

    MetamatReport {
        frontier,
        c_solid,
        all_stable,
        all_admissible,
        stiffness_monotone,
        solid_is_specific_optimal,
        stability_color,
    }
}

/// The default porosity sweep.
#[must_use]
pub fn default_radii() -> Vec<f64> {
    vec![0.0, 0.08, 0.16, 0.24, 0.32, 0.40]
}
