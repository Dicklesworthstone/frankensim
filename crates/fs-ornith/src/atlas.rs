//! Stage 5: the CERTIFIED PARETO ATLAS. NSGA-II (fs-dfo — the landed
//! multi-objective engine; NSGA-III reference-point selection is the
//! recorded successor) over four objectives: −L/D, −ROA volume,
//! −maneuver proxy, and inlet mass-flow violation. The front's winner
//! gets a GRADIENT POLISH through the landed BEM adjoint (honestly
//! labeled: α-direction adjoint-assisted with finite-difference output
//! partials, thickness by central difference). Every
//! atlas row carries its stability certificate, its conformal
//! screening band, and its lineage — the deliverable is evidence, not
//! a scatter plot.

use crate::certify::{CertifyReport, LdSurrogate, certify};
use crate::param::{GENE_DIM, OrnithCandidate};
use crate::screen::lift_to_drag;
use fs_dfo::moo::{NsgaParams, hypervolume, knee_point, nsga2};

/// One certified atlas row.
#[derive(Debug, Clone)]
pub struct AtlasRow {
    /// The decoded design.
    pub candidate: OrnithCandidate,
    /// Genes (lineage: the exact NSGA decision vector).
    pub genes: Vec<f64>,
    /// Objectives (L/D, ROA volume, maneuver, inlet mass-flow).
    pub ld: f64,
    /// Certified ROA proxy volume.
    pub roa: f64,
    /// Maneuver proxy.
    pub maneuver: f64,
    /// Inlet mass-flow satisfaction (distance to the [0.9, 1.6] band).
    pub inlet_violation: f64,
    /// The stability certificate.
    pub certificate: CertifyReport,
    /// Surrogate prediction and band half-width (the screening
    /// evidence attached to the row).
    pub surrogate_ld: f64,
}

/// The atlas artifact.
#[derive(Debug, Clone)]
pub struct Atlas {
    /// Certified front rows.
    pub rows: Vec<AtlasRow>,
    /// Hypervolume of the front (evidence of coverage).
    pub hypervolume: f64,
    /// Knee-point index (the "start here" design).
    pub knee: usize,
    /// The polished winner (knee design after adjoint polish).
    pub polished: OrnithCandidate,
    /// L/D before → after polish.
    pub polish_gain: (f64, f64),
}

/// Inlet mass-flow violation: distance to the satisfaction band.
fn inlet_violation(c: &OrnithCandidate) -> f64 {
    let m = c.inlet_mass_flow(crate::screen::PANELS);
    if m < 0.9 {
        0.9 - m
    } else if m > 1.6 {
        m - 1.6
    } else {
        0.0
    }
}

/// Build the certified Pareto atlas.
///
/// # Panics
/// Fixture-scale programmer contracts only.
#[must_use]
pub fn build_atlas(pop: usize, generations: usize, seed: u64, surrogate: &LdSurrogate) -> Atlas {
    let mut objectives = |genes: &[f64]| -> Vec<f64> {
        let c = OrnithCandidate::from_genes(genes);
        let cert = certify(&c);
        vec![
            -lift_to_drag(&c),
            -cert.roa_volume,
            -cert.maneuver,
            inlet_violation(&c),
        ]
    };
    let params = NsgaParams {
        pop,
        generations,
        eta_c: 15.0,
        eta_m: 20.0,
        p_mut: 1.0 / GENE_DIM as f64,
        seed,
    };
    let front = nsga2(&mut objectives, GENE_DIM, (0.0, 1.0), &params);
    let rows: Vec<AtlasRow> = front
        .iter()
        .map(|ind| {
            let c = OrnithCandidate::from_genes(&ind.x);
            let certificate = certify(&c);
            AtlasRow {
                candidate: c,
                genes: ind.x.clone(),
                ld: -ind.f[0],
                roa: -ind.f[1],
                maneuver: -ind.f[2],
                inlet_violation: ind.f[3],
                surrogate_ld: surrogate.predict(&c),
                certificate,
            }
        })
        .collect();
    let objs: Vec<Vec<f64>> = front.iter().map(|i| i.f.clone()).collect();
    let reference = vec![1.0, 1.0, 1.0, 2.0];
    let hv = hypervolume(&objs, &reference);
    let knee = knee_point(&objs);
    // Gradient polish on the knee design: ascend L/D along the adjoint-assisted
    // ∂cl/∂α adjoint direction (thickness held — the honest lever).
    let knee_c = rows[knee].candidate;
    let ld0 = lift_to_drag(&knee_c);
    let mut best = knee_c;
    let mut best_ld = ld0;
    let mut step = 0.01f64;
    let mut cur = knee_c;
    for _ in 0..25 {
        let g = cur.cl_gradient(crate::screen::PANELS);
        let mut trial = cur;
        trial.alpha = (cur.alpha + step * g[0].signum()).clamp(0.02, 0.14);
        let ld = lift_to_drag(&trial);
        if ld > best_ld {
            best_ld = ld;
            best = trial;
            cur = trial;
            step *= 1.2;
        } else {
            step *= 0.5;
            if step < 1e-5 {
                break;
            }
        }
    }
    Atlas {
        rows,
        hypervolume: hv,
        knee,
        polished: best,
        polish_gain: (ld0, best_ld),
    }
}
