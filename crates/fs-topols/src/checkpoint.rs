//! Exact durable checkpoints for the single-load level-set optimizer.
//!
//! The fixed-run optimizer remains unchanged. This module reuses the same
//! corrected numerical evolution in a segment primitive that accepts an explicit
//! global iteration ordinal and current augmented-Lagrange multiplier, allowing
//! durable continuation without replaying prior geometry updates.

use crate::fim::{RedistanceAudit, redistance};
use crate::gridsdf::GridSdf;
use crate::optimize::{Cantilever, OptimizeReport, OptimizeSettings, material_volume};
use crate::topder::{NucleationEvent, nucleate, topological_derivative};
use crate::veloext::extend_velocity;
use crate::weno::{Velocity, advect, build_band};
use fs_cutfem::{
    BoundaryTraction, CutElasticity, CutElasticitySolution, CutFemError, DesignBoxEdge, EdgeBand,
    MAX_PLANE_STRAIN_STIFFNESS_RATIO, Quadtree,
};
use fs_material::IsotropicElastic;
use std::fmt::Write as _;

const MATERIAL_STRAIN_LIMIT: f64 = 1.0;
const SOLVER_TOL: f64 = 1e-12;
const SOLVER_MAX_ITERS: usize = 60_000;

fn fnv(phi: &GridSdf) -> u64 {
    let mut h = 0xcbf2_9ce4_8422_2325u64;
    for value in phi.nodes() {
        for byte in value.to_bits().to_le_bytes() {
            h ^= u64::from(byte);
            h = h.wrapping_mul(0x0000_0100_0000_01b3);
        }
    }
    h
}

fn invalid_input(what: impl Into<String>) -> CutFemError {
    CutFemError::InvalidElasticityInput { what: what.into() }
}

fn cantilever_support(fixture: Cantilever) -> Result<EdgeBand, CutFemError> {
    if !(fixture.load.is_finite() && fixture.load > 0.0) {
        return Err(invalid_input(format!(
            "cantilever load magnitude {} must be finite and strictly positive",
            fixture.load
        )));
    }
    EdgeBand::new(DesignBoxEdge::Right, 0.5 - fixture.band, 0.5 + fixture.band).map_err(|error| {
        invalid_input(format!(
            "cantilever load-band half-width {} must be finite and lie in [0, 0.5]: {error}",
            fixture.band
        ))
    })
}

fn retain_cantilever_load_pad(phi: &mut GridSdf, support: EdgeBand) -> usize {
    let n = phi.n();
    let h = phi.h();
    let x_min = 1.0 - 2.0 * h;
    let x_max = 1.0 + h;
    let y_min = support.start() - h;
    let y_max = support.end() + h;
    let mut changed = 0usize;
    for j in 0..=n {
        for i in n.saturating_sub(2)..=n {
            let [x, y] = phi.pos(i, j);
            let pad = (x_min - x).max(x - x_max).max(y_min - y).max(y - y_max);
            if pad <= 0.0 && pad < phi.node(i, j) {
                *phi.node_mut(i, j) = pad;
                changed += 1;
            }
        }
    }
    changed
}

fn validated_plane_strain_material(
    settings: OptimizeSettings,
) -> Result<(IsotropicElastic, f64, f64), CutFemError> {
    if !(settings.youngs.is_finite() && settings.youngs > 0.0) {
        return Err(invalid_input(format!(
            "optimizer Young's modulus {} must be finite and positive",
            settings.youngs
        )));
    }
    if !(settings.poisson.is_finite() && settings.poisson > -1.0 && settings.poisson < 0.5) {
        return Err(invalid_input(format!(
            "optimizer Poisson ratio {} must lie in (-1, 0.5)",
            settings.poisson
        )));
    }
    let material = IsotropicElastic::new(settings.youngs, settings.poisson, MATERIAL_STRAIN_LIMIT)
        .map_err(|error| invalid_input(format!("optimizer material card was refused: {error}")))?;
    let (lambda, mu) = material.lame();
    let bulk_2d = lambda + mu;
    let stiffness_ratio = (lambda + 2.0 * mu) / mu;
    if !(lambda.is_finite()
        && mu.is_finite()
        && mu > 0.0
        && bulk_2d.is_finite()
        && bulk_2d > 0.0
        && stiffness_ratio.is_finite())
    {
        return Err(invalid_input(
            "optimizer material does not define a finite coercive plane-strain law",
        ));
    }
    if stiffness_ratio > MAX_PLANE_STRAIN_STIFFNESS_RATIO {
        return Err(invalid_input(format!(
            "optimizer plane-strain stiffness ratio (lambda + 2*mu)/mu = {stiffness_ratio} exceeds the certified limit {MAX_PLANE_STRAIN_STIFFNESS_RATIO}"
        )));
    }
    Ok((material, lambda, mu))
}

fn mass_stiffness(n: usize) -> (fs_sparse::Csr, fs_sparse::Csr) {
    #[allow(clippy::cast_precision_loss)]
    let h = 1.0 / n as f64;
    let stride = n + 1;
    let nn = stride * stride;
    let mut mc = fs_sparse::Coo::new(nn, nn);
    let mut kc = fs_sparse::Coo::new(nn, nn);
    let ke = [
        [2.0 / 3.0, -1.0 / 6.0, -1.0 / 3.0, -1.0 / 6.0],
        [-1.0 / 6.0, 2.0 / 3.0, -1.0 / 6.0, -1.0 / 3.0],
        [-1.0 / 3.0, -1.0 / 6.0, 2.0 / 3.0, -1.0 / 6.0],
        [-1.0 / 6.0, -1.0 / 3.0, -1.0 / 6.0, 2.0 / 3.0],
    ];
    for cj in 0..n {
        for ci in 0..n {
            let ids = [
                ci + cj * stride,
                ci + 1 + cj * stride,
                ci + 1 + (cj + 1) * stride,
                ci + (cj + 1) * stride,
            ];
            for (a, &ia) in ids.iter().enumerate() {
                mc.push(ia, ia, 0.25 * h * h);
                for (b, &ib) in ids.iter().enumerate() {
                    kc.push(ia, ib, ke[a][b]);
                }
            }
        }
    }
    (mc.assemble(), kc.assemble())
}

#[path = "optimize/engine_checkpoint.rs"]
mod engine;
#[path = "optimize/stateful.rs"]
mod stateful;

pub use stateful::OptimizeCheckpoint;
