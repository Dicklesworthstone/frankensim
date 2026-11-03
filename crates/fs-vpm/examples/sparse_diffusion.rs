//! Run with cargo run -p fs-vpm --example sparse_diffusion.
//! Reuse one compact PSE graph for two successive fixed-site Gaussian-vortex
//! diffusion runs. Compare the physical field to the continuum heat solution
//! and report actual pair work, without a machine-speed or accuracy certificate.

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_math::det;
use fs_vpm::diffusion::{PreparationBudget, PreparedPse, PseDiffusion};
use fs_vpm::{VortexParticle, VpmBudget, VpmError, total_circulation};
use std::f64::consts::PI;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let side = 41_usize;
    let spacing_m = 0.15;
    let initial_width_m2 = 1.0;
    let circulation_m2_s = 1.0;
    // Declared scenario coefficient, not a measured fluid specification.
    let nu_m2_s = 0.01;
    let support_radius_m = 4.0 * spacing_m;
    let dt_s = 0.005;
    let steps_per_run = 50;
    let area_m2 = spacing_m * spacing_m;
    let mut particles = Vec::with_capacity(side * side);
    for y in 0..side {
        for x in 0..side {
            let pos = [(x as f64 - 20.0) * spacing_m, (y as f64 - 20.0) * spacing_m];
            let r2 = pos[0] * pos[0] + pos[1] * pos[1];
            let omega = circulation_m2_s / (PI * initial_width_m2)
                * det::exp(-r2 / initial_width_m2);
            particles.push(VortexParticle::new(pos, area_m2 * omega));
        }
    }
    // Zero-strength support would also be retained. No particles or Gaussian
    // tails are discarded by a strength threshold.
    let areas = vec![area_m2; particles.len()];
    let n = particles.len();
    let max_edges = 32 * n;
    let law = PseDiffusion::compact_wendland(nu_m2_s, support_radius_m)?;
    let preparation_budget = PreparationBudget::new(
        n,
        max_edges,
        1000 * n as u64,
        5000 * n as u64,
        2 * 1024 * 1024,
    );
    let run_budget = VpmBudget::new(
        n,
        steps_per_run,
        2 * max_edges as u64 * steps_per_run as u64,
        1024 * 1024,
    );
    let pool = ArenaPool::new(ArenaConfig::default());
    let gate = CancelGate::new_clock_free();
    let (result, edges, preparation, pair_visits) = pool.scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 0,
                kernel_id: 8,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        let prepared =
            PreparedPse::prepare_with_cx(&cx, &particles, &areas, law, preparation_budget)?;
        let first =
            prepared.diffuse_with_cx(&cx, &particles, dt_s, steps_per_run, run_budget)?;
        // Geometry is unchanged by diffusion; the same immutable graph admits
        // the evolved circulation inventory without repeating spatial searches.
        let second =
            prepared.diffuse_with_cx(&cx, &first.particles, dt_s, steps_per_run, run_budget)?;
        let pair_visits = first.pair_evaluations + second.pair_evaluations;
        Ok::<_, VpmError>((
            second,
            prepared.edges().len(),
            *prepared.preparation(),
            pair_visits,
        ))
    })?;

    let time_s = dt_s * (2 * steps_per_run) as f64;
    let continuum_width_m2 = initial_width_m2 + 4.0 * nu_m2_s * time_s;
    let moment = |state: &[VortexParticle]| {
        state
            .iter()
            .map(|p| p.circulation * (p.pos[0] * p.pos[0] + p.pos[1] * p.pos[1]))
            .sum::<f64>()
            / total_circulation(state)
    };
    let (mut error2, mut norm2) = (0.0, 0.0);
    for particle in &result.particles {
        let r2 = particle.pos[0] * particle.pos[0] + particle.pos[1] * particle.pos[1];
        let exact = circulation_m2_s / (PI * continuum_width_m2)
            * det::exp(-r2 / continuum_width_m2);
        let difference = particle.circulation / area_m2 - exact;
        error2 += area_m2 * difference * difference;
        norm2 += area_m2 * exact * exact;
    }
    let dense_pair_visits = n as u64 * (n - 1) as u64 * (2 * steps_per_run) as u64;
    println!("particles: {n}; supported unordered edges: {edges}");
    println!("declared viscosity: {nu_m2_s:.6e} m2/s; elapsed: {time_s:.6e} s");
    println!(
        "preparation: {} candidate visits; {} logical units; {} payload bytes",
        preparation.candidate_visits, preparation.work_units, preparation.peak_live_bytes
    );
    println!(
        "circulation: {:.12e} -> {:.12e} m2/s",
        total_circulation(&particles), total_circulation(&result.particles)
    );
    println!(
        "mean radius squared: {:.8e} -> {:.8e} m2; continuum width: {continuum_width_m2:.8e} m2",
        moment(&particles), moment(&result.particles)
    );
    println!("relative Gaussian heat-field L2 error: {:.8e}", det::sqrt(error2 / norm2));
    println!(
        "diffusion pair visits: {pair_visits}; dense compact equivalent: {dense_pair_visits}"
    );
    println!("final run live payload: {} bytes", result.peak_live_bytes);
    Ok(())
}
