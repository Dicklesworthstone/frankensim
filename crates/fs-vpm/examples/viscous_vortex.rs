//! Run with `cargo run -p fs-vpm --example viscous_vortex`.
//! A declared SI viscosity broadens a Gaussian vortex through real PSE/RK4 flow.

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_vpm::{
    VortexParticle, VpmBudget,
    diffusion::{PseDiffusion, simulate_viscous_with_cx},
    total_circulation,
};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let cells = 21;
    let spacing_m = 0.0008;
    let initial_width_m2 = 0.002_f64.powi(2);
    let circulation_m2_s = 1e-6;
    // Explicit scenario coefficient, not a measured material specification.
    let nu_m2_s = 1.5e-5;
    let dt_s = 0.001;
    let steps = 10;
    let law = PseDiffusion::new(nu_m2_s, 0.0016)?;
    let mut particles = Vec::with_capacity(cells * cells);
    for y in 0..cells {
        for x in 0..cells {
            let position = [(x as f64 - 10.0) * spacing_m, (y as f64 - 10.0) * spacing_m];
            let r2 = position[0] * position[0] + position[1] * position[1];
            let vorticity = circulation_m2_s / (std::f64::consts::PI * initial_width_m2)
                * fs_math::det::exp(-r2 / initial_width_m2);
            particles.push(VortexParticle::new(
                position,
                vorticity * spacing_m * spacing_m,
            ));
        }
    }
    let areas = vec![spacing_m * spacing_m; particles.len()];
    let pool = ArenaPool::new(ArenaConfig::default());
    let gate = CancelGate::new_clock_free();
    let result = pool.scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 0,
                kernel_id: 7,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        simulate_viscous_with_cx(
            &cx,
            &particles,
            &areas,
            dt_s,
            steps,
            0.0004,
            law,
            VpmBudget::new(particles.len(), steps, 20_000_000, 200_000),
        )
    })?;
    let circulation = total_circulation(&result.particles);
    let width_m2 = result
        .particles
        .iter()
        .map(|p| p.circulation * (p.pos[0] * p.pos[0] + p.pos[1] * p.pos[1]))
        .sum::<f64>()
        / circulation;
    let enstrophy = |state: &[VortexParticle]| {
        state
            .iter()
            .zip(&areas)
            .map(|(p, a)| 0.5 * p.circulation * p.circulation / a)
            .sum::<f64>()
    };
    println!("declared kinematic viscosity: {nu_m2_s:.6e} m²/s");
    println!("elapsed: {:.6e} s", dt_s * steps as f64);
    println!(
        "circulation: {:.12e} -> {circulation:.12e} m²/s",
        total_circulation(&particles)
    );
    println!("mean radius squared: {initial_width_m2:.8e} -> {width_m2:.8e} m²");
    println!(
        "Lamb–Oseen continuum width: {:.8e} m²",
        initial_width_m2 + 4.0 * nu_m2_s * dt_s * steps as f64
    );
    println!(
        "enstrophy: {:.8e} -> {:.8e} m²/s²",
        enstrophy(&particles),
        enstrophy(&result.particles)
    );
    println!(
        "pair visits: {}; logical work: {}; live payload: {} bytes",
        result.pair_evaluations, result.work_units, result.peak_live_bytes
    );
    Ok(())
}
