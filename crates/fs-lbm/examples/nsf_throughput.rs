//! Complete-step NSF timing fixture; no recorder, turbulence or audio claim.
//! Run in release mode: nsf_throughput <nx> <ny> <steps> <workers> [workers...].
//! The 50 µm dry-air box is periodic, all fluid and one cell deep. Timing excludes
//! construction, two warmup steps and hashing, but includes every allocation,
//! material query, coupling pass and energy ledger inside `NsfGrid3::step`.
use fs_lbm::d3q27::D3q27;
use fs_lbm::nsf3::{MAX_QUALIFIED_MACH, NsfGrid3, NsfSpec, WallThermal};
use fs_material::gas::GasSpec;
use std::time::Instant;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args = std::env::args()
        .skip(1)
        .map(|s| s.parse::<usize>())
        .collect::<Result<Vec<_>, _>>()?;
    if args.len() < 4 {
        return Err("usage: nsf_throughput <nx> <ny> <steps> <workers> [workers...]".into());
    }
    let (nx, ny, steps) = (args[0], args[1], args[2]);
    let cells = nx.checked_mul(ny).ok_or("cell count overflow")?;
    if nx < 8 || ny < 8 || cells > 10_000_000 || steps == 0 {
        return Err("require nx,ny >= 8, cells <= 10000000 and steps > 0".into());
    }
    if args[3..].iter().any(|&n| n == 0 || n > 128) {
        return Err("worker counts must lie in 1..=128".into());
    }
    let mut expected_hash = None;
    for &workers in &args[3..] {
        println!(
            "NSF_THROUGHPUT begin nx={nx} ny={ny} nz=1 cells={cells} workers={workers} steps={steps} warmup_steps=2 dx_m=0.00005 gas=dry_air_ussa1976 temperature_k=293.15 pressure_pa=101325 periodic=true all_fluid=true"
        );
        let spec = NsfSpec {
            dims: [nx, ny, 1],
            periodic: [true; 3],
            gas: GasSpec::dry_air_ussa1976(),
            temperature_k: 293.15,
            pressure_pa: 101_325.0,
            dx_m: 5e-5,
            temperature_window_k: [273.15, 333.15],
            maximum_mach: MAX_QUALIFIED_MACH,
            maximum_temperature_ratio: 0.2,
            bulk_to_shear_ratio: 0.73,
            bulk_ratio_source: "Cramer 2012, Phys. Fluids 24, 066102: N2 at 293 K",
            maximum_cells: 10_000_000,
            wall: WallThermal::Adiabatic,
            workers,
            transport_points: 33,
        };
        let mut grid = NsfGrid3::<D3q27>::new(
            spec,
            |_| false,
            |x| {
                let px = std::f64::consts::TAU * (x[0] as f64 + 0.5) / nx as f64;
                let py = std::f64::consts::TAU * (x[1] as f64 + 0.5) / ny as f64;
                // A smooth evolving field exercises stress and thermal transport.
                // These dimensionless amplitudes define a timing fixture only.
                let t_ratio = 1.0 + 1e-4 * fs_math::det::cos(px) * fs_math::det::cos(py);
                (
                    1.0 / t_ratio,
                    [
                        0.005 * fs_math::det::sin(py),
                        -0.005 * fs_math::det::sin(px),
                        0.0,
                    ],
                    293.15 * t_ratio,
                )
            },
        )?;
        for _ in 0..2 {
            grid.step()?;
        }
        let start = Instant::now();
        let mut max_residual = 0.0_f64;
        for _ in 0..steps {
            let ledger = grid.step()?;
            if !ledger.residual_rel.is_finite() || ledger.residual_rel > 1e-12 {
                return Err(format!("energy refusal: residual_rel={}", ledger.residual_rel).into());
            }
            max_residual = max_residual.max(ledger.residual_rel);
        }
        let elapsed_s = start.elapsed().as_secs_f64();
        let hash = grid.state_hash();
        if let Some(expected) = expected_hash {
            if hash != expected {
                return Err(format!(
                    "worker-count replay mismatch: {expected:016x} != {hash:016x}"
                )
                .into());
            }
        }
        expected_hash = Some(hash);
        let cell_steps_s = cells as f64 * steps as f64 / elapsed_s;
        println!(
            "NSF_THROUGHPUT result cells={cells} workers={workers} steps={steps} elapsed_s={elapsed_s:.9} complete_cell_steps_s={cell_steps_s:.6} dt_s={:.12e} max_energy_residual_rel={max_residual:e} state_hash={hash:016x} recorder_claim=false",
            grid.law().dt_s
        );
    }
    Ok(())
}
