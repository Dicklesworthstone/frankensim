//! Periodic D3Q27 Taylor–Green evolution with explicit exploratory SGS.
//! Re = U/(k nu), k = 2*pi on the unit box; U = 1 defines normalization.
//! Keep Re and lattice velocity fixed and double steps when doubling grid
//! size to compare the same convective time at finer physical filter width.
//! This is an isothermal low-Mach experiment, not a validated turbulent DNS.

use fs_exec::{CancelGate, Cx};
use fs_lbm::d3q27::{
    CollisionConfig27, CollisionModel27, D3q27, Equilibrium27, Rates27, equilibrium27,
};
use fs_lbm::lattice3::{Grid3, GridBudget3};
use fs_lbm::sgs3::{SgsConfig, SgsModel, SgsUse};
use fs_math::det;
use std::io::{BufWriter, Write};

const HELP: &str = "Periodic Taylor-Green experiment; CSV goes to stdout.
Run: cargo run --release -p fs-lbm --example turbulent_taylor_green -- [options]
  --grid N                 cells per side, 4..128 (default 24)
  --steps N                steps, 0..1000000 (default 120; 0 reports initial state)
  --reynolds R             Re = U/(k nu), k=2*pi, U=1 (default 1600)
  --velocity U             lattice velocity, 0 < U <= 0.08 (default 0.04)
  --model M                none | smagorinsky | wale (default none)
  --coefficient C          nonnegative (defaults: none=0, smagorinsky=0.12, wale=0.5)
  --max-eddy-ratio R        maximum nu_sgs/nu, refusal not clipping (default 1000)
  --help                   show this help
Example refinement at matched convective time (repeat with --model none):
  --grid 16 --steps 80  --reynolds 1600 --velocity 0.04 --model wale --coefficient 0.5
  --grid 32 --steps 160 --reynolds 1600 --velocity 0.04 --model wale --coefficient 0.5
  --grid 64 --steps 320 --reynolds 1600 --velocity 0.04 --model wale --coefficient 0.5
Compare the CSV dissipation ratio and energy histories; convergence is measured,
not presumed. A coarse high-Re run is not a DNS claim.\n";

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let (mut n, mut steps) = (24_usize, 120_usize);
    let (mut re, mut speed) = (1600.0_f64, 0.04_f64);
    let (mut label, mut coefficient, mut cap) = (String::from("none"), None, 1000.0_f64);
    let mut args = std::env::args().skip(1);
    while let Some(arg) = args.next() {
        if arg == "--help" {
            print!("{HELP}");
            return Ok(());
        }
        let value = args
            .next()
            .ok_or_else(|| format!("missing value for {arg}; use --help"))?;
        match arg.as_str() {
            "--grid" => n = value.parse()?,
            "--steps" => steps = value.parse()?,
            "--reynolds" => re = value.parse()?,
            "--velocity" => speed = value.parse()?,
            "--model" => label = value,
            "--coefficient" => coefficient = Some(value.parse::<f64>()?),
            "--max-eddy-ratio" => cap = value.parse()?,
            _ => return Err(format!("unknown option {arg}; use --help").into()),
        }
    }
    if !(4..=128).contains(&n) || steps > 1_000_000 {
        return Err("require grid in 4..128 and steps in 0..1000000".into());
    }
    if !(re.is_finite() && re > 0.0 && speed.is_finite() && speed > 0.0 && speed <= 0.08) {
        return Err("require finite Re > 0 and 0 < lattice velocity <= 0.08".into());
    }
    let coefficient = coefficient.unwrap_or(match label.as_str() {
        "none" => 0.0,
        "smagorinsky" => 0.12,
        "wale" => 0.5,
        _ => return Err("model must be none, smagorinsky, or wale".into()),
    });
    let model = match label.as_str() {
        "none" if coefficient == 0.0 => SgsModel::Smagorinsky { cs: 0.0 },
        "smagorinsky" => SgsModel::Smagorinsky { cs: coefficient },
        "wale" => SgsModel::Wale { cw: coefficient },
        _ => return Err("the none control requires coefficient=0".into()),
    };
    let config = SgsConfig {
        model,
        filter_width_lu: 1.0,
        maximum_eddy_to_molecular: cap,
        usage: SgsUse::Exploration,
    };
    config.validate()?;
    let k = core::f64::consts::TAU;
    let cells = n * n * n;
    let nu = speed * n as f64 / (k * re);
    let tau = 0.5 + 3.0 * nu;
    let mut rates = Rates27 { omega: [1.0; 10] };
    rates.omega[0] = 1.0 / tau;
    let collision = CollisionConfig27 {
        model: CollisionModel27::Cumulant { rates },
        maximum_mach: 0.2,
    };
    collision.validate()?;
    let mut grid = Grid3::<D3q27, 27>::try_new(
        [n; 3],
        collision,
        [0.0; 3],
        [true; 3],
        GridBudget3 {
            maximum_cells: cells,
            maximum_mach: 0.2,
        },
        |_| false,
        |p| {
            let [x, y, z] = p.map(|i| k * (i as f64 + 0.5) / n as f64);
            let u = [
                speed * det::sin(x) * det::cos(y) * det::cos(z),
                -speed * det::cos(x) * det::sin(y) * det::cos(z),
                0.0,
            ];
            // Incompressible TGV pressure initializes rho through p=cs² rho
            // (cs²=1/3), reducing the uniform-pressure acoustic startup error.
            let rho = 1.0
                + 3.0 * speed * speed / 16.0
                    * (det::cos(2.0 * x) + det::cos(2.0 * y))
                    * (det::cos(2.0 * z) + 2.0);
            equilibrium27(rho, u, Equilibrium27::Factorized)
        },
    )?;
    let mut workspace = grid.sgs_workspace()?;
    let gate = CancelGate::new_clock_free();
    let arenas = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    let mut output = BufWriter::new(std::io::stdout().lock());
    writeln!(
        output,
        "# grid={n} steps={steps} Re={re} U_lu={speed} model={label} coefficient={coefficient} max_eddy_ratio={cap}"
    )?;
    writeln!(
        output,
        "# unit_periodic_box=1 U_reference=1 rho_reference=1 k={k} nu_reference={} dx={} dt={} dt_star={} tau={tau}",
        1.0 / (k * re),
        1.0 / n as f64,
        speed / n as f64,
        k * speed / n as f64
    )?;
    writeln!(
        output,
        "# collision=cumulant non_shear_rates=1 maximum_mach=0.2 filter_width_cells=1 mode=exploration"
    )?;
    writeln!(
        output,
        "# time_star=k*U*t; E_star=<rho*u^2>/(2*rho_ref*U^2); eps_star=<2*rho*nu*Sdev:Sdev>/(rho_ref*U^3*k)"
    )?;
    writeln!(
        output,
        "# instantaneous shear diagnostics exclude bulk/numerical dissipation; equilibrium startup; no DNS-validation claim"
    )?;
    writeln!(
        output,
        "step,time_star,kinetic_energy_star,resolved_dissipation_star,sgs_dissipation_star,sgs_to_resolved,max_eddy_to_molecular,mean_density"
    )?;
    arenas.scope(|arena| -> Result<(), Box<dyn std::error::Error>> {
        let cx = Cx::new(
            &gate,
            arena,
            fs_exec::StreamKey {
                seed: 0,
                kernel_id: 136,
                tile: 0,
                iteration: 0,
            },
            fs_exec::Budget::INFINITE,
            fs_exec::ExecMode::Deterministic,
        );
        for step in 0..=steps {
            // Both diagnostics describe the SAME current state. Preparation
            // is reused by collision; there is no duplicate SGS sweep.
            let report = grid.prepare_sgs(&cx, config, None, &mut workspace)?;
            let (mut energy, mut mass) = (0.0, 0.0);
            for cell in 0..cells {
                let (rho, u) = grid.macroscopic(cell)?;
                energy += 0.5 * rho * u.iter().map(|v| v * v).sum::<f64>();
                mass += rho;
            }
            let dissipation_scale = n as f64 / (cells as f64 * speed * speed * speed * k);
            writeln!(
                output,
                "{step},{:.12e},{:.12e},{:.12e},{:.12e},{:.12e},{:.12e},{:.12e}",
                step as f64 * k * speed / n as f64,
                energy / (cells as f64 * speed * speed),
                report.resolved_dissipation_lu * dissipation_scale,
                report.sgs_dissipation_lu * dissipation_scale,
                report.dissipation_ratio(),
                report.maximum_eddy_to_molecular,
                mass / cells as f64
            )?;
            if step < steps {
                grid.step_serial_with_collisions(&cx, workspace.collisions())?;
            }
        }
        Ok(())
    })?;
    output.flush()?;
    Ok(())
}
