//! Forced slit jet after de la Cuadra, Vergez & Fabre 2007 (arXiv 0811.0333,
//! Fig. 10). A 2-D jet leaves a long straight channel of height h into a
//! periodic strip and is forced by a uniform transverse acoustic pressure
//! gradient. Writes per-column jet centrelines as
//! `frankensim-jet-centreline-v1` (two estimators), the exit profile and a
//! run log with the cost.
//!
//! Inputs are dimensionless: Re = U h / nu, St = omega h / U (the paper's
//! definition, from its Fig. 10 titles), v_ac / U and lengths in h. The
//! declared h and U only label the output in SI. The lattice is D3Q19
//! (nz = 4, periodic, z-uniform) with a reduced-cumulant collision in
//! convective scaling: isothermal, weakly compressible NS at lattice Mach
//! sqrt(3) U_lat. Inflow is a uniform plug at the channel entry, or with
//! `plenum_h > 0` a uniform plug across a full-width cavity that feeds the
//! channel through a square-edged entry (the paper's Fig. 1 blows from a
//! cavity at p_f), adjusted
//! before the forcing until the centreline speed one h past the exit,
//! averaged over two periods of the channel's lattice quarter-wave mode,
//! equals U with that mode rung down, then frozen. `control=0` holds the
//! plug instead (a diagnostic, not a gate run). The outlet is the
//! regularized pressure face, not a nonreflecting one; the forcing vanishes
//! on the two layers it reads.
//!
//! Usage: `d3q19_forced_jet out=DIR [key=value ...]`; keys and defaults are
//! listed in `Config::parse` (the defaults are the paper's Re 500 case).

use std::f64::consts::PI;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::time::Instant;

use fs_exec::{CancelGate, TilePool};
use fs_lbm::d3q19::{
    BoundaryGrid3, BoundarySpec3, CollisionModel3, Face3, FaceBoundary3, InletFlux3, equilibrium3,
};
use fs_math::det;

struct Config {
    re: f64,
    st: f64,
    vu: f64,
    h_cells: usize,
    channel_h: f64,
    plenum_h: f64,
    entry_r_h: f64,
    free_h: f64,
    width_h: f64,
    u_lat: f64,
    h_m: f64,
    u_m_s: f64,
    periods: f64,
    samples_per_period: usize,
    settle_h: f64,
    column_step_h: f64,
    column_max_h: f64,
    control: bool,
    ring_max: f64,
    workers: usize,
    out: String,
}

impl Config {
    fn parse() -> Result<Config, String> {
        let mut c = Config {
            re: 500.0,
            st: 0.199958,
            vu: 0.005,
            h_cells: 20,
            channel_h: 27.0,
            plenum_h: 0.0,
            entry_r_h: 0.0,
            free_h: 30.0,
            width_h: 16.0,
            u_lat: 0.05,
            h_m: 1e-3,
            u_m_s: 7.487,
            periods: 8.0,
            samples_per_period: 50,
            settle_h: f64::NAN,
            column_step_h: 0.25,
            column_max_h: f64::NAN,
            control: true,
            ring_max: 0.01,
            workers: 8,
            out: String::new(),
        };
        for arg in std::env::args().skip(1) {
            let (key, value) = arg
                .split_once('=')
                .ok_or_else(|| format!("expected key=value, got {arg}"))?;
            let number = || value.parse::<f64>().map_err(|e| format!("{key}: {e}"));
            let count = || value.parse::<usize>().map_err(|e| format!("{key}: {e}"));
            match key {
                "re" => c.re = number()?,
                "st" => c.st = number()?,
                "vu" => c.vu = number()?,
                "h_cells" => c.h_cells = count()?,
                "channel_h" => c.channel_h = number()?,
                "plenum_h" => c.plenum_h = number()?,
                "entry_r_h" => c.entry_r_h = number()?,
                "free_h" => c.free_h = number()?,
                "width_h" => c.width_h = number()?,
                "u_lat" => c.u_lat = number()?,
                "h_m" => c.h_m = number()?,
                "u_m_s" => c.u_m_s = number()?,
                "periods" => c.periods = number()?,
                "samples_per_period" => c.samples_per_period = count()?,
                "settle_h" => c.settle_h = number()?,
                "column_step_h" => c.column_step_h = number()?,
                "column_max_h" => c.column_max_h = number()?,
                "control" => c.control = count()? != 0,
                "ring_max" => c.ring_max = number()?,
                "workers" => c.workers = count()?,
                "out" => c.out = value.to_string(),
                _ => return Err(format!("unknown key {key}")),
            }
        }
        if c.out.is_empty() {
            return Err("out=DIR is required".into());
        }
        if c.h_cells < 4 || !c.h_cells.is_multiple_of(2) {
            return Err("h_cells must be even and at least 4".into());
        }
        if c.entry_r_h > 0.0
            && (c.plenum_h <= 0.0
                || c.entry_r_h > c.channel_h
                || 1.0 + 2.0 * c.entry_r_h >= c.width_h)
        {
            return Err(
                "entry_r_h needs a cavity, a channel at least r long and walls wider than r".into(),
            );
        }
        // The starting vortex leaves the free region at about U/2.
        if c.settle_h.is_nan() {
            c.settle_h = 2.0 * c.free_h;
        }
        if c.column_max_h.is_nan() {
            c.column_max_h = (c.free_h - 3.0).min(25.0);
        }
        Ok(c)
    }
}

/// Streamwise velocity of one column at z = 0 (the flow is z-uniform).
fn column(grid: &BoundaryGrid3, x: usize) -> Vec<f64> {
    let ny = grid.dimensions()[1];
    (0..ny)
        .map(|y| {
            if grid.is_solid(x, y, 0) {
                0.0
            } else {
                grid.velocity(x, y, 0)[0]
            }
        })
        .collect()
}

/// Peak of a periodic column by a parabola through the largest sample:
/// (position of the peak in cell-centre coordinates, peak value, index).
fn peak(u: &[f64]) -> (f64, f64, usize) {
    let n = u.len();
    let j = (0..n)
        .max_by(|&a, &b| u[a].total_cmp(&u[b]))
        .expect("nonempty column");
    let (low, mid, high) = (u[(j + n - 1) % n], u[j], u[(j + 1) % n]);
    let curvature = low - 2.0 * mid + high;
    let shift = if curvature < 0.0 {
        0.5 * (low - high) / curvature
    } else {
        0.0
    };
    let value = mid - 0.25 * (low - high) * shift;
    (j as f64 + 0.5 + shift, value, j)
}

/// Centroid of u^2 over the contiguous region around the peak where
/// u >= 0.25 u_max, in cell-centre coordinates (periodic, unwrapped about j).
fn centroid(u: &[f64], j: usize, u_max: f64) -> f64 {
    let n = u.len() as isize;
    let at = |k: isize| u[((j as isize + k).rem_euclid(n)) as usize];
    let floor = 0.25 * u_max;
    let mut low = 0isize;
    while low > -n / 2 && at(low - 1) >= floor {
        low -= 1;
    }
    let mut high = 0isize;
    while high < n / 2 && at(high + 1) >= floor {
        high += 1;
    }
    let (mut moment, mut weight) = (0.0, 0.0);
    for k in low..=high {
        let w = at(k) * at(k);
        moment += k as f64 * w;
        weight += w;
    }
    j as f64 + 0.5 + moment / weight
}

/// FNV-1a over the bits of every fluid cell's populations: a replay receipt
/// for comparing builds and hosts, not a physical quantity.
fn state_digest(grid: &BoundaryGrid3) -> u64 {
    let [nx, ny, nz] = grid.dimensions();
    let mut hash = 0xcbf2_9ce4_8422_2325_u64;
    for z in 0..nz {
        for y in 0..ny {
            for x in 0..nx {
                if grid.is_solid(x, y, z) {
                    continue;
                }
                for value in grid.populations(x, y, z) {
                    for byte in value.to_bits().to_le_bytes() {
                        hash = (hash ^ u64::from(byte)).wrapping_mul(0x0100_0000_01b3);
                    }
                }
            }
        }
    }
    hash
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let c = Config::parse()?;
    let h = c.h_cells;
    let round4 = |n: usize| n.div_ceil(4) * 4;
    let channel = (c.channel_h * h as f64).round() as usize;
    let plenum = (c.plenum_h * h as f64).round() as usize;
    // First free column: the channel exit plane.
    let exit = plenum + channel;
    let nx = round4(exit + (c.free_h * h as f64).round() as usize);
    let ny = round4((c.width_h * h as f64).round() as usize);
    let nz = 4;
    let axis = (ny / 2) as f64;
    let nu = c.u_lat * h as f64 / c.re;
    let rate = 1.0 / (3.0 * nu + 0.5);
    let omega = c.st * c.u_lat / h as f64;
    let v_ac = c.vu * c.u_lat;
    let dx = c.h_m / h as f64;
    let dt = c.u_lat * dx / c.u_m_s;
    let f_hz = omega / (2.0 * PI * dt);
    let period_steps = 2.0 * PI / omega;

    std::fs::create_dir_all(&c.out)?;
    let mut log = BufWriter::new(File::create(format!("{}/run.log", c.out))?);
    macro_rules! say {
        ($($t:tt)*) => {{
            let line = format!($($t)*);
            println!("{line}");
            writeln!(log, "{line}")?;
            log.flush()?;
        }};
    }
    say!(
        "d3q19_forced_jet re {} st(omega h/U) {} v_ac/U {} h_cells {} channel_h {} plenum_h {} entry_r_h {} free_h {} width_h {} u_lat {} ring_max {} workers {}",
        c.re,
        c.st,
        c.vu,
        h,
        c.channel_h,
        c.plenum_h,
        c.entry_r_h,
        c.free_h,
        c.width_h,
        c.u_lat,
        c.ring_max,
        c.workers
    );
    say!(
        "lattice nx {nx} ny {ny} nz {nz} cells {} nu_lat {nu:.6e} rate2 {rate:.9} dx_m {dx:.6e} dt_s {dt:.6e} f_hz {f_hz:.6} period_steps {period_steps:.1}",
        nx * ny * nz
    );

    // u_in is the channel's bulk speed. A cavity carries the same flux over
    // its full width ny, so its plug moves at u_in h / ny.
    let mut u_in = 2.0 * c.u_lat / 3.0;
    let face_speed = |u: f64| {
        if plenum > 0 {
            u * h as f64 / ny as f64
        } else {
            u
        }
    };
    let periodic = FaceBoundary3::Periodic;
    let spec = BoundarySpec3::new([
        FaceBoundary3::Velocity {
            velocity: [face_speed(u_in), 0.0, 0.0],
        },
        FaceBoundary3::Pressure { density: 1.0 },
        periodic,
        periodic,
        periodic,
        periodic,
    ]);
    let model = CollisionModel3::ReducedCumulant {
        second_order_rate: rate,
        third_order_rate: 1.0,
        fourth_order_rate: 1.0,
    };
    let mut grid = BoundaryGrid3::with_collision_model(nx, ny, nz, model, [0.0; 3], spec);
    let half = h as f64 / 2.0;
    // The channel's entry edges are quarter circles of radius r tangent to
    // the cavity wall and the channel wall (r = 0: square edges); only the
    // sign is sampled.
    let r = c.entry_r_h * h as f64;
    grid.voxelize_sdf(|p| {
        if p[0] >= plenum as f64 && p[0] < exit as f64 {
            let d = plenum as f64 + r - p[0];
            let flare = if d > 0.0 {
                r - (r * r - d * d).sqrt()
            } else {
                0.0
            };
            half + flare - (p[1] - axis).abs()
        } else {
            1.0
        }
    });
    // Developed Poiseuille in the channel with its own pressure gradient
    // (rho c_s^2 falls by 12 nu U_mean / h^2 per cell, U_mean = 2U/3), the
    // cavity's plug at the entry density, rest elsewhere. An initial state
    // only, which the startup washes out; the gradient keeps the
    // inlet-outlet lattice resonator from being kicked.
    let gradient = 24.0 * nu * c.u_lat / (h as f64 * h as f64);
    for z in 0..nz {
        for y in 0..ny {
            for x in 0..exit {
                if !grid.is_solid(x, y, z) {
                    let x_c = (x.max(plenum) - plenum) as f64;
                    let rho = 1.0 + gradient * (channel as f64 - x_c - 0.5);
                    let u = if x < plenum {
                        face_speed(u_in)
                    } else {
                        let s = (y as f64 + 0.5 - axis) / half;
                        c.u_lat * (1.0 - s * s).max(0.0)
                    };
                    let f = equilibrium3(rho, [u, 0.0, 0.0]);
                    grid.set_populations(x, y, z, &f);
                }
            }
        }
    }
    // Uniform transverse pressure-gradient force density over the free
    // region, zero in the channel and on the two layers the outlet reads.
    // Held amplitude a = 2 sin(omega/2) v_ac with scale cos(omega (m + 1/2))
    // gives the free fluid a raw transverse velocity v_ac sin(omega m) after
    // m forced steps (Guo adds the held force each step).
    let amplitude = 2.0 * det::sin(0.5 * omega) * v_ac;
    grid.set_force_field(|[x, _, _]| {
        if x >= exit && x + 3 <= nx {
            [0.0, amplitude, 0.0]
        } else {
            [0.0; 3]
        }
    });
    grid.set_force_field_scale(0.0);

    let gate = CancelGate::new();
    let set_inflow = |grid: &mut BoundaryGrid3, u: f64| {
        grid.set_velocity_profile(Face3::XMin, InletFlux3::Unconstrained, &gate, |_| {
            Ok([face_speed(u), 0.0, 0.0])
        })
    };
    set_inflow(&mut grid, u_in)?;
    let fluid_cells = (0..nz)
        .flat_map(|z| (0..ny).flat_map(move |y| (0..nx).map(move |x| [x, y, z])))
        .filter(|&[x, y, z]| !grid.is_solid(x, y, z))
        .count();
    let cells = nx * ny * nz;

    // Exit centreline speed one h past the exit plane (between two cells).
    let probe = |grid: &BoundaryGrid3| {
        let a = peak(&column(grid, exit + h - 1)).1;
        let b = peak(&column(grid, exit + h)).1;
        0.5 * (a + b)
    };

    let step_columns: Vec<usize> = {
        let stride = ((c.column_step_h * h as f64).round() as usize).max(1);
        let last = exit + (c.column_max_h * h as f64).round() as usize;
        (exit..=last.min(nx - 3)).step_by(stride).collect()
    };
    let distances: Vec<f64> = step_columns
        .iter()
        .map(|&x| (x as f64 + 0.5 - exit as f64) * dx)
        .collect();

    let pool = TilePool::for_host(c.workers, 0);
    let wall = Instant::now();
    let result: Result<(), Box<dyn std::error::Error>> = pool.with_parked_crew_local(|parked| {
        let mut steps = 0usize;
        let step = |grid: &mut BoundaryGrid3, steps: &mut usize| -> Result<(), String> {
            grid.step_pooled(parked, &gate)
                .map_err(|e| format!("step {}: {e}", *steps + 1))?;
            *steps += 1;
            Ok(())
        };

        // Startup: settle the operating point with the forcing off. The
        // cavity and channel are a lattice quarter-wave duct (rigid inlet,
        // open exit) of period about 4 (plenum + channel + h) / c_s steps, and
        // a controller reading single samples pumps it. Each window averages
        // the probe over two such periods; the inflow then moves by half the
        // mean error.
        let duct = 4.0 * (exit + h) as f64 * 3f64.sqrt();
        let stride = ((2.0 * duct / 40.0).round() as usize).max(1);
        let window = 40 * stride;
        let settle = (c.settle_h * h as f64 / c.u_lat).round() as usize;
        let limit = (4 * settle).max(20 * window);
        let mut good = 0usize;
        // Every startup step's probe speed, to tell the channel's lattice
        // mode from jet flapping or outlet feedback when the ring persists.
        let mut trace = BufWriter::new(File::create(format!("{}/startup-probe.txt", c.out))?);
        writeln!(trace, "# step u_exit_over_U (probe one h past the exit), dt_s {dt:.6e}")?;
        loop {
            let (mut sum, mut low, mut high) = (0.0, f64::INFINITY, f64::NEG_INFINITY);
            for k in 1..=window {
                step(&mut grid, &mut steps)?;
                writeln!(trace, "{steps} {:.6}", probe(&grid) / c.u_lat)?;
                if k % stride == 0 {
                    let u = probe(&grid) / c.u_lat;
                    if !u.is_finite() {
                        return Err("exit speed is not finite".into());
                    }
                    sum += u;
                    low = low.min(u);
                    high = high.max(u);
                }
            }
            let mean = sum / 40.0;
            let ring = 0.5 * (high - low) / mean;
            // ring_max above the default 0.01 freezes a persistent oscillation
            // on its mean; the ring is still printed for every window.
            let settled = ring < c.ring_max && (!c.control || (mean - 1.0).abs() < 0.005);
            good = if settled { good + 1 } else { 0 };
            let elapsed = wall.elapsed().as_secs_f64();
            say!(
                "startup step {steps} t_ms {:.4} u_in/U {:.6} u_exit/U mean {mean:.6} ring {ring:.6} window {window} wall_s {elapsed:.1} mlups {:.3} digest {:016x}",
                steps as f64 * dt * 1e3,
                u_in / c.u_lat,
                (cells * steps) as f64 / elapsed / 1e6,
                state_digest(&grid)
            );
            trace.flush()?;
            if good >= 2 && steps >= settle {
                break;
            }
            if steps > limit {
                return Err(format!("operating point unsettled after {steps} steps").into());
            }
            if c.control {
                u_in *= (1.0 + 0.5 * (1.0 / mean - 1.0)).clamp(0.95, 1.05);
                set_inflow(&mut grid, u_in)?;
            }
        }
        if !c.control {
            say!("control off: the operating point is not U, so this is not a gate run");
        }
        let startup_steps = steps;
        say!(
            "operating point frozen: u_in/U {:.6} u_exit/U {:.6} after {startup_steps} steps ({:.4} ms) digest {:016x}",
            u_in / c.u_lat,
            probe(&grid) / c.u_lat,
            startup_steps as f64 * dt * 1e3,
            state_digest(&grid)
        );
        {
            let mut out = BufWriter::new(File::create(format!("{}/exit-profile.txt", c.out))?);
            writeln!(out, "# exit profile one h past the exit, u_x/U averaged over the two straddling columns")?;
            writeln!(out, "# y_m (from the channel axis) u_over_U")?;
            let (a, b) = (column(&grid, exit + h - 1), column(&grid, exit + h));
            for y in 0..ny {
                writeln!(
                    out,
                    "{:.6e} {:.9}",
                    (y as f64 + 0.5 - axis) * dx,
                    0.5 * (a[y] + b[y]) / c.u_lat
                )?;
            }
        }

        let header = |estimator: &str| -> Result<BufWriter<File>, std::io::Error> {
            let mut out = BufWriter::new(File::create(format!(
                "{}/centreline-{estimator}.txt",
                c.out
            ))?);
            writeln!(out, "# frankensim-jet-centreline-v1")?;
            writeln!(
                out,
                "# h_m {} u_b_m_s {} f_hz {f_hz:.9} v_ac_m_s {:.9e} estimator {estimator}",
                c.h_m,
                c.u_m_s,
                c.vu * c.u_m_s
            )?;
            writeln!(out, "# forcing v_ac(t) = v_ac_m_s * sin(2*pi*f*t)")?;
            write!(out, "# x_m")?;
            for d in &distances {
                write!(out, " {d:.6e}")?;
            }
            writeln!(out)?;
            Ok(out)
        };
        let mut argmax = header("argmax-ux-parabolic")?;
        let mut momentum = header("momentum-centroid")?;
        let record = |grid: &BoundaryGrid3,
                      m: usize,
                      argmax: &mut BufWriter<File>,
                      momentum: &mut BufWriter<File>|
         -> Result<(), std::io::Error> {
            let t = m as f64 * dt;
            write!(argmax, "{t:.9e}")?;
            write!(momentum, "{t:.9e}")?;
            for &x in &step_columns {
                let u = column(grid, x);
                let (y, value, j) = peak(&u);
                write!(argmax, " {:.6e}", (y - axis) * dx)?;
                write!(momentum, " {:.6e}", (centroid(&u, j, value) - axis) * dx)?;
            }
            writeln!(argmax)?;
            writeln!(momentum)?;
            Ok(())
        };

        // Forced phase.
        let forced = (c.periods * period_steps).round() as usize;
        let sample = ((period_steps / c.samples_per_period as f64).round() as usize).max(1);
        let report = (forced / 20).max(1);
        record(&grid, 0, &mut argmax, &mut momentum)?;
        for m in 0..forced {
            grid.set_force_field_scale(det::cos(omega * (m as f64 + 0.5)));
            step(&mut grid, &mut steps)?;
            if (m + 1) % sample == 0 {
                record(&grid, m + 1, &mut argmax, &mut momentum)?;
            }
            if (m + 1) % report == 0 {
                argmax.flush()?;
                momentum.flush()?;
                let elapsed = wall.elapsed().as_secs_f64();
                say!(
                    "forced step {} of {forced} t_forced_ms {:.4} u_exit/U {:.6} wall_s {elapsed:.1} mlups {:.3} digest {:016x}",
                    m + 1,
                    (m + 1) as f64 * dt * 1e3,
                    probe(&grid) / c.u_lat,
                    (cells * steps) as f64 / elapsed / 1e6,
                    state_digest(&grid)
                );
            }
        }
        argmax.flush()?;
        momentum.flush()?;
        let elapsed = wall.elapsed().as_secs_f64();
        let flow_ms = steps as f64 * dt * 1e3;
        say!(
            "cost cells {cells} fluid_cells {fluid_cells} dt_s {dt:.6e} steps {steps} flow_ms {flow_ms:.4} wall_s {elapsed:.1} wall_s_per_flow_ms {:.1} mlups {:.3} workers {} digest {:016x}",
            elapsed / flow_ms,
            (cells * steps) as f64 / elapsed / 1e6,
            c.workers,
            state_digest(&grid)
        );
        Ok(())
    });
    result
}
