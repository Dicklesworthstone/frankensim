//! Native Apple C ABI over FrankenSim's real bounded laboratory kernels.

use std::cell::{Cell, RefCell};
use std::panic::{AssertUnwindSafe, catch_unwind};

const SCHEMA_VERSION: f64 = 1.0;
const HEADER_LEN: usize = 6;
const MAX_RESULT_VALUES: usize = 2_000_000;

#[repr(u32)]
#[derive(Clone, Copy)]
enum Shape {
    Signal = 0,
    Grid = 1,
    GridFrames = 2,
    XyzPath = 3,
    Triangles = 4,
    Campaign = 5,
}

thread_local! {
    static RESULT: RefCell<Vec<f64>> = const { RefCell::new(Vec::new()) };
    static LAST_ERROR: Cell<i32> = const { Cell::new(0) };
}

fn packet(
    id: u32,
    shape: Shape,
    width: usize,
    height: usize,
    frames: usize,
    payload: Vec<f64>,
) -> Vec<f64> {
    if payload.len() > MAX_RESULT_VALUES {
        return Vec::new();
    }
    let mut out = Vec::with_capacity(HEADER_LEN + payload.len());
    out.extend_from_slice(&[
        SCHEMA_VERSION,
        id as f64,
        shape as u32 as f64,
        width as f64,
        height as f64,
        frames as f64,
    ]);
    out.extend(payload);
    out
}

fn grid_packet(id: u32, n: usize, frames: usize, payload: Vec<f64>) -> Vec<f64> {
    packet(
        id,
        if frames == 1 {
            Shape::Grid
        } else {
            Shape::GridFrames
        },
        n,
        n,
        frames,
        payload,
    )
}

/// Mirror `fs_wasm::wave2d_frames`' admitted FFT grid so the packet metadata
/// describes the payload the kernel actually returns. The wave kernel rounds
/// non-power-of-two requests upward; reporting the pre-admission request made
/// the 48-cell quality tier claim a 48x48 field while carrying 64x64 frames.
fn admitted_wave_grid(n_in: usize) -> usize {
    n_in.clamp(8, 128).next_power_of_two().min(128)
}

fn signal_packet(id: u32, payload: Vec<f64>) -> Vec<f64> {
    packet(id, Shape::Signal, payload.len(), 1, 1, payload)
}

fn campaign_packet(id: u32, payload: Vec<f64>) -> Vec<f64> {
    packet(id, Shape::Campaign, payload.len(), 1, 1, payload)
}

fn navier_speed_frames(payload: Vec<f64>) -> Vec<f64> {
    let Some((&grid, rest)) = payload.split_first() else {
        return Vec::new();
    };
    let Some((&frames, fields)) = rest.split_first() else {
        return Vec::new();
    };
    let grid = grid as usize;
    let frames = frames as usize;
    let field_len = grid.saturating_mul(grid);
    let frame_len = field_len.saturating_mul(2);
    if grid != 20 || fields.len() < frames.saturating_mul(frame_len) {
        return Vec::new();
    }
    let mut speed = Vec::with_capacity(frames.saturating_mul(field_len));
    for frame in fields.chunks_exact(frame_len).take(frames) {
        speed.extend_from_slice(&frame[..field_len]);
    }
    speed
}

fn run(id: u32, quality: f64, seed: u32) -> Option<Vec<f64>> {
    let q = quality.clamp(0.0, 1.0);
    let grid = if q < 0.34 {
        32
    } else if q < 0.75 {
        48
    } else {
        64
    };
    let frames = if q < 0.34 {
        8
    } else if q < 0.75 {
        14
    } else {
        20
    };
    let signal = if q < 0.34 {
        96
    } else if q < 0.75 {
        160
    } else {
        256
    };

    let result = match id {
        0 => grid_packet(id, grid, frames, fs_wasm::heat_frames(grid, frames, 5)),
        1 => signal_packet(
            id,
            fs_wasm::orr_sommerfeld_curve(1.02, 40, 4_000.0, 8_000.0, signal),
        ),
        2 => signal_packet(id, fs_wasm::chebyshev_fit(seed % 4, signal)),
        3 => signal_packet(id, fs_wasm::taylor_bound(0.25, 0.55, 12)),
        4 => signal_packet(id, fs_wasm::autodiff_derivatives(-3.0, 3.0, signal)),
        5 => signal_packet(id, fs_wasm::randomized_svd(36, 8, seed)),
        6 => signal_packet(id, fs_wasm::fft_power_spectrum(512, seed)),
        7 => signal_packet(id, fs_wasm::laplacian_modes(36, 5)),
        8 => signal_packet(id, fs_wasm::qmc_vs_mc(14, seed)),
        9 => signal_packet(id, fs_wasm::robust_hull(8)),
        10 => packet(
            id,
            Shape::GridFrames,
            52,
            30,
            frames,
            fs_wasm::topopt_frames(52, 30, frames, 0.42),
        ),
        11 => grid_packet(id, 40, 40, fs_wasm::sdf_volume(40, seed % 4, 0.35)),
        12 => packet(
            id,
            Shape::Triangles,
            0,
            0,
            1,
            fs_wasm::marching_cubes(30, seed % 4, 0.0),
        ),
        13 => packet(
            id,
            Shape::XyzPath,
            0,
            0,
            1,
            fs_wasm::lorenz_points(5_000, 0.006, 28.0),
        ),
        14 => {
            let wave_grid = admitted_wave_grid(grid);
            grid_packet(
                id,
                wave_grid,
                frames,
                fs_wasm::wave2d_frames(wave_grid, frames, 3),
            )
        }
        15 => grid_packet(id, grid, frames, fs_wasm::fluid_frames(grid, frames)),
        16 => grid_packet(
            id,
            grid,
            frames,
            fs_wasm::gray_scott_frames(grid, frames, 0.037, 0.060),
        ),
        17 => grid_packet(
            id,
            96,
            1,
            fs_wasm::mandelbrot_certified(96, 96, -0.55, 0.0, 1.5, 120),
        ),
        18 => {
            let mut values = fs_wasm::ga_motor_orbit(28, 80);
            if values.len() >= 2 {
                values.drain(0..2);
            }
            packet(id, Shape::XyzPath, 28, 1, 80, values)
        }
        19 => signal_packet(id, fs_wasm::symplectic_vs_euler(1_200, 0.015)),
        20 => signal_packet(id, fs_wasm::hodge_decomposition(seed % 3)),
        21 => {
            let values = navier_speed_frames(fs_wasm::navier_stokes_cavity(5, 5, 100.0, 2));
            packet(id, Shape::GridFrames, 20, 20, 5, values)
        }
        22 => signal_packet(id, fs_wasm::gp_regression(9, signal)),
        23 => signal_packet(id, fs_wasm::cmaes_trace(seed, 36)),
        24 => signal_packet(id, fs_wasm::optimal_transport(32, 0.05)),
        25 => signal_packet(id, fs_wasm::cyclic_symmetry(24, 0.8)),
        26 => signal_packet(id, fs_wasm::krylov_convergence(16, 300)),
        27 => signal_packet(id, fs_wasm::cutfem_quadtree(3, 6, 0.31)),
        28 => signal_packet(id, fs_wasm::ffd_deform(14, 4, 0.65, seed % 3)),
        29 => signal_packet(id, fs_wasm::betti_shapes(seed % 3)),
        30 => campaign_packet(id, fs_wasm::proofrobust(0.95, 0.35, 51)),
        31 => campaign_packet(id, fs_wasm::metamatcert(10, 8, 0.35)),
        32 => campaign_packet(id, fs_wasm::fluttercert(0.2, 2.8, 80)),
        33 => campaign_packet(id, fs_wasm::schedule_campaign(9.0, 0.84, 0.11)),
        34 => campaign_packet(id, fs_wasm::trusspath(4, 3, 1.0e-4)),
        35 => campaign_packet(id, fs_wasm::sensorforge(0.08, 8, 0.25)),
        36 => campaign_packet(id, fs_wasm::neuroshape(7.0, 2.3, 0.4)),
        37 => campaign_packet(id, fs_wasm::grammarforge(0.12, 0.03)),
        38 => campaign_packet(id, fs_wasm::anytimebo(30, 0.05, 0.05)),
        39 => campaign_packet(id, fs_wasm::flowcert(2_500, 0.08)),
        40 => campaign_packet(id, fs_wasm::run_ornithoid(seed)),
        41 => campaign_packet(id, fs_wasm::run_vessel(650)),
        42 => campaign_packet(id, fs_wasm::run_frame(seed)),
        43 => signal_packet(id, fs_wasm::run_instrument_reed(2_800.0 + quality * 4_000.0)),
        _ => return None,
    };
    Some(result)
}

#[unsafe(no_mangle)]
pub extern "C" fn frankensim_apple_schema_version() -> u32 {
    SCHEMA_VERSION as u32
}

#[unsafe(no_mangle)]
pub extern "C" fn frankensim_apple_run(id: u32, quality: f64, seed: u32) -> u64 {
    LAST_ERROR.with(|error| error.set(0));
    let attempted = catch_unwind(AssertUnwindSafe(|| run(id, quality, seed)));
    let (values, error) = match attempted {
        Ok(Some(values)) if values.len() >= HEADER_LEN => (values, 0),
        Ok(Some(_)) => (Vec::new(), 3),
        Ok(None) => (Vec::new(), 1),
        Err(_) => (Vec::new(), 2),
    };
    let len = values.len() as u64;
    RESULT.with(|result| *result.borrow_mut() = values);
    LAST_ERROR.with(|last| last.set(error));
    len
}

#[unsafe(no_mangle)]
pub extern "C" fn frankensim_apple_result_len() -> u64 {
    RESULT.with(|result| result.borrow().len() as u64)
}

#[unsafe(no_mangle)]
pub extern "C" fn frankensim_apple_result_value(index: u64) -> f64 {
    RESULT.with(|result| {
        result
            .borrow()
            .get(index as usize)
            .copied()
            .unwrap_or(f64::NAN)
    })
}

#[unsafe(no_mangle)]
pub extern "C" fn frankensim_apple_last_error() -> i32 {
    LAST_ERROR.with(Cell::get)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn packet_header_is_stable() {
        assert!(frankensim_apple_run(13, 0.0, 7) > HEADER_LEN as u64);
        assert_eq!(frankensim_apple_result_value(0), SCHEMA_VERSION);
        assert_eq!(frankensim_apple_result_value(1), 13.0);
        assert_eq!(
            frankensim_apple_result_value(2),
            Shape::XyzPath as u32 as f64
        );
        assert_eq!(frankensim_apple_last_error(), 0);
    }

    #[test]
    fn unknown_id_refuses_and_clears_previous_result() {
        assert!(frankensim_apple_run(0, 0.0, 1) > 0);
        assert_eq!(frankensim_apple_run(u32::MAX, 0.0, 1), 0);
        assert_eq!(frankensim_apple_result_len(), 0);
        assert_eq!(frankensim_apple_last_error(), 1);
    }

    #[test]
    fn structured_native_visuals_preserve_their_declared_shapes() {
        let volume = run(11, 0.0, 1).expect("SDF slice packet");
        assert_eq!(volume[2], Shape::GridFrames as u32 as f64);
        assert_eq!(&volume[3..6], &[40.0, 40.0, 40.0]);
        assert_eq!(volume.len(), HEADER_LEN + 40 * 40 * 40);

        let surface = run(12, 0.0, 1).expect("surface packet");
        assert_eq!(surface[2], Shape::Triangles as u32 as f64);
        let triangles = surface[HEADER_LEN] as usize;
        assert_eq!(surface.len(), HEADER_LEN + 1 + triangles * 18);

        let cavity = run(21, 0.0, 1).expect("cavity packet");
        assert_eq!(cavity[2], Shape::GridFrames as u32 as f64);
        assert_eq!(&cavity[3..6], &[20.0, 20.0, 5.0]);
        assert_eq!(cavity.len(), HEADER_LEN + 20 * 20 * 5);
    }

    #[test]
    fn spectral_wave_metadata_matches_the_fft_admitted_grid() {
        let medium = run(14, 0.55, 1).expect("medium spectral wave packet");
        assert_eq!(medium[2], Shape::GridFrames as u32 as f64);
        assert_eq!(&medium[3..6], &[64.0, 64.0, 14.0]);
        assert_eq!(medium.len(), HEADER_LEN + 64 * 64 * 14);

        let low = run(14, 0.2, 1).expect("low spectral wave packet");
        assert_eq!(&low[3..6], &[32.0, 32.0, 8.0]);
        assert_eq!(low.len(), HEADER_LEN + 32 * 32 * 8);
    }

    #[test]
    fn every_public_catalog_entry_returns_a_bounded_packet() {
        for id in 0..44 {
            let result = run(id, 0.0, 0x5EED).unwrap_or_else(|| panic!("catalog id {id}"));
            assert!(result.len() >= HEADER_LEN, "catalog id {id}");
            assert!(
                result.len() <= HEADER_LEN + MAX_RESULT_VALUES,
                "catalog id {id}"
            );
            assert_eq!(result[0], SCHEMA_VERSION, "catalog id {id}");
            assert_eq!(result[1], id as f64, "catalog id {id}");
            assert!(
                result[HEADER_LEN..].iter().any(|value| value.is_finite()),
                "catalog id {id}"
            );
        }
    }
}
