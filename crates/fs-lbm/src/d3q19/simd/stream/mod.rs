//! Bit-neutral SIMD dispatch for the Duct pull-stream move.
//!
//! The kernel processes one four-cell x-row at a time. Directions without an
//! x shift are contiguous moves; x-moving directions shift a contiguous row
//! and insert the one neighboring-tile or wall-bounce edge value. No floating
//! arithmetic is performed, so every population bit is preserved exactly.

use std::sync::OnceLock;

#[cfg(test)]
use super::super::D3Q19_BIT_SEMANTICS_VERSION;
#[cfg(any(test, miri, not(target_arch = "aarch64")))]
use super::super::{E3, OPP3};
use super::super::{Q3, TILE, Tile};

#[cfg(all(target_arch = "aarch64", not(miri)))]
mod neon;
#[cfg(all(target_arch = "x86_64", not(miri)))]
mod x86;

pub(super) type StreamKernel = fn(&[Vec<Tile>; Q3], &mut [Vec<Tile>; Q3], usize, usize, usize);

/// Effective implementation selected for the Duct pull-stream move.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum D3q19StreamSimdTier {
    /// Portable scalar twin.
    Scalar,
    /// AArch64 NEON, two cells per vector move.
    Neon,
    /// x86-64 AVX2, four cells per vector move.
    Avx2,
}

impl D3q19StreamSimdTier {
    /// Stable lowercase name for receipts, ledger rows, and tune keys.
    #[must_use]
    pub const fn name(self) -> &'static str {
        match self {
            Self::Scalar => "scalar",
            Self::Neon => "neon",
            Self::Avx2 => "avx2",
        }
    }
}

struct StreamDispatch {
    tier: D3q19StreamSimdTier,
    kernel: StreamKernel,
}

static STREAM_DISPATCH: OnceLock<StreamDispatch> = OnceLock::new();

/// Tier selected once for the Duct stream operation.
#[must_use]
pub fn d3q19_stream_simd_tier() -> D3q19StreamSimdTier {
    dispatch().tier
}

pub(super) fn stream_duct(
    post: &[Vec<Tile>; Q3],
    populations: &mut [Vec<Tile>; Q3],
    nx: usize,
    ny: usize,
    nz: usize,
) {
    validate_layout(post, populations, nx, ny, nz);
    (dispatch().kernel)(post, populations, nx, ny, nz);
}

fn validate_layout(
    post: &[Vec<Tile>; Q3],
    populations: &[Vec<Tile>; Q3],
    nx: usize,
    ny: usize,
    nz: usize,
) {
    assert!(
        nx > 0
            && ny > 0
            && nz > 0
            && nx.is_multiple_of(TILE)
            && ny.is_multiple_of(TILE)
            && nz.is_multiple_of(TILE),
        "D3Q19 stream dimensions must be positive multiples of {TILE}"
    );
    let tiles = (nx / TILE)
        .checked_mul(ny / TILE)
        .and_then(|count| count.checked_mul(nz / TILE))
        .expect("D3Q19 stream tile count overflow");
    assert!(
        post.iter().all(|field| field.len() == tiles)
            && populations.iter().all(|field| field.len() == tiles),
        "D3Q19 stream field tile-count mismatch"
    );
}

#[inline]
pub(super) const fn tile_index(ntx: usize, nty: usize, tx: usize, ty: usize, tz: usize) -> usize {
    (tz * nty + ty) * ntx + tx
}

#[inline]
pub(super) const fn row_lane(y: usize, z: usize) -> usize {
    (z % TILE * TILE + y % TILE) * TILE
}

#[inline]
pub(super) const fn pull_source(
    coordinate: usize,
    extent: usize,
    velocity: i32,
    periodic: bool,
) -> Option<usize> {
    match velocity {
        1 if coordinate == 0 => {
            if periodic {
                Some(extent - 1)
            } else {
                None
            }
        }
        1 => Some(coordinate - 1),
        -1 if coordinate + 1 == extent => {
            if periodic {
                Some(0)
            } else {
                None
            }
        }
        -1 => Some(coordinate + 1),
        0 => Some(coordinate),
        _ => unreachable!(),
    }
}

#[cfg(any(test, miri, not(target_arch = "aarch64")))]
pub(super) fn scalar_kernel(
    post: &[Vec<Tile>; Q3],
    populations: &mut [Vec<Tile>; Q3],
    nx: usize,
    ny: usize,
    nz: usize,
) {
    let (ntx, nty) = (nx / TILE, ny / TILE);
    for z in 0..nz {
        for y in 0..ny {
            for x in 0..nx {
                let dtile = tile_index(ntx, nty, x / TILE, y / TILE, z / TILE);
                let dlane = row_lane(y, z) + x % TILE;
                for direction in 0..Q3 {
                    let (ex, ey, ez) = E3[direction];
                    let sx = pull_source(x, nx, ex, false);
                    let sy = pull_source(y, ny, ey, false);
                    let sz = pull_source(z, nz, ez, true).expect("periodic z source");
                    populations[direction][dtile].0[dlane] = match (sx, sy) {
                        (Some(sx), Some(sy)) => {
                            let stile = tile_index(ntx, nty, sx / TILE, sy / TILE, sz / TILE);
                            let slane = row_lane(sy, sz) + sx % TILE;
                            post[direction][stile].0[slane]
                        }
                        _ => post[OPP3[direction]][dtile].0[dlane],
                    };
                }
            }
        }
    }
}

fn dispatch() -> &'static StreamDispatch {
    STREAM_DISPATCH.get_or_init(build_dispatch)
}

#[cfg(miri)]
fn build_dispatch() -> StreamDispatch {
    StreamDispatch {
        tier: D3q19StreamSimdTier::Scalar,
        kernel: scalar_kernel,
    }
}

#[cfg(all(not(miri), target_arch = "aarch64"))]
fn build_dispatch() -> StreamDispatch {
    StreamDispatch {
        tier: D3q19StreamSimdTier::Neon,
        kernel: neon::selected_kernel,
    }
}

#[cfg(all(not(miri), target_arch = "x86_64"))]
fn build_dispatch() -> StreamDispatch {
    let (kernel, tier) = x86::select_kernel();
    StreamDispatch { tier, kernel }
}

#[cfg(all(not(miri), not(any(target_arch = "aarch64", target_arch = "x86_64"))))]
fn build_dispatch() -> StreamDispatch {
    StreamDispatch {
        tier: D3q19StreamSimdTier::Scalar,
        kernel: scalar_kernel,
    }
}

#[cfg(any(test, target_arch = "x86_64"))]
pub(super) const fn selected_x86_tier_for(avx2_available: bool) -> D3q19StreamSimdTier {
    if avx2_available {
        D3q19StreamSimdTier::Avx2
    } else {
        D3q19StreamSimdTier::Scalar
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const RECEIPT_SCHEMA: &str = "frankensim-d3q19-stream-simd-v1";

    fn address(nx: usize, ny: usize, x: usize, y: usize, z: usize) -> (usize, usize) {
        let tile = tile_index(nx / TILE, ny / TILE, x / TILE, y / TILE, z / TILE);
        (tile, row_lane(y, z) + x % TILE)
    }

    fn seeded_fields(nx: usize, ny: usize, nz: usize, seed: u64) -> [Vec<Tile>; Q3] {
        let tiles = (nx / TILE) * (ny / TILE) * (nz / TILE);
        let mut fields = core::array::from_fn(|_| vec![Tile::filled(0.0); tiles]);
        let cells = nx * ny * nz;
        assert!((cells as u64) * (Q3 as u64) < 1_u64 << 32);
        for (direction, field) in fields.iter_mut().enumerate() {
            for z in 0..nz {
                for y in 0..ny {
                    for x in 0..nx {
                        let identity = direction * cells + (z * ny + y) * nx + x;
                        let mantissa = ((seed & 0x000f_ffff) << 32) | identity as u64;
                        let (tile, lane) = address(nx, ny, x, y, z);
                        field[tile].0[lane] = f64::from_bits(0x3ff0_0000_0000_0000 | mantissa);
                    }
                }
            }
        }
        fields
    }

    fn blank_fields(tiles: usize) -> [Vec<Tile>; Q3] {
        core::array::from_fn(|_| vec![Tile::filled(f64::from_bits(0x7ff8_dead_beef_0001)); tiles])
    }

    fn first_divergence(
        actual: &[Vec<Tile>; Q3],
        expected: &[Vec<Tile>; Q3],
        nx: usize,
        ny: usize,
        nz: usize,
    ) -> Option<(usize, usize, usize, usize, usize, usize)> {
        (0..nz).find_map(|z| {
            (0..ny).find_map(|y| {
                (0..nx).find_map(|x| {
                    let (tile, lane) = address(nx, ny, x, y, z);
                    (0..Q3)
                        .find(|&direction| {
                            actual[direction][tile].0[lane].to_bits()
                                != expected[direction][tile].0[lane].to_bits()
                        })
                        .map(|direction| (x, y, z, direction, tile, lane))
                })
            })
        })
    }

    fn assert_route(
        actual: &[Vec<Tile>; Q3],
        post: &[Vec<Tile>; Q3],
        nx: usize,
        ny: usize,
        destination: (usize, usize, usize, usize),
        source: (usize, usize, usize, usize),
    ) {
        let (dx, dy, dz, direction) = destination;
        let (sx, sy, sz, source_direction) = source;
        let (dtile, dlane) = address(nx, ny, dx, dy, dz);
        let (stile, slane) = address(nx, ny, sx, sy, sz);
        assert_eq!(
            actual[direction][dtile].0[dlane].to_bits(),
            post[source_direction][stile].0[slane].to_bits(),
            "D3Q19 stream route ({dx},{dy},{dz},q{direction}) must read ({sx},{sy},{sz},q{source_direction})"
        );
    }

    #[test]
    fn active_stream_is_bitwise_to_scalar_stencil() {
        let cases = [(4, 4, 4, 0x4be2_5101), (8, 12, 16, 0x4be2_5102)];
        for (nx, ny, nz, seed) in cases {
            let post = seeded_fields(nx, ny, nz, seed);
            let tiles = post[0].len();
            let mut expected = blank_fields(tiles);
            scalar_kernel(&post, &mut expected, nx, ny, nz);
            let mut actual = blank_fields(tiles);
            stream_duct(&post, &mut actual, nx, ny, nz);
            if let Some((x, y, z, direction, tile, lane)) =
                first_divergence(&actual, &expected, nx, ny, nz)
            {
                println!(
                    "{{\"schema\":\"{RECEIPT_SCHEMA}\",\"suite\":\"fs-lbm/d3q19-stream-simd\",\"semantics_version\":{D3Q19_BIT_SEMANTICS_VERSION},\"fixture_version\":1,\"layout\":\"q-major-tile-major-lane64-v1\",\"seed\":\"0x{seed:016x}\",\"dims\":[{nx},{ny},{nz}],\"tile_dims\":[{},{},{}],\"stream_tier\":\"{}\",\"first_divergence\":{{\"destination\":{{\"xyz\":[{x},{y},{z}],\"direction\":{direction},\"tile\":{tile},\"lane\":{lane}}},\"expected_bits\":\"0x{:016x}\",\"actual_bits\":\"0x{:016x}\"}},\"verdict\":\"fail\"}}",
                    nx / TILE,
                    ny / TILE,
                    nz / TILE,
                    d3q19_stream_simd_tier().name(),
                    expected[direction][tile].0[lane].to_bits(),
                    actual[direction][tile].0[lane].to_bits(),
                );
                panic!(
                    "D3Q19 stream SIMD diverged at direction {direction}, tile {tile}, lane {lane}"
                );
            }
            if (nx, ny, nz) == (8, 12, 16) {
                // Independent anchors cover all three tile crossings, periodic
                // z, the xy corner, and wall precedence over a simultaneous
                // z-wrap source.
                assert_route(&actual, &post, nx, ny, (4, 5, 6, 1), (3, 5, 6, 1));
                assert_route(&actual, &post, nx, ny, (5, 4, 6, 3), (5, 3, 6, 3));
                assert_route(&actual, &post, nx, ny, (5, 6, 4, 5), (5, 6, 3, 5));
                assert_route(&actual, &post, nx, ny, (5, 6, 0, 5), (5, 6, 15, 5));
                assert_route(&actual, &post, nx, ny, (0, 0, 7, 7), (0, 0, 7, 8));
                assert_route(&actual, &post, nx, ny, (0, 5, 0, 11), (0, 5, 0, 12));
            }
            println!(
                "{{\"schema\":\"{RECEIPT_SCHEMA}\",\"suite\":\"fs-lbm/d3q19-stream-simd\",\"semantics_version\":{D3Q19_BIT_SEMANTICS_VERSION},\"fixture_version\":1,\"layout\":\"q-major-tile-major-lane64-v1\",\"seed\":\"0x{seed:016x}\",\"dims\":[{nx},{ny},{nz}],\"tile_dims\":[{},{},{}],\"directions\":{Q3},\"cells\":{},\"populations\":{},\"route_anchors\":{},\"stream_tier\":\"{}\",\"first_divergence\":null,\"verdict\":\"pass\"}}",
                nx / TILE,
                ny / TILE,
                nz / TILE,
                nx * ny * nz,
                nx * ny * nz * Q3,
                if (nx, ny, nz) == (8, 12, 16) { 6 } else { 0 },
                d3q19_stream_simd_tier().name(),
            );
        }
    }

    #[test]
    fn stream_dispatch_is_one_shot_and_fail_closed() {
        assert_eq!(selected_x86_tier_for(true), D3q19StreamSimdTier::Avx2);
        assert_eq!(selected_x86_tier_for(false), D3q19StreamSimdTier::Scalar);
        #[cfg(all(target_arch = "aarch64", not(miri)))]
        assert_eq!(d3q19_stream_simd_tier(), D3q19StreamSimdTier::Neon);
        #[cfg(all(target_arch = "x86_64", not(miri)))]
        assert_eq!(
            d3q19_stream_simd_tier(),
            selected_x86_tier_for(std::arch::is_x86_feature_detected!("avx2"))
        );
        #[cfg(miri)]
        assert_eq!(d3q19_stream_simd_tier(), D3q19StreamSimdTier::Scalar);
        assert!(std::ptr::eq(dispatch(), dispatch()));
    }
}
