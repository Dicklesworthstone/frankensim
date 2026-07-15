//! Bit-neutral SIMD dispatch for the frozen Duct axial-z BGK collision.
//!
//! SIMD lanes are independent cells inside one 4x4x4 SoA tile. Direction
//! reductions remain scalar and direction-ascending during admission, while
//! the equilibrium, Guo forcing, and relaxation expressions run lane-wise.
//! The general-force boundary-grid expression is intentionally not routed
//! here because it has a distinct frozen arithmetic surface.

use std::sync::OnceLock;

use super::{CollisionError3, CollisionModel3, E3, Q3, TILE_CELLS};
#[cfg(any(test, miri, not(target_arch = "aarch64")))]
use super::{W3, equilibrium3};

#[cfg(all(target_arch = "aarch64", not(miri)))]
mod neon;
mod stream;
#[cfg(all(target_arch = "x86_64", not(miri)))]
mod x86;

pub(super) use stream::stream_duct;
pub use stream::{D3q19StreamSimdTier, d3q19_stream_simd_tier};

pub(super) type TileInput<'a> = [&'a [f64; TILE_CELLS]; Q3];
pub(super) type TileOutput<'a> = [&'a mut [f64; TILE_CELLS]; Q3];

type TileKernel = for<'input, 'output> fn(
    &TileInput<'input>,
    &mut TileOutput<'output>,
    &MacroscopicTile,
    f64,
    f64,
) -> Result<(), CollisionError3>;

/// Effective implementation selected for the Duct axial-z BGK tile kernel.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum D3q19BgkSimdTier {
    /// Portable scalar twin.
    Scalar,
    /// AArch64 NEON, two independent cells per vector.
    Neon,
    /// x86-64 AVX2, four independent cells per vector.
    Avx2,
}

impl D3q19BgkSimdTier {
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

struct Dispatch {
    tier: D3q19BgkSimdTier,
    kernel: TileKernel,
}

pub(super) struct MacroscopicTile {
    pub(super) rho: [f64; TILE_CELLS],
    pub(super) velocity: [[f64; TILE_CELLS]; 3],
}

static DISPATCH: OnceLock<Dispatch> = OnceLock::new();

/// Tier selected once for the process, suitable for ledger/performance keys.
#[must_use]
pub fn d3q19_bgk_simd_tier() -> D3q19BgkSimdTier {
    dispatch().tier
}

pub(super) fn collide_bgk_axial_z_tile(
    input: &TileInput<'_>,
    output: &mut TileOutput<'_>,
    tau: f64,
    gz: f64,
) -> Result<(), CollisionError3> {
    CollisionModel3::Bgk { tau }.validate()?;
    if !gz.is_finite() {
        return Err(CollisionError3::NonFiniteForce { axis: 2, value: gz });
    }
    let macros = prepare_macroscopic_tile(input, gz)?;
    (dispatch().kernel)(input, output, &macros, tau, gz)
}

fn prepare_macroscopic_tile(
    input: &TileInput<'_>,
    gz: f64,
) -> Result<MacroscopicTile, CollisionError3> {
    let force = [0.0, 0.0, gz];
    let mut rho = [0.0; TILE_CELLS];
    let mut velocity = [[0.0; TILE_CELLS]; 3];
    for lane in 0..TILE_CELLS {
        let mut lane_rho = 0.0;
        let mut momentum = [0.0; 3];
        for direction in 0..Q3 {
            let value = input[direction][lane];
            if !value.is_finite() {
                return Err(CollisionError3::NonFinitePopulation { direction, value });
            }
            lane_rho += value;
            momentum[0] += f64::from(E3[direction].0) * value;
            momentum[1] += f64::from(E3[direction].1) * value;
            momentum[2] += f64::from(E3[direction].2) * value;
        }
        if !(lane_rho.is_finite() && lane_rho > 0.0) {
            return Err(CollisionError3::NonPositiveDensity { rho: lane_rho });
        }
        rho[lane] = lane_rho;
        for axis in 0..3 {
            let value = (momentum[axis] + 0.5 * force[axis]) / lane_rho;
            if !value.is_finite() {
                return Err(CollisionError3::NonFiniteVelocity { axis, value });
            }
            velocity[axis][lane] = value;
        }
    }
    Ok(MacroscopicTile { rho, velocity })
}

#[cfg(any(test, miri, not(target_arch = "aarch64")))]
pub(super) fn scalar_kernel(
    input: &TileInput<'_>,
    output: &mut TileOutput<'_>,
    macros: &MacroscopicTile,
    tau: f64,
    gz: f64,
) -> Result<(), CollisionError3> {
    let coefficient = 1.0 - 0.5 / tau;
    for lane in 0..TILE_CELLS {
        let velocity = [
            macros.velocity[0][lane],
            macros.velocity[1][lane],
            macros.velocity[2][lane],
        ];
        let equilibrium = equilibrium3(macros.rho[lane], velocity);
        for direction in 0..Q3 {
            let e = [
                f64::from(E3[direction].0),
                f64::from(E3[direction].1),
                f64::from(E3[direction].2),
            ];
            let eu = e[0] * velocity[0] + e[1] * velocity[1] + e[2] * velocity[2];
            let forcing =
                coefficient * W3[direction] * (3.0 * (e[2] - velocity[2]) + 9.0 * eu * e[2]) * gz;
            let value = input[direction][lane]
                + (equilibrium[direction] - input[direction][lane]) / tau
                + forcing;
            if !value.is_finite() {
                return Err(CollisionError3::NonFiniteOutput { direction, value });
            }
            output[direction][lane] = value;
        }
    }
    Ok(())
}

#[cfg(all(not(miri), any(target_arch = "aarch64", target_arch = "x86_64")))]
pub(super) fn validate_output(output: &TileOutput<'_>) -> Result<(), CollisionError3> {
    for lane in 0..TILE_CELLS {
        for direction in 0..Q3 {
            let value = output[direction][lane];
            if !value.is_finite() {
                return Err(CollisionError3::NonFiniteOutput { direction, value });
            }
        }
    }
    Ok(())
}

fn dispatch() -> &'static Dispatch {
    DISPATCH.get_or_init(build_dispatch)
}

#[cfg(miri)]
fn build_dispatch() -> Dispatch {
    Dispatch {
        tier: D3q19BgkSimdTier::Scalar,
        kernel: scalar_kernel,
    }
}

#[cfg(all(not(miri), target_arch = "aarch64"))]
fn build_dispatch() -> Dispatch {
    Dispatch {
        tier: D3q19BgkSimdTier::Neon,
        kernel: neon::selected_kernel,
    }
}

#[cfg(all(not(miri), target_arch = "x86_64"))]
fn build_dispatch() -> Dispatch {
    let (kernel, tier) = x86::select_kernel();
    Dispatch { tier, kernel }
}

#[cfg(all(not(miri), not(any(target_arch = "aarch64", target_arch = "x86_64"))))]
fn build_dispatch() -> Dispatch {
    Dispatch {
        tier: D3q19BgkSimdTier::Scalar,
        kernel: scalar_kernel,
    }
}

#[cfg(any(test, target_arch = "x86_64"))]
pub(super) const fn selected_x86_tier_for(
    avx2_available: bool,
    fma_available: bool,
) -> D3q19BgkSimdTier {
    if avx2_available && fma_available {
        D3q19BgkSimdTier::Avx2
    } else {
        D3q19BgkSimdTier::Scalar
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const RECEIPT_SCHEMA: &str = "frankensim-d3q19-bgk-simd-v1";

    fn first_divergence(
        actual: &[[f64; TILE_CELLS]; Q3],
        expected: &[[f64; TILE_CELLS]; Q3],
    ) -> Option<(usize, usize)> {
        (0..TILE_CELLS).find_map(|lane| {
            (0..Q3)
                .find(|&direction| {
                    actual[direction][lane].to_bits() != expected[direction][lane].to_bits()
                })
                .map(|direction| (lane, direction))
        })
    }

    fn json_string(value: &str) -> String {
        format!("\"{}\"", value.replace('\\', "\\\\").replace('"', "\\\""))
    }

    const fn collision_error_kind(error: &CollisionError3) -> &'static str {
        match error {
            CollisionError3::InvalidRelaxationTime { .. } => "invalid_relaxation_time",
            CollisionError3::InvalidMomentRelaxationRate { .. } => "invalid_moment_relaxation_rate",
            CollisionError3::NonFiniteForce { .. } => "non_finite_force",
            CollisionError3::NonFinitePopulation { .. } => "non_finite_population",
            CollisionError3::NonPositiveDensity { .. } => "non_positive_density",
            CollisionError3::NonFiniteVelocity { .. } => "non_finite_velocity",
            CollisionError3::CentralMomentForceUnsupported { .. } => {
                "central_moment_force_unsupported"
            }
            CollisionError3::ReducedCumulantForceUnsupported { .. } => {
                "reduced_cumulant_force_unsupported"
            }
            CollisionError3::SingularCentralMomentTransform { .. } => {
                "singular_central_moment_transform"
            }
            CollisionError3::NonFiniteOutput { .. } => "non_finite_output",
        }
    }

    #[allow(clippy::too_many_arguments)] // replay-complete retained receipt
    fn emit_divergence(
        case: &str,
        seed: u64,
        tau: f64,
        gz: f64,
        tier: &str,
        input: &[[f64; TILE_CELLS]; Q3],
        expected: &[[f64; TILE_CELLS]; Q3],
        actual: &[[f64; TILE_CELLS]; Q3],
        lane: usize,
        direction: usize,
        kernel_error: Option<&CollisionError3>,
    ) {
        let population_bits = (0..Q3)
            .map(|q| format!("0x{:016x}", input[q][lane].to_bits()))
            .collect::<Vec<_>>();
        let kernel_error_kind = kernel_error
            .map(collision_error_kind)
            .map_or_else(|| "null".to_owned(), json_string);
        let kernel_error = kernel_error
            .map(ToString::to_string)
            .map_or_else(|| "null".to_owned(), |error| json_string(&error));
        println!(
            "{{\"schema\":\"{RECEIPT_SCHEMA}\",\"suite\":\"fs-lbm/d3q19-bgk-simd\",\"case\":\"{case}\",\"semantics_version\":{},\"layout\":\"q-major-soa-64-v1\",\"seed\":\"0x{seed:016x}\",\"tau_bits\":\"0x{:016x}\",\"gz_bits\":\"0x{:016x}\",\"tier\":\"{tier}\",\"kernel_error_kind\":{kernel_error_kind},\"kernel_error\":{kernel_error},\"first_divergence\":{{\"lane\":{lane},\"direction\":{direction},\"population_bits\":{population_bits:?},\"expected_bits\":\"0x{:016x}\",\"actual_bits\":\"0x{:016x}\"}},\"verdict\":\"fail\"}}",
            super::super::D3Q19_BIT_SEMANTICS_VERSION,
            tau.to_bits(),
            gz.to_bits(),
            expected[direction][lane].to_bits(),
            actual[direction][lane].to_bits(),
        );
    }

    fn emit_kernel_error(
        case: &str,
        seed: u64,
        tau: f64,
        gz: f64,
        tier: &str,
        error: &CollisionError3,
    ) {
        let error_kind = json_string(collision_error_kind(error));
        let error = json_string(&error.to_string());
        println!(
            "{{\"schema\":\"{RECEIPT_SCHEMA}\",\"suite\":\"fs-lbm/d3q19-bgk-simd\",\"case\":\"{case}\",\"semantics_version\":{},\"layout\":\"q-major-soa-64-v1\",\"seed\":\"0x{seed:016x}\",\"tau_bits\":\"0x{:016x}\",\"gz_bits\":\"0x{:016x}\",\"tier\":\"{tier}\",\"kernel_error_kind\":{error_kind},\"kernel_error\":{error},\"first_divergence\":null,\"verdict\":\"fail\"}}",
            super::super::D3Q19_BIT_SEMANTICS_VERSION,
            tau.to_bits(),
            gz.to_bits(),
        );
    }

    fn emit_pass(case: &str, seed: u64, tau: f64, gz: f64, tier: &str) {
        println!(
            "{{\"schema\":\"{RECEIPT_SCHEMA}\",\"suite\":\"fs-lbm/d3q19-bgk-simd\",\"case\":\"{case}\",\"semantics_version\":{},\"layout\":\"q-major-soa-64-v1\",\"seed\":\"0x{seed:016x}\",\"tau_bits\":\"0x{:016x}\",\"gz_bits\":\"0x{:016x}\",\"cells\":{TILE_CELLS},\"tier\":\"{tier}\",\"kernel_error_kind\":null,\"kernel_error\":null,\"first_divergence\":null,\"verdict\":\"pass\"}}",
            super::super::D3Q19_BIT_SEMANTICS_VERSION,
            tau.to_bits(),
            gz.to_bits(),
        );
    }

    fn seeded_tile(seed: u64) -> [[f64; TILE_CELLS]; Q3] {
        let mut state = seed;
        let mut tile = [[0.0; TILE_CELLS]; Q3];
        for lane in 0..TILE_CELLS {
            state = state
                .wrapping_mul(6_364_136_223_846_793_005)
                .wrapping_add(1_442_695_040_888_963_407);
            let signed = ((state >> 11) as f64 / (1_u64 << 53) as f64) * 2.0 - 1.0;
            let rho = 1.0 + signed * 0.02;
            let velocity = [signed * 0.015, -signed * 0.01, signed * 0.02];
            let equilibrium = equilibrium3(rho, velocity);
            for direction in 0..Q3 {
                let perturbation = (direction as f64 - 9.0) * (lane as f64 - 31.5) * 1e-10;
                tile[direction][lane] = equilibrium[direction] + perturbation;
            }
        }
        tile
    }

    fn resting_tile() -> [[f64; TILE_CELLS]; Q3] {
        let equilibrium = equilibrium3(1.0, [0.0; 3]);
        core::array::from_fn(|direction| [equilibrium[direction]; TILE_CELLS])
    }

    #[test]
    fn scalar_tile_twin_matches_frozen_axial_cell_authority() {
        let case = "scalar-twin-vs-frozen-axial";
        let seed = 0x4be2_0001;
        let tau = 0.83;
        let gz = 2.5e-6;
        let input = seeded_tile(seed);
        let input_refs = input.each_ref();
        let macros = prepare_macroscopic_tile(&input_refs, gz).expect("admitted tile");
        let mut actual = [[0.0; TILE_CELLS]; Q3];
        scalar_kernel(&input_refs, &mut actual.each_mut(), &macros, tau, gz)
            .expect("admitted tile");

        let mut expected = [[0.0; TILE_CELLS]; Q3];
        for lane in 0..TILE_CELLS {
            let populations = core::array::from_fn(|direction| input[direction][lane]);
            let cell =
                super::super::collide_axial_z_cell3(populations, CollisionModel3::Bgk { tau }, gz)
                    .expect("admitted cell");
            for direction in 0..Q3 {
                expected[direction][lane] = cell[direction];
            }
        }
        if let Some((lane, direction)) = first_divergence(&actual, &expected) {
            emit_divergence(
                case,
                seed,
                tau,
                gz,
                "scalar-twin",
                &input,
                &expected,
                &actual,
                lane,
                direction,
                None,
            );
            panic!("D3Q19 BGK scalar twin diverged at lane {lane}, direction {direction}");
        }
        emit_pass(case, seed, tau, gz, "scalar-twin");
    }

    #[test]
    fn active_tile_kernel_is_bitwise_to_scalar_twin() {
        let cases = [
            ("rest-zero-force", 0x4be2_0000, 0.8, 0.0, true),
            ("near-tau-floor", 0x4be2_0013, 0.500_001, 1e-7, false),
            ("negative-axial-force", 0x4be2_0041, 1.17, -3e-6, false),
        ];
        for (case, seed, tau, gz, is_resting) in cases {
            let input = if is_resting {
                resting_tile()
            } else {
                seeded_tile(seed)
            };
            let input_refs = input.each_ref();
            let macros = prepare_macroscopic_tile(&input_refs, gz).expect("admitted tile");
            let mut expected = [[0.0; TILE_CELLS]; Q3];
            scalar_kernel(&input_refs, &mut expected.each_mut(), &macros, tau, gz)
                .expect("admitted tile");
            let mut actual = [[0.0; TILE_CELLS]; Q3];
            let result = collide_bgk_axial_z_tile(&input_refs, &mut actual.each_mut(), tau, gz);
            if let Some((lane, direction)) = first_divergence(&actual, &expected) {
                emit_divergence(
                    case,
                    seed,
                    tau,
                    gz,
                    d3q19_bgk_simd_tier().name(),
                    &input,
                    &expected,
                    &actual,
                    lane,
                    direction,
                    result.as_ref().err(),
                );
                panic!("D3Q19 BGK SIMD first divergence at lane {lane}, direction {direction}");
            }
            if let Err(error) = result {
                emit_kernel_error(case, seed, tau, gz, d3q19_bgk_simd_tier().name(), &error);
                panic!("D3Q19 BGK SIMD returned {error}");
            }
            emit_pass(case, seed, tau, gz, d3q19_bgk_simd_tier().name());
        }
    }

    #[test]
    fn dispatch_selector_is_fail_closed_and_operation_specific() {
        assert_eq!(selected_x86_tier_for(true, true), D3q19BgkSimdTier::Avx2);
        assert_eq!(selected_x86_tier_for(true, false), D3q19BgkSimdTier::Scalar);
        assert_eq!(selected_x86_tier_for(false, true), D3q19BgkSimdTier::Scalar);
        assert_eq!(
            selected_x86_tier_for(false, false),
            D3q19BgkSimdTier::Scalar
        );
        #[cfg(all(target_arch = "aarch64", not(miri)))]
        assert_eq!(d3q19_bgk_simd_tier(), D3q19BgkSimdTier::Neon);
        #[cfg(all(target_arch = "x86_64", not(miri)))]
        assert_eq!(
            d3q19_bgk_simd_tier(),
            selected_x86_tier_for(
                std::arch::is_x86_feature_detected!("avx2"),
                std::arch::is_x86_feature_detected!("fma")
            )
        );
        #[cfg(miri)]
        assert_eq!(d3q19_bgk_simd_tier(), D3q19BgkSimdTier::Scalar);
        assert!(std::ptr::eq(dispatch(), dispatch()));
    }
}
