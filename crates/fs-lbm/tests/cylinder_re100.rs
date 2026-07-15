//! G2 circular-cylinder crossflow battery for the D2Q9 core.
//!
//! The release-scale flow run is ignored in the default fast suite and must be
//! invoked explicitly in release mode. The FFT estimator, blockage
//! extrapolation, and geometry contracts remain ordinary tests.

use fs_fft::RealFft;
use fs_lbm::core2::VelocityPressureX2;
use fs_lbm::{Cell, Grid, equilibrium, plan_scaling};
use fs_math::det;

const REYNOLDS: f64 = 100.0;
const DIAMETER: usize = 10;
const INLET_SPEED: f64 = 0.1;
const STREAMWISE_DIAMETERS: usize = 32;
const UPSTREAM_DIAMETERS: usize = 8;
const WARMUP_STEPS: usize = 8_192;
const SPECTRUM_SAMPLES: usize = 16_384;
const MAX_STROUHAL_BIN_WIDTH: f64 = 0.01;

#[derive(Debug, Clone, Copy, PartialEq)]
struct StrouhalReceipt {
    sample_count: usize,
    peak_bin: usize,
    frequency: f64,
    strouhal: f64,
    delta_strouhal: f64,
    peak_power: f64,
}

#[derive(Debug, Clone, Copy, PartialEq)]
enum SpectrumRefusal {
    InvalidLength,
    NonFiniteHistory,
    InvalidScale,
    NoSignal,
    Underresolved { delta_strouhal: f64 },
    NyquistPeak,
    InvalidPeak,
}

#[allow(clippy::cast_precision_loss)]
fn estimate_strouhal(
    history: &[f64],
    sample_dt: f64,
    diameter: f64,
    reference_velocity: f64,
) -> Result<StrouhalReceipt, SpectrumRefusal> {
    let n = history.len();
    if n < 8 || !n.is_power_of_two() {
        return Err(SpectrumRefusal::InvalidLength);
    }
    if history.iter().any(|sample| !sample.is_finite()) {
        return Err(SpectrumRefusal::NonFiniteHistory);
    }
    if !sample_dt.is_finite()
        || sample_dt <= 0.0
        || !diameter.is_finite()
        || diameter <= 0.0
        || !reference_velocity.is_finite()
        || reference_velocity <= 0.0
    {
        return Err(SpectrumRefusal::InvalidScale);
    }

    let delta_strouhal = diameter / (n as f64 * sample_dt * reference_velocity);
    if !delta_strouhal.is_finite() || delta_strouhal > MAX_STROUHAL_BIN_WIDTH {
        return Err(SpectrumRefusal::Underresolved { delta_strouhal });
    }
    let mean = history.iter().sum::<f64>() / n as f64;
    let centered_energy = history
        .iter()
        .map(|sample| (sample - mean) * (sample - mean))
        .sum::<f64>();
    if !centered_energy.is_finite() || centered_energy <= 0.0 {
        return Err(SpectrumRefusal::NoSignal);
    }

    let denominator = (n - 1) as f64;
    let mut windowed = Vec::with_capacity(n);
    for (index, sample) in history.iter().enumerate() {
        let phase = 2.0 * core::f64::consts::PI * index as f64 / denominator;
        let hann = 0.5 - 0.5 * det::cos(phase);
        windowed.push((sample - mean) * hann);
    }
    let spectrum = RealFft::new(n).forward(&windowed);
    let peak_bin = (1..spectrum.len()).fold(1, |best, candidate| {
        if spectrum[candidate].norm_sq() > spectrum[best].norm_sq() {
            candidate
        } else {
            best
        }
    });
    if peak_bin + 1 == spectrum.len() {
        return Err(SpectrumRefusal::NyquistPeak);
    }
    let peak_power = spectrum[peak_bin].norm_sq();
    if !peak_power.is_finite() || peak_power <= 0.0 {
        return Err(SpectrumRefusal::InvalidPeak);
    }
    let frequency = peak_bin as f64 / (n as f64 * sample_dt);
    let strouhal = frequency * diameter / reference_velocity;
    if !frequency.is_finite() || !strouhal.is_finite() || frequency <= 0.0 {
        return Err(SpectrumRefusal::InvalidPeak);
    }
    Ok(StrouhalReceipt {
        sample_count: n,
        peak_bin,
        frequency,
        strouhal,
        delta_strouhal,
        peak_power,
    })
}

#[allow(clippy::cast_precision_loss)]
fn mark_cylinder(grid: &mut Grid, center: [usize; 2], diameter: usize) -> Vec<bool> {
    assert!(diameter >= 4 && diameter.is_multiple_of(2));
    let mut mask = vec![false; grid.flags.len()];
    let diameter_sq = i64::try_from(
        diameter
            .checked_mul(diameter)
            .expect("fixture diameter square fits usize"),
    )
    .expect("fixture diameter square fits i64");
    let center_x = i64::try_from(center[0]).expect("fixture x center fits i64");
    let center_y = i64::try_from(center[1]).expect("fixture y center fits i64");
    for y in 0..grid.ny {
        for x in 0..grid.nx {
            // Cell centers have odd doubled coordinates; the cylinder center
            // lies on an even doubled coordinate. For an even diameter this
            // produces exactly `diameter` occupied rows without privileging a
            // center cell.
            let doubled_x = 2 * i64::try_from(x).expect("fixture x fits i64") + 1;
            let doubled_y = 2 * i64::try_from(y).expect("fixture y fits i64") + 1;
            let dx = doubled_x - 2 * center_x;
            let dy = doubled_y - 2 * center_y;
            if dx * dx + dy * dy <= diameter_sq {
                let index = grid.idx(x, y);
                grid.flags[index] = Cell::Wall;
                mask[index] = true;
            }
        }
    }
    mask
}

#[derive(Debug, Clone, Copy)]
struct CylinderRun {
    span_diameters: usize,
    blockage: f64,
    tau: f64,
    mean_reference_density: f64,
    mean_cd: f64,
    strouhal: StrouhalReceipt,
    measured_links: usize,
}

#[allow(clippy::cast_precision_loss)]
fn run_cylinder(span_diameters: usize) -> CylinderRun {
    assert!(span_diameters >= 12 && span_diameters.is_multiple_of(2));
    let ny = span_diameters * DIAMETER;
    let nx = STREAMWISE_DIAMETERS * DIAMETER + 1;
    let center = [UPSTREAM_DIAMETERS * DIAMETER, ny / 2];
    let scaling = plan_scaling(REYNOLDS, DIAMETER as f64, INLET_SPEED);
    assert!(
        scaling.stable,
        "the declared Re=100 scaling must be admitted"
    );
    let mut grid = Grid::uniform(nx, ny, scaling.tau);
    grid.periodic_x = false;
    grid.f.fill(equilibrium(1.0, INLET_SPEED, 0.0));
    let measured_walls = mark_cylinder(&mut grid, center, DIAMETER);
    let seed = grid.idx(center[0] + DIAMETER, center[1] + DIAMETER / 2 + 1);
    assert_eq!(grid.flags[seed], Cell::Fluid);
    grid.f[seed] = equilibrium(1.0, INLET_SPEED, 1e-4);

    let boundary = VelocityPressureX2::new([INLET_SPEED, 0.0], 1.0);
    let mut scratch = Vec::new();
    for _ in 0..WARMUP_STEPS {
        let _ = grid.step_velocity_pressure_x_with_wall_momentum(
            &mut scratch,
            boundary,
            &measured_walls,
        );
    }

    let mut lift = Vec::with_capacity(SPECTRUM_SAMPLES);
    let mut drag_sum = 0.0;
    let mut inlet_density_sum = 0.0;
    let mut measured_links = None;
    for _ in 0..SPECTRUM_SAMPLES {
        let receipt = grid.step_velocity_pressure_x_with_wall_momentum(
            &mut scratch,
            boundary,
            &measured_walls,
        );
        assert!(receipt.wall_impulse.into_iter().all(f64::is_finite));
        if let Some(expected) = measured_links {
            assert_eq!(receipt.measured_links, expected);
        } else {
            assert!(receipt.measured_links > 0);
            measured_links = Some(receipt.measured_links);
        }
        drag_sum += receipt.wall_impulse[0];
        lift.push(receipt.wall_impulse[1]);
        inlet_density_sum += (0..grid.ny)
            .map(|y| grid.moments(grid.idx(0, y)).rho)
            .sum::<f64>()
            / grid.ny as f64;
    }
    let mean_drag = drag_sum / SPECTRUM_SAMPLES as f64;
    let mean_reference_density = inlet_density_sum / SPECTRUM_SAMPLES as f64;
    let mean_cd =
        2.0 * mean_drag / (mean_reference_density * INLET_SPEED * INLET_SPEED * DIAMETER as f64);
    let strouhal = estimate_strouhal(&lift, 1.0, DIAMETER as f64, INLET_SPEED)
        .expect("the release lift history must resolve a finite shedding peak");
    CylinderRun {
        span_diameters,
        blockage: 1.0 / span_diameters as f64,
        tau: scaling.tau,
        mean_reference_density,
        mean_cd,
        strouhal,
        measured_links: measured_links.expect("cylinder has measured links"),
    }
}

fn zero_blockage_linear(
    narrow_blockage: f64,
    narrow_value: f64,
    wide_blockage: f64,
    wide_value: f64,
) -> f64 {
    assert!(
        narrow_blockage.is_finite()
            && wide_blockage.is_finite()
            && narrow_blockage > wide_blockage
            && wide_blockage > 0.0
            && narrow_value.is_finite()
            && wide_value.is_finite()
    );
    (narrow_blockage * wide_value - wide_blockage * narrow_value)
        / (narrow_blockage - wide_blockage)
}

#[test]
fn lift_fft_estimator_resolves_an_exact_bin_and_refuses_bad_evidence() {
    const N: usize = 1_024;
    const BIN: usize = 16;
    #[allow(clippy::cast_precision_loss)]
    let signal: Vec<f64> = (0..N)
        .map(|index| det::sin(2.0 * core::f64::consts::PI * BIN as f64 * index as f64 / N as f64))
        .collect();
    let receipt = estimate_strouhal(&signal, 0.5, 2.0, 0.5).expect("resolved exact-bin signal");
    assert_eq!(receipt.sample_count, N);
    assert_eq!(receipt.peak_bin, BIN);
    assert!((receipt.frequency - 0.031_25).abs() < 1e-15);
    assert!((receipt.strouhal - 0.125).abs() < 1e-15);
    assert!((receipt.delta_strouhal - 0.007_812_5).abs() < 1e-15);
    assert!(receipt.peak_power.is_finite() && receipt.peak_power > 0.0);

    let shifted_scaled: Vec<f64> = signal.iter().map(|sample| 4.0 + 3.0 * sample).collect();
    let transformed =
        estimate_strouhal(&shifted_scaled, 0.5, 2.0, 0.5).expect("affine signal transform");
    assert_eq!(transformed.peak_bin, receipt.peak_bin);
    assert_eq!(transformed.frequency.to_bits(), receipt.frequency.to_bits());
    assert_eq!(transformed.strouhal.to_bits(), receipt.strouhal.to_bits());
    assert_eq!(
        estimate_strouhal(&signal, 0.5, 2.0, 0.5),
        estimate_strouhal(&signal, 0.5, 2.0, 0.5)
    );

    assert_eq!(
        estimate_strouhal(&[0.0; N], 0.5, 2.0, 0.5),
        Err(SpectrumRefusal::NoSignal)
    );
    let mut non_finite = signal.clone();
    non_finite[7] = f64::NAN;
    assert_eq!(
        estimate_strouhal(&non_finite, 0.5, 2.0, 0.5),
        Err(SpectrumRefusal::NonFiniteHistory)
    );
    assert!(matches!(
        estimate_strouhal(&signal, 1.0, 4.0, 0.1),
        Err(SpectrumRefusal::Underresolved { .. })
    ));
    assert_eq!(
        estimate_strouhal(&signal, 0.0, 2.0, 0.5),
        Err(SpectrumRefusal::InvalidScale)
    );
    let nyquist: Vec<f64> = (0..N)
        .map(|index| if index.is_multiple_of(2) { 1.0 } else { -1.0 })
        .collect();
    assert_eq!(
        estimate_strouhal(&nyquist, 0.5, 2.0, 0.5),
        Err(SpectrumRefusal::NyquistPeak)
    );
    assert_eq!(
        estimate_strouhal(&signal[..1_000], 0.5, 2.0, 0.5),
        Err(SpectrumRefusal::InvalidLength)
    );
}

#[test]
fn cylinder_mask_and_linear_blockage_extrapolation_are_pinned() {
    let mut grid = Grid::uniform(40, 30, 0.8);
    let mask = mark_cylinder(&mut grid, [15, 15], DIAMETER);
    let occupied_rows = (0..grid.ny)
        .filter(|&y| (0..grid.nx).any(|x| mask[grid.idx(x, y)]))
        .count();
    assert_eq!(occupied_rows, DIAMETER);
    assert_eq!(
        mask.iter().filter(|measured| **measured).count(),
        grid.flags
            .iter()
            .filter(|flag| **flag == Cell::Wall)
            .count()
    );

    let q0 = 1.325;
    let slope = 0.6;
    let narrow_blockage = 1.0 / 12.0;
    let wide_blockage = 1.0 / 16.0;
    let extrapolated = zero_blockage_linear(
        narrow_blockage,
        q0 + slope * narrow_blockage,
        wide_blockage,
        q0 + slope * wide_blockage,
    );
    assert!((extrapolated - q0).abs() < 1e-14);
}

/// Release invocation:
/// `cargo test -p fs-lbm --release --test cylinder_re100 -- --ignored --nocapture`
#[test]
#[ignore = "release-scale G2: two Re=100 domains, force history, and lift FFT"]
fn lbm_109_cylinder_re100_cd_and_strouhal() {
    // Primary references:
    // - Roshko, NACA TR-1191: St=0.212(1-21.2/Re), so St(100)=0.167056;
    //   the report says the fit is about 1% accurate and its tunnel data were
    //   blockage-corrected. https://ntrs.nasa.gov/citations/19930092207
    // - Posdziech & Grundmann (2007), doi:10.1016/j.jfluidstructs.2006.09.004:
    //   domain-asymptotic 2D Re=100 mean Cd=1.325, St=0.1644.
    // - Behr et al. (1995), doi:10.1016/0045-7825(94)00736-7: at least 8D
    //   cylinder-to-lateral-boundary distance for ~1% domain stability.
    // - Maskell, ARC R&M 3400: closed-tunnel wake blockage is leading-order
    //   linear in frontal-area ratio but its coefficient is geometry/base-
    //   pressure dependent. We therefore estimate the zero-blockage intercept
    //   from two widths instead of assuming a universal coefficient.
    let narrow = run_cylinder(12);
    let wide = run_cylinder(16);
    let corrected_cd =
        zero_blockage_linear(narrow.blockage, narrow.mean_cd, wide.blockage, wide.mean_cd);
    let corrected_st = zero_blockage_linear(
        narrow.blockage,
        narrow.strouhal.strouhal,
        wide.blockage,
        wide.strouhal.strouhal,
    );
    let pass = (1.25..=1.45).contains(&corrected_cd)
        && (0.155..=0.175).contains(&corrected_st)
        && corrected_cd.is_finite()
        && corrected_st.is_finite();
    println!(
        "{{\"test\":\"lbm-109\",\"model\":\"D2Q9-BGK-stair-step-cylinder-periodic-y\",\
         \"reynolds\":{REYNOLDS},\"diameter_lu\":{DIAMETER},\"u_lattice\":{INLET_SPEED},\
         \"warmup_steps\":{WARMUP_STEPS},\"spectrum_samples\":{SPECTRUM_SAMPLES},\
         \"narrow\":{{\"span_diameters\":{},\"blockage\":{:.6},\"tau\":{:.6},\"mean_reference_density\":{:.9},\"mean_cd\":{:.6},\"st\":{:.6},\"delta_st\":{:.6},\"links\":{}}},\
         \"wide\":{{\"span_diameters\":{},\"blockage\":{:.6},\"tau\":{:.6},\"mean_reference_density\":{:.9},\"mean_cd\":{:.6},\"st\":{:.6},\"delta_st\":{:.6},\"links\":{}}},\
         \"blockage_model\":\"two-width linear-beta intercept\",\"corrected_cd\":{corrected_cd:.6},\"corrected_st\":{corrected_st:.6},\
         \"cd_envelope\":[1.25,1.45],\"st_envelope\":[0.155,0.175],\"verdict\":\"{}\"}}",
        narrow.span_diameters,
        narrow.blockage,
        narrow.tau,
        narrow.mean_reference_density,
        narrow.mean_cd,
        narrow.strouhal.strouhal,
        narrow.strouhal.delta_strouhal,
        narrow.measured_links,
        wide.span_diameters,
        wide.blockage,
        wide.tau,
        wide.mean_reference_density,
        wide.mean_cd,
        wide.strouhal.strouhal,
        wide.strouhal.delta_strouhal,
        wide.measured_links,
        if pass { "pass" } else { "fail" }
    );
    assert!(pass, "Re=100 cylinder envelope failed");
}
