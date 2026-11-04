//! The forced jet's opt-in base development (fs_aeroac::forced_jet::BaseSpreading, MB 97): the
//! base widens as `b(x) = b0 beta(x)` with the centre speed `beta^(-1/2)`, quasi-parallel, held past
//! the station. The parallel base is the default and its arithmetic is unchanged; a spread base
//! whose width never leaves `b0` runs the spread arithmetic and returns the parallel jet bit for
//! bit; the laws are their closed forms; a laminar Bickley base at a very high Reynolds number
//! tends to the parallel jet; a wide base moves the response; the measured exit stays in phase at
//! the exit (the exit is the unspread profile); and every law, hold or band outside its domain is
//! refused.

use fs_aeroac::forced_jet::{
    BaseSpreading, ExitReceptivity, ForcedJet, ForcedJetResolution, JetStepResponse,
    bickley_step_response_with_base, bickley_step_response_with_exit,
};
use fs_math::det;

const WORKERS: usize = 4;
/// A station inside the coarse box, 6 exit half-widths from the exit. A numerical choice.
const STATION: f64 = 6.0;
/// Fletcher's `db/dx = 2 tan 5 deg` (JASA 60, 1976, eq. 5 and Table I, with `b = W`).
fn fletcher_slope() -> f64 {
    2.0 * det::tan(5.0_f64.to_radians())
}

/// Numerics, not physics: the coarse row of `tests/forced_jet_exit.rs` (`spacing 0.25`, `window
/// 40`, 256 cells, Courant 0.8, Re 600), with the band half-height given per case.
fn coarse(half_height: f64, reynolds: f64) -> ForcedJetResolution {
    ForcedJetResolution {
        spacing: 0.25,
        streamwise_cells: 256,
        half_height,
        window: 40.0,
        courant: 0.8,
        reynolds,
        workers: WORKERS,
    }
}

fn bits(step: &JetStepResponse) -> Vec<u64> {
    step.values().iter().map(|v| v.to_bits()).collect()
}

#[test]
fn parallel_is_the_default_and_unchanged() {
    assert_eq!(BaseSpreading::default(), BaseSpreading::Parallel);
    let res = coarse(8.0, 600.0);
    for exit in [
        ExitReceptivity::Kutta,
        ExitReceptivity::Measured {
            displacement_per_velocity: 4.0 / 3.0,
        },
    ] {
        let before = bickley_step_response_with_exit(STATION, &res, exit).unwrap();
        let named =
            bickley_step_response_with_base(STATION, &res, exit, BaseSpreading::Parallel).unwrap();
        assert_eq!(bits(&before), bits(&named), "{exit:?}");
        assert_eq!(named.spreading(), BaseSpreading::Parallel);
        assert_eq!(named.station_width_ratio(), 1.0);
    }
}

/// `beta = 1 + 1e-300 x` rounds to exactly 1 at every node, so the spread arrays hold the exit
/// profile at every column and the centre speed is exactly 1: the spread arithmetic (per-column
/// `U`, `U''` and centre speed) must then give the parallel jet's bits.
#[test]
fn a_spread_base_that_never_widens_is_the_parallel_jet_bit_for_bit() {
    let res = coarse(8.0, 600.0);
    let flat = BaseSpreading::Linear {
        width_slope: 1e-300,
    };
    assert_eq!(flat.width_ratio(res.length(), res.reynolds), 1.0);
    for exit in [
        ExitReceptivity::Kutta,
        ExitReceptivity::Measured {
            displacement_per_velocity: 4.0 / 3.0,
        },
    ] {
        let parallel = bickley_step_response_with_exit(STATION, &res, exit).unwrap();
        let spread = bickley_step_response_with_base(STATION, &res, exit, flat).unwrap();
        assert_eq!(bits(&parallel), bits(&spread), "{exit:?}");
    }
}

#[test]
fn the_laws_widths_are_their_closed_forms() {
    let re = 600.0;
    let bickley = BaseSpreading::LaminarBickley;
    assert_eq!(bickley.virtual_origin(re), Some(100.0));
    for x in [0.0, 1.0, 6.0, 30.0] {
        let beta = bickley.width_ratio(x, re);
        // beta^(3/2) = 1 + 6 x / Re.
        let closed = 1.0 + 6.0 * x / re;
        assert!(
            (det::pow(beta, 1.5) / closed - 1.0).abs() < 1e-14,
            "x {x}: beta {beta}"
        );
        // The centreline momentum balance u_c du_c/dx = -2 nu u_c / b^2 in exit units:
        // da/dx = -2 / (Re beta^2), a = beta^(-1/2), by a central difference.
        let d = 1e-4;
        let a = |x: f64| 1.0 / det::sqrt(bickley.width_ratio(x, re));
        let slope = (a(x + d) - a((x - d).max(0.0))) / (x + d - (x - d).max(0.0));
        let law = -2.0 / (re * beta * beta);
        assert!(
            (slope / law - 1.0).abs() < 1e-5,
            "x {x}: da/dx {slope} against {law}"
        );
    }
    let s = fletcher_slope();
    assert!((s - 0.174_977).abs() < 1e-6, "2 tan 5 deg = {s}");
    let linear = BaseSpreading::Linear { width_slope: s };
    assert_eq!(linear.width_ratio(4.0, re), 1.0 + s * 4.0);
    assert_eq!(linear.virtual_origin(re), Some(1.0 / s));
    assert_eq!(BaseSpreading::Parallel.virtual_origin(re), None);
}

/// The laminar base widens by `6 x / Re`: at Re 1e8 it is within 4e-7 of the exit's over the
/// station, so the response must be the parallel jet's to well within 1e-4 at w 0.2, yet not bit
/// for bit (the spread arithmetic ran). Tolerances are test choices, not physics.
#[test]
fn a_very_high_reynolds_bickley_base_tends_to_the_parallel_jet() {
    let res = coarse(8.0, 1e8);
    let parallel = bickley_step_response_with_exit(STATION, &res, ExitReceptivity::Kutta).unwrap();
    let spread = bickley_step_response_with_base(
        STATION,
        &res,
        ExitReceptivity::Kutta,
        BaseSpreading::LaminarBickley,
    )
    .unwrap();
    let (gp, gs) = (
        parallel.transfer(0.2).unwrap(),
        spread.transfer(0.2).unwrap(),
    );
    let rel = (gs - gp).abs() / gp.abs();
    println!(
        "Re 1e8: |G| parallel {:.9}, spread {:.9}, |dG|/|G| {rel:.3e}",
        gp.abs(),
        gs.abs()
    );
    assert!(rel < 1e-4, "{rel}");
    assert_ne!(bits(&parallel), bits(&spread));
}

/// Fletcher's slope doubles the width by the station (beta 2.05 at 6 b0): the response at w 0.3
/// must move by more than 5 % from the parallel jet's. A non-vacuity check, not a measured value:
/// it fails if the spread base never reaches the solve.
#[test]
fn a_wide_base_moves_the_response() {
    let wide = BaseSpreading::Linear {
        width_slope: fletcher_slope(),
    };
    let beta = wide.width_ratio(STATION, 600.0);
    let half = (4.0 * beta / 0.25).ceil() * 0.25;
    let res = coarse(half.max(8.0), 600.0);
    let parallel = bickley_step_response_with_exit(STATION, &res, ExitReceptivity::Kutta).unwrap();
    let spread =
        bickley_step_response_with_base(STATION, &res, ExitReceptivity::Kutta, wide).unwrap();
    assert!((spread.station_width_ratio() - beta).abs() < 1e-15);
    let (gp, gs) = (
        parallel.transfer(0.3).unwrap(),
        spread.transfer(0.3).unwrap(),
    );
    let rel = (gs - gp).abs() / gp.abs();
    println!(
        "Fletcher slope, beta(X) {beta:.4}: |G(0.3)| parallel {:.6}, spread {:.6}, |dG|/|G| {rel:.4}",
        gp.abs(),
        gs.abs()
    );
    assert!(rel > 0.05, "{rel}");
}

/// The exit is the unspread profile under every base, so the measured exit still issues the jet
/// displaced by `R v_ac` in phase: near the exit (0.5 b0) at w 0.01, `|G / R - 1| < 1 %` and
/// `|arg G| < 0.01 rad`, the tolerances of `tests/forced_jet_exit.rs`.
#[test]
fn the_measured_exit_stays_in_phase_at_the_exit_under_spreading() {
    let r = 4.0 / 3.0;
    let wide = BaseSpreading::Linear {
        width_slope: fletcher_slope(),
    };
    let res = coarse(8.0, 600.0);
    let g = bickley_step_response_with_base(
        0.5,
        &res,
        ExitReceptivity::Measured {
            displacement_per_velocity: r,
        },
        wide,
    )
    .unwrap()
    .transfer(0.01)
    .unwrap();
    let (mag, arg) = (g.abs(), det::atan2(g.im, g.re));
    println!(
        "near-exit G(0.01) under Fletcher spreading = {mag:.6} at {arg:+.6} rad against R {r:.6}"
    );
    assert!((mag / r - 1.0).abs() < 0.01, "{mag}");
    assert!(arg.abs() < 0.01, "{arg}");
}

#[test]
fn laws_holds_and_bands_outside_their_domain_are_refused() {
    let res = coarse(8.0, 600.0);
    let kutta = ExitReceptivity::Kutta;
    let inviscid = coarse(8.0, f64::INFINITY);
    assert!(
        bickley_step_response_with_base(STATION, &inviscid, kutta, BaseSpreading::LaminarBickley)
            .is_err()
    );
    for slope in [0.0, -0.1, f64::NAN, f64::INFINITY] {
        let law = BaseSpreading::Linear { width_slope: slope };
        assert!(
            bickley_step_response_with_base(STATION, &res, kutta, law).is_err(),
            "{slope}"
        );
    }
    // Fletcher's slope at 20 b0 needs a band of 4 x 4.5 b0; 8 b0 is refused, 18 b0 admitted.
    let wide = BaseSpreading::Linear {
        width_slope: fletcher_slope(),
    };
    assert!(bickley_step_response_with_base(20.0, &res, kutta, wide).is_err());
    assert!(ForcedJet::with_base(coarse(18.0, 600.0), kutta, wide, 20.0).is_ok());
    for hold in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(
            ForcedJet::with_base(res, kutta, wide, hold).is_err(),
            "{hold}"
        );
    }
    // The parallel base ignores the hold, as before.
    assert!(ForcedJet::with_base(res, kutta, BaseSpreading::Parallel, f64::INFINITY).is_ok());
}
