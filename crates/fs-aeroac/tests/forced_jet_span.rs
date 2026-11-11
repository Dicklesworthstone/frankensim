//! The forced jet's station span (fs_aeroac::forced_jet::bickley_step_span, MB 23643): ONE
//! integration records every grid node of `[near, far]`, and `response_at` rebuilds the step
//! response at any station inside it by `displacement_at`'s own arithmetic. A span of one station is
//! the single-station response bit for bit under the parallel and the laminar Bickley bases and
//! both exits; a parallel base never reads
//! the station, so every station of a parallel span is its own single-station integration bit for
//! bit; a spread base held past `far` gives the single-station response at `far` bit for bit and
//! moves the near station only slightly; stations outside their domain are refused.

use fs_aeroac::forced_jet::{
    BaseSpreading, ExitReceptivity, ForcedJetResolution, JetStepResponse,
    bickley_step_response_with_base, bickley_step_span,
};

const WORKERS: usize = 4;
/// Numerical choices inside the coarse box: a span of one exit half-width around 6.
const NEAR: f64 = 5.5;
const FAR: f64 = 6.5;

/// Numerics, not physics: the coarse row of `tests/forced_jet_spreading.rs` (`spacing 0.25`,
/// `window 40`, 256 cells, Courant 0.8, half-height 8, Re 600).
fn coarse() -> ForcedJetResolution {
    ForcedJetResolution {
        spacing: 0.25,
        streamwise_cells: 256,
        half_height: 8.0,
        window: 40.0,
        courant: 0.8,
        reynolds: 600.0,
        workers: WORKERS,
    }
}

fn exits() -> [ExitReceptivity; 2] {
    [
        ExitReceptivity::Kutta,
        ExitReceptivity::Measured {
            displacement_per_velocity: 4.0 / 3.0,
        },
    ]
}

fn bits(step: &JetStepResponse) -> Vec<u64> {
    step.values().iter().map(|v| v.to_bits()).collect()
}

/// The values, the integrals the hat weights read (at the records and between them) and the
/// metadata are the single-station response's.
fn assert_same(span: &JetStepResponse, single: &JetStepResponse, what: &str) {
    assert_eq!(bits(span), bits(single), "{what}: values");
    assert_eq!(span.station().to_bits(), single.station().to_bits(), "{what}: station");
    assert_eq!(span.time_step().to_bits(), single.time_step().to_bits(), "{what}: time step");
    assert_eq!(span.spreading(), single.spreading(), "{what}: spreading");
    assert_eq!(span.exit(), single.exit(), "{what}: exit");
    for k in 0..=2 * span.values().len() {
        let tau = 0.5 * k as f64 * span.time_step();
        assert_eq!(span.integral_at(tau).to_bits(), single.integral_at(tau).to_bits(), "{what}: integral at {tau}");
    }
}

#[test]
fn a_span_of_one_station_is_the_single_station_response_bit_for_bit() {
    let res = coarse();
    for spreading in [BaseSpreading::Parallel, BaseSpreading::LaminarBickley] {
        for exit in exits() {
            for station in [NEAR, 6.0, 6.37] {
                let span = bickley_step_span(station, station, &res, exit, spreading).unwrap();
                let single = bickley_step_response_with_base(station, &res, exit, spreading).unwrap();
                assert_same(&span.response_at(station).unwrap(), &single, &format!("{spreading:?} {exit:?} {station}"));
            }
        }
    }
}

#[test]
fn every_station_of_a_parallel_span_is_its_own_integration_bit_for_bit() {
    let res = coarse();
    for exit in exits() {
        let span = bickley_step_span(NEAR, FAR, &res, exit, BaseSpreading::Parallel).unwrap();
        assert_eq!(span.stations(), (NEAR, FAR));
        for station in [NEAR, 5.62, 6.0, 6.37, FAR] {
            let single = bickley_step_response_with_base(station, &res, exit, BaseSpreading::Parallel).unwrap();
            assert_same(&span.response_at(station).unwrap(), &single, &format!("{exit:?} {station}"));
        }
    }
}

/// The base is held past `far`, so `far` reads the single-station base exactly. A nearer station
/// reads a base that keeps spreading out to `far`: the difference is real and small.
#[test]
fn a_spread_span_is_exact_at_its_far_end_and_close_inside() {
    let res = coarse();
    for exit in exits() {
        let span = bickley_step_span(NEAR, FAR, &res, exit, BaseSpreading::LaminarBickley).unwrap();
        let far = bickley_step_response_with_base(FAR, &res, exit, BaseSpreading::LaminarBickley).unwrap();
        assert_same(&span.response_at(FAR).unwrap(), &far, &format!("{exit:?} far"));
        let near = bickley_step_response_with_base(NEAR, &res, exit, BaseSpreading::LaminarBickley).unwrap();
        let read = span.response_at(NEAR).unwrap();
        let rel = (read.held_gain() - near.held_gain()).abs() / near.held_gain().abs();
        let g_read = read.transfer(0.3).unwrap();
        let g_near = near.transfer(0.3).unwrap();
        let rel_g = (g_read - g_near).abs() / g_near.abs();
        println!("{exit:?}: near station held gain {:.9e} (span) vs {:.9e} (single), relative {rel:.3e}; |dG(0.3)|/|G| {rel_g:.3e}",
            read.held_gain(), near.held_gain());
        assert_ne!(bits(&read), bits(&near), "{exit:?}: the spread base differs past the near station");
        assert!(rel < 1e-2 && rel_g < 1e-2, "{exit:?}: {rel:e} {rel_g:e}");
    }
}

#[test]
fn stations_outside_their_domain_are_refused() {
    let res = coarse();
    let (exit, base) = (ExitReceptivity::Kutta, BaseSpreading::Parallel);
    for (near, far) in [(FAR, NEAR), (0.0, FAR), (-1.0, FAR), (f64::NAN, FAR), (NEAR, f64::INFINITY), (0.1, FAR), (NEAR, 1e6)] {
        assert!(bickley_step_span(near, far, &res, exit, base).is_err(), "({near}, {far})");
    }
    let span = bickley_step_span(NEAR, FAR, &res, exit, base).unwrap();
    for station in [NEAR - 1e-9, FAR + 1e-9, f64::NAN] {
        assert!(span.response_at(station).is_err(), "{station}");
    }
}
