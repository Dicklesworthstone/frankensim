//! The pair-mean frequency sampling of the multi-port radiation law, `modal_loading::PairMeanLaw`,
//! on a planar baffled Rayleigh bake of three modes: each entry is the reciprocal part of the row at
//! its pair's mean frequency, an owner's block is dropped, the band statement is the largest change
//! across each mode's half-power band, and the law refuses what it cannot sample.
//! Run with `cargo test -p fs-bem --test pair_mean_loading`.

// Exact float comparisons here assert entries carried unchanged from a bake row.
#![allow(clippy::float_cmp)]

use std::f64::consts::TAU;

use fs_bem::helmholtz::Medium;
use fs_bem::helmholtz::rayleigh_modal::AreaSample;
use fs_bem::modal_loading::{
    LoadingRoute, ModalLoadingBake, ModalLoadingError, PairMeanLaw, SurfaceModeTraces,
    modal_radiation_loading,
};
use fs_blake3::hash_domain;

/// A bar's first mode, its tube's first mode 1.9 Hz above it, and a mode two octaves up.
const HZ: [f64; 3] = [329.63, 331.53, 1_318.5];
/// Half-power half-widths [rad/s]: a slow bar mode, a fast tube mode, a slow upper mode.
const HALF_BANDS: [f64; 3] = [1.2, 25.0, 4.0];

fn omegas() -> [f64; 3] {
    HZ.map(|f| TAU * f)
}

/// Twelve cells of 2 cm^2 on a 4 x 3 grid: a piston, a rocking mode and a mode with a node line,
/// placed off the grid's centre lines so that no mutual entry vanishes by symmetry.
fn cells() -> Vec<AreaSample> {
    (0..12)
        .map(|i| {
            let (x, y) = (f64::from(i % 4), f64::from(i / 4));
            AreaSample {
                position_m: [0.02 * x, 0.025 * y],
                area_m2: 2.0e-4,
                mode_shape: vec![1.0, x - 1.2, 0.4 * (y - 0.8) * (x + 0.5) - 0.3],
            }
        })
        .collect()
}

#[expect(
    clippy::result_large_err,
    reason = "ModalLoadingError is the bake's refusal type; boxing it is outside this test"
)]
fn bake(samples: &[AreaSample], omegas: &[f64]) -> Result<ModalLoadingBake, ModalLoadingError> {
    modal_radiation_loading(
        SurfaceModeTraces::PlanarBaffled { samples },
        LoadingRoute::PlanarBaffledRayleigh,
        omegas,
        Medium::air(),
        hash_domain("analytic-fixture", "pair-mean loading".as_bytes()),
        1_000_000,
        || false,
    )
}

/// `(-Im / w, Re)` of `(Z_rs + Z_sr) / 2` at exactly `omega`, read from the bake by hand.
fn by_hand(bake: &ModalLoadingBake, omega: f64, r: usize, s: usize) -> (f64, f64) {
    let row = bake.at_omega(omega).expect("declared");
    assert_eq!(row.omega(), omega);
    let z = row.impedance();
    let (a, b) = (z[r * 3 + s], z[s * 3 + r]);
    (-(0.5 * (a.im + b.im)) / omega, 0.5 * (a.re + b.re))
}

#[test]
fn each_entry_is_the_reciprocal_part_at_its_pairs_mean_frequency() {
    let w = omegas();
    let law = PairMeanLaw::new(&w, &HALF_BANDS, &[None; 3]).expect("law");
    let frequencies = law.frequencies();
    // Six pair means and six band edges, none equal here, ascending.
    assert_eq!(frequencies.len(), 12);
    assert!(frequencies.windows(2).all(|pair| pair[0] < pair[1]));
    let bake = bake(&cells(), &frequencies).expect("bake");
    let loading = law.sample(&bake).expect("sample");
    for r in 0..3 {
        for s in 0..3 {
            let mean = (w[r] + w[s]) / 2.0;
            assert_eq!(law.pair_mean(r, s), mean);
            let (mass, real) = by_hand(&bake, mean, r, s);
            assert_eq!(loading.added_mass[r * 3 + s], mass, "M_a ({r}, {s})");
            assert_eq!(loading.resistance[r * 3 + s], real, "R ({r}, {s})");
            assert_eq!(
                loading.added_mass[r * 3 + s].to_bits(),
                loading.added_mass[s * 3 + r].to_bits()
            );
        }
    }
    // A baffled piston's own loading is a positive mass and a positive resistance (exp(-i w t):
    // Z = rho c S (R_1 - i X_1) with R_1 and X_1 positive), so a sign slip in the law shows.
    assert!(
        loading.added_mass[0] > 0.0 && loading.resistance[0] > 0.0,
        "piston: M_a {} R {}",
        loading.added_mass[0],
        loading.resistance[0]
    );
    // The loading is real and changes with frequency, so a wrong sample would show: the piston's
    // resistance at the upper mode's frequency is far above its value at its own.
    let (_, at_own) = by_hand(&bake, w[0], 0, 0);
    let (_, at_upper) = by_hand(&bake, w[2], 0, 0);
    assert!(
        at_own > 0.0 && at_upper > 4.0 * at_own,
        "{at_own} -> {at_upper}"
    );
}

#[test]
fn the_modes_of_one_owner_lose_their_shared_entries() {
    let w = omegas();
    let free = PairMeanLaw::new(&w, &HALF_BANDS, &[None; 3]).expect("law");
    let owned = PairMeanLaw::new(&w, &HALF_BANDS, &[None, Some(7), Some(7)]).expect("law");
    assert!(owned.shared(1, 2) && owned.shared(2, 2) && !owned.shared(0, 1) && !owned.shared(0, 0));
    let bake = bake(&cells(), &free.frequencies()).expect("bake");
    let (all, kept) = (
        free.sample(&bake).expect("free"),
        owned.sample(&bake).expect("owned"),
    );
    for r in 0..3 {
        for s in 0..3 {
            let k = r * 3 + s;
            if r > 0 && s > 0 {
                assert_eq!(
                    (kept.added_mass[k], kept.resistance[k]),
                    (0.0, 0.0),
                    "({r}, {s})"
                );
            } else {
                assert_eq!(
                    (kept.added_mass[k], kept.resistance[k]),
                    (all.added_mass[k], all.resistance[k])
                );
                assert!(all.resistance[k] != 0.0, "({r}, {s}) carries loading");
            }
        }
    }
}

#[test]
fn the_band_statement_is_the_largest_change_across_each_band() {
    let w = omegas();
    let owners = [None, Some(1), Some(1)];
    let law = PairMeanLaw::new(&w, &HALF_BANDS, &owners).expect("law");
    let bake = bake(&cells(), &law.frequencies()).expect("bake");
    let loading = law.sample(&bake).expect("sample");
    for r in 0..3 {
        let (mut mass, mut real, mut scale_mass, mut scale_real) =
            (0.0_f64, 0.0_f64, 0.0_f64, 0.0_f64);
        for edge in [w[r] - HALF_BANDS[r], w[r] + HALF_BANDS[r]] {
            for s in 0..3 {
                if r > 0 && s > 0 {
                    continue;
                }
                let (at_edge_mass, at_edge_real) = by_hand(&bake, edge, r, s);
                let (kept_mass, kept_real) = by_hand(&bake, (w[r] + w[s]) / 2.0, r, s);
                mass = mass.max((at_edge_mass - kept_mass).abs());
                real = real.max((at_edge_real - kept_real).abs());
                scale_mass = scale_mass.max(kept_mass.abs());
                scale_real = scale_real.max(kept_real.abs());
            }
        }
        let band = loading.band[r];
        assert_eq!(band.mode, r);
        assert_eq!((band.added_mass, band.resistance), (mass, real), "mode {r}");
        assert_eq!(
            (band.largest_added_mass, band.largest_resistance),
            (scale_mass, scale_real)
        );
        assert!(
            band.resistance > 0.0,
            "mode {r}: the loading moves across the band"
        );
    }
}

#[test]
fn equal_pair_means_are_declared_once() {
    let w = TAU * 440.0;
    let law = PairMeanLaw::new(&[w, w], &[0.0, 3.0], &[None, None]).expect("law");
    // A zero half-width puts both edges of mode 0 on its pair mean.
    assert_eq!(law.frequencies(), vec![w - 3.0, w, w + 3.0]);
}

#[test]
fn the_law_refuses_what_it_cannot_sample() {
    let w = omegas();
    let invalid = |result: Result<PairMeanLaw, ModalLoadingError>| {
        assert!(
            matches!(result, Err(ModalLoadingError::Invalid(_))),
            "{result:?}"
        );
    };
    invalid(PairMeanLaw::new(&[], &[], &[]));
    invalid(PairMeanLaw::new(&w, &HALF_BANDS[..2], &[None; 3]));
    invalid(PairMeanLaw::new(&w, &HALF_BANDS, &[None; 2]));
    for bad in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        invalid(PairMeanLaw::new(&[bad], &[0.0], &[None]));
    }
    for bad in [-0.5, w[0], 2.0 * w[0], f64::NAN, f64::INFINITY] {
        invalid(PairMeanLaw::new(&[w[0]], &[bad], &[None]));
    }
    let law = PairMeanLaw::new(&w, &HALF_BANDS, &[None; 3]).expect("law");
    // A bake of the pair means alone lacks the band edges.
    let mut means: Vec<f64> = (0..3)
        .flat_map(|r| (r..3).map(move |s| (r, s)))
        .map(|(r, s)| law.pair_mean(r, s))
        .collect();
    means.sort_by(f64::total_cmp);
    let short = bake(&cells(), &means).expect("bake");
    assert!(matches!(
        law.sample(&short),
        Err(ModalLoadingError::Invalid(_))
    ));
    // A bake of two modes does not hold the law's three.
    let two: Vec<AreaSample> = cells()
        .into_iter()
        .map(|mut cell| {
            cell.mode_shape.truncate(2);
            cell
        })
        .collect();
    let other = bake(&two, &law.frequencies()).expect("bake");
    assert!(matches!(
        law.sample(&other),
        Err(ModalLoadingError::Invalid(_))
    ));
}
