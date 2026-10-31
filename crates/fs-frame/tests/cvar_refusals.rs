//! G1/G3: a frame design cannot hide a failed history or failed catalog recheck.

use fs_frame::cvar::{FrameCvarError, try_cvar_mass_min, try_ensemble_cvar, try_losses};
use fs_frame::{StoryFrame, StoryParams};
use fs_qty::{Dims, QtyAny};
use fs_scenario::ensemble::{SpectrumModel, StochasticEnsemble};

fn motions(duration_s: f64, omega_g: f64) -> StochasticEnsemble {
    StochasticEnsemble {
        name: "catalog-resonance-refusal".into(),
        seed: 29,
        members: 3,
        duration: QtyAny::new(duration_s, Dims([0, 0, 1, 0, 0, 0])),
        dt: QtyAny::new(0.02, Dims([0, 0, 1, 0, 0, 0])),
        model: SpectrumModel::KanaiTajimi {
            // Stay near the linear response while resolving section-dependent
            // resonance. These are numerical motion/model fixtures.
            s0: 1e-12,
            omega_g: QtyAny::new(omega_g, Dims([0, 0, -1, 0, 0, 0])),
            zeta_g: 0.02,
        },
    }
}

#[test]
fn invalid_or_overflowing_sections_cannot_become_zero_drift_designs() {
    let ensemble = motions(0.08, 12.5);
    let base = StoryParams::default();
    for scale in [f64::NAN, f64::INFINITY, 0.0, -1.0, f64::MAX] {
        assert!(matches!(
            try_losses(&ensemble, base, scale),
            Err(FrameCvarError::History { member: 0, .. })
        ));
    }
    assert!(matches!(
        try_losses(&ensemble, StoryParams { h: f64::NAN, ..base }, 1.0),
        Err(FrameCvarError::History { member: 0, .. })
    ));
}

#[test]
fn catalog_snap_into_resonance_must_pass_its_own_risk_limit() {
    let base = StoryParams::default();
    let omega = (StoryFrame::new(base).initial_stiffness() / base.mass).sqrt();
    let ensemble = motions(8.0, omega);
    let beta = 0.5;
    let low = try_ensemble_cvar(&ensemble, base, 0.25, beta).unwrap();
    let high = try_ensemble_cvar(&ensemble, base, 4.0, beta).unwrap();
    let resonant = try_ensemble_cvar(&ensemble, base, 1.0, beta).unwrap();
    let limit = 1.05 * low.max(high);
    assert!(
        resonant > limit && low > 0.0 && high > 0.0,
        "fixture must resolve resonance: low={low}, high={high}, snapped={resonant}"
    );
    // Both bisection endpoints pass. Snapping their accepted lower endpoint
    // to scale 1 moves the actual structure into the motion's resonance band.
    match try_cvar_mass_min(&ensemble, base, beta, limit, &[1.0]) {
        Err(FrameCvarError::CatalogRiskExceeded {
            scale,
            cvar_observed,
            limit: observed_limit,
        }) => {
            assert_eq!(scale, 1.0);
            assert_eq!(cvar_observed, resonant);
            assert_eq!(observed_limit, limit);
        }
        _ => panic!("an independently infeasible catalog design was published"),
    }
}
