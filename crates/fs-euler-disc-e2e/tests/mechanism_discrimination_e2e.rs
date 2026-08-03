//! G0/G2 checks of declared reduced-decay equations through the public API.
//!
//! These tests check numerical implementation of caller-declared dry-contour
//! and Bildsten energy-only equations. They do not establish an emergent
//! physical mechanism, experimental agreement, target ranking, or any Mould
//! outcome.

use fs_euler_disc_e2e::{
    BILDSTEN_PUBLISHED_POWER_COEFFICIENT, BildstenBoundaryLayerChannel, ChannelCrossoverDiagnostic,
    DryContourChannel, REDUCED_DECAY_MODEL_ID, ReducedDecayInput, ReducedDecayTerminal,
    STANDARD_GRAVITY_M_PER_S2, channel_crossover_diagnostic, run_reduced_decay,
};
use fs_tribo::{InputAuthority, InterfaceMedium, InterfaceSystemRef};

const MASS_KG: f64 = 0.12;
const RADIUS_M: f64 = 0.038;
const INITIAL_THETA_RAD: f64 = 0.08;
const CUTOFF_THETA_RAD: f64 = 0.001;
const TIMESTEP_S: f64 = 1.0e-5;

fn dry_channel(force_n: f64) -> DryContourChannel {
    DryContourChannel {
        interface: InterfaceSystemRef::new(
            "mechanism-discrimination/disc->support",
            "mechanism-discrimination/synthetic-history",
            "synthetic/mechanism-discrimination-v1",
            InputAuthority::SyntheticFixture,
            InterfaceMedium::Dry,
        )
        .expect("declared synthetic dry interface is valid"),
        normal_force_n: MASS_KG * STANDARD_GRAVITY_M_PER_S2,
        contour_force_n: force_n,
    }
}

fn bildsten_channel() -> BildstenBoundaryLayerChannel {
    BildstenBoundaryLayerChannel {
        source_id: "synthetic/bildsten-energy-only-v1".to_owned(),
        density_kg_per_m3: 1.2,
        dynamic_viscosity_pa_s: 1.8e-5,
        dimensionless_prefactor: 1.0,
    }
}

fn input(dry_force_n: Option<f64>, with_bildsten: bool) -> ReducedDecayInput {
    ReducedDecayInput {
        mass_kg: MASS_KG,
        radius_m: RADIUS_M,
        gravity_m_per_s2: STANDARD_GRAVITY_M_PER_S2,
        initial_theta_rad: INITIAL_THETA_RAD,
        validity_cutoff_theta_rad: CUTOFF_THETA_RAD,
        timestep_s: TIMESTEP_S,
        maximum_steps: 120_000,
        small_angle_oracle_source_id: "synthetic/small-angle-oracle-v1".to_owned(),
        dry_contour: dry_force_n.map(dry_channel),
        bildsten_boundary_layer: with_bildsten.then(bildsten_channel),
    }
}

fn assert_monotone_accounting(label: &str, input: &ReducedDecayInput) {
    let run = run_reduced_decay(input).expect("declared reduced equation run");
    assert_eq!(
        run.terminal,
        ReducedDecayTerminal::ValidityCutoff,
        "{label}: expected positive-cutoff termination, got {:?}",
        run.terminal
    );
    assert!(
        run.energy_closure_residual_j.abs() < 1.0e-11,
        "{label}: encoded P*dt closure residual={} J",
        run.energy_closure_residual_j
    );
    for pair in run.samples.windows(2) {
        let previous = pair[0];
        let next = pair[1];
        assert!(
            next.energy_j <= previous.energy_j,
            "{label}: energy rose from {} to {} J at t={} s",
            previous.energy_j,
            next.energy_j,
            next.time_s
        );
        assert!(
            next.work.total_j() >= previous.work.total_j(),
            "{label}: total work fell from {} to {} J at t={} s",
            previous.work.total_j(),
            next.work.total_j(),
            next.time_s
        );
        assert!(
            next.work.dry_contour_j >= previous.work.dry_contour_j,
            "{label}: dry work fell from {} to {} J",
            previous.work.dry_contour_j,
            next.work.dry_contour_j
        );
        assert!(
            next.work.bildsten_boundary_layer_j >= previous.work.bildsten_boundary_layer_j,
            "{label}: Bildsten work fell from {} to {} J",
            previous.work.bildsten_boundary_layer_j,
            next.work.bildsten_boundary_layer_j
        );
    }
    let final_sample = run
        .final_sample()
        .expect("successful run retains final sample");
    if input.dry_contour.is_some() {
        assert!(
            final_sample.work.dry_contour_j > 0.0,
            "{label}: declared dry channel accumulated {} J",
            final_sample.work.dry_contour_j
        );
    } else {
        assert_eq!(
            final_sample.work.dry_contour_j, 0.0,
            "{label}: dry ablation leaked work"
        );
    }
    if input.bildsten_boundary_layer.is_some() {
        assert!(
            final_sample.work.bildsten_boundary_layer_j > 0.0,
            "{label}: declared Bildsten channel accumulated {} J",
            final_sample.work.bildsten_boundary_layer_j
        );
    } else {
        assert_eq!(
            final_sample.work.bildsten_boundary_layer_j, 0.0,
            "{label}: Bildsten ablation leaked work"
        );
    }
}

#[test]
fn g0_declared_channels_are_separated_and_accounted_monotonically() {
    assert_monotone_accounting("dry-only", &input(Some(0.001), false));
    assert_monotone_accounting("bildsten-air-only", &input(None, true));
    assert_monotone_accounting("combined", &input(Some(0.001), true));
}

#[test]
fn g0_reduced_runs_are_deterministic_and_retain_no_claim_provenance() {
    let scenario = input(Some(0.001), true);
    let first = run_reduced_decay(&scenario).expect("first deterministic run");
    let second = run_reduced_decay(&scenario).expect("second deterministic run");
    assert_eq!(
        first, second,
        "same declared input must retain identical samples and work"
    );

    let provenance = &first.provenance;
    assert_eq!(provenance.model_id, REDUCED_DECAY_MODEL_ID);
    assert_eq!(
        provenance.small_angle_oracle_source_id,
        "synthetic/small-angle-oracle-v1"
    );
    assert_eq!(provenance.model_authority, "numerical-reference-only");
    assert_eq!(provenance.physical_validation, "not-claimed");
    assert_eq!(provenance.cancellation_capability, "not-implemented");
    assert_eq!(
        provenance.dry_source_id.as_deref(),
        Some("synthetic/mechanism-discrimination-v1")
    );
    assert_eq!(
        provenance.bildsten_source_id.as_deref(),
        Some("synthetic/bildsten-energy-only-v1")
    );
    assert_eq!(provenance.bildsten_multiplier_authority, "caller-declared");
    assert_eq!(
        provenance.dry_authority,
        Some(InputAuthority::SyntheticFixture),
        "synthetic fixture provenance must remain synthetic rather than experimental authority"
    );
}

#[test]
fn g2_crossover_diagnostic_closes_the_declared_power_equation() {
    let scenario = input(Some(0.001), true);
    let diagnostic = channel_crossover_diagnostic(&scenario).expect("crossover diagnostic");
    let theta_rad = match diagnostic {
        ChannelCrossoverDiagnostic::AtInclination { theta_rad } => theta_rad,
        other => panic!("expected declared-model crossover, got {other:?}"),
    };
    let dry_power_w = scenario
        .dry_contour
        .as_ref()
        .expect("dry channel")
        .contour_force_n
        * scenario.radius_m
        * theta_rad.cos()
        * scenario
            .omega_rad_s(theta_rad)
            .expect("declared precession equation");
    let air = scenario
        .bildsten_boundary_layer
        .as_ref()
        .expect("Bildsten channel");
    let bildsten_power_w = air.dimensionless_prefactor
        * BILDSTEN_PUBLISHED_POWER_COEFFICIENT
        * (air.dynamic_viscosity_pa_s * air.density_kg_per_m3).sqrt()
        * scenario.gravity_m_per_s2.powf(1.25)
        * scenario.radius_m.powf(2.75)
        * theta_rad.powf(-1.25);
    let residual_w = dry_power_w - bildsten_power_w;
    let scale_w = dry_power_w.abs().max(bildsten_power_w.abs());
    assert!(
        residual_w.abs() <= 1.0e-10 * scale_w,
        "declared crossover residual={} W (dry={} W, Bildsten={} W, theta={} rad)",
        residual_w,
        dry_power_w,
        bildsten_power_w,
        theta_rad
    );
    assert!(
        (CUTOFF_THETA_RAD..=INITIAL_THETA_RAD).contains(&theta_rad),
        "crossover theta={} rad lies outside declared interval [{}, {}] rad",
        theta_rad,
        CUTOFF_THETA_RAD,
        INITIAL_THETA_RAD
    );
}
