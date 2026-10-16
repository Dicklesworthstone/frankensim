//! Numerical self-consistency checks for the public reduced-decay API.
//!
//! Exponent recovery and `P * dt` energy closure are internal consistency
//! checks of the stated reduced equations. They are not independent emergence,
//! experimental agreement, or physical validation.

use fs_euler_disc_e2e::reduced_decay::{
    THORNE_2026_AMBIENT_AIR_DENSITY_KG_PER_M3, THORNE_2026_DECLARED_AIR_VISCOSITY_PA_S,
    THORNE_2026_FITTED_ROLLING_COEFFICIENT, THORNE_2026_SOURCE_ID,
    THORNE_2026_STEEL_DISC_DIAMETER_M, THORNE_2026_STEEL_DISC_FILLET_RADIUS_M,
    THORNE_2026_STEEL_DISC_MASS_KG, THORNE_2026_STEEL_DISC_THICKNESS_M,
    THORNE_2026_VACUUM_AIR_DENSITY_KG_PER_M3, Thorne2026SteelGlassBenchmark,
    run_thorne_2026_steel_glass_benchmark, thorne_2026_channel_crossover_diagnostic,
    thorne_2026_refinement_evidence,
};
use fs_euler_disc_e2e::specimen::DiscProfileSpec;
use fs_euler_disc_e2e::{
    BildstenBoundaryLayerChannel, ChannelCrossoverDiagnostic, ChannelCrossoverNotComparable,
    DryContourChannel, ReducedDecayError, ReducedDecayInput, ReducedDecayTerminal,
    STANDARD_GRAVITY_M_PER_S2, channel_crossover_diagnostic, refinement_evidence,
    run_reduced_decay, structured_runner_output,
};
use fs_rep_frep::SquatDiscEdgeTreatment;
use fs_tribo::{InputAuthority, InterfaceMedium, InterfaceSystemRef};

fn dry_channel(force_n: f64) -> DryContourChannel {
    DryContourChannel {
        interface: InterfaceSystemRef::new(
            "reduced-decay/test-disc->test-support",
            "reduced-decay/test-history",
            "synthetic/numerical-reference",
            InputAuthority::SyntheticFixture,
            InterfaceMedium::Dry,
        )
        .expect("fixed test interface is valid"),
        normal_force_n: 1.0,
        contour_force_n: force_n,
    }
}

fn air_channel(density_kg_per_m3: f64) -> BildstenBoundaryLayerChannel {
    BildstenBoundaryLayerChannel {
        source_id: "synthetic/bildsten-energy-only".to_owned(),
        density_kg_per_m3,
        dynamic_viscosity_pa_s: 1.8e-5,
        dimensionless_prefactor: 1.0,
    }
}

fn input(dry: Option<f64>, air: Option<f64>) -> ReducedDecayInput {
    ReducedDecayInput {
        mass_kg: 0.12,
        radius_m: 0.038,
        gravity_m_per_s2: STANDARD_GRAVITY_M_PER_S2,
        initial_theta_rad: 0.05,
        validity_cutoff_theta_rad: 0.001,
        timestep_s: 1.0e-5,
        maximum_steps: 100_000,
        small_angle_oracle_source_id: "synthetic/small-angle-oracle-v1".to_owned(),
        dry_contour: dry.map(dry_channel),
        bildsten_boundary_layer: air.map(air_channel),
    }
}

fn exponent_from_terminal_times(low: f64, high: f64, low_time_s: f64, high_time_s: f64) -> f64 {
    (high / low).ln() / (high_time_s / low_time_s).ln()
}

#[test]
fn e2e_reduced_decay_contact_only_recovers_two_thirds_reference_exponent() {
    let mut low = input(Some(0.002), None);
    low.initial_theta_rad = 0.03;
    low.validity_cutoff_theta_rad = 1.0e-5;
    let mut high = low.clone();
    high.initial_theta_rad = 0.06;
    let low_run = run_reduced_decay(&low).expect("contact-only low run");
    let high_run = run_reduced_decay(&high).expect("contact-only high run");
    let exponent = exponent_from_terminal_times(
        low.initial_theta_rad,
        high.initial_theta_rad,
        low_run.final_sample().expect("low final sample").time_s,
        high_run.final_sample().expect("high final sample").time_s,
    );
    assert!((exponent - 2.0 / 3.0).abs() < 0.01, "got {exponent}");
    assert_eq!(low_run.terminal, ReducedDecayTerminal::ValidityCutoff);
    assert_eq!(
        low_run
            .final_sample()
            .expect("low final sample")
            .powers
            .bildsten_boundary_layer_w,
        0.0
    );
}

#[test]
fn e2e_reduced_decay_air_only_recovers_four_ninths_reference_exponent() {
    let mut low = input(None, Some(1.2));
    low.initial_theta_rad = 0.03;
    low.validity_cutoff_theta_rad = 1.0e-5;
    low.maximum_steps = 200_000;
    let mut high = low.clone();
    high.initial_theta_rad = 0.06;
    let low_run = run_reduced_decay(&low).expect("air-only low run");
    let high_run = run_reduced_decay(&high).expect("air-only high run");
    let exponent = exponent_from_terminal_times(
        low.initial_theta_rad,
        high.initial_theta_rad,
        low_run.final_sample().expect("low final sample").time_s,
        high_run.final_sample().expect("high final sample").time_s,
    );
    assert!((exponent - 4.0 / 9.0).abs() < 0.01, "got {exponent}");
    assert_eq!(low_run.terminal, ReducedDecayTerminal::ValidityCutoff);
    assert_eq!(high_run.terminal, ReducedDecayTerminal::ValidityCutoff);
    assert_eq!(
        low_run
            .final_sample()
            .expect("low final sample")
            .powers
            .dry_contour_w,
        0.0
    );
}

#[test]
fn e2e_reduced_decay_combined_channels_cross_over_without_blending_identity() {
    let mut combined = input(Some(0.002), Some(1.2));
    combined.initial_theta_rad = 0.08;
    let run = run_reduced_decay(&combined).expect("combined run");
    let start = run.samples.first().expect("initial sample");
    let end = run
        .samples
        .get(run.samples.len() - 2)
        .expect("pre-cutoff sample");
    assert!(start.powers.dry_contour_w > start.powers.bildsten_boundary_layer_w);
    assert!(end.powers.bildsten_boundary_layer_w > end.powers.dry_contour_w);
    assert!(run.final_sample().expect("final sample").work.dry_contour_j > 0.0);
    assert!(
        run.final_sample()
            .expect("final sample")
            .work
            .bildsten_boundary_layer_j
            > 0.0
    );
    assert!(run.energy_closure_residual_j.abs() < 1.0e-12);
    assert!(matches!(
        channel_crossover_diagnostic(&combined),
        Ok(ChannelCrossoverDiagnostic::AtInclination { theta_rad })
            if theta_rad > combined.validity_cutoff_theta_rad
                && theta_rad < combined.initial_theta_rad
    ));
}

#[test]
fn e2e_reduced_decay_crossover_diagnostic_handles_ablations_and_no_crossing() {
    assert_eq!(
        channel_crossover_diagnostic(&input(None, Some(1.2))).expect("air-only diagnostic"),
        ChannelCrossoverDiagnostic::NotComparable {
            reason: ChannelCrossoverNotComparable::MissingDryContour,
        }
    );
    assert_eq!(
        channel_crossover_diagnostic(&input(Some(0.002), None)).expect("dry-only diagnostic"),
        ChannelCrossoverDiagnostic::NotComparable {
            reason: ChannelCrossoverNotComparable::MissingBildstenBoundaryLayer,
        }
    );
    let mut dry_dominant = input(Some(1.0), Some(1.2));
    dry_dominant
        .bildsten_boundary_layer
        .as_mut()
        .expect("Bildsten channel")
        .dimensionless_prefactor = 1.0e-8;
    assert_eq!(
        channel_crossover_diagnostic(&dry_dominant).expect("no-crossing diagnostic"),
        ChannelCrossoverDiagnostic::NoneWithinInterval,
    );
}

#[test]
fn e2e_reduced_decay_ablation_isolates_real_channels() {
    let contact = run_reduced_decay(&input(Some(0.002), None)).expect("contact run");
    let air = run_reduced_decay(&input(None, Some(1.2))).expect("air run");
    assert_eq!(
        contact
            .final_sample()
            .expect("contact final sample")
            .work
            .bildsten_boundary_layer_j,
        0.0
    );
    assert_eq!(
        air.final_sample()
            .expect("air final sample")
            .work
            .dry_contour_j,
        0.0
    );
    assert!(
        contact
            .final_sample()
            .expect("contact final sample")
            .work
            .dry_contour_j
            > 0.0
    );
    assert!(
        air.final_sample()
            .expect("air final sample")
            .work
            .bildsten_boundary_layer_j
            > 0.0
    );
}

#[test]
fn e2e_reduced_decay_mass_radius_and_gas_density_scaling_are_explicit() {
    let base = input(None, Some(1.0));
    let base_sample = run_reduced_decay(&base).expect("base run").samples[0];
    let mut mass_double = base.clone();
    mass_double.mass_kg *= 2.0;
    let mass_sample = run_reduced_decay(&mass_double).expect("mass run").samples[0];
    assert!((mass_sample.energy_j / base_sample.energy_j - 2.0).abs() < 1.0e-12);
    assert!((mass_sample.powers.total_w() / base_sample.powers.total_w() - 1.0).abs() < 1.0e-12);

    let mut radius_double = base.clone();
    radius_double.radius_m *= 2.0;
    let radius_sample = run_reduced_decay(&radius_double)
        .expect("radius run")
        .samples[0];
    assert!(
        (radius_sample.omega_rad_s / base_sample.omega_rad_s - 1.0 / 2.0_f64.sqrt()).abs()
            < 1.0e-12
    );
    assert!(
        (radius_sample.powers.bildsten_boundary_layer_w
            / base_sample.powers.bildsten_boundary_layer_w
            - 2.0_f64.powf(11.0 / 4.0))
        .abs()
            < 1.0e-10
    );

    let mut density_quadruple = base.clone();
    density_quadruple
        .bildsten_boundary_layer
        .as_mut()
        .expect("air present")
        .density_kg_per_m3 *= 4.0;
    let density_sample = run_reduced_decay(&density_quadruple)
        .expect("density run")
        .samples[0];
    assert!(
        (density_sample.powers.bildsten_boundary_layer_w
            / base_sample.powers.bildsten_boundary_layer_w
            - 2.0)
            .abs()
            < 1.0e-12
    );
}

#[test]
fn e2e_reduced_decay_multiplier_one_is_the_published_bildsten_law() {
    let scenario = input(None, Some(1.2));
    let sample = run_reduced_decay(&scenario)
        .expect("published-law run")
        .samples[0];
    let channel = scenario
        .bildsten_boundary_layer
        .as_ref()
        .expect("Bildsten channel");
    let expected_w = 4.0
        * (channel.dynamic_viscosity_pa_s * channel.density_kg_per_m3).sqrt()
        * scenario.gravity_m_per_s2.powf(1.25)
        * scenario.radius_m.powf(2.75)
        * scenario.initial_theta_rad.powf(-1.25);
    assert!((sample.powers.bildsten_boundary_layer_w / expected_w - 1.0).abs() < 1.0e-12);
}

#[test]
fn e2e_reduced_decay_timestep_refinement_is_bounded_and_closes_energy() {
    let evidence = refinement_evidence(&input(Some(0.002), Some(1.2))).expect("refinement");
    assert_eq!(
        evidence.coarse.terminal,
        ReducedDecayTerminal::ValidityCutoff
    );
    assert_eq!(evidence.fine.terminal, ReducedDecayTerminal::ValidityCutoff);
    assert!(evidence.terminal_time_difference_s < 5.0e-5);
    assert!(evidence.total_work_difference_j < 1.0e-9);
    assert!(evidence.coarse.energy_closure_residual_j.abs() < 1.0e-12);
    assert!(evidence.fine.energy_closure_residual_j.abs() < 1.0e-12);
}

#[test]
fn e2e_reduced_decay_refinement_refuses_non_cutoff_terminal_evidence() {
    let mut step_limited = input(Some(0.002), Some(1.2));
    step_limited.maximum_steps = 1;
    assert!(matches!(
        refinement_evidence(&step_limited),
        Err(ReducedDecayError::RefinementIncompleteTerminal {
            terminal: ReducedDecayTerminal::StepBudgetExhausted
        })
    ));
}

#[test]
fn e2e_reduced_decay_invalid_inputs_and_cutoff_are_structured() {
    let mut invalid = input(Some(0.002), None);
    invalid.mass_kg = 0.0;
    assert!(matches!(
        run_reduced_decay(&invalid),
        Err(ReducedDecayError::InvalidInput { field: "mass_kg" })
    ));
    let mut below_cutoff = input(Some(0.002), None);
    below_cutoff.initial_theta_rad = below_cutoff.validity_cutoff_theta_rad;
    assert!(matches!(
        run_reduced_decay(&below_cutoff),
        Err(ReducedDecayError::InitialStateOutsideValidity { .. })
    ));
    let mut outside_small_angle = input(Some(0.002), None);
    outside_small_angle.initial_theta_rad = 0.21;
    assert!(matches!(
        run_reduced_decay(&outside_small_angle),
        Err(ReducedDecayError::InvalidInput {
            field: "initial_theta_rad_small_angle"
        })
    ));
    let no_channel = input(None, None);
    assert!(matches!(
        run_reduced_decay(&no_channel),
        Err(ReducedDecayError::NoActiveChannel)
    ));
    let mut zero_multiplier = input(None, Some(1.2));
    zero_multiplier
        .bildsten_boundary_layer
        .as_mut()
        .expect("Bildsten channel")
        .dimensionless_prefactor = 0.0;
    assert!(matches!(
        run_reduced_decay(&zero_multiplier),
        Err(ReducedDecayError::InvalidInput {
            field: "bildsten.dimensionless_prefactor"
        })
    ));
    let mut zero_contour_force = input(Some(0.002), None);
    zero_contour_force
        .dry_contour
        .as_mut()
        .expect("dry channel")
        .contour_force_n = 0.0;
    assert!(matches!(
        run_reduced_decay(&zero_contour_force),
        Err(ReducedDecayError::InvalidInput {
            field: "dry_contour.contour_force_n"
        })
    ));
    let mut unsupported_unloaded_contour = input(Some(0.002), None);
    unsupported_unloaded_contour
        .dry_contour
        .as_mut()
        .expect("dry channel")
        .normal_force_n = 0.0;
    assert!(matches!(
        run_reduced_decay(&unsupported_unloaded_contour),
        Err(ReducedDecayError::InvalidInput {
            field: "dry_contour.normal_force_n"
        })
    ));
}

#[test]
fn e2e_reduced_decay_runner_output_is_structured_and_deterministic() {
    let scenario = input(Some(0.002), Some(1.2));
    let run = run_reduced_decay(&scenario).expect("runner run");
    let refinement = refinement_evidence(&scenario).expect("runner refinement");
    let first = structured_runner_output(&run, &refinement).expect("structured output");
    let second = structured_runner_output(&run, &refinement).expect("structured output");
    assert_eq!(first, second);
    assert!(
        first.starts_with("schema=reduced-decay-v1 model_id=euler-disc-small-angle-late-stage-v2 ")
    );
    assert!(first.contains("model_authority=numerical-reference-only"));
    assert!(first.contains("physical_validation=not-claimed"));
    assert!(first.contains("small_angle_oracle_source_id=synthetic/small-angle-oracle-v1"));
    assert!(first.contains("dry_authority=Some(SyntheticFixture)"));
    assert!(first.contains("bildsten_source_id=synthetic/bildsten-energy-only"));
    assert!(first.contains("terminal=ValidityCutoff"));
    assert!(first.contains("dry_work_j="));
    assert!(first.contains("bildsten_work_j="));
    assert!(first.contains("closure_residual_j="));
    assert!(first.contains("\nrefinement_terminal_time_difference_s="));
    assert!(first.ends_with("evidence_scope=numerical-self-consistency-only"));
}

#[test]
fn e2e_reduced_decay_distinguishes_oracle_and_bildsten_sources() {
    let nominal = ReducedDecayInput::nominal_reference().expect("nominal reference");
    assert_eq!(
        nominal.small_angle_oracle_source_id,
        "analytic/euler-disc-small-angle-oracle-v1"
    );
    assert_eq!(
        nominal
            .bildsten_boundary_layer
            .as_ref()
            .expect("Bildsten channel")
            .source_id,
        "doi:10.1103/PhysRevE.66.056309"
    );
    let mut blank_oracle = input(Some(0.002), None);
    blank_oracle.small_angle_oracle_source_id = " \t".to_owned();
    assert!(matches!(
        run_reduced_decay(&blank_oracle),
        Err(ReducedDecayError::MissingIdentity {
            field: "small_angle_oracle_source_id"
        })
    ));
}

#[test]
fn e2e_thorne_2026_benchmark_binds_reported_specimen_and_direct_power_laws() {
    let benchmark = Thorne2026SteelGlassBenchmark::ambient().expect("ambient benchmark");
    assert_eq!(benchmark.specimen.source_id, THORNE_2026_SOURCE_ID);
    assert_eq!(
        benchmark.specimen.diameter_m.to_bits(),
        THORNE_2026_STEEL_DISC_DIAMETER_M.to_bits()
    );
    assert_eq!(
        benchmark.specimen.thickness_m.to_bits(),
        THORNE_2026_STEEL_DISC_THICKNESS_M.to_bits()
    );
    assert_eq!(
        benchmark.specimen.mass_kg.to_bits(),
        THORNE_2026_STEEL_DISC_MASS_KG.to_bits()
    );
    assert_eq!(
        benchmark.specimen.outer_fillet_radius_m.to_bits(),
        THORNE_2026_STEEL_DISC_FILLET_RADIUS_M.to_bits()
    );
    assert!(matches!(
        benchmark.specimen.profile_spec(),
        DiscProfileSpec::SolidCylinder {
            outer_radius_m,
            thickness_m,
            edge_treatment: SquatDiscEdgeTreatment::CircularFillet { radius },
        } if outer_radius_m.to_bits() == (0.5 * THORNE_2026_STEEL_DISC_DIAMETER_M).to_bits()
            && thickness_m.to_bits() == THORNE_2026_STEEL_DISC_THICKNESS_M.to_bits()
            && radius.to_bits() == THORNE_2026_STEEL_DISC_FILLET_RADIUS_M.to_bits()
    ));

    let run = run_thorne_2026_steel_glass_benchmark(&benchmark).expect("benchmark run");
    let initial = run.samples.first().expect("initial sample");
    let radius_m = 0.5 * THORNE_2026_STEEL_DISC_DIAMETER_M;
    let expected_rolling_w = THORNE_2026_FITTED_ROLLING_COEFFICIENT
        * THORNE_2026_STEEL_DISC_MASS_KG
        * STANDARD_GRAVITY_M_PER_S2
        * radius_m
        * initial.theta_rad.cos()
        * initial.omega_rad_s;
    assert_eq!(
        initial.powers.published_rolling_w.to_bits(),
        expected_rolling_w.to_bits()
    );
    assert_ne!(
        initial.powers.published_rolling_w.to_bits(),
        (expected_rolling_w / initial.theta_rad.cos()).to_bits(),
        "the source-bound Eq. 5 channel must retain its explicit cos(theta) factor"
    );
    let expected_air_w = 4.0
        * STANDARD_GRAVITY_M_PER_S2.powf(1.25)
        * radius_m.powf(2.75)
        * (THORNE_2026_DECLARED_AIR_VISCOSITY_PA_S * THORNE_2026_AMBIENT_AIR_DENSITY_KG_PER_M3)
            .sqrt()
        * initial.theta_rad.powf(-1.25);
    assert!((initial.powers.bildsten_boundary_layer_w / expected_air_w - 1.0).abs() < 1.0e-12);
    assert_eq!(initial.powers.dry_contour_w, 0.0);
    assert_eq!(
        run.provenance.model_authority,
        "literature-calibrated-analytical"
    );
    assert_eq!(
        run.provenance.physical_validation,
        "no-raw-trajectory-or-full-fsi-validation-claimed"
    );
    assert_eq!(run.terminal, ReducedDecayTerminal::ValidityCutoff);
    assert!(
        run.final_sample().expect("cutoff sample").time_s > 8.0,
        "the source run must cover the complete eight-second presentation horizon"
    );
}

#[test]
fn e2e_thorne_2026_crossover_and_ambient_vacuum_direction_match_the_declared_model() {
    let ambient = Thorne2026SteelGlassBenchmark::ambient().expect("ambient benchmark");
    let vacuum = Thorne2026SteelGlassBenchmark::partial_vacuum().expect("vacuum benchmark");
    let crossover =
        thorne_2026_channel_crossover_diagnostic(&ambient).expect("published-channel crossover");
    assert!(matches!(
        crossover,
        ChannelCrossoverDiagnostic::AtInclination { theta_rad }
            if (theta_rad - 0.03).abs() < 0.01
    ));

    let ambient_run =
        run_thorne_2026_steel_glass_benchmark(&ambient).expect("ambient benchmark run");
    let vacuum_run = run_thorne_2026_steel_glass_benchmark(&vacuum).expect("vacuum benchmark run");
    let ambient_initial = ambient_run.samples.first().expect("ambient initial sample");
    let vacuum_initial = vacuum_run.samples.first().expect("vacuum initial sample");
    assert_eq!(
        ambient_initial.powers.published_rolling_w.to_bits(),
        vacuum_initial.powers.published_rolling_w.to_bits()
    );
    assert!(
        ambient_initial.powers.bildsten_boundary_layer_w
            > vacuum_initial.powers.bildsten_boundary_layer_w
    );
    assert!(
        ambient_run.final_sample().expect("ambient final").time_s
            < vacuum_run.final_sample().expect("vacuum final").time_s,
        "lower gas density must lengthen the declared model's decay from the same initial state"
    );
    let density_ratio =
        THORNE_2026_AMBIENT_AIR_DENSITY_KG_PER_M3 / THORNE_2026_VACUUM_AIR_DENSITY_KG_PER_M3;
    let power_ratio = ambient_initial.powers.bildsten_boundary_layer_w
        / vacuum_initial.powers.bildsten_boundary_layer_w;
    assert!((power_ratio - density_ratio.sqrt()).abs() < 1.0e-12);
}

#[test]
fn e2e_thorne_2026_refinement_preserves_both_published_channels() {
    let benchmark = Thorne2026SteelGlassBenchmark::ambient().expect("ambient benchmark");
    let evidence = thorne_2026_refinement_evidence(&benchmark).expect("benchmark refinement");
    assert_eq!(
        evidence.coarse.terminal,
        ReducedDecayTerminal::ValidityCutoff
    );
    assert_eq!(evidence.fine.terminal, ReducedDecayTerminal::ValidityCutoff);
    assert_eq!(
        evidence.fine.parameters.timestep_s.to_bits(),
        (0.5 * evidence.coarse.parameters.timestep_s).to_bits()
    );
    let coarse = evidence.coarse.final_sample().expect("coarse cutoff");
    let fine = evidence.fine.final_sample().expect("fine cutoff");
    assert!(coarse.work.published_rolling_j > 0.0);
    assert!(coarse.work.bildsten_boundary_layer_j > 0.0);
    assert!(fine.work.published_rolling_j > 0.0);
    assert!(fine.work.bildsten_boundary_layer_j > 0.0);
    assert!(evidence.terminal_time_difference_s < 1.0e-4);
    assert!(evidence.total_work_difference_j < 1.0e-12);
    assert!((coarse.work.published_rolling_j - fine.work.published_rolling_j).abs() < 3.0e-7);
    assert!(
        (coarse.work.bildsten_boundary_layer_j - fine.work.bildsten_boundary_layer_j).abs()
            < 3.0e-7
    );
    assert_eq!(
        evidence.coarse.provenance.model_authority,
        "literature-calibrated-analytical"
    );
    assert_eq!(
        evidence
            .fine
            .provenance
            .published_rolling_source_id
            .as_deref(),
        Some(THORNE_2026_SOURCE_ID)
    );
}
