//! Numerical self-consistency checks for the public reduced-decay API.
//!
//! Exponent recovery and `P * dt` energy closure are internal consistency
//! checks of the stated reduced equations. They are not independent emergence,
//! experimental agreement, or physical validation.

use fs_euler_disc_e2e::{
    BildstenBoundaryLayerChannel, DryContourChannel, ReducedDecayError, ReducedDecayInput,
    ReducedDecayTerminal, STANDARD_GRAVITY_M_PER_S2, refinement_evidence, run_reduced_decay,
    structured_runner_output,
};
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
        first.starts_with("schema=reduced-decay-v1 model_id=euler-disc-small-angle-late-stage-v1 ")
    );
    assert!(first.contains("model_authority=numerical-reference-only"));
    assert!(first.contains("physical_validation=not-claimed"));
    assert!(first.contains("dry_authority=SyntheticFixture"));
    assert!(first.contains("bildsten_source_id=synthetic/bildsten-energy-only"));
    assert!(first.contains("terminal=ValidityCutoff"));
    assert!(first.contains("dry_work_j="));
    assert!(first.contains("bildsten_work_j="));
    assert!(first.contains("closure_residual_j="));
    assert!(first.contains("\nrefinement_terminal_time_difference_s="));
    assert!(first.ends_with("evidence_scope=numerical-self-consistency-only"));
}
