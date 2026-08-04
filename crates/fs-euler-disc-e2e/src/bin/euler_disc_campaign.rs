//! Deterministic production JSONL campaign for the committed Euler-disc rungs.
//!
//! Contact and reduced-decay records are deliberately registered only once the
//! owning production modules have landed on `main`; this runner never compiles
//! uncommitted source through a path inclusion.

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_euler_disc_e2e::convergence::{
    CalibrationEvidenceKind, CalibrationReadinessError, CalibrationReadinessInput,
    CensorAwareDurationOrdering, CensorAwareRankingRefusal, ConvergenceScales, DeclaredEvidence,
    HorizonContinuationPolicy, ObservedOrder, RefinementMode, RunOutcome, ThreeRungConvergence,
    admit_calibration_readiness, analyse_three_rung_convergence, classify_outcome,
    compare_censor_aware_durations, missing_calibration_evidence,
};
use fs_euler_disc_e2e::coupled_runner::{
    CoupledChannelFactors, CoupledControls, CoupledInitialState, CoupledNumericalRefusalReason,
    CoupledRun, run_closed_profile_reduced,
};
use fs_euler_disc_e2e::specimen::{DiscProfileSpec, ResolvedDiscProfile};
use fs_euler_disc_e2e::{
    BaseGeometryScope, BaseResponseInput, BaselineRunOutput, ChannelCrossoverDiagnostic,
    ChannelCrossoverNotComparable, ContactDiscGeometry, ContactDynamicsInput, ContactLoadScope,
    ContactTermination, LevelSupportInput, MovingContactLoad, ProfileContactDynamicsInput,
    SquatDiscInput, channel_crossover_diagnostic, refine_profile_timestep_by_two,
    refine_reduced_base_response, refinement_evidence, run_ideal_conservative_baseline,
    run_profile_contact_dynamics, run_reduced_decay, small_angle_rolling_profile_initializer,
};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_flux::{
    ApplicabilityEnvelope, BodyKinematics, ClosedRange, ContributionFamily, CorrelationIdentity,
    CorrelationUncertainty, DiscGeometry, DiscPose, EdgeFlow, FormDrag, GasProperties,
    GasPropertyCard, ReducedAeroComponents, ReducedAeroInput, ReducedAeroModel,
    RotationalSkinFriction, SurfaceRoughness, Vec3, WorkWindow,
};
use fs_mbd::Vec3 as MbdVec3;
use fs_rep_frep::{AxisymmetricChart, SquatDiscEdgeTreatment};
use fs_solid::{
    AssemblyBudget, DampingModel, ShellIdentity, ShellMaterial, ShellNode, ShellPlate, ShellSupport,
};
use fs_tribo::{InputAuthority, InterfaceMedium, InterfaceSystemRef};

const SCHEMA: &str = "euler-disc-campaign-jsonl-v1";
const SI_UNITS: &str = "SI:m,kg,s,N,J,Pa,rad";
const NO_PHYSICAL_VALIDATION: &str =
    "numerical-slice-only:no-physical-validation-or-target-ranking";
const DECLARED_CX_BUDGET_NO_CLAIM: &str =
    "declared-cx-contract-metadata-not-enforced-quota;hard-model-step-bounds-only";
// Exact decimal identity for `0x4555_4c45_5243_414d`; retain the decimal form
// because it is emitted in replay records.
const CAMPAIGN_SEED: u64 = 4_995_983_222_254_027_085;
const CAMPAIGN_SEED_DEC: &str = "4995983222254027085";
const EXECUTION_POLL_QUOTA: u32 = 10_000;
const EXECUTION_COST_QUOTA: u64 = 100_000;

fn main() {
    let result = match std::env::args().nth(1).as_deref() {
        None => run_campaign().and_then(|_| run_closed_campaign()),
        Some("--closed-only") => run_closed_campaign(),
        Some(_) => Err("usage: euler_disc_campaign [--closed-only]".to_string()),
    };
    if let Err(error) = result {
        eprintln!("euler-disc campaign refusal: {error}");
        std::process::exit(1);
    }
}

fn run_closed_campaign() -> Result<(), String> {
    with_closed_context(run_closed_campaign_with_context)
}

#[derive(Clone, Copy)]
struct ClosedCase {
    name: &'static str,
    profile_kind: &'static str,
    profile: DiscProfileSpec,
    density_kg_per_m3: f64,
    base_stiffness_scale: f64,
    rolling_resistance_m: f64,
    gas_rotational_damping_n_m_s: f64,
    gas_translation_damping_n_s_per_m: f64,
}

struct ContinuedRun {
    run: CoupledRun,
    last_sample: fs_euler_disc_e2e::coupled_runner::CoupledSample,
    declared_horizon_s: f64,
    continuation_count: u32,
}

fn run_closed_campaign_with_context(cx: &Cx<'_>) -> Result<(), String> {
    const RADIUS: f64 = 0.038;
    const THICKNESS: f64 = 0.006;
    const INNER_RATIO: f64 = 0.65;
    const DENSITY: f64 = 2_680.0;
    let sharp = DiscProfileSpec::SolidCylinder {
        outer_radius_m: RADIUS,
        thickness_m: THICKNESS,
        edge_treatment: SquatDiscEdgeTreatment::Sharp,
    };
    let cases = [
        ClosedCase {
            name: "solid",
            profile_kind: "solid-cylinder-sharp",
            profile: sharp,
            density_kg_per_m3: DENSITY,
            base_stiffness_scale: 1.0,
            rolling_resistance_m: 4.0e-5,
            gas_rotational_damping_n_m_s: 2.0e-7,
            gas_translation_damping_n_s_per_m: 4.0e-4,
        },
        ClosedCase {
            name: "solid-fillet-1mm",
            profile_kind: "solid-cylinder-circular-fillet",
            profile: DiscProfileSpec::SolidCylinder {
                outer_radius_m: RADIUS,
                thickness_m: THICKNESS,
                edge_treatment: SquatDiscEdgeTreatment::CircularFillet { radius: 0.001 },
            },
            density_kg_per_m3: DENSITY,
            ..cases_defaults(sharp)
        },
        ClosedCase {
            name: "solid-chamfer-1mm",
            profile_kind: "solid-cylinder-linear-chamfer",
            profile: DiscProfileSpec::ChamferedCylinder {
                outer_radius_m: RADIUS,
                thickness_m: THICKNESS,
                chamfer_radial_m: 0.001,
                chamfer_axial_m: 0.001,
            },
            density_kg_per_m3: DENSITY,
            ..cases_defaults(sharp)
        },
        ClosedCase {
            name: "ring-fixed-density",
            profile_kind: "annular-cylinder",
            profile: DiscProfileSpec::AnnularCylinder {
                outer_radius_m: RADIUS,
                inner_radius_m: INNER_RATIO * RADIUS,
                thickness_m: THICKNESS,
            },
            density_kg_per_m3: DENSITY,
            ..cases_defaults(sharp)
        },
        ClosedCase {
            name: "ring-equal-mass",
            profile_kind: "annular-cylinder-equal-mass-control",
            profile: DiscProfileSpec::AnnularCylinder {
                outer_radius_m: RADIUS,
                inner_radius_m: INNER_RATIO * RADIUS,
                thickness_m: THICKNESS,
            },
            density_kg_per_m3: DENSITY / (1.0 - INNER_RATIO * INNER_RATIO),
            ..cases_defaults(sharp)
        },
        ClosedCase {
            name: "symmetric-tapered",
            profile_kind: "symmetric-double-frustum",
            profile: DiscProfileSpec::SymmetricTapered {
                outer_radius_m: RADIUS,
                face_radius_m: 0.012,
                thickness_m: THICKNESS,
            },
            density_kg_per_m3: DENSITY,
            ..cases_defaults(sharp)
        },
        ClosedCase {
            name: "large-rim",
            profile_kind: "solid-cylinder-sharp",
            profile: DiscProfileSpec::SolidCylinder {
                outer_radius_m: 0.052,
                thickness_m: THICKNESS,
                edge_treatment: SquatDiscEdgeTreatment::Sharp,
            },
            density_kg_per_m3: DENSITY,
            ..cases_defaults(sharp)
        },
        ClosedCase {
            name: "dense-material",
            profile_kind: "solid-cylinder-sharp",
            profile: sharp,
            density_kg_per_m3: 7_850.0,
            ..cases_defaults(sharp)
        },
        ClosedCase {
            name: "compliant-base",
            profile_kind: "solid-cylinder-sharp",
            profile: sharp,
            density_kg_per_m3: DENSITY,
            base_stiffness_scale: 0.35,
            ..cases_defaults(sharp)
        },
        ClosedCase {
            name: "solid-no-gas",
            profile_kind: "solid-cylinder-sharp",
            profile: sharp,
            density_kg_per_m3: DENSITY,
            gas_rotational_damping_n_m_s: 0.0,
            gas_translation_damping_n_s_per_m: 0.0,
            ..cases_defaults(sharp)
        },
        ClosedCase {
            name: "solid-no-rolling",
            profile_kind: "solid-cylinder-sharp",
            profile: sharp,
            density_kg_per_m3: DENSITY,
            rolling_resistance_m: 0.0,
            ..cases_defaults(sharp)
        },
    ];
    let controls = CoupledControls {
        timestep_s: 2.0e-5,
        maximum_steps: 100_000,
        terminal_inclination_rad: 0.002,
        reimpact_limit: 128,
    };
    let initial = CoupledInitialState {
        inclination_rad: 0.08,
        precession_rad_per_s: 16.0,
        spin_rad_per_s: 120.0,
    };
    let continuation = HorizonContinuationPolicy {
        initial_horizon_s: 2.0,
        maximum_horizon_s: 8.0,
        multiplier: 2.0,
        maximum_extensions: 2,
    };
    let mut records = Vec::new();
    let mut resolved_solid = None;
    let mut resolved_equal_mass_ring = None;
    for case in cases {
        let profile = case
            .profile
            .resolve(case.density_kg_per_m3, cx)
            .map_err(|error| format!("{} profile: {error}", case.name))?;
        let channels = channels_for(case);
        let continued =
            run_with_declared_continuation(&profile, channels, controls, initial, continuation, cx)
                .map_err(|error| format!("{}: {error}", case.name))?;
        let outcome = classify_outcome(&continued.run)
            .map_err(|error| format!("{} outcome: {error}", case.name))?;
        records.push(closed_case_record(
            case, &profile, controls, initial, &continued, outcome,
        )?);
        if case.name == "solid" {
            resolved_solid = Some((profile, channels));
        } else if case.name == "ring-equal-mass" {
            resolved_equal_mass_ring = Some((profile, channels));
        }
    }
    let (solid_profile, solid_channels) =
        resolved_solid.ok_or_else(|| "solid convergence anchor missing".to_owned())?;
    let (ring_profile, ring_channels) = resolved_equal_mass_ring
        .ok_or_else(|| "equal-mass ring convergence anchor missing".to_owned())?;
    records.push(convergence_record(
        &solid_profile,
        solid_channels,
        initial,
        cx,
    )?);
    records.push(ranking_convergence_record(
        RankingRefinementInput {
            solid_profile: &solid_profile,
            solid_channels,
            ring_profile: &ring_profile,
            ring_channels,
            initial,
        },
        cx,
    )?);
    records.push(calibration_readiness_record()?);

    let payload = records.join("\n");
    let digest = fs_blake3::hash_domain(
        "org.frankensim.euler-disc-campaign-jsonl.v3",
        payload.as_bytes(),
    );
    let record_count = records.len();
    for record in records {
        println!("{record}");
    }
    println!(
        "{{\"schema\":\"euler-disc-campaign-jsonl-v3\",\"scenario\":\"campaign-complete\",\"model\":\"closed-time-evolving-profile-native-reduced-euler-disc\",\"authority\":\"integration-local\",\"units\":\"SI:m,kg,s,N,J,rad\",\"terminal\":\"completed\",\"record_count\":{},\"digest_blake3\":\"{}\",\"no_claim\":\"no-physical-validation-or-video-ranking;right-censored-cases-are-not-durations\"}}",
        record_count,
        digest.to_hex()
    );
    Ok(())
}

const fn cases_defaults(profile: DiscProfileSpec) -> ClosedCase {
    ClosedCase {
        name: "unused-default-name",
        profile_kind: "unused-default-profile-kind",
        profile,
        density_kg_per_m3: 2_680.0,
        base_stiffness_scale: 1.0,
        rolling_resistance_m: 4.0e-5,
        gas_rotational_damping_n_m_s: 2.0e-7,
        gas_translation_damping_n_s_per_m: 4.0e-4,
    }
}

fn channels_for(case: ClosedCase) -> CoupledChannelFactors {
    CoupledChannelFactors {
        gravity_m_per_s2: 9.806_65,
        sliding_friction_coefficient: 0.42,
        rolling_resistance_m: case.rolling_resistance_m,
        contact_stiffness_n_per_m: 8.0e4,
        contact_damping_n_s_per_m: 3.0,
        base_effective_mass_kg: 0.25,
        base_stiffness_n_per_m: 4.0e4 * case.base_stiffness_scale,
        base_damping_n_s_per_m: 4.0,
        gas_rotational_damping_n_m_s: case.gas_rotational_damping_n_m_s,
        gas_translation_damping_n_s_per_m: case.gas_translation_damping_n_s_per_m,
    }
}

fn run_with_declared_continuation(
    profile: &ResolvedDiscProfile,
    channels: CoupledChannelFactors,
    controls: CoupledControls,
    initial: CoupledInitialState,
    policy: HorizonContinuationPolicy,
    cx: &Cx<'_>,
) -> Result<ContinuedRun, String> {
    policy
        .validate()
        .map_err(|error| format!("invalid continuation policy: {error}"))?;
    let mut declared_horizon_s = policy.initial_horizon_s;
    let mut continuation_count = 0;
    let mut restart = None;
    let mut last_sample = None;
    loop {
        let completed_time_s = restart.as_ref().map_or(
            0.0,
            |checkpoint: &fs_euler_disc_e2e::coupled_runner::CoupledCheckpoint| checkpoint.time_s,
        );
        let remaining_s = declared_horizon_s - completed_time_s;
        if !remaining_s.is_finite() || remaining_s <= 0.0 {
            return Err("continuation horizon did not advance".to_owned());
        }
        let step_count_f64 = (remaining_s / controls.timestep_s).round();
        if !(1.0..=f64::from(u32::MAX)).contains(&step_count_f64) {
            return Err("continuation step count exceeded the runner budget".to_owned());
        }
        let segment_controls = CoupledControls {
            maximum_steps: step_count_f64 as u32,
            ..controls
        };
        let run =
            run_closed_profile_reduced(profile, channels, segment_controls, initial, restart, cx)
                .map_err(|error| error.to_string())?;
        if let Some(sample) = run.samples.last() {
            last_sample = Some(sample.clone());
        }
        let outcome = classify_outcome(&run).map_err(|error| error.to_string())?;
        let next = policy
            .next_horizon_s(outcome, declared_horizon_s, continuation_count)
            .map_err(|error| error.to_string())?;
        let Some(next_horizon_s) = next else {
            return Ok(ContinuedRun {
                run,
                last_sample: last_sample
                    .ok_or_else(|| "trajectory retained no committed sample".to_owned())?,
                declared_horizon_s,
                continuation_count,
            });
        };
        restart = Some(run.checkpoint);
        declared_horizon_s = next_horizon_s;
        continuation_count += 1;
    }
}

fn closed_case_record(
    case: ClosedCase,
    profile: &ResolvedDiscProfile,
    controls: CoupledControls,
    initial: CoupledInitialState,
    continued: &ContinuedRun,
    outcome: RunOutcome,
) -> Result<String, String> {
    let run = &continued.run;
    let last = &continued.last_sample;
    let (outcome_kind, observed_physical_terminal, retained_time_s, refusal_reason) = match outcome
    {
        RunOutcome::PhysicalTerminal { event_time_s, .. } => {
            ("physical-terminal-inclination", true, event_time_s, "none")
        }
        RunOutcome::RightCensored { censor_time_s } => {
            ("right-censored", false, censor_time_s, "none")
        }
        RunOutcome::NumericalRefusal {
            last_valid_time_s,
            reason,
        } => (
            "numerical-refusal",
            false,
            last_valid_time_s,
            numerical_refusal_reason_name(reason),
        ),
    };
    let properties = profile.mass_properties;
    let channels = channels_for(case);
    let support_feature = last
        .support_source_feature
        .ok_or_else(|| format!("{}: profile support feature was not retained", case.name))?;
    Ok(format!(
        "{{\"schema\":\"euler-disc-campaign-jsonl-v3\",\"scenario\":\"closed-reduced-{}\",\"model\":\"fs-mbd-profile-native-reduced-coupled-runner\",\"authority\":\"numerical-slice-only\",\"units\":\"SI:m,kg,s,N,J,rad\",\"profile\":{{\"kind\":\"{}\",\"identity_u64_dec\":\"{}\",\"support_source_feature\":{},\"mass_and_support_same_chart\":true}},\"inputs\":{{\"input_units\":\"SI:kg,m,s,N,J,rad\",\"mass_kg\":{:.17e},\"radius_m\":{:.17e},\"thickness_m\":{:.17e},\"density_kg_m3\":{:.17e},\"gravity_m_per_s2\":{:.17e},\"transverse_inertia_kg_m2\":{:.17e},\"axial_inertia_kg_m2\":{:.17e},\"timestep_s\":{:.17e},\"initial_segment_maximum_steps\":{},\"maximum_steps\":{},\"initial_horizon_s\":{:.17e},\"maximum_horizon_s\":{:.17e},\"declared_final_horizon_s\":{:.17e},\"continuation_count\":{},\"terminal_inclination_rad\":{:.17e},\"reimpact_limit\":{},\"initial_inclination_rad\":{:.17e},\"initial_precession_rad_per_s\":{:.17e},\"initial_spin_rad_per_s\":{:.17e},\"sliding_friction_coefficient\":{:.17e},\"rolling_resistance_m\":{:.17e},\"base_effective_mass_kg\":{:.17e},\"base_stiffness_n_per_m\":{:.17e},\"base_damping_n_s_per_m\":{:.17e},\"contact_stiffness_n_per_m\":{:.17e},\"contact_damping_n_s_per_m\":{:.17e},\"gas_rotational_damping_n_m_s\":{:.17e},\"gas_translation_damping_n_s_per_m\":{:.17e}}},\"terminal\":\"{}\",\"outcome\":{{\"kind\":\"{}\",\"observed_physical_terminal\":{},\"retained_time_s\":{:.17e},\"numerical_refusal_reason\":\"{}\"}},\"qoi\":{{\"retained_time_s\":{:.17e},\"inclination_rad\":{:.17e},\"precession_rad_per_s\":{:.17e},\"spin_rad_per_s\":{:.17e},\"precession_acceleration_rad_per_s2\":{:.17e},\"reimpact_count\":{}}},\"channel_work_j\":{{\"gravity\":{:.17e},\"contact\":{:.17e},\"rolling\":{:.17e},\"base\":{:.17e},\"gas\":{:.17e}}},\"last_step_channel_work_j\":{{\"gravity\":{:.17e},\"contact\":{:.17e},\"rolling\":{:.17e},\"base\":{:.17e},\"gas\":{:.17e}}},\"energy\":{{\"initial_total_j\":{:.17e},\"final_total_j\":{:.17e},\"defect_j\":{:.17e},\"relative_defect\":{:.17e}}},\"applicability\":\"{}\",\"model_disagreement\":\"{}\",\"no_claim\":\"no-physical-validation-or-video-ranking;right-censored-time-is-not-duration\"}}",
        case.name,
        case.profile_kind,
        profile.identity.0,
        support_feature,
        properties.mass,
        profile.dimensions.outer_radius_m,
        profile.dimensions.thickness_m,
        profile.density_kg_per_m3,
        channels.gravity_m_per_s2,
        properties.principal_inertia.transverse,
        properties.principal_inertia.axial,
        controls.timestep_s,
        controls.maximum_steps,
        (8.0 / controls.timestep_s).round() as u32,
        2.0,
        8.0,
        continued.declared_horizon_s,
        continued.continuation_count,
        controls.terminal_inclination_rad,
        controls.reimpact_limit,
        initial.inclination_rad,
        initial.precession_rad_per_s,
        initial.spin_rad_per_s,
        channels.sliding_friction_coefficient,
        channels.rolling_resistance_m,
        channels.base_effective_mass_kg,
        channels.base_stiffness_n_per_m,
        channels.base_damping_n_s_per_m,
        channels.contact_stiffness_n_per_m,
        channels.contact_damping_n_s_per_m,
        channels.gas_rotational_damping_n_m_s,
        channels.gas_translation_damping_n_s_per_m,
        outcome_kind,
        outcome_kind,
        observed_physical_terminal,
        retained_time_s,
        refusal_reason,
        last.time_s,
        last.inclination_rad,
        last.precession_rad_per_s,
        last.spin_rad_per_s,
        last.precession_acceleration_rad_per_s2,
        run.checkpoint.reimpact_count,
        run.checkpoint.accumulated_channel_work_j[0],
        run.checkpoint.accumulated_channel_work_j[1],
        run.checkpoint.accumulated_channel_work_j[2],
        run.checkpoint.accumulated_channel_work_j[3],
        run.checkpoint.accumulated_channel_work_j[4],
        last.channels.gravity.work_j,
        last.channels.contact.work_j,
        last.channels.rolling.work_j,
        last.channels.base.work_j,
        last.channels.gas.work_j,
        run.checkpoint.initial_total_energy_j,
        last.mechanical_energy_j,
        last.energy_defect_j,
        last.energy_defect_j.abs()
            / run
                .checkpoint
                .initial_total_energy_j
                .abs()
                .max(f64::MIN_POSITIVE),
        run.applicability,
        run.model_disagreement,
    ))
}

fn convergence_record(
    profile: &ResolvedDiscProfile,
    channels: CoupledChannelFactors,
    initial: CoupledInitialState,
    cx: &Cx<'_>,
) -> Result<String, String> {
    const HORIZON_S: f64 = 0.2;
    let run = |timestep_s: f64| {
        let steps = (HORIZON_S / timestep_s).round() as u32;
        run_closed_profile_reduced(
            profile,
            channels,
            CoupledControls {
                timestep_s,
                maximum_steps: steps,
                terminal_inclination_rad: 0.002,
                reimpact_limit: 128,
            },
            initial,
            None,
            cx,
        )
        .map_err(|error| error.to_string())
    };
    let coarse = run(4.0e-5)?;
    let fine = run(2.0e-5)?;
    let reference = run(1.0e-5)?;
    let energy_scale = reference
        .checkpoint
        .initial_total_energy_j
        .abs()
        .max(f64::MIN_POSITIVE);
    let receipt = analyse_three_rung_convergence(ThreeRungConvergence {
        coarse: &coarse,
        fine: &fine,
        reference: &reference,
        coarse_timestep_s: 4.0e-5,
        fine_timestep_s: 2.0e-5,
        reference_timestep_s: 1.0e-5,
        mode: RefinementMode::Eventful {
            reason: "unilateral contact and support-feature switching".to_owned(),
        },
        scales: ConvergenceScales {
            inclination_rad: initial.inclination_rad,
            precession_rad_per_s: initial.precession_rad_per_s.abs().max(1.0),
            spin_rad_per_s: initial.spin_rad_per_s.abs().max(1.0),
            work_j: energy_scale,
            energy_j: energy_scale,
        },
    })
    .map_err(|error| error.to_string())?;
    let fine_reference_qoi_linf = receipt
        .fine_reference_qoi
        .inclination
        .max(receipt.fine_reference_qoi.precession)
        .max(receipt.fine_reference_qoi.spin);
    let fine_reference_work_linf = receipt
        .fine_reference_work_energy
        .channel_work
        .into_iter()
        .fold(0.0_f64, f64::max);
    let declared_delta_limit = 5.0e-3;
    let within_declared_delta_band = fine_reference_qoi_linf <= declared_delta_limit
        && fine_reference_work_linf <= declared_delta_limit
        && receipt.fine_reference_work_energy.energy_defect <= declared_delta_limit;
    let order_status = match receipt.observed_order {
        ObservedOrder::Available { .. } => "available",
        ObservedOrder::NotApplicable { .. } => "withheld-eventful-mode",
    };
    Ok(format!(
        "{{\"schema\":\"euler-disc-campaign-jsonl-v3\",\"scenario\":\"closed-reduced-solid-timestep-convergence\",\"model\":\"h-h2-h4-fixed-horizon-profile-native-runner\",\"authority\":\"numerical-convergence-evidence-only\",\"units\":\"dimensionless-normalized-deltas\",\"terminal\":\"analysis-complete\",\"horizon_s\":{HORIZON_S:.17e},\"timesteps_s\":{{\"h\":{:.17e},\"h2\":{:.17e},\"h4\":{:.17e}}},\"terminal_class_agreement\":{},\"fine_reference_qoi\":{{\"inclination\":{:.17e},\"precession\":{:.17e},\"spin\":{:.17e}}},\"fine_reference_qoi_linf\":{fine_reference_qoi_linf:.17e},\"fine_reference_work_linf\":{fine_reference_work_linf:.17e},\"fine_reference_energy_defect\":{:.17e},\"declared_normalized_delta_limit\":{declared_delta_limit:.17e},\"within_declared_delta_band\":{},\"observed_order\":\"{}\",\"no_claim\":\"fixed-horizon-numerical-sensitivity-only;eventful-mode-withholds-smooth-order;not-terminal-time-or-physical-validation\"}}",
        4.0e-5,
        2.0e-5,
        1.0e-5,
        receipt.terminal_class_agreement,
        receipt.fine_reference_qoi.inclination,
        receipt.fine_reference_qoi.precession,
        receipt.fine_reference_qoi.spin,
        receipt.fine_reference_work_energy.energy_defect,
        within_declared_delta_band,
        order_status,
    ))
}

#[derive(Clone, Copy)]
struct RankingRefinementInput<'profile> {
    solid_profile: &'profile ResolvedDiscProfile,
    solid_channels: CoupledChannelFactors,
    ring_profile: &'profile ResolvedDiscProfile,
    ring_channels: CoupledChannelFactors,
    initial: CoupledInitialState,
}

#[derive(Clone, Copy)]
struct RankingRung {
    timestep_s: f64,
    solid: RunOutcome,
    ring: RunOutcome,
    ordering: Result<CensorAwareDurationOrdering, CensorAwareRankingRefusal>,
}

/// Refines the actual ranking claim, rather than inferring it from a
/// phase-sensitive endpoint state. Every rung uses the same two-second
/// observation window. A ring event before the solid censor bound can prove a
/// strict ordering without inventing a terminal time for the solid.
fn ranking_convergence_record(
    input: RankingRefinementInput<'_>,
    cx: &Cx<'_>,
) -> Result<String, String> {
    const HORIZON_S: f64 = 2.0;
    const TIMESTEPS_S: [f64; 3] = [8.0e-5, 4.0e-5, 2.0e-5];
    const EVENT_TIME_RELATIVE_LIMIT: f64 = 5.0e-3;

    let mut rungs = Vec::with_capacity(TIMESTEPS_S.len());
    for timestep_s in TIMESTEPS_S {
        let maximum_steps = (HORIZON_S / timestep_s).round() as u32;
        let controls = CoupledControls {
            timestep_s,
            maximum_steps,
            terminal_inclination_rad: 0.002,
            reimpact_limit: 128,
        };
        let solid_run = run_closed_profile_reduced(
            input.solid_profile,
            input.solid_channels,
            controls,
            input.initial,
            None,
            cx,
        )
        .map_err(|error| format!("solid ranking refinement at dt={timestep_s}: {error}"))?;
        let ring_run = run_closed_profile_reduced(
            input.ring_profile,
            input.ring_channels,
            controls,
            input.initial,
            None,
            cx,
        )
        .map_err(|error| format!("ring ranking refinement at dt={timestep_s}: {error}"))?;
        let solid = classify_outcome(&solid_run)
            .map_err(|error| format!("solid ranking outcome at dt={timestep_s}: {error}"))?;
        let ring = classify_outcome(&ring_run)
            .map_err(|error| format!("ring ranking outcome at dt={timestep_s}: {error}"))?;
        let ordering = compare_censor_aware_durations(ring, solid);
        rungs.push(RankingRung {
            timestep_s,
            solid,
            ring,
            ordering,
        });
    }
    let [coarse, fine, reference] = rungs.as_slice() else {
        return Err("ranking refinement did not retain exactly three rungs".to_owned());
    };
    let coarse_fine_event_delta = event_time_delta(coarse.ring, fine.ring);
    let fine_reference_event_delta = event_time_delta(fine.ring, reference.ring);
    let ordering_agreement =
        coarse.ordering == fine.ordering && fine.ordering == reference.ordering;
    let ring_shorter_bound_proven = ordering_agreement
        && reference.ordering == Ok(CensorAwareDurationOrdering::ProvenLeftShorter);
    let event_time_within_declared_band = fine_reference_event_delta
        .map(|(_, relative)| relative <= EVENT_TIME_RELATIVE_LIMIT)
        .unwrap_or(false);
    let ranking_numerically_supported =
        ring_shorter_bound_proven && event_time_within_declared_band;

    Ok(format!(
        "{{\"schema\":\"euler-disc-campaign-jsonl-v3\",\"scenario\":\"equal-mass-ring-vs-solid-ranking-convergence\",\"model\":\"three-rung-common-window-censor-aware-profile-native-runner\",\"authority\":\"numerical-convergence-evidence-only\",\"units\":\"SI:s\",\"terminal\":\"analysis-complete\",\"common_horizon_s\":{HORIZON_S:.17e},\"timesteps_s\":{{\"h\":{:.17e},\"h2\":{:.17e},\"h4\":{:.17e}}},\"rungs\":{{\"h\":{},\"h2\":{},\"h4\":{}}},\"ordering_agreement\":{},\"ring_shorter_than_solid_bound_proven_at_all_rungs\":{},\"ring_event_time_coarse_fine\":{},\"ring_event_time_fine_reference\":{},\"declared_event_time_relative_limit\":{EVENT_TIME_RELATIVE_LIMIT:.17e},\"ring_event_time_within_declared_band\":{},\"ranking_numerically_supported\":{},\"no_claim\":\"numerical-ranking-of-declared-reduced-model-only;solid-censor-is-a-lower-bound-not-a-duration;not-experimental-calibration-or-video-validation\"}}",
        coarse.timestep_s,
        fine.timestep_s,
        reference.timestep_s,
        ranking_rung_json(*coarse),
        ranking_rung_json(*fine),
        ranking_rung_json(*reference),
        ordering_agreement,
        ring_shorter_bound_proven,
        event_delta_json(coarse_fine_event_delta),
        event_delta_json(fine_reference_event_delta),
        event_time_within_declared_band,
        ranking_numerically_supported,
    ))
}

fn ranking_rung_json(rung: RankingRung) -> String {
    let (ring_kind, ring_time_s, ring_refusal_reason) = outcome_kind_and_time(rung.ring);
    let (solid_kind, solid_time_s, solid_refusal_reason) = outcome_kind_and_time(rung.solid);
    format!(
        "{{\"ring_outcome\":\"{ring_kind}\",\"ring_retained_time_s\":{ring_time_s:.17e},\"ring_numerical_refusal_reason\":\"{ring_refusal_reason}\",\"solid_outcome\":\"{solid_kind}\",\"solid_retained_time_s\":{solid_time_s:.17e},\"solid_numerical_refusal_reason\":\"{solid_refusal_reason}\",\"censor_aware_ordering\":\"{}\"}}",
        censor_ordering_name(rung.ordering),
    )
}

fn outcome_kind_and_time(outcome: RunOutcome) -> (&'static str, f64, &'static str) {
    match outcome {
        RunOutcome::PhysicalTerminal { event_time_s, .. } => {
            ("physical-terminal-inclination", event_time_s, "none")
        }
        RunOutcome::RightCensored { censor_time_s } => ("right-censored", censor_time_s, "none"),
        RunOutcome::NumericalRefusal {
            last_valid_time_s,
            reason,
        } => (
            "numerical-refusal",
            last_valid_time_s,
            numerical_refusal_reason_name(reason),
        ),
    }
}

const fn numerical_refusal_reason_name(reason: CoupledNumericalRefusalReason) -> &'static str {
    match reason {
        CoupledNumericalRefusalReason::ReimpactLimitExceeded => "reimpact-limit-exceeded",
        CoupledNumericalRefusalReason::NonFiniteEnergyOrBaseState => {
            "non-finite-energy-or-base-state"
        }
    }
}

const fn censor_ordering_name(
    ordering: Result<CensorAwareDurationOrdering, CensorAwareRankingRefusal>,
) -> &'static str {
    match ordering {
        Ok(CensorAwareDurationOrdering::ProvenLeftShorter) => "ring-proven-shorter",
        Ok(CensorAwareDurationOrdering::EqualObserved) => "equal-observed",
        Ok(CensorAwareDurationOrdering::ProvenLeftLonger) => "ring-proven-longer",
        Ok(CensorAwareDurationOrdering::Indeterminate) => "indeterminate",
        Err(CensorAwareRankingRefusal::LeftNumericalRefusal) => "ring-numerical-refusal",
        Err(CensorAwareRankingRefusal::RightNumericalRefusal) => "solid-numerical-refusal",
        Err(CensorAwareRankingRefusal::InvalidLeftPhysicalTerminalTime) => {
            "ring-invalid-physical-terminal-time"
        }
        Err(CensorAwareRankingRefusal::InvalidRightPhysicalTerminalTime) => {
            "solid-invalid-physical-terminal-time"
        }
        Err(CensorAwareRankingRefusal::InvalidLeftCensorTime) => "ring-invalid-censor-time",
        Err(CensorAwareRankingRefusal::InvalidRightCensorTime) => "solid-invalid-censor-time",
        Err(CensorAwareRankingRefusal::DifferentPhysicalTerminalKinds) => {
            "different-physical-terminal-kinds"
        }
    }
}

fn event_time_delta(left: RunOutcome, right: RunOutcome) -> Option<(f64, f64)> {
    let (Some(left_time_s), Some(right_time_s)) = (
        left.observed_terminal_time_s(),
        right.observed_terminal_time_s(),
    ) else {
        return None;
    };
    let absolute_s = (left_time_s - right_time_s).abs();
    Some((
        absolute_s,
        absolute_s / right_time_s.abs().max(f64::MIN_POSITIVE),
    ))
}

fn event_delta_json(delta: Option<(f64, f64)>) -> String {
    match delta {
        Some((absolute_s, relative)) => format!(
            "{{\"available\":true,\"absolute_s\":{absolute_s:.17e},\"relative\":{relative:.17e}}}"
        ),
        None => "{\"available\":false}".to_owned(),
    }
}

fn calibration_readiness_record() -> Result<String, String> {
    let input = CalibrationReadinessInput {
        specimen: DeclaredEvidence::Missing,
        rig: DeclaredEvidence::Missing,
        instrument: DeclaredEvidence::Missing,
        raw_observations: DeclaredEvidence::Missing,
        observation_covariance: DeclaredEvidence::Missing,
        calibration_partition: DeclaredEvidence::Missing,
        blind_holdout: DeclaredEvidence::Missing,
    };
    let missing = missing_calibration_evidence(&input);
    let expected = [
        CalibrationEvidenceKind::Specimen,
        CalibrationEvidenceKind::Rig,
        CalibrationEvidenceKind::Instrument,
        CalibrationEvidenceKind::RawObservations,
        CalibrationEvidenceKind::ObservationCovariance,
        CalibrationEvidenceKind::CalibrationPartition,
        CalibrationEvidenceKind::BlindHoldout,
    ];
    if missing.as_slice() != expected {
        return Err(
            "calibration no-data record omitted or reordered a typed prerequisite".to_owned(),
        );
    }
    if !matches!(
        admit_calibration_readiness(input),
        Err(CalibrationReadinessError::MissingEvidence {
            kind: CalibrationEvidenceKind::Specimen,
        })
    ) {
        return Err("calibration readiness unexpectedly admitted missing evidence".to_owned());
    }
    let missing_evidence = missing
        .into_iter()
        .map(CalibrationEvidenceKind::slug)
        .collect::<Vec<_>>()
        .join(",");
    Ok(format!(
        "{{\"schema\":\"euler-disc-campaign-jsonl-v3\",\"scenario\":\"physical-calibration-readiness\",\"model\":\"structural-evidence-admission\",\"authority\":\"no-data\",\"units\":\"not-applicable\",\"terminal\":\"no-data\",\"missing_evidence\":\"{missing_evidence}\",\"synthetic_substitution\":false,\"target_ordering_fit\":false,\"no_claim\":\"no-calibration-or-physical-validation-without-retained-independent-evidence\"}}"
    ))
}

fn with_closed_context<R>(
    operation: impl FnOnce(&Cx<'_>) -> Result<R, String>,
) -> Result<R, String> {
    let gate = CancelGate::new_clock_free();
    let pool = ArenaPool::new(ArenaConfig::default());
    pool.scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: CAMPAIGN_SEED,
                kernel_id: 3,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        operation(&cx)
    })
}

fn run_campaign() -> Result<(), String> {
    // This is the physical squat specimen used by the geometry and exterior
    // screening rungs.  It is intentionally distinct from the ultra-thin
    // analytic oracle below.
    let radius_m = 0.038;
    let thickness_m = 0.006;
    let density_kg_m3 = 2_680.0;
    let sharp = specimen(
        radius_m,
        thickness_m,
        density_kg_m3,
        SquatDiscEdgeTreatment::Sharp,
    )?;
    let filleted = specimen(
        radius_m,
        thickness_m,
        density_kg_m3,
        SquatDiscEdgeTreatment::CircularFillet { radius: 0.001 },
    )?;
    let filleted_equal_mass = specimen(
        radius_m,
        thickness_m,
        density_kg_m3 * sharp.properties.mass / filleted.properties.mass,
        SquatDiscEdgeTreatment::CircularFillet { radius: 0.001 },
    )?;
    let contact = contact_snapshot(
        &sharp,
        &filleted,
        &filleted_equal_mass,
        radius_m,
        thickness_m,
    )?;
    let mut records = Vec::new();
    records.push(geometry_record(
        "geometry-sharp-squat-disc",
        "sharp",
        radius_m,
        thickness_m,
        0.0,
        density_kg_m3,
        sharp.properties,
    )?);
    records.push(geometry_record(
        "geometry-filleted-squat-disc",
        "circular-fillet",
        radius_m,
        thickness_m,
        0.001,
        density_kg_m3,
        filleted.properties,
    )?);

    // The conservative oracle's thin-disc relationship is not admitted for
    // the squat specimen.  Keep its distinct ultra-thin nominal geometry and
    // report it as an encoded analytic reference law.
    records.push(baseline_record(run_ideal_conservative_baseline(
        SquatDiscInput::nominal(),
    ))?);
    records.push(contact.record.clone());

    let base_input = nominal_base_input(contact.normal_reaction_n);
    let base = fs_euler_disc_e2e::run_reduced_base_response(&base_input)
        .map_err(|error| format!("base response: {error}"))?;
    let refinement = refine_reduced_base_response(&base_input)
        .map_err(|error| format!("base refinement: {error}"))?;
    let base_terminal = base
        .final_sample()
        .ok_or_else(|| "base response returned no samples".to_owned())?;
    records.push(format!(
        concat!(
            "{{\"schema\":\"{SCHEMA}\",\"scenario\":\"reduced-flexible-base\",",
            "\"model\":\"fs-solid-flat-single-patch-one-mode\",",
            "\"source\":\"{}\",\"authority\":\"simulation-only-local-operator-evidence\",",
            "\"units\":\"{SI_UNITS}\",\"inputs\":{{\"normal_force_n\":{:.17e}}},",
            "\"campaign_seed_manifest_ref\":\"campaign-complete\",\"budget\":{{\"timestep_s\":{:.17e},\"steps\":{},\"declared_cx_poll_quota\":{EXECUTION_POLL_QUOTA},\"declared_cx_cost_quota\":{EXECUTION_COST_QUOTA}}},",
            "\"terminal\":\"time-horizon-reached\",",
            "\"powers_w\":{{}},\"loads_n\":{{\"modal_force_n\":{:.17e}}},",
            "\"work_j\":{{\"damping\":{:.17e},\"external\":{:.17e}}},",
            "\"residual\":{{\"energy_j\":{:.17e},\"refinement_displacement_m\":{:.17e},",
            "\"refinement_elastic_energy_j\":{:.17e},\"refinement_normalized_energy\":{:.17e},",
            "\"normalized_energy_closure\":{:.17e},\"reduced_solve_scaled\":{:.17e},",
            "\"reduced_solve_scaled_limit\":{:.17e}}},\"no_claim\":\"{NO_PHYSICAL_VALIDATION};",
            "no-resolved-contact-or-curved-shell;{DECLARED_CX_BUDGET_NO_CLAIM}\"}}"
        ),
        base_input.plate.identity.source_id,
        base_input.load.normal_force_n,
        base_input.timestep_s,
        base_input.steps,
        base_terminal.modal_force_n,
        base_terminal.damping_work_j,
        base_terminal.external_work_j,
        base.energy_closure_residual_j,
        refinement.terminal_displacement_difference_m,
        refinement.terminal_elastic_energy_difference_j,
        refinement.terminal_normalized_energy_difference,
        base.normalized_energy_closure_residual,
        base.diagnostics.reduced_solve_scaled_residual,
        base.diagnostics.reduced_solve_scaled_residual_limit,
        SCHEMA = SCHEMA,
        SI_UNITS = SI_UNITS,
        EXECUTION_POLL_QUOTA = EXECUTION_POLL_QUOTA,
        EXECUTION_COST_QUOTA = EXECUTION_COST_QUOTA,
        NO_PHYSICAL_VALIDATION = NO_PHYSICAL_VALIDATION,
        DECLARED_CX_BUDGET_NO_CLAIM = DECLARED_CX_BUDGET_NO_CLAIM,
    ));

    records.push(reduced_decay_record(
        "contour-only-decay",
        filleted.properties.mass,
        radius_m,
        contact.normal_reaction_n,
        true,
        false,
    )?);
    records.push(reduced_decay_record(
        "boundary-layer-only-decay",
        filleted.properties.mass,
        radius_m,
        contact.normal_reaction_n,
        false,
        true,
    )?);
    records.push(reduced_decay_record(
        "combined-decay",
        filleted.properties.mass,
        radius_m,
        contact.normal_reaction_n,
        true,
        true,
    )?);
    records.push(flux_probe_record(
        radius_m,
        thickness_m,
        filleted.properties.mass,
        filleted.properties.principal_inertia.transverse,
        filleted.properties.principal_inertia.axial,
        contact.final_state.clone(),
    )?);

    let payload = records.join("\n");
    let digest_domain = "org.frankensim.euler-disc-campaign-jsonl.v1";
    let digest = fs_blake3::hash_domain(digest_domain, payload.as_bytes());
    let manifest = format!(
        "{{\"schema\":\"{SCHEMA}\",\"scenario\":\"campaign-complete\",\"model\":\"committed-production-rungs\",\"source\":\"campaign/committed-production-rungs\",\"authority\":\"integration-local\",\"units\":\"{SI_UNITS}\",\"campaign_seed_u64_dec\":\"{CAMPAIGN_SEED_DEC}\",\"budget\":{{\"record_count\":{},\"declared_cx_poll_quota\":{EXECUTION_POLL_QUOTA},\"declared_cx_cost_quota\":{EXECUTION_COST_QUOTA}}},\"terminal\":\"completed\",\"powers_w\":{{}},\"work_j\":{{}},\"residual\":{{\"record_count\":{}}},\"no_claim\":\"{NO_PHYSICAL_VALIDATION};{DECLARED_CX_BUDGET_NO_CLAIM}\",\"digest_domain\":\"{}\",\"digest_scope\":\"preceding-data-records-LF-joined-no-trailing-LF\",\"digest_blake3\":\"{}\"}}",
        records.len(),
        records.len(),
        digest_domain,
        digest.to_hex(),
    );
    for record in records {
        println!("{record}");
    }
    println!("{manifest}");
    Ok(())
}

struct Specimen {
    chart: AxisymmetricChart,
    properties: fs_rep_frep::AxisymmetricMassProperties,
    density_kg_m3: f64,
}

fn specimen(
    radius_m: f64,
    thickness_m: f64,
    density_kg_m3: f64,
    edge: SquatDiscEdgeTreatment,
) -> Result<Specimen, String> {
    let chart = AxisymmetricChart::squat_disc(radius_m, thickness_m, edge)
        .map_err(|error| format!("squat-disc chart: {error:?}"))?;
    let properties = with_context(|cx| {
        chart
            .mass_properties(density_kg_m3, cx)
            .map_err(|error| format!("mass properties: {error:?}"))
    })?;
    Ok(Specimen {
        chart,
        properties,
        density_kg_m3,
    })
}

fn with_context<R>(operation: impl FnOnce(&Cx<'_>) -> Result<R, String>) -> Result<R, String> {
    let gate = CancelGate::new_clock_free();
    let pool = ArenaPool::new(ArenaConfig::default());
    pool.scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: CAMPAIGN_SEED,
                kernel_id: 1,
                tile: 0,
                iteration: 0,
            },
            Budget::new()
                .with_poll_quota(EXECUTION_POLL_QUOTA)
                .with_cost_quota(EXECUTION_COST_QUOTA),
            ExecMode::Deterministic,
        );
        operation(&cx)
    })
}

fn geometry_record(
    scenario: &str,
    edge: &str,
    radius_m: f64,
    thickness_m: f64,
    edge_radius_m: f64,
    density_kg_m3: f64,
    properties: fs_rep_frep::AxisymmetricMassProperties,
) -> Result<String, String> {
    if !(properties.mass.is_finite() && properties.mass > 0.0) {
        return Err("geometry mass was not positive and finite".to_owned());
    }
    Ok(format!(
        concat!(
            "{{\"schema\":\"{SCHEMA}\",\"scenario\":\"{}\",",
            "\"model\":\"frep-axisymmetric-line-arc-mass-properties\",",
            "\"source\":\"frep/axisymmetric-line-arc\",",
            "\"authority\":\"non-authoritative-roundoff-telemetry\",\"units\":\"{SI_UNITS}\",\"campaign_seed_manifest_ref\":\"campaign-complete\",",
            "\"inputs\":{{\"edge\":\"{}\",\"radius_m\":{:.17e},\"thickness_m\":{:.17e},",
            "\"edge_radius_m\":{:.17e},\"density_kg_m3\":{:.17e}}},\"budget\":{{\"exec_mode\":\"deterministic\",",
            "\"seed_u64_dec\":\"{CAMPAIGN_SEED_DEC}\",\"declared_cx_poll_quota\":{EXECUTION_POLL_QUOTA},\"declared_cx_cost_quota\":{EXECUTION_COST_QUOTA}}},\"terminal\":\"completed\",",
            "\"powers_w\":{{}},\"work_j\":{{}},",
            "\"residual\":{{\"roundoff_volume_term_scale\":{:.17e}}},",
            "\"mass_kg\":{:.17e},\"inertia_transverse_kg_m2\":{:.17e},",
            "\"inertia_axial_kg_m2\":{:.17e},\"no_claim\":\"{NO_PHYSICAL_VALIDATION};",
            "density-is-declared-not-material-calibration-or-certified-mass-bound;{DECLARED_CX_BUDGET_NO_CLAIM}\"}}"
        ),
        scenario,
        edge,
        radius_m,
        thickness_m,
        edge_radius_m,
        density_kg_m3,
        properties.roundoff_diagnostics.volume_term_scale,
        properties.mass,
        properties.principal_inertia.transverse,
        properties.principal_inertia.axial,
        SCHEMA = SCHEMA,
        SI_UNITS = SI_UNITS,
        CAMPAIGN_SEED_DEC = CAMPAIGN_SEED_DEC,
        EXECUTION_POLL_QUOTA = EXECUTION_POLL_QUOTA,
        EXECUTION_COST_QUOTA = EXECUTION_COST_QUOTA,
        NO_PHYSICAL_VALIDATION = NO_PHYSICAL_VALIDATION,
        DECLARED_CX_BUDGET_NO_CLAIM = DECLARED_CX_BUDGET_NO_CLAIM,
    ))
}

fn baseline_record(output: BaselineRunOutput) -> Result<String, String> {
    let BaselineRunOutput::Completed(trajectory) = output else {
        return Err("conservative baseline admission refused".to_owned());
    };
    let sample = trajectory
        .samples
        .last()
        .ok_or_else(|| "baseline returned no samples".to_owned())?;
    Ok(format!(
        concat!(
            "{{\"schema\":\"{SCHEMA}\",\"scenario\":\"conservative-steady-oracle\",",
            "\"model\":\"{}\",\"source\":\"{}\",",
            "\"authority\":\"caller-declared-thin-disc-oracle\",\"units\":\"{SI_UNITS}\",\"campaign_seed_manifest_ref\":\"campaign-complete\",",
            "\"inputs\":{{\"radius_m\":{:.17e},\"thickness_m\":{:.17e},\"mass_kg\":{:.17e},\"inclination_rad\":{:.17e}}},",
            "\"budget\":{{\"timestep_s\":{:.17e},\"steps\":{}}},\"terminal\":\"{}\",",
            "\"powers_w\":{{\"ideal_support\":0.0}},\"work_j\":{{\"ideal_support\":{:.17e}}},",
            "\"residual\":{{\"energy_j\":{:.17e},\"precession_s_inv2\":{:.17e},",
            "\"small_angle_energy_j\":{:.17e}}},\"no_claim\":\"{NO_PHYSICAL_VALIDATION};",
            "encoded-analytic-reference-law;not-squat-specimen-mechanics\"}}"
        ),
        trajectory.model_id,
        trajectory.input.interface_source_id,
        trajectory.input.radius_m,
        trajectory.input.thickness_m,
        trajectory.input.mass_kg,
        trajectory.input.inclination_from_vertical_rad,
        trajectory.input.step_seconds,
        trajectory.input.steps,
        baseline_terminal_label(trajectory.terminal),
        sample.energy.ideal_support_work_j,
        sample.energy.residual_from_initial_j,
        trajectory.equilibrium.precession_balance_residual_s_inv2,
        trajectory.equilibrium.small_angle_energy_residual_j,
        SCHEMA = SCHEMA,
        SI_UNITS = SI_UNITS,
        NO_PHYSICAL_VALIDATION = NO_PHYSICAL_VALIDATION,
    ))
}

#[derive(Clone)]
struct ContactSnapshot {
    normal_reaction_n: f64,
    final_state: fs_mbd::RigidBodyState,
    record: String,
}

struct ContactCase {
    label: &'static str,
    density_kg_m3: f64,
    properties: fs_rep_frep::AxisymmetricMassProperties,
    initializer: fs_euler_disc_e2e::ProfileRollingInitializer,
    input: ProfileContactDynamicsInput,
    run: fs_euler_disc_e2e::ContactDynamicsRun,
    refinement: fs_euler_disc_e2e::TimestepRefinement,
    mean_reaction_n: f64,
    final_reaction_n: f64,
    peak_reaction_n: f64,
    max_required_static_mu: f64,
    summed_mechanical_balance_residual_j: f64,
    sum_abs_mechanical_balance_residual_j: f64,
    max_abs_mechanical_balance_residual_j: f64,
    mechanical_energy_scale_j: f64,
}

fn contact_snapshot(
    sharp: &Specimen,
    filleted: &Specimen,
    filleted_equal_mass: &Specimen,
    radius_m: f64,
    thickness_m: f64,
) -> Result<ContactSnapshot, String> {
    let sharp_case = run_contact_case("sharp-fixed-density", sharp, radius_m, thickness_m)?;
    let filleted_case = run_contact_case("fillet-fixed-density", filleted, radius_m, thickness_m)?;
    let equal_mass_case = run_contact_case(
        "fillet-equal-mass",
        filleted_equal_mass,
        radius_m,
        thickness_m,
    )?;
    let normal_reaction_n = filleted_case.mean_reaction_n;
    let record = format!(
        concat!(
            "{{\"schema\":\"{SCHEMA}\",\"scenario\":\"dynamic-unilateral-contact\",",
            "\"model\":\"axisymmetric-profile-unilateral-sticking-contact\",",
            "\"source\":\"frep/matched-profile-contact-comparison\",",
            "\"authority\":\"numerical-reference-only\",\"units\":\"{SI_UNITS}\",\"campaign_seed_manifest_ref\":\"campaign-complete\",",
            "\"inputs\":{{\"primary_case\":\"fillet-fixed-density\",\"radius_m\":{:.17e},\"thickness_m\":{:.17e},\"edge_radius_m\":1.0e-3,\"rolling_inclination_rad\":5.0e-2,\"static_friction_coefficient\":100.0,\"interface_system_id\":\"campaign/squat-disc->plane\",\"interface_history_id\":\"campaign/shared-matched-profile-contact-history-v1\",\"dry_interface_authority\":\"caller-declared-dry-interface\"}},",
            "\"budget\":{{\"timestep_s\":1.0e-6,\"maximum_steps\":8,\"declared_cx_poll_quota\":{EXECUTION_POLL_QUOTA},\"declared_cx_cost_quota\":{EXECUTION_COST_QUOTA}}},\"terminal\":\"horizon-reached\",",
            "\"cases\":{{\"sharp_fixed_density\":{},\"fillet_fixed_density\":{},\"fillet_equal_mass\":{}}},",
            "\"deltas\":{{\"fillet_fixed_density_minus_sharp_fixed_density\":{},\"fillet_equal_mass_minus_sharp_fixed_density\":{}}},",
            "\"loads_n\":{{\"normal_reaction_mean_n\":{:.17e}}},",
            "\"powers_w\":{{\"gravity\":{:.17e},\"contact_impulse\":{:.17e}}},",
            "\"work_j\":{{}},",
            "\"residual\":{{\"energy_j\":{:.17e},\"sum_abs_energy_j\":{:.17e},\"max_abs_energy_j\":{:.17e},\"energy_scale_j\":{:.17e},\"energy_residual_limit_j\":{:.17e},\"sum_abs_energy_limit_j\":{:.17e},\"max_abs_energy_limit_j\":{:.17e},\"refinement_position_m\":{:.17e},\"refinement_position_limit_m\":{:.17e},\"refinement_energy_j\":{:.17e},\"refinement_energy_limit_j\":{:.17e}}},",
            "\"no_claim\":\"{NO_PHYSICAL_VALIDATION};one-way-not-closed-coupling;mu-100-branch-isolation-control;matched-geometry-comparison-not-outcome-ranking-or-physical-validation;{DECLARED_CX_BUDGET_NO_CLAIM}\"}}"
        ),
        radius_m,
        thickness_m,
        contact_case_json(&sharp_case)?,
        contact_case_json(&filleted_case)?,
        contact_case_json(&equal_mass_case)?,
        contact_case_delta_json(&filleted_case, &sharp_case),
        contact_case_delta_json(&equal_mass_case, &sharp_case),
        normal_reaction_n,
        average_gravity_power(&filleted_case),
        average_contact_power(&filleted_case),
        filleted_case.summed_mechanical_balance_residual_j,
        filleted_case.sum_abs_mechanical_balance_residual_j,
        filleted_case.max_abs_mechanical_balance_residual_j,
        filleted_case.mechanical_energy_scale_j,
        filleted_case.mechanical_energy_scale_j * 1.0e-4,
        filleted_case.mechanical_energy_scale_j * 1.0e-3,
        filleted_case.mechanical_energy_scale_j * 1.0e-4,
        filleted_case.refinement.final_position_difference_m,
        radius_m * 1.0e-1,
        filleted_case
            .refinement
            .final_mechanical_energy_difference_j,
        filleted_case.mechanical_energy_scale_j * 1.0e-1,
        SCHEMA = SCHEMA,
        SI_UNITS = SI_UNITS,
        EXECUTION_POLL_QUOTA = EXECUTION_POLL_QUOTA,
        EXECUTION_COST_QUOTA = EXECUTION_COST_QUOTA,
        NO_PHYSICAL_VALIDATION = NO_PHYSICAL_VALIDATION,
        DECLARED_CX_BUDGET_NO_CLAIM = DECLARED_CX_BUDGET_NO_CLAIM,
    );
    Ok(ContactSnapshot {
        normal_reaction_n,
        final_state: filleted_case.run.final_state,
        record,
    })
}

fn run_contact_case(
    label: &'static str,
    specimen: &Specimen,
    radius_m: f64,
    thickness_m: f64,
) -> Result<ContactCase, String> {
    let initializer = with_context(|cx| {
        small_angle_rolling_profile_initializer(
            &specimen.chart,
            specimen.density_kg_m3,
            0.05,
            9.806_65,
            cx,
        )
        .map_err(|error| format!("{label} rolling initializer: {error}"))
    })?;
    let controls = ContactDynamicsInput {
        geometry: ContactDiscGeometry {
            radius_m,
            thickness_m,
            mass_kg: specimen.properties.mass,
        },
        initial_state: initializer.state,
        gravity_m_per_s2: 9.806_65,
        // Deliberately high branch-isolation control, not a measured friction coefficient.
        static_friction_coefficient: 100.0,
        interface: InterfaceSystemRef::new(
            "campaign/squat-disc->plane",
            "campaign/shared-matched-profile-contact-history-v1",
            "campaign/caller-declared-dry-interface-v1",
            InputAuthority::CallerDeclared,
            InterfaceMedium::Dry,
        )
        .map_err(|error| format!("{label} contact interface: {error}"))?,
        timestep_s: 1.0e-6,
        maximum_steps: 8,
        contact_tolerance_m: 1.0e-9,
        maximum_initial_penetration_m: 1.0e-10,
        release_speed_tolerance_m_per_s: 1.0e-8,
    };
    let input = ProfileContactDynamicsInput {
        chart: specimen.chart.clone(),
        density_kg_per_m3: specimen.density_kg_m3,
        controls,
    };
    let timestep_s = input.controls.timestep_s;
    let run = with_context(|cx| {
        run_profile_contact_dynamics(&input, cx)
            .map_err(|error| format!("{label} profile dynamic contact: {error}"))
    })?;
    if !matches!(run.termination, ContactTermination::HorizonReached) {
        return Err(format!(
            "{label} contact did not reach its required horizon"
        ));
    }
    let refinement = with_context(|cx| {
        refine_profile_timestep_by_two(&input, cx)
            .map_err(|error| format!("{label} profile refinement: {error}"))
    })?;
    let mut impulse_sum = 0.0;
    let mut peak_reaction_n = 0.0_f64;
    let mut max_required_static_mu = 0.0_f64;
    let mut summed_mechanical_balance_residual_j = 0.0;
    let mut sum_abs_mechanical_balance_residual_j = 0.0;
    let mut max_abs_mechanical_balance_residual_j = 0.0_f64;
    for step in &run.steps {
        if !(step.normal_impulse_ns.is_finite() && step.normal_impulse_ns > 0.0) {
            return Err(format!("{label} has zero or invalid normal impulse"));
        }
        impulse_sum += step.normal_impulse_ns;
        peak_reaction_n = peak_reaction_n.max(step.normal_impulse_ns / timestep_s);
        max_required_static_mu = max_required_static_mu
            .max(step.stick.required_tangential_impulse_ns / step.stick.normal_impulse_ns);
        let residual_j = step.energy.mechanical_balance_residual_j;
        summed_mechanical_balance_residual_j += residual_j;
        sum_abs_mechanical_balance_residual_j += residual_j.abs();
        max_abs_mechanical_balance_residual_j =
            max_abs_mechanical_balance_residual_j.max(residual_j.abs());
    }
    let count = run.steps.len();
    if count == 0 {
        return Err(format!("{label} contact retained no steps"));
    }
    let final_reaction_n = run
        .steps
        .last()
        .ok_or_else(|| format!("{label} contact retained no terminal step"))?
        .normal_impulse_ns
        / timestep_s;
    Ok(ContactCase {
        label,
        density_kg_m3: specimen.density_kg_m3,
        properties: specimen.properties,
        initializer,
        input,
        run,
        refinement,
        mean_reaction_n: impulse_sum / (count as f64 * timestep_s),
        final_reaction_n,
        peak_reaction_n,
        max_required_static_mu,
        summed_mechanical_balance_residual_j,
        sum_abs_mechanical_balance_residual_j,
        max_abs_mechanical_balance_residual_j,
        mechanical_energy_scale_j: (specimen.properties.mass * 9.806_65 * radius_m).max(1.0e-12),
    })
}

fn contact_case_json(case: &ContactCase) -> Result<String, String> {
    validate_contact_case(case)?;
    let support = case.initializer.contact.contact.radius_world_m;
    let velocity = case.initializer.initial_contact_velocity_world_m_per_s;
    Ok(format!(
        "{{\"label\":\"{}\",\"density_kg_m3\":{:.17e},\"mass_kg\":{:.17e},\"inertia_transverse_kg_m2\":{:.17e},\"inertia_axial_kg_m2\":{:.17e},\"support_vector_world_m\":{{\"x\":{:.17e},\"y\":{:.17e},\"z\":{:.17e}}},\"center_height_m\":{:.17e},\"initial_material_contact_speed_m_per_s\":{:.17e},\"mean_reaction_n\":{:.17e},\"final_reaction_n\":{:.17e},\"peak_reaction_n\":{:.17e},\"max_required_static_mu\":{:.17e},\"summed_mechanical_balance_residual_j\":{:.17e},\"sum_abs_mechanical_balance_residual_j\":{:.17e},\"max_abs_mechanical_balance_residual_j\":{:.17e},\"mechanical_energy_scale_j\":{:.17e},\"energy_residual_limit_j\":{:.17e},\"sum_abs_energy_limit_j\":{:.17e},\"max_abs_energy_limit_j\":{:.17e},\"refinement_position_m\":{:.17e},\"refinement_energy_j\":{:.17e},\"quarter_step_refinement\":{},\"terminal\":\"horizon-reached\"}}",
        case.label,
        case.density_kg_m3,
        case.properties.mass,
        case.properties.principal_inertia.transverse,
        case.properties.principal_inertia.axial,
        support.x,
        support.y,
        support.z,
        case.initializer.state.pose().position_world().z,
        (velocity.x.mul_add(
            velocity.x,
            velocity.y.mul_add(velocity.y, velocity.z * velocity.z)
        ))
        .sqrt(),
        case.mean_reaction_n,
        case.final_reaction_n,
        case.peak_reaction_n,
        case.max_required_static_mu,
        case.summed_mechanical_balance_residual_j,
        case.sum_abs_mechanical_balance_residual_j,
        case.max_abs_mechanical_balance_residual_j,
        case.mechanical_energy_scale_j,
        case.mechanical_energy_scale_j * 1.0e-4,
        case.mechanical_energy_scale_j * 1.0e-3,
        case.mechanical_energy_scale_j * 1.0e-4,
        case.refinement.final_position_difference_m,
        case.refinement.final_mechanical_energy_difference_j,
        quarter_step_refinement_json(case)?,
    ))
}

fn validate_contact_case(case: &ContactCase) -> Result<(), String> {
    for (rung, run) in [
        ("primary", &case.run),
        ("coarse", &case.refinement.coarse),
        ("half-step", &case.refinement.fine),
        ("quarter-step", &case.refinement.reference),
    ] {
        if !matches!(run.termination, ContactTermination::HorizonReached) {
            return Err(format!(
                "{} {rung} contact did not reach its required horizon",
                case.label
            ));
        }
    }

    let support = case.initializer.contact.contact.radius_world_m;
    let velocity = case.initializer.initial_contact_velocity_world_m_per_s;
    let refinement = &case.refinement;
    let finite_values = [
        ("density_kg_m3", case.density_kg_m3),
        ("mass_kg", case.properties.mass),
        (
            "inertia_transverse_kg_m2",
            case.properties.principal_inertia.transverse,
        ),
        (
            "inertia_axial_kg_m2",
            case.properties.principal_inertia.axial,
        ),
        ("support_x_m", support.x),
        ("support_y_m", support.y),
        ("support_z_m", support.z),
        (
            "center_height_m",
            case.initializer.state.pose().position_world().z,
        ),
        ("initial_contact_velocity_x_m_per_s", velocity.x),
        ("initial_contact_velocity_y_m_per_s", velocity.y),
        ("initial_contact_velocity_z_m_per_s", velocity.z),
        ("mean_reaction_n", case.mean_reaction_n),
        ("final_reaction_n", case.final_reaction_n),
        ("peak_reaction_n", case.peak_reaction_n),
        ("max_required_static_mu", case.max_required_static_mu),
        (
            "summed_mechanical_balance_residual_j",
            case.summed_mechanical_balance_residual_j,
        ),
        (
            "sum_abs_mechanical_balance_residual_j",
            case.sum_abs_mechanical_balance_residual_j,
        ),
        (
            "max_abs_mechanical_balance_residual_j",
            case.max_abs_mechanical_balance_residual_j,
        ),
        ("mechanical_energy_scale_j", case.mechanical_energy_scale_j),
        (
            "refinement_position_m",
            refinement.final_position_difference_m,
        ),
        (
            "refinement_energy_j",
            refinement.final_mechanical_energy_difference_j,
        ),
        (
            "coarse_reference_position_m",
            refinement.coarse_reference_position_difference_m,
        ),
        (
            "fine_reference_position_m",
            refinement.fine_reference_position_difference_m,
        ),
        (
            "coarse_reference_linear_momentum",
            refinement.coarse_reference_linear_momentum_difference_kg_m_per_s,
        ),
        (
            "fine_reference_linear_momentum",
            refinement.fine_reference_linear_momentum_difference_kg_m_per_s,
        ),
        (
            "coarse_reference_angular_momentum",
            refinement.coarse_reference_angular_momentum_difference_kg_m2_per_s,
        ),
        (
            "fine_reference_angular_momentum",
            refinement.fine_reference_angular_momentum_difference_kg_m2_per_s,
        ),
        (
            "coarse_reference_orientation_rad",
            refinement.coarse_reference_orientation_angle_rad,
        ),
        (
            "fine_reference_orientation_rad",
            refinement.fine_reference_orientation_angle_rad,
        ),
    ];
    for (field, value) in finite_values {
        if !value.is_finite() {
            return Err(format!("{} emitted non-finite {field}", case.label));
        }
    }

    let energy_limit_j = case.mechanical_energy_scale_j * 1.0e-4;
    let sum_abs_limit_j = case.mechanical_energy_scale_j * 1.0e-3;
    let position_limit_m = case.input.controls.geometry.radius_m * 1.0e-1;
    let refinement_energy_limit_j = case.mechanical_energy_scale_j * 1.0e-1;
    if !(case.mechanical_energy_scale_j > 0.0
        && case.mean_reaction_n > 0.0
        && case.final_reaction_n > 0.0
        && case.peak_reaction_n > 0.0
        && case.max_required_static_mu >= 0.0)
    {
        return Err(format!(
            "{} contact produced a non-positive scale or reaction",
            case.label
        ));
    }
    if case.summed_mechanical_balance_residual_j.abs() > energy_limit_j
        || case.sum_abs_mechanical_balance_residual_j > sum_abs_limit_j
        || case.max_abs_mechanical_balance_residual_j > energy_limit_j
    {
        return Err(format!(
            "{} contact energy residual exceeded its declared bound",
            case.label
        ));
    }
    if refinement.final_position_difference_m > position_limit_m
        || refinement.final_mechanical_energy_difference_j.abs() > refinement_energy_limit_j
    {
        return Err(format!(
            "{} contact timestep refinement exceeded its declared bound",
            case.label
        ));
    }
    if !(refinement.position_refinement_improved
        && refinement.linear_momentum_refinement_improved
        && refinement.angular_momentum_refinement_improved
        && refinement.orientation_refinement_improved)
    {
        return Err(format!(
            "{} half-step contact endpoint did not improve against the quarter-step reference",
            case.label
        ));
    }

    let reference_energy_j = terminal_mechanical_energy_j(&refinement.reference, "quarter-step")?;
    let coarse_energy_error_j =
        (terminal_mechanical_energy_j(&refinement.coarse, "coarse")? - reference_energy_j).abs();
    let fine_energy_error_j =
        (terminal_mechanical_energy_j(&refinement.fine, "half-step")? - reference_energy_j).abs();
    if !(reference_energy_j.is_finite()
        && coarse_energy_error_j.is_finite()
        && fine_energy_error_j.is_finite()
        && fine_energy_error_j <= coarse_energy_error_j)
    {
        return Err(format!(
            "{} half-step contact energy did not improve against the quarter-step reference",
            case.label
        ));
    }
    Ok(())
}

fn terminal_mechanical_energy_j(
    run: &fs_euler_disc_e2e::ContactDynamicsRun,
    rung: &str,
) -> Result<f64, String> {
    run.steps
        .last()
        .map(|step| step.energy.mechanical_energy_after_j)
        .ok_or_else(|| format!("{rung} refinement retained no terminal energy step"))
}

fn quarter_step_refinement_json(case: &ContactCase) -> Result<String, String> {
    let refinement = &case.refinement;
    let reference_energy_j = terminal_mechanical_energy_j(&refinement.reference, "quarter-step")?;
    let coarse_energy_error_j =
        (terminal_mechanical_energy_j(&refinement.coarse, "coarse")? - reference_energy_j).abs();
    let fine_energy_error_j =
        (terminal_mechanical_energy_j(&refinement.fine, "half-step")? - reference_energy_j).abs();
    Ok(format!(
        concat!(
            "{{\"coarse_position_error_m\":{:.17e},\"fine_position_error_m\":{:.17e},",
            "\"coarse_linear_momentum_error_kg_m_per_s\":{:.17e},\"fine_linear_momentum_error_kg_m_per_s\":{:.17e},",
            "\"coarse_angular_momentum_error_kg_m2_per_s\":{:.17e},\"fine_angular_momentum_error_kg_m2_per_s\":{:.17e},",
            "\"coarse_orientation_error_rad\":{:.17e},\"fine_orientation_error_rad\":{:.17e},",
            "\"coarse_energy_error_j\":{:.17e},\"fine_energy_error_j\":{:.17e},",
            "\"position_refinement_improved\":{},\"linear_momentum_refinement_improved\":{},",
            "\"angular_momentum_refinement_improved\":{},\"orientation_refinement_improved\":{},\"energy_refinement_improved\":{}}}"
        ),
        refinement.coarse_reference_position_difference_m,
        refinement.fine_reference_position_difference_m,
        refinement.coarse_reference_linear_momentum_difference_kg_m_per_s,
        refinement.fine_reference_linear_momentum_difference_kg_m_per_s,
        refinement.coarse_reference_angular_momentum_difference_kg_m2_per_s,
        refinement.fine_reference_angular_momentum_difference_kg_m2_per_s,
        refinement.coarse_reference_orientation_angle_rad,
        refinement.fine_reference_orientation_angle_rad,
        coarse_energy_error_j,
        fine_energy_error_j,
        refinement.position_refinement_improved,
        refinement.linear_momentum_refinement_improved,
        refinement.angular_momentum_refinement_improved,
        refinement.orientation_refinement_improved,
        fine_energy_error_j <= coarse_energy_error_j,
    ))
}

fn contact_case_delta_json(lhs: &ContactCase, rhs: &ContactCase) -> String {
    format!(
        "{{\"mean_reaction_n\":{:.17e},\"final_reaction_n\":{:.17e},\"peak_reaction_n\":{:.17e},\"max_required_static_mu\":{:.17e},\"summed_mechanical_balance_residual_j\":{:.17e},\"sum_abs_mechanical_balance_residual_j\":{:.17e},\"max_abs_mechanical_balance_residual_j\":{:.17e},\"mechanical_energy_scale_j\":{:.17e},\"refinement_position_m\":{:.17e},\"refinement_energy_j\":{:.17e}}}",
        lhs.mean_reaction_n - rhs.mean_reaction_n,
        lhs.final_reaction_n - rhs.final_reaction_n,
        lhs.peak_reaction_n - rhs.peak_reaction_n,
        lhs.max_required_static_mu - rhs.max_required_static_mu,
        lhs.summed_mechanical_balance_residual_j - rhs.summed_mechanical_balance_residual_j,
        lhs.sum_abs_mechanical_balance_residual_j - rhs.sum_abs_mechanical_balance_residual_j,
        lhs.max_abs_mechanical_balance_residual_j - rhs.max_abs_mechanical_balance_residual_j,
        lhs.mechanical_energy_scale_j - rhs.mechanical_energy_scale_j,
        lhs.refinement.final_position_difference_m - rhs.refinement.final_position_difference_m,
        lhs.refinement.final_mechanical_energy_difference_j
            - rhs.refinement.final_mechanical_energy_difference_j,
    )
}

fn average_gravity_power(case: &ContactCase) -> f64 {
    case.run
        .steps
        .iter()
        .map(|step| step.energy.gravity_work_j)
        .sum::<f64>()
        / (case.run.steps.len() as f64 * case.input.controls.timestep_s)
}

fn average_contact_power(case: &ContactCase) -> f64 {
    case.run
        .steps
        .iter()
        .map(|step| step.energy.contact_impulse_work_estimate_j)
        .sum::<f64>()
        / (case.run.steps.len() as f64 * case.input.controls.timestep_s)
}

fn baseline_terminal_label(terminal: fs_euler_disc_e2e::BaselineTerminal) -> &'static str {
    match terminal {
        fs_euler_disc_e2e::BaselineTerminal::TimeHorizonReached { .. } => "time-horizon-reached",
        fs_euler_disc_e2e::BaselineTerminal::StaticFrictionCapacityExceeded { .. } => {
            "static-friction-capacity-exceeded"
        }
    }
}

fn nominal_base_input(normal_force_n: f64) -> BaseResponseInput {
    let support = ShellSupport {
        node_indices: [0, 1, 2],
        normal: [0.0, 0.0, 1.0],
    };
    BaseResponseInput {
        plate: ShellPlate {
            nodes: vec![
                ShellNode {
                    position_m: [-0.10, -0.08, 0.0],
                },
                ShellNode {
                    position_m: [0.10, -0.08, 0.0],
                },
                ShellNode {
                    position_m: [0.0, 0.12, 0.0],
                },
                ShellNode {
                    position_m: [0.0, 0.0, 0.0],
                },
            ],
            triangles: vec![[0, 1, 3], [1, 2, 3], [2, 0, 3]],
            thickness_m: 0.004,
            material: ShellMaterial {
                youngs_modulus_pa: 70.0e9,
                poisson_ratio: 0.33,
                density_kg_m3: 2_700.0,
            },
            identity: ShellIdentity {
                model_id: "campaign/reduced-base-response-v1".to_owned(),
                source_id: "campaign/flat-tripod-plate-v1".to_owned(),
                state_id: "initial".to_owned(),
            },
            support: Some(support),
            damping: DampingModel::Rayleigh {
                mass_proportional_per_s: 0.2,
                stiffness_proportional_s: 1.0e-6,
            },
            budget: AssemblyBudget::default(),
        },
        level_support: LevelSupportInput {
            support,
            level_normal: [0.0, 0.0, 1.0],
            maximum_tilt_rad: 1.0e-6,
        },
        geometry_scope: BaseGeometryScope::FlatSinglePatch,
        contact_scope: ContactLoadScope::NodalNormalLoad,
        load: MovingContactLoad {
            start_node: 3,
            end_node: 1,
            normal_force_n,
        },
        initial_modal_displacement_m: 0.0,
        initial_modal_velocity_m_per_s: 0.0,
        timestep_s: 1.0e-6,
        steps: 100,
    }
}

fn reduced_decay_record(
    scenario: &str,
    mass_kg: f64,
    radius_m: f64,
    normal_force_n: f64,
    include_dry: bool,
    include_boundary_layer: bool,
) -> Result<String, String> {
    let mut input = fs_euler_disc_e2e::ReducedDecayInput::nominal_reference()
        .map_err(|error| format!("reduced-decay input: {error}"))?;
    input.mass_kg = mass_kg;
    input.radius_m = radius_m;
    // Keep the same one-second declared horizon as the nominal 1e-5/100k
    // reference while bounding retained trajectory samples for a real campaign.
    input.timestep_s = 1.0e-4;
    input.maximum_steps = 10_000;
    input.dry_contour = include_dry.then(|| input.dry_contour.clone()).flatten();
    input.bildsten_boundary_layer = include_boundary_layer
        .then(|| input.bildsten_boundary_layer.clone())
        .flatten();
    if let Some(channel) = input.dry_contour.as_mut() {
        channel.normal_force_n = normal_force_n;
        channel.contour_force_n = 0.0;
    }
    let run = run_reduced_decay(&input).map_err(|error| format!("reduced decay: {error}"))?;
    let refinement = refinement_evidence(&input)
        .map_err(|error| format!("reduced decay refinement: {error}"))?;
    let final_sample = run
        .final_sample()
        .map_err(|error| format!("reduced decay final sample: {error}"))?;
    let energy_scale_j =
        (input.mass_kg * input.gravity_m_per_s2 * input.radius_m * input.initial_theta_rad)
            .max(1.0e-12);
    let closure_limit_j = energy_scale_j * 1.0e-4;
    let refinement_work_limit_j = energy_scale_j * 1.0e-1;
    let refinement_time_limit_s = final_sample.time_s.max(input.timestep_s) * 1.0e-1;
    if !matches!(
        run.terminal,
        fs_euler_disc_e2e::ReducedDecayTerminal::ValidityCutoff
    ) {
        return Err(format!(
            "{scenario} exhausted its step budget before the validity cutoff"
        ));
    }
    for (field, value) in [
        ("final_time_s", final_sample.time_s),
        ("final_theta_rad", final_sample.theta_rad),
        ("final_omega_rad_s", final_sample.omega_rad_s),
        ("final_energy_j", final_sample.energy_j),
        ("dry_power_w", final_sample.powers.dry_contour_w),
        (
            "boundary_layer_power_w",
            final_sample.powers.bildsten_boundary_layer_w,
        ),
        ("dry_work_j", final_sample.work.dry_contour_j),
        (
            "boundary_layer_work_j",
            final_sample.work.bildsten_boundary_layer_j,
        ),
        ("energy_closure_residual_j", run.energy_closure_residual_j),
        (
            "refinement_terminal_time_difference_s",
            refinement.terminal_time_difference_s,
        ),
        (
            "refinement_total_work_difference_j",
            refinement.total_work_difference_j,
        ),
    ] {
        if !value.is_finite() {
            return Err(format!("{scenario} emitted non-finite {field}"));
        }
    }
    if run.energy_closure_residual_j.abs() > closure_limit_j
        || refinement.terminal_time_difference_s > refinement_time_limit_s
        || refinement.total_work_difference_j.abs() > refinement_work_limit_j
    {
        return Err(format!(
            "{scenario} exceeded its declared closure or refinement bound"
        ));
    }
    let dry_source = run.provenance.dry_source_id.as_deref().unwrap_or("none");
    let bildsten_source = run
        .provenance
        .bildsten_source_id
        .as_deref()
        .unwrap_or("none");
    let primary_source = if dry_source != "none" {
        dry_source
    } else {
        bildsten_source
    };
    let terminal = "validity-cutoff";
    let crossover = crossover_receipt(&input, &run)?;
    Ok(format!(
        concat!(
            "{{\"schema\":\"{SCHEMA}\",\"scenario\":\"{}\",\"model\":\"{}\",",
            "\"source\":\"{}\",\"authority\":\"{}\",\"units\":\"{SI_UNITS}\",\"campaign_seed_manifest_ref\":\"campaign-complete\",",
            "\"sources\":{{\"small_angle_oracle\":\"{}\",\"dry_interface\":\"{}\",\"dry\":\"{}\",\"dry_contour_scaling\":\"caller-declared/normal-scaled-contour-force-v1\",\"bildsten\":\"{}\"}},",
            "\"inputs\":{{\"mass_kg\":{:.17e},\"radius_m\":{:.17e},\"initial_theta_rad\":{:.17e},\"validity_cutoff_theta_rad\":{:.17e},\"gravity_m_per_s2\":{:.17e},\"dry_normal_force_n\":{:.17e},\"contour_force_n\":{:.17e},\"contour_force_per_normal_force\":{:.17e},\"rho_kg_m3\":{:.17e},\"mu_pa_s\":{:.17e},\"prefactor\":{:.17e}}},",
            "\"budget\":{{\"timestep_s\":{:.17e},\"maximum_steps\":{}}},\"terminal\":\"{}\",",
            "\"powers_w\":{{\"dry_contour\":{:.17e},\"bildsten_boundary_layer\":{:.17e}}},",
            "\"work_j\":{{\"dry_contour\":{:.17e},\"bildsten_boundary_layer\":{:.17e}}},",
            "\"final\":{{\"time_s\":{:.17e},\"theta_rad\":{:.17e},\"omega_rad_s\":{:.17e},\"energy_j\":{:.17e}}},",
            "\"crossover\":{},",
            "\"residual\":{{\"energy_j\":{:.17e},\"energy_scale_j\":{:.17e},\"energy_residual_limit_j\":{:.17e},\"refinement_time_s\":{:.17e},\"refinement_time_limit_s\":{:.17e},",
            "\"refinement_work_j\":{:.17e},\"refinement_work_limit_j\":{:.17e}}},\"no_claim\":\"{};{};{};",
            "numerical-self-consistency-only\"}}"
        ),
        scenario,
        run.provenance.model_id,
        primary_source,
        run.provenance.model_authority,
        run.provenance.small_angle_oracle_source_id,
        run.provenance
            .dry_interface_system_id
            .as_deref()
            .unwrap_or("none"),
        dry_source,
        bildsten_source,
        input.mass_kg,
        input.radius_m,
        input.initial_theta_rad,
        input.validity_cutoff_theta_rad,
        input.gravity_m_per_s2,
        input
            .dry_contour
            .as_ref()
            .map_or(0.0, |channel| channel.normal_force_n),
        input
            .dry_contour
            .as_ref()
            .map_or(0.0, |channel| channel.contour_force_n),
        0.0,
        input
            .bildsten_boundary_layer
            .as_ref()
            .map_or(0.0, |channel| channel.density_kg_per_m3),
        input
            .bildsten_boundary_layer
            .as_ref()
            .map_or(0.0, |channel| channel.dynamic_viscosity_pa_s),
        input
            .bildsten_boundary_layer
            .as_ref()
            .map_or(0.0, |channel| channel.dimensionless_prefactor),
        input.timestep_s,
        input.maximum_steps,
        terminal,
        final_sample.powers.dry_contour_w,
        final_sample.powers.bildsten_boundary_layer_w,
        final_sample.work.dry_contour_j,
        final_sample.work.bildsten_boundary_layer_j,
        final_sample.time_s,
        final_sample.theta_rad,
        final_sample.omega_rad_s,
        final_sample.energy_j,
        crossover,
        run.energy_closure_residual_j,
        energy_scale_j,
        closure_limit_j,
        refinement.terminal_time_difference_s,
        refinement_time_limit_s,
        refinement.total_work_difference_j,
        refinement_work_limit_j,
        NO_PHYSICAL_VALIDATION,
        run.provenance.physical_validation,
        run.provenance.cancellation_capability,
        SCHEMA = SCHEMA,
        SI_UNITS = SI_UNITS,
    ))
}

fn crossover_receipt(
    input: &fs_euler_disc_e2e::ReducedDecayInput,
    run: &fs_euler_disc_e2e::ReducedDecayRun,
) -> Result<String, String> {
    match channel_crossover_diagnostic(input)
        .map_err(|error| format!("reduced decay crossover: {error}"))?
    {
        ChannelCrossoverDiagnostic::AtInclination { theta_rad } => {
            let time_s = run
                .samples
                .windows(2)
                .find_map(|window| {
                    let before = window[0];
                    let after = window[1];
                    if before.theta_rad >= theta_rad && theta_rad >= after.theta_rad {
                        let fraction =
                            (before.theta_rad - theta_rad) / (before.theta_rad - after.theta_rad);
                        Some(before.time_s + fraction * (after.time_s - before.time_s))
                    } else {
                        None
                    }
                });
            Ok(match time_s {
                Some(time_s) => format!(
                    "{{\"class\":\"encoded-power-law-crossover\",\"reason\":\"both-channels-declared\",\"theta_rad\":{theta_rad:.17e},\"time_s\":{time_s:.17e}}}"
                ),
                None => format!(
                    "{{\"class\":\"encoded-power-law-crossover\",\"reason\":\"both-channels-declared\",\"theta_rad\":{theta_rad:.17e},\"time_status\":\"outside-retained-trajectory\"}}"
                ),
            })
        }
        ChannelCrossoverDiagnostic::NoneWithinInterval => Ok(
            "{\"class\":\"none\",\"reason\":\"no-encoded-power-crossing-within-validity-envelope\",\"theta_status\":\"none-within-validity-envelope\",\"time_status\":\"not-applicable\"}".to_owned(),
        ),
        ChannelCrossoverDiagnostic::NotComparable { reason } => Ok(format!(
            "{{\"class\":\"not-comparable\",\"reason\":\"{}\",\"theta_status\":\"unavailable\",\"time_status\":\"unavailable\"}}",
            match reason {
                ChannelCrossoverNotComparable::MissingDryContour => "missing-dry-contour",
                ChannelCrossoverNotComparable::MissingBildstenBoundaryLayer => {
                    "missing-bildsten-boundary-layer"
                }
            }
        )),
    }
}

fn flux_probe_record(
    radius_m: f64,
    thickness_m: f64,
    mass_kg: f64,
    transverse_inertia_kg_m2: f64,
    axial_inertia_kg_m2: f64,
    contact_state: fs_mbd::RigidBodyState,
) -> Result<String, String> {
    if !(transverse_inertia_kg_m2.is_finite()
        && transverse_inertia_kg_m2 > 0.0
        && axial_inertia_kg_m2.is_finite()
        && axial_inertia_kg_m2 > 0.0)
    {
        return Err("flux profile inertia was not positive and finite".to_owned());
    }
    let position = contact_state.pose().position_world();
    let linear_velocity = contact_state.linear_momentum_world().scale(1.0 / mass_kg);
    let angular_momentum_body = contact_state.angular_momentum_body();
    let angular_velocity_body = MbdVec3::new(
        angular_momentum_body.x / transverse_inertia_kg_m2,
        angular_momentum_body.y / transverse_inertia_kg_m2,
        angular_momentum_body.z / axial_inertia_kg_m2,
    );
    let orientation = contact_state.pose().orientation();
    let angular_velocity_world = orientation.rotate_body_to_world(angular_velocity_body);
    let normal_world = orientation.rotate_body_to_world(MbdVec3::new(0.0, 0.0, 1.0));
    let range = |minimum, maximum| {
        ClosedRange::try_new(minimum, maximum).map_err(|error| format!("flux range: {error}"))
    };
    let model = ReducedAeroModel::try_new(
        CorrelationIdentity::try_new(
            "campaign.reduced-exterior-wrench",
            "v1",
            "campaign/caller-declared-screening-correlation-v1",
        )
        .map_err(|error| format!("flux correlation: {error}"))?,
        ApplicabilityEnvelope {
            translational_reynolds: range(0.0, 1.0e9)?,
            rotational_reynolds: range(0.0, 1.0e9)?,
            relative_roughness: range(0.0, 0.1)?,
            maximum_tip_mach: 0.8,
        },
        CorrelationUncertainty {
            source_id: "campaign/caller-declared-screening-correlation-uncertainty-v1".to_owned(),
            coefficient_relative_half_width: 0.15,
        },
        ReducedAeroComponents {
            form_drag: Some(FormDrag { coefficient: 1.1 }),
            rotational_skin_friction: Some(RotationalSkinFriction { coefficient: 0.01 }),
            edge_flow: Some(EdgeFlow { coefficient: 0.04 }),
            ..ReducedAeroComponents::default()
        },
        &[
            ContributionFamily::TranslationalFormDrag,
            ContributionFamily::RotationalSkinFriction,
            ContributionFamily::EdgeFlow,
        ],
    )
    .map_err(|error| format!("flux model: {error}"))?;
    let input = ReducedAeroInput {
        world_frame_id: "campaign/world-inertial-v1".to_owned(),
        geometry: DiscGeometry {
            radius_m,
            exterior_thickness_m: thickness_m,
        },
        pose: DiscPose::try_new(Vec3::new(normal_world.x, normal_world.y, normal_world.z))
            .map_err(|error| format!("flux pose: {error}"))?,
        kinematics: BodyKinematics {
            reference_point_world_m: Vec3::new(position.x, position.y, position.z),
            linear_velocity_world_m_per_s: Vec3::new(
                linear_velocity.x,
                linear_velocity.y,
                linear_velocity.z,
            ),
            angular_velocity_world_rad_per_s: Vec3::new(
                angular_velocity_world.x,
                angular_velocity_world.y,
                angular_velocity_world.z,
            ),
        },
        gas: GasProperties::try_from(GasPropertyCard {
            source_id: "campaign/air-card-v1".to_owned(),
            density_kg_per_m3: Some(1.2),
            dynamic_viscosity_pa_s: Some(1.8e-5),
            speed_of_sound_m_per_s: Some(343.0),
            velocity_world_m_per_s: Vec3::ZERO,
        })
        .map_err(|error| format!("flux gas: {error}"))?,
        roughness: SurfaceRoughness {
            source_id: "campaign/exterior-roughness-card-v1".to_owned(),
            height_m: 1.0e-6,
        },
    };
    let wrench = model
        .evaluate(&input)
        .map_err(|error| format!("flux probe: {error}"))?;
    if wrench.receipt.relative_power_w > 1.0e-10 {
        return Err("exterior wrench violated passive relative-power identity".to_owned());
    }
    let work = WorkWindow::default()
        .record_once(1, 1.0e-4, &wrench)
        .map_err(|error| format!("flux work: {error}"))?;
    Ok(format!(
        concat!(
            "{{\"schema\":\"{SCHEMA}\",\"scenario\":\"reduced-exterior-wrench-passivity\",",
            "\"model\":\"{}\",\"source\":\"{}\",\"authority\":\"estimate-only-caller-declared-screening\",",
            "\"units\":\"{SI_UNITS}\",\"campaign_seed_manifest_ref\":\"campaign-complete\",\"inputs\":{{\"radius_m\":{:.17e},\"thickness_m\":{:.17e},\"profile_transverse_inertia_kg_m2\":{:.17e},\"profile_axial_inertia_kg_m2\":{:.17e},\"snapshot_position_z_m\":{:.17e},\"snapshot_linear_speed_m_per_s\":{:.17e},\"snapshot_angular_speed_rad_s\":{:.17e},\"gas_source_id\":\"campaign/air-card-v1\",\"roughness_source_id\":\"campaign/exterior-roughness-card-v1\",\"correlation_authority\":\"caller-declared-screening\"}},",
            "\"budget\":{{\"work_window_s\":1.0e-4}},\"terminal\":\"passivity-checked\",",
            "\"powers_w\":{{\"relative_dissipation\":{:.17e},\"body\":{:.17e}}},",
            "\"work_j\":{{\"relative_dissipation\":{:.17e},\"body\":{:.17e}}},",
            "\"residual\":{{\"passivity_w\":{:.17e}}},\"no_claim\":\"{NO_PHYSICAL_VALIDATION};",
            "one-way-contact-state-snapshot;screening-correlation-not-cfd-or-gas-film\"}}"
        ),
        wrench.correlation.id,
        wrench.correlation.source_id,
        radius_m,
        thickness_m,
        transverse_inertia_kg_m2,
        axial_inertia_kg_m2,
        position.z,
        (linear_velocity.x.mul_add(
            linear_velocity.x,
            linear_velocity
                .y
                .mul_add(linear_velocity.y, linear_velocity.z * linear_velocity.z),
        ))
        .sqrt(),
        (angular_velocity_world.x.mul_add(
            angular_velocity_world.x,
            angular_velocity_world.y.mul_add(
                angular_velocity_world.y,
                angular_velocity_world.z * angular_velocity_world.z,
            ),
        ))
        .sqrt(),
        wrench.receipt.dissipated_relative_power_w,
        wrench.receipt.body_power_w,
        work.relative_dissipation_j,
        work.body_work_j,
        wrench.receipt.relative_power_w,
        SCHEMA = SCHEMA,
        SI_UNITS = SI_UNITS,
        NO_PHYSICAL_VALIDATION = NO_PHYSICAL_VALIDATION,
    ))
}
