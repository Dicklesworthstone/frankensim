//! G0/G1/G3 checks for the declared independent-sector tilted gas-film adapter.
//!
//! These verify code-level limits of declared equations only. They do not
//! establish an air-dominance mechanism, cross-sector flow, or any video/device outcome.

use fs_euler_disc_e2e::air::{
    AirFilmDiscretization, AirFilmError, AirFilmIdentity, AirVec3, ContactExclusion,
    PrescribedPlaneBase, TILTED_DISC_GAS_FILM_ADAPTER_ID, TiltedDiscAirFilmInput,
    TiltedDiscKinematics, sample_tilted_disc_gap, solve_tilted_disc_air_film,
};
use fs_flux::{
    GasFilmApplicability, GasFilmBoundaryTopology, GasFilmBudget, GasFilmInputAuthority,
    GasFilmUncertainty, IsothermalIdealGas, RoughnessPolicy, SlipPolicy,
};

fn fixture() -> TiltedDiscAirFilmInput {
    let pressure = 101_325.0;
    let gas_constant = 287.05;
    let temperature = 300.0;
    TiltedDiscAirFilmInput {
        identity: AirFilmIdentity {
            case_id: "air-film-synthetic-fixture-v1".to_owned(),
            adapter_model_id: TILTED_DISC_GAS_FILM_ADAPTER_ID.to_owned(),
            frame_id: "fixture-world-z-up".to_owned(),
            base_motion_id: "prescribed-horizontal-plane-v1".to_owned(),
            gas_species_id: "synthetic-air".to_owned(),
            eos_id: "synthetic-isothermal-ideal-gas-v1".to_owned(),
            viscosity_source_id: "synthetic-viscosity-v1".to_owned(),
            thermal_model_id: "synthetic-isothermal-v1".to_owned(),
            configuration_id: "synthetic-air-film-config-v1".to_owned(),
            deterministic_seed: 17,
            authority: GasFilmInputAuthority::SyntheticFixture,
        },
        disc_radius_m: 1.0e-2,
        disc_half_thickness_m: 1.0e-4,
        disc: TiltedDiscKinematics {
            center_world_m: AirVec3::new(0.0, 0.0, 1.0e-3),
            normal_away_from_base_world: AirVec3::new(0.0, 0.0, 1.0),
            center_velocity_world_m_per_s: AirVec3::ZERO,
            angular_velocity_world_rad_per_s: AirVec3::ZERO,
        },
        base: PrescribedPlaneBase {
            height_m: 0.0,
            vertical_velocity_m_per_s: 0.0,
        },
        discretization: AirFilmDiscretization {
            azimuthal_sectors: 4,
            radial_cells: 4,
        },
        contact_exclusion: ContactExclusion {
            handoff_gap_m: 1.0e-6,
        },
        gas: IsothermalIdealGas {
            specific_gas_constant_j_kg_k: gas_constant,
            temperature_k: temperature,
            dynamic_viscosity_pa_s: 1.8e-5,
            declared_density_kg_m3: pressure / (gas_constant * temperature),
            declared_specific_enthalpy_j_kg: 300_000.0,
        },
        boundary: GasFilmBoundaryTopology::Sealed,
        slip_policy: SlipPolicy::NoSlipContinuum {
            source_id: "continuum-no-slip-v1".to_owned(),
        },
        roughness_policy: RoughnessPolicy::ResolvedSmooth {
            source_id: "synthetic-smooth-wall-v1".to_owned(),
            maximum_roughness_m: 1.0e-8,
        },
        applicability: GasFilmApplicability {
            mean_free_path_m: 65.0e-9,
            maximum_knudsen_number: 0.01,
            maximum_gap_slope: 0.2,
            speed_of_sound_m_per_s: 347.0,
            maximum_mach_number: 0.3,
        },
        uncertainty: GasFilmUncertainty {
            viscosity_relative_bound: 0.0,
            gap_relative_bound: 0.0,
            pressure_relative_bound: 0.0,
        },
        initial_absolute_pressure_pa: pressure,
        gauge_reference_absolute_pressure_pa: pressure,
        timestep_s: 1.0e-6,
        budget: GasFilmBudget {
            maximum_iterations: 2_000,
            mass_residual_tolerance_kg_m2_s: 1.0e-9,
            relaxation: 0.8,
        },
    }
}

#[test]
fn g1_parallel_plate_squeeze_has_positive_pressure_work_and_mass_closure() {
    let mut input = fixture();
    input.disc.center_velocity_world_m_per_s.z = -1.0e-3;
    let step = solve_tilted_disc_air_film(&input, None).expect("parallel squeeze admitted");
    assert!(
        step.receipt.wall_power_to_gas_w > 0.0,
        "wall power={}",
        step.receipt.wall_power_to_gas_w
    );
    assert!(
        step.receipt.wrench.force_world_n.z > 0.0,
        "pressure force={:?}",
        step.receipt.wrench.force_world_n
    );
    assert!(
        step.receipt.mass_closure_residual_kg_per_s.abs() < 1.0e-12,
        "closure={}",
        step.receipt.mass_closure_residual_kg_per_s
    );
}

#[test]
fn g1_parallel_plate_radial_shear_has_positive_heat_and_reaction_force() {
    let mut input = fixture();
    input.disc.center_velocity_world_m_per_s.x = 1.0;
    let step = solve_tilted_disc_air_film(&input, None).expect("radial-strip shear admitted");
    assert!(
        step.receipt.viscous_heat_w >= 0.0,
        "heat={}",
        step.receipt.viscous_heat_w
    );
    assert!(step.receipt.wrench.force_world_n.x.is_finite());
    assert!(step.receipt.wrench.force_world_n.y.is_finite());
}

#[test]
fn g3_tilted_sector_permutation_preserves_sorted_gaps_and_vertical_load() {
    let mut first = fixture();
    first.disc.normal_away_from_base_world = AirVec3::new(0.02, 0.0, (1.0 - 0.0004_f64).sqrt());
    let mut rotated = first.clone();
    rotated.disc.normal_away_from_base_world = AirVec3::new(0.0, 0.02, (1.0 - 0.0004_f64).sqrt());
    let a = solve_tilted_disc_air_film(&first, None).expect("first tilt admitted");
    let b = solve_tilted_disc_air_film(&rotated, None).expect("quarter-turn tilt admitted");
    let mut gaps_a = a
        .samples
        .iter()
        .map(|sample| sample.gap_m)
        .collect::<Vec<_>>();
    let mut gaps_b = b
        .samples
        .iter()
        .map(|sample| sample.gap_m)
        .collect::<Vec<_>>();
    gaps_a.sort_by(f64::total_cmp);
    gaps_b.sort_by(f64::total_cmp);
    assert_eq!(
        gaps_a, gaps_b,
        "sector permutation must preserve sampled gap multiset"
    );
    assert!(
        (a.receipt.wrench.force_world_n.z - b.receipt.wrench.force_world_n.z).abs() < 1.0e-10,
        "vertical loads: {} vs {}",
        a.receipt.wrench.force_world_n.z,
        b.receipt.wrench.force_world_n.z
    );
}

#[test]
fn g1_half_thickness_moves_horizontal_and_tilted_face_gap_exactly() {
    let input = fixture();
    let horizontal = sample_tilted_disc_gap(&input, 0, 0).expect("horizontal sample");
    assert_eq!(
        horizontal.gap_m,
        input.disc.center_world_m.z - input.disc_half_thickness_m
    );
    let mut tilted = input.clone();
    tilted.disc.normal_away_from_base_world = AirVec3::new(0.02, 0.0, (1.0 - 0.0004_f64).sqrt());
    let sample = sample_tilted_disc_gap(&tilted, 0, 3).expect("tilted sample");
    assert!(
        sample.lever_arm_world_m.z < 0.0,
        "gas-facing arm must include negative face offset"
    );
    assert_eq!(
        sample.gap_m,
        tilted.disc.center_world_m.z + sample.lever_arm_world_m.z - tilted.base.height_m
    );
}

#[test]
fn g1_normal_axis_spin_has_circumferential_heat_and_opposing_torque_without_radial_drive() {
    let mut input = fixture();
    input.disc.angular_velocity_world_rad_per_s = AirVec3::new(0.0, 0.0, 40.0);
    let sample = sample_tilted_disc_gap(&input, 0, 3).expect("spin sample");
    assert!(
        sample.radial_relative_velocity_m_per_s.abs() < 1.0e-14,
        "normal-axis spin has no radial strip drive: {}",
        sample.radial_relative_velocity_m_per_s
    );
    assert!(sample.circumferential_relative_velocity_m_per_s.abs() > 0.0);
    let step = solve_tilted_disc_air_film(&input, None).expect("spin admitted");
    assert!(step.receipt.circumferential_couette_heat_w > 0.0);
    assert!(
        step.receipt.wrench.moment_about_com_world_n_m.z < 0.0,
        "gas torque must oppose positive spin"
    );
}

#[test]
fn g0_equal_area_surrogate_reports_exact_negative_quarter_first_moment_defect() {
    let input = fixture();
    let step = solve_tilted_disc_air_film(&input, None).expect("surrogate admitted");
    let exact = 2.0 * core::f64::consts::PI * input.disc_radius_m.powi(3) / 3.0;
    let expected = -0.25 * exact;
    assert_eq!(
        step.receipt.equal_area_strip_width_m,
        0.5 * input.disc_radius_m
            * (2.0 * core::f64::consts::PI / input.discretization.azimuthal_sectors as f64)
    );
    assert!(
        (step.receipt.signed_first_radial_moment_discrepancy_m3 - expected).abs() < 1.0e-20,
        "reported={} expected={expected:e}",
        step.receipt.signed_first_radial_moment_discrepancy_m3
    );
}

#[test]
fn g0_action_work_sign_and_checkpoint_replay_are_deterministic() {
    let mut input = fixture();
    input.disc.center_velocity_world_m_per_s.z = -2.0e-3;
    let first = solve_tilted_disc_air_film(&input, None).expect("first admitted");
    let mut next = input.clone();
    next.disc.center_world_m.z += input.disc.center_velocity_world_m_per_s.z * input.timestep_s;
    let replay_a = solve_tilted_disc_air_film(&next, Some(&first.checkpoint)).expect("replay a");
    let replay_b = solve_tilted_disc_air_film(&next, Some(&first.checkpoint)).expect("replay b");
    assert_eq!(
        replay_a, replay_b,
        "identity-bound sector state replays deterministically"
    );
    assert!(
        first.receipt.wrench.force_world_n.z * input.disc.center_velocity_world_m_per_s.z < 0.0,
        "gas force must oppose squeeze velocity"
    );
}

#[test]
fn g0_changing_tilt_advances_nonuniform_gap_profile_across_checkpoint() {
    let mut first_input = fixture();
    let normal = AirVec3::new(0.02, 0.0, (1.0 - 0.0004_f64).sqrt());
    let angular_velocity = AirVec3::new(0.0, 5.0, 0.0);
    first_input.disc.normal_away_from_base_world = normal;
    first_input.disc.angular_velocity_world_rad_per_s = angular_velocity;
    let first = solve_tilted_disc_air_film(&first_input, None).expect("initial tilted film step");

    let mut second_input = first_input.clone();
    let candidate = AirVec3::new(
        normal.x + angular_velocity.y * normal.z * first_input.timestep_s,
        normal.y,
        normal.z - angular_velocity.y * normal.x * first_input.timestep_s,
    );
    let norm =
        (candidate.x * candidate.x + candidate.y * candidate.y + candidate.z * candidate.z).sqrt();
    second_input.disc.normal_away_from_base_world =
        AirVec3::new(candidate.x / norm, candidate.y / norm, candidate.z / norm);
    let second = solve_tilted_disc_air_film(&second_input, Some(&first.checkpoint))
        .expect("changing wedge slope must restart");
    assert!(
        second
            .checkpoint
            .sectors
            .iter()
            .all(|sector| sector.step_index == 2)
    );
    assert_eq!(
        second,
        solve_tilted_disc_air_film(&second_input, Some(&first.checkpoint))
            .expect("changing-slope replay")
    );
}

#[test]
fn g0_contact_and_boundary_topology_are_explicit() {
    let mut contact = fixture();
    contact.disc.center_world_m.z = contact.disc_half_thickness_m + 1.1e-6;
    contact.disc.normal_away_from_base_world =
        AirVec3::new(2.0e-5, 0.0, (1.0 - 4.0e-10_f64).sqrt());
    let step =
        solve_tilted_disc_air_film(&contact, None).expect("outer-rim contact handoff admitted");
    assert!(step.samples.iter().any(|sample| sample.excluded_by_contact));

    let mut vented = fixture();
    vented.boundary = GasFilmBoundaryTopology::Vented {
        cell_index: 1,
        absolute_pressure_pa: 90_000.0,
    };
    let vented_step = solve_tilted_disc_air_film(&vented, None).expect("vented strips admitted");
    assert!(vented_step.receipt.outward_enthalpy_flux_w.is_finite());

    let mut open = fixture();
    open.boundary = GasFilmBoundaryTopology::Open {
        left_absolute_pressure_pa: 110_000.0,
        right_absolute_pressure_pa: 90_000.0,
    };
    open.budget.mass_residual_tolerance_kg_m2_s = 1.0e-7;
    let open_step = solve_tilted_disc_air_film(&open, None).expect("open strips admitted");
    assert!(open_step.receipt.outward_mass_flux_kg_per_s.is_finite());
}

#[test]
fn g3_radius_and_pressure_scaling_are_declared_model_limits() {
    let base = solve_tilted_disc_air_film(&fixture(), None).expect("base admitted");
    let mut doubled = fixture();
    doubled.initial_absolute_pressure_pa *= 2.0;
    doubled.gas.declared_density_kg_m3 *= 2.0;
    let scaled = solve_tilted_disc_air_film(&doubled, None).expect("pressure-scaled admitted");
    assert!(
        (scaled.receipt.gas_mass_kg - 2.0 * base.receipt.gas_mass_kg).abs() < 1.0e-14,
        "mass scaling: {} vs {}",
        scaled.receipt.gas_mass_kg,
        base.receipt.gas_mass_kg
    );
}

#[test]
fn hostile_invalid_contact_rarefied_roughness_and_topology_refuse() {
    let mut rarefied = fixture();
    rarefied.applicability.mean_free_path_m = 1.0e-5;
    assert!(matches!(
        solve_tilted_disc_air_film(&rarefied, None),
        Err(AirFilmError::GasFilmRefusal { .. })
    ));

    let mut rough = fixture();
    rough.roughness_policy = RoughnessPolicy::Unresolved {
        source_id: "no-roughness-card".to_owned(),
    };
    assert!(matches!(
        solve_tilted_disc_air_film(&rough, None),
        Err(AirFilmError::GasFilmRefusal { .. })
    ));

    let mut molecular = fixture();
    molecular.slip_policy = SlipPolicy::RarefiedSlipRequested {
        source_id: "molecular-request".to_owned(),
    };
    assert!(matches!(
        solve_tilted_disc_air_film(&molecular, None),
        Err(AirFilmError::GasFilmRefusal { .. })
    ));

    let mut invalid = fixture();
    invalid.disc.center_world_m.z = f64::NAN;
    assert_eq!(
        solve_tilted_disc_air_film(&invalid, None),
        Err(AirFilmError::InvalidInput {
            field: "disc.center_world_m"
        })
    );

    let valid = solve_tilted_disc_air_film(&fixture(), None).expect("checkpoint fixture");
    let mut other = fixture();
    other.identity.case_id = "different-synthetic-case".to_owned();
    assert_eq!(
        solve_tilted_disc_air_film(&other, Some(&valid.checkpoint)),
        Err(AirFilmError::CheckpointMismatch { field: "case_id" })
    );

    let mut provenance_mutation = fixture();
    provenance_mutation.identity.viscosity_source_id = "other-viscosity-card".to_owned();
    provenance_mutation.identity.configuration_id = "different-config".to_owned();
    assert_eq!(
        solve_tilted_disc_air_film(&provenance_mutation, Some(&valid.checkpoint)),
        Err(AirFilmError::CheckpointMismatch {
            field: "configuration_id"
        })
    );
}
