//! G0/G1/G3 checks for the bounded isothermal gas-film foundation.

use fs_flux::{
    ContactExclusionMask, GasFilmApplicability, GasFilmBoundaryTopology, GasFilmBudget,
    GasFilmError, GasFilmGrid1d, GasFilmIdentity, GasFilmInput, GasFilmInputAuthority,
    GasFilmUncertainty, IsothermalIdealGas, MovingWallInput, RoughnessPolicy, SlipPolicy,
    isothermal_compressible_reynolds_model_id, solve_isothermal_gas_film_1d,
};

fn fixture(cells: usize) -> GasFilmInput {
    let pressure = 101_325.0;
    let gas_constant = 287.05;
    let temperature = 300.0;
    GasFilmInput {
        identity: GasFilmIdentity {
            case_id: "gas-film-synthetic-fixture-v1".to_owned(),
            model_id: isothermal_compressible_reynolds_model_id().to_owned(),
            gas_species_id: "synthetic-air".to_owned(),
            eos_id: "ideal-gas-isothermal-v1".to_owned(),
            viscosity_source_id: "synthetic-viscosity-v1".to_owned(),
            thermal_model_id: "isothermal-v1".to_owned(),
            frame_id: "fixture-line-x".to_owned(),
            deterministic_seed: 7,
            authority: GasFilmInputAuthority::SyntheticFixture,
        },
        gas: IsothermalIdealGas {
            specific_gas_constant_j_kg_k: gas_constant,
            temperature_k: temperature,
            dynamic_viscosity_pa_s: 1.8e-5,
            declared_density_kg_m3: pressure / (gas_constant * temperature),
            declared_specific_enthalpy_j_kg: 300_000.0,
        },
        grid: GasFilmGrid1d {
            length_m: 1.0e-2,
            gap_m: vec![10.0e-6; cells],
            contact_exclusion: ContactExclusionMask {
                excluded: vec![false; cells],
            },
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
            maximum_gap_slope: 0.1,
            speed_of_sound_m_per_s: 347.0,
            maximum_mach_number: 0.3,
        },
        uncertainty: GasFilmUncertainty {
            viscosity_relative_bound: 0.0,
            gap_relative_bound: 0.0,
            pressure_relative_bound: 0.0,
        },
        wall_motion: MovingWallInput {
            lower_tangential_velocity_m_per_s: 0.0,
            upper_tangential_velocity_m_per_s: 0.0,
            gap_rate_m_per_s: 0.0,
        },
        initial_absolute_pressure_pa: pressure,
        timestep_s: 1.0e-5,
        budget: GasFilmBudget {
            maximum_iterations: 2_000,
            mass_residual_tolerance_kg_m2_s: 1.0e-9,
            relaxation: 0.8,
        },
    }
}

fn active_pressure(step: &fs_flux::GasFilmStep) -> Vec<f64> {
    step.absolute_pressure_pa
        .iter()
        .map(|pressure| pressure.expect("fixture active cell"))
        .collect()
}

#[test]
fn g0_uniform_sealed_equilibrium_is_exact_and_replayable() {
    let input = fixture(8);
    let first = solve_isothermal_gas_film_1d(&input, None).expect("uniform equilibrium admitted");
    let second = solve_isothermal_gas_film_1d(&input, None).expect("same replay admitted");
    assert_eq!(first, second, "same synthetic request replays bit-for-bit");
    for pressure in active_pressure(&first) {
        assert_eq!(pressure, input.initial_absolute_pressure_pa);
    }
    assert_eq!(first.receipt.max_mass_residual_kg_m2_s, 0.0);
    assert_eq!(first.receipt.mass_closure_residual_kg_per_m_s, 0.0);
    assert_eq!(first.receipt.input_uncertainty, input.uncertainty);
}

#[test]
fn g1_planar_couette_is_uniform_and_has_analytic_shear_and_heat() {
    let mut input = fixture(6);
    input.wall_motion.upper_tangential_velocity_m_per_s = 2.0;
    let step = solve_isothermal_gas_film_1d(&input, None).expect("Couette limit admitted");
    let expected_shear = input.gas.dynamic_viscosity_pa_s * 2.0 / input.grid.gap_m[0];
    for pressure in active_pressure(&step) {
        assert!((pressure - input.initial_absolute_pressure_pa).abs() < 1.0e-9);
    }
    for shear in step.receipt.upper_wall_shear_pa {
        assert!((shear.expect("active") - expected_shear).abs() < 1.0e-12);
    }
    let expected_power = expected_shear * 2.0 * input.grid.length_m;
    assert!((step.receipt.wall_power_to_gas_w_per_m - expected_power).abs() < 1.0e-12);
    assert_eq!(
        step.receipt.wall_power_to_gas_w_per_m,
        step.receipt.viscous_heat_w_per_m
    );
}

#[test]
fn g3_uniform_couette_spatial_refinement_preserves_the_analytic_limit() {
    let mut coarse = fixture(4);
    coarse.wall_motion.upper_tangential_velocity_m_per_s = 1.25;
    let mut fine = fixture(32);
    fine.wall_motion.upper_tangential_velocity_m_per_s = 1.25;
    let coarse_step = solve_isothermal_gas_film_1d(&coarse, None).expect("coarse admitted");
    let fine_step = solve_isothermal_gas_film_1d(&fine, None).expect("fine admitted");
    let expected = coarse.gas.dynamic_viscosity_pa_s * 1.25 / coarse.grid.gap_m[0];
    assert!(
        (coarse_step.receipt.upper_wall_shear_pa[0].expect("active") - expected).abs() < 1.0e-12,
        "coarse shear={} expected={expected}",
        coarse_step.receipt.upper_wall_shear_pa[0].expect("active")
    );
    assert!(
        (fine_step.receipt.upper_wall_shear_pa[0].expect("active") - expected).abs() < 1.0e-12,
        "fine shear={} expected={expected}",
        fine_step.receipt.upper_wall_shear_pa[0].expect("active")
    );
}

#[test]
fn g1_uniform_sealed_squeeze_preserves_mass_and_matches_uniform_limit() {
    let mut input = fixture(10);
    input.wall_motion.gap_rate_m_per_s = -1.0e-3;
    let step = solve_isothermal_gas_film_1d(&input, None).expect("uniform squeeze admitted");
    let old_gap = input.grid.gap_m[0] - input.wall_motion.gap_rate_m_per_s * input.timestep_s;
    let expected_pressure = input.initial_absolute_pressure_pa * old_gap / input.grid.gap_m[0];
    for pressure in active_pressure(&step) {
        assert!(
            (pressure - expected_pressure).abs() / expected_pressure < 1.0e-10,
            "pressure={pressure:e} expected={expected_pressure:e}"
        );
    }
    assert!(step.receipt.mass_closure_residual_kg_per_m_s.abs() < 1.0e-11);
    assert_eq!(step.receipt.left_boundary_outward_mass_flux_kg_per_m_s, 0.0);
    assert_eq!(
        step.receipt.right_boundary_outward_mass_flux_kg_per_m_s,
        0.0
    );
    let expected_normal_power =
        -expected_pressure * input.wall_motion.gap_rate_m_per_s * input.grid.length_m;
    assert!(
        (step.receipt.normal_gap_power_to_gas_w_per_m - expected_normal_power).abs()
            / expected_normal_power
            < 1.0e-10,
        "normal power={} expected={expected_normal_power:e}",
        step.receipt.normal_gap_power_to_gas_w_per_m
    );
    assert_eq!(
        step.receipt.wall_power_to_gas_w_per_m,
        step.receipt.normal_gap_power_to_gas_w_per_m
    );
}

#[test]
fn g1_uniform_squeeze_time_refinement_preserves_the_exact_limit() {
    let mut coarse = fixture(8);
    coarse.wall_motion.gap_rate_m_per_s = -1.0e-3;
    let mut fine = coarse.clone();
    fine.timestep_s *= 0.5;
    let coarse_step = solve_isothermal_gas_film_1d(&coarse, None).expect("coarse admitted");
    let fine_step = solve_isothermal_gas_film_1d(&fine, None).expect("fine admitted");
    for (input, step) in [(&coarse, &coarse_step), (&fine, &fine_step)] {
        let old_gap = input.grid.gap_m[0] - input.wall_motion.gap_rate_m_per_s * input.timestep_s;
        let expected = input.initial_absolute_pressure_pa * old_gap / input.grid.gap_m[0];
        let actual = active_pressure(step)[0];
        assert!(
            (actual - expected).abs() / expected < 1.0e-10,
            "dt={} pressure={actual:e} expected={expected:e}",
            input.timestep_s
        );
    }
}

#[test]
fn g1_open_poiseuille_has_monotone_pressure_and_boundary_mass_closure() {
    let mut input = fixture(32);
    input.boundary = GasFilmBoundaryTopology::Open {
        left_absolute_pressure_pa: 120_000.0,
        right_absolute_pressure_pa: 100_000.0,
    };
    input.budget.maximum_iterations = 20_000;
    input.budget.mass_residual_tolerance_kg_m2_s = 1.0e-7;
    let step = solve_isothermal_gas_film_1d(&input, None).expect("Poiseuille limit admitted");
    let pressure = active_pressure(&step);
    assert!(pressure.windows(2).all(|pair| pair[0] > pair[1]));
    assert!(step.receipt.left_boundary_outward_mass_flux_kg_per_m_s < 0.0);
    assert!(step.receipt.right_boundary_outward_mass_flux_kg_per_m_s > 0.0);
    assert!(step.receipt.mass_closure_residual_kg_per_m_s.abs() < 5.0e-9);
}

#[test]
fn g0_vented_and_open_boundaries_close_storage_with_named_fluxes() {
    let mut vented = fixture(7);
    vented.boundary = GasFilmBoundaryTopology::Vented {
        cell_index: 3,
        absolute_pressure_pa: 90_000.0,
    };
    let step = solve_isothermal_gas_film_1d(&vented, None).expect("vented topology admitted");
    assert_eq!(step.absolute_pressure_pa[3], Some(90_000.0));
    assert!(step.receipt.mass_closure_residual_kg_per_m_s.abs() < 1.0e-9);
    assert!(step.receipt.vent_outward_mass_flux_kg_per_m_s.is_finite());
    assert_eq!(
        step.receipt.vent_outward_enthalpy_flux_w_per_m,
        step.receipt.vent_outward_mass_flux_kg_per_m_s * vented.gas.declared_specific_enthalpy_j_kg
    );
}

#[test]
fn g3_reversing_relative_wall_motion_preserves_pressure_and_reverses_shear() {
    let mut forward = fixture(8);
    forward.wall_motion.lower_tangential_velocity_m_per_s = -1.0;
    forward.wall_motion.upper_tangential_velocity_m_per_s = 2.0;
    let mut reversed = forward.clone();
    reversed.wall_motion.lower_tangential_velocity_m_per_s = 1.0;
    reversed.wall_motion.upper_tangential_velocity_m_per_s = -2.0;
    let positive = solve_isothermal_gas_film_1d(&forward, None).expect("forward admitted");
    let negative = solve_isothermal_gas_film_1d(&reversed, None).expect("reversed admitted");
    assert_eq!(active_pressure(&positive), active_pressure(&negative));
    assert_eq!(
        positive.receipt.upper_wall_shear_pa[0].expect("active"),
        -negative.receipt.upper_wall_shear_pa[0].expect("active")
    );
    assert_eq!(
        positive.receipt.wall_power_to_gas_w_per_m,
        negative.receipt.wall_power_to_gas_w_per_m
    );
}

#[test]
fn g3_pressure_and_density_scaling_hold_for_uniform_sealed_state() {
    let base = fixture(5);
    let mut scaled = base.clone();
    scaled.initial_absolute_pressure_pa *= 2.0;
    scaled.gas.declared_density_kg_m3 *= 2.0;
    let base_step = solve_isothermal_gas_film_1d(&base, None).expect("base admitted");
    let scaled_step = solve_isothermal_gas_film_1d(&scaled, None).expect("scaled admitted");
    assert_eq!(
        active_pressure(&scaled_step),
        active_pressure(&base_step)
            .into_iter()
            .map(|value| 2.0 * value)
            .collect::<Vec<_>>()
    );
    assert_eq!(
        scaled_step.receipt.gas_mass_kg_per_m,
        2.0 * base_step.receipt.gas_mass_kg_per_m
    );
}

#[test]
fn g0_contact_exclusion_is_explicit_and_never_fills_a_gap() {
    let mut input = fixture(6);
    input.grid.contact_exclusion.excluded[5] = true;
    let step =
        solve_isothermal_gas_film_1d(&input, None).expect("static terminal exclusion admitted");
    assert_eq!(step.absolute_pressure_pa[5], None);
    assert_eq!(step.receipt.upper_wall_shear_pa[5], None);
    assert_eq!(step.checkpoint.active_cells, 5);

    input.grid.contact_exclusion.excluded = vec![false, true, false, false, false, false];
    assert_eq!(
        solve_isothermal_gas_film_1d(&input, None),
        Err(GasFilmError::TopologyChangeUnavailable)
    );
}

#[test]
fn hostile_rarefaction_negative_state_slip_and_roughness_refuse() {
    let mut rarefied = fixture(4);
    rarefied.grid.gap_m[0] = 1.0e-7;
    assert!(matches!(
        solve_isothermal_gas_film_1d(&rarefied, None),
        Err(GasFilmError::Unavailable {
            reason: "rarefied-knudsen-outside-continuum-envelope"
        })
    ));

    let mut negative = fixture(4);
    negative.initial_absolute_pressure_pa = -1.0;
    assert_eq!(
        solve_isothermal_gas_film_1d(&negative, None),
        Err(GasFilmError::InvalidInput {
            field: "initial_absolute_pressure_pa"
        })
    );

    let mut slip = fixture(4);
    slip.slip_policy = SlipPolicy::RarefiedSlipRequested {
        source_id: "unadmitted-slip-card".to_owned(),
    };
    assert!(matches!(
        solve_isothermal_gas_film_1d(&slip, None),
        Err(GasFilmError::Unavailable {
            reason: "slip-or-rarefied-model-not-implemented"
        })
    ));

    let mut roughness = fixture(4);
    roughness.roughness_policy = RoughnessPolicy::Unresolved {
        source_id: "unknown-roughness".to_owned(),
    };
    assert!(matches!(
        solve_isothermal_gas_film_1d(&roughness, None),
        Err(GasFilmError::Unavailable {
            reason: "roughness-model-not-admitted"
        })
    ));

    let oversized = fixture(4_097);
    assert!(matches!(
        solve_isothermal_gas_film_1d(&oversized, None),
        Err(GasFilmError::Unavailable {
            reason: "gas-film-grid-exceeds-bounded-cell-cap"
        })
    ));
}

#[test]
fn checkpoint_replay_and_budget_refusal_are_deterministic() {
    let input = fixture(8);
    let first = solve_isothermal_gas_film_1d(&input, None).expect("first step");
    let checkpoint_before_refusal = first.checkpoint.clone();
    let replay_a = solve_isothermal_gas_film_1d(&input, Some(&first.checkpoint)).expect("replay a");
    let replay_b = solve_isothermal_gas_film_1d(&input, Some(&first.checkpoint)).expect("replay b");
    assert_eq!(replay_a, replay_b);

    let mut mismatched_case = input.clone();
    mismatched_case.identity.case_id = "different-synthetic-case".to_owned();
    assert_eq!(
        solve_isothermal_gas_film_1d(&mismatched_case, Some(&first.checkpoint)),
        Err(GasFilmError::CheckpointMismatch { field: "case_id" })
    );

    let mut constrained = fixture(16);
    constrained.boundary = GasFilmBoundaryTopology::Open {
        left_absolute_pressure_pa: 130_000.0,
        right_absolute_pressure_pa: 90_000.0,
    };
    constrained.budget.maximum_iterations = 1;
    constrained.budget.mass_residual_tolerance_kg_m2_s = 1.0e-30;
    assert!(matches!(
        solve_isothermal_gas_film_1d(&constrained, None),
        Err(GasFilmError::IterationBudgetExceeded { .. })
    ));
    assert_eq!(
        first.checkpoint, checkpoint_before_refusal,
        "a refused solve cannot mutate a published checkpoint"
    );
}
