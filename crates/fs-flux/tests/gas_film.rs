//! G0/G1/G3 checks for the bounded isothermal gas-film foundation.

use fs_flux::{
    AxisymmetricCapillaryVent, AxisymmetricGasFilmReservoir, CapillaryVent, ContactExclusionMask,
    GasFilmApplicability, GasFilmBoundaryTopology, GasFilmBudget, GasFilmError, GasFilmGrid1d,
    GasFilmIdentity, GasFilmInput, GasFilmInputAuthority, GasFilmReservoir, GasFilmUncertainty,
    IsothermalIdealGas, MovingWallInput, RoughnessPolicy, SlipPolicy,
    isothermal_compressible_reynolds_model_id, solve_isothermal_axisymmetric_gas_film_1d,
    solve_isothermal_gas_film_1d,
};

fn radial_connection_fixture(
    cells: usize,
    pressure: f64,
    count: u32,
    volume: Option<f64>,
) -> GasFilmInput {
    let mut input = fixture(cells);
    let capillary = AxisymmetricCapillaryVent {
        radius_m: 1.0e-5,
        length_m: 1.0e-2,
        capillary_count: count,
    };
    input.boundary = match volume {
        None => GasFilmBoundaryTopology::AxisymmetricCapillaryVented {
            cell_index: cells - 1,
            absolute_pressure_pa: pressure,
            capillary,
        },
        Some(volume_m3) => GasFilmBoundaryTopology::AxisymmetricCapillaryReservoir {
            cell_index: cells - 1,
            capillary,
            reservoir: AxisymmetricGasFilmReservoir {
                volume_m3,
                initial_absolute_pressure_pa: pressure,
            },
        },
    };
    input.budget.maximum_iterations = 12;
    input.budget.relaxation = 1.0;
    input.budget.mass_residual_tolerance_kg_m2_s = 1.0e-12;
    input
}

#[test]
fn g1_radial_capillary_counts_total_tubes_and_conserves_each_annulus() {
    let pi = std::f64::consts::PI;
    for inner in [0.0, 0.02] {
        for bath in [100_325.0, 102_325.0] {
            for count in [1, 4] {
                for radius in [1.0e-5, 2.0e-5] {
                    let mut input = radial_connection_fixture(2, bath, count, None);
                    input.timestep_s = 1.0e-4;
                    if let GasFilmBoundaryTopology::AxisymmetricCapillaryVented {
                        capillary, ..
                    } = &mut input.boundary
                    {
                        capillary.radius_m = radius;
                    }
                    let step =
                        solve_isothermal_axisymmetric_gas_film_1d(&input, inner, None).unwrap();
                    let p = &step.checkpoint.active_absolute_pressure_pa;
                    let rt = input.gas.specific_gas_constant_j_kg_k * input.gas.temperature_k;
                    let mu = input.gas.dynamic_viscosity_pa_s;
                    let h = input.grid.gap_m[0];
                    let dr = input.grid.length_m / 2.0;
                    let expected_flow =
                        f64::from(count) * pi * radius.powi(4) * (p[1] * p[1] - bath * bath)
                            / (16.0 * mu * 0.01 * rt);
                    let flow = step.receipt.vent_outward_mass_flow_kg_per_s;
                    assert!((flow / expected_flow - 1.0).abs() < 1.0e-12);
                    assert!(flow * (input.initial_absolute_pressure_pa - bath) > 0.0);
                    // Integrate the radial continuum resistance between sample
                    // radii, independently of the implementation's face metric.
                    let internal = pi * h.powi(3) * (p[0] * p[0] - p[1] * p[1])
                        / (12.0 * mu * rt * ((inner + 1.5 * dr) / (inner + 0.5 * dr)).ln());
                    for cell in 0..2 {
                        let r0 = inner + cell as f64 * dr;
                        let area = pi * ((r0 + dr).powi(2) - r0.powi(2));
                        let storage = (p[cell] - input.initial_absolute_pressure_pa) * h * area
                            / (rt * input.timestep_s);
                        let outward = if cell == 0 { internal } else { flow - internal };
                        assert!((storage + outward).abs() < 3.0e-16);
                    }
                    let vent = step.receipt.capillary_vent.unwrap();
                    let expected_work = expected_flow * rt * (p[1] / bath).ln();
                    assert!((vent.viscous_dissipation_w / expected_work - 1.0).abs() < 1.0e-12);
                    assert!(vent.viscous_dissipation_w >= 0.0);
                    let expected_re = 2.0 * flow.abs() / (f64::from(count) * pi * radius * mu);
                    assert!((vent.reynolds_number / expected_re - 1.0).abs() < 1.0e-14);
                    assert_eq!(
                        step.receipt.vent_outward_enthalpy_power_w,
                        flow * input.gas.declared_specific_enthalpy_j_kg
                    );
                    assert_eq!(step.receipt.inner_boundary_outward_mass_flow_kg_per_s, 0.0);
                    assert_eq!(step.receipt.outer_boundary_outward_mass_flow_kg_per_s, 0.0);
                    assert_eq!(step.checkpoint.capillary_reservoir_mass_kg, None);
                    assert_eq!(step.checkpoint.capillary_reservoir_mass_kg_per_m, None);
                }
            }
        }
    }
}

#[test]
fn g1_radial_finite_reservoir_equilibrates_actual_total_volume_and_mass() {
    for inner in [0.0, 0.007] {
        for volume in [5.0e-9, 2.0e-8] {
            for initial_reservoir_pressure in [100_325.0, 102_325.0] {
                let mut input =
                    radial_connection_fixture(4, initial_reservoir_pressure, 128, Some(volume));
                // The physical outer boundary is at the final active circular
                // face. Excluded material contributes neither film volume nor load.
                input.grid.contact_exclusion.excluded[3] = true;
                if let GasFilmBoundaryTopology::AxisymmetricCapillaryReservoir {
                    cell_index, ..
                } = &mut input.boundary
                {
                    *cell_index = 2;
                }
                input.timestep_s = 0.01;
                input.budget.mass_residual_tolerance_kg_m2_s = 1.0e-13;
                let outer = inner + 0.75 * input.grid.length_m;
                let area = std::f64::consts::PI * (outer * outer - inner * inner);
                let film_volume = area * input.grid.gap_m[0];
                let rt = input.gas.specific_gas_constant_j_kg_k * input.gas.temperature_k;
                let initial_mass = (input.initial_absolute_pressure_pa * film_volume
                    + initial_reservoir_pressure * volume)
                    / rt;
                let equilibrium = initial_mass * rt / (film_volume + volume);
                let mut checkpoint = None;
                for index in 0..100 {
                    let step = solve_isothermal_axisymmetric_gas_film_1d(
                        &input,
                        inner,
                        checkpoint.as_ref(),
                    )
                    .unwrap();
                    let reservoir = step.receipt.finite_reservoir.unwrap();
                    assert!(
                        (reservoir.absolute_pressure_pa - reservoir.gas_mass_kg * rt / volume)
                            .abs()
                            < 1.0e-9
                    );
                    assert!(
                        (step.receipt.gas_mass_kg + reservoir.gas_mass_kg - initial_mass).abs()
                            < 2.0e-16
                    );
                    assert!(
                        reservoir.combined_mass_closure_residual_kg_per_s.abs() < 2.0e-16,
                        "inner={inner}, volume={volume}, step={index}, closure={:e}",
                        reservoir.combined_mass_closure_residual_kg_per_s
                    );
                    let flow = step.receipt.vent_outward_mass_flow_kg_per_s;
                    assert_eq!(
                        reservoir.inward_enthalpy_power_w,
                        step.receipt.vent_outward_enthalpy_power_w
                    );
                    assert_eq!(reservoir.required_isothermal_heat_to_gas_w, -rt * flow);
                    let stored_energy_rate = reservoir.storage_rate_kg_per_s
                        * (input.gas.declared_specific_enthalpy_j_kg - rt);
                    assert!(
                        (stored_energy_rate
                            - reservoir.inward_enthalpy_power_w
                            - reservoir.required_isothermal_heat_to_gas_w)
                            .abs()
                            < 1.0e-13
                    );
                    assert_eq!(
                        step.checkpoint.capillary_reservoir_mass_kg,
                        Some(reservoir.gas_mass_kg)
                    );
                    assert_eq!(step.checkpoint.capillary_reservoir_mass_kg_per_m, None);
                    assert_eq!(step.absolute_pressure_pa[3], None);
                    if index == 99 {
                        for pressure in &step.checkpoint.active_absolute_pressure_pa {
                            assert!((pressure - equilibrium).abs() < 0.01);
                        }
                        assert!((reservoir.absolute_pressure_pa - equilibrium).abs() < 0.01);
                    }
                    checkpoint = Some(step.checkpoint);
                }
            }
        }
    }
}

#[test]
fn g1_radial_reservoir_time_refinement_matches_independent_three_volume_ode() {
    for inner in [0.0, 0.01] {
        let mut input = radial_connection_fixture(2, 100_325.0, 128, Some(2.0e-9));
        let pi = std::f64::consts::PI;
        let h = input.grid.gap_m[0];
        let dr = input.grid.length_m / 2.0;
        let rt = input.gas.specific_gas_constant_j_kg_k * input.gas.temperature_k;
        let mu = input.gas.dynamic_viscosity_pa_s;
        let capacity: [f64; 3] = [
            pi * ((inner + dr).powi(2) - inner.powi(2)) * h / rt,
            pi * ((inner + 2.0 * dr).powi(2) - (inner + dr).powi(2)) * h / rt,
            2.0e-9 / rt,
        ];
        let radial =
            pi * h.powi(3) / (12.0 * mu * rt * ((inner + 1.5 * dr) / (inner + 0.5 * dr)).ln());
        let tube = 128.0 * pi * 1.0e-5_f64.powi(4) / (16.0 * mu * 0.01 * rt);
        let rate = |p: [f64; 3]| {
            let internal = radial * (p[0] * p[0] - p[1] * p[1]);
            let exchange = tube * (p[1] * p[1] - p[2] * p[2]);
            [
                -internal / capacity[0],
                (internal - exchange) / capacity[1],
                exchange / capacity[2],
            ]
        };
        let duration = 0.01;
        let dt = duration / 16_384.0;
        let mut reference = [
            input.initial_absolute_pressure_pa,
            input.initial_absolute_pressure_pa,
            100_325.0,
        ];
        for _ in 0..16_384 {
            let a = rate(reference);
            let b = rate(core::array::from_fn(|i| reference[i] + 0.5 * dt * a[i]));
            let c = rate(core::array::from_fn(|i| reference[i] + 0.5 * dt * b[i]));
            let d = rate(core::array::from_fn(|i| reference[i] + dt * c[i]));
            reference = core::array::from_fn(|i| {
                reference[i] + dt * (a[i] + 2.0 * b[i] + 2.0 * c[i] + d[i]) / 6.0
            });
        }
        let mut errors = Vec::new();
        for count in [16, 32, 64] {
            input.timestep_s = duration / f64::from(count);
            let mut checkpoint = None;
            let mut actual = [0.0; 3];
            for _ in 0..count {
                let step =
                    solve_isothermal_axisymmetric_gas_film_1d(&input, inner, checkpoint.as_ref())
                        .unwrap();
                actual = [
                    step.checkpoint.active_absolute_pressure_pa[0],
                    step.checkpoint.active_absolute_pressure_pa[1],
                    step.receipt.finite_reservoir.unwrap().absolute_pressure_pa,
                ];
                checkpoint = Some(step.checkpoint);
            }
            errors.push(
                actual
                    .into_iter()
                    .zip(reference)
                    .map(|(p, q)| (p - q).abs())
                    .fold(0.0_f64, f64::max),
            );
        }
        for pair in errors.windows(2) {
            assert!(
                pair[0] / pair[1] > 1.8 && pair[0] / pair[1] < 2.1,
                "radial reservoir time convergence: inner={inner}, errors={errors:?}"
            );
        }
    }
}

#[test]
fn g1_radial_reservoir_handles_stiff_capacity_and_fixed_bath_limit() {
    for inner in [0.0, 0.01] {
        for volume in [1.0e-18, 1.0e-24] {
            let mut input = radial_connection_fixture(2, 100_325.0, 64, Some(volume));
            input.timestep_s = 1.0e-3;
            let step = solve_isothermal_axisymmetric_gas_film_1d(&input, inner, None).unwrap();
            let reservoir = step.receipt.finite_reservoir.unwrap();
            assert!(
                (reservoir.absolute_pressure_pa - step.checkpoint.active_absolute_pressure_pa[1])
                    .abs()
                    < 1.0e-3
            );
            assert!(reservoir.gas_mass_kg > 0.0);
            assert!(step.receipt.vent_outward_mass_flow_kg_per_s > 0.0);
            assert!(step.receipt.capillary_vent.unwrap().viscous_dissipation_w >= 0.0);
            if volume < 1.0e-20 {
                assert!(reservoir.gas_mass_kg < 1.0e-20);
            }
        }
        let mut bath = radial_connection_fixture(2, 100_325.0, 64, None);
        bath.timestep_s = 1.0e-3;
        let reference = solve_isothermal_axisymmetric_gas_film_1d(&bath, inner, None).unwrap();
        let mut errors = Vec::new();
        for volume in [1.0e-7, 1.0e-6, 1.0e-5] {
            let mut input = radial_connection_fixture(2, 100_325.0, 64, Some(volume));
            input.timestep_s = bath.timestep_s;
            let step = solve_isothermal_axisymmetric_gas_film_1d(&input, inner, None).unwrap();
            let error = step.checkpoint.active_absolute_pressure_pa[1]
                - reference.checkpoint.active_absolute_pressure_pa[1];
            assert!(error > 0.0);
            errors.push(error);
        }
        for pair in errors.windows(2) {
            assert!(
                pair[0] / pair[1] > 9.0 && pair[0] / pair[1] < 11.0,
                "radial bath limit: {errors:?}"
            );
        }
    }
}

#[test]
fn g3_radial_reservoir_replays_squeeze_and_expansion_with_total_inventory() {
    let mut input = radial_connection_fixture(4, 101_325.0, 128, Some(3.0e-9));
    let inner = 0.002;
    input.timestep_s = 1.0e-3;
    input.wall_motion.gap_rate_m_per_s = -1.0e-5;
    let outer = inner + input.grid.length_m;
    let area = std::f64::consts::PI * (outer * outer - inner * inner);
    let rt = input.gas.specific_gas_constant_j_kg_k * input.gas.temperature_k;
    let old_gap = input.grid.gap_m[0] - input.timestep_s * input.wall_motion.gap_rate_m_per_s;
    let mass = input.initial_absolute_pressure_pa * (area * old_gap + 3.0e-9) / rt;
    let mut checkpoint = None;
    let mut filling = false;
    let mut draining = false;
    for index in 0..30 {
        let step =
            solve_isothermal_axisymmetric_gas_film_1d(&input, inner, checkpoint.as_ref()).unwrap();
        if index == 15 {
            assert_eq!(
                step,
                solve_isothermal_axisymmetric_gas_film_1d(&input, inner, checkpoint.as_ref())
                    .unwrap()
            );
        }
        let reservoir = step.receipt.finite_reservoir.unwrap();
        assert!((step.receipt.gas_mass_kg + reservoir.gas_mass_kg - mass).abs() < 1.0e-17);
        filling |= step.receipt.vent_outward_mass_flow_kg_per_s > 0.0;
        draining |= step.receipt.vent_outward_mass_flow_kg_per_s < 0.0;
        checkpoint = Some(step.checkpoint);
        if index == 9 {
            input.wall_motion.gap_rate_m_per_s = 2.0e-5;
        }
        for gap in &mut input.grid.gap_m {
            *gap += input.timestep_s * input.wall_motion.gap_rate_m_per_s;
        }
    }
    assert!(filling && draining);
}

#[test]
fn g4_radial_connections_reject_wrong_units_geometry_and_restart_inventory() {
    let input = radial_connection_fixture(2, 100_325.0, 32, Some(3.0e-9));
    let step = solve_isothermal_axisymmetric_gas_film_1d(&input, 0.0, None).unwrap();
    for total_mass in [None, Some(0.0), Some(-1.0), Some(f64::NAN)] {
        let mut invalid = step.checkpoint.clone();
        invalid.capillary_reservoir_mass_kg = total_mass;
        assert!(solve_isothermal_axisymmetric_gas_film_1d(&input, 0.0, Some(&invalid)).is_err());
    }
    let mut wrong_basis = step.checkpoint.clone();
    wrong_basis.capillary_reservoir_mass_kg_per_m = wrong_basis.capillary_reservoir_mass_kg;
    assert!(matches!(
        solve_isothermal_axisymmetric_gas_film_1d(&input, 0.0, Some(&wrong_basis)),
        Err(GasFilmError::CheckpointMismatch {
            field: "capillary_reservoir_mass_kg_per_m"
        })
    ));
    for change in 0..4 {
        let mut changed = input.clone();
        if let GasFilmBoundaryTopology::AxisymmetricCapillaryReservoir {
            cell_index,
            capillary,
            reservoir,
        } = &mut changed.boundary
        {
            match change {
                0 => capillary.capillary_count *= 2,
                1 => reservoir.volume_m3 *= 2.0,
                2 => reservoir.initial_absolute_pressure_pa += 1.0,
                _ => *cell_index = 0,
            }
        }
        assert!(matches!(
            solve_isothermal_axisymmetric_gas_film_1d(&changed, 0.0, Some(&step.checkpoint)),
            Err(GasFilmError::CheckpointMismatch {
                field: "configuration_fingerprint"
            })
        ));
    }
    for boundary in [
        input.boundary,
        radial_connection_fixture(2, 100_325.0, 32, None).boundary,
    ] {
        let mut changed = input.clone();
        changed.boundary = boundary;
        assert_eq!(
            solve_isothermal_gas_film_1d(&changed, None),
            Err(GasFilmError::Unavailable {
                reason: "cartesian-total-unit-capillary-connection"
            })
        );
    }
    let mut changed = input.clone();
    if let GasFilmBoundaryTopology::AxisymmetricCapillaryReservoir { capillary, .. } =
        &mut changed.boundary
    {
        capillary.capillary_count = 0;
    }
    assert!(matches!(
        solve_isothermal_axisymmetric_gas_film_1d(&changed, 0.0, None),
        Err(GasFilmError::InvalidInput {
            field: "capillary.capillary_count"
        })
    ));
    let mut equilibrium = radial_connection_fixture(2, 101_325.0, 32, Some(3.0e-9));
    equilibrium.budget.maximum_iterations = 1;
    let parked = solve_isothermal_axisymmetric_gas_film_1d(&equilibrium, 0.0, None).unwrap();
    let mut squeezed = equilibrium.clone();
    squeezed.wall_motion.gap_rate_m_per_s = -1.0e-4;
    for gap in &mut squeezed.grid.gap_m {
        *gap += squeezed.timestep_s * squeezed.wall_motion.gap_rate_m_per_s;
    }
    assert!(matches!(
        solve_isothermal_axisymmetric_gas_film_1d(&squeezed, 0.0, Some(&parked.checkpoint)),
        Err(GasFilmError::IterationBudgetExceeded { .. })
    ));
    let resumed =
        solve_isothermal_axisymmetric_gas_film_1d(&equilibrium, 0.0, Some(&parked.checkpoint))
            .unwrap();
    assert_eq!(
        resumed.checkpoint.capillary_reservoir_mass_kg,
        parked.checkpoint.capillary_reservoir_mass_kg
    );
}

#[test]
fn g1_stiff_radial_drainage_reaches_exact_endpoint_in_bounded_newton_updates() {
    let mut heat_errors = Vec::new();
    for cells in [32, 128, 256] {
        let mut input = fixture(cells);
        let inner: f64 = 0.005;
        let left: f64 = 120_000.0;
        let right: f64 = 100_000.0;
        let log_span = ((inner + input.grid.length_m) / inner).ln();
        let h = input.grid.gap_m[0];
        let mu = input.gas.dynamic_viscosity_pa_s;
        let rt = input.gas.specific_gas_constant_j_kg_k * input.gas.temperature_k;
        input.boundary = GasFilmBoundaryTopology::Open {
            left_absolute_pressure_pa: left,
            right_absolute_pressure_pa: right,
        };
        // The time step is much longer than lateral pressure relaxation, so
        // independent cell updates are diffusion-stiff. Keep the existing mass
        // tolerance and require a bounded number of full nonlinear corrections.
        input.timestep_s = 1.0;
        input.budget.maximum_iterations = 12;
        input.budget.relaxation = 1.0;
        let exact_pressure = (0..cells)
            .map(|cell| {
                let radius = inner + (cell as f64 + 0.5) * input.grid.length_m / cells as f64;
                (left.powi(2) + (right.powi(2) - left.powi(2)) * (radius / inner).ln() / log_span)
                    .sqrt()
            })
            .collect::<Vec<_>>();
        // Manufacture the endpoint through prescribed gap motion, not the
        // numerical flux or pressure initial guess. Initial pressure is uniform;
        // old_gap=h*p_exact/p_initial makes the endpoint storage exactly zero.
        // The independent continuum logarithmic profile then has constant flow.
        input.wall_motion.gap_rate_profile_m_per_s = Some(
            exact_pressure
                .iter()
                .map(|pressure| {
                    (h - h * pressure / input.initial_absolute_pressure_pa) / input.timestep_s
                })
                .collect(),
        );
        let step = solve_isothermal_axisymmetric_gas_film_1d(&input, inner, None)
            .expect("stiff radial pressure endpoint within the fixed Newton budget");
        let pressure_error = step
            .checkpoint
            .active_absolute_pressure_pa
            .iter()
            .zip(&exact_pressure)
            .map(|(actual, expected)| (actual / expected - 1.0).abs())
            .fold(0.0_f64, f64::max);
        let flow = std::f64::consts::PI * h.powi(3) * (left.powi(2) - right.powi(2))
            / (12.0 * mu * rt * log_span);
        let exact_heat = flow * rt * (left / right).ln();
        let heat_error = (step.receipt.viscous_heat_w / exact_heat - 1.0).abs();
        eprintln!(
            "stiff radial {cells} cells: {} Newton corrections, residual {:.3e}, pressure error {:.3e}, heat error {:.3e}",
            step.receipt.iterations,
            step.receipt.max_mass_residual_kg_m2_s,
            pressure_error,
            heat_error
        );
        assert!(step.receipt.iterations <= 8);
        assert!(pressure_error < 1.0e-9, "{pressure_error:e}");
        assert!(
            (step.receipt.outer_boundary_outward_mass_flow_kg_per_s / flow - 1.0).abs() < 1.0e-9
        );
        assert!(
            (-step.receipt.inner_boundary_outward_mass_flow_kg_per_s / flow - 1.0).abs() < 1.0e-9
        );
        assert!(step.receipt.mass_closure_residual_kg_per_s.abs() < 1.0e-14);
        assert!(
            step.receipt.max_mass_residual_kg_m2_s <= input.budget.mass_residual_tolerance_kg_m2_s
        );
        heat_errors.push(heat_error);
    }
    assert!(
        heat_errors
            .windows(2)
            .all(|errors| errors[1] < errors[0] / 3.0),
        "{heat_errors:?}"
    );
    assert!(heat_errors[2] < 1.0e-6, "{heat_errors:?}");
}

#[test]
fn g1_axisymmetric_sealed_squeeze_integrates_actual_annular_mass_load_and_work() {
    let mut disk_masses = Vec::new();
    for inner in [0.0, 0.005] {
        for span in [0.01, 0.02] {
            let mut input = fixture(8);
            input.grid.length_m = span;
            input.grid.contact_exclusion.excluded[6..].fill(true);
            input.grid.gap_m[6..].fill(0.0);
            input.wall_motion.gap_rate_m_per_s = -1.0e-4;
            input.timestep_s = 1.0e-4;
            input.budget.mass_residual_tolerance_kg_m2_s = 1.0e-12;
            let step = solve_isothermal_axisymmetric_gas_film_1d(&input, inner, None).unwrap();
            assert_eq!(
                step,
                solve_isothermal_axisymmetric_gas_film_1d(&input, inner, None).unwrap()
            );
            let outer = inner + span * 6.0 / 8.0;
            let area = std::f64::consts::PI * (outer * outer - inner * inner);
            let h = input.grid.gap_m[0];
            let old_h = h - input.wall_motion.gap_rate_m_per_s * input.timestep_s;
            let p = input.initial_absolute_pressure_pa * old_h / h;
            let rt = input.gas.specific_gas_constant_j_kg_k * input.gas.temperature_k;
            let expected_mass = input.initial_absolute_pressure_pa * old_h * area / rt;
            let receipt = &step.receipt;
            assert!((receipt.gas_mass_kg / expected_mass - 1.0).abs() < 1.0e-9);
            assert!((receipt.absolute_pressure_load_n / (p * area) - 1.0).abs() < 1.0e-9);
            assert!(
                (receipt.gauge_pressure_load_n / ((p - input.initial_absolute_pressure_pa) * area)
                    - 1.0)
                    .abs()
                    < 1.0e-6
            );
            assert!(
                (receipt.normal_gap_power_to_gas_w
                    / (-p * input.wall_motion.gap_rate_m_per_s * area)
                    - 1.0)
                    .abs()
                    < 1.0e-9
            );
            assert!(receipt.mass_closure_residual_kg_per_s.abs() < 1.0e-14);
            assert_eq!(receipt.inner_boundary_outward_mass_flow_kg_per_s, 0.0);
            assert_eq!(receipt.outer_boundary_outward_mass_flow_kg_per_s, 0.0);
            assert!(step.absolute_pressure_pa[6..].iter().all(Option::is_none));
            assert!(step.checkpoint.capillary_reservoir_mass_kg_per_m.is_none());
            if inner == 0.0 {
                disk_masses.push(receipt.gas_mass_kg);
            }
        }
    }
    assert!((disk_masses[1] / disk_masses[0] - 4.0).abs() < 1.0e-9);
}

#[test]
fn g1_axisymmetric_annular_log_pressure_flux_and_heat_with_maxwell_slip() {
    let inner: f64 = 0.005;
    for slip in [false, true] {
        for (left, right) in [
            (120_000.0_f64, 100_000.0_f64),
            (100_000.0_f64, 120_000.0_f64),
        ] {
            let mut heat_errors = Vec::new();
            for cells in [8, 32] {
                let mut input = fixture(cells);
                if slip {
                    input.slip_policy = SlipPolicy::MaxwellFirstOrder {
                        source_id: "synthetic-radial-maxwell-v1".into(),
                        tangential_momentum_accommodation: 0.8,
                    };
                }
                input.boundary = GasFilmBoundaryTopology::Open {
                    left_absolute_pressure_pa: right,
                    right_absolute_pressure_pa: right,
                };
                let mut state = solve_isothermal_axisymmetric_gas_film_1d(&input, inner, None)
                    .unwrap()
                    .checkpoint;
                input.boundary = GasFilmBoundaryTopology::Open {
                    left_absolute_pressure_pa: left,
                    right_absolute_pressure_pa: right,
                };
                let h = input.grid.gap_m[0];
                let mu = input.gas.dynamic_viscosity_pa_s;
                let rt = input.gas.specific_gas_constant_j_kg_k * input.gas.temperature_k;
                let shift = if slip {
                    6.0 * 1.5
                        * input.applicability.mean_free_path_m
                        * input.initial_absolute_pressure_pa
                        / h
                } else {
                    0.0
                };
                let log_span = ((inner + input.grid.length_m) / inner).ln();
                // Exact steady continuum potential: (p+6*b*p/h)^2 is affine
                // in ln(r), including pressure-dependent Maxwell slip.
                state.active_absolute_pressure_pa = (0..cells)
                    .map(|cell| {
                        let radius =
                            inner + (cell as f64 + 0.5) * input.grid.length_m / cells as f64;
                        let x = (radius / inner).ln() / log_span;
                        ((left + shift).powi(2)
                            + ((right + shift).powi(2) - (left + shift).powi(2)) * x)
                            .sqrt()
                            - shift
                    })
                    .collect();
                let expected_pressure = state.active_absolute_pressure_pa.clone();
                let step =
                    solve_isothermal_axisymmetric_gas_film_1d(&input, inner, Some(&state)).unwrap();
                let flow = std::f64::consts::PI
                    * h.powi(3)
                    * (left - right)
                    * (left + right + 2.0 * shift)
                    / (12.0 * mu * rt * log_span);
                let receipt = &step.receipt;
                assert!(
                    (receipt.outer_boundary_outward_mass_flow_kg_per_s / flow - 1.0).abs()
                        < 1.0e-10
                );
                assert!(
                    (-receipt.inner_boundary_outward_mass_flow_kg_per_s / flow - 1.0).abs()
                        < 1.0e-10
                );
                for (actual, expected) in step
                    .checkpoint
                    .active_absolute_pressure_pa
                    .iter()
                    .zip(expected_pressure)
                {
                    assert!((actual / expected - 1.0).abs() < 1.0e-12);
                }
                assert!(
                    receipt
                        .gas_on_wall_radial_traction_pa
                        .iter()
                        .all(|value| value.unwrap() * (left - right) > 0.0)
                );
                assert!(
                    (receipt.pressure_gradient_power_to_gas_w / receipt.viscous_heat_w - 1.0).abs()
                        < 1.0e-12
                );
                assert_eq!(receipt.normal_gap_power_to_gas_w, 0.0);
                assert!(receipt.mass_closure_residual_kg_per_s.abs() < 1.0e-14);
                assert_eq!(
                    receipt.outer_boundary_outward_enthalpy_power_w,
                    receipt.outer_boundary_outward_mass_flow_kg_per_s
                        * input.gas.declared_specific_enthalpy_j_kg
                );
                if slip {
                    assert!(receipt.wall_slip_dissipation_w > 0.0);
                    assert!(receipt.wall_slip_dissipation_w < receipt.viscous_heat_w);
                } else {
                    assert_eq!(receipt.wall_slip_dissipation_w, 0.0);
                }
                let exact_heat = flow * rt * (left / right).ln();
                heat_errors.push((receipt.viscous_heat_w / exact_heat - 1.0).abs());
            }
            assert!(
                heat_errors[1] < 1.0e-4 && heat_errors[1] < heat_errors[0] / 8.0,
                "{heat_errors:?}"
            );
        }
    }
}

#[test]
fn g1_open_rim_disk_squeeze_and_expansion_refine_to_compressible_reynolds_solution() {
    for central_gauge_pressure in [1_000.0, -1_000.0] {
        let mut pressure_errors = Vec::new();
        let mut load_errors = Vec::new();
        for cells in [8, 16, 32] {
            let mut input = fixture(cells);
            input.timestep_s = 2.0e-3;
            input.budget.maximum_iterations = 30_000;
            input.budget.mass_residual_tolerance_kg_m2_s = 1.0e-11;
            let ambient = input.initial_absolute_pressure_pa;
            input.boundary = GasFilmBoundaryTopology::MixedEnds {
                left_absolute_pressure_pa: None,
                right_absolute_pressure_pa: Some(ambient),
            };
            let mut old = solve_isothermal_axisymmetric_gas_film_1d(&input, 0.0, None)
                .unwrap()
                .checkpoint;
            let radius = input.grid.length_m;
            let h = input.grid.gap_m[0];
            let mu = input.gas.dynamic_viscosity_pa_s;
            let c = central_gauge_pressure / radius.powi(2);
            let mut rates = Vec::new();
            // Manufacture a backward-Euler endpoint from the CONTINUUM radial
            // equation, without using any numerical face flux. For p=pa+C(R²-r²),
            // (1/r)d(r*p*p_r)/dr = -4*C*p+4*C²*r². Holding p across the
            // step therefore requires h_dot=-h³*C/(3*mu)*(1-C*r²/p).
            for cell in 0..cells {
                let r = (cell as f64 + 0.5) * radius / cells as f64;
                let p = ambient + c * (radius.powi(2) - r.powi(2));
                let rate = -h.powi(3) * c / (3.0 * mu) * (1.0 - c * r.powi(2) / p);
                old.active_absolute_pressure_pa[cell] = p;
                old.active_gap_m[cell] = h - rate * input.timestep_s;
                rates.push(rate);
            }
            input.wall_motion.gap_rate_profile_m_per_s = Some(rates);
            let step = solve_isothermal_axisymmetric_gas_film_1d(&input, 0.0, Some(&old)).unwrap();
            let max_error = step
                .checkpoint
                .active_absolute_pressure_pa
                .iter()
                .zip(&old.active_absolute_pressure_pa)
                .map(|(actual, expected)| (actual - expected).abs())
                .fold(0.0_f64, f64::max);
            pressure_errors.push(max_error / central_gauge_pressure.abs());
            let exact_load = 0.5 * std::f64::consts::PI * central_gauge_pressure * radius.powi(2);
            load_errors.push((step.receipt.gauge_pressure_load_n / exact_load - 1.0).abs());
            let receipt = &step.receipt;
            assert_eq!(receipt.inner_boundary_outward_mass_flow_kg_per_s, 0.0);
            assert!(
                receipt.outer_boundary_outward_mass_flow_kg_per_s * central_gauge_pressure > 0.0
            );
            assert!(receipt.normal_gap_power_to_gas_w * central_gauge_pressure > 0.0);
            assert!(receipt.mass_closure_residual_kg_per_s.abs() < 1.0e-14);
            assert!(receipt.viscous_heat_w > 0.0);
            assert!(
                (receipt.pressure_gradient_power_to_gas_w / receipt.viscous_heat_w - 1.0).abs()
                    < 1.0e-12
            );
        }
        assert!(pressure_errors[2] < 0.005, "{pressure_errors:?}");
        assert!(
            pressure_errors
                .windows(2)
                .all(|pair| pair[1] < pair[0] / 2.5),
            "{pressure_errors:?}"
        );
        assert!(load_errors[2] < 0.002, "{load_errors:?}");
        assert!(
            load_errors.windows(2).all(|pair| pair[1] < pair[0] / 3.0),
            "{load_errors:?}"
        );
    }
}

#[test]
fn g0_axisymmetric_restart_binds_geometry_and_refuses_unsupported_physics() {
    let mut input = fixture(8);
    input.boundary = GasFilmBoundaryTopology::MixedEnds {
        left_absolute_pressure_pa: None,
        right_absolute_pressure_pa: Some(input.initial_absolute_pressure_pa),
    };
    let first = solve_isothermal_axisymmetric_gas_film_1d(&input, 0.0, None).unwrap();
    for inner in [0.001, -1.0, f64::NAN] {
        assert!(
            solve_isothermal_axisymmetric_gas_film_1d(&input, inner, Some(&first.checkpoint))
                .is_err()
        );
    }
    assert!(matches!(
        solve_isothermal_gas_film_1d(&input, Some(&first.checkpoint)),
        Err(GasFilmError::CheckpointMismatch {
            field: "configuration_fingerprint"
        })
    ));
    let cartesian = solve_isothermal_gas_film_1d(&input, None).unwrap();
    assert!(
        solve_isothermal_axisymmetric_gas_film_1d(&input, 0.0, Some(&cartesian.checkpoint))
            .is_err()
    );
    input.wall_motion.gap_rate_profile_m_per_s = Some(
        (0..8)
            .map(|cell| -1.0e-5 * (1.0 - (cell as f64 + 0.5).powi(2) / 64.0))
            .collect(),
    );
    for (gap, rate) in input
        .grid
        .gap_m
        .iter_mut()
        .zip(input.wall_motion.gap_rate_profile_m_per_s.as_ref().unwrap())
    {
        *gap += rate * input.timestep_s;
    }
    let evolved =
        solve_isothermal_axisymmetric_gas_film_1d(&input, 0.0, Some(&first.checkpoint)).unwrap();
    assert_eq!(
        evolved,
        solve_isothermal_axisymmetric_gas_film_1d(&input, 0.0, Some(&first.checkpoint)).unwrap()
    );
    assert_eq!(evolved.checkpoint.active_gap_m, input.grid.gap_m);
    input.boundary = GasFilmBoundaryTopology::Open {
        left_absolute_pressure_pa: input.initial_absolute_pressure_pa,
        right_absolute_pressure_pa: input.initial_absolute_pressure_pa,
    };
    assert_eq!(
        solve_isothermal_axisymmetric_gas_film_1d(&input, 0.0, None),
        Err(GasFilmError::Unavailable {
            reason: "axisymmetric-axis-must-be-sealed"
        })
    );
    input.boundary = GasFilmBoundaryTopology::Sealed;
    input.wall_motion.upper_tangential_velocity_m_per_s = 0.1;
    assert_eq!(
        solve_isothermal_axisymmetric_gas_film_1d(&input, 0.0, None),
        Err(GasFilmError::Unavailable {
            reason: "axisymmetric-in-plane-wall-motion"
        })
    );
    assert_eq!(
        solve_isothermal_axisymmetric_gas_film_1d(&capillary_fixture(8, 101_325.0, 1.0), 0.0, None),
        Err(GasFilmError::Unavailable {
            reason: "axisymmetric-per-width-capillary-connection"
        })
    );
    assert_eq!(
        solve_isothermal_axisymmetric_gas_film_1d(
            &reservoir_fixture(8, 101_325.0, 1.0e-7),
            0.0,
            None
        ),
        Err(GasFilmError::Unavailable {
            reason: "axisymmetric-per-width-capillary-connection"
        })
    );
    input = fixture(8);
    input.boundary = GasFilmBoundaryTopology::MixedEnds {
        left_absolute_pressure_pa: None,
        right_absolute_pressure_pa: Some(1_000.0),
    };
    assert_eq!(
        solve_isothermal_axisymmetric_gas_film_1d(&input, 0.0, None),
        Err(GasFilmError::Unavailable {
            reason: "rarefied-knudsen-outside-continuum-envelope"
        })
    );
}

#[test]
fn g1_mixed_cartesian_ends_retain_wall_backpressure_and_mechanical_balance() {
    for left_open in [false, true] {
        let mut input = maxwell_fixture(8, 0.8);
        input.boundary = GasFilmBoundaryTopology::MixedEnds {
            left_absolute_pressure_pa: left_open.then_some(input.initial_absolute_pressure_pa),
            right_absolute_pressure_pa: (!left_open).then_some(input.initial_absolute_pressure_pa),
        };
        input.wall_motion.lower_tangential_velocity_m_per_s = 0.1;
        input.wall_motion.upper_tangential_velocity_m_per_s = 0.1;
        let step = solve_isothermal_gas_film_1d(&input, None).unwrap();
        if left_open {
            assert_eq!(
                step.receipt.right_boundary_outward_mass_flux_kg_per_m_s,
                0.0
            );
            assert!(step.receipt.left_boundary_outward_mass_flux_kg_per_m_s < 0.0);
        } else {
            assert_eq!(step.receipt.left_boundary_outward_mass_flux_kg_per_m_s, 0.0);
            assert!(step.receipt.right_boundary_outward_mass_flux_kg_per_m_s > 0.0);
        }
        assert!(step.receipt.mass_closure_residual_kg_per_m_s.abs() < 1.0e-11);
        assert_mechanical_balance(&input, &step);
    }
}

fn capillary_fixture(cells: usize, bath_pressure: f64, count_per_m: f64) -> GasFilmInput {
    let mut input = fixture(cells);
    input.boundary = GasFilmBoundaryTopology::CapillaryVented {
        cell_index: 0,
        absolute_pressure_pa: bath_pressure,
        capillary: CapillaryVent {
            radius_m: 1.0e-5,
            length_m: 1.0e-2,
            capillaries_per_m: count_per_m,
        },
    };
    input.budget.mass_residual_tolerance_kg_m2_s = 1.0e-12;
    input
}

fn reservoir_fixture(cells: usize, pressure: f64, volume: f64) -> GasFilmInput {
    let mut input = capillary_fixture(cells, pressure, 10_000.0);
    let GasFilmBoundaryTopology::CapillaryVented {
        cell_index,
        capillary,
        ..
    } = input.boundary
    else {
        unreachable!()
    };
    input.boundary = GasFilmBoundaryTopology::CapillaryReservoir {
        cell_index,
        capillary,
        reservoir: GasFilmReservoir {
            volume_per_width_m2: volume,
            initial_absolute_pressure_pa: pressure,
        },
    };
    input
}

#[test]
fn g1_finite_reservoir_backpressure_equilibrates_with_conservative_gas_storage() {
    for initial_reservoir_pressure in [100_325.0, 102_325.0] {
        let mut input = reservoir_fixture(2, initial_reservoir_pressure, 1.0e-7);
        input.timestep_s = 2.5e-3;
        let rt = input.gas.specific_gas_constant_j_kg_k * input.gas.temperature_k;
        let film_volume = input.grid.gap_m[0] * input.grid.length_m;
        let reservoir_volume = 1.0e-7;
        let initial_mass = (input.initial_absolute_pressure_pa * film_volume
            + initial_reservoir_pressure * reservoir_volume)
            / rt;
        let equilibrium = initial_mass * rt / (film_volume + reservoir_volume);
        let mut checkpoint = None;
        let mut last = None;
        for _ in 0..100 {
            let step = solve_isothermal_gas_film_1d(&input, checkpoint.as_ref()).unwrap();
            let reservoir = step.receipt.finite_reservoir.unwrap();
            assert!(reservoir.gas_mass_kg_per_m > 0.0);
            assert!(
                (reservoir.absolute_pressure_pa
                    - reservoir.gas_mass_kg_per_m * rt / reservoir_volume)
                    .abs()
                    < 1.0e-9
            );
            assert!(
                (step.receipt.gas_mass_kg_per_m + reservoir.gas_mass_kg_per_m - initial_mass).abs()
                    < 3.0e-14
            );
            assert!(reservoir.combined_mass_closure_residual_kg_per_m_s.abs() < 2.0e-14);
            assert_eq!(
                reservoir.inward_enthalpy_flux_w_per_m,
                step.receipt.vent_outward_enthalpy_flux_w_per_m
            );
            let flow = step.receipt.vent_outward_mass_flux_kg_per_m_s;
            assert_eq!(
                reservoir.required_isothermal_heat_to_gas_w_per_m,
                -rt * flow
            );
            // U=m*(h-RT) at this fixed ideal-gas temperature. Reconstruct the
            // reservoir's first law independently from its represented storage.
            let energy_storage = reservoir.storage_rate_kg_per_m_s
                * (input.gas.declared_specific_enthalpy_j_kg - rt);
            let energy_input = reservoir.inward_enthalpy_flux_w_per_m
                + reservoir.required_isothermal_heat_to_gas_w_per_m;
            assert!((energy_storage - energy_input).abs() < 5.0e-9);
            if last.is_none() {
                let direction = input.initial_absolute_pressure_pa - initial_reservoir_pressure;
                assert!(flow * direction > 0.0);
                assert!(
                    (reservoir.absolute_pressure_pa - initial_reservoir_pressure) * direction > 0.0
                );
                assert!(reservoir.required_isothermal_heat_to_gas_w_per_m * flow < 0.0);
            }
            checkpoint = Some(step.checkpoint.clone());
            last = Some(step);
        }
        let last = last.unwrap();
        for pressure in &last.checkpoint.active_absolute_pressure_pa {
            assert!((pressure - equilibrium).abs() < 0.01);
        }
        assert!(
            (last.receipt.finite_reservoir.unwrap().absolute_pressure_pa - equilibrium).abs()
                < 0.01
        );
    }
}

#[test]
fn g1_finite_reservoir_pressure_time_refinement_matches_three_storage_ode() {
    let mut input = reservoir_fixture(2, 100_325.0, 1.0e-7);
    let rt = input.gas.specific_gas_constant_j_kg_k * input.gas.temperature_k;
    let dx = input.grid.length_m / 2.0;
    let h = input.grid.gap_m[0];
    let mu = input.gas.dynamic_viscosity_pa_s;
    let film_storage = h * dx / rt;
    let reservoir_storage = 1.0e-7 / rt;
    let film = h.powi(3) / (24.0 * mu * rt * dx);
    let capillary = 10_000.0 * std::f64::consts::PI * 1.0e-5_f64.powi(4) / (16.0 * mu * 0.01 * rt);
    let rate = |p: [f64; 3]| {
        let internal = film * (p[0] - p[1]) * (p[0] + p[1]);
        let exchange = capillary * (p[0] - p[2]) * (p[0] + p[2]);
        [
            -(internal + exchange) / film_storage,
            internal / film_storage,
            exchange / reservoir_storage,
        ]
    };
    let duration = 0.005;
    let dt = duration / 8_192.0;
    let mut reference = [
        input.initial_absolute_pressure_pa,
        input.initial_absolute_pressure_pa,
        100_325.0,
    ];
    for _ in 0..8_192 {
        let a = rate(reference);
        let b = rate(core::array::from_fn(|i| reference[i] + 0.5 * dt * a[i]));
        let c = rate(core::array::from_fn(|i| reference[i] + 0.5 * dt * b[i]));
        let d = rate(core::array::from_fn(|i| reference[i] + dt * c[i]));
        reference = core::array::from_fn(|i| {
            reference[i] + dt * (a[i] + 2.0 * b[i] + 2.0 * c[i] + d[i]) / 6.0
        });
    }
    let mut errors = Vec::new();
    for count in [16, 32, 64] {
        input.timestep_s = duration / f64::from(count);
        let mut checkpoint = None;
        let mut last = None;
        for _ in 0..count {
            let step = solve_isothermal_gas_film_1d(&input, checkpoint.as_ref()).unwrap();
            checkpoint = Some(step.checkpoint.clone());
            last = Some(step);
        }
        let last = last.unwrap();
        let actual = [
            last.checkpoint.active_absolute_pressure_pa[0],
            last.checkpoint.active_absolute_pressure_pa[1],
            last.receipt.finite_reservoir.unwrap().absolute_pressure_pa,
        ];
        errors.push(
            actual
                .into_iter()
                .zip(reference)
                .map(|(p, r)| (p - r).abs())
                .fold(0.0_f64, f64::max),
        );
    }
    for pair in errors.windows(2) {
        assert!(
            pair[0] / pair[1] > 1.8 && pair[0] / pair[1] < 2.1,
            "finite reservoir time convergence: {errors:?}"
        );
    }
}

#[test]
fn g1_finite_reservoir_handles_stiff_small_volume_and_approaches_fixed_bath_limit() {
    for volume in [1.0e-14, 1.0e-24] {
        let mut input = reservoir_fixture(2, 100_325.0, volume);
        input.timestep_s = 1.0e-3;
        let step = solve_isothermal_gas_film_1d(&input, None).unwrap();
        let reservoir = step.receipt.finite_reservoir.unwrap();
        let film_pressure = step.checkpoint.active_absolute_pressure_pa[0];
        assert!((reservoir.absolute_pressure_pa - film_pressure).abs() < 1.0e-3);
        assert!(reservoir.gas_mass_kg_per_m > 0.0);
        assert!(step.receipt.vent_outward_mass_flux_kg_per_m_s > 0.0);
        assert!(
            step.receipt
                .capillary_vent
                .unwrap()
                .viscous_dissipation_w_per_m
                >= 0.0
        );
        if volume < 1.0e-20 {
            assert!(
                reservoir.gas_mass_kg_per_m < 1.0e-20,
                "no artificial mass/pressure floor"
            );
        }
    }
    let mut infinite = capillary_fixture(2, 100_325.0, 10_000.0);
    infinite.timestep_s = 1.0e-3;
    let reference = solve_isothermal_gas_film_1d(&infinite, None).unwrap();
    let mut errors = Vec::new();
    for volume in [1.0e-5, 1.0e-4, 1.0e-3] {
        let mut input = reservoir_fixture(2, 100_325.0, volume);
        input.timestep_s = infinite.timestep_s;
        let step = solve_isothermal_gas_film_1d(&input, None).unwrap();
        let finite = step.checkpoint.active_absolute_pressure_pa[0];
        let bath = reference.checkpoint.active_absolute_pressure_pa[0];
        assert!(finite > bath, "finite backpressure reduces discharge");
        errors.push(finite - bath);
    }
    for pair in errors.windows(2) {
        assert!(
            pair[0] / pair[1] > 9.0 && pair[0] / pair[1] < 11.0,
            "fixed-bath volume limit: {errors:?}"
        );
    }
}

#[test]
fn g3_finite_reservoir_tracks_squeeze_expansion_and_replays_without_resetting_inventory() {
    let mut input = reservoir_fixture(4, 101_325.0, 1.0e-7);
    input.timestep_s = 1.0e-3;
    input.wall_motion.gap_rate_m_per_s = -1.0e-5;
    let rt = input.gas.specific_gas_constant_j_kg_k * input.gas.temperature_k;
    let old_gap = input.grid.gap_m[0] - input.wall_motion.gap_rate_m_per_s * input.timestep_s;
    let initial_mass =
        input.initial_absolute_pressure_pa * (old_gap * input.grid.length_m + 1.0e-7) / rt;
    let mut checkpoint = None;
    let mut filled = false;
    let mut drained = false;
    for index in 0..30 {
        let step = solve_isothermal_gas_film_1d(&input, checkpoint.as_ref()).unwrap();
        if index == 15 {
            assert_eq!(
                step,
                solve_isothermal_gas_film_1d(&input, checkpoint.as_ref()).unwrap()
            );
        }
        let reservoir = step.receipt.finite_reservoir.unwrap();
        assert!(
            (step.receipt.gas_mass_kg_per_m + reservoir.gas_mass_kg_per_m - initial_mass).abs()
                < 3.0e-14
        );
        filled |= step.receipt.vent_outward_mass_flux_kg_per_m_s > 0.0;
        drained |= step.receipt.vent_outward_mass_flux_kg_per_m_s < 0.0;
        checkpoint = Some(step.checkpoint);
        if index == 9 {
            input.wall_motion.gap_rate_m_per_s = 2.0e-5;
        }
        for gap in &mut input.grid.gap_m {
            *gap += input.timestep_s * input.wall_motion.gap_rate_m_per_s;
        }
    }
    assert!(filled && drained);
}

#[test]
fn g4_finite_reservoir_refuses_invalid_or_missing_storage_and_preserves_restart() {
    let input = reservoir_fixture(2, 100_325.0, 1.0e-7);
    let first = solve_isothermal_gas_film_1d(&input, None).unwrap();
    let saved = first.checkpoint.clone();
    for value in [None, Some(0.0), Some(-1.0), Some(f64::NAN)] {
        let mut invalid = saved.clone();
        invalid.capillary_reservoir_mass_kg_per_m = value;
        assert!(solve_isothermal_gas_film_1d(&input, Some(&invalid)).is_err());
    }
    let mut changed = input.clone();
    if let GasFilmBoundaryTopology::CapillaryReservoir { reservoir, .. } = &mut changed.boundary {
        reservoir.volume_per_width_m2 *= 2.0;
    }
    assert!(matches!(
        solve_isothermal_gas_film_1d(&changed, Some(&saved)),
        Err(GasFilmError::CheckpointMismatch {
            field: "configuration_fingerprint"
        })
    ));
    if let GasFilmBoundaryTopology::CapillaryReservoir { reservoir, .. } = &mut changed.boundary {
        reservoir.volume_per_width_m2 = 0.0;
    }
    assert!(matches!(
        solve_isothermal_gas_film_1d(&changed, None),
        Err(GasFilmError::InvalidInput { .. })
    ));
    let mut equilibrium = reservoir_fixture(2, 101_325.0, 1.0e-7);
    equilibrium.budget.maximum_iterations = 1;
    let parked = solve_isothermal_gas_film_1d(&equilibrium, None).unwrap();
    let mut squeezed = equilibrium.clone();
    squeezed.wall_motion.gap_rate_m_per_s = -1.0e-5;
    for gap in &mut squeezed.grid.gap_m {
        *gap += squeezed.timestep_s * squeezed.wall_motion.gap_rate_m_per_s;
    }
    assert!(matches!(
        solve_isothermal_gas_film_1d(&squeezed, Some(&parked.checkpoint)),
        Err(GasFilmError::IterationBudgetExceeded { .. })
    ));
    assert_eq!(saved, first.checkpoint);
    let resumed = solve_isothermal_gas_film_1d(&equilibrium, Some(&parked.checkpoint)).unwrap();
    assert_eq!(
        resumed.checkpoint.capillary_reservoir_mass_kg_per_m,
        parked.checkpoint.capillary_reservoir_mass_kg_per_m
    );
}

#[test]
fn g1_capillary_vent_retains_pressure_and_closes_signed_mass_work_and_enthalpy() {
    for bath in [100_325.0, 102_325.0] {
        let mut coefficients = Vec::new();
        for radius in [1.0e-5, 2.0e-5] {
            let mut input = capillary_fixture(4, bath, 5_000.0);
            input.timestep_s = 5.0e-4;
            let GasFilmBoundaryTopology::CapillaryVented { capillary, .. } = &mut input.boundary
            else {
                unreachable!()
            };
            capillary.radius_m = radius;
            let tube = *capillary;
            let step = solve_isothermal_gas_film_1d(&input, None).unwrap();
            let p = step.absolute_pressure_pa[0].unwrap();
            assert!(p > bath.min(input.initial_absolute_pressure_pa));
            assert!(p < bath.max(input.initial_absolute_pressure_pa));
            assert!(
                (p - bath).abs() > 100.0,
                "finite vent must not pin its cell"
            );
            let rt = input.gas.specific_gas_constant_j_kg_k * input.gas.temperature_k;
            let coefficient = tube.capillaries_per_m * std::f64::consts::PI * radius.powi(4)
                / (16.0 * input.gas.dynamic_viscosity_pa_s * tube.length_m * rt);
            let expected = coefficient * (p * p - bath * bath);
            let flow = step.receipt.vent_outward_mass_flux_kg_per_m_s;
            assert!((flow - expected).abs() < expected.abs() * 1.0e-12);
            coefficients.push(flow / ((p - bath) * (p + bath)));
            assert!(flow * (input.initial_absolute_pressure_pa - bath) > 0.0);
            assert!(step.receipt.mass_closure_residual_kg_per_m_s.abs() < 2.0e-14);
            assert_eq!(step.receipt.left_boundary_outward_mass_flux_kg_per_m_s, 0.0);
            assert_eq!(
                step.receipt.right_boundary_outward_mass_flux_kg_per_m_s,
                0.0
            );
            assert_eq!(
                step.receipt.vent_outward_enthalpy_flux_w_per_m,
                flow * input.gas.declared_specific_enthalpy_j_kg
            );
            let report = step.receipt.capillary_vent.unwrap();
            let pressure_work = flow * rt * (p / bath).ln();
            assert!(report.viscous_dissipation_w_per_m > 0.0);
            assert!(
                (report.viscous_dissipation_w_per_m - pressure_work).abs()
                    < pressure_work * 2.0e-12
            );
            let peak = 2.0 * flow.abs() * rt
                / (tube.capillaries_per_m * std::f64::consts::PI * radius * radius * p.min(bath));
            assert!(
                (report.maximum_gas_mach_number
                    - peak / input.applicability.speed_of_sound_m_per_s)
                    .abs()
                    < 1.0e-16
            );
            // Capillary work is outside the planar film: adding it to either
            // side of this identity would count the connection twice.
            assert!(
                (step.receipt.viscous_heat_w_per_m
                    - step.receipt.pressure_gradient_power_to_gas_w_per_m)
                    .abs()
                    < 1.0e-15
            );
        }
        assert!((coefficients[1] / coefficients[0] - 16.0).abs() < 1.0e-12);
    }
}

#[test]
fn g1_capillary_pressure_relaxation_converges_to_independent_conservative_ode() {
    let mut input = capillary_fixture(2, 100_325.0, 10_000.0);
    let rt = input.gas.specific_gas_constant_j_kg_k * input.gas.temperature_k;
    let dx = input.grid.length_m / 2.0;
    let h = input.grid.gap_m[0];
    let mu = input.gas.dynamic_viscosity_pa_s;
    let storage = h * dx / rt;
    let film = h.powi(3) / (24.0 * mu * rt * dx);
    let vent = 10_000.0 * std::f64::consts::PI * 1.0e-5_f64.powi(4) / (16.0 * mu * 1.0e-2 * rt);
    let bath = 100_325.0;
    let rate = |p: [f64; 2]| {
        let transfer = film * (p[0] - p[1]) * (p[0] + p[1]);
        let discharge = vent * (p[0] - bath) * (p[0] + bath);
        [-(transfer + discharge) / storage, transfer / storage]
    };
    // Independently integrate the two physical storage balances with RK4.
    // Production uses backward Euler and a coupled nonlinear pressure solve.
    let duration = 5.0e-3;
    let dt = duration / 8_192.0;
    let mut reference = [input.initial_absolute_pressure_pa; 2];
    for _ in 0..8_192 {
        let a = rate(reference);
        let b = rate(core::array::from_fn(|i| reference[i] + 0.5 * dt * a[i]));
        let c = rate(core::array::from_fn(|i| reference[i] + 0.5 * dt * b[i]));
        let d = rate(core::array::from_fn(|i| reference[i] + dt * c[i]));
        reference = core::array::from_fn(|i| {
            reference[i] + dt * (a[i] + 2.0 * b[i] + 2.0 * c[i] + d[i]) / 6.0
        });
    }
    let mut errors = Vec::new();
    for count in [16, 32, 64] {
        input.timestep_s = duration / f64::from(count);
        let mut checkpoint = None;
        let mut exchanged = 0.0;
        let mut last = None;
        for index in 0..count {
            let step = solve_isothermal_gas_film_1d(&input, checkpoint.as_ref()).unwrap();
            if index == count / 2 {
                assert_eq!(
                    step,
                    solve_isothermal_gas_film_1d(&input, checkpoint.as_ref()).unwrap()
                );
            }
            exchanged += input.timestep_s * step.receipt.vent_outward_mass_flux_kg_per_m_s;
            checkpoint = Some(step.checkpoint.clone());
            last = Some(step);
        }
        let last = last.unwrap();
        let initial_mass = input.initial_absolute_pressure_pa * h * input.grid.length_m / rt;
        assert!((last.receipt.gas_mass_kg_per_m + exchanged - initial_mass).abs() < 1.0e-15);
        let error = last
            .checkpoint
            .active_absolute_pressure_pa
            .iter()
            .zip(reference)
            .map(|(p, expected)| (p - expected).abs())
            .fold(0.0_f64, f64::max);
        errors.push(error);
    }
    for pair in errors.windows(2) {
        let ratio = pair[0] / pair[1];
        assert!(
            ratio > 1.8 && ratio < 2.1,
            "first-order pressure convergence: {errors:?}"
        );
    }
}

#[test]
fn g1_restricted_vent_interpolates_sealed_and_ideal_vent_squeeze_pressure() {
    let mut sealed = fixture(4);
    sealed.grid.gap_m.fill(1.0e-6);
    sealed.wall_motion.gap_rate_m_per_s = -1.0e-5;
    sealed.timestep_s = 1.0e-3;
    sealed.slip_policy = SlipPolicy::MaxwellFirstOrder {
        source_id: "declared-film-only-accommodation".into(),
        tangential_momentum_accommodation: 0.9,
    };
    sealed.applicability.maximum_knudsen_number = 0.1;
    sealed.budget.mass_residual_tolerance_kg_m2_s = 1.0e-12;
    let sealed_step = solve_isothermal_gas_film_1d(&sealed, None).unwrap();
    let bath = sealed.initial_absolute_pressure_pa;
    let mut pinned = sealed.clone();
    pinned.boundary = GasFilmBoundaryTopology::Vented {
        cell_index: 0,
        absolute_pressure_pa: bath,
    };
    let pinned_step = solve_isothermal_gas_film_1d(&pinned, None).unwrap();
    let mut previous = sealed_step.checkpoint.active_absolute_pressure_pa[0];
    for count in [1.0, 100.0, 20_000.0] {
        let mut finite = sealed.clone();
        finite.boundary = capillary_fixture(4, bath, count).boundary;
        let step = solve_isothermal_gas_film_1d(&finite, None).unwrap();
        let p = step.checkpoint.active_absolute_pressure_pa[0];
        assert!(p < previous && p > bath);
        assert!(step.receipt.vent_outward_mass_flux_kg_per_m_s > 0.0);
        assert!(
            step.receipt.gauge_pressure_load_n_per_m
                > pinned_step.receipt.gauge_pressure_load_n_per_m
        );
        assert!(
            step.receipt.gauge_pressure_load_n_per_m
                < sealed_step.receipt.gauge_pressure_load_n_per_m
        );
        assert!(step.receipt.mass_closure_residual_kg_per_m_s.abs() < 2.0e-14);
        previous = p;
    }
    assert!(
        previous - bath < 0.1 * (sealed_step.checkpoint.active_absolute_pressure_pa[0] - bath),
        "restricted vent pressure={previous}, bath={bath}, sealed pressure={}",
        sealed_step.checkpoint.active_absolute_pressure_pa[0]
    );
}

#[test]
fn g4_capillary_admits_its_own_no_slip_geometry_and_endpoint_flow_regime() {
    let input = capillary_fixture(4, 100_325.0, 1_000.0);
    for field in [0, 1, 2] {
        let mut bad = input.clone();
        let GasFilmBoundaryTopology::CapillaryVented { capillary, .. } = &mut bad.boundary else {
            unreachable!()
        };
        match field {
            0 => capillary.radius_m = 0.0,
            1 => capillary.length_m = f64::NAN,
            _ => capillary.capillaries_per_m = -1.0,
        }
        assert!(matches!(
            solve_isothermal_gas_film_1d(&bad, None),
            Err(GasFilmError::InvalidInput { .. })
        ));
    }
    let mut short = input.clone();
    let GasFilmBoundaryTopology::CapillaryVented { capillary, .. } = &mut short.boundary else {
        unreachable!()
    };
    capillary.length_m = 10.0 * capillary.radius_m;
    assert_eq!(
        solve_isothermal_gas_film_1d(&short, None),
        Err(GasFilmError::Unavailable {
            reason: "capillary-outside-long-tube-envelope",
        })
    );
    let mut fast = capillary_fixture(4, 80_000.0, 1_000.0);
    assert_eq!(
        solve_isothermal_gas_film_1d(&fast, None),
        Err(GasFilmError::Unavailable {
            reason: "capillary-reynolds-outside-creeping-envelope",
        })
    );
    fast.boundary = capillary_fixture(4, 10_000.0, 1_000.0).boundary;
    fast.slip_policy = SlipPolicy::MaxwellFirstOrder {
        source_id: "film-accommodation-does-not-admit-tube-slip".into(),
        tangential_momentum_accommodation: 0.9,
    };
    fast.applicability.maximum_knudsen_number = 0.1;
    assert_eq!(
        solve_isothermal_gas_film_1d(&fast, None),
        Err(GasFilmError::Unavailable {
            reason: "capillary-knudsen-outside-no-slip-envelope",
        })
    );
    // A mean-speed gate would admit this; the centreline speed is twice as large.
    let mut mach = input.clone();
    let p = mach.initial_absolute_pressure_pa;
    let bath = 100_325.0;
    let mean_speed = 1.0e-5_f64.powi(2) * (p * p - bath * bath)
        / (16.0 * mach.gas.dynamic_viscosity_pa_s * 0.01 * bath);
    mach.applicability.maximum_mach_number =
        1.5 * mean_speed / mach.applicability.speed_of_sound_m_per_s;
    assert_eq!(
        solve_isothermal_gas_film_1d(&mach, None),
        Err(GasFilmError::Unavailable {
            reason: "capillary-mach-outside-quasisteady-envelope",
        })
    );
}

#[test]
fn g3_capillary_restart_keeps_geometry_but_allows_reversed_bath_pressure() {
    let input = capillary_fixture(4, 100_325.0, 1_000.0);
    let first = solve_isothermal_gas_film_1d(&input, None).unwrap();
    let saved = first.checkpoint.clone();
    let mut reverse = input.clone();
    let GasFilmBoundaryTopology::CapillaryVented {
        absolute_pressure_pa,
        ..
    } = &mut reverse.boundary
    else {
        unreachable!()
    };
    *absolute_pressure_pa = 102_325.0;
    let second = solve_isothermal_gas_film_1d(&reverse, Some(&saved)).unwrap();
    assert!(first.receipt.vent_outward_mass_flux_kg_per_m_s > 0.0);
    assert!(second.receipt.vent_outward_mass_flux_kg_per_m_s < 0.0);
    assert_eq!(
        second,
        solve_isothermal_gas_film_1d(&reverse, Some(&saved)).unwrap()
    );
    for field in [0, 1, 2, 3] {
        let mut changed = input.clone();
        let GasFilmBoundaryTopology::CapillaryVented {
            cell_index,
            capillary,
            ..
        } = &mut changed.boundary
        else {
            unreachable!()
        };
        match field {
            0 => *cell_index = 1,
            1 => capillary.radius_m *= 1.1,
            2 => capillary.length_m *= 1.1,
            _ => capillary.capillaries_per_m *= 1.1,
        }
        assert!(matches!(
            solve_isothermal_gas_film_1d(&changed, Some(&saved)),
            Err(GasFilmError::CheckpointMismatch {
                field: "configuration_fingerprint"
            })
        ));
    }
    assert_eq!(saved, first.checkpoint);
}

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
            gap_rate_profile_m_per_s: None,
        },
        initial_absolute_pressure_pa: pressure,
        gauge_reference_absolute_pressure_pa: pressure,
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
    input.boundary = GasFilmBoundaryTopology::Open {
        left_absolute_pressure_pa: 101_325.0,
        right_absolute_pressure_pa: 101_325.0,
    };
    input.wall_motion.upper_tangential_velocity_m_per_s = 2.0;
    let step = solve_isothermal_gas_film_1d(&input, None).expect("Couette limit admitted");
    let expected_shear = input.gas.dynamic_viscosity_pa_s * 2.0 / input.grid.gap_m[0];
    for pressure in active_pressure(&step) {
        assert!((pressure - input.initial_absolute_pressure_pa).abs() < 1.0e-9);
    }
    for shear in &step.receipt.upper_wall_shear_pa {
        assert!((shear.expect("active") - expected_shear).abs() < 1.0e-12);
    }
    for index in 0..input.grid.gap_m.len() {
        let wall_on_gas = step.receipt.upper_wall_on_gas_tangential_traction_pa[index]
            .expect("active upper-wall traction");
        assert_eq!(
            step.receipt.gas_on_upper_wall_tangential_traction_pa[index],
            Some(-wall_on_gas)
        );
        let gas_on_wall = step.receipt.gas_on_upper_wall_normal_traction_pa[index]
            .expect("active gas pressure traction");
        assert_eq!(
            step.receipt.upper_wall_on_gas_normal_traction_pa[index],
            Some(-gas_on_wall)
        );
        assert_eq!(step.receipt.upper_wall_shear_pa[index], Some(wall_on_gas));
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
    coarse.boundary = GasFilmBoundaryTopology::Open {
        left_absolute_pressure_pa: coarse.initial_absolute_pressure_pa,
        right_absolute_pressure_pa: coarse.initial_absolute_pressure_pa,
    };
    let mut fine = fixture(32);
    fine.wall_motion.upper_tangential_velocity_m_per_s = 1.25;
    fine.boundary = coarse.boundary.clone();
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
            (pressure - expected_pressure).abs() / expected_pressure < 1.0e-8,
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
            < 1.0e-8,
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
            (actual - expected).abs() / expected < 1.0e-8,
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
    let mirrored_negative_pressure: Vec<f64> =
        active_pressure(&negative).into_iter().rev().collect();
    assert_eq!(active_pressure(&positive), mirrored_negative_pressure);
    for (forward_shear, reversed_shear) in positive
        .receipt
        .upper_wall_shear_pa
        .iter()
        .zip(negative.receipt.upper_wall_shear_pa.iter().rev())
    {
        assert_eq!(
            forward_shear.expect("active"),
            -reversed_shear.expect("active")
        );
    }
    assert!(
        (positive.receipt.wall_power_to_gas_w_per_m - negative.receipt.wall_power_to_gas_w_per_m)
            .abs()
            < 1.0e-12
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
    assert_eq!(
        step.receipt.gas_on_lower_wall_tangential_traction_pa[5],
        None
    );
    assert_eq!(
        step.receipt.lower_wall_on_gas_tangential_traction_pa[5],
        None
    );
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
fn pressure_scaled_knudsen_checks_reservoirs_and_vent() {
    for boundary in [
        GasFilmBoundaryTopology::Open {
            left_absolute_pressure_pa: 50_000.0,
            right_absolute_pressure_pa: 101_325.0,
        },
        GasFilmBoundaryTopology::Open {
            left_absolute_pressure_pa: 101_325.0,
            right_absolute_pressure_pa: 50_000.0,
        },
        GasFilmBoundaryTopology::Vented {
            cell_index: 2,
            absolute_pressure_pa: 50_000.0,
        },
    ] {
        let mut input = fixture(4);
        input.boundary = boundary;
        assert_eq!(
            solve_isothermal_gas_film_1d(&input, None),
            Err(GasFilmError::Unavailable {
                reason: "rarefied-knudsen-outside-continuum-envelope"
            })
        );
    }

    let mut input = fixture(4);
    input.boundary = GasFilmBoundaryTopology::Open {
        left_absolute_pressure_pa: 110_000.0,
        right_absolute_pressure_pa: 80_000.0,
    };
    let step = solve_isothermal_gas_film_1d(&input, None).expect("continuum reservoir");
    let expected = input.applicability.mean_free_path_m / input.grid.gap_m[3]
        * input.initial_absolute_pressure_pa
        / 80_000.0;
    assert!((step.receipt.maximum_knudsen_number / expected - 1.0).abs() < 1.0e-12);
}

#[test]
fn sealed_gap_changes_scale_mean_free_path_with_pressure() {
    for new_gap in [5.0e-6, 20.0e-6] {
        let mut input = fixture(4);
        let old_gap = input.grid.gap_m[0];
        input.grid.gap_m.fill(new_gap);
        input.wall_motion.gap_rate_m_per_s = (new_gap - old_gap) / input.timestep_s;
        let step = solve_isothermal_gas_film_1d(&input, None)
            .expect("isothermal compression/expansion preserves the continuum regime");
        // In a sealed uniform film, p*h is conserved. Both lambda and h scale
        // together, so Kn stays fixed even when a pressure-blind check would
        // reject the compressed gap or understate the expanded gap's Kn.
        let expected_knudsen = input.applicability.mean_free_path_m / old_gap;
        let expected_pressure = input.initial_absolute_pressure_pa * old_gap / new_gap;
        assert!((step.receipt.maximum_knudsen_number / expected_knudsen - 1.0).abs() < 1.0e-8);
        for pressure in active_pressure(&step) {
            assert!((pressure / expected_pressure - 1.0).abs() < 1.0e-8);
        }
    }

    let input = fixture(4);
    let mut checkpoint = solve_isothermal_gas_film_1d(&input, None)
        .expect("continuum checkpoint")
        .checkpoint;
    checkpoint.active_absolute_pressure_pa[1] = 50_000.0;
    let before = checkpoint.clone();
    assert_eq!(
        solve_isothermal_gas_film_1d(&input, Some(&checkpoint)),
        Err(GasFilmError::Unavailable {
            reason: "rarefied-knudsen-outside-continuum-envelope"
        })
    );
    assert_eq!(checkpoint, before);
}

#[test]
fn converged_pressure_must_remain_inside_the_continuum_envelope() {
    let mut input = fixture(2);
    input.timestep_s = 4.0e-5;
    input.wall_motion.lower_tangential_velocity_m_per_s = 100.0;
    input.wall_motion.upper_tangential_velocity_m_per_s = 100.0;
    // The initial state and wall speeds are admitted. Couette transport moves
    // mass toward the sealed right end, rarefying the left cell during the step.
    assert_eq!(
        solve_isothermal_gas_film_1d(&input, None),
        Err(GasFilmError::Unavailable {
            reason: "rarefied-knudsen-outside-continuum-envelope"
        })
    );
}

#[test]
fn stationary_walls_do_not_bypass_pressure_driven_mach_limit() {
    let mut input = fixture(2);
    input.boundary = GasFilmBoundaryTopology::Open {
        left_absolute_pressure_pa: 101_325.0,
        right_absolute_pressure_pa: 101_325.0,
    };
    let mut checkpoint = solve_isothermal_gas_film_1d(&input, None)
        .expect("quiescent checkpoint")
        .checkpoint;
    let left: f64 = 2.0e6;
    let right: f64 = 101_325.0;
    input.boundary = GasFilmBoundaryTopology::Open {
        left_absolute_pressure_pa: left,
        right_absolute_pressure_pa: right,
    };
    // The exact steady isothermal profile has p^2 linear in x. It satisfies
    // the mass balance but its Poiseuille speed exceeds the declared Mach cap.
    checkpoint.active_absolute_pressure_pa = [0.25, 0.75]
        .map(|x| (left * left + (right * right - left * left) * x).sqrt())
        .to_vec();
    let before = checkpoint.clone();
    assert_eq!(
        solve_isothermal_gas_film_1d(&input, Some(&checkpoint)),
        Err(GasFilmError::Unavailable {
            reason: "gas-mach-outside-declared-envelope"
        })
    );
    assert_eq!(checkpoint, before);
}

#[test]
fn couette_poiseuille_mach_uses_off_center_peak_in_both_flow_directions() {
    for direction in [1.0, -1.0] {
        let mut input = fixture(2);
        input.boundary = GasFilmBoundaryTopology::Open {
            left_absolute_pressure_pa: if direction > 0.0 {
                120_000.0
            } else {
                100_000.0
            },
            right_absolute_pressure_pa: if direction > 0.0 {
                100_000.0
            } else {
                120_000.0
            },
        };
        input.wall_motion.lower_tangential_velocity_m_per_s = direction;
        input.wall_motion.upper_tangential_velocity_m_per_s = 3.0 * direction;
        let step = solve_isothermal_gas_film_1d(&input, None).expect("subsonic pressure flow");
        let inlet = if direction > 0.0 { 0 } else { 1 };
        let pressure = step.absolute_pressure_pa[inlet].expect("active inlet cell");
        let gradient_magnitude = (120_000.0 - pressure) / (0.25 * input.grid.length_m);
        let curvature = input.grid.gap_m[inlet].powi(2) * gradient_magnitude
            / (2.0 * input.gas.dynamic_viscosity_pa_s);
        // Completing the square in the quadratic velocity profile gives the
        // interior peak: |U_mean| + |A|/4 + (U_upper-U_lower)^2/(4|A|).
        let midpoint_speed = 2.0 + curvature / 4.0;
        let expected_peak_speed = midpoint_speed + 1.0 / curvature;
        let actual_peak_speed =
            step.receipt.maximum_gas_mach_number * input.applicability.speed_of_sound_m_per_s;
        assert!(
            curvature > 2.0,
            "the profile extremum must be inside the gap"
        );
        assert!((actual_peak_speed / expected_peak_speed - 1.0).abs() < 1.0e-12);

        // A cap between the midpoint and the true peak must refuse even though
        // both wall velocities and the entire mean-flow estimate pass.
        input.applicability.maximum_mach_number = 0.5 * (midpoint_speed + expected_peak_speed)
            / input.applicability.speed_of_sound_m_per_s;
        assert_eq!(
            solve_isothermal_gas_film_1d(&input, None),
            Err(GasFilmError::Unavailable {
                reason: "gas-mach-outside-declared-envelope"
            })
        );
    }
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

#[test]
fn checkpoint_fingerprint_rejects_each_semantic_mutation_group() {
    let input = fixture(8);
    let checkpoint = solve_isothermal_gas_film_1d(&input, None)
        .expect("base checkpoint")
        .checkpoint;
    let mut variants = Vec::new();
    let mut identity = input.clone();
    identity.identity.frame_id = "other-fixture-line".to_owned();
    variants.push(identity);
    let mut gas = input.clone();
    gas.gas.dynamic_viscosity_pa_s *= 1.1;
    variants.push(gas);
    let mut topology = input.clone();
    topology.boundary = GasFilmBoundaryTopology::Open {
        left_absolute_pressure_pa: 101_325.0,
        right_absolute_pressure_pa: 101_325.0,
    };
    variants.push(topology);
    let mut slip_and_roughness = input.clone();
    slip_and_roughness.slip_policy = SlipPolicy::NoSlipContinuum {
        source_id: "other-no-slip-source".to_owned(),
    };
    slip_and_roughness.roughness_policy = RoughnessPolicy::ResolvedSmooth {
        source_id: "other-roughness-source".to_owned(),
        maximum_roughness_m: 2.0e-8,
    };
    variants.push(slip_and_roughness);
    let mut applicability = input.clone();
    applicability.applicability.maximum_gap_slope = 0.2;
    applicability.uncertainty.gap_relative_bound = 0.01;
    variants.push(applicability);
    let mut gas_variant = input.clone();
    gas_variant.identity.gas_species_id = "mutated-gas-species-v2".to_owned();
    variants.push(gas_variant);
    let mut initial_and_gauge = input.clone();
    initial_and_gauge.gauge_reference_absolute_pressure_pa = 99_000.0;
    variants.push(initial_and_gauge);
    let mut timestep_and_budget = input.clone();
    timestep_and_budget.timestep_s *= 2.0;
    timestep_and_budget.budget.relaxation = 0.7;
    variants.push(timestep_and_budget);
    for changed in variants {
        assert_eq!(
            solve_isothermal_gas_film_1d(&changed, Some(&checkpoint)),
            Err(GasFilmError::CheckpointMismatch {
                field: "configuration_fingerprint"
            })
        );
    }
}

#[test]
fn restart_admits_physical_gap_evolution_and_changed_tangential_speed() {
    let mut first_input = fixture(8);
    first_input.wall_motion.gap_rate_m_per_s = -1.0e-4;
    let first = solve_isothermal_gas_film_1d(&first_input, None).expect("first squeeze step");
    let mut second_input = first_input.clone();
    second_input.grid.gap_m = first_input
        .grid
        .gap_m
        .iter()
        .map(|gap| gap + second_input.wall_motion.gap_rate_m_per_s * second_input.timestep_s)
        .collect();
    second_input.wall_motion.upper_tangential_velocity_m_per_s = 0.5;
    let second = solve_isothermal_gas_film_1d(&second_input, Some(&first.checkpoint))
        .expect("consistent second squeeze step");
    assert_eq!(second.checkpoint.step_index, 2);
}

#[test]
fn restart_admits_nonuniform_gap_evolution_and_uses_cellwise_pressure_work() {
    let mut first_input = fixture(8);
    let rates = (0..8)
        .map(|index| -2.0e-5 - index as f64 * 1.0e-6)
        .collect::<Vec<_>>();
    first_input.wall_motion.gap_rate_profile_m_per_s = Some(rates.clone());
    let first =
        solve_isothermal_gas_film_1d(&first_input, None).expect("first changing-wedge step");

    let spacing = first_input.grid.length_m / first_input.grid.gap_m.len() as f64;
    let reconstructed_normal_power = first
        .absolute_pressure_pa
        .iter()
        .zip(&rates)
        .map(|(pressure, rate)| -pressure.expect("active pressure") * rate * spacing)
        .sum::<f64>();
    assert!(
        (first.receipt.normal_gap_power_to_gas_w_per_m - reconstructed_normal_power).abs()
            <= 1.0e-12 * reconstructed_normal_power.abs().max(1.0)
    );

    let mut second_input = first_input.clone();
    second_input.grid.gap_m = first_input
        .grid
        .gap_m
        .iter()
        .zip(&rates)
        .map(|(gap, rate)| gap + rate * second_input.timestep_s)
        .collect();
    let second = solve_isothermal_gas_film_1d(&second_input, Some(&first.checkpoint))
        .expect("nonuniform physical restart");
    assert_eq!(second.checkpoint.step_index, 2);
}

#[test]
fn per_cell_gap_rate_requires_one_finite_value_per_geometric_cell() {
    let mut wrong_length = fixture(8);
    wrong_length.wall_motion.gap_rate_profile_m_per_s = Some(vec![0.0; 7]);
    assert_eq!(
        solve_isothermal_gas_film_1d(&wrong_length, None),
        Err(GasFilmError::InvalidInput {
            field: "wall_motion.gap_rate_profile_m_per_s.length"
        })
    );

    let mut nonfinite = fixture(8);
    let mut rates = vec![0.0; 8];
    rates[3] = f64::NAN;
    nonfinite.wall_motion.gap_rate_profile_m_per_s = Some(rates);
    assert_eq!(
        solve_isothermal_gas_film_1d(&nonfinite, None),
        Err(GasFilmError::InvalidInput {
            field: "wall_motion.gap_rate_profile_m_per_s"
        })
    );
}

#[test]
fn restart_rejects_inconsistent_gap_evolution_without_mutating_checkpoint() {
    let mut input = fixture(8);
    input.wall_motion.gap_rate_m_per_s = -1.0e-4;
    let first = solve_isothermal_gas_film_1d(&input, None).expect("first squeeze step");
    let before = first.checkpoint.clone();
    let mut inconsistent = input.clone();
    inconsistent.grid.gap_m[0] += 1.0e-7;
    assert_eq!(
        solve_isothermal_gas_film_1d(&inconsistent, Some(&first.checkpoint)),
        Err(GasFilmError::CheckpointMismatch {
            field: "active_gap_evolution"
        })
    );
    assert_eq!(first.checkpoint, before);
}

#[test]
fn gauge_reference_and_mass_flux_sign_receipts_are_independent() {
    let mut uniform = fixture(6);
    uniform.gauge_reference_absolute_pressure_pa = 100_000.0;
    let uniform_step = solve_isothermal_gas_film_1d(&uniform, None).expect("uniform gauge case");
    assert!(
        (uniform_step.receipt.gauge_pressure_load_n_per_m
            - (101_325.0 - 100_000.0) * uniform.grid.length_m)
            .abs()
            < 1.0e-12
    );

    let mut vented = fixture(8);
    vented.boundary = GasFilmBoundaryTopology::Vented {
        cell_index: 3,
        absolute_pressure_pa: 90_000.0,
    };
    let step = solve_isothermal_gas_film_1d(&vented, None).expect("vented sign case");
    assert!(step.receipt.vent_outward_mass_flux_kg_per_m_s > 0.0);
    let independently_summed = step.receipt.storage_rate_kg_per_m_s
        + step.receipt.left_boundary_outward_mass_flux_kg_per_m_s
        + step.receipt.right_boundary_outward_mass_flux_kg_per_m_s
        + step.receipt.vent_outward_mass_flux_kg_per_m_s;
    assert!((independently_summed - step.receipt.mass_closure_residual_kg_per_m_s).abs() < 1.0e-12);
}

fn assert_mechanical_balance(input: &GasFilmInput, step: &fs_flux::GasFilmStep) {
    let spacing = input.grid.length_m / input.grid.gap_m.len() as f64;
    let mut tangential_work = 0.0;
    for cell in 0..step.checkpoint.active_cells {
        let upper = step.receipt.upper_wall_on_gas_tangential_traction_pa[cell].expect("active");
        let lower = step.receipt.lower_wall_on_gas_tangential_traction_pa[cell].expect("active");
        assert_eq!(
            step.receipt.gas_on_upper_wall_tangential_traction_pa[cell],
            Some(-upper)
        );
        assert_eq!(
            step.receipt.gas_on_lower_wall_tangential_traction_pa[cell],
            Some(-lower)
        );
        tangential_work += (upper * input.wall_motion.upper_tangential_velocity_m_per_s
            + lower * input.wall_motion.lower_tangential_velocity_m_per_s)
            * spacing;
    }
    let wall = step.receipt.wall_power_to_gas_w_per_m;
    let pressure = step.receipt.pressure_gradient_power_to_gas_w_per_m;
    let normal = step.receipt.normal_gap_power_to_gas_w_per_m;
    let heat = step.receipt.viscous_heat_w_per_m;
    let tolerance = 1.0e-11 * heat.abs().max(1.0e-10);
    assert!((wall - normal - tangential_work).abs() <= tolerance);
    assert!((tangential_work + pressure - heat).abs() <= tolerance);
    assert!(heat >= 0.0);
}

#[test]
fn g1_steady_compressible_poiseuille_has_full_shear_and_convergent_heat() {
    let mut relative_errors = Vec::new();
    for cells in [8, 32] {
        let mut input = fixture(cells);
        let left: f64 = 120_000.0;
        let right: f64 = 100_000.0;
        input.boundary = GasFilmBoundaryTopology::Open {
            left_absolute_pressure_pa: left,
            right_absolute_pressure_pa: right,
        };
        input.budget.maximum_iterations = 20_000;
        // Seed the exact steady ideal-gas solution p^2(x), independent of the
        // iteration. The FV face flux is exact for this profile and constant h.
        let mut state = solve_isothermal_gas_film_1d(&input, None)
            .expect("open checkpoint")
            .checkpoint;
        state.active_absolute_pressure_pa = (0..cells)
            .map(|cell| {
                let x_over_length = (cell as f64 + 0.5) / cells as f64;
                (left * left + (right * right - left * left) * x_over_length).sqrt()
            })
            .collect();
        let step = solve_isothermal_gas_film_1d(&input, Some(&state))
            .expect("exact steady Poiseuille restart");
        let h = input.grid.gap_m[0];
        let mu = input.gas.dynamic_viscosity_pa_s;
        let rt = input.gas.specific_gas_constant_j_kg_k * input.gas.temperature_k;
        let expected_mass =
            h.powi(3) * (left * left - right * right) / (24.0 * mu * rt * input.grid.length_m);
        assert!(
            (step.receipt.right_boundary_outward_mass_flux_kg_per_m_s / expected_mass - 1.0).abs()
                < 1.0e-10
        );
        assert!(
            (-step.receipt.left_boundary_outward_mass_flux_kg_per_m_s / expected_mass - 1.0).abs()
                < 1.0e-10
        );
        let spacing = input.grid.length_m / cells as f64;
        for tractions in [
            &step.receipt.upper_wall_on_gas_tangential_traction_pa,
            &step.receipt.lower_wall_on_gas_tangential_traction_pa,
        ] {
            assert!(tractions.iter().all(|value| value.expect("active") < 0.0));
            let force = tractions
                .iter()
                .map(|value| value.expect("active") * spacing)
                .sum::<f64>();
            let expected_force = 0.5 * h * (right - left);
            assert!((force / expected_force - 1.0).abs() < 1.0e-12);
        }
        // Integral of h^3 p_x^2/(12 mu) for p^2 linear in x:
        // heat = mass_flux * R*T * log(p_left/p_right).
        let expected_heat = expected_mass * rt * (left / right).ln();
        let relative_error = (step.receipt.viscous_heat_w_per_m / expected_heat - 1.0).abs();
        assert!(step.receipt.viscous_heat_w_per_m > 0.0);
        assert_eq!(step.receipt.wall_power_to_gas_w_per_m, 0.0);
        assert_mechanical_balance(&input, &step);
        relative_errors.push(relative_error);
    }
    assert!(relative_errors[1] < 3.0e-6, "{relative_errors:?}");
    assert!(
        relative_errors[1] < relative_errors[0] / 10.0,
        "{relative_errors:?}"
    );
}

#[test]
fn g1_sealed_comoving_walls_have_backpressure_shear_and_dissipation() {
    let mut input = fixture(8);
    let velocity = 1.0;
    input.wall_motion.lower_tangential_velocity_m_per_s = velocity;
    input.wall_motion.upper_tangential_velocity_m_per_s = velocity;
    let h = input.grid.gap_m[0];
    let mu = input.gas.dynamic_viscosity_pa_s;
    let gradient = 12.0 * mu * velocity / h.powi(2);
    let mut state = solve_isothermal_gas_film_1d(&input, None)
        .expect("sealed checkpoint")
        .checkpoint;
    // This linear pressure cancels Couette transport at every face, including
    // the two sealed boundary half-cells, despite zero relative wall speed.
    state.active_absolute_pressure_pa = (0..input.grid.gap_m.len())
        .map(|cell| {
            let x = (cell as f64 + 0.5) * input.grid.length_m / input.grid.gap_m.len() as f64;
            input.initial_absolute_pressure_pa + gradient * (x - 0.5 * input.grid.length_m)
        })
        .collect();
    let step = solve_isothermal_gas_film_1d(&input, Some(&state)).expect("sealed steady state");
    let expected_shear = 6.0 * mu * velocity / h;
    for tractions in [
        &step.receipt.upper_wall_on_gas_tangential_traction_pa,
        &step.receipt.lower_wall_on_gas_tangential_traction_pa,
    ] {
        for traction in tractions {
            assert!((traction.expect("active") / expected_shear - 1.0).abs() < 1.0e-11);
        }
    }
    let expected_heat = 12.0 * mu * velocity.powi(2) * input.grid.length_m / h;
    assert!((step.receipt.viscous_heat_w_per_m / expected_heat - 1.0).abs() < 1.0e-11);
    assert!((step.receipt.wall_power_to_gas_w_per_m / expected_heat - 1.0).abs() < 1.0e-11);
    assert!(step.receipt.pressure_gradient_power_to_gas_w_per_m.abs() < 1.0e-12);
    assert_mechanical_balance(&input, &step);
}

#[test]
fn g1_pressure_flow_can_drive_both_moving_walls() {
    let mut input = fixture(8);
    input.boundary = GasFilmBoundaryTopology::Open {
        left_absolute_pressure_pa: 120_000.0,
        right_absolute_pressure_pa: 100_000.0,
    };
    input.wall_motion.lower_tangential_velocity_m_per_s = 1.0;
    input.wall_motion.upper_tangential_velocity_m_per_s = 1.0;
    let step =
        solve_isothermal_gas_film_1d(&input, None).expect("pressure-driven translating walls");
    // Pressure traction integrates to h*(p_right-p_left); both walls translate
    // at one m/s, so gas delivers mechanical work to the walls.
    let expected_wall_work = input.grid.gap_m[0] * (100_000.0 - 120_000.0);
    assert!((step.receipt.wall_power_to_gas_w_per_m / expected_wall_work - 1.0).abs() < 1.0e-12);
    assert!(step.receipt.wall_power_to_gas_w_per_m < 0.0);
    assert!(
        step.receipt.pressure_gradient_power_to_gas_w_per_m > step.receipt.viscous_heat_w_per_m
    );
    assert_mechanical_balance(&input, &step);
}

#[test]
fn g1_opposing_pressure_gradients_do_not_cancel_viscous_heat() {
    let input = fixture(3);
    let mut state = solve_isothermal_gas_film_1d(&input, None)
        .expect("sealed checkpoint")
        .checkpoint;
    state.active_absolute_pressure_pa = vec![100_000.0, 120_000.0, 100_000.0];
    let step = solve_isothermal_gas_film_1d(&input, Some(&state)).expect("symmetric pressure peak");
    let p = active_pressure(&step);
    assert_eq!(p[0], p[2]);
    assert!(p[1] > p[0]);
    let spacing = input.grid.length_m / 3.0;
    let gradient = (p[1] - p[0]) / spacing;
    let expected_heat = 2.0 * input.grid.gap_m[0].powi(3) * gradient.powi(2) * spacing
        / (12.0 * input.gas.dynamic_viscosity_pa_s);
    assert_eq!(step.receipt.upper_wall_shear_pa[1], Some(0.0));
    assert!(step.receipt.viscous_heat_w_per_m > 0.0);
    assert!((step.receipt.viscous_heat_w_per_m / expected_heat - 1.0).abs() < 1.0e-12);
    assert_mechanical_balance(&input, &step);
}

fn maxwell_fixture(cells: usize, accommodation: f64) -> GasFilmInput {
    let mut input = fixture(cells);
    // Lower pressure for the same declared isothermal gas moves this fixture
    // into the slip regime: lambda*p stays fixed, and Kn is about 0.033.
    input.applicability.mean_free_path_m *= input.initial_absolute_pressure_pa / 20_000.0;
    input.initial_absolute_pressure_pa = 20_000.0;
    input.gauge_reference_absolute_pressure_pa = 20_000.0;
    input.gas.declared_density_kg_m3 = input.initial_absolute_pressure_pa
        / (input.gas.specific_gas_constant_j_kg_k * input.gas.temperature_k);
    input.applicability.maximum_knudsen_number = 0.1;
    input.slip_policy = SlipPolicy::MaxwellFirstOrder {
        source_id: "synthetic-shared-wall-maxwell-v1".to_owned(),
        tangential_momentum_accommodation: accommodation,
    };
    input
}

#[test]
fn maxwell_couette_has_analytic_shear_and_bulk_plus_wall_dissipation() {
    for sigma in [1.0, 0.5] {
        let mut input = maxwell_fixture(6, sigma);
        input.boundary = GasFilmBoundaryTopology::Open {
            left_absolute_pressure_pa: input.initial_absolute_pressure_pa,
            right_absolute_pressure_pa: input.initial_absolute_pressure_pa,
        };
        let speed = 2.0;
        input.wall_motion.upper_tangential_velocity_m_per_s = speed;
        let step = solve_isothermal_gas_film_1d(&input, None).expect("Maxwell Couette flow");
        let h = input.grid.gap_m[0];
        let mu = input.gas.dynamic_viscosity_pa_s;
        let b = (2.0 - sigma) / sigma * input.applicability.mean_free_path_m;
        let effective_gap = h + 2.0 * b;
        let shear = mu * speed / effective_gap;
        let total_heat = shear * speed * input.grid.length_m;
        let bulk_heat = mu * (speed / effective_gap).powi(2) * h * input.grid.length_m;
        let slip_heat = 2.0 * mu * b * (speed / effective_gap).powi(2) * input.grid.length_m;
        for pressure in active_pressure(&step) {
            assert!((pressure / input.initial_absolute_pressure_pa - 1.0).abs() < 1.0e-12);
        }
        for cell in 0..input.grid.gap_m.len() {
            assert!(
                (step.receipt.upper_wall_shear_pa[cell].unwrap() / shear - 1.0).abs() < 1.0e-12
            );
            assert!(
                (step.receipt.lower_wall_on_gas_tangential_traction_pa[cell].unwrap() / -shear
                    - 1.0)
                    .abs()
                    < 1.0e-12
            );
        }
        assert!((step.receipt.viscous_heat_w_per_m / total_heat - 1.0).abs() < 1.0e-12);
        assert!((step.receipt.wall_slip_dissipation_w_per_m / slip_heat - 1.0).abs() < 1.0e-12);
        assert!(
            ((step.receipt.viscous_heat_w_per_m - step.receipt.wall_slip_dissipation_w_per_m)
                / bulk_heat
                - 1.0)
                .abs()
                < 1.0e-12
        );
        let peak_gas_speed = speed * (h + b) / effective_gap;
        assert!(
            (step.receipt.maximum_gas_mach_number * input.applicability.speed_of_sound_m_per_s
                / peak_gas_speed
                - 1.0)
                .abs()
                < 1.0e-12
        );
        assert_mechanical_balance(&input, &step);
    }
}

#[test]
fn maxwell_poiseuille_matches_pressure_dependent_slip_flux_and_heat_limit() {
    let sigma = 0.8;
    let left: f64 = 30_000.0;
    let right: f64 = 20_000.0;
    let mut heat_errors = Vec::new();
    for cells in [8, 32] {
        let mut input = maxwell_fixture(cells, sigma);
        input.boundary = GasFilmBoundaryTopology::Open {
            left_absolute_pressure_pa: right,
            right_absolute_pressure_pa: right,
        };
        let mut checkpoint = solve_isothermal_gas_film_1d(&input, None)
            .expect("quiescent slip checkpoint")
            .checkpoint;
        input.boundary = GasFilmBoundaryTopology::Open {
            left_absolute_pressure_pa: left,
            right_absolute_pressure_pa: right,
        };
        let h = input.grid.gap_m[0];
        let mu = input.gas.dynamic_viscosity_pa_s;
        let rt = input.gas.specific_gas_constant_j_kg_k * input.gas.temperature_k;
        let bp = (2.0 - sigma) / sigma
            * input.applicability.mean_free_path_m
            * input.initial_absolute_pressure_pa;
        let shift = 6.0 * bp / h;
        // With b(p)*p constant, the exact steady solution has (p+shift)^2
        // linear in x, not p^2 and not the constant-slip-length solution.
        checkpoint.active_absolute_pressure_pa = (0..cells)
            .map(|cell| {
                let x = (cell as f64 + 0.5) / cells as f64;
                ((left + shift).powi(2) + ((right + shift).powi(2) - (left + shift).powi(2)) * x)
                    .sqrt()
                    - shift
            })
            .collect();
        let expected_pressure = checkpoint.active_absolute_pressure_pa.clone();
        let step = solve_isothermal_gas_film_1d(&input, Some(&checkpoint))
            .expect("steady compressible Maxwell Poiseuille flow");
        let expected_mass = h.powi(3) * (left - right) * (left + right + 2.0 * shift)
            / (24.0 * mu * rt * input.grid.length_m);
        assert!(
            (step.receipt.right_boundary_outward_mass_flux_kg_per_m_s / expected_mass - 1.0).abs()
                < 1.0e-10
        );
        assert!(
            (-step.receipt.left_boundary_outward_mass_flux_kg_per_m_s / expected_mass - 1.0).abs()
                < 1.0e-10
        );
        for (actual, expected) in active_pressure(&step).iter().zip(expected_pressure) {
            assert!((actual / expected - 1.0).abs() < 1.0e-12);
        }
        // Pressure stress remains Newtonian. Increased throughput is not an
        // invented effective viscosity or a change to the imposed pressure drop.
        let force = step
            .receipt
            .upper_wall_on_gas_tangential_traction_pa
            .iter()
            .map(|traction| traction.unwrap() * input.grid.length_m / cells as f64)
            .sum::<f64>();
        assert!((force / (0.5 * h * (right - left)) - 1.0).abs() < 1.0e-12);
        let exact_total_heat = expected_mass * rt * (left / right).ln();
        heat_errors.push((step.receipt.viscous_heat_w_per_m / exact_total_heat - 1.0).abs());
        assert!(step.receipt.wall_slip_dissipation_w_per_m > 0.0);
        assert!(step.receipt.wall_slip_dissipation_w_per_m < step.receipt.viscous_heat_w_per_m);
        assert_mechanical_balance(&input, &step);
    }
    assert!(heat_errors[1] < 2.0e-5, "{heat_errors:?}");
    assert!(heat_errors[1] < heat_errors[0] / 10.0, "{heat_errors:?}");
}

#[test]
fn maxwell_squeeze_wedge_and_boundaries_preserve_mass_and_mechanical_power() {
    for boundary in [
        GasFilmBoundaryTopology::Sealed,
        GasFilmBoundaryTopology::Open {
            left_absolute_pressure_pa: 20_000.0,
            right_absolute_pressure_pa: 15_000.0,
        },
        GasFilmBoundaryTopology::Vented {
            cell_index: 5,
            absolute_pressure_pa: 16_000.0,
        },
    ] {
        let mut input = maxwell_fixture(12, 0.7);
        input.boundary = boundary;
        input.grid.gap_m = (0..12).map(|cell| 10.0e-6 + cell as f64 * 0.1e-6).collect();
        input.wall_motion.lower_tangential_velocity_m_per_s = -0.5;
        input.wall_motion.upper_tangential_velocity_m_per_s = 1.5;
        let rates = (0..12)
            .map(|cell| -1.0e-4 - cell as f64 * 1.0e-5)
            .collect::<Vec<_>>();
        input.wall_motion.gap_rate_profile_m_per_s = Some(rates.clone());
        let first = solve_isothermal_gas_film_1d(&input, None).expect("slipping wedge step");
        assert_mechanical_balance(&input, &first);
        assert!(first.receipt.mass_closure_residual_kg_per_m_s.abs() < 2.0e-11);
        for (gap, rate) in input.grid.gap_m.iter_mut().zip(rates) {
            *gap += rate * input.timestep_s;
        }
        let second = solve_isothermal_gas_film_1d(&input, Some(&first.checkpoint))
            .expect("slipping wedge restart");
        assert_mechanical_balance(&input, &second);
        assert!(second.receipt.mass_closure_residual_kg_per_m_s.abs() < 2.0e-11);
        assert!(second.receipt.wall_slip_dissipation_w_per_m > 0.0);
        assert_eq!(second.checkpoint.step_index, 2);
    }
}

#[test]
fn maxwell_admission_requires_physical_accommodation_and_first_order_knudsen() {
    for sigma in [0.0, -0.1, 1.01, f64::NAN, f64::INFINITY] {
        assert_eq!(
            solve_isothermal_gas_film_1d(&maxwell_fixture(4, sigma), None),
            Err(GasFilmError::InvalidInput {
                field: "slip_policy.tangential_momentum_accommodation"
            })
        );
    }
    let mut input = maxwell_fixture(4, 1.0);
    input.applicability.maximum_knudsen_number = 1.0;
    input.boundary = GasFilmBoundaryTopology::Open {
        left_absolute_pressure_pa: 20_000.0,
        right_absolute_pressure_pa: 5_000.0,
    };
    assert_eq!(
        solve_isothermal_gas_film_1d(&input, None),
        Err(GasFilmError::Unavailable {
            reason: "maxwell-slip-knudsen-outside-first-order-envelope"
        })
    );
    input.boundary = GasFilmBoundaryTopology::Sealed;
    input.applicability.mean_free_path_m = 0.11 * input.grid.gap_m[0];
    assert_eq!(
        solve_isothermal_gas_film_1d(&input, None),
        Err(GasFilmError::Unavailable {
            reason: "maxwell-slip-knudsen-outside-first-order-envelope"
        })
    );
}

#[test]
fn maxwell_checkpoint_binds_wall_law_and_replays_without_mutation() {
    let input = maxwell_fixture(4, 0.8);
    let step = solve_isothermal_gas_film_1d(&input, None).expect("Maxwell state");
    let replay = solve_isothermal_gas_film_1d(&input, Some(&step.checkpoint)).expect("replay");
    assert_eq!(
        replay,
        solve_isothermal_gas_film_1d(&input, Some(&step.checkpoint)).unwrap()
    );
    for policy in [
        SlipPolicy::MaxwellFirstOrder {
            source_id: "synthetic-shared-wall-maxwell-v1".to_owned(),
            tangential_momentum_accommodation: 0.6,
        },
        SlipPolicy::NoSlipContinuum {
            source_id: "synthetic-shared-wall-maxwell-v1".to_owned(),
        },
    ] {
        let mut changed = input.clone();
        changed.slip_policy = policy;
        assert_eq!(
            solve_isothermal_gas_film_1d(&changed, Some(&step.checkpoint)),
            Err(GasFilmError::CheckpointMismatch {
                field: "configuration_fingerprint"
            })
        );
    }
}

#[test]
fn maxwell_mach_uses_slipping_gas_profile_and_off_center_extremum() {
    let sigma = 0.8;
    for direction in [1.0, -1.0] {
        let mut input = maxwell_fixture(2, sigma);
        input.boundary = GasFilmBoundaryTopology::Open {
            left_absolute_pressure_pa: if direction > 0.0 { 30_000.0 } else { 20_000.0 },
            right_absolute_pressure_pa: if direction > 0.0 { 20_000.0 } else { 30_000.0 },
        };
        input.wall_motion.lower_tangential_velocity_m_per_s = direction;
        input.wall_motion.upper_tangential_velocity_m_per_s = 3.0 * direction;
        let step =
            solve_isothermal_gas_film_1d(&input, None).expect("Maxwell pressure/Couette flow");
        let inlet = if direction > 0.0 { 0 } else { 1 };
        let pressure = step.absolute_pressure_pa[inlet].unwrap();
        let face_pressure = 0.5 * (30_000.0 + pressure);
        let h = input.grid.gap_m[inlet];
        let beta = (2.0 - sigma) / sigma
            * input.applicability.mean_free_path_m
            * input.initial_absolute_pressure_pa
            / (face_pressure * h);
        let gradient = (30_000.0 - pressure) / (0.25 * input.grid.length_m);
        let curvature = h.powi(2) * gradient / (2.0 * input.gas.dynamic_viscosity_pa_s);
        let gas_velocity_difference = 2.0 / (1.0 + 2.0 * beta);
        let midpoint_speed = 2.0 + curvature * (0.25 + beta);
        let peak_speed = midpoint_speed + gas_velocity_difference.powi(2) / (4.0 * curvature);
        assert!(curvature > gas_velocity_difference);
        assert!(
            (step.receipt.maximum_gas_mach_number * input.applicability.speed_of_sound_m_per_s
                / peak_speed
                - 1.0)
                .abs()
                < 1.0e-12
        );
        input.applicability.maximum_mach_number =
            (midpoint_speed + peak_speed) / (2.0 * input.applicability.speed_of_sound_m_per_s);
        assert_eq!(
            solve_isothermal_gas_film_1d(&input, None),
            Err(GasFilmError::Unavailable {
                reason: "gas-mach-outside-declared-envelope"
            })
        );
    }
}
