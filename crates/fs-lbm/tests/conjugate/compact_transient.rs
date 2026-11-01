//! Physical junction storage in compact components, independent of voxel aliases.

use fs_exec::CancelGate;
use fs_lbm::Face3;
use fs_lbm::conjugate::{
    ChtError, CompactComponent, EnergyConfig, FlowField, FluidProperties, SolidMaterial,
    ThermalFace, ThermalSetup, TransientConfig, Voxel, VoxelDomain, march_energy,
};

fn unit_fluid() -> FluidProperties {
    FluidProperties {
        density_kg_m3: 1.0,
        specific_heat_j_kg_k: 1.0,
        conductivity_w_m_k: 1.0,
        kinematic_viscosity_m2_s: 1.0,
    }
}

fn controls(dt: f64, steps: usize) -> TransientConfig {
    TransientConfig {
        time_step_s: dt,
        steps,
        energy: EnergyConfig::default(),
    }
}

fn pulse(time: f64) -> f64 {
    if time <= 1.5 { 1.0 } else { 0.0 }
}

fn reservoir_case() -> (VoxelDomain, [SolidMaterial; 2], ThermalSetup) {
    let domain =
        VoxelDomain::graded_from_fn([vec![0.25, 0.75], vec![1.0], vec![1.0, 0.2, 1.0]], |p| {
            Voxel::Solid(u16::from(p[2] > 1.0 && p[2] < 1.2))
        })
        .unwrap();
    let solids = [
        SolidMaterial::new("reservoir", 1.0).with_heat_capacity(2.0),
        // Neither package volume nor material capacity supplies junction storage.
        SolidMaterial::new("collapsed package", 1e-9),
    ];
    let mut setup = ThermalSetup::new([ThermalFace::Adiabatic; 6]);
    for z in [0, 2] {
        for x in 0..2 {
            setup.fixed_temperature.push((domain.index(x, 0, z), 300.0));
        }
    }
    setup.compact_components.push(CompactComponent {
        lo: [0, 0, 1],
        hi: [2, 1, 2],
        board_face: Face3::ZMin,
        power_w: 4.0,
        junction_to_case_k_w: 1.0,
        junction_to_board_k_w: 1.0,
        heat_capacity_j_k: Some(6.0),
    });
    (domain, solids, setup)
}

#[test]
fn compact_junction_matches_rc_pulse_and_counts_unequal_face_aliases_once() {
    let (domain, solids, setup) = reservoir_case();
    let dt = 0.75;
    let solution = march_energy(
        &domain,
        &unit_fluid(),
        &solids,
        &FlowField::quiescent(&domain),
        &setup,
        &vec![300.0; domain.cell_count()],
        pulse,
        &controls(dt, 6),
        None,
        &CancelGate::new(),
    )
    .unwrap();
    // Each whole branch has R=1+0.5=1.5 K/W. Unequal faces split its
    // conductance 1:3; the two alias cells still own only one 6 J/K capacity.
    let conductance = 4.0 / 3.0;
    let mut rise = 0.0;
    let mut outflow = 0.0;
    for record in &solution.records {
        let old = rise;
        let power = 4.0 * pulse(record.time_s);
        rise = (6.0 / dt * old + power) / (6.0 / dt + conductance);
        assert!((record.max_temperature_k - (300.0 + rise)).abs() < 1e-8);
        assert!((record.stored_energy_change_j - 6.0 * (rise - old)).abs() < 1e-8);
        assert!((record.source_j - power * dt).abs() < 1e-12);
        assert!((record.boundary_outflow_j - dt * conductance * rise).abs() < 1e-8);
        assert!(record.closure_j.abs() < 1e-8, "{record:?}");
        outflow += record.boundary_outflow_j;
    }
    let junction = solution.temperature[domain.index(0, 0, 1)];
    assert!((junction - (300.0 + rise)).abs() < 1e-8);
    assert_eq!(junction, solution.temperature[domain.index(1, 0, 1)]);
    assert!((6.0 * rise + outflow - 6.0).abs() < 1e-8);
}

#[test]
fn closed_compact_exchange_conserves_energy_and_converges_in_time() {
    let domain = VoxelDomain::from_fn(1, 1, 3, 1.0, |p| {
        Voxel::Solid(u16::from(p[2] > 1.0 && p[2] < 2.0))
    })
    .unwrap();
    let solids = [
        SolidMaterial::new("outside masses", 1.0).with_heat_capacity(2.0),
        SolidMaterial::new("junction", 1.0),
    ];
    let mut setup = ThermalSetup::new([ThermalFace::Adiabatic; 6]);
    setup.compact_components.push(CompactComponent {
        lo: [0, 0, 1],
        hi: [1, 1, 2],
        board_face: Face3::ZMin,
        power_w: 0.0,
        junction_to_case_k_w: 1.0,
        junction_to_board_k_w: 1.0,
        heat_capacity_j_k: Some(6.0),
    });
    // Symmetric outside temperatures reduce this to capacities 6 and 4 J/K,
    // joined by G=4/3 W/K. The weighted mean is 360 K and the difference
    // decays at lambda=G*(1/6+1/4)=5/9 per second.
    let initial = [300.0, 400.0, 300.0];
    let duration = 3.0;
    let lambda = 5.0 / 9.0;
    let exact = 360.0 + 40.0 * f64::exp(-lambda * duration);
    let mut errors = Vec::new();
    for steps in [1, 2, 4, 8, 16] {
        let dt = duration / steps as f64;
        let solution = march_energy(
            &domain,
            &unit_fluid(),
            &solids,
            &FlowField::quiescent(&domain),
            &setup,
            &initial,
            |_| 1.0,
            &controls(dt, steps),
            None,
            &CancelGate::new(),
        )
        .unwrap();
        let factor = (1.0 + lambda * dt).powi(-(steps as i32));
        assert!((solution.temperature[1] - (360.0 + 40.0 * factor)).abs() < 1e-8);
        for outside in [0, 2] {
            assert!((solution.temperature[outside] - (360.0 - 60.0 * factor)).abs() < 1e-8);
        }
        let total = 2.0 * solution.temperature[0]
            + 6.0 * solution.temperature[1]
            + 2.0 * solution.temperature[2];
        assert!((total - 3600.0).abs() < 1e-7);
        for record in &solution.records {
            assert!(record.stored_energy_change_j.abs() < 1e-8);
            assert_eq!(record.source_j, 0.0);
            assert_eq!(record.boundary_outflow_j, 0.0);
            assert!(record.closure_j.abs() < 1e-8, "{record:?}");
        }
        errors.push((solution.temperature[1] - exact).abs());
    }
    assert!(errors.windows(2).all(|pair| pair[1] < pair[0]));
    let ratio = errors[3] / errors[4];
    assert!(
        (1.8..2.1).contains(&ratio),
        "BE refinement errors: {errors:?}"
    );
}

fn nonlinear_branch_heat(junction: f64, outside: f64, half_width: f64, area: f64, r: f64) -> f64 {
    // k(T)=1+0.01*(T-300), Q=(Tj-Tsurface)/R and
    // integral(Toutside..Tsurface,k dT)=half_width*Q/area.
    // Solve the resulting quadratic for Tsurface-Toutside independently.
    let beta = half_width / (area * r);
    let linear = 1.0 + 0.01 * (outside - 300.0) + beta;
    let drive = beta * (junction - outside);
    let drop = 2.0 * drive / (linear + (linear * linear + 0.02 * drive).sqrt());
    (junction - outside - drop) / r
}

#[test]
fn compact_storage_and_nonlinear_half_cells_share_the_implicit_endpoint() {
    let domain =
        VoxelDomain::graded_from_fn([vec![0.2, 0.4], vec![0.3], vec![0.5, 0.2, 0.8]], |p| {
            Voxel::Solid(u16::from(p[2] > 0.5 && p[2] < 0.7))
        })
        .unwrap();
    let solids = [
        SolidMaterial::new("nonlinear outside", 1.0)
            .with_heat_capacity(2.0)
            .with_conductivity_table(&[
                (250.0, 0.5),
                (290.0, 0.9),
                (330.0, 1.3),
                (370.0, 1.7),
                (450.0, 2.5),
            ]),
        SolidMaterial::new("collapsed package", 1e-9),
    ];
    let (capacity, dt, old) = (7.0, 3.0, 350.0);
    let area = (0.2 + 0.4) * 0.3;
    for target in [280.0, 400.0] {
        let board = nonlinear_branch_heat(target, 300.0, 0.25, area, 2.0);
        let case = nonlinear_branch_heat(target, 310.0, 0.4, area, 5.0);
        let stored = capacity * (target - old);
        let power = stored / dt + board + case;
        let mut setup = ThermalSetup::new([ThermalFace::Adiabatic; 6]);
        let mut initial = vec![old; domain.cell_count()];
        for (z, temperature) in [(0, 300.0), (2, 310.0)] {
            for x in 0..2 {
                let cell = domain.index(x, 0, z);
                initial[cell] = temperature;
                setup.fixed_temperature.push((cell, temperature));
            }
        }
        setup.compact_components.push(CompactComponent {
            lo: [0, 0, 1],
            hi: [2, 1, 2],
            board_face: Face3::ZMin,
            power_w: power,
            junction_to_case_k_w: 5.0,
            junction_to_board_k_w: 2.0,
            heat_capacity_j_k: Some(capacity),
        });
        let solution = march_energy(
            &domain,
            &unit_fluid(),
            &solids,
            &FlowField::quiescent(&domain),
            &setup,
            &initial,
            |_| 1.0,
            &controls(dt, 1),
            None,
            &CancelGate::new(),
        )
        .unwrap();
        for x in 0..2 {
            assert!((solution.temperature[domain.index(x, 0, 1)] - target).abs() < 1e-7);
        }
        let record = solution.records[0];
        assert!((record.stored_energy_change_j - stored).abs() < 1e-7);
        assert!((record.source_j - power * dt).abs() < 1e-10);
        assert!((record.boundary_outflow_j - (board + case) * dt).abs() < 1e-7);
        assert!(record.closure_j.abs() < 1e-7, "target={target}: {record:?}");
    }
}

#[test]
fn compact_transients_require_capacity_and_one_unambiguous_initial_junction_state() {
    let (domain, solids, mut setup) = reservoir_case();
    let flow = FlowField::quiescent(&domain);
    let initial = vec![300.0; domain.cell_count()];
    for capacity in [
        None,
        Some(0.0),
        Some(-1.0),
        Some(f64::NAN),
        Some(f64::INFINITY),
    ] {
        setup.compact_components[0].heat_capacity_j_k = capacity;
        let result = march_energy(
            &domain,
            &unit_fluid(),
            &solids,
            &flow,
            &setup,
            &initial,
            |_| 1.0,
            &controls(0.75, 1),
            None,
            &CancelGate::new(),
        );
        assert!(
            matches!(result, Err(ChtError::InvalidInput { .. })),
            "capacity={capacity:?}"
        );
    }
    setup.compact_components[0].heat_capacity_j_k = Some(6.0);
    let mut ambiguous = initial;
    ambiguous[domain.index(1, 0, 1)] = 301.0;
    let result = march_energy(
        &domain,
        &unit_fluid(),
        &solids,
        &flow,
        &setup,
        &ambiguous,
        |_| 1.0,
        &controls(0.75, 1),
        None,
        &CancelGate::new(),
    );
    assert!(matches!(result, Err(ChtError::InvalidInput { .. })));
}

#[test]
fn compact_pulse_restart_reproduces_the_march_after_a_cancelled_attempt() {
    let (domain, solids, setup) = reservoir_case();
    let flow = FlowField::quiescent(&domain);
    let initial = vec![300.0; domain.cell_count()];
    let cancelled = CancelGate::new();
    let interrupted = march_energy(
        &domain,
        &unit_fluid(),
        &solids,
        &flow,
        &setup,
        &initial,
        |time| {
            // Cancel during the third step, after two complete physical steps.
            if time >= 2.25 {
                cancelled.request();
            }
            pulse(time)
        },
        &controls(0.75, 6),
        None,
        &cancelled,
    );
    assert!(matches!(interrupted, Err(ChtError::Cancelled)));
    let whole = march_energy(
        &domain,
        &unit_fluid(),
        &solids,
        &flow,
        &setup,
        &initial,
        pulse,
        &controls(0.75, 6),
        None,
        &CancelGate::new(),
    )
    .unwrap();
    let first = march_energy(
        &domain,
        &unit_fluid(),
        &solids,
        &flow,
        &setup,
        &initial,
        pulse,
        &controls(0.75, 2),
        None,
        &CancelGate::new(),
    )
    .unwrap();
    assert_eq!(
        first.temperature[domain.index(0, 0, 1)],
        first.temperature[domain.index(1, 0, 1)]
    );
    let rest = march_energy(
        &domain,
        &unit_fluid(),
        &solids,
        &flow,
        &setup,
        &first.temperature,
        |time| pulse(time + 1.5),
        &controls(0.75, 4),
        None,
        &CancelGate::new(),
    )
    .unwrap();
    for (a, b) in whole.temperature.iter().zip(&rest.temperature) {
        assert!((a - b).abs() < 1e-10);
    }
    for (a, b) in whole.records[2..].iter().zip(&rest.records) {
        assert!((a.time_s - (b.time_s + 1.5)).abs() < 1e-12);
        assert!((a.stored_energy_change_j - b.stored_energy_change_j).abs() < 1e-10);
        assert!((a.boundary_outflow_j - b.boundary_outflow_j).abs() < 1e-10);
        assert_eq!(a.source_j, b.source_j);
        assert!(b.closure_j.abs() < 1e-8, "{b:?}");
    }
}

fn bind_tabulated_junction(setup: &mut ThermalSetup) {
    use fs_lbm::conjugate::{CompactHeatStorage, CompactHeatStorageLaw};
    use fs_material::thermal::transport::storage::TabulatedHeatStorage;

    setup.compact_components[0].heat_capacity_j_k = None;
    setup.compact_heat_storage = vec![CompactHeatStorage {
        component: 0,
        law: CompactHeatStorageLaw::Tabulated(
            TabulatedHeatStorage::try_new(&[
                [290.0, 1.0], [300.0, 2.0], [306.0, 5.0], [320.0, 33.0],
            ]).unwrap(),
        ),
    }];
    let controls = &mut setup.compact_storage_controls;
    controls.maximum_iterations = 100;
    controls.temperature_tolerance_k = 1e-10;
    controls.energy_absolute_tolerance_j = 1e-9;
    controls.energy_relative_tolerance = 1e-12;
    controls.temperature_relaxation = 0.5;
}

fn pulse_junction_energy(rise: f64) -> f64 {
    // Independent integral from 300 K; all tested endpoints are >=300 K.
    if rise <= 6.0 {
        2.0 * rise + 0.25 * rise * rise
    } else {
        let above = rise - 6.0;
        21.0 + 5.0 * above + above * above
    }
}

#[test]
fn nonlinear_compact_pulse_crosses_knots_both_ways_and_conserves_closed_energy() {
    let (domain, solids, mut setup) = reservoir_case();
    setup.fixed_temperature.clear();
    setup.compact_components[0].power_w = 40.0;
    bind_tabulated_junction(&mut setup);
    let initial = vec![300.0; domain.cell_count()];
    let flow = FlowField::quiescent(&domain);
    let mut config = controls(1.0, 8);
    config.energy.tolerance = 1e-13;
    let run = march_energy(
        &domain, &unit_fluid(), &solids, &flow, &setup, &initial,
        |time| if time <= 1.0 { 1.0 } else { 0.0 },
        &config, None, &CancelGate::new(),
    ).unwrap();

    // Each outside layer owns 2 J/K and connects through G=2/3 W/K.
    // Eliminate its BE endpoint, then solve the remaining scalar equation
    // using the independent piecewise energy integral above.
    let g = 2.0 / 3.0;
    let (mut outside, mut junction) = (0.0_f64, 0.0_f64);
    for (step, record) in run.records.iter().enumerate() {
        let power = if step == 0 { 40.0 } else { 0.0 };
        let (old_outside, old_junction) = (outside, junction);
        let (mut lo, mut hi) = (0.0, 20.0);
        for _ in 0..90 {
            let j = 0.5 * (lo + hi);
            let o = (2.0 * old_outside + g * j) / (2.0 + g);
            let residual = pulse_junction_energy(j) - pulse_junction_energy(old_junction)
                + 2.0 * g * (j - o) - power;
            if residual > 0.0 { hi = j; } else { lo = j; }
        }
        junction = 0.5 * (lo + hi);
        outside = (2.0 * old_outside + g * junction) / (2.0 + g);
        assert!((record.max_temperature_k - 300.0 - junction).abs() < 1e-8);
        let actual_j = record.max_temperature_k - 300.0;
        let actual_o = (2.0 * old_outside + g * actual_j) / (2.0 + g);
        let local = pulse_junction_energy(actual_j) - pulse_junction_energy(old_junction)
            + 2.0 * g * (actual_j - actual_o) - power;
        assert!(local.abs() < 2e-8, "step={step}, local={local}");
        assert!((pulse_junction_energy(junction) + 4.0 * outside - 40.0).abs() < 1e-10);
        assert!((record.stored_energy_change_j - power).abs() < 2e-8);
        assert!((record.source_j - power).abs() < 1e-12);
        assert_eq!(record.boundary_outflow_j, 0.0);
        assert!(record.closure_j.abs() < 2e-8, "{record:?}");
    }
    assert!(run.records[0].max_temperature_k > 306.0);
    assert!(run.records[4].max_temperature_k < 306.0);
    for (cell, &temperature) in run.temperature.iter().enumerate() {
        let expected = if domain.coords(cell)[2] == 1 { junction } else { outside };
        assert!((temperature - 300.0 - expected).abs() < 1e-8);
    }
    assert_eq!(run.temperature[domain.index(0, 0, 1)], run.temperature[domain.index(1, 0, 1)]);

    // The public resolved-field step shares the same nonlinear owner.
    // Both compact aliases have zero independently supplied cell capacity.
    let capacities: Vec<_> = (0..domain.cell_count()).map(|cell| {
        if domain.coords(cell)[2] == 1 { 0.0 } else { 2.0 * domain.volume(cell) }
    }).collect();
    let (step, heat) = config.energy.step_with_conductivity(
        &domain, &unit_fluid(), &solids, &flow, &setup,
        &vec![[1.0; 3]; domain.cell_count()], &capacities, &initial, 1.0, &CancelGate::new(),
    ).unwrap();
    let first_junction = 3.0 + 22.0_f64.sqrt();
    for (cell, &temperature) in step.temperature.iter().enumerate() {
        let expected = if domain.coords(cell)[2] == 1 { first_junction } else { first_junction / 4.0 };
        assert!((temperature - 300.0 - expected).abs() < 1e-8);
    }
    assert!((heat.stored_energy_change_j - 40.0).abs() < 2e-8);
    assert!(heat.closure_j.abs() < 2e-8);
}

#[test]
fn nonlinear_compact_storage_refuses_overlap_extrapolation_and_unclosed_steps() {
    let (domain, solids, mut setup) = reservoir_case();
    bind_tabulated_junction(&mut setup);
    let initial = vec![300.0; domain.cell_count()];
    let flow = FlowField::quiescent(&domain);
    let mut config = controls(1.0, 1);
    config.energy.tolerance = 1e-13;
    let run = |setup: &ThermalSetup, previous: &[f64]| march_energy(
        &domain, &unit_fluid(), &solids, &flow, setup, previous,
        |_| 1.0, &config, None, &CancelGate::new(),
    );
    setup.compact_components[0].heat_capacity_j_k = Some(6.0);
    assert!(matches!(run(&setup, &initial), Err(ChtError::InvalidInput { .. })));
    setup.compact_components[0].heat_capacity_j_k = None;
    let duplicate = setup.compact_heat_storage[0].clone();
    setup.compact_heat_storage.push(duplicate);
    assert!(run(&setup, &initial).is_err());
    setup.compact_heat_storage.pop();
    setup.compact_heat_storage[0].component = 1;
    assert!(run(&setup, &initial).is_err());
    setup.compact_heat_storage[0].component = 0;
    let mut unsupported_initial = initial.clone();
    for x in 0..2 { unsupported_initial[domain.index(x, 0, 1)] = 289.0; }
    assert!(run(&setup, &unsupported_initial).is_err());

    // Even at the upper endpoint, storage plus branch outflow cannot absorb
    // 1000 J. A clamped 320 K endpoint must not be returned as success.
    setup.compact_components[0].power_w = 1000.0;
    assert!(run(&setup, &initial).is_err());
    setup.compact_components[0].power_w = 40.0;
    setup.compact_storage_controls.maximum_iterations = 1;
    setup.compact_storage_controls.temperature_tolerance_k = 1000.0;
    // A loose temperature gate cannot waive the local energy defect.
    assert!(run(&setup, &initial).is_err());
    setup.compact_storage_controls.maximum_iterations = 100;
    setup.compact_storage_controls.temperature_tolerance_k = 1e-10;
    let retry = run(&setup, &initial).unwrap();
    assert!(retry.records[0].closure_j.abs() < 2e-8);
    assert!((306.0..312.0).contains(&retry.records[0].max_temperature_k));
}

#[test]
fn compact_tabulated_storage_recovers_from_an_unsupported_linear_predictor() {
    use fs_lbm::conjugate::{CompactHeatStorage, CompactHeatStorageLaw};
    use fs_material::thermal::transport::storage::TabulatedHeatStorage;
    let knots = [[300.0, 1.0], [399.0, 1.0], [400.0, 1000.0]];
    let conductance = 4.0 / 3.0;
    let expected_temperature = 399.5;
    // 99 K at 1 J/K, then half of the steep final affine segment.
    let expected_storage = 99.0 + 0.5 * (1.0 + 500.5) * 0.5;
    let expected_outflow = conductance * (expected_temperature - 300.0);
    let power = expected_storage + expected_outflow;
    assert!(300.0 + power / (1.0 + conductance) > 400.0);
    for relaxation in [1.0, 0.1] {
        let (domain, solids, mut setup) = reservoir_case();
        setup.compact_components[0].heat_capacity_j_k = None;
        setup.compact_components[0].power_w = power;
        setup.compact_heat_storage = vec![CompactHeatStorage {
            component: 0,
            law: CompactHeatStorageLaw::Tabulated(TabulatedHeatStorage::try_new(&knots).unwrap()),
        }];
        setup.compact_storage_controls.maximum_iterations = 100;
        setup.compact_storage_controls.temperature_relaxation = relaxation;
        let solution = march_energy(
            &domain, &unit_fluid(), &solids, &FlowField::quiescent(&domain),
            &setup, &vec![300.0; domain.cell_count()], |_| 1.0,
            &controls(1.0, 1), None, &CancelGate::new(),
        ).unwrap();
        let junction = solution.temperature[domain.index(0, 0, 1)];
        assert!((junction - expected_temperature).abs() < 1e-8,
            "relaxation={relaxation}, junction={junction}");
        assert_eq!(junction, solution.temperature[domain.index(1, 0, 1)]);
        let above = junction - 399.0;
        let actual_storage = 99.0 + above + 499.5 * above * above;
        let actual_outflow = conductance * (junction - 300.0);
        assert!((actual_storage + actual_outflow - power).abs() < 1e-6);
        assert_eq!(solution.records.len(), 1);
        let record = solution.records[0];
        assert!((record.stored_energy_change_j - expected_storage).abs() < 1e-6, "{record:?}");
        assert!((record.boundary_outflow_j - expected_outflow).abs() < 1e-6, "{record:?}");
        assert!((record.source_j - power).abs() < 1e-10);
        assert!(record.closure_j.abs() < 1e-6, "{record:?}");
    }
}

#[test]
fn compact_tabulated_storage_resolves_each_junction_despite_global_cancellation() {
    use fs_lbm::conjugate::{CompactHeatStorage, CompactHeatStorageLaw};
    use fs_material::thermal::transport::storage::TabulatedHeatStorage;
    let (_, solids, mut setup) = reservoir_case();
    let template = setup.compact_components[0];
    let domain = VoxelDomain::from_fn(2, 1, 3, 1.0, |p| {
        Voxel::Solid(u16::from(p[2] > 1.0 && p[2] < 2.0))
    }).unwrap();
    setup.fixed_temperature.clear();
    setup.compact_components.clear();
    for x in 0..2 {
        for z in [0, 2] { setup.fixed_temperature.push((domain.index(x, 0, z), 300.0)); }
        let mut component = template;
        component.lo = [x, 0, 1];
        component.hi = [x + 1, 1, 2];
        component.power_w = if x == 0 { 0.001 } else { -0.001 };
        component.heat_capacity_j_k = None;
        setup.compact_components.push(component);
        setup.compact_heat_storage.push(CompactHeatStorage {
            component: x,
            law: CompactHeatStorageLaw::Tabulated(
                TabulatedHeatStorage::try_new(&[[299.0, 6.0], [301.0, 6.0]]).unwrap(),
            ),
        });
    }
    let mut config = controls(1.0, 1);
    // The unchanged warm start passes this scaled linear tolerance.
    // Per-junction physical gates must tighten the same inner solver.
    config.energy.tolerance = 1e-6;
    let solution = march_energy(
        &domain, &unit_fluid(), &solids, &FlowField::quiescent(&domain),
        &setup, &vec![300.0; domain.cell_count()], |_| 1.0,
        &config, None, &CancelGate::new(),
    ).unwrap();
    let capacity = 6.0;
    let conductance = 4.0 / 3.0;
    for x in 0..2 {
        let power = setup.compact_components[x].power_w;
        let rise = solution.temperature[domain.index(x, 0, 1)] - 300.0;
        let exact_rise = power / (capacity + conductance);
        assert!((rise - exact_rise).abs() < 1e-9, "component={x}, rise={rise}");
        let physical_residual = capacity * rise + conductance * rise - power;
        assert!(physical_residual.abs() < 1e-8, "component={x}, residual={physical_residual}");
    }
    let record = solution.records[0];
    assert_eq!(record.source_j, 0.0);
    assert!(record.stored_energy_change_j.abs() < 1e-8);
    assert!(record.boundary_outflow_j.abs() < 1e-8);
    assert!(record.closure_j.abs() < 1e-8, "{record:?}");
}

#[test]
fn nonlinear_compact_support_boundary_allows_other_junction_trials_to_advance() {
    use fs_lbm::conjugate::{CompactHeatStorage, CompactHeatStorageLaw};
    use fs_material::thermal::transport::storage::TabulatedHeatStorage;
    // An ordinary conductive cell separates the two compact boxes.
    const K: f64 = 1e12;
    const CONNECTOR_CAPACITY: f64 = 1e-9;
    const OUTER_RESISTANCE: f64 = 1e12;
    let domain = VoxelDomain::from_fn(1, 1, 5, 1.0, |p| {
        let package = (1.0..2.0).contains(&p[2]) || (3.0..4.0).contains(&p[2]);
        Voxel::Solid(u16::from(package))
    }).unwrap();
    let solids = [
        SolidMaterial::new("conductive surroundings", K).with_heat_capacity(CONNECTOR_CAPACITY),
        SolidMaterial::new("collapsed packages", 1.0),
    ];
    let mut setup = ThermalSetup::new([ThermalFace::Adiabatic; 6]);
    setup.fixed_temperature = vec![(domain.index(0, 0, 0), 300.0), (domain.index(0, 0, 4), 300.0)];
    setup.compact_components = vec![
        CompactComponent {
            lo: [0, 0, 1], hi: [1, 1, 2], board_face: Face3::ZMin, power_w: 30.0,
            junction_to_case_k_w: 0.5, junction_to_board_k_w: OUTER_RESISTANCE,
            heat_capacity_j_k: None,
        },
        CompactComponent {
            lo: [0, 0, 3], hi: [1, 1, 4], board_face: Face3::ZMin, power_w: 0.0,
            junction_to_case_k_w: OUTER_RESISTANCE, junction_to_board_k_w: 0.5,
            heat_capacity_j_k: None,
        },
    ];
    setup.compact_heat_storage = vec![
        CompactHeatStorage {
            component: 0,
            law: CompactHeatStorageLaw::Tabulated(
                TabulatedHeatStorage::try_new(&[[300.0, 1.0], [400.0, 1.0]]).unwrap(),
            ),
        },
        CompactHeatStorage {
            component: 1,
            law: CompactHeatStorageLaw::Tabulated(
                TabulatedHeatStorage::try_new(&[[350.0, 1.0], [400.0, 501.0]]).unwrap(),
            ),
        },
    ];
    setup.compact_storage_controls.temperature_tolerance_k = 1e-10;
    setup.compact_storage_controls.energy_relative_tolerance = 1e-12;
    let initial = vec![300.0, 400.0, 375.0, 350.0, 300.0];
    let mut config = controls(1.0, 1);
    config.energy.tolerance = 1e-13;

    // Eliminate the finite-capacity connector and first junction. The
    // remaining scalar energy equation has storage y+5*y² for junction 2.
    let g = 1.0 / (0.5 + 0.5 / K);
    let outer_g = 1.0 / (OUTER_RESISTANCE + 0.5 / K);
    let connector_diagonal = CONNECTOR_CAPACITY + 2.0 * g;
    let coupling = g * g / connector_diagonal;
    let first_diagonal = 1.0 + g + outer_g - coupling;
    let first_rhs = 400.0 + 30.0 + 300.0 * outer_g
        + g * CONNECTOR_CAPACITY * 375.0 / connector_diagonal;
    let first_from_second = |second: f64| (first_rhs + coupling * second) / first_diagonal;
    let connector_from_junctions = |first: f64, second: f64| {
        (CONNECTOR_CAPACITY * 375.0 + g * (first + second)) / connector_diagonal
    };
    let second_initial_rhs = 350.0 + 300.0 * outer_g
        + g * CONNECTOR_CAPACITY * 375.0 / connector_diagonal;
    let frozen_first = (first_rhs * first_diagonal + coupling * second_initial_rhs)
        / (first_diagonal * first_diagonal - coupling * coupling);
    assert!(frozen_first > 403.0);
    let (mut low, mut high) = (350.0_f64, 400.0_f64);
    for _ in 0..90 {
        let second = 0.5 * (low + high);
        let first = first_from_second(second);
        let connector = connector_from_junctions(first, second);
        let rise = second - 350.0;
        let residual = rise + 5.0 * rise * rise
            + g * (second - connector) + outer_g * (second - 300.0);
        if residual > 0.0 { high = second; } else { low = second; }
    }
    let expected_second = 0.5 * (low + high);
    let expected_first = first_from_second(expected_second);
    let expected_connector = connector_from_junctions(expected_first, expected_second);
    assert!((391.0..392.0).contains(&expected_first));
    assert!((352.0..353.0).contains(&expected_second));
    let run = march_energy(
        &domain, &unit_fluid(), &solids, &FlowField::quiescent(&domain),
        &setup, &initial, |_| 1.0, &config, None, &CancelGate::new(),
    ).unwrap();
    let first = run.temperature[domain.index(0, 0, 1)];
    let connector = run.temperature[domain.index(0, 0, 2)];
    let second = run.temperature[domain.index(0, 0, 3)];
    assert!((first - expected_first).abs() < 1e-8);
    assert!((connector - expected_connector).abs() < 1e-8);
    assert!((second - expected_second).abs() < 1e-8);
    assert!((300.0..=400.0).contains(&first));
    assert!((350.0..=400.0).contains(&second));
    assert!(first < initial[1]);
    let first_storage = first - 400.0;
    let second_rise = second - 350.0;
    let second_storage = second_rise + 5.0 * second_rise * second_rise;
    let first_defect = first_storage + g * (first - connector) + outer_g * (first - 300.0) - 30.0;
    let second_defect = second_storage + g * (second - connector) + outer_g * (second - 300.0);
    assert!(first_defect.abs() < 2e-8, "{first_defect}");
    assert!(second_defect.abs() < 2e-8, "{second_defect}");
    let total_storage = first_storage + second_storage + CONNECTOR_CAPACITY * (connector - 375.0);
    let reservoir_outflow = outer_g * ((first - 300.0) + (second - 300.0));
    let heat = &run.records[0];
    assert!((heat.stored_energy_change_j - total_storage).abs() < 2e-8);
    assert!((heat.boundary_outflow_j - reservoir_outflow).abs() < 2e-8);
    assert!((heat.source_j - 30.0).abs() < 1e-12);
    assert!(heat.closure_j.abs() < 2e-8);
}
