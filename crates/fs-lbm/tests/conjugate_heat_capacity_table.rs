//! G1/G3: ordinary solid storage integrates the authored volumetric C(T).
//! The independent energy polynomial below crosses original source knots;
//! endpoint-frozen capacity does not satisfy these balances.

use fs_exec::CancelGate;
use fs_lbm::Face3;
use fs_lbm::conjugate::{
    ChtError, CompactComponent, EnergyConfig, FaceContactResistance, FlowField, FluidProperties,
    FvBoundary, RadiationConfig, STEFAN_BOLTZMANN, SimpleConfig, SolidMaterial, ThermalFace,
    ThermalSetup, TimeScheme, TransientConfig, UnsteadyConfig, Voxel, VoxelDomain, march_conjugate,
    march_energy,
};

const CAPACITY: [[f64; 2]; 4] = [[280.0, 1.0], [300.0, 2.0], [306.0, 5.0], [320.0, 33.0]];

fn solid() -> SolidMaterial {
    SolidMaterial::new("manufactured variable-capacity solid", 1.0)
        .with_heat_capacity_table(&CAPACITY)
        .unwrap()
}

fn fluid() -> FluidProperties {
    FluidProperties {
        density_kg_m3: 1.0,
        specific_heat_j_kg_k: 1.0,
        conductivity_w_m_k: 1.0,
        kinematic_viscosity_m2_s: 1.0,
    }
}

fn setup() -> ThermalSetup {
    let mut setup = ThermalSetup::new([ThermalFace::Adiabatic; 6]);
    setup.compact_storage_controls.maximum_iterations = 200;
    setup.compact_storage_controls.temperature_tolerance_k = 1e-10;
    setup.compact_storage_controls.energy_absolute_tolerance_j = 1e-9;
    setup.compact_storage_controls.energy_relative_tolerance = 1e-12;
    setup
}

fn controls(dt: f64, steps: usize) -> TransientConfig {
    TransientConfig {
        time_step_s: dt,
        steps,
        energy: EnergyConfig {
            tolerance: 1e-13,
            ..EnergyConfig::default()
        },
    }
}

// Independent integral of the manufactured Cvol(T), relative to 300 K.
fn energy_density(t: f64) -> f64 {
    let rise = t - 300.0;
    if rise < 0.0 {
        2.0 * rise + 0.025 * rise * rise
    } else if rise <= 6.0 {
        2.0 * rise + 0.25 * rise * rise
    } else {
        let above = rise - 6.0;
        21.0 + 5.0 * above + above * above
    }
}

fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() < 2e-8,
        "{actual:.15e} != {expected:.15e}"
    );
}

fn root(mut lo: f64, mut hi: f64, residual: impl Fn(f64) -> f64) -> f64 {
    for _ in 0..90 {
        let mid = 0.5 * (lo + hi);
        if residual(mid) > 0.0 {
            hi = mid;
        } else {
            lo = mid;
        }
    }
    0.5 * (lo + hi)
}

#[test]
fn isolated_heat_and_cooling_cross_original_capacity_knots() {
    let domain = VoxelDomain::from_fn(1, 1, 1, 1.0, |_| Voxel::Solid(0)).unwrap();
    let flow = FlowField::quiescent(&domain);
    // Reflecting the same source about 300 K also exercises steeply declining
    // positive capacity, where unrestricted fixed-point iterations can diverge.
    for reflected in [false, true] {
        let points: Vec<_> = if reflected {
            CAPACITY
                .iter()
                .rev()
                .map(|[t, c]| [600.0 - t, *c])
                .collect()
        } else {
            CAPACITY.to_vec()
        };
        let solids = [SolidMaterial::new("manufactured source", 1.0)
            .with_heat_capacity_table(&points)
            .unwrap()];
        let energy = |t| {
            if reflected {
                -energy_density(600.0 - t)
            } else {
                energy_density(t)
            }
        };
        let mut initial = 300.0;
        for unreflected in [312.0, 296.0, 306.0] {
            let target = if reflected {
                600.0 - unreflected
            } else {
                unreflected
            };
            let heat = energy(target) - energy(initial);
            let mut setup = setup();
            setup.power_w = vec![heat / 0.5];
            let result = march_energy(
                &domain,
                &fluid(),
                &solids,
                &flow,
                &setup,
                &[initial],
                |_| 1.0,
                &controls(0.5, 1),
                None,
                &CancelGate::new(),
            )
            .unwrap();
            close(result.temperature[0], target);
            close(result.records[0].stored_energy_change_j, heat);
            close(result.records[0].source_j, heat);
            close(result.records[0].boundary_outflow_j, 0.0);
            close(result.records[0].closure_j, 0.0);
            initial = result.temperature[0];
        }
    }
}

#[test]
fn graded_cell_volume_scales_energy_without_changing_temperature_response() {
    let domain = VoxelDomain::graded_from_fn([vec![0.1, 0.2, 0.7], vec![0.5], vec![2.0]], |_| {
        Voxel::Solid(0)
    })
    .unwrap();
    let mut setup = setup();
    let target = 312.0;
    let energy = energy_density(target);
    setup.power_w = (0..domain.cell_count())
        .map(|cell| energy * domain.volume(cell) / 0.25)
        .collect();
    let result = march_energy(
        &domain,
        &fluid(),
        &[solid()],
        &FlowField::quiescent(&domain),
        &setup,
        &[300.0; 3],
        |_| 1.0,
        &controls(0.25, 1),
        None,
        &CancelGate::new(),
    )
    .unwrap();
    for temperature in result.temperature {
        close(temperature, target);
    }
    close(result.records[0].stored_energy_change_j, energy);
    close(result.records[0].source_j, energy);
    close(result.records[0].closure_j, 0.0);
}

#[test]
fn closed_unequal_capacity_contact_exchange_conserves_exact_sensible_energy() {
    let domain = VoxelDomain::graded_from_fn([vec![0.25, 0.75], vec![1.0], vec![1.0]], |p| {
        Voxel::Solid(u16::from(p[0] > 0.25))
    })
    .unwrap();
    let solids = [
        solid(),
        SolidMaterial::new("constant-capacity receiver", 1.0).with_heat_capacity(4.0),
    ];
    let mut setup = setup();
    setup.face_contacts.push(FaceContactResistance {
        cells: [0, 1],
        resistance_m2_k_w: 0.4,
    });
    let result = march_energy(
        &domain,
        &fluid(),
        &solids,
        &FlowField::quiescent(&domain),
        &setup,
        &[312.0, 300.0],
        |_| 1.0,
        &controls(1.0, 16),
        None,
        &CancelGate::new(),
    )
    .unwrap();
    // G=A/(dx_hot/(2 k_hot)+R_contact+dx_cold/(2 k_cold)).
    let g = 1.0 / (0.125 + 0.4 + 0.375);
    let (mut hot, mut cold) = (312.0, 300.0);
    for record in &result.records {
        let old_hot = hot;
        let old_cold = cold;
        hot = root(cold, hot, |t| {
            let receiver = (3.0 * old_cold + g * t) / (3.0 + g);
            0.25 * (energy_density(t) - energy_density(old_hot)) + g * (t - receiver)
        });
        cold = (3.0 * old_cold + g * hot) / (3.0 + g);
        close(record.max_temperature_k, hot);
        close(record.stored_energy_change_j, 0.0);
        close(record.source_j, 0.0);
        close(record.boundary_outflow_j, 0.0);
        close(record.closure_j, 0.0);
    }
    assert!(hot < 306.0, "the cooling path must cross the Cp knot");
    close(result.temperature[0], hot);
    close(result.temperature[1], cold);
    close(
        0.25 * (energy_density(hot) - energy_density(312.0)) + 3.0 * (cold - 300.0),
        0.0,
    );
}

#[test]
fn coupled_momentum_schemes_reuse_the_same_nonlinear_storage_endpoint() {
    let domain = VoxelDomain::from_fn(2, 1, 1, 1.0, |p| {
        if p[0] < 1.0 {
            Voxel::Fluid
        } else {
            Voxel::Solid(0)
        }
    })
    .unwrap();
    let mut setup = setup();
    setup.fixed_temperature.push((0, 300.0));
    setup.power_w = vec![0.0, 4.0];
    let solids = [solid()];
    let schedule = |time: f64| if time <= 0.5 { 1.0 } else { 0.0 };
    let config = controls(0.5, 12);
    let frozen = march_energy(
        &domain,
        &fluid(),
        &solids,
        &FlowField::quiescent(&domain),
        &setup,
        &[300.0, 312.0],
        schedule,
        &config,
        None,
        &CancelGate::new(),
    )
    .unwrap();
    for scheme in [TimeScheme::BackwardEuler, TimeScheme::Bdf2] {
        let mut clock = UnsteadyConfig::new(config.time_step_s, config.steps);
        clock.scheme = scheme;
        clock.average_from_step = 1;
        let result = march_conjugate(
            &domain,
            &fluid(),
            &solids,
            &setup,
            &SimpleConfig::new([FvBoundary::wall(); 6]),
            &clock,
            None,
            &[300.0, 312.0],
            |_| 1.0,
            schedule,
            &config.energy,
            None,
            &CancelGate::new(),
        )
        .unwrap();
        for (actual, expected) in result.temperature.iter().zip(&frozen.temperature) {
            close(*actual, *expected);
        }
        let mut old = 312.0;
        for (step, record) in result.energy_records.iter().enumerate() {
            let power = 4.0 * schedule(record.time_s);
            let expected = root(300.0, 320.0, |t| {
                energy_density(t) - energy_density(old) + 0.5 * (t - 300.0 - power)
            });
            close(record.max_temperature_k, expected);
            close(
                record.stored_energy_change_j,
                energy_density(expected) - energy_density(old),
            );
            close(record.boundary_outflow_j, 0.5 * (expected - 300.0));
            close(record.closure_j, 0.0);
            close(
                record.stored_energy_change_j,
                frozen.records[step].stored_energy_change_j,
            );
            assert_eq!(result.flow.records[step].kinetic_energy_j, 0.0);
            old = expected;
        }
    }
}

#[test]
fn fixed_temperature_jump_books_the_integrated_capacity_as_thermostat_heat() {
    let domain = VoxelDomain::from_fn(1, 1, 1, 1.0, |_| Voxel::Solid(0)).unwrap();
    let mut setup = setup();
    setup.fixed_temperature.push((0, 312.0));
    let result = march_energy(
        &domain,
        &fluid(),
        &[solid()],
        &FlowField::quiescent(&domain),
        &setup,
        &[300.0],
        |_| 1.0,
        &controls(0.5, 2),
        None,
        &CancelGate::new(),
    )
    .unwrap();
    close(result.temperature[0], 312.0);
    close(result.records[0].stored_energy_change_j, 87.0);
    close(result.records[0].source_j, 0.0);
    close(result.records[0].boundary_outflow_j, -87.0);
    close(result.records[0].closure_j, 0.0);
    close(result.records[1].stored_energy_change_j, 0.0);
    close(result.records[1].boundary_outflow_j, 0.0);
}

#[test]
fn collapsed_compact_voxels_do_not_acquire_a_second_capacity() {
    let domain = VoxelDomain::from_fn(1, 1, 3, 1.0, |_| Voxel::Solid(0)).unwrap();
    let mut setup = setup();
    setup.fixed_temperature = vec![(0, 300.0), (2, 300.0)];
    setup.compact_components.push(CompactComponent {
        lo: [0, 0, 1],
        hi: [1, 1, 2],
        board_face: Face3::ZMin,
        power_w: 4.0,
        junction_to_case_k_w: 1.0,
        junction_to_board_k_w: 1.0,
        heat_capacity_j_k: Some(6.0),
    });
    let result = march_energy(
        &domain,
        &fluid(),
        &[solid()],
        &FlowField::quiescent(&domain),
        &setup,
        &[300.0; 3],
        |_| 1.0,
        &controls(1.0, 1),
        None,
        &CancelGate::new(),
    )
    .unwrap();
    let rise = 4.0 / (6.0 + 4.0 / 3.0);
    close(result.temperature[1], 300.0 + rise);
    close(result.records[0].stored_energy_change_j, 6.0 * rise);
    close(result.records[0].boundary_outflow_j, 4.0 / 3.0 * rise);
    close(result.records[0].closure_j, 0.0);
}

#[test]
fn radiation_and_capacity_integrals_converge_at_one_physical_endpoint() {
    let domain = VoxelDomain::from_fn(1, 1, 2, 1.0, |p| {
        if p[2] < 1.0 {
            Voxel::Solid(0)
        } else {
            Voxel::Fluid
        }
    })
    .unwrap();
    let fluid = FluidProperties {
        conductivity_w_m_k: 1e-12,
        ..fluid()
    };
    let mut setup = setup();
    setup.fixed_temperature.push((1, 300.0));
    let mut surroundings = [Some(300.0); 6];
    surroundings[4] = None;
    let mut radiation = RadiationConfig::new(vec![0.8], surroundings);
    radiation.surface_exchange = false;
    radiation.max_iterations = 200;
    let dt = 2.0;
    let coefficient = 0.8 * STEFAN_BOLTZMANN;
    let conductive_g = 1.0 / (0.5 + 0.5 / fluid.conductivity_w_m_k);
    let expected = root(300.0, 312.0, |t| {
        energy_density(t) - energy_density(312.0)
            + dt * (coefficient * (t.powi(4) - 300.0_f64.powi(4)) + conductive_g * (t - 300.0))
    });
    let result = march_energy(
        &domain,
        &fluid,
        &[solid()],
        &FlowField::quiescent(&domain),
        &setup,
        &[312.0, 300.0],
        |_| 1.0,
        &controls(dt, 1),
        Some(&radiation),
        &CancelGate::new(),
    )
    .unwrap();
    close(result.temperature[0], expected);
    let record = result.records[0];
    close(
        record.stored_energy_change_j,
        energy_density(expected) - energy_density(312.0),
    );
    close(
        record.radiated_w,
        coefficient * (expected.powi(4) - 300.0_f64.powi(4)),
    );
    close(record.boundary_outflow_j, -record.stored_energy_change_j);
    close(record.closure_j, 0.0);
}

#[test]
fn malformed_competing_and_out_of_domain_storage_refuse() {
    for points in [
        vec![],
        vec![[300.0, 2.0]],
        vec![[300.0, 2.0], [300.0, 3.0]],
        vec![[310.0, 2.0], [300.0, 3.0]],
        vec![[300.0, 0.0], [310.0, 3.0]],
        vec![[300.0, -2.0], [310.0, 3.0]],
        vec![[0.0, 2.0], [310.0, 3.0]],
        vec![[300.0, f64::NAN], [310.0, 3.0]],
        vec![[300.0, 2.0], [f64::INFINITY, 3.0]],
        (0..4097).map(|i| [300.0 + i as f64, 2.0]).collect(),
    ] {
        assert!(
            SolidMaterial::new("invalid", 1.0)
                .with_heat_capacity_table(&points)
                .is_err()
        );
    }
    assert!(
        SolidMaterial::new("competing", 1.0)
            .with_heat_capacity(2.0)
            .with_heat_capacity_table(&CAPACITY)
            .is_err()
    );
    let domain = VoxelDomain::from_fn(1, 1, 1, 1.0, |_| Voxel::Solid(0)).unwrap();
    let flow = FlowField::quiescent(&domain);
    let mut setup = setup();
    let run = |solid: SolidMaterial, setup: &ThermalSetup, initial: f64| {
        march_energy(
            &domain,
            &fluid(),
            &[solid],
            &flow,
            setup,
            &[initial],
            |_| 1.0,
            &controls(1.0, 1),
            None,
            &CancelGate::new(),
        )
    };
    assert!(matches!(
        run(solid().with_heat_capacity(2.0), &setup, 300.0),
        Err(ChtError::InvalidInput { .. })
    ));
    for initial in [279.0, 321.0] {
        assert!(run(solid(), &setup, initial).is_err());
    }
    setup.power_w = vec![1e6];
    assert!(
        run(solid(), &setup, 300.0).is_err(),
        "unsupported heating must not clamp to 320 K"
    );
    setup.power_w = vec![87.0];
    setup.compact_storage_controls.maximum_iterations = 1;
    setup.compact_storage_controls.temperature_tolerance_k = 1e3;
    assert!(
        run(solid(), &setup, 300.0).is_err(),
        "a loose temperature gate cannot waive the energy equation"
    );
}
