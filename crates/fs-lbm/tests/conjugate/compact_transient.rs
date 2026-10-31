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
