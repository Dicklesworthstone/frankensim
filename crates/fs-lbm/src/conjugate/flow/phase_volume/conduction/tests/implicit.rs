//! G1/G3/G4: independent stiff-grid, nonlinear-source and native-flow limits.

use super::*;
use crate::conjugate::flow::ThermalVolumeMomentumInput;
use crate::d3q19::Face3;

fn solve_controls() -> ImplicitThermalVolumeSolveControls {
    ImplicitThermalVolumeSolveControls {
        cell_tolerance_j: 1e-11,
        ..ImplicitThermalVolumeSolveControls::default()
    }
}

fn bath(
    cell: usize,
    face: Face3,
    temperature_k: f64,
    resistance_m2_k_w: f64,
) -> ThermalVolumeBathFace {
    ThermalVolumeBathFace {
        cell,
        face,
        temperature_k,
        resistance_m2_k_w,
    }
}

#[test]
fn stiff_cartesian_pair_cycle_and_cube_match_backward_euler_without_explicit_cfl() {
    let material = card("T", &[(300.0, 4.0), (500.0, 4.0)]);
    let caloric = curve(
        &material,
        &[(0.0, 300.0, 1.0, 1.0), (400.0, 500.0, 1.0, 1.0)],
    );
    let source = source(&material, &caloric, [0.0, 400.0]);
    let gate = CancelGate::new_clock_free();
    for axes in 1..=3 {
        let domain = VoxelDomain::new(
            2,
            if axes >= 2 { 2 } else { 1 },
            if axes == 3 { 2 } else { 1 },
            1.0,
        )
        .unwrap();
        let n = domain.cell_count();
        let field = FlowField::quiescent(&domain);
        let mass = vec![1.0; n];
        let old: Vec<_> = (0..n)
            .map(|cell| {
                if cell.count_ones() % 2 == 0 {
                    100.0
                } else {
                    300.0
                }
            })
            .collect();
        let states: Vec<_> = old
            .iter()
            .map(|&h| Some(caloric.state_at_specific_enthalpy(h).unwrap()))
            .collect();
        let zero = vec![0.0; n];
        let data = input(&mass, &states, &zero);
        let limits = controls(8.0);
        invalid(
            field.advance_isotropic_conductive_thermal_volume(
                &domain,
                [FlowFace::Wall; 6],
                &source,
                data,
                limits,
                &gate,
            ),
            "mixing weight exceeds one",
        );
        let mut workspace = ImplicitThermalVolumeWorkspace::new(n).unwrap();
        let result = field
            .advance_isotropic_implicit_conductive_thermal_volume_with_baths(
                &domain,
                [FlowFace::Wall; 6],
                &source,
                &[],
                data,
                limits,
                solve_controls(),
                &mut workspace,
                &gate,
            )
            .unwrap();
        // Every neighbor has opposite checkerboard sign. M=1, Cp=2, G=4,
        // so each of d incident links gives eigenvalue 4 and amplitude
        // a_BE=100/(1+4*d*dt), at explicit mixing weight 16*d.
        for cell in 0..n {
            let expected = 200.0 + (old[cell] - 200.0) / (1.0 + 32.0 * f64::from(axes));
            close(
                result.thermal.volume.transport.specific_enthalpy_j_kg[cell],
                expected,
            );
            close(
                result.conductive_enthalpy_increment_j[cell],
                expected - old[cell],
            );
            close(result.thermal.volume.transport.mass_kg[cell], 1.0);
        }
        assert_eq!(result.report.conductive_faces, axes as usize * n / 2);
        assert_eq!(result.report.solve.energy.substeps, 1);
        close(result.report.solve.energy.external_heat_j, 0.0);
        close(result.report.conduction_balance_defect_j, 0.0);
        assert!(result.report.endpoint_heat_defect_j <= 1e-9);
        let replay = field
            .advance_isotropic_implicit_conductive_thermal_volume_with_baths(
                &domain,
                [FlowFace::Wall; 6],
                &source,
                &[],
                data,
                limits,
                solve_controls(),
                &mut workspace,
                &gate,
            )
            .unwrap();
        assert_eq!(result, replay);
        assert_eq!(field, FlowField::quiescent(&domain));
    }
}

#[test]
fn decimal_graded_interior_keeps_exact_native_face_coverage() {
    let material = card("T", &[(300.0, 4.0), (500.0, 4.0)]);
    let caloric = curve(
        &material,
        &[(0.0, 300.0, 1.0, 1.0), (400.0, 500.0, 1.0, 1.0)],
    );
    let source = source(&material, &caloric, [0.0, 400.0]);
    let domain = VoxelDomain::graded([vec![0.1; 3], vec![0.2; 3], vec![0.7; 3]]).unwrap();
    let n = domain.cell_count();
    let mass: Vec<_> = (0..n).map(|cell| domain.volume(cell)).collect();
    let states = vec![Some(caloric.state_at_specific_enthalpy(100.0).unwrap()); n];
    let result = FlowField::quiescent(&domain)
        .advance_isotropic_implicit_conductive_thermal_volume_with_baths(
            &domain,
            [FlowFace::Wall; 6],
            &source,
            &[],
            input(&mass, &states, &vec![0.0; n]),
            controls(10.0),
            solve_controls(),
            &mut ImplicitThermalVolumeWorkspace::new(n).unwrap(),
            &CancelGate::new_clock_free(),
        )
        .unwrap();
    // For the central cell the native six-link sum is 0.46, whereas
    // 2*(Ax+Ay+Az) rounds to 0.45999999999999996. Keep actual occupied
    // area before adding exposed faces, without masking overlap with a tolerance.
    assert_eq!(result.report.conductive_faces, 54);
    assert_eq!(result.conductive_enthalpy_increment_j, vec![0.0; n]);
    assert_eq!(result.thermal.volume.transport.mass_kg, mass);
    assert_eq!(result.thermal.volume.phase_states, states);
}

#[test]
fn nonlinear_source_and_series_bath_use_the_same_heated_endpoint() {
    // T=300+h, k=1+h. At the exact h*=2 endpoint, k*=3 and
    // G=1/(.5/3+.5)=1.5; dt=.25 gives Q_bath=.75, authored Q=.25.
    let material = card("T", &[(300.0, 1.0), (304.0, 5.0)]);
    let caloric = curve(&material, &[(0.0, 300.0, 1.0, 1.0), (4.0, 304.0, 1.0, 1.0)]);
    let source = source(&material, &caloric, [0.0, 4.0]);
    let domain = VoxelDomain::new(1, 1, 1, 1.0).unwrap();
    let field = FlowField::quiescent(&domain);
    let states = [Some(caloric.state_at_specific_enthalpy(1.0).unwrap())];
    let mut workspace = ImplicitThermalVolumeWorkspace::new(1).unwrap();
    let result = field
        .advance_isotropic_implicit_conductive_thermal_volume_with_baths(
            &domain,
            [FlowFace::Wall; 6],
            &source,
            &[bath(0, Face3::XMin, 304.0, 0.5)],
            input(&[1.0], &states, &[0.25]),
            controls(0.25),
            solve_controls(),
            &mut workspace,
            &CancelGate::new_clock_free(),
        )
        .unwrap();
    close(
        result.thermal.volume.transport.specific_enthalpy_j_kg[0],
        2.0,
    );
    close(
        result.thermal.volume.phase_states[0]
            .unwrap()
            .temperature_k(),
        302.0,
    );
    close(result.conductive_enthalpy_increment_j[0], 0.75);
    close(result.report.boundary_heat_j[0], 0.75);
    close(result.report.solve.energy.external_heat_j, 1.0);
    close(
        result.thermal.volume.report.requested_enthalpy_increment_j,
        1.0,
    );
    assert_eq!(
        result.report.conductivity_source_identity,
        source.identity()
    );
    assert!(result.report.solve.entropy.has_unclassified_sources);
    assert_eq!(states[0].unwrap().specific_enthalpy_j_kg(), 1.0);
}

#[test]
fn phase_dependent_conductivity_is_implicit_along_a_true_latent_plateau() {
    let material = card("f", &[(0.0, 2.0), (1.0, 6.0)]);
    let caloric = curve(
        &material,
        &[
            (-100.0, 300.0, 0.0, 1.0),
            (0.0, 400.0, 0.0, 1.0),
            (100.0, 400.0, 1.0, 1.0),
            (200.0, 500.0, 1.0, 1.0),
        ],
    );
    let source = source(&material, &caloric, [0.0, 100.0]);
    let domain = VoxelDomain::new(1, 1, 1, 1.0).unwrap();
    let field = FlowField::quiescent(&domain);
    let states = [Some(caloric.state_at_specific_enthalpy(25.0).unwrap())];
    let mut workspace = ImplicitThermalVolumeWorkspace::new(1).unwrap();
    let result = field
        .advance_isotropic_implicit_conductive_thermal_volume_with_baths(
            &domain,
            [FlowFace::Wall; 6],
            &source,
            &[bath(0, Face3::XMax, 410.0, 0.0)],
            input(&[1.0], &states, &[0.0]),
            controls(0.5),
            solve_controls(),
            &mut workspace,
            &CancelGate::new_clock_free(),
        )
        .unwrap();
    // k=2+.04*h and T=400: h*=25+.5*(4+.08*h*)*10 gives h*=75.
    let phase = result.thermal.volume.phase_states[0].unwrap();
    close(phase.specific_enthalpy_j_kg(), 75.0);
    close(phase.temperature_k(), 400.0);
    close(phase.liquid_mass_fraction(), 0.75);
    close(result.conductive_enthalpy_increment_j[0], 50.0);
    close(
        source
            .state_at_specific_enthalpy(75.0)
            .unwrap()
            .conductivity_w_m_k(),
        5.0,
    );
    assert_eq!(result.report.solve.energy.substeps, 1);
}

#[test]
fn distinct_graded_wall_baths_keep_opposed_heat_and_zero_solid_storage() {
    let material = card("T", &[(300.0, 4.0), (500.0, 4.0)]);
    let caloric = curve(
        &material,
        &[(0.0, 300.0, 1.0, 1.0), (400.0, 500.0, 1.0, 1.0)],
    );
    let source = source(&material, &caloric, [0.0, 400.0]);
    let mut domain = VoxelDomain::graded([vec![2.0, 7.0], vec![3.0], vec![4.0]]).unwrap();
    domain.set(1, 0, 0, Voxel::Solid(0));
    let field = FlowField::quiescent(&domain);
    let states = [
        Some(caloric.state_at_specific_enthalpy(100.0).unwrap()),
        None,
    ];
    let baths = [
        bath(0, Face3::XMax, 360.0, 0.5),
        bath(0, Face3::YMin, 340.0, 0.125),
        bath(0, Face3::ZMin, 370.0, 0.0),
    ];
    let run = |baths: &[ThermalVolumeBathFace]| {
        field
            .advance_isotropic_implicit_conductive_thermal_volume_with_baths(
                &domain,
                [FlowFace::Wall; 6],
                &source,
                baths,
                input(&[24.0, 0.0], &states, &[0.0; 2]),
                controls(12.0),
                solve_controls(),
                &mut ImplicitThermalVolumeWorkspace::new(1).unwrap(),
                &CancelGate::new_clock_free(),
            )
            .unwrap()
    };
    let result = run(&baths);
    // Cp*M=48, G=[16,16,12]. At dt=12, T*=350+12*240/(48+12*44)=355.
    // Individual signed transfers are +960, -2880 and +2160 J, net +240 J.
    close(
        result.thermal.volume.phase_states[0]
            .unwrap()
            .temperature_k(),
        355.0,
    );
    close(
        result.thermal.volume.transport.specific_enthalpy_j_kg[0],
        110.0,
    );
    assert_eq!(result.thermal.volume.transport.mass_kg, [24.0, 0.0]);
    assert_eq!(result.thermal.volume.phase_states[1], None);
    close(result.conductive_enthalpy_increment_j[0], 240.0);
    assert_eq!(result.conductive_enthalpy_increment_j[1], 0.0);
    for (side, expected) in [(1, 960.0), (2, -2880.0), (4, 2160.0)] {
        close(result.report.boundary_heat_j[side], expected);
    }
    close(result.report.absolute_boundary_heat_j, 6000.0);
    close(result.report.solve.energy.external_heat_j, 240.0);
    assert_eq!(result.report.boundary_faces, 3);
    assert_eq!(result.report.conductive_faces, 0);
    assert_eq!(result.report, run(&[baths[2], baths[0], baths[1]]).report);
}

#[test]
fn implicit_heat_drives_original_density_flow_and_momentum_once() {
    // v=1+h/2, T=300+h, k=1. Initially M=.5, h=2. With dt=3/8,
    // h*=[2+3/8,2-3/8] solves the original two-cell/bath BE system,
    // supplying Q=[3/16,-3/16] and driving actual volume flow 1/4 m3/s.
    let material = card("T", &[(300.0, 1.0), (304.0, 1.0)]);
    let caloric = curve(
        &material,
        &[(0.0, 300.0, 1.0, 1.0), (4.0, 304.0, 1.0, 1.0 / 3.0)],
    );
    let source = source(&material, &caloric, [0.0, 4.0]);
    let domain = VoxelDomain::new(2, 1, 1, 1.0).unwrap();
    let field = FlowField::quiescent(&domain);
    let states = [Some(caloric.state_at_specific_enthalpy(2.0).unwrap()); 2];
    let zero = [[0.0; 3]; 2];
    let data = ThermalVolumeMomentumInput {
        volume: input(&[0.5; 2], &states, &[0.0; 2]),
        momentum_kg_m_s: &zero,
        acceleration_m_s2: [0.0; 3],
        force_density_n_m3: &zero,
        boundary_velocity_m_s: [None; 6],
    };
    let limits = super::momentum::coupled_controls(0.375);
    let gate = CancelGate::new_clock_free();
    let result = field
        .advance_isotropic_implicit_conductive_thermal_volume_momentum_with_baths(
            &domain,
            [FlowFace::Wall; 6],
            None,
            &source,
            &[
                bath(0, Face3::XMin, 303.0, 0.0),
                bath(1, Face3::XMax, 301.0, 0.0),
            ],
            data,
            limits,
            solve_controls(),
            &mut ImplicitThermalVolumeWorkspace::new(2).unwrap(),
            &gate,
        )
        .unwrap();
    let expected = field
        .advance_piecewise_thermal_volume_momentum(
            &domain,
            [FlowFace::Wall; 6],
            &caloric,
            ThermalVolumeMomentumInput {
                volume: input(&[0.5; 2], &states, &[0.1875, -0.1875]),
                ..data
            },
            limits.evolution,
            &gate,
        )
        .unwrap();
    for cell in 0..2 {
        close(
            result.conductive_enthalpy_increment_j[cell],
            [0.1875, -0.1875][cell],
        );
        close(
            result.evolution.thermal.volume.transport.mass_kg[cell],
            expected.thermal.volume.transport.mass_kg[cell],
        );
        close(
            result
                .evolution
                .thermal
                .volume
                .transport
                .specific_enthalpy_j_kg[cell],
            expected.thermal.volume.transport.specific_enthalpy_j_kg[cell],
        );
        for axis in 0..3 {
            close(
                result.evolution.momentum.momentum_kg_m_s[cell][axis],
                expected.momentum.momentum_kg_m_s[cell][axis],
            );
        }
    }
    close(
        result
            .evolution
            .thermal
            .volume
            .transport
            .projection
            .field
            .flux_x(1, 0, 0),
        0.25,
    );
    close(
        result
            .evolution
            .thermal
            .volume
            .report
            .material_enthalpy_change_j,
        0.0,
    );
    assert!(
        result
            .evolution
            .thermal
            .volume
            .report
            .required_source_pressure_work_j
            > 0.0
    );
    assert_eq!(field, FlowField::quiescent(&domain));
    assert_eq!(states[0].unwrap().specific_enthalpy_j_kg(), 2.0);
}

#[test]
fn source_work_and_cancel_refusals_publish_no_state_and_workspace_replays() {
    let material = card("T", &[(300.0, 4.0), (500.0, 4.0)]);
    let caloric = curve(
        &material,
        &[(0.0, 300.0, 1.0, 1.0), (400.0, 500.0, 1.0, 1.0)],
    );
    let source = source(&material, &caloric, [0.0, 400.0]);
    let domain = VoxelDomain::new(1, 1, 1, 1.0).unwrap();
    let field = FlowField::quiescent(&domain);
    let states = [Some(caloric.state_at_specific_enthalpy(100.0).unwrap())];
    let mut workspace = ImplicitThermalVolumeWorkspace::new(1).unwrap();
    let gate = CancelGate::new_clock_free();
    let run = |temperature, solver, workspace: &mut _, gate: &_| {
        field.advance_isotropic_implicit_conductive_thermal_volume_with_baths(
            &domain,
            [FlowFace::Wall; 6],
            &source,
            &[bath(0, Face3::XMin, temperature, 0.0)],
            input(&[1.0], &states, &[0.0]),
            controls(1.0),
            solver,
            workspace,
            gate,
        )
    };
    let mut bounded = solve_controls();
    bounded.maximum_residual_evaluations = 16;
    assert!(run(1000.0, bounded, &mut workspace, &gate).is_err());
    bounded.maximum_residual_evaluations = 1;
    assert!(run(400.0, bounded, &mut workspace, &gate).is_err());
    gate.request();
    assert!(matches!(
        run(400.0, solve_controls(), &mut workspace, &gate),
        Err(ChtError::Cancelled)
    ));
    let fresh = CancelGate::new_clock_free();
    let accepted = run(400.0, solve_controls(), &mut workspace, &fresh).unwrap();
    close(
        accepted.thermal.volume.phase_states[0]
            .unwrap()
            .temperature_k(),
        390.0,
    );
    assert_eq!(
        accepted,
        run(400.0, solve_controls(), &mut workspace, &fresh).unwrap()
    );
    assert_eq!(field, FlowField::quiescent(&domain));
    assert_eq!(states[0].unwrap().specific_enthalpy_j_kg(), 100.0);

    let insulating = card("T", &[(300.0, 0.0), (500.0, 0.0)]);
    let chart = curve(
        &insulating,
        &[(0.0, 300.0, 1.0, 1.0), (400.0, 500.0, 1.0, 1.0)],
    );
    let source = super::source(&insulating, &chart, [0.0, 400.0]);
    let state = [Some(chart.state_at_specific_enthalpy(100.0).unwrap())];
    let insulated = field
        .advance_isotropic_implicit_conductive_thermal_volume_with_baths(
            &domain,
            [FlowFace::Wall; 6],
            &source,
            &[bath(0, Face3::XMin, 1000.0, 0.0)],
            input(&[1.0], &state, &[2.0]),
            controls(1.0),
            solve_controls(),
            &mut workspace,
            &fresh,
        )
        .unwrap();
    assert_eq!(insulated.conductive_enthalpy_increment_j, [0.0]);
    close(
        insulated.thermal.volume.transport.specific_enthalpy_j_kg[0],
        102.0,
    );
    close(
        insulated
            .thermal
            .volume
            .report
            .requested_enthalpy_increment_j,
        2.0,
    );
}
