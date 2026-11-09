//! Pure kinematic source observation and transactional nonlinear SGS refusal.

use super::*;
use crate::conjugate::FvSgsConfig;
use crate::sgs3::{SgsModel, SgsUse};

fn config() -> FvSgsConfig {
    FvSgsConfig {
        model: SgsModel::Smagorinsky { cs: 0.5 },
        filter_width_factor: 1.0,
        maximum_eddy_to_molecular: 0.1,
        usage: SgsUse::Exploration,
    }
}

#[test]
fn unit_rate_observer_hides_and_restores_cell_and_mouth_sgs_even_on_refusal() {
    let domain = VoxelDomain::graded([vec![0.4], vec![0.2, 0.6], vec![0.3]]).unwrap();
    let mut boundaries = [FvBoundary::Symmetry; 6];
    boundaries[0] = FvBoundary::Outlet;
    boundaries[1] = FvBoundary::Outlet;
    let settings = SimpleConfig::new(boundaries);
    let gate = CancelGate::new_clock_free();
    let mut solver = Solver::with_molecular_coefficients(&domain, 1.5, 0.03, &settings);
    solver.set_material_transport(&[0.03, 0.11], None, &gate).unwrap();
    let mass = [1.5 * domain.volume(0), 9.0 * domain.volume(1)];
    solver.set_sgs_density(&mass, &gate).unwrap();
    solver.molecular_bulk_viscosity = vec![0.05, 0.2];
    solver.pressure_work_axis = Some(0);
    solver.start_time(0.1, TimeScheme::BackwardEuler);
    solver.enable_sgs(
        FvSgsConfig { maximum_eddy_to_molecular: 1e6, ..config() },
        &gate,
    ).unwrap();
    let (minus, plus) = ([0, 0, 0], [0, 1, 0]);
    let tangent = solver.cell_face(1, minus, true);
    solver.vel[1][tangent] = 0.3;
    for (high, velocity) in [(false, 0.28), (true, 0.44)] {
        let face = solver.cell_face(0, plus, high);
        solver.vel[0][face] = velocity;
    }
    solver.update_sgs(&gate).unwrap();
    let field = solver.field();
    let expected = physical_shear_rates(&domain, &field, boundaries, Some(0), 10, &gate).unwrap();
    let molecular = solver.molecular_viscosity.clone();
    let eddy = solver.eddy_viscosity.clone();
    let molecular_storage = solver.molecular_viscosity.as_ptr();
    let eddy_storage = solver.eddy_viscosity.as_ptr();
    let mouth = [false, true].map(|high| {
        solver.pressure_mouth_sgs_eddy(1, minus, plus, high).unwrap()
    });
    assert!(mouth.iter().all(|nu| *nu > 0.0));
    let mut unit = vec![1.0; 2];
    let restored = |solver: &Solver<'_>, unit: &[f64]| {
        assert!(solver.has_sgs());
        assert_eq!(solver.molecular_viscosity, molecular);
        assert_eq!(solver.eddy_viscosity, eddy);
        assert_eq!(solver.molecular_viscosity.as_ptr(), molecular_storage);
        assert_eq!(solver.eddy_viscosity.as_ptr(), eddy_storage);
        assert_eq!(unit, &[1.0; 2]);
        assert_eq!(
            [false, true].map(|high| solver
                .pressure_mouth_sgs_eddy(1, minus, plus, high).unwrap()),
            mouth,
        );
    };
    let observed = observe_molecular_rate(&mut solver, &mut unit, &gate).unwrap();
    for (&rate, &reference) in observed.iter().zip(&expected) {
        // The independent sampler round-trips native velocity through Q/A on
        // graded faces. Preserve that physical check without requiring exact
        // equality after the separate floating-point representation.
        assert!((rate - reference).abs() <= 1e-12 * reference.abs().max(1.0));
    }
    restored(&solver, &unit);
    assert_eq!(solver.field(), field);
    let cancelled = CancelGate::new_clock_free();
    cancelled.request();
    assert!(matches!(
        observe_molecular_rate(&mut solver, &mut unit, &cancelled),
        Err(ChtError::Cancelled),
    ));
    restored(&solver, &unit);
    assert_eq!(solver.field(), field);
    let bad_face = solver.cell_face(0, minus, false);
    solver.vel[0][bad_face] = f64::NAN;
    assert!(observe_molecular_rate(&mut solver, &mut unit, &gate).is_err());
    restored(&solver, &unit);
}

fn two_cell_step(
    gate: &CancelGate,
    resolve: impl FnMut(&[f64]) -> Result<(Vec<f64>, Vec<f64>), ChtError>,
) -> Result<NativeMaterialRheologyStep, ChtError> {
    let domain = VoxelDomain::new(2, 1, 1, 1.0).unwrap();
    let field = FlowField::from_face_arrays(
        &domain, vec![0.0, 1.0, 0.0], vec![0.0; 4], vec![0.0; 4],
    );
    kick_with_source_and_sgs(
        &domain,
        NativeMaterialRheologyInput {
            field: &field,
            mass_kg: &[2.0; 2],
            momentum_kg_m_s: &[[2.0, 0.0, 0.0]; 2],
            face_mass_kg: [vec![0.0, 2.0, 0.0], vec![0.0; 4], vec![0.0; 4]],
            boundaries: [FvBoundary::Symmetry; 6],
            pressure_mouth_axis: None,
            duration_s: 0.1,
        },
        config(),
        MaterialViscosityControls {
            max_cells: 10,
            max_sweeps: 200,
            relative_tolerance: 1e-12,
            maximum_impulse_defect_kg_m_s: 1e-8,
            maximum_angular_impulse_defect_kg_m2_s: 1e-8,
            maximum_energy_defect_j: 1e-7,
        },
        gate,
        resolve,
    )
}

#[test]
fn first_exact_source_owns_sgs_cap_and_late_source_cap_or_cancel_returns_no_step() {
    let gate = CancelGate::new_clock_free();
    let mut calls = 0;
    let accepted = two_cell_step(&gate, |rates| {
        calls += 1;
        if calls == 1 {
            assert!((rates[0] - (4.0_f64 / 3.0).sqrt()).abs() < 1e-12);
            assert_eq!(rates[0], rates[1]);
        }
        Ok((vec![20.0; 2], vec![0.4; 2]))
    }).unwrap();
    // At the initial velocity rho*nu_t = sqrt(2)/2. The actual mu=20 passes
    // the .1 cap; the artificial initialization seed mu=rho=2 would refuse.
    assert_eq!(accepted.source_updates, calls);
    assert_eq!(accepted.source_updates, accepted.viscosity.report.sweeps + 1);
    assert!(accepted.viscosity.report.sgs_transfer_j > 0.0);
    let mut calls = 0;
    assert!(matches!(
        two_cell_step(&gate, |_| {
            calls += 1;
            if calls > 1 {
                return Err(refuse("late original material source refusal"));
            }
            Ok((vec![20.0; 2], vec![0.4; 2]))
        }),
        Err(ChtError::InvalidInput { field: "flow.material_viscosity", .. }),
    ));
    assert_eq!(calls, 2);
    let mut calls = 0;
    assert!(matches!(
        two_cell_step(&gate, |_| {
            calls += 1;
            Ok((vec![if calls == 1 { 20.0 } else { 0.001 }; 2], vec![0.4; 2]))
        }),
        Err(ChtError::InvalidInput { field: "simple.sgs", .. }),
    ));
    assert_eq!(calls, 2);
    let cancelled = CancelGate::new_clock_free();
    let mut calls = 0;
    assert!(matches!(
        two_cell_step(&cancelled, |_| {
            calls += 1;
            if calls > 1 {
                cancelled.request();
                return Err(refuse("simultaneous late source refusal"));
            }
            Ok((vec![20.0; 2], vec![0.4; 2]))
        }),
        Err(ChtError::Cancelled),
    ));
    assert_eq!(calls, 2);
    let replay = two_cell_step(&gate, |_| Ok((vec![20.0; 2], vec![0.4; 2]))).unwrap();
    assert_eq!(replay.viscosity.report, accepted.viscosity.report);
    assert_eq!(replay.viscosity.momentum_kg_m_s, accepted.viscosity.momentum_kg_m_s);
    assert_eq!(replay.endpoint_shear_rate_per_s, accepted.endpoint_shear_rate_per_s);
    assert_eq!(replay.source_updates, accepted.source_updates);
}
