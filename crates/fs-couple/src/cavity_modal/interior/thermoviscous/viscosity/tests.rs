use super::super::{ThermoviscousTimeLimits, source_medium};
use super::*;
use fs_feec::acoustics::thermoviscous::moving_wall::MovingWall;
use fs_feec::acoustics::thermoviscous::{
    ThermoviscousBoundary, ThermoviscousFe, ThermoviscousLimits, ThermoviscousThermodynamics,
};
use fs_material::gas::{GasSpec, GasState};
use fs_math::c64::C64;
use fs_modal::time::{ModalAcousticMode, ModalAcousticTimeBudget, ModalAcousticTimeModel};
use fs_rep_mesh::TetComplex;

fn model<'c>(mesh: &'c TetComplex, points: &[[f64; 3]]) -> ThermoviscousInterior<'c> {
    let gas = GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 101325.0).unwrap();
    let medium = source_medium(&gas, 0.0).unwrap();
    let mut shape = vec![[0.0; 3]; mesh.vertex_count + mesh.edges.len()];
    shape[..mesh.vertex_count].fill([100.0, 0.0, 0.0]);
    let structure = ModalAcousticTimeModel::try_new_with_free_coordinates(
        48_000,
        vec![
            ModalAcousticMode {
                angular_frequency_rad_s: 0.0,
                damping_ratio: 0.0,
                pressure_per_modal_velocity: C64::ZERO,
            };
            2
        ],
        ModalAcousticTimeBudget::audible_reference(),
    )
    .unwrap();
    let fe = ThermoviscousFe::new_with_moving_wall(
        mesh,
        points,
        medium,
        ThermoviscousThermodynamics::ideal_gas(medium).unwrap(),
        &ThermoviscousBoundary {
            no_slip_faces: mesh.faces.clone(),
            isothermal_faces: vec![],
            ports: vec![],
        },
        MovingWall {
            faces: mesh.faces.clone(),
            velocity_shapes: vec![shape.clone(), shape],
            structure,
        },
        ThermoviscousLimits::default(),
        &mut || false,
    )
    .unwrap();
    ThermoviscousInterior::from_gas(fe, gas, [10.0, 2000.0]).unwrap()
}

#[test]
fn original_material_factor_reaches_harmonic_and_causal_heat_independently_of_reference() {
    let mesh = TetComplex::from_tets(4, vec![[0, 1, 2, 3]]);
    let points = [
        [0.0; 3],
        [0.002, 0.0, 0.0],
        [0.0, 0.002, 0.0],
        [0.0, 0.0, 0.002],
    ];
    let gate = CancelGate::new_clock_free();
    let matrix = vec![2000.0, -2000.0, -2000.0, 2000.0];
    let frequency = 100.0;
    let exact = C64::ONE / C64::new(4000.0, std::f64::consts::TAU * frequency);
    for resistance in [1e-120, 7.0, 1e120] {
        let viscosity =
            ModalPlateViscosity::from_psd_modal_matrix(matrix.clone(), 2, resistance, &gate)
                .unwrap();
        let original = model(&mesh, &points);
        let medium = original.inner().medium();
        let source = original
            .with_structural_viscosity(&viscosity, &gate)
            .unwrap();
        assert_eq!(source.inner().medium(), medium);
        for (&actual, &expected) in source
            .inner()
            .moving_wall_damping_matrix_per_s()
            .unwrap()
            .iter()
            .zip(&matrix)
        {
            assert!((actual - expected).abs() < 2e-12 * expected.abs());
        }
        let field = source
            .response_with_wall_forces(frequency, &[], &[C64::ONE, -C64::ONE], &mut || false)
            .unwrap();
        assert!(
            (field.moving_wall.modal_velocity_m_sqrt_kg_per_s[0] - exact).abs()
                < 2e-9 * exact.abs()
        );
        let expected_heat = 4000.0 * exact.norm_sq();
        assert!((field.power.structural_damping_w / expected_heat - 1.0).abs() < 2e-9);
        let prepared = source
            .prepare_time(1e-4, ThermoviscousTimeLimits::default(), &mut || false)
            .unwrap();
        let step = prepared
            .step(&prepared.zero_state(), &[], &[1.0, -1.0], &mut || false)
            .unwrap();
        assert!(step.work.structural_damping_heat_j > 0.0);
        assert!(step.work.balance_residual_j().abs() < 2e-8 * step.work.supplied_j());
    }
    let wrong_basis = ModalPlateViscosity::from_psd_modal_matrix(vec![1.0], 1, 1.0, &gate).unwrap();
    assert!(
        model(&mesh, &points)
            .with_structural_viscosity(&wrong_basis, &gate)
            .is_err()
    );
}
