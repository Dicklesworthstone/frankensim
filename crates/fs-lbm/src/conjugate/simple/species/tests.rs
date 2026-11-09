//! G0/G3 coupling checks on accepted mass and the native momentum closure.

use super::*;
use crate::conjugate::{
    DensityProjectionControls, FlowFace, FlowField, FluidProperties, FvBoundary,
    FvSgsConfig, MaterialAdvectionControls, MaterialAdvectionInput, MaterialVolumeSource,
    SimpleConfig, UnsteadyConfig, simple_unsteady_with_sgs,
};
use crate::sgs3::{SgsModel, SgsUse};

fn controls() -> SpeciesDiffusionControls<3> {
    SpeciesDiffusionControls {
        substeps: 4,
        maximum_pair_updates: 1 << 16,
        maximum_storage_bytes: 1 << 20,
        maximum_mass_defect_kg: [1e-11; 3],
        maximum_partition_defect_kg: 1e-11,
    }
}

fn carrier(domain: &VoxelDomain, field: &FlowField, mass: &[f64]) -> MaterialAdvectionStep {
    field.project_and_advect_material(
        domain,
        [FlowFace::Wall; 6],
        MaterialAdvectionInput {
            mass_kg: mass,
            specific_enthalpy_j_kg: &vec![10.0; mass.len()],
            volume_sources: &vec![MaterialVolumeSource::None; mass.len()],
            boundary_inflow: [None; 6],
        },
        MaterialAdvectionControls {
            projection: DensityProjectionControls {
                duration_s: 0.01,
                relative_tolerance: 1e-12,
                maximum_volume_defect_m3_s: 1e-10,
                maximum_energy_defect_j: 1e-9,
            },
            max_cells: mass.len(),
            max_connections: 6 * mass.len(),
            max_iterations: 200,
            relative_tolerance: 1e-12,
            maximum_mass_defect_kg: 1e-11,
            maximum_enthalpy_defect_j: 1e-9,
        },
        &CancelGate::new_clock_free(),
    ).unwrap()
}

fn fluid() -> FluidProperties {
    FluidProperties {
        density_kg_m3: 1.0,
        specific_heat_j_kg_k: 2.0,
        conductivity_w_m_k: 0.02,
        kinematic_viscosity_m2_s: 0.02,
    }
}

fn quiescent_flow(domain: &VoxelDomain, eddy: Vec<f64>) -> FvFlow {
    let fluid = fluid();
    let config = SimpleConfig::new([FvBoundary::wall(); 6]);
    let solver = super::super::Solver::new(domain, &fluid, &config);
    let mut flow = solver.finish_projected(&fluid, 0, (0.0, 0.0));
    flow.eddy_viscosity_m2_s = eddy;
    flow
}

#[test]
fn native_transient_sgs_drives_conservative_three_species_mixing() {
    let domain = VoxelDomain::new(3, 3, 2, 0.25).unwrap();
    let mut faces = [FvBoundary::wall(); 6];
    faces[3] = FvBoundary::Wall { velocity: [0.2, 0.0, 0.0] };
    let config = SimpleConfig::new(faces);
    let mut march = UnsteadyConfig::new(0.02, 2);
    march.inner_iterations = 400;
    march.inner_tolerance = 1e-8;
    let gate = CancelGate::new_clock_free();
    let flow = simple_unsteady_with_sgs(
        &domain, &fluid(), &config, &march,
        FvSgsConfig {
            model: SgsModel::Smagorinsky { cs: 0.2 },
            filter_width_factor: 1.0,
            maximum_eddy_to_molecular: 100.0,
            usage: SgsUse::Exploration,
        },
        |_| 1.0, &gate,
    ).unwrap().flow;
    assert!(flow.eddy_viscosity_m2_s.iter().any(|&nu| nu > 0.0));
    let mass: Vec<_> = (0..domain.cell_count()).map(|c| domain.volume(c)).collect();
    let accepted = carrier(&domain, &flow.field, &mass);
    let initial: Vec<_> = accepted.mass_kg.iter().enumerate().map(|(c, &m)| {
        let a = if domain.coords(c)[0] == 0 { 0.7 } else { 0.1 };
        [a * m, 0.2 * m, (0.8 - a) * m]
    }).collect();
    let result = flow.diffuse_turbulent_species(
        &domain, &accepted, &initial, 0.7, 1.0, controls(), &gate,
    ).unwrap();
    let d: Vec<_> = flow.eddy_viscosity_m2_s.iter().map(|nu| nu / 0.7).collect();
    let reference = accepted.diffuse_species(
        &domain, &initial, &d, 1.0, controls(), &gate,
    ).unwrap();
    assert_eq!(result, reference);
    assert_ne!(result.species_mass_kg, initial);
    for k in 0..3 {
        assert!(result.report.mass_change_kg[k].abs() < 1e-13);
    }
    for (c, values) in result.species_mass_kg.iter().enumerate() {
        assert!(values.iter().all(|v| v.is_finite() && *v >= 0.0));
        assert!((values.iter().sum::<f64>() - accepted.mass_kg[c]).abs() < 1e-13);
    }
    // A component uniform in mass fraction must not follow the gradients of
    // the other components or acquire an invented composition source.
    for (c, values) in result.species_mass_kg.iter().enumerate() {
        assert!((values[1] - initial[c][1]).abs() < 1e-15);
    }
}

#[test]
fn schmidt_uses_momentum_eddy_field_independently_of_heat_wall_policy() {
    let domain = VoxelDomain::graded([vec![0.5, 1.5], vec![2.0], vec![1.0]]).unwrap();
    let mut flow = quiescent_flow(&domain, vec![0.03, 0.12]);
    let accepted = carrier(&domain, &flow.field, &[1.0, 3.0]);
    let initial = [[0.7, 0.2, 0.1], [0.3, 0.9, 1.8]];
    let gate = CancelGate::new_clock_free();
    let run = |flow: &FvFlow, sc| flow.diffuse_turbulent_species(
        &domain, &accepted, &initial, sc, 0.7, controls(), &gate,
    ).unwrap();
    let original = run(&flow, 0.5);
    flow.heat_eddy_viscosity_m2_s.fill(9.0);
    flow.turbulent_prandtl = 4.0;
    assert_eq!(run(&flow, 0.5), original);
    for eddy in &mut flow.eddy_viscosity_m2_s { *eddy *= 2.0; }
    assert_eq!(run(&flow, 1.0), original);
    assert!(run(&flow, 0.5).species_mass_kg[0][0] < original.species_mass_kg[0][0]);
    assert_eq!(accepted.mass_kg, [1.0, 3.0]);
    assert_eq!(accepted.specific_enthalpy_j_kg, [10.0; 2]);
}

#[test]
fn zero_turbulence_is_exact_and_invalid_input_never_changes_state() {
    let domain = VoxelDomain::new(2, 1, 1, 1.0).unwrap();
    let mut flow = quiescent_flow(&domain, vec![0.0; 2]);
    let accepted = carrier(&domain, &flow.field, &[1.0; 2]);
    let before = accepted.clone();
    let initial = [[0.7, 0.2, 0.1], [0.1, 0.3, 0.6]];
    let gate = CancelGate::new_clock_free();
    let disabled = flow.diffuse_turbulent_species(
        &domain, &accepted, &initial, 0.7, 100.0, controls(), &gate,
    ).unwrap();
    assert_eq!(disabled.species_mass_kg, initial);
    assert_eq!(disabled.report.pair_updates, 0);
    for schmidt in [0.0, -1.0, f64::INFINITY, f64::NAN] {
        assert!(flow.diffuse_turbulent_species(
            &domain, &accepted, &initial, schmidt, 1.0, controls(), &gate,
        ).is_err());
    }
    let mut limited = controls();
    limited.maximum_storage_bytes = 15;
    assert!(flow.diffuse_turbulent_species(
        &domain, &accepted, &initial, 0.7, 1.0, limited, &gate,
    ).is_err());
    for eddy in [-0.1, f64::NAN, f64::INFINITY] {
        flow.eddy_viscosity_m2_s[0] = eddy;
        assert!(flow.diffuse_turbulent_species(
            &domain, &accepted, &initial, 0.7, 1.0, controls(), &gate,
        ).is_err());
    }
    flow.eddy_viscosity_m2_s[0] = f64::MIN_POSITIVE;
    assert!(flow.diffuse_turbulent_species(
        &domain, &accepted, &initial, f64::MAX, 1.0, controls(), &gate,
    ).is_err());
    flow.eddy_viscosity_m2_s[0] = 0.0;
    gate.request();
    assert!(matches!(flow.diffuse_turbulent_species(
        &domain, &accepted, &initial, 0.7, 1.0, controls(), &gate,
    ), Err(ChtError::Cancelled)));
    assert_eq!(accepted, before);
}
