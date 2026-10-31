//! Authored compact-junction storage in the native coupled flow/thermal march.

use fs_exec::CancelGate;
use fs_lbm::Face3;
use fs_lbm::conjugate::{
    CompactComponent, EnergyConfig, FluidProperties, FvBoundary, SimpleConfig, SolidMaterial,
    ThermalFace, ThermalSetup, TimeScheme, UnsteadyConfig, Voxel, VoxelDomain, march_conjugate,
};

#[test]
fn coupled_compact_pulse_keeps_be_storage_for_both_momentum_schemes() {
    // The side columns connect the fluid above and below the package. Its
    // two collapsed cells have unequal face areas, but share one junction.
    let domain = VoxelDomain::graded_from_fn(
        [vec![0.5, 0.25, 0.75, 0.5], vec![1.0], vec![1.0, 0.2, 1.0]],
        |p| {
            if p[0] > 0.5 && p[0] < 1.5 && p[2] > 1.0 && p[2] < 1.2 {
                Voxel::Solid(0)
            } else {
                Voxel::Fluid
            }
        },
    )
    .unwrap();
    let fluid = FluidProperties {
        density_kg_m3: 1.0,
        specific_heat_j_kg_k: 1.0,
        conductivity_w_m_k: 1.0,
        kinematic_viscosity_m2_s: 1.0,
    };
    // Collapsed material contributes neither bulk conductance nor capacity.
    let solids = [SolidMaterial::new("collapsed package", 1e-9)];
    let mut setup = ThermalSetup::new([ThermalFace::Adiabatic; 6]);
    for cell in 0..domain.cell_count() {
        if domain.is_fluid(cell) {
            setup.fixed_temperature.push((cell, 300.0));
        }
    }
    let capacity = 6.0;
    let power = 4.0;
    setup.compact_components.push(CompactComponent {
        lo: [1, 0, 1],
        hi: [3, 1, 2],
        board_face: Face3::ZMin,
        power_w: power,
        junction_to_case_k_w: 1.0,
        junction_to_board_k_w: 1.0,
        heat_capacity_j_k: Some(capacity),
    });
    let flow = SimpleConfig::new([FvBoundary::wall(); 6]);
    let initial = vec![300.0; domain.cell_count()];
    let schedule = |time: f64| if time <= 1.5 { 1.0 } else { 0.0 };
    let (dt, steps) = (0.75, 6);
    // Each branch has R = 1 + 0.5/(k*A) = 1.5 K/W. The 1:3 face
    // areas split its conductance without duplicating the 6 J/K storage.
    let conductance = 4.0 / 3.0;
    for scheme in [TimeScheme::BackwardEuler, TimeScheme::Bdf2] {
        let mut controls = UnsteadyConfig::new(dt, steps);
        controls.scheme = scheme;
        controls.inner_tolerance = 1e-10;
        controls.average_from_step = 1;
        let solution = march_conjugate(
            &domain,
            &fluid,
            &solids,
            &setup,
            &flow,
            &controls,
            None,
            &initial,
            |_| 1.0,
            schedule,
            &EnergyConfig::default(),
            None,
            &CancelGate::new(),
        )
        .unwrap();
        assert_eq!(solution.energy_records.len(), steps);
        assert_eq!(solution.flow.records.len(), steps);
        assert_eq!(solution.flow.averaged_steps, steps);
        let (mut rise, mut summed_rise, mut supplied, mut outflow) = (0.0, 0.0, 0.0, 0.0);
        for (step, record) in solution.energy_records.iter().enumerate() {
            let time = (step + 1) as f64 * dt;
            let old_rise = rise;
            let scheduled_power = power * schedule(time);
            // Thermal storage remains backward Euler even when momentum uses BDF2.
            rise = (capacity / dt * old_rise + scheduled_power) / (capacity / dt + conductance);
            assert_eq!(record.time_s, time);
            assert!(
                (record.max_temperature_k - (300.0 + rise)).abs() < 1e-8,
                "{scheme:?}, step {}: {record:?}",
                step + 1
            );
            assert!((record.stored_energy_change_j - capacity * (rise - old_rise)).abs() < 1e-8);
            assert!((record.source_j - scheduled_power * dt).abs() < 1e-12);
            assert!((record.boundary_outflow_j - conductance * rise * dt).abs() < 1e-8);
            assert!(record.closure_j.abs() < 1e-8, "{scheme:?}: {record:?}");
            assert!(solution.flow.records[step].kinetic_energy_j.abs() < 1e-24);
            supplied += scheduled_power * dt;
            outflow += record.boundary_outflow_j;
            summed_rise += rise;
            assert!((capacity * rise + outflow - supplied).abs() < 1e-8);
        }
        for cell in 0..domain.cell_count() {
            let (expected, expected_mean) = if domain.is_fluid(cell) {
                (300.0, 300.0)
            } else {
                (300.0 + rise, 300.0 + summed_rise / steps as f64)
            };
            assert!((solution.temperature[cell] - expected).abs() < 1e-8);
            assert!((solution.mean_temperature[cell] - expected_mean).abs() < 1e-8);
        }
        assert_eq!(supplied, 6.0);
        assert_eq!(setup.compact_components[0].power_w, power);
    }
}
