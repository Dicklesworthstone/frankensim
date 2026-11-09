//! Applied magnetic fields retain their Faraday drive through finite heating.

use super::*;
use fs_couple::circuit_enthalpy::{
    contacts::{heat::MaterialHeatNetwork, material::ElectricalContactNetwork},
    material::MaterialNetwork,
    spatial_induction::heat_assembly_with_boundary_candidate,
};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};

fn uniform_z_trace(problem: &SpatialVectorMagnetostatics, bz: f64) -> VectorMagneticBoundary {
    // Exact edge integral of A = Bz * (-y, x, 0) / 2.
    let edge_potential: Vec<_> = problem
        .mesh()
        .edges
        .iter()
        .map(|&[u, v]| {
            let p = problem.positions_m()[u as usize];
            let q = problem.positions_m()[v as usize];
            0.5 * bz * (p[0] * q[1] - p[1] * q[0])
        })
        .collect();
    VectorMagneticBoundary {
        potential: problem
            .exterior_faces()
            .into_iter()
            .map(|face| VectorPotentialFace {
                face,
                edge_integrals_wb: problem
                    .face_edges(face)
                    .unwrap()
                    .map(|edge| edge_potential[edge]),
            })
            .collect(),
    }
}

#[test]
fn applied_field_heats_the_original_sources_and_agrees_with_complete_heat_assembly() {
    let basis = ElectricalPropertyBasis::Resistivity;
    let (rho, alpha, dt) = (0.2, 0.005, 0.1);
    let card = source::absolute_affine_magnetic_card(rho, alpha, 1.0, 0.0, basis);
    let curve = curve(&card);
    let thermal = thermal(&curve, 0.0);
    let problem = problem(&card, basis);
    let state = initial(&problem, basis, 0.0, 1.0);
    let context = source::absolute_point(300.0);
    let cells = [MaterialCell {
        cell: CellId(7),
        material_id: MaterialId(9),
        curve: curve.curve(),
        point: ThermalMaterialPoint {
            material: &card,
            context: &context,
            temperature_axis: "T",
            liquid_fraction_axis: None,
        },
        conductivity: None,
    }];
    let regions = regions();
    let bindings = laws(&card, &problem, 5.0);
    let network = SpatialMagneticThermalNetwork {
        cells: &cells,
        links: &[],
        magnetic_regions: &regions,
        volume_tolerance_m3: 1e-12,
    };
    let boundary = uniform_z_trace(&problem, 0.4);
    let current = [[0.0; 3]; 6];
    let velocity = [[0.0; 3]; 8];
    let candidate = network
        .step_with_boundary_candidate(
            &thermal,
            dt,
            &state,
            &current,
            &velocity,
            &boundary,
            &bindings,
            thermal_controls(),
            field_controls(),
            assembly(),
            &[0.0],
            || false,
        )
        .unwrap();
    let temperature = candidate.thermal().temperature_k(CellId(7)).unwrap();
    let report = candidate.report();
    let step = &report.coupled.electrical;
    assert!(temperature > 300.0);
    assert!(step.heat_j > 0.0);
    assert!(step.magnetic_boundary_work_j > 0.0);
    near(step.impressed_work_j, 0.0, 1e-12);
    near(step.mechanical_work_j, 0.0, 1e-12);
    near(report.external_heat_j, 0.0, 1e-12);
    near(report.coupled.source.reversible_heat_j, 0.0, 1e-12);
    near(report.coupled.deposited_joule_j, step.heat_j, 2e-8);
    // The manufactured caloric source has mass 1 kg and capacity 4 J/(kg K).
    // Applied-field work supplies BOTH magnetic storage and actual Joule heat.
    near(4.0 * (temperature - 300.0), step.heat_j, 2e-8);
    near(
        4.0 * (temperature - 300.0) + step.stored_after_j - step.stored_before_j,
        step.magnetic_boundary_work_j,
        2e-8,
    );
    near(report.coupled.total_energy_residual_j, 0.0, 2e-8);
    near(report.total_energy_residual_j, 0.0, 2e-8);
    for conductor in candidate.field().conductors().iter().flatten() {
        near(
            conductor.resistivity_ohm_m(),
            rho * (1.0 + alpha * (temperature - 300.0)),
            2e-12,
        );
    }

    // The complete provider must solve the same field/temperature residual,
    // even when there are no optional thermal links, contacts or radiation.
    let bulk = MaterialNetwork {
        cells: &cells,
        resistors: &[],
        links: &[],
    };
    let electrical = ElectricalContactNetwork {
        bulk: &bulk,
        contacts: &[],
        electrical_bases: &[],
    };
    let heat_network = MaterialHeatNetwork {
        electrical: &electrical,
        thermal_contacts: &[],
        radiation: None,
    };
    let gate = CancelGate::new();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 47,
                kernel_id: 3,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        let prepared = heat_network.prepare(&cx, &thermal).unwrap();
        let complete = heat_assembly_with_boundary_candidate(
            &prepared,
            &cx,
            &thermal,
            dt,
            &state,
            &current,
            &velocity,
            &boundary,
            &regions,
            1e-12,
            &bindings,
            thermal_controls(),
            field_controls(),
            assembly(),
            &[0.0],
        )
        .unwrap();
        near(
            complete.thermal().temperature_k(CellId(7)).unwrap(),
            temperature,
            2e-9,
        );
        for (&a, &b) in complete
            .field()
            .potential_wb()
            .iter()
            .zip(candidate.field().potential_wb())
        {
            near(a, b, 2e-10);
        }
        let complete_report = complete.report();
        let complete_step = &complete_report.coupled.electrical;
        near(complete_step.heat_j, step.heat_j, 2e-9);
        near(
            complete_step.magnetic_boundary_work_j,
            step.magnetic_boundary_work_j,
            2e-9,
        );
        near(complete_report.external_heat_j, 0.0, 1e-12);
        near(complete_report.total_energy_residual_j, 0.0, 2e-8);
    });

    // A later ordinary step holds the newly accepted magnetic trace rather
    // than returning to the original zero boundary or repeating its ramp.
    let hold = network
        .step_candidate(
            candidate.thermal(),
            2.0 * dt,
            candidate.field(),
            &current,
            &velocity,
            &bindings,
            thermal_controls(),
            field_controls(),
            assembly(),
            &[0.0],
            || false,
        )
        .unwrap();
    let held = &hold.report().coupled.electrical;
    assert_eq!(held.magnetic_boundary_work_j, 0.0);
    near(
        4.0 * (hold.thermal().temperature_k(CellId(7)).unwrap() - temperature)
            + held.stored_after_j
            - held.stored_before_j,
        0.0,
        2e-8,
    );
    for prescribed in &boundary.potential {
        for (edge, expected) in problem
            .face_edges(prescribed.face)
            .unwrap()
            .into_iter()
            .zip(prescribed.edge_integrals_wb)
        {
            assert_eq!(hold.field().potential_wb()[edge], expected);
        }
    }
    assert_eq!(hold.field().time_s(), 2.0 * dt);
    assert_eq!(state.time_s(), 0.0);
    assert!(state.potential_wb().iter().all(|&a| a == 0.0));
    near(thermal.temperature_k(CellId(7)).unwrap(), 300.0, 1e-12);
}
