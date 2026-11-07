//! Source gas storage, thermal compaction and native transport share one endpoint.

use super::*;
use fs_couple::iqn_ils::{IqnIlsConfig, driver::CouplingMethod};
use fs_material::thermochemical_gas::energy::entropy::GasEntropySolve;
use fs_material::thermochemical_gas::energy::parcel::ClosedGasParcel;
use fs_solid::hyper3::pressure::cavity::{AdiabaticCavity, ClosedTetCavity};
use fs_solid::hyper3::pressure::{ClosedTetPressureSurface, TetPressure};
use std::sync::Arc;

#[path = "../../../fs-solid/tests/hyper3_cavity/gas_fixture.rs"]
mod gas_fixture;

fn subtract(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    std::array::from_fn(|i| a[i] - b[i])
}

fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

fn volume(nodes: &[[f64; 3]], tet: [usize; 4]) -> f64 {
    let a = subtract(nodes[tet[1]], nodes[tet[0]]);
    let normal = cross(
        subtract(nodes[tet[2]], nodes[tet[0]]),
        subtract(nodes[tet[3]], nodes[tet[0]]),
    );
    a.iter().zip(normal).map(|(a, b)| a * b).sum::<f64>() / 6.0
}

/// A conforming hollow tetrahedron: outer V = 4.5 L^3, gas V = L^3/6.
fn shell() -> (Vec<[f64; 3]>, Vec<[usize; 4]>) {
    let nodes = [
        [0.0, 0.0, 0.0],
        [3.0, 0.0, 0.0],
        [0.0, 3.0, 0.0],
        [0.0, 0.0, 3.0],
        [0.5, 0.5, 0.5],
        [1.5, 0.5, 0.5],
        [0.5, 1.5, 0.5],
        [0.5, 0.5, 1.5],
    ]
    .map(|node| node.map(|x| x * L))
    .to_vec();
    let mut tets = Vec::new();
    for [a, b, c] in [[0, 1, 2], [0, 1, 3], [0, 2, 3], [1, 2, 3]] {
        tets.push([a + 4, b + 4, c + 4, c]);
        tets.push([a + 4, b + 4, b, c]);
        tets.push([a + 4, a, b, c]);
    }
    for tet in &mut tets {
        if volume(&nodes, *tet) < 0.0 {
            tet.swap(0, 1);
        }
    }
    (nodes, tets)
}

fn original_thermal<'a>(
    fixture: &'a Fixture,
    nodes: &[[f64; 3]],
    tets: &[[usize; 4]],
) -> DepositionSimulation<'a> {
    let events = tets
        .iter()
        .enumerate()
        .map(|(i, &tet)| {
            let volume_m3 = volume(nodes, tet);
            let start_m = std::array::from_fn(|axis| {
                tet.iter()
                    .map(|&n| nodes[n][axis])
                    .fold(f64::INFINITY, f64::min)
            });
            let end_m = std::array::from_fn(|axis| {
                tet.iter()
                    .map(|&n| nodes[n][axis])
                    .fold(f64::NEG_INFINITY, f64::max)
            });
            let surface_area_m2 = [[0, 1, 2], [0, 1, 3], [0, 2, 3], [1, 2, 3]]
                .iter()
                .map(|&[a, b, c]| {
                    let normal = cross(
                        subtract(nodes[tet[b]], nodes[tet[a]]),
                        subtract(nodes[tet[c]], nodes[tet[a]]),
                    );
                    0.5 * normal.iter().map(|x| x * x).sum::<f64>().sqrt()
                })
                .sum();
            TimedEvent {
                time_s: 0.0,
                event: Event::Birth(DepositedCell {
                    id: CellId(100 + i as u64),
                    material_id: MaterialId(0),
                    segment_id: i as u64,
                    body_id: 7,
                    mass_kg: RHO * volume_m3,
                    specific_h_j_kg: 0.0,
                    geometry: CellGeometry {
                        start_m,
                        end_m,
                        volume_m3,
                        surface_area_m2,
                    },
                    boundary: Boundary {
                        ambient_k: 300.0,
                        htc_w_m2_k: 0.0,
                        reservoir: None,
                        source_w: 0.0,
                    },
                }),
            }
        })
        .collect();
    DepositionSimulation::new(
        events,
        &[(MaterialId(0), fixture.caloric.curve() as &dyn ThermalLaw)],
        DepositionBudget {
            max_cells: tets.len(),
            max_links: 0,
            max_events: tets.len(),
            max_substeps: 1000,
        },
    )
    .unwrap()
    .advance_to(0.0, thermal_config(), || Ok(()))
    .unwrap()
    .0
}

fn bind_shell<'a>(
    fixture: &'a Fixture,
    thermal: &DepositionSimulation<'_>,
    nodes: &'a [[f64; 3]],
    tets: &'a [[usize; 4]],
    cx: &Cx<'_>,
) -> FiniteThermalPorousSolid<'a> {
    FiniteThermalPorousSolid::try_new(
        thermal,
        nodes,
        tets,
        (0..tets.len())
            .map(|i| FiniteThermalPorousCell {
                cell: CellId(100 + i as u64),
                query: fixture.query(),
            })
            .collect(),
        TetAssemblyBudget::standard(),
        cx,
    )
    .unwrap()
}

fn boundaries(
    nodes: &[[f64; 3]],
    tets: &[[usize; 4]],
    cx: &Cx<'_>,
) -> (Arc<ClosedTetPressureSurface>, ClosedTetCavity) {
    let surface = Arc::new(
        ClosedTetPressureSurface::new(
            nodes,
            tets,
            TetAssemblyBudget::standard(),
            4 * tets.len(),
            cx,
        )
        .unwrap(),
    );
    let faces: Vec<_> = surface
        .faces()
        .iter()
        .enumerate()
        .filter_map(|(i, face)| face.iter().all(|node| (4..8).contains(node)).then_some(i))
        .collect();
    assert_eq!(faces.len(), 4);
    let cavity = ClosedTetCavity::new(Arc::clone(&surface), &faces, cx).unwrap();
    (surface, cavity)
}

fn binding<'a>(
    cavity: &'a ClosedTetCavity,
    reference: &'a ClosedGasParcel,
    exterior: f64,
) -> AdiabaticCavity<'a> {
    AdiabaticCavity {
        cavity,
        reference,
        exterior_pressure_pa: exterior,
        controls: GasEntropySolve {
            minimum_temperature_k: 200.0,
            maximum_temperature_k: 600.0,
            entropy_tolerance_j_kg_k: 1e-10,
            maximum_iterations: 128,
        },
    }
}

fn affine(nodes: &[[f64; 3]], stretch: f64) -> Vec<(usize, f64)> {
    nodes
        .iter()
        .flatten()
        .enumerate()
        .map(|(i, x)| (i, (stretch - 1.0) * x))
        .collect()
}

fn cavity_settings() -> HyperTetSettings {
    HyperTetSettings {
        max_newton_iterations: 80,
        ..settings()
    }
}

/// Independent source entropy integral for Cv/R = 2.5 + slope*T.
fn gas_temperature(stretch: f64, slope: f64) -> f64 {
    let (mut lo, mut hi) = (300.0, 300.0 * stretch.powf(-1.2));
    for _ in 0..100 {
        let t = 0.5 * (lo + hi);
        if 2.5 * (t / 300.0).ln() + slope * (t - 300.0) + 3.0 * stretch.ln() > 0.0 {
            hi = t;
        } else {
            lo = t;
        }
    }
    0.5 * (lo + hi)
}

fn gas_energy(gas: &ClosedGasParcel, slope: f64, before: f64, after: f64) -> f64 {
    gas.mass_kg() * fs_thermochem::UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K / 0.028
        * (after - before)
        * (2.5 + 0.5 * slope * (after + before))
}

#[test]
fn gas_pressure_drives_a_free_wall_and_material_heat_feedback_in_both_methods() {
    context(|cx| {
        let fixture = Fixture::new(CP, 0);
        let (nodes, tets) = shell();
        let thermal = original_thermal(&fixture, &nodes, &tets);
        let solid = bind_shell(&fixture, &thermal, &nodes, &tets, cx);
        let (surface, cavity) = boundaries(&nodes, &tets, cx);
        let stretch: f64 = 0.94;
        let temperature = adiabatic_temperature(stretch, CP);
        let expected = hydrostatic(stretch, temperature);
        let pressure = -expected.first_piola_pa / stretch.powi(2);
        let exterior = 1e5;
        let gas = gas_fixture::gas(
            cavity.reference_volume_m3(),
            (pressure + exterior) * stretch.powf(4.2),
            0.0,
            1e24,
        );
        let original_gas = gas.clone();
        let sources = [binding(&cavity, &gas, exterior)];
        let fields: Vec<_> = (0..surface.faces().len())
            .map(|face| {
                if cavity.solid_face_indices().contains(&face) {
                    0.0
                } else {
                    pressure
                }
            })
            .collect();
        // The x displacement of an inner-wall node is actually unknown. Its
        // pressure traction comes only from the source gas, not the outer load.
        let supports: Vec<_> = affine(&nodes, stretch)
            .into_iter()
            .filter(|(dof, _)| *dof != 15)
            .collect();
        let solid_volume = (13.0 / 3.0) * L.powi(3);
        let outer_volume = 4.5 * L.powi(3);
        let final_gas_volume = cavity.reference_volume_m3() * stretch.powi(3);
        let gas_t = 300.0 * stretch.powf(-1.2);
        let gas_u = gas_energy(&gas, 0.0, 300.0, gas_t);
        for method in [
            CouplingMethod::RelaxedPicard,
            CouplingMethod::IqnIls(IqnIlsConfig::default()),
        ] {
            let next = solid
                .advance_with_cavities_and_coupling(
                    &thermal,
                    &(),
                    1.0,
                    &supports,
                    &[],
                    [0.0; 3],
                    &sources,
                    Some((&surface, TetPressure::PerFace(&fields))),
                    cavity_settings(),
                    controls(),
                    thermal_config(),
                    config(),
                    method,
                    cx,
                )
                .unwrap();
            let mechanics = next.solid().mechanics();
            let report = next.report();
            let energy = &report.material_energy;
            assert!(report.equilibrium.solid.steps[0].linear_iterations > 0);
            close(
                mechanics.displacement_m()[5][0],
                (stretch - 1.0) * nodes[5][0],
                2e-10,
            );
            for (i, state) in mechanics.element_states().iter().enumerate() {
                close(state.alpha, expected.alpha, 3e-8);
                close(
                    mechanics.element_plastic_jacobians()[i],
                    expected.plastic_jacobian,
                    3e-8,
                );
                close(mechanics.temperatures_k()[i], temperature, 3e-6);
                close(
                    next.thermal()
                        .temperature_k(CellId(100 + i as u64))
                        .unwrap(),
                    temperature,
                    3e-6,
                );
                close(
                    next.thermal().cells()[&CellId(100 + i as u64)].mass_kg,
                    RHO * volume(&nodes, tets[i]),
                    0.0,
                );
                assert!(std::ptr::eq(
                    mechanics.families()[i].caloric_source(),
                    &fixture.caloric
                ));
            }
            let accepted_gas = &report.equilibrium.cavities[0];
            close(accepted_gas.volume_m3, final_gas_volume, 1e-13);
            close(accepted_gas.net_pressure_pa, pressure, 0.2);
            close(accepted_gas.gas.state().state().temperature, gas_t, 2e-6);
            assert_eq!(
                accepted_gas.gas.mass_kg().to_bits(),
                gas.mass_kg().to_bits()
            );
            assert_eq!(accepted_gas.gas.model(), gas.model());
            close(energy.gas_stored_energy_change_j, gas_u, 2e-7);
            close(
                energy.exterior_work_j,
                exterior * (cavity.reference_volume_m3() - final_gas_volume),
                2e-9,
            );
            let pressure_work = pressure * outer_volume * (1.0 - stretch.powi(3));
            close(
                energy.pressure.as_ref().unwrap().work_j,
                pressure_work,
                5e-7,
            );
            let heat = solid_volume * expected.caloric_heat_j_m3();
            close(energy.material.caloric_heat_j, heat, 5e-7);
            close(report.thermal_audit.stored_energy_change_j, heat, 5e-7);
            close(report.applied_heat_j, heat, 5e-7);
            close(report.external_heat_j, 0.0, config().energy_tolerance_j);
            close(
                energy.material.endpoint_work_defect_j,
                solid_volume * expected.endpoint_work_defect_j_m3,
                5e-7,
            );
            assert!(energy.external_work_defect_j.abs() > 1.0);
            close(
                energy.material.dead_load_endpoint_work_j
                    + energy.material.support_endpoint_work_j
                    + pressure_work
                    + energy.exterior_work_j,
                energy.material.internal_energy_change_j
                    + report.thermal_audit.stored_energy_change_j
                    + gas_u
                    + energy.external_work_defect_j,
                1e-6,
            );
            close(report.thermal_residual_j, 0.0, config().energy_tolerance_j);
            close(report.combined_residual_j, 0.0, config().energy_tolerance_j);
        }
        assert_eq!(thermal.time_s(), 0.0);
        assert!(
            solid
                .mechanics()
                .element_states()
                .iter()
                .all(|state| state.alpha == 0.0)
        );
        assert_eq!(gas, original_gas);
    });
}

#[test]
fn prescribed_cavity_motion_keeps_heat_inertia_and_original_sources_atomic() {
    context(|cx| {
        let fixture = Fixture::new(CP, 0);
        let (nodes, tets) = shell();
        let thermal = original_thermal(&fixture, &nodes, &tets);
        let solid = bind_shell(&fixture, &thermal, &nodes, &tets, cx)
            .with_inertia(&vec![[0.001, -0.0005, 0.0003]; nodes.len()], cx)
            .unwrap();
        let (_, cavity) = boundaries(&nodes, &tets, cx);
        let gas = gas_fixture::gas(cavity.reference_volume_m3(), 5e6, 0.001, 1e24);
        let original_gas = gas.clone();
        let sources = [binding(&cavity, &gas, 2e5)];
        let peak: f64 = 0.94;
        let native = solid
            .advance(
                &thermal,
                &(),
                1.0,
                &affine(&nodes, peak),
                &[],
                [0.0; 3],
                cavity_settings(),
                controls(),
                thermal_config(),
                config(),
                cx,
            )
            .unwrap();
        let empty: [AdiabaticCavity<'_>; 0] = [];
        let no_cavity = solid
            .advance_with_cavities(
                &thermal,
                &(),
                1.0,
                &affine(&nodes, peak),
                &[],
                [0.0; 3],
                &empty,
                None,
                cavity_settings(),
                controls(),
                thermal_config(),
                config(),
                cx,
            )
            .unwrap();
        assert_eq!(no_cavity.thermal().cells(), native.thermal().cells());
        assert_eq!(
            no_cavity.solid().mechanics().element_states(),
            native.solid().mechanics().element_states()
        );
        assert!(no_cavity.report().equilibrium.cavities.is_empty());
        assert_eq!(
            no_cavity
                .report()
                .material_energy
                .gas_stored_energy_change_j,
            0.0
        );
        let first = solid
            .advance_with_cavities(
                &thermal,
                &(),
                1.0,
                &affine(&nodes, peak),
                &[],
                [0.0; 3],
                &sources,
                None,
                cavity_settings(),
                controls(),
                thermal_config(),
                config(),
                cx,
            )
            .unwrap();
        // Prescribed identical motion has identical material heat despite gas
        // compression, gas energy datum and different support work.
        assert_eq!(
            first.solid().mechanics().element_states(),
            native.solid().mechanics().element_states()
        );
        assert_eq!(
            first.solid().mechanics().velocity_mps(),
            native.solid().mechanics().velocity_mps()
        );
        assert_eq!(first.thermal().cells(), native.thermal().cells());
        close(
            first.report().material_energy.material.caloric_heat_j,
            native.report().material_energy.caloric_heat_j,
            1e-12,
        );
        let (thermal, solid, (), first_report) = first.into_parts();
        let old_cells = thermal.cells().clone();
        let old_states = solid.mechanics().element_states().to_vec();
        let old_velocity = solid.mechanics().velocity_mps().unwrap().to_vec();
        let old_displacement = solid.mechanics().displacement_m().to_vec();
        let old_temperatures = solid.mechanics().temperatures_k().to_vec();
        let old_internal = solid.mechanics().internal_energy_j();
        let stretch: f64 = 0.9;
        let observer = LateRefusal {
            completed_substeps: 0,
        };
        let refused = solid.advance_with_cavities(
            &thermal,
            &observer,
            2.0,
            &affine(&nodes, stretch),
            &[],
            [0.0; 3],
            &sources,
            None,
            cavity_settings(),
            controls(),
            thermal_config(),
            config(),
            cx,
        );
        assert!(matches!(
            refused,
            Err(FiniteThermoplasticError::Observer(
                DepositionError::Cancelled
            ))
        ));
        let short = solid.advance_with_cavities(
            &thermal,
            &(),
            2.0,
            &affine(&nodes, stretch),
            &[],
            [0.0; 3],
            &sources,
            None,
            cavity_settings(),
            controls(),
            thermal_config(),
            FiniteThermoplasticConfig {
                maximum_iterations: 1,
                ..config()
            },
            cx,
        );
        assert!(matches!(
            short,
            Err(FiniteThermoplasticError::NotConverged { iterations: 1, .. })
        ));
        // An endpoint parcel cannot replace the original undeformed binding.
        assert!(
            solid
                .advance_with_cavities(
                    &thermal,
                    &(),
                    2.0,
                    &affine(&nodes, stretch),
                    &[],
                    [0.0; 3],
                    &[binding(
                        &cavity,
                        &first_report.equilibrium.cavities[0].gas,
                        2e5
                    )],
                    None,
                    cavity_settings(),
                    controls(),
                    thermal_config(),
                    config(),
                    cx,
                )
                .is_err()
        );
        assert_eq!(observer.completed_substeps, 0);
        assert_eq!(thermal.cells(), &old_cells);
        assert_eq!(solid.mechanics().element_states(), old_states);
        assert_eq!(solid.mechanics().velocity_mps().unwrap(), old_velocity);
        assert_eq!(solid.mechanics().displacement_m(), old_displacement);
        assert_eq!(solid.mechanics().temperatures_k(), old_temperatures);
        assert_eq!(solid.mechanics().internal_energy_j(), old_internal);
        assert_eq!(thermal.time_s(), 1.0);
        assert_eq!(solid.time_s(), 1.0);
        let next = solid
            .advance_with_cavities(
                &thermal,
                &(),
                2.0,
                &affine(&nodes, stretch),
                &[],
                [0.0; 3],
                &sources,
                None,
                cavity_settings(),
                controls(),
                thermal_config(),
                config(),
                cx,
            )
            .unwrap();
        let report = next.report();
        let energy = &report.material_energy;
        let inertia = report.inertia.as_ref().unwrap();
        let before_t = gas_temperature(peak, 0.001);
        let after_t = gas_temperature(stretch, 0.001);
        close(
            report.equilibrium.cavities[0]
                .gas
                .state()
                .state()
                .temperature,
            after_t,
            2e-7,
        );
        close(
            energy.gas_stored_energy_change_j,
            gas_energy(&gas, 0.001, before_t, after_t),
            2e-9,
        );
        close(
            first_report.material_energy.gas_stored_energy_change_j
                + energy.gas_stored_energy_change_j,
            gas_energy(&gas, 0.001, 300.0, after_t),
            3e-9,
        );
        close(
            energy.exterior_work_j,
            2e5 * cavity.reference_volume_m3() * (peak.powi(3) - stretch.powi(3)),
            2e-12,
        );
        for (old, new) in old_states
            .iter()
            .zip(next.solid().mechanics().element_states())
        {
            assert!(new.alpha > old.alpha);
        }
        assert!(inertia.numerical_loss_j > 0.0);
        close(
            inertia.reference_mass_kg,
            RHO * (13.0 / 3.0) * L.powi(3),
            1e-17,
        );
        close(
            energy.material.dead_load_endpoint_work_j
                + energy.material.support_endpoint_work_j
                + energy.exterior_work_j,
            energy.material.internal_energy_change_j
                + report.thermal_audit.stored_energy_change_j
                + energy.gas_stored_energy_change_j
                + inertia.kinetic_change_j
                + inertia.numerical_loss_j
                + energy.external_work_defect_j,
            3e-9,
        );
        close(report.thermal_residual_j, 0.0, config().energy_tolerance_j);
        close(report.combined_residual_j, 0.0, config().energy_tolerance_j);
        assert_eq!(gas, original_gas);
        assert_eq!(thermal.cells(), &old_cells);
        assert_eq!(solid.mechanics().element_states(), old_states);
    });
}
