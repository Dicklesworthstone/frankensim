//! Thermal compaction and sealed gas storage share original mechanics, never heat.

use super::*;
use fs_material::thermochemical_gas::energy::entropy::GasEntropySolve;
use fs_material::thermochemical_gas::energy::parcel::ClosedGasParcel;
use fs_solid::hyper3::HyperTetError;
use fs_solid::hyper3::plastic::FinitePlasticTetError;
use fs_solid::hyper3::plastic::porous::thermal::{
    FiniteThermalPorousCavitiesStep, FiniteThermalPorousTetError,
};
use fs_solid::hyper3::pressure::cavity::{AdiabaticCavity, ClosedTetCavity};
use fs_solid::hyper3::pressure::{ClosedTetPressureSurface, TetPressure};
use std::sync::Arc;

const RHO: f64 = 6.0;

fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    std::array::from_fn(|i| a[i] - b[i])
}
fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a.iter().zip(b).map(|(a, b)| a * b).sum()
}

/// Conforming hollow tetrahedron: material volume 13/3 and gas volume 1/6.
fn shell() -> (Vec<[f64; 3]>, Vec<[usize; 4]>) {
    let nodes = vec![
        [0.0, 0.0, 0.0],
        [3.0, 0.0, 0.0],
        [0.0, 3.0, 0.0],
        [0.0, 0.0, 3.0],
        [0.5, 0.5, 0.5],
        [1.5, 0.5, 0.5],
        [0.5, 1.5, 0.5],
        [0.5, 0.5, 1.5],
    ];
    let mut tets = Vec::new();
    for [a, b, c] in [[0, 1, 2], [0, 1, 3], [0, 2, 3], [1, 2, 3]] {
        tets.push([a + 4, b + 4, c + 4, c]);
        tets.push([a + 4, b + 4, b, c]);
        tets.push([a + 4, a, b, c]);
    }
    for tet in &mut tets {
        if dot(
            sub(nodes[tet[1]], nodes[tet[0]]),
            cross(
                sub(nodes[tet[2]], nodes[tet[0]]),
                sub(nodes[tet[3]], nodes[tet[0]]),
            ),
        ) < 0.0
        {
            tet.swap(0, 1);
        }
    }
    (nodes, tets)
}
fn connected_shells() -> (Vec<[f64; 3]>, Vec<[usize; 4]>) {
    let (mut nodes, mut tets) = shell();
    let reflected = [0, 1, 2, 8, 9, 10, 11, 12];
    for node in [3, 4, 5, 6, 7] {
        let [x, y, z] = nodes[node];
        nodes.push([x, y, -z]);
    }
    for index in 0..12 {
        let mut tet = tets[index].map(|node| reflected[node]);
        tet.swap(0, 1);
        tets.push(tet);
    }
    (nodes, tets)
}
fn surface(nodes: &[[f64; 3]], tets: &[[usize; 4]], cx: &Cx<'_>) -> Arc<ClosedTetPressureSurface> {
    Arc::new(
        ClosedTetPressureSurface::new(
            nodes,
            tets,
            TetAssemblyBudget::standard(),
            4 * tets.len(),
            cx,
        )
        .unwrap(),
    )
}
fn cavity(surface: &Arc<ClosedTetPressureSurface>, first: usize, cx: &Cx<'_>) -> ClosedTetCavity {
    let faces: Vec<_> = surface
        .faces()
        .iter()
        .enumerate()
        .filter_map(|(i, face)| {
            face.iter()
                .all(|node| (first..first + 4).contains(node))
                .then_some(i)
        })
        .collect();
    assert_eq!(faces.len(), 4);
    ClosedTetCavity::new(Arc::clone(surface), &faces, cx).unwrap()
}
fn thermal_solid<'a>(
    nodes: &'a [[f64; 3]],
    tets: &'a [[usize; 4]],
    cx: &Cx<'_>,
) -> FiniteThermalPorousTetSolid<'a> {
    FiniteThermalPorousTetSolid::try_new(
        nodes,
        tets,
        tets.iter().map(|_| family()).collect(),
        &vec![300.0; tets.len()],
        TetAssemblyBudget::standard(),
        cx,
    )
    .unwrap()
}
fn binding<'a>(
    cavity: &'a ClosedTetCavity,
    gas: &'a ClosedGasParcel,
    exterior: f64,
) -> AdiabaticCavity<'a> {
    AdiabaticCavity {
        cavity,
        reference: gas,
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
        .map(|(dof, x)| (dof, (stretch - 1.0) * x))
        .collect()
}

/// Independent consistent P1 mass action, M_ab = rho*V*(1+delta_ab)/20.
fn mass_loads(
    nodes: &[[f64; 3]],
    tets: &[[usize; 4]],
    acceleration: &[[f64; 3]],
) -> Vec<(usize, f64)> {
    let mut force = vec![[0.0; 3]; nodes.len()];
    for tet in tets {
        let volume = dot(
            sub(nodes[tet[1]], nodes[tet[0]]),
            cross(
                sub(nodes[tet[2]], nodes[tet[0]]),
                sub(nodes[tet[3]], nodes[tet[0]]),
            ),
        ) / 6.0;
        assert!(volume > 0.0);
        let sum: [f64; 3] = std::array::from_fn(|i| tet.iter().map(|&n| acceleration[n][i]).sum());
        for &node in tet {
            for i in 0..3 {
                force[node][i] += RHO * volume / 20.0 * (sum[i] + acceleration[node][i]);
            }
        }
    }
    force.into_iter().flatten().enumerate().collect()
}
/// Independent entropy integral for the authored Cv/R=2.5+slope*T.
fn gas_temperature(slope: f64, volume_ratio: f64) -> f64 {
    let (mut low, mut high) = (300.0, 300.0 * volume_ratio.powf(-0.4));
    for _ in 0..80 {
        let mid = 0.5 * (low + high);
        let entropy = 2.5 * (mid / 300.0).ln() + slope * (mid - 300.0) + volume_ratio.ln();
        if entropy < 0.0 {
            low = mid;
        } else {
            high = mid;
        }
    }
    0.5 * (low + high)
}
fn gas_energy(gas: &ClosedGasParcel, slope: f64, before: f64, after: f64) -> f64 {
    gas.mass_kg() * fs_thermochem::UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K / 0.028
        * (after - before)
        * (2.5 + 0.5 * slope * (after + before))
}
fn material_account(
    result: &FiniteThermalPorousCavitiesStep<'_>,
    old: Hydro,
    new: Hydro,
    volume: f64,
    tolerance: f64,
) {
    let [work, reversible, defect, coefficients] = new.energy(old);
    let e = &result.step.energy;
    close(e.internal_endpoint_work_j, volume * work, tolerance);
    close(
        e.internal_energy_change_j,
        volume * (new.internal - old.internal),
        tolerance,
    );
    close(
        e.free_energy_change_j,
        volume * (new.free - old.free),
        tolerance,
    );
    close(e.plastic_dissipation_j, volume * new.heat, tolerance);
    close(e.reversible_heat_j, volume * reversible, tolerance);
    close(
        e.caloric_heat_j,
        volume * (new.heat + reversible),
        tolerance,
    );
    close(
        e.coefficient_storage_change_j,
        volume * coefficients,
        tolerance,
    );
    close(e.endpoint_work_defect_j, volume * defect, tolerance);
}
fn complete_balance(result: &FiniteThermalPorousCavitiesStep<'_>, tolerance: f64) {
    let e = &result.step.energy;
    let held = result.step.pressure.as_ref().map_or(0.0, |p| p.work_j);
    let correction = result
        .step
        .pressure
        .as_ref()
        .map_or(0.0, |p| p.endpoint_work_correction_j);
    let gas = result
        .cavity_energy
        .iter()
        .map(|w| w.gas_stored_energy_change_j)
        .sum::<f64>();
    let exterior = result
        .cavity_energy
        .iter()
        .map(|w| w.exterior_work_j)
        .sum::<f64>();
    let (kinetic, numerical) = result
        .step
        .inertia
        .as_ref()
        .map_or((0.0, 0.0), |i| (i.kinetic_change_j, i.numerical_loss_j));
    close(
        e.dead_load_endpoint_work_j + e.support_endpoint_work_j + held + exterior,
        e.internal_energy_change_j
            + e.caloric_heat_j
            + gas
            + kinetic
            + numerical
            + result.external_work_defect_j
            + e.balance_error_j,
        tolerance,
    );
    close(
        result.external_work_defect_j,
        e.endpoint_work_defect_j
            + correction
            + result
                .cavity_energy
                .iter()
                .map(|w| w.endpoint_work_correction_j)
                .sum::<f64>(),
        tolerance,
    );
    close(e.balance_error_j, 0.0, tolerance);
}

#[test]
fn free_heated_shell_solves_sealed_gas_and_consistent_inertia_together() {
    with_cx(|cx, _| {
        let (mut nodes, tets) = shell();
        for point in &mut nodes {
            for coordinate in point {
                *coordinate -= 0.75;
            }
        }
        let surface = surface(&nodes, &tets, cx);
        let cavity = cavity(&surface, 4, cx);
        let initial = thermal_solid(&nodes, &tets, cx)
            .with_inertia(&vec![RHO; tets.len()], &vec![[0.0; 3]; nodes.len()], cx)
            .unwrap();
        let stretch = 0.97_f64;
        let dt = 0.25;
        let old = Hydro::at(0.0, 300.0, None);
        let expected = Hydro::at(stretch.ln(), 420.0, Some(old));
        let pressure = -expected.piola / stretch.powi(2);
        let exterior = 50.0;
        let gas = gas_fixture::gas(
            cavity.reference_volume_m3(),
            (pressure + exterior) * stretch.powf(4.2),
            0.0,
            1e24,
        );
        let original = gas.clone();
        let sources = [binding(&cavity, &gas, exterior)];
        let face_pressure: Vec<_> = (0..surface.faces().len())
            .map(|face| {
                if cavity.solid_face_indices().contains(&face) {
                    0.0
                } else {
                    pressure
                }
            })
            .collect();
        // At the target, the inner gas and outer current-area pressure balance
        // the analytic thermal stress. Independent M*a supplies the nonzero
        // acceleration. Every displacement remains free, including the gas wall.
        let acceleration: Vec<_> = nodes
            .iter()
            .map(|x| x.map(|v| (stretch - 1.0) * v / (dt * dt)))
            .collect();
        let loads = mass_loads(&nodes, &tets, &acceleration);
        let result = initial
            .advance_with_cavities(
                &vec![420.0; tets.len()],
                &[],
                &loads,
                [0.0; 3],
                dt,
                &sources,
                Some((surface.as_ref(), TetPressure::PerFace(&face_pressure))),
                HyperTetSettings {
                    max_newton_iterations: 80,
                    ..settings()
                },
                controls(),
                cx,
            )
            .unwrap();
        assert!(result.step.equilibrium.reactions_n.is_empty());
        assert!(result.step.equilibrium.steps[0].linear_iterations > 0);
        for (node, x) in nodes.iter().enumerate() {
            for i in 0..3 {
                close(
                    result.step.solid.displacement_m()[node][i],
                    (stretch - 1.0) * x[i],
                    5e-7,
                );
                close(
                    result.step.solid.velocity_mps().unwrap()[node][i],
                    (stretch - 1.0) * x[i] / dt,
                    2e-6,
                );
            }
        }
        for cell in 0..tets.len() {
            close(
                result.step.solid.element_states()[cell].alpha,
                expected.alpha,
                5e-7,
            );
            close(
                result.step.solid.element_jacobians()[cell],
                stretch.powi(3),
                1e-6,
            );
            close(
                result.step.solid.element_plastic_jacobians()[cell],
                expected.log_jp.exp(),
                1e-6,
            );
        }
        let material_volume = 13.0 / 3.0;
        close(
            result.step.solid.reference_mass_kg().unwrap(),
            RHO * material_volume,
            1e-12,
        );
        for density in result
            .step
            .solid
            .current_density_kg_m3(&vec![RHO; tets.len()], cx)
            .unwrap()
        {
            close(density, RHO / stretch.powi(3), 2e-5);
        }
        assert_eq!(
            result.step.solid.reference_density_kg_m3(),
            Some(vec![RHO; tets.len()].as_slice())
        );
        assert_eq!(result.step.solid.temperatures_k(), vec![420.0; tets.len()]);
        assert_eq!(result.step.solid.families(), initial.families());
        material_account(&result, old, expected, material_volume, 3e-4);
        // The centered hollow tetrahedron has integral |X|² dV = 363/80.
        let kinetic = 0.5 * RHO * (363.0 / 80.0) * (stretch - 1.0).powi(2) / (dt * dt);
        let inertia = result.step.inertia.as_ref().unwrap();
        close(inertia.kinetic_before_j, 0.0, 0.0);
        close(inertia.kinetic_after_j, kinetic, 8e-6);
        close(inertia.kinetic_change_j, kinetic, 8e-6);
        close(inertia.numerical_loss_j, kinetic, 8e-6);
        close(inertia.inertial_work_j, 2.0 * kinetic, 1e-5);
        close(
            result.step.energy.dead_load_endpoint_work_j,
            2.0 * kinetic,
            1e-5,
        );
        assert!(kinetic > 0.19);
        let gas_t = 300.0 * stretch.powf(-1.2);
        close(
            result.cavities[0].gas.state().state().temperature,
            gas_t,
            2e-4,
        );
        assert!(gas_t < 320.0 && gas_t != expected.t);
        close(result.cavities[0].net_pressure_pa, pressure, 3e-4);
        close(
            result.cavities[0].volume_m3,
            cavity.reference_volume_m3() * stretch.powi(3),
            5e-7,
        );
        let work = result.cavity_energy[0];
        // The enormous energy datum defeats subtracting absolute stored energies.
        close(
            work.gas_stored_energy_change_j,
            gas_energy(&gas, 0.0, 300.0, gas_t),
            5e-5,
        );
        close(work.gas_work_j, -work.gas_stored_energy_change_j, 0.0);
        close(
            work.exterior_work_j,
            exterior * cavity.reference_volume_m3() * (1.0 - stretch.powi(3)),
            3e-5,
        );
        close(
            work.endpoint_work_j,
            3.0 * pressure * cavity.reference_volume_m3() * stretch.powi(2) * (stretch - 1.0),
            5e-5,
        );
        close(
            result.step.pressure.as_ref().unwrap().work_j,
            pressure * 4.5 * (1.0 - stretch.powi(3)),
            3e-4,
        );
        complete_balance(&result, 2e-6);
        assert_eq!(
            result.cavities[0].gas.mass_kg().to_bits(),
            gas.mass_kg().to_bits()
        );
        assert_eq!(result.cavities[0].gas.model(), gas.model());
        assert_eq!(gas, original);
        assert!(
            initial
                .element_states()
                .iter()
                .all(|state| state.alpha == 0.0)
        );
        assert_eq!(initial.temperatures_k(), vec![300.0; tets.len()]);
        assert_eq!(initial.velocity_mps().unwrap(), vec![[0.0; 3]; nodes.len()]);
    });
}

#[test]
fn heating_at_fixed_cavity_volume_retains_independent_gases_and_atomic_history() {
    with_cx(|cx, gate| {
        let (nodes, tets) = connected_shells();
        let surface = surface(&nodes, &tets, cx);
        let cavities = [cavity(&surface, 4, cx), cavity(&surface, 9, cx)];
        let slopes = [0.0, 0.002];
        let pressures = [100.0, 250.0];
        let exteriors = [20.0, 30.0];
        let gases = [
            gas_fixture::gas(
                cavities[0].reference_volume_m3(),
                pressures[0],
                slopes[0],
                0.0,
            ),
            gas_fixture::gas(
                cavities[1].reference_volume_m3(),
                pressures[1],
                slopes[1],
                1e24,
            ),
        ];
        let original_gases = gases.clone();
        // Reverse source order on a connected solid: neither gases nor their
        // distinct heat-capacity laws may be merged or reordered.
        let sources = [
            binding(&cavities[1], &gases[1], exteriors[1]),
            binding(&cavities[0], &gases[0], exteriors[0]),
        ];
        let initial = thermal_solid(&nodes, &tets, cx);
        let stretch = 0.97_f64;
        let supports = affine(&nodes, stretch);
        let virgin = Hydro::at(0.0, 300.0, None);
        let compressed = Hydro::at(stretch.ln(), 400.0, Some(virgin));
        let first = initial
            .advance_with_cavities(
                &vec![400.0; tets.len()],
                &supports,
                &[],
                [0.0; 3],
                0.2,
                &sources,
                None,
                settings(),
                controls(),
                cx,
            )
            .unwrap();
        let volume = 26.0 / 3.0;
        material_account(&first, virgin, compressed, volume, 3e-6);
        complete_balance(&first, 2e-6);
        assert!(first.step.inertia.is_none());
        let accepted = &first.step.solid;
        let original_history = accepted.element_states().to_vec();
        let original_u = accepted.displacement_m().to_vec();
        let original_t = accepted.temperatures_k().to_vec();
        let original_internal = accepted.internal_energy_j();
        let original_free = accepted.free_energy_j();
        let heated = Hydro::at(stretch.ln(), 450.0, Some(compressed));
        assert!(heated.alpha > compressed.alpha && heated.log_jp < compressed.log_jp);
        let second = accepted
            .advance_with_cavities(
                &vec![450.0; tets.len()],
                &supports,
                &[],
                [0.0; 3],
                0.2,
                &sources,
                None,
                settings(),
                controls(),
                cx,
            )
            .unwrap();
        material_account(&second, compressed, heated, volume, 3e-6);
        complete_balance(&second, 2e-6);
        assert!(second.step.energy.caloric_heat_j > 70.0);
        for cell in 0..tets.len() {
            close(
                second.step.solid.element_states()[cell].alpha,
                heated.alpha,
                2e-8,
            );
            close(
                second.step.solid.element_plastic_jacobians()[cell],
                heated.log_jp.exp(),
                2e-8,
            );
            close(
                second.step.solid.element_jacobians()[cell],
                stretch.powi(3),
                2e-12,
            );
        }
        for (slot, source_index) in [1, 0].into_iter().enumerate() {
            let gas_t = gas_temperature(slopes[source_index], stretch.powi(3));
            close(
                first.cavities[slot].gas.state().state().temperature,
                gas_t,
                2e-8,
            );
            close(
                first.cavity_energy[slot].gas_stored_energy_change_j,
                gas_energy(&gases[source_index], slopes[source_index], 300.0, gas_t),
                2e-8,
            );
            assert_eq!(second.cavities[slot].gas, first.cavities[slot].gas);
            close(
                second.cavities[slot].volume_m3,
                cavities[source_index].reference_volume_m3() * stretch.powi(3),
                2e-12,
            );
            close(
                second.cavity_energy[slot].gas_stored_energy_change_j,
                0.0,
                1e-12,
            );
            close(second.cavity_energy[slot].exterior_work_j, 0.0, 1e-12);
            close(second.cavity_energy[slot].endpoint_work_j, 0.0, 1e-12);
            close(
                second.cavity_energy[slot].endpoint_work_correction_j,
                0.0,
                1e-12,
            );
            assert_eq!(
                second.cavities[slot].gas.mass_kg().to_bits(),
                gases[source_index].mass_kg().to_bits()
            );
            assert_eq!(
                second.cavities[slot].gas.model(),
                gases[source_index].model()
            );
        }
        // Release partway toward the hot stress-free plastic stretch. Actual
        // gas expansion follows total geometry; retained plastic volume stays.
        let release_log = 0.5 * (stretch.ln() + thermal_log(450.0) + heated.log_jp / 3.0);
        let released = Hydro::at(release_log, 450.0, Some(heated));
        assert!(released.stretch > stretch);
        close(released.heat, 0.0, 0.0);
        let third = second
            .step
            .solid
            .advance_with_cavities(
                &vec![450.0; tets.len()],
                &affine(&nodes, released.stretch),
                &[],
                [0.0; 3],
                0.2,
                &sources,
                None,
                settings(),
                controls(),
                cx,
            )
            .unwrap();
        material_account(&third, heated, released, volume, 3e-6);
        complete_balance(&third, 2e-6);
        close(third.step.energy.plastic_dissipation_j, 0.0, 2e-8);
        for &j in third.step.solid.element_plastic_jacobians() {
            close(j, heated.log_jp.exp(), 2e-8);
        }
        for (slot, source_index) in [1, 0].into_iter().enumerate() {
            let gas_t0 = gas_temperature(slopes[source_index], stretch.powi(3));
            let gas_t1 = gas_temperature(slopes[source_index], released.stretch.powi(3));
            let work = third.cavity_energy[slot];
            close(
                third.cavities[slot].gas.state().state().temperature,
                gas_t1,
                2e-8,
            );
            close(
                work.gas_stored_energy_change_j,
                gas_energy(&gases[source_index], slopes[source_index], gas_t0, gas_t1),
                3e-8,
            );
            close(
                first.cavity_energy[slot].gas_stored_energy_change_j
                    + second.cavity_energy[slot].gas_stored_energy_change_j
                    + work.gas_stored_energy_change_j,
                gas_energy(&gases[source_index], slopes[source_index], 300.0, gas_t1),
                3e-8,
            );
            close(
                work.exterior_work_j,
                exteriors[source_index]
                    * cavities[source_index].reference_volume_m3()
                    * (stretch.powi(3) - released.stretch.powi(3)),
                2e-10,
            );
            let p_after = pressures[source_index] * gas_t1 / 300.0 / released.stretch.powi(3)
                - exteriors[source_index];
            close(
                work.endpoint_work_j,
                3.0 * p_after
                    * cavities[source_index].reference_volume_m3()
                    * released.stretch.powi(2)
                    * (released.stretch - stretch),
                3e-8,
            );
            assert!(work.gas_stored_energy_change_j < 0.0 && work.gas_work_j > 0.0);
        }

        let strict = FiniteThermoplasticControls {
            material: FinitePlasticStepControls {
                maximum_endpoint_work_defect_j: 0.0,
                ..controls().material
            },
            ..controls()
        };
        assert!(matches!(
            accepted.advance_with_cavities(
                &vec![450.0; tets.len()],
                &supports,
                &[],
                [0.0; 3],
                0.2,
                &sources,
                None,
                settings(),
                strict,
                cx,
            ),
            Err(FiniteThermalPorousTetError::Solid(
                FinitePlasticTetError::EnergyGate {
                    quantity: "thermal porous cavity/pressure/constitutive endpoint integration defect",
                    ..
                }
            ))
        ));
        let retry = accepted
            .advance_with_cavities(
                &vec![450.0; tets.len()],
                &supports,
                &[],
                [0.0; 3],
                0.2,
                &sources,
                None,
                settings(),
                controls(),
                cx,
            )
            .unwrap();
        assert_eq!(
            retry.step.solid.element_states(),
            second.step.solid.element_states()
        );
        assert_eq!(retry.step.energy, second.step.energy);
        assert_eq!(retry.cavity_energy, second.cavity_energy);
        assert_eq!(retry.external_work_defect_j, second.external_work_defect_j);
        let empty = accepted
            .advance_with_cavities::<ClosedGasParcel>(
                &vec![450.0; tets.len()],
                &supports,
                &[],
                [0.0; 3],
                0.2,
                &[],
                None,
                settings(),
                controls(),
                cx,
            )
            .unwrap();
        let plain = accepted
            .advance(
                &vec![450.0; tets.len()],
                &supports,
                &[],
                [0.0; 3],
                0.2,
                settings(),
                controls(),
                cx,
            )
            .unwrap();
        assert_eq!(
            empty.step.solid.element_states(),
            plain.solid.element_states()
        );
        assert_eq!(empty.step.energy, plain.energy);
        assert_eq!(empty.step.element_energy, plain.element_energy);
        assert!(empty.cavities.is_empty() && empty.cavity_energy.is_empty());
        close(
            empty.external_work_defect_j,
            plain.energy.endpoint_work_defect_j,
            0.0,
        );
        // A deformed output parcel is never a substitute undeformed reference.
        let rebased = [binding(&cavities[1], &first.cavities[0].gas, exteriors[1])];
        assert!(
            accepted
                .advance_with_cavities(
                    &vec![450.0; tets.len()],
                    &supports,
                    &[],
                    [0.0; 3],
                    0.2,
                    &rebased,
                    None,
                    settings(),
                    controls(),
                    cx,
                )
                .is_err()
        );
        gate.request();
        assert!(matches!(
            accepted.advance_with_cavities(
                &vec![450.0; tets.len()],
                &supports,
                &[],
                [0.0; 3],
                0.2,
                &sources,
                None,
                settings(),
                controls(),
                cx,
            ),
            Err(FiniteThermalPorousTetError::Solid(
                FinitePlasticTetError::Native(HyperTetError::Cancelled)
            ))
        ));
        assert_eq!(accepted.element_states(), original_history);
        assert_eq!(accepted.displacement_m(), original_u);
        assert_eq!(accepted.temperatures_k(), original_t);
        assert_eq!(accepted.internal_energy_j(), original_internal);
        assert_eq!(accepted.free_energy_j(), original_free);
        assert_eq!(gases, original_gases);
    });
}
