//! Sealed source-gas work and retained foam compaction through one native solve.
//! These synthetic coefficients test physical identities, not measured accuracy.

#[path = "hyper_porous_cavity/dynamics.rs"]
mod dynamics;

use std::sync::Arc;

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_material::plastic::finite::porous::{
    FinitePorous, FinitePorousLimits, FinitePorousParameters,
};
use fs_material::thermochemical_gas::energy::entropy::GasEntropySolve;
use fs_material::thermochemical_gas::energy::parcel::ClosedGasParcel;
use fs_solid::hyper3::plastic::porous::{FinitePorousTetError, FinitePorousTetSolid};
use fs_solid::hyper3::plastic::{
    FinitePlasticLawField, FinitePlasticStepControls, FinitePlasticTetError,
};
use fs_solid::hyper3::pressure::cavity::{AdiabaticCavity, ClosedTetCavity};
use fs_solid::hyper3::pressure::{ClosedTetPressureSurface, TetPressure};
use fs_solid::hyper3::{HyperTetError, HyperTetSettings};
use fs_solid::linear3::TetAssemblyBudget;

#[path = "hyper3_cavity/gas_fixture.rs"]
mod gas_fixture;

const MU: f64 = 80.0;
const K: f64 = 160.0;
const Y: f64 = 2.0;
const H: f64 = 8.0;
const A: f64 = 1.2;
const FIXED: [(usize, f64); 6] = [(0, 0.0), (1, 0.0), (2, 0.0), (4, 0.0), (5, 0.0), (8, 0.0)];

fn with_cx(f: impl FnOnce(&Cx<'_>, &CancelGate)) {
    let gate = CancelGate::new();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        f(
            &Cx::new(
                &gate,
                arena,
                StreamKey {
                    seed: 47,
                    kernel_id: 601,
                    tile: 0,
                    iteration: 0,
                },
                Budget::INFINITE,
                ExecMode::Deterministic,
            ),
            &gate,
        );
    });
}

fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        actual.is_finite() && expected.is_finite() && (actual - expected).abs() <= tolerance,
        "{actual:e} != {expected:e}; tolerance={tolerance:e}",
    );
}

fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    std::array::from_fn(|i| a[i] - b[i])
}

fn triple(a: [f64; 3], b: [f64; 3], c: [f64; 3]) -> f64 {
    a[0] * (b[1] * c[2] - b[2] * c[1])
        + a[1] * (b[2] * c[0] - b[0] * c[2])
        + a[2] * (b[0] * c[1] - b[1] * c[0])
}

/// One conforming hollow tetrahedron, V_outer=9/2, V_gas=1/6.
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
        if triple(
            sub(nodes[tet[1]], nodes[tet[0]]),
            sub(nodes[tet[2]], nodes[tet[0]]),
            sub(nodes[tet[3]], nodes[tet[0]]),
        ) < 0.0
        {
            tet.swap(0, 1);
        }
    }
    (nodes, tets)
}

/// Reflect the conforming shell across its outer z=0 triangle. The chambers
/// share solid vertices and a face, but never share boundary faces or gas.
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

fn material() -> FinitePorous {
    FinitePorous::new(
        FinitePorousParameters {
            shear_modulus_pa: MU,
            bulk_modulus_pa: K,
            yield_stress_pa: Y,
            hardening_modulus_pa: H,
            pressure_sensitivity: A,
        },
        FinitePorousLimits {
            minimum_jacobian: 0.2,
            maximum_elastic_log_strain: 0.5,
            maximum_plastic_log_strain: 2.0,
            maximum_accumulated_plastic_strain: 3.0,
        },
    )
    .unwrap()
}

fn solid<'a>(
    nodes: &'a [[f64; 3]],
    tets: &'a [[usize; 4]],
    material: &'a FinitePorous,
    cx: &Cx<'_>,
) -> FinitePorousTetSolid<'a> {
    FinitePorousTetSolid::try_new(
        nodes,
        tets,
        FinitePlasticLawField::Uniform(material),
        TetAssemblyBudget::standard(),
        cx,
    )
    .unwrap()
}

fn controls() -> FinitePlasticStepControls {
    FinitePlasticStepControls {
        equilibrium_work_tolerance_j: 2e-6,
        maximum_endpoint_work_defect_j: 20.0,
        yield_tolerance_pa: 2e-6,
    }
}

fn gas_controls() -> GasEntropySolve {
    GasEntropySolve {
        minimum_temperature_k: 200.0,
        maximum_temperature_k: 600.0,
        entropy_tolerance_j_kg_k: 1e-10,
        maximum_iterations: 128,
    }
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
        controls: gas_controls(),
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

/// Independent hydrostatic Deshpande--Fleck/Hencky return.
fn hydrostatic(stretch: f64) -> (f64, f64, f64) {
    let beta = A / (1.0 + A * A / 9.0).sqrt();
    let alpha = (-3.0 * K * beta * stretch.ln() - Y) / (H + K * beta * beta);
    let plastic_j = (-beta * alpha).exp();
    let tau = -(Y + H * alpha) / beta;
    (alpha, plastic_j, tau)
}

/// Independent entropy integral of the authored Cv/R = 2.5 + slope*T.
fn temperature(slope: f64, volume_ratio: f64) -> f64 {
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

#[test]
fn g1_sealed_cavity_and_outer_pressure_drive_native_permanent_compaction() {
    with_cx(|cx, _| {
        let (nodes, tets) = shell();
        let surface = surface(&nodes, &tets, cx);
        let cavity = cavity(&surface, 4, cx);
        let material = material();
        let initial = solid(&nodes, &tets, &material, cx);
        let stretch = 0.94_f64;
        let (alpha, plastic_j, tau) = hydrostatic(stretch);
        let pressure = -tau / stretch.powi(3);
        let exterior = 0.5;
        // At the target every solid boundary has the same compressive traction.
        // The inner pressure evolves from its original sealed source.
        let gas = gas_fixture::gas(
            cavity.reference_volume_m3(),
            (pressure + exterior) * stretch.powf(4.2),
            0.0,
            1e24,
        );
        let original = gas.clone();
        let fields: Vec<_> = (0..surface.faces().len())
            .map(|face| {
                if cavity.solid_face_indices().contains(&face) {
                    0.0
                } else {
                    pressure
                }
            })
            .collect();
        let result = initial
            .advance_with_cavities(
                &FIXED,
                &[],
                [0.0; 3],
                &[binding(&cavity, &gas, exterior)],
                Some((&surface, TetPressure::PerFace(&fields))),
                HyperTetSettings {
                    max_newton_iterations: 80,
                    ..HyperTetSettings::default()
                },
                controls(),
                cx,
            )
            .unwrap();
        assert!(result.step.equilibrium.steps[0].linear_iterations > 0);
        for (&actual, &x) in result
            .step
            .solid
            .displacement_m()
            .iter()
            .flatten()
            .zip(nodes.iter().flatten())
        {
            near(actual, (stretch - 1.0) * x, 8e-7);
        }
        for (&j, state) in result
            .step
            .solid
            .element_plastic_jacobians()
            .iter()
            .zip(result.step.solid.element_states())
        {
            near(j, plastic_j, 8e-7);
            near(state.alpha, alpha, 8e-7);
        }
        near(
            result.step.energy.plastic_dissipation_j,
            (13.0 / 3.0) * Y * alpha,
            2e-6,
        );
        let volume = cavity.reference_volume_m3() * stretch.powi(3);
        let gas_temperature = 300.0 * stretch.powf(-1.2);
        near(result.cavities[0].volume_m3, volume, 5e-7);
        near(
            result.cavities[0].gas.state().state().temperature,
            gas_temperature,
            8e-5,
        );
        near(result.cavities[0].net_pressure_pa, pressure, 8e-6);
        let work = result.cavity_energy[0];
        near(
            work.gas_stored_energy_change_j,
            gas_energy(&gas, 0.0, 300.0, gas_temperature),
            2e-6,
        );
        near(
            work.exterior_work_j,
            exterior * (cavity.reference_volume_m3() - volume),
            2e-7,
        );
        near(work.work_j, work.gas_work_j + work.exterior_work_j, 1e-13);
        near(
            work.endpoint_work_j,
            pressure * 3.0 * volume * (stretch - 1.0) / stretch,
            2e-6,
        );
        near(
            result.external_work_defect_j + result.step.energy.equilibrium_work_residual_j,
            result.step.energy.dead_load_endpoint_work_j
                + result.step.energy.support_endpoint_work_j
                + result.step.pressure.as_ref().unwrap().work_j
                + work.work_j
                - (result.step.energy.stored_after_j - result.step.energy.stored_before_j)
                - result.step.energy.plastic_dissipation_j,
            3e-6,
        );
        assert!(
            initial
                .element_states()
                .iter()
                .all(|state| state.alpha == 0.0)
        );
        assert_eq!(gas, original);
    });
}

#[test]
fn g1_original_cavity_isentrope_survives_retained_compaction_and_elastic_unloading() {
    with_cx(|cx, _| {
        let (nodes, tets) = shell();
        let surface = surface(&nodes, &tets, cx);
        let cavity = cavity(&surface, 4, cx);
        let material = material();
        let initial = solid(&nodes, &tets, &material, cx);
        let gas = gas_fixture::gas(cavity.reference_volume_m3(), 5.0, 0.001, 1e24);
        let original = gas.clone();
        let source = [binding(&cavity, &gas, 1.0)];
        let peak = 0.9_f64;
        let first = initial
            .advance_with_cavities(
                &affine(&nodes, peak),
                &[],
                [0.0; 3],
                &source,
                None,
                HyperTetSettings::default(),
                controls(),
                cx,
            )
            .unwrap();
        let unloaded = first.step.solid.element_plastic_jacobians()[0].cbrt();
        let stretch = (peak * unloaded).sqrt();
        let old_states = first.step.solid.element_states().to_vec();
        let second = first
            .step
            .solid
            .advance_with_cavities(
                &affine(&nodes, stretch),
                &[],
                [0.0; 3],
                &source,
                None,
                HyperTetSettings::default(),
                controls(),
                cx,
            )
            .unwrap();
        let t1 = temperature(0.001, peak.powi(3));
        let t2 = temperature(0.001, stretch.powi(3));
        near(second.cavities[0].gas.state().state().temperature, t2, 2e-7);
        near(
            second.cavity_energy[0].gas_stored_energy_change_j,
            gas_energy(&gas, 0.001, t1, t2),
            2e-9,
        );
        near(
            first.cavity_energy[0].gas_stored_energy_change_j
                + second.cavity_energy[0].gas_stored_energy_change_j,
            gas_energy(&gas, 0.001, 300.0, t2),
            3e-9,
        );
        near(second.step.energy.plastic_dissipation_j, 0.0, 2e-9);
        for (before, after) in old_states.iter().zip(second.step.solid.element_states()) {
            near(after.alpha, before.alpha, 2e-9);
        }
        assert_eq!(first.step.solid.element_states(), old_states);
        assert_eq!(
            second.cavities[0].gas.mass_kg().to_bits(),
            gas.mass_kg().to_bits()
        );
        assert_eq!(second.cavities[0].gas.model(), gas.model());
        assert_eq!(gas, original);
        // A deformed parcel must never silently become an undeformed reference.
        assert!(
            first
                .step
                .solid
                .advance_with_cavities(
                    &affine(&nodes, stretch),
                    &[],
                    [0.0; 3],
                    &[binding(&cavity, &first.cavities[0].gas, 1.0)],
                    None,
                    HyperTetSettings::default(),
                    controls(),
                    cx,
                )
                .is_err()
        );
        assert_eq!(first.step.solid.element_states(), old_states);
    });
}

#[test]
fn g1_connected_chambers_keep_independent_inventories_work_and_source_order() {
    with_cx(|cx, _| {
        let (nodes, tets) = connected_shells();
        let surface = surface(&nodes, &tets, cx);
        let cavities = [cavity(&surface, 4, cx), cavity(&surface, 9, cx)];
        let gases = [
            gas_fixture::gas(cavities[0].reference_volume_m3(), 5.0, 0.0, 0.0),
            gas_fixture::gas(cavities[1].reference_volume_m3(), 7.0, 0.002, 1e24),
        ];
        let original = gases.clone();
        let slopes = [0.0, 0.002];
        let exterior = [0.5, 1.5];
        let material = material();
        let initial = solid(&nodes, &tets, &material, cx);
        let stretch = 0.96_f64;
        let mut forward = None;
        for order in [[0, 1], [1, 0]] {
            let sources = order.map(|i| binding(&cavities[i], &gases[i], exterior[i]));
            let result = initial
                .advance_with_cavities(
                    &affine(&nodes, stretch),
                    &[],
                    [0.0; 3],
                    &sources,
                    None,
                    HyperTetSettings::default(),
                    controls(),
                    cx,
                )
                .unwrap();
            if let Some((stored, loss, defect)) = forward {
                near(result.step.energy.stored_after_j, stored, 0.0);
                near(result.step.energy.plastic_dissipation_j, loss, 0.0);
                near(result.external_work_defect_j, defect, 2e-12);
            } else {
                forward = Some((
                    result.step.energy.stored_after_j,
                    result.step.energy.plastic_dissipation_j,
                    result.external_work_defect_j,
                ));
            }
            for (slot, &i) in order.iter().enumerate() {
                let gas = &result.cavities[slot].gas;
                let t = temperature(slopes[i], stretch.powi(3));
                near(gas.state().state().temperature, t, 2e-7);
                assert_eq!(gas.mass_kg().to_bits(), gases[i].mass_kg().to_bits());
                assert_eq!(gas.model(), gases[i].model());
                near(
                    result.cavity_energy[slot].gas_stored_energy_change_j,
                    gas_energy(&gases[i], slopes[i], 300.0, t),
                    2e-9,
                );
                near(
                    result.cavity_energy[slot].exterior_work_j,
                    exterior[i] * cavities[i].reference_volume_m3() * (1.0 - stretch.powi(3)),
                    2e-12,
                );
            }
        }
        assert_eq!(gases, original);
    });
}

#[test]
fn g3_refused_cavity_work_sources_and_cancellation_leave_both_owners_retryable() {
    with_cx(|cx, gate| {
        let (nodes, tets) = shell();
        let surface = surface(&nodes, &tets, cx);
        let cavity = cavity(&surface, 4, cx);
        let material = material();
        let initial = solid(&nodes, &tets, &material, cx);
        let gas = gas_fixture::gas(cavity.reference_volume_m3(), 5.0, 0.0, 0.0);
        let original = gas.clone();
        let source = binding(&cavity, &gas, 0.5);
        let target = affine(&nodes, 0.9);
        let old_states = initial.element_states().to_vec();
        let old_u = initial.displacement_m().to_vec();
        let mut strict = controls();
        strict.maximum_endpoint_work_defect_j = 0.0;
        let refused = initial
            .advance_with_cavities(
                &target,
                &[],
                [0.0; 3],
                &[source],
                None,
                HyperTetSettings::default(),
                strict,
                cx,
            )
            .unwrap_err();
        assert!(matches!(
            refused,
            FinitePorousTetError::Solid(FinitePlasticTetError::EnergyGate { .. })
        ));
        // Duplicate inventory binding cannot count one surface twice.
        assert!(
            initial
                .advance_with_cavities(
                    &target,
                    &[],
                    [0.0; 3],
                    &[source, source],
                    None,
                    HyperTetSettings::default(),
                    controls(),
                    cx,
                )
                .is_err()
        );
        let mut narrow = source;
        narrow.controls.maximum_temperature_k = 301.0;
        assert!(
            initial
                .advance_with_cavities(
                    &target,
                    &[],
                    [0.0; 3],
                    &[narrow],
                    None,
                    HyperTetSettings::default(),
                    controls(),
                    cx,
                )
                .is_err()
        );
        let retry = initial
            .advance_with_cavities(
                &target,
                &[],
                [0.0; 3],
                &[source],
                None,
                HyperTetSettings::default(),
                controls(),
                cx,
            )
            .unwrap();
        assert!(retry.step.energy.plastic_dissipation_j > 0.0);
        // The empty route must keep the existing native endpoint bit-for-bit.
        let empty = initial
            .advance_with_cavities::<ClosedGasParcel>(
                &target,
                &[],
                [0.0; 3],
                &[],
                None,
                HyperTetSettings::default(),
                controls(),
                cx,
            )
            .unwrap();
        let plain = initial
            .advance(
                &target,
                &[],
                [0.0; 3],
                HyperTetSettings::default(),
                controls(),
                cx,
            )
            .unwrap();
        assert_eq!(empty.step.equilibrium, plain.equilibrium);
        assert_eq!(empty.step.energy, plain.energy);
        assert!(empty.cavities.is_empty() && empty.cavity_energy.is_empty());
        gate.request();
        let cancelled = initial
            .advance_with_cavities(
                &target,
                &[],
                [0.0; 3],
                &[source],
                None,
                HyperTetSettings::default(),
                controls(),
                cx,
            )
            .unwrap_err();
        assert!(matches!(
            cancelled,
            FinitePorousTetError::Solid(FinitePlasticTetError::Native(HyperTetError::Cancelled))
        ));
        assert_eq!(initial.element_states(), old_states);
        assert_eq!(initial.displacement_m(), old_u);
        assert_eq!(gas, original);
    });
}
