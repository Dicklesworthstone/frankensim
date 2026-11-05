//! Independent constitutive power through the native thermal/material owner.
//! These explicitly declared synthetic laws are analytical fixtures, not
//! empirical claims about a named liquid or gas.

use fs_exec::CancelGate;
use fs_lbm::d3q19::{CollisionModel3, E3, W3, equilibrium3, update_tau3};
use fs_lbm::material::transport::MaterialTransportError;
use fs_lbm::material::transport::rheology::{RheologyPoint, RheologySolve};
use fs_lbm::material::{MaterialSlab3, MaterialThermalSlab3, plan_scaling_from_material};
use fs_lbm::rheology::Rheology;
use fs_lbm::thermal::sources::HeatSourceLimits;
use fs_material::MaterialError;
use fs_material::fluid_transport::FluidTransportState;
use fs_material::fluid_transport::rheology::CarreauYasuda;

const DX: f64 = 0.01;
const DT: f64 = 0.01;
const CELLS: usize = 16;
const CELL: usize = 6;

fn fluid(mu: f64) -> FluidTransportState {
    FluidTransportState::try_new(1000.0, 10.0, 10.0, mu).unwrap()
}
fn viscosity(temperature: f64, rate: f64) -> f64 {
    // Source-declared reciprocal-temperature shift of an original Carreau law.
    CarreauYasuda::new(1.0, 0.2, 0.5, 0.6, 2.0, [0.0, 1e5])
        .unwrap()
        .viscosity_pa_s(rate)
        .unwrap()
        * (310.0 / temperature)
}
fn resolve(
    point: RheologyPoint,
    _: &FluidTransportState,
) -> Result<FluidTransportState, MaterialError> {
    Ok(fluid(viscosity(
        point.temperature_k,
        point.shear_rate_per_s,
    )))
}
fn controls() -> RheologySolve {
    RheologySolve {
        maximum_iterations: 100,
        relaxation: 0.8,
        absolute_relaxation_tolerance: 1e-14,
        relative_relaxation_tolerance: 1e-13,
    }
}
fn limits() -> HeatSourceLimits {
    HeatSourceLimits {
        max_cells: CELLS,
        temperature_range_k: [290.0, 330.0],
        balance_tolerance_j: 1e-10,
    }
}
fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(actual.is_finite() && expected.is_finite());
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:e} != {expected:e} (tol {tolerance:e})"
    );
}
fn bulk_tau(zeta: f64) -> f64 {
    0.5 + 4.5 * (zeta / 1000.0) * DT / (DX * DX)
}
fn slab() -> MaterialSlab3 {
    let mut slab = plan_scaling_from_material(&fluid(1.0), DX, DT, 0.1)
        .unwrap()
        .thermal_slab3(MaterialThermalSlab3 {
            nx: 2,
            ny: 1,
            nz: 6,
            bottom_temperature_k: 310.0,
            top_temperature_k: 300.0,
            gravity_m_s2: 0.0,
            expansion_per_k: 0.0,
            max_cells: CELLS,
        })
        .unwrap();
    slab.lattice
        .set_buoyancy_direction([2.0, -1.0, 3.0])
        .unwrap();
    slab.lattice.grid.set_trt_magic(3.0 / 16.0).unwrap();
    slab
}

// Construct an independent physical symmetric strain and Guo half-source.
// No production stress estimator or collision creates this manufactured input.
fn manufacture(
    slab: &mut MaterialSlab3,
    strain: [[f64; 3]; 3],
    zeta: f64,
    rho: f64,
) -> (f64, f64, [[f64; 3]; CELLS]) {
    let temperature = slab
        .temperature_si(slab.lattice.gpop[CELL].iter().sum())
        .unwrap();
    let trace: f64 = (0..3).map(|a| strain[a][a]).sum();
    let gamma = (2.0 * strain.iter().flatten().map(|d| d * d).sum::<f64>()).sqrt();
    let mu = viscosity(temperature, gamma);
    let tau = 0.5 + 3.0 * (mu / 1000.0) * DT / (DX * DX);
    let deviatoric: [[f64; 3]; 3] = core::array::from_fn(|a| {
        core::array::from_fn(|b| strain[a][b] - if a == b { trace / 3.0 } else { 0.0 })
    });
    let pi: [[f64; 3]; 3] = core::array::from_fn(|a| {
        core::array::from_fn(|b| {
            -(2.0 / 3.0)
                * rho
                * DT
                * (tau * deviatoric[a][b]
                    + if a == b {
                        bulk_tau(zeta) * trace / 3.0
                    } else {
                        0.0
                    })
        })
    });
    let trace_pi: f64 = (0..3).map(|a| pi[a][a]).sum();
    let velocity = [0.03, -0.015, 0.0];
    let force = [0.00012, -0.00009, 0.00004];
    let equilibrium = equilibrium3(rho, velocity);
    slab.lattice.grid.f[CELL] = core::array::from_fn(|q| {
        let e = [f64::from(E3[q].0), f64::from(E3[q].1), f64::from(E3[q].2)];
        let eu: f64 = e.iter().zip(velocity).map(|(e, u)| e * u).sum();
        let cpc: f64 = (0..3)
            .flat_map(|a| (0..3).map(move |b| e[a] * pi[a][b] * e[b]))
            .sum();
        let guo: f64 = (0..3)
            .map(|a| (3.0 * (e[a] - velocity[a]) + 9.0 * eu * e[a]) * force[a])
            .sum();
        equilibrium[q] + 4.5 * W3[q] * (cpc - trace_pi / 3.0) - 0.5 * W3[q] * guo
    });
    slab.lattice.grid.fext[CELL] = force;
    slab.lattice.grid.tau[CELL] = 1.4; // Warm start is deliberately not the material answer.
    let power =
        2.0 * mu * deviatoric.iter().flatten().map(|d| d * d).sum::<f64>() + zeta * trace * trace;
    let mut physical_forces = [[0.0; 3]; CELLS];
    physical_forces[CELL] = force.map(|f| f * 1000.0 * DX / (DT * DT));
    (gamma, power, physical_forces)
}

#[test]
fn independent_shear_dilation_and_mixed_heat_match_native_source_transport() {
    for strain in [
        [[0.0, 0.8, 0.0], [0.8, 0.0, -0.3], [0.0, -0.3, 0.0]],
        [[0.6, 0.0, 0.0], [0.0, 0.6, 0.0], [0.0, 0.0, 0.6]],
        [[1.2, 0.3, -0.2], [0.3, -0.3, 0.1], [-0.2, 0.1, 0.6]],
    ] {
        for (zeta, rho) in [(0.2, 0.8), (1.6, 1.2)] {
            let mut actual = slab();
            let mut reference = slab();
            let direction = actual.lattice.buoyancy_direction();
            let (gamma, power, forces) = manufacture(&mut actual, strain, zeta, rho);
            manufacture(&mut reference, strain, zeta, rho);
            let report = actual
                .step_with_material_rheology_and_bulk_viscous_heating(
                    &[0.0; CELLS],
                    &forces,
                    zeta,
                    1.3,
                    limits(),
                    controls(),
                    &CancelGate::new(),
                    resolve,
                )
                .unwrap();
            near(report.viscous.power_w_m3[CELL], power, 2e-11);
            near(
                report.viscous.viscous_heat_j,
                power * DX.powi(3) * DT,
                2e-19,
            );
            near(report.viscous.source_partition_defect_j, 0.0, 1e-19);
            let applied = report
                .transport
                .applied
                .iter()
                .find(|s| s.cell == CELL)
                .unwrap();
            near(applied.material.point().shear_rate_per_s, gamma, 2e-12);
            near(
                applied.transport.dynamic_viscosity_pa_s(),
                viscosity(applied.temperature_k, gamma),
                2e-13,
            );
            assert_eq!(actual.lattice.buoyancy_direction(), direction);
            assert_eq!(actual.lattice.grid.trt_magic(), None);
            near(
                actual
                    .lattice
                    .grid
                    .collision_model_at(CELL)
                    .kinematic_bulk_viscosity(),
                (zeta / 1000.0) * DT / (DX * DX),
                2e-16,
            );
            reference
                .lattice
                .grid
                .set_central_moment_bulk_collision(1.0 / bulk_tau(zeta), 1.3)
                .unwrap();
            let mut prescribed = [0.0; CELLS];
            prescribed[CELL] = power;
            let original = reference
                .step_with_material_rheology(
                    &prescribed,
                    &forces,
                    limits(),
                    controls(),
                    &CancelGate::new(),
                    resolve,
                )
                .unwrap();
            for (a, b) in actual
                .lattice
                .grid
                .f
                .iter()
                .flatten()
                .zip(reference.lattice.grid.f.iter().flatten())
            {
                near(*a, *b, 2e-15);
            }
            for (a, b) in actual
                .lattice
                .gpop
                .iter()
                .flatten()
                .zip(reference.lattice.gpop.iter().flatten())
            {
                near(*a, *b, 2e-14);
            }
            near(
                report.transport.heat.stored_change_j,
                original.heat.stored_change_j,
                1e-13,
            );
        }
    }
}

#[test]
fn repeated_temperature_and_shear_queries_retain_physical_bulk_collision() {
    let mut slab = slab();
    let (_, _, forces) = manufacture(
        &mut slab,
        [[0.4, 0.8, 0.0], [0.8, -0.2, 0.0], [0.0, 0.0, 0.1]],
        0.7,
        1.0,
    );
    let mut previous_endpoint = None;
    for zeta in [0.7, 0.7, 1.4, 0.3] {
        let report = slab
            .step_with_material_rheology_and_bulk_viscous_heating(
                &[0.0; CELLS],
                &forces,
                zeta,
                1.1,
                limits(),
                controls(),
                &CancelGate::new(),
                resolve,
            )
            .unwrap();
        let applied = report
            .transport
            .applied
            .iter()
            .find(|s| s.cell == CELL)
            .unwrap();
        if let Some(temperature) = previous_endpoint {
            assert_eq!(applied.temperature_k, temperature);
        }
        for value in report
            .transport
            .applied
            .iter()
            .chain(&report.transport.endpoint)
        {
            let point = value.material.point();
            assert_eq!(point.temperature_k, value.temperature_k);
            assert_eq!(
                value.transport.dynamic_viscosity_pa_s(),
                viscosity(point.temperature_k, point.shear_rate_per_s)
            );
            assert_eq!(value.material.source(), &value.transport);
        }
        previous_endpoint = Some(
            report
                .transport
                .endpoint
                .iter()
                .find(|s| s.cell == CELL)
                .unwrap()
                .temperature_k,
        );
        near(
            slab.lattice
                .grid
                .collision_model_at(CELL)
                .kinematic_bulk_viscosity(),
            (zeta / 1000.0) * DT / (DX * DX),
            2e-16,
        );
        assert!(report.viscous.viscous_heat_j > 0.0);
    }
}

fn accepted_bits(slab: &MaterialSlab3) -> Vec<u64> {
    slab.lattice
        .grid
        .f
        .iter()
        .flatten()
        .chain(slab.lattice.gpop.iter().flatten())
        .chain(&slab.lattice.grid.tau)
        .chain(slab.lattice.grid.fext.iter().flatten())
        .copied()
        .chain(slab.lattice.buoyancy_direction())
        .map(f64::to_bits)
        .collect()
}

#[test]
fn bulk_source_refusals_and_cancellation_preserve_state_and_selector() {
    let mut slab = slab();
    let (_, _, forces) = manufacture(
        &mut slab,
        [[0.8, 0.2, 0.0], [0.2, 0.2, 0.0], [0.0, 0.0, 0.2]],
        0.7,
        1.0,
    );
    let before = accepted_bits(&slab);
    let collision = slab.lattice.grid.collision_model_at(CELL);
    for (zeta, higher) in [
        (0.0, 1.0),
        (-1.0, 1.0),
        (f64::NAN, 1.0),
        (f64::INFINITY, 1.0),
        (f64::MIN_POSITIVE, 1.0),
        (0.7, 2.0),
    ] {
        let result = slab
            .step_with_material_rheology_and_bulk_viscous_heating::<FluidTransportState>(
                &[0.0; CELLS],
                &forces,
                zeta,
                higher,
                limits(),
                controls(),
                &CancelGate::new(),
                |_, _| panic!("invalid bulk must refuse before material queries"),
            );
        assert!(matches!(result, Err(MaterialTransportError::Invalid(_))));
        assert_eq!(accepted_bits(&slab), before);
        assert_eq!(slab.lattice.grid.collision_model_at(CELL), collision);
    }
    let mut reached_last_initial_cell = false;
    let result = slab.step_with_material_rheology_and_bulk_viscous_heating(
        &[0.0; CELLS],
        &forces,
        0.7,
        1.2,
        limits(),
        controls(),
        &CancelGate::new(),
        |point, reference| {
            // Source traversal returns to the first cell only after every
            // initial sample and the complete tentative native heat step.
            if point.cell == CELLS - 3 {
                reached_last_initial_cell = true;
            }
            if point.cell == 2 && reached_last_initial_cell {
                return Err(MaterialError::State {
                    what: "endpoint domain".into(),
                });
            }
            Ok(*reference)
        },
    );
    assert!(matches!(result, Err(MaterialTransportError::State { .. })));
    assert_eq!(accepted_bits(&slab), before);
    assert_eq!(slab.lattice.grid.collision_model_at(CELL), collision);
    let gate = CancelGate::new();
    gate.request();
    assert!(matches!(
        slab.step_with_material_rheology_and_bulk_viscous_heating(
            &[0.0; CELLS],
            &forces,
            0.7,
            1.2,
            limits(),
            controls(),
            &gate,
            resolve,
        ),
        Err(MaterialTransportError::Cancelled)
    ));
    assert_eq!(accepted_bits(&slab), before);
    slab.lattice
        .grid
        .set_central_moment_bulk_collision(1.0 / bulk_tau(0.7), 1.2)
        .unwrap();
    let numerical = accepted_bits(&slab);
    assert!(matches!(
        slab.step_with_material_rheology_and_viscous_heating(
            &[0.0; CELLS],
            &forces,
            limits(),
            controls(),
            &CancelGate::new(),
            resolve,
        ),
        Err(MaterialTransportError::Invalid(
            "independent bulk relaxation requires explicitly sourced bulk heating"
        ))
    ));
    assert_eq!(accepted_bits(&slab), numerical);
}

#[test]
fn plate_material_retuning_keeps_trace_rate_and_rejects_bad_selection_atomically() {
    let mut slab = slab();
    let before = accepted_bits(&slab);
    assert!(
        slab.lattice
            .grid
            .set_central_moment_bulk_collision(2.0, 1.0)
            .is_err()
    );
    assert_eq!(accepted_bits(&slab), before);
    assert_eq!(slab.lattice.grid.trt_magic(), Some(3.0 / 16.0));
    slab.lattice
        .grid
        .set_central_moment_bulk_collision(1.0 / bulk_tau(0.7), 1.2)
        .unwrap();
    let bulk = slab.lattice.grid.central_moment_bulk_rates();
    let stats = update_tau3(&mut slab.lattice.grid, Rheology::Newtonian { nu: 0.07 });
    assert_eq!(stats.floored, 0);
    assert_eq!(stats.capped, 0);
    assert_eq!(slab.lattice.grid.central_moment_bulk_rates(), bulk);
    assert_eq!(
        slab.lattice.grid.collision_model_at(CELL),
        CollisionModel3::CentralMomentBulk {
            shear_rate: 1.0 / (0.5 + 0.07 / fs_lbm::CS2),
            bulk_rate: bulk.unwrap()[0],
            higher_order_rate: 1.2,
        }
    );
    slab.lattice.grid.set_bgk_collision().unwrap();
    assert_eq!(slab.lattice.grid.central_moment_bulk_rates(), None);
    slab.lattice
        .grid
        .set_central_moment_bulk_collision(1.0, 1.2)
        .unwrap();
    slab.lattice.grid.set_trt_magic(0.2).unwrap();
    assert_eq!(slab.lattice.grid.central_moment_bulk_rates(), None);
    assert_eq!(slab.lattice.grid.trt_magic(), Some(0.2));
}
