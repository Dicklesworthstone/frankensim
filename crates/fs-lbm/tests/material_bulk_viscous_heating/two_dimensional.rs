//! Planar physical bulk viscosity, actual trace collision and shared heat balance.
use super::{
    CELL, CELLS, DT, DX, controls, fluid, limits, local_bulk, local_material, near, viscosity,
};
use fs_exec::CancelGate;
use fs_lbm::core2::CollisionModel2;
use fs_lbm::material::transport::MaterialTransportError;
use fs_lbm::material::transport::rheology::heating::BulkViscousMaterial;
use fs_lbm::material::{MaterialSlab, MaterialThermalSlab, plan_scaling_from_material};
use fs_material::fluid_transport::{
    FluidTransportState,
    viscous::{DeformationRate, NewtonianViscosity},
};

const E: [[f64; 2]; 9] = [
    [0.0, 0.0],
    [1.0, 0.0],
    [0.0, 1.0],
    [-1.0, 0.0],
    [0.0, -1.0],
    [1.0, 1.0],
    [-1.0, 1.0],
    [-1.0, -1.0],
    [1.0, -1.0],
];
const W: [f64; 9] = [
    4.0 / 9.0,
    1.0 / 9.0,
    1.0 / 9.0,
    1.0 / 9.0,
    1.0 / 9.0,
    1.0 / 36.0,
    1.0 / 36.0,
    1.0 / 36.0,
    1.0 / 36.0,
];

fn slab() -> MaterialSlab {
    plan_scaling_from_material(&fluid(1.0), DX, DT, 0.1)
        .unwrap()
        .thermal_slab(MaterialThermalSlab {
            nx: 2,
            ny: 6,
            bottom_temperature_k: 310.0,
            top_temperature_k: 300.0,
            gravity_m_s2: 0.0,
            expansion_per_k: 0.0,
            max_cells: CELLS,
        })
        .unwrap()
}
fn trace_tau(mu: f64, zeta: f64) -> f64 {
    // Physical 3-D zeta becomes zeta+mu/3 in the in-plane stress trace.
    0.5 + 3.0 * ((zeta + mu / 3.0) / 1000.0) * DT / (DX * DX)
}
fn rate(strain: [[f64; 2]; 2]) -> DeformationRate {
    DeformationRate::from_tensor([strain[0][0], strain[1][1], 0.0, strain[0][1], 0.0, 0.0]).unwrap()
}
fn populations(
    strain: [[f64; 2]; 2],
    rho: f64,
    u: [f64; 2],
    force: [f64; 2],
    shear_tau: f64,
    bulk_tau: f64,
) -> [f64; 9] {
    let equilibrium = fs_lbm::equilibrium(rho, u[0], u[1]);
    let mean = 0.5 * (strain[0][0] + strain[1][1]);
    core::array::from_fn(|q| {
        let mut contraction = 0.0;
        for a in 0..2 {
            for b in 0..2 {
                let trace = if a == b { mean } else { 0.0 };
                let stress = -(2.0 / 3.0)
                    * rho
                    * DT
                    * (shear_tau * (strain[a][b] - trace) + bulk_tau * trace)
                    - 0.5 * (u[a] * force[b] + force[a] * u[b]);
                contraction += (E[q][a] * E[q][b] - if a == b { 1.0 / 3.0 } else { 0.0 }) * stress;
            }
        }
        let first = -0.5 * (E[q][0] * force[0] + E[q][1] * force[1]);
        equilibrium[q] + W[q] * (3.0 * first + 4.5 * contraction)
    })
}
fn manufacture(
    slab: &mut MaterialSlab,
    strain: [[f64; 2]; 2],
    rho: f64,
) -> (f64, f64, [[f64; 2]; CELLS]) {
    let deformation = rate(strain);
    let gamma = deformation.shear_rate_per_s();
    let t = slab
        .temperature_si(slab.lattice.gpop[CELL].iter().sum())
        .unwrap();
    let mu = viscosity(t, gamma);
    let zeta = local_bulk(CELL, t, gamma);
    let tau_s = 0.5 + 3.0 * mu * DT / (1000.0 * DX * DX);
    let force = [0.00012, -0.00009];
    slab.lattice.grid.f[CELL] = populations(
        strain,
        rho,
        [0.03, -0.015],
        force,
        tau_s,
        trace_tau(mu, zeta),
    );
    slab.lattice.grid.fext[CELL] = force;
    slab.lattice.grid.tau[CELL] = 1.4;
    let mut forces = [[0.0; 2]; CELLS];
    forces[CELL] = force.map(|f| f * 1000.0 * DX / (DT * DT));
    let power = NewtonianViscosity::new(mu, zeta)
        .unwrap()
        .evaluate(deformation)
        .unwrap()
        .dissipation_w_m3();
    (gamma, power, forces)
}
fn bits(slab: &MaterialSlab) -> Vec<u64> {
    slab.lattice
        .grid
        .f
        .iter()
        .flatten()
        .chain(slab.lattice.gpop.iter().flatten())
        .chain(&slab.lattice.grid.tau)
        .chain(slab.lattice.grid.fext.iter().flatten())
        .chain(slab.lattice.grid.central_moment_bulk_tau.iter().flatten())
        .copied()
        .map(f64::to_bits)
        .collect()
}

#[test]
fn planar_local_sources_recover_zero_shear_trace_and_mixed_power_on_native_collision() {
    for strain in [
        [[0.0, 0.0], [0.0, 0.0]],
        [[0.0, 0.8], [0.8, 0.0]],
        [[0.6, 0.0], [0.0, 0.6]],
        [[1.2, 0.3], [0.3, -0.3]],
    ] {
        for rho in [0.8, 1.2] {
            let mut actual = slab();
            let mut reference = slab();
            let (gamma, power, forces) = manufacture(&mut actual, strain, rho);
            manufacture(&mut reference, strain, rho);
            let report = actual
                .step_with_local_bulk_rheology_and_viscous_heating(
                    &[0.0; CELLS],
                    &forces,
                    DX,
                    limits(),
                    controls(),
                    &CancelGate::new_clock_free(),
                    local_material,
                )
                .unwrap();
            let applied = report
                .transport
                .applied
                .iter()
                .find(|s| s.cell == CELL)
                .unwrap();
            near(applied.material.point().shear_rate_per_s, gamma, 3e-12);
            near(
                applied.transport.dynamic_viscosity_pa_s(),
                viscosity(applied.temperature_k, gamma),
                3e-13,
            );
            near(report.viscous.power_w_m3[CELL], power, 3e-11);
            near(
                report.viscous.viscous_heat_j,
                power * DX.powi(3) * DT,
                3e-19,
            );
            near(report.viscous.source_partition_defect_j, 0.0, 2e-19);
            for value in &report.transport.applied {
                let source = value.material.source();
                assert_eq!(source.source(), &value.transport);
                near(
                    actual
                        .lattice
                        .grid
                        .central_moment_bulk_tau_at(value.cell)
                        .unwrap(),
                    trace_tau(
                        value.transport.dynamic_viscosity_pa_s(),
                        source.bulk_viscosity_pa_s(),
                    ),
                    3e-15,
                );
            }
            // Pure planar expansion has gamma=|theta|/sqrt(3), not zero.
            if strain[0][0] == 0.6 {
                near(gamma, 1.2 / 3.0_f64.sqrt(), 2e-15);
                assert!(power > 0.0);
            }
            // Independently install actual physical trace coefficients and
            // prescribe the tensor power through the existing unheated owner.
            let mut taus = [1.0; CELLS];
            for cell in 2..CELLS - 2 {
                let t = reference
                    .temperature_si(reference.lattice.gpop[cell].iter().sum())
                    .unwrap();
                let g = if cell == CELL { gamma } else { 0.0 };
                taus[cell] = trace_tau(viscosity(t, g), local_bulk(cell, t, g));
            }
            reference
                .lattice
                .grid
                .set_central_moment_bulk_field(&taus)
                .unwrap();
            let mut heat = [0.0; CELLS];
            heat[CELL] = power;
            let direct = reference
                .step_with_material_rheology(
                    &heat,
                    &forces,
                    DX,
                    limits(),
                    controls(),
                    &CancelGate::new_clock_free(),
                    super::resolve,
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
                near(*a, *b, 3e-14);
            }
            for (a, b) in actual
                .lattice
                .gpop
                .iter()
                .flatten()
                .zip(reference.lattice.gpop.iter().flatten())
            {
                near(*a, *b, 3e-14);
            }
            near(report.transport.heat.source_j, direct.heat.source_j, 3e-19);
        }
    }
}

#[test]
fn constant_and_zero_physical_bulk_keep_the_shear_dependent_planar_trace() {
    for zeta in [0.0, 0.7] {
        let mut slab = slab();
        let strains = [[[0.6, 0.0], [0.0, 0.6]], [[1.2, 0.3], [0.3, -0.3]]];
        let mut expected = [(0.0, 0.0); 2];
        for (offset, strain) in strains.into_iter().enumerate() {
            let cell = CELL + offset;
            let deformation = rate(strain);
            let gamma = deformation.shear_rate_per_s();
            let t = slab
                .temperature_si(slab.lattice.gpop[cell].iter().sum())
                .unwrap();
            let mu = viscosity(t, gamma);
            let tau_s = 0.5 + 3.0 * mu * DT / (1000.0 * DX * DX);
            slab.lattice.grid.f[cell] =
                populations(strain, 1.0, [0.0; 2], [0.0; 2], tau_s, trace_tau(mu, zeta));
            slab.lattice.grid.tau[cell] = 1.4;
            expected[offset] = (
                gamma,
                NewtonianViscosity::new(mu, zeta)
                    .unwrap()
                    .evaluate(deformation)
                    .unwrap()
                    .dissipation_w_m3(),
            );
        }
        let report = slab
            .step_with_material_rheology_and_bulk_viscous_heating(
                &[0.0; CELLS],
                &[[0.0; 2]; CELLS],
                DX,
                zeta,
                limits(),
                controls(),
                &CancelGate::new_clock_free(),
                super::resolve,
            )
            .unwrap();
        for (offset, (gamma, power)) in expected.into_iter().enumerate() {
            let cell = CELL + offset;
            let sample = report
                .transport
                .applied
                .iter()
                .find(|v| v.cell == cell)
                .unwrap();
            near(sample.material.point().shear_rate_per_s, gamma, 3e-12);
            assert_eq!(sample.material.source().bulk_viscosity_pa_s(), zeta);
            near(report.viscous.power_w_m3[cell], power, 3e-11);
            near(
                slab.lattice.grid.central_moment_bulk_tau_at(cell).unwrap(),
                trace_tau(sample.transport.dynamic_viscosity_pa_s(), zeta),
                3e-15,
            );
        }
        assert_ne!(
            slab.lattice.grid.central_moment_bulk_tau_at(CELL),
            slab.lattice.grid.central_moment_bulk_tau_at(CELL + 1)
        );
    }
}

#[test]
fn local_planar_feedback_retains_actual_source_epochs_and_trace_field() {
    let mut slab = slab();
    let (_, _, forces) = manufacture(&mut slab, [[0.8, 0.3], [0.3, -0.1]], 1.0);
    let mut previous_temperature = None;
    for _ in 0..3 {
        let report = slab
            .step_with_local_bulk_rheology_and_viscous_heating(
                &[0.0; CELLS],
                &forces,
                2.0 * DX,
                limits(),
                controls(),
                &CancelGate::new(),
                local_material,
            )
            .unwrap();
        let applied = report
            .transport
            .applied
            .iter()
            .find(|v| v.cell == CELL)
            .unwrap();
        if let Some(t) = previous_temperature {
            assert_eq!(applied.temperature_k, t);
        }
        for value in report
            .transport
            .applied
            .iter()
            .chain(&report.transport.endpoint)
        {
            let point = value.material.point();
            let source = value.material.source();
            assert_eq!(source.source(), &value.transport);
            assert_eq!(
                source.bulk_viscosity_pa_s(),
                local_bulk(point.cell, point.temperature_k, point.shear_rate_per_s)
            );
            assert_eq!(
                value.transport.dynamic_viscosity_pa_s(),
                viscosity(point.temperature_k, point.shear_rate_per_s)
            );
        }
        for value in &report.transport.applied {
            near(
                slab.lattice
                    .grid
                    .central_moment_bulk_tau_at(value.cell)
                    .unwrap(),
                trace_tau(
                    value.transport.dynamic_viscosity_pa_s(),
                    value.material.source().bulk_viscosity_pa_s(),
                ),
                3e-15,
            );
        }
        previous_temperature = Some(
            report
                .transport
                .endpoint
                .iter()
                .find(|v| v.cell == CELL)
                .unwrap()
                .temperature_k,
        );
        assert!(report.viscous.viscous_heat_j > 0.0);
    }
}

#[test]
fn planar_trace_field_reuses_uniform_collision_and_has_the_actual_local_decay() {
    let mut uniform = fs_lbm::Grid::uniform(3, 4, 0.9);
    uniform.f[CELL] = populations([[0.5, 0.0], [0.0, 0.5]], 1.0, [0.0; 2], [0.0; 2], 0.9, 1.3);
    uniform.collision = CollisionModel2::CentralMomentBulk { bulk_tau: 1.3 };
    let mut local = uniform.clone();
    local.set_central_moment_bulk_field(&[1.3; 12]).unwrap();
    let mut a = Vec::new();
    let mut b = Vec::new();
    uniform.collide_into(&mut a);
    local.collide_into(&mut b);
    assert_eq!(a, b);
    let trace = |f: &[f64; 9]| -> f64 {
        f.iter()
            .zip(E)
            .map(|(p, e)| p * (e[0] * e[0] + e[1] * e[1]))
            .sum::<f64>()
            - (2.0 / 3.0) * f.iter().sum::<f64>()
    };
    let before_trace = trace(&local.f[CELL]);
    let mut taus = [1.3; 12];
    taus[CELL] = 0.7;
    local.set_central_moment_bulk_field(&taus).unwrap();
    local.collide_into(&mut b);
    near(
        trace(&b[CELL]),
        (1.0 - 1.0 / taus[CELL]) * before_trace,
        2e-15,
    );
    assert_ne!(a[CELL], b[CELL]);
    let selected = local.central_moment_bulk_tau.clone();
    assert!(local.set_central_moment_bulk_field(&[1.0; 11]).is_err());
    taus[CELL] = f64::NAN;
    assert!(local.set_central_moment_bulk_field(&taus).is_err());
    assert_eq!(local.central_moment_bulk_tau, selected);
}

#[test]
fn planar_bulk_refusals_preserve_populations_heat_and_the_accepted_trace_field() {
    let mut slab = slab();
    let (_, _, forces) = manufacture(&mut slab, [[0.8, 0.2], [0.2, 0.2]], 1.0);
    slab.step_with_local_bulk_rheology_and_viscous_heating(
        &[0.0; CELLS],
        &forces,
        DX,
        limits(),
        controls(),
        &CancelGate::new(),
        local_material,
    )
    .unwrap();
    let before = bits(&slab);
    let collision = slab.lattice.grid.collision;
    for zeta in [-1.0, f64::NAN, f64::INFINITY, f64::MIN_POSITIVE] {
        let result = slab
            .step_with_material_rheology_and_bulk_viscous_heating::<FluidTransportState>(
                &[0.0; CELLS],
                &forces,
                DX,
                zeta,
                limits(),
                controls(),
                &CancelGate::new(),
                |_, _| panic!("invalid bulk must refuse before a source query"),
            );
        assert!(matches!(result, Err(MaterialTransportError::Invalid(_))));
        assert_eq!(bits(&slab), before);
    }
    for endpoint in [false, true] {
        let mut at_last = false;
        let result = slab.step_with_local_bulk_rheology_and_viscous_heating(
            &[0.0; CELLS],
            &forces,
            DX,
            limits(),
            controls(),
            &CancelGate::new(),
            |point, retained| {
                if point.cell == CELLS - 3 {
                    at_last = true;
                }
                if (!endpoint && point.cell == CELLS - 3)
                    || (endpoint && at_last && point.cell == 2)
                {
                    return BulkViscousMaterial::try_new(*retained, -1.0);
                }
                local_material(point, retained)
            },
        );
        assert!(matches!(result, Err(MaterialTransportError::State { .. })));
        assert_eq!(bits(&slab), before);
    }
    let mut one = controls();
    one.maximum_iterations = 1;
    assert!(matches!(
        slab.step_with_local_bulk_rheology_and_viscous_heating(
            &[0.0; CELLS],
            &forces,
            DX,
            limits(),
            one,
            &CancelGate::new(),
            local_material,
        ),
        Err(MaterialTransportError::RheologyConvergence { .. })
    ));
    assert_eq!(bits(&slab), before);
    let gate = CancelGate::new();
    assert!(matches!(
        slab.step_with_local_bulk_rheology_and_viscous_heating(
            &[0.0; CELLS],
            &forces,
            DX,
            limits(),
            controls(),
            &gate,
            |point, retained| {
                if point.cell == CELLS - 3 {
                    gate.request();
                }
                local_material(point, retained)
            },
        ),
        Err(MaterialTransportError::Cancelled)
    ));
    assert_eq!(bits(&slab), before);
    assert_eq!(slab.lattice.grid.collision, collision);
    assert!(matches!(
        slab.step_with_material_rheology_and_viscous_heating(
            &[0.0; CELLS],
            &forces,
            DX,
            limits(),
            controls(),
            &CancelGate::new(),
            super::resolve,
        ),
        Err(MaterialTransportError::Invalid(
            "independent bulk relaxation requires explicitly sourced bulk heating"
        ))
    ));
    assert_eq!(bits(&slab), before);
}
