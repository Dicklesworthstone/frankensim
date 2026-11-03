//! G1/G3/G4 thermal-plastic memory, impulse work and native continuation.
//! Synthetic coefficients; inertia oracles use original simplex/tensor moments.
use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_material::plastic::finite::thermal::{
    AffineThermalFiniteJ2, FiniteJ2TemperatureSlopes, FiniteJ2ThermalParameters,
};
use fs_material::plastic::finite::{FiniteJ2Limits, FiniteJ2Parameters};
use fs_solid::hyper3::plastic::FinitePlasticStepControls;
use fs_solid::hyper3::plastic::hex::dynamics::FiniteThermoplasticHexMotion;
use fs_solid::hyper3::plastic::hex::{FiniteHexFormulation, FiniteThermoplasticHexSolid};
use fs_solid::hyper3::plastic::thermal::{
    FiniteThermoplasticControls, FiniteThermoplasticTetSolid, SolidImpulseControls,
};
use fs_solid::plastic3::IncrementSettings;
use fs_solid::plastic3::hex::HexAssemblyBudget;
use fs_solid::{HyperTetSettings, TetAssemblyBudget};
use fs_solver::CheckedCgConfig;

const TET_NODES: [[f64; 3]; 4] = [[0.0; 3], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
const TETS: [[usize; 4]; 1] = [[0, 1, 2, 3]];
const HEX_NODES: [[f64; 3]; 8] = [
    [0.0, 0.0, 0.0],
    [1.0, 0.0, 0.0],
    [1.0, 1.0, 0.0],
    [0.0, 1.0, 0.0],
    [0.0, 0.0, 1.0],
    [1.0, 0.0, 1.0],
    [1.0, 1.0, 1.0],
    [0.0, 1.0, 1.0],
];
const HEX: [[usize; 8]; 1] = [[0, 1, 2, 3, 4, 5, 6, 7]];
const HOT: f64 = 450.0;
const RHO: [f64; 1] = [6.0];
const PLASTIC: f64 = (3.0 * 8000.0 * 0.2 - 200.0) / (3.0 * 8000.0 + 800.0);

fn with_cx<R>(cancelled: bool, f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate = CancelGate::new_clock_free();
    if cancelled {
        gate.request();
    }
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        f(&Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 41,
                kernel_id: 179,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        ))
    })
}

fn family() -> AffineThermalFiniteJ2 {
    AffineThermalFiniteJ2::new(
        FiniteJ2ThermalParameters {
            reference_temperature_k: 300.0,
            temperature_range_k: [250.0, 650.0],
            reference: FiniteJ2Parameters {
                shear_modulus_pa: 8000.0,
                bulk_modulus_pa: 16000.0,
                yield_stress_pa: 200.0,
                hardening_modulus_pa: 800.0,
            },
            slopes: FiniteJ2TemperatureSlopes {
                shear_modulus_pa_k: -6.0,
                bulk_modulus_pa_k: -8.0,
                yield_stress_pa_k: -0.4,
                hardening_modulus_pa_k: -1.0,
            },
            logarithmic_expansion_per_k: 3e-4,
            logarithmic_expansion_slope_per_k2: 1e-7,
        },
        FiniteJ2Limits {
            minimum_jacobian: 0.25,
            maximum_elastic_log_strain: 0.7,
            maximum_plastic_log_strain: 2.0,
            maximum_accumulated_plastic_strain: 3.0,
        },
    )
    .unwrap()
}

fn controls() -> FiniteThermoplasticControls {
    FiniteThermoplasticControls {
        material: FinitePlasticStepControls {
            equilibrium_work_tolerance_j: 2e-6,
            maximum_endpoint_work_defect_j: 800.0,
            yield_tolerance_pa: 1e-5,
        },
        maximum_numerical_kinetic_loss_j: 1.0,
    }
}
fn impulse_controls() -> SolidImpulseControls {
    SolidImpulseControls {
        solve: CheckedCgConfig {
            tolerance: 1e-13,
            max_iterations: 64,
            max_corrections: 2,
        },
        maximum_workspace_components: 4096,
        linear_tolerance_kg_m_s: 1e-11,
        angular_tolerance_kg_m2_s: 1e-11,
        work_tolerance_j: 1e-11,
    }
}
fn tet_settings() -> HyperTetSettings {
    HyperTetSettings {
        load_steps: 1,
        relative_tolerance: 1e-11,
        ..Default::default()
    }
}
fn hex_settings() -> IncrementSettings {
    IncrementSettings {
        relative_tolerance: 1e-11,
        ..Default::default()
    }
}
fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:.16e} != {expected:.16e}"
    );
}
fn cross(x: [f64; 3], p: [f64; 3]) -> [f64; 3] {
    [
        x[1] * p[2] - x[2] * p[1],
        x[2] * p[0] - x[0] * p[2],
        x[0] * p[1] - x[1] * p[0],
    ]
}
fn shape(nodes: &[[f64; 3]], strain: f64, temperature: f64) -> Vec<(usize, f64)> {
    let d = temperature - 300.0;
    let thermal = 3e-4 * d + 0.5e-7 * d * d;
    let stretch = [
        (thermal + strain).exp(),
        (thermal - 0.5 * strain).exp(),
        (thermal - 0.5 * strain).exp(),
    ];
    (0..3 * nodes.len())
        .map(|dof| (dof, (stretch[dof % 3] - 1.0) * nodes[dof / 3][dof % 3]))
        .collect()
}

fn hot_tet(cx: &Cx<'_>) -> FiniteThermoplasticTetSolid<'static> {
    let mut solid = FiniteThermoplasticTetSolid::try_new(
        &TET_NODES,
        &TETS,
        vec![family()],
        &[300.0],
        TetAssemblyBudget::standard(),
        cx,
    )
    .unwrap();
    // Coaxial plastic loading, then unloading to the stress-free plastic shape
    // at a different temperature. Both hardening and thermal storage remain.
    for (strain, temperature) in [(0.2, 300.0), (PLASTIC, HOT)] {
        solid = solid
            .advance(
                &[temperature],
                &shape(&TET_NODES, strain, temperature),
                &[],
                [0.0; 3],
                0.05,
                tet_settings(),
                controls(),
                cx,
            )
            .unwrap()
            .solid;
    }
    close(solid.element_states()[0].alpha, PLASTIC, 1e-11);
    close(
        solid.free_energy_j(),
        0.5 * 650.0 * PLASTIC * PLASTIC / 6.0,
        1e-8,
    );
    close(
        solid.internal_energy_j(),
        0.5 * 1100.0 * PLASTIC * PLASTIC / 6.0,
        1e-8,
    );
    solid
}

fn hot_hex(formulation: FiniteHexFormulation, cx: &Cx<'_>) -> FiniteThermoplasticHexSolid<'static> {
    let mut solid = FiniteThermoplasticHexSolid::try_new(
        &HEX_NODES,
        &HEX,
        vec![family()],
        &[300.0],
        formulation,
        HexAssemblyBudget {
            maximum_tangent_entries: 768,
            ..Default::default()
        },
        cx,
    )
    .unwrap();
    for (strain, temperature) in [(0.2, 300.0), (PLASTIC, HOT)] {
        solid = solid
            .advance(
                &[temperature],
                &shape(&HEX_NODES, strain, temperature),
                &[],
                0.05,
                hex_settings(),
                controls().material,
                cx,
            )
            .unwrap()
            .solid;
    }
    for state in solid.quadrature_states() {
        close(state.alpha, PLASTIC, 1e-11);
    }
    close(solid.free_energy_j(), 0.5 * 650.0 * PLASTIC * PLASTIC, 1e-7);
    close(
        solid.internal_energy_j(),
        0.5 * 1100.0 * PLASTIC * PLASTIC,
        1e-7,
    );
    solid
}

#[test]
fn g1_hot_plastic_tet_impulse_retains_thermal_reference_storage_and_native_motion() {
    with_cx(false, |cx| {
        let solid = hot_tet(cx);
        assert!(
            solid
                .impulse_candidate(&[], &[(0, 0.1)], [0.0; 3], impulse_controls(), cx)
                .is_err()
        );
        let v0 = [0.2, -0.1, 0.05];
        let delta = [0.03, -0.02, 0.04];
        let expected_v: [f64; 3] = std::array::from_fn(|c| v0[c] + delta[c]);
        let solid = solid.with_inertia(&RHO, &[v0; 4], cx).unwrap();
        let impulses: Vec<_> = (0..12).map(|d| (d, 0.25 * delta[d % 3])).collect();
        let origin = [-0.2, 0.3, 0.4];
        let first_moment: [f64; 3] = std::array::from_fn(|c| {
            (0..4)
                .map(|i| 0.25 * (TET_NODES[i][c] + solid.displacement_m()[i][c] - origin[c]))
                .sum()
        });
        let event = solid
            .impulse_candidate(&[], &impulses, origin, impulse_controls(), cx)
            .unwrap();
        for velocity in event.solid.velocity_mps().unwrap() {
            for c in 0..3 {
                close(velocity[c], expected_v[c], 1e-13);
            }
        }
        for c in 0..3 {
            close(event.impulse.before.linear_kg_m_s[c], v0[c], 1e-13);
            close(event.impulse.after.linear_kg_m_s[c], expected_v[c], 1e-13);
            close(
                event.impulse.applied.angular_kg_m2_s[c],
                cross(first_moment, delta)[c],
                1e-13,
            );
            close(event.impulse.angular_residual_kg_m2_s[c], 0.0, 1e-13);
        }
        close(
            event.impulse.kinetic_change_j,
            (0..3).map(|c| (v0[c] + 0.5 * delta[c]) * delta[c]).sum(),
            1e-14,
        );
        close(event.impulse.work_residual_j, 0.0, 1e-14);
        assert_eq!(event.solid.temperatures_k(), &[HOT]);
        assert!(std::ptr::eq(
            event.solid.families().as_ptr(),
            solid.families().as_ptr()
        ));
        assert_eq!(event.solid.displacement_m(), solid.displacement_m());
        assert_eq!(
            event.solid.element_deformation(),
            solid.element_deformation()
        );
        assert_eq!(event.solid.element_states(), solid.element_states());
        assert_eq!(
            event.solid.element_first_piola_pa(),
            solid.element_first_piola_pa()
        );
        assert_eq!(event.solid.free_energy_j(), solid.free_energy_j());
        assert_eq!(event.solid.internal_energy_j(), solid.internal_energy_j());
        assert_eq!(event.solid.reference_mass_kg(), Some(1.0));
        assert_eq!(event.solid.reference_density_kg_m3(), Some(RHO.as_slice()));
        assert_eq!(solid.velocity_mps().unwrap(), [v0; 4]);
        with_cx(true, |stopped| {
            assert!(
                solid
                    .impulse_candidate(&[], &impulses, origin, impulse_controls(), stopped)
                    .is_err()
            );
        });
        assert_eq!(solid.velocity_mps().unwrap(), [v0; 4]);
        let dt = 0.02;
        let next = event
            .solid
            .advance(
                &[HOT],
                &[],
                &[],
                [0.0; 3],
                dt,
                tet_settings(),
                controls(),
                cx,
            )
            .unwrap();
        assert_eq!(next.solid.element_states(), solid.element_states());
        assert_eq!(next.solid.temperatures_k(), &[HOT]);
        assert_eq!(next.solid.families(), solid.families());
        for i in 0..4 {
            for c in 0..3 {
                close(
                    next.solid.displacement_m()[i][c],
                    solid.displacement_m()[i][c] + dt * expected_v[c],
                    1e-10,
                );
                close(
                    next.solid.velocity_mps().unwrap()[i][c],
                    expected_v[c],
                    1e-8,
                );
            }
        }
        close(
            next.solid.internal_energy_j(),
            solid.internal_energy_j(),
            1e-7,
        );
        close(next.energy.caloric_heat_j, 0.0, 1e-7);
        close(next.energy.plastic_dissipation_j, 0.0, 0.0);
        close(
            next.inertia.as_ref().unwrap().kinetic_before_j,
            event.impulse.kinetic_after_j,
            1e-13,
        );
        close(next.energy.balance_error_j, 0.0, 2e-6);
    });
}

// Exact tensor-product integral on the original unit brick; current hot volume
// and the mean-dilatation constitutive projection do not enter reference mass.
fn hex_mass(i: usize, j: usize) -> f64 {
    RHO[0]
        * (0..3)
            .map(|c| {
                if HEX_NODES[i][c] == HEX_NODES[j][c] {
                    1.0 / 3.0
                } else {
                    1.0 / 6.0
                }
            })
            .product::<f64>()
}

#[test]
fn g1_hot_hex_impulse_keeps_all_quadrature_memory_and_moving_support_work() {
    with_cx(false, |cx| {
        for formulation in [
            FiniteHexFormulation::Standard,
            FiniteHexFormulation::MeanDilatation,
        ] {
            let v0 = [[0.0, 0.1, 0.0]; 8];
            let motion =
                FiniteThermoplasticHexMotion::try_new(hot_hex(formulation, cx), &RHO, &v0, cx)
                    .unwrap();
            let old = motion.solid();
            let held: Vec<_> = (0..24).filter(|&d| d != 4).collect();
            let origin = [0.4, -0.3, 0.2];
            let impulses = [(1, 0.04), (4, 0.03)];
            let delta = 0.03 / hex_mass(1, 1);
            let event = motion
                .impulse_candidate(&held, &impulses, origin, impulse_controls(), cx)
                .unwrap();
            let mut support = [0.0; 24];
            for &(d, value) in &event.impulse.support_impulses_kg_m_s {
                support[d] = value;
            }
            let mut total_torque = [0.0; 3];
            for i in 0..8 {
                let dp = [0.0, hex_mass(i, 1) * delta, 0.0];
                let torque = cross(
                    std::array::from_fn(|c| {
                        HEX_NODES[i][c] + old.displacement_m()[i][c] - origin[c]
                    }),
                    dp,
                );
                for c in 0..3 {
                    close(
                        event.motion.velocity_mps()[i][c],
                        v0[i][c] + if i == 1 && c == 1 { delta } else { 0.0 },
                        1e-13,
                    );
                    let expected_support = if i == 1 {
                        0.0
                    } else {
                        dp[c] - if i == 0 && c == 1 { 0.04 } else { 0.0 }
                    };
                    close(support[3 * i + c], expected_support, 1e-13);
                    total_torque[c] += torque[c];
                }
            }
            for c in 0..3 {
                close(
                    event.impulse.applied.angular_kg_m2_s[c]
                        + event.impulse.support.angular_kg_m2_s[c],
                    total_torque[c],
                    1e-13,
                );
                close(event.impulse.angular_residual_kg_m2_s[c], 0.0, 1e-13);
            }
            close(
                event.impulse.kinetic_change_j,
                0.1 * (RHO[0] / 8.0) * delta + 0.5 * hex_mass(1, 1) * delta * delta,
                1e-14,
            );
            close(
                event.impulse.applied_work_j,
                0.04 * 0.1 + 0.03 * (0.1 + 0.5 * delta),
                1e-14,
            );
            close(
                event.impulse.support_work_j,
                0.1 * ((RHO[0] / 8.0 - hex_mass(1, 1)) * delta - 0.04),
                1e-14,
            );
            close(event.impulse.work_residual_j, 0.0, 1e-14);
            let new = event.motion.solid();
            assert_eq!(new.temperatures_k(), &[HOT]);
            assert!(std::ptr::eq(
                new.families().as_ptr(),
                old.families().as_ptr()
            ));
            assert_eq!(new.formulation(), formulation);
            assert_eq!(new.displacement_m(), old.displacement_m());
            assert_eq!(new.quadrature_deformation(), old.quadrature_deformation());
            assert_eq!(
                new.quadrature_assumed_deformation(),
                old.quadrature_assumed_deformation()
            );
            assert_eq!(new.quadrature_states(), old.quadrature_states());
            assert_eq!(
                new.quadrature_first_piola_pa(),
                old.quadrature_first_piola_pa()
            );
            assert_eq!(new.free_energy_j(), old.free_energy_j());
            assert_eq!(new.internal_energy_j(), old.internal_energy_j());
            close(event.motion.reference_mass_kg(), RHO[0], 1e-13);
            assert_eq!(event.motion.reference_density_kg_m3(), RHO);
            assert_eq!(motion.velocity_mps(), v0);
            let mut limited = impulse_controls();
            limited.maximum_workspace_components = 0;
            assert!(
                motion
                    .impulse_candidate(&held, &impulses, origin, limited, cx)
                    .is_err()
            );
            assert_eq!(motion.velocity_mps(), v0);
            let dt = 0.001;
            let supports: Vec<_> = held
                .iter()
                .map(|&d| {
                    (
                        d,
                        old.displacement_m()[d / 3][d % 3] + dt * v0[d / 3][d % 3],
                    )
                })
                .collect();
            let next = event
                .motion
                .advance(&[HOT], &supports, &[], dt, hex_settings(), controls(), cx)
                .unwrap();
            close(
                next.inertia.kinetic_before_j,
                event.impulse.kinetic_after_j,
                1e-13,
            );
            assert_eq!(
                next.motion.solid().quadrature_states(),
                old.quadrature_states()
            );
            assert_eq!(next.motion.solid().temperatures_k(), &[HOT]);
            assert_eq!(next.motion.solid().formulation(), formulation);
            close(
                next.material_energy.internal_energy_change_j,
                next.motion.solid().internal_energy_j() - old.internal_energy_j(),
                1e-7,
            );
            close(next.balance_error_j, 0.0, 2e-6);
            close(next.motion.reference_mass_kg(), RHO[0], 1e-13);
        }
    });
}
