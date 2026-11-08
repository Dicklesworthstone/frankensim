//! G1/G3/G4: thermally biased transient memories retain their physical DC state.
//! Uniform fields give independent capacitor, Debye and Maxwell oracles. These
//! authored SI sources make no experimental or nonlinear-material claim.
use super::*;
use crate::harmonic::thermal_memory_source;
use crate::transient_load::resistor;
use fs_couple::piezoelectric::spatial::transient::load::ElectrodeDrive as Drive;

const COUPLING: f64 = 12.0;
const EXPANSION: f64 = 2e-6 * 10.0;
const POLARIZATION: f64 = -3e-5 * 10.0;
const INCREMENT: f64 = 0.8 * EPSILON[2];
const FRACTION: f64 = 0.5;
const TAU_P: f64 = 1e-6;
const TAU_M: f64 = 2e-7;
const AREA: f64 = LENGTH * LENGTH;
const VOLUME: f64 = AREA * LENGTH;
const MEMORY_BUDGET: usize = 256;

fn axial(cube: &Cube, strain: f64) -> Vec<[f64; 3]> {
    cube.nodes
        .iter()
        .map(|p| [0.0, 0.0, strain * p[2]])
        .collect()
}

fn fixed(displacement: &[[f64; 3]]) -> Vec<(usize, f64)> {
    displacement
        .iter()
        .enumerate()
        .flat_map(|(node, value)| (0..3).map(move |axis| (3 * node + axis, value[axis])))
        .collect()
}

#[test]
fn g1_g3_hot_relaxed_voltage_and_charge_states_keep_polarization_without_startup_heat() {
    with_cx(|cx, _| {
        let cube = cube(1);
        let n = cube.nodes.len();
        let (sample, loss) = thermal_memory_source([0.0, f64::MAX]);
        let model = Model::assemble_thermal_materials(
            &cube.nodes,
            &cube.tets,
            vec![sample],
            FRAME,
            budget(),
            cx,
        )
        .unwrap();
        let (zero, charge) = zero_source(n);
        let strain = 3e-5;
        let displacement = axial(&cube, strain);
        let supports = fixed(&displacement);
        let voltage = 4.0;
        let field = -voltage / LENGTH;
        let electric_displacement =
            COUPLING * strain + (EPSILON[2] + INCREMENT) * field + POLARIZATION;
        let terminal_charge = -AREA * electric_displacement;
        let stress = STIFFNESS[2] * (strain - EXPANSION) - COUPLING * field;
        let memory = 0.5 * VOLUME * INCREMENT * field * field;
        let stored = 0.5
            * VOLUME
            * (STIFFNESS[2] * (strain - EXPANSION).powi(2)
                + (EPSILON[2] + INCREMENT) * field * field);
        // Uniform/per-cell laws and a large common voltage must recover the
        // same physical state, including the affine clamped-pyroelectric charge.
        for per_cell in [false, true] {
            for (gauge, charged) in [(0.0, false), (2.0_f64.powi(40), true)] {
                let top = if charged {
                    Condition::Charge(terminal_charge)
                } else {
                    Condition::Voltage(gauge + voltage)
                };
                let laws = if per_cell {
                    vec![loss.clone(); cube.tets.len()]
                } else {
                    vec![loss.clone()]
                };
                let body = Body::new_with_equilibrium_relaxation(
                    &model,
                    &supports,
                    &electrodes(&cube, gauge, top),
                    &displacement,
                    &zero,
                    &charge,
                    laws,
                    MEMORY_BUDGET,
                    controls(),
                    cx,
                )
                .unwrap();
                let state = body.state();
                assert_eq!(state.displacement_m(), &displacement);
                assert_eq!(state.velocity_m_s(), &zero);
                near(state.electrodes()[1].voltage_v, gauge + voltage, 2e-8, 0.0);
                near(state.electrodes()[1].charge_c, terminal_charge, 2e-18, 3e-8);
                near(state.material_memory_energy_j(), memory, 2e-20, 3e-8);
                near(state.stored_energy_j(), stored, 2e-17, 3e-8);
                let initial = body.material_relaxation_steps().unwrap();
                for step in initial {
                    assert_eq!(step.mechanical_heat_j, 0.0);
                    assert_eq!(step.dielectric_heat_j, 0.0);
                    assert_eq!(step.memory_energy_change_j, 0.0);
                    near(
                        step.memory_energy_j[0],
                        memory / cube.tets.len() as f64,
                        2e-20,
                        3e-8,
                    );
                    assert_eq!(step.memory_energy_j[0], step.memory_energy_j[1]);
                }
                let candidate = body
                    .step_candidate(
                        0.7 * TAU_P,
                        &[Condition::Voltage(gauge), top],
                        &zero,
                        &charge,
                        controls(),
                        cx,
                    )
                    .unwrap();
                let state = candidate.body().state();
                for cell in state.cells() {
                    near(cell.strain_mandel[2], strain, 2e-12, 3e-8);
                    near(cell.response.electric_field_v_m[2], field, 2e-5, 3e-8);
                    near(cell.response.stress_mandel_pa[2], stress, 0.05, 3e-8);
                    near(
                        cell.response.displacement_c_m2[2],
                        electric_displacement,
                        2e-12,
                        3e-8,
                    );
                }
                near(state.material_memory_energy_j(), memory, 2e-20, 3e-8);
                near(state.total_energy_j(), stored, 2e-17, 3e-8);
                near(candidate.step().material_heat_j, 0.0, 1e-24, 0.0);
                near(candidate.step().electrical_work_j, 0.0, 1e-20, 0.0);
                near(candidate.step().energy_residual_j, 0.0, 1e-20, 0.0);
                near(
                    cube.top
                        .iter()
                        .map(|&i| candidate.midpoint_support_force_n()[i][2])
                        .sum(),
                    AREA * stress,
                    1e-6,
                    3e-8,
                );
                assert_eq!(body.state().time_s(), 0.0);
            }
        }
    });
}

#[test]
fn g1_g4_relaxed_initialization_uses_its_own_field_limit_and_source_budget() {
    with_cx(|cx, gate| {
        let cube = cube(1);
        let n = cube.nodes.len();
        let (sample, loss) = thermal_memory_source([0.0, f64::MAX]);
        let model = Model::assemble_thermal_materials(
            &cube.nodes,
            &cube.tets,
            vec![sample],
            FRAME,
            budget(),
            cx,
        )
        .unwrap();
        let (zero, charge) = zero_source(n);
        let supports = fixed(&zero);
        let field = 1.5e6;
        let terminal_charge = -AREA * ((EPSILON[2] + INCREMENT) * field + POLARIZATION);
        let terminals = electrodes(&cube, 0.0, Condition::Charge(terminal_charge));
        assert!(field < controls().maximum_electric_field_v_m);
        assert!(
            field * (EPSILON[2] + INCREMENT) / EPSILON[2] > controls().maximum_electric_field_v_m
        );
        // An instantaneous-epsilon preliminary solve would incorrectly reject
        // this fully relaxed charged state. The direct constructor must not do it.
        assert!(
            Body::new(
                &model,
                &supports,
                &terminals,
                &zero,
                &zero,
                &charge,
                controls(),
                cx
            )
            .is_err()
        );
        let variables = loss.internal_variables() * cube.tets.len();
        assert!(
            Body::new_with_equilibrium_relaxation(
                &model,
                &supports,
                &terminals,
                &zero,
                &zero,
                &charge,
                vec![loss.clone()],
                variables - 1,
                controls(),
                cx,
            )
            .is_err()
        );
        let body = Body::new_with_equilibrium_relaxation(
            &model,
            &supports,
            &terminals,
            &zero,
            &zero,
            &charge,
            vec![loss.clone()],
            variables,
            controls(),
            cx,
        )
        .unwrap();
        near(
            body.state().electrodes()[1].voltage_v,
            -LENGTH * field,
            2e-5,
            3e-8,
        );
        near(
            body.state().material_memory_energy_j(),
            0.5 * VOLUME * INCREMENT * field * field,
            2e-16,
            3e-8,
        );
        let (sample, finite_band) = thermal_memory_source([0.0, 8e6]);
        let bounded = Model::assemble_thermal_materials(
            &cube.nodes,
            &cube.tets,
            vec![sample],
            FRAME,
            budget(),
            cx,
        )
        .unwrap();
        assert!(
            Body::new_with_equilibrium_relaxation(
                &bounded,
                &supports,
                &terminals,
                &zero,
                &zero,
                &charge,
                vec![finite_band],
                MEMORY_BUDGET,
                controls(),
                cx,
            )
            .is_err()
        );
        gate.request();
        assert!(matches!(
            Body::new_with_equilibrium_relaxation(
                &model,
                &supports,
                &terminals,
                &zero,
                &zero,
                &charge,
                vec![loss],
                MEMORY_BUDGET,
                controls(),
                cx,
            ),
            Err(TetPiezoelectricError::Cancelled)
        ));
    });
}

#[test]
fn g1_hot_equilibrated_polarization_discharges_through_the_native_resistor_solve() {
    with_cx(|cx, _| {
        let cube = cube(1);
        let n = cube.nodes.len();
        let (sample, loss) = thermal_memory_source([0.0, f64::MAX]);
        let model = Model::assemble_thermal_materials(
            &cube.nodes,
            &cube.tets,
            vec![sample],
            FRAME,
            budget(),
            cx,
        )
        .unwrap();
        let (zero, charge) = zero_source(n);
        let c = EPSILON[2] * LENGTH;
        let cp = INCREMENT * LENGTH;
        let q_reference = -AREA * POLARIZATION;
        let resistance = 5e4;
        let voltage = 20.0;
        let mut body = Body::new_with_equilibrium_relaxation(
            &model,
            &fixed(&zero),
            &electrodes(
                &cube,
                0.0,
                Condition::Charge((c + cp) * voltage + q_reference),
            ),
            &zero,
            &zero,
            &charge,
            vec![loss],
            MEMORY_BUDGET,
            controls(),
            cx,
        )
        .unwrap();
        let elastic = 0.5 * VOLUME * STIFFNESS[2] * EXPANSION * EXPANSION;
        let initial = elastic + 0.5 * (c + cp) * voltage * voltage;
        let dt = 0.4 * TAU_P;
        let h = 0.5 * dt;
        let (mut v, mut p) = (voltage, cp * voltage);
        let (mut resistor_heat, mut ceramic_heat) = (0.0, 0.0);
        for _ in 0..24 {
            // C Vdot + pdot = -V/R; tau pdot+p=Cp V. The fixed
            // pyroelectric charge shifts Q but cannot create current or heat.
            let vm = (c * v / h + p / (TAU_P + h)) / (c / h + cp / (TAU_P + h) + 1.0 / resistance);
            let pm = (TAU_P * p + h * cp * vm) / (TAU_P + h);
            let next_v = 2.0 * vm - v;
            let next_p = 2.0 * pm - p;
            let dielectric_heat = TAU_P * (next_p - p).powi(2) / (cp * dt);
            let joule_heat = dt * vm * vm / resistance;
            let memory = 0.5 * next_p * next_p / cp;
            let (candidate, loaded) = body
                .step_with_loads_candidate(
                    dt,
                    &[Drive::Voltage(0.0), Drive::Current(0.0)],
                    &[resistor(resistance, [1, 0])],
                    &zero,
                    &charge,
                    controls(),
                    cx,
                )
                .unwrap();
            let state = candidate.body().state();
            near(state.electrodes()[1].voltage_v, next_v, 2e-9, 3e-8);
            near(
                state.electrodes()[1].charge_c,
                c * next_v + next_p + q_reference,
                2e-18,
                3e-8,
            );
            near(state.material_memory_energy_j(), memory, 2e-20, 3e-8);
            near(
                state.stored_energy_j(),
                elastic + 0.5 * c * next_v * next_v + memory,
                2e-17,
                3e-8,
            );
            near(
                candidate.step().dielectric_material_heat_j,
                dielectric_heat,
                2e-20,
                3e-8,
            );
            near(candidate.step().mechanical_material_heat_j, 0.0, 1e-24, 0.0);
            near(loaded.joule_heat_j, joule_heat, 2e-20, 3e-8);
            near(loaded.electrical_source_work_j, 0.0, 1e-20, 0.0);
            near(loaded.energy_residual_j, 0.0, 1e-19, 0.0);
            resistor_heat += joule_heat;
            ceramic_heat += dielectric_heat;
            body = candidate.into_parts().0;
            (v, p) = (next_v, next_p);
        }
        assert!(resistor_heat > 0.0);
        assert!(ceramic_heat > 0.0);
        near(
            body.state().total_energy_j() + resistor_heat + ceramic_heat,
            initial,
            2e-17,
            3e-8,
        );
    });
}

#[test]
fn g1_fixed_thermal_preload_and_maxwell_support_motion_keep_separate_storage() {
    with_cx(|cx, _| {
        let cube = cube(1);
        let n = cube.nodes.len();
        let (sample, loss) = thermal_memory_source([0.0, f64::MAX]);
        let model = Model::assemble_thermal_materials(
            &cube.nodes,
            &cube.tets,
            vec![sample],
            FRAME,
            budget(),
            cx,
        )
        .unwrap();
        let (zero, charge) = zero_source(n);
        let mut body = Body::new(
            &model,
            &fixed(&zero),
            &electrodes(&cube, 0.0, Condition::Voltage(0.0)),
            &zero,
            &zero,
            &charge,
            controls(),
            cx,
        )
        .unwrap()
        .with_relaxation(vec![loss], MEMORY_BUDGET, controls(), cx)
        .unwrap();
        assert_eq!(body.state().material_memory_energy_j(), 0.0);
        let dt = 0.2 * TAU_M;
        let a = TAU_M / (TAU_M + 0.5 * dt);
        let branch_stiffness = FRACTION * STIFFNESS[2];
        let (mut strain, mut rate, mut stress_memory) = (0.0, 0.0, 0.0);
        for next_rate in [150.0, 150.0, 0.0, -75.0, 0.0, 0.0] {
            let increment = 0.5 * dt * (rate + next_rate);
            let next_strain = strain + increment;
            // At fixed T, delta(S-alpha deltaT)=deltaS. A thermal preload
            // is already in the base energy; it is not a fresh Maxwell drive.
            let mid_memory = a * (stress_memory + 0.5 * branch_stiffness * increment);
            let next_memory = 2.0 * mid_memory - stress_memory;
            let memory = 0.5 * VOLUME * next_memory * next_memory / branch_stiffness;
            let heat = VOLUME * dt * mid_memory * mid_memory / (branch_stiffness * TAU_M);
            let kinetic = 7500.0 * VOLUME * LENGTH * LENGTH * next_rate * next_rate / 6.0;
            let midpoint_stress =
                STIFFNESS[2] * (0.5 * (strain + next_strain) - EXPANSION) + mid_memory;
            let before_kinetic = 7500.0 * VOLUME * LENGTH * LENGTH * rate * rate / 6.0;
            let work = VOLUME * midpoint_stress * increment + kinetic - before_kinetic;
            let candidate = body
                .step_with_support_motion_candidate(
                    dt,
                    &fixed(&axial(&cube, next_rate)),
                    &[Condition::Voltage(0.0); 2],
                    &zero,
                    &charge,
                    controls(),
                    cx,
                )
                .unwrap();
            let state = candidate.body().state();
            for cell in state.cells() {
                near(cell.strain_mandel[2], next_strain, 2e-12, 3e-8);
                near(
                    cell.response.stress_mandel_pa[2],
                    STIFFNESS[2] * (next_strain - EXPANSION) + next_memory,
                    0.05,
                    3e-8,
                );
                near(
                    cell.response.displacement_c_m2[2],
                    COUPLING * next_strain + POLARIZATION,
                    2e-12,
                    3e-8,
                );
            }
            near(state.material_memory_energy_j(), memory, 2e-20, 3e-8);
            near(
                state.stored_energy_j(),
                0.5 * VOLUME * STIFFNESS[2] * (next_strain - EXPANSION).powi(2) + memory,
                2e-17,
                3e-8,
            );
            near(state.kinetic_energy_j(), kinetic, 2e-20, 3e-8);
            near(
                candidate.step().mechanical_material_heat_j,
                heat,
                2e-20,
                3e-8,
            );
            near(candidate.step().dielectric_material_heat_j, 0.0, 1e-24, 0.0);
            near(candidate.step().electrical_work_j, 0.0, 1e-20, 0.0);
            near(candidate.step().support_work_j, work, 2e-17, 3e-8);
            near(candidate.step().energy_residual_j, 0.0, 1e-19, 0.0);
            body = candidate.into_parts().0;
            (strain, rate, stress_memory) = (next_strain, next_rate, next_memory);
        }
    });
}
