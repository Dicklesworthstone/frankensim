//! G1/G3/G4: independent peak-phasor circuit and mechanical response oracles.
use super::*;
use fs_couple::piezoelectric::spatial::harmonic::{
    HarmonicElectrode, HarmonicElectrodeCondition as HCondition, HarmonicElectrodeLoad as HLoad,
    HarmonicPiezoelectricControls as HControls,
};
use fs_math::c64::C64;

#[path = "harmonic/relaxation.rs"]
pub(super) mod relaxation;

#[path = "harmonic/acoustic.rs"]
mod acoustic;

#[path = "harmonic/operating_point.rs"]
mod operating_point;

const ZERO: C64 = C64::new(0.0, 0.0);

fn hc() -> HControls {
    HControls {
        field: controls(),
        power_tolerance_w: 1e-9,
        maximum_real_coordinates: 512,
    }
}
fn scale(z: C64, a: f64) -> C64 {
    C64::new(a * z.re, a * z.im)
}
fn add(a: C64, b: C64) -> C64 {
    C64::new(a.re + b.re, a.im + b.im)
}
fn mul(a: C64, b: C64) -> C64 {
    C64::new(a.re * b.re - a.im * b.im, a.re * b.im + a.im * b.re)
}
fn div(a: C64, b: C64) -> C64 {
    scale(mul(a, b.conj()), 1.0 / b.norm_sq())
}
fn cnear(a: C64, b: C64, absolute: f64, relative: f64) {
    let error = (a.re - b.re).hypot(a.im - b.im);
    let magnitude = a.re.hypot(a.im).max(b.re.hypot(b.im));
    assert!(
        error <= absolute + relative * magnitude,
        "{a:?} != {b:?}, error {error:e}"
    );
}
fn helectrodes(cube: &Cube, reference: C64, top: HCondition) -> [HarmonicElectrode<'_>; 2] {
    [
        HarmonicElectrode {
            nodes: &cube.bottom,
            condition: HCondition::Voltage(reference),
        },
        HarmonicElectrode {
            nodes: &cube.top,
            condition: top,
        },
    ]
}
fn complex_supports(fixed: &[(usize, f64)]) -> Vec<(usize, C64)> {
    fixed
        .iter()
        .map(|&(dof, value)| (dof, C64::new(value, 0.0)))
        .collect()
}

#[test]
fn g1_clamped_current_drive_has_exact_rc_phase_charge_and_average_power() {
    with_cx(|cx, _| {
        let cube = cube(1);
        let n = cube.nodes.len();
        let materials = [material([0.0, 0.0, 12.0, 0.0, 0.0, 0.0])];
        let model =
            Model::assemble(&cube.nodes, &cube.tets, &materials, FRAME, budget(), cx).unwrap();
        let fixed: Vec<_> = (0..3 * n).map(|dof| (dof, ZERO)).collect();
        let resistance = 1e5;
        let capacitance = EPSILON[2] * LENGTH;
        let omega = 1.0 / (resistance * capacitance);
        let current = C64::new(3e-5, -4e-5);
        let load = HLoad::Resistive(super::transient_load::resistor(resistance, [1, 0]));
        let response = model
            .solve_harmonic(
                omega,
                &fixed,
                &helectrodes(&cube, ZERO, HCondition::Current(current)),
                &[load],
                &vec![[ZERO; 3]; n],
                &vec![ZERO; n],
                hc(),
                cx,
            )
            .unwrap();
        // With exp(-i omega t), the clamped ceramic is a parallel capacitor:
        // I=(1/R-i omega C)V. Peak amplitudes contribute one half to power.
        let voltage = div(current, C64::new(1.0 / resistance, -omega * capacitance));
        cnear(response.electrodes[1].voltage_v, voltage, 1e-10, 2e-9);
        cnear(
            response.electrodes[1].charge_c,
            scale(voltage, capacitance),
            1e-20,
            2e-9,
        );
        cnear(
            response.electrodes[1].source_current_a,
            current,
            1e-14,
            2e-9,
        );
        cnear(
            response.electrodes[0].source_current_a,
            scale(current, -1.0),
            1e-14,
            2e-9,
        );
        let heat = 0.5 * voltage.norm_sq() / resistance;
        near(response.average_load_heat_w, heat, 1e-14, 2e-9);
        near(response.average_electrical_power_w, heat, 1e-14, 2e-9);
        near(response.average_mechanical_power_w, 0.0, 1e-20, 0.0);
        near(
            response.average_stored_energy_j,
            0.25 * capacitance * voltage.norm_sq(),
            1e-20,
            2e-9,
        );
        assert_eq!(response.average_kinetic_energy_j, 0.0);
        near(response.power_balance_residual_w, 0.0, 1e-13, 0.0);
        for cell in &response.cells {
            cnear(
                cell.electric_field_v_m[2],
                scale(voltage, -1.0 / LENGTH),
                1e-7,
                2e-9,
            );
            cnear(
                cell.displacement_c_m2[2],
                scale(voltage, -EPSILON[2] / LENGTH),
                1e-17,
                2e-9,
            );
        }
    });
}

#[test]
fn g1_clamped_reactive_shunts_match_lc_and_rl_parallel_admittance() {
    use super::transient_load::reactive::{capacitor, inductor};
    use fs_couple::piezoelectric::spatial::transient::load::reactive::ReactiveElectrodeLoad;

    with_cx(|cx, _| {
        let cube = cube(1);
        let n = cube.nodes.len();
        let materials = [material([0.0, 0.0, 12.0, 0.0, 0.0, 0.0])];
        let model =
            Model::assemble(&cube.nodes, &cube.tets, &materials, FRAME, budget(), cx).unwrap();
        let fixed: Vec<_> = (0..3 * n).map(|dof| (dof, ZERO)).collect();
        let ceramic_c = EPSILON[2] * LENGTH;
        let extra_c = 2.0 * ceramic_c;
        let l = 2e-3;
        let current = C64::new(1e-5, 2e-5);
        // Stay away from the lossless electrical pole while testing both
        // inductor and external-capacitor phasor signs independently.
        let omega = 0.7 / (l * (ceramic_c + extra_c)).sqrt();
        for resistance in [None, Some(1e4)] {
            let loads = [
                HLoad::Reactive(ReactiveElectrodeLoad::SeriesRl {
                    electrodes: [1, 0],
                    inductor: inductor(l),
                    conductor: resistance
                        .map(|r| super::transient_load::resistor(r, [1, 0]).conductor),
                }),
                HLoad::Reactive(ReactiveElectrodeLoad::Capacitor {
                    electrodes: [1, 0],
                    capacitor: capacitor(extra_c),
                }),
            ];
            let winding_admittance = div(
                C64::new(1.0, 0.0),
                C64::new(resistance.unwrap_or(0.0), -omega * l),
            );
            let admittance = add(
                winding_admittance,
                C64::new(0.0, -omega * (ceramic_c + extra_c)),
            );
            let voltage = div(current, admittance);
            let winding_current = mul(winding_admittance, voltage);
            let capacitor_current = mul(C64::new(0.0, -omega * extra_c), voltage);
            let response = model
                .solve_harmonic(
                    omega,
                    &fixed,
                    &helectrodes(&cube, ZERO, HCondition::Current(current)),
                    &loads,
                    &vec![[ZERO; 3]; n],
                    &vec![ZERO; n],
                    hc(),
                    cx,
                )
                .unwrap();
            cnear(response.electrodes[1].voltage_v, voltage, 1e-10, 2e-9);
            cnear(response.loads[0].current_a, winding_current, 1e-14, 2e-9);
            cnear(response.loads[1].current_a, capacitor_current, 1e-14, 2e-9);
            let heat = 0.5 * resistance.unwrap_or(0.0) * winding_current.norm_sq();
            near(response.average_load_heat_w, heat, 1e-14, 2e-9);
            near(response.average_electrical_power_w, heat, 1e-14, 2e-9);
            assert_eq!(response.loads[1].average_heat_w, 0.0);
            if resistance.is_none() {
                assert_eq!(response.loads[0].average_heat_w, 0.0);
            }
        }
    });
}

#[test]
fn g1_low_frequency_direct_and_converse_responses_recover_open_short_limits() {
    with_cx(|cx, _| {
        let cube = cube(1);
        let n = cube.nodes.len();
        let e = 10.0;
        let materials = [material([0.0, 0.0, e, 0.0, 0.0, 0.0])];
        let model =
            Model::assemble(&cube.nodes, &cube.tets, &materials, FRAME, budget(), cx).unwrap();
        let fixed = complex_supports(&cube.supports);
        let stress = 1e6;
        let force: Vec<_> = axial_traction(&cube, stress)
            .into_iter()
            .map(|f| f.map(|value| C64::new(value, 0.0)))
            .collect();
        // Inertia is retained, but omega L / c is below 4e-6 here.
        let omega = 10.0;
        for (condition, modulus) in [
            (HCondition::Voltage(ZERO), STIFFNESS[2]),
            (HCondition::Current(ZERO), STIFFNESS[2] + e * e / EPSILON[2]),
        ] {
            let response = model
                .solve_harmonic(
                    omega,
                    &fixed,
                    &helectrodes(&cube, ZERO, condition),
                    &[],
                    &force,
                    &vec![ZERO; n],
                    hc(),
                    cx,
                )
                .unwrap();
            let strain = stress / modulus;
            for &node in &cube.top {
                cnear(
                    response.displacement_m[node][2],
                    C64::new(strain * LENGTH, 0.0),
                    2e-15,
                    2e-7,
                );
            }
            match condition {
                HCondition::Current(_) => {
                    cnear(
                        response.electrodes[1].voltage_v,
                        C64::new(e * strain * LENGTH / EPSILON[2], 0.0),
                        2e-7,
                        2e-7,
                    );
                    cnear(response.electrodes[1].charge_c, ZERO, 1e-18, 0.0);
                }
                HCondition::Voltage(_) => {
                    cnear(
                        response.electrodes[1].charge_c,
                        C64::new(-e * strain * LENGTH * LENGTH, 0.0),
                        1e-18,
                        2e-7,
                    );
                }
            }
        }
        let voltage = C64::new(-10.0, 6.0);
        let actuator = model
            .solve_harmonic(
                omega,
                &fixed,
                &helectrodes(&cube, ZERO, HCondition::Voltage(voltage)),
                &[],
                &vec![[ZERO; 3]; n],
                &vec![ZERO; n],
                hc(),
                cx,
            )
            .unwrap();
        for &node in &cube.top {
            cnear(
                actuator.displacement_m[node][2],
                scale(voltage, -e / STIFFNESS[2]),
                2e-15,
                2e-7,
            );
        }
        let free_capacitance = (EPSILON[2] + e * e / STIFFNESS[2]) * LENGTH;
        cnear(
            actuator.electrodes[1].charge_c,
            scale(voltage, free_capacitance),
            1e-18,
            2e-7,
        );
    });
}

#[test]
fn g1_loaded_mechanical_resonance_matches_scalar_dynamic_stiffness() {
    with_cx(|cx, _| {
        let cube = cube(1);
        let n = cube.nodes.len();
        let materials = [material([0.0, 0.0, 12.0, 0.0, 0.0, 0.0])];
        let model =
            Model::assemble(&cube.nodes, &cube.tets, &materials, FRAME, budget(), cx).unwrap();
        let node = *cube.top.last().unwrap();
        let dof = 3 * node + 2;
        let fixed: Vec<_> = (0..3 * n)
            .filter(|&i| i != dof)
            .map(|i| (i, ZERO))
            .collect();
        let diagonal = |matrix: &fs_sparse::Csr| {
            let (cols, values) = matrix.row(dof);
            values[cols.iter().position(|&i| i == dof).unwrap()]
        };
        let mass = diagonal(&model.solid().mass);
        let stiffness = diagonal(&model.solid().stiffness);
        let (cols, values) = model.coupling().row(dof);
        let coupling = cols
            .iter()
            .zip(values)
            .filter(|(col, _)| cube.top.contains(col))
            .map(|(_, value)| value)
            .sum::<f64>();
        let capacitance = EPSILON[2] * LENGTH;
        let resonance = (stiffness / mass).sqrt();
        let resistance = 1.0 / (resonance * capacitance);
        let loads = [HLoad::Resistive(super::transient_load::resistor(
            resistance,
            [1, 0],
        ))];
        let force = C64::new(0.02, 0.01);
        let mut forces = vec![[ZERO; 3]; n];
        forces[node][2] = force;
        assert!(coupling.abs() > 0.0);
        for ratio in [0.4, 1.0, 1.2] {
            let omega = ratio * resonance;
            let admittance = C64::new(1.0 / resistance, -omega * capacitance);
            // The native G sign gives q=C V-G u and (K-w^2 M)u+G V=F.
            // Eliminating I=0 yields V/u=-i w G/(1/R-i w C).
            let voltage_per_motion = div(C64::new(0.0, -omega * coupling), admittance);
            let dynamic_stiffness = add(
                C64::new(stiffness - omega * omega * mass, 0.0),
                scale(voltage_per_motion, coupling),
            );
            let motion = div(force, dynamic_stiffness);
            let voltage = mul(voltage_per_motion, motion);
            let velocity = mul(C64::new(0.0, -omega), motion);
            let response = model
                .solve_harmonic(
                    omega,
                    &fixed,
                    &helectrodes(&cube, ZERO, HCondition::Current(ZERO)),
                    &loads,
                    &forces,
                    &vec![ZERO; n],
                    hc(),
                    cx,
                )
                .unwrap();
            assert!(response.maximum_force_residual_n <= hc().field.force_tolerance_n);
            assert!(response.maximum_charge_residual_c <= hc().field.charge_tolerance_c);
            cnear(response.displacement_m[node][2], motion, 1e-18, 3e-8);
            cnear(response.velocity_m_s[node][2], velocity, 1e-12, 3e-8);
            cnear(response.electrodes[1].voltage_v, voltage, 1e-9, 3e-8);
            let heat = 0.5 * voltage.norm_sq() / resistance;
            let supplied = 0.5 * mul(force, velocity.conj()).re;
            near(response.average_load_heat_w, heat, 1e-12, 3e-8);
            near(response.average_mechanical_power_w, supplied, 1e-12, 3e-8);
            near(supplied, heat, 1e-14, 3e-10);
            near(response.average_electrical_power_w, 0.0, 1e-12, 0.0);
            near(
                response.average_kinetic_energy_j,
                0.25 * mass * velocity.norm_sq(),
                1e-20,
                3e-8,
            );
            near(
                response.average_stored_energy_j,
                0.25 * (stiffness * motion.norm_sq() + capacitance * voltage.norm_sq()),
                1e-20,
                3e-8,
            );
        }
    });
}

#[test]
fn g3_complex_voltage_gauge_and_global_phase_preserve_physical_response() {
    with_cx(|cx, _| {
        let cube = cube(1);
        let n = cube.nodes.len();
        let materials = [material([0.0, 0.0, 12.0, 0.0, 0.0, 0.0])];
        let model =
            Model::assemble(&cube.nodes, &cube.tets, &materials, FRAME, budget(), cx).unwrap();
        let fixed: Vec<_> = (0..3 * n).map(|i| (i, ZERO)).collect();
        let loads = [HLoad::Resistive(super::transient_load::resistor(
            1e5,
            [1, 0],
        ))];
        let current = C64::new(2e-5, 3e-5);
        let phase = C64::new(0.6, 0.8);
        let shift = C64::new(73.0, -19.0);
        let solve = |reference, current| {
            model
                .solve_harmonic(
                    1e6,
                    &fixed,
                    &helectrodes(&cube, reference, HCondition::Current(current)),
                    &loads,
                    &vec![[ZERO; 3]; n],
                    &vec![ZERO; n],
                    hc(),
                    cx,
                )
                .unwrap()
        };
        let base = solve(ZERO, current);
        let transformed = solve(shift, mul(phase, current));
        for response in [&base, &transformed] {
            assert!(response.maximum_force_residual_n <= hc().field.force_tolerance_n);
            assert!(response.maximum_charge_residual_c <= hc().field.charge_tolerance_c);
        }
        for (a, b) in base.voltage_v.iter().zip(&transformed.voltage_v) {
            cnear(*b, add(shift, mul(phase, *a)), 1e-10, 2e-9);
        }
        for (a, b) in base.cells.iter().zip(&transformed.cells) {
            for axis in 0..3 {
                cnear(
                    b.electric_field_v_m[axis],
                    mul(phase, a.electric_field_v_m[axis]),
                    1e-7,
                    2e-9,
                );
            }
        }
        for (a, b) in base.electrodes.iter().zip(&transformed.electrodes) {
            cnear(b.charge_c, mul(phase, a.charge_c), 1e-20, 2e-9);
            cnear(
                b.source_current_a,
                mul(phase, a.source_current_a),
                1e-14,
                2e-9,
            );
        }
        near(
            transformed.average_load_heat_w,
            base.average_load_heat_w,
            1e-13,
            2e-9,
        );
        near(
            transformed.average_electrical_power_w,
            base.average_electrical_power_w,
            1e-13,
            2e-9,
        );
        near(
            transformed.average_stored_energy_j,
            base.average_stored_energy_j,
            1e-20,
            2e-9,
        );
    });
}

#[test]
fn g4_amplitude_limits_cover_the_whole_cycle_and_cancellation_refuses() {
    with_cx(|cx, gate| {
        let cube = cube(1);
        let n = cube.nodes.len();
        let materials = [material([0.0, 0.0, 12.0, 0.0, 0.0, 0.0])];
        let model =
            Model::assemble(&cube.nodes, &cube.tets, &materials, FRAME, budget(), cx).unwrap();
        let mut limits = hc();
        limits.field.maximum_strain_norm = 1e-5;
        let amplitude = 0.8 * limits.field.maximum_strain_norm;
        // Orthogonal real/imaginary strains trace a circle of radius amplitude;
        // their Frobenius sum is NOT the largest instantaneous strain.
        let circular: Vec<_> = cube
            .nodes
            .iter()
            .enumerate()
            .flat_map(|(node, p)| {
                [
                    (3 * node, C64::new(amplitude * p[0], 0.0)),
                    (3 * node + 1, ZERO),
                    (3 * node + 2, C64::new(0.0, amplitude * p[2])),
                ]
            })
            .collect();
        let linear: Vec<_> = cube
            .nodes
            .iter()
            .enumerate()
            .flat_map(|(node, p)| {
                [
                    (3 * node, ZERO),
                    (3 * node + 1, ZERO),
                    (3 * node + 2, C64::new(amplitude * p[2], amplitude * p[2])),
                ]
            })
            .collect();
        let solve = |fixed: &[(usize, C64)], voltage| {
            model.solve_harmonic(
                1e4,
                fixed,
                &helectrodes(&cube, ZERO, HCondition::Voltage(voltage)),
                &[],
                &vec![[ZERO; 3]; n],
                &vec![ZERO; n],
                limits,
                cx,
            )
        };
        assert!(solve(&circular, ZERO).is_ok());
        // Each quadrature separately lies within the envelope, but their
        // collinear sum reaches sqrt(2)*amplitude over the physical cycle.
        assert!(matches!(
            solve(&linear, ZERO),
            Err(TetPiezoelectricError::Invalid(_))
        ));
        let fixed: Vec<_> = (0..3 * n).map(|i| (i, ZERO)).collect();
        let voltage = 0.8 * limits.field.maximum_electric_field_v_m * LENGTH;
        assert!(matches!(
            solve(&fixed, C64::new(voltage, voltage)),
            Err(TetPiezoelectricError::Invalid(_))
        ));
        assert!(solve(&circular, ZERO).is_ok());
        gate.request();
        assert!(matches!(
            solve(&circular, ZERO),
            Err(TetPiezoelectricError::Cancelled)
        ));
    });
}
