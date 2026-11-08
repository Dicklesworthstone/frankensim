//! G1/G3/G4: opening a voltage source retains the actual charged material state.
//! The operation has no duration, impulse or work. Subsequent redistribution
//! follows the same native constitutive/circuit memories, including thermal bias.
use super::*;
use crate::harmonic::thermal_memory_source;
use fs_couple::piezoelectric::spatial::transient::load::ElectrodeDrive as Drive;
use fs_couple::piezoelectric::spatial::transient::load::reactive::{
    ReactiveElectrodeLoad as Reactive, ReactiveElectrodeState as Storage,
    ReactivePiezoelectricTransient as Shunted,
};
use fs_material::piezoelectric::tensor::relaxation::RelaxingTensorPiezoelectric;

const INCREMENT: f64 = 0.8 * EPSILON[2];
const TAU: f64 = 1e-6;
const AREA: f64 = LENGTH * LENGTH;
const VOLUME: f64 = AREA * LENGTH;
const EXPANSION: f64 = 2e-5;
const POLARIZATION: f64 = -3e-4;

fn thermal_model<'a>(cube: &'a Cube, cx: &Cx<'_>) -> (Model<'a>, RelaxingTensorPiezoelectric) {
    let (sample, law) = thermal_memory_source([0.0, f64::MAX]);
    let model = Model::assemble_thermal_materials(
        &cube.nodes,
        &cube.tets,
        vec![sample],
        FRAME,
        budget(),
        cx,
    )
    .unwrap();
    (model, law)
}

fn supports(values: &[[f64; 3]]) -> Vec<(usize, f64)> {
    values
        .iter()
        .enumerate()
        .flat_map(|(node, value)| (0..3).map(move |axis| (3 * node + axis, value[axis])))
        .collect()
}

fn unchanged(before: &Body<'_, '_>, after: &Body<'_, '_>) {
    let (a, b) = (before.state(), after.state());
    assert_eq!(a.time_s(), b.time_s());
    assert_eq!(a.displacement_m(), b.displacement_m());
    assert_eq!(a.velocity_m_s(), b.velocity_m_s());
    assert_eq!(a.voltage_v(), b.voltage_v());
    assert_eq!(a.boundary_charge_c(), b.boundary_charge_c());
    assert_eq!(a.material_memory_energy_j(), b.material_memory_energy_j());
    assert_eq!(a.total_energy_j(), b.total_energy_j());
    for (a, b) in a.cells().iter().zip(b.cells()) {
        assert_eq!(a.strain_mandel, b.strain_mandel);
        assert_eq!(a.response, b.response);
    }
    for (a, b) in before
        .material_relaxation_steps()
        .unwrap()
        .iter()
        .zip(after.material_relaxation_steps().unwrap())
    {
        assert_eq!(a.memory_energy_j, b.memory_energy_j);
        assert_eq!(a.memory_energy_change_j, b.memory_energy_change_j);
        assert_eq!(a.mechanical_heat_j, b.mechanical_heat_j);
        assert_eq!(a.dielectric_heat_j, b.dielectric_heat_j);
    }
}

#[test]
fn g1_released_hot_voltage_source_retains_partial_polarization_and_total_free_charge() {
    with_cx(|cx, _| {
        let cube = cube(1);
        let (model, law) = thermal_model(&cube, cx);
        let (zero, charge) = zero_source(cube.nodes.len());
        let c = EPSILON[2] * LENGTH;
        let cp = INCREMENT * LENGTH;
        let voltage = 20.0;
        let q_reference = -AREA * POLARIZATION;
        let body = Body::new(
            &model,
            &supports(&zero),
            &electrodes(&cube, 0.0, Condition::Voltage(voltage)),
            &zero,
            &zero,
            &charge,
            controls(),
            cx,
        )
        .unwrap()
        .with_relaxation(vec![law], 256, controls(), cx)
        .unwrap();
        let prepare_dt = 0.5 * TAU;
        let prepared = body
            .step_candidate(
                prepare_dt,
                &[Condition::Voltage(0.0), Condition::Voltage(voltage)],
                &zero,
                &charge,
                controls(),
                cx,
            )
            .unwrap()
            .into_parts()
            .0;
        // Voltage is held for a finite time, not until dielectric equilibrium.
        let mut polarization = cp * voltage * prepare_dt / (TAU + 0.5 * prepare_dt);
        assert!(polarization > 0.0 && polarization < cp * voltage);
        near(
            prepared.state().material_memory_energy_j(),
            0.5 * polarization * polarization / cp,
            2e-20,
            3e-8,
        );
        let free_charge = prepared.state().electrodes()[1].charge_c;
        near(
            free_charge,
            c * voltage + polarization + q_reference,
            2e-18,
            3e-8,
        );
        let mut released = prepared.release_electrodes_candidate(&[1], cx).unwrap();
        unchanged(&prepared, &released);
        let initial = released.state().total_energy_j();
        let elastic = 0.5 * VOLUME * STIFFNESS[2] * EXPANSION * EXPANSION;
        let mut v = voltage;
        let dt = 0.4 * TAU;
        let h = 0.5 * dt;
        let mut total_heat = 0.0;
        for _ in 0..12 {
            // C Vdot + pdot=0; tau pdot+p=Cp V. The retained p cannot
            // be reconstructed as zero or Cp V when the source opens.
            let vm = (c * v / h + polarization / (TAU + h)) / (c / h + cp / (TAU + h));
            let pm = (TAU * polarization + h * cp * vm) / (TAU + h);
            let next_v = 2.0 * vm - v;
            let next_p = 2.0 * pm - polarization;
            let heat = TAU * (next_p - polarization).powi(2) / (cp * dt);
            let memory = 0.5 * next_p * next_p / cp;
            let candidate = released
                .step_candidate(
                    dt,
                    &[Condition::Voltage(0.0), Condition::Charge(free_charge)],
                    &zero,
                    &charge,
                    controls(),
                    cx,
                )
                .unwrap();
            let state = candidate.body().state();
            near(state.electrodes()[1].voltage_v, next_v, 2e-9, 3e-8);
            near(state.electrodes()[1].charge_c, free_charge, 2e-18, 3e-8);
            near(state.material_memory_energy_j(), memory, 2e-20, 3e-8);
            near(
                state.stored_energy_j(),
                elastic + 0.5 * c * next_v * next_v + memory,
                2e-17,
                3e-8,
            );
            near(
                candidate.step().dielectric_material_heat_j,
                heat,
                2e-20,
                3e-8,
            );
            near(candidate.step().mechanical_material_heat_j, 0.0, 1e-24, 0.0);
            near(candidate.step().electrical_work_j, 0.0, 1e-20, 0.0);
            near(candidate.step().energy_residual_j, 0.0, 1e-19, 0.0);
            total_heat += heat;
            released = candidate.into_parts().0;
            (v, polarization) = (next_v, next_p);
        }
        assert!(total_heat > 0.0);
        near(
            released.state().total_energy_j() + total_heat,
            initial,
            2e-17,
            3e-8,
        );
        assert_eq!(prepared.state().time_s(), prepare_dt);
        assert_eq!(prepared.state().electrodes()[1].voltage_v, voltage);
    });
}

#[test]
fn g3_releasing_the_old_ground_preserves_fields_below_a_large_common_gauge() {
    with_cx(|cx, _| {
        let cube = cube(2);
        let (model, law) = thermal_model(&cube, cx);
        let (zero, mut charge) = zero_source(cube.nodes.len());
        let center = cube
            .nodes
            .iter()
            .position(|p| p.iter().all(|&x| x == 0.5 * LENGTH))
            .unwrap();
        charge[center] = 1.3e-13;
        let variables = law.internal_variables() * cube.tets.len();
        let initialize = |gauge| {
            let body = Body::new(
                &model,
                &supports(&zero),
                &electrodes(&cube, gauge, Condition::Voltage(gauge + 4.0)),
                &zero,
                &zero,
                &charge,
                controls(),
                cx,
            )
            .unwrap()
            .with_relaxation(vec![law.clone()], variables, controls(), cx)
            .unwrap();
            body.step_candidate(
                0.3 * TAU,
                &[Condition::Voltage(gauge), Condition::Voltage(gauge + 4.0)],
                &zero,
                &charge,
                controls(),
                cx,
            )
            .unwrap()
            .into_parts()
            .0
        };
        let plain = initialize(0.0);
        let gauge = 2.0_f64.powi(40);
        let shifted = initialize(gauge);
        // Interior volume charge makes nodal differences non-binary. Rebase
        // the original relative potentials; subtracting rounded absolute
        // nodal voltages would silently corrupt the next electric field.
        let plain_charge = plain.state().electrodes()[0].charge_c;
        let shifted_charge = shifted.state().electrodes()[0].charge_c;
        let mut a = plain.release_electrodes_candidate(&[0], cx).unwrap();
        let mut b = shifted.release_electrodes_candidate(&[0], cx).unwrap();
        unchanged(&plain, &a);
        unchanged(&shifted, &b);
        for _ in 0..3 {
            let next_a = a
                .step_candidate(
                    0.2 * TAU,
                    &[Condition::Charge(plain_charge), Condition::Voltage(4.0)],
                    &zero,
                    &charge,
                    controls(),
                    cx,
                )
                .unwrap();
            let next_b = b
                .step_candidate(
                    0.2 * TAU,
                    &[
                        Condition::Charge(shifted_charge),
                        Condition::Voltage(gauge + 4.0),
                    ],
                    &zero,
                    &charge,
                    controls(),
                    cx,
                )
                .unwrap();
            let (sa, sb) = (next_a.body().state(), next_b.body().state());
            assert_eq!(sa.displacement_m(), sb.displacement_m());
            assert_eq!(sa.velocity_m_s(), sb.velocity_m_s());
            near(
                sa.material_memory_energy_j(),
                sb.material_memory_energy_j(),
                2e-20,
                3e-8,
            );
            near(sa.total_energy_j(), sb.total_energy_j(), 2e-17, 3e-8);
            near(
                next_a.step().material_heat_j,
                next_b.step().material_heat_j,
                2e-20,
                3e-8,
            );
            for (a, b) in sa.cells().iter().zip(sb.cells()) {
                for axis in 0..3 {
                    near(
                        a.response.electric_field_v_m[axis],
                        b.response.electric_field_v_m[axis],
                        2e-5,
                        3e-8,
                    );
                    near(
                        a.response.displacement_c_m2[axis],
                        b.response.displacement_c_m2[axis],
                        2e-12,
                        3e-8,
                    );
                }
                for axis in 0..6 {
                    near(
                        a.response.stress_mandel_pa[axis],
                        b.response.stress_mandel_pa[axis],
                        0.05,
                        3e-8,
                    );
                }
            }
            near(sa.electrodes()[0].charge_c, plain_charge, 2e-18, 3e-8);
            near(sb.electrodes()[0].charge_c, shifted_charge, 2e-18, 3e-8);
            near(next_a.step().electrical_work_j, 0.0, 1e-20, 0.0);
            near(next_b.step().electrical_work_j, 0.0, 1e-20, 0.0);
            a = next_a.into_parts().0;
            b = next_b.into_parts().0;
        }
    });
}

#[test]
fn g3_g4_release_preserves_moving_supports_and_refuses_invalid_or_cancelled_changes() {
    with_cx(|cx, gate| {
        let cube = cube(1);
        let (model, law) = thermal_model(&cube, cx);
        let n = cube.nodes.len();
        let (zero, charge) = zero_source(n);
        let body = Body::new(
            &model,
            &supports(&zero),
            &electrodes(&cube, 0.0, Condition::Voltage(4.0)),
            &zero,
            &zero,
            &charge,
            controls(),
            cx,
        )
        .unwrap()
        .with_relaxation(vec![law], 256, controls(), cx)
        .unwrap();
        let dt = 1e-7;
        let velocity = vec![[0.01, -0.02, 0.03]; n];
        let driven = supports(&velocity);
        let moving = body
            .step_with_support_motion_candidate(
                dt,
                &driven,
                &[Condition::Voltage(0.0), Condition::Voltage(4.0)],
                &zero,
                &charge,
                controls(),
                cx,
            )
            .unwrap()
            .into_parts()
            .0;
        for indices in [&[][..], &[2][..], &[1, 1][..], &[0, 1][..]] {
            assert!(moving.release_electrodes_candidate(indices, cx).is_err());
        }
        let released = moving.release_electrodes_candidate(&[1], cx).unwrap();
        unchanged(&moving, &released);
        assert!(released.release_electrodes_candidate(&[1], cx).is_err());
        assert!(released.release_electrodes_candidate(&[0], cx).is_err());
        let q = released.state().electrodes()[1].charge_c;
        let candidate = released
            .step_with_support_motion_candidate(
                dt,
                &driven,
                &[Condition::Voltage(0.0), Condition::Charge(q)],
                &zero,
                &charge,
                controls(),
                cx,
            )
            .unwrap();
        for node in 0..n {
            for axis in 0..3 {
                near(
                    candidate.body().state().displacement_m()[node][axis],
                    moving.state().displacement_m()[node][axis] + dt * velocity[node][axis],
                    2e-18,
                    3e-8,
                );
                near(
                    candidate.body().state().velocity_m_s()[node][axis],
                    velocity[node][axis],
                    2e-12,
                    3e-8,
                );
            }
        }
        unchanged(&moving, &released);
        gate.request();
        assert!(matches!(
            moving.release_electrodes_candidate(&[1], cx),
            Err(TetPiezoelectricError::Cancelled)
        ));
        assert_eq!(moving.state().time_s(), dt);
        assert_eq!(moving.state().velocity_m_s(), &velocity);
    });
}

#[test]
fn g1_reactive_release_preserves_external_flux_charge_and_total_energy() {
    use crate::transient_load::reactive::{capacitor, inductor};
    with_cx(|cx, _| {
        let cube = cube(1);
        let (model, law) = thermal_model(&cube, cx);
        let (zero, charge) = zero_source(cube.nodes.len());
        let voltage = 4.0;
        let body = Body::new_with_equilibrium_relaxation(
            &model,
            &supports(&zero),
            &electrodes(&cube, 0.0, Condition::Voltage(voltage)),
            &zero,
            &zero,
            &charge,
            vec![law],
            256,
            controls(),
            cx,
        )
        .unwrap();
        let capacitance = 2.0 * EPSILON[2] * LENGTH;
        let inductance = 2e-3;
        let flux = inductance * 1e-3;
        let capacitor_charge = capacitance * voltage;
        let external = 0.5 * flux * flux / inductance
            + 0.5 * capacitor_charge * capacitor_charge / capacitance;
        let initial = body.state().total_energy_j() + external;
        let coupled = Shunted::new(
            body,
            vec![
                Reactive::SeriesRl {
                    electrodes: [1, 0],
                    inductor: inductor(inductance),
                    conductor: None,
                },
                Reactive::Capacitor {
                    electrodes: [1, 0],
                    capacitor: capacitor(capacitance),
                },
            ],
            vec![
                Storage::FluxLinkage(flux),
                Storage::Charge(capacitor_charge),
            ],
            vec![],
            controls(),
            cx,
        )
        .unwrap();
        let released = coupled.release_electrodes_candidate(&[1], cx).unwrap();
        unchanged(coupled.body(), released.body());
        assert_eq!(coupled.states(), released.states());
        assert_eq!(coupled.loads().len(), released.loads().len());
        let candidate = released
            .step_candidate(
                1e-8,
                &[Drive::Voltage(0.0), Drive::Current(0.0)],
                &zero,
                &charge,
                controls(),
                cx,
            )
            .unwrap();
        let next = candidate.model();
        let next_flux = match next.states()[0] {
            Storage::FluxLinkage(value) => value,
            _ => panic!("original winding changed storage kind"),
        };
        let next_charge = match next.states()[1] {
            Storage::Charge(value) => value,
            _ => panic!("original capacitor changed storage kind"),
        };
        let next_external = 0.5 * next_flux * next_flux / inductance
            + 0.5 * next_charge * next_charge / capacitance;
        assert_ne!(next_flux, flux);
        assert_ne!(next_charge, capacitor_charge);
        near(candidate.step().reactive_energy_j[0], external, 2e-20, 3e-8);
        near(
            candidate.step().reactive_energy_j[1],
            next_external,
            2e-20,
            3e-8,
        );
        near(
            next.body().state().total_energy_j()
                + next_external
                + candidate.step().body.material_heat_j,
            initial,
            2e-17,
            3e-8,
        );
        near(candidate.step().electrical_source_work_j, 0.0, 1e-20, 0.0);
        assert_eq!(candidate.step().joule_heat_j, 0.0);
        near(candidate.step().energy_residual_j, 0.0, 1e-19, 0.0);
        assert_eq!(
            coupled.states(),
            &[
                Storage::FluxLinkage(flux),
                Storage::Charge(capacitor_charge)
            ]
        );
        assert_eq!(coupled.body().state().time_s(), 0.0);
    });
}
