//! Physical midpoint oracles for retained Maxwell and polarization memories.
//! These reuse the explicit synthetic C0/epsilon_infinity source roles; no
//! experimental loss, thermal receiver or broad-band material fit is inferred.
use super::*;
use crate::harmonic::relaxation::fixture;
use crate::transient_load::{reactive::inductor, resistor};
use fs_couple::piezoelectric::spatial::transient::load::{
    ElectrodeDrive as Drive,
    reactive::{
        ReactiveElectrodeLoad as Reactive, ReactiveElectrodeState as State,
        ReactivePiezoelectricTransient as Shunted,
    },
};
use fs_material::piezoelectric::tensor::relaxation::RelaxingTensorPiezoelectric;

const MEMORY_BUDGET: usize = 256;

struct OneCoordinate {
    node: usize,
    fixed: Vec<(usize, f64)>,
    stiffness: f64,
    mass: f64,
    coupling: f64,
}

fn one_coordinate(model: &Model<'_>, cube: &Cube) -> OneCoordinate {
    let node = *cube.top.last().unwrap();
    let dof = 3 * node + 2;
    let diagonal = |matrix: &fs_sparse::Csr| {
        let (cols, values) = matrix.row(dof);
        values[cols.iter().position(|&i| i == dof).unwrap()]
    };
    let (cols, values) = model.coupling().row(dof);
    OneCoordinate {
        node,
        fixed: (0..3 * cube.nodes.len())
            .filter(|&i| i != dof)
            .map(|i| (i, 0.0))
            .collect(),
        stiffness: diagonal(&model.solid().stiffness),
        mass: diagonal(&model.solid().mass),
        coupling: cols
            .iter()
            .zip(values)
            .filter(|(col, _)| cube.top.contains(col))
            .map(|(_, value)| value)
            .sum(),
    }
}

// Preserve the shared fixture's actual coefficient declarations, changing only
// one auxiliary source's applicability. A frequency-only fit is admissible in
// harmonic response but cannot certify an arbitrary transient history.
fn finite_band_fixture() -> (TensorPiezoelectric, RelaxingTensorPiezoelectric) {
    let (reference, memory) = fixture(None, Some((0.8 * EPSILON[2], 1e-6)), FRAME);
    let mut claims = ClaimSet::new();
    for property in reference
        .material()
        .properties()
        .iter()
        .chain(memory.material().properties())
    {
        let name = property.requirement().name();
        let dims = property.requirement().dims();
        let mut validity = ValidityDomain::unconstrained().with("T", 290.0, 310.0);
        if name == "memory_tau_p" {
            validity = validity.with("omega", 10.0, 1000.0);
        }
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value: PropertyValue::Scalar {
                    value: property.value_si(),
                    dims,
                },
                validity,
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::ConstantWithinValidity,
                observations: vec![],
                provenance: Provenance {
                    source: "authored finite-band polarization fixture".into(),
                    license: "CC0-1.0".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic reciprocal ceramic with bounded loss evidence".into(),
            phase: "solid".into(),
            process: "explicit finite-band applicability".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    let base = TensorPiezoelectric::resolve(
        &card,
        &QueryPoint::new()
            .with("T", 300.0)
            .unwrap()
            .with("omega", 100.0)
            .unwrap(),
        MaterialPropertySelection::SingleClaimOnly,
        FRAME,
    )
    .unwrap();
    let law = RelaxingTensorPiezoelectric::resolve(
        base.clone(),
        &card,
        vec![],
        memory.dielectric_properties().to_vec(),
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    (base, law)
}

#[test]
fn g1_clamped_polarization_retains_current_driven_charge_storage_and_heat() {
    with_cx(|cx, _| {
        let cube = cube(1);
        let n = cube.nodes.len();
        let (zero, charge) = zero_source(n);
        let fixed: Vec<_> = (0..3 * n).map(|i| (i, 0.0)).collect();
        let increment = 0.8 * EPSILON[2];
        let tau = 1e-6;
        let (base, loss) = fixture(None, Some((increment, tau)), FRAME);
        let materials = [base];
        let model =
            Model::assemble(&cube.nodes, &cube.tets, &materials, FRAME, budget(), cx).unwrap();
        let mut body = Body::new(
            &model,
            &fixed,
            &electrodes(&cube, 0.0, Condition::Charge(0.0)),
            &zero,
            &zero,
            &charge,
            controls(),
            cx,
        )
        .unwrap()
        .with_relaxation(vec![loss], MEMORY_BUDGET, controls(), cx)
        .unwrap();
        let c = EPSILON[2] * LENGTH;
        let cp = increment * LENGTH;
        let dt = 0.4 * tau;
        let h = 0.5 * dt;
        let (mut voltage, mut polarization, mut free_charge) = (0.0_f64, 0.0_f64, 0.0_f64);
        let (mut supplied, mut heat_sum) = (0.0, 0.0);
        assert_eq!(body.state().material_memory_energy_j(), 0.0);
        for current in [4e-5, 4e-5, 4e-5, 0.0, 0.0, -2e-5, 0.0] {
            // Positive terminal polarization charge P obeys tau Pdot+P=Cp V.
            // I=C Vdot+Pdot fixes V_mid without using a material candidate.
            let vm =
                (c * voltage / h + current + polarization / (tau + h)) / (c / h + cp / (tau + h));
            let pm = (tau * polarization + h * cp * vm) / (tau + h);
            let next_v = 2.0 * vm - voltage;
            let next_p = 2.0 * pm - polarization;
            let memory = 0.5 * next_p * next_p / cp;
            let heat = tau * (next_p - polarization).powi(2) / (cp * dt);
            let energy = 0.5 * c * next_v * next_v + memory;
            let previous_time = body.state().time_s();
            let previous_memory = body.state().material_memory_energy_j();
            let (candidate, loaded) = body
                .step_with_loads_candidate(
                    dt,
                    &[Drive::Voltage(0.0), Drive::Current(current)],
                    &[],
                    &zero,
                    &charge,
                    controls(),
                    cx,
                )
                .unwrap();
            let step = candidate.step();
            let state = candidate.body().state();
            free_charge += dt * current;
            near(state.electrodes()[1].voltage_v, next_v, 1e-10, 3e-8);
            near(state.electrodes()[1].charge_c, free_charge, 1e-20, 3e-8);
            near(state.stored_energy_j(), energy, 1e-20, 3e-8);
            near(state.material_memory_energy_j(), memory, 1e-20, 3e-8);
            near(
                step.material_memory_energy_j[0],
                previous_memory,
                1e-20,
                3e-8,
            );
            near(step.material_memory_energy_j[1], memory, 1e-20, 3e-8);
            near(step.material_heat_j, heat, 1e-20, 3e-8);
            near(step.dielectric_material_heat_j, heat, 1e-20, 3e-8);
            assert_eq!(step.mechanical_material_heat_j, 0.0);
            assert_eq!(loaded.joule_heat_j, 0.0);
            near(
                loaded.electrical_source_work_j,
                dt * current * vm,
                1e-20,
                3e-8,
            );
            near(step.energy_residual_j, 0.0, 1e-20, 0.0);
            near(loaded.energy_residual_j, 0.0, 1e-20, 0.0);
            let cell_steps = candidate.body().material_relaxation_steps().unwrap();
            near(
                cell_steps.iter().map(|s| s.dielectric_heat_j).sum(),
                heat,
                1e-20,
                3e-8,
            );
            near(
                cell_steps.iter().map(|s| s.memory_energy_j[1]).sum(),
                memory,
                1e-20,
                3e-8,
            );
            for cell in state.cells() {
                near(
                    cell.response.electric_field_v_m[2],
                    -next_v / LENGTH,
                    1e-7,
                    3e-8,
                );
                near(
                    cell.response.displacement_c_m2[2],
                    -free_charge / (LENGTH * LENGTH),
                    1e-17,
                    3e-8,
                );
            }
            // A candidate includes fresh material histories; the accepted one
            // remains byte-for-byte at its previous physical time and storage.
            assert_eq!(body.state().time_s(), previous_time);
            assert_eq!(body.state().material_memory_energy_j(), previous_memory);
            supplied += dt * current * vm;
            heat_sum += heat;
            body = candidate.into_parts().0;
            (voltage, polarization) = (next_v, next_p);
        }
        near(
            body.state().total_energy_j() + heat_sum,
            supplied,
            1e-19,
            3e-8,
        );
    });
}

#[test]
fn g1_open_circuit_maxwell_motion_matches_an_independent_retained_state_midpoint_oracle() {
    with_cx(|cx, _| {
        let cube = cube(1);
        let n = cube.nodes.len();
        let (zero, charge) = zero_source(n);
        let fraction = 0.4;
        let tau = 2e-7;
        let (base, loss) = fixture(Some((fraction, tau)), None, FRAME);
        let materials = [base];
        let model =
            Model::assemble(&cube.nodes, &cube.tets, &materials, FRAME, budget(), cx).unwrap();
        let coordinate = one_coordinate(&model, &cube);
        let capacitance = EPSILON[2] * LENGTH;
        let k = coordinate.stiffness + coordinate.coupling.powi(2) / capacitance;
        let km = fraction * coordinate.stiffness;
        let amplitude = 1e-8;
        let mut initial_u = zero.clone();
        initial_u[coordinate.node][2] = amplitude;
        let plain = Body::new(
            &model,
            &coordinate.fixed,
            &electrodes(&cube, 0.0, Condition::Charge(0.0)),
            &initial_u,
            &zero,
            &charge,
            controls(),
            cx,
        )
        .unwrap();
        let shape: Vec<_> = plain
            .state()
            .cells()
            .iter()
            .map(|cell| cell.strain_mandel.map(|x| x / amplitude))
            .collect();
        let initial_energy = plain.state().total_energy_j();
        let initial_voltages = plain.state().voltage_v().to_vec();
        let mut body = plain
            .with_relaxation(vec![loss], MEMORY_BUDGET, controls(), cx)
            .unwrap();
        assert_eq!(body.state().material_memory_energy_j(), 0.0);
        assert_eq!(body.state().total_energy_j(), initial_energy);
        assert_eq!(body.state().voltage_v(), initial_voltages);
        let dt = 0.2 * (coordinate.mass / k).sqrt();
        let h = 0.5 * dt;
        let a = tau / (tau + h);
        let (mut u, mut v, mut memory_force) = (amplitude, 0.0, 0.0);
        let mut heat_sum = 0.0;
        for _ in 0..32 {
            // M vdot + k_open u + w = 0; wdot=km udot-w/tau.
            // Open-circuit voltage follows V=G u/C at every endpoint.
            let vm = (coordinate.mass * v - h * k * u - h * a * memory_force)
                / (coordinate.mass + h * h * (k + a * km));
            let wm = a * (memory_force + h * km * vm);
            let next_u = u + dt * vm;
            let next_v = 2.0 * vm - v;
            let next_w = 2.0 * wm - memory_force;
            let memory = 0.5 * next_w * next_w / km;
            let heat = dt * wm * wm / (km * tau);
            let candidate = body
                .step_candidate(
                    dt,
                    &[Condition::Voltage(0.0), Condition::Charge(0.0)],
                    &zero,
                    &charge,
                    controls(),
                    cx,
                )
                .unwrap();
            let state = candidate.body().state();
            let step = candidate.step();
            near(
                state.displacement_m()[coordinate.node][2],
                next_u,
                1e-18,
                3e-8,
            );
            near(
                state.velocity_m_s()[coordinate.node][2],
                next_v,
                1e-11,
                3e-8,
            );
            near(
                state.electrodes()[1].voltage_v,
                coordinate.coupling * next_u / capacitance,
                1e-10,
                3e-8,
            );
            near(state.electrodes()[1].charge_c, 0.0, 1e-18, 0.0);
            near(state.material_memory_energy_j(), memory, 1e-20, 3e-8);
            near(
                state.stored_energy_j(),
                0.5 * k * next_u * next_u + memory,
                1e-20,
                3e-8,
            );
            near(
                state.kinetic_energy_j(),
                0.5 * coordinate.mass * next_v * next_v,
                1e-20,
                3e-8,
            );
            near(step.mechanical_material_heat_j, heat, 1e-20, 3e-8);
            near(step.material_heat_j, heat, 1e-20, 3e-8);
            assert_eq!(step.dielectric_material_heat_j, 0.0);
            assert_eq!(step.mechanical_work_j, 0.0);
            near(step.electrical_work_j, 0.0, 1e-20, 0.0);
            near(step.energy_residual_j, 0.0, 1e-20, 0.0);
            for (cell, shape) in state.cells().iter().zip(&shape) {
                for axis in 0..6 {
                    let electric = if axis == 2 {
                        12.0 * cell.response.electric_field_v_m[2]
                    } else {
                        0.0
                    };
                    let stress = STIFFNESS[axis] * cell.strain_mandel[axis] - electric
                        + fraction * STIFFNESS[axis] * shape[axis] * next_w / km;
                    near(cell.response.stress_mandel_pa[axis], stress, 1e-7, 3e-8);
                }
            }
            heat_sum += heat;
            body = candidate.into_parts().0;
            (u, v, memory_force) = (next_u, next_v, next_w);
        }
        assert!(heat_sum > 0.0);
        near(
            body.state().total_energy_j() + heat_sum,
            initial_energy,
            1e-19,
            3e-8,
        );
    });
}

#[test]
fn g1_intrinsic_and_rl_shunt_memories_share_one_trajectory_and_all_heat_accounts() {
    with_cx(|cx, _| {
        let cube = cube(1);
        let n = cube.nodes.len();
        let (zero, charge) = zero_source(n);
        let (fraction, tau_m) = (0.3, 2e-7);
        let (increment, tau_p) = (0.6 * EPSILON[2], 4e-7);
        let (base, loss) = fixture(Some((fraction, tau_m)), Some((increment, tau_p)), FRAME);
        let materials = [base];
        let model =
            Model::assemble(&cube.nodes, &cube.tets, &materials, FRAME, budget(), cx).unwrap();
        let coordinate = one_coordinate(&model, &cube);
        let (c, cp) = (EPSILON[2] * LENGTH, increment * LENGTH);
        let km = fraction * coordinate.stiffness;
        let (l, rp, rs) = (2e-3, 2e4, 3e3);
        let (mut u, mut v, mut voltage) = (2e-9, 0.0, 5.0);
        let (mut w, mut p, mut flux) = (0.0_f64, 0.0_f64, l * 1e-4);
        let mut initial_u = zero.clone();
        initial_u[coordinate.node][2] = u;
        let body = Body::new(
            &model,
            &coordinate.fixed,
            &electrodes(
                &cube,
                0.0,
                Condition::Charge(c * voltage - coordinate.coupling * u),
            ),
            &initial_u,
            &zero,
            &charge,
            controls(),
            cx,
        )
        .unwrap()
        .with_relaxation(vec![loss], MEMORY_BUDGET, controls(), cx)
        .unwrap();
        let initial_energy = body.state().total_energy_j() + 0.5 * flux * flux / l;
        near(body.state().electrodes()[1].voltage_v, voltage, 1e-10, 3e-8);
        assert_eq!(body.state().material_memory_energy_j(), 0.0);
        let mut coupled = Shunted::new(
            body,
            vec![Reactive::SeriesRl {
                electrodes: [1, 0],
                inductor: inductor(l),
                conductor: Some(resistor(rs, [1, 0]).conductor),
            }],
            vec![State::FluxLinkage(flux)],
            vec![resistor(rp, [1, 0])],
            controls(),
            cx,
        )
        .unwrap();
        let dt = 3e-8;
        let h = 0.5 * dt;
        let am = tau_m / (tau_m + h);
        let (mut heat_sum, mut supplied) = (0.0, 0.0);
        for current in [
            5e-5, 5e-5, 5e-5, 0.0, 0.0, 0.0, -2e-5, -2e-5, -2e-5, 0.0, 0.0, 0.0,
        ] {
            // Independently eliminate w_mid, P_mid, and I_L_mid. The two
            // remaining equations are A u_mid+G V_mid=f and
            // -G u_mid/h+B V_mid=j, retaining all five physical state variables.
            let a = coordinate.mass / (h * h) + coordinate.stiffness + am * km;
            let b = c / h + cp / (tau_p + h) + 1.0 / rp + h / (l + h * rs);
            let f = coordinate.mass * (u + h * v) / (h * h) + am * km * u - am * w;
            let j = current + c * voltage / h - coordinate.coupling * u / h + p / (tau_p + h)
                - flux / (l + h * rs);
            let determinant = a * b + coordinate.coupling.powi(2) / h;
            let um = (f * b - coordinate.coupling * j) / determinant;
            let vm = (a * j + coordinate.coupling * f / h) / determinant;
            let wm = am * (w + km * (um - u));
            let pm = (tau_p * p + h * cp * vm) / (tau_p + h);
            let im = (flux + h * vm) / (l + h * rs);
            let (next_u, next_v, next_voltage) =
                (2.0 * um - u, 2.0 * (um - u) / h - v, 2.0 * vm - voltage);
            let (next_w, next_p, next_flux) = (2.0 * wm - w, 2.0 * pm - p, 2.0 * l * im - flux);
            let memory_before = 0.5 * w * w / km + 0.5 * p * p / cp;
            let memory_after = 0.5 * next_w * next_w / km + 0.5 * next_p * next_p / cp;
            let mechanical_heat = dt * wm * wm / (km * tau_m);
            let dielectric_heat = tau_p * (next_p - p).powi(2) / (cp * dt);
            let external_heat = dt * (vm * vm / rp + rs * im * im);
            let body_energy = 0.5 * coordinate.mass * next_v * next_v
                + 0.5 * coordinate.stiffness * next_u * next_u
                + 0.5 * c * next_voltage * next_voltage
                + memory_after;
            let previous_time = coupled.body().state().time_s();
            let previous_flux = coupled.states()[0];
            let candidate = coupled
                .step_candidate(
                    dt,
                    &[Drive::Voltage(0.0), Drive::Current(current)],
                    &zero,
                    &charge,
                    controls(),
                    cx,
                )
                .unwrap();
            let state = candidate.model().body().state();
            let step = candidate.step();
            near(
                state.displacement_m()[coordinate.node][2],
                next_u,
                1e-18,
                3e-8,
            );
            near(
                state.velocity_m_s()[coordinate.node][2],
                next_v,
                1e-11,
                3e-8,
            );
            near(state.electrodes()[1].voltage_v, next_voltage, 1e-10, 3e-8);
            near(
                state.electrodes()[1].charge_c,
                c * next_voltage - coordinate.coupling * next_u + next_p,
                1e-20,
                3e-8,
            );
            let State::FluxLinkage(actual_flux) = candidate.model().states()[0] else {
                panic!("winding state changed kind")
            };
            near(actual_flux, next_flux, 1e-18, 3e-8);
            near(step.reactive[0].midpoint_current_a, im, 1e-13, 3e-8);
            near(step.resistors[0].midpoint_current_a, vm / rp, 1e-13, 3e-8);
            near(state.material_memory_energy_j(), memory_after, 1e-20, 3e-8);
            near(state.total_energy_j(), body_energy, 1e-20, 3e-8);
            near(
                step.body.material_memory_energy_j[0],
                memory_before,
                1e-20,
                3e-8,
            );
            near(
                step.body.material_memory_energy_j[1],
                memory_after,
                1e-20,
                3e-8,
            );
            near(
                step.body.mechanical_material_heat_j,
                mechanical_heat,
                1e-20,
                3e-8,
            );
            near(
                step.body.dielectric_material_heat_j,
                dielectric_heat,
                1e-20,
                3e-8,
            );
            near(
                step.body.material_heat_j,
                mechanical_heat + dielectric_heat,
                1e-20,
                3e-8,
            );
            near(step.joule_heat_j, external_heat, 1e-20, 3e-8);
            near(
                step.reactive_energy_j[1],
                0.5 * next_flux * next_flux / l,
                1e-20,
                3e-8,
            );
            near(
                step.electrical_source_work_j,
                dt * current * vm,
                1e-20,
                3e-8,
            );
            near(
                step.body.electrical_work_j,
                vm * (c * (next_voltage - voltage) - coordinate.coupling * (next_u - u) + next_p
                    - p),
                1e-20,
                3e-8,
            );
            near(step.body.energy_residual_j, 0.0, 1e-20, 0.0);
            near(step.energy_residual_j, 0.0, 1e-20, 0.0);
            assert_eq!(coupled.body().state().time_s(), previous_time);
            near(
                coupled.body().state().material_memory_energy_j(),
                memory_before,
                1e-20,
                3e-8,
            );
            assert_eq!(coupled.states()[0], previous_flux);
            heat_sum += mechanical_heat + dielectric_heat + external_heat;
            supplied += dt * current * vm;
            coupled = candidate.into_parts().0;
            (u, v, voltage, w, p, flux) = (next_u, next_v, next_voltage, next_w, next_p, next_flux);
        }
        near(
            coupled.body().state().total_energy_j() + 0.5 * flux * flux / l + heat_sum,
            initial_energy + supplied,
            1e-19,
            3e-8,
        );
    });
}

#[test]
fn g4_transient_memory_refuses_bounded_evidence_and_preserves_history_on_failure() {
    with_cx(|cx, gate| {
        let cube = cube(1);
        let n = cube.nodes.len();
        let (zero, charge) = zero_source(n);
        let fixed: Vec<_> = (0..3 * n).map(|i| (i, 0.0)).collect();
        let pole = Some((0.8 * EPSILON[2], 1e-6));
        let (base, loss) = fixture(None, pole, FRAME);
        let (_, foreign_source) = fixture(None, Some((0.7 * EPSILON[2], 1e-6)), FRAME);
        let (_, foreign_frame) = fixture(None, pole, ContentHash([81; 32]));
        let materials = [base];
        let model =
            Model::assemble(&cube.nodes, &cube.tets, &materials, FRAME, budget(), cx).unwrap();
        let plain = Body::new(
            &model,
            &fixed,
            &electrodes(&cube, 0.0, Condition::Charge(0.0)),
            &zero,
            &zero,
            &charge,
            controls(),
            cx,
        )
        .unwrap();
        assert!(
            plain
                .clone()
                .with_relaxation(vec![loss.clone()], 0, controls(), cx)
                .is_err()
        );
        assert!(
            plain
                .clone()
                .with_relaxation(vec![foreign_source], MEMORY_BUDGET, controls(), cx)
                .is_err()
        );
        assert!(
            plain
                .clone()
                .with_relaxation(vec![foreign_frame], MEMORY_BUDGET, controls(), cx)
                .is_err()
        );
        let (bounded_base, bounded_loss) = finite_band_fixture();
        assert_eq!(bounded_loss.frequency_band_rad_s(), [10.0, 1000.0]);
        assert!(bounded_loss.frequency_coefficients(100.0).is_ok());
        let bounded_materials = [bounded_base];
        let bounded_model = Model::assemble(
            &cube.nodes,
            &cube.tets,
            &bounded_materials,
            FRAME,
            budget(),
            cx,
        )
        .unwrap();
        let bounded_body = Body::new(
            &bounded_model,
            &fixed,
            &electrodes(&cube, 0.0, Condition::Charge(0.0)),
            &zero,
            &zero,
            &charge,
            controls(),
            cx,
        )
        .unwrap();
        assert!(
            bounded_body
                .with_relaxation(vec![bounded_loss], MEMORY_BUDGET, controls(), cx)
                .is_err()
        );
        assert_eq!(plain.state().time_s(), 0.0);
        assert_eq!(plain.state().material_memory_energy_j(), 0.0);
        let body = plain
            .with_relaxation(vec![loss], MEMORY_BUDGET, controls(), cx)
            .unwrap()
            .step_candidate(
                2e-7,
                &[Condition::Voltage(0.0), Condition::Charge(2e-11)],
                &zero,
                &charge,
                controls(),
                cx,
            )
            .unwrap()
            .into_parts()
            .0;
        let before_time = body.state().time_s();
        let before_voltage = body.state().voltage_v().to_vec();
        let before_memory = body.state().material_memory_energy_j();
        let before_energy = body.state().total_energy_j();
        let before_d: Vec<_> = body
            .state()
            .cells()
            .iter()
            .map(|c| c.response.displacement_c_m2)
            .collect();
        assert!(before_memory > 0.0);
        assert!(
            body.step_candidate(
                2e-7,
                &[Condition::Voltage(0.0), Condition::Charge(1e-5)],
                &zero,
                &charge,
                controls(),
                cx
            )
            .is_err()
        );
        assert_eq!(body.state().time_s(), before_time);
        assert_eq!(body.state().voltage_v(), before_voltage);
        assert_eq!(body.state().material_memory_energy_j(), before_memory);
        assert_eq!(body.state().total_energy_j(), before_energy);
        assert_eq!(
            body.state()
                .cells()
                .iter()
                .map(|c| c.response.displacement_c_m2)
                .collect::<Vec<_>>(),
            before_d
        );
        // A valid candidate still starts from the same accepted polarization.
        assert!(
            body.step_candidate(
                2e-7,
                &[Condition::Voltage(0.0), Condition::Charge(2e-11)],
                &zero,
                &charge,
                controls(),
                cx
            )
            .is_ok()
        );
        gate.request();
        assert!(matches!(
            body.step_candidate(
                2e-7,
                &[Condition::Voltage(0.0), Condition::Charge(2e-11)],
                &zero,
                &charge,
                controls(),
                cx
            ),
            Err(TetPiezoelectricError::Cancelled)
        ));
        assert_eq!(body.state().time_s(), before_time);
        assert_eq!(body.state().material_memory_energy_j(), before_memory);
        assert_eq!(body.state().total_energy_j(), before_energy);
    });
}
