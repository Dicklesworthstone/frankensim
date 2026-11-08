//! G1/G3 fixed-temperature thermoelastic and pyroelectric spatial responses.
//! Synthetic constant SI coefficients give independent affine-field oracles;
//! no measured ceramic or inferred pyroelectric convention is claimed.
use super::*;
use fs_couple::piezoelectric::spatial::transient::TetPiezoelectricTransient as Body;
use fs_material::piezoelectric::{
    tensor::relaxation::RelaxingTensorPiezoelectric, thermal::ThermalTensorPiezoelectric,
};
use fs_qty::{
    QuantitySpec,
    semantic::{QuantityKind, SemanticType, ValueForm},
};

const COUPLING: f64 = 12.0;
const ALPHA: f64 = 2e-6;
const PYRO: f64 = -3e-5;
const REFERENCE: f64 = 300.0;
const HOT: f64 = 310.0;

fn temperature() -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(
        QuantityKind::AbsoluteTemperature,
        ValueForm::Static,
    ))
}
fn point(t: f64) -> QueryPoint {
    QueryPoint::new()
        .with_quantity("T", temperature(), t)
        .unwrap()
}
pub(super) fn source() -> (MaterialCard, ThermalTensorPiezoelectric) {
    let tensor = material([0.0, 0.0, COUPLING, 0.0, 0.0, 0.0]);
    let mut claims = ClaimSet::new();
    let mut add = |name: &str, dims, value| {
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value: PropertyValue::Scalar { value, dims },
                validity: ValidityDomain::unconstrained().with_quantity(
                    "T",
                    temperature(),
                    280.0,
                    340.0,
                ),
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::ConstantWithinValidity,
                observations: vec![],
                provenance: Provenance {
                    source: "authored clamped-pyroelectric spatial fixture".into(),
                    license: "CC0-1.0".into(),
                    artifact: None,
                },
            })
            .unwrap();
    };
    for property in tensor.material().properties() {
        add(
            property.requirement().name(),
            property.requirement().dims(),
            property.value_si(),
        );
    }
    for requirement in ThermalTensorPiezoelectric::thermal_bias_requirements().unwrap() {
        let name = requirement.name();
        let value = match name {
            "piezoelectric_thermal_expansion_mandel_3" => ALPHA,
            "piezoelectric_pyroelectric_clamped_3" => PYRO,
            _ => 0.0,
        };
        add(name, requirement.dims(), value);
    }
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic thermopyroelectric ceramic".into(),
            phase: "solid".into(),
            process: "explicit strain and charge reference".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    let tensor = TensorPiezoelectric::resolve(
        &card,
        &point(REFERENCE),
        MaterialPropertySelection::SingleClaimOnly,
        FRAME,
    )
    .unwrap();
    let law = ThermalTensorPiezoelectric::new(&card, &tensor, "T", [280.0, 340.0])
        .unwrap()
        .with_thermal_bias(REFERENCE, MaterialPropertySelection::SingleClaimOnly)
        .unwrap();
    (card, law)
}

#[test]
fn g1_clamped_and_free_thermal_strain_and_pyroelectric_boundaries() {
    with_cx(|cx, _| {
        let cube = cube(1);
        let n = cube.nodes.len();
        let (_, law) = source();
        let sample = law.at_temperature(HOT).unwrap();
        let identity = sample.thermal_bias_identity();
        let model = Model::assemble_thermal_materials(
            &cube.nodes,
            &cube.tets,
            vec![sample],
            FRAME,
            budget(),
            cx,
        )
        .unwrap();
        assert_eq!(
            model.thermal_materials().unwrap()[0].thermal_bias_identity(),
            identity
        );
        let fixed = (0..3 * n).map(|dof| (dof, 0.0)).collect::<Vec<_>>();
        let zero = vec![[0.0; 3]; n];
        let charge = vec![0.0; n];
        let a = ALPHA * (HOT - REFERENCE);
        let p = PYRO * (HOT - REFERENCE);
        let c = STIFFNESS[2];
        let epsilon = EPSILON[2];
        let open_strain = (c * a - COUPLING * p / epsilon) / (c + COUPLING * COUPLING / epsilon);
        let cases = [
            (fixed.as_slice(), Condition::Voltage(0.0), 0.0, 0.0),
            (fixed.as_slice(), Condition::Charge(0.0), 0.0, -p / epsilon),
            (cube.supports.as_slice(), Condition::Voltage(0.0), a, 0.0),
            (
                cube.supports.as_slice(),
                Condition::Charge(0.0),
                open_strain,
                -(COUPLING * open_strain + p) / epsilon,
            ),
        ];
        for (supports, condition, strain, field) in cases {
            let solved = model
                .solve(
                    supports,
                    &electrodes(&cube, 0.0, condition),
                    &zero,
                    &charge,
                    controls(),
                    cx,
                )
                .unwrap();
            let stress = c * (strain - a) - COUPLING * field;
            let displacement = COUPLING * strain + epsilon * field + p;
            for cell in &solved.cells {
                near(cell.strain_mandel[2], strain, 2e-12, 3e-8);
                near(cell.response.electric_field_v_m[2], field, 1e-5, 3e-8);
                near(cell.response.stress_mandel_pa[2], stress, 0.05, 3e-8);
                near(
                    cell.response.displacement_c_m2[2],
                    displacement,
                    2e-12,
                    3e-8,
                );
            }
            near(solved.electrodes[1].voltage_v, -field * LENGTH, 2e-8, 3e-8);
            near(
                solved.electrodes[1].charge_c,
                -displacement * LENGTH.powi(2),
                2e-18,
                3e-8,
            );
            let volume = LENGTH.powi(3);
            near(
                solved.stored_energy_j,
                0.5 * volume * (c * (strain - a).powi(2) + epsilon * field * field),
                1e-17,
                3e-8,
            );
            near(
                solved.reference_energy_j,
                0.5 * volume * c * a * a,
                1e-18,
                3e-8,
            );
            near(
                solved.mechanical_ramp_work_j,
                0.5 * volume * strain * (stress - c * a),
                1e-17,
                3e-8,
            );
            near(
                solved.electrical_ramp_work_j,
                0.5 * volume * field * (displacement - p),
                1e-17,
                3e-8,
            );
            near(solved.energy_residual_j, 0.0, 1e-16, 0.0);
        }
    });
}

#[test]
fn g1_biased_midpoint_oscillates_about_the_thermal_equilibrium() {
    with_cx(|cx, _| {
        let cube = cube(1);
        let n = cube.nodes.len();
        let (_, law) = source();
        let model = Model::assemble_thermal_materials(
            &cube.nodes,
            &cube.tets,
            vec![law.at_temperature(HOT).unwrap()],
            FRAME,
            budget(),
            cx,
        )
        .unwrap();
        let node = *cube.top.last().unwrap();
        let dof = 3 * node + 2;
        let fixed = (0..3 * n)
            .filter(|&i| i != dof)
            .map(|i| (i, 0.0))
            .collect::<Vec<_>>();
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
        let capacitance = cube
            .top
            .iter()
            .map(|&row| {
                let (cols, values) = model.dielectric().stiffness().row(row);
                cols.iter()
                    .zip(values)
                    .filter(|(col, _)| cube.top.contains(col))
                    .map(|(_, value)| value)
                    .sum::<f64>()
            })
            .sum::<f64>();
        // Independent weak integral of the initial thermal stress; the clamped
        // polarization charge is its exact planar surface integral.
        let reference_force = cube
            .tets
            .iter()
            .enumerate()
            .filter_map(|(cell, tet)| {
                tet.iter().position(|&id| id == node).map(|local| {
                    -STIFFNESS[2]
                        * ALPHA
                        * (HOT - REFERENCE)
                        * model.dielectric().volumes_m3()[cell]
                        * model.dielectric().gradients()[cell][local][2]
                })
            })
            .sum::<f64>();
        let reference_charge = -PYRO * (HOT - REFERENCE) * LENGTH.powi(2);
        let effective_stiffness = stiffness + coupling * coupling / capacitance;
        let equilibrium =
            (coupling * reference_charge / capacitance - reference_force) / effective_stiffness;
        let omega = (effective_stiffness / mass).sqrt();
        let amplitude = 1e-9;
        let zero = vec![[0.0; 3]; n];
        let charge = vec![0.0; n];
        let mut initial = zero.clone();
        initial[node][2] = equilibrium + amplitude;
        let mut body = Body::new(
            &model,
            &fixed,
            &electrodes(&cube, 0.0, Condition::Charge(0.0)),
            &initial,
            &zero,
            &charge,
            controls(),
            cx,
        )
        .unwrap();
        let energy = body.state().total_energy_j();
        for _ in 0..32 {
            let (next, step) = body
                .step_candidate(
                    0.7 / (32.0 * omega),
                    &[Condition::Voltage(0.0), Condition::Charge(0.0)],
                    &zero,
                    &charge,
                    controls(),
                    cx,
                )
                .unwrap()
                .into_parts();
            near(step.energy_residual_j, 0.0, 1e-16, 0.0);
            body = next;
        }
        near(
            (body.state().displacement_m()[node][2] - equilibrium) / amplitude,
            0.7_f64.cos(),
            4e-5,
            0.0,
        );
        near(
            body.state().velocity_m_s()[node][2] / (amplitude * omega),
            -0.7_f64.sin(),
            4e-5,
            0.0,
        );
        near(
            body.state().electrodes()[1].voltage_v,
            (coupling * body.state().displacement_m()[node][2] - reference_charge) / capacitance,
            2e-7,
            3e-8,
        );
        near(body.state().total_energy_j(), energy, 1e-16, 3e-8);
    });
}

#[test]
fn g3_zero_bias_is_exact_and_nonzero_bias_refuses_phasors_and_unreferenced_memories() {
    use fs_couple::piezoelectric::spatial::harmonic::{
        HarmonicElectrode, HarmonicElectrodeCondition as HCondition, HarmonicPiezoelectricControls,
    };
    use fs_math::c64::C64;
    with_cx(|cx, _| {
        let cube = cube(1);
        let n = cube.nodes.len();
        let (card, law) = source();
        let zero = vec![[0.0; 3]; n];
        let charge = vec![0.0; n];
        let sample = law.at_temperature(REFERENCE).unwrap();
        let tangents = [sample.material().clone()];
        let original =
            Model::assemble(&cube.nodes, &cube.tets, &tangents, FRAME, budget(), cx).unwrap();
        let unbiased = Model::assemble_thermal_materials(
            &cube.nodes,
            &cube.tets,
            vec![sample],
            FRAME,
            budget(),
            cx,
        )
        .unwrap();
        let electrodes = electrodes(&cube, 0.0, Condition::Voltage(-20.0));
        let a = original
            .solve(&cube.supports, &electrodes, &zero, &charge, controls(), cx)
            .unwrap();
        let b = unbiased
            .solve(&cube.supports, &electrodes, &zero, &charge, controls(), cx)
            .unwrap();
        assert_eq!(a.displacement_m, b.displacement_m);
        assert_eq!(a.voltage_v, b.voltage_v);
        assert_eq!(a.stored_energy_j, b.stored_energy_j);
        assert_eq!(a.mechanical_ramp_work_j, b.mechanical_ramp_work_j);
        assert_eq!(a.electrical_ramp_work_j, b.electrical_ramp_work_j);
        assert_eq!(b.reference_energy_j, 0.0);
        let mut final_state = Vec::new();
        for model in [&original, &unbiased] {
            let body = Body::new(
                model,
                &cube.supports,
                &electrodes,
                &zero,
                &zero,
                &charge,
                controls(),
                cx,
            )
            .unwrap();
            let (next, step) = body
                .step_candidate(
                    1e-8,
                    &[Condition::Voltage(0.0), Condition::Voltage(-20.0)],
                    &zero,
                    &charge,
                    controls(),
                    cx,
                )
                .unwrap()
                .into_parts();
            final_state.push((
                next.state().displacement_m().to_vec(),
                next.state().voltage_v().to_vec(),
                next.state().total_energy_j(),
                step.energy_residual_j,
            ));
        }
        assert_eq!(final_state[0], final_state[1]);
        let hot_sample = law.at_temperature(HOT).unwrap();
        let no_memory = RelaxingTensorPiezoelectric::resolve(
            hot_sample.material().clone(),
            &card,
            vec![],
            vec![],
            MaterialPropertySelection::SingleClaimOnly,
        )
        .unwrap();
        let hot = Model::assemble_thermal_materials(
            &cube.nodes,
            &cube.tets,
            vec![hot_sample],
            FRAME,
            budget(),
            cx,
        )
        .unwrap();
        let zero_phasor = C64::new(0.0, 0.0);
        let supports = cube
            .supports
            .iter()
            .map(|&(dof, value)| (dof, C64::new(value, 0.0)))
            .collect::<Vec<_>>();
        let harmonic_electrodes = [
            HarmonicElectrode {
                nodes: &cube.bottom,
                condition: HCondition::Voltage(zero_phasor),
            },
            HarmonicElectrode {
                nodes: &cube.top,
                condition: HCondition::Voltage(zero_phasor),
            },
        ];
        assert!(matches!(
            hot.solve_harmonic(
                1e5,
                &supports,
                &harmonic_electrodes,
                &[],
                &vec![[zero_phasor; 3]; n],
                &vec![zero_phasor; n],
                HarmonicPiezoelectricControls {
                    field: controls(),
                    power_tolerance_w: 1e-9,
                    maximum_real_coordinates: 512
                },
                cx
            ),
            Err(TetPiezoelectricError::Invalid(_))
        ));
        let body = Body::new(
            &hot,
            &cube.supports,
            &electrodes,
            &zero,
            &zero,
            &charge,
            controls(),
            cx,
        )
        .unwrap();
        assert!(matches!(
            body.with_relaxation(vec![no_memory], 0, controls(), cx),
            Err(TetPiezoelectricError::Invalid(
                "intrinsic memories require a declared thermal strain and polarization reference"
            ))
        ));
    });
}
