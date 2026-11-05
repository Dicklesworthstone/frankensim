//! Independent finite-body Maxwell/Debye oracles using authored SI sources.
//! The base roles are deliberately relaxed C^E and frozen epsilon^S. These
//! fixtures do not infer a loss spectrum from a named ceramic or a measured Q.
use super::*;
use fs_material::field_storage::ELECTRIC_PERMITTIVITY_DIMS;
use fs_material::piezoelectric::tensor::relaxation::{
    DirectedDebyeProperties, RelaxingTensorPiezoelectric,
};
use fs_material::visco_tensor::source::TensorBranchProperties;
use fs_qty::{Dims, Pressure, Time};

pub(crate) fn fixture(
    mechanical: Option<(f64, f64)>,
    dielectric: Option<(f64, f64)>,
    frame: ContentHash,
) -> (TensorPiezoelectric, RelaxingTensorPiezoelectric) {
    let reference = material([0.0, 0.0, 12.0, 0.0, 0.0, 0.0]);
    let mut values: Vec<_> = reference
        .material()
        .properties()
        .iter()
        .map(|property| {
            (
                property.requirement().name().to_owned(),
                property.requirement().dims(),
                property.value_si(),
            )
        })
        .collect();
    let branches = mechanical
        .map(|(fraction, tau)| {
            let mut names = Vec::with_capacity(21);
            for (row, &stiffness) in STIFFNESS.iter().enumerate() {
                for col in row..6 {
                    let name = format!("memory_c_{}{}", row + 1, col + 1);
                    let value = if row == col {
                        fraction * stiffness
                    } else {
                        0.0
                    };
                    values.push((name.clone(), Pressure::DIMS, value));
                    names.push(name);
                }
            }
            values.push(("memory_tau_m".into(), Time::DIMS, tau));
            TensorBranchProperties {
                stiffness_upper: names.try_into().unwrap(),
                relaxation_time: "memory_tau_m".into(),
            }
        })
        .into_iter()
        .collect();
    let poles = dielectric
        .map(|(increment, tau)| {
            values.push((
                "memory_delta_epsilon".into(),
                ELECTRIC_PERMITTIVITY_DIMS,
                increment,
            ));
            values.push(("memory_tau_p".into(), Time::DIMS, tau));
            let direction = ["memory_nx", "memory_ny", "memory_nz"].map(str::to_owned);
            for (name, value) in direction.iter().zip([0.0, 0.0, 1.0]) {
                values.push((name.clone(), Dims::NONE, value));
            }
            DirectedDebyeProperties {
                permittivity_increment: "memory_delta_epsilon".into(),
                relaxation_time: "memory_tau_p".into(),
                direction,
            }
        })
        .into_iter()
        .collect();
    let mut claims = ClaimSet::new();
    for (name, dims, value) in values {
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value: PropertyValue::Scalar { value, dims },
                validity: ValidityDomain::unconstrained().with("T", 290.0, 310.0),
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::ConstantWithinValidity,
                observations: vec![],
                provenance: Provenance {
                    source:
                        "authored relaxed stiffness, frozen permittivity and passive memory fixture"
                            .into(),
                    license: "CC0-1.0".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic reciprocal ceramic with memory".into(),
            phase: "solid".into(),
            process: "explicit Maxwell and directed Debye source coefficients".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    // Resolve both owners from this same COMPLETE card, including loss roles.
    let base = TensorPiezoelectric::resolve(
        &card,
        &QueryPoint::new().with("T", 300.0).unwrap(),
        MaterialPropertySelection::SingleClaimOnly,
        frame,
    )
    .unwrap();
    let loss = RelaxingTensorPiezoelectric::resolve(
        base.clone(),
        &card,
        branches,
        poles,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    (base, loss)
}

struct OneCoordinate {
    node: usize,
    fixed: Vec<(usize, C64)>,
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
            .map(|i| (i, ZERO))
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

// Closed scalar expressions, independent of the material harmonic evaluator.
fn maxwell_stiffness(k: f64, fraction: f64, omega_tau: f64) -> C64 {
    let denominator = 1.0 + omega_tau * omega_tau;
    C64::new(
        k * (1.0 + fraction * omega_tau * omega_tau / denominator),
        -k * fraction * omega_tau / denominator,
    )
}

fn debye_capacitance(increment: f64, omega_tau: f64) -> C64 {
    let c = increment * LENGTH / (1.0 + omega_tau * omega_tau);
    C64::new(EPSILON[2] * LENGTH + c, c * omega_tau)
}

#[test]
fn g1_clamped_debye_current_drive_matches_capacitance_phase_storage_and_heat() {
    with_cx(|cx, _| {
        let cube = cube(1);
        let n = cube.nodes.len();
        let increment = 0.8 * EPSILON[2];
        let tau = 1e-6;
        let (base, loss) = fixture(None, Some((increment, tau)), FRAME);
        let materials = [base];
        let losses = [loss];
        let model =
            Model::assemble(&cube.nodes, &cube.tets, &materials, FRAME, budget(), cx).unwrap();
        let fixed: Vec<_> = (0..3 * n).map(|dof| (dof, ZERO)).collect();
        let current = C64::new(3e-5, -4e-5);
        for x in [0.05, 1.0, 20.0] {
            let omega = x / tau;
            let capacitance = debye_capacitance(increment, x);
            let voltage = div(current, mul(C64::new(0.0, -omega), capacitance));
            let heat = 0.5 * omega * capacitance.im * voltage.norm_sq();
            let response = model
                .solve_harmonic_with_relaxation(
                    omega,
                    &losses,
                    &fixed,
                    &helectrodes(&cube, ZERO, HCondition::Current(current)),
                    &[],
                    &vec![[ZERO; 3]; n],
                    &vec![ZERO; n],
                    hc(),
                    cx,
                )
                .unwrap();
            cnear(response.electrodes[1].voltage_v, voltage, 1e-10, 3e-9);
            cnear(
                response.electrodes[1].charge_c,
                mul(capacitance, voltage),
                1e-20,
                3e-9,
            );
            cnear(
                response.electrodes[1].source_current_a,
                current,
                1e-14,
                3e-9,
            );
            assert!(heat > 0.0);
            near(response.average_material_heat_w, heat, 1e-13, 3e-9);
            near(response.average_electrical_power_w, heat, 1e-13, 3e-9);
            assert_eq!(response.average_load_heat_w, 0.0);
            assert_eq!(response.average_mechanical_power_w, 0.0);
            near(
                response.average_stored_energy_j,
                0.25 * capacitance.re * voltage.norm_sq(),
                1e-20,
                3e-9,
            );
            let memory = 0.25 * increment * LENGTH * voltage.norm_sq() / (1.0 + x * x);
            near(
                response
                    .cells
                    .iter()
                    .map(|c| c.average_memory_energy_j)
                    .sum(),
                memory,
                1e-20,
                3e-9,
            );
            for cell in &response.cells {
                assert_eq!(cell.average_mechanical_heat_w, 0.0);
                assert!(cell.average_dielectric_heat_w > 0.0);
                cnear(
                    cell.electric_field_v_m[2],
                    scale(voltage, -1.0 / LENGTH),
                    1e-7,
                    3e-9,
                );
                cnear(
                    cell.displacement_c_m2[2],
                    scale(mul(capacitance, voltage), -1.0 / (LENGTH * LENGTH)),
                    1e-17,
                    3e-9,
                );
            }
            near(
                response
                    .cells
                    .iter()
                    .map(|c| c.average_dielectric_heat_w)
                    .sum(),
                heat,
                1e-13,
                3e-9,
            );
            near(response.power_balance_residual_w, 0.0, 1e-13, 0.0);
        }
    });
}

#[test]
fn g1_maxwell_memory_resolves_the_mechanical_pole_and_accounts_for_dissipation() {
    with_cx(|cx, _| {
        let cube = cube(1);
        let n = cube.nodes.len();
        let fraction = 0.5;
        let tau = 2e-7;
        let (base, loss) = fixture(Some((fraction, tau)), None, FRAME);
        let materials = [base];
        let losses = [loss];
        let model =
            Model::assemble(&cube.nodes, &cube.tets, &materials, FRAME, budget(), cx).unwrap();
        let coordinate = one_coordinate(&model, &cube);
        // Solve Re(K*(omega))-omega^2 M=0 analytically. The remaining
        // imaginary stiffness comes entirely from the declared Maxwell branch.
        let a = coordinate.stiffness * tau * tau / coordinate.mass;
        let b = 1.0 - a * (1.0 + fraction);
        let x_squared = 0.5 * ((b * b + 4.0 * a).sqrt() - b);
        let x = x_squared.sqrt();
        let omega = x / tau;
        let stiffness = maxwell_stiffness(coordinate.stiffness, fraction, x);
        let dynamic = add(stiffness, C64::new(-omega * omega * coordinate.mass, 0.0));
        assert!(dynamic.re.abs() <= 1e-12 * coordinate.stiffness);
        assert!(dynamic.im < 0.0);
        let force = C64::new(0.02, 0.01);
        let motion = div(force, dynamic);
        let velocity = mul(C64::new(0.0, -omega), motion);
        let heat = -0.5 * omega * stiffness.im * motion.norm_sq();
        let mut forces = vec![[ZERO; 3]; n];
        forces[coordinate.node][2] = force;
        let response = model
            .solve_harmonic_with_relaxation(
                omega,
                &losses,
                &coordinate.fixed,
                &helectrodes(&cube, ZERO, HCondition::Voltage(ZERO)),
                &[],
                &forces,
                &vec![ZERO; n],
                hc(),
                cx,
            )
            .unwrap();
        cnear(
            response.displacement_m[coordinate.node][2],
            motion,
            1e-18,
            3e-8,
        );
        cnear(
            response.velocity_m_s[coordinate.node][2],
            velocity,
            1e-12,
            3e-8,
        );
        cnear(
            response.electrodes[1].charge_c,
            scale(motion, -coordinate.coupling),
            1e-20,
            3e-8,
        );
        near(response.average_material_heat_w, heat, 1e-13, 3e-8);
        near(
            response.average_mechanical_power_w,
            0.5 * mul(force, velocity.conj()).re,
            1e-13,
            3e-8,
        );
        near(response.average_mechanical_power_w, heat, 1e-13, 3e-8);
        near(
            response.average_stored_energy_j,
            0.25 * stiffness.re * motion.norm_sq(),
            1e-20,
            3e-8,
        );
        near(
            response.average_kinetic_energy_j,
            0.25 * coordinate.mass * velocity.norm_sq(),
            1e-20,
            3e-8,
        );
        let memory = 0.25 * fraction * coordinate.stiffness * x_squared / (1.0 + x_squared)
            * motion.norm_sq();
        near(
            response
                .cells
                .iter()
                .map(|c| c.average_memory_energy_j)
                .sum(),
            memory,
            1e-20,
            3e-8,
        );
        near(
            response
                .cells
                .iter()
                .map(|c| c.average_mechanical_heat_w)
                .sum(),
            heat,
            1e-13,
            3e-8,
        );
        assert!(
            response
                .cells
                .iter()
                .all(|c| c.average_dielectric_heat_w == 0.0)
        );
        assert_eq!(response.average_load_heat_w, 0.0);
        assert_eq!(response.average_electrical_power_w, 0.0);
        near(response.power_balance_residual_w, 0.0, 1e-13, 0.0);
    });
}

#[test]
fn g1_combined_intrinsic_loss_and_resistor_match_the_coupled_scalar_power_balance() {
    with_cx(|cx, _| {
        let cube = cube(1);
        let n = cube.nodes.len();
        let fraction = 0.4;
        let tau_m = 2e-7;
        let increment = 0.75 * EPSILON[2];
        let tau_p = 4e-7;
        let (base, loss) = fixture(Some((fraction, tau_m)), Some((increment, tau_p)), FRAME);
        let materials = [base];
        let losses = [loss];
        let model =
            Model::assemble(&cube.nodes, &cube.tets, &materials, FRAME, budget(), cx).unwrap();
        let coordinate = one_coordinate(&model, &cube);
        let omega = 0.8 * (coordinate.stiffness / coordinate.mass).sqrt();
        let resistance = 1.0 / (omega * EPSILON[2] * LENGTH);
        let stiffness = maxwell_stiffness(coordinate.stiffness, fraction, omega * tau_m);
        let capacitance = debye_capacitance(increment, omega * tau_p);
        let admittance = add(
            C64::new(1.0 / resistance, 0.0),
            mul(C64::new(0.0, -omega), capacitance),
        );
        let voltage_per_motion = div(C64::new(0.0, -omega * coordinate.coupling), admittance);
        let dynamic = add(
            add(stiffness, C64::new(-omega * omega * coordinate.mass, 0.0)),
            scale(voltage_per_motion, coordinate.coupling),
        );
        let force = C64::new(0.02, 0.01);
        let motion = div(force, dynamic);
        let voltage = mul(voltage_per_motion, motion);
        let mechanical_heat = -0.5 * omega * stiffness.im * motion.norm_sq();
        let dielectric_heat = 0.5 * omega * capacitance.im * voltage.norm_sq();
        let resistor_heat = 0.5 * voltage.norm_sq() / resistance;
        assert!(mechanical_heat > 0.0 && dielectric_heat > 0.0 && resistor_heat > 0.0);
        let mut forces = vec![[ZERO; 3]; n];
        forces[coordinate.node][2] = force;
        let response = model
            .solve_harmonic_with_relaxation(
                omega,
                &losses,
                &coordinate.fixed,
                &helectrodes(&cube, ZERO, HCondition::Current(ZERO)),
                &[HLoad::Resistive(super::super::transient_load::resistor(
                    resistance,
                    [1, 0],
                ))],
                &forces,
                &vec![ZERO; n],
                hc(),
                cx,
            )
            .unwrap();
        cnear(
            response.displacement_m[coordinate.node][2],
            motion,
            1e-18,
            3e-8,
        );
        cnear(response.electrodes[1].voltage_v, voltage, 1e-10, 3e-8);
        cnear(response.electrodes[1].source_current_a, ZERO, 1e-13, 0.0);
        cnear(
            response.electrodes[1].charge_c,
            add(
                mul(capacitance, voltage),
                scale(motion, -coordinate.coupling),
            ),
            1e-20,
            3e-8,
        );
        near(
            response
                .cells
                .iter()
                .map(|c| c.average_mechanical_heat_w)
                .sum(),
            mechanical_heat,
            1e-13,
            3e-8,
        );
        near(
            response
                .cells
                .iter()
                .map(|c| c.average_dielectric_heat_w)
                .sum(),
            dielectric_heat,
            1e-13,
            3e-8,
        );
        near(
            response.average_material_heat_w,
            mechanical_heat + dielectric_heat,
            1e-13,
            3e-8,
        );
        near(response.average_load_heat_w, resistor_heat, 1e-13, 3e-8);
        near(
            response.average_mechanical_power_w,
            mechanical_heat + dielectric_heat + resistor_heat,
            1e-13,
            3e-8,
        );
        near(response.average_electrical_power_w, 0.0, 1e-13, 0.0);
        near(
            response.average_stored_energy_j,
            0.25 * (stiffness.re * motion.norm_sq() + capacitance.re * voltage.norm_sq()),
            1e-20,
            3e-8,
        );
        let memory = 0.25
            * ((stiffness.re - coordinate.stiffness) * motion.norm_sq()
                + (capacitance.re - EPSILON[2] * LENGTH) * voltage.norm_sq());
        near(
            response
                .cells
                .iter()
                .map(|c| c.average_memory_energy_j)
                .sum(),
            memory,
            1e-20,
            3e-8,
        );
        near(response.power_balance_residual_w, 0.0, 1e-13, 0.0);
    });
}

#[test]
fn g3_empty_relaxation_preserves_the_existing_harmonic_solution_exactly() {
    with_cx(|cx, _| {
        let cube = cube(1);
        let n = cube.nodes.len();
        let (base, loss) = fixture(None, None, FRAME);
        let materials = [base];
        let losses = [loss];
        let model =
            Model::assemble(&cube.nodes, &cube.tets, &materials, FRAME, budget(), cx).unwrap();
        let coordinate = one_coordinate(&model, &cube);
        let omega = 0.7 * (coordinate.stiffness / coordinate.mass).sqrt();
        let electrodes = helectrodes(&cube, ZERO, HCondition::Current(C64::new(2e-5, -3e-5)));
        let loads = [HLoad::Resistive(super::super::transient_load::resistor(
            1e4,
            [1, 0],
        ))];
        let mut forces = vec![[ZERO; 3]; n];
        forces[coordinate.node][2] = C64::new(0.02, 0.01);
        let old = model
            .solve_harmonic(
                omega,
                &coordinate.fixed,
                &electrodes,
                &loads,
                &forces,
                &vec![ZERO; n],
                hc(),
                cx,
            )
            .unwrap();
        let new = model
            .solve_harmonic_with_relaxation(
                omega,
                &losses,
                &coordinate.fixed,
                &electrodes,
                &loads,
                &forces,
                &vec![ZERO; n],
                hc(),
                cx,
            )
            .unwrap();
        assert_eq!(new.displacement_m, old.displacement_m);
        assert_eq!(new.velocity_m_s, old.velocity_m_s);
        assert_eq!(new.voltage_v, old.voltage_v);
        assert_eq!(new.support_force_n, old.support_force_n);
        for (actual, expected) in new.electrodes.iter().zip(&old.electrodes) {
            assert_eq!(actual.voltage_v, expected.voltage_v);
            assert_eq!(actual.charge_c, expected.charge_c);
            assert_eq!(actual.source_current_a, expected.source_current_a);
        }
        for (actual, expected) in new.cells.iter().zip(&old.cells) {
            assert_eq!(actual.strain_mandel, expected.strain_mandel);
            assert_eq!(actual.stress_mandel_pa, expected.stress_mandel_pa);
            assert_eq!(actual.electric_field_v_m, expected.electric_field_v_m);
            assert_eq!(actual.displacement_c_m2, expected.displacement_c_m2);
            assert_eq!(actual.average_memory_energy_j, 0.0);
            assert_eq!(actual.average_mechanical_heat_w, 0.0);
            assert_eq!(actual.average_dielectric_heat_w, 0.0);
        }
        assert_eq!(new.average_stored_energy_j, old.average_stored_energy_j);
        assert_eq!(new.average_kinetic_energy_j, old.average_kinetic_energy_j);
        assert_eq!(
            new.average_mechanical_power_w,
            old.average_mechanical_power_w
        );
        assert_eq!(
            new.average_electrical_power_w,
            old.average_electrical_power_w
        );
        assert_eq!(new.average_load_heat_w, old.average_load_heat_w);
        assert_eq!(new.average_material_heat_w, 0.0);
        assert_eq!(new.power_balance_residual_w, old.power_balance_residual_w);
    });
}

#[test]
fn g4_relaxation_rejects_foreign_sources_frames_cycle_overruns_and_cancellation() {
    with_cx(|cx, gate| {
        let cube = cube(1);
        let n = cube.nodes.len();
        let poles = Some((0.8 * EPSILON[2], 1e-6));
        let (base, loss) = fixture(Some((0.4, 2e-7)), poles, FRAME);
        // Identical base coefficients with a distinct original card still do
        // not authorize replacing the model's complete material source.
        let (_, foreign_source) = fixture(Some((0.5, 2e-7)), poles, FRAME);
        let (_, foreign_frame) = fixture(Some((0.4, 2e-7)), poles, ContentHash([81; 32]));
        let materials = [base];
        let losses = [loss];
        let model =
            Model::assemble(&cube.nodes, &cube.tets, &materials, FRAME, budget(), cx).unwrap();
        let fixed: Vec<_> = (0..3 * n).map(|i| (i, ZERO)).collect();
        let solve = |laws: &[RelaxingTensorPiezoelectric], voltage| {
            model.solve_harmonic_with_relaxation(
                1e6,
                laws,
                &fixed,
                &helectrodes(&cube, ZERO, HCondition::Voltage(voltage)),
                &[],
                &vec![[ZERO; 3]; n],
                &vec![ZERO; n],
                hc(),
                cx,
            )
        };
        for laws in [&[][..], &[foreign_source][..], &[foreign_frame][..]] {
            assert!(matches!(
                solve(laws, C64::new(1.0, 0.0)),
                Err(TetPiezoelectricError::Invalid(_))
            ));
        }
        let quadrature = 0.8 * hc().field.maximum_electric_field_v_m * LENGTH;
        assert!(matches!(
            solve(&losses, C64::new(quadrature, quadrature)),
            Err(TetPiezoelectricError::Invalid(_))
        ));
        assert!(solve(&losses, C64::new(1.0, 0.0)).is_ok());
        gate.request();
        assert!(matches!(
            solve(&losses, C64::new(1.0, 0.0)),
            Err(TetPiezoelectricError::Cancelled)
        ));
    });
}
