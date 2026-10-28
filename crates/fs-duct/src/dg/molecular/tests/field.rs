//! G0/G1: pointwise molecular reservoirs retain the original field mass and work.
use super::*;
use crate::profile::gas_field::tests::field;
use fs_matdb::QueryPoint;
use fs_material::gas_field::AxialGasField;

fn material(active: bool) -> RelaxationAbsorption {
    RelaxationAbsorption::new(
        ["analytic-field-a", "analytic-field-b"]
            .into_iter()
            .enumerate()
            .map(|(i, name)| SpeciesRelaxation {
                species: SpeciesId::new(name).unwrap(),
                processes: if active && i == 0 {
                    vec![RelaxationProcess {
                        reservoir: InternalReservoir::Vibrational(vec![VibrationalMode {
                            theta_k: 900.0,
                            degeneracy: 1,
                        }]),
                        rate: RelaxationRateLaw::PartnerIndependent {
                            p_tau_ref_pa_s: 2.0,
                            t_ref_k: 300.0,
                            exponent: 1.2,
                            t_min_k: 200.0,
                            t_max_k: 500.0,
                        },
                    }]
                } else {
                    vec![]
                },
            })
            .collect(),
    )
    .unwrap()
}

fn evaluated(field: &AxialGasField, law: &RelaxationAbsorption, x: f64) -> RelaxationState {
    law.evaluate(&RelaxingMedium::from_thermochemical(&field.state_at(x).unwrap()).unwrap())
        .unwrap()
}

fn shape(field: &AxialGasField, cone: bool) -> Duct {
    Duct {
        segments: vec![if cone {
            Segment::Cone {
                inlet_radius: 0.001,
                outlet_radius: 0.002,
                length: field.length_m(),
            }
        } else {
            Segment::Cylinder {
                radius: 0.001,
                length: field.length_m(),
            }
        }],
    }
}

fn resolution(points: usize) -> DgWallResolution {
    DgWallResolution {
        poles_per_field: 8,
        quadrature_points: points,
        band_rad_s: [50.0, 2000.0].map(|v| core::f64::consts::TAU * v),
        maximum_relative_error: 0.05,
    }
}

fn base(field: &AxialGasField, cone: bool, elements: usize, points: usize) -> DgDuct {
    DgDuct::new_with_gas_field(
        &shape(field, cone),
        field,
        DgMesh { elements, order: 1 },
        DgBoundary::Port(PortId(7)),
        DgBoundary::Closed,
    )
    .unwrap()
    .with_spatial_viscothermal(resolution(points))
    .unwrap()
}

// Physical weak compliance from an independent 64-point mass integral and the
// declared local heat equation. It does not consume retained projections, mass
// roots, reservoir coefficients, or the pressure mass used by production.
fn physical_pressure_oracle(
    field: &AxialGasField,
    law: &RelaxationAbsorption,
    duct: &DgDuct,
    omega: f64,
    points: usize,
    cone: bool,
) -> Vec<C64> {
    let m = duct.field_dofs;
    let local = duct.mesh.order + 1;
    let s = C64::new(0.0, -omega);
    let mut matrix = vec![C64::ZERO; m * m];
    for (e, element) in duct.elements.iter().enumerate() {
        let jac = 0.5 * element.length;
        for (count, correction) in [(64, false), (points, true)] {
            let (nodes, weights) = gauss_legendre_rule(count);
            for (&xi, &weight) in nodes.iter().zip(&weights) {
                let x = element.x_left + jac * (xi + 1.0);
                let radius = 0.001
                    * if cone {
                        1.0 + x / field.length_m()
                    } else {
                        1.0
                    };
                let state = evaluated(field, law, x);
                let gas = state.gas();
                let frozen = 1.0 / (gas.density * state.frozen_sound_speed_m_s().powi(2));
                let value = if correction {
                    let cp = state.processes().iter().fold(
                        C64::from_re(state.frozen_cv_j_kg_k() + gas.specific_gas_constant),
                        |sum, p| {
                            sum + C64::from_re(p.heat_capacity_j_kg_k)
                                / (C64::ONE + s.scale(p.relaxation_time_s))
                        },
                    );
                    let r = resolution(points);
                    let bessel = fs_phs::ZkDiffusive::new(
                        radius,
                        gas.density,
                        gas.dynamic_viscosity,
                        gas.gamma,
                        gas.prandtl,
                        r.poles_per_field,
                        r.band_rad_s,
                        r.maximum_relative_error,
                    )
                    .unwrap();
                    let heat = bessel.shunt_terms().iter().fold(C64::ZERO, |sum, t| {
                        let w = t.gain_per_s / (t.rate_per_s * (gas.gamma - 1.0));
                        sum + s.scale(w)
                            / (s * cp + C64::from_re(gas.specific_heat_cp * t.rate_per_s))
                    });
                    (C64::ONE - heat.scale(gas.specific_gas_constant))
                        .scale(1.0 / (gas.density * gas.specific_gas_constant * gas.temperature))
                        - C64::from_re(frozen)
                } else {
                    C64::from_re(frozen)
                };
                let coefficient =
                    (s * value).scale(jac * weight * core::f64::consts::PI * radius * radius);
                let phi = duct.lgl.eval_basis(xi);
                for i in 0..local {
                    for j in 0..local {
                        let k = (e * local + i) * m + e * local + j;
                        matrix[k] = matrix[k] + coefficient.scale(phi[i] * phi[j]);
                    }
                }
            }
        }
    }
    matrix
}

fn physical_pressure(duct: &DgDuct, omega: f64) -> Vec<C64> {
    let normalized = super::spatial::pressure_schur(duct, omega);
    let m = duct.field_dofs;
    let local = duct.mesh.order + 1;
    let mut out = vec![C64::ZERO; m * m];
    for (e, element) in duct.elements.iter().enumerate() {
        for i in 0..local {
            for j in 0..local {
                let mut value = C64::ZERO;
                for k in 0..local {
                    for l in 0..local {
                        value = value
                            + normalized[(e * local + k) * m + e * local + l].scale(
                                element.mp_sqrt[i * local + k] * element.mp_sqrt[l * local + j],
                            );
                    }
                }
                out[(e * local + i) * m + e * local + j] = value;
            }
        }
    }
    out
}

fn relative(a: &[C64], b: &[C64]) -> f64 {
    a.iter()
        .zip(b)
        .map(|(a, b)| (*a - *b).abs())
        .fold(0.0, f64::max)
        / b.iter().map(|v| v.abs()).fold(0.0, f64::max)
}

#[test]
fn field_molecular_schur_uses_local_material_and_exact_frozen_mass() {
    let field = field([270.0, 370.0], [0.95, 0.05], 101325.0);
    let law = material(true);
    for points in [2, 4, 8] {
        let original = base(&field, true, 1, points);
        let total = original.molecular_relaxation_field_dofs(&law).unwrap();
        assert_eq!(total, original.total_dofs() + points * 8);
        let duct = original
            .with_molecular_relaxation_field(law.clone())
            .unwrap();
        assert_eq!(duct.total_dofs(), total);
        assert!(duct.molecular_relaxation().is_none());
        assert_eq!(duct.molecular_relaxation_field(), Some(&law));
        for omega in [350.0, 1500.0, 7000.0] {
            let actual = physical_pressure(&duct, omega);
            let expected = physical_pressure_oracle(&field, &law, &duct, omega, points, true);
            assert!(
                relative(&actual, &expected) < 3e-10,
                "{points} {omega}: {}",
                relative(&actual, &expected)
            );
        }
    }
}

#[test]
fn field_molecular_quadrature_refines_the_original_local_heat_equation() {
    let field = field([270.0, 370.0], [0.95, 0.05], 101325.0);
    let law = material(true);
    let mut errors = Vec::new();
    for points in [2, 4, 8] {
        let duct = base(&field, true, 1, points)
            .with_molecular_relaxation_field(law.clone())
            .unwrap();
        let actual = physical_pressure(&duct, 3000.0);
        let expected = physical_pressure_oracle(&field, &law, &duct, 3000.0, 64, true);
        errors.push(relative(&actual, &expected));
    }
    assert!(
        errors[1] < errors[0] * 0.1 && errors[2] < errors[1] * 0.1,
        "{errors:?}"
    );
    assert!(errors[2] < 1e-7, "{errors:?}");
}

#[test]
fn field_molecular_shared_loss_ports_and_structural_wall_remain_physical() {
    let field = field([290.0, 340.0], [0.9, 0.1], 101325.0);
    let law = material(true);
    let original = base(&field, false, 2, 3);
    let junction = [DgSidePort {
        id: PortId(11),
        downstream_element: 1,
    }];
    let a = original
        .clone()
        .with_side_ports(&junction)
        .unwrap()
        .with_molecular_relaxation_field(law.clone())
        .unwrap();
    let b = original
        .with_molecular_relaxation_field(law)
        .unwrap()
        .with_side_ports(&junction)
        .unwrap();
    for ((id_a, a), (id_b, b)) in a.ports.iter().zip(&b.ports) {
        assert_eq!(id_a, id_b);
        let scale = b.iter().map(|v| v.abs()).fold(0.0, f64::max);
        assert!(a.iter().zip(b).all(|(a, b)| (a - b).abs() < 5e-13 * scale));
    }
    let duct = a
        .with_compliant_wall(fs_phs::WallPin {
            surface_density: 0.2,
            stiffness_per_area: 4e7,
            resistance: 300.0,
        })
        .unwrap();
    let mut state = duct.zero_state();
    for (i, value) in state.tilde.iter_mut().enumerate() {
        *value = 1e-4 * det::sin(0.37 * i as f64);
    }
    let n = duct.total_dofs();
    let derivative: f64 = (0..n)
        .map(|i| {
            state.tilde[i]
                * (0..n)
                    .map(|j| duct.l_tilde[i * n + j] * state.tilde[j])
                    .sum::<f64>()
        })
        .sum();
    let wall = duct.wall_dissipated_power_w(&state).unwrap();
    let molecular = duct.molecular_dissipated_power_w(&state).unwrap();
    assert!(wall > 0.0 && molecular > 0.0);
    assert!((derivative + wall + molecular).abs() < 3e-11 * (wall + molecular));
    let response = duct
        .harmonic_response(
            3000.0,
            &[(PortId(7), C64::from_re(1e-8))],
            1_000_000_000,
            1e-9,
            &mut || false,
        )
        .unwrap();
    assert!(response.power_residual_w.abs() < 2e-8 * response.input_power_w.abs());
    assert!(response.molecular_dissipation_w > 0.0);
}

#[test]
fn field_molecular_uniform_and_zero_laws_keep_original_operators() {
    let field = field([300.0; 2], [0.8; 2], 101325.0);
    let law = material(true);
    let state = evaluated(&field, &law, 0.0);
    let bare = DgDuct::new_with_gas_field(
        &shape(&field, false),
        &field,
        DgMesh {
            elements: 1,
            order: 1,
        },
        DgBoundary::Closed,
        DgBoundary::Closed,
    )
    .unwrap();
    for base in [
        bare.clone(),
        bare.clone()
            .with_viscothermal(8, [50.0, 2000.0], 0.05)
            .unwrap(),
        bare.with_spatial_viscothermal(resolution(3)).unwrap(),
    ] {
        let expected = base
            .clone()
            .with_molecular_relaxation(state.clone())
            .unwrap();
        assert_eq!(
            base.molecular_relaxation_field_dofs(&law).unwrap(),
            expected.total_dofs()
        );
        let actual = base.with_molecular_relaxation_field(law.clone()).unwrap();
        assert_eq!(actual.l_tilde, expected.l_tilde);
        assert_eq!(actual.ports, expected.ports);
        assert_eq!(actual.operator_identity(), expected.operator_identity());
    }
    let varying = crate::profile::gas_field::tests::field([290.0, 340.0], [0.9, 0.1], 101325.0);
    let base = base(&varying, true, 1, 3);
    let expected = base.clone();
    assert_eq!(
        base.molecular_relaxation_field_dofs(&material(false))
            .unwrap(),
        base.total_dofs()
    );
    let actual = base
        .with_molecular_relaxation_field(material(false))
        .unwrap();
    assert_eq!(actual.l_tilde, expected.l_tilde);
    assert_eq!(actual.operator_identity(), expected.operator_identity());
}

#[test]
fn field_molecular_counts_local_zero_reservoirs_and_refuses_invalid_domains() {
    let source = field([300.0; 2], [0.0, 0.8], 101325.0);
    let length = source.length_m();
    let field = AxialGasField::new(
        vec![0.0, length / 3.0, 2.0 * length / 3.0, length],
        101325.0,
        vec![
            source.primitive_states()[0].clone(),
            source.primitive_states()[0].clone(),
            source.primitive_states()[1].clone(),
            source.primitive_states()[1].clone(),
        ],
        source.model().clone(),
        QueryPoint::new(),
        source.phase(),
    )
    .unwrap();
    let original = base(&field, false, 3, 2);
    let expected = original.total_dofs() + 2 * 2 * 8;
    assert_eq!(
        original
            .molecular_relaxation_field_dofs(&material(true))
            .unwrap(),
        expected
    );
    let duct = original
        .clone()
        .with_molecular_relaxation_field(material(true))
        .unwrap();
    assert_eq!(duct.total_dofs(), expected);
    let memory = duct.molecular.as_ref().unwrap();
    assert_eq!(
        memory.modes.iter().filter(|m| m.rates.is_empty()).count(),
        2 * 8
    );
    assert!(memory.modes.iter().all(|m| m.heat_rate > 0.0));
    assert!(
        duct.clone()
            .with_molecular_relaxation_field(material(true))
            .is_err()
    );
    assert!(
        original
            .clone()
            .with_molecular_relaxation_field_controlled(material(true), &mut || true)
            .is_err()
    );
    let bare = DgDuct::new_with_gas_field(
        &shape(&field, false),
        &field,
        DgMesh {
            elements: 3,
            order: 1,
        },
        DgBoundary::Closed,
        DgBoundary::Closed,
    )
    .unwrap();
    assert!(
        bare.with_molecular_relaxation_field(material(true))
            .is_err()
    );
    let hot = crate::profile::gas_field::tests::field([300.0, 550.0], [0.8; 2], 101325.0);
    assert!(
        base(&hot, false, 1, 2)
            .molecular_relaxation_field_dofs(&material(true))
            .is_err()
    );
    let mut live = duct;
    let mut state = live.zero_state();
    let before = state.clone();
    assert!(
        live.rebind_spatial_viscothermal_gas_field(&mut state, &field, 0.05, 1.0, &mut || false)
            .is_err()
    );
    assert_eq!(state, before);
}
