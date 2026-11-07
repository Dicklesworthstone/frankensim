//! Compact original calorics in the native spatial phase and conductivity solve.

use super::*;
use fs_conduction::transient::enthalpy::EnthalpyRadiation;
use fs_conduction::transient::enthalpy::material::MaterialEnthalpyError;
use fs_conduction::transient::enthalpy::material::melting::ExactMaterialMeltingSource;
use fs_conduction::transient::enthalpy::sdirk::EnthalpySdirk2Input;

fn resolve_capacity(
    card: &MaterialCard,
    closure: CoexistenceConductivity,
    datum: f64,
    maximum_chord: Option<f64>,
    fusion_rounding: f64,
) -> Result<MaterialMeltingTransport, MaterialEnthalpyError> {
    let cp = PropertyKey::new("cp", CP);
    let rho = PropertyKey::new("rho", RHO);
    let k = PropertyKey::new("k", K);
    let tm = PropertyKey::with_quantity("Tm", absolute_temperature());
    let fusion = PropertyKey::new("fusion", ENERGY);
    let points = [
        point(300.0, Some(0.0), 100_000.0),
        point(350.0, Some(0.0), 100_000.0),
        point(350.0, Some(1.0), 100_000.0),
        point(450.0, Some(1.0), 100_000.0),
        point(350.0, None, 100_000.0),
    ];
    let branch = |index| HeatCapacitySourceQuery {
        heat_capacity: &cp,
        density: &rho,
        temperature_axis: "T",
        lower: &points[index],
        upper: &points[index + 1],
        selection: SelectionPolicy::SingleClaimOnly,
    };
    let caloric = MeltingSourceQuery {
        solid: branch(0),
        liquid: branch(2),
        melting_point: &tm,
        latent_heat_fusion: &fusion,
        transition: &points[4],
        transition_selection: SelectionPolicy::SingleClaimOnly,
        phase_selector: Some("phase"),
    };
    match maximum_chord {
        Some(maximum_error) => MaterialMeltingTransport::from_isotropic_card(
            card,
            MaterialMeltingSource {
                caloric,
                conductivity: [&k, &k],
                reference_phase: SolidLiquidPhase::Solid,
                reference: &points[0],
                lower_specific_enthalpy_j_kg: datum,
                maximum_interpolation_error_j_kg: maximum_error,
                coexistence: closure,
            },
        ),
        None => MaterialMeltingTransport::from_isotropic_card_exact(
            card,
            ExactMaterialMeltingSource {
                caloric,
                conductivity: [&k, &k],
                reference_phase: SolidLiquidPhase::Solid,
                reference: &points[0],
                lower_specific_enthalpy_j_kg: datum,
                maximum_fusion_rounding_error_j_kg: fusion_rounding,
                coexistence: closure,
            },
        ),
    }
}

// Independent original Cp integrals: solid h=2*x+0.02*x^2, liquid
// h=350+4*y+0.02*y^2. Fusion occupies [150,350] at 350 K.
fn oracle(h: f64, closure: CoexistenceConductivity) -> (f64, f64, f64, f64) {
    if h < 150.0 {
        let x = h / (1.0 + (1.0 + 0.02 * h).sqrt());
        return (300.0 + x, 0.0, 2.0 + 0.04 * x, 0.04 / (2.0 + 0.04 * x));
    }
    if h >= 350.0 {
        let heat = h - 350.0;
        let y = heat / (2.0 + (4.0 + 0.02 * heat).sqrt());
        return (350.0 + y, 1.0, 8.0 + 0.04 * y, 0.04 / (4.0 + 0.04 * y));
    }
    let f = (h - 150.0) / 200.0;
    let (k, derivative) = match closure {
        CoexistenceConductivity::Parallel => (4.0 + 8.0 * f / (1.0 + f), 0.04 / (1.0 + f).powi(2)),
        CoexistenceConductivity::Series => (4.0 * (1.0 + f), 0.02),
    };
    (350.0, f, k, derivative)
}

#[test]
fn compact_sources_keep_original_cp_inverse_and_latent_conductivity() {
    let card = card_with_capacities(true, false, [[2.0, 4.0], [4.0, 8.0]]);
    for closure in [
        CoexistenceConductivity::Parallel,
        CoexistenceConductivity::Series,
    ] {
        let refined = resolve_capacity(&card, closure, 0.0, Some(0.05), 0.0).unwrap();
        assert!(refined.caloric().curve().knots().len() > 4);
        for datum in [0.0, 2.0_f64.powi(30)] {
            let material = resolve_capacity(&card, closure, datum, None, 0.0).unwrap();
            let curve = material.caloric().curve();
            assert_eq!(curve.knots().len(), 4);
            assert_eq!(material.reference_density_kg_m3(), 2.0);
            assert_eq!(material.card_identity(), card.content_hash());
            assert_eq!(
                material.caloric().plateau_enthalpy_bounds_j_kg(),
                [datum + 150.0, datum + 350.0]
            );
            for h in [37.5, 150.0, 250.0, 350.0, 500.0, 950.0] {
                let expected = oracle(h, closure);
                let state = curve.state_at_specific_enthalpy(datum + h).unwrap();
                let transport = material
                    .conductivity_at_specific_enthalpy(datum + h)
                    .unwrap();
                near(state.temperature_k(), expected.0, 2e-12);
                near(state.liquid_mass_fraction(), expected.1, 2e-14);
                near(transport.conductivity_w_m_k, expected.2, 2e-12);
                near(transport.derivative_per_j_kg, expected.3, 2e-13);
                near(
                    refined
                        .caloric()
                        .curve()
                        .temperature_at_specific_enthalpy(h)
                        .unwrap(),
                    state.temperature_k(),
                    2e-11,
                );
            }
        }
    }
}

fn independent_residual(h: &[f64], old: &[f64], closure: CoexistenceConductivity) -> [f64; 4] {
    let values: Vec<_> = h.iter().map(|&value| oracle(value, closure)).collect();
    let k = values.iter().map(|value| value.2 / 4.0).sum::<f64>();
    let t: Vec<_> = values.iter().map(|value| value.0).collect();
    let laplace = [
        3.0 * t[0] - t[1] - t[2] - t[3],
        t[1] - t[0],
        t[2] - t[0],
        t[3] - t[0],
    ];
    std::array::from_fn(|i| {
        MASS * (h[i] - old[i])
            + DT * (k * laplace[i] / 6.0 - (Q.iter().sum::<f64>() + Q[i]) / 120.0)
    })
}

#[test]
fn compact_nonlinear_cp_reaches_native_spatial_balance_and_jacobian() {
    let card = card_with_capacities(true, false, [[2.0, 4.0], [4.0, 8.0]]);
    let mesh = mesh();
    let boundary = ThermalBoundaryBuilder::new(&mesh)
        .adiabatic_remainder()
        .finish()
        .unwrap();
    let q = ScalarField::nodal("source-exact phase heat", 4, Q.to_vec()).unwrap();
    let old = [60.0, 210.0, 290.0, 500.0];
    for closure in [
        CoexistenceConductivity::Parallel,
        CoexistenceConductivity::Series,
    ] {
        let material = resolve_capacity(&card, closure, 0.0, None, 0.0).unwrap();
        with_cx(|cx| {
            let solver = material.prepare_uniform(cx, &mesh, budget()).unwrap();
            let response = solver
                .linearize_step(
                    cx,
                    material.problem(&mesh, &boundary, &q),
                    None,
                    &old,
                    DT,
                    config(),
                )
                .unwrap();
            let next = &response.primal().specific_enthalpy_j_kg;
            for value in independent_residual(next, &old, closure) {
                near(value, 0.0, 3e-10);
            }
            near(
                response.primal().stored_energy_change_j,
                DT * Q.iter().sum::<f64>() / 24.0,
                1e-10,
            );
            for column in 0..4 {
                let mut direction = [0.0; 4];
                direction[column] = 1.0;
                let action = response.apply_jacobian(cx, &direction).unwrap();
                let mut plus = next.clone();
                let mut minus = next.clone();
                plus[column] += 0.001;
                minus[column] -= 0.001;
                let rp = independent_residual(&plus, &old, closure);
                let rm = independent_residual(&minus, &old, closure);
                for row in 0..4 {
                    near(action[row], (rp[row] - rm[row]) / 0.002, 4e-9);
                }
            }
        });
    }
}

#[test]
fn compact_original_cp_heats_and_freezes_with_fixed_mass_in_both_native_integrators() {
    let card = card_with_capacities(true, false, [[2.0, 4.0], [4.0, 8.0]]);
    let material =
        resolve_capacity(&card, CoexistenceConductivity::Parallel, 0.0, None, 0.0).unwrap();
    let mesh = mesh();
    let boundary = ThermalBoundaryBuilder::new(&mesh)
        .adiabatic_remainder()
        .finish()
        .unwrap();
    with_cx(|cx| {
        let solver = material.prepare_uniform(cx, &mesh, budget()).unwrap();
        let masses = solver.reference_nodal_masses_kg().to_vec();
        for sdirk in [false, true] {
            let mut h = vec![60.0; 4];
            for power in [100.0, -100.0] {
                let q = ScalarField::uniform("held heating or cooling", power).unwrap();
                for _ in 0..8 {
                    let expected_h = h[0] + power / 2.0;
                    let problem = material.problem(&mesh, &boundary, &q);
                    let (next, temperature, stored) = if sdirk {
                        let input = EnthalpySdirk2Input {
                            problem,
                            interfaces: None,
                            radiation: EnthalpyRadiation::None,
                            contacts: None,
                        };
                        let step = solver
                            .advance_sdirk2(cx, [input; 2], &h, 1.0, config())
                            .unwrap();
                        (
                            step.second_stage.specific_enthalpy_j_kg,
                            step.second_stage.temperature,
                            step.stored_energy_change_j,
                        )
                    } else {
                        let step = solver
                            .advance(cx, problem, None, &h, 1.0, config())
                            .unwrap();
                        (
                            step.specific_enthalpy_j_kg,
                            step.temperature,
                            step.stored_energy_change_j,
                        )
                    };
                    for vertex in 0..4 {
                        near(next[vertex], expected_h, 2e-8);
                        near(
                            temperature[vertex],
                            oracle(expected_h, material.coexistence()).0,
                            2e-8,
                        );
                    }
                    near(stored, power / 6.0, 1e-9);
                    h = next;
                }
            }
            for value in h {
                near(value, 60.0, 3e-8);
            }
        }
        assert_eq!(solver.reference_nodal_masses_kg(), masses);
    });
}

#[test]
fn compact_ingress_retains_source_refusals_and_real_conductivity_corners() {
    let card = card_with_capacities(true, true, [[2.0, 4.0], [4.0, 8.0]]);
    for rounding in [-1.0, f64::NAN, f64::INFINITY] {
        assert!(
            resolve_capacity(
                &card,
                CoexistenceConductivity::Parallel,
                0.0,
                None,
                rounding
            )
            .is_err()
        );
    }
    assert!(
        resolve_capacity(
            &card_with_capacities(false, false, [[2.0, 4.0], [4.0, 8.0]]),
            CoexistenceConductivity::Parallel,
            0.0,
            None,
            0.0,
        )
        .is_err()
    );
    let material =
        resolve_capacity(&card, CoexistenceConductivity::Parallel, 0.0, None, 0.0).unwrap();
    assert_eq!(material.caloric().curve().knots().len(), 4);
    near(
        material.conductivity_sources()[0].eval(325.0).unwrap(),
        3.5,
        0.0,
    );
    let mesh = mesh();
    let boundary = ThermalBoundaryBuilder::new(&mesh)
        .adiabatic_remainder()
        .finish()
        .unwrap();
    let q = ScalarField::uniform("equilibrium", 0.0).unwrap();
    with_cx(|cx| {
        let solver = material.prepare_uniform(cx, &mesh, budget()).unwrap();
        let result = solver.linearize_step(
            cx,
            material.problem(&mesh, &boundary, &q),
            None,
            &[62.5; 4],
            DT,
            config(),
        );
        assert!(matches!(result, Err(EnthalpyAdjointError::MaterialKink)));
    });
}
