//! Native P1 entropy against independent caloric and face-matrix identities.
#[path = "entropy/robin.rs"]
mod robin;
use super::*;
use fs_conduction::transient::backward_euler::{ArealHeatCapacity, SurfaceHeatCapacity};
use fs_conduction::transient::enthalpy::entropy::{EnthalpyEntropyBalance, EnthalpyEntropyError};

fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        actual.is_finite() && expected.is_finite() && (actual - expected).abs() <= tolerance,
        "{actual:.17e} != {expected:.17e} (allowance {tolerance:e})"
    );
}

fn sensible_h(temperature: &[f64]) -> Vec<f64> {
    temperature.iter().map(|t| 5.0 * (t - 280.0)).collect()
}

#[test]
fn closed_sensible_step_matches_exact_star_graph_and_original_trajectory() {
    let mesh = tetrahedron();
    let boundary = adiabatic(&mesh);
    let material = ConductivityModel::isotropic_declared(2.0).unwrap();
    let source = ScalarField::Uniform(0.0);
    let chart = sensible_curve();
    let old_t = [400.0, 300.0, 300.0, 300.0];
    let old_h = sensible_h(&old_t);
    let saved = bits(&old_h);
    with_cx(|cx| {
        let model =
            EnthalpyBackwardEuler::uniform(cx, &mesh, &chart, REFERENCE_DENSITY, budget()).unwrap();
        let p = problem(&mesh, &boundary, &material, &source);
        let original = model.advance(cx, p, None, &old_h, 0.25, config()).unwrap();
        let (step, entropy) = model
            .advance_with_entropy(cx, p, None, &old_h, 0.25, config(), 1e-10)
            .unwrap();
        assert_eq!(
            bits(&step.specific_enthalpy_j_kg),
            bits(&original.specific_enthalpy_j_kg)
        );
        assert_eq!(step.newton, original.newton);
        assert_eq!(bits(&old_h), saved);
        let capacity = REFERENCE_DENSITY * VOLUME * 5.0 / 4.0;
        let conductance = 2.0 * VOLUME;
        let gap = 100.0 / (1.0 + 4.0 * conductance * 0.25 / capacity);
        let expected = [
            325.0 + 0.75 * gap,
            325.0 - 0.25 * gap,
            325.0 - 0.25 * gap,
            325.0 - 0.25 * gap,
        ];
        let storage = capacity
            * expected
                .iter()
                .zip(old_t)
                .map(|(t, old)| (t / old).ln())
                .sum::<f64>();
        let end_energy_over_t = capacity
            * expected
                .iter()
                .zip(old_t)
                .map(|(t, old)| (t - old) / t)
                .sum::<f64>();
        let generation = 3.0 * 0.25 * conductance * gap * gap / (expected[0] * expected[1]);
        for (t, reference) in step.temperature.iter().zip(expected) {
            near(*t, reference, 1e-8);
        }
        near(entropy.thermal_storage_change_j_k, storage, 1e-10);
        near(entropy.conduction_generation_j_k, generation, 1e-10);
        near(
            entropy.temporal_defect_j_k,
            storage - end_energy_over_t,
            1e-10,
        );
        near(entropy.algebraic_residual_j_k, 0.0, 1e-10);
        near(entropy.balance_residual_j_k, 0.0, 1e-12);
        assert!(entropy.conduction_generation_j_k > 0.0);
        assert!(entropy.temporal_defect_j_k > 0.0);
        assert!(!entropy.has_unclassified_inputs);
        assert_eq!(entropy.dirichlet_entropy_input_j_k, 0.0);
        assert_eq!(entropy.robin_entropy_input_j_k, 0.0);
        near(step.stored_energy_change_j, 0.0, 1e-9);
    });
}

#[test]
fn film_entropy_uses_sensible_capacity_without_counting_host_latent_heat_twice() {
    let mesh = tetrahedron();
    let boundary = adiabatic(&mesh);
    let material = ConductivityModel::isotropic_declared(2.0).unwrap();
    let chart = latent_curve();
    let zero = ScalarField::Uniform(0.0);
    let source = ScalarField::Uniform(100.0);
    let coating = SurfaceHeatCapacity::uniform(0, ArealHeatCapacity::declared(2.0).unwrap());
    let film_capacity = 2.0 * mesh.boundary()[0].area;
    let mass = REFERENCE_DENSITY * VOLUME;
    with_cx(|cx| {
        let model = EnthalpyBackwardEuler::uniform(cx, &mesh, &chart, REFERENCE_DENSITY, budget())
            .unwrap()
            .with_surface_heat_capacity(cx, &[coating])
            .unwrap();
        // Host Cp=2 on this sensible segment; film stores Cfilm*dT.
        let (endpoint, account) = model
            .evaluate_endpoint_with_entropy(
                cx,
                problem(&mesh, &boundary, &material, &zero),
                None,
                &[40.0; 4],
                &[80.0; 4],
                0.1,
                config(),
            )
            .unwrap();
        let ds = (mass * 2.0 + film_capacity) * (340.0_f64 / 320.0).ln();
        near(account.thermal_storage_change_j_k, ds, 1e-12);
        near(
            endpoint.stored_energy_change_j,
            mass * 40.0 + film_capacity * 20.0,
            1e-10,
        );
        near(account.balance_residual_j_k, 0.0, 1e-12);
        assert!(account.algebraic_residual_j_k > 0.0);
        // A native step entirely on T=350 K stores host latent heat only.
        let (step, account) = model
            .advance_with_entropy(
                cx,
                problem(&mesh, &boundary, &material, &source),
                None,
                &[150.0; 4],
                0.1,
                config(),
                1e-10,
            )
            .unwrap();
        for h in &step.specific_enthalpy_j_kg {
            near(*h, 155.0, 1e-8);
        }
        for t in &step.temperature {
            assert_eq!(*t, 350.0);
        }
        near(step.stored_energy_change_j, mass * 5.0, 1e-9);
        near(
            account.thermal_storage_change_j_k,
            mass * 5.0 / 350.0,
            1e-11,
        );
        near(
            account.source_at_node_j_k,
            100.0 * VOLUME * 0.1 / 350.0,
            1e-12,
        );
        near(account.temporal_defect_j_k, 0.0, 1e-12);
        near(account.conduction_generation_j_k, 0.0, 1e-12);
        assert!(account.has_unclassified_inputs);
    });
}

#[test]
fn robin_faces_use_their_actual_baths_and_the_original_consistent_flux() {
    let mesh = tetrahedron();
    let first = mesh.boundary()[0].vertices;
    let second = mesh.boundary()[1].vertices;
    let boundary = ThermalBoundaryBuilder::new(&mesh)
        .region(
            "cold bath",
            |face| face.vertices == first,
            ThermalBc::robin(2.0, 300.0).unwrap(),
        )
        .unwrap()
        .region(
            "hot bath",
            |face| face.vertices == second,
            ThermalBc::robin(3.0, 450.0).unwrap(),
        )
        .unwrap()
        .adiabatic_remainder()
        .finish()
        .unwrap();
    let material = ConductivityModel::isotropic_declared(2.0).unwrap();
    let zero = ScalarField::Uniform(0.0);
    let chart = sensible_curve();
    let t = [320.0, 350.0, 380.0, 400.0];
    let h = sensible_h(&t);
    let dt = 0.3;
    let mut bath = 0.0;
    let mut generation = 0.0;
    for (slot, htc, reference) in [(0, 2.0, 300.0), (1, 3.0, 450.0)] {
        let face = &mesh.boundary()[slot];
        for (a, &va) in face.vertices.iter().enumerate() {
            let mut outward = 0.0;
            for (b, &vb) in face.vertices.iter().enumerate() {
                let weight = htc * face.area / 12.0 * if a == b { 2.0 } else { 1.0 };
                outward += weight * (t[vb as usize] - reference);
            }
            bath -= dt * outward / reference;
            generation += dt * outward * (1.0 / reference - 1.0 / t[va as usize]);
        }
    }
    with_cx(|cx| {
        let model =
            EnthalpyBackwardEuler::uniform(cx, &mesh, &chart, REFERENCE_DENSITY, budget()).unwrap();
        let (_, account) = model
            .evaluate_endpoint_with_entropy(
                cx,
                problem(&mesh, &boundary, &material, &zero),
                None,
                &h,
                &h,
                dt,
                config(),
            )
            .unwrap();
        near(account.robin_entropy_input_j_k, bath, 1e-12);
        near(account.robin_generation_j_k, generation, 1e-12);
        let bulk = dt
            * 2.0
            * VOLUME
            * t[1..]
                .iter()
                .map(|tj| (t[0] - tj).powi(2) / (t[0] * tj))
                .sum::<f64>();
        near(account.conduction_generation_j_k, bulk, 1e-12);
        near(account.thermal_storage_change_j_k, 0.0, 0.0);
        near(account.temporal_defect_j_k, 0.0, 0.0);
        near(account.balance_residual_j_k, 0.0, 1e-12);
        assert!(!account.has_unclassified_inputs);
    });
}

#[test]
fn positive_anisotropic_tensor_can_refuse_a_negative_discrete_entropy_step() {
    let mesh = tetrahedron();
    let boundary = adiabatic(&mesh);
    // SPD eigenvalues are 1-c, 1, 1+c, but the FE edge 1--2 has negative g.
    let c = 0.9999;
    let material =
        ConductivityModel::constant_tensor([[1.0, c, 0.0], [c, 1.0, 0.0], [0.0, 0.0, 1.0]])
            .unwrap();
    let zero = ScalarField::Uniform(0.0);
    let chart = sensible_curve();
    let middle = (300.0_f64 * 480.0).sqrt();
    let t = [middle, 300.0, 480.0, middle];
    let h = sensible_h(&t);
    let saved = bits(&h);
    let rate = (1.0 + c)
        * VOLUME
        * ((middle - 300.0).powi(2) / (middle * 300.0)
            + (middle - 480.0).powi(2) / (middle * 480.0))
        - c * VOLUME * 180.0_f64.powi(2) / (300.0 * 480.0);
    assert!(rate < 0.0);
    with_cx(|cx| {
        let model =
            EnthalpyBackwardEuler::uniform(cx, &mesh, &chart, REFERENCE_DENSITY, budget()).unwrap();
        let p = problem(&mesh, &boundary, &material, &zero);
        let (_, account) = model
            .evaluate_endpoint_with_entropy(cx, p, None, &h, &h, 0.001, config())
            .unwrap();
        near(account.conduction_generation_j_k, 0.001 * rate, 1e-12);
        near(account.balance_residual_j_k, 0.0, 1e-12);
        assert!(matches!(
            account.admit(1e-10),
            Err(EnthalpyEntropyError::NegativeProduction {
                account: "bulk/surface/contact",
                ..
            })
        ));
        // The unmodified native energy solve can accept this state; the
        // explicitly requested second-law gate must not publish it.
        let ordinary = model.advance(cx, p, None, &h, 0.001, config()).unwrap();
        near(ordinary.energy_residual_j, 0.0, 1e-9);
        assert!(matches!(
            model.advance_with_entropy(cx, p, None, &h, 0.001, config(), 1e-10,),
            Err(EnthalpyEntropyError::NegativeProduction {
                account: "bulk/surface/contact",
                ..
            })
        ));
        assert_eq!(bits(&h), saved);
    });
}

#[test]
fn consistent_robin_negative_production_is_retained_without_lumping_the_face() {
    let mesh = tetrahedron();
    let face = &mesh.boundary()[0];
    let selected = face.vertices;
    let htc = 12.0 / face.area; // The consistent 3x3 matrix is I + 11^T W/K.
    let boundary = ThermalBoundaryBuilder::new(&mesh)
        .region(
            "actual 1000 K bath",
            |candidate| candidate.vertices == selected,
            ThermalBc::robin(htc, 1000.0).unwrap(),
        )
        .unwrap()
        .adiabatic_remainder()
        .finish()
        .unwrap();
    let material = ConductivityModel::isotropic_declared(0.001).unwrap();
    let zero = ScalarField::Uniform(0.0);
    let chart = EquilibriumEnthalpyPhaseCurve::try_single_phase(
        ContentHash([0x73; 32]),
        SolidLiquidPhase::Solid,
        [(0.0, 100.0), (34500.0, 7000.0)]
            .into_iter()
            .map(|(h, t)| EnthalpyPhaseKnot {
                specific_enthalpy_j_kg: h,
                temperature_k: t,
                liquid_mass_fraction: 0.0,
                bulk_density_kg_m3: REFERENCE_DENSITY,
            })
            .collect(),
    )
    .unwrap();
    let mut t = [300.0; 4];
    t[selected[0] as usize] = 6000.0;
    let h: Vec<_> = t.iter().map(|t| 5.0 * (t - 100.0)).collect();
    with_cx(|cx| {
        let model =
            EnthalpyBackwardEuler::uniform(cx, &mesh, &chart, REFERENCE_DENSITY, budget()).unwrap();
        let (_, account) = model
            .evaluate_endpoint_with_entropy(
                cx,
                problem(&mesh, &boundary, &material, &zero),
                None,
                &h,
                &h,
                0.3,
                config(),
            )
            .unwrap();
        // q=[8600,2900,2900] W at T=[6000,300,300] K. Positive total
        // outward heat does not imply nonnegative nodal consistent-face entropy.
        near(account.robin_entropy_input_j_k, -4.32, 1e-11);
        near(account.robin_generation_j_k, -1.91, 1e-11);
        near(account.balance_residual_j_k, 0.0, 1e-11);
        assert!(account.conduction_generation_j_k >= -1e-12);
        assert!(matches!(
            account.admit(1e-10),
            Err(EnthalpyEntropyError::NegativeProduction {
                account: "Robin",
                ..
            })
        ));
    });
}

#[test]
fn nonpositive_robin_baths_nonfinite_accounts_and_cancellation_refuse() {
    let mesh = tetrahedron();
    let boundary = ThermalBoundaryBuilder::new(&mesh)
        .region(
            "nonpositive bath has no thermal entropy",
            |_| true,
            ThermalBc::Robin {
                htc: ScalarField::Uniform(1.0),
                t_ref: ScalarField::Nodal(vec![0.0, 310.0, 320.0, 330.0]),
            },
        )
        .unwrap()
        .finish()
        .unwrap();
    let material = ConductivityModel::isotropic_declared(2.0).unwrap();
    let zero = ScalarField::Uniform(0.0);
    let chart = sensible_curve();
    with_cx(|cx| {
        let model =
            EnthalpyBackwardEuler::uniform(cx, &mesh, &chart, REFERENCE_DENSITY, budget()).unwrap();
        let p = problem(&mesh, &boundary, &material, &zero);
        assert!(matches!(
            model.evaluate_endpoint_with_entropy(
                cx,
                p,
                None,
                &[100.0; 4],
                &[100.0; 4],
                1.0,
                config(),
            ),
            Err(EnthalpyError::InvalidInput(
                "entropy requires positive actual Robin bath temperatures"
            ))
        ));
        assert!(matches!(
            model.advance_with_entropy(cx, p, None, &[100.0; 4], 1.0, config(), -1.0,),
            Err(EnthalpyEntropyError::Physical(EnthalpyError::InvalidInput(
                _
            )))
        ));
        with_cancelled_cx(|cancelled| {
            assert!(matches!(
                model.advance_with_entropy(cancelled, p, None, &[100.0; 4], 1.0, config(), 1e-9,),
                Err(EnthalpyEntropyError::Physical(EnthalpyError::Conduction(
                    ConductionError::Cancelled { .. }
                )))
            ));
        });
    });
    for account in [
        EnthalpyEntropyBalance {
            temporal_defect_j_k: f64::NAN,
            ..Default::default()
        },
        EnthalpyEntropyBalance {
            algebraic_residual_j_k: f64::INFINITY,
            ..Default::default()
        },
    ] {
        assert!(account.admit(1e-9).is_err());
    }
}

#[test]
fn nodal_sources_and_fluxes_retain_consistent_moments_and_unclassified_inputs() {
    let mesh = tetrahedron();
    let face = &mesh.boundary()[0];
    let selected = face.vertices;
    let flux = [12.0, -8.0, 7.0, -3.0];
    let boundary = ThermalBoundaryBuilder::new(&mesh)
        .region(
            "prescribed heat only",
            |candidate| candidate.vertices == selected,
            ThermalBc::Neumann {
                outward_flux: ScalarField::Nodal(flux.to_vec()),
            },
        )
        .unwrap()
        .adiabatic_remainder()
        .finish()
        .unwrap();
    let material = ConductivityModel::isotropic_declared(2.0).unwrap();
    let source_values = [100.0, -40.0, 30.0, -90.0]; // Zero integral, nonzero local source entropy.
    let source = ScalarField::Nodal(source_values.to_vec());
    let chart = sensible_curve();
    let t = [310.0, 340.0, 380.0, 450.0];
    let h = sensible_h(&t);
    let dt = 0.4;
    let mut source_entropy = 0.0;
    for i in 0..4 {
        let mut load = 0.0;
        for (j, &q) in source_values.iter().enumerate() {
            load += VOLUME / 20.0 * if i == j { 2.0 * q } else { q };
        }
        source_entropy += dt * load / t[i];
    }
    let mut flux_entropy = 0.0;
    for (a, &va) in selected.iter().enumerate() {
        let mut out = 0.0;
        for (b, &vb) in selected.iter().enumerate() {
            out += face.area / 12.0
                * if a == b {
                    2.0 * flux[vb as usize]
                } else {
                    flux[vb as usize]
                };
        }
        flux_entropy -= dt * out / t[va as usize];
    }
    with_cx(|cx| {
        let model =
            EnthalpyBackwardEuler::uniform(cx, &mesh, &chart, REFERENCE_DENSITY, budget()).unwrap();
        let (_, account) = model
            .evaluate_endpoint_with_entropy(
                cx,
                problem(&mesh, &boundary, &material, &source),
                None,
                &h,
                &h,
                dt,
                config(),
            )
            .unwrap();
        near(account.source_at_node_j_k, source_entropy, 1e-12);
        near(account.neumann_at_node_j_k, flux_entropy, 1e-12);
        near(account.balance_residual_j_k, 0.0, 1e-12);
        assert!(account.has_unclassified_inputs);
        assert_eq!(account.dirichlet_entropy_input_j_k, 0.0);
        assert_eq!(account.robin_entropy_input_j_k, 0.0);
    });
}

#[test]
fn nonlinear_film_entropy_keeps_original_piecewise_integrals_and_nodal_inventory() {
    use fs_conduction::transient::enthalpy::SurfaceHeatStorage;
    use fs_material::thermal::transport::storage::TabulatedHeatStorage;

    let mesh = tetrahedron();
    let boundary = ThermalBoundaryBuilder::new(&mesh)
        .region("400 K bath", |_| true, ThermalBc::dirichlet(400.0).unwrap())
        .unwrap()
        .finish()
        .unwrap();
    let material = ConductivityModel::isotropic_declared(2.0).unwrap();
    let source = ScalarField::Uniform(0.0);
    let chart = sensible_curve();
    // Total face inventory: C=.01*T below 350 K; C=.02*T-3.5 above it.
    let nonlinear =
        TabulatedHeatStorage::try_new(&[[280.0, 2.8], [350.0, 3.5], [450.0, 5.5]]).unwrap();
    let constant = SurfaceHeatCapacity::uniform(0, ArealHeatCapacity::declared(2.0).unwrap());
    let constant_capacity = 2.0 * mesh.boundary()[0].area;
    let old_h = sensible_h(&[300.0; 4]);
    let old_bits = bits(&old_h);
    with_cx(|cx| {
        let model = EnthalpyBackwardEuler::uniform(cx, &mesh, &chart, REFERENCE_DENSITY, budget())
            .unwrap()
            .with_surface_heat_capacity(cx, &[constant])
            .unwrap()
            .with_surface_heat_storage(cx, &[SurfaceHeatStorage::tabulated(0, &nonlinear)])
            .unwrap();
        assert!(std::ptr::eq(
            model.surface_heat_storage()[0].tabulated_storage().unwrap(),
            &nonlinear,
        ));
        let p = problem(&mesh, &boundary, &material, &source);
        let (step, entropy) = model
            .advance_with_entropy(cx, p, None, &old_h, 0.5, config(), 1e-10)
            .unwrap();
        // Integrating each side of the original 350 K corner gives these
        // independent exact totals, before distributing one third per face node.
        let film_energy = 162.5 + 200.0;
        let film_entropy = 1.5 - 3.5 * (400.0_f64 / 350.0).ln();
        let capacity = REFERENCE_DENSITY * VOLUME * 5.0 + constant_capacity;
        let energy = capacity * 100.0 + film_energy;
        let storage_entropy = capacity * (400.0_f64 / 300.0).ln() + film_entropy;
        for &temperature in &step.temperature {
            near(temperature, 400.0, 1e-12);
        }
        near(step.stored_energy_change_j, energy, 1e-9);
        near(step.dirichlet_in_w * 0.5, energy, 1e-9);
        near(entropy.thermal_storage_change_j_k, storage_entropy, 1e-11);
        near(entropy.dirichlet_entropy_input_j_k, energy / 400.0, 1e-11);
        near(
            entropy.temporal_defect_j_k,
            storage_entropy - energy / 400.0,
            1e-11,
        );
        near(entropy.conduction_generation_j_k, 0.0, 1e-12);
        near(entropy.algebraic_residual_j_k, 0.0, 0.0);
        near(entropy.balance_residual_j_k, 0.0, 1e-11);
        assert!(!entropy.has_unclassified_inputs);
        assert_eq!(bits(&old_h), old_bits);

        // Host support extends to 500 K; the original film stops at 450 K.
        let outside = ThermalBoundaryBuilder::new(&mesh)
            .region("475 K bath", |_| true, ThermalBc::dirichlet(475.0).unwrap())
            .unwrap()
            .finish()
            .unwrap();
        assert!(
            model
                .advance_with_entropy(
                    cx,
                    problem(&mesh, &outside, &material, &source),
                    None,
                    &old_h,
                    0.5,
                    config(),
                    1e-10,
                )
                .is_err()
        );
        assert_eq!(bits(&old_h), old_bits);
    });
}

#[test]
fn fully_prescribed_entropy_counts_boundary_storage_and_real_sources_once() {
    let mesh = tetrahedron();
    let material = ConductivityModel::isotropic_declared(2.0).unwrap();
    let source = ScalarField::Uniform(12.0);
    let chart = sensible_curve();
    let capacity = REFERENCE_DENSITY * VOLUME * 5.0;
    let dt = 0.5;
    let supplied = 12.0 * VOLUME * dt;
    with_cx(|cx| {
        let model =
            EnthalpyBackwardEuler::uniform(cx, &mesh, &chart, REFERENCE_DENSITY, budget()).unwrap();
        for (old_t, new_t) in [(300.0_f64, 400.0_f64), (400.0_f64, 300.0_f64)] {
            let boundary = ThermalBoundaryBuilder::new(&mesh)
                .region(
                    "prescribed bath",
                    |_| true,
                    ThermalBc::dirichlet(new_t).unwrap(),
                )
                .unwrap()
                .finish()
                .unwrap();
            let old_h = sensible_h(&[old_t; 4]);
            let saved = bits(&old_h);
            let p = problem(&mesh, &boundary, &material, &source);
            let (step, entropy) = model
                .advance_with_entropy(cx, p, None, &old_h, dt, config(), 1e-10)
                .unwrap();
            let energy = capacity * (new_t - old_t);
            let storage = capacity * (new_t / old_t).ln();
            assert_eq!(step.newton.iterations, 0);
            near(step.stored_energy_change_j, energy, 1e-10);
            near(step.dirichlet_in_w * dt, energy - supplied, 1e-10);
            near(entropy.thermal_storage_change_j_k, storage, 1e-12);
            near(
                entropy.dirichlet_entropy_input_j_k,
                (energy - supplied) / new_t,
                1e-12,
            );
            near(entropy.source_at_node_j_k, supplied / new_t, 1e-12);
            near(entropy.temporal_defect_j_k, storage - energy / new_t, 1e-12);
            near(entropy.conduction_generation_j_k, 0.0, 1e-12);
            near(entropy.algebraic_residual_j_k, 0.0, 0.0);
            near(entropy.balance_residual_j_k, 0.0, 1e-12);
            assert!(entropy.temporal_defect_j_k > 0.0);
            assert!(entropy.has_unclassified_inputs);
            let (endpoint, measured) = model
                .evaluate_endpoint_with_entropy(
                    cx,
                    p,
                    None,
                    &old_h,
                    &step.specific_enthalpy_j_kg,
                    dt,
                    config(),
                )
                .unwrap();
            assert_eq!(endpoint.free_dofs, 0);
            assert_eq!(endpoint.residual_norm_j, 0.0);
            assert_eq!(measured, entropy);
            assert_eq!(bits(&old_h), saved);
        }
    });
}

#[test]
fn heterogeneous_entropy_uses_each_original_chart_and_frozen_reference_mass() {
    use fs_conduction::transient::enthalpy::heterogeneous::{
        HeterogeneousEnthalpyBackwardEuler, ReferenceEnthalpyMaterial,
    };

    let mesh = ConductionMesh::new(
        TetComplex::from_tets(8, vec![[0, 1, 2, 3], [4, 5, 6, 7]]),
        vec![
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
            [2.0, 0.0, 0.0],
            [3.0, 0.0, 0.0],
            [2.0, 1.0, 0.0],
            [2.0, 0.0, 1.0],
        ],
    )
    .unwrap();
    let first = sensible_curve();
    let second = EquilibriumEnthalpyPhaseCurve::try_single_phase(
        ContentHash([0x74; 32]),
        SolidLiquidPhase::Solid,
        [(0.0, 300.0), (800.0, 500.0)]
            .into_iter()
            .map(|(h, t)| EnthalpyPhaseKnot {
                specific_enthalpy_j_kg: h,
                temperature_k: t,
                liquid_mass_fraction: 0.0,
                // Deliberately not the separately declared frozen density.
                bulk_density_kg_m3: 9.0,
            })
            .collect(),
    )
    .unwrap();
    let materials = [
        ReferenceEnthalpyMaterial {
            curve: &first,
            reference_density_kg_m3: 2.0,
        },
        ReferenceEnthalpyMaterial {
            curve: &second,
            reference_density_kg_m3: 6.0,
        },
    ];
    let boundary = adiabatic(&mesh);
    let conductivity = ConductivityModel::isotropic_declared(2.0).unwrap();
    let source = ScalarField::Uniform(12.0);
    let old_h = [150.0, 150.0, 150.0, 150.0, 240.0, 240.0, 240.0, 240.0];
    let saved = bits(&old_h);
    with_cx(|cx| {
        let model = HeterogeneousEnthalpyBackwardEuler::new(
            cx,
            &mesh,
            &materials,
            &[0, 1],
            EnthalpyBudget {
                max_vertices: 8,
                max_elements: 2,
            },
        )
        .unwrap();
        assert!(std::ptr::eq(model.phase_curve_at(0).unwrap(), &first));
        assert!(std::ptr::eq(model.phase_curve_at(7).unwrap(), &second));
        let p = problem(&mesh, &boundary, &conductivity, &source);
        let (step, entropy) = model
            .advance_with_entropy(cx, p, None, &old_h, 0.5, config(), 1e-10)
            .unwrap();
        let (first_t, second_t) = (310.6_f64, 360.25_f64);
        for vertex in 0..8 {
            let (expected_h, expected_t, expected_mass) = if vertex < 4 {
                (153.0, first_t, 1.0 / 12.0)
            } else {
                (241.0, second_t, 0.25)
            };
            near(step.specific_enthalpy_j_kg[vertex], expected_h, 1e-8);
            near(step.temperature[vertex], expected_t, 1e-9);
            near(
                model.reference_nodal_masses_kg()[vertex],
                expected_mass,
                1e-15,
            );
        }
        // Each disconnected tet receives exactly one joule. Its caloric
        // entropy uses its own Cp and reference mass, never chart bulk density.
        let storage = (2.0 * VOLUME * 5.0) * (first_t / 310.0).ln()
            + (6.0 * VOLUME * 4.0) * (second_t / 360.0).ln();
        let source_entropy = 1.0 / first_t + 1.0 / second_t;
        near(step.stored_energy_change_j, 2.0, 1e-9);
        near(entropy.thermal_storage_change_j_k, storage, 1e-11);
        near(entropy.source_at_node_j_k, source_entropy, 1e-11);
        near(entropy.temporal_defect_j_k, storage - source_entropy, 1e-11);
        near(entropy.conduction_generation_j_k, 0.0, 1e-11);
        near(entropy.algebraic_residual_j_k, 0.0, 1e-10);
        near(entropy.balance_residual_j_k, 0.0, 1e-11);
        assert!(entropy.has_unclassified_inputs);
        let (endpoint, measured) = model
            .evaluate_endpoint_with_entropy(
                cx,
                p,
                None,
                &old_h,
                &step.specific_enthalpy_j_kg,
                0.5,
                config(),
            )
            .unwrap();
        near(endpoint.stored_energy_change_j, 2.0, 1e-9);
        assert_eq!(measured, entropy);
        assert_eq!(bits(&old_h), saved);
    });
}
