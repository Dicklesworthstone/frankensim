//! Source surfaces and prescribed baths share the original caloric history.
use super::*;
use fs_conduction::ThermalBc;

fn surface_card() -> InterfaceSystemCard {
    let dims = fs_conduction::AREA_SPECIFIC_THERMAL_RESISTANCE_DIMS;
    let mut claims = ClaimSet::new();
    claims
        .insert_claim(PropertyClaim {
            key: PropertyKey::new(
                fs_conduction::AREA_SPECIFIC_THERMAL_RESISTANCE_PROPERTY,
                dims,
            ),
            value: PropertyValue::Surface2D {
                axes: ["T".into(), "normal_pressure".into()],
                axis_dims: [
                    fs_conduction::TEMPERATURE_DIMS,
                    fs_qty::Dims([-1, 1, -2, 0, 0, 0]),
                ],
                samples: [
                    [280.0, 0.5e6],
                    [380.0, 0.5e6],
                    [280.0, 1.5e6],
                    [380.0, 1.5e6],
                ]
                .into_iter()
                .map(|p| (p, 0.1 + 0.006 * (p[0] - 280.0) + 1e-8 * (p[1] - 1e6)))
                .collect(),
                triangles: vec![[0, 1, 2], [1, 3, 2]],
                dims,
            },
            validity: ValidityDomain::unconstrained()
                .with("T", 280.0, 380.0)
                .with("normal_pressure", 0.5e6, 1.5e6),
            uncertainty: UncertaintyModel::Unstated,
            interpolation: InterpolationPolicy::LinearInside,
            observations: vec![],
            provenance: Provenance {
                source: "manufactured R(T,p) surface, not measurements".into(),
                license: "internal-test-use".into(),
                artifact: None,
            },
        })
        .unwrap();
    let side = |name: &str| SurfaceSpec {
        material: MaterialStateId {
            chemistry: name.into(),
            phase: "solid".into(),
            process: "fixture".into(),
            revision: 0,
        },
        texture_frame: "fixed reference normal".into(),
    };
    InterfaceSystemCard::assemble(
        side("A"),
        side("B"),
        SystemContext {
            medium: "dry".into(),
            third_body: None,
            environment: "vacuum".into(),
            history: "unaged".into(),
        },
        claims,
        vec![],
    )
    .unwrap()
}

#[test]
fn g1_source_surface_fixed_point_retains_latent_storage_and_independent_roots() {
    context(&CancelGate::new_clock_free(), |cx| {
        let mesh = mesh();
        let bc = ThermalBoundaryBuilder::new(&mesh)
            .adiabatic_remainder()
            .finish()
            .unwrap();
        let card = surface_card();
        let contacts =
            TemperatureContacts::new(cx, &mesh, &bc, vec![], vec![law(&mesh, &card)], 1).unwrap();
        let cold = chart(false);
        let hot = chart(true);
        let engine = HeterogeneousEnthalpyBackwardEuler::new(
            cx,
            &mesh,
            &[
                ReferenceEnthalpyMaterial {
                    curve: &cold,
                    reference_density_kg_m3: 1.0,
                },
                ReferenceEnthalpyMaterial {
                    curve: &hot,
                    reference_density_kg_m3: 2.0,
                },
            ],
            &[0, 1],
            EnthalpyBudget {
                max_vertices: 8,
                max_elements: 2,
            },
        )
        .unwrap();
        let material = ConductivityModel::isotropic(
            ConductivityTable::declared_curve(vec![(280.0, 1.6), (380.0, 3.6)]).unwrap(),
        );
        let source = ScalarField::Nodal(vec![0.0, 0.0, 0.0, 0.0, 60.0, 60.0, 60.0, 60.0]);
        let problem = ConductionProblem {
            mesh: &mesh,
            boundary: &bc,
            material: &material,
            element_materials: None,
            source: &source,
        };
        let initial = [80.0, 80.0, 80.0, 80.0, 450.0, 450.0, 450.0, 450.0];
        let mut old = initial.to_vec();
        for (index, dt) in [0.15, 0.35, 0.5].into_iter().enumerate() {
            let before = old.clone();
            let accepted = contacts
                .advance_enthalpy(cx, &engine, problem, &old, dt, config(), coupling())
                .unwrap();
            assert_eq!(old, before);
            assert!(accepted.iterations > 1);
            assert!(accepted.nonlinear_iterations > accepted.conduction.newton.iterations);
            assert!(
                accepted.krylov_iterations
                    >= accepted
                        .conduction
                        .newton
                        .history
                        .iter()
                        .map(|row| row.linear_iterations)
                        .sum()
            );
            assert!(accepted.physical.residual_norm_j <= coupling().residual_tolerance_j);
            close(accepted.physical.stored_energy_change_j, 10.0 * dt, 1e-8);
            for (&actual, &expected) in accepted
                .conduction
                .specific_enthalpy_j_kg
                .iter()
                .zip(&EXPECTED_H[index])
            {
                close(actual, expected, 3e-6);
            }
            let state = &accepted.physical_contacts.states()[0];
            close(
                state.resistance.value_m2_k_per_w(),
                0.1 + 0.006 * (state.interface_temperature_k - 280.0),
                2e-14,
            );
            assert_eq!(state.resistance.card_identity(), card.content_hash());
            old = accepted.conduction.specific_enthalpy_j_kg;
        }
        let stored: f64 = old
            .iter()
            .zip(initial)
            .zip(engine.reference_nodal_masses_kg())
            .map(|((&h, h0), &m)| m * (h - h0))
            .sum();
        close(stored, 10.0, 3e-8);
    });
}

fn linear_chart() -> EquilibriumEnthalpyPhaseCurve {
    EquilibriumEnthalpyPhaseCurve::try_new(
        ContentHash([91; 32]),
        [
            (0.0, 280.0, 0.0),
            (1000.0, 380.0, 0.0),
            (1100.0, 380.0, 1.0),
        ]
        .into_iter()
        .map(|(h, t, fraction)| EnthalpyPhaseKnot {
            specific_enthalpy_j_kg: h,
            temperature_k: t,
            liquid_mass_fraction: fraction,
            bulk_density_kg_m3: 1.0,
        })
        .collect(),
    )
    .unwrap()
}

fn bath(mesh: &ConductionMesh, temperature: f64) -> fs_conduction::ThermalBoundary {
    ThermalBoundaryBuilder::new(mesh)
        .region(
            "bath",
            |face| {
                let mut vertices = face.vertices;
                vertices.sort_unstable();
                vertices == [0, 1, 3]
            },
            ThermalBc::dirichlet(temperature).unwrap(),
        )
        .unwrap()
        .adiabatic_remainder()
        .finish()
        .unwrap()
}

// Independent 8x8 temperature system: k V grad(Na).grad(Nb), rho V Cp/4,
// and the consistent triangular contact blocks A(1+delta_ab)/(12 R).
// The scalar R(Tmean) root is bracketed separately from the product's h-Newton
// and relaxed contact iterations. No native assembled operator enters this oracle.
fn oracle(old: [f64; 8], bath: f64, dt: f64) -> ([f64; 8], f64) {
    let solve = |r: f64| {
        let mut a = [[0.0; 8]; 8];
        let mut b = [0.0; 8];
        for side in 0..2 {
            let s = if side == 0 { -1.0 } else { 1.0 };
            let g = [
                [-1.0, -1.0, -s],
                [1.0, 0.0, 0.0],
                [0.0, 1.0, 0.0],
                [0.0, 0.0, s],
            ];
            let c = 10.0 * (side + 1) as f64 / 24.0;
            for i in 0..4 {
                let row = 4 * side + i;
                a[row][row] = c;
                b[row] = c * old[row];
                for j in 0..4 {
                    a[row][4 * side + j] +=
                        dt * 2.0 / 6.0 * g[i].iter().zip(g[j]).map(|(&x, y)| x * y).sum::<f64>();
                }
            }
        }
        for i in 0..3 {
            for j in 0..3 {
                let c = dt * 0.5 * if i == j { 2.0 } else { 1.0 } / (12.0 * r);
                a[i][j] += c;
                a[4 + i][4 + j] += c;
                a[i][4 + j] -= c;
                a[4 + i][j] -= c;
            }
        }
        let physical_a = a;
        let physical_b = b;
        for i in [0, 1, 3] {
            a[i] = [0.0; 8];
            a[i][i] = 1.0;
            b[i] = bath;
        }
        for k in 0..8 {
            let p = (k..8)
                .max_by(|&i, &j| a[i][k].abs().total_cmp(&a[j][k].abs()))
                .unwrap();
            a.swap(k, p);
            b.swap(k, p);
            let pivot = a[k][k];
            for j in k..8 {
                a[k][j] /= pivot;
            }
            b[k] /= pivot;
            for i in 0..8 {
                if i != k {
                    let factor = a[i][k];
                    for j in k..8 {
                        a[i][j] -= factor * a[k][j];
                    }
                    b[i] -= factor * b[k];
                }
            }
        }
        let reaction = [0, 1, 3]
            .into_iter()
            .map(|i| {
                (physical_a[i]
                    .iter()
                    .zip(b)
                    .map(|(&x, t)| x * t)
                    .sum::<f64>()
                    - physical_b[i])
                    / dt
            })
            .sum();
        (b, reaction)
    };
    let (mut lo, mut hi) = (0.1, 0.7);
    for _ in 0..80 {
        let mid = (lo + hi) * 0.5;
        let (t, _) = solve(mid);
        let mean = [0, 1, 2, 4, 5, 6].into_iter().map(|i| t[i]).sum::<f64>() / 6.0;
        if 0.1 + 0.006 * (mean - 280.0) > mid {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    solve((lo + hi) * 0.5)
}

#[test]
fn g1_g3_surface_contact_bath_jump_counts_original_prescribed_storage_once() {
    context(&CancelGate::new_clock_free(), |cx| {
        let mesh = mesh();
        let card = surface_card();
        let curve = linear_chart();
        let initial_bc = bath(&mesh, 315.0);
        let contacts =
            TemperatureContacts::new(cx, &mesh, &initial_bc, vec![], vec![law(&mesh, &card)], 1)
                .unwrap();
        let engine = HeterogeneousEnthalpyBackwardEuler::new(
            cx,
            &mesh,
            &[
                ReferenceEnthalpyMaterial {
                    curve: &curve,
                    reference_density_kg_m3: 1.0,
                },
                ReferenceEnthalpyMaterial {
                    curve: &curve,
                    reference_density_kg_m3: 2.0,
                },
            ],
            &[0, 1],
            EnthalpyBudget {
                max_vertices: 8,
                max_elements: 2,
            },
        )
        .unwrap();
        let material = ConductivityModel::isotropic_declared(2.0).unwrap();
        let source = ScalarField::Uniform(0.0);
        let mut old_t = [300.0, 300.0, 300.0, 300.0, 360.0, 360.0, 360.0, 360.0];
        for (temperature, dt) in [(315.0, 0.2), (300.0, 0.3)] {
            let old = old_t.map(|t| 10.0 * (t - 280.0));
            let unchanged = old;
            let bc = bath(&mesh, temperature);
            let problem = ConductionProblem {
                mesh: &mesh,
                boundary: &bc,
                material: &material,
                element_materials: None,
                source: &source,
            };
            let accepted = contacts
                .advance_enthalpy(cx, &engine, problem, &old, dt, config(), coupling())
                .unwrap();
            assert_eq!(old, unchanged);
            let (expected, reaction) = oracle(old_t, temperature, dt);
            for (&actual, expected) in accepted.conduction.temperature.iter().zip(expected) {
                close(actual, expected, 2e-7);
            }
            close(accepted.physical.dirichlet_in_w, reaction, 2e-7);
            let prescribed_storage_w: f64 = [0, 1, 3]
                .into_iter()
                .map(|i| (10.0 / 24.0) * (temperature - old_t[i]) / dt)
                .sum();
            assert!(prescribed_storage_w.abs() > 50.0);
            let stored: f64 = accepted
                .conduction
                .specific_enthalpy_j_kg
                .iter()
                .zip(old)
                .zip(engine.reference_nodal_masses_kg())
                .map(|((&h, h0), &m)| m * (h - h0))
                .sum();
            close(stored, dt * reaction, 5e-8);
            close(accepted.physical.stored_energy_change_j, stored, 1e-10);
            close(accepted.physical.energy_residual_j, 0.0, 1e-8);
            assert!(accepted.physical.residual_norm_j <= coupling().residual_tolerance_j);
            assert!(accepted.iterations > 1);
            old_t.copy_from_slice(&accepted.conduction.temperature);
        }
    });
}

#[test]
fn g4_surface_contact_keeps_latent_prescription_and_work_refusals_atomic() {
    context(&CancelGate::new_clock_free(), |cx| {
        let mesh = mesh();
        let card = surface_card();
        let cold = chart(false);
        let hot = chart(true);
        let bc = bath(&mesh, 310.0);
        let contacts =
            TemperatureContacts::new(cx, &mesh, &bc, vec![], vec![law(&mesh, &card)], 1).unwrap();
        let engine = HeterogeneousEnthalpyBackwardEuler::new(
            cx,
            &mesh,
            &[
                ReferenceEnthalpyMaterial {
                    curve: &cold,
                    reference_density_kg_m3: 1.0,
                },
                ReferenceEnthalpyMaterial {
                    curve: &hot,
                    reference_density_kg_m3: 2.0,
                },
            ],
            &[0, 1],
            EnthalpyBudget {
                max_vertices: 8,
                max_elements: 2,
            },
        )
        .unwrap();
        let material = ConductivityModel::isotropic_declared(2.0).unwrap();
        let source = ScalarField::Uniform(0.0);
        let old = [80.0, 80.0, 80.0, 80.0, 450.0, 450.0, 450.0, 450.0];
        let before = old;
        let problem = ConductionProblem {
            mesh: &mesh,
            boundary: &bc,
            material: &material,
            element_materials: None,
            source: &source,
        };
        assert!(matches!(
            contacts.advance_enthalpy(cx, &engine, problem, &old, 0.2, config(), coupling()),
            Err(EnthalpyError::AmbiguousPrescribedTemperature { .. })
        ));
        let bc = bath(&mesh, 305.0);
        let problem = ConductionProblem {
            boundary: &bc,
            ..problem
        };
        assert!(
            contacts
                .advance_enthalpy(
                    cx,
                    &engine,
                    problem,
                    &old,
                    0.2,
                    config(),
                    ContactIterationConfig {
                        max_iterations: 1,
                        ..coupling()
                    }
                )
                .is_err()
        );
        let cancelled = CancelGate::new_clock_free();
        cancelled.request();
        context(&cancelled, |cancelled_cx| {
            assert!(matches!(
                contacts.advance_enthalpy(
                    cancelled_cx,
                    &engine,
                    problem,
                    &old,
                    0.2,
                    config(),
                    coupling()
                ),
                Err(EnthalpyError::Conduction(
                    fs_conduction::ConductionError::Cancelled { .. }
                ))
            ))
        });
        assert_eq!(old, before);
    });
}
