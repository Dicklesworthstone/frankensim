use super::*;

struct Law {
    plateau: bool,
    nonlinear: bool,
}
impl ThermalLaw for Law {
    fn temperature_k(&self, h: f64) -> Result<f64, DepositionError> {
        if !h.is_finite() || !(100.0..=1500.0).contains(&h) {
            return Err(DepositionError::MaterialDomain);
        }
        Ok(if !self.plateau || h < 300.0 {
            h
        } else if h <= 500.0 {
            300.0
        } else {
            h - 200.0
        })
    }
    fn maximum_temperature_slope(&self) -> f64 {
        1.0
    }
    fn specific_entropy_change_j_kg_k(
        &self,
        a: f64,
        b: f64,
    ) -> Result<Option<f64>, DepositionError> {
        let primitive = |h| -> Result<f64, DepositionError> {
            let t = self.temperature_k(h)?;
            Ok(if !self.plateau || h < 300.0 {
                t.ln()
            } else if h <= 500.0 {
                300.0_f64.ln() + (h - 300.0) / 300.0
            } else {
                t.ln() + 200.0 / 300.0
            })
        };
        Ok(Some(primitive(b)? - primitive(a)?))
    }
}
impl ConductivityLaw for Law {
    fn conductivity_w_m_k(&self, h: f64) -> Result<f64, DepositionError> {
        self.temperature_k(h)?;
        Ok(if self.nonlinear {
            1.0 + 0.002 * (h - 100.0)
        } else {
            1.0
        })
    }
}
impl ImplicitConductionLaw for Law {
    fn caloric_law(&self) -> &dyn ThermalLaw {
        self
    }
    fn conductivity_law(&self) -> &dyn ConductivityLaw {
        self
    }
    fn enthalpy_bounds_j_kg(&self) -> [f64; 2] {
        [100.0, 1500.0]
    }
    fn state_at_specific_enthalpy(
        &self,
        h: f64,
    ) -> Result<ImplicitConductionState, DepositionError> {
        Ok(ImplicitConductionState {
            temperature_k: self.temperature_k(h)?,
            temperature_slope_kg_k_j: if self.plateau && (300.0..=500.0).contains(&h) {
                0.0
            } else {
                1.0
            },
            conductivity_w_m_k: self.conductivity_w_m_k(h)?,
            conductivity_slope_kg_w_m_k_j: if self.nonlinear { 0.002 } else { 0.0 },
        })
    }
}
fn id(i: usize) -> CellId {
    CellId(7 + 3 * i as u64)
}
fn setup<'a>(
    law: &'a Law,
    h: &[f64],
    mass: &[f64],
    edges: &[(usize, usize, f64)],
) -> (DepositionSimulation<'a>, ConductionObserver<'a>) {
    let mut events: Vec<_> = h
        .iter()
        .enumerate()
        .map(|(i, &h)| TimedEvent {
            time_s: 0.0,
            event: Event::Birth(DepositedCell {
                id: id(i),
                material_id: MaterialId(1),
                body_id: i as u64,
                segment_id: 0,
                mass_kg: mass[i],
                specific_h_j_kg: h,
                geometry: CellGeometry {
                    start_m: [i as f64, 0.0, 0.0],
                    end_m: [i as f64 + 1.0, 0.0, 0.0],
                    volume_m3: mass[i],
                    surface_area_m2: 10.0,
                },
                boundary: Boundary {
                    ambient_k: 300.0,
                    htc_w_m2_k: 0.0,
                    source_w: 0.0,
                    reservoir: None,
                },
            }),
        })
        .collect();
    let mut geometry = Vec::new();
    for (i, &(a, b, resistance)) in edges.iter().enumerate() {
        events.push(TimedEvent {
            time_s: 0.0,
            event: Event::Contact(ThermalLink {
                id: i as u64,
                a: id(a),
                b: id(b),
                area_m2: 1.0,
                conductance_w_k: 0.0,
            }),
        });
        geometry.push(ConductiveLink {
            link_id: i as u64,
            half_distances_m: [0.5; 2],
            contact_resistance_m2_k_w: resistance,
        });
    }
    let initial = DepositionSimulation::new(
        events,
        &[(MaterialId(1), law)],
        DepositionBudget {
            max_cells: 1024,
            max_links: 4096,
            max_events: 8192,
            max_substeps: 100_000,
        },
    )
    .unwrap();
    let (thermal, _) = initial
        .advance_to(
            0.0,
            StepConfig {
                max_dt_s: 1.0,
                energy_tolerance_j: 0.0,
            },
            || Ok(()),
        )
        .unwrap();
    let observer =
        ConductionObserver::new(&thermal, &[(MaterialId(1), law)], &geometry, &[], || Ok(()))
            .unwrap();
    (thermal, observer)
}
fn graph<'a>(
    thermal: &DepositionSimulation<'_>,
    observer: &ConductionObserver<'a>,
    law: &'a Law,
) -> ImplicitConductionGraph<'a> {
    ImplicitConductionGraph::new(
        thermal,
        observer,
        &[(MaterialId(1), law)],
        ImplicitConductionLinearControls::default(),
        || Ok(()),
    )
    .unwrap()
}
fn controls(dt: f64) -> ImplicitConductionControls {
    ImplicitConductionControls {
        maximum_dt_s: dt,
        cell_tolerance_j: 1e-9,
        total_tolerance_j: 1e-8,
        ..Default::default()
    }
}
fn near(a: f64, b: f64, tolerance: f64) {
    assert!((a - b).abs() <= tolerance, "{a:.15e} != {b:.15e}");
}

#[test]
fn stiff_star_matches_independent_unequal_capacity_endpoint_and_conserves_heat() {
    let law = Law {
        plateau: false,
        nonlinear: false,
    };
    let h = [600.0, 200.0, 300.0, 400.0];
    let masses = [3.0, 2.0, 1.0, 4.0];
    let edges = [(0, 1, 0.5), (0, 2, 1.0), (0, 3, 0.0)];
    let (thermal, observer) = setup(&law, &h, &masses, &edges);
    assert!(
        ImplicitConductionChain::new(&thermal, &observer, &[(MaterialId(1), &law)], || Ok(()))
            .is_err()
    );
    let graph = graph(&thermal, &observer, &law);
    let mut workspace = ImplicitConductionWorkspace::new(4).unwrap();
    let dt = 100.0;
    let (next, report) = graph
        .advance_to(&thermal, dt, controls(dt), &mut workspace, || Ok(()))
        .unwrap();
    // Eliminate the three scalar leaf balances analytically, independently of
    // the graph Jacobian/Krylov implementation. Heat capacities are m * cp.
    let mut numerator = masses[0] * h[0];
    let mut denominator = masses[0];
    for &(a, b, resistance) in &edges {
        assert_eq!(a, 0);
        let dg = dt / (1.0 + resistance);
        let weight = dg * masses[b] / (masses[b] + dg);
        numerator += weight * h[b];
        denominator += weight;
    }
    let center = numerator / denominator;
    near(next.temperature_k(id(0)).unwrap(), center, 1e-9);
    for &(_, b, resistance) in &edges {
        let dg = dt / (1.0 + resistance);
        let expected = (masses[b] * h[b] + dg * center) / (masses[b] + dg);
        near(next.temperature_k(id(b)).unwrap(), expected, 1e-9);
    }
    assert_eq!(report.energy.substeps, 1);
    assert!(report.linear_iterations > 0);
    near(report.energy.external_heat_j, 0.0, 0.0);
    near(report.energy.residual_j, 0.0, 1e-8);
    assert!(report.entropy.conduction_generation_j_k > 0.0);
    assert_eq!(thermal.time_s(), 0.0);
    for link in next.links().values() {
        assert_eq!(link.conductance_w_k, 0.0);
    }
}

#[test]
fn stiff_branched_network_with_parallel_contacts_needs_one_krylov_iteration() {
    let law = Law {
        plateau: false,
        nonlinear: false,
    };
    let count = 49;
    let dt = 1000.0;
    let masses: Vec<_> = (0..count)
        .map(|i| 0.5 + 0.25 * ((5 * i) % 11) as f64)
        .collect();
    let expected: Vec<_> = (0..count)
        .map(|i| 400.0 + (((7 * i) % 13) as f64 - 6.0) / dt)
        .collect();
    let mut edges = Vec::new();
    // Cell IDs do not follow a path. Every seventh branch has an additional
    // physical contact, including the opposite stored endpoint orientation.
    for i in 1..count {
        edges.push(((i - 1) / 2, i, 0.25 * (i % 5) as f64));
        if i % 7 == 0 {
            edges.push((i, (i - 1) / 2, 0.5));
        }
    }
    // G1 manufactured backward-Euler endpoint from independent pairwise heat
    // balances, not from the production Jacobian or its factorization.
    let mut initial = expected.clone();
    for &(a, b, resistance) in &edges {
        let heat = dt * (expected[a] - expected[b]) / (1.0 + resistance);
        initial[a] += heat / masses[a];
        initial[b] -= heat / masses[b];
    }
    let (thermal, observer) = setup(&law, &initial, &masses, &edges);
    let graph = ImplicitConductionGraph::new(
        &thermal,
        &observer,
        &[(MaterialId(1), &law)],
        ImplicitConductionLinearControls {
            maximum_iterations: 1,
            ..Default::default()
        },
        || Ok(()),
    )
    .unwrap();
    let (next, report) = graph
        .advance_to(
            &thermal,
            dt,
            ImplicitConductionControls {
                minimum_dt_s: dt,
                maximum_halvings: 0,
                cell_tolerance_j: 2e-8,
                ..controls(dt)
            },
            &mut ImplicitConductionWorkspace::new(count).unwrap(),
            || Ok(()),
        )
        .unwrap();
    for (i, &temperature) in expected.iter().enumerate() {
        near(next.temperature_k(id(i)).unwrap(), temperature, 2e-9);
    }
    assert_eq!(report.linear_iterations, 1);
    assert_eq!(report.energy.substeps, 1);
    assert_eq!(report.rejected_intervals, 0);
    near(report.energy.external_heat_j, 0.0, 0.0);
    near(report.energy.residual_j, 0.0, 1e-8);
    assert!(report.entropy.conduction_generation_j_k > 0.0);
    assert_eq!(thermal.time_s(), 0.0);
}

#[test]
fn stiff_cyclic_conduction_uses_sparse_ilu_with_original_energy_acceptance() {
    let law = Law {
        plateau: false,
        nonlinear: false,
    };
    let count = 49;
    let dt = 1000.0;
    let masses: Vec<_> = (0..count)
        .map(|i| 0.5 + 0.25 * ((5 * i) % 11) as f64)
        .collect();
    let expected: Vec<_> = (0..count)
        .map(|i| 400.0 + (((7 * i) % 13) as f64 - 6.0) / dt)
        .collect();
    let mut edges = Vec::new();
    for i in 0..count {
        edges.push((i, (i + 1) % count, 0.25 * (i % 5) as f64));
        if i % 7 == 0 {
            edges.push(((i + 1) % count, i, 0.5));
        }
    }
    let mut initial = expected.clone();
    for &(a, b, resistance) in &edges {
        let heat = dt * (expected[a] - expected[b]) / (1.0 + resistance);
        initial[a] += heat / masses[a];
        initial[b] -= heat / masses[b];
    }
    let (thermal, observer) = setup(&law, &initial, &masses, &edges);
    let graph = ImplicitConductionGraph::new(
        &thermal,
        &observer,
        &[(MaterialId(1), &law)],
        ImplicitConductionLinearControls {
            maximum_iterations: 8,
            ..Default::default()
        },
        || Ok(()),
    )
    .unwrap();
    let (next, report) = graph
        .advance_to(
            &thermal,
            dt,
            ImplicitConductionControls {
                minimum_dt_s: dt,
                maximum_halvings: 0,
                cell_tolerance_j: 2e-8,
                ..controls(dt)
            },
            &mut ImplicitConductionWorkspace::new(count).unwrap(),
            || Ok(()),
        )
        .unwrap();
    for (i, &temperature) in expected.iter().enumerate() {
        near(next.temperature_k(id(i)).unwrap(), temperature, 2e-9);
    }
    assert!(report.linear_iterations <= 8);
    assert_eq!(report.energy.substeps, 1);
    assert_eq!(report.rejected_intervals, 0);
    near(report.energy.external_heat_j, 0.0, 0.0);
    near(report.energy.residual_j, 0.0, 1e-8);
    assert!(report.entropy.conduction_generation_j_k > 0.0);
    assert_eq!(thermal.time_s(), 0.0);
}

#[test]
fn unsafe_sparse_pivots_retain_the_original_general_krylov_solve() {
    // Invertible nonsymmetric-capable system requiring a pivot: neither the
    // forest nor zero-fill factor may turn its unsafe first pivot into a
    // refusal of the original linear equations.
    let mut workspace = ImplicitConductionWorkspace::new(2).unwrap();
    workspace.diagonal.copy_from_slice(&[0.0, 2.0]);
    workspace.graph_entries = Some(vec![(0, 1, 1.0), (1, 0, 1.0)]);
    workspace.rhs.copy_from_slice(&[2.0, 5.0]);
    solve(
        &mut workspace,
        ImplicitConductionLinearControls::default(),
        &mut empty_report(),
        &mut || Ok(()),
    )
    .unwrap();
    near(workspace.direction[0], 1.0, 1e-12);
    near(workspace.direction[1], 2.0, 1e-12);
}

#[test]
fn branched_latent_transport_and_cyclic_nonlinear_jacobian_keep_original_calorics() {
    let law = Law {
        plateau: true,
        nonlinear: true,
    };
    let (thermal, observer) = setup(
        &law,
        &[900.0, 350.0, 400.0, 350.0, 700.0],
        &[0.5, 5.0, 5.0, 5.0, 1.0],
        &[
            (0, 1, 0.0),
            (0, 2, 0.2),
            (0, 3, 0.0),
            (1, 2, 0.0),
            (1, 2, 1.0),
        ],
    );
    let graph = graph(&thermal, &observer, &law);
    let mut workspace = ImplicitConductionWorkspace::new(5).unwrap();
    let (next, report) = graph
        .advance_to(&thermal, 10.0, controls(10.0), &mut workspace, || Ok(()))
        .unwrap();
    for i in 1..=3 {
        near(next.temperature_k(id(i)).unwrap(), 300.0, 0.0);
        assert!(next.cells()[&id(i)].specific_h_j_kg > thermal.cells()[&id(i)].specific_h_j_kg);
    }
    assert!(next.temperature_k(id(0)).unwrap() < 320.0);
    assert_eq!(next.cells()[&id(4)], thermal.cells()[&id(4)]);
    near(report.energy.residual_j, 0.0, 1e-8);

    // Unequal masses and evolving latent conductivity make J nonsymmetric.
    // Check its physical entries against the ORIGINAL endpoint residual.
    let old_patch = thermal.contact_patch(0).unwrap();
    let mut thermal = thermal
        .replace_contact_patches(
            &[ContactPatchUpdate {
                id: 0,
                expected: Some(old_patch),
                replacement: Some(ContactPatch {
                    covered_area_m2: [0.25, 0.75],
                    ..old_patch
                }),
            }],
            || Ok(()),
        )
        .unwrap();
    thermal.cells.get_mut(&id(0)).unwrap().boundary.htc_w_m2_k = 0.03;
    // Activate both parallel 1--2 paths across a real temperature difference,
    // with one side on the latent plateau and the other in sensible storage.
    thermal.cells.get_mut(&id(2)).unwrap().specific_h_j_kg = 250.0;
    let step = BackwardEulerStep::new(&thermal, 0.125, || Ok(())).unwrap();
    let mut work = empty_report();
    let energy = [1.0, 2.0, -1.0, 3.0, 0.0];
    let trial = graph
        .inner
        .evaluate(
            &step,
            &energy,
            &BTreeMap::new(),
            &[],
            false,
            controls(1.0),
            &mut work,
            &mut || Ok(()),
        )
        .unwrap();
    graph
        .inner
        .jacobian(
            &thermal,
            trial.payload(),
            step.duration_s(),
            &mut workspace,
            &mut || Ok(()),
        )
        .unwrap();
    let mut checkpoint = || Ok(());
    let checkpoint: &mut dyn FnMut() -> Result<(), DepositionError> = &mut checkpoint;
    let checkpoint = RefCell::new(checkpoint);
    let error = RefCell::new(None);
    let operator = GraphOperator {
        diagonal: &workspace.diagonal,
        entries: workspace.graph_entries.as_deref().unwrap(),
        preconditioner: None,
        checkpoint: &checkpoint,
        error: &error,
    };
    let v = [1.0, -2.0, 0.5, 4.0, -3.0];
    let w = [2.0, 1.0, -3.0, 0.25, 4.0];
    let mut jv = [0.0; 5];
    let mut jtw = [0.0; 5];
    LinearOp::apply(&operator, &v, &mut jv);
    operator.apply_transpose(&w, &mut jtw);
    near(
        jv.iter().zip(w).map(|(a, b)| a * b).sum(),
        v.iter().zip(jtw).map(|(a, b)| a * b).sum(),
        1e-12,
    );
    assert!(error.borrow().is_none());
    for column in 0..5 {
        let mut plus = energy;
        let mut minus = energy;
        let epsilon = 1e-4;
        plus[column] += epsilon;
        minus[column] -= epsilon;
        let plus = graph
            .inner
            .evaluate(
                &step,
                &plus,
                &BTreeMap::new(),
                &[],
                false,
                controls(1.0),
                &mut work,
                &mut || Ok(()),
            )
            .unwrap();
        let minus = graph
            .inner
            .evaluate(
                &step,
                &minus,
                &BTreeMap::new(),
                &[],
                false,
                controls(1.0),
                &mut work,
                &mut || Ok(()),
            )
            .unwrap();
        for row in 0..5 {
            let mut analytic = if row == column {
                workspace.diagonal[row]
            } else {
                0.0
            };
            for &(r, c, value) in workspace.graph_entries.as_ref().unwrap() {
                if r == row && c == column {
                    analytic += value;
                }
            }
            near(
                analytic,
                (plus.residual_j()[row] - minus.residual_j()[row]) / (2.0 * epsilon),
                2e-8,
            );
        }
    }
}

#[test]
fn adaptive_star_refines_toward_exact_diffusion_mode() {
    let law = Law {
        plateau: false,
        nonlinear: false,
    };
    let (thermal, observer) = setup(
        &law,
        &[400.0, 300.0, 300.0, 300.0],
        &[1.0; 4],
        &[(0, 1, 0.0), (0, 2, 0.0), (0, 3, 0.0)],
    );
    let graph = graph(&thermal, &observer, &law);
    let exact = 325.0 + 75.0 * (-2.0_f64).exp();
    let run = |tolerance| {
        let policy = adaptive::AdaptiveImplicitConductionControls {
            step: controls(0.5),
            initial_dt_s: 0.5,
            enthalpy_absolute_tolerance_j_kg: tolerance,
            temperature_absolute_tolerance_k: tolerance,
            relative_tolerance: 0.0,
            ..Default::default()
        };
        graph
            .advance_adaptive_to(
                &thermal,
                0.5,
                policy,
                &mut ImplicitConductionWorkspace::new(4).unwrap(),
                || Ok(()),
            )
            .unwrap()
    };
    let (coarse, coarse_report) = run(2.0);
    let (fine, fine_report) = run(0.02);
    let coarse_error = (coarse.temperature_k(id(0)).unwrap() - exact).abs();
    let fine_error = (fine.temperature_k(id(0)).unwrap() - exact).abs();
    assert!(
        fine_error < coarse_error / 3.0,
        "{coarse_error} -> {fine_error}"
    );
    assert!(fine_report.solve.energy.substeps > coarse_report.solve.energy.substeps);
    assert!(fine_report.accuracy_rejections > 0);
    assert!(fine_report.solve.linear_iterations > 0);
    near(fine_report.solve.energy.residual_j, 0.0, 1e-8);
    // Coarse and fine candidates share one linear-work ceiling; rejected time
    // comparisons must not restart that counter and evade the caller's budget.
    let limited = ImplicitConductionGraph::new(
        &thermal,
        &observer,
        &[(MaterialId(1), &law)],
        ImplicitConductionLinearControls {
            maximum_iterations: 4,
            ..Default::default()
        },
        || Ok(()),
    )
    .unwrap();
    let result = limited.advance_adaptive_to(
        &thermal,
        0.5,
        adaptive::AdaptiveImplicitConductionControls {
            step: controls(0.5),
            initial_dt_s: 0.5,
            enthalpy_absolute_tolerance_j_kg: 0.02,
            temperature_absolute_tolerance_k: 0.02,
            relative_tolerance: 0.0,
            ..Default::default()
        },
        &mut ImplicitConductionWorkspace::new(4).unwrap(),
        || Ok(()),
    );
    assert!(matches!(
        result,
        Err(adaptive::AdaptiveImplicitConductionError::Implicit(
            ImplicitConductionError::Thermal(DepositionError::Budget(
                "conduction Krylov iterations"
            ))
        ))
    ));
}

#[test]
fn krylov_work_limits_and_cancellation_publish_no_candidate() {
    let law = Law {
        plateau: false,
        nonlinear: true,
    };
    let (thermal, observer) = setup(
        &law,
        &[700.0, 200.0, 350.0, 400.0],
        &[1.0, 2.0, 3.0, 4.0],
        &[(0, 1, 0.0), (0, 2, 0.0), (0, 3, 0.0)],
    );
    let limited = ImplicitConductionGraph::new(
        &thermal,
        &observer,
        &[(MaterialId(1), &law)],
        ImplicitConductionLinearControls {
            maximum_iterations: 1,
            ..Default::default()
        },
        || Ok(()),
    )
    .unwrap();
    let initial = thermal.cells().clone();
    assert!(matches!(
        limited.advance_to(
            &thermal,
            5.0,
            controls(5.0),
            &mut ImplicitConductionWorkspace::new(4).unwrap(),
            || Ok(())
        ),
        Err(ImplicitConductionError::Thermal(DepositionError::Budget(
            "conduction Krylov iterations"
        )))
    ));
    let graph = graph(&thermal, &observer, &law);
    let mut polls = 0;
    let result = graph.advance_to(
        &thermal,
        5.0,
        controls(5.0),
        &mut ImplicitConductionWorkspace::new(4).unwrap(),
        || {
            polls += 1;
            if polls >= 100 {
                Err(DepositionError::Cancelled)
            } else {
                Ok(())
            }
        },
    );
    assert!(matches!(
        result,
        Err(ImplicitConductionError::Thermal(DepositionError::Cancelled))
    ));
    assert_eq!(thermal.cells(), &initial);
    assert_eq!(thermal.time_s(), 0.0);
}

#[test]
fn cyclic_anisotropic_diffusion_retains_the_original_heat_equations() {
    let side = 12;
    let n = side * side;
    let law = Law {
        plateau: false,
        nonlinear: false,
    };
    let mut initial = Vec::new();
    let mut edges = Vec::new();
    for y in 0..side {
        for x in 0..side {
            let row = y * side + x;
            let cx = (core::f64::consts::PI * (x as f64 + 0.5) / side as f64).cos();
            let cy = (core::f64::consts::PI * (y as f64 + 0.5) / side as f64).cos();
            initial.push(300.0 + 50.0 * cx + 20.0 * cy);
            if x + 1 < side {
                edges.push((row, row + 1, 0.0));
            }
            // Series contact resistance 999 gives conductance 1/1000.
            if y + 1 < side {
                edges.push((row, row + side, 999.0));
            }
        }
    }
    let (thermal, observer) = setup(&law, &initial, &vec![1.0; n], &edges);
    let graph = graph(&thermal, &observer, &law);
    let dt = 1000.0;
    let (next, report) = graph
        .advance_to(
            &thermal,
            dt,
            controls(dt),
            &mut ImplicitConductionWorkspace::new(n).unwrap(),
            || Ok(()),
        )
        .unwrap();
    // Cell-centred cosine modes are exact eigenvectors of the reflecting
    // finite-volume grid. This checks the actual backward-Euler endpoint, not
    // an independently reassembled copy of the production matrix.
    let lambda = 4.0 * (core::f64::consts::PI / (2.0 * side as f64)).sin().powi(2);
    for y in 0..side {
        for x in 0..side {
            let cx = (core::f64::consts::PI * (x as f64 + 0.5) / side as f64).cos();
            let cy = (core::f64::consts::PI * (y as f64 + 0.5) / side as f64).cos();
            let expected =
                300.0 + 50.0 * cx / (1.0 + dt * lambda) + 20.0 * cy / (1.0 + 0.001 * dt * lambda);
            near(
                next.temperature_k(id(y * side + x)).unwrap(),
                expected,
                2e-8,
            );
        }
    }
    near(report.energy.residual_j, 0.0, 1e-8);
    assert_eq!(report.rejected_intervals, 0);
    assert_eq!(report.energy.substeps, 1);
}

#[test]
fn cyclic_incomplete_factor_reduces_iterations_without_changing_true_residual() {
    let side = 16;
    let n = side * side;
    let mut diagonal = vec![1.0; n];
    let mut entries = Vec::new();
    for y in 0..side {
        for x in 0..side {
            let a = y * side + x;
            for (b, conductance) in [
                ((x + 1 < side).then_some(a + 1), 1000.0),
                ((y + 1 < side).then_some(a + side), 1.0),
            ] {
                if let Some(b) = b {
                    diagonal[a] += conductance;
                    diagonal[b] += conductance;
                    // Parallel half-faces exercise original contribution
                    // coalescing solely inside the approximate inverse.
                    for _ in 0..2 {
                        entries.push((a, b, -0.5 * conductance));
                        entries.push((b, a, -0.5 * conductance));
                    }
                }
            }
        }
    }
    let expected: Vec<_> = (0..n).map(|i| (i * 37 % 103) as f64 / 103.0).collect();
    let mut rhs: Vec<_> = diagonal.iter().zip(&expected).map(|(d, x)| d * x).collect();
    for &(row, col, value) in &entries {
        rhs[row] += value * expected[col];
    }
    let solve_with = |use_ilu| {
        let ilu = if use_ilu {
            Some(
                incomplete_factor(&diagonal, &entries, &mut || Ok(()))
                    .unwrap()
                    .unwrap(),
            )
        } else {
            None
        };
        let mut poll = || Ok(());
        let poll: &mut dyn FnMut() -> Result<(), DepositionError> = &mut poll;
        let checkpoint = RefCell::new(poll);
        let error = RefCell::new(None);
        let operator = GraphOperator {
            diagonal: &diagonal,
            entries: &entries,
            preconditioner: ilu.map(SparsePreconditioner::Incomplete),
            checkpoint: &checkpoint,
            error: &error,
        };
        let mut state = FgmresState::new(&rhs, 32);
        let result = state.run_controlled(&operator, &operator, &rhs, 1e-11, 64, &mut || false);
        assert!(result.report.converged, "{result:?}");
        assert!(error.borrow().is_none());
        for (&actual, &target) in state.x.iter().zip(&expected) {
            near(actual, target, 1e-8);
        }
        (state.iters, state.rel_residual())
    };
    let (diagonal_iterations, diagonal_residual) = solve_with(false);
    let (ilu_iterations, ilu_residual) = solve_with(true);
    assert!(
        ilu_iterations * 2 < diagonal_iterations,
        "diagonal {diagonal_iterations} ({diagonal_residual}), ILU {ilu_iterations} ({ilu_residual})"
    );
}

#[test]
fn incomplete_pivot_failure_falls_back_and_cancellation_is_preserved() {
    // Invertible original matrix with a zero first incomplete pivot.
    assert!(
        incomplete_factor(&[0.0, 2.0], &[(0, 1, 1.0), (1, 0, 1.0)], &mut || Ok(()))
            .unwrap()
            .is_none()
    );
    let mut calls = 0;
    let result = incomplete_factor(
        &[4.0; 3],
        &[(0, 1, -1.0), (1, 2, -1.0), (2, 0, -1.0)],
        &mut || {
            calls += 1;
            if calls == 4 {
                Err(DepositionError::Cancelled)
            } else {
                Ok(())
            }
        },
    );
    assert!(matches!(result, Err(DepositionError::Cancelled)));
}

use super::super::coupled_source::{ImplicitCoupledHeatSourceLaw, ImplicitCoupledHeatSourceState};
use super::super::source::{HeatSourceSchedule, HeatSourceScheduleKnot};

/// Stationary two-resistor circuit: the same current heats both original
/// caloric bodies. A resistance increase changes BOTH Joule powers.
struct SeriesResistanceSource<'a> {
    calorics: [&'a dyn ThermalLaw; 2],
    materials: [&'a dyn ImplicitConductionLaw; 2],
    resistance_ohm_at_300k: [f64; 2],
    resistance_slope_ohm_k: [f64; 2],
    voltage_v: f64,
    dependencies: Vec<(usize, usize)>,
    drive: Option<&'a HeatSourceSchedule>,
}
fn series_source(
    law: &Law,
    voltage_v: f64,
    resistance_ohm_at_300k: [f64; 2],
    resistance_slope_ohm_k: [f64; 2],
) -> SeriesResistanceSource<'_> {
    SeriesResistanceSource {
        calorics: [law, law],
        materials: [law, law],
        resistance_ohm_at_300k,
        resistance_slope_ohm_k,
        voltage_v,
        dependencies: vec![(0, 0), (0, 1), (1, 0), (1, 1)],
        drive: None,
    }
}
impl ImplicitCoupledHeatSourceLaw for SeriesResistanceSource<'_> {
    fn caloric_laws(&self) -> &[&dyn ThermalLaw] {
        &self.calorics
    }
    fn power_dependencies(&self) -> &[(usize, usize)] {
        &self.dependencies
    }
    fn state_at_specific_enthalpies_and_time(
        &self,
        h: &[f64],
        time_s: f64,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<ImplicitCoupledHeatSourceState, DepositionError> {
        checkpoint()?;
        if h.len() != 2 || !time_s.is_finite() || time_s < 0.0 {
            return Err(DepositionError::MaterialDomain);
        }
        let states = [
            self.materials[0].state_at_specific_enthalpy(h[0])?,
            self.materials[1].state_at_specific_enthalpy(h[1])?,
        ];
        let mut resistance = [0.0; 2];
        let mut slope = [0.0; 2];
        for i in 0..2 {
            resistance[i] = self.resistance_ohm_at_300k[i]
                + self.resistance_slope_ohm_k[i] * (states[i].temperature_k - 300.0);
            slope[i] = self.resistance_slope_ohm_k[i] * states[i].temperature_slope_kg_k_j;
            if !positive(resistance[i]) {
                return Err(DepositionError::MaterialDomain);
            }
        }
        let multiplier = self
            .drive
            .map_or(Ok(1.0), |drive| drive.multiplier_at(time_s))?;
        let total = resistance[0] + resistance[1];
        let current = self.voltage_v / total;
        let current_squared = multiplier * current * current;
        Ok(ImplicitCoupledHeatSourceState {
            powers_w: resistance.iter().map(|r| current_squared * r).collect(),
            power_enthalpy_derivatives_w_kg_j: self
                .dependencies
                .iter()
                .map(|&(row, col)| {
                    let direct = if row == col { slope[col] } else { 0.0 };
                    current_squared * (direct - 2.0 * resistance[row] * slope[col] / total)
                })
                .collect(),
        })
    }
    fn next_transition_time_s(&self, after_s: f64) -> Option<f64> {
        self.drive
            .and_then(|drive| drive.next_transition_time_s(after_s))
    }
}

#[test]
fn coupled_joule_feedback_and_conduction_match_independent_current_balance() {
    let law = Law {
        plateau: false,
        nonlinear: false,
    };
    let masses = [2.0, 3.0];
    let (mut thermal, observer) = setup(&law, &[300.0, 320.0], &masses, &[(0, 1, 0.0)]);
    thermal.cells.get_mut(&id(1)).unwrap().boundary.source_w = 3.0;
    let source = series_source(&law, 12.0, [2.0, 3.0], [0.01, 0.02]);
    let chain =
        ImplicitConductionChain::new(&thermal, &observer, &[(MaterialId(1), &law)], || Ok(()))
            .unwrap()
            .with_coupled_heat_source(&thermal, &[id(0), id(1)], &source, 4, || Ok(()))
            .unwrap();
    let sparse = graph(&thermal, &observer, &law)
        .with_coupled_heat_source(&thermal, &[id(0), id(1)], &source, 4, || Ok(()))
        .unwrap();
    let dt = 0.5;
    // Independently eliminate the two BE heat balances at a trial I^2=q.
    // The remaining scalar equation is I*(R0(T0)+R1(T1))=12 V.
    // G=1 W/K, cp=1 J/(kg K), and the original cell-1 source contributes 3 W.
    let temperatures = |q: f64| {
        let a = masses[0] + dt - dt * q * 0.01;
        let d = masses[1] + dt - dt * q * 0.02;
        let b0 = dt * q * 2.0;
        let b1 = masses[1] * 20.0 + dt * (q * 3.0 + 3.0);
        let determinant = a * d - dt * dt;
        [
            (d * b0 + dt * b1) / determinant,
            (dt * b0 + a * b1) / determinant,
        ]
    };
    let mut lower = 0.0;
    let mut upper = (12.0_f64 / 5.0).powi(2);
    for _ in 0..100 {
        let q = f64::midpoint(lower, upper);
        let x = temperatures(q);
        if q.sqrt() * (5.0 + 0.01 * x[0] + 0.02 * x[1]) < 12.0 {
            lower = q;
        } else {
            upper = q;
        }
    }
    let q = f64::midpoint(lower, upper);
    let expected = temperatures(q);
    let mut workspace = ImplicitConductionWorkspace::new(2).unwrap();
    let direct_result = chain
        .advance_to(&thermal, dt, controls(dt), &mut workspace, || Ok(()))
        .unwrap();
    let sparse_result = sparse
        .advance_to(&thermal, dt, controls(dt), &mut workspace, || Ok(()))
        .unwrap();
    for (next, report) in [direct_result, sparse_result] {
        for i in 0..2 {
            near(
                next.temperature_k(id(i)).unwrap(),
                300.0 + expected[i],
                1e-8,
            );
        }
        let storage = masses[0] * expected[0] + masses[1] * (expected[1] - 20.0);
        // The circuit's thermal credit is exactly VI, not VI plus I^2R again.
        near(
            report.energy.external_heat_j,
            dt * (12.0 * q.sqrt() + 3.0),
            1e-8,
        );
        near(report.energy.stored_energy_change_j, storage, 1e-8);
        near(report.energy.residual_j, 0.0, 1e-8);
        assert_eq!(report.energy.substeps, 1);
        assert!(report.entropy.conduction_generation_j_k > 0.0);
        assert!(report.entropy.has_unclassified_sources);
        near(next.cells()[&id(1)].boundary.source_w, 3.0, 0.0);
    }
    near(thermal.time_s(), 0.0, 0.0);
    near(thermal.temperature_k(id(0)).unwrap(), 300.0, 0.0);
}

#[test]
fn nonlocal_joule_tangents_reach_original_graph_action_and_transpose() {
    let law = Law {
        plateau: false,
        nonlinear: true,
    };
    let (thermal, observer) = setup(
        &law,
        &[360.0, 330.0, 410.0],
        &[1.5, 2.0, 3.0],
        &[(0, 1, 0.0), (1, 2, 0.0)],
    );
    let source = series_source(&law, 12.0, [2.0, 3.0], [0.01, 0.02]);
    // Deliberately reverse the source coordinate order. The electrical
    // coupling skips the middle thermal cell and has unequal column masses.
    let bindings = [id(2), id(0)];
    let sparse = graph(&thermal, &observer, &law)
        .with_coupled_heat_source(&thermal, &bindings, &source, 4, || Ok(()))
        .unwrap();
    let chain =
        ImplicitConductionChain::new(&thermal, &observer, &[(MaterialId(1), &law)], || Ok(()))
            .unwrap();
    assert!(matches!(
        chain.with_coupled_heat_source(&thermal, &bindings, &source, 4, || Ok(())),
        Err(ImplicitConductionError::Thermal(
            DepositionError::InvalidInput(
                "non-neighbor heat-source feedback requires the implicit graph solver",
            )
        ))
    ));

    let dt = 0.37;
    let energy = [-4.0, 2.0, 6.0];
    let step = BackwardEulerStep::new(&thermal, dt, || Ok(())).unwrap();
    let mut work = empty_report();
    let trial = sparse
        .inner
        .evaluate(
            &step,
            &energy,
            &BTreeMap::new(),
            &[],
            false,
            controls(dt),
            &mut work,
            &mut || Ok(()),
        )
        .unwrap();
    let mut workspace = ImplicitConductionWorkspace::new(3).unwrap();
    sparse
        .inner
        .jacobian(
            &thermal,
            trial.payload(),
            dt,
            &mut workspace,
            &mut || Ok(()),
        )
        .unwrap();
    let entries = workspace.graph_entries.as_deref().unwrap();
    assert!(
        entries
            .iter()
            .any(|&(row, col, value)| row == 0 && col == 2 && value > 0.0)
    );
    assert!(
        entries
            .iter()
            .any(|&(row, col, value)| row == 2 && col == 0 && value > 0.0)
    );

    let mut poll = || Ok(());
    let poll: &mut dyn FnMut() -> Result<(), DepositionError> = &mut poll;
    let checkpoint = RefCell::new(poll);
    let error = RefCell::new(None);
    let operator = GraphOperator {
        diagonal: &workspace.diagonal,
        entries,
        preconditioner: None,
        checkpoint: &checkpoint,
        error: &error,
    };
    let mut difference = [[0.0; 3]; 3];
    for column in 0..3 {
        let mut plus = energy;
        let mut minus = energy;
        plus[column] += 1e-3;
        minus[column] -= 1e-3;
        let plus = sparse
            .inner
            .evaluate(
                &step,
                &plus,
                &BTreeMap::new(),
                &[],
                false,
                controls(dt),
                &mut work,
                &mut || Ok(()),
            )
            .unwrap();
        let minus = sparse
            .inner
            .evaluate(
                &step,
                &minus,
                &BTreeMap::new(),
                &[],
                false,
                controls(dt),
                &mut work,
                &mut || Ok(()),
            )
            .unwrap();
        let mut basis = [0.0; 3];
        basis[column] = 1.0;
        let mut action = [0.0; 3];
        LinearOp::apply(&operator, &basis, &mut action);
        for row in 0..3 {
            difference[row][column] = (plus.residual_j()[row] - minus.residual_j()[row]) / 0.002;
            near(action[row], difference[row][column], 5e-8);
        }
    }
    let v = [0.3, -1.2, 2.0];
    let w = [2.0, -0.5, 0.2];
    let mut jv = [0.0; 3];
    let mut jtw = [0.0; 3];
    LinearOp::apply(&operator, &v, &mut jv);
    operator.apply_transpose(&w, &mut jtw);
    for i in 0..3 {
        near(jv[i], (0..3).map(|j| difference[i][j] * v[j]).sum(), 1e-7);
        near(jtw[i], (0..3).map(|j| difference[j][i] * w[j]).sum(), 1e-7);
    }
    near(
        jv.iter().zip(w).map(|(a, b)| a * b).sum(),
        v.iter().zip(jtw).map(|(a, b)| a * b).sum(),
        1e-12,
    );
    assert!(error.borrow().is_none());

    // A declared source dependency is admitted once before any stage. Work
    // and shape errors do not become changing or silently truncated patterns.
    assert!(
        graph(&thermal, &observer, &law)
            .with_coupled_heat_source(&thermal, &bindings, &source, 3, || Ok(()),)
            .is_err()
    );
    assert!(
        graph(&thermal, &observer, &law)
            .with_coupled_heat_source(&thermal, &[id(0), id(0)], &source, 4, || Ok(()),)
            .is_err()
    );
    for pattern in [vec![(0, 1), (0, 1)], vec![(0, 2)]] {
        let invalid = SeriesResistanceSource {
            dependencies: pattern,
            ..series_source(&law, 12.0, [2.0, 3.0], [0.01, 0.02])
        };
        assert!(
            graph(&thermal, &observer, &law)
                .with_coupled_heat_source(&thermal, &bindings, &invalid, 4, || Ok(()),)
                .is_err()
        );
    }
}

#[test]
fn coupled_generation_uses_actual_stage_times_and_adaptive_source_knots() {
    let law = Law {
        plateau: false,
        nonlinear: false,
    };
    let (thermal, observer) = setup(&law, &[300.0, 300.0], &[1.0, 1.0], &[(0, 1, 0.0)]);
    let drive = HeatSourceSchedule::new(
        [(0.0, 0.0), (0.25, 1.0), (1.0, 0.0), (1.5, 0.0)]
            .into_iter()
            .map(|(time_s, power_multiplier)| HeatSourceScheduleKnot {
                time_s,
                power_multiplier,
            })
            .collect(),
        4,
        || Ok(()),
    )
    .unwrap();
    let beta = 0.01;
    let voltage = 10.0;
    let source = SeriesResistanceSource {
        drive: Some(&drive),
        ..series_source(&law, voltage, [1.0, 1.0], [beta, beta])
    };
    let chain =
        ImplicitConductionChain::new(&thermal, &observer, &[(MaterialId(1), &law)], || Ok(()))
            .unwrap()
            .with_coupled_heat_source(&thermal, &[id(0), id(1)], &source, 4, || Ok(()))
            .unwrap();
    let sparse = graph(&thermal, &observer, &law)
        .with_coupled_heat_source(&thermal, &[id(0), id(1)], &source, 4, || Ok(()))
        .unwrap();
    let mut workspace = ImplicitConductionWorkspace::new(2).unwrap();

    // For identical bodies T0=T1=300+x, Joule power per body is
    // V^2*drive(t)/(4*(1+beta*x)). Independently integrate:
    // x + beta*x^2/2 = V^2/4 * integral(drive dt), whose pulse area is 1/2.
    let exact = 300.0 + ((1.0_f64 + beta * voltage * voltage * 0.5 / 2.0).sqrt() - 1.0) / beta;
    let direct_result = chain
        .advance_sdirk2_to(&thermal, 1.5, controls(2.0), &mut workspace, || Ok(()))
        .unwrap();
    let sparse_result = sparse
        .advance_sdirk2_to(&thermal, 1.5, controls(2.0), &mut workspace, || Ok(()))
        .unwrap();
    for (next, report) in [direct_result, sparse_result] {
        assert_eq!(report.energy.substeps, 3); // .25, 1.0, then the declared tail.
        let t = next.temperature_k(id(0)).unwrap();
        assert!(t > 300.0);
        near(next.temperature_k(id(1)).unwrap(), t, 1e-8);
        near(report.energy.external_heat_j, 2.0 * (t - 300.0), 1e-8);
        near(report.energy.residual_j, 0.0, 1e-8);
        assert!(report.entropy.has_unclassified_sources);
    }
    let mut errors = Vec::new();
    for dt in [0.125, 0.0625, 0.03125] {
        let (next, report) = chain
            .advance_sdirk2_to(&thermal, 1.5, controls(dt), &mut workspace, || Ok(()))
            .unwrap();
        let t = next.temperature_k(id(0)).unwrap();
        errors.push((t - exact).abs());
        near(report.energy.external_heat_j, 2.0 * (t - 300.0), 1e-8);
        near(report.energy.residual_j, 0.0, 1e-8);
    }
    assert!(errors[0] / errors[1] > 3.7, "{errors:?}");
    assert!(errors[1] / errors[2] > 3.7, "{errors:?}");

    let adaptive = super::super::adaptive::AdaptiveImplicitConductionControls {
        step: controls(2.0),
        initial_dt_s: 2.0,
        enthalpy_absolute_tolerance_j_kg: 1e-5,
        temperature_absolute_tolerance_k: 1e-5,
        relative_tolerance: 0.0,
        ..Default::default()
    };
    for second_order in [false, true] {
        let (next, report) = if second_order {
            sparse.advance_sdirk2_adaptive_to(&thermal, 1.5, adaptive, &mut workspace, || Ok(()))
        } else {
            chain.advance_adaptive_to(&thermal, 1.5, adaptive, &mut workspace, || Ok(()))
        }
        .unwrap();
        let t = next.temperature_k(id(0)).unwrap();
        near(t, exact, if second_order { 1e-3 } else { 0.02 });
        near(next.temperature_k(id(1)).unwrap(), t, 1e-8);
        near(report.solve.energy.external_heat_j, 2.0 * (t - 300.0), 1e-8);
        near(report.solve.energy.residual_j, 0.0, 1e-8);
        assert!(report.maximum_accepted_normalized_error <= 1.0);
        assert!(report.accuracy_rejections > 0);
    }
    // The same source refuses after its explicit time domain. No already
    // computed heating prefix or coupled-source state is published.
    assert!(matches!(
        sparse.advance_to(&thermal, 2.0, controls(2.0), &mut workspace, || Ok(())),
        Err(ImplicitConductionError::Thermal(
            DepositionError::InvalidInput(
                "heat-source query is outside the prescribed time domain",
            )
        ))
    ));
    near(thermal.time_s(), 0.0, 0.0);
    near(thermal.temperature_k(id(0)).unwrap(), 300.0, 0.0);
    near(thermal.temperature_k(id(1)).unwrap(), 300.0, 0.0);
}

struct ConstantCellCooling<'a>(&'a dyn ThermalLaw);
impl super::super::source::ImplicitHeatSourceLaw for ConstantCellCooling<'_> {
    fn caloric_law(&self) -> &dyn ThermalLaw {
        self.0
    }
    fn state_at_specific_enthalpy(
        &self,
        h: f64,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<super::super::source::ImplicitHeatSourceState, DepositionError> {
        checkpoint()?;
        self.0.temperature_k(h)?;
        Ok(super::super::source::ImplicitHeatSourceState {
            power_w: -4.0,
            power_enthalpy_derivative_w_kg_j: 0.0,
        })
    }
}

#[test]
fn cancelling_source_owners_retain_entropy_flags_through_native_stages() {
    let law = Law {
        plateau: false,
        nonlinear: false,
    };
    let (mut thermal, observer) = setup(&law, &[300.0, 300.0], &[1.0, 1.0], &[(0, 1, 0.0)]);
    for cell in thermal.cells.values_mut() {
        cell.boundary.source_w = -4.0;
    }
    // V=4 V across two 1-ohm resistors gives exactly 4 W in each body.
    // Authored cooling cancels that power bit-for-bit, but both owners act.
    let circuit = series_source(&law, 4.0, [1.0, 1.0], [0.0, 0.0]);
    let solver = graph(&thermal, &observer, &law)
        .with_coupled_heat_source(&thermal, &[id(0), id(1)], &circuit, 4, || Ok(()))
        .unwrap();
    let mut workspace = ImplicitConductionWorkspace::new(2).unwrap();
    let (_, report) = solver
        .advance_to(&thermal, 1.0, controls(1.0), &mut workspace, || Ok(()))
        .unwrap();
    near(report.energy.external_heat_j, 0.0, 0.0);
    near(report.energy.stored_energy_change_j, 0.0, 0.0);
    near(report.entropy.unclassified_source_at_cell_j_k, 0.0, 0.0);
    assert!(report.entropy.has_unclassified_sources);

    let mut clean = thermal.clone();
    for cell in clean.cells.values_mut() {
        cell.boundary.source_w = 0.0;
    }
    let ramp = HeatSourceSchedule::new(
        vec![
            HeatSourceScheduleKnot {
                time_s: 0.0,
                power_multiplier: 1.0,
            },
            HeatSourceScheduleKnot {
                time_s: 1.0,
                power_multiplier: 0.0,
            },
        ],
        2,
        || Ok(()),
    )
    .unwrap();
    let cooling = ConstantCellCooling(&law);
    let ramped_cooling = super::super::source::ScheduledHeatSource::new(&cooling, &ramp);
    let ramped_circuit = SeriesResistanceSource {
        drive: Some(&ramp),
        ..series_source(&law, 4.0, [1.0, 1.0], [0.0, 0.0])
    };
    let solver = graph(&clean, &observer, &law)
        .with_heat_sources(
            &clean,
            &[(id(0), &ramped_cooling), (id(1), &ramped_cooling)],
            || Ok(()),
        )
        .unwrap()
        .with_coupled_heat_source(&clean, &[id(0), id(1)], &ramped_circuit, 4, || Ok(()))
        .unwrap();
    let (next, report) = solver
        .advance_sdirk2_to(&clean, 1.0, controls(1.0), &mut workspace, || Ok(()))
        .unwrap();
    // Local/coupled powers cancel throughout. Both are ZERO at the second
    // stage, so only correctly retained first-stage ownership can set this flag.
    assert_eq!(report.energy.substeps, 1);
    near(report.energy.external_heat_j, 0.0, 0.0);
    near(report.energy.stored_energy_change_j, 0.0, 0.0);
    near(report.entropy.unclassified_source_at_cell_j_k, 0.0, 0.0);
    assert!(report.entropy.has_unclassified_sources);
    for i in 0..2 {
        near(next.temperature_k(id(i)).unwrap(), 300.0, 0.0);
    }
}

/// Two cells in contact over 1 m2 of the 10 m2 each has, the contact
/// authored with 0.7 W/K, each exchanging with the 300 K ambient through a
/// film of `htc` on what the contact leaves exposed; and the binding that
/// names the contact as carried outside.
fn carried_outside(
    law: &Law,
    h: [f64; 2],
    mass: [f64; 2],
    htc: f64,
) -> (DepositionSimulation<'_>, ConductionObserver<'_>) {
    let (mut thermal, _) = setup(law, &h, &mass, &[(0, 1, 0.0)]);
    thermal.links.get_mut(&0).unwrap().conductance_w_k = 0.7;
    for i in 0..2 {
        thermal.cells.get_mut(&id(i)).unwrap().boundary.htc_w_m2_k = htc;
    }
    let observer = ConductionObserver::new(&thermal, &[(MaterialId(1), law)], &[], &[], || Ok(()))
        .unwrap()
        .with_links_carried_outside(&thermal, &[0], || Ok(()))
        .unwrap();
    (thermal, observer)
}

/// `observer` with `links` also named as carried outside.
fn named_on<'a>(
    thermal: &DepositionSimulation<'_>,
    observer: &ConductionObserver<'a>,
    links: &[u64],
) -> Result<ConductionObserver<'a>, DepositionError> {
    observer.clone().with_links_carried_outside(thermal, links, || Ok(()))
}

/// The implicit graph on `observer`, or its refusal.
fn solve_on<'a>(
    thermal: &DepositionSimulation<'_>,
    observer: &ConductionObserver<'a>,
    law: &'a Law,
) -> Result<ImplicitConductionGraph<'a>, ImplicitConductionError> {
    ImplicitConductionGraph::new(
        thermal,
        observer,
        &[(MaterialId(1), law)],
        ImplicitConductionLinearControls::default(),
        || Ok(()),
    )
}

/// A contact law and a link law that refuse to be consulted: attaching
/// either to a link carried outside must fail before it is asked anything.
struct NeverAsked<'a>(&'a Law);
const ASKED: DepositionError =
    DepositionError::InvalidInput("a law on a link carried outside was consulted");
impl ImplicitContactLaw for NeverAsked<'_> {
    fn caloric_laws(&self) -> [&dyn ThermalLaw; 2] {
        [self.0; 2]
    }
    fn enthalpy_bounds_j_kg(&self) -> [[f64; 2]; 2] {
        [[100.0, 1500.0]; 2]
    }
    fn state_at_specific_enthalpies(
        &self,
        _: [f64; 2],
    ) -> Result<ImplicitContactState, DepositionError> {
        Err(ASKED)
    }
}
impl ImplicitLinkLaw for NeverAsked<'_> {
    fn caloric_laws(&self) -> [&dyn ThermalLaw; 2] {
        [self.0; 2]
    }
    fn state_at_specific_enthalpies(
        &self,
        _: [f64; 2],
        _: ConductiveLink,
        _: ContactPatch,
        _: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<ImplicitLinkState, DepositionError> {
        Err(ASKED)
    }
}

#[test]
fn a_link_carried_outside_conducts_nothing_and_keeps_the_area_it_covers() {
    let law = Law {
        plateau: false,
        nonlinear: false,
    };
    let (h, mass, htc, dt) = ([600.0, 200.0], [2.0, 3.0], 0.3, 10.0);
    let config = StepConfig {
        max_dt_s: 1.0,
        energy_tolerance_j: 1e-9,
    };
    let enthalpy = |state: &DepositionSimulation<'_>, i: usize| state.cells()[&id(i)].specific_h_j_kg;
    let mut workspace = ImplicitConductionWorkspace::new(2).unwrap();

    // Insulated cells 400 K apart. What the binding hands the explicit
    // timeline for the contact is zero, not a small number.
    let (still, observer) = carried_outside(&law, h, mass, 0.0);
    let mut asked = observer.clone();
    let handed = asked.conductances(&still, &mut || Ok(())).unwrap();
    assert_eq!(handed.len(), 1);
    assert_eq!(handed[&0].to_bits(), 0.0_f64.to_bits());
    // The authored 0.7 W/K conducts; with the binding nothing moves.
    let (authored, _) = still.advance_to(dt, config, || Ok(())).unwrap();
    assert!(authored.temperature_k(id(0)).unwrap() < 500.0);
    let (held, _, audit) = still.advance_observed(&observer, dt, config, || Ok(())).unwrap();
    for (i, start) in h.into_iter().enumerate() {
        assert_eq!(enthalpy(&held, i).to_bits(), start.to_bits());
    }
    near(audit.stored_energy_change_j, 0.0, 0.0);
    // The implicit solves, with tolerances of zero: a conductance of any
    // size would leave a residual at the start that nothing is allowed to.
    let exact = ImplicitConductionControls {
        cell_tolerance_j: 0.0,
        total_tolerance_j: 0.0,
        ..controls(dt)
    };
    let solver = graph(&still, &observer, &law);
    let euler = solver.advance_to(&still, dt, exact, &mut workspace, || Ok(()));
    let staged = solver.advance_sdirk2_to(&still, dt, exact, &mut workspace, || Ok(()));
    for (next, report) in [euler.unwrap(), staged.unwrap()] {
        for (i, start) in h.into_iter().enumerate() {
            assert_eq!(enthalpy(&next, i).to_bits(), start.to_bits());
        }
        assert_eq!(report.energy.substeps, 1);
        near(report.maximum_cell_residual_j, 0.0, 0.0);
        near(report.energy.stored_energy_change_j, 0.0, 0.0);
        near(report.energy.residual_j, 0.0, 0.0);
        near(report.entropy.conduction_generation_j_k, 0.0, 0.0);
        // The contact is still the owner's, with what was authored on it.
        near(next.links()[&0].conductance_w_k, 0.7, 0.0);
        for side_m2 in next.contact_patch(0).unwrap().covered_area_m2 {
            near(side_m2, 1.0, 0.0);
        }
    }

    // With a film each cell cools or warms alone through the 9 m2 the
    // contact leaves exposed: one backward-Euler step of
    // m (T1 - T0) = dt htc A (300 - T1), which differs by kelvins from the
    // step through all 10 m2.
    let (thermal, observer) = carried_outside(&law, h, mass, htc);
    let solver = graph(&thermal, &observer, &law);
    let (next, report) = solver
        .advance_to(&thermal, dt, controls(dt), &mut workspace, || Ok(()))
        .unwrap();
    let mut leaving_j = 0.0;
    for (i, (start, kg)) in h.into_iter().zip(mass).enumerate() {
        let alone = |area_m2: f64| {
            let r = dt * htc * area_m2;
            (kg * start + r * 300.0) / (kg + r)
        };
        let after = next.temperature_k(id(i)).unwrap();
        near(after, alone(9.0), 1e-9);
        assert!((alone(9.0) - alone(10.0)).abs() > 0.5, "{} vs {}", alone(9.0), alone(10.0));
        leaving_j += kg * (after - start);
    }
    assert_eq!(report.energy.substeps, 1);
    // The owner's external heat is the two cells' stored change. (Heat
    // between the cells would cancel in this sum; that there was none is
    // shown by each cell's own closed form above.)
    near(report.energy.external_heat_j, leaving_j, 1e-9);
    // The explicit timeline with the binding is the one whose contact was
    // authored without conductance, to the bit.
    let mut twin = thermal.clone();
    twin.links.get_mut(&0).unwrap().conductance_w_k = 0.0;
    let (reference, _) = twin.advance_to(dt, config, || Ok(())).unwrap();
    let (observed, _, _) = thermal.advance_observed(&observer, dt, config, || Ok(())).unwrap();
    for (i, start) in h.into_iter().enumerate() {
        assert_eq!(enthalpy(&observed, i).to_bits(), enthalpy(&reference, i).to_bits());
        assert!((enthalpy(&observed, i) - start).abs() > 1.0);
    }
}

#[test]
fn the_tangent_of_a_link_carried_outside_has_no_entry_and_the_exposed_area() {
    let law = Law {
        plateau: false,
        nonlinear: false,
    };
    let (mass, htc, dt) = ([2.0, 3.0], 0.3, 10.0);
    let (thermal, observer) = carried_outside(&law, [600.0, 200.0], mass, htc);
    let solver = graph(&thermal, &observer, &law);
    let step = BackwardEulerStep::new(&thermal, dt, || Ok(())).unwrap();
    let mut work = empty_report();
    let energy = [-40.0, 25.0];
    let mut residual = |energy: &[f64; 2]| {
        solver
            .inner
            .evaluate(
                &step,
                energy,
                &BTreeMap::new(),
                &[],
                false,
                controls(dt),
                &mut work,
                &mut || Ok(()),
            )
            .unwrap()
    };
    let trial = residual(&energy);
    let mut workspace = ImplicitConductionWorkspace::new(2).unwrap();
    solver
        .inner
        .jacobian(&thermal, trial.payload(), step.duration_s(), &mut workspace, &mut || Ok(()))
        .unwrap();
    // No entry joins the two cells.
    assert!(workspace.graph_entries.as_ref().unwrap().is_empty());
    let epsilon = 1e-3;
    for column in 0..2 {
        let (mut plus, mut minus) = (energy, energy);
        plus[column] += epsilon;
        minus[column] -= epsilon;
        let (plus, minus) = (residual(&plus), residual(&minus));
        for row in 0..2 {
            let analytic = if row == column { workspace.diagonal[row] } else { 0.0 };
            near(
                analytic,
                (plus.residual_j()[row] - minus.residual_j()[row]) / (2.0 * epsilon),
                1e-8,
            );
        }
    }
    // Each diagonal is 1 + dt htc A / m with the 9 m2 left exposed; with
    // the contact's area not taken off they would be 16 and 11.
    near(workspace.diagonal[0], 1.0 + dt * htc * 9.0 / mass[0], 1e-12);
    near(workspace.diagonal[1], 1.0 + dt * htc * 9.0 / mass[1], 1e-12);
}

#[test]
fn a_contact_has_one_owner_a_conductance_or_a_carrier_outside() {
    let law = Law {
        plateau: false,
        nonlinear: false,
    };
    let h = [600.0, 200.0, 300.0];
    let (thermal, bound) = setup(&law, &h, &[2.0, 3.0, 5.0], &[(0, 1, 0.0), (1, 2, 0.0)]);
    let why = |result: Result<ConductionObserver<'_>, DepositionError>| match result {
        Err(DepositionError::InvalidInput(why)) => why,
        Err(other) => panic!("{other:?}"),
        Ok(_) => panic!("admitted"),
    };
    let name = |observer, links: &[u64]| named_on(&thermal, observer, links);
    let solve = |observer| solve_on(&thermal, observer, &law);
    let second = ConductiveLink {
        link_id: 1,
        half_distances_m: [0.5; 2],
        contact_resistance_m2_k_w: 0.0,
    };
    let partial =
        ConductionObserver::new(&thermal, &[(MaterialId(1), &law)], &[second], &[], || Ok(()))
            .unwrap();
    let also_bound = "a link carried outside is also bound to a conductance";
    assert_eq!(why(name(&bound, &[0])), also_bound);
    assert_eq!(why(name(&partial, &[1])), also_bound);
    assert_eq!(why(name(&partial, &[9])), "inactive link carried outside");
    assert_eq!(why(name(&partial, &[0, 0])), "repeated link carried outside");
    let named = name(&partial, &[0]).unwrap();
    assert_eq!(why(name(&named, &[0])), "repeated link carried outside");
    // Naming none changes nothing.
    assert!(name(&partial, &[]).is_ok());

    // A contact neither bound nor named is still refused by the implicit
    // solve, and the direct chain refuses a named one: it has no link to
    // put between that pair.
    let laws: [(MaterialId, &dyn ImplicitConductionLaw); 1] = [(MaterialId(1), &law)];
    assert!(matches!(
        solve(&partial).err(),
        Some(ImplicitConductionError::Invalid(
            "thermal graph bindings changed or direct-chain topology is invalid"
        ))
    ));
    assert!(ImplicitConductionChain::new(&thermal, &bound, &laws, || Ok(())).is_ok());
    assert!(matches!(
        ImplicitConductionChain::new(&thermal, &named, &laws, || Ok(())).err(),
        Some(ImplicitConductionError::Invalid(
            "each neighboring cell pair needs exactly one link"
        ))
    ));

    // The graph takes it: the bound contact conducts and the named one
    // does not, so the first cell keeps its enthalpy to the bit while the
    // other two exchange.
    let solver = solve(&named).unwrap();
    let mut workspace = ImplicitConductionWorkspace::new(3).unwrap();
    let (next, _) = solver
        .advance_to(&thermal, 5.0, controls(5.0), &mut workspace, || Ok(()))
        .unwrap();
    assert_eq!(next.cells()[&id(0)].specific_h_j_kg.to_bits(), h[0].to_bits());
    assert!(next.temperature_k(id(1)).unwrap() > h[1] + 10.0);
    assert!(next.temperature_k(id(2)).unwrap() < h[2] - 5.0);

    // A law attached to the named contact would be a second owner. The
    // attachment looks for a bound link and refuses before asking the law;
    // on the bound contact the same law is asked, and its own answer comes
    // back.
    let never = NeverAsked(&law);
    assert!(matches!(
        solve(&named).unwrap().with_contacts(&thermal, &[(0, &never)], || Ok(())).err(),
        Some(ImplicitConductionError::Invalid("inactive implicit contact link"))
    ));
    assert!(matches!(
        solve(&named).unwrap().with_link_laws(&thermal, &[(0, &never)], || Ok(())).err(),
        Some(ImplicitConductionError::Invalid("inactive implicit flux link"))
    ));
    assert_eq!(
        solve(&named).unwrap().with_link_laws(&thermal, &[(1, &never)], || Ok(())).err(),
        Some(ImplicitConductionError::Thermal(ASKED))
    );

    // Planar layers and current metrics, in either order.
    let beside = "links carried outside are not admitted beside planar layers or current metrics";
    let layer = PlanarPhaseLayer {
        cell_id: id(2),
        footprint_area_m2: 5.0,
        density_kg_m3: 1.0,
    };
    assert_eq!(
        why(named.clone().with_planar_phase_layers(&thermal, &[layer], || Ok(()))),
        "planar phase geometry is not admitted beside links carried outside"
    );
    let layered = partial
        .clone()
        .with_planar_phase_layers(&thermal, &[layer], || Ok(()))
        .unwrap();
    assert_eq!(why(name(&layered, &[0])), beside);
    assert_eq!(
        why(named.with_current_geometry(&thermal, CurrentConductionGeometry::default(), || Ok(()))),
        "current conduction metrics are not admitted beside links carried outside"
    );
    let metrics = CurrentConductionGeometry {
        links: [(
            1,
            CurrentConductiveLink {
                area_m2: 1.0,
                half_distances_m: [0.5; 2],
            },
        )]
        .into(),
        ..CurrentConductionGeometry::default()
    };
    let current = partial.with_current_geometry(&thermal, metrics, || Ok(())).unwrap();
    assert_eq!(why(name(&current, &[0])), beside);
}

#[test]
fn an_event_on_a_link_carried_outside_is_refused_like_one_on_a_bound_link() {
    let law = Law {
        plateau: false,
        nonlinear: false,
    };
    let config = StepConfig {
        max_dt_s: 1.0,
        energy_tolerance_j: 1e-9,
    };
    let (thermal, observer) = carried_outside(&law, [600.0, 200.0], [2.0, 3.0], 0.0);
    let link = thermal.links()[&0];
    // The contact separating at 1 s, on the fixture's own schedule: the
    // explicit timeline asks the binding, which refuses.
    let mut scheduled = thermal.clone();
    let mut events = (*scheduled.schedule).clone();
    events.push(TimedEvent {
        time_s: 1.0,
        event: Event::Separation(0),
    });
    scheduled.schedule = std::rc::Rc::new(events);
    assert_eq!(
        scheduled
            .advance_observed(&observer, 2.0, config, || Ok(()))
            .err(),
        Some(DepositionError::InvalidInput("event removes bound conduction face"))
    );
    // The owner itself refuses a second contact with a link's identifier
    // ("contact requires distinct active cells and a unique ID"), so a
    // changed contact is put to the binding's hook directly, as the owner
    // of a candidate state would.
    let heard = |event: Event| observer.clone().event(&event, &thermal);
    assert_eq!(
        heard(Event::Contact(ThermalLink { area_m2: 0.5, ..link })),
        Err(DepositionError::InvalidInput("event changes bound conduction face"))
    );
    assert_eq!(
        heard(Event::Contact(ThermalLink { a: link.b, b: link.a, ..link })),
        Err(DepositionError::InvalidInput("event changes bound conduction face"))
    );
    // The same cells and area with another conductance is the same contact.
    assert_eq!(
        heard(Event::Contact(ThermalLink {
            conductance_w_k: 3.0,
            ..link
        })),
        Ok(())
    );
}
