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
        .evaluate(&step, &energy, false, controls(1.0), &mut work, &mut || {
            Ok(())
        })
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
                false,
                controls(1.0),
                &mut work,
                &mut || Ok(()),
            )
            .unwrap();
        let minus = graph
            .inner
            .evaluate(&step, &minus, false, controls(1.0), &mut work, &mut || {
                Ok(())
            })
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
