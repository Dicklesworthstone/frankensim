use super::*;

struct Law {
    nonlinear_k: bool,
    plateau: bool,
}
impl ThermalLaw for Law {
    fn temperature_k(&self, h: f64) -> Result<f64, DepositionError> {
        if !h.is_finite() || !(100.0..=1000.0).contains(&h) {
            return Err(DepositionError::MaterialDomain);
        }
        Ok(if self.plateau {
            if h < 300.0 {
                h
            } else if h <= 500.0 {
                300.0
            } else {
                h - 200.0
            }
        } else {
            h
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
        let primitive = |h| {
            let t = self.temperature_k(h)?;
            Ok::<_, DepositionError>(if !self.plateau || h < 300.0 {
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
        Ok(if self.nonlinear_k {
            self.temperature_k(h)? / 100.0
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
        [100.0, 1000.0]
    }
    fn state_at_specific_enthalpy(
        &self,
        h: f64,
    ) -> Result<ImplicitConductionState, DepositionError> {
        let slope = if self.plateau && (300.0..=500.0).contains(&h) {
            0.0
        } else {
            1.0
        };
        Ok(ImplicitConductionState {
            temperature_k: self.temperature_k(h)?,
            temperature_slope_kg_k_j: slope,
            conductivity_w_m_k: self.conductivity_w_m_k(h)?,
            conductivity_slope_kg_w_m_k_j: if self.nonlinear_k { slope / 100.0 } else { 0.0 },
        })
    }
}

fn boundary(reservoir: Option<f64>) -> Boundary {
    Boundary {
        ambient_k: 300.0,
        htc_w_m2_k: 0.0,
        source_w: 0.0,
        reservoir: reservoir.map(|temperature_k| Reservoir {
            area_m2: 1.0,
            conductance_w_k: 0.0,
            temperature_k,
        }),
    }
}
pub(super) fn setup<'a, L: ThermalLaw + ConductivityLaw>(
    law: &'a L,
    h: &[f64],
    reservoir: Option<f64>,
    extra: Vec<TimedEvent>,
) -> (DepositionSimulation<'a>, ConductionObserver<'a>) {
    let mut events = Vec::new();
    for (i, &specific_h_j_kg) in h.iter().enumerate() {
        events.push(TimedEvent {
            time_s: 0.0,
            event: Event::Birth(DepositedCell {
                id: CellId(i as u64),
                material_id: MaterialId(0),
                body_id: 0,
                segment_id: 0,
                mass_kg: 1.0,
                specific_h_j_kg,
                geometry: CellGeometry {
                    start_m: [i as f64, 0.0, 0.0],
                    end_m: [i as f64 + 1.0, 0.0, 0.0],
                    volume_m3: 1.0,
                    surface_area_m2: 2.0,
                },
                boundary: boundary(if i == 0 { reservoir } else { None }),
            }),
        });
    }
    let mut links = Vec::new();
    for i in 0..h.len() - 1 {
        events.push(TimedEvent {
            time_s: 0.0,
            event: Event::Contact(ThermalLink {
                id: i as u64,
                a: CellId(i as u64),
                b: CellId(i as u64 + 1),
                area_m2: 1.0,
                conductance_w_k: 0.0,
            }),
        });
        links.push(ConductiveLink {
            link_id: i as u64,
            half_distances_m: [0.5; 2],
            contact_resistance_m2_k_w: 0.0,
        });
    }
    events.extend(extra);
    let thermal = DepositionSimulation::new(
        events,
        &[(MaterialId(0), law)],
        DepositionBudget {
            max_cells: 1024,
            max_links: 1024,
            max_events: 4096,
            max_substeps: 100_000,
        },
    )
    .unwrap();
    let (thermal, _) = thermal
        .advance_to(
            0.0,
            StepConfig {
                max_dt_s: 1.0,
                energy_tolerance_j: 0.0,
            },
            || Ok(()),
        )
        .unwrap();
    let reservoirs = reservoir.map(|_| ConductiveReservoir {
        cell_id: CellId(0),
        half_distance_m: 1.0,
        surface_resistance_m2_k_w: 0.0,
    });
    let observer = ConductionObserver::new(
        &thermal,
        &[(MaterialId(0), law)],
        &links,
        &reservoirs.into_iter().collect::<Vec<_>>(),
        || Ok(()),
    )
    .unwrap();
    (thermal, observer)
}
fn controls(dt: f64) -> ImplicitConductionControls {
    ImplicitConductionControls {
        maximum_dt_s: dt,
        cell_tolerance_j: 1e-10,
        total_tolerance_j: 1e-9,
        ..ImplicitConductionControls::default()
    }
}
pub(super) fn near(a: f64, b: f64, tolerance: f64) {
    assert!((a - b).abs() <= tolerance, "{a:.15e} != {b:.15e}");
}

#[test]
fn stiff_pair_and_latent_crossing_use_original_storage_and_local_residual_gates() {
    let law = Law {
        nonlinear_k: false,
        plateau: false,
    };
    let (thermal, observer) = setup(&law, &[400.0, 200.0], None, vec![]);
    let chain =
        ImplicitConductionChain::new(&thermal, &observer, &[(MaterialId(0), &law)], || Ok(()))
            .unwrap();
    let mut workspace = ImplicitConductionWorkspace::new(2).unwrap();
    let (next, report) = chain
        .advance_to(&thermal, 2.0, controls(2.0), &mut workspace, || Ok(()))
        .unwrap();
    near(next.temperature_k(CellId(0)).unwrap(), 320.0, 1e-12);
    near(next.temperature_k(CellId(1)).unwrap(), 280.0, 1e-12);
    near(report.energy.external_heat_j, 0.0, 0.0);
    near(report.energy.residual_j, 0.0, 1e-12);
    assert_eq!(report.energy.substeps, 1);
    assert!(report.entropy.thermal_storage_change_j_k > 0.0);
    assert!(report.maximum_cell_residual_j <= 1e-10);
    assert_eq!(thermal.temperature_k(CellId(0)).unwrap(), 400.0);

    let phase = Law {
        nonlinear_k: false,
        plateau: true,
    };
    let (thermal, observer) = setup(&phase, &[500.0], Some(250.0), vec![]);
    let chain =
        ImplicitConductionChain::new(&thermal, &observer, &[(MaterialId(0), &phase)], || Ok(()))
            .unwrap();
    let mut workspace = ImplicitConductionWorkspace::new(1).unwrap();
    let (next, report) = chain
        .advance_to(&thermal, 10.0, controls(10.0), &mut workspace, || Ok(()))
        .unwrap();
    near(
        next.cells()[&CellId(0)].specific_h_j_kg,
        3000.0 / 11.0,
        1e-10,
    );
    near(report.energy.external_heat_j, 3000.0 / 11.0 - 500.0, 1e-9);
    assert_eq!(report.energy.substeps, 1);
    assert!(!report.entropy.has_unclassified_sources);
}

#[test]
fn analytic_nonlinear_tangent_matches_original_residual_and_time_refinement() {
    let law = Law {
        nonlinear_k: true,
        plateau: false,
    };
    let (thermal, observer) = setup(&law, &[400.0, 350.0, 320.0], Some(300.0), vec![]);
    let chain =
        ImplicitConductionChain::new(&thermal, &observer, &[(MaterialId(0), &law)], || Ok(()))
            .unwrap();
    let step = BackwardEulerStep::new(&thermal, 0.7, || Ok(())).unwrap();
    let energy = [-10.0, 3.0, 2.0];
    let mut report = empty_report();
    let trial = chain
        .evaluate(
            &step,
            &energy,
            false,
            controls(1.0),
            &mut report,
            &mut || Ok(()),
        )
        .unwrap();
    let mut workspace = ImplicitConductionWorkspace::new(3).unwrap();
    chain
        .jacobian(&thermal, trial.payload(), 0.7, &mut workspace, &mut || {
            Ok(())
        })
        .unwrap();
    for column in 0..3 {
        let mut plus = energy;
        let mut minus = energy;
        plus[column] += 1e-3;
        minus[column] -= 1e-3;
        let a = chain
            .evaluate(&step, &plus, false, controls(1.0), &mut report, &mut || {
                Ok(())
            })
            .unwrap();
        let b = chain
            .evaluate(
                &step,
                &minus,
                false,
                controls(1.0),
                &mut report,
                &mut || Ok(()),
            )
            .unwrap();
        for row in 0..3 {
            let analytic = if row == column {
                workspace.diagonal[row]
            } else if row == column + 1 {
                workspace.lower[column]
            } else if column == row + 1 {
                workspace.upper[row]
            } else {
                0.0
            };
            near(
                analytic,
                (a.residual_j()[row] - b.residual_j()[row]) / 0.002,
                1e-8,
            );
        }
    }
    let (thermal, observer) = setup(&law, &[400.0], Some(300.0), vec![]);
    let chain =
        ImplicitConductionChain::new(&thermal, &observer, &[(MaterialId(0), &law)], || Ok(()))
            .unwrap();
    let mut workspace = ImplicitConductionWorkspace::new(1).unwrap();
    let (next, report) = chain
        .advance_to(&thermal, 1.0, controls(1.0), &mut workspace, || Ok(()))
        .unwrap();
    near(
        next.temperature_k(CellId(0)).unwrap(),
        100.0 + 50_000.0_f64.sqrt(),
        1e-10,
    );
    assert_eq!(report.energy.substeps, 1);
    let exact = 300.0 / (1.0 - 0.25 * (-3.0_f64).exp());
    let mut errors = Vec::new();
    for dt in [0.1, 0.05, 0.025] {
        let (next, _) = chain
            .advance_to(&thermal, 1.0, controls(dt), &mut workspace, || Ok(()))
            .unwrap();
        errors.push((next.temperature_k(CellId(0)).unwrap() - exact).abs());
    }
    assert!(
        errors[1] < 0.6 * errors[0] && errors[2] < 0.6 * errors[1],
        "{errors:?}"
    );
}

#[test]
fn scheduled_baths_split_exactly_and_failures_do_not_publish_partial_state() {
    let law = Law {
        nonlinear_k: false,
        plateau: false,
    };
    let event = TimedEvent {
        time_s: 1.0,
        event: Event::Boundary(CellId(0), boundary(Some(500.0))),
    };
    let (thermal, observer) = setup(&law, &[400.0], Some(200.0), vec![event]);
    let chain =
        ImplicitConductionChain::new(&thermal, &observer, &[(MaterialId(0), &law)], || Ok(()))
            .unwrap();
    let mut workspace = ImplicitConductionWorkspace::new(1).unwrap();
    let (next, report) = chain
        .advance_to(&thermal, 2.0, controls(2.0), &mut workspace, || Ok(()))
        .unwrap();
    // First implicit interval cools to300; second heats to400.
    near(next.temperature_k(CellId(0)).unwrap(), 400.0, 1e-12);
    assert_eq!(report.energy.substeps, 2);
    assert_eq!(report.energy.events, 1);
    near(report.energy.external_heat_j, 0.0, 1e-12);
    assert!(
        report.entropy.thermal_storage_change_j_k - report.entropy.reservoir_entropy_input_j_k
            > 0.0
    );
    assert_eq!(
        next.cells()[&CellId(0)]
            .boundary
            .reservoir
            .unwrap()
            .conductance_w_k,
        0.0
    );
    let mut checks = 0;
    let cancelled = chain.advance_to(&thermal, 2.0, controls(2.0), &mut workspace, || {
        checks += 1;
        if checks > 30 {
            Err(DepositionError::Cancelled)
        } else {
            Ok(())
        }
    });
    assert!(matches!(
        cancelled,
        Err(ImplicitConductionError::Thermal(DepositionError::Cancelled))
    ));
    let mut limited = controls(2.0);
    limited.maximum_residual_evaluations = 1;
    assert!(matches!(
        chain.advance_to(&thermal, 2.0, limited, &mut workspace, || Ok(())),
        Err(ImplicitConductionError::Thermal(DepositionError::Budget(_)))
    ));
    assert_eq!(thermal.time_s(), 0.0);
    assert_eq!(thermal.temperature_k(CellId(0)).unwrap(), 400.0);
    assert_eq!(
        thermal.cells()[&CellId(0)]
            .boundary
            .reservoir
            .unwrap()
            .temperature_k,
        200.0
    );
}

/// A latent plateau with a falling conductivity, as occurs when a conducting
/// solid melts into a less conducting liquid. Its large-step residual need not
/// be monotone, even though both calorics and Fourier transport are passive.
struct FallingConductivity;
impl ThermalLaw for FallingConductivity {
    fn temperature_k(&self, h: f64) -> Result<f64, DepositionError> {
        Law {
            nonlinear_k: false,
            plateau: true,
        }
        .temperature_k(h)
    }
    fn maximum_temperature_slope(&self) -> f64 {
        1.0
    }
    fn specific_entropy_change_j_kg_k(
        &self,
        a: f64,
        b: f64,
    ) -> Result<Option<f64>, DepositionError> {
        Law {
            nonlinear_k: false,
            plateau: true,
        }
        .specific_entropy_change_j_kg_k(a, b)
    }
}
impl ConductivityLaw for FallingConductivity {
    fn conductivity_w_m_k(&self, h: f64) -> Result<f64, DepositionError> {
        if !(100.0..=550.0).contains(&h) {
            return Err(DepositionError::MaterialDomain);
        }
        Ok((600.0 - h) / 200.0)
    }
}
impl ImplicitConductionLaw for FallingConductivity {
    fn caloric_law(&self) -> &dyn ThermalLaw {
        self
    }
    fn conductivity_law(&self) -> &dyn ConductivityLaw {
        self
    }
    fn enthalpy_bounds_j_kg(&self) -> [f64; 2] {
        [100.0, 550.0]
    }
    fn state_at_specific_enthalpy(
        &self,
        h: f64,
    ) -> Result<ImplicitConductionState, DepositionError> {
        Ok(ImplicitConductionState {
            temperature_k: self.temperature_k(h)?,
            temperature_slope_kg_k_j: if (300.0..=500.0).contains(&h) {
                0.0
            } else {
                1.0
            },
            conductivity_w_m_k: self.conductivity_w_m_k(h)?,
            conductivity_slope_kg_w_m_k_j: -0.005,
        })
    }
}

#[test]
fn phase_dependent_transport_reduces_a_singular_interval_transactionally() {
    let law = FallingConductivity;
    let (thermal, observer) = setup(&law, &[400.0], Some(250.0), vec![]);
    let chain =
        ImplicitConductionChain::new(&thermal, &observer, &[(MaterialId(0), &law)], || Ok(()))
            .unwrap();
    let mut workspace = ImplicitConductionWorkspace::new(1).unwrap();
    let mut no_retry = controls(4.0);
    no_retry.maximum_halvings = 0;
    assert!(matches!(
        chain.advance_to(&thermal, 4.0, no_retry, &mut workspace, || Ok(())),
        Err(ImplicitConductionError::Linear(
            TridiagonalError::Singular { index: 0 }
        ))
    ));
    let (next, report) = chain
        .advance_to(&thermal, 4.0, controls(4.0), &mut workspace, || Ok(()))
        .unwrap();
    // At dt=4, dR/de=1-dt/4=0 on the starting latent branch. The dt=2
    // admissible sensible root solves h²-950h+150000+100*h_previous=0.
    let advance_exact = |h: f64| 475.0 - (75_625.0 - 100.0 * h).sqrt();
    let exact = advance_exact(advance_exact(400.0));
    near(next.cells()[&CellId(0)].specific_h_j_kg, exact, 1e-10);
    near(report.energy.external_heat_j, exact - 400.0, 1e-9);
    assert_eq!(report.energy.substeps, 2);
    assert_eq!(report.rejected_intervals, 1);
    assert!(report.maximum_cell_residual_j <= controls(4.0).cell_tolerance_j);
    assert_eq!(thermal.time_s(), 0.0);
    assert_eq!(thermal.cells()[&CellId(0)].specific_h_j_kg, 400.0);
}

#[test]
fn long_chain_resolves_heat_without_unrepresentable_far_field_increments() {
    let law = Law {
        nonlinear_k: false,
        plateau: false,
    };
    let (thermal, observer) = setup(&law, &[400.0; 96], Some(200.0), vec![]);
    let chain =
        ImplicitConductionChain::new(&thermal, &observer, &[(MaterialId(0), &law)], || Ok(()))
            .unwrap();
    let mut workspace = ImplicitConductionWorkspace::new(96).unwrap();
    let (next, report) = chain
        .advance_to(&thermal, 0.1, controls(0.1), &mut workspace, || Ok(()))
        .unwrap();
    // The uniform nearest-neighbor resolvent decays geometrically. Reflection
    // from cell95 is below 1e-200 here, while the first cells resolve actual heat.
    let ratio = 0.2 / (1.2 + 1.4_f64.sqrt());
    for i in 0..8 {
        near(
            next.temperature_k(CellId(i)).unwrap(),
            400.0 - 200.0 * ratio.powi(i as i32 + 1),
            1e-12,
        );
    }
    assert_eq!(next.cells()[&CellId(95)].specific_h_j_kg, 400.0);
    assert_eq!(report.energy.substeps, 1);
    assert_eq!(report.rejected_intervals, 0);
    assert_eq!(report.nonlinear_iterations, 1);
    assert!(report.maximum_cell_residual_j <= controls(0.1).cell_tolerance_j);
    near(
        report.energy.residual_j,
        0.0,
        controls(0.1).total_tolerance_j,
    );
    near(
        report.energy.external_heat_j,
        0.1 * (200.0 - next.temperature_k(CellId(0)).unwrap()),
        1e-12,
    );
    assert_eq!(thermal.cells()[&CellId(0)].specific_h_j_kg, 400.0);
}
