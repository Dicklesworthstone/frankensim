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
            &BTreeMap::new(),
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
            .evaluate(
                &step,
                &plus,
                &BTreeMap::new(),
                false,
                controls(1.0),
                &mut report,
                &mut || Ok(()),
            )
            .unwrap();
        let b = chain
            .evaluate(
                &step,
                &minus,
                &BTreeMap::new(),
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

struct IntegratedLinearK<'a>(&'a Law);
impl ImplicitLinkLaw for IntegratedLinearK<'_> {
    fn caloric_laws(&self) -> [&dyn ThermalLaw; 2] { [self.0; 2] }
    fn state_at_specific_enthalpies(
        &self, h: [f64; 2], geometry: ConductiveLink, patch: ContactPatch,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<ImplicitLinkState, DepositionError> {
        checkpoint()?;
        if geometry.contact_resistance_m2_k_w != 0.0 {
            return Err(DepositionError::InvalidInput("this test law needs perfect contact"));
        }
        let a = self.0.state_at_specific_enthalpy(h[0])?;
        let b = self.0.state_at_specific_enthalpy(h[1])?;
        let scale = patch.link.area_m2 / geometry.half_distances_m.iter().sum::<f64>();
        // Independent Kirchhoff integral of k(T)=T/100 through one homogeneous link.
        Ok(ImplicitLinkState {
            conductance_w_k: scale * (a.temperature_k + b.temperature_k) / 200.0,
            power_enthalpy_derivatives_w_kg_j: [
                scale * a.temperature_k / 100.0 * a.temperature_slope_kg_k_j,
                -scale * b.temperature_k / 100.0 * b.temperature_slope_kg_k_j,
            ],
        })
    }
}

#[test]
fn integrated_link_flux_drives_chain_and_graph_with_exact_internal_heat() {
    let law = Law { nonlinear_k: true, plateau: false };
    let (thermal, observer) = setup(&law, &[400.0, 200.0], None, vec![]);
    let integrated = IntegratedLinearK(&law);
    let native = ImplicitConductionChain::new(
        &thermal, &observer, &[(MaterialId(0), &law)], || Ok(()),
    ).unwrap();
    near(native.endpoint(&thermal, &mut || Ok(())).unwrap().link_values[&0], 8.0 / 3.0, 1e-14);
    let chain = native.with_link_laws(&thermal, &[(0, &integrated)], || Ok(())).unwrap();
    near(chain.endpoint(&thermal, &mut || Ok(())).unwrap().link_values[&0], 3.0, 1e-14);
    let graph = ImplicitConductionGraph::new(
        &thermal, &observer, &[(MaterialId(0), &law)],
        ImplicitConductionLinearControls::default(), || Ok(()),
    ).unwrap().with_link_laws(&thermal, &[(0, &integrated)], || Ok(())).unwrap();
    let mut workspace = ImplicitConductionWorkspace::new(2).unwrap();
    let (a, report) = chain.advance_to(&thermal, 2.0, controls(2.0), &mut workspace, || Ok(())).unwrap();
    let (b, graph_report) = graph.advance_to(&thermal, 2.0, controls(2.0), &mut workspace, || Ok(())).unwrap();
    // Ta+Tb=600 K is conserved, so the integrated G is always 3 W/K.
    // One backward-Euler interval therefore has this independent exact result.
    for (id, expected) in [(CellId(0), 300.0 + 100.0 / 13.0), (CellId(1), 300.0 - 100.0 / 13.0)] {
        near(a.temperature_k(id).unwrap(), expected, 1e-10);
        near(b.temperature_k(id).unwrap(), expected, 1e-10);
    }
    for report in [report, graph_report] {
        assert_eq!(report.energy.substeps, 1);
        assert_eq!(report.energy.external_heat_j, 0.0);
        assert!(report.energy.residual_j.abs() < 1e-9);
        assert!(report.entropy.conduction_generation_j_k > 0.0);
        assert!(report.maximum_cell_residual_j <= 1e-10);
    }
    assert_eq!(thermal.time_s(), 0.0);
    assert_eq!(thermal.temperature_k(CellId(0)).unwrap(), 400.0);
}

#[test]
fn integrated_link_tangent_matches_original_energy_residual_with_unequal_masses() {
    let law = Law { nonlinear_k: true, plateau: false };
    let (mut thermal, mut observer) = setup(&law, &[400.0, 320.0, 250.0], None, vec![]);
    for (id, mass) in [(CellId(1), 2.0), (CellId(2), 0.5)] {
        thermal.cells.get_mut(&id).unwrap().mass_kg = mass;
        observer.cells.get_mut(&id).unwrap().mass_kg = mass;
    }
    let integrated = IntegratedLinearK(&law);
    let chain = ImplicitConductionChain::new(
        &thermal, &observer, &[(MaterialId(0), &law)], || Ok(()),
    ).unwrap().with_link_laws(&thermal, &[(0, &integrated)], || Ok(())).unwrap();
    // One integrated link and one unchanged bulk-series link share each row.
    let step = BackwardEulerStep::new(&thermal, 0.3, || Ok(())).unwrap();
    let energy = [-7.0, 4.0, 3.0];
    let mut report = empty_report();
    let trial = chain.evaluate(&step, &energy, &BTreeMap::new(), false, controls(1.0),
        &mut report, &mut || Ok(())).unwrap();
    let mut workspace = ImplicitConductionWorkspace::new(3).unwrap();
    chain.jacobian(&thermal, trial.payload(), 0.3, &mut workspace, &mut || Ok(())).unwrap();
    for column in 0..3 {
        let mut plus = energy; let mut minus = energy;
        plus[column] += 1e-3; minus[column] -= 1e-3;
        let a = chain.evaluate(&step, &plus, &BTreeMap::new(), false, controls(1.0),
            &mut report, &mut || Ok(())).unwrap();
        let b = chain.evaluate(&step, &minus, &BTreeMap::new(), false, controls(1.0),
            &mut report, &mut || Ok(())).unwrap();
        for row in 0..3 {
            let derivative = if row == column { workspace.diagonal[row] }
                else if row == column + 1 { workspace.lower[column] }
                else if column == row + 1 { workspace.upper[row] }
                else { 0.0 };
            near(derivative, (a.residual_j()[row] - b.residual_j()[row]) / 0.002, 1e-8);
        }
    }
}

#[test]
fn integrated_link_attachment_and_late_failure_preserve_original_history() {
    let law = Law { nonlinear_k: true, plateau: false };
    let other = Law { nonlinear_k: true, plateau: false };
    let (thermal, observer) = setup(&law, &[400.0, 200.0, 300.0], None, vec![]);
    let integrated = IntegratedLinearK(&law);
    let foreign = IntegratedLinearK(&other);
    let make = || ImplicitConductionChain::new(
        &thermal, &observer, &[(MaterialId(0), &law)], || Ok(()),
    ).unwrap();
    assert!(make().with_link_laws(&thermal, &[(0, &foreign)], || Ok(())).is_err());
    assert!(make().with_link_laws(&thermal, &[(99, &integrated)], || Ok(())).is_err());
    assert!(make().with_link_laws(&thermal, &[(0, &integrated), (0, &integrated)], || Ok(())).is_err());
    assert!(make().with_link_laws(&thermal, &[], || Ok(())).is_err());
    struct Contact<'a>(&'a Law);
    impl ImplicitContactLaw for Contact<'_> {
        fn caloric_laws(&self) -> [&dyn ThermalLaw; 2] { [self.0; 2] }
        fn enthalpy_bounds_j_kg(&self) -> [[f64; 2]; 2] { [[100.0, 1000.0]; 2] }
        fn state_at_specific_enthalpies(&self, _: [f64; 2]) -> Result<ImplicitContactState, DepositionError> {
            Ok(ImplicitContactState { resistance_m2_k_w: 1.0,
                resistance_enthalpy_derivatives_m2_kg_k_w_j: [0.0; 2] })
        }
    }
    let contact = Contact(&law);
    assert!(make().with_contacts(&thermal, &[(0, &contact)], || Ok(())).unwrap()
        .with_link_laws(&thermal, &[(0, &integrated)], || Ok(())).is_err());
    assert!(make().with_link_laws(&thermal, &[(0, &integrated)], || Ok(())).unwrap()
        .with_contacts(&thermal, &[(0, &contact)], || Ok(())).is_err());
    let chain = make().with_link_laws(&thermal, &[(0, &integrated)], || Ok(())).unwrap();
    let mut workspace = ImplicitConductionWorkspace::new(3).unwrap();
    let mut polls = 0;
    let (expected, _) = chain.advance_to(&thermal, 0.2, controls(0.1), &mut workspace,
        || { polls += 1; Ok(()) }).unwrap();
    let mut calls = 0;
    assert!(chain.advance_to(&thermal, 0.2, controls(0.1), &mut workspace,
        || { calls += 1; if calls == polls { Err(DepositionError::Cancelled) } else { Ok(()) } }).is_err());
    assert_eq!(thermal.time_s(), 0.0);
    assert_eq!(thermal.temperature_k(CellId(0)).unwrap(), 400.0);
    let (retry, _) = chain.advance_to(&thermal, 0.2, controls(0.1), &mut workspace, || Ok(())).unwrap();
    for id in chain.cell_ids() { assert_eq!(retry.cells()[id], expected.cells()[id]); }
    for state in [
        ImplicitLinkState { conductance_w_k: -1.0, power_enthalpy_derivatives_w_kg_j: [1.0, -1.0] },
        ImplicitLinkState { conductance_w_k: 1.0, power_enthalpy_derivatives_w_kg_j: [f64::NAN, -1.0] },
    ] { assert!(state.admit().is_err()); }
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

#[derive(Clone)]
struct Forcing {
    time_s: f64,
    constant_w: f64,
    interval_w: f64,
    ramp_w_s: f64,
    transitions: Vec<f64>,
    maximum_dt_s: f64,
    maximum_evaluations: usize,
    evaluations: usize,
    begins: usize,
    events: usize,
    energy_j: f64,
    pending: Option<(f64, f64, f64)>,
    accepted: Vec<(f64, f64, f64)>,
    unsupported: u8,
    refuse_at_s: f64,
}
impl Default for Forcing {
    fn default() -> Self {
        Self {
            time_s: 0.0,
            constant_w: 0.0,
            interval_w: 0.0,
            ramp_w_s: 0.0,
            transitions: vec![],
            maximum_dt_s: f64::MAX,
            maximum_evaluations: usize::MAX,
            evaluations: 0,
            begins: 0,
            events: 0,
            energy_j: 0.0,
            pending: None,
            accepted: vec![],
            unsupported: 0,
            refuse_at_s: f64::INFINITY,
        }
    }
}
impl DepositionObserver for Forcing {
    type Error = DepositionError;
    fn begin(&mut self) -> Result<(), Self::Error> {
        self.begins += 1;
        self.evaluations = 0;
        self.pending = None;
        Ok(())
    }
    fn requires_step_interval(&self) -> bool {
        true
    }
    fn next_transition_time_s(
        &mut self,
        _: &DepositionSimulation<'_>,
        _: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<Option<f64>, Self::Error> {
        Ok(self
            .transitions
            .iter()
            .copied()
            .find(|&time| time > self.time_s))
    }
    fn sources(
        &mut self,
        thermal: &DepositionSimulation<'_>,
        _: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<BTreeMap<CellId, f64>, Self::Error> {
        assert_eq!(thermal.time_s(), self.time_s);
        self.pending = None;
        Ok(BTreeMap::from([(CellId(0), self.constant_w)]))
    }
    fn conductances(
        &mut self,
        _: &DepositionSimulation<'_>,
        _: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<BTreeMap<u64, f64>, Self::Error> {
        Ok(if self.unsupported == 1 {
            BTreeMap::from([(0, 1.0)])
        } else {
            BTreeMap::new()
        })
    }
    fn reservoir_conductances(
        &mut self,
        _: &DepositionSimulation<'_>,
        _: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<BTreeMap<CellId, f64>, Self::Error> {
        Ok(if self.unsupported == 2 {
            BTreeMap::from([(CellId(0), 1.0)])
        } else {
            BTreeMap::new()
        })
    }
    fn heat_transfers(
        &mut self,
        _: &DepositionSimulation<'_>,
        _: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<Vec<HeatTransfer>, Self::Error> {
        Ok(if self.unsupported == 3 {
            vec![HeatTransfer {
                a: CellId(0),
                b: CellId(1),
                power_w: 1.0,
            }]
        } else {
            vec![]
        })
    }
    fn limit_timestep_s(
        &mut self,
        _: &DepositionSimulation<'_>,
        proposed_s: f64,
        _: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<f64, Self::Error> {
        Ok(proposed_s.min(self.maximum_dt_s))
    }
    fn step_sources(
        &mut self,
        thermal: &DepositionSimulation<'_>,
        start: f64,
        dt: f64,
        end: f64,
        _: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<BTreeMap<CellId, f64>, Self::Error> {
        assert_eq!(thermal.time_s(), self.time_s);
        assert_eq!(start, self.time_s);
        assert_eq!(dt, end - start);
        if self.evaluations >= self.maximum_evaluations {
            return Err(DepositionError::Budget("forcing evaluations"));
        }
        self.evaluations += 1;
        let mean = self.interval_w + self.ramp_w_s * f64::midpoint(start, end);
        self.pending = Some((end, dt, mean));
        Ok(BTreeMap::from([(CellId(0), mean)]))
    }
    fn event(&mut self, _: &Event, thermal: &DepositionSimulation<'_>) -> Result<(), Self::Error> {
        assert_eq!(self.time_s, thermal.time_s());
        self.events += 1;
        Ok(())
    }
    fn step(
        &mut self,
        start: f64,
        dt: f64,
        old: &BTreeMap<CellId, f64>,
        thermal: &DepositionSimulation<'_>,
        _: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), Self::Error> {
        let (end, duration, mean) = self.pending.unwrap();
        assert_eq!(self.time_s, start);
        assert_eq!(dt, duration);
        assert_eq!(end, thermal.time_s());
        assert!(old.contains_key(&CellId(0)));
        if end >= self.refuse_at_s {
            return Err(DepositionError::MaterialDomain);
        }
        self.energy_j += (self.constant_w + mean) * dt;
        self.accepted.push((start, end, mean));
        self.pending = None;
        self.time_s = end;
        Ok(())
    }
}

#[test]
fn observed_stiff_sources_use_native_chain_graph_and_robin_endpoints() {
    let law = Law {
        nonlinear_k: false,
        plateau: false,
    };
    let (thermal, conduction) = setup(&law, &[400.0, 200.0], None, vec![]);
    let chain =
        ImplicitConductionChain::new(&thermal, &conduction, &[(MaterialId(0), &law)], || Ok(()))
            .unwrap();
    let graph = ImplicitConductionGraph::new(
        &thermal,
        &conduction,
        &[(MaterialId(0), &law)],
        ImplicitConductionLinearControls::default(),
        || Ok(()),
    )
    .unwrap();
    let forcing = Forcing {
        constant_w: 3.0,
        interval_w: 7.0,
        ..Forcing::default()
    };
    let mut workspace = ImplicitConductionWorkspace::new(2).unwrap();
    let (next, source, report) = chain
        .advance_observed_to(
            &thermal,
            &forcing,
            10.0,
            controls(10.0),
            &mut workspace,
            || Ok(()),
        )
        .unwrap();
    // C=G=1, dt=10 exceeds the native explicit monotone limit by 12.5x.
    let sum = 600.0 + 100.0;
    let difference = (200.0 + 100.0) / 21.0;
    near(
        next.temperature_k(CellId(0)).unwrap(),
        (sum + difference) * 0.5,
        1e-10,
    );
    near(
        next.temperature_k(CellId(1)).unwrap(),
        (sum - difference) * 0.5,
        1e-10,
    );
    near(report.energy.external_heat_j, 100.0, 1e-10);
    near(source.energy_j, 100.0, 0.0);
    assert_eq!(report.energy.substeps, 1);
    assert_eq!(source.evaluations, 1);
    assert_eq!(source.begins, 1);
    assert!(report.entropy.has_unclassified_sources);
    let (graph_next, graph_source, graph_report) = graph
        .advance_observed_to(
            &thermal,
            &forcing,
            10.0,
            controls(10.0),
            &mut workspace,
            || Ok(()),
        )
        .unwrap();
    for id in [CellId(0), CellId(1)] {
        near(
            graph_next.temperature_k(id).unwrap(),
            next.temperature_k(id).unwrap(),
            1e-10,
        );
    }
    near(graph_source.energy_j, source.energy_j, 0.0);
    near(graph_report.energy.residual_j, 0.0, 1e-9);

    let mut robin = boundary(None);
    robin.htc_w_m2_k = 2.0;
    robin.source_w = 2.0;
    let (thermal, conduction) = setup(
        &law,
        &[400.0],
        None,
        vec![TimedEvent {
            time_s: 0.0,
            event: Event::Boundary(CellId(0), robin),
        }],
    );
    let chain =
        ImplicitConductionChain::new(&thermal, &conduction, &[(MaterialId(0), &law)], || Ok(()))
            .unwrap();
    let forcing = Forcing {
        constant_w: 3.0,
        interval_w: 5.0,
        ..Forcing::default()
    };
    let mut workspace = ImplicitConductionWorkspace::new(1).unwrap();
    let (next, source, report) = chain
        .advance_observed_to(
            &thermal,
            &forcing,
            10.0,
            controls(10.0),
            &mut workspace,
            || Ok(()),
        )
        .unwrap();
    let expected = (400.0 + 10.0 * (4.0 * 300.0 + 10.0)) / 41.0;
    near(next.temperature_k(CellId(0)).unwrap(), expected, 1e-10);
    near(report.energy.external_heat_j, expected - 400.0, 1e-9);
    near(source.energy_j, 80.0, 0.0);
    assert_eq!(next.cells()[&CellId(0)].boundary, robin);
    assert_eq!(thermal.time_s(), 0.0);
    assert_eq!(forcing.energy_j, 0.0);
}

#[test]
fn observed_ramp_retries_recompute_means_without_resetting_source_work() {
    let law = FallingConductivity;
    let (thermal, conduction) = setup(&law, &[400.0], Some(250.0), vec![]);
    let chain =
        ImplicitConductionChain::new(&thermal, &conduction, &[(MaterialId(0), &law)], || Ok(()))
            .unwrap();
    let source = Forcing {
        ramp_w_s: 1.0,
        maximum_evaluations: 3,
        ..Forcing::default()
    };
    let mut workspace = ImplicitConductionWorkspace::new(1).unwrap();
    let (next, accepted, report) = chain
        .advance_observed_to(
            &thermal,
            &source,
            4.0,
            controls(4.0),
            &mut workspace,
            || Ok(()),
        )
        .unwrap();
    // The dt=4 latent-branch Jacobian is singular; each dt=2 root includes
    // its own ramp mean. Only the two accepted intervals contribute joules.
    let root = |h: f64, mean: f64| 475.0 - (75_625.0 - 100.0 * h - 200.0 * mean).sqrt();
    let expected = root(root(400.0, 1.0), 3.0);
    near(next.cells()[&CellId(0)].specific_h_j_kg, expected, 1e-10);
    near(report.energy.external_heat_j, expected - 400.0, 1e-9);
    assert_eq!(report.rejected_intervals, 1);
    assert_eq!(accepted.evaluations, 3);
    assert_eq!(accepted.begins, 1);
    assert_eq!(accepted.accepted, [(0.0, 2.0, 1.0), (2.0, 4.0, 3.0)]);
    assert_eq!(accepted.energy_j, 8.0);
    let limited = Forcing {
        maximum_evaluations: 2,
        ..source.clone()
    };
    assert!(matches!(
        chain.advance_observed_to(
            &thermal,
            &limited,
            4.0,
            controls(4.0),
            &mut workspace,
            || Ok(())
        ),
        Err(ImplicitConductionObservedError::Observer(
            DepositionError::Budget("forcing evaluations")
        ))
    ));
    assert_eq!(thermal.time_s(), 0.0);
    assert_eq!(thermal.cells()[&CellId(0)].specific_h_j_kg, 400.0);
    assert_eq!(limited.energy_j, 0.0);
    assert_eq!(limited.evaluations, 0);
    assert_eq!(limited.begins, 0);
}

#[test]
fn observed_transition_events_source_limit_and_restart_keep_exact_history() {
    let law = Law {
        nonlinear_k: false,
        plateau: false,
    };
    let mut changed = boundary(None);
    changed.source_w = 2.0;
    let (thermal, conduction) = setup(
        &law,
        &[300.0],
        None,
        vec![TimedEvent {
            time_s: 1.0,
            event: Event::Boundary(CellId(0), changed),
        }],
    );
    let chain =
        ImplicitConductionChain::new(&thermal, &conduction, &[(MaterialId(0), &law)], || Ok(()))
            .unwrap();
    let source = Forcing {
        ramp_w_s: 1.0,
        transitions: vec![1.0, 3.0],
        maximum_dt_s: 2.0,
        ..Forcing::default()
    };
    let mut workspace = ImplicitConductionWorkspace::new(1).unwrap();
    let (next, accepted, report) = chain
        .advance_observed_to(
            &thermal,
            &source,
            4.0,
            controls(10.0),
            &mut workspace,
            || Ok(()),
        )
        .unwrap();
    assert_eq!(
        accepted.accepted,
        [(0.0, 1.0, 0.5), (1.0, 3.0, 2.0), (3.0, 4.0, 3.5)]
    );
    assert_eq!(accepted.energy_j, 8.0);
    assert_eq!(accepted.events, 1);
    assert_eq!(accepted.begins, 1);
    assert_eq!(report.energy.events, 1);
    near(next.temperature_k(CellId(0)).unwrap(), 314.0, 1e-10);
    near(report.energy.external_heat_j, 14.0, 1e-10);
    let (half, half_source, _) = chain
        .advance_observed_to(
            &thermal,
            &source,
            1.0,
            controls(10.0),
            &mut workspace,
            || Ok(()),
        )
        .unwrap();
    let (restarted, restarted_source, _) = chain
        .advance_observed_to(
            &half,
            &half_source,
            4.0,
            controls(10.0),
            &mut workspace,
            || Ok(()),
        )
        .unwrap();
    assert_eq!(restarted.cells(), next.cells());
    assert_eq!(restarted_source.accepted, accepted.accepted);
    assert_eq!(restarted_source.energy_j, accepted.energy_j);
    assert_eq!(restarted_source.events, accepted.events);
    assert_eq!(restarted_source.begins, 2);

    // A separate no-transition interval proves that the source's own bound is
    // honored even though the implicit conductive bound has been removed.
    let limited = Forcing {
        maximum_dt_s: 0.25,
        ..Forcing::default()
    };
    let (_, limited, report) = chain
        .advance_observed_to(
            &thermal,
            &limited,
            0.5,
            controls(10.0),
            &mut workspace,
            || Ok(()),
        )
        .unwrap();
    assert_eq!(report.energy.substeps, 2);
    assert_eq!(limited.accepted, [(0.0, 0.25, 0.0), (0.25, 0.5, 0.0)]);
}

#[test]
fn observed_refusals_preserve_native_state_and_source_history() {
    let law = Law {
        nonlinear_k: false,
        plateau: false,
    };
    let (thermal, conduction) = setup(&law, &[300.0, 300.0], None, vec![]);
    let chain =
        ImplicitConductionChain::new(&thermal, &conduction, &[(MaterialId(0), &law)], || Ok(()))
            .unwrap();
    let mut workspace = ImplicitConductionWorkspace::new(2).unwrap();
    let source = Forcing {
        constant_w: 1.0,
        refuse_at_s: 2.0,
        ..Forcing::default()
    };
    assert!(matches!(
        chain.advance_observed_to(
            &thermal,
            &source,
            2.0,
            controls(1.0),
            &mut workspace,
            || Ok(())
        ),
        Err(ImplicitConductionObservedError::Observer(
            DepositionError::MaterialDomain
        ))
    ));
    for unsupported in 1..=3 {
        let input = Forcing {
            unsupported,
            ..source.clone()
        };
        assert!(matches!(
            chain.advance_observed_to(&thermal, &input, 1.0, controls(1.0), &mut workspace, || Ok(
                ()
            )),
            Err(ImplicitConductionObservedError::Implicit(
                ImplicitConductionError::Invalid(
                    "observed implicit conduction accepts source-only observers"
                )
            ))
        ));
    }
    let mut calls = 0;
    assert!(matches!(
        chain.advance_observed_to(
            &thermal,
            &source,
            1.0,
            controls(1.0),
            &mut workspace,
            || {
                calls += 1;
                if calls == 40 {
                    Err(DepositionError::Cancelled)
                } else {
                    Ok(())
                }
            }
        ),
        Err(ImplicitConductionObservedError::Implicit(
            ImplicitConductionError::Thermal(DepositionError::Cancelled)
        ))
    ));
    assert_eq!(thermal.time_s(), 0.0);
    assert_eq!(thermal.temperature_k(CellId(0)).unwrap(), 300.0);
    assert_eq!(thermal.cells()[&CellId(0)].boundary.source_w, 0.0);
    assert_eq!(source.time_s, 0.0);
    assert_eq!(source.energy_j, 0.0);
    assert_eq!(source.evaluations, 0);
    assert_eq!(source.begins, 0);
}

#[path = "tests/ambient.rs"]
mod ambient;

