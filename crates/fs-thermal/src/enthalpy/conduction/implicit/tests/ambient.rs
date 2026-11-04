//! One caloric cell can conduct through two separately owned physical baths.
use super::super::boundary::{
    ImplicitBoundaryGeometry, ImplicitBoundaryLaw, ImplicitBoundaryState,
};
use super::*;
use crate::observers::{CoupledObservers, ObserverCouplingError};

// Independent exact integral for the existing test law k(T)=T/100. The
// production material adapter uses its original piecewise-linear source owner.
struct IntegratedBoundary<'a>(&'a Law);
impl ImplicitBoundaryLaw for IntegratedBoundary<'_> {
    fn caloric_law(&self) -> &dyn ThermalLaw {
        self.0
    }

    fn state_at_specific_enthalpy(
        &self,
        h: f64,
        bath: f64,
        geometry: ImplicitBoundaryGeometry,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<ImplicitBoundaryState, DepositionError> {
        checkpoint()?;
        if geometry.surface_resistance_m2_k_w != 0.0 {
            return Err(DepositionError::MaterialDomain);
        }
        let state = self.0.state_at_specific_enthalpy(h)?;
        let metric = geometry.area_m2 / geometry.half_distance_m;
        Ok(ImplicitBoundaryState {
            conductance_w_k: metric * (state.temperature_k + bath) / 200.0,
            power_enthalpy_derivative_w_kg_j: metric * state.temperature_k / 100.0
                * state.temperature_slope_kg_k_j,
        })
    }
}

#[test]
fn integrated_boundary_fluxes_share_exact_energy_jacobian_and_bath_entropy() {
    let law = Law {
        nonlinear_k: true,
        plateau: false,
    };
    let exact = IntegratedBoundary(&law);
    let boundary = baths(700.0, 300.0);
    let (mut thermal, _) = setup(
        &law,
        &[500.0],
        Some(300.0),
        vec![TimedEvent {
            time_s: 0.0,
            event: Event::Boundary(CellId(0), boundary),
        }],
    );
    let mass = 2.5;
    thermal.cells.get_mut(&CellId(0)).unwrap().mass_kg = mass;
    let conduction = ConductionObserver::new(
        &thermal,
        &[(MaterialId(0), &law)],
        &[],
        &[ConductiveReservoir {
            cell_id: CellId(0),
            half_distance_m: 1.0,
            surface_resistance_m2_k_w: 0.0,
        }],
        || Ok(()),
    )
    .unwrap()
    .with_ambient_faces(&thermal, &[ambient(0.0)], || Ok(()))
    .unwrap();
    let chain =
        ImplicitConductionChain::new(&thermal, &conduction, &[(MaterialId(0), &law)], || Ok(()))
            .unwrap()
            .with_boundary_laws(
                &thermal,
                &[(CellId(0), &exact)],
                &[(CellId(0), &exact)],
                || Ok(()),
            )
            .unwrap();
    let graph = ImplicitConductionGraph::new(
        &thermal,
        &conduction,
        &[(MaterialId(0), &law)],
        ImplicitConductionLinearControls::default(),
        || Ok(()),
    )
    .unwrap()
    .with_boundary_laws(
        &thermal,
        &[(CellId(0), &exact)],
        &[(CellId(0), &exact)],
        || Ok(()),
    )
    .unwrap();
    let dt = 0.25;
    let mut workspace = ImplicitConductionWorkspace::new(1).unwrap();
    let endpoint = chain.endpoint(&thermal, &mut || Ok(())).unwrap();
    chain
        .jacobian(&thermal, &endpoint, dt, &mut workspace, &mut || Ok(()))
        .unwrap();
    near(
        workspace.diagonal[0],
        1.0 + dt * 3.0 * 500.0 / (100.0 * mass),
        1e-14,
    );
    // m(T-T0)=dt*[.01*(Ta²-T²)+.005*(Tr²-T²)], solved independently.
    let rhs = mass * 500.0 + dt * (0.01 * 700.0_f64.powi(2) + 0.005 * 300.0_f64.powi(2));
    let quadratic = dt * 0.015;
    let expected = 2.0 * rhs / (mass + (mass * mass + 4.0 * quadratic * rhs).sqrt());
    for result in [
        chain.advance_to(&thermal, dt, controls(dt), &mut workspace, || Ok(())),
        graph.advance_to(&thermal, dt, controls(dt), &mut workspace, || Ok(())),
    ] {
        let (next, report) = result.unwrap();
        near(next.temperature_k(CellId(0)).unwrap(), expected, 1e-10);
        let ambient_heat = dt * 0.01 * (700.0_f64.powi(2) - expected * expected);
        let reservoir_heat = dt * 0.005 * (300.0_f64.powi(2) - expected * expected);
        near(
            report.energy.external_heat_j,
            ambient_heat + reservoir_heat,
            1e-9,
        );
        near(
            report.energy.stored_energy_change_j,
            mass * (expected - 500.0),
            1e-9,
        );
        near(report.energy.residual_j, 0.0, 1e-9);
        near(
            report.entropy.reservoir_entropy_input_j_k,
            ambient_heat / 700.0 + reservoir_heat / 300.0,
            1e-11,
        );
        assert_eq!(next.cells[&CellId(0)].boundary, boundary);
        assert_eq!(report.energy.substeps, 1);
    }
}

#[test]
fn integrated_boundary_binding_refuses_other_calorics_and_duplicate_channels() {
    let law = Law {
        nonlinear_k: true,
        plateau: false,
    };
    let different = Law {
        nonlinear_k: true,
        plateau: false,
    };
    let exact = IntegratedBoundary(&law);
    let unrelated = IntegratedBoundary(&different);
    let (thermal, conduction) = setup(&law, &[500.0], Some(300.0), vec![]);
    let build = || {
        ImplicitConductionChain::new(&thermal, &conduction, &[(MaterialId(0), &law)], || Ok(()))
            .unwrap()
    };
    assert!(
        build()
            .with_boundary_laws(&thermal, &[(CellId(0), &unrelated)], &[], || Ok(()))
            .is_err()
    );
    assert!(
        build()
            .with_boundary_laws(
                &thermal,
                &[(CellId(0), &exact), (CellId(0), &exact)],
                &[],
                || Ok(())
            )
            .is_err()
    );
    assert!(
        build()
            .with_boundary_laws(&thermal, &[], &[(CellId(0), &exact)], || Ok(()))
            .is_err()
    );
    assert!(
        build()
            .with_boundary_laws(&thermal, &[(CellId(9), &exact)], &[], || Ok(()))
            .is_err()
    );
    let chain = build()
        .with_boundary_laws(&thermal, &[(CellId(0), &exact)], &[], || Ok(()))
        .unwrap();
    assert!(
        chain
            .with_boundary_laws(&thermal, &[(CellId(0), &exact)], &[], || Ok(()))
            .is_err()
    );
    assert_eq!(thermal.time_s(), 0.0);
    assert_eq!(thermal.cells[&CellId(0)].specific_h_j_kg, 500.0);
}

fn baths(ambient_k: f64, reservoir_k: f64) -> Boundary {
    Boundary {
        ambient_k,
        // Deliberately unrelated initialization slots: selected faces replace
        // these values with their actual material/film conductance.
        htc_w_m2_k: 999.0,
        source_w: 0.0,
        reservoir: Some(Reservoir {
            area_m2: 1.0,
            conductance_w_k: 999.0,
            temperature_k: reservoir_k,
        }),
    }
}

fn ambient(surface_resistance_m2_k_w: f64) -> ConductiveAmbient {
    ConductiveAmbient {
        cell_id: CellId(0),
        half_distance_m: 0.5,
        surface_resistance_m2_k_w,
    }
}

#[test]
fn explicit_two_bath_faces_retain_actual_entropy_events_and_single_ownership() {
    let law = Law {
        nonlinear_k: false,
        plateau: false,
    };
    let (thermal, conduction) = setup(
        &law,
        &[500.0],
        Some(300.0),
        vec![
            TimedEvent {
                time_s: 0.0,
                event: Event::Boundary(CellId(0), baths(700.0, 300.0)),
            },
            TimedEvent {
                time_s: 0.05,
                event: Event::Boundary(CellId(0), baths(300.0, 700.0)),
            },
        ],
    );
    let conduction = conduction
        .with_ambient_faces(&thermal, &[ambient(0.0)], || Ok(()))
        .unwrap();
    let power = thermal
        .power_balance_observed(&conduction, || Ok(()))
        .unwrap();
    near(power.external_power_w, 200.0, 1e-12);
    let step = StepConfig {
        max_dt_s: 0.05,
        energy_tolerance_j: 1e-9,
    };
    let (next, _, heat, entropy) = thermal
        .advance_observed_with_entropy(&conduction, 0.1, step, || Ok(()))
        .unwrap();
    // Native explicit updates: 500 + .05*(2*200-200)=510,
    // then 510 + .05*(2*(-210)+190)=498.5.
    near(next.temperature_k(CellId(0)).unwrap(), 498.5, 1e-12);
    near(heat.external_heat_j, -1.5, 1e-12);
    assert_eq!(heat.substeps, 2);
    near(
        entropy.reservoir_entropy_input_j_k,
        20.0 / 700.0 - 10.0 / 300.0 - 21.0 / 300.0 + 9.5 / 700.0,
        1e-13,
    );
    let production = 20.0 * (1.0 / 500.0 - 1.0 / 700.0)
        - 10.0 * (1.0 / 500.0 - 1.0 / 300.0)
        - 21.0 * (1.0 / 510.0 - 1.0 / 300.0)
        + 9.5 * (1.0 / 510.0 - 1.0 / 700.0);
    near(entropy.conduction_generation_j_k, production, 1e-13);
    assert!(!entropy.has_unclassified_sources);
    assert_eq!(next.cells[&CellId(0)].boundary, baths(300.0, 700.0));
    assert!(matches!(
        thermal.advance_observed(
            &CoupledObservers::new(conduction.clone(), conduction.clone()),
            0.05,
            step,
            || Ok(())
        ),
        Err(ObserverCouplingError::ConflictingAmbientConductance(
            CellId(0)
        ))
    ));
    let mut polls = 0;
    thermal
        .advance_observed(&conduction, 0.05, step, || {
            polls += 1;
            Ok(())
        })
        .unwrap();
    let mut calls = 0;
    assert!(matches!(
        thermal.advance_observed(&conduction, 0.05, step, || {
            calls += 1;
            if calls == polls {
                Err(DepositionError::Cancelled)
            } else {
                Ok(())
            }
        }),
        Err(DepositionError::Cancelled)
    ));
    near(thermal.temperature_k(CellId(0)).unwrap(), 500.0, 0.0);
    assert_eq!(thermal.time_s(), 0.0);
}

// Synthetic phase-sensitive source: conductivity changes on the REAL latent
// plateau even though temperature and dT/dh remain fixed.
struct PhaseTransport {
    conducting: bool,
}
impl ThermalLaw for PhaseTransport {
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
impl ConductivityLaw for PhaseTransport {
    fn conductivity_w_m_k(&self, h: f64) -> Result<f64, DepositionError> {
        self.temperature_k(h)?;
        Ok(if self.conducting { h / 100.0 } else { 0.0 })
    }
}
impl ImplicitConductionLaw for PhaseTransport {
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
        Ok(ImplicitConductionState {
            temperature_k: self.temperature_k(h)?,
            temperature_slope_kg_k_j: if (300.0..=500.0).contains(&h) {
                0.0
            } else {
                1.0
            },
            conductivity_w_m_k: self.conductivity_w_m_k(h)?,
            conductivity_slope_kg_w_m_k_j: if self.conducting { 0.01 } else { 0.0 },
        })
    }
}

#[test]
fn chain_and_graph_two_bath_latent_heat_use_actual_dk_dh_and_exact_insulation() {
    let dt = 0.01;
    for conducting in [true, false] {
        for sign in [1.0, -1.0] {
            let law = PhaseTransport { conducting };
            let boundary = baths(300.0 + sign * 200.0, 300.0 - sign * 100.0);
            let (thermal, conduction) = setup(
                &law,
                &[400.0],
                Some(boundary.reservoir.unwrap().temperature_k),
                vec![TimedEvent {
                    time_s: 0.0,
                    event: Event::Boundary(CellId(0), boundary),
                }],
            );
            let conduction = conduction
                .with_ambient_faces(&thermal, &[ambient(0.0)], || Ok(()))
                .unwrap();
            let chain = ImplicitConductionChain::new(
                &thermal,
                &conduction,
                &[(MaterialId(0), &law)],
                || Ok(()),
            )
            .unwrap();
            let graph = ImplicitConductionGraph::new(
                &thermal,
                &conduction,
                &[(MaterialId(0), &law)],
                ImplicitConductionLinearControls::default(),
                || Ok(()),
            )
            .unwrap();
            let mut workspace = ImplicitConductionWorkspace::new(1).unwrap();
            let endpoint = chain.endpoint(&thermal, &mut || Ok(())).unwrap();
            chain
                .jacobian(&thermal, &endpoint, dt, &mut workspace, &mut || Ok(()))
                .unwrap();
            let rate = if conducting { sign * 3.0 } else { 0.0 };
            near(workspace.diagonal[0], 1.0 - dt * rate, 1e-14);
            // Gambient=2h/100, Greservoir=h/100; T=300 on the plateau.
            // dE/dt = sign*3h, so backward Euler has this closed form.
            let h = 400.0 / (1.0 - dt * rate);
            for result in [
                chain.advance_to(&thermal, dt, controls(dt), &mut workspace, || Ok(())),
                graph.advance_to(&thermal, dt, controls(dt), &mut workspace, || Ok(())),
            ] {
                let (next, report) = result.unwrap();
                near(next.cells[&CellId(0)].specific_h_j_kg, h, 1e-9);
                near(next.temperature_k(CellId(0)).unwrap(), 300.0, 0.0);
                near(report.energy.external_heat_j, h - 400.0, 1e-9);
                let qa = if conducting { dt * sign * 4.0 * h } else { 0.0 };
                let qr = if conducting { -dt * sign * h } else { 0.0 };
                near(
                    report.entropy.reservoir_entropy_input_j_k,
                    qa / boundary.ambient_k + qr / boundary.reservoir.unwrap().temperature_k,
                    1e-11,
                );
                near(
                    report.entropy.thermal_storage_change_j_k,
                    (h - 400.0) / 300.0,
                    1e-11,
                );
                near(report.energy.residual_j, 0.0, 1e-9);
                assert!(!report.entropy.has_unclassified_sources);
                assert_eq!(next.cells[&CellId(0)].boundary, boundary);
            }
            if !conducting {
                let (next, _, heat, entropy) = thermal
                    .advance_observed_with_entropy(
                        &conduction,
                        dt,
                        StepConfig {
                            max_dt_s: dt,
                            energy_tolerance_j: 1e-10,
                        },
                        || Ok(()),
                    )
                    .unwrap();
                near(next.cells[&CellId(0)].specific_h_j_kg, 400.0, 0.0);
                near(heat.external_heat_j, 0.0, 0.0);
                near(entropy.reservoir_entropy_input_j_k, 0.0, 0.0);
            } else {
                let mut limited = controls(dt);
                limited.maximum_residual_evaluations = 1;
                assert!(
                    chain
                        .advance_to(&thermal, dt, limited, &mut workspace, || Ok(()))
                        .is_err()
                );
                near(thermal.cells[&CellId(0)].specific_h_j_kg, 400.0, 0.0);
                assert_eq!(thermal.time_s(), 0.0);
            }
        }
    }
}

#[test]
fn finite_ambient_film_and_material_derivatives_close_nonlinear_cell_energy() {
    let law = Law {
        nonlinear_k: true,
        plateau: false,
    };
    let boundary = baths(600.0, 200.0);
    let (thermal, conduction) = setup(
        &law,
        &[500.0],
        Some(200.0),
        vec![TimedEvent {
            time_s: 0.0,
            event: Event::Boundary(CellId(0), boundary),
        }],
    );
    let conduction = conduction
        .with_ambient_faces(&thermal, &[ambient(2.0)], || Ok(()))
        .unwrap();
    let chain =
        ImplicitConductionChain::new(&thermal, &conduction, &[(MaterialId(0), &law)], || Ok(()))
            .unwrap();
    let dt = 0.01;
    let power = |t: f64| t / (50.0 + 2.0 * t) * (600.0 - t) + t / 100.0 * (200.0 - t);
    let residual = |t| t - 500.0 - dt * power(t);
    let (mut lo, mut hi) = (400.0, 500.0);
    assert!(residual(lo) < 0.0 && residual(hi) > 0.0);
    for _ in 0..80 {
        let mid = (lo + hi) * 0.5;
        if residual(mid) < 0.0 {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    let expected = (lo + hi) * 0.5;
    let mut workspace = ImplicitConductionWorkspace::new(1).unwrap();
    let endpoint = chain.endpoint(&thermal, &mut || Ok(())).unwrap();
    chain
        .jacobian(&thermal, &endpoint, dt, &mut workspace, &mut || Ok(()))
        .unwrap();
    let dg = 50.0 / (1050.0 * 1050.0);
    let derivative = dg * 100.0 - 500.0 / 1050.0 + 0.01 * (-300.0) - 5.0;
    near(workspace.diagonal[0], 1.0 - dt * derivative, 1e-14);
    let (next, report) = chain
        .advance_to(&thermal, dt, controls(dt), &mut workspace, || Ok(()))
        .unwrap();
    near(next.temperature_k(CellId(0)).unwrap(), expected, 1e-9);
    near(report.energy.external_heat_j, dt * power(expected), 1e-9);
    near(report.energy.stored_energy_change_j, expected - 500.0, 1e-9);
    near(report.energy.residual_j, 0.0, 1e-9);
}

#[derive(Clone)]
struct AmbientOnly;
impl DepositionObserver for AmbientOnly {
    type Error = DepositionError;
    fn ambient_conductances(
        &mut self,
        _: &DepositionSimulation<'_>,
        _: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<BTreeMap<CellId, f64>, Self::Error> {
        Ok(BTreeMap::from([(CellId(0), 1.0)]))
    }
}
struct UnsupportedConductivity;
impl ConductivityLaw for UnsupportedConductivity {
    fn conductivity_w_m_k(&self, _: f64) -> Result<f64, DepositionError> {
        Err(DepositionError::MaterialDomain)
    }
}

#[test]
fn ambient_binding_source_coverage_budget_and_cancellation_refuse_without_history() {
    let law = Law {
        nonlinear_k: false,
        plateau: false,
    };
    let (thermal, conduction) = setup(&law, &[300.0], Some(300.0), vec![]);
    let original = thermal.cells.clone();
    for face in [
        ConductiveAmbient {
            half_distance_m: 0.0,
            ..ambient(0.0)
        },
        ConductiveAmbient {
            surface_resistance_m2_k_w: f64::NAN,
            ..ambient(0.0)
        },
        ConductiveAmbient {
            cell_id: CellId(88),
            ..ambient(0.0)
        },
    ] {
        assert!(
            conduction
                .clone()
                .with_ambient_faces(&thermal, &[face], || Ok(()))
                .is_err()
        );
    }
    assert!(
        conduction
            .clone()
            .with_ambient_faces(&thermal, &[ambient(0.0), ambient(0.0)], || Ok(()))
            .is_err()
    );
    assert!(matches!(
        conduction
            .clone()
            .with_ambient_faces(&thermal, &[ambient(0.0)], || Err(
                DepositionError::Cancelled
            )),
        Err(DepositionError::Cancelled)
    ));
    let unsupported = UnsupportedConductivity;
    let empty =
        ConductionObserver::new(&thermal, &[(MaterialId(0), &unsupported)], &[], &[], || {
            Ok(())
        })
        .unwrap();
    assert!(matches!(
        empty.with_ambient_faces(&thermal, &[ambient(0.0)], || Ok(())),
        Err(DepositionError::MaterialDomain)
    ));

    let chain =
        ImplicitConductionChain::new(&thermal, &conduction, &[(MaterialId(0), &law)], || Ok(()))
            .unwrap();
    let mut workspace = ImplicitConductionWorkspace::new(1).unwrap();
    assert!(
        chain
            .advance_observed_to(
                &thermal,
                &AmbientOnly,
                0.1,
                controls(0.1),
                &mut workspace,
                || Ok(())
            )
            .is_err()
    );
    assert_eq!(thermal.cells, original);

    // This observer owns only ambient conduction, so changing an unowned
    // reservoir's coverage must still refuse at the event endpoint itself.
    let mut changed = boundary(Some(300.0));
    changed.reservoir.as_mut().unwrap().area_m2 = 1.5;
    let (thermal, _) = setup(
        &law,
        &[300.0],
        Some(300.0),
        vec![TimedEvent {
            time_s: 0.1,
            event: Event::Boundary(CellId(0), changed),
        }],
    );
    let only_ambient =
        ConductionObserver::new(&thermal, &[(MaterialId(0), &law)], &[], &[], || Ok(()))
            .unwrap()
            .with_ambient_faces(&thermal, &[ambient(0.0)], || Ok(()))
            .unwrap();
    assert!(
        thermal
            .advance_observed(
                &only_ambient,
                0.1,
                StepConfig {
                    max_dt_s: 0.1,
                    energy_tolerance_j: 1e-9
                },
                || Ok(())
            )
            .is_err()
    );
    near(thermal.cells[&CellId(0)].specific_h_j_kg, 300.0, 0.0);
    assert_eq!(thermal.time_s(), 0.0);
}
