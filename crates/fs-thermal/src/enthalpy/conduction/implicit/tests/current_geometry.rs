//! Frozen mechanical metrics share the native implicit heat and entropy owners.
use super::super::boundary::{
    BathTemperatureKnot, BathTemperatureSchedule, GrayAmbientBoundary, GrayAmbientControls,
    GrayTwoBathBoundary,
};
use super::super::link::{GrayRadiativeControls, GrayRadiativeLink};
use super::*;

fn current_pair<'a>(
    law: &'a Law,
) -> (
    DepositionSimulation<'a>,
    ConductionObserver<'a>,
    ConductionObserver<'a>,
) {
    let (mut thermal, _) = setup(law, &[600.0, 300.0], None, vec![]);
    thermal.cells.get_mut(&CellId(0)).unwrap().mass_kg = 2.0;
    thermal.cells.get_mut(&CellId(1)).unwrap().mass_kg = 3.0;
    let original = ConductionObserver::new(
        &thermal,
        &[(MaterialId(0), law)],
        &[ConductiveLink {
            link_id: 0,
            half_distances_m: [0.5; 2],
            contact_resistance_m2_k_w: 0.125,
        }],
        &[],
        || Ok(()),
    )
    .unwrap();
    let current = original
        .with_current_geometry(
            &thermal,
            CurrentConductionGeometry {
                links: [(
                    0,
                    CurrentConductiveLink {
                        area_m2: 2.5,
                        half_distances_m: [0.1, 0.4],
                    },
                )]
                .into(),
                ..Default::default()
            },
            || Ok(()),
        )
        .unwrap();
    (thermal, original, current)
}

fn unchanged_inventory(initial: &DepositionSimulation<'_>, next: &DepositionSimulation<'_>) {
    for (&id, old) in initial.cells() {
        let new = &next.cells()[&id];
        assert_eq!(new.mass_kg, old.mass_kg);
        assert_eq!(new.geometry, old.geometry);
        assert_eq!(new.material_id, old.material_id);
        assert_eq!(new.body_id, old.body_id);
        assert_eq!(new.segment_id, old.segment_id);
        assert_eq!(new.boundary, old.boundary);
    }
    assert_eq!(next.links(), initial.links());
}

fn check_jacobian(
    chain: &ImplicitConductionChain<'_>,
    thermal: &DepositionSimulation<'_>,
    energy: &[f64],
    dt: f64,
) {
    let step = BackwardEulerStep::new(thermal, dt, || Ok(())).unwrap();
    let mut report = empty_report();
    let trial = chain
        .evaluate(
            &step,
            energy,
            &BTreeMap::new(),
            &[],
            false,
            controls(dt),
            &mut report,
            &mut || Ok(()),
        )
        .unwrap();
    let mut workspace = ImplicitConductionWorkspace::new(energy.len()).unwrap();
    chain
        .jacobian(thermal, trial.payload(), dt, &mut workspace, &mut || Ok(()))
        .unwrap();
    for column in 0..energy.len() {
        let mut plus = energy.to_vec();
        let mut minus = energy.to_vec();
        plus[column] += 1e-3;
        minus[column] -= 1e-3;
        let a = chain
            .evaluate(
                &step,
                &plus,
                &BTreeMap::new(),
                &[],
                false,
                controls(dt),
                &mut report,
                &mut || Ok(()),
            )
            .unwrap();
        let b = chain
            .evaluate(
                &step,
                &minus,
                &BTreeMap::new(),
                &[],
                false,
                controls(dt),
                &mut report,
                &mut || Ok(()),
            )
            .unwrap();
        for row in 0..energy.len() {
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
                2e-8,
            );
        }
    }
}

#[test]
fn current_faces_drive_implicit_two_body_relaxation_and_original_inventory() {
    let law = Law {
        nonlinear_k: false,
        plateau: false,
    };
    let (thermal, original, current) = current_pair(&law);
    let chain =
        ImplicitConductionChain::new(&thermal, &current, &[(MaterialId(0), &law)], || Ok(()))
            .unwrap();
    let graph = ImplicitConductionGraph::new(
        &thermal,
        &current,
        &[(MaterialId(0), &law)],
        ImplicitConductionLinearControls::default(),
        || Ok(()),
    )
    .unwrap();
    let mut workspace = ImplicitConductionWorkspace::new(2).unwrap();
    // Physical current G=2.5/(0.1+0.125+0.4)=4 W/K, capacities 2 and 3 J/K.
    // The weighted mean remains 420 K; only the temperature difference decays.
    let conductance = 4.0;
    let decay_rate = conductance * (0.5 + 1.0 / 3.0);
    let gamma = 1.0 - core::f64::consts::FRAC_1_SQRT_2;
    for second_order in [false, true] {
        let dt = if second_order { 0.3 } else { 2.0 };
        let z = dt * decay_rate;
        let amplification = if second_order {
            (1.0 - (1.0 - 2.0 * gamma) * z) / (1.0 + gamma * z).powi(2)
        } else {
            1.0 / (1.0 + z)
        };
        let expected = [420.0 + 180.0 * amplification, 420.0 - 120.0 * amplification];
        let production = |t: [f64; 2]| conductance * (t[0] - t[1]).powi(2) / (t[0] * t[1]);
        let generation = if second_order {
            let first = [
                420.0 + 180.0 / (1.0 + gamma * z),
                420.0 - 120.0 / (1.0 + gamma * z),
            ];
            dt * ((1.0 - gamma) * production(first) + gamma * production(expected))
        } else {
            dt * production(expected)
        };
        for sparse in [false, true] {
            let (next, report) = match (sparse, second_order) {
                (false, false) => {
                    chain.advance_to(&thermal, dt, controls(dt), &mut workspace, || Ok(()))
                }
                (true, false) => {
                    graph.advance_to(&thermal, dt, controls(dt), &mut workspace, || Ok(()))
                }
                (false, true) => {
                    chain.advance_sdirk2_to(&thermal, dt, controls(dt), &mut workspace, || Ok(()))
                }
                (true, true) => {
                    graph.advance_sdirk2_to(&thermal, dt, controls(dt), &mut workspace, || Ok(()))
                }
            }
            .unwrap();
            for (i, temperature) in expected.into_iter().enumerate() {
                near(
                    next.temperature_k(CellId(i as u64)).unwrap(),
                    temperature,
                    1e-10,
                );
            }
            near(report.energy.external_heat_j, 0.0, 0.0);
            near(report.energy.stored_energy_change_j, 0.0, 1e-9);
            near(report.entropy.conduction_generation_j_k, generation, 1e-11);
            near(
                report.entropy.thermal_storage_change_j_k,
                2.0 * (expected[0] / 600.0).ln() + 3.0 * (expected[1] / 300.0).ln(),
                1e-12,
            );
            assert_eq!(report.energy.substeps, 1);
            assert!(report.maximum_cell_residual_j <= 1e-10);
            unchanged_inventory(&thermal, &next);
        }
    }
    // The same refreshed conductance needs several native explicit CFL steps.
    let (_, _, explicit) = thermal
        .advance_observed(
            &current,
            2.0,
            StepConfig {
                max_dt_s: 2.0,
                energy_tolerance_j: 1e-9,
            },
            || Ok(()),
        )
        .unwrap();
    assert!(explicit.substeps > 1);
    assert!(original.current_geometry().is_none());
    assert_eq!(thermal.time_s(), 0.0);
    assert_eq!(thermal.temperature_k(CellId(0)).unwrap(), 600.0);
}

#[test]
fn current_graph_faces_retain_asymmetric_reference_coverage_and_stiff_balance() {
    let law = Law {
        nonlinear_k: false,
        plateau: false,
    };
    let (mut thermal, _) = setup(&law, &[300.0, 600.0, 700.0], None, vec![]);
    for (i, mass) in [2.0, 1.0, 3.0].into_iter().enumerate() {
        thermal.cells.get_mut(&CellId(i as u64)).unwrap().mass_kg = mass;
    }
    // A branched hub, with distinct reference footprints on both sides.
    thermal.links.get_mut(&1).unwrap().a = CellId(0);
    thermal.covered_areas.insert(0, [0.4, 0.7]);
    thermal.covered_areas.insert(1, [0.6, 0.9]);
    let original = ConductionObserver::new(
        &thermal,
        &[(MaterialId(0), &law)],
        &[
            ConductiveLink {
                link_id: 0,
                half_distances_m: [0.5; 2],
                contact_resistance_m2_k_w: 0.5,
            },
            ConductiveLink {
                link_id: 1,
                half_distances_m: [0.5; 2],
                contact_resistance_m2_k_w: 1.0,
            },
        ],
        &[],
        || Ok(()),
    )
    .unwrap();
    let current = original
        .with_current_geometry(
            &thermal,
            CurrentConductionGeometry {
                links: [
                    (
                        0,
                        CurrentConductiveLink {
                            area_m2: 4.0,
                            half_distances_m: [0.125, 0.375],
                        },
                    ),
                    (
                        1,
                        CurrentConductiveLink {
                            area_m2: 3.0,
                            half_distances_m: [0.4, 0.6],
                        },
                    ),
                ]
                .into(),
                ..Default::default()
            },
            || Ok(()),
        )
        .unwrap();
    let graph = ImplicitConductionGraph::new(
        &thermal,
        &current,
        &[(MaterialId(0), &law)],
        ImplicitConductionLinearControls::default(),
        || Ok(()),
    )
    .unwrap();
    let dt = 10.0;
    // Independently eliminate each leaf from m_i(T_i-T_i0)=dt*G_i(T_hub-T_i).
    let g = [4.0, 1.5];
    let weights = [
        dt * g[0] / (1.0 + dt * g[0]),
        dt * g[1] * 3.0 / (3.0 + dt * g[1]),
    ];
    let hub =
        (2.0 * 300.0 + weights[0] * 600.0 + weights[1] * 700.0) / (2.0 + weights[0] + weights[1]);
    let expected = [
        hub,
        (600.0 + dt * g[0] * hub) / (1.0 + dt * g[0]),
        (3.0 * 700.0 + dt * g[1] * hub) / (3.0 + dt * g[1]),
    ];
    let mut workspace = ImplicitConductionWorkspace::new(3).unwrap();
    let (next, report) = graph
        .advance_to(&thermal, dt, controls(dt), &mut workspace, || Ok(()))
        .unwrap();
    for (i, temperature) in expected.into_iter().enumerate() {
        near(
            next.temperature_k(CellId(i as u64)).unwrap(),
            temperature,
            2e-10,
        );
    }
    assert_eq!(report.energy.substeps, 1);
    near(report.energy.external_heat_j, 0.0, 0.0);
    near(report.energy.stored_energy_change_j, 0.0, 1e-9);
    assert_eq!(next.covered_areas, thermal.covered_areas);
    unchanged_inventory(&thermal, &next);

    let mut stale = thermal.clone();
    stale.covered_areas.get_mut(&0).unwrap()[0] += 0.01;
    assert!(
        graph
            .advance_to(&stale, dt, controls(dt), &mut workspace, || Ok(()))
            .is_err()
    );
    assert!(
        ImplicitConductionGraph::new(
            &stale,
            &current,
            &[(MaterialId(0), &law)],
            ImplicitConductionLinearControls::default(),
            || Ok(()),
        )
        .is_err()
    );
    let mut changed_mass = thermal.clone();
    changed_mass.cells.get_mut(&CellId(0)).unwrap().mass_kg *= 2.0;
    assert!(
        ImplicitConductionGraph::new(
            &changed_mass,
            &current,
            &[(MaterialId(0), &law)],
            ImplicitConductionLinearControls::default(),
            || Ok(()),
        )
        .is_err()
    );
}

struct Resistance<'a>(&'a Law);
impl ImplicitContactLaw for Resistance<'_> {
    fn caloric_laws(&self) -> [&dyn ThermalLaw; 2] {
        [self.0; 2]
    }
    fn enthalpy_bounds_j_kg(&self) -> [[f64; 2]; 2] {
        [[100.0, 1000.0]; 2]
    }
    fn state_at_specific_enthalpies(
        &self,
        h: [f64; 2],
    ) -> Result<ImplicitContactState, DepositionError> {
        Ok(ImplicitContactState {
            resistance_m2_k_w: 0.001 * h[0] + 0.0005 * h[1],
            resistance_enthalpy_derivatives_m2_kg_k_w_j: [0.001, 0.0005],
        })
    }
}

#[test]
fn current_metrics_enter_nonlinear_contact_bath_and_scheduled_bath_jacobians() {
    let law = Law {
        nonlinear_k: true,
        plateau: false,
    };
    let resistance = Resistance(&law);
    let (mut thermal, _) = setup(&law, &[550.0, 400.0, 300.0], Some(700.0), vec![]);
    for (i, mass) in [2.0, 3.0, 0.5].into_iter().enumerate() {
        thermal.cells.get_mut(&CellId(i as u64)).unwrap().mass_kg = mass;
    }
    thermal
        .cells
        .get_mut(&CellId(2))
        .unwrap()
        .boundary
        .ambient_k = 250.0;
    let current = ConductionObserver::new(
        &thermal,
        &[(MaterialId(0), &law)],
        &[
            ConductiveLink {
                link_id: 0,
                half_distances_m: [0.5; 2],
                contact_resistance_m2_k_w: 0.0,
            },
            ConductiveLink {
                link_id: 1,
                half_distances_m: [0.5; 2],
                contact_resistance_m2_k_w: 0.0,
            },
        ],
        &[ConductiveReservoir {
            cell_id: CellId(0),
            half_distance_m: 0.5,
            surface_resistance_m2_k_w: 0.2,
        }],
        || Ok(()),
    )
    .unwrap()
    .with_ambient_faces(
        &thermal,
        &[ConductiveAmbient {
            cell_id: CellId(2),
            half_distance_m: 1.0,
            surface_resistance_m2_k_w: 0.1,
        }],
        || Ok(()),
    )
    .unwrap()
    .with_current_geometry(
        &thermal,
        CurrentConductionGeometry {
            links: [
                (
                    0,
                    CurrentConductiveLink {
                        area_m2: 3.0,
                        half_distances_m: [0.15, 0.9],
                    },
                ),
                (
                    1,
                    CurrentConductiveLink {
                        area_m2: 0.8,
                        half_distances_m: [0.6, 0.1],
                    },
                ),
            ]
            .into(),
            reservoirs: [(
                CellId(0),
                CurrentConductiveBoundary {
                    area_m2: 1.8,
                    half_distance_m: 0.2,
                },
            )]
            .into(),
            ambient_faces: [(
                CellId(2),
                CurrentConductiveBoundary {
                    area_m2: 0.4,
                    half_distance_m: 1.3,
                },
            )]
            .into(),
        },
        || Ok(()),
    )
    .unwrap();
    let chain =
        ImplicitConductionChain::new(&thermal, &current, &[(MaterialId(0), &law)], || Ok(()))
            .unwrap()
            .with_contacts(&thermal, &[(0, &resistance)], || Ok(()))
            .unwrap();
    let endpoint = chain.endpoint(&thermal, &mut || Ok(())).unwrap();
    near(
        endpoint.link_values[&0],
        3.0 / (0.15 / 5.5 + 0.75 + 0.9 / 4.0),
        1e-14,
    );
    near(
        endpoint.link_values[&1],
        0.8 / (0.6 / 4.0 + 0.1 / 3.0),
        1e-14,
    );
    near(
        endpoint.boundary_values.reservoir_w_k[&CellId(0)],
        1.8 / (0.2 / 5.5 + 0.2),
        1e-14,
    );
    near(
        endpoint.boundary_values.ambient_w_k[&CellId(2)],
        0.4 / (1.3 / 3.0 + 0.1),
        1e-14,
    );
    check_jacobian(&chain, &thermal, &[-15.0, 8.0, 6.0], 0.23);

    let history = |a, b| {
        BathTemperatureSchedule::new(
            vec![
                BathTemperatureKnot {
                    time_s: 0.0,
                    temperature_k: a,
                },
                BathTemperatureKnot {
                    time_s: 1.0,
                    temperature_k: b,
                },
            ],
            2,
            || Ok(()),
        )
        .unwrap()
    };
    let hot = history(700.0, 650.0);
    let cold = history(250.0, 400.0);
    let scheduled = chain
        .with_bath_temperature_histories(
            &thermal,
            &[(CellId(0), &hot)],
            &[(CellId(2), &cold)],
            || Ok(()),
        )
        .unwrap();
    check_jacobian(&scheduled, &thermal, &[-15.0, 8.0, 6.0], 0.23);
}

#[test]
fn current_metrics_preserve_the_finite_insulating_onset_tangent() {
    let law = ConductivityOnset;
    let (thermal, observer) = setup(&law, &[300.0, 400.0], None, vec![]);
    let current = observer
        .with_current_geometry(
            &thermal,
            CurrentConductionGeometry {
                links: [(
                    0,
                    CurrentConductiveLink {
                        area_m2: 2.5,
                        half_distances_m: [0.1, 0.4],
                    },
                )]
                .into(),
                ..Default::default()
            },
            || Ok(()),
        )
        .unwrap();
    let chain =
        ImplicitConductionChain::new(&thermal, &current, &[(MaterialId(0), &law)], || Ok(()))
            .unwrap();
    let endpoint = chain.endpoint(&thermal, &mut || Ok(())).unwrap();
    let mut workspace = ImplicitConductionWorkspace::new(2).unwrap();
    chain
        .jacobian(&thermal, &endpoint, 0.01, &mut workspace, &mut || Ok(()))
        .unwrap();
    // G=0 but dG/dh_a=A*k'_a/d_a=0.25. The 100 K difference gives dQ/dh_a=25.
    near(endpoint.link_values[&0], 0.0, 0.0);
    near(workspace.diagonal[0], 0.75, 1e-15);
    near(workspace.lower[0], 0.25, 1e-15);
    near(workspace.upper[0], 0.0, 0.0);
    near(workspace.diagonal[1], 1.0, 0.0);
}

// The legacy fixture deliberately retains the default reference-geometry
// contract. This wrapper explicitly declares that its exact integral consumes
// the current area/distances passed to every flux and tangent query.
struct CurrentIntegratedLink<'a>(IntegratedLinearK<'a>);
impl ImplicitLinkLaw for CurrentIntegratedLink<'_> {
    fn caloric_laws(&self) -> [&dyn ThermalLaw; 2] {
        self.0.caloric_laws()
    }
    fn supports_current_geometry(&self) -> bool {
        true
    }
    fn state_at_specific_enthalpies(
        &self,
        h: [f64; 2],
        geometry: ConductiveLink,
        patch: ContactPatch,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<ImplicitLinkState, DepositionError> {
        self.0
            .state_at_specific_enthalpies(h, geometry, patch, checkpoint)
    }
}

#[test]
fn source_integrated_current_link_matches_kirchhoff_relaxation() {
    let law = Law {
        nonlinear_k: true,
        plateau: false,
    };
    let (thermal, observer) = setup(&law, &[400.0, 200.0], None, vec![]);
    let current = observer
        .with_current_geometry(
            &thermal,
            CurrentConductionGeometry {
                links: [(
                    0,
                    CurrentConductiveLink {
                        area_m2: 2.5,
                        half_distances_m: [0.125, 0.375],
                    },
                )]
                .into(),
                ..Default::default()
            },
            || Ok(()),
        )
        .unwrap();
    let pinned = IntegratedLinearK(&law);
    assert!(
        ImplicitConductionChain::new(&thermal, &current, &[(MaterialId(0), &law)], || Ok(()),)
            .unwrap()
            .with_link_laws(&thermal, &[(0, &pinned)], || Ok(()))
            .is_err()
    );
    let integrated = CurrentIntegratedLink(IntegratedLinearK(&law));
    let chain =
        ImplicitConductionChain::new(&thermal, &current, &[(MaterialId(0), &law)], || Ok(()))
            .unwrap()
            .with_link_laws(&thermal, &[(0, &integrated)], || Ok(()))
            .unwrap();
    let graph = ImplicitConductionGraph::new(
        &thermal,
        &current,
        &[(MaterialId(0), &law)],
        ImplicitConductionLinearControls::default(),
        || Ok(()),
    )
    .unwrap()
    .with_link_laws(&thermal, &[(0, &integrated)], || Ok(()))
    .unwrap();
    // Integral k(T)=T/100: G=A*(Ta+Tb)/(200*(da+db))=15 W/K.
    // Conservation fixes Ta+Tb=600 K, so one BE interval has an exact linear mode.
    near(
        chain
            .endpoint(&thermal, &mut || Ok(()))
            .unwrap()
            .link_values[&0],
        15.0,
        1e-14,
    );
    check_jacobian(&chain, &thermal, &[-9.0, 6.0], 0.5);
    let mut workspace = ImplicitConductionWorkspace::new(2).unwrap();
    for result in [
        chain.advance_to(&thermal, 0.5, controls(0.5), &mut workspace, || Ok(())),
        graph.advance_to(&thermal, 0.5, controls(0.5), &mut workspace, || Ok(())),
    ] {
        let (next, report) = result.unwrap();
        near(next.temperature_k(CellId(0)).unwrap(), 306.25, 1e-10);
        near(next.temperature_k(CellId(1)).unwrap(), 293.75, 1e-10);
        assert_eq!(report.energy.substeps, 1);
        near(report.energy.external_heat_j, 0.0, 0.0);
        near(report.energy.stored_energy_change_j, 0.0, 1e-9);
        unchanged_inventory(&thermal, &next);
    }
    // Native gray links opt in, but retain their equal-area/full-view domain.
    // Current half-distances can change independently of the original coverage;
    // an unsupported change of radiating area is still refused at the source.
    let gray = GrayRadiativeLink::new(
        [&law, &law],
        [0.0; 2],
        2.0,
        [0.0; 2],
        GrayRadiativeControls::default(),
    )
    .unwrap();
    assert!(
        ImplicitConductionChain::new(&thermal, &current, &[(MaterialId(0), &law)], || Ok(()),)
            .unwrap()
            .with_link_laws(&thermal, &[(0, &gray)], || Ok(()))
            .is_err()
    );
    let mut metrics = current.current_geometry().unwrap().clone();
    metrics.links.get_mut(&0).unwrap().area_m2 = 1.0;
    let admitted = observer
        .with_current_geometry(&thermal, metrics, || Ok(()))
        .unwrap();
    let gray_chain =
        ImplicitConductionChain::new(&thermal, &admitted, &[(MaterialId(0), &law)], || Ok(()))
            .unwrap()
            .with_link_laws(&thermal, &[(0, &gray)], || Ok(()))
            .unwrap();
    near(
        gray_chain
            .endpoint(&thermal, &mut || Ok(()))
            .unwrap()
            .link_values[&0],
        1.0 / (0.125 / 4.0 + 0.5 + 0.375 / 2.0),
        1e-12,
    );
    check_jacobian(&gray_chain, &thermal, &[-9.0, 6.0], 0.5);
}

struct IntegratedBoundary<'a>(&'a Law);
impl ImplicitBoundaryLaw for IntegratedBoundary<'_> {
    fn supports_current_geometry(&self) -> bool {
        true
    }

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
            return Err(DepositionError::InvalidInput(
                "the exact test boundary requires perfect surface contact",
            ));
        }
        let state = self.0.state_at_specific_enthalpy(h)?;
        let scale = geometry.area_m2 / geometry.half_distance_m;
        Ok(ImplicitBoundaryState {
            conductance_w_k: scale * (state.temperature_k + bath) / 200.0,
            power_enthalpy_derivative_w_kg_j: scale
                * state.temperature_k
                * state.temperature_slope_kg_k_j
                / 100.0,
        })
    }
}

struct PinnedBoundary<'a>(&'a Law);
impl ImplicitBoundaryLaw for PinnedBoundary<'_> {
    fn caloric_law(&self) -> &dyn ThermalLaw {
        self.0
    }
    fn state_at_specific_enthalpy(
        &self,
        _: f64,
        _: f64,
        _: ImplicitBoundaryGeometry,
        _: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<ImplicitBoundaryState, DepositionError> {
        panic!("a reference-pinned boundary must refuse before its flux query")
    }
}
impl ImplicitBoundaryExchangeLaw for PinnedBoundary<'_> {
    fn caloric_law(&self) -> &dyn ThermalLaw {
        self.0
    }
    fn state_at_specific_enthalpy(
        &self,
        _: f64,
        _: f64,
        _: ImplicitBoundaryGeometry,
        _: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<ImplicitBoundaryExchangeState, DepositionError> {
        panic!("a reference-pinned exchange must refuse before its flux query")
    }
}

#[test]
fn source_integrated_current_baths_retain_separate_heat_and_entropy() {
    let law = Law {
        nonlinear_k: true,
        plateau: false,
    };
    let (mut thermal, _) = setup(&law, &[500.0], Some(300.0), vec![]);
    let cell = thermal.cells.get_mut(&CellId(0)).unwrap();
    cell.mass_kg = 2.5;
    cell.boundary.ambient_k = 700.0;
    let current = ConductionObserver::new(
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
    .with_ambient_faces(
        &thermal,
        &[ConductiveAmbient {
            cell_id: CellId(0),
            half_distance_m: 0.5,
            surface_resistance_m2_k_w: 0.0,
        }],
        || Ok(()),
    )
    .unwrap()
    .with_current_geometry(
        &thermal,
        CurrentConductionGeometry {
            links: BTreeMap::new(),
            reservoirs: [(
                CellId(0),
                CurrentConductiveBoundary {
                    area_m2: 1.5,
                    half_distance_m: 0.25,
                },
            )]
            .into(),
            ambient_faces: [(
                CellId(0),
                CurrentConductiveBoundary {
                    area_m2: 0.4,
                    half_distance_m: 0.8,
                },
            )]
            .into(),
        },
        || Ok(()),
    )
    .unwrap();
    let fresh = || {
        ImplicitConductionChain::new(&thermal, &current, &[(MaterialId(0), &law)], || Ok(()))
            .unwrap()
    };
    let pinned = PinnedBoundary(&law);
    assert!(
        fresh()
            .with_boundary_laws(&thermal, &[(CellId(0), &pinned)], &[], || Ok(()),)
            .is_err()
    );
    assert!(
        fresh()
            .with_boundary_exchanges(&thermal, &[], &[(CellId(0), &pinned)], || Ok(()),)
            .is_err()
    );
    let integrated = IntegratedBoundary(&law);
    let chain =
        ImplicitConductionChain::new(&thermal, &current, &[(MaterialId(0), &law)], || Ok(()))
            .unwrap()
            .with_boundary_laws(
                &thermal,
                &[(CellId(0), &integrated)],
                &[(CellId(0), &integrated)],
                || Ok(()),
            )
            .unwrap();
    let graph = ImplicitConductionGraph::new(
        &thermal,
        &current,
        &[(MaterialId(0), &law)],
        ImplicitConductionLinearControls::default(),
        || Ok(()),
    )
    .unwrap()
    .with_boundary_laws(
        &thermal,
        &[(CellId(0), &integrated)],
        &[(CellId(0), &integrated)],
        || Ok(()),
    )
    .unwrap();
    let dt = 0.25;
    let mass = 2.5;
    // Kirchhoff integration yields independent Q_in=(A/d)*(Tb²-T²)/200.
    let metric = [6.0, 0.5];
    let rhs =
        mass * 500.0 + dt * (metric[0] * 300.0_f64.powi(2) + metric[1] * 700.0_f64.powi(2)) / 200.0;
    let quadratic = dt * (metric[0] + metric[1]) / 200.0;
    let expected = 2.0 * rhs / (mass + (mass * mass + 4.0 * quadratic * rhs).sqrt());
    check_jacobian(&chain, &thermal, &[-30.0], dt);
    let mut workspace = ImplicitConductionWorkspace::new(1).unwrap();
    for result in [
        chain.advance_to(&thermal, dt, controls(dt), &mut workspace, || Ok(())),
        graph.advance_to(&thermal, dt, controls(dt), &mut workspace, || Ok(())),
    ] {
        let (next, report) = result.unwrap();
        let reservoir_heat = dt * metric[0] * (300.0_f64.powi(2) - expected * expected) / 200.0;
        let ambient_heat = dt * metric[1] * (700.0_f64.powi(2) - expected * expected) / 200.0;
        near(next.temperature_k(CellId(0)).unwrap(), expected, 1e-10);
        near(
            report.energy.external_heat_j,
            reservoir_heat + ambient_heat,
            1e-9,
        );
        near(
            report.energy.stored_energy_change_j,
            mass * (expected - 500.0),
            1e-9,
        );
        near(
            report.entropy.reservoir_entropy_input_j_k,
            reservoir_heat / 300.0 + ambient_heat / 700.0,
            1e-11,
        );
        assert_eq!(report.energy.substeps, 1);
        unchanged_inventory(&thermal, &next);
    }
    // Native single- and two-bath gray closures consume those same measures.
    // Zero emissivity gives an independent exact bulk-plus-convection limit.
    let gray = GrayAmbientBoundary::new(&law, 0.0, 2.0, GrayAmbientControls::default()).unwrap();
    let split =
        GrayTwoBathBoundary::new(&law, 0.0, 2.0, 800.0, GrayAmbientControls::default()).unwrap();
    let gray_chain = fresh()
        .with_boundary_laws(&thermal, &[], &[(CellId(0), &gray)], || Ok(()))
        .unwrap()
        .with_boundary_exchanges(&thermal, &[(CellId(0), &split)], &[], || Ok(()))
        .unwrap();
    let endpoint = gray_chain.endpoint(&thermal, &mut || Ok(())).unwrap();
    near(
        endpoint.ambient_states[&CellId(0)].conductance_w_k,
        0.4 / (0.8 / 5.0 + 0.5),
        1e-12,
    );
    let exchange = endpoint.reservoir_exchange_states[&CellId(0)];
    near(
        exchange.baths[0].outward_power_w,
        1.5 * (500.0 - 300.0) / (0.25 / 5.0 + 0.5),
        1e-10,
    );
    near(exchange.baths[1].outward_power_w, 0.0, 0.0);
    check_jacobian(&gray_chain, &thermal, &[-30.0], dt);
}

#[test]
fn current_geometry_late_cancellation_retries_the_same_sdirk_interval() {
    let law = Law {
        nonlinear_k: false,
        plateau: false,
    };
    let (thermal, original, current) = current_pair(&law);
    let chain =
        ImplicitConductionChain::new(&thermal, &current, &[(MaterialId(0), &law)], || Ok(()))
            .unwrap();
    let mut workspace = ImplicitConductionWorkspace::new(2).unwrap();
    let mut polls = 0;
    let (expected, _) = chain
        .advance_sdirk2_to(&thermal, 0.6, controls(0.2), &mut workspace, || {
            polls += 1;
            Ok(())
        })
        .unwrap();
    let mut calls = 0;
    assert!(
        chain
            .advance_sdirk2_to(&thermal, 0.6, controls(0.2), &mut workspace, || {
                calls += 1;
                if calls == polls {
                    Err(DepositionError::Cancelled)
                } else {
                    Ok(())
                }
            },)
            .is_err()
    );
    assert_eq!(calls, polls);
    let (retry, _) = chain
        .advance_sdirk2_to(&thermal, 0.6, controls(0.2), &mut workspace, || Ok(()))
        .unwrap();
    assert_eq!(retry.cells(), expected.cells());
    assert_eq!(thermal.time_s(), 0.0);
    assert_eq!(thermal.temperature_k(CellId(0)).unwrap(), 600.0);
    assert!(original.current_geometry().is_none());
    unchanged_inventory(&thermal, &retry);
}
