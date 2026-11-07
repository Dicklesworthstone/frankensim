//! Conserved inventory and original density in the nonlinear heat geometry.
use super::tests::near;
use super::*;

#[derive(Clone, Copy)]
struct Expanding {
    variable_density: bool,
    nonlinear_k: bool,
    density_tangent: bool,
}
static EXPANDING: Expanding = Expanding {
    variable_density: true,
    nonlinear_k: false,
    density_tangent: true,
};
impl ThermalLaw for Expanding {
    fn temperature_k(&self, h: f64) -> Result<f64, DepositionError> {
        if !(100.0..=1000.0).contains(&h) {
            return Err(DepositionError::MaterialDomain);
        }
        Ok(h)
    }
    fn maximum_temperature_slope(&self) -> f64 {
        1.0
    }
    fn bulk_density_kg_m3(&self, h: f64) -> Result<Option<f64>, DepositionError> {
        let temperature = self.temperature_k(h)?;
        Ok(Some(if self.variable_density {
            300.0 / temperature
        } else {
            1.0
        }))
    }
    fn specific_entropy_change_j_kg_k(
        &self,
        a: f64,
        b: f64,
    ) -> Result<Option<f64>, DepositionError> {
        Ok(Some((self.temperature_k(b)? / self.temperature_k(a)?).ln()))
    }
}
impl ConductivityLaw for Expanding {
    fn conductivity_w_m_k(&self, h: f64) -> Result<f64, DepositionError> {
        let temperature = self.temperature_k(h)?;
        Ok(if self.nonlinear_k {
            1.0 + temperature / 300.0
        } else {
            1.0
        })
    }
}
impl ImplicitConductionLaw for Expanding {
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
            temperature_slope_kg_k_j: 1.0,
            conductivity_w_m_k: self.conductivity_w_m_k(h)?,
            conductivity_slope_kg_w_m_k_j: if self.nonlinear_k { 1.0 / 300.0 } else { 0.0 },
        })
    }
    fn specific_volume_derivative_at_specific_enthalpy(
        &self,
        h: f64,
    ) -> Result<Option<f64>, DepositionError> {
        self.temperature_k(h)?;
        Ok(self.density_tangent.then_some(if self.variable_density {
            1.0 / 300.0
        } else {
            0.0
        }))
    }
}

fn setup<'a>(
    law: &'a Expanding,
    enthalpies: &[f64],
    baths: Option<(f64, Option<f64>)>,
) -> (DepositionSimulation<'a>, ConductionObserver<'a>) {
    let mut events = Vec::new();
    let mut distances = Vec::new();
    for (i, &h) in enthalpies.iter().enumerate() {
        let mass_kg = (i + 1) as f64;
        let volume_m3 = mass_kg / law.bulk_density_kg_m3(h).unwrap().unwrap();
        distances.push(0.5 * volume_m3);
        events.push(TimedEvent {
            time_s: 0.0,
            event: Event::Birth(DepositedCell {
                id: CellId(i as u64),
                material_id: MaterialId(0),
                body_id: 0,
                segment_id: 0,
                mass_kg,
                specific_h_j_kg: h,
                geometry: CellGeometry {
                    start_m: [i as f64, 0.0, 0.0],
                    end_m: [i as f64 + 1.0, 0.0, 0.0],
                    volume_m3,
                    surface_area_m2: 2.0,
                },
                boundary: Boundary {
                    ambient_k: baths.and_then(|b| b.1).unwrap_or(300.0),
                    htc_w_m2_k: 0.0,
                    source_w: 0.0,
                    reservoir: if i == 0 {
                        baths.map(|(temperature_k, ambient)| Reservoir {
                            area_m2: if ambient.is_some() && enthalpies.len() > 1 {
                                0.5
                            } else {
                                1.0
                            },
                            conductance_w_k: 0.0,
                            temperature_k,
                        })
                    } else {
                        None
                    },
                },
            }),
        });
    }
    let mut links = Vec::new();
    for i in 0..enthalpies.len() - 1 {
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
            half_distances_m: [distances[i], distances[i + 1]],
            contact_resistance_m2_k_w: if law.nonlinear_k { 0.2 } else { 0.0 },
        });
    }
    let thermal = DepositionSimulation::new(
        events,
        &[(MaterialId(0), law)],
        DepositionBudget {
            max_cells: 8,
            max_links: 8,
            max_events: 32,
            max_substeps: 1000,
        },
    )
    .unwrap()
    .advance_to(
        0.0,
        StepConfig {
            max_dt_s: 1.0,
            energy_tolerance_j: 1e-9,
        },
        || Ok(()),
    )
    .unwrap()
    .0;
    let reservoirs = baths.map(|_| ConductiveReservoir {
        cell_id: CellId(0),
        half_distance_m: distances[0],
        surface_resistance_m2_k_w: if law.nonlinear_k { 0.3 } else { 0.0 },
    });
    let mut observer = ConductionObserver::new(
        &thermal,
        &[(MaterialId(0), law)],
        &links,
        &reservoirs.into_iter().collect::<Vec<_>>(),
        || Ok(()),
    )
    .unwrap();
    if baths.and_then(|b| b.1).is_some() {
        observer = observer
            .with_ambient_faces(
                &thermal,
                &[ConductiveAmbient {
                    cell_id: CellId(0),
                    half_distance_m: distances[0],
                    surface_resistance_m2_k_w: 0.4,
                }],
                || Ok(()),
            )
            .unwrap();
    }
    (thermal, observer)
}

fn layers<'a>(
    thermal: &DepositionSimulation<'_>,
    observer: &ConductionObserver<'a>,
) -> ConductionObserver<'a> {
    observer
        .clone()
        .with_material_density_layers(
            thermal,
            &thermal
                .cells
                .keys()
                .map(|&cell_id| PlanarDensityLayer {
                    cell_id,
                    footprint_area_m2: 1.0,
                })
                .collect::<Vec<_>>(),
            || Ok(()),
        )
        .unwrap()
}

fn controls(dt: f64) -> ImplicitConductionControls {
    ImplicitConductionControls {
        maximum_dt_s: dt,
        maximum_halvings: 0,
        cell_tolerance_j: 1e-10,
        total_tolerance_j: 1e-9,
        ..ImplicitConductionControls::default()
    }
}

#[test]
fn expanding_and_contracting_layers_match_the_implicit_energy_equation() {
    for bath in [600.0, 150.0] {
        let (initial, observer) = setup(&EXPANDING, &[300.0], Some((bath, None)));
        let observer = layers(&initial, &observer);
        let chain = ImplicitConductionChain::new(
            &initial,
            &observer,
            &[(MaterialId(0), &EXPANDING)],
            || Ok(()),
        )
        .unwrap();
        let (next, report) = chain
            .advance_to(
                &initial,
                1.0,
                controls(1.0),
                &mut ImplicitConductionWorkspace::new(1).unwrap(),
                || Ok(()),
            )
            .unwrap();
        // d=T/600, G=600/T. The original BE balance is
        // T-300 = (600/T)*(bath-T), with no mass or extra p*dV source.
        let expected = 0.5 * ((90_000.0 + 2400.0 * bath).sqrt() - 300.0);
        near(next.temperature_k(CellId(0)).unwrap(), expected, 1e-10);
        near(report.energy.external_heat_j, expected - 300.0, 1e-10);
        near(report.energy.residual_j, 0.0, 1e-10);
        assert_eq!(report.energy.substeps, 1);
        assert_eq!(report.rejected_intervals, 0);
        assert_eq!(next.cells[&CellId(0)].mass_kg, 1.0);
        assert_eq!(
            next.cells[&CellId(0)].geometry,
            initial.cells[&CellId(0)].geometry
        );
        assert_eq!(initial.temperature_k(CellId(0)).unwrap(), 300.0);
    }
}

#[test]
fn density_and_conductivity_jacobian_matches_link_and_both_bath_residuals() {
    let law = Expanding {
        nonlinear_k: true,
        ..EXPANDING
    };
    let (initial, observer) = setup(&law, &[360.0, 450.0, 600.0], Some((500.0, Some(450.0))));
    let observer = layers(&initial, &observer);
    let chain =
        ImplicitConductionChain::new(&initial, &observer, &[(MaterialId(0), &law)], || Ok(()))
            .unwrap();
    let step = BackwardEulerStep::new(&initial, 0.7, || Ok(())).unwrap();
    let energy = [20.0, -30.0, 10.0];
    let mut report = empty_report();
    let mut evaluate = |energy: &[f64]| {
        chain
            .evaluate(
                &step,
                energy,
                &BTreeMap::new(),
                &[],
                false,
                controls(1.0),
                &mut report,
                &mut || Ok(()),
            )
            .unwrap()
    };
    let trial = evaluate(&energy);
    let mut workspace = ImplicitConductionWorkspace::new(3).unwrap();
    chain
        .jacobian(&initial, trial.payload(), 0.7, &mut workspace, &mut || {
            Ok(())
        })
        .unwrap();
    for column in 0..3 {
        let mut plus = energy;
        let mut minus = energy;
        plus[column] += 1e-3;
        minus[column] -= 1e-3;
        let a = evaluate(&plus);
        let b = evaluate(&minus);
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
                2e-8,
            );
        }
    }
}

#[test]
fn closed_density_layers_conserve_mass_and_energy_in_chain_graph_and_sdirk() {
    let (initial, observer) = setup(&EXPANDING, &[500.0, 200.0], None);
    let observer = layers(&initial, &observer);
    let chain =
        ImplicitConductionChain::new(&initial, &observer, &[(MaterialId(0), &EXPANDING)], || {
            Ok(())
        })
        .unwrap();
    let graph = ImplicitConductionGraph::new(
        &initial,
        &observer,
        &[(MaterialId(0), &EXPANDING)],
        ImplicitConductionLinearControls::default(),
        || Ok(()),
    )
    .unwrap();
    let mut workspace = ImplicitConductionWorkspace::new(2).unwrap();
    let euler = chain
        .advance_to(&initial, 20.0, controls(20.0), &mut workspace, || Ok(()))
        .unwrap();
    // m=[1,2], rho=300/T: conservation keeps d1+d2=1.5, hence G=2/3.
    // The temperature difference has exactly unit decay rate.
    let gamma = 1.0 - 1.0 / 2.0_f64.sqrt();
    let sdirk_factor = (2.0 * gamma) / (1.0 + gamma).powi(2);
    let chain_sdirk = chain
        .advance_sdirk2_to(&initial, 1.0, controls(1.0), &mut workspace, || Ok(()))
        .unwrap();
    let graph_sdirk = graph
        .advance_sdirk2_to(&initial, 1.0, controls(1.0), &mut workspace, || Ok(()))
        .unwrap();
    for ((next, report), factor) in [
        (euler, 1.0 / 21.0),
        (chain_sdirk, sdirk_factor),
        (graph_sdirk, sdirk_factor),
    ] {
        near(
            next.temperature_k(CellId(0)).unwrap(),
            300.0 + 200.0 * factor,
            1e-9,
        );
        near(
            next.temperature_k(CellId(1)).unwrap(),
            300.0 - 100.0 * factor,
            1e-9,
        );
        near(next.energy_j(), initial.energy_j(), 1e-9);
        near(report.energy.external_heat_j, 0.0, 0.0);
        near(report.energy.residual_j, 0.0, 1e-9);
        for (&id, cell) in next.cells() {
            assert_eq!(cell.mass_kg, initial.cells[&id].mass_kg);
            assert_eq!(cell.geometry, initial.cells[&id].geometry);
        }
    }
}

#[test]
fn constant_density_reduces_to_fixed_geometry_and_changed_mass_or_missing_tangent_refuse() {
    let law = Expanding {
        variable_density: false,
        ..EXPANDING
    };
    let (initial, fixed) = setup(&law, &[500.0, 200.0], None);
    let sourced = layers(&initial, &fixed);
    let constant = fixed
        .clone()
        .with_planar_phase_layers(
            &initial,
            &[
                PlanarPhaseLayer {
                    cell_id: CellId(0),
                    density_kg_m3: 1.0,
                    footprint_area_m2: 1.0,
                },
                PlanarPhaseLayer {
                    cell_id: CellId(1),
                    density_kg_m3: 1.0,
                    footprint_area_m2: 1.0,
                },
            ],
            || Ok(()),
        )
        .unwrap();
    let mut solutions = Vec::new();
    for observer in [&fixed, &sourced, &constant] {
        let chain =
            ImplicitConductionChain::new(&initial, observer, &[(MaterialId(0), &law)], || Ok(()))
                .unwrap();
        let mut workspace = ImplicitConductionWorkspace::new(2).unwrap();
        let (next, _) = chain
            .advance_to(&initial, 1.0, controls(1.0), &mut workspace, || Ok(()))
            .unwrap();
        solutions.push(next.cells().clone());
        let mut changed = initial.clone();
        changed.cells.get_mut(&CellId(0)).unwrap().mass_kg = 2.0;
        assert!(matches!(
            chain.advance_to(&changed, 1.0, controls(1.0), &mut workspace, || Ok(())),
            Err(ImplicitConductionError::Invalid(
                "implicit reference cell binding changed"
            ))
        ));
    }
    assert_eq!(solutions[0], solutions[1]);
    assert_eq!(solutions[0], solutions[2]);
    let missing = Expanding {
        density_tangent: false,
        ..EXPANDING
    };
    let (initial, observer) = setup(&missing, &[300.0], Some((600.0, None)));
    assert!(
        ImplicitConductionChain::new(&initial, &observer, &[(MaterialId(0), &missing)], || Ok(()))
            .is_ok()
    );
    let observer = layers(&initial, &observer);
    assert!(matches!(
        ImplicitConductionChain::new(&initial, &observer, &[(MaterialId(0), &missing)], || Ok(())),
        Err(ImplicitConductionError::Invalid(
            "implicit material-density layers require an original specific-volume tangent"
        ))
    ));
}

impl ImplicitLinkLaw for Expanding {
    fn caloric_laws(&self) -> [&dyn ThermalLaw; 2] {
        [self, self]
    }
    fn state_at_specific_enthalpies(
        &self,
        _: [f64; 2],
        _: ConductiveLink,
        _: ContactPatch,
        _: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<ImplicitLinkState, DepositionError> {
        panic!("fixed-geometry flux must refuse before evaluation")
    }
}
impl ImplicitBoundaryLaw for Expanding {
    fn caloric_law(&self) -> &dyn ThermalLaw {
        self
    }
    fn state_at_specific_enthalpy(
        &self,
        _: f64,
        _: f64,
        _: ImplicitBoundaryGeometry,
        _: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<ImplicitBoundaryState, DepositionError> {
        panic!("fixed-geometry flux must refuse before evaluation")
    }
}
impl ImplicitBoundaryExchangeLaw for Expanding {
    fn caloric_law(&self) -> &dyn ThermalLaw {
        self
    }
    fn state_at_specific_enthalpy(
        &self,
        _: f64,
        _: f64,
        _: ImplicitBoundaryGeometry,
        _: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<ImplicitBoundaryExchangeState, DepositionError> {
        panic!("fixed-geometry flux must refuse before evaluation")
    }
}

#[test]
fn integrated_fluxes_refuse_density_feedback_without_geometry_tangents() {
    let (initial, observer) = setup(&EXPANDING, &[300.0, 400.0], Some((600.0, Some(450.0))));
    let observer = layers(&initial, &observer);
    let build = || {
        ImplicitConductionChain::new(&initial, &observer, &[(MaterialId(0), &EXPANDING)], || {
            Ok(())
        })
        .unwrap()
    };
    let refused = |result: Result<ImplicitConductionChain<'_>, ImplicitConductionError>| {
        assert!(matches!(
            result,
            Err(ImplicitConductionError::Invalid(
                "integrated flux laws do not supply material-density geometry derivatives"
            ))
        ));
    };
    refused(build().with_link_laws(&initial, &[(0, &EXPANDING)], || Ok(())));
    refused(build().with_boundary_laws(&initial, &[(CellId(0), &EXPANDING)], &[], || Ok(())));
    refused(build().with_boundary_laws(&initial, &[], &[(CellId(0), &EXPANDING)], || Ok(())));
    refused(build().with_boundary_exchanges(&initial, &[(CellId(0), &EXPANDING)], &[], || Ok(())));
    refused(build().with_boundary_exchanges(&initial, &[], &[(CellId(0), &EXPANDING)], || Ok(())));
}
