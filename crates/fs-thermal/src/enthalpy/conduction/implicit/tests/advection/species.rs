//! Frozen mixture compositions share the original native heat/flow solve.
use super::*;
use crate::enthalpy::advection::species::{CaloricComposition, InletSpecies};

struct MixtureLaw {
    names: Vec<String>,
    fractions: Vec<f64>,
    datum: f64,
    reference: &'static str,
}

impl MixtureLaw {
    fn new(fractions: Vec<f64>, datum: f64) -> Self {
        Self {
            names: vec!["A".into(), "B".into()],
            fractions,
            datum,
            reference: "fixed-species-formation-reference",
        }
    }
    fn cp(&self) -> f64 {
        2.0 * self.fractions[0] + 6.0 * self.fractions[1]
    }
    fn formation(&self) -> f64 {
        1000.0 * self.fractions[0] - 400.0 * self.fractions[1]
    }
    fn h(&self, t: f64) -> f64 {
        self.datum + self.formation() + self.cp() * t
    }
}

impl ThermalLaw for MixtureLaw {
    fn temperature_k(&self, h: f64) -> Result<f64, DepositionError> {
        let t = ((h - self.datum) - self.formation()) / self.cp();
        if t.is_finite() && (100.0..=1000.0).contains(&t) {
            Ok(t)
        } else {
            Err(DepositionError::MaterialDomain)
        }
    }
    fn maximum_temperature_slope(&self) -> f64 {
        1.0 / self.cp()
    }
    fn advected_composition(&self) -> Option<CaloricComposition<'_>> {
        Some(CaloricComposition {
            species: &self.names,
            mass_fractions: &self.fractions,
            enthalpy_reference: self.reference,
        })
    }
    fn specific_entropy_change_j_kg_k(
        &self,
        a: f64,
        b: f64,
    ) -> Result<Option<f64>, DepositionError> {
        // Thermal storage at this fixed composition; no mixing entropy.
        Ok(Some(
            self.cp() * (self.temperature_k(b)? / self.temperature_k(a)?).ln(),
        ))
    }
}

impl ConductivityLaw for MixtureLaw {
    fn conductivity_w_m_k(&self, h: f64) -> Result<f64, DepositionError> {
        self.temperature_k(h)?;
        Ok(1.0)
    }
}

impl ImplicitConductionLaw for MixtureLaw {
    fn caloric_law(&self) -> &dyn ThermalLaw {
        self
    }
    fn conductivity_law(&self) -> &dyn ConductivityLaw {
        self
    }
    fn enthalpy_bounds_j_kg(&self) -> [f64; 2] {
        [self.h(100.0), self.h(1000.0)]
    }
    fn state_at_specific_enthalpy(
        &self,
        h: f64,
    ) -> Result<ImplicitConductionState, DepositionError> {
        Ok(ImplicitConductionState {
            temperature_k: self.temperature_k(h)?,
            temperature_slope_kg_k_j: 1.0 / self.cp(),
            conductivity_w_m_k: 1.0,
            conductivity_slope_kg_w_m_k_j: 0.0,
        })
    }
}

fn bindings(laws: &[MixtureLaw]) -> Vec<(MaterialId, &dyn ImplicitConductionLaw)> {
    laws.iter()
        .enumerate()
        .map(|(i, law)| (MaterialId(i as u64), law as &dyn ImplicitConductionLaw))
        .collect()
}

fn prepare_mixtures<'a>(
    laws: &'a [MixtureLaw],
    temperatures: &[f64],
    masses: &[f64],
    conducting: bool,
) -> (DepositionSimulation<'a>, ConductionObserver<'a>) {
    let mut events: Vec<_> = temperatures
        .iter()
        .enumerate()
        .map(|(i, &t)| TimedEvent {
            time_s: 0.0,
            event: Event::Birth(DepositedCell {
                id: CellId(i as u64),
                material_id: MaterialId(i as u64),
                segment_id: i as u64,
                body_id: 0,
                mass_kg: masses[i],
                specific_h_j_kg: laws[i].h(t),
                geometry: CellGeometry {
                    start_m: [i as f64, 0.0, 0.0],
                    end_m: [i as f64 + 1.0, 0.0, 0.0],
                    volume_m3: masses[i],
                    surface_area_m2: 2.0,
                },
                boundary: Boundary {
                    ambient_k: 300.0,
                    htc_w_m2_k: 0.0,
                    reservoir: None,
                    source_w: 0.0,
                },
            }),
        })
        .collect();
    let links = if conducting {
        events.push(TimedEvent {
            time_s: 0.0,
            event: Event::Contact(ThermalLink {
                id: 0,
                a: CellId(0),
                b: CellId(1),
                area_m2: 1.0,
                conductance_w_k: 0.0,
            }),
        });
        vec![ConductiveLink {
            link_id: 0,
            half_distances_m: [0.5; 2],
            contact_resistance_m2_k_w: 0.0,
        }]
    } else {
        Vec::new()
    };
    let calorics: Vec<_> = laws
        .iter()
        .enumerate()
        .map(|(i, law)| (MaterialId(i as u64), law as &dyn ThermalLaw))
        .collect();
    let conductivity: Vec<_> = laws
        .iter()
        .enumerate()
        .map(|(i, law)| (MaterialId(i as u64), law as &dyn ConductivityLaw))
        .collect();
    let budget = DepositionBudget {
        max_cells: temperatures.len(),
        max_links: 4,
        max_events: events.len(),
        max_substeps: 1024,
    };
    let thermal = DepositionSimulation::new(events, &calorics, budget)
        .unwrap()
        .advance_to(
            0.0,
            StepConfig {
                max_dt_s: 1.0,
                energy_tolerance_j: 1e-8,
            },
            || Ok(()),
        )
        .unwrap()
        .0;
    let conduction =
        ConductionObserver::new(&thermal, &conductivity, &links, &[], || Ok(())).unwrap();
    (thermal, conduction)
}

fn single_mixer(datum: f64) -> (Vec<MixtureLaw>, Advected<()>) {
    let laws = vec![
        MixtureLaw::new(vec![0.25, 0.75], datum),
        MixtureLaw::new(vec![1.0, 0.0], datum),
        MixtureLaw::new(vec![0.0, 1.0], datum),
    ];
    let flow = Advected::new((), vec![], 0)
        .unwrap()
        .with_boundaries(
            vec![
                MassBoundary::Inlet {
                    id: 10,
                    cell: CellId(0),
                    kg_s: 0.25,
                    material_id: MaterialId(1),
                    specific_h_j_kg: laws[1].h(400.0),
                },
                MassBoundary::Inlet {
                    id: 11,
                    cell: CellId(0),
                    kg_s: 0.75,
                    material_id: MaterialId(2),
                    specific_h_j_kg: laws[2].h(300.0),
                },
                MassBoundary::Outlet {
                    id: 12,
                    cell: CellId(0),
                    kg_s: 1.0,
                },
            ],
            3,
        )
        .unwrap()
        .with_stationary_species(
            vec!["A".into(), "B".into()],
            vec![
                InletSpecies {
                    boundary_id: 10,
                    mass_fractions: vec![1.0, 0.0],
                },
                InletSpecies {
                    boundary_id: 11,
                    mass_fractions: vec![0.0, 1.0],
                },
            ],
            1,
            2,
        )
        .unwrap();
    (laws, flow)
}

#[test]
fn stationary_mixed_inlets_use_native_chain_above_cfl_and_preserve_enthalpy_datum() {
    for datum in [0.0, 1_099_511_627_776.0] {
        let (laws, flow) = single_mixer(datum);
        let (thermal, conduction) = prepare_mixtures(&laws, &[280.0], &[2.0], false);
        let system =
            ImplicitConductionChain::new(&thermal, &conduction, &bindings(&laws), || Ok(()))
                .unwrap()
                .with_implicit_advection(&thermal, &flow, || Ok(()))
                .unwrap();
        let (next, report) = system
            .advance_to(
                &thermal,
                18.0,
                ImplicitConductionControls {
                    maximum_halvings: 0,
                    ..controls(18.0)
                },
                &mut ImplicitConductionWorkspace::new(1).unwrap(),
                || Ok(()),
            )
            .unwrap();
        // Inlet h_mix=.25*1800+.75*1400=1500. At lambda=dt*q/m=9,
        // BE gives h=(1350+9*1500)/10=1485 and T=(1485+50)/5=307.
        near(
            next.cells()[&CellId(0)].specific_h_j_kg - datum,
            1485.0,
            1e-10,
        );
        near(next.temperature_k(CellId(0)).unwrap(), 307.0, 1e-10);
        near(report.energy.external_heat_j, 270.0, 1e-9);
        near(report.energy.stored_energy_change_j, 270.0, 1e-9);
        near(report.energy.residual_j, 0.0, 1e-9);
        assert_eq!(report.energy.substeps, 1);
        assert_eq!(report.rejected_intervals, 0);
        assert!(report.entropy.advection.is_none());
        assert!(report.entropy.has_unclassified_sources);
        assert_eq!(next.cells()[&CellId(0)].mass_kg, 2.0);
        assert_eq!(
            flow.stationary_species().unwrap().cells()[&CellId(0)],
            vec![0.25, 0.75]
        );
        assert_eq!(flow.audit().substeps, 0);
        assert_eq!(thermal.time_s(), 0.0);
    }
}

#[test]
fn inlet_fed_mixture_recirculation_shares_conduction_and_actual_sdirk_stage_fluxes() {
    let a = MixtureLaw::new(vec![1.0, 0.0], 0.0);
    let b = MixtureLaw::new(vec![0.0, 1.0], 0.0);
    let inlet = [a.h(400.0), b.h(300.0)];
    let flow = Advected::new(
        (),
        vec![
            MassFlow {
                id: 10,
                a: CellId(0),
                b: CellId(1),
                kg_s: 2.0,
            },
            MassFlow {
                id: 11,
                a: CellId(1),
                b: CellId(0),
                kg_s: 1.0,
            },
        ],
        2,
    )
    .unwrap()
    .with_boundaries(
        vec![
            MassBoundary::Inlet {
                id: 20,
                cell: CellId(0),
                kg_s: 1.0,
                material_id: MaterialId(2),
                specific_h_j_kg: inlet[0],
            },
            MassBoundary::Inlet {
                id: 21,
                cell: CellId(1),
                kg_s: 1.0,
                material_id: MaterialId(3),
                specific_h_j_kg: inlet[1],
            },
            MassBoundary::Outlet {
                id: 22,
                cell: CellId(1),
                kg_s: 2.0,
            },
        ],
        5,
    )
    .unwrap()
    .with_stationary_species(
        vec!["A".into(), "B".into()],
        vec![
            InletSpecies {
                boundary_id: 20,
                mass_fractions: vec![1.0, 0.0],
            },
            InletSpecies {
                boundary_id: 21,
                mass_fractions: vec![0.0, 1.0],
            },
        ],
        2,
        2,
    )
    .unwrap();
    let compositions = flow.stationary_species().unwrap().cells();
    near(compositions[&CellId(0)][0], 0.75, 1e-14);
    near(compositions[&CellId(1)][0], 0.5, 1e-14);
    let laws = vec![
        MixtureLaw::new(compositions[&CellId(0)].clone(), 0.0),
        MixtureLaw::new(compositions[&CellId(1)].clone(), 0.0),
        a,
        b,
    ];
    let masses = [1.0, 2.0];
    let (thermal, conduction) = prepare_mixtures(&laws, &[280.0, 360.0], &masses, true);
    let system = ImplicitConductionGraph::new(
        &thermal,
        &conduction,
        &bindings(&laws),
        Default::default(),
        || Ok(()),
    )
    .unwrap()
    .with_implicit_advection(&thermal, &flow, || Ok(()))
    .unwrap();
    let old = [laws[0].h(280.0), laws[1].h(360.0)];
    // Independent two-cell balance: convection plus G=1 W/K conduction.
    // No temperature is advected between unlike Cp/formation-enthalpy laws.
    let slope = [
        [-2.0 - 1.0 / laws[0].cp(), 1.0 + 1.0 / laws[1].cp()],
        [2.0 + 1.0 / laws[0].cp(), -3.0 - 1.0 / laws[1].cp()],
    ];
    let offset = laws[0].formation() / laws[0].cp() - laws[1].formation() / laws[1].cp();
    let forcing = [inlet[0] + offset, inlet[1] - offset];
    let rates = |h: [f64; 2]| {
        [
            slope[0][0] * h[0] + slope[0][1] * h[1] + forcing[0],
            slope[1][0] * h[0] + slope[1][1] * h[1] + forcing[1],
        ]
    };
    let solve = |span: f64, history: [f64; 2]| {
        let a = masses[0] - span * slope[0][0];
        let b = -span * slope[0][1];
        let c = -span * slope[1][0];
        let d = masses[1] - span * slope[1][1];
        let rhs = [
            history[0] + span * forcing[0],
            history[1] + span * forcing[1],
        ];
        let determinant = a * d - b * c;
        [
            (d * rhs[0] - b * rhs[1]) / determinant,
            (a * rhs[1] - c * rhs[0]) / determinant,
        ]
    };
    let dt = 2.0; // Courant numbers 4 and 3.
    let history = [masses[0] * old[0], masses[1] * old[1]];
    for second_order in [false, true] {
        let (expected, external) = if second_order {
            let gamma = 1.0 - 0.5_f64.sqrt();
            let stage = solve(gamma * dt, history);
            let power = rates(stage);
            let next = solve(
                gamma * dt,
                [
                    history[0] + (1.0 - gamma) * dt * power[0],
                    history[1] + (1.0 - gamma) * dt * power[1],
                ],
            );
            let outlet = (1.0 - gamma) * stage[1] + gamma * next[1];
            (next, dt * (inlet[0] + inlet[1] - 2.0 * outlet))
        } else {
            let next = solve(dt, history);
            (next, dt * (inlet[0] + inlet[1] - 2.0 * next[1]))
        };
        let mut workspace = ImplicitConductionWorkspace::new(2).unwrap();
        let control = ImplicitConductionControls {
            maximum_halvings: 0,
            ..controls(dt)
        };
        let (next, report) = if second_order {
            system.advance_sdirk2_to(&thermal, dt, control, &mut workspace, || Ok(()))
        } else {
            system.advance_to(&thermal, dt, control, &mut workspace, || Ok(()))
        }
        .unwrap();
        for (i, &h) in expected.iter().enumerate() {
            near(next.cells()[&CellId(i as u64)].specific_h_j_kg, h, 1e-8);
            near(
                next.temperature_k(CellId(i as u64)).unwrap(),
                laws[i].temperature_k(h).unwrap(),
                1e-8,
            );
        }
        near(report.energy.external_heat_j, external, 1e-8);
        near(report.energy.stored_energy_change_j, external, 1e-8);
        near(report.energy.residual_j, 0.0, 1e-9);
        assert_eq!(report.energy.substeps, 1);
        assert_eq!(report.rejected_intervals, 0);
        assert!(report.entropy.advection.is_none());
    }
    assert_eq!(flow.audit().substeps, 0);
    assert_eq!(thermal.time_s(), 0.0);
}

#[test]
fn implicit_stationary_species_keep_caloric_bindings_and_refuse_chemical_entropy_claims() {
    let (laws, flow) = single_mixer(0.0);
    let (thermal, conduction) = prepare_mixtures(&laws, &[280.0], &[2.0], false);
    assert!(matches!(
        flow.endpoint_transport_with_entropy(&thermal, || Ok(())),
        Err(DepositionError::InvalidInput(
            "thermal advection entropy does not cover stationary species mixing"
        ))
    ));
    let system = ImplicitConductionChain::new(&thermal, &conduction, &bindings(&laws), || Ok(()))
        .unwrap()
        .with_implicit_advection(&thermal, &flow, || Ok(()))
        .unwrap();
    assert!(matches!(
        system.with_implicit_advection_entropy(&thermal, || Ok(())),
        Err(ImplicitConductionError::Thermal(
            DepositionError::InvalidInput(
                "thermal advection entropy does not cover stationary species mixing",
            )
        ))
    ));
    for wrong_reference in [false, true] {
        let (mut laws, flow) = single_mixer(0.0);
        if wrong_reference {
            laws[0].reference = "different-formation-reference";
        } else {
            laws[0].fractions = vec![0.5, 0.5];
        }
        let (thermal, conduction) = prepare_mixtures(&laws, &[280.0], &[2.0], false);
        let system =
            ImplicitConductionChain::new(&thermal, &conduction, &bindings(&laws), || Ok(()))
                .unwrap();
        assert!(
            system
                .with_implicit_advection(&thermal, &flow, || Ok(()))
                .is_err()
        );
        assert_eq!(thermal.time_s(), 0.0);
        assert_eq!(flow.audit().substeps, 0);
    }
}

#[test]
fn scheduled_pump_retains_frozen_species_and_original_inlet_calorics() {
    use crate::enthalpy::advection::schedule::{AdvectionSchedule, AdvectionScheduleKnot};
    for datum in [0.0, 1_099_511_627_776.0] {
        let (laws, flow) = single_mixer(datum);
        let (thermal, conduction) = prepare_mixtures(&laws, &[280.0], &[2.0], false);
        let knots = vec![
            AdvectionScheduleKnot {
                time_s: 0.0, flow_multiplier: 0.0,
                inlet_enthalpies_j_kg: BTreeMap::from([
                    (10, laws[1].h(400.0)), (11, laws[2].h(300.0)),
                ]),
            },
            AdvectionScheduleKnot {
                time_s: 1.0, flow_multiplier: 2.0,
                inlet_enthalpies_j_kg: BTreeMap::from([
                    (10, laws[1].h(600.0)), (11, laws[2].h(200.0)),
                ]),
            },
        ];
        let schedule = AdvectionSchedule::new(knots.clone(), 2, 2, || Ok(())).unwrap();
        let stopped = schedule.flow_at(&flow, 0.0, &mut || Ok(())).unwrap().0;
        assert_eq!(stopped.stationary_species().unwrap().cells(),
            flow.stationary_species().unwrap().cells());
        assert!(stopped.endpoint_transport(&thermal, || Ok(())).unwrap().heat.bulk.sources_w.is_empty());
        let system = ImplicitConductionChain::new(
            &thermal, &conduction, &bindings(&laws), || Ok(()),
        ).unwrap().with_scheduled_implicit_advection(&thermal, &flow, &schedule, || Ok(())).unwrap();
        let (next, report) = system.advance_to(
            &thermal, 1.0, controls(1.0), &mut ImplicitConductionWorkspace::new(1).unwrap(), || Ok(()),
        ).unwrap();
        // Final original inlet enthalpies mix to h=1150+datum; BE with
        // mass=2 and total mdot=2 gives h=(1350+1150)/2 = 1250+datum.
        near(next.temperature_k(CellId(0)).unwrap(), 260.0, 1e-10);
        near(report.energy.external_heat_j, -200.0, 1e-9);
        near(report.energy.residual_j, 0.0, 1e-9);
        assert!(report.entropy.has_unclassified_sources);
        assert!(report.entropy.advection.is_none());
        assert_eq!(flow.audit().substeps, 0);
        assert_eq!(flow.stationary_species().unwrap().cells()[&CellId(0)], vec![0.25, 0.75]);

        // This value is inside the receiving mixture's h-domain but outside
        // species A's original inlet law. A schedule must reject the latter.
        let mut bad_knots = knots;
        bad_knots[1].inlet_enthalpies_j_kg.insert(10, laws[1].h(1050.0));
        let bad = AdvectionSchedule::new(bad_knots, 2, 2, || Ok(())).unwrap();
        assert!(matches!(ImplicitConductionChain::new(
            &thermal, &conduction, &bindings(&laws), || Ok(()),
        ).unwrap().with_scheduled_implicit_advection(&thermal, &flow, &bad, || Ok(())),
            Err(ImplicitConductionError::Thermal(DepositionError::MaterialDomain))));
        assert_eq!(thermal.time_s(), 0.0);
    }
}
