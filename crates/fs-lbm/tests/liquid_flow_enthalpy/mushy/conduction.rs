//! Phase-dependent heat transport and momentum use one accepted h history.
use super::*;
use fs_lbm::conjugate::enthalpy::advection::feedback::MaterialFlowError;
use fs_material::fluid_transport::FluidTransportState;
use fs_material::phase::conductivity::LiquidMassFractionConductivity;

fn thermal_curve(card: &MaterialCard) -> EquilibriumEnthalpyPhaseCurve {
    EquilibriumEnthalpyPhaseCurve::try_new(
        card.content_hash(),
        vec![
            EnthalpyPhaseKnot {
                specific_enthalpy_j_kg: 0.0,
                temperature_k: 500.0,
                liquid_mass_fraction: 0.0,
                bulk_density_kg_m3: 1.0,
            },
            EnthalpyPhaseKnot {
                specific_enthalpy_j_kg: 100.0,
                temperature_k: 500.0,
                liquid_mass_fraction: 1.0,
                bulk_density_kg_m3: 1.0,
            },
            EnthalpyPhaseKnot {
                specific_enthalpy_j_kg: 300.0,
                temperature_k: 600.0,
                liquid_mass_fraction: 1.0,
                bulk_density_kg_m3: 1.0,
            },
        ],
    )
    .unwrap()
}

fn conducting_initial<'a>(
    domain: &VoxelDomain,
    curve: &'a EquilibriumEnthalpyPhaseCurve,
) -> DepositionSimulation<'a> {
    let seed = initial(domain, curve, 0.0);
    let h = [25.0, 75.0, 150.0, 170.0];
    let mut events: Vec<_> = seed
        .cells()
        .values()
        .cloned()
        .map(|mut cell| {
            cell.specific_h_j_kg = h[cell.id.0 as usize];
            TimedEvent {
                time_s: 0.0,
                event: Event::Birth(cell),
            }
        })
        .collect();
    for c in 0..3 {
        events.push(TimedEvent {
            time_s: 0.0,
            event: Event::Contact(ThermalLink {
                id: 100 + c,
                a: CellId(c),
                b: CellId(c + 1),
                area_m2: 1.0,
                // The actual source epoch must replace this authored value.
                conductance_w_k: 0.001,
            }),
        });
    }
    DepositionSimulation::new(
        events,
        &[(MAT, curve as &dyn ThermalLaw)],
        DepositionBudget {
            max_cells: 4,
            max_links: 3,
            max_events: 7,
            max_substeps: 1000,
        },
    )
    .unwrap()
    .advance_to(0.0, config(), || Ok(()))
    .unwrap()
    .0
}

fn temperature(h: f64) -> f64 {
    if h <= 100.0 { 500.0 } else { 500.0 + (h - 100.0) / 2.0 }
}

fn fraction(h: f64) -> f64 {
    if h <= 100.0 { h / 100.0 } else { 1.0 }
}

#[test]
fn phase_data_changes_native_conduction_and_darcy_feedback_at_constant_melting_temperature() {
    let card = card(claims(1.0, 1.0, false));
    let liquid = liquid(&card);
    let curve = thermal_curve(&card);
    let domain = VoxelDomain::new(4, 1, 1, 1.0).unwrap();
    let config_flow = momentum();
    let gate = CancelGate::new();
    let ids: Vec<_> = (0..4).map(|c| Some(CellId(c))).collect();
    let old = conducting_initial(&domain, &curve);
    let mobility = MushyPhaseLaw::try_new(&curve, 1.0, 1.0, 1.0, 0.1).unwrap();
    let mut outputs = Vec::new();
    for solid in [1.0, 4.0] {
        let law = LiquidMassFractionConductivity::declared(solid, 1.0).unwrap();
        let observer = MaterialFlowFeedback::try_new(
            &old, &domain, input(&ids, &curve), &config_flow,
            Heater::default(), &gate, 100,
            |_, _: &DepositedCell, t| {
                liquid.at_temperature(Temperature::new(t)).map_err(|e| {
                    fs_material::MaterialError::State { what: e.to_string() }
                })
            },
        ).unwrap()
            .with_mushy_mobility(&mobility).unwrap()
            .with_fluid_conduction(&old, &[100, 101, 102]).unwrap()
            .with_phase_conductivity(&curve, law).unwrap();
        // These two cells have identical T, but their different latent h must
        // produce distinct half-cell conductivities when the phase law varies.
        assert_eq!(old.temperature_k(CellId(0)).unwrap(), 500.0);
        assert_eq!(old.temperature_k(CellId(1)).unwrap(), 500.0);
        let mut source = observer.clone();
        source.sources(&old, &mut || Ok(())).unwrap();
        let conductances = source.conductances(&old, &mut || Ok(())).unwrap();
        let k0 = 0.01 * (0.75 * solid + 0.25);
        let k1 = 0.01 * (0.25 * solid + 0.75);
        near(conductances[&100], 1.0 / (0.5 / k0 + 0.5 / k1));
        let (next, accepted, audit) = old.advance_observed(
            &observer, 0.1,
            StepConfig { max_dt_s: 0.0125, energy_tolerance_j: 1e-9 },
            || Ok(()),
        ).unwrap();
        let mut h = [25.0, 75.0, 150.0, 170.0];
        let mut external_j = 0.0;
        let mut last_q = 0.0;
        for &dt in &accepted.inner().durations {
            // Independent series Darcy/thermal resistance and donor-cell total
            // enthalpy recurrence, using the actual admitted native step sizes.
            let phi = h.map(fraction);
            let drag = phi.map(|f| (1.0 - f).powi(2) / f.powi(3));
            let q = 4.0 / drag.iter().sum::<f64>();
            let k = phi.map(|f| 0.01 * ((1.0 - f) * solid + f));
            let t = h.map(temperature);
            let previous = h;
            for c in 0..4 {
                let upstream = if c == 0 { 50.0 } else { previous[c - 1] };
                h[c] += dt * (q * (upstream - previous[c]) + if c == 1 { 100.0 } else { 0.0 });
            }
            for c in 0..3 {
                let g = 1.0 / (0.5 / k[c] + 0.5 / k[c + 1]);
                let heat = dt * g * (t[c + 1] - t[c]);
                h[c] += heat;
                h[c + 1] -= heat;
            }
            external_j += dt * (100.0 + q * (50.0 - previous[3]));
            last_q = q;
        }
        for (c, expected) in h.into_iter().enumerate() {
            near(next.cells()[&CellId(c as u64)].specific_h_j_kg, expected);
        }
        for c in [CellId(0), CellId(1)] {
            assert_eq!(next.temperature_k(c).unwrap(), 500.0);
            assert_ne!(next.cells()[&c].specific_h_j_kg, old.cells()[&c].specific_h_j_kg);
        }
        near(accepted.applied().unwrap().flow.field.flux_x(2, 0, 0), last_q);
        near(next.energy_j() - old.energy_j(), external_j);
        near(audit.external_heat_j, external_j);
        near(next.mass_kg(), old.mass_kg());
        assert_eq!(accepted.phase_conductivity_law(), Some(law));
        assert_eq!(observer.momentum_solves(), 0);
        outputs.push((h[1], last_q));
    }
    assert!(outputs[1].0 > outputs[0].0);
    assert_ne!(outputs[1].1, outputs[0].1);
}

#[test]
fn identity_phase_law_preserves_the_unmodified_coupled_trajectory() {
    let card = card(claims(1.0, 1.0, false));
    let liquid = liquid(&card);
    let curve = thermal_curve(&card);
    let domain = VoxelDomain::new(4, 1, 1, 1.0).unwrap();
    let flow = momentum();
    let gate = CancelGate::new();
    let ids: Vec<_> = (0..4).map(|c| Some(CellId(c))).collect();
    let old = conducting_initial(&domain, &curve);
    let mobility = MushyPhaseLaw::try_new(&curve, 1.0, 1.0, 1.0, 0.1).unwrap();
    let mut outputs = Vec::new();
    for identity in [false, true] {
        let mut observer = MaterialFlowFeedback::try_new(
            &old, &domain, input(&ids, &curve), &flow, Heater::default(), &gate, 100,
            |_, _: &DepositedCell, t| liquid.at_temperature(Temperature::new(t))
                .map_err(|e| fs_material::MaterialError::State { what: e.to_string() }),
        ).unwrap().with_mushy_mobility(&mobility).unwrap()
            .with_fluid_conduction(&old, &[100, 101, 102]).unwrap();
        if identity {
            observer = observer.with_phase_conductivity(
                &curve, LiquidMassFractionConductivity::declared(1.0, 1.0).unwrap(),
            ).unwrap();
        }
        let (next, accepted, _) = old.advance_observed(
            &observer, 0.1,
            StepConfig { max_dt_s: 0.0125, energy_tolerance_j: 1e-9 },
            || Ok(()),
        ).unwrap();
        outputs.push((
            next.cells().values().map(|c| c.specific_h_j_kg.to_bits()).collect::<Vec<_>>(),
            accepted.applied().cloned(), accepted.inner().clone(), accepted.audit(),
        ));
    }
    assert_eq!(outputs[0], outputs[1]);
}

#[test]
fn foreign_chart_endpoint_overflow_and_cancellation_preserve_accepted_heat_and_flow() {
    let card = card(claims(1.0, 1.0, false));
    let curve = thermal_curve(&card);
    let copy = curve.clone();
    let domain = VoxelDomain::new(4, 1, 1, 1.0).unwrap();
    let flow = momentum();
    let gate = CancelGate::new();
    let ids: Vec<_> = (0..4).map(|c| Some(CellId(c))).collect();
    let old = initial(&domain, &curve, 0.0);
    let mobility = MushyPhaseLaw::try_new(&curve, 1.0, 1.0, 1.0, 0.1).unwrap();
    let law = LiquidMassFractionConductivity::declared(2.0, 2.0).unwrap();
    let observer = MaterialFlowFeedback::try_new(
        &old, &domain, input(&ids, &curve), &flow, Heater::default(), &gate, 100,
        |_, cell: &DepositedCell, _| {
            // A bounded arithmetic-refusal fixture, not a physical material.
            // The base state is representable; only the new phase-scaled
            // endpoint conductivity overflows after the first valid heat step.
            let k = if cell.specific_h_j_kg > 50.0 { f64::MAX } else { 1.0 };
            FluidTransportState::try_new(1.0, 2.0, k, 1.0)
        },
    ).unwrap().with_mushy_mobility(&mobility).unwrap();
    assert!(observer.clone().with_phase_conductivity(&copy, law).is_err());
    let observer = observer.with_phase_conductivity(&curve, law).unwrap();
    assert!(matches!(
        old.advance_observed(&observer, 0.01, config(), || Ok(())),
        Err(MaterialFlowError::Thermal(DepositionError::InvalidInput(
            "phase-scaled fluid conductivity is not representable"
        )))
    ));
    assert_eq!(old.time_s(), 0.0);
    near(old.energy_j(), 200.0);
    assert_eq!(observer.momentum_solves(), 0);
    assert!(observer.applied().is_none());
    assert_eq!(observer.inner(), &Heater::default());
    gate.request();
    assert!(matches!(
        old.advance_observed(&observer, 0.01, config(), || Ok(())),
        Err(MaterialFlowError::Thermal(DepositionError::Cancelled))
    ));
    assert_eq!(observer.inner(), &Heater::default());
    near(old.energy_j(), 200.0);
}
