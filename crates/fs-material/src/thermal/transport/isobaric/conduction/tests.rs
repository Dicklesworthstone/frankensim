use super::*;
use crate::thermal::transport::{SPECIFIC_HEAT, isobaric::IsobaricHeatCapacityQuery};
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialStateId, PropertyClaim, PropertyValue, Provenance,
    QueryPoint, SelectionPolicy, UncertaintyModel,
};
use fs_qty::{Pressure, Temperature};
use fs_thermal::enthalpy::{
    Boundary, CellGeometry, CellId, DepositedCell, DepositionBudget, DepositionSimulation, Event,
    MaterialId, StepConfig, ThermalLink, TimedEvent,
    conduction::{ConductionObserver, ConductiveLink},
    conduction::implicit::{
        ImplicitConductionChain, ImplicitConductionControls, ImplicitConductionWorkspace,
    },
};

struct Fixture {
    card: MaterialCard,
    cp: PropertyKey,
    conductivity: PropertyKey,
    lower: QueryPoint,
    upper: QueryPoint,
}

const SELECTION: ClaimSelection = ClaimSelection::Policy(SelectionPolicy::SingleClaimOnly);

impl Fixture {
    fn new(conductivity_pressure: Option<f64>, conductivities: [f64; 4]) -> Self {
        let temperature = QuantitySpec::semantic(SemanticType::new(
            QuantityKind::AbsoluteTemperature,
            ValueForm::Static,
        ));
        let pressure = QuantitySpec::dimensional(Pressure::DIMS);
        let cp = PropertyKey::new("isobaric_cp", SPECIFIC_HEAT);
        let conductivity = PropertyKey::new("scalar_conductivity", CONDUCTIVITY);
        let validity = || {
            ValidityDomain::unconstrained().with_quantity("T", temperature, 300.0, 700.0)
        };
        let claim = |key: PropertyKey, dims, knots, validity| PropertyClaim {
            value: PropertyValue::Curve {
                abscissa: "T".into(),
                abscissa_dims: Temperature::DIMS,
                knots,
                dims,
            },
            key,
            validity,
            interpolation: InterpolationPolicy::LinearInside,
            uncertainty: UncertaintyModel::Unstated,
            observations: vec![],
            provenance: Provenance {
                source: "synthetic piecewise-affine Cp and conductivity fixture".into(),
                license: "test-only".into(),
                artifact: None,
            },
        };
        let mut claims = ClaimSet::new();
        claims
            .insert_claim(claim(
                cp.clone(),
                SPECIFIC_HEAT,
                vec![(300.0, 1000.0), (500.0, 1400.0), (700.0, 2000.0)],
                validity().with_quantity("p", pressure, 100000.0, 100000.0),
            ))
            .unwrap();
        let k_validity = if let Some(p) = conductivity_pressure {
            validity().with_quantity("p", pressure, p, p)
        } else {
            validity()
        };
        claims
            .insert_claim(claim(
                conductivity.clone(),
                CONDUCTIVITY,
                [300.0, 425.0, 550.0, 700.0]
                    .into_iter()
                    .zip(conductivities)
                    .collect(),
                k_validity,
            ))
            .unwrap();
        let card = MaterialCard::assemble(
            MaterialStateId {
                chemistry: "synthetic thermal material".into(),
                phase: "single phase".into(),
                process: "numerical fixture".into(),
                revision: 0,
            },
            claims,
            vec![],
        )
        .unwrap();
        let point = |t| {
            QueryPoint::new()
                .with_quantity("T", temperature, t)
                .unwrap()
                .with_quantity("p", pressure, 100000.0)
                .unwrap()
        };
        Self {
            card,
            cp,
            conductivity,
            lower: point(300.0),
            upper: point(700.0),
        }
    }

    fn caloric(&self) -> IsobaricHeatCapacityLaw {
        IsobaricHeatCapacityLaw::try_new(
            &self.card,
            IsobaricHeatCapacityQuery {
                specific_heat_cp: &self.cp,
                temperature_axis: "T",
                pressure_axis: "p",
                lower: &self.lower,
                upper: &self.upper,
                selection: SELECTION,
            },
            300.0,
        )
        .unwrap()
    }
}

fn fixture() -> Fixture {
    Fixture::new(Some(100000.0), [10.0, 15.0, 12.0, 18.0])
}

// Independent closed-form primitive and inverse of the authored fixture Cp.
fn enthalpy(t: f64) -> f64 {
    if t <= 500.0 {
        let x = t - 300.0;
        1000.0 * x + x * x
    } else {
        let x = t - 500.0;
        240000.0 + 1400.0 * x + 1.5 * x * x
    }
}

fn temperature(h: f64) -> f64 {
    if h <= 240000.0 {
        300.0 + 2.0 * h / (1000.0 + (1000000.0 + 4.0 * h).sqrt())
    } else {
        let delta = h - 240000.0;
        500.0 + 2.0 * delta / (1400.0 + (1960000.0 + 6.0 * delta).sqrt())
    }
}

fn conductivity(t: f64) -> (f64, f64) {
    if t < 425.0 {
        (10.0 + 0.04 * (t - 300.0), 0.04)
    } else if t < 550.0 {
        (15.0 - 0.024 * (t - 425.0), -0.024)
    } else {
        (12.0 + 0.04 * (t - 550.0), 0.04)
    }
}

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "actual={actual:.17e}, expected={expected:.17e}, tolerance={tolerance:.3e}"
    );
}

#[test]
fn source_exact_state_and_tangents_preserve_independent_cp_and_conductivity_knots() {
    let fixture = fixture();
    let caloric = fixture.caloric();
    let law = IsobaricConductionLaw::try_new(
        &fixture.card,
        &caloric,
        &fixture.conductivity,
        SELECTION,
    )
    .unwrap();
    assert!(core::ptr::eq(law.caloric(), &caloric));
    assert_eq!(
        law.source_line()
            .knots
            .iter()
            .map(|(t, _)| *t)
            .collect::<Vec<_>>(),
        vec![300.0, 425.0, 550.0, 700.0],
    );
    assert_eq!(law.enthalpy_bounds_j_kg(), [0.0, 580000.0]);
    for t in [300.0, 317.0, 425.0, 499.0, 500.0, 550.0, 610.0, 700.0] {
        let h = enthalpy(t);
        let cp = if t <= 500.0 {
            1000.0 + 2.0 * (t - 300.0)
        } else {
            1400.0 + 3.0 * (t - 500.0)
        };
        let (k, dk_dt) = conductivity(t);
        let state = law.state_at_specific_enthalpy(h).unwrap();
        close(state.temperature_k, t, 2e-12);
        close(state.temperature_slope_kg_k_j, 1.0 / cp, 2e-18);
        close(state.conductivity_w_m_k, k, 2e-13);
        close(state.conductivity_slope_kg_w_m_k_j, dk_dt / cp, 2e-18);
        close(law.conductivity_w_m_k(h).unwrap(), k, 2e-13);
        let observables = law.time_observables(h).unwrap();
        close(observables.temperature_k, t, 2e-12);
        assert_eq!(observables.phase_mass_fraction, None);
    }
    for h in [-1.0, 580001.0, f64::NAN, f64::INFINITY] {
        assert!(law.state_at_specific_enthalpy(h).is_err());
    }
}

#[test]
fn analytic_tangents_match_constitutive_finite_differences_away_from_corners() {
    let fixture = fixture();
    let caloric = fixture.caloric();
    let law = IsobaricConductionLaw::try_new(
        &fixture.card,
        &caloric,
        &fixture.conductivity,
        SELECTION,
    )
    .unwrap();
    for t in [350.0, 470.0, 620.0] {
        let h = enthalpy(t);
        let delta = 0.1;
        let state = law.state_at_specific_enthalpy(h).unwrap();
        let lower = law.state_at_specific_enthalpy(h - delta).unwrap();
        let upper = law.state_at_specific_enthalpy(h + delta).unwrap();
        close(
            state.temperature_slope_kg_k_j,
            (upper.temperature_k - lower.temperature_k) / (2.0 * delta),
            1e-12,
        );
        close(
            state.conductivity_slope_kg_w_m_k_j,
            (upper.conductivity_w_m_k - lower.conductivity_w_m_k) / (2.0 * delta),
            5e-14,
        );
    }
}

#[test]
fn source_binding_rejects_changed_card_pressure_and_nonpositive_conductivity() {
    let original = fixture();
    let caloric = original.caloric();
    let changed = Fixture::new(Some(100000.0), [11.0, 15.0, 12.0, 18.0]);
    assert!(
        IsobaricConductionLaw::try_new(
            &changed.card,
            &caloric,
            &changed.conductivity,
            SELECTION,
        )
        .is_err()
    );
    for bad in [
        Fixture::new(None, [10.0, 15.0, 12.0, 18.0]),
        Fixture::new(Some(200000.0), [10.0, 15.0, 12.0, 18.0]),
        Fixture::new(Some(100000.0), [10.0, 0.0, 12.0, 18.0]),
    ] {
        let caloric = bad.caloric();
        assert!(
            IsobaricConductionLaw::try_new(
                &bad.card,
                &caloric,
                &bad.conductivity,
                SELECTION,
            )
            .is_err()
        );
    }
}

#[test]
fn native_implicit_exchange_matches_independent_nonlinear_heat_root_and_conserves_energy() {
    let fixture = fixture();
    let caloric = fixture.caloric();
    let law = IsobaricConductionLaw::try_new(
        &fixture.card,
        &caloric,
        &fixture.conductivity,
        SELECTION,
    )
    .unwrap();
    let material = MaterialId(1);
    let boundary = Boundary {
        ambient_k: 300.0,
        htc_w_m2_k: 0.0,
        reservoir: None,
        source_w: 0.0,
    };
    let cell = |id, mass_kg, t| DepositedCell {
        id: CellId(id),
        segment_id: id,
        body_id: 1,
        material_id: material,
        mass_kg,
        specific_h_j_kg: enthalpy(t),
        geometry: CellGeometry {
            start_m: [0.2 * id as f64, 0.0, 0.0],
            end_m: [0.2 * (id + 1) as f64, 0.0, 0.0],
            volume_m3: mass_kg / 1000.0,
            surface_area_m2: 2.0,
        },
        boundary,
    };
    let schedule = vec![
        TimedEvent {
            time_s: 0.0,
            event: Event::Birth(cell(1, 1.0, 300.0)),
        },
        TimedEvent {
            time_s: 0.0,
            event: Event::Birth(cell(2, 2.0, 650.0)),
        },
        TimedEvent {
            time_s: 0.0,
            event: Event::Contact(ThermalLink {
                id: 1,
                a: CellId(1),
                b: CellId(2),
                area_m2: 1.0,
                conductance_w_k: 0.0,
            }),
        },
    ];
    let thermal = DepositionSimulation::new(
        schedule,
        &[(material, &caloric as &dyn ThermalLaw)],
        DepositionBudget {
            max_cells: 2,
            max_links: 1,
            max_events: 3,
            max_substeps: 32,
        },
    )
    .unwrap();
    let (thermal, _) = thermal
        .advance_to(
            0.0,
            StepConfig {
                max_dt_s: 1.0,
                energy_tolerance_j: 1e-7,
            },
            || Ok(()),
        )
        .unwrap();
    let conduction = ConductionObserver::new(
        &thermal,
        &[(material, &law as &dyn ConductivityLaw)],
        &[ConductiveLink {
            link_id: 1,
            half_distances_m: [0.1, 0.1],
            contact_resistance_m2_k_w: 0.0,
        }],
        &[],
        || Ok(()),
    )
    .unwrap();
    let implicit = ImplicitConductionChain::new(
        &thermal,
        &conduction,
        &[(material, &law as &dyn ImplicitConductionLaw)],
        || Ok(()),
    )
    .unwrap();
    let dt = 50.0;
    let (next, report) = implicit
        .advance_to(
            &thermal,
            dt,
            ImplicitConductionControls {
                maximum_dt_s: dt,
                cell_tolerance_j: 1e-7,
                total_tolerance_j: 1e-7,
                ..ImplicitConductionControls::default()
            },
            &mut ImplicitConductionWorkspace::new(2).unwrap(),
            || Ok(()),
        )
        .unwrap();

    // Backward Euler's heat Q satisfies one scalar equation independently of
    // the production Newton solver, caloric inverse and conductivity evaluator.
    let hot_initial_h = enthalpy(650.0);
    let residual = |q: f64| {
        let cold_t = temperature(q);
        let hot_t = temperature(hot_initial_h - q / 2.0);
        let g = 1.0 / (0.1 / conductivity(cold_t).0 + 0.1 / conductivity(hot_t).0);
        q - dt * g * (hot_t - cold_t)
    };
    let (mut low, mut high) = (0.0, hot_initial_h * 2.0 / 3.0);
    assert!(residual(low) < 0.0 && residual(high) > 0.0);
    for _ in 0..80 {
        let mid = (low + high) / 2.0;
        if residual(mid) < 0.0 {
            low = mid;
        } else {
            high = mid;
        }
    }
    let q = (low + high) / 2.0;
    let cold = &next.cells()[&CellId(1)];
    let hot = &next.cells()[&CellId(2)];
    close(cold.specific_h_j_kg, q, 2e-7);
    close(hot.specific_h_j_kg, hot_initial_h - q / 2.0, 2e-7);
    assert!(temperature(cold.specific_h_j_kg) > 500.0);
    close(
        cold.mass_kg * cold.specific_h_j_kg + hot.mass_kg * hot.specific_h_j_kg,
        2.0 * hot_initial_h,
        2e-7,
    );
    assert_eq!(report.energy.substeps, 1);
    assert_eq!(report.rejected_intervals, 0);
    assert!(report.energy.residual_j.abs() < 1e-7);
    assert!(report.entropy.thermal_storage_change_j_k > 0.0);
    assert!(report.entropy.conduction_generation_j_k > 0.0);
    // Advancing the returned candidate must not mutate the original timeline.
    close(thermal.cells()[&CellId(1)].specific_h_j_kg, 0.0, 0.0);
    close(
        thermal.cells()[&CellId(2)].specific_h_j_kg,
        hot_initial_h,
        0.0,
    );
}
