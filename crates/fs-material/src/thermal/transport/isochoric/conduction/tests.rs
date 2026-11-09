use super::*;
use crate::thermal::transport::{
    CONDUCTIVITY, DENSITY, SPECIFIC_HEAT, isochoric::IsochoricHeatCapacityQuery,
};
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialStateId, PropertyClaim, PropertyValue, Provenance,
    QueryPoint, SelectionPolicy, UncertaintyModel,
};
use fs_qty::{
    Pressure, QuantitySpec, Temperature,
    semantic::{QuantityKind, SemanticType, ValueForm},
};
use fs_thermal::enthalpy::{
    Boundary, CellGeometry, CellId, DepositedCell, DepositionBudget, DepositionSimulation, Event,
    MaterialId, StepConfig, ThermalLink, ThermalStorageBasis, TimedEvent,
    conduction::implicit::{
        ImplicitConductionChain, ImplicitConductionControls, ImplicitConductionGraph,
        ImplicitConductionLinearControls, ImplicitConductionWorkspace,
        adaptive::AdaptiveImplicitConductionControls,
    },
    conduction::{ConductionObserver, ConductiveLink},
};

const RHO: f64 = 1000.0;
const AREA: f64 = 0.01;
const MATERIAL: MaterialId = MaterialId(1);
const SELECTION: ClaimSelection = ClaimSelection::Policy(SelectionPolicy::SingleClaimOnly);

struct Fixture {
    card: MaterialCard,
    cv: PropertyKey,
    k: PropertyKey,
    lower: QueryPoint,
    upper: QueryPoint,
}

impl Fixture {
    fn new(
        capacities: [f64; 3],
        conductivities: [f64; 4],
        density: Option<f64>,
        pressure: bool,
    ) -> Self {
        let temperature = QuantitySpec::semantic(SemanticType::new(
            QuantityKind::AbsoluteTemperature,
            ValueForm::Static,
        ));
        let rho = QuantitySpec::dimensional(DENSITY);
        let cv = PropertyKey::new("isochoric_cv", SPECIFIC_HEAT);
        let k = PropertyKey::new("scalar_conductivity", CONDUCTIVITY);
        let validity =
            || ValidityDomain::unconstrained().with_quantity("T", temperature, 300.0, 700.0);
        let claim = |key: PropertyKey, dims, knots, validity| PropertyClaim {
            key,
            value: PropertyValue::Curve {
                abscissa: "T".into(),
                abscissa_dims: Temperature::DIMS,
                knots,
                dims,
            },
            validity,
            interpolation: InterpolationPolicy::LinearInside,
            uncertainty: UncertaintyModel::Unstated,
            observations: vec![],
            provenance: Provenance {
                source: "synthetic fixed-density Cv and conductivity".into(),
                license: "test-only".into(),
                artifact: None,
            },
        };
        let mut claims = ClaimSet::new();
        claims
            .insert_claim(claim(
                cv.clone(),
                SPECIFIC_HEAT,
                [300.0, 500.0, 700.0].into_iter().zip(capacities).collect(),
                validity().with_quantity("rho", rho, RHO, RHO),
            ))
            .unwrap();
        let mut k_validity = validity();
        if let Some(value) = density {
            k_validity = k_validity.with_quantity("rho", rho, value, value);
        }
        if pressure {
            k_validity = k_validity.with_quantity(
                "p",
                QuantitySpec::dimensional(Pressure::DIMS),
                100000.0,
                100000.0,
            );
        }
        claims
            .insert_claim(claim(
                k.clone(),
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
                process: "fixed-volume numerical fixture".into(),
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
                .with_quantity("rho", rho, RHO)
                .unwrap()
        };
        Self {
            card,
            cv,
            k,
            lower: point(300.0),
            upper: point(700.0),
        }
    }

    fn nonlinear() -> Self {
        Self::new(
            [500.0, 700.0, 1000.0],
            [10.0, 15.0, 12.0, 18.0],
            Some(RHO),
            false,
        )
    }

    fn caloric(&self) -> IsochoricHeatCapacityLaw {
        IsochoricHeatCapacityLaw::try_new(
            &self.card,
            IsochoricHeatCapacityQuery {
                specific_heat_cv: &self.cv,
                temperature_axis: "T",
                density_axis: "rho",
                lower: &self.lower,
                upper: &self.upper,
                selection: SELECTION,
            },
            300.0,
        )
        .unwrap()
    }

    fn conduction<'a>(&self, caloric: &'a IsochoricHeatCapacityLaw) -> IsochoricConductionLaw<'a> {
        IsochoricConductionLaw::try_new(&self.card, caloric, &self.k, SELECTION).unwrap()
    }
}

// Independent piecewise quadratic primitive/inverse for the nonlinear fixture.
fn energy(t: f64) -> f64 {
    if t <= 500.0 {
        let x = t - 300.0;
        500.0 * x + 0.5 * x * x
    } else {
        let x = t - 500.0;
        120000.0 + 700.0 * x + 0.75 * x * x
    }
}

fn temperature(u: f64) -> f64 {
    if u <= 120000.0 {
        300.0 + 2.0 * u / (500.0 + (250000.0 + 2.0 * u).sqrt())
    } else {
        let delta = u - 120000.0;
        500.0 + 2.0 * delta / (700.0 + (490000.0 + 3.0 * delta).sqrt())
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

fn cell(law: &IsochoricHeatCapacityLaw, id: usize, mass_kg: f64, t: f64) -> DepositedCell {
    let length = mass_kg / (RHO * AREA);
    DepositedCell {
        id: CellId(id as u64),
        segment_id: id as u64,
        body_id: 1,
        material_id: MATERIAL,
        mass_kg,
        specific_h_j_kg: law.specific_internal_energy_j_kg(t).unwrap(),
        geometry: CellGeometry {
            start_m: [id as f64, 0.0, 0.0],
            end_m: [id as f64 + length, 0.0, 0.0],
            volume_m3: mass_kg / RHO,
            surface_area_m2: 2.0 * AREA + 4.0 * AREA.sqrt() * length,
        },
        boundary: Boundary {
            ambient_k: 300.0,
            htc_w_m2_k: 0.0,
            reservoir: None,
            source_w: 0.0,
        },
    }
}

fn budget() -> DepositionBudget {
    DepositionBudget {
        max_cells: 8,
        max_links: 8,
        max_events: 16,
        max_substeps: 4096,
    }
}

fn timeline<'a>(
    caloric: &'a IsochoricHeatCapacityLaw,
    values: &[(f64, f64)],
    edges: &[(usize, usize)],
) -> (DepositionSimulation<'a>, Vec<ConductiveLink>) {
    let mut schedule: Vec<_> = values
        .iter()
        .enumerate()
        .map(|(i, &(mass, t))| TimedEvent {
            time_s: 0.0,
            event: Event::Birth(cell(caloric, i, mass, t)),
        })
        .collect();
    let mut links = Vec::new();
    for (i, &(a, b)) in edges.iter().enumerate() {
        schedule.push(TimedEvent {
            time_s: 0.0,
            event: Event::Contact(ThermalLink {
                id: i as u64,
                a: CellId(a as u64),
                b: CellId(b as u64),
                area_m2: AREA,
                conductance_w_k: 0.0,
            }),
        });
        links.push(ConductiveLink {
            link_id: i as u64,
            half_distances_m: [
                values[a].0 / (2.0 * RHO * AREA),
                values[b].0 / (2.0 * RHO * AREA),
            ],
            contact_resistance_m2_k_w: 0.0,
        });
    }
    let thermal = DepositionSimulation::new(
        schedule,
        &[(MATERIAL, caloric as &dyn ThermalLaw)],
        budget(),
    )
    .unwrap()
    .advance_to(
        0.0,
        StepConfig {
            max_dt_s: 1.0,
            energy_tolerance_j: 1e-7,
        },
        || Ok(()),
    )
    .unwrap()
    .0;
    (thermal, links)
}

#[test]
fn exact_cv_state_retains_internal_energy_and_independent_conductivity_knots() {
    let fixture = Fixture::nonlinear();
    let caloric = fixture.caloric();
    let law = fixture.conduction(&caloric);
    assert!(core::ptr::eq(law.caloric(), &caloric));
    assert_eq!(
        law.caloric_law().storage_basis(),
        ThermalStorageBasis::InternalEnergy
    );
    assert_eq!(law.enthalpy_bounds_j_kg(), [0.0, 290000.0]);
    assert_eq!(
        law.source_line()
            .knots
            .iter()
            .map(|(t, _)| *t)
            .collect::<Vec<_>>(),
        vec![300.0, 425.0, 550.0, 700.0]
    );
    for t in [300.0, 317.0, 425.0, 499.0, 500.0, 550.0, 610.0, 700.0] {
        let u = energy(t);
        let cv = if t <= 500.0 {
            500.0 + t - 300.0
        } else {
            700.0 + 1.5 * (t - 500.0)
        };
        let (k, dk_dt) = conductivity(t);
        let state = law.state_at_specific_internal_energy(u).unwrap();
        close(state.temperature_k, t, 2e-12);
        close(state.temperature_slope_kg_k_j, 1.0 / cv, 3e-18);
        close(state.conductivity_w_m_k, k, 2e-13);
        close(state.conductivity_slope_kg_w_m_k_j, dk_dt / cv, 3e-18);
        assert_eq!(law.state_at_specific_enthalpy(u).unwrap(), state);
        close(law.conductivity_w_m_k(u).unwrap(), k, 2e-13);
        assert_eq!(law.time_observables(u).unwrap().phase_mass_fraction, None);
    }
    for t in [350.0, 470.0, 620.0] {
        let u = energy(t);
        let state = law.state_at_specific_internal_energy(u).unwrap();
        let a = law.state_at_specific_internal_energy(u - 0.1).unwrap();
        let b = law.state_at_specific_internal_energy(u + 0.1).unwrap();
        close(
            state.temperature_slope_kg_k_j,
            (b.temperature_k - a.temperature_k) / 0.2,
            1e-12,
        );
        close(
            state.conductivity_slope_kg_w_m_k_j,
            (b.conductivity_w_m_k - a.conductivity_w_m_k) / 0.2,
            5e-14,
        );
    }
    for u in [-1.0, 290001.0, f64::NAN, f64::INFINITY] {
        assert!(law.state_at_specific_internal_energy(u).is_err());
    }
}

#[test]
fn original_card_density_and_cell_mass_volume_cannot_be_substituted() {
    let fixture = Fixture::nonlinear();
    let caloric = fixture.caloric();
    let changed = Fixture::new(
        [500.0, 700.0, 1000.0],
        [11.0, 15.0, 12.0, 18.0],
        Some(RHO),
        false,
    );
    assert!(
        IsochoricConductionLaw::try_new(&changed.card, &caloric, &changed.k, SELECTION).is_err()
    );
    for (rho, pressure, k) in [
        (None, false, [10.0, 15.0, 12.0, 18.0]),
        (Some(2.0 * RHO), false, [10.0, 15.0, 12.0, 18.0]),
        (Some(RHO), true, [10.0, 15.0, 12.0, 18.0]),
        (Some(RHO), false, [10.0, 0.0, 12.0, 18.0]),
    ] {
        let bad = Fixture::new([500.0, 700.0, 1000.0], k, rho, pressure);
        let cv = bad.caloric();
        assert!(IsochoricConductionLaw::try_new(&bad.card, &cv, &bad.k, SELECTION).is_err());
    }
    let law = fixture.conduction(&caloric);
    let mut wrong_cell = cell(&caloric, 0, 1.0, 450.0);
    wrong_cell.geometry.volume_m3 *= 1.01;
    assert!(law.caloric_law().validate_cell(&wrong_cell).is_err());
    let thermal = DepositionSimulation::new(
        vec![TimedEvent {
            time_s: 0.0,
            event: Event::Birth(wrong_cell),
        }],
        &[(MATERIAL, law.caloric_law())],
        budget(),
    )
    .unwrap();
    assert!(
        thermal
            .advance_to(
                0.0,
                StepConfig {
                    max_dt_s: 1.0,
                    energy_tolerance_j: 1e-7
                },
                || Ok(()),
            )
            .is_err()
    );
    assert!(thermal.cells().is_empty());
}

#[test]
fn nonlinear_cv_exchange_matches_independent_heat_root_and_closes_internal_energy() {
    let fixture = Fixture::nonlinear();
    let caloric = fixture.caloric();
    let law = fixture.conduction(&caloric);
    let (thermal, links) = timeline(&caloric, &[(1.0, 300.0), (2.0, 650.0)], &[(0, 1)]);
    let conduction = ConductionObserver::new(
        &thermal,
        &[(MATERIAL, &law as &dyn ConductivityLaw)],
        &links,
        &[],
        || Ok(()),
    )
    .unwrap();
    let implicit = ImplicitConductionChain::new(
        &thermal,
        &conduction,
        &[(MATERIAL, &law as &dyn ImplicitConductionLaw)],
        || Ok(()),
    )
    .unwrap();
    let dt = 2000.0;
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
    let initial_hot = energy(650.0);
    let residual = |q: f64| {
        let cold_t = temperature(q);
        let hot_t = temperature(initial_hot - q / 2.0);
        let g = AREA / (0.05 / conductivity(cold_t).0 + 0.1 / conductivity(hot_t).0);
        q - dt * g * (hot_t - cold_t)
    };
    let (mut low, mut high) = (0.0, initial_hot * 2.0 / 3.0);
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
    let cold_u = next.cells()[&CellId(0)].specific_h_j_kg;
    let hot_u = next.cells()[&CellId(1)].specific_h_j_kg;
    close(cold_u, q, 2e-7);
    close(hot_u, initial_hot - q / 2.0, 2e-7);
    assert!(temperature(cold_u) > 500.0); // crosses Cv and k breakpoints
    close(cold_u + 2.0 * hot_u, 2.0 * initial_hot, 2e-7);
    assert_eq!(
        next.material_storage_basis(MATERIAL).unwrap(),
        ThermalStorageBasis::InternalEnergy
    );
    assert_eq!(report.energy.substeps, 1);
    assert!(report.energy.residual_j.abs() <= 1e-7);
    assert!(report.entropy.thermal_storage_change_j_k > 0.0);
    assert!(report.entropy.conduction_generation_j_k > 0.0);
    close(thermal.cells()[&CellId(0)].specific_h_j_kg, 0.0, 0.0);
    close(
        thermal.cells()[&CellId(1)].specific_h_j_kg,
        initial_hot,
        1e-10,
    );
}

#[test]
fn branched_cv_graph_uses_native_adaptive_sdirk2_and_matches_analytic_decay() {
    let fixture = Fixture::new([500.0; 3], [10.0; 4], Some(RHO), false);
    let caloric = fixture.caloric();
    let law = fixture.conduction(&caloric);
    let (thermal, links) = timeline(
        &caloric,
        &[(1.0, 300.0), (1.0, 600.0), (1.0, 600.0), (1.0, 600.0)],
        &[(0, 1), (0, 2), (0, 3)],
    );
    let conduction = ConductionObserver::new(
        &thermal,
        &[(MATERIAL, &law as &dyn ConductivityLaw)],
        &links,
        &[],
        || Ok(()),
    )
    .unwrap();
    let graph = ImplicitConductionGraph::new(
        &thermal,
        &conduction,
        &[(MATERIAL, &law as &dyn ImplicitConductionLaw)],
        ImplicitConductionLinearControls::default(),
        || Ok(()),
    )
    .unwrap();
    let time = 125.0;
    let (next, report) = graph
        .advance_sdirk2_adaptive_to(
            &thermal,
            time,
            AdaptiveImplicitConductionControls {
                initial_dt_s: time,
                step: ImplicitConductionControls {
                    maximum_dt_s: time,
                    cell_tolerance_j: 1e-7,
                    total_tolerance_j: 1e-7,
                    ..ImplicitConductionControls::default()
                },
                enthalpy_absolute_tolerance_j_kg: 5e-4,
                temperature_absolute_tolerance_k: 1e-6,
                relative_tolerance: 0.0,
                ..AdaptiveImplicitConductionControls::default()
            },
            &mut ImplicitConductionWorkspace::new(4).unwrap(),
            || Ok(()),
        )
        .unwrap();
    // Each edge has G = A*k/(0.05+0.05) = 1 W/K. Symmetry leaves one
    // decay mode with rate 4*G/(m*Cv), around the conserved mean 525 K.
    let decay = (-4.0 * time / 500.0).exp();
    close(
        next.temperature_k(CellId(0)).unwrap(),
        525.0 - 225.0 * decay,
        1e-3,
    );
    for id in [CellId(1), CellId(2), CellId(3)] {
        close(next.temperature_k(id).unwrap(), 525.0 + 75.0 * decay, 1e-3);
    }
    close(next.energy_j(), thermal.energy_j(), 1e-6);
    assert_eq!(
        next.material_storage_basis(MATERIAL).unwrap(),
        ThermalStorageBasis::InternalEnergy
    );
    assert!(report.accepted_comparisons > 0);
    assert!(report.solve.linear_iterations > 0);
    assert!(report.maximum_accepted_normalized_error <= 1.0);
    assert!(report.solve.entropy.thermal_storage_change_j_k > 0.0);
}

#[test]
fn g1_direct_temperature_conductivity_retains_source_values_and_support() {
    let fixture = Fixture::nonlinear();
    let caloric = fixture.caloric();
    let law = fixture.conduction(&caloric);
    for t in [300.0, 350.0, 425.0, 500.0, 550.0, 625.0, 700.0] {
        let direct = law.conductivity_at_temperature_w_m_k(t).unwrap();
        close(direct, conductivity(t).0, 2e-13);
        let u = caloric.specific_internal_energy_j_kg(t).unwrap();
        close(
            direct,
            law.state_at_specific_internal_energy(u)
                .unwrap()
                .conductivity_w_m_k,
            2e-13,
        );
    }
    for t in [
        300.0_f64.next_down(),
        700.0_f64.next_up(),
        f64::NAN,
        f64::INFINITY,
    ] {
        assert!(law.conductivity_at_temperature_w_m_k(t).is_err());
    }
}
