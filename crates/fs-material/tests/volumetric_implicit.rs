//! G0/G1/G3: original volumetric capacity on native implicit conduction.

use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSelection, ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim,
    PropertyKey, PropertyValue, Provenance, QueryPoint, SelectionPolicy, UncertaintyModel,
};
use fs_material::thermal::transport::volumetric::{
    PreparedVolumetricThermalTransport, VOLUMETRIC_HEAT_CAPACITY_DIMS, VolumetricThermalPathQuery,
    timeline::FixedVolumeThermalLaw,
};
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_qty::{Dims, QuantitySpec, Temperature};
use fs_thermal::enthalpy::{
    Boundary, CellGeometry, CellId, DepositedCell, DepositionBudget, DepositionError,
    DepositionSimulation, Event, MaterialId, StepConfig, ThermalLaw, ThermalLink,
    ThermalStorageBasis, TimedEvent,
    conduction::implicit::{
        ImplicitConductionChain, ImplicitConductionControls, ImplicitConductionGraph,
        ImplicitConductionLaw, ImplicitConductionLinearControls, ImplicitConductionWorkspace,
    },
    conduction::{ConductionObserver, ConductiveLink, ConductivityLaw},
};

fn path() -> PreparedVolumetricThermalTransport {
    path_with_capacity(vec![(300.0, 1000.0), (500.0, 1400.0), (700.0, 2000.0)])
}

fn path_with_capacity(capacity_knots: Vec<(f64, f64)>) -> PreparedVolumetricThermalTransport {
    let temperature = QuantitySpec::semantic(SemanticType::new(
        QuantityKind::AbsoluteTemperature,
        ValueForm::Static,
    ));
    let capacity = PropertyKey::new("volumetric-capacity", VOLUMETRIC_HEAT_CAPACITY_DIMS);
    let conductivity_dims = Dims([1, 1, -3, -1, 0, 0]);
    let conductivity = PropertyKey::new("directional-conductivity", conductivity_dims);
    let mut claims = ClaimSet::new();
    for (key, dims, knots) in [
        (
            capacity.clone(),
            VOLUMETRIC_HEAT_CAPACITY_DIMS,
            capacity_knots,
        ),
        (
            conductivity.clone(),
            conductivity_dims,
            vec![(300.0, 10.0), (425.0, 15.0), (550.0, 12.0), (700.0, 18.0)],
        ),
    ] {
        claims
            .insert_claim(PropertyClaim {
                key,
                value: PropertyValue::Curve {
                    abscissa: "T".into(),
                    abscissa_dims: Temperature::DIMS,
                    knots,
                    dims,
                },
                validity: ValidityDomain::unconstrained().with_quantity(
                    "T",
                    temperature,
                    300.0,
                    700.0,
                ),
                interpolation: InterpolationPolicy::LinearInside,
                uncertainty: UncertaintyModel::Unstated,
                observations: vec![],
                provenance: Provenance {
                    source: "synthetic independent Cvol and conductivity curves".into(),
                    license: "test-only".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic fixed-volume body".into(),
            phase: "solid".into(),
            process: "declared sensible storage".into(),
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
    };
    PreparedVolumetricThermalTransport::try_new(
        &card,
        VolumetricThermalPathQuery {
            volumetric_heat_capacity: &capacity,
            conductivity: &conductivity,
            temperature_axis: "T",
            lower: &point(300.0),
            upper: &point(700.0),
            selections: [ClaimSelection::Policy(SelectionPolicy::SingleClaimOnly); 2],
        },
    )
    .unwrap()
}

// Independent primitive/inverse, in J/m3 with datum at 300 K.
fn energy(t: f64) -> f64 {
    if t <= 500.0 {
        let x = t - 300.0;
        1000.0 * x + x * x
    } else {
        let x = t - 500.0;
        240000.0 + 1400.0 * x + 1.5 * x * x
    }
}
fn temperature(e: f64) -> f64 {
    if e <= 240000.0 {
        300.0 + 2.0 * e / (1000.0 + (1000000.0 + 4.0 * e).sqrt())
    } else {
        let x = e - 240000.0;
        500.0 + 2.0 * x / (1400.0 + (1960000.0 + 6.0 * x).sqrt())
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
fn native_volumetric_tangents_keep_source_knots_and_storage_basis() {
    let path = path();
    let law = FixedVolumeThermalLaw::try_new(&path, 1.0, 300.0).unwrap();
    assert_eq!(
        law.storage_basis(),
        ThermalStorageBasis::DeclaredSensibleStorage
    );
    assert!(std::ptr::addr_eq(
        law.caloric_law() as *const dyn ThermalLaw,
        &law as *const FixedVolumeThermalLaw<'_>
    ));
    assert_eq!(law.enthalpy_bounds_j_kg(), [0.0, 580000.0]);
    assert_eq!(
        path.temperatures_k().collect::<Vec<_>>(),
        vec![300.0, 425.0, 500.0, 550.0, 700.0]
    );
    for t in [300.0, 350.0, 425.0, 470.0, 500.0, 550.0, 620.0, 700.0] {
        let state = law.conduction_state_at_specific_storage(energy(t)).unwrap();
        let capacity = if t <= 500.0 {
            1000.0 + 2.0 * (t - 300.0)
        } else {
            1400.0 + 3.0 * (t - 500.0)
        };
        let (k, dk) = conductivity(t);
        close(state.temperature_k, t, 2e-12);
        close(state.temperature_slope_kg_k_j, 1.0 / capacity, 2e-18);
        close(state.conductivity_w_m_k, k, 2e-13);
        close(state.conductivity_slope_kg_w_m_k_j, dk / capacity, 2e-18);
        close(law.conductivity_w_m_k(energy(t)).unwrap(), k, 2e-13);
        assert_eq!(
            law.time_observables(energy(t)).unwrap().phase_mass_fraction,
            None
        );
    }
    for t in [350.0, 470.0, 620.0] {
        let h = energy(t);
        let center = law.conduction_state_at_specific_storage(h).unwrap();
        let a = law.conduction_state_at_specific_storage(h - 0.1).unwrap();
        let b = law.conduction_state_at_specific_storage(h + 0.1).unwrap();
        close(
            center.temperature_slope_kg_k_j,
            (b.temperature_k - a.temperature_k) / 0.2,
            1e-12,
        );
        close(
            center.conductivity_slope_kg_w_m_k_j,
            (b.conductivity_w_m_k - a.conductivity_w_m_k) / 0.2,
            5e-14,
        );
    }
    for h in [-1.0, 580001.0, f64::NAN, f64::INFINITY] {
        assert!(law.conduction_state_at_specific_storage(h).is_err());
    }
}

#[test]
fn capacity_only_breakpoints_cannot_change_the_original_conductivity_tangent() {
    let next = 300.0_f64.next_up();
    let path = path_with_capacity(vec![(300.0, 1000.0), (next, 1000.0), (700.0, 2000.0)]);
    let law = FixedVolumeThermalLaw::try_new(&path, 1.0, 300.0).unwrap();
    let state = law.conduction_state_at_specific_storage(0.0).unwrap();
    assert_eq!(state.conductivity_w_m_k, 10.0);
    close(state.conductivity_slope_kg_w_m_k_j, 0.04 / 1000.0, 1e-19);
}

fn cell(law: &FixedVolumeThermalLaw<'_>, id: u64, volume: f64, t: f64) -> DepositedCell {
    DepositedCell {
        id: CellId(id),
        segment_id: id,
        body_id: 1,
        material_id: MaterialId(1),
        mass_kg: law.reference_mass_kg(volume).unwrap(),
        specific_h_j_kg: law.specific_storage_j_kg(t).unwrap(),
        geometry: CellGeometry {
            start_m: [id as f64, 0.0, 0.0],
            end_m: [id as f64 + 1.0, 0.0, 0.0],
            volume_m3: volume,
            surface_area_m2: 2.0,
        },
        boundary: Boundary {
            ambient_k: 300.0,
            htc_w_m2_k: 0.0,
            reservoir: None,
            source_w: 0.0,
        },
    }
}

#[test]
fn native_chain_and_graph_conserve_original_volumetric_storage_across_capacity_knots() {
    let path = path();
    let dt = 50.0;
    let residual = |q: f64| {
        let cold = temperature(q);
        let hot = temperature(energy(650.0) - q / 2.0);
        let g = 1.0 / (0.1 / conductivity(cold).0 + 0.1 / conductivity(hot).0);
        q - dt * g * (hot - cold)
    };
    let (mut low, mut high) = (0.0, 2.0 * energy(650.0) / 3.0);
    assert!(residual(low) < 0.0 && residual(high) > 0.0);
    for _ in 0..80 {
        let mid = 0.5 * (low + high);
        if residual(mid) < 0.0 {
            low = mid;
        } else {
            high = mid;
        }
    }
    let expected_heat = 0.5 * (low + high);
    for (density, datum, graph) in [(2.0, 300.0, false), (100.0, 450.0, true)] {
        let law = FixedVolumeThermalLaw::try_new(&path, density, datum)
            .unwrap()
            .with_sensible_entropy();
        let initial = [cell(&law, 1, 1.0, 300.0), cell(&law, 2, 2.0, 650.0)];
        let mut wrong_mass = initial[0].clone();
        wrong_mass.mass_kg *= 1.01;
        assert!(law.validate_cell(&wrong_mass).is_err());
        let thermal = DepositionSimulation::new(
            vec![
                TimedEvent {
                    time_s: 0.0,
                    event: Event::Birth(initial[0].clone()),
                },
                TimedEvent {
                    time_s: 0.0,
                    event: Event::Birth(initial[1].clone()),
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
            ],
            &[(MaterialId(1), &law as &dyn ThermalLaw)],
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
            &[(MaterialId(1), &law as &dyn ConductivityLaw)],
            &[ConductiveLink {
                link_id: 1,
                half_distances_m: [0.1, 0.1],
                contact_resistance_m2_k_w: 0.0,
            }],
            &[],
            || Ok(()),
        )
        .unwrap();
        let laws = [(MaterialId(1), &law as &dyn ImplicitConductionLaw)];
        let controls = ImplicitConductionControls {
            maximum_dt_s: dt,
            cell_tolerance_j: 1e-7,
            total_tolerance_j: 1e-7,
            ..ImplicitConductionControls::default()
        };
        let mut workspace = ImplicitConductionWorkspace::new(2).unwrap();
        let (next, report) = if graph {
            let implicit = ImplicitConductionGraph::new(
                &thermal,
                &conduction,
                &laws,
                ImplicitConductionLinearControls::default(),
                || Ok(()),
            )
            .unwrap();
            implicit
                .advance_to(&thermal, dt, controls, &mut workspace, || Ok(()))
                .unwrap()
        } else {
            let implicit =
                ImplicitConductionChain::new(&thermal, &conduction, &laws, || Ok(())).unwrap();
            assert!(
                implicit
                    .advance_to(&thermal, dt, controls, &mut workspace, || Err(
                        DepositionError::Cancelled
                    ))
                    .is_err()
            );
            implicit
                .advance_to(&thermal, dt, controls, &mut workspace, || Ok(()))
                .unwrap()
        };
        let cold = &next.cells()[&CellId(1)];
        let hot = &next.cells()[&CellId(2)];
        close(
            cold.mass_kg * (cold.specific_h_j_kg - initial[0].specific_h_j_kg),
            expected_heat,
            2e-7,
        );
        close(
            hot.mass_kg * (hot.specific_h_j_kg - initial[1].specific_h_j_kg),
            -expected_heat,
            2e-7,
        );
        close(
            law.temperature_k(cold.specific_h_j_kg).unwrap(),
            temperature(expected_heat),
            2e-9,
        );
        assert!(law.temperature_k(cold.specific_h_j_kg).unwrap() > 500.0);
        assert!(report.energy.residual_j.abs() < 1e-7);
        assert_eq!(report.energy.substeps, 1);
        assert_eq!(report.rejected_intervals, 0);
        assert!(report.entropy.thermal_storage_change_j_k > 0.0);
        assert!(report.entropy.conduction_generation_j_k > 0.0);
        assert_eq!(thermal.cells()[&CellId(1)], initial[0]);
        assert_eq!(thermal.cells()[&CellId(2)], initial[1]);
    }
}
