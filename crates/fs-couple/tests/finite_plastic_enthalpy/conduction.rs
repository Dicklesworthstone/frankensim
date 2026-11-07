//! Actual normal transport follows the finite mechanical proposal. Material
//! values are synthetic; the scalar thermal recurrence is independent of the
//! coupled implementation. The parent point-law fixture supplies forming work.

use super::*;
use fs_solid::hyper3::plastic::thermal::FiniteThermoplasticTetSolid;
use fs_thermal::enthalpy::conduction::{
    ConductionObserver, ConductiveReservoir, ConductivityLaw, CurrentConductionGeometry,
    CurrentConductiveBoundary,
};

const AREA: f64 = 0.5 * L * L;
const DISTANCE: f64 = 0.25 * L;
const CONDUCTIVITY: f64 = 0.5;
const BATH: f64 = 340.0;

struct SyntheticConductivity;
impl ConductivityLaw for SyntheticConductivity {
    fn conductivity_w_m_k(&self, h: f64) -> Result<f64, DepositionError> {
        if h.is_finite() && (-20.0 * CP..=120.0 * CP).contains(&h) {
            Ok(CONDUCTIVITY)
        } else {
            Err(DepositionError::MaterialDomain)
        }
    }
}
static CONDUCTOR: SyntheticConductivity = SyntheticConductivity;

fn normal_metric(
    solid: &FiniteThermoplasticTetSolid<'_, FiniteThermoplasticMaterial<'_>>,
    cx: &Cx<'_>,
) -> Result<CurrentConductionGeometry, DepositionError> {
    cx.checkpoint().map_err(|_| DepositionError::Cancelled)?;
    let f = solid.element_deformation()[0];
    // This authored one-cell reduction has only a face normal to reference x.
    // Diagonal stretches preserve its normal path; a skewed path needs another
    // transport discretization instead of silently applying this formula.
    if [1, 2, 3, 5, 6, 7].into_iter().any(|i| f[i].abs() > 1.0e-12)
        || [0, 4, 8]
            .into_iter()
            .any(|i| !f[i].is_finite() || f[i] <= 0.0)
    {
        return Err(DepositionError::InvalidInput(
            "fixture needs positive normal stretches",
        ));
    }
    Ok(CurrentConductionGeometry {
        links: BTreeMap::new(),
        reservoirs: [(
            CELL,
            CurrentConductiveBoundary {
                area_m2: AREA * f[4] * f[8],
                half_distance_m: DISTANCE * f[0],
            },
        )]
        .into(),
        ambient_faces: BTreeMap::new(),
    })
}

fn source_power(expected: &Endpoint) -> f64 {
    let f = expected.deformation;
    let g = CONDUCTIVITY * AREA * f[4] * f[8] / (DISTANCE * f[0]);
    let capacity = RHO * V * CP;
    // Four actual explicit quarter-second steps with endpoint geometry held.
    // T4=Tb+r^4*(T0-Tb)+(P_total/G)*(1-r^4), with r=1-G*dt/C.
    let r4 = (1.0 - g * 0.25 / capacity).powi(4);
    let total_w =
        g * (expected.temperature_k - BATH - r4 * (expected.old_temperature_k - BATH)) / (1.0 - r4);
    total_w - V * expected.increment.caloric_heat_j_m3()
}

fn bath_boundary(expected: &Endpoint) -> Boundary {
    Boundary {
        ambient_k: 300.0,
        htc_w_m2_k: 0.0,
        reservoir: Some(Reservoir {
            area_m2: AREA,
            conductance_w_k: 0.0,
            temperature_k: BATH,
        }),
        source_w: source_power(expected),
    }
}

fn bath_timeline<'a>(fixture: &'a Fixture, endpoints: &[Endpoint]) -> DepositionSimulation<'a> {
    let mut events = vec![TimedEvent {
        time_s: 0.0,
        event: Event::Birth(DepositedCell {
            id: CELL,
            material_id: MaterialId(0),
            segment_id: 1,
            body_id: 7,
            mass_kg: RHO * V,
            specific_h_j_kg: 0.0,
            geometry: CellGeometry {
                start_m: [0.0; 3],
                end_m: [L; 3],
                volume_m3: V,
                surface_area_m2: 3.0 * L * L,
            },
            boundary: bath_boundary(&endpoints[0]),
        }),
    }];
    for (i, expected) in endpoints.iter().enumerate().skip(1) {
        events.push(TimedEvent {
            time_s: i as f64,
            event: Event::Boundary(CELL, bath_boundary(expected)),
        });
    }
    DepositionSimulation::new(
        events,
        &[(MaterialId(0), fixture.caloric.curve() as &dyn ThermalLaw)],
        DepositionBudget {
            max_cells: 1,
            max_links: 0,
            max_events: 8,
            max_substeps: 1000,
        },
    )
    .unwrap()
    .advance_to(0.0, thermal_config(), || Ok(()))
    .unwrap()
    .0
}

#[test]
fn solved_current_metrics_close_load_heat_reversal_and_geometry_refusal_retry() {
    context(|cx| {
        let fixture = Fixture::new(CP, 703);
        let endpoints = path(&fixture, false);
        let mut thermal = bath_timeline(&fixture, &endpoints);
        let mut solid = bind(&fixture, &thermal, false, cx);
        let mut conduction = ConductionObserver::new(
            &thermal,
            &[(MaterialId(0), &CONDUCTOR)],
            &[],
            &[ConductiveReservoir {
                cell_id: CELL,
                half_distance_m: DISTANCE,
                surface_resistance_m2_k_w: 0.0,
            }],
            || Ok(()),
        )
        .unwrap();
        let first = &endpoints[0];
        let accepted_cells = thermal.cells().clone();
        let accepted_history = solid.mechanics().element_states().to_vec();
        // Loose thermal gates cannot conceal a geometrically inconsistent
        // endpoint. The strict geometry gate refuses this first proposal.
        let refusal = solid.advance_with_deformed_conduction(
            &thermal,
            &conduction,
            first.time_s,
            &first.prescribed(),
            &[(3, first.free_force_n)],
            [0.0; 3],
            settings(),
            controls(),
            thermal_config(),
            FiniteThermoplasticConfig {
                maximum_iterations: 1,
                temperature_tolerance_k: 100.0,
                energy_tolerance_j: 20.0,
                enthalpy_relaxation: 1.0,
            },
            CouplingMethod::RelaxedPicard,
            0.0,
            normal_metric,
            cx,
        );
        assert!(
            matches!(refusal, Err(FiniteThermoplasticError::TransportNotConverged {
            geometry_relative_defect, ..
        }) if geometry_relative_defect > 0.0)
        );
        assert_eq!(thermal.cells(), &accepted_cells);
        assert_eq!(solid.mechanics().element_states(), accepted_history);
        assert_eq!(thermal.time_s(), 0.0);
        assert!(conduction.current_geometry().is_none());

        // A reference-area calculation on the same physics misses the target;
        // the changed heat transfer must affect the actual caloric solve.
        let fixed = solid
            .advance(
                &thermal,
                &conduction,
                first.time_s,
                &first.prescribed(),
                &[(3, first.free_force_n)],
                [0.0; 3],
                settings(),
                controls(),
                thermal_config(),
                config(),
                cx,
            )
            .unwrap();
        assert!((fixed.thermal().temperature_k(CELL).unwrap() - first.temperature_k).abs() > 0.1);

        for expected in &endpoints {
            let (candidate, geometry_defect) = solid
                .advance_with_deformed_conduction(
                    &thermal,
                    &conduction,
                    expected.time_s,
                    &expected.prescribed(),
                    &[(3, expected.free_force_n)],
                    [0.0; 3],
                    settings(),
                    controls(),
                    thermal_config(),
                    config(),
                    CouplingMethod::RelaxedPicard,
                    1.0e-10,
                    normal_metric,
                    cx,
                )
                .unwrap();
            let actual = candidate.solid().mechanics();
            let report = candidate.report();
            let temperature = candidate.thermal().temperature_k(CELL).unwrap();
            assert_eq!(temperature, actual.temperatures_k()[0]);
            close(temperature, expected.temperature_k, 2.0e-7);
            for (actual, wanted) in actual.element_deformation()[0]
                .iter()
                .zip(expected.deformation)
            {
                close(*actual, wanted, 2.0e-9);
            }
            close(
                actual.element_states()[0].alpha,
                expected.increment.mechanical.state.alpha,
                3.0e-9,
            );
            assert!(actual.element_states()[0].alpha > 0.14);
            assert!(std::ptr::eq(
                actual.families()[0].caloric_source(),
                &fixture.caloric
            ));
            assert_eq!(
                candidate.thermal().cells()[&CELL].mass_kg.to_bits(),
                (RHO * V).to_bits()
            );
            assert_eq!(
                candidate.thermal().cells()[&CELL].geometry,
                accepted_cells[&CELL].geometry
            );
            assert!(geometry_defect <= 1.0e-10);
            let used = candidate.observer().current_geometry().unwrap();
            let actual_metrics = normal_metric(actual, cx).unwrap();
            close(
                geometry_defect,
                used.relative_difference(&actual_metrics, || Ok(()))
                    .unwrap(),
                0.0,
            );
            close(report.external_heat_j, expected.boundary_heat_j, 2.0e-8);
            close(
                report.thermal_audit.stored_energy_change_j,
                RHO * V * CP * (expected.temperature_k - expected.old_temperature_k),
                2.0e-8,
            );
            close(report.thermal_residual_j, 0.0, config().energy_tolerance_j);
            close(report.combined_residual_j, 0.0, config().energy_tolerance_j);
            assert_eq!(report.thermal_audit.substeps, 4);
            assert!(report.iterations > 1);
            (thermal, solid, conduction, _) = candidate.into_parts();
        }
        assert!(solid.mechanics().element_states()[0].alpha > 0.3);
        assert!(solid.mechanics().element_states()[0].inverse_plastic_cauchy_green[0] > 1.0);
    });
}
