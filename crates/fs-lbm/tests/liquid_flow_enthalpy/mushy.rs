//! Actual h-dependent mixture mobility on the original native heat timeline.
use super::*;
use fs_lbm::conjugate::enthalpy::advection::FlowEnthalpyInput;
use fs_lbm::conjugate::enthalpy::advection::feedback::MaterialFlowFeedback;
use fs_material::phase::mushy::MushyPhaseLaw;
use fs_material::phase::{EnthalpyPhaseKnot, EquilibriumEnthalpyPhaseCurve};

#[path = "mushy/heating.rs"]
mod heating;

#[path = "mushy/conduction.rs"]
mod conduction;

#[path = "mushy/electrical.rs"]
mod electrical;

#[path = "mushy/continuity.rs"]
mod continuity;

#[path = "mushy/occupancy.rs"]
mod occupancy;

#[path = "mushy/release.rs"]
mod release;

#[path = "mushy/unsteady.rs"]
mod unsteady;

fn phase_curve(card: &MaterialCard, rho: [f64; 2]) -> EquilibriumEnthalpyPhaseCurve {
    EquilibriumEnthalpyPhaseCurve::try_new(
        card.content_hash(),
        vec![
            EnthalpyPhaseKnot {
                specific_enthalpy_j_kg: 0.0,
                temperature_k: 500.0,
                liquid_mass_fraction: 0.0,
                bulk_density_kg_m3: rho[0],
            },
            EnthalpyPhaseKnot {
                specific_enthalpy_j_kg: 100.0,
                temperature_k: 500.0,
                liquid_mass_fraction: 1.0,
                bulk_density_kg_m3: rho[1],
            },
        ],
    )
    .unwrap()
}

fn initial<'a>(
    domain: &VoxelDomain,
    caloric: &'a EquilibriumEnthalpyPhaseCurve,
    source_w: f64,
) -> DepositionSimulation<'a> {
    let events = (0..domain.cell_count())
        .map(|c| {
            let [x, y, z] = domain.coords(c);
            let center = domain.center(x, y, z);
            let mut start = center;
            let mut end = center;
            start[0] -= 0.5;
            end[0] += 0.5;
            TimedEvent {
                time_s: 0.0,
                event: Event::Birth(DepositedCell {
                    id: CellId(c as u64),
                    segment_id: c as u64,
                    body_id: 1,
                    material_id: MAT,
                    mass_kg: 1.0,
                    specific_h_j_kg: 50.0,
                    geometry: CellGeometry {
                        start_m: start,
                        end_m: end,
                        volume_m3: 1.0,
                        surface_area_m2: 6.0,
                    },
                    boundary: Boundary {
                        ambient_k: 500.0,
                        htc_w_m2_k: 0.0,
                        reservoir: None,
                        source_w: if c == 1 { source_w } else { 0.0 },
                    },
                }),
            }
        })
        .collect();
    DepositionSimulation::new(
        events,
        &[(MAT, caloric as &dyn ThermalLaw)],
        DepositionBudget {
            max_cells: 4,
            max_links: 64,
            max_events: 4,
            max_substeps: 1000,
        },
    )
    .unwrap()
    .advance_to(0.0, config(), || Ok(()))
    .unwrap()
    .0
}

fn input<'a>(
    ids: &'a [Option<CellId>],
    curve: &'a EquilibriumEnthalpyPhaseCurve,
) -> FlowEnthalpyInput<'a> {
    FlowEnthalpyInput {
        cell_ids: ids,
        material_id: MAT,
        caloric: curve,
        inlet_specific_h_j_kg: [Some(50.0); 6],
        max_cells: 4,
        max_connections: 64,
    }
}

fn momentum() -> SimpleConfig {
    let mut faces = [FvBoundary::Symmetry; 6];
    faces[0] = FvBoundary::Outlet;
    faces[1] = FvBoundary::Outlet;
    let mut config = SimpleConfig::new(faces);
    config.tolerance = 1e-13;
    config.pressure_tolerance = 1e-12;
    config.internal_fans.push(InternalFan {
        patch: FacePatch {
            axis: 0,
            index: 2,
            lo: [0, 0, 0],
            hi: [4, 1, 1],
        },
        blows_positive: true,
        curve: FanCurve::new(&[(0.0, 4.0), (100.0, 4.0)]).unwrap(),
    });
    config
}

#[test]
fn latent_plateau_changes_native_momentum_and_energy_by_material_data_only() {
    let card = card(claims(1.0, 1.0, false));
    let liquid = liquid(&card);
    let curve = phase_curve(&card, [1.0; 2]);
    let domain = VoxelDomain::new(4, 1, 1, 1.0).unwrap();
    let config_flow = momentum();
    let gate = CancelGate::new();
    let ids: Vec<_> = (0..4).map(|c| Some(CellId(c))).collect();
    let old = initial(&domain, &curve, 0.0);
    let mut outputs = Vec::new();
    for scale in [1.0, 2.0] {
        let mobility = MushyPhaseLaw::try_new(&curve, 1.0, 1.0, scale, 0.1).unwrap();
        let resolve = |_: usize, _: &DepositedCell, t: f64| {
            liquid.at_temperature(Temperature::new(t)).map_err(|error| {
                fs_material::MaterialError::State {
                    what: error.to_string(),
                }
            })
        };
        let observer = MaterialFlowFeedback::try_new(
            &old,
            &domain,
            input(&ids, &curve),
            &config_flow,
            Heater::default(),
            &gate,
            100,
            resolve,
        )
        .unwrap()
        .with_mushy_mobility(&mobility)
        .unwrap();
        let (next, accepted, audit) = old
            .advance_observed(
                &observer,
                0.1,
                StepConfig {
                    max_dt_s: 0.0125,
                    energy_tolerance_j: 1e-9,
                },
                || Ok(()),
            )
            .unwrap();
        let mut h = [50.0_f64; 4];
        let mut boundary_j = 0.0;
        let mut last_drag = [0.0; 4];
        let mut last_q = 0.0;
        for &dt in &accepted.inner().durations {
            // Independent pressure/series-Darcy and total-h donor-cell balance.
            // Temperature stays exactly 500K, so no T-only update can pass this.
            last_drag = h.map(|h| {
                let phi = h / 100.0;
                (1.0 - phi).powi(2) / (scale * phi.powi(3))
            });
            last_q = 4.0 / last_drag.iter().sum::<f64>();
            boundary_j += dt * last_q * (50.0 - h[3]);
            let previous = h;
            for c in 0..4 {
                let upstream = if c == 0 { 50.0 } else { previous[c - 1] };
                h[c] += dt * (if c == 1 { 100.0 } else { 0.0 } + last_q * (upstream - previous[c]));
            }
        }
        assert!(accepted.momentum_solves() > 2);
        assert_eq!(accepted.momentum_solves(), accepted.inner().durations.len());
        let applied = accepted.applied().unwrap();
        near(applied.flow.field.flux_x(2, 0, 0), last_q);
        assert!(last_q > 0.5 * scale);
        for (actual, expected) in applied
            .mushy_drag_pa_s_m2
            .as_ref()
            .unwrap()
            .iter()
            .zip(last_drag)
        {
            near(*actual, expected);
        }
        for sample in accepted.endpoint() {
            near(next.cells()[&sample.cell].specific_h_j_kg, h[sample.voxel]);
            assert_eq!(sample.temperature_k, 500.0);
            assert_eq!(
                sample.material.material().card_identity(),
                card.content_hash()
            );
            assert_eq!(
                sample
                    .material
                    .material()
                    .property(DYNAMIC_VISCOSITY_PROPERTY)
                    .unwrap()
                    .answer()
                    .receipt
                    .selected,
                liquid
                    .reference_state()
                    .material()
                    .property(DYNAMIC_VISCOSITY_PROPERTY)
                    .unwrap()
                    .answer()
                    .receipt
                    .selected
            );
        }
        near(accepted.audit().external_enthalpy_j, boundary_j);
        near(next.energy_j() - old.energy_j(), 10.0 + boundary_j);
        near(audit.external_heat_j, 10.0 + boundary_j);
        near(next.mass_kg(), old.mass_kg());
        assert_eq!(observer.momentum_solves(), 0);
        outputs.push((last_q, h[1]));
    }
    assert!(outputs[1].0 > outputs[0].0);
    assert!(outputs[1].1 < outputs[0].1);
}

#[test]
fn phase_volume_change_wrong_chart_and_solid_locking_refuse_without_state_publication() {
    let card = card(claims(1.0, 1.0, false));
    let liquid = liquid(&card);
    let curve = phase_curve(&card, [1.0; 2]);
    let clone = curve.clone();
    let unequal = phase_curve(&card, [1.0, 2.0]);
    let wrong_chart = MushyPhaseLaw::try_new(&clone, 1.0, 1.0, 1.0, 0.2).unwrap();
    let shrinkage = MushyPhaseLaw::try_new(&unequal, 1.0, 2.0, 1.0, 0.2).unwrap();
    let mobility = MushyPhaseLaw::try_new(&curve, 1.0, 1.0, 1.0, 0.2).unwrap();
    let domain = VoxelDomain::new(4, 1, 1, 1.0).unwrap();
    let config_flow = momentum();
    let gate = CancelGate::new();
    let ids: Vec<_> = (0..4).map(|c| Some(CellId(c))).collect();
    let old = initial(&domain, &curve, -400.0);
    let observer = MaterialFlowFeedback::try_new(
        &old,
        &domain,
        input(&ids, &curve),
        &config_flow,
        (),
        &gate,
        100,
        |_, _: &DepositedCell, t| {
            liquid.at_temperature(Temperature::new(t)).map_err(|error| {
                fs_material::MaterialError::State {
                    what: error.to_string(),
                }
            })
        },
    )
    .unwrap();
    assert!(observer.clone().with_mushy_mobility(&wrong_chart).is_err());
    assert!(observer.clone().with_mushy_mobility(&shrinkage).is_err());
    let observer = observer.with_mushy_mobility(&mobility).unwrap();
    assert!(
        old.advance_observed(
            &observer,
            0.1,
            StepConfig {
                max_dt_s: 0.1,
                energy_tolerance_j: 1e-9
            },
            || Ok(())
        )
        .is_err()
    );
    assert_eq!(old.time_s(), 0.0);
    assert_eq!(old.cells()[&CellId(1)].specific_h_j_kg, 50.0);
    assert_eq!(observer.momentum_solves(), 0);
    assert!(observer.applied().is_none());
}
