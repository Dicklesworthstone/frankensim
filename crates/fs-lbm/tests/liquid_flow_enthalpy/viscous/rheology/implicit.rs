//! G1/G3/G4: source material and native flow close on the accepted thermal endpoint.
use super::*;
use fs_lbm::conjugate::enthalpy::advection::FlowEnthalpyInput;
use fs_lbm::conjugate::enthalpy::advection::feedback::MaterialFlowFeedback;
use fs_lbm::conjugate::enthalpy::advection::feedback::implicit::{
    ImplicitFlowControls, ImplicitFlowError,
};
use fs_lbm::conjugate::enthalpy::advection::feedback::rheology::QuasisteadyRheologyControls;
use fs_material::phase::mushy::MushyPhaseLaw;
use fs_material::phase::{EnthalpyPhaseKnot, EquilibriumEnthalpyPhaseCurve};

#[path = "implicit/adaptive.rs"]
mod adaptive;
#[path = "implicit/unsteady.rs"]
mod unsteady;

fn implicit_controls() -> ImplicitFlowControls {
    ImplicitFlowControls {
        max_unknowns: 16,
        max_newton_iterations: 32,
        difference_step_j: 1e-4,
        cell_energy_tolerance_j: 1e-7,
        total_energy_tolerance_j: 4e-7,
        newton: fs_solver::NewtonKrylovConfig {
            absolute_tolerance: 2e-8,
            relative_tolerance: 0.0,
            linear_restart: 8,
            ..Default::default()
        },
    }
}

fn rheology_controls() -> QuasisteadyRheologyControls {
    QuasisteadyRheologyControls {
        max_iterations: 160,
        absolute_rate_tolerance_per_s: 1e-13,
        relative_rate_tolerance: 1e-12,
        relative_transport_tolerance: 1e-12,
        relaxation: 0.5,
    }
}

fn thermal_state<'a>(
    domain: &VoxelDomain,
    caloric: &'a dyn ThermalLaw,
    h: &[f64],
    links: &[ThermalLink],
) -> DepositionSimulation<'a> {
    let mut events = (0..domain.cell_count())
        .map(|c| {
            let [x, y, z] = domain.coords(c);
            let mut start = domain.center(x, y, z);
            let mut end = start;
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
                    specific_h_j_kg: h[c],
                    geometry: CellGeometry {
                        start_m: start,
                        end_m: end,
                        volume_m3: 1.0,
                        surface_area_m2: 6.0,
                    },
                    boundary: Boundary {
                        ambient_k: 300.0,
                        htc_w_m2_k: 0.0,
                        reservoir: None,
                        source_w: 0.0,
                    },
                }),
            }
        })
        .collect::<Vec<_>>();
    events.extend(links.iter().cloned().map(|link| TimedEvent {
        time_s: 0.0,
        event: Event::Contact(link),
    }));
    DepositionSimulation::new(
        events,
        &[(MAT, caloric)],
        DepositionBudget {
            max_cells: domain.cell_count(),
            max_links: 64,
            max_events: domain.cell_count() + links.len(),
            max_substeps: 1000,
        },
    )
    .unwrap()
    .advance_to(0.0, thermal_step(), || Ok(()))
    .unwrap()
    .0
}

fn pressure_channel() -> SimpleConfig {
    let mut faces = [FvBoundary::Symmetry; 6];
    faces[0] = FvBoundary::Pressure {
        pressure_pa: 100030.0,
    };
    faces[1] = FvBoundary::Pressure {
        pressure_pa: 100000.0,
    };
    faces[2] = FvBoundary::wall();
    faces[3] = FvBoundary::wall();
    let mut config = SimpleConfig::new(faces);
    config.tolerance = 1e-13;
    config.pressure_tolerance = 1e-12;
    config.momentum_tolerance = 1e-10;
    config
}

#[test]
fn thermal_and_rate_material_feedback_matches_one_independent_coupled_root() {
    let domain = VoxelDomain::new(2, 2, 1, 1.0).unwrap();
    let config = closed();
    let gate = CancelGate::new();
    let ids: Vec<_> = (0..4).map(|c| Some(CellId(c))).collect();
    let dt = 2.0;
    let mut speeds = Vec::new();
    for (temperature_slope, rate_slope) in [(0.0, 0.0), (40.0, 0.0), (40.0, 100.0)] {
        let card = rheological_card(temperature_slope, rate_slope, 2.0);
        let liquid = law(&card);
        let caloric = caloric(&card);
        let old = initial(&domain, 1.0, &caloric);
        let source = LiquidFlowRheology::try_new(input(&ids, &liquid, &caloric), "rate").unwrap();
        let flow = source
            .quasisteady_feedback(
                &old,
                &domain,
                &config,
                (),
                &gate,
                8192,
                vec![0.1; 4],
                rheology_controls(),
            )
            .unwrap()
            .with_body_force_heating(force(), controls())
            .unwrap();
        // Exact discrete circulation: mu*q=50/3; each cell receives 100*q W.
        // BE gives h=100*dt*q and T-300=dt*q. Together with gamma=sqrt(6)*q,
        // [(s_T*dt+s_gamma*sqrt(6))*q+100]*q=50/3.
        let a = temperature_slope * dt + rate_slope * 6.0_f64.sqrt();
        let q = (100.0 / 3.0) / (100.0 + (10000.0 + (200.0 / 3.0) * a).sqrt());
        let expected_h = 100.0 * dt * q;
        let (next, report) = flow
            .step_implicit(&old, dt, implicit_controls(), &[10.0; 4], || Ok(()))
            .unwrap();
        near(report.endpoint.flow.field.flux_x(1, 1, 0), q);
        near(report.endpoint.time_s, dt);
        for sample in &report.endpoint.samples {
            near(sample.specific_enthalpy_j_kg, expected_h);
            near(next.cells()[&sample.cell].specific_h_j_kg, expected_h);
            near(sample.temperature_k, 300.0 + dt * q);
            near(sample.material.dynamic_viscosity_pa_s(), 100.0 + a * q);
            near(
                sample.material.thermal_conductivity_w_m_k(),
                0.01 + 0.05 * 6.0_f64.sqrt() * q,
            );
            assert_eq!(
                sample.material.material().card_identity(),
                card.content_hash()
            );
        }
        near(report.mechanical.supplied_body_work_j, dt * 400.0 * q);
        near(report.deposited_heat_j, dt * 400.0 * q);
        near(
            report.thermal.stored_energy_change_j,
            report.deposited_heat_j,
        );
        near(report.net_stream_enthalpy_j, 0.0);
        near(report.other_thermal_input_j, 0.0);
        assert!(report.combined_residual_j.abs() <= implicit_controls().total_energy_tolerance_j);
        assert!(
            report.nonlinear.converged
                && report.momentum_solves > report.rheology.unwrap().iterations
        );
        assert!(flow.applied().is_none());
        assert_eq!(flow.momentum_solves(), 0);
        assert_eq!(old.energy_j(), 0.0);
        speeds.push(q);
    }
    assert!(speeds[0] > speeds[1] && speeds[1] > speeds[2]);
}

#[test]
fn pressure_work_advection_and_source_conduction_share_the_implicit_endpoint() {
    let domain = VoxelDomain::new(3, 2, 1, 1.0).unwrap();
    let config = pressure_channel();
    let gate = CancelGate::new();
    let ids: Vec<_> = (0..6).map(|c| Some(CellId(c))).collect();
    for (mu, k) in [(100.0, 0.01), (200.0, 2.0)] {
        let original = source(mu, false);
        let mut claims = ClaimSet::new();
        for (_, original) in original.claims().claims_ordered() {
            let mut claim = original.clone();
            if claim.key.name() == THERMAL_CONDUCTIVITY_PROPERTY {
                claim.value = PropertyValue::Scalar {
                    value: k,
                    dims: THERMAL_CONDUCTIVITY_DIMS,
                };
            }
            claims.insert_claim(claim).unwrap();
        }
        let card = card(claims);
        let liquid = liquid(&card);
        let caloric = liquid_caloric(&card);
        let link = ThermalLink {
            id: 100,
            a: CellId(1),
            b: CellId(4),
            area_m2: 1.0,
            conductance_w_k: 999.0,
        };
        let old_h = [100.0, 100.0, 100.0, 0.0, 0.0, 0.0];
        let old = thermal_state(&domain, caloric.curve(), &old_h, &[link]);
        let flow = input(&ids, &liquid, &caloric)
            .feedback(&old, &domain, &config, (), &gate, 512)
            .unwrap()
            .with_pressure_work_heating(controls())
            .unwrap()
            .with_fluid_conduction(&old, &[100])
            .unwrap();
        // Independent donor-chain elimination with one cross-row Fourier link.
        // Each mouth carries u=5/mu, each cell receives 10u W, and G=k.
        let dt = 4.0;
        let u = 5.0 / mu;
        let lambda = dt * u;
        let eta = dt * k / 100.0;
        let heat = dt * 10.0 * u;
        let hot0 = (100.0 + heat) / (1.0 + lambda);
        let cold0 = heat / (1.0 + lambda);
        let sum1 = (100.0 + 2.0 * heat + lambda * (hot0 + cold0)) / (1.0 + lambda);
        let difference1 = (100.0 + lambda * (hot0 - cold0)) / (1.0 + lambda + 2.0 * eta);
        let hot1 = 0.5 * (sum1 + difference1);
        let cold1 = 0.5 * (sum1 - difference1);
        let hot2 = (100.0 + heat + lambda * hot1) / (1.0 + lambda);
        let cold2 = (heat + lambda * cold1) / (1.0 + lambda);
        let expected = [hot0, hot1, hot2, cold0, cold1, cold2];
        let (next, report) = flow
            .step_implicit(&old, dt, implicit_controls(), &[1.0; 6], || Ok(()))
            .unwrap();
        for (i, expected) in expected.into_iter().enumerate() {
            near(next.cells()[&CellId(i as u64)].specific_h_j_kg, expected);
            near(report.endpoint.samples[i].specific_enthalpy_j_kg, expected);
            near(
                report.endpoint.samples[i]
                    .material
                    .thermal_conductivity_w_m_k(),
                k,
            );
        }
        near(report.mechanical.supplied_pressure_work_j, dt * 60.0 * u);
        near(report.deposited_heat_j, dt * 60.0 * u);
        near(report.net_stream_enthalpy_j, -dt * u * (hot2 + cold2));
        near(
            next.energy_j() - old.energy_j(),
            report.net_stream_enthalpy_j + report.deposited_heat_j,
        );
        assert!(report.combined_residual_j.abs() <= implicit_controls().total_energy_tolerance_j);
        assert_eq!(next.links()[&100].conductance_w_k, 999.0);
        for cell in next.cells().values() {
            assert_eq!(cell.boundary.source_w, 0.0);
        }
    }
}

#[test]
fn latent_storage_and_endpoint_mobility_close_without_a_temperature_change() {
    let domain = VoxelDomain::new(2, 2, 1, 1.0).unwrap();
    let config = closed();
    let gate = CancelGate::new();
    let ids: Vec<_> = (0..4).map(|c| Some(CellId(c))).collect();
    let card = source(100.0, false);
    let liquid = liquid(&card);
    let curve = EquilibriumEnthalpyPhaseCurve::try_new(
        card.content_hash(),
        [(0.0, 0.0), (100.0, 1.0)]
            .map(|(h, f)| EnthalpyPhaseKnot {
                specific_enthalpy_j_kg: h,
                temperature_k: 500.0,
                liquid_mass_fraction: f,
                bulk_density_kg_m3: 1.0,
            })
            .to_vec(),
    )
    .unwrap();
    let old = thermal_state(&domain, &curve, &[50.0; 4], &[]);
    let mut speeds = Vec::new();
    for permeability in [0.5, 2.0] {
        let mobility = MushyPhaseLaw::try_new(&curve, 1.0, 1.0, permeability, 0.2).unwrap();
        let flow = MaterialFlowFeedback::try_new(
            &old,
            &domain,
            FlowEnthalpyInput {
                cell_ids: &ids,
                material_id: MAT,
                caloric: &curve,
                inlet_specific_h_j_kg: [Some(50.0); 6],
                max_cells: 4,
                max_connections: 64,
            },
            &config,
            (),
            &gate,
            512,
            |_, _: &DepositedCell, t| {
                liquid.at_temperature(Temperature::new(t)).map_err(|e| {
                    fs_material::MaterialError::State {
                        what: e.to_string(),
                    }
                })
            },
        )
        .unwrap()
        .with_body_force_heating(force(), controls())
        .unwrap()
        .with_mushy_mobility(&mobility)
        .unwrap();
        let dt = 0.1;
        let drag = |h: f64| {
            let phi = h / 100.0;
            100.0 * (1.0 - phi).powi(2) / (permeability * phi.powi(3))
        };
        // Independent scalar equation: h-50=dt*100q and
        // q=100/(600+mu/K(h)). Its unique small-step root lies in [50,52].
        let (mut lo, mut hi) = (50.0, 52.0);
        for _ in 0..64 {
            let h = 0.5 * (lo + hi);
            if h - 50.0 > dt * 10000.0 / (600.0 + drag(h)) {
                hi = h;
            } else {
                lo = h;
            }
        }
        let h = 0.5 * (lo + hi);
        let q = 100.0 / (600.0 + drag(h));
        let (next, report) = flow
            .step_implicit(&old, dt, implicit_controls(), &[1.0; 4], || Ok(()))
            .unwrap();
        near(report.endpoint.flow.field.flux_x(1, 1, 0), q);
        for sample in &report.endpoint.samples {
            near(sample.temperature_k, 500.0);
            near(sample.specific_enthalpy_j_kg, h);
            near(next.cells()[&sample.cell].specific_h_j_kg, h);
            near(
                report.endpoint.mushy_drag_pa_s_m2.as_ref().unwrap()[sample.voxel],
                drag(h),
            );
        }
        near(report.deposited_heat_j, 4.0 * (h - 50.0));
        near(report.mechanical.supplied_body_work_j, dt * 400.0 * q);
        near(report.mechanical.darcy_heat_j, dt * 4.0 * drag(h) * q * q);
        near(report.mechanical.molecular_heat_j, dt * 2400.0 * q * q);
        near(next.mass_kg(), old.mass_kg());
        assert!(report.endpoint.mushy_drag_pa_s_m2.as_ref().unwrap()[0] < drag(50.0));
        speeds.push(q);
    }
    assert!(speeds[1] > speeds[0]);
}

#[test]
fn constant_material_keeps_the_original_flow_and_physical_heat_accounts() {
    let domain = VoxelDomain::new(2, 2, 1, 1.0).unwrap();
    let config = closed();
    let gate = CancelGate::new();
    let ids: Vec<_> = (0..4).map(|c| Some(CellId(c))).collect();
    let card = source(100.0, false);
    let liquid = liquid(&card);
    let caloric = liquid_caloric(&card);
    let old = initial(&domain, 1.0, &caloric);
    let flow = input(&ids, &liquid, &caloric)
        .feedback(&old, &domain, &config, (), &gate, 512)
        .unwrap()
        .with_body_force_heating(force(), controls())
        .unwrap();
    let dt = 0.125;
    let (explicit, history, _) = old
        .advance_observed(&flow, dt, thermal_step(), || Ok(()))
        .unwrap();
    for guess in [0.0, 1.0] {
        let (implicit, report) = flow
            .step_implicit(&old, dt, implicit_controls(), &[guess; 4], || Ok(()))
            .unwrap();
        assert_eq!(report.endpoint.flow, history.applied().unwrap().flow);
        assert_eq!(report.mechanical, history.viscous_heating_audit().unwrap());
        for (id, cell) in implicit.cells() {
            near(cell.specific_h_j_kg, explicit.cells()[id].specific_h_j_kg);
        }
        near(report.deposited_heat_j, dt * 200.0 / 3.0);
        near(
            report.thermal.stored_energy_change_j,
            report.deposited_heat_j,
        );
        assert!(report.momentum_solves >= 2);
    }
}

#[test]
fn shared_work_source_and_late_cancellation_refuse_without_advancing_any_history() {
    let domain = VoxelDomain::new(2, 2, 1, 1.0).unwrap();
    let config = closed();
    let gate = CancelGate::new();
    let ids: Vec<_> = (0..4).map(|c| Some(CellId(c))).collect();
    let card = source(100.0, true);
    let liquid = liquid(&card);
    let caloric = liquid_caloric(&card);
    let old = initial(&domain, 1.0, &caloric);
    let make = |budget| {
        input(&ids, &liquid, &caloric)
            .feedback(&old, &domain, &config, (), &gate, budget)
            .unwrap()
            .with_body_force_heating(force(), controls())
            .unwrap()
    };
    let limited = make(1);
    let dt = 0.5;
    assert!(matches!(
        limited.step_implicit(&old, dt, implicit_controls(), &[1.0; 4], || Ok(())),
        Err(ImplicitFlowError::Flow(MaterialFlowError::Thermal(
            DepositionError::Budget(_)
        )))
    ));
    let flow = make(512);
    let one = ImplicitFlowControls {
        max_newton_iterations: 1,
        ..implicit_controls()
    };
    assert!(matches!(
        flow.step_implicit(&old, dt, one, &[1.0; 4], || Ok(())),
        Err(ImplicitFlowError::NotConverged(_))
    ));
    let mut polls = 0;
    let (complete, first) = flow
        .step_implicit(&old, dt, implicit_controls(), &[1.0; 4], || {
            polls += 1;
            Ok(())
        })
        .unwrap();
    assert!(polls > 2 && first.momentum_solves > 1);
    let mut repeat = 0;
    assert!(matches!(
        flow.step_implicit(&old, dt, implicit_controls(), &[1.0; 4], || {
            repeat += 1;
            if repeat == polls - 1 {
                Err(DepositionError::Cancelled)
            } else {
                Ok(())
            }
        }),
        Err(ImplicitFlowError::Flow(MaterialFlowError::Thermal(
            DepositionError::Cancelled
        )))
    ));
    let (retry, repeated) = flow
        .step_implicit(&old, dt, implicit_controls(), &[1.0; 4], || Ok(()))
        .unwrap();
    assert_eq!(retry.cells(), complete.cells());
    assert_eq!(repeated.endpoint, first.endpoint);
    assert_eq!(repeated.momentum_solves, first.momentum_solves);
    assert_eq!(old.time_s(), 0.0);
    assert_eq!(old.energy_j(), 0.0);
    assert_eq!(flow.momentum_solves(), 0);
    assert!(flow.applied().is_none());
    assert_eq!(
        flow.viscous_heating_audit(),
        Some(ViscousHeatingAudit::default())
    );

    // The initial iterate is source-supported; its actual heated endpoint is
    // not. An unsupported Newton trial is a physical refusal, not convergence.
    let mut claims = ClaimSet::new();
    for (_, original) in card.claims().claims_ordered() {
        let mut claim = original.clone();
        if claim.key.name() == DYNAMIC_VISCOSITY_PROPERTY {
            claim.validity = claim.validity.with_quantity("T", absolute(), 300.0, 300.05);
        }
        claims.insert_claim(claim).unwrap();
    }
    let bounded_card = crate::card(claims);
    let bounded_liquid = crate::liquid(&bounded_card);
    let bounded_caloric = liquid_caloric(&bounded_card);
    let bounded_old = initial(&domain, 1.0, &bounded_caloric);
    let bounded = input(&ids, &bounded_liquid, &bounded_caloric)
        .feedback(&bounded_old, &domain, &config, (), &gate, 512)
        .unwrap()
        .with_body_force_heating(force(), controls())
        .unwrap();
    assert!(matches!(
        bounded.step_implicit(&bounded_old, dt, implicit_controls(), &[1.0; 4], || Ok(())),
        Err(ImplicitFlowError::Flow(MaterialFlowError::Material { .. }))
    ));
    assert_eq!(bounded_old.time_s(), 0.0);
    assert_eq!(bounded_old.energy_j(), 0.0);
    assert!(bounded.applied().is_none());
}
