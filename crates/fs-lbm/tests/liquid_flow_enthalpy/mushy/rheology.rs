//! The original phase timeline consumes measured-form T/rate source tables.
//! Manufactured coefficients isolate native strain, latent heat and dissipation.
use super::*;
use fs_lbm::conjugate::ViscousHeatControls;
use fs_lbm::conjugate::enthalpy::advection::feedback::{
    PhaseFlowTransport, rheology::FlowRheologySource,
};
use fs_material::MaterialError;
use fs_material::phase::transport::{
    PhaseTransportSourceQuery, ResolvedPhaseTransportState,
    resolve_phase_transport_state_with_shear_rate,
};
use fs_material::thermal::source::ThermalMaterialPoint;
use fs_qty::Frequency;

pub(super) fn source(rate_slope: f64, max_rate: f64) -> MaterialCard {
    let original = super::continuity::source(100.0);
    let mut claims = ClaimSet::new();
    for (_, claim) in original.claims().claims_ordered() {
        let mut claim = claim.clone();
        let viscosity = claim.key.name() == DYNAMIC_VISCOSITY_PROPERTY;
        if viscosity || claim.key.name() == THERMAL_CONDUCTIVITY_PROPERTY {
            claim.value = PropertyValue::Surface2D {
                axes: ["T".into(), "rate".into()],
                axis_dims: [Temperature::DIMS, Frequency::DIMS],
                samples: [
                    [450.0, 0.0],
                    [600.0, 0.0],
                    [450.0, max_rate],
                    [600.0, max_rate],
                ]
                .into_iter()
                .map(|p| {
                    (
                        p,
                        if viscosity {
                            100.0 - 0.05 * (p[0] - 550.0) + rate_slope * p[1]
                        } else {
                            0.2 + 0.001 * (p[0] - 500.0) + 0.05 * p[1]
                        },
                    )
                })
                .collect(),
                triangles: vec![[0, 1, 2], [1, 3, 2]],
                dims: claim.key.dims(),
            };
            claim.interpolation = InterpolationPolicy::LinearInside;
        }
        claim.validity = claim.validity.with_quantity(
            "rate",
            QuantitySpec::dimensional(Frequency::DIMS),
            0.0,
            max_rate,
        );
        claims.insert_claim(claim).unwrap();
    }
    MaterialCard::assemble(original.id().clone(), claims, vec![]).unwrap()
}

pub(super) fn resolve_at(
    card: &MaterialCard,
    law: &MushyPhaseLaw<'_>,
    h: f64,
    rate: f64,
) -> Result<ResolvedPhaseTransportState, MaterialError> {
    let context = point(500.0, 100000.0)
        .with("liquid", 0.0)
        .unwrap()
        .with_quantity("rate", QuantitySpec::dimensional(Frequency::DIMS), 0.0)
        .unwrap();
    resolve_phase_transport_state_with_shear_rate(
        law,
        h,
        PhaseTransportSourceQuery {
            source: ThermalMaterialPoint {
                material: card,
                context: &context,
                temperature_axis: "T",
                liquid_fraction_axis: Some("liquid"),
            },
            viscosity_basis: LiquidViscosityBasis::Dynamic,
            extra_requirements: &[],
            selection: MaterialPropertySelection::SingleClaimOnly,
        },
        "rate",
        rate,
    )
}

fn cold<'a>(
    domain: &VoxelDomain,
    curve: &'a EquilibriumEnthalpyPhaseCurve,
) -> DepositionSimulation<'a> {
    let seed = super::initial(domain, curve, 0.0);
    let mut events: Vec<_> = seed
        .cells()
        .values()
        .cloned()
        .map(|mut cell| {
            cell.specific_h_j_kg = -20.0;
            cell.boundary.source_w = 2000.0;
            TimedEvent {
                time_s: 0.0,
                event: Event::Birth(cell),
            }
        })
        .collect();
    for cell in seed.cells().values() {
        events.push(TimedEvent {
            time_s: 0.11,
            event: Event::Boundary(cell.id, cell.boundary),
        });
    }
    DepositionSimulation::new(
        events,
        &[(MAT, curve as &dyn ThermalLaw)],
        DepositionBudget {
            max_cells: 4,
            max_links: 4,
            max_events: 8,
            max_substeps: 64,
        },
    )
    .unwrap()
    .advance_to(0.0, step(), || Ok(()))
    .unwrap()
    .0
}
fn closed() -> SimpleConfig {
    let mut faces = [FvBoundary::wall(); 6];
    faces[4] = FvBoundary::Symmetry;
    faces[5] = FvBoundary::Symmetry;
    let mut config = SimpleConfig::new(faces);
    config.tolerance = 1e-13;
    config.pressure_tolerance = 1e-12;
    config.momentum_tolerance = 1e-3;
    config
}
fn controls() -> ViscousHeatControls {
    ViscousHeatControls {
        maximum_reynolds: 0.1,
        power_tolerance_w: 1e-8,
    }
}
fn step() -> StepConfig {
    StepConfig {
        max_dt_s: 0.2,
        energy_tolerance_j: 1e-8,
    }
}
fn forces() -> Vec<[f64; 3]> {
    vec![
        [-100.0, 100.0, 0.0],
        [-100.0, -100.0, 0.0],
        [100.0, 100.0, 0.0],
        [100.0, -100.0, 0.0],
    ]
}
#[allow(clippy::too_many_arguments)]
fn bound<'a>(
    thermal: &DepositionSimulation<'_>,
    domain: &'a VoxelDomain,
    ids: &'a [Option<CellId>],
    config: &'a SimpleConfig,
    gate: &'a CancelGate,
    card: &'a MaterialCard,
    law: &'a MushyPhaseLaw<'a>,
    source: &'a dyn FlowRheologySource<ResolvedPhaseTransportState>,
) -> MaterialFlowFeedback<
    'a,
    ResolvedPhaseTransportState,
    impl Clone
    + FnMut(usize, &DepositedCell, f64) -> Result<ResolvedPhaseTransportState, MaterialError>
    + 'a,
    (),
    PhaseFlowTransport,
> {
    MaterialFlowFeedback::try_new_phase(
        thermal,
        domain,
        FlowEnthalpyInput {
            cell_ids: ids,
            material_id: MAT,
            caloric: law.curve(),
            inlet_specific_h_j_kg: [None; 6],
            max_cells: 4,
            max_connections: 32,
        },
        config,
        (),
        gate,
        16,
        law,
        move |_, cell: &DepositedCell, _| resolve_at(card, law, cell.specific_h_j_kg, 0.0),
    )
    .unwrap()
    .with_unsteady_viscous_heating(controls())
    .unwrap()
    .with_native_rheology(thermal, source)
    .unwrap()
}

#[test]
fn real_phase_shear_tables_drive_native_motion_and_heat_after_latent_melting() {
    let domain = VoxelDomain::new(2, 2, 1, 1.0).unwrap();
    let ids: Vec<_> = (0..4).map(|c| Some(CellId(c))).collect();
    let config = closed();
    let gate = CancelGate::new_clock_free();
    let mut outcomes = Vec::new();
    for slope in [0.0, 20.0, -20.0] {
        let card = source(slope, 2.0);
        let curve = super::continuity::chart(&card);
        let law = MushyPhaseLaw::try_new_with_immobile_solid(&curve, 1.0, 1.0, 1.0).unwrap();
        let original = cold(&domain, &curve);
        let source = |_: usize, cell: &DepositedCell, _: f64, rate: f64| {
            resolve_at(&card, &law, cell.specific_h_j_kg, rate)
        };
        let history = bound(
            &original, &domain, &ids, &config, &gate, &card, &law, &source,
        );
        assert_eq!(history.endpoint_shear_rate_per_s().unwrap(), &[0.0; 4]);
        assert!(history.endpoint().iter().all(|sample| {
            sample.material.mobile_transport().is_none()
                && sample
                    .material
                    .material()
                    .property(DYNAMIC_VISCOSITY_PROPERTY)
                    .is_none()
        }));
        let (molten, history, heat) = original
            .advance_observed(&history, 0.11, step(), || Ok(()))
            .unwrap();
        near(heat.external_heat_j, 4.0 * 220.0);
        assert_eq!(history.unsteady_time_s(), Some(0.11));
        for sample in history.endpoint() {
            near(sample.specific_enthalpy_j_kg, 200.0);
            near(sample.temperature_k, 550.0);
        }
        let history = history
            .with_unsteady_body_force_heating(forces(), controls())
            .unwrap();
        let (warm, first, heat) = molten
            .advance_observed(&history, 0.21, step(), || Ok(()))
            .unwrap();
        // Four free MAC face masses and symmetric-stress/no-slip dissipation
        // reduce independently to (4/dt + 24 mu) q = 400, gamma² = 6 q².
        let dt = 0.1;
        let q1 = 400.0 / (4.0 / dt + 24.0 * 100.0);
        let rate = 6.0_f64.sqrt() * q1;
        let energy1 = dt * 24.0 * 100.0 * q1 * q1;
        let t1 = 550.0 + energy1 / (4.0 * 2.0);
        let mu1 = 100.0 - 0.05 * (t1 - 550.0) + slope * rate;
        near(first.unsteady_flow().unwrap().field.flux_x(1, 1, 0), q1);
        near(heat.external_heat_j, energy1);
        near(warm.energy_j() - molten.energy_j(), energy1);
        for sample in first.endpoint() {
            near(
                first.endpoint_shear_rate_per_s().unwrap()[sample.voxel],
                rate,
            );
            near(sample.temperature_k, t1);
            near(
                sample
                    .transport
                    .mobile_transport()
                    .unwrap()
                    .dynamic_viscosity_pa_s(),
                mu1,
            );
            near(
                sample.transport.thermal_conductivity_w_m_k(),
                0.2 + 0.001 * (t1 - 500.0) + 0.05 * rate,
            );
            assert_eq!(
                sample.material.material().card_identity(),
                card.content_hash()
            );
            near(sample.material.material().query().axes()["rate"], rate);
        }
        let (next, second, heat) = warm
            .advance_observed(&first, 0.31, step(), || Ok(()))
            .unwrap();
        let q2 = (4.0 * q1 / dt + 400.0) / (4.0 / dt + 24.0 * mu1);
        let energy2 = dt * 24.0 * mu1 * q2 * q2;
        let mechanical = second.unsteady_mechanical_step().unwrap();
        near(second.unsteady_flow().unwrap().field.flux_x(1, 1, 0), q2);
        near(heat.external_heat_j, energy2);
        near(next.energy_j() - warm.energy_j(), energy2);
        near(mechanical.power.molecular_heat_w, energy2 / dt);
        near(mechanical.power.darcy_heat_w, 0.0);
        near(mechanical.final_kinetic_energy_j, 2.0 * q2 * q2);
        near(mechanical.backward_euler_loss_j, 2.0 * (q2 - q1).powi(2));
        near(
            dt * mechanical.power.supplied_body_power_w,
            energy2 + 2.0 * (q2 * q2 - q1 * q1) + mechanical.backward_euler_loss_j,
        );
        assert_eq!(history.unsteady_time_s(), Some(0.11));
        near(next.mass_kg(), original.mass_kg());
        outcomes.push(q2);
    }
    assert!(outcomes[1] < outcomes[0] && outcomes[2] > outcomes[0]);
}

#[test]
fn unsupported_actual_phase_shear_rolls_back_heat_native_motion_and_clock() {
    let card = source(20.0, 0.2);
    let curve = super::continuity::chart(&card);
    let law = MushyPhaseLaw::try_new_with_immobile_solid(&curve, 1.0, 1.0, 1.0).unwrap();
    let domain = VoxelDomain::new(2, 2, 1, 1.0).unwrap();
    let ids: Vec<_> = (0..4).map(|c| Some(CellId(c))).collect();
    let config = closed();
    let gate = CancelGate::new_clock_free();
    let original = cold(&domain, &curve);
    let source = |_: usize, cell: &DepositedCell, _: f64, rate: f64| {
        resolve_at(&card, &law, cell.specific_h_j_kg, rate)
    };
    let history = bound(
        &original, &domain, &ids, &config, &gate, &card, &law, &source,
    );
    let (molten, history, _) = original
        .advance_observed(&history, 0.11, step(), || Ok(()))
        .unwrap();
    let history = history
        .with_unsteady_body_force_heating(forces(), controls())
        .unwrap();
    let before = (
        molten.cells().clone(),
        history.unsteady_flow().cloned(),
        history.unsteady_heating_audit(),
    );
    // The accepted-start rate is zero and supported; the solved endpoint's
    // sqrt(6)*400/(40+2400) exceeds 0.2/s and must refuse the whole interval.
    assert!(
        molten
            .advance_observed(&history, 0.21, step(), || Ok(()))
            .is_err()
    );
    assert_eq!(molten.cells(), &before.0);
    assert_eq!(history.unsteady_flow(), before.1.as_ref());
    assert_eq!(history.unsteady_heating_audit(), before.2);
    assert_eq!(history.unsteady_time_s(), Some(0.11));
    assert_eq!(history.endpoint_shear_rate_per_s().unwrap(), &[0.0; 4]);
    // Refining from the same accepted prefix falls inside the SAME source.
    let (_, resumed, _) = molten
        .advance_observed(&history, 0.111, step(), || Ok(()))
        .unwrap();
    assert!(
        resumed
            .endpoint_shear_rate_per_s()
            .unwrap()
            .iter()
            .all(|rate| *rate > 0.0 && *rate < 0.2)
    );
}
