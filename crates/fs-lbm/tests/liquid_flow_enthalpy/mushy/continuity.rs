//! Original current, thermal and momentum owners survive exact solid endpoints.
//! Manufactured source domains deliberately do not define liquid Cp or mu at phi=0.
use super::*;
use fs_lbm::conjugate::enthalpy::advection::feedback::electrical::{
    MaterialFlowElectricalError, advance_joule_enthalpy_with_material_flow,
};
use fs_lbm::conjugate::enthalpy::advection::feedback::{FlowMaterialSample, PhaseFlowTransport};
use fs_lbm::conjugate::enthalpy::{JouleEnthalpyError, JouleEnthalpyInput};
use fs_lbm::conjugate::{ChtError, ElectricSetup, Electrode, ElectrodeDrive, ViscousHeatControls};
use fs_material::MaterialError;
use fs_material::conductor::ElectricalPropertyBasis;
use fs_material::phase::transport::{
    PhaseTransportSourceQuery, ResolvedPhaseTransportState, resolve_phase_transport_state,
};
use fs_material::thermal::source::ThermalMaterialPoint;

#[path = "continuity/electrothermal_momentum.rs"]
mod electrothermal_momentum;

#[path = "continuity/anchored_mhd.rs"]
mod anchored_mhd;

pub(super) fn source(mu: f64) -> MaterialCard {
    let mut claims = ClaimSet::new();
    for (name, dims, value) in [
        (DENSITY_PROPERTY, Density::DIMS, 1.0),
        (
            SPECIFIC_HEAT_CAPACITY_PROPERTY,
            SPECIFIC_HEAT_CAPACITY_DIMS,
            2.0,
        ),
        (DYNAMIC_VISCOSITY_PROPERTY, DYNAMIC_VISCOSITY_DIMS, mu),
        (
            THERMAL_CONDUCTIVITY_PROPERTY,
            THERMAL_CONDUCTIVITY_DIMS,
            0.1,
        ),
        (
            ElectricalPropertyBasis::Resistivity.property_name(),
            ElectricalPropertyBasis::Resistivity.dimensions(),
            1.0,
        ),
    ] {
        let mut claim = scalar(name, dims, value);
        let liquid_only =
            name == SPECIFIC_HEAT_CAPACITY_PROPERTY || name == DYNAMIC_VISCOSITY_PROPERTY;
        claim.validity = claim.validity.with(
            "liquid",
            if liquid_only { f64::from_bits(1) } else { 0.0 },
            1.0,
        );
        if name == THERMAL_CONDUCTIVITY_PROPERTY
            || name == ElectricalPropertyBasis::Resistivity.property_name()
        {
            claim.value = PropertyValue::Curve {
                abscissa: "liquid".into(),
                abscissa_dims: Dims([0; 6]),
                knots: vec![(0.0, value), (1.0, 2.0 * value)],
                dims,
            };
            claim.interpolation = InterpolationPolicy::LinearInside;
        }
        claims.insert_claim(claim).unwrap();
    }
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: "manufactured phase-change material".into(),
            phase: "solid-liquid".into(),
            process: "numerical coupling fixture".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap()
}

pub(super) fn chart(card: &MaterialCard) -> EquilibriumEnthalpyPhaseCurve {
    EquilibriumEnthalpyPhaseCurve::try_new(
        card.content_hash(),
        [
            (-100.0, 450.0, 0.0),
            (0.0, 500.0, 0.0),
            (100.0, 500.0, 1.0),
            (300.0, 600.0, 1.0),
        ]
        .into_iter()
        .map(|(h, t, phi)| EnthalpyPhaseKnot {
            specific_enthalpy_j_kg: h,
            temperature_k: t,
            liquid_mass_fraction: phi,
            bulk_density_kg_m3: 1.0,
        })
        .collect(),
    )
    .unwrap()
}

pub(super) fn resolve(
    card: &MaterialCard,
    law: &MushyPhaseLaw<'_>,
    h: f64,
) -> Result<ResolvedPhaseTransportState, MaterialError> {
    let context = point(500.0, 100000.0).with("liquid", 0.0).unwrap();
    let extras = [ElectricalPropertyBasis::Resistivity.requirement().unwrap()];
    resolve_phase_transport_state(
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
            extra_requirements: &extras,
            selection: MaterialPropertySelection::SingleClaimOnly,
        },
    )
}

fn initial<'a>(
    domain: &VoxelDomain,
    curve: &'a EquilibriumEnthalpyPhaseCurve,
) -> DepositionSimulation<'a> {
    let mut events = Vec::new();
    for voxel in 0..domain.cell_count() {
        let [x, y, z] = domain.coords(voxel);
        let center = domain.center(x, y, z);
        let mut start = center;
        let mut end = center;
        start[0] -= 0.5;
        end[0] += 0.5;
        events.push(TimedEvent {
            time_s: 0.0,
            event: Event::Birth(DepositedCell {
                id: CellId(voxel as u64),
                segment_id: voxel as u64,
                body_id: 1,
                material_id: MAT,
                mass_kg: 1.0,
                specific_h_j_kg: -20.0,
                geometry: CellGeometry {
                    start_m: start,
                    end_m: end,
                    volume_m3: 1.0,
                    surface_area_m2: 6.0,
                },
                boundary: Boundary {
                    ambient_k: 490.0,
                    htc_w_m2_k: 0.0,
                    reservoir: None,
                    source_w: 0.0,
                },
            }),
        });
    }
    let mut link = 0;
    for y in 0..3 {
        for x in 0..3 {
            for [xx, yy] in [[x + 1, y], [x, y + 1]] {
                if xx < 3 && yy < 3 {
                    events.push(TimedEvent {
                        time_s: 0.0,
                        event: Event::Contact(ThermalLink {
                            id: link,
                            a: CellId((3 * y + x) as u64),
                            b: CellId((3 * yy + xx) as u64),
                            area_m2: 1.0,
                            conductance_w_k: 0.0,
                        }),
                    });
                    link += 1;
                }
            }
        }
    }
    DepositionSimulation::new(
        events,
        &[(MAT, curve as &dyn ThermalLaw)],
        DepositionBudget {
            max_cells: 9,
            max_links: 12,
            max_events: 21,
            max_substeps: 2000,
        },
    )
    .unwrap()
    .advance_to(0.0, step(), || Ok(()))
    .unwrap()
    .0
}

#[derive(Clone, Default)]
struct CycleHeat;
impl DepositionObserver for CycleHeat {
    type Error = DepositionError;
    fn sources(
        &mut self,
        thermal: &DepositionSimulation<'_>,
        _: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<BTreeMap<CellId, f64>, Self::Error> {
        // Each drive change is also an explicit caller-owned coupling endpoint.
        let watts = if thermal.time_s() < 1.5 {
            100.0
        } else if thermal.time_s() < 1.6 {
            0.0
        } else {
            -200.0
        };
        Ok(thermal.cells().keys().map(|&id| (id, watts)).collect())
    }
    fn event(&mut self, _: &Event, _: &DepositionSimulation<'_>) -> Result<(), Self::Error> {
        Ok(())
    }
    fn step(
        &mut self,
        _: f64,
        _: f64,
        _: &BTreeMap<CellId, f64>,
        _: &DepositionSimulation<'_>,
        _: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), Self::Error> {
        Ok(())
    }
}

fn step() -> StepConfig {
    StepConfig {
        max_dt_s: 0.025,
        energy_tolerance_j: 1e-8,
    }
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
fn setup(voltage: f64) -> ElectricSetup {
    ElectricSetup {
        conductors: vec![],
        electrodes: vec![
            Electrode {
                lo: [0, 0, 0],
                hi: [1, 3, 1],
                drive: ElectrodeDrive::Voltage(voltage),
            },
            Electrode {
                lo: [2, 0, 0],
                hi: [3, 3, 1],
                drive: ElectrodeDrive::Voltage(0.0),
            },
        ],
    }
}
fn current_input<'a>(
    domain: &'a VoxelDomain,
    setup: &'a ElectricSetup,
    ids: &'a [CellId],
) -> JouleEnthalpyInput<'a> {
    JouleEnthalpyInput {
        domain,
        setup,
        contacts: &[],
        cell_ids: ids,
        max_cells: 9,
        max_current_solves: 1000,
        power_tolerance_w: 1e-8,
        energy_tolerance_j: 1e-8,
    }
}
fn resistivity(
    sample: &FlowMaterialSample<ResolvedPhaseTransportState, PhaseFlowTransport>,
) -> Result<Option<f64>, ChtError> {
    Ok(Some(
        sample
            .material
            .material()
            .property(ElectricalPropertyBasis::Resistivity.property_name())
            .unwrap()
            .value_si(),
    ))
}

#[test]
fn the_same_source_and_cells_conduct_current_melt_flow_and_freeze_without_liquid_queries_on_solids()
{
    let domain = VoxelDomain::new(3, 3, 1, 1.0).unwrap();
    let ids = std::array::from_fn::<_, 9, _>(|c| CellId(c as u64));
    let flow_ids = ids.map(Some);
    let gate = CancelGate::new();
    let config = closed();
    let on = setup(1.0);
    let off = setup(0.0);
    let links: Vec<_> = (0..12).collect();
    let mut fluxes = Vec::new();
    for mu in [100.0, 200.0] {
        let card = source(mu);
        let curve = chart(&card);
        let law = MushyPhaseLaw::try_new_with_immobile_solid(&curve, 1.0, 1.0, 1.0).unwrap();
        let old = initial(&domain, &curve);
        let initial_energy = old.energy_j();
        let input = FlowEnthalpyInput {
            cell_ids: &flow_ids,
            material_id: MAT,
            caloric: &curve,
            inlet_specific_h_j_kg: [None; 6],
            max_cells: 9,
            max_connections: 64,
        };
        let observer = MaterialFlowFeedback::try_new_phase(
            &old,
            &domain,
            input,
            &config,
            CycleHeat,
            &gate,
            1000,
            &law,
            |_, cell: &DepositedCell, _| resolve(&card, &law, cell.specific_h_j_kg),
        )
        .unwrap()
        .with_fluid_conduction(&old, &links)
        .unwrap()
        .with_viscous_heating(controls())
        .unwrap();
        assert_eq!(observer.endpoint().len(), 9);
        for sample in observer.endpoint() {
            assert!(sample.transport.mobile_transport().is_none());
            assert!(
                sample
                    .material
                    .material()
                    .property(DYNAMIC_VISCOSITY_PROPERTY)
                    .is_none()
            );
            assert!(
                sample
                    .material
                    .material()
                    .property(SPECIFIC_HEAT_CAPACITY_PROPERTY)
                    .is_none()
            );
            assert!(
                sample
                    .material
                    .material()
                    .property(ElectricalPropertyBasis::Resistivity.property_name())
                    .is_some()
            );
        }
        let (solid, solid_flow, solid_report) = advance_joule_enthalpy_with_material_flow(
            &old,
            0.1,
            step(),
            current_input(&domain, &on, &ids),
            &observer,
            resistivity,
        )
        .unwrap();
        assert!(
            solid_flow
                .endpoint()
                .iter()
                .all(|s| s.transport.mobile_transport().is_none())
        );
        assert!(
            solid_flow
                .applied()
                .unwrap()
                .flow
                .velocity_m_s
                .iter()
                .flatten()
                .all(|&u| u == 0.0)
        );
        assert!(solid_report.electrode_energy_j > 0.0);
        assert!(solid_report.endpoint.electric.electrodes[0].current_a > 0.0);
        near(
            solid.energy_j() - initial_energy,
            90.0 + solid_report.electrode_energy_j,
        );

        let (melt, melted_flow, melt_report) = advance_joule_enthalpy_with_material_flow(
            &solid,
            1.5,
            step(),
            current_input(&domain, &on, &ids),
            &solid_flow,
            resistivity,
        )
        .unwrap();
        for sample in melted_flow.endpoint() {
            assert_eq!(sample.material.phase_state().liquid_mass_fraction(), 1.0);
            assert!(sample.transport.mobile_transport().is_some());
            assert_eq!(
                sample.material.material().card_identity(),
                card.content_hash()
            );
        }
        near(
            melt.energy_j() - initial_energy,
            1350.0 + solid_report.electrode_energy_j + melt_report.electrode_energy_j,
        );

        // Stir only the already mobile region. No force is silently discarded
        // on an immobile skeleton, and no momentum/thermal owner is replaced.
        let force = vec![
            [-1.0, 1.0, 0.0],
            [-1.0, 0.0, 0.0],
            [-1.0, -1.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0; 3],
            [0.0, -1.0, 0.0],
            [1.0, 1.0, 0.0],
            [1.0, 0.0, 0.0],
            [1.0, -1.0, 0.0],
        ];
        let driven = melted_flow
            .with_body_force_heating(force, controls())
            .unwrap();
        let (moving, moving_flow, moving_report) = advance_joule_enthalpy_with_material_flow(
            &melt,
            1.6,
            step(),
            current_input(&domain, &on, &ids),
            &driven,
            resistivity,
        )
        .unwrap();
        let maximum = moving_flow
            .applied()
            .unwrap()
            .flow
            .velocity_m_s
            .iter()
            .flatten()
            .map(|v| v.abs())
            .fold(0.0_f64, f64::max);
        assert!(maximum > 0.0);
        let work = moving_flow.viscous_heating_audit().unwrap();
        assert!(work.generated_heat_j > 0.0);
        near(work.generated_heat_j, work.supplied_body_work_j);
        near(
            moving.energy_j() - melt.energy_j(),
            moving_report.electrode_energy_j + work.generated_heat_j,
        );
        fluxes.push(maximum);

        let cooling = moving_flow
            .with_body_force_heating(vec![[0.0; 3]; 9], controls())
            .unwrap();
        let (frozen, final_flow, final_report) = advance_joule_enthalpy_with_material_flow(
            &moving,
            2.6,
            step(),
            current_input(&domain, &off, &ids),
            &cooling,
            resistivity,
        )
        .unwrap();
        for sample in final_flow.endpoint() {
            assert!(sample.material.phase_state().is_fully_solid());
            assert!(sample.transport.mobile_transport().is_none());
            assert!(
                sample
                    .material
                    .material()
                    .property(DYNAMIC_VISCOSITY_PROPERTY)
                    .is_none()
            );
            assert!(
                sample
                    .material
                    .material()
                    .property(SPECIFIC_HEAT_CAPACITY_PROPERTY)
                    .is_none()
            );
        }
        assert!(
            final_flow
                .applied()
                .unwrap()
                .flow
                .velocity_m_s
                .iter()
                .flatten()
                .all(|&u| u == 0.0)
        );
        assert_eq!(final_report.electrode_energy_j, 0.0);
        near(frozen.energy_j() - moving.energy_j(), -1800.0);
        near(frozen.mass_kg(), old.mass_kg());
        assert_eq!(frozen.cells().len(), old.cells().len());
        assert_eq!(old.time_s(), 0.0);
        assert!(observer.applied().is_none());

        // Frozen does not mean electrically disconnected. Re-enable current
        // below the fusion endpoint and preserve the same thermal cells.
        let (powered_solid, still_solid, report) = advance_joule_enthalpy_with_material_flow(
            &frozen,
            2.625,
            step(),
            current_input(&domain, &on, &ids),
            &final_flow,
            resistivity,
        )
        .unwrap();
        assert!(report.electrode_energy_j > 0.0);
        assert!(
            still_solid
                .endpoint()
                .iter()
                .all(|s| s.transport.mobile_transport().is_none())
        );
        near(
            powered_solid.energy_j() - frozen.energy_j(),
            -45.0 + report.electrode_energy_j,
        );
    }
    // Only the material card changed; the mechanical law and drive did not.
    near(fluxes[0], 2.0 * fluxes[1]);
}

#[test]
fn endpoint_current_refusal_publishes_no_part_of_a_solid_to_mobile_transition() {
    let card = source(100.0);
    let curve = chart(&card);
    let law = MushyPhaseLaw::try_new_with_immobile_solid(&curve, 1.0, 1.0, 1.0).unwrap();
    let domain = VoxelDomain::new(3, 3, 1, 1.0).unwrap();
    let old = initial(&domain, &curve);
    let ids = std::array::from_fn::<_, 9, _>(|c| CellId(c as u64));
    let flow_ids = ids.map(Some);
    let gate = CancelGate::new();
    let config = closed();
    let setup = setup(1.0);
    let observer = MaterialFlowFeedback::try_new_phase(
        &old,
        &domain,
        FlowEnthalpyInput {
            cell_ids: &flow_ids,
            material_id: MAT,
            caloric: &curve,
            inlet_specific_h_j_kg: [None; 6],
            max_cells: 9,
            max_connections: 64,
        },
        &config,
        CycleHeat,
        &gate,
        1000,
        &law,
        |_, cell: &DepositedCell, _| resolve(&card, &law, cell.specific_h_j_kg),
    )
    .unwrap();
    let mut limited = current_input(&domain, &setup, &ids);
    limited.max_current_solves = 1;
    let one = StepConfig {
        max_dt_s: 0.3,
        energy_tolerance_j: 1e-8,
    };
    assert!(matches!(
        advance_joule_enthalpy_with_material_flow(&old, 0.3, one, limited, &observer, resistivity),
        Err(MaterialFlowElectricalError::Electrical(
            JouleEnthalpyError::CurrentBudget
        ))
    ));
    assert_eq!(old.time_s(), 0.0);
    assert!(old.cells().values().all(|c| c.specific_h_j_kg == -20.0));
    assert!(observer.applied().is_none());
    assert_eq!(observer.momentum_solves(), 0);
    let (a, flow_a, report_a) = advance_joule_enthalpy_with_material_flow(
        &old,
        0.3,
        one,
        current_input(&domain, &setup, &ids),
        &observer,
        resistivity,
    )
    .unwrap();
    let (b, flow_b, report_b) = advance_joule_enthalpy_with_material_flow(
        &old,
        0.3,
        one,
        current_input(&domain, &setup, &ids),
        &observer,
        resistivity,
    )
    .unwrap();
    assert_eq!(a.cells(), b.cells());
    assert_eq!(flow_a.endpoint(), flow_b.endpoint());
    assert_eq!(report_a, report_b);
    assert!(
        flow_a
            .endpoint()
            .iter()
            .all(|s| s.transport.mobile_transport().is_some())
    );
}
