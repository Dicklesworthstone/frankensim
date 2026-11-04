//! G1/G3/G4: original paired conduction, physical drag heat and source bounds.

use super::*;
use crate::conjugate::flow::{
    CoexistenceVolumeControls, CoexistenceVolumeInput, DensityForceControls,
    DensityMomentumReconstructionControls, DensityProjectionControls, MaterialAdvectionControls,
    MaterialDensityMomentumControls, MaterialMomentumControls, PiecewiseThermalVolumeControls,
    ThermalVolumeControls, ThermalVolumeMomentumControls,
};
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSelection, ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim,
    PropertyKey, PropertyValue, Provenance, QueryPoint, SelectionPolicy, UncertaintyModel,
};
use fs_material::liquid::{KINEMATIC_VISCOSITY_DIMS, KINEMATIC_VISCOSITY_PROPERTY};
use fs_material::liquid_temperature::LiquidViscosityBasis;
use fs_material::phase::transport::{
    PhaseTransportSourceQuery, resolve_phase_transport_state_for_volume_change,
};
use fs_material::phase::{EnthalpyPhaseKnot, EquilibriumEnthalpyPhaseCurve};
use fs_material::state_point::{
    DENSITY_PROPERTY, MaterialPropertySelection, SPECIFIC_HEAT_CAPACITY_DIMS,
    SPECIFIC_HEAT_CAPACITY_PROPERTY, THERMAL_CONDUCTIVITY_PROPERTY,
};
use fs_material::thermal::source::ThermalMaterialPoint;
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_qty::{Density, Dims, QuantitySpec};

fn quantity(kind: QuantityKind) -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(kind, ValueForm::Static))
}

fn fixture() -> (MaterialCard, EquilibriumEnthalpyPhaseCurve) {
    let mut claims = ClaimSet::new();
    for (name, spec, value) in [
        (
            DENSITY_PROPERTY,
            QuantitySpec::dimensional(Density::DIMS),
            0.0,
        ),
        (
            THERMAL_CONDUCTIVITY_PROPERTY,
            quantity(QuantityKind::ThermalConductivity),
            0.1,
        ),
        (
            SPECIFIC_HEAT_CAPACITY_PROPERTY,
            QuantitySpec::dimensional(SPECIFIC_HEAT_CAPACITY_DIMS),
            1.0,
        ),
        (
            KINEMATIC_VISCOSITY_PROPERTY,
            QuantitySpec::dimensional(KINEMATIC_VISCOSITY_DIMS),
            0.5,
        ),
    ] {
        let density = name == DENSITY_PROPERTY;
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::with_quantity(name, spec),
                // Explicit manufactured identity coordinate, not measured data.
                // Power-of-two interpolation retains arbitrary current chart rho
                // exactly; a linear rho(f) table would change the mixture model.
                value: if density {
                    PropertyValue::Curve {
                        abscissa: "manufactured_mass_concentration".into(),
                        abscissa_dims: Density::DIMS,
                        knots: vec![(0.0, 0.0), (8.0, 8.0)],
                        dims: Density::DIMS,
                    }
                } else {
                    PropertyValue::Scalar {
                        value,
                        dims: spec.dims(),
                    }
                },
                validity: ValidityDomain::unconstrained()
                    .with_quantity(
                        "T",
                        quantity(QuantityKind::AbsoluteTemperature),
                        400.0,
                        600.0,
                    )
                    .with_quantity("liquid", QuantitySpec::dimensional(Dims([0; 6])), 0.0, 1.0)
                    .with_quantity(
                        "manufactured_mass_concentration",
                        quantity(QuantityKind::MassConcentration),
                        1.0,
                        8.0,
                    ),
                uncertainty: UncertaintyModel::Unstated,
                interpolation: if density {
                    InterpolationPolicy::LinearInside
                } else {
                    InterpolationPolicy::ConstantWithinValidity
                },
                observations: vec![],
                provenance: Provenance {
                    source: "MANUFACTURED conductive Darcy/phase fixture; not measurements".into(),
                    license: "test-only".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "manufactured phase mixture".into(),
            phase: "solid/liquid".into(),
            process: "original conductive heat and coupled density verification".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    let curve = EquilibriumEnthalpyPhaseCurve::try_new(
        card.content_hash(),
        [
            (-100.0, 400.0, 0.0, 4.0),
            (0.0, 500.0, 0.0, 4.0),
            (100.0, 500.0, 1.0, 2.0),
            (200.0, 600.0, 1.0, 2.0),
        ]
        .into_iter()
        .map(|(h, t, liquid, rho)| EnthalpyPhaseKnot {
            specific_enthalpy_j_kg: h,
            temperature_k: t,
            liquid_mass_fraction: liquid,
            bulk_density_kg_m3: rho,
        })
        .collect(),
    )
    .unwrap();
    (card, curve)
}

fn context(density: f64) -> QueryPoint {
    QueryPoint::new()
        .with_quantity("T", quantity(QuantityKind::AbsoluteTemperature), 550.0)
        .unwrap()
        .with_quantity("liquid", QuantitySpec::dimensional(Dims([0; 6])), 0.9)
        .unwrap()
        .with_quantity(
            "manufactured_mass_concentration",
            quantity(QuantityKind::MassConcentration),
            density,
        )
        .unwrap()
        .with_quantity("pressure", quantity(QuantityKind::Pressure), 100000.0)
        .unwrap()
}

fn conductivity<'a>(
    card: &MaterialCard,
    curve: &'a EquilibriumEnthalpyPhaseCurve,
    bounds: [f64; 2],
) -> PreparedPhaseConductivity<'a> {
    let context = context(3.0);
    PreparedPhaseConductivity::new(
        curve,
        ThermalMaterialPoint {
            material: card,
            context: &context,
            temperature_axis: "T",
            liquid_fraction_axis: Some("liquid"),
        },
        &PropertyKey::with_quantity(
            THERMAL_CONDUCTIVITY_PROPERTY,
            quantity(QuantityKind::ThermalConductivity),
        ),
        ClaimSelection::Policy(SelectionPolicy::SingleClaimOnly),
        bounds,
        || Ok(()),
    )
    .unwrap()
}

fn source(
    card: &MaterialCard,
    law: &MushyPhaseLaw<'_>,
    phase: &EquilibriumPhaseState,
) -> Result<ResolvedPhaseTransportState, ChtError> {
    let context = context(phase.bulk_density_kg_m3());
    resolve_phase_transport_state_for_volume_change(
        law,
        phase.specific_enthalpy_j_kg(),
        PhaseTransportSourceQuery {
            source: ThermalMaterialPoint {
                material: card,
                context: &context,
                temperature_axis: "T",
                liquid_fraction_axis: Some("liquid"),
            },
            viscosity_basis: LiquidViscosityBasis::Kinematic,
            extra_requirements: &[],
            selection: MaterialPropertySelection::SingleClaimOnly,
        },
    )
    .map_err(|error| refuse(error.to_string()))
}

fn controls(dt: f64) -> ConductiveResistiveThermalVolumeMomentumControls {
    let projection = DensityProjectionControls {
        duration_s: dt,
        relative_tolerance: 1e-12,
        maximum_volume_defect_m3_s: 1e-10,
        maximum_energy_defect_j: 1e-7,
    };
    ConductiveResistiveThermalVolumeMomentumControls {
        resistance: ResistiveThermalVolumeMomentumControls {
            thermal: ThermalVolumeMomentumControls {
                volume: PiecewiseThermalVolumeControls {
                    thermal: ThermalVolumeControls {
                        volume: CoexistenceVolumeControls {
                            transport: MaterialAdvectionControls {
                                projection,
                                max_cells: 16,
                                max_connections: 128,
                                max_iterations: 200,
                                relative_tolerance: 1e-12,
                                maximum_mass_defect_kg: 1e-10,
                                maximum_enthalpy_defect_j: 1e-7,
                            },
                            maximum_cell_volume_defect_m3: 1e-10,
                        },
                        max_volume_iterations: 64,
                        source_relaxation: 1.0,
                    },
                    max_mixing_intervals: 4096,
                },
                mechanics: MaterialDensityMomentumControls {
                    reconstruction: DensityMomentumReconstructionControls {
                        max_cells: 16,
                        maximum_mass_defect_kg: 1e-10,
                        maximum_momentum_defect_kg_m_s: 1e-9,
                        maximum_energy_defect_j: 1e-7,
                    },
                    forcing: DensityForceControls {
                        max_cells: 16,
                        maximum_impulse_defect_kg_m_s: 1e-9,
                        maximum_energy_defect_j: 1e-7,
                    },
                    projection,
                    maximum_feedback_momentum_defect_kg_m_s: 1e-9,
                    maximum_feedback_angular_impulse_defect_kg_m2_s: 1e-9,
                    maximum_feedback_energy_defect_j: 1e-7,
                },
                momentum: MaterialMomentumControls {
                    max_iterations: 200,
                    relative_tolerance: 1e-12,
                    maximum_momentum_defect_kg_m_s: 1e-9,
                    maximum_energy_defect_j: 1e-7,
                },
            },
            max_coupling_iterations: 64,
            heat_relaxation: 1.0,
            maximum_cell_heat_defect_j: 1e-11,
        },
        max_conductive_faces: 32,
    }
}

fn input<'a>(
    mass: &'a [f64],
    phases: &'a [Option<EquilibriumPhaseState>],
    heat: &'a [f64],
    momentum: &'a [[f64; 3]],
    force: &'a [[f64; 3]],
) -> ThermalVolumeMomentumInput<'a> {
    ThermalVolumeMomentumInput {
        volume: CoexistenceVolumeInput {
            mass_kg: mass,
            phase_states: phases,
            net_enthalpy_increment_j: heat,
            boundary_inflow: [None; 6],
        },
        momentum_kg_m_s: momentum,
        acceleration_m_s2: [0.0; 3],
        force_density_n_m3: force,
        boundary_velocity_m_s: [None; 6],
    }
}

fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() < 5e-10 * expected.abs().max(1.0),
        "{actual:.17e} != {expected:.17e}"
    );
}

#[test]
fn paired_original_heat_drives_mushy_density_flow_and_one_physical_drag_account() {
    let (card, curve) = fixture();
    let conductivity = conductivity(&card, &curve, [-100.0, 200.0]);
    let law = MushyPhaseLaw::try_new_with_immobile_solid(&curve, 4.0, 2.0, 1.0).unwrap();
    let domain = VoxelDomain::new(2, 1, 1, 1.0).unwrap();
    let field = FlowField::quiescent(&domain);
    let mut faces = [FlowFace::Wall; 6];
    faces[1] = FlowFace::Free;
    let mass = [8.0 / 3.0, 2.0];
    let phases = [50.0, 150.0].map(|h| Some(curve.state_at_specific_enthalpy(h).unwrap()));
    let zero = [[0.0; 3]; 2];
    let authored = [0.0; 2];
    let data = input(&mass, &phases, &authored, &zero, &zero);
    let result = field
        .advance_isotropic_conductive_phase_resistance_thermal_volume_momentum(
            &domain,
            faces,
            &conductivity,
            &law,
            data,
            controls(0.25),
            &CancelGate::new_clock_free(),
            |_, phase| source(&card, &law, phase),
        )
        .unwrap();

    // k=.1, half widths .5+.5 and A=1 give G=.1. Original temperatures
    // [500,550] give paired Q=dt G delta_T=1.25 J, once for all drag trials.
    assert_eq!(result.conductive_enthalpy_increment_j, [1.25, -1.25]);
    assert_eq!(
        result.conduction.conductivity_source_identity,
        conductivity.identity()
    );
    assert_eq!(result.conduction.conductive_faces, 1);
    close(result.conduction.exchanged_heat_j, 1.25);
    close(
        result.conduction.maximum_enthalpy_mixing_weight,
        1.0 / 160.0,
    );
    close(result.conduction.conduction_balance_defect_j, 0.0);
    let exact_conduction = [1.25, -1.25];
    let direct = field
        .advance_resolved_phase_resistance_thermal_volume_momentum(
            &domain,
            faces,
            &law,
            ThermalVolumeMomentumInput {
                volume: CoexistenceVolumeInput {
                    net_enthalpy_increment_j: &exact_conduction,
                    ..data.volume
                },
                ..data
            },
            controls(0.25).resistance,
            &CancelGate::new_clock_free(),
            |_, phase| source(&card, &law, phase),
        )
        .unwrap();
    assert_eq!(result.evolution, direct);
    let native = &result.evolution.evolution;
    let volume = &native.thermal.volume;
    let transport = &volume.transport;
    assert!(native.coupling.coupling_iterations > 1);
    assert!(native.resistance_pa_s_m2[0] > 0.0);
    assert_eq!(native.resistance_pa_s_m2[1], 0.0);
    assert!(native.resistance.cell_resistance_heat_j[0] > 1e-6);
    assert_eq!(native.resistance.cell_resistance_heat_j[1], 0.0);
    // Cold-cell v=1/4+h/400. It has no incoming route, so its physical
    // expansion source and implicit outgoing mass have independent formulas.
    let q = (1.25 + native.caloric_resistance_heat_j[0]) / 100.0;
    close(transport.projection.field.flux_x(1, 0, 0), q);
    close(transport.mass_kg[0], mass[0] / (1.0 + 0.25 * q));
    close(
        transport.specific_enthalpy_j_kg[0],
        50.0 + (1.25 + native.caloric_resistance_heat_j[0]) / mass[0],
    );
    close(
        native.resistance.cell_resistance_heat_j[0],
        0.25 * 0.5 * native.resistance_pa_s_m2[0] * q * q,
    );
    close(transport.mass_kg[1], 2.0);
    assert!(transport.projection.field.flux_x(2, 0, 0) > q);
    assert_eq!(transport.report.mass_routes, 2);
    close(transport.report.incoming_mass_kg, 0.0);
    close(
        transport.report.outgoing_mass_kg,
        mass.iter().sum::<f64>() - transport.mass_kg.iter().sum::<f64>(),
    );
    for (cell, source) in result.evolution.endpoint_sources.iter().enumerate() {
        let phase = volume.phase_states[cell].unwrap();
        close(transport.mass_kg[cell] / phase.bulk_density_kg_m3(), 1.0);
        let source = source.as_ref().unwrap();
        assert_eq!(*source.phase_state(), phase);
        assert_eq!(source.source_density_kg_m3(), phase.bulk_density_kg_m3());
        assert_eq!(source.material().card_identity(), card.content_hash());
    }
    close(
        volume.report.material_enthalpy_change_j,
        native.resistance.resistance_heat_j + transport.report.carried_enthalpy_j,
    );
    close(
        native.resistance.kinetic_energy_change_j
            + native.resistance.backward_euler_loss_j
            + native.resistance.resistance_heat_j,
        volume.report.required_source_pressure_work_j,
    );
    assert!(
        result.conduction.absolute_enthalpy_defect_j + native.resistance.absolute_energy_defect_j
            <= controls(0.25)
                .resistance
                .thermal
                .mechanics
                .forcing
                .maximum_energy_defect_j
    );
    assert_eq!(native.momentum.report.represented_impulse_kg_m_s, [0.0; 3]);

    let x = [None, None, Some(1e12)];
    let transverse = [None; 4];
    let shifted = field.advance_isotropic_conductive_phase_resistance_thermal_volume_momentum_with_pressure_profile(
        &domain, faces, StaggeredBoundaryPressure { x_pa: &x, y_pa: &transverse, z_pa: &transverse },
        &conductivity, &law, data, controls(0.25), &CancelGate::new_clock_free(),
        |_, phase| source(&card, &law, phase)).unwrap();
    assert_eq!(shifted.conduction, result.conduction);
    assert_eq!(
        shifted.evolution.endpoint_sources,
        result.evolution.endpoint_sources
    );
    assert_eq!(shifted.evolution.evolution.momentum, native.momentum);
    assert_eq!(shifted.evolution.evolution.resistance, native.resistance);
    assert_eq!(
        shifted.evolution.evolution.thermal.volume.phase_states,
        volume.phase_states
    );
    assert_eq!(field, FlowField::quiescent(&domain));
}

#[test]
fn original_source_stability_and_drag_heated_intermediate_bounds_gate_publication() {
    let (card, curve) = fixture();
    let original_source = conductivity(&card, &curve, [-100.0, 200.0]);
    let law = MushyPhaseLaw::try_new_with_immobile_solid(&curve, 4.0, 2.0, 1.0).unwrap();
    let domain = VoxelDomain::new(2, 1, 1, 1.0).unwrap();
    let field = FlowField::quiescent(&domain);
    let mass = [8.0 / 3.0, 2.0];
    let phases = [50.0, 150.0].map(|h| Some(curve.state_at_specific_enthalpy(h).unwrap()));
    let zero = [[0.0; 3]; 2];
    let authored = [0.0; 2];
    let data = input(&mass, &phases, &authored, &zero, &zero);
    let mut faces = [FlowFace::Wall; 6];
    faces[1] = FlowFace::Free;
    let mut queries = 0;
    let unstable = field.advance_isotropic_conductive_phase_resistance_thermal_volume_momentum(
        &domain,
        faces,
        &original_source,
        &law,
        data,
        controls(64.0),
        &CancelGate::new_clock_free(),
        |_, phase| {
            queries += 1;
            source(&card, &law, phase)
        },
    );
    assert!(
        matches!(unstable, Err(ChtError::InvalidInput { reason, .. })
        if reason.contains("explicit conductive enthalpy mixing weight exceeds one"))
    );
    assert_eq!(queries, 0);
    let mut knots = curve.knots().to_vec();
    knots[0].temperature_k = 399.0;
    let other_curve = EquilibriumEnthalpyPhaseCurve::try_new(card.content_hash(), knots).unwrap();
    let other_law =
        MushyPhaseLaw::try_new_with_immobile_solid(&other_curve, 4.0, 2.0, 1.0).unwrap();
    let mismatched = field.advance_isotropic_conductive_phase_resistance_thermal_volume_momentum(
        &domain,
        faces,
        &original_source,
        &other_law,
        data,
        controls(0.25),
        &CancelGate::new_clock_free(),
        |_, phase| {
            queries += 1;
            source(&card, &other_law, phase)
        },
    );
    assert!(
        matches!(mismatched, Err(ChtError::InvalidInput { reason, .. })
        if reason.contains("same original chart and material"))
    );
    assert_eq!(queries, 0);
    let gate = CancelGate::new_clock_free();
    gate.request();
    assert!(matches!(
        field.advance_isotropic_conductive_phase_resistance_thermal_volume_momentum(
            &domain,
            faces,
            &original_source,
            &law,
            data,
            controls(0.25),
            &gate,
            |_, phase| source(&card, &law, phase)
        ),
        Err(ChtError::Cancelled)
    ));
    assert_eq!(field, FlowField::quiescent(&domain));
    assert_eq!(mass, [8.0 / 3.0, 2.0]);

    // Throughflow mixes the drag-heated original inventory back into the
    // source range. Its valid returned phase cannot license the hotter
    // pre-transport caloric state; that state needs its own bound check.
    let single = VoxelDomain::new(1, 1, 1, 1.0).unwrap();
    let mut through_faces = [FlowFace::Wall; 6];
    through_faces[0] = FlowFace::Fixed;
    through_faces[1] = FlowFace::Free;
    let through =
        FlowField::from_face_velocity(
            &single,
            through_faces,
            |axis, _| if axis == 0 { 1.0 } else { 0.0 },
        )
        .unwrap();
    let original_through = through.clone();
    let one_mass = [8.0 / 3.0];
    let one_phase = [Some(curve.state_at_specific_enthalpy(50.0).unwrap())];
    let one_momentum = [[8.0 / 3.0, 0.0, 0.0]];
    let one_zero = [[0.0; 3]];
    let one_heat = [0.0];
    let mut through_input = input(&one_mass, &one_phase, &one_heat, &one_momentum, &one_zero);
    through_input.volume.boundary_inflow[0] = Some(curve.state_at_specific_enthalpy(49.0).unwrap());
    through_input.boundary_velocity_m_s[0] = Some([1.0, 0.0, 0.0]);
    let direct = through
        .advance_resolved_phase_resistance_thermal_volume_momentum(
            &single,
            through_faces,
            &law,
            through_input,
            controls(0.25).resistance,
            &CancelGate::new_clock_free(),
            |_, phase| source(&card, &law, phase),
        )
        .unwrap();
    let bounds = [49.0, 50.02];
    let caloric_h = 50.0 + direct.evolution.caloric_resistance_heat_j[0] / one_mass[0];
    assert!(caloric_h > bounds[1]);
    assert!(
        (bounds[0]..=bounds[1]).contains(
            &direct
                .evolution
                .thermal
                .volume
                .transport
                .specific_enthalpy_j_kg[0]
        )
    );
    let bounded_source = conductivity(&card, &curve, bounds);
    let rejected = through.advance_isotropic_conductive_phase_resistance_thermal_volume_momentum(
        &single,
        through_faces,
        &bounded_source,
        &law,
        through_input,
        controls(0.25),
        &CancelGate::new_clock_free(),
        |_, phase| source(&card, &law, phase),
    );
    assert!(
        matches!(rejected, Err(ChtError::InvalidInput { reason, .. })
        if reason.contains("drag-heated caloric state leaves the original conductivity source"))
    );
    let retried = through
        .advance_isotropic_conductive_phase_resistance_thermal_volume_momentum(
            &single,
            through_faces,
            &original_source,
            &law,
            through_input,
            controls(0.25),
            &CancelGate::new_clock_free(),
            |_, phase| source(&card, &law, phase),
        )
        .unwrap();
    assert_eq!(retried.evolution, direct);
    assert_eq!(retried.conductive_enthalpy_increment_j, [0.0]);
    assert_eq!(through, original_through);
    assert_eq!(one_mass, [8.0 / 3.0]);
    assert_eq!(
        one_phase,
        [Some(curve.state_at_specific_enthalpy(50.0).unwrap())]
    );
    assert_eq!(one_momentum, [[8.0 / 3.0, 0.0, 0.0]]);
}
