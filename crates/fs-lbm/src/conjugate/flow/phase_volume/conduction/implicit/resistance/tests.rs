//! Independent latent-cell bath/drag balance on the real source and flow owners.

use super::*;
use crate::conjugate::flow::{
    CoexistenceVolumeControls, DensityForceControls, DensityMomentumReconstructionControls,
    DensityProjectionControls, MaterialAdvectionControls, MaterialDensityMomentumControls,
    MaterialMomentumControls, ResistiveThermalVolumeMomentumControls, ThermalVolumeControls,
    ThermalVolumeMomentumControls,
};
use crate::d3q19::Face3;
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
    fixture_with_densities([2.0; 2])
}

fn fixture_with_densities(densities: [f64; 2]) -> (MaterialCard, EquilibriumEnthalpyPhaseCurve) {
    let mut claims = ClaimSet::new();
    for (name, spec, value) in [
        (
            DENSITY_PROPERTY,
            QuantitySpec::dimensional(Density::DIMS),
            2.0,
        ),
        (
            THERMAL_CONDUCTIVITY_PROPERTY,
            quantity(QuantityKind::ThermalConductivity),
            0.0,
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
        let conductivity = name == THERMAL_CONDUCTIVITY_PROPERTY;
        let density = name == DENSITY_PROPERTY;
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::with_quantity(name, spec),
                value: if density {
                    // Exact manufactured chart-density query, not an
                    // approximate linear mass-fraction density mixture.
                    PropertyValue::Curve {
                        abscissa: "manufactured_mass_concentration".into(),
                        abscissa_dims: Density::DIMS,
                        knots: vec![(0.0, 0.0), (8.0, 8.0)],
                        dims: Density::DIMS,
                    }
                } else if conductivity {
                    PropertyValue::Curve {
                        abscissa: "liquid".into(),
                        abscissa_dims: Dims([0; 6]),
                        knots: vec![(0.0, 0.125), (1.0, 0.375)],
                        dims: spec.dims(),
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
                interpolation: if conductivity || density {
                    InterpolationPolicy::LinearInside
                } else {
                    InterpolationPolicy::ConstantWithinValidity
                },
                observations: vec![],
                provenance: Provenance {
                    source: "MANUFACTURED latent implicit bath/drag balance; not measurements"
                        .into(),
                    license: "test-only".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "manufactured phase-density mixture".into(),
            phase: "solid/liquid".into(),
            process: "coupled implicit Fourier and Darcy calorimetry".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    let curve = EquilibriumEnthalpyPhaseCurve::try_new(
        card.content_hash(),
        [
            (-100.0, 400.0, 0.0),
            (0.0, 500.0, 0.0),
            (100.0, 500.0, 1.0),
            (200.0, 600.0, 1.0),
        ]
        .into_iter()
        .map(|(h, t, f)| EnthalpyPhaseKnot {
            specific_enthalpy_j_kg: h,
            temperature_k: t,
            liquid_mass_fraction: f,
            bulk_density_kg_m3: if f == 0.0 { densities[0] } else { densities[1] },
        })
        .collect(),
    )
    .unwrap();
    (card, curve)
}

fn context() -> QueryPoint {
    QueryPoint::new()
        .with_quantity("T", quantity(QuantityKind::AbsoluteTemperature), 500.0)
        .unwrap()
        .with_quantity("liquid", QuantitySpec::dimensional(Dims([0; 6])), 0.5)
        .unwrap()
        .with_quantity("pressure", quantity(QuantityKind::Pressure), 100000.0)
        .unwrap()
        .with_quantity(
            "manufactured_mass_concentration",
            quantity(QuantityKind::MassConcentration),
            2.0,
        )
        .unwrap()
}

fn source(
    card: &MaterialCard,
    law: &MushyPhaseLaw<'_>,
    phase: &EquilibriumPhaseState,
) -> Result<ResolvedPhaseTransportState, ChtError> {
    let context = context()
        .with_quantity(
            "manufactured_mass_concentration",
            quantity(QuantityKind::MassConcentration),
            phase.bulk_density_kg_m3(),
        )
        .unwrap();
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
    .map_err(|e| refuse(e.to_string()))
}

fn controls() -> ConductiveResistiveThermalVolumeMomentumControls {
    let projection = DensityProjectionControls {
        duration_s: 0.25,
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
                                max_cells: 1,
                                max_connections: 16,
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
                        max_cells: 1,
                        maximum_mass_defect_kg: 1e-10,
                        maximum_momentum_defect_kg_m_s: 1e-9,
                        maximum_energy_defect_j: 1e-7,
                    },
                    forcing: DensityForceControls {
                        max_cells: 1,
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
        max_conductive_faces: 1,
    }
}

fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() < 2e-9 * expected.abs().max(1.0),
        "{actual:.17e} != {expected:.17e}"
    );
}

mod viscous;

#[test]
fn latent_bath_drag_and_phase_share_one_accepted_energy_and_momentum_endpoint() {
    let (card, curve) = fixture();
    let law = MushyPhaseLaw::try_new_with_immobile_solid(&curve, 2.0, 2.0, 1.0).unwrap();
    let context = context();
    let conductivity = PreparedPhaseConductivity::new(
        &curve,
        ThermalMaterialPoint {
            material: &card,
            context: &context,
            temperature_axis: "T",
            liquid_fraction_axis: Some("liquid"),
        },
        &PropertyKey::with_quantity(
            THERMAL_CONDUCTIVITY_PROPERTY,
            quantity(QuantityKind::ThermalConductivity),
        ),
        ClaimSelection::Policy(SelectionPolicy::SingleClaimOnly),
        [40.0, 90.0],
        || Ok(()),
    )
    .unwrap();
    let domain = VoxelDomain::new(1, 1, 1, 1.0).unwrap();
    let field = FlowField::quiescent(&domain);
    let mut faces = [FlowFace::Wall; 6];
    faces[0] = FlowFace::Free;
    faces[1] = FlowFace::Free;
    let mass = [2.0];
    let phases = [Some(curve.state_at_specific_enthalpy(50.0).unwrap())];
    let momentum = [[0.0; 3]];
    // At f*=3/4, equal phase densities give phi=f. With nu=1/2 and rho=2,
    // mu=1, K0=1, beta*=mu(1-phi)^2/phi^3=4/27. The free-face momentum law
    // rho*u*=dt*(F-beta*u*) gives u*=1 when F=8+4/27 and dt=1/4.
    let beta = 4.0 / 27.0;
    let force = [[8.0 + beta, 0.0, 0.0]];
    // k(f*)=1/8+f*/4=5/16, A=1, d=1/2, Tb-T=16: Qbath=5/2.
    // Physical drag heat is dt*beta*u*^2*V=1/27; only these physical terms
    // join authored heat in M*(75-50)=50. BE kinetic loss stays mechanical.
    let authored = [50.0 - 2.5 - 1.0 / 27.0];
    let mut data = ThermalVolumeMomentumInput {
        volume: CoexistenceVolumeInput {
            mass_kg: &mass,
            phase_states: &phases,
            net_enthalpy_increment_j: &authored,
            boundary_inflow: [None; 6],
        },
        momentum_kg_m_s: &momentum,
        acceleration_m_s2: [0.0; 3],
        force_density_n_m3: &force,
        boundary_velocity_m_s: [None; 6],
    };
    // Inflow at h*=75 leaves the manufactured caloric endpoint unchanged
    // by conservative transport, while exercising real nonzero mass routes.
    data.volume.boundary_inflow[0] = Some(curve.state_at_specific_enthalpy(75.0).unwrap());
    data.boundary_velocity_m_s[0] = Some([1.0, 0.0, 0.0]);
    let baths = [ThermalVolumeBathFace {
        cell: 0,
        face: Face3::YMin,
        temperature_k: 516.0,
        resistance_m2_k_w: 0.0,
    }];
    let solve = ImplicitThermalVolumeSolveControls {
        cell_tolerance_j: 1e-12,
        ..ImplicitThermalVolumeSolveControls::default()
    };
    let mut workspace = ImplicitThermalVolumeWorkspace::new(1).unwrap();
    let mut queries = 0;
    let result = field
        .advance_isotropic_implicit_conductive_phase_resistance_thermal_volume_momentum_with_baths(
            &domain,
            faces,
            None,
            &conductivity,
            &baths,
            &law,
            data,
            controls(),
            solve,
            &mut workspace,
            &CancelGate::new_clock_free(),
            |_, phase| {
                queries += 1;
                source(&card, &law, phase)
            },
        )
        .unwrap();
    let native = &result.evolution.evolution;
    let volume = &native.thermal.volume;
    let transport = &volume.transport;
    close(
        volume.phase_states[0].unwrap().specific_enthalpy_j_kg(),
        75.0,
    );
    close(transport.mass_kg[0], 2.0);
    close(transport.projection.field.flux_x(0, 0, 0), 1.0);
    close(transport.projection.field.flux_x(1, 0, 0), 1.0);
    close(transport.report.incoming_mass_kg, 0.5);
    close(transport.report.outgoing_mass_kg, 0.5);
    close(transport.report.carried_enthalpy_j, 0.0);
    close(native.momentum.momentum_kg_m_s[0][0], 2.0);
    assert_eq!(native.momentum.report.represented_impulse_kg_m_s, [0.0; 3]);
    close(native.resistance_pa_s_m2[0], beta);
    close(native.caloric_resistance_heat_j[0], 1.0 / 27.0);
    close(native.resistance.resistance_heat_j, 1.0 / 27.0);
    close(native.resistance.kinetic_energy_change_j, 1.0);
    close(native.resistance.backward_euler_loss_j, 1.0);
    close(native.resistance.body_endpoint_work_j, 2.0 + 1.0 / 27.0);
    close(native.resistance.resistance_impulse_kg_m_s[0], -1.0 / 27.0);
    close(result.conductive_enthalpy_increment_j[0], 2.5);
    close(result.report.boundary_heat_j[2], 2.5);
    close(result.report.absolute_boundary_heat_j, 2.5);
    close(
        volume.report.material_enthalpy_change_j,
        authored[0] + native.resistance.resistance_heat_j + result.report.boundary_heat_j[2],
    );
    assert_eq!(result.report.conductive_faces, 0);
    assert_eq!(result.report.boundary_faces, 1);
    assert_eq!(
        result.report.conductivity_source_identity,
        conductivity.identity()
    );
    let endpoint = result.evolution.endpoint_sources[0].as_ref().unwrap();
    assert_eq!(*endpoint.phase_state(), volume.phase_states[0].unwrap());
    assert_eq!(endpoint.source_density_kg_m3(), 2.0);
    assert!(native.coupling.coupling_iterations > 2);
    assert_eq!(result.work.heat_trials, native.coupling.coupling_iterations);
    assert_eq!(queries, result.work.heat_trials + 1);
    assert!(result.work.residual_evaluations > result.report.solve.residual_evaluations);
    assert!(
        result.report.absolute_enthalpy_defect_j + native.resistance.absolute_energy_defect_j
            <= 1e-7
    );

    // The same native budgets include rejected thermal guesses. A ceiling
    // sufficient for only the accepted heat solve cannot fund all drag trials.
    let exhausted = ImplicitThermalVolumeSolveControls {
        maximum_residual_evaluations: result.report.solve.residual_evaluations,
        ..solve
    };
    assert!(field
        .advance_isotropic_implicit_conductive_phase_resistance_thermal_volume_momentum_with_baths(
            &domain, faces, None, &conductivity, &baths, &law, data,
            controls(), exhausted, &mut workspace, &CancelGate::new_clock_free(),
            |_, phase| source(&card, &law, phase),
        ).is_err());
    let late_gate = CancelGate::new_clock_free();
    let mut late_queries = 0;
    let cancelled = field
        .advance_isotropic_implicit_conductive_phase_resistance_thermal_volume_momentum_with_baths(
            &domain,
            faces,
            None,
            &conductivity,
            &baths,
            &law,
            data,
            controls(),
            solve,
            &mut workspace,
            &late_gate,
            |_, phase| {
                late_queries += 1;
                if late_queries == 2 {
                    late_gate.request();
                }
                source(&card, &law, phase)
            },
        );
    assert!(matches!(cancelled, Err(ChtError::Cancelled)));
    assert_eq!(late_queries, 2);
    let source_refusal = field
        .advance_isotropic_implicit_conductive_phase_resistance_thermal_volume_momentum_with_baths(
            &domain,
            faces,
            None,
            &conductivity,
            &baths,
            &law,
            data,
            controls(),
            solve,
            &mut workspace,
            &CancelGate::new_clock_free(),
            |_, phase| {
                if phase.specific_enthalpy_j_kg() > 74.9 {
                    return Err(refuse("manufactured endpoint source bound"));
                }
                source(&card, &law, phase)
            },
        );
    assert!(
        matches!(source_refusal, Err(ChtError::InvalidInput { reason, .. })
        if reason.contains("manufactured endpoint source bound"))
    );
    let retry = field
        .advance_isotropic_implicit_conductive_phase_resistance_thermal_volume_momentum_with_baths(
            &domain,
            faces,
            None,
            &conductivity,
            &baths,
            &law,
            data,
            controls(),
            solve,
            &mut workspace,
            &CancelGate::new_clock_free(),
            |_, phase| source(&card, &law, phase),
        )
        .unwrap();
    assert_eq!(retry, result);
    assert_eq!(field, FlowField::quiescent(&domain));
    assert_eq!(mass, [2.0]);
    assert_eq!(momentum, [[0.0; 3]]);
    assert_eq!(
        phases,
        [Some(curve.state_at_specific_enthalpy(50.0).unwrap())]
    );
}

#[test]
fn finite_wall_is_reheated_from_original_inventory_in_every_phase_drag_trial() {
    use crate::conjugate::{FaceContactResistance, Voxel};

    let (card, curve) = fixture();
    let law = MushyPhaseLaw::try_new_with_immobile_solid(&curve, 2.0, 2.0, 1.0).unwrap();
    let context = context();
    let prepare = |bounds| {
        PreparedPhaseConductivity::new(
            &curve,
            ThermalMaterialPoint {
                material: &card,
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
    };
    let conductivity = prepare([40.0, 90.0]);
    let wall_conductivity = prepare([110.0, 130.0]);
    let mut domain = VoxelDomain::new(1, 2, 1, 1.0).unwrap();
    domain.set(0, 1, 0, Voxel::Solid(0));
    let field = FlowField::quiescent(&domain);
    let mut faces = [FlowFace::Wall; 6];
    faces[0] = FlowFace::Free;
    faces[1] = FlowFace::Free;
    let mass = [2.0, 0.0];
    let phases = [Some(curve.state_at_specific_enthalpy(50.0).unwrap()), None];
    let momentum = [[0.0; 3]; 2];
    let beta = 4.0 / 27.0;
    let force = [[8.0 + beta, 0.0, 0.0], [0.0; 3]];
    let authored = [49.0 - 1.0 / 27.0, 0.0];
    let mut data = ThermalVolumeMomentumInput {
        volume: CoexistenceVolumeInput {
            mass_kg: &mass,
            phase_states: &phases,
            net_enthalpy_increment_j: &authored,
            boundary_inflow: [None; 6],
        },
        momentum_kg_m_s: &momentum,
        acceleration_m_s2: [0.0; 3],
        force_density_n_m3: &force,
        boundary_velocity_m_s: [None; 6],
    };
    data.volume.boundary_inflow[0] = Some(curve.state_at_specific_enthalpy(75.0).unwrap());
    data.boundary_velocity_m_s[0] = Some([1.0, 0.0, 0.0]);
    let solids = [ThermalVolumeSolidCell {
        cell: 1,
        conductivity: &wall_conductivity,
        mass_kg: 2.0,
        phase: curve.state_at_specific_enthalpy(116.5).unwrap(),
        net_enthalpy_increment_j: 0.0,
    }];
    let contacts = [ThermalVolumeContact::Resistance(FaceContactResistance {
        cells: [0, 1],
        resistance_m2_k_w: 16.0 / 15.0,
    })];
    // At f_fluid*=3/4, k_f=5/16 and k_s=3/8. The actual half cells and
    // contact give R=8/5+4/3+16/15=4 K/W. At T_wall*=516, T_f*=500,
    // dt=1/4 the wall supplies exactly 1 J, so M_s*(116-116.5)=-1 J.
    // The fluid gains 50 J from Q_authored=49-1/27, Q_drag=1/27 and Q_wall=1.
    // Its manufactured transported h=75 still gives beta=4/27 and u=1.
    let mut limits = controls();
    limits
        .resistance
        .thermal
        .volume
        .thermal
        .volume
        .transport
        .max_cells = 2;
    limits.resistance.thermal.mechanics.reconstruction.max_cells = 2;
    limits.resistance.thermal.mechanics.forcing.max_cells = 2;
    let solve = ImplicitThermalVolumeSolveControls {
        cell_tolerance_j: 1e-12,
        ..ImplicitThermalVolumeSolveControls::default()
    };
    let mut workspace = ImplicitThermalVolumeWorkspace::new(2).unwrap();
    let run = |workspace: &mut ImplicitThermalVolumeWorkspace, gate: &CancelGate| {
        field
            .advance_isotropic_implicit_conjugate_phase_resistance_thermal_volume_momentum_with_baths(
                &domain, faces, None, &conductivity, ThermalVolumeBoundaries::default(), &solids, &contacts,
                &law, data, limits, solve, workspace, gate,
                |_, phase| source(&card, &law, phase),
            )
    };
    let result = run(&mut workspace, &CancelGate::new_clock_free()).unwrap();
    let native = &result.evolution.evolution;
    let volume = &native.thermal.volume;
    close(
        volume.phase_states[0].unwrap().specific_enthalpy_j_kg(),
        75.0,
    );
    assert_eq!(volume.phase_states[1], None);
    assert_eq!(volume.transport.mass_kg, vec![2.0, 0.0]);
    close(native.momentum.momentum_kg_m_s[0][0], 2.0);
    assert_eq!(native.momentum.momentum_kg_m_s[1], [0.0; 3]);
    close(native.resistance_pa_s_m2[0], beta);
    close(native.resistance.resistance_heat_j, 1.0 / 27.0);
    close(native.resistance.backward_euler_loss_j, 1.0);
    close(result.conductive_enthalpy_increment_j[0], 1.0);
    assert_eq!(result.conductive_enthalpy_increment_j[1], 0.0);
    assert_eq!(result.solid_endpoints.len(), 1);
    let wall = &result.solid_endpoints[0];
    assert_eq!(wall.cell, 1);
    close(wall.mass_kg, 2.0);
    close(wall.phase.specific_enthalpy_j_kg(), 116.0);
    close(wall.phase.temperature_k(), 516.0);
    close(result.report.solid_enthalpy_change_j, -1.0);
    close(result.report.solid_to_fluid_heat_j, 1.0);
    assert_eq!(result.report.solid_authored_heat_j, 0.0);
    assert_eq!(result.report.conduction.boundary_faces, 0);
    close(
        volume.report.material_enthalpy_change_j + result.report.solid_enthalpy_change_j,
        authored[0] + native.resistance.resistance_heat_j,
    );
    assert!(native.coupling.coupling_iterations > 2);
    assert_eq!(result.work.heat_trials, native.coupling.coupling_iterations);
    assert!(result.work.residual_evaluations > result.report.conduction.solve.residual_evaluations);
    assert!(result.report.conduction.absolute_enthalpy_defect_j <= 1e-7);

    // Cancellation after at least one solved heat trial must not warm-start
    // physical wall history. Scratch reuse then reproduces the original result.
    let gate = CancelGate::new_clock_free();
    let mut queries = 0;
    let cancelled = field
        .advance_isotropic_implicit_conjugate_phase_resistance_thermal_volume_momentum_with_baths(
            &domain,
            faces,
            None,
            &conductivity,
            ThermalVolumeBoundaries::default(),
            &solids,
            &contacts,
            &law,
            data,
            limits,
            solve,
            &mut workspace,
            &gate,
            |_, phase| {
                queries += 1;
                if queries == 2 {
                    gate.request();
                }
                source(&card, &law, phase)
            },
        );
    assert!(matches!(cancelled, Err(ChtError::Cancelled)));
    assert_eq!(queries, 2);
    assert_eq!(
        run(&mut workspace, &CancelGate::new_clock_free()).unwrap(),
        result
    );
    close(solids[0].phase.specific_enthalpy_j_kg(), 116.5);
    assert_eq!(mass, [2.0, 0.0]);
    assert_eq!(momentum, [[0.0; 3]; 2]);
    assert_eq!(field, FlowField::quiescent(&domain));
}

