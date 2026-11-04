//! G1/G3/G4: phase/source feedback, physical heat and actual occupied volume.

use super::*;
use crate::conjugate::{
    CoexistenceVolumeControls, CoexistenceVolumeInput, DensityForceControls,
    DensityMomentumReconstructionControls, DensityProjectionControls, MaterialAdvectionControls,
    MaterialDensityMomentumControls, MaterialMomentumControls, PiecewiseThermalVolumeControls,
    ThermalVolumeControls, ThermalVolumeMomentumControls,
};
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
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

fn temperature() -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(
        QuantityKind::AbsoluteTemperature,
        ValueForm::Static,
    ))
}

fn concentration() -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(
        QuantityKind::MassConcentration,
        ValueForm::Static,
    ))
}

fn fixture(densities: [f64; 2]) -> (MaterialCard, EquilibriumEnthalpyPhaseCurve) {
    let mut claims = ClaimSet::new();
    for (name, quantity, value) in [
        (
            DENSITY_PROPERTY,
            QuantitySpec::dimensional(Density::DIMS),
            0.0,
        ),
        (
            THERMAL_CONDUCTIVITY_PROPERTY,
            QuantitySpec::semantic(SemanticType::new(
                QuantityKind::ThermalConductivity,
                ValueForm::Static,
            )),
            1.0,
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
                key: PropertyKey::with_quantity(name, quantity),
                // Manufactured density response on an explicitly supplied
                // typed mass-concentration coordinate. The 0..8 identity
                // interpolation is exact; no approximate mass-fraction table
                // stands in for the harmonic unequal-density mixture law.
                // This is a test source, not measured material evidence.
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
                        dims: quantity.dims(),
                    }
                },
                validity: ValidityDomain::unconstrained()
                    .with_quantity("T", temperature(), 400.0, 600.0)
                    .with_quantity("liquid", QuantitySpec::dimensional(Dims([0; 6])), 0.0, 1.0)
                    .with_quantity("manufactured_mass_concentration", concentration(), 1.0, 8.0),
                uncertainty: UncertaintyModel::Unstated,
                interpolation: if density {
                    InterpolationPolicy::LinearInside
                } else {
                    InterpolationPolicy::ConstantWithinValidity
                },
                observations: vec![],
                provenance: Provenance {
                    source: "MANUFACTURED coupled Darcy/phase fixture; not measurements".into(),
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
            process: "coupled Darcy heat and density verification".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    let curve = EquilibriumEnthalpyPhaseCurve::try_new(
        card.content_hash(),
        [
            (-100.0, 400.0, 0.0, densities[0]),
            (0.0, 500.0, 0.0, densities[0]),
            (100.0, 500.0, 1.0, densities[1]),
            (200.0, 600.0, 1.0, densities[1]),
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

fn source(
    card: &MaterialCard,
    law: &MushyPhaseLaw<'_>,
    phase: &EquilibriumPhaseState,
) -> Result<ResolvedPhaseTransportState, ChtError> {
    let context = QueryPoint::new()
        .with_quantity("T", temperature(), 550.0)
        .unwrap()
        .with_quantity("liquid", QuantitySpec::dimensional(Dims([0; 6])), 0.9)
        .unwrap()
        .with_quantity(
            "manufactured_mass_concentration",
            concentration(),
            phase.bulk_density_kg_m3(),
        )
        .unwrap()
        .with_quantity(
            "pressure",
            QuantitySpec::semantic(SemanticType::new(QuantityKind::Pressure, ValueForm::Static)),
            100000.0,
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
    .map_err(|error| refuse(error.to_string()))
}

fn controls() -> ResistiveThermalVolumeMomentumControls {
    let projection = DensityProjectionControls {
        duration_s: 0.5,
        relative_tolerance: 1e-13,
        maximum_volume_defect_m3_s: 1e-12,
        maximum_energy_defect_j: 1e-8,
    };
    ResistiveThermalVolumeMomentumControls {
        thermal: ThermalVolumeMomentumControls {
            volume: PiecewiseThermalVolumeControls {
                thermal: ThermalVolumeControls {
                    volume: CoexistenceVolumeControls {
                        transport: MaterialAdvectionControls {
                            projection,
                            max_cells: 16,
                            max_connections: 128,
                            max_iterations: 200,
                            relative_tolerance: 1e-13,
                            maximum_mass_defect_kg: 1e-10,
                            maximum_enthalpy_defect_j: 1e-8,
                        },
                        maximum_cell_volume_defect_m3: 1e-10,
                    },
                    max_volume_iterations: 16,
                    source_relaxation: 1.0,
                },
                max_mixing_intervals: 1024,
            },
            mechanics: MaterialDensityMomentumControls {
                reconstruction: DensityMomentumReconstructionControls {
                    max_cells: 16,
                    maximum_mass_defect_kg: 1e-10,
                    maximum_momentum_defect_kg_m_s: 1e-9,
                    maximum_energy_defect_j: 1e-8,
                },
                forcing: DensityForceControls {
                    max_cells: 16,
                    maximum_impulse_defect_kg_m_s: 1e-9,
                    maximum_energy_defect_j: 1e-8,
                },
                projection,
                maximum_feedback_momentum_defect_kg_m_s: 1e-9,
                maximum_feedback_angular_impulse_defect_kg_m2_s: 1e-9,
                maximum_feedback_energy_defect_j: 1e-8,
            },
            momentum: MaterialMomentumControls {
                max_iterations: 200,
                relative_tolerance: 1e-13,
                maximum_momentum_defect_kg_m_s: 1e-9,
                maximum_energy_defect_j: 1e-8,
            },
        },
        max_coupling_iterations: 64,
        heat_relaxation: 1.0,
        maximum_cell_heat_defect_j: 1e-12,
    }
}

fn data<'a>(
    mass: &'a [f64],
    phase: &'a [Option<EquilibriumPhaseState>],
    heat: &'a [f64],
    momentum: &'a [[f64; 3]],
) -> ThermalVolumeMomentumInput<'a> {
    ThermalVolumeMomentumInput {
        volume: CoexistenceVolumeInput {
            mass_kg: mass,
            phase_states: phase,
            net_enthalpy_increment_j: heat,
            boundary_inflow: [None; 6],
        },
        momentum_kg_m_s: momentum,
        acceleration_m_s2: [0.0; 3],
        force_density_n_m3: &[[0.0; 3]],
        boundary_velocity_m_s: [None; 6],
    }
}

fn outlet() -> [FlowFace; 6] {
    let mut faces = [FlowFace::Wall; 6];
    faces[1] = FlowFace::Free;
    faces
}

fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 5e-9 * expected.abs().max(1.0),
        "{actual:.17e} != {expected:.17e}"
    );
}

#[test]
fn endpoint_phase_density_and_resistance_close_the_same_heat_volume_and_momentum_step() {
    let (card, curve) = fixture([4.0, 2.0]);
    let law = MushyPhaseLaw::try_new_with_immobile_solid(&curve, 4.0, 2.0, 1.0).unwrap();
    let domain = VoxelDomain::new(1, 1, 1, 1.0).unwrap();
    let field = FlowField::quiescent(&domain);
    let phase = [Some(curve.state_at_specific_enthalpy(50.0).unwrap())];
    let mass = [8.0 / 3.0];
    let external = [80.0 / 3.0 - 1.0 / 1215.0];
    let mut queries = 0;
    let result = field
        .advance_resolved_phase_resistance_thermal_volume_momentum(
            &domain,
            outlet(),
            &law,
            data(&mass, &phase, &external, &[[0.0; 3]]),
            controls(),
            &CancelGate::new_clock_free(),
            |_, phase| {
                queries += 1;
                source(&card, &law, phase)
            },
        )
        .unwrap();
    let step = &result.evolution;
    // On the latent interval v(h)=1/4+h/400. At h1=60, rho1=5/2,
    // phi_V=3/4, mu=rho1*nu=5/4 and beta=mu*(1-phi)^2/phi^3=5/27.
    // M0*(h1-h0)=E+H=80/3; Q=(E+H)/(400 dt)=2/15.
    // H=dt*(beta V/2)*Q^2=1/1215. Native inertia remains M0/2=4/3.
    assert!(step.coupling.coupling_iterations > 1);
    assert_eq!(queries, 1 + step.coupling.coupling_iterations);
    let returned = step.thermal.volume.phase_states[0].unwrap();
    close(returned.specific_enthalpy_j_kg(), 60.0);
    close(returned.bulk_density_kg_m3(), 2.5);
    close(step.thermal.volume.transport.mass_kg[0], 2.5);
    close(
        step.thermal
            .volume
            .transport
            .projection
            .field
            .flux_x(1, 0, 0),
        2.0 / 15.0,
    );
    close(
        step.thermal
            .volume
            .transport
            .projection
            .pressure_correction_pa[0],
        149.0 / 405.0,
    );
    close(step.resistance_pa_s_m2[0], 5.0 / 27.0);
    close(step.resistance.resistance_heat_j, 1.0 / 1215.0);
    close(step.caloric_resistance_heat_j[0], 1.0 / 1215.0);
    close(step.resistance.kinetic_energy_after_j, 8.0 / 675.0);
    close(step.resistance.backward_euler_loss_j, 8.0 / 675.0);
    close(step.momentum.momentum_kg_m_s[0][0], 1.0 / 6.0);
    assert!(step.coupling.maximum_cell_heat_defect_j <= controls().maximum_cell_heat_defect_j);
    let source = result.endpoint_sources[0].as_ref().unwrap();
    assert_eq!(*source.phase_state(), returned);
    assert_eq!(source.mobility().law_identity(), law.identity());
    close(source.source_density_kg_m3(), 2.5);
    close(source.mobility().volumes().liquid, 0.75);
    close(
        source.mobile_transport().unwrap().dynamic_viscosity_pa_s(),
        1.25,
    );
    for property in source.material().properties() {
        let point = &property.answer().receipt.query_point;
        assert!(point.contains(&("T".into(), returned.temperature_k())));
        assert!(point.contains(&("liquid".into(), returned.liquid_mass_fraction())));
        assert!(point.contains(&(
            "manufactured_mass_concentration".into(),
            returned.bulk_density_kg_m3()
        )));
        assert!(point.contains(&("pressure".into(), 100000.0)));
    }
    assert_eq!(mass, [8.0 / 3.0]);
    assert_eq!(phase[0].unwrap().specific_enthalpy_j_kg(), 50.0);
}

#[test]
fn wrong_sources_and_immobile_states_refuse_before_an_unchanged_retry() {
    let (card, curve) = fixture([2.0; 2]);
    let law = MushyPhaseLaw::try_new_with_immobile_solid(&curve, 2.0, 2.0, 1.0).unwrap();
    let wrong_law = MushyPhaseLaw::try_new_with_immobile_solid(&curve, 2.0, 2.0, 2.0).unwrap();
    let domain = VoxelDomain::new(1, 1, 1, 1.0).unwrap();
    let field = FlowField::quiescent(&domain);
    let original = field.clone();
    let phase = [Some(curve.state_at_specific_enthalpy(50.0).unwrap())];
    let input = data(&[2.0], &phase, &[0.0], &[[0.0; 3]]);
    let run = |input, gate: &CancelGate| {
        field.advance_resolved_phase_resistance_thermal_volume_momentum(
            &domain,
            outlet(),
            &law,
            input,
            controls(),
            gate,
            |_, phase| source(&card, &law, phase),
        )
    };
    let expected = run(input, &CancelGate::new_clock_free()).unwrap();
    assert!(
        field
            .advance_resolved_phase_resistance_thermal_volume_momentum(
                &domain,
                outlet(),
                &law,
                input,
                controls(),
                &CancelGate::new_clock_free(),
                |_, phase| source(&card, &wrong_law, phase),
            )
            .is_err()
    );
    let wrong_phase = curve.state_at_specific_enthalpy(70.0).unwrap();
    assert!(
        field
            .advance_resolved_phase_resistance_thermal_volume_momentum(
                &domain,
                outlet(),
                &law,
                input,
                controls(),
                &CancelGate::new_clock_free(),
                |_, _| source(&card, &law, &wrong_phase),
            )
            .is_err()
    );
    let solid = [Some(curve.state_at_specific_enthalpy(0.0).unwrap())];
    let solid_input = ThermalVolumeMomentumInput {
        volume: CoexistenceVolumeInput {
            phase_states: &solid,
            ..input.volume
        },
        ..input
    };
    let mut queried = false;
    assert!(
        field
            .advance_resolved_phase_resistance_thermal_volume_momentum(
                &domain,
                outlet(),
                &law,
                solid_input,
                controls(),
                &CancelGate::new_clock_free(),
                |_, phase| {
                    queried = true;
                    source(&card, &law, phase)
                },
            )
            .is_err()
    );
    assert!(!queried);
    let mut incoming_solid = input;
    incoming_solid.volume.boundary_inflow[1] = solid[0];
    assert!(
        field
            .advance_resolved_phase_resistance_thermal_volume_momentum(
                &domain,
                outlet(),
                &law,
                incoming_solid,
                controls(),
                &CancelGate::new_clock_free(),
                |_, phase| {
                    queried = true;
                    source(&card, &law, phase)
                },
            )
            .is_err()
    );
    assert!(!queried);
    let gate = CancelGate::new_clock_free();
    gate.request();
    assert!(matches!(run(input, &gate), Err(ChtError::Cancelled)));
    assert_eq!(field, original);
    let retry = run(input, &CancelGate::new_clock_free()).unwrap();
    assert_eq!(retry, expected);
}

#[test]
fn mobile_inflow_cannot_hide_an_immobile_accepted_caloric_intermediate() {
    let (card, curve) = fixture([2.0; 2]);
    let law = MushyPhaseLaw::try_new_with_immobile_solid(&curve, 2.0, 2.0, 1.0).unwrap();
    let domain = VoxelDomain::new(1, 1, 1, 1.0).unwrap();
    let mut faces = [FlowFace::Wall; 6];
    faces[0] = FlowFace::Fixed;
    faces[1] = FlowFace::Fixed;
    let field =
        FlowField::from_face_velocity(&domain, faces, |axis, _| if axis == 0 { 1.0 } else { 0.0 })
            .unwrap();
    let original = field.clone();
    let phase = [Some(curve.state_at_specific_enthalpy(50.0).unwrap())];
    let mut input = data(&[2.0], &phase, &[-151.0], &[[2.0, 0.0, 0.0]]);
    input.volume.boundary_inflow[0] = Some(curve.state_at_specific_enthalpy(200.0).unwrap());
    input.boundary_velocity_m_s[0] = Some([1.0, 0.0, 0.0]);
    let mut queried_h = Vec::new();
    let error = field
        .advance_resolved_phase_resistance_thermal_volume_momentum(
            &domain,
            faces,
            &law,
            input,
            controls(),
            &CancelGate::new_clock_free(),
            |_, phase| {
                queried_h.push(phase.specific_enthalpy_j_kg());
                source(&card, &law, phase)
            },
        )
        .unwrap_err();
    // The converged mobile endpoint is h=50, beta=2 and physical drag H=1.
    // Before one unit of hot incoming mass mixes, h_star=50+(-151+1)/2=-25.
    // That intermediate is solid and has no admitted mobile momentum owner.
    assert!(
        matches!(error, ChtError::InvalidInput { reason, .. } if reason.contains("exact solid phase"))
    );
    assert!(queried_h.len() > 2);
    assert!(queried_h.into_iter().all(|h| h > 0.0));
    assert_eq!(field, original);
    assert_eq!(phase[0].unwrap().specific_enthalpy_j_kg(), 50.0);
}
