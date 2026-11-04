//! Original phase/source provenance and actual-mass stationary resistance.

use super::*;
use crate::conjugate::flow::{
    DensityForceControls, DensityMomentumReconstructionControls, DensityProjectionControls,
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
    resolve_phase_transport_state_with_density_departure,
};
use fs_material::phase::{EnthalpyPhaseKnot, EquilibriumEnthalpyPhaseCurve};
use fs_material::state_point::{
    DENSITY_PROPERTY, MaterialPropertySelection, SPECIFIC_HEAT_CAPACITY_DIMS,
    SPECIFIC_HEAT_CAPACITY_PROPERTY, THERMAL_CONDUCTIVITY_PROPERTY,
};
use fs_material::thermal::source::ThermalMaterialPoint;
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_qty::{Density, Dims, QuantitySpec};

fn absolute_temperature() -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(
        QuantityKind::AbsoluteTemperature,
        ValueForm::Static,
    ))
}

fn fixture(
    densities: [f64; 2],
    source_density_scale: f64,
    kinematic_viscosity: f64,
) -> (MaterialCard, EquilibriumEnthalpyPhaseCurve) {
    let mixture_density = (0.5 / densities[0] + 0.5 / densities[1]).recip();
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
            3.0,
        ),
        (
            KINEMATIC_VISCOSITY_PROPERTY,
            QuantitySpec::dimensional(KINEMATIC_VISCOSITY_DIMS),
            kinematic_viscosity,
        ),
    ] {
        let density = name == DENSITY_PROPERTY;
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::with_quantity(name, quantity),
                value: if density {
                    PropertyValue::Curve {
                        abscissa: "liquid".into(),
                        abscissa_dims: Dims([0; 6]),
                        knots: vec![
                            (0.0, densities[0] * source_density_scale),
                            (0.5, mixture_density * source_density_scale),
                            (1.0, densities[1] * source_density_scale),
                        ],
                        dims: Density::DIMS,
                    }
                } else {
                    PropertyValue::Scalar {
                        value,
                        dims: quantity.dims(),
                    }
                },
                validity: ValidityDomain::unconstrained()
                    .with_quantity("T", absolute_temperature(), 400.0, 600.0)
                    .with_quantity("liquid", QuantitySpec::dimensional(Dims([0; 6])), 0.0, 1.0),
                uncertainty: UncertaintyModel::Unstated,
                interpolation: if density {
                    InterpolationPolicy::LinearInside
                } else {
                    InterpolationPolicy::ConstantWithinValidity
                },
                observations: vec![],
                provenance: Provenance {
                    source: "MANUFACTURED phase resistance fixture; not measurements".into(),
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
            process: "actual-density resistance verification".into(),
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
        .map(
            |(h, temperature_k, liquid_mass_fraction, bulk_density_kg_m3)| EnthalpyPhaseKnot {
                specific_enthalpy_j_kg: h,
                temperature_k,
                liquid_mass_fraction,
                bulk_density_kg_m3,
            },
        )
        .collect(),
    )
    .unwrap();
    (card, curve)
}

fn resolve(
    card: &MaterialCard,
    law: &MushyPhaseLaw<'_>,
    h: f64,
    departure: Option<f64>,
) -> Result<ResolvedPhaseTransportState, ChtError> {
    let context = QueryPoint::new()
        .with_quantity("T", absolute_temperature(), 550.0)
        .unwrap()
        .with_quantity("liquid", QuantitySpec::dimensional(Dims([0; 6])), 0.9)
        .unwrap()
        .with_quantity(
            "pressure",
            QuantitySpec::semantic(SemanticType::new(QuantityKind::Pressure, ValueForm::Static)),
            100000.0,
        )
        .unwrap();
    let query = PhaseTransportSourceQuery {
        source: ThermalMaterialPoint {
            material: card,
            context: &context,
            temperature_axis: "T",
            liquid_fraction_axis: Some("liquid"),
        },
        viscosity_basis: LiquidViscosityBasis::Kinematic,
        extra_requirements: &[],
        selection: MaterialPropertySelection::SingleClaimOnly,
    };
    match departure {
        Some(limit) => resolve_phase_transport_state_with_density_departure(law, h, query, limit),
        None => resolve_phase_transport_state_for_volume_change(law, h, query),
    }
    .map_err(|error| ChtError::InvalidInput {
        field: "flow.phase_resistance_test_source",
        reason: error.to_string(),
    })
}

fn controls() -> PhaseResistanceMomentumControls {
    PhaseResistanceMomentumControls {
        mechanics: MaterialDensityMomentumControls {
            reconstruction: DensityMomentumReconstructionControls {
                max_cells: 100,
                maximum_mass_defect_kg: 1e-10,
                maximum_momentum_defect_kg_m_s: 1e-9,
                maximum_energy_defect_j: 1e-7,
            },
            forcing: DensityForceControls {
                max_cells: 100,
                maximum_impulse_defect_kg_m_s: 1e-9,
                maximum_energy_defect_j: 1e-7,
            },
            projection: DensityProjectionControls {
                duration_s: 0.5,
                relative_tolerance: 1e-12,
                maximum_volume_defect_m3_s: 1e-10,
                maximum_energy_defect_j: 1e-7,
            },
            maximum_feedback_momentum_defect_kg_m_s: 1e-8,
            maximum_feedback_angular_impulse_defect_kg_m2_s: 1e-8,
            maximum_feedback_energy_defect_j: 1e-7,
        },
        maximum_cell_volume_defect_m3: 1e-10,
    }
}

fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() < 4e-10 * expected.abs().max(1.0),
        "{actual:e} != {expected:e}"
    );
}

fn open_x() -> [FlowFace; 6] {
    [
        FlowFace::Free,
        FlowFace::Free,
        FlowFace::Wall,
        FlowFace::Wall,
        FlowFace::Wall,
        FlowFace::Wall,
    ]
}

fn input<'a>(
    mass: &'a [f64],
    phases: &'a [Option<EquilibriumPhaseState>],
    momentum: &'a [[f64; 3]],
    force: &'a [[f64; 3]],
    targets: &'a [f64],
) -> PhaseResistanceMomentumInput<'a> {
    PhaseResistanceMomentumInput {
        mechanics: MaterialDensityMomentumInput {
            mass_kg: mass,
            momentum_kg_m_s: momentum,
            acceleration_m_s2: [0.0; 3],
            force_density_n_m3: force,
            target_outflow_m3_s: targets,
        },
        phase_states: phases,
    }
}

#[allow(clippy::too_many_arguments)]
fn run(
    field: &FlowField,
    domain: &VoxelDomain,
    law: &MushyPhaseLaw<'_>,
    input: PhaseResistanceMomentumInput<'_>,
    controls: PhaseResistanceMomentumControls,
    gate: &CancelGate,
    resolve: impl FnMut(usize, &EquilibriumPhaseState) -> Result<ResolvedPhaseTransportState, ChtError>,
) -> Result<ResolvedPhaseResistanceMomentumStep, ChtError> {
    field.force_project_material_momentum_with_phase_resistance(
        domain,
        open_x(),
        law,
        input,
        controls,
        gate,
        resolve,
    )
}

#[test]
fn unequal_phase_density_resolves_volume_fraction_viscosity_and_native_impulse_once() {
    let (card, curve) = fixture([4.0, 2.0], 1.0, 0.5);
    let law = MushyPhaseLaw::try_new_with_immobile_solid(&curve, 4.0, 2.0, 1.0).unwrap();
    let domain = VoxelDomain::graded([vec![1.0, 2.0], vec![1.0], vec![1.0]]).unwrap();
    let field = FlowField::quiescent(&domain);
    let phases = [
        Some(curve.state_at_specific_enthalpy(50.0).unwrap()),
        Some(curve.state_at_specific_enthalpy(150.0).unwrap()),
    ];
    let mass = [8.0 / 3.0, 4.0];
    let zero = [[0.0; 3]; 2];
    let data = input(&mass, &phases, &zero, &zero, &[0.0; 2]);
    let mut calls = [0; 2];
    let result = field
        .force_project_material_momentum_with_phase_resistance_and_reservoir_pressure(
            &domain,
            open_x(),
            [Some(4.0), Some(0.0), None, None, None, None],
            &law,
            data,
            controls(),
            &CancelGate::new_clock_free(),
            |cell, phase| {
                calls[cell] += 1;
                resolve(&card, &law, phase.specific_enthalpy_j_kg(), None)
            },
        )
        .unwrap();
    assert_eq!(calls, [1, 1]);
    let mixed = result.original_sources[0].as_ref().unwrap();
    close(mixed.phase_state().liquid_mass_fraction(), 0.5);
    close(mixed.mobility().volumes().liquid, 2.0 / 3.0);
    close(mixed.source_density_kg_m3(), 8.0 / 3.0);
    close(
        mixed.mobile_transport().unwrap().dynamic_viscosity_pa_s(),
        4.0 / 3.0,
    );
    close(result.resistance_pa_s_m2[0], 0.5);
    assert_eq!(result.resistance_pa_s_m2[1], 0.0);
    for (cell, source) in result.original_sources.iter().enumerate() {
        let source = source.as_ref().unwrap();
        let phase = phases[cell].unwrap();
        assert_eq!(*source.phase_state(), phase);
        assert_eq!(source.mobility().law_identity(), law.identity());
        assert_eq!(source.material().card_identity(), card.content_hash());
        for property in source.material().properties() {
            let point = &property.answer().receipt.query_point;
            assert!(point.contains(&("T".into(), phase.temperature_k())));
            assert!(point.contains(&("liquid".into(), phase.liquid_mass_fraction())));
            assert!(point.contains(&("pressure".into(), 100000.0)));
        }
    }
    // The series column has actual M=20/3 and integrated B=1/2.
    // (M + dt B)u = dt delta_p gives u=24/83, not a reference-density speed.
    let speed = 24.0 / 83.0;
    for cell in 0..2 {
        close(
            result.evolution.momentum_kg_m_s[cell][0],
            mass[cell] * speed,
        );
        close(
            result.evolution.density_kg_m3[cell],
            mass[cell] / domain.volume(cell),
        );
    }
    let native = &result.evolution.evolution.report;
    close(native.pressure_impulse_kg_m_s[0], 2.0);
    close(native.resistance_impulse_kg_m_s[0], -0.25 * speed);
    close(native.cell_resistance_heat_j[0], 0.25 * speed * speed);
    assert_eq!(native.cell_resistance_heat_j[1], 0.0);
    close(native.resistance_heat_j, 0.25 * speed * speed);
    close(result.evolution.report.composed_energy_defect_j, 0.0);
    assert_eq!(phases[0].unwrap().specific_enthalpy_j_kg(), 50.0);
}

#[test]
fn liquid_zero_drag_preserves_actual_mass_and_rebuilds_stale_face_guesses() {
    let (card, curve) = fixture([4.0, 2.0], 1.0, 0.5);
    let law = MushyPhaseLaw::try_new_with_immobile_solid(&curve, 4.0, 2.0, 1.0).unwrap();
    let domain = VoxelDomain::new(1, 1, 1, 1.0).unwrap();
    let field = FlowField::quiescent(&domain);
    let mass = [2.125];
    let momentum = [[1.0625, 0.0, 0.0]];
    let phases = [Some(curve.state_at_specific_enthalpy(150.0).unwrap())];
    let mut data = input(&mass, &phases, &momentum, &[[0.0; 3]], &[0.0]);
    data.mechanics.acceleration_m_s2 = [2.0, 0.0, 0.0];
    let mut admitted = controls();
    admitted.maximum_cell_volume_defect_m3 = 0.0625;
    let result = run(
        &field,
        &domain,
        &law,
        data,
        admitted,
        &CancelGate::new_clock_free(),
        |_, phase| resolve(&card, &law, phase.specific_enthalpy_j_kg(), None),
    )
    .unwrap();
    let original = field
        .force_project_material_momentum(
            &domain,
            open_x(),
            data.mechanics,
            admitted.mechanics,
            &CancelGate::new_clock_free(),
        )
        .unwrap();
    assert_eq!(result.resistance_pa_s_m2, [0.0]);
    assert_eq!(result.evolution.density_kg_m3, [2.125]);
    assert_eq!(
        result.original_sources[0]
            .as_ref()
            .unwrap()
            .source_density_kg_m3(),
        2.0
    );
    assert_eq!(result.evolution.momentum_kg_m_s, original.momentum_kg_m_s);
    close(result.evolution.momentum_kg_m_s[0][0], 2.125 * 1.5);
    assert_eq!(result.evolution.evolution.report.resistance_heat_j, 0.0);
    let stale =
        FlowField::from_face_velocity(
            &domain,
            open_x(),
            |axis, _| if axis == 0 { 17.0 } else { 0.0 },
        )
        .unwrap();
    let rebuilt = run(
        &stale,
        &domain,
        &law,
        data,
        admitted,
        &CancelGate::new_clock_free(),
        |_, phase| resolve(&card, &law, phase.specific_enthalpy_j_kg(), None),
    )
    .unwrap();
    assert_eq!(rebuilt.evolution, result.evolution);
    assert_eq!(rebuilt.original_sources, result.original_sources);
}

#[test]
fn phase_source_and_support_refusals_preserve_an_unchanged_retry() {
    let (card, curve) = fixture([4.0, 2.0], 1.0, 0.5);
    let law = MushyPhaseLaw::try_new_with_immobile_solid(&curve, 4.0, 2.0, 1.0).unwrap();
    let other_law = MushyPhaseLaw::try_new_with_immobile_solid(&curve, 4.0, 2.0, 2.0).unwrap();
    let domain = VoxelDomain::new(1, 1, 1, 1.0).unwrap();
    let field = FlowField::quiescent(&domain);
    let original = field.clone();
    let mass = [8.0 / 3.0];
    let momentum = [[8.0 / 3.0, 0.0, 0.0]];
    let phases = [Some(curve.state_at_specific_enthalpy(50.0).unwrap())];
    let data = input(&mass, &phases, &momentum, &[[0.0; 3]], &[0.0]);
    let expected = run(
        &field,
        &domain,
        &law,
        data,
        controls(),
        &CancelGate::new_clock_free(),
        |_, phase| resolve(&card, &law, phase.specific_enthalpy_j_kg(), None),
    )
    .unwrap();
    assert!(
        run(
            &field,
            &domain,
            &law,
            data,
            controls(),
            &CancelGate::new_clock_free(),
            |_, phase| resolve(&card, &other_law, phase.specific_enthalpy_j_kg(), None)
        )
        .is_err()
    );
    assert!(
        run(
            &field,
            &domain,
            &law,
            data,
            controls(),
            &CancelGate::new_clock_free(),
            |_, _| resolve(&card, &law, 150.0, None)
        )
        .is_err()
    );

    let solid = [Some(curve.state_at_specific_enthalpy(-50.0).unwrap())];
    let solid_data = input(&[4.0], &solid, &[[0.0; 3]], &[[0.0; 3]], &[0.0]);
    let mut queried_solid = false;
    assert!(
        run(
            &field,
            &domain,
            &law,
            solid_data,
            controls(),
            &CancelGate::new_clock_free(),
            |_, phase| {
                queried_solid = true;
                resolve(&card, &law, phase.specific_enthalpy_j_kg(), None)
            }
        )
        .is_err()
    );
    assert!(!queried_solid);
    let mut bad_volume = data;
    bad_volume.mechanics.mass_kg = &[3.0];
    assert!(
        run(
            &field,
            &domain,
            &law,
            bad_volume,
            controls(),
            &CancelGate::new_clock_free(),
            |_, phase| resolve(&card, &law, phase.specific_enthalpy_j_kg(), None)
        )
        .is_err()
    );

    // A valid explicit reference-density source is not actual-density proof.
    let (departure_card, departure_curve) = fixture([2.0; 2], 0.95, 0.5);
    let departure_law =
        MushyPhaseLaw::try_new_with_immobile_solid(&departure_curve, 2.0, 2.0, 1.0).unwrap();
    let departure_phase = [Some(
        departure_curve.state_at_specific_enthalpy(50.0).unwrap(),
    )];
    let departure_data = input(&[2.0], &departure_phase, &[[0.0; 3]], &[[0.0; 3]], &[0.0]);
    let snapshot = resolve(&departure_card, &departure_law, 50.0, Some(0.1)).unwrap();
    close(snapshot.source_density_kg_m3(), 1.9);
    assert!(
        run(
            &field,
            &domain,
            &departure_law,
            departure_data,
            controls(),
            &CancelGate::new_clock_free(),
            |_, _| Ok(snapshot.clone())
        )
        .is_err()
    );

    for (nu, permeability) in [(1e-30, 1e300), (1e20, 1e-300)] {
        let (card, curve) = fixture([4.0, 2.0], 1.0, nu);
        let law =
            MushyPhaseLaw::try_new_with_immobile_solid(&curve, 4.0, 2.0, permeability).unwrap();
        let phases = [Some(curve.state_at_specific_enthalpy(50.0).unwrap())];
        let data = input(&mass, &phases, &momentum, &[[0.0; 3]], &[0.0]);
        // Both operands resolve normally; only their beta product is outside
        // the representable positive range, so this reaches the adapter gate.
        let snapshot = resolve(&card, &law, 50.0, None).unwrap();
        let inverse = law
            .at_specific_enthalpy(50.0)
            .unwrap()
            .inverse_permeability_per_m2;
        assert!(inverse.is_finite() && inverse > 0.0);
        let beta = snapshot
            .mobile_transport()
            .unwrap()
            .dynamic_viscosity_pa_s()
            * inverse;
        if nu < 1.0 {
            assert_eq!(beta, 0.0);
        } else {
            assert!(beta.is_infinite());
        }
        assert!(
            run(
                &field,
                &domain,
                &law,
                data,
                controls(),
                &CancelGate::new_clock_free(),
                |_, _| Ok(snapshot.clone())
            )
            .is_err()
        );
    }
    let gate = CancelGate::new_clock_free();
    gate.request();
    assert!(matches!(
        run(
            &field,
            &domain,
            &law,
            data,
            controls(),
            &gate,
            |_, phase| resolve(&card, &law, phase.specific_enthalpy_j_kg(), None)
        ),
        Err(ChtError::Cancelled)
    ));
    assert_eq!(field, original);
    assert_eq!(mass, [8.0 / 3.0]);
    assert_eq!(momentum, [[8.0 / 3.0, 0.0, 0.0]]);
    let retry = run(
        &field,
        &domain,
        &law,
        data,
        controls(),
        &CancelGate::new_clock_free(),
        |_, phase| resolve(&card, &law, phase.specific_enthalpy_j_kg(), None),
    )
    .unwrap();
    assert_eq!(retry.evolution, expected.evolution);
    assert_eq!(retry.original_sources, expected.original_sources);
}
