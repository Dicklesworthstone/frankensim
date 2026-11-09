//! Coupled native molecular rate and caloric response with a distinct pore law.

use super::*;
use crate::conjugate::flow::phase_volume::momentum::{
    admit_mobile_phase, resolve_endpoint_rheology_sources, source_resistance,
};
use fs_material::liquid::{DYNAMIC_VISCOSITY_DIMS, DYNAMIC_VISCOSITY_PROPERTY};
use fs_material::phase::transport::resolve_phase_transport_state_with_shear_for_volume_change;
use fs_qty::{Frequency, Temperature};

fn rate_fixture(
    molecular_scale: f64,
    darcy_viscosity: f64,
) -> (MaterialCard, EquilibriumEnthalpyPhaseCurve) {
    let densities = [4.0, 2.0];
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
            DYNAMIC_VISCOSITY_PROPERTY,
            QuantitySpec::dimensional(DYNAMIC_VISCOSITY_DIMS),
            darcy_viscosity,
        ),
        (
            KINEMATIC_VISCOSITY_PROPERTY,
            QuantitySpec::dimensional(KINEMATIC_VISCOSITY_DIMS),
            0.5,
        ),
    ] {
        let conductivity = name == THERMAL_CONDUCTIVITY_PROPERTY;
        let density = name == DENSITY_PROPERTY;
        let molecular = name == KINEMATIC_VISCOSITY_PROPERTY;
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
                } else if molecular {
                    // Independent manufactured molecular nu(T, gamma), with
                    // actual source density supplying its dynamic conversion.
                    let nu = |t: f64, rate: f64| {
                        molecular_scale * (0.05 + 0.002 * (t - 500.0) + 0.3 * rate)
                    };
                    PropertyValue::Surface2D {
                        axes: ["T".into(), "rate".into()],
                        axis_dims: [Temperature::DIMS, Frequency::DIMS],
                        samples: vec![
                            ([500.0, 0.0], nu(500.0, 0.0)),
                            ([520.0, 0.0], nu(520.0, 0.0)),
                            ([500.0, 4.0], nu(500.0, 4.0)),
                            ([520.0, 4.0], nu(520.0, 4.0)),
                        ],
                        triangles: vec![[0, 1, 2], [1, 3, 2]],
                        dims: spec.dims(),
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
                validity: {
                    let validity = ValidityDomain::unconstrained()
                    .with_quantity(
                        "T",
                        quantity(QuantityKind::AbsoluteTemperature),
                        if molecular { 500.0 } else { 400.0 },
                        if molecular { 520.0 } else { 600.0 },
                    )
                    .with_quantity("liquid", QuantitySpec::dimensional(Dims([0; 6])), 0.0, 1.0)
                    .with_quantity(
                        "manufactured_mass_concentration",
                        quantity(QuantityKind::MassConcentration),
                        1.0,
                        8.0,
                    );
                    if molecular {
                        validity.with_quantity("rate", QuantitySpec::dimensional(Frequency::DIMS), 0.0, 4.0)
                    } else {
                        validity
                    }
                },
                uncertainty: UncertaintyModel::Unstated,
                interpolation: if conductivity || density || molecular {
                    InterpolationPolicy::LinearInside
                } else {
                    InterpolationPolicy::ConstantWithinValidity
                },
                observations: vec![],
                provenance: Provenance {
                    source: "MANUFACTURED independent molecular T/rate and pore-response laws; not measurements"
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
            process: "nonisothermal mushy native rheology and independent pore drag".into(),
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
            (100.0, 520.0, 1.0),
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

/// The scalar dynamic pore response is explicitly authored independently of
/// the molecular kinematic T/rate surface. Both retain this chart and material;
/// no macroscopic rate is fabricated or forwarded to the Darcy law.
fn pore_source(
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
            viscosity_basis: LiquidViscosityBasis::Dynamic,
            extra_requirements: &[],
            selection: MaterialPropertySelection::SingleClaimOnly,
        },
    )
    .map_err(|error| refuse(error.to_string()))
}

fn molecular_source(
    card: &MaterialCard,
    law: &MushyPhaseLaw<'_>,
    phase: &EquilibriumPhaseState,
    rate: f64,
) -> Result<(ResolvedPhaseTransportState, f64), ChtError> {
    let context = context()
        .with_quantity(
            "manufactured_mass_concentration",
            quantity(QuantityKind::MassConcentration),
            phase.bulk_density_kg_m3(),
        )
        .unwrap()
        .with_quantity("rate", QuantitySpec::dimensional(Frequency::DIMS), rate)
        .unwrap();
    let source = resolve_phase_transport_state_with_shear_for_volume_change(
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
        "rate",
        rate,
    )
    .map_err(|error| refuse(error.to_string()))?;
    Ok((
        source,
        0.08 + 0.001 * (phase.temperature_k() - 500.0) + 0.2 * rate,
    ))
}

#[test]
fn mobile_molecular_evidence_does_not_require_an_unused_darcy_product() {
    let (card, curve) = rate_fixture(1e200, 1e-200);
    let law = MushyPhaseLaw::try_new_with_immobile_solid(&curve, 4.0, 2.0, 1e-200).unwrap();
    let phase = curve.state_at_specific_enthalpy(50.0).unwrap();
    let domain = VoxelDomain::new(1, 1, 1, 1.0).unwrap();
    let mut calls = 0;
    let admitted = resolve_endpoint_rheology_sources(
        &domain,
        &curve,
        Some(&law),
        |_| Ok(phase),
        &[0.5],
        &mut |_, phase, rate| molecular_source(&card, &law, phase, rate),
        &mut calls,
        &CancelGate::new_clock_free(),
    )
    .unwrap();
    assert_eq!(calls, 1);
    assert!(admitted.shear_pa_s[0].is_finite() && admitted.shear_pa_s[0] > 1e198);
    let mobile = admit_mobile_phase(&law, &phase).unwrap();
    assert!(!(admitted.shear_pa_s[0] * mobile.inverse_permeability_per_m2).is_finite());
    // The actual Darcy consumer still refuses a nonrepresentable beta.
    assert!(source_resistance(mobile, admitted.receipts[0].as_ref().unwrap()).is_err());
    // The existing thermal-only shared-source entry still consumes its real
    // Darcy viscosity and must retain the original beta representability gate.
    let mut thermal_work = ThermalViscosityFeedbackReport::default();
    assert!(crate::conjugate::flow::phase_volume::conduction::implicit::viscous::feedback::resolve_thermal_sources(
        &domain, &curve, Some(&law), |_| Ok(phase),
        &mut |_, phase| molecular_source(&card, &law, phase, 0.5),
        &mut thermal_work, &CancelGate::new_clock_free(),
    ).is_err());
    let pore = pore_source(&card, &law, &phase).unwrap();
    let beta = source_resistance(mobile, &pore).unwrap();
    assert!(beta.is_finite() && beta > 0.0);
    close(beta, 1e-200 * mobile.inverse_permeability_per_m2);
}

#[test]
fn joint_thermal_rate_mushy_flow_retains_original_mass_and_independent_darcy_heat() {
    let (card, curve) = rate_fixture(1.0, 0.75);
    let law = MushyPhaseLaw::try_new_with_immobile_solid(&curve, 4.0, 2.0, 1.0).unwrap();
    let context = context();
    let conductivity = prepare(&card, &curve, &context, [40.0, 90.0]);
    let wall_conductivity = prepare(&card, &curve, &context, [110.0, 130.0]);
    let mut domain = VoxelDomain::new(1, 2, 1, 1.0).unwrap();
    domain.set(0, 1, 0, Voxel::Solid(0));
    let field = FlowField::quiescent(&domain);
    let initial = curve.state_at_specific_enthalpy(50.0).unwrap();
    let phases = [Some(initial), None];
    let mass = [initial.bulk_density_kg_m3(), 0.0];
    let momentum = [[0.1 * mass[0], 0.0, 0.0], [0.0; 3]];
    let heat = [mass[0] * 10.0, 0.0];
    let force = [[0.0; 3]; 2];
    let mut input = ThermalVolumeMomentumInput {
        volume: CoexistenceVolumeInput {
            mass_kg: &mass,
            phase_states: &phases,
            net_enthalpy_increment_j: &heat,
            boundary_inflow: [None; 6],
        },
        momentum_kg_m_s: &momentum,
        acceleration_m_s2: [0.0; 3],
        force_density_n_m3: &force,
        boundary_velocity_m_s: [None; 6],
    };
    let mut faces = [FlowFace::Wall; 6];
    let mut boundaries = [FvBoundary::Symmetry; 6];
    for side in 0..2 {
        faces[side] = FlowFace::Free;
        boundaries[side] = FvBoundary::Outlet;
        boundaries[2 + side] = FvBoundary::wall();
        input.volume.boundary_inflow[side] = Some(initial);
        input.boundary_velocity_m_s[side] = Some([0.0; 3]);
    }
    let px = [Some(0.0), Some(0.0), None, None];
    let py = [None; 3];
    let pz = [None; 4];
    let reservoir = Some((
        faces,
        StaggeredBoundaryPressure {
            x_pa: &px,
            y_pa: &py,
            z_pa: &pz,
        },
    ));
    let solids = [ThermalVolumeSolidCell {
        cell: 1,
        conductivity: &wall_conductivity,
        mass_kg: 2.0,
        phase: curve.state_at_specific_enthalpy(116.5).unwrap(),
        net_enthalpy_increment_j: 0.0,
    }];
    let contact_resistance = 16.0 / 15.0;
    let contacts = [ThermalVolumeContact::Resistance(FaceContactResistance {
        cells: [0, 1],
        resistance_m2_k_w: contact_resistance,
    })];
    let limits = limits(2, 1);
    let solve = ImplicitThermalVolumeSolveControls {
        cell_tolerance_j: 1e-12,
        ..ImplicitThermalVolumeSolveControls::default()
    };
    let mut workspace = ImplicitThermalVolumeWorkspace::new(2).unwrap();
    let molecular = MaterialViscosityControls {
        max_sweeps: 40_000,
        ..molecular_controls()
    };
    let feedback = ThermalViscosityFeedbackControls {
        maximum_trials: 1024,
        relaxation: 1.0,
        maximum_heat_defect_j: 1e-11,
    };
    let run = |workspace: &mut ImplicitThermalVolumeWorkspace,
               solve,
               molecular,
               feedback,
               gate: &CancelGate,
               mode: u8| {
        let mut molecular_calls = 0;
        let mut darcy_calls = 0;
        let result = field
            .advance_isotropic_implicit_conjugate_thermorheological_phase_resistance_thermal_volume_momentum_with_baths(
                &domain, reservoir, &conductivity, ThermalVolumeBoundaries::default(),
                &solids, &contacts, &law, input, boundaries, limits, molecular,
                solve, feedback, workspace, gate,
                |_, phase, rate| {
                    molecular_calls += 1;
                    let heated = phase.specific_enthalpy_j_kg() > 51.0;
                    if mode == 1 && heated { gate.request(); }
                    let source_phase = if mode == 2 && heated { &initial } else { phase };
                    molecular_source(&card, &law, source_phase, rate)
                },
                |_, phase| {
                    darcy_calls += 1;
                    if mode == 3 && phase.specific_enthalpy_j_kg() > 51.0 {
                        return Err(refuse("manufactured late pore-source refusal"));
                    }
                    pore_source(&card, &law, phase)
                },
            );
        (result, molecular_calls, darcy_calls)
    };
    let (result, molecular_calls, darcy_calls) = run(
        &mut workspace,
        solve,
        molecular,
        feedback,
        &CancelGate::new_clock_free(),
        0,
    );
    let result = result.unwrap();
    let thermal = &result.evolution;
    let coupled = &thermal.evolution;
    let combined = &coupled.evolution;
    let native = &combined.evolution.evolution;
    let volume = &native.thermal.volume;
    let accepted = combined.viscosity_sources[0].as_ref().unwrap();
    let caloric = accepted.phase_state();
    let rate = result.endpoint_shear_rate_per_s[0];
    assert!(rate > 0.0);
    assert!(caloric.temperature_k() > initial.temperature_k());
    assert!(caloric.liquid_mass_fraction() > initial.liquid_mass_fraction());
    assert!(caloric.liquid_mass_fraction() < 1.0);
    assert!(caloric.bulk_density_kg_m3() < mass[0]);
    let nu = 0.05 + 0.002 * (caloric.temperature_k() - 500.0) + 0.3 * rate;
    let mu = accepted
        .mobile_transport()
        .unwrap()
        .dynamic_viscosity_pa_s();
    close(mu, caloric.bulk_density_kg_m3() * nu);
    close(
        thermal.bulk_viscosity_pa_s[0],
        0.08 + 0.001 * (caloric.temperature_k() - 500.0) + 0.2 * rate,
    );
    let pore = combined.evolution.endpoint_sources[0].as_ref().unwrap();
    close(
        pore.mobile_transport().unwrap().dynamic_viscosity_pa_s(),
        0.75,
    );
    assert!((mu - 0.75).abs() > 1e-3);
    assert_eq!(
        thermal.feedback.source_evaluations,
        molecular_calls + darcy_calls
    );
    assert!(result.rheology_source_updates > thermal.feedback.molecular_solves);
    assert!(thermal.feedback.trials > native.coupling.coupling_iterations);
    assert!(combined.viscosity.relative_residual <= molecular.relative_tolerance);
    assert!(combined.viscosity.molecular_heat_j > 0.0);
    assert!(native.resistance.resistance_heat_j > 0.0);

    // Close the same nonlinear native rows independently at the accepted exact
    // thermal phase. Original mass is authoritative despite its lower hot rho.
    let prepare_reference = |mass: &[f64]| {
        field.prepare_material_rheology_with_source(
            &domain,
            faces,
            Some(0),
            mass,
            &momentum,
            boundaries,
            limits.resistance.thermal.mechanics.reconstruction,
            0.25,
            molecular,
            &CancelGate::new_clock_free(),
            |rates| {
                let (source, bulk) = molecular_source(&card, &law, caloric, rates[0])?;
                Ok((
                    vec![
                        source.mobile_transport().unwrap().dynamic_viscosity_pa_s(),
                        0.0,
                    ],
                    vec![bulk, 0.0],
                ))
            },
        )
    };
    let reference = prepare_reference(&mass).unwrap();
    assert_eq!(combined.viscosity, reference.viscosity.report);
    assert_eq!(
        combined.viscosity_reconstruction,
        reference.viscosity.reconstruction
    );
    assert_eq!(
        result.endpoint_shear_rate_per_s,
        reference.endpoint_shear_rate_per_s
    );
    let wrong_mass = [caloric.bulk_density_kg_m3(), 0.0];
    let wrong = prepare_reference(&wrong_mass).unwrap();
    assert!(
        (wrong.viscosity.report.molecular_heat_j - combined.viscosity.molecular_heat_j).abs()
            > 1e-7
    );

    let used_heat = heat[0]
        + thermal.used_molecular_heat_j[0]
        + native.caloric_resistance_heat_j[0]
        + coupled.conductive_enthalpy_increment_j[0];
    close(caloric.specific_enthalpy_j_kg(), 50.0 + used_heat / mass[0]);
    close(volume.report.requested_enthalpy_increment_j, used_heat);
    close(volume.volume_source_m3_s[0], used_heat / (400.0 * 0.25));
    let endpoint = volume.phase_states[0].unwrap();
    close(
        volume.transport.mass_kg[0] / endpoint.bulk_density_kg_m3(),
        1.0,
    );
    let physical_heat =
        heat[0] + combined.viscosity.molecular_heat_j + native.resistance.resistance_heat_j;
    close(
        volume.transport.mass_kg[0] * endpoint.specific_enthalpy_j_kg() - mass[0] * 50.0
            + coupled.report.solid_enthalpy_change_j,
        physical_heat + volume.transport.report.carried_enthalpy_j,
    );
    let wall = &coupled.solid_endpoints[0];
    close(
        wall.mass_kg * (wall.phase.specific_enthalpy_j_kg() - 116.5),
        -coupled.conductive_enthalpy_increment_j[0],
    );
    close(
        volume.transport.mass_kg[0] - mass[0],
        volume.transport.report.incoming_mass_kg - volume.transport.report.outgoing_mass_kg,
    );
    assert!(thermal.feedback.absolute_heat_coupling_defect_j <= feedback.maximum_heat_defect_j);
    assert_eq!(
        coupled.report.conduction.absolute_enthalpy_defect_j,
        combined.absolute_enthalpy_defect_j
    );
    assert!(combined.absolute_energy_defect_j <= molecular.maximum_energy_defect_j);
    assert_eq!(result.endpoint_shear_rate_per_s[1], 0.0);
    assert!(combined.viscosity_sources[1].is_none());

    let (cancelled, _, _) = run(
        &mut workspace,
        solve,
        molecular,
        feedback,
        &CancelGate::new_clock_free(),
        1,
    );
    assert!(matches!(cancelled, Err(ChtError::Cancelled)));
    for mode in [2, 3] {
        assert!(
            run(
                &mut workspace,
                solve,
                molecular,
                feedback,
                &CancelGate::new_clock_free(),
                mode
            )
            .0
            .is_err()
        );
    }
    let short_heat = ImplicitThermalVolumeSolveControls {
        maximum_residual_evaluations: thermal.feedback.heat_residual_evaluations - 1,
        ..solve
    };
    assert!(
        run(
            &mut workspace,
            short_heat,
            molecular,
            feedback,
            &CancelGate::new_clock_free(),
            0
        )
        .0
        .is_err()
    );
    let short_molecular = MaterialViscosityControls {
        max_sweeps: thermal.feedback.molecular_sweeps - 1,
        ..molecular
    };
    assert!(
        run(
            &mut workspace,
            solve,
            short_molecular,
            feedback,
            &CancelGate::new_clock_free(),
            0
        )
        .0
        .is_err()
    );
    let short_trials = ThermalViscosityFeedbackControls {
        maximum_trials: thermal.feedback.trials - 1,
        ..feedback
    };
    assert!(
        run(
            &mut workspace,
            solve,
            molecular,
            short_trials,
            &CancelGate::new_clock_free(),
            0
        )
        .0
        .is_err()
    );
    let (retry, m_calls, d_calls) = run(
        &mut workspace,
        solve,
        molecular,
        feedback,
        &CancelGate::new_clock_free(),
        0,
    );
    assert_eq!(retry.unwrap(), result);
    assert_eq!((m_calls, d_calls), (molecular_calls, darcy_calls));
}
