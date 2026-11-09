//! Original nonlinear conduction, finite wall heat and native expansion work.

use super::super::implicit::{
    ImplicitHeatWorkThermalVolumeControls, ImplicitHeatWorkThermalVolumeStep,
};
use super::*;
use crate::conjugate::domain::Voxel;
use crate::conjugate::{
    ConductiveThermalVolumeMomentumControls, FaceContactResistance,
    ImplicitThermalVolumeSolveControls, ImplicitThermalVolumeWorkspace, ThermalVolumeBathFace,
    ThermalVolumeBoundaries, ThermalVolumeContact, ThermalVolumeSolidCell,
};
use crate::d3q19::Face3;
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSelection, ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim,
    PropertyKey, PropertyValue, Provenance, QueryPoint, SelectionPolicy, UncertaintyModel,
};
use fs_material::state_point::THERMAL_CONDUCTIVITY_DIMS as KD;
use fs_material::thermal::conduction::phase::PreparedPhaseConductivity;
use fs_material::thermal::source::ThermalMaterialPoint;
use fs_qty::{Dims, QuantitySpec, Temperature};

fn card(axis: &str, values: &[(f64, f64)]) -> MaterialCard {
    let mut claims = ClaimSet::new();
    claims
        .insert_claim(PropertyClaim {
            key: PropertyKey::new("k", KD),
            value: PropertyValue::Curve {
                abscissa: axis.into(),
                abscissa_dims: if axis == "T" {
                    Temperature::DIMS
                } else {
                    Dims::NONE
                },
                knots: values.to_vec(),
                dims: KD,
            },
            validity: ValidityDomain::unconstrained()
                .with_quantity(
                    "T",
                    QuantitySpec::dimensional(Temperature::DIMS),
                    300.0,
                    800.0,
                )
                .with_quantity("f", QuantitySpec::dimensional(Dims::NONE), 0.0, 1.0),
            interpolation: InterpolationPolicy::LinearInside,
            uncertainty: UncertaintyModel::Unstated,
            observations: vec![],
            provenance: Provenance {
                source: "manufactured conductive density-volume fixture".into(),
                license: "CC0-1.0".into(),
                artifact: None,
            },
        })
        .unwrap();
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: "manufactured thermal-volume material".into(),
            phase: "declared original caloric chart".into(),
            process: "numerical fixture".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap()
}

fn curve(card: &MaterialCard, values: &[(f64, f64, f64, f64)]) -> EquilibriumEnthalpyPhaseCurve {
    let knots = values
        .iter()
        .map(|&(h, t, f, rho)| EnthalpyPhaseKnot {
            specific_enthalpy_j_kg: h,
            temperature_k: t,
            liquid_mass_fraction: f,
            bulk_density_kg_m3: rho,
        })
        .collect();
    if values.iter().all(|value| value.2 == 1.0) {
        EquilibriumEnthalpyPhaseCurve::try_single_phase(
            card.content_hash(),
            SolidLiquidPhase::Liquid,
            knots,
        )
        .unwrap()
    } else if values.iter().all(|value| value.2 == 0.0) {
        EquilibriumEnthalpyPhaseCurve::try_single_phase(
            card.content_hash(),
            SolidLiquidPhase::Solid,
            knots,
        )
        .unwrap()
    } else {
        EquilibriumEnthalpyPhaseCurve::try_new(card.content_hash(), knots).unwrap()
    }
}

fn source<'a>(
    card: &MaterialCard,
    curve: &'a EquilibriumEnthalpyPhaseCurve,
    bounds: [f64; 2],
) -> PreparedPhaseConductivity<'a> {
    let context = QueryPoint::new()
        .with_quantity("T", QuantitySpec::dimensional(Temperature::DIMS), 300.0)
        .unwrap()
        .with_quantity("f", QuantitySpec::dimensional(Dims::NONE), 0.0)
        .unwrap();
    PreparedPhaseConductivity::new(
        curve,
        ThermalMaterialPoint {
            material: card,
            context: &context,
            temperature_axis: "T",
            liquid_fraction_axis: Some("f"),
        },
        &PropertyKey::new("k", KD),
        ClaimSelection::Policy(SelectionPolicy::SingleClaimOnly),
        bounds,
        || Ok(()),
    )
    .unwrap()
}

fn limits() -> ImplicitHeatWorkThermalVolumeControls {
    ImplicitHeatWorkThermalVolumeControls {
        heat_work: controls(),
        maximum_conductive_faces: 100,
    }
}

fn solve_limits() -> ImplicitThermalVolumeSolveControls {
    let mut solve = ImplicitThermalVolumeSolveControls {
        maximum_residual_evaluations: 4000,
        cell_tolerance_j: 1e-12,
        ..ImplicitThermalVolumeSolveControls::default()
    };
    solve.linear.maximum_iterations = 4000;
    solve
}

fn check_global(step: &ImplicitHeatWorkThermalVolumeStep, external: f64) {
    let volume = &step.evolution.evolution.thermal.volume;
    close(step.external_heat_j, external);
    close(
        volume.report.material_enthalpy_change_j - volume.transport.report.carried_enthalpy_j
            + step.represented_solid_enthalpy_change_j
            + volume.transport.projection.report.pressure_source_work_j,
        external,
    );
    assert!(step.total_heat_defect_j.abs() <= controls().maximum_total_heat_defect_j);
    assert!(step.absolute_local_heat_defect_j <= controls().maximum_absolute_local_heat_defect_j);
}

#[test]
fn original_endpoint_conductivity_and_nonuniform_work_close_heating_and_cooling() {
    for datum in [0.0, 1024.0] {
        let material = card("T", &[(300.0, 1.0), (340.0, 1.8)]);
        let caloric = curve(
            &material,
            &[
                (datum, 300.0, 1.0, 2.0),
                (datum + 40.0, 340.0, 1.0, 1.0 / 0.9),
            ],
        );
        let conductivity = source(&material, &caloric, [datum, datum + 40.0]);
        let h = [datum + 10.0, datum + 20.0];
        let states = h.map(|h| Some(caloric.state_at_specific_enthalpy(h).unwrap()));
        let mass = states.map(|s| s.unwrap().bulk_density_kg_m3());
        let zero = [[0.0; 3]; 2];
        for sign in [-1.0, 1.0] {
            let net = [sign, 2.0 * sign];
            let h_star = core::array::from_fn::<_, 2, _>(|i| h[i] + net[i] / mass[i]);
            let source = net.map(|q| 0.01 * q / DT);
            let flux = [source[0], source[0] + source[1]];
            let p1 = 0.5 * mass[1] * flux[1] / DT;
            let pressure = [p1 + 0.5 * (mass[0] + mass[1]) * flux[0] / DT, p1];
            let work = core::array::from_fn::<_, 2, _>(|i| DT * pressure[i] * source[i]);
            let k = h_star.map(|h| 1.0 + 0.02 * (h - datum));
            let heat01 = DT * (h_star[1] - h_star[0]) / (0.5 / k[0] + 0.5 / k[1]);
            let conductive = [heat01, -heat01];
            let external = core::array::from_fn::<_, 2, _>(|i| net[i] + work[i] - conductive[i]);
            let frozen_heat = DT * 10.0 / (0.5 / 1.2 + 0.5 / 1.4);
            assert!((heat01 - frozen_heat).abs() > 0.05);
            for axis in 0..3 {
                let mut shape = [1; 3];
                shape[axis] = 2;
                let domain = VoxelDomain::new(shape[0], shape[1], shape[2], 1.0).unwrap();
                let field = FlowField::quiescent(&domain);
                let mut faces = [FlowFace::Wall; 6];
                let outlet = 2 * axis + 1;
                faces[outlet] = FlowFace::Free;
                let result = field
                    .advance_isotropic_implicit_heat_work_thermal_volume_momentum_with_baths(
                        &domain,
                        faces,
                        &conductivity,
                        ThermalVolumeBoundaries::default(),
                        &[],
                        &[],
                        input(&mass, &states, &external, &zero, &zero, outlet),
                        &external,
                        limits(),
                        solve_limits(),
                        &mut ImplicitThermalVolumeWorkspace::new(2).unwrap(),
                        &CancelGate::new_clock_free(),
                    )
                    .unwrap();
                let volume = &result.evolution.evolution.thermal.volume;
                for c in 0..2 {
                    let heat = result.evolution.conductive_enthalpy_increment_j[c];
                    close(heat, conductive[c]);
                    close(
                        result.authored_net_enthalpy_increment_j[c],
                        external[c] - work[c],
                    );
                    close(
                        h[c] - datum
                            + (result.authored_net_enthalpy_increment_j[c] + heat) / mass[c],
                        h_star[c] - datum,
                    );
                    close(result.cell_pressure_work_j[c], work[c]);
                    close(volume.volume_source_m3_s[c], source[c]);
                    close(
                        volume.transport.projection.pressure_correction_pa[c],
                        pressure[c],
                    );
                    let returned = volume.phase_states[c].unwrap();
                    close(
                        volume.transport.mass_kg[c] / returned.bulk_density_kg_m3(),
                        1.0,
                    );
                    assert_eq!(returned.phase_curve_identity(), caloric.identity());
                }
                check_global(&result, external.iter().sum());
                assert_eq!(
                    result
                        .evolution
                        .report
                        .conduction
                        .conductivity_source_identity,
                    conductivity.identity()
                );
                assert!(result.work.heat_trials > 1);
                assert_eq!(
                    result.work.heat_trials,
                    result.heat_work.heat_work_iterations
                );
                assert!(
                    result.work.residual_evaluations
                        > result
                            .evolution
                            .report
                            .conduction
                            .solve
                            .residual_evaluations
                );
                assert!(volume.transport.projection.report.projection_loss_j > 0.0);
                // Net authored graph input includes work. Physical external
                // heat is independently returned, without heating any loss.
                close(
                    result
                        .evolution
                        .report
                        .conduction
                        .solve
                        .energy
                        .external_heat_j,
                    external.iter().sum::<f64>() - work.iter().sum::<f64>(),
                );
                assert_eq!(field, FlowField::quiescent(&domain));
            }
        }
    }
}

#[test]
fn finite_sensible_and_latent_walls_share_endpoint_heat_and_pay_native_work() {
    let fluid_card = card("T", &[(300.0, 1.0), (400.0, 11.0)]);
    let fluid_curve = curve(
        &fluid_card,
        &[(0.0, 300.0, 1.0, 2.0), (100.0, 400.0, 1.0, 2.0 / 3.0)],
    );
    let fluid = source(&fluid_card, &fluid_curve, [0.0, 100.0]);
    let mut domain = VoxelDomain::graded([vec![1.0, 2.0], vec![1.0], vec![1.0]]).unwrap();
    domain.set(1, 0, 0, Voxel::Solid(0));
    let field = FlowField::quiescent(&domain);
    let states = [
        Some(fluid_curve.state_at_specific_enthalpy(10.0).unwrap()),
        None,
    ];
    let mass = [states[0].unwrap().bulk_density_kg_m3(), 0.0];
    let zero = [[0.0; 3]; 2];
    let authored = [0.0; 2];
    let mut faces = [FlowFace::Wall; 6];
    faces[0] = FlowFace::Free;
    // One-cell half-face inertia: h*=11, Delta H=M, S=.01*M/dt,
    // p=(M/2)*S/dt. Original wall/bath source resistances below
    // manufacture exactly this coupled endpoint, including mechanical work.
    let volume_source = 0.01 * mass[0] / DT;
    let pressure = 0.5 * mass[0] * volume_source / DT;
    let work = DT * pressure * volume_source;
    let into_fluid = mass[0] + work;
    for latent in [false, true] {
        let solid_card = if latent {
            card("f", &[(0.0, 1.0), (1.0, 2.0)])
        } else {
            card("T", &[(300.0, 1.0), (400.0, 3.0)])
        };
        let solid_curve = if latent {
            curve(
                &solid_card,
                &[
                    (0.0, 300.0, 0.0, 2.0),
                    (50.0, 350.0, 0.0, 2.0),
                    (150.0, 350.0, 1.0, 2.0),
                    (200.0, 400.0, 1.0, 2.0),
                ],
            )
        } else {
            curve(
                &solid_card,
                &[(0.0, 300.0, 0.0, 2.0), (200.0, 400.0, 0.0, 2.0)],
            )
        };
        let solid = source(&solid_card, &solid_curve, [0.0, 200.0]);
        let (old_h, target_h, target_k, bath_t) = if latent {
            (100.0, 99.99, 1.4999, 400.0)
        } else {
            (102.0, 100.0, 2.0, 300.0)
        };
        let solid_change = 4.0 * (target_h - old_h);
        let solid_authored = if latent { 0.25 } else { 1.25 };
        let bath_heat = into_fluid + solid_change - solid_authored;
        let conductance = into_fluid / (DT * (350.0 - 311.0));
        let contact_r = 1.0 / conductance - 0.5 / 2.1 - 1.0 / target_k;
        let bath_r = DT * (bath_t - 350.0) / bath_heat - 1.0 / target_k;
        assert!(contact_r > 0.0 && bath_r > 0.0);
        let walls = [ThermalVolumeSolidCell {
            cell: 1,
            conductivity: &solid,
            mass_kg: 4.0,
            phase: solid_curve.state_at_specific_enthalpy(old_h).unwrap(),
            net_enthalpy_increment_j: solid_authored,
        }];
        let contacts = [ThermalVolumeContact::Resistance(FaceContactResistance {
            cells: [0, 1],
            resistance_m2_k_w: contact_r,
        })];
        let baths = [ThermalVolumeBathFace {
            cell: 1,
            face: Face3::XMax,
            temperature_k: bath_t,
            resistance_m2_k_w: bath_r,
        }];
        let mut data = input(&mass, &states, &authored, &zero, &zero, 0);
        data.volume.boundary_inflow[0] = states[0];
        let result = field
            .advance_isotropic_implicit_heat_work_thermal_volume_momentum_with_baths(
                &domain,
                faces,
                &fluid,
                (&baths[..]).into(),
                &walls,
                &contacts,
                data,
                &authored,
                limits(),
                solve_limits(),
                &mut ImplicitThermalVolumeWorkspace::new(2).unwrap(),
                &CancelGate::new_clock_free(),
            )
            .unwrap();
        let returned = &result.evolution.solid_endpoints[0];
        close(returned.phase.specific_enthalpy_j_kg(), target_h);
        close(returned.phase.temperature_k(), 350.0);
        if latent {
            close(returned.phase.liquid_mass_fraction(), 0.4999);
        }
        assert_eq!(returned.mass_kg, 4.0);
        assert_eq!(returned.conductivity_source_identity, solid.identity());
        assert_eq!(
            returned.phase.phase_curve_identity(),
            solid_curve.identity()
        );
        close(
            result.evolution.conductive_enthalpy_increment_j[0],
            into_fluid,
        );
        assert_eq!(result.evolution.conductive_enthalpy_increment_j[1], 0.0);
        close(result.cell_pressure_work_j[0], work);
        close(result.authored_net_enthalpy_increment_j[0], -work);
        close(result.represented_solid_enthalpy_change_j, solid_change);
        close(
            returned.conductive_enthalpy_increment_j,
            solid_change - solid_authored,
        );
        close(
            result.evolution.report.solid_authored_heat_j,
            solid_authored,
        );
        close(
            result.evolution.report.conduction.bath_exchanges[0].heat_j,
            bath_heat,
        );
        assert_eq!(result.evolution.report.conduction.bath_exchanges.len(), 1);
        close(result.heat_work.supplied_heat_j, into_fluid);
        close(
            result
                .evolution
                .report
                .conduction
                .solve
                .energy
                .external_heat_j,
            bath_heat + solid_authored - work,
        );
        check_global(&result, bath_heat + solid_authored);
        assert!(result.work.heat_trials > 1);
        assert_eq!(walls[0].phase.specific_enthalpy_j_kg(), old_h);
        assert_eq!(field, FlowField::quiescent(&domain));
    }
}

#[test]
fn zero_expansion_and_zero_conduction_recover_the_original_native_owners() {
    let domain = VoxelDomain::new(2, 1, 1, 1.0).unwrap();
    let field = FlowField::quiescent(&domain);
    let mut faces = [FlowFace::Wall; 6];
    faces[1] = FlowFace::Free;
    let zero = [[0.0; 3]; 2];
    let supplied = [1.0, 2.0];
    let gate = CancelGate::new_clock_free();
    for beta in [0.0, 0.01] {
        let k = if beta == 0.0 { 1.0 } else { 0.0 };
        let material = card("T", &[(300.0, k), (340.0, k)]);
        let caloric = curve(
            &material,
            &[
                (0.0, 300.0, 1.0, 2.0),
                (40.0, 340.0, 1.0, 1.0 / (0.5 + beta * 40.0)),
            ],
        );
        let conductivity = source(&material, &caloric, [0.0, 40.0]);
        let states = [10.0, 20.0].map(|h| Some(caloric.state_at_specific_enthalpy(h).unwrap()));
        let mass = states.map(|s| s.unwrap().bulk_density_kg_m3());
        let data = input(&mass, &states, &supplied, &zero, &zero, 1);
        let mut workspace = ImplicitThermalVolumeWorkspace::new(2).unwrap();
        let actual = field
            .advance_isotropic_implicit_heat_work_thermal_volume_momentum_with_baths(
                &domain,
                faces,
                &conductivity,
                ThermalVolumeBoundaries::default(),
                &[],
                &[],
                data,
                &supplied,
                limits(),
                solve_limits(),
                &mut workspace,
                &gate,
            )
            .unwrap();
        if beta == 0.0 {
            let expected = field
                .advance_isotropic_implicit_conjugate_thermal_volume_momentum_with_baths(
                    &domain,
                    faces,
                    None,
                    &conductivity,
                    ThermalVolumeBoundaries::default(),
                    &[],
                    &[],
                    data,
                    ConductiveThermalVolumeMomentumControls {
                        evolution: controls().native,
                        max_conductive_faces: 100,
                    },
                    solve_limits(),
                    &mut workspace,
                    &gate,
                )
                .unwrap();
            assert_eq!(actual.evolution, expected);
            assert_eq!(actual.work.heat_trials, 1);
            assert_eq!(actual.heat_work.source_pressure_work_j, 0.0);
        } else {
            let expected = field
                .advance_supplied_heat_thermal_volume_momentum(
                    &domain,
                    faces,
                    &caloric,
                    data,
                    &supplied,
                    controls(),
                    &gate,
                )
                .unwrap();
            assert_eq!(actual.evolution.evolution, expected.evolution);
            assert_eq!(actual.heat_work, expected.report);
            assert_eq!(
                actual.authored_net_enthalpy_increment_j,
                expected.net_enthalpy_increment_j
            );
        }
        check_global(&actual, 3.0);
    }
}

#[test]
fn rejected_outer_trials_consume_thermal_work_and_publish_no_state() {
    let material = card("T", &[(300.0, 1.0), (340.0, 1.8)]);
    let caloric = curve(
        &material,
        &[(0.0, 300.0, 1.0, 2.0), (40.0, 340.0, 1.0, 1.0 / 0.9)],
    );
    let conductivity = source(&material, &caloric, [0.0, 40.0]);
    let domain = VoxelDomain::new(2, 1, 1, 1.0).unwrap();
    let field = FlowField::quiescent(&domain);
    let states = [10.0, 20.0].map(|h| Some(caloric.state_at_specific_enthalpy(h).unwrap()));
    let mass = states.map(|s| s.unwrap().bulk_density_kg_m3());
    let zero = [[0.0; 3]; 2];
    let supplied = [0.0; 2];
    let data = input(&mass, &states, &supplied, &zero, &zero, 1);
    let mut faces = [FlowFace::Wall; 6];
    faces[1] = FlowFace::Free;
    let gate = CancelGate::new_clock_free();
    let mut workspace = ImplicitThermalVolumeWorkspace::new(2).unwrap();
    let first = field
        .advance_isotropic_implicit_conjugate_thermal_volume_momentum_with_baths(
            &domain,
            faces,
            None,
            &conductivity,
            ThermalVolumeBoundaries::default(),
            &[],
            &[],
            data,
            ConductiveThermalVolumeMomentumControls {
                evolution: controls().native,
                max_conductive_faces: 100,
            },
            solve_limits(),
            &mut workspace,
            &gate,
        )
        .unwrap();
    let run = |data, limits, solve, workspace: &mut _, gate: &_| {
        field.advance_isotropic_implicit_heat_work_thermal_volume_momentum_with_baths(
            &domain,
            faces,
            &conductivity,
            ThermalVolumeBoundaries::default(),
            &[],
            &[],
            data,
            &supplied,
            limits,
            solve,
            workspace,
            gate,
        )
    };
    let expected = run(data, limits(), solve_limits(), &mut workspace, &gate).unwrap();
    assert!(expected.work.heat_trials > 1);
    assert!(expected.heat_work.source_pressure_work_j > 0.0);
    check_global(&expected, 0.0);
    let mut short_solve = solve_limits();
    short_solve.maximum_residual_evaluations = first.report.conduction.solve.residual_evaluations;
    assert!(matches!(
        run(data, limits(), short_solve, &mut workspace, &gate),
        Err(ChtError::InvalidInput { reason, .. }) if reason.contains("whole-call thermal work")
    ));
    let mut short_outer = limits();
    short_outer.heat_work.maximum_heat_work_iterations = 1;
    assert!(run(data, short_outer, solve_limits(), &mut workspace, &gate).is_err());
    let unsupported = [1e6; 2];
    assert!(
        run(
            input(&mass, &states, &unsupported, &zero, &zero, 1),
            limits(),
            solve_limits(),
            &mut workspace,
            &gate
        )
        .is_err()
    );
    let cancelled = CancelGate::new_clock_free();
    cancelled.request();
    assert!(matches!(
        run(data, limits(), solve_limits(), &mut workspace, &cancelled),
        Err(ChtError::Cancelled)
    ));
    assert_eq!(
        run(data, limits(), solve_limits(), &mut workspace, &gate).unwrap(),
        expected
    );
    assert_eq!(field, FlowField::quiescent(&domain));
    assert_eq!(
        states,
        [10.0, 20.0].map(|h| Some(caloric.state_at_specific_enthalpy(h).unwrap()))
    );
}
