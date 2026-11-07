//! Thermal viscosity must close on the native heat endpoint and original history.

mod rheology;

use super::*;
use crate::conjugate::flow::ThermalViscosityFeedbackControls;
use fs_material::liquid::{KINEMATIC_VISCOSITY_DIMS, KINEMATIC_VISCOSITY_PROPERTY as NU_KEY};
use fs_material::liquid_temperature::LiquidViscosityBasis;
use fs_material::phase::HeatCapacitySourceQuery;
use fs_material::phase::melting::{MeltingSourceQuery, ResolvedMeltingCurve};
use fs_material::phase::mushy::MushyPhaseLaw;
use fs_material::phase::transport::{
    PhaseTransportSourceQuery, ResolvedPhaseTransportState,
    resolve_phase_transport_state_with_shear_for_volume_change,
};
use fs_material::state_point::{
    MaterialPropertySelection, SPECIFIC_HEAT_CAPACITY_DIMS,
    SPECIFIC_HEAT_CAPACITY_PROPERTY as CP_KEY, THERMAL_CONDUCTIVITY_PROPERTY as K_KEY,
};
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_qty::{Density, Frequency};

fn absolute_temperature() -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(
        QuantityKind::AbsoluteTemperature,
        ValueForm::Static,
    ))
}

fn material_point(temperature: f64) -> QueryPoint {
    QueryPoint::new()
        .with_quantity("T", absolute_temperature(), temperature)
        .unwrap()
        .with_quantity("rate", QuantitySpec::dimensional(Frequency::DIMS), 0.0)
        .unwrap()
}

/// Manufactured, finite-support original source: rho=2, Cp=1, k=2,
/// T=300+h, mu(h)=1/8+6(h-1), with NO viscosity extrapolation.
fn thermal_material() -> (MaterialCard, ResolvedMeltingCurve) {
    thermal_material_with_rate(false)
}

fn thermal_material_with_rate(rate_dependent: bool) -> (MaterialCard, ResolvedMeltingCurve) {
    let cp = QuantitySpec::dimensional(SPECIFIC_HEAT_CAPACITY_DIMS);
    let rho = QuantitySpec::dimensional(Density::DIMS);
    let latent = QuantitySpec::dimensional(Dims([2, 0, -2, 0, 0, 0]));
    let nu = QuantitySpec::dimensional(KINEMATIC_VISCOSITY_DIMS);
    let k = QuantitySpec::semantic(SemanticType::new(
        QuantityKind::ThermalConductivity,
        ValueForm::Static,
    ));
    let mut claims = ClaimSet::new();
    for (name, quantity, range, values) in [
        ("solid_cp", cp, [299.0, 300.0], [1.0; 2]),
        ("solid_density", rho, [299.0, 300.0], [4.0; 2]),
        (CP_KEY, cp, [300.0, 302.0], [1.0; 2]),
        ("density", rho, [300.0, 302.0], [2.0; 2]),
        (NU_KEY, nu, [301.0, 301.1], [1.0 / 16.0, 29.0 / 80.0]),
        (K_KEY, k, [299.0, 302.0], [2.0; 2]),
        (
            "melting_point",
            absolute_temperature(),
            [299.0, 302.0],
            [300.0; 2],
        ),
        ("fusion", latent, [299.0, 302.0], [1.0; 2]),
    ] {
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::with_quantity(name, quantity),
                value: if rate_dependent && name == NU_KEY {
                    // rho=2; mu=1/8+4*(T-301)+(sqrt(3)/8)*rate.
                    let nu = |t: f64, rate: f64| {
                        1.0 / 16.0 + 2.0 * (t - 301.0) + 3.0_f64.sqrt() * rate / 16.0
                    };
                    PropertyValue::Surface2D {
                        axes: ["T".into(), "rate".into()],
                        axis_dims: [Temperature::DIMS, Frequency::DIMS],
                        samples: vec![
                            ([range[0], 0.0], nu(range[0], 0.0)),
                            ([range[1], 0.0], nu(range[1], 0.0)),
                            ([range[0], 4.0], nu(range[0], 4.0)),
                            ([range[1], 4.0], nu(range[1], 4.0)),
                        ],
                        triangles: vec![[0, 1, 2], [1, 3, 2]],
                        dims: quantity.dims(),
                    }
                } else {
                    PropertyValue::Curve {
                        abscissa: "T".into(),
                        abscissa_dims: Temperature::DIMS,
                        knots: vec![(range[0], values[0]), (range[1], values[1])],
                        dims: quantity.dims(),
                    }
                },
                validity: ValidityDomain::unconstrained()
                    .with_quantity("T", absolute_temperature(), range[0], range[1])
                    .with_quantity("rate", QuantitySpec::dimensional(Frequency::DIMS), 0.0, 4.0),
                interpolation: InterpolationPolicy::LinearInside,
                uncertainty: UncertaintyModel::Unstated,
                observations: vec![],
                provenance: Provenance {
                    source: "MANUFACTURED finite-support thermoviscous oracle, not measurements"
                        .into(),
                    license: "test-only".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    let material = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "manufactured pure substance".into(),
            phase: "solid/liquid".into(),
            process: "thermal-viscosity endpoint regression".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    let solid_cp = PropertyKey::with_quantity("solid_cp", cp);
    let solid_density = PropertyKey::with_quantity("solid_density", rho);
    let liquid_cp = PropertyKey::with_quantity(CP_KEY, cp);
    let liquid_density = PropertyKey::with_quantity("density", rho);
    let melting_point = PropertyKey::with_quantity("melting_point", absolute_temperature());
    let fusion = PropertyKey::with_quantity("fusion", latent);
    let points = [
        material_point(299.0),
        material_point(300.0),
        material_point(302.0),
    ];
    let melting = ResolvedMeltingCurve::try_new(
        &material,
        MeltingSourceQuery {
            solid: HeatCapacitySourceQuery {
                heat_capacity: &solid_cp,
                density: &solid_density,
                temperature_axis: "T",
                lower: &points[0],
                upper: &points[1],
                selection: SelectionPolicy::SingleClaimOnly,
            },
            liquid: HeatCapacitySourceQuery {
                heat_capacity: &liquid_cp,
                density: &liquid_density,
                temperature_axis: "T",
                lower: &points[1],
                upper: &points[2],
                selection: SelectionPolicy::SingleClaimOnly,
            },
            melting_point: &melting_point,
            latent_heat_fusion: &fusion,
            transition: &points[1],
            transition_selection: SelectionPolicy::SingleClaimOnly,
            phase_selector: None,
        },
        -2.0,
        1e-12,
    )
    .unwrap();
    (material, melting)
}

fn thermal_conductivity<'a>(
    material: &MaterialCard,
    caloric: &'a EquilibriumEnthalpyPhaseCurve,
) -> PreparedPhaseConductivity<'a> {
    PreparedPhaseConductivity::new(
        caloric,
        ThermalMaterialPoint {
            material,
            context: &material_point(301.0),
            temperature_axis: "T",
            liquid_fraction_axis: None,
        },
        &PropertyKey::with_quantity(
            K_KEY,
            QuantitySpec::semantic(SemanticType::new(
                QuantityKind::ThermalConductivity,
                ValueForm::Static,
            )),
        ),
        ClaimSelection::Policy(SelectionPolicy::SingleClaimOnly),
        [0.0, 2.0],
        || Ok(()),
    )
    .unwrap()
}

fn thermal_response(
    material: &MaterialCard,
    law: &MushyPhaseLaw<'_>,
    h: f64,
) -> Result<(ResolvedPhaseTransportState, f64), ChtError> {
    let source = resolve_phase_transport_state_with_shear_for_volume_change(
        law,
        h,
        PhaseTransportSourceQuery {
            source: ThermalMaterialPoint {
                material,
                context: &material_point(301.0),
                temperature_axis: "T",
                liquid_fraction_axis: None,
            },
            viscosity_basis: LiquidViscosityBasis::Kinematic,
            extra_requirements: &[],
            selection: MaterialPropertySelection::SingleClaimOnly,
        },
        "rate",
        0.0,
    )
    .map_err(|error| refuse(error.to_string()))?;
    // Independent bulk response; no Stokes hypothesis or shear-rate feedback.
    Ok((source, 1.0 / 3.0 + 8.0 * (h - 1.0)))
}

fn feedback_controls() -> ThermalViscosityFeedbackControls {
    ThermalViscosityFeedbackControls {
        maximum_trials: 100,
        relaxation: 0.5,
        maximum_heat_defect_j: 1e-11,
    }
}

fn molecular_controls() -> MaterialViscosityControls {
    MaterialViscosityControls {
        max_cells: 100,
        max_sweeps: 20_000,
        relative_tolerance: 1e-12,
        maximum_impulse_defect_kg_m_s: 1e-8,
        maximum_angular_impulse_defect_kg_m2_s: 1e-8,
        maximum_energy_defect_j: 1e-8,
    }
}

fn thermal_solve_controls() -> ImplicitThermalVolumeSolveControls {
    let mut solve = solve_controls();
    solve.cell_tolerance_j = 1e-12;
    solve.maximum_residual_evaluations = 20_000;
    solve.linear.maximum_iterations = 20_000;
    solve
}

fn near(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() < 2e-9 * expected.abs().max(1.0),
        "{actual:.17e} != {expected:.17e}",
    );
}

#[test]
fn sourced_shear_and_bulk_close_on_exact_fluid_and_finite_wall_endpoint() {
    let (material, melting) = thermal_material();
    let caloric = melting.curve();
    let law = MushyPhaseLaw::try_from_melting_curve(&melting, 1e-8).unwrap();
    let conductivity = thermal_conductivity(&material, caloric);
    for finite_wall in [false, true] {
        let n = if finite_wall { 4 } else { 3 };
        let mut domain = VoxelDomain::new(n, 1, 1, 1.0).unwrap();
        if finite_wall {
            domain.set(3, 0, 0, Voxel::Solid(0));
        }
        let field = FlowField::quiescent(&domain);
        let mut mass = vec![2.0; n];
        let mut phases = vec![Some(caloric.state_at_specific_enthalpy(1.0).unwrap()); n];
        let mut momentum = vec![[0.0; 3]; n];
        momentum[0][0] = 2.0;
        momentum[2][0] = -2.0;
        let mut authored = vec![0.0; n];
        let force = vec![[0.0; 3]; n];
        let mut solids = vec![];
        if finite_wall {
            mass[3] = 0.0;
            phases[3] = None;
            authored[2] = -1.0 / 16.0;
            solids.push(ThermalVolumeSolidCell {
                cell: 3,
                conductivity: &conductivity,
                mass_kg: 2.0,
                phase: caloric.state_at_specific_enthalpy(9.0 / 8.0).unwrap(),
                net_enthalpy_increment_j: 0.0,
            });
        }
        let calls = Cell::new(0);
        let result = field
            .advance_isotropic_implicit_conjugate_thermoviscous_thermal_volume_momentum_with_baths(
                &domain,
                None,
                &conductivity,
                ThermalVolumeBoundaries::default(),
                &solids,
                &[],
                ThermalVolumeMomentumInput {
                    volume: input(&mass, &phases, &authored),
                    momentum_kg_m_s: &momentum,
                    acceleration_m_s2: [0.0; 3],
                    force_density_n_m3: &force,
                    boundary_velocity_m_s: [None; 6],
                },
                [FvBoundary::Symmetry; 6],
                super::super::super::momentum::coupled_controls(0.5),
                molecular_controls(),
                thermal_solve_controls(),
                feedback_controls(),
                &mut ImplicitThermalVolumeWorkspace::new(n).unwrap(),
                &CancelGate::new_clock_free(),
                |cell, phase| {
                    assert!(cell < 3);
                    calls.set(calls.get() + 1);
                    thermal_response(&material, &law, phase.specific_enthalpy_j_kg())
                },
            )
            .unwrap();
        // The two original free MAC faces have velocities +1/2,-1/2.
        // With lambda=4*mu/3+zeta, their common endpoint magnitude b obeys
        // b=(1/2)/(1+lambda_edge/4+lambda_center/2).
        // At h=[33/32,17/16,33/32], lambda=[1,3/2,1], so b=1/4.
        // dt=1/2 then gives Q=[1/32,3/16,1/32], NOT the frozen-law heat.
        // Each finite-volume BE row is 2*(h_i-1)=Q_i+sum_neighbors(h_j-h_i).
        // The held wall adds +1/16 J to the right fluid; its authored -1/16
        // keeps the same fluid solution while the wall cools from 9/8 to 35/32.
        let expected_h = [33.0 / 32.0, 17.0 / 16.0, 33.0 / 32.0];
        let expected_mu = [5.0 / 16.0, 0.5, 5.0 / 16.0];
        let expected_bulk = [7.0 / 12.0, 5.0 / 6.0, 7.0 / 12.0];
        let expected_heat = [1.0 / 32.0, 3.0 / 16.0, 1.0 / 32.0];
        let viscous = &result.evolution.evolution.viscosity;
        let volume = &result.evolution.evolution.evolution.thermal.volume;
        let mut signed_defect = 0.0;
        let mut absolute_defect = 0.0;
        for cell in 0..3 {
            let source = result.endpoint_sources[cell].as_ref().unwrap();
            near(
                source.phase_state().specific_enthalpy_j_kg(),
                expected_h[cell],
            );
            near(
                source.transport_state().unwrap().dynamic_viscosity_pa_s(),
                expected_mu[cell],
            );
            near(result.bulk_viscosity_pa_s[cell], expected_bulk[cell]);
            // Receipts must describe the same native temperature as the
            // caloric endpoint, not a fabricated density or post-advection law.
            for property in source.material().properties() {
                let temperature = property
                    .answer()
                    .receipt
                    .query_point
                    .iter()
                    .find(|(axis, _)| axis == "T")
                    .unwrap()
                    .1;
                near(temperature, 300.0 + expected_h[cell]);
            }
            near(viscous.cell_molecular_heat_j[cell], expected_heat[cell]);
            near(result.used_molecular_heat_j[cell], expected_heat[cell]);
            near(
                volume.phase_states[cell].unwrap().specific_enthalpy_j_kg(),
                expected_h[cell],
            );
            near(volume.transport.mass_kg[cell], 2.0);
            near(
                result.evolution.conductive_enthalpy_increment_j[cell],
                2.0 * (expected_h[cell] - 1.0) - authored[cell] - expected_heat[cell],
            );
            let difference =
                result.used_molecular_heat_j[cell] - viscous.cell_molecular_heat_j[cell];
            signed_defect += difference;
            absolute_defect += difference.abs();
        }
        near(viscous.molecular_heat_j, 1.0 / 4.0);
        near(viscous.backward_euler_loss_j, 1.0 / 8.0);
        let report = &result.evolution.report;
        near(
            report.conduction.solve.energy.external_heat_j,
            1.0 / 4.0 + authored.iter().sum::<f64>(),
        );
        near(
            volume.report.material_enthalpy_change_j + report.solid_enthalpy_change_j,
            1.0 / 4.0 + authored.iter().sum::<f64>(),
        );
        if finite_wall {
            near(
                result.evolution.solid_endpoints[0]
                    .phase
                    .specific_enthalpy_j_kg(),
                35.0 / 32.0,
            );
            near(report.solid_to_fluid_heat_j, 1.0 / 16.0);
            near(report.solid_enthalpy_change_j, -1.0 / 16.0);
            assert!(result.endpoint_sources[3].is_none());
            assert_eq!(result.bulk_viscosity_pa_s[3], 0.0);
            assert_eq!(result.used_molecular_heat_j[3], 0.0);
            assert_eq!(solids[0].phase.specific_enthalpy_j_kg(), 9.0 / 8.0);
        }
        assert!(result.feedback.trials > 1);
        assert_eq!(calls.get(), 3 * (result.feedback.trials + 1));
        assert_eq!(result.feedback.source_evaluations, calls.get());
        assert!(result.feedback.molecular_solves > 1);
        assert!(result.feedback.molecular_sweeps <= molecular_controls().max_sweeps);
        assert!(
            result.feedback.heat_residual_evaluations
                <= thermal_solve_controls().maximum_residual_evaluations
        );
        assert!(result.feedback.absolute_heat_coupling_defect_j <= 1e-11);
        assert!((result.feedback.heat_coupling_defect_j - signed_defect).abs() < 1e-15);
        assert!((result.feedback.absolute_heat_coupling_defect_j - absolute_defect).abs() < 1e-15);
        assert!(result.feedback.absolute_enthalpy_defect_j <= 1e-9);
        assert_eq!(field, FlowField::quiescent(&domain));
        assert_eq!(phases[0].unwrap().specific_enthalpy_j_kg(), 1.0);
        assert_eq!(momentum[0], [2.0, 0.0, 0.0]);
    }
}

#[test]
fn refused_feedback_preserves_original_history_and_reuses_workspace_for_retry() {
    let (material, melting) = thermal_material();
    let caloric = melting.curve();
    let law = MushyPhaseLaw::try_from_melting_curve(&melting, 1e-8).unwrap();
    let conductivity = thermal_conductivity(&material, caloric);
    let mut domain = VoxelDomain::new(4, 1, 1, 1.0).unwrap();
    domain.set(3, 0, 0, Voxel::Solid(0));
    let field = FlowField::quiescent(&domain);
    let initial_phase = caloric.state_at_specific_enthalpy(1.0).unwrap();
    let mass = [2.0, 2.0, 2.0, 0.0];
    let phases = [
        Some(initial_phase),
        Some(initial_phase),
        Some(initial_phase),
        None,
    ];
    let momentum = [[2.0, 0.0, 0.0], [0.0; 3], [-2.0, 0.0, 0.0], [0.0; 3]];
    let authored = [0.0, 0.0, -1.0 / 16.0, 0.0];
    let forces = [[0.0; 3]; 4];
    let solids = [ThermalVolumeSolidCell {
        cell: 3,
        conductivity: &conductivity,
        mass_kg: 2.0,
        phase: caloric.state_at_specific_enthalpy(9.0 / 8.0).unwrap(),
        net_enthalpy_increment_j: 0.0,
    }];
    #[derive(Clone, Copy, PartialEq)]
    enum SourceFailure {
        None,
        StaleEndpoint,
        CancelEndpoint,
    }
    let run = |feedback: ThermalViscosityFeedbackControls,
               viscosity: MaterialViscosityControls,
               solve: ImplicitThermalVolumeSolveControls,
               heat: &[f64],
               failure: SourceFailure,
               calls: &Cell<usize>,
               workspace: &mut ImplicitThermalVolumeWorkspace,
               gate: &CancelGate| {
        field.advance_isotropic_implicit_conjugate_thermoviscous_thermal_volume_momentum_with_baths(
            &domain,
            None,
            &conductivity,
            ThermalVolumeBoundaries::default(),
            &solids,
            &[],
            ThermalVolumeMomentumInput {
                volume: input(&mass, &phases, heat),
                momentum_kg_m_s: &momentum,
                acceleration_m_s2: [0.0; 3],
                force_density_n_m3: &forces,
                boundary_velocity_m_s: [None; 6],
            },
            [FvBoundary::Symmetry; 6],
            super::super::super::momentum::coupled_controls(0.5),
            viscosity,
            solve,
            feedback,
            workspace,
            gate,
            |_, phase| {
                calls.set(calls.get() + 1);
                if calls.get() > 3 && failure == SourceFailure::CancelEndpoint {
                    gate.request();
                    return Err(refuse("resolver error after cancellation"));
                }
                let h = phase.specific_enthalpy_j_kg()
                    + if calls.get() > 3 && failure == SourceFailure::StaleEndpoint {
                        0.001
                    } else {
                        0.0
                    };
                thermal_response(&material, &law, h)
            },
        )
    };
    let mut workspace = ImplicitThermalVolumeWorkspace::new(4).unwrap();
    let gate = CancelGate::new_clock_free();
    let calls = Cell::new(0);
    let expected = run(
        feedback_controls(),
        molecular_controls(),
        thermal_solve_controls(),
        &authored,
        SourceFailure::None,
        &calls,
        &mut workspace,
        &gate,
    )
    .unwrap();
    calls.set(0);
    let refused = run(
        ThermalViscosityFeedbackControls {
            maximum_trials: 1,
            ..feedback_controls()
        },
        molecular_controls(),
        thermal_solve_controls(),
        &authored,
        SourceFailure::None,
        &calls,
        &mut workspace,
        &gate,
    )
    .unwrap_err();
    assert!(
        refused
            .to_string()
            .contains("did not converge within its trial budget")
    );
    assert_eq!(calls.get(), 6);

    // Each ceiling is just below the measured aggregate work, but comfortably
    // above one trial's work: resetting the budget every trial would pass.
    assert!(expected.feedback.heat_residual_evaluations > 1);
    let mut bounded_heat = thermal_solve_controls();
    bounded_heat.maximum_residual_evaluations = expected.feedback.heat_residual_evaluations - 1;
    calls.set(0);
    assert!(
        run(
            feedback_controls(),
            molecular_controls(),
            bounded_heat,
            &authored,
            SourceFailure::None,
            &calls,
            &mut workspace,
            &gate,
        )
        .is_err()
    );
    assert!(calls.get() > 3);
    assert!(expected.feedback.molecular_sweeps > 1);
    let mut bounded_molecular = molecular_controls();
    bounded_molecular.max_sweeps = expected.feedback.molecular_sweeps - 1;
    calls.set(0);
    assert!(
        run(
            feedback_controls(),
            bounded_molecular,
            thermal_solve_controls(),
            &authored,
            SourceFailure::None,
            &calls,
            &mut workspace,
            &gate,
        )
        .is_err()
    );
    assert!(calls.get() > 3);

    calls.set(0);
    let stale = run(
        feedback_controls(),
        molecular_controls(),
        thermal_solve_controls(),
        &authored,
        SourceFailure::StaleEndpoint,
        &calls,
        &mut workspace,
        &gate,
    )
    .unwrap_err();
    assert!(
        stale
            .to_string()
            .contains("must match the exact endpoint phase and density")
    );
    assert_eq!(calls.get(), 4);
    // The caloric chart still admits this endpoint; only original viscosity
    // support ends at 301.1 K. Refusal must preserve the bounded source law.
    calls.set(0);
    assert!(
        run(
            feedback_controls(),
            molecular_controls(),
            thermal_solve_controls(),
            &[0.0, 1.0, -1.0 / 16.0, 0.0],
            SourceFailure::None,
            &calls,
            &mut workspace,
            &gate,
        )
        .is_err()
    );
    assert!(calls.get() > 3);

    let cancelled = CancelGate::new_clock_free();
    calls.set(0);
    assert!(matches!(
        run(
            feedback_controls(),
            molecular_controls(),
            thermal_solve_controls(),
            &authored,
            SourceFailure::CancelEndpoint,
            &calls,
            &mut workspace,
            &cancelled,
        ),
        Err(ChtError::Cancelled)
    ));
    assert_eq!(calls.get(), 4);
    assert_eq!(field, FlowField::quiescent(&domain));
    assert_eq!(mass, [2.0, 2.0, 2.0, 0.0]);
    assert_eq!(
        phases,
        [
            Some(initial_phase),
            Some(initial_phase),
            Some(initial_phase),
            None
        ]
    );
    assert_eq!(
        momentum,
        [[2.0, 0.0, 0.0], [0.0; 3], [-2.0, 0.0, 0.0], [0.0; 3]]
    );
    assert_eq!(authored, [0.0, 0.0, -1.0 / 16.0, 0.0]);
    assert_eq!(solids[0].phase.specific_enthalpy_j_kg(), 9.0 / 8.0);
    calls.set(0);
    let retry = run(
        feedback_controls(),
        molecular_controls(),
        thermal_solve_controls(),
        &authored,
        SourceFailure::None,
        &calls,
        &mut workspace,
        &gate,
    )
    .unwrap();
    assert_eq!(retry, expected);
    assert_eq!(calls.get(), retry.feedback.source_evaluations);
}
