//! The same h/T/shear source reaches actual VOF flow, heat and conduction.
use super::super::rheology::{resolve_at, source as shear_card};
use super::*;
use fs_lbm::conjugate::enthalpy::advection::feedback::rheology::FlowRheologySource;
use fs_lbm::material::free_surface::enthalpy::conduction::FreeSurfaceConductionFace;
use fs_lbm::material::free_surface::enthalpy::rheology::RheologySolve;
use fs_material::liquid::surface::ResolvedLiquidSurfaceState;
use fs_material::phase::conductivity::LiquidMassFractionConductivity;
use fs_qty::Frequency;

fn solve() -> RheologySolve {
    RheologySolve {
        maximum_iterations: 128,
        relaxation: 0.5,
        absolute_relaxation_tolerance: 1e-13,
        relative_relaxation_tolerance: 1e-13,
    }
}

#[allow(clippy::too_many_arguments)]
fn history<'a>(
    thermal: &DepositionSimulation<'_>,
    domain: &'a VoxelDomain,
    ids: &'a [Option<CellId>],
    config: &'a SimpleConfig,
    gate: &'a CancelGate,
    card: &'a MaterialCard,
    mushy: &'a MushyPhaseLaw<'a>,
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
            caloric: mushy.curve(),
            inlet_specific_h_j_kg: [None; 6],
            max_cells: 64,
            max_connections: 1024,
        },
        config,
        (),
        gate,
        256,
        mushy,
        move |_, cell: &DepositedCell, _| resolve_at(card, mushy, cell.specific_h_j_kg, 0.0),
    )
    .unwrap()
    .with_unsteady_viscous_heating(viscous())
    .unwrap()
    .with_native_rheology(thermal, source)
    .unwrap()
}

#[test]
fn molten_native_shear_response_survives_release_and_drives_vof_heat_and_motion() {
    let mut outcomes = Vec::new();
    for slope in [0.0, 500.0] {
        let card = shear_card(slope, 2.0);
        let curve = chart(&card);
        let mushy = MushyPhaseLaw::try_new_with_immobile_solid(&curve, 1.0, 1.0, 1.0).unwrap();
        let domain = domain();
        let ids = ids(&domain);
        let config = config();
        let gate = CancelGate::new_clock_free();
        let original = initial(&domain, &curve, -100.0, 3000.0, 0.1);
        let source = |_: usize, cell: &DepositedCell, _: f64, rate: f64| {
            resolve_at(&card, &mushy, cell.specific_h_j_kg, rate)
        };
        let bound = history(
            &original, &domain, &ids, &config, &gate, &card, &mushy, &source,
        );
        let (molten, bound, _) = original
            .advance_observed(&bound, 0.1, step(), || Ok(()))
            .unwrap();
        near(molten.energy_j() - original.energy_j(), 8.0 * 300.0);
        let force = (0..domain.cell_count())
            .map(|c| {
                if !domain.is_fluid(c) {
                    return [0.0; 3];
                }
                let [x, y, z] = domain.coords(c);
                let p = domain.center(x, y, z);
                [-100.0 * (p[1] - 2.0), 100.0 * (p[0] - 2.0), 0.0]
            })
            .collect();
        let bound = bound
            .with_unsteady_body_force_heating(force, viscous())
            .unwrap();
        let (thermal, bound, _) = molten
            .advance_observed(&bound, 0.11, step(), || Ok(()))
            .unwrap();
        let speed = bound.unsteady_flow().unwrap().field.flux_x(2, 1, 0);
        assert!(speed > 0.0);
        assert!(
            bound
                .endpoint_shear_rate_per_s()
                .unwrap()
                .iter()
                .any(|rate| *rate > 0.0)
        );
        // Deliberately use a different REFERENCE rate. The transfer must query
        // each actual MAC rate and the receiver must use its own solved stress.
        let point = point_liquid(550.0)
            .with_quantity("rate", QuantitySpec::dimensional(Frequency::DIMS), 1.0)
            .unwrap();
        let liquid_law = LiquidTemperatureLaw::try_new(
            &card,
            &point,
            "T",
            LiquidViscosityBasis::Dynamic,
            MaterialPropertySelection::SingleClaimOnly,
        )
        .unwrap();
        let surface = ResolvedLiquidSurfaceState::try_for_shear_rate_evolution(
            liquid_law.reference_state().material().clone(),
            LiquidViscosityBasis::Dynamic,
            "T",
            "rate",
        )
        .unwrap();
        with_cx(&gate, |cx| {
            let scaling = plan_free_surface_from_material(&surface, 1.0, 0.0001, 0.2).unwrap();
            let mut allowance = limits();
            allowance.maximum_kinetic_projection_deficit_j = 0.1;
            let before = (
                thermal.cells().clone(),
                bound.unsteady_flow().unwrap().clone(),
            );
            assert!(
                bound
                    .release_fully_liquid(
                        &thermal,
                        scaling,
                        receiving([0.0; 3]),
                        &liquid_law,
                        &curve,
                        StationaryConstraintRelease::Withdraw,
                        allowance,
                        FreeSurfaceEnthalpyControls::default(),
                        cx
                    )
                    .is_err()
            );
            assert!(
                bound
                    .release_fully_liquid_with_rheology(
                        &thermal,
                        scaling,
                        receiving([0.0; 3]),
                        &liquid_law,
                        &curve,
                        StationaryConstraintRelease::Withdraw,
                        allowance,
                        FreeSurfaceEnthalpyControls::default(),
                        RheologySolve {
                            maximum_iterations: 0,
                            ..solve()
                        },
                        cx
                    )
                    .is_err()
            );
            // An extra phase multiplier cannot silently vanish at release.
            let enhanced = bound
                .clone()
                .with_phase_conductivity(
                    &curve,
                    LiquidMassFractionConductivity::declared(4.0, 2.0).unwrap(),
                )
                .unwrap();
            assert!(
                matches!(enhanced.release_fully_liquid_with_rheology(&thermal, scaling,
                receiving([0.0; 3]), &liquid_law, &curve, StationaryConstraintRelease::Withdraw,
                allowance, FreeSurfaceEnthalpyControls::default(), solve(), cx),
                Err(LiquidReleaseError::Invalid(message)) if message.contains("phase-conductivity"))
            );
            let (mut released, report) = bound
                .release_fully_liquid_with_rheology(
                    &thermal,
                    scaling,
                    receiving([0.0; 3]),
                    &liquid_law,
                    &curve,
                    StationaryConstraintRelease::Withdraw,
                    allowance,
                    FreeSurfaceEnthalpyControls::default(),
                    solve(),
                    cx,
                )
                .unwrap();
            near(report.source.kinetic_energy_j, 4.0 * speed * speed);
            near(report.target.kinetic_energy_j, 2.0 * speed * speed);
            near(report.kinetic_projection_deficit_j, 2.0 * speed * speed);
            near(report.target_enthalpy_j, thermal.energy_j());
            near(report.target.angular_momentum_kg_m2_s[2], 4.0 * speed);
            near(report.source.angular_momentum_kg_m2_s[2], 4.0 * speed);
            assert_eq!(released.physical_time_s(), thermal.time_s());
            assert_eq!(released.accepted_steps(), 0);
            let first = released
                .step_with_heat(
                    &[FreeSurfaceCellHeat {
                        coordinates: [1, 1, 0],
                        energy_j: 1.0,
                    }],
                    &gate,
                )
                .unwrap();
            near(first.enthalpy_after_j - thermal.energy_j(), 1.0);
            let face = FreeSurfaceConductionFace {
                cells: [[1, 1, 0], [2, 1, 0]],
                area_m2: 1.0,
                distances_m: [0.5; 2],
            };
            let mut total_heat = 0.0;
            for _ in 0..4 {
                let step = released
                    .step_with_conduction_and_viscous_heating(
                        &[face],
                        &[],
                        FreeSurfaceConductionControls {
                            maximum_faces: 1,
                            maximum_substeps: 64,
                            maximum_dt_s: 0.0001,
                            energy_tolerance_j: 1e-8,
                        },
                        &gate,
                    )
                    .unwrap();
                let mut independently_summed_heat = 0.0;
                for sample in &step.applied {
                    let state = sample.material.material().query();
                    let rate = sample.shear_rate_per_s;
                    near(state.axes()["rate"], rate);
                    let mu = 100.0 - 0.05 * (state.axes()["T"] - 550.0) + slope * rate;
                    near(sample.material.dynamic_viscosity_pa_s(), mu);
                    independently_summed_heat +=
                        mu * rate * rate * sample.liquid_volume_m3 * 0.0001;
                }
                assert!(independently_summed_heat > 0.0);
                near(
                    step.conduction.flow.viscous_heat_j,
                    independently_summed_heat,
                );
                let applied = &step.conduction.last_applied[0];
                let k = applied.materials.each_ref().map(|material| {
                    let q = material.material().query();
                    0.2 + 0.001 * (q.axes()["T"] - 500.0) + 0.05 * q.axes()["rate"]
                });
                near(applied.conductance_w_k, 1.0 / (0.5 / k[0] + 0.5 / k[1]));
                near(step.conduction.flow.external_heat_j, 0.0);
                total_heat += independently_summed_heat;
            }
            near(
                released.total_enthalpy_j().unwrap(),
                thermal.energy_j() + 1.0 + total_heat,
            );
            assert_eq!(released.accepted_steps(), 5);
            close(released.physical_time_s(), thermal.time_s() + 0.0005, 1e-15);
            let at = [1, 1, 0];
            let response = released.rheology_at(at).unwrap();
            let material = released.material_at(at).unwrap();
            assert_eq!(material.material().card_identity(), card.content_hash());
            near(
                response.collision_tau,
                0.5 + 3.0 * material.dynamic_viscosity_pa_s() * 0.0001,
            );
            assert_eq!(
                material
                    .material()
                    .property(DYNAMIC_VISCOSITY_PROPERTY)
                    .unwrap()
                    .answer()
                    .receipt
                    .selected,
                liquid_law
                    .reference_state()
                    .material()
                    .property(DYNAMIC_VISCOSITY_PROPERTY)
                    .unwrap()
                    .answer()
                    .receipt
                    .selected
            );
            let velocity = released
                .flow()
                .sample_at_si(1, 1, 0)
                .unwrap()
                .unwrap()
                .velocity_m_s;
            outcomes.push(velocity);
            let rejected = (
                released.physical_time_s(),
                released.total_enthalpy_j().unwrap(),
                released.flow().lattice().grid().state_bits(),
            );
            assert!(
                released
                    .step_with_heat(
                        &[FreeSurfaceCellHeat {
                            coordinates: at,
                            energy_j: 1e6,
                        }],
                        &gate
                    )
                    .is_err()
            );
            assert_eq!(released.physical_time_s(), rejected.0);
            assert_eq!(released.total_enthalpy_j().unwrap(), rejected.1);
            assert_eq!(released.flow().lattice().grid().state_bits(), rejected.2);
            assert_eq!(thermal.cells(), &before.0);
            assert_eq!(bound.unsteady_flow().unwrap(), &before.1);
        });
    }
    assert!(
        outcomes[0]
            .iter()
            .zip(outcomes[1])
            .any(|(a, b)| (a - b).abs() > 1e-12)
    );
}
