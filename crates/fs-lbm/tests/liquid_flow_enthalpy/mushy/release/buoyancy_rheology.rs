//! One original h/T/rate/density law survives phase flow and complete-melt release.
use super::super::rheology::source as shear_card;
use super::*;
use fs_lbm::material::free_surface::enthalpy::conduction::FreeSurfaceConductionFace;
use fs_lbm::material::free_surface::enthalpy::rheology::RheologySolve;
use fs_lbm::material::transport::boussinesq::BoussinesqConfig;
use fs_material::fluid_transport::FluidTransport;
use fs_material::liquid::surface::ResolvedLiquidSurfaceState;
use fs_material::phase::transport::{
    PhaseTransportSourceQuery, resolve_phase_transport_state_with_shear_rate_and_density_departure,
};
use fs_material::thermal::source::ThermalMaterialPoint;
use fs_qty::Frequency;

fn original_source(slope: f64) -> MaterialCard {
    let original = shear_card(slope, 2.0);
    let mut claims = ClaimSet::new();
    for (_, claim) in original.claims().claims_ordered() {
        let mut claim = claim.clone();
        if claim.key.name() == DENSITY_PROPERTY {
            claim.value = PropertyValue::Curve {
                abscissa: "T".into(),
                abscissa_dims: Temperature::DIMS,
                knots: vec![(450.0, 1.01), (500.0, 1.0), (600.0, 0.98)],
                dims: Density::DIMS,
            };
            claim.interpolation = InterpolationPolicy::LinearInside;
        }
        claims.insert_claim(claim).unwrap();
    }
    MaterialCard::assemble(original.id().clone(), claims, vec![]).unwrap()
}

fn solve() -> RheologySolve {
    RheologySolve {
        maximum_iterations: 128,
        relaxation: 0.5,
        absolute_relaxation_tolerance: 1e-13,
        relative_relaxation_tolerance: 1e-13,
    }
}

#[test]
fn buoyant_shear_dependent_melting_releases_original_mass_heat_and_motion_to_native_vof() {
    let mut native_speeds = Vec::new();
    let mut released_speeds = Vec::new();
    for slope in [0.0, 500.0] {
        let card = original_source(slope);
        let curve = chart(&card);
        let mushy = MushyPhaseLaw::try_new_with_immobile_solid(&curve, 1.0, 1.0, 1.0).unwrap();
        let domain = domain();
        let ids = ids(&domain);
        let config = config();
        let gate = CancelGate::new_clock_free();
        // The reference rate deliberately differs from initial and live MAC rates.
        let context = point_liquid(500.0)
            .with_quantity("rate", QuantitySpec::dimensional(Frequency::DIMS), 1.0)
            .unwrap();
        let liquid = LiquidTemperatureLaw::try_new(
            &card,
            &context,
            "T",
            LiquidViscosityBasis::Dynamic,
            MaterialPropertySelection::SingleClaimOnly,
        )
        .unwrap();
        let reference = liquid.reference_state().transport_state().unwrap();
        let surface = ResolvedLiquidSurfaceState::try_for_shear_rate_evolution(
            liquid.reference_state().material().clone(),
            LiquidViscosityBasis::Dynamic,
            "T",
            "rate",
        )
        .unwrap();
        let gravity = [0.2, 0.0, 0.0];
        let bounds = BoussinesqConfig {
            gravity_m_s2: gravity,
            max_relative_density_departure: 0.03,
            max_relative_specific_heat_departure: 0.0,
        };
        let raw_source = |cell: &DepositedCell, rate: f64| {
            resolve_phase_transport_state_with_shear_rate_and_density_departure(
                &mushy,
                cell.specific_h_j_kg,
                PhaseTransportSourceQuery {
                    source: ThermalMaterialPoint {
                        material: &card,
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
                bounds.max_relative_density_departure,
            )
        };
        let original = initial(&domain, &curve, -100.0, 2000.0, 0.15);
        let history = MaterialFlowFeedback::try_new_phase_boussinesq(
            &original,
            &domain,
            FlowEnthalpyInput {
                cell_ids: &ids,
                material_id: MAT,
                caloric: &curve,
                inlet_specific_h_j_kg: [None; 6],
                max_cells: 64,
                max_connections: 1024,
            },
            &config,
            (),
            &gate,
            256,
            &mushy,
            reference,
            bounds,
            |_, cell: &DepositedCell, _| raw_source(cell, 0.0),
        )
        .unwrap();
        let projection = history.endpoint()[0].material.clone();
        let source = |_: usize, cell: &DepositedCell, _: f64, rate: f64| {
            projection.at_source(raw_source(cell, rate)?)
        };
        let history = history
            .with_unsteady_viscous_heating(viscous())
            .unwrap()
            .with_native_rheology(&original, &source)
            .unwrap();
        assert!(history.endpoint().iter().all(|sample| {
            sample.material.material().mobile_transport().is_none()
                && sample
                    .material
                    .material()
                    .material()
                    .property(DYNAMIC_VISCOSITY_PROPERTY)
                    .is_none()
        }));
        let (plateau, history, _) = original
            .advance_observed(&history, 0.075, step(), || Ok(()))
            .unwrap();
        for sample in history.endpoint() {
            near(sample.specific_enthalpy_j_kg, 50.0);
            near(sample.temperature_k, 500.0);
            near(
                sample
                    .material
                    .material()
                    .phase_state()
                    .liquid_mass_fraction(),
                0.5,
            );
            assert_eq!(sample.material.reference(), reference);
            assert_eq!(sample.material.config(), bounds);
        }
        let (molten, history, _) = plateau
            .advance_observed(&history, 0.15, step(), || Ok(()))
            .unwrap();
        close(molten.energy_j() - original.energy_j(), 8.0 * 300.0, 1e-8);
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
        let history = history
            .with_unsteady_body_force_heating(force, viscous())
            .unwrap();
        let (warm, history, first_heat) = molten
            .advance_observed(&history, 0.16, step(), || Ok(()))
            .unwrap();
        assert!(first_heat.external_heat_j > 0.0);
        near(
            warm.energy_j() - molten.energy_j(),
            first_heat.external_heat_j,
        );
        // A second interval consumes the actual native rate from the first.
        let (thermal, history, _) = warm
            .advance_observed(&history, 0.17, step(), || Ok(()))
            .unwrap();
        assert!(
            history
                .endpoint_shear_rate_per_s()
                .unwrap()
                .iter()
                .any(|rate| *rate > 0.0)
        );
        for sample in history.endpoint() {
            let rate = history.endpoint_shear_rate_per_s().unwrap()[sample.voxel];
            let original = sample.material.material();
            near(original.material().query().axes()["rate"], rate);
            near(
                original.source_density_kg_m3(),
                1.0 - 0.0002 * (sample.temperature_k - 500.0),
            );
            near(
                original.transport_state().unwrap().dynamic_viscosity_pa_s(),
                100.0 - 0.05 * (sample.temperature_k - 550.0) + slope * rate,
            );
            assert_eq!(sample.material.reference(), reference);
            assert_eq!(sample.material.config(), bounds);
        }
        native_speeds.push(history.unsteady_flow().unwrap().field.flux_x(2, 1, 0));
        let before = (
            thermal.cells().clone(),
            history.unsteady_flow().unwrap().clone(),
        );
        with_cx(&gate, |cx| {
            let scaling = plan_free_surface_from_material(&surface, 1.0, 0.0001, 0.2).unwrap();
            let mut allowance = limits();
            allowance.maximum_kinetic_projection_deficit_j = 0.1;
            assert!(
                history
                    .release_fully_liquid(
                        &thermal,
                        scaling,
                        receiving(gravity),
                        &liquid,
                        &curve,
                        StationaryConstraintRelease::Withdraw,
                        allowance,
                        FreeSurfaceEnthalpyControls::default(),
                        cx,
                    )
                    .is_err()
            );
            assert!(
                history
                    .release_fully_liquid_with_rheology(
                        &thermal,
                        scaling,
                        receiving([0.1, 0.0, 0.0]),
                        &liquid,
                        &curve,
                        StationaryConstraintRelease::Withdraw,
                        allowance,
                        FreeSurfaceEnthalpyControls::default(),
                        solve(),
                        cx,
                    )
                    .is_err()
            );
            assert!(
                history
                    .release_fully_liquid_with_rheology(
                        &thermal,
                        scaling,
                        receiving(gravity),
                        &liquid,
                        &curve,
                        StationaryConstraintRelease::Withdraw,
                        allowance,
                        FreeSurfaceEnthalpyControls::default(),
                        RheologySolve {
                            maximum_iterations: 0,
                            ..solve()
                        },
                        cx,
                    )
                    .is_err()
            );
            let (mut released, report) = history
                .release_fully_liquid_with_rheology(
                    &thermal,
                    scaling,
                    receiving(gravity),
                    &liquid,
                    &curve,
                    StationaryConstraintRelease::Withdraw,
                    allowance,
                    FreeSurfaceEnthalpyControls::default(),
                    solve(),
                    cx,
                )
                .unwrap();
            close(report.source.mass_kg, 8.0, 1e-12);
            close(report.target.mass_kg, 8.0, 1e-12);
            close(report.source_enthalpy_j, thermal.energy_j(), 1e-9);
            close(report.target_enthalpy_j, thermal.energy_j(), 1e-9);
            for axis in 0..3 {
                close(
                    report.source.momentum_kg_m_s[axis],
                    report.target.momentum_kg_m_s[axis],
                    1e-10,
                );
                close(
                    report.source.angular_momentum_kg_m2_s[axis],
                    report.target.angular_momentum_kg_m2_s[axis],
                    1e-10,
                );
            }
            assert!(report.kinetic_projection_deficit_j >= -allowance.kinetic_roundoff_tolerance_j);
            near(
                report.source.kinetic_energy_j - report.target.kinetic_energy_j,
                report.kinetic_projection_deficit_j,
            );
            assert_eq!(released.physical_time_s(), thermal.time_s());
            assert_eq!(released.accepted_steps(), 0);
            let heating: Vec<_> = ids
                .iter()
                .enumerate()
                .filter_map(|(voxel, id)| {
                    id.map(|_| FreeSurfaceCellHeat {
                        coordinates: domain.coords(voxel),
                        energy_j: 20.0,
                    })
                })
                .collect();
            let heated = released.step_with_heat(&heating, &gate).unwrap();
            close(heated.enthalpy_after_j - thermal.energy_j(), 160.0, 1e-7);
            close(heated.mass_after_kg, 8.0, 1e-11);
            let step = released
                .step_with_conduction_and_viscous_heating(
                    &[FreeSurfaceConductionFace {
                        cells: [[1, 1, 0], [2, 1, 0]],
                        area_m2: 1.0,
                        distances_m: [0.5; 2],
                    }],
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
            let mut expected_heat = 0.0;
            for sample in &step.applied {
                let query = sample.material.material().query();
                let temperature = query.axes()["T"];
                let rate = sample.shear_rate_per_s;
                let mu = 100.0 - 0.05 * (temperature - 550.0) + slope * rate;
                near(query.axes()["rate"], rate);
                near(sample.material.dynamic_viscosity_pa_s(), mu);
                near(
                    sample.material.density_kg_m3(),
                    1.0 - 0.0002 * (temperature - 500.0),
                );
                expected_heat += mu * rate * rate * sample.liquid_volume_m3 * 0.0001;
            }
            assert!(expected_heat > 0.0);
            near(step.conduction.flow.viscous_heat_j, expected_heat);
            near(
                released.total_enthalpy_j().unwrap(),
                thermal.energy_j() + 160.0 + expected_heat,
            );
            let at = [1, 1, 0];
            let material = released.material_at(at).unwrap();
            let rate = released.rheology_at(at).unwrap();
            assert_eq!(material.material().card_identity(), card.content_hash());
            assert!(material.density_kg_m3() < 1.0);
            // Collision uses mu/rho0, never mu divided by the departing source rho.
            near(
                rate.collision_tau,
                0.5 + 3.0 * material.dynamic_viscosity_pa_s() * 0.0001,
            );
            released_speeds.push(
                released
                    .flow()
                    .sample_at_si(1, 1, 0)
                    .unwrap()
                    .unwrap()
                    .velocity_m_s,
            );
        });
        assert_eq!(thermal.cells(), &before.0);
        assert_eq!(history.unsteady_flow().unwrap(), &before.1);
    }
    assert!((native_speeds[0] - native_speeds[1]).abs() > 1e-9);
    assert!(
        released_speeds[0]
            .iter()
            .zip(released_speeds[1])
            .any(|(a, b)| (a - b).abs() > 1e-9)
    );
}
