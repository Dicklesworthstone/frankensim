//! G1/G3: original finite ductile work, reciprocal caloric heat and retained motion.
//! These explicitly authored constants and low heat capacity are manufactured
//! numerical fixtures, not experimental metal/fracture or heat-partition data.
use super::*;
use fs_couple::solid_enthalpy::finite_ductile::{
    DuctileHeatPartition, FiniteDuctileThermalCell, FiniteDuctileThermalSolid,
    ReferenceCaloricMaterial, ReferenceCaloricQuery,
};
use fs_material::ductile::evolution::DuctileJ2Settings;
use fs_material::ductile::source::{
    DuctileSourceRequest, FRACTURE_PLASTIC_STRAIN_PROPERTY, PLASTIC_STRAIN_RATE_AXIS,
    STRESS_TRIAXIALITY_AXIS, resolve_ductile_locus,
};
use fs_material::ductile::{DuctileDomain, JohnsonCookInitiation, JohnsonCookParameters};
use fs_material::plastic::energy::EnergyTolerances;
use fs_material::plastic::finite::ductile::FiniteDuctileJ2;
use fs_material::plastic::finite::{FiniteJ2, FiniteJ2Limits, FiniteJ2Parameters};
use fs_solid::brittle::TetCrackBandGeometry;
use fs_solid::ductile::{DuctileStepConditions, DuctileTemperatures};
use fs_solid::hyper3::plastic::ductile::FiniteDuctileTetSolid;
use fs_solid::hyper3::plastic::ductile::dynamics::FiniteDuctileMotionControls;
use fs_solid::loads3::TetLoadGeometry3;
use fs_solid::plastic3::IncrementSettings;

const MU: f64 = 80.0e6;
const BULK: f64 = 160.0e6;
const YIELD: f64 = 2.0e6;
const GP: f64 = 400.0;
const LOW_CP: f64 = 2.0;
const ONSET: f64 = 0.008;
const RATE: f64 = 0.04;
const TEMPERATURE: f64 = -0.6;
const PARTITION: DuctileHeatPartition = DuctileHeatPartition {
    plastic_dissipation_heat_fraction: 0.75,
    damage_release_heat_fraction: 0.25,
};

fn native() -> FiniteJ2 {
    FiniteJ2::new(
        FiniteJ2Parameters {
            shear_modulus_pa: MU,
            bulk_modulus_pa: BULK,
            yield_stress_pa: YIELD,
            hardening_modulus_pa: 0.0,
        },
        FiniteJ2Limits {
            minimum_jacobian: 0.25,
            maximum_elastic_log_strain: 0.5,
            maximum_plastic_log_strain: 2.0,
            maximum_accumulated_plastic_strain: 3.0,
        },
    )
    .unwrap()
}

fn damage_settings(width: f64) -> DuctileJ2Settings {
    DuctileJ2Settings {
        characteristic_length_m: width,
        post_initiation_plastic_work_j_m2: GP,
        maximum_damage_increment: 1.0,
        energy_tolerances: EnergyTolerances {
            energy_j_m3: 1e-6,
            yield_pa: 1e-4,
        },
    }
}

fn material(width: f64) -> FiniteDuctileJ2<JohnsonCookInitiation> {
    // Explicit numerical mechanical law. The separate original card supplies
    // Cp/rho only; its unrelated affine thermal model is not mechanical authority.
    FiniteDuctileJ2::new(
        native(),
        JohnsonCookInitiation::new(
            JohnsonCookParameters {
                d1: ONSET,
                d2: 0.0,
                d3: 0.0,
                d4: RATE,
                d5: TEMPERATURE,
                reference_strain_rate_s_inv: 1.0,
                transition_temperature_k: 280.0,
                melting_temperature_k: 420.0,
            },
            DuctileDomain {
                triaxiality: [-20.0, 20.0],
                strain_rate_s_inv: [1e-9, 1e6],
                temperature_k: [280.0, 420.0],
            },
        )
        .unwrap(),
        damage_settings(width),
    )
    .unwrap()
}

fn reference(fixture: &Fixture) -> ReferenceCaloricMaterial<'_> {
    ReferenceCaloricMaterial::try_new(
        ReferenceCaloricQuery {
            card: &fixture.card,
            caloric: &fixture.caloric,
            context: source_context(300.0),
            temperature_axis: "temperature".into(),
            caloric_reference: CaloricReference::StressFreeZeroPressure,
        },
        0.0,
        300.0,
        [280.0, 420.0],
    )
    .unwrap()
}

fn sources(
    fixture: &Fixture,
    heat_partition: DuctileHeatPartition,
) -> Vec<FiniteDuctileThermalCell<'_>> {
    vec![FiniteDuctileThermalCell {
        cell: CELL,
        caloric: reference(fixture),
        heat_partition,
    }]
}

fn thermal(fixture: &Fixture, mass_factor: f64, h: f64) -> DepositionSimulation<'_> {
    DepositionSimulation::new(
        vec![TimedEvent {
            time_s: 0.0,
            event: Event::Birth(DepositedCell {
                id: CELL,
                material_id: MaterialId(0),
                segment_id: 1,
                body_id: 7,
                mass_kg: mass_factor * RHO * V,
                specific_h_j_kg: h,
                geometry: CellGeometry {
                    start_m: [0.0; 3],
                    end_m: [L; 3],
                    volume_m3: V,
                    surface_area_m2: 3.0 * L * L,
                },
                boundary: boundary(0.0),
            }),
        }],
        &[(MaterialId(0), fixture.caloric.curve() as &dyn ThermalLaw)],
        DepositionBudget {
            max_cells: 1,
            max_links: 0,
            max_events: 8,
            max_substeps: 1000,
        },
    )
    .unwrap()
    .advance_to(0.0, thermal_config(), || Ok(()))
    .unwrap()
    .0
}

fn increments() -> IncrementSettings {
    IncrementSettings {
        relative_tolerance: 1e-11,
        ..IncrementSettings::default()
    }
}

fn motion_controls() -> FiniteDuctileMotionControls {
    FiniteDuctileMotionControls {
        material: FinitePlasticStepControls {
            equilibrium_work_tolerance_j: 1e-8,
            maximum_endpoint_work_defect_j: 1.0,
            yield_tolerance_pa: 1e-4,
        },
        maximum_numerical_kinetic_loss_j: 1.0,
    }
}

fn prescribed(lambda: f64) -> Vec<(usize, f64)> {
    (0..12)
        .map(|dof| (dof, if dof == 3 { L * (lambda - 1.0) } else { 0.0 }))
        .collect()
}

fn alpha(lambda: f64) -> f64 {
    (2.0 * MU * lambda.ln() - YIELD) / (3.0 * MU)
}
fn active_storage(lambda: f64) -> f64 {
    0.5 * BULK * lambda.ln().powi(2) + YIELD * YIELD / (6.0 * MU)
}
fn integrity(lambda: f64, dt: f64, temperature: f64, width: f64) -> f64 {
    let alpha = alpha(lambda);
    let onset = ONSET
        * (1.0 + RATE * (alpha / dt).ln())
        * (1.0 + TEMPERATURE * (temperature - 280.0) / 140.0);
    (-width * YIELD * (alpha - onset).max(0.0) / GP).exp()
}
fn piola(lambda: f64, integrity: f64) -> f64 {
    integrity * (BULK * lambda.ln() + 2.0 * YIELD / 3.0) / lambda
}

fn adiabatic_first_temperature(lambda: f64, dt: f64, width: f64) -> f64 {
    // Independent scalar conservation root: no calls to a native material trial,
    // coupled solver or thermal step when manufacturing this expected endpoint.
    let mut lo = 300.0;
    let mut hi = 420.0;
    for _ in 0..80 {
        let t = 0.5 * (lo + hi);
        let heat = PARTITION.plastic_dissipation_heat_fraction
            * integrity(lambda, dt, t, width)
            * YIELD
            * alpha(lambda);
        if RHO * LOW_CP * (t - 300.0) < heat {
            lo = t;
        } else {
            hi = t;
        }
    }
    0.5 * (lo + hi)
}

#[derive(Clone, Debug, Default, PartialEq)]
struct Trace {
    accepted_endpoints: Vec<f64>,
    refuse: bool,
}
impl DepositionObserver for Trace {
    type Error = DepositionError;
    fn event(&mut self, _: &Event, _: &DepositionSimulation<'_>) -> Result<(), Self::Error> {
        Ok(())
    }
    fn step(
        &mut self,
        _: f64,
        _: f64,
        _: &BTreeMap<CellId, f64>,
        thermal: &DepositionSimulation<'_>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), Self::Error> {
        checkpoint()?;
        self.accepted_endpoints.push(thermal.time_s());
        if self.refuse {
            Err(DepositionError::Arithmetic)
        } else {
            Ok(())
        }
    }
}

#[test]
fn heated_first_crossing_solves_independent_reciprocal_onset_with_picard_and_iqn() {
    context(|cx| {
        let fixture = Fixture::new(LOW_CP, 81);
        let thermal = thermal(&fixture, 1.0, 0.0);
        let geometry =
            TetLoadGeometry3::new(&NODES, &TETS, TetAssemblyBudget::standard(), 4, cx).unwrap();
        let bands = [TetCrackBandGeometry::from_element(&geometry, 0, [1.0, 0.0, 0.0]).unwrap()];
        let width = bands[0].band_width_m();
        let material = material(width);
        let materials = [&material];
        let native = FiniteDuctileTetSolid::try_new(
            &geometry,
            &materials,
            &bands,
            TetAssemblyBudget::standard(),
            cx,
        )
        .unwrap();
        let body =
            FiniteDuctileThermalSolid::try_new(native, &thermal, sources(&fixture, PARTITION), cx)
                .unwrap();
        let lambda = 1.06;
        let expected_t = adiabatic_first_temperature(lambda, 1.0, width);
        let g = integrity(lambda, 1.0, expected_t, width);
        let plastic = V * g * YIELD * alpha(lambda);
        let stored = V * g * active_storage(lambda);
        let work = V * piola(lambda, g) * (lambda - 1.0);
        let initial_cells = thermal.cells().clone();
        let trace = Trace::default();
        for method in [
            CouplingMethod::RelaxedPicard,
            CouplingMethod::IqnIls(IqnIlsConfig::default()),
        ] {
            let result = body
                .advance_with_coupling(
                    &thermal,
                    &trace,
                    1.0,
                    &prescribed(lambda),
                    &[],
                    [0.0; 3],
                    increments(),
                    motion_controls(),
                    thermal_config(),
                    config(),
                    method,
                    cx,
                )
                .unwrap();
            let report = result.report();
            let energy = &report.material_energy;
            let mechanics = result.solid().mechanics();
            close(
                result.thermal().temperature_k(CELL).unwrap(),
                expected_t,
                1e-7,
            );
            close(
                mechanics.element_states()[0].native_state().alpha,
                alpha(lambda),
                1e-11,
            );
            close(mechanics.element_states()[0].integrity(), g, 1e-9);
            close(
                mechanics.element_first_piola_pa()[0][0],
                piola(lambda, g),
                2e-3,
            );
            close(
                report.equilibrium.reactions_n[3],
                V / L * piola(lambda, g),
                1e-7,
            );
            close(energy.mechanical.stored_after_j, stored, 2e-11);
            close(energy.mechanical.plastic_dissipation_j, plastic, 2e-11);
            close(energy.mechanical.damage_release_j, 0.0, 0.0);
            close(
                energy.mechanical.endpoint_work_defect_j,
                work - stored - plastic,
                2e-11,
            );
            close(
                energy.plastic_heat_j,
                PARTITION.plastic_dissipation_heat_fraction * plastic,
                2e-11,
            );
            close(
                energy.nonthermal_plastic_j,
                (1.0 - PARTITION.plastic_dissipation_heat_fraction) * plastic,
                2e-11,
            );
            close(
                report.applied_heat_j,
                RHO * V * LOW_CP * (expected_t - 300.0),
                5e-11,
            );
            close(report.external_heat_j, 0.0, 1e-12);
            close(report.combined_residual_j, 0.0, config().energy_tolerance_j);
            close(report.thermal_residual_j, 0.0, config().energy_tolerance_j);
            assert!(report.iterations > 1);
            assert_eq!(
                result.observer().accepted_endpoints,
                vec![0.25, 0.5, 0.75, 1.0]
            );
            assert!(std::ptr::eq(
                result.solid().sources()[0].caloric.caloric_source(),
                &fixture.caloric
            ));
            assert_eq!(result.thermal().cells()[&CELL].mass_kg, RHO * V);
        }
        assert!(expected_t > 310.0);
        assert!((g - integrity(lambda, 1.0, 300.0, width)).abs() > 0.004);
        assert!(
            (work - stored - plastic).abs() > PARTITION.plastic_dissipation_heat_fraction * plastic
        );
        assert_eq!(
            body.mechanics().element_states()[0],
            material.initial_state()
        );
        assert_eq!(thermal.cells(), &initial_cells);
        assert_eq!(trace, Trace::default());
        assert_eq!(body.time_s(), 0.0);
    });
}

#[test]
fn later_damage_releases_old_storage_with_explicit_partition_and_no_held_replay() {
    context(|cx| {
        let fixture = Fixture::new(LOW_CP, 82);
        let thermal = thermal(&fixture, 1.0, 0.0);
        let geometry =
            TetLoadGeometry3::new(&NODES, &TETS, TetAssemblyBudget::standard(), 4, cx).unwrap();
        let bands = [TetCrackBandGeometry::from_element(&geometry, 0, [1.0, 0.0, 0.0]).unwrap()];
        let width = bands[0].band_width_m();
        let material = material(width);
        let materials = [&material];
        let native = FiniteDuctileTetSolid::try_new(
            &geometry,
            &materials,
            &bands,
            TetAssemblyBudget::standard(),
            cx,
        )
        .unwrap();
        let body =
            FiniteDuctileThermalSolid::try_new(native, &thermal, sources(&fixture, PARTITION), cx)
                .unwrap();
        let lambda1 = 1.06;
        let first = body
            .advance(
                &thermal,
                &Trace::default(),
                1.0,
                &prescribed(lambda1),
                &[],
                [0.0; 3],
                increments(),
                motion_controls(),
                thermal_config(),
                config(),
                cx,
            )
            .unwrap();
        let (thermal, body, trace, first_report) = first.into_parts();
        let t1 = adiabatic_first_temperature(lambda1, 1.0, width);
        let g1 = integrity(lambda1, 1.0, t1, width);
        let lambda2 = 1.08;
        let delta_alpha = alpha(lambda2) - alpha(lambda1);
        let g2 = g1 * (-width * YIELD * delta_alpha / GP).exp();
        let plastic = V * g2 * YIELD * delta_alpha;
        let damage = V * (g1 - g2) * active_storage(lambda1);
        let heat = PARTITION.plastic_dissipation_heat_fraction * plastic
            + PARTITION.damage_release_heat_fraction * damage;
        let second = body
            .advance(
                &thermal,
                &trace,
                2.0,
                &prescribed(lambda2),
                &[],
                [0.0; 3],
                increments(),
                motion_controls(),
                thermal_config(),
                config(),
                cx,
            )
            .unwrap();
        let energy = &second.report().material_energy;
        assert!(damage > 0.005);
        close(energy.mechanical.plastic_dissipation_j, plastic, 2e-11);
        close(energy.mechanical.damage_release_j, damage, 2e-11);
        close(
            energy.damage_heat_j,
            PARTITION.damage_release_heat_fraction * damage,
            2e-11,
        );
        close(
            energy.nonthermal_damage_j,
            (1.0 - PARTITION.damage_release_heat_fraction) * damage,
            2e-11,
        );
        close(second.report().applied_heat_j, heat, 5e-11);
        close(
            second.thermal().temperature_k(CELL).unwrap(),
            t1 + heat / (RHO * V * LOW_CP),
            1e-7,
        );
        close(
            second.solid().mechanics().element_states()[0].integrity(),
            g2,
            1e-9,
        );
        close(
            second.solid().nonthermal_energy_j()[0].plastic_j,
            first_report.material_energy.nonthermal_plastic_j + energy.nonthermal_plastic_j,
            1e-13,
        );
        close(
            second.solid().nonthermal_energy_j()[0].damage_j,
            energy.nonthermal_damage_j,
            1e-13,
        );
        close(
            second.report().combined_residual_j,
            0.0,
            config().energy_tolerance_j,
        );
        let (thermal, body, trace, _) = second.into_parts();
        let held = body
            .advance(
                &thermal,
                &trace,
                3.0,
                &prescribed(lambda2),
                &[],
                [0.0; 3],
                increments(),
                motion_controls(),
                thermal_config(),
                config(),
                cx,
            )
            .unwrap();
        assert_eq!(
            held.solid().mechanics().element_states(),
            body.mechanics().element_states()
        );
        assert_eq!(
            held.solid().nonthermal_energy_j(),
            body.nonthermal_energy_j()
        );
        close(held.report().applied_heat_j, 0.0, 0.0);
        close(
            held.report()
                .material_energy
                .mechanical
                .plastic_dissipation_j,
            0.0,
            0.0,
        );
        close(
            held.report().material_energy.mechanical.damage_release_j,
            0.0,
            0.0,
        );
        assert_eq!(held.thermal().cells(), thermal.cells());
        assert_eq!(held.observer().accepted_endpoints.len(), 12);
        assert_eq!(trace.accepted_endpoints.len(), 8);
    });
}

#[test]
fn actual_free_coordinate_solves_rate_temperature_damage_and_original_consistent_inertia() {
    context(|cx| {
        let fixture = Fixture::new(LOW_CP, 83);
        let thermal = thermal(&fixture, 1.0, 0.0);
        let geometry =
            TetLoadGeometry3::new(&NODES, &TETS, TetAssemblyBudget::standard(), 4, cx).unwrap();
        let bands = [TetCrackBandGeometry::from_element(&geometry, 0, [1.0, 0.0, 0.0]).unwrap()];
        let width = bands[0].band_width_m();
        let material = material(width);
        let materials = [&material];
        let supports: Vec<_> = (0..12)
            .filter(|dof| *dof != 3)
            .map(|dof| (dof, 0.0))
            .collect();
        let lambda = 1.05;
        let displacement = L * (lambda - 1.0);
        let effective_mass = RHO * V / 10.0;
        let mut damage = Vec::new();
        for dt in [1e-5, 1.5e-5] {
            let native = FiniteDuctileTetSolid::try_new(
                &geometry,
                &materials,
                &bands,
                TetAssemblyBudget::standard(),
                cx,
            )
            .unwrap();
            let body = FiniteDuctileThermalSolid::try_new(
                native,
                &thermal,
                sources(&fixture, PARTITION),
                cx,
            )
            .unwrap()
            .with_inertia(&[[0.0; 3]; 4], cx)
            .unwrap();
            let expected_t = adiabatic_first_temperature(lambda, dt, width);
            let g = integrity(lambda, dt, expected_t, width);
            let acceleration = displacement / (dt * dt);
            let force = V / L * piola(lambda, g) + effective_mass * acceleration;
            let result = body
                .advance(
                    &thermal,
                    &(),
                    dt,
                    &supports,
                    &[(3, force)],
                    [0.0; 3],
                    increments(),
                    motion_controls(),
                    thermal_config(),
                    config(),
                    cx,
                )
                .unwrap();
            let report = result.report();
            let inertia = report.inertia.as_ref().unwrap();
            let kinetic = 0.5 * effective_mass * (displacement / dt).powi(2);
            let plastic = V * g * YIELD * alpha(lambda);
            close(
                result.thermal().temperature_k(CELL).unwrap(),
                expected_t,
                1e-7,
            );
            close(
                result.solid().mechanics().displacement_m()[1][0],
                displacement,
                1e-11,
            );
            close(
                result.solid().velocity_mps().unwrap()[1][0],
                displacement / dt,
                1e-6,
            );
            close(
                result.solid().mechanics().element_states()[0].integrity(),
                g,
                1e-9,
            );
            close(result.solid().reference_mass_kg().unwrap(), RHO * V, 1e-18);
            close(inertia.reference_mass_kg, RHO * V, 1e-18);
            close(inertia.kinetic_change_j, kinetic, 1e-9);
            close(inertia.numerical_loss_j, kinetic, 1e-9);
            close(inertia.inertial_work_j, 2.0 * kinetic, 1e-9);
            close(
                report.material_energy.mechanical.dead_load_endpoint_work_j,
                force * displacement,
                1e-9,
            );
            close(
                report.material_energy.mechanical.support_endpoint_work_j,
                0.0,
                0.0,
            );
            close(
                report.applied_heat_j,
                PARTITION.plastic_dissipation_heat_fraction * plastic,
                5e-11,
            );
            assert!(inertia.numerical_loss_j > report.applied_heat_j);
            assert!(report.equilibrium.residual_history.len() > 1);
            close(report.combined_residual_j, 0.0, config().energy_tolerance_j);
            let reaction = supports
                .iter()
                .zip(&report.equilibrium.reactions_n)
                .find(|((dof, _), _)| *dof == 6)
                .unwrap()
                .1;
            close(*reaction, RHO * V / 20.0 * acceleration, 1e-6);
            damage.push(result.solid().mechanics().element_states()[0].damage());
            let (thermal, body, (), _) = result.into_parts();
            // Arresting velocity through fixed supports is numerical kinetic
            // loss with zero displacement work, not another material heat pulse.
            let held_supports: Vec<_> = body
                .mechanics()
                .displacement_m()
                .iter()
                .flatten()
                .copied()
                .enumerate()
                .collect();
            let stopped = body
                .advance(
                    &thermal,
                    &(),
                    2.0 * dt,
                    &held_supports,
                    &[],
                    [0.0; 3],
                    increments(),
                    motion_controls(),
                    thermal_config(),
                    config(),
                    cx,
                )
                .unwrap();
            close(stopped.report().applied_heat_j, 0.0, 1e-13);
            close(
                stopped.report().inertia.as_ref().unwrap().numerical_loss_j,
                kinetic,
                1e-9,
            );
            close(
                stopped.report().inertia.as_ref().unwrap().kinetic_change_j,
                -kinetic,
                1e-9,
            );
            close(
                stopped.report().combined_residual_j,
                0.0,
                config().energy_tolerance_j,
            );
            assert_eq!(
                stopped.solid().nonthermal_energy_j(),
                body.nonthermal_energy_j()
            );
            close(
                stopped.thermal().temperature_k(CELL).unwrap(),
                expected_t,
                1e-7,
            );
        }
        assert!((damage[0] - damage[1]).abs() > 0.0005);
    });
}

#[test]
fn refused_candidates_preserve_original_source_caloric_mass_history_and_observer() {
    context(|cx| {
        let fixture = Fixture::new(LOW_CP, 84);
        let thermal = thermal(&fixture, 1.0, 0.0);
        let geometry =
            TetLoadGeometry3::new(&NODES, &TETS, TetAssemblyBudget::standard(), 4, cx).unwrap();
        let bands = [TetCrackBandGeometry::from_element(&geometry, 0, [1.0, 0.0, 0.0]).unwrap()];
        let material = material(bands[0].band_width_m());
        let materials = [&material];
        let make_native = || {
            FiniteDuctileTetSolid::try_new(
                &geometry,
                &materials,
                &bands,
                TetAssemblyBudget::standard(),
                cx,
            )
            .unwrap()
        };
        let body = FiniteDuctileThermalSolid::try_new(
            make_native(),
            &thermal,
            sources(&fixture, PARTITION),
            cx,
        )
        .unwrap();
        let original_cells = thermal.cells().clone();
        let trace = Trace::default();
        let request = |simulation: &DepositionSimulation<'_>, observer: &Trace, policy| {
            body.advance(
                simulation,
                observer,
                1.0,
                &prescribed(1.06),
                &[],
                [0.0; 3],
                increments(),
                motion_controls(),
                thermal_config(),
                policy,
                cx,
            )
            .map(|candidate| candidate.report().applied_heat_j)
        };
        let one_trial = FiniteThermoplasticConfig {
            maximum_iterations: 1,
            ..config()
        };
        assert!(matches!(
            request(&thermal, &trace, one_trial),
            Err(FiniteThermoplasticError::NotConverged { .. })
        ));
        let refusing = Trace {
            refuse: true,
            ..Trace::default()
        };
        assert!(matches!(
            request(&thermal, &refusing, config()),
            Err(FiniteThermoplasticError::Observer(
                DepositionError::Arithmetic
            ))
        ));
        assert!(refusing.accepted_endpoints.is_empty());
        let replaced_enthalpy = self::thermal(&fixture, 1.0, 0.5);
        assert!(matches!(
            request(&replaced_enthalpy, &trace, config()),
            Err(FiniteThermoplasticError::Invalid(_))
        ));
        let wrong_mass = self::thermal(&fixture, 2.0, 0.0);
        assert!(matches!(
            FiniteDuctileThermalSolid::try_new(
                make_native(),
                &wrong_mass,
                sources(&fixture, PARTITION),
                cx
            ),
            Err(FiniteThermoplasticError::Invalid(_))
        ));
        let duplicate_chart = Fixture::new(LOW_CP, 84);
        assert_eq!(
            duplicate_chart.card.content_hash(),
            fixture.card.content_hash()
        );
        assert!(matches!(
            FiniteDuctileThermalSolid::try_new(
                make_native(),
                &thermal,
                sources(&duplicate_chart, PARTITION),
                cx
            ),
            Err(FiniteThermoplasticError::Invalid(_))
        ));
        let unaccounted = make_native()
            .advance(
                &prescribed(1.06),
                &[],
                [0.0; 3],
                DuctileStepConditions {
                    duration_s: 1.0,
                    temperatures_k: DuctileTemperatures::Uniform(300.0),
                },
                increments(),
                motion_controls().material,
                cx,
            )
            .unwrap();
        assert!(matches!(
            FiniteDuctileThermalSolid::try_new(
                unaccounted.solid,
                &thermal,
                sources(&fixture, PARTITION),
                cx
            ),
            Err(FiniteThermoplasticError::Invalid(_))
        ));
        let accepted_heat = request(&thermal, &trace, config()).unwrap();
        close(
            accepted_heat,
            request(&thermal, &trace, config()).unwrap(),
            0.0,
        );
        assert!(accepted_heat > 0.0);
        assert_eq!(
            body.mechanics().element_states()[0],
            material.initial_state()
        );
        assert_eq!(body.nonthermal_energy_j()[0].plastic_j, 0.0);
        assert_eq!(body.nonthermal_energy_j()[0].damage_j, 0.0);
        assert_eq!(thermal.cells(), &original_cells);
        assert_eq!(trace, Trace::default());
        assert_eq!(body.time_s(), 0.0);
    });
}

fn sourced_fixture(identity: u32) -> Fixture {
    let mut claims = ClaimSet::new();
    claims
        .insert_claim(scalar(CP_NAME, CP_DIMS, LOW_CP))
        .unwrap();
    claims
        .insert_claim(scalar(DENSITY_PROPERTY, Density::DIMS, RHO))
        .unwrap();
    claims
        .insert_claim(PropertyClaim {
            key: PropertyKey::new(FRACTURE_PLASTIC_STRAIN_PROPERTY, Dims::NONE),
            value: PropertyValue::Curve {
                abscissa: STRESS_TRIAXIALITY_AXIS.into(),
                abscissa_dims: Dims::NONE,
                knots: vec![(-20.0, ONSET), (20.0, ONSET)],
                dims: Dims::NONE,
            },
            validity: ValidityDomain::unconstrained()
                .with_quantity("temperature", temperature_quantity(), 280.0, 420.0)
                .with_quantity(
                    STRESS_TRIAXIALITY_AXIS,
                    QuantitySpec::dimensional(Dims::NONE),
                    -20.0,
                    20.0,
                )
                .with_quantity(
                    PLASTIC_STRAIN_RATE_AXIS,
                    QuantitySpec::dimensional(Dims([0, 0, -1, 0, 0, 0])),
                    1e-9,
                    1e6,
                ),
            interpolation: InterpolationPolicy::TabulatedOnly,
            uncertainty: UncertaintyModel::Unstated,
            observations: vec![],
            provenance: provenance(),
        })
        .unwrap();
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: format!(
                "synthetic same-card caloric and fixed-temperature ductility {identity}"
            ),
            phase: "solid".into(),
            process: "numerical regression".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    let caloric = ResolvedHeatCapacityCurve::try_new(
        &card,
        HeatCapacitySourceQuery {
            heat_capacity: &PropertyKey::new(CP_NAME, CP_DIMS),
            density: &PropertyKey::new(DENSITY_PROPERTY, Density::DIMS),
            temperature_axis: "temperature",
            lower: &source_context(280.0),
            upper: &source_context(420.0),
            selection: SelectionPolicy::SingleClaimOnly,
        },
        SolidLiquidPhase::Solid,
        -20.0 * LOW_CP,
        1e-6,
    )
    .unwrap();
    Fixture { card, caloric }
}

#[test]
fn original_fixed_temperature_source_refuses_heating_foreign_cards_and_cancellation() {
    let gate = CancelGate::new_clock_free();
    fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default()).scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 3,
                kernel_id: 0x44554354,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        let fixture = sourced_fixture(85);
        let thermal = thermal(&fixture, 1.0, 0.0);
        let geometry =
            TetLoadGeometry3::new(&NODES, &TETS, TetAssemblyBudget::standard(), 4, &cx).unwrap();
        let bands = [TetCrackBandGeometry::from_element(&geometry, 0, [1.0, 0.0, 0.0]).unwrap()];
        let point = source_context(300.0)
            .with_quantity(
                STRESS_TRIAXIALITY_AXIS,
                QuantitySpec::dimensional(Dims::NONE),
                -20.0,
            )
            .unwrap()
            .with_quantity(
                PLASTIC_STRAIN_RATE_AXIS,
                QuantitySpec::dimensional(Dims([0, 0, -1, 0, 0, 0])),
                0.01,
            )
            .unwrap();
        let locus = resolve_ductile_locus(
            &fixture.card,
            &point,
            &DuctileSourceRequest {
                temperature_axis: "temperature".into(),
                maximum_points: 2,
            },
            fs_matdb::ClaimSelection::Policy(SelectionPolicy::SingleClaimOnly),
        )
        .unwrap();
        let material =
            FiniteDuctileJ2::new(native(), locus, damage_settings(bands[0].band_width_m()))
                .unwrap();
        let materials = [&material];
        let make_native = || {
            FiniteDuctileTetSolid::try_new(
                &geometry,
                &materials,
                &bands,
                TetAssemblyBudget::standard(),
                &cx,
            )
            .unwrap()
        };
        let retained = DuctileHeatPartition {
            plastic_dissipation_heat_fraction: 0.0,
            damage_release_heat_fraction: 0.0,
        };
        let body = FiniteDuctileThermalSolid::try_new(
            make_native(),
            &thermal,
            sources(&fixture, retained),
            &cx,
        )
        .unwrap();
        let accepted = body
            .advance(
                &thermal,
                &(),
                1.0,
                &prescribed(1.06),
                &[],
                [0.0; 3],
                increments(),
                motion_controls(),
                thermal_config(),
                config(),
                &cx,
            )
            .unwrap();
        assert!(accepted.solid().mechanics().element_states()[0].damage() > 0.1);
        close(accepted.thermal().temperature_k(CELL).unwrap(), 300.0, 0.0);
        close(accepted.report().applied_heat_j, 0.0, 0.0);
        close(
            accepted.solid().nonthermal_energy_j()[0].plastic_j,
            accepted
                .report()
                .material_energy
                .mechanical
                .plastic_dissipation_j,
            0.0,
        );
        close(
            accepted.report().combined_residual_j,
            0.0,
            config().energy_tolerance_j,
        );
        let heated = FiniteDuctileThermalSolid::try_new(
            make_native(),
            &thermal,
            sources(&fixture, PARTITION),
            &cx,
        )
        .unwrap();
        assert!(matches!(
            heated.advance(
                &thermal,
                &(),
                1.0,
                &prescribed(1.06),
                &[],
                [0.0; 3],
                increments(),
                motion_controls(),
                thermal_config(),
                config(),
                &cx
            ),
            Err(FiniteThermoplasticError::DuctileSolid(_))
        ));
        let foreign = sourced_fixture(86);
        let foreign_thermal = self::thermal(&foreign, 1.0, 0.0);
        // Chart, inventory and numerical Cp/rho all match their own source;
        // refusal is the original ductile calibration's distinct card identity.
        assert!(matches!(
            FiniteDuctileThermalSolid::try_new(
                make_native(),
                &foreign_thermal,
                sources(&foreign, retained),
                &cx
            ),
            Err(FiniteThermoplasticError::Invalid(_))
        ));
        gate.request();
        assert!(matches!(
            body.advance(
                &thermal,
                &(),
                1.0,
                &prescribed(1.06),
                &[],
                [0.0; 3],
                increments(),
                motion_controls(),
                thermal_config(),
                config(),
                &cx
            ),
            Err(FiniteThermoplasticError::Cancelled)
        ));
        assert_eq!(
            body.mechanics().element_states()[0],
            material.initial_state()
        );
        assert_eq!(
            heated.mechanics().element_states()[0],
            material.initial_state()
        );
        assert_eq!(thermal.cells()[&CELL].specific_h_j_kg, 0.0);
        assert_eq!(thermal.time_s(), 0.0);
    });
}
