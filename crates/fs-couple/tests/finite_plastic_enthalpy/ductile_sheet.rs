//! G1/G3: progressive sheet heat, reciprocal initiation, native quadrature and
//! retained membrane motion. Synthetic low heat capacities magnify feedback;
//! these fixtures do not qualify an alloy or measured heat partition.

use super::*;
use fs_couple::solid_enthalpy::finite_ductile::sheet::FiniteDuctileThermalSheet;
use fs_couple::solid_enthalpy::finite_ductile::{
    DuctileHeatPartition, FiniteDuctileThermalCell, ReferenceCaloricMaterial, ReferenceCaloricQuery,
};
use fs_material::ductile::evolution::DuctileJ2Settings;
use fs_material::ductile::{DuctileDomain, JohnsonCookInitiation, JohnsonCookParameters};
use fs_material::hyper::plane_stress::PlaneStressSettings;
use fs_material::plastic::energy::EnergyTolerances;
use fs_material::plastic::finite::ductile::FiniteDuctileJ2;
use fs_material::plastic::finite::rate::FiniteOverstressParameters;
use fs_material::plastic::finite::{FiniteJ2, FiniteJ2Limits, FiniteJ2Parameters};
use fs_solid::ductile::{DuctileStepConditions, DuctileTemperatures};
use fs_solid::hyper2d::plastic::SheetAssemblyBudget;
use fs_solid::hyper2d::plastic::ductile::dynamics::FiniteDuctileSheetMotionControls;
use fs_solid::hyper2d::plastic::ductile::{FiniteDuctileSheet, SheetCrackBandGeometry};
use fs_solid::{IncrementSettings, Mesh2};

const MU: f64 = 80e6;
const BULK: f64 = 160e6;
const YIELD: f64 = 2e6;
const GP: f64 = 400.0;
const LOW_CP: f64 = 2.0;
const ONSET: f64 = 0.008;
const RATE: f64 = 0.04;
const TEMPERATURE: f64 = -0.6;
const PARTITION: DuctileHeatPartition = DuctileHeatPartition {
    plastic_dissipation_heat_fraction: 0.75,
    damage_release_heat_fraction: 0.25,
};
const ZERO_HEAT: DuctileHeatPartition = DuctileHeatPartition {
    plastic_dissipation_heat_fraction: 0.0,
    damage_release_heat_fraction: 0.0,
};
const C: f64 = 1.0 / (6.0 * MU) + 2.0 / (9.0 * BULK);

fn material(width: f64, viscosity: f64) -> FiniteDuctileJ2<JohnsonCookInitiation> {
    let native = FiniteJ2::new(
        FiniteJ2Parameters {
            shear_modulus_pa: MU,
            bulk_modulus_pa: BULK,
            yield_stress_pa: YIELD,
            hardening_modulus_pa: 0.0,
        },
        FiniteJ2Limits {
            minimum_jacobian: 0.2,
            maximum_elastic_log_strain: 0.5,
            maximum_plastic_log_strain: 2.0,
            maximum_accumulated_plastic_strain: 3.0,
        },
    )
    .unwrap();
    let initiation = JohnsonCookInitiation::new(
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
            triaxiality: [-50.0, 50.0],
            strain_rate_s_inv: [1e-9, 1e6],
            temperature_k: [280.0, 420.0],
        },
    )
    .unwrap();
    let material = FiniteDuctileJ2::new(
        native,
        initiation,
        DuctileJ2Settings {
            characteristic_length_m: width,
            post_initiation_plastic_work_j_m2: GP,
            maximum_damage_increment: 1.0,
            energy_tolerances: EnergyTolerances {
                energy_j_m3: 1e-6,
                yield_pa: 1e-4,
            },
        },
    )
    .unwrap();
    if viscosity == 0.0 {
        material
    } else {
        material
            .with_overstress(FiniteOverstressParameters {
                viscosity_pa_s: viscosity,
                stress_scale_pa: YIELD,
                exponent: 1.0,
            })
            .unwrap()
    }
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
fn sources<'a>(
    fixture: &'a Fixture,
    partitions: &[DuctileHeatPartition],
) -> Vec<FiniteDuctileThermalCell<'a>> {
    partitions
        .iter()
        .enumerate()
        .map(|(i, &heat_partition)| FiniteDuctileThermalCell {
            cell: CellId(CELL.0 + i as u64),
            caloric: reference(fixture),
            heat_partition,
        })
        .collect()
}
fn thermal<'a>(fixture: &'a Fixture, volumes: &[f64]) -> DepositionSimulation<'a> {
    let events = volumes
        .iter()
        .enumerate()
        .map(|(i, &volume)| TimedEvent {
            time_s: 0.0,
            event: Event::Birth(DepositedCell {
                id: CellId(CELL.0 + i as u64),
                material_id: MaterialId(0),
                segment_id: i as u64,
                body_id: 7,
                mass_kg: RHO * volume,
                specific_h_j_kg: 0.0,
                geometry: CellGeometry {
                    start_m: [i as f64 * L, 0.0, 0.0],
                    end_m: [(i + 1) as f64 * L, L, 0.001],
                    volume_m3: volume,
                    surface_area_m2: 2.0 * L * L,
                },
                boundary: boundary(0.0),
            }),
        })
        .collect();
    DepositionSimulation::new(
        events,
        &[(MaterialId(0), fixture.caloric.curve() as &dyn ThermalLaw)],
        DepositionBudget {
            max_cells: volumes.len(),
            max_links: 0,
            max_events: volumes.len() + 2,
            max_substeps: 1000,
        },
    )
    .unwrap()
    .advance_to(0.0, thermal_config(), || Ok(()))
    .unwrap()
    .0
}
fn bands(mesh: &Mesh2, cx: &Cx<'_>) -> Vec<SheetCrackBandGeometry> {
    (0..mesh.elems.len())
        .map(|e| SheetCrackBandGeometry::from_element(mesh, e, [1.0, 0.0], cx).unwrap())
        .collect()
}
fn sheet<'g, 'm>(
    mesh: &'g Mesh2,
    models: &'m [&'m FiniteDuctileJ2<JohnsonCookInitiation>],
    thickness: &'g [f64],
    bands: &'g [SheetCrackBandGeometry],
    cx: &Cx<'_>,
) -> FiniteDuctileSheet<'g, 'm, JohnsonCookInitiation> {
    FiniteDuctileSheet::try_new(
        mesh,
        models,
        thickness,
        bands,
        PlaneStressSettings {
            relative_tolerance: 1e-12,
            ..PlaneStressSettings::default()
        },
        SheetAssemblyBudget {
            maximum_nodes: 20,
            maximum_elements: 8,
            maximum_free_dofs: 40,
            minimum_scaled_jacobian: 1e-4,
        },
        cx,
    )
    .unwrap()
}
fn motion_controls() -> FiniteDuctileSheetMotionControls {
    FiniteDuctileSheetMotionControls {
        material: FinitePlasticStepControls {
            equilibrium_work_tolerance_j: 1e-8,
            maximum_endpoint_work_defect_j: 1.0,
            yield_tolerance_pa: 1e-4,
        },
        maximum_numerical_kinetic_loss_j: 1.0,
    }
}
fn increments() -> IncrementSettings {
    IncrementSettings {
        relative_tolerance: 1e-11,
        ..IncrementSettings::default()
    }
}
fn supports(mesh: &Mesh2, stretch: f64) -> Vec<(usize, f64)> {
    mesh.nodes
        .iter()
        .enumerate()
        .flat_map(|(i, p)| {
            [
                (2 * i, (stretch - 1.0) * p[0]),
                (2 * i + 1, (stretch - 1.0) * p[1]),
            ]
        })
        .collect()
}
fn alpha(stretch: f64) -> f64 {
    2.0 * (stretch.ln() - C * YIELD)
}
fn integrity(stretch: f64, dt: f64, temperature: f64, width: f64) -> f64 {
    let alpha = alpha(stretch);
    let onset = ONSET
        * (1.0 + RATE * (alpha / dt).ln())
        * (1.0 + TEMPERATURE * (temperature - 280.0) / 140.0);
    (-width * YIELD * (alpha - onset).max(0.0) / GP).exp()
}
fn adiabatic_temperature(stretch: f64, dt: f64, width: f64) -> f64 {
    // Independent plane-stress Hencky equations and scalar caloric balance.
    // No native material return or feedback solve manufactures this reference.
    let (mut low, mut high) = (300.0, 420.0);
    for _ in 0..80 {
        let temperature = 0.5 * (low + high);
        let heat = PARTITION.plastic_dissipation_heat_fraction
            * integrity(stretch, dt, temperature, width)
            * YIELD
            * alpha(stretch);
        if RHO * LOW_CP * (temperature - 300.0) < heat {
            low = temperature;
        } else {
            high = temperature;
        }
    }
    0.5 * (low + high)
}

#[derive(Clone, Debug, Default, PartialEq)]
struct Trace {
    endpoints: Vec<f64>,
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
        self.endpoints.push(thermal.time_s());
        if self.refuse {
            Err(DepositionError::Arithmetic)
        } else {
            Ok(())
        }
    }
}

#[test]
fn independent_reciprocal_sheet_onset_heat_and_thickness_use_one_history_for_picard_and_iqn() {
    context(|cx| {
        let fixture = Fixture::new(LOW_CP, 201);
        let stretch = 1.025;
        for mesh in [
            Mesh2 {
                nodes: vec![[0.0, 0.0], [L, 0.0], [0.0, L]],
                elems: vec![vec![0, 1, 2]],
                patches: vec![],
            },
            Mesh2::quads(L, L, 1, 1),
        ] {
            let bands = bands(&mesh, cx);
            let model = material(bands[0].band_width_m(), 0.0);
            let models = [&model];
            let thickness = [0.001];
            let native = sheet(&mesh, &models, &thickness, &bands, cx);
            let volume = native.reference_volume_m3(0).unwrap();
            let thermal = thermal(&fixture, &[volume]);
            let original_cells = thermal.cells().clone();
            let body = FiniteDuctileThermalSheet::try_new(
                native,
                &thermal,
                sources(&fixture, &[PARTITION]),
                cx,
            )
            .unwrap();
            let t = adiabatic_temperature(stretch, 1.0, bands[0].band_width_m());
            let g = integrity(stretch, 1.0, t, bands[0].band_width_m());
            let plastic = volume * g * YIELD * alpha(stretch);
            let storage = volume * g * C * YIELD * YIELD;
            for method in [
                CouplingMethod::RelaxedPicard,
                CouplingMethod::IqnIls(IqnIlsConfig::default()),
            ] {
                let result = body
                    .advance_with_coupling(
                        &thermal,
                        &Trace::default(),
                        1.0,
                        &supports(&mesh, stretch),
                        &[],
                        increments(),
                        motion_controls(),
                        thermal_config(),
                        config(),
                        method,
                        cx,
                    )
                    .unwrap();
                let energy = &result.report().material_energy;
                close(result.thermal().temperature_k(CELL).unwrap(), t, 1e-7);
                assert_eq!(
                    result.sheet().mechanics().quadrature_states().len(),
                    if mesh.elems[0].len() == 3 { 1 } else { 4 }
                );
                for state in result.sheet().mechanics().quadrature_states() {
                    close(state.native_state().alpha, alpha(stretch), 2e-10);
                    close(state.integrity(), g, 2e-9);
                    close(state.first_piola_pa()[0], g * YIELD / stretch, 2e-3);
                    close(
                        state.deformation()[8],
                        (-2.0 * stretch.ln() + 2.0 * YIELD / (3.0 * BULK)).exp(),
                        2e-10,
                    );
                }
                close(energy.mechanical.plastic_dissipation_j, plastic, 2e-10);
                close(energy.mechanical.stored_after_j, storage, 2e-10);
                assert_eq!(energy.mechanical.damage_release_j, 0.0);
                close(
                    energy.plastic_heat_j,
                    PARTITION.plastic_dissipation_heat_fraction * plastic,
                    2e-10,
                );
                close(
                    energy.nonthermal_plastic_j,
                    (1.0 - PARTITION.plastic_dissipation_heat_fraction) * plastic,
                    2e-10,
                );
                close(
                    result.report().applied_heat_j,
                    RHO * volume * LOW_CP * (t - 300.0),
                    2e-10,
                );
                close(
                    result.report().combined_residual_j,
                    0.0,
                    config().energy_tolerance_j,
                );
                assert!(result.report().iterations > 1);
                assert_eq!(result.observer().endpoints, vec![0.25, 0.5, 0.75, 1.0]);
                assert!(std::ptr::eq(
                    result.sheet().sources()[0].caloric.caloric_source(),
                    &fixture.caloric
                ));
                assert_eq!(result.thermal().cells()[&CELL].mass_kg, RHO * volume);
            }
            assert!(t > 305.0);
            assert!((g - integrity(stretch, 1.0, 300.0, bands[0].band_width_m())).abs() > 1e-4);
            assert!(
                body.mechanics()
                    .quadrature_states()
                    .iter()
                    .all(|s| s.native_state().alpha == 0.0 && s.damage() == 0.0)
            );
            assert_eq!(thermal.cells(), &original_cells);
            assert_eq!(body.time_s(), 0.0);
        }
    });
}

#[test]
fn warped_quad_cells_keep_pointwise_rate_heat_and_nonthermal_inventory_without_heating_from_inertia()
 {
    context(|cx| {
        let fixture = Fixture::new(LOW_CP, 202);
        let mut mesh = Mesh2::quads(2.0 * L, L, 2, 1);
        let middle = mesh.nodes.iter().position(|p| *p == [L, L]).unwrap();
        mesh.nodes[middle] = [0.9 * L, 1.1 * L];
        let bands = bands(&mesh, cx);
        let models: Vec<_> = bands
            .iter()
            .map(|b| material(b.band_width_m(), 20e6))
            .collect();
        let refs: Vec<_> = models.iter().collect();
        let thickness = [0.001, 0.0025];
        let native = sheet(&mesh, &refs, &thickness, &bands, cx);
        let volumes: Vec<_> = (0..2)
            .map(|e| native.reference_volume_m3(e).unwrap())
            .collect();
        let thermal = thermal(&fixture, &volumes);
        let partitions = [ZERO_HEAT, PARTITION];
        let body = FiniteDuctileThermalSheet::try_new(
            native,
            &thermal,
            sources(&fixture, &partitions),
            cx,
        )
        .unwrap()
        .with_inertia(&vec![[0.2, -0.1]; mesh.nodes.len()], 1 << 20, cx)
        .unwrap();
        let target: Vec<_> = mesh
            .nodes
            .iter()
            .map(|&[x, y]| [0.025 * x + 0.015 * x * y / L, 0.015 * y])
            .collect();
        let prescribed: Vec<_> = target.as_flattened().iter().copied().enumerate().collect();
        let result = body
            .advance(
                &thermal,
                &Trace::default(),
                0.5,
                &prescribed,
                &[],
                increments(),
                motion_controls(),
                thermal_config(),
                config(),
                cx,
            )
            .unwrap();
        assert_eq!(result.sheet().mechanics().quadrature_states().len(), 8);
        let mut point_temperatures = Vec::new();
        for e in 0..2 {
            point_temperatures.extend(std::iter::repeat_n(
                result
                    .thermal()
                    .temperature_k(CellId(CELL.0 + e as u64))
                    .unwrap(),
                4,
            ));
        }
        let reference = body
            .mechanics()
            .advance(
                &prescribed,
                &[],
                DuctileStepConditions {
                    duration_s: 0.5,
                    temperatures_k: DuctileTemperatures::PerPoint(&point_temperatures),
                },
                increments(),
                motion_controls().material,
                cx,
            )
            .unwrap();
        for e in 0..2 {
            let mut plastic = 0.0;
            let mut damage = 0.0;
            for q in 0..4 {
                let volume = body.mechanics().quadrature_volume_m3(e, q).unwrap();
                let account = reference.point_energy[4 * e + q];
                plastic += volume * account.plastic_dissipation_j_m3;
                damage += volume * account.damage_release_j_m3;
                let accepted = &result.sheet().mechanics().quadrature_states()[4 * e + q];
                let expected = &reference.sheet.quadrature_states()[4 * e + q];
                close(
                    accepted.native_state().alpha,
                    expected.native_state().alpha,
                    2e-10,
                );
                close(accepted.damage(), expected.damage(), 2e-9);
            }
            let cell = result.report().material_energy.cells[e];
            close(
                cell.plastic_heat_j,
                partitions[e].plastic_dissipation_heat_fraction * plastic,
                2e-10,
            );
            close(
                cell.damage_heat_j,
                partitions[e].damage_release_heat_fraction * damage,
                2e-10,
            );
            close(
                cell.nonthermal_plastic_j,
                (1.0 - partitions[e].plastic_dissipation_heat_fraction) * plastic,
                2e-10,
            );
            assert_eq!(
                result.thermal().cells()[&CellId(CELL.0 + e as u64)].mass_kg,
                RHO * volumes[e]
            );
        }
        let v0 = body.mechanics().quadrature_volume_m3(0, 0).unwrap();
        assert!(
            (1..4)
                .any(|q| (body.mechanics().quadrature_volume_m3(0, q).unwrap() - v0).abs() > 1e-12)
        );
        let states = result.sheet().mechanics().quadrature_states();
        assert!((1..4).any(|q| {
            (states[q].native_state().alpha - states[0].native_state().alpha).abs() > 1e-5
        }));
        close(result.thermal().temperature_k(CELL).unwrap(), 300.0, 1e-10);
        assert!(result.thermal().temperature_k(CellId(CELL.0 + 1)).unwrap() > 300.1);
        assert!(result.report().inertia.as_ref().unwrap().numerical_loss_j > 1e-7);
        assert!(result.sheet().nonthermal_energy_j()[0].plastic_j > 0.0);
        close(
            result.report().combined_residual_j,
            0.0,
            config().energy_tolerance_j,
        );
        let (thermal, body, observer, first_report) = result.into_parts();
        let first_alpha: Vec<_> = body
            .mechanics()
            .quadrature_states()
            .iter()
            .map(|s| s.native_state().alpha)
            .collect();
        let first_temperature = thermal.temperature_k(CellId(CELL.0 + 1)).unwrap();
        let held = body
            .advance(
                &thermal,
                &observer,
                1.0,
                &prescribed,
                &[],
                increments(),
                motion_controls(),
                thermal_config(),
                config(),
                cx,
            )
            .unwrap();
        assert!(
            held.sheet()
                .mechanics()
                .quadrature_states()
                .iter()
                .zip(first_alpha)
                .all(|(s, old)| s.native_state().alpha > old)
        );
        assert!(held.report().material_energy.mechanical.damage_release_j > 0.0);
        assert!(held.report().material_energy.damage_heat_j > 0.0);
        close(held.thermal().temperature_k(CELL).unwrap(), 300.0, 1e-10);
        assert!(held.thermal().temperature_k(CellId(CELL.0 + 1)).unwrap() > first_temperature);
        assert_eq!(held.report().inertia.as_ref().unwrap().kinetic_after_j, 0.0);
        assert_eq!(held.report().inertia.as_ref().unwrap().inertial_work_j, 0.0);
        assert_eq!(held.sheet().reference_mass_kg(), body.reference_mass_kg());
        for e in 0..2 {
            close(
                held.sheet().nonthermal_energy_j()[e].plastic_j,
                first_report.material_energy.cells[e].nonthermal_plastic_j
                    + held.report().material_energy.cells[e].nonthermal_plastic_j,
                2e-12,
            );
            close(
                held.sheet().nonthermal_energy_j()[e].damage_j,
                first_report.material_energy.cells[e].nonthermal_damage_j
                    + held.report().material_energy.cells[e].nonthermal_damage_j,
                2e-12,
            );
        }
        assert_eq!(held.observer().endpoints, vec![0.25, 0.5, 0.75, 1.0]);
        assert_eq!(body.time_s(), 0.5);
    });
}

#[test]
fn source_clock_loss_history_and_failed_feedback_leave_original_sheet_thermal_and_observer_unchanged()
 {
    context(|cx| {
        let fixture = Fixture::new(LOW_CP, 203);
        let other = Fixture::new(LOW_CP, 204);
        let mesh = Mesh2::quads(L, L, 1, 1);
        let bands = bands(&mesh, cx);
        let model = material(bands[0].band_width_m(), 0.0);
        let models = [&model];
        let thickness = [0.001];
        let native = sheet(&mesh, &models, &thickness, &bands, cx);
        let volume = native.reference_volume_m3(0).unwrap();
        let thermal = thermal(&fixture, &[volume]);
        let original_cells = thermal.cells().clone();
        assert!(
            FiniteDuctileThermalSheet::try_new(
                native.clone(),
                &thermal,
                sources(&other, &[PARTITION]),
                cx
            )
            .is_err()
        );
        assert!(FiniteDuctileThermalSheet::try_new(native.clone(), &thermal, vec![], cx).is_err());
        let later = thermal
            .advance_to(0.5, thermal_config(), || Ok(()))
            .unwrap()
            .0;
        assert!(
            FiniteDuctileThermalSheet::try_new(
                native.clone(),
                &later,
                sources(&fixture, &[PARTITION]),
                cx
            )
            .is_err()
        );
        let prescribed = supports(&mesh, 1.025);
        let loaded = native
            .advance(
                &prescribed,
                &[],
                DuctileStepConditions {
                    duration_s: 0.5,
                    temperatures_k: DuctileTemperatures::Uniform(300.0),
                },
                increments(),
                motion_controls().material,
                cx,
            )
            .unwrap();
        assert!(
            FiniteDuctileThermalSheet::try_new(
                loaded.sheet,
                &later,
                sources(&fixture, &[PARTITION]),
                cx
            )
            .is_err()
        );
        let body = FiniteDuctileThermalSheet::try_new(
            native,
            &thermal,
            sources(&fixture, &[PARTITION]),
            cx,
        )
        .unwrap();
        let observer = Trace::default();
        assert!(matches!(
            body.advance(
                &thermal,
                &observer,
                1.0,
                &prescribed,
                &[],
                increments(),
                motion_controls(),
                thermal_config(),
                FiniteThermoplasticConfig {
                    maximum_iterations: 1,
                    ..config()
                },
                cx
            ),
            Err(FiniteThermoplasticError::NotConverged { .. })
        ));
        let refusing = Trace {
            refuse: true,
            ..Trace::default()
        };
        assert!(matches!(
            body.advance(
                &thermal,
                &refusing,
                1.0,
                &prescribed,
                &[],
                increments(),
                motion_controls(),
                thermal_config(),
                config(),
                cx
            ),
            Err(FiniteThermoplasticError::Observer(_))
        ));
        assert_eq!(observer, Trace::default());
        assert!(refusing.endpoints.is_empty());
        assert_eq!(thermal.cells(), &original_cells);
        assert_eq!(body.time_s(), 0.0);
        assert!(
            body.mechanics()
                .quadrature_states()
                .iter()
                .all(|s| s.native_state().alpha == 0.0 && s.damage() == 0.0)
        );
    });
}
