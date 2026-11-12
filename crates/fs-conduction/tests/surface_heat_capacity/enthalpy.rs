//! G1: sensible film energy stays separate from the host's exact latent heat.

use super::*;
use fs_blake3::ContentHash;
use fs_conduction::transient::enthalpy::heterogeneous::EnthalpyRadiation;
use fs_conduction::transient::enthalpy::sdirk::{EnthalpySdirk2Input, SDIRK2_GAMMA as G};
use fs_conduction::transient::enthalpy::{
    EnthalpyBackwardEuler, EnthalpyBudget, EnthalpyError, EnthalpyStepConfig, SurfaceHeatStorage,
};
use fs_exec::Cx;
use fs_material::phase::{EnthalpyPhaseKnot, EquilibriumEnthalpyPhaseCurve};
use fs_material::thermal::transport::storage::{
    HeatStorageAmount, HeatStorageComponentQuery, PreparedHeatStorage, TabulatedHeatStorage,
};
use fs_solver::NewtonKrylovConfig;

const VOLUME: f64 = 8.0 / 3.0;
const DENSITY: f64 = 0.375;

fn chart() -> EquilibriumEnthalpyPhaseCurve {
    // One kilogram of host: Cp_s=2, latent heat=60, Cp_l=3 in SI units.
    EquilibriumEnthalpyPhaseCurve::try_new(
        ContentHash([0x67; 32]),
        [
            (0.0, 280.0, 0.0),
            (40.0, 300.0, 0.0),
            (100.0, 300.0, 1.0),
            (220.0, 340.0, 1.0),
        ]
        .map(|(h, t, liquid)| EnthalpyPhaseKnot {
            specific_enthalpy_j_kg: h,
            temperature_k: t,
            liquid_mass_fraction: liquid,
            bulk_density_kg_m3: DENSITY,
        })
        .to_vec(),
    )
    .unwrap()
}

fn engine<'m, 'c>(
    cx: &Cx<'_>,
    mesh: &'m ConductionMesh,
    chart: &'c EquilibriumEnthalpyPhaseCurve,
) -> EnthalpyBackwardEuler<'m, 'c> {
    // Every face has area 2 sqrt(3). Four films each store 1 J/K, giving
    // C_film=4 J/K and exactly 1 J/K at each of the four original nodes.
    let capacity = areal(1.0 / (2.0 * 3.0_f64.sqrt()));
    let films: Vec<_> = (0..mesh.boundary().len())
        .map(|face| SurfaceHeatCapacity::uniform(face, capacity))
        .collect();
    EnthalpyBackwardEuler::uniform(
        cx,
        mesh,
        chart,
        DENSITY,
        EnthalpyBudget {
            max_vertices: 4,
            max_elements: 1,
        },
    )
    .unwrap()
    .with_surface_heat_capacity(cx, &films)
    .unwrap()
}

fn controls() -> EnthalpyStepConfig {
    EnthalpyStepConfig {
        newton: NewtonKrylovConfig {
            absolute_tolerance: 1e-10,
            relative_tolerance: 1e-12,
            linear_restart: 4,
            max_linear_cycles: 8,
            forcing_minimum: 1e-12,
            forcing_maximum: 1e-4,
            ..NewtonKrylovConfig::default()
        },
        max_newton_iterations: 32,
        energy_tolerance_j: 1e-8,
    }
}

fn problem<'a>(
    mesh: &'a ConductionMesh,
    boundary: &'a fs_conduction::ThermalBoundary,
    material: &'a ConductivityModel,
    source: &'a ScalarField,
) -> ConductionProblem<'a> {
    ConductionProblem {
        mesh,
        boundary,
        material,
        element_materials: None,
        source,
    }
}

fn input(problem: ConductionProblem<'_>) -> EnthalpySdirk2Input<'_> {
    EnthalpySdirk2Input {
        problem,
        interfaces: None,
        radiation: EnthalpyRadiation::None,
        contacts: None,
    }
}

#[test]
fn uniform_heating_stores_film_sensible_heat_without_inflating_host_latent_heat() {
    with_cx(|cx| {
        let mesh = regular_tetrahedron();
        let boundary = coated_tetrahedron_boundary(cx, &mesh);
        let chart = chart();
        let engine = engine(cx, &mesh, &chart);
        let material = ConductivityModel::isotropic_declared(2.0).unwrap();
        let source = ScalarField::Uniform(3.75); // q V = 10 W.
        let p = problem(&mesh, &boundary, &material, &source);
        // Symmetric geometry, storage and forcing keep T uniform, so Fourier
        // transport vanishes. Integrate the physical inventory E=h+4*T:
        // dT/dt=10/(2+4), then dh/dt=10 on the exact latent plateau,
        // then dT/dt=10/(3+4). No film latent heat is supplied or inferred.
        for sdirk in [false, true] {
            let mut h = vec![20.0; 4];
            let mut temperature = vec![290.0; 4];
            let mut previous_time = 0.0;
            let mut accumulated = 0.0;
            for (time, expected_h, expected_t, expected_liquid) in [
                (3.0, 30.0, 295.0, 0.0),
                (6.0, 40.0, 300.0, 0.0),
                (9.0, 70.0, 300.0, 0.5),
                (12.0, 100.0, 300.0, 1.0),
                (19.0, 130.0, 310.0, 1.0),
            ] {
                let dt = time - previous_time;
                let (next_h, next_t, liquid, stored, defect) = if sdirk {
                    let step = engine
                        .advance_sdirk2(cx, [input(p); 2], &h, dt, controls())
                        .unwrap();
                    // Independently check the first physical stage inventory,
                    // rather than admitting a correct final sum alone.
                    for i in 0..4 {
                        let stage_storage = 0.25
                            * (step.first_stage.specific_enthalpy_j_kg[i] - h[i])
                            + (step.first_stage.temperature[i] - temperature[i]);
                        close(stage_storage, 2.5 * G * dt, 2e-8);
                    }
                    close(step.source_j, 10.0 * dt, 1e-10);
                    close(step.dirichlet_in_j, 0.0, 1e-12);
                    assert!(
                        step.nodal_residual_j
                            .iter()
                            .all(|r| r.abs() < controls().energy_tolerance_j)
                    );
                    (
                        step.second_stage.specific_enthalpy_j_kg,
                        step.second_stage.temperature,
                        step.second_stage.liquid_mass_fraction,
                        step.stored_energy_change_j,
                        step.energy_residual_j,
                    )
                } else {
                    let step = engine.advance(cx, p, None, &h, dt, controls()).unwrap();
                    close(step.source_w, 10.0, 1e-12);
                    close(step.dirichlet_in_w, 0.0, 1e-12);
                    (
                        step.specific_enthalpy_j_kg,
                        step.temperature,
                        step.liquid_mass_fraction,
                        step.stored_energy_change_j,
                        step.energy_residual_j,
                    )
                };
                for i in 0..4 {
                    close(next_h[i], expected_h, 3e-8);
                    close(next_t[i], expected_t, 2e-8);
                    close(liquid[i], expected_liquid, 1e-9);
                }
                close(stored, 10.0 * dt, 3e-8);
                close(defect, 0.0, controls().energy_tolerance_j);
                accumulated += stored;
                // At t=6 and t=12 the films have the same temperature: all
                // 60 J supplied between them belongs to the host's latent heat.
                close(
                    accumulated,
                    expected_h - 20.0 + 4.0 * (expected_t - 290.0),
                    8e-8,
                );
                h = next_h;
                temperature = next_t;
                previous_time = time;
            }
            close(accumulated, 190.0, 8e-8); // 110 J host + 80 J film.
        }
    });
}

#[test]
fn film_history_adjoints_use_the_old_and_new_caloric_slopes() {
    with_cx(|cx| {
        let mesh = regular_tetrahedron();
        let boundary = coated_tetrahedron_boundary(cx, &mesh);
        let chart = chart();
        let engine = engine(cx, &mesh, &chart);
        let material = ConductivityModel::isotropic_declared(2.0).unwrap();
        let source = ScalarField::Uniform(3.75);
        let p = problem(&mesh, &boundary, &material, &source);
        let linear = LinearConfig {
            tolerance: 1e-11,
            max_iterations: 64,
            restart: 4,
        };
        // Initial h=20 is solid (dT/dh=1/2), while h=70 is latent
        // (dT/dh=0). Both first SDIRK stages are latent; both final stages
        // are liquid (dT/dh=1/3), so all three storage slopes matter.
        // dT_new/dh_old = (M+Cfilm*dT_old/dh)/(M*Cp_liquid+Cfilm).
        for (old_h, dt, expected_t, expected_history) in [
            (20.0, 30.0, 2280.0 / 7.0, 3.0 / 7.0),
            (70.0, 10.0, 310.0, 1.0 / 7.0),
        ] {
            let old = [old_h; 4];
            let accepted = engine.advance(cx, p, None, &old, dt, controls()).unwrap();
            let response = engine
                .linearize_accepted(cx, p, None, &old, dt, controls(), accepted)
                .unwrap();
            let seed = response.temperature_pullback(cx, &[0.25; 4]).unwrap();
            let gradient = response.pullback(cx, &seed, linear).unwrap();
            for i in 0..4 {
                close(response.primal().temperature[i], expected_t, 2e-8);
                close(
                    gradient.previous_specific_enthalpy[i],
                    expected_history / 4.0,
                    2e-9,
                );
            }
            close(
                gradient.source_density.iter().sum(),
                dt * VOLUME / 7.0,
                2e-8,
            );

            let response = engine
                .linearize_sdirk2(cx, [input(p); 2], &old, dt, controls(), 0)
                .unwrap();
            let seeds = response
                .temperature_pullback(cx, [&[0.0; 4], &[0.25; 4]])
                .unwrap();
            let gradient = response
                .pullback_stages(cx, [&seeds[0], &seeds[1]], linear)
                .unwrap();
            for i in 0..4 {
                close(response.primal().first_stage.temperature[i], 300.0, 2e-8);
                close(
                    response.primal().second_stage.temperature[i],
                    expected_t,
                    2e-8,
                );
                close(
                    gradient.previous_specific_enthalpy[i],
                    expected_history / 4.0,
                    2e-9,
                );
            }
            for (stage, weight) in [1.0 - G, G].into_iter().enumerate() {
                close(
                    gradient.source_density[stage].iter().sum(),
                    weight * dt * VOLUME / 7.0,
                    2e-8,
                );
            }
        }
    });
}

#[test]
fn prescribed_baths_pay_for_the_same_film_storage_once_in_each_integrator() {
    with_cx(|cx| {
        let mesh = regular_tetrahedron();
        let chart = chart();
        let engine = engine(cx, &mesh, &chart);
        let material = ConductivityModel::isotropic_declared(2.0).unwrap();
        let baths = [295.0, 310.0].map(|t| {
            ThermalBoundaryBuilder::new(&mesh)
                .region("bath", |_| true, ThermalBc::dirichlet(t).unwrap())
                .unwrap()
                .finish()
                .unwrap()
        });
        let sources = [ScalarField::Uniform(3.75), ScalarField::Uniform(-1.875)];
        let old = [20.0; 4];
        let dt = 0.4;
        let expected_storage = 130.0 - 20.0 + 4.0 * (310.0 - 290.0);
        let step = engine
            .advance(
                cx,
                problem(&mesh, &baths[1], &material, &sources[0]),
                None,
                &old,
                dt,
                controls(),
            )
            .unwrap();
        close(step.stored_energy_change_j, expected_storage, 1e-9);
        close(dt * step.dirichlet_in_w, expected_storage - dt * 10.0, 1e-9);
        close(step.energy_residual_j, 0.0, 1e-9);

        let stages = [0, 1].map(|i| input(problem(&mesh, &baths[i], &material, &sources[i])));
        let step = engine
            .advance_sdirk2(cx, stages, &old, dt, controls())
            .unwrap();
        let source_heat = dt * ((1.0 - G) * 10.0 - G * 5.0);
        close(step.stored_energy_change_j, expected_storage, 1e-9);
        close(step.source_j, source_heat, 1e-10);
        close(step.dirichlet_in_j, expected_storage - source_heat, 1e-9);
        close(step.energy_residual_j, 0.0, 1e-9);
        assert_eq!(step.first_stage.temperature, [295.0; 4]);
        assert_eq!(step.second_stage.temperature, [310.0; 4]);
        assert_eq!(old, [20.0; 4]);
    });
}

fn bare_engine<'m, 'c>(
    cx: &Cx<'_>,
    mesh: &'m ConductionMesh,
    chart: &'c EquilibriumEnthalpyPhaseCurve,
) -> EnthalpyBackwardEuler<'m, 'c> {
    EnthalpyBackwardEuler::uniform(
        cx,
        mesh,
        chart,
        DENSITY,
        EnthalpyBudget {
            max_vertices: 4,
            max_elements: 1,
        },
    )
    .unwrap()
}

fn nonlinear_engine<'m, 'c>(
    cx: &Cx<'_>,
    mesh: &'m ConductionMesh,
    chart: &'c EquilibriumEnthalpyPhaseCurve,
    law: &'c TabulatedHeatStorage,
) -> EnthalpyBackwardEuler<'m, 'c> {
    let films: Vec<_> = (0..mesh.boundary().len())
        .map(|face| SurfaceHeatStorage::tabulated(face, law))
        .collect();
    bare_engine(cx, mesh, chart)
        .with_surface_heat_storage(cx, &films)
        .unwrap()
}

fn nonlinear_table() -> TabulatedHeatStorage {
    // Each face's complete inventory: C(T)=1+0.02*(T-300) J/K.
    TabulatedHeatStorage::try_new(&[[280.0, 0.6], [300.0, 1.0], [340.0, 1.8]]).unwrap()
}

// Independent closed-form integration of the four physical face inventories.
fn nonlinear_film_energy(old_t: f64, new_t: f64) -> f64 {
    (new_t - old_t) * (4.0 + 0.04 * (new_t + old_t - 600.0))
}

#[test]
fn nonlinear_films_integrate_calorics_through_host_melting_in_both_integrators() {
    with_cx(|cx| {
        let mesh = regular_tetrahedron();
        let boundary = coated_tetrahedron_boundary(cx, &mesh);
        let chart = chart();
        let law = nonlinear_table();
        let engine = nonlinear_engine(cx, &mesh, &chart, &law);
        let material = ConductivityModel::isotropic_declared(2.0).unwrap();
        let source = ScalarField::Uniform(3.75); // 10 W into one kilogram.
        let p = problem(&mesh, &boundary, &material, &source);
        assert!(std::ptr::eq(
            engine.surface_heat_storage()[0]
                .tabulated_storage()
                .unwrap(),
            &law
        ));
        for sdirk in [false, true] {
            let mut old = vec![20.0; 4];
            let mut previous_t = 290.0;
            let mut previous_time = 0.0;
            // Each time is the 10 W source's cumulative inventory from
            // (h=20, T=290): (h-20) + nonlinear_film_energy(290, T), over 10 W.
            // To (40, 300): 20 + 10*(4 - 0.4) = 56 J, t = 5.6 (not 5.8, which
            // asks the step to store both 29 J and 31 J below); the plateau
            // then takes 30 J per 3 s, and (130, 310) needs 110 + 80 = 190 J.
            for (time, expected_h, expected_t) in [
                (2.7, 30.0, 295.0),
                (5.6, 40.0, 300.0),
                (8.6, 70.0, 300.0),
                (11.6, 100.0, 300.0),
                (19.0, 130.0, 310.0),
                (31.6, 175.0, 325.0),
            ] {
                let dt = time - previous_time;
                let (next_h, next_t, stored, defect) = if sdirk {
                    let step = engine
                        .advance_sdirk2(cx, [input(p); 2], &old, dt, controls())
                        .unwrap();
                    (
                        step.second_stage.specific_enthalpy_j_kg,
                        step.second_stage.temperature,
                        step.stored_energy_change_j,
                        step.energy_residual_j,
                    )
                } else {
                    let step = engine.advance(cx, p, None, &old, dt, controls()).unwrap();
                    (
                        step.specific_enthalpy_j_kg,
                        step.temperature,
                        step.stored_energy_change_j,
                        step.energy_residual_j,
                    )
                };
                for (&h, &t) in next_h.iter().zip(&next_t) {
                    close(h, expected_h, 3e-8);
                    close(t, expected_t, 2e-8);
                }
                close(
                    stored,
                    expected_h - old[0] + nonlinear_film_energy(previous_t, expected_t),
                    5e-8,
                );
                close(stored, dt * 10.0, 5e-8);
                close(defect, 0.0, controls().energy_tolerance_j);
                old = next_h;
                previous_t = expected_t;
                previous_time = time;
            }
        }
    });
}

#[test]
fn nonlinear_film_adjoints_use_separate_old_stage_and_endpoint_capacities() {
    with_cx(|cx| {
        let mesh = regular_tetrahedron();
        let boundary = coated_tetrahedron_boundary(cx, &mesh);
        let chart = chart();
        let law = nonlinear_table();
        let engine = nonlinear_engine(cx, &mesh, &chart, &law);
        let material = ConductivityModel::isotropic_declared(2.0).unwrap();
        let source = ScalarField::Uniform(3.75);
        let p = problem(&mesh, &boundary, &material, &source);
        let linear = LinearConfig {
            tolerance: 1e-11,
            max_iterations: 64,
            restart: 4,
        };
        // Final T=310: total dE/dT = Cp_liquid + 4*C_face(310) = 7.8.
        // Old dE/dh is 2.6 on the solid branch, or 1 on the latent plateau.
        for (old_h, dt, old_slope) in [(20.0, 19.0, 2.6), (70.0, 10.4, 1.0)] {
            let old = [old_h; 4];
            let accepted = engine.advance(cx, p, None, &old, dt, controls()).unwrap();
            let response = engine
                .linearize_accepted(cx, p, None, &old, dt, controls(), accepted)
                .unwrap();
            let seed = response.temperature_pullback(cx, &[0.25; 4]).unwrap();
            let gradient = response.pullback(cx, &seed, linear).unwrap();
            for i in 0..4 {
                close(response.primal().temperature[i], 310.0, 2e-8);
                close(
                    gradient.previous_specific_enthalpy[i],
                    old_slope / (4.0 * 7.8),
                    2e-9,
                );
            }
            close(
                gradient.source_density.iter().sum(),
                dt * VOLUME / 7.8,
                2e-8,
            );
            let epsilon = 1e-3;
            let lower = engine
                .advance(cx, p, None, &[old_h - epsilon; 4], dt, controls())
                .unwrap();
            let upper = engine
                .advance(cx, p, None, &[old_h + epsilon; 4], dt, controls())
                .unwrap();
            close(
                (upper.temperature[0] - lower.temperature[0]) / (2.0 * epsilon),
                old_slope / 7.8,
                2e-7,
            );

            let response = engine
                .linearize_sdirk2(cx, [input(p); 2], &old, dt, controls(), 0)
                .unwrap();
            let seeds = response
                .temperature_pullback(cx, [&[0.0; 4], &[0.25; 4]])
                .unwrap();
            let gradient = response
                .pullback_stages(cx, [&seeds[0], &seeds[1]], linear)
                .unwrap();
            for i in 0..4 {
                close(response.primal().second_stage.temperature[i], 310.0, 2e-8);
                close(
                    gradient.previous_specific_enthalpy[i],
                    old_slope / (4.0 * 7.8),
                    2e-9,
                );
            }
            for (stage, weight) in [1.0 - G, G].into_iter().enumerate() {
                close(
                    gradient.source_density[stage].iter().sum(),
                    weight * dt * VOLUME / 7.8,
                    2e-8,
                );
            }
        }
    });
}

fn prepared_film() -> PreparedHeatStorage {
    let capacity = PropertyKey::new("film-specific-capacity", SPECIFIC_HEAT_CAPACITY_DIMS);
    let temperature = fs_qty::QuantitySpec::dimensional(fs_qty::Temperature::DIMS);
    let mut claims = ClaimSet::new();
    claims
        .insert_claim(PropertyClaim {
            key: capacity.clone(),
            value: PropertyValue::Curve {
                abscissa: "T".into(),
                abscissa_dims: fs_qty::Temperature::DIMS,
                knots: vec![(285.0, 1.4), (315.0, 2.6)],
                dims: SPECIFIC_HEAT_CAPACITY_DIMS,
            },
            validity: ValidityDomain::unconstrained().with_quantity("T", temperature, 285.0, 315.0),
            uncertainty: UncertaintyModel::Unstated,
            interpolation: InterpolationPolicy::LinearInside,
            observations: Vec::new(),
            provenance: Provenance {
                source: "manufactured film; not measured material data".into(),
                license: "CC0-1.0".into(),
                artifact: None,
            },
        })
        .unwrap();
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic film".into(),
            phase: "solid".into(),
            process: "fixed mass".into(),
            revision: 0,
        },
        claims,
        Vec::new(),
    )
    .unwrap();
    let point = |t| {
        QueryPoint::new()
            .with_quantity("T", temperature, t)
            .unwrap()
    };
    PreparedHeatStorage::try_new(&[HeatStorageComponentQuery {
        id: 7,
        card: &card,
        capacity: &capacity,
        amount: HeatStorageAmount::MassKg(0.5),
        temperature_axis: "T",
        lower: &point(285.0),
        upper: &point(315.0),
        selection: fs_matdb::ClaimSelection::Policy(fs_matdb::SelectionPolicy::SingleClaimOnly),
    }])
    .unwrap()
}

#[test]
fn original_film_sources_bound_histories_candidates_and_prescribed_bath_heat() {
    with_cx(|cx| {
        let mesh = regular_tetrahedron();
        let chart = chart();
        let law = prepared_film();
        let faces: Vec<_> = (0..4)
            .map(|face| SurfaceHeatStorage::prepared(face, &law))
            .collect();
        let engine = bare_engine(cx, &mesh, &chart)
            .with_surface_heat_storage(cx, &faces)
            .unwrap();
        assert!(std::ptr::eq(
            engine.surface_heat_storage()[0].prepared_storage().unwrap(),
            &law
        ));
        assert_eq!(law.components()[0].amount(), HeatStorageAmount::MassKg(0.5));
        let material = ConductivityModel::isotropic_declared(2.0).unwrap();
        let source = ScalarField::Uniform(3.75);
        let baths = [295.0, 310.0, 320.0].map(|t| {
            ThermalBoundaryBuilder::new(&mesh)
                .region("bath", |_| true, ThermalBc::dirichlet(t).unwrap())
                .unwrap()
                .finish()
                .unwrap()
        });
        let old = [20.0; 4];
        let dt = 0.4;
        let step = engine
            .advance(
                cx,
                problem(&mesh, &baths[1], &material, &source),
                None,
                &old,
                dt,
                controls(),
            )
            .unwrap();
        close(step.stored_energy_change_j, 190.0, 1e-9);
        close(dt * step.dirichlet_in_w, 190.0 - dt * 10.0, 1e-9);
        let stages = [0, 1].map(|i| input(problem(&mesh, &baths[i], &material, &source)));
        let step = engine
            .advance_sdirk2(cx, stages, &old, dt, controls())
            .unwrap();
        close(step.stored_energy_change_j, 190.0, 1e-9);
        close(step.dirichlet_in_j, 190.0 - dt * 10.0, 1e-9);
        close(step.energy_residual_j, 0.0, 1e-9);
        for (history, bath) in [([0.0; 4], 1), (old, 2)] {
            assert!(matches!(
                engine.advance(
                    cx,
                    problem(&mesh, &baths[bath], &material, &source),
                    None,
                    &history,
                    dt,
                    controls()
                ),
                Err(EnthalpyError::Conduction(
                    ConductionError::OutsideTemperatureSpan { .. }
                ))
            ));
        }
        assert_eq!(old, [20.0; 4]);
        assert_eq!(law.components()[0].amount(), HeatStorageAmount::MassKg(0.5));
    });
}

#[test]
fn constant_nonlinear_films_match_the_existing_path_and_heterogeneous_owner() {
    use fs_conduction::transient::enthalpy::heterogeneous::{
        HeterogeneousEnthalpyBackwardEuler, ReferenceEnthalpyMaterial,
    };
    with_cx(|cx| {
        let mesh = regular_tetrahedron();
        let chart = chart();
        let boundary = coated_tetrahedron_boundary(cx, &mesh);
        let law = TabulatedHeatStorage::try_new(&[[280.0, 1.0], [340.0, 1.0]]).unwrap();
        let original = engine(cx, &mesh, &chart);
        let nonlinear = nonlinear_engine(cx, &mesh, &chart, &law);
        let faces: Vec<_> = (0..4)
            .map(|face| SurfaceHeatStorage::tabulated(face, &law))
            .collect();
        let heterogeneous = HeterogeneousEnthalpyBackwardEuler::new(
            cx,
            &mesh,
            &[ReferenceEnthalpyMaterial {
                curve: &chart,
                reference_density_kg_m3: DENSITY,
            }],
            &[0],
            EnthalpyBudget {
                max_vertices: 4,
                max_elements: 1,
            },
        )
        .unwrap()
        .with_surface_heat_storage(cx, &faces)
        .unwrap();
        let material = ConductivityModel::isotropic_declared(2.0).unwrap();
        let source = ScalarField::Uniform(0.0);
        let p = problem(&mesh, &boundary, &material, &source);
        let old = [21.0, 25.0, 29.0, 33.0];
        let reference = original
            .advance(cx, p, None, &old, 0.2, controls())
            .unwrap();
        for result in [
            nonlinear
                .advance(cx, p, None, &old, 0.2, controls())
                .unwrap(),
            heterogeneous
                .advance(cx, p, None, &old, 0.2, controls())
                .unwrap(),
        ] {
            for (actual, expected) in result
                .specific_enthalpy_j_kg
                .iter()
                .zip(&reference.specific_enthalpy_j_kg)
            {
                close(*actual, *expected, 1e-8);
            }
            close(result.stored_energy_change_j, 0.0, 1e-8);
        }
        assert!(std::ptr::eq(
            heterogeneous.surface_heat_storage()[0]
                .tabulated_storage()
                .unwrap(),
            &law
        ));
    });
}

#[test]
fn nonlinear_film_admission_and_two_sided_storage_derivatives_refuse_at_support_edges() {
    use fs_conduction::transient::enthalpy::adjoint::EnthalpyAdjointError;
    with_cx(|cx| {
        let mesh = regular_tetrahedron();
        let chart = chart();
        let law = TabulatedHeatStorage::try_new(&[[290.0, 0.8], [310.0, 1.2]]).unwrap();
        let face = SurfaceHeatStorage::tabulated(0, &law);
        assert!(
            bare_engine(cx, &mesh, &chart)
                .with_surface_heat_storage(cx, &[face, face])
                .is_err()
        );
        assert!(
            bare_engine(cx, &mesh, &chart)
                .with_surface_heat_storage(cx, &[SurfaceHeatStorage::tabulated(4, &law)])
                .is_err()
        );
        let outside = TabulatedHeatStorage::try_new(&[[350.0, 1.0], [360.0, 1.0]]).unwrap();
        assert!(
            bare_engine(cx, &mesh, &chart)
                .with_surface_heat_storage(cx, &[SurfaceHeatStorage::tabulated(0, &outside)])
                .is_err()
        );
        with_cancelled_cx(|cancelled| {
            assert!(
                bare_engine(cx, &mesh, &chart)
                    .with_surface_heat_storage(cancelled, &[face])
                    .is_err()
            );
        });
        let engine = nonlinear_engine(cx, &mesh, &chart, &law);
        let boundary = coated_tetrahedron_boundary(cx, &mesh);
        let material = ConductivityModel::isotropic_declared(2.0).unwrap();
        let source = ScalarField::Uniform(3.75);
        let p = problem(&mesh, &boundary, &material, &source);
        let old = [20.0; 4]; // Valid primal history, but T=290 is the film support edge.
        let step = engine.advance(cx, p, None, &old, 1.0, controls()).unwrap();
        assert!(matches!(
            engine.linearize_accepted(cx, p, None, &old, 1.0, controls(), step),
            Err(EnthalpyAdjointError::FilmStorageEndpoint { .. })
        ));
        let plateau = TabulatedHeatStorage::try_new(&[[300.0, 1.0], [310.0, 1.2]]).unwrap();
        let engine = nonlinear_engine(cx, &mesh, &chart, &plateau);
        let old = [70.0; 4];
        let step = engine.advance(cx, p, None, &old, 1.0, controls()).unwrap();
        let response = engine
            .linearize_accepted(cx, p, None, &old, 1.0, controls(), step)
            .unwrap();
        assert_eq!(
            response.temperature_pullback(cx, &[0.25; 4]).unwrap(),
            vec![0.0; 4]
        );
        assert_eq!(old, [70.0; 4]);
    });
}

#[test]
fn nonlinear_inventory_profiles_retain_exact_p1_moments_and_original_amounts() {
    let law = prepared_film();
    let original = SurfaceHeatStorage::prepared(2, &law);
    // Integral(N_i*N_j) on a triangle is A/6 on the diagonal and A/12
    // elsewhere; total P1 inventory is A*sum(a)/3. These are the resulting
    // row fractions for a=[0,1,3], independently evaluated as rational values.
    let expected = [1.0 / 4.0, 5.0 / 16.0, 7.0 / 16.0];
    let tiny = f64::from_bits(1);
    for profile in [
        [0.0, 1.0, 3.0],
        [0.0, tiny, 3.0 * tiny],
        [0.0, 2.0_f64.powi(1020), 3.0 * 2.0_f64.powi(1020)],
    ] {
        let face = original.with_inventory_profile(profile).unwrap();
        assert_eq!(face.boundary_face(), 2);
        assert!(std::ptr::eq(face.prepared_storage().unwrap(), &law));
        for (actual, expected) in face.nodal_inventory_fractions().into_iter().zip(expected) {
            close(actual, expected, 2e-16);
        }
        close(face.nodal_inventory_fractions().iter().sum(), 1.0, 2e-16);
    }
    let reversed = original.with_inventory_profile([3.0, 0.0, 1.0]).unwrap();
    for (actual, expected) in reversed.nodal_inventory_fractions().into_iter().zip([
        expected[2],
        expected[0],
        expected[1],
    ]) {
        close(actual, expected, 2e-16);
    }
    for value in [tiny, 1.0, f64::MAX] {
        assert_eq!(
            original
                .with_inventory_profile([value; 3])
                .unwrap()
                .nodal_inventory_fractions(),
            original.nodal_inventory_fractions()
        );
    }
    for profile in [
        [0.0; 3],
        [-1.0, 1.0, 1.0],
        [f64::INFINITY, 1.0, 1.0],
        [f64::NAN, 1.0, 1.0],
    ] {
        assert!(original.with_inventory_profile(profile).is_err());
    }
    assert_eq!(law.components()[0].amount(), HeatStorageAmount::MassKg(0.5));
}

#[test]
fn profiled_original_films_charge_nonuniform_baths_and_entropy_in_both_integrators() {
    with_cx(|cx| {
        let mesh = regular_tetrahedron();
        let chart = chart();
        let law = prepared_film(); // Total C(T)=0.02*T-5 J/K, original mass=0.5 kg.
        let face = SurfaceHeatStorage::prepared(0, &law)
            .with_inventory_profile([0.0, 1.0, 3.0])
            .unwrap();
        let engine = bare_engine(cx, &mesh, &chart)
            .with_surface_heat_storage(cx, &[face])
            .unwrap();
        let mut shares = [0.0; 4];
        for (vertex, fraction) in
            mesh.boundary()[0]
                .vertices
                .into_iter()
                .zip([1.0 / 4.0, 5.0 / 16.0, 7.0 / 16.0])
        {
            shares[vertex as usize] = fraction;
        }
        let old_t = [289.0, 290.0, 291.0, 292.0];
        let new_t = [294.0, 296.0, 295.0, 297.0];
        let old_h = old_t.map(|t| 2.0 * (t - 280.0));
        let new_h = new_t.map(|t| 2.0 * (t - 280.0));
        let middle = std::array::from_fn(|i| 0.5 * (old_t[i] + new_t[i]));
        let baths = [middle, new_t].map(|temperatures: [f64; 4]| {
            ThermalBoundaryBuilder::new(&mesh)
                .region(
                    "bath",
                    |_| true,
                    ThermalBc::Dirichlet {
                        temperature: ScalarField::Nodal(temperatures.to_vec()),
                    },
                )
                .unwrap()
                .finish()
                .unwrap()
        });
        let expected_energy: f64 = (0..4)
            .map(|i| {
                let delta = new_t[i] - old_t[i];
                0.5 * delta + shares[i] * delta * (0.01 * (new_t[i] + old_t[i]) - 5.0)
            })
            .sum();
        let expected_entropy: f64 = (0..4)
            .map(|i| {
                let logarithm = (new_t[i] / old_t[i]).ln();
                0.5 * logarithm + shares[i] * (0.02 * (new_t[i] - old_t[i]) - 5.0 * logarithm)
            })
            .sum();
        close(
            engine
                .thermal_entropy_change_j_k(cx, &old_h, &new_h)
                .unwrap(),
            expected_entropy,
            2e-14,
        );
        let material = ConductivityModel::isotropic_declared(2.0).unwrap();
        let source = ScalarField::Uniform(0.0);
        let dt = 0.4;
        let be = engine
            .advance(
                cx,
                problem(&mesh, &baths[1], &material, &source),
                None,
                &old_h,
                dt,
                controls(),
            )
            .unwrap();
        close(be.stored_energy_change_j, expected_energy, 2e-10);
        close(dt * be.dirichlet_in_w, expected_energy, 2e-10);
        let stages = [0, 1].map(|i| input(problem(&mesh, &baths[i], &material, &source)));
        let sdirk = engine
            .advance_sdirk2(cx, stages, &old_h, dt, controls())
            .unwrap();
        close(sdirk.stored_energy_change_j, expected_energy, 2e-10);
        close(sdirk.dirichlet_in_j, expected_energy, 2e-10);
        close(sdirk.energy_residual_j, 0.0, 2e-10);
        for i in 0..4 {
            close(be.temperature[i], new_t[i], 2e-12);
            close(sdirk.second_stage.temperature[i], new_t[i], 2e-12);
        }
        assert!(std::ptr::eq(
            engine.surface_heat_storage()[0].prepared_storage().unwrap(),
            &law
        ));
    });
}

#[test]
fn profiled_nonlinear_storage_drives_the_original_transient_and_history_adjoint() {
    with_cx(|cx| {
        let mesh = regular_tetrahedron();
        let chart = chart();
        let law = nonlinear_table();
        let uniform = SurfaceHeatStorage::tabulated(0, &law);
        let film = uniform.with_inventory_profile([0.0, 1.0, 3.0]).unwrap();
        let engine = bare_engine(cx, &mesh, &chart)
            .with_surface_heat_storage(cx, &[film])
            .unwrap();
        let mut shares = [0.0; 4];
        for (vertex, fraction) in
            mesh.boundary()[0]
                .vertices
                .into_iter()
                .zip([1.0 / 4.0, 5.0 / 16.0, 7.0 / 16.0])
        {
            shares[vertex as usize] = fraction;
        }
        let boundary = ThermalBoundaryBuilder::new(&mesh)
            .adiabatic_remainder()
            .finish()
            .unwrap();
        let material = ConductivityModel::isotropic_declared(2.0).unwrap();
        let source = ScalarField::Uniform(0.0);
        let p = problem(&mesh, &boundary, &material, &source);
        let (target, dt): ([f64; 4], f64) = ([292.0, 293.0, 294.0, 295.0], 0.1);
        // Independent regular-tet K*T = (4/3)*(T-mean(T)). Each host row has
        // mass=1/4 and Cp=2. Invert its quadratic film-plus-host energy to
        // manufacture the OLD state for this explicitly specified endpoint.
        let old = std::array::from_fn::<_, 4, _>(|i| {
            let energy = dt * (4.0 / 3.0) * (target[i] - 293.5);
            let linear = 0.5 + shares[i] * (0.02 * target[i] - 5.0);
            let quadratic = 0.01 * shares[i];
            let delta =
                2.0 * energy / (linear + (linear * linear + 4.0 * quadratic * energy).sqrt());
            2.0 * (target[i] + delta - 280.0)
        });
        let step = engine.advance(cx, p, None, &old, dt, controls()).unwrap();
        for i in 0..4 {
            close(step.temperature[i], target[i], 2e-8);
        }
        close(step.stored_energy_change_j, 0.0, 2e-9);
        close(step.energy_residual_j, 0.0, 2e-9);
        let response = engine
            .linearize_accepted(cx, p, None, &old, dt, controls(), step)
            .unwrap();
        let seed = response
            .temperature_pullback(cx, &[1.0, 0.0, 0.0, 0.0])
            .unwrap();
        let gradient = response
            .pullback(
                cx,
                &seed,
                LinearConfig {
                    tolerance: 1e-11,
                    max_iterations: 64,
                    restart: 4,
                },
            )
            .unwrap();
        for i in 0..4 {
            let mut lower = old;
            let mut upper = old;
            lower[i] -= 1e-3;
            upper[i] += 1e-3;
            let a = engine.advance(cx, p, None, &lower, dt, controls()).unwrap();
            let b = engine.advance(cx, p, None, &upper, dt, controls()).unwrap();
            close(
                gradient.previous_specific_enthalpy[i],
                (b.temperature[0] - a.temperature[0]) / 2e-3,
                3e-7,
            );
        }
        // Explicit uniform profiles must keep the old native arithmetic too.
        let results = [
            uniform,
            uniform.with_inventory_profile([f64::MAX; 3]).unwrap(),
        ]
        .map(|face| {
            bare_engine(cx, &mesh, &chart)
                .with_surface_heat_storage(cx, &[face])
                .unwrap()
                .advance(cx, p, None, &old, dt, controls())
                .unwrap()
        });
        assert_eq!(
            results[0].specific_enthalpy_j_kg,
            results[1].specific_enthalpy_j_kg
        );
        assert_eq!(
            results[0].stored_energy_change_j,
            results[1].stored_energy_change_j
        );
    });
}
