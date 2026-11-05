//! G1: sensible film energy stays separate from the host's exact latent heat.

use super::*;
use fs_blake3::ContentHash;
use fs_conduction::transient::enthalpy::heterogeneous::EnthalpyRadiation;
use fs_conduction::transient::enthalpy::sdirk::{EnthalpySdirk2Input, SDIRK2_GAMMA as G};
use fs_conduction::transient::enthalpy::{
    EnthalpyBackwardEuler, EnthalpyBudget, EnthalpyStepConfig,
};
use fs_exec::Cx;
use fs_material::phase::{EnthalpyPhaseKnot, EquilibriumEnthalpyPhaseCurve};
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
