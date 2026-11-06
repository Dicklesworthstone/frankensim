use super::*;
use fs_material::creep::plastic::{NortonJ2, NortonJ2State};
use fs_material::creep::{CreepState, NortonCreep, NortonCreepParameters, StrainHardeningCreep};
use fs_material::tensor::{contract, von_mises};
use fs_solid::plastic_shell::dynamics::PlasticShellDensity;

const PRIMARY_STRAIN: f64 = 8e-4;
const REFERENCE_RATE: f64 = 2e-3;

fn creep(duration_s: f64, reference_rate_s: f64) -> NortonCreep {
    NortonCreep::new(
        IsotropicElastic::new(E, NU, 0.05).unwrap(),
        NortonCreepParameters {
            reference_stress_pa: 1.0,
            reference_rate_s,
            stress_exponent: 1.0,
            activation_energy_j_mol: 0.0,
            reference_temperature_k: 300.0,
            temperature_range_k: [300.0, 600.0],
        },
        300.0,
        duration_s,
    )
    .unwrap()
    .with_strain_hardening(StrainHardeningCreep {
        reference_strain: PRIMARY_STRAIN,
        exponent: 1.0,
        residual_rate_fraction: 0.0,
    })
    .unwrap()
}

fn triangle() -> ShellMesh {
    ShellMesh::new(
        vec![[0.0; 3], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]],
        vec![[0, 1, 2]],
    )
    .unwrap()
}

// Positive root of d² + b d - c = 0, independent of the native return map.
fn quadratic(b: f64, c: f64) -> f64 {
    2.0 * c / (b + (b * b + 4.0 * c).sqrt())
}

#[test]
fn primary_shell_creep_matches_constant_stress_relaxation_and_permanent_extension() {
    with_cx(|_, cx| {
        let mesh = triangle();
        let dt = 0.25;
        let stress = 0.8;
        let volume = 0.5 * H;
        let material = section(creep(dt, REFERENCE_RATE));
        // Only axial extension and Poisson contraction are free. The original
        // CST membrane is exactly uniaxial; no section strain is imposed.
        let fixed: Vec<_> = (0..18)
            .filter(|d| *d != 6 && *d != 13)
            .map(|d| (d, 0.0))
            .collect();
        let loads = [(6, volume * stress)];
        let owner = problem(&mesh, &material, &fixed, &loads, &[1.0], 0.37);
        let mut previous: Option<SmallStrainShellSolution<CreepState>> = None;
        let mut alpha = 0.0;
        let mut previous_increment = f64::INFINITY;
        for interval in 0..4 {
            // dc = dt A sigma / (1 + (alpha_old + dc)/a).
            let increment = quadratic(
                PRIMARY_STRAIN + alpha,
                dt * REFERENCE_RATE * stress * PRIMARY_STRAIN,
            );
            alpha += increment;
            let step = owner
                .solve_with_energy(previous.as_ref(), gates(), cx)
                .unwrap();
            near(step.solution.coordinates[1][0], stress / E + alpha, 2e-12);
            near(
                step.solution.coordinates[2][1],
                -NU * stress / E - 0.5 * alpha,
                2e-12,
            );
            for state in &step.solution.section_states {
                for point in &state.points {
                    near(point.history.equivalent_creep_strain, alpha, 2e-12);
                    near(point.history.creep_strain[0], alpha, 2e-12);
                    near(point.history.creep_strain[1], -0.5 * alpha, 2e-12);
                    near(point.history.creep_strain[2], -0.5 * alpha, 2e-12);
                }
            }
            let stored = volume * stress * stress / (2.0 * E);
            near(step.energy.final_state.elastic_j, stored, 2e-12);
            near(step.energy.final_state.hardening_j, 0.0, 0.0);
            near(
                step.energy.physical_dissipation_j,
                volume * stress * increment,
                2e-12,
            );
            near(
                step.energy.material_numerical_dissipation_j,
                if interval == 0 { stored } else { 0.0 },
                2e-12,
            );
            near(step.energy.balance_residual_j, 0.0, 1e-9);
            assert!(increment < previous_increment);
            previous_increment = increment;
            if interval == 1 {
                let saved = previous.clone();
                let repeated = owner
                    .solve_with_energy(previous.as_ref(), gates(), cx)
                    .unwrap();
                assert_eq!(step.solution, repeated.solution);
                assert_eq!(step.energy, repeated.energy);
                assert_eq!(previous, saved);
            }
            previous = Some(step.solution);
        }

        // Hold the actual loaded length while the transverse contraction stays
        // free. The stress now relaxes through the same accumulated history.
        let loaded = previous.unwrap();
        let mut clamped = fixed.clone();
        clamped.push((6, loaded.coordinates[1][0]));
        let k = dt * REFERENCE_RATE * PRIMARY_STRAIN;
        let increment = quadratic(PRIMARY_STRAIN + alpha + k * E, k * stress);
        let relaxed_stress = stress - E * increment;
        alpha += increment;
        let relaxed = problem(&mesh, &material, &clamped, &[], &[1.0], 0.37)
            .solve_with_energy(Some(&loaded), gates(), cx)
            .unwrap();
        for resultants in &relaxed.solution.section_resultants {
            near(resultants[0], H * relaxed_stress, 2e-11);
            near(resultants[1], 0.0, 2e-11);
        }
        near(
            relaxed.energy.physical_dissipation_j,
            volume * relaxed_stress * increment,
            2e-12,
        );
        near(relaxed.energy.support_endpoint_work_j, 0.0, 2e-12);
        assert!(relaxed_stress > 0.0 && relaxed_stress < stress);
        let unloaded = problem(&mesh, &material, &fixed, &[], &[1.0], 0.37)
            .solve_with_energy(Some(&relaxed.solution), gates(), cx)
            .unwrap();
        near(unloaded.solution.coordinates[1][0], alpha, 2e-12);
        near(unloaded.solution.coordinates[2][1], -0.5 * alpha, 2e-12);
        near(unloaded.energy.physical_dissipation_j, 0.0, 2e-12);
        near(unloaded.energy.final_state.elastic_j, 0.0, 2e-12);
    });
}

#[test]
fn creeping_plastic_cantilever_keeps_distinct_histories_and_native_energy_through_unloading() {
    with_cx(|_, cx| {
        let (mesh, fixed, loads) = cantilever(2, 1);
        let law = NortonJ2::new(
            creep(0.25, 1e-4),
            J2Plasticity::new(IsotropicElastic::new(E, NU, 0.05).unwrap(), 1.0, 100.0).unwrap(),
        )
        .unwrap();
        let material = section(&law);
        let mut previous: Option<SmallStrainShellSolution<NortonJ2State>> = None;
        let mut cumulative_creep_heat = 0.0;
        let mut cumulative_plastic_heat = 0.0;
        let mut peak_rotation = 0.0;
        for (interval, factor) in [0.5, 1.0, 1.0, 0.0].into_iter().enumerate() {
            let step = problem(&mesh, &material, &fixed, &loads, &[factor], 0.37)
                .solve_with_energy(previous.as_ref(), gates(), cx)
                .unwrap();
            let mut creep_heat = 0.0;
            let mut plastic_heat = 0.0;
            let mut elastic_storage = 0.0;
            let mut hardening_storage = 0.0;
            for e in 0..mesh.tris.len() {
                for q in 0..3 {
                    for (i, point) in material.points().iter().enumerate() {
                        let after = &step.solution.section_states[3 * e + q].points[i];
                        let before = previous
                            .as_ref()
                            .map(|old| &old.section_states[3 * e + q].points[i].history);
                        let dc = after.history.creep.equivalent_creep_strain
                            - before.map_or(0.0, |s| s.creep.equivalent_creep_strain);
                        let dp =
                            after.history.plastic.alpha - before.map_or(0.0, |s| s.plastic.alpha);
                        assert!(dc >= 0.0 && dp >= 0.0);
                        let elastic_strain = core::array::from_fn(|j| {
                            after.strain[j]
                                - after.history.creep.creep_strain[j]
                                - after.history.plastic.plastic_strain[j]
                        });
                        let stress = law.creep().elastic().stress(&elastic_strain, &());
                        let weight = mesh.facet(e).unwrap().area_m2 / 3.0 * point.weight_m;
                        creep_heat += weight * von_mises(&stress) * dc;
                        plastic_heat += weight * law.plastic().yield_stress * dp;
                        elastic_storage += weight * 0.5 * contract(&elastic_strain, &stress);
                        hardening_storage += weight
                            * 0.5
                            * law.plastic().hardening
                            * after.history.plastic.alpha.powi(2);
                    }
                }
            }
            near(
                step.energy.physical_dissipation_j,
                creep_heat + plastic_heat,
                2e-12,
            );
            near(step.energy.final_state.elastic_j, elastic_storage, 2e-12);
            near(
                step.energy.final_state.hardening_j,
                hardening_storage,
                2e-12,
            );
            near(step.energy.final_state.kinematic_hardening_j, 0.0, 0.0);
            near(step.energy.balance_residual_j, 0.0, 1e-9);
            let tip_rotation: f64 = mesh
                .nodes
                .iter()
                .zip(&step.solution.coordinates)
                .filter(|(p, _)| p[0] == 1.0)
                .map(|(_, q)| q[4])
                .sum::<f64>()
                / 2.0;
            if interval == 1 {
                peak_rotation = tip_rotation;
                assert!(hardening_storage > 0.0);
            } else if interval == 2 {
                assert!(tip_rotation > peak_rotation);
                assert!(creep_heat > 0.0);
            } else if interval == 3 {
                assert!(tip_rotation > 1e-5 && tip_rotation < peak_rotation);
            }
            cumulative_creep_heat += creep_heat;
            cumulative_plastic_heat += plastic_heat;
            previous = Some(step.solution);
        }
        assert!(cumulative_creep_heat > 0.0 && cumulative_plastic_heat > 0.0);
    });
}

#[test]
fn creep_bending_motion_keeps_its_physical_clock_and_separate_inertial_loss() {
    with_cx(|_, cx| {
        let mesh = triangle();
        let dt = 0.25;
        let density = 2.0;
        let curvature = 0.02;
        let law = creep(dt, REFERENCE_RATE);
        let material = section(&law);
        let prescribed: Vec<_> = mesh
            .nodes
            .iter()
            .flat_map(|&[x, _, _]| [0.0, 0.0, 0.5 * curvature * x * x, 0.0, -curvature * x, 0.0])
            .enumerate()
            .collect();
        let owner = problem(&mesh, &material, &prescribed, &[], &[1.0], 0.37);
        let reference = owner
            .prepare_reference_inertia(PlasticShellDensity::Uniform(density), cx)
            .unwrap();
        let velocity = [[0.0; 6]; 3];
        let error = owner
            .advance_backward_euler_with_energy(None, &velocity, &reference, 2.0 * dt, gates(), cx)
            .unwrap_err();
        assert!(error.to_string().contains("durations must agree"));
        let step = owner
            .advance_backward_euler_with_energy(None, &velocity, &reference, dt, gates(), cx)
            .unwrap();
        // Only corner (1,0) moves; integral N1² dA = A/6. The physical
        // director inertia is rho*h³/12, independent of constitutive creep.
        let kinetic = 0.5
            * (density * H * 0.5 / 6.0)
            * ((0.5 * curvature / dt).powi(2) + H * H / 12.0 * (curvature / dt).powi(2));
        near(step.inertia.physical.kinetic_after_j, kinetic, 2e-15);
        near(step.inertia.physical.numerical_loss_j, kinetic, 2e-15);
        near(step.inertia.drilling_gauge.kinetic_after_j, 0.0, 1e-25);
        let mut physical_heat = 0.0;
        for state in &step.solution.section_states {
            for (point, after) in material.points().iter().zip(&state.points) {
                let elastic =
                    core::array::from_fn(|i| after.strain[i] - after.history.creep_strain[i]);
                let stress = law.elastic().stress(&elastic, &());
                physical_heat += (0.5 / 3.0)
                    * point.weight_m
                    * von_mises(&stress)
                    * after.history.equivalent_creep_strain;
            }
        }
        near(step.energy.physical_dissipation_j, physical_heat, 2e-12);
        assert!(physical_heat > 0.0 && step.energy.material_numerical_dissipation_j > 0.0);
        near(step.balance_residual_j(), 0.0, 1e-9);
    });
}
