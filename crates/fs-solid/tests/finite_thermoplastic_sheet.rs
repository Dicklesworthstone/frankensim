//! G0/G1/G3/G4 checks of finite thermal sheet deformation and native heat accounts.
//! All coefficients are synthetic; these are not experimental formability claims.

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_material::hyper::plane_stress::PlaneStressSettings;
use fs_material::plastic::finite::thermal::{
    AffineThermalFiniteJ2, FiniteJ2TemperatureSlopes, FiniteJ2ThermalOverstress,
    FiniteJ2ThermalParameters,
};
use fs_material::plastic::finite::{FiniteJ2Limits, FiniteJ2Parameters};
use fs_solid::hyper2d::plastic::SheetAssemblyBudget;
use fs_solid::hyper2d::plastic::thermal::{FiniteThermoplasticSheet, FiniteThermoplasticSheetStep};
use fs_solid::hyper3::plastic::thermal::FiniteThermoplasticElementEnergy;
use fs_solid::hyper3::plastic::{FinitePlasticLawField, FinitePlasticStepControls};
use fs_solid::mesh2::Mesh2;
use fs_solid::plastic3::IncrementSettings;

const REFERENCE_T: f64 = 300.0;
const THICKNESS: f64 = 0.01;

fn with_cx<R>(f: impl FnOnce(&Cx<'_>, &CancelGate) -> R) -> R {
    let gate = CancelGate::new();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        f(
            &Cx::new(
                &gate,
                arena,
                StreamKey {
                    seed: 47,
                    kernel_id: 176,
                    tile: 0,
                    iteration: 0,
                },
                Budget::INFINITE,
                ExecMode::Deterministic,
            ),
            &gate,
        )
    })
}

fn parameters() -> FiniteJ2ThermalParameters {
    FiniteJ2ThermalParameters {
        reference_temperature_k: REFERENCE_T,
        temperature_range_k: [250.0, 650.0],
        reference: FiniteJ2Parameters {
            shear_modulus_pa: 8000.0,
            bulk_modulus_pa: 16000.0,
            yield_stress_pa: 200.0,
            hardening_modulus_pa: 800.0,
        },
        slopes: FiniteJ2TemperatureSlopes {
            shear_modulus_pa_k: -6.0,
            bulk_modulus_pa_k: -8.0,
            yield_stress_pa_k: -0.4,
            hardening_modulus_pa_k: -1.0,
        },
        logarithmic_expansion_per_k: 3e-4,
        logarithmic_expansion_slope_per_k2: 1e-7,
    }
}

fn family(p: FiniteJ2ThermalParameters) -> AffineThermalFiniteJ2 {
    AffineThermalFiniteJ2::new(
        p,
        FiniteJ2Limits {
            minimum_jacobian: 0.25,
            maximum_elastic_log_strain: 0.7,
            maximum_plastic_log_strain: 2.0,
            maximum_accumulated_plastic_strain: 3.0,
        },
    )
    .unwrap()
}

fn budget() -> SheetAssemblyBudget {
    SheetAssemblyBudget {
        maximum_nodes: 6,
        maximum_elements: 2,
        maximum_free_dofs: 12,
        minimum_scaled_jacobian: 1e-4,
    }
}

fn controls() -> FinitePlasticStepControls {
    FinitePlasticStepControls {
        equilibrium_work_tolerance_j: 1e-7,
        maximum_endpoint_work_defect_j: 10.0,
        yield_tolerance_pa: 1e-5,
    }
}

fn settings() -> IncrementSettings {
    IncrementSettings {
        relative_tolerance: 1e-11,
        ..IncrementSettings::default()
    }
}

fn local() -> PlaneStressSettings {
    PlaneStressSettings {
        relative_tolerance: 1e-13,
        ..PlaneStressSettings::default()
    }
}

fn unit_mesh(triangles: bool) -> Mesh2 {
    if triangles {
        Mesh2::triangles(1.0, 1.0, 1, 1)
    } else {
        Mesh2::quads(1.0, 1.0, 1, 1)
    }
}

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "actual {actual:.16e}, expected {expected:.16e}, tolerance {tolerance:e}"
    );
}

// Hold axial stretch; the two upper-edge y coordinates solve lateral contraction.
fn axial_supports(mesh: &Mesh2, stretch: f64) -> Vec<(usize, f64)> {
    let mut result = Vec::new();
    for (n, &[x, y]) in mesh.nodes.iter().enumerate() {
        result.push((2 * n, (stretch - 1.0) * x));
        if y == 0.0 {
            result.push((2 * n + 1, 0.0));
        }
    }
    result
}

fn clamp(displacement: &[[f64; 2]]) -> Vec<(usize, f64)> {
    displacement
        .as_flattened()
        .iter()
        .copied()
        .enumerate()
        .collect()
}

// Uniform-temperature unit-square fixtures have equal reference-volume weights.
fn check_native_accounts(
    law: &AffineThermalFiniteJ2,
    old: &FiniteThermoplasticSheet<'_>,
    step: &FiniteThermoplasticSheetStep<'_>,
    new_temperature: f64,
    dt: f64,
) {
    let count = step.quadrature_energy.len();
    let per_element = count / old.temperatures_k().len();
    let reference_points: Vec<_> = (0..count)
        .map(|q| (q / per_element, THICKNESS / count as f64))
        .collect();
    check_native_field_accounts(
        FinitePlasticLawField::Uniform(law),
        old,
        step,
        &reference_points,
        dt,
    );
    assert!(
        step.sheet
            .temperatures_k()
            .iter()
            .all(|&t| t == new_temperature)
    );
}

// Recompute each account through its original 3-D thermal owner, using its
// element's accepted temperature, original volume and one physical duration.
// This catches temperature averaging, history advances inside Newton and
// accidental deposition of a kinetic potential as physical heat.
fn check_native_field_accounts(
    families: FinitePlasticLawField<'_, AffineThermalFiniteJ2>,
    old: &FiniteThermoplasticSheet<'_>,
    step: &FiniteThermoplasticSheetStep<'_>,
    reference_points: &[(usize, f64)],
    dt: f64,
) {
    let new = &step.sheet;
    let reports = &step.quadrature_energy;
    assert_eq!(reports.len(), old.quadrature_states().len());
    assert_eq!(reports.len(), new.quadrature_states().len());
    assert_eq!(reports.len(), reference_points.len());
    for (q, report) in reports.iter().enumerate() {
        let (element, _) = reference_points[q];
        let law = match families {
            FinitePlasticLawField::Uniform(law) => law,
            FinitePlasticLawField::PerElement(laws) => &laws[element],
        };
        let native = law
            .checked_increment_for_step(
                &old.quadrature_deformation()[q],
                &old.quadrature_states()[q],
                old.temperatures_k()[element],
                &new.quadrature_deformation()[q],
                new.temperatures_k()[element],
                controls().yield_tolerance_pa,
                dt,
            )
            .unwrap();
        close(
            new.quadrature_states()[q].alpha,
            native.mechanical.state.alpha,
            1e-10,
        );
        close(
            new.quadrature_states()[q].creep_alpha,
            native.mechanical.state.creep_alpha,
            1e-10,
        );
        for (a, b) in new.quadrature_states()[q]
            .inverse_plastic_cauchy_green
            .iter()
            .zip(native.mechanical.state.inverse_plastic_cauchy_green)
        {
            close(*a, b, 1e-9);
        }
        for (a, b) in new.quadrature_first_piola_pa()[q]
            .iter()
            .zip(native.mechanical.first_piola_pa)
        {
            close(*a, b, 2e-6);
        }
        for (actual, expected) in [
            (report.free_before_j_m3, native.old.free_energy_j_m3),
            (report.free_after_j_m3, native.new.free_energy_j_m3),
            (report.internal_before_j_m3, native.old.internal_energy_j_m3),
            (report.internal_after_j_m3, native.new.internal_energy_j_m3),
            (
                report.plastic_dissipation_j_m3,
                native.mechanical.plastic_dissipation_j_m3,
            ),
            (
                report.creep_dissipation_j_m3,
                native.mechanical.creep_dissipation_j_m3,
            ),
            (report.reversible_heat_j_m3, native.reversible_heat_j_m3),
            (
                report.coefficient_storage_change_j_m3,
                native.coefficient_storage_change_j_m3,
            ),
            (
                report.thermoelastic_storage_change_j_m3,
                native.thermoelastic_storage_change_j_m3,
            ),
            (
                report.endpoint_work_defect_j_m3,
                native.endpoint_work_defect_j_m3,
            ),
            (
                report.incremental_potential_j_m3,
                native.mechanical.incremental_potential_j_m3,
            ),
            (
                report.internal_endpoint_work_j_m3,
                native.mechanical_work_j_m3,
            ),
            (
                report.mechanical_coordinate_work_j_m3,
                native.mechanical_coordinate_work_j_m3,
            ),
            (
                report.thermal_deformation_work_j_m3,
                native.thermal_deformation_work_j_m3,
            ),
            (report.balance_error_j_m3, native.balance_error_j_m3),
        ] {
            close(actual, expected, 2e-6);
        }
        for transverse in [2, 5, 8] {
            close(new.quadrature_first_piola_pa()[q][transverse], 0.0, 2e-6);
        }
    }
    let integral = |field: fn(&FiniteThermoplasticElementEnergy) -> f64| {
        reports
            .iter()
            .zip(reference_points)
            .map(|(report, &(_, volume))| volume * field(report))
            .sum::<f64>()
    };
    close(new.free_energy_j(), integral(|r| r.free_after_j_m3), 1e-8);
    close(
        new.internal_energy_j(),
        integral(|r| r.internal_after_j_m3),
        1e-8,
    );
    for (actual, expected) in [
        (
            step.energy.free_energy_change_j,
            integral(|r| r.free_after_j_m3 - r.free_before_j_m3),
        ),
        (
            step.energy.internal_energy_change_j,
            integral(FiniteThermoplasticElementEnergy::internal_energy_change_j_m3),
        ),
        (
            step.energy.plastic_dissipation_j,
            integral(|r| r.plastic_dissipation_j_m3),
        ),
        (
            step.energy.creep_dissipation_j,
            integral(|r| r.creep_dissipation_j_m3),
        ),
        (
            step.energy.reversible_heat_j,
            integral(|r| r.reversible_heat_j_m3),
        ),
        (
            step.energy.caloric_heat_j,
            integral(FiniteThermoplasticElementEnergy::caloric_heat_j_m3),
        ),
        (
            step.energy.coefficient_storage_change_j,
            integral(|r| r.coefficient_storage_change_j_m3),
        ),
        (
            step.energy.thermoelastic_storage_change_j,
            integral(|r| r.thermoelastic_storage_change_j_m3),
        ),
        (
            step.energy.endpoint_work_defect_j,
            integral(|r| r.endpoint_work_defect_j_m3),
        ),
        (
            step.energy.incremental_potential_j,
            integral(|r| r.incremental_potential_j_m3),
        ),
        (
            step.energy.internal_endpoint_work_j,
            integral(|r| r.internal_endpoint_work_j_m3),
        ),
    ] {
        close(actual, expected, 1e-8);
    }
    let external = step.energy.dead_load_endpoint_work_j + step.energy.support_endpoint_work_j;
    close(
        external
            - step.energy.internal_energy_change_j
            - step.energy.caloric_heat_j
            - step.energy.endpoint_work_defect_j,
        0.0,
        controls().equilibrium_work_tolerance_j,
    );
    close(
        step.energy.equilibrium_work_residual_j,
        0.0,
        controls().equilibrium_work_tolerance_j,
    );
    close(
        step.energy.balance_error_j,
        0.0,
        controls().equilibrium_work_tolerance_j,
    );
    close(new.physical_time_s(), old.physical_time_s() + dt, 1e-12);
}

#[test]
fn free_thermal_expansion_and_cooling_solve_thickness_without_plastic_memory() {
    with_cx(|cx, _| {
        let law = family(parameters());
        for triangular in [false, true] {
            let mesh = unit_mesh(triangular);
            let thickness = vec![THICKNESS; mesh.elems.len()];
            let temperatures = vec![REFERENCE_T; mesh.elems.len()];
            let initial = FiniteThermoplasticSheet::try_new(
                &mesh,
                FinitePlasticLawField::Uniform(&law),
                &thickness,
                &temperatures,
                local(),
                budget(),
                cx,
            )
            .unwrap();
            let delta: f64 = 200.0;
            let theta = (3e-4 * delta + 0.5e-7 * delta * delta).exp();
            // Prescribe stress-free motion except at the top-right corner, whose
            // two coordinates actually equilibrate. All reactions must vanish.
            let supports: Vec<_> = mesh
                .nodes
                .iter()
                .enumerate()
                .filter(|(_, p)| **p != [1.0, 1.0])
                .flat_map(|(n, p)| {
                    [
                        (2 * n, (theta - 1.0) * p[0]),
                        (2 * n + 1, (theta - 1.0) * p[1]),
                    ]
                })
                .collect();
            let step = initial
                .advance(
                    &vec![500.0; mesh.elems.len()],
                    &supports,
                    &[],
                    0.2,
                    settings(),
                    controls(),
                    cx,
                )
                .unwrap();
            check_native_accounts(&law, &initial, &step, 500.0, 0.2);
            assert!(step.equilibrium.residual_history.len() > 1);
            for (x, u) in mesh.nodes.iter().zip(step.sheet.displacement_m()) {
                for i in 0..2 {
                    close(u[i], (theta - 1.0) * x[i], 2e-9);
                }
            }
            for &reaction in &step.equilibrium.reactions_n {
                close(reaction, 0.0, 2e-7);
            }
            for (f, h) in step
                .sheet
                .quadrature_deformation()
                .iter()
                .zip(step.sheet.quadrature_thickness_m())
            {
                for i in [0, 4, 8] {
                    close(f[i], theta, 2e-9);
                }
                close(h, THICKNESS * theta, 2e-11);
            }
            assert!(
                step.sheet
                    .quadrature_states()
                    .iter()
                    .all(|s| *s == law.initial_state())
            );
            close(step.sheet.free_energy_j(), 0.0, 1e-12);
            close(step.energy.plastic_dissipation_j, 0.0, 0.0);
            let FinitePlasticLawField::Uniform(retained) = step.sheet.families() else {
                panic!("the original uniform thermal family was replaced");
            };
            assert!(std::ptr::eq(retained, &law));
            let cool_supports: Vec<_> = supports.iter().map(|&(dof, _)| (dof, 0.0)).collect();
            let cooled = step
                .sheet
                .advance(
                    &temperatures,
                    &cool_supports,
                    &[],
                    0.3,
                    settings(),
                    controls(),
                    cx,
                )
                .unwrap();
            check_native_accounts(&law, &step.sheet, &cooled, REFERENCE_T, 0.3);
            for u in cooled.sheet.displacement_m() {
                for &v in u {
                    close(v, 0.0, 2e-9);
                }
            }
            for h in cooled.sheet.quadrature_thickness_m() {
                close(h, THICKNESS, 2e-11);
            }
            assert!(
                cooled
                    .sheet
                    .quadrature_states()
                    .iter()
                    .all(|s| *s == law.initial_state())
            );
            close(initial.physical_time_s(), 0.0, 0.0);
        }
    });
}

#[test]
fn heating_a_clamped_preloaded_sheet_rechecks_history_at_old_temperature_and_yields() {
    with_cx(|cx, _| {
        let mut p = parameters();
        p.logarithmic_expansion_per_k = 0.0;
        p.logarithmic_expansion_slope_per_k2 = 0.0;
        let law = family(p);
        for triangular in [false, true] {
            let mesh = unit_mesh(triangular);
            let thickness = vec![THICKNESS; mesh.elems.len()];
            let temperatures = vec![REFERENCE_T; mesh.elems.len()];
            let initial = FiniteThermoplasticSheet::try_new(
                &mesh,
                FinitePlasticLawField::Uniform(&law),
                &thickness,
                &temperatures,
                local(),
                budget(),
                cx,
            )
            .unwrap();
            let loaded = initial
                .advance(
                    &temperatures,
                    &axial_supports(&mesh, 1.12),
                    &[],
                    0.1,
                    settings(),
                    controls(),
                    cx,
                )
                .unwrap();
            let old = &loaded.sheet;
            assert!(old.quadrature_states().iter().all(|s| s.alpha > 0.05));
            for (f, state) in old
                .quadrature_deformation()
                .iter()
                .zip(old.quadrature_states())
            {
                law.validate_committed(f, state, REFERENCE_T, controls().yield_tolerance_pa)
                    .unwrap();
                assert!(
                    law.validate_committed(f, state, 600.0, controls().yield_tolerance_pa)
                        .is_err()
                );
            }
            let before = old.quadrature_states().to_vec();
            let heated = old
                .advance(
                    &vec![600.0; mesh.elems.len()],
                    &clamp(old.displacement_m()),
                    &[],
                    0.2,
                    settings(),
                    controls(),
                    cx,
                )
                .unwrap();
            check_native_accounts(&law, old, &heated, 600.0, 0.2);
            assert_eq!(heated.sheet.displacement_m(), old.displacement_m());
            for q in 0..before.len() {
                assert!(heated.sheet.quadrature_states()[q].alpha > before[q].alpha);
                assert!(
                    heated.sheet.quadrature_first_piola_pa()[q][0]
                        < old.quadrature_first_piola_pa()[q][0]
                );
            }
            assert!(heated.energy.plastic_dissipation_j > 0.0);
            assert!(heated.energy.coefficient_storage_change_j < 0.0);
            close(heated.energy.support_endpoint_work_j, 0.0, 1e-12);
            close(heated.energy.dead_load_endpoint_work_j, 0.0, 0.0);
            assert_eq!(old.quadrature_states(), before);
            assert_eq!(old.temperatures_k(), temperatures);
        }

        // Adjacent clamped cells receive independent temperatures and original
        // thicknesses. Heating one cannot contaminate the cold cell's history.
        let mesh = Mesh2::quads(1.0, 1.0, 2, 1);
        let mut neighbor = parameters();
        neighbor.reference_temperature_k = 310.0;
        let families = [family(parameters()), family(neighbor)];
        let thickness = [THICKNESS, 2.0 * THICKNESS];
        let initial = FiniteThermoplasticSheet::try_new(
            &mesh,
            FinitePlasticLawField::PerElement(&families),
            &thickness,
            &[REFERENCE_T, 310.0],
            local(),
            budget(),
            cx,
        )
        .unwrap();
        let step = initial
            .advance(
                &[500.0, 310.0],
                &clamp(&vec![[0.0; 2]; mesh.nodes.len()]),
                &[],
                0.1,
                settings(),
                controls(),
                cx,
            )
            .unwrap();
        assert_eq!(step.sheet.temperatures_k(), &[500.0, 310.0]);
        let points: Vec<_> = (0..8)
            .map(|q| {
                if q < 4 {
                    (0, THICKNESS / 8.0)
                } else {
                    (1, THICKNESS / 4.0)
                }
            })
            .collect();
        check_native_field_accounts(
            FinitePlasticLawField::PerElement(&families),
            &initial,
            &step,
            &points,
            0.1,
        );
        for q in 0..4 {
            assert!(step.sheet.quadrature_states()[q].alpha > 0.0);
            assert!(step.sheet.quadrature_first_piola_pa()[q][0] < 0.0);
            assert!(step.sheet.quadrature_first_piola_pa()[q][4] < 0.0);
        }
        for q in 4..8 {
            assert_eq!(
                step.sheet.quadrature_states()[q],
                families[1].initial_state()
            );
            for &stress in &step.sheet.quadrature_first_piola_pa()[q] {
                close(stress, 0.0, 1e-8);
            }
            close(step.sheet.quadrature_thickness_m()[q], thickness[1], 1e-11);
        }
        let FinitePlasticLawField::PerElement(retained) = step.sheet.families() else {
            panic!("element thermal families were homogenized");
        };
        assert!(std::ptr::eq(retained, families.as_slice()));
        assert_eq!(initial.temperatures_k(), &[REFERENCE_T, 310.0]);
    });
}

#[test]
fn arrhenius_sheet_relaxation_uses_actual_duration_and_exact_native_heat_at_each_point() {
    with_cx(|cx, _| {
        let mut p = parameters();
        p.slopes = FiniteJ2TemperatureSlopes {
            shear_modulus_pa_k: 0.0,
            bulk_modulus_pa_k: 0.0,
            yield_stress_pa_k: 0.0,
            hardening_modulus_pa_k: 0.0,
        };
        p.logarithmic_expansion_per_k = 0.0;
        p.logarithmic_expansion_slope_per_k2 = 0.0;
        let law = family(p)
            .with_overstress(FiniteJ2ThermalOverstress {
                reference_viscosity_pa_s: 4000.0,
                stress_scale_pa: 100.0,
                exponent: 1.0,
                activation_temperature_k: 3000.0,
            })
            .unwrap();
        for triangular in [false, true] {
            let mesh = unit_mesh(triangular);
            let thickness = vec![THICKNESS; mesh.elems.len()];
            let temperatures = vec![REFERENCE_T; mesh.elems.len()];
            let initial = FiniteThermoplasticSheet::try_new(
                &mesh,
                FinitePlasticLawField::Uniform(&law),
                &thickness,
                &temperatures,
                local(),
                budget(),
                cx,
            )
            .unwrap();
            let supports = axial_supports(&mesh, 1.16);
            let loaded = initial
                .advance(
                    &temperatures,
                    &supports,
                    &[],
                    0.04,
                    settings(),
                    controls(),
                    cx,
                )
                .unwrap();
            check_native_accounts(&law, &initial, &loaded, REFERENCE_T, 0.04);
            assert!(loaded.equilibrium.residual_history.len() > 1);
            let old = &loaded.sheet;
            let before = old.quadrature_states().to_vec();
            let mut flow = Vec::new();
            for (dt, temperature) in [(0.02, REFERENCE_T), (0.4, REFERENCE_T), (0.02, 420.0)] {
                let step = old
                    .advance(
                        &vec![temperature; mesh.elems.len()],
                        &supports,
                        &[],
                        dt,
                        settings(),
                        controls(),
                        cx,
                    )
                    .unwrap();
                check_native_accounts(&law, old, &step, temperature, dt);
                let increment: f64 = step
                    .sheet
                    .quadrature_states()
                    .iter()
                    .zip(&before)
                    .map(|(new, old)| new.alpha - old.alpha)
                    .sum();
                assert!(increment > 0.0);
                for (new, old) in step
                    .sheet
                    .quadrature_first_piola_pa()
                    .iter()
                    .zip(old.quadrature_first_piola_pa())
                {
                    assert!(new[0] < old[0]);
                }
                assert!(step.energy.plastic_dissipation_j > 0.0);
                close(step.energy.reversible_heat_j, 0.0, 1e-10);
                close(
                    step.energy.caloric_heat_j,
                    step.energy.plastic_dissipation_j,
                    1e-10,
                );
                assert!(
                    step.energy.incremental_potential_j
                        < step.sheet.free_energy_j() + step.energy.plastic_dissipation_j - 1e-8
                );
                flow.push(increment);
            }
            assert!(
                flow[1] > flow[0],
                "a longer physical hold must relax further"
            );
            assert!(
                flow[2] > flow[0],
                "Arrhenius activation must accelerate the same hold"
            );
            assert_eq!(old.quadrature_states(), before);
            assert_eq!(old.temperatures_k(), temperatures);
            close(old.physical_time_s(), 0.04, 1e-12);
        }
    });
}

#[test]
fn reference_temperatures_refusals_and_cancellation_preserve_every_accepted_sheet_field() {
    with_cx(|cx, gate| {
        let mesh = unit_mesh(false);
        let law = family(parameters());
        let thickness = [THICKNESS];
        for temperatures in [vec![], vec![301.0], vec![f64::NAN]] {
            assert!(
                FiniteThermoplasticSheet::try_new(
                    &mesh,
                    FinitePlasticLawField::Uniform(&law),
                    &thickness,
                    &temperatures,
                    local(),
                    budget(),
                    cx
                )
                .is_err()
            );
        }
        // Each family has its own stress-free reference: a common mesh-wide
        // reference cannot silently replace a heterogeneous source binding.
        let triangle_mesh = unit_mesh(true);
        let mut second = parameters();
        second.reference_temperature_k = 310.0;
        let families = [family(parameters()), family(second)];
        let two_thicknesses = [THICKNESS; 2];
        assert!(
            FiniteThermoplasticSheet::try_new(
                &triangle_mesh,
                FinitePlasticLawField::PerElement(&families),
                &two_thicknesses,
                &[300.0, 300.0],
                local(),
                budget(),
                cx
            )
            .is_err()
        );
        assert!(
            FiniteThermoplasticSheet::try_new(
                &triangle_mesh,
                FinitePlasticLawField::PerElement(&families),
                &two_thicknesses,
                &[300.0, 310.0],
                local(),
                budget(),
                cx
            )
            .is_ok()
        );
        let initial = FiniteThermoplasticSheet::try_new(
            &mesh,
            FinitePlasticLawField::Uniform(&law),
            &thickness,
            &[REFERENCE_T],
            local(),
            budget(),
            cx,
        )
        .unwrap();
        let sheet = initial
            .advance(
                &[350.0],
                &axial_supports(&mesh, 1.08),
                &[],
                0.1,
                settings(),
                controls(),
                cx,
            )
            .unwrap()
            .sheet;
        let u = sheet.displacement_m().to_vec();
        let states = sheet.quadrature_states().to_vec();
        let f = sheet.quadrature_deformation().to_vec();
        let piola = sheet.quadrature_first_piola_pa().to_vec();
        let h = sheet.quadrature_thickness_m();
        let free = sheet.free_energy_j();
        let internal = sheet.internal_energy_j();
        let unchanged = || {
            assert_eq!(sheet.displacement_m(), u);
            assert_eq!(sheet.quadrature_states(), states);
            assert_eq!(sheet.quadrature_deformation(), f);
            assert_eq!(sheet.quadrature_first_piola_pa(), piola);
            assert_eq!(sheet.quadrature_thickness_m(), h);
            assert_eq!(sheet.temperatures_k(), &[350.0]);
            assert_eq!(sheet.free_energy_j(), free);
            assert_eq!(sheet.internal_energy_j(), internal);
            assert_eq!(sheet.physical_time_s(), 0.1);
        };
        let supports = clamp(&u);
        for temperatures in [vec![], vec![700.0], vec![f64::NAN]] {
            assert!(
                sheet
                    .advance(
                        &temperatures,
                        &supports,
                        &[],
                        0.1,
                        settings(),
                        controls(),
                        cx
                    )
                    .is_err()
            );
            unchanged();
        }
        for dt in [0.0, -0.1, f64::NAN, f64::INFINITY] {
            assert!(
                sheet
                    .advance(&[350.0], &supports, &[], dt, settings(), controls(), cx)
                    .is_err()
            );
            unchanged();
        }
        gate.request();
        assert!(
            sheet
                .advance(&[400.0], &supports, &[], 0.1, settings(), controls(), cx)
                .is_err()
        );
        unchanged();
    });
}
