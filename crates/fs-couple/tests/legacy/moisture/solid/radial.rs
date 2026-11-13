//! Authored radial-to-native-solid checks; these are not measured wood data.
use super::*;
use fs_couple::moisture::radial::{MoistureRadial, MoistureRadialGeometry, RadialMoistureSymmetry};

#[path = "radial/solution.rs"]
mod finite_solution;

fn radial_specimen(
    material: MaterialCard,
    radii: Vec<f64>,
    contents: Vec<f64>,
    cx: &Cx<'_>,
) -> MoistureRadial {
    MoistureRadial::new(
        material,
        point(0.1),
        "moisture_diffusivity",
        MaterialPropertySelection::SingleClaimOnly,
        MoistureRadialGeometry {
            symmetry: RadialMoistureSymmetry::Cylinder { length_m: 0.1 },
            faces_m: radii,
            maximum_moisture: 0.3,
            maximum_cells: 32,
        },
        contents,
        cx,
    )
    .unwrap()
}

fn resolution() -> MoistureRadialSolidMesh {
    MoistureRadialSolidMesh {
        angular_divisions: 8,
        axial_divisions: 1,
        maximum_relative_shell_volume_error: 0.11,
    }
}

fn nonlinear() -> MoistureNonlinearConfig {
    MoistureNonlinearConfig {
        relative_tolerance: 1e-12,
        maximum_iterations: 128,
        relaxation: 1.0,
    }
}

fn radial_binding(
    moisture: &MoistureRadial,
    fixed: Vec<usize>,
    cx: &Cx<'_>,
) -> MoistureRadialSolid {
    MoistureRadialSolid::new(
        moisture,
        resolution(),
        SolidStiffnessModel::Isotropic,
        MaterialPropertySelection::SingleClaimOnly,
        MaterialPropertySelection::SingleClaimOnly,
        axes(),
        fixed,
        controls(),
        cx,
    )
    .unwrap()
}

#[test]
fn radial_uniform_swelling_reproduces_affine_shape_and_exact_curved_wet_inventory() {
    with_cx(|cx| {
        let radii = vec![0.0, 0.02, 0.04];
        let dry = radial_specimen(card([0.01; 2], 0.3), radii.clone(), vec![0.1; 2], cx);
        let wet = radial_specimen(card([0.01; 2], 0.3), radii, vec![0.2; 2], cx);
        // Eliminate rigid motion without restraining diagonal free strain:
        // center xyz, first +x ring node yz, and first +y ring node z.
        let binding = radial_binding(&dry, vec![0, 1, 2, 4, 5, 11], cx);
        let state = binding.resolve(&wet, cx).unwrap();
        let strain = [0.001, 0.004, 0.002];
        for (x, u) in binding
            .reference_nodes_m()
            .iter()
            .zip(state.solution().displacement_m())
        {
            for k in 0..3 {
                near(u[k], strain[k] * x[k], 1e-10);
            }
        }
        for element in state.stress().elements() {
            for stress in element.stress_mandel_pa() {
                near(*stress, 0.0, 0.1);
            }
        }
        let curved_volume = core::f64::consts::PI * 0.04_f64.powi(2) * 0.1;
        let polygon_fraction = 8.0 * core::f64::consts::FRAC_PI_4.sin() / core::f64::consts::TAU;
        let expected_volume =
            curved_volume * polygon_fraction * strain.iter().map(|s| 1.0 + s).product::<f64>();
        near(
            state.geometry().element_volumes_m3().iter().sum(),
            expected_volume,
            1e-13,
        );
        near(
            binding.maximum_relative_shell_volume_error(),
            1.0 - polygon_fraction,
            2e-14,
        );
        // The polygonal tessellation does not lose the curved shell's water.
        near(state.wet_mass_kg(), 500.0 * curved_volume * 1.2, 1e-12);
        near(
            state.tangent().unwrap().pencil().total_mass_kg,
            wet.wet_mass_kg(),
            1e-12,
        );
        near(
            state
                .current_density_kg_m3()
                .iter()
                .zip(state.geometry().element_volumes_m3())
                .map(|(rho, volume)| rho * volume)
                .sum(),
            wet.wet_mass_kg(),
            1e-12,
        );
        assert_eq!(state.elastic_states().len(), 2);
        assert_eq!(state.materials().len(), binding.tetrahedra().len());
    });
}

#[test]
fn radial_diffusion_gradient_changes_compatible_stress_and_native_modal_pencil() {
    with_cx(|cx| {
        let moisture = radial_specimen(
            card([0.01; 2], 0.3),
            vec![0.01, 0.02, 0.04],
            vec![0.1; 2],
            cx,
        );
        // Three rings x eight angles x two planes. One actual generalized DOF
        // gives an independent exact K00/M00 eigenvalue, without mode fitting.
        let free = 141;
        let binding = radial_binding(&moisture, (0..144).filter(|d| *d != free).collect(), cx);
        let before = binding.resolve(&moisture, cx).unwrap();
        let step = binding
            .advance(
                &moisture,
                0.01,
                &[MoistureFace::Sealed, MoistureFace::PrescribedDryBasis(0.2)],
                config(),
                nonlinear(),
                cx,
            )
            .unwrap();
        assert!(step.moisture.moisture()[1] > step.moisture.moisture()[0]);
        assert!(step.moisture.moisture()[0] > 0.1);
        let eigenvalue = |state: &MoistureSolidState| {
            let p = state.tangent().unwrap().pencil();
            assert_eq!(p.free_dofs.as_slice(), &[free]);
            p.stiffness.row(0).1[0] / p.mass.row(0).1[0]
        };
        assert!(eigenvalue(&step.solid) > 0.0);
        assert!(eigenvalue(&step.solid) < eigenvalue(&before));
        assert!(step.solid.stress().elastic_energy_j() > 0.0);
        near(step.solid.wet_mass_kg(), step.moisture.wet_mass_kg(), 1e-12);
        near(
            step.solid.tangent().unwrap().pencil().total_mass_kg,
            step.moisture.wet_mass_kg(),
            1e-12,
        );
        assert_eq!(step.solid.time_s(), step.moisture.time_s());
        assert_eq!(moisture.time_s(), 0.0);
        assert_eq!(moisture.moisture(), &[0.1; 2]);
    });
}

#[test]
fn radial_solid_refuses_unresolved_geometry_foreign_sources_and_late_material_extrapolation() {
    with_cx(|cx| {
        let limited = radial_specimen(card([0.01; 2], 0.12), vec![0.0, 0.04], vec![0.1], cx);
        let fixed = vec![0, 1, 2, 4, 5, 11];
        let binding = radial_binding(&limited, fixed.clone(), cx);
        let mut fine_geometry_required = resolution();
        fine_geometry_required.maximum_relative_shell_volume_error = 0.01;
        assert!(
            MoistureRadialSolid::new(
                &limited,
                fine_geometry_required,
                SolidStiffnessModel::Isotropic,
                MaterialPropertySelection::SingleClaimOnly,
                MaterialPropertySelection::SingleClaimOnly,
                axes(),
                fixed,
                controls(),
                cx,
            )
            .is_err()
        );
        let foreign = radial_specimen(card([0.01; 2], 0.3), vec![0.0, 0.05], vec![0.1], cx);
        assert!(binding.resolve(&foreign, cx).is_err());
        let boundary = [MoistureFace::Sealed, MoistureFace::PrescribedDryBasis(0.2)];
        let water_only = limited
            .advance(0.5, &boundary, config(), nonlinear(), cx)
            .unwrap()
            .0;
        assert!(water_only.moisture()[0] > 0.12);
        assert!(
            binding
                .advance(&limited, 0.5, &boundary, config(), nonlinear(), cx)
                .is_err()
        );
        assert_eq!(limited.moisture(), &[0.1]);
        assert_eq!(limited.time_s(), 0.0);
    });
}

#[test]
fn finite_radial_vapor_drives_native_shape_and_preserves_both_inputs_on_late_refusal() {
    use crate::atmosphere::{atmosphere as gas, controls as air_controls, sorption};
    with_cx(|cx| {
        let air = gas(0.8, 1.0);
        let saved_air = air.clone();
        let original = radial_specimen(card([0.01; 2], 0.3), vec![0.0, 0.04], vec![0.1], cx);
        let binding = radial_binding(&original, vec![0, 1, 2, 4, 5, 11], cx);
        let boundary = [None, Some(sorption(0.25))];
        let (candidate, solid) = binding
            .advance_finite_atmosphere(
                &air,
                &original,
                0.25,
                &boundary,
                config(),
                nonlinear(),
                air_controls(),
                cx,
            )
            .unwrap();
        assert!(candidate.step().water_to_solid_kg > 0.0);
        assert!(
            candidate.atmosphere().relative_humidity().unwrap() < air.relative_humidity().unwrap()
        );
        near(
            candidate.atmosphere().gas().mass_kg() + solid.wet_mass_kg(),
            air.gas().mass_kg() + original.wet_mass_kg(),
            3e-12,
        );
        near(
            solid.tangent().unwrap().pencil().total_mass_kg,
            candidate.moisture().wet_mass_kg(),
            1e-12,
        );
        let change = candidate.moisture().moisture()[0] - original.moisture()[0];
        for (x, u) in binding
            .reference_nodes_m()
            .iter()
            .zip(solid.solution().displacement_m())
        {
            for k in 0..3 {
                near(u[k], [0.01, 0.04, 0.02][k] * change * x[k], 1e-10);
            }
        }
        assert_eq!(solid.time_s(), candidate.moisture().time_s());
        let limited = radial_specimen(card([0.01; 2], 0.105), vec![0.0, 0.04], vec![0.1], cx);
        let limited_binding = radial_binding(&limited, vec![0, 1, 2, 4, 5, 11], cx);
        let water_only = air
            .advance_radial(
                &limited,
                0.25,
                &boundary,
                config(),
                nonlinear(),
                air_controls(),
                cx,
            )
            .unwrap();
        assert!(water_only.moisture().moisture()[0] > 0.105);
        assert!(
            limited_binding
                .advance_finite_atmosphere(
                    &air,
                    &limited,
                    0.25,
                    &boundary,
                    config(),
                    nonlinear(),
                    air_controls(),
                    cx,
                )
                .is_err()
        );
        assert_eq!(air, saved_air);
        assert_eq!(original.moisture(), &[0.1]);
        assert_eq!(original.time_s(), 0.0);
        assert_eq!(limited.moisture(), &[0.1]);
        assert_eq!(limited.time_s(), 0.0);
    });
}
