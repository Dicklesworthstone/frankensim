//! Forward contact steps at source knots, without weakening physical adjoints.
use super::*;

fn corner_card() -> InterfaceSystemCard {
    card_with_knots(vec![(256.0, 0.25), (320.0, 0.5), (384.0, 1.0)])
}

#[test]
fn g1_equilibrium_source_kink_has_the_exact_two_sided_contact_flux_tangent() {
    context(&CancelGate::new_clock_free(), |cx| {
        let mesh = mesh();
        let bc = ThermalBoundaryBuilder::new(&mesh).adiabatic_remainder().finish().unwrap();
        let card = corner_card();
        let contacts = contacts(cx, &mesh, &bc, &card, false);
        let temperature = [320.0; 8];
        assert!(contacts.linearize_temperature(cx, &bc, &temperature, 1).is_err());
        let linear = contacts.linearize_temperature_for_newton(cx, &bc, &temperature, 1).unwrap();
        let direction = [1.0, -2.0, 3.0, 0.0, -0.5, 1.0, -1.0, 0.0];
        let tangent = linear.apply(cx, &direction).unwrap();
        // Exact area=1/2 consistent triangle mass: M_ij=(1+delta_ij)/24.
        // At equilibrium d(1/R)/dT multiplies zero jump, so the full physical
        // flux has this unique derivative even though R itself has a corner.
        let jumps: [f64; 3] = std::array::from_fn(|i| direction[i] - direction[i + 4]);
        let mut expected = [0.0; 8];
        for i in 0..3 {
            expected[i] = (jumps.iter().sum::<f64>() + jumps[i]) / (24.0 * 0.5);
            expected[i + 4] = -expected[i];
        }
        for i in 0..8 {
            near(linear.nodal_outflow_w()[i], 0.0, 1e-12);
            near(tangent[i], expected[i], 1e-13);
        }
        near(tangent.iter().sum(), 0.0, 1e-13);
        for epsilon in [-1e-5, 1e-5] {
            let trial: [f64; 8] = std::array::from_fn(|i| temperature[i] + epsilon * direction[i]);
            let fresh = contacts.linearize_temperature_for_newton(cx, &bc, &trial, 1).unwrap();
            for i in 0..8 {
                near((fresh.nodal_outflow_w()[i] - linear.nodal_outflow_w()[i]) / epsilon,
                    expected[i], 1e-7);
            }
        }
    });
}

#[test]
fn g3_primal_contact_slopes_use_declared_segments_and_refuse_extrapolation() {
    context(&CancelGate::new_clock_free(), |cx| {
        let mesh = mesh();
        let bc = ThermalBoundaryBuilder::new(&mesh).adiabatic_remainder().finish().unwrap();
        let card = corner_card();
        let contacts = contacts(cx, &mesh, &bc, &card, false);
        for mean in [256.0, 319.0, 320.0, 321.0, 384.0] {
            let mut temperature = [mean - 10.0; 8];
            temperature[4..].fill(mean + 10.0);
            let linear = contacts.linearize_temperature_for_newton(cx, &bc, &temperature, 1).unwrap();
            let (resistance, slope) = if mean < 320.0 {
                (0.25 + (mean - 256.0) / 256.0, 1.0 / 256.0)
            } else {
                (0.5 + (mean - 320.0) / 128.0, 1.0 / 128.0)
            };
            // A uniform shift leaves every temperature jump unchanged. Only
            // the independently known source slope changes the contact load.
            let expected = slope / resistance * 20.0 / (6.0 * resistance);
            let action = linear.apply(cx, &[1.0; 8]).unwrap();
            for (actual, expected) in action.iter().zip([
                expected, expected, expected, 0.0, -expected, -expected, -expected, 0.0,
            ]) {
                near(*actual, expected, 2e-13);
            }
            if [256.0, 320.0, 384.0].contains(&mean) {
                assert!(contacts.linearize_temperature(cx, &bc, &temperature, 1).is_err());
            }
        }
        for temperature in [255.0, 385.0] {
            assert!(contacts.linearize_temperature_for_newton(cx, &bc, &[temperature; 8], 1).is_err());
        }
    });
}

#[test]
fn g0_primal_contact_keeps_surface_gradient_corners_strict() {
    context(&CancelGate::new_clock_free(), |cx| {
        let mesh = mesh();
        let bc = ThermalBoundaryBuilder::new(&mesh).adiabatic_remainder().finish().unwrap();
        let card = card_with_response(PropertyValue::Surface2D {
            axes: ["TA".into(), "TB".into()], axis_dims: [fs_conduction::TEMPERATURE_DIMS; 2],
            samples: vec![([256.0, 256.0], 0.25), ([384.0, 256.0], 0.5),
                ([256.0, 384.0], 0.5), ([384.0, 384.0], 1.0)],
            triangles: vec![[0, 1, 2], [1, 3, 2]],
            dims: fs_conduction::AREA_SPECIFIC_THERMAL_RESISTANCE_DIMS,
        }, two_sided_validity());
        let contacts = two_sided_contacts(cx, &mesh, &bc, &card, &two_sided_point());
        assert!(contacts.evaluate(cx, &bc, &[320.0; 8]).is_ok());
        assert!(contacts.linearize_temperature_for_newton(cx, &bc, &[320.0; 8], 1).is_err());
    });
}

#[test]
fn g1_capacity_newton_admits_source_corners_and_inward_ends_without_extrapolation() {
    use fs_conduction::{ConductionProblem, ConductivityModel, LinearConfig, ScalarField};
    use fs_conduction::interface::temperature::ContactIterationConfig;
    use fs_conduction::transient::VolumetricHeatCapacity;
    use fs_conduction::transient::backward_euler::{BackwardEuler, NonlinearStepConfig, StepConfig};

    context(&CancelGate::new_clock_free(), |cx| {
        let mesh = mesh();
        let bc = ThermalBoundaryBuilder::new(&mesh).adiabatic_remainder().finish().unwrap();
        let card = corner_card();
        let contacts = contacts(cx, &mesh, &bc, &card, false);
        let capacities = [VolumetricHeatCapacity::declared(10.0).unwrap(); 2];
        let engine = BackwardEuler::per_element(cx, &mesh, &capacities).unwrap();
        // The bulk material has no temperature-domain restriction. Any outward
        // refusal below must therefore come from the original contact source.
        let material = ConductivityModel::isotropic_declared(2.0).unwrap();
        let config = StepConfig {
            linear: LinearConfig { tolerance: 1e-12, max_iterations: 256, restart: 8 },
            energy_tolerance_j: 1e-8,
        };
        let nonlinear = NonlinearStepConfig {
            residual_rtol: 1e-12, residual_atol_j: 1e-11, ..Default::default()
        };
        let coupling = ContactIterationConfig {
            max_iterations: 1, relaxation: 0.5, residual_tolerance_j: 1e-8,
        };
        for (initial_temperature, volumetric_heat, expected_resistance) in [
            (320.0, 10.0, 0.5078125), (320.0, -10.0, 0.49609375),
            (256.0, 10.0, 0.25390625), (384.0, -10.0, 0.9921875),
        ] {
            let source = ScalarField::Uniform(volumetric_heat);
            let problem = ConductionProblem { mesh: &mesh, boundary: &bc, material: &material,
                element_materials: None, source: &source };
            let old = [initial_temperature; 8];
            let before = old.map(f64::to_bits);
            let initial = contacts.evaluate(cx, &bc, &old).unwrap();
            assert!(contacts.linearize_temperature(cx, &bc, &old, 1).is_err(),
                "forward admission must not weaken the strict endpoint derivative");
            let next = contacts.advance_capacity_newton(
                cx, &engine, problem, &old, 1.0, config, nonlinear, coupling, 1,
            ).unwrap();
            let expected_temperature = initial_temperature + volumetric_heat / 10.0;
            for &temperature in &next.conduction.temperature {
                near(temperature, expected_temperature, 2e-9);
            }
            assert_eq!(next.iterations, 1, "one complete coupled thermal response");
            assert!(next.nonlinear_updates > 0 && next.krylov_iterations > 0);
            // The original two unit tetrahedra have total volume 1/3. Their
            // uniform capacity/source gives delta T = dt*f/c_v exactly; no
            // contact heat is external input and no trial advances old storage.
            near(next.physical.source_w, volumetric_heat / 3.0, 1e-12);
            near(next.physical.stored_energy_change_j, volumetric_heat / 3.0, 1e-8);
            near(next.physical.energy_residual_j, 0.0, config.energy_tolerance_j);
            assert!(next.physical.residual_norm_j <= coupling.residual_tolerance_j);
            assert_eq!(next.physical_contacts.states(), next.applied_contacts.states());
            for state in next.physical_contacts.states() {
                near(state.interface_temperature_k, expected_temperature, 2e-9);
                near(state.resistance.value_m2_k_per_w(), expected_resistance, 2e-11);
                assert_eq!(state.resistance.card_identity(), card.content_hash());
                assert_eq!(state.resistance.receipt().selected,
                    initial.states()[0].resistance.receipt().selected);
            }
            assert_eq!(old.map(f64::to_bits), before);
        }
        for (initial_temperature, outward_heat) in [(256.0, -10.0), (384.0, 10.0)] {
            let old = [initial_temperature; 8];
            let before = old.map(f64::to_bits);
            let source = ScalarField::Uniform(outward_heat);
            let problem = ConductionProblem { mesh: &mesh, boundary: &bc, material: &material,
                element_materials: None, source: &source };
            assert!(contacts.evaluate(cx, &bc, &old).is_ok());
            assert!(matches!(contacts.advance_capacity_newton(
                cx, &engine, problem, &old, 1.0, config, nonlinear, coupling, 1,
            ), Err(ConductionError::LineSearchFailed { iteration: 0, backtracks, .. })
                if backtracks == nonlinear.line_search.max_backtracks),
                "every outward trial must remain outside the pinned source domain");
            assert_eq!(old.map(f64::to_bits), before);

            // The refused attempt leaves the same original history usable for
            // an inward retry through the unchanged physical owners.
            let inward = ScalarField::Uniform(-outward_heat);
            let retry = contacts.advance_capacity_newton(
                cx, &engine, ConductionProblem { source: &inward, ..problem }, &old,
                1.0, config, nonlinear, coupling, 1,
            ).unwrap();
            for &temperature in &retry.conduction.temperature {
                near(temperature, initial_temperature - outward_heat / 10.0, 2e-9);
            }
            assert_eq!(old.map(f64::to_bits), before);
        }
    });
}

#[test]
fn g1_enthalpy_newton_heats_and_cools_from_a_source_corner_and_inward_from_its_ends() {
    use fs_blake3::ContentHash;
    use fs_conduction::{ConductionProblem, ConductivityModel, ScalarField};
    use fs_conduction::interface::temperature::enthalpy::ContactEnthalpyNewtonConfig;
    use fs_conduction::transient::enthalpy::{EnthalpyBudget, EnthalpyStepConfig};
    use fs_conduction::transient::enthalpy::heterogeneous::{
        HeterogeneousEnthalpyBackwardEuler, ReferenceEnthalpyMaterial,
    };
    use fs_material::phase::{EnthalpyPhaseKnot, EquilibriumEnthalpyPhaseCurve};
    use fs_solver::NewtonKrylovConfig;

    context(&CancelGate::new_clock_free(), |cx| {
        let mesh = mesh();
        let bc = ThermalBoundaryBuilder::new(&mesh).adiabatic_remainder().finish().unwrap();
        let card = corner_card();
        let contacts = contacts(cx, &mesh, &bc, &card, false);
        // The caloric chart extends beyond contact support, so an outward
        // trial must be refused by the pinned contact source, not this chart.
        let chart = EquilibriumEnthalpyPhaseCurve::try_new(ContentHash([181; 32]),
            [(0.0, 240.0), (1600.0, 400.0)].into_iter().map(|(h, temperature)| EnthalpyPhaseKnot {
                specific_enthalpy_j_kg: h, temperature_k: temperature,
                liquid_mass_fraction: 0.0, bulk_density_kg_m3: 1.0,
            }).collect()).unwrap();
        let engine = HeterogeneousEnthalpyBackwardEuler::new(cx, &mesh,
            &[ReferenceEnthalpyMaterial { curve: &chart, reference_density_kg_m3: 1.0 }],
            &[0, 0], EnthalpyBudget { max_vertices: 8, max_elements: 2 }).unwrap();
        let material = ConductivityModel::isotropic_declared(2.0).unwrap();
        let config = ContactEnthalpyNewtonConfig {
            step: EnthalpyStepConfig {
                newton: NewtonKrylovConfig { relative_tolerance: 0.0, absolute_tolerance: 1e-11,
                    linear_restart: 8, max_linear_cycles: 16,
                    forcing_minimum: 1e-12, forcing_maximum: 1e-3, ..Default::default() },
                max_newton_iterations: 32, energy_tolerance_j: 1e-8,
            },
            max_residual_evaluations: 256, max_face_pairs: 1, residual_tolerance_j: 1e-8,
        };
        for (initial_temperature, volumetric_heat) in [
            (320.0, 10.0), (320.0, -10.0), (256.0, 10.0), (384.0, -10.0),
        ] {
            let source = ScalarField::Uniform(volumetric_heat);
            let problem = ConductionProblem { mesh: &mesh, boundary: &bc, material: &material,
                element_materials: None, source: &source };
            let old = [10.0 * (initial_temperature - 240.0); 8];
            let before = old;
            let next = contacts.advance_enthalpy_newton(cx, &engine, problem, &old, 1.0, config).unwrap();
            let expected_temperature = initial_temperature + volumetric_heat / 10.0;
            for (&h, &temperature) in next.conduction.specific_enthalpy_j_kg.iter()
                .zip(&next.conduction.temperature)
            {
                near(h, 10.0 * (expected_temperature - 240.0), 2e-8);
                near(temperature, expected_temperature, 2e-9);
            }
            // Two exact unit tetrahedra have total volume/mass 1/3. Contact
            // transfer is internal and uniform heating has no spatial flux.
            near(next.physical.stored_energy_change_j, volumetric_heat / 3.0, 1e-8);
            near(next.physical.energy_residual_j, 0.0, 1e-8);
            near(next.physical_contacts.states()[0].interface_temperature_k, expected_temperature, 2e-9);
            assert!(next.physical.residual_norm_j <= config.residual_tolerance_j);
            assert_eq!(old, before);
        }
        for (initial_temperature, volumetric_heat) in [(256.0, -10.0), (384.0, 10.0)] {
            let source = ScalarField::Uniform(volumetric_heat);
            let problem = ConductionProblem { mesh: &mesh, boundary: &bc, material: &material,
                element_materials: None, source: &source };
            let old = [10.0 * (initial_temperature - 240.0); 8];
            let before = old;
            assert!(contacts.advance_enthalpy_newton(cx, &engine, problem, &old, 1.0, config).is_err());
            assert_eq!(old, before);
        }
    });
}
