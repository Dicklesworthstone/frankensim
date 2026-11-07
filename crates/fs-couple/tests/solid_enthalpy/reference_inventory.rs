//! G0/G1: caloric and inertial owners must share the local reference inventory.
use super::*;
use fs_solid::loads3::transient::events::MaterialRebindConfig3;

#[test]
fn finite_heat_rejects_foreign_nodal_mass_but_accepts_equivalent_uniform_inventory() {
    with_cx(|cx| {
        let (card, caloric) = material_with_modulus_endpoint(60., 60., 1e-5, 10., 0., false);
        let thermal = initial(caloric.curve(), &[(30, 20_000.)], 1.);
        let context = point(300.);
        let model = SolidEnthalpy::try_new(
            &thermal, &NODES, &TETS, vec![source(30, &card, caloric.curve(), &context)],
            TetAssemblyBudget::standard(), 0.1, config().mesh_frame, cx,
        ).unwrap();
        let loads = model.geometry().assemble(&[], VolumeLoad3::None, &[], cx).unwrap();
        let start = model.resolve(&thermal, cx).unwrap();
        let heated = thermal.advance_to(0.01, thermal_config(), || Ok(())).unwrap().0;
        let next = model.resolve(&heated, cx).unwrap();
        let fixed: Vec<_> = (0..12).filter(|d| *d != 3).collect();
        let before = start.problem(&fixed);
        let after = next.problem(&fixed);
        let zero = [[0.; 3]; 4];
        let default_run = loads.begin_transient(
            &before, Some(start.thermal_strains()), 0., &zero, &zero, Default::default(), cx,
        ).unwrap();
        let (expected, _) = next.advance_transient_candidate(
            &start, &default_run, &after, &loads, cx,
        ).unwrap();
        for (density, mass) in [([12.; 4], 2.), ([3., 5., 7., 9.], 1.), ([6.; 4], 1.)] {
            let run = loads.begin_transient_with_nodal_reference_density(
                &before, &density, Some(start.thermal_strains()), 0., &zero, &zero,
                Default::default(), cx,
            ).unwrap();
            near(run.reference_mass_kg(), mass, 1e-12);
            let original = run.state().clone();
            let result = next.advance_transient_candidate(&start, &run, &after, &loads, cx);
            if density == [6.; 4] {
                let (accepted, report) = result.unwrap();
                near(accepted.reference_mass_kg(), 1., 1e-12);
                near(accepted.state().displacement[0], expected.state().displacement[0], 1e-14);
                near(accepted.state().velocity[0], expected.state().velocity[0], 1e-12);
                near(report.physical_balance_error_j, 0., 1e-12);
            } else {
                // The graded field has the same total mass but different local
                // inertia; a global mass-only comparison would accept it.
                assert!(matches!(result, Err(SolidEnthalpyError::Invalid(_))));
            }
            assert_eq!(run.state(), &original);
        }
        assert_eq!(thermal.cells()[&CellId(30)].specific_h_j_kg, 0.);
        near(heated.cells()[&CellId(30)].specific_h_j_kg, 200., 1e-9);
    });
}

#[test]
fn material_replacement_cannot_retain_an_override_from_a_different_caloric_inventory() {
    with_cx(|cx| {
        let (card, caloric) = material_with_density(60., 60., 1e-5, 10., 0., false, 6.);
        let (replacement_card, replacement_caloric) =
            material_with_density(60., 60., 1e-5, 10., 0., false, 12.);
        let thermal = initial(caloric.curve(), &[(30, 0.)], 1.);
        let replacement_heat = initial(replacement_caloric.curve(), &[(30, 0.)], 2.);
        let context = point(300.);
        let model = SolidEnthalpy::try_new(
            &thermal, &NODES, &TETS, vec![source(30, &card, caloric.curve(), &context)],
            TetAssemblyBudget::standard(), 0.1, config().mesh_frame, cx,
        ).unwrap();
        let replacement = SolidEnthalpy::try_new(
            &replacement_heat, &NODES, &TETS,
            vec![source(30, &replacement_card, replacement_caloric.curve(), &context)],
            TetAssemblyBudget::standard(), 0.1, config().mesh_frame, cx,
        ).unwrap();
        let loads = model.geometry().assemble(&[], VolumeLoad3::None, &[], cx).unwrap();
        let replacement_loads = replacement.geometry().assemble(&[], VolumeLoad3::None, &[], cx).unwrap();
        let start = model.resolve(&thermal, cx).unwrap();
        let next = replacement.resolve(&replacement_heat, cx).unwrap();
        let fixed: Vec<_> = (0..12).filter(|d| *d != 3).collect();
        let before = start.problem(&fixed);
        let after = next.problem(&fixed);
        let zero = [[0.; 3]; 4];
        let run = loads.begin_transient_with_nodal_reference_density(
            &before, &[6.; 4], Some(start.thermal_strains()), 0., &zero, &zero,
            Default::default(), cx,
        ).unwrap();
        let controls = MaterialRebindConfig3 {
            maximum_absolute_work_j: 1.,
            absolute_momentum_tolerance_kg_m_s: 1e-12,
            maximum_strain: 0.1,
        };
        let original = run.state().clone();
        assert!(matches!(
            next.rebind_transient_candidate(&start, &run, &after, &replacement_loads, controls, cx),
            Err(SolidEnthalpyError::Invalid(_))
        ));
        assert_eq!(run.state(), &original);
        let element_density_run = loads.begin_transient(
            &before, Some(start.thermal_strains()), 0., &zero, &zero, Default::default(), cx,
        ).unwrap();
        let (accepted, report) = next.rebind_transient_candidate(
            &start, &element_density_run, &after, &replacement_loads, controls, cx,
        ).unwrap();
        near(accepted.reference_mass_kg(), 2., 1e-12);
        near(report.mechanical.reference_mass_kg[0], 1., 1e-12);
        near(report.mechanical.reference_mass_kg[1], 2., 1e-12);
    });
}
