//! G0/G1: native spatial heat and full solid motion share local reference mass.
use super::*;

#[test]
fn spatial_heat_and_motion_admit_only_the_same_elementwise_reference_inventory() {
    with_cx(|cx| {
        let (card, caloric) = fixture(2., 60., 1e-5, 10., 500., 45.);
        let mesh = mesh();
        let native = native(&mesh, caloric.curve(), 6., cx);
        let context = point(300.);
        let model = bind(&native, source(&card, caloric.curve(), &context), cx).unwrap();
        let geometry = model.geometry(cx).unwrap();
        let loads = geometry.assemble(&[], VolumeLoad3::None, &[], cx).unwrap();
        let boundary = boundary(&mesh);
        let source = ScalarField::uniform("volumetric heating", 120_000.).unwrap();
        let old = [0.; 4];
        let start = model.resolve_at(&old, 0., cx).unwrap();
        let fixed: Vec<_> = (0..12).filter(|d| *d != 3).collect();
        let problem = start.problem(&fixed);
        let zero = [[0.; 3]; 4];
        let default_run = loads.begin_transient(
            &problem, Some(start.thermal_strains()), 0., &zero, &zero,
            Default::default(), cx,
        ).unwrap();
        let expected = model.advance_transient(
            &old, 0.01, &boundary, None, &source, thermal_config(), EnthalpyRadiation::None,
            &default_run, &loads, cx,
        ).unwrap();
        near(expected.thermal.conduction().stored_energy_change_j, 200., 1e-8);
        let next_problem = expected.state.problem(&fixed);
        for (density, mass) in [([12.; 4], 2.), ([3., 5., 7., 9.], 1.), ([6.; 4], 1.)] {
            let run = loads.begin_transient_with_nodal_reference_density(
                &problem, &density, Some(start.thermal_strains()), 0., &zero, &zero,
                Default::default(), cx,
            ).unwrap();
            near(run.reference_mass_kg(), mass, 1e-12);
            let original = run.state().clone();
            let direct = expected.state.advance_transient_candidate(
                &start, &run, &next_problem, &loads, cx,
            );
            let integrated = model.advance_transient(
                &old, 0.01, &boundary, None, &source, thermal_config(), EnthalpyRadiation::None,
                &run, &loads, cx,
            );
            if density == [6.; 4] {
                let (direct, _) = direct.unwrap();
                let integrated = integrated.unwrap();
                near(direct.reference_mass_kg(), 1., 1e-12);
                for actual in [direct.state(), integrated.motion.state()] {
                    near(actual.displacement[0], expected.motion.state().displacement[0], 1e-14);
                    near(actual.velocity[0], expected.motion.state().velocity[0], 1e-12);
                    assert_eq!(actual.time_s, 0.01);
                }
                assert_eq!(
                    integrated.thermal.conduction().specific_enthalpy_j_kg,
                    expected.thermal.conduction().specific_enthalpy_j_kg,
                );
            } else {
                // [3,5,7,9] has mean density 6 and total mass 1 kg, but
                // distinct nodal inertia cannot share the scalar caloric cells.
                assert!(matches!(direct, Err(SpatialSolidEnthalpyError::Invalid(_))));
                assert!(matches!(integrated, Err(SpatialSolidEnthalpyError::Invalid(_))));
            }
            assert_eq!(run.state(), &original);
        }
        assert_eq!(native.reference_nodal_masses_kg(), &[0.25; 4]);
        assert_eq!(old.map(f64::to_bits), [0.; 4].map(f64::to_bits));
    });
}
