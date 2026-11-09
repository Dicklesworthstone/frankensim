use super::*;

#[test]
fn g1_complete_original_nodal_tones_keep_element_heat_and_refine_exactly() {
    with_cx(|cx, gate| {
        let mut nodes = CORNERS.to_vec();
        nodes.extend(CORNERS.map(|mut p| {
            p[0] += 2.0;
            p
        }));
        let cells = [[0, 1, 2, 3], [4, 5, 6, 7]];
        // Every remaining original free coordinate survives; no modal basis.
        let fixed: Vec<_> = (0..24).filter(|&i| i != 11 && i != 23).collect();
        let materials = [material(0.0002, 0.3), material(0.0001, 0.5)];
        let source = ViscoElements::new(
            SolidMeshRef::Tet4 {
                nodes_m: &nodes,
                elements: &cells,
            },
            TetViscoField::PerElement(&materials),
            &fixed,
            TetAssemblyBudget::standard(),
            ViscoStateBudget {
                max_internal_variables: 1000,
            },
            cx,
        )
        .unwrap();
        let frequencies = [500.0, 2200.0];
        let make = |h| {
            PreparedTensorSolidBody::new_tones(
                source
                    .project_basis(&identity(2), projection_budget(), cx)
                    .unwrap(),
                h,
                &frequencies,
                1.0,
                controls(),
                gate,
            )
            .unwrap()
        };
        let whole = make(0.0006);
        let half = make(0.0003);
        let mut a = whole
            .history(
                &[1e-6, -0.5e-6],
                &[0.01, -0.005],
                InitialMemory::Instantaneous,
                gate,
            )
            .unwrap();
        let mut b = half
            .history(
                &[1e-6, -0.5e-6],
                &[0.01, -0.005],
                InitialMemory::Instantaneous,
                gate,
            )
            .unwrap();
        let mut cosine = vec![vec![[0.0; 3]; 8]; 2];
        let mut sine = cosine.clone();
        cosine[0][3][2] = 0.1;
        cosine[0][7][2] = -0.07;
        cosine[1][3][2] = -0.05;
        sine[0][7][2] = 0.03;
        sine[1][3][2] = 0.09;
        sine[1][7][2] = 0.06;
        // Original support forces remain present and have exactly zero motion.
        cosine[0][0] = [1000.0; 3];
        sine[1][4] = [-900.0; 3];
        let full = a.advance_tones(&cosine, &sine, gate).unwrap();
        let left = b.advance_tones(&cosine, &sine, gate).unwrap();
        let (next_cosine, next_sine) = half.advance_force_phase(&cosine, &sine, gate).unwrap();
        let right = b.advance_tones(&next_cosine, &next_sine, gate).unwrap();
        for (x, y) in a
            .material_history()
            .state()
            .iter()
            .zip(b.material_history().state())
        {
            assert!((x - y).abs() < 3e-11, "full={x:e}, half={y:e}");
        }
        for group in 0..2 {
            let expected =
                left.motion.element_heat_j()[group] + right.motion.element_heat_j()[group];
            assert!(full.motion.element_heat_j()[group] > 0.0);
            assert!((full.motion.element_heat_j()[group] - expected).abs() < 1e-10);
        }
        assert!(
            (full.motion.material.energy.supplied_work_j
                - left.motion.material.energy.supplied_work_j
                - right.motion.material.energy.supplied_work_j)
                .abs()
                < 1e-10
        );
        assert!(full.motion.nodal_work_error_j.abs() < 1e-12);
        for moments in full
            .cosine_nodal_velocity_moment_mps
            .iter()
            .chain(&full.sine_nodal_velocity_moment_mps)
        {
            for &d in &fixed {
                assert_eq!(moments[d / 3][d % 3], 0.0);
            }
        }
        let before = a.material_history().state().to_vec();
        let heat = a.material_history().dissipated_heat_j();
        let mut bad = cosine.clone();
        bad[0].pop();
        assert!(a.advance_tones(&bad, &sine, gate).is_err());
        assert_eq!(a.material_history().state(), before);
        assert_eq!(a.material_history().dissipated_heat_j(), heat);
    });
}
