//! One physical clock, refreshed coupled stages, and original closed balances.

use crate::miscibility::multicomponent::{
    advection::{MulticomponentMiscibilityAdvection, MulticomponentMiscibilityAdvectionFace},
    nonisothermal::{
        MulticomponentNonisothermalError, MulticomponentNonisothermalField,
        MulticomponentNonisothermalSolve,
        advection::MulticomponentNonisothermalAdvectionControls,
        tests::{Case, enthalpy, entropy, near},
    },
    tests::matrix,
};

fn opposing(
    field: &MulticomponentNonisothermalField,
    rate: f64,
) -> MulticomponentMiscibilityAdvection {
    field
        .network()
        .prepare_advection(
            vec![MulticomponentMiscibilityAdvectionFace {
                forward_m3_s: rate,
                reverse_m3_s: rate,
            }],
            0.0,
            || false,
        )
        .unwrap()
}

#[test]
fn same_interval_transport_refreshes_full_diffusion_after_stirring_and_closes_both_receipts() {
    let case = Case {
        a: matrix(3, &[(0, 1, 0.2), (0, 2, -0.1), (1, 2, 0.3)]),
        b: matrix(3, &[(0, 1, 120.0), (0, 2, -60.0), (1, 2, 80.0)]),
        conductivity: 1.2,
        heat: [500.0, -200.0, 0.0],
        ..Case::default()
    };
    let initial = case.field();
    let saved = initial.clone();
    let flow = opposing(&initial, 0.0003);
    let dt = 0.04;
    let advection_controls = MulticomponentNonisothermalAdvectionControls::default();
    let controls = MulticomponentNonisothermalSolve::default();
    let advected = initial
        .advect_candidate(&flow, dt, advection_controls, controls, || false)
        .unwrap();
    let refreshed = advected
        .field()
        .step_candidate(dt, controls, || false)
        .unwrap();
    let unrefreshed = initial.step_candidate(dt, controls, || false).unwrap();
    let mut polls = 0usize;
    let candidate = initial
        .transport_candidate(&flow, dt, advection_controls, controls, || {
            polls += 1;
            false
        })
        .unwrap();
    assert_eq!(candidate.field().elapsed_s(), dt);
    assert_eq!(candidate.step().advection.duration_s, dt);
    assert_eq!(candidate.step().diffusion.duration_s, dt);
    assert_eq!(refreshed.field().elapsed_s(), 2.0 * dt);
    assert_eq!(
        candidate.field().component_moles(),
        refreshed.field().component_moles()
    );
    assert_eq!(
        candidate.field().temperatures(),
        refreshed.field().temperatures()
    );
    let mut naive_error = 0.0_f64;
    for s in 0..3 {
        near(
            candidate.step().transferred_mol[0][s],
            advected.step().transferred_mol[0][s] + refreshed.step().transferred_mol[0][s],
            2e-15,
        );
        let naive = advected.field().component_moles()[s]
            + unrefreshed.field().component_moles()[s]
            - initial.component_moles()[s];
        naive_error = naive_error.max((candidate.field().component_moles()[s] - naive).abs());
        near(
            candidate.field().component_moles()[s] - initial.component_moles()[s],
            -candidate.step().transferred_mol[0][s],
            3e-15,
        );
        near(
            candidate.field().component_moles()[3 + s] - initial.component_moles()[3 + s],
            candidate.step().transferred_mol[0][s],
            3e-15,
        );
    }
    assert!(naive_error > 1e-8);
    near(
        candidate.step().transferred_energy_j[0],
        advected.step().transferred_energy_j[0] + refreshed.step().transferred_energy_j[0],
        2e-10,
    );
    let mut energy_change = 0.0;
    let mut entropy_change = 0.0;
    for cell in 0..2 {
        let old = initial.cell_component_moles(cell).unwrap();
        let new = candidate.field().cell_component_moles(cell).unwrap();
        let t = candidate.field().temperatures()[cell].value();
        let dh = enthalpy(new, t, &case) - enthalpy(old, case.temperatures[cell], &case);
        energy_change += dh;
        entropy_change += entropy(new, t, &case) - entropy(old, case.temperatures[cell], &case);
        near(
            dh,
            if cell == 0 {
                -candidate.step().transferred_energy_j[0]
            } else {
                candidate.step().transferred_energy_j[0]
            },
            1e-8,
        );
    }
    near(energy_change, 0.0, 1e-8);
    near(candidate.step().energy_residual_j, energy_change, 1e-8);
    near(
        candidate.step().entropy_change_j_per_k,
        entropy_change,
        1e-9,
    );
    assert!(entropy_change > 0.0);
    assert!(polls > 3);
    let mut seen = 0usize;
    assert!(matches!(
        initial.transport_candidate(&flow, dt, advection_controls, controls, || {
            seen += 1;
            seen == polls
        }),
        Err(MulticomponentNonisothermalError::Cancelled)
    ));
    assert_eq!(initial, saved);
    assert_eq!(
        candidate,
        initial
            .transport_candidate(&flow, dt, advection_controls, controls, || false)
            .unwrap()
    );
}

#[test]
fn zero_flow_and_zero_molecular_response_preserve_each_original_stage_and_clock() {
    let initial = Case {
        conductivity: 0.7,
        heat: [200.0, -100.0, 0.0],
        ..Case::default()
    }
    .field();
    let flow = opposing(&initial, 0.0);
    let dt = 0.05;
    let controls = MulticomponentNonisothermalSolve::default();
    let molecular = initial.step_candidate(dt, controls, || false).unwrap();
    let combined = initial
        .transport_candidate(&flow, dt, Default::default(), controls, || false)
        .unwrap();
    assert_eq!(combined.field(), molecular.field());
    assert_eq!(combined.step().diffusion, *molecular.step());
    assert_eq!(
        combined.step().transferred_mol,
        molecular.step().transferred_mol
    );
    assert_eq!(
        combined.step().transferred_energy_j,
        molecular.step().transferred_energy_j
    );

    // Equal compositions exclude molecular diffusion without pretending that
    // the strictly positive authored D matrix is zero. With k=Q=0, unequal
    // temperatures travel solely with the original parcels.
    let initial = Case {
        volumes: vec![0.001; 2],
        fractions: vec![[0.25, 0.25, 0.5]; 2],
        ..Case::default()
    }
    .field();
    let flow = opposing(&initial, 0.0005);
    let dt = 0.5;
    let advected = initial
        .advect_candidate(&flow, dt, Default::default(), controls, || false)
        .unwrap();
    let combined = initial
        .transport_candidate(&flow, dt, Default::default(), controls, || false)
        .unwrap();
    assert_eq!(combined.field(), advected.field());
    assert_eq!(
        combined.step().diffusion.transferred_mol,
        vec![vec![0.0; 3]]
    );
    assert_eq!(combined.step().diffusion.transferred_energy_j, vec![0.0]);
    assert_eq!(
        combined.step().transferred_mol,
        advected.step().transferred_mol
    );
    assert_eq!(
        combined.step().transferred_energy_j,
        advected.step().transferred_energy_j
    );
    near(combined.field().temperatures()[0].value(), 375.0, 1e-9);
    near(combined.field().temperatures()[1].value(), 525.0, 1e-9);
    assert_eq!(combined.field().elapsed_s(), dt);
}
