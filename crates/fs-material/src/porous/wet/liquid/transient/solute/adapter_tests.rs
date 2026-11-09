use super::super::{LiquidSaturationCell, LiquidSaturationStepInput};
use super::*;
use crate::porous::wet::liquid::{WetPoreLiquidLaw, WetPoreLiquidSupport};

#[test]
fn g1_sorption_consumes_the_original_hydraulic_interval_and_refuses_a_different_clock() {
    let law = WetPoreLiquidLaw::try_new(
        0.1,
        0.1,
        0.5,
        0.5,
        10_000.0,
        WetPoreLiquidSupport {
            liquid_saturation: (0.1, 0.9),
            temperature_k: (280.0, 330.0),
            liquid_pressure_pa: (50_000.0, 200_000.0),
            maximum_superficial_speed_m_s: 0.1,
        },
        "synthetic aqueous sorption fixture",
    )
    .unwrap();
    let original = LiquidSaturationColumn::new(
        vec![LiquidSaturationCell {
            thickness_m: 0.01,
            porosity: 0.4,
            intrinsic_permeability_m2: 2e-10,
            law,
        }],
        vec![0.5],
    )
    .unwrap();
    let p = super::super::pressures(original.states(), 150_000.0).unwrap()[0].unwrap();
    let water = original
        .step_candidate(
            LiquidSaturationStepInput {
                duration_s: 1.0,
                temperature_k: 300.0,
                gas_pressure_pa: 150_000.0,
                dynamic_viscosity_pa_s: 1e-3,
                body_force_pressure_gradient_pa_m: 0.0,
                boundary_liquid_pressure_pa: [p, p],
            },
            || false,
        )
        .unwrap();
    let controls = DiluteSoluteControls {
        absolute_equation_tolerance_mol_m2: [1e-14],
        relative_equation_tolerance: 1e-12,
        amount_balance_tolerance_mol_m2: [1e-14],
        water_balance_tolerance_m: 1e-14,
        maximum_cell_species: 1,
        maximum_storage_bytes: 64 * 1024,
    };
    let run = |column: &LiquidSaturationColumn| {
        water.transport_sorbing_dilute_solutes(
            column,
            &[[0.006]],
            &[[0.0]],
            &[[LinearSoluteSorption::Kinetic {
                capacity_m: 0.004,
                rate_s_inverse: 1.0,
            }]],
            &[[0.01]],
            &[1.0],
            [DiluteSoluteBoundary::Outflow; 2],
            controls,
            || false,
        )
    };
    let before = water.clone();
    let result = run(&original).unwrap();
    assert!((result.concentration_mol_m3[0][0] - 1.5).abs() < 1e-14);
    assert!((result.aqueous_amount_mol_m2[0][0] - 0.003).abs() < 1e-15);
    assert!((result.sorbed_amount_mol_m2[0][0] - 0.003).abs() < 1e-15);
    assert_eq!(result.report.total.interval_s, [0.0, 1.0]);
    assert!(result.report.total.amount_balance_defect_mol_m2[0].abs() < 1e-15);
    assert_eq!(water, before);
    let mut wrong = original.clone();
    wrong.elapsed_s = 1.0;
    assert!(matches!(run(&wrong), Err(DiluteSoluteError::Invalid(_))));
    assert_eq!(original.elapsed_s(), 0.0);
}
