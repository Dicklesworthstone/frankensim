//! Original NASA sources on composition-selective material streams.

use super::*;
use fs_thermal::enthalpy::advection::{
    Advected, MassBoundary, SelectiveOutlet, SelectiveOutletLaw, species::InletSpecies,
};

fn law(source: &ThermochemicalGasModel, fractions: Vec<f64>) -> GasEnthalpyLaw {
    GasEnthalpyLaw::new(
        source,
        Composition::new(CompositionBasis::MassFraction, fractions).unwrap(),
        point(320.0, 100000.0),
        controls(),
    )
    .unwrap()
}

fn selective_flow(
    inlet: &GasEnthalpyLaw,
    selected: &GasEnthalpyLaw,
    selected_t: f64,
) -> Advected<()> {
    let composition = selected.advected_composition().unwrap();
    Advected::new((), vec![], 0)
        .unwrap()
        .with_boundaries(
            vec![
                MassBoundary::Inlet {
                    id: 1,
                    cell: CellId(1),
                    kg_s: 1.0,
                    material_id: MaterialId(2),
                    specific_h_j_kg: inlet.enthalpy_at_temperature(320.0).unwrap(),
                },
                MassBoundary::Outlet {
                    id: 2,
                    cell: CellId(1),
                    kg_s: 0.75,
                },
                MassBoundary::Outlet {
                    id: 3,
                    cell: CellId(1),
                    kg_s: 0.25,
                },
            ],
            3,
        )
        .unwrap()
        .with_selective_outlets(
            vec![SelectiveOutlet {
                boundary_id: 3,
                mass_fractions: composition.mass_fractions.to_vec(),
                specific_h_j_kg: selected.enthalpy_at_temperature(selected_t).unwrap(),
                enthalpy_reference: composition.enthalpy_reference.to_owned(),
            }],
            2,
        )
        .unwrap()
        .with_stationary_species(
            vec!["A".into(), "B".into()],
            vec![InletSpecies {
                boundary_id: 1,
                mass_fractions: vec![0.5, 0.5],
            }],
            1,
            2,
        )
        .unwrap()
}

// Independent ideal-mixture expression for the synthetic NASA source:
// Cp_A/R=3.5+0.001*T, Cp_B/R=2.5, M=[0.028,0.004], p=p0.
fn reduced_entropy(fractions: [f64; 2], t: f64) -> f64 {
    let amounts = [fractions[0] / 0.028, fractions[1] / 0.004];
    let total = amounts.iter().sum::<f64>();
    let thermal = [
        3.5 * (t / 280.0).ln() + 0.001 * (t - 280.0),
        2.5 * (t / 280.0).ln(),
    ];
    UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K
        * amounts
            .into_iter()
            .zip(thermal)
            .filter(|(n, _)| *n > 0.0)
            .map(|(n, s)| n * (s - (n / total).ln()))
            .sum::<f64>()
}

fn specific_h(fractions: [f64; 2], t: f64) -> f64 {
    UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K
        * (fractions[0] / 0.028 * (3.5 * t + 0.0005 * t * t) + fractions[1] / 0.004 * 2.5 * t)
}

#[test]
fn g1_source_selective_stream_retains_separation_entropy_under_large_reference_shifts() {
    let boundary_rate = reduced_entropy([0.5, 0.5], 320.0)
        - 0.75 * reduced_entropy([1.0 / 3.0, 2.0 / 3.0], 320.0)
        - 0.25 * reduced_entropy([1.0, 0.0], 320.0);
    assert!(boundary_rate > 0.0);
    for offset in [0.0, 2.0_f64.powi(50), -2.0_f64.powi(50)] {
        let source = entropy_shifted_source(offset);
        let remaining = law(&source, vec![1.0 / 3.0, 2.0 / 3.0]);
        let inlet = law(&source, vec![0.5, 0.5]);
        let selected = law(&source, vec![1.0, 0.0]);
        let h0 = remaining.enthalpy_at_temperature(320.0).unwrap();
        let initial = thermal_timeline(
            &[
                (MaterialId(1), &remaining),
                (MaterialId(2), &inlet),
                (MaterialId(3), &selected),
            ],
            h0,
            0.0,
        );
        let unbound = selective_flow(&inlet, &selected, 320.0);
        assert!(
            unbound
                .endpoint_transport_with_entropy(&initial, || Ok(()))
                .is_err()
        );
        let flow = unbound
            .clone()
            .with_selective_outlet_laws(vec![SelectiveOutletLaw {
                boundary_id: 3,
                material_id: MaterialId(3),
            }])
            .unwrap();
        let (_, rates) = flow
            .endpoint_transport_with_entropy(&initial, || Ok(()))
            .unwrap();
        assert!((rates.boundary_entropy_input_w_k - boundary_rate).abs() < 2e-9);
        // Isothermal separation has no enthalpy change but reduces mixing
        // entropy. Keep that deficit visible: this is no work-free separator.
        assert!((rates.upwind_mixing_w_k + boundary_rate).abs() < 2e-9);
        assert!(rates.upwind_mixing_w_k < 0.0);
        assert!(rates.balance_residual_w_k.abs() < 1e-9);
        let (next, _, energy, entropy) = initial
            .advance_advected_with_entropy(&flow, 0.05, thermal_config(), || Ok(()))
            .unwrap();
        let carried = entropy.advection.unwrap();
        assert!((carried.boundary_entropy_input_j_k - 0.05 * boundary_rate).abs() < 2e-9);
        assert!((carried.upwind_mixing_j_k + 0.05 * boundary_rate).abs() < 2e-9);
        assert!(entropy.thermal_storage_change_j_k.abs() < 1e-9);
        assert!(energy.external_heat_j.abs() < 1e-8);
        assert!((next.cells()[&CellId(1)].specific_h_j_kg - h0).abs() < 1e-8);
        let wrong = unbound
            .with_selective_outlet_laws(vec![SelectiveOutletLaw {
                boundary_id: 3,
                material_id: MaterialId(2),
            }])
            .unwrap();
        assert!(
            wrong
                .endpoint_transport_with_entropy(&initial, || Ok(()))
                .is_err()
        );
        assert!(matches!(
            initial.advance_advected_with_entropy(&flow, 0.05, thermal_config(), || Err(
                DepositionError::Cancelled
            )),
            Err(DepositionError::Cancelled)
        ));
        assert_eq!(initial.cells()[&CellId(1)].specific_h_j_kg, h0);
    }
}

#[test]
fn g1_source_selective_stream_uses_its_own_caloric_state_and_conserves_heat_once() {
    let source = model(vec![0.5, 0.5], 0.001);
    let remaining_y = [1.0 / 3.0, 2.0 / 3.0];
    let remaining = law(&source, remaining_y.to_vec());
    let inlet = law(&source, vec![0.5, 0.5]);
    let selected = law(&source, vec![1.0, 0.0]);
    let h0 = remaining.enthalpy_at_temperature(320.0).unwrap();
    let initial = thermal_timeline(
        &[
            (MaterialId(1), &remaining),
            (MaterialId(2), &inlet),
            (MaterialId(3), &selected),
        ],
        h0,
        0.0,
    );
    for selected_t in [300.0, 340.0] {
        let flow = selective_flow(&inlet, &selected, selected_t)
            .with_selective_outlet_laws(vec![SelectiveOutletLaw {
                boundary_id: 3,
                material_id: MaterialId(3),
            }])
            .unwrap();
        let power = specific_h([0.5, 0.5], 320.0)
            - 0.75 * specific_h(remaining_y, 320.0)
            - 0.25 * specific_h([1.0, 0.0], selected_t);
        let boundary_rate = reduced_entropy([0.5, 0.5], 320.0)
            - 0.75 * reduced_entropy(remaining_y, 320.0)
            - 0.25 * reduced_entropy([1.0, 0.0], selected_t);
        let (transport, rates) = flow
            .endpoint_transport_with_entropy(&initial, || Ok(()))
            .unwrap();
        assert!((transport.heat.bulk.sources_w[&CellId(1)] - power).abs() < 1e-8);
        assert!((rates.boundary_entropy_input_w_k - boundary_rate).abs() < 2e-9);
        assert!((rates.upwind_mixing_w_k - (power / 320.0 - boundary_rate)).abs() < 2e-9);
        let dt = 0.01;
        let (next, _, energy, entropy) = initial
            .advance_advected_with_entropy(&flow, dt, thermal_config(), || Ok(()))
            .unwrap();
        assert!((next.cells()[&CellId(1)].specific_h_j_kg - h0 - dt * power).abs() < 1e-8);
        assert!((energy.external_heat_j - dt * power).abs() < 1e-8);
        let t = next.temperature_k(CellId(1)).unwrap();
        let storage = reduced_entropy(remaining_y, t) - reduced_entropy(remaining_y, 320.0);
        assert!((entropy.thermal_storage_change_j_k - storage).abs() < 2e-9);
        assert!((entropy.temporal_defect_j_k - (storage - dt * power / 320.0)).abs() < 2e-9);
    }
}
