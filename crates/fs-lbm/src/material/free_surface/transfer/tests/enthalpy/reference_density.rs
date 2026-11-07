//! The original reference inventory owns both rho*h and cubic mechanical work.

use super::super::reference_density::{compatible_reference, density_controls};
use super::*;

#[test]
fn deformed_original_density_retains_spatial_enthalpy_through_native_motion() {
    let card = thermal_card(1000.0, 2e-5);
    let law = LiquidTemperatureLaw::try_new(
        &card,
        &point(),
        "T",
        LiquidViscosityBasis::Dynamic,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    let surface = resolve_liquid_surface_state(
        &card,
        &point(),
        LiquidViscosityBasis::Dynamic,
        LiquidCapillarity::Neglected,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    let caloric = chart(&card, 1000.0);
    let (p, t, v) = cube([0.0; 3]);
    let (reference, density) = compatible_reference(&p, &t, 1000.0);
    let h: Vec<_> = p.iter().copied().map(affine_h).collect();
    let unchanged = (reference.clone(), density.clone(), h.clone());
    let gate = CancelGate::new_clock_free();
    with_cx(&gate, |cx| {
        let scaling = plan_free_surface_from_material(&surface, 0.001, 0.001, 0.1).unwrap();
        let (mut thermal, report) = scaling
            .free_surface3_from_reference_tetrahedra_with_enthalpy(
                domain(),
                [0.0; 3],
                source(&p, &t, &v, &h),
                TetrahedralLiquidReferenceDensity {
                    positions_m: &reference,
                    corner_density_kg_m3: &density,
                },
                density_controls(),
                &law,
                &caloric,
                thermal_limits(),
                FreeSurfaceEnthalpyControls::default(),
                cx,
                &gate,
            )
            .unwrap();
        let mass = 1000.0 * 0.004_f64.powi(3);
        let energy = mass * affine_h([0.004; 3]);
        close(report.source_enthalpy_j, energy, 2e-12);
        close(report.target_enthalpy_j, energy, 2e-12);
        close(report.enthalpy_defect_j, 0.0, 2e-12);
        let inventory = report.mechanical.reference_density.unwrap();
        close(inventory.reference.mass, mass, 2e-18);
        close(
            report.mechanical.target.mass_kg,
            inventory.reference.mass,
            2e-18,
        );
        let actual_volume: f64 = report
            .mechanical
            .components
            .iter()
            .flat_map(|piece| &piece.target_cells)
            .map(|&[x, y, z]| {
                thermal
                    .flow()
                    .sample_at_si(x, y, z)
                    .unwrap()
                    .unwrap()
                    .fill_fraction
                    * 0.001_f64.powi(3)
            })
            .sum();
        close(inventory.native_occupied_volume_m3, actual_volume, 2e-22);
        assert!(report.mechanical.kinetic_energy_change_j > 0.0);
        for z in 2..6 {
            for y in 2..6 {
                for x in 2..6 {
                    close(
                        thermal.specific_enthalpy_at([x, y, z]).unwrap(),
                        affine_h([x, y, z].map(|i| (i as f64 + 0.5) * 0.001)),
                        3e-8,
                    );
                }
            }
        }
        let step = thermal.step(&gate).unwrap();
        close(step.mass_after_kg, mass, 1e-14);
        close(step.enthalpy_before_j, energy, 2e-12);
        close(step.enthalpy_after_j, energy, 1e-9);
        close(step.external_heat_j, 0.0, 0.0);
        close(step.viscous_heat_j, 0.0, 0.0);
        close(step.balance_defect_j, 0.0, 1e-9);

        let incompatible: Vec<_> = density
            .iter()
            .map(|rho| rho.map(|value| value * 1.02))
            .collect();
        assert!(
            scaling
                .free_surface3_from_reference_tetrahedra_with_enthalpy(
                    domain(),
                    [0.0; 3],
                    source(&p, &t, &v, &h),
                    TetrahedralLiquidReferenceDensity {
                        positions_m: &reference,
                        corner_density_kg_m3: &incompatible,
                    },
                    density_controls(),
                    &law,
                    &caloric,
                    thermal_limits(),
                    FreeSurfaceEnthalpyControls::default(),
                    cx,
                    &gate,
                )
                .is_err()
        );
    });
    assert_eq!((reference, density, h), unchanged);
}
