//! Several measured faces exchange through one actual native liquid inventory.

use super::*;
use fs_couple::liquid_gas::evaporation::multiple::{
    LiquidGasHeatSurface, advance_liquid_surfaces_with_gas_and_heat,
};
use fs_material::evaporation::gas::GasEvaporationError;
use fs_material::evaporation::gas::chamber::heat::GasLiquidChamberHeatSolve;
use fs_thermal::evaporation::{InterfaceKinetics, KineticEvaporation};

const CORNER: [usize; 3] = [3, 6, 6];

fn corner_faces<'a>(
    law: &'a dyn EvaporationLaw,
    binding: &'a GasLiquidBinding,
) -> [LiquidGasEvaporationSurface<'a>; 3] {
    [
        ([4, 6, 6], 1e-6, 0.006),
        ([3, 7, 6], 0.6e-6, 0.02),
        ([3, 6, 7], 0.9e-6, 0.065),
    ]
    .map(
        |(gas, area_m2, mass_transfer_coefficient_m_s)| LiquidGasEvaporationSurface {
            liquid: CORNER,
            gas,
            area_m2,
            mass_transfer_coefficient_m_s,
            law,
            binding,
        },
    )
}

fn face_controls() -> DistributedLiquidGasSolve {
    let mut controls = joint_controls();
    controls.chamber.maximum_surfaces = 3;
    controls
}

fn requested_masses(
    faces: &[LiquidGasEvaporationSurface<'_>],
    liquid_t: f64,
    gas_t: f64,
    total_pressure: f64,
    vapor_pressure: f64,
) -> Vec<f64> {
    faces
        .iter()
        .map(|face| {
            DT * face.area_m2
                * surface_flux(
                    face.law,
                    liquid_t,
                    GasFilm {
                        total_pressure_pa: total_pressure,
                        vapor_partial_pressure_pa: vapor_pressure,
                        temperature_k: gas_t,
                        mass_transfer_coefficient_m_s: face.mass_transfer_coefficient_m_s,
                    },
                )
                .unwrap()
                .mass_flux_kg_m2_s
        })
        .collect()
}

#[test]
fn g1_three_native_faces_sum_heterogeneous_kinetic_rates_on_one_inventory_and_clock() {
    with_liquid_interfaces(0.0, &[(CORNER, 0.4, 300.0)], |liquid, gate| {
        let source = source(-13.4e6);
        let binding = binding(&source, 0.0);
        let law = KineticEvaporation::try_new(
            &Phase,
            InterfaceKinetics::Schrage {
                accommodation_coefficient: 1e-4,
            },
        )
        .unwrap();
        let faces = corner_faces(&law, &binding);
        let area: f64 = faces.iter().map(|face| face.area_m2).sum();
        assert!(area > liquid.flow().scaling().cell_size_m().powi(2));
        for vapor_pressure in [2000.0, 6000.0] {
            let gas = gas(&source, 300.0, vapor_pressure);
            let expected: f64 = requested_masses(&faces, 300.0, 300.0, 100000.0, vapor_pressure)
                .into_iter()
                .sum();
            let next = advance_liquid_surfaces_with_gas(
                liquid,
                &gas,
                &source,
                &faces,
                face_controls(),
                gate,
            )
            .unwrap();
            let report = next.report();
            assert_eq!(report.surfaces.len(), 1);
            assert_eq!(report.mass_exchange.samples.len(), 1);
            let row = &report.surfaces[0];
            assert_eq!(row.coordinates, CORNER);
            assert_eq!(row.phase.id, 0);
            assert_eq!(row.phase.transferred_mass_kg > 0.0, vapor_pressure < 3536.8);
            near(row.phase.transferred_mass_kg, expected, 1e-18);
            near(row.phase.liquid.area_m2, area, 1e-20);
            assert_eq!(row.native.mass_change_kg, -row.phase.transferred_mass_kg);
            assert_eq!(next.liquid().accepted_steps(), 1);
            assert_eq!(next.liquid().physical_time_s(), DT);
            complete_inventory(liquid, &gas, &binding, &next);

            // Averaging k before the kinetic resistance is a different model.
            let mut averaged = faces[0];
            averaged.area_m2 = area;
            averaged.mass_transfer_coefficient_m_s = faces
                .iter()
                .map(|face| face.area_m2 * face.mass_transfer_coefficient_m_s)
                .sum::<f64>()
                / area;
            let wrong = requested_masses(&[averaged], 300.0, 300.0, 100000.0, vapor_pressure)[0];
            assert!((wrong - expected).abs() > 0.01 * expected.abs());

            let reversed = advance_liquid_surfaces_with_gas(
                liquid,
                &gas,
                &source,
                &[faces[2], faces[0], faces[1]],
                face_controls(),
                gate,
            )
            .unwrap();
            assert_eq!(next.gas(), reversed.gas());
            assert_eq!(report.exchange, reversed.report().exchange);
            assert_eq!(report.mass_exchange, reversed.report().mass_exchange);
            assert_eq!(
                next.liquid().flow().lattice().grid().state_bits(),
                reversed.liquid().flow().lattice().grid().state_bits(),
            );
        }
        assert_eq!(liquid.accepted_steps(), 0);
        assert_eq!(liquid.physical_time_s(), 0.0);
    });
}

#[test]
fn g1_two_native_faces_use_actual_heated_temperature_and_one_reciprocal_heat_account() {
    with_liquid_interfaces(0.0, &[(CORNER, 0.4, 295.0)], |liquid, gate| {
        let source = source(-13.4e6);
        let gas = gas(&source, 300.0, 3500.0);
        let binding = binding(&source, 0.0);
        let law = KineticEvaporation::try_new(
            &Phase,
            InterfaceKinetics::HertzKnudsen {
                accommodation_coefficient: 1e-4,
            },
        )
        .unwrap();
        let faces = corner_faces(&law, &binding);
        let heated = [
            LiquidGasHeatSurface {
                surface: faces[0],
                heat_transfer_coefficient_w_m2_k: 500.0,
                heat_controls: heat_controls(),
            },
            LiquidGasHeatSurface {
                surface: faces[1],
                heat_transfer_coefficient_w_m2_k: 1000.0,
                heat_controls: heat_controls(),
            },
        ];
        let next = advance_liquid_surfaces_with_gas_and_heat(
            liquid,
            &gas,
            &source,
            &heated,
            face_controls(),
            GasLiquidChamberHeatSolve {
                maximum_iterations: 128,
            },
            gate,
        )
        .unwrap();
        let report = next.report();
        assert_eq!(report.surfaces.len(), 1);
        assert_eq!(report.mass_exchange.samples.len(), 1);
        let row = &report.surfaces[0];
        let heat = row.heat.as_ref().unwrap();
        let gas_heat = report.heat.as_ref().unwrap();
        let liquid_t = row.phase.surface_temperature_k[0];
        let gas_t = report.phase.gas_temperature_k[0];
        assert!(liquid_t > 295.0);
        assert!(gas_t < 300.0);
        assert_eq!(heat.temperature_k[1], liquid_t);
        assert_eq!(gas_heat.gas_temperature_k[1], gas_t);
        let conductance: f64 = heated
            .iter()
            .map(|face| face.surface.area_m2 * face.heat_transfer_coefficient_w_m2_k)
            .sum();
        near(
            heat.heat_into_liquid_j,
            DT * conductance * (gas_t - liquid_t),
            1e-12,
        );
        near(heat.heat_into_liquid_j, gas_heat.heat_into_liquid_j, 1e-14);
        let expected: f64 = requested_masses(
            &faces[..2],
            liquid_t,
            gas_t,
            100000.0 * gas_t / 300.0,
            row.phase.vapor_partial_pressure_pa[0],
        )
        .into_iter()
        .sum();
        near(row.phase.transferred_mass_kg, expected, 1e-18);
        assert_eq!(next.liquid().accepted_steps(), 1);
        complete_inventory(liquid, &gas, &binding, &next);

        let reversed = advance_liquid_surfaces_with_gas_and_heat(
            liquid,
            &gas,
            &source,
            &[heated[1], heated[0]],
            face_controls(),
            GasLiquidChamberHeatSolve {
                maximum_iterations: 128,
            },
            gate,
        )
        .unwrap();
        assert_eq!(next.gas(), reversed.gas());
        assert_eq!(report.exchange, reversed.report().exchange);
        assert_eq!(report.mass_exchange, reversed.report().mass_exchange);
    });
}

#[test]
fn g4_native_faces_share_the_original_liquid_donor_allowance() {
    with_liquid_interfaces(0.0, &[(CORNER, 0.4, 300.0)], |liquid, gate| {
        let source = source(-13.4e6);
        let gas = gas(&source, 300.0, 2000.0);
        let binding = binding(&source, 0.0);
        let faces = corner_faces(&Phase, &binding);
        let requested = requested_masses(&faces, 300.0, 300.0, 100000.0, 2000.0);
        let largest = requested.iter().copied().fold(0.0_f64, f64::max);
        let initial = liquid.project_interface_mass_exchange(CORNER, 0.0).unwrap();
        let allowed = 1.05 * largest;
        assert!(requested.iter().all(|&mass| mass < allowed));
        assert!(requested.iter().sum::<f64>() > allowed);
        let mut controls = face_controls();
        controls.chamber.surface.maximum_inventory_fraction = allowed / initial.initial_mass_kg;
        let bits = liquid.flow().lattice().grid().state_bits();
        let original_gas = gas.clone();
        assert!(matches!(
            advance_liquid_surfaces_with_gas(liquid, &gas, &source, &faces, controls, gate),
            Err(EvaporatingLiquidGasError::Phase(
                GasEvaporationError::StepSizeLimit { .. }
            ))
        ));
        assert_eq!(gas, original_gas);
        assert_eq!(liquid.flow().lattice().grid().state_bits(), bits);
        assert_eq!(liquid.accepted_steps(), 0);
    });
}

#[test]
fn g0_multiface_geometry_source_binding_and_heat_control_mismatches_refuse_atomically() {
    with_liquid_interfaces(0.0, &[(CORNER, 0.4, 300.0)], |liquid, gate| {
        let source = source(-13.4e6);
        let gas = gas(&source, 300.0, 2000.0);
        let binding = binding(&source, 0.0);
        let first = KineticEvaporation::try_new(&Phase, InterfaceKinetics::Equilibrium).unwrap();
        let other = KineticEvaporation::try_new(&Phase, InterfaceKinetics::Equilibrium).unwrap();
        let faces = corner_faces(&first, &binding);
        let mut foreign_binding = binding.clone();
        foreign_binding.liquid_enthalpy_offset_j_kg += 1.0;
        let bits = liquid.flow().lattice().grid().state_bits();
        let original_h = liquid.total_enthalpy_j().unwrap();
        let original_gas = gas.clone();
        for change in 0..6 {
            let mut rows = [faces[0], faces[1]];
            match change {
                0 => rows[1] = rows[0],
                1 => rows[1].gas = [2, 6, 6],
                2 => rows[1].gas = [4, 7, 6],
                3 => rows[1].area_m2 = 2e-6,
                4 => rows[1].binding = &foreign_binding,
                _ => rows[1].law = &other,
            }
            assert!(
                advance_liquid_surfaces_with_gas(
                    liquid,
                    &gas,
                    &source,
                    &rows,
                    face_controls(),
                    gate,
                )
                .is_err()
            );
        }
        let mut heated = faces[..2]
            .iter()
            .map(|&surface| LiquidGasHeatSurface {
                surface,
                heat_transfer_coefficient_w_m2_k: 10.0,
                heat_controls: heat_controls(),
            })
            .collect::<Vec<_>>();
        heated[1].heat_controls.temperature_tolerance_k *= 2.0;
        assert!(
            advance_liquid_surfaces_with_gas_and_heat(
                liquid,
                &gas,
                &source,
                &heated,
                face_controls(),
                GasLiquidChamberHeatSolve {
                    maximum_iterations: 128
                },
                gate,
            )
            .is_err()
        );
        assert_eq!(gas, original_gas);
        assert_eq!(liquid.flow().lattice().grid().state_bits(), bits);
        assert_eq!(liquid.total_enthalpy_j().unwrap(), original_h);
        assert_eq!(liquid.accepted_steps(), 0);
        assert_eq!(liquid.physical_time_s(), 0.0);
    });
}
