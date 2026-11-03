//! Native finite C/L storage and the full 3-D terminal field share one heat endpoint.
use super::*;
use fs_couple::circuit_enthalpy::conduction::spatial::circuit::SpatialConductionCircuit;
use fs_phs::circuit::{Branch, CircuitGraph, MutualInductance};

fn response_controls() -> SpatialConductanceControls {
    SpatialConductanceControls {
        maximum_terminals: 3,
        maximum_field_solves: 2,
        maximum_response_values: 1000,
        maximum_dense_work: 1000,
        reference_voltage_v: 0.5,
        coefficient_tolerance_s: 1e-10,
    }
}

#[test]
fn native_capacitor_inductor_and_source_exchange_actual_multiterminal_joule_work() {
    for sign in [1.0, -1.0] {
        let (a, ca) = material(1.0, 1.0, 0.002);
        let (b, cb) = material(3.0, 2.0, 0.0);
        let cards = [a, b];
        let curves = [ca, cb];
        let context = point(300.0);
        let key = PropertyKey::new("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS);
        let cells = bindings(&cards, &curves, &context, &key);
        let mut thermal = initial(&curves, 999.0, None);
        let conductors = cards.each_ref().map(|card| {
            resolve_conductor_state_with_basis(
                card,
                &context,
                MaterialPropertySelection::SingleClaimOnly,
                ElectricalPropertyBasis::Resistivity,
            )
            .unwrap()
        });
        let (mut positions, mut tets) = slab(&[0.0, 0.5, 1.0]);
        let (other, other_tets) = slab(&[1.0, 1.5, 2.0]);
        let offset = positions.len() as u32;
        positions.extend(other);
        tets.extend(other_tets.into_iter().map(|t| t.map(|v| v + offset)));
        let field = SpatialConduction::new(
            positions,
            tets,
            (0..24).map(|i| conductors[i / 12].clone()).collect(),
            limits(),
            &mut || false,
        )
        .unwrap();
        let faces = [0.0, 1.0, 2.0].map(|x| {
            field
                .exterior_faces()
                .filter(|&f| {
                    field.mesh().faces[f].map(|v| field.positions_m()[v as usize][0]) == [x; 3]
                })
                .collect::<Vec<_>>()
        });
        let regions = (0..24)
            .map(|tetrahedron| SpatialConductionRegion {
                tetrahedron,
                cell: cells[tetrahedron / 12].cell,
            })
            .collect::<Vec<_>>();
        let network = SpatialConductionThermalNetwork {
            cells: &cells,
            links: &LINKS,
            regions: &regions,
            volume_tolerance_m3: 1e-12,
        };
        // C=1 F at the middle electrode; L=2 H at the right electrode.
        // The spatial body alone supplies both dissipation paths.
        let graph = CircuitGraph {
            node_count: 3,
            branches: vec![
                (1, 0, Branch::Capacitor { farads: 1.0 }),
                (2, 0, Branch::Inductor { henries: 2.0 }),
                (1, 0, Branch::CurrentSource { port: 0 }),
            ],
            transformers: vec![],
        };
        let mut charge = sign * 4.0;
        let mut flux = sign * 0.5;
        let ports = [sign];
        let mut circuit = SpatialConductionCircuit::new(
            field,
            graph,
            faces.to_vec(),
            vec![0, 1, 2],
            0.0,
            vec![flux],
            vec![charge],
            &ports,
            response_controls(),
            field_controls(),
            || false,
        )
        .unwrap();
        let mut theta = [0.0; 2];
        let dt = 0.5;
        for interval in 1..=2 {
            // Independent midpoint RLC elimination at each proposed final RA:
            // (1+h/RA)*v+h*i=q_old+h*u, -h/2*v+(1+3h/2)*i=phi_old/2.
            // Heat then obeys two original BE balances, Cthermal=[1,2], G=1.
            let endpoint = |d_a: f64| {
                let rho = 1.0 + 0.002 * d_a;
                let h = 0.5 * dt;
                let a = 1.0 + h / rho;
                let d = 1.0 + 1.5 * h;
                let determinant = a * d + 0.5 * h * h;
                let rhs = charge + h * ports[0];
                let v = (d * rhs - h * flux / 2.0) / determinant;
                let i = (0.5 * h * rhs + a * flux / 2.0) / determinant;
                let q_a = dt * v * v / rho;
                let q_b = dt * 3.0 * i * i;
                let d_b = (2.0 * theta[1] + q_b + dt * d_a) / (2.0 + dt);
                (
                    v,
                    i,
                    q_a,
                    q_b,
                    d_b,
                    (1.0 + dt) * d_a - dt * d_b - theta[0] - q_a,
                )
            };
            let (mut lo, mut hi) = (0.0, 50.0);
            assert!(endpoint(lo).5 < 0.0 && endpoint(hi).5 > 0.0);
            for _ in 0..90 {
                let mid = 0.5 * (lo + hi);
                if endpoint(mid).5 < 0.0 {
                    lo = mid;
                } else {
                    hi = mid;
                }
            }
            let d_a = 0.5 * (lo + hi);
            let (v, i, q_a, q_b, d_b, _) = endpoint(d_a);
            let final_charge = 2.0 * v - charge;
            let final_flux = 4.0 * i - flux;
            let before_thermal = thermal.energy_j();
            let next = network
                .step_circuit_candidate(
                    &thermal,
                    interval as f64 * dt,
                    &circuit,
                    &ports,
                    &[],
                    response_controls(),
                    limits(),
                    field_controls(),
                    config(),
                    &[0.0; 2],
                    || false,
                )
                .unwrap();
            near(
                next.thermal().temperature_k(CellId(30)).unwrap(),
                300.0 + d_a,
            );
            near(
                next.thermal().temperature_k(CellId(10)).unwrap(),
                300.0 + d_b,
            );
            near(next.circuit().charges_c()[0], final_charge);
            near(next.circuit().fluxes_wb()[0], final_flux);
            let step = &next.report().electrical;
            near(step.terminal_voltage_v[1], v);
            near(step.terminal_voltage_v[2], v - 3.0 * i);
            near(step.terminal_current_a[0], -v / (1.0 + 0.002 * d_a));
            near(step.terminal_current_a[2], -i);
            near(
                step.spatial.field.cell_power_w[..12].iter().sum::<f64>() * dt,
                q_a,
            );
            near(
                step.spatial.field.cell_power_w[12..].iter().sum::<f64>() * dt,
                q_b,
            );
            near(step.native.dissipated, q_a + q_b);
            near(step.spatial.heat_j, q_a + q_b);
            near(step.spatial.supplied_work_j, q_a + q_b);
            near(step.supplied_work_j, dt * ports[0] * v);
            let storage_change = 0.5 * (final_charge * final_charge - charge * charge)
                + 0.25 * (final_flux * final_flux - flux * flux);
            near(step.stored_change_j, storage_change);
            near(step.native.delta_h, storage_change);
            near(step.energy_residual_j, 0.0);
            near(next.thermal().energy_j() - before_thermal, q_a + q_b);
            near(next.report().total_energy_residual_j, 0.0);
            assert_eq!(next.report().circuit_trials % 3, 0);
            assert!(next.report().circuit_trials >= 3);
            if interval == 1 {
                let mut invalid = config();
                invalid.circuit.power_tolerance_w = f64::NAN;
                assert!(
                    network
                        .step_circuit_candidate(
                            &thermal,
                            dt,
                            &circuit,
                            &ports,
                            &[],
                            response_controls(),
                            limits(),
                            field_controls(),
                            invalid,
                            &[0.0; 2],
                            || false,
                        )
                        .is_err()
                );
                let mut short = config();
                short.maximum_circuit_trials = 2;
                assert!(matches!(
                    network.step_circuit_candidate(
                        &thermal,
                        dt,
                        &circuit,
                        &ports,
                        &[],
                        response_controls(),
                        limits(),
                        field_controls(),
                        short,
                        &[0.0; 2],
                        || false,
                    ),
                    Err(CircuitEnthalpyError::CircuitBudget)
                ));
                assert!(matches!(
                    network.step_circuit_candidate(
                        &thermal,
                        dt,
                        &circuit,
                        &ports,
                        &[],
                        response_controls(),
                        limits(),
                        field_controls(),
                        config(),
                        &[0.0; 2],
                        || true,
                    ),
                    Err(CircuitEnthalpyError::Cancelled)
                ));
                near(circuit.charges_c()[0], charge);
                near(circuit.fluxes_wb()[0], flux);
                assert_eq!(thermal.time_s(), 0.0);
            }
            theta = [d_a, d_b];
            charge = final_charge;
            flux = final_flux;
            let (accepted_heat, accepted_circuit, _) = next.into_parts();
            thermal = accepted_heat;
            circuit = accepted_circuit;
        }
    }
}

#[test]
fn transformer_leakage_power_heats_the_original_spatial_body_and_continues() {
    for orientation in [1.0, -1.0] {
        let (a, ca) = material(1.0, 1.0, 0.0);
        let (b, cb) = material(1.0, 1.0, 0.0);
        let cards = [a, b];
        let curves = [ca, cb];
        let context = point(300.0);
        let key = PropertyKey::new("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS);
        let cells = bindings(&cards, &curves, &context, &key);
        let mut thermal = initial(&curves, 999.0, None);
        let conductors = cards.each_ref().map(|card| {
            resolve_conductor_state_with_basis(
                card,
                &context,
                MaterialPropertySelection::SingleClaimOnly,
                ElectricalPropertyBasis::Resistivity,
            )
            .unwrap()
        });
        // Two equal 1 m³ thermal cells make one 2 ohm spatial conductor.
        let (positions, tets) = slab(&[0.0, 0.5, 1.0, 1.5, 2.0]);
        let field = SpatialConduction::new(
            positions,
            tets,
            (0..24).map(|i| conductors[i / 12].clone()).collect(),
            limits(),
            &mut || false,
        )
        .unwrap();
        let faces = electrodes(&field, 2.0);
        let regions = (0..24)
            .map(|tetrahedron| SpatialConductionRegion {
                tetrahedron,
                cell: cells[tetrahedron / 12].cell,
            })
            .collect::<Vec<_>>();
        let network = SpatialConductionThermalNetwork {
            cells: &cells,
            links: &LINKS,
            regions: &regions,
            volume_tolerance_m3: 1e-12,
        };
        let graph = CircuitGraph {
            node_count: 3,
            branches: vec![
                (1, 0, Branch::VoltageSource { port: 0 }),
                (1, 0, Branch::Inductor { henries: 4.0 }),
                if orientation > 0.0 {
                    (2, 0, Branch::Inductor { henries: 1.0 })
                } else {
                    (0, 2, Branch::Inductor { henries: 1.0 })
                },
            ],
            transformers: vec![],
        };
        let links = vec![MutualInductance {
            branches: (1, 2),
            henries: orientation * 1.5,
        }];
        let mut circuit = SpatialConductionCircuit::new_with_couplings(
            field,
            graph,
            faces.to_vec(),
            vec![2, 0],
            Vec::new(),
            links.clone(),
            0.0,
            vec![0.0; 2],
            vec![],
            &[1.0],
            response_controls(),
            field_controls(),
            || false,
        )
        .unwrap();
        let (mut phi1, mut i2, mut heat, mut supply) = (0.0, 0.0, 0.0, 0.0);
        let dt = 0.2;
        for interval in 1..=3 {
            // V1=1, M=1.5: (L2-M²/L1) di2/dt = -2 i2 - M/L1.
            // Independent midpoint solution retains finite leakage L_eff=7/16.
            let h = dt / 2.0;
            let leakage = 7.0 / 16.0;
            let mid_i2 = (i2 - h * 0.375 / leakage) / (1.0 + h * 2.0 / leakage);
            let mid_i1 = (phi1 + h - 1.5 * mid_i2) / 4.0;
            let step_heat = dt * 2.0 * mid_i2 * mid_i2;
            phi1 += dt;
            i2 = 2.0 * mid_i2 - i2;
            let i1 = (phi1 - 1.5 * i2) / 4.0;
            let phi2 = 1.5 * i1 + i2;
            heat += step_heat;
            supply += dt * mid_i1;
            let next = network
                .step_circuit_candidate(
                    &thermal,
                    interval as f64 * dt,
                    &circuit,
                    &[1.0],
                    &[],
                    response_controls(),
                    limits(),
                    field_controls(),
                    config(),
                    &[0.0; 2],
                    || false,
                )
                .unwrap();
            let step = &next.report().electrical;
            near(step.terminal_voltage_v[0], -2.0 * mid_i2);
            near(step.terminal_current_a[0], -mid_i2);
            near(step.spatial.field.power_w, 2.0 * mid_i2 * mid_i2);
            near(step.heat_j, step_heat);
            near(step.native.dissipated, step_heat);
            near(step.supplied_work_j, dt * mid_i1);
            near(step.energy_residual_j, 0.0);
            near(next.report().total_energy_residual_j, 0.0);
            near(next.circuit().fluxes_wb()[0], phi1);
            near(next.circuit().fluxes_wb()[1], orientation * phi2);
            assert_eq!(next.circuit().mutual_inductance(), links);
            for cell in [CellId(30), CellId(10)] {
                near(
                    next.thermal().temperature_k(cell).unwrap(),
                    300.0 + heat / 2.0,
                );
            }
            let stored = 0.5 * (4.0 * i1 * i1 + 3.0 * i1 * i2 + i2 * i2);
            near(step.stored_after_j, stored);
            near(supply, stored + heat);
            let (accepted_heat, accepted_circuit, _) = next.into_parts();
            thermal = accepted_heat;
            circuit = accepted_circuit;
        }
        assert!(heat > 0.0);
    }
}
