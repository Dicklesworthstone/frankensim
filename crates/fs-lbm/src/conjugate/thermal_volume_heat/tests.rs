//! Original-chart nonuniform heat, pressure work and conservative inventories.

#[cfg(feature = "implicit-enthalpy-coupling")]
mod implicit;

use super::*;
use crate::conjugate::{
    CoexistenceVolumeControls, CoexistenceVolumeInput, DensityForceControls,
    DensityMomentumReconstructionControls, DensityProjectionControls, MaterialAdvectionControls,
    MaterialDensityMomentumControls, MaterialMomentumControls, PiecewiseThermalVolumeControls,
    ThermalVolumeControls,
};
use fs_blake3::ContentHash;
use fs_material::phase::{EnthalpyPhaseKnot, EquilibriumPhaseState, SolidLiquidPhase};

const DT: f64 = 0.1;

fn controls() -> SuppliedHeatThermalVolumeControls {
    let projection = DensityProjectionControls {
        duration_s: DT,
        relative_tolerance: 1e-12,
        maximum_volume_defect_m3_s: 1e-10,
        maximum_energy_defect_j: 1e-8,
    };
    SuppliedHeatThermalVolumeControls {
        native: ThermalVolumeMomentumControls {
            volume: PiecewiseThermalVolumeControls {
                thermal: ThermalVolumeControls {
                    volume: CoexistenceVolumeControls {
                        transport: MaterialAdvectionControls {
                            projection,
                            max_cells: 100,
                            max_connections: 1000,
                            max_iterations: 1000,
                            relative_tolerance: 1e-12,
                            maximum_mass_defect_kg: 1e-10,
                            maximum_enthalpy_defect_j: 1e-9,
                        },
                        maximum_cell_volume_defect_m3: 1e-10,
                    },
                    max_volume_iterations: 40,
                    source_relaxation: 0.8,
                },
                max_mixing_intervals: 4000,
            },
            mechanics: MaterialDensityMomentumControls {
                reconstruction: DensityMomentumReconstructionControls {
                    max_cells: 100,
                    maximum_mass_defect_kg: 1e-10,
                    maximum_momentum_defect_kg_m_s: 1e-9,
                    maximum_energy_defect_j: 1e-8,
                },
                forcing: DensityForceControls {
                    max_cells: 100,
                    maximum_impulse_defect_kg_m_s: 1e-9,
                    maximum_energy_defect_j: 1e-8,
                },
                projection,
                maximum_feedback_momentum_defect_kg_m_s: 1e-9,
                maximum_feedback_angular_impulse_defect_kg_m2_s: 1e-9,
                maximum_feedback_energy_defect_j: 1e-8,
            },
            momentum: MaterialMomentumControls {
                max_iterations: 1000,
                relative_tolerance: 1e-12,
                maximum_momentum_defect_kg_m_s: 1e-9,
                maximum_energy_defect_j: 1e-8,
            },
        },
        maximum_heat_work_iterations: 80,
        work_relaxation: 0.5,
        maximum_absolute_local_heat_defect_j: 5e-10,
        maximum_total_heat_defect_j: 5e-9,
    }
}

fn chart(beta: f64, latent: bool) -> EquilibriumEnthalpyPhaseCurve {
    let samples = if latent {
        vec![
            (-10.0, 300.0, 0.0),
            (0.0, 400.0, 0.0),
            (40.0, 400.0, 1.0),
            (50.0, 500.0, 1.0),
        ]
    } else {
        vec![(0.0, 300.0, 1.0), (40.0, 340.0, 1.0)]
    };
    let knots = samples
        .into_iter()
        .map(|(h, t, f)| EnthalpyPhaseKnot {
            specific_enthalpy_j_kg: h,
            temperature_k: t,
            liquid_mass_fraction: f,
            bulk_density_kg_m3: 1.0 / (0.5 + beta * h),
        })
        .collect();
    if latent {
        EquilibriumEnthalpyPhaseCurve::try_new(ContentHash([0x93; 32]), knots).unwrap()
    } else {
        EquilibriumEnthalpyPhaseCurve::try_single_phase(
            ContentHash([0x93; 32]),
            SolidLiquidPhase::Liquid,
            knots,
        )
        .unwrap()
    }
}

fn input<'a>(
    mass: &'a [f64],
    states: &'a [Option<EquilibriumPhaseState>],
    guess: &'a [f64],
    momentum: &'a [[f64; 3]],
    force: &'a [[f64; 3]],
    outlet: usize,
) -> ThermalVolumeMomentumInput<'a> {
    let mut boundary_inflow = [None; 6];
    boundary_inflow[outlet] = states.last().copied().flatten();
    let mut boundary_velocity_m_s = [None; 6];
    boundary_velocity_m_s[outlet] = Some([0.0; 3]);
    ThermalVolumeMomentumInput {
        volume: CoexistenceVolumeInput {
            mass_kg: mass,
            phase_states: states,
            net_enthalpy_increment_j: guess,
            boundary_inflow,
        },
        momentum_kg_m_s: momentum,
        acceleration_m_s2: [0.0; 3],
        force_density_n_m3: force,
        boundary_velocity_m_s,
    }
}

fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 5e-8 * expected.abs().max(1.0),
        "{actual:.17e} != {expected:.17e}"
    );
}

#[test]
fn g1_spatial_heating_cooling_sensible_and_latent_expansion_pay_their_native_work() {
    for (beta, latent, datum) in [
        (0.01, false, 0.0),
        (-0.005, false, 0.0),
        (0.01, true, 0.0),
        (0.01, false, 1024.0),
    ] {
        let original_curve = chart(beta, latent);
        let curve = if datum == 0.0 {
            original_curve
        } else {
            let mut knots = original_curve.knots().to_vec();
            for knot in &mut knots {
                knot.specific_enthalpy_j_kg += datum;
            }
            EquilibriumEnthalpyPhaseCurve::try_single_phase(
                original_curve.material_card_identity(),
                SolidLiquidPhase::Liquid,
                knots,
            )
            .unwrap()
        };
        let h = [10.0 + datum, 20.0 + datum];
        let states = h.map(|h| Some(curve.state_at_specific_enthalpy(h).unwrap()));
        let mass = states.map(|state| state.unwrap().bulk_density_kg_m3());
        for sign in [-1.0, 1.0] {
            let net = [sign, 2.0 * sign];
            // Independent two-cell staggered pressure solution, with one
            // half-cell open reservoir and no initial motion/body force.
            let source = net.map(|q| beta * q / DT);
            let flux = [source[0], source[0] + source[1]];
            let p1 = 0.5 * mass[1] * flux[1] / DT;
            let pressure = [p1 + 0.5 * (mass[0] + mass[1]) * flux[0] / DT, p1];
            let work = core::array::from_fn::<_, 2, _>(|i| DT * pressure[i] * source[i]);
            let heat = core::array::from_fn::<_, 2, _>(|i| net[i] + work[i]);
            let h_star = core::array::from_fn::<_, 2, _>(|i| h[i] + net[i] / mass[i]);
            let kicked_velocity = [0.5 * flux[0], 0.5 * (flux[0] + flux[1])];
            let (expected_mass, expected_h, expected_p) = if flux[0] > 0.0 {
                let m0 = mass[0] / (1.0 + DT * flux[0]);
                let incoming = DT * flux[0] * m0;
                let m1 = (mass[1] + incoming) / (1.0 + DT * flux[1]);
                let h1 = (mass[1] * h_star[1] + incoming * h_star[0]) / (mass[1] + incoming);
                let p0 = m0 * kicked_velocity[0];
                let p1 = (mass[1] * kicked_velocity[1] + incoming * kicked_velocity[0])
                    / (1.0 + DT * flux[1]);
                ([m0, m1], [h_star[0], h1], [p0, p1])
            } else {
                let incoming = -DT * flux[1] * mass[1];
                let m1 = (mass[1] + incoming) / (1.0 - DT * flux[0]);
                let h1 = (mass[1] * h_star[1] + incoming * h[1]) / (mass[1] + incoming);
                let transfer = -DT * flux[0] * m1;
                let m0 = mass[0] + transfer;
                let h0 = (mass[0] * h_star[0] + transfer * h1) / m0;
                let p1 = mass[1] * kicked_velocity[1] / (1.0 - DT * flux[0]);
                let p0 = mass[0] * kicked_velocity[0] - DT * flux[0] * p1;
                ([m0, m1], [h0, h1], [p0, p1])
            };
            for axis in 0..3 {
                let mut shape = [1; 3];
                shape[axis] = 2;
                let domain = VoxelDomain::new(shape[0], shape[1], shape[2], 1.0).unwrap();
                let field = FlowField::quiescent(&domain);
                let original = field.clone();
                let outlet = 2 * axis + 1;
                let mut faces = [FlowFace::Wall; 6];
                faces[outlet] = FlowFace::Free;
                let zero = [[0.0; 3]; 2];
                let step = field
                    .advance_supplied_heat_thermal_volume_momentum(
                        &domain,
                        faces,
                        &curve,
                        input(&mass, &states, &heat, &zero, &zero, outlet),
                        &heat,
                        controls(),
                        &CancelGate::new_clock_free(),
                    )
                    .unwrap();
                let volume = &step.evolution.thermal.volume;
                let transport = &volume.transport;
                assert!(step.report.heat_work_iterations > 1);
                assert!(step.report.absolute_local_heat_defect_j <= 5e-10);
                assert!(step.report.total_heat_defect_j.abs() <= 5e-9);
                assert!(step.report.pressure_work_partition_defect_j.abs() < 1e-12);
                for c in 0..2 {
                    close(step.net_enthalpy_increment_j[c], net[c]);
                    close(step.cell_pressure_work_j[c], work[c]);
                    close(volume.volume_source_m3_s[c], source[c]);
                    close(transport.projection.pressure_correction_pa[c], pressure[c]);
                    close(transport.mass_kg[c], expected_mass[c]);
                    close(
                        transport.specific_enthalpy_j_kg[c] - datum,
                        expected_h[c] - datum,
                    );
                    close(
                        step.evolution.momentum.momentum_kg_m_s[c][axis],
                        expected_p[c],
                    );
                    let phase = volume.phase_states[c].unwrap();
                    close(transport.mass_kg[c] / phase.bulk_density_kg_m3(), 1.0);
                    assert_eq!(phase.phase_curve_identity(), curve.identity());
                    if latent {
                        assert_eq!(phase.temperature_k(), 400.0);
                    }
                }
                let carried = if flux[1] > 0.0 {
                    -DT * flux[1] * expected_mass[1] * expected_h[1]
                } else {
                    -DT * flux[1] * mass[1] * h[1]
                };
                close(transport.report.carried_enthalpy_j, carried);
                close(
                    expected_mass.iter().sum::<f64>() - mass.iter().sum::<f64>(),
                    transport.report.incoming_mass_kg - transport.report.outgoing_mass_kg,
                );
                close(
                    volume.report.material_enthalpy_change_j - carried + work.iter().sum::<f64>(),
                    heat.iter().sum(),
                );
                assert!(transport.projection.report.projection_loss_j > 0.0);
                // Pressure projection loss already has its own native channel.
                // It is not an additional source of heat or a second debit.
                close(step.report.source_pressure_work_j, work.iter().sum());
                assert_eq!(field, original);
            }
        }
    }
}

#[test]
fn g1_original_knot_mixing_volume_also_pays_its_pressure_work() {
    let curve = EquilibriumEnthalpyPhaseCurve::try_single_phase(
        ContentHash([0x94; 32]),
        SolidLiquidPhase::Liquid,
        [(0.0, 300.0, 2.0), (20.0, 400.0, 1.5), (40.0, 500.0, 1.0)]
            .map(|(h, t, rho)| EnthalpyPhaseKnot {
                specific_enthalpy_j_kg: h,
                temperature_k: t,
                liquid_mass_fraction: 1.0,
                bulk_density_kg_m3: rho,
            })
            .to_vec(),
    )
    .unwrap();
    let states = [10.0, 30.0].map(|h| Some(curve.state_at_specific_enthalpy(h).unwrap()));
    let mass = states.map(|state| state.unwrap().bulk_density_kg_m3());
    let domain = VoxelDomain::new(2, 1, 1, 1.0).unwrap();
    let field = FlowField::quiescent(&domain);
    let mut faces = [FlowFace::Wall; 6];
    faces[1] = FlowFace::Free;
    let zero = [[0.0; 3]; 2];
    let heat = [0.5, 0.3];
    let step = field
        .advance_supplied_heat_thermal_volume_momentum(
            &domain,
            faces,
            &curve,
            input(&mass, &states, &heat, &zero, &zero, 1),
            &heat,
            controls(),
            &CancelGate::new_clock_free(),
        )
        .unwrap();
    let thermal = &step.evolution.thermal;
    let pressure = &thermal.volume.transport.projection.pressure_correction_pa;
    let mixing_work: f64 = (0..2)
        .map(|i| DT * pressure[i] * thermal.mixing_volume_source_m3_s[i])
        .sum();
    assert!(mixing_work.abs() > 1e-5);
    assert!(thermal.source_iterations > 1);
    assert!(step.report.volume_source_iterations > step.report.heat_work_iterations);
    for c in 0..2 {
        let state = states[c].unwrap();
        let next = curve
            .advance_specific_energy(state, step.net_enthalpy_increment_j[c] / mass[c])
            .unwrap();
        let represented =
            mass[c] * (next.specific_enthalpy_j_kg() - state.specific_enthalpy_j_kg());
        close(
            represented + DT * pressure[c] * thermal.volume.volume_source_m3_s[c],
            heat[c],
        );
        close(
            thermal.volume.transport.mass_kg[c]
                / thermal.volume.phase_states[c].unwrap().bulk_density_kg_m3(),
            1.0,
        );
    }
    assert!(step.report.absolute_local_heat_defect_j <= 5e-10);
    assert!(step.report.total_heat_defect_j.abs() <= 5e-9);
}

#[test]
fn g3_zero_expansion_recovers_native_heat_transport_exactly_and_continuation_keeps_momentum() {
    let curve = chart(0.0, false);
    let domain = VoxelDomain::new(2, 1, 1, 1.0).unwrap();
    let field = FlowField::quiescent(&domain);
    let mut faces = [FlowFace::Wall; 6];
    faces[1] = FlowFace::Free;
    let states = [10.0, 20.0].map(|h| Some(curve.state_at_specific_enthalpy(h).unwrap()));
    let mass = [2.0; 2];
    let zero = [[0.0; 3]; 2];
    let momentum = [[0.2, 0.0, 0.0], [0.4, 0.0, 0.0]];
    let heat = [0.5, -0.25];
    let gate = CancelGate::new_clock_free();
    let data = input(&mass, &states, &heat, &momentum, &zero, 1);
    let expected = field
        .advance_piecewise_thermal_volume_momentum(
            &domain,
            faces,
            &curve,
            data,
            controls().native,
            &gate,
        )
        .unwrap();
    let actual = field
        .advance_supplied_heat_thermal_volume_momentum(
            &domain,
            faces,
            &curve,
            data,
            &heat,
            controls(),
            &gate,
        )
        .unwrap();
    assert_eq!(actual.evolution, expected);
    assert_eq!(actual.report.heat_work_iterations, 1);
    assert_eq!(actual.report.source_pressure_work_j, 0.0);
    assert!(
        actual
            .evolution
            .thermal
            .volume
            .transport
            .projection
            .report
            .projection_loss_j
            > 0.0
    );
    let v = &actual.evolution.thermal.volume;
    let next_heat = [0.0; 2];
    let next = input(
        &v.transport.mass_kg,
        &v.phase_states,
        &next_heat,
        &actual.evolution.momentum.momentum_kg_m_s,
        &zero,
        1,
    );
    let continued = v
        .transport
        .projection
        .field
        .advance_supplied_heat_thermal_volume_momentum(
            &domain,
            faces,
            &curve,
            next,
            &next_heat,
            controls(),
            &gate,
        )
        .unwrap();
    assert_eq!(continued.report.source_pressure_work_j, 0.0);
    for c in 0..2 {
        close(
            continued.evolution.thermal.volume.transport.mass_kg[c],
            v.transport.mass_kg[c],
        );
        close(
            continued
                .evolution
                .thermal
                .volume
                .transport
                .specific_enthalpy_j_kg[c],
            v.transport.specific_enthalpy_j_kg[c],
        );
    }
}

#[test]
fn g4_unsupported_source_unanchored_work_and_cancel_refuse_without_publication() {
    let curve = chart(0.01, false);
    let domain = VoxelDomain::new(2, 1, 1, 1.0).unwrap();
    let field = FlowField::quiescent(&domain);
    let original = field.clone();
    let states = [10.0, 20.0].map(|h| Some(curve.state_at_specific_enthalpy(h).unwrap()));
    let mass = states.map(|state| state.unwrap().bulk_density_kg_m3());
    let zero = [[0.0; 3]; 2];
    let heat = [1.0, 2.0];
    let data = input(&mass, &states, &heat, &zero, &zero, 1);
    let mut faces = [FlowFace::Wall; 6];
    faces[1] = FlowFace::Free;
    let gate = CancelGate::new_clock_free();
    let run = |data, heat: &[f64], limits, gate: &CancelGate| {
        field.advance_supplied_heat_thermal_volume_momentum(
            &domain, faces, &curve, data, heat, limits, gate,
        )
    };
    let expected = run(data, &heat, controls(), &gate).unwrap();
    let mut short = controls();
    short.maximum_heat_work_iterations = 1;
    assert!(run(data, &heat, short, &gate).is_err());
    let unsupported = [1e6; 2];
    assert!(
        run(
            input(&mass, &states, &unsupported, &zero, &zero, 1),
            &heat,
            controls(),
            &gate
        )
        .is_err()
    );
    assert!(run(data, &[f64::NAN, 2.0], controls(), &gate).is_err());
    let mut late = controls();
    late.native.momentum.max_iterations = 0;
    assert!(run(data, &heat, late, &gate).is_err());
    let no_heat = [0.0; 2];
    assert!(
        matches!(field.advance_supplied_heat_thermal_volume_momentum(
        &domain, [FlowFace::Wall; 6], &curve,
        input(&mass, &states, &no_heat, &zero, &zero, 1), &no_heat, controls(), &gate,
    ), Err(ChtError::InvalidInput { field: "flow.supplied_heat_thermal_volume", reason })
        if reason.contains("reference-pressure Free opening"))
    );
    assert_eq!(run(data, &heat, controls(), &gate).unwrap(), expected);
    let cancelled = CancelGate::new_clock_free();
    cancelled.request();
    assert!(matches!(
        run(data, &heat, controls(), &cancelled),
        Err(ChtError::Cancelled)
    ));
    assert_eq!(field, original);
    assert_eq!(heat, [1.0, 2.0]);
    assert_eq!(
        states,
        [10.0, 20.0].map(|h| Some(curve.state_at_specific_enthalpy(h).unwrap()))
    );
}
