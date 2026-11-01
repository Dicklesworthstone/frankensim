use super::*;
use fs_bem::near_field::FirstOrder;
use fs_bem::observation_transfers::{HarmonicField, HarmonicWorkBudget};
use fs_duct::network::{
    CoupledSourceBank, NetworkCoherenceBudget, NetworkLoad, NetworkNode, NetworkSourceQuantity,
    PortExcitation, PortSource, network_drive_source_bank, network_drive_sources,
};

mod angular;

const FIELD_BUDGET: HarmonicWorkBudget = HarmonicWorkBudget {
    maximum_products: 100_000,
};

fn two_source_model() -> RadiatingNetwork {
    let mut network = graph();
    network.nodes.push(NetworkNode {
        id: "second-inlet".into(),
        boundary: NodeBoundary::Port {
            label: "aux".into(),
        },
    });
    // Identical source arms give an independent antisymmetric-mode reference,
    // even though the full exterior opening matrix remains mutually coupled.
    let mut arm = network.edges[0].clone();
    arm.id = "second-stem".into();
    arm.from = 4;
    network.edges.push(arm);
    let (surface, patches) = exterior();
    let driving = NetworkDrivingPoint::new_coupled(
        network,
        "inlet".into(),
        vec![PortTermination {
            label: "aux".into(),
            load: NetworkLoad::Closed,
        }],
        patches.iter().map(|p| p.port.clone()).collect(),
        schema(),
    )
    .unwrap();
    RadiatingNetwork::new(
        driving,
        surface,
        patches,
        gas(),
        receivers(),
        0.001,
        ContentHash([12; 32]),
        budget(),
        &CancelGate::new_clock_free(),
    )
    .unwrap()
}

fn bank(omega: f64) -> CoupledSourceBank {
    let (self_z, mutual) = (C64::new(2.0e6, -1.0e6), C64::new(-0.4e6, 0.3e6));
    CoupledSourceBank {
        impedance: CoupledTermination {
            omega,
            labels: vec!["inlet".into(), "aux".into()],
            z: vec![self_z, mutual, mutual, self_z],
            source_id: "original-coupled-generators".into(),
        },
        open_circuit_pressures_pa: vec![C64::new(0.8, -0.2), C64::new(-0.3, 0.7)],
    }
}

fn real_close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() < 3e-9 * expected.abs().max(1e-16),
        "actual {actual:e}, expected {expected:e}"
    );
}

// A freshly assembled physical panel velocity is independent of the retained
// unit-opening source bank and its composed source-to-receiver transfers.
fn direct_fields(model: &RadiatingNetwork, omega: f64, outward: &[C64]) -> Vec<HarmonicField> {
    let mut velocity = vec![C64::ZERO; model.surface().areas().len()];
    for ((patch, &area), &flow) in model
        .patches()
        .iter()
        .zip(model.opening_areas_m2())
        .zip(outward)
    {
        for &panel in &patch.panels {
            velocity[panel] = flow.scale(area.recip());
        }
    }
    let medium = Medium::from(&model.gas());
    let solution = helmholtz::solve_radiation(
        model.surface(),
        omega / medium.sound_speed,
        medium,
        &velocity,
        Formulation::BurtonMiller,
    )
    .unwrap();
    let euler = C64::new(0.0, -1.0 / (omega * medium.density));
    laplace::observe_frequency(model.surface(), &solution, model.receivers(), 0.001)
        .unwrap()
        .into_iter()
        .map(|field| HarmonicField {
            pressure_pa: field.pressure_pa,
            particle_velocity_m_s: field.gradient_pa_per_m.map(|gradient| euler * gradient),
        })
        .collect()
}

#[test]
fn full_generator_bank_matches_independent_boundary_fields_and_signed_work() {
    let model = two_source_model();
    let omega = core::f64::consts::TAU * 1200.0;
    let bank = bank(omega);
    let gate = CancelGate::new_clock_free();
    let prepared = model
        .prepare_sources(omega, LossModel::Lossless, &[], Some(&bank), &[], &gate)
        .unwrap();
    assert_eq!(
        prepared
            .source_coordinates()
            .iter()
            .map(|s| s.label.as_str())
            .collect::<Vec<_>>(),
        ["aux", "inlet"]
    );
    assert!(
        prepared
            .source_coordinates()
            .iter()
            .all(|s| s.quantity == NetworkSourceQuantity::OpenCircuitPressure)
    );
    assert_eq!(
        prepared.boundary().source_id,
        format!("fs-bem:network-openings:{}", model.boundary_hash())
    );
    let mut result = prepared.response_buffer();
    prepared
        .drive_into(prepared.default_excitations(), &mut result, &gate)
        .unwrap();
    let network = result.network().unwrap();
    let direct = network_drive_source_bank(
        model.driving_point().network(),
        &model.gas(),
        omega,
        LossModel::Lossless,
        &[],
        &bank,
        &[],
        Some(prepared.boundary()),
    )
    .unwrap();
    for (&a, &b) in network
        .fields
        .port_pressures
        .iter()
        .zip(&direct.fields.port_pressures)
    {
        close(a, b);
    }
    for (&a, &b) in network
        .fields
        .port_inflows
        .iter()
        .zip(&direct.fields.port_inflows)
    {
        close(a, b);
    }
    let observed = direct_fields(&model, omega, result.outward_flows_m3_s().unwrap());
    for (a, b) in result.receiver_fields().unwrap().iter().zip(&observed) {
        close(a.pressure_pa, b.pressure_pa);
        for i in 0..3 {
            close(a.particle_velocity_m_s[i], b.particle_velocity_m_s[i]);
        }
    }
    let exterior = result.exterior_surface_power_w().unwrap();
    real_close(network.fields.external_load_power_w, exterior);
    real_close(
        network.generator_power_w,
        network.source_dissipation_w + exterior,
    );
    assert!(network.source_dissipation_w > 0.0 && exterior > 0.0);
    let address = result.receiver_fields().unwrap().as_ptr();
    prepared
        .drive_into(
            &[C64::new(1.0, 0.4), C64::new(-0.2, 0.3)],
            &mut result,
            &gate,
        )
        .unwrap();
    assert_eq!(address, result.receiver_fields().unwrap().as_ptr());
}

#[test]
fn antisymmetric_generators_cancel_radiation_but_retain_analytic_internal_loss() {
    let model = two_source_model();
    let omega = core::f64::consts::TAU * 1200.0;
    let bank = bank(omega);
    let gate = CancelGate::new_clock_free();
    let prepared = model
        .prepare_sources(omega, LossModel::Lossless, &[], Some(&bank), &[], &gate)
        .unwrap();
    let mut common = prepared.response_buffer();
    let mut opposed = prepared.response_buffer();
    prepared
        .drive_into(&[C64::ONE; 2], &mut common, &gate)
        .unwrap();
    prepared
        .drive_into(&[C64::ONE, -C64::ONE], &mut opposed, &gate)
        .unwrap();
    let gas = model.gas();
    let zc = gas.density * gas.sound_speed / (core::f64::consts::PI * 0.006_f64.powi(2));
    let shorted_arm = C64::new(0.0, -zc * (omega * 0.13 / gas.sound_speed).tan());
    let differential = bank.impedance.z[0] - bank.impedance.z[1];
    let flow = C64::ONE / (differential + shorted_arm);
    let expected_loss = differential.re * flow.norm_sq();
    let network = opposed.network().unwrap();
    let aux = prepared.ports().iter().position(|p| p == "aux").unwrap();
    close(network.fields.port_inflows[aux], flow);
    real_close(network.source_dissipation_w, expected_loss);
    real_close(network.generator_power_w, expected_loss);
    assert!(expected_loss > 0.0);
    for (cancelled, bright) in opposed
        .receiver_fields()
        .unwrap()
        .iter()
        .zip(common.receiver_fields().unwrap())
    {
        assert!(cancelled.pressure_pa.abs() < 1e-10 * bright.pressure_pa.abs());
    }
    assert!(
        opposed.exterior_surface_power_w().unwrap().abs()
            < 1e-18 * common.exterior_surface_power_w().unwrap()
    );
    let mut statistics = prepared.statistics_buffer();
    prepared
        .independent_into(
            &[1.0; 2],
            NetworkCoherenceBudget::default(),
            FIELD_BUDGET,
            &mut statistics,
            &gate,
        )
        .unwrap();
    for (actual, bright) in statistics
        .receiver_fields()
        .unwrap()
        .iter()
        .zip(common.receiver_fields().unwrap())
    {
        // Independent covariance is half common plus half differential. The
        // differential mode is dark outside but dissipates inside the sources.
        real_close(
            actual.mean_square_pressure_pa2,
            0.25 * bright.pressure_pa.norm_sq(),
        );
    }
    let work = statistics.network().unwrap().work().unwrap();
    real_close(
        work.source_dissipation_w,
        0.5 * (common.network().unwrap().source_dissipation_w + expected_loss),
    );
    real_close(
        work.generator_power_w,
        work.source_dissipation_w + work.external_load_power_w,
    );
    let reference = statistics.receiver_fields().unwrap().to_vec();
    prepared
        .coherence_into(
            &[C64::ONE, C64::ZERO, C64::ZERO, C64::ONE],
            2,
            NetworkCoherenceBudget::default(),
            FIELD_BUDGET,
            &mut statistics,
            &gate,
        )
        .unwrap();
    assert_eq!(statistics.receiver_fields().unwrap(), reference);
    let mut microphones = prepared.microphones(&[FirstOrder::default(); 2]).unwrap();
    let independent = microphones
        .independent(&[1.0; 2], FIELD_BUDGET, &gate)
        .unwrap();
    for i in 0..2 {
        for j in 0..2 {
            close(
                independent[i * 2 + j],
                (common.receiver_fields().unwrap()[i].pressure_pa
                    * common.receiver_fields().unwrap()[j].pressure_pa.conj())
                .scale(0.25),
            );
        }
    }
}

#[test]
fn mixed_pressure_flow_covariance_reaches_full_near_field_and_microphone_cross_spectra() {
    let model = two_source_model();
    let omega = core::f64::consts::TAU * 1200.0;
    let sources = [
        PortSource {
            label: "inlet".into(),
            excitation: PortExcitation::VolumeFlow(C64::new(1e-6, -0.4e-6)),
        },
        PortSource {
            label: "aux".into(),
            excitation: PortExcitation::Pressure {
                open_circuit_pressure_pa: C64::new(0.6, 0.2),
                impedance_pa_s_m3: C64::new(2e6, -1e6),
            },
        },
    ];
    let gate = CancelGate::new_clock_free();
    let prepared = model
        .prepare_sources(omega, LossModel::Bessel, &sources, None, &[], &gate)
        .unwrap();
    assert_eq!(
        prepared.source_coordinates()[0].quantity,
        NetworkSourceQuantity::OpenCircuitPressure
    );
    assert_eq!(
        prepared.source_coordinates()[1].quantity,
        NetworkSourceQuantity::VolumeFlow
    );
    let direct = network_drive_sources(
        model.driving_point().network(),
        &model.gas(),
        omega,
        LossModel::Bessel,
        &sources,
        &[],
        Some(prepared.boundary()),
    )
    .unwrap();
    let mut coherent = prepared.response_buffer();
    prepared
        .drive_into(prepared.default_excitations(), &mut coherent, &gate)
        .unwrap();
    for (&a, &b) in coherent
        .network()
        .unwrap()
        .fields
        .port_inflows
        .iter()
        .zip(&direct.fields.port_inflows)
    {
        close(a, b);
    }
    // Nontrivial imaginary cross covariance in Pa*m³/s: the source coordinates
    // are neither normalized into watts nor replaced by independent openings.
    let factors = [
        C64::new(0.7, -0.2),
        C64::new(1.0e-6, 0.5e-6),
        C64::new(-0.1, 0.4),
        C64::new(-0.3e-6, 0.2e-6),
    ];
    let mut basis_fields = Vec::new();
    for source in 0..2 {
        let mut vector = [C64::ZERO; 2];
        vector[source] = C64::ONE;
        prepared.drive_into(&vector, &mut coherent, &gate).unwrap();
        basis_fields.push(coherent.receiver_fields().unwrap().to_vec());
    }
    let covariance: Vec<_> = (0..2)
        .flat_map(|i| {
            (0..2).map(move |j| {
                factors[i] * factors[j].conj() + factors[2 + i] * factors[2 + j].conj()
            })
        })
        .collect();
    let contract = |left: &[C64], right: &[C64]| {
        (0..2)
            .flat_map(|i| (0..2).map(move |j| (i, j)))
            .fold(C64::ZERO, |value, (i, j)| {
                value + (left[i] * covariance[i * 2 + j] * right[j].conj()).scale(0.5)
            })
    };
    let mut statistics = prepared.statistics_buffer();
    prepared
        .coherence_into(
            &factors,
            2,
            NetworkCoherenceBudget::default(),
            FIELD_BUDGET,
            &mut statistics,
            &gate,
        )
        .unwrap();
    for (receiver, actual) in statistics.receiver_fields().unwrap().iter().enumerate() {
        let pressure: Vec<_> = basis_fields
            .iter()
            .map(|s| s[receiver].pressure_pa)
            .collect();
        real_close(
            actual.mean_square_pressure_pa2,
            contract(&pressure, &pressure).re,
        );
        for axis in 0..3 {
            let velocity: Vec<_> = basis_fields
                .iter()
                .map(|s| s[receiver].particle_velocity_m_s[axis])
                .collect();
            close(
                actual.complex_intensity_w_m2[axis],
                contract(&pressure, &velocity),
            );
            for other in 0..3 {
                let other_velocity: Vec<_> = basis_fields
                    .iter()
                    .map(|s| s[receiver].particle_velocity_m_s[other])
                    .collect();
                close(
                    actual.velocity_cross_m2_s2[axis][other],
                    contract(&velocity, &other_velocity),
                );
            }
        }
    }
    let patterns = [
        FirstOrder::new(0.5, [1.0, 0.0, 0.0]).unwrap(),
        FirstOrder::new(0.2, [0.0, 1.0, 0.0]).unwrap(),
    ];
    let microphone_columns: Vec<Vec<_>> = patterns
        .iter()
        .enumerate()
        .map(|(receiver, &pattern)| {
            basis_fields
                .iter()
                .map(|s| {
                    pattern
                        .observe(
                            s[receiver].pressure_pa,
                            s[receiver].particle_velocity_m_s,
                            Medium::from(&model.gas()),
                        )
                        .unwrap()
                })
                .collect()
        })
        .collect();
    let mut microphones = prepared.microphones(&patterns).unwrap();
    let spectrum = microphones
        .coherence(&factors, 2, FIELD_BUDGET, &gate)
        .unwrap();
    for i in 0..2 {
        for j in 0..2 {
            close(
                spectrum[i * 2 + j],
                contract(&microphone_columns[i], &microphone_columns[j]),
            );
        }
    }
    assert!(spectrum[1].im.abs() > 1e-12 * spectrum[0].re);
    let work = statistics.network().unwrap().work().unwrap();
    real_close(
        work.generator_power_w,
        work.source_dissipation_w + work.network_power_w + work.external_load_power_w,
    );
    real_close(
        statistics.exterior_surface_power_w().unwrap(),
        work.external_load_power_w,
    );
}

#[test]
fn original_single_source_parity_and_transactional_refusals() {
    let model = model(gas(), receivers());
    let gate = CancelGate::new_clock_free();
    let sources = [PortSource {
        label: "inlet".into(),
        excitation: PortExcitation::VolumeFlow(C64::new(1e-6, -2e-6)),
    }];
    let omega = core::f64::consts::TAU * 1200.0;
    let legacy = model
        .solve(
            omega,
            LossModel::Lossless,
            match sources[0].excitation {
                PortExcitation::VolumeFlow(x) => x,
                _ => unreachable!(),
            },
            &gate,
        )
        .unwrap();
    let prepared = model
        .prepare_sources(omega, LossModel::Lossless, &sources, None, &[], &gate)
        .unwrap();
    let mut result = prepared.response_buffer();
    prepared
        .drive_into(prepared.default_excitations(), &mut result, &gate)
        .unwrap();
    for (new, old) in result
        .receiver_fields()
        .unwrap()
        .iter()
        .zip(legacy.receiver_fields)
    {
        close(new.pressure_pa, old.pressure_pa);
    }
    real_close(
        result.exterior_surface_power_w().unwrap(),
        legacy.exterior_surface_power_w,
    );
    let other = model
        .prepare_sources(omega, LossModel::Lossless, &sources, None, &[], &gate)
        .unwrap();
    assert!(
        other
            .drive_into(other.default_excitations(), &mut result, &gate)
            .is_err()
    );
    assert!(result.network().is_none() && result.receiver_fields().is_none());
    let stopped = CancelGate::new_clock_free();
    stopped.request();
    assert!(matches!(
        prepared.drive_into(prepared.default_excitations(), &mut result, &stopped),
        Err(ExteriorNetworkError::Cancelled)
    ));
    let mut statistics = prepared.statistics_buffer();
    for (rank, factors, budget) in [
        (33, vec![C64::ONE; 33], FIELD_BUDGET),
        (1, vec![C64::new(f64::NAN, 0.0)], FIELD_BUDGET),
        (
            1,
            vec![C64::ONE],
            HarmonicWorkBudget {
                maximum_products: 1,
            },
        ),
    ] {
        assert!(
            prepared
                .coherence_into(
                    &factors,
                    rank,
                    NetworkCoherenceBudget::default(),
                    budget,
                    &mut statistics,
                    &gate
                )
                .is_err()
        );
        assert!(
            !statistics.is_valid()
                && statistics.network().is_none()
                && statistics.receiver_fields().is_none()
        );
    }
    assert!(
        model
            .prepare_sources(
                omega,
                LossModel::Lossless,
                &sources,
                None,
                &[PortTermination {
                    label: "right".into(),
                    load: NetworkLoad::Closed
                }],
                &gate
            )
            .is_err()
    );
    assert!(
        model
            .prepare_sources(omega, LossModel::Lossless, &sources, None, &[], &stopped)
            .is_err()
    );
}

#[test]
fn loaded_source_observations_survive_exact_duct_open_port_pole() {
    let mut model = two_source_model();
    // Equal half-wave arms give an open-port pole and a finite loaded response.
    // This exterior frequency stays inside the coarse fixture's physical band.
    for edge in &mut model.driving.network.edges {
        edge.duct.segments = vec![Segment::Cylinder {
            radius: 0.006,
            length: 0.2,
        }];
    }
    let omega = core::f64::consts::PI * model.gas().sound_speed / 0.2;
    assert!(matches!(
        fs_duct::network::network_port_impedance(
            model.driving_point().network(),
            &model.gas(),
            omega,
            LossModel::Lossless
        ),
        Err(DuctError::Singular)
    ));
    let bank = bank(omega);
    let gate = CancelGate::new_clock_free();
    let prepared = model
        .prepare_sources(omega, LossModel::Lossless, &[], Some(&bank), &[], &gate)
        .unwrap();
    let mut coherent = prepared.response_buffer();
    prepared
        .drive_into(prepared.default_excitations(), &mut coherent, &gate)
        .unwrap();
    let direct = direct_fields(&model, omega, coherent.outward_flows_m3_s().unwrap());
    for (actual, expected) in coherent.receiver_fields().unwrap().iter().zip(direct) {
        close(actual.pressure_pa, expected.pressure_pa);
        for i in 0..3 {
            close(
                actual.particle_velocity_m_s[i],
                expected.particle_velocity_m_s[i],
            );
        }
    }
    let mut statistics = prepared.statistics_buffer();
    prepared
        .coherence_into(
            prepared.default_excitations(),
            1,
            NetworkCoherenceBudget::default(),
            FIELD_BUDGET,
            &mut statistics,
            &gate,
        )
        .unwrap();
    for (actual, expected) in statistics
        .receiver_fields()
        .unwrap()
        .iter()
        .zip(coherent.receiver_fields().unwrap())
    {
        real_close(
            actual.mean_square_pressure_pa2,
            0.5 * expected.pressure_pa.norm_sq(),
        );
    }
    let work = statistics.network().unwrap().work().unwrap();
    real_close(
        work.generator_power_w,
        work.source_dissipation_w + work.external_load_power_w,
    );
}
