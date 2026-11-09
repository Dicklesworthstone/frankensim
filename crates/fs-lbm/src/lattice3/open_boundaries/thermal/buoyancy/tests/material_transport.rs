use super::*;
use crate::lattice3::ViscousLattice3;
use crate::material::plan_scaling_from_material;
use fs_material::{
    MaterialError,
    fluid_transport::{FluidTransport, FluidTransportState},
};
use std::sync::{
    Arc,
    atomic::{AtomicUsize, Ordering},
};

const PHYSICAL_RHO: f64 = 1000.0;
const PHYSICAL_CP: f64 = 4000.0;
const DX: f64 = 0.001;
const DT: f64 = 0.05;

fn material_config() -> BoussinesqOpenMaterialConfig3 {
    BoussinesqOpenMaterialConfig3 {
        temperature_at_zero_k: 300.0,
        temperature_scale_k: 20.0,
        maximum_relative_density_departure: 0.01,
        maximum_relative_specific_heat_departure: 0.02,
        maximum_workspace_bytes: 1024 * 1024,
        maximum_report_bytes: 1024 * 1024,
    }
}

fn material_at(temperature_k: f64) -> FluidTransportState {
    // Independent, explicit test constitutive law. Variable density/Cp must
    // be admitted but must not replace the native reference storage or force.
    let delta = temperature_k - 310.0;
    FluidTransportState::try_new(
        PHYSICAL_RHO * (1.0 - 0.0001 * delta),
        PHYSICAL_CP * (1.0 + 0.0002 * delta),
        2.0 * (1.0 + 0.006 * delta),
        0.002 * (1.0 - 0.004 * delta),
    )
    .unwrap()
}

#[derive(Debug)]
struct TracedMaterial {
    cell: usize,
    temperature_k: f64,
    ordinal: usize,
    state: FluidTransportState,
    reads: Arc<AtomicUsize>,
}
impl FluidTransport for TracedMaterial {
    fn transport_state(&self) -> Result<FluidTransportState, MaterialError> {
        self.reads.fetch_add(1, Ordering::SeqCst);
        Ok(self.state)
    }
}

fn check_explicit_material_fields<
    L: EquilibriumLattice3<Q> + ViscousLattice3<Q>,
    const Q: usize,
>(
    collision: L::Collision,
    temperature_scale_k: f64,
) {
    let mut actual = periodic::<L, Q>(collision, 0.5);
    let mut explicit = periodic::<L, Q>(collision, 0.5);
    let reference = material_at(310.0);
    let scaling = plan_scaling_from_material(&reference, DX, DT, 0.0).unwrap();
    let config = BoussinesqOpenMaterialConfig3 {
        temperature_scale_k,
        ..material_config()
    };
    let mut workspace = actual
        .material_transport_workspace(&scaling, config)
        .unwrap();
    let cells = actual.thermal().populations().len();
    let heat: Vec<_> = (0..cells)
        .map(|cell| 0.001 * (1 + cell % 4) as f64)
        .collect();
    let reads = Arc::new(AtomicUsize::new(0));
    let mut largest_tau_spread = 0.0_f64;
    for _ in 0..4 {
        let mut collisions = Vec::new();
        let mut thermal = Vec::new();
        let mut opening_k = Vec::new();
        for cell in 0..cells {
            let kelvin =
                300.0 + temperature_scale_k * explicit.thermal().temperature(cell).unwrap();
            let source = material_at(kelvin);
            let units = (DT / DX) / DX;
            // Independent physical conversion, using reference rho and Cp
            // rather than the variable values returned by the source.
            let momentum_tau =
                0.5 + (source.dynamic_viscosity_pa_s() / PHYSICAL_RHO) * units / (1.0 / 3.0);
            let thermal_tau = 0.5
                + ((source.thermal_conductivity_w_m_k() / PHYSICAL_RHO) / PHYSICAL_CP) * units
                    / 0.25;
            collisions.push(L::with_shear_relaxation(collision, momentum_tau).unwrap());
            thermal.push(thermal_tau);
            opening_k.push(kelvin);
        }
        let expected = explicit
            .step_with_poll(
                Some(&collisions),
                &heat,
                Some(&thermal),
                [0.0; 3],
                &mut || Ok(()),
            )
            .unwrap();
        let mut ordinal = 0;
        reads.store(0, Ordering::SeqCst);
        let report = actual
            .advance_material_with_poll(
                None,
                None,
                &heat,
                &mut workspace,
                [0.0; 3],
                &mut |cell, temperature_k, retained| {
                    assert_eq!(*retained, reference);
                    let source = TracedMaterial {
                        cell,
                        temperature_k,
                        ordinal,
                        state: material_at(temperature_k),
                        reads: Arc::clone(&reads),
                    };
                    ordinal += 1;
                    Ok(source)
                },
                &mut || Ok(()),
            )
            .unwrap();
        assert_eq!(ordinal, 2 * cells);
        assert_eq!(reads.load(Ordering::SeqCst), 2 * cells);
        assert_eq!(report.transport, expected);
        assert_eq!(accepted(&actual), accepted(&explicit));
        assert_eq!(workspace.thermal_relaxation(), thermal);
        for (epoch, samples) in [&report.applied, &report.endpoint].into_iter().enumerate() {
            assert_eq!(samples.len(), cells);
            for sample in samples {
                let source = sample.material.material();
                assert_eq!(source.cell, sample.cell);
                assert_eq!(source.ordinal, epoch * cells + sample.cell);
                assert_eq!(
                    source.temperature_k.to_bits(),
                    sample.temperature_k.to_bits()
                );
                assert_eq!(sample.material.resolved_transport(), &source.state);
                assert_eq!(sample.transport.density_kg_m3(), PHYSICAL_RHO);
                assert_eq!(sample.transport.specific_heat_j_kg_k(), PHYSICAL_CP);
                let expected_k = if epoch == 0 {
                    opening_k[sample.cell]
                } else {
                    300.0 + temperature_scale_k * actual.thermal().temperature(sample.cell).unwrap()
                };
                assert_eq!(sample.temperature_k.to_bits(), expected_k.to_bits());
            }
        }
        let range = workspace
            .thermal_relaxation()
            .iter()
            .fold([f64::INFINITY, f64::NEG_INFINITY], |[lo, hi], &tau| {
                [lo.min(tau), hi.max(tau)]
            });
        largest_tau_spread = largest_tau_spread.max(range[1] - range[0]);
        near(
            report.transport.heat.stored_change,
            heat.iter().sum(),
            2e-13,
        );
        near(report.transport.heat.combined_defect, 0.0, 2e-13);
    }
    assert!(largest_tau_spread > 1e-6, "{largest_tau_spread:e}");
}

#[test]
fn material_queries_and_local_rates_drive_both_native_lattices_without_rebinding_storage() {
    check_explicit_material_fields::<D3q19, 19>(CollisionModel3::Bgk { tau: 0.9 }, 20.0);
    check_explicit_material_fields::<D3q27, 27>(
        CollisionConfig27 {
            model: CollisionModel27::Trt {
                tau: 0.9,
                magic: 3.0 / 16.0,
                equilibrium: Equilibrium27::Factorized,
            },
            maximum_mach: 0.2,
        },
        -20.0,
    );
}

#[test]
fn endpoint_material_refusal_and_final_cancellation_preserve_background_walls_and_exact_retry() {
    let mut actual = channel(BETA);
    let mut reference = channel(BETA);
    let heat = vec![0.0004; 72];
    for coupled in [&mut actual, &mut reference] {
        advance(coupled, &heat);
        coupled.set_finite_wall_heat(2, &vec![0.0003; 18]).unwrap();
    }
    let physical_reference = material_at(310.0);
    let scaling = plan_scaling_from_material(&physical_reference, DX, DT, 0.0).unwrap();
    assert!(
        actual
            .material_transport_workspace(
                &scaling,
                BoussinesqOpenMaterialConfig3 {
                    temperature_at_zero_k: 0.0,
                    ..material_config()
                }
            )
            .is_err()
    );
    let mut workspace = actual
        .material_transport_workspace(&scaling, material_config())
        .unwrap();
    let mut reference_workspace = reference
        .material_transport_workspace(&scaling, material_config())
        .unwrap();
    let before = accepted(&actual);
    let energy_before = total_energy(&reference);
    let opening_k: Vec<_> = (0..72)
        .map(|cell| 300.0 + 20.0 * actual.thermal().temperature(cell).unwrap())
        .collect();
    let mut background: Vec<_> = (0..72)
        .map(|cell| {
            [
                BACKGROUND[0] + 1e-6 * (1 + cell % 3) as f64,
                BACKGROUND[1] - 1e-6,
                BACKGROUND[2],
            ]
        })
        .collect();
    let original_endpoint = background.clone();
    let mut reference_background = background.clone();
    let mut polls = 0;
    let expected = reference
        .advance_material_with_poll(
            Some(&mut reference_background),
            None,
            &heat,
            &mut reference_workspace,
            [0.0; 3],
            &mut |_, temperature, _| Ok(material_at(temperature)),
            &mut || {
                polls += 1;
                Ok(())
            },
        )
        .unwrap();
    let changed = expected
        .endpoint
        .iter()
        .rev()
        .find(|sample| (sample.temperature_k - opening_k[sample.cell]).abs() > 1e-6)
        .expect("the sourced, conducting wall changes an actual fluid temperature");
    let refusal_cell = changed.cell;
    let admitted_change = 0.5 * (changed.temperature_k - opening_k[refusal_cell]).abs();
    let refused = actual.advance_material_with_poll(
        Some(&mut background),
        None,
        &heat,
        &mut workspace,
        [0.0; 3],
        &mut |cell, temperature, _| {
            if cell == refusal_cell && (temperature - opening_k[cell]).abs() > admitted_change {
                Err(MaterialError::State {
                    what: "outside this source's temperature support".into(),
                })
            } else {
                Ok(material_at(temperature))
            }
        },
        &mut || Ok(()),
    );
    assert!(
        matches!(refused, Err(BoussinesqOpenMaterialError3::Material(
        crate::material::transport::MaterialTransportError::State { cell, .. }
    )) if cell == refusal_cell)
    );
    assert_eq!(background, original_endpoint);
    assert_eq!(accepted(&actual), before);

    let mut calls = 0;
    let refused = actual.advance_material_with_poll(
        Some(&mut background),
        None,
        &heat,
        &mut workspace,
        [0.0; 3],
        &mut |_, temperature, _| Ok(material_at(temperature)),
        &mut || {
            calls += 1;
            if calls == polls {
                Err(GridError3::Cancelled)
            } else {
                Ok(())
            }
        },
    );
    assert!(matches!(
        refused,
        Err(BoussinesqOpenMaterialError3::Transport(
            BoussinesqOpenError3::Flow(GridError3::Cancelled)
        ))
    ));
    assert_eq!(calls, polls);
    assert_eq!(background, original_endpoint);
    assert_eq!(accepted(&actual), before);
    let retried = actual
        .advance_material_with_poll(
            Some(&mut background),
            None,
            &heat,
            &mut workspace,
            [0.0; 3],
            &mut |_, temperature, _| Ok(material_at(temperature)),
            &mut || Ok(()),
        )
        .unwrap();
    assert_eq!(retried.transport, expected.transport);
    assert_eq!(retried.applied, expected.applied);
    assert_eq!(retried.endpoint, expected.endpoint);
    assert_eq!(background, reference_background);
    assert_eq!(accepted(&actual), accepted(&reference));
    let report = retried.transport.heat;
    let wall_source: f64 = report
        .walls
        .into_iter()
        .flatten()
        .filter_map(|wall| wall.solid)
        .map(|solid| solid.source)
        .sum();
    near(
        total_energy(&actual) - energy_before,
        heat.iter().sum::<f64>()
            + wall_source
            + report
                .reservoirs
                .iter()
                .map(|reservoir| reservoir.heat_in)
                .sum::<f64>(),
        2e-11,
    );
    for cell in 0..72 {
        let theta = actual.thermal().temperature(cell).unwrap();
        for axis in 0..3 {
            near(
                actual.flow().forces()[cell][axis],
                original_endpoint[cell][axis]
                    - REFERENCE_RHO * BETA * (theta - 0.5) * GRAVITY[axis],
                2e-15,
            );
        }
    }
}

mod sgs;
