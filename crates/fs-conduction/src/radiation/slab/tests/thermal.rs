use super::*;
use fs_thermal::enthalpy::backward_euler::BackwardEulerStep;
use fs_thermal::enthalpy::{
    Boundary, CellGeometry, CellId, DepositedCell, DepositionBudget, DepositionError,
    DepositionObserver, DepositionSimulation, Event, MaterialId, StepConfig, ThermalLaw,
    ThermalLink, TimedEvent,
};

const SIGMA: f64 = crate::radiation::STEFAN_BOLTZMANN_W_M2_K4;

mod cells;
mod spectral;

struct Caloric {
    capacity: f64,
    datum: f64,
    latent: f64,
}
impl Caloric {
    fn energy(&self, temperature: f64) -> f64 {
        self.datum
            + self.capacity * temperature
            + if temperature > 500.0 {
                self.latent
            } else {
                0.0
            }
    }
}
impl ThermalLaw for Caloric {
    fn temperature_k(&self, energy: f64) -> Result<f64, DepositionError> {
        let e = energy - self.datum;
        let melt = self.capacity * 500.0;
        let t = if e < melt {
            e / self.capacity
        } else if e <= melt + self.latent {
            500.0
        } else {
            (e - self.latent) / self.capacity
        };
        if t.is_finite() && (100.0..=10_000.0).contains(&t) {
            Ok(t)
        } else {
            Err(DepositionError::MaterialDomain)
        }
    }
    fn maximum_temperature_slope(&self) -> f64 {
        1.0 / self.capacity
    }
}

fn exitance(t: f64) -> f64 {
    (SIGMA * (t * t)) * (t * t)
}
fn config() -> StepConfig {
    StepConfig {
        max_dt_s: 1.0,
        energy_tolerance_j: 1e-7,
    }
}
fn bindings(law: &Caloric, n: usize) -> Vec<(CellId, &dyn ThermalLaw)> {
    (0..n)
        .map(|i| (CellId(i as u64), law as &dyn ThermalLaw))
        .collect()
}
fn state<'a>(
    law: &'a Caloric,
    slab: &AbsorbingSlab,
    energies: &[f64],
    area: f64,
    contact: f64,
    max_links: usize,
) -> DepositionSimulation<'a> {
    let mut schedule: Vec<_> = energies
        .iter()
        .enumerate()
        .map(|(i, &h)| TimedEvent {
            time_s: 0.0,
            event: Event::Birth(DepositedCell {
                id: CellId(i as u64),
                segment_id: i as u64 + 20,
                body_id: 19,
                material_id: MaterialId(7),
                mass_kg: 1.0,
                specific_h_j_kg: h,
                geometry: CellGeometry {
                    start_m: [i as f64, 0.0, 0.0],
                    end_m: [i as f64 + 1.0, 0.0, 0.0],
                    volume_m3: area * slab.layers()[i].thickness_m(),
                    surface_area_m2: 6.0,
                },
                boundary: Boundary {
                    ambient_k: 300.0,
                    htc_w_m2_k: 0.0,
                    reservoir: None,
                    source_w: 0.0,
                },
            }),
        })
        .collect();
    if contact > 0.0 {
        schedule.push(TimedEvent {
            time_s: 0.0,
            event: Event::Contact(ThermalLink {
                id: 31,
                a: CellId(0),
                b: CellId(1),
                area_m2: area,
                conductance_w_k: contact,
            }),
        });
    }
    DepositionSimulation::new(
        schedule,
        &[(MaterialId(7), law as &dyn ThermalLaw)],
        DepositionBudget {
            max_cells: 8,
            max_links,
            max_events: 10,
            max_substeps: 100_000,
        },
    )
    .unwrap()
    .advance_to(0.0, config(), || Ok(()))
    .unwrap()
    .0
}

#[test]
fn g1_gray_irradiation_heats_original_latent_storage_without_inventing_temperature() {
    let slab = slab(&[1000.0]);
    let law = Caloric {
        capacity: 1000.0,
        datum: 1e9,
        latent: 20_000.0,
    };
    let old = state(&law, &slab, &[law.energy(500.0) + 10_000.0], 0.01, 0.0, 0);
    let source = GraySlabThermal::new_full_spectrum(
        &slab,
        0.01,
        [2.0 * exitance(500.0); 2],
        &bindings(&law, 1),
        &old,
        || Ok(()),
    )
    .unwrap();
    let (next, _, audit) = old
        .advance_observed(&source, 0.05, config(), || Ok(()))
        .unwrap();
    let expected = 0.05 * 0.02 * exitance(500.0);
    close(audit.external_heat_j, expected, 1e-12);
    close(audit.stored_energy_change_j, expected, 1e-7);
    assert_eq!(next.temperature_k(CellId(0)).unwrap(), 500.0);
    assert!(next.cells()[&CellId(0)].specific_h_j_kg > old.cells()[&CellId(0)].specific_h_j_kg);
    assert_eq!(
        next.cells()[&CellId(0)].geometry,
        old.cells()[&CellId(0)].geometry
    );
    assert_eq!(
        next.cells()[&CellId(0)].boundary,
        old.cells()[&CellId(0)].boundary
    );
    assert_eq!(source.slab().layers().as_ptr(), slab.layers().as_ptr());
    assert_eq!(old.time_s(), 0.0);
}

#[test]
fn g3_gray_trapped_exchange_joins_native_conduction_without_external_heat_or_coverage() {
    let slab = refractive(&[1000.0, 1000.0], 1.5);
    let law = Caloric {
        capacity: 1000.0,
        datum: 0.0,
        latent: 0.0,
    };
    let old = state(
        &law,
        &slab,
        &[law.energy(900.0), law.energy(300.0)],
        0.01,
        0.5,
        1,
    );
    let mut source = GraySlabThermal::new_full_spectrum(
        &slab,
        0.01,
        [exitance(900.0), exitance(300.0)],
        &bindings(&law, 2),
        &old,
        || Ok(()),
    )
    .unwrap();
    let transport = source.endpoint_transport(&old, || Ok(())).unwrap();
    assert!(transport.heat.bulk.sources_w.is_empty());
    assert_eq!(transport.transfers.len(), 1);
    let power = 0.01 * 2.25 * (exitance(900.0) - exitance(300.0));
    close(transport.transfers[0].power_w, power, 1e-11);
    // Native conduction remains on its original contact; radiation consumes no
    // extra patch and contributes no external energy even above the n=1 limit.
    let dt = 1e-4;
    let (next, _, audit) = old
        .advance_observed(&source, dt, config(), || Ok(()))
        .unwrap();
    let energy = dt * (power + 0.5 * 600.0);
    close(
        next.cells()[&CellId(0)].specific_h_j_kg,
        law.energy(900.0) - energy,
        2e-10,
    );
    close(
        next.cells()[&CellId(1)].specific_h_j_kg,
        law.energy(300.0) + energy,
        2e-10,
    );
    assert_eq!(audit.external_heat_j, 0.0);
    assert_eq!(next.links(), old.links());
    let proposed = 0.25;
    let limited = source
        .limit_timestep_s(&old, proposed, &mut || Ok(()))
        .unwrap();
    let row = slab.face_coupling()[0].iter().sum::<f64>() + 2.25;
    let decay = 0.01 * row * 4.0 * SIGMA * 3000.0_f64.powi(3) / 1000.0;
    close(limited, proposed / (1.0 + proposed * decay / 0.8), 1e-15);
    assert!(limited < proposed);
}

#[test]
fn g0_gray_equilibrium_and_nearby_temperatures_retain_detailed_balance() {
    let slab = refractive(&[0.2, 0.0, 2.0], 1.7);
    let law = Caloric {
        capacity: 1.0,
        datum: 0.0,
        latent: 0.0,
    };
    let old = state(&law, &slab, &[500.0; 3], 0.001, 0.0, 3);
    let source = GraySlabThermal::new_full_spectrum(
        &slab,
        0.001,
        [exitance(500.0); 2],
        &bindings(&law, 3),
        &old,
        || Ok(()),
    )
    .unwrap();
    let transport = source.endpoint_transport(&old, || Ok(())).unwrap();
    assert!(transport.heat.bulk.sources_w.is_empty());
    assert!(transport.transfers.is_empty());
    let (next, _, audit) = old
        .advance_observed(&source, 0.2, config(), || Ok(()))
        .unwrap();
    assert_eq!(next.cells(), old.cells());
    assert_eq!(audit.external_heat_j, 0.0);
    let slab = refractive(&[1000.0, 1000.0], 1.5);
    let b = f64::from_bits(500.0_f64.to_bits() - 1);
    let old = state(&law, &slab, &[500.0, b], 0.01, 0.0, 1);
    let source = GraySlabThermal::new_full_spectrum(
        &slab,
        0.01,
        [exitance(500.0), exitance(b)],
        &bindings(&law, 2),
        &old,
        || Ok(()),
    )
    .unwrap();
    let actual = source
        .endpoint_transport(&old, || Ok(()))
        .unwrap()
        .transfers[0]
        .power_w;
    let first_order = 0.01 * 2.25 * 4.0 * SIGMA * 500.0_f64.powi(3) * (500.0 - b);
    assert!(actual > 0.0);
    assert!((actual / first_order - 1.0).abs() < 8.0 * f64::EPSILON);
}

#[test]
fn g1_gray_endpoint_transport_drives_original_backward_euler_energy_residual() {
    let slab = refractive(&[1000.0, 1000.0], 1.5);
    let law = Caloric {
        capacity: 1000.0,
        datum: 0.0,
        latent: 0.0,
    };
    let [a, b] = [550.0, 350.0];
    let irradiation = [exitance(a), exitance(b)];
    let energy = 0.01 * 2.25 * (irradiation[0] - irradiation[1]) + 0.5 * (a - b);
    // Independently author the initial energy for this known dt=1 endpoint.
    // The original residual must combine nonlinear radiation with native G.
    let old = state(
        &law,
        &slab,
        &[law.energy(a) + energy, law.energy(b) - energy],
        0.01,
        0.5,
        1,
    );
    let source = GraySlabThermal::new_full_spectrum(
        &slab,
        0.01,
        irradiation,
        &bindings(&law, 2),
        &old,
        || Ok(()),
    )
    .unwrap();
    let step = BackwardEulerStep::new(&old, 1.0, || Ok(())).unwrap();
    let trial = step
        .evaluate_with_transport(
            &[-energy, energy],
            |candidate| source.endpoint_transport(candidate, || Ok(())),
            || Ok(()),
        )
        .unwrap();
    assert!(trial.maximum_cell_residual_j() < 1e-8);
    let (next, (), audit) = trial.accept(1e-8, 1e-8, || Ok(())).unwrap();
    close(next.temperature_k(CellId(0)).unwrap(), a, 1e-12);
    close(next.temperature_k(CellId(1)).unwrap(), b, 1e-12);
    assert_eq!(audit.external_heat_j, 0.0);
    assert_eq!(next.links(), old.links());
    assert_eq!(old.time_s(), 0.0);
}

#[test]
fn g1_gray_immersed_thermal_equilibrium_heat_and_stiffness_use_exterior_n_squared() {
    let slab = AbsorbingSlab::new_in_medium(
        vec![
            AbsorbingLayer::new(1.0, 1000.0, [100.0, 3000.0], "full-spectrum absorption")
                .unwrap()
                .with_refractive_index(2.0, "gray material index")
                .unwrap(),
        ],
        2.0,
        "gray exterior index",
        SlabQuadrature::default(),
        || false,
    )
    .unwrap();
    let law = Caloric {
        capacity: 1000.0,
        datum: 0.0,
        latent: 0.0,
    };
    let old = state(&law, &slab, &[law.energy(500.0)], 0.01, 0.0, 0);
    let equilibrium = GraySlabThermal::new_full_spectrum(
        &slab,
        0.01,
        [4.0 * exitance(500.0); 2],
        &bindings(&law, 1),
        &old,
        || Ok(()),
    )
    .unwrap();
    assert!(
        equilibrium
            .endpoint_transport(&old, || Ok(()))
            .unwrap()
            .heat
            .bulk
            .sources_w
            .is_empty()
    );
    let (unchanged, _, audit) = old
        .advance_observed(&equilibrium, 0.1, config(), || Ok(()))
        .unwrap();
    assert_eq!(unchanged.cells(), old.cells());
    assert_eq!(audit.external_heat_j, 0.0);

    let mut illuminated = GraySlabThermal::new_full_spectrum(
        &slab,
        0.01,
        [4.0 * exitance(600.0), 4.0 * exitance(500.0)],
        &bindings(&law, 1),
        &old,
        || Ok(()),
    )
    .unwrap();
    // Matched-index opaque slab: each face has absorptance one, while its
    // reciprocal escape is n_exterior² σT⁴. No Fresnel/quadrature oracle needed.
    let expected_power = 0.01 * 4.0 * (exitance(600.0) - exitance(500.0));
    let transport = illuminated.endpoint_transport(&old, || Ok(())).unwrap();
    close(
        transport.heat.bulk.sources_w[&CellId(0)],
        expected_power,
        1e-12,
    );
    assert!(transport.transfers.is_empty());
    let dt = 1e-4;
    let (heated, _, audit) = old
        .advance_observed(&illuminated, dt, config(), || Ok(()))
        .unwrap();
    close(audit.external_heat_j, dt * expected_power, 1e-13);
    close(
        heated.cells()[&CellId(0)].specific_h_j_kg - old.cells()[&CellId(0)].specific_h_j_kg,
        dt * expected_power,
        1e-10,
    );
    let proposed = 0.25;
    let limited = illuminated
        .limit_timestep_s(&old, proposed, &mut || Ok(()))
        .unwrap();
    let decay = 0.01 * 2.0 * 4.0 * 4.0 * SIGMA * 3000.0_f64.powi(3) / 1000.0;
    close(limited, proposed / (1.0 + proposed * decay / 0.8), 1e-15);
    assert_eq!(
        illuminated.slab().exterior_index_source(),
        Some("gray exterior index")
    );
}

#[test]
fn g2_gray_native_binding_support_budget_and_cancellation_refuse_transactionally() {
    let slab = slab(&[1.0, 1.0]);
    let law = Caloric {
        capacity: 1000.0,
        datum: 0.0,
        latent: 0.0,
    };
    let copy = Caloric {
        capacity: 1000.0,
        datum: 0.0,
        latent: 0.0,
    };
    let old = state(&law, &slab, &[300_000.0, 400_000.0], 0.01, 0.0, 1);
    let make = |cells: &[(CellId, &dyn ThermalLaw)], area, thermal: &DepositionSimulation<'_>| {
        GraySlabThermal::new_full_spectrum(&slab, area, [0.0; 2], cells, thermal, || Ok(()))
            .map(|_| ())
    };
    assert!(make(&bindings(&copy, 2), 0.01, &old).is_err());
    assert!(make(&[(CellId(0), &law), (CellId(0), &law)], 0.01, &old).is_err());
    assert!(make(&bindings(&law, 2), 0.02, &old).is_err());
    let tiny_budget = state(&law, &slab, &[300_000.0, 400_000.0], 0.01, 0.0, 0);
    assert!(matches!(
        make(&bindings(&law, 2), 0.01, &tiny_budget),
        Err(DepositionError::Budget(_))
    ));
    let unsupported = state(&law, &slab, &[3_500_000.0, 400_000.0], 0.01, 0.0, 1);
    assert_eq!(
        make(&bindings(&law, 2), 0.01, &unsupported),
        Err(DepositionError::MaterialDomain)
    );
    let source =
        GraySlabThermal::new_full_spectrum(&slab, 0.01, [0.0; 2], &bindings(&law, 2), &old, || {
            Ok(())
        })
        .unwrap();
    let changed = state(&law, &slab, &[300_000.0, 400_000.0], 0.02, 0.0, 1);
    assert!(source.endpoint_transport(&changed, || Ok(())).is_err());
    let saved = old.cells().clone();
    let mut seen = 0;
    assert!(matches!(
        old.advance_observed(&source, 0.1, config(), || {
            seen += 1;
            if seen == 8 {
                Err(DepositionError::Cancelled)
            } else {
                Ok(())
            }
        }),
        Err(DepositionError::Cancelled)
    ));
    assert_eq!(old.cells(), &saved);
    assert_eq!(old.time_s(), 0.0);
    assert!(source.endpoint_transport(&old, || Ok(())).is_ok());
}

#[test]
fn g1_gray_asymmetric_ports_keep_each_n_squared_in_heat_and_stiffness() {
    let slab = super::asymmetric::between(&[1000.0], &[1.5], [1.0, 1.7]);
    let law = Caloric {
        capacity: 1000.0,
        datum: 0.0,
        latent: 0.0,
    };
    let area = 0.01;
    let old = state(&law, &slab, &[law.energy(500.0)], area, 0.0, 0);
    let n_squared = [1.0, 1.7 * 1.7];
    let bath = n_squared.map(|n2| n2 * exitance(500.0));
    let equilibrium =
        GraySlabThermal::new_full_spectrum(&slab, area, bath, &bindings(&law, 1), &old, || Ok(()))
            .unwrap();
    assert!(
        equilibrium
            .endpoint_transport(&old, || Ok(()))
            .unwrap()
            .heat
            .bulk
            .sources_w
            .is_empty()
    );
    let h = [bath[0] + 100.0, bath[1] + 300.0];
    let mut source =
        GraySlabThermal::new_full_spectrum(&slab, area, h, &bindings(&law, 1), &old, || Ok(()))
            .unwrap();
    let face = slab.face_coupling()[0];
    let power = area * (face[0] * 100.0 + face[1] * 300.0);
    let transport = source.endpoint_transport(&old, || Ok(())).unwrap();
    close(transport.heat.bulk.sources_w[&CellId(0)], power, 2e-12);
    let (_, _, audit) = old
        .advance_observed(&source, 1e-4, config(), || Ok(()))
        .unwrap();
    close(audit.external_heat_j, 1e-4 * power, 2e-14);
    let proposed = 0.25;
    let limited = source
        .limit_timestep_s(&old, proposed, &mut || Ok(()))
        .unwrap();
    let coupling = area * (face[0] * n_squared[0] + face[1] * n_squared[1]);
    let decay = coupling * 4.0 * SIGMA * 3000.0_f64.powi(3) / 1000.0;
    close(limited, proposed / (1.0 + proposed * decay / 0.8), 2e-15);
}
