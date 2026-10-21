//! G1/G3/G4: wall work melts an admitted stationary Darcy mixture.
use super::*;
use fs_lbm::conjugate::enthalpy::advection::feedback::{MaterialFlowError, ViscousHeatingAudit};
use fs_lbm::conjugate::{
    FluidProperties, ViscousHeatControls, creeping_flow_with_darcy_heat,
    creeping_flow_with_viscous_heat,
};

fn stokes() -> SimpleConfig {
    let mut faces = [FvBoundary::wall(); 6];
    faces[3] = FvBoundary::Wall {
        velocity: [1.0, 0.0, 0.0],
    };
    faces[4] = FvBoundary::Symmetry;
    faces[5] = FvBoundary::Symmetry;
    let mut config = SimpleConfig::new(faces);
    config.tolerance = 1e-13;
    config.pressure_tolerance = 1e-12;
    config.momentum_tolerance = 1e-3;
    config
}

fn controls() -> ViscousHeatControls {
    ViscousHeatControls {
        maximum_reynolds: 0.1,
        power_tolerance_w: 1e-7,
    }
}

fn thermal() -> StepConfig {
    // The original harmonic advection CFL remains active. This large proposal
    // lets the explicit short advance endpoint select the accepted interval.
    StepConfig {
        max_dt_s: 10.0,
        energy_tolerance_j: 1e-8,
    }
}

#[test]
fn exact_stokes_darcy_work_uses_the_native_half_cells() {
    let domain = VoxelDomain::new(2, 2, 1, 1.0).unwrap();
    let fluid = FluidProperties {
        density_kg_m3: 1.0,
        specific_heat_j_kg_k: 2.0,
        conductivity_w_m_k: 0.01,
        kinematic_viscosity_m2_s: 100.0,
    };
    let config = stokes();
    let gate = CancelGate::new();
    for beta in [0.0, 200.0, 800.0] {
        let (flow, power) = creeping_flow_with_darcy_heat(
            &domain,
            &fluid,
            &[100.0; 4],
            &[beta; 4],
            &config,
            controls(),
            &gate,
        )
        .unwrap();
        // Independent constrained four-face solution: all nonzero fluxes are
        // +/-a. D=mu(24a²-4Ua+2U²)+4 beta a² and P=2mu U(U-a).
        let a = 100.0 / (1200.0 + 2.0 * beta);
        near(flow.field.flux_x(1, 1, 0), a);
        near(flow.field.flux_x(1, 0, 0), -a);
        near(power.wall_tangential_force_n[3][0], 200.0 * (1.0 - a));
        near(power.supplied_wall_power_w, 200.0 * (1.0 - a));
        near(power.darcy_heat_w, 4.0 * beta * a * a);
        near(
            power.molecular_heat_w,
            100.0 * (24.0 * a * a - 4.0 * a + 2.0),
        );
        near(power.heat_w, power.molecular_heat_w + power.darcy_heat_w);
        near(power.heat_w, power.supplied_wall_power_w);
        for c in 0..4 {
            let molecular = if c < 2 {
                600.0 * a * a
            } else {
                500.0 * a * a + 100.0 * (1.0 - a).powi(2)
            };
            near(power.cell_darcy_heat_w.as_ref().unwrap()[c], beta * a * a);
            near(power.cell_heat_w[c], molecular + beta * a * a);
        }
        near(power.pressure_constraint_power_w, 0.0);
        near(power.residual_w, 0.0);
    }
    let (_, unchanged) =
        creeping_flow_with_viscous_heat(&domain, &fluid, &[100.0; 4], &config, controls(), &gate)
            .unwrap();
    assert!(unchanged.cell_darcy_heat_w.is_none());
    assert_eq!(unchanged.darcy_heat_w, 0.0);
    near(unchanged.heat_w, unchanged.molecular_heat_w);
    assert!(
        creeping_flow_with_darcy_heat(
            &domain,
            &fluid,
            &[100.0; 4],
            &[-1.0; 4],
            &config,
            controls(),
            &gate
        )
        .is_err()
    );
}

#[test]
fn permeability_data_and_latent_heat_drive_the_next_flow_without_temperature_change() {
    let card = card(claims(1.0, 100.0, false));
    let liquid = liquid(&card);
    let curve = phase_curve(&card, [1.0; 2]);
    let domain = VoxelDomain::new(2, 2, 1, 1.0).unwrap();
    let config = stokes();
    let gate = CancelGate::new();
    let ids = [
        Some(CellId(0)),
        Some(CellId(1)),
        Some(CellId(2)),
        Some(CellId(3)),
    ];
    let old = initial(&domain, &curve, 0.0);
    let mut results = Vec::new();
    for permeability in [0.5, 2.0] {
        let mobility = MushyPhaseLaw::try_new(&curve, 1.0, 1.0, permeability, 0.2).unwrap();
        let base = MaterialFlowFeedback::try_new(
            &old,
            &domain,
            input(&ids, &curve),
            &config,
            (),
            &gate,
            20,
            |_, _: &DepositedCell, t| {
                liquid.at_temperature(Temperature::new(t)).map_err(|error| {
                    fs_material::MaterialError::State {
                        what: error.to_string(),
                    }
                })
            },
        )
        .unwrap();
        // Both attachment orders retain the same declared owners.
        assert!(
            base.clone()
                .with_viscous_heating(controls())
                .unwrap()
                .with_mushy_mobility(&mobility)
                .is_ok()
        );
        let observer = base
            .with_mushy_mobility(&mobility)
            .unwrap()
            .with_viscous_heating(controls())
            .unwrap();
        let (warm, first, first_thermal) = old
            .advance_observed(&observer, 0.1, thermal(), || Ok(()))
            .unwrap();
        assert_eq!(first.momentum_solves(), 1);
        let beta0 = 200.0 / permeability;
        let a0 = 100.0 / (1200.0 + 2.0 * beta0);
        let mut h = [50.0; 4];
        for c in 0..4 {
            let shear = if c < 2 {
                600.0 * a0 * a0
            } else {
                500.0 * a0 * a0 + 100.0 * (1.0 - a0).powi(2)
            };
            h[c] += 0.1 * (shear + beta0 * a0 * a0);
            near(warm.cells()[&CellId(c as u64)].specific_h_j_kg, h[c]);
            assert_eq!(warm.temperature_k(CellId(c as u64)).unwrap(), 500.0);
        }
        let first_work = first.viscous_heating_audit().unwrap();
        near(first_work.darcy_heat_j, 0.1 * 4.0 * beta0 * a0 * a0);
        near(first_work.generated_heat_j, 0.1 * 200.0 * (1.0 - a0));
        near(first_thermal.external_heat_j, first_work.generated_heat_j);
        let (next, accepted, thermal_audit) = warm
            .advance_observed(&first, 0.2, thermal(), || Ok(()))
            .unwrap();
        let applied = accepted.applied().unwrap();
        near(applied.time_s, 0.1);
        let beta = h.map(|h| {
            let phi = h / 100.0;
            100.0 * (1.0 - phi).powi(2) / (permeability * phi.powi(3))
        });
        // Heterogeneous beta contributes sum(beta_i)*a² to the four-face
        // dissipation. This reconstructs the next solve independently of its code.
        let a1 = 200.0 / (2400.0 + beta.iter().sum::<f64>());
        assert!(a1 > a0);
        near(applied.flow.field.flux_x(1, 1, 0), a1);
        let before = h;
        let donors = [1, 3, 0, 2];
        for c in 0..4 {
            near(applied.mushy_drag_pa_s_m2.as_ref().unwrap()[c], beta[c]);
            let shear = if c < 2 {
                600.0 * a1 * a1
            } else {
                500.0 * a1 * a1 + 100.0 * (1.0 - a1).powi(2)
            };
            h[c] += 0.1 * (shear + beta[c] * a1 * a1 + a1 * (before[donors[c]] - before[c]));
            near(next.cells()[&CellId(c as u64)].specific_h_j_kg, h[c]);
            assert_eq!(next.temperature_k(CellId(c as u64)).unwrap(), 500.0);
            near(
                mobility.at_specific_enthalpy(h[c]).unwrap().volumes.liquid,
                h[c] / 100.0,
            );
        }
        for sample in accepted.endpoint() {
            assert_eq!(
                sample.material.material().card_identity(),
                card.content_hash()
            );
            assert_eq!(
                sample
                    .material
                    .material()
                    .property(DYNAMIC_VISCOSITY_PROPERTY)
                    .unwrap()
                    .answer()
                    .receipt
                    .selected,
                liquid
                    .reference_state()
                    .material()
                    .property(DYNAMIC_VISCOSITY_PROPERTY)
                    .unwrap()
                    .answer()
                    .receipt
                    .selected
            );
        }
        let work = accepted.viscous_heating_audit().unwrap();
        near(work.duration_s, 0.1);
        near(work.darcy_heat_j, 0.1 * beta.iter().sum::<f64>() * a1 * a1);
        near(
            work.generated_heat_j,
            work.molecular_heat_j + work.darcy_heat_j,
        );
        near(work.supplied_wall_work_j, work.generated_heat_j);
        near(thermal_audit.external_heat_j, work.generated_heat_j);
        near(
            next.energy_j() - old.energy_j(),
            first_work.supplied_wall_work_j + work.supplied_wall_work_j,
        );
        near(accepted.audit().external_enthalpy_j, 0.0);
        near(work.residual_j, 0.0);
        near(next.mass_kg(), old.mass_kg());
        results.push((a1, next.energy_j()));
    }
    assert!(results[1].0 > results[0].0);
    assert!(results[1].1 < results[0].1);
}

#[derive(Clone)]
struct CancelSecond<'a> {
    gate: &'a CancelGate,
    steps: usize,
}
impl DepositionObserver for CancelSecond<'_> {
    type Error = DepositionError;
    fn event(&mut self, _: &Event, _: &DepositionSimulation<'_>) -> Result<(), Self::Error> {
        Ok(())
    }
    fn step(
        &mut self,
        _: f64,
        _: f64,
        _: &BTreeMap<CellId, f64>,
        _: &DepositionSimulation<'_>,
        _: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), Self::Error> {
        self.steps += 1;
        if self.steps == 2 {
            self.gate.request();
        }
        Ok(())
    }
}

#[test]
fn source_phase_refusal_and_late_cancel_publish_no_darcy_heat() {
    let card = card(claims(1.0, 100.0, false));
    let liquid = liquid(&card);
    let curve = phase_curve(&card, [1.0; 2]);
    let mobility = MushyPhaseLaw::try_new(&curve, 1.0, 1.0, 1.0, 0.2).unwrap();
    let domain = VoxelDomain::new(2, 2, 1, 1.0).unwrap();
    let config = stokes();
    let gate = CancelGate::new();
    let ids = [
        Some(CellId(0)),
        Some(CellId(1)),
        Some(CellId(2)),
        Some(CellId(3)),
    ];
    let resolve = |_, _: &DepositedCell, t| {
        liquid.at_temperature(Temperature::new(t)).map_err(|error| {
            fs_material::MaterialError::State {
                what: error.to_string(),
            }
        })
    };
    let old = initial(&domain, &curve, -400.0);
    let observer = MaterialFlowFeedback::try_new(
        &old,
        &domain,
        input(&ids, &curve),
        &config,
        (),
        &gate,
        20,
        resolve,
    )
    .unwrap()
    .with_mushy_mobility(&mobility)
    .unwrap()
    .with_viscous_heating(controls())
    .unwrap();
    // The caloric endpoint still exists, but its cooled liquid fraction drops
    // below the ORIGINAL admitted mobility support, after candidate heating.
    assert!(matches!(
        old.advance_observed(&observer, 0.1, thermal(), || Ok(())),
        Err(MaterialFlowError::Material { .. })
    ));
    assert_eq!(
        observer.viscous_heating_audit(),
        Some(ViscousHeatingAudit::default())
    );
    assert!(observer.applied().is_none());
    let (retry, accepted, audit) = old
        .advance_observed(&observer, 0.01, thermal(), || Ok(()))
        .unwrap();
    let (again, second, _) = old
        .advance_observed(&observer, 0.01, thermal(), || Ok(()))
        .unwrap();
    assert_eq!(retry.cells(), again.cells());
    assert_eq!(
        accepted.viscous_heating_audit(),
        second.viscous_heating_audit()
    );
    near(
        audit.external_heat_j,
        accepted.viscous_heating_audit().unwrap().generated_heat_j - 4.0,
    );
    near(old.time_s(), 0.0);
    near(old.cells()[&CellId(1)].specific_h_j_kg, 50.0);
    let warm = initial(&domain, &curve, 0.0);
    let cancel = MaterialFlowFeedback::try_new(
        &warm,
        &domain,
        input(&ids, &curve),
        &config,
        CancelSecond {
            gate: &gate,
            steps: 0,
        },
        &gate,
        20,
        resolve,
    )
    .unwrap()
    .with_mushy_mobility(&mobility)
    .unwrap()
    .with_viscous_heating(controls())
    .unwrap();
    assert!(matches!(
        warm.advance_observed(
            &cancel,
            0.15,
            StepConfig {
                max_dt_s: 0.06,
                ..thermal()
            },
            || Ok(())
        ),
        Err(MaterialFlowError::Thermal(DepositionError::Cancelled))
    ));
    assert_eq!(cancel.inner().steps, 0);
    assert_eq!(
        cancel.viscous_heating_audit(),
        Some(ViscousHeatingAudit::default())
    );
    near(warm.energy_j(), 200.0);
    near(warm.time_s(), 0.0);
}
