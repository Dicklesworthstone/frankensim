//! Phase occupancy is the native fixed-wall topology, while h cells persist.
use super::continuity::{chart, resolve, source};
use super::*;
use fs_lbm::conjugate::ViscousHeatControls;
use fs_lbm::conjugate::enthalpy::advection::feedback::MaterialFlowError;

#[path = "occupancy/unsteady.rs"]
mod unsteady;

#[path = "occupancy/environment.rs"]
mod environment;

#[path = "occupancy/locking.rs"]
mod locking;

fn geometry<'a>(
    domain: &VoxelDomain,
    curve: &'a EquilibriumEnthalpyPhaseCurve,
    enthalpy: [f64; 9],
) -> DepositionSimulation<'a> {
    let events = (0..9)
        .map(|voxel| {
            let [x, y, z] = domain.coords(voxel);
            let mut start = domain.center(x, y, z);
            let mut end = start;
            start[0] -= 0.5;
            end[0] += 0.5;
            TimedEvent {
                time_s: 0.0,
                event: Event::Birth(DepositedCell {
                    id: CellId(voxel as u64),
                    segment_id: voxel as u64,
                    body_id: 1,
                    material_id: MAT,
                    mass_kg: 1.0,
                    specific_h_j_kg: enthalpy[voxel],
                    geometry: CellGeometry {
                        start_m: start,
                        end_m: end,
                        volume_m3: 1.0,
                        surface_area_m2: 6.0,
                    },
                    boundary: Boundary {
                        ambient_k: 500.0,
                        htc_w_m2_k: 0.0,
                        reservoir: None,
                        source_w: 0.0,
                    },
                }),
            }
        })
        .collect();
    DepositionSimulation::new(
        events,
        &[(MAT, curve as &dyn ThermalLaw)],
        DepositionBudget {
            max_cells: 9,
            max_links: 32,
            max_events: 9,
            max_substeps: 20,
        },
    )
    .unwrap()
    .advance_to(0.0, step(0.1), || Ok(()))
    .unwrap()
    .0
}

fn closed() -> SimpleConfig {
    let mut faces = [FvBoundary::wall(); 6];
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
        power_tolerance_w: 1e-8,
    }
}
fn step(max_dt_s: f64) -> StepConfig {
    StepConfig {
        max_dt_s,
        energy_tolerance_j: 1e-8,
    }
}
fn input<'a>(
    ids: &'a [Option<CellId>],
    curve: &'a EquilibriumEnthalpyPhaseCurve,
) -> FlowEnthalpyInput<'a> {
    FlowEnthalpyInput {
        cell_ids: ids,
        material_id: MAT,
        caloric: curve,
        inlet_specific_h_j_kg: [None; 6],
        max_cells: 9,
        max_connections: 64,
    }
}

#[derive(Clone)]
struct CellPower {
    cell: CellId,
    watts: f64,
}
impl DepositionObserver for CellPower {
    type Error = DepositionError;
    fn sources(
        &mut self,
        _: &DepositionSimulation<'_>,
        _: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<BTreeMap<CellId, f64>, Self::Error> {
        Ok(BTreeMap::from([(self.cell, self.watts)]))
    }
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
        Ok(())
    }
}

fn corner_force() -> Vec<[f64; 3]> {
    let mut force = vec![[0.0; 3]; 9];
    force[0] = [-800.0, 800.0, 0.0];
    force
}

#[test]
fn a_mushy_ring_uses_exact_solid_faces_then_reopens_them_when_the_same_cell_melts() {
    let domain = VoxelDomain::new(3, 3, 1, 1.0).unwrap();
    let ids = std::array::from_fn::<_, 9, _>(|c| Some(CellId(c as u64)));
    let gate = CancelGate::new();
    let config = closed();
    for mu in [100.0, 200.0] {
        let card = source(mu);
        let curve = chart(&card);
        let law = MushyPhaseLaw::try_new_with_immobile_solid(&curve, 1.0, 1.0, 1.0).unwrap();
        let mut h = [50.0; 9];
        h[4] = -20.0;
        let old = geometry(&domain, &curve, h);
        let observer = MaterialFlowFeedback::try_new_phase(
            &old,
            &domain,
            input(&ids, &curve),
            &config,
            CellPower {
                cell: CellId(4),
                watts: 2000.0,
            },
            &gate,
            10,
            &law,
            |_, cell: &DepositedCell, _| resolve(&card, &law, cell.specific_h_j_kg),
        )
        .unwrap()
        .with_body_force_heating(corner_force(), controls())
        .unwrap();

        let (melted, accepted, audit) = old
            .advance_observed(&observer, 0.1, step(0.1), || Ok(()))
            .unwrap();
        assert_eq!(accepted.momentum_solves(), 1);
        let epoch = accepted.applied().unwrap();
        assert_eq!(epoch.samples.len(), 9);
        assert!(epoch.samples[4].transport.mobile_transport().is_none());
        assert_eq!(epoch.flow.velocity_m_s[4], [0.0; 3]);
        for x in 1..3 {
            assert_eq!(epoch.flow.field.flux_x(x, 1, 0), 0.0);
        }
        for y in 1..3 {
            assert_eq!(epoch.flow.field.flux_y(1, y, 0), 0.0);
        }

        // Independent one-coordinate balance around the immobile center.
        // Eight fluid cells form one cycle q. Full normal stress + outer and
        // inner no-slip walls give 48 mu q². At phi=1/2, beta=2mu and eight MAC
        // dual volumes give 8 beta q². Two half-cell forces supply 800 q.
        let beta = 2.0 * mu;
        let q = 800.0 / (48.0 * mu + 8.0 * beta);
        for x in 1..3 {
            near(epoch.flow.field.flux_x(x, 0, 0), -q);
            near(epoch.flow.field.flux_x(x, 2, 0), q);
        }
        for y in 1..3 {
            near(epoch.flow.field.flux_y(0, y, 0), q);
            near(epoch.flow.field.flux_y(2, y, 0), -q);
        }
        let power = epoch.viscous_power.as_ref().unwrap();
        near(power.molecular_heat_w, 48.0 * mu * q * q);
        near(power.darcy_heat_w, 8.0 * beta * q * q);
        near(power.supplied_body_power_w, 800.0 * q);
        near(power.heat_w, power.supplied_body_power_w);
        assert_eq!(power.cell_heat_w[4], 0.0);
        assert_eq!(power.cell_body_power_w.as_ref().unwrap()[4], 0.0);
        assert_eq!(power.cell_darcy_heat_w.as_ref().unwrap()[4], 0.0);
        near(melted.energy_j() - old.energy_j(), 200.0 + 0.1 * 800.0 * q);
        near(audit.external_heat_j, 200.0 + 0.1 * 800.0 * q);
        near(melted.cells()[&CellId(4)].specific_h_j_kg, 180.0);
        assert!(
            accepted.endpoint()[4]
                .transport
                .mobile_transport()
                .is_some()
        );
        assert_eq!(melted.cells().len(), 9);
        near(melted.mass_kg(), old.mass_kg());

        // The next applied source epoch uses the now-liquid center. The
        // asymmetric corner drive distinguishes an open center from the old
        // obstacle; it is not an interpolation of the old velocity field.
        let (next, reopened, _) = melted
            .advance_observed(&accepted, 0.11, step(0.01), || Ok(()))
            .unwrap();
        let applied = reopened.applied().unwrap();
        assert!(
            applied
                .samples
                .iter()
                .all(|sample| sample.transport.mobile_transport().is_some())
        );
        let center_flux = [
            applied.flow.field.flux_x(1, 1, 0),
            applied.flow.field.flux_x(2, 1, 0),
            applied.flow.field.flux_y(1, 1, 0),
            applied.flow.field.flux_y(1, 2, 0),
        ];
        assert!(center_flux.iter().any(|q| q.abs() > 1e-9));
        near(next.mass_kg(), old.mass_kg());
        assert!(observer.applied().is_none());
        assert_eq!(old.cells()[&CellId(4)].specific_h_j_kg, -20.0);
    }
}

#[test]
fn freezing_a_forced_cell_refuses_atomically_instead_of_discarding_its_actuator() {
    let card = source(100.0);
    let curve = chart(&card);
    let law = MushyPhaseLaw::try_new_with_immobile_solid(&curve, 1.0, 1.0, 1.0).unwrap();
    let domain = VoxelDomain::new(3, 3, 1, 1.0).unwrap();
    let old = geometry(&domain, &curve, [150.0; 9]);
    let ids = std::array::from_fn::<_, 9, _>(|c| Some(CellId(c as u64)));
    let gate = CancelGate::new();
    let config = closed();
    let observer = MaterialFlowFeedback::try_new_phase(
        &old,
        &domain,
        input(&ids, &curve),
        &config,
        CellPower {
            cell: CellId(0),
            watts: -2000.0,
        },
        &gate,
        10,
        &law,
        |_, cell: &DepositedCell, _| resolve(&card, &law, cell.specific_h_j_kg),
    )
    .unwrap()
    .with_body_force_heating(corner_force(), controls())
    .unwrap();
    let result = old.advance_observed(&observer, 0.1, step(0.1), || Ok(()));
    assert!(matches!(result, Err(MaterialFlowError::Thermal(_))));
    assert_eq!(old.time_s(), 0.0);
    assert!(
        old.cells()
            .values()
            .all(|cell| cell.specific_h_j_kg == 150.0)
    );
    assert!(observer.applied().is_none());
    assert_eq!(observer.momentum_solves(), 0);
}

#[test]
fn all_solid_advance_has_no_momentum_solve_and_rejects_incompatible_drives() {
    let card = source(100.0);
    let curve = chart(&card);
    let law = MushyPhaseLaw::try_new_with_immobile_solid(&curve, 1.0, 1.0, 1.0).unwrap();
    let domain = VoxelDomain::new(3, 3, 1, 1.0).unwrap();
    let old = geometry(&domain, &curve, [-20.0; 9]);
    let ids = std::array::from_fn::<_, 9, _>(|c| Some(CellId(c as u64)));
    let gate = CancelGate::new();
    let config = closed();
    let observer = MaterialFlowFeedback::try_new_phase(
        &old,
        &domain,
        input(&ids, &curve),
        &config,
        CellPower {
            cell: CellId(4),
            watts: 1.0,
        },
        &gate,
        1,
        &law,
        |_, cell: &DepositedCell, _| resolve(&card, &law, cell.specific_h_j_kg),
    )
    .unwrap()
    .with_viscous_heating(controls())
    .unwrap();
    let (next, accepted, audit) = old
        .advance_observed(&observer, 0.1, step(0.025), || Ok(()))
        .unwrap();
    assert_eq!(accepted.momentum_solves(), 0);
    assert_eq!(accepted.applied().unwrap().flow.report.iterations, 0);
    assert!(
        accepted
            .applied()
            .unwrap()
            .flow
            .velocity_m_s
            .iter()
            .flatten()
            .all(|&v| v == 0.0)
    );
    assert_eq!(
        accepted.viscous_heating_audit().unwrap().generated_heat_j,
        0.0
    );
    near(next.energy_j() - old.energy_j(), 0.1);
    near(audit.external_heat_j, 0.1);
    assert!(
        observer
            .clone()
            .with_body_force_heating(corner_force(), controls())
            .is_err()
    );
    for face in [
        FvBoundary::Wall {
            velocity: [1.0, 0.0, 0.0],
        },
        FvBoundary::Outlet,
    ] {
        let mut drive = closed();
        drive.faces[3] = face;
        assert!(
            MaterialFlowFeedback::try_new_phase(
                &old,
                &domain,
                input(&ids, &curve),
                &drive,
                (),
                &gate,
                1,
                &law,
                |_, cell: &DepositedCell, _| resolve(&card, &law, cell.specific_h_j_kg),
            )
            .is_err()
        );
    }
}
