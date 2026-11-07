//! G0/G1/G3: original density changes resistance without manufacturing mass.
use super::*;

struct Expanding<'a> {
    queries: &'a Cell<usize>,
    maximum_h: f64,
    constant_density: bool,
}

impl ThermalLaw for Expanding<'_> {
    fn temperature_k(&self, h: f64) -> Result<f64, DepositionError> {
        SENSIBLE.temperature_k(h)
    }
    fn maximum_temperature_slope(&self) -> f64 {
        1.0
    }
    fn bulk_density_kg_m3(&self, h: f64) -> Result<Option<f64>, DepositionError> {
        self.queries.set(self.queries.get() + 1);
        if !positive(h) || h > self.maximum_h {
            return Err(DepositionError::MaterialDomain);
        }
        Ok(Some(if self.constant_density {
            1.0
        } else {
            300.0 / h
        }))
    }
}

fn density_layer() -> PlanarDensityLayer {
    PlanarDensityLayer {
        cell_id: CellId(0),
        footprint_area_m2: 1.0,
    }
}

#[test]
fn heating_refreshes_normal_resistance_and_preserves_original_mass_and_storage() {
    let queries = Cell::new(0);
    let law = Expanding {
        queries: &queries,
        maximum_h: 1000.0,
        constant_density: false,
    };
    let initial = simulation(cells(&[300.0], 1.0, Some(600.0)), &law, vec![], 100);
    let fixed = ConductionObserver::new(
        &initial,
        &[(MaterialId(1), &UNIT_K)],
        &[],
        &[face(1.0, 0.0)],
        || Ok(()),
    )
    .unwrap();
    let expanding = fixed
        .clone()
        .with_material_density_layers(&initial, &[density_layer()], || Ok(()))
        .unwrap();
    let (heated, observer, audit) = initial
        .advance_observed(&expanding, 0.2, config(0.1), || Ok(()))
        .unwrap();
    // rho=300/T, d=T/600: the two actual steps use G=2 then G=5/3.
    // A stale density instead gives T=408 K after the same two steps.
    near(heated.temperature_k(CellId(0)).unwrap(), 400.0, 1e-12);
    let (stale, _, _) = initial
        .advance_observed(&fixed, 0.2, config(0.1), || Ok(()))
        .unwrap();
    near(stale.temperature_k(CellId(0)).unwrap(), 408.0, 1e-12);
    assert_eq!(heated.cells()[&CellId(0)].mass_kg, 1.0);
    assert_eq!(
        heated.cells()[&CellId(0)].geometry,
        initial.cells()[&CellId(0)].geometry
    );
    near(heated.energy_j() - initial.energy_j(), 100.0, 1e-12);
    near(audit.residual_j, 0.0, 1e-12);
    near(
        observer.last_planar_half_distance_m(CellId(0)).unwrap(),
        0.6,
        1e-15,
    );
    near(
        heated
            .power_balance_observed(&observer, || Ok(()))
            .unwrap()
            .external_power_w,
        300.0,
        1e-12,
    );
}

#[test]
fn one_density_evaluation_serves_reservoir_and_ambient_faces() {
    let queries = Cell::new(0);
    let law = Expanding {
        queries: &queries,
        maximum_h: 1000.0,
        constant_density: false,
    };
    let mut deposited = cells(&[300.0], 1.0, Some(500.0));
    deposited[0].boundary.ambient_k = 450.0;
    let initial = simulation(deposited, &law, vec![], 100);
    let observer = ConductionObserver::new(
        &initial,
        &[(MaterialId(1), &UNIT_K)],
        &[],
        &[face(1.0, 0.0)],
        || Ok(()),
    )
    .unwrap()
    .with_ambient_faces(
        &initial,
        &[ConductiveAmbient {
            cell_id: CellId(0),
            half_distance_m: 0.5,
            surface_resistance_m2_k_w: 0.0,
        }],
        || Ok(()),
    )
    .unwrap()
    .with_material_density_layers(&initial, &[density_layer()], || Ok(()))
    .unwrap();
    queries.set(0);
    let power = initial
        .power_balance_observed(&observer, || Ok(()))
        .unwrap();
    assert_eq!(queries.get(), 1, "query density per cell, not per face");
    near(power.external_power_w, 700.0, 1e-12);
    let (heated, observer, _) = initial
        .advance_observed(&observer, 0.1, config(0.1), || Ok(()))
        .unwrap();
    near(heated.temperature_k(CellId(0)).unwrap(), 370.0, 1e-12);
    queries.set(0);
    let power = heated.power_balance_observed(&observer, || Ok(())).unwrap();
    near(
        power.external_power_w,
        (600.0 / 370.0) * (130.0 + 80.0),
        1e-12,
    );
    assert_eq!(queries.get(), 1);
}

#[test]
fn constant_original_density_reproduces_the_existing_phase_layer() {
    let queries = Cell::new(0);
    let law = Expanding {
        queries: &queries,
        maximum_h: 1000.0,
        constant_density: true,
    };
    let initial = simulation(cells(&[300.0, 400.0], 1.0, None), &law, vec![], 100);
    let fixed = ConductionObserver::new(
        &initial,
        &[(MaterialId(1), &UNIT_K)],
        &link_geometry(2, 1.0),
        &[],
        || Ok(()),
    )
    .unwrap();
    let sourced = fixed
        .clone()
        .with_material_density_layers(&initial, &[density_layer()], || Ok(()))
        .unwrap();
    let explicit = fixed
        .with_planar_phase_layers(
            &initial,
            &[PlanarPhaseLayer {
                cell_id: CellId(0),
                footprint_area_m2: 1.0,
                density_kg_m3: 1.0,
            }],
            || Ok(()),
        )
        .unwrap();
    let (a, _, _) = initial
        .advance_observed(&sourced, 0.3, config(0.1), || Ok(()))
        .unwrap();
    let (b, _, _) = initial
        .advance_observed(&explicit, 0.3, config(0.1), || Ok(()))
        .unwrap();
    assert_eq!(a.cells(), b.cells());
    assert_eq!(a.energy_j(), initial.energy_j());
}

#[test]
fn absent_density_and_endpoint_domain_exit_refuse_without_publishing() {
    let unavailable = simulation(cells(&[300.0], 1.0, Some(600.0)), &SENSIBLE, vec![], 100);
    let observer = ConductionObserver::new(
        &unavailable,
        &[(MaterialId(1), &UNIT_K)],
        &[],
        &[face(1.0, 0.0)],
        || Ok(()),
    )
    .unwrap();
    assert!(
        observer
            .with_material_density_layers(&unavailable, &[density_layer()], || Ok(()))
            .is_err()
    );

    let queries = Cell::new(0);
    let law = Expanding {
        queries: &queries,
        maximum_h: 350.0,
        constant_density: false,
    };
    let initial = simulation(cells(&[300.0], 1.0, Some(600.0)), &law, vec![], 100);
    let observer = ConductionObserver::new(
        &initial,
        &[(MaterialId(1), &UNIT_K)],
        &[],
        &[face(1.0, 0.0)],
        || Ok(()),
    )
    .unwrap()
    .with_material_density_layers(&initial, &[density_layer()], || Ok(()))
    .unwrap();
    // The final caloric T=360 is supported, but its density is not. There is
    // no next substep to discover that invalid endpoint: acceptance must refuse.
    assert!(
        initial
            .advance_observed(&observer, 0.1, config(0.1), || Ok(()))
            .is_err()
    );
    assert_eq!(initial.temperature_k(CellId(0)).unwrap(), 300.0);
    assert_eq!(observer.last_planar_half_distance_m(CellId(0)), Some(0.5));
    assert!(
        initial
            .advance_observed(&observer, 0.05, config(0.1), || Ok(()))
            .is_ok()
    );
    assert!(
        observer
            .clone()
            .with_material_density_layers(&initial, &[density_layer()], || Ok(()))
            .is_err()
    );
    assert!(
        observer
            .clone()
            .with_material_density_layers(&initial, &[], || Err(DepositionError::Cancelled))
            .is_err()
    );
}
