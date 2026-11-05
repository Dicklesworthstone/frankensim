//! Surface convection changes the original latent heat and native flow history.
use super::*;
use fs_evidence::{Evidence, ModelEvidence, NumericalCertificate, ProvenanceHash, StatisticalCertificate};
use fs_lbm::conjugate::enthalpy::advection::feedback::environment::{
    FiniteGasConvection, FiniteGasConvectionContact, FiniteGasConvectionQuery,
    FiniteGasHeatTransferCoefficient, advance_finite_gas_with_material_convection,
};
use fs_material::thermal::conduction::LayerConductionConfig;
use fs_material::MaterialError;
use std::cell::RefCell;

#[derive(Debug)]
struct Film {
    h: f64,
    minimum_bulk_k: Option<f64>,
    variable: bool,
    calls: RefCell<Vec<[f64; 3]>>,
}
impl Film {
    fn new(h: f64) -> Self {
        Self { h, minimum_bulk_k: None, variable: false, calls: RefCell::new(Vec::new()) }
    }
}
impl FiniteGasConvection for Film {
    fn coefficient(&self, query: FiniteGasConvectionQuery<'_>)
        -> Result<Option<Evidence<FiniteGasHeatTransferCoefficient>>, MaterialError>
    {
        let bulk = query.bulk.state();
        let film = query.film.state();
        assert_eq!(query.contact.outward_normal, [0.0, 0.0, 1.0]);
        assert_eq!(bulk.pressure.to_bits(), film.pressure.to_bits());
        near(film.temperature, 0.5 * query.surface_temperature_k + 0.5 * bulk.temperature);
        self.calls.borrow_mut().push([
            query.surface_temperature_k, query.condensed_temperature_k, bulk.pressure,
        ]);
        if self.minimum_bulk_k.is_some_and(|limit| bulk.temperature < limit) {
            return Err(MaterialError::State { what: "manufactured late convection source refusal".into() });
        }
        if query.surface_temperature_k == bulk.temperature { return Ok(None); }
        let h = if self.variable { self.h * query.surface_temperature_k / 500.0 } else { self.h };
        Ok(Some(Evidence {
            value: FiniteGasHeatTransferCoefficient::new(h),
            qoi: h,
            numerical: NumericalCertificate::exact(h),
            statistical: StatisticalCertificate::None,
            model: ModelEvidence::none(),
            sensitivity: Default::default(),
            provenance: ProvenanceHash(0x5048_4153_4546_494c),
            adjoint_ref: None,
        }))
    }
}
fn surface_controls() -> LayerConductionConfig {
    LayerConductionConfig { max_iterations: 96, absolute_heat_tolerance_w: 1e-11, relative_heat_tolerance: 1e-12 }
}
fn convective_contact(id: u64) -> FiniteGasConvectionContact {
    FiniteGasConvectionContact {
        contact: FiniteGasContact { distances_m: [0.1, 0.0], ..contact(id) },
        outward_normal: [0.0, 0.0, 1.0],
    }
}

#[test]
fn actual_surface_and_nasa_energy_close_in_both_heat_directions_including_latent_storage() {
    let card = source(100.0);
    let curve = chart(&card);
    let law = MushyPhaseLaw::try_new_with_immobile_solid(&curve, 1.0, 1.0, 1.0).unwrap();
    let domain = VoxelDomain::new(3, 3, 1, 1.0).unwrap();
    let ids = std::array::from_fn::<_, 9, _>(|c| Some(CellId(c as u64)));
    let gate = CancelGate::new();
    let config = closed();
    for (initial_h, cell_t, k, gas_t) in [(-20.0, 490.0, 0.1, 600.0), (300.0, 600.0, 0.2, 350.0)] {
        let old = geometry(&domain, &curve, [initial_h; 9]);
        let flow = MaterialFlowFeedback::try_new_phase(
            &old, &domain, input(&ids, &curve), &config, (), &gate, 16, &law,
            |_, cell: &DepositedCell, _| resolve(&card, &law, cell.specific_h_j_kg),
        ).unwrap().with_unsteady_momentum().unwrap()
            .with_unsteady_viscous_heating(controls()).unwrap();
        let gas = IsochoricGasParcel::new(
            gas_source::model(12000.0, 250000.0, 100.0),
            gas_source::point(gas_t, 100000.0), 0.1,
        ).unwrap();
        let film = Film::new(5.0);
        let contacts: Vec<_> = (0..9).map(convective_contact).collect();
        let dt = 0.25;
        let (next, accepted, changed, report) = advance_finite_gas_with_material_convection(
            &old, dt, step(1.0), &flow, &gas, &contacts, &film,
            gas_controls(), surface_controls(),
        ).unwrap();
        // Independent two-resistance heat and NASA cv integral; no product
        // surface, conductance or reported rate is used to construct this oracle.
        let conductance = 1.0 / (0.1 / k + 1.0 / film.h);
        let q_cell = dt * conductance * (gas_t - cell_t);
        let q = 9.0 * q_cell;
        let target = 2.0 * gas_t + 0.0005 * gas_t * gas_t - q / gas.mass_kg() * 0.028 / R;
        let expected_gas = ((4.0 + 0.002 * target).sqrt() - 2.0) / 0.001;
        let expected_surface = ((k / 0.1) * cell_t + film.h * gas_t) / (k / 0.1 + film.h);
        assert_eq!(report.substeps, 1);
        near(changed.state().state().temperature, expected_gas);
        near(changed.state().state().pressure, 100000.0 * expected_gas / gas_t);
        near(report.heat_into_condensed_j, q);
        near(report.gas_internal_energy_change_j, -q);
        near(next.energy_j() - old.energy_j(), q);
        near(report.other_thermal_external_heat_j, 0.0);
        assert!(report.combined_thermal_balance_defect_j.abs() < 1e-7);
        let applied = &report.last_applied[0];
        let boundary = applied.convection.as_ref().unwrap();
        near(boundary.surface_temperature_k, expected_surface);
        near(applied.conductance_w_k, conductance);
        near(boundary.condensed_power_to_gas_w + applied.heat_into_condensed_w,
            boundary.surface_heat_residual_w);
        for (sample, endpoint) in accepted.endpoint().iter().zip(&report.endpoint) {
            near(sample.specific_enthalpy_j_kg, initial_h + q_cell);
            let boundary = endpoint.convection.as_ref().unwrap();
            near(endpoint.film.state().temperature,
                0.5 * boundary.surface_temperature_k + 0.5 * expected_gas);
            near(endpoint.film.state().pressure, changed.state().state().pressure);
            if initial_h < 0.0 {
                near(sample.temperature_k, 500.0);
                assert!(sample.transport.mobile_transport().is_some());
                near(endpoint.condensed_conductivity_w_m_k, 0.1 * (1.0 + (initial_h + q_cell) / 100.0));
            }
        }
        assert!(film.calls.borrow().iter().any(|q| q[0] != q[1]));
        assert!(film.calls.borrow().iter().any(|q| q[2] != 100000.0));
        assert_eq!(flow.unsteady_time_s(), Some(0.0));
        assert_eq!(old.time_s(), 0.0);
        assert_eq!(gas.state().state().temperature, gas_t);
    }
}

#[test]
fn convection_strength_opens_original_phase_faces_and_changes_native_momentum() {
    let card = source(100.0);
    let curve = chart(&card);
    let law = MushyPhaseLaw::try_new_with_immobile_solid(&curve, 1.0, 1.0, 1.0).unwrap();
    let domain = VoxelDomain::new(3, 3, 1, 1.0).unwrap();
    let ids = std::array::from_fn::<_, 9, _>(|c| Some(CellId(c as u64)));
    let gate = CancelGate::new();
    let config = closed();
    let mut enthalpy = [50.0; 9];
    enthalpy[4] = -20.0;
    let old = geometry(&domain, &curve, enthalpy);
    let flow = MaterialFlowFeedback::try_new_phase(
        &old, &domain, input(&ids, &curve), &config, (), &gate, 20, &law,
        |_, cell: &DepositedCell, _| resolve(&card, &law, cell.specific_h_j_kg),
    ).unwrap().with_unsteady_body_force_heating(corner_force(), controls()).unwrap();
    let gas = IsochoricGasParcel::new(gas_source::model(0.0, 250000.0, 100.0),
        gas_source::point(600.0, 100000.0), 1.0).unwrap();
    let mut outputs = Vec::new();
    for h in [0.01, 50.0] {
        let film = Film::new(h);
        let (next, accepted, _, report) = advance_finite_gas_with_material_convection(
            &old, 0.3, step(0.3), &flow, &gas, &[convective_contact(4)], &film,
            gas_controls(), surface_controls(),
        ).unwrap();
        let field = &accepted.unsteady_flow().unwrap().field;
        let flux = [field.flux_x(1, 1, 0), field.flux_x(2, 1, 0),
            field.flux_y(1, 1, 0), field.flux_y(1, 2, 0)]
            .into_iter().map(f64::abs).fold(0.0, f64::max);
        outputs.push((accepted.endpoint()[4].transport.mobile_transport().is_some(), flux,
            next.cells()[&CellId(4)].specific_h_j_kg));
        near(report.other_thermal_external_heat_j, accepted.unsteady_heating_audit().unwrap().generated_heat_j);
        assert!(report.combined_thermal_balance_defect_j.abs() < 1e-7);
    }
    assert!(!outputs[0].0);
    assert_eq!(outputs[0].1, 0.0);
    assert!(outputs[1].0);
    assert!(outputs[1].1 > 1e-10);
    assert!(outputs[1].2 > outputs[0].2);
}

#[test]
fn geometry_surface_budget_and_late_provider_failures_preserve_every_input() {
    let card = source(100.0);
    let curve = chart(&card);
    let law = MushyPhaseLaw::try_new_with_immobile_solid(&curve, 1.0, 1.0, 1.0).unwrap();
    let domain = VoxelDomain::new(3, 3, 1, 1.0).unwrap();
    let ids = std::array::from_fn::<_, 9, _>(|c| Some(CellId(c as u64)));
    let gate = CancelGate::new();
    let config = closed();
    let old = geometry(&domain, &curve, [-20.0; 9]);
    let flow = MaterialFlowFeedback::try_new_phase(
        &old, &domain, input(&ids, &curve), &config, (), &gate, 20, &law,
        |_, cell: &DepositedCell, _| resolve(&card, &law, cell.specific_h_j_kg),
    ).unwrap().with_unsteady_momentum().unwrap();
    let gas = IsochoricGasParcel::new(gas_source::model(0.0, 250000.0, 100.0),
        gas_source::point(600.0, 100000.0), 0.1).unwrap();
    let original = gas.clone();
    let film = Film::new(50.0);
    let late = Film { minimum_bulk_k: Some(599.99), ..Film::new(50.0) };
    let nonlinear = Film { variable: true, ..Film::new(50.0) };
    let mut bad_film = convective_contact(0);
    bad_film.contact.distances_m[1] = 0.1;
    let mut bad_normal = convective_contact(0);
    bad_normal.outward_normal = [0.0; 3];
    let exhausted = LayerConductionConfig {
        max_iterations: 1, absolute_heat_tolerance_w: 0.0, relative_heat_tolerance: 0.0,
    };
    for (index, (contact, provider, controls)) in [
        (bad_film, &film, surface_controls()),
        (bad_normal, &film, surface_controls()),
        (convective_contact(0), &nonlinear, exhausted),
        (convective_contact(0), &late, surface_controls()),
    ].into_iter().enumerate() {
        let error = advance_finite_gas_with_material_convection(
            &old, 0.25, step(1.0), &flow, &gas, &[contact], provider, gas_controls(), controls,
        ).err().expect("the physical source or surface budget must refuse");
        if index == 2 { assert!(matches!(error, FiniteGasFlowError::SurfaceBalance { .. })); }
        if index == 3 { assert!(matches!(error, FiniteGasFlowError::Flow(MaterialFlowError::Material { .. }))); }
        assert_eq!(gas, original);
        assert_eq!(old.time_s(), 0.0);
        assert!(old.cells().values().all(|cell| cell.specific_h_j_kg == -20.0));
        assert_eq!(flow.unsteady_time_s(), Some(0.0));
    }
    gate.request();
    let cancelled = advance_finite_gas_with_material_convection(
        &old, 0.25, step(1.0), &flow, &gas, &[convective_contact(0)], &film,
        gas_controls(), surface_controls(),
    );
    assert!(matches!(cancelled,
        Err(FiniteGasFlowError::Flow(MaterialFlowError::Thermal(DepositionError::Cancelled)))));
    assert_eq!(gas, original);
}

#[test]
fn equilibrium_keeps_latent_fraction_and_gas_energy_without_inventing_a_coefficient() {
    let card = source(100.0);
    let curve = chart(&card);
    let law = MushyPhaseLaw::try_new_with_immobile_solid(&curve, 1.0, 1.0, 1.0).unwrap();
    let domain = VoxelDomain::new(3, 3, 1, 1.0).unwrap();
    let ids = std::array::from_fn::<_, 9, _>(|c| Some(CellId(c as u64)));
    let gate = CancelGate::new();
    let config = closed();
    let old = geometry(&domain, &curve, [50.0; 9]);
    let flow = MaterialFlowFeedback::try_new_phase(
        &old, &domain, input(&ids, &curve), &config, (), &gate, 20, &law,
        |_, cell: &DepositedCell, _| resolve(&card, &law, cell.specific_h_j_kg),
    ).unwrap().with_unsteady_momentum().unwrap();
    let gas = IsochoricGasParcel::new(gas_source::model(0.0, 250000.0, 100.0),
        gas_source::point(500.0, 100000.0), 0.1).unwrap();
    let film = Film::new(50.0);
    let (next, _, unchanged, report) = advance_finite_gas_with_material_convection(
        &old, 0.01, step(0.01), &flow, &gas, &[convective_contact(0)], &film,
        gas_controls(), surface_controls(),
    ).unwrap();
    assert_eq!(unchanged, gas);
    assert_eq!(report.heat_into_condensed_j, 0.0);
    assert!(next.cells().values().all(|cell| cell.specific_h_j_kg == 50.0));
    let sample = &report.endpoint[0];
    assert_eq!(sample.conductance_w_k, 0.0);
    assert_eq!(sample.convection.as_ref().unwrap().surface_temperature_k, 500.0);
    assert_eq!(sample.convection.as_ref().unwrap().coefficient, None);
}
