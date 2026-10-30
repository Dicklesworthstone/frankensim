//! Independent discrete MHD: one shared source drives current, force and heat.
use super::*;
use fs_lbm::conjugate::enthalpy::advection::FlowEnthalpyInput;
use fs_lbm::conjugate::enthalpy::advection::feedback::mhd::{
    MaterialFlowMhdError, MhdCurrentBoundary, MhdEnthalpyInput,
    advance_mhd_enthalpy_with_material_flow,
};
use fs_lbm::conjugate::enthalpy::advection::feedback::{FlowMaterialSample, MaterialFlowFeedback};
use fs_lbm::conjugate::enthalpy::{JouleEnthalpyError, JouleEnthalpyInput};
use fs_lbm::conjugate::{ChtError, ElectricSetup, Electrode, ElectrodeDrive, ViscousHeatControls};
use fs_material::MaterialError;
use fs_material::conductor::{ElectricalPropertyBasis, ResolvedConductorState};
use fs_material::fluid_transport::{FluidTransport, FluidTransportState};
use fs_material::state_point::{
    ScalarAdmissibility, ScalarPropertyRequirement, resolve_material_state_point,
};
use std::cell::{Cell, RefCell};

const MU: f64 = 100.0;
const DT: f64 = 0.1;

fn source(resistivity: f64, temperature_feedback: bool) -> MaterialCard {
    let mut claims = ClaimSet::new();
    for (name, dims, value) in [
        (DENSITY_PROPERTY, Density::DIMS, 1.0),
        (
            SPECIFIC_HEAT_CAPACITY_PROPERTY,
            SPECIFIC_HEAT_CAPACITY_DIMS,
            100.0,
        ),
        (DYNAMIC_VISCOSITY_PROPERTY, DYNAMIC_VISCOSITY_DIMS, MU),
        (
            THERMAL_CONDUCTIVITY_PROPERTY,
            THERMAL_CONDUCTIVITY_DIMS,
            0.01,
        ),
        (
            ElectricalPropertyBasis::Resistivity.property_name(),
            ElectricalPropertyBasis::Resistivity.dimensions(),
            resistivity,
        ),
    ] {
        let mut claim = scalar(name, dims, value);
        if temperature_feedback && name == ElectricalPropertyBasis::Resistivity.property_name() {
            claim.value = PropertyValue::Curve {
                abscissa: "T".into(),
                abscissa_dims: Temperature::DIMS,
                knots: vec![
                    (300.0, resistivity),
                    (301.0, 2.0 * resistivity),
                    (810.0, 5.0 * resistivity),
                ],
                dims,
            };
            claim.interpolation = InterpolationPolicy::LinearInside;
        }
        claims.insert_claim(claim).unwrap();
    }
    card(claims)
}

#[derive(Clone, Debug, PartialEq)]
struct ConductiveLiquid {
    conductor: ResolvedConductorState,
    transport: FluidTransportState,
    query: usize,
}
impl FluidTransport for ConductiveLiquid {
    fn transport_state(&self) -> Result<FluidTransportState, MaterialError> {
        Ok(self.transport)
    }
}
fn source_error(error: impl core::fmt::Display) -> MaterialError {
    MaterialError::State {
        what: error.to_string(),
    }
}
fn resolve(
    card: &MaterialCard,
    temperature: f64,
    query: usize,
) -> Result<ConductiveLiquid, MaterialError> {
    let requirements = [
        ScalarPropertyRequirement::try_new(
            DENSITY_PROPERTY,
            Density::DIMS,
            ScalarAdmissibility::StrictlyPositive,
        )
        .unwrap(),
        ScalarPropertyRequirement::try_new(
            SPECIFIC_HEAT_CAPACITY_PROPERTY,
            SPECIFIC_HEAT_CAPACITY_DIMS,
            ScalarAdmissibility::StrictlyPositive,
        )
        .unwrap(),
        ScalarPropertyRequirement::try_new(
            DYNAMIC_VISCOSITY_PROPERTY,
            DYNAMIC_VISCOSITY_DIMS,
            ScalarAdmissibility::StrictlyPositive,
        )
        .unwrap(),
        ScalarPropertyRequirement::try_with_quantity(
            THERMAL_CONDUCTIVITY_PROPERTY,
            QuantitySpec::semantic(SemanticType::new(
                QuantityKind::ThermalConductivity,
                ValueForm::Static,
            )),
            ScalarAdmissibility::StrictlyPositive,
        )
        .unwrap(),
        ElectricalPropertyBasis::Resistivity.requirement().unwrap(),
    ];
    let material = resolve_material_state_point(
        card,
        &point(temperature, 100000.0),
        &requirements,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .map_err(source_error)?;
    let value = |name| material.property(name).unwrap().value_si();
    let transport = FluidTransportState::try_new(
        value(DENSITY_PROPERTY),
        value(SPECIFIC_HEAT_CAPACITY_PROPERTY),
        value(THERMAL_CONDUCTIVITY_PROPERTY),
        value(DYNAMIC_VISCOSITY_PROPERTY),
    )?;
    let conductor = ResolvedConductorState::try_from_material_with_basis(
        material,
        ElectricalPropertyBasis::Resistivity,
    )
    .map_err(source_error)?;
    Ok(ConductiveLiquid {
        conductor,
        transport,
        query,
    })
}
fn resistivity(sample: &FlowMaterialSample<ConductiveLiquid>) -> Result<Option<f64>, ChtError> {
    Ok(Some(sample.material.conductor.resistivity_ohm_m()))
}
fn closed_box() -> SimpleConfig {
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
fn thermal_step() -> StepConfig {
    StepConfig {
        max_dt_s: 1.0,
        energy_tolerance_j: 1e-8,
    }
}
fn setup(voltage: f64) -> ElectricSetup {
    ElectricSetup {
        conductors: vec![],
        electrodes: vec![
            Electrode {
                lo: [0, 0, 0],
                hi: [1, 2, 1],
                drive: ElectrodeDrive::Voltage(voltage),
            },
            Electrode {
                lo: [2, 0, 0],
                hi: [3, 2, 1],
                drive: ElectrodeDrive::Voltage(0.0),
            },
        ],
    }
}
fn electrical_input<'a>(
    domain: &'a VoxelDomain,
    setup: &'a ElectricSetup,
    ids: &'a [CellId],
) -> JouleEnthalpyInput<'a> {
    JouleEnthalpyInput {
        domain,
        setup,
        contacts: &[],
        cell_ids: ids,
        max_cells: ids.len(),
        max_current_solves: 300,
        power_tolerance_w: 1e-9,
        energy_tolerance_j: 1e-8,
    }
}
fn motion(field: &[[f64; 3]]) -> MhdEnthalpyInput<'_> {
    MhdEnthalpyInput {
        magnetic_flux_density_t: field,
        current_boundary: MhdCurrentBoundary::Electrodes,
        max_coupling_iterations: 40,
        relaxation: 1.0,
        velocity_tolerance_m_s: 1e-12,
    }
}
fn flow_input<'a>(
    ids: &'a [Option<CellId>],
    caloric: &'a ResolvedHeatCapacityCurve,
) -> FlowEnthalpyInput<'a> {
    FlowEnthalpyInput {
        cell_ids: ids,
        material_id: MAT,
        caloric: caloric.curve(),
        inlet_specific_h_j_kg: [None; 6],
        max_cells: ids.len(),
        max_connections: 64,
    }
}

#[derive(Clone, Default, Debug, PartialEq)]
struct Intervals {
    prepared: Vec<f64>,
    accepted: Vec<f64>,
}
impl DepositionObserver for Intervals {
    type Error = DepositionError;
    fn step_sources(
        &mut self,
        _: &DepositionSimulation<'_>,
        _: f64,
        dt: f64,
        _: f64,
        _: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<BTreeMap<CellId, f64>, Self::Error> {
        self.prepared.push(dt);
        Ok(BTreeMap::new())
    }
    fn step(
        &mut self,
        _: f64,
        dt: f64,
        _: &BTreeMap<CellId, f64>,
        _: &DepositionSimulation<'_>,
        _: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), Self::Error> {
        self.accepted.push(dt);
        Ok(())
    }
}

fn amplitude(old: f64, dt: f64, voltage: f64, b: f64, rho_e: f64, external: f64) -> f64 {
    // Symmetry reduces this 3x2 MAC box to the middle y-face amplitude q.
    // y faces=(-q/2,q,-q/2); lower x faces=(q/2,-q/2), upper x faces reversed.
    // Direct face masses give 2K=(5/2)q². Normal stress contributes 12mu q²,
    // reciprocal interior shear mu q²/2 and wall shear 3mu q², hence D=15.5mu q².
    // Both electrode rows have I=(V+Bq/2)/rho_e and W_L=-B I q.
    (2.5 * old / dt + external - b * voltage / rho_e)
        / (2.5 / dt + 15.5 * MU + b * b / (2.0 * rho_e))
}

#[test]
fn reciprocal_material_mhd_matches_independent_mac_current_work_and_heat() {
    let domain = VoxelDomain::new(3, 2, 1, 1.0).unwrap();
    let ids: Vec<_> = (0..6).map(CellId).collect();
    let flow_ids: Vec<_> = ids.iter().copied().map(Some).collect();
    let gate = CancelGate::new();
    let config = closed_box();
    let drive = setup(-1.0);
    for (rho_e, b, external) in [
        (1.0, 10.0, 0.0),
        (2.0, 10.0, 0.0),
        (1.0, -10.0, 0.0),
        (1.0, 0.0, 2.0),
    ] {
        let card = source(rho_e, false);
        let caloric = liquid_caloric(&card);
        let old = initial(&domain, 1.0, &caloric);
        let fields = vec![[0.0, 0.0, b]; 6];
        let queries = Cell::new(0);
        let reads = RefCell::new(Vec::new());
        let mut forces = vec![[0.0; 3]; 6];
        forces[1][1] = external;
        forces[4][1] = external;
        let observer = MaterialFlowFeedback::try_new(
            &old,
            &domain,
            flow_input(&flow_ids, &caloric),
            &config,
            Intervals::default(),
            &gate,
            100,
            |_, _, t| {
                queries.set(queries.get() + 1);
                resolve(&card, t, queries.get())
            },
        )
        .unwrap()
        .with_unsteady_body_force_heating(forces, controls())
        .unwrap();
        let (next, accepted, report) = advance_mhd_enthalpy_with_material_flow(
            &old,
            DT,
            thermal_step(),
            electrical_input(&domain, &drive, &ids),
            motion(&fields),
            &observer,
            |sample| {
                reads
                    .borrow_mut()
                    .push((sample.voxel, sample.material.query));
                resistivity(sample)
            },
        )
        .unwrap();
        let q = amplitude(0.0, DT, -1.0, b, rho_e, external);
        let current = (-1.0 + b * q / 2.0) / rho_e;
        let native = accepted.unsteady_flow().unwrap();
        near(native.field.flux_y(1, 1, 0), q);
        near(native.field.flux_y(0, 1, 0), -q / 2.0);
        near(native.field.flux_y(2, 1, 0), -q / 2.0);
        near(native.field.flux_x(1, 0, 0), q / 2.0);
        near(native.field.flux_x(1, 1, 0), -q / 2.0);
        near(report.electrode_energy_j, -2.0 * current * DT);
        near(report.joule_energy_j, 2.0 * rho_e * current * current * DT);
        near(report.lorentz_work_j, -b * current * q * DT);
        near(report.other_body_work_j, external * q * DT);
        near(report.mechanical.generated_heat_j, 15.5 * MU * q * q * DT);
        near(report.mechanical.final_kinetic_energy_j, 1.25 * q * q);
        near(report.mechanical.backward_euler_loss_j, 1.25 * q * q);
        near(
            next.energy_j() - old.energy_j(),
            report.joule_energy_j + report.mechanical.generated_heat_j,
        );
        near(report.combined_residual_j, 0.0);
        near(report.electrical_defect_j, 0.0);
        near(report.mechanical_transfer_defect_j, 0.0);
        near(report.native_body_transfer_defect_j, 0.0);
        near(report.heat_transfer_defect_j, 0.0);
        near(report.other_thermal_input_j, 0.0);
        assert!(report.momentum_trials > 1);
        assert_eq!(report.current_solves, report.momentum_trials + 1);
        assert_eq!(accepted.momentum_solves(), report.momentum_trials);
        assert_eq!(report.thermal.substeps, 1);
        assert_eq!(report.mechanical.substeps, 1);
        assert_eq!(accepted.inner().prepared, [DT]);
        assert_eq!(accepted.inner().accepted, [DT]);
        // Six constructor, six actual start and six endpoint material queries:
        // current Picard iterations read those retained source objects only.
        assert_eq!(queries.get(), 18);
        for (voxel, query) in reads.borrow().iter() {
            let start = accepted.applied().unwrap().samples[*voxel].material.query;
            let end = accepted.endpoint()[*voxel].material.query;
            assert!(*query == start || *query == end);
        }
        for sample in accepted.endpoint() {
            assert_eq!(
                sample.temperature_k.to_bits(),
                next.temperature_k(sample.cell).unwrap().to_bits()
            );
        }
        for c in [0, 2, 3, 5] {
            assert_eq!(
                report
                    .last_applied
                    .as_ref()
                    .unwrap()
                    .solution
                    .electric
                    .joule_w[c],
                0.0
            );
            assert_eq!(
                report
                    .last_applied
                    .as_ref()
                    .unwrap()
                    .solution
                    .lorentz_force_density_n_m3[c],
                [0.0; 3]
            );
        }
        assert_eq!(observer.unsteady_time_s(), Some(0.0));
        assert_eq!(observer.inner(), &Intervals::default());
        assert_eq!(old.energy_j(), 0.0);
    }
}

#[test]
fn zero_voltage_brakes_retained_motion_and_heats_from_kinetic_energy() {
    let domain = VoxelDomain::new(3, 2, 1, 1.0).unwrap();
    let ids: Vec<_> = (0..6).map(CellId).collect();
    let flow_ids: Vec<_> = ids.iter().copied().map(Some).collect();
    let card = source(1.0, false);
    let caloric = liquid_caloric(&card);
    let old = initial(&domain, 1.0, &caloric);
    let config = closed_box();
    let gate = CancelGate::new();
    let fields = vec![[0.0, 0.0, 10.0]; 6];
    let drive = setup(-1.0);
    let observer = MaterialFlowFeedback::try_new(
        &old,
        &domain,
        flow_input(&flow_ids, &caloric),
        &config,
        Intervals::default(),
        &gate,
        100,
        |_, _, t| resolve(&card, t, 0),
    )
    .unwrap()
    .with_unsteady_viscous_heating(controls())
    .unwrap();
    let (spun, moving, first) = advance_mhd_enthalpy_with_material_flow(
        &old,
        DT,
        thermal_step(),
        electrical_input(&domain, &drive, &ids),
        motion(&fields),
        &observer,
        resistivity,
    )
    .unwrap();
    let short = setup(0.0);
    let (next, stopped, report) = advance_mhd_enthalpy_with_material_flow(
        &spun,
        2.0 * DT,
        thermal_step(),
        electrical_input(&domain, &short, &ids),
        motion(&fields),
        &moving,
        resistivity,
    )
    .unwrap();
    let q0 = amplitude(0.0, DT, -1.0, 10.0, 1.0, 0.0);
    let q = amplitude(q0, DT, 0.0, 10.0, 1.0, 0.0);
    let current = 5.0 * q;
    near(stopped.unsteady_flow().unwrap().field.flux_y(1, 1, 0), q);
    assert!(q > 0.0 && q < q0);
    assert!(report.lorentz_work_j < 0.0);
    assert!(report.joule_energy_j > 0.0);
    near(report.electrode_energy_j, 0.0);
    near(report.other_body_work_j, 0.0);
    near(report.joule_energy_j, 2.0 * current * current * DT);
    near(report.lorentz_work_j, -report.joule_energy_j);
    near(
        report.mechanical.initial_kinetic_energy_j,
        first.mechanical.final_kinetic_energy_j,
    );
    near(report.mechanical.final_kinetic_energy_j, 1.25 * q * q);
    near(
        report.mechanical.backward_euler_loss_j,
        1.25 * (q - q0) * (q - q0),
    );
    near(
        next.energy_j() - spun.energy_j(),
        report.joule_energy_j + report.mechanical.generated_heat_j,
    );
    near(
        first.mechanical.final_kinetic_energy_j - report.mechanical.final_kinetic_energy_j,
        next.energy_j() - spun.energy_j() + report.mechanical.backward_euler_loss_j,
    );
    near(report.combined_residual_j, 0.0);
    assert_eq!(stopped.inner().prepared.len(), 2);
    assert_eq!(stopped.inner().accepted.len(), 2);
    assert_eq!(moving.inner().accepted, [DT]);
    near(moving.unsteady_flow().unwrap().field.flux_y(1, 1, 0), q0);
}

#[test]
fn shared_temperature_feedback_changes_current_and_lorentz_without_changing_physics() {
    let domain = VoxelDomain::new(3, 2, 1, 1.0).unwrap();
    let ids: Vec<_> = (0..6).map(CellId).collect();
    let flow_ids: Vec<_> = ids.iter().copied().map(Some).collect();
    let card = source(1.0, true);
    let caloric = liquid_caloric(&card);
    let old = initial(&domain, 1.0, &caloric);
    let config = closed_box();
    let gate = CancelGate::new();
    let fields = vec![[0.0, 0.0, 10.0]; 6];
    let drive = setup(-1.0);
    let observer = MaterialFlowFeedback::try_new(
        &old,
        &domain,
        flow_input(&flow_ids, &caloric),
        &config,
        (),
        &gate,
        100,
        |_, _, t| resolve(&card, t, 0),
    )
    .unwrap()
    .with_unsteady_viscous_heating(controls())
    .unwrap();
    let (warm, moving, first) = advance_mhd_enthalpy_with_material_flow(
        &old,
        DT,
        thermal_step(),
        electrical_input(&domain, &drive, &ids),
        motion(&fields),
        &observer,
        resistivity,
    )
    .unwrap();
    let (next, accepted, second) = advance_mhd_enthalpy_with_material_flow(
        &warm,
        2.0 * DT,
        thermal_step(),
        electrical_input(&domain, &drive, &ids),
        motion(&fields),
        &moving,
        resistivity,
    )
    .unwrap();
    let r0 = 1.0 + warm.temperature_k(CellId(1)).unwrap() - 300.0;
    let r1 = 1.0 + warm.temperature_k(CellId(4)).unwrap() - 300.0;
    assert!(r0 > 1.0 && r1 > 1.0);
    let q0 = amplitude(0.0, DT, -1.0, 10.0, 1.0, 0.0);
    let g = 1.0 / r0 + 1.0 / r1;
    let q = (2.5 * q0 / DT + 5.0 * g) / (2.5 / DT + 15.5 * MU + 25.0 * g);
    let current_sum = (-1.0 + 5.0 * q) * g;
    near(accepted.unsteady_flow().unwrap().field.flux_y(1, 1, 0), q);
    near(second.electrode_energy_j, -current_sum * DT);
    near(second.joule_energy_j, (-1.0 + 5.0 * q).powi(2) * g * DT);
    let applied = second.last_applied.as_ref().unwrap();
    near(applied.resistivity_ohm_m[1].unwrap(), r0);
    near(applied.resistivity_ohm_m[4].unwrap(), r1);
    assert_eq!(applied.temperature_k, first.endpoint.temperature_k);
    assert_eq!(
        applied.specific_enthalpy_j_kg,
        first.endpoint.specific_enthalpy_j_kg
    );
    for c in [1, 4] {
        near(
            second.endpoint.resistivity_ohm_m[c].unwrap(),
            1.0 + next.temperature_k(CellId(c as u64)).unwrap() - 300.0,
        );
        assert!(second.endpoint.resistivity_ohm_m[c] > applied.resistivity_ohm_m[c]);
    }
    near(second.combined_residual_j, 0.0);
}

#[test]
fn late_current_budget_and_nonconvergent_transfer_publish_no_history() {
    let domain = VoxelDomain::new(3, 2, 1, 1.0).unwrap();
    let ids: Vec<_> = (0..6).map(CellId).collect();
    let flow_ids: Vec<_> = ids.iter().copied().map(Some).collect();
    let card = source(1.0, false);
    let caloric = liquid_caloric(&card);
    let old = initial(&domain, 1.0, &caloric);
    let config = closed_box();
    let gate = CancelGate::new();
    let fields = vec![[0.0, 0.0, 10.0]; 6];
    let drive = setup(-1.0);
    let observer = MaterialFlowFeedback::try_new(
        &old,
        &domain,
        flow_input(&flow_ids, &caloric),
        &config,
        Intervals::default(),
        &gate,
        100,
        |_, _, t| resolve(&card, t, 0),
    )
    .unwrap()
    .with_unsteady_viscous_heating(controls())
    .unwrap();
    let (_, _, report) = advance_mhd_enthalpy_with_material_flow(
        &old,
        DT,
        thermal_step(),
        electrical_input(&domain, &drive, &ids),
        motion(&fields),
        &observer,
        resistivity,
    )
    .unwrap();
    let mut limited = electrical_input(&domain, &drive, &ids);
    limited.max_current_solves = report.current_solves - 1; // only endpoint support fails
    assert!(matches!(
        advance_mhd_enthalpy_with_material_flow(
            &old,
            DT,
            thermal_step(),
            limited,
            motion(&fields),
            &observer,
            resistivity
        ),
        Err(MaterialFlowMhdError::Electrical(
            JouleEnthalpyError::CurrentBudget
        ))
    ));
    let mut one_trial = motion(&fields);
    one_trial.max_coupling_iterations = 1;
    assert!(matches!(
        advance_mhd_enthalpy_with_material_flow(
            &old,
            DT,
            thermal_step(),
            electrical_input(&domain, &drive, &ids),
            one_trial,
            &observer,
            resistivity
        ),
        Err(MaterialFlowMhdError::CouplingDidNotConverge { .. })
    ));
    assert_eq!(observer.inner(), &Intervals::default());
    assert_eq!(observer.unsteady_time_s(), Some(0.0));
    assert_eq!(observer.unsteady_heating_audit().unwrap().substeps, 0);
    assert!(
        observer
            .unsteady_flow()
            .unwrap()
            .velocity_m_s
            .iter()
            .all(|u| *u == [0.0; 3])
    );
    assert_eq!(old.time_s(), 0.0);
    assert_eq!(old.energy_j(), 0.0);
}

mod interfaces;
