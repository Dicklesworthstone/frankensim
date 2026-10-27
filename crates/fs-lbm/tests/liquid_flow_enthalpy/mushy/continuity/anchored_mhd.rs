//! Anchored phase MHD keeps current, force, heat and original material cells.
use super::*;
use fs_lbm::conjugate::enthalpy::advection::feedback::MaterialFlowError;
use fs_lbm::conjugate::enthalpy::advection::feedback::mhd::{
    MaterialFlowMhdError, MhdCurrentBoundary, MhdEnthalpyInput,
    advance_mhd_enthalpy_with_anchored_phase, advance_mhd_enthalpy_with_material_flow,
};

fn motion(field: &[[f64; 3]]) -> MhdEnthalpyInput<'_> {
    MhdEnthalpyInput {
        magnetic_flux_density_t: field,
        current_boundary: MhdCurrentBoundary::Electrodes,
        max_coupling_iterations: 40,
        relaxation: 1.0,
        velocity_tolerance_m_s: 1e-12,
    }
}

fn interval() -> StepConfig {
    StepConfig {
        max_dt_s: 0.1,
        energy_tolerance_j: 1e-8,
    }
}

fn flow_input<'a>(
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

// The same phase-dependent card data as the parent fixture, with electrical
// resistivity also varied. Liquid-only properties are absent in exact solids.
fn electrical_source(mu: f64, rho_e: f64) -> MaterialCard {
    let mut claims = ClaimSet::new();
    for (name, dims, value) in [
        (DENSITY_PROPERTY, Density::DIMS, 1.0),
        (SPECIFIC_HEAT_CAPACITY_PROPERTY, SPECIFIC_HEAT_CAPACITY_DIMS, 2.0),
        (DYNAMIC_VISCOSITY_PROPERTY, DYNAMIC_VISCOSITY_DIMS, mu),
        (THERMAL_CONDUCTIVITY_PROPERTY, THERMAL_CONDUCTIVITY_DIMS, 0.1),
        (
            ElectricalPropertyBasis::Resistivity.property_name(),
            ElectricalPropertyBasis::Resistivity.dimensions(),
            rho_e,
        ),
    ] {
        let mut claim = scalar(name, dims, value);
        let liquid_only =
            name == SPECIFIC_HEAT_CAPACITY_PROPERTY || name == DYNAMIC_VISCOSITY_PROPERTY;
        claim.validity = claim.validity.with(
            "liquid",
            if liquid_only { f64::from_bits(1) } else { 0.0 },
            1.0,
        );
        if name == THERMAL_CONDUCTIVITY_PROPERTY
            || name == ElectricalPropertyBasis::Resistivity.property_name()
        {
            claim.value = PropertyValue::Curve {
                abscissa: "liquid".into(),
                abscissa_dims: Dims([0; 6]),
                knots: vec![(0.0, value), (1.0, 2.0 * value)],
                dims,
            };
            claim.interpolation = InterpolationPolicy::LinearInside;
        }
        claims.insert_claim(claim).unwrap();
    }
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: "manufactured phase-change conductor".into(),
            phase: "solid-liquid".into(),
            process: "anchored MHD coupling fixture".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap()
}

#[derive(Clone, Default, Debug, PartialEq)]
struct HeatThenCool {
    hot_w: f64,
    cold_w: f64,
    accepted: Vec<f64>,
}
impl HeatThenCool {
    fn melting() -> Self {
        Self { hot_w: 1500.0, cold_w: -2000.0, accepted: Vec::new() }
    }
}
impl DepositionObserver for HeatThenCool {
    type Error = DepositionError;

    fn event(&mut self, _: &Event, _: &DepositionSimulation<'_>) -> Result<(), Self::Error> {
        Ok(())
    }

    fn sources(
        &mut self,
        thermal: &DepositionSimulation<'_>,
        _: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<BTreeMap<CellId, f64>, Self::Error> {
        // Changes in this external drive coincide with explicit advance targets.
        let watts = if thermal.time_s() < 0.1 {
            self.hot_w
        } else if thermal.time_s() < 0.2 {
            0.0
        } else {
            self.cold_w
        };
        Ok(thermal.cells().keys().map(|&id| (id, watts)).collect())
    }

    fn step(
        &mut self,
        _: f64,
        dt_s: f64,
        _: &BTreeMap<CellId, f64>,
        _: &DepositionSimulation<'_>,
        _: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), Self::Error> {
        self.accepted.push(dt_s);
        Ok(())
    }
}

#[test]
fn exact_solid_current_heats_and_reports_its_equal_opposite_lorentz_support_load() {
    let domain = VoxelDomain::new(3, 3, 1, 1.0).unwrap();
    let ids = std::array::from_fn::<_, 9, _>(|c| CellId(c as u64));
    let flow_ids = ids.map(Some);
    let gate = CancelGate::new();
    let config = closed();
    let drive = setup(10.0);
    for (rho_e, b) in [(1.0, 5.0), (2.0, 5.0), (1.0, -5.0)] {
        let card = electrical_source(100.0, rho_e);
        let curve = chart(&card);
        let law = MushyPhaseLaw::try_new_with_immobile_solid(&curve, 1.0, 1.0, 1.0).unwrap();
        let old = initial(&domain, &curve);
        let fields = vec![[0.0, 0.0, b]; 9];
        let observer = MaterialFlowFeedback::try_new_phase(
            &old,
            &domain,
            flow_input(&flow_ids, &curve),
            &config,
            HeatThenCool::default(),
            &gate,
            1,
            &law,
            |_, cell: &DepositedCell, _| resolve(&card, &law, cell.specific_h_j_kg),
        )
        .unwrap()
        .with_unsteady_viscous_heating(controls())
        .unwrap();
        // Existing callers must still choose this additional mechanical boundary.
        assert!(advance_mhd_enthalpy_with_material_flow(
            &old, 0.1, interval(), current_input(&domain, &drive, &ids),
            motion(&fields), &observer, resistivity,
        ).is_err());
        let (next, accepted, report) = advance_mhd_enthalpy_with_anchored_phase(
            &old, 0.1, interval(), current_input(&domain, &drive, &ids),
            motion(&fields), &observer, resistivity,
        )
        .unwrap();

        // Three independent unit-length/unit-area bulk resistors sit between
        // ideal electrode columns. R=rho_e, I=10/rho_e, JxB_y=-I*B.
        let current = 10.0 / rho_e;
        let joule = 3.0 * 100.0 / rho_e * 0.1;
        near(next.energy_j() - old.energy_j(), joule);
        near(report.electrode_energy_j, joule);
        near(report.joule_energy_j, joule);
        near(report.deposited_joule_j, joule);
        near(report.combined_residual_j, 0.0);
        assert_eq!(report.momentum_trials, 0);
        assert_eq!(accepted.momentum_solves(), 0);
        assert_eq!(report.current_solves, 2);
        assert_eq!(report.thermal.substeps, 1);
        assert_eq!(report.mechanical.substeps, 1);
        assert_eq!(report.mechanical.generated_heat_j, 0.0);
        assert_eq!(report.mechanical.final_kinetic_energy_j, 0.0);
        assert_eq!(report.mechanical.backward_euler_loss_j, 0.0);
        assert_eq!(report.lorentz_work_j, 0.0);
        assert_eq!(report.sampled_lorentz_work_j, 0.0);
        assert_eq!(report.frozen_support_work_j, 0.0);
        let applied = report.last_applied.as_ref().unwrap();
        assert_eq!(applied.frozen_lorentz_reaction_n.len(), 9);
        assert_eq!(report.frozen_lorentz_reaction_impulse_n_s.len(), 9);
        for c in 0..9 {
            let expected = if c % 3 == 1 { [0.0, current * b, 0.0] } else { [0.0; 3] };
            for axis in 0..3 {
                near(applied.frozen_lorentz_reaction_n[c][axis], expected[axis]);
                near(report.endpoint.frozen_lorentz_reaction_n[c][axis], expected[axis]);
                near(report.frozen_lorentz_reaction_impulse_n_s[c][axis], 0.1 * expected[axis]);
                near(
                    applied.solution.lorentz_force_density_n_m3[c][axis] * domain.volume(c)
                        + applied.frozen_lorentz_reaction_n[c][axis],
                    0.0,
                );
            }
            near(
                next.cells()[&ids[c]].specific_h_j_kg,
                -20.0 + if c % 3 == 1 { 10.0 / rho_e } else { 0.0 },
            );
            let sample = &accepted.endpoint()[c];
            assert_eq!(sample.cell, ids[c]);
            assert_eq!(sample.voxel, c);
            assert_eq!(sample.material.material().card_identity(), card.content_hash());
            assert!(sample.transport.mobile_transport().is_none());
            assert!(sample.material.material().property(DYNAMIC_VISCOSITY_PROPERTY).is_none());
            assert!(sample.material.material().property(SPECIFIC_HEAT_CAPACITY_PROPERTY).is_none());
            assert_eq!(accepted.unsteady_flow().unwrap().velocity_m_s[c], [0.0; 3]);
        }
        assert_eq!(accepted.inner().accepted, [0.1]);
        assert_eq!(observer.inner(), &HeatThenCool::default());
        assert_eq!(observer.unsteady_time_s(), Some(0.0));
        assert!(old.cells().values().all(|cell| cell.specific_h_j_kg == -20.0));
        near(next.mass_kg(), old.mass_kg());
    }
}

#[test]
fn shared_heat_releases_anchors_into_material_dependent_flow_but_never_captures_momentum() {
    let domain = VoxelDomain::new(3, 3, 1, 1.0).unwrap();
    let ids = std::array::from_fn::<_, 9, _>(|c| CellId(c as u64));
    let flow_ids = ids.map(Some);
    let gate = CancelGate::new();
    let config = closed();
    let drive = setup(1.0);
    let off = setup(0.0);
    let fields = vec![[0.0, 0.0, 5.0]; 9];
    let zero_fields = vec![[0.0; 3]; 9];
    let mut speeds = Vec::new();
    for mu in [100.0, 200.0] {
        let card = source(mu);
        let curve = chart(&card);
        let law = MushyPhaseLaw::try_new_with_immobile_solid(&curve, 1.0, 1.0, 1.0).unwrap();
        let old = initial(&domain, &curve);
        let observer = MaterialFlowFeedback::try_new_phase(
            &old, &domain, flow_input(&flow_ids, &curve), &config,
            HeatThenCool::melting(), &gate, 100, &law,
            |_, cell: &DepositedCell, _| resolve(&card, &law, cell.specific_h_j_kg),
        )
        .unwrap()
        .with_unsteady_viscous_heating(controls())
        .unwrap();
        let (melted, melting_flow, first) = advance_mhd_enthalpy_with_anchored_phase(
            &old, 0.1, interval(), current_input(&domain, &drive, &ids),
            motion(&fields), &observer, resistivity,
        )
        .unwrap();
        near(melted.energy_j() - old.energy_j(), 1350.3);
        near(first.other_thermal_input_j, 1350.0);
        near(first.joule_energy_j, 0.3);
        near(first.combined_residual_j, 0.0);
        assert_eq!(first.momentum_trials, 0);
        assert_eq!(first.frozen_support_work_j, 0.0);
        for c in 0..9 {
            near(melted.cells()[&ids[c]].specific_h_j_kg,
                130.0 + if c % 3 == 1 { 0.1 } else { 0.0 });
            assert_eq!(melting_flow.endpoint()[c].material.phase_state().liquid_mass_fraction(), 1.0);
            assert!(melting_flow.endpoint()[c].transport.mobile_transport().is_some());
            assert_eq!(first.endpoint.frozen_lorentz_reaction_n[c], [0.0; 3]);
            near(first.frozen_lorentz_reaction_impulse_n_s[c][1],
                if c % 3 == 1 { 0.5 } else { 0.0 });
            assert_eq!(melting_flow.unsteady_flow().unwrap().velocity_m_s[c], [0.0; 3]);
        }

        let (moving, moving_flow, second) = advance_mhd_enthalpy_with_anchored_phase(
            &melted, 0.2, interval(), current_input(&domain, &drive, &ids),
            motion(&fields), &melting_flow, resistivity,
        )
        .unwrap();
        let maximum = moving_flow.unsteady_flow().unwrap().velocity_m_s
            .iter().flatten().map(|v| v.abs()).fold(0.0_f64, f64::max);
        assert!(maximum > 0.0);
        speeds.push(maximum);
        assert!(second.momentum_trials > 1);
        assert_eq!(second.current_solves, second.momentum_trials + 1);
        assert_eq!(moving_flow.momentum_solves(), second.momentum_trials);
        assert!(second.lorentz_work_j > 0.0);
        assert!(second.mechanical.final_kinetic_energy_j > 0.0);
        assert!(second.mechanical.generated_heat_j > 0.0);
        assert!(second.frozen_lorentz_reaction_impulse_n_s.iter().flatten().all(|&p| p == 0.0));
        assert_eq!(second.frozen_support_work_j, 0.0);
        near(second.combined_residual_j, 0.0);
        near(second.electrical_defect_j, 0.0);
        near(second.mechanical_transfer_defect_j, 0.0);
        near(second.native_body_transfer_defect_j, 0.0);
        near(second.other_thermal_input_j, 0.0);
        near(moving.energy_j() - melted.energy_j(),
            second.joule_energy_j + second.mechanical.generated_heat_j);
        for sample in moving_flow.endpoint() {
            assert_eq!(sample.cell, ids[sample.voxel]);
            assert_eq!(sample.material.material().card_identity(), card.content_hash());
            near(sample.transport.mobile_transport().unwrap().dynamic_viscosity_pa_s(), mu);
        }

        // Switching drive and imposed field is an explicit new static-field
        // advance. Cooling cannot silently erase the retained moving faces.
        let retained_native = moving_flow.unsteady_flow().unwrap().clone();
        let retained_mechanical = moving_flow.unsteady_mechanical_step().unwrap().clone();
        let retained_ledger = moving_flow.unsteady_heating_audit();
        let retained_inner = moving_flow.inner().clone();
        let retained_cells = moving.cells().clone();
        let refused = advance_mhd_enthalpy_with_anchored_phase(
            &moving, 0.3, interval(), current_input(&domain, &off, &ids),
            motion(&zero_fields), &moving_flow, resistivity,
        );
        assert!(matches!(refused,
            Err(MaterialFlowMhdError::Flow(MaterialFlowError::Momentum(
                ChtError::InvalidInput { field: "simple.momentum_state", reason }
            ))) if reason.contains("freezing a moving native face")
                && reason.contains("kinetic-energy handoff")
        ));
        assert_eq!(moving_flow.unsteady_flow(), Some(&retained_native));
        assert_eq!(moving_flow.unsteady_mechanical_step(), Some(&retained_mechanical));
        assert_eq!(moving_flow.unsteady_heating_audit(), retained_ledger);
        assert_eq!(moving_flow.inner(), &retained_inner);
        assert_eq!(moving_flow.unsteady_time_s(), Some(0.2));
        assert_eq!(moving.cells(), &retained_cells);
        assert_eq!(moving.time_s(), 0.2);
        assert!(observer.applied().is_none());
        assert_eq!(observer.unsteady_time_s(), Some(0.0));
    }
    // Only source viscosity changed; the same electrical and momentum owners
    // respond without any material-specific branch in the physics code.
    assert!(speeds[1] < speeds[0]);
}

#[test]
fn final_current_budget_refusal_rolls_back_phase_change_and_support_impulse() {
    let card = source(100.0);
    let curve = chart(&card);
    let law = MushyPhaseLaw::try_new_with_immobile_solid(&curve, 1.0, 1.0, 1.0).unwrap();
    let domain = VoxelDomain::new(3, 3, 1, 1.0).unwrap();
    let old = initial(&domain, &curve);
    let ids = std::array::from_fn::<_, 9, _>(|c| CellId(c as u64));
    let flow_ids = ids.map(Some);
    let gate = CancelGate::new();
    let config = closed();
    let drive = setup(1.0);
    let fields = vec![[0.0, 0.0, 5.0]; 9];
    let observer = MaterialFlowFeedback::try_new_phase(
        &old, &domain, flow_input(&flow_ids, &curve), &config,
        HeatThenCool::melting(), &gate, 1, &law,
        |_, cell: &DepositedCell, _| resolve(&card, &law, cell.specific_h_j_kg),
    )
    .unwrap()
    .with_unsteady_viscous_heating(controls())
    .unwrap();
    let mut limited = current_input(&domain, &drive, &ids);
    limited.max_current_solves = 1;
    assert!(matches!(
        advance_mhd_enthalpy_with_anchored_phase(
            &old, 0.1, interval(), limited, motion(&fields), &observer, resistivity,
        ),
        Err(MaterialFlowMhdError::Electrical(JouleEnthalpyError::CurrentBudget))
    ));
    assert_eq!(old.time_s(), 0.0);
    assert!(old.cells().values().all(|cell| cell.specific_h_j_kg == -20.0));
    assert_eq!(observer.unsteady_time_s(), Some(0.0));
    assert_eq!(observer.momentum_solves(), 0);
    assert!(observer.applied().is_none());
    assert!(observer.unsteady_mechanical_step().is_none());
    assert_eq!(observer.unsteady_heating_audit().unwrap().substeps, 0);
    assert_eq!(observer.inner(), &HeatThenCool::melting());
    assert!(observer.endpoint().iter().all(|s| s.transport.mobile_transport().is_none()));
    let (next, accepted, report) = advance_mhd_enthalpy_with_anchored_phase(
        &old, 0.1, interval(), current_input(&domain, &drive, &ids),
        motion(&fields), &observer, resistivity,
    )
    .unwrap();
    near(next.energy_j() - old.energy_j(), 1350.3);
    assert_eq!(report.current_solves, 2);
    assert_eq!(report.momentum_trials, 0);
    assert_eq!(report.frozen_support_work_j, 0.0);
    assert_eq!(accepted.inner().accepted, [0.1]);
    for c in 0..9 {
        assert!(accepted.endpoint()[c].transport.mobile_transport().is_some());
        assert_eq!(report.endpoint.frozen_lorentz_reaction_n[c], [0.0; 3]);
        near(report.frozen_lorentz_reaction_impulse_n_s[c][1],
            if c % 3 == 1 { 0.5 } else { 0.0 });
    }
    // A successful retry from the same unmodified input publishes the same
    // thermal/native and support accounts, with no leaked failed-trial impulse.
    let (repeated, repeated_flow, repeated_report) = advance_mhd_enthalpy_with_anchored_phase(
        &old, 0.1, interval(), current_input(&domain, &drive, &ids),
        motion(&fields), &observer, resistivity,
    )
    .unwrap();
    assert_eq!(next.cells(), repeated.cells());
    assert_eq!(accepted.endpoint(), repeated_flow.endpoint());
    assert_eq!(report, repeated_report);
}

#[derive(Clone)]
struct MeltRing;
impl DepositionObserver for MeltRing {
    type Error = DepositionError;
    fn event(&mut self, _: &Event, _: &DepositionSimulation<'_>) -> Result<(), Self::Error> {
        Ok(())
    }
    fn sources(
        &mut self,
        thermal: &DepositionSimulation<'_>,
        _: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<BTreeMap<CellId, f64>, Self::Error> {
        Ok(thermal.cells().keys().map(|&id| {
            let watts = if id != CellId(4) && thermal.time_s() < 0.1 { 1500.0 } else { 0.0 };
            (id, watts)
        }).collect())
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

#[test]
fn mixed_phase_mhd_retains_solid_reaction_and_mobile_native_force_work() {
    let card = source(100.0);
    let curve = chart(&card);
    let law = MushyPhaseLaw::try_new_with_immobile_solid(&curve, 1.0, 1.0, 1.0).unwrap();
    let domain = VoxelDomain::new(3, 3, 1, 1.0).unwrap();
    let old = initial(&domain, &curve);
    let ids = std::array::from_fn::<_, 9, _>(|c| CellId(c as u64));
    let flow_ids = ids.map(Some);
    let gate = CancelGate::new();
    let config = closed();
    let drive = setup(1.0);
    let fields = vec![[0.0, 0.0, 5.0]; 9];
    let observer = MaterialFlowFeedback::try_new_phase(
        &old, &domain, flow_input(&flow_ids, &curve), &config,
        MeltRing, &gate, 100, &law,
        |_, cell: &DepositedCell, _| resolve(&card, &law, cell.specific_h_j_kg),
    )
    .unwrap()
    .with_unsteady_viscous_heating(controls())
    .unwrap();
    let (mixed, phase, _) = advance_mhd_enthalpy_with_anchored_phase(
        &old, 0.1, interval(), current_input(&domain, &drive, &ids),
        motion(&fields), &observer, resistivity,
    )
    .unwrap();
    for c in 0..9 {
        assert_eq!(phase.endpoint()[c].transport.mobile_transport().is_none(), c == 4);
    }
    // Only the mobile ring receives this actuator. Retain the native momentum
    // history and the conducting central solid; its force cannot drive fluid.
    let forces = vec![
        [-1.0, 1.0, 0.0], [-1.0, 0.0, 0.0], [-1.0, -1.0, 0.0],
        [0.0, 1.0, 0.0], [0.0; 3], [0.0, -1.0, 0.0],
        [1.0, 1.0, 0.0], [1.0, 0.0, 0.0], [1.0, -1.0, 0.0],
    ];
    let driven = phase.with_unsteady_body_force_heating(forces, controls()).unwrap();
    let (next, accepted, report) = advance_mhd_enthalpy_with_anchored_phase(
        &mixed, 0.2, interval(), current_input(&domain, &drive, &ids),
        motion(&fields), &driven, resistivity,
    )
    .unwrap();
    let native = accepted.unsteady_flow().unwrap();
    assert!(native.velocity_m_s.iter().flatten().any(|v| v.abs() > 0.0));
    assert_eq!(native.velocity_m_s[4], [0.0; 3]);
    assert_eq!(native.field.flux_x(1, 1, 0), 0.0);
    assert_eq!(native.field.flux_x(2, 1, 0), 0.0);
    assert_eq!(native.field.flux_y(1, 1, 0), 0.0);
    assert_eq!(native.field.flux_y(1, 2, 0), 0.0);
    assert!(report.other_body_work_j > 0.0);
    assert!(report.mechanical.final_kinetic_energy_j > 0.0);
    assert!(report.mechanical.generated_heat_j > 0.0);
    assert_eq!(report.thermal.substeps, 1);
    assert_eq!(report.frozen_support_work_j, 0.0);
    near(report.combined_residual_j, 0.0);
    near(report.native_body_transfer_defect_j, 0.0);
    near(report.mechanical_transfer_defect_j, 0.0);
    near(report.electrical_defect_j, 0.0);
    near(next.energy_j() - mixed.energy_j(),
        report.joule_energy_j + report.mechanical.generated_heat_j);
    let applied = report.last_applied.as_ref().unwrap();
    let dt = report.last_duration_s.unwrap();
    assert!(applied.frozen_lorentz_reaction_n[4].iter().any(|f| f.abs() > 0.0));
    assert!(applied.solution.electric.joule_w[4] > 0.0);
    for c in 0..9 {
        assert_eq!(accepted.endpoint()[c].transport.mobile_transport().is_none(), c == 4);
        for axis in 0..3 {
            let reaction = applied.frozen_lorentz_reaction_n[c][axis];
            near(reaction, if c == 4 {
                -applied.solution.lorentz_force_density_n_m3[c][axis] * domain.volume(c)
            } else { 0.0 });
            near(report.frozen_lorentz_reaction_impulse_n_s[c][axis], reaction * dt);
        }
    }
    assert!(accepted.endpoint()[4].material.material()
        .property(DYNAMIC_VISCOSITY_PROPERTY).is_none());
    assert_eq!(driven.unsteady_time_s(), Some(0.1));
    assert!(driven.unsteady_flow().unwrap().velocity_m_s.iter().flatten().all(|&v| v == 0.0));
}
