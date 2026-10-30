//! Synthetic interface sources against the independent 3x2 MAC/current balance.
use super::*;
use fs_lbm::conjugate::ElectricalContact;
use fs_lbm::conjugate::enthalpy::advection::feedback::mhd::interfaces::advance_mhd_enthalpy_with_material_interfaces;
use fs_lbm::conjugate::enthalpy::transport::interfaces::{
    ElectricalInterfaceResponse, ElectricalInterfaceState,
};
use fs_lbm::conjugate::material_coupling::interfaces::MaterialContact;
use fs_matdb::{InterfaceSystemCard, SurfaceSpec, SystemContext};
use fs_material::electrical_interface::{
    ELECTRICAL_CONTACT_RESISTANCE_DIMS, ELECTRICAL_CONTACT_RESISTANCE_PROPERTY,
    resolve_electrical_interface_state,
};

fn surface(material: &ConductiveLiquid, voxel: usize) -> SurfaceSpec {
    SurfaceSpec {
        material: material.conductor.material().material().clone(),
        texture_frame: if voxel % 3 == 1 {
            "specimen side"
        } else {
            "electrode side"
        }
        .into(),
    }
}
fn interface_source(
    material: &ConductiveLiquid,
    resistance: f64,
    slope: f64,
    upper: f64,
) -> InterfaceSystemCard {
    let dims = ELECTRICAL_CONTACT_RESISTANCE_DIMS;
    let mut claim = scalar(ELECTRICAL_CONTACT_RESISTANCE_PROPERTY, dims, resistance);
    claim.validity = claim.validity.with_quantity("T", absolute(), 300.0, upper);
    if slope != 0.0 {
        claim.value = PropertyValue::Curve {
            abscissa: "T".into(),
            abscissa_dims: Temperature::DIMS,
            knots: vec![
                (300.0, resistance),
                (upper, resistance + slope * (upper - 300.0)),
            ],
            dims,
        };
        claim.interpolation = InterpolationPolicy::LinearInside;
    }
    let mut claims = ClaimSet::new();
    claims.insert_claim(claim).unwrap();
    InterfaceSystemCard::assemble(
        surface(material, 1),
        surface(material, 0),
        SystemContext {
            medium: "synthetic".into(),
            third_body: None,
            environment: "numerical fixture".into(),
            history: "explicit switching law".into(),
        },
        claims,
        vec![],
    )
    .unwrap()
}

// Deliberately not Clone: original constitutive evidence must survive without
// being cloned on each rejected momentum trial.
#[derive(Debug)]
struct Receipt {
    contact: MaterialContact,
    bulk_queries: [usize; 2],
    temperature_k: [f64; 2],
}
fn response(
    card: &InterfaceSystemCard,
    template: &ElectricalContact,
    sides: [&FlowMaterialSample<ConductiveLiquid>; 2],
    open_above_k: Option<f64>,
) -> Result<ElectricalInterfaceState<Receipt>, JouleEnthalpyError> {
    // The declared fixture law uses the first/specimen temperature, not an
    // inferred average. Resolve full source support even when explicitly open.
    let source = resolve_electrical_interface_state(
        card,
        &point(sides[0].temperature_k, 100000.0),
        MaterialPropertySelection::SingleClaimOnly,
    )
    .map_err(|e| ChtError::InvalidInput {
        field: "mhd.interface_source",
        reason: e.to_string(),
    })?;
    let actual = sides.map(|s| surface(&s.material, s.voxel));
    let contact = MaterialContact::try_new(source, template.cells, [&actual[0], &actual[1]], 0.25)?;
    let face = *contact.face();
    let receipt = Receipt {
        contact,
        bulk_queries: sides.map(|s| s.material.query),
        temperature_k: sides.map(|s| s.temperature_k),
    };
    if open_above_k.is_some_and(|limit| sides[0].temperature_k > limit) {
        Ok(ElectricalInterfaceState::Open(receipt))
    } else {
        Ok(ElectricalInterfaceResponse::try_new(
            receipt,
            face.resistance_ohm_m2,
            face.heat_fraction_to_first,
        )?
        .into())
    }
}
fn contacts() -> [ElectricalContact; 2] {
    // Reverse spatial order deliberately: heat fraction and source surfaces
    // must follow the authored roles, not sorted voxel indices.
    [[1, 0], [4, 3]].map(|cells| ElectricalContact {
        cells,
        resistance_ohm_m2: 999.0,
        heat_fraction_to_first: 1.0,
    })
}

#[test]
fn source_contact_changes_native_current_force_heat_and_endpoint_feedback() {
    let domain = VoxelDomain::new(3, 2, 1, 1.0).unwrap();
    let ids: Vec<_> = (0..6).map(CellId).collect();
    let flow_ids: Vec<_> = ids.iter().copied().map(Some).collect();
    let card = source(1.0, false);
    let caloric = liquid_caloric(&card);
    let old = initial(&domain, 1.0, &caloric);
    let config = closed_box();
    let gate = CancelGate::new();
    let drive = setup(-1.0);
    let fields = vec![[0.0, 0.0, 10.0]; 6];
    let faces = contacts();
    let mut currents = Vec::new();
    for (rc, slope) in [(0.5, 0.0), (2.0, 0.0), (1.0, 0.5)] {
        let interface = interface_source(&resolve(&card, 300.0, 0).unwrap(), rc, slope, 810.0);
        let queries = Cell::new(0);
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
        .with_unsteady_viscous_heating(controls())
        .unwrap();
        let mut input = electrical_input(&domain, &drive, &ids);
        input.contacts = &faces;
        let (next, accepted, report) = advance_mhd_enthalpy_with_material_interfaces(
            &old,
            DT,
            thermal_step(),
            input,
            motion(&fields),
            &observer,
            resistivity,
            |_, template, sides| response(&interface, template, sides, None),
        )
        .unwrap();
        // Each row has one unit bulk resistance plus the source interface.
        // Interfaces add no motional EMF or mechanical work of their own.
        let total_r = 1.0 + rc;
        let q = amplitude(0.0, DT, -1.0, 10.0, total_r, 0.0);
        let current = (-1.0 + 10.0 * q / 2.0) / total_r;
        let mhd = &report.mhd;
        let applied = mhd.last_applied.as_ref().unwrap();
        near(accepted.unsteady_flow().unwrap().field.flux_y(1, 1, 0), q);
        near(mhd.electrode_energy_j, -2.0 * current * DT);
        near(mhd.joule_energy_j, 2.0 * total_r * current * current * DT);
        near(mhd.lorentz_work_j, -10.0 * current * q * DT);
        near(mhd.mechanical.generated_heat_j, 15.5 * MU * q * q * DT);
        for c in [0, 3] {
            near(
                applied.solution.electric.joule_w[c],
                0.75 * rc * current * current,
            );
        }
        for c in [1, 4] {
            near(
                applied.solution.electric.joule_w[c],
                (1.0 + 0.25 * rc) * current * current,
            );
        }
        for c in [2, 5] {
            near(applied.solution.electric.joule_w[c], 0.0);
        }
        near(
            next.energy_j() - old.energy_j(),
            mhd.joule_energy_j + mhd.mechanical.generated_heat_j,
        );
        near(mhd.combined_residual_j, 0.0);
        near(mhd.heat_transfer_defect_j, 0.0);
        near(mhd.mechanical_transfer_defect_j, 0.0);
        assert_eq!(accepted.inner().accepted, [DT]);
        assert_eq!(queries.get(), 18);
        let last = report.last_applied.as_ref().unwrap();
        near(last.time_s, applied.time_s);
        near(report.endpoint.time_s, mhd.endpoint.time_s);
        for (index, face) in faces.iter().enumerate() {
            let a = last.contacts[index].response.material();
            let e = report.endpoint.contacts[index].response.material();
            assert_eq!(
                a.bulk_queries,
                face.cells
                    .map(|c| accepted.applied().unwrap().samples[c].material.query)
            );
            assert_eq!(
                e.bulk_queries,
                face.cells.map(|c| accepted.endpoint()[c].material.query)
            );
            assert_eq!(
                e.temperature_k,
                face.cells.map(|c| next.temperature_k(ids[c]).unwrap())
            );
            assert_eq!(
                e.contact.source().resolved().card_identity(),
                interface.content_hash()
            );
            assert_eq!(
                e.contact.source().resolved().state_coordinate("T"),
                Some(e.temperature_k[0])
            );
            near(
                e.contact.source().resistance_ohm_m2(),
                rc + slope * (e.temperature_k[0] - 300.0),
            );
            near(
                last.contacts[index].applied().unwrap().resistance_ohm_m2,
                rc,
            );
        }
        let endpoint_r = report.endpoint.contacts[0]
            .applied()
            .unwrap()
            .resistance_ohm_m2;
        near(
            mhd.endpoint.solution.electric.electrodes[0].current_a,
            2.0 * (-1.0 + 10.0 * q / 2.0) / (1.0 + endpoint_r),
        );
        if slope > 0.0 {
            assert!(endpoint_r > rc);
            let (_, continued, feedback) = advance_mhd_enthalpy_with_material_interfaces(
                &next,
                2.0 * DT,
                thermal_step(),
                input,
                motion(&fields),
                &accepted,
                resistivity,
                |_, template, sides| response(&interface, template, sides, None),
            )
            .unwrap();
            near(
                continued.unsteady_flow().unwrap().field.flux_y(1, 1, 0),
                amplitude(q, DT, -1.0, 10.0, 1.0 + endpoint_r, 0.0),
            );
            near(
                feedback.last_applied.as_ref().unwrap().contacts[0]
                    .applied()
                    .unwrap()
                    .resistance_ohm_m2,
                endpoint_r,
            );
            near(feedback.mhd.combined_residual_j, 0.0);
        }
        currents.push(current.abs());
    }
    assert!(currents[0] > currents[1]);
    assert_eq!(old.time_s(), 0.0);
    near(old.energy_j(), 0.0);
}

#[test]
fn opening_stops_current_and_force_but_retains_native_flow_and_thermal_history() {
    let domain = VoxelDomain::new(3, 2, 1, 1.0).unwrap();
    let ids: Vec<_> = (0..6).map(CellId).collect();
    let flow_ids: Vec<_> = ids.iter().copied().map(Some).collect();
    let card = source(1.0, false);
    let caloric = liquid_caloric(&card);
    let old = initial(&domain, 1.0, &caloric);
    let config = closed_box();
    let gate = CancelGate::new();
    let drive = setup(-1.0);
    let fields = vec![[0.0, 0.0, 10.0]; 6];
    let faces = contacts();
    let material = resolve(&card, 300.0, 0).unwrap();
    let interface = interface_source(&material, 1.0, 0.0, 810.0);
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
    let mut input = electrical_input(&domain, &drive, &ids);
    input.contacts = &faces;
    // Explicit synthetic opening criterion; missing source support is never Open.
    let law = |_: usize,
               template: &ElectricalContact,
               sides: [&FlowMaterialSample<ConductiveLiquid>; 2]| {
        response(&interface, template, sides, Some(300.00001))
    };
    let (warm, moving, first) = advance_mhd_enthalpy_with_material_interfaces(
        &old,
        DT,
        thermal_step(),
        input,
        motion(&fields),
        &observer,
        resistivity,
        law,
    )
    .unwrap();
    assert!(
        first
            .last_applied
            .as_ref()
            .unwrap()
            .contacts
            .iter()
            .all(|s| !s.response.is_open())
    );
    assert!(first.endpoint.contacts.iter().all(|s| s.response.is_open()));
    near(first.mhd.endpoint.solution.electric.total_joule_w, 0.0);
    let q1 = moving.unsteady_flow().unwrap().field.flux_y(1, 1, 0);
    let (coasting, free, second) = advance_mhd_enthalpy_with_material_interfaces(
        &warm,
        2.0 * DT,
        thermal_step(),
        input,
        motion(&fields),
        &moving,
        resistivity,
        law,
    )
    .unwrap();
    let q2 = amplitude(q1, DT, 0.0, 0.0, 1.0, 0.0);
    near(free.unsteady_flow().unwrap().field.flux_y(1, 1, 0), q2);
    assert!(q2 > 0.0 && q2 < q1);
    near(second.mhd.joule_energy_j, 0.0);
    near(second.mhd.electrode_energy_j, 0.0);
    near(second.mhd.lorentz_work_j, 0.0);
    near(
        coasting.energy_j() - warm.energy_j(),
        second.mhd.mechanical.generated_heat_j,
    );
    near(second.mhd.combined_residual_j, 0.0);
    assert_eq!(coasting.cells().len(), warm.cells().len());
    assert_eq!(free.inner().accepted, [DT, DT]);
    assert!(
        second
            .last_applied
            .as_ref()
            .unwrap()
            .contacts
            .iter()
            .all(|s| s.response.is_open())
    );
    let (_, restarted, third) = advance_mhd_enthalpy_with_material_interfaces(
        &coasting,
        3.0 * DT,
        thermal_step(),
        input,
        motion(&fields),
        &free,
        resistivity,
        |_, template, sides| response(&interface, template, sides, None),
    )
    .unwrap();
    near(
        restarted.unsteady_flow().unwrap().field.flux_y(1, 1, 0),
        amplitude(q2, DT, -1.0, 10.0, 2.0, 0.0),
    );
    assert!(third.mhd.joule_energy_j > 0.0);
    assert!(
        third
            .last_applied
            .as_ref()
            .unwrap()
            .contacts
            .iter()
            .all(|s| !s.response.is_open())
    );
    // Source support fails only AFTER positive heating. The candidate flow and
    // contact history must not escape, and an opening cannot conceal the refusal.
    let narrow = interface_source(&material, 1.0, 0.0, 300.000001);
    let failed = advance_mhd_enthalpy_with_material_interfaces(
        &old,
        DT,
        thermal_step(),
        input,
        motion(&fields),
        &observer,
        resistivity,
        |_, template, sides| response(&narrow, template, sides, Some(300.00001)),
    );
    assert!(matches!(
        failed,
        Err(MaterialFlowMhdError::Electrical(
            JouleEnthalpyError::Electric(ChtError::InvalidInput {
                field: "mhd.interface_source",
                ..
            })
        ))
    ));
    assert_eq!(observer.unsteady_time_s(), Some(0.0));
    assert_eq!(observer.inner(), &Intervals::default());
    assert_eq!(old.time_s(), 0.0);
    near(old.energy_j(), 0.0);
}
