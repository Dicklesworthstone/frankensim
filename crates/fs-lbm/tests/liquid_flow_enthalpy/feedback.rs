//! Data-only liquid temperature feedback through the existing physical owners.
use super::*;
use fs_lbm::conjugate::FlowResistance;
use fs_lbm::conjugate::enthalpy::advection::feedback::MaterialFlowError;

fn feedback_card(scale: f64) -> MaterialCard {
    feedback_card_with_conductivity(
        scale,
        scalar(
            THERMAL_CONDUCTIVITY_PROPERTY,
            THERMAL_CONDUCTIVITY_DIMS,
            0.01,
        ),
    )
}

fn feedback_card_with_conductivity(scale: f64, conductivity: PropertyClaim) -> MaterialCard {
    let mut claims = ClaimSet::new();
    for (name, dims, value) in [
        (DENSITY_PROPERTY, Density::DIMS, 1.0),
        (
            SPECIFIC_HEAT_CAPACITY_PROPERTY,
            SPECIFIC_HEAT_CAPACITY_DIMS,
            2.0,
        ),
    ] {
        claims.insert_claim(scalar(name, dims, value)).unwrap();
    }
    claims.insert_claim(conductivity).unwrap();
    let mut viscosity = scalar(DYNAMIC_VISCOSITY_PROPERTY, DYNAMIC_VISCOSITY_DIMS, scale);
    viscosity.value = PropertyValue::Curve {
        abscissa: "T".into(),
        abscissa_dims: Temperature::DIMS,
        knots: vec![(300.0, scale), (340.0, 0.2 * scale), (810.0, 0.1 * scale)],
        dims: DYNAMIC_VISCOSITY_DIMS,
    };
    viscosity.interpolation = InterpolationPolicy::LinearInside;
    claims.insert_claim(viscosity).unwrap();
    card(claims)
}

fn porous_config() -> SimpleConfig {
    let mut faces = [FvBoundary::Symmetry; 6];
    faces[0] = FvBoundary::Outlet;
    faces[1] = FvBoundary::Outlet;
    let mut config = SimpleConfig::new(faces);
    config.tolerance = 1e-13;
    config.pressure_tolerance = 1e-12;
    config.resistances.push(FlowResistance::Volume {
        lo: [0, 0, 0],
        hi: [4, 1, 1],
        permeability_m2: [1.0; 3],
        inertial_per_m: [0.0; 3],
    });
    config.internal_fans.push(InternalFan {
        patch: FacePatch {
            axis: 0,
            index: 2,
            lo: [0, 0, 0],
            hi: [4, 1, 1],
        },
        blows_positive: true,
        curve: FanCurve::new(&[(0.0, 4.0), (100.0, 4.0)]).unwrap(),
    });
    config
}

fn viscosity(scale: f64, h: f64) -> f64 {
    let t = 300.0 + h / 2.0;
    if t <= 340.0 {
        scale * (1.0 - 0.8 * (t - 300.0) / 40.0)
    } else {
        scale * (0.2 - 0.1 * (t - 340.0) / 470.0)
    }
}

#[test]
fn source_temperature_changes_actual_flow_each_substep_and_closes_enthalpy() {
    let domain = VoxelDomain::new(4, 1, 1, 1.0).unwrap();
    let config_flow = porous_config();
    let gate = CancelGate::new();
    let ids: Vec<_> = (0..4).map(|c| Some(CellId(c))).collect();
    let mut outputs = Vec::new();
    for scale in [1.0, 2.0] {
        let card = feedback_card(scale);
        let liquid = liquid(&card);
        let caloric = liquid_caloric(&card);
        let old = initial(&domain, 1.0, &caloric);
        let observer = input(&ids, &liquid, &caloric)
            .feedback(&old, &domain, &config_flow, Heater::default(), &gate, 100)
            .unwrap();
        let (next, accepted, audit) = old
            .advance_observed(&observer, 1.0, config(), || Ok(()))
            .unwrap();
        let mut h = [0.0; 4];
        let mut exported = 0.0;
        let mut mass_through = 0.0;
        let mut last_q = 0.0;
        for &dt in &accepted.inner().durations {
            // Independent Darcy series-resistance + donor-cell energy balance.
            // The native MAC solver owns the actual flow used by the test.
            let q = 4.0 / h.iter().map(|&h| viscosity(scale, h)).sum::<f64>();
            exported += dt * q * h[3];
            mass_through += dt * q;
            let old_h = h;
            for c in 0..4 {
                let upstream = if c == 0 { 0.0 } else { old_h[c - 1] };
                let heater = if c == 1 { 100.0 } else { 0.0 };
                h[c] += dt * (heater + q * (upstream - old_h[c]));
            }
            last_q = q;
        }
        assert!(accepted.momentum_solves() > 2);
        assert_eq!(accepted.momentum_solves(), accepted.inner().durations.len());
        assert_eq!(accepted.audit().substeps, accepted.inner().durations.len());
        let applied = accepted.applied().unwrap();
        near(applied.flow.field.flux_x(2, 0, 0), last_q);
        assert!(last_q > 1.0 / scale);
        assert!(applied.time_s < next.time_s());
        for (c, expected) in h.iter().enumerate() {
            near(next.cells()[&CellId(c as u64)].specific_h_j_kg, *expected);
        }
        for sample in accepted.endpoint() {
            near(sample.specific_enthalpy_j_kg, h[sample.voxel]);
            near(sample.temperature_k, 300.0 + h[sample.voxel] / 2.0);
            near(
                sample.transport.dynamic_viscosity_pa_s(),
                viscosity(scale, h[sample.voxel]),
            );
            let receipt = &sample
                .material
                .material()
                .property(DYNAMIC_VISCOSITY_PROPERTY)
                .unwrap()
                .answer()
                .receipt;
            assert_eq!(
                receipt.selected,
                liquid
                    .reference_state()
                    .material()
                    .property(DYNAMIC_VISCOSITY_PROPERTY)
                    .unwrap()
                    .answer()
                    .receipt
                    .selected
            );
            assert!(receipt.query_point.contains(&("pressure".into(), 100000.0)));
            assert!(receipt.query_point.contains(&("lot".into(), 7.0)));
        }
        near(accepted.audit().inlet_mass_kg, mass_through);
        near(accepted.audit().outlet_mass_kg, mass_through);
        near(accepted.audit().external_enthalpy_j, -exported);
        near(next.energy_j() - old.energy_j(), 100.0 - exported);
        near(
            audit.external_heat_j,
            accepted.inner().supplied_j - exported,
        );
        near(next.mass_kg(), old.mass_kg());
        assert_eq!(observer.momentum_solves(), 0);
        assert!(observer.applied().is_none());
        assert_eq!(observer.inner(), &Heater::default());
        outputs.push((last_q, h[1]));
    }
    assert!(outputs[1].0 < outputs[0].0);
    assert!(outputs[1].1 > outputs[0].1);
}

#[test]
fn fixed_viscosity_feedback_matches_existing_frozen_owner() {
    let card = card(claims(1.0, 1.0, false));
    let liquid = liquid(&card);
    let caloric = liquid_caloric(&card);
    let domain = VoxelDomain::new(4, 1, 1, 1.0).unwrap();
    let config_flow = porous_config();
    let gate = CancelGate::new();
    let ids: Vec<_> = (0..4).map(|c| Some(CellId(c))).collect();
    let old = initial(&domain, 1.0, &caloric);
    let live = input(&ids, &liquid, &caloric)
        .feedback(&old, &domain, &config_flow, Heater::default(), &gate, 100)
        .unwrap();
    let fluid = MaterialChtFluid::try_new(liquid.reference_state()).unwrap();
    let flow = fluid.solve_flow(&domain, &config_flow, &gate).unwrap();
    let frozen = flow
        .advect_liquid_enthalpy(
            &old,
            input(&ids, &liquid, &caloric),
            Heater::default(),
            &gate,
        )
        .unwrap();
    let (a, a_observer, _) = old
        .advance_observed(&live, 1.0, config(), || Ok(()))
        .unwrap();
    let (b, b_observer, _) = old
        .advance_observed(&frozen, 1.0, config(), || Ok(()))
        .unwrap();
    for &id in ids.iter().flatten() {
        near(
            a.cells()[&id].specific_h_j_kg,
            b.cells()[&id].specific_h_j_kg,
        );
    }
    near(
        a_observer.audit().external_enthalpy_j,
        b_observer.advection().audit().external_enthalpy_j,
    );
}

#[test]
fn source_support_budget_inner_history_and_late_cancellation_are_atomic() {
    let card = card(claims(1.0, 1.0, true));
    let liquid = liquid(&card);
    let caloric = liquid_caloric(&card);
    let domain = VoxelDomain::new(4, 1, 1, 1.0).unwrap();
    let config_flow = porous_config();
    let gate = CancelGate::new();
    let ids: Vec<_> = (0..4).map(|c| Some(CellId(c))).collect();
    let old = initial(&domain, 1.0, &caloric);
    let observer = input(&ids, &liquid, &caloric)
        .feedback(&old, &domain, &config_flow, Heater::default(), &gate, 100)
        .unwrap();
    assert!(matches!(
        old.advance_observed(&observer, 1.0, config(), || Ok(())),
        Err(MaterialFlowError::Material { .. })
    ));
    let limited = input(&ids, &liquid, &caloric)
        .feedback(&old, &domain, &config_flow, Heater::default(), &gate, 1)
        .unwrap();
    assert!(matches!(
        old.advance_observed(&limited, 1.0, config(), || Ok(())),
        Err(MaterialFlowError::Thermal(DepositionError::Budget(
            "live material momentum solves"
        )))
    ));
    let history = Heater {
        reject_step: Some(1),
        ..Heater::default()
    };
    let rejecting = input(&ids, &liquid, &caloric)
        .feedback(&old, &domain, &config_flow, history.clone(), &gate, 100)
        .unwrap();
    assert!(matches!(
        old.advance_observed(&rejecting, 1.0, config(), || Ok(())),
        Err(MaterialFlowError::Inner(DepositionError::InvalidInput(_)))
    ));
    let mut checks = 0;
    let result = old.advance_observed(&observer, 1.0, config(), || {
        checks += 1;
        if checks == 30 {
            Err(DepositionError::Cancelled)
        } else {
            Ok(())
        }
    });
    assert!(matches!(
        result,
        Err(MaterialFlowError::Thermal(DepositionError::Cancelled))
            | Err(MaterialFlowError::Inner(DepositionError::Cancelled))
    ));
    assert_eq!(observer.momentum_solves(), 0);
    assert!(observer.applied().is_none());
    assert_eq!(observer.audit().substeps, 0);
    assert_eq!(observer.inner(), &Heater::default());
    assert_eq!(rejecting.inner(), &history);
    near(old.time_s(), 0.0);
    near(old.energy_j(), 0.0);
    near(old.mass_kg(), 4.0);
}

#[test]
fn live_binding_keeps_original_caloric_context_and_borrowed_law() {
    let card = feedback_card(1.0);
    let liquid = liquid(&card);
    let caloric = liquid_caloric(&card);
    let pressure = source_caloric(
        &card,
        200000.0,
        SPECIFIC_HEAT_CAPACITY_PROPERTY,
        DENSITY_PROPERTY,
        SolidLiquidPhase::Liquid,
    );
    let equal_curve = caloric.clone();
    let domain = VoxelDomain::new(4, 1, 1, 1.0).unwrap();
    let config_flow = porous_config();
    let gate = CancelGate::new();
    let ids: Vec<_> = (0..4).map(|c| Some(CellId(c))).collect();
    let old = initial(&domain, 1.0, &caloric);
    for other in [&pressure, &equal_curve] {
        assert!(
            input(&ids, &liquid, other)
                .feedback(&old, &domain, &config_flow, (), &gate, 100,)
                .is_err()
        );
    }
}

fn conducting_card(k_scale: f64) -> MaterialCard {
    let mut k = scalar(
        THERMAL_CONDUCTIVITY_PROPERTY,
        THERMAL_CONDUCTIVITY_DIMS,
        k_scale,
    );
    k.value = PropertyValue::Curve {
        abscissa: "T".into(),
        abscissa_dims: Temperature::DIMS,
        knots: vec![(300.0, k_scale), (810.0, 6.1 * k_scale)],
        dims: THERMAL_CONDUCTIVITY_DIMS,
    };
    k.interpolation = InterpolationPolicy::LinearInside;
    feedback_card_with_conductivity(1.0, k)
}

fn conductivity(k_scale: f64, h: f64) -> f64 {
    k_scale * (1.0 + h / 200.0)
}

// Same four-cell flow fixture, now with physical full-face thermal contacts.
fn conducting_initial<'a>(
    domain: &VoxelDomain,
    caloric: &'a ResolvedHeatCapacityCurve,
    h: [f64; 4],
) -> DepositionSimulation<'a> {
    let seed = initial(domain, 1.0, caloric);
    let mut events: Vec<_> = seed
        .cells()
        .values()
        .cloned()
        .map(|mut cell| {
            let c = cell.id.0 as usize;
            cell.specific_h_j_kg = h[c];
            cell.geometry.surface_area_m2 =
                2.0 * (0..3).map(|axis| domain.face_area(c, axis)).sum::<f64>();
            TimedEvent {
                time_s: 0.0,
                event: Event::Birth(cell),
            }
        })
        .collect();
    for c in 0..3 {
        events.push(TimedEvent {
            time_s: 0.0,
            event: Event::Contact(ThermalLink {
                id: 100 + c as u64,
                a: CellId(c as u64),
                b: CellId(c as u64 + 1),
                area_m2: domain.face_area(c, 0),
                conductance_w_k: 0.001,
            }),
        });
    }
    let max_events = events.len();
    DepositionSimulation::new(
        events,
        &[(MAT, caloric.curve() as &dyn ThermalLaw)],
        DepositionBudget {
            max_cells: 4,
            max_links: 64,
            max_events,
            max_substeps: 10000,
        },
    )
    .unwrap()
    .advance_to(0.0, config(), || Ok(()))
    .unwrap()
    .0
}

#[test]
fn data_only_conductivity_changes_native_heat_feedback_and_combined_cfl() {
    let domain = VoxelDomain::new(4, 1, 1, 1.0).unwrap();
    let config_flow = porous_config();
    let gate = CancelGate::new();
    let ids: Vec<_> = (0..4).map(|c| Some(CellId(c))).collect();
    let mut outputs = Vec::new();
    for k_scale in [0.25, 40.0] {
        let card = conducting_card(k_scale);
        let liquid = liquid(&card);
        let caloric = liquid_caloric(&card);
        let old = conducting_initial(&domain, &caloric, [0.0; 4]);
        let observer = input(&ids, &liquid, &caloric)
            .feedback(&old, &domain, &config_flow, Heater::default(), &gate, 1000)
            .unwrap()
            .with_fluid_conduction(&old, &[100, 101, 102])
            .unwrap();
        let (next, accepted, audit) = old
            .advance_observed(&observer, 1.0, config(), || Ok(()))
            .unwrap();
        let mut h = [0.0; 4];
        let mut exported = 0.0;
        let mut last_q = 0.0;
        for &dt in &accepted.inner().durations {
            // Independent series Darcy resistance, series thermal resistance,
            // and donor-cell energy balance at the ACTUAL native step lengths.
            let q = 4.0 / h.iter().map(|&h| viscosity(1.0, h)).sum::<f64>();
            let k = h.map(|h| conductivity(k_scale, h));
            let g = [0, 1, 2].map(|c| 2.0 / (1.0 / k[c] + 1.0 / k[c + 1]));
            let old_h = h;
            exported += dt * q * old_h[3];
            for c in 0..4 {
                let upstream = if c == 0 { 0.0 } else { old_h[c - 1] };
                let heater = if c == 1 { 100.0 } else { 0.0 };
                let left = if c > 0 {
                    g[c - 1] * (old_h[c - 1] - old_h[c]) / 2.0
                } else {
                    0.0
                };
                let right = if c < 3 {
                    g[c] * (old_h[c + 1] - old_h[c]) / 2.0
                } else {
                    0.0
                };
                h[c] += dt * (heater + q * (upstream - old_h[c]) + left + right);
            }
            last_q = q;
        }
        for (c, expected) in h.iter().enumerate() {
            near(next.cells()[&CellId(c as u64)].specific_h_j_kg, *expected);
        }
        let applied = accepted.applied().unwrap();
        near(applied.flow.field.flux_x(2, 0, 0), last_q);
        for sample in &applied.samples {
            near(
                sample.transport.thermal_conductivity_w_m_k(),
                conductivity(k_scale, sample.specific_enthalpy_j_kg),
            );
            let receipt = &sample
                .material
                .material()
                .property(THERMAL_CONDUCTIVITY_PROPERTY)
                .unwrap()
                .answer()
                .receipt;
            assert_eq!(
                receipt.selected,
                liquid
                    .reference_state()
                    .material()
                    .property(THERMAL_CONDUCTIVITY_PROPERTY)
                    .unwrap()
                    .answer()
                    .receipt
                    .selected
            );
            assert!(receipt.query_point.contains(&("pressure".into(), 100000.0)));
            assert!(receipt.query_point.contains(&("lot".into(), 7.0)));
        }
        near(accepted.audit().external_enthalpy_j, -exported);
        near(audit.external_heat_j, 100.0 - exported);
        near(next.energy_j() - old.energy_j(), audit.external_heat_j);
        near(next.mass_kg(), old.mass_kg());
        assert_eq!(
            next.links(),
            old.links(),
            "step coefficients must not replace authored links"
        );
        assert_eq!(accepted.momentum_solves(), accepted.inner().durations.len());
        assert!(accepted.audit().maximum_courant <= 0.8 + 1e-12);
        assert_eq!(observer.inner(), &Heater::default());
        assert!(observer.applied().is_none());
        outputs.push((last_q, h[1], accepted.inner().durations.len()));
    }
    // Higher k spreads heater energy to the outlet sooner. That changes both
    // stored h and the temperature-dependent viscosity driving the next flow.
    assert!(outputs[1].0 < outputs[0].0);
    assert!(outputs[1].1 < outputs[0].1);
    assert!(
        outputs[1].2 > outputs[0].2,
        "conduction must participate in native CFL"
    );
}

#[test]
fn graded_fluid_faces_use_both_half_cell_resistances_and_actual_cell_mass() {
    let domain = VoxelDomain::graded([vec![0.5, 1.0, 1.5, 1.0], vec![2.0], vec![3.0]]).unwrap();
    let mut config_flow = porous_config();
    config_flow.internal_fans.clear();
    let gate = CancelGate::new();
    let ids: Vec<_> = (0..4).map(|c| Some(CellId(c))).collect();
    let card = conducting_card(10.0);
    let liquid = liquid(&card);
    let caloric = liquid_caloric(&card);
    let old = conducting_initial(&domain, &caloric, [40.0, 0.0, 0.0, 0.0]);
    let observer = input(&ids, &liquid, &caloric)
        .feedback(&old, &domain, &config_flow, (), &gate, 4)
        .unwrap()
        .with_fluid_conduction(&old, &[100, 101, 102])
        .unwrap();
    let dt = 0.001;
    let (next, accepted, audit) = old
        .advance_observed(
            &observer,
            dt,
            StepConfig {
                max_dt_s: dt,
                energy_tolerance_j: 1e-9,
            },
            || Ok(()),
        )
        .unwrap();
    // A=6 m2, half-widths 0.25 and 0.5 m, k(T)=12 and 10 W/(m K).
    let transferred = dt * 6.0 / (0.25 / 12.0 + 0.5 / 10.0) * 20.0;
    near(
        next.cells()[&CellId(0)].specific_h_j_kg,
        40.0 - transferred / 3.0,
    );
    near(next.cells()[&CellId(1)].specific_h_j_kg, transferred / 6.0);
    near(next.cells()[&CellId(2)].specific_h_j_kg, 0.0);
    near(next.cells()[&CellId(3)].specific_h_j_kg, 0.0);
    near(next.energy_j(), old.energy_j());
    near(audit.external_heat_j, 0.0);
    near(accepted.audit().external_enthalpy_j, 0.0);
    near(accepted.applied().unwrap().flow.field.flux_x(2, 0, 0), 0.0);
    assert_eq!(audit.substeps, 1);
    assert_eq!(next.links(), old.links());
}

#[derive(Clone)]
struct ConductanceOwner;
impl DepositionObserver for ConductanceOwner {
    type Error = DepositionError;
    fn conductances(
        &mut self,
        _: &DepositionSimulation<'_>,
        _: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<BTreeMap<u64, f64>, Self::Error> {
        Ok(BTreeMap::from([(100, 0.0)]))
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

#[test]
fn conductive_faces_refuse_partial_parallel_competing_and_stale_ownership() {
    let domain = VoxelDomain::new(4, 1, 1, 1.0).unwrap();
    let config_flow = porous_config();
    let gate = CancelGate::new();
    let ids: Vec<_> = (0..4).map(|c| Some(CellId(c))).collect();
    let card = conducting_card(1.0);
    let liquid = liquid(&card);
    let caloric = liquid_caloric(&card);
    let old = conducting_initial(&domain, &caloric, [0.0; 4]);
    let observer = input(&ids, &liquid, &caloric)
        .feedback(&old, &domain, &config_flow, (), &gate, 100)
        .unwrap()
        .with_fluid_conduction(&old, &[100, 101, 102])
        .unwrap();
    assert!(
        observer
            .clone()
            .with_fluid_conduction(&old, &[100, 100])
            .is_err()
    );

    // Nominal area alone cannot admit an asymmetric native contact patch.
    let patch = old.contact_patch(100).unwrap();
    let mut partial = patch;
    partial.covered_area_m2[1] *= 0.5;
    let partial = old
        .replace_contact_patches(
            &[ContactPatchUpdate {
                id: 100,
                expected: Some(patch),
                replacement: Some(partial),
            }],
            || Ok(()),
        )
        .unwrap();
    assert!(matches!(
        partial.advance_observed(&observer, 0.1, config(), || Ok(())),
        Err(MaterialFlowError::Thermal(DepositionError::InvalidInput(
            "fluid conduction must retain both full-face coverages"
        )))
    ));

    // An unselected parallel link would still add heat to the same physical face.
    let mut extra = patch.link;
    extra.id = 200;
    extra.area_m2 = 0.25;
    let parallel = old
        .replace_contacts(
            &[ContactUpdate {
                id: 200,
                expected: None,
                replacement: Some(extra),
            }],
            || Ok(()),
        )
        .unwrap();
    assert!(matches!(
        parallel.advance_observed(&observer, 0.1, config(), || Ok(())),
        Err(MaterialFlowError::Thermal(DepositionError::InvalidInput(
            "parallel thermal links duplicate an owned fluid face"
        )))
    ));

    let competing = input(&ids, &liquid, &caloric)
        .feedback(&old, &domain, &config_flow, ConductanceOwner, &gate, 100)
        .unwrap()
        .with_fluid_conduction(&old, &[100, 101, 102])
        .unwrap();
    assert!(matches!(
        old.advance_observed(&competing, 0.1, config(), || Ok(())),
        Err(MaterialFlowError::Thermal(DepositionError::InvalidInput(
            "multiple owners supplied one fluid-face conductance"
        )))
    ));

    // Same clock and geometry do not make a different enthalpy field this epoch.
    let other = conducting_initial(&domain, &caloric, [1.0, 0.0, 0.0, 0.0]);
    let mut pending = observer.clone();
    pending.begin().unwrap();
    pending.sources(&old, &mut || Ok(())).unwrap();
    assert!(pending.conductances(&old, &mut || Ok(())).is_ok());
    assert!(matches!(
        pending.conductances(&other, &mut || Ok(())),
        Err(MaterialFlowError::Thermal(DepositionError::InvalidInput(
            "fluid conduction must use the actual applied enthalpy state"
        )))
    ));
    assert_eq!(observer.momentum_solves(), 0);
    assert!(observer.applied().is_none());
    assert_eq!(competing.momentum_solves(), 0);
    assert_eq!(old.contact_patch(100), Some(patch));
    near(old.time_s(), 0.0);
    near(old.energy_j(), 0.0);
}
