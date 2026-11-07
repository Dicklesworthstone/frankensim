//! Synthetic nonlinear damping/enthalpy checks; no measured material claims.
use fs_couple::membrane_thermal::{
    ElementThermalCell, MembraneThermalConfig, MembraneThermalError, advance_candidate,
};
use fs_evidence::ValidityDomain;
use fs_exec::CancelGate;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::plane_stress_viscosity::{
    MaterialPlane, PLANE_STRESS_VISCOSITY_PROPERTIES as KEYS, PlaneStressViscosity,
};
use fs_material::sheet::{
    AREAL_DENSITY, MEMBRANE_DIMS, ResolvedSheet, THICKNESS, thermal::SheetCaloric,
};
use fs_material::state_point::{
    MaterialPropertySelection, SPECIFIC_HEAT_CAPACITY_DIMS, SPECIFIC_HEAT_CAPACITY_PROPERTY as CP,
    ScalarAdmissibility, ScalarPropertyRequirement, resolve_material_state_point,
};
use fs_qty::{Dims, DynViscosity, QuantitySpec,
    semantic::{MoistureBasis, QuantityKind, SemanticType, ValueForm}};
use fs_solid::membrane::{
    MembraneBudget, MembraneDomain, MembranePatch,
    dynamics::{
        MembraneMotion,
        viscosity::{DampedMembraneMotion, MembraneDampingBudget},
    },
    viscosity::MembraneViscosity,
};
use fs_solid::mesh2::Mesh2;
use fs_thermal::enthalpy::{
    Boundary, CellGeometry, CellId, DepositedCell, DepositionBudget, DepositionError,
    DepositionObserver, DepositionSimulation, Event, MaterialId, StepConfig, ThermalLaw,
    TimedEvent,
};
use std::collections::BTreeMap;

fn fixture(
    variable_cp: bool,
    density: &[[f64; 3]],
    revision: u32,
) -> (MembranePatch, Vec<PlaneStressViscosity>, SheetCaloric) {
    fixture_at_moisture(variable_cp, density, revision, None)
}

fn fixture_at_moisture(
    variable_cp: bool,
    density: &[[f64; 3]],
    revision: u32,
    moisture: Option<f64>,
) -> (MembranePatch, Vec<PlaneStressViscosity>, SheetCaloric) {
    let mut properties = vec![
        (AREAL_DENSITY, 6.0, Dims([-2, 1, 0, 0, 0, 0])),
        (THICKNESS, 0.002, Dims([1, 0, 0, 0, 0, 0])),
        ("warp_stiffness", 1200.0, MEMBRANE_DIMS),
        ("weft_stiffness", 600.0, MEMBRANE_DIMS),
        ("poisson_warp_weft", 0.0, Dims::NONE),
        ("shear_stiffness", 250.0, MEMBRANE_DIMS),
    ];
    properties.extend(
        KEYS.into_iter()
            .zip([20_000.0, 20_000.0, 20_000.0, 0.0, 0.0, 0.0])
            .map(|(key, value)| (key, value, DynViscosity::DIMS)),
    );
    let mut claims = ClaimSet::new();
    let mut requirements = Vec::new();
    let provenance = Provenance {
        source: "synthetic sheet damping and calorics; no specimen authority".into(),
        license: "test-only".into(),
        artifact: None,
    };
    for (key, value, dims) in properties {
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(key, dims),
                value: PropertyValue::Scalar { value, dims },
                interpolation: InterpolationPolicy::ConstantWithinValidity,
                validity: ValidityDomain::unconstrained()
                    .with("T", 299.0, 301.0)
                    .with("omega", 0.0, 1000.0),
                uncertainty: UncertaintyModel::Unstated,
                observations: vec![],
                provenance: provenance.clone(),
            })
            .unwrap();
        requirements.push(
            ScalarPropertyRequirement::try_new(key, dims, ScalarAdmissibility::Finite).unwrap(),
        );
    }
    let moisture_quantity = QuantitySpec::semantic(SemanticType::new(
        QuantityKind::Moisture(MoistureBasis::DryMass), ValueForm::Static,
    ));
    let cp_validity = ValidityDomain::unconstrained().with("T", 299.0, 301.0);
    let cp_validity = if moisture.is_some() {
        cp_validity.with_quantity("water_content", moisture_quantity, 0.0, 0.2)
    } else {
        cp_validity.with("omega", 0.0, 1000.0)
    };
    claims
        .insert_claim(PropertyClaim {
            key: PropertyKey::new(CP, SPECIFIC_HEAT_CAPACITY_DIMS),
            value: if moisture.is_some() {
                PropertyValue::Surface2D {
                    axes: ["T".into(), "water_content".into()],
                    axis_dims: [fs_qty::Temperature::DIMS, Dims::NONE],
                    samples: vec![
                        ([299.0, 0.0], 0.019), ([300.0, 0.0], 0.020), ([301.0, 0.0], 0.021),
                        ([299.0, 0.2], 0.038), ([300.0, 0.2], 0.040), ([301.0, 0.2], 0.042),
                    ],
                    triangles: vec![[0, 1, 3], [1, 4, 3], [1, 2, 4], [2, 5, 4]],
                    dims: SPECIFIC_HEAT_CAPACITY_DIMS,
                }
            } else if variable_cp {
                PropertyValue::Curve {
                    abscissa: "T".into(),
                    abscissa_dims: Dims([0, 0, 0, 1, 0, 0]),
                    knots: vec![(299.0, 0.019), (300.0, 0.020), (301.0, 0.021)],
                    dims: SPECIFIC_HEAT_CAPACITY_DIMS,
                }
            } else {
                PropertyValue::Scalar {
                    value: 1.0,
                    dims: SPECIFIC_HEAT_CAPACITY_DIMS,
                }
            },
            interpolation: if variable_cp || moisture.is_some() {
                InterpolationPolicy::LinearInside
            } else {
                InterpolationPolicy::ConstantWithinValidity
            },
            validity: cp_validity,
            uncertainty: UncertaintyModel::Unstated,
            observations: vec![],
            provenance,
        })
        .unwrap();
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic membrane".into(),
            phase: "solid".into(),
            process: "damping heat test".into(),
            revision,
        },
        claims,
        vec![],
    )
    .unwrap();
    let mut query = QueryPoint::new()
        .with("T", 300.0)
        .unwrap()
        .with("omega", 0.0)
        .unwrap();
    if let Some(moisture) = moisture {
        query = query.with_quantity("water_content", moisture_quantity, moisture).unwrap();
    }
    let source = resolve_material_state_point(
        &card,
        &query,
        &requirements,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    let law =
        PlaneStressViscosity::from_material(source.clone(), MaterialPlane::Axes12, KEYS).unwrap();
    let sheet = ResolvedSheet::try_from_material(source).unwrap();
    let cp = resolve_material_state_point(
        &card,
        &query,
        &[ScalarPropertyRequirement::try_new(
            CP,
            SPECIFIC_HEAT_CAPACITY_DIMS,
            ScalarAdmissibility::StrictlyPositive,
        )
        .unwrap()],
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap()
    .property(CP)
    .unwrap()
    .answer()
    .receipt
    .selected;
    let caloric =
        SheetCaloric::try_new(sheet.clone(), &card, cp, [299.0, 301.0], 0.0, 1e-9).unwrap();
    let count = density.len();
    let mesh = Mesh2 {
        nodes: (0..count)
            .flat_map(|i| {
                let x = 2.0 * i as f64;
                [[x, 0.0], [x + 1.0, 0.0], [x, 1.0]]
            })
            .collect(),
        elems: (0..count)
            .map(|i| vec![3 * i, 3 * i + 1, 3 * i + 2])
            .collect(),
        patches: vec![],
    };
    let patch = MembranePatch::try_new(
        mesh,
        vec![sheet; count],
        vec![0.0; count],
        MembraneDomain {
            minimum_area_ratio: 0.01,
            maximum_abs_strain: 0.5,
        },
        MembraneBudget {
            max_nodes: 3 * count,
            max_triangles: count,
        },
        || false,
    )
    .unwrap()
    .with_reference_areal_density(density.to_vec(), || false)
    .unwrap();
    (patch, vec![law; count], caloric)
}

fn motion<'m, 'p>(
    viscosity: &'m MembraneViscosity<'p>,
    states: &[[f64; 2]],
) -> DampedMembraneMotion<'m, 'p> {
    let mut u = vec![[0.0; 3]; 3 * states.len()];
    let mut v = u.clone();
    let free: Vec<_> = states
        .iter()
        .enumerate()
        .map(|(i, &[z, velocity])| {
            u[3 * i + 1] = [0.08, 0.0, z];
            v[3 * i + 1][2] = velocity;
            9 * i + 5
        })
        .collect();
    let fixed: Vec<_> = (0..3 * u.len())
        .filter(|d| !free.contains(d))
        .map(|d| (d, u[d / 3][d % 3]))
        .collect();
    let native = MembraneMotion::new(viscosity.patch(), &fixed, &u, &v, 0.0, || false).unwrap();
    DampedMembraneMotion::new(
        viscosity,
        native,
        MembraneDampingBudget {
            max_source_rows: 6 * states.len(),
            max_mass_solves_per_step: 0,
            max_step_s: 0.05,
            max_energy_defect_j: 1.0,
        },
        || false,
    )
    .unwrap()
}

fn loads(values: &[f64]) -> Vec<[f64; 3]> {
    let mut result = vec![[0.0; 3]; 3 * values.len()];
    for (i, &load) in values.iter().enumerate() {
        result[3 * i + 1][2] = load;
    }
    result
}

fn boundary(source_w: f64) -> Boundary {
    Boundary {
        ambient_k: 300.0,
        htc_w_m2_k: 0.0,
        reservoir: None,
        source_w,
    }
}

fn thermal<'a>(
    caloric: &'a SheetCaloric,
    inventory: &[(u64, f64, f64)],
    extra: Vec<TimedEvent>,
) -> DepositionSimulation<'a> {
    let initial_h = caloric
        .state_at_temperature(300.0)
        .unwrap()
        .specific_enthalpy_j_kg();
    let mut schedule: Vec<_> = inventory
        .iter()
        .enumerate()
        .map(|(i, &(id, mass, volume))| TimedEvent {
            time_s: 0.0,
            event: Event::Birth(DepositedCell {
                id: CellId(id),
                segment_id: id,
                body_id: 7,
                material_id: MaterialId(0),
                mass_kg: mass,
                specific_h_j_kg: initial_h,
                geometry: CellGeometry {
                    start_m: [2.0 * i as f64, 0.0, 0.0],
                    end_m: [2.0 * i as f64 + 1.0, 0.0, 0.0],
                    volume_m3: volume,
                    surface_area_m2: 1.0,
                },
                boundary: boundary(0.0),
            }),
        })
        .collect();
    schedule.extend(extra);
    DepositionSimulation::new(
        schedule,
        &[(MaterialId(0), caloric.curve() as &dyn ThermalLaw)],
        DepositionBudget {
            max_cells: 4,
            max_links: 4,
            max_events: 16,
            max_substeps: 64,
        },
    )
    .unwrap()
    .advance_to(
        0.0,
        StepConfig {
            max_dt_s: 0.04,
            energy_tolerance_j: 1e-10,
        },
        || Ok(()),
    )
    .unwrap()
    .0
}

fn cells<'a>(caloric: &'a SheetCaloric, ids: &[u64]) -> Vec<ElementThermalCell<'a>> {
    ids.iter()
        .map(|&id| ElementThermalCell {
            cell: CellId(id),
            caloric,
            temperature_interval_k: [299.0, 301.0],
        })
        .collect()
}

fn config() -> MembraneThermalConfig {
    MembraneThermalConfig {
        thermal: StepConfig {
            max_dt_s: 0.04,
            energy_tolerance_j: 1e-10,
        },
        energy_tolerance_j: 0.1,
    }
}

fn near(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance,
        "{a:e} != {b:e}, tolerance {tolerance:e}"
    );
}

struct ScalarStep {
    z: f64,
    v: f64,
    heat: f64,
    work: f64,
    defect: f64,
}

// Independent authored-sheet reduction, not a mechanical or caloric API:
// U=300(.0832+z²/2)²; f=49.92*z+300*z³; C(z)=20*z².
fn scalar_step(z0: f64, v0: f64, mass: f64, dt: f64, load: f64) -> ScalarStep {
    let energy = |z: f64| 300.0 * (0.0832 + 0.5 * z * z).powi(2);
    let force = |z: f64| 49.92 * z + 300.0 * z.powi(3);
    let decay = |z: f64, v: f64| {
        let exponent = -20.0 * z * z * dt / mass;
        (
            v * (0.5 * exponent).exp(),
            -0.5 * mass * v * v * exponent.exp_m1(),
        )
    };
    let (va, qa) = decay(z0, v0);
    let half = va + 0.5 * dt * (load - force(z0)) / mass;
    let z = z0 + dt * half;
    let vb = half + 0.5 * dt * (load - force(z)) / mass;
    let (v, qb) = decay(z, vb);
    let work = load * (z - z0);
    ScalarStep {
        z,
        v,
        heat: qa + qb,
        work,
        defect: energy(z) - energy(z0) + 0.5 * mass * (vb - va) * (vb + va) - work,
    }
}

#[test]
fn variable_heat_capacity_uses_original_enthalpy_inversion_and_full_facet_mass() {
    let (p, laws, caloric) = fixture(true, &[[6.0; 3]], 0);
    let viscosity = MembraneViscosity::new(&p, laws, (0.0, 1000.0), || false).unwrap();
    let original = motion(&viscosity, &[[0.2, -0.3]]);
    let heat = thermal(&caloric, &[(30, 3.0, 0.001)], vec![]);
    let expected = scalar_step(0.2, -0.3, 1.0, 0.04, 0.7);
    let next = advance_candidate(
        &original,
        &heat,
        &(),
        0.04,
        &loads(&[0.7]),
        &cells(&caloric, &[30]),
        config(),
        &CancelGate::new_clock_free(),
    )
    .unwrap();
    near(
        next.motion().motion().displacement_m()[1][2],
        expected.z,
        2e-12,
    );
    near(
        next.motion().motion().velocity_m_s()[1][2],
        expected.v,
        2e-12,
    );
    // The free mechanical lumped mass is 1 kg, but all 3 kg receive this facet's heat.
    let delta_h = expected.heat / 3.0;
    let delta_t = 2.0 * delta_h / (0.02 + (0.02_f64.powi(2) + 0.002 * delta_h).sqrt());
    let temperature = next.thermal().temperature_k(CellId(30)).unwrap();
    // Cp(T)=.02+.001*(T-300); its integral is .02*dT+.0005*dT².
    // 1e-9 J/kg chart interpolation bounds T error by 1e-9/.019 K.
    near(temperature, 300.0 + delta_t, 6e-8);
    assert!((temperature - (300.0 + delta_h / 0.02)).abs() > 3e-4);
    near(
        next.thermal().cells()[&CellId(30)].specific_h_j_kg
            - heat.cells()[&CellId(30)].specific_h_j_kg,
        delta_h,
        2e-12,
    );
    let report = next.report();
    near(report.element_heat_j[0], expected.heat, 2e-12);
    near(report.material_heat_j, expected.heat, 2e-12);
    near(report.deposited_heat_j, expected.heat, 2e-12);
    near(report.thermal.stored_energy_change_j, expected.heat, 2e-12);
    near(report.mechanical_work_j, expected.work, 2e-12);
    near(report.mechanical_numerical_defect_j, expected.defect, 2e-12);
    near(report.mechanical_residual_j, expected.defect, 2e-12);
    near(report.combined_residual_j, expected.defect, 2e-12);
    assert!(expected.defect.abs() > 1e-3);
    near(
        report.absolute_mechanical_defect_j,
        expected.defect.abs(),
        2e-12,
    );
    near(report.heat_handoff_residual_j, 0.0, 2e-12);
    near(report.external_thermal_heat_j, 0.0, 2e-12);
    near(
        report.support_impulse_n_s[2],
        expected.v + 0.3 - 0.04 * 0.7,
        2e-12,
    );
    assert_eq!(original.motion().time_s(), 0.0);
    assert_eq!(heat.time_s(), 0.0);
}

#[test]
fn g1_held_moisture_changes_native_membrane_heat_response_without_changing_mechanics() {
    let expected = scalar_step(0.2, -0.3, 1.0, 0.04, 0.7);
    let mut temperatures = Vec::new();
    let mut card_identity = None;
    for (moisture, capacity_scale) in [(0.0, 1.0), (0.2, 2.0)] {
        // Same complete material card, mechanics, geometry, reference mass and
        // forcing. Only the held source coordinate selects a different Cp(T).
        let (p, laws, caloric) = fixture_at_moisture(true, &[[6.0; 3]], 0, Some(moisture));
        let identity = caloric.sheet().material().card_identity();
        if let Some(original) = card_identity { assert_eq!(identity, original); }
        card_identity = Some(identity);
        for sample in caloric.heat_capacity_samples() {
            assert_eq!(sample.query().axes()["water_content"], moisture);
            assert_eq!(sample.card_identity(), identity);
        }
        let viscosity = MembraneViscosity::new(&p, laws, (0.0, 1000.0), || false).unwrap();
        let original = motion(&viscosity, &[[0.2, -0.3]]);
        let heat = thermal(&caloric, &[(30, 3.0, 0.001)], vec![]);
        let next = advance_candidate(
            &original, &heat, &(), 0.04, &loads(&[0.7]), &cells(&caloric, &[30]),
            config(), &CancelGate::new_clock_free(),
        ).unwrap();
        near(next.motion().motion().displacement_m()[1][2], expected.z, 2e-12);
        near(next.motion().motion().velocity_m_s()[1][2], expected.v, 2e-12);
        near(next.report().material_heat_j, expected.heat, 2e-12);
        near(next.report().thermal.stored_energy_change_j, expected.heat, 2e-12);
        near(next.report().heat_handoff_residual_j, 0.0, 2e-12);
        let delta_h = expected.heat / 3.0;
        let normalized_h = delta_h / capacity_scale;
        let delta_t = 2.0 * normalized_h
            / (0.02 + (0.02_f64.powi(2) + 0.002 * normalized_h).sqrt());
        let temperature = next.thermal().temperature_k(CellId(30)).unwrap();
        near(temperature, 300.0 + delta_t, 6e-8);
        let cell = &next.thermal().cells()[&CellId(30)];
        near(cell.specific_h_j_kg - heat.cells()[&CellId(30)].specific_h_j_kg,
            delta_h, 2e-12);
        assert_eq!(cell.mass_kg, 3.0, "held conditioning does not simulate water uptake");
        assert_eq!(cell.geometry, heat.cells()[&CellId(30)].geometry);
        temperatures.push(temperature);
    }
    assert!(temperatures[0] > temperatures[1] && temperatures[1] > 300.0);
}

#[test]
fn facet_routes_preserve_original_p1_inventory_and_do_not_sort_cell_ids() {
    let (p, laws, caloric) = fixture(false, &[[6.0; 3], [6.0, 9.0, 12.0]], 0);
    assert_eq!(
        p.reference_triangle_areas_m2().collect::<Vec<_>>(),
        [0.5, 0.5]
    );
    near(p.reference_mass_kg(), 7.5, 1e-14);
    let viscosity = MembraneViscosity::new(&p, laws, (0.0, 1000.0), || false).unwrap();
    let original = motion(&viscosity, &[[0.2, -0.3], [0.1, 0.7]]);
    // Original facet order deliberately opposes BTreeMap/CellId order.
    let heat = thermal(&caloric, &[(30, 3.0, 0.001), (10, 4.5, 0.001)], vec![]);
    let bindings = cells(&caloric, &[30, 10]);
    let next = advance_candidate(
        &original,
        &heat,
        &(),
        0.04,
        &loads(&[0.7, -0.4]),
        &bindings,
        config(),
        &CancelGate::new_clock_free(),
    )
    .unwrap();
    let expected = [
        scalar_step(0.2, -0.3, 1.0, 0.04, 0.7),
        // Integral rho*N_1 dA=.5*(6+2*9+12)/12=1.5, total mass=4.5.
        scalar_step(0.1, 0.7, 1.5, 0.04, -0.4),
    ];
    for (e, (id, mass)) in [(30, 3.0), (10, 4.5)].into_iter().enumerate() {
        near(next.report().element_heat_j[e], expected[e].heat, 2e-12);
        near(
            next.motion().motion().displacement_m()[3 * e + 1][2],
            expected[e].z,
            2e-12,
        );
        near(
            next.thermal().temperature_k(CellId(id)).unwrap(),
            300.0 + expected[e].heat / mass,
            2e-11,
        );
        let before = &heat.cells()[&CellId(id)];
        let after = &next.thermal().cells()[&CellId(id)];
        assert_eq!(after.mass_kg, before.mass_kg);
        assert_eq!(after.geometry, before.geometry);
        near(
            mass * (after.specific_h_j_kg - before.specific_h_j_kg),
            expected[e].heat,
            2e-12,
        );
    }
    near(
        next.report().material_heat_j,
        expected.iter().map(|s| s.heat).sum(),
        2e-12,
    );
    let duplicated = cells(&caloric, &[30, 30]);
    assert!(
        advance_candidate(
            &original,
            &heat,
            &(),
            0.04,
            &loads(&[0.7, -0.4]),
            &duplicated,
            config(),
            &CancelGate::new_clock_free()
        )
        .is_err()
    );
}

#[derive(Clone, Default)]
struct Observer<'g> {
    intervals: Vec<(f64, f64, f64)>,
    source_calls: usize,
    events: usize,
    steps: usize,
    fail_after: Option<f64>,
    cancel_after: Option<(&'g CancelGate, f64)>,
}
impl DepositionObserver for Observer<'_> {
    type Error = DepositionError;
    fn sources(
        &mut self,
        _: &DepositionSimulation<'_>,
        _: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<BTreeMap<CellId, f64>, Self::Error> {
        self.source_calls += 1;
        Ok(BTreeMap::from([(CellId(30), 0.001)]))
    }
    fn limit_timestep_s(
        &mut self,
        _: &DepositionSimulation<'_>,
        proposed: f64,
        _: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<f64, Self::Error> {
        Ok(proposed.min(0.015625))
    }
    fn requires_step_interval(&self) -> bool {
        true
    }
    fn step_sources(
        &mut self,
        _: &DepositionSimulation<'_>,
        start: f64,
        dt: f64,
        end: f64,
        _: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<BTreeMap<CellId, f64>, Self::Error> {
        self.intervals.push((start, dt, end));
        Ok(BTreeMap::from([(CellId(30), 0.002)]))
    }
    fn event(&mut self, _: &Event, _: &DepositionSimulation<'_>) -> Result<(), Self::Error> {
        self.events += 1;
        Ok(())
    }
    fn step(
        &mut self,
        _: f64,
        _: f64,
        _: &BTreeMap<CellId, f64>,
        thermal: &DepositionSimulation<'_>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), Self::Error> {
        self.steps += 1;
        if self.fail_after.is_some_and(|time| thermal.time_s() > time) {
            return Err(DepositionError::InvalidInput(
                "synthetic late observer refusal",
            ));
        }
        if let Some((gate, time)) = self.cancel_after {
            if thermal.time_s() > time {
                gate.request();
                checkpoint()?;
            }
        }
        Ok(())
    }
}

#[test]
fn thermal_event_cuts_and_observer_intervals_drive_the_same_continuing_motion() {
    let (p, laws, caloric) = fixture(false, &[[6.0; 3]], 0);
    let viscosity = MembraneViscosity::new(&p, laws, (0.0, 1000.0), || false).unwrap();
    let mut original = motion(&viscosity, &[[0.2, -0.3]]);
    let switch = 0.0234375;
    let mut heat = thermal(
        &caloric,
        &[(30, 3.0, 0.001)],
        vec![TimedEvent {
            time_s: switch,
            event: Event::Boundary(CellId(30), boundary(0.004)),
        }],
    );
    let mut observer = Observer::default();
    let mut z = 0.2;
    let mut v = -0.3;
    let mut total_heat = 0.0;
    let mut thermal_h = heat.cells()[&CellId(30)].specific_h_j_kg;
    for target in [0.0625, 0.09375] {
        let before = observer.intervals.len();
        let start = heat.time_s();
        let next = advance_candidate(
            &original,
            &heat,
            &observer,
            target,
            &loads(&[0.7]),
            &cells(&caloric, &[30]),
            config(),
            &CancelGate::new_clock_free(),
        )
        .unwrap();
        let intervals = &next.observer().intervals[before..];
        assert_eq!(intervals.len(), if before == 0 { 5 } else { 2 });
        let mut interval_heat = 0.0;
        let mut defects = 0.0;
        for &(time, dt, end) in intervals {
            assert_eq!(time + dt, end);
            let expected = scalar_step(z, v, 1.0, dt, 0.7);
            let outside = 0.003 + if time >= switch { 0.004 } else { 0.0 };
            thermal_h += (expected.heat + outside * dt) / 3.0;
            interval_heat += expected.heat;
            defects += expected.defect;
            z = expected.z;
            v = expected.v;
        }
        total_heat += interval_heat;
        near(next.motion().motion().displacement_m()[1][2], z, 2e-12);
        near(next.motion().motion().velocity_m_s()[1][2], v, 2e-12);
        near(next.motion().cumulative_dissipated_j(), total_heat, 2e-12);
        near(
            next.thermal().cells()[&CellId(30)].specific_h_j_kg,
            thermal_h,
            2e-12,
        );
        near(next.report().element_heat_j[0], interval_heat, 2e-12);
        near(next.report().mechanical_numerical_defect_j, defects, 2e-12);
        let external = 0.003 * (target - start) + 0.004 * (target - start.max(switch));
        near(next.report().external_thermal_heat_j, external, 2e-12);
        assert_eq!(next.motion().motion().time_s(), target);
        assert_eq!(next.thermal().time_s(), target);
        assert_eq!(next.observer().steps, next.observer().intervals.len());
        assert_eq!(next.observer().source_calls, next.observer().steps);
        assert_eq!(next.observer().events, 1);
        assert!(
            next.observer()
                .intervals
                .iter()
                .any(|&(time, _, _)| time == switch)
        );
        let (m, t, o, _) = next.into_parts();
        original = m;
        heat = t;
        observer = o;
    }
    assert!(heat.temperature_k(CellId(30)).unwrap() > 300.0);
}

#[derive(Debug, PartialEq)]
struct MechanicalSnapshot {
    displacement: Vec<[f64; 3]>,
    velocity: Vec<[f64; 3]>,
    time: f64,
    energy: [f64; 2],
    element_heat: Vec<f64>,
    cumulative_heat: f64,
}
fn snapshot(m: &DampedMembraneMotion<'_, '_>) -> MechanicalSnapshot {
    MechanicalSnapshot {
        displacement: m.motion().displacement_m().to_vec(),
        velocity: m.motion().velocity_m_s().to_vec(),
        time: m.motion().time_s(),
        energy: [m.motion().strain_energy_j(), m.motion().kinetic_energy_j()],
        element_heat: m.element_dissipated_j().to_vec(),
        cumulative_heat: m.cumulative_dissipated_j(),
    }
}

#[test]
fn late_temperature_observer_and_cancellation_refusals_preserve_both_histories() {
    let (p, laws, caloric) = fixture(false, &[[6.0; 3]], 0);
    let viscosity = MembraneViscosity::new(&p, laws, (0.0, 1000.0), || false).unwrap();
    let original = motion(&viscosity, &[[0.2, -0.3]]);
    let heat = thermal(&caloric, &[(30, 3.0, 0.001)], vec![]);
    let saved = snapshot(&original);
    let saved_cells = heat.cells().clone();
    let mut limited = cells(&caloric, &[30]);
    limited[0].temperature_interval_k = [299.0, 300.0 + 1e-10];
    assert!(
        advance_candidate(
            &original,
            &heat,
            &(),
            0.04,
            &loads(&[0.7]),
            &limited,
            config(),
            &CancelGate::new_clock_free()
        )
        .is_err()
    );
    let failing = Observer {
        fail_after: Some(0.02),
        ..Observer::default()
    };
    assert!(
        advance_candidate(
            &original,
            &heat,
            &failing,
            0.04,
            &loads(&[0.7]),
            &cells(&caloric, &[30]),
            config(),
            &CancelGate::new_clock_free()
        )
        .is_err()
    );
    assert!(failing.intervals.is_empty());
    assert_eq!(failing.steps, 0);
    let gate = CancelGate::new_clock_free();
    let cancelling = Observer {
        cancel_after: Some((&gate, 0.02)),
        ..Observer::default()
    };
    assert!(
        advance_candidate(
            &original,
            &heat,
            &cancelling,
            0.04,
            &loads(&[0.7]),
            &cells(&caloric, &[30]),
            config(),
            &gate
        )
        .is_err()
    );
    assert!(gate.is_requested());
    assert!(cancelling.intervals.is_empty());
    assert_eq!(snapshot(&original), saved);
    assert_eq!(heat.cells(), &saved_cells);
    assert_eq!(heat.time_s(), 0.0);
    let retry = advance_candidate(
        &original,
        &heat,
        &(),
        0.04,
        &loads(&[0.7]),
        &cells(&caloric, &[30]),
        config(),
        &CancelGate::new_clock_free(),
    )
    .unwrap();
    assert!(retry.report().deposited_heat_j > 0.0);
}

#[test]
fn exact_source_law_inventory_clock_and_numerical_energy_allowances_are_required() {
    let (p, laws, caloric) = fixture(false, &[[6.0; 3]], 0);
    let viscosity = MembraneViscosity::new(&p, laws, (0.0, 1000.0), || false).unwrap();
    let original = motion(&viscosity, &[[0.2, -0.3]]);
    let saved = snapshot(&original);
    let heat = thermal(&caloric, &[(30, 3.0, 0.001)], vec![]);
    let same_numbers = caloric.clone();
    let gate = CancelGate::new_clock_free();
    assert!(
        advance_candidate(
            &original,
            &heat,
            &(),
            0.04,
            &loads(&[0.7]),
            &cells(&same_numbers, &[30]),
            config(),
            &gate
        )
        .is_err()
    );
    let (_, _, foreign) = fixture(false, &[[6.0; 3]], 1);
    let foreign_heat = thermal(&foreign, &[(30, 3.0, 0.001)], vec![]);
    assert!(
        advance_candidate(
            &original,
            &foreign_heat,
            &(),
            0.04,
            &loads(&[0.7]),
            &cells(&foreign, &[30]),
            config(),
            &gate
        )
        .is_err()
    );
    for (mass, volume) in [(1.0, 0.001), (3.0, 0.002)] {
        let wrong = thermal(&caloric, &[(30, mass, volume)], vec![]);
        assert!(
            advance_candidate(
                &original,
                &wrong,
                &(),
                0.04,
                &loads(&[0.7]),
                &cells(&caloric, &[30]),
                config(),
                &gate
            )
            .is_err()
        );
    }
    assert!(
        advance_candidate(
            &original,
            &heat,
            &(),
            0.04,
            &loads(&[0.7]),
            &[],
            config(),
            &gate
        )
        .is_err()
    );
    let strict = MembraneThermalConfig {
        energy_tolerance_j: 1e-8,
        ..config()
    };
    assert!(matches!(
        advance_candidate(
            &original,
            &heat,
            &(),
            0.04,
            &loads(&[0.7]),
            &cells(&caloric, &[30]),
            strict,
            &gate
        ),
        Err(MembraneThermalError::Energy { .. })
    ));
    let mut later = original.clone();
    later.step(0.005, &loads(&[0.7]), || false).unwrap();
    assert!(
        advance_candidate(
            &later,
            &heat,
            &(),
            0.04,
            &loads(&[0.7]),
            &cells(&caloric, &[30]),
            config(),
            &gate
        )
        .is_err()
    );
    assert_eq!(snapshot(&original), saved);
    assert_eq!(heat.time_s(), 0.0);
}
