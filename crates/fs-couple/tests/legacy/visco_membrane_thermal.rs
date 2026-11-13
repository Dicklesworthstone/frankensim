//! Synthetic finite-strain Maxwell heat checks; no measured material claims.
#[path = "visco_membrane_thermal/temperature.rs"]
mod temperature;

use fs_couple::membrane_thermal::{
    ElementThermalCell, MembraneThermalConfig, MembraneThermalError,
    visco::{ViscoMembraneThermalControls, advance_candidate},
};
use fs_evidence::ValidityDomain;
use fs_exec::CancelGate;
use fs_matdb::{
    ClaimSet, ConstitutiveModelCard, InitialStatePolicy, InterpolationPolicy, LawId, LawParameter,
    MaterialCard, MaterialStateId, PropertyClaim, PropertyKey, PropertyValue, Provenance,
    QueryPoint, UncertaintyModel,
};
use fs_material::hyper::plane_stress::PlaneStressSettings;
use fs_material::sheet::{
    AREAL_DENSITY, MembraneLaw, SheetCapabilities, THICKNESS, resolve_sheet, thermal::SheetCaloric,
};
use fs_material::state_point::{
    MaterialPropertySelection, SPECIFIC_HEAT_CAPACITY_DIMS, SPECIFIC_HEAT_CAPACITY_PROPERTY as CP,
    ScalarAdmissibility, ScalarPropertyRequirement, resolve_material_state_point,
};
use fs_material::visco_hyper::{
    BranchEnergy, FiniteStrainBranch, ViscoHyperelastic, membrane::ViscoMembraneState,
};
use fs_qty::Dims;
use fs_solid::membrane::{
    MembraneBudget, MembraneDomain, MembranePatch, MembraneSolveSettings,
    inertia::MembraneMassImage,
    visco::{
        FiniteViscoMembrane, FiniteViscoStepControls, ViscoelasticLawField,
        dynamics::{FiniteViscoMembraneMotion, FiniteViscoMotionControls},
    },
};
use fs_solid::mesh2::Mesh2;
use fs_solver::nonlinear::NewtonKrylovConfig;
use fs_thermal::enthalpy::{
    Boundary, CellGeometry, CellId, DepositedCell, DepositionBudget, DepositionError,
    DepositionObserver, DepositionSimulation, Event, MaterialId, StepConfig, ThermalLaw,
    TimedEvent,
};
use std::collections::BTreeMap;

const AREA_M2: f64 = 0.005;
const THICKNESS_M: f64 = 0.02;

fn local() -> PlaneStressSettings {
    PlaneStressSettings {
        relative_tolerance: 2e-13,
        ..Default::default()
    }
}

fn fixture(
    variable_cp: bool,
    density: &[[f64; 3]],
    relaxation_s: &[f64],
    revision: u32,
) -> (MembranePatch, Vec<ViscoHyperelastic>, SheetCaloric) {
    assert_eq!(density.len(), relaxation_s.len());
    let provenance = Provenance {
        source: "synthetic finite-strain membrane and calorics; no specimen authority".into(),
        license: "test-only".into(),
        artifact: None,
    };
    let pa = Dims([-1, 1, -2, 0, 0, 0]);
    let model = ConstitutiveModelCard {
        law: LawId("mooney-rivlin-decoupled".into()),
        law_version: 1,
        parameters: [
            ("c10", 10.0, pa),
            ("c01", 2.0, pa),
            ("kappa", 200.0, pa),
            ("stretch_limit", 4.0, Dims::NONE),
        ]
        .into_iter()
        .map(|(name, value, dims)| (name.into(), LawParameter { value, dims }))
        .collect(),
        state_schema_version: 1,
        initial_state: InitialStatePolicy::ZeroInternalState,
        validity: ValidityDomain::unconstrained().with("T", 299.0, 301.0),
        sources: vec![],
        provenance: provenance.clone(),
    };
    let mut claims = ClaimSet::new();
    for (name, value, dims) in [
        (AREAL_DENSITY, 0.8, Dims([-2, 1, 0, 0, 0, 0])),
        (THICKNESS, THICKNESS_M, Dims([1, 0, 0, 0, 0, 0])),
    ] {
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value: PropertyValue::Scalar { value, dims },
                interpolation: InterpolationPolicy::ConstantWithinValidity,
                validity: ValidityDomain::unconstrained().with("T", 299.0, 301.0),
                uncertainty: UncertaintyModel::Unstated,
                observations: vec![],
                provenance: provenance.clone(),
            })
            .unwrap();
    }
    claims
        .insert_claim(PropertyClaim {
            key: PropertyKey::new(CP, SPECIFIC_HEAT_CAPACITY_DIMS),
            value: if variable_cp {
                PropertyValue::Curve {
                    abscissa: "T".into(),
                    abscissa_dims: Dims([0, 0, 0, 1, 0, 0]),
                    knots: vec![(299.0, 0.015), (300.0, 0.020), (301.0, 0.025)],
                    dims: SPECIFIC_HEAT_CAPACITY_DIMS,
                }
            } else {
                PropertyValue::Scalar {
                    value: 1.0,
                    dims: SPECIFIC_HEAT_CAPACITY_DIMS,
                }
            },
            interpolation: if variable_cp {
                InterpolationPolicy::LinearInside
            } else {
                InterpolationPolicy::ConstantWithinValidity
            },
            validity: ValidityDomain::unconstrained().with("T", 299.0, 301.0),
            uncertainty: UncertaintyModel::Unstated,
            observations: vec![],
            provenance,
        })
        .unwrap();
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic membrane".into(),
            phase: "solid".into(),
            process: "finite-strain relaxation heat test".into(),
            revision,
        },
        claims,
        vec![model],
    )
    .unwrap();
    let query = QueryPoint::new().with("T", 300.0).unwrap();
    let sheet = resolve_sheet(
        &card,
        &query,
        SheetCapabilities {
            membrane: false,
            bending: false,
        },
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap()
    .with_incompressible_hyperelastic_membrane(&card, 0, local())
    .unwrap();
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
    let MembraneLaw::Hyperelastic(equilibrium) = sheet.membrane_law().unwrap() else {
        unreachable!()
    };
    let laws = relaxation_s
        .iter()
        .map(|&tau_s| {
            ViscoHyperelastic::new(
                equilibrium.model().law().clone(),
                vec![FiniteStrainBranch {
                    energy: BranchEnergy::NeoHookean { mu_pa: 12.0 },
                    tau_s,
                }],
            )
            .unwrap()
        })
        .collect();
    let count = density.len();
    let patch = MembranePatch::try_new(
        Mesh2 {
            nodes: (0..count)
                .flat_map(|e| {
                    let x = 0.2 * e as f64;
                    [[x, 0.0], [x + 0.1, 0.0], [x, 0.1]]
                })
                .collect(),
            elems: (0..count)
                .map(|e| vec![3 * e, 3 * e + 1, 3 * e + 2])
                .collect(),
            patches: vec![],
        },
        vec![sheet; count],
        vec![0.0; count],
        MembraneDomain {
            minimum_area_ratio: 0.01,
            maximum_abs_strain: 2.0,
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
    (patch, laws, caloric)
}

fn motion<'a>(
    patch: &'a MembranePatch,
    laws: &'a [ViscoHyperelastic],
    strained: bool,
) -> FiniteViscoMembraneMotion<'a> {
    let gradient = [[0.08, 0.03], [0.01, 0.04], [0.03, -0.02]];
    let u: Vec<_> = patch
        .mesh()
        .nodes
        .iter()
        .map(|x| {
            core::array::from_fn(|i| {
                if strained {
                    gradient[i][0] * x[0] + gradient[i][1] * x[1]
                } else {
                    0.0
                }
            })
        })
        .collect();
    let membrane = FiniteViscoMembrane::try_new_from(
        patch,
        ViscoelasticLawField::PerElement(laws),
        &u,
        || false,
    )
    .unwrap();
    FiniteViscoMembraneMotion::try_new(
        membrane,
        &vec![[0.0; 3]; u.len()],
        MembraneMassImage::Consistent,
        || false,
    )
    .unwrap()
}

fn held(motion: &FiniteViscoMembraneMotion<'_>) -> Vec<(usize, f64)> {
    motion
        .membrane()
        .displacement_m()
        .iter()
        .flatten()
        .copied()
        .enumerate()
        .collect()
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
    inventory: &[(u64, f64)],
    extra: Vec<TimedEvent>,
) -> DepositionSimulation<'a> {
    let initial_h = caloric
        .state_at_temperature(300.0)
        .unwrap()
        .specific_enthalpy_j_kg();
    let mut schedule: Vec<_> = inventory
        .iter()
        .enumerate()
        .map(|(e, &(id, mass_kg))| TimedEvent {
            time_s: 0.0,
            event: Event::Birth(DepositedCell {
                id: CellId(id),
                segment_id: id,
                body_id: 7,
                material_id: MaterialId(0),
                mass_kg,
                specific_h_j_kg: initial_h,
                geometry: CellGeometry {
                    start_m: [0.2 * e as f64, 0.0, 0.0],
                    end_m: [0.2 * e as f64 + 0.1, 0.0, 0.0],
                    volume_m3: AREA_M2 * THICKNESS_M,
                    surface_area_m2: 0.01,
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
    .advance_to(0.0, config().thermal, || Ok(()))
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

fn mechanics() -> ViscoMembraneThermalControls {
    ViscoMembraneThermalControls {
        solve: MembraneSolveSettings {
            newton: NewtonKrylovConfig {
                relative_tolerance: 1e-10,
                absolute_tolerance: 1e-11,
                linear_restart: 16,
                ..Default::default()
            },
            max_iterations: 60,
            max_evaluations: 4096,
            max_krylov_restart: 64,
        },
        motion: FiniteViscoMotionControls {
            material: FiniteViscoStepControls {
                equilibrium_work_tolerance_j: 1e-9,
                maximum_endpoint_work_defect_j: 10.0,
            },
            maximum_numerical_kinetic_loss_j: 1.0,
        },
        maximum_step_s: 0.04,
        pressure_pa: 0.0,
    }
}

fn config() -> MembraneThermalConfig {
    MembraneThermalConfig {
        thermal: StepConfig {
            max_dt_s: 0.04,
            energy_tolerance_j: 1e-10,
        },
        // Explicitly admit the reported BE work defect without heating from it.
        energy_tolerance_j: 0.1,
    }
}

fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(actual.is_finite() && expected.is_finite());
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:e} != {expected:e}, tolerance {tolerance:e}"
    );
}

#[test]
fn regional_finite_strain_relaxation_heats_original_mass_with_variable_cp() {
    let density = [[0.8; 3], [1.0, 1.6, 2.2]];
    let (patch, laws, caloric) = fixture(true, &density, &[0.015, 0.5], 0);
    let original = motion(&patch, &laws, true);
    let fixed = held(&original);
    // Triangle order deliberately opposes CellId order; the second P1 density
    // varies within its facet. Thermal receivers own each complete facet mass.
    let ids = [90, 10];
    let masses = density.map(|rho| AREA_M2 * rho.iter().sum::<f64>() / 3.0);
    let heat = thermal(
        &caloric,
        &[(ids[0], masses[0]), (ids[1], masses[1])],
        vec![],
    );
    let dt = 0.04;
    let expected: Vec<_> = original
        .membrane()
        .triangle_states()
        .iter()
        .zip(&laws)
        .map(|(state, law)| {
            let ViscoMembraneState::Incompressible(state) = state else {
                unreachable!()
            };
            let f = state.deformation_gradient;
            let response = law
                .step_incompressible_plane_stress_response(
                    state,
                    &[f[0], f[1], f[3], f[4]],
                    dt,
                    None,
                    local(),
                )
                .unwrap();
            AREA_M2 * THICKNESS_M * response.material.dissipated_pa
        })
        .collect();
    assert!(expected[0] > expected[1] && expected[1] > 0.0);
    let next = advance_candidate(
        &original,
        &heat,
        &(),
        dt,
        &vec![[0.0; 3]; patch.mesh().nodes.len()],
        &fixed,
        &cells(&caloric, &ids),
        mechanics(),
        config(),
        &CancelGate::new_clock_free(),
    )
    .unwrap();
    for e in 0..ids.len() {
        let before = &heat.cells()[&CellId(ids[e])];
        let after = &next.thermal().cells()[&CellId(ids[e])];
        let delta_h = expected[e] / masses[e];
        // Independently integrate Cp(T)=.02+.005*(T-300) and invert the
        // resulting .02*dT+.0025*dT². No constant-Cp shortcut is valid here.
        let delta_t = 2.0 * delta_h / (0.02 + (0.02_f64.powi(2) + 0.01 * delta_h).sqrt());
        near(next.report().element_heat_j[e], expected[e], 2e-12);
        near(
            after.specific_h_j_kg - before.specific_h_j_kg,
            delta_h,
            2e-10,
        );
        let temperature = next.thermal().temperature_k(CellId(ids[e])).unwrap();
        // The source chart's 1e-9 J/kg error bounds T error by 1e-9/.015 K.
        near(temperature, 300.0 + delta_t, 8e-8);
        assert!((temperature - (300.0 + delta_h / 0.02)).abs() > 1e-5);
        assert_eq!(after.mass_kg, before.mass_kg);
        assert_eq!(after.geometry, before.geometry);
    }
    let total: f64 = expected.iter().sum();
    let report = next.report();
    near(report.material_heat_j, total, 2e-12);
    near(report.deposited_heat_j, total, 2e-12);
    near(report.thermal.stored_energy_change_j, total, 2e-12);
    near(report.mechanical_energy_change_j, -total, 2e-12);
    near(report.mechanical_work_j, 0.0, 2e-12);
    near(report.numerical_kinetic_loss_j, 0.0, 2e-12);
    near(report.mechanical_numerical_defect_j, 0.0, 2e-12);
    near(report.heat_handoff_residual_j, 0.0, 2e-12);
    near(report.external_thermal_heat_j, 0.0, 2e-12);
    near(report.combined_residual_j, 0.0, 2e-12);
    assert_eq!(
        next.motion().membrane().displacement_m(),
        original.membrane().displacement_m()
    );
    assert_eq!(next.motion().membrane().time_s(), dt);
    assert_eq!(next.thermal().time_s(), dt);
    assert_eq!(original.membrane().time_s(), 0.0);
    assert_eq!(heat.time_s(), 0.0);
}

#[test]
fn uniform_acceleration_reports_be_kinetic_loss_without_turning_it_into_heat() {
    let (patch, laws, caloric) = fixture(false, &[[0.8; 3]], &[0.25], 0);
    let original = motion(&patch, &laws, false);
    let heat = thermal(&caloric, &[(30, AREA_M2 * 0.8)], vec![]);
    let acceleration = [1.2, -0.6, 0.2];
    let load = patch
        .reference_mass_apply(
            MembraneMassImage::Consistent,
            &vec![acceleration; patch.mesh().nodes.len()],
            || false,
        )
        .unwrap();
    let dt = 0.04;
    let next = advance_candidate(
        &original,
        &heat,
        &(),
        dt,
        &load,
        &[],
        &cells(&caloric, &[30]),
        mechanics(),
        config(),
        &CancelGate::new_clock_free(),
    )
    .unwrap();
    for (u, v) in next
        .motion()
        .membrane()
        .displacement_m()
        .iter()
        .zip(next.motion().velocity_mps())
    {
        for d in 0..3 {
            near(u[d], acceleration[d] * dt * dt, 2e-10);
            near(v[d], acceleration[d] * dt, 2e-9);
        }
    }
    let loss = 0.5 * (AREA_M2 * 0.8) * acceleration.iter().map(|a| a * a).sum::<f64>() * dt * dt;
    let report = next.report();
    assert!(report.numerical_kinetic_loss_j > 0.0);
    near(report.numerical_kinetic_loss_j, loss, 2e-12);
    near(report.mechanical_energy_change_j, loss, 2e-12);
    near(report.mechanical_work_j, 2.0 * loss, 2e-12);
    near(report.mechanical_numerical_defect_j, -loss, 2e-12);
    near(report.mechanical_residual_j, -loss, 2e-12);
    near(report.absolute_mechanical_defect_j, loss, 2e-12);
    near(report.combined_residual_j, -loss, 2e-12);
    near(report.element_heat_j[0], 0.0, 2e-13);
    near(report.material_heat_j, 0.0, 2e-13);
    near(report.deposited_heat_j, 0.0, 2e-13);
    near(report.thermal.stored_energy_change_j, 0.0, 2e-13);
    near(
        next.thermal().temperature_k(CellId(30)).unwrap(),
        300.0,
        2e-10,
    );
}

#[derive(Clone, Debug, Default, PartialEq)]
struct Observer {
    intervals: Vec<(f64, f64, f64)>,
    events: usize,
    steps: usize,
    fail_after: Option<f64>,
}

impl DepositionObserver for Observer {
    type Error = DepositionError;

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
        Ok(BTreeMap::new())
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
        _: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), Self::Error> {
        self.steps += 1;
        if self.fail_after.is_some_and(|time| thermal.time_s() > time) {
            return Err(DepositionError::InvalidInput(
                "synthetic late observer refusal",
            ));
        }
        Ok(())
    }
}

#[test]
fn thermal_events_and_observer_intervals_advance_the_same_retained_maxwell_history() {
    let (patch, laws, caloric) = fixture(false, &[[0.8; 3]], &[0.05], 0);
    let mut original = motion(&patch, &laws, true);
    let fixed = held(&original);
    let load = vec![[0.0; 3]; patch.mesh().nodes.len()];
    let switch = 0.0234375;
    let source_w = 0.0002;
    let mass = AREA_M2 * 0.8;
    let mut heat = thermal(
        &caloric,
        &[(30, mass)],
        vec![TimedEvent {
            time_s: switch,
            event: Event::Boundary(CellId(30), boundary(source_w)),
        }],
    );
    let mut observer = Observer::default();
    let controls = ViscoMembraneThermalControls {
        maximum_step_s: 0.03125,
        ..mechanics()
    };
    let mut expected_h = heat.cells()[&CellId(30)].specific_h_j_kg;
    for target in [0.0625, 0.09375] {
        let before = observer.intervals.len();
        let start = heat.time_s();
        let next = advance_candidate(
            &original,
            &heat,
            &observer,
            target,
            &load,
            &fixed,
            &cells(&caloric, &[30]),
            controls,
            config(),
            &CancelGate::new_clock_free(),
        )
        .unwrap();
        let intervals = &next.observer().intervals[before..];
        assert_eq!(intervals.len(), if before == 0 { 5 } else { 2 });
        let mut replay = None;
        let mut interval_heat = 0.0;
        for &(time, dt, end) in intervals {
            assert_eq!(time + dt, end);
            let accepted = replay.as_ref().unwrap_or(&original);
            assert_eq!(accepted.membrane().time_s(), time);
            // Replay the actual native thermal intervals, never Newton trials
            // or one coarse call-wide material interval.
            let step = accepted
                .advance(dt, &load, &fixed, controls.solve, controls.motion, || false)
                .unwrap();
            interval_heat += step.material_energy.dissipated_j;
            replay = Some(step.motion);
        }
        let replay = replay.unwrap();
        assert_eq!(
            next.motion().membrane().triangle_states(),
            replay.membrane().triangle_states()
        );
        near(
            next.motion().membrane().stored_energy_j(),
            replay.membrane().stored_energy_j(),
            2e-12,
        );
        let external = source_w * (target - start.max(switch));
        expected_h += (interval_heat + external) / mass;
        near(next.report().element_heat_j[0], interval_heat, 2e-12);
        near(next.report().material_heat_j, interval_heat, 2e-12);
        near(next.report().deposited_heat_j, interval_heat, 2e-12);
        near(next.report().external_thermal_heat_j, external, 2e-12);
        near(
            next.thermal().cells()[&CellId(30)].specific_h_j_kg,
            expected_h,
            2e-10,
        );
        assert_eq!(next.motion().membrane().time_s(), target);
        assert_eq!(next.thermal().time_s(), target);
        assert_eq!(next.observer().steps, next.observer().intervals.len());
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
    storage: f64,
    states: Vec<ViscoMembraneState>,
}

fn snapshot(motion: &FiniteViscoMembraneMotion<'_>) -> MechanicalSnapshot {
    MechanicalSnapshot {
        displacement: motion.membrane().displacement_m().to_vec(),
        velocity: motion.velocity_mps().to_vec(),
        time: motion.membrane().time_s(),
        storage: motion.membrane().stored_energy_j(),
        states: motion.membrane().triangle_states().to_vec(),
    }
}

#[test]
fn source_inventory_and_late_observer_refusals_preserve_all_accepted_owners() {
    let (patch, laws, caloric) = fixture(false, &[[0.8; 3]], &[0.05], 0);
    let original = motion(&patch, &laws, true);
    let fixed = held(&original);
    let load = vec![[0.0; 3]; patch.mesh().nodes.len()];
    let heat = thermal(&caloric, &[(30, AREA_M2 * 0.8)], vec![]);
    let saved_motion = snapshot(&original);
    let saved_cells = heat.cells().clone();
    let observer = Observer {
        fail_after: Some(0.02),
        ..Observer::default()
    };
    let saved_observer = observer.clone();
    let gate = CancelGate::new_clock_free();
    let (_, _, foreign) = fixture(false, &[[0.8; 3]], &[0.05], 1);
    let foreign_heat = thermal(&foreign, &[(30, AREA_M2 * 0.8)], vec![]);
    let saved_foreign_cells = foreign_heat.cells().clone();
    assert!(matches!(
        advance_candidate(
            &original,
            &foreign_heat,
            &observer,
            0.04,
            &load,
            &fixed,
            &cells(&foreign, &[30]),
            mechanics(),
            config(),
            &gate,
        ),
        Err(MembraneThermalError::Invalid(_))
    ));
    assert_eq!(foreign_heat.cells(), &saved_foreign_cells);
    assert_eq!(foreign_heat.time_s(), 0.0);
    let wrong_mass = thermal(&caloric, &[(30, 0.5 * AREA_M2 * 0.8)], vec![]);
    let saved_wrong_cells = wrong_mass.cells().clone();
    assert!(matches!(
        advance_candidate(
            &original,
            &wrong_mass,
            &observer,
            0.04,
            &load,
            &fixed,
            &cells(&caloric, &[30]),
            mechanics(),
            config(),
            &gate,
        ),
        Err(MembraneThermalError::Invalid(_))
    ));
    assert_eq!(wrong_mass.cells(), &saved_wrong_cells);
    assert_eq!(wrong_mass.time_s(), 0.0);
    assert!(matches!(
        advance_candidate(
            &original,
            &heat,
            &observer,
            0.04,
            &load,
            &fixed,
            &cells(&caloric, &[30]),
            mechanics(),
            config(),
            &gate,
        ),
        Err(MembraneThermalError::Observer(
            DepositionError::InvalidInput("synthetic late observer refusal")
        ))
    ));
    assert_eq!(snapshot(&original), saved_motion);
    assert_eq!(heat.cells(), &saved_cells);
    assert_eq!(heat.time_s(), 0.0);
    assert_eq!(observer, saved_observer);
    let retry = advance_candidate(
        &original,
        &heat,
        &(),
        0.04,
        &load,
        &fixed,
        &cells(&caloric, &[30]),
        mechanics(),
        config(),
        &gate,
    )
    .unwrap();
    assert!(retry.report().deposited_heat_j > 0.0);
}
