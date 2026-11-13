//! Complete material heat on the same native rod potentials and motion history.
use super::*;
use fs_couple::circuit_enthalpy::{CircuitEnthalpyError, CircuitThermalCoefficients,
    contacts::{ContactThermalCoefficients, material::ElectricalContactNetwork,
        heat::{HeatAssembly, HeatNetworkError, HeatNetworkSample, MaterialHeatNetwork, PreparedHeatNetwork}}};
use fs_couple::rod_motion::enthalpy::{
    RodAngularViscosity, RodBendingViscosity, RodHeatError, RodTorsionalViscosity,
    RodTransverseViscosity,
};
use fs_conduction::transient::deposition::{
    contact::{ContactCell, MaterialContactLink},
    radiation::{MaterialRadiation, RadiatingCell},
};
use fs_conduction::{GrayDiffuseEnclosure, RadiationSurface, SurfaceEmissivity, ViewFactorMatrix,
    AREA_SPECIFIC_THERMAL_RESISTANCE_DIMS, AREA_SPECIFIC_THERMAL_RESISTANCE_PROPERTY,
    STEFAN_BOLTZMANN_W_M2_K4};
use fs_exec::Cx;
use fs_matdb::{ClaimSelection, InterfaceSystemCard, SurfaceSpec, SystemContext};
use fs_material::thermal::interface::ThermalInterfacePoint;
use fs_solid::rod::material_field::motion::thermal::RodSegmentStrainState;
use std::cell::Cell;
use std::collections::BTreeMap;

fn with_cx(gate: &CancelGate, run: impl FnOnce(&Cx<'_>)) {
    fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default()).scope(|arena| {
        let cx = Cx::new(gate, arena,
            fs_exec::StreamKey { seed: 917, kernel_id: 63, tile: 0, iteration: 0 },
            fs_exec::Budget::INFINITE, fs_exec::ExecMode::Deterministic);
        run(&cx);
    });
}

fn interface(cards: [&MaterialCard; 2], resistance: f64) -> InterfaceSystemCard {
    let mut claims = ClaimSet::new();
    claims.insert_claim(PropertyClaim {
        key: PropertyKey::new(AREA_SPECIFIC_THERMAL_RESISTANCE_PROPERTY,
            AREA_SPECIFIC_THERMAL_RESISTANCE_DIMS),
        value: PropertyValue::Scalar { value: resistance, dims: AREA_SPECIFIC_THERMAL_RESISTANCE_DIMS },
        validity: ValidityDomain::unconstrained().with("Ta", 200.0, 500.0).with("Tb", 200.0, 500.0),
        interpolation: InterpolationPolicy::ConstantWithinValidity,
        uncertainty: UncertaintyModel::Unstated, observations: vec![],
        provenance: Provenance { source: "manufactured native rod contact".into(),
            license: "test-only".into(), artifact: None },
    }).unwrap();
    InterfaceSystemCard::assemble(
        SurfaceSpec { material: cards[0].id().clone(), texture_frame: "a".into() },
        SurfaceSpec { material: cards[1].id().clone(), texture_frame: "b".into() },
        SystemContext { medium: "manufactured".into(), third_body: None,
            environment: "fixed reference geometry".into(), history: "original".into() },
        claims, vec![]).unwrap()
}

#[test]
fn complete_contact_radiation_heat_refreshes_native_rod_stiffness_and_thermal_work() {
    let specs = [Spec::default(), Spec { shear_scale: 1.4, cp: 0.15, ..Spec::default() }];
    let initial = [310.0, 350.0];
    let strain = 0.012;
    let dt = 0.001;
    let mut final_forces = Vec::new();
    for resistance in [0.2, 2.0] {
        fixture(&specs, 1, &initial, &[0.0; 2], true, |source, mut thermal, network| {
            let gate = CancelGate::new_clock_free();
            let motion = axial(source, &initial, strain, true);
            let pose = motion.pose().clone();
            let inertia = motion.model().inertia().to_vec();
            let mut owner = RodThermoelasticSystem::try_new(motion, source, &thermal,
                network, &destinations(2, 1), 1e-11, &gate).unwrap();
            let electrical = ElectricalContactNetwork { bulk: network, contacts: &[], electrical_bases: &[] };
            let card = interface([network.cells[0].point.material, network.cells[1].point.material], resistance);
            let context = QueryPoint::new().with("Ta", initial[0]).unwrap().with("Tb", initial[1]).unwrap();
            let contacts = [MaterialContactLink {
                link: 17,
                sides: [0, 1].map(|i| ContactCell { cell: CellId(i as u64),
                    material_id: network.cells[i].material_id, curve: network.cells[i].curve }),
                point: ThermalInterfacePoint { interface: &card,
                    materials: [network.cells[0].point.material, network.cells[1].point.material],
                    context: &context, temperature_axes: ["Ta", "Tb"],
                    liquid_fraction_axes: [None, None], mean_temperature_axis: None },
                selection: ClaimSelection::Policy(SelectionPolicy::SingleClaimOnly),
            }];
            let radiating = [0, 1].map(|i| RadiatingCell { surface: i, cell: CellId(i as u64),
                material_id: network.cells[i].material_id, curve: network.cells[i].curve,
                material: network.cells[i].point, selection: MaterialPropertySelection::SingleClaimOnly });
            let (complex, positions) = fs_conduction::fixtures::box_grid([1, 1, 1], [0.1; 3]);
            let mesh = fs_conduction::mesh::ConductionMesh::new(complex, positions).unwrap();
            let surfaces = (0..2).map(|i| {
                let emissivity = SurfaceEmissivity::from_thermal_state("rod test", network.cells[i].point,
                    network.cells[i].curve.state_at_specific_enthalpy(thermal.cells()[&CellId(i as u64)].specific_h_j_kg).unwrap(),
                    MaterialPropertySelection::SingleClaimOnly).unwrap().0;
                RadiationSurface::new(&mesh, format!("face {i}"),
                    |face| fs_conduction::fixtures::on_box_face(face.centroid[0], 0.1 * i as f64), emissivity).unwrap()
            }).collect();
            let enclosure = GrayDiffuseEnclosure::new(surfaces,
                ViewFactorMatrix::infinite_parallel_plates(0.01).unwrap()).unwrap();
            let radiation = MaterialRadiation { enclosure: &enclosure, cells: &radiating,
                power_tolerance_w: 1e-10, radiosity_tolerance_w_m2: 1e-8 };
            let heat = MaterialHeatNetwork { electrical: &electrical, thermal_contacts: &contacts,
                radiation: Some(&radiation) };
            with_cx(&gate, |cx| {
                let prepared = heat.prepare(cx, &thermal).unwrap();
                let g = 0.01 / (0.02 / specs[0].conductivity + 0.03 / specs[1].conductivity + resistance);
                let radiation_factor = STEFAN_BOLTZMANN_W_M2_K4 * 0.01 / (2.0 / 0.6 - 1.0);
                let potential = |i: usize, t| {
                    let (a, second) = geometry(i, 2);
                    independent_potential(specs[i], a, second, t, [strain, 0.0, 0.0], [0.0; 3])
                };
                let capacity: [f64; 2] = core::array::from_fn(|i| RHO * geometry(i, 2).0 * LENGTH * specs[i].cp);
                let initial_total = owner.mechanical_internal_energy_j().unwrap() + total_heat(&thermal);
                let loads = [RodNodalWrench::default(); 2];
                let spatial = PreparedRodLoads::new(&loads, &[[0.0; 3]; 2], &[], &gate).unwrap();
                let mut old = initial;
                for step in 1..=2 {
                    let mut expected = old;
                    // Independent two-cell Newton equations for caloric U plus
                    // the original closed-form mechanical Psi+TS potential.
                    for _ in 0..20 {
                        let p = [potential(0, expected[0]), potential(1, expected[1])];
                        let q = dt * (g * (expected[1] - expected[0])
                            + radiation_factor * (expected[1].powi(4) - expected[0].powi(4)));
                        let r = [capacity[0] * (expected[0] - old[0]) + LENGTH * (p[0].2 - potential(0, old[0]).2) - q,
                            capacity[1] * (expected[1] - old[1]) + LENGTH * (p[1].2 - potential(1, old[1]).2) + q];
                        let h = expected.map(|t| dt * (g + 4.0 * radiation_factor * t.powi(3)));
                        let d = [capacity[0] + LENGTH * p[0].3 + h[0], capacity[1] + LENGTH * p[1].3 + h[1]];
                        let determinant = d[0] * d[1] - h[0] * h[1];
                        expected[0] -= (d[1] * r[0] + h[1] * r[1]) / determinant;
                        expected[1] -= (d[0] * r[1] + h[0] * r[0]) / determinant;
                    }
                    let report = if step == 1 {
                        owner.step_heat_assembly(&mut thermal, dt, &loads, &prepared, cx, None,
                            controls(dt), solve(), &[0.0; 2], &gate)
                    } else {
                        owner.step_spatial_heat_assembly(&mut thermal, 2.0 * dt, &spatial, &prepared, cx, None,
                            controls(dt), solve(), &[0.0; 2], &gate)
                    }.unwrap();
                    for i in 0..2 {
                        near(thermal.temperature_k(CellId(i as u64)).unwrap(), expected[i], 2e-7);
                        near(report.layers[i].states[1].internal_energy_j, LENGTH * potential(i, expected[i]).2, 2e-10);
                    }
                    near(report.source.thermal_contacts[0].interface.value_m2_k_per_w(), resistance, 1e-14);
                    near(report.source.thermal_contacts[0].conductance_w_k, g, 1e-12);
                    assert!(report.source.radiation.is_some());
                    assert!(report.entropy.radiation_exchange_entropy_j_k > 0.0);
                    assert!(!report.entropy.has_unclassified_external_sources);
                    near(report.closed_radiation_net_j, 0.0, 1e-12);
                    near(report.external_thermal_heat_j, 0.0, 3e-11);
                    near(report.total_energy_residual_j, 0.0, 5e-10);
                    near(report.entropy.accounting_residual_j_k, 0.0, 2e-12);
                    assert_eq!(owner.mechanical().pose(), &pose);
                    assert_eq!(owner.mechanical().model().inertia(), inertia);
                    assert_eq!(owner.mechanical().accepted_steps(), step);
                    old = expected;
                }
                near(owner.mechanical_internal_energy_j().unwrap() + total_heat(&thermal), initial_total, 5e-10);
                let expected_force: f64 = (0..2).map(|i| {
                    let young = YOUNG * specs[i].young_scale + specs[i].young_slope * (old[i] - 300.0);
                    -young * geometry(i, 2).0 * (strain - specs[i].alpha * (old[i] - 300.0))
                }).sum();
                let force = owner.mechanical().model().evaluate(owner.mechanical().pose(), &gate).unwrap().force_n[1][0];
                near(force, expected_force, 2e-10);
                final_forces.push(force);
            });
        });
    }
    assert!((final_forces[0] - final_forces[1]).abs() > 1e-6);
}

#[derive(Debug)]
enum SourceError { Heat(HeatNetworkError), Refused }
impl From<HeatNetworkError> for SourceError {
    fn from(error: HeatNetworkError) -> Self { Self::Heat(error) }
}
struct Receipt { heat: HeatNetworkSample, conductance_w_k: f64 }
impl AsRef<HeatNetworkSample> for Receipt {
    fn as_ref(&self) -> &HeatNetworkSample { &self.heat }
}
struct FiniteHeat<'p, 'n, 'a> {
    prepared: &'p PreparedHeatNetwork<'n, 'a>,
    conductance_w_k: f64,
    calls: Cell<usize>,
    refuse: Cell<bool>,
    cancel: Cell<bool>,
    cancel_gate: &'p CancelGate,
}
impl HeatAssembly for FiniteHeat<'_, '_, '_> {
    type Sample = Receipt;
    type Error = SourceError;
    fn prepared_heat_network(&self) -> &PreparedHeatNetwork<'_, '_> { self.prepared }
    fn resolve_heat(&self, cx: &Cx<'_>, thermal: &DepositionSimulation<'_>)
        -> Result<ContactThermalCoefficients<Self::Sample>, CircuitEnthalpyError<Self::Error>> {
        self.calls.set(self.calls.get() + 1);
        if self.refuse.get() { return Err(CircuitEnthalpyError::Source(SourceError::Refused)); }
        let resolved = self.prepared.resolve_with_link_conductances(cx, thermal,
            &BTreeMap::from([(17, self.conductance_w_k)]))
            .map_err(|error| CircuitEnthalpyError::Source(SourceError::Heat(error)))?;
        if self.cancel.get() { self.cancel_gate.request(); }
        let bulk = resolved.bulk;
        Ok(ContactThermalCoefficients { contacts: resolved.contacts, bulk: CircuitThermalCoefficients {
            resistors: bulk.resistors, sources_w: bulk.sources_w, conductances_w_k: bulk.conductances_w_k,
            source: Receipt { heat: bulk.source, conductance_w_k: self.conductance_w_k },
        } })
    }
}

#[test]
fn finite_heat_assembly_changes_continuing_rod_motion_and_keeps_both_angular_heats_once() {
    finite_heat_assembly_case(false);
}

#[test]
fn finite_heat_assembly_cools_native_rod_while_depositing_transverse_heat_once() {
    finite_heat_assembly_case(true);
}

fn finite_heat_assembly_case(include_transverse: bool) {
    let spec = Spec { viscosity: Some(ViscositySpec { value_pa_s: 3.0, slope_pa_s_k: 0.01,
        lower_k: 200.0, upper_k: 500.0, upper_rad_s: f64::MAX }),
        bending_viscosity: Some(ViscositySpec { value_pa_s: 8.0, slope_pa_s_k: 0.02,
            lower_k: 200.0, upper_k: 500.0, upper_rad_s: f64::MAX }), ..Spec::default() };
    let old = [310.0, 350.0];
    let dt = 0.001;
    let mut final_motion = Vec::new();
    for conductance in [0.0, 0.02] {
        fixture(&[spec; 2], 1, &old, &[0.0; 2], true, |source, mut thermal, original| {
            let gate = CancelGate::new_clock_free();
            let cx_gate = CancelGate::new_clock_free();
            let network = MaterialNetwork { cells: original.cells, resistors: &[], links: &[] };
            // Admission certifies the extra link; its cold coefficient must not
            // replace the conductance supplied by the live assembly below.
            let extra = BTreeMap::from([(17, 0.0)]);
            let rod = source.resolve_at(&old, &gate).unwrap().with_supports(&[RodSupport {
                node: 1, position_m: [None, Some(0.0), Some(0.0)], frame: None,
            }], &gate).unwrap();
            let model = Arc::new(rod.into_motion_model(RodRootBoundary::Clamped, &gate).unwrap());
            let mut pose = model.initial_pose();
            pose.positions_m[1][0] = LENGTH * 1.012;
            let motion = RodMotion::new(model, pose, &[[0.0; 3]; 2], &[[0.0; 3], [0.2, 0.3, -0.25]], 0.0, &gate).unwrap();
            let mut owner = RodThermoelasticSystem::try_new_with_link_conductances(motion, source,
                &thermal, &network, &extra, &destinations(2, 1), 1e-11, &gate).unwrap();
            let electrical = ElectricalContactNetwork { bulk: &network, contacts: &[], electrical_bases: &[] };
            let heat = MaterialHeatNetwork { electrical: &electrical, thermal_contacts: &[], radiation: None };
            with_cx(&cx_gate, |cx| {
                let prepared = heat.prepare_with_link_conductances(cx, &thermal, &extra).unwrap();
                let mut assembly = FiniteHeat { prepared: &prepared, conductance_w_k: conductance,
                    calls: Cell::new(0), refuse: Cell::new(false), cancel: Cell::new(false),
                    cancel_gate: &cx_gate };
                let selections = vec![MaterialPropertySelection::SingleClaimOnly; 2];
                let viscosity = RodAngularViscosity {
                    torsional: Some(RodTorsionalViscosity { property: "shear_viscosity", selections: &selections }),
                    bending: Some(RodBendingViscosity { property: "bending_viscosity", selections: &selections }),
                };
                let transverse_viscosity = RodTransverseViscosity {
                    property: "shear_viscosity", selections: &selections,
                };
                let initial_total = owner.mechanical_internal_energy_j().unwrap() + total_heat(&thermal);
                let loads = [RodNodalWrench::default(); 2];
                let spatial = PreparedRodLoads::new(&loads, &[[0.0; 3]; 2], &[], &gate).unwrap();
                // A cold additional link is never silently frozen by an ordinary step.
                assert!(owner.step(&mut thermal, dt, &loads, controls(dt), solve(), &[0.0; 2], &gate).is_err());
                assert_eq!(owner.mechanical().accepted_steps(), 0);
                let mut trials = 0;
                let mut temperatures = old;
                let mut defects = 0.0;
                for step in 1..=2 {
                    let active_conductance = step as f64 * conductance;
                    assembly.conductance_w_k = active_conductance;
                    let before = owner.mechanical().model()
                        .segment_strain_states(owner.mechanical().pose(), &gate).unwrap()[0];
                    let report = match (include_transverse, step == 1) {
                        (true, true) => owner.step_heat_assembly_with_transverse(
                            &mut thermal, dt, &loads, &assembly, cx, viscosity,
                            transverse_viscosity, controls(dt), solve(), &[0.0; 2], &gate),
                        (true, false) => owner.step_spatial_heat_assembly_with_transverse(
                            &mut thermal, 2.0 * dt, &spatial, &assembly, cx, viscosity,
                            transverse_viscosity, controls(dt), solve(), &[0.0; 2], &gate),
                        (false, true) => owner.step_heat_assembly(&mut thermal, dt, &loads,
                            &assembly, cx, viscosity, controls(dt), solve(), &[0.0; 2], &gate),
                        (false, false) => owner.step_spatial_heat_assembly(&mut thermal, 2.0 * dt, &spatial,
                            &assembly, cx, viscosity, controls(dt), solve(), &[0.0; 2], &gate),
                    }.unwrap();
                    trials += report.material_trials;
                    assert_eq!(assembly.calls.get(), trials);
                    assert_eq!(report.source.conductance_w_k, active_conductance);
                    let torsional = report.torsional.as_ref().unwrap();
                    let bending = report.bending.as_ref().unwrap();
                    assert!(torsional.mechanical_heat_j > 0.0 && bending.mechanical_heat_j > 0.0);
                    near(torsional.deposited_heat_j, torsional.mechanical_heat_j, 1e-17);
                    near(bending.deposited_heat_j, bending.mechanical_heat_j, 1e-17);
                    assert_eq!(report.transverse.is_some(), include_transverse);
                    let transverse_total = if let Some(transverse) = &report.transverse {
                        assert!(transverse.mechanical_heat_j > 0.0);
                        near(transverse.deposited_heat_j, transverse.mechanical_heat_j, 1e-17);
                        transverse.mechanical_heat_j
                    } else { 0.0 };
                    let next: [f64; 2] = core::array::from_fn(|i| thermal.temperature_k(CellId(i as u64)).unwrap());
                    let after = owner.mechanical().model()
                        .segment_strain_states(owner.mechanical().pose(), &gate).unwrap()[0];
                    let eta_s = temperatures.map(|t| 3.0 + 0.01 * (t - 300.0));
                    let eta_b = temperatures.map(|t| 8.0 + 0.02 * (t - 300.0));
                    let torsional_weight: f64 = (0..2).map(|i| eta_s[i] * 2.0 * geometry(i, 2).1).sum();
                    let bending_weight: f64 = (0..2).map(|i| eta_b[i] * geometry(i, 2).1).sum();
                    let transverse_weight: f64 = (0..2).map(|i| 0.8 * eta_s[i] * geometry(i, 2).0).sum();
                    let mut caloric_entropy = 0.0;
                    let mut mechanical_entropy = 0.0;
                    let mut torsional_entropy = 0.0;
                    let mut bending_entropy = 0.0;
                    let mut transverse_entropy = 0.0;
                    for i in 0..2 {
                        let (area, second) = geometry(i, 2);
                        let capacity = RHO * LENGTH * area * spec.cp;
                        let potential = |t, strain: RodSegmentStrainState| {
                            independent_potential(spec, area, second, t, strain.gamma, strain.curvature)
                        };
                        let old_p = potential(temperatures[i], before);
                        let frozen_p = potential(temperatures[i], after);
                        let new_p = potential(next[i], after);
                        let reversible = LENGTH * (frozen_p.0 - new_p.0
                            - next[i] * new_p.1 + temperatures[i] * old_p.1);
                        let shear_heat = torsional.mechanical_heat_j * eta_s[i] * 2.0 * second / torsional_weight;
                        let bend_heat = bending.mechanical_heat_j * eta_b[i] * second / bending_weight;
                        let transverse_heat = transverse_total * (0.8 * eta_s[i] * area) / transverse_weight;
                        near(torsional.layers[i].viscosity_pa_s, eta_s[i], 1e-12);
                        near(bending.layers[i].viscosity_pa_s, eta_b[i], 1e-12);
                        near(torsional.layers[i].heat_j, shear_heat, 1e-17);
                        near(bending.layers[i].heat_j, bend_heat, 1e-17);
                        if let Some(transverse) = &report.transverse {
                            near(transverse.layers[i].viscosity_pa_s, eta_s[i], 1e-12);
                            near(transverse.layers[i].shear_correction, 0.8, 0.0);
                            near(transverse.layers[i].heat_j, transverse_heat, 1e-17);
                        }
                        near(report.layers[i].reversible_heat_j, reversible, 1e-12);
                        // Independent finite-cell balance: transport and reversible
                        // work remain distinct from each actual intrinsic loss.
                        let transport = dt * active_conductance * (next[1 - i] - next[i]);
                        near(capacity * (next[i] - temperatures[i]),
                            transport + reversible + shear_heat + bend_heat + transverse_heat, 3e-11);
                        caloric_entropy += capacity * (next[i] / temperatures[i]).ln();
                        mechanical_entropy += LENGTH * (new_p.1 - old_p.1);
                        torsional_entropy += shear_heat / next[i];
                        bending_entropy += bend_heat / next[i];
                        transverse_entropy += transverse_heat / next[i];
                    }
                    near(report.entropy.caloric_storage_change_j_k, caloric_entropy, 2e-12);
                    near(report.entropy.mechanical_storage_change_j_k, mechanical_entropy, 2e-12);
                    near(report.entropy.torsional_generation_j_k, torsional_entropy, 2e-15);
                    near(report.entropy.bending_generation_j_k, bending_entropy, 2e-15);
                    near(report.entropy.transverse_generation_j_k, transverse_entropy, 2e-15);
                    near(report.entropy.conduction_generation_j_k,
                        dt * active_conductance * (next[0] - next[1]).powi(2) / (next[0] * next[1]), 2e-12);
                    assert!(!report.entropy.has_unclassified_external_sources);
                    near(report.external_thermal_heat_j, 0.0, 5e-11);
                    near(report.coupling_energy_residual_j, 0.0, 5e-10);
                    near(report.entropy.accounting_residual_j_k, 0.0, 2e-12);
                    defects += report.mechanical.integration_defect_j;
                    temperatures = next;
                }
                assert!((owner.mechanical().pose().positions_m[1][0] - LENGTH * 1.012).abs() > 1e-7);
                // Native integration error stays separate from physical heat.
                near(owner.mechanical_internal_energy_j().unwrap() + total_heat(&thermal) - initial_total,
                    defects, 1e-9);
                let saved = (owner.mechanical().pose().clone(), owner.mechanical().linear_momenta_kg_m_s().to_vec(),
                    owner.mechanical().body_spins_kg_m2_s().to_vec(), owner.layers().to_vec(), thermal.cells().clone());
                assembly.refuse.set(true);
                let refused = if include_transverse {
                    owner.step_heat_assembly_with_transverse(&mut thermal, 3.0 * dt, &loads, &assembly, cx,
                        viscosity, transverse_viscosity, controls(dt), solve(), &[0.0; 2], &gate)
                } else {
                    owner.step_heat_assembly(&mut thermal, 3.0 * dt, &loads, &assembly, cx,
                        viscosity, controls(dt), solve(), &[0.0; 2], &gate)
                };
                assert!(matches!(refused, Err(RodHeatError::Heat(SourceError::Refused))));
                assembly.refuse.set(false);
                assembly.cancel.set(true);
                let cancelled = if include_transverse {
                    owner.step_spatial_heat_assembly_with_transverse(&mut thermal, 3.0 * dt, &spatial, &assembly, cx,
                        viscosity, transverse_viscosity, controls(dt), solve(), &[0.0; 2], &gate)
                } else {
                    owner.step_spatial_heat_assembly(&mut thermal, 3.0 * dt, &spatial, &assembly, cx,
                        viscosity, controls(dt), solve(), &[0.0; 2], &gate)
                };
                assert!(matches!(cancelled, Err(RodHeatError::Rod(CircuitEnthalpyError::Cancelled))));
                assert_eq!(owner.mechanical().pose(), &saved.0);
                assert_eq!(owner.mechanical().linear_momenta_kg_m_s(), saved.1);
                assert_eq!(owner.mechanical().body_spins_kg_m2_s(), saved.2);
                assert_eq!(owner.layers(), saved.3);
                assert_eq!(thermal.cells(), &saved.4);
                assert_eq!(thermal.time_s(), 2.0 * dt);
                assert_eq!(owner.mechanical().accepted_steps(), 2);
                final_motion.push((owner.mechanical().pose().positions_m[1][0], owner.mechanical().body_spins_kg_m2_s()[1]));
            });
        });
    }
    assert!((final_motion[0].0 - final_motion[1].0).abs() > 1e-12);
    for axis in 0..3 {
        assert!((final_motion[0].1[axis] - final_motion[1].1[axis]).abs() > 1e-12);
    }
}
