//! Native inertia, finite gas and complete contact heat share one endpoint.
//! The independent scalar oracle solves a constrained P1 tetrahedron; all
//! source values are manufactured and no production solver is reused.
use super::*;
use fs_couple::solid_enthalpy::plastic::dynamics::ThermoplasticDynamics;

const STEP: f64 = 2e-6;
const HEAT_CAPACITY: f64 = 5.0;
const CONDUCTIVITY: f64 = 400.0;
const VISCOSITY: f64 = 2000.0;
const INTERFACE_RESISTANCE: f64 = 5e-6;
const SHEAR: f64 = E / 2.4;
const CONSTRAINED_MODULUS: f64 = E * 0.8 / (1.2 * 0.6);

fn moving_fixture(density: f64) -> Fixture {
    let mut claims = ClaimSet::new();
    for claim in [
        scalar(YOUNG_MODULUS_PROPERTY, Pressure::DIMS, E),
        scalar(POISSON_RATIO_PROPERTY, Dims([0; 6]), 0.2),
        affine(YIELD_STRESS_PROPERTY, Pressure::DIMS, YIELD, -SOFTENING),
        scalar(BASIS.property_name(), Pressure::DIMS, H),
        affine(
            OVERSTRESS_VISCOSITY_PROPERTY,
            OVERSTRESS_VISCOSITY_DIMS,
            VISCOSITY,
            -20.0,
        ),
        scalar(DENSITY_PROPERTY, Density::DIMS, density),
        scalar(CP, CP_DIMS, HEAT_CAPACITY),
        scalar(
            "thermal_conductivity",
            THERMAL_CONDUCTIVITY_DIMS,
            CONDUCTIVITY,
        ),
    ] {
        claims.insert_claim(claim).unwrap();
    }
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "manufactured inertial gas-cooled member".into(),
            phase: "solid".into(),
            process: "independent constrained-tet reference".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    let lower = QueryPoint::new().with("T", 280.0).unwrap();
    let upper = QueryPoint::new().with("T", 320.0).unwrap();
    let caloric = ResolvedHeatCapacityCurve::try_new(
        &card,
        HeatCapacitySourceQuery {
            heat_capacity: &PropertyKey::new(CP, CP_DIMS),
            density: &PropertyKey::new(DENSITY_PROPERTY, Density::DIMS),
            temperature_axis: "T",
            lower: &lower,
            upper: &upper,
            selection: SelectionPolicy::SingleClaimOnly,
        },
        SolidLiquidPhase::Solid,
        -20.0 * HEAT_CAPACITY,
        1e-5,
    )
    .unwrap();
    Fixture {
        card,
        caloric,
        expansion: false,
    }
}

#[derive(Clone, Copy)]
struct MovingReference {
    strain: f64,
    plastic: f64,
    alpha: f64,
    velocity: f64,
    temperature: f64,
    passive: f64,
    gas: f64,
}
impl MovingReference {
    fn initial() -> Self {
        Self {
            strain: 0.0,
            plastic: 0.0,
            alpha: 0.0,
            velocity: 30.0,
            temperature: 300.0,
            passive: 300.0,
            gas: 310.0,
        }
    }
    fn free_energy(self) -> f64 {
        0.5 * CONSTRAINED_MODULUS * self.strain * self.strain
            - 2.0 * SHEAR * self.strain * self.plastic
            + 1.5 * SHEAR * self.plastic * self.plastic
            + 0.5 * H * self.alpha * self.alpha
    }
    fn stress(self) -> f64 {
        CONSTRAINED_MODULUS * self.strain - 2.0 * SHEAR * self.plastic
    }
    fn next_at_temperature(self, t: f64, density: f64, g_contact: f64) -> Self {
        let mass = density * V / 10.0;
        let eta = VISCOSITY - 20.0 * (t - 300.0);
        let state = |strain: f64| {
            let trial = 2.0 * SHEAR * strain - 3.0 * SHEAR * self.plastic;
            let da = ((trial.abs() - (YIELD - SOFTENING * (t - 300.0)) - H * self.alpha)
                / (3.0 * SHEAR + H + eta / STEP))
                .max(0.0);
            let plastic = self.plastic + trial.signum() * da;
            let stress = CONSTRAINED_MODULUS * strain - 2.0 * SHEAR * plastic;
            let residual = mass * (L * (strain - self.strain) - STEP * self.velocity)
                / (STEP * STEP)
                + stress * V / L;
            (residual, plastic, self.alpha + da)
        };
        let (mut lo, mut hi) = (-0.03, 0.03);
        assert!(state(lo).0 < 0.0 && state(hi).0 > 0.0);
        for _ in 0..80 {
            let mid = 0.5 * (lo + hi);
            if state(mid).0 > 0.0 {
                hi = mid;
            } else {
                lo = mid;
            }
        }
        let strain = 0.5 * (lo + hi);
        let (_, plastic, alpha) = state(strain);
        let capacity = density * V * HEAT_CAPACITY;
        let g_gas = AREA / (1e-4 / 0.03 + 0.005 / CONDUCTIVITY);
        Self {
            strain,
            plastic,
            alpha,
            velocity: L * (strain - self.strain) / STEP,
            temperature: t,
            passive: (capacity * self.passive + STEP * g_contact * t)
                / (capacity + STEP * g_contact),
            gas: (GAS_CAPACITY * self.gas + STEP * g_gas * t) / (GAS_CAPACITY + STEP * g_gas),
        }
    }
    fn physical_heat_from(self, old: Self) -> f64 {
        let da = self.alpha - old.alpha;
        V * ((YIELD - SOFTENING * (self.temperature - 300.0)) * da
            + (VISCOSITY - 20.0 * (self.temperature - 300.0)) / STEP * da * da)
    }
    fn next(self, density: f64, g_contact: f64) -> Self {
        let capacity = density * V * HEAT_CAPACITY;
        let residual = |t| {
            let new = self.next_at_temperature(t, density, g_contact);
            capacity * (new.temperature - self.temperature)
                + capacity * (new.passive - self.passive)
                + GAS_CAPACITY * (new.gas - self.gas)
                - new.physical_heat_from(self)
        };
        let (mut lo, mut hi) = (280.0, 320.0);
        assert!(residual(lo) < 0.0 && residual(hi) > 0.0);
        for _ in 0..80 {
            let mid = 0.5 * (lo + hi);
            if residual(mid) > 0.0 {
                hi = mid;
            } else {
                lo = mid;
            }
        }
        self.next_at_temperature(0.5 * (lo + hi), density, g_contact)
    }
}

#[test]
fn finite_gas_complete_heat_continues_native_inertia_and_material_history() {
    context(|cx| {
        let fixed: Vec<_> = (0..12).filter(|d| *d != 3).map(|d| (d, 0.0)).collect();
        let mut responses = Vec::new();
        // The second run changes only source density. Both the reference mass
        // and each original caloric inventory must follow that material data.
        for density in [RHO, 2.0 * RHO] {
            let fixture = moving_fixture(density);
            let mut thermal = initial(&[&fixture, &fixture], &[(0, 1, 999.0)], &[0.0; 2]);
            let bound = bind(&thermal, &[&fixture], &NODES, &TETS, cx)
                .with_linear_overstress(&[Some(MaterialPropertySelection::SingleClaimOnly)], cx)
                .unwrap();
            let mut velocity = vec![[0.0; 3]; 4];
            velocity[1][0] = 30.0;
            let mut moving = ThermoplasticDynamics::try_new(bound.clone(), velocity, cx).unwrap();
            let mut gas = gas();
            let conductivity = PropertyKey::new("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS);
            let selections = [Some((
                conductivity.clone(),
                MaterialPropertySelection::SingleClaimOnly,
            ))];
            let mut cells = bound.thermal_material_cells(&selections).unwrap();
            let point = cells[0].point;
            cells.push(MaterialCell {
                cell: CellId(101),
                material_id: MaterialId(1),
                curve: fixture.caloric.curve(),
                point,
                conductivity: Some((&conductivity, MaterialPropertySelection::SingleClaimOnly)),
            });
            let links = [ConductionLink {
                link: 0,
                distance_a_m: 0.005,
                distance_b_m: 0.005,
            }];
            let bulk = MaterialNetwork {
                cells: &cells,
                resistors: &[],
                links: &links,
            };
            let electrical = ElectricalContactNetwork {
                bulk: &bulk,
                contacts: &[],
                electrical_bases: &[],
            };
            let mut claims = ClaimSet::new();
            claims
                .insert_claim(scalar(CONTACT_R, CONTACT_DIMS, INTERFACE_RESISTANCE))
                .unwrap();
            let interface = InterfaceSystemCard::assemble(
                SurfaceSpec {
                    material: fixture.card.id().clone(),
                    texture_frame: "a".into(),
                },
                SurfaceSpec {
                    material: fixture.card.id().clone(),
                    texture_frame: "b".into(),
                },
                SystemContext {
                    medium: "synthetic fixed contact".into(),
                    third_body: None,
                    environment: "test".into(),
                    history: "fixed".into(),
                },
                claims,
                vec![],
            )
            .unwrap();
            let point = QueryPoint::new()
                .with("Ta", 300.0)
                .unwrap()
                .with("Tb", 300.0)
                .unwrap();
            let contacts = [MaterialContactLink {
                link: 0,
                sides: [0, 1].map(|i| ContactCell {
                    cell: CellId(100 + i),
                    material_id: MaterialId(i),
                    curve: fixture.caloric.curve(),
                }),
                point: ThermalInterfacePoint {
                    interface: &interface,
                    materials: [&fixture.card, &fixture.card],
                    context: &point,
                    temperature_axes: ["Ta", "Tb"],
                    liquid_fraction_axes: [None, None],
                    mean_temperature_axis: None,
                },
                selection: ClaimSelection::Policy(SelectionPolicy::SingleClaimOnly),
            }];
            let heat = MaterialHeatNetwork {
                electrical: &electrical,
                thermal_contacts: &contacts,
                radiation: None,
            };
            let prepared = heat.prepare(cx, &thermal).unwrap();
            let g_contact = 0.5 * L * L / (0.01 / CONDUCTIVITY + INTERFACE_RESISTANCE);
            let mut reference = MovingReference::initial();
            let mut solve = controls();
            solve.maximum_cells = 2;
            solve.newton.linear_restart = 3;
            solve.newton.absolute_tolerance = 1e-12;
            solve.difference_step_j = 1e-8;
            solve.cell_energy_tolerance_j = 1e-11;
            solve.thermal_energy_tolerance_j = 1e-11;
            solve.gas_energy_tolerance_j = 1e-11;
            for step in 1..=3 {
                let expected = reference.next(density, g_contact);
                let old_solution = moving.solid().solution().clone();
                let old_velocity = moving.velocity_mps().to_vec();
                let old_thermal = thermal.cells().clone();
                let candidate = moving
                    .advance_with_isochoric_gas_heat_network(
                        &thermal,
                        &gas,
                        step as f64 * STEP,
                        0.0,
                        &[contact()],
                        &prepared,
                        &fixed,
                        &[],
                        IncrementSettings::default(),
                        1e-3,
                        solve,
                        &[0.0; 3],
                        cx,
                    )
                    .unwrap();
                let report = candidate.report();
                let material = &report.material;
                let inertia = report.inertia;
                let solution = candidate.dynamics().solid().solution();
                close(solution.displacement_m[1][0], L * expected.strain, 3e-11);
                close(
                    candidate.dynamics().velocity_mps()[1][0],
                    expected.velocity,
                    3e-5,
                );
                close(solution.element_states[0].alpha, expected.alpha, 3e-9);
                close(
                    candidate.thermal().temperature_k(CellId(100)).unwrap(),
                    expected.temperature,
                    3e-6,
                );
                close(
                    candidate.thermal().temperature_k(CellId(101)).unwrap(),
                    expected.passive,
                    3e-6,
                );
                close(
                    candidate.gas().state().state().temperature,
                    expected.gas,
                    3e-6,
                );
                let mass = density * V / 10.0;
                let before = 0.5 * mass * reference.velocity * reference.velocity;
                let after = 0.5 * mass * expected.velocity * expected.velocity;
                let inertial_loss = 0.5 * mass * (expected.velocity - reference.velocity).powi(2);
                let stored = V * (expected.free_energy() - reference.free_energy());
                let physical = expected.physical_heat_from(reference);
                let material_loss = V * expected.stress() * (expected.strain - reference.strain)
                    - stored
                    - physical;
                close(inertia.reference_mass_kg, density * V, 1e-18);
                close(inertia.kinetic_before_j, before, 1e-8);
                close(inertia.kinetic_after_j, after, 1e-8);
                close(inertia.kinetic_change_j, after - before, 1e-8);
                close(inertia.numerical_loss_j, inertial_loss, 1e-8);
                close(material.physical_heat_j, physical, 1e-8);
                close(material.stored_mechanical_change_j, stored, 1e-8);
                close(material.numerical_loss_j, material_loss, 1e-8);
                assert!(inertia.numerical_loss_j > 0.0);
                assert!(material.numerical_loss_j > 0.0);
                close(material.reversible_heat_j, 0.0, 1e-12);
                close(material.mechanical_work_j, 0.0, 1e-12);
                close(material.thermal_gas.combined_energy_residual_j, 0.0, 1e-8);
                // Independent caloric and mechanical identities distinguish
                // actual material heat from BOTH strictly positive BE losses.
                let capacity = density * V * HEAT_CAPACITY;
                let caloric = capacity * (expected.temperature - reference.temperature)
                    + capacity * (expected.passive - reference.passive)
                    + GAS_CAPACITY * (expected.gas - reference.gas);
                close(caloric, physical, 1e-10);
                close(
                    caloric + stored + after - before + material_loss + inertial_loss,
                    0.0,
                    1e-10,
                );
                let source = material.thermal_gas.transport.as_ref().unwrap();
                assert_eq!(source.thermal_contacts.len(), 1);
                let contact = &source.thermal_contacts[0];
                close(contact.bulk_conductance_w_k, 2.0, 1e-12);
                close(contact.conductance_w_k, g_contact, 1e-12);
                close(
                    contact.phases[0].temperature_k(),
                    expected.temperature,
                    3e-6,
                );
                close(contact.phases[1].temperature_k(), expected.passive, 3e-6);
                interface
                    .claims()
                    .verify_receipt(contact.interface.receipt())
                    .unwrap();
                assert_eq!(moving.solid().solution(), &old_solution);
                assert_eq!(moving.velocity_mps(), old_velocity);
                assert_eq!(thermal.cells(), &old_thermal);
                close(gas.state().state().temperature, reference.gas, 3e-6);
                let (next_thermal, next_gas, next_motion, _) = candidate.into_parts();
                thermal = next_thermal;
                gas = next_gas;
                moving = next_motion;
                reference = expected;
            }
            assert!(reference.alpha > 0.0);
            assert!(reference.passive > 300.0);
            responses.push(reference);
        }
        assert!(responses[1].strain > responses[0].strain);
        assert!(responses[1].alpha > responses[0].alpha);
        assert!((responses[1].temperature - responses[0].temperature).abs() > 0.05);
        assert!((responses[1].gas - responses[0].gas).abs() > 1e-8);
    });
}

#[test]
fn bulk_finite_gas_keeps_rigid_velocity_and_refusals_leave_all_inputs_unchanged() {
    context(|cx| {
        let fixture = moving_fixture(RHO);
        let thermal = initial(&[&fixture], &[], &[0.0]);
        let gas = gas();
        let moving = ThermoplasticDynamics::try_new(
            bind(&thermal, &[&fixture], &NODES, &TETS, cx),
            vec![[0.1, 0.0, 0.0]; 4],
            cx,
        )
        .unwrap();
        let old_solution = moving.solid().solution().clone();
        let old_thermal = thermal.cells().clone();
        let duration = 0.001;
        let advance = |external_heat, cx: &Cx<'_>| {
            moving.advance_with_isochoric_gas(
                &thermal,
                &gas,
                duration,
                external_heat,
                &[contact()],
                None,
                &[],
                &[],
                IncrementSettings::default(),
                1e-3,
                controls(),
                &[0.0; 2],
                cx,
            )
        };
        let candidate = advance(0.0, cx).unwrap();
        let capacity = RHO * V * HEAT_CAPACITY;
        let g = AREA / (1e-4 / 0.03 + 0.005 / CONDUCTIVITY);
        let heat =
            duration * g * 10.0 / (1.0 + duration * g * (capacity.recip() + GAS_CAPACITY.recip()));
        close(
            candidate.thermal().temperature_k(CellId(100)).unwrap(),
            300.0 + heat / capacity,
            3e-6,
        );
        close(
            candidate.gas().state().state().temperature,
            310.0 - heat / GAS_CAPACITY,
            3e-6,
        );
        for (u, v) in candidate
            .dynamics()
            .solid()
            .solution()
            .displacement_m
            .iter()
            .zip(candidate.dynamics().velocity_mps())
        {
            close(u[0], 0.1 * duration, 1e-10);
            close(v[0], 0.1, 1e-9);
        }
        close(candidate.report().inertia.kinetic_change_j, 0.0, 1e-10);
        close(candidate.report().inertia.numerical_loss_j, 0.0, 1e-10);
        close(candidate.report().material.physical_heat_j, 0.0, 1e-10);
        close(
            candidate
                .report()
                .material
                .thermal_gas
                .combined_energy_residual_j,
            0.0,
            1e-8,
        );
        assert!(advance(10.0, cx).is_err());
        let gate = CancelGate::new();
        gate.request();
        fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default()).scope(|arena| {
            let cancelled = Cx::new(
                &gate,
                arena,
                StreamKey {
                    seed: 1,
                    kernel_id: 2,
                    tile: 0,
                    iteration: 0,
                },
                Budget::INFINITE,
                ExecMode::Deterministic,
            );
            assert!(advance(0.0, &cancelled).is_err());
        });
        assert_eq!(moving.solid().solution(), &old_solution);
        assert_eq!(moving.velocity_mps(), &[[0.1, 0.0, 0.0]; 4]);
        assert_eq!(thermal.cells(), &old_thermal);
        assert_eq!(moving.time_s(), 0.0);
        assert_eq!(thermal.time_s(), 0.0);
        close(gas.state().state().temperature, 310.0, 1e-12);
    });
}
