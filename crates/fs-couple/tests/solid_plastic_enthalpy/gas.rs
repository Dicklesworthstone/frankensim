//! Finite gas changes actual cyclic creep; mechanical heat returns to that gas.
mod dynamics;

use super::*;
use fs_conduction::transient::deposition::contact::{ContactCell, MaterialContactLink};
use fs_conduction::{
    AREA_SPECIFIC_THERMAL_RESISTANCE_DIMS as CONTACT_DIMS,
    AREA_SPECIFIC_THERMAL_RESISTANCE_PROPERTY as CONTACT_R,
};
use fs_couple::circuit_enthalpy::contacts::heat::MaterialHeatNetwork;
use fs_couple::circuit_enthalpy::contacts::material::ElectricalContactNetwork;
use fs_couple::circuit_enthalpy::material::{ConductionLink, MaterialCell, MaterialNetwork};
use fs_couple::solid_enthalpy::gas::{
    GasConductionGeometry, SolidGasContact, SolidGasError, SolidGasSolve,
};
use fs_matdb::{ClaimSelection, InterfaceSystemCard, SurfaceSpec, SystemContext};
use fs_material::plastic::rate::{OVERSTRESS_VISCOSITY_DIMS, OVERSTRESS_VISCOSITY_PROPERTY};
use fs_material::plastic::source::KINEMATIC_HARDENING_MODULUS_PROPERTY;
use fs_material::state_point::THERMAL_CONDUCTIVITY_DIMS;
use fs_material::thermal::interface::ThermalInterfacePoint;
use fs_material::thermochemical_gas::{
    ThermochemicalGasModel,
    energy::{GasEnergySolve, IsochoricGasParcel},
};
use fs_qty::QuantitySpec;
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_solver::NewtonKrylovConfig;
use fs_thermochem::mixture::FrozenIdealGasMixtureModelV1;
use fs_thermochem::{Composition, CompositionBasis};

#[allow(dead_code)]
#[path = "../material_exchange/gas.rs"]
mod gas_fixture;

const DT: f64 = 0.1;
const HK: f64 = 2e8;
const HK_SLOPE: f64 = -5e6;
const CP_VALUE: f64 = 100.0;
const GAS_VOLUME: f64 = 2e-5;
const GAS_CAPACITY: f64 = 2.5 * 100000.0 * GAS_VOLUME / 310.0;
const AREA: f64 = 2e-4;

fn affine(name: &str, dims: Dims, value: f64, slope: f64) -> PropertyClaim {
    let mut claim = scalar(name, dims, value);
    claim.value = PropertyValue::Curve {
        abscissa: "T".into(),
        abscissa_dims: Temperature::DIMS,
        knots: [280.0, 300.0, 320.0]
            .map(|t| (t, value + slope * (t - 300.0)))
            .to_vec(),
        dims,
    };
    claim.interpolation = InterpolationPolicy::LinearInside;
    claim.validity = ValidityDomain::unconstrained().with("T", 280.0, 320.0);
    claim
}
fn fixture(k: f64) -> Fixture {
    let mut claims = ClaimSet::new();
    for claim in [
        scalar(YOUNG_MODULUS_PROPERTY, Pressure::DIMS, E),
        scalar(POISSON_RATIO_PROPERTY, Dims([0; 6]), 0.2),
        affine(YIELD_STRESS_PROPERTY, Pressure::DIMS, YIELD, -SOFTENING),
        scalar(BASIS.property_name(), Pressure::DIMS, H),
        affine(
            KINEMATIC_HARDENING_MODULUS_PROPERTY,
            Pressure::DIMS,
            HK,
            HK_SLOPE,
        ),
        affine(
            OVERSTRESS_VISCOSITY_PROPERTY,
            OVERSTRESS_VISCOSITY_DIMS,
            1e7,
            -1e5,
        ),
        scalar(DENSITY_PROPERTY, Density::DIMS, RHO),
        scalar(CP, CP_DIMS, CP_VALUE),
        scalar("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS, k),
    ] {
        claims.insert_claim(claim).unwrap();
    }
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "manufactured heated viscoplastic member".into(),
            phase: "solid".into(),
            process: "synthetic cyclic reference".into(),
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
        -20.0 * CP_VALUE,
        1e-5,
    )
    .unwrap();
    Fixture {
        card,
        caloric,
        expansion: false,
    }
}
fn gas() -> IsochoricGasParcel {
    let (species, transport) = gas_fixture::component("heated-creep-reference", 0.028, 0.0);
    let model = ThermochemicalGasModel::new(
        FrozenIdealGasMixtureModelV1::new(
            vec![species],
            Composition::new(CompositionBasis::MoleFraction, vec![1.0]).unwrap(),
        )
        .unwrap(),
        vec![transport],
        None,
    )
    .unwrap();
    let typed = |kind| QuantitySpec::semantic(SemanticType::new(kind, ValueForm::Static));
    let point = QueryPoint::new()
        .with_quantity(
            "temperature",
            typed(QuantityKind::AbsoluteTemperature),
            310.0,
        )
        .unwrap()
        .with_quantity("pressure", typed(QuantityKind::Pressure), 100000.0)
        .unwrap();
    IsochoricGasParcel::new(model, point, GAS_VOLUME).unwrap()
}
fn controls() -> SolidGasSolve {
    SolidGasSolve {
        maximum_cells: 1,
        maximum_contacts: 1,
        maximum_trials: 300,
        maximum_newton_iterations: 30,
        difference_step_j: 1e-6,
        cell_energy_tolerance_j: 1e-9,
        thermal_energy_tolerance_j: 1e-9,
        gas_energy_tolerance_j: 1e-9,
        total_energy_tolerance_j: 1e-8,
        maximum_gas_temperature_increment_k: 20.0,
        gas: GasEnergySolve {
            minimum_temperature_k: 280.0,
            maximum_temperature_k: 320.0,
            energy_tolerance_j_kg: 1e-9,
            maximum_iterations: 128,
        },
        newton: NewtonKrylovConfig {
            absolute_tolerance: 1e-10,
            relative_tolerance: 1e-14,
            forcing_minimum: 1e-4,
            linear_restart: 2,
            ..NewtonKrylovConfig::default()
        },
    }
}
fn contact() -> SolidGasContact {
    SolidGasContact {
        id: 7,
        cell: CellId(100),
        geometry: GasConductionGeometry {
            area_m2: AREA,
            gas_distance_m: 1e-4,
            cell_distance_m: 0.005,
        },
        conductivity: PropertyKey::new("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS),
        selection: MaterialPropertySelection::SingleClaimOnly,
    }
}

// Independent uniaxial viscous consistency and two finite-reservoir balances.
// Gas T is eliminated analytically; only the scalar solid T root remains.
#[derive(Clone, Copy)]
struct Reference {
    t: f64,
    tg: f64,
    sigma: f64,
    p: f64,
    a: f64,
    hk: f64,
}
impl Reference {
    fn initial() -> Self {
        Self {
            t: 300.0,
            tg: 310.0,
            sigma: 0.0,
            p: 0.0,
            a: 0.0,
            hk: HK,
        }
    }
    fn next(self, t: f64, sigma: f64, g: f64) -> Self {
        let hk = HK + HK_SLOPE * (t - 300.0);
        let relative = sigma - hk * self.p;
        let da = ((relative.abs() - (YIELD - SOFTENING * (t - 300.0)) - H * self.a)
            / (H + hk + (1e7 - 1e5 * (t - 300.0)) / DT))
            .max(0.0);
        Self {
            t,
            tg: (GAS_CAPACITY * self.tg + DT * g * t) / (GAS_CAPACITY + DT * g),
            sigma,
            p: self.p + relative.signum() * da,
            a: self.a + da,
            hk,
        }
    }
    fn entropy(self) -> f64 {
        -0.5 * HK_SLOPE * self.p * self.p
    }
    fn free_energy(self) -> f64 {
        0.5 * self.sigma * self.sigma / E
            + 0.5 * H * self.a * self.a
            + 0.5 * self.hk * self.p * self.p
    }
    fn strain(self) -> f64 {
        self.sigma / E + self.p
    }
}
fn energies(old: Reference, new: Reference) -> [f64; 4] {
    let da = new.a - old.a;
    let physical =
        (YIELD - SOFTENING * (new.t - 300.0)) * da + (1e7 - 1e5 * (new.t - 300.0)) / DT * da * da;
    let coefficient = 0.5 * (new.hk - old.hk) * old.p * old.p;
    let entropy_storage = new.t * new.entropy() - old.t * old.entropy();
    [
        physical,
        coefficient,
        entropy_storage,
        -coefficient - entropy_storage,
    ]
}
fn endpoint(old: Reference, sigma: f64, g: f64) -> Reference {
    let residual = |t| {
        let new = old.next(t, sigma, g);
        let [physical, _, _, reversible] = energies(old, new);
        RHO * V * CP_VALUE * (t - old.t) + GAS_CAPACITY * (new.tg - old.tg)
            - V * (physical + reversible)
    };
    let (mut lo, mut hi) = (old.t, 319.0);
    assert!(residual(lo) < 0.0 && residual(hi) > 0.0);
    for _ in 0..80 {
        let mid = 0.5 * (lo + hi);
        if residual(mid) > 0.0 {
            hi = mid;
        } else {
            lo = mid;
        }
    }
    old.next(0.5 * (lo + hi), sigma, g)
}

#[test]
fn finite_gas_and_material_transport_change_native_cyclic_creep_with_reciprocal_energy() {
    context(|cx| {
        let mut responses = Vec::new();
        for k in [1.0, 4.0] {
            let fixture = fixture(k);
            let mut thermal = initial(&[&fixture], &[], &[0.0]);
            let mut gas = gas();
            let mut query = fixture.query();
            query.hardening_basis = HardeningPropertyBasis::MixedPlasticStrainModuli;
            let mut solid = ThermoplasticSolid::try_new(
                &thermal,
                &NODES,
                &TETS,
                vec![ThermoplasticCell {
                    cell: CellId(100),
                    query,
                }],
                TetAssemblyBudget::standard(),
                cx,
            )
            .unwrap()
            .with_linear_overstress(&[Some(MaterialPropertySelection::SingleClaimOnly)], cx)
            .unwrap();
            let g = AREA / (1e-4 / 0.03 + 0.005 / k);
            let mut old = Reference::initial();
            for (i, sigma) in [3.4e6, -3.4e6, 3.4e6].into_iter().enumerate() {
                let original_solution = solid.solution().clone();
                let original_cells = thermal.cells().clone();
                let original_gas = gas.clone();
                let candidate = solid
                    .advance_with_isochoric_gas(
                        &thermal,
                        &gas,
                        (i + 1) as f64 * DT,
                        0.0,
                        &[contact()],
                        None,
                        &ROLLERS,
                        &[(3, sigma * V / L)],
                        IncrementSettings::default(),
                        1e-3,
                        controls(),
                        &[0.0; 2],
                        cx,
                    )
                    .unwrap();
                let expected = endpoint(old, sigma, g);
                close(
                    candidate.thermal().temperature_k(CellId(100)).unwrap(),
                    expected.t,
                    3e-6,
                );
                close(
                    candidate.gas().state().state().temperature,
                    expected.tg,
                    3e-6,
                );
                let solution = candidate.solid().solution();
                close(solution.element_states[0].alpha, expected.a, 3e-9);
                close(
                    solution.element_states[0].plastic_strain[0],
                    expected.p,
                    3e-9,
                );
                close(solution.element_stress[0][0], sigma, 3e-3);
                close(solution.displacement_m[1][0], L * expected.strain(), 3e-11);
                let report = candidate.report();
                let [physical, coefficient, ts, reversible] = energies(old, expected);
                close(report.physical_heat_j, V * physical, 1e-8);
                close(report.coefficient_storage_change_j, V * coefficient, 1e-8);
                close(report.thermoelastic_storage_change_j, V * ts, 1e-8);
                close(report.reversible_heat_j, V * reversible, 1e-8);
                close(
                    report.stored_mechanical_change_j,
                    V * (expected.free_energy() - old.free_energy()),
                    1e-8,
                );
                close(
                    report.mechanical_work_j,
                    V * sigma * (expected.strain() - old.strain()),
                    1e-8,
                );
                let (ds, dp, da) = (
                    expected.sigma - old.sigma,
                    expected.p - old.p,
                    expected.a - old.a,
                );
                close(
                    report.numerical_loss_j,
                    V * (0.5 * ds * ds / E + 0.5 * H * da * da + 0.5 * expected.hk * dp * dp),
                    1e-8,
                );
                close(
                    report.thermal_gas.solid_to_gas_heat_j,
                    GAS_CAPACITY * (expected.tg - old.tg),
                    1e-8,
                );
                close(report.thermal_gas.combined_energy_residual_j, 0.0, 1e-8);
                assert_eq!(report.thermal_gas.external_thermal_heat_j, 0.0);
                assert_eq!(report.thermal_gas.external_gas_heat_j, 0.0);
                let patch = &report.thermal_gas.contacts[0];
                close(patch.conductance_w_k, g, 1e-12);
                let claim = patch.source.property("thermal_conductivity").unwrap();
                fixture
                    .card
                    .claims()
                    .verify_receipt(&claim.answer().receipt)
                    .unwrap();
                assert_eq!(candidate.gas().mass_kg(), gas.mass_kg());
                assert_eq!(candidate.thermal().cells()[&CellId(100)].mass_kg, RHO * V);
                assert_eq!(solid.solution(), &original_solution);
                assert_eq!(thermal.cells(), &original_cells);
                assert_eq!(gas.state(), original_gas.state());
                // A failed joint Newton trial must not consume creep time or gas heat.
                let limited = SolidGasSolve {
                    maximum_trials: 1,
                    ..controls()
                };
                assert!(matches!(
                    solid.advance_with_isochoric_gas(
                        &thermal,
                        &gas,
                        (i + 1) as f64 * DT,
                        0.0,
                        &[contact()],
                        None,
                        &ROLLERS,
                        &[(3, sigma * V / L)],
                        IncrementSettings::default(),
                        1e-3,
                        limited,
                        &[0.0; 2],
                        cx
                    ),
                    Err(SolidGasError::TrialBudget)
                ));
                assert_eq!(solid.solution(), &original_solution);
                assert_eq!(thermal.cells(), &original_cells);
                assert_eq!(gas.state(), original_gas.state());
                let (next_thermal, next_gas, next_solid, _) = candidate.into_parts();
                thermal = next_thermal;
                gas = next_gas;
                solid = next_solid;
                old = expected;
            }
            responses.push(old);
        }
        assert!(responses[1].t - responses[0].t > 0.8);
        assert!(responses[1].a - responses[0].a > 5e-4);
        assert!(responses[1].p - responses[0].p > 7e-5);
    });
}

#[test]
fn complete_heat_contact_changes_loaded_creep_and_closes_three_finite_inventories() {
    context(|cx| {
        let fixture = fixture(4.0);
        let thermal = initial(&[&fixture, &fixture], &[(0, 1, 999.0)], &[0.0; 2]);
        let gas = gas();
        let mut query = fixture.query();
        query.hardening_basis = HardeningPropertyBasis::MixedPlasticStrainModuli;
        let solid = ThermoplasticSolid::try_new(
            &thermal,
            &NODES,
            &TETS,
            vec![ThermoplasticCell {
                cell: CellId(100),
                query,
            }],
            TetAssemblyBudget::standard(),
            cx,
        )
        .unwrap()
        .with_linear_overstress(&[Some(MaterialPropertySelection::SingleClaimOnly)], cx)
        .unwrap();
        let conductivity = PropertyKey::new("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS);
        let selections = [Some((
            conductivity.clone(),
            MaterialPropertySelection::SingleClaimOnly,
        ))];
        let mut cells = solid.thermal_material_cells(&selections).unwrap();
        let point = cells[0].point;
        // The second original caloric cell is a finite, unloaded heat inventory.
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
        let sigma = 3.4e6;
        let capacity = RHO * V * CP_VALUE;
        let g_gas = AREA / (1e-4 / 0.03 + 0.005 / 4.0);
        let old = Reference::initial();
        let mut responses = Vec::new();
        for resistance in [0.0005, 0.05] {
            let mut claims = ClaimSet::new();
            claims
                .insert_claim(scalar(CONTACT_R, CONTACT_DIMS, resistance))
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
            let mut solve = controls();
            solve.maximum_cells = 2;
            solve.newton.linear_restart = 3;
            let candidate = solid
                .advance_with_isochoric_gas_heat_network(
                    &thermal,
                    &gas,
                    DT,
                    0.0,
                    &[contact()],
                    &prepared,
                    &ROLLERS,
                    &[(3, sigma * V / L)],
                    IncrementSettings::default(),
                    1e-3,
                    solve,
                    &[0.0; 3],
                    cx,
                )
                .unwrap();
            // Independent series law and finite-reservoir elimination. The
            // authored 999 W/K is deliberately incompatible with either case.
            let g_contact = 0.5 * L * L / (0.005 / 4.0 + resistance + 0.005 / 4.0);
            let passive = |t| (capacity * 300.0 + DT * g_contact * t) / (capacity + DT * g_contact);
            let residual = |t| {
                let new = old.next(t, sigma, g_gas);
                let [physical, _, _, reversible] = energies(old, new);
                capacity * (t - old.t)
                    + capacity * (passive(t) - 300.0)
                    + GAS_CAPACITY * (new.tg - old.tg)
                    - V * (physical + reversible)
            };
            let (mut lo, mut hi) = (300.0, 319.0);
            assert!(residual(lo) < 0.0 && residual(hi) > 0.0);
            for _ in 0..80 {
                let mid = 0.5 * (lo + hi);
                if residual(mid) > 0.0 {
                    hi = mid;
                } else {
                    lo = mid;
                }
            }
            let expected = old.next(0.5 * (lo + hi), sigma, g_gas);
            close(
                candidate.thermal().temperature_k(CellId(100)).unwrap(),
                expected.t,
                3e-6,
            );
            close(
                candidate.thermal().temperature_k(CellId(101)).unwrap(),
                passive(expected.t),
                3e-6,
            );
            close(
                candidate.gas().state().state().temperature,
                expected.tg,
                3e-6,
            );
            close(
                candidate.solid().solution().element_states[0].alpha,
                expected.a,
                3e-9,
            );
            close(
                candidate.solid().solution().displacement_m[1][0],
                L * expected.strain(),
                3e-11,
            );
            let report = candidate.report();
            let source = report.thermal_gas.transport.as_ref().unwrap();
            assert_eq!(source.thermal_contacts.len(), 1);
            let contact = &source.thermal_contacts[0];
            close(contact.bulk_conductance_w_k, 0.02, 1e-12);
            close(contact.conductance_w_k, g_contact, 1e-12);
            close(contact.phases[0].temperature_k(), expected.t, 3e-6);
            close(contact.phases[1].temperature_k(), passive(expected.t), 3e-6);
            interface
                .claims()
                .verify_receipt(contact.interface.receipt())
                .unwrap();
            let [physical, _, ts, reversible] = energies(old, expected);
            close(report.physical_heat_j, V * physical, 1e-8);
            close(report.reversible_heat_j, V * reversible, 1e-8);
            close(report.thermoelastic_storage_change_j, V * ts, 1e-8);
            close(
                report.thermal_gas.solid_to_gas_heat_j,
                GAS_CAPACITY * (expected.tg - old.tg),
                1e-8,
            );
            close(report.thermal_gas.combined_energy_residual_j, 0.0, 1e-8);
            assert_eq!(report.thermal_gas.external_thermal_heat_j, 0.0);
            assert_eq!(report.thermal_gas.closed_radiation_net_j, 0.0);
            responses.push(expected);
        }
        // Changing only the interface card alters the actual return-map history.
        assert!(responses[1].t - responses[0].t > 0.05);
        assert!(responses[1].a - responses[0].a > 1e-5);
        assert_eq!(solid.time_s(), 0.0);
        assert_eq!(thermal.time_s(), 0.0);
        close(gas.state().state().temperature, 310.0, 1e-12);
    });
}
