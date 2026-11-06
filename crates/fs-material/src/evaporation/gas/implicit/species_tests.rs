//! Distinct pure liquids share real multicomponent NASA9 gas mass and energy.
use super::*;
use crate::evaporation::gas::implicit::adaptive::{
    AdaptiveImplicitGasControls, advance_deposition_with_gas_adaptive,
};
use crate::evaporation::gas::implicit::timeline::advance_deposition_with_gas_timeline;
use crate::thermochemical_gas::tests::{point, species_with_energy_reference, transport};
use fs_thermal::evaporation::{EvaporationError, EvaporationLaw, EvaporationProperties};
use fs_thermochem::{Composition, CompositionBasis, mixture::FrozenIdealGasMixtureModelV1};

const AREA: f64 = 0.001;
const VOLUME: f64 = 0.001;
const DT: f64 = 20.0;

#[derive(Clone, Copy)]
struct PureLiquid {
    name: &'static str,
    molar_mass: f64,
    cp_r: f64,
    slope: f64,
    latent: f64,
    saturation_at_300: f64,
}
static LIQUIDS: [PureLiquid; 2] = [
    PureLiquid {
        name: "B",
        molar_mass: 0.018015,
        cp_r: 4.0,
        slope: 0.001,
        latent: 200_000.0,
        saturation_at_300: 1500.0,
    },
    PureLiquid {
        name: "C",
        molar_mass: 0.03204,
        cp_r: 5.0,
        slope: 0.002,
        latent: 300_000.0,
        saturation_at_300: 2000.0,
    },
];
impl PureLiquid {
    fn id(&self) -> SpeciesId {
        SpeciesId::new(self.name).unwrap()
    }
    fn r(&self) -> f64 {
        fs_thermochem::UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K / self.molar_mass
    }
    fn vapor_h(&self, t: f64) -> f64 {
        self.r() * (self.cp_r * t + 0.5 * self.slope * t * t)
    }
    fn h(&self, t: f64) -> f64 {
        self.vapor_h(t) - self.latent
    }
    fn u(&self, t: f64) -> f64 {
        self.vapor_h(t) - self.r() * t
    }
    fn dh(&self, a: f64, b: f64) -> f64 {
        self.r() * (b - a) * (self.cp_r + 0.5 * self.slope * (a + b))
    }
    fn du(&self, a: f64, b: f64) -> f64 {
        self.dh(a, b) - self.r() * (b - a)
    }
    fn psat(&self, t: f64) -> f64 {
        self.saturation_at_300 * fs_math::det::exp(self.latent / self.r() * (1.0 / 300.0 - 1.0 / t))
    }
    fn unit_rate(&self, t: f64, tg: f64, pressure: f64, partial: f64) -> f64 {
        pressure / (self.r() * tg) * ((pressure - partial) / (pressure - self.psat(t))).ln()
    }
}
// These are adapters of the same analytic NASA9 coefficients used by the real
// gas source. No gas EOS, mixture calorics, or production phase owner is mocked.
impl ThermalLaw for PureLiquid {
    fn temperature_k(&self, h: f64) -> Result<f64, DepositionError> {
        let y = (h + self.latent) / self.r();
        let t = 2.0 * y / (self.cp_r + (self.cp_r * self.cp_r + 2.0 * self.slope * y).sqrt());
        if (200.0..=900.0).contains(&t) {
            Ok(t)
        } else {
            Err(DepositionError::MaterialDomain)
        }
    }
    fn maximum_temperature_slope(&self) -> f64 {
        1.0 / (self.r() * (self.cp_r + 200.0 * self.slope))
    }
}
impl EvaporationLaw for PureLiquid {
    fn temperature_range_k(&self) -> [f64; 2] {
        [200.0, 900.0]
    }
    fn properties(&self, t: f64) -> Result<EvaporationProperties, EvaporationError> {
        let p = self.psat(t);
        Ok(EvaporationProperties {
            saturation_pressure_pa: p,
            saturation_pressure_derivative_pa_k: p * self.latent / (self.r() * t * t),
            latent_heat_j_kg: self.latent,
            liquid_density_kg_m3: 1000.0,
            molar_mass_kg_mol: self.molar_mass,
        })
    }
}

fn sources(references: [f64; 2]) -> (ThermochemicalGasModel, [GasLiquidBinding; 2]) {
    let catalog = ThermochemicalGasModel::new(
        FrozenIdealGasMixtureModelV1::new(
            vec![
                species_with_energy_reference("A", 0.02897, 3.5, 0.0, 0.0),
                species_with_energy_reference(
                    "B",
                    LIQUIDS[0].molar_mass,
                    LIQUIDS[0].cp_r,
                    LIQUIDS[0].slope,
                    references[0],
                ),
                species_with_energy_reference(
                    "C",
                    LIQUIDS[1].molar_mass,
                    LIQUIDS[1].cp_r,
                    LIQUIDS[1].slope,
                    references[1],
                ),
            ],
            Composition::new(CompositionBasis::MoleFraction, vec![0.9, 0.05, 0.05]).unwrap(),
        )
        .unwrap(),
        vec![
            transport("A", 1.7e-5, 0.025),
            transport("B", 2e-5, 0.03),
            transport("C", 1.5e-5, 0.02),
        ],
        Some(&LIQUIDS[0].id()),
    )
    .unwrap();
    let bindings = std::array::from_fn(|i| GasLiquidBinding {
        species: LIQUIDS[i].id(),
        liquid_enthalpy_offset_j_kg: references[i],
    });
    (catalog, bindings)
}
fn mixture_gas(
    catalog: &ThermochemicalGasModel,
    temperature: f64,
    partials: [f64; 2],
) -> ClosedGasParcel {
    let [b, c] = partials.map(|p| p / 100_000.0);
    ClosedGasParcel::new(
        catalog
            .with_composition(
                Composition::new(CompositionBasis::MoleFraction, vec![1.0 - b - c, b, c]).unwrap(),
            )
            .unwrap(),
        point(temperature, 100_000.0)
            .with("local-lot", 37.0)
            .unwrap(),
        VOLUME,
    )
    .unwrap()
}
fn pure_cell(id: u64, species: usize, temperature: f64, source_w: f64) -> DepositedCell {
    let mut cell = cell(id, temperature, source_w, AREA);
    cell.material_id = MaterialId(species as u64 + 1);
    cell.specific_h_j_kg = LIQUIDS[species].h(temperature);
    cell
}
fn scheduled(events: Vec<TimedEvent>) -> DepositionSimulation<'static> {
    DepositionSimulation::new(
        events,
        &[(MaterialId(1), &LIQUIDS[0]), (MaterialId(2), &LIQUIDS[1])],
        DepositionBudget {
            max_cells: 4,
            max_links: 2,
            max_events: 16,
            max_substeps: 1024,
        },
    )
    .unwrap()
}
fn events_at<'m>(thermal: &DepositionSimulation<'m>, time_s: f64) -> DepositionSimulation<'m> {
    thermal
        .advance_to(
            time_s,
            StepConfig {
                max_dt_s: 0.5,
                energy_tolerance_j: 1e-6,
            },
            || Ok(()),
        )
        .unwrap()
        .0
}
fn born(cells: Vec<DepositedCell>, contact: bool) -> DepositionSimulation<'static> {
    let mut events: Vec<_> = cells
        .into_iter()
        .map(|cell| TimedEvent {
            time_s: 0.0,
            event: Event::Birth(cell),
        })
        .collect();
    if contact {
        events.push(TimedEvent {
            time_s: 0.0,
            event: Event::Contact(ThermalLink {
                id: 9,
                a: CellId(1),
                b: CellId(2),
                area_m2: 0.0002,
                conductance_w_k: 2.0,
            }),
        });
    }
    events_at(&scheduled(events), 0.0)
}
fn bound<'a>(
    id: u64,
    species: usize,
    bindings: &'a [GasLiquidBinding; 2],
    mass_coefficient: f64,
    heat_coefficient: f64,
) -> ImplicitGasSurface<'a> {
    ImplicitGasSurface {
        cell_id: CellId(id),
        law: &LIQUIDS[species],
        binding: &bindings[species],
        pure_liquid_enthalpy_range_j_kg: [LIQUIDS[species].h(200.0), LIQUIDS[species].h(900.0)],
        mass_transfer_coefficient_m_s: mass_coefficient,
        heat_transfer_coefficient_w_m2_k: heat_coefficient,
    }
}
fn canonical_gas_change(initial: &ClosedGasParcel, final_gas: &ClosedGasParcel) -> f64 {
    let a = initial.state().state().temperature;
    let b = final_gas.state().state().temperature;
    initial.species_mass_kg(&carrier()).unwrap() * CARRIER_R * 2.5 * (b - a)
        + LIQUIDS
            .iter()
            .map(|law| {
                let old = initial.species_mass_kg(&law.id()).unwrap();
                let next = final_gas.species_mass_kg(&law.id()).unwrap();
                next * law.du(a, b) + (next - old) * law.u(a)
            })
            .sum::<f64>()
}

struct Manufactured {
    thermal: DepositionSimulation<'static>,
    gas: ClosedGasParcel,
    mass: f64,
    gas_temperature: f64,
    coefficients: [f64; 2],
    supplied_heat: [f64; 2],
    sensible_heat: [f64; 2],
    external_heat: f64,
}
fn manufactured(catalog: &ThermochemicalGasModel) -> Manufactured {
    let gas = mixture_gas(catalog, 320.0, [0.0, 5000.0]);
    let d = 2.0f64.powi(-20);
    let [b, c] = LIQUIDS;
    let carrier_mass = gas.species_mass_kg(&carrier()).unwrap();
    let final_mass = [d, gas.species_mass_kg(&c.id()).unwrap() - d];
    let k = [20.0, 10.0].map(|htc| htc * AREA * DT);
    // ΔU = d*h_B(360) - d*h_C(Tg) - Σ K_i*(Tg-T_i).
    // Unlike equal-opposite SAME-species streams, d*(R_B-R_C)*T0 remains.
    let linear = carrier_mass * CARRIER_R * 2.5
        + (0..2)
            .map(|i| {
                final_mass[i] * LIQUIDS[i].r() * (LIQUIDS[i].cp_r - 1.0 + LIQUIDS[i].slope * 320.0)
            })
            .sum::<f64>()
        + d * c.r() * (c.cp_r + c.slope * 320.0)
        + k[0]
        + k[1];
    let quadratic = 0.5
        * ((0..2)
            .map(|i| final_mass[i] * LIQUIDS[i].r() * LIQUIDS[i].slope)
            .sum::<f64>()
            + d * c.r() * c.slope);
    let rhs = d * ((b.r() - c.r()) * 320.0 + b.dh(320.0, 360.0)) + k[0] * 40.0 - k[1] * 40.0;
    let gas_temperature = 320.0 + quadratic_increment(linear, quadratic, rhs);
    let pressure = (carrier_mass * CARRIER_R + final_mass[0] * b.r() + final_mass[1] * c.r())
        * gas_temperature
        / VOLUME;
    let temperatures = [360.0, 280.0];
    let signed = [d, -d];
    let coefficients = std::array::from_fn(|i| {
        signed[i]
            / (AREA
                * DT
                * LIQUIDS[i].unit_rate(
                    temperatures[i],
                    gas_temperature,
                    pressure,
                    final_mass[i] * LIQUIDS[i].r() * gas_temperature / VOLUME,
                ))
    });
    assert!(coefficients.iter().all(|coefficient| *coefficient > 0.0));
    let supplied_heat = [
        LIQUID_MASS * b.dh(355.0, 360.0) + d * b.latent,
        LIQUID_MASS * c.dh(285.0, 280.0) - d * (c.latent + c.dh(280.0, gas_temperature)),
    ];
    let sensible_heat = std::array::from_fn(|i| k[i] * (gas_temperature - temperatures[i]));
    let contact_w = 2.0 * (temperatures[1] - temperatures[0]);
    let source_w = [
        (supplied_heat[0] - sensible_heat[0]) / DT - contact_w,
        (supplied_heat[1] - sensible_heat[1]) / DT + contact_w,
    ];
    let cells = [355.0, 285.0]
        .into_iter()
        .enumerate()
        .map(|(i, t)| {
            let mut cell = pure_cell(i as u64 + 1, i, t, source_w[i]);
            cell.geometry.surface_area_m2 += 0.0002;
            cell
        })
        .collect();
    Manufactured {
        thermal: born(cells, true),
        gas,
        mass: d,
        gas_temperature,
        coefficients,
        supplied_heat,
        sensible_heat,
        external_heat: DT * (source_w[0] + source_w[1]),
    }
}

#[test]
fn equal_opposite_distinct_vapors_retain_each_latent_gap_upstream_caloric_and_pressure_work() {
    let mut reference: Option<[f64; 5]> = None;
    for references in [[0.0, 0.0], [1e12, -5e11]] {
        let (catalog, bindings) = sources(references);
        let fixture = manufactured(&catalog);
        let mut config = controls();
        config.phase.latent_compatibility_tolerance_j_kg = 0.02;
        let mut previous = None;
        for order in [[0, 1], [1, 0]] {
            let surfaces = order.map(|i| {
                bound(
                    i as u64 + 1,
                    i,
                    &bindings,
                    fixture.coefficients[i],
                    [20.0, 10.0][i],
                )
            });
            let candidate = solve_deposition_with_gas(
                &fixture.thermal,
                &surfaces,
                &fixture.gas,
                &catalog,
                DT,
                config,
                |_, liquids, gas| Ok(coefficients((liquids.clone(), gas.species_inventory_kg()?))),
                || Ok(()),
            )
            .unwrap();
            let report = candidate.report();
            let tg = candidate.gas().state().state().temperature;
            near(tg, fixture.gas_temperature, 2e-8);
            for (i, temperature) in [360.0, 280.0].into_iter().enumerate() {
                let id = CellId(i as u64 + 1);
                let dm = [fixture.mass, -fixture.mass][i];
                let phase = &report.thermal.thermal.surfaces[&id];
                near(
                    candidate.thermal().temperature_k(id).unwrap(),
                    temperature,
                    2e-8,
                );
                near(
                    report.species_transferred_mass_kg[&LIQUIDS[i].id()],
                    dm,
                    2e-14,
                );
                near(
                    candidate.gas().species_mass_kg(&LIQUIDS[i].id()).unwrap()
                        - fixture.gas.species_mass_kg(&LIQUIDS[i].id()).unwrap(),
                    dm,
                    2e-14,
                );
                near(phase.liquid.mass_kg, LIQUID_MASS - dm, 2e-14);
                near(phase.latent_energy_j, dm * LIQUIDS[i].latent, 2e-8);
                near(
                    phase.vapor_energy_j,
                    dm * LIQUIDS[i].vapor_h(if dm < 0.0 { tg } else { temperature }),
                    2e-8,
                );
                near(phase.supplied_heat_j, fixture.supplied_heat[i], 2e-7);
                near(
                    report.sensible_heat_into_liquid_j[&id],
                    fixture.sensible_heat[i],
                    2e-8,
                );
                near(
                    candidate.payload().0[&id].mass_kg,
                    phase.liquid.mass_kg,
                    1e-14,
                );
            }
            near(candidate.gas().mass_kg(), fixture.gas.mass_kg(), 1e-14);
            near(
                candidate.gas().species_mass_kg(&carrier()).unwrap(),
                fixture.gas.species_mass_kg(&carrier()).unwrap(),
                1e-14,
            );
            near(
                report.species_transferred_mass_kg.values().sum(),
                0.0,
                2e-14,
            );
            let latent: f64 = report
                .thermal
                .thermal
                .surfaces
                .values()
                .map(|surface| surface.latent_energy_j)
                .sum();
            near(
                latent,
                fixture.mass * (LIQUIDS[0].latent - LIQUIDS[1].latent),
                2e-8,
            );
            assert!(
                latent.abs() > 0.09,
                "zero net mass does not cancel distinct latent heats"
            );
            near(
                candidate.thermal().energy_j() - fixture.thermal.energy_j()
                    + canonical_gas_change(&fixture.gas, candidate.gas()),
                fixture.external_heat,
                3e-7,
            );
            near(report.external_heat_j, fixture.external_heat, 2e-8);
            near(report.energy_residual_j, 0.0, 1e-7);
            assert_eq!(candidate.gas().point().axes().get("local-lot"), Some(&37.0));
            assert!(std::ptr::eq(candidate.catalog(), &catalog));
            let actual = [
                candidate.thermal().temperature_k(CellId(1)).unwrap(),
                candidate.thermal().temperature_k(CellId(2)).unwrap(),
                tg,
                candidate.gas().species_mass_kg(&LIQUIDS[0].id()).unwrap(),
                candidate.gas().species_mass_kg(&LIQUIDS[1].id()).unwrap(),
            ];
            if let Some(reference) = &reference {
                for i in 0..5 {
                    near(actual[i], reference[i], if i < 3 { 2e-8 } else { 1e-13 });
                }
            } else {
                reference = Some(actual);
            }
            let physical = (candidate.thermal().cells().clone(), candidate.gas().clone());
            if let Some(previous) = &previous {
                assert_eq!(&physical, previous);
            } else {
                previous = Some(physical);
            }
        }
    }
}

#[test]
fn two_absent_vapors_enter_from_their_own_liquid_without_seed_or_shared_latent_constant() {
    let (catalog, bindings) = sources([0.0, 0.0]);
    let gas = mixture_gas(&catalog, 300.0, [0.0, 0.0]);
    let thermal = born(
        vec![pure_cell(1, 0, 300.0, 0.0), pure_cell(2, 1, 300.0, 0.0)],
        false,
    );
    let surfaces = [
        bound(1, 0, &bindings, 0.02, 0.0),
        bound(2, 1, &bindings, 0.02, 0.0),
    ];
    let candidate = solve_deposition_with_gas(
        &thermal,
        &surfaces,
        &gas,
        &catalog,
        1.0,
        controls(),
        empty_heat,
        || Ok(()),
    )
    .unwrap();
    let tg = candidate.gas().state().state().temperature;
    let pressure = candidate.gas().state().state().pressure;
    for i in 0..2 {
        let id = CellId(i as u64 + 1);
        let law = LIQUIDS[i];
        let dm = candidate.report().species_transferred_mass_kg[&law.id()];
        assert_eq!(gas.species_mass_kg(&law.id()).unwrap(), 0.0);
        assert!(dm > 0.0);
        near(
            candidate.gas().species_mass_kg(&law.id()).unwrap(),
            dm,
            1e-14,
        );
        near(
            LIQUID_MASS - candidate.thermal().cells()[&id].mass_kg,
            dm,
            1e-14,
        );
        let t = candidate.thermal().temperature_k(id).unwrap();
        assert!(t < 300.0);
        near(
            dm,
            AREA * 0.02
                * law.unit_rate(
                    t,
                    tg,
                    pressure,
                    candidate
                        .gas()
                        .species_partial_pressure_pa(&law.id())
                        .unwrap(),
                ),
            1e-14,
        );
        near(
            candidate.report().thermal.thermal.surfaces[&id].latent_energy_j,
            dm * law.latent,
            1e-8,
        );
    }
    near(
        candidate.thermal().energy_j() - thermal.energy_j()
            + canonical_gas_change(&gas, candidate.gas()),
        0.0,
        1e-7,
    );
    near(
        candidate.gas().species_mass_kg(&carrier()).unwrap(),
        gas.species_mass_kg(&carrier()).unwrap(),
        1e-14,
    );
    assert_eq!(candidate.gas().point().axes().get("local-lot"), Some(&37.0));
}

#[test]
fn each_original_species_donor_bounds_gross_condensation_even_with_abundant_or_evaporating_other_vapor()
 {
    let (catalog, bindings) = sources([0.0, 0.0]);
    for (abundant_b, cold_patches) in [(true, 1usize), (true, 2), (false, 1)] {
        let gas = mixture_gas(
            &catalog,
            300.0,
            [if abundant_b { 30_000.0 } else { 0.0 }, 5000.0],
        );
        let gas_snapshot = gas.clone();
        let withdrawal = 0.36 * gas.species_mass_kg(&LIQUIDS[1].id()).unwrap();
        let cold_rate = -LIQUIDS[1].unit_rate(280.0, 300.0, 100_000.0, 5000.0);
        let cold_coefficient = withdrawal / (cold_patches as f64 * 10.0 * AREA * cold_rate);
        let hot_coefficient = if abundant_b {
            0.0
        } else {
            withdrawal / (10.0 * AREA * LIQUIDS[0].unit_rate(320.0, 300.0, 100_000.0, 0.0))
        };
        let mut cells = vec![pure_cell(1, 0, 320.0, 0.0)];
        let mut surfaces = vec![bound(1, 0, &bindings, hot_coefficient, 0.0)];
        for i in 0..cold_patches {
            cells.push(pure_cell(i as u64 + 2, 1, 280.0, 0.0));
            surfaces.push(bound(i as u64 + 2, 1, &bindings, cold_coefficient, 0.0));
        }
        let thermal = born(cells, false);
        let snapshot = thermal.clone();
        let calls = Cell::new(0usize);
        let error = solve_deposition_with_gas(
            &thermal,
            &surfaces,
            &gas,
            &catalog,
            10.0,
            controls(),
            |_, liquids, _| {
                calls.set(calls.get() + 1);
                let condensed: f64 = (0..cold_patches)
                    .map(|i| liquids[&CellId(i as u64 + 2)].mass_kg - LIQUID_MASS)
                    .sum();
                near(condensed, withdrawal, 1e-14);
                if !abundant_b {
                    near(LIQUID_MASS - liquids[&CellId(1)].mass_kg, condensed, 1e-14);
                }
                Ok(coefficients(()))
            },
            || Ok(()),
        )
        .err()
        .expect("the original C donor is insufficient under its declared fraction ceiling");
        assert!(matches!(
            error,
            ImplicitGasError::Phase(GasEvaporationError::StepSizeLimit { .. })
        ));
        assert_eq!(
            calls.get(),
            1,
            "the initial physical proposal must refuse before a later iterate can hide the donor violation"
        );
        unchanged(&thermal, &snapshot, &gas, &gas_snapshot);
    }
}

struct CountPhase<'a> {
    law: &'a PureLiquid,
    queries: Cell<usize>,
}
impl EvaporationLaw for CountPhase<'_> {
    fn temperature_range_k(&self) -> [f64; 2] {
        self.queries.set(self.queries.get() + 1);
        self.law.temperature_range_k()
    }
    fn properties(&self, t: f64) -> Result<EvaporationProperties, EvaporationError> {
        self.queries.set(self.queries.get() + 1);
        self.law.properties(t)
    }
}
fn adaptive_policy() -> AdaptiveImplicitGasControls {
    let mut solve = controls();
    solve.thermal.maximum_iterations = 10_000;
    solve.thermal.maximum_evaluations = 100_000;
    solve.thermal.newton.absolute_tolerance = 1e-10;
    let mut config = AdaptiveImplicitGasControls::new(solve);
    config.initial_dt_s = 0.5;
    config.maximum_dt_s = 0.5;
    config.minimum_dt_s = 1e-5;
    config.maximum_comparisons = 1024;
    config.maximum_gas_iterations = 100_000;
    config.enthalpy_absolute_tolerance_j_kg = 2000.0;
    config.temperature_absolute_tolerance_k = 1.0;
    config.mass_absolute_tolerance_kg = 1e-5;
    config.relative_tolerance = 0.0;
    config
}

#[test]
fn scheduled_second_volatile_activates_on_the_native_clock_and_keeps_each_species_account() {
    let (catalog, bindings) = sources([0.0, 0.0]);
    let gas = mixture_gas(&catalog, 320.0, [0.0, 5000.0]);
    let original = scheduled(vec![
        TimedEvent {
            time_s: 0.0,
            event: Event::Birth(pure_cell(1, 0, 300.0, 2.0)),
        },
        TimedEvent {
            time_s: 1.0,
            event: Event::Birth(pure_cell(2, 1, 280.0, 1.0)),
        },
        TimedEvent {
            time_s: 1.0,
            event: Event::Contact(ThermalLink {
                id: 9,
                a: CellId(1),
                b: CellId(2),
                area_m2: 0.0002,
                conductance_w_k: 0.1,
            }),
        },
    ]);
    let spy = CountPhase {
        law: &LIQUIDS[1],
        queries: Cell::new(0),
    };
    let mut second = bound(2, 1, &bindings, 0.02, 20.0);
    second.law = &spy;
    let surfaces = [second, bound(1, 0, &bindings, 0.02, 20.0)];
    let before = advance_deposition_with_gas_timeline(
        &original,
        &surfaces,
        &gas,
        &catalog,
        0.5,
        adaptive_policy(),
        empty_heat,
        || Ok(()),
    )
    .unwrap();
    assert_eq!(spy.queries.get(), 0);
    near(
        before.gas().species_mass_kg(&LIQUIDS[1].id()).unwrap(),
        gas.species_mass_kg(&LIQUIDS[1].id()).unwrap(),
        1e-14,
    );
    assert_eq!(
        before
            .report()
            .species_transferred_mass_kg
            .get(&LIQUIDS[1].id())
            .copied()
            .unwrap_or(0.0),
        0.0
    );
    let candidate = advance_deposition_with_gas_timeline(
        &original,
        &surfaces,
        &gas,
        &catalog,
        2.0,
        adaptive_policy(),
        |thermal, liquids, _| {
            let ids: Vec<_> = liquids.keys().copied().collect();
            assert_eq!(
                ids,
                if thermal.time_s() <= 1.0 {
                    vec![CellId(1)]
                } else {
                    vec![CellId(1), CellId(2)]
                }
            );
            Ok(coefficients((thermal.time_s(), ids)))
        },
        || Ok(()),
    )
    .unwrap();
    assert!(spy.queries.get() > 0);
    let initial = events_at(&original, 0.0);
    let first = advance_deposition_with_gas_adaptive(
        &initial,
        &[surfaces[1]],
        &gas,
        &catalog,
        1.0,
        adaptive_policy(),
        empty_heat,
        || Ok(()),
    )
    .unwrap();
    let at_birth = events_at(first.thermal(), 1.0);
    let final_part = advance_deposition_with_gas_adaptive(
        &at_birth,
        &surfaces,
        first.gas(),
        &catalog,
        2.0,
        adaptive_policy(),
        empty_heat,
        || Ok(()),
    )
    .unwrap();
    assert_eq!(candidate.thermal().cells(), final_part.thermal().cells());
    assert_eq!(candidate.thermal().links(), final_part.thermal().links());
    assert_eq!(candidate.gas(), final_part.gas());
    assert_eq!(
        candidate.payload(),
        Some(&(2.0, vec![CellId(1), CellId(2)]))
    );
    assert_eq!(candidate.report().events, 3);
    near(
        candidate.report().incoming_mass_kg,
        2.0 * LIQUID_MASS,
        1e-16,
    );
    near(candidate.report().external_heat_j, 5.0, 1e-9);
    for i in 0..2 {
        let species = LIQUIDS[i].id();
        let first_mass = first
            .report()
            .species_transferred_mass_kg
            .get(&species)
            .copied()
            .unwrap_or(0.0);
        let dm = candidate.report().species_transferred_mass_kg[&species];
        near(
            dm,
            first_mass + final_part.report().species_transferred_mass_kg[&species],
            1e-15,
        );
        near(
            candidate.gas().species_mass_kg(&species).unwrap()
                - gas.species_mass_kg(&species).unwrap(),
            dm,
            1e-14,
        );
        near(
            LIQUID_MASS - candidate.thermal().cells()[&CellId(i as u64 + 1)].mass_kg,
            dm,
            1e-14,
        );
        near(
            candidate.report().species_mass_residual_kg[&species],
            0.0,
            1e-14,
        );
    }
    assert!(candidate.report().species_transferred_mass_kg[&LIQUIDS[0].id()] > 0.0);
    assert!(candidate.report().species_transferred_mass_kg[&LIQUIDS[1].id()] < 0.0);
    near(
        candidate.thermal().energy_j() - candidate.report().incoming_energy_j
            + canonical_gas_change(&gas, candidate.gas()),
        5.0,
        1e-7,
    );
    near(candidate.report().mass_residual_kg, 0.0, 1e-14);
    assert_eq!(candidate.gas().point().axes().get("local-lot"), Some(&37.0));
}
