//! Two finite stores, independent surface balances and original time quadrature.
use super::super::contact::{ImplicitContactLaw, ImplicitContactState};
use super::super::link::{GrayRadiativeControls, GrayRadiativeLink};
use super::*;
use std::cell::Cell;

mod surface_sources;

const SIGMA: f64 = 5.670_374_419e-8;

struct Material {
    caloric: Law,
    conductivity_intercept: f64,
    conductivity_slope: f64,
}
impl ThermalLaw for Material {
    fn temperature_k(&self, h: f64) -> Result<f64, DepositionError> {
        self.caloric.temperature_k(h)
    }
    fn maximum_temperature_slope(&self) -> f64 {
        1.0
    }
    fn specific_entropy_change_j_kg_k(
        &self,
        a: f64,
        b: f64,
    ) -> Result<Option<f64>, DepositionError> {
        self.caloric.specific_entropy_change_j_kg_k(a, b)
    }
}
impl ConductivityLaw for Material {
    fn conductivity_w_m_k(&self, h: f64) -> Result<f64, DepositionError> {
        self.temperature_k(h)?;
        Ok(self.conductivity_intercept + self.conductivity_slope * h)
    }
}
impl ImplicitConductionLaw for Material {
    fn caloric_law(&self) -> &dyn ThermalLaw {
        self
    }
    fn conductivity_law(&self) -> &dyn ConductivityLaw {
        self
    }
    fn enthalpy_bounds_j_kg(&self) -> [f64; 2] {
        [100.0, 1000.0]
    }
    fn state_at_specific_enthalpy(
        &self,
        h: f64,
    ) -> Result<ImplicitConductionState, DepositionError> {
        let mut state = self.caloric.state_at_specific_enthalpy(h)?;
        state.conductivity_w_m_k = self.conductivity_w_m_k(h)?;
        state.conductivity_slope_kg_w_m_k_j = self.conductivity_slope;
        Ok(state)
    }
}
fn material(plateau: bool, k0: f64, slope: f64) -> Material {
    Material {
        caloric: Law {
            nonlinear_k: false,
            plateau,
        },
        conductivity_intercept: k0,
        conductivity_slope: slope,
    }
}
fn patch(area: f64) -> ContactPatch {
    ContactPatch::symmetric(ThermalLink {
        id: 7,
        a: CellId(0),
        b: CellId(1),
        area_m2: area,
        // The selected native constitutive law replaces this authored value.
        conductance_w_k: 12345.0,
    })
}
fn geometry(distances: [f64; 2]) -> ConductiveLink {
    ConductiveLink {
        link_id: 7,
        half_distances_m: distances,
        contact_resistance_m2_k_w: 0.0,
    }
}

// Independent oracle: bisect the signed physical flux, use absolute Kelvin and
// direct fourth powers. No production surface fraction, secant or tangent.
fn flux(t: [f64; 2], resistance: [f64; 2], eps: [f64; 2], htc: f64) -> f64 {
    let effective = if eps.contains(&0.0) {
        0.0
    } else {
        1.0 / (1.0 / eps[0] + 1.0 / eps[1] - 1.0)
    };
    let bound = (t[0] - t[1]) / (resistance[0] + resistance[1]);
    let (mut lo, mut hi) = (bound.min(0.0), bound.max(0.0));
    for _ in 0..90 {
        let q = 0.5 * (lo + hi);
        let surfaces = [t[0] - resistance[0] * q, t[1] + resistance[1] * q];
        let gap = htc * (surfaces[0] - surfaces[1])
            + effective * SIGMA * (surfaces[0].powi(4) - surfaces[1].powi(4));
        if q > gap {
            hi = q;
        } else {
            lo = q;
        }
    }
    0.5 * (lo + hi)
}

#[test]
fn two_surface_radiation_retains_each_coating_half_cell_and_reversal() {
    let materials = [material(false, 2.0, 0.0), material(false, 4.0, 0.0)];
    let (eps, htc, area) = ([0.8, 0.6], 12.0, 1.7);
    let surfaces = [600.0_f64, 400.0_f64];
    let bodies = [800.0, 300.0];
    let effective = 1.0 / (1.0 / eps[0] + 1.0 / eps[1] - 1.0);
    let radiated = effective * SIGMA * (surfaces[0].powi(4) - surfaces[1].powi(4));
    let conducted = htc * (surfaces[0] - surfaces[1]);
    let q = radiated + conducted;
    let resistance = [(bodies[0] - surfaces[0]) / q, (surfaces[1] - bodies[1]) / q];
    let coatings = [0.25 * resistance[0], 0.4 * resistance[1]];
    let distances = [
        2.0 * (resistance[0] - coatings[0]),
        4.0 * (resistance[1] - coatings[1]),
    ];
    let exchange = GrayRadiativeLink::new(
        [&materials[0], &materials[1]],
        eps,
        htc,
        coatings,
        GrayRadiativeControls::default(),
    )
    .unwrap();
    let result = exchange
        .evaluate(bodies, geometry(distances), patch(area), &mut || Ok(()))
        .unwrap();
    for side in 0..2 {
        near(
            result.surface_temperatures_k.unwrap()[side],
            surfaces[side],
            2e-10,
        );
    }
    near(result.power_w, area * q, 2e-8);
    near(result.gap_power_w, area * conducted, 2e-8);
    near(result.radiation_power_w, area * radiated, 2e-8);
    near(result.surface_residual_w, 0.0, 2e-8);
    let slopes = [
        htc + 4.0 * effective * SIGMA * surfaces[0].powi(3),
        htc + 4.0 * effective * SIGMA * surfaces[1].powi(3),
    ];
    let denominator = 1.0 + slopes[0] * resistance[0] + slopes[1] * resistance[1];
    near(
        result.link.power_enthalpy_derivatives_w_kg_j[0],
        area * slopes[0] / denominator,
        1e-10,
    );
    near(
        result.link.power_enthalpy_derivatives_w_kg_j[1],
        -area * slopes[1] / denominator,
        1e-10,
    );
    assert!(result.iterations < 16);

    let reversed = GrayRadiativeLink::new(
        [&materials[1], &materials[0]],
        [eps[1], eps[0]],
        htc,
        [coatings[1], coatings[0]],
        GrayRadiativeControls::default(),
    )
    .unwrap()
    .evaluate(
        [bodies[1], bodies[0]],
        geometry([distances[1], distances[0]]),
        patch(area),
        &mut || Ok(()),
    )
    .unwrap();
    near(reversed.power_w, -result.power_w, 2e-8);
    near(
        reversed.link.conductance_w_k,
        result.link.conductance_w_k,
        1e-11,
    );
    for side in 0..2 {
        near(
            reversed.surface_temperatures_k.unwrap()[side],
            surfaces[1 - side],
            2e-10,
        );
        near(
            reversed.link.power_enthalpy_derivatives_w_kg_j[side],
            -result.link.power_enthalpy_derivatives_w_kg_j[1 - side],
            1e-10,
        );
    }
}

#[test]
fn exact_enthalpy_tangents_retain_latent_conductivity_and_near_equilibrium_heat() {
    let materials = [material(true, 0.4, 0.004), material(false, 0.2, 0.001)];
    let exchange = GrayRadiativeLink::new(
        [&materials[0], &materials[1]],
        [0.7, 0.9],
        5.0,
        [0.01, 0.03],
        GrayRadiativeControls::default(),
    )
    .unwrap();
    for h in [[400.0, 650.0], [700.0, 250.0]] {
        let evaluate = |h| {
            exchange
                .evaluate(h, geometry([0.1, 0.04]), patch(0.8), &mut || Ok(()))
                .unwrap()
        };
        let response = evaluate(h);
        for side in 0..2 {
            let mut plus = h;
            let mut minus = h;
            plus[side] += 0.001;
            minus[side] -= 0.001;
            let finite_difference = (evaluate(plus).power_w - evaluate(minus).power_w) / 0.002;
            near(
                response.link.power_enthalpy_derivatives_w_kg_j[side],
                finite_difference,
                2e-7,
            );
        }
        if h[0] == 400.0 {
            assert_eq!(
                materials[0]
                    .state_at_specific_enthalpy(h[0])
                    .unwrap()
                    .temperature_slope_kg_k_j,
                0.0
            );
            assert!(response.link.power_enthalpy_derivatives_w_kg_j[0].abs() > 1e-4);
        }
    }
    let sensible = material(false, 2.0, 0.0);
    let equal = GrayRadiativeLink::new(
        [&sensible; 2],
        [1.0; 2],
        3.0,
        [0.01, 0.03],
        GrayRadiativeControls::default(),
    )
    .unwrap();
    let t = 512.0_f64;
    let surface_g = 3.0 + 4.0 * SIGMA * t.powi(3);
    let expected_g = 0.8 / (0.1 / 2.0 + 0.04 / 2.0 + 0.04 + 1.0 / surface_g);
    for hot in [t, t.next_up()] {
        let response = equal
            .evaluate([hot, t], geometry([0.1, 0.04]), patch(0.8), &mut || Ok(()))
            .unwrap();
        near(response.link.conductance_w_k, expected_g, 2e-13);
        near(response.power_w, expected_g * (hot - t), 2e-26);
        near(
            response.link.power_enthalpy_derivatives_w_kg_j[0],
            expected_g,
            2e-13,
        );
        near(
            response.link.power_enthalpy_derivatives_w_kg_j[1],
            -expected_g,
            2e-13,
        );
    }
}

#[test]
fn gap_only_radiation_only_and_insulating_limits_keep_one_finite_body_transfer() {
    let conductive = material(false, 2.0, 0.0);
    let insulating = material(false, 0.0, 0.0);
    let (bodies, distances, coatings, area) = ([700.0, 300.0], [0.1, 0.2], [0.02, 0.03], 0.8);
    for (eps, htc) in [([0.0, 0.8], 5.0), ([0.8, 0.7], 0.0), ([0.0; 2], 0.0)] {
        let law = GrayRadiativeLink::new(
            [&conductive; 2],
            eps,
            htc,
            coatings,
            GrayRadiativeControls::default(),
        )
        .unwrap();
        let response = law
            .evaluate(bodies, geometry(distances), patch(area), &mut || Ok(()))
            .unwrap();
        let expected = area * flux(bodies, [0.07, 0.13], eps, htc);
        near(response.power_w, expected, 2e-8);
        near(
            response.power_w,
            response.gap_power_w + response.radiation_power_w,
            2e-8,
        );
        if eps[0] == 0.0 {
            assert_eq!(response.radiation_power_w, 0.0);
        }
        if htc == 0.0 {
            assert_eq!(response.gap_power_w, 0.0);
        }
    }
    for laws in [
        [&insulating, &conductive],
        [&conductive, &insulating],
        [&insulating; 2],
    ] {
        let law = GrayRadiativeLink::new(
            laws.map(|law| law as &dyn ImplicitConductionLaw),
            [0.8; 2],
            5.0,
            coatings,
            GrayRadiativeControls::default(),
        )
        .unwrap();
        let response = law
            .evaluate(bodies, geometry(distances), patch(area), &mut || Ok(()))
            .unwrap();
        assert_eq!(response.power_w, 0.0);
        assert_eq!(response.link.power_enthalpy_derivatives_w_kg_j, [0.0; 2]);
        match (
            laws[0].conductivity_intercept,
            laws[1].conductivity_intercept,
        ) {
            (0.0, 0.0) => assert_eq!(response.surface_temperatures_k, None),
            (0.0, _) => assert_eq!(response.surface_temperatures_k, Some([bodies[1]; 2])),
            _ => assert_eq!(response.surface_temperatures_k, Some([bodies[0]; 2])),
        }
        assert_eq!(
            law.evaluate(
                [1100.0, 300.0],
                geometry(distances),
                patch(area),
                &mut || Ok(())
            ),
            Err(DepositionError::MaterialDomain)
        );
    }
}

fn pair<'a>(
    materials: [&'a Material; 2],
    h: [f64; 2],
    masses: [f64; 2],
    face: ConductiveLink,
    area: f64,
) -> (DepositionSimulation<'a>, ConductionObserver<'a>) {
    let mut events = Vec::new();
    for side in 0..2 {
        events.push(TimedEvent {
            time_s: 0.0,
            event: Event::Birth(DepositedCell {
                id: CellId(side as u64),
                material_id: MaterialId(side as u64),
                body_id: side as u64,
                segment_id: 0,
                mass_kg: masses[side],
                specific_h_j_kg: h[side],
                geometry: CellGeometry {
                    start_m: [side as f64, 0.0, 0.0],
                    end_m: [side as f64 + 1.0, 0.0, 0.0],
                    volume_m3: 1.0,
                    surface_area_m2: area,
                },
                // Fully occupied faces must not leak to this authored ambient.
                boundary: Boundary {
                    ambient_k: 200.0,
                    htc_w_m2_k: 999.0,
                    source_w: 0.0,
                    reservoir: None,
                },
            }),
        });
    }
    events.push(TimedEvent {
        time_s: 0.0,
        event: Event::Contact(patch(area).link),
    });
    let thermal = DepositionSimulation::new(
        events,
        &[(MaterialId(0), materials[0]), (MaterialId(1), materials[1])],
        DepositionBudget {
            max_cells: 2,
            max_links: 1,
            max_events: 3,
            max_substeps: 1000,
        },
    )
    .unwrap()
    .advance_to(
        0.0,
        StepConfig {
            max_dt_s: 1.0,
            energy_tolerance_j: 0.0,
        },
        || Ok(()),
    )
    .unwrap()
    .0;
    let observer = ConductionObserver::new(
        &thermal,
        &[(MaterialId(0), materials[0]), (MaterialId(1), materials[1])],
        &[face],
        &[],
        || Ok(()),
    )
    .unwrap();
    (thermal, observer)
}

struct Reference<'a> {
    materials: [&'a Material; 2],
    masses: [f64; 2],
    distances: [f64; 2],
    coatings: [f64; 2],
    eps: [f64; 2],
    htc: f64,
    gap: Option<&'a Gap<'a>>,
    area: f64,
}
impl Reference<'_> {
    fn temperatures(&self, h: [f64; 2]) -> [f64; 2] {
        [
            self.materials[0].temperature_k(h[0]).unwrap(),
            self.materials[1].temperature_k(h[1]).unwrap(),
        ]
    }
    fn power(&self, h: [f64; 2]) -> f64 {
        let resistance = [
            self.distances[0] / self.materials[0].conductivity_w_m_k(h[0]).unwrap()
                + self.coatings[0],
            self.distances[1] / self.materials[1].conductivity_w_m_k(h[1]).unwrap()
                + self.coatings[1],
        ];
        let htc = self.gap.map_or(self.htc, |gap| 1.0 / gap.resistance(h));
        self.area * flux(self.temperatures(h), resistance, self.eps, htc)
    }
    // Reduce the independent endpoint balance to one transferred energy E,
    // preserving the unequal reference masses and both original caloric laws.
    fn stage(&self, base: [f64; 2], dt: f64) -> [f64; 2] {
        let mut lo = (self.masses[0] * (base[0] - 1000.0)).max(self.masses[1] * (100.0 - base[1]));
        let mut hi = (self.masses[0] * (base[0] - 100.0)).min(self.masses[1] * (1000.0 - base[1]));
        for _ in 0..90 {
            let energy = 0.5 * (lo + hi);
            let h = [
                base[0] - energy / self.masses[0],
                base[1] + energy / self.masses[1],
            ];
            if energy > dt * self.power(h) {
                hi = energy;
            } else {
                lo = energy;
            }
        }
        let energy = 0.5 * (lo + hi);
        [
            base[0] - energy / self.masses[0],
            base[1] + energy / self.masses[1],
        ]
    }
    fn generation(&self, h: [f64; 2]) -> f64 {
        let t = self.temperatures(h);
        self.power(h) * (1.0 / t[1] - 1.0 / t[0])
    }
}

#[test]
fn finite_unequal_mass_bodies_use_original_be_storage_and_entropy_on_chain_and_graph() {
    for plateau in [false, true] {
        let materials = [material(false, 0.8, 0.002), material(plateau, 0.5, 0.002)];
        let laws = [&materials[0], &materials[1]];
        let reference = Reference {
            materials: laws,
            masses: [2.0, 3.0],
            distances: [0.04, 0.08],
            coatings: [0.01, 0.02],
            eps: [0.8, 0.6],
            htc: 6.0,
            gap: None,
            area: 0.7,
        };
        let initial = [800.0, if plateau { 400.0 } else { 250.0 }];
        let dt = 0.05;
        let (thermal, observer) = pair(
            laws,
            initial,
            reference.masses,
            geometry(reference.distances),
            reference.area,
        );
        let exchange = GrayRadiativeLink::new(
            laws.map(|law| law as &dyn ImplicitConductionLaw),
            reference.eps,
            reference.htc,
            reference.coatings,
            GrayRadiativeControls::default(),
        )
        .unwrap();
        let materials: [(MaterialId, &dyn ImplicitConductionLaw); 2] =
            [(MaterialId(0), laws[0]), (MaterialId(1), laws[1])];
        let chain = ImplicitConductionChain::new(&thermal, &observer, &materials, || Ok(()))
            .unwrap()
            .with_link_laws(&thermal, &[(7, &exchange)], || Ok(()))
            .unwrap();
        let graph = ImplicitConductionGraph::new(
            &thermal,
            &observer,
            &materials,
            ImplicitConductionLinearControls::default(),
            || Ok(()),
        )
        .unwrap()
        .with_link_laws(&thermal, &[(7, &exchange)], || Ok(()))
        .unwrap();
        let expected = reference.stage(initial, dt);
        let expected_generation = dt * reference.generation(expected);
        let expected_entropy: f64 = (0..2)
            .map(|i| {
                reference.masses[i]
                    * laws[i]
                        .specific_entropy_change_j_kg_k(initial[i], expected[i])
                        .unwrap()
                        .unwrap()
            })
            .sum();
        let mut workspace = ImplicitConductionWorkspace::new(2).unwrap();
        for result in [
            chain.advance_to(&thermal, dt, controls(dt), &mut workspace, || Ok(())),
            graph.advance_to(&thermal, dt, controls(dt), &mut workspace, || Ok(())),
        ] {
            let (next, report) = result.unwrap();
            for side in 0..2 {
                near(
                    next.cells()[&CellId(side as u64)].specific_h_j_kg,
                    expected[side],
                    2e-9,
                );
            }
            near(report.energy.external_heat_j, 0.0, 0.0);
            near(report.energy.residual_j, 0.0, 1e-9);
            near(report.entropy.reservoir_entropy_input_j_k, 0.0, 0.0);
            near(
                report.entropy.conduction_generation_j_k,
                expected_generation,
                2e-10,
            );
            near(
                report.entropy.thermal_storage_change_j_k,
                expected_entropy,
                2e-10,
            );
            assert!(report.entropy.conduction_generation_j_k > 0.0);
            assert!(!report.entropy.has_unclassified_sources);
            assert_eq!(report.energy.substeps, 1);
            assert!(next.cells()[&CellId(0)].specific_h_j_kg < initial[0]);
            assert!(next.cells()[&CellId(1)].specific_h_j_kg > initial[1]);
            if plateau {
                assert_eq!(next.temperature_k(CellId(1)).unwrap(), 300.0);
            }
        }
        assert_eq!(thermal.cells()[&CellId(0)].specific_h_j_kg, initial[0]);
        assert_eq!(thermal.cells()[&CellId(1)].specific_h_j_kg, initial[1]);
    }
}

#[test]
fn finite_bodies_use_original_sdirk_stages_for_internal_heat_and_entropy() {
    for plateau in [false, true] {
        let materials = [material(false, 0.8, 0.002), material(plateau, 0.5, 0.002)];
        let laws = [&materials[0], &materials[1]];
        let reference = Reference {
            materials: laws,
            masses: [2.0, 3.0],
            distances: [0.04, 0.08],
            coatings: [0.01, 0.02],
            eps: [0.8, 0.6],
            htc: 6.0,
            gap: None,
            area: 0.7,
        };
        let initial = [800.0, if plateau { 400.0 } else { 250.0 }];
        let dt = 0.05;
        let (thermal, observer) = pair(
            laws,
            initial,
            reference.masses,
            geometry(reference.distances),
            reference.area,
        );
        let exchange = GrayRadiativeLink::new(
            laws.map(|law| law as &dyn ImplicitConductionLaw),
            reference.eps,
            reference.htc,
            reference.coatings,
            GrayRadiativeControls::default(),
        )
        .unwrap();
        let materials: [(MaterialId, &dyn ImplicitConductionLaw); 2] =
            [(MaterialId(0), laws[0]), (MaterialId(1), laws[1])];
        let chain = ImplicitConductionChain::new(&thermal, &observer, &materials, || Ok(()))
            .unwrap()
            .with_link_laws(&thermal, &[(7, &exchange)], || Ok(()))
            .unwrap();
        let graph = ImplicitConductionGraph::new(
            &thermal,
            &observer,
            &materials,
            ImplicitConductionLinearControls::default(),
            || Ok(()),
        )
        .unwrap()
        .with_link_laws(&thermal, &[(7, &exchange)], || Ok(()))
        .unwrap();
        let gamma = 1.0 - 1.0 / 2.0_f64.sqrt();
        let first = reference.stage(initial, gamma * dt);
        let earlier_energy = (1.0 - gamma) * dt * reference.power(first);
        let base = [
            initial[0] - earlier_energy / reference.masses[0],
            initial[1] + earlier_energy / reference.masses[1],
        ];
        let expected = reference.stage(base, gamma * dt);
        let expected_generation = dt
            * ((1.0 - gamma) * reference.generation(first)
                + gamma * reference.generation(expected));
        let mut workspace = ImplicitConductionWorkspace::new(2).unwrap();
        for result in [
            chain.advance_sdirk2_to(&thermal, dt, controls(dt), &mut workspace, || Ok(())),
            graph.advance_sdirk2_to(&thermal, dt, controls(dt), &mut workspace, || Ok(())),
        ] {
            let (next, report) = result.unwrap();
            for side in 0..2 {
                near(
                    next.cells()[&CellId(side as u64)].specific_h_j_kg,
                    expected[side],
                    2e-9,
                );
            }
            near(report.energy.external_heat_j, 0.0, 0.0);
            near(report.energy.residual_j, 0.0, 1e-9);
            near(report.entropy.reservoir_entropy_input_j_k, 0.0, 0.0);
            near(
                report.entropy.conduction_generation_j_k,
                expected_generation,
                2e-10,
            );
            assert!(report.entropy.conduction_generation_j_k > 0.0);
            assert!(!report.entropy.has_unclassified_sources);
            assert_eq!(report.energy.substeps, 1);
        }
    }
}

#[test]
fn gray_links_refuse_ambiguous_faces_source_extrapolation_exhaustion_and_cancellation() {
    let law = material(false, 2.0, 0.0);
    let mut control = GrayRadiativeControls::default();
    assert!(GrayRadiativeLink::new([&law; 2], [-0.1, 0.5], 0.0, [0.0; 2], control).is_err());
    assert!(GrayRadiativeLink::new([&law; 2], [0.5; 2], -1.0, [0.0; 2], control).is_err());
    assert!(GrayRadiativeLink::new([&law; 2], [0.5; 2], 1.0, [-1.0, 0.0], control).is_err());
    let exchange = GrayRadiativeLink::new([&law; 2], [1.0; 2], 0.0, [0.0; 2], control).unwrap();
    let face = geometry([0.1, 0.2]);
    let mut ambiguous = face;
    ambiguous.contact_resistance_m2_k_w = 0.1;
    assert!(
        exchange
            .evaluate([800.0, 300.0], ambiguous, patch(1.0), &mut || Ok(()))
            .is_err()
    );
    let mut unequal = patch(1.0);
    unequal.covered_area_m2[1] = 2.0;
    assert!(
        exchange
            .evaluate([800.0, 300.0], face, unequal, &mut || Ok(()))
            .is_err()
    );
    assert_eq!(
        exchange.evaluate([99.0, 300.0], face, patch(1.0), &mut || Ok(())),
        Err(DepositionError::MaterialDomain)
    );
    let mut calls = 0;
    assert_eq!(
        exchange.evaluate([800.0, 300.0], face, patch(1.0), &mut || {
            calls += 1;
            if calls >= 5 {
                Err(DepositionError::Cancelled)
            } else {
                Ok(())
            }
        }),
        Err(DepositionError::Cancelled)
    );
    assert_eq!(calls, 5);
    control.maximum_iterations = 1;
    control.absolute_power_tolerance_w = 0.0;
    control.relative_power_tolerance = 0.0;
    let exhausted = GrayRadiativeLink::new([&law; 2], [1.0; 2], 0.0, [0.0; 2], control).unwrap();
    assert!(
        exhausted
            .evaluate([800.0, 300.0], face, patch(1.0), &mut || Ok(()))
            .is_err()
    );

    let other = material(false, 2.0, 0.0);
    let (thermal, observer) = pair([&law; 2], [800.0, 300.0], [2.0, 3.0], face, 1.0);
    let wrong = GrayRadiativeLink::new(
        [&law, &other],
        [0.8; 2],
        1.0,
        [0.0; 2],
        GrayRadiativeControls::default(),
    )
    .unwrap();
    let chain = ImplicitConductionChain::new(
        &thermal,
        &observer,
        &[(MaterialId(0), &law), (MaterialId(1), &law)],
        || Ok(()),
    )
    .unwrap();
    assert!(
        chain
            .with_link_laws(&thermal, &[(7, &wrong)], || Ok(()))
            .is_err()
    );
    assert_eq!(thermal.cells()[&CellId(0)].specific_h_j_kg, 800.0);
    assert_eq!(thermal.cells()[&CellId(1)].specific_h_j_kg, 300.0);
}

// Explicit original bulk-enthalpy resistance, including a phase-sensitive
// contribution when dT/dh is zero. Counters observe source-query ownership.
struct Gap<'a> {
    materials: [&'a Material; 2],
    bounds: [[f64; 2]; 2],
    offset: f64,
    slopes: [f64; 2],
    queries: Cell<usize>,
    last_h: Cell<[f64; 2]>,
}
impl<'a> Gap<'a> {
    fn new(materials: [&'a Material; 2]) -> Self {
        Self {
            materials,
            bounds: [[100.0, 1000.0]; 2],
            offset: 0.01,
            slopes: [1e-4, -1e-5],
            queries: Cell::new(0),
            last_h: Cell::new([0.0; 2]),
        }
    }
    fn resistance(&self, h: [f64; 2]) -> f64 {
        self.offset + self.slopes[0] * h[0] + self.slopes[1] * h[1]
    }
}
impl ImplicitContactLaw for Gap<'_> {
    fn caloric_laws(&self) -> [&dyn ThermalLaw; 2] {
        [self.materials[0], self.materials[1]]
    }
    fn enthalpy_bounds_j_kg(&self) -> [[f64; 2]; 2] {
        self.bounds
    }
    fn state_at_specific_enthalpies(
        &self,
        h: [f64; 2],
    ) -> Result<ImplicitContactState, DepositionError> {
        self.queries.set(self.queries.get() + 1);
        self.last_h.set(h);
        if h.into_iter()
            .zip(self.bounds)
            .any(|(h, [lo, hi])| !(lo..=hi).contains(&h))
        {
            return Err(DepositionError::MaterialDomain);
        }
        Ok(ImplicitContactState {
            resistance_m2_k_w: self.resistance(h),
            resistance_enthalpy_derivatives_m2_kg_k_w_j: self.slopes,
        })
    }
}

#[test]
fn source_gap_shares_radiating_surfaces_and_preserves_sensible_and_latent_tangents() {
    let materials = [material(true, 2.0, 0.0), material(false, 3.0, 0.0)];
    let laws = [&materials[0], &materials[1]];
    let gap = Gap::new(laws);
    let (distances, coatings, area, eps) = ([0.1, 0.07], [0.01, 0.02], 0.8, [0.7, 0.9]);
    let resistance = [
        distances[0] / 2.0 + coatings[0],
        distances[1] / 3.0 + coatings[1],
    ];
    let exchange = GrayRadiativeLink::new(
        laws.map(|m| m as &dyn ImplicitConductionLaw),
        eps,
        0.0,
        coatings,
        GrayRadiativeControls::default(),
    )
    .unwrap()
    .with_gap_resistance(&gap)
    .unwrap();
    let oracle = |h: [f64; 2]| {
        area * flux(
            [
                laws[0].temperature_k(h[0]).unwrap(),
                laws[1].temperature_k(h[1]).unwrap(),
            ],
            resistance,
            eps,
            1.0 / gap.resistance(h),
        )
    };
    for h in [[400.0, 650.0], [700.0, 250.0]] {
        let before = gap.queries.get();
        let response = exchange
            .evaluate(h, geometry(distances), patch(area), &mut || Ok(()))
            .unwrap();
        assert_eq!(gap.queries.get(), before + 1);
        assert_eq!(gap.last_h.get(), h);
        assert!(response.iterations > 1);
        let q = oracle(h) / area;
        let t = [
            laws[0].temperature_k(h[0]).unwrap(),
            laws[1].temperature_k(h[1]).unwrap(),
        ];
        let surfaces = [t[0] - resistance[0] * q, t[1] + resistance[1] * q];
        near(response.power_w, area * q, 2e-8);
        for side in 0..2 {
            near(
                response.surface_temperatures_k.unwrap()[side],
                surfaces[side],
                2e-10,
            );
            let (mut plus, mut minus) = (h, h);
            plus[side] += 0.001;
            minus[side] -= 0.001;
            near(
                response.link.power_enthalpy_derivatives_w_kg_j[side],
                (oracle(plus) - oracle(minus)) / 0.002,
                2e-7,
            );
        }
        near(
            response.gap_power_w,
            area * (surfaces[0] - surfaces[1]) / gap.resistance(h),
            2e-8,
        );
        near(
            response.power_w,
            response.gap_power_w + response.radiation_power_w,
            2e-8,
        );
        if h[0] == 400.0 {
            assert_eq!(
                laws[0]
                    .state_at_specific_enthalpy(h[0])
                    .unwrap()
                    .temperature_slope_kg_k_j,
                0.0
            );
            // k is constant here: this derivative can ONLY come from R'(h).
            assert!(response.link.power_enthalpy_derivatives_w_kg_j[0].abs() > 1e-4);
        }
    }
}

#[test]
fn variable_gap_updates_original_be_and_sdirk_heat_and_entropy_on_chain_and_graph() {
    for plateau in [false, true] {
        let materials = [material(false, 0.8, 0.002), material(plateau, 0.5, 0.002)];
        let laws = [&materials[0], &materials[1]];
        let gap = Gap::new(laws);
        let reference = Reference {
            materials: laws,
            masses: [2.0, 3.0],
            distances: [0.04, 0.08],
            coatings: [0.01, 0.02],
            eps: [0.8, 0.6],
            htc: 0.0,
            gap: Some(&gap),
            area: 0.7,
        };
        let initial = [800.0, if plateau { 400.0 } else { 250.0 }];
        let dt = 0.05;
        let (thermal, observer) = pair(
            laws,
            initial,
            reference.masses,
            geometry(reference.distances),
            reference.area,
        );
        let exchange = GrayRadiativeLink::new(
            laws.map(|m| m as &dyn ImplicitConductionLaw),
            reference.eps,
            0.0,
            reference.coatings,
            GrayRadiativeControls::default(),
        )
        .unwrap()
        .with_gap_resistance(&gap)
        .unwrap();
        let bindings: [(MaterialId, &dyn ImplicitConductionLaw); 2] =
            [(MaterialId(0), laws[0]), (MaterialId(1), laws[1])];
        let chain = ImplicitConductionChain::new(&thermal, &observer, &bindings, || Ok(()))
            .unwrap()
            .with_link_laws(&thermal, &[(7, &exchange)], || Ok(()))
            .unwrap();
        let graph = ImplicitConductionGraph::new(
            &thermal,
            &observer,
            &bindings,
            ImplicitConductionLinearControls::default(),
            || Ok(()),
        )
        .unwrap()
        .with_link_laws(&thermal, &[(7, &exchange)], || Ok(()))
        .unwrap();
        let mut workspace = ImplicitConductionWorkspace::new(2).unwrap();
        for sdirk in [false, true] {
            let (expected, generation) = if sdirk {
                let gamma = 1.0 - 1.0 / 2.0_f64.sqrt();
                let first = reference.stage(initial, gamma * dt);
                let earlier = (1.0 - gamma) * dt * reference.power(first);
                let base = [
                    initial[0] - earlier / reference.masses[0],
                    initial[1] + earlier / reference.masses[1],
                ];
                let expected = reference.stage(base, gamma * dt);
                (
                    expected,
                    dt * ((1.0 - gamma) * reference.generation(first)
                        + gamma * reference.generation(expected)),
                )
            } else {
                let expected = reference.stage(initial, dt);
                (expected, dt * reference.generation(expected))
            };
            let results = if sdirk {
                [
                    chain.advance_sdirk2_to(&thermal, dt, controls(dt), &mut workspace, || Ok(())),
                    graph.advance_sdirk2_to(&thermal, dt, controls(dt), &mut workspace, || Ok(())),
                ]
            } else {
                [
                    chain.advance_to(&thermal, dt, controls(dt), &mut workspace, || Ok(())),
                    graph.advance_to(&thermal, dt, controls(dt), &mut workspace, || Ok(())),
                ]
            };
            for result in results {
                let (next, report) = result.unwrap();
                for side in 0..2 {
                    near(
                        next.cells()[&CellId(side as u64)].specific_h_j_kg,
                        expected[side],
                        2e-9,
                    );
                }
                near(report.energy.external_heat_j, 0.0, 0.0);
                near(report.energy.residual_j, 0.0, 1e-9);
                near(report.entropy.reservoir_entropy_input_j_k, 0.0, 0.0);
                near(report.entropy.conduction_generation_j_k, generation, 2e-10);
                assert!(report.entropy.conduction_generation_j_k > 0.0);
                assert!(!report.entropy.has_unclassified_sources);
                assert_eq!(report.energy.substeps, 1);
                if plateau {
                    assert_eq!(next.temperature_k(CellId(1)).unwrap(), 300.0);
                }
            }
        }
        assert_eq!(thermal.cells()[&CellId(0)].specific_h_j_kg, initial[0]);
        assert_eq!(thermal.cells()[&CellId(1)].specific_h_j_kg, initial[1]);
    }
}

#[test]
fn zero_emissivity_recovers_original_source_contact_without_an_extra_bulk_resistance() {
    let materials = [material(false, 0.8, 0.002), material(true, 0.5, 0.002)];
    let laws = [&materials[0], &materials[1]];
    let gap = Gap::new(laws);
    let (initial, masses, distances, area, dt) =
        ([800.0, 400.0], [2.0, 3.0], [0.04, 0.08], 0.7, 0.05);
    let exchange = GrayRadiativeLink::new(
        laws.map(|m| m as &dyn ImplicitConductionLaw),
        [0.0; 2],
        0.0,
        [0.0; 2],
        GrayRadiativeControls::default(),
    )
    .unwrap()
    .with_gap_resistance(&gap)
    .unwrap();
    let response = exchange
        .evaluate(initial, geometry(distances), patch(area), &mut || Ok(()))
        .unwrap();
    let resistance = distances[0] / laws[0].conductivity_w_m_k(initial[0]).unwrap()
        + gap.resistance(initial)
        + distances[1] / laws[1].conductivity_w_m_k(initial[1]).unwrap();
    near(response.link.conductance_w_k, area / resistance, 1e-12);
    assert_eq!(response.radiation_power_w, 0.0);
    let (thermal, observer) = pair(laws, initial, masses, geometry(distances), area);
    let bindings: [(MaterialId, &dyn ImplicitConductionLaw); 2] =
        [(MaterialId(0), laws[0]), (MaterialId(1), laws[1])];
    let legacy = ImplicitConductionChain::new(&thermal, &observer, &bindings, || Ok(()))
        .unwrap()
        .with_contacts(&thermal, &[(7, &gap)], || Ok(()))
        .unwrap();
    let composed = ImplicitConductionChain::new(&thermal, &observer, &bindings, || Ok(()))
        .unwrap()
        .with_link_laws(&thermal, &[(7, &exchange)], || Ok(()))
        .unwrap();
    let mut workspace = ImplicitConductionWorkspace::new(2).unwrap();
    let (expected, old_report) = legacy
        .advance_to(&thermal, dt, controls(dt), &mut workspace, || Ok(()))
        .unwrap();
    let (next, report) = composed
        .advance_to(&thermal, dt, controls(dt), &mut workspace, || Ok(()))
        .unwrap();
    for side in 0..2 {
        near(
            next.cells()[&CellId(side as u64)].specific_h_j_kg,
            expected.cells()[&CellId(side as u64)].specific_h_j_kg,
            2e-9,
        );
    }
    near(report.energy.residual_j, old_report.energy.residual_j, 1e-9);
    near(
        report.entropy.conduction_generation_j_k,
        old_report.entropy.conduction_generation_j_k,
        2e-10,
    );
    assert!(
        legacy
            .with_link_laws(&thermal, &[(7, &exchange)], || Ok(()))
            .is_err()
    );
    assert!(
        composed
            .with_contacts(&thermal, &[(7, &gap)], || Ok(()))
            .is_err()
    );
}

#[test]
fn gap_source_identity_bounds_finite_resistance_and_cancellation_refuse_without_publication() {
    let materials = [material(false, 2.0, 0.0), material(false, 3.0, 0.0)];
    let laws = [&materials[0], &materials[1]];
    let gap = Gap::new(laws);
    let make = |htc| {
        GrayRadiativeLink::new(
            laws.map(|m| m as &dyn ImplicitConductionLaw),
            [0.8; 2],
            htc,
            [0.0; 2],
            GrayRadiativeControls::default(),
        )
        .unwrap()
    };
    assert!(make(1.0).with_gap_resistance(&gap).is_err());
    assert!(
        make(0.0)
            .with_gap_resistance(&gap)
            .unwrap()
            .with_gap_resistance(&gap)
            .is_err()
    );
    let alien = material(false, 2.0, 0.0);
    let wrong = Gap::new([&alien, laws[1]]);
    assert!(make(0.0).with_gap_resistance(&wrong).is_err());
    let reversed = Gap::new([laws[1], laws[0]]);
    assert!(make(0.0).with_gap_resistance(&reversed).is_err());
    let mut narrow = Gap::new(laws);
    narrow.bounds[0] = [500.0, 900.0];
    let bounded = make(0.0).with_gap_resistance(&narrow).unwrap();
    assert_eq!(
        bounded.evaluate([400.0, 300.0], geometry([0.1; 2]), patch(1.0), &mut || Ok(
            ()
        )),
        Err(DepositionError::MaterialDomain)
    );
    assert_eq!(narrow.queries.get(), 0);
    for resistance in [0.0, -0.1, f64::NAN] {
        let invalid = Gap {
            offset: resistance,
            slopes: [0.0; 2],
            ..Gap::new(laws)
        };
        let law = make(0.0).with_gap_resistance(&invalid).unwrap();
        assert_eq!(
            law.evaluate([800.0, 300.0], geometry([0.1; 2]), patch(1.0), &mut || Ok(
                ()
            )),
            Err(DepositionError::MaterialDomain)
        );
    }
    let exchange = make(0.0).with_gap_resistance(&gap).unwrap();
    let mut calls = 0;
    assert_eq!(
        exchange.evaluate([800.0, 300.0], geometry([0.1; 2]), patch(1.0), &mut || {
            calls += 1;
            if calls == 5 {
                Err(DepositionError::Cancelled)
            } else {
                Ok(())
            }
        }),
        Err(DepositionError::Cancelled)
    );
    assert_eq!(gap.queries.get(), 1);
    assert_eq!(calls, 5);

    let (thermal, observer) = pair(laws, [800.0, 300.0], [2.0, 3.0], geometry([0.1; 2]), 1.0);
    let bindings: [(MaterialId, &dyn ImplicitConductionLaw); 2] =
        [(MaterialId(0), laws[0]), (MaterialId(1), laws[1])];
    let chain = ImplicitConductionChain::new(&thermal, &observer, &bindings, || Ok(()))
        .unwrap()
        .with_link_laws(&thermal, &[(7, &exchange)], || Ok(()))
        .unwrap();
    let mut workspace = ImplicitConductionWorkspace::new(2).unwrap();
    assert!(
        chain
            .advance_to(&thermal, 0.05, controls(0.05), &mut workspace, || Err(
                DepositionError::Cancelled
            ))
            .is_err()
    );
    assert_eq!(thermal.cells()[&CellId(0)].specific_h_j_kg, 800.0);
    assert_eq!(thermal.cells()[&CellId(1)].specific_h_j_kg, 300.0);
}
