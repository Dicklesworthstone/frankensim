//! Coupled reservoir intervals, independently resolved throughput and refusals.
use super::super::enthalpy::adaptive::{
    AdaptiveCoupledBoundarySpeciesEnthalpyDiffusionControls,
    AdaptiveSpeciesEnthalpyDiffusionControls, SpeciesEnthalpyDiffusionModel,
};
use super::*;
use std::cell::Cell;

#[derive(Clone, Copy)]
enum CaloricLaw {
    Affine([f64; 2]),
    Cubic,
    CompetingFaces,
}

struct ReservoirModel {
    law: CaloricLaw,
    diffusivity: f64,
    preparations: Cell<usize>,
    boundary_queries: Cell<usize>,
    fail_at_boundary: usize,
    cancel_at_boundary: usize,
}

#[derive(Clone, Debug, PartialEq)]
struct ReservoirState {
    temperature: Vec<f64>,
    diffusivity: Vec<f64>,
    fraction: Vec<f64>,
}

impl ReservoirModel {
    fn new(law: CaloricLaw, diffusivity: f64) -> Self {
        Self {
            law,
            diffusivity,
            preparations: Cell::new(0),
            boundary_queries: Cell::new(0),
            fail_at_boundary: usize::MAX,
            cancel_at_boundary: usize::MAX,
        }
    }
}

impl SpeciesEnthalpyDiffusionModel<2> for ReservoirModel {
    type Prepared = ReservoirState;
    fn prepare(
        &self,
        domain: &VoxelDomain,
        carrier: &MaterialAdvectionStep,
        species: &[[f64; 2]],
        h: &[f64],
        gate: &CancelGate,
    ) -> Result<Self::Prepared, ChtError> {
        poll(gate)?;
        self.preparations.set(self.preparations.get() + 1);
        let mut temperature = Vec::new();
        let mut diffusivity = Vec::new();
        let mut fraction = Vec::new();
        for c in 0..domain.cell_count() {
            poll_every(gate, c)?;
            if !domain.is_fluid(c) {
                temperature.push(0.0);
                diffusivity.push(0.0);
                fraction.push(0.0);
                continue;
            }
            let y = species[c][0] / carrier.mass_kg[c];
            let t = match self.law {
                CaloricLaw::Affine(partials) => {
                    // Binary caloric chart h=h_B,0+Y_A*(h_A,0-h_B,0)
                    // +(T-300). Anchor the common datum to the fixed carrier,
                    // so roundoff in a mass partition cannot amplify it into
                    // a fictitious temperature. Retain chemical low parts.
                    let (difference, difference_tail) =
                        fs_math::eft::two_sum(partials[0], -partials[1]);
                    let (product, product_tail) = fs_math::eft::two_prod(y, difference);
                    let (extra, extra_tail) = fs_math::eft::two_prod(y, difference_tail);
                    let (sum, tail) = fs_math::eft::two_sum(partials[1], product);
                    300.0 + (h[c] - (sum + (tail + product_tail + extra + extra_tail)))
                }
                CaloricLaw::Cubic => 300.0 + h[c] - 4.0 * y.powi(3),
                // Positive Cp=1+Y_A, h=(1+Y_A)*(T-300)+10Y_A.
                CaloricLaw::CompetingFaces => 300.0 + (h[c] - 10.0 * y) / (1.0 + y),
            };
            if !t.is_finite() || t <= 0.0 {
                return Err(refuse("original reservoir caloric law is out of domain"));
            }
            temperature.push(t);
            diffusivity.push(self.diffusivity);
            fraction.push(y);
        }
        Ok(ReservoirState {
            temperature,
            diffusivity,
            fraction,
        })
    }
    fn temperature_k<'p>(&self, p: &'p Self::Prepared) -> &'p [f64] {
        &p.temperature
    }
    fn diffusivity_m2_s<'p>(&self, p: &'p Self::Prepared) -> &'p [f64] {
        &p.diffusivity
    }
    fn face_partial_enthalpy_j_kg(
        &self,
        _: &Self::Prepared,
        _: &VoxelDomain,
        _: [usize; 2],
        _: usize,
        gate: &CancelGate,
    ) -> Result<[f64; 2], ChtError> {
        poll(gate)?;
        match self.law {
            CaloricLaw::Affine(partials) => Ok(partials),
            _ => Err(refuse(
                "this authored nonlinear fixture has no interior face",
            )),
        }
    }
    fn boundary_partial_enthalpy_j_kg(
        &self,
        p: &Self::Prepared,
        _: &VoxelDomain,
        index: usize,
        boundary: &SpeciesDiffusionBoundary<2>,
        gate: &CancelGate,
    ) -> Result<[f64; 2], ChtError> {
        let calls = self.boundary_queries.get() + 1;
        self.boundary_queries.set(calls);
        if calls == self.cancel_at_boundary {
            gate.request();
        }
        poll(gate)?;
        if calls == self.fail_at_boundary {
            return Err(refuse("authored adaptive boundary source refused"));
        }
        Ok(match self.law {
            CaloricLaw::Affine(partials) => partials,
            CaloricLaw::Cubic => {
                // Original partial derivatives of h=T-300+4Y_A^3.
                let y = p.fraction[boundary.cell];
                let t = p.temperature[boundary.cell] - 300.0;
                [
                    t + 4.0 * (3.0 * y * y - 2.0 * y.powi(3)),
                    t - 8.0 * y.powi(3),
                ]
            }
            CaloricLaw::CompetingFaces => {
                assert_eq!((boundary.cell, boundary.face), (0, index));
                // Same source species Cp=[2,1], formations=[10,0];
                // independently authored physical face temperatures300/310K.
                if index == 0 {
                    [10.0, 0.0]
                } else {
                    [30.0, 10.0]
                }
            }
        })
    }
    fn maximum_prepared_bytes(&self, cells: usize) -> Option<usize> {
        cells
            .checked_mul(24)?
            .checked_add(core::mem::size_of::<ReservoirState>())
    }
    fn caloric_defect_j(&self, _: &Self::Prepared) -> Result<f64, ChtError> {
        Ok(0.0)
    }
}

fn reservoir(cell: usize, face: usize, fraction: f64) -> SpeciesDiffusionBoundary<2> {
    SpeciesDiffusionBoundary {
        cell,
        face,
        reservoir_mass_fractions: [fraction, 1.0 - fraction],
        film_resistance_m2_s_kg: 0.0,
    }
}

fn open_controls() -> AdaptiveCoupledBoundarySpeciesEnthalpyDiffusionControls<2> {
    AdaptiveCoupledBoundarySpeciesEnthalpyDiffusionControls {
        integrator: AdaptiveSpeciesEnthalpyDiffusionControls {
            solver: controls(1),
            initial_substeps: 1,
            maximum_refinements: 12,
            maximum_source_evaluations: 100_000,
            absolute_temperature_tolerance_k: 1e-4,
            absolute_specific_enthalpy_tolerance_j_kg: 1e-4,
            absolute_mass_fraction_tolerance: [1e-8; 2],
            relative_mass_fraction_tolerance: 0.0,
            maximum_enthalpy_defect_j: 1e-7,
        },
        absolute_boundary_transfer_tolerance_kg: [1e-5; 2],
        absolute_boundary_enthalpy_tolerance_j: 1e-4,
    }
}

fn reason(error: ChtError) -> String {
    match error {
        ChtError::InvalidInput { reason, .. } => reason,
        other => panic!("{other:?}"),
    }
}

#[test]
fn adaptive_reservoir_matches_full_interval_film_oracle_and_empty_boundary_parity() {
    let domain = VoxelDomain::new(1, 1, 1, 1.0).unwrap();
    let state = carrier(&domain, &[2.0]);
    let mut boundary = reservoir(0, 1, 0.25);
    boundary.film_resistance_m2_s_kg = 0.5;
    let model = ReservoirModel::new(CaloricLaw::Affine([8.0, 0.0]), 0.2);
    let gate = CancelGate::new_clock_free();
    let got = state
        .diffuse_species_with_enthalpy_and_boundaries_coupled_adaptive(
            &domain,
            &[[1.6, 0.4]],
            &[6.4],
            &[boundary],
            1.3,
            &model,
            open_controls(),
            &gate,
        )
        .unwrap();
    let conductance = 1.0 / (0.5 / (2.0 * 0.2) + 0.5);
    let dm = 2.0 * (0.25 - 0.8) * -(-conductance * 1.3 / 2.0_f64).exp_m1();
    near(got.step.species_mass_kg[0][0], 1.6 + dm, 2e-14);
    near(got.step.specific_enthalpy_j_kg[0], 6.4 + 4.0 * dm, 2e-13);
    near(got.step.temperature_k[0], 300.0, 2e-12);
    near(got.boundary.mass_into_domain_kg[0], dm, 2e-14);
    near(got.thermal_boundary.enthalpy_into_domain_j, 8.0 * dm, 2e-13);
    near(got.boundary.faces[0].conductance_kg_s, conductance, 1e-14);
    assert_eq!(got.step.report.accepted_substeps, 2);
    assert_eq!(got.step.state.temperature, got.step.temperature_k);

    let domain = VoxelDomain::new(2, 1, 1, 1.0).unwrap();
    let state = carrier(&domain, &[1.0; 2]);
    let species = [[0.8, 0.2], [0.2, 0.8]];
    let h = [6.4, 1.6];
    let closed = state
        .diffuse_species_with_enthalpy_adaptive(
            &domain,
            &species,
            &h,
            1.0,
            &model,
            open_controls().integrator,
            &gate,
        )
        .unwrap();
    let open = state
        .diffuse_species_with_enthalpy_and_boundaries_coupled_adaptive(
            &domain,
            &species,
            &h,
            &[],
            1.0,
            &model,
            open_controls(),
            &gate,
        )
        .unwrap();
    assert_eq!(open.step, closed);
    assert!(open.boundary.faces.is_empty());
    assert!(open.thermal_boundary.faces.is_empty());
    assert_eq!(open.estimated_boundary_transfer_error_kg, [0.0; 2]);
    assert_eq!(open.estimated_boundary_enthalpy_error_j, 0.0);
}

#[test]
fn adaptive_reservoir_throughput_refines_independently_of_nearly_stationary_fields() {
    let domain = VoxelDomain::new(1, 1, 1, 1.0).unwrap();
    let state = carrier(&domain, &[1.0]);
    let boundaries = [reservoir(0, 0, 1.0), reservoir(0, 1, 0.0)];
    let model = ReservoirModel::new(CaloricLaw::Affine([8.0, 0.0]), 0.25);
    let gate = CancelGate::new_clock_free();
    let mut relaxed = open_controls();
    relaxed.integrator.absolute_temperature_tolerance_k = 1.0;
    relaxed.integrator.absolute_specific_enthalpy_tolerance_j_kg = 8.0;
    relaxed.integrator.absolute_mass_fraction_tolerance = [1.0; 2];
    relaxed.absolute_boundary_transfer_tolerance_kg = [1e6; 2];
    relaxed.absolute_boundary_enthalpy_tolerance_j = 1e6;
    let run = |control| {
        state
            .diffuse_species_with_enthalpy_and_boundaries_coupled_adaptive(
                &domain,
                &[[0.5, 0.5]],
                &[4.0],
                &boundaries,
                4.0,
                &model,
                control,
                &gate,
            )
            .unwrap()
    };
    let coarse = run(relaxed);
    let mut mass_control = relaxed;
    mass_control.absolute_boundary_transfer_tolerance_kg = [2e-5; 2];
    let mut heat_control = relaxed;
    heat_control.absolute_boundary_enthalpy_tolerance_j = 1.6e-4;
    for fine in [run(mass_control), run(heat_control)] {
        assert!(fine.step.report.accepted_substeps > coarse.step.report.accepted_substeps);
        // Independent continuous law: dY/dt=G(1-2Y), so initial Y=.5
        // remains EXACTLY stationary. G=.5kg/s and duration4s imply 1kg A
        // enters face0 and exits face1, together with 8J in each direction.
        near(fine.step.species_mass_kg[0][0], 0.5, 2e-5);
        near(fine.step.temperature_k[0], 300.0, 1e-10);
        near(fine.boundary.faces[0].mass_into_cell_kg[0], 1.0, 2e-5);
        near(fine.boundary.faces[1].mass_into_cell_kg[0], -1.0, 2e-5);
        near(fine.boundary.faces[0].incoming_mass_kg[0], 1.0, 2e-5);
        near(fine.boundary.faces[0].outgoing_mass_kg[1], 1.0, 2e-5);
        near(fine.boundary.faces[1].outgoing_mass_kg[0], 1.0, 2e-5);
        near(fine.boundary.faces[1].incoming_mass_kg[1], 1.0, 2e-5);
        near(
            fine.thermal_boundary.positive_enthalpy_into_domain_j,
            8.0,
            1.6e-4,
        );
        near(
            fine.thermal_boundary.negative_enthalpy_into_domain_j,
            8.0,
            1.6e-4,
        );
        near(fine.thermal_boundary.enthalpy_into_domain_j, 0.0, 1.6e-4);
        assert!(
            fine.estimated_boundary_transfer_error_kg
                .iter()
                .all(|&v| v <= 2e-5)
        );
        assert!(fine.estimated_boundary_enthalpy_error_j <= 1.6e-4);
    }
}

#[test]
fn adaptive_reservoir_external_energy_tails_survive_zero_joule_admission() {
    let domain = VoxelDomain::new(1, 1, 1, 1.0).unwrap();
    let state = carrier(&domain, &[1.0]);
    let formation = 2.0_f64.powi(54);
    let model = ReservoirModel::new(CaloricLaw::Affine([-formation, formation + 4.0]), 1e30);
    let mut control = open_controls();
    control.integrator.maximum_enthalpy_defect_j = 0.0;
    control.integrator.solver.maximum_mass_defect_kg = [0.0; 2];
    control.integrator.solver.maximum_partition_defect_kg = 0.0;
    let got = state
        .diffuse_species_with_enthalpy_and_boundaries_coupled_adaptive(
            &domain,
            &[[1.0, 0.0]],
            &[-formation],
            &[reservoir(0, 0, 0.5)],
            1.0,
            &model,
            control,
            &CancelGate::new_clock_free(),
        )
        .unwrap();
    assert_eq!(got.step.specific_enthalpy_j_kg, [2.0]);
    assert_eq!(got.step.temperature_k, [300.0]);
    assert_eq!(got.thermal_boundary.enthalpy_into_domain_j, formation);
    assert_eq!(got.step.report.interval_enthalpy_defect_j, 0.0);
    assert_eq!(got.step.report.absolute_enthalpy_defect_j, 0.0);
    assert_eq!(got.step.report.accepted_substeps, 2);
}

#[test]
fn coupled_temporal_heat_error_survives_identical_large_datum_endpoint_readbacks() {
    let domain = VoxelDomain::new(3, 1, 1, 1.0).unwrap();
    let state = carrier(&domain, &[1.0; 3]);
    let datum = 2.0_f64.powi(80);
    let contrast = 2.0_f64.powi(40);
    let perturbation = 2.0_f64.powi(-20);
    let fractions = [0.5 + perturbation, 0.5, 0.5 - perturbation];
    let species = fractions.map(|y| [y, 1.0 - y]);
    let h = [datum + 0.5 * contrast; 3];
    let model = ReservoirModel::new(CaloricLaw::Affine([datum + contrast, datum]), 0.2);
    let gate = CancelGate::new_clock_free();
    let mut loose = open_controls();
    loose.integrator.absolute_temperature_tolerance_k = 1e9;
    loose.integrator.absolute_specific_enthalpy_tolerance_j_kg = 1e9;
    loose.integrator.absolute_mass_fraction_tolerance = [1.0; 2];
    // Deliberately admit the inability to store sub-ulp heat at this huge
    // datum. Physical representation error remains separate from time error.
    loose.integrator.maximum_enthalpy_defect_j = 1e9;
    let mut tight = loose;
    tight.integrator.absolute_specific_enthalpy_tolerance_j_kg = 1.0;
    let run = |control| {
        state
            .diffuse_species_with_enthalpy_and_boundaries_coupled_adaptive(
                &domain,
                &species,
                &h,
                &[],
                1.0,
                &model,
                control,
                &gate,
            )
            .unwrap()
    };
    let coarse = run(loose);
    let fine = run(tight);
    assert_eq!(coarse.step.specific_enthalpy_j_kg, h);
    assert_eq!(fine.step.specific_enthalpy_j_kg, h);
    assert!(coarse.step.report.estimated_specific_enthalpy_error_j_kg > 1.0);
    assert!(fine.step.report.estimated_specific_enthalpy_error_j_kg <= 1.0);
    assert!(fine.step.report.accepted_substeps > coarse.step.report.accepted_substeps);
    // Independent three-cell eigenmode: [eps,0,-eps] decays at rateD.
    // Its actual transported heat is resolvable even though all h readbacks
    // remain identical to their initial datum at every temporal resolution.
    let heat = contrast * perturbation * (-0.2_f64).exp_m1();
    near(fine.step.diffusive_enthalpy_change_j[0], heat, 3.0);
    near(fine.step.diffusive_enthalpy_change_j[1], 0.0, 3.0);
    near(fine.step.diffusive_enthalpy_change_j[2], -heat, 3.0);
    assert!(fine.step.report.absolute_enthalpy_defect_j > 0.0);
}

#[test]
fn concurrent_frozen_boundary_api_resolves_competition_and_retains_isolated_fast_path() {
    use super::super::enthalpy::boundary::adaptive::BoundarySpeciesEnthalpyAccuracy;
    let domain = VoxelDomain::new(1, 1, 1, 1.0).unwrap();
    let state = carrier(&domain, &[1.0]);
    let boundaries = [reservoir(0, 0, 0.0), reservoir(0, 1, 0.0)];
    let gate = CancelGate::new_clock_free();
    let accuracy = BoundarySpeciesEnthalpyAccuracy {
        absolute_boundary_transfer_tolerance_kg: [1e-5; 2],
        absolute_specific_enthalpy_tolerance_j_kg: 1e-4,
        absolute_boundary_enthalpy_tolerance_j: 1e-4,
    };
    let mut limited = AdaptiveSpeciesDiffusionControls {
        solver: controls(1),
        absolute_mass_fraction_tolerance: [1e-8; 2],
        relative_tolerance: 0.0,
        maximum_refinements: 1,
    };
    let run = |diffusivity, controls| {
        state.diffuse_species_with_enthalpy_and_boundaries_adaptive(
            &domain,
            &[[1.0, 0.0]],
            &[10.0],
            &[diffusivity],
            &boundaries,
            1.0,
            |_, _| unreachable!(),
            |index| {
                Ok(if index == 0 {
                    [10.0, 0.0]
                } else {
                    [30.0, 10.0]
                })
            },
            controls,
            accuracy,
            1e-7,
            &gate,
        )
    };
    assert!(run(1000.0, limited).is_err());
    limited.maximum_refinements = 16;
    let got = run(2.0, limited).unwrap();
    let dm = -0.5 * (1.0 - (-8.0_f64).exp());
    near(got.step.boundary.faces[0].mass_into_cell_kg[0], dm, 3e-5);
    near(got.step.boundary.faces[1].mass_into_cell_kg[0], dm, 3e-5);
    near(
        got.step.thermal_boundary.faces[0].enthalpy_into_cell_j,
        10.0 * dm,
        3e-4,
    );
    near(
        got.step.thermal_boundary.faces[1].enthalpy_into_cell_j,
        20.0 * dm,
        3e-4,
    );
    let isolated = state
        .diffuse_species_with_enthalpy_and_boundaries_adaptive(
            &domain,
            &[[1.0, 0.0]],
            &[10.0],
            &[1000.0],
            &[reservoir(0, 0, 0.0)],
            1.0,
            |_, _| unreachable!(),
            |_| Ok([10.0, 0.0]),
            limited,
            accuracy,
            0.0,
            &gate,
        )
        .unwrap();
    assert_eq!(isolated.refinements, 1);
    assert_eq!(isolated.step.step.specific_enthalpy_j_kg, [0.0]);
    assert_eq!(isolated.step.thermal_boundary.enthalpy_into_domain_j, -10.0);
}

#[test]
fn adaptive_reservoir_changing_partial_law_cannot_alias_a_stiff_thermal_transient() {
    let domain = VoxelDomain::new(1, 1, 1, 1.0).unwrap();
    let state = carrier(&domain, &[1.0]);
    let boundaries = [reservoir(0, 0, 0.5)];
    let gate = CancelGate::new_clock_free();
    for inner in [1, 4096] {
        let model = ReservoirModel::new(CaloricLaw::Cubic, 1000.0);
        let mut control = open_controls();
        control.integrator.solver.substeps = inner;
        control.integrator.maximum_refinements = 2;
        let error = state
            .diffuse_species_with_enthalpy_and_boundaries_coupled_adaptive(
                &domain,
                &[[1.0, 0.0]],
                &[4.0],
                &boundaries,
                1.0,
                &model,
                control,
                &gate,
            )
            .unwrap_err();
        assert!(reason(error).contains("relaxation"));
    }
    let model = ReservoirModel::new(CaloricLaw::Cubic, 4.0);
    let got = state
        .diffuse_species_with_enthalpy_and_boundaries_coupled_adaptive(
            &domain,
            &[[1.0, 0.0]],
            &[4.0],
            &boundaries,
            1.0,
            &model,
            open_controls(),
            &gate,
        )
        .unwrap();
    // Y=.5+.5exp(-2Dt), dh=12Y²dY: h=4Y³ and T=300 for the
    // complete continuous interval. A saturated unresolved midpoint gives302K.
    let y = 0.5 + 0.5 * (-8.0_f64).exp();
    near(got.step.species_mass_kg[0][0], y, 2e-12);
    near(got.step.specific_enthalpy_j_kg[0], 4.0 * y.powi(3), 3e-4);
    near(got.step.temperature_k[0], 300.0, 3e-4);
    assert!(got.step.report.maximum_changing_source_relaxation <= 0.5);
}

#[test]
fn adaptive_competing_reservoir_maps_cannot_alias_frozen_face_heat_assignment() {
    let domain = VoxelDomain::new(1, 1, 1, 1.0).unwrap();
    let state = carrier(&domain, &[1.0]);
    let boundaries = [reservoir(0, 0, 0.0), reservoir(0, 1, 0.0)];
    let gate = CancelGate::new_clock_free();
    let model = ReservoirModel::new(CaloricLaw::CompetingFaces, 1000.0);
    let mut limited = open_controls();
    limited.integrator.maximum_refinements = 1;
    let error = state
        .diffuse_species_with_enthalpy_and_boundaries_coupled_adaptive(
            &domain,
            &[[1.0, 0.0]],
            &[10.0],
            &boundaries,
            1.0,
            &model,
            limited,
            &gate,
        )
        .unwrap_err();
    assert!(reason(error).contains("relaxation"));
    // A single active reservoir per cell can also compete with an interior
    // map. Exact long-time A exits these equal-D faces in proportions3/4,1/4;
    // saturated boundary-before-interior sweeps incorrectly assign all to0.
    let pair_domain = VoxelDomain::new(2, 1, 1, 1.0).unwrap();
    let pair = carrier(&pair_domain, &[1.0; 2]);
    let pair_model = ReservoirModel::new(CaloricLaw::Affine([8.0, 0.0]), 1000.0);
    let error = pair
        .diffuse_species_with_enthalpy_and_boundaries_coupled_adaptive(
            &pair_domain,
            &[[1.0, 0.0], [0.0, 1.0]],
            &[8.0, 0.0],
            &[reservoir(0, 0, 0.0), reservoir(1, 1, 0.0)],
            1.0,
            &pair_model,
            limited,
            &gate,
        )
        .unwrap_err();
    assert!(reason(error).contains("relaxation"));

    let model = ReservoirModel::new(CaloricLaw::CompetingFaces, 2.0);
    let got = state
        .diffuse_species_with_enthalpy_and_boundaries_coupled_adaptive(
            &domain,
            &[[1.0, 0.0]],
            &[10.0],
            &boundaries,
            1.0,
            &model,
            open_controls(),
            &gate,
        )
        .unwrap();
    // Two equal G=2D reservoirs each receive half the lost A. Original
    // partial contrasts10/20 give Q=−15(1−exp(−4Dt)), h=−5+15exp(−4Dt).
    let y = (-8.0_f64).exp();
    let dm = -0.5 * (1.0 - y);
    near(got.step.species_mass_kg[0][0], y, 2e-12);
    near(got.boundary.faces[0].mass_into_cell_kg[0], dm, 3e-5);
    near(got.boundary.faces[1].mass_into_cell_kg[0], dm, 3e-5);
    near(
        got.thermal_boundary.faces[0].enthalpy_into_cell_j,
        10.0 * dm,
        3e-4,
    );
    near(
        got.thermal_boundary.faces[1].enthalpy_into_cell_j,
        20.0 * dm,
        3e-4,
    );
    near(got.step.specific_enthalpy_j_kg[0], -5.0 + 15.0 * y, 3e-4);
    near(
        got.step.temperature_k[0],
        295.0 + 10.0 * y / (1.0 + y),
        3e-4,
    );
    assert!(got.step.report.maximum_competing_boundary_relaxation <= 0.5);

    // Explicitly resolved inner sweeps need no unnecessary outer-property
    // resolution for these unchanged sources. Nonzero reservoir traces avoid
    // asking the original positive map to represent an underflowing species.
    let model = ReservoirModel::new(CaloricLaw::CompetingFaces, 1000.0);
    let mut resolved = open_controls();
    resolved.integrator.solver.substeps = 4096;
    resolved.integrator.maximum_refinements = 1;
    resolved.integrator.absolute_temperature_tolerance_k = 1.0;
    resolved
        .integrator
        .absolute_specific_enthalpy_tolerance_j_kg = 1.0;
    resolved.absolute_boundary_transfer_tolerance_kg = [0.1; 2];
    resolved.absolute_boundary_enthalpy_tolerance_j = 1.0;
    let got = state
        .diffuse_species_with_enthalpy_and_boundaries_coupled_adaptive(
            &domain,
            &[[1.0, 0.0]],
            &[10.0],
            &[reservoir(0, 0, 0.25), reservoir(0, 1, 0.25)],
            1.0,
            &model,
            resolved,
            &gate,
        )
        .unwrap();
    assert_eq!(got.step.report.accepted_substeps, 2);
    assert_eq!(got.step.report.maximum_changing_source_relaxation, 0.0);
    assert!(got.step.report.maximum_competing_boundary_relaxation <= 0.5);
    near(got.step.specific_enthalpy_j_kg[0], -1.25, 0.05);
}

#[test]
fn adaptive_reservoir_source_and_total_resource_refusals_leave_inputs_unchanged() {
    let domain = VoxelDomain::new(1, 1, 1, 1.0).unwrap();
    let state = carrier(&domain, &[1.0]);
    let species = [[0.75, 0.25]];
    let h = [6.0];
    let boundaries = [reservoir(0, 0, 0.25)];
    let gate = CancelGate::new_clock_free();
    let model = ReservoirModel::new(CaloricLaw::Affine([8.0, 0.0]), 0.2);
    let accepted = state
        .diffuse_species_with_enthalpy_and_boundaries_coupled_adaptive(
            &domain,
            &species,
            &h,
            &boundaries,
            1.0,
            &model,
            open_controls(),
            &gate,
        )
        .unwrap();
    for budget in 0..3 {
        let model = ReservoirModel::new(CaloricLaw::Affine([8.0, 0.0]), 0.2);
        let mut control = open_controls();
        let expected = match budget {
            0 => {
                control.integrator.solver.maximum_storage_bytes = 0;
                "storage"
            }
            1 => {
                control.integrator.solver.maximum_pair_updates =
                    accepted.step.report.total_pair_updates - 1;
                "pair"
            }
            _ => {
                control.integrator.maximum_source_evaluations =
                    accepted.step.report.source_preparations
                        + accepted.step.report.face_source_evaluations
                        - 1;
                "source"
            }
        };
        let error = state
            .diffuse_species_with_enthalpy_and_boundaries_coupled_adaptive(
                &domain,
                &species,
                &h,
                &boundaries,
                1.0,
                &model,
                control,
                &gate,
            )
            .unwrap_err();
        assert!(reason(error).contains(expected));
        if budget == 0 {
            assert_eq!(model.preparations.get(), 0);
        } else {
            assert!(model.boundary_queries.get() > 0);
        }
    }
    let mut model = ReservoirModel::new(CaloricLaw::Affine([8.0, 0.0]), 0.2);
    model.fail_at_boundary = 2;
    let error = state
        .diffuse_species_with_enthalpy_and_boundaries_coupled_adaptive(
            &domain,
            &species,
            &h,
            &boundaries,
            1.0,
            &model,
            open_controls(),
            &gate,
        )
        .unwrap_err();
    assert!(reason(error).contains("authored adaptive boundary source refused"));
    let mut model = ReservoirModel::new(CaloricLaw::Affine([8.0, 0.0]), 0.2);
    model.cancel_at_boundary = 2;
    assert!(
        state
            .diffuse_species_with_enthalpy_and_boundaries_coupled_adaptive(
                &domain,
                &species,
                &h,
                &boundaries,
                1.0,
                &model,
                open_controls(),
                &gate,
            )
            .is_err()
    );
    assert_eq!(model.boundary_queries.get(), 2);
    assert_eq!(state.mass_kg, [1.0]);
    assert_eq!(species, [[0.75, 0.25]]);
    assert_eq!(h, [6.0]);
}
