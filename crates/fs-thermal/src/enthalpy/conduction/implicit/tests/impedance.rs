//! G1/G3/G4: original caloric and conductive laws, analytical impedance and
//! direct native time-domain refinement. Synthetic fixtures are not validation.
use super::super::harmonic::{ThermalImpedanceBudget, ThermalImpedanceError};
use super::super::source::{ImplicitHeatSourceLaw, ImplicitHeatSourceState};
use super::*;
use fs_math::c64::C64;

fn budget() -> ThermalImpedanceBudget {
    ThermalImpedanceBudget {
        maximum_cells: 16,
        maximum_dense_bytes: 100_000,
        maximum_dense_work: 100_000,
        equilibrium_tolerance_w: 1e-9,
        relative_residual_tolerance: 1e-12,
    }
}
fn complex_near(actual: C64, expected: C64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:?} != {expected:?}, tolerance {tolerance:e}"
    );
}

struct Caloric;
impl ThermalLaw for Caloric {
    fn temperature_k(&self, h: f64) -> Result<f64, DepositionError> {
        if !(300.0..=3000.0).contains(&h) {
            return Err(DepositionError::MaterialDomain);
        }
        Ok(h / 3.0)
    }
    fn maximum_temperature_slope(&self) -> f64 {
        1.0 / 3.0
    }
}
impl ConductivityLaw for Caloric {
    fn conductivity_w_m_k(&self, _: f64) -> Result<f64, DepositionError> {
        Ok(2.0)
    }
}
impl ImplicitConductionLaw for Caloric {
    fn caloric_law(&self) -> &dyn ThermalLaw {
        self
    }
    fn conductivity_law(&self) -> &dyn ConductivityLaw {
        self
    }
    fn enthalpy_bounds_j_kg(&self) -> [f64; 2] {
        [300.0, 3000.0]
    }
    fn state_at_specific_enthalpy(
        &self,
        h: f64,
    ) -> Result<ImplicitConductionState, DepositionError> {
        Ok(ImplicitConductionState {
            temperature_k: self.temperature_k(h)?,
            temperature_slope_kg_k_j: 1.0 / 3.0,
            conductivity_w_m_k: 2.0,
            conductivity_slope_kg_w_m_k_j: 0.0,
        })
    }
}

#[test]
fn thermal_impedance_uses_original_capacity_and_physical_contact_resistance() {
    let law = Caloric;
    let (mut thermal, _) = setup(&law, &[900.0], Some(300.0), Vec::new());
    thermal.cells.get_mut(&CellId(0)).unwrap().mass_kg = 2.0;
    // A=1, d/k=1/2 and R_surface=1/2: G=1 W/K; m=2, cp=3 gives C=6 J/K.
    let conduction = ConductionObserver::new(
        &thermal,
        &[(MaterialId(0), &law)],
        &[],
        &[ConductiveReservoir {
            cell_id: CellId(0),
            half_distance_m: 1.0,
            surface_resistance_m2_k_w: 0.5,
        }],
        || Ok(()),
    )
    .unwrap();
    let chain =
        ImplicitConductionChain::new(&thermal, &conduction, &[(MaterialId(0), &law)], || Ok(()))
            .unwrap();
    let model = chain
        .prepare_thermal_impedance(&thermal, budget(), || Ok(()))
        .unwrap();
    assert_eq!(model.cell_ids(), &[CellId(0)]);
    assert_eq!(model.operating_temperatures_k(), &[300.0]);
    assert_eq!(model.operating_time_s(), thermal.time_s());
    for omega in [0.0, 1e-10, 0.7, 1e4] {
        let factor = model.factor(omega, || Ok(())).unwrap();
        let expected = C64::new(1.0, -6.0 * omega).recip();
        let z = factor.impedance_column(CellId(0), || Ok(())).unwrap()[0];
        complex_near(z, expected, 2e-15);
        let response = factor.response(&[C64::new(2.0, -0.7)], || Ok(())).unwrap();
        complex_near(
            response.temperature_k[0],
            expected * C64::new(2.0, -0.7),
            4e-15,
        );
        complex_near(
            response.energy_j[0],
            response.temperature_k[0].scale(6.0),
            4e-15,
        );
        assert!(response.relative_residual < 1e-14);
    }
    assert_eq!(thermal.time_s(), 0.0);
    assert_eq!(thermal.cells()[&CellId(0)].specific_h_j_kg, 900.0);
}

#[test]
fn thermal_mutual_impedance_preserves_coupled_diffusion_and_passive_quadratic_form() {
    let law = Law {
        nonlinear_k: false,
        plateau: false,
    };
    let (thermal, conduction) = setup(&law, &[300.0; 2], Some(300.0), Vec::new());
    let graph = ImplicitConductionGraph::new(
        &thermal,
        &conduction,
        &[(MaterialId(0), &law)],
        ImplicitConductionLinearControls::default(),
        || Ok(()),
    )
    .unwrap();
    let model = graph
        .prepare_thermal_impedance(&thermal, budget(), || Ok(()))
        .unwrap();
    for omega in [0.0, 0.1, 1.0, 10.0, 1e3] {
        let factor = model.factor(omega, || Ok(())).unwrap();
        let a = C64::new(2.0, -omega);
        let d = C64::new(1.0, -omega);
        let determinant = a * d - C64::ONE;
        let first = factor.impedance_column(CellId(0), || Ok(())).unwrap();
        let second = factor.impedance_column(CellId(1), || Ok(())).unwrap();
        complex_near(first[0], d / determinant, 2e-14);
        complex_near(first[1], C64::ONE / determinant, 2e-14);
        complex_near(second[1], a / determinant, 2e-14);
        complex_near(second[0], first[1], 2e-14);
        let p = [C64::new(0.3, -0.6), C64::new(0.7, 0.2)];
        let response = factor.response(&p, || Ok(())).unwrap();
        let t = &response.temperature_k;
        // W*K, not physical heat: exact positive-real identity of this passive
        // fixed-coefficient graph, with G_bath=G_link=1 W/K.
        let quadratic = (p[0].conj() * t[0] + p[1].conj() * t[1]).re;
        let conduction = t[0].norm_sq() + (t[0] - t[1]).norm_sq();
        assert!(quadratic > 0.0);
        near(quadratic, conduction, 2e-14);
    }
}

struct Source<'a> {
    law: &'a Law,
    sinusoidal: bool,
}
impl ImplicitHeatSourceLaw for Source<'_> {
    fn caloric_law(&self) -> &dyn ThermalLaw {
        self.law
    }
    fn state_at_specific_enthalpy(
        &self,
        h: f64,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<ImplicitHeatSourceState, DepositionError> {
        self.state_at_specific_enthalpy_and_time(h, 0.0, checkpoint)
    }
    fn state_at_specific_enthalpy_and_time(
        &self,
        h: f64,
        time: f64,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<ImplicitHeatSourceState, DepositionError> {
        checkpoint()?;
        Ok(if self.sinusoidal {
            ImplicitHeatSourceState {
                power_w: 0.2 * time.sin(),
                power_enthalpy_derivative_w_kg_j: 0.0,
            }
        } else {
            ImplicitHeatSourceState {
                power_w: 1800.0 + 2.0 * (h - 600.0),
                power_enthalpy_derivative_w_kg_j: 2.0,
            }
        })
    }
}

#[test]
fn thermal_impedance_differentiates_original_nonlinear_heat_and_source_feedback() {
    let law = Law {
        nonlinear_k: true,
        plateau: false,
    };
    let source = Source {
        law: &law,
        sinusoidal: false,
    };
    let (thermal, conduction) = setup(&law, &[600.0], Some(300.0), Vec::new());
    let chain =
        ImplicitConductionChain::new(&thermal, &conduction, &[(MaterialId(0), &law)], || Ok(()))
            .unwrap();
    assert!(matches!(
        chain.prepare_thermal_impedance(&thermal, budget(), || Ok(())),
        Err(ThermalImpedanceError::NotEquilibrium { .. })
    ));
    let chain = chain
        .with_heat_sources(&thermal, &[(CellId(0), &source)], || Ok(()))
        .unwrap();
    let model = chain
        .prepare_thermal_impedance(&thermal, budget(), || Ok(()))
        .unwrap();
    // Q_out=(T/100)*(T-300), so dQ/dT=9 at T=600. Source dP/dT=2.
    // Freezing just the secant conductance would incorrectly give 6-2=4.
    let z = model
        .factor(0.4, || Ok(()))
        .unwrap()
        .impedance_column(CellId(0), || Ok(()))
        .unwrap()[0];
    complex_near(z, C64::new(7.0, -0.4).recip(), 2e-15);
}

#[test]
fn thermal_impedance_agrees_with_refined_native_two_cell_time_response() {
    let law = Law {
        nonlinear_k: false,
        plateau: false,
    };
    let source = Source {
        law: &law,
        sinusoidal: true,
    };
    let (thermal, conduction) = setup(&law, &[300.0; 2], Some(300.0), Vec::new());
    let chain =
        ImplicitConductionChain::new(&thermal, &conduction, &[(MaterialId(0), &law)], || Ok(()))
            .unwrap()
            .with_heat_sources(&thermal, &[(CellId(1), &source)], || Ok(()))
            .unwrap();
    let model = chain
        .prepare_thermal_impedance(&thermal, budget(), || Ok(()))
        .unwrap();
    let response = model
        .factor(1.0, || Ok(()))
        .unwrap()
        .response(&[C64::ZERO, C64::new(0.0, 0.2)], || Ok(()))
        .unwrap();
    let time = 40.0_f64; // Decaying homogeneous mode is below 1e-7 K here.
    let phase = C64::new(time.cos(), -time.sin());
    let mut errors = Vec::new();
    for dt in [0.1, 0.05, 0.025] {
        let mut workspace = ImplicitConductionWorkspace::new(2).unwrap();
        let (next, _) = chain
            .advance_to(&thermal, time, controls(dt), &mut workspace, || Ok(()))
            .unwrap();
        let error = (0..2)
            .map(|i| {
                (next.temperature_k(CellId(i as u64)).unwrap()
                    - 300.0
                    - (response.temperature_k[i] * phase).re)
                    .abs()
            })
            .fold(0.0, f64::max);
        errors.push(error);
    }
    assert!(
        errors[1] < 0.6 * errors[0] && errors[2] < 0.6 * errors[1],
        "{errors:?}"
    );
}

#[test]
fn thermal_impedance_preserves_latent_energy_and_refuses_ungrounded_dc() {
    for plateau in [false, true] {
        let law = Law {
            nonlinear_k: false,
            plateau,
        };
        let h = if plateau { 400.0 } else { 300.0 };
        let (thermal, conduction) = setup(&law, &[h; 2], None, Vec::new());
        let chain = ImplicitConductionChain::new(
            &thermal,
            &conduction,
            &[(MaterialId(0), &law)],
            || Ok(()),
        )
        .unwrap();
        let model = chain
            .prepare_thermal_impedance(&thermal, budget(), || Ok(()))
            .unwrap();
        assert!(model.factor(0.0, || Ok(())).is_err());
        let factor = model.factor(2.0, || Ok(())).unwrap();
        let response = factor.response(&[C64::ONE; 2], || Ok(())).unwrap();
        for i in 0..2 {
            complex_near(response.energy_j[i], C64::new(0.0, 0.5), 1e-15);
            if plateau {
                assert_eq!(response.temperature_k[i], C64::ZERO);
            }
        }
        assert!(factor.impedance_column(CellId(99), || Ok(())).is_err());
        assert!(
            factor
                .response(&[C64::from_re(f64::NAN); 2], || Ok(()))
                .is_err()
        );
        let mut calls = 0;
        assert!(
            factor
                .response(&[C64::ONE; 2], || {
                    calls += 1;
                    if calls == 6 {
                        Err(DepositionError::Cancelled)
                    } else {
                        Ok(())
                    }
                })
                .is_err()
        );
        assert_eq!(
            factor.response(&[C64::ONE; 2], || Ok(())).unwrap(),
            response
        );
        assert_eq!(thermal.time_s(), 0.0);
    }
}
