//! Prescribed ghost states for the hybrid NSF lattice and energy equation.
//!
//! A ghost centre lies one cell spacing outside the adjacent fluid centre;
//! their common face is halfway between them. Incoming kinetic populations use
//! non-equilibrium extrapolation (Guo, Zheng & Shi, Chinese Physics 11 (2002),
//! 366–374, doi:10.1088/1009-1963/11/4/310). The thermal-pressure force and total
//! energy flux use the SAME ghost state. This is a boundary primitive, not an
//! acoustic non-reflection condition. In particular a prescribed state does not
//! preserve the outgoing characteristic of an impinging sound wave.

use super::{NsfError, NsfGrid3, THETA0, check_state};
use crate::lattice3::Lattice3;

/// Physical state at an exterior ghost centre, held until explicitly changed.
///
/// This specifies incoming populations, not an exact Dirichlet condition on
/// the adjacent fluid centre. Values must come from the caller's physical
/// boundary or a coupled solver; no pressure, temperature or velocity defaults.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct NsfOpenState {
    /// Mass density [kg/m³].
    pub density_kg_m3: f64,
    /// Velocity in the grid's Cartesian frame [m/s].
    pub velocity_m_s: [f64; 3],
    /// Absolute temperature [K].
    pub temperature_k: f64,
}

#[derive(Clone, Copy)]
pub(super) struct OpenGhost {
    pub rho: f64,
    pub u: [f64; 3],
    pub theta: f64,
    pub energy: f64,
}

impl<L: Lattice3<27>> NsfGrid3<L> {
    /// Replace one non-periodic face with a prescribed exterior state.
    /// `None` restores its stationary wall. Solid cells remain solid.
    ///
    /// This first boundary rung admits one open axis with periodic transverse
    /// domain edges. Internal solids may describe channels/nozzles, but may not
    /// intersect the open plane: solid/open diagonal-link corners are unqualified.
    /// It deliberately refuses intersecting open/wall domain edges rather than
    /// choosing a corner rule by axis order. It does not claim transparency,
    /// characteristic coupling, jet validation or a converged open-boundary order.
    ///
    /// Updating the ghost state restarts force/energy-flux history while
    /// preserving the current fluid's density, physical velocity and energy.
    /// A refusal leaves the boundary and all published state unchanged.
    ///
    /// # Errors
    /// Invalid axis/topology, a moving face, manufactured sources, or an
    /// exterior state outside the same material/Mach/temperature bounds as the grid.
    pub fn set_open_face(
        &mut self,
        axis: usize,
        high: bool,
        state: Option<NsfOpenState>,
    ) -> Result<(), NsfError> {
        if self.sgs.is_some() {
            return Err(NsfError::Input(
                "NSF SGS does not support prescribed-reservoir faces",
            ));
        }
        if self.open.is_some() {
            return Err(NsfError::Input(
                "prescribed reservoir and characteristic faces are exclusive",
            ));
        }
        if axis >= 3 || self.spec.periodic[axis] || self.spec.dims[axis] < 2 {
            return Err(NsfError::Input(
                "open face requires a non-periodic axis with at least two cells",
            ));
        }
        if (0..3).any(|a| a != axis && !self.spec.periodic[a]) {
            return Err(NsfError::Input(
                "open-face transverse domain edges must be periodic; use interior solids for walls",
            ));
        }
        let plane = if high { self.spec.dims[axis] - 1 } else { 0 };
        if state.is_some()
            && self
                .solid
                .iter()
                .enumerate()
                .any(|(i, solid)| *solid && self.at(i)[axis] == plane)
        {
            return Err(NsfError::Input(
                "open plane must be free of solids; solid/open junctions are unqualified",
            ));
        }
        let side = usize::from(high);
        if self.edge_velocity[axis][side] != [0.0; 3] || self.sources.is_some() {
            return Err(NsfError::Input(
                "open face cannot combine with a moving face or manufactured sources",
            ));
        }
        let ghost = state.map(|s| self.admit_open_state(s)).transpose()?;
        let previous = self.reservoir_faces[axis][side];
        let (rho, mom) = self.macro_fields(&self.f);
        let u = self.velocity_field(&rho, &mom, &self.force);
        self.reservoir_faces[axis][side] = ghost;
        let force = self.force_field(&rho, &self.theta, None);
        let mut populations = self.f.clone();
        for (i, f) in populations.iter_mut().enumerate() {
            if self.solid[i] {
                continue;
            }
            // Preserve the whole physical distribution, including its stress,
            // when changing the thermal-pressure force of a reservoir face.
            super::rebase_force::<L>(
                f,
                u[i],
                std::array::from_fn(|a| self.force[i][a] - force[i][a]),
            );
        }
        let sig = self.stress_field(&populations, &force, &rho, &self.theta, &self.material, &[]);
        let history = self.flux_divergence(&self.e, &rho, &u, &self.theta, &sig, &[], None);
        let (div, wall, open) = match history {
            Ok(value) => value,
            Err(error) => {
                self.reservoir_faces[axis][side] = previous;
                return Err(error);
            }
        };
        self.f = populations;
        self.force = force;
        self.div = div;
        self.wall_out = wall;
        self.open_out = open;
        Ok(())
    }

    fn admit_open_state(&self, s: NsfOpenState) -> Result<OpenGhost, NsfError> {
        let rho = s.density_kg_m3 / self.law.rho_ref_kg_m3;
        let u = s.velocity_m_s.map(|v| v * self.law.dt_s / self.law.dx_m);
        let theta = self.law.theta_of_kelvin(s.temperature_k);
        if !(rho.is_finite()
            && rho > 0.0
            && u.iter().all(|v| v.is_finite())
            && theta.is_finite()
            && theta > 0.0)
        {
            return Err(NsfError::Input(
                "open-face density, velocity and temperature must be finite and physical",
            ));
        }
        let [low_k, high_k] = self.spec.temperature_window_k;
        if s.temperature_k < low_k || s.temperature_k > high_k {
            return Err(NsfError::Input(
                "open-face temperature outside the declared run window",
            ));
        }
        let material = self.spec.gas.resolve(&self.law, rho, theta)?;
        let gas = self
            .spec
            .gas
            .at(s.temperature_k, self.law.pressure_pa(rho, theta))?;
        check_state(&self.spec, gas.gamma, self.n, rho, u, theta)?;
        let internal = material
            .as_ref()
            .map_or(theta / (self.law.gamma - 1.0), |m| {
                m.sensible_energy_lattice
            });
        let energy = rho * (internal + 0.5 * u.iter().map(|v| v * v).sum::<f64>());
        if !energy.is_finite() {
            return Err(NsfError::Input("open-face energy is unrepresentable"));
        }
        Ok(OpenGhost {
            rho,
            u,
            theta,
            energy,
        })
    }

    pub(super) fn open_at(&self, i: usize, axis: usize, direction: i64) -> Option<OpenGhost> {
        if self.shifted(self.at(i), axis, direction).is_none() {
            self.reservoir_faces[axis][usize::from(direction > 0)]
        } else {
            None
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::d3q27::D3q27;
    use crate::nsf3::{MAX_QUALIFIED_MACH, NsfSpec, WallThermal};
    use fs_material::gas::{GasSpec, GasState};

    fn spec(axis: usize, workers: usize) -> NsfSpec {
        let mut dims = [3; 3];
        dims[axis] = 24;
        let mut periodic = [true; 3];
        periodic[axis] = false;
        NsfSpec {
            dims,
            periodic,
            gas: GasSpec::dry_air_ussa1976(),
            temperature_k: 293.15,
            pressure_pa: 101325.0,
            dx_m: 5e-5,
            temperature_window_k: [273.15, 333.15],
            maximum_mach: MAX_QUALIFIED_MACH,
            maximum_temperature_ratio: 0.2,
            bulk_to_shear_ratio: 0.73,
            bulk_ratio_source: "Cramer 2012, Phys. Fluids 24 066102, N2 at 293 K",
            maximum_cells: 216,
            wall: WallThermal::Adiabatic,
            workers,
            transport_points: 33,
        }
    }

    fn state(velocity_m_s: [f64; 3], temperature_k: f64) -> NsfOpenState {
        let gas = GasState::try_new(&GasSpec::dry_air_ussa1976(), temperature_k, 101325.0).unwrap();
        NsfOpenState {
            density_kg_m3: gas.density,
            velocity_m_s,
            temperature_k,
        }
    }

    fn uniform(axis: usize, workers: usize, s: NsfOpenState) -> NsfGrid3<D3q27> {
        let spec = spec(axis, workers);
        let reference = GasState::try_new(&spec.gas, spec.temperature_k, spec.pressure_pa).unwrap();
        let speed_unit = (3.0 * reference.specific_gas_constant * reference.temperature).sqrt();
        NsfGrid3::new(
            spec,
            |_| false,
            |_| {
                (
                    s.density_kg_m3 / reference.density,
                    s.velocity_m_s.map(|u| u / speed_unit),
                    s.temperature_k,
                )
            },
        )
        .unwrap()
    }

    #[test]
    fn g0_open_faces_preserve_uniform_throughflow_all_axes() {
        for axis in 0..3 {
            for speed in [-10.0, 10.0] {
                let mut velocity = [0.0; 3];
                velocity[axis] = speed;
                let s = state(velocity, 303.15);
                let mut grid = uniform(axis, 1, s);
                let expected = grid.cell([1; 3]);
                for high in [false, true] {
                    grid.set_open_face(axis, high, Some(s)).unwrap();
                }
                let initial_energy = grid.total_energy_j();
                for _ in 0..32 {
                    let ledger = grid.step().unwrap();
                    assert!(ledger.residual_rel < 1e-12);
                }
                let mut max_error: f64 = 0.0;
                for i in 0..grid.n {
                    let (rho, u, temperature) = grid.cell(grid.at(i));
                    max_error = max_error
                        .max((rho - expected.0).abs())
                        .max((temperature / expected.2 - 1.0).abs());
                    for a in 0..3 {
                        max_error = max_error.max((u[a] - expected.1[a]).abs());
                    }
                }
                assert!(
                    max_error < 1e-11,
                    "axis={axis} speed={speed} error={max_error:e}"
                );
                assert!((grid.total_energy_j() / initial_energy - 1.0).abs() < 1e-12);
                // A sealed domain cannot satisfy this throughflow fixed point.
                let mut sealed = uniform(axis, 1, s);
                sealed.step().unwrap();
                assert!((sealed.cell([0; 3]).1[axis] - expected.1[axis]).abs() > 1e-4);
                println!("NSF_OPEN uniform axis={axis} speed_m_s={speed} error={max_error:e}");
            }
        }
    }

    #[test]
    fn g0_open_energy_flux_is_separate_from_wall_heat() {
        let cold = state([8.0, 0.0, 0.0], 293.15);
        let hot = state([8.0, 0.0, 0.0], 294.15);
        let mut grid = uniform(0, 1, cold);
        grid.set_open_face(0, false, Some(hot)).unwrap();
        grid.set_open_face(0, true, Some(cold)).unwrap();
        let before = grid.total_energy_j();
        for _ in 0..64 {
            let ledger = grid.step().unwrap();
            assert_eq!(ledger.wall_heat_j, 0.0);
            assert!(ledger.residual_rel < 1e-12, "{:?}", ledger);
        }
        let change = grid.total_energy_j() - before;
        assert!(change > 0.0);
        assert!(grid.cell([0, 1, 1]).2 > cold.temperature_k);
        let incoming = grid.cumulative_ledger().open_energy_j;
        assert!((change - incoming).abs() / before < 1e-12);
        println!(
            "NSF_OPEN heat change_j={change:e} open_energy_j={incoming:e} inlet_k={}",
            grid.cell([0, 1, 1]).2
        );
    }

    #[test]
    fn g0_open_refusals_and_updates_preserve_published_state() {
        let s = state([0.0; 3], 293.15);
        let mut grid = uniform(0, 1, s);
        let before = grid.state_hash();
        for axis in [1, 2, 3] {
            assert!(grid.set_open_face(axis, false, Some(s)).is_err());
        }
        for bad in [
            NsfOpenState {
                density_kg_m3: 0.0,
                ..s
            },
            NsfOpenState {
                density_kg_m3: f64::NAN,
                ..s
            },
            NsfOpenState {
                temperature_k: 400.0,
                ..s
            },
            NsfOpenState {
                temperature_k: f64::NAN,
                ..s
            },
            NsfOpenState {
                velocity_m_s: [1000.0, 0.0, 0.0],
                ..s
            },
            NsfOpenState {
                velocity_m_s: [f64::NAN, 0.0, 0.0],
                ..s
            },
        ] {
            assert!(grid.set_open_face(0, false, Some(bad)).is_err());
            assert_eq!(grid.state_hash(), before);
            assert!(grid.reservoir_faces[0][0].is_none());
        }
        let cells: Vec<_> = (0..grid.n).map(|i| grid.cell(grid.at(i))).collect();
        grid.set_open_face(0, false, Some(state([0.0; 3], 303.15)))
            .unwrap();
        for (i, expected) in cells.iter().enumerate() {
            let actual = grid.cell(grid.at(i));
            assert!((actual.0 - expected.0).abs() < 1e-14);
            assert!((actual.2 - expected.2).abs() < 1e-12);
            for a in 0..3 {
                assert!((actual.1[a] - expected.1[a]).abs() < 1e-14);
            }
        }
        assert!(grid.set_edge_velocity(0, false, [0.0, 1e-3, 0.0]).is_err());
        grid.set_open_face(0, false, None).unwrap();
        grid.set_edge_velocity(0, false, [0.0, 1e-3, 0.0]).unwrap();
        assert!(grid.set_open_face(0, false, Some(s)).is_err());
        let mut corner =
            NsfGrid3::<D3q27>::new(spec(0, 1), |x| x == [0; 3], |_| (1.0, [0.0; 3], 293.15))
                .unwrap();
        let before = corner.state_hash();
        assert!(corner.set_open_face(0, false, Some(s)).is_err());
        assert_eq!(corner.state_hash(), before);
        assert!(corner.reservoir_faces[0][0].is_none());
    }

    #[test]
    fn g0_reservoir_rebinding_preserves_physical_stress_and_energy() {
        use crate::nsf3::{Q, neq_second};

        for axis in 0..3 {
            let s = state([8.0, -4.0, 3.0], 293.15);
            let mut grid = uniform(axis, 1, s);
            // Nonzero shear and bulk stresses must survive the rebase too.
            // These Hermite coefficients have zero mass and first moment.
            let physical = [2e-4, -1e-4, 3e-4, 1.5e-4, -2.5e-4, 5e-5];
            let trace = physical[0] + physical[1] + physical[2];
            for f in &mut grid.f {
                for q in 0..Q {
                    let [x, y, z] = D3q27::E[q].map(f64::from);
                    let cpc = x * x * physical[0]
                        + y * y * physical[1]
                        + z * z * physical[2]
                        + 2.0 * (x * y * physical[3] + x * z * physical[4] + y * z * physical[5]);
                    f[q] += D3q27::W[q] * (cpc - THETA0 * trace) / (2.0 * THETA0.powi(2));
                }
            }
            let moments = |grid: &NsfGrid3<D3q27>| -> Vec<[f64; 6]> {
                (0..grid.n)
                    .map(|i| {
                        let (rho, u, _) = grid.cell(grid.at(i));
                        let eq = super::super::equilibrium::<D3q27>(rho, u, THETA0);
                        let mut second = neq_second::<D3q27>(&grid.f[i], &eq);
                        let f = grid.force[i];
                        second[0] += u[0] * f[0];
                        second[1] += u[1] * f[1];
                        second[2] += u[2] * f[2];
                        second[3] += 0.5 * (u[0] * f[1] + u[1] * f[0]);
                        second[4] += 0.5 * (u[0] * f[2] + u[2] * f[0]);
                        second[5] += 0.5 * (u[1] * f[2] + u[2] * f[1]);
                        second
                    })
                    .collect()
            };
            let before = moments(&grid);
            let cells: Vec<_> = (0..grid.n).map(|i| grid.cell(grid.at(i))).collect();
            let energy = grid.total_energy_j();
            for ghost in [
                Some(state(s.velocity_m_s, 303.15)),
                Some(state(s.velocity_m_s, 283.15)),
                None,
            ] {
                grid.set_open_face(axis, false, ghost).unwrap();
                let after = moments(&grid);
                for i in 0..grid.n {
                    for component in 0..6 {
                        assert!(
                            (after[i][component] - before[i][component]).abs() < 2e-14,
                            "axis={axis} cell={i} component={component}: {:?} -> {:?}",
                            before[i],
                            after[i]
                        );
                    }
                    let cell = grid.cell(grid.at(i));
                    assert!((cell.0 - cells[i].0).abs() < 2e-14);
                    assert_eq!(cell.2, cells[i].2);
                    for a in 0..3 {
                        assert!((cell.1[a] - cells[i].1[a]).abs() < 2e-14);
                    }
                }
                assert_eq!(grid.total_energy_j(), energy);
                if ghost.is_some() {
                    assert!(grid.force.iter().any(|f| f[axis].abs() > 1e-5));
                }
            }
        }
    }

    #[test]
    fn g0_local_density_controls_viscosity_admission_and_step_rollback() {
        use crate::nsf3::MAX_VERIFIED_LATTICE_VISCOSITY;

        let s = state([0.0; 3], 293.15);
        let mut grid = uniform(0, 1, s);
        let minimum_density = grid.law.mu_lattice(THETA0) / MAX_VERIFIED_LATTICE_VISCOSITY;
        let admitted = 1.01 * minimum_density;
        let refused = 0.99 * minimum_density;
        assert!(matches!(
            NsfGrid3::<D3q27>::new(spec(0, 1), |_| false, |_| (refused, [0.0; 3], 293.15)),
            Err(NsfError::Viscosity { .. })
        ));
        NsfGrid3::<D3q27>::new(spec(0, 1), |_| false, |_| (admitted, [0.0; 3], 293.15)).unwrap();

        // The same source law must reject a dilute prescribed reservoir
        // before publishing its boundary or any half-force rebase.
        let before = (grid.state_hash(), grid.total_energy_j());
        assert!(matches!(
            grid.set_open_face(
                0,
                false,
                Some(NsfOpenState {
                    density_kg_m3: refused * grid.law.rho_ref_kg_m3,
                    ..s
                })
            ),
            Err(NsfError::Viscosity { .. })
        ));
        assert_eq!((grid.state_hash(), grid.total_energy_j()), before);
        assert!(grid.reservoir_faces[0][0].is_none());

        // A constant-temperature manufactured rarefaction crosses the local
        // envelope after streaming. Match removed energy to removed mass so
        // that a temperature-domain check cannot mask the viscosity defect.
        let removed = refused - admitted;
        let cv_theta = THETA0 / (grid.law.gamma - 1.0);
        let mut periodic = spec(0, 1);
        periodic.periodic = [true; 3];
        let mut evolving = NsfGrid3::<D3q27>::new_with_sources(
            periodic,
            |_| false,
            |_| (admitted, [0.0; 3], 293.15),
            Some(Box::new(move |_, _| {
                (removed, [0.0; 3], removed * cv_theta)
            })),
        )
        .unwrap();
        let before = (
            evolving.state_hash(),
            evolving.steps(),
            evolving.cumulative_ledger(),
        );
        assert!(matches!(evolving.step(), Err(NsfError::Viscosity { .. })));
        assert_eq!(
            (
                evolving.state_hash(),
                evolving.steps(),
                evolving.cumulative_ledger()
            ),
            before
        );
    }

    #[test]
    fn g0_prescribed_and_characteristic_faces_are_mutually_exclusive() {
        use crate::nsf3::{NsfCharacteristicAxis, NsfCharacteristicBoundary};
        let s = state([0.0; 3], 293.15);
        let face = NsfCharacteristicBoundary {
            reference_pressure_pa: 101325.0,
            reference_temperature_k: s.temperature_k,
            reference_velocity_m_s: s.velocity_m_s,
            incoming_characteristic_m_s: 0.0,
        };
        let axis = NsfCharacteristicAxis {
            axis: 0,
            faces: [face; 2],
            maximum_relative_perturbation: 0.05,
        };
        let mut prescribed = uniform(0, 1, s);
        prescribed.set_open_face(0, false, Some(s)).unwrap();
        assert!(prescribed.with_characteristic_axis(axis).is_err());
        let mut characteristic = uniform(0, 1, s).with_characteristic_axis(axis).unwrap();
        let before = characteristic.state_hash();
        let traces = characteristic.open_boundary_traces(false);
        assert!(characteristic.set_open_face(0, false, Some(s)).is_err());
        assert!(characteristic.set_open_face(0, true, None).is_err());
        assert_eq!(characteristic.state_hash(), before);
        assert_eq!(characteristic.open_boundary_traces(false), traces);
    }

    #[test]
    fn g5_open_face_evolution_is_worker_count_invariant() {
        let s = state([8.0, 0.0, 0.0], 293.15);
        let hot = state([8.0, 0.0, 0.0], 294.15);
        let mut hashes = Vec::new();
        for workers in [1, 3] {
            let mut grid = uniform(0, workers, s);
            grid.set_open_face(0, false, Some(hot)).unwrap();
            grid.set_open_face(0, true, Some(s)).unwrap();
            for _ in 0..8 {
                grid.step().unwrap();
            }
            hashes.push(grid.state_hash());
        }
        assert_eq!(hashes[0], hashes[1]);
    }

    #[test]
    fn g2_prescribed_open_face_reflection_diagnostic() {
        // Fixed physical pulse and domain; only dx changes. Prediction: a
        // prescribed exterior state is not a characteristic-matched outlet.
        // This measures its reflection; the test does NOT assert transparency.
        for n in [64, 128, 256] {
            let mut s = spec(0, 1);
            s.dims = [n, 1, 1];
            s.maximum_cells = n;
            s.dx_m = 0.02 / n as f64;
            let gas = GasState::try_new(&s.gas, s.temperature_k, s.pressure_pa).unwrap();
            let c = (gas.gamma * THETA0).sqrt();
            let amplitude = 1e-6;
            let mut grid = NsfGrid3::<D3q27>::new(
                s,
                |_| false,
                |x| {
                    let position = (x[0] as f64 + 0.5) / n as f64;
                    let a = amplitude * fs_math::det::exp(-0.5 * ((position - 0.3) / 0.05).powi(2));
                    (
                        1.0 + a,
                        [c * a, 0.0, 0.0],
                        293.15 * (1.0 + (gas.gamma - 1.0) * a),
                    )
                },
            )
            .unwrap();
            for high in [false, true] {
                grid.set_open_face(0, high, Some(state([0.0; 3], 293.15)))
                    .unwrap();
            }
            let mut incident: f64 = 0.0;
            let mut reflected: f64 = 0.0;
            for step in 0..(1.5 * n as f64 / c).ceil() as usize {
                grid.step().unwrap();
                let (rho, u, temperature) = grid.cell([n / 2, 0, 0]);
                let pressure = rho * grid.law.theta_of_kelvin(temperature) - THETA0;
                let travel = step as f64 * c / n as f64;
                if travel < 0.5 {
                    incident = incident.max((pressure + c * u[0]).abs());
                }
                if travel > 0.8 {
                    reflected = reflected.max((pressure - c * u[0]).abs());
                }
            }
            assert!(incident > amplitude * 0.5 * c * c);
            assert!(reflected.is_finite());
            println!(
                "NSF_OPEN prescribed_reflection cells={n} dx_m={} amplitude_ratio={:e} transparency_qualified=false",
                grid.law.dx_m,
                reflected / incident
            );
        }
    }
}
