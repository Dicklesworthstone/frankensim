//! G1/G4: original stress operator, physical vorticity production and reuse.

use super::*;
use fs_gas::multicomponent::dissipation::DissipationError;

#[test]
fn g1_viscous_vortex_rate_matches_discrete_decay_and_refines() {
    let lengths = [1.1, 0.7];
    let (amplitude, density, mu) = (0.08, 1.3, 0.025);
    let mut errors = Vec::new();
    for n in [8, 16, 32] {
        let h = lengths.map(|l| l / n as f64);
        let k = lengths.map(|l| std::f64::consts::TAU / l);
        let lambda = h.map(|dx| 4.0 * (std::f64::consts::PI / n as f64).sin().powi(2) / dx.powi(2));
        let derivative = [0, 1].map(|a| (k[a] * h[a]).sin() / h[a]);
        let grid = CartesianGrid::new(n, n, 1, h[0], h[1], 0.4).unwrap();
        let mut flow = field(gas(0.0), grid, PERIODIC, |x| {
            [
                density,
                -amplitude * (k[1] * x[1]).sin(),
                amplitude * (k[0] * x[0]).sin(),
                0.0,
                2.0,
            ]
        });
        let mut op =
            MolecularTransport3d::new(&flow, law(flow.gas(), mu, 0.017, 0.0), || false).unwrap();
        let mut work = op
            .vorticity_workspace(op.required_vorticity_workspace_bytes().unwrap())
            .unwrap();
        let values = op.vorticity(&flow, &mut work, || false).unwrap();
        let mut error = 0.0;
        let mut norm = 0.0;
        for (i, value) in values.iter().enumerate() {
            let x = grid.cell_center(i).unwrap();
            let cos = [0, 1].map(|a| (k[a] * x[a]).cos());
            let exact_discrete = -mu / density
                * amplitude
                * (lambda[0] * derivative[0] * cos[0] + lambda[1] * derivative[1] * cos[1]);
            let continuum =
                -mu / density * amplitude * (k[0].powi(3) * cos[0] + k[1].powi(3) * cos[1]);
            close(value.stress[0], 0.0, 3e-12);
            close(value.stress[1], 0.0, 3e-12);
            close(value.stress[2], exact_discrete, 3e-11);
            error += (value.stress[2] - continuum).powi(2);
            norm += continuum.powi(2);
        }
        errors.push((error / norm).sqrt());

        // The actual native molecular SSP-RK2 update must damp these same
        // vorticity modes, including anisotropic physical wave numbers.
        let dt = 1e-5;
        let report = op
            .advance(
                &mut flow,
                Time::new(dt),
                MolecularControl::default(),
                || false,
            )
            .unwrap();
        assert_eq!(report.accepted_steps, 1);
        let values = op.vorticity(&flow, &mut work, || false).unwrap();
        let decay = lambda.map(|v| {
            let z = mu / density * v * dt;
            1.0 - z + 0.5 * z * z
        });
        for (i, value) in values.iter().enumerate() {
            let x = grid.cell_center(i).unwrap();
            let expected = amplitude
                * (decay[0] * derivative[0] * (k[0] * x[0]).cos()
                    + decay[1] * derivative[1] * (k[1] * x[1]).cos());
            close(value.inviscid.vorticity[2], expected, 3e-12);
        }
    }
    for pair in errors.windows(2) {
        assert!(
            pair[1] < 0.27 * pair[0],
            "vortex curl refinement: {errors:?}"
        );
    }
}

struct TemperatureViscosity(ConstantTransport);

impl TransportLaw for TemperatureViscosity {
    fn source(&self) -> &IdealGasEquilibriumV1 {
        self.0.source()
    }

    fn evaluate(
        &mut self,
        temperature: Temperature,
        pressure: Pressure,
        fractions: &[f64],
        binary: &mut [f64],
        cancelled: &mut dyn FnMut() -> bool,
    ) -> Result<TransportCoefficients, DissipationError> {
        let mut value = self
            .0
            .evaluate(temperature, pressure, fractions, binary, cancelled)?;
        value.dynamic_viscosity =
            DynViscosity::new(value.dynamic_viscosity.value() * temperature.value());
        Ok(value)
    }
}

#[test]
fn g1_original_variable_viscosity_and_density_generate_irrotational_curl() {
    let open = AxisBoundaries::Ends {
        lower: BoundaryCondition::Transmissive,
        upper: BoundaryCondition::Transmissive,
    };
    let boundaries = Boundaries3d {
        x: open,
        y: open,
        z: AxisBoundaries::Periodic,
    };
    let (strain, rho0, density_slope, quadratic_mu) = (0.15, 0.8, 0.4, 0.01);
    let mut errors = Vec::new();
    for n in [8, 16, 32] {
        let h = 1.0 / n as f64;
        let grid = CartesianGrid::new(n, 5, 1, h, 0.19, 0.4).unwrap();
        // u=(b*y,b*x,0) has identically zero curl and divergence.
        // mu=0.02+0.01*x^2, rho=rho0+c*x, T=2+x^2.
        // Both pressure and density depend only on x: baroclinicity is zero.
        let flow = field(gas(0.4), grid, boundaries, |x| {
            let rho = rho0 + density_slope * x[0];
            [
                rho,
                strain * x[1],
                strain * x[0],
                0.0,
                rho * (2.0 + x[0] * x[0]),
            ]
        });
        let original = TemperatureViscosity(law(flow.gas(), quadratic_mu, 0.007, 0.0));
        let mut op = MolecularTransport3d::new(&flow, original, || false).unwrap();
        let mut work = op.vorticity_workspace(usize::MAX).unwrap();
        let cell = grid.index([n / 2, 2, 0]);
        let value = op.vorticity(&flow, &mut work, || false).unwrap()[cell];
        let x = grid.cell_center(cell).unwrap()[0];
        let rho = rho0 + density_slope * x;
        // div(tau)=(0,4*b*k*x,0). This independent rational derivative
        // distinguishes curl(div(tau)/rho) from nu*laplacian(omega)=0.
        let numerator = 4.0 * strain * quadratic_mu * rho0;
        let continuum = numerator / (rho * rho);
        let exact_discrete = numerator / (rho * rho - (density_slope * h).powi(2));
        for v in value
            .inviscid
            .vorticity
            .into_iter()
            .chain(value.inviscid.baroclinic)
        {
            close(v, 0.0, 2e-11);
        }
        close(value.stress[0], 0.0, 2e-11);
        close(value.stress[1], 0.0, 2e-11);
        close(value.stress[2], exact_discrete, 2e-10);
        assert!(value.stress[2] > 0.001);
        close(
            value.resolved_material_rate()[2],
            value.inviscid.inviscid_material_rate()[2] + exact_discrete,
            2e-10,
        );
        errors.push((value.stress[2] - continuum).abs());
    }
    for pair in errors.windows(2) {
        assert!(
            pair[1] < 0.29 * pair[0],
            "variable-stress refinement: {errors:?}"
        );
    }
}

#[test]
fn g1_stationary_no_slip_and_navier_walls_create_their_actual_vorticity() {
    let (h, mu, velocity) = (0.2_f64, 0.03, 0.2);
    let grid = CartesianGrid::new(1, 5, 1, 0.4, h, 0.6).unwrap();
    let flow = field(
        gas(0.0),
        grid,
        Boundaries3d {
            y: WALL,
            ..PERIODIC
        },
        |_| [1.0, velocity, 0.0, 0.0, 2.0],
    );
    let mut slip_default =
        MolecularTransport3d::new(&flow, law(flow.gas(), mu, 0.0, 0.0), || false).unwrap();
    let mut work = slip_default.vorticity_workspace(usize::MAX).unwrap();
    for value in slip_default.vorticity(&flow, &mut work, || false).unwrap() {
        assert_eq!(value.stress, [0.0; 3]);
    }
    for slip in [0.0, 0.13] {
        let wall = if slip == 0.0 {
            MolecularWall3d::NoSlip {
                velocity: [0.0; 3],
                thermal: WallThermalCondition::Adiabatic,
            }
        } else {
            MolecularWall3d::NavierSlip {
                velocity: [0.0; 3],
                thermal: WallThermalCondition::Adiabatic,
                slip_length_m: slip,
                fluid_heat_fraction: 1.0,
            }
        };
        let mut op = MolecularTransport3d::new(&flow, law(flow.gas(), mu, 0.0, 0.0), || false)
            .unwrap()
            .with_wall(BoundarySide::YMin, wall)
            .unwrap()
            .with_wall(BoundarySide::YMax, wall)
            .unwrap();
        let factor = mu * velocity / (h * h * (0.5 * h + slip));
        for (value, coefficient) in op
            .vorticity(&flow, &mut work, || false)
            .unwrap()
            .iter()
            .zip([-1.5, -0.5, 0.0, 0.5, 1.5])
        {
            close(value.stress[0], 0.0, 2e-13);
            close(value.stress[1], 0.0, 2e-13);
            close(value.stress[2], coefficient * factor, 2e-13);
        }
    }
}

#[test]
fn g4_coupled_query_caps_cancellation_and_reuse_preserve_accepted_state() {
    let flow = cube(4, |x| {
        let s = (std::f64::consts::TAU * x[0]).sin();
        [1.0, 0.0, 0.1 * s, 0.0, 2.0]
    });
    let source = law(flow.gas(), 0.02, 0.003, 0.01);
    let mut owner = MolecularEuler3d::new(flow, source, || false).unwrap();
    owner
        .step(Time::new(1e-5), MolecularStepControl::default(), || false)
        .unwrap();
    let bytes = owner.required_vorticity_workspace_bytes().unwrap();
    assert!(owner.vorticity_workspace(bytes - 1).is_err());
    let mut work = owner.vorticity_workspace(bytes).unwrap();
    let capacity = work.workspace_bytes();
    let states = owner.flow().states().to_vec();
    let time = owner.flow().time();
    let report = *owner.last_step().unwrap();
    let mut polls = 0;
    let expected = owner
        .vorticity(&mut work, || {
            polls += 1;
            false
        })
        .unwrap()
        .to_vec();
    for stop in [1, polls / 2, polls] {
        let mut seen = 0;
        assert_eq!(
            owner.vorticity(&mut work, || {
                seen += 1;
                seen == stop
            }),
            Err(MolecularError::Cancelled)
        );
        assert_eq!(owner.flow().states(), states);
        assert_eq!(owner.flow().time(), time);
        assert_eq!(owner.last_step(), Some(&report));
        assert_eq!(owner.vorticity(&mut work, || false).unwrap(), expected);
        assert_eq!(work.workspace_bytes(), capacity);
    }
}
