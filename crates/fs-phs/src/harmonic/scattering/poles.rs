//! Natural modes with outgoing matched boundaries and original physical loss.
use super::{C64, HarmonicError, ScatteringModel, finite, poll};
use crate::harmonic::{
    BeynOptions, HarmonicPole, HarmonicPoleError, HarmonicSpectrumError, PoleOptions,
    ResonanceRegion,
};
use fs_la::nonlinear_eig::NonlinearEigError;

pub mod response;

/// One outgoing boundary account for a mode proportional to exp(s*t).
/// Quadratic values are complex-envelope powers at the returned arbitrary
/// amplitude, without the factor 1/2 used for peak forced harmonic averages.
#[derive(Clone, Debug, PartialEq)]
pub struct DecayWavePort {
    pub original_port: usize,
    /// Positive output/input reference, with the same units as in wave response.
    pub reference: f64,
    /// b=y: there is no incident wave in a free outgoing mode.
    pub outgoing: C64,
    /// Actual original input u=-y/Z.
    pub physical_input: C64,
    /// |b|²/Z, power leaving through this boundary; not internal material heat.
    pub outgoing_power_w: f64,
    /// Re(conj(u)*y), power into the original physical model, evaluated separately.
    pub physical_power_w: f64,
    /// Outgoing + physical power. No residual is folded into a loss channel.
    pub power_residual_w: f64,
}

/// A mode of the fully closed original descriptor pencil, with separate heat
/// and outgoing-wave accounts. A finite mode need not be visible at a port.
#[derive(Clone, Debug, PartialEq)]
pub struct ScatteringPole {
    /// Laplace s [1/s], with Re(s)<0 for decay. Unstable values are not clipped.
    pub laplace_pole_per_s: C64,
    /// omega=i*s under exp(-i*omega*t).
    pub angular_frequency_rad_s: C64,
    /// Every original differential and algebraic coordinate, unit numerical norm.
    pub state: Vec<C64>,
    pub effort: Vec<C64>,
    pub port_output: Vec<C64>,
    /// Actual inputs at all original ports; unselected ordinary inputs are zero.
    pub physical_inputs: Vec<C64>,
    pub waves: Vec<DecayWavePort>,
    /// x^H Q x/2. Algebraic constraints and boundary references store no energy.
    pub stored_energy_j: f64,
    /// Original Re(e^H R e + e^H T x). For equilibrium linearizations, T gives
    /// signed incremental work, not the total heat of the biased physical system.
    pub internal_dissipated_power_w: f64,
    /// Sum of independently computed |y|²/Z at the selected boundaries.
    pub outgoing_power_w: f64,
    /// Internal dissipative work plus outward boundary exchange.
    pub total_decay_power_w: f64,
    /// 2*Re(s)*H + internal loss + outgoing power.
    pub decay_power_residual_w: f64,
    /// 2*Re(s)*H + internal loss - sum Re(conj(u)*y), independently evaluated.
    pub physical_power_residual_w: f64,
    /// Maximum of closed-pencil and original actual-input equation residuals.
    pub relative_residual: f64,
    pub iterations: usize,
}

/// Numerically resolved finite outgoing modes in an authored Laplace region.
/// Shared contour count/rank/quadrature evidence is not an interval proof of
/// completeness; all source bands and spatial model limits remain in force.
#[derive(Clone, Debug)]
pub struct ScatteringSpectrum {
    pub poles: Vec<ScatteringPole>,
    pub unrefined_laplace_poles_per_s: Vec<C64>,
    pub argument_count: usize,
    pub resolved_rank: usize,
    pub quadrature_nodes: usize,
    pub argument_integral: C64,
    pub coarse_argument_integral: C64,
    pub moment_relative_change: f64,
    pub minimum_reciprocal_condition: f64,
}

fn pole_error(error: HarmonicError) -> HarmonicPoleError {
    match error {
        HarmonicError::Cancelled => HarmonicPoleError::Solve(NonlinearEigError::Cancelled),
        HarmonicError::Invalid { what } => HarmonicPoleError::Invalid { what },
        HarmonicError::Residual {
            relative,
            tolerance,
        } => HarmonicPoleError::Residual {
            relative,
            tolerance,
        },
        HarmonicError::Singular { .. } => HarmonicPoleError::Invalid {
            what: "original outgoing-mode observation cannot factor an operator",
        },
    }
}

impl ScatteringModel {
    // Additional wave observation storage and work, on top of the existing
    // physical/numerical owner's resources. Original mode vectors move without
    // copying. Arithmetic uses u128 before admission and narrowing.
    fn pole_observation_resources(&self, modes: usize) -> (u128, u128) {
        let (n, m, k, modes) = (
            self.closed.n as u128,
            self.closed.m as u128,
            self.ports.len() as u128,
            modes as u128,
        );
        (
            modes * (16 * m + 160 * k + 256),
            modes * (32 * n * n + 16 * n * m + 32 * k + 8 * m),
        )
    }

    /// Refine one finite natural mode with zero incident waves at every selected
    /// port, retaining every original material, load and constraint coordinate.
    /// Ordinary source inputs are homogeneous. Closure precedes the solve, so
    /// no bare component impedance or singular descriptor mass is inverted.
    ///
    /// The shared original-state Newton owner solves the complete pencil. This
    /// observer restores actual physical boundary inputs and checks the original
    /// unclosed equations independently. Boundary decay is not added to R or
    /// reported as material heat. No accepted time history or time step is used.
    /// The supplied options include additional outgoing-wave observation costs.
    pub fn resonance_pole(
        &self,
        laplace_seed_per_s: C64,
        state_seed: &[C64],
        options: &PoleOptions,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<ScatteringPole, HarmonicPoleError> {
        poll(&mut cancelled).map_err(pole_error)?;
        let (bytes, work) = self.pole_observation_resources(1);
        if bytes > options.max_dense_bytes as u128 || work > options.max_work as u128 {
            return Err(HarmonicPoleError::Invalid {
                what: "insufficient outgoing-mode observation resources",
            });
        }
        let mut numerical = *options;
        numerical.max_dense_bytes -= bytes as usize;
        numerical.max_work -= work as u64;
        let mode = self.closed.resonance_pole(
            laplace_seed_per_s,
            state_seed,
            &numerical,
            &mut cancelled,
        )?;
        self.observe_outgoing_pole(mode, &mut cancelled)
    }

    /// Discover finite natural modes without individual seeds. The shared Beyn
    /// owner retains its rank, count, conditioning and nested-quadrature gates;
    /// every accepted mode additionally passes the original physical rows with
    /// actual outgoing-boundary inputs. Dark internal modes are retained too.
    /// Any failed observation refuses the entire result, leaving AC factors and
    /// accepted runtime history unchanged. Options also bound all wave accounts.
    pub fn resonance_spectrum(
        &self,
        region: ResonanceRegion,
        options: &BeynOptions,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<ScatteringSpectrum, HarmonicSpectrumError> {
        poll(&mut cancelled).map_err(|e| HarmonicSpectrumError::PhysicalPole(pole_error(e)))?;
        let (bytes, work) = self.pole_observation_resources(options.probe_columns);
        if bytes > options.max_dense_bytes as u128 || work > options.max_work as u128 {
            return Err(HarmonicSpectrumError::Invalid {
                what: "insufficient outgoing-spectrum observation resources",
            });
        }
        let mut numerical = *options;
        numerical.max_dense_bytes -= bytes as usize;
        numerical.max_work -= work as u64;
        let spectrum = self
            .closed
            .resonance_spectrum(region, &numerical, &mut cancelled)?;
        let mut poles = Vec::with_capacity(spectrum.poles.len());
        for mode in spectrum.poles {
            poles.push(
                self.observe_outgoing_pole(mode, &mut cancelled)
                    .map_err(HarmonicSpectrumError::PhysicalPole)?,
            );
        }
        poll(&mut cancelled).map_err(|e| HarmonicSpectrumError::PhysicalPole(pole_error(e)))?;
        Ok(ScatteringSpectrum {
            poles,
            unrefined_laplace_poles_per_s: spectrum.unrefined_laplace_poles_per_s,
            argument_count: spectrum.argument_count,
            resolved_rank: spectrum.resolved_rank,
            quadrature_nodes: spectrum.quadrature_nodes,
            argument_integral: spectrum.argument_integral,
            coarse_argument_integral: spectrum.coarse_argument_integral,
            moment_relative_change: spectrum.moment_relative_change,
            minimum_reciprocal_condition: spectrum.minimum_reciprocal_condition,
        })
    }

    fn observe_outgoing_pole(
        &self,
        mode: HarmonicPole,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<ScatteringPole, HarmonicPoleError> {
        let mut physical_inputs = vec![C64::ZERO; self.closed.m];
        let mut waves = Vec::with_capacity(self.ports.len());
        let (mut outgoing_power_w, mut physical_power_w) = (0.0, 0.0);
        for &(p, z) in &self.ports {
            poll(cancelled).map_err(pole_error)?;
            let output = mode.port_output[p];
            let input = output.scale(-1.0 / z);
            physical_inputs[p] = input;
            let outgoing = output.norm_sq() / z;
            let physical = (input.conj() * output).re;
            let residual = outgoing + physical;
            if !finite(input)
                || [outgoing, physical, residual]
                    .iter()
                    .any(|x| !x.is_finite())
            {
                return Err(HarmonicPoleError::Invalid {
                    what: "outgoing-mode boundary power is not representable",
                });
            }
            outgoing_power_w += outgoing;
            physical_power_w += physical;
            waves.push(DecayWavePort {
                original_port: p,
                reference: z,
                outgoing: output,
                physical_input: input,
                outgoing_power_w: outgoing,
                physical_power_w: physical,
                power_residual_w: residual,
            });
        }
        let total_decay_power_w = mode.dissipated_power_w + outgoing_power_w;
        let energy_rate = 2.0 * mode.laplace_pole_per_s.re * mode.stored_energy_j;
        let decay_power_residual_w = energy_rate + total_decay_power_w;
        let physical_power_residual_w = energy_rate + mode.dissipated_power_w - physical_power_w;
        if [
            outgoing_power_w,
            total_decay_power_w,
            decay_power_residual_w,
            physical_power_residual_w,
        ]
        .iter()
        .any(|x| !x.is_finite())
        {
            return Err(HarmonicPoleError::Invalid {
                what: "outgoing-mode decay balance is not representable",
            });
        }
        let relative_residual = mode.relative_residual.max(
            self.original_residual(
                mode.laplace_pole_per_s,
                &mode.state,
                &mode.effort,
                &physical_inputs,
                cancelled,
            )
            .map_err(pole_error)?,
        );
        poll(cancelled).map_err(pole_error)?;
        Ok(ScatteringPole {
            laplace_pole_per_s: mode.laplace_pole_per_s,
            angular_frequency_rad_s: mode.angular_frequency_rad_s,
            state: mode.state,
            effort: mode.effort,
            port_output: mode.port_output,
            physical_inputs,
            waves,
            stored_energy_j: mode.stored_energy_j,
            internal_dissipated_power_w: mode.dissipated_power_w,
            outgoing_power_w,
            total_decay_power_w,
            decay_power_residual_w,
            physical_power_residual_w,
            relative_residual,
            iterations: mode.iterations,
        })
    }
}
