//! Admit and retain the native energy and momentum transport choices.

use fs_lbm::conjugate::{ConvectionScheme, EnergyConfig, MomentumReconstruction, SimpleConfig};
use std::fmt::Write as _;

use super::{J, Result, bad, quote};

pub(super) const HELP: &str = "\nConvection: solver.energy_convection selects \"power-law\", \"upwind\",\n\"exponential\" or \"enthalpy-minmod\". Source fluid_enthalpy requires\n\"upwind\" (its default) or \"enthalpy-minmod\"; other energy paths default\nto \"power-law\" and also accept \"enthalpy-minmod\" under their declared\nconstant Cp, with a converged native temperature reconstruction. solver.momentum_convection\nselects \"power-law\" (default), \"upwind\", \"exponential\" or \"minmod\".\nMomentum minmod uses native limited face velocities with an implicit\nupwind operator. Exponential fitting resolves constant-coefficient 1-D\nconvection-diffusion between neighboring cells; it is not a general\nhigher-order multidimensional method. Limited reconstruction retains\nfirst-order treatment where the native boundary stencil is incomplete.\nSpatial choices do not change the selected time integrator or its order.\n";

#[derive(Clone, Copy)]
pub(super) struct Selection {
    energy: ConvectionScheme,
    momentum: ConvectionScheme,
    momentum_reconstruction: MomentumReconstruction,
}

impl Selection {
    pub(super) fn parse(solver: Option<&J>, source_enthalpy: bool) -> Result<Self> {
        if solver.is_some_and(|solver| solver.as_object().is_none()) {
            return Err(bad("solver must be an object"));
        }
        let choice = |key| -> Result<Option<&str>> {
            solver
                .and_then(|solver| solver.get(key))
                .map(|value| {
                    value
                        .as_str()
                        .ok_or_else(|| bad(format!("solver.{key} must be a string")))
                })
                .transpose()
        };
        let energy = match choice("energy_convection")? {
            None if source_enthalpy => ConvectionScheme::Upwind,
            None | Some("power-law") => ConvectionScheme::PowerLaw,
            Some("upwind") => ConvectionScheme::Upwind,
            Some("exponential") => ConvectionScheme::Exponential,
            Some("enthalpy-minmod") => ConvectionScheme::EnthalpyMinmod,
            Some(_) => {
                return Err(bad(
                    "solver.energy_convection must be power-law, upwind, exponential or enthalpy-minmod",
                ));
            }
        };
        if source_enthalpy
            && !matches!(
                energy,
                ConvectionScheme::Upwind | ConvectionScheme::EnthalpyMinmod
            )
        {
            return Err(bad(
                "fluid_enthalpy requires solver.energy_convection upwind or enthalpy-minmod",
            ));
        }
        // Without a source caloric field, the native owner converges the same
        // limited transport under the declared reference Cp.
        let (momentum, momentum_reconstruction) = match choice("momentum_convection")? {
            None | Some("power-law") => (
                ConvectionScheme::PowerLaw,
                MomentumReconstruction::DonorCell,
            ),
            Some("upwind") => (ConvectionScheme::Upwind, MomentumReconstruction::DonorCell),
            Some("exponential") => (
                ConvectionScheme::Exponential,
                MomentumReconstruction::DonorCell,
            ),
            Some("minmod") => (ConvectionScheme::Upwind, MomentumReconstruction::Minmod),
            Some(_) => {
                return Err(bad(
                    "solver.momentum_convection must be power-law, upwind, exponential or minmod; enthalpy-minmod is an energy scheme",
                ));
            }
        };
        Ok(Self {
            energy,
            momentum,
            momentum_reconstruction,
        })
    }

    pub(super) fn energy_config(self) -> EnergyConfig {
        EnergyConfig {
            scheme: self.energy,
            ..EnergyConfig::default()
        }
    }

    pub(super) fn configure_momentum(self, config: &mut SimpleConfig) {
        config.scheme = self.momentum;
        config.momentum_reconstruction = self.momentum_reconstruction;
    }

    pub(super) fn energy_name(self) -> &'static str {
        name(self.energy)
    }

    fn momentum_name(self) -> &'static str {
        if self.momentum_reconstruction == MomentumReconstruction::Minmod {
            "minmod"
        } else {
            name(self.momentum)
        }
    }

    pub(super) fn append_json(self, out: &mut String) {
        let _ = write!(
            out,
            ",\"convection\":{{\"energy\":{},\"momentum\":{},\"energy_scope\":{},\"momentum_scope\":{}}}",
            quote(self.energy_name()),
            quote(self.momentum_name()),
            quote(scope(self.energy)),
            quote(
                if self.momentum_reconstruction == MomentumReconstruction::Minmod {
                    "conservative limited face-velocity reconstruction; incomplete boundary stencils retain donor transport; no general kinetic-energy stability or multidimensional maximum-principle claim"
                } else {
                    scope(self.momentum)
                }
            ),
        );
    }

    pub(super) fn append_text(self, out: &mut String) {
        let _ = writeln!(
            out,
            "energy_convection={}\nmomentum_convection={}",
            self.energy_name(),
            self.momentum_name()
        );
    }
}

fn name(scheme: ConvectionScheme) -> &'static str {
    match scheme {
        ConvectionScheme::PowerLaw => "power-law",
        ConvectionScheme::Upwind => "upwind",
        ConvectionScheme::Exponential => "exponential",
        ConvectionScheme::EnthalpyMinmod => "enthalpy-minmod",
    }
}

fn scope(scheme: ConvectionScheme) -> &'static str {
    match scheme {
        ConvectionScheme::PowerLaw => {
            "power-law fitting; first-order donor transport at high cell Peclet number"
        }
        ConvectionScheme::Upwind => "first-order donor transport with molecular diffusion",
        ConvectionScheme::Exponential => {
            "constant-coefficient 1-D exponential fitting; no general higher-order multidimensional claim"
        }
        ConvectionScheme::EnthalpyMinmod => {
            "conservative limited enthalpy reconstruction; native temperature paths use declared constant Cp and converged face corrections; incomplete boundary stencils retain donor transport; no general multidimensional maximum-principle claim"
        }
    }
}
