//! Authored exact DG cylinder; geometry and Bessel memory lower once into the
//! original fs-duct owner, then enter the existing physical aperture runtime.
use super::*;
use fs_duct::dg::{
    DgBoundary, DgDuct, DgMesh, MAX_DG_ELEMENTS, MAX_DG_ORDER, MAX_DG_PROPAGATOR_DOFS, PortId,
};
use fs_duct::{Duct, Segment};

pub(super) struct DgInput {
    length_m: f64,
    radius_m: f64,
    outlet: DgBoundary,
    mesh: DgMesh,
    poles: usize,
    band_hz: [f64; 2],
    maximum_relative_error: f64,
    exponential_tolerance: f64,
    maximum_dofs: usize,
}

impl DgInput {
    pub(super) fn read(row: &mut Row<'_>) -> Result<Self, PlateValveInputError> {
        let length_m = row.scalar()?;
        let radius_m = row.scalar()?;
        let outlet = match row.word()? {
            "open" => DgBoundary::IdealOpen,
            "closed" => DgBoundary::Closed,
            "matched" => DgBoundary::Matched,
            _ => {
                return Err(bad(
                    row.line,
                    "dg_tube outlet must be open, closed or matched; no radiation model is implied",
                ));
            }
        };
        let mesh = DgMesh {
            elements: row.count(MAX_DG_ELEMENTS)?,
            order: row.count(MAX_DG_ORDER)?,
        };
        let poles = row.count(fs_phs::MAX_ZK_DIFFUSIVE_POLES)?;
        let band_hz = [row.scalar()?, row.scalar()?];
        let maximum_relative_error = row.scalar()?;
        let exponential_tolerance = row.scalar()?;
        let maximum_dofs = row.count(MAX_DG_PROPAGATOR_DOFS)?;
        if length_m <= 0.0
            || radius_m <= 0.0
            || mesh.elements == 0
            || mesh.order == 0
            || poles == 0
            || exponential_tolerance <= 0.0
            || exponential_tolerance >= 1.0
        {
            return Err(bad(
                row.line,
                "dg_tube needs positive geometry, mesh, poles and an exponential tolerance below one",
            ));
        }
        Ok(Self {
            length_m,
            radius_m,
            outlet,
            mesh,
            poles,
            band_hz,
            maximum_relative_error,
            exponential_tolerance,
            maximum_dofs,
        })
    }

    pub(super) fn prepare(
        self,
        gas: &GasState,
        dt: f64,
        gate: &CancelGate,
    ) -> Result<ApertureDgLoad, PlateValveInputError> {
        checkpoint(gate)?;
        // Preflight the actual state count before the original owner allocates
        // dense matrices. gamma=1 omits the exact-zero thermal operator.
        let branches = self.poles + if gas.gamma == 1.0 { 0 } else { self.poles };
        let dofs = self.mesh.elements * (self.mesh.order + 1) * (2 + branches);
        if dofs > self.maximum_dofs {
            return Err(bad(
                0,
                "dg_tube acoustic and Bessel memory states exceed the declared DOF budget",
            ));
        }
        let geometry = Duct {
            segments: vec![Segment::Cylinder {
                radius: self.radius_m,
                length: self.length_m,
            }],
        };
        let owner = DgDuct::new(
            &geometry,
            gas,
            self.mesh,
            DgBoundary::Port(PortId(0)),
            self.outlet,
        )
        .and_then(|duct| {
            duct.with_viscothermal(self.poles, self.band_hz, self.maximum_relative_error)
        })
        .map_err(|error| PlateValveInputError::Physics(AcousticRealizeError::Duct(error)))?;
        checkpoint(gate)?;
        let load =
            ApertureDgLoad::new(owner, dt, self.exponential_tolerance, self.maximum_dofs)
                .map_err(PlateValveInputError::Physics)?;
        checkpoint(gate)?;
        Ok(load)
    }
}
