# El Centro 1940 NS fixture provenance

## Record identity

- Fixture: `el_centro_1940_ns.at2`
- Upstream description: El Centro 1940 North South Component, Peknold version
- Upstream repository: `peer-open-source/opensees-gallery`
- Pinned revision: `fd3958fde8a3d0a350321b3fbdd3f415ee16e2a2`
- Pinned source: <https://github.com/peer-open-source/opensees-gallery/blob/fd3958fde8a3d0a350321b3fbdd3f415ee16e2a2/content/en/examples/Example3/elCentro.AT2>
- Upstream SHA-256: `e096f0458ae7565ac8e9fac80eda018e4715a357ea4c14504ce1cfa6556c53a6`
- Upstream byte count: 16,031
- Repository SHA-256: `e9ae5a8c2163c28a2e52f29f8b363e2f19a3eac400e61ff73c547f8eb309b6ab`
- Repository byte count: 16,032; the sole normalization is one terminal LF
- Header contract: 1,559 uniformly spaced samples, `dt = 0.02 s`, values in `g`
- Peak absolute stored sample: `0.31882 g`
- Test conversion: exactly `9.80665 m/s²` per `g`
- FrankenSim timing convention: source row 0 is designated initial forcing;
  checked step-end integration advances rows 1 through 1,558
- Station/event context: <https://www.strongmotioncenter.org/vdc/scripts/event.plx?evt=88>

This pins the OpenSees Gallery's Peknold byte stream. It is not represented as
the byte-identical CESMD corrected trace; different processing histories produce
different PGA values and must remain separately identified artifacts. The
fixture is redistributed under the upstream repository's BSD-3-Clause notice
below; no separate file-specific public-domain declaration is claimed.

## Public comparison artifact

The pinned OpenSees Example 3 model is a one-story, one-bay reinforced-concrete
frame with force-based fiber-section columns and the same Peknold motion:

- Model: <https://github.com/peer-open-source/opensees-gallery/blob/fd3958fde8a3d0a350321b3fbdd3f415ee16e2a2/content/en/examples/Example3/portal.tcl>
- Model SHA-256: `b5b6e08cec39d6b1826a9d27778c3aa1f76adfb115d7ba66fb6eaf1bc7c6dd02`
- Committed output: <https://github.com/peer-open-source/opensees-gallery/blob/fd3958fde8a3d0a350321b3fbdd3f415ee16e2a2/content/en/examples/Example3/out/disp.out>
- Output SHA-256: `b1f6f933b71e543e85deea342aca54d3389e0435fd6bccdf9e0da595490efdef`
- Committed peak roof displacement: `2.47158 in = 0.062778132 m`
- Model story height: `144 in = 3.6576 m`
- Committed peak drift ratio: `0.01716375`

The model records envelope roof displacement and one column's section force. It
does not record total dynamic support reaction. A section force is not total
base reaction, and independently peaking nodal accelerations cannot be summed
into one. Therefore a public total-base-reaction oracle is unavailable from this artifact.

## FrankenSim comparison authority

`fs-frame::StoryFrame` is materially non-equivalent to Example 3: it is a
two-column concentrated-base-hinge story with a 3.0 m default height, 280,000 kg
mass, hard-coded 0.5 m by 0.35 m section fixture, no elastic girder, no gravity
preload, and no P-Delta transform. Example 3 uses 3.6576 m force-based distributed
columns, a different fiber section and mass, an elastic girder, gravity preload,
and P-Delta kinematics.

Accordingly, the current battery may prove record integrity and produce an
auditable FrankenSim displacement/restoring-shear response, but it must not gate
that response against the OpenSees displacement as though the models were the
same. The external displacement is a pinned diagnostic reference. A numerical
acceptance band requires an admitted model mapping; a public total-base-
reaction band additionally requires an authoritative reaction artifact.

## Upstream license notice

BSD 3-Clause License

Copyright (c) 2024, Structural Artificial Intelligence Research Lab

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.
3. Neither the name of the copyright holder nor the names of its contributors
   may be used to endorse or promote products derived from this software without
   specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
