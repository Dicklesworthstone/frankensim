# fs-kamen-wasm browser contract

`kamen_cluster_step(state_index)` maps the six public teaching states for US
5,701,965 onto the nominal geometry printed in Table 1 and the discrete poses
shown in Figures 39–42, then delegates rigid wheel-axis transforms and
horizontal-support gap evaluation to
`fs-mbd::tri_wheel_cluster::step_tri_wheel_stair_contact`.

State indices are ordered as follows:

0. ground support
1. fore-aft balance
2. stair start
3. weight transfer
4. climb
5. transition gate

The success envelope reports SI geometry, pose, three wheel centers, signed
vertical gaps, contact flags, and the generic owner/boundary strings. The
boundary is deliberately narrow: three rigid equal wheels, one planar carrier,
horizontal level/tread contacts, ideal sharp stair corners, and no force,
friction, tire compliance, impact, motor, controller, sensor, or riser-side
contact calculation. Unsupported or penetrating poses are typed refusals.
