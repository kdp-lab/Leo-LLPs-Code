# LCTuple configurations

Configurations for the [`LCTuple`](https://github.com/MuonColliderSoft/LCTuple) Marlin processor that produces flat `ROOT` trees and takes `LCIO` files as input.  
Several configurations are available for different stages:
* `mcp.xml` - before GEANT4 simulation stage, takes only `MCParticle` objects as input;
* `sim.xml` - after GEANT4 simulation stage, takes only `SimHit` objects as input;
* `digi.xml` - after digitisation stage, takes only `RecHit` objects as input;
