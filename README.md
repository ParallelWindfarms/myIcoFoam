---
title: icoFoam with custom integrator
author: Johan Hidding (Netherlands eScience Center)
---

This shows howto implement a custom `ddt` solver in a OpenFOAM project. I've copied the Euler integrator and disassembled it to a literate (markdown) file. The sources of this Euler integrator are in the `lit` folder.

You need to have [`enTangleD` installed](https://jhidding.github.io/enTangleD) to work with the literate source.

# Building

Assuming you have OpenFOAM installed,

- Activate the OpenFOAM environment by sourcing `${FOAM_INST_DIR}/etc/bashrc`.
- Run `entangled` on `lit/euler.md`, from the project root.
- Compile using `wmake`
  
# Running

The `elbow` example is included here. To run it, first generate the mesh using `fluentMeshToFoam elbow.msh`, then run `myIcoFoam`.
