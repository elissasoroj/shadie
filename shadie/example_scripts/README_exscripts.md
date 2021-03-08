# READ ME: Example Scripts Folder

This folder contains example SLiM scripts for different scenarios. These are the kinds of scripts that `shadie` needs to be able to produce

Note: SLiM3 does not use a standard computer language. It uses `Eidos`, a custom computer language that was written just for SLiM. Text editors like `sublime` or `vscode` will not recognize `Eidos`, but the code can be visualized pretty well using `Java` syntax highlighting.

SLiM3 scripts are organized into two parts, the `initialize()` callbacks, which define a number of parameters for the simulation. Some of these parameters are required and many are optional. Once the callbacks are all defined, the simulation begins with the `Eidos events`

```Java
// Simple simuulation of neutral evolution in panmictic population

//initialize() callbacks - required callbacks only
initialize() {
	initializeMutationRate(1e-7);
	initializeMutationType("m1", 0.5, "f", 0.0);
	initializeGenomicElementType("g1", m1, 1.0);
	initializeGenomicElement(g1, 0, 99999);
	initializeRecombinationRate(1e-8);
}

//Eidos events
1 { sim.addSubpop("p1", 500); }
10000 late() { sim.outputFull(); }
```
