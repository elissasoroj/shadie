### Goal of the project: Is it clear to you from the proposal.md how the goal can be accomplished using Python and the specified packages?

Yes, I understand it. The one thing that is not totally clear to me is which of the functions of SLiM3, which seems extremely complex, are being preserved, a decision which I assume was made based on what is important for understanding the biphasic lifecycle of plants, which is the type of modeling that shadie is supposed to help with. I understand why it would be beneficial to have a python wrapper in terms of improving ease of use and simplifying so that only the parameters that need to be changed are accessible, but I'm not sure from the proposal who the target audience is, and, as a result, which are the parameters that are particularly important to them. 

### The Data: Is it clear to you from the proposal.md what the data for this project is, or will look like?
To some extent. I know the input is a phylogeny, along with some other editable parameters- I am not totally sure if the user needs to provide a toytree phylogeny, or if shadie will use toytree to build their phylogeny based on other inputs.  

### Does the current code include a proper skeleton (pseudocode) for starting this project?

The code seems very organized and most of the steps seem to be in place, though a lot of it is a little over my head! Next steps would be using the inputs from the toytree that have been read by the code to write a file that SLiM3 can take as input. 

### Code contributions/ideas: 
In the jupyter notebook I added, I showed a proposal for how to make a simple graphical representation of a gene using toyplot, using a single bar stacked bar graph aligned along the x-axis. I also created a basic function to generate some random genes, which can be replaced with actual data once you have it. I'm not sure what format SLiM3 would provide this information in, but the demo I made takes a list as input. 