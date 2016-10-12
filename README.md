# Chlamy-Cheminformatics
Implementation of the ECFP algorithm for use in a chemical genomics project in the Jonikas Lab, Stanford Dept. of Plant Biology.

## Background
### ECFP
The ECFP algorithm was first published by David Rogers and Mathew Hahn as a means of generating numerical fingerprints from the graphical structures of molecules. Since the numerical features of a fingerprint are directly encoded by structural elements of the molecule that fingerprint comes from, cheminformatic algorithms can use them to approximate the structural similarity of two molecules without wasting time on computing subgraph isomorphisms (who _wouldn't_ want to avoid solving millions of NP-complete problems?) A description of the algorithm can be found at http://pubs.acs.org/doi/abs/10.1021/ci100050t.

### Why Fingeprinting?
When I worked at the Jonikas lab during the summer of 2016, I was helping out with a chemical genomics project in the green alga _Chlamydomonas reinhardtii_ (aka Chlamy). Chlamy has several metabolic pathways that are of great interest to plant biologists--including [a particularly efficient CO2-concentrating mechanism](https://en.wikipedia.org/wiki/Pyrenoid) that could dramatically increase crop yield if it were engineered into higher plants--and the lab wanted to characterize these pathways by screening a library of single-gene-knockout Chlamy mutants with a series of drugs and analyzing how different groups of mutants responded to different families of molecules.

As a precursor to this project, I helped screen the effects of over 1000 small molecules on the growth of wild-type Chlamy colonies. These experiments were intended to help us figure out how much of each drug we would need to use in order to get a response out of the algae in the larger mutant screen, but the data also presented us with an opportunity to try and identify particular chemical substructures that were common amongst growth-inhibiting molecules--which would give us hints about how some of these drugs were actually retarding colony growth.

That was our motivation for clustering molecules by structural similarity--and it just happens that ECFPs are great for that purpose.

Now, software packages exist to compute similarity coefficients between molecules using the ECFP algorithm. But none of the vendors we contacted got back to us before we needed the clustering data--so I read up on the Rogers/Hahn paper and took a crack at implementing the algorithm myself.

## The Code
The main program